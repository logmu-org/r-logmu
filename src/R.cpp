// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#include <cpp11.hpp>
#include <algorithm>
#include <cmath>
#include "vec_ops/vec_ops.hpp"
#include "mortality/convert_m_to_mu.hpp"
#include "mortality/lookup_log_mu.hpp"

[[cpp11::register]]
bool cpp_all_finite(cpp11::doubles x)
{
  for (double val : x)
  {
    if (!std::isfinite(val)) { return false; }
  }

  return true;
}

[[cpp11::register]]
cpp11::doubles cpp_matrix_q_to_log_mu(cpp11::doubles q_matrix)
{
  R_xlen_t n = q_matrix.size();
  cpp11::writable::doubles result(n);

  if (n > 0)
  {
    auto matrix_dims = cpp11::as_cpp<cpp11::integers>(q_matrix.attr("dim"));

    // R is column major.
    // So what look like q[x,t] in R looks like q[t,x] in C++.
    // Individual dimensions are guaranteed by R to be `int`,
    // although total length in theory could be more that 2^32-1.
    // Given that we will need to multiple dimensions, it's safest to work with
    // `ptrdiff_t` everywhere.
    auto n_x = static_cast<ptrdiff_t>(matrix_dims[0]);
    auto n_t = static_cast<ptrdiff_t>(matrix_dims[1]);

    const double* p_q = REAL(q_matrix);
    double* p =  REAL(result);

    // 1. m = -log1p(-q)
    tier::m_from_q_V_V(n, p_q, p);

    // 2. mu = convexity adjusted m
    // NB C-style array argument order!
    mortality::convert_tx_array_from_m_to_mu(p, n_t, n_x);

    // 3. log mu = log mu
    tier::log_V_V(n, p, p);
  }

  return result;
}

[[cpp11::register]]
// Vectorised over individuals as well as over time: `b_clicks` and `t_clicks`
// recycle against each other, so one birth against many times reads a lifeline
// and many births against one time reads a population at an instant --
// `simulate()` does the latter, click by click.
//
// THE ONE DEFINITION in `mortality/lookup_log_mu.hpp`, which the veil engine compiles too.
//
// THE TABLE'S DIMENSIONS ARE NOT `n_x` AND `n_t`. They were until 2026-10-28, and the inner
// `n_t` shadowed the length of `t_clicks`, so the recycling below cycled through as many instants
// as the table has periods: the wrong instants when there were more of them, and a read past the
// end of `t_clicks` when there were fewer, as in `simulate()`'s one instant against many births.
cpp11::doubles cpp_slow_lookup_log_mu(
  cpp11::doubles log_mu_matrix,
  int x0_clicks, int t0_clicks, // xt order when interfacing with R
  cpp11::integers b_clicks, cpp11::integers t_clicks
)
{
  const R_xlen_t n_b = b_clicks.size();
  const R_xlen_t n_t = t_clicks.size();
  const R_xlen_t n = (n_b == 0 || n_t == 0) ? 0 : std::max(n_b, n_t);
  cpp11::writable::doubles result(n);

  if (n > 0)
  {
    auto matrix_dims = cpp11::as_cpp<cpp11::integers>(log_mu_matrix.attr("dim"));

    const int age_count = matrix_dims[0];
    const int period_count = matrix_dims[1];

    const double* p_log_mu = REAL(log_mu_matrix);
    double* p_result =  REAL(result);

    for (R_xlen_t i = 0; i < n; ++i)
    {
      p_result[i] = mortality::interpolatedLogMu(
        p_log_mu, period_count, age_count, // tx order when inside C++
        t0_clicks, x0_clicks,
        b_clicks[i % n_b], t_clicks[i % n_t]);
    }
  }

  return result;
}

[[cpp11::register]]
// One birth along many instants through the veil engine's bulk lookup, in the active SIMD tier.
// FOR THE TESTS: it is the only way to compare the bulk path with `cpp_slow_lookup_log_mu()` value
// by value, since the engine itself only ever reports integrals.
cpp11::doubles cpp_table_log_mu_series(
  cpp11::doubles log_mu_matrix,
  int x0_clicks, int t0_clicks,
  int b_clicks, cpp11::integers t_clicks
)
{
  auto matrix_dims = cpp11::as_cpp<cpp11::integers>(log_mu_matrix.attr("dim"));
  if (matrix_dims.size() != 2 || matrix_dims[0] < 1 || matrix_dims[1] < 1)
  {
    cpp11::stop("`log_mu_matrix` must be a non-empty matrix.");
  }
  if (!isValidDatey(b_clicks))
  {
    cpp11::stop("`b_clicks` must be a valid datey.");
  }
  for (int t : t_clicks)
  {
    if (!isValidDatey(t)) { cpp11::stop("`t_clicks` must all be valid dateys."); }
  }

  const R_xlen_t n = t_clicks.size();
  cpp11::writable::doubles result(n);

  const tier::annual_log_mu_table table{
    REAL(log_mu_matrix), matrix_dims[1], matrix_dims[0], t0_clicks, x0_clicks};
  tier::table_log_mu_SV_V(n, table, b_clicks, INTEGER(t_clicks), REAL(result));

  return result;
}
