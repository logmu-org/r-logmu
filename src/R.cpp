// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#include <cpp11.hpp>
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
cpp11::doubles cpp_slow_lookup_log_mu(
  cpp11::doubles log_mu_matrix,
  int x0_clicks, int t0_clicks, // xt order when interfacing with R
  int b_clicks, cpp11::integers t_clicks
)
{
  R_xlen_t n = t_clicks.size();
  cpp11::writable::doubles result(n);

  if (n > 0)
  {
    auto matrix_dims = cpp11::as_cpp<cpp11::integers>(log_mu_matrix.attr("dim"));

    int n_x = matrix_dims[0];
    int n_t = matrix_dims[1];

    const double* p_log_mu =  REAL(log_mu_matrix);
    double* p_result =  REAL(result);

    for (R_xlen_t i = 0; i < n; ++i)
    {
      p_result[i] = mortality::get_interpolated_log_mu_bt_from_annual_tx_array(
        p_log_mu, n_t, n_x, // tx order when inside C++
        t0_clicks, x0_clicks,
        b_clicks, t_clicks[i]
      );

    }
  }

  return result;
}
