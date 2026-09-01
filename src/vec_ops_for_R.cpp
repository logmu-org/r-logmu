// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

// cpp11-registered R entry points. Compiled at baseline (no SIMD headers here);
// the actual SIMD work lives behind the tier:: dispatch functions.
//
// NOTE: these [[cpp11::register]] signatures must be literal (the cpp11 code
// generator scans source text and does NOT run the preprocessor), so they
// cannot be produced by a macro. The shared body lives in R_vec_unary().

#include <cpp11.hpp>
#include "vec_ops/vec_ops.hpp"
#include <algorithm>   // std::min
#include <cstddef>
#include <cmath>     // std::pow, std::log, std::exp, std::isfinite

[[cpp11::register]] cpp11::strings cpp_vec_active_tier()
{
  cpp11::writable::strings s(1);
  s[0] = tier::active_tier();
  return s;
}

[[cpp11::register]] int cpp_vec_active_lanes()
{
  return tier::active_lanes();
}

#define R_UNARY(NAME)                                          \
R_xlen_t n = x.size();                                         \
cpp11::writable::doubles result(n);                            \
if (n > 0) { tier::NAME##_V_V(n, REAL(x), REAL(result)); }     \
return result;

#define R_BINARY(NAME)                                         \
R_xlen_t n_x = x.size();                                       \
R_xlen_t n_y = y.size();                                       \
R_xlen_t n = std::max(n_x, n_y);                               \
cpp11::writable::doubles result(n);                            \
if (n > 0)                                                     \
{                                                              \
  if (n_x == n_y)    { tier::NAME##_VV_V(n, REAL(x), REAL(y), REAL(result)); }\
  else if (n_y == 1) { tier::NAME##_VS_V(n, REAL(x),    y[0], REAL(result)); }\
  else if (n_x == 1) { tier::NAME##_SV_V(n,    x[0], REAL(y), REAL(result)); }\
  else { cpp11::stop("Lengths of `x` and `y` must be equal or 1."); }\
}                                                              \
return result;

// Scan for [[cpp11::register]] takes place *before* macro expansion.
// So we can't replace all of this with macros.

[[cpp11::register]] cpp11::doubles cpp_vec_neg  (cpp11::doubles x){ R_UNARY(neg) }
[[cpp11::register]] cpp11::doubles cpp_vec_exp  (cpp11::doubles x){ R_UNARY(exp) }
[[cpp11::register]] cpp11::doubles cpp_vec_expm1(cpp11::doubles x){ R_UNARY(expm1) }
[[cpp11::register]] cpp11::doubles cpp_vec_log  (cpp11::doubles x){ R_UNARY(log) }
[[cpp11::register]] cpp11::doubles cpp_vec_log1p(cpp11::doubles x){ R_UNARY(log1p) }
[[cpp11::register]] cpp11::doubles cpp_vec_m_from_q(cpp11::doubles x){ R_UNARY(m_from_q) }

/*
// This approach is ~17% slower
[[cpp11::register]] cpp11::doubles cpp_vec_m_from_q2(cpp11::doubles x)
{
  R_xlen_t n = x.size();
  cpp11::writable::doubles result(n);
  if (n > 0)
  {
    const double* p_x = REAL(x);
    double* p_result =  REAL(result);
    tier::neg_V_V(n, p_x, p_result);
    tier::log1p_V_V(n, p_result, p_result);
    tier::neg_V_V(n, p_result, p_result);
  }
  return result;
}
*/

[[cpp11::register]] cpp11::doubles cpp_vec_add(cpp11::doubles x, cpp11::doubles y){ R_BINARY(add) }
[[cpp11::register]] cpp11::doubles cpp_vec_sub(cpp11::doubles x, cpp11::doubles y){ R_BINARY(sub) }
[[cpp11::register]] cpp11::doubles cpp_vec_mul(cpp11::doubles x, cpp11::doubles y){ R_BINARY(mul) }
[[cpp11::register]] cpp11::doubles cpp_vec_div(cpp11::doubles x, cpp11::doubles y){ R_BINARY(div) }
// POW IS IMPLEMENTED HERE AND NOWHERE ELSE, and does not go through the SIMD
// tier system at all. Three reasons, all found the hard way over 2026-08-31
// and 2026-09-01.
//
// EVE'S POW CRASHED. On the AVX2 tier it segfaulted on a GitHub Windows
// runner -- inside `pow_VV_V`, confirmed by a gdb backtrace after exit code
// 139 -- and could not be reproduced on any machine here, on either tier, with
// or without gctorture.
//
// EVE'S POW WAS ALSO WRONG. Edge-case tests written before any change failed
// SIXTEEN IEEE cases on the WORKING AVX-512 tier: `pow(NaN, Inf)` gave Inf,
// `pow(-2, 0.5)` was not NaN. The old tests covered bases in (0, 1) and
// exponents in (-1, 1), so none of it had ever been exercised.
//
// AND POW HAS NO BUSINESS IN THE TIER SYSTEM. Only one of its three shapes is
// vectorised at all, so a dispatch table entry bought nothing. `std::pow` is
// correctly rounded and gets the whole IEEE special-case ladder right for
// free -- `x^0` is 1 for every x including NaN, `1^y` is 1 for every y
// including Inf, a negative base keeps its sign for integer exponents. The
// engine already evaluates `Op::Pow` with `std::pow`, so the two now agree.
//
// THE ONE VECTORISED SHAPE is a scalar base, which is what a projection factor
// looks like: `1.02 ^ t` over a vector of durations. Tim's idea, and it works
// because the guard is one check OUTSIDE the loop and `log x` is computed once
// in full scalar precision. With `x` finite, strictly positive and not 1, every
// special exponent then falls out correctly:
//
//     y = 0     ->  0 * finite = 0  ->  exp(0) = 1      correct
//     y = +Inf  ->  +-Inf           ->  Inf or 0        correct
//     y = -Inf  ->  -+Inf           ->  0 or Inf        correct
//     y = NaN   ->  NaN             ->  NaN             correct
//
// `x == 1` is filled directly: `log 1` is 0 and `Inf * 0` is NaN, where
// `pow(1, Inf)` is 1. That path still uses the vectorised `mul` and `exp`
// kernels, which are exercised everywhere and have never been implicated.
//
// THE COST is accuracy on that one path: relative error is about
// |y log x| * epsilon, since error in the exponent becomes relative error after
// `exp`. A few ulp for ordinary arguments. The tests require exactness for the
// special values and tolerance elsewhere.
[[cpp11::register]] cpp11::doubles cpp_vec_pow(cpp11::doubles x, cpp11::doubles y)
{
  R_xlen_t n_x = x.size();
  R_xlen_t n_y = y.size();
  R_xlen_t n = std::max(n_x, n_y);

  cpp11::writable::doubles result(n);
  if (n == 0) { return result; }

  if (n_x != n_y && n_x != 1 && n_y != 1)
  {
    cpp11::stop("Lengths of `x` and `y` must be equal or 1.");
  }

  const double* p_x = REAL(x);
  const double* p_y = REAL(y);
  double* p_result = REAL(result);

  if (n_x == n_y)
  {
    for (R_xlen_t i = 0; i < n; ++i) { p_result[i] = std::pow(p_x[i], p_y[i]); }
  }
  else if (n_y == 1)
  {
    const double exponent = p_y[0];
    for (R_xlen_t i = 0; i < n; ++i) { p_result[i] = std::pow(p_x[i], exponent); }
  }
  else
  {
    const double base = p_x[0];

    if (base == 1.0)
    {
      for (R_xlen_t i = 0; i < n; ++i) { p_result[i] = 1.0; }
    }
    else if (std::isfinite(base) && base > 0.0)
    {
      tier::mul_VS_V(n, p_y, std::log(base), p_result);
      tier::exp_V_V(n, p_result, p_result);
    }
    else
    {
      for (R_xlen_t i = 0; i < n; ++i) { p_result[i] = std::pow(base, p_y[i]); }
    }
  }

  return result;
}
[[cpp11::register]] cpp11::doubles cpp_vec_max(cpp11::doubles x, cpp11::doubles y){ R_BINARY(max) }
[[cpp11::register]] cpp11::doubles cpp_vec_min(cpp11::doubles x, cpp11::doubles y){ R_BINARY(min) }

[[cpp11::register]] cpp11::doubles cpp_vec_clamp(cpp11::doubles x, double min, double max)
{
  R_xlen_t n = x.size();
  cpp11::writable::doubles result(n);
  if (n > 0) { tier::clamp_VSS_V(n, REAL(x), min, max, REAL(result)); }
  return result;
}
