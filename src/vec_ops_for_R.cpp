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
[[cpp11::register]] cpp11::doubles cpp_vec_pow(cpp11::doubles x, cpp11::doubles y){ R_BINARY(pow) }
[[cpp11::register]] cpp11::doubles cpp_vec_max(cpp11::doubles x, cpp11::doubles y){ R_BINARY(max) }
[[cpp11::register]] cpp11::doubles cpp_vec_min(cpp11::doubles x, cpp11::doubles y){ R_BINARY(min) }

[[cpp11::register]] cpp11::doubles cpp_vec_clamp(cpp11::doubles x, double min, double max)
{
  R_xlen_t n = x.size();
  cpp11::writable::doubles result(n);
  if (n > 0) { tier::clamp_VSS_V(n, REAL(x), min, max, REAL(result)); }
  return result;
}
