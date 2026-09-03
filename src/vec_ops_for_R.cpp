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
// tier dispatch table. It has ONE vectorised shape -- a scalar base over a
// vector of exponents -- which calls the `mul` and `exp` kernels directly.
//
// EVE'S POW WAS WRONG, AND THAT IS WHY IT IS GONE. Edge-case tests written
// before any change failed SIXTEEN IEEE cases on the AVX-512 tier:
// `pow(NaN, Inf)` gave Inf, `pow(-2, 0.5)` was not NaN. The old tests covered
// bases in (0, 1) and exponents in (-1, 1), so none of it had ever run.
// `std::pow` is correctly rounded and gets the whole IEEE ladder right for
// free -- `x^0` is 1 for every x including NaN, `1^y` is 1 for every y
// including Inf, a negative base keeps its sign for integer exponents. The
// engine already evaluates `Op::Pow` with `std::pow`, so the two agree.
//
// EVE'S POW DID NOT CRASH. THAT WAS A MISATTRIBUTION, corrected 2026-09-02.
// The Windows segfault was `vmovapd %ymm3,(%rcx)` inside
// `tier::avx2::exp_V_V`: an aligned 256-bit store to a stack address, where
// Windows guarantees only 16-byte alignment, so it faulted about half the
// time. The gdb backtrace that appeared to name `pow_VV_V` was resolving to
// the nearest EXPORT rather than the real function, which is what symbol
// lookup inside a stripped DLL does. `SIMD_SAFE_STACK` in `Makevars.win`
// fixes it, and `test-simd-stack-alignment.R` guards it.
//
// SO THE SCALAR-BASE FAST PATH IS BACK, AND IT WRITES IN PLACE -- the same
// buffer passed as both `const double* arg` and `double* result`. That is
// sound, and was never shown to be otherwise. No `restrict` appears anywhere
// in `src/vec_ops/`, so the compiler is told nothing about non-overlap; both
// parameters are `double*`, so no strict-aliasing question arises; and every
// lane loads and stores the same index. It was removed on 2026-09-01 purely on
// the crash diagnosis above, which was wrong.
//
// THE GUARD IS x FINITE, POSITIVE AND NOT 1. Inside it every special exponent
// falls out of exp(y log x) with no further cases. log 2 is positive, so
// y = Inf gives exp(Inf) = Inf and y = -Inf gives exp(-Inf) = 0; log 0.5 is
// negative and reverses both; y = 0 gives exp(0) = 1; NaN propagates. Outside
// the guard -- a negative, zero, infinite or NaN base -- `std::pow` handles it.
//
// `x == 1` NEEDS ITS OWN BRANCH, and it is not an optimisation. log 1 is 0 and
// Inf * 0 is NaN, so any exp(y log x) formulation returns NaN where
// `pow(1, Inf)` must be 1.
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
      // y * log(base) into the result buffer, then exp over that same buffer.
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
