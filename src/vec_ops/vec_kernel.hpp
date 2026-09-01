// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#ifndef VEC_KERNEL_HPP
#define VEC_KERNEL_HPP

#include "vec_ops.hpp" // tier::vec_size
#include <cstddef>
#include <cmath>     // std::pow, std::log, std::isfinite
#include <eve/wide.hpp>
#include <eve/module/core.hpp>
#include <eve/module/math.hpp>

namespace tier
{

template <typename Op>
inline void V_V(vec_size n, const double* arg, double* result, Op op)
{
  using simd = eve::wide<double>;
  constexpr vec_size width = static_cast<vec_size>(simd::size());

  vec_size i = 0;
  for (; i + width <= n; i += width)
  {
    // Unaligned loads and stores
    simd v(arg + i);
    eve::store(op(v), result + i);
  }
  for (; i < n; ++i) { result[i] = op(arg[i]); }
}

template <typename Op0, typename Op1>
inline void V_Vx2(vec_size n, const double* arg, double* result, Op0 op0, Op1 op1)
{
  using simd = eve::wide<double>;
  constexpr vec_size width = static_cast<vec_size>(simd::size());

  vec_size i = 0;
  for (; i + width <= n; i += width)
  {
    // Unaligned loads and stores
    simd v(arg + i);
    eve::store(op1(op0(v)), result + i);
  }
  for (; i < n; ++i) { result[i] = op1(op0(arg[i])); }
}

template <typename Op0, typename Op1, typename Op2>
inline void V_Vx3(vec_size n, const double* arg, double* result, Op0 op0, Op1 op1, Op2 op2)
{
  using simd = eve::wide<double>;
  constexpr vec_size width = static_cast<vec_size>(simd::size());

  vec_size i = 0;
  for (; i + width <= n; i += width)
  {
    // Unaligned loads and stores
    simd v(arg + i);
    eve::store(op2(op1(op0(v))), result + i);
  }
  for (; i < n; ++i) { result[i] = op2(op1(op0(arg[i]))); }
}

template <typename Op>
inline void VV_V(vec_size n, const double* arg0, const double* arg1, double* result, Op op)
{
  using simd = eve::wide<double>;
  constexpr vec_size width = static_cast<vec_size>(simd::size());

  vec_size i = 0;
  for (; i + width <= n; i += width)
  {
    // Unaligned loads and stores
    simd v0(arg0 + i);
    simd v1(arg1 + i);
    eve::store(op(v0, v1), result + i);
  }
  for (; i < n; ++i) { result[i] = op(arg0[i], arg1[i]); }
}
template <typename Op>
inline void VS_V(vec_size n, const double* arg0, double arg1, double* result, Op op)
{
  using simd = eve::wide<double>;
  constexpr vec_size width = static_cast<vec_size>(simd::size());

  eve::wide<double> v1(arg1);
  vec_size i = 0;
  for (; i + width <= n; i += width)
  {
    // Unaligned loads and stores
    simd v0(arg0 + i);
    eve::store(op(v0, v1), result + i);
  }
  for (; i < n; ++i) { result[i] = op(arg0[i], arg1); }
}
template <typename Op>
inline void SV_V(vec_size n, double arg0, const double* arg1, double* result, Op op)
{
  using simd = eve::wide<double>;
  constexpr vec_size width = static_cast<vec_size>(simd::size());

  eve::wide<double> v0(arg0);
  vec_size i = 0;
  for (; i + width <= n; i += width)
  {
    // Unaligned loads and stores
    simd v1(arg1 + i);
    eve::store(op(v0, v1), result + i);
  }
  for (; i < n; ++i) { result[i] = op(arg0, arg1[i]); }
}


template <typename Op>
inline void VSS_V(vec_size n, const double* arg0, double arg1, double arg2, double* result, Op op)
{
  using simd = eve::wide<double>;
  constexpr vec_size width = static_cast<vec_size>(simd::size());

  eve::wide<double> v1(arg1);
  eve::wide<double> v2(arg2);
  vec_size i = 0;
  for (; i + width <= n; i += width)
  {
    // Unaligned loads and stores
    simd v0(arg0 + i);
    eve::store(op(v0, v1, v2), result + i);
  }
  for (; i < n; ++i) { result[i] = op(arg0[i], arg1, arg2); }
}

} // namespace tier


#define DEFINE_WRAPPER1(NAME) \
void NAME##_V_V  (vec_size n, const double* arg, double* result){ ::tier::V_V(n,arg,result,eve::NAME); }

#define DEFINE_WRAPPER1x1(NAME, OP) \
void NAME##_V_V  (vec_size n, const double* arg, double* result){ ::tier::V_V(n,arg,result,eve::OP); }

#define DEFINE_WRAPPER1x2(NAME, OP0, OP1) \
void NAME##_V_V  (vec_size n, const double* arg, double* result){ ::tier::V_Vx3(n,arg,result,eve::OP0,eve::OP1); }

#define DEFINE_WRAPPER1x3(NAME, OP0, OP1, OP2) \
void NAME##_V_V  (vec_size n, const double* arg, double* result){ ::tier::V_Vx3(n,arg,result,eve::OP0,eve::OP1,eve::OP2); }

#define DEFINE_WRAPPER2(NAME)                                                                                                       \
void NAME##_VV_V(vec_size n, const double* arg0, const double* arg1, double* result){ ::tier::VV_V(n,arg0,arg1,result,eve::NAME); } \
void NAME##_VS_V(vec_size n, const double* arg0,       double  arg1, double* result){ ::tier::VS_V(n,arg0,arg1,result,eve::NAME); } \
void NAME##_SV_V(vec_size n,       double  arg0, const double* arg1, double* result){ ::tier::SV_V(n,arg0,arg1,result,eve::NAME); }


// POW IS NOT EVE'S, AND IS THE ONE OPERATION HERE THAT IS NOT VECTORISED
// THROUGHOUT. Two reasons, both found the hard way on 2026-08-31.
//
// IT CRASHED. `eve::pow` on the AVX2 tier segfaulted on a GitHub Windows
// runner -- inside `pow_VV_V`, on exactly one full four-lane vector, confirmed
// by a gdb backtrace. It could not be reproduced on any machine here, on either
// tier, with or without `gctorture`.
//
// IT WAS ALSO WRONG. Edge-case tests written before this change failed 16 ways
// on the WORKING AVX-512 tier: `pow(NaN, Inf)` gave Inf, `pow(-2, 0.5)` was not
// NaN, and more. The old tests only covered bases in (0, 1) and exponents in
// (-1, 1), so none of it had ever been exercised.
//
// `std::pow` is correctly rounded and gets the whole IEEE 754 special-case
// ladder right for nothing -- `x^0` is 1 for every x including NaN, `1^y` is 1
// for every y including Inf, a negative base keeps its sign for integer
// exponents. Reimplementing that vectorised is a great deal of branchless work
// for an operation the engine does not even use: `Interpreter.hpp` evaluates
// `Op::Pow` with `std::pow` already.
//
// THE ONE VECTORISED CASE is a scalar base, which is the shape a projection
// factor takes -- `1.02 ^ t` over a vector of durations. Tim's idea, and it is
// sound precisely because the guard is one check OUTSIDE the loop and `log x`
// is computed once in full scalar precision. With `x` finite, strictly positive
// and not 1, every special exponent then falls out correctly:
//
//     y = 0     ->  0 * finite = 0     ->  exp(0)  = 1     correct
//     y = +Inf  ->  +-Inf              ->  Inf or 0        correct
//     y = -Inf  ->  -+Inf              ->  0 or Inf        correct
//     y = NaN   ->  NaN                ->  NaN             correct
//
// `x == 1` must be excluded because `log 1` is 0 and `Inf * 0` is NaN, where
// `pow(1, Inf)` is 1. It is filled directly instead.
//
// THE COST is accuracy: relative error is about |y log x| * epsilon, because
// error in the exponent becomes relative error after `exp`. A few ulp for
// ordinary arguments, a few hundred near the edge of `exp`'s range. The tests
// therefore require exactness for the special values and tolerance elsewhere.
#define DEFINE_POW()                                                          \
void pow_VV_V(vec_size n, const double* arg0, const double* arg1, double* result) \
{                                                                             \
  for (vec_size i = 0; i < n; ++i) { result[i] = std::pow(arg0[i], arg1[i]); }\
}                                                                             \
void pow_VS_V(vec_size n, const double* arg0, double arg1, double* result)    \
{                                                                             \
  for (vec_size i = 0; i < n; ++i) { result[i] = std::pow(arg0[i], arg1); }   \
}                                                                             \
void pow_SV_V(vec_size n, double arg0, const double* arg1, double* result)    \
{                                                                             \
  if (arg0 == 1.0)                                                            \
  {                                                                           \
    for (vec_size i = 0; i < n; ++i) { result[i] = 1.0; }                     \
    return;                                                                   \
  }                                                                           \
                                                                              \
  if (std::isfinite(arg0) && arg0 > 0.0)                                      \
  {                                                                           \
    const double log_arg0 = std::log(arg0);                                   \
    ::tier::VS_V(n, arg1, log_arg0, result, eve::mul);                        \
    ::tier::V_V(n, result, result, eve::exp);                                 \
    return;                                                                   \
  }                                                                           \
                                                                              \
  for (vec_size i = 0; i < n; ++i) { result[i] = std::pow(arg0, arg1[i]); }   \
}

#define DEFINE_TIER()                                          \
int simd_lanes(){ return static_cast<int>(eve::wide<double>::size()); } \
\
DEFINE_WRAPPER1x1(neg,minus)                                   \
DEFINE_WRAPPER1(exp)                                           \
DEFINE_WRAPPER1(expm1)                                         \
DEFINE_WRAPPER1(log)                                           \
DEFINE_WRAPPER1(log1p)                                         \
DEFINE_WRAPPER1x3(m_from_q,minus,log1p,minus)                  \
\
DEFINE_WRAPPER2(add)                                           \
DEFINE_WRAPPER2(sub)                                           \
DEFINE_WRAPPER2(mul)                                           \
DEFINE_WRAPPER2(div)                                           \
DEFINE_POW()                                                   \
DEFINE_WRAPPER2(min)                                           \
DEFINE_WRAPPER2(max)                                           \
\
void clamp_VSS_V(vec_size n, const double* arg0, double arg1, double arg2, double* result){ ::tier::VSS_V(n,arg0,arg1,arg2,result,eve::clamp); }

#endif // VEC_KERNEL_HPP
