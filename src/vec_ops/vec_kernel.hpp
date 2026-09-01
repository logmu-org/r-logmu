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
DEFINE_WRAPPER2(min)                                           \
DEFINE_WRAPPER2(max)                                           \
\
void clamp_VSS_V(vec_size n, const double* arg0, double arg1, double arg2, double* result){ ::tier::VSS_V(n,arg0,arg1,arg2,result,eve::clamp); }

#endif // VEC_KERNEL_HPP
