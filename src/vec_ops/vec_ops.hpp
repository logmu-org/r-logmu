// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#ifndef VEC_OPS_HPP
#define VEC_OPS_HPP

#include <cstddef>   // std::ptrdiff_t

// Is this an x86 target where AVX2 / AVX-512 tiers make sense?
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#  define IS_X86_TARGET 1
#else
#  define IS_X86_TARGET 0
#endif

namespace tier
{

// Length type for the vectorised kernels. std::ptrdiff_t is the R-free
// equivalent of R's R_xlen_t (itself ptrdiff_t on 64-bit R) and matches
// Python's Py_ssize_t / NumPy npy_intp, so each binding passes its native
// length type straight in.
using vec_size = std::ptrdiff_t;

#define DECLARE_WRAPPER1(NAME) \
void NAME##_V_V(vec_size, const double*, double*);

#define DECLARE_WRAPPER2(NAME) \
void NAME##_VV_V(vec_size, const double*, const double*, double*); \
void NAME##_VS_V(vec_size, const double*,       double , double*); \
void NAME##_SV_V(vec_size,       double , const double*, double*);

#define DECLARE_WRAPPER3(NAME) \
void NAME##_VVV_V(vec_size, const double*, const double*, const double*, double*); \
void NAME##_VVS_V(vec_size, const double*, const double*,       double , double*); \
void NAME##_VSV_V(vec_size, const double*,       double , const double*, double*); \
void NAME##_VSS_V(vec_size, const double*,       double ,       double , double*); \
void NAME##_SVV_V(vec_size,       double , const double*, const double*, double*); \
void NAME##_SVS_V(vec_size,       double , const double*,       double , double*); \
void NAME##_SSV_V(vec_size,       double ,       double , const double*, double*);


#define DECLARE_TIER(NS)                                       \
namespace NS                                                   \
{                                                              \
  int simd_lanes();                                            \
  DECLARE_WRAPPER1(neg)                                        \
  DECLARE_WRAPPER1(exp)                                        \
  DECLARE_WRAPPER1(expm1)                                      \
  DECLARE_WRAPPER1(log)                                        \
  DECLARE_WRAPPER1(log1p)                                      \
  DECLARE_WRAPPER1(m_from_q)                                   \
\
  DECLARE_WRAPPER2(add)                                        \
  DECLARE_WRAPPER2(sub)                                        \
  DECLARE_WRAPPER2(mul)                                        \
  DECLARE_WRAPPER2(div)                                        \
  DECLARE_WRAPPER2(min)                                        \
  DECLARE_WRAPPER2(max)                                        \
\
  void clamp_VSS_V(vec_size, const double*, double, double, double*);\
}

DECLARE_TIER(baseline)
#if IS_X86_TARGET
  DECLARE_TIER(avx2)
  DECLARE_TIER(avx512)
#endif

// Public entry points. Resolved once at load time to the best available tier
// (see dispatch.cpp).

// Name of tier selected ("avx512" / "avx2" / "baseline").
const char* active_tier();

// SIMD lane count of selected tier (2 baseline, 4 avx2, 8 avx512).
int active_lanes();

DECLARE_WRAPPER1(neg)
DECLARE_WRAPPER1(exp)
DECLARE_WRAPPER1(expm1)
DECLARE_WRAPPER1(log)
DECLARE_WRAPPER1(log1p)
DECLARE_WRAPPER1(m_from_q)

DECLARE_WRAPPER2(add)
DECLARE_WRAPPER2(sub)
DECLARE_WRAPPER2(mul)
DECLARE_WRAPPER2(div)
DECLARE_WRAPPER2(min)
DECLARE_WRAPPER2(max)

void clamp_VSS_V(vec_size, const double*, double, double, double*);

} // namespace tier

#endif // VEC_OPS_HPP
