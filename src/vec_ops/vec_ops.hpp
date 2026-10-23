// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#ifndef VEC_OPS_HPP
#define VEC_OPS_HPP

#include <cstddef>   // std::ptrdiff_t

// ARE THE AVX2 / AVX-512 TIERS PART OF THIS BUILD?
//
// THE BUILD DECIDES, NOT THIS HEADER. `configure` probes the compiler and picks
// the objects to compile; it passes that decision here as
// LOGMU_BUILD_SIMD_TIERS. Makevars.win sets it to 1, since Windows always
// builds all three.
//
// This used to make its own judgement from __x86_64__, independently of
// configure. The two could disagree -- and when configure's probe came back
// empty in a sanitizer container, they did: the AVX objects were never
// compiled while these declarations and the dispatch table still referred to
// them, so the package failed to load with
//
//     undefined symbol: _ZN4tier4avx210simd_lanesEv
//
// The fallback below is for a build that does not run configure at all.
#if defined(LOGMU_BUILD_SIMD_TIERS)
#  define IS_X86_TARGET LOGMU_BUILD_SIMD_TIERS
#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
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

// An annual mortality table as the lookup reads it: log mu in `tx` order, periods first, then ages,
// with the instant and age of the first cell in clicks.
struct annual_log_mu_table
{
  const double* log_mu;
  int period_count;
  int age_count;
  int t0_clicks;
  int x0_clicks;
};

// SCALAR FORMS, ONE PER TIER, AND THE REASON THEY ARE HERE AT ALL. A scalar `exp` compiled outside
// the tiers would be a different implementation from the vector kernel, or the same EVE code built
// with different ISA flags, and either way the two could disagree in the last place. Compiled inside
// each tier and chosen by the same dispatch, the scalar and vector forms of a tier run the same
// arithmetic and give identical answers, so whether a value took the scalar or the vector path can
// never show in a result.
#define DECLARE_SCALARS                                                \
double exp_S_S(double);                                                \
double log_S_S(double);                                                \
double table_log_mu_SS_S(const annual_log_mu_table&, int b_clicks, int t_clicks);

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
  DECLARE_SCALARS                                              \
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

DECLARE_SCALARS

// FOR TESTS ONLY: one named tier's scalar and vector forms, bypassing selection, so that every tier
// this CPU can run is checked in one process instead of one per `LOGMU_TIER` setting. Null where
// this build or this CPU lacks the tier; the CPU checks are the ones selection makes.
struct tier_forms
{
  double (*exp_S_S)(double);
  double (*log_S_S)(double);
  double (*table_log_mu_SS_S)(const annual_log_mu_table&, int, int);
  void (*exp_V_V)(vec_size, const double*, double*);
  void (*log_V_V)(vec_size, const double*, double*);
};
const tier_forms* forms_of_tier(const char* name);

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
