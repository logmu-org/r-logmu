// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#include "vec_ops.hpp"
#include <cstddef>
#include <cstring>

namespace tier
{

namespace {

using V_V_t = void (*)(vec_size, const double*, double*);

using VV_V_t = void (*)(vec_size, const double*, const double*, double*);
using VS_V_t = void (*)(vec_size, const double*,       double , double*);
using SV_V_t = void (*)(vec_size,       double , const double*, double*);

using VVV_V_t = void (*)(vec_size, const double*, const double*, const double*, double*);
using VVS_V_t = void (*)(vec_size, const double*, const double*,       double , double*);
using VSV_V_t = void (*)(vec_size, const double*,       double , const double*, double*);
using VSS_V_t = void (*)(vec_size, const double*,       double ,       double , double*);
using SVV_V_t = void (*)(vec_size,       double , const double*, const double*, double*);
using SVS_V_t = void (*)(vec_size,       double , const double*,       double , double*);
using SSV_V_t = void (*)(vec_size,       double ,       double , const double*, double*);

struct tier_table
{
  const char* name;
  int (*lanes)();

  V_V_t neg_V_V, exp_V_V, expm1_V_V, log_V_V, log1p_V_V, m_from_q_V_V;

  VV_V_t add_VV_V, sub_VV_V, mul_VV_V, div_VV_V, pow_VV_V, max_VV_V, min_VV_V;
  VS_V_t add_VS_V, sub_VS_V, mul_VS_V, div_VS_V, pow_VS_V, max_VS_V, min_VS_V;
  SV_V_t add_SV_V, sub_SV_V, mul_SV_V, div_SV_V, pow_SV_V, max_SV_V, min_SV_V;

  VSS_V_t clamp_VSS_V;
};

#define LOGMU_TIER_TABLE(NS)  \
{                             \
  #NS                         \
  , &NS::simd_lanes           \
  , &NS::neg_V_V, &NS::exp_V_V, &NS::expm1_V_V, &NS::log_V_V, &NS::log1p_V_V, &NS::m_from_q_V_V  \
  , &NS::add_VV_V, &NS::sub_VV_V, &NS::mul_VV_V, &NS::div_VV_V, &NS::pow_VV_V, &NS::max_VV_V, &NS::min_VV_V \
  , &NS::add_VS_V, &NS::sub_VS_V, &NS::mul_VS_V, &NS::div_VS_V, &NS::pow_VS_V, &NS::max_VS_V, &NS::min_VS_V \
  , &NS::add_SV_V, &NS::sub_SV_V, &NS::mul_SV_V, &NS::div_SV_V, &NS::pow_SV_V, &NS::max_SV_V, &NS::min_SV_V \
  , &NS::clamp_VSS_V \
}

const tier_table& select_tier()
{
  static const tier_table baseline_tbl = LOGMU_TIER_TABLE(baseline);

#if IS_X86_TARGET && (defined(__GNUC__) || defined(__clang__))
  static const tier_table avx2_tbl   = LOGMU_TIER_TABLE(avx2);
  static const tier_table avx512_tbl = LOGMU_TIER_TABLE(avx512);

  __builtin_cpu_init();
  // __builtin_cpu_supports() checks both CPUID feature bits and (for AVX
  // families) OS support via XGETBV, so a CPU with AVX-512 silicon but no OS
  // ZMM-save support correctly falls back.
  if (__builtin_cpu_supports("avx512f")) return avx512_tbl;
  if (__builtin_cpu_supports("avx2"))    return avx2_tbl;
#endif
  return baseline_tbl;
}

// Resolved exactly once (first call); the reference is cheap to re-read.
const tier_table& T = select_tier();

} // namespace

const char* active_tier() { return T.name; }
int active_lanes() { return static_cast<int>(T.lanes()); }

#define DEFINE_WRAPPER1(NAME) \
void NAME##_V_V(vec_size n, const double* a0, double* r) { T.NAME##_V_V(n, a0, r); }

DEFINE_WRAPPER1(neg)
DEFINE_WRAPPER1(exp)
DEFINE_WRAPPER1(expm1)
DEFINE_WRAPPER1(log)
DEFINE_WRAPPER1(log1p)
DEFINE_WRAPPER1(m_from_q)

#define DEFINE_WRAPPER2(NAME) \
void NAME##_VV_V(vec_size n, const double* a0, const double* a1, double* r) { T.NAME##_VV_V(n, a0, a1, r); } \
void NAME##_VS_V(vec_size n, const double* a0,       double  a1, double* r) { T.NAME##_VS_V(n, a0, a1, r); } \
void NAME##_SV_V(vec_size n,        double a0, const double* a1, double* r) { T.NAME##_SV_V(n, a0, a1, r); }

DEFINE_WRAPPER2(add)
DEFINE_WRAPPER2(sub)
DEFINE_WRAPPER2(mul)
DEFINE_WRAPPER2(div)
DEFINE_WRAPPER2(pow)
DEFINE_WRAPPER2(min)
DEFINE_WRAPPER2(max)

#define DEFINE_WRAPPER3(NAME) \
void NAME##_VVV_V(vec_size n, const double* a0, const double* a1, const double* a2, double* r) { T.NAME##_VVV_V(n, a0, a1, a2, r); } \
void NAME##_VVS_V(vec_size n, const double* a0, const double* a1,       double  a2, double* r) { T.NAME##_VVS_V(n, a0, a1, a2, r); } \
void NAME##_VSV_V(vec_size n, const double* a0,       double  a1, const double* a2, double* r) { T.NAME##_VSV_V(n, a0, a1, a2, r); } \
void NAME##_VSS_V(vec_size n, const double* a0,       double  a1,       double  a2, double* r) { T.NAME##_VSS_V(n, a0, a1, a2, r); } \
void NAME##_SVV_V(vec_size n,       double  a0, const double* a1, const double* a2, double* r) { T.NAME##_SVV_V(n, a0, a1, a2, r); } \
void NAME##_SVS_V(vec_size n,       double  a0, const double* a1,       double  a2, double* r) { T.NAME##_SVS_V(n, a0, a1, a2, r); } \
void NAME##_SSV_V(vec_size n,       double  a0,       double  a1, const double* a2, double* r) { T.NAME##_SSV_V(n, a0, a1, a2, r); }

void clamp_VSS_V(vec_size n, const double* a0,       double  a1,       double  a2, double* r) { T.clamp_VSS_V(n, a0, a1, a2, r); } \
//DEFINE_WRAPPER3(clamp)

} // namespace tier
