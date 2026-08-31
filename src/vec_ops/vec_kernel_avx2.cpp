// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#include "vec_ops.hpp"   // pulls in IS_X86_TARGET

#if IS_X86_TARGET
#include "vec_kernel.hpp"

namespace tier::avx2
{
DEFINE_TIER()
}

#endif
