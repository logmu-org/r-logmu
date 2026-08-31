// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <algorithm>

namespace mortality
{

namespace
{
  // n must be non-negative!
  inline int T(int n) { return (n * (n + 1)) >> 1; }
}

// Step 1. Clamp `ib`
inline int clamp_ib(int ib, int ix_hi, int it_hi)
{
  return std::clamp(ib, -ix_hi, it_hi);
}

// Step 2. Get valid range for `it`. This is [min_it, max_it].
inline int min_it(int ib) { return std::max(0, ib); }
inline int max_it(int ib, int ix_hi, int it_hi) { return std::min(it_hi, ib + ix_hi); }

// Step 3. Given clamped `ib` get the offset for the first valid `it`.
// This function actually gets the offset given any valid `ib` and `it`.
inline int offset_for_valid(int ib, int it, int ix_hi, int it_hi)
{
  return T(std::min(it_hi, ib + ix_hi))
    - T(std::max(0, ib))
    + (it_hi + 1) * std::max(0, ib + ix_hi - it_hi)
    + it;
}

// Remap what in R is an x-t matrix (but in C++ looks like a t-x matrix)
// to a b-t matrix.
inline void remap(const double* values_tx, double* values_bt, int ix_hi, int it_hi)
{
  /*
  for(int ib = -ix_hi; ib <= it_hi; ++ib)
  {
    int it_min = std::max(0, ib);
    int it_max = std::min(it_hi, ix_hi + ib);
    for(int it = it_min; it <= it_max; ++it)
    {
      int ix = it - ib;
      *values_bt++ = values_tx[ix + it * (ix_hi + 1)];
    }
  }
   */

  const int stride = ix_hi + 1;
  const int stride_plus_1 = stride + 1;
  for (int ib = -ix_hi; ib <= it_hi; ++ib)
  {
    const int it_min = std::max(0, ib);
    const int it_max = std::min(it_hi, ix_hi + ib);
    const double* p_source = values_tx + (it_min - ib) + it_min * stride;
    for (int it = it_min; it <= it_max; ++it, p_source += stride_plus_1)
    {
      *values_bt++ = *p_source;
    }
  }
}

} // namespace mortality
