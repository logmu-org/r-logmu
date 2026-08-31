// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <cmath>
#include <algorithm>

namespace mortality
{

//const double MachineEpsilon = 2.2204460492503131E-16;
//const double MachineUnitRoundOff = MachineEpsilon / 2.0;
//const double QCeiling = 1.0 - MachineUnitRoundOff;
//const double QFloor = MachineUnitRoundOff;

// Converts the interior of a vector of m = −log(1 − q) rates to μ rates
// in place with an arbitrary stride allowing for convexity.
inline void convert_vector_m_to_mu(double* p, ptrdiff_t count, ptrdiff_t stride)
{
  if (count <= 2) { return; }

  // Formula is
  //
  //  m ≈ μ + β/24
  //
  // where β is d²μ/dt²
  //
  // Hence
  //
  //  μ ≈ m − (m₊ − 2m + m₋) / 24
  //
  //  Worst proportional CMI base table convexity I've seen is
  //  8.6% for PFA08 at age 91. That would be an 0.4% adjustment.
  //  So using a floor of 50% of m to ensure positivity is highly unlikely to be
  //  hit in practice.

  const double one_24th = 1.0 / 24.0;

  double m_lo = *p;
  p += stride;
  double m = *p;

  for (ptrdiff_t i = 2; i < count; ++i)
  {
    auto p_mu = p;

    p += stride;
    double m_hi = *p;

    double convexity_adj = (m_hi - 2.0 * m + m_lo) * one_24th;
    *p_mu = std::max(m - convexity_adj, 0.5 * m);

    m_lo = m;
    m = m_hi;
  }
}

// Converts the interior and age edges of a block of m[t,x]
// to μ[t,x] in place
// allowing for convexity.
inline void convert_tx_array_from_m_to_mu(double* p, ptrdiff_t n_t, ptrdiff_t n_x)
{
  if (n_x <= 2 || n_t <= 0) { return; }

  // □ ● ● ● ● ● ● □
  // □ \ \ \ \ \ \ □
  // □ \ \ \ \ \ \ □
  // □ \ \ \ \ \ \ □
  // □ ● ● ● ● ● ● □

  auto stride = n_x + 1;

  if (n_t > 2)
  {
    // Diagonals starting on the top row
    auto p_diagonal = p;
    for (ptrdiff_t i = 0; i < n_x - 2; ++i)
    {
      auto count = std::min(n_t, n_x - i);
      convert_vector_m_to_mu(p_diagonal, count, stride);
      ++p_diagonal;
    }
    // Diagonals starting on the left column
    // (omitting the main diagonal, which was handled above)
    p_diagonal = p + n_x;
    for (ptrdiff_t i = 1; i < n_t - 2; ++i)
    {
      auto count = std::min(n_x, n_t - i);
      convert_vector_m_to_mu(p_diagonal, count, stride);
      p_diagonal += n_x;
    }
  }

  // Do the age edges *after* the interior
  // (because the interior op depends on the edges).

  // Row 0:
  convert_vector_m_to_mu(p, n_x, 1);

  if (n_t >= 2)
  {
    // Row n_t - 1:
    convert_vector_m_to_mu(p + (n_t - 1) * n_x, n_x, 1);
  }
}

} // namespace mortality

