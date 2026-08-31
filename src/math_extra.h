// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstdlib>

struct EuclideanDivMod
{
  int div;
  int mod;
};

// Assumes divisor > 0
inline EuclideanDivMod divmod(int numerator, int divisor)
{
  std::div_t trunc = std::div(numerator, divisor);

  // If remainder is negative, shift it forward by 'divisor' and decrement the quotient
  int is_neg = (trunc.rem < 0);

  // Branchless addition of 'divisor' if remainder is negative
  return {trunc.quot - is_neg, trunc.rem + (divisor & -is_neg)};
}
