// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#include <cpp11.hpp>
#include <cmath>

inline double deviance_residual(double A, double E, double V)
{
  // THE QUOTIENT COMES FIRST, AND THE ZERO TEST IS ON THE QUOTIENT, NOT ON `A`.
  //
  // Testing `A_over_E` rather than `A` catches A/E underflowing to zero as well
  // as A being zero. They are the same limit -- the bracket A*log(A/E) - (A - E)
  // tends to E either way -- but an underflowed quotient reaching the general
  // branch below would compute 0 * log(0), i.e. 0 * -inf, and return NaN.
  //
  // Forming `alpha` from the rounded quotient is also load-bearing: `A_over_E`
  // lies within a factor of two of 1 on the branch that uses `alpha`, so the
  // subtraction is exact by Sterbenz, and that exactness is what lets the
  // general branch cancel gracefully. Rewriting these lines as
  // `alpha = (A - E) / E` costs about three orders of magnitude and fails no
  // test.
  double A_over_E = A / E;

  if (A_over_E == 0.0)
  {
    return -E * std::sqrt(2.0 / V);
  }

  double alpha = A_over_E - 1.0;

  if (std::abs(alpha) > 0.0015)
  {
    double T1 = A_over_E * std::log(A_over_E) - alpha;
    double T2 = E * std::sqrt(2.0 * T1 / V);
    return std::copysign(T2, alpha);
  }
  else
  {
    const double MinusOneSixth = -1.0 / 6.0;
    const double FiveSeventySeconds = 5.0 / 72.0;
    const double MinusEightyThreeOver2160 = -83.0 / 2160.0;
    double quadratic = FiveSeventySeconds + alpha * MinusEightyThreeOver2160;
    double polynomial = alpha * (1.0 + alpha * (MinusOneSixth + alpha * quadratic));
    return E / sqrt(V) * polynomial;
  }
}

[[cpp11::register]]
cpp11::doubles cpp_deviance_residual(cpp11::doubles A, cpp11::doubles E, cpp11::doubles V)
{
  R_xlen_t n = A.size();

  if (E.size() != n || V.size() != n) { cpp11::stop("A, E and V must be the same size."); }

  cpp11::writable::doubles result(n);

  for (R_xlen_t i = 0; i < n; ++i)
  {
    result[i] = deviance_residual(A[i], E[i], V[i]);
  }

  return result;
}

[[cpp11::register]]
void cpp_validate_aev(cpp11::doubles A, cpp11::doubles E, cpp11::doubles V)
{
  R_xlen_t n = A.size();

  if (E.size() != n || V.size() != n) { cpp11::stop("A, E and V must be the same size."); }

  for (R_xlen_t i = 0; i < n; ++i)
  {
    double A_i = A[i];
    double E_i = E[i];
    double V_i = V[i];

    if (A_i < 0.0) { cpp11::stop("A cannot be negative."); }
    if (E_i < 0.0) { cpp11::stop("E cannot be negative."); }
    if (V_i < 0.0) { cpp11::stop("V cannot be negative."); }

    if (E_i > 0.0 && V_i == 0.0) { cpp11::stop("E cannot be non-zero if V is zero."); }
    if (V_i > 0.0 && E_i == 0.0) { cpp11::stop("V cannot be non-zero if E is zero."); }
    if (A_i > 0.0 && E_i == 0.0) { cpp11::stop("A cannot be non-zero if E and V are zero."); }
  }
}

