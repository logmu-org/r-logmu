// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cmath>
#include <cstdint>

const int ValidDateStartYear = 1000;
const int ValidDateEndYear = 3000;
const int ValidDurationMaxYears = 2000;

const int ClicksPerYear = 534360;
const int ClicksPerHalfYear = ClicksPerYear / 2;
const int ClicksPerDay366 = ClicksPerYear / 366;
const int ClicksPerDay365 = ClicksPerYear / 365;

const double YearsPerClick = 1.0 / (double)ClicksPerYear;

const int ValidDateStartClicks = ClicksPerYear * ValidDateStartYear;
const int ValidDateEndClicks = ClicksPerYear * ValidDateEndYear;
const int ValidDurationMaxClicks = ClicksPerYear * ValidDurationMaxYears;

// Both take int64 because that is what a click is inside the veil engine, and a 32-bit click from a
// column widens to it without changing any answer. R's NA integer is the most negative int, so it
// widens to something far below the start of the calendar and is reported invalid, which is what
// makes these double as the missing-value test.
inline bool isValidDatey(int64_t clicks)
{
  return clicks >= ValidDateStartClicks && clicks <= ValidDateEndClicks;
}
inline bool isValidDurationy(int64_t clicks)
{
  // Don't use abs() -- NA_INTEGER is pathological case
  return clicks >= -ValidDurationMaxClicks && clicks <= ValidDurationMaxClicks;
}

// The canonical conversions between clicks and years. These two definitions are what a click and a
// year MEAN to each other, so every part of the system must use them rather than open-coding the
// arithmetic: the engine, the compile-time rewrites that reason about comparisons, and the R
// reference evaluator. Two of them computing a value differently is how a comparison ends up
// selecting different rows in R than in C++.
//
// Clicks to years is a multiply by the reciprocal, not a divide. A divide is not rewritten to a
// multiply without -ffast-math, which is unavailable here, and it is several times slower in
// vectorised code. It costs nothing measurable: whole years, ages, months and quarters all convert
// exactly, and only sub-month values can sit one unit in the last place from a correctly rounded
// divide. That exactness is a property of ClicksPerYear's particular value rather than a guarantee
// of the arithmetic, so it is pinned by tests -- see test-veil_conversion.R.
// Takes int64 because that is what a click is inside the engine. No value changes: a 32-bit click
// widens exactly, and every click is far inside the range a double represents exactly, so this is the
// same multiplication it always was.
inline double yearsFromClicks(int64_t clicks)
{
  return static_cast<double>(clicks) * YearsPerClick;
}

// Banker's rounding. WHEREVER ROUNDING IS DEFINITIONAL -- where it decides what a value means
// rather than merely how it is computed -- this is the one to use, so that logmu and datey can
// never disagree about a boundary. It is also IEEE 754-2019's roundToIntegralTiesToEven, which is
// what veil's `round` operator has to be.
//
// Taken from the datey package; keep the two in step. The only difference is the non-finite guard
// below, which datey needs as well: it cannot change any answer datey can reach, since a click is
// always finite, but leaving the two spellings different is how they start to drift.
//
// Deliberately not std::nearbyint or std::rint: both follow the CURRENT rounding mode, which is
// round-to-nearest-even only until something in the process calls fesetround. That is ambient
// state another package or a BLAS can change underneath us, and it is per-thread once the engine
// runs on more than one. std::round is no good either -- it rounds half away from zero, which is
// not what the specification says.
inline double roundBankers(double x)
{
  // An infinity rounds to itself and a NaN to itself, per IEEE 754-2019 roundToIntegralTiesToEven,
  // and both have to be taken out BEFORE the arithmetic below. std::fmod of an infinity is a NaN,
  // so the final sum would otherwise turn an infinity silently into a NaN. Within the datey
  // framework that could not arise -- a click is always finite -- but this is also veil's `round`
  // operator, whose argument is whatever a user expression produced, `1 / 0` included.

  if (!std::isfinite(x)) { return x; }

  double trunc_x = std::trunc(x);

  //abs_x=ABS(x-trunc_x)
  double abs_x = std::abs(x - trunc_x);

  //=IF(abs_x<0.5,trunc_x,IF(abs_x>0.5,trunc_x+SIGN(x),trunc_x+fmod_trunc_x_2))
  if (abs_x < 0.5) { return trunc_x; }
  if (abs_x > 0.5) { return trunc_x + std::copysign(1.0, x); }
  return trunc_x + std::fmod(trunc_x, 2.0);
}

// Years to clicks, per the datey specification: scale by ClicksPerYear, round half to even, then
// take the integer. It returns a double so the caller can range-check before narrowing, because
// converting an out-of-range double to int is undefined behaviour. The two conversions round-trip
// exactly -- yearsFromClicks then clicksFromYears returns every valid click unchanged, with around
// a millionfold margin -- which is what lets a rewrite ask whether a year value lands on a whole
// click by converting it and converting it back.
//
// This is a boundary and compile-time conversion, not an inner-loop one, so spelling the rounding
// out costs nothing worth measuring.
inline double clicksFromYears(double years)
{
  return roundBankers(years * ClicksPerYear);
}
