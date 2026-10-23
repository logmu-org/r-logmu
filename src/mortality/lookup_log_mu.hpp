// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

// THE DEFINITION OF A TABLE'S LOG MU AT ONE BIRTH AND ONE INSTANT.
//
// INCLUDED ONLY BY THE TIER KERNELS (`vec_ops/vec_kernel.hpp`), and reached from everywhere else
// through `tier::table_log_mu_SS_S`. R's `log_mu()` and the veil engine therefore run the same
// compiled arithmetic on any one machine, and the bulk path, when it comes, is written from the
// same steps. Including it anywhere else would give that caller a second, differently compiled
// copy of the arithmetic, and the two would no longer be guaranteed to agree.
//
// AN ANONYMOUS NAMESPACE, DELIBERATELY. Each tier's translation unit is compiled with its own ISA
// flags, so an ordinary `inline` function here would be one symbol with three different bodies, and
// the linker keeps whichever it sees first. A baseline CPU could then run the AVX-512 body. Internal
// linkage gives each tier its own copy.
//
// COHORT FIRST, THEN PERIOD, EACH STEP AN `eve::fma`. The lattice is cohort-period, so a person's
// whole exposure lies between two cohort columns fixed for life: blending those first leaves one
// interpolation in period per instant, which is the shape the bulk path needs. `eve::fma` is a fused
// multiply-add on a tier that has the instruction and a separate multiply and add on one that does
// not. Either way the scalar and vector forms of a tier do the same thing, which is the property that
// matters; the tiers may differ from each other by an ulp, which does not. `std::lerp` is not used:
// its guarantees at t = 1 and on monotonicity buy nothing for a weight in [0, 1), its branches get in
// the way of lanes, and the standard fixes only its properties, so libraries may round it differently.
// `std::fma` is not used either: compiled without `-mfma`, Rtools turns it into a call to a library
// `fma` that is not correctly rounded.
//
// AN INDEX PAST THE END IS CLAMPED to the edge of the table, age and period independently. An
// individual older than the table's last age is read at that last age rather than refused, which is
// consistent with `end_age` being a check the caller makes rather than an input here.

#include <algorithm>
#include <eve/module/core.hpp>
#include "datey.h"   // ClicksPerYear

namespace mortality
{
namespace
{

struct WholeYearsAndClicks
{
	int years;
	int clicks;
};

// Floor division by a year, so a remainder is never negative. Kept here rather than shared with
// `math_extra.h` for the linkage reason in the header comment.
inline WholeYearsAndClicks splitClicks(int clicks)
{
	int years = clicks / ClicksPerYear;
	int remainder = clicks % ClicksPerYear;
	if (remainder < 0)
	{
		remainder += ClicksPerYear;
		--years;
	}
	return {years, remainder};
}

inline double valueFromTxArray(
	const double* logMu, int periodCount, int ageCount,
	int periodIndex, int ageIndex)
{
	periodIndex = std::clamp(periodIndex, 0, periodCount - 1);
	ageIndex = std::clamp(ageIndex, 0, ageCount - 1);

	// The caller has checked that periodCount * ageCount <= 2^31 - 1, so this cannot overflow.
	return logMu[periodIndex * ageCount + ageIndex];
}

// Parameters are `int`s because R matrix dimensions are `int`s and a `datey` is an `int` under the
// covers. `logMu` is in `tx` order: periods first, then ages.
inline double interpolatedLogMu(
	const double* logMu, int periodCount, int ageCount,
	int t0Clicks, int x0Clicks,
	int birthClicks, int tClicks)
{
	const int b0Clicks = t0Clicks - x0Clicks;

	// Cohort `cohort` is the older of the two either side of this birth, `cohort + 1` the younger.
	const WholeYearsAndClicks birth = splitClicks(birthClicks - b0Clicks);
	const int cohort = birth.years;
	const double cohortWeight = birth.clicks / static_cast<double>(ClicksPerYear);

	const WholeYearsAndClicks instant = splitClicks(tClicks - t0Clicks);
	const int period = instant.years;
	const double periodWeight = instant.clicks / static_cast<double>(ClicksPerYear);

	// Age is period minus cohort.
	const double olderNow = valueFromTxArray(logMu, periodCount, ageCount, period, period - cohort);
	const double youngerNow = valueFromTxArray(logMu, periodCount, ageCount, period, period - cohort - 1);
	const double olderNext = valueFromTxArray(logMu, periodCount, ageCount, period + 1, period + 1 - cohort);
	const double youngerNext = valueFromTxArray(logMu, periodCount, ageCount, period + 1, period - cohort);

	const double now = eve::fma(cohortWeight, youngerNow - olderNow, olderNow);
	const double next = eve::fma(cohortWeight, youngerNext - olderNext, olderNext);

	return eve::fma(periodWeight, next - now, now);
}

} // namespace
} // namespace mortality
