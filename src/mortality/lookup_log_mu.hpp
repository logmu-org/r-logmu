// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

// THE DEFINITION OF A TABLE'S LOG MU AT ONE BIRTH AND ONE INSTANT, AND ALONG ONE LIFELINE.
//
// ONE SOURCE FOR EVERY PATH. R's `log_mu()` and the veil engine's bulk lookup both compile this
// header, so they run the same algorithm. They are compiled in different translation units, some
// with FMA and some without, so they may differ by an ulp or two, which is the agreed tolerance
// between paths. A second, hand-written copy of the lookup anywhere else would break the one-source
// rule even if it agreed today.
//
// THE TWO FORMS ARE BUILT FROM THE SAME TWO PIECES, a cohort blend at one period and an
// interpolation in period between two blends. `interpolatedLogMu` does both for one instant;
// `interpolatedLogMuSeries` blends once per period an individual spans and then does one
// interpolation per instant. Compiled with the same flags, the two give bit-identical answers,
// because the series form only keeps blends the single form would recompute.
//
// AN ANONYMOUS NAMESPACE, DELIBERATELY. The bulk lookup is compiled into each SIMD tier's
// translation unit with that tier's ISA flags, so an ordinary `inline` function here would be one
// symbol with several different bodies, and the linker keeps whichever it sees first. A baseline CPU
// could then run the AVX-512 body. Internal linkage gives each translation unit its own copy.
//
// COHORT FIRST, THEN PERIOD, EACH STEP AN `eve::fma`. The lattice is cohort-period, so a person's
// whole exposure lies between two cohort columns fixed for life: blending those first leaves one
// interpolation in period per instant, which is the shape the bulk path needs. `eve::fma` is a fused
// multiply-add where the target has the instruction and a separate multiply and add where it does
// not. `std::lerp` is not used: its guarantees at t = 1 and on monotonicity buy nothing for a weight
// in [0, 1), its branches get in the way of lanes, and the standard fixes only its properties, so
// libraries may round it differently. `std::fma` is not used either: compiled without `-mfma`, Rtools
// turns it into a call to a library `fma` that is not correctly rounded.
//
// AN INDEX PAST THE END IS CLAMPED to the edge of the table, age and period independently. An
// individual older than the table's last age is read at that last age rather than refused, which is
// consistent with `end_age` being a check the caller makes rather than an input here.

#include <algorithm>
#include <cstddef>
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

// Where a birth sits between the two cohort columns either side of it. `cohort` is the older of the
// two, `cohort + 1` the younger, and `weight` runs from 0 at the older to 1 at the younger.
struct CohortPosition
{
	int cohort;
	double weight;
};

inline CohortPosition cohortPositionOf(int t0Clicks, int x0Clicks, int birthClicks)
{
	const WholeYearsAndClicks birth = splitClicks(birthClicks - (t0Clicks - x0Clicks));
	return {birth.years, birth.clicks / static_cast<double>(ClicksPerYear)};
}

// The two cohort columns blended at one whole period. Age is period minus cohort.
inline double cohortBlendAt(
	const double* logMu, int periodCount, int ageCount,
	const CohortPosition& position, int period)
{
	const double older = valueFromTxArray(logMu, periodCount, ageCount, period, period - position.cohort);
	const double younger = valueFromTxArray(logMu, periodCount, ageCount, period, period - position.cohort - 1);
	return eve::fma(position.weight, younger - older, older);
}

inline double periodWeightOf(const WholeYearsAndClicks& instant)
{
	return instant.clicks / static_cast<double>(ClicksPerYear);
}

inline double periodBlend(double now, double next, double periodWeight)
{
	return eve::fma(periodWeight, next - now, now);
}

// Parameters are `int`s because R matrix dimensions are `int`s and a `datey` is an `int` under the
// covers. `logMu` is in `tx` order: periods first, then ages.
inline double interpolatedLogMu(
	const double* logMu, int periodCount, int ageCount,
	int t0Clicks, int x0Clicks,
	int birthClicks, int tClicks)
{
	const CohortPosition position = cohortPositionOf(t0Clicks, x0Clicks, birthClicks);
	const WholeYearsAndClicks instant = splitClicks(tClicks - t0Clicks);

	return periodBlend(
		cohortBlendAt(logMu, periodCount, ageCount, position, instant.years),
		cohortBlendAt(logMu, periodCount, ageCount, position, instant.years + 1),
		periodWeightOf(instant));
}

// How many cohort blends the series form keeps at once. A five-year investigation spans six or seven
// periods; a longer one, or instants out of order, refills rather than failing.
constexpr int BlendsHeld = 32;

// log mu for one birth at `count` instants, written to `result`.
//
// THE COHORT BLEND IS DONE ONCE PER PERIOD rather than twice per instant, which is the whole point:
// a birth fixes both cohort columns for life, so along one lifeline each instant needs only its
// period's two blends and one interpolation between them. The blends are kept for a run of periods
// starting at the first instant's, sized from the last instant's on the assumption that instants
// increase, which every caller's do.
//
// THE ASSUMPTION AFFECTS SPEED, NEVER THE ANSWER. An instant whose two periods are not both held
// refills from that instant's period, so instants out of order, or a lifeline longer than
// `BlendsHeld` periods, come out exactly as they would one at a time.
inline void interpolatedLogMuSeries(
	const double* logMu, int periodCount, int ageCount,
	int t0Clicks, int x0Clicks,
	int birthClicks, std::ptrdiff_t count, const int* tClicks, double* result)
{
	if (count <= 0)
	{
		return;
	}

	const CohortPosition position = cohortPositionOf(t0Clicks, x0Clicks, birthClicks);
	const int lastPeriod = splitClicks(tClicks[count - 1] - t0Clicks).years;

	double blends[BlendsHeld];
	int heldFrom = 0;
	int heldCount = 0;

	for (std::ptrdiff_t index = 0; index < count; ++index)
	{
		const WholeYearsAndClicks instant = splitClicks(tClicks[index] - t0Clicks);
		const int period = instant.years;

		// Both `period` and `period + 1` must be held.
		if (period < heldFrom || period + 1 >= heldFrom + heldCount)
		{
			heldFrom = period;
			heldCount = std::clamp(lastPeriod - period + 2, 2, BlendsHeld);
			for (int offset = 0; offset < heldCount; ++offset)
			{
				blends[offset] = cohortBlendAt(logMu, periodCount, ageCount, position, heldFrom + offset);
			}
		}

		const int offset = period - heldFrom;
		result[index] = periodBlend(blends[offset], blends[offset + 1], periodWeightOf(instant));
	}
}

} // namespace
} // namespace mortality
