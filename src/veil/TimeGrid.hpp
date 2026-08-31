// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>
#include "datey.h" // ClicksPerYear, yearsFromClicks, isValidDatey

namespace veil
{

// The time vector for one individual: where the integration samples their exposure.
//
// The integral is a rectangular midpoint rule over uniform intervals of length dt, with a final
// interval of dt or less to reach the end of the exposure. For an exposure running from nu to tau,
// the sample points are
//
//     nu + dt/2,  nu + 3dt/2,  ...,  (nu + n*dt + tau) / 2
//
// where n is the number of whole intervals that fit. The last of those is the midpoint of the short
// final interval. If the individual died, tau itself is appended as one further slot, because the
// value wanted at the moment of death usually travels through the same calculation as the integrand
// and there is no sense computing it twice. The integral does not read that slot.
//
// EVERYTHING HERE IS IN CLICKS, and the arithmetic is exact integer arithmetic. That is the whole
// point of holding it this way: the sample points of two individuals with the same exposure are
// bit-identical, and a slot never drifts off the annual grid. Conversion to years happens once, at
// the point where a slot is read into a double vector.

// The integration intervals the framework allows, as a whole number of clicks. Each divides
// ClicksPerYear exactly, and each is itself even -- which is what lets the midpoint of a full
// interval land on a whole click, with no rounding anywhere.
//
// THAT GUARANTEE IS WHY ONLY THE FINAL INTERVAL NEEDS A ROUNDING RULE (Tim, 2026-07-28). It is a
// property of ClicksPerYear's particular value rather than of the arithmetic, so it is asserted
// here: change ClicksPerYear, or add an interval that leaves an odd number of clicks, and the build
// fails rather than the midpoints quietly starting to round.
//
// WHY THESE FOUR AND NO OTHERS (Tim, 2026-08-06). ClicksPerYear factorises as 2^3 x 3 x 5 x 61 x 73.
// Requiring half an interval to be a whole number of clicks means the denominator must divide
// ClicksPerYear / 2 = 2^2 x 3 x 5 x 61 x 73. The primes 61 and 73 are in ClicksPerYear so that a day
// is a whole number of clicks in both a common and a leap year (534360/365 = 1464 and
// 534360/366 = 1460, and see ClicksPerDay365/366 in datey.h). They are no use as an integration
// interval, because A DAY IS NOT A FIXED DURATION -- 1/365 of a year is not a day in a leap year --
// so the fine end is closed rather than cut off. Dropping them leaves the divisors of 60, of which
// 1, 4, 12 and 60 are the intuitive chain: a year, a quarter, a month, and a fifth of a month.
//
// HALF A YEAR IS DELIBERATELY ABSENT, for two reasons. Cost is linear in the number of intervals
// while the midpoint rule's error is quadratic, so a quarter costs twice what a half costs and buys
// four times the accuracy; there is no regime where a half is the right pick. And it would be the
// only entry breaking the refinement chain: midpoint samples sit at (2k+1)/2n of a year, so a
// coarse grid's samples survive in a finer one exactly when the RATIO OF DENOMINATORS IS ODD.
// 4 -> 12 and 12 -> 60 are odd, so refining the time scale keeps every sample already taken and adds
// more between them. Only the step off the year is even. Adding 2 would put an even ratio inside the
// fine end of the set.
constexpr int timeScaleClickOptions[] = {
    ClicksPerYear / 1,  // a year
    ClicksPerYear / 4,  // a quarter
    ClicksPerYear / 12, // a month
    ClicksPerYear / 60, // a fifth of a month
};

constexpr bool everyTimeScaleIsWholeClicks() noexcept
{
    for (const int clicks : timeScaleClickOptions)
    {
        if (clicks <= 0) { return false; }
        if (ClicksPerYear % clicks != 0) { return false; }
        if (clicks % 2 != 0) { return false; } // The half-interval too.
    }
    return true;
}

static_assert(everyTimeScaleIsWholeClicks(),
              "An integration interval, or half of one, is not a whole number of clicks.");

// The default is a quarter of a year: short enough that a mortality table's annual cells are
// sampled sensibly, long enough that a three to five year exposure needs only twelve to twenty
// slots.
constexpr int defaultTimeScaleClicks = ClicksPerYear / 4;

// THE TIME SCALE CROSSES THE BOUNDARY AS CLICKS, not as a count of intervals per year, because
// clicks are what ExposureColumns::deltaTClicks wants -- so this validates rather than converts.
//
// THE GUARD BELONGS HERE AND NOT ONLY IN R. The R entry points are called positionally, so a call
// site left behind by the rename would otherwise pass a bare `4` and integrate on a four-click step
// -- about seven minutes -- returning an answer that is wrong but entirely plausible. Refusing
// anything outside the permitted set turns that into a loud failure.
inline int validateTimeScaleClicks(int clicks)
{
    for (const int option : timeScaleClickOptions)
    {
        if (option == clicks) { return clicks; }
    }
    throw std::runtime_error(
        "veil: `time_scale` must be 1, 1/4, 1/12 or 1/60 of a year.");
}

// Half of a click count, rounded half to even.
//
// USED FOR THE FINAL INTERVAL'S MIDPOINT AND NOWHERE ELSE. Every full interval has an exact
// midpoint, because half of dt is itself a whole number of clicks for all four permitted intervals
// -- see the assertion above. The final, short interval is the one case where the midpoint,
// half of `nu + n*dt + tau`, can fall between two clicks.
//
// RULED BY TIM, 2026-07-28, after the three options were costed:
//   - a click grid with banker's halving, which is this;
//   - a half-click grid, which would make every sample exact and need no rounding at all, REJECTED
//     because a unit nobody expects surprises the reader more than the rounding costs;
//   - a plain `x >> 1`, REJECTED because it always shifts the sample the same way and so biases
//     systematically, where banker's lets the errors cancel across a portfolio.
// The residual is half a click -- about 29 seconds -- on one interval per individual, whatever the
// exposure length or dt. It is around three thousand times smaller than the midpoint rule's own
// discretisation error, which no choice here affects.
//
// NB the `halve` snippet in the specification sits under "ARCHIVED -- IGNORE FOR NOW", so it is not
// the authority for any of this; the ruling above is.
inline int64_t halveClicks(int64_t value) noexcept
{
    return (value + ((value >> 1) & 1)) >> 1;
}

// One individual's time vector, described rather than materialised.
//
// THERE IS NO ARRAY OF SAMPLE POINTS, and there should not be. All but two of the slots lie on an
// arithmetic progression -- `firstSlotClicks + slot * deltaTClicks` -- so a consumer generates them
// with one add per slot instead of one load per slot. The two that do not are the midpoint of the
// short final interval and, for a death, tau itself, and both are single values that live here.
//
// That also suits the SIMD kernel this becomes: a constant stride over the full slots needs no memory
// traffic at all, and the two odd slots peel off outside the loop, which is the shape one would want
// regardless. Reproducibility is untouched by the choice, because the arithmetic is integer --
// `firstSlotClicks + slot * deltaTClicks` is the same value however it is computed, so a vectorised
// form with per-lane offsets is bit-identical to the scalar loop. That is a reason to keep the grid
// in clicks rather than years.
struct TimeGrid final
{
    // Slots the integral reads. The first `fullSlots` carry the full dt; anything beyond that is
    // the single final slot, carrying `finalWidthYears`.
    int integrationSlots = 0;
    int fullSlots = 0;
    double finalWidthYears = 0.0;

    // The progression the full slots lie on.
    int64_t firstSlotClicks = 0;
    int64_t deltaTClicks = 0;

    // The midpoint of the short final interval, read only when `integrationSlots > fullSlots`.
    int64_t finalSlotClicks = 0;

    // Set when the individual died, in which case the last slot of all is tau.
    bool hasDeathSlot = false;
    int64_t deathSlotClicks = 0;

    int slotCount() const noexcept { return this->integrationSlots + (this->hasDeathSlot ? 1 : 0); }

    // The sample point at one slot. Written for a consumer that walks slots in order, which every
    // consumer does; the two comparisons fall away in a kernel that peels the tail off instead.
    int64_t clicksAt(int slot) const noexcept
    {
        if (slot < this->fullSlots)
        {
            return this->firstSlotClicks + static_cast<int64_t>(slot) * this->deltaTClicks;
        }
        return slot < this->integrationSlots ? this->finalSlotClicks : this->deathSlotClicks;
    }
};

// Describes the sample points for one exposure. Nothing is allocated: the grid is a handful of
// scalars, so there is no buffer to reuse between individuals and none to size.
inline TimeGrid buildTimeGrid(int64_t startClicks, int64_t endClicks, bool died, int64_t deltaTClicks)
{
    // MISSING BOUNDS ARE CAUGHT BEFORE THE ARITHMETIC, and this is not a formality. R's NA integer
    // is the most negative int, so an NA start against a real end passes the ordering test below and
    // describes an exposure running back four millennia -- tens of thousands of slots, and a
    // contribution that is enormous rather than absent. That is the one failure that survives
    // eyeballing the output, so it is refused by name here.
    //
    // An individual the include empties never reaches this: the interpreter builds an empty grid
    // instead. So this fires on data that is genuinely unusable, not on data that is merely excluded.
    if (!isValidDatey(startClicks) || !isValidDatey(endClicks))
    {
        throw std::runtime_error("veil: an exposure bound is missing, or outside the calendar.");
    }
    if (endClicks <= startClicks)
    {
        // An empty exposure has no integral to take and no death to value. The caller is meant to
        // skip such an individual entirely, so reaching here is a fault rather than a data question.
        throw std::runtime_error("veil: an exposure must end after it starts.");
    }
    if (deltaTClicks <= 0)
    {
        throw std::runtime_error("veil: the integration interval must be a positive number of clicks.");
    }

    const int64_t span = endClicks - startClicks;
    const int64_t wholeIntervals = span / deltaTClicks;
    const int64_t remainder = span - wholeIntervals * deltaTClicks;

    TimeGrid grid;
    grid.fullSlots = static_cast<int>(wholeIntervals);
    grid.hasDeathSlot = died;
    grid.deltaTClicks = deltaTClicks;
    grid.firstSlotClicks = startClicks + deltaTClicks / 2;
    grid.deathSlotClicks = endClicks;

    // A final interval of zero length would contribute nothing, so it is left out rather than
    // carried as a slot weighted by zero. That is the common case whenever an exposure is a whole
    // number of intervals long, which a calendar-year experience period usually is.
    grid.integrationSlots = grid.fullSlots + (remainder > 0 ? 1 : 0);
    grid.finalWidthYears = remainder > 0 ? yearsFromClicks(remainder) : 0.0;

    if (remainder > 0)
    {
        const int64_t finalIntervalStart = startClicks + wholeIntervals * deltaTClicks;
        grid.finalSlotClicks = halveClicks(finalIntervalStart + endClicks);
    }

    return grid;
}

// The integral of a sampled integrand: dt times the sum of the full-width slots, plus the final
// slot weighted by however much of an interval was left. The death slot, if there is one, sits past
// `integrationSlots` and is not read.
//
// The sum runs in slot order and nothing reorders it, so the answer is a function of the data alone.
// That matters more once this runs under SIMD and across threads -- see the reproducible-reduction
// note, since addition is not associative and IEEE 754 does not help here.
inline double integrateSlots(const double* values, const TimeGrid& grid) noexcept
{
    double total = 0.0;
    for (int slot = 0; slot < grid.fullSlots; ++slot) { total += values[slot]; }
    total *= yearsFromClicks(grid.deltaTClicks);

    if (grid.integrationSlots > grid.fullSlots)
    {
        total += values[grid.fullSlots] * grid.finalWidthYears;
    }
    return total;
}

} // namespace veil
