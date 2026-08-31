// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>
#include "datey.h" // isValidDatey
#include "mortality/lookup_log_mu.hpp"

namespace veil
{

// An age-period table of log mu. The values are held in the source matrix's own column-major order,
// so `logMu[age + period * ageCount]`, with `age` counted from `x0` and `period` from `t0`.
//
// `x0` is a durationy (the youngest age) and `t0` a datey (the earliest period), both in clicks, and
// both denote the MIDDLE of their year rather than its start.
struct MortalityTable final
{
    std::vector<double> logMu;
    uint32_t ageCount = 0;
    uint32_t periodCount = 0;
    int x0Clicks = 0;
    int t0Clicks = 0;

    bool isWellFormed() const noexcept
    {
        return this->logMu.size() == static_cast<size_t>(this->ageCount) * this->periodCount
            && this->ageCount > 0
            && this->periodCount > 0;
    }
};

// A table one block reads, identified by its position in that block's own list. Block-local for the
// same reason an OperandId is: a block is the unit that gets handed to a worker.
using TableId = uint32_t;

constexpr TableId invalidTableId = static_cast<TableId>(-1);

// log mu for one individual at one instant.
//
// THE SEMANTICS ARE NOT VEIL'S TO CHOOSE. `log_mu()` on a `mortality_table` in R already answers
// this question, by way of `mortality::get_interpolated_log_mu_bt_from_annual_tx_array`, and veil
// has to give the same answer for the same table -- so it calls that same function rather than
// working the lookup out a second time. Two properties of it matter here, because both are easy to
// assume wrongly:
//
//   - IT INTERPOLATES, bilinearly, and on the COHORT-PERIOD lattice rather than the age-period one.
//     So log mu varies continuously along an individual's exposure and is NOT constant within an
//     annual cell. Anything that wants to exploit annual cells has to reckon with that first.
//   - AN INDEX PAST THE END IS CLAMPED to the edge of the table, age and period independently. An
//     individual older than the table's last age is read at that last age rather than refused,
//     which is consistent with `end_age` being a check the caller makes rather than an input here.
//
// ARGUMENT ORDER IS THE TRAP. The shared function counts in `tx` order -- periods first, then ages
// -- where the table's own fields read the other way round. This wrapper is the one place that
// translation happens, so there is only one line to get wrong rather than one per call site.
// CLICKS ARE 64-BIT INSIDE THE ENGINE AND THE SHARED LOOKUP IS `int`-TYPED, so the narrowing happens
// here and is stated rather than left implicit. It is guarded rather than merely cast because of how
// a bad value would surface: an out-of-calendar date puts the lattice index far outside the table,
// the lookup CLAMPS to the edge, and what comes back is a perfectly plausible log mu for the wrong
// individual. A missing birth date is the way that actually happens, R's NA integer being the most
// negative int. The two comparisons cost nothing against the divmod, four loads and three lerps below.
inline double logMuAt(const MortalityTable& table, int64_t birthClicks, int64_t tClicks)
{
    if (!isValidDatey(birthClicks))
    {
        throw std::runtime_error("veil: a birth date is missing, or outside the calendar, so log mu "
                                 "cannot be looked up.");
    }
    if (!isValidDatey(tClicks))
    {
        throw std::runtime_error("veil: a log mu lookup was asked for an instant outside the calendar.");
    }

    return mortality::get_interpolated_log_mu_bt_from_annual_tx_array(
        table.logMu.data(),
        static_cast<int>(table.periodCount), static_cast<int>(table.ageCount),
        table.t0Clicks, table.x0Clicks,
        static_cast<int>(birthClicks), static_cast<int>(tClicks));
}

} // namespace veil
