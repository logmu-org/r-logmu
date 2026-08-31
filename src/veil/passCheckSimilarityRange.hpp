// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
#include "veil/AevRecipe.hpp"
#include "veil/Node.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/passPropagateIntervals.hpp"

namespace veil
{

// The similarity range check. A similarity s is a proportion, so it belongs in [0, 1], and a
// distance d = -log s therefore cannot be negative. Neither bound can be ENFORCED: both are
// expressions of the record and of time, so their values are not known until the walk, and checking
// each evaluated value would put a comparison in the innermost loop to police something the user
// has already been told about.
//
// WHAT CAN BE DONE IS ANALYTIC, AND THIS IS IT. Interval propagation already works out the range
// each node can take, from literals and scanned column ranges upwards, so a violation that is
// certain BEFORE ANY RECORD IS READ can be refused at compile time. `similarity = 2`,
// `distance = -1`, and a similarity that is a column the scan says runs from 5 to 200, are all
// settled here for nothing.
//
// IT REFUSES ONLY WHAT IS ALWAYS WRONG -- an interval lying WHOLLY outside the bound -- and that
// restraint is deliberate rather than timid. Intervals are conservative: they widen through
// arithmetic and are unknown wherever a rule is missing. Refusing anything that MIGHT violate would
// reject `distance = (2025 - .t) / 10`, the ordinary decay kernel, because `.t` carries the whole
// representable range of dates rather than the span of the data, so the interval permits a negative
// distance that the exposure never reaches. A check that fires on correct code is worse than one
// that stays quiet on incorrect code, because only the first makes a user distrust it.
//
// UNITS. Intervals are in their own node's units -- clicks for a datey or durationy, plain numbers
// for a double -- so this compares nothing until it has established that the node is a double. A
// non-double here would already have failed typing, since the recipe multiplies it by the weight,
// which is why this is an assertion rather than a diagnostic.

// Which bound was broken, so the caller can say it in the user's own words rather than in terms of
// the factor the engine built.
enum class SimilarityViolation : unsigned char
{
    None,
    Above,  // a similarity above 1, i.e. a negative distance
    Below,  // a similarity below 0, which a distance cannot produce
};

inline SimilarityViolation passCheckSimilarityRange(
    const Tree& tree,
    const std::vector<Interval>& intervals,
    NodeId similarity,
    SimilarityForm form)
{
    if (similarity == invalidNodeId) { return SimilarityViolation::None; }
    if (similarity >= static_cast<NodeId>(intervals.size())) { return SimilarityViolation::None; }

    const Node& node = tree.at(similarity);
    if (!node.type.has_value() || node.type->type != Type::Double)
    {
        throw std::runtime_error("veil: a similarity or distance must be a number.");
    }

    const Interval& range = intervals[similarity];

    // A distance is the negative log, so the bound reads the other way round: d < 0 is s > 1, and no
    // distance can produce a similarity below zero at all, since exp is positive everywhere.
    if (form == SimilarityForm::Distance)
    {
        return range.hi < 0.0 ? SimilarityViolation::Above : SimilarityViolation::None;
    }

    if (range.hi < 0.0) { return SimilarityViolation::Below; }
    if (range.lo > 1.0) { return SimilarityViolation::Above; }
    return SimilarityViolation::None;
}

} // namespace veil
