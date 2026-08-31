// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <optional>
#include <variant>
#include <vector>
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"
#include "veil/passPropagateIntervals.hpp" // Interval

namespace veil
{

// The interval-based comparison folding pass. It rewrites a comparison to a static bool whenever the
// two operands' propagated intervals settle the answer for every row and every time-point -- so
// `.i$age > 200` folds to false against a column that never exceeds 60, and so does `.x < 0`, which
// no pass reading a bare column could reach.
//
// It reads the intervals passPropagateIntervals worked out and takes them as an argument rather than
// computing them, keeping the analysis and the rewrite separately readable.
//
// This SUPERSEDES the earlier range-based fold, which compared a bare field's scanned range against a
// literal. A field's interval IS its scanned range and a literal's is a single point, so everything
// that pass folded is folded here, plus derived expressions and field-against-field.
//
// UNITS. Both operands must already measure in the same unit. Comparison narrowing is what puts them
// there, converting a plain number of years into a click threshold, and it now does so for any
// click-typed operand rather than only a bare field. The guard here is deliberate duplication: if
// narrowing ever declines -- an inexact equality, an overflowing threshold -- this refuses to fold
// rather than comparing clicks against years, which is how a datey comparison once folded to exactly
// the wrong answer.
//
// Folding does not re-run the interval analysis, so a fold that would let an enclosing expression
// narrow further is not chased in this pass. Running the two alternately to a fixed point is a
// possible refinement, deliberately not taken yet.

namespace detail
{

// Clicks and plain numbers are different units. Datey, durationy and category all count in their own
// integer space, so an operand of one of those types lines up only with the same type; everything in
// the plain double / bool family shares one unit.
inline bool comparableUnits(const TypeFull& a, const TypeFull& b) noexcept
{
    const auto isCounted = [](Type t) noexcept
    {
        return t == Type::Datey || t == Type::Durationy || t == Type::Category;
    };

    if (isCounted(a.type) || isCounted(b.type)) { return a.type == b.type; }
    return (a.type == Type::Double || a.type == Type::Bool)
        && (b.type == Type::Double || b.type == Type::Bool);
}

inline bool isPoint(const Interval& x) noexcept { return x.lo == x.hi && std::isfinite(x.lo); }

// Whether the two intervals cannot overlap at all, which is what settles an equality.
inline bool disjoint(const Interval& a, const Interval& b) noexcept
{
    return a.hi < b.lo || a.lo > b.hi;
}

// The settled answer of `a OP b` over every pair drawn from the two intervals, or empty when both
// outcomes remain possible. Each test asks whether the operator holds across the whole of both
// ranges, which is sound because the true values are contained in them.
inline std::optional<bool> intervalAnswer(Op op, const Interval& a, const Interval& b) noexcept
{
    switch (op)
    {
        case Op::Lt:
            if (a.hi < b.lo) { return true; }
            if (a.lo >= b.hi) { return false; }
            return std::nullopt;
        case Op::Le:
            if (a.hi <= b.lo) { return true; }
            if (a.lo > b.hi) { return false; }
            return std::nullopt;
        case Op::Gt:
            if (a.lo > b.hi) { return true; }
            if (a.hi <= b.lo) { return false; }
            return std::nullopt;
        case Op::Ge:
            if (a.lo >= b.hi) { return true; }
            if (a.hi < b.lo) { return false; }
            return std::nullopt;
        case Op::Eq:
            if (disjoint(a, b)) { return false; }
            // Only two intervals pinned to the same single value can force equality.
            if (isPoint(a) && isPoint(b) && a.lo == b.lo) { return true; }
            return std::nullopt;
        case Op::Ne:
            if (disjoint(a, b)) { return true; }
            if (isPoint(a) && isPoint(b) && a.lo == b.lo) { return false; }
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

} // namespace detail

inline void passFoldIntervalComparisons(Tree& tree, const std::vector<Interval>& intervals)
{
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        Node& node = tree.at(id);
        if (!isCall(node)) { continue; }

        const CallPayload& call = std::get<CallPayload>(node.payload);
        if (!isComparisonOp(call.op)) { continue; }

        const NodeId leftId = call.args[0];
        const NodeId rightId = call.args[1];
        if (leftId >= intervals.size() || rightId >= intervals.size()) { continue; }

        const std::optional<TypeFull>& leftType = tree.at(leftId).type;
        const std::optional<TypeFull>& rightType = tree.at(rightId).type;
        if (!leftType.has_value() || !rightType.has_value()) { continue; }
        if (!detail::comparableUnits(*leftType, *rightType)) { continue; }

        const std::optional<bool> answer =
            detail::intervalAnswer(call.op, intervals[leftId], intervals[rightId]);
        if (answer.has_value())
        {
            node.payload.emplace<LitPayload>(LitPayload{LiteralValue(*answer), TypeFull::createBool()});
        }
    }
}

} // namespace veil
