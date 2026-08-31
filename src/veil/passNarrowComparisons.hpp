// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cmath>
#include <limits>
#include <optional>
#include "datey.h" // ClicksPerYear
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// The comparison-narrowing pass. By the datey package's own operator rules (see the
// datey-ops-coupling reference, https://r-datey.logmu.org/reference/ops.html), a datey/durationy
// value compared against a plain number is evaluated by converting it to YEARS and comparing there
// -- `.i$birth > 1990` means "birth's year exceeds 1990", not "birth's click count exceeds 1990". A
// datey/durationy is a whole number of clicks, so that comparison has an integer form: convert the
// number to clicks once, round it in the direction the operator wants (`>`/`>=` up, `<`/`<=` down),
// and compare click against click.
//
// THIS IS THE ONE PLACE YEARS BECOME CLICKS. Everything downstream then compares like with like,
// which is why the folding passes need no unit arithmetic of their own -- only a guard that refuses
// a comparison this pass did not normalise. It applies to ANY click-typed operand, not just a bare
// field, so a derived expression is covered too: `.x` is `.t - .b`, which the type pass types as a
// durationy, so `.x < 65` narrows exactly as `.i$birth < 1990` does. Coercion deliberately leaves
// comparison operands alone, so a click-typed operand is still visible here rather than buried
// under a `to_double`.
//
// EXACTNESS: EVERY comparison is narrowed only when the number lands on a whole click. Given that,
// the rewrite is exact for all six operators, and the argument is short. Clicks to years is
// strictly increasing -- consecutive clicks are about 1.9e-6 years apart, millions of times the
// spacing of doubles at these magnitudes -- so if some click c0 converts to exactly the number,
// then `t < y` holds for exactly the clicks below c0, `t == y` for exactly c0, and so on. No
// rounding enters any of those.
//
// Without that precondition a threshold still exists mathematically, but finding it means
// reasoning about how the scaled product rounded and how the conversion rounds either side of the
// boundary, and the year-space and click-space answers can then disagree about the single click on
// the boundary. Declining costs almost nothing: whole years, ages, months and quarters are all
// whole clicks, which is what expressions are written in.
//
// Data-independent: unlike the folding passes this needs no column scan, only the node types set by
// type annotation, so it runs on type/literal shape alone.
//
// Left for a later pass: narrowing a plain (non-click) double operand by its `allIntegral` scan
// flag -- that is data-dependent, unlike this pass's pure type-driven rewrite.

namespace detail
{

inline bool isClickType(Type t) noexcept { return t == Type::Datey || t == Type::Durationy; }

inline std::optional<int> toIntChecked(double v) noexcept
{
    if (v < static_cast<double>(std::numeric_limits<int>::min())
        || v > static_cast<double>(std::numeric_limits<int>::max()))
    {
        return std::nullopt;
    }
    return static_cast<int>(v);
}

// What `clicks OP years` narrows to: the operator to compare with, and the click threshold to
// compare against. Empty when the rewrite is declined, which happens when the number is not a whole
// click and when a threshold would overflow `int`. The comparison is then left exactly as it was,
// and the unit guard in passFoldIntervalComparisons refuses to fold it rather than mixing years
// with clicks.
struct Narrowed final
{
    Op op;
    int thresholdClicks;
};

// Whether the number is exactly the value of some click, and if so which. Asked by converting it to
// clicks and straight back: the two conversions round-trip every valid click exactly, so a number
// survives the journey unchanged precisely when a click denotes it. This also finds the right
// candidate, since a number denoted by c0 scales to within a millionth of c0 and so rounds to it.
inline std::optional<int> wholeClick(double years)
{
    const std::optional<int> candidate = toIntChecked(clicksFromYears(years));
    if (!candidate.has_value()) { return std::nullopt; }
    if (yearsFromClicks(*candidate) != years) { return std::nullopt; }
    return candidate;
}

inline std::optional<Narrowed> narrow(Op op, double years)
{
    // Every operator needs the number to land on a whole click; see the exactness note above.
    const std::optional<int> c0 = wholeClick(years);
    if (!c0.has_value()) { return std::nullopt; }

    // With `t` a click and `y` the click c0, each comparison has an exact integer form. The two
    // strict operators shift by a click, which is why they need their own range check.
    switch (op)
    {
        case Op::Eq: return Narrowed{Op::Eq, *c0};
        case Op::Ne: return Narrowed{Op::Ne, *c0};
        case Op::Le: return Narrowed{Op::Le, *c0}; // t <= y  <=>  t <= c0
        case Op::Ge: return Narrowed{Op::Ge, *c0}; // t >= y  <=>  t >= c0
        case Op::Lt: // t <  y  <=>  t <= c0 - 1
        {
            const std::optional<int> below = toIntChecked(static_cast<double>(*c0) - 1.0);
            if (!below.has_value()) { return std::nullopt; }
            return Narrowed{Op::Le, *below};
        }
        case Op::Gt: // t >  y  <=>  t >= c0 + 1
        {
            const std::optional<int> above = toIntChecked(static_cast<double>(*c0) + 1.0);
            if (!above.has_value()) { return std::nullopt; }
            return Narrowed{Op::Ge, *above};
        }
        default:
            return std::nullopt;
    }
}

} // namespace detail

inline void passNarrowComparisons(Tree& tree)
{
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        // Everything this iteration needs is copied out into plain values before anything is added
        // to the tree. Adding a node can reallocate the arena, so no reference into it may be held
        // across the call that builds the threshold literal.
        bool leftIsClicks = false;
        std::optional<TypeFull> clickType; // TypeFull holds const members, so it is emplaced, not assigned.
        Op op = Op::Eq;
        double years = 0.0;
        {
            const Node& node = tree.at(id);
            if (!isCall(node)) { continue; }

            const CallPayload& call = std::get<CallPayload>(node.payload);
            if (!isComparisonOp(call.op)) { continue; }

            const Node& left = tree.at(call.args[0]);
            const Node& right = tree.at(call.args[1]);

            // Any click-typed operand qualifies, whether it is a field, a literal or a derived
            // expression such as `.t - .b`. Exactly one side must be clicks: if both are, they
            // already share a unit and there is nothing to convert.
            leftIsClicks = left.type.has_value() && detail::isClickType(left.type->type);
            const bool rightIsClicks = right.type.has_value() && detail::isClickType(right.type->type);
            if (leftIsClicks == rightIsClicks) { continue; }

            const Node& litNode = leftIsClicks ? right : left;
            if (!isLit(litNode)) { continue; }

            const LitPayload& lit = std::get<LitPayload>(litNode.payload);
            if (lit.type.type != Type::Double) { continue; } // Already a click literal -- nothing to narrow.
            const double* value = std::get_if<double>(&lit.value);
            if (value == nullptr) { continue; }

            clickType.emplace(*(leftIsClicks ? left : right).type);
            op = leftIsClicks ? call.op : reverseComparison(call.op);
            years = *value;
        }

        const std::optional<detail::Narrowed> narrowed = detail::narrow(op, years);
        if (!narrowed.has_value()) { continue; } // Declined: leave the comparison exactly as it was.

        // Build a NEW literal holding the click threshold and point this comparison at it, rather
        // than rewriting the existing one where it lies. Common sub-expression elimination makes the
        // tree a DAG, and a shared literal may be read as a plain number elsewhere; converting it in
        // place would silently change its value for every other reader. Appending is safe mid-walk
        // because the new node is a literal, which the loop skips.
        const NodeId thresholdId = tree.buildLitInt(narrowed->thresholdClicks, *clickType);

        // Re-fetch after the append: the reference taken above may no longer be valid.
        CallPayload& call = std::get<CallPayload>(tree.at(id).payload);
        call.args[leftIsClicks ? 1 : 0] = thresholdId;
        // Reverse back to the original argument order when the literal sits on the left, so the
        // operator still reads correctly against the operand positions.
        call.op = leftIsClicks ? narrowed->op : reverseComparison(narrowed->op);
    }
}

} // namespace veil
