// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <variant>
#include <vector>
#include "datey.h" // ClicksPerYear, ValidDateStartClicks, ValidDateEndClicks
#include "veil/ColumnScan.hpp" // TypeWithConstraints
#include "veil/Node.hpp"
#include "veil/ObjStore.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"
#include "veil/TypeSpecificConstraint.hpp"

namespace veil
{

// The interval propagation pass. It works out, for every node, the range of values that node can
// take, from the literals and scanned column ranges at the leaves upwards through the operators. It
// ANNOTATES ONLY -- nothing is rewritten here. The folding that consumes these intervals is a
// separate pass, which is what lets this one be read and tested as a single question: what can this
// expression evaluate to?
//
// It is what lets a later fold settle a comparison over a DERIVED expression, which range folding
// cannot: `.x < 65` is a subtraction, not a bare field, so only a propagated interval can decide it.
//
// UNITS -- READ THIS BEFORE CONSUMING AN INTERVAL. Every interval is expressed in its own node's
// units, matching the type annotation pass: a datey or durationy node's interval counts CLICKS, a
// double node's counts plain numbers. `to_double` is the one operator that rescales, dividing by
// ClicksPerYear. A node's `type` is what says which unit its interval is in, so a consumer must
// check that type before comparing an interval against anything.
//
// This matters because the obvious next use is exactly where it bites. `.x` is `.t - .b`, a
// subtraction of two dateys, so its interval is a durationy in CLICKS -- while the `65` in
// `.x < 65` is a plain number of YEARS. Comparing the two directly is off by a factor of 534360,
// which is the very mistake that made range folding answer `.i$birth > 1990` with true. Comparison
// narrowing rescales the literal for a bare field; a fold over derived expressions will have to do
// the same for itself.
//
// TRAVERSAL. Node ids do NOT reliably run children-before-parents by the time this pass runs: the
// coercion pass appends its `to_double` nodes at the END of the arena and repoints parents at them,
// so a parent can hold an argument with a larger id than its own. So this walks post-order from the
// root and memoises, rather than sweeping ids forwards. Memoisation also means a node shared by
// several parents -- the DAG that common sub-expression elimination produces -- is costed once.
//
// An unknown value is the whole real line. That is the safe default, and every rule below either
// tightens it soundly or leaves it alone, so a rule that is missing or too hard simply yields no
// optimisation rather than a wrong one.

struct Interval final
{
    double lo;
    double hi;

    static constexpr Interval unknown() noexcept
    {
        return Interval{-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    }

    static constexpr Interval point(double value) noexcept { return Interval{value, value}; }

    static constexpr Interval bounds(double lo, double hi) noexcept { return Interval{lo, hi}; }

    // The unit interval a logical takes, read as 0 or 1.
    static constexpr Interval boolean() noexcept { return Interval{0.0, 1.0}; }

    bool isKnown() const noexcept { return std::isfinite(this->lo) || std::isfinite(this->hi); }
    bool isBounded() const noexcept { return std::isfinite(this->lo) && std::isfinite(this->hi); }
};

namespace detail
{

// Guards against the indeterminate forms interval arithmetic runs into at the edges -- an infinite
// bound multiplied by a zero one, most obviously. Any NaN endpoint collapses to unknown.
inline Interval sanitised(double lo, double hi) noexcept
{
    if (std::isnan(lo) || std::isnan(hi)) { return Interval::unknown(); }
    if (lo > hi) { return Interval::unknown(); } // Should not arise; treated as "no information".
    return Interval{lo, hi};
}

inline Interval negate(const Interval& x) noexcept { return sanitised(-x.hi, -x.lo); }

inline Interval absolute(const Interval& x) noexcept
{
    if (x.lo >= 0.0) { return x; }
    if (x.hi <= 0.0) { return negate(x); }
    return sanitised(0.0, std::max(-x.lo, x.hi)); // Spans zero.
}

inline Interval add(const Interval& a, const Interval& b) noexcept
{
    return sanitised(a.lo + b.lo, a.hi + b.hi);
}

inline Interval subtract(const Interval& a, const Interval& b) noexcept
{
    return sanitised(a.lo - b.hi, a.hi - b.lo);
}

inline Interval multiply(const Interval& a, const Interval& b) noexcept
{
    const double products[4] = {a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    double lo = products[0];
    double hi = products[0];
    for (const double product : products)
    {
        if (std::isnan(product)) { return Interval::unknown(); } // An infinity met a zero.
        lo = std::min(lo, product);
        hi = std::max(hi, product);
    }
    return sanitised(lo, hi);
}

// Division, which needs a guard none of the operators above do.
//
// A divisor SPANNING zero puts a gap in the result: [1, 2] / [-1, 1] is (-inf, -1] together with
// [1, +inf), which one interval cannot express. That much is obvious.
//
// A divisor merely TOUCHING zero from BELOW is the trap, because the endpoint arithmetic looks like
// it copes. Take [1, 2] / [-2, 0]: the four quotients are -0.5, +inf, -1 and +inf, so the extremes
// would say [-1, +inf] -- yet a divisor of -0.2 yields -5, and one of -0.001 yields -1000, neither
// of which that bound contains. Not merely loose: WRONG, in the one direction a bound must never be,
// since it omits values the engine can produce. The cause is that an interval endpoint written 0 is
// a POSITIVE zero, so dividing by it gives +inf where the values just inside the interval tend to
// -inf. An interval carries no sign for its endpoints, so there is nothing to test.
//
// Touching zero from ABOVE would in fact be safe -- a / +0 gives the correctly signed infinity
// there, so the corners are sound -- but the guard refuses it too. An asymmetric test (>= on one
// side, < on the other) buys a half-bounded interval in a rare case and invites exactly the kind of
// sign mistake this comment exists to warn about.
//
// With the divisor held strictly to one side of zero the quotient is monotonic in each argument
// across the whole box, so the extremes sit at its corners -- the same reasoning `multiply` rests on.
inline Interval divide(const Interval& a, const Interval& b) noexcept
{
    if (!(b.lo > 0.0 || b.hi < 0.0)) { return Interval::unknown(); }

    const double quotients[4] = {a.lo / b.lo, a.lo / b.hi, a.hi / b.lo, a.hi / b.hi};
    double lo = quotients[0];
    double hi = quotients[0];
    for (const double quotient : quotients)
    {
        if (std::isnan(quotient)) { return Interval::unknown(); } // An infinity divided by an infinity.
        lo = std::min(lo, quotient);
        hi = std::max(hi, quotient);
    }
    return sanitised(lo, hi);
}

inline Interval smallest(const Interval& a, const Interval& b) noexcept
{
    return sanitised(std::min(a.lo, b.lo), std::min(a.hi, b.hi));
}

inline Interval largest(const Interval& a, const Interval& b) noexcept
{
    return sanitised(std::max(a.lo, b.lo), std::max(a.hi, b.hi));
}

// The smallest interval containing both -- what a value that could come from either one lies in.
inline Interval spanning(const Interval& a, const Interval& b) noexcept
{
    return sanitised(std::min(a.lo, b.lo), std::max(a.hi, b.hi));
}

// clamp(x, low, high) never leaves the bounds it was given, whatever x does. Deliberately the loose
// reading -- it ignores what x contributes -- because that is sound for every x and needs no case
// analysis of the three intervals against each other.
inline Interval clamped(const Interval& low, const Interval& high) noexcept
{
    return sanitised(low.lo, high.hi);
}

// A monotonically increasing function maps endpoints to endpoints.
template <typename Function>
inline Interval increasing(const Interval& x, Function function) noexcept
{
    return sanitised(function(x.lo), function(x.hi));
}

// The interval a scanned column's values lie in, in the column's own units. A column with NAs gives
// nothing away: the same caution the folding passes take, so no rule has to define what an interval
// containing a missing value would mean.
inline Interval columnInterval(const TypeWithConstraints& constraints) noexcept
{
    if (constraints.hasNAs || !constraints.hasValues) { return Interval::unknown(); }
    if (const auto* d = std::get_if<DoubleConstraint>(&constraints.specific))
    {
        return sanitised(d->min, d->max);
    }
    if (const auto* n = std::get_if<IntConstraint>(&constraints.specific))
    {
        return sanitised(static_cast<double>(n->min), static_cast<double>(n->max));
    }
    if (const auto* b = std::get_if<BoolConstraint>(&constraints.specific))
    {
        return sanitised(b->min ? 1.0 : 0.0, b->max ? 1.0 : 0.0);
    }
    return Interval::unknown(); // Text carries no numeric range.
}

inline Interval literalInterval(const LitPayload& lit) noexcept
{
    if (const auto* b = std::get_if<bool>(&lit.value)) { return Interval::point(*b ? 1.0 : 0.0); }
    if (const auto* d = std::get_if<double>(&lit.value)) { return Interval::point(*d); }
    if (const auto* i = std::get_if<int>(&lit.value)) { return Interval::point(static_cast<double>(*i)); }
    return Interval::unknown(); // Text.
}

inline Interval objInterval(const ObjPayload& payload, const ObjStore& objs) noexcept
{
    const Obj& obj = objs.at(payload.obj);
    if (const auto* constant = std::get_if<MortalityConst>(&obj)) { return Interval::point(constant->logMu); }
    return Interval::unknown(); // A table's values and an include's interval are not bounded here.
}

} // namespace detail

// Computes an interval for every node reachable from the root. Nodes that are not reached keep the
// unknown interval. `columnConstraints` is indexed by ColumnId, as the folding passes take it.
// `timeInterval` is what `.t` can be, IN CLICKS. The default is the representable calendar, which is
// all a caller with no exposure knows -- the diagnostic entry points are in that position. A caller
// that has the data can do much better: every sample point of the integration lies within the
// individual's own exposure, and an include only ever narrows that, so across a whole dataset `.t`
// is confined to the first exposure start and the last exposure end. That is a span of years rather
// than of millennia, and it is what lets a range test over a derived expression settle.
inline std::vector<Interval> passPropagateIntervals(
    const Tree& tree,
    const std::vector<std::optional<TypeWithConstraints>>& columnConstraints,
    const ObjStore& objs,
    Interval timeInterval = Interval::bounds(static_cast<double>(ValidDateStartClicks),
                                             static_cast<double>(ValidDateEndClicks)))
{
    std::vector<Interval> intervals(tree.size(), Interval::unknown());
    std::vector<char> settled(tree.size(), 0);

    // Explicit post-order walk. Each node is pushed twice: once to have its children scheduled, and
    // once more, after them, to be computed.
    struct Frame final
    {
        NodeId id;
        bool childrenScheduled;
    };

    if (tree.size() == 0) { return intervals; }

    // ONE WALK OVER EVERY ROOT, sharing the settled map. A node reached from two roots is computed
    // once, which is the whole reason the outputs of a calculation live in one arena.
    // SEEDED FROM EVERY NODE, NOT JUST THE ROOTS. An include's gate is a real expression that has to
    // be analysed like any other, but it is not an output, so it is not a root and walking down from
    // the roots would step straight past it. Sweeping the arena reaches it, and reaches anything else
    // an earlier pass left behind; the post-order below still visits each node exactly once, so the
    // cost is the same walk in a different order.
    std::vector<Frame> stack;
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        stack.push_back(Frame{id, false});
    }

    while (!stack.empty())
    {
        const Frame frame = stack.back();
        stack.pop_back();

        if (settled[frame.id] != 0) { continue; }

        const Node& node = tree.at(frame.id);

        if (!frame.childrenScheduled && isCall(node))
        {
            stack.push_back(Frame{frame.id, true});
            for (const NodeId argId : std::get<CallPayload>(node.payload).args)
            {
                if (settled[argId] == 0) { stack.push_back(Frame{argId, false}); }
            }
            continue;
        }

        Interval result = Interval::unknown();

        if (isLit(node))
        {
            result = detail::literalInterval(std::get<LitPayload>(node.payload));
        }
        else if (isField(node))
        {
            const ColumnId column = std::get<FieldPayload>(node.payload).column;
            if (column < columnConstraints.size() && columnConstraints[column].has_value())
            {
                result = detail::columnInterval(*columnConstraints[column]);
            }
        }
        else if (isTime(node))
        {
            result = timeInterval;
        }
        else if (isObj(node))
        {
            result = detail::objInterval(std::get<ObjPayload>(node.payload), objs);
        }
        else
        {
            const CallPayload& call = std::get<CallPayload>(node.payload);
            const auto argument = [&](size_t position) { return intervals[call.args[position]]; };

            switch (call.op)
            {
                case Op::Pos: result = argument(0); break;
                case Op::Neg: result = detail::negate(argument(0)); break;
                case Op::Abs: result = detail::absolute(argument(0)); break;
                case Op::Add: result = detail::add(argument(0), argument(1)); break;
                case Op::Sub: result = detail::subtract(argument(0), argument(1)); break;
                case Op::Mul: result = detail::multiply(argument(0), argument(1)); break;
                case Op::RDiv: result = detail::divide(argument(0), argument(1)); break;
                case Op::Min: result = detail::smallest(argument(0), argument(1)); break;
                case Op::Max: result = detail::largest(argument(0), argument(1)); break;
                case Op::Clamp: result = detail::clamped(argument(1), argument(2)); break;

                // A selection yields one branch or the other, so it spans both. The condition, which
                // is argument 0, contributes no value of its own.
                case Op::Select:
                    result = detail::spanning(argument(1), argument(2));
                    break;

                case Op::ToDouble:
                {
                    // The one rescaling operator: clicks become years. A logical already reads as
                    // 0 or 1 and needs no conversion.
                    const Interval source = argument(0);
                    const std::optional<TypeFull>& sourceType = tree.at(call.args[0]).type;
                    if (sourceType.has_value()
                        && (sourceType->type == Type::Datey || sourceType->type == Type::Durationy))
                    {
                        const double perClick = 1.0 / static_cast<double>(ClicksPerYear);
                        result = detail::sanitised(source.lo * perClick, source.hi * perClick);
                    }
                    else if (sourceType.has_value() && sourceType->type == Type::Bool)
                    {
                        result = source;
                    }
                    break;
                }

                case Op::Exp:
                    result = detail::increasing(argument(0), [](double v) { return std::exp(v); });
                    break;
                case Op::Sqrt:
                {
                    const Interval x = argument(0);
                    if (x.lo >= 0.0) { result = detail::increasing(x, [](double v) { return std::sqrt(v); }); }
                    break;
                }
                case Op::Log:
                {
                    const Interval x = argument(0);
                    if (x.lo > 0.0) { result = detail::increasing(x, [](double v) { return std::log(v); }); }
                    break;
                }
                case Op::Log10:
                {
                    const Interval x = argument(0);
                    if (x.lo > 0.0) { result = detail::increasing(x, [](double v) { return std::log10(v); }); }
                    break;
                }
                case Op::Log1p:
                {
                    const Interval x = argument(0);
                    if (x.lo > -1.0) { result = detail::increasing(x, [](double v) { return std::log1p(v); }); }
                    break;
                }
                case Op::Expm1:
                    result = detail::increasing(argument(0), [](double v) { return std::expm1(v); });
                    break;
                case Op::Floor:
                    result = detail::increasing(argument(0), [](double v) { return std::floor(v); });
                    break;
                case Op::Ceiling:
                    result = detail::increasing(argument(0), [](double v) { return std::ceil(v); });
                    break;
                case Op::Trunc:
                    result = detail::increasing(argument(0), [](double v) { return std::trunc(v); });
                    break;
                case Op::Round:
                    // The engine now implements `round`, and it rounds half to even, so this bound
                    // is the rounding itself rather than the floor-and-ceiling bracket it used to
                    // be. Worth up to one at each end: over an age of 30 to 60, `round(.i$age *
                    // 0.09)` was [2, 6] and is now [3, 5].
                    //
                    // Sound because rounding is monotonic: lo <= x <= hi gives round(lo) <=
                    // round(x) <= round(hi). CORRECT, rather than merely sound, because the engine
                    // rounds the same way -- Op::Round in Interpreter.hpp calls this very function,
                    // so the two cannot drift. An infinite bound rounds to itself, which is
                    // roundBankers' own behaviour and is what keeps a half-bounded interval from
                    // collapsing to unknown.
                    result = detail::increasing(argument(0), [](double v) { return roundBankers(v); });
                    break;
                case Op::Sign:
                    result = detail::increasing(argument(0), [](double v) { return (v > 0.0) - (v < 0.0); });
                    break;

                // Sine and cosine are not monotonic, but they are bounded, which is all a range test
                // usually needs.
                case Op::Sin:
                case Op::Cos:
                    result = Interval::bounds(-1.0, 1.0);
                    break;

                default:
                    // Comparisons, logicals and the NA test all answer 0 or 1. Everything else --
                    // powers, the floored division pair, the time-vector and finalisation operators
                    // -- stays unknown until there is a reason to tighten it.
                    if (opInfo(call.op).resultRule == ResultRule::AlwaysBool) { result = Interval::boolean(); }
                    break;
            }
        }

        intervals[frame.id] = result;
        settled[frame.id] = 1;
    }

    return intervals;
}

} // namespace veil
