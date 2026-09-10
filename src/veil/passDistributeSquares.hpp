// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <optional>
#include <variant>
#include <vector>
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/TypeFull.hpp"
#include "veil/TypeSpecificConstraint.hpp"
#include "veil/ZeroOrOneTest.hpp"

namespace veil
{

// `(a * b) * (a * b)` becomes `(a * a) * (b * b)`, but ONLY where one of the factors can be nothing
// but zero or one.
//
// WHAT IT IS FOR. The commonest covariate shape in a fitted model is an indicator times a shape --
// five regions each carrying the same age curve, `X^(j) = I^(j) phi`. The information matrix wants
// `X^(j) X^(j)` on its diagonal, and written as it stands that is `(I phi)(I phi)`: a square of a
// product, where `passFoldIndicatorSquares` fires only on a square of a node that is ITSELF zero or
// one, and `I phi` is not. The `I^2` is there and is invisible.
//
// Distributing the square exposes it. `(I*I) * (phi*phi)` has `I*I` as a bare square of an
// indicator, which the fold then collapses to `I`, leaving `I * (phi*phi)`. `I` is time-invariant,
// so the hoist lifts it clean out of the integral -- and what remains, `integrate(mu w phi^2)`, is
// the SAME NODE for every j. Sharing then makes k diagonal integrals into one.
//
// So three passes that already existed do the work; this one only puts the expression in a shape
// where they can see it.
//
// IT MUST RUN BEFORE THE HOIST, which is what actually banks the saving, and the hoist runs before
// sharing. That is the whole of its placement.
//
// WHY THE INDICATOR CONDITION IS NOT AN OPTIMISATION BUT A GUARD. Without a factor that folds,
// distributing buys nothing at all -- no square collapses, no factor becomes time-invariant that was
// not before -- and it is not free: `(ab)(ab)` and `(aa)(bb)` group the multiplications differently
// and round differently, exactly as the hoist does. Rewriting every squared product in the package
// to gain nothing on almost all of them would be a poor trade. Firing only where a fold follows
// keeps the rewrite where it pays.
//
// AND IT BOUNDS THE ONE REAL NUMERICAL DIFFERENCE. `b*b` can overflow to infinity where `(ab)^2`
// would not, if `a` is zero on that record and `|b|` exceeds about 1.3e154; the old form gives zero
// and the new one gives `0 * inf`, a NaN. That needs a covariate at a scale no experience analysis
// has, and the package already squares user quantities directly -- an AEV's V is `mu * (w * w)`.
// Worth knowing about rather than worth guarding against, and it is why this note exists.
//
// EACH FACTOR IS COERCED ON ITS WAY INTO THE SQUARE, and that is not tidiness. `.x` is an age,
// which is a `durationy`, and `durationy * durationy` is a product datey does not define -- so
// `(I * .x)` squared is perfectly legal while `.x * .x` is not. Since `I * .x` is the single most
// common covariate there is, declining click-backed factors would decline the whole reason this
// pass exists. Wrapping each factor in `ToDouble` first is what the fit recipe already does one
// level up, it reads a click-backed value as years, and it preserves the value exactly:
// `(I .x)^2` and `I^2 * years(.x)^2` are the same number.
//
// IT APPENDS AND REPOINTS; IT NEVER REWRITES A NODE IN PLACE. The recipe already references one
// covariate from many roots, so changing a node's meaning would corrupt every other user of it.
// Building new nodes of the value the parent already had is always safe.
//
// TWO THINGS HERE HAVE NO TEST, and that is recorded rather than papered over -- see
// `tests/testthat/test-veil_distribute_squares.R`, where the rest of this is witnessed.
//
//   * THE COERCION CHANGES NO ANSWER TODAY. Dropping it was measured and moved nothing, on either
//     operand: a duration on the time vector already holds years, and `passLowerToBlock` converts a
//     scalar click-backed operand when a double is asked of it. What the coercion buys is that the
//     nodes appended here are WELL TYPED IN THEIR OWN RIGHT rather than correct by the lowering's
//     leniency -- which matters because the types below are written by hand, so nothing else will
//     ever check them.
//   * REPOINTING THE ROOTS is unreachable from any real calculation, where every root is an
//     `integrate` or a `died_value` and never a bare square. It is here because
//     `passFoldIndicatorSquares` does the same and a rewrite that repointed parents but not roots
//     would be a trap for whoever next builds a tree whose root is an expression.

namespace detail
{

// The two factors of a `Mul`, seen through any coercion wrapping the product itself. Returns false
// for anything that is not a two-argument `Mul`.
inline bool multiplyFactors(const Tree& tree, NodeId id, NodeId& left, NodeId& right)
{
    const Node& node = tree.at(id);
    if (!isCall(node)) { return false; }

    const CallPayload& call = std::get<CallPayload>(node.payload);
    if (call.op == Op::ToDouble && call.args.size() == 1)
    {
        return multiplyFactors(tree, call.args[0], left, right);
    }
    if (call.op != Op::Mul || call.args.size() != 2) { return false; }

    left = call.args[0];
    right = call.args[1];
    return true;
}

// Whether a factor can go into a square at all: a plain double already can, and a click-backed one
// can once it is read as years. Anything else -- text, an interval, an untyped node -- cannot.
inline bool canBeSquared(const Tree& tree, NodeId id)
{
    const std::optional<TypeFull>& type = tree.at(id).type;
    if (!type.has_value()) { return false; }
    return type->type == Type::Double || type->type == Type::Datey
        || type->type == Type::Durationy;
}

inline bool needsCoercion(const Tree& tree, NodeId id)
{
    const std::optional<TypeFull>& type = tree.at(id).type;
    return type.has_value() && type->type != Type::Double;
}

} // namespace detail

// Rewrites qualifying squares and returns how many were rewritten.
inline size_t passDistributeSquares(
    Tree& tree, const std::vector<std::optional<TypeWithConstraints>>& columnConstraints)
{
    if (tree.size() == 0) { return 0; }

    ZeroOrOneTest isZeroOrOne(tree, columnConstraints);

    // Worked out over the tree AS IT STANDS, before anything is appended, so the scan cannot trip
    // over its own new nodes -- and a Node reference is never held across an append, which would
    // dangle the moment the arena grows.
    std::vector<NodeId> replacement(tree.size(), invalidNodeId);
    const NodeId originalSize = static_cast<NodeId>(tree.size());

    for (NodeId id = 0; id < originalSize; ++id)
    {
        NodeId squared = invalidNodeId;
        {
            const Node& node = tree.at(id);
            if (!isCall(node)) { continue; }

            const CallPayload& call = std::get<CallPayload>(node.payload);
            if (call.op != Op::Mul || call.args.size() != 2 || call.args[0] != call.args[1])
            {
                continue;
            }
            squared = call.args[0];
        }

        // A square of something ALREADY zero or one is `passFoldIndicatorSquares`'s to collapse
        // outright, and collapsing beats distributing. Left alone.
        if (isZeroOrOne(squared)) { continue; }

        NodeId left = invalidNodeId;
        NodeId right = invalidNodeId;
        if (!detail::multiplyFactors(tree, squared, left, right)) { continue; }

        // `a * a` distributed over itself is `(a*a) * (a*a)`, which is the same expression again.
        if (left == right) { continue; }

        // Exactly the payoff condition: one factor folds, so the rewrite exposes something.
        //
        // A LITERAL DOES NOT COUNT, even though 0 and 1 are the most zero-or-one values there are.
        // `(a * 1)` squared distributes to `(a*a) * (1*1)`, which folds to `(a*a) * 1` -- a constant
        // split off from a product, which constant folding and the hoist already deal with far more
        // cheaply than appending five nodes to say it. The point of this pass is to expose an
        // indicator that VARIES BETWEEN INDIVIDUALS, because that is the one the hoist can lift out
        // to leave a shared integral behind.
        const auto qualifies = [&](NodeId id) {
            return isZeroOrOne(id) && !std::holds_alternative<LitPayload>(tree.at(id).payload);
        };
        if (!qualifies(left) && !qualifies(right)) { continue; }

        if (!detail::canBeSquared(tree, left) || !detail::canBeSquared(tree, right)) { continue; }

        // NOTHING IS APPENDED UNTIL EVERY TEST HAS PASSED, so a declined square leaves the tree
        // exactly as it was rather than littering it with nodes nobody references.
        const bool coerceLeft = detail::needsCoercion(tree, left);
        const bool coerceRight = detail::needsCoercion(tree, right);

        const NodeId leftDouble = coerceLeft ? tree.buildCall(Op::ToDouble, {left}) : left;
        const NodeId rightDouble = coerceRight ? tree.buildCall(Op::ToDouble, {right}) : right;
        const NodeId leftSquared = tree.buildCall(Op::Mul, {leftDouble, leftDouble});
        const NodeId rightSquared = tree.buildCall(Op::Mul, {rightDouble, rightDouble});
        const NodeId product = tree.buildCall(Op::Mul, {leftSquared, rightSquared});

        // TYPED HERE RATHER THAN BY RE-RUNNING ANNOTATION. Every node built above is a double by
        // construction -- a coercion to double, or a product of two of them -- so there is nothing
        // to derive. Re-annotating the whole tree at this point would also have to redo the text
        // encoding and the coercions that have already been applied.
        for (const NodeId built : {leftDouble, rightDouble, leftSquared, rightSquared, product})
        {
            tree.at(built).type.emplace(TypeFull::createDouble());
        }

        replacement[id] = product;
    }

    size_t rewritten = 0;
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        if (!isCall(tree.at(id))) { continue; }
        for (NodeId& argId : std::get<CallPayload>(tree.at(id).payload).args)
        {
            if (argId < originalSize && replacement[argId] != invalidNodeId)
            {
                argId = replacement[argId];
                ++rewritten;
            }
        }
    }

    // A root may be one of these squares in its own right, which no parent would reach.
    std::vector<NodeId> roots = tree.roots();
    for (NodeId& rootId : roots)
    {
        if (rootId != invalidNodeId && rootId < originalSize
            && replacement[rootId] != invalidNodeId)
        {
            rootId = replacement[rootId];
            ++rewritten;
        }
    }
    tree.setRoots(std::move(roots));

    return rewritten;
}

} // namespace veil
