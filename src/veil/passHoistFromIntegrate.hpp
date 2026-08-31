// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <variant>
#include <vector>
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// Takes the time-invariant factors out of an integral: `integrate(c * x)` becomes `c * integrate(x)`
// whenever `c` does not vary over time.
//
// THIS IS EXACT IN REAL ARITHMETIC, AND IT IS NOT A NUMERICAL TRADE. Integration is linear, and a
// time-invariant factor is by definition constant across the whole of one individual's exposure, so
// it comes out of the integral with nothing lost.
//
// IN FLOATING POINT IT IS NOT BIT-IDENTICAL, and that is worth stating plainly so nobody spends an
// afternoon on it. `sum(c * x_i) * dt` and `c * (sum(x_i) * dt)` distribute the multiplication
// differently, so they round differently and can differ in the last place. Nothing about
// reproducibility is affected: the rewrite is a deterministic function of the tree, so every
// platform and every run hoists the same factors and gets the same answer. What it does mean is that
// a hoisted result need not match an unhoisted reference to the last ULP, which is a property of the
// rewrite rather than a fault in either.
//
// WHY IT IS WORTH MORE THAN IT LOOKS. A weight is usually constant per individual -- a pension
// amount, a lives indicator -- or a constant times a simple function of time. For an AEV that turns
//
//     E = integrate(mu * w)        into    E = w * integrate(mu)
//     V = integrate(mu * (w * w))  into    V = (w * w) * integrate(mu)
//
// which is one vector multiply and one reduction instead of two of each, because the two integrals
// are now the same node and sharing merges them. The weight never enters the time vector at all.
//
// This is also why the AEV recipe writes V as `mu * (w * w)` rather than reusing E's `mu * w`:
// reusing it would look like a saving and would pin `w` inside the integral, out of reach of this.
//
// SCOPE. A factor is peeled only when it is time-invariant AND the other side of the product still
// varies -- hoisting everything would leave an integral of nothing to integrate, which the existing
// broadcast path already handles better. Only `Mul` is walked: addition does not distribute over an
// integral without knowing the exposure length, which is a per-individual value rather than a factor,
// and division would need the same treatment on the reciprocal to be worth having.
//
// It runs AFTER the folds and BEFORE sharing, and it needs time-varying tags taken after the last
// pass that appended a node. It appends nodes itself, so the tags must be taken again afterwards --
// which the pipeline does anyway.

// Rewrites integrals in place and returns how many factors were taken out.
//
// `timeVarying` must cover the tree as it stands on entry. Nodes this pass appends are not in it, and
// are not consulted: every one of them is a product or an integral built from operands already
// classified.
inline size_t passHoistFromIntegrate(Tree& tree, const std::vector<char>& timeVarying)
{
    if (timeVarying.size() != tree.size())
    {
        throw std::runtime_error("veil: hoisting needs time-varying tags covering the whole tree.");
    }

    const size_t originalCount = tree.size();
    size_t hoisted = 0;

    auto isInvariantNumber = [&](NodeId id) {
        if (id >= timeVarying.size() || timeVarying[id] != 0) { return false; }
        const std::optional<TypeFull>& type = tree.at(id).type;
        return type.has_value() && type->type == Type::Double;
    };

    for (NodeId id = 0; id < static_cast<NodeId>(originalCount); ++id)
    {
        {
            const Node& node = tree.at(id);
            if (!isCall(node)) { continue; }
            const CallPayload& call = std::get<CallPayload>(node.payload);
            if (call.op != Op::Integrate || call.args.size() != 1) { continue; }
        }

        // Peel invariant factors off a chain of products, outermost first. Nothing is written until
        // the walk is done, so an integral that yields nothing is left exactly as it was.
        std::vector<NodeId> factors;
        NodeId inner = std::get<CallPayload>(tree.at(id).payload).args[0];

        while (true)
        {
            const Node& product = tree.at(inner);
            if (!isCall(product)) { break; }

            const CallPayload& call = std::get<CallPayload>(product.payload);
            if (call.op != Op::Mul || call.args.size() != 2) { break; }

            const NodeId left = call.args[0];
            const NodeId right = call.args[1];

            // The side that stays must still vary, or there is no integral left worth taking.
            if (isInvariantNumber(left) && right < timeVarying.size() && timeVarying[right] != 0)
            {
                factors.push_back(left);
                inner = right;
                continue;
            }
            if (isInvariantNumber(right) && left < timeVarying.size() && timeVarying[left] != 0)
            {
                factors.push_back(right);
                inner = left;
                continue;
            }
            break;
        }

        if (factors.empty()) { continue; }

        // Build the new integral and the product chain around it. Every append can reallocate the
        // arena, so nothing holds a reference to a node across one; ids are re-fetched instead.
        const TypeFull number = TypeFull::createDouble();

        NodeId current = tree.buildCall(Op::Integrate, {inner});
        tree.at(current).type.emplace(number);

        for (size_t position = factors.size(); position-- > 1;)
        {
            current = tree.buildCall(Op::Mul, {factors[position], current});
            tree.at(current).type.emplace(number);
        }

        // The original integral node BECOMES the outermost product, so every parent pointing at it
        // keeps pointing at the right thing without anything having to find them.
        //
        // THIS IS AN IN-PLACE REWRITE, AND IT IS SAFE FOR A STRONGER REASON THAN "INTEGRALS ARE NOT
        // SHARED YET". It would be safe even for a node with two parents, because it is VALUE-
        // PRESERVING: the node computes exactly what it computed before, by a different route. The
        // DAG rule forbids changing what a shared node MEANS, which is a different thing. Resting on
        // how the AEV recipe happens to build its integrals would be an assumption about a caller;
        // this rests on the rewrite itself.
        Node& original = tree.at(id);
        original.payload = CallPayload{Op::Mul, {factors[0], current}};
        original.type.emplace(number);

        hoisted += factors.size();
    }

    return hoisted;
}

} // namespace veil
