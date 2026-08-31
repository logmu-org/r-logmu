// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <variant>
#include <vector>
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"

namespace veil
{

// The time-varying tagging analysis. It marks each node as time-varying (a value that changes across
// the integration time-points -- a vector) or time-invariant (a per-individual scalar). Time enters
// through the `.t` pronoun and the vector-source ops (vector_log_mu and friends); a finalise op
// collapses a vector back to a scalar; everything else is time-varying exactly when one operand is.
//
// Non-mutating and re-runnable, like the type annotation pass. Its result feeds passCheckTimeInvariance
// and, later, lowering's scalar-versus-vector operand split.

namespace detail
{

// A memoised post-order walk. The coercion pass has appended nodes, so the children-before-parents
// ordering no longer holds and a forward scan would read an operand before computing it; recursion
// with a visited flag copes with any ordering and with the shared nodes a DAG carries.
inline char visitTimeVarying(const Tree& tree, NodeId id, std::vector<char>& flags, std::vector<char>& done)
{
    if (done[id]) { return flags[id]; }

    const Node& node = tree.at(id);
    char result = 0;

    if (isTime(node))
    {
        result = 1;
    }
    else if (isCall(node))
    {
        const CallPayload& call = std::get<CallPayload>(node.payload);
        if (opInfo(call.op).resultRule == ResultRule::AlwaysVector)
        {
            result = 1; // A vector source, whatever its arguments.
        }
        else if (call.op == Op::Integrate || call.op == Op::DiedValue)
        {
            result = 0; // A finalise collapses a vector to a scalar.
        }
        else
        {
            for (const NodeId argId : call.args)
            {
                if (visitTimeVarying(tree, argId, flags, done)) { result = 1; break; }
            }
        }
    }
    // A lit, a field and an obj node are per-individual scalars.

    flags[id] = result;
    done[id] = 1;
    return result;
}

} // namespace detail

// Returns a flag per node: 1 if the node is time-varying (a vector), 0 if time-invariant (a scalar).
inline std::vector<char> passTagTimeVarying(const Tree& tree)
{
    std::vector<char> flags(tree.size(), 0);
    std::vector<char> done(tree.size(), 0);
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        detail::visitTimeVarying(tree, id, flags, done);
    }
    return flags;
}

} // namespace veil
