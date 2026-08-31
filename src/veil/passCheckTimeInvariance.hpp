// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <stdexcept>
#include <variant>
#include <vector>
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"

namespace veil
{

// The select-condition time-invariance check. A select's condition must not depend on time, so the
// branch is chosen once per individual rather than once per time-point. The R front end already
// refuses a time-dependent condition; this is the same guarantee for the language-independent core,
// which a Python or C# front end gets for free.
//
// It reads the flags passTagTimeVarying produced and takes them as an argument rather than computing
// them, so the tagging and the rule it enforces stay separately readable and separately testable.

inline void passCheckTimeInvariance(const Tree& tree, const std::vector<char>& timeVarying)
{
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        const Node& node = tree.at(id);
        if (!isCall(node)) { continue; }

        const CallPayload& call = std::get<CallPayload>(node.payload);
        if (call.op == Op::Select && timeVarying[call.args[0]])
        {
            throw std::runtime_error("veil: a `select` condition must not depend on time.");
        }
    }
}

} // namespace veil
