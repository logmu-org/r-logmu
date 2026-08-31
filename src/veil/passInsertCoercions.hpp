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
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// The coercion insertion pass. It makes the implicit numeric widening the annotation pass recorded
// explicit: any value-operator operand that is not already double is wrapped in a ToDouble node, so
// every value operator is left with operands of one matching type.
//
// Only value-producing arithmetic and select branches are coerced. Comparisons are left as they are
// -- a comparison is a bool result however it executes, and its widen-versus-narrow choice belongs to
// the narrowing pass. ToDouble is source-aware (datey / durationy read as years, bool as 0 / 1); this
// pass only tags the operand, and lowering does the conversion.
//
// It reads the types the annotation pass wrote and never recomputes them.

namespace detail
{

// Whether operand `position` of a call must be a double, given the call's own resolved type. A call
// whose result is not double is either pure datey / durationy arithmetic or a non-value operator, and
// in both cases nothing widens.
inline bool positionWantsDouble(const OpInfo& info, Op op, const TypeFull& resultType, size_t position)
{
    if (resultType.type != Type::Double) { return false; }
    if (info.category == OpCategory::Arithmetic) { return true; }
    if (op == Op::Select) { return position != 0; } // The condition is a bool test, not a value.
    return false; // Comparisons, logicals, conversions, time-vector and finalise ops are left alone.
}

struct Coercion final
{
    NodeId parent;
    size_t position;
    NodeId operand;
};

} // namespace detail

inline void passInsertCoercions(Tree& tree)
{
    const size_t originalCount = tree.size();

    // Decide over the original nodes only, so appending ToDouble nodes cannot disturb the walk and no
    // inserted node is ever examined.
    std::vector<detail::Coercion> coercions;
    for (NodeId id = 0; id < static_cast<NodeId>(originalCount); ++id)
    {
        const Node& node = tree.at(id);
        if (!isCall(node)) { continue; }

        const CallPayload& call = std::get<CallPayload>(node.payload);
        const OpInfo& info = opInfo(call.op);
        const TypeFull resultType = node.type.value(); // A call always carries a value type.

        for (size_t position = 0; position < call.args.size(); ++position)
        {
            if (!detail::positionWantsDouble(info, call.op, resultType, position)) { continue; }

            const NodeId operand = call.args[position];
            const std::optional<TypeFull>& operandType = tree.at(operand).type;
            if (operandType.has_value() && operandType->type != Type::Double)
            {
                coercions.push_back({id, position, operand});
            }
        }
    }

    // Insert. Each decision appends a ToDouble node and repoints its parent's argument at it. The
    // recorded operand ids stay valid because existing nodes are never moved, only appended to.
    for (const detail::Coercion& coercion : coercions)
    {
        const NodeId coerced = tree.buildCall(Op::ToDouble, {coercion.operand});
        tree.at(coerced).type.emplace(TypeFull::createDouble());
        Node& parent = tree.at(coercion.parent);
        std::get<CallPayload>(parent.payload).args[coercion.position] = coerced;
    }
}

} // namespace veil
