// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cmath>
#include <optional>
#include <variant>
#include <vector>
#include "veil/ColumnScan.hpp" // TypeWithConstraints
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"
#include "veil/TypeSpecificConstraint.hpp"

namespace veil
{

// The constant-folding pass. Using the constraints the column scan read off the data, it rewrites in
// place two kinds of node whose value the data settles: a field of a constant, NA-free column becomes
// that constant, and an is_na whose answer is settled becomes true or false. It mutates the tree but
// does not recompute types -- a folded literal keeps the node's type. Range-based comparison folding
// is a separate pass -- see passFoldIntervalComparisons.hpp.

namespace detail
{

// The constant value of a column whose non-NA values are all one value, and which has no NAs. Empty
// otherwise -- a column with NAs is not a clean constant.
inline std::optional<LiteralValue> constantLiteral(const TypeWithConstraints& constraints)
{
    if (!constraints.valuesAreConstant || constraints.hasNAs || !constraints.hasValues)
    {
        return std::nullopt;
    }
    if (const auto* d = std::get_if<DoubleConstraint>(&constraints.specific)) { return LiteralValue(d->min); }
    if (const auto* n = std::get_if<IntConstraint>(&constraints.specific)) { return LiteralValue(n->min); }
    if (const auto* b = std::get_if<BoolConstraint>(&constraints.specific)) { return LiteralValue(b->min); }
    return std::nullopt;
}

// The settled answer of is_na(arg), or empty when it depends on the row. A literal is never NA unless
// it is a NaN double; a field is never NA when its column has no NAs, and always NA when its column
// is all NA.
inline std::optional<bool> isNaAnswer(
    const Tree& tree,
    NodeId argId,
    const std::vector<std::optional<TypeWithConstraints>>& constraints)
{
    const Node& arg = tree.at(argId);

    if (isLit(arg))
    {
        const LitPayload& lit = std::get<LitPayload>(arg.payload);
        if (const auto* d = std::get_if<double>(&lit.value)) { return std::isnan(*d); }
        return false;
    }

    if (isField(arg))
    {
        const ColumnId column = std::get<FieldPayload>(arg.payload).column;
        if (column < constraints.size() && constraints[column].has_value())
        {
            const TypeWithConstraints& c = *constraints[column];
            if (!c.hasNAs) { return false; }
            if (!c.hasValues) { return true; }
        }
    }

    return std::nullopt;
}

} // namespace detail

inline void passFoldConstants(Tree& tree,
                          const std::vector<std::optional<TypeWithConstraints>>& columnConstraints)
{
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        Node& node = tree.at(id);

        if (isField(node))
        {
            const ColumnId column = std::get<FieldPayload>(node.payload).column;
            if (column < columnConstraints.size() && columnConstraints[column].has_value())
            {
                if (const std::optional<LiteralValue> value = detail::constantLiteral(*columnConstraints[column]))
                {
                    // The payload variant is not assignable (a LitPayload holds a const-member
                    // TypeFull), so emplace the new alternative in place.
                    node.payload.emplace<LitPayload>(LitPayload{*value, node.type.value()});
                }
            }
            continue;
        }

        if (isCall(node) && std::get<CallPayload>(node.payload).op == Op::IsNa)
        {
            const NodeId arg = std::get<CallPayload>(node.payload).args[0];
            if (const std::optional<bool> answer = detail::isNaAnswer(tree, arg, columnConstraints))
            {
                node.payload.emplace<LitPayload>(LitPayload{LiteralValue(*answer), TypeFull::createBool()});
            }
        }
    }
}

} // namespace veil
