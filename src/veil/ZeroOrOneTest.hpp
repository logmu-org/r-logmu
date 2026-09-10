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
#include "veil/ColumnScan.hpp"
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/TypeFull.hpp"
#include "veil/TypeSpecificConstraint.hpp"

namespace veil
{

// EXTRACTED SO TWO PASSES CAN SHARE ONE DEFINITION of what "only zero or one" means, rather than two
// that have to agree. `passFoldIndicatorSquares` asks it whether `x * x` may become `x`;
// `passDistributeSquares` asks it whether distributing a square over a product would expose such an
// x at all. A second, looser copy of this test in either place would be a silently wrong answer, not
// a missed optimisation.
//
// AN INTERVAL OF [0, 1] IS NOT ENOUGH, and this is the trap worth naming: a half squares to a
// quarter. The fact needed is "only zero or one".

// Whether a node can only ever be zero or one. Conservative: anything it cannot show, it declines.
class ZeroOrOneTest final
{
public:
    ZeroOrOneTest(
        const Tree& tree,
        const std::vector<std::optional<TypeWithConstraints>>& columnConstraints)
        : tree(tree), columnConstraints(columnConstraints), answers(tree.size(), unknown) {}

    bool operator()(NodeId id) { return this->test(id, 0); }

private:
    static constexpr char unknown = 0;
    static constexpr char yes = 1;
    static constexpr char no = 2;

    // The depth bound stops a cycle from looping for ever. A tree cannot hold one, but this is
    // cheaper than proving that here and it costs nothing on any expression a person writes.
    bool test(NodeId id, int depth)
    {
        if (depth > 32) { return false; }
        if (this->answers[id] != unknown) { return this->answers[id] == yes; }

        const bool answer = this->compute(id, depth);
        this->answers[id] = answer ? yes : no;
        return answer;
    }

    bool compute(NodeId id, int depth)
    {
        const Node& node = this->tree.at(id);

        // A logical is zero or one by what it is.
        if (node.type.has_value() && node.type->type == Type::Bool) { return true; }

        if (const auto* lit = std::get_if<LitPayload>(&node.payload))
        {
            if (const auto* value = std::get_if<double>(&lit->value))
            {
                return *value == 0.0 || *value == 1.0;
            }
            return false;
        }

        if (const auto* field = std::get_if<FieldPayload>(&node.payload))
        {
            return this->columnIsZeroOrOne(field->column);
        }

        if (!isCall(node)) { return false; }
        const CallPayload& call = std::get<CallPayload>(node.payload);

        // A comparison or a logical operator answers zero or one whatever its arguments are.
        if (opInfo(call.op).resultRule == ResultRule::AlwaysBool) { return true; }

        switch (call.op)
        {
            // The conversion a logical takes on its way into arithmetic.
            case Op::ToDouble:
                return call.args.size() == 1 && this->test(call.args[0], depth + 1);

            // Both branches, since either may be the one taken.
            case Op::Select:
                return call.args.size() == 3 && this->test(call.args[1], depth + 1)
                    && this->test(call.args[2], depth + 1);

            // Zero and one are closed under multiplication, and under min and max.
            case Op::Mul:
            case Op::Min:
            case Op::Max:
                return call.args.size() == 2 && this->test(call.args[0], depth + 1)
                    && this->test(call.args[1], depth + 1);

            default:
                return false;
        }
    }

    bool columnIsZeroOrOne(ColumnId column) const
    {
        if (column >= this->columnConstraints.size()) { return false; }
        const std::optional<TypeWithConstraints>& constraint = this->columnConstraints[column];
        if (!constraint.has_value() || constraint->hasNAs || !constraint->hasValues) { return false; }

        if (std::holds_alternative<BoolConstraint>(constraint->specific)) { return true; }

        const auto* range = std::get_if<DoubleConstraint>(&constraint->specific);
        if (range == nullptr) { return false; }

        // Integral AND inside [0, 1]: either alone would admit a half or a two.
        return range->allIntegral && range->min >= 0.0 && range->max <= 1.0;
    }

    const Tree& tree;
    const std::vector<std::optional<TypeWithConstraints>>& columnConstraints;
    std::vector<char> answers;
};

} // namespace veil
