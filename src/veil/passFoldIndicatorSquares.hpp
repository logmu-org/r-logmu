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

// `x * x` becomes `x` wherever `x` can only be zero or one, since both square to themselves.
//
// WHAT IT IS ACTUALLY FOR. An AEV's second moment is `integrate(mu * (w * w))`, and for an indicator
// weight -- a lives count, a subset flag -- that is the same integral as E. Folding the square makes
// V's expression identical to E's, and the sharing pass then makes V free rather than merely cheap.
// So V = E falls out of an algebraic rewrite plus sharing, with nothing in the recipe knowing that
// the case exists. Nothing here mentions AEV, and a log-likelihood weighted the same way gets it too.
//
// AN INTERVAL OF [0, 1] IS NOT ENOUGH, and this is the trap worth naming: a half squares to a
// quarter. The fact needed is "only zero or one", which comes from one of
//
//   - a logical, or anything the op table marks AlwaysBool -- a comparison, a logical operator;
//   - `to_double` of one of those, which is how a logical reaches arithmetic after coercion;
//   - a literal zero or one;
//   - a DOUBLE COLUMN THE SCAN FOUND HOLDING ONLY ZEROS AND ONES, which is the case R cannot prove
//     for itself and the reason this is worth doing in the engine at all;
//   - a select whose branches are both such a thing, or a product of two of them.
//
// A COLUMN WITH MISSING VALUES IS REFUSED, though the arithmetic would survive one: a NaN squares to
// a NaN, so the rewrite holds there too. It is refused because "only zero or one" should mean what it
// says -- a column whose values include NA is not a column of zeros and ones, and leaning on NaN
// algebra to make the claim true would be the sort of reasoning that stops being true when someone
// changes something else.
//
// WHERE IT RUNS: AFTER SHARING, AND SHARING RUNS AGAIN AFTER IT. Both halves matter. It needs sharing
// first because `w * w` on a logical weight arrives as a product of two SEPARATE `to_double` nodes --
// coercion inserts one per argument -- and only sharing makes them the same node, which is what lets
// this pass see a square rather than a product of two different things. And it needs sharing again
// afterwards because folding V's square leaves V's expression identical to E's, which is a new
// redundancy that did not exist when sharing last ran. That is the concrete case for CSE running more
// than once, rather than a general principle.
//
// IT REPOINTS PARENTS; IT DOES NOT REWRITE THE PRODUCT NODE. Sharing has already made the tree a DAG,
// so changing a node's meaning in place would corrupt every other user of it. Repointing a parent at
// a DIFFERENT node of EQUAL VALUE is a different thing and is always safe -- the parent computes what
// it always computed. The product node simply stops being referenced.

namespace detail
{

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

} // namespace detail

// Repoints every reference to `x * x` at `x`, where `x` can only be zero or one. Returns how many
// references moved.
inline size_t passFoldIndicatorSquares(
    Tree& tree, const std::vector<std::optional<TypeWithConstraints>>& columnConstraints)
{
    if (tree.size() == 0) { return 0; }

    detail::ZeroOrOneTest isZeroOrOne(tree, columnConstraints);

    // What each node should be read as instead of itself. Worked out once, so a square referenced
    // from several places is settled the same way everywhere.
    std::vector<NodeId> replacement(tree.size(), invalidNodeId);
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        const Node& node = tree.at(id);
        if (!isCall(node)) { continue; }

        const CallPayload& call = std::get<CallPayload>(node.payload);
        if (call.op != Op::Mul || call.args.size() != 2 || call.args[0] != call.args[1]) { continue; }

        if (isZeroOrOne(call.args[0])) { replacement[id] = call.args[0]; }
    }

    // FOLLOWED THROUGH TO THE END, because a replacement can itself have one: `(x*x)*(x*x)` is a
    // square of `x*x`, which is a square of `x`. Repointing one step would leave the parent reading
    // `x*x` -- the right VALUE, so no answer changes, but not the same node as `x`, and being the
    // same node is the entire purpose. V would go on being paid for.
    //
    // The walk terminates because a replacement is always an ARGUMENT of the node it replaces, so
    // each step moves strictly down the graph.
    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        NodeId target = replacement[id];
        while (target != invalidNodeId && replacement[target] != invalidNodeId)
        {
            target = replacement[target];
        }
        replacement[id] = target;
    }

    size_t moved = 0;

    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        if (!isCall(tree.at(id))) { continue; }
        CallPayload& call = std::get<CallPayload>(tree.at(id).payload);
        for (NodeId& argId : call.args)
        {
            if (replacement[argId] != invalidNodeId)
            {
                argId = replacement[argId];
                ++moved;
            }
        }
    }

    std::vector<NodeId> roots = tree.roots();
    for (NodeId& rootId : roots)
    {
        if (rootId != invalidNodeId && replacement[rootId] != invalidNodeId)
        {
            rootId = replacement[rootId];
            ++moved;
        }
    }
    tree.setRoots(std::move(roots));

    return moved;
}

} // namespace veil
