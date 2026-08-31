// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// Common sub-expression sharing: the tree becomes a DAG, so a value computed in two places is
// computed once.
//
// WHAT IT IS FOR, AND WHAT IT IS NOT FOR. Sharing that a caller constructs deliberately -- an AEV
// recipe splicing the same mortality node into both E and V -- needs nothing from this pass, because
// lowering already memoises by NodeId and would compute it once anyway. This pass earns its keep on
// sharing NOBODY ARRANGED: the same condition written into an include and a weight, a batch of
// several calculations over one mortality, a sub-expression a user simply wrote out twice. That is
// also why it is a general pass over the whole forest rather than anything a recipe calls.
//
// IT WILL NOT RESCUE A BADLY ASSOCIATED EXPRESSION. `(a*b)*c` and `a*(b*c)` are genuinely different
// subtrees and stay that way; no reassociation happens here. Anything wanting the same value out of
// two expressions has to write them the same way round.
//
// WHERE IT RUNS, AND WHY THAT IS LAST. Every pass before it may rewrite a node in place -- constant
// folding replaces a call with a literal, interval folding settles a comparison. Once a node has two
// parents, an in-place rewrite corrupts every other user of it, which is exactly the fault found in
// passNarrowComparisons on 2026-07-27. So this pass runs after all of them, and the rule it
// establishes is:
//
//     ONCE THE TREE IS A DAG, NO PASS MAY MUTATE A NODE IN PLACE.
//
// A pass that must change something builds a new node and repoints the parent, which is what
// passNarrowComparisons was corrected to do. Running this pass a SECOND time is safe and cheap --
// it is idempotent and only ever reduces the node count -- so once there are rewrites that expose
// fresh redundancy, a second run after them is the right answer rather than a reason to move this
// one earlier.
//
// The memoising passes need nothing: tagging, interval propagation and lowering all settle a shared
// node once and reuse the answer.

namespace detail
{

// What makes two nodes the same node. Everything that can distinguish their VALUE, and nothing that
// cannot: two nodes agreeing on all of this compute the same thing for every individual at every
// instant, so either may stand for both.
//
// The key is built from ALREADY-CANONICAL argument ids, so structural equality falls out of a single
// bottom-up pass -- if the arguments have been shared, equal argument ids mean equal sub-expressions.
//
// A TYPE MISMATCH KEEPS TWO NODES APART even when they look identical otherwise. Two literals holding
// the integer 1 are a different value if one is a durationy and the other a category code, and the
// annotation pass has already settled which each is.
// EVERY PART OF A TYPE IS KEYED, not just its base. `max` carries a category's cardinality, and
// while two nodes agreeing on everything else here cannot today disagree about it -- annotation being
// a deterministic function of the structure this key already covers -- that is a property of another
// pass, held nowhere near this one. Keying it outright costs a tuple field and means the merge no
// longer depends on a theorem nobody restates.
//
// A LITERAL'S OWN TYPE IS A SEPARATE FIELD from the node's, rather than the two being packed into one
// by arithmetic. They can differ, both must distinguish, and encoding a pair as `a * 1000 + b` is
// exactly the shape a collision hides in even where the arithmetic happens to be safe.
//
// A DOUBLE LITERAL IS KEYED BY ITS BITS, NOT ITS VALUE, AND THAT IS LOAD-BEARING. `std::map` orders
// with `<`, and `<` on doubles is NOT a strict weak ordering once NaN is in play: NaN is unordered
// against everything, so `!(NaN < x) && !(x < NaN)` makes the tuple comparison call NaN EQUIVALENT to
// whatever double it is compared against on the way down the tree. That is undefined behaviour by the
// standard's own precondition, and it had a visible wrong answer attached (found 2026-08-20): a
// literal NaN merged with the AEV recipe's literal 1.0, so `mortality = NaN` computed exp(1) per unit
// of exposure and returned a perfectly finite E. Bits give a total order and NaN merges only with the
// identical bit pattern.
//
// Two consequences, both in the safe direction. `+0.0` and `-0.0` stop sharing, which is right rather
// than merely harmless -- they differ under division, so they are different values. Distinct NaN
// payloads (R's NA_real_ against a plain NaN) stop sharing too. Both cost a node and nothing else;
// this pass may always share LESS than it could, never more.
using ShareKey = std::tuple<int,                 // which sort of node
                            int,                 // op, or the payload discriminator for a literal
                            std::string,         // a text literal's value, or a field's column
                            std::uint64_t,       // a numeric literal's value, AS BITS -- see above
                            long long,           // an integer or boolean literal, an obj id, a column
                            int,                 // the node's own type
                            int,                 // the node's own type's max, a category's cardinality
                            int,                 // a literal payload's type
                            int,                 // a literal payload's type's max
                            std::vector<NodeId>>; // canonical arguments

inline int typeTag(const std::optional<TypeFull>& type)
{
    return type.has_value() ? static_cast<int>(type->type) : -1;
}

inline int typeMaxTag(const std::optional<TypeFull>& type)
{
    return type.has_value() ? type->max : -1;
}

inline ShareKey shareKeyOf(const Node& node, std::vector<NodeId> canonicalArgs)
{
    const int selfType = typeTag(node.type);
    const int selfMax = typeMaxTag(node.type);

    if (const auto* lit = std::get_if<LitPayload>(&node.payload))
    {
        // The literal's own alternative is part of the key: a double 1.0 and an integer 1 are not
        // interchangeable even where both would print the same.
        const int which = static_cast<int>(lit->value.index());
        std::string text;
        std::uint64_t numberBits = 0;
        long long integer = 0;

        if (const auto* s = std::get_if<std::string>(&lit->value)) { text = *s; }
        else if (const auto* d = std::get_if<double>(&lit->value)) { numberBits = std::bit_cast<std::uint64_t>(*d); }
        else if (const auto* i = std::get_if<int>(&lit->value)) { integer = *i; }
        else if (const auto* b = std::get_if<bool>(&lit->value)) { integer = *b ? 1 : 0; }

        return ShareKey{0, which, std::move(text), numberBits, integer, selfType, selfMax,
                        static_cast<int>(lit->type.type), lit->type.max, {}};
    }

    if (const auto* field = std::get_if<FieldPayload>(&node.payload))
    {
        return ShareKey{1, 0, field->name, 0, static_cast<long long>(field->column),
                        selfType, selfMax, -1, -1, {}};
    }

    if (isTime(node)) { return ShareKey{2, 0, {}, 0, 0, selfType, selfMax, -1, -1, {}}; }

    if (const auto* obj = std::get_if<ObjPayload>(&node.payload))
    {
        // Two references to the SAME stored object share; two structurally identical objects stored
        // separately do not, since nothing here can see inside one.
        return ShareKey{3, 0, {}, 0, static_cast<long long>(obj->obj),
                        selfType, selfMax, -1, -1, {}};
    }

    const CallPayload& call = std::get<CallPayload>(node.payload);
    return ShareKey{4, static_cast<int>(call.op), {}, 0, 0, selfType, selfMax, -1, -1,
                    std::move(canonicalArgs)};
}

} // namespace detail

// Rewrites the forest so that structurally identical sub-expressions become one node, and returns
// how many nodes stopped being reachable as a result.
//
// Nothing is deleted. The arena keeps every node it ever held, because a NodeId is an index into it
// and shrinking it would invalidate every id anything else is holding. The unreachable ones simply
// stop being walked, and lowering only ever visits what a root reaches.
inline size_t passShareCommonSubtrees(Tree& tree)
{
    if (tree.size() == 0) { return 0; }

    std::vector<NodeId> canonical(tree.size(), invalidNodeId);
    std::map<detail::ShareKey, NodeId> seen;

    // Post-order, from every root, over a graph that may already share nodes. `canonical` doubles as
    // the visited marker, so a node reached twice is settled once.
    struct Frame final
    {
        NodeId id;
        bool childrenScheduled;
    };

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

        if (canonical[frame.id] != invalidNodeId) { continue; }

        const Node& node = tree.at(frame.id);

        if (!frame.childrenScheduled && isCall(node))
        {
            stack.push_back(Frame{frame.id, true});
            for (const NodeId argId : std::get<CallPayload>(node.payload).args)
            {
                if (canonical[argId] == invalidNodeId) { stack.push_back(Frame{argId, false}); }
            }
            continue;
        }

        // Point this node's arguments at their canonical selves BEFORE keying it, so that two
        // sub-expressions which have just been merged make their parents equal too. This is the one
        // in-place edit the pass makes, and it is safe because a node is rewritten exactly once --
        // whereas the node it might be MERGED INTO is never touched.
        std::vector<NodeId> args;
        if (isCall(node))
        {
            CallPayload& call = std::get<CallPayload>(tree.at(frame.id).payload);
            for (NodeId& argId : call.args)
            {
                const NodeId target = canonical[argId];
                if (target == invalidNodeId)
                {
                    throw std::runtime_error("veil: sharing reached a parent before its argument.");
                }
                argId = target;
            }
            args = call.args;
        }

        detail::ShareKey key = detail::shareKeyOf(tree.at(frame.id), std::move(args));
        const auto found = seen.find(key);
        if (found != seen.end())
        {
            canonical[frame.id] = found->second;
        }
        else
        {
            seen.emplace(std::move(key), frame.id);
            canonical[frame.id] = frame.id;
        }
    }

    // A root can itself be shared -- two outputs that turn out to be the same value, which is exactly
    // what an AEV's V and E become when the weight is an indicator.
    std::vector<NodeId> roots;
    roots.reserve(tree.roots().size());
    size_t shared = 0;
    for (const NodeId rootId : tree.roots())
    {
        roots.push_back(rootId == invalidNodeId ? rootId : canonical[rootId]);
    }
    tree.setRoots(std::move(roots));

    for (NodeId id = 0; id < static_cast<NodeId>(tree.size()); ++id)
    {
        if (canonical[id] != invalidNodeId && canonical[id] != id) { ++shared; }
    }
    return shared;
}

} // namespace veil
