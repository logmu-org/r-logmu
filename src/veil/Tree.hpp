// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// The arena holding one calculation's tree, and the builder API that front ends call to construct
// it. The builders take plain C++ values only -- no front-end object, and in the R case no SEXP,
// ever reaches this far. That makes this API, rather than any wire format, the boundary a Python or
// C# front end reuses: each walks its own source tree and calls the same builders.
//
// A tree is built once per calculation over a small number of nodes, so it favours clarity over
// layout tricks. It is rewritten in place before lowering, which is why `at` hands back a mutable
// reference. Common sub-expression elimination during rewriting lets a node be referenced by more
// than one parent, so a rewritten tree is a DAG; nothing here assumes otherwise.
class Tree final
{
public:
    NodeId add(Node node)
    {
        this->nodes.push_back(std::move(node));
        return static_cast<NodeId>(this->nodes.size() - 1);
    }

    Node& at(NodeId id) { return this->nodes[id]; }
    const Node& at(NodeId id) const { return this->nodes[id]; }

    size_t size() const noexcept { return this->nodes.size(); }

    // A TREE HAS AS MANY ROOTS AS THE CALCULATION HAS OUTPUTS, so it is a forest rather than a tree
    // whenever more than one value is wanted from one pass over the data. An AEV has three roots, a
    // log-likelihood two, its gradient one per parameter and its second differential one per pair.
    // They share a single arena, which is what lets them share sub-expressions -- the mortality is
    // computed once however many outputs read it.
    const std::vector<NodeId>& roots() const noexcept { return this->rootList; }

    void addRoot(NodeId id) { this->rootList.push_back(id); }

    void setRoot(NodeId id)
    {
        this->rootList.clear();
        this->rootList.push_back(id);
    }

    // Replaces the roots wholesale, which only a pass that rewrites the graph has cause to do: common
    // sub-expression sharing can make one root the canonical node for another.
    void setRoots(std::vector<NodeId> ids) { this->rootList = std::move(ids); }

    // The single root, for the callers and the debug entry points that genuinely have one. Asking a
    // forest for "the" root is a mistake worth naming rather than answering with its first.
    NodeId root() const
    {
        if (this->rootList.size() != 1)
        {
            throw std::runtime_error("veil: this tree does not have exactly one root.");
        }
        return this->rootList.front();
    }

    // Builders. A `call` builds its children first, then itself.

    NodeId buildLitBool(bool value)
    {
        return this->add(Node{LitPayload{value, TypeFull::createBool()}, TypeFull::createBool(), {}});
    }

    NodeId buildLitDouble(double value)
    {
        return this->add(Node{LitPayload{value, TypeFull::createDouble()}, TypeFull::createDouble(), {}});
    }

    // Covers datey, durationy and category alike; `type` says which.
    NodeId buildLitInt(int value, TypeFull type)
    {
        return this->add(Node{LitPayload{value, type}, type, {}});
    }

    NodeId buildLitText(std::string value)
    {
        return this->add(Node{LitPayload{std::move(value), TypeFull::createText()}, TypeFull::createText(), {}});
    }

    NodeId buildField(ColumnId column, std::string name)
    {
        return this->add(Node{FieldPayload{column, std::move(name)}, {}, {}});
    }

    NodeId buildTime()
    {
        return this->add(Node{TimePayload{}, {}, {}});
    }

    NodeId buildCall(Op op, std::vector<NodeId> args)
    {
        return this->add(Node{CallPayload{op, std::move(args)}, {}, {}});
    }

    NodeId buildObj(ObjId obj)
    {
        return this->add(Node{ObjPayload{obj}, {}, {}});
    }

private:
    std::vector<Node> nodes;
    std::vector<NodeId> rootList;
};

} // namespace veil
