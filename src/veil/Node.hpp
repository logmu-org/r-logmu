// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <variant>
#include <vector>
#include "veil/ColumnView.hpp" // for ColumnId
#include "veil/Op.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

using NodeId = uint32_t;
using ObjId = uint32_t;

constexpr NodeId invalidNodeId = static_cast<NodeId>(-1);

// A folded constant. `int` carries datey, durationy and category clicks alike, so the accompanying
// `TypeFull` on `LitPayload` is what distinguishes them.
using LiteralValue = std::variant<bool, double, int, std::string>;

struct LitPayload final
{
    LiteralValue value;
    TypeFull type;
};

// A fact access. The name is resolved to a column index when the tree is built, against the
// scanned columns; the name is retained for diagnostics only.
struct FieldPayload final
{
    ColumnId column = 0;
    std::string name;
};

// The time pronoun `.t`.
struct TimePayload final
{
};

struct CallPayload final
{
    Op op = Op::Pos;
    std::vector<NodeId> args;
};

// An opaque concept object spliced by reference -- a mortality or an include -- already lowered to
// the parameters the core needs, so no front-end object crosses.
struct ObjPayload final
{
    ObjId obj = 0;
};

using NodePayload = std::variant<LitPayload, FieldPayload, TimePayload, CallPayload, ObjPayload>;

// Nodes are deliberately mutable: the tree is rewritten in place before lowering (constraint
// folding, algebraic simplification, idiom recognition), so the const-member style used across the
// type headers would work against us here.
struct Node final
{
    NodePayload payload;

    // Filled in by type resolution; empty until then.
    std::optional<TypeFull> type;

    // Filled in by interval propagation during tree rewriting; empty until then. Held as an id into
    // the rewriting pass's own constraint store rather than inline, so the node stays small.
    std::optional<uint32_t> constraint;
};

inline bool isCall(const Node& node) noexcept { return std::holds_alternative<CallPayload>(node.payload); }
inline bool isLit(const Node& node) noexcept { return std::holds_alternative<LitPayload>(node.payload); }
inline bool isField(const Node& node) noexcept { return std::holds_alternative<FieldPayload>(node.payload); }
inline bool isTime(const Node& node) noexcept { return std::holds_alternative<TimePayload>(node.payload); }
inline bool isObj(const Node& node) noexcept { return std::holds_alternative<ObjPayload>(node.payload); }

} // namespace veil
