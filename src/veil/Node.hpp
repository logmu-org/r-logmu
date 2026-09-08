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
#include "veil/Operand.hpp" // for ParamId, which the block owns
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

// A FITTED PARAMETER -- one beta of a model the fitter is iterating on.
//
// IT IS NOT A LITERAL, AND THAT IS THE WHOLE POINT. A parameter is constant across every individual
// but changes between runs, so it lowers to a `ConstantBinding` like any other constant. What it
// must NOT do is behave like a literal in the tree:
//
//   - `passFoldConstants` would bake the STARTING value into a folded product, so the block would be
//     compiled for one beta and then asked to answer for another;
//   - `passShareCommonSubtrees` keys a double literal BY ITS BITS and merges equal ones -- and a fit
//     starts every coefficient at zero by default, so all of them would collapse onto ONE operand
//     and setting beta_1 would move beta_2. A silent wrong answer in the default configuration at
//     every k >= 2.
//
// So it is its own leaf, opaque to both. Being a leaf rather than a call it is time-invariant by
// construction, which `passTagTimeVarying` gets right by falling through to zero.
struct ParamPayload final
{
    ParamId param = 0;
};

using NodePayload =
    std::variant<LitPayload, FieldPayload, TimePayload, CallPayload, ObjPayload, ParamPayload>;

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
inline bool isParam(const Node& node) noexcept { return std::holds_alternative<ParamPayload>(node.payload); }

} // namespace veil
