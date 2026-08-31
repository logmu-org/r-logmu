// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstdint>
#include <variant>
#include <vector>
#include "veil/MortalityTable.hpp"
#include "veil/Node.hpp"

namespace veil
{

// The concept objects a front end may splice into a tree, reduced to the parameters the core needs.
// No front-end object reaches this far: the binding reads a mortality or an include out into one of
// these before the tree ever refers to it.
//
// MortalityTable has a header of its own because a lowered block holds one directly, and the block
// is the boundary the tree stops at -- so it must not have to include the tree's own types to say so.

// A mortality that is the same log mu everywhere. It lowers to a scalar, and is broadcast only
// where a time vector is actually wanted.
struct MortalityConst final
{
    double logMu = 0.0;
};

// An include is a conjunction of terms. Absolute and offset terms contribute interval bounds;
// a gate contributes a time-invariant condition which, when false, empties the include for that
// individual. The geometry crosses as bounds and expressions, never as a front-end object.

// `t` lies in [from, to), with both bounds datey clicks.
struct AbsoluteTerm final
{
    int fromClicks = 0;
    int toClicks = 0;
};

// `t - offset` lies in [from, to), where the offset is always SUBTRACTED from t, so the durationy
// bounds are never negative. `age(from, to)` is this with the offset expression `.i$birth`.
//
// THE OFFSET IS AN EXPRESSION, NOT A COLUMN, already ingested into the tree exactly as a gate is. A
// bare field is the common case and lowers to a column read, but nothing here requires one:
// `.t - min(.i$entry, .i$retirement)` is as good an origin as `.t - .i$entry`, and a user who has
// written the second is entitled to be surprised that the first was ever refused.
//
// It must be time-invariant and datey-valued, both checked at lowering. Time-invariance is not a
// style rule: the offset is what DECIDES where the exposure starts, so one that sampled time would
// need the grid that it is itself being used to build.
struct OffsetTerm final
{
    NodeId offset = invalidNodeId;
    int fromClicks = 0;
    int toClicks = 0;
};

// A time-invariant indicator, already ingested into the tree.
struct GateTerm final
{
    NodeId ast = invalidNodeId;
};

using IncludeTerm = std::variant<AbsoluteTerm, OffsetTerm, GateTerm>;

struct Include final
{
    std::vector<IncludeTerm> terms;
};

using Obj = std::variant<MortalityConst, MortalityTable, Include>;

// Objects are owned per calculation, alongside the tree that refers to them by ObjId.
class ObjStore final
{
public:
    ObjId add(Obj obj)
    {
        this->objs.push_back(std::move(obj));
        return static_cast<ObjId>(this->objs.size() - 1);
    }

    const Obj& at(ObjId id) const { return this->objs[id]; }
    size_t size() const noexcept { return this->objs.size(); }

private:
    std::vector<Obj> objs;
};

} // namespace veil
