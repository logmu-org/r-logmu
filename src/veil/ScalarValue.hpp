// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstdint>
#include <stdexcept>
#include <variant>
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// The value a scalar operand holds while a block runs.
//
// Three storage forms cover every veil base type: `double` for a number, `int64_t` for the
// click-backed types (datey and durationy) and for a category or text code, and `bool` for a logical.
// Which form a given operand holds follows from its TypeFull and is fixed before anything runs, so
// the variant is not a run-time decision -- it is a check that lowering assigned the operand the type
// its instructions go on to treat it as having.
//
// CLICKS ARE 64-BIT INSIDE THE ENGINE, THOUGH THEY ARRIVE AND LEAVE AS 32 (Tim, 2026-07-28). A datey
// or durationy column is R's own 32-bit integer and is widened once on load; the only way back to 32
// bits is an output crossing to R, which is checked there. In between, nothing can overflow and so
// nothing checks.
//
// WHY IT CANNOT OVERFLOW, rather than merely being unlikely to: click arithmetic has NO
// MULTIPLICATION. The click ops are Pos, Neg, Abs, Add, Sub, Min, Max and Clamp, and of those only
// the first four produce a new value, so every click expression is a sum of inputs with coefficients
// of plus or minus one. The largest input magnitude is a valid datey, about 1.6e9, against an int64
// range of 9.2e18 -- so overflow needs billions of chained terms, where an expression tree has tens
// of nodes. It costs nothing: scalar 64-bit integer arithmetic is the same instruction at the same
// latency as 32-bit on any current target, and clicks never vectorise, so no SIMD width is given up.
using ScalarValue = std::variant<double, int64_t, bool>;

// The form an operand of `type` holds, zeroed. A register file starts out this way, so every read
// finds the alternative the operand's type promises whether or not an instruction has written to it.
inline ScalarValue zeroValueOf(const TypeFull& type)
{
    switch (type.type)
    {
        case Type::Bool: return ScalarValue(false);
        case Type::Double: return ScalarValue(0.0);
        case Type::Datey:
        case Type::Durationy:
        case Type::Category:
        case Type::Text: return ScalarValue(int64_t{0});
        case Type::DateyInterval: break; // An include resolves to bounds and gates, never to a value.
    }
    throw std::runtime_error("veil: a datey_interval cannot be held in an operand.");
}

} // namespace veil
