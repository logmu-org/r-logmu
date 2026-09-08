// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstdint>
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// Operands are block-local: an id indexes the arena of the block that declared it, never a global
// table. Two concurrent invocations of a block, or two iterations of a host-driven solve, each hold
// their own arena, so ids from different invocations must never be compared or exchanged.
using OperandId = uint32_t;

constexpr OperandId invalidOperandId = static_cast<OperandId>(-1);

// WHICH COEFFICIENT OF A FITTED MODEL, by its position in the model's covariate list. It lives here
// rather than with the tree's own ids because the BLOCK owns parameters -- the tree merely has a
// leaf that names one, exactly as it has a leaf naming a `ColumnId` the data layer owns.
using ParamId = uint32_t;

enum class OperandShape : uint8_t
{
    Scalar = 1,
    Vector = 2,
};

// An operand is either a scalar of one of the veil base types, or the `double` time vector -- which
// is the only vectorised type. The base type is fixed for the life of the operand.
class Operand final
{
public:
    const TypeFull type;
    const OperandShape shape;

private:
    Operand(TypeFull type, OperandShape shape)
        : type(type), shape(shape) {}

public:
    static Operand createScalar(TypeFull type) { return Operand(type, OperandShape::Scalar); }

    // The time vector is always `double`; there is no other vectorised base type.
    static Operand createVector() { return Operand(TypeFull::createDouble(), OperandShape::Vector); }

    bool isVector() const noexcept { return this->shape == OperandShape::Vector; }
    bool isScalar() const noexcept { return this->shape == OperandShape::Scalar; }

    bool operator==(const Operand& other) const = default;
};

} // namespace veil
