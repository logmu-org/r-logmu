// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include "veil/Op.hpp"
#include "veil/Operand.hpp"

namespace veil
{

// No veil op takes more than three arguments. The op table's own well-formedness check refuses a row
// with a higher arity, so this bound cannot drift away from the table.
constexpr size_t MaxArgCount = 3;

// An instruction's `parameter` when it has none, which is all but a handful of them.
constexpr uint32_t invalidParameter = static_cast<uint32_t>(-1);

// One three-address instruction: an op moniker, the operands it reads, and the single operand it
// assigns. That is the whole of veil's executable form -- there is no nesting here, and an operand is
// referred to by index into the block that declared it rather than by a pointer to a sub-expression.
//
// `argCount` is the op's arity, so the unused tail of `args` is never read. It is carried rather than
// looked up from the op table on every access because an interpreter reads it once per instruction.
//
// `parameter` is a COMPILE-TIME parameter of the instruction: something the op needs that is not a
// value flowing through a register, and so cannot be an operand. Today the only op that carries one
// is `vector_log_mu`, where it is the TableId of the mortality table the instruction reads -- a table
// is fixed when the block is built, is far too big to sit in a register, and is shared by every
// individual the block runs over. Everything else leaves it invalid, and an op that does not declare
// a parameter must never read one.
struct Instruction final
{
    Op op = Op::Pos;
    OperandId result = invalidOperandId;
    std::array<OperandId, MaxArgCount> args = {invalidOperandId, invalidOperandId, invalidOperandId};
    uint8_t argCount = 0;
    uint32_t parameter = invalidParameter;
};

// Builds an instruction from the operands it reads, in order. Written this way so that a caller
// cannot silently disagree with `argCount` about how many of `args` are meaningful.
inline Instruction makeInstruction(Op op, OperandId result, OperandId a)
{
    return Instruction{op, result, {a, invalidOperandId, invalidOperandId}, 1};
}

inline Instruction makeInstruction(Op op, OperandId result, OperandId a, OperandId b)
{
    return Instruction{op, result, {a, b, invalidOperandId}, 2};
}

inline Instruction makeInstruction(Op op, OperandId result, OperandId a, OperandId b, OperandId c)
{
    return Instruction{op, result, {a, b, c}, 3};
}

// An instruction that reads one operand and one compile-time parameter. Kept separate from the
// overloads above, rather than given a defaulted argument, so that a parameter can never be passed
// where an operand was meant -- both are a bare uint32_t, and the compiler would not object.
inline Instruction makeParameterisedInstruction(
    Op op,
    OperandId result,
    uint32_t parameter,
    OperandId a)
{
    return Instruction{op, result, {a, invalidOperandId, invalidOperandId}, 1, parameter};
}

} // namespace veil
