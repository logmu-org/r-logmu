// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>
#include "veil/Block.hpp"
#include "veil/Instruction.hpp"
#include "veil/Operand.hpp"

namespace veil
{

// Where each time-vector operand is written and where it is last read.
//
// ANNOTATES NOTHING AND REWRITES NOTHING. It answers one question about a block so that a separate
// pass can act on the answer, which is how the interval propagation is arranged too.
//
// WHY THIS IS A SCAN AND NOT A DATAFLOW ANALYSIS. A block is straight-line code -- no jumps, no
// loops, no phi nodes -- and every operand an instruction assigns is assigned exactly once. So an
// operand is live from the instruction that writes it to the last instruction that reads it, and
// that whole range is one contiguous interval discovered by looking at each instruction once. This
// is the reason the IL is a linear three-address form rather than a tree.
//
// SCALARS ARE DELIBERATELY NOT ANALYSED. A scalar register is eight bytes or so and there are a few
// dozen of them; a time-vector buffer is one double per slot per operand, and the whole point of
// working out the ranges is to stop a block needing one of those for every temporary it ever names.

// The live range of one time-vector operand, in body instruction indices.
//
// `lastUsedAt` equals `definedAt` for a vector nothing reads -- the two are not distinguishable
// here, and it is dead-code elimination's business to notice, not this pass's.
struct VectorLiveRange final
{
    size_t definedAt = 0;
    size_t lastUsedAt = 0;
    bool defined = false; // False for every scalar operand.
};

// Works out the live range of every time-vector operand in `block`, indexed by OperandId.
//
// The guards here are the ones a mis-ordered or mis-built block would trip, and each is checked
// rather than assumed because the whole of buffer allocation rests on them:
//
//   - a vector operand never appears in the PROLOGUE, which runs before the grid exists and so has
//     no slots to write into;
//   - a vector operand is never a block OUTPUT, an individual's answer being a single value;
//   - a vector operand is never bound to a constant or a column, so nothing writes it but an
//     instruction;
//   - each is assigned exactly ONCE, which is what makes its live range a single interval;
//   - each is READ only after it is written, which is the fault a pass reordering instructions
//     would introduce.
inline std::vector<VectorLiveRange> passComputeLiveness(const Block& block)
{
    const size_t operandCount = block.operandCount();
    std::vector<VectorLiveRange> ranges(operandCount);

    const auto isVector = [&block](OperandId id) { return block.operandAt(id).isVector(); };

    for (const Instruction& instruction : block.prologue())
    {
        if (instruction.result != invalidOperandId && isVector(instruction.result))
        {
            throw std::runtime_error("veil: a time vector was written in the prologue, which runs "
                                     "before the time grid is built.");
        }
        for (uint8_t i = 0; i < instruction.argCount; ++i)
        {
            if (isVector(instruction.args[i]))
            {
                throw std::runtime_error("veil: the prologue reads a time vector, which does not "
                                         "exist until the grid is built.");
            }
        }
    }

    for (const OperandId output : block.outputs())
    {
        if (isVector(output))
        {
            throw std::runtime_error("veil: a block output varies over time, so it is not one value "
                                     "per individual.");
        }
    }

    for (const ConstantBinding& constant : block.constants())
    {
        if (isVector(constant.operand))
        {
            throw std::runtime_error("veil: a time vector cannot hold a constant binding.");
        }
    }
    for (const ColumnBinding& column : block.columns())
    {
        if (isVector(column.operand))
        {
            throw std::runtime_error("veil: a time vector cannot be bound to a column.");
        }
    }

    const std::vector<Instruction>& body = block.body();
    for (size_t index = 0; index < body.size(); ++index)
    {
        const Instruction& instruction = body[index];

        // Arguments first: an operand read by the very instruction that writes it would be a
        // self-reference, which single assignment rules out, and taking the reads first is what
        // makes that show up as a use before definition rather than passing unnoticed.
        for (uint8_t i = 0; i < instruction.argCount; ++i)
        {
            const OperandId argument = instruction.args[i];
            if (!isVector(argument)) { continue; }

            VectorLiveRange& range = ranges.at(argument);
            if (!range.defined)
            {
                throw std::runtime_error("veil: a time vector is read before anything writes it.");
            }
            range.lastUsedAt = index;
        }

        const OperandId result = instruction.result;
        if (result == invalidOperandId || !isVector(result)) { continue; }

        VectorLiveRange& range = ranges.at(result);
        if (range.defined)
        {
            throw std::runtime_error("veil: a time vector is written twice, so its live range is "
                                     "not one interval.");
        }
        range.defined = true;
        range.definedAt = index;
        range.lastUsedAt = index;
    }

    return ranges;
}

} // namespace veil
