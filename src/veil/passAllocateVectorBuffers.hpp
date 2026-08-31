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
#include "veil/passComputeLiveness.hpp"

namespace veil
{

// Gives every time-vector operand a physical buffer, reusing one whose last reader has already run.
//
// LINEAR SCAN, WHICH IS ALL A STRAIGHT-LINE BLOCK NEEDS. Register allocation is graph colouring in
// general because control flow makes live ranges into arbitrary sets. A veil block has no control
// flow, so each range is one interval on one axis, and walking the instructions in order while
// keeping a free list assigns the minimum number of buffers with no search at all.
//
// WHAT IT IS WORTH. An AEV's inner loop holds a handful of vectors at once but names more than that
// over its length, and the ones it names are the ones the interpreter was allocating. The prize
// grows with the output count rather than with the expression: a log-likelihood Hessian has
// p(p+1)/2 outputs over one grid, so at p = 10 the virtual form wants fifty-odd buffers where the
// live count at any point is a few. Keeping the working set in L1 is the whole game there.
//
// THE ONE SUBTLETY, AND IT IS THE ALIASING QUESTION. A buffer is released only once the instruction
// that last read it has FINISHED, so a result never lands in a buffer one of its own arguments is
// still using. Letting a result overwrite an argument would work for the elementwise operations --
// slot i is read before slot i is written -- and would break for anything that reads a slot other
// than the one it writes. Rather than depend on which ops are elementwise, and on that never
// changing, the allocator gives up at most one buffer and the question does not arise. Note this
// costs nothing asymptotically: a chain still ping-pongs between two buffers.

// Assigns buffers and records the assignment on the block. Returns how many buffers it needed.
//
// `ranges` must be this block's own liveness, as passComputeLiveness computes it. The two are kept
// apart because the analysis answers a question about the block and this decides what to do with
// the answer, which is the same division as interval propagation and the fold that consumes it.
inline size_t passAllocateVectorBuffers(Block& block, const std::vector<VectorLiveRange>& ranges)
{
    const size_t operandCount = block.operandCount();
    if (ranges.size() != operandCount)
    {
        throw std::runtime_error("veil: the liveness handed to buffer allocation is not this "
                                 "block's own.");
    }

    std::vector<BufferId> bufferOf(operandCount, invalidBufferId);

    // For each buffer, the first instruction index at which it may be handed to something else.
    std::vector<size_t> freeFrom;

    const std::vector<Instruction>& body = block.body();
    for (size_t index = 0; index < body.size(); ++index)
    {
        const OperandId result = body[index].result;
        if (result == invalidOperandId || !block.operandAt(result).isVector()) { continue; }

        const VectorLiveRange& range = ranges[result];
        if (!range.defined || range.definedAt != index)
        {
            throw std::runtime_error("veil: the liveness disagrees with the body about where a time "
                                     "vector is written.");
        }

        // Lowest free buffer first, so the assignment is a function of the block alone and a test
        // can state what it should be rather than only that it is small enough.
        BufferId chosen = invalidBufferId;
        for (size_t candidate = 0; candidate < freeFrom.size(); ++candidate)
        {
            if (freeFrom[candidate] <= index)
            {
                chosen = static_cast<BufferId>(candidate);
                break;
            }
        }
        if (chosen == invalidBufferId)
        {
            chosen = static_cast<BufferId>(freeFrom.size());
            freeFrom.push_back(0);
        }

        bufferOf[result] = chosen;
        freeFrom[chosen] = range.lastUsedAt + 1;
    }

    const size_t count = freeFrom.size();
    block.assignVectorBuffers(std::move(bufferOf), count);
    return count;
}

} // namespace veil
