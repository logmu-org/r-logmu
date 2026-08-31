// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <vector>
#include "veil/Block.hpp"
#include "veil/Instruction.hpp"
#include "veil/Operand.hpp"

namespace veil
{

// Finds instructions that write a result nothing ever reads.
//
// IT FINDS NOTHING TODAY, AND THAT IS THE POINT OF RUNNING IT. Lowering is demand-driven -- it walks
// down from the roots and emits an instruction only because something asked for its value -- so a
// block as lowered contains no dead code at all. Measured across the whole test suite rather than
// argued: not one instruction in any block any entry point built.
//
// SO WHY HAVE IT. Because that property is an invariant of lowering rather than a law of nature, and
// the first TAC-level rewrite will break it deliberately. Fusing E's and V's reductions over one
// grid, or giving if/else arms that omit their instructions, both leave results behind that nothing
// reads. This is the pass that notices, and a test asserting the count is zero is what will report
// the day it stops being.
//
// REMOVAL IS NOT BUILT, on purpose. A pass that deletes instructions and has never had one to delete
// is untested code holding a knife; the analysis is the half that can be checked today. Whoever adds
// the first rewrite that leaves something dead should take the indices below and erase them, and
// they will have a case to test it against.
//
// WHAT COUNTS AS READ. A block's outputs, and the three or four operands the time grid is built from
// -- the exposure bounds, the death flag, and whether the individual is included at all. Those last
// are read by the host between the prologue and the body rather than by any instruction, so a walk
// that only looked at instruction arguments would call the whole include clip dead.

struct DeadInstructions final
{
    std::vector<size_t> body;
    std::vector<size_t> prologue;

    size_t count() const noexcept { return this->body.size() + this->prologue.size(); }
};

inline DeadInstructions passFindDeadInstructions(const Block& block)
{
    std::vector<char> read(block.operandCount(), 0);

    for (const OperandId output : block.outputs()) { read.at(output) = 1; }

    if (const std::optional<TimeGridBinding>& grid = block.timeGrid())
    {
        read.at(grid->start) = 1;
        read.at(grid->end) = 1;
        read.at(grid->died) = 1;
        if (grid->included != invalidOperandId) { read.at(grid->included) = 1; }
    }

    DeadInstructions dead;

    // Backwards, and the BODY BEFORE THE PROLOGUE, because that is the order they run in: a prologue
    // instruction can only be needed by something later, which includes the whole body.
    const auto sweep = [&read](const std::vector<Instruction>& list, std::vector<size_t>& found)
    {
        for (size_t remaining = list.size(); remaining-- > 0;)
        {
            const Instruction& instruction = list[remaining];
            if (instruction.result != invalidOperandId && read.at(instruction.result) == 0)
            {
                found.push_back(remaining);
                continue;
            }
            for (uint8_t i = 0; i < instruction.argCount; ++i) { read.at(instruction.args[i]) = 1; }
        }
    };

    sweep(block.body(), dead.body);
    sweep(block.prologue(), dead.prologue);
    return dead;
}

} // namespace veil
