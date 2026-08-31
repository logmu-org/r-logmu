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
#include "veil/Op.hpp"
#include "veil/Operand.hpp"

namespace veil
{

// Works out, for every time-vector operand, which slots of the grid are ever read from it.
//
// ANNOTATES NOTHING. It answers the question and a separate pass records the answer on the block,
// the same division as liveness and buffer allocation.
//
// WHY THERE IS ANYTHING TO WORK OUT. The grid is not one thing. `integrate` reads the midpoint
// slots and stops before the death slot; `died_value` reads the death slot and nothing else. So a
// chain feeding only one of the two is computed where the other looks for no reason at all. For an
// A/E/V that chain is A's whole integrand, and the saving is every slot but one -- or every slot
// including that one, for an individual who did not die.
//
// THE ANALYSIS IS ONE BACKWARD PASS, for the same reason liveness is one forward pass: the block is
// straight-line single-assignment code, so an operand's readers all sit after the instruction that
// writes it, and walking from the end means every reader has been seen by the time its argument is
// reached.

namespace detail
{

// Whether slot i of an op's result depends only on slot i of its arguments.
//
// EVERY OP VEIL HAS TODAY IS SLOT-ALIGNED EXCEPT THE TWO FINALISING ONES, and that is what makes
// the demand propagate straight through an expression: if only the death slot of a product is
// wanted, only the death slot of each factor is wanted.
//
// THE OP THAT WILL BREAK THIS IS ALREADY NAMED IN THE SPECIFICATION. `accumulate`, the integral
// from nu up to t, is deferred to liability valuation, and its slot i reads every slot before i --
// so wanting one slot of its result means wanting all the earlier slots of its argument. A gather,
// a shift or a running total is the same shape. Anything of that kind must be handled by name here
// rather than falling into the general case, which is why this asks the question explicitly instead
// of assuming.
inline bool isSlotAligned(Op op)
{
    switch (opInfo(op).category)
    {
        case OpCategory::Finalise: return false;
        case OpCategory::Arithmetic:
        case OpCategory::Comparison:
        case OpCategory::Logical:
        case OpCategory::Selection:
        case OpCategory::Conversion:
        case OpCategory::TimeVector: return true;
    }
    return false;
}

} // namespace detail

// The slots each operand is read at, indexed by OperandId. A scalar operand's entry is meaningless
// and left empty.
inline std::vector<SlotDemand> passComputeSlotDemand(const Block& block)
{
    const size_t operandCount = block.operandCount();
    std::vector<SlotDemand> demand(operandCount);

    const auto isVector = [&block](OperandId id) { return block.operandAt(id).isVector(); };

    const std::vector<Instruction>& body = block.body();
    for (size_t remaining = body.size(); remaining-- > 0;)
    {
        const Instruction& instruction = body[remaining];
        const bool resultIsVector =
            instruction.result != invalidOperandId && isVector(instruction.result);

        if (!detail::isSlotAligned(instruction.op))
        {
            // A finalising op collapses the vector to one number, so it is where demand STARTS
            // rather than something demand passes through.
            if (resultIsVector)
            {
                throw std::runtime_error("veil: an op that collapses the time vector was given a "
                                         "time vector to write into.");
            }
            if (instruction.argCount == 0 || !isVector(instruction.args[0]))
            {
                throw std::runtime_error("veil: an op that collapses the time vector was not given "
                                         "one to read.");
            }

            SlotDemand& argument = demand.at(instruction.args[0]);
            if (instruction.op == Op::Integrate) { argument.midpoints = true; }
            else if (instruction.op == Op::DiedValue) { argument.death = true; }
            else
            {
                throw std::runtime_error("veil: this op collapses the time vector but does not say "
                                         "which slots it reads.");
            }
            continue;
        }

        if (!resultIsVector)
        {
            // A scalar instruction reading a time vector would be reading one value out of many
            // without saying which, and only the finalising ops are allowed to do that.
            for (uint8_t i = 0; i < instruction.argCount; ++i)
            {
                if (isVector(instruction.args[i]))
                {
                    throw std::runtime_error("veil: a scalar instruction reads a time vector.");
                }
            }
            continue;
        }

        const SlotDemand wanted = demand.at(instruction.result);
        if (!wanted.any()) { continue; } // Nothing reads this, so nothing it reads is needed for it.

        for (uint8_t i = 0; i < instruction.argCount; ++i)
        {
            const OperandId argument = instruction.args[i];
            if (isVector(argument)) { demand.at(argument).merge(wanted); }
        }
    }

    return demand;
}

} // namespace veil
