// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>
#include "veil/ColumnView.hpp" // ColumnId
#include "veil/Instruction.hpp"
#include "veil/MortalityTable.hpp"
#include "veil/Operand.hpp"
#include "veil/ScalarValue.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// The physical store a time-vector operand computes into: one double per slot of the individual's
// exposure. Several operands share one buffer whenever their live ranges do not overlap, which is
// what stops a block needing a buffer for every temporary it names.
//
// A BufferId is not an OperandId and the two must never be swapped. An OperandId says WHAT a value
// is, and stays fixed for the life of the block; a BufferId says WHERE it sits while it is needed,
// and is a decision the allocation pass makes.
using BufferId = uint32_t;

constexpr BufferId invalidBufferId = static_cast<BufferId>(-1);

// Which parts of the time vector a value is actually read at.
//
// The grid has two halves and they are read by different operators: `integrate` reads the midpoint
// slots and never the death slot, `died_value` reads the death slot and never the midpoints. So a
// value feeding only one of them needs computing only where that one looks.
//
// THE CASE THIS EXISTS FOR is the A of an A/E/V. Its integrand reaches `died_value` alone, so of the
// twelve to twenty slots an exposure is cut into, exactly one is ever read -- and for a survivor,
// none at all, since there is no death slot to read. Survivors are the great majority of any
// portfolio, so the whole A contribution disappears for most individuals without a branch being
// written anywhere.
//
// Nothing demanded at all means the value is never read. That is dead code, and it costs nothing to
// skip because the demand already says so.
struct SlotDemand final
{
    bool midpoints = false;
    bool death = false;

    bool any() const noexcept { return this->midpoints || this->death; }

    void merge(const SlotDemand& other) noexcept
    {
        this->midpoints = this->midpoints || other.midpoints;
        this->death = this->death || other.death;
    }

    bool operator==(const SlotDemand& other) const = default;
};

// What an unanalysed block assumes: every slot of every vector might be read.
constexpr SlotDemand everySlot = SlotDemand{true, true};

// A constant, written into its operand once before any individual is processed.
struct ConstantBinding final
{
    OperandId operand = invalidOperandId;
    ScalarValue value = 0.0;
};

// An operand the host refills from a data column, once for each individual in turn.
struct ColumnBinding final
{
    OperandId operand = invalidOperandId;
    ColumnId column = 0;
};

// Where a block's time vector comes from. The three operands are per-individual scalars -- the
// exposure start, its end, and whether the individual died -- so the grid is rebuilt per individual
// from values already in the register file, with no separate loading path.
//
// WITHOUT AN INCLUDE they are the three exposure columns as loaded. WITH ONE they are temporaries the
// prologue computed, holding the exposure after it has been clipped to the include's interval, and
// `died` is the RECOMPUTED flag rather than the column. Nothing here can tell the two cases apart,
// which is the point: the grid is built from whatever those registers hold.
//
// `included` is the answer to "does this individual contribute at all" -- every gate held, no offset
// was missing, and the clipped interval is not empty. Invalid means nothing could have made it false,
// so the individual is always included. When it is false the grid is EMPTY rather than built, which
// is what makes an excluded individual integrate to zero rather than being an error.
//
// `deltaTClicks` is fixed for the whole calculation, which is what lets the number of slots follow
// from the exposure alone.
struct TimeGridBinding final
{
    OperandId start = invalidOperandId;
    OperandId end = invalidOperandId;
    OperandId died = invalidOperandId;
    int deltaTClicks = 0;
    OperandId included = invalidOperandId;
};

// One veil block: a table of typed operands, the values that flow into them from outside, the
// three-address body, and the operands the host reads back when the body has run.
//
// OPERANDS ARE BLOCK-LOCAL. An OperandId indexes this block's own table and means nothing against
// any other block, or against another concurrent invocation of this one. That is what lets a worker
// thread, or a host-driven solve iterating on a parameter, hold its own register file without any
// shared operand table to collide over.
//
// The three sorts of operand differ only in where their value comes from. A constant is written once;
// a column-bound operand is rewritten per individual by the host, which is what makes the body a
// function of the individual without the body itself knowing how a column is stored; everything else
// is a temporary an instruction assigns. Nothing here distinguishes a scalar from a time vector --
// that is the Operand's own shape -- so this same structure carries the vector phase unchanged.
class Block final
{
public:
    OperandId addConstant(const TypeFull& type, ScalarValue value)
    {
        const OperandId id = this->addOperand(type);
        this->constantList.push_back(ConstantBinding{id, value});
        return id;
    }

    OperandId addColumn(const TypeFull& type, ColumnId column)
    {
        const OperandId id = this->addOperand(type);
        this->columnList.push_back(ColumnBinding{id, column});
        return id;
    }

    OperandId addTemp(const TypeFull& type) { return this->addOperand(type); }

    // A time-varying temporary: one value per slot of the individual's time vector. Always double,
    // whatever the tree said the expression's type was -- a datey that varies over time reaches a
    // vector as a number of years, because there is no vectorised integer type and no need for one.
    OperandId addVectorTemp()
    {
        this->operandList.push_back(Operand::createVector());
        return static_cast<OperandId>(this->operandList.size() - 1);
    }

    // Takes a copy of a mortality table for an instruction to read, and returns how to name it.
    //
    // The block OWNS its tables rather than pointing at whatever the tree was built from. A block is
    // meant to be a self-contained executable artefact -- that is the whole point of the TAC being
    // the language-independence boundary -- and a borrowed table would tie its lifetime to an
    // ObjStore that the host is free to drop. The copy costs a few tens of kilobytes once when the
    // block is built, against a table read for every slot of every individual.
    TableId addTable(MortalityTable table)
    {
        this->tableList.push_back(std::move(table));
        return static_cast<TableId>(this->tableList.size() - 1);
    }

    // Bounds-checked, for the same reason operandAt is: a table id that is not this block's own is a
    // lowering fault, and reading past the list is how it would surface as a wrong log mu instead.
    const MortalityTable& tableAt(TableId id) const { return this->tableList.at(id); }

    size_t tableCount() const noexcept { return this->tableList.size(); }

    // Declares that this block samples time, and says where the exposure comes from: three scalar
    // operands already bound to the start, end and died columns. The host builds the grid from them
    // once per individual, before the body runs.
    void useTimeGrid(
        OperandId start,
        OperandId end,
        OperandId died,
        int deltaTClicks,
        OperandId included = invalidOperandId)
    {
        this->grid = TimeGridBinding{start, end, died, deltaTClicks, included};
    }

    const std::optional<TimeGridBinding>& timeGrid() const noexcept { return this->grid; }

    void emit(const Instruction& instruction) { this->instructionList.push_back(instruction); }

    // The prologue runs once per individual, after the columns are loaded and BEFORE the time grid is
    // built. It exists because the include has to resolve first: it is what decides where the
    // exposure starts and ends, so it cannot be part of a body that only runs once the grid says how
    // many slots there are. Only the include lowers into it; everything else is body.
    void emitPrologue(const Instruction& instruction) { this->prologueList.push_back(instruction); }

    const std::vector<Instruction>& prologue() const noexcept { return this->prologueList; }

    void addOutput(OperandId operand) { this->outputList.push_back(operand); }

    // Records which physical buffer each time-vector operand computes into. This is the block's
    // ALLOCATED form: the body is unchanged and still names operands, but the host now knows how
    // few buffers it has to find. Until this is called the block is in its VIRTUAL form, where
    // every vector operand has a buffer to itself -- which is what the accessors below answer, so
    // a block that never goes through allocation still runs, it just costs more memory.
    //
    // Validated here rather than trusted, because a mapping that is wrong by one produces a
    // plausible wrong number rather than a failure: a buffer holding the wrong operand's values is
    // still a buffer full of doubles.
    void assignVectorBuffers(std::vector<BufferId> mapping, size_t count)
    {
        if (mapping.size() != this->operandList.size())
        {
            throw std::runtime_error("veil: a buffer assignment does not cover every operand.");
        }
        for (size_t id = 0; id < mapping.size(); ++id)
        {
            const bool isVector = this->operandList[id].isVector();
            if (isVector && mapping[id] >= count)
            {
                throw std::runtime_error("veil: a time vector was assigned no buffer.");
            }
            if (!isVector && mapping[id] != invalidBufferId)
            {
                throw std::runtime_error("veil: a scalar operand was assigned a time-vector buffer.");
            }
        }
        this->vectorBufferList = std::move(mapping);
        this->vectorBufferTotal = count;
    }

    // The buffer an operand's time vector lives in. Before allocation has run this is the operand's
    // own id, so the two forms differ in how much memory they need and in nothing else.
    BufferId vectorBufferOf(OperandId id) const
    {
        if (this->vectorBufferList.empty()) { return id; }
        return this->vectorBufferList.at(id);
    }

    size_t vectorBufferCount() const noexcept
    {
        return this->vectorBufferList.empty() ? this->operandList.size() : this->vectorBufferTotal;
    }

    // Records where each time-vector operand is actually read, so the host computes it nowhere else.
    // Until this is called every vector is taken to be read everywhere, which is the safe answer and
    // the one the block behaved by before the analysis existed.
    void assignSlotDemand(std::vector<SlotDemand> demand)
    {
        if (demand.size() != this->operandList.size())
        {
            throw std::runtime_error("veil: a slot demand does not cover every operand.");
        }
        this->slotDemandList = std::move(demand);
    }

    SlotDemand slotDemandOf(OperandId id) const
    {
        if (this->slotDemandList.empty()) { return everySlot; }
        return this->slotDemandList.at(id);
    }

    // Bounds-checked: an operand id that is not this block's own is a lowering fault, and reading
    // past the table is how it would otherwise show up as a wrong answer rather than an error.
    const Operand& operandAt(OperandId id) const { return this->operandList.at(id); }

    size_t operandCount() const noexcept { return this->operandList.size(); }

    const std::vector<ConstantBinding>& constants() const noexcept { return this->constantList; }
    const std::vector<ColumnBinding>& columns() const noexcept { return this->columnList; }
    const std::vector<Instruction>& body() const noexcept { return this->instructionList; }
    const std::vector<OperandId>& outputs() const noexcept { return this->outputList; }

private:
    OperandId addOperand(const TypeFull& type)
    {
        this->operandList.push_back(Operand::createScalar(type));
        return static_cast<OperandId>(this->operandList.size() - 1);
    }

    std::vector<Operand> operandList;
    std::vector<MortalityTable> tableList;
    std::vector<ConstantBinding> constantList;
    std::vector<ColumnBinding> columnList;
    std::vector<Instruction> prologueList;
    std::vector<Instruction> instructionList;
    std::vector<OperandId> outputList;
    std::optional<TimeGridBinding> grid;

    // Both empty until their passes have run, which is what the virtual form is.
    std::vector<BufferId> vectorBufferList;
    size_t vectorBufferTotal = 0;
    std::vector<SlotDemand> slotDemandList;
};

} // namespace veil
