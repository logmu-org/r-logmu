// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>
#include "veil/Block.hpp"
#include "veil/ColumnView.hpp" // ColumnId
#include "veil/Instruction.hpp"
#include "veil/MortalityTable.hpp"
#include "veil/Node.hpp"
#include "veil/ObjStore.hpp"
#include "veil/Op.hpp"
#include "veil/Operand.hpp"
#include "veil/ScalarValue.hpp"
#include "veil/TimeGrid.hpp"
#include "veil/Tree.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// Lowering: the tree becomes a block of three-address instructions.
//
// This is the step that ends the value-level phase. Everything above it asks what to compute and is
// natural on a tree, where a rewrite matches the shape of a sub-expression; everything below it asks
// how to lay the computation out and is natural on a linear body, where liveness is a backward scan.
//
// TRAVERSAL IS POST-ORDER FROM THE ROOT, NOT AN ID SWEEP. The coercion pass appends its ToDouble
// nodes at the end of the arena and repoints the parents, so a node can hold an argument with a
// larger id than its own; only the passes that run before coercion may assume children come first.
// The memo doubles as common-sub-expression sharing: a node reached twice, as `%in%` desugaring and
// CSE both produce, is lowered once and its operand reused.
//
// SCALAR VERSUS VECTOR IS DECIDED BY THE TIME-VARYING TAGS, NOT BY THE TYPES. A node the tagging
// pass marked time-varying becomes a vector operand; everything else is a per-individual scalar.
//
// A TIME-VARYING VALUE IS ALWAYS A DOUBLE VECTOR HOLDING YEARS. There is no vectorised integer type
// and no need for one, so a datey or durationy that varies over time reaches its operand already
// converted. That is why `.t` lowers straight to vector_t rather than to clicks followed by a
// conversion, and why `.t - .b` lowers to vector_durn, which does the subtraction in exact clicks
// and converts once. It is also why ToDouble applied to something already vectorised emits nothing:
// the conversion happened at the source. Tim's specification names this wrinkle directly -- datey
// and durationy arithmetic has to be resolved while it is still scalar.
//
// WHAT IT ASSUMES, AND CHECKS FOR ITSELF. Coercion insertion has already made every arithmetic
// widening explicit, so an arithmetic instruction's operands all share the result's storage. The one
// mixture lowering has to settle itself is a comparison, which the coercion pass deliberately leaves
// alone: `.i$birth == 1990.5` can still arrive as a datey against a plain number when narrowing
// declined it. Rather than let the interpreter guess, lowering inserts the ToDouble -- reading the
// clicks as years, per the datey rules -- so that by the time an instruction exists its operands
// agree. Everything else is asserted rather than assumed: a mismatch throws here, where the tree is
// still readable, instead of surfacing as a wrong number.

// Where a block that samples time gets its exposure from, and how finely it samples. Supplied by the
// caller because the core does not know what a column is called; the binding resolves the names.
struct ExposureColumns final
{
    ColumnId start = 0;
    ColumnId end = 0;
    ColumnId died = 0;
    int deltaTClicks = 0;
};

namespace detail
{

// How an operand is held at run time. Two operands of the same class can be handed to one
// instruction; two of different classes cannot, and lowering must have reconciled them first.
enum class StorageClass
{
    Number,  // double
    Click,   // datey, durationy -- integer clicks on the annual grid
    Logical, // bool
    Code,    // category, text -- an integer code, which is not a number and does not convert to one
};

inline StorageClass storageOf(const TypeFull& type)
{
    switch (type.type)
    {
        case Type::Double: return StorageClass::Number;
        case Type::Datey:
        case Type::Durationy: return StorageClass::Click;
        case Type::Bool: return StorageClass::Logical;
        case Type::Category:
        case Type::Text: return StorageClass::Code;
        case Type::DateyInterval: break;
    }
    throw std::runtime_error("veil: a datey_interval has no operand storage.");
}

inline std::string opName(Op op) { return std::string(moniker(op)); }

// The value of a literal node, in the form its type says the operand holds. A text literal has
// nowhere to go: passEncodeText has already turned every one of them into a category index or folded
// the comparison it sat in, so one reaching here is a pass that stopped firing rather than a value.
inline ScalarValue literalValue(const LitPayload& lit)
{
    const StorageClass storage = storageOf(lit.type);

    if (const auto* b = std::get_if<bool>(&lit.value))
    {
        if (storage == StorageClass::Logical) { return ScalarValue(*b); }
    }
    else if (const auto* d = std::get_if<double>(&lit.value))
    {
        if (storage == StorageClass::Number) { return ScalarValue(*d); }
    }
    else if (const auto* i = std::get_if<int>(&lit.value))
    {
        if (storage == StorageClass::Click || storage == StorageClass::Code) { return ScalarValue(*i); }
    }
    else
    {
        throw std::runtime_error("veil: a text literal reached lowering unencoded.");
    }

    throw std::runtime_error("veil: a literal's value does not match the type it carries.");
}

// Lowers one tree into one block. Held as a struct so the recursion carries the memo and the block
// without threading them through every call.
struct Lowerer final
{
    const Tree& tree;
    const ObjStore& objs;
    const std::vector<char>& timeVarying;
    std::optional<ExposureColumns> exposure;
    std::optional<ObjId> includeObj;
    Block block;
    std::vector<OperandId> operandOfNode;

    // Set while the include is being lowered, which is the only thing that belongs in the prologue.
    bool intoPrologue = false;

    // Set while the include is being resolved, so that anything inside it reaching back for the time
    // grid is caught as the circularity it is rather than recursing until the stack runs out.
    bool resolvingInclude = false;

    // Whether anything can empty this individual, and where that answer is held. Empty when nothing
    // can, which is the case for every block that has no include.
    std::optional<OperandId> includedFlag;

    Lowerer(
        const Tree& tree,
        const ObjStore& objs,
        const std::vector<char>& timeVarying,
        const std::optional<ExposureColumns>& exposure,
        std::optional<ObjId> includeObj)
        : tree(tree), objs(objs), timeVarying(timeVarying), exposure(exposure),
            includeObj(includeObj), operandOfNode(tree.size(), invalidOperandId) {}

    // Every instruction goes through here so that a gate's own expression, which lowers by the same
    // ordinary paths as anything else, lands in the prologue along with the clip that reads it.
    void emit(const Instruction& instruction)
    {
        if (this->intoPrologue) { this->block.emitPrologue(instruction); }
        else { this->block.emit(instruction); }
    }

    bool isVector(OperandId operand) const { return this->block.operandAt(operand).isVector(); }

    // A scalar temporary written by one instruction, which is how every step of the clip is held.
    OperandId emitScalar(Op op, const TypeFull& type, OperandId a, OperandId b)
    {
        const OperandId result = this->block.addTemp(type);
        this->emit(makeInstruction(op, result, a, b));
        return result;
    }

    // An existing operand already bound to `column`, so two references to the same field read it
    // once. Distinct field nodes for one column are common: `.x` expands to `.t - .b` wherever it
    // appears, and the memo alone would not catch that.
    OperandId columnOperand(ColumnId column, const TypeFull& type)
    {
        for (const ColumnBinding& binding : this->block.columns())
        {
            if (binding.column == column)
            {
                // ONE COLUMN, ONE TYPE. An assertion on an invariant, NOT a guard against anything a
                // user can currently write, and worth being exact about which.
                //
                // Every request now takes its type from `columnTypes` and from nothing else. The
                // exposure columns are resolved by `exposureColumn`, which refuses a column of the
                // wrong type before it hands the id over; a field node carries the type
                // passAnnotateTypes read out of that same vector. Two requests for one column
                // therefore cannot disagree, so no input reaches this throw.
                //
                // IT USED TO BE REACHABLE, and what changed is instructive. An include's offset term
                // held a column id and asked for it as a datey outright, without consulting what the
                // column held -- so naming a durationy column tripped exactly this. Offsets are
                // expressions now, checked on the lowered operand, and that caller is gone.
                //
                // KEPT ANYWAY, because being wrong here is silent rather than loud: a durationy
                // operand would satisfy every check downstream, both being clicks, and loadRecord
                // validates a column against its FIRST binding rather than against this request. It
                // costs one comparison per column per block.
                //
                // NO TEST COVERS IT AND NONE CAN BE WRITTEN while the invariant above holds. Do not
                // invent one that appears to; it would be testing a path that does not exist. What
                // would make it reachable again is a new caller asserting a type against a column id
                // independently, which is the shape the offset term had.
                if (!(this->block.operandAt(binding.operand).type == type))
                {
                    throw std::runtime_error("veil: one column is referred to as two different types.");
                }
                return binding.operand;
            }
        }
        return this->block.addColumn(type, column);
    }

    // Conjoins one more condition onto the running "is this individual included at all" flag,
    // starting it off if this is the first. Kept as a running temp rather than a fold over a list so
    // that a block with no gates and no offsets carries no flag at all.
    void conjoinIncluded(std::optional<OperandId>& included, OperandId condition)
    {
        included = included.has_value()
            ? this->emitScalar(Op::And, TypeFull::createBool(), *included, condition)
            : condition;
    }

    // One interval term of the include, narrowing the running bounds.
    //
    // The bounds arrive as CLICKS on the same grid the exposure uses, so the whole clip is exact
    // integer arithmetic and `max`/`min` on clicks need none of the NaN and signed-zero care their
    // floating-point counterparts do.
    void lowerIntervalTerm(
        const IncludeTerm& term,
        OperandId& nu,
        OperandId& tau,
        std::optional<OperandId>& included)
    {
        const TypeFull datey = TypeFull::createDatey();

        if (const auto* absolute = std::get_if<AbsoluteTerm>(&term))
        {
            const OperandId from = this->block.addConstant(datey, ScalarValue(absolute->fromClicks));
            const OperandId to = this->block.addConstant(datey, ScalarValue(absolute->toClicks));
            nu = this->emitScalar(Op::Max, datey, nu, from);
            tau = this->emitScalar(Op::Min, datey, tau, to);
            return;
        }

        const OffsetTerm& offset = std::get<OffsetTerm>(term);

        // CHECKED BEFORE THE OFFSET IS LOWERED, for the reason spelled out on the gate below: an
        // expression that sampled time would reach ensureTimeGrid from inside the include resolution
        // that has not yet produced a grid, and recurse until the stack ran out. Here it is also
        // circular in meaning -- the offset is what decides where the exposure starts.
        if (offset.offset >= this->timeVarying.size())
        {
            throw std::runtime_error("veil: an include's offset is not covered by the time-varying "
                                     "tags; re-run the tagging pass.");
        }
        if (this->timeVarying[offset.offset] != 0)
        {
            throw std::runtime_error("veil: an include's offset must not depend on time.");
        }

        const OperandId field = this->lower(offset.offset);

        // The bounds are durationy clicks added to this, so anything but a datey would be adding
        // years to a number and calling the result a date.
        if (this->isVector(field) || this->block.operandAt(field).type.type != Type::Datey)
        {
            throw std::runtime_error("veil: an include's offset must be a datey.");
        }

        // A MISSING OFFSET IS NO EXPOSURE. That is what `period_included()` does in R, and it is
        // documented user-facing behaviour -- "an offset field that is NA yields the empty interval"
        // -- so veil has to give the same answer.
        //
        // It is tested for EXPLICITLY rather than left to fall out of the arithmetic. R's NA integer
        // is the most negative int, so an NA offset plus a non-negative bound stays in range and
        // gives a very negative click; tau would then fall below nu and the individual would drop out
        // anyway. That is the right answer for the wrong reason -- it holds only because of what the
        // sentinel's value happens to be, and because an offset bound is never negative. The test
        // says what is meant, and costs three instructions once per individual, none per slot.
        const OperandId missing = this->block.addTemp(TypeFull::createBool());
        this->emit(makeInstruction(Op::IsNa, missing, field));
        const OperandId present = this->block.addTemp(TypeFull::createBool());
        this->emit(makeInstruction(Op::Not, present, missing));
        this->conjoinIncluded(included, present);

        const OperandId zero = this->block.addConstant(datey, ScalarValue(0));
        const OperandId safeField = this->block.addTemp(datey);
        this->emit(makeInstruction(Op::Select, safeField, missing, zero, field));

        const OperandId from = this->block.addConstant(datey, ScalarValue(offset.fromClicks));
        const OperandId to = this->block.addConstant(datey, ScalarValue(offset.toClicks));
        nu = this->emitScalar(Op::Max, datey, nu, this->emitScalar(Op::Add, datey, safeField, from));
        tau = this->emitScalar(Op::Min, datey, tau, this->emitScalar(Op::Add, datey, safeField, to));
    }

    // The include, resolved into the clipped exposure the grid is then built from.
    //
    // This mirrors `period_included()` in R, which is the package's own definition of what an include
    // means, and then applies the clip the specification sets out: nu' = max(nu, from),
    // tau' = min(tau, to), with the individual counting as a death only if the clip did not cut the
    // end off. R starts its bounds at the whole representable calendar; here they start at the
    // exposure itself, which is the same answer once the two are intersected and saves two constants.
    void lowerInclude(const Include& include, OperandId& nu, OperandId& tau, OperandId& died)
    {
        const TypeFull boolean = TypeFull::createBool();
        std::optional<OperandId> included;

        const OperandId exposureEnd = tau;

        for (const IncludeTerm& term : include.terms)
        {
            if (const auto* gate = std::get_if<GateTerm>(&term))
            {
                // CHECKED BEFORE THE GATE IS LOWERED, AND THE ORDER IS THE WHOLE POINT. Lowering a
                // gate that samples time reaches lowerTime, which calls ensureTimeGrid -- from
                // inside the include resolution that has not yet produced a grid. That recurses
                // until the stack runs out, and control never returns to the check below, so this
                // one has to come first.
                //
                // It is not merely a restriction on what a user may write. An include RESOLVES to a
                // single clopen interval, and it is what decides where the exposure starts and ends,
                // so a gate that needed the sample points to be evaluated would be circular.
                if (gate->ast >= this->timeVarying.size())
                {
                    throw std::runtime_error("veil: an include's gate is not covered by the "
                                             "time-varying tags; re-run the tagging pass.");
                }
                if (this->timeVarying[gate->ast] != 0)
                {
                    throw std::runtime_error("veil: an include's gate must not depend on time.");
                }

                const OperandId condition = this->lower(gate->ast);
                if (this->isVector(condition)
                    || storageOf(this->block.operandAt(condition).type) != StorageClass::Logical)
                {
                    throw std::runtime_error("veil: an include's gate must be a per-individual logical.");
                }
                this->conjoinIncluded(included, condition);
                continue;
            }
            this->lowerIntervalTerm(term, nu, tau, included);
        }

        // An empty clipped interval contributes nothing, exactly as a false gate does, so the two
        // reach the grid by the same route rather than the grid having to test the bounds itself.
        this->conjoinIncluded(included, this->emitScalar(Op::Lt, boolean, nu, tau));

        // Death is at the end of the exposure, so it survives the clip only if the end did.
        died = this->emitScalar(
            Op::And, boolean, died,
            this->emitScalar(Op::Eq, boolean, tau, exposureEnd));

        this->includedFlag = included;
    }

    // Binds the exposure the first time anything samples time. Idempotent: the three columns dedupe
    // through columnOperand, and the grid is only described once.
    void ensureTimeGrid()
    {
        if (this->block.timeGrid().has_value()) { return; }

        // REACHED FROM INSIDE THE INCLUDE ITSELF, which cannot be answered: the include is what
        // decides where the exposure starts and ends, so anything within it that wants the sample
        // points is asking for a grid that its own result determines. The time-varying tag check on
        // a gate catches most of this, but NOT an `integrate` or `died_value`, which collapse to a
        // per-individual scalar and so are not tagged time-varying while still needing the grid.
        if (this->resolvingInclude)
        {
            throw std::runtime_error("veil: an include cannot sample time -- it is what decides the "
                                     "exposure that time is sampled over.");
        }

        if (!this->exposure.has_value())
        {
            throw std::runtime_error("veil: this expression samples time, but no exposure was supplied.");
        }

        const ExposureColumns& columns = *this->exposure;
        OperandId start = this->columnOperand(columns.start, TypeFull::createDatey());
        OperandId end = this->columnOperand(columns.end, TypeFull::createDatey());
        OperandId died = this->columnOperand(columns.died, TypeFull::createBool());

        if (this->includeObj.has_value())
        {
            const Obj& obj = this->objs.at(*this->includeObj);
            const Include* include = std::get_if<Include>(&obj);
            if (include == nullptr)
            {
                throw std::runtime_error("veil: the object supplied as an include is not one.");
            }

            // THE PROLOGUE MUST BE BUILT BEFORE ANY BODY INSTRUCTION EXISTS. The memo is shared, so
            // a gate reaching a node the body had already lowered would leave the prologue reading a
            // register the body has not written yet -- and the prologue runs first. Lowering resolves
            // an include eagerly to make that impossible; this restates it as a check rather than
            // trusting the caller to keep doing so.
            if (!this->block.body().empty())
            {
                throw std::runtime_error("veil: an include must be resolved before the body is lowered.");
            }

            // Everything the clip emits belongs to the prologue, because the grid it produces is
            // what the body's very first vector instruction depends on.
            this->intoPrologue = true;
            this->resolvingInclude = true;
            this->lowerInclude(*include, start, end, died);
            this->resolvingInclude = false;
            this->intoPrologue = false;
        }

        this->block.useTimeGrid(
            start, end, died, columns.deltaTClicks,
            this->includedFlag.value_or(invalidOperandId));
    }

    // Reads a click-backed or logical operand as a number of years, or as 0 and 1. Anything already
    // a number -- which includes every vector -- is handed straight back.
    OperandId asNumber(OperandId operand)
    {
        if (this->isVector(operand)) { return operand; }

        const TypeFull sourceType = this->block.operandAt(operand).type;
        const StorageClass storage = storageOf(sourceType);
        if (storage == StorageClass::Number) { return operand; }
        if (storage == StorageClass::Code)
        {
            throw std::runtime_error("veil: a category or text code is not a number.");
        }

        const OperandId result = this->block.addTemp(TypeFull::createDouble());
        this->emit(makeInstruction(Op::ToDouble, result, operand));
        return result;
    }

    // The guard that a scalar instruction passes before it is emitted. It restates, as a check, what
    // the passes above are supposed to have arranged -- so that if one of them stops firing, the
    // failure is a message naming the op rather than an interpreter reading a register as the wrong
    // sort of thing.
    void requireScalarOperands(Op op, const TypeFull& resultType, const std::vector<OperandId>& args) const
    {
        const StorageClass resultStorage = storageOf(resultType);
        const OpInfo& info = opInfo(op);

        // Checked before anything indexes `args`, so a malformed call is a message rather than a
        // read past the end of the operand list.
        if (args.size() != info.arity)
        {
            throw std::runtime_error("veil: `" + opName(op) + "` was given the wrong number of arguments.");
        }

        auto storageAt = [&](size_t position) {
            return storageOf(this->block.operandAt(args[position]).type);
        };
        auto require = [&](bool condition) {
            if (!condition)
            {
                throw std::runtime_error("veil: `" + opName(op) + "` was given operands it cannot take.");
            }
        };

        switch (info.category)
        {
            case OpCategory::Arithmetic:
                // Coercion insertion has widened everything that needed widening, so the operands
                // and the result agree; a number never meets a raw click here.
                require(resultStorage == StorageClass::Number || resultStorage == StorageClass::Click);
                for (size_t i = 0; i < args.size(); ++i) { require(storageAt(i) == resultStorage); }
                return;

            case OpCategory::Logical:
                require(resultStorage == StorageClass::Logical);
                for (size_t i = 0; i < args.size(); ++i) { require(storageAt(i) == StorageClass::Logical); }
                return;

            case OpCategory::Comparison:
                require(resultStorage == StorageClass::Logical);
                if (op == Op::IsNa)
                {
                    require(storageAt(0) == StorageClass::Number || storageAt(0) == StorageClass::Click);
                    return;
                }
                // Both sides are read the same way, which asNumber has already arranged where they
                // started out differently.
                require(storageAt(0) == storageAt(1));
                if (storageAt(0) == StorageClass::Code)
                {
                    // TWO SEPARATE RULES, and only the second of them ever lifted. Ordering a code is
                    // refused permanently -- a code says which value, not where it sits in a
                    // sequence -- whereas equality was refused only until a literal could be encoded
                    // to a code, which passEncodeText now does. Keeping them apart is what let the
                    // second lift without touching the first.
                    if (isOrderingComparison(op))
                    {
                        throw std::runtime_error(
                            "veil: `" + opName(op) + "` cannot order a text or category value.");
                    }
                }
                return;

            case OpCategory::Selection:
                require(storageAt(0) == StorageClass::Logical);
                require(storageAt(1) == resultStorage && storageAt(2) == resultStorage);
                return;

            case OpCategory::Conversion:
                if (op != Op::ToDouble)
                {
                    throw std::runtime_error("veil: `" + opName(op) + "` cannot be lowered as a scalar.");
                }
                require(resultStorage == StorageClass::Number);
                require(storageAt(0) != StorageClass::Code);
                return;

            case OpCategory::TimeVector:
            case OpCategory::Finalise:
                break;
        }

        throw std::runtime_error("veil: `" + opName(op) + "` cannot be lowered as a scalar.");
    }

    void emitByArity(Op op, OperandId result, const std::vector<OperandId>& args)
    {
        switch (args.size())
        {
            case 1: this->emit(makeInstruction(op, result, args[0])); return;
            case 2: this->emit(makeInstruction(op, result, args[0], args[1])); return;
            case 3: this->emit(makeInstruction(op, result, args[0], args[1], args[2])); return;
            default:
                throw std::runtime_error("veil: `" + opName(op) + "` has an arity veil cannot lower.");
        }
    }

    // `.t` itself: the sample points, read as years.
    OperandId lowerTime()
    {
        this->ensureTimeGrid();
        const OperandId result = this->block.addVectorTemp();
        this->emit(Instruction{Op::VectorT, result, {invalidOperandId, invalidOperandId, invalidOperandId}, 0});
        return result;
    }

    // `.t - offset`, which is `.x` whenever the offset is birth. Worth its own op rather than a
    // vector_t followed by a subtraction, because the subtraction then happens in exact clicks and
    // only the difference is converted -- one rounding instead of two, on the value ages come from.
    OperandId lowerDuration(NodeId offsetNode)
    {
        this->ensureTimeGrid();
        const OperandId offset = this->lower(offsetNode);
        if (this->isVector(offset) || storageOf(this->block.operandAt(offset).type) != StorageClass::Click)
        {
            throw std::runtime_error("veil: subtracting from `.t` needs a per-individual date.");
        }

        const OperandId result = this->block.addVectorTemp();
        this->emit(makeInstruction(Op::VectorDurn, result, offset));
        return result;
    }

    // Collapses a vector back to a scalar. A scalar argument is broadcast first, so `integrate(1)`
    // means what it should -- the exposure length -- rather than being refused.
    OperandId lowerFinalise(Op op, const std::vector<NodeId>& args)
    {
        if (args.size() != 1)
        {
            throw std::runtime_error("veil: `" + opName(op) + "` takes one argument.");
        }

        this->ensureTimeGrid();
        OperandId value = this->asNumber(this->lower(args[0]));
        if (!this->isVector(value))
        {
            const OperandId broadcast = this->block.addVectorTemp();
            this->emit(makeInstruction(Op::Broadcast, broadcast, value));
            value = broadcast;
        }

        const OperandId result = this->block.addTemp(TypeFull::createDouble());
        this->emit(makeInstruction(op, result, value));
        return result;
    }

    // A call whose result varies over time. Every operand is read as a number, because a vector holds
    // doubles and nothing else; a scalar operand stays scalar and the interpreter broadcasts it,
    // which is why no explicit Broadcast is needed here.
    OperandId lowerVectorCall(Op op, const std::vector<NodeId>& argNodes)
    {
        const OpInfo& info = opInfo(op);
        if (!info.vectorisable)
        {
            throw std::runtime_error("veil: `" + opName(op) + "` is not defined over the time vector.");
        }
        if (argNodes.size() != info.arity)
        {
            throw std::runtime_error("veil: `" + opName(op) + "` was given the wrong number of arguments.");
        }

        std::vector<OperandId> args;
        args.reserve(argNodes.size());
        for (size_t position = 0; position < argNodes.size(); ++position)
        {
            const OperandId operand = this->lower(argNodes[position]);

            // A select's condition is a per-individual logical and stays one; the tree phase has
            // already checked that it does not vary over time, which is what makes the choice once
            // per individual rather than once per slot.
            if (op == Op::Select && position == 0)
            {
                if (this->isVector(operand)
                    || storageOf(this->block.operandAt(operand).type) != StorageClass::Logical)
                {
                    throw std::runtime_error("veil: `select` needs a per-individual logical condition.");
                }
                args.push_back(operand);
                continue;
            }

            args.push_back(this->asNumber(operand));
        }

        const OperandId result = this->block.addVectorTemp();
        this->emitByArity(op, result, args);
        return result;
    }

    // log mu read from an age-period table over the time vector: the one time-vector op that arrives
    // in a tree rather than being built here.
    //
    // THE TABLE IS NOT AN OPERAND, and its node must not go anywhere near `lower`. It carries no
    // value type -- the annotation pass leaves a table untyped on purpose, precisely so that a table
    // reaching any other position is caught -- and it is a compile-time parameter of the instruction
    // rather than something that flows through a register. So it is read straight out of the ObjStore
    // and copied into the block, and only `birth` is lowered.
    OperandId lowerVectorLogMu(const std::vector<NodeId>& argNodes)
    {
        if (argNodes.size() != 2)
        {
            throw std::runtime_error("veil: `vector_log_mu` takes a table and a birth date.");
        }

        // Restated here rather than assumed from the annotation pass, per the rule that a pass
        // guards its own preconditions: if that check ever stops firing, this is a message naming
        // the op instead of a variant access reading the wrong alternative.
        const Node& tableNode = this->tree.at(argNodes[0]);
        if (!isObj(tableNode))
        {
            throw std::runtime_error("veil: `vector_log_mu` expects a mortality table as its first argument.");
        }
        const Obj& obj = this->objs.at(std::get<ObjPayload>(tableNode.payload).obj);
        const MortalityTable* table = std::get_if<MortalityTable>(&obj);
        if (table == nullptr)
        {
            throw std::runtime_error("veil: `vector_log_mu` expects a mortality table as its first argument.");
        }
        if (!table->isWellFormed())
        {
            throw std::runtime_error("veil: a mortality table's log mu does not match its dimensions.");
        }

        this->ensureTimeGrid();

        // Birth is per-individual and stays in clicks. It must not become a vector or a number: the
        // lookup indexes a cohort lattice by exact integer arithmetic, and a birth read as years
        // would put the interpolation weights very slightly off.
        const OperandId birth = this->lower(argNodes[1]);
        if (this->isVector(birth)
            || storageOf(this->block.operandAt(birth).type) != StorageClass::Click)
        {
            throw std::runtime_error("veil: `vector_log_mu` needs a per-individual birth date.");
        }

        const TableId tableId = this->block.addTable(*table);
        const OperandId result = this->block.addVectorTemp();
        this->emit(makeParameterisedInstruction(Op::VectorLogMu, result, tableId, birth));
        return result;
    }

    OperandId lowerCall(NodeId id, const Node& node, const CallPayload& call)
    {
        // Stated once here rather than at each `node.type.value()` below. An untyped call means a
        // pass appended a node without annotating it -- the sort of slip a `std::optional<TypeFull>`
        // assignment error already caused once -- and without this it surfaces as the standard
        // library's own "bad optional access" with nothing naming the op.
        if (!node.type.has_value())
        {
            throw std::runtime_error("veil: `" + opName(call.op) + "` carries no type; a pass added "
                                     "a node without annotating it.");
        }

        // The two ops that turn a vector back into a per-individual number.
        if (call.op == Op::Integrate || call.op == Op::DiedValue)
        {
            return this->lowerFinalise(call.op, call.args);
        }

        if (call.op == Op::VectorLogMu)
        {
            return this->lowerVectorLogMu(call.args);
        }

        // `.t - offset` is recognised before its arguments are lowered, so that the subtraction
        // stays in clicks rather than becoming a vector subtraction after the conversion.
        //
        // The result must itself be click-typed, which is what distinguishes `.t - .b` from
        // `.t - 5`. The latter is a date less a number of years, so it is already a double by the
        // datey rules, and there is no exact integer subtraction to preserve.
        if (call.op == Op::Sub && call.args.size() == 2 && isTime(this->tree.at(call.args[0]))
            && !this->timeVarying[call.args[1]]
            && storageOf(node.type.value()) == StorageClass::Click)
        {
            return this->lowerDuration(call.args[1]);
        }

        // A conversion applied to something already vectorised has nothing left to do: a vector
        // holds years already.
        if (call.op == Op::ToDouble && call.args.size() == 1 && this->timeVarying[call.args[0]])
        {
            return this->asNumber(this->lower(call.args[0]));
        }

        if (this->timeVarying[id]) { return this->lowerVectorCall(call.op, call.args); }

        std::vector<OperandId> args;
        args.reserve(call.args.size());
        for (const NodeId argId : call.args) { args.push_back(this->lower(argId)); }

        // A comparison is the one place a mixture survives the passes above: the coercion pass
        // leaves comparisons alone on purpose, so a click can still meet a plain number when
        // narrowing declined to rewrite the threshold. Read both sides as years and compare those,
        // which is the datey package's own rule.
        if (isComparisonOp(call.op) && args.size() == 2
            && storageOf(this->block.operandAt(args[0]).type)
                   != storageOf(this->block.operandAt(args[1]).type))
        {
            args[0] = this->asNumber(args[0]);
            args[1] = this->asNumber(args[1]);
        }

        const TypeFull resultType = node.type.value();
        this->requireScalarOperands(call.op, resultType, args);

        const OperandId result = this->block.addTemp(resultType);
        this->emitByArity(call.op, result, args);
        return result;
    }

    OperandId lower(NodeId id)
    {
        if (this->operandOfNode[id] != invalidOperandId) { return this->operandOfNode[id]; }

        const Node& node = this->tree.at(id);

        if (isObj(node))
        {
            throw std::runtime_error("veil: a mortality or include object cannot yet be lowered.");
        }

        OperandId operand = invalidOperandId;
        if (isTime(node))
        {
            operand = this->lowerTime();
        }
        else if (isLit(node))
        {
            const LitPayload& lit = std::get<LitPayload>(node.payload);
            operand = this->block.addConstant(lit.type, literalValue(lit));
        }
        else if (isField(node))
        {
            const FieldPayload& field = std::get<FieldPayload>(node.payload);
            operand = this->columnOperand(field.column, node.type.value());
        }
        else
        {
            operand = this->lowerCall(id, node, std::get<CallPayload>(node.payload));
        }

        this->operandOfNode[id] = operand;
        return operand;
    }
};

} // namespace detail

// Lowers a rewritten tree to a block whose single output is the value of the root. `timeVarying` is
// what passTagTimeVarying returned for this tree, and must have been taken after the last pass that
// added a node -- a short one would leave the trailing nodes unchecked, which is exactly how a time
// vector would slip into a scalar register.
//
// `objs` holds the concept objects the tree refers to by ObjId. Lowering needs it because a mortality
// table is copied into the block as a parameter of the instruction that reads it, rather than being
// left behind a reference to a store the block does not own.
//
// `exposure` says where the time vector comes from. Leaving it empty is how a caller states that the
// expression must not sample time; anything that does is then refused by name rather than reading a
// grid that was never built.
inline Block passLowerToBlock(
    const Tree& tree,
    const ObjStore& objs,
    const std::vector<char>& timeVarying,
    const std::optional<ExposureColumns>& exposure = std::nullopt,
    std::optional<ObjId> includeObj = std::nullopt)
{
    if (tree.roots().empty())
    {
        throw std::runtime_error("veil: a tree with no root cannot be lowered.");
    }
    if (timeVarying.size() != tree.size())
    {
        throw std::runtime_error("veil: the time-varying tags do not cover the tree; re-run the "
                                 "tagging pass after anything that adds a node.");
    }

    detail::Lowerer lowerer(tree, objs, timeVarying, exposure, includeObj);

    // An include is resolved EAGERLY, before the body, rather than when something first samples time.
    // Two reasons, and the second is the one that matters: the clip decides the exposure, so it has
    // to be settled before anything reads the grid; and the memo lowering shares between the two
    // would otherwise let the prologue reuse an operand the body defines, which runs later.
    if (includeObj.has_value()) { lowerer.ensureTimeGrid(); }

    // ONE OUTPUT PER ROOT, in the order the roots were given, so a caller reads its results back by
    // the position it put them in. The roots share one memo, which is what makes a sub-expression two
    // outputs have in common cost one instruction rather than two -- the mortality shared by an AEV's
    // E and V being the case that matters.
    for (const NodeId rootId : tree.roots())
    {
        const OperandId result = lowerer.lower(rootId);

        if (lowerer.block.operandAt(result).isVector())
        {
            throw std::runtime_error("veil: this expression varies over time, so it has no single "
                                     "value per individual; integrate it, or take its value at death.");
        }

        lowerer.block.addOutput(result);
    }

    return std::move(lowerer.block);
}

} // namespace veil
