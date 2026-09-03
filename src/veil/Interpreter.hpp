// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>
#include "datey.h" // yearsFromClicks, roundBankers, isValidDatey, isValidDurationy
#include "vec_ops/vec_ops.hpp" // tier::exp_V_V, tier::log_V_V -- R-free, so the core stays R-free
#include "veil/Block.hpp"
#include "veil/ColumnView.hpp"
#include "veil/Instruction.hpp"
#include "veil/MortalityTable.hpp"
#include "veil/Op.hpp"
#include "veil/Operand.hpp"
#include "veil/ScalarValue.hpp"
#include "veil/TimeGrid.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// The veil interpreter: a register file, and a switch that runs one block over it.
//
// A switch over op monikers is the right shape here rather than a compromise. Dispatch costs the
// same whatever an instruction goes on to do, while the work an instruction does grows with the time
// vector, so the overhead amortises across the slots of an individual's exposure. The scalar stretch
// this file currently covers is the part with nothing to amortise against, and it is also the part
// that is a handful of operations per individual.
//
// SCOPE: THE SCALAR SPINE, matching passLowerToBlock. Time-vector and finalise ops are refused.
//
// IEEE 754 IS FOLLOWED EVERYWHERE, AND THE EDITION IS 754-2019. logmu's arithmetic has to mean the
// same thing on every platform, in every R build, and through every front end, and the standard is
// the only thing that makes each operation mean one fixed thing. So where R and the standard part
// company, the standard wins -- and the two places that look like departures from IEEE are not:
//
//   - A comparison is never missing. veil draws no distinction between R's NA_real_ and any other
//     NaN, so `NaN == NaN` is false and `NaN != NaN` is true where R would give NA. That is the
//     unordered rule of 754, and the veil specification's own.
//   - `min` and `max` PROPAGATE a NaN. This is 754-2019 `minimum` / `maximum`, not a concession to
//     R. The 2008 edition's `minNum` / `maxNum` -- which is what C's fmin / fmax implement --
//     DISCARDED a NaN, and 2019 removed them, because discarding is not associative in the presence
//     of a NaN and so breaks reductions and any reassociation a vectoriser wants to do. The
//     discarding behaviour survives in 2019 under the separate names `minimumNumber` /
//     `maximumNumber`, which veil does not use.
//
// Rounding is banker's rounding, from datey.h, because rounding here is definitional -- it decides
// what a value means, and logmu and datey must never disagree about a boundary.
//
// NOTE FOR THE VECTOR PHASE: EVE ships `minimumNumber` (as `eve::min[eve::numeric]`) but NOT
// `minimum`. Undecorated `eve::min` lowers to the raw MINPD, which returns its second operand on any
// unordered or zero-tied pair and so is order-dependent; `eve::min[eve::pedantic]` gets the signed
// zeros right but propagates a NaN only from its first argument. The composition that matches this
// file is `if_else(is_unordered(a, b), nan, min[pedantic](a, b))`. It costs about a nanosecond an
// element, which lowering can avoid where the scan and the intervals prove no NaN and no zero
// crossing can reach the operands.

namespace detail
{

inline std::string interpreterOpName(Op op) { return std::string(moniker(op)); }

[[noreturn]] inline void opNotRunnable(Op op, const char* why)
{
    throw std::runtime_error("veil: `" + interpreterOpName(op) + "` " + why);
}

inline double quietNaN() noexcept { return std::numeric_limits<double>::quiet_NaN(); }

// IEEE 754-2019 `minimum` and `maximum` (clause 9.6). Two rules beyond the obvious comparison:
//
//   - a NaN in either operand yields a NaN, symmetrically;
//   - the zeros are ordered, so `minimum(+0, -0)` is -0 and `maximum(+0, -0)` is +0.
//
// The second is what a plain `a < b ? a : b` gets wrong, and only in one argument order: +0 and -0
// compare equal, so the comparison is false and the second operand comes back whichever way round
// they were given. It is observable through anything that reads the sign of a zero -- a division
// yielding an infinity is the usual route.
//
// `a == b` is true for a genuine tie as well as for the zeros, and returning either operand is
// correct there, so the zero branch needs no test for zero-ness.
inline double ieeeMinimum(double a, double b) noexcept
{
    if (std::isnan(a) || std::isnan(b)) { return quietNaN(); }
    if (a == b) { return std::signbit(a) ? a : b; }
    return a < b ? a : b;
}

inline double ieeeMaximum(double a, double b) noexcept
{
    if (std::isnan(a) || std::isnan(b)) { return quietNaN(); }
    if (a == b) { return std::signbit(a) ? b : a; }
    return a > b ? a : b;
}

// Back to the 32 bits R holds a datey in. THE ONLY PLACE THIS IS NEEDED IS AN OUTPUT CROSSING TO R:
// inside the engine a click is 64-bit and cannot overflow, so nothing on the way checks. A result
// outside the range is a real error and this is where someone can act on it.
inline int narrowClicksForOutput(int64_t value)
{
    if (value < static_cast<int64_t>(std::numeric_limits<int>::min())
        || value > static_cast<int64_t>(std::numeric_limits<int>::max()))
    {
        throw std::runtime_error("veil: a datey/durationy result is outside the range R can hold.");
    }
    return static_cast<int>(value);
}

inline double unaryNumber(Op op, double x)
{
    switch (op)
    {
        case Op::Pos: return x;
        case Op::Neg: return -x;
        case Op::Abs: return std::abs(x);
        case Op::Sqrt: return std::sqrt(x);
        case Op::Exp: return std::exp(x);
        case Op::Log: return std::log(x);
        case Op::Log10: return std::log10(x);
        case Op::Expm1: return std::expm1(x);
        case Op::Log1p: return std::log1p(x);
        case Op::Sin: return std::sin(x);
        case Op::Cos: return std::cos(x);
        case Op::Floor: return std::floor(x);
        case Op::Ceiling: return std::ceil(x);
        case Op::Round: return roundBankers(x);
        case Op::Trunc: return std::trunc(x);
        // R's sign: minus one, zero or one, and a NaN back for a NaN, which falls out of both
        // comparisons being false.
        case Op::Sign: return x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : x);
        default: break;
    }
    opNotRunnable(op, "is not a unary numeric operation.");
}

inline double binaryNumber(Op op, double a, double b)
{
    switch (op)
    {
        case Op::Add: return a + b;
        case Op::Sub: return a - b;
        case Op::Mul: return a * b;
        case Op::RDiv: return a / b;
        case Op::Pow: return std::pow(a, b);
        case Op::Min: return ieeeMinimum(a, b);
        case Op::Max: return ieeeMaximum(a, b);
        // R's `%/%` and `%%` are floored division, which is what actuarial banding wants: an age
        // band is the same width either side of zero.
        case Op::IDiv: return std::floor(a / b);
        case Op::IMod: return a - std::floor(a / b) * b;
        default: break;
    }
    opNotRunnable(op, "is not a binary numeric operation.");
}

inline double ternaryNumber(Op op, double a, double b, double c)
{
    switch (op)
    {
        // clamp is veil's own, so its edge cases are a choice rather than a standard. It is defined
        // as the composition, which is the only way it does not become a third set of NaN and
        // signed-zero rules to keep in step with the other two.
        case Op::Clamp: return ieeeMinimum(ieeeMaximum(a, b), c);
        case Op::Fma: return std::fma(a, b, c);
        default: break;
    }
    opNotRunnable(op, "is not a ternary numeric operation.");
}

// The click arithmetic the type rules actually admit. Multiplying, dividing or taking the modulus of
// two dates is undefined and the annotation pass has already refused it, so anything else reaching
// here is a lowering fault rather than a user error.
//
// min, max and clamp need none of the care their double counterparts do: an integer has no NaN and
// no signed zero, so the plain comparison IS 754-2019 `minimum` here. Worth knowing because the
// include clip -- max(nu, from), min(tau, to) -- is exactly this case, so the compliant form costs
// nothing where min and max are used most.
//
// NOTHING HERE CHECKS FOR OVERFLOW, and nothing needs to: the operations below are the whole of click
// arithmetic, there is no multiplication among them, and only the first four make a new value. So any
// click expression is a sum of 32-bit inputs with coefficients of plus or minus one, and int64 has
// room for billions of terms. See the note on ScalarValue.
inline int64_t clickArithmetic(Op op, const std::array<int64_t, MaxArgCount>& args)
{
    switch (op)
    {
        case Op::Pos: return args[0];
        case Op::Neg: return -args[0];
        case Op::Abs: return args[0] < 0 ? -args[0] : args[0];
        case Op::Add: return args[0] + args[1];
        case Op::Sub: return args[0] - args[1];
        case Op::Min: return args[0] < args[1] ? args[0] : args[1];
        case Op::Max: return args[0] > args[1] ? args[0] : args[1];
        case Op::Clamp: return args[0] < args[1] ? args[1] : (args[0] > args[2] ? args[2] : args[0]);
        default: break;
    }
    opNotRunnable(op, "is not defined on datey or durationy values.");
}

inline bool logical(Op op, const std::array<bool, MaxArgCount>& args)
{
    switch (op)
    {
        case Op::Not: return !args[0];
        case Op::And: return args[0] && args[1];
        case Op::Or: return args[0] || args[1];
        case Op::Xor: return args[0] != args[1];
        default: break;
    }
    opNotRunnable(op, "is not a logical operation.");
}

template <typename T>
inline bool compareValues(Op op, const T& a, const T& b)
{
    switch (op)
    {
        case Op::Eq: return a == b;
        case Op::Ne: return a != b;
        case Op::Lt: return a < b;
        case Op::Le: return a <= b;
        case Op::Gt: return a > b;
        case Op::Ge: return a >= b;
        default: break;
    }
    opNotRunnable(op, "is not a comparison.");
}

// One value read out of a column. The view checks its own tag before handing back a pointer, so the
// only thing to get right here is which accessor a type asks for.
inline ScalarValue valueFromColumn(const ColumnView& column, size_t record)
{
    switch (column.type.type)
    {
        case Type::Double: return ScalarValue(column.doubles()[record]);
        case Type::Datey:
        case Type::Durationy:
        case Type::Category: return ScalarValue(static_cast<int64_t>(column.ints()[record]));
        case Type::Bool: return ScalarValue(column.bools()[record]);
        case Type::Text:
        case Type::DateyInterval: break;
    }
    throw std::runtime_error("veil: a column of this type cannot be read into an operand.");
}

} // namespace detail

class Interpreter final
{
public:
    // Sizes the register file from the block's operand table and writes the constants, which do not
    // change from one individual to the next. LIFETIME: the block must outlive the interpreter.
    explicit Interpreter(const Block& block)
        : block(block)
    {
        this->registers.reserve(block.operandCount());

        // One buffer per PHYSICAL slot, not per operand: where allocation has run, operands whose
        // live ranges do not overlap share a buffer, and this is the number that says how many the
        // block actually needs at once.
        this->buffers.resize(block.vectorBufferCount());
        for (OperandId id = 0; id < static_cast<OperandId>(block.operandCount()); ++id)
        {
            // A vector operand's scalar register is never read; it is filled anyway so that the two
            // tables stay the same length and an OperandId indexes both.
            this->registers.push_back(zeroValueOf(block.operandAt(id).type));
        }

        for (const ConstantBinding& constant : block.constants())
        {
            if (constant.value.index() != this->registers.at(constant.operand).index())
            {
                throw std::runtime_error("veil: a constant does not match the type of its operand.");
            }
            this->registers[constant.operand] = constant.value;
        }
    }

    // Writes each column-bound operand from one individual's row, and resolves the include.
    // `columns` is indexed by the same ColumnId the block's field references carry, and may hold a
    // null where the host did not read a column; a block that binds one of those is an error rather
    // than a zero.
    //
    // ANSWERS WHETHER THIS INDIVIDUAL CONTRIBUTES AT ALL. False means the include emptied them, and
    // the caller should not run the body: their contribution is the identity, which is zero for a
    // sum, and adding a zero is the same arithmetic as adding nothing. Skipping is therefore free of
    // any effect on the answer and saves the whole body, which for an excluded individual is the
    // entire cost of them.
    //
    // The grid is cleared before answering false, so `slotCount` describes this individual rather
    // than whoever came before them.
    bool loadRecord(const std::vector<const ColumnView*>& columns, size_t record)
    {
        for (const ColumnBinding& binding : this->block.columns())
        {
            const ColumnView* column = binding.column < columns.size() ? columns[binding.column] : nullptr;
            if (column == nullptr)
            {
                throw std::runtime_error("veil: the block reads a column the host did not supply.");
            }
            if (record >= column->count)
            {
                throw std::runtime_error("veil: the block was run past the end of a column.");
            }
            if (!(column->type == this->block.operandAt(binding.operand).type))
            {
                throw std::runtime_error("veil: a column's type does not match the operand it feeds.");
            }
            this->registers[binding.operand] = detail::valueFromColumn(*column, record);
        }

        // The prologue resolves the include, and so decides what the exposure actually is. It runs
        // between loading the columns and building the grid because it reads the first and produces
        // what the second is built from.
        for (const Instruction& instruction : this->block.prologue()) { this->execute(instruction); }

        // The time vector is rebuilt from the exposure this individual just loaded. The slot buffer
        // is a member and only ever resized, so after the first few individuals it stops allocating.
        if (const std::optional<TimeGridBinding>& binding = this->block.timeGrid())
        {
            // An excluded individual is reported rather than run. An empty grid would answer
            // correctly too -- every vector would hold nothing, the integral would sum nothing and
            // there would be no death slot -- but it walks the body to reach that conclusion.
            // buildTimeGrid keeps refusing an inverted exposure as the fault it is, because an
            // exclusion never reaches it.
            if (binding->included != invalidOperandId && !this->logicalAt(binding->included))
            {
                this->grid = TimeGrid{};
                return false;
            }

            this->grid = buildTimeGrid(this->clicksAt(binding->start), this->clicksAt(binding->end),
                                       this->logicalAt(binding->died), binding->deltaTClicks);
        }
        return true;
    }

    // How many slots this individual's time vector holds, the death slot included. Meaningful only
    // after loadRecord, and only for a block that samples time.
    int slotCount() const noexcept { return this->grid.slotCount(); }

    // Every slot of every time vector this interpreter has filled since it was built, added up over
    // instructions and over individuals.
    //
    // THIS IS THE ONLY THING THAT SEES SLOT DEMAND WORKING. Narrowing where a value is computed
    // cannot change an answer -- that is the whole claim -- so no result tells you whether it
    // happened. This does, and it is also the measurement the design asks the proof of concept to
    // take before anyone argues about interpreter dispatch: it is the count that dispatch amortises
    // against. One addition per instruction, not per slot.
    size_t slotEvaluations() const noexcept { return this->slotsFilled; }

    // Runs the body once, in order. There are no jumps: the block is straight-line code.
    void run()
    {
        for (const Instruction& instruction : this->block.body()) { this->execute(instruction); }
    }

    double outputAsNumber(size_t index) const { return this->numberAt(this->outputOperand(index)); }
    int outputAsClicks(size_t index) const
    {
        return detail::narrowClicksForOutput(this->clicksAt(this->outputOperand(index)));
    }
    bool outputAsLogical(size_t index) const { return this->logicalAt(this->outputOperand(index)); }

private:
    OperandId outputOperand(size_t index) const
    {
        const std::vector<OperandId>& outputs = this->block.outputs();
        if (index >= outputs.size())
        {
            throw std::runtime_error("veil: the block has no output at that position.");
        }
        return outputs[index];
    }

    double numberAt(OperandId id) const
    {
        if (const double* value = std::get_if<double>(&this->registers.at(id))) { return *value; }
        throw std::runtime_error("veil: an operand was read as a number but does not hold one.");
    }

    int64_t clicksAt(OperandId id) const
    {
        if (const int64_t* value = std::get_if<int64_t>(&this->registers.at(id))) { return *value; }
        throw std::runtime_error("veil: an operand was read as clicks but does not hold them.");
    }

    bool logicalAt(OperandId id) const
    {
        if (const bool* value = std::get_if<bool>(&this->registers.at(id))) { return *value; }
        throw std::runtime_error("veil: an operand was read as a logical but does not hold one.");
    }

    // Storing checks the alternative rather than replacing it, so an instruction writing the wrong
    // sort of value is caught where it happens instead of at whatever reads the register next.
    template <typename T>
    void store(OperandId id, T value)
    {
        ScalarValue& target = this->registers.at(id);
        if (!std::holds_alternative<T>(target))
        {
            throw std::runtime_error("veil: an instruction wrote a value of the wrong sort into an operand.");
        }
        target = value;
    }

    void execute(const Instruction& instruction)
    {
        if (this->block.operandAt(instruction.result).isVector())
        {
            this->executeVector(instruction);
            return;
        }
        if (instruction.op == Op::Integrate || instruction.op == Op::DiedValue)
        {
            this->executeFinalise(instruction);
            return;
        }

        switch (opInfo(instruction.op).category)
        {
            case OpCategory::Arithmetic: this->executeArithmetic(instruction); return;
            case OpCategory::Logical: this->executeLogical(instruction); return;
            case OpCategory::Comparison: this->executeComparison(instruction); return;
            case OpCategory::Selection: this->executeSelection(instruction); return;
            case OpCategory::Conversion: this->executeConversion(instruction); return;
            case OpCategory::TimeVector:
            case OpCategory::Finalise: break;
        }
        detail::opNotRunnable(instruction.op, "needs the time vector, which this interpreter does not run.");
    }

    // One argument of a vector instruction, resolved once rather than per slot: either a run of
    // values, or a single one to be read at every slot. That is what makes an explicit Broadcast
    // unnecessary for a mixed scalar-and-vector operation, exactly as the specification says.
    struct Lane final
    {
        const double* values = nullptr; // Null when the argument is a per-individual scalar.
        double scalar = 0.0;

        double at(int slot) const noexcept { return this->values != nullptr ? this->values[slot] : this->scalar; }
    };

    Lane laneOf(OperandId operand) const
    {
        if (this->block.operandAt(operand).isVector())
        {
            return Lane{this->bufferOf(operand).data(), 0.0};
        }
        return Lane{nullptr, this->numberAt(operand)};
    }

    // The physical buffer an operand's time vector lives in. The indirection is the whole of what
    // buffer allocation costs at run time, and it is one bounds-checked array read per vector
    // instruction rather than per slot.
    const std::vector<double>& bufferOf(OperandId operand) const
    {
        return this->buffers.at(this->block.vectorBufferOf(operand));
    }

    std::vector<double>& bufferOf(OperandId operand)
    {
        return this->buffers.at(this->block.vectorBufferOf(operand));
    }

    void executeVector(const Instruction& instruction)
    {
        // WHICH SLOTS THIS VALUE IS ACTUALLY READ AT. A value feeding only `died_value` is wanted at
        // the death slot and nowhere else, and for an individual who did not die there is no death
        // slot, so the range is empty and the instruction does nothing at all. That is the delta
        // gate the design asks for -- it costs a comparison, not a branch in the block, and it is
        // why survivors never walk A's integrand.
        //
        // Slots outside the range keep whatever the buffer last held. Nothing reads them: that is
        // exactly what the demand analysis established.
        const SlotDemand demand = this->block.slotDemandOf(instruction.result);
        const int first = demand.midpoints ? 0 : this->grid.integrationSlots;
        const int slots = demand.death ? this->grid.slotCount() : this->grid.integrationSlots;

        // Sized before the arguments are resolved: a Lane holds a pointer into a buffer, and a
        // resize could move one. The allocator never gives a result the buffer of one of its own
        // arguments, so the two can never be the same store, but the ordering costs nothing and
        // does not rely on that being true.
        //
        // Sized to the WHOLE grid rather than to the range written, so that a slot index means the
        // same thing in every buffer and the death slot sits where its reader looks for it.
        std::vector<double>& out = this->bufferOf(instruction.result);
        out.resize(static_cast<size_t>(this->grid.slotCount()));

        this->slotsFilled += slots > first ? static_cast<size_t>(slots - first) : 0;

        switch (instruction.op)
        {
            case Op::VectorT:
                // The sample points themselves, as years.
                for (int slot = first; slot < slots; ++slot)
                {
                    out[static_cast<size_t>(slot)] = yearsFromClicks(this->grid.clicksAt(slot));
                }
                return;

            case Op::VectorDurn:
            {
                // A duration from a fixed date -- an age, when the date is birth. The subtraction is
                // done in clicks and only the difference is converted, so this is one rounding
                // rather than two and the ages of two individuals born the same day agree exactly.
                const int64_t offset = this->clicksAt(instruction.args[0]);
                for (int slot = first; slot < slots; ++slot)
                {
                    out[static_cast<size_t>(slot)] =
                        yearsFromClicks(this->grid.clicksAt(slot) - offset);
                }
                return;
            }

            case Op::VectorLogMu:
            {
                // The table is a parameter of the instruction, not an operand: it is fixed when the
                // block is built and shared by every individual, so it never enters a register.
                //
                // Checked by name because `parameter` is an untyped uint32_t whose meaning is the
                // op's alone: an instruction lowered without one would otherwise index the table
                // list with -1 and surface as the standard library's own out-of-range message.
                if (instruction.parameter == invalidParameter)
                {
                    detail::opNotRunnable(instruction.op, "was lowered without its table parameter.");
                }
                const MortalityTable& table = this->block.tableAt(instruction.parameter);
                const int64_t birth = this->clicksAt(instruction.args[0]);

                // Read at the sample point in CLICKS. The lookup does its own integer arithmetic to
                // find the surrounding lattice points and the weights between them, so handing it
                // years would round twice and shift the interpolation.
                for (int slot = first; slot < slots; ++slot)
                {
                    out[static_cast<size_t>(slot)] =
                        logMuAt(table, birth, this->grid.clicksAt(slot));
                }
                return;
            }

            case Op::Broadcast:
            {
                const double value = this->numberAt(instruction.args[0]);
                for (int slot = first; slot < slots; ++slot) { out[static_cast<size_t>(slot)] = value; }
                return;
            }

            case Op::Select:
            {
                // The condition does not vary over time, so the choice is made once and then the
                // chosen branch is copied. Both branches were computed; skipping one is the vector
                // phase's structured if/else, which is a lowering matter rather than an interpreter
                // one.
                const Lane chosen = this->laneOf(this->logicalAt(instruction.args[0])
                                                     ? instruction.args[1] : instruction.args[2]);
                for (int slot = first; slot < slots; ++slot) { out[static_cast<size_t>(slot)] = chosen.at(slot); }
                return;
            }

            default:
                break;
        }

        const Lane a = instruction.argCount > 0 ? this->laneOf(instruction.args[0]) : Lane{};
        const Lane b = instruction.argCount > 1 ? this->laneOf(instruction.args[1]) : Lane{};
        const Lane c = instruction.argCount > 2 ? this->laneOf(instruction.args[2]) : Lane{};

        if (this->runKernel(instruction, a, first, slots, out)) { return; }

        const OpCategory category = opInfo(instruction.op).category;
        for (int slot = first; slot < slots; ++slot)
        {
            double value = 0.0;
            if (category == OpCategory::Comparison)
            {
                // A vector holds doubles and nothing else, so a comparison over time answers 1 or 0
                // rather than a logical. That is what an indicator weight becomes.
                const bool answer = instruction.op == Op::IsNa
                    ? std::isnan(a.at(slot))
                    : detail::compareValues(instruction.op, a.at(slot), b.at(slot));
                value = answer ? 1.0 : 0.0;
            }
            else
            {
                switch (instruction.argCount)
                {
                    case 1: value = detail::unaryNumber(instruction.op, a.at(slot)); break;
                    case 2: value = detail::binaryNumber(instruction.op, a.at(slot), b.at(slot)); break;
                    case 3: value = detail::ternaryNumber(instruction.op, a.at(slot), b.at(slot), c.at(slot)); break;
                    default: detail::opNotRunnable(instruction.op, "has no vector form.");
                }
            }
            out[static_cast<size_t>(slot)] = value;
        }
    }

    // BELOW THIS MANY SLOTS THE SCALAR LOOP WINS. A kernel call goes through the tier's function
    // pointer, and under a handful of elements that indirection costs more than the lanes save.
    // One threshold for every op rather than one per op: the crossing point is a property of the
    // call, not of the arithmetic behind it, and a per-op table would be tuning noise pretending
    // to be a decision.
    static constexpr int MinimumKernelSlots = 4;

    // THE SIMD KERNELS, for the ops where lanes pay for themselves. Answers true when it has
    // written the result, false when the caller should fall through to the scalar loop.
    //
    // COVERAGE HERE IS A PERFORMANCE QUESTION, NEVER A CORRECTNESS ONE. The scalar loop handles
    // every op and remains the definition; an op missing from the switch below is slow, not wrong.
    // So ops are added one at a time as a profile asks for them, never generated wholesale.
    //
    // Only `exp` and `log` are routed today, and the profile is why. A five-year monthly A/E over
    // an age-period table spends about 72% of its time on the exponential and the pointwise
    // arithmetic around it, against 16% on the table lookup and 12% on the scalar spine. One op
    // is therefore most of the calculation.
    //
    // THE ULP PROMISE IS WEAKER FOR THESE TWO THAN FOR ARITHMETIC, and that is not a concession
    // made here. Neither `exp` nor `log` is correctly rounded by any mainstream library -- IEEE
    // 754 recommends it for the transcendentals rather than requiring it -- so the kernel and the
    // scalar fallback agree to within an ulp rather than exactly. Measured against this platform's
    // `std::exp`, about one value in ten differs in the last place. Two further differences are
    // real and confined to the extremes of the range: the kernel returns an infinity from about
    // 709.44 upwards where `std::exp` still has finite room to 709.78, and it returns zero below
    // about -708.4 where `std::exp` continues into the subnormals. A log-mortality reaching either
    // band describes a hazard of 10^308 or 10^-308 a year, so neither is attainable from data.
    //
    // A SCALAR ARGUMENT IS LEFT TO THE LOOP. `Op::Broadcast` is what puts a per-individual value
    // into a buffer, so a scalar lane arriving here means the whole result is one repeated value,
    // and a kernel over it would evaluate the same exponential once per slot.
    //
    // ALIASING IS NOT A CONCERN. The allocator never gives an instruction's result the buffer of
    // one of its own arguments, so the two ranges are distinct; and were they ever the same, these
    // are element-wise unary kernels, for which writing over the input in place is well defined.
    bool runKernel(const Instruction& instruction, const Lane& a, int first, int slots, std::vector<double>& out) const
    {
        const int kernelSlots = slots - first;
        if (a.values == nullptr || kernelSlots < MinimumKernelSlots) { return false; }

        const tier::vec_size count = static_cast<tier::vec_size>(kernelSlots);
        const double* const from = a.values + first;
        double* const to = out.data() + first;

        switch (instruction.op)
        {
            case Op::Exp: tier::exp_V_V(count, from, to); return true;
            case Op::Log: tier::log_V_V(count, from, to); return true;
            default: return false;
        }
    }

    // Collapsing the time vector back to one number per individual.
    void executeFinalise(const Instruction& instruction)
    {
        if (!this->block.operandAt(instruction.args[0]).isVector())
        {
            detail::opNotRunnable(instruction.op, "needs a value that varies over time.");
        }
        const std::vector<double>& values = this->bufferOf(instruction.args[0]);

        if (instruction.op == Op::Integrate)
        {
            // The interval comes from the grid rather than from the binding: the grid is what was
            // actually built for this individual, so there is one source of truth for the width.
            this->store(instruction.result, integrateSlots(values.data(), this->grid));
            return;
        }

        // The value at the moment of death, which is nothing at all for a survivor -- and survivors
        // are the great majority, which is why lowering will eventually gate this whole branch.
        this->store(
            instruction.result,
            this->grid.hasDeathSlot ? values[static_cast<size_t>(this->grid.integrationSlots)] : 0.0);
    }

    void executeArithmetic(const Instruction& instruction)
    {
        // The result's type says which algebra applies: a plain number, or the exact click integers.
        // Lowering has already checked that the operands agree with it.
        if (this->block.operandAt(instruction.result).type.type == Type::Double)
        {
            switch (instruction.argCount)
            {
                case 1:
                    this->store(
                        instruction.result,
                        detail::unaryNumber(instruction.op, this->numberAt(instruction.args[0])));
                    return;
                case 2:
                    this->store(
                        instruction.result,
                        detail::binaryNumber(
                            instruction.op,
                            this->numberAt(instruction.args[0]),
                            this->numberAt(instruction.args[1])));
                    return;
                default:
                    this->store(
                        instruction.result,
                        detail::ternaryNumber(
                            instruction.op,
                            this->numberAt(instruction.args[0]),
                            this->numberAt(instruction.args[1]),
                            this->numberAt(instruction.args[2])));
                    return;
            }
        }

        std::array<int64_t, MaxArgCount> args = {0, 0, 0};
        for (uint8_t i = 0; i < instruction.argCount; ++i) { args[i] = this->clicksAt(instruction.args[i]); }
        this->store(instruction.result, detail::clickArithmetic(instruction.op, args));
    }

    void executeLogical(const Instruction& instruction)
    {
        std::array<bool, MaxArgCount> args = {false, false, false};
        for (uint8_t i = 0; i < instruction.argCount; ++i) { args[i] = this->logicalAt(instruction.args[i]); }
        this->store(instruction.result, detail::logical(instruction.op, args));
    }

    void executeComparison(const Instruction& instruction)
    {
        const TypeFull argType = this->block.operandAt(instruction.args[0]).type;

        if (instruction.op == Op::IsNa)
        {
            this->store(instruction.result, this->isMissing(argType, instruction.args[0]));
            return;
        }

        switch (argType.type)
        {
            case Type::Double:
                this->store(
                    instruction.result,
                    detail::compareValues(
                        instruction.op,
                        this->numberAt(instruction.args[0]),
                        this->numberAt(instruction.args[1])));
                return;
            case Type::Datey:
            case Type::Durationy:
                this->store(
                    instruction.result,
                    detail::compareValues(
                        instruction.op,
                        this->clicksAt(instruction.args[0]),
                        this->clicksAt(instruction.args[1])));
                return;
            case Type::Bool:
                this->store(
                    instruction.result,
                    detail::compareValues(
                        instruction.op,
                        this->logicalAt(instruction.args[0]),
                        this->logicalAt(instruction.args[1])));
                return;
            case Type::Category:
            {
                // A plain integer equality: every code names a string, and a category column carrying
                // an NA was refused when it was read, so there is no missing value here to give one
                // code special behaviour.
                //
                // The op is restated rather than assumed. Ordering a category is refused twice above
                // this, on type and again at lowering, and this is the third place that would give a
                // plausible wrong answer instead of a complaint if both ever stopped firing.
                if (instruction.op != Op::Eq && instruction.op != Op::Ne)
                {
                    detail::opNotRunnable(instruction.op, "cannot order a category.");
                }
                this->store(
                    instruction.result,
                    detail::compareValues(
                        instruction.op,
                        this->clicksAt(instruction.args[0]),
                        this->clicksAt(instruction.args[1])));
                return;
            }
            default:
                break;
        }
        detail::opNotRunnable(instruction.op, "cannot compare values of this type yet.");
    }

    // What counts as missing is per base type: a double is any NaN, and a click-backed value is one
    // outside the calendar the datey framework can represent, which is also where R's integer NA
    // sentinel lands.
    bool isMissing(const TypeFull& type, OperandId operand) const
    {
        switch (type.type)
        {
            case Type::Double: return std::isnan(this->numberAt(operand));
            case Type::Datey: return !isValidDatey(this->clicksAt(operand));
            case Type::Durationy: return !isValidDurationy(this->clicksAt(operand));
            default: break;
        }
        throw std::runtime_error("veil: `is_na` does not apply to this type.");
    }

    // A scalar select copies whichever branch the condition picked. Both branches were computed --
    // the branch-omitting form of if/else belongs to the vector phase, where skipping a branch saves
    // a whole time vector's work rather than one register move.
    void executeSelection(const Instruction& instruction)
    {
        const OperandId chosen = this->logicalAt(instruction.args[0]) ? instruction.args[1] : instruction.args[2];
        ScalarValue& target = this->registers.at(instruction.result);
        const ScalarValue& source = this->registers.at(chosen);
        if (target.index() != source.index())
        {
            throw std::runtime_error("veil: a select's branch does not match the type of its result.");
        }
        target = source;
    }

    // ToDouble is source-aware, which is why the operand carries its type: a click-backed value reads
    // as a number of years, and a logical as zero or one.
    void executeConversion(const Instruction& instruction)
    {
        if (instruction.op != Op::ToDouble)
        {
            detail::opNotRunnable(instruction.op, "needs the time vector, which this interpreter does not run.");
        }

        const TypeFull sourceType = this->block.operandAt(instruction.args[0]).type;
        switch (sourceType.type)
        {
            case Type::Double:
                this->store(instruction.result, this->numberAt(instruction.args[0]));
                return;
            case Type::Datey:
            case Type::Durationy:
                this->store(instruction.result, yearsFromClicks(this->clicksAt(instruction.args[0])));
                return;
            case Type::Bool:
                this->store(instruction.result, this->logicalAt(instruction.args[0]) ? 1.0 : 0.0);
                return;
            default:
                break;
        }
        throw std::runtime_error("veil: a value of this type does not convert to a number.");
    }

    const Block& block;

    // The register file, indexed by OperandId, and the time-vector buffers, indexed by BufferId.
    // A vector operand's scalar register is never read and a scalar operand has no buffer, which is
    // what the operand's shape decides. A buffer is resized rather than reallocated per individual,
    // so the allocations settle after the first few.
    //
    // WHICH BUFFER AN OPERAND USES IS THE BLOCK'S DECISION, not this file's: passAllocateVectorBuffers
    // works out where live ranges do not overlap and several operands can share one store. An
    // unallocated block gives every vector operand a buffer of its own, so this runs either way.
    std::vector<ScalarValue> registers;
    std::vector<std::vector<double>> buffers;

    TimeGrid grid;
    size_t slotsFilled = 0;
};

} // namespace veil
