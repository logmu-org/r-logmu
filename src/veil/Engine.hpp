// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <vector>
#include "veil/Block.hpp"
#include "veil/ColumnView.hpp"
#include "veil/Interpreter.hpp"
#include "veil/RecordChunk.hpp"
#include "veil/ThreadPool.hpp"
#include "veil/Type.hpp"

namespace veil
{

// Running a block over a whole dataset and adding up what comes out.
//
// THE ANSWER IS THE TOTAL, AND IT IS COMPUTED HERE. An A/E/V is one A, one E and one V, not three
// vectors for the host to add up: the summation is part of the calculation, it belongs in the
// platform-independent core so that every front end gets the same numbers, and a per-individual
// array would cost one double per output per life for no purpose.
//
// THE ORDER OF THE SUM IS FIXED BY THE DATA ALONE. Individuals are added within a chunk in record
// order, and the chunk partials are folded in chunk order, so nothing about scheduling, thread count
// or machine load can reach the result. See RecordChunk.hpp for why that is arranged this way and
// what it costs. Threading changes only which worker walks which chunk -- not a single addition, and
// so not a single digit. That is what makes the threading test the strongest one available: the same
// data at any thread count, compared exactly rather than to a tolerance.
//
// EVERY REPORTED NUMBER IS ACCUMULATED PER CHUNK, not just the totals. `recordsIncluded` and
// `slotEvaluations` are counts, so integer addition would let them be summed in any order at all --
// but they are diagnostics that tests assert on, and a diagnostic that is reproducible for a
// different reason from the answer is one an unrelated change can quietly break. They fold in chunk
// order alongside the doubles, for one size_t per chunk.
//
// EVERY OUTPUT IS A SUM. That holds for A, E and V, for a log-likelihood and for its first and second
// derivatives. A maximum or a minimum would accumulate just as safely, since both are associative --
// and safely in veil in particular, because Interpreter.hpp uses IEEE 754-2019 `minimum` and
// `maximum` rather than C's `fmin` and `fmax`, whose 2008 ancestors discarded NaN and so were not
// associative. Something like a count of distinct values would not fit, and is meant to be a
// non-threaded data check rather than a calculation output.

// What one run of a block over a dataset produces.
struct CalculationResult final
{
    // One per block output, in output order. This is the answer.
    std::vector<double> totals;

    // DIAGNOSTIC ONLY, and empty unless asked for. `records * outputCount`, record-major, so
    // individual i's output j is at `i * outputCount + j`. An excluded individual reads zero.
    //
    // This exists for the tests rather than for a caller: the analytic oracles for the midpoint rule
    // are written per individual, and a total that disagrees with them says only that something is
    // wrong. It is also the oracle for the accumulation itself, since a plain sum of these must
    // agree with the chunked total, which catches a chunk folded twice or a boundary out by one in a
    // way an analytic total would not. It goes away when the real batch entry point lands.
    std::vector<double> contributions;

    size_t chunks = 0;          // How many partials were folded.
    size_t recordsIncluded = 0; // Individuals that actually contributed.
    size_t slotEvaluations = 0; // Time-vector slots filled, over every instruction and individual.

    // True when the host asked to stop before every chunk had run. THE TOTALS ARE THEN MEANINGLESS --
    // some chunks contributed and some did not -- so a caller must raise rather than report them. It
    // is not an error in itself, which is why it is a field rather than an exception: the host asked.
    bool interrupted = false;
};

// What one chunk contributed, beside its output partials.
struct ChunkTally final
{
    size_t recordsIncluded = 0;
    size_t slotEvaluations = 0;
};

// Runs ONE block over ONE chunk of the records, adding what it finds to `partial` and `tally`.
//
// THIS IS THE ONLY PLACE A RECORD IS WALKED, and both entry points below call it. They differ in how
// they number their tasks and in how they fold the partials afterwards; what happens inside a chunk
// is the same work either way, and a copy of it in each would be two things to keep in step.
//
// `partial` points at this chunk's own `outputCount` accumulators, and `contributions` is null unless
// the caller wants the diagnostic per-individual values. Neither is read by anyone else while this
// runs -- see the entry points for why that is true by construction rather than by convention.
inline void runOneChunk(
    const Block& block,
    const std::vector<const ColumnView*>& columns,
    size_t records,
    size_t chunk,
    double* const partial,
    ChunkTally& tally,
    std::vector<double>* const contributions)
{
    const size_t outputCount = block.outputs().size();
    const RecordChunk range = chunkOf(records, chunk);

    Interpreter interpreter(block);
    for (size_t record = range.startIndex; record < range.endIndex; ++record)
    {
        // An excluded individual contributes the identity, so there is nothing to add and nothing to
        // run.
        if (!interpreter.loadRecord(columns, record)) { continue; }

        interpreter.run();
        ++tally.recordsIncluded;

        for (size_t output = 0; output < outputCount; ++output)
        {
            const double value = interpreter.outputAsNumber(output);
            partial[output] += value;
            if (contributions != nullptr)
            {
                (*contributions)[record * outputCount + output] = value;
            }
        }
    }

    // Cumulative since this interpreter was built, which is this chunk and no other.
    tally.slotEvaluations = interpreter.slotEvaluations();
}

// Runs `block` over `records` individuals and returns the totals.
//
// `columns` is indexed by the ColumnId the block's field references carry, exactly as
// Interpreter::loadRecord expects.
//
// `threadCount` follows resolveThreadCount: 1 runs on the calling thread and spawns nothing, 0 asks
// for as many threads as the machine reports, anything else is taken literally. IT CANNOT MOVE AN
// ANSWER -- see the note on summation order above -- which is the only reason it is safe to expose.
//
// `hostPoll` is called only by the calling thread, only between chunks, and answers true when the
// host wants to stop. It may be empty. See ThreadPool.hpp for why the core takes a callback here
// rather than knowing what an interrupt is.
//
// NOTHING BELOW IS SHARED MUTABLE STATE, and that is not luck. Operands are block-local by design;
// the register file, the time grid and the vector buffers all live inside an Interpreter, one of
// which is built per chunk; the Block is const for the whole run and owns its mortality table by
// copy; and the column views are read-only pointers into the host's memory. The two places a chunk
// writes -- its slice of `partials` and its own records' slots in `contributions` -- are disjoint by
// construction, because chunks partition the records.
//
// THE INTERPRETER IS BUILT PER CHUNK RATHER THAN PER WORKER. It costs a register file and a buffer
// resize, amortised over RecordsPerChunk individuals, which is nothing; and it keeps the pool free of
// any per-worker plumbing, so a task index is all a body ever needs. If profiling ever says
// otherwise, the change is local to this function.
inline CalculationResult runCalculation(
    const Block& block,
    const std::vector<const ColumnView*>& columns,
    size_t records,
    bool keepContributions,
    size_t threadCount = 1,
    const std::function<bool()>& hostPoll = {})
{
    const std::vector<OperandId>& outputs = block.outputs();
    const size_t outputCount = outputs.size();

    // Guarded rather than assumed, because accumulating a click count or a logical would answer with
    // a plausible number rather than failing. A block whose outputs are not numbers is not a
    // calculation, and the per-individual entry points are where it should be run.
    for (const OperandId output : outputs)
    {
        if (block.operandAt(output).type.type != Type::Double)
        {
            throw std::runtime_error("veil: only a numeric output can be accumulated over a "
                                     "dataset.");
        }
    }

    CalculationResult result;
    result.totals.assign(outputCount, 0.0);
    result.chunks = chunkCount(records);
    if (keepContributions) { result.contributions.assign(records * outputCount, 0.0); }

    // One accumulator per chunk, not one per worker. That is the whole reproducibility argument, and
    // it costs a few doubles per chunk.
    std::vector<double> partials(result.chunks * outputCount, 0.0);
    std::vector<ChunkTally> tallies(result.chunks);

    // ONE CHUNK'S WORK, which is `runOneChunk` above -- all this task numbering has to say is which
    // chunk.
    const auto runChunk = [&](size_t chunk)
    {
        runOneChunk(
            block,
            columns,
            records,
            chunk,
            partials.data() + chunk * outputCount,
            tallies[chunk],
            keepContributions ? &result.contributions : nullptr);
    };

    const ParallelRunOutcome outcome =
        runInParallel(result.chunks, threadCount, hostPoll, runChunk);

    result.interrupted = !outcome.completed;
    if (result.interrupted) { return result; }

    // Folded in CHUNK ORDER. Written as the outer loop over chunks so that the order is visible in
    // the shape of the code rather than only in a comment.
    for (size_t chunk = 0; chunk < result.chunks; ++chunk)
    {
        for (size_t output = 0; output < outputCount; ++output)
        {
            result.totals[output] += partials[chunk * outputCount + output];
        }
        result.recordsIncluded += tallies[chunk].recordsIncluded;
        result.slotEvaluations += tallies[chunk].slotEvaluations;
    }

    return result;
}

// Runs several blocks over ONE dataset and returns one result per block, in block order.
//
// THE UNIT OF DISPATCH IS (BLOCK, CHUNK), AND THAT IS THE POINT OF THIS FUNCTION. Op-level
// parallelism is the primary axis -- a user runs several A/E/Vs and several fits together, so there
// are usually more operations than cores -- but until now a batch lowered to one block and there was
// nothing to schedule across. Flattening the (block, chunk) pairs into one task index is the whole
// mechanism, and `runInParallel` already takes exactly that shape, so the pool needs no change.
//
// EVERY BLOCK SEES THE SAME DATASET, so every block divides into the same chunks and the flattening
// is plain division: task `t` is block `t / chunks`, chunk `t % chunks`. Several datasets would need
// a per-block task offset instead, and that is a later slice.
//
// SCHEDULING STILL CANNOT MOVE AN ANSWER, for the same reason one level down. Each (block, chunk)
// pair owns its own partial, and the partials are folded per block in chunk order once every task has
// run, so which worker took which pair reaches nothing. The bit-identity test therefore applies to a
// batch exactly as it applies to a single calculation.
//
// AN EXCEPTION SURFACES BY TASK INDEX, which here means the lowest-numbered block and then its
// lowest-numbered chunk. That is a property of the caller's own numbering rather than of the
// schedule, so a failing batch reports the same error every time.
//
// NOTHING BELOW IS SHARED MUTABLE STATE. The argument is the one on `runCalculation`, extended by one
// dimension: two tasks with the same block write disjoint record ranges, and two tasks with different
// blocks write different result objects entirely. No vector is resized once the run has started.
//
// AN INTERRUPT VOIDS THE WHOLE BATCH rather than the calculation that happened to be running. Some
// chunks contributed and some did not, and that is true of every block at once, so every result is
// marked and the caller must raise.
inline std::vector<CalculationResult> runCalculations(
    const std::vector<const Block*>& blocks,
    const std::vector<const ColumnView*>& columns,
    size_t records,
    bool keepContributions,
    size_t threadCount = 1,
    const std::function<bool()>& hostPoll = {})
{
    const size_t blockCount = blocks.size();
    std::vector<CalculationResult> results(blockCount);
    if (blockCount == 0) { return results; }

    // Checked for EVERY block before any thread is spawned, so that a batch containing one bad block
    // fails the same way whether or not the others would have run.
    for (const Block* const block : blocks)
    {
        if (block == nullptr) { throw std::runtime_error("veil: a null block cannot be run."); }

        for (const OperandId output : block->outputs())
        {
            if (block->operandAt(output).type.type != Type::Double)
            {
                throw std::runtime_error("veil: only a numeric output can be accumulated over a "
                                         "dataset.");
            }
        }
    }

    const size_t chunks = chunkCount(records);

    // One accumulator per (block, chunk), which is the reproducibility argument of `runCalculation`
    // with a block index in front of it.
    std::vector<std::vector<double>> partials(blockCount);
    std::vector<std::vector<ChunkTally>> tallies(blockCount);

    for (size_t blockIndex = 0; blockIndex < blockCount; ++blockIndex)
    {
        const size_t outputCount = blocks[blockIndex]->outputs().size();
        CalculationResult& result = results[blockIndex];

        result.totals.assign(outputCount, 0.0);
        result.chunks = chunks;
        if (keepContributions) { result.contributions.assign(records * outputCount, 0.0); }

        partials[blockIndex].assign(chunks * outputCount, 0.0);
        tallies[blockIndex].resize(chunks);
    }

    // ONE BLOCK'S WORK OVER ONE CHUNK, which is the same `runOneChunk` the single-block path uses.
    // Only the task numbering is new. The division is safe because the body runs only when there is
    // at least one task, and a dataset with no records has no chunks and so none.
    const auto runBlockChunk = [&](size_t task)
    {
        const size_t blockIndex = task / chunks;
        const size_t chunk = task % chunks;

        const Block& block = *blocks[blockIndex];
        CalculationResult& result = results[blockIndex];

        runOneChunk(
            block,
            columns,
            records,
            chunk,
            partials[blockIndex].data() + chunk * block.outputs().size(),
            tallies[blockIndex][chunk],
            keepContributions ? &result.contributions : nullptr);
    };

    const ParallelRunOutcome outcome =
        runInParallel(blockCount * chunks, threadCount, hostPoll, runBlockChunk);

    if (!outcome.completed)
    {
        for (CalculationResult& result : results) { result.interrupted = true; }
        return results;
    }

    // Folded per block in CHUNK ORDER, so that each result is bit-identical to the same block run on
    // its own.
    for (size_t blockIndex = 0; blockIndex < blockCount; ++blockIndex)
    {
        CalculationResult& result = results[blockIndex];
        const size_t outputCount = blocks[blockIndex]->outputs().size();

        for (size_t chunk = 0; chunk < chunks; ++chunk)
        {
            for (size_t output = 0; output < outputCount; ++output)
            {
                result.totals[output] += partials[blockIndex][chunk * outputCount + output];
            }

            result.recordsIncluded += tallies[blockIndex][chunk].recordsIncluded;
            result.slotEvaluations += tallies[blockIndex][chunk].slotEvaluations;
        }
    }

    return results;
}

} // namespace veil
