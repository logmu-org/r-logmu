// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>

namespace veil
{

// How a dataset is divided into chunks for accumulation.
//
// WHY A DATASET IS DIVIDED AT ALL. A calculation's answer is a sum over individuals, and floating
// point addition is not associative, so a sum is only as reproducible as its order. Workers pulling
// work off a queue and accumulating per THREAD would give a total that depends on which chunk landed
// on which worker, and so on timing -- the same binary, on the same machine, over the same data,
// disagreeing with itself between runs. That is the one thing that must not happen. Accumulating per
// CHUNK and folding the chunk partials in chunk order removes the schedule from the answer
// altogether, while leaving the scheduling itself free to be as dynamic as it likes.
//
// SO THE DIVISION HAPPENS WHETHER OR NOT ANYTHING IS THREADED. A single-threaded run walks the same
// chunks in the same order and folds the same partials, which is what lets threading be added later
// without moving a number, and what lets a user turn threading off to check something and get the
// same answer back.
//
// THE PRIMARY UNIT OF PARALLELISM IS THE OPERATION, NOT THE CHUNK (Tim, 2026-07-29). A user
// typically runs several A/E/Vs and several model fits together, so there are usually more
// operations than cores, and each operation summing its own records is both simpler and enough.
// Chunking is mitigation for the case that leaves out: a few calculations over a very large dataset,
// and the single model fit, whose iteration runs one operation at a time.

// The largest number of records in a chunk.
//
// A MAXIMUM RATHER THAN A MINIMUM: the division below rounds the chunk COUNT up, so 10,001 records
// become two chunks of about five thousand rather than one chunk of 10,001.
//
// THIS NUMBER IS PART OF WHAT AN ANSWER MEANS. A total is the sum of the chunk partials in chunk
// order, so changing it moves the last digits of every result. It is a compile-time constant so that
// changing it is a deliberate change to the code rather than something that can vary with the
// machine, and it must never be derived from core count or available memory.
//
// STILL TO BE TUNED ON TEST DATA, and it may well want to differ between kinds of calculation, since
// a log-likelihood with p parameters costs far more per individual than an A/E/V.
constexpr size_t RecordsPerChunk = 10000;

// How many chunks a dataset of `records` individuals divides into. None at all when there are no
// records, so a driving loop simply does not run.
constexpr size_t chunkCount(size_t records) noexcept
{
    return (records + RecordsPerChunk - 1) / RecordsPerChunk;
}

// A half-open range of record indices.
struct RecordChunk final
{
    size_t startIndex = 0;
    size_t endIndex = 0;

    constexpr size_t size() const noexcept { return this->endIndex - this->startIndex; }
};

// The records belonging to one chunk.
//
// THE PROPORTIONAL FORM, chosen deliberately over the more common remainder-to-the-front version
// that numpy's `array_split` and the MPI scatter helpers use. Coverage is exact by construction --
// chunk i's `endIndex` is the same expression as chunk i+1's `startIndex` -- so there is no gap or
// overlap to reason about or to test for, and no remainder arithmetic to get wrong. Sizes still
// differ by at most one, because consecutive terms of `floor(i * records / chunks)` differ by either
// the floor or the ceiling of `records / chunks`.
//
// The two forms are NOT interchangeable: 25,000 records in three chunks divide 8333 / 8333 / 8334
// here, and 8334 / 8333 / 8333 the other way. Same sizes, different membership, different last
// digits in the total.
//
// Equal parts matter most when there are fewest of them. Fifteen thousand records divide into 7500
// and 7500 rather than 10000 and 5000, which is the difference between two workers finishing
// together and one of them waiting.
//
// GUIDED OR GEOMETRIC CHUNKING IS REFUSED, though it is what a parallel-loop textbook would suggest:
// it balances better for less overhead, and it makes the partition a function of the schedule, which
// is precisely what this exists to prevent.
//
// An index past the end answers an EMPTY range rather than reading off the end, which also covers a
// dataset of no records at all -- the arithmetic below would divide by zero for it.
constexpr RecordChunk chunkOf(size_t records, size_t index) noexcept
{
    const size_t chunks = chunkCount(records);
    if (index >= chunks) { return RecordChunk{records, records}; }

    return RecordChunk{index * records / chunks, (index + 1) * records / chunks};
}

// The division rule is part of what an answer means, so it is pinned at compile time in the same
// spirit as `everyTimeScaleIsWholeClicks()`: a change to it fails the build rather than quietly
// moving the numbers.
static_assert(chunkCount(0) == 0);
static_assert(chunkCount(1) == 1);
static_assert(chunkCount(RecordsPerChunk) == 1);
static_assert(chunkCount(RecordsPerChunk + 1) == 2);

// Tim's worked example.
static_assert(chunkOf(15000, 0).startIndex == 0 && chunkOf(15000, 0).size() == 7500);
static_assert(chunkOf(15000, 1).startIndex == 7500 && chunkOf(15000, 1).endIndex == 15000);

// An uneven division: the parts differ by one, and the extra record falls at the END. This is what
// separates the proportional form from the remainder-to-the-front form, which would put it first --
// see test-veil_chunks.R, which pins the same distinction from R over more record counts.
static_assert(chunkOf(25000, 0).size() == 8333);
static_assert(chunkOf(25000, 1).size() == 8333);
static_assert(chunkOf(25000, 2).size() == 8334);
static_assert(chunkOf(25000, 2).endIndex == 25000);

// An odd record count over two chunks, which the other form would divide 7501 then 7500.
static_assert(chunkOf(15001, 0).size() == 7500);
static_assert(chunkOf(15001, 1).size() == 7501);

// Empty, rather than undefined, off the end and for an empty dataset.
static_assert(chunkOf(15000, 2).size() == 0);
static_assert(chunkOf(0, 0).size() == 0);

} // namespace veil
