# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests that threading cannot move an answer.
#
# THIS IS THE STRONGEST TEST AVAILABLE TO US, and it exists only because the chunking landed
# single-threaded first. The engine already produced the final total before any thread was involved,
# so threading is a pure scheduling change: a different worker walks a chunk, and nothing else. That
# lets the assertion be EXACT EQUALITY rather than a tolerance -- `expect_identical` on doubles is
# bitwise -- which no analytic oracle could give us, since an analytic answer agrees to a tolerance
# and would happily absorb a genuine reordering.
#
# WHAT EACH ASSERTION HERE CAN AND CANNOT CATCH:
#
#   - The bit-identity comparisons catch a total that depends on the schedule, which is what
#     accumulating per thread rather than per chunk would produce. They CANNOT catch threading
#     silently not happening: a fall back to one thread passes them trivially.
#   - `threads_used` is the only observable that a requested thread count was honoured. It witnesses
#     configuration, not execution.
#   - The repeated-run test catches a genuine data race, which would show as answers that vary
#     between runs of the same request.
#
# THE TWO-CORE LIMIT SPLITS THIS FILE. A CRAN check runs on two cores and a package may not use more
# than that in its own tests, so every count above two -- and `0`, which asks for as many as the
# machine reports -- sits behind `skip_on_cran()`. What runs everywhere is the two-thread case, and
# that is not a token: two workers over four chunks is already enough for a schedule to decide the
# order of a sum, which is precisely the thing that must not matter. The higher counts buy contention
# on the task counter rather than a different claim.
#
# THERE IS DELIBERATELY NO "DID WORK ACTUALLY LAND ON MORE THAN ONE THREAD" ASSERTION. The honest
# versions of it -- counting how many chunks the calling thread took, or how many distinct thread ids
# ran a task -- are all timing-dependent: on a fast enough machine the calling thread can finish every
# chunk before a worker starts, which is correct behaviour and would fail the test. A flaky test that
# pins a scheduling accident is worse than no test, and the correctness claim is the one that matters.
#
# THE INTERRUPT PATH IS NOT TESTED HERE. It needs a real user interrupt part-way through a running
# calculation, which testthat cannot raise against itself without something far more elaborate than
# the thing it would be checking.

clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

records_per_chunk <- cpp_veil_record_chunks(1L)$records_per_chunk

# FOUR CHUNKS, so that there is more work than there are threads in every case below and a scheduler
# has something to get wrong. Individuals differ in exposure length and in weight, so a chunk boundary
# landing in the wrong place cannot be masked by them all being identical.
n <- 3L * records_per_chunk + 1L
start <- rep(2010L * clicks_per_year, n)
quarters <- rep(c(4L, 8L, 12L, 16L), length.out = n)

wide <- list(
  birth     = datey::datey(1940 + (seq_len(n) %% 30)),
  amount    = 100 + (seq_len(n) %% 7) * 250,
  E2R_start = datey_clicks(start),
  E2R_end   = datey_clicks(start + quarters * quarter),
  E2R_died  = (seq_len(n) %% 5L) == 0L
)

run_at <- function(threads, columns = wide) {
  cpp_veil_aev(it_obj(mortality_const(log_mu = -3.2)), it_ast(~ .i$amount), columns,
               quarter_scale, NULL, no_overdispersion, threads)
}

single <- run_at(1L)

# EXACT, not approximate. A tolerance here would defeat the purpose: the claim is that the arithmetic
# is the same arithmetic, not that it is close. Per individual as well as in total, so that a failure
# localises to a chunk rather than leaving a total that is merely wrong by an unknown amount.
expect_matches_single_threaded <- function(threaded) {
  expect_identical(threaded$A, single$A)
  expect_identical(threaded$E, single$E)
  expect_identical(threaded$V, single$V)

  expect_identical(threaded$contributions$A, single$contributions$A)
  expect_identical(threaded$contributions$E, single$contributions$E)
  expect_identical(threaded$contributions$V, single$contributions$V)

  # The diagnostics fold per chunk for the same reason the totals do, so they are exact as well.
  expect_identical(threaded$chunk_count, single$chunk_count)
  expect_identical(threaded$records_included, single$records_included)
  expect_identical(threaded$slot_evaluations, single$slot_evaluations)
}

test_that("the answer is bit-identical on two threads", {
  # Runs everywhere, CRAN included: two is the limit a check allows, and two workers over four chunks
  # already lets the schedule decide the order of the sum.
  expect_equal(single$chunk_count, 4L)
  expect_matches_single_threaded(run_at(2L))
})

test_that("the answer is bit-identical at higher thread counts", {
  skip_on_cran()

  for (threads in c(3L, 4L, 8L)) {
    expect_matches_single_threaded(run_at(threads))
  }
})

test_that("a thread count is honoured", {
  # The only observable that a request was not quietly ignored -- see the note at the top of the file.
  expect_equal(run_at(1L)$threads_used, 1L)
  expect_equal(run_at(2L)$threads_used, 2L)
})

test_that("zero asks the machine", {
  skip_on_cran()

  # Zero means "as many as the machine reports", which is an explicit opt-in and never a default --
  # and is why this cannot run under a check. hardware_concurrency is a hint and may answer 0, so the
  # floor is guarded rather than assumed; that guard is the whole content of this assertion, since the
  # number itself belongs to the machine.
  expect_gte(run_at(0L)$threads_used, 1L)
})

test_that("repeating the same threaded request gives the same answer every time", {
  skip_on_cran()

  # WHAT CATCHES A GENUINE RACE. Two workers writing the same accumulator, or a shared interpreter
  # register file, would show as answers that differ between runs of an identical request rather than
  # as a wrong answer -- and a single comparison against the single-threaded run could pass by luck.
  # Eight threads over four chunks so that workers contend for the counter.
  runs <- lapply(seq_len(10L), function(i) run_at(8L))

  for (run in runs) {
    expect_identical(run$A, runs[[1]]$A)
    expect_identical(run$E, runs[[1]]$E)
    expect_identical(run$V, runs[[1]]$V)
    expect_identical(run$records_included, runs[[1]]$records_included)
    expect_identical(run$slot_evaluations, runs[[1]]$slot_evaluations)
  }
})

test_that("more threads than chunks is harmless", {
  # SAFE UNDER A CHECK DESPITE THE SIXTEEN, and that is the point of the test rather than an accident:
  # a single chunk cannot be shared out, so runInParallel takes the branch that spawns no thread at
  # all. It is worth pinning because it is a different branch, and because a small dataset with a
  # large thread count is exactly what an interactive user will try first.
  small <- lapply(wide, function(column) column[seq_len(50L)])

  many <- run_at(16L, small)
  one <- run_at(1L, small)

  expect_equal(one$chunk_count, 1L)
  expect_identical(many$A, one$A)
  expect_identical(many$E, one$E)
  expect_identical(many$V, one$V)
  expect_identical(many$contributions$A, one$contributions$A)
})

test_that("a negative thread count is refused", {
  # Refused rather than clamped: a negative count is a caller's mistake, and a size_t cast would turn
  # it into an enormous one. Refused before anything is spawned, so this is safe under a check.
  expect_error(run_at(-1L), "cannot be negative")
})
