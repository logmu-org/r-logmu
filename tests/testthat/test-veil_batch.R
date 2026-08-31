# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests that running several specifications together gives each of them its own answer, unchanged.
#
# THIS IS THE OP-LEVEL COUNTERPART OF test-veil_threads.R, and it earns the same exactness for the
# same reason: a batch is a scheduling change and nothing else. Each (block, chunk) pair owns its own
# partial and the partials fold per block in chunk order, so batching must not move a single digit
# either from batching itself or from the thread count it is run at. Both claims are therefore
# `expect_identical` rather than a tolerance.
#
# WHAT EACH ASSERTION HERE CAN AND CANNOT CATCH:
#
#   - Comparing each entry against the same specification run ALONE through cpp_veil_aev catches the
#     block index being mixed up, a partial written to the wrong block's accumulator, and a fold that
#     reads the wrong slice. It cannot catch batching silently running one block repeatedly, which is
#     why the three specifications below deliberately produce three different answers.
#   - The thread-count comparison catches a total that depends on the schedule once there is more than
#     one block to schedule, which is the case that did not exist before this entry point.
#   - `task_count` is the only observable that the (block, chunk) flattening happened at all. It
#     witnesses the task numbering, not that anything ran in parallel.
#
# THE SPECIFICATIONS MUST DISAGREE WITH EACH OTHER. If all three asked the same question, every
# assertion below would pass just as happily with the block index hard-wired to zero. They differ in
# mortality, in weight and in include, so each of the three ways a spec can vary is covered.

clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

records_per_chunk <- cpp_veil_record_chunks(1L)$records_per_chunk

# FOUR CHUNKS over three blocks, so the flat task list has twelve entries and there is more work than
# threads in every case below.
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

# Three questions that genuinely differ: a weighted A/E/V, an unweighted one at a different level of
# mortality, and a weighted one restricted to an age band.
specs <- list(
  list(
    mortality      = it_obj(mortality_const(log_mu = -3.2)),
    weight         = it_ast(~ .i$amount),
    include        = NULL,
    overdispersion = no_overdispersion
  ),
  list(
    mortality      = it_obj(mortality_const(log_mu = -2.5)),
    weight         = NULL,
    include        = NULL,
    overdispersion = no_overdispersion
  ),
  list(
    mortality      = it_obj(mortality_const(log_mu = -3.2)),
    weight         = it_ast(~ .i$amount),
    include        = age(30, 60),
    overdispersion = no_overdispersion
  )
)

alone <- function(spec, threads = 1L, columns = wide) {
  cpp_veil_aev(spec$mortality, spec$weight, columns, quarter_scale, spec$include,
               spec$overdispersion, threads)
}

batch_at <- function(threads, columns = wide, keep_contributions = TRUE) {
  cpp_veil_run(specs, columns, quarter_scale, keep_contributions, threads)
}

test_that("the three specifications ask genuinely different questions", {
  # THE PRECONDITION FOR EVERY OTHER TEST IN THIS FILE. If these ever coincide, the comparisons below
  # stop being able to see a block index that is wrong, and they would keep passing while doing so.
  singles <- lapply(specs, alone)

  expect_false(singles[[1]]$E == singles[[2]]$E)
  expect_false(singles[[1]]$E == singles[[3]]$E)
  expect_false(singles[[2]]$E == singles[[3]]$E)

  # The third is a strict subset of the first's exposure, which is the include doing its work.
  expect_lt(singles[[3]]$records_included, singles[[1]]$records_included)
})

test_that("each batch entry is bit-identical to the same specification run alone", {
  batch <- batch_at(1L)

  expect_equal(batch$spec_count, 3L)
  expect_equal(length(batch$results), 3L)

  for (i in seq_along(specs)) {
    single <- alone(specs[[i]])
    entry <- batch$results[[i]]

    expect_identical(entry$A, single$A)
    expect_identical(entry$E, single$E)
    expect_identical(entry$V, single$V)

    # Per individual as well as in total, so a failure localises to a chunk rather than leaving a
    # total that is merely wrong by an unknown amount.
    expect_identical(entry$contributions$A, single$contributions$A)
    expect_identical(entry$contributions$E, single$contributions$E)
    expect_identical(entry$contributions$V, single$contributions$V)

    expect_identical(entry$records_included, single$records_included)
    expect_identical(entry$slot_evaluations, single$slot_evaluations)
  }
})

test_that("the batch is bit-identical on two threads", {
  # Runs everywhere, CRAN included. Two workers over twelve tasks is already enough for the schedule
  # to decide the order of every one of the three sums.
  single <- batch_at(1L)
  threaded <- batch_at(2L)

  expect_equal(single$chunk_count, 4L)
  expect_equal(single$task_count, 12L)

  for (i in seq_along(specs)) {
    expect_identical(threaded$results[[i]]$A, single$results[[i]]$A)
    expect_identical(threaded$results[[i]]$E, single$results[[i]]$E)
    expect_identical(threaded$results[[i]]$V, single$results[[i]]$V)
    expect_identical(threaded$results[[i]]$contributions$E, single$results[[i]]$contributions$E)
    expect_identical(threaded$results[[i]]$records_included, single$results[[i]]$records_included)
    expect_identical(threaded$results[[i]]$slot_evaluations, single$results[[i]]$slot_evaluations)
  }
})

test_that("the batch is bit-identical at higher thread counts", {
  skip_on_cran()

  single <- batch_at(1L)

  for (threads in c(3L, 4L, 8L)) {
    threaded <- batch_at(threads)
    for (i in seq_along(specs)) {
      expect_identical(threaded$results[[i]]$A, single$results[[i]]$A)
      expect_identical(threaded$results[[i]]$E, single$results[[i]]$E)
      expect_identical(threaded$results[[i]]$V, single$results[[i]]$V)
    }
  }
})

test_that("repeating the same threaded batch gives the same answer every time", {
  skip_on_cran()

  # WHAT CATCHES A GENUINE RACE across blocks specifically: two tasks from different blocks sharing an
  # accumulator would show as answers that vary between runs of an identical request rather than as a
  # wrong answer, and a single comparison could pass by luck. Eight threads over twelve tasks so that
  # workers contend for the counter.
  runs <- lapply(seq_len(10L), function(i) batch_at(8L))

  for (run in runs) {
    for (i in seq_along(specs)) {
      expect_identical(run$results[[i]]$A, runs[[1]]$results[[i]]$A)
      expect_identical(run$results[[i]]$E, runs[[1]]$results[[i]]$E)
      expect_identical(run$results[[i]]$V, runs[[1]]$results[[i]]$V)
    }
  }
})

test_that("the task list is the product of blocks and chunks", {
  # The only observable that the flattening happened. A batch of three over four chunks is twelve
  # tasks, and it is the task count -- not the chunk count -- that the pool is handed.
  batch <- batch_at(1L)
  expect_equal(batch$task_count, batch$spec_count * batch$chunk_count)
})

test_that("a thread count is honoured", {
  expect_equal(batch_at(1L)$threads_used, 1L)
  expect_equal(batch_at(2L)$threads_used, 2L)
})

test_that("contributions can be turned off", {
  # They cost one double per output per individual per specification, so a real batch will not want
  # them. Off means absent rather than zeroed, so that nothing reads them by accident.
  batch <- batch_at(1L, keep_contributions = FALSE)

  expect_equal(length(batch$results[[1]]$contributions$A), 0L)
  expect_equal(length(batch$results[[1]]$contributions$E), 0L)

  # The totals are unaffected by whether the diagnostic was kept.
  expect_identical(batch$results[[1]]$A, batch_at(1L)$results[[1]]$A)
})

test_that("more threads than tasks is harmless", {
  # A small dataset with a large thread count is exactly what an interactive user tries first. Safe
  # under a check despite the sixteen: one chunk per block over three blocks is three tasks, which
  # still spawns, so this exercises the ordinary path rather than the no-thread branch.
  small <- lapply(wide, function(column) column[seq_len(50L)])

  many <- batch_at(16L, small)
  one <- batch_at(1L, small)

  expect_equal(one$chunk_count, 1L)
  expect_equal(one$task_count, 3L)
  for (i in seq_along(specs)) {
    expect_identical(many$results[[i]]$A, one$results[[i]]$A)
    expect_identical(many$results[[i]]$E, one$results[[i]]$E)
  }
})

test_that("overdispersion is read per specification, not once for the batch", {
  # A batch entry may legitimately disagree with its neighbours about overdispersion -- that is what
  # lets a user override it for one entry. THE SOLE WITNESS that the value is read per spec: hoist it
  # out of the loop and both entries below would come back with the same V, which is exactly what a
  # batch-wide reading would produce and exactly what nobody would question.
  pair <- list(
    specs[[1]],
    list(
      mortality      = specs[[1]]$mortality,
      weight         = specs[[1]]$weight,
      include        = specs[[1]]$include,
      overdispersion = 4
    )
  )
  batch <- cpp_veil_run(pair, wide, quarter_scale, FALSE, 1L)

  # Same question but for the dispersion, so A and E must be bit-identical and only V may move.
  expect_identical(batch$results[[2]]$A, batch$results[[1]]$A)
  expect_identical(batch$results[[2]]$E, batch$results[[1]]$E)
  expect_equal(batch$results[[2]]$V, batch$results[[1]]$V * 4, tolerance = 1e-14)
})

test_that("an empty batch runs nothing and answers nothing", {
  batch <- cpp_veil_run(list(), wide, quarter_scale, TRUE, 1L)

  expect_equal(batch$spec_count, 0L)
  expect_equal(length(batch$results), 0L)
  expect_equal(batch$task_count, 0L)
})

test_that("a specification without a mortality is refused, and says which one", {
  # Named by position, because a batch is built by a loop and "one of them is wrong" is not actionable.
  bad <- list(specs[[1]], list(weight = it_ast(~ .i$amount)))
  expect_error(cpp_veil_run(bad, wide, quarter_scale, TRUE, 1L), "Specification 2")
})

test_that("a specification without an overdispersion is refused, and says which one", {
  # OVERDISPERSION IS MANDATORY AND HAS NO DEFAULT ANYWHERE, so a spec that omits it must fail rather
  # than fall back to one. A silent default of 1 is the whole failure this rule exists to prevent:
  # every confidence interval and residual an aev reports would come back too tight by sqrt(Omega),
  # and nothing about the numbers would look wrong.
  bad <- list(specs[[1]], specs[[2]][c("mortality", "weight", "include")])
  expect_error(cpp_veil_run(bad, wide, quarter_scale, TRUE, 1L), "Specification 2")
})

test_that("a non-positive overdispersion is refused", {
  # NaN alongside zero and the negatives: the check is written as `!(x > 0)` precisely so that a NaN
  # fails it rather than slipping through a `x <= 0` test, which NaN passes.
  for (bad_value in c(0, -1, NaN)) {
    bad <- list(list(
      mortality      = it_obj(mortality_const(log_mu = -3.2)),
      overdispersion = bad_value
    ))
    expect_error(cpp_veil_run(bad, wide, quarter_scale, TRUE, 1L), "must be a positive number")
  }
})

test_that("one bad specification fails the batch before any of it runs", {
  # Compiled up front, so the failure does not depend on where the bad one sits in the list.
  bad_first <- list(list(mortality = it_obj(mortality_const(log_mu = -3.2)),
                         include = "not an include",
                         overdispersion = no_overdispersion),
                    specs[[1]])
  expect_error(cpp_veil_run(bad_first, wide, quarter_scale, TRUE, 1L),
               "must be an `include` object")
})

test_that("a negative thread count is refused", {
  # Refused rather than clamped, before anything is compiled or spawned.
  expect_error(batch_at(-1L), "cannot be negative")
})
