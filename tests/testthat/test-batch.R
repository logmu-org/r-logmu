# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests for `batch()`.
#
# THE STRONGEST ORACLE IS THE STANDALONE CALL. A batched analysis must return
# exactly what the same analysis returns on its own -- identically, not nearly.
# Everything `batch()` does is scheduling, so any difference in the answer is a
# defect by definition, and that comparison is what most of these assert.

data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940)),
    pension   = c(5000, 12000, 30000, 8000),
    male      = c(TRUE, FALSE, TRUE, FALSE),
    E2R_start = datey::datey(rep(2015, 4)),
    E2R_end   = datey::datey(c(2020, 2020, 2018, 2020)),
    E2R_died  = c(FALSE, FALSE, TRUE, FALSE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

other_data <- exp_data(
  list(
    birth     = datey::datey(1950),
    E2R_start = datey::datey(2015),
    E2R_end   = datey::datey(2020),
    E2R_died  = FALSE
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

light <- mortality_const(log_mu = -4.5)
heavy <- mortality_const(log_mu = -4.0)

# Varying with age, so that the integration interval actually matters. A
# constant mortality integrates exactly at every time scale and would hide any
# difference between them.
ageing <- mortality(-10 + 0.1 * .x)

test_that("a batch returns a plain named list in the order written", {
  b <- batch(
    .exp_data = data, .overdispersion = 2,
    second = aev(mortality = heavy),
    first  = aev(mortality = light)
  )

  expect_type(b, "list")
  expect_identical(names(b), c("second", "first"))
  expect_null(attributes(b)[["class"]])
  expect_true(all(vapply(b, is_aev, logical(1L))))
})

test_that("a batched analysis is identical to the same analysis alone", {
  b <- batch(
    .exp_data = data, .overdispersion = 2, .weight = .i$pension,
    plain  = aev(mortality = light),
    banded = aev(mortality = heavy,
                 breakdown = bands(.i$pension, thresholds = 10000))
  )

  expect_identical(
    b$plain,
    aev(data, mortality = light, weight = .i$pension, overdispersion = 2)
  )
  expect_identical(
    b$banded,
    aev(data, mortality = heavy, weight = .i$pension, overdispersion = 2,
        breakdown = bands(.i$pension, thresholds = 10000))
  )
})

test_that("labels survive the trip through a batch", {
  b <- batch(
    .exp_data = data, .overdispersion = 1,
    by_amount = aev(mortality = light,
                    breakdown = bands(.i$pension, thresholds = 10000))
  )

  expect_identical(names(b$by_amount), c("< 10000", ">= 10000"))
  expect_identical(group_names(b$by_amount), rep("pension", 2L))
})

# ---- settings as defaults --------------------------------------------------

test_that("a batch setting applies where an analysis has none, and loses where it has", {
  b <- batch(
    .exp_data = data, .overdispersion = 2,
    inherited = aev(mortality = light),
    own       = aev(mortality = light, overdispersion = 1)
  )

  # V carries the overdispersion, so it is what tells the two apart: E is the
  # same run either way and V scales with it.
  expect_equal(b$inherited$E, b$own$E)
  expect_equal(b$inherited$V, 2 * b$own$V)
})

test_that("overdispersion is still required when neither level gives one", {
  expect_error(
    batch(.exp_data = data, a = aev(mortality = light)),
    "overdispersion"
  )
})

test_that("a weight default is captured in the caller's frame, not the batch's", {
  # `amounts` exists only here. Captured anywhere but the caller's frame it
  # would not be found at all.
  amounts <- variable(.i$pension)

  b <- batch(
    .exp_data = data, .overdispersion = 1, .weight = amounts,
    weighted = aev(mortality = light)
  )

  expect_identical(
    b$weighted,
    aev(data, mortality = light, weight = .i$pension, overdispersion = 1)
  )
})

test_that("include and breakdown default from the batch too", {
  b <- batch(
    .exp_data = data, .overdispersion = 1,
    .include = include(.i$male),
    .breakdown = bands(.i$pension, thresholds = 10000),
    males = aev(mortality = light)
  )

  expect_identical(
    b$males,
    aev(data, mortality = light, overdispersion = 1,
        include = include(.i$male),
        breakdown = bands(.i$pension, thresholds = 10000))
  )
})

test_that("exp_data may be given per analysis instead of once", {
  b <- batch(.overdispersion = 1, a = aev(data, mortality = light))
  expect_identical(b$a, aev(data, mortality = light, overdispersion = 1))
})

# ---- one dataset, several time scales --------------------------------------

test_that("a batch refuses a second experience dataset", {
  expect_error(
    batch(.overdispersion = 1,
          a = aev(data, mortality = light),
          b = aev(other_data, mortality = light)),
    "one experience dataset"
  )
})

test_that("analyses may differ in time_scale, each still matching its own run", {
  # The time grid is a property of the RUN, so these cannot share a crossing.
  # The point of the test is that the results are right regardless.
  b <- batch(
    .exp_data = data, .overdispersion = 1,
    yearly  = aev(mortality = ageing, time_scale = 1),
    monthly = aev(mortality = ageing, time_scale = 1 / 12)
  )

  expect_identical(b$yearly,
                   aev(data, mortality = ageing, overdispersion = 1, time_scale = 1))
  expect_identical(b$monthly,
                   aev(data, mortality = ageing, overdispersion = 1, time_scale = 1 / 12))

  # And they really are different answers, or the test above proves nothing.
  expect_false(isTRUE(all.equal(b$yearly$E, b$monthly$E)))
})

test_that("a threads argument inside a batched call is ignored, not refused", {
  b <- batch(
    .exp_data = data, .overdispersion = 1, .threads = 1,
    a = aev(mortality = light, threads = 2)
  )
  expect_identical(b$a, aev(data, mortality = light, overdispersion = 1))
})

# ---- what a batch refuses --------------------------------------------------

test_that("every analysis must be named, and no name may begin with a dot", {
  expect_error(batch(.exp_data = data, .overdispersion = 1, aev(mortality = light)),
               "must be named")
  expect_error(batch(.exp_data = data, .overdispersion = 1, a = aev(mortality = light),
                     a = aev(mortality = heavy)),
               "own name")

  # A dotted name is the reservation that lets settings be added later without
  # breaking anyone's script.
  expect_error(batch(.exp_data = data, .overdispersion = 1,
                     .mine = aev(mortality = light)),
               "leading dot")
})

test_that("a batch needs at least one analysis", {
  expect_error(batch(.exp_data = data, .overdispersion = 1), "at least one")
})

test_that("an element that is not a logmu analysis is refused, and a missing dot is named", {
  expect_error(batch(.exp_data = data, .overdispersion = 1, a = 1 + 1),
               "not a call to a logmu analysis")

  # The mistake this exists for: a setting written without its dot becomes an
  # element, and the name says what was meant.
  expect_error(batch(.exp_data = data, overdispersion = 2, a = aev(mortality = light)),
               "Did you mean `\\.overdispersion`")
})

test_that("an analysis may not refer to another in the same batch", {
  # Nothing binds `first` here, so without the check this would either fail with
  # R's own "object not found" or, worse, silently use a stale object of that
  # name from the workspace.
  expect_error(
    batch(.exp_data = data, .overdispersion = 1,
          first  = aev(mortality = light),
          second = aev(mortality = light, weight = first$E)),
    "another analysis in the same"
  )
})

test_that("an analysis with no mortality is refused by name", {
  expect_error(
    batch(.exp_data = data, .overdispersion = 1, missing_one = aev(include = include(.i$male))),
    "`missing_one` has no `mortality`"
  )
})

test_that("an analysis with no data anywhere is refused by name", {
  expect_error(
    batch(.overdispersion = 1, nowhere = aev(mortality = light)),
    "`nowhere` has no `exp_data`"
  )
})
