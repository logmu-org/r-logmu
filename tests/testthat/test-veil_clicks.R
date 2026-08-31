# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests that clicks are 64-bit inside the engine.
#
# A datey and a durationy each arrive from R as a 32-bit integer, and each is valid across its own
# declared range -- dates over [1000, 3000] and durations over +/- 2000 years. Those ranges overlap
# enough that a perfectly valid date plus a perfectly valid duration can exceed what 32 bits hold:
# year 2900 is 1,549,644,000 clicks and 1500 years is 801,540,000, and their sum passes INT_MAX. That
# used to be refused mid-calculation. It is now simply computed, because a click is an int64 from the
# moment it is loaded until the moment it goes back to R.
#
# The check did not disappear, it MOVED. It now sits at the one boundary where 32 bits are forced --
# an output crossing back to R -- which is where the answer is meant to be a datey R can hold, and
# where whoever asked for it can do something about it.

clicks_per_year <- 534360L
int_max <- 2147483647

cols <- list(
  birth = datey::datey(c(2900, 2000, 1950)),
  gap   = datey::durationy(c(1500, 10, 20))
)

birth_clicks <- unclass(cols$birth)
gap_clicks <- unclass(cols$gap)

test_that("the first individual really does overflow 32 bits", {
  # If this stops being true the tests below stop testing anything, so it is asserted rather than
  # assumed -- and asserted against R's own arithmetic in doubles, which has room for it.
  expect_gt(as.double(birth_clicks[1]) + as.double(gap_clicks[1]), int_max)
})

test_that("an intermediate past 32 bits is computed, not refused", {
  # `(birth + gap) - gap` is birth again. The intermediate passes INT_MAX for the first individual,
  # so this is the case that used to throw; the answer it gives now is the one it always should have.
  res <- cpp_veil_eval(it_ast(~ (.i$birth + .i$gap) - .i$gap), cols)

  expect_equal(res$type, "datey")
  expect_equal(res$values, birth_clicks)
})

test_that("a wide intermediate does not disturb a narrower one", {
  # The same expression for individuals whose intermediate fits comfortably, so a change of width
  # cannot have altered the ordinary case.
  res <- cpp_veil_eval(it_ast(~ (.i$birth + .i$gap) - .i$gap), cols)
  expect_equal(res$values[2:3], birth_clicks[2:3])
})

test_that("subtraction across the calendar is exact", {
  res <- cpp_veil_eval(it_ast(~ .i$birth - .i$gap), cols)

  expect_equal(res$type, "datey")
  expect_equal(res$values, birth_clicks - gap_clicks)
})

test_that("a result R cannot hold is refused at the boundary, naming what went wrong", {
  expect_error(
    cpp_veil_eval(it_ast(~ .i$birth + .i$gap), cols),
    "outside the range R can hold"
  )
})

test_that("min and max still select without any arithmetic at all", {
  # These cannot overflow whatever the width, because they return one of their arguments. Worth a
  # test because the include clip is built from them.
  res <- cpp_veil_eval(it_ast(~ min(.i$birth, .i$birth - .i$gap)), cols)
  expect_equal(res$values, pmin(birth_clicks, birth_clicks - gap_clicks))
})
