# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("constructors build includes", {
  expect_true(is_include(period(2010, 2020)))
  expect_true(is_include(age(65, 95)))
  expect_true(is_include(band(.t - .i$joined, 0, 5)))
  expect_true(is_logmu_function(period(2010, 2020)))
})

test_that("an absolute period resolves to its calendar bounds", {
  expect_equal(
    period_included(period(2010, 2020), .i = list()),
    datey::datey_interval(datey::datey(2010), datey::datey(2020))
  )
})

test_that("age() offsets the interval from birth", {
  expect_equal(
    period_included(age(65, 95), .i = list(birth = datey::datey(1950))),
    datey::datey_interval(datey::datey(2015), datey::datey(2045))
  )
})

test_that("banding `.t` minus a field offsets from that field", {
  expect_equal(
    period_included(band(.t - .i$joined, 0, 5), .i = list(joined = datey::datey(2018))),
    datey::datey_interval(datey::datey(2018), datey::datey(2023))
  )
})

test_that("band edges must increase", {
  expect_error(period(2020, 2010), "strictly increasing")
  expect_error(age(95, 65), "strictly increasing")
  expect_error(band(.t - .i$joined, 5, 0), "strictly increasing")
})

test_that("a banded variable must be increasing in `.t` at unit slope", {
  # The rule that keeps a band a single clopen interval. A decreasing expression
  # would flip it to (a, b], so adjacent bands would double-count or drop an instant.
  expect_error(band(2020 - .t, 0, 5), "time-invariant")
  expect_error(band(.t * 2, 0, 5), "time-invariant")
  expect_error(band(.t - .t, 0, 5), "time-invariant")
})

test_that("`.x` and `.t - .b` band the same thing as age()", {
  .i <- list(birth = datey::datey(1950))
  expected <- period_included(age(65, 95), .i = .i)
  expect_equal(period_included(band(.x, 65, 95), .i = .i), expected)
  expect_equal(period_included(band(.t - .b, 65, 95), .i = .i), expected)
})

# THE ORIGIN A DURATION IS MEASURED FROM IS AS COMPUTED AS ANYTHING ELSE. A user who has written
# `.t - .i$joined` has no reason to expect `.t - min(...)` to be refused, so the offset is any
# time-invariant datey expression rather than a bare field.
test_that("an offset may be a computed datey, not only a field", {
  .i <- list(entry = datey::datey(2018), retirement = datey::datey(2012))

  expect_equal(
    period_included(band(.t - min(.i$entry, .i$retirement), 0, 5), .i = .i),
    datey::datey_interval(datey::datey(2012), datey::datey(2017))
  )
})

test_that("a computed offset is evaluated per individual, not fixed once", {
  # The same include against two individuals whose `min` falls on different fields. If the
  # expression were resolved to one field when the include was built, one of these would be wrong.
  s <- band(.t - min(.i$entry, .i$retirement), 0, 5)

  expect_equal(
    period_included(s, .i = list(entry = datey::datey(2010), retirement = datey::datey(2018))),
    datey::datey_interval(datey::datey(2010), datey::datey(2015))
  )
  expect_equal(
    period_included(s, .i = list(entry = datey::datey(2018), retirement = datey::datey(2010))),
    datey::datey_interval(datey::datey(2010), datey::datey(2015))
  )
})

test_that("a computed offset prints as what was written", {
  shown <- function(x) paste(capture.output(print(x)), collapse = "\n")

  # `age` is the one origin with a name of its own; everything else renders as the expression, so a
  # computed offset stays readable rather than becoming a column index or a node id.
  expect_match(shown(age(65, 95)), "age [65 yr, 95 yr)", fixed = TRUE)
  expect_match(shown(band(.t - .i$joined, 0, 5)), "since .i$joined", fixed = TRUE)
  expect_match(shown(band(.t - min(.i$entry, .i$retirement), 0, 5)),
               "since min(.i$entry, .i$retirement)", fixed = TRUE)
})

test_that("an offset include needs its field, of datey type", {
  s <- band(.t - .i$joined, 0, 5)
  expect_error(period_included(s, .i = list()), "was not supplied")
  expect_error(period_included(s, .i = list(joined = 2018)), "must be a single `datey`")
})

test_that("includes are immutable", {
  p <- period(2010, 2020)
  expect_error({ p$from <- 1 }, "immutable")
})

test_that("period_included.default errors for non-includes", {
  expect_error(period_included(42, .i = list()), "not implemented")
})
