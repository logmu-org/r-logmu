# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the canonical conversions between clicks and years.
#
# Clicks to years is a multiply by the reciprocal rather than a divide, because a divide is not
# rewritten to a multiply without -ffast-math and is several times slower in vectorised code. What
# justifies that choice is that it agrees with a correctly rounded divide at every granularity an
# expression is written in -- whole years, ages, months, quarters -- and only differs, by one unit in
# the last place, below a month. That agreement follows from the particular value of ClicksPerYear
# rather than from the arithmetic, so these tests are what notices if that value ever changes.
#
# Comparison narrowing rewrites a comparison into click space only when a number lands on a whole
# click, and decides that by converting and converting back, so the round-trip here is the property
# that pass rests on.

cpy <- 534360

test_that("clicks per year has the divisors the framework needs", {
  expect_equal(cpy %% 4, 0)   # quarters
  expect_equal(cpy %% 12, 0)  # months
  expect_equal(cpy %% 60, 0)  # the finest integration step
  expect_equal(cpy %% 365, 0) # a day in a common year
  expect_equal(cpy %% 366, 0) # a day in a leap year
})

test_that("whole years across the representable calendar convert exactly", {
  years <- 1000:3000
  expect_identical(cpp_veil_to_years(as.integer(years * cpy)), as.double(years))
})

test_that("whole ages convert exactly, so age last birthday is right on the birthday", {
  ages <- 0:120
  converted <- cpp_veil_to_years(as.integer(ages * cpy))
  expect_identical(converted, as.double(ages))
  expect_identical(floor(converted), as.double(ages)) # `floor(.x)` is age last birthday
})

# Below a whole year the true value is generally not representable, so the property to hold the
# conversion to is that it agrees with a correctly rounded divide -- not that it equals whatever
# expression R would write, which rounds twice and can land elsewhere.

test_that("month boundaries agree with a correctly rounded divide", {
  # Built by integer arithmetic: (a + m/12) * cpy would be inexact, and as.integer truncates.
  clicks <- as.integer(outer(0:120, 0:11, function(a, m) a * cpy + m * (cpy / 12)))
  expect_identical(cpp_veil_to_years(clicks), clicks / cpy)
})

test_that("quarter boundaries agree with a correctly rounded divide", {
  clicks <- as.integer(outer(0:120, 0:3, function(a, q) a * cpy + q * (cpy / 4)))
  expect_identical(cpp_veil_to_years(clicks), clicks / cpy)
})

test_that("the two conversions round-trip every valid click", {
  set.seed(4)
  lo <- 1000 * cpy
  hi <- 3000 * cpy
  clicks <- as.integer(c(lo, lo + 1, hi - 1, hi, lo + sample.int(hi - lo, 200000)))
  expect_identical(cpp_veil_to_clicks(cpp_veil_to_years(clicks)), clicks)
})

test_that("years to clicks rounds half to even, as the datey spec requires", {
  # Exactly halfway between clicks: banker's rounding takes the even one each time, on both sides of
  # zero. The rounding is spelled out in C++ rather than delegated to nearbyint, which would follow
  # whatever rounding mode the process happens to be in.
  expect_identical(cpp_veil_to_clicks(c(0.5, 1.5, 2.5, 3.5) / cpy), c(0L, 2L, 2L, 4L))
  expect_identical(cpp_veil_to_clicks(c(-0.5, -1.5, -2.5, -3.5) / cpy), c(0L, -2L, -2L, -4L))
})

test_that("values either side of a tie round to the nearer click", {
  expect_identical(cpp_veil_to_clicks(c(2.4, 2.6, -2.4, -2.6) / cpy), c(2L, 3L, -2L, -3L))
})

test_that("a value that cannot be represented as clicks reports NA rather than overflowing", {
  expect_identical(cpp_veil_to_clicks(c(1e9, -1e9, Inf, -Inf, NaN, NA_real_)), rep(NA_integer_, 6))
})
