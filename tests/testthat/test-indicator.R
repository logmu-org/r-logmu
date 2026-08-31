# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("an indicator is both a variable and an include", {
  male <- indicator(.i$sex == "male")
  expect_true(is_indicator(male))
  expect_true(is_variable(male))
  expect_true(is_static_variable(male))
  expect_true(is_include(male))
  expect_true(is_logmu_function(male))
})

test_that("an indicator must not depend on time", {
  expect_error(indicator(.t > .b), "must not depend on time")
})

test_that("logical_value evaluates to TRUE/FALSE", {
  male <- indicator(.i$sex == "male")
  expect_true(logical_value(male, .i = list(sex = "male")))
  expect_false(logical_value(male, .i = list(sex = "female")))
})

test_that("indicators may be built from logical, integer or double 0/1", {
  flag <- indicator(.i$flag)
  expect_true(logical_value(flag, .i = list(flag = 1L)))
  expect_false(logical_value(flag, .i = list(flag = 0)))
  expect_true(logical_value(flag, .i = list(flag = TRUE)))
  # a non-0/1 value is rejected at evaluation
  expect_error(logical_value(flag, .i = list(flag = 2)), "0 or 1")
})

test_that("a constant non-0/1 indicator is rejected at construction", {
  expect_error(indicator(2), "0 or 1")
  expect_silent(indicator(1))
  expect_silent(indicator(TRUE))
})

test_that("an indicator coerces to an include interval", {
  male <- indicator(.i$sex == "male")
  expect_equal(
    period_included(male, .i = list(sex = "male")),
    datey::datey_interval(datey::datey(datey::valid_years_start),
                          datey::datey(datey::valid_years_end))
  )
  expect_true(is.na(period_included(male, .i = list(sex = "female"))))
})

test_that("value() (the variable side) returns the raw evaluation", {
  male <- indicator(.i$sex == "male")
  expect_true(value(male, .i = list(sex = "male")))
})

test_that("an indicator is immutable", {
  male <- indicator(.i$sex == "male")
  expect_error({ male$ast <- 1 }, "immutable")
})
