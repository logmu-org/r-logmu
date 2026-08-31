# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("indicator & indicator stays an indicator (logical AND)", {
  male <- indicator(.i$sex == "male")
  ret  <- indicator(.i$category == "ret")
  mr   <- male & ret
  expect_true(is_indicator(mr))
  expect_true(logical_value(mr, .i = list(sex = "male", category = "ret")))
  expect_false(logical_value(mr, .i = list(sex = "male", category = "dep")))
  expect_false(logical_value(mr, .i = list(sex = "female", category = "ret")))
})

test_that("indicator & include is an include gated by the indicator", {
  s <- indicator(.i$sex == "male") & age(65, 95)
  expect_true(is_include(s))
  expect_false(is_indicator(s))
  expect_equal(
    period_included(s, .i = list(sex = "male", birth = datey::datey(1950))),
    datey::datey_interval(datey::datey(2015), datey::datey(2045))
  )
  # excluded individual -> empty interval
  expect_true(is.na(period_included(s, .i = list(sex = "female", birth = datey::datey(1950)))))
})

test_that("include & include intersects the intervals", {
  s <- age(65, 95) & period(2010, 2040)
  expect_equal(
    period_included(s, .i = list(birth = datey::datey(1950))),
    datey::datey_interval(datey::datey(2015), datey::datey(2040))
  )
})

test_that("a non-overlapping intersection is empty", {
  s <- age(65, 95) & period(1990, 2000)   # [2015,2045) vs [1990,2000)
  expect_true(is.na(period_included(s, .i = list(birth = datey::datey(1950)))))
})

test_that("chained & narrows across indicators and intervals", {
  male <- indicator(.i$sex == "male")
  ret  <- indicator(.i$category == "ret")
  s <- male & ret & age(65, 95)
  expect_true(is_include(s))
  expect_equal(
    period_included(s, .i = list(sex = "male", category = "ret", birth = datey::datey(1950))),
    datey::datey_interval(datey::datey(2015), datey::datey(2045))
  )
  expect_true(is.na(period_included(s, .i = list(sex = "male", category = "dep", birth = datey::datey(1950)))))
})

test_that("a ~ formula operand is coerced to an indicator", {
  s <- age(65, 95) & (~ .i$sex == "male")
  expect_true(is_include(s))
  expect_equal(
    period_included(s, .i = list(sex = "male", birth = datey::datey(1950))),
    datey::datey_interval(datey::datey(2015), datey::datey(2045))
  )
})

test_that("a time-dependent formula cannot be intersected", {
  expect_error(age(65, 95) & (~ .t > .b), "time-invariant")
})

test_that("operators other than & are rejected at the object level", {
  male <- indicator(.i$sex == "male")
  ret  <- indicator(.i$category == "ret")
  expect_error(male | ret, "not defined")
  expect_error(period(2010, 2020) + period(2010, 2020), "not defined")
})
