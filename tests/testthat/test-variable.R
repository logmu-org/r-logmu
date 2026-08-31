# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("variable() narrows to the most specific provable type", {
  # uses .t -> general variable
  gv <- variable(.i$pension * (.t - .b))
  expect_true(is_variable(gv))
  expect_false(is_static_variable(gv))
  expect_false(is_indicator(gv))

  # time-invariant, non-logical -> static_variable
  sv <- variable(.i$pension * 2)
  expect_true(is_static_variable(sv))
  expect_false(is_indicator(sv))
  expect_true(is_logmu_function(sv))

  # time-invariant, logical -> indicator (hence also an include)
  iv <- variable(.i$pension > 0)
  expect_true(is_indicator(iv))
  expect_true(is_include(iv))

  # static_variable() asserts no .t but still narrows to indicator when logical
  expect_true(is_static_variable(static_variable(.i$pension)))
  expect_true(is_indicator(static_variable(.i$sex == "male")))
})

test_that("static_variable() rejects time dependence", {
  expect_error(static_variable(.t - .b), "must not depend on time")
  expect_error(static_variable(.x), "must not depend on time")
})

test_that("value() evaluates against facts and time", {
  expect_equal(
    value(static_variable(.i$pension), .i = list(pension = 1000)),
    1000
  )
  expect_equal(
    value(variable(.i$pension * 2), .i = list(pension = 1000)),
    2000
  )
  # a comparison evaluates to logical
  expect_equal(
    value(variable(.i$pension > 0), .i = list(pension = 1000)),
    TRUE
  )
})

test_that("value() needs .t when the variable uses time", {
  v <- variable(.t)
  expect_error(value(v, .i = list()), "uses `.t`")
})

test_that("value() errors on a missing field", {
  expect_error(
    value(static_variable(.i$missing), .i = list(pension = 1)),
    "was not supplied"
  )
})

test_that("a variable is immutable", {
  v <- variable(.i$pension)
  expect_error({ v$ast <- 1 }, "immutable")
  expect_error({ v[[1]] <- 1 }, "immutable")
})

test_that("value.default errors for non-variables", {
  expect_error(value(42, .i = list()), "not implemented")
})
