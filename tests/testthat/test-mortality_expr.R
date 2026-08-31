# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("mortality() scales a mortality (additive in log mu)", {
  base <- mortality_const(log_mu = -4)
  .i <- list(birth = datey::datey(1950))
  .t <- datey::datey(2020)

  m <- mortality(base + 0.05)
  expect_true(is_mortality(m))
  expect_s3_class(m, "mortality_expr")
  expect_equal(log_mu(m, .i, .t), -3.95)
})

test_that("a bare mortality reference is returned unchanged", {
  base <- mortality_const(log_mu = -4)
  expect_identical(mortality(base), base)
})

test_that("mortality() selects between mortalities", {
  hi <- mortality_const(log_mu = -3)
  lo <- mortality_const(log_mu = -5)
  m <- mortality(ifelse(.i$smoker, hi, lo))
  b <- datey::datey(1950)
  t <- datey::datey(2020)
  expect_equal(log_mu(m, list(birth = b, smoker = TRUE),  t), -3)
  expect_equal(log_mu(m, list(birth = b, smoker = FALSE), t), -5)
})

test_that("an inline mortality constructor is read as a leaf", {
  m <- mortality(ifelse(.i$smoker, mortality_const(log_mu = -3), mortality_const(log_mu = -5)))
  b <- datey::datey(1950)
  t <- datey::datey(2020)
  expect_equal(log_mu(m, list(birth = b, smoker = TRUE), t), -3)
})

test_that("a closed-form log-mu formula needs no leaf", {
  m <- mortality(.i$log_mu_value)
  expect_equal(log_mu(m, list(log_mu_value = -4), datey::datey(2020)), -4)
})

test_that("mortality() rejects non-mortality leaves", {
  inc <- age(65, 95)
  expect_error(
    mortality(ifelse(.i$smoker, inc, mortality_const(log_mu = -4))),
    "only reference"
  )
})
