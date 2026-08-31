# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("simple test of mortality_const", {

  q <- 0.01
  mu <- -log1p(-q)
  log_mu <- log(mu)
  const_q <- mortality_const(q = q)
  const_mu <- mortality_const(mu = mu)
  const_log_mu <- mortality_const(log_mu = log_mu)

  expect_equal(const_q, const_log_mu)
  expect_equal(const_mu, const_log_mu)

  test_on_const <- function(mortality, .b, .t) {
    .b <- datey::datey(.b)
    .t <- datey::datey(.t)
    .i <- list(birth = .b)
    actual <- log_mu(mortality, .i, .t)
    expect_equal(actual, rep_len(log_mu, length(.t)))
  }
  test <- function(.b, .t) {
    test_on_const(const_q, .b, .t)
    test_on_const(const_mu, .b, .t)
    test_on_const(const_log_mu, .b, .t)
  }

  test(1950, 2020)
  test(1900, c(2000, 2020))
})
