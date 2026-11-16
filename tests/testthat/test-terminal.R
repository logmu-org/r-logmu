# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

rates <- matrix(c(0.01, 0.02, 0.03, 0.011, 0.021, 0.031), nrow = 3L, ncol = 2L)
births <- list(birth = datey::datey(c(1950, 1961.5)))
none <- datey::datey(datey::valid_years_end)

test_that("a table's terminal age is the row after the last, plus half a year for q", {
  expect_identical(
    terminal_age(mortality_table(x_start = 60, t_start = 2020, mu = rates)),
    datey::durationy(63)
  )
  expect_identical(
    terminal_age(mortality_table(x_start = 60, t_start = 2020, log_mu = log(rates))),
    datey::durationy(63)
  )
  expect_identical(
    terminal_age(mortality_table(x_start = 60.5, t_start = 2020.5, q = rates)),
    datey::durationy(64)
  )
})

test_that("terminal_age() refuses anything but a table and points to terminal()", {
  pointer <- "use `terminal\\(x, \\.i\\)`"
  expect_error(terminal_age(mortality_const(log_mu = -4)), pointer)
  expect_error(terminal_age(gompertz_mortality()), pointer)
  base <- mortality_const(log_mu = -4)
  expect_error(terminal_age(mortality(base + 0.1)), pointer)
  expect_error(terminal_age(42), pointer)
})

test_that("a table's terminal date is each birth plus its terminal age", {
  table <- mortality_table(x_start = 60, t_start = 2020, mu = rates)
  result <- terminal(table, births)
  expect_true(datey::is_datey(result))
  expect_identical(result, datey::datey(c(2013, 2024.5)))
})

test_that("a mortality with no terminal age gives the last valid date, once per individual", {
  expect_identical(terminal(mortality_const(log_mu = -4), births), rep(none, 2L))
  expect_identical(terminal(gompertz_mortality(), births), rep(none, 2L))
  base <- mortality_const(log_mu = -4)
  expect_identical(terminal(mortality(base + 0.05), births), rep(none, 2L))
})

test_that("an expression that refers to a table is refused for now", {
  table <- mortality_table(x_start = 60, t_start = 2020, mu = rates)
  base <- mortality_const(log_mu = -4)
  expect_error(terminal(mortality(table + 0.1), births), "does not yet support")
  expect_error(terminal(mortality(ifelse(.i$male, table, base)), births), "does not yet support")
})

test_that("terminal() needs a mortality and a valid .i", {
  expect_error(terminal(42, births), "not implemented")
  expect_error(terminal(mortality_const(log_mu = -4), list(age = 1)), "birth")
})
