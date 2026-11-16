# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

values_of <- function(x) as.vector(unclass(x))

test_that("base_table() holds log mu from log_mu, mu or q", {
  log_mu <- c(-5, -4.9, -4.7, -4.6)
  from_log_mu <- base_table(x_start = 60, as_at = 2020, log_mu = log_mu)
  expect_identical(values_of(from_log_mu), log_mu)
  expect_identical(attr(from_log_mu, "x_start"), datey::durationy(60))
  expect_identical(as_at(from_log_mu), datey::datey(2020))
  expect_s3_class(from_log_mu, "base_table")
  expect_false(is_mortality(from_log_mu))
  expect_false(is_logmu_function(from_log_mu))

  from_mu <- base_table(x_start = 60, as_at = 2020, mu = exp(log_mu))
  expect_equal(values_of(from_mu), log_mu)
})

test_that("base_table() converts q along age with convexity, edges unadjusted", {
  q <- c(0.010, 0.012, 0.0145, 0.017)
  m <- -log1p(-q)
  mu <- c(
    m[1],
    m[2] - (m[3] - 2 * m[2] + m[1]) / 24,
    m[3] - (m[4] - 2 * m[3] + m[2]) / 24,
    m[4]
  )
  base <- base_table(x_start = 60.5, as_at = 2020.5, q = q)
  expect_equal(values_of(base), log(mu), tolerance = 1e-14)

  # The same conversion as a one-period `mortality_table()`.
  table <- mortality_table(x_start = 60.5, t_start = 2020.5, q = matrix(q))
  expect_identical(values_of(base), values_of(table))
})

test_that("base_table() keeps no q", {
  base <- base_table(x_start = 60, as_at = 2020, q = c(0.01, 0.012, 0.0145))
  expect_identical(
    sort(names(attributes(base))),
    c("as_at", "class", "terminal_age", "x_start")
  )
})

test_that("base_table() terminal age follows mortality_table()", {
  rates <- c(0.01, 0.012, 0.0145)
  expect_identical(
    terminal_age(base_table(x_start = 60, as_at = 2020, log_mu = log(rates))),
    datey::durationy(63)
  )
  expect_identical(
    terminal_age(base_table(x_start = 60, as_at = 2020, mu = rates)),
    datey::durationy(63)
  )
  expect_identical(
    terminal_age(base_table(x_start = 60.5, as_at = 2020.5, q = rates)),
    datey::durationy(64)
  )
})

test_that("base_table() refuses bad arguments", {
  expect_error(base_table(60, 2020), "At least one of q, mu or log_mu")
  expect_error(base_table(60, 2020, q = 0.01, log_mu = -4), "Only one of")
  expect_error(base_table(60, 2020, mu = 0.01, log_mu = -4), "Only one of")
  expect_error(base_table(60, 2020, log_mu = matrix(-4)), "log_mu must be a vector of doubles")
  expect_error(base_table(60, 2020, q = 1L), "q must be a vector of doubles")
  expect_error(base_table(60, 2020, log_mu = c(-4, Inf)), "log_mu must be a vector of finite values")
  expect_error(base_table(60, 2020, mu = c(0.01, 0)), "log\\(mu\\) must all be finite")
  expect_error(base_table(60, 2020, q = c(0.01, 1)), "implied by q must all be finite")
  expect_error(base_table("60", 2020, log_mu = -4), "x_start")
  expect_error(base_table(60, "2020", log_mu = -4), "as_at")
  expect_error(base_table(60, 2020, log_mu = -4, name = " x"), "valid name")
})

test_that("a base_table is immutable and prints its grid", {
  base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9), name = "S4")
  expect_error(base[1] <- 0, "immutable")
  expect_error(base[[1]] <- 0, "immutable")
  expect_output(print(base), "<base_table: S4> 2 ages from 60 yr, as at 2020-01-01.0")
})

test_that("mortality_table() turns a base_table into a mortality that is constant in time", {
  log_mu <- c(-5, -4.9, -4.7)
  base <- base_table(x_start = 60, as_at = 2020, log_mu = log_mu, name = "S4")
  table <- mortality_table(base)

  expect_true(is_mortality(table))
  expect_identical(
    table,
    mortality_table(x_start = 60, t_start = 2020, log_mu = matrix(log_mu), name = "S4")
  )
  expect_identical(terminal_age(table), terminal_age(base))

  .i <- list(birth = datey::datey(1959))
  expect_identical(log_mu(table, .i, datey::datey(2020)), log_mu[2])
  expect_identical(log_mu(table, list(birth = datey::datey(1909)), datey::datey(1970)), log_mu[2])
  expect_identical(log_mu(table, list(birth = datey::datey(2009)), datey::datey(2070)), log_mu[2])

  expect_identical(name(mortality_table(base, name = "other")), "other")
  expect_error(mortality_table(base, 2020), "either a `base_table` alone")
  expect_error(mortality_table(base, log_mu = matrix(-4)), "either a `base_table` alone")
  expect_error(mortality_table(base, name = ""), "valid name")
})

test_that("mu_scaled() multiplies mu", {
  log_mu <- c(-5, -4.9, -4.7)
  base <- base_table(x_start = 60, as_at = 2020, log_mu = log_mu, name = "S4PFA")
  scaled <- mu_scaled(base, 1.1)
  expect_identical(values_of(scaled), log_mu + log(1.1))
  expect_identical(attr(scaled, "x_start"), attr(base, "x_start"))
  expect_identical(as_at(scaled), as_at(base))
  expect_identical(terminal_age(scaled), terminal_age(base))
  expect_identical(name(scaled), "110%\u03bc S4PFA")
})

test_that("q_scaled() multiplies q read without convexity", {
  log_mu <- c(-5, -4.9, -2, 0.5)
  base <- base_table(x_start = 60.5, as_at = 2020.5, log_mu = log_mu, name = "S4PFA")
  scaled <- q_scaled(base, 1.125)
  expect_identical(
    values_of(scaled),
    log(-log1p(-(1.125 * -expm1(-exp(log_mu)))))
  )
  expect_identical(terminal_age(scaled), terminal_age(base))
  expect_identical(as_at(scaled), as_at(base))
  expect_identical(name(scaled), "112.5%q S4PFA")

  # A factor of one gives the table back, there being no convexity step.
  expect_equal(values_of(q_scaled(base, 1)), log_mu, tolerance = 1e-14)
})

test_that("q_scaled() caps q at 1 - u and leaves a small q alone", {
  u <- .Machine$double.eps / 2
  base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-40, 3))
  scaled <- q_scaled(base, 1.1)
  expect_identical(
    values_of(scaled),
    log(-log1p(-c(1.1 * -expm1(-exp(-40)), 1 - u)))
  )
  expect_true(all(is.finite(values_of(scaled))))
  expect_equal(values_of(q_scaled(base, 1))[1], -40, tolerance = 1e-14)

  # Only a q that underflows to zero is refused.
  tiny <- base_table(x_start = 60, as_at = 2020, log_mu = -700)
  expect_error(q_scaled(tiny, 1e-30), "scaled log mu must all be finite")
})

test_that("scaling names and refusals", {
  base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9))
  expect_null(name(q_scaled(base, 1.1)))
  expect_null(name(mu_scaled(base, 1.1)))
  named <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9), name = "S4PFA")
  expect_identical(name(q_scaled(named, 0.95)), "95%q S4PFA")
  expect_identical(name(q_scaled(q_scaled(named, 1.1), 0.9)), "90%q 110%q S4PFA")

  expect_error(q_scaled(base, 0), "single positive number")
  expect_error(mu_scaled(base, -1), "single positive number")
  expect_error(q_scaled(base, c(1, 2)), "single positive number")
  expect_error(q_scaled(base, NA_real_), "single positive number")
  table <- mortality_table(base)
  expect_error(q_scaled(table, 1.1), "must be a `base_table`")
})

test_that("terminal_age() accepts a base_table and names it in its refusal", {
  expect_error(terminal_age(mortality_const(log_mu = -4)), "`mortality_table` or a `base_table`")
})

test_that("as_at() reads a base_table or a val_data and nothing else", {
  base <- base_table(x_start = 60, as_at = 2020.5, log_mu = c(-5, -4.9))
  expect_identical(as_at(base), datey::datey(2020.5))

  valuation <- val_data(list(birth = datey::datey(c(1950, 1960))), datey::datey(2025))
  expect_identical(as_at(valuation), datey::datey(2025))

  expect_error(as_at(mortality_table(base)), "`val_data` or a `base_table`")
  expect_error(as_at(42), "`val_data` or a `base_table`")
})
