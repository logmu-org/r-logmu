# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

values_of <- function(x) as.vector(unclass(x))

# Log projection factors with age, period and cohort variation, the last so
# that interpolating in cohort-period and in age-period space differ.
make_log_factor <- function(ages, periods) {
  outer(ages, periods, function(x, t) -0.015 * t + 0.0004 * x * t + 0.01 * sin(t - x))
}

u <- .Machine$double.eps / 2

test_that("projection tables hold their grid and refuse bad arguments", {
  log_factor <- make_log_factor(0:3, 0:10)
  q_projection <- q_projection_table(x_start = 60, t_start = 2015, log_factor, name = "CMI")
  log_mu_projection <- log_mu_projection_table(x_start = 60, t_start = 2015, log_factor)

  expect_s3_class(q_projection, c("q_projection_table", "projection_table"), exact = TRUE)
  expect_s3_class(log_mu_projection, c("log_mu_projection_table", "projection_table"), exact = TRUE)
  expect_identical(values_of(q_projection), as.vector(log_factor))
  expect_identical(attr(q_projection, "x_start"), datey::durationy(60))
  expect_identical(attr(q_projection, "t_start"), datey::datey(2015))
  expect_identical(name(q_projection), "CMI")
  expect_false(is_mortality(q_projection))
  expect_false(is_logmu_function(q_projection))

  expect_output(print(q_projection), "<q_projection_table: CMI> 4 ages from 60 yr, 11 years from 2015-01-01.0")
  expect_error(q_projection[1] <- 0, "immutable")
  expect_error(q_projection[[1]] <- 0, "immutable")
  expect_error(terminal_age(q_projection), "`mortality_table` or a `base_table`")

  expect_error(q_projection_table(60, 2015, c(0, 1)), "log_factor must be a matrix of finite values")
  expect_error(q_projection_table(60, 2015, matrix(c(0, NA))), "log_factor must be a matrix of finite values")
  expect_error(log_mu_projection_table(60, 2015, matrix(1L)), "log_factor must be a matrix of finite values")
  expect_error(log_mu_projection_table("60", 2015, log_factor), "x_start")
  expect_error(log_mu_projection_table(60, "2015", log_factor), "t_start")
  expect_error(log_mu_projection_table(60, 2015, log_factor, name = ""), "valid name")
})

test_that("a log mu projection on aligned grids adds the change in P", {
  log_mu_base <- c(-5, -4.9, -4.7, -4.6)
  base <- base_table(x_start = 60, as_at = 2020, log_mu = log_mu_base)
  log_factor <- make_log_factor(0:3, 0:10)
  projection <- log_mu_projection_table(x_start = 60, t_start = 2015, log_factor)

  table <- base + projection

  expect_true(is_mortality(table))
  expect_identical(dim(table), c(4L, 11L))
  expect_identical(attr(table, "x_start"), datey::durationy(60))
  expect_identical(attr(table, "t_start"), datey::datey(2015))
  expect_identical(terminal_age(table), terminal_age(base))
  expect_identical(values_of(table), as.vector(log_mu_base + (log_factor - log_factor[, 6])))
  expect_identical(projection + base, table)
})

test_that("a q projection on aligned grids multiplies q read without convexity", {
  log_mu_base <- c(-5, -4.9, -4.7, -4.6)
  base <- base_table(x_start = 60, as_at = 2020, log_mu = log_mu_base)
  log_factor <- make_log_factor(0:3, 0:10)
  projection <- q_projection_table(x_start = 60, t_start = 2015, log_factor)

  q <- -expm1(-exp(log_mu_base)) * exp(log_factor - log_factor[, 6])
  expect_identical(values_of(base + projection), as.vector(log(-log1p(-q))))
})

test_that("a projection with no change in P gives back the base table", {
  log_mu_base <- c(-5, -4.9, -4.7, -4.6)
  base <- base_table(x_start = 60, as_at = 2020, log_mu = log_mu_base)
  flat <- matrix(0, nrow = 4, ncol = 11)

  expect_identical(
    values_of(base + log_mu_projection_table(60, 2015, flat)),
    rep(log_mu_base, 11)
  )
  expect_equal(
    values_of(base + q_projection_table(60, 2015, flat)),
    rep(log_mu_base, 11),
    tolerance = 1e-14
  )
})

test_that("a table made into a projection and added to its own base period gives the table back", {
  log_mu <- -4.8 + make_log_factor(0:3, 0:10)
  table <- mortality_table(x_start = 60, t_start = 2015, log_mu = log_mu)
  base <- base_table(x_start = 60, as_at = 2020, log_mu = log_mu[, 6])

  from_log_mu <- base + as_log_mu_projection_table(table)
  from_q <- base + as_q_projection_table(table)

  for (result in list(from_log_mu, from_q)) {
    expect_identical(attr(result, "x_start"), attr(table, "x_start"))
    expect_identical(attr(result, "t_start"), attr(table, "t_start"))
    expect_equal(values_of(result), as.vector(log_mu), tolerance = 1e-14)
  }
})

test_that("as_q_projection_table() and as_log_mu_projection_table() set P from a table", {
  log_mu <- -4.8 + make_log_factor(0:3, 0:10)
  table <- mortality_table(x_start = 60.5, t_start = 2015.5, log_mu = log_mu, name = "T")

  q_projection <- as_q_projection_table(table)
  expect_s3_class(q_projection, "q_projection_table")
  expect_identical(values_of(q_projection), as.vector(log(-expm1(-exp(log_mu)))))
  expect_identical(dim(q_projection), dim(log_mu))
  expect_identical(attr(q_projection, "x_start"), datey::durationy(60.5))
  expect_identical(attr(q_projection, "t_start"), datey::datey(2015.5))
  expect_identical(name(q_projection), "T")
  expect_null(attr(q_projection, "terminal_age"))

  log_mu_projection <- as_log_mu_projection_table(table)
  expect_s3_class(log_mu_projection, "log_mu_projection_table")
  expect_identical(values_of(log_mu_projection), as.vector(log_mu))

  base <- base_table(x_start = 60, as_at = 2020, log_mu = -4)
  expect_error(as_q_projection_table(base), "must be a `mortality_table`")
  expect_error(as_log_mu_projection_table(log_mu), "must be a `mortality_table`")
  expect_error(
    as_q_projection_table(mortality_table(60, 2020, log_mu = matrix(-800))),
    "log q implied by `x` must all be finite"
  )
})

test_that("P at an as-at date between periods is interpolated in cohort-period space", {
  log_mu_base <- c(-5, -4.9, -4.7, -4.6)
  base <- base_table(x_start = 60, as_at = 2020.25, log_mu = log_mu_base)
  log_factor <- make_log_factor(55:70, 2015:2025)
  projection <- log_mu_projection_table(x_start = 55, t_start = 2015, log_factor)

  table <- base + projection
  expect_identical(attr(table, "x_start"), datey::durationy(55))
  expect_identical(attr(table, "t_start"), datey::datey(2014.25))
  expect_identical(dim(table), c(16L, 12L))

  # Column 7 is the as-at date, where nothing has changed.
  values <- unclass(table)
  expect_identical(as.vector(values[6:9, 7]), log_mu_base)

  # P a fraction f of a year after period `t`, at age `x`, both whole years
  # on the projection's grid: the cohort blend at `t` and at `t + 1`, then
  # the period blend.
  f <- 0.25
  P <- function(x, t) log_factor[x - 54, t - 2014]
  between <- function(x, t) {
    (1 - f) * ((1 - f) * P(x, t) + f * P(x - 1, t)) +
      f * ((1 - f) * P(x + 1, t + 1) + f * P(x, t + 1))
  }
  for (x in 60:63) {
    expected <- log_mu_base[x - 59] + between(x, 2021) - between(x, 2020)
    expect_equal(values[x - 54, 8], expected, tolerance = 1e-14)
  }

  # Interpolating in age-period space instead would give something else.
  in_period <- function(x, t) (1 - f) * P(x, t) + f * P(x, t + 1)
  expect_gt(abs(between(61, 2020) - in_period(61, 2020)), 1e-5)
})

test_that("the grid is the base table's, extended to cover the projection", {
  log_mu_base <- c(-5, -4.9, -4.7, -4.6)
  base <- base_table(x_start = 60, as_at = 2020.5, log_mu = log_mu_base)
  projection <- log_mu_projection_table(
    x_start = 50.5,
    t_start = 2000,
    make_log_factor(0:20, 0:40)
  )

  table <- base + projection
  expect_identical(attr(table, "x_start"), datey::durationy(50))
  expect_identical(attr(table, "t_start"), datey::datey(1999.5))
  expect_identical(dim(table), c(22L, 42L))
  expect_identical(terminal_age(table), datey::durationy(64))

  # Beyond its own ages the base table is held at its edges.
  values <- unclass(table)
  expect_identical(values[, 22], c(rep(log_mu_base[1], 10), log_mu_base, rep(log_mu_base[4], 8)))

  # A base table that already covers a projection is not extended.
  wide <- base_table(x_start = 40, as_at = 2020.5, log_mu = rep(-5, 40))
  narrow <- log_mu_projection_table(60, 2010, make_log_factor(0:5, 0:20))
  expect_identical(dim(wide + narrow), c(40L, 22L))
  expect_identical(attr(wide + narrow, "t_start"), datey::datey(2009.5))
})

test_that("a projected q is capped at 1 - u and a small q is left alone", {
  base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-40, 2.5))
  projection <- q_projection_table(60, 2020, rbind(c(0, -10), c(0, 10)))
  values <- unclass(base + projection)
  expect_identical(
    as.vector(values[, 2]),
    log(-log1p(-c(-expm1(-exp(-40)) * exp(-10), 1 - u)))
  )
  expect_lt(values[1, 2], -49)
  # With no change in P, q is read and converted back. Near q = 1 that loses
  # a few digits (mu = 12 here), so the base comes back closely, not exactly.
  expect_equal(values[, 1], c(-40, 2.5), tolerance = 1e-12)

  # Only a q that underflows to zero is refused.
  tiny <- base_table(x_start = 60, as_at = 2020, log_mu = -700)
  expect_error(
    tiny + q_projection_table(60, 2020, matrix(c(0, -100), nrow = 1)),
    "projected log mu must all be finite"
  )
})

test_that("the as-at date must lie within the projection's periods", {
  base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9))
  log_factor <- make_log_factor(0:1, 0:10)
  expect_error(
    base + log_mu_projection_table(60, 2021, log_factor),
    "as-at date \\(2020-01-01.0\\) must lie within the projection's periods \\(2021-01-01.0 to 2031-01-01.0\\)"
  )
  expect_error(base + log_mu_projection_table(60, 2009, log_factor), "must lie within")
  expect_no_error(base + log_mu_projection_table(60, 2010, log_factor))
  expect_no_error(base + log_mu_projection_table(60, 2020, log_factor))
})

test_that("the result is named after both parts when both are named", {
  log_factor <- make_log_factor(0:1, 0:10)
  base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9))
  named_base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9), name = "S4PFA")
  projection <- q_projection_table(60, 2015, log_factor)
  named_projection <- q_projection_table(60, 2015, log_factor, name = "CMI_2025(0.015)")

  expect_identical(name(named_base + named_projection), "S4PFA + CMI_2025(0.015)")
  expect_null(name(base + named_projection))
  expect_null(name(named_base + projection))
})

test_that("the result is an ordinary mortality_table", {
  base <- base_table(x_start = 60.5, as_at = 2020.5, q = c(0.01, 0.012, 0.0145))
  projection <- q_projection_table(60.5, 2015.5, make_log_factor(0:2, 0:10))
  table <- base + projection
  .i <- list(birth = datey::datey(1950))
  expect_identical(terminal(table, .i), datey::datey(1950) + terminal_age(base))
  expect_length(log_mu(table, .i, datey::datey(c(2015, 2020, 2030))), 3L)
})

test_that("a base table and a projection can only be added to each other", {
  log_factor <- make_log_factor(0:1, 0:10)
  base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9))
  projection <- log_mu_projection_table(60, 2015, log_factor)
  table <- mortality_table(base)
  message <- "can only be added to each other"

  expect_error(base + 1, message)
  expect_error(1 + projection, message)
  expect_error(base - projection, message)
  expect_error(base * projection, message)
  expect_error(projection + projection, message)
  expect_error(base + base, message)
  expect_error(-base, message)

  # Against an operand with an `Ops` method of its own, ours is chosen, so
  # there is neither R's `Incompatible methods` warning nor arithmetic on the
  # numbers underneath.
  for (pair in list(list(base, table), list(table, projection), list(base, datey::datey(2020)))) {
    warned <- FALSE
    expect_error(
      withCallingHandlers(
        pair[[1]] + pair[[2]],
        warning = function(w) {
          warned <<- TRUE
          invokeRestart("muffleWarning")
        }
      ),
      message
    )
    expect_false(warned)
  }
})
