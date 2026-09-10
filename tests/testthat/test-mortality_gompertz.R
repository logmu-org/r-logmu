# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `gompertz_mortality()` is stored in its COHORT form but specified in the
# age-period form, so the load-bearing question is whether the rewrite is
# faithful. Two independent oracles answer it:
#
#   1. the same law written out by hand as an age-period pronoun expression,
#      which must agree in R and through the engine; and
#   2. the closed-form integral of `exp(c + s(t - t_0))`, which the engine's
#      numerical integration must converge to as `time_scale` falls.
#
# The cohort origin `b_0 = t_0 - x_0` is what makes the rewrite hold, and both
# oracles fail if it is wrong.

log_mu_0 <- -4.6
slope_x <- 0.10
x_0 <- 70
slope_t <- -0.02
t_0 <- 2020

gompertz <- gompertz_mortality(
  log_mu_0 = log_mu_0, slope_x = slope_x, x_0 = x_0,
  slope_t = slope_t, t_0 = t_0
)

# The age-period form, written out independently of the constructor.
age_period_log_mu <- function(birth, t) {
  x <- as.numeric(t) - as.numeric(birth)
  log_mu_0 + slope_x * (x - x_0) + slope_t * (as.numeric(t) - t_0)
}

test_that("gompertz_mortality() builds a mortality", {
  expect_true(is_mortality(gompertz))
  expect_s3_class(gompertz, "mortality_expr")
})

test_that("the defaults are the recorded standard basis", {
  # log mu(x, t) = -3.8 + 0.1 (x - 75) - 0.01 (t - 2020)
  #
  # Read at four points rather than compared against a call that repeats the
  # same five numbers, so each default is pinned by an independent arithmetic
  # statement. Between them the four fix all five: any change to `x_0` or `t_0`
  # moves the first, and the last two separate the slopes from each other.
  standard <- gompertz_mortality()

  # Age 75 in 2020, so sitting on both origins.
  expect_equal(
    log_mu(standard, list(birth = datey::datey(1945)), datey::datey(2020)),
    -3.8
  )
  # Age 65 in 2020: age moves alone.
  expect_equal(
    log_mu(standard, list(birth = datey::datey(1955)), datey::datey(2020)),
    -3.8 - 0.1 * 10
  )
  # Age 75 in 2010: time moves alone.
  expect_equal(
    log_mu(standard, list(birth = datey::datey(1935)), datey::datey(2010)),
    -3.8 + 0.01 * 10
  )
  # Age 85 in 2030: ten years of ageing against ten of improvement.
  expect_equal(
    log_mu(standard, list(birth = datey::datey(1945)), datey::datey(2030)),
    -3.8 + 0.1 * 10 - 0.01 * 10
  )
})

test_that("log mu is log_mu_0 at age x_0 and time t_0", {
  # Born t_0 - x_0, so exactly age x_0 at t_0.
  expect_equal(
    log_mu(gompertz, list(birth = datey::datey(t_0 - x_0)), datey::datey(t_0)),
    log_mu_0
  )
})

test_that("the cohort form reproduces the age-period form", {
  # A grid that moves age and time independently, and off whole years, so a
  # confusion of the two origins cannot pass by coincidence.
  times <- datey::datey(c(2000, 2010.5, 2020, 2031.25))
  for (birth in c(1930, 1945, 1950.25, 1962)) {
    expect_equal(
      log_mu(gompertz, list(birth = datey::datey(birth)), times),
      age_period_log_mu(birth, as.numeric(times))
    )
  }
})

test_that("each slope moves log mu in the right direction", {
  b <- list(birth = datey::datey(1950))

  # Age: a year older is `slope_x` higher, holding the cohort fixed means
  # moving forward a year, which also picks up `slope_t`.
  older <- log_mu(gompertz, list(birth = datey::datey(1949)), datey::datey(2020))
  expect_equal(older - log_mu(gompertz, b, datey::datey(2020)), slope_x)

  # Time at a fixed age: born a year later and observed a year later.
  later <- log_mu(gompertz, list(birth = datey::datey(1951)), datey::datey(2021))
  expect_equal(later - log_mu(gompertz, b, datey::datey(2020)), slope_t)
})

test_that("slopes that cancel leave an individual's log mu flat in time", {
  # slope_x + slope_t == 0, so ageing and improvement exactly offset.
  flat <- gompertz_mortality(log_mu_0 = -4, slope_x = 0.1, x_0 = 70,
                             slope_t = -0.1, t_0 = 2020)
  got <- log_mu(flat, list(birth = datey::datey(1950)), datey::datey(c(2000, 2020, 2040)))
  expect_equal(got, rep(-4, 3L))
})

test_that("x_0 and t_0 accept datey/durationy as well as numbers", {
  as_objects <- gompertz_mortality(
    log_mu_0 = log_mu_0, slope_x = slope_x, x_0 = datey::durationy(x_0),
    slope_t = slope_t, t_0 = datey::datey(t_0)
  )
  expect_identical(as_objects$ast, gompertz$ast)
})

test_that("a gompertz splices into a pronoun expression", {
  scaled <- mortality(gompertz + 0.05)
  expect_equal(
    log_mu(scaled, list(birth = datey::datey(1950)), datey::datey(2020)),
    log_mu_0 + 0.05
  )
})

test_that("an optional name is carried", {
  named <- gompertz_mortality(log_mu_0 = -4, slope_x = 0.1, x_0 = 70,
                              slope_t = -0.02, t_0 = 2020, name = "Test basis")
  # `exact = TRUE` matters: a `mortality_expr` is a list, so a partial match
  # would find its `names` attribute instead.
  expect_identical(attr(named, "name", exact = TRUE), "Test basis")
  expect_null(attr(gompertz, "name", exact = TRUE))
})

test_that("gompertz_mortality() rejects parameters it cannot use", {
  ok <- list(log_mu_0 = -4, slope_x = 0.1, x_0 = 70, slope_t = -0.02, t_0 = 2020)
  bad <- function(...) do.call(gompertz_mortality, utils::modifyList(ok, list(...)))

  # Match the whole message, not just the parameter name: the name alone also
  # appears in the parser's own error, so it would pass with no check at all.
  expect_error(bad(log_mu_0 = NA_real_), "log_mu_0 must be a finite numeric scalar")
  expect_error(bad(log_mu_0 = c(-4, -5)), "log_mu_0 must be a finite numeric scalar")
  expect_error(bad(slope_x = Inf), "slope_x must be a finite numeric scalar")
  expect_error(bad(slope_t = NaN), "slope_t must be a finite numeric scalar")
  expect_error(bad(x_0 = "seventy"), "`x_0` must be a valid `durationy`")
  expect_error(bad(t_0 = NA), "`t_0` must be a valid `datey`")
  expect_error(bad(name = 42), "`name` must be a valid name")

  # Two finite slopes can still overflow when added.
  expect_error(bad(slope_x = 1e308, slope_t = 1e308),
               "slope_x + slope_t must be finite", fixed = TRUE)
})

# ---- through the engine ----------------------------------------------------

gompertz_data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
    E2R_start = datey::datey(rep(2015, 5)),
    E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
    E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

# The closed form of E for this law, integrated exactly rather than on a grid.
# For one individual log mu is `c + s(t - t_0)` with `c` fixed by birth, so
# `E = exp(c)/s * [exp(s(t2 - t_0)) - exp(s(t1 - t_0))]`.
gompertz_E_exact <- function() {
  birth <- c(1945, 1950, 1955, 1940, 1948)
  start <- rep(2015, 5)
  end <- c(2020, 2020, 2018, 2020, 2019)
  s <- slope_x + slope_t
  c_i <- log_mu_0 - slope_x * (birth - (t_0 - x_0))
  sum(exp(c_i) / s * (exp(s * (end - t_0)) - exp(s * (start - t_0))))
}

test_that("the engine agrees with the same law written age-period", {
  age_period <- mortality(-4.6 + 0.10 * (.x - 70) - 0.02 * (.t - 2020))
  basis <- settings(overdispersion = 1, time_scale = 1 / 12)

  from_gompertz <- aev(gompertz_data, mortality = gompertz, settings = basis)
  from_expression <- aev(gompertz_data, mortality = age_period, settings = basis)

  expect_equal(from_gompertz$A, from_expression$A)
  expect_equal(from_gompertz$E, from_expression$E, tolerance = 1e-12)
  expect_equal(from_gompertz$V, from_expression$V, tolerance = 1e-12)
})

test_that("the engine converges to the closed-form integral", {
  exact <- gompertz_E_exact()

  coarse <- aev(gompertz_data, mortality = gompertz,
                settings = settings(overdispersion = 1, time_scale = 1 / 4))
  fine <- aev(gompertz_data, mortality = gompertz,
              settings = settings(overdispersion = 1, time_scale = 1 / 60))

  expect_equal(fine$E, exact, tolerance = 1e-5)

  # The residual is quadrature error, so a finer grid must be strictly closer.
  expect_lt(abs(fine$E - exact), abs(coarse$E - exact))
})
