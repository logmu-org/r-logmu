# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests `vector_log_mu`: an age-period table read over the time vector.
#
# THERE ARE TWO ORACLES HERE, AND NEITHER IS A SECOND COPY OF THE LOOKUP.
#
#   1. `log_mu()` on the same table. The package already answers "what is log mu for this individual
#      at this instant", and veil has to agree with it exactly -- the same C++ function serves both,
#      so this checks the plumbing rather than the arithmetic: that the right table reaches the right
#      instruction, that birth stays in clicks, and that the sample points are what the grid says.
#      It is applied through `died_value`, whose slot is the exposure end -- a click R knows exactly,
#      with no need to work out where any other sample point fell.
#
#   2. The analytic integral. The lookup interpolates bilinearly on the cohort-period lattice, so a
#      table whose log mu is affine in the age and period indices produces a log mu that is affine in
#      time along any one cohort, with no knot at the year boundaries. The midpoint rule is exact for
#      a linear integrand, so the expected integral can be written down.
#
# EXPOSURES ARE WHOLE QUARTERS on purpose, so there is no short final interval and therefore none of
# the half-click midpoint rounding that test-veil_time.R has to allow for. What is being tested here
# is the lookup, not the grid.

clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

# The table's own extent. Wide enough that none of the individuals below is ever clamped, which
# matters because clamping is silent -- it is asserted deliberately in its own test instead.
x0_years <- 55
t0_years <- 2005
age_count <- 40L
period_count <- 20L

birth_years <- c(1940, 1945, 1950)
start_clicks <- c(2010, 2011, 2010) * clicks_per_year
end_clicks <- start_clicks + c(12L, 8L, 4L) * quarter

# Everyone dies, so the value at death is available for all three rather than just one.
cols <- list(
  birth     = datey::datey(birth_years),
  E2R_start = datey_clicks(start_clicks),
  E2R_end   = datey_clicks(end_clicks),
  E2R_died  = c(TRUE, TRUE, TRUE)
)

exposure_years <- (end_clicks - start_clicks) / clicks_per_year

constant_table <- function(value) {
  mortality_table(
    x0 = x0_years, t0 = t0_years,
    log_mu = matrix(value, nrow = age_count, ncol = period_count)
  )
}

# log mu = base + age_slope * age_index + period_slope * period_index, on the tabulated lattice.
# Bilinear interpolation reproduces an affine function exactly, so this is also log mu everywhere
# between the lattice points.
affine_table <- function(base, age_slope, period_slope) {
  age_index <- seq_len(age_count) - 1L
  period_index <- seq_len(period_count) - 1L
  log_mu <- outer(age_index, period_index, function(x, t) base + age_slope * x + period_slope * t)
  mortality_table(x0 = x0_years, t0 = t0_years, log_mu = log_mu)
}

# Where each individual's cohort sits on the lattice, in years from the table's own origin. The
# earliest cohort the table covers is the one aged x0 at t0.
cohort_index <- birth_years - (t0_years - x0_years)

# log mu along one cohort, as a function of calendar time in years. Substituting age = t - b into the
# affine table leaves a straight line in t whose slope is the sum of the two tabulated slopes.
affine_log_mu <- function(base, age_slope, period_slope, i, t_years) {
  base - age_slope * cohort_index[i] + (age_slope + period_slope) * (t_years - t0_years)
}

# The integral of a linear integrand, with the two squares kept apart rather than subtracted: over
# calendar time they are around four million and their difference a few thousand, so the obvious
# spelling would lose most of its precision to cancellation.
analytic_linear <- function(from, to) (to - from) * (to + from) / 2

integrate_obj <- function(obj, time_scale = quarter_scale) {
  cpp_veil_integrate(it_obj(obj), cols, time_scale, NULL)
}

test_that("a mortality table lowers to vector_log_mu over the time vector", {
  res <- integrate_obj(constant_table(-3.5))

  expect_true("vector_log_mu" %in% res$monikers)
  expect_true("integrate" %in% res$monikers)

  # Nothing converts the sample points to years and nothing broadcasts: the table is read at the
  # slots directly, so the body is the source op and the reduction and no more.
  expect_false("vector_t" %in% res$monikers)
  expect_false("broadcast" %in% res$monikers)
})

test_that("a constant table integrates to log mu times the exposure", {
  value <- -3.5
  res <- integrate_obj(constant_table(value))

  expect_equal(res$integral, value * exposure_years, tolerance = 1e-12)
  expect_equal(res$died_value, rep(value, 3), tolerance = 1e-12)
})

test_that("the value at death agrees with log_mu() on the same table", {
  tbl <- affine_table(base = -4.2, age_slope = 0.09, period_slope = -0.015)
  res <- integrate_obj(tbl)

  expected <- vapply(
    seq_along(birth_years),
    function(i) log_mu(tbl, list(birth = cols$birth[i]), datey_clicks(end_clicks[i])),
    double(1)
  )

  expect_equal(res$died_value, expected, tolerance = 1e-12)

  # And the same values worked out from the table's definition rather than from its own lookup, so
  # the two routes to the answer are genuinely independent.
  by_hand <- vapply(
    seq_along(birth_years),
    function(i) affine_log_mu(-4.2, 0.09, -0.015, i, end_clicks[i] / clicks_per_year),
    double(1)
  )
  expect_equal(res$died_value, by_hand, tolerance = 1e-12)
})

test_that("an affine table integrates to its analytic value", {
  base <- -4.2
  age_slope <- 0.09
  period_slope <- -0.015
  res <- integrate_obj(affine_table(base, age_slope, period_slope))

  from_years <- start_clicks / clicks_per_year
  to_years <- end_clicks / clicks_per_year

  # The integral of `intercept + slope * (t - t0)` over the exposure.
  intercept <- base - age_slope * cohort_index
  slope <- age_slope + period_slope
  analytic <- intercept * (to_years - from_years) +
    slope * analytic_linear(from_years - t0_years, to_years - t0_years)

  expect_equal(res$integral, analytic, tolerance = 1e-12)
})

test_that("log mu is exact at a finer interval too, since the integrand is linear", {
  tbl <- affine_table(base = -4.2, age_slope = 0.09, period_slope = -0.015)

  monthly <- integrate_obj(tbl, time_scale = month_scale)
  quarterly <- integrate_obj(tbl, time_scale = quarter_scale)

  expect_equal(monthly$integral, quarterly$integral, tolerance = 1e-12)
  expect_true(all(monthly$slot_count > quarterly$slot_count))
})

test_that("a table composes with arithmetic in a mortality expression", {
  value <- -3.5
  # Bound to a name first: `mortality()` captures its argument unevaluated, and the folding step
  # resolves a name to the object behind it rather than evaluating a call in the middle of the tree.
  base_table <- constant_table(value)
  res <- integrate_obj(mortality(base_table + 0.25))

  expect_true("vector_log_mu" %in% res$monikers)
  expect_equal(res$integral, (value + 0.25) * exposure_years, tolerance = 1e-12)
})

test_that("an age past the end of the table is clamped, as log_mu() clamps it", {
  # Born long before the table's earliest cohort, so the age index runs off the end throughout.
  old_cols <- cols
  old_cols$birth <- datey::datey(rep(1900, 3))

  tbl <- affine_table(base = -4.2, age_slope = 0.09, period_slope = -0.015)
  res <- cpp_veil_integrate(it_obj(tbl), old_cols, quarter_scale, NULL)

  expected <- vapply(
    seq_len(3),
    function(i) log_mu(tbl, list(birth = old_cols$birth[i]), datey_clicks(end_clicks[i])),
    double(1)
  )
  expect_equal(res$died_value, expected, tolerance = 1e-12)
})

test_that("a table on its own has no single value per individual", {
  expect_error(
    cpp_veil_eval(it_obj(constant_table(-3.5)), cols),
    "varies over time"
  )
})
