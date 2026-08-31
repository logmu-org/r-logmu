# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the time vector: the sample points an individual's exposure is cut into, the integral over
# them, and the value taken at death.
#
# THE ORACLE IS THE ANALYTIC INTEGRAL, not a second implementation of the midpoint rule. The
# rectangular midpoint rule is exact for an integrand that is constant or linear in time, on every
# subinterval and so on any composition of them, which means the expected answer for those cases can
# be written down: the integral of 1 is the exposure length, and the integral of age is the
# difference of two halved squares. Anything a re-implementation would share a mistake with is
# therefore excluded.
#
# WHERE EXACTNESS STOPS. Sample points are whole clicks, which keeps them reproducible and lets an
# age be an exact integer subtraction. Half of the integration interval is itself a whole number of
# clicks, so every full interval has an exact midpoint. The FINAL, short interval need not: its
# midpoint is half of `nu + n*dt + tau`, and when that sum is odd the halving moves the sample half a
# click -- about 29 seconds. The resulting error in the integral is at most half a click times the
# width of that final interval, it applies to one interval however many there are, and the test
# below states it as a bound rather than pretending it is not there.

clicks_per_year <- 534360
quarter <- clicks_per_year / 4

# Clicks in, clicks out: the exposures are built by integer arithmetic so that no test value depends
# on a year-to-click conversion rounding the way the test happens to expect.
years_of <- function(clicks) clicks * (1 / clicks_per_year)

start_clicks <- c(2010, 2010, 2011) * clicks_per_year
end_clicks <- c(
  2013 * clicks_per_year,                                    # exactly 12 quarters
  2010 * clicks_per_year + 12 * quarter + clicks_per_year / 10, # 12 quarters and a short tail
  2011 * clicks_per_year + clicks_per_year / 8               # shorter than one interval
)

cols <- list(
  birth     = datey::datey(c(1960, 1965, 1975)),
  E2R_start = structure(as.integer(start_clicks), class = class(datey::datey(2010))),
  E2R_end   = structure(as.integer(end_clicks), class = class(datey::datey(2010))),
  E2R_died  = c(FALSE, TRUE, FALSE),
  weight    = c(1.5, 2.5, 0.5)
)

birth_clicks <- unclass(cols$birth)
exposure_years <- years_of(end_clicks - start_clicks)

# The half-click bound on the final interval, per individual. Zero where the exposure is a whole
# number of intervals, because then there is no final interval at all.
final_width_years <- years_of((end_clicks - start_clicks) %% quarter)
midpoint_bound <- 0.5 * (1 / clicks_per_year) * final_width_years

# The integral of a linear integrand, written so the two squares are not subtracted from each other.
# Over calendar time they are both around four million and their difference is a few thousand, so the
# obvious spelling loses most of its precision to cancellation and would make a worse oracle than the
# thing it is checking.
analytic_linear <- function(from, to) (to - from) * (to + from) / 2

# Asserts an integral against its analytic value, allowing for the half-click shift of the final
# sample point. `slope` is the derivative of the integrand, since the shift feeds through by it.
expect_matches_analytic <- function(got, analytic, slope = 1) {
  bound <- abs(slope) * midpoint_bound + 1e-11 * (1 + abs(analytic))
  expect_true(all(abs(got - analytic) <= bound))
}

integrate_of <- function(expr, time_scale = quarter_scale) {
  cpp_veil_integrate(expr, cols, time_scale, NULL)
}

test_that("the exposure is cut into the slots the specification describes", {
  res <- integrate_of(it_ast(~ 1))
  # Twelve whole quarters; then twelve plus a short tail plus a death slot; then a single short
  # interval on its own.
  expect_equal(res$slot_count, c(12L, 14L, 1L))
  expect_equal(res$delta_t_clicks, quarter)
})

test_that("a finer interval cuts the same exposure into more slots", {
  # Monthly: 36 whole months; then 37 whole plus a short tail plus the death slot; then one whole
  # month plus a tail. Yearly: 3; then 3 plus a tail plus death; then a single short interval.
  expect_equal(integrate_of(it_ast(~ 1), month_scale)$slot_count, c(36L, 39L, 2L))
  expect_equal(integrate_of(it_ast(~ 1), year_scale)$slot_count, c(3L, 5L, 1L))
})

test_that("an unsupported interval is refused", {
  # A CLICK COUNT THAT IS NOT ONE OF THE FOUR, and 7 in particular: it is what a call site left
  # behind by the rename to clicks would pass, having meant seven intervals a year. Seven clicks is
  # about two minutes, so without this guard the answer would be wrong and entirely plausible.
  expect_error(integrate_of(it_ast(~ 1), 7L), "1, 1/4, 1/12 or 1/60")

  # The old spelling is refused for the same reason, which is the point of the guard.
  expect_error(integrate_of(it_ast(~ 1), 4L), "1, 1/4, 1/12 or 1/60")
})

test_that("integrating one gives the length of the exposure", {
  # Exact: the midpoint rule reproduces a constant integrand on every interval, and the widths sum to
  # the exposure however the intervals fall.
  expect_equal(integrate_of(it_ast(~ 1))$integral, exposure_years)
  expect_equal(integrate_of(it_ast(~ 1), fine_scale)$integral, exposure_years)
})

test_that("integrating a constant scales the exposure by it", {
  expect_equal(integrate_of(it_ast(~ exp(-4.5)))$integral, exp(-4.5) * exposure_years)
  expect_equal(integrate_of(it_ast(~ .i$weight))$integral, cols$weight * exposure_years)
})

test_that("integrating age matches the analytic integral", {
  # The integral of (t - b) from nu to tau is the difference of the halved squares of the ages at
  # each end. Linear in t, so the midpoint rule is exact but for the final-interval bound.
  from <- years_of(start_clicks - birth_clicks)
  to <- years_of(end_clicks - birth_clicks)

  res <- integrate_of(it_ast(~ .x))
  expect_matches_analytic(res$integral, analytic_linear(from, to))

  # The third individual is the one whose final interval has an odd click span, so it is the one that
  # actually uses the bound. Left as an assertion so that a change making it exact, or making it very
  # much worse, both show up here.
  expect_true(abs(res$integral[[3]] - analytic_linear(from, to)[[3]]) > 0)
  expect_equal(res$integral[[1]], analytic_linear(from, to)[[1]])

  # `.x` is one operation, not a time vector followed by a subtraction: the subtraction stays in
  # exact clicks and only the difference is converted.
  expect_equal(res$monikers, c("vector_durn", "integrate"))
})

test_that("integrating time itself matches the analytic integral", {
  from <- years_of(start_clicks)
  to <- years_of(end_clicks)
  res <- integrate_of(it_ast(~ .t))
  expect_matches_analytic(res$integral, analytic_linear(from, to))
  expect_equal(res$monikers, c("vector_t", "integrate"))
})

test_that("a scalar mixes into a time-varying expression without a broadcast", {
  # weight does not vary over time and age does, so the product is per-slot with the weight read at
  # every slot. No explicit broadcast is emitted; the interpreter reads the scalar in place.
  #
  # The weight is also the integrand's slope, so it scales the half-click bound with it.
  from <- years_of(start_clicks - birth_clicks)
  to <- years_of(end_clicks - birth_clicks)
  res <- integrate_of(it_ast(~ .x * .i$weight))
  expect_matches_analytic(res$integral, cols$weight * analytic_linear(from, to), slope = cols$weight)
  expect_false("vector" %in% res$monikers)
})

test_that("a constant integrand is broadcast, because there is nothing to integrate over otherwise", {
  # Nothing here varies over time, so the value has to be spread across the slots before the integral
  # can read it. That is the one place the explicit broadcast op appears.
  res <- integrate_of(it_ast(~ .i$weight))
  expect_equal(res$monikers, c("vector", "integrate"))
})

test_that("a comparison over time answers one or zero at each slot", {
  # Every individual is over 20 throughout, so the indicator is 1 at every slot and its integral is
  # the whole exposure. An indicator weight is exactly this shape.
  expect_equal(integrate_of(it_ast(~ .x > 20))$integral, exposure_years)
  # And nobody is over 200, so that indicator integrates to nothing.
  expect_equal(integrate_of(it_ast(~ .x > 200))$integral, rep(0, 3))
})

test_that("an indicator that turns on partway through integrates to the part it covers", {
  # The first individual is 50 at the start of 2010 and 53 at the end of 2013, so `.x > 51` covers
  # the exposure from age 51 onwards. Quarterly slots put the boundary on a slot edge here, so the
  # midpoint rule gets it exactly right rather than approximately.
  res <- integrate_of(it_ast(~ .x > 51))
  age_at_start <- years_of(start_clicks - birth_clicks)
  age_at_end <- years_of(end_clicks - birth_clicks)
  covered <- pmax(0, pmin(age_at_end, Inf) - pmax(age_at_start, 51))
  expect_equal(res$integral[[1]], covered[[1]], tolerance = 1e-9)
})

test_that("the value at death is the integrand at the end of the exposure, and nothing otherwise", {
  # Only the second individual died. Their age at death is exact -- it is a whole-click subtraction,
  # with no midpoint involved.
  res <- integrate_of(it_ast(~ .x))
  age_at_end <- years_of(end_clicks - birth_clicks)
  expect_equal(res$died_value, c(0, age_at_end[[2]], 0))

  res <- integrate_of(it_ast(~ 1))
  expect_equal(res$died_value, c(0, 1, 0))
})

test_that("a date less a number of years is not the exact click subtraction", {
  # `.t - .b` is datey minus datey and stays in clicks; `.t - 5` is a date less a number of YEARS,
  # which the datey rules already make a double. The second must not be lowered as a duration.
  expect_equal(integrate_of(it_ast(~ .t - .b))$monikers, c("vector_durn", "integrate"))
  expect_equal(integrate_of(it_ast(~ .t - 5))$monikers, c("vector_t", "sub", "integrate"))
})

test_that("integrating needs the exposure columns", {
  bare <- list(birth = cols$birth, weight = cols$weight)
  expect_error(cpp_veil_integrate(it_ast(~ 1), bare, quarter_scale, NULL), "E2R_start")
})

test_that("a time-varying expression has no single value per individual", {
  # The scalar entry point is handed the exposure, because these columns carry one, so it gets as far
  # as discovering the real problem: the expression produces a value per slot, not per individual.
  expect_error(cpp_veil_eval(it_ast(~ .t), cols), "no single value per individual")
  expect_error(cpp_veil_eval(it_ast(~ .x), cols), "no single value per individual")

  # A time-invariant expression over the same data is unaffected.
  expect_equal(cpp_veil_eval(it_ast(~ .i$weight * 2), cols)$values, cols$weight * 2)
})

test_that("a logical operator over time is refused rather than guessed at", {
  # The op table says the logical operators are not vectorisable, and a vector holds doubles, so
  # there is nowhere for a per-slot logical to live. Refused by name rather than half-handled.
  expect_error(integrate_of(it_ast(~ !(.x > 20))), "not defined over the time vector")
})
