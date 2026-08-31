# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the include clip: the exposure narrowed to the interval an include resolves to, and the death
# flag recomputed to match.
#
# THE ORACLE IS `period_included()`, the package's own plain-R resolution of an include. Where the
# exposure is wide enough to contain the whole resolved interval, the integral of 1 is that interval's
# length, so `durationy(period_included(inc, .i))` gives the expected answer with no arithmetic
# restated on this side at all. Where the clip actually bites, the specification's own formula is
# used -- nu' = max(nu, from), tau' = min(tau, to) -- because that is the thing being tested.
#
# `integrate(~ 1)` is the exposure length and `died_value(~ 1)` is 1 for a death and 0 otherwise, so
# between them they report exactly what the clip did to each end.

clicks_per_year <- 534360L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

birth_years <- c(1940, 1945, 1950)

# Deliberately wide, so that an include sitting inside it is not clipped by the exposure and the
# integral reports the include's own length.
wide_start <- rep(1990, 3) * clicks_per_year
wide_end <- rep(2060, 3) * clicks_per_year

cols <- list(
  birth      = datey::datey(birth_years),
  retirement = datey::datey(c(2005, NA, 2012)),
  smoker     = c(TRUE, FALSE, TRUE),
  E2R_start  = datey_clicks(wide_start),
  E2R_end    = datey_clicks(wide_end),
  E2R_died   = c(TRUE, TRUE, TRUE)
)

facts <- function(i) {
  list(birth = cols$birth[i], retirement = cols$retirement[i], smoker = cols$smoker[i])
}

# The length of the interval an include resolves to, in years, and zero where it resolves to nothing.
resolved_years <- function(inc) {
  vapply(seq_along(birth_years), function(i) {
    interval <- period_included(inc, facts(i))
    if (is.na(interval)) 0 else unclass(datey::durationy(interval)) / clicks_per_year
  }, double(1))
}

integrate_one <- function(include = NULL, columns = cols) {
  cpp_veil_integrate(it_ast(~ 1), columns, quarter_scale, include)
}

test_that("no include leaves the exposure alone", {
  res <- integrate_one(NULL)

  expect_equal(res$integral, (wide_end - wide_start) / clicks_per_year, tolerance = 1e-12)
  expect_equal(res$died_value, rep(1, 3), tolerance = 1e-12)
})

test_that("an absolute period clips to its own bounds", {
  inc <- period(2010, 2020)
  res <- integrate_one(inc)

  expect_equal(res$integral, resolved_years(inc), tolerance = 1e-12)
  expect_equal(res$integral, rep(10, 3), tolerance = 1e-12)
})

test_that("an age band clips per individual", {
  inc <- age(65, 95)
  res <- integrate_one(inc)

  expect_equal(res$integral, resolved_years(inc), tolerance = 1e-12)
  expect_equal(res$integral, rep(30, 3), tolerance = 1e-12)
})

test_that("an offset band uses the field it names", {
  inc <- band(.t - .i$retirement, 0, 5)
  res <- integrate_one(inc)

  # The middle individual has no retirement date, so contributes nothing -- and gets there without
  # the arithmetic ever touching R's NA integer, which would overflow rather than mean anything.
  expect_equal(res$integral, resolved_years(inc), tolerance = 1e-12)
  expect_equal(res$integral, c(5, 0, 5), tolerance = 1e-12)
})

# AN OFFSET IS AN EXPRESSION, NOT A COLUMN. These two run it through the engine rather than through
# `period_included()` alone, since the R oracle would agree with itself whatever the boundary did.
test_that("a computed offset agrees with the bare field it reduces to", {
  # Every retirement here is after the birth, so `max` of the two IS the retirement. The computed
  # form must therefore give exactly what naming the field gives.
  computed <- band(.t - max(.i$birth, .i$retirement), 0, 5)

  expect_equal(integrate_one(computed)$integral,
               integrate_one(band(.t - .i$retirement, 0, 5))$integral,
               tolerance = 1e-12)
  expect_equal(integrate_one(computed)$integral, resolved_years(computed), tolerance = 1e-12)
})

test_that("a computed offset is really computed, not reduced to one of its fields", {
  # `min` of birth and retirement is the birth by value -- but it is MISSING wherever retirement is,
  # because a missing value poisons the min. So this must differ from banding birth directly, which
  # is what would happen if the expression were being quietly resolved to a single column.
  # Ages 60-65 rather than 0-5, so that every band lands inside the exposure window and a zero means
  # the offset was missing rather than that the interval fell outside the data.
  computed <- band(.t - min(.i$birth, .i$retirement), 60, 65)

  expect_equal(integrate_one(computed)$integral, c(5, 0, 5), tolerance = 1e-12)
  expect_equal(integrate_one(age(60, 65))$integral, c(5, 5, 5), tolerance = 1e-12)
  expect_equal(integrate_one(computed)$integral, resolved_years(computed), tolerance = 1e-12)
})

test_that("intersecting includes applies both", {
  inc <- age(65, 95) & period(2010, 2020)
  res <- integrate_one(inc)

  expect_equal(res$integral, resolved_years(inc), tolerance = 1e-12)

  # Born 1940 is 70 to 80 across the period, so all ten years count; born 1950 turns 65 in 2015, so
  # only the second half does.
  expect_equal(res$integral, c(10, 10, 5), tolerance = 1e-12)
})

test_that("a gate empties an individual entirely", {
  inc <- period(2010, 2020) & ~ .i$smoker
  res <- integrate_one(inc)

  expect_equal(res$integral, resolved_years(inc), tolerance = 1e-12)
  expect_equal(res$integral, c(10, 0, 10), tolerance = 1e-12)

  # Nobody counts as a death here, gate or no gate: the include ends in 2020 and the exposure runs to
  # 2060, so the clip cut the end off for all three.
  expect_equal(res$died_value, rep(0, 3), tolerance = 1e-12)
})

test_that("an include disjoint from the exposure contributes nothing rather than erroring", {
  narrow <- cols
  narrow$E2R_start <- datey_clicks(rep(2030, 3) * clicks_per_year)
  narrow$E2R_end <- datey_clicks(rep(2040, 3) * clicks_per_year)

  res <- integrate_one(period(2010, 2020), narrow)

  expect_equal(res$integral, rep(0, 3), tolerance = 1e-12)
  expect_equal(res$died_value, rep(0, 3), tolerance = 1e-12)
  expect_equal(res$slot_count, rep(0L, 3))
})

test_that("the clip bites at whichever end the exposure does not cover", {
  # Exposure 2015 to 2025 against an include of 2010 to 2020: the start comes from the exposure, the
  # end from the include, so the death at the exposure's end is cut off.
  narrow <- cols
  narrow$E2R_start <- datey_clicks(rep(2015, 3) * clicks_per_year)
  narrow$E2R_end <- datey_clicks(rep(2025, 3) * clicks_per_year)

  res <- integrate_one(period(2010, 2020), narrow)

  expect_equal(res$integral, rep(5, 3), tolerance = 1e-12)
  expect_equal(res$died_value, rep(0, 3), tolerance = 1e-12)
})

test_that("a death survives a clip that does not cut the end off", {
  # The include ends after the exposure does, so tau is the exposure's own end and the death counts.
  narrow <- cols
  narrow$E2R_start <- datey_clicks(rep(2015, 3) * clicks_per_year)
  narrow$E2R_end <- datey_clicks(rep(2018, 3) * clicks_per_year)

  res <- integrate_one(period(2010, 2020), narrow)

  expect_equal(res$integral, rep(3, 3), tolerance = 1e-12)
  expect_equal(res$died_value, rep(1, 3), tolerance = 1e-12)
})

test_that("a survivor stays a survivor whatever the clip does", {
  alive <- cols
  alive$E2R_died <- c(FALSE, FALSE, FALSE)

  res <- integrate_one(period(2010, 2020), alive)

  expect_equal(res$integral, rep(10, 3), tolerance = 1e-12)
  expect_equal(res$died_value, rep(0, 3), tolerance = 1e-12)
})

test_that("the clip applies to a real integrand, not just to a constant", {
  # Age integrates over the clipped window rather than the whole exposure, which is the point: the
  # clip moves where the time vector starts and ends, so everything sampled on it moves with it.
  inc <- period(2010, 2020)
  res <- cpp_veil_integrate(it_ast(~ .x), cols, quarter_scale, inc)

  from <- 2010 - birth_years
  to <- 2020 - birth_years
  expect_equal(res$integral, (to - from) * (to + from) / 2, tolerance = 1e-9)
})

test_that("something that is not an include is refused", {
  expect_error(integrate_one(mortality_const(log_mu = -4.5)), "must be an `include`")
})
