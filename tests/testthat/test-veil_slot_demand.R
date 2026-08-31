# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests slot demand: which parts of the time vector a value is actually read at, and so where it
# needs computing.
#
# The grid is not one thing. `integrate` reads the midpoint slots and stops before the death slot;
# `died_value` reads the death slot and nothing else. A value feeding only the second is wanted at
# one slot out of the twelve to twenty an exposure is cut into -- and at none at all for an
# individual who did not die, which is most of them.
#
# THE FAILURE THIS GUARDS AGAINST IS STALE DATA, not a wrong sum. When a survivor's death-only
# instruction is skipped, its buffer still holds whatever the last individual left there. If
# `died_value` were ever to read it the answer would be the PREVIOUS individual's, which is a
# plausible number rather than an obvious fault. So the fixture puts a death before a survivor, and
# gives them weights far enough apart that one showing up in the other's place could not be missed.

clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

start_clicks <- c(2010, 2010, 2010) * clicks_per_year
end_clicks <- start_clicks + c(12L, 8L, 4L) * quarter

# The first individual dies, the second does not, and their weights are nothing like each other.
cols <- list(
  birth     = datey::datey(c(1940, 1945, 1950)),
  amount    = c(1000, 2500, 400),
  E2R_start = datey_clicks(start_clicks),
  E2R_end   = datey_clicks(end_clicks),
  E2R_died  = c(TRUE, FALSE, TRUE)
)

exposure_years <- (end_clicks - start_clicks) / clicks_per_year
age_to <- (end_clicks - unclass(cols$birth)) / clicks_per_year
age_from <- (start_clicks - unclass(cols$birth)) / clicks_per_year
integral_of_age <- (age_to - age_from) * (age_to + age_from) / 2

log_mu_value <- -3.2
mu_value <- exp(log_mu_value)

flat_table <- mortality_table(
  x0 = 55, t0 = 2005,
  log_mu = matrix(log_mu_value, nrow = 40L, ncol = 20L)
)

test_that("A's integrand is wanted at the death slot alone", {
  res <- cpp_veil_aev(it_obj(flat_table), it_ast(~ .i$amount), cols, quarter_scale, NULL, no_overdispersion, 1L)

  # One of the three time vectors -- the one A reads -- is never wanted at a midpoint. The other two
  # are log mu over the grid and its exponential, which E integrates.
  expect_equal(res$vector_operand_count, 3L)
  expect_equal(res$death_only_count, 1L)

  # THE SLOT COUNT IS WHAT SHOWS THE NARROWING ACTUALLY HAPPENED, and nothing else can: computing a
  # slot nobody reads changes no answer, so only the work done tells you whether it was skipped.
  # Exposures are 12, 8 and 4 quarters, and the first and third individuals died, so the grids hold
  # 13, 8 and 5 slots. Log mu and its exponential are wanted at the midpoints alone -- 12, 8 and 4
  # -- and A's weight at the death slot alone, which is one slot for a death and none for the
  # survivor. That is 25 + 16 + 9. Filling every vector over every slot would be 3 * 26 = 78.
  expect_equal(res$slot_evaluations, 50L)

  # A IS THE WEIGHT FOR THE ONES WHO DIED AND NOTHING FOR THE SURVIVOR. The middle individual is the
  # one that matters: their death-only instruction never ran, so a zero here says the buffer left
  # over from the individual before was not read.
  expect_equal(res$contributions$A, cols$amount * as.double(cols$E2R_died), tolerance = 1e-12)
  expect_equal(res$contributions$A[2], 0)

  expect_equal(res$contributions$E, mu_value * cols$amount * exposure_years, tolerance = 1e-12)
  expect_equal(res$contributions$V, mu_value * cols$amount^2 * exposure_years, tolerance = 1e-12)
})

test_that("a weight A and E both read is wanted everywhere", {
  # `.x * 1` reaches `died_value` for A and the integrand for E, and sharing has made those the same
  # operand. So it is genuinely wanted at both ends and nothing narrows -- which is the right answer,
  # not a missed opportunity: computing it twice would cost the same work and more instructions.
  res <- cpp_veil_aev(it_obj(flat_table), it_ast(~ .x * 1), cols, quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(res$death_only_count, 0L)

  # A is the age at death, which is where the exposure ends for the two who died.
  expect_equal(res$contributions$A, ifelse(cols$E2R_died, age_to, 0), tolerance = 1e-12)
  expect_equal(res$contributions$E, mu_value * integral_of_age, tolerance = 1e-12)
})

test_that("a survivor's death value is nothing however the individual before them died", {
  # The same calculation run over the individuals in both orders. Slot demand means the buffer A
  # reads is written only for a death, so if anything were reading it unwritten the answer for a
  # given individual would depend on who came before them. It must not.
  reversed <- lapply(cols, rev)

  forwards <- cpp_veil_aev(it_obj(flat_table), it_ast(~ .i$amount), cols, quarter_scale, NULL, no_overdispersion, 1L)
  backwards <- cpp_veil_aev(it_obj(flat_table), it_ast(~ .i$amount), reversed, quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(rev(backwards$contributions$A), forwards$contributions$A, tolerance = 1e-12)
  expect_equal(rev(backwards$contributions$E), forwards$contributions$E, tolerance = 1e-12)
  expect_equal(rev(backwards$contributions$V), forwards$contributions$V, tolerance = 1e-12)
})

test_that("an integral wants the midpoints and not the death slot", {
  # Nothing here reads a death slot at all, so no vector is death-only and the death slot of a
  # deceased individual's grid is never filled. The integral is over age, which is linear in time
  # and so exact under the midpoint rule.
  res <- cpp_veil_integrate(it_ast(~ .x), cols, quarter_scale, NULL)

  expect_equal(res$death_only_count, 0L)

  # One vector, the age, filled at the 12, 8 and 4 midpoints and not at the two death slots. Filling
  # the whole grid would be 26.
  expect_equal(res$slot_evaluations, 24L)

  expect_equal(res$integral, integral_of_age, tolerance = 1e-12)

  # The same expression taken at death rather than integrated, which is the other half of the grid.
  expect_equal(res$died_value, ifelse(cols$E2R_died, age_to, 0), tolerance = 1e-12)
})
