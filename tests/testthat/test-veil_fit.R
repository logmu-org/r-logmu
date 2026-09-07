# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the fit recipe: everything one Newton-Raphson iteration needs, out of one
# pass over the data, at the beta already folded into the mortality.
#
#     A                = died_value(w log mu)
#     E                = integrate(mu w)
#     score_actual[j]  = died_value(w X_j)
#     score_expected[j]= integrate(mu w X_j)
#     ew_xx[j,l]       = integrate(mu w X_j X_l)
#     ew2_xx[j,l]      = integrate(mu w^2 X_j X_l)
#
# THE ORACLES ARE MOSTLY OTHER RECIPES, WHICH IS THE POINT OF TESTING THIS SLICE
# ON ITS OWN. With no terms at all the block reduces to the log-likelihood, and
# its E must agree with the E an A/E computes over the same data -- to the BIT,
# since it is the same integrand reached by a different assembly. With a single
# constant term of one, every family collapses onto an A/E's own three answers.
# Only where a term varies over time is an oracle worked out on paper here, and
# then it is the midpoint sum rather than the true integral, because the
# midpoint rule is what the engine promises.
#
# Exposures are whole quarters, so there is no short final interval.

clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

birth_years <- c(1940, 1945, 1950)
birth_clicks <- birth_years * clicks_per_year
start_clicks <- c(2010, 2010, 2010) * clicks_per_year
end_clicks <- start_clicks + c(12L, 8L, 4L) * quarter

cols <- list(
  birth     = datey::datey(birth_years),
  amount    = c(1000, 2500, 400),
  male      = c(TRUE, FALSE, TRUE),
  E2R_start = datey_clicks(start_clicks),
  E2R_end   = datey_clicks(end_clicks),
  E2R_died  = c(TRUE, FALSE, TRUE)
)

exposure_years <- (end_clicks - start_clicks) / clicks_per_year
died <- cols$E2R_died
quarters_of <- (end_clicks - start_clicks) %/% quarter

log_mu_value <- -3.2
mu_value <- exp(log_mu_value)
constant_mortality <- mortality_const(log_mu = log_mu_value)

# `mortality` is the BASE, and the recipe adds `sum_j beta_j X_j` to it. Beta defaults to zeros, so
# an oracle written against the base mortality alone still holds: `0 * X` is exactly zero for a
# finite covariate and adding it changes no bit.
fit <- function(terms = list(), mortality = it_obj(constant_mortality), weight = NULL,
                include = NULL, time_scale = quarter_scale, columns = cols,
                keep_contributions = TRUE, threads = 1L,
                beta = rep(0, length(terms))) {
  cpp_veil_fit(mortality, terms, as.double(beta), weight, NULL, NULL, columns, time_scale, include,
               keep_contributions, threads)
}

aev <- function(weight = NULL, include = NULL, time_scale = quarter_scale, columns = cols) {
  cpp_veil_aev(it_obj(constant_mortality), weight, columns, time_scale, include,
               no_overdispersion, 1L)
}

# The midpoints of an individual's quarterly grid, as ages in years. The engine
# builds exactly these, and test-veil_time.R is what holds it to them; here they
# are the independent oracle for a term that varies over time.
midpoint_ages <- function(individual) {
  steps <- seq_len(quarters_of[[individual]]) - 1L
  midpoints <- start_clicks[[individual]] + steps * quarter + quarter / 2
  (midpoints - birth_clicks[[individual]]) / clicks_per_year
}

test_that("with no terms the fit reduces to the log-likelihood, and E is an A/E's E to the bit", {
  # THE SHARPEST ORACLE AVAILABLE, and it is exact rather than approximate: E is
  # `integrate(mu * w)` in both recipes, so two different assemblies of the same
  # integrand must produce the identical double, not merely a close one.
  res <- fit(weight = it_ast(~ .i$amount))
  reference <- aev(weight = it_ast(~ .i$amount))

  expect_identical(res$E, reference$E)
  expect_identical(res$term_count, 0L)

  # Two roots and nothing else.
  expect_identical(res$output_count, 2L)

  # A is the weight times log mu at the death, and nothing for a survivor -- so
  # it is the A/E's own A scaled by log mu, which is constant here.
  expect_equal(res$A, reference$A * log_mu_value, tolerance = 1e-12)
})

test_that("a single constant term of one collapses every family onto an A/E", {
  # X = 1 makes the score's expected part, the information and the A/E's own E
  # the same integral, and the second-moment triangle the A/E's V.
  res <- fit(terms = list(it_ast(~ 1)), weight = it_ast(~ .i$amount))
  reference <- aev(weight = it_ast(~ .i$amount))

  expect_identical(res$term_count, 1L)
  # 2 + 2k + k(k+1) with k = 1.
  expect_identical(res$output_count, 6L)

  expect_equal(res$score_actual, reference$A, tolerance = 1e-12)
  expect_equal(res$score_expected, reference$E, tolerance = 1e-12)

  # THE SHARING IS THE ASSERTION HERE, not merely the arithmetic. `X * X` folds
  # to one, the integrand becomes the same node as E's, and the two outputs end
  # up reading one operand -- so this is bit-identical or the fold did not fire.
  expect_identical(res$ew_xx, res$E)
  expect_equal(res$ew2_xx, reference$V, tolerance = 1e-12)
})

test_that("the packed triangle is upper, row-major, and symmetric in its two indices", {
  # Constant terms of 2 and 3, so the six integrals are the plain E scaled by
  # 4, 6 and 9 -- and the middle one is the only place a transposed index would
  # not show up as an obviously wrong magnitude.
  res <- fit(terms = list(it_ast(~ 2), it_ast(~ 3)), weight = it_ast(~ .i$amount))
  reference <- aev(weight = it_ast(~ .i$amount))

  expect_identical(res$term_count, 2L)
  # 2 + 2k + k(k+1) with k = 2.
  expect_identical(res$output_count, 12L)
  expect_length(res$ew_xx, 3L)

  expect_equal(res$ew_xx, reference$E * c(4, 6, 9), tolerance = 1e-12)
  expect_equal(res$score_expected, reference$E * c(2, 3), tolerance = 1e-12)
  expect_equal(res$score_actual, reference$A * c(2, 3), tolerance = 1e-12)
})

test_that("three terms pin the packing order, which two terms cannot", {
  # AT k = 2 A ROW-MAJOR AND A COLUMN-MAJOR TRIANGLE ARE THE SAME SEQUENCE, so
  # the two-term test above cannot see a transposed packing at all. Three terms
  # is the smallest case where the two orders differ, and they differ in exactly
  # one adjacent pair -- row-major gives (0,2) before (1,1), column-major the
  # reverse.
  res <- fit(terms = list(it_ast(~ 2), it_ast(~ 3), it_ast(~ 5)),
             weight = it_ast(~ .i$amount))
  reference <- aev(weight = it_ast(~ .i$amount))

  expect_length(res$ew_xx, 6L)
  # (0,0) (0,1) (0,2) (1,1) (1,2) (2,2), so 4, 6, 10, 9, 15, 25 -- and NOT
  # 4, 6, 9, 10, 15, 25, which is what column-major would give.
  expect_equal(res$ew_xx, reference$E * c(4, 6, 10, 9, 15, 25), tolerance = 1e-12)
})

test_that("time-invariant covariates cost no vector work however many there are", {
  # THE SOLE WITNESS TO THE HOIST, and the reason the integrands are spelled with
  # `mu` innermost. Where the covariates and the weight do not vary over an
  # individual's exposure, every one of the k(k+1)/2 information integrals should
  # peel down to a scalar multiple of the SAME `integrate(mu)` -- so the vector
  # work must not grow with the number of parameters, even though the number of
  # outputs does.
  #
  # THE MORTALITY HAS TO VARY WITH TIME FOR THIS TO MEAN ANYTHING. Against a
  # constant mortality the whole integrand is time-invariant, and the hoist
  # DECLINES it: lifting the last factor out would need the exposure length,
  # which is a per-individual value rather than a factor. That is why a Gompertz
  # expression is used here and not `constant_mortality`.
  # THE THREE TERMS MUST BE DISTINCT, and the first version of this test got that
  # wrong. Repeating one term three times leaves the vector work flat whatever
  # the hoist does, because the sharing pass merges the identical covariates and
  # the identical products -- the test passed with `passHoistFromIntegrate`
  # switched off entirely, which is what found it.
  gompertz <- it_ast(~ -10.5 + 0.09 * .x)
  distinct <- list(it_ast(~ .i$male), it_ast(~ .i$amount), it_ast(~ .b))

  one <- fit(terms = distinct[1], mortality = gompertz,
             weight = it_ast(~ .i$amount), keep_contributions = FALSE)
  three <- fit(terms = distinct, mortality = gompertz,
               weight = it_ast(~ .i$amount), keep_contributions = FALSE)

  expect_identical(one$output_count, 6L)
  expect_identical(three$output_count, 20L)

  # The outputs more than trebled; the VECTOR work barely moved -- 106 slot
  # evaluations to 110 as measured, against 662 with the hoist switched off.
  # A ratio rather than an equality, because the extra four are real: the
  # score's actual part is a `died_value`, which needs a vector to read at the
  # death slot, so a time-invariant covariate gets broadcast into one.
  expect_lt(three$slot_evaluations, one$slot_evaluations * 1.2)

  # INSTRUCTIONS GROW, AND THAT IS THE HOIST WORKING RATHER THAN FAILING. Every
  # factor peeled out of an integral becomes a scalar multiply outside it, so
  # trading vector work for scalar work is precisely the trade being made. An
  # equality here would be asserting the opposite of what is wanted.
  expect_gt(three$instruction_count, one$instruction_count)
})

test_that("the weight is linear in the information and squared in the second moment", {
  # The one asymmetry between the two triangles, and the reason there are two.
  res <- fit(terms = list(it_ast(~ 1)), weight = it_ast(~ .i$amount))

  amounts <- cols$amount
  expect_equal(res$ew_xx, sum(mu_value * amounts * exposure_years), tolerance = 1e-12)
  expect_equal(res$ew2_xx, sum(mu_value * amounts^2 * exposure_years), tolerance = 1e-12)
})

test_that("an indicator weight collapses the two triangles onto each other", {
  # `w * w` simplifies to `w` where the scan proves the weight is zero or one,
  # so the second-moment triangle becomes the same operand as the information.
  # Bit-identical, since it is one integral read twice.
  res <- fit(terms = list(it_ast(~ 1)), weight = it_ast(~ .i$male))

  expect_identical(res$ew2_xx, res$ew_xx)
})

test_that("a term that varies over time is integrated at the midpoints", {
  # `.x` is age, which moves across the exposure, so nothing hoists and the
  # oracle is the midpoint sum rather than the true integral. The diagonal
  # carries age squared, where the midpoint rule is NOT exact -- which is why
  # the expected value is built the same way the engine promises to build it.
  res <- fit(terms = list(it_ast(~ .x)))

  dt <- 1 / 4
  expected_score <- sum(vapply(seq_along(birth_years), function(i) {
    mu_value * sum(midpoint_ages(i)) * dt
  }, numeric(1)))
  expected_information <- sum(vapply(seq_along(birth_years), function(i) {
    mu_value * sum(midpoint_ages(i)^2) * dt
  }, numeric(1)))

  expect_equal(res$score_expected, expected_score, tolerance = 1e-10)
  expect_equal(res$ew_xx, expected_information, tolerance = 1e-10)

  # The actual part reads the term at the moment of death, which is the exposure
  # end, and is nothing at all for a survivor.
  ages_at_end <- (end_clicks - birth_clicks) / clicks_per_year
  expect_equal(res$score_actual, sum(ages_at_end[died]), tolerance = 1e-10)
})

test_that("every total is the sum of its own per-individual contributions", {
  # The oracle for the accumulation itself rather than for any one integrand: a
  # chunk folded twice, or a boundary out by one, moves a total without moving
  # any individual's value.
  res <- fit(terms = list(it_ast(~ 1), it_ast(~ .x)), weight = it_ast(~ .i$amount))

  contributions <- res$contributions
  expect_identical(dim(contributions), c(res$output_count, length(birth_years)))

  totals <- rowSums(contributions)
  expect_equal(totals[[1]], res$A, tolerance = 1e-12)
  expect_equal(totals[[2]], res$E, tolerance = 1e-12)
  expect_equal(totals[3:4], res$score_actual, tolerance = 1e-12)
  expect_equal(totals[5:6], res$score_expected, tolerance = 1e-12)
  expect_equal(totals[7:9], res$ew_xx, tolerance = 1e-12)
  expect_equal(totals[10:12], res$ew2_xx, tolerance = 1e-12)
})

test_that("an include narrows the fit to the individuals it keeps", {
  kept <- include(.i$male)
  res <- fit(terms = list(it_ast(~ 1)), weight = it_ast(~ .i$amount), include = kept)
  reference <- aev(weight = it_ast(~ .i$amount), include = kept)

  expect_identical(res$records_included, 2L)
  expect_identical(res$E, reference$E)
  expect_equal(res$score_expected, reference$E, tolerance = 1e-12)
})

test_that("the thread count cannot move a single digit", {
  # The strongest test the engine offers, and it is exact rather than to a
  # tolerance: chunk partials fold in chunk order whoever computed them.
  terms <- list(it_ast(~ 1), it_ast(~ .x), it_ast(~ .i$male))
  weight <- it_ast(~ .i$amount)

  one <- fit(terms = terms, weight = weight, threads = 1L, keep_contributions = FALSE)
  many <- fit(terms = terms, weight = weight, threads = 4L, keep_contributions = FALSE)

  expect_identical(one$A, many$A)
  expect_identical(one$E, many$E)
  expect_identical(one$score_actual, many$score_actual)
  expect_identical(one$score_expected, many$score_expected)
  expect_identical(one$ew_xx, many$ew_xx)
  expect_identical(one$ew2_xx, many$ew2_xx)
})

test_that("contributions are off unless asked for", {
  # They cost one double per output per individual, and a real fit has
  # 2 + 2k + k(k+1) outputs -- 462 of them at twenty terms.
  res <- fit(terms = list(it_ast(~ 1)), keep_contributions = FALSE)

  expect_length(res$contributions, 0L)
})

test_that("beta reaches the accumulation through the mortality it was folded into", {
  # THE SLICE'S OWN CONTRACT. There is no beta argument: a caller varies beta by
  # varying log mu, and the recipe is not told which part of it came from the
  # model. Doubling mu must therefore double every integral and leave the
  # actual-death families alone apart from the log.
  shifted <- mortality_const(log_mu = log_mu_value + log(2))

  plain <- fit(terms = list(it_ast(~ 1)), weight = it_ast(~ .i$amount))
  doubled <- fit(terms = list(it_ast(~ 1)), weight = it_ast(~ .i$amount),
                 mortality = it_obj(shifted))

  expect_equal(doubled$E, plain$E * 2, tolerance = 1e-12)
  expect_equal(doubled$score_expected, plain$score_expected * 2, tolerance = 1e-12)
  expect_equal(doubled$ew_xx, plain$ew_xx * 2, tolerance = 1e-12)

  # The score's actual part does not see the mortality at all.
  expect_identical(doubled$score_actual, plain$score_actual)
})
