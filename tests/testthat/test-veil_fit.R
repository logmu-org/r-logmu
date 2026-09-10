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
  result <- cpp_veil_fit(mortality, terms, list(as.double(beta)), weight, NULL, NULL, columns,
                         time_scale, include, keep_contributions, threads)
  # Flattened to one run, so a test that cares about one beta reads the answers directly. The
  # multi-beta shape is what `fit_at()` below exercises.
  c(result[setdiff(names(result), "runs")], result$runs[[1]])
}

# The same block run at SEVERAL betas, which is the shape the Newton loop needs: compile once, set
# the coefficients, run again.
fit_at <- function(betas, terms, ...) {
  cpp_veil_fit(it_obj(constant_mortality), terms, lapply(betas, as.double), NULL, NULL, NULL,
               cols, quarter_scale, NULL, FALSE, 1L)
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

# THE PARAMETER LEAF.
#
# A coefficient is constant across every individual and changes between runs, so it lowers to a
# `ConstantBinding` the block can rewrite. What it must never be is a LITERAL: `passFoldConstants`
# would bake the starting value into a folded product, and `passShareCommonSubtrees` keys a double
# literal by its bits and merges equal ones -- and a fit starts every coefficient at zero, so as
# literals they would all collapse onto one operand and setting one would move the rest.

test_that("each coefficient gets its own parameter slot, however equal their values", {
  # THE MERGE IS INVISIBLE IN THE ARITHMETIC UNTIL A BETA MOVES, which is why this is asserted
  # directly rather than inferred. Three terms all starting at zero is exactly the case that would
  # collapse to one slot were the coefficients literals.
  res <- fit_at(list(c(0, 0, 0)), list(it_ast(~ 1), it_ast(~ 2), it_ast(~ 3)))
  expect_identical(res$term_count, 3L)
  expect_identical(res$parameter_count, 3L)

  # And with no terms there is nothing to parameterise.
  expect_identical(fit_at(list(numeric(0)), list())$parameter_count, 0L)
})

test_that("a coefficient reaches the answer, and only its own term", {
  terms <- list(it_ast(~ 1), it_ast(~ 1))

  # `log mu = -3.2 + b1 * 1 + b2 * 1`, so E scales by exp(b1 + b2). Moving ONE coefficient must
  # scale by exp of that one alone -- if the two slots had merged, setting the first would move both
  # and the factor would be exp(2 * 0.25) instead.
  res <- fit_at(list(c(0, 0), c(0.25, 0), c(0, 0.25), c(0.25, 0.25)), terms)
  base <- res$runs[[1]]$E

  expect_equal(res$runs[[2]]$E, base * exp(0.25), tolerance = 1e-12)
  expect_equal(res$runs[[3]]$E, base * exp(0.25), tolerance = 1e-12)
  expect_equal(res$runs[[4]]$E, base * exp(0.50), tolerance = 1e-12)
})

test_that("a parameter answers exactly what the same value written into the mortality does", {
  # THE ORACLE IS THE OTHER MECHANISM. Setting a coefficient to 0.5 against a constant term of one
  # must give precisely what a mortality of `log_mu + 0.5` gives with no coefficient at all -- the
  # same arithmetic reached two entirely different ways.
  viaParameter <- fit_at(list(0.5), list(it_ast(~ 1)))$runs[[1]]

  shifted <- mortality_const(log_mu = log_mu_value + 0.5)
  viaMortality <- fit(terms = list(), mortality = it_obj(shifted), keep_contributions = FALSE)

  expect_equal(viaParameter$E, viaMortality$E, tolerance = 1e-12)
  expect_equal(viaParameter$A, viaMortality$A, tolerance = 1e-12)
})

test_that("running the same block again at the same beta gives bit-identical answers", {
  # NO STATE MAY LEAK BETWEEN RUNS. The interpreter writes the constants into its registers when it
  # is built for a chunk, so a second run at the same coefficients must reproduce the first exactly
  # -- and coming back to a beta after visiting another must reproduce it too, which is what a
  # backtracking line search does every time it rejects a step.
  res <- fit_at(list(c(0.3, -0.2), c(1.1, 0.7), c(0.3, -0.2)), list(it_ast(~ 1), it_ast(~ 2)))

  expect_identical(res$runs[[1]]$E, res$runs[[3]]$E)
  expect_identical(res$runs[[1]]$ew_xx, res$runs[[3]]$ew_xx)
  expect_identical(res$runs[[1]]$score_expected, res$runs[[3]]$score_expected)
  expect_false(identical(res$runs[[1]]$E, res$runs[[2]]$E))
})

test_that("one compiled block serves every beta", {
  # The structure is settled at compile time and nothing in the loop may move it. Were the block
  # recompiled per beta, folding could give a different shape for different starting values.
  many <- fit_at(list(0, 0.5, -0.5, 2), list(it_ast(~ 1)))
  expect_identical(length(many$runs), 4L)
  expect_identical(many$parameter_count, 1L)

  one <- fit_at(list(0.5), list(it_ast(~ 1)))
  expect_identical(many$instruction_count, one$instruction_count)
  expect_identical(many$output_count, one$output_count)
  expect_identical(many$runs[[2]]$E, one$runs[[1]]$E)
})

test_that("a beta of the wrong length is refused", {
  expect_error(fit_at(list(c(0, 0)), list(it_ast(~ 1))), "one value for each model term")
})

# THE NEWTON-RAPHSON LOOP.
#
# `cpp_veil_fit_run` compiles one block and iterates on it, setting the coefficients between runs.
# Nothing is recompiled, which is what the parameter leaf above exists for.

fit_run <- function(terms = list(it_ast(~ 1)), mortality = it_obj(constant_mortality),
                    weight = it_ast(~ .i$amount), include = NULL, columns = cols,
                    start = rep(0, length(terms)), max_iterations = 25L, tolerance = 1e-6,
                    armijo = 1e-4, max_halvings = 30, overdispersion = no_overdispersion,
                    z_scale = 1, disjoint = FALSE, threads = 1L) {
  cpp_veil_fit_run(mortality, terms, weight, NULL, NULL, columns, quarter_scale, include,
                   as.double(start), as.integer(max_iterations), tolerance, armijo,
                   max_halvings, overdispersion, z_scale, disjoint, as.integer(threads))
}

test_that("a single constant covariate has an exact answer, and the fit finds it", {
  # log mu = log mu_ref + beta, so L is maximised where exp(beta) Ew = Aw and
  #
  #     beta_hat = log(Aw / Ew)
  #
  # taken at beta = 0 -- which is precisely the A and E of an A/E on the reference mortality. An
  # analytic oracle from a different recipe, not a second run of this one.
  reference <- aev(weight = it_ast(~ .i$amount))
  expected <- log(reference$A / reference$E)

  res <- fit_run()

  expect_identical(res$status, "converged")

  # THE BOUND COMES FROM THE STOPPING RULE, not from a round number. Stopping when the gain still
  # available is below eps means the loss in L is at most eps, and that loss is `d^2 I / 2` where I
  # is the information -- so `d <= sqrt(2 eps / I)`. Here I is Ew at the optimum, which is Aw, about
  # 1400, giving about 4e-5. Measured: 8e-10.
  #
  # NOT IN UNITS OF THE REPORTED STANDARD ERROR, which is a different quantity. `sqrt(2 eps)` in
  # standard errors holds only where J = I -- unweighted or indicator-weighted data. With amount
  # weights the sandwich variance here is 1.33 against an information inverse of 7e-4, so the two
  # yardsticks differ by a factor of forty.
  expect_lt(abs(res$beta - expected), sqrt(2 * 1e-6 / reference$A))

  expect_lte(res$predicted_gain, 1e-6)
  expect_true(is.finite(res$log_likelihood))
})

test_that("the answer does not depend on where beta started", {
  # The objective is strictly concave for a non-negative weight, so there is ONE maximum and any
  # start must reach it. Starting far out is also what exercises the damping: a full Newton step
  # from beta = 4 overshoots badly.
  base <- fit_run()
  target <- base$beta
  standardError <- sqrt(base$variance)

  # Each run is within `sqrt(2 eps / I)` of the maximum, so two are within twice that of each other.
  # `standardError` is deliberately NOT the yardstick -- see the note above on why the sandwich
  # variance and the information inverse are different quantities here.
  reference <- aev(weight = it_ast(~ .i$amount))
  bound <- 2 * sqrt(2 * 1e-6 / reference$A)

  for (start in list(0, 2, -2, 4, 8, -8)) {
    res <- fit_run(start = start)
    expect_identical(res$status, "converged")
    expect_lt(abs(res$beta - target), bound)
  }
  expect_true(standardError > 0)
})

test_that("the variance and the penalty match what the A/E says they must be", {
  # For the intercept-only model both collapse onto the reference A/E:
  #
  #     Var(beta_hat) = Omega Ew^2(beta_hat) / Aw^2  =  V / (A E)
  #
  # because Ew(beta_hat) = Aw at the maximum and V already carries Omega, which therefore cancels.
  reference <- aev(weight = it_ast(~ .i$amount))
  res <- fit_run()

  expect_identical(res$status, "converged")
  expect_equal(res$variance, reference$V / (reference$A * reference$E), tolerance = 1e-8)
})

test_that("one parameter costs exactly one on the L scale for an indicator weight", {
  # p = tr(J I^-1) / Z, and where w^2 = w the second moment IS the first, so the trace is the number
  # of parameters and Z is one. This is the fact the whole `L_tolerance` convention rests on -- that a
  # tolerance in units of L means "a fraction of one parameter's worth".
  res <- fit_run(weight = it_ast(~ .i$male), z_scale = 1)

  expect_identical(res$status, "converged")
  expect_equal(res$penalty, 1, tolerance = 1e-10)
})

test_that("an unidentifiable model fails and names the dependency", {
  # THE SAME COVARIATE TWICE. beta_1 + beta_2 is determined but neither one is, so there is no
  # maximum -- and a fit returns a fit or it fails.
  res <- fit_run(terms = list(it_ast(~ 1), it_ast(~ 1)))

  expect_identical(res$status, "not_identifiable")

  # Detected at the SECOND covariate, because the first is fine on its own.
  expect_identical(res$failed_parameter, 1L)

  # And the diagnosis names the earlier one: covariate 2 is 1.0 times covariate 1.
  expect_equal(res$dependency, 1, tolerance = 1e-8)
})

test_that("a covariate with no exposure of its own is refused rather than fitted", {
  # `.i$male & !.i$male` is identically zero, so it contributes nothing anywhere and its own
  # information diagonal is zero. Nothing can be estimated from it.
  res <- fit_run(terms = list(it_ast(~ .i$male & !.i$male)))
  expect_identical(res$status, "not_identifiable")
  expect_identical(res$failed_parameter, 0L)
  expect_identical(length(res$dependency), 0L)
})

test_that("running out of iterations is a failure, not a quiet answer", {
  res <- fit_run(start = 6, max_iterations = 1L)
  expect_identical(res$status, "did_not_converge")

  # The gain still available is reported, so a caller can see how far off it was.
  expect_true(is.finite(res$predicted_gain))
  expect_gt(res$predicted_gain, 1e-6)
})

test_that("damping fires below the optimum and not above it", {
  # THE ASYMMETRY IS `exp`, AND IT IS THE OPPOSITE WAY ROUND FROM "far away means damping". A step
  # that RAISES mortality overshoots, because the exponential grows faster than the quadratic model
  # predicts; a step that LOWERS it undershoots, so the full step is always an improvement and
  # nothing is ever halved. Measured on this data, where the reference table is out by a factor of
  # four so even zero is below the optimum:
  #
  #     start  iterations  walks  halvings
  #        -8           6     20        14
  #         0           5      6         1
  #        +4           7      7         0
  #        +8          11     11         0
  #
  # Note the trade: from below it converges in FEWER iterations but more walks, because from above
  # the Newton step is bounded by about one per iteration.
  below <- fit_run(start = -5)
  above <- fit_run(start = 4)

  expect_identical(below$status, "converged")
  expect_identical(above$status, "converged")

  # A walk per iteration, plus one for the starting point, and one more for every rejected trial.
  expect_gt(below$evaluations, below$iterations)
  expect_identical(above$evaluations, above$iterations)
})

test_that("the block is compiled once however many iterations run", {
  # One parameter slot per term, whatever the iteration count -- the loop never recompiles, so this
  # cannot drift with the number of walks.
  res <- fit_run(terms = list(it_ast(~ 1), it_ast(~ .i$male)), start = c(0, 0))
  expect_identical(res$term_count, 2L)
  expect_identical(res$parameter_count, 2L)
})

test_that("armijo at or above one half is refused", {
  # One half is exactly what a full Newton step delivers on a quadratic, so anything at or above it
  # rejects the full step near the optimum and the iteration can never finish.
  expect_error(fit_run(armijo = 0.5), "strictly between 0 and 0.5")
  expect_error(fit_run(armijo = 0), "strictly between 0 and 0.5")
})

test_that("the predicted gain is the gain actually available, not twice it", {
  # THE HALF IS NOT DECORATION. Along the Newton direction the linear term promises `grad . step`
  # and the curvature gives half of it back, so half is what a full step actually delivers. Take one
  # step from close to the optimum, where the quadratic approximation is good, and the realised gain
  # in L must match `predicted_gain` -- not half it, and not twice it.
  start <- 1.3

  atStart <- cpp_veil_fit(it_obj(constant_mortality), list(it_ast(~ 1)), list(start),
                          it_ast(~ .i$amount), NULL, NULL, cols, quarter_scale, NULL, FALSE, 1L)
  before <- atStart$runs[[1]]$A - atStart$runs[[1]]$E

  stepped <- fit_run(start = start, max_iterations = 1L)
  expect_identical(stepped$status, "did_not_converge")

  # ONE WALK FOR THE STARTING POINT PLUS ONE FOR THE STEP. `evaluations == iterations` holds only
  # when the last iteration CONVERGED, because a converging iteration returns without stepping. Here
  # it stepped and then ran out of budget, so there is one more walk than iterations -- and nothing
  # was halved, which is what matters for the comparison below.
  expect_identical(stepped$evaluations, stepped$iterations + 1L)

  # THE CUBIC REMAINDER IS WHY THIS IS NOT EXACT. The quadratic model ignores the third-order term,
  # which over a step of about 0.11 in beta is a few per cent -- and it is optimistic, since the step
  # raises mortality and the exponential outruns the quadratic. Measured: 7.89 realised against 8.21
  # predicted. Dropping the half would predict 16.4, which no remainder explains.
  realised <- stepped$log_likelihood - before
  expect_equal(realised, stepped$predicted_gain, tolerance = 0.05)
})

test_that("the Armijo constant is actually used", {
  # At 1e-4 the condition is nearly inert, which is the point -- it rejects only a step that has
  # overshot. Wound up close to its ceiling of one half it becomes demanding, and a start below the
  # optimum then needs more halvings. Were the test written as "any improvement at all", the
  # constant would make no difference whatever and these two would agree.
  slack <- fit_run(start = -5, armijo = 1e-4)
  strict <- fit_run(start = -5, armijo = 0.49)

  expect_identical(slack$status, "converged")
  expect_identical(strict$status, "converged")
  expect_gt(strict$evaluations, slack$evaluations)
})

test_that("Omega and Z land where the maths says and nowhere else", {
  # A TIGHT TOLERANCE, so all three stop at effectively the same beta. Omega and Z reach the
  # CONVERGENCE TEST as well as the reported quantities, so at the ordinary tolerance the three runs
  # stop at slightly different points and the ratios below carry that wobble rather than the scaling
  # being tested.
  plain <- fit_run(tolerance = 1e-12)
  dispersed <- fit_run(tolerance = 1e-12, overdispersion = 2)
  scaled <- fit_run(tolerance = 1e-12, z_scale = 2)

  # Var(beta_hat) = Omega A^-1 B A^-1, so overdispersion scales it and Z does not.
  expect_equal(dispersed$variance, plain$variance * 2, tolerance = 1e-8)
  expect_equal(scaled$variance, plain$variance, tolerance = 1e-8)

  # p = tr(B A^-1) / Z, so Z scales it and Omega cancels out of it entirely.
  expect_equal(dispersed$penalty, plain$penalty, tolerance = 1e-8)
  expect_equal(scaled$penalty, plain$penalty / 2, tolerance = 1e-8)

  # L = (Aw log mu - Ew) / (Omega Z), so both scale it.
  expect_equal(dispersed$log_likelihood, plain$log_likelihood / 2, tolerance = 1e-8)
  expect_equal(scaled$log_likelihood, plain$log_likelihood / 2, tolerance = 1e-8)
})

test_that("Omega and Z reach the convergence test", {
  # `predicted_gain` is on the L scale, so dividing L by a large Z makes the same step look
  # insignificant and the fit stops sooner. Were the scale left out of that one expression, Z would
  # change what is reported and nothing about when the loop ends.
  tight <- fit_run(start = -5, z_scale = 1)
  loose <- fit_run(start = -5, z_scale = 1e6)

  expect_identical(tight$status, "converged")
  expect_identical(loose$status, "converged")
  expect_lt(loose$iterations, tight$iterations)
})

# THE DISJOINTNESS PRE-FLIGHT WALK.
#
# `disjoint()` is an assertion the user makes and R cannot check. Omitting the off-diagonal
# integrals is a COMPILE-TIME decision, so the claim has to be verified against the data BEFORE the
# fit block is built -- which is why this is the one place a second block and a second crossing are
# unavoidable.
#
# THE ARITHMETIC IS UNCHANGED, and that is the strongest thing to assert: the omitted integrals are
# identically zero, so a verified assertion must give BIT-IDENTICAL answers by a shorter route.
# Nothing numeric can therefore witness the omission, and `off_diagonals_omitted` and `output_count`
# exist for exactly that reason.

# Three groups by date of birth -- the fixture's births are 1940, 1945 and 1950, one record each.
birth_groups <- list(
  it_ast(~ .b < datey::datey(1943)),
  it_ast(~ .b >= datey::datey(1943) & .b < datey::datey(1948)),
  it_ast(~ .b >= datey::datey(1948))
)

test_that("a verified assertion drops the off-diagonal outputs and nothing else", {
  full <- fit_run(terms = birth_groups, disjoint = FALSE)
  lean <- fit_run(terms = birth_groups, disjoint = TRUE)

  expect_false(full$off_diagonals_omitted)
  expect_true(lean$off_diagonals_omitted)

  # 2 + 2k + 2 * triangle against 2 + 2k + 2k.
  expect_identical(full$output_count, 2L + 2L * 3L + 2L * 6L)
  expect_identical(lean$output_count, 2L + 2L * 3L + 2L * 3L)
})

test_that("a verified assertion gives bit-identical answers", {
  full <- fit_run(terms = birth_groups, disjoint = FALSE)
  lean <- fit_run(terms = birth_groups, disjoint = TRUE)

  expect_identical(lean$status, full$status)
  expect_identical(lean$beta, full$beta)
  expect_identical(lean$log_likelihood, full$log_likelihood)
  expect_identical(lean$penalty, full$penalty)
  expect_identical(lean$iterations, full$iterations)

  # The variance carries the off-diagonals back as EXACT zeros, which is what they are.
  expect_identical(lean$variance, full$variance)
  expect_equal(lean$variance[c(2, 3, 5)], rep(0, 3))
})

test_that("a false assertion is refused, naming the pair and the overlap", {
  # `.i$male` is TRUE on records 1 and 3; `.i$amount > 500` on records 1 and 2. They overlap on
  # record 1, so the claim is false and the block must not be compiled.
  overlapping <- list(it_ast(~ .i$male), it_ast(~ .i$amount > 500))

  expect_error(fit_run(terms = overlapping, disjoint = TRUE),
               "terms 1 and 2 overlap")
  expect_error(fit_run(terms = overlapping, disjoint = TRUE),
               "years of exposure")

  # And without the claim it is an ordinary, correct fit.
  expect_identical(fit_run(terms = overlapping, disjoint = FALSE)$status, "converged")
})

test_that("the check names the pair it found, not the first pair", {
  # Terms 1 and 2 are exclusive; 1 and 3 are not. The report must reach the third pair.
  mixed <- list(it_ast(~ .i$male), it_ast(~ !.i$male), it_ast(~ .i$amount > 500))
  expect_error(fit_run(terms = mixed, disjoint = TRUE), "terms 1 and 3 overlap")
})

test_that("the check only has to hold where the fit integrates", {
  # The same overlapping pair, but an include that keeps only record 2, where they do not overlap.
  # An individual excluded contributes zero to every off-diagonal integral whatever groups they are
  # in, so refusing this fit would refuse a sound one.
  overlapping <- list(it_ast(~ .i$male), it_ast(~ .i$amount > 500))

  expect_error(fit_run(terms = overlapping, disjoint = TRUE), "overlap")

  # RETURNING AT ALL IS THE ASSERTION. The check refuses by throwing, so a result -- of any status --
  # is proof that it accepted. What the fit then makes of a population this small is a separate
  # question and not what this test is about.
  narrowed <- fit_run(terms = overlapping, include = include(.i$amount > 2000), disjoint = TRUE)
  expect_true(narrowed$off_diagonals_omitted)
  expect_identical(narrowed$term_count, 2L)
})

test_that("a NaN convicts rather than acquitting", {
  # Every comparison against NaN is false, so a test written `total != 0` would let one THROUGH --
  # and an assertion that cannot be evaluated has not been verified. Same shape and same reason as
  # the Cholesky pivot test.
  not_a_number <- list(it_ast(~ .i$male), it_ast(~ log(0 - .i$amount)))
  expect_error(suppressWarnings(fit_run(terms = not_a_number, disjoint = TRUE)),
               "terms 1 and 2 overlap")
})

test_that("fewer than two terms has nothing to check and is not refused", {
  # No pairs, so no block to run at all -- and a single term is trivially exclusive of nothing.
  one <- fit_run(terms = list(it_ast(~ 1)), disjoint = TRUE)
  expect_identical(one$status, "converged")
  expect_true(one$off_diagonals_omitted)
  expect_identical(one$output_count, 2L + 2L + 2L)

  none <- fit_run(terms = list(), disjoint = TRUE)
  expect_identical(none$output_count, 2L)
})

# The two witnesses below exist because the disable-and-recheck found their guards unwitnessed: a
# suite of indicator terms can see neither the absolute value nor the coercion, since indicators are
# already non-negative and already coerce cleanly.

test_that("a signed term cannot cancel its way past the check", {
  # Exposures are 3, 2 and 1 years. These two terms have products 2, -3 and 0, so the raw
  # exposure-weighted sum is 2*3 - 3*2 = 0 EXACTLY -- a false assertion that would sail through on a
  # sum of products. Taking the magnitude first gives 12 and convicts it.
  #
  # Non-negativity is precisely what nobody has proved about a term: it is the same unknown that
  # makes this walk necessary at all.
  signed_cols <- modifyList(cols, list(left = c(2, -3, 0), right = c(1, 1, 0)))
  signed_terms <- list(it_ast(~ .i$left), it_ast(~ .i$right))

  expect_error(fit_run(terms = signed_terms, columns = signed_cols, disjoint = TRUE),
               "terms 1 and 2 overlap")
})

test_that("a duration-valued term is checked rather than refused", {
  # `durationy * durationy` is a product datey does not define and veil rightly refuses, so the
  # check has to coerce every term exactly as the fit recipe does. These two never overlap, so the
  # assertion holds -- and without the coercion it would not get as far as saying so.
  duration_cols <- modifyList(cols, list(
    service = datey::durationy(c(0, 5, 0)),
    leave   = datey::durationy(c(3, 0, 0))))
  duration_terms <- list(it_ast(~ .i$service), it_ast(~ .i$leave))

  checked <- fit_run(terms = duration_terms, columns = duration_cols, disjoint = TRUE)
  expect_true(checked$off_diagonals_omitted)
  expect_identical(checked$term_count, 2L)
})
