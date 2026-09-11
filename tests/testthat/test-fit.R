# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `fit()` is where the model syntax and the Newton-Raphson loop meet. The loop
# itself is tested in `test-veil_fit.R` against the engine entry point; what is
# tested here is everything R does around it -- the model, the Z scale, the
# result, the fitted mortality and the failure messages.
#
# THE ORACLES ARE ANALYTIC AND COME FROM `aev()`, which is a different recipe:
#
#   * one constant covariate       beta_hat = log(Aw / Ew) on the reference
#   * disjoint indicators          the same, group by group, with a DIAGONAL
#                                  variance of 1/A_j -- and at three terms that
#                                  is what witnesses the packed-triangle
#                                  convention, which two terms cannot see
#   * the fitted mortality         an A/E on it must come out at exactly one,
#                                  because the score is zero at the optimum

reference <- mortality_const(log_mu = -4)
basis <- settings(overdispersion = 1)

fit_data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
    pension   = c(5000, 12000, 30000, 8000, 15000),
    male      = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    E2R_start = datey::datey(rep(2015, 5)),
    E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
    E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

# Three groups, each with a death, so every coefficient is finite.
group_data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948, 1952, 1938, 1960, 1943)),
    group     = c("a", "a", "a", "b", "b", "b", "c", "c", "c"),
    pension   = c(5000, 12000, 30000, 8000, 15000, 22000, 4000, 9000, 40000),
    E2R_start = datey::datey(rep(2015, 9)),
    E2R_end   = datey::datey(c(2020, 2018, 2020, 2019, 2020, 2020, 2020, 2017, 2020)),
    E2R_died  = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

intercept_model <- model(ref_mortality = reference, covariates = covariates(level = 1))

# The bound the stopping rule promises.
#
# The loop stops when the gain still available in L falls below `L_tolerance`, and
# that gain is `lambda^2 / (2 Omega Z)`. Converting to a distance in beta: the
# loss in L* from being off by d is `d^2 I* / 2` with `I* = Omega^-1 Ew`, so
#
#     d <= sqrt(2 L_tolerance Omega^2 Z / Ew)
#
# and `Ew` at the optimum is `Aw`. SO OMEGA AND Z DO MOVE WHERE THE LOOP STOPS,
# even though neither moves the maximum -- which is why two fits differing only
# in Z are NOT equal to the bit, and each is checked against the analytic answer
# within its own bound instead.
#
# NOT IN UNITS OF THE REPORTED STANDARD ERROR, which is a different quantity
# once the weight is not an indicator.
beta_bound <- function(A, L_tolerance = 1e-6, overdispersion = 1, Z = 1) {
  sqrt(2 * L_tolerance * overdispersion^2 * Z / A)
}

# ---- model() ---------------------------------------------------------------

test_that("model() holds a reference mortality and its covariates", {
  m <- model(ref_mortality = reference, covariates = covariates(a = 1, b = .x))
  expect_true(is_model(m))
  expect_identical(length(m$covariates), 2L)
  expect_identical(m$ref_mortality$kind, "obj")
})

test_that("model() takes a pronoun expression as its reference", {
  m <- model(ref_mortality = ~ -10.5 + 0.09 * .x, covariates = covariates(level = 1))
  expect_true(is_model(m))
  # Anything held fixed lives in the reference, which is already an offset, so
  # fixed effects need no feature of their own.
  expect_true(it_uses_t(m$ref_mortality))
})

test_that("model() refuses what it cannot fit", {
  expect_error(model(covariates = covariates(level = 1)), "`ref_mortality` is required")
  expect_error(model(ref_mortality = reference), "`covariates` is required")
  expect_error(model(ref_mortality = reference, covariates = variable(.x)),
               "must be a `covariates` object")
  expect_error(model(ref_mortality = reference, covariates = covariates()),
               "at least one covariate")
  expect_error(model(ref_mortality = age(65, 95), covariates = covariates(level = 1)),
               "only reference `mortality` objects")
})

# ---- the analytic oracle ---------------------------------------------------

test_that("one constant covariate finds log(A/E) on the reference", {
  # log mu = log mu_ref + beta, so L is maximised where exp(beta) Ew = Aw and
  # beta_hat = log(Aw / Ew) taken at beta = 0 -- which is exactly the A and the
  # E of an A/E on the reference mortality. A different recipe supplies it.
  a <- aev(fit_data, mortality = reference, settings = basis)
  f <- fit(fit_data, intercept_model, settings = basis)

  expect_true(is_fit(f))
  expect_lt(abs(unname(f$beta) - log(a$A / a$E)), beta_bound(a$A))
  expect_lte(f$predicted_gain, 1e-6)
})

test_that("the fitted mortality is a mortality that the data is neutral against", {
  # THE SCORE IS ZERO AT THE OPTIMUM, so `Ew` on the fitted mortality equals
  # `Aw`: an A/E on it must come out at one. That checks the coefficient, the
  # rebuilt AST and the mortality wrapper in a single statement, and it goes
  # back through the engine rather than through R arithmetic.
  f <- fit(fit_data, intercept_model, settings = basis)

  expect_true(is_mortality(f$mortality))
  fitted <- f$mortality
  a <- aev(fit_data, mortality = fitted, settings = basis)
  expect_equal(a$A / a$E, 1, tolerance = 1e-5)

  # And it says the same thing as the arithmetic, read directly.
  expect_equal(
    log_mu(fitted, list(birth = datey::datey(1950)), datey::datey(2018)),
    -4 + unname(f$beta)
  )
})

test_that("a fitted model serves as the reference of the next one", {
  # Layered fitting comes for nothing once the coefficients are numbers. Fitting
  # the same covariate again against the fitted mortality has nothing left to
  # find, so the second coefficient must be zero.
  first <- fit(fit_data, intercept_model, settings = basis)
  again <- model(ref_mortality = first$mortality, covariates = covariates(level = 1))
  second <- fit(fit_data, again, settings = basis)

  a <- aev(fit_data, mortality = reference, settings = basis)
  expect_lt(abs(unname(second$beta)), 2 * beta_bound(a$A))
})

test_that("the variance and the penalty are what the A/E says they must be", {
  # For the intercept-only model both collapse onto the reference A/E:
  # Var(beta_hat) = Omega Ew^2 / Aw^2 = V / (A E), because Ew = Aw at the
  # maximum and V already carries Omega. With a lives weight that is 1 / A.
  a <- aev(fit_data, mortality = reference, settings = basis)
  f <- fit(fit_data, intercept_model, settings = basis)

  expect_equal(unname(f$variance[1, 1]), 1 / a$A, tolerance = 2 * beta_bound(a$A))

  # `p = k` EXACTLY where `w^2 = w`, so for lives, and then L_P is -AIC/2.
  expect_equal(f$penalty, 1)
  expect_equal(f$penalised_log_likelihood, f$log_likelihood - 1)
})

test_that("a time-varying covariate fits, and both score equations come out zero", {
  # AGE IS THE COMMONEST COVARIATE THERE IS -- Gompertz is `log mu = a + b x` --
  # and `.x` is a `durationy`, so this is the path where `X_j X_l` would ask for
  # `durationy * durationy` if the recipe did not coerce every term once.
  #
  # THE ORACLE IS BOTH SCORE EQUATIONS AT ONCE, and each is an A/E on the FITTED
  # mortality with a different weight. At the optimum every term's score is
  # zero, `Aw X_j = Ew X_j`, and an A/E weighted by `X_j` is exactly that ratio.
  # So the level term says an unweighted A/E is one, and the age term says an
  # age-weighted A/E is one. Nothing here recomputes what the fit computed.
  two_terms <- model(ref_mortality = reference, covariates = covariates(level = 1, age = .x))
  f <- fit(fit_data, two_terms, settings = basis)

  expect_identical(length(f$beta), 2L)
  expect_true(all(is.finite(f$beta)))

  fitted <- f$mortality
  level_score <- aev(fit_data, mortality = fitted, settings = basis)
  # `.x * 1` RATHER THAN `.x`, and the reason is worth stating: an A/E also forms
  # `V = Omega Ew^2`, and `durationy * durationy` is undefined, so a bare
  # duration cannot be an A/E weight. `fit()` itself has no such trouble -- its
  # recipe coerces every term through `ToDouble` once, which is exactly why age
  # works as a covariate.
  age_score <- aev(fit_data, mortality = fitted, weight = .x * 1, settings = basis)

  expect_equal(level_score$A / level_score$E, 1, tolerance = 1e-4)
  expect_equal(age_score$A / age_score$E, 1, tolerance = 1e-4)

  # Two terms, so the variance has a genuine off-diagonal: age and level are
  # informative about each other, unlike the disjoint indicators below.
  expect_identical(dim(f$variance), c(2L, 2L))
  expect_false(f$variance[1, 2] == 0)
  expect_identical(f$variance[1, 2], f$variance[2, 1])
})

# ---- three terms, which is the smallest case that sees the packing ---------

test_that("disjoint indicators fit group by group, with a diagonal variance", {
  # Each group is its own intercept-only problem, so beta_j = log(A_j / E_j) and
  # the terms do not inform each other at all: the information is diagonal and
  # so is the variance, at 1 / A_j for a lives weight.
  #
  # THREE TERMS IS THE POINT. At two, the row-major and column-major upper
  # packings are the same sequence, so a transposed unpacking cannot be seen.
  terms <- covariates(one = .i$group == "a", two = .i$group == "b", three = .i$group == "c")
  f <- fit(group_data, model(ref_mortality = reference, covariates = terms),
           settings = basis)

  groups <- c("a", "b", "c")
  per_group <- lapply(1:3, function(g) {
    aev(group_data, mortality = reference,
        include = include(.i$group == groups[[g]]), settings = basis)
  })

  for (g in 1:3) {
    a <- per_group[[g]]
    expect_lt(abs(unname(f$beta[[g]]) - log(a$A / a$E)), beta_bound(a$A))

    # The variance is read AT beta_hat, so it carries the same first-order error
    # the stopping rule allows in beta.
    expect_equal(unname(f$variance[g, g]), 1 / a$A, tolerance = 2 * beta_bound(a$A))
  }

  off <- f$variance[upper.tri(f$variance)]
  expect_equal(off, rep(0, 3L))
  expect_equal(f$penalty, 3)
  expect_identical(dimnames(f$variance), list(c("one", "two", "three"),
                                              c("one", "two", "three")))
})

test_that("disjoint() gives a BIT-IDENTICAL fit, by a shorter route", {
  # `disjoint()` asserts something extra about the same terms, and the engine
  # verifies it against the data before compiling. What it then omits are
  # integrals that are identically zero, so the answer cannot move by so much as
  # a bit -- which is the strongest statement available and the reason the
  # structural witnesses in `test-veil_fit.R` exist at all.
  plain <- covariates(one = .i$group == "a", two = .i$group == "b", three = .i$group == "c")
  asserted <- disjoint(.i$group == "a", .i$group == "b", .i$group == "c")

  a <- fit(group_data, model(ref_mortality = reference, covariates = plain), settings = basis)
  b <- fit(group_data, model(ref_mortality = reference, covariates = asserted), settings = basis)

  expect_identical(unname(a$beta), unname(b$beta))
  expect_identical(a$log_likelihood, b$log_likelihood)
  expect_identical(a$penalty, b$penalty)
  expect_identical(unname(a$variance), unname(b$variance))
})

test_that("a false disjoint() assertion stops the fit", {
  # R can only see that a term does not use `.t`; whether two terms overlap is a
  # property of the data, so the engine walks it once before compiling. A wrong
  # claim is an error rather than a quiet fallback to the full triangle: the
  # right answer computed anyway would leave the user believing something about
  # their data that is not true.
  overlapping <- disjoint(.i$group == "a", .i$group != "b")
  m <- model(ref_mortality = reference, covariates = overlapping)

  expect_error(fit(group_data, m, settings = basis), "no individual is in more than one")
  expect_error(fit(group_data, m, settings = basis), "terms 1 and 2 overlap")

  # The same terms without the claim are an ordinary, correct fit.
  honest <- covariates(.i$group == "a", .i$group != "b")
  expect_true(is_fit(fit(group_data, model(ref_mortality = reference, covariates = honest),
                         settings = basis)))
})

test_that("disjointness survives multiplication by a shape, and is still verified", {
  # `I_j I_l = 0` implies `I_j phi * I_l phi = 0`, so the assertion carries
  # through `* shape` -- and the engine checks the terms it is actually given.
  shaped <- disjoint(.i$group == "a", .i$group == "b", .i$group == "c") * variable(.x)
  expect_true(is_disjoint(shaped))

  f <- fit(group_data, model(ref_mortality = reference, covariates = shaped), settings = basis)
  expect_true(is_fit(f))

  plain <- covariates(.i$group == "a", .i$group == "b", .i$group == "c") * variable(.x)
  g <- fit(group_data, model(ref_mortality = reference, covariates = plain), settings = basis)
  expect_identical(unname(f$beta), unname(g$beta))
})

# ---- the Z scale -----------------------------------------------------------

test_that("the Z scale is one where `w^2 = w` can be shown", {
  # No weight, the literal 1, and an expression an operator has forced to be
  # logical. All three give `Z = Ew^2/Ew = 1` exactly.
  expect_identical(fit(fit_data, intercept_model, settings = basis)$Z, 1)
  expect_identical(fit(fit_data, intercept_model, weight = 1, settings = basis)$Z, 1)
  expect_identical(
    fit(fit_data, intercept_model, weight = .i$pension > 10000, settings = basis)$Z,
    1
  )

  # A bare logical column is an indicator in fact and unprovable from the
  # expression, because a column's type is unknown until the data arrives. So it
  # is measured rather than assumed -- and the measurement gives 1 anyway.
  expect_equal(fit(fit_data, intercept_model, weight = .i$male, settings = basis)$Z, 1)
})

test_that("a fit makes exactly one pass outside the iteration, whatever the weight", {
  # ONE `cpp_veil_run()`, always: the header's population numbers and the
  # test-mortality A/E are two blocks of it. Z is then read from that run rather
  # than costing a run of its own, so the weighted case pays what it always did
  # and the unweighted case pays one pass for a header.
  count <- function(...) {
    runs <- 0L
    original <- cpp_veil_run
    local_mocked_bindings(
      cpp_veil_run = function(...) {
        runs <<- runs + 1L
        original(...)
      }
    )
    fit(fit_data, intercept_model, settings = basis, ...)
    runs
  }

  expect_identical(count(), 1L)
  expect_identical(count(weight = 1), 1L)
  expect_identical(count(weight = .i$pension), 1L)
  expect_identical(count(Z = 3, weight = .i$pension), 1L)
})

test_that("the engine makes `Ew^2 = Ew` exactly where `w^2 = w`", {
  # THIS IS WHAT THE `weight_squares_to_itself()` SHORTCUT NOW RESTS ON, and it is
  # the engine's property rather than R's: `Ew` and `Ew^2` accumulate identical
  # per-record values in the same order, so their ratio is 1 to the bit and not
  # merely to a tolerance. The shortcut used to earn its place by avoiding a pass;
  # now that the diagnostics run regardless, the only thing left for it to protect
  # is exactness -- and IF THIS TEST EVER FAILS the shortcut is load-bearing again.
  columns <- exp_data_columns(fit_data)
  clicks <- time_scale_clicks(1 / 4)
  test_ast <- it_obj(default_mortality())

  measured <- function(weight_ast) {
    d <- fit_diagnostics(test_ast, weight_ast, NULL, columns, clicks, 1L)
    d$Ew2 / d$Ew
  }

  expect_identical(measured(NULL), 1)
  expect_identical(measured(it_capture(quote(1), environment())), 1)
  expect_identical(measured(it_capture(quote(.i$pension > 10000), environment())), 1)

  # The control: with an amount weight it is nowhere near 1, so the comparison
  # above is not vacuous.
  expect_gt(measured(it_capture(quote(.i$pension), environment())), 100)
})

test_that("an amount weight uses the default mortality rather than refusing", {
  # It used to be an error. The default is fixed and shared, which is what Z
  # needs -- one yardstick beats a good one.
  f <- fit(fit_data, intercept_model, weight = .i$pension, settings = basis)
  named <- fit(fit_data, intercept_model, weight = .i$pension,
               test_mortality = default_mortality(), settings = basis)

  # THE SIGNATURE DEFAULT IS NEVER EVALUATED, so the default is named twice --
  # once for `?fit` and once in the code. This is what stops them drifting.
  expect_identical(f, named)
  expect_false(isTRUE(all.equal(f$Z, 1)))
})

test_that("an explicit Z wins even where one would have been provable", {
  # THE ORDER OF THE TWO CHECKS. `Z` is read before the `w^2 = w` shortcut,
  # so a caller forcing a common scale across a set of runs gets the number they
  # asked for rather than the 1 this particular run could have proved.
  f <- fit(fit_data, intercept_model, Z = 4, settings = basis)
  expect_identical(f$Z, 4)
  expect_equal(f$penalty, 1 / 4)
})

test_that("a test mortality and an explicit Z are alternatives, not a pair", {
  expect_error(
    fit(fit_data, intercept_model, weight = .i$pension,
        test_mortality = gompertz_mortality(), Z = 2, settings = basis),
    "Give either `test_mortality` or `Z`, not both"
  )
})

test_that("Z from a test mortality is V/E on a run pinned at overdispersion one", {
  # THE PIN MATTERS. The engine returns `V = Omega Ew^2` for the Omega of the
  # run that produced it, so only at Omega = 1 is `V/E` equal to `Ew^2/Ew`.
  # Using the fit's own Omega would scale L and p together -- invisible in any
  # ranking, and wrong in the calibration that says a change of one is
  # significant.
  test_mortality <- gompertz_mortality()
  scale <- aev(fit_data, mortality = test_mortality, weight = .i$pension,
               settings = settings(overdispersion = 1))
  expected <- scale$V / scale$E

  f <- fit(fit_data, intercept_model, weight = .i$pension, test_mortality = test_mortality,
           settings = settings(overdispersion = 3))

  expect_equal(f$Z, expected)

  # Not the same as reading it off a run at the fit's own overdispersion.
  wrong <- aev(fit_data, mortality = test_mortality, weight = .i$pension,
               settings = settings(overdispersion = 3))
  expect_false(isTRUE(all.equal(f$Z, wrong$V / wrong$E)))
})

test_that("Z as a number is taken as given", {
  f <- fit(fit_data, intercept_model, weight = .i$pension, Z = 2500, settings = basis)
  expect_identical(f$Z, 2500)
})

test_that("Z rescales L and p but never beta", {
  # L = Z^-1 L*, and p = Z^-1 tr(J I^-1), so both scale together and the
  # argument of the maximum does not move at all.
  one <- fit(fit_data, intercept_model, weight = .i$pension, Z = 1000, settings = basis)
  two <- fit(fit_data, intercept_model, weight = .i$pension, Z = 2000, settings = basis)

  # Z scales the convergence test as well, so the two stop at slightly different
  # points. Each is checked against the analytic maximum within its own bound.
  weighted <- aev(fit_data, mortality = reference, weight = .i$pension, settings = basis)
  target <- log(weighted$A / weighted$E)
  expect_lt(abs(unname(one$beta) - target), beta_bound(weighted$A, Z = 1000))
  expect_lt(abs(unname(two$beta) - target), beta_bound(weighted$A, Z = 2000))

  expect_equal(two$log_likelihood, one$log_likelihood / 2, tolerance = 1e-6)
  expect_equal(two$penalty, one$penalty / 2, tolerance = 1e-6)
})

test_that("a Z that is not a usable scale is refused", {
  # `Z` is an ordinary value rather than a pronoun expression, so all three
  # are now caught by the same check rather than the last one falling to the
  # parser.
  for (bad in list(0, -1, c(1, 2), NA_real_, Inf, "1")) {
    expect_error(fit(fit_data, intercept_model, Z = bad, settings = basis),
                 "single positive finite number")
  }
})

# ---- overdispersion --------------------------------------------------------

test_that("overdispersion scales L and the variance but never beta", {
  # L* = Omega^-1 (Aw log mu - Ew), so Omega does not move the maximum; it
  # scales the likelihood, and it scales the sandwich variance once.
  one <- fit(fit_data, intercept_model, settings = settings(overdispersion = 1))
  two <- fit(fit_data, intercept_model, settings = settings(overdispersion = 2))

  # Omega does not move the maximum, but it does move where the loop stops:
  # the convergence test divides the gain still available by `Omega Z`, so a
  # larger Omega stops sooner. Each is checked against the analytic answer
  # within its own bound rather than against the other to the bit.
  a <- aev(fit_data, mortality = reference, settings = basis)
  target <- log(a$A / a$E)
  expect_lt(abs(unname(one$beta) - target), beta_bound(a$A))
  expect_lt(abs(unname(two$beta) - target), beta_bound(a$A, overdispersion = 2))
  expect_equal(two$log_likelihood, one$log_likelihood / 2, tolerance = 1e-5)
  expect_equal(two$variance[1, 1], one$variance[1, 1] * 2, tolerance = 1e-5)

  # Omega cancels in `J I^-1`, so the penalty does not move.
  expect_equal(two$penalty, one$penalty, tolerance = 1e-5)
  expect_identical(two$overdispersion, 2)
})

test_that("overdispersion is required, here as everywhere", {
  expect_error(fit(fit_data, intercept_model), "`overdispersion` is required")
})

# ---- the iteration ---------------------------------------------------------

test_that("the answer does not depend on where beta started", {
  # The objective is strictly concave for a non-negative weight, so there is one
  # maximum and every start must reach it.
  a <- aev(fit_data, mortality = reference, settings = basis)
  target <- log(a$A / a$E)

  for (start in list(0, 3, -3, 6)) {
    f <- fit(fit_data, intercept_model, start = start, settings = basis)
    expect_lt(abs(unname(f$beta) - target), beta_bound(a$A))
  }
})

test_that("a start of the wrong shape is refused", {
  expect_error(fit(fit_data, intercept_model, start = c(0, 0), settings = basis),
               "must be 1 finite number, one for each model term")
  expect_error(fit(fit_data, intercept_model, start = NA_real_, settings = basis),
               "must be 1 finite number")
})

test_that("a tighter tolerance costs iterations and buys accuracy", {
  loose <- fit(fit_data, intercept_model, L_tolerance = 1e-2, settings = basis)
  tight <- fit(fit_data, intercept_model, L_tolerance = 1e-10, settings = basis)

  expect_lte(loose$iterations, tight$iterations)
  expect_lt(tight$predicted_gain, loose$predicted_gain)
  expect_lte(tight$predicted_gain, 1e-10)
})

test_that("the numerical-method arguments are passed through and checked", {
  # The ceiling on armijo is one half: that is exactly what a full Newton step
  # delivers on a quadratic, so anything at or above it rejects the full step
  # near the optimum and the iteration cannot finish.
  expect_error(fit(fit_data, intercept_model, armijo = 0.5, settings = basis),
               "strictly between 0 and 0.5")
  expect_error(fit(fit_data, intercept_model, armijo = 0, settings = basis),
               "strictly between 0 and 0.5")
  expect_error(fit(fit_data, intercept_model, max_iterations = 0, settings = basis),
               "at least one")
})

test_that("a fit that runs out of iterations fails rather than returning", {
  # A fit returns a fit or it fails: not converged is also not a fit.
  expect_error(
    fit(fit_data, intercept_model, start = 6, max_iterations = 2, settings = basis),
    "did not converge in 2 iterations"
  )
})

# ---- failure ---------------------------------------------------------------

test_that("a collinear model fails, naming the term and what it duplicates", {
  # THE PIVOT NAMES WHERE THE DEPENDENCY WAS DETECTED, NOT THE CULPRIT: with
  # terms 1 and 3 identical it is term 3 that cannot be added, because 1 and 2
  # are still independent between them. The message must therefore say what
  # term 3 is a combination of, never that term 3 is the bad one.
  terms <- covariates(level = 1, male = .i$male, again = 1)
  m <- model(ref_mortality = reference, covariates = terms)

  expect_error(fit(fit_data, m, settings = basis), "not identifiable")
  expect_error(fit(fit_data, m, settings = basis), "Term 3 \\(again\\)")
  expect_error(fit(fit_data, m, settings = basis), "linear combination of the terms before it")

  # And it names the actual dependency: `again` is one times `level`.
  expect_error(fit(fit_data, m, settings = basis), "\\[1\\] level")
})

test_that("include narrows the population the fit sees", {
  restricted <- fit(group_data, intercept_model,
                    include = include(.i$group == "a"), settings = basis)
  a <- aev(group_data, mortality = reference,
           include = include(.i$group == "a"), settings = basis)

  expect_lt(abs(unname(restricted$beta) - log(a$A / a$E)), beta_bound(a$A))

  whole <- fit(group_data, intercept_model, settings = basis)
  expect_false(isTRUE(all.equal(unname(restricted$beta), unname(whole$beta))))
})

test_that("fit() refuses a breakdown-shaped include", {
  expect_error(
    fit(group_data, intercept_model,
        include = includes(include(.i$group == "a"), include(.i$group == "b")),
        settings = basis),
    "takes a single `include`"
  )
})

# ---- the result ------------------------------------------------------------

test_that("coef() and vcov() give the fitted coefficients and their variance", {
  terms <- covariates(one = .i$group == "a", two = .i$group == "b", three = .i$group == "c")
  f <- fit(group_data, model(ref_mortality = reference, covariates = terms), settings = basis)

  expect_identical(coef(f), f$beta)
  expect_identical(names(coef(f)), c("one", "two", "three"))
  expect_identical(vcov(f), f$variance)
  expect_identical(dim(vcov(f)), c(3L, 3L))
  expect_equal(vcov(f), t(vcov(f)))
})

test_that("an unnamed covariate is labelled by the expression the user wrote", {
  # Naming covariates is not usual practice, so the label falls back to what
  # they actually wrote, which identifies the term better than `X1` would.
  f <- fit(fit_data, model(ref_mortality = reference, covariates = covariates(.i$male)),
           settings = basis)
  expect_identical(names(coef(f)), ".i$male")
})

test_that("a fit prints its estimates, its scale and its cost", {
  f <- fit(fit_data, intercept_model, settings = basis)
  printed <- paste(capture.output(print(f)), collapse = "\n")

  expect_match(printed, "<fit: 1 term")
  expect_match(printed, "estimate")
  expect_match(printed, "std_error")
  # A lone fit shows the three likelihood numbers unshifted, with `k`.
  expect_match(printed, "L -[0-9.]+   p [0-9.]+   L_P -[0-9.]+   k 1")
  # The Z scale is reported because it is the unit the numbers are quoted in.
  expect_match(printed, "Z 1")
  # And the header, which is shared with a comparison.
  expect_match(printed, "deaths, ")
  expect_match(printed, "test mortality  default_mortality")
})

test_that("a model prints its reference and its terms", {
  printed <- paste(capture.output(print(intercept_model)), collapse = "\n")
  expect_match(printed, "<model: 1 term>")
  expect_match(printed, "reference: <mortality_const>")
  expect_match(printed, "level")
})
