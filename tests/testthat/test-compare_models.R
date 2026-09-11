# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `compare_models()` computes no statistic of its own: every number in its table
# comes from a `fit()`. So the oracles are mostly identities against `fit()`
# itself, which is what pins the claim that a candidate in a comparison is fitted
# exactly as it would have been alone:
#
#   * against fit()        every field of every candidate, to the bit, when Z is
#                          given as a number
#   * p == terms           for an indicator weight, where w^2 = w, and then
#                          L_P = -AIC/2
#   * Omega                p is free of it and L scales as 1/Omega, which is what
#                          makes a larger Omega favour simpler models
#   * one Z                two candidates with different reference mortalities
#                          still get the same Z, because it comes from the test
#                          mortality and never from the candidate

reference <- mortality_const(log_mu = -4)
basis <- settings(overdispersion = 1)

# Three groups, each with a death, so a group split is identifiable.
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

level_model <- model(ref_mortality = reference, covariates = covariates(level = 1))
group_model <- model(
  ref_mortality = reference,
  covariates = covariates(a = .i$group == "a",
                          b = .i$group == "b",
                          c = .i$group == "c")
)

# Term 3 duplicates term 1, so the information matrix is singular and this
# candidate cannot be fitted at all.
collinear_model <- model(
  ref_mortality = reference,
  covariates = covariates(level = 1, group_a = .i$group == "a", again = 1)
)

two <- function(...) {
  compare_models(group_data, models = list(level = level_model, group = group_model),
                 settings = basis, ...)
}

# ---- the shape of the answer -----------------------------------------------

test_that("compare_models() returns a comparison over every candidate", {
  comparison <- two()

  expect_true(is_model_comparison(comparison))
  expect_identical(nrow(comparison$table), 2L)
  expect_identical(length(comparison$fits), 2L)
  expect_length(comparison$failures, 0L)
  expect_identical(comparison$Z, 1)
  expect_identical(comparison$overdispersion, 1)

  expect_identical(
    names(comparison$table),
    c("k", "log_likelihood", "penalty", "penalised_log_likelihood", "status")
  )
  expect_true(all(vapply(comparison$fits, is_fit, logical(1L))))
  expect_identical(sort(rownames(comparison$table)), c("group", "level"))
  expect_identical(rownames(comparison$table), names(comparison$fits))
})

test_that("the table ranks on penalised log-likelihood, best first", {
  comparison <- two()
  penalised <- comparison$table$penalised_log_likelihood

  # Descending, because larger is better -- this is the AIC without the -2.
  expect_true(all(diff(penalised) <= 0))
  expect_identical(penalised[[1L]], max(penalised))

  # THE TABLE KEEPS THE RAW VALUES. Only `print()` deducts the maximum, which is
  # why no argument is needed to get at them.
  expect_true(all(penalised < 0))
  expect_equal(penalised, comparison$table$log_likelihood - comparison$table$penalty)
})

test_that("the fits stay alongside the rows they rank", {
  # WORST CANDIDATE FIRST, DELIBERATELY. With the winner already in position one
  # the ranking is the identity and this test cannot see a reorder at all -- which
  # is what it did before the input order was turned round, and it passed with the
  # reordering of `fits` removed entirely.
  comparison <- compare_models(group_data,
                               models = list(group = group_model, level = level_model),
                               settings = basis)
  expect_identical(rownames(comparison$table), c("level", "group"))

  # THE FITS ARE REORDERED WITH THE TABLE. Reading `fits[[i]]` against
  # `table[i, ]` is the whole convenience of the object, so the two orders must
  # agree row by row rather than merely contain the same names.
  expect_identical(names(comparison$fits), rownames(comparison$table))
  for (at in seq_len(nrow(comparison$table))) {
    expect_identical(comparison$fits[[at]]$penalised_log_likelihood,
                     comparison$table$penalised_log_likelihood[[at]])
    expect_identical(length(comparison$fits[[at]]$beta),
                     comparison$table$k[[at]])
  }
})

test_that("best_fit() returns the winner", {
  comparison <- two()
  best <- best_fit(comparison)

  expect_true(is_fit(best))
  expect_identical(best$penalised_log_likelihood,
                   max(comparison$table$penalised_log_likelihood))
  expect_identical(best, comparison$fits[[1L]])
})

test_that("best_fit() refuses anything that is not a comparison", {
  expect_error(best_fit(fit(group_data, level_model, settings = basis)),
               "must be a `model_comparison`")
})

# ---- the same answer as fit() ----------------------------------------------

test_that("a candidate is fitted exactly as it would have been alone", {
  comparison <- two()
  alone <- list(level = fit(group_data, level_model, settings = basis),
                group = fit(group_data, group_model, settings = basis))

  # TO THE BIT, and every field. `compare_models()` computes nothing itself, so
  # anything that differs here is a difference in how the candidate was fitted.
  for (label in names(alone)) {
    expect_identical(comparison$fits[[label]], alone[[label]])
  }
})

test_that("the penalty is the term count for an indicator weight", {
  # w^2 = w makes Z one and p exactly dim(beta), so L_P is then -AIC/2 and a gap
  # of one is one parameter's worth.
  comparison <- compare_models(
    group_data,
    models = list(level = level_model, group = group_model),
    weight = .i$group != "z",
    settings = basis
  )

  expect_equal(comparison$table$penalty, as.double(comparison$table$k))
  expect_equal(comparison$table$penalised_log_likelihood,
               comparison$table$log_likelihood - comparison$table$k)
})

test_that("overdispersion scales the fit term and leaves the penalty alone", {
  # p = Z^-1 tr(J I^-1) with I and J both carrying Omega^-1, which cancels; L
  # carries it once. So a larger Omega shrinks the gain from a richer model
  # against an unchanged penalty, and simpler models win -- the QAIC behaviour,
  # and the reason overdispersion is the user's to supply.
  #
  # NOT TO THE BIT: Omega moves where the loop stops, though not where the
  # maximum is, so a tight tolerance and a tolerance-sized comparison.
  one <- two(L_tolerance = 1e-12)
  ten <- compare_models(group_data,
                        models = list(level = level_model, group = group_model),
                        overdispersion = 10, L_tolerance = 1e-12)

  expect_equal(ten$table$penalty, one$table$penalty)
  expect_equal(ten$table$log_likelihood, one$table$log_likelihood / 10)

  # And the consequence: the richer candidate's advantage in L is divided by ten
  # while the extra two parameters still cost two.
  advantage <- function(comparison) {
    table <- comparison$table
    table["group", "log_likelihood"] - table["level", "log_likelihood"]
  }
  expect_equal(advantage(ten), advantage(one) / 10)
})

# ---- one Z, shared -----------------------------------------------------------

test_that("Z from a test mortality is the same for every candidate", {
  # AND IT IS NEVER THE CANDIDATE'S OWN REFERENCE. These two candidates sit on
  # different reference mortalities; a Z read from the reference would differ
  # between them and destroy the comparability it exists to provide.
  other <- model(ref_mortality = mortality_const(log_mu = -2),
                 covariates = covariates(level = 1))
  comparison <- compare_models(group_data,
                               models = list(low = level_model, high = other),
                               weight = .i$pension,
                               test_mortality = mortality_const(log_mu = -3),
                               settings = basis)

  expected <- fit(group_data, level_model, weight = .i$pension,
                  test_mortality = mortality_const(log_mu = -3), settings = basis)$Z

  expect_equal(comparison$Z, expected)
  expect_equal(comparison$fits[["low"]]$Z, expected)
  expect_equal(comparison$fits[["high"]]$Z, expected)
  expect_false(isTRUE(all.equal(expected, 1)))

  # AND THE SHARED Z IS THE ONE THEY WERE FITTED AT. The `$Z` field is written by
  # the comparison rather than read back from the run, so it agrees with itself
  # however Z was obtained; only this identity against `fit()` can tell a shared Z
  # from one taken per candidate from its own reference mortality.
  expect_identical(comparison$fits[["low"]],
                   fit(group_data, level_model, weight = .i$pension,
                       test_mortality = mortality_const(log_mu = -3), settings = basis))
})

test_that("Z as a number reaches every candidate", {
  comparison <- compare_models(group_data,
                               models = list(level = level_model, group = group_model),
                               weight = .i$pension, Z = 2.5, settings = basis)

  expect_identical(comparison$Z, 2.5)
  expect_true(all(vapply(comparison$fits, function(x) identical(x$Z, 2.5), logical(1L))))

  # AND IT REACHES THE ENGINE, not merely the reported field. Z divides L and p
  # alike, so a Z that was recorded but not passed on would leave every number in
  # the table on the wrong scale while `$Z` still read 2.5. Only a comparison with
  # `fit()` at the same Z can see that, and it holds to the bit.
  expect_identical(comparison$fits[["level"]],
                   fit(group_data, level_model, weight = .i$pension, Z = 2.5,
                       settings = basis))
})

test_that("a weighted comparison falls back on the default mortality, once", {
  # It used to refuse. Now it measures Z on `default_mortality()`, and it does so
  # a SINGLE time however many candidates there are, which is the whole
  # mechanical reason this function exists.
  runs <- 0L
  original <- cpp_veil_run
  local_mocked_bindings(
    cpp_veil_run = function(...) {
      runs <<- runs + 1L
      original(...)
    }
  )

  comparison <- compare_models(group_data,
                               models = list(level = level_model, group = group_model,
                                             again = level_model),
                               weight = .i$pension, settings = basis)

  expect_identical(runs, 1L)
  expect_equal(comparison$Z,
               fit(group_data, level_model, weight = .i$pension, settings = basis)$Z)
  expect_false(isTRUE(all.equal(comparison$Z, 1)))
})

test_that("the diagnostic pass happens once, whatever the weight", {
  # ONE PASS, NOT ONE PER CANDIDATE, and not one per diagnostic either: the
  # population numbers and the test-mortality A/E are two blocks of a single
  # `cpp_veil_run()`, because the engine crosses the data once for the whole list.
  count <- function(...) {
    runs <- 0L
    original <- cpp_veil_run
    local_mocked_bindings(
      cpp_veil_run = function(...) {
        runs <<- runs + 1L
        original(...)
      }
    )
    compare_models(group_data,
                   models = list(level = level_model, group = group_model,
                                 again = level_model),
                   settings = basis, ...)
    runs
  }

  expect_identical(count(), 1L)
  expect_identical(count(weight = .i$pension), 1L)
})

test_that("the header reports the population it fitted on", {
  comparison <- two()
  d <- comparison$diagnostics

  # Nine lives, five of whom die, over an exposure this A/E measures directly:
  # `mu = 1` makes E the person-years and A the death count.
  expect_identical(d$deaths, 5)
  expect_equal(d$exposure, sum(c(5, 3, 5, 4, 5, 5, 5, 2, 5)))

  # The test-mortality block, and Z read from it.
  expect_equal(d$Ew2 / d$Ew, 1)
  expect_identical(comparison$Z, 1)

  shown <- capture.output(print(comparison))
  expect_true(any(grepl("5 deaths, 39 years of exposure", shown)))
  expect_true(any(grepl("test mortality  default_mortality\\(\\)", shown)))
  expect_true(any(grepl("weight          lives", shown)))
  expect_true(any(grepl("include         all records", shown)))
})

test_that("both printed likelihood columns shift by the same constant", {
  # THE ONE PROPERTY THE SHIFT EXISTS FOR: `L - L_P = p` still holds row by row,
  # so the arithmetic can be checked on the face of the table. Shifting only `L_P`,
  # or shifting `L` by its own maximum instead, both break it and both were SILENT
  # until this test existed.
  comparison <- two()
  shown <- comparison_shown_table(comparison)
  best <- max(comparison$table$penalised_log_likelihood)

  expect_equal(shown$L - shown$L_P, shown$p)
  expect_identical(shown$L_P[[1L]], 0)
  expect_equal(shown$L_P, comparison$table$penalised_log_likelihood - best)
  expect_equal(shown$L, comparison$table$log_likelihood - best)
  expect_identical(names(shown), c("L_P", "L", "p", "k"))

  # Not `L` shifted by its own maximum, which for these candidates differs.
  expect_false(isTRUE(all.equal(best, max(comparison$table$log_likelihood))))
})

test_that("the header's interval carries the user's overdispersion", {
  # Z MUST COME FROM A RUN PINNED AT OMEGA = 1 or it is not `Ew^2/Ew`, but the
  # interval is a statement about the DATA and quoting it at 1 would understate it
  # by sqrt(Omega). The two numbers therefore disagree about Omega on purpose, and
  # dropping it from the interval was SILENT until this test existed.
  for (omega in c(1, 4)) {
    comparison <- compare_models(group_data, models = list(level = level_model),
                                 weight = .i$pension, overdispersion = omega)
    d <- comparison$diagnostics
    expected <- format(sqrt(omega * d$Ew2) / d$Ew, digits = 5L)

    shown <- capture.output(print(comparison))
    line <- grep("log A/E", shown, value = TRUE)
    expect_length(line, 1L)
    expect_true(grepl(expected, line, fixed = TRUE))
    expect_true(grepl(sprintf("1 sd at overdispersion %s", omega), line, fixed = TRUE))
  }

  # And at Omega = 4 the interval is twice the Omega = 1 one, while Z is unmoved.
  one <- compare_models(group_data, models = list(level = level_model),
                        weight = .i$pension, overdispersion = 1)
  four <- compare_models(group_data, models = list(level = level_model),
                         weight = .i$pension, overdispersion = 4)
  expect_equal(four$Z, one$Z)
})

test_that("the header names the include, the weight and the test mortality given", {
  comparison <- compare_models(group_data, models = list(level = level_model),
                               include = include(.i$group != "c"),
                               weight = .i$pension,
                               test_mortality = gompertz_mortality(slope_x = 0.2),
                               settings = basis)
  shown <- capture.output(print(comparison))

  # WHAT THE USER WROTE, which identifies a mortality where the object cannot:
  # `it_deparse()` renders an obj leaf as `<mortality_expr>`.
  expect_true(any(grepl("gompertz_mortality\\(slope_x = 0.2\\)", shown)))
  expect_true(any(grepl("include\\(.i\\$group != \"c\"\\)", shown, fixed = FALSE)))
  expect_true(any(grepl("weight          .i\\$pension", shown)))

  # The include narrows the population the header describes, not just the fit.
  expect_lt(comparison$diagnostics$deaths, 5)
})

# ---- a failing candidate ---------------------------------------------------

test_that("a candidate that cannot be fitted does not destroy the comparison", {
  # `fit()` raises, because a fit returns a fit or it fails. A comparison must
  # not: one bad candidate out of fifty cannot be allowed to take the other
  # forty-nine with it.
  expect_error(fit(group_data, collinear_model, settings = basis), "not identifiable")

  # ASSIGNED INSIDE `expect_warning()`, which returns the condition rather than
  # the value of the expression.
  expect_warning(
    comparison <- compare_models(group_data,
                                 models = list(level = level_model,
                                               broken = collinear_model,
                                               group = group_model),
                                 settings = basis),
    "1 of 3 candidates could not be fitted"
  )

  expect_identical(nrow(comparison$table), 3L)
  expect_identical(comparison$table["broken", "status"], "not_identifiable")
  expect_identical(comparison$table["broken", "penalised_log_likelihood"], NA_real_)
  expect_identical(comparison$table["broken", "log_likelihood"], NA_real_)
  expect_identical(comparison$table["broken", "penalty"], NA_real_)

  # LAST, whatever its status, because `order(na.last = TRUE)` puts NA at the end.
  expect_identical(rownames(comparison$table)[[3L]], "broken")
  expect_null(comparison$fits[["broken"]])

  # The other two are unaffected, to the bit.
  expect_identical(comparison$fits[["level"]],
                   fit(group_data, level_model, settings = basis))
  expect_true(is_fit(comparison$fits[["group"]]))
  expect_true(is_fit(best_fit(comparison)))
})

test_that("the failure carries fit()'s diagnosis, not just its status", {
  comparison <- suppressWarnings(
    compare_models(group_data,
                   models = list(level = level_model, broken = collinear_model),
                   settings = basis)
  )

  expect_identical(names(comparison$failures), "broken")
  # The whole message `fit()` would have raised, which names the term and what it
  # duplicates. Reporting only `not_identifiable` would throw that away.
  expect_match(comparison$failures[["broken"]], "Term 3 \\(again\\)")
  expect_match(comparison$failures[["broken"]], "\\[1\\] level")
})

test_that("a comparison in which nothing converged has no best fit", {
  comparison <- suppressWarnings(
    compare_models(group_data, models = list(broken = collinear_model), settings = basis)
  )

  expect_identical(comparison$table$penalised_log_likelihood, NA_real_)
  expect_error(best_fit(comparison), "No candidate converged")

  # And printing it must not fail for want of a maximum to deduct.
  expect_output(print(comparison), "1 candidate, 0 converged")
})

test_that("a failure to converge is a candidate failure too", {
  # A collinear candidate never reaches a log-likelihood; one that runs out of
  # iterations does, and is a different status down a different path.
  comparison <- suppressWarnings(
    compare_models(group_data,
                   models = list(level = level_model, group = group_model),
                   max_iterations = 1, settings = basis)
  )

  expect_identical(comparison$table["level", "status"], "did_not_converge")
  # Singular, which is also the only test of that pluralisation.
  expect_match(comparison$failures[["level"]], "did not converge in 1 iteration;")
})

# ---- the candidate list ----------------------------------------------------

test_that("compare_models() refuses a list that is not models", {
  expect_error(compare_models(group_data, settings = basis), "`models` is required")
  expect_error(compare_models(group_data, models = level_model, settings = basis),
               "must be a list of `model` objects")
  expect_error(compare_models(group_data, models = list(), settings = basis),
               "nothing to compare")
  expect_error(compare_models(group_data, models = list(level_model, reference),
                              settings = basis),
               "`models\\[\\[2\\]\\]` is not a `model`")
})

test_that("candidates with no name are labelled by position", {
  comparison <- compare_models(group_data, models = list(level_model, group_model),
                              settings = basis)
  expect_identical(sort(rownames(comparison$table)), c("model_1", "model_2"))
})

test_that("a named and an unnamed candidate can be mixed", {
  comparison <- compare_models(group_data, models = list(level = level_model, group_model),
                              settings = basis)
  expect_identical(sort(rownames(comparison$table)), c("level", "model_2"))
})

test_that("candidates may not share a name", {
  # The name indexes `fits`, so a duplicate would make one of the two candidates
  # unreachable.
  expect_error(
    compare_models(group_data, models = list(m = level_model, m = group_model),
                   settings = basis),
    "`m` is used more than once"
  )
})

test_that("a comparison of one candidate is allowed", {
  # Degenerate, but a comparison is often built programmatically and refusing it
  # would break a loop that happens to produce one candidate.
  comparison <- compare_models(group_data, models = list(only = level_model),
                               settings = basis)
  expect_identical(nrow(comparison$table), 1L)
  expect_true(is_fit(best_fit(comparison)))
})

# ---- the shared arguments --------------------------------------------------

test_that("include narrows the population every candidate sees", {
  comparison <- compare_models(group_data, models = list(level = level_model),
                               include = include(.i$group == "a"), settings = basis)
  expect_identical(comparison$fits[["level"]],
                   fit(group_data, level_model, include = include(.i$group == "a"),
                       settings = basis))
})

test_that("compare_models() takes one include, not a breakdown", {
  expect_error(
    compare_models(group_data, models = list(level = level_model),
                   include = includes(include(.i$group == "a"), include(.i$group == "b")),
                   settings = basis),
    "takes a single `include`"
  )
})

test_that("overdispersion is required, as everywhere else", {
  expect_error(compare_models(group_data, models = list(level = level_model)),
               "`overdispersion` is required")
})

test_that("there is no `start` argument, and every candidate begins at zero", {
  # DELIBERATELY ABSENT (Tim, 2026-09-10). Newton-Raphson with relaxation gets
  # there from the reference mortality, and a per-candidate `start` could only be
  # a list, since candidates differ in term count. A user who needs the iteration
  # to begin somewhere else translates `log mu^ref` instead.
  expect_false("start" %in% names(formals(compare_models)))
  expect_true("start" %in% names(formals(fit)))

  # And the engine is reached with zeros, which is the same answer `fit()` gives
  # with its own default -- to the bit.
  comparison <- two()
  expect_identical(comparison$fits[["level"]],
                   fit(group_data, level_model, start = 0, settings = basis))
})

test_that("a bad numerical setting is the caller's error, not a candidate failure", {
  # Nothing to record against a candidate here: the method itself is wrong, so
  # every candidate would fail the same way and raising is right.
  expect_error(
    compare_models(group_data, models = list(level = level_model),
                   armijo = 0.5, settings = basis),
    "strictly between 0 and 0.5"
  )
})

# ---- printing --------------------------------------------------------------

test_that("print() shows the ranking, the scale and any diagnosis", {
  comparison <- two()
  shown <- capture.output(print(comparison))

  expect_match(shown[[1L]], "<model_comparison: 2 candidates, 2 converged>")
  # The printed columns are Tim's four, not the table's field names.
  expect_true(any(grepl("^ *L_P +L +p +k$", shown)))
  expect_true(any(grepl("^  max L_P", shown)))
  expect_true(any(grepl("^  overdispersion  1$", shown)))

  # THE STATUS COLUMN IS NOT SHOWN WHEN THERE IS NOTHING TO SAY. It would read
  # `converged` all the way down and wrap the table. It is still in `$table`.
  expect_false(any(grepl("converged", shown[-1L])))
  expect_true("status" %in% names(comparison$table))

  # Returns its argument invisibly, so `x` still prints once at the console.
  expect_identical(capture.output(returned <- print(comparison)), shown)
  expect_identical(returned, comparison)
})

test_that("print() counts the failures and explains them", {
  comparison <- suppressWarnings(
    compare_models(group_data,
                   models = list(level = level_model, broken = collinear_model),
                   settings = basis)
  )
  shown <- capture.output(print(comparison))

  expect_match(shown[[1L]], "<model_comparison: 2 candidates, 1 converged>")
  expect_true(any(grepl("^broken: The model is not identifiable", shown)))

  # And here the status column IS shown, because it now distinguishes the rows.
  expect_true(any(grepl("not_identifiable", shown)))
})
