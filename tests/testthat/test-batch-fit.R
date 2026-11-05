# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests for `fit()` and `compare_models()` inside `batch()`.
#
# THE ORACLE IS THE LONE CALL, as it is for `aev()` in `test-batch.R`. A batched
# fit must be `identical()` to the same fit made on its own, because everything
# `batch()` does is scheduling. The weight is an amount wherever it can be, so
# that Z is not 1 and a Z obtained the wrong way would show.

reference <- mortality_const(log_mu = -4)

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

# Term 3 duplicates term 1, so the information matrix is singular.
collinear_model <- model(
  ref_mortality = reference,
  covariates = covariates(level = 1, group_a = .i$group == "a", again = 1)
)

# ---- identity with the lone call ---------------------------------------------

test_that("a batched fit is identical to the same fit alone", {
  b <- batch(
    .exp_data = group_data, .overdispersion = 2, .weight = .i$pension,
    level = fit(model = level_model, include = include(.i$group != "c")),
    group = fit(model = group_model)
  )

  expect_identical(
    b$level,
    fit(group_data, level_model, weight = .i$pension, overdispersion = 2,
        include = include(.i$group != "c"))
  )
  expect_identical(
    b$group,
    fit(group_data, group_model, weight = .i$pension, overdispersion = 2)
  )

  # Z really is away from 1, or the identity says nothing about how it was got.
  expect_false(isTRUE(all.equal(b$level$Z, 1)))
})

test_that("a batched comparison is identical to the same comparison alone", {
  candidates <- list(level = level_model, group = group_model)

  b <- batch(
    .exp_data = group_data, .overdispersion = 1, .weight = .i$pension,
    ranked = compare_models(models = candidates)
  )

  expect_identical(
    b$ranked,
    compare_models(group_data, candidates, weight = .i$pension, overdispersion = 1)
  )
})

test_that("A/Es, fits and comparisons share one crossing for their population passes", {
  real_run <- cpp_veil_run
  real_fits <- cpp_veil_fit_runs
  crossings <- 0L
  fit_calls <- 0L
  models_fitted <- 0L
  local_mocked_bindings(
    cpp_veil_run = function(...) {
      crossings <<- crossings + 1L
      real_run(...)
    },
    cpp_veil_fit_runs = function(fits, ...) {
      fit_calls <<- fit_calls + 1L
      models_fitted <<- models_fitted + length(fits)
      real_fits(fits, ...)
    }
  )

  b <- batch(
    .exp_data = group_data, .overdispersion = 1, .weight = .i$pension,
    crude  = aev(mortality = reference),
    level  = fit(model = level_model),
    ranked = compare_models(models = list(level = level_model, group = group_model))
  )

  # One pass over the data for everything that is not an iteration, and ONE
  # engine call iterating every model together -- the lone fit and both of the
  # comparison's candidates.
  expect_identical(crossings, 1L)
  expect_identical(fit_calls, 1L)
  expect_identical(models_fitted, 3L)

  expect_true(is_aev(b$crude))
  expect_true(is_fit(b$level))
  expect_true(is_model_comparison(b$ranked))
})

test_that("fits at different time scales each match their own lone call", {
  # Varying with age, so that the integration interval actually matters.
  ageing_mortality <- mortality(-10 + 0.1 * .x)
  ageing <- model(ref_mortality = ageing_mortality, covariates = covariates(level = 1))

  b <- batch(
    .exp_data = group_data, .overdispersion = 1,
    yearly  = fit(model = ageing, time_scale = 1),
    monthly = fit(model = ageing, time_scale = 1 / 12)
  )

  expect_identical(b$yearly, fit(group_data, ageing, overdispersion = 1, time_scale = 1))
  expect_identical(b$monthly, fit(group_data, ageing, overdispersion = 1, time_scale = 1 / 12))
  expect_false(isTRUE(all.equal(b$yearly$beta, b$monthly$beta)))
})

test_that("a threads argument inside a batched fit is ignored, not refused", {
  b <- batch(
    .exp_data = group_data, .overdispersion = 1, .threads = 1,
    level = fit(model = level_model, threads = 2)
  )
  expect_identical(b$level, fit(group_data, level_model, overdispersion = 1))
})

# ---- the new defaults --------------------------------------------------------

test_that(".Z applies where a fit names neither half of its pair", {
  b <- batch(
    .exp_data = group_data, .overdispersion = 1, .weight = .i$pension, .Z = 2.5,
    given  = fit(model = level_model),
    own    = fit(model = level_model, test_mortality = reference),
    own_Z  = fit(model = level_model, Z = 4)
  )

  expect_identical(b$given, fit(group_data, level_model, weight = .i$pension,
                                overdispersion = 1, Z = 2.5))

  # Naming the test mortality is naming the scale, so the batch's `.Z` must not
  # also arrive -- which would be the pair the lone call refuses.
  expect_identical(b$own, fit(group_data, level_model, weight = .i$pension,
                              overdispersion = 1, test_mortality = reference))
  expect_identical(b$own_Z$Z, 4)
})

test_that(".test_mortality applies to fits and comparisons alike", {
  b <- batch(
    .exp_data = group_data, .overdispersion = 1, .weight = .i$pension,
    .test_mortality = reference,
    level  = fit(model = level_model),
    ranked = compare_models(models = list(level = level_model))
  )

  expect_identical(b$level, fit(group_data, level_model, weight = .i$pension,
                                overdispersion = 1, test_mortality = reference))
  expect_identical(b$ranked$Z, b$level$Z)
})

test_that("the method's defaults come from the batch, and an analysis's own value wins", {
  b <- suppressWarnings(batch(
    .exp_data = group_data, .overdispersion = 1, .max_iterations = 1,
    starved = fit(model = group_model),
    own     = fit(model = group_model, max_iterations = 25)
  ))

  expect_s3_class(b$starved, "fit_failure")
  expect_identical(b$starved$status, "did_not_converge")
  expect_identical(b$own, fit(group_data, group_model, overdispersion = 1))
})

test_that("a batch may not give both halves of a pair", {
  expect_error(
    batch(.exp_data = group_data, .overdispersion = 1, .Z = 2, .test_mortality = reference,
          level = fit(model = level_model)),
    "Give either `.test_mortality` or `.Z`, not both"
  )
  expect_error(
    batch(.exp_data = group_data, .overdispersion = 1, .val_relevance = 1, .val_distance = 0,
          level = fit(model = level_model)),
    "Give either `.val_relevance` or `.val_distance`, not both"
  )
})

test_that("an analysis giving both halves of a pair is refused by name", {
  expect_error(
    batch(.exp_data = group_data, .overdispersion = 1,
          both = fit(model = level_model, Z = 2, test_mortality = reference)),
    "`both` gives both `test_mortality` and `Z`"
  )
})

# ---- failure -----------------------------------------------------------------

test_that("a lone fit that fails raises a `fit_failure` carrying its status", {
  failure <- tryCatch(fit(group_data, collinear_model, overdispersion = 1),
                      fit_failure = function(failure) failure)

  expect_s3_class(failure, "fit_failure")
  expect_s3_class(failure, "error")
  expect_identical(failure$status, "not_identifiable")
  expect_match(conditionMessage(failure), "not identifiable")
})

test_that("a failed fit is held on its element and the rest of the batch survives", {
  expect_warning(
    b <- batch(
      .exp_data = group_data, .overdispersion = 1,
      level  = fit(model = level_model),
      broken = fit(model = collinear_model),
      crude  = aev(mortality = reference)
    ),
    "1 of 3 analyses in the batch could not be fitted: `broken` \\(not_identifiable\\)"
  )

  expect_identical(b$level, fit(group_data, level_model, overdispersion = 1))
  expect_true(is_aev(b$crude))

  # THE HELD ERROR IS THE ONE THE LONE CALL RAISES, to the object, so `stop()`
  # on it is that call's failure exactly.
  alone <- tryCatch(fit(group_data, collinear_model, overdispersion = 1),
                    fit_failure = function(failure) failure)
  expect_identical(b$broken, alone)
  expect_error(stop(b$broken), conditionMessage(alone), fixed = TRUE)
})

test_that("a comparison with a failing candidate still comes back as a comparison", {
  expect_warning(
    b <- batch(
      .exp_data = group_data, .overdispersion = 1,
      ranked = compare_models(models = list(level = level_model, broken = collinear_model))
    ),
    "could not be fitted"
  )

  # The candidate's failure belongs to the comparison, which records it; the
  # element itself succeeded.
  expect_true(is_model_comparison(b$ranked))
  expect_identical(b$ranked$table["broken", "status"], "not_identifiable")
})

test_that("an error that is not a failed fit still stops the batch", {
  expect_error(
    batch(.exp_data = group_data, .overdispersion = 1,
          level = fit(model = level_model),
          wrong = fit(model = level_model, Z = -1)),
    "`Z` must be a single positive finite number"
  )
})

# ---- what a batched fit refuses ----------------------------------------------

test_that("a fit or a comparison with no model is refused by name", {
  expect_error(
    batch(.exp_data = group_data, .overdispersion = 1, nothing = fit()),
    "`nothing` has no `model`"
  )
  expect_error(
    batch(.exp_data = group_data, .overdispersion = 1, nothing = compare_models()),
    "`nothing` has no `models`"
  )
})

test_that("a batched fit with no data anywhere is refused by name", {
  expect_error(
    batch(.overdispersion = 1, nowhere = fit(model = level_model)),
    "`nowhere` has no `exp_data`"
  )
})

test_that("a batched fit may not use another element's result", {
  expect_error(
    batch(.exp_data = group_data, .overdispersion = 1,
          first  = fit(model = level_model),
          second = fit(model = model(ref_mortality = first$mortality,
                                     covariates = covariates(level = 1)))),
    "another analysis in the same"
  )
})
