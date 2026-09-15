# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `diagnose()` is an `aev()` per candidate per breakdown in one engine pass, plus
# each candidate's summary. THE ORACLE IS `aev()` ITSELF: every A/E a diagnosis
# holds must be identical to the bit to the one `aev()` gives when asked directly
# with the same arguments, because both build the same plan.

reference <- mortality_const(log_mu = -4)
basis <- settings(overdispersion = 1)

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

level_model <- model(reference, covariates(level = 1))
group_model <- model(reference, covariates(a = .i$group == "a", b = .i$group == "b",
                                           c = .i$group == "c"))

by_group <- categories(.i$group, "a", "b", "c")
by_amount <- bands(.i$pension, thresholds = c(10000, 20000))

test_that("a lone fit is diagnosed as one candidate over one ungrouped A/E", {
  f <- fit(group_data, level_model, settings = basis)
  d <- diagnose(group_data, f)

  expect_true(is_diagnosis(d))
  expect_identical(names(d$candidates), "fit")
  expect_identical(d$breakdowns, "all")
  expect_identical(d$candidates$fit$summary, summary(f))
  expect_null(d$pv_impact)

  # THE ORACLE: the same A/E asked for directly.
  direct <- aev(group_data, mortality = f$mortality, settings = basis)
  expect_identical(d$candidates$fit$aev$all, direct)
})

test_that("every converged candidate is diagnosed, best first, on every breakdown", {
  comparison <- compare_models(group_data, list(level = level_model, group = group_model),
                               settings = basis)
  d <- diagnose(group_data, comparison, breakdown = list(group = by_group, amount = by_amount))

  expect_identical(names(d$candidates), rownames(comparison$table))
  expect_identical(d$breakdowns, c("group", "amount"))

  for (label in names(d$candidates)) {
    candidate <- d$candidates[[label]]
    expect_identical(candidate$summary, summary(comparison$fits[[label]]))
    expect_identical(names(candidate$aev), c("group", "amount"))
    expect_identical(
      candidate$aev$group,
      aev(group_data, mortality = comparison$fits[[label]]$mortality, breakdown = by_group,
          settings = basis)
    )
    expect_identical(
      candidate$aev$amount,
      aev(group_data, mortality = comparison$fits[[label]]$mortality, breakdown = by_amount,
          settings = basis)
    )
  }

  # The labels travel: an aev from a breakdown carries the breakdown's names.
  expect_identical(names(d$candidates$group$aev$group), c("a", "b", "c"))
})

test_that("the whole diagnosis is one pass over the data", {
  comparison <- compare_models(group_data, list(level = level_model, group = group_model),
                               settings = basis)
  runs <- 0L
  original <- cpp_veil_run
  local_mocked_bindings(
    cpp_veil_run = function(...) {
      runs <<- runs + 1L
      original(...)
    }
  )
  d <- diagnose(group_data, comparison, breakdown = list(by_group, by_amount))
  expect_identical(runs, 1L)
  # Two candidates on two breakdowns is four A/Es, and they all came from that pass.
  expect_identical(sum(lengths(lapply(d$candidates, function(c) c$aev))), 4L)
})

test_that("the fit's own include, weight and overdispersion are the defaults", {
  f <- fit(group_data, level_model, include = include(.i$group != "c"),
           weight = .i$pension, overdispersion = 2.5)
  d <- diagnose(group_data, f, breakdown = by_group)

  expect_identical(d$overdispersion, 2.5)
  expect_identical(d$labels$include, f$labels$include)
  expect_identical(d$labels$weight, ".i$pension")
  expect_identical(
    d$candidates$fit$aev[[1L]],
    aev(group_data, mortality = f$mortality, include = include(.i$group != "c"),
        weight = .i$pension, breakdown = by_group, overdispersion = 2.5)
  )

  # GIVEN AS NULL MEANS NONE, which `missing()` tells apart from not given.
  lives <- diagnose(group_data, f, breakdown = by_group, weight = NULL, include = NULL,
                    overdispersion = 1)
  expect_identical(lives$labels$weight, "lives")
  expect_identical(lives$labels$include, "all records")
  expect_identical(lives$overdispersion, 1)
  expect_identical(
    lives$candidates$fit$aev[[1L]],
    aev(group_data, mortality = f$mortality, breakdown = by_group, overdispersion = 1)
  )
  expect_false(identical(lives$candidates$fit$aev[[1L]], d$candidates$fit$aev[[1L]]))

  # And an override is what was written.
  other <- diagnose(group_data, f, weight = .i$pension / 1000)
  expect_identical(other$labels$weight, "(.i$pension / 1000)")
})

test_that("breakdowns are labelled by name, else by what was written", {
  f <- fit(group_data, level_model, settings = basis)

  expect_identical(diagnose(group_data, f, breakdown = by_group)$breakdowns, "by_group")
  expect_identical(diagnose(group_data, f, breakdown = list(by_group, by_amount))$breakdowns,
                   c("by_group", "by_amount"))
  expect_identical(diagnose(group_data, f, breakdown = list(g = by_group, by_amount))$breakdowns,
                   c("g", "by_amount"))
  expect_identical(
    diagnose(group_data, f, breakdown = categories(.i$group, "a", "b"))$breakdowns,
    'categories(.i$group, "a", "b")'
  )
  expect_error(diagnose(group_data, f, breakdown = 3), "must be an `include`")
})

test_that("a failed candidate is left out, and nothing at all is an error", {
  broken <- model(reference, covariates(a = .i$group == "a", also_a = .i$group == "a"))
  comparison <- suppressWarnings(
    compare_models(group_data, list(level = level_model, broken = broken), settings = basis)
  )
  d <- diagnose(group_data, comparison)
  expect_identical(names(d$candidates), "level")

  nothing <- suppressWarnings(compare_models(group_data, list(broken = broken), settings = basis))
  expect_error(diagnose(group_data, nothing), "No candidate")
  expect_error(diagnose(group_data, level_model), "must be a `fit` or a `model_comparison`")
  expect_error(diagnose(group_data), "`x` is required")
})

test_that("a diagnosis prints each candidate's tables and each A/E", {
  comparison <- compare_models(group_data, list(level = level_model, group = group_model),
                               settings = basis)
  d <- diagnose(group_data, comparison, breakdown = list(group = by_group))
  printed <- capture.output(print(d))

  expect_match(printed[[1L]], "<diagnosis: 2 candidates, 1 breakdown>")
  expect_identical(grep("^---- ", printed, value = TRUE),
                   paste("----", rownames(comparison$table)))
  expect_identical(sum(printed == "A/E by group"), 2L)
  # The three-term candidate has a correlation to show; the one-term one has not.
  expect_identical(sum(printed == "correlation"), 1L)
  expect_true(any(grepl("^  overdispersion  1$", printed)))
})

test_that("plot() draws one chart per candidate per breakdown, in a grid", {
  comparison <- compare_models(group_data, list(level = level_model, group = group_model),
                               settings = basis)
  d <- diagnose(group_data, comparison, breakdown = list(group = by_group, amount = by_amount))

  charts <- diagnosis_charts(d)
  expect_length(charts, 4L)
  expect_true(all(vapply(charts, ggplot2::is_ggplot, logical(1L))))
  # Row-major: the first row is the best candidate across both breakdowns.
  titles <- vapply(charts, function(chart) chart$labels$title, character(1L))
  subtitles <- vapply(charts, function(chart) chart$labels$subtitle, character(1L))
  expect_identical(titles, rep(rownames(comparison$table), each = 2L))
  expect_identical(subtitles, rep(c("group", "amount"), times = 2L))

  # Anything else goes to every chart; the two names `plot()` owns are refused.
  wide <- diagnosis_charts(d, log_range = c(-0.6, 0.6), log_step = 0.2)
  expect_length(wide, 4L)
  expect_error(diagnosis_charts(d, title = "x"), "set by `plot()`", fixed = TRUE)

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_identical(plot(d), d)
})
