# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `measure_overdispersion()` is a Pearson statistic over the cells of a
# breakdown. THE ORACLES:
#
#   * `aev()` itself, at an overdispersion of one, from which the statistic is a
#     line of arithmetic;
#   * data simulated from the truth, where deaths are independent and the true
#     overdispersion is one;
#   * the same data with every life duplicated k times, where the true
#     overdispersion is k -- and where, cell by cell, (A - E)^2 scales by k^2
#     and V by k, so the measure scales by exactly k for ANY weight.

od_population <- function(n, seed = 1) {
  set.seed(seed)
  data.frame(
    birth = datey::datey(2015 - sample(60:85, n, replace = TRUE)),
    pension = round(exp(rnorm(n, log(12000), 0.7)))
  )
}

od_truth <- gompertz_mortality(log_mu_0 = -3.6, slope_t = 0)

od_data <- function(n, seed = 1) {
  simulate(od_truth, seed = seed, lives = od_population(n, seed), exp_start = 2015, exp_end = 2020)
}

# Every life repeated `times` times: perfectly dependent deaths.
od_duplicate <- function(data, times) {
  columns <- unclass(data)
  attributes(columns) <- list(names = names(columns))
  exp_data(
    lapply(columns, function(column) rep(column, times)),
    exp_start = exp_start(data),
    exp_end = exp_end(data))
}

od_small <- od_data(5000)
od_large <- od_data(20000)
od_cells <- ages(60, 90, by = 5)

test_that("the measure is the Pearson statistic of aev() over its degrees of freedom", {
  check <- function(measure, cells, k) {
    values <- unclass(cells)
    expected <- sum((values$A - values$E)^2 / values$V) / (length(cells) - k)
    expect_equal(measure$overdispersion, expected)
    expect_identical(measure$cells, length(cells))
  }

  check(
    measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 0),
    aev(od_small, od_truth, breakdown = od_cells, overdispersion = 1),
    0L)

  # Every pronoun argument passes through as it does in `aev()`: a weight, a
  # population and a relevance, the last because V already carries it.
  check(
    measure_overdispersion(
      od_small,
      od_truth,
      include = include(.i$pension > 5000),
      weight = .i$pension,
      val_distance = (2020 - .t) / 10,
      breakdown = od_cells,
      parameter_count = 2),
    aev(
      od_small,
      od_truth,
      include = include(.i$pension > 5000),
      weight = .i$pension,
      val_distance = (2020 - .t) / 10,
      breakdown = od_cells,
      overdispersion = 1),
    2L)
})

test_that("the overdispersion a `settings` carries is not used", {
  plain <- measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 0)
  given <- measure_overdispersion(
    od_small,
    od_truth,
    breakdown = od_cells,
    parameter_count = 0,
    settings = settings(overdispersion = 3))
  expect_identical(given, plain)
})

test_that("on data simulated from the truth the measure is near one", {
  measure <- measure_overdispersion(
    od_large,
    od_truth,
    breakdown = ages(60, 90, by = 1),
    parameter_count = 0)
  expect_identical(measure$cells, 30L)
  expect_lt(measure$lower_95pc, 1)
  expect_gt(measure$upper_95pc, 1)
  expect_lt(abs(measure$overdispersion - 1), 0.2)
})

test_that("a fitted model is measured with its parameters taken off", {
  reference <- gompertz_mortality(log_mu_0 = -3.3, slope_t = 0)
  f <- fit(od_large, model(reference, covariates(level = 1)), overdispersion = 1)
  measure <- measure_overdispersion(
    od_large,
    f$mortality,
    breakdown = ages(60, 90, by = 1),
    parameter_count = length(coef(f)))
  expect_identical(measure$degrees_of_freedom, 29L)
  expect_lt(measure$lower_95pc, 1)
  expect_gt(measure$upper_95pc, 1)
})

test_that("duplicating every life k times multiplies the measure by exactly k, for any weight", {
  tripled <- od_duplicate(od_small, 3L)

  lives <- measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 0)
  lives_3 <- measure_overdispersion(tripled, od_truth, breakdown = od_cells, parameter_count = 0)
  expect_equal(lives_3$overdispersion, 3 * lives$overdispersion)

  weighted <- suppressWarnings(measure_overdispersion(
    od_small,
    od_truth,
    weight = .i$pension,
    breakdown = od_cells,
    parameter_count = 0))
  weighted_3 <- suppressWarnings(measure_overdispersion(
    tripled,
    od_truth,
    weight = .i$pension,
    breakdown = od_cells,
    parameter_count = 0))
  expect_equal(weighted_3$overdispersion, 3 * weighted$overdispersion)

  # And so on the large dataset, where the truth is one, three is recovered.
  large_3 <- measure_overdispersion(
    od_duplicate(od_large, 3L),
    od_truth,
    breakdown = ages(60, 90, by = 1),
    parameter_count = 0)
  expect_lt(large_3$lower_95pc, 3)
  expect_gt(large_3$upper_95pc, 3)
})

test_that("the parameter count comes off the cell count", {
  none <- measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 0)
  two <- measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 2)
  expect_identical(none$cells, 6L)
  expect_identical(two$degrees_of_freedom, 4L)
  expect_identical(two$parameter_count, 2L)
  expect_equal(two$overdispersion * 4, none$overdispersion * 6)
})

test_that("a cell with nothing in it is not counted", {
  covered <- measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 1)
  beyond <- measure_overdispersion(
    od_small,
    od_truth,
    breakdown = ages(60, 120, by = 5),
    parameter_count = 1)
  expect_identical(beyond$cells, 6L)
  expect_identical(beyond, covered)
})

test_that("a mortality that underflows against real deaths gives an infinite measure", {
  measure <- measure_overdispersion(
    od_small,
    mortality_const(log_mu = -800),
    breakdown = od_cells,
    parameter_count = 0)
  expect_identical(measure$overdispersion, Inf)
})

test_that("cells with few effective expected deaths draw a warning naming them", {
  few <- od_data(300)
  expect_warning(
    measure_overdispersion(few, od_truth, breakdown = od_cells, parameter_count = 0),
    "fewer than 5 effective expected deaths \\(.*60-65.*\\)")
  expect_no_warning(
    measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 0))
})

test_that("the warning counts effective deaths, which a weight reduces", {
  # Lives: plenty in every cell. The same cells under a skewed weight hold far
  # fewer effective deaths, and the youngest cell drops below five.
  sparse <- od_data(1200)
  expect_no_warning(
    measure_overdispersion(sparse, od_truth, breakdown = od_cells, parameter_count = 0))
  expect_warning(
    measure_overdispersion(
      sparse,
      od_truth,
      weight = .i$pension^3,
      breakdown = od_cells,
      parameter_count = 0),
    "effective expected deaths")
})

test_that("no degrees of freedom left is an error", {
  expect_error(
    measure_overdispersion(od_small, od_truth, parameter_count = 1),
    "no degrees of freedom left.*1 cell with exposure, less 1 fitted parameter")
  expect_error(
    measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 6),
    "6 cells with exposure, less 6 fitted parameters")

  # One cell and nothing fitted leaves one degree of freedom, which is allowed.
  expect_identical(
    measure_overdispersion(od_small, od_truth, parameter_count = 0)$degrees_of_freedom,
    1L)
})

test_that("`parameter_count` is required and must be a whole number", {
  expect_error(
    measure_overdispersion(od_small, od_truth, breakdown = od_cells),
    "`parameter_count` is required")
  for (bad in list(-1, 1.5, NA_real_, Inf, c(1, 2), "1", TRUE)) {
    expect_error(
      measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = bad),
      "must be a single whole number")
  }
  expect_error(measure_overdispersion(od_small, parameter_count = 0), "`mortality` is required")
})

test_that("the interval and the standard error are the chi-squared ones", {
  measure <- measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 1)
  degrees <- 5
  expect_equal(measure$std_error, measure$overdispersion * sqrt(2 / degrees))
  expect_equal(measure$lower_95pc, measure$overdispersion * degrees / qchisq(0.975, degrees))
  expect_equal(measure$upper_95pc, measure$overdispersion * degrees / qchisq(0.025, degrees))
  expect_lt(measure$lower_95pc, measure$overdispersion)
  expect_gt(measure$upper_95pc, measure$overdispersion)
})

test_that("print() names every figure it shows", {
  measure <- measure_overdispersion(od_small, od_truth, breakdown = od_cells, parameter_count = 1)
  output <- capture.output(result <- print(measure))
  expect_identical(result, measure)
  expect_match(output[1], "<overdispersion measure>", fixed = TRUE)
  shown <- formatC(measure$overdispersion, format = "f", digits = 3L)
  expect_true(any(grepl(shown, output, fixed = TRUE)))
  expect_true(any(grepl("(1 sd)", output, fixed = TRUE)))
  expect_true(any(grepl("95% interval", output, fixed = TRUE)))
  expect_true(any(grepl("5 (6 cells less 1 parameter)", output, fixed = TRUE)))
})
