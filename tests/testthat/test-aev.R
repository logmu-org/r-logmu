# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `aev()` is where the R front end and the engine finally meet: a breakdown of G
# elements becomes G specifications run in one crossing, and the two levels of
# naming travel from the `includes` to the result.
#
# THE STRONGEST ORACLE HERE IS THE DISJOINT SUM. On a breakdown deliberately
# built to be disjoint and to cover everybody, the elements must add up to the
# ungrouped answer. That is an oracle for OUR fixtures, not a rule logmu
# enforces -- a breakdown is not a partition and elements may overlap.

data <- exp_data(
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

# The same five individuals with a `log_mu` column bolted on, for the underflow
# tests below. Rebuilt rather than modified: an `exp_data` validates its columns
# on construction, and reaching into one with `$<-` would go round that.
pensions <- c(5000, 12000, 30000, 8000, 15000)

with_log_mu <- function(log_mu) {
  exp_data(
    list(
      birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
      pension   = pensions,
      log_mu    = log_mu,
      E2R_start = datey::datey(rep(2015, 5)),
      E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
      E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )
}

flat <- mortality_const(log_mu = -4)
basis <- settings(overdispersion = 2)

test_that("an ungrouped aev is a single record", {
  res <- aev(data, mortality = flat, settings = basis)

  expect_true(is_aev(res))
  expect_equal(length(res), 1L)
  expect_null(names(res))
  expect_null(group_names(res))

  # Two deaths in the data, unweighted, so A counts them.
  expect_equal(res$A, 2)
  expect_gt(res$E, 0)
})

test_that("a breakdown gives one record per element, carrying both label levels", {
  res <- aev(data,
             mortality = flat,
             breakdown = bands(.i$pension, thresholds = c(10000, 20000)),
             settings  = basis)

  expect_equal(length(res), 3L)
  expect_identical(names(res), c("< 10000", "10000-20000", ">= 20000"))
  expect_identical(group_names(res), rep("pension", 3L))

  # ONE aev, not a list of them -- the type is vectorised and that was the
  # reason for having a single `aev()` rather than an `aev_multi()`.
  expect_true(is_aev(res))
})

test_that("disjoint covering groups sum to the ungrouped answer", {
  # The oracle. Thresholds chosen so every record falls in exactly one band.
  grouped <- aev(data, mortality = flat, weight = .i$pension,
                 breakdown = bands(.i$pension, thresholds = c(10000, 20000)),
                 settings = basis)
  whole <- aev(data, mortality = flat, weight = .i$pension, settings = basis)

  expect_equal(sum(grouped$A), whole$A, tolerance = 1e-12)
  expect_equal(sum(grouped$E), whole$E, tolerance = 1e-12)
  expect_equal(sum(grouped$V), whole$V, tolerance = 1e-12)
})

test_that("an age breakdown sums the same way", {
  grouped <- aev(data, mortality = flat,
                 breakdown = ages(60, 85, by = 5), settings = basis)
  whole <- aev(data, mortality = flat, settings = basis)

  # Every individual is between 60 and 85 throughout, so the bands cover the
  # whole exposure. If that stops being true this fails, which is the point.
  expect_equal(sum(grouped$E), whole$E, tolerance = 1e-12)
  expect_identical(group_names(grouped), rep("age", 5L))
})

test_that("include narrows the population and breakdown divides what is left", {
  males <- include(.i$male)

  # Complementary populations add back to everybody.
  m <- aev(data, mortality = flat, include = males, settings = basis)
  f <- aev(data, mortality = flat, include = include(!.i$male), settings = basis)
  whole <- aev(data, mortality = flat, settings = basis)
  expect_equal(m$E + f$E, whole$E, tolerance = 1e-12)

  # EVERY BREAKDOWN ELEMENT IS INTERSECTED WITH THE POPULATION. Without that the
  # groups below would sum to the whole portfolio rather than to the males.
  split_males <- aev(data, mortality = flat, include = males,
                     breakdown = bands(.i$pension, thresholds = c(10000, 20000)),
                     settings = basis)
  expect_equal(sum(split_males$E), m$E, tolerance = 1e-12)
  expect_lt(sum(split_males$E), whole$E)
})

test_that("an indicator works as the population include", {
  # `include(...)` RETURNS AN INDICATOR by design, and an indicator carries an
  # `ast` rather than `terms`. It is the ordinary shape a population include
  # arrives in, and it used to fail with "expected 'list' actual 'NULL'".
  males <- include(.i$male)
  expect_true(is_indicator(males))

  by_indicator <- aev(data, mortality = flat, include = males, settings = basis)
  # Same question written as a band-style gate, which always worked.
  by_terms <- aev(data, mortality = flat,
                  include = males & age(0, 200), settings = basis)

  expect_equal(by_indicator$E, by_terms$E, tolerance = 1e-12)
  expect_gt(by_indicator$E, 0)
})

test_that("a weight may be an expression, a variable object, or absent", {
  by_expression <- aev(data, mortality = flat, weight = .i$pension, settings = basis)

  amounts <- static_variable(.i$pension)
  by_object <- aev(data, mortality = flat, weight = amounts, settings = basis)

  expect_equal(by_object$E, by_expression$E, tolerance = 1e-14)
  expect_equal(by_object$V, by_expression$V, tolerance = 1e-14)

  # No weight is a count of lives, so E and V agree up to the overdispersion.
  lives <- aev(data, mortality = flat, overdispersion = 1)
  expect_equal(lives$V, lives$E, tolerance = 1e-12)
})

test_that("an indicator weight counts the same exposure as the matching subset", {
  # An indicator is BOTH a variable and an include, and for a time-invariant
  # {0,1} function the two readings coincide: weighting by it and restricting to
  # it select the same exposure. That is why one type can be both.
  #
  # NOT A WITNESS TO THE OBJ-LEAF LOWERING ORDER. Checked by experiment: swapping
  # `variable` and `include` in `ingestObj` leaves the whole suite passing,
  # because a single-gate include lowers back to that gate's value. The order
  # there is the correct reading rather than an observable one.
  males <- indicator(.i$male)

  as_weight <- aev(data, mortality = flat, weight = males, overdispersion = 1)
  as_population <- aev(data, mortality = flat, include = males, overdispersion = 1)

  expect_equal(as_weight$E, as_population$E, tolerance = 1e-12)
  expect_lt(as_weight$E, aev(data, mortality = flat, overdispersion = 1)$E)
})

test_that("a mortality may be written as an expression over mortality objects", {
  plain <- aev(data, mortality = flat, overdispersion = 1)
  raised <- aev(data, mortality = flat + log(2), overdispersion = 1)

  # Adding in log mu is multiplying mu, so E doubles and A does not move.
  expect_equal(raised$E, plain$E * 2, tolerance = 1e-12)
  expect_equal(raised$A, plain$A)
})

test_that("overdispersion scales V alone and the call beats the settings", {
  plain <- aev(data, mortality = flat, weight = .i$pension, overdispersion = 1)
  dispersed <- aev(data, mortality = flat, weight = .i$pension, overdispersion = 4)

  expect_identical(dispersed$A, plain$A)
  expect_identical(dispersed$E, plain$E)
  expect_equal(dispersed$V, plain$V * 4, tolerance = 1e-14)

  # A direct argument overrides the settings object rather than being ignored.
  overridden <- aev(data, mortality = flat, weight = .i$pension,
                    settings = basis, overdispersion = 4)
  expect_equal(overridden$V, dispersed$V, tolerance = 1e-14)
})

test_that("the time scale can be set directly or through settings", {
  monthly_settings <- settings(overdispersion = 1, time_scale = 1 / 12)

  by_settings <- aev(data, mortality = flat, settings = monthly_settings)
  by_argument <- aev(data, mortality = flat, overdispersion = 1, time_scale = 1 / 12)
  quarterly <- aev(data, mortality = flat, overdispersion = 1)

  expect_equal(by_settings$E, by_argument$E, tolerance = 1e-14)

  # A constant mortality integrates exactly at any scale, so the totals agree --
  # which means this pins that both routes were ACCEPTED, not that either bit.
  expect_equal(by_settings$E, quarterly$E, tolerance = 1e-9)

  expect_error(aev(data, mortality = flat, overdispersion = 1, time_scale = 0.5),
               "must be one of")
})

test_that("the thread count cannot move an answer", {
  single <- aev(data, mortality = flat, weight = .i$pension,
                breakdown = bands(.i$pension, thresholds = c(10000, 20000)),
                settings = basis, threads = 1L)
  many <- aev(data, mortality = flat, weight = .i$pension,
              breakdown = bands(.i$pension, thresholds = c(10000, 20000)),
              settings = basis, threads = 4L)

  expect_identical(single$A, many$A)
  expect_identical(single$E, many$E)
  expect_identical(single$V, many$V)
})

test_that("a single include is accepted as a breakdown of one", {
  one_band <- aev(data, mortality = flat, breakdown = age(60, 85), settings = basis)

  expect_equal(length(one_band), 1L)
  expect_identical(group_names(one_band), "age")
  expect_identical(names(one_band), "60-85")
})

test_that("required arguments are required and the wrong types are refused", {
  expect_error(aev(data, mortality = flat), "`overdispersion` is required")
  expect_error(aev(data, settings = basis), "`mortality` is required")

  # `include` is scalar; several go in `breakdown`. The signatures differ so a
  # mix-up errors rather than quietly analysing the wrong thing.
  expect_error(aev(data, mortality = flat, include = ages(60, 80, by = 10), settings = basis),
               "single `include`")
  expect_error(aev(data, mortality = flat, include = "male", settings = basis),
               "must be an `include`")
  expect_error(aev(data, mortality = flat, breakdown = "amounts", settings = basis),
               "must be an `include` or an `includes`")

  expect_error(aev(list(a = 1), mortality = flat, settings = basis), "must be `exp_data`")
  expect_error(aev(data, mortality = flat, settings = list(overdispersion = 2)),
               "must be a `settings`")
})

test_that("a mortality that underflows to zero returns rather than throwing", {
  # `exp()` REACHES EXACTLY ZERO at a log mortality of about -746, which is well
  # inside the range a table sentinel or a wandering fit would produce. E and V
  # are then zero while A still counts the deaths that happened -- the triple
  # `create_aev()` refuses, arrived at by arithmetic rather than by typing.
  res <- aev(data, mortality = -800, settings = basis)

  expect_true(is_aev(res))
  expect_equal(res$A, 2)
  expect_equal(res$E, 0)
  expect_identical(res$A_over_E, Inf)
})

test_that("one underflowing band does not take the other bands with it", {
  # THE COST OF THE OLD CHECK, MADE VISIBLE. Validation was vectorised over the
  # whole breakdown, so a single degenerate band discarded every other band's
  # answer alongside it.
  sentinel <- with_log_mu(ifelse(pensions >= 20000, -800, -4))

  res <- aev(sentinel,
             mortality = .i$log_mu,
             breakdown = bands(.i$pension, thresholds = c(10000, 20000)),
             settings  = basis)

  expect_equal(length(res), 3L)
  expect_true(all(is.finite(res$E[1:2])))
  expect_equal(res$E[[3L]], 0)
})

test_that("an underflowing element does not take the rest of a batch with it", {
  # `batch()` finalises in a plain loop with no error handling, so an error out of
  # any one element aborted every other element's already-computed answer.
  sentinel <- with_log_mu(rep(-800, length(pensions)))

  b <- batch(.exp_data = sentinel, .settings = basis,
             ordinary  = aev(mortality = -4),
             underflow = aev(mortality = .i$log_mu),
             also_fine = aev(mortality = -5))

  expect_named(b, c("ordinary", "underflow", "also_fine"))
  expect_true(all(vapply(b, is_aev, logical(1L))))
  expect_true(is.finite(b$ordinary$E))
  expect_equal(b$underflow$E, 0)
  expect_true(is.finite(b$also_fine$E))
})
