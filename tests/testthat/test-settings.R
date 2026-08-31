# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `settings` carries the two choices an analysis needs that are not part of the
# question being asked. Both have a trap in them:
#
#   - overdispersion is MANDATORY with no default anywhere, because a default is
#     an invitation to ignore it and ignoring it understates every confidence
#     interval by sqrt(Omega) without looking wrong;
#   - a time scale is checked as WHOLE CLICKS, which is what lets 1/12 be
#     accepted exactly despite being unrepresentable as a double.

clicks_per_year <- unclass(datey::durationy(1))

test_that("overdispersion is required and has no default", {
  # The error names the argument rather than complaining about a missing value
  # somewhere downstream, since "which argument" is the whole content of it.
  expect_error(settings(), "`overdispersion` is required")
  expect_error(settings(time_scale = 1 / 12), "`overdispersion` is required")
})

test_that("overdispersion must be a single positive number", {
  expect_error(settings(overdispersion = 0), "positive")
  expect_error(settings(overdispersion = -1), "positive")
  expect_error(settings(overdispersion = NaN), "positive")
  expect_error(settings(overdispersion = NA_real_), "positive")
  expect_error(settings(overdispersion = Inf), "positive")
  expect_error(settings(overdispersion = c(2, 3)), "single")
  expect_error(settings(overdispersion = "2"), "positive")

  expect_equal(settings(overdispersion = 2)$overdispersion, 2)
  expect_equal(settings(overdispersion = 2.5)$overdispersion, 2.5)
})

test_that("the four permitted time scales are accepted, as numbers or durationy", {
  for (per_year in c(1, 4, 12, 60)) {
    expected <- clicks_per_year %/% per_year

    expect_equal(settings(2, time_scale = 1 / per_year)$time_scale_clicks, expected)
    expect_equal(
      settings(2, time_scale = datey::durationy(1 / per_year))$time_scale_clicks,
      expected)
  }
})

test_that("a month lands on a whole click, so the check is exact integer work", {
  # THE PROPERTY THE WHOLE CHECK RESTS ON. A month is a whole number of clicks
  # and `durationy()` gets there from a double, so the comparison against the
  # permitted set is integer and needs no tolerance anywhere. If datey ever stops
  # storing whole clicks, this fails before anything subtler does.
  expect_identical(settings(2, time_scale = 1 / 12)$time_scale_clicks, 44530L)
  expect_identical(clicks_per_year %% 44530L, 0L)
  expect_identical(44530L %% 2L, 0L) # The half-interval too -- midpoints need it.
})

test_that("a truncated fraction is forgiven to half a click and no further", {
  # A consequence of routing a double through `durationy()` rather than a choice.
  # It is worth pinning because it is the only tolerance anywhere in this path:
  # about half a minute, so a user writing out a fraction by hand is fine and a
  # genuinely different number is not.
  expect_identical(settings(2, time_scale = 0.083333)$time_scale_clicks, 44530L)
  expect_error(settings(2, time_scale = 0.08333), "must be one of")
})

test_that("the default time scale is a quarter", {
  expect_equal(settings(overdispersion = 2)$time_scale_clicks, clicks_per_year %/% 4L)
})

test_that("a time scale outside the permitted set is refused", {
  # Half a year and a day both satisfy the whole-clicks constraint the engine
  # asserts, so they are excluded by the CHOSEN SET rather than by arithmetic --
  # which is exactly why they need refusing explicitly.
  expect_error(settings(2, time_scale = 0.5), "must be one of")
  expect_error(settings(2, time_scale = 1 / 365), "must be one of")
  expect_error(settings(2, time_scale = 1 / 2), "must be one of")
  expect_error(settings(2, time_scale = 0), "must be one of")

  # A number of INTERVALS is the old spelling, and 4 years is not what anyone
  # writing it meant.
  expect_error(settings(2, time_scale = 4), "must be one of")

  expect_error(settings(2, time_scale = "quarterly"), "durationy")
  expect_error(settings(2, time_scale = c(0.25, 1)), "durationy")
})

test_that("the error names the values that are permitted", {
  # Derived from the engine's own list, so the message cannot describe a set that
  # is no longer the set being enforced.
  expect_error(settings(2, time_scale = 0.5), "1, 1/4, 1/12, 1/60", fixed = TRUE)
})

test_that("the permitted set comes from the engine, not from a copy in R", {
  # A second spelling in R is somewhere for the two to drift, and the failure it
  # hides is one-sided and silent.
  expect_identical(
    sort(cpp_veil_time_scales()),
    sort(as.integer(clicks_per_year %/% c(1L, 4L, 12L, 60L))))
})

test_that("settings are immutable", {
  s <- settings(overdispersion = 2)
  expect_error(s$overdispersion <- 3, "invalid")
  expect_error(s[["time_scale_clicks"]] <- 1L, "invalid")
})

test_that("settings print both choices", {
  expect_output(print(settings(2, time_scale = 1 / 12)), "overdispersion: 2")
  expect_output(print(settings(2, time_scale = 1 / 12)), "1/12 year")
  expect_output(print(settings(2, time_scale = 1)), "1 year")
})

test_that("is_settings recognises one and only one thing", {
  expect_true(is_settings(settings(overdispersion = 2)))
  expect_false(is_settings(list(overdispersion = 2)))
  expect_false(is_settings(2))
  expect_false(is_settings(NULL))
})
