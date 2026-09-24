# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `ensure_is_valid_name()` guards every optional `name` in the UI. Its message
# promised "not blank and no leading, trailing or consecutive spaces" before the
# body checked any of it, so each promise gets its own witness here.

test_that("a valid name is a single non-blank string with tidy spacing", {
  expect_true(is_valid_name("S4PMA"))
  expect_true(is_valid_name("men, pensioners"))
  expect_true(is_valid_name("a"))
  expect_silent(ensure_is_valid_name("men, pensioners"))
})

test_that("a name may not be blank", {
  expect_false(is_valid_name(""))
  expect_false(is_valid_name(" "))
  expect_false(is_valid_name("   "))
  expect_error(ensure_is_valid_name(""), "must be a valid name")
})

test_that("a name may not have leading, trailing or consecutive spaces", {
  expect_false(is_valid_name(" men"))
  expect_false(is_valid_name("men "))
  expect_false(is_valid_name("men  pensioners"))
  expect_error(ensure_is_valid_name(" men"), "must be a valid name")
  expect_error(ensure_is_valid_name("men "), "must be a valid name")
  expect_error(ensure_is_valid_name("men  pensioners"), "must be a valid name")
})

test_that("a name may not hold a control character anywhere", {
  # An interior tab is neither leading, trailing nor consecutive, so it needs
  # its own rule: it would break the printing that lines names up in columns.
  expect_false(is_valid_name("men\tpensioners"))
  expect_false(is_valid_name("men\npensioners"))
  expect_false(is_valid_name("men\r"))
  expect_error(ensure_is_valid_name("men\tpensioners"), "must be a valid name")
})

test_that("a name must be a single plain string", {
  expect_false(is_valid_name(42))
  expect_false(is_valid_name(NA_character_))
  expect_false(is_valid_name(c("a", "b")))
  expect_false(is_valid_name(character(0)))
  expect_false(is_valid_name(factor("a")))
  expect_error(ensure_is_valid_name(42), "`42` must be a valid name")
})

test_that("the mortality constructors refuse a blank name", {
  # The exported path, so the guard is known to be wired in and not just to exist.
  expect_error(mortality_const(log_mu = -4, name = ""), "`name` must be a valid name")
  expect_error(gompertz_mortality(name = "men "), "`name` must be a valid name")
  expect_error(mortality_table(x0 = datey::durationy(60), t0 = datey::datey(2000),
                               log_mu = matrix(-4, 2, 2), name = "  "),
               "`name` must be a valid name")
  expect_silent(mortality_const(log_mu = -4, name = "flat"))
})
