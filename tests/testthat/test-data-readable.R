# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("all recognised column types are accepted", {
  b <- datey::datey(1950)
  cols <- list(
    birth = c(b, b),
    dbl   = c(1.0, 2.0),
    int   = c(1L, 2L),
    lgl   = c(TRUE, FALSE),
    chr   = c("a", "b"),
    fct   = factor(c("x", "y")),
    dur   = datey::durationy(c(1, 2))
  )
  expect_s3_class(val_data(cols, as_at = datey::datey(2025)), "val_data")
})

test_that("unreadable column types follow the on_unreadable policy", {
  b <- datey::datey(1950)
  cols <- list(birth = c(b, b), when = as.Date(c("2020-01-01", "2021-01-01")))

  expect_error(val_data(cols, as_at = datey::datey(2025)), "cannot read")
  expect_warning(val_data(cols, as_at = datey::datey(2025), on_unreadable = "warn"), "cannot read")
  expect_s3_class(
    suppressWarnings(val_data(cols, as_at = datey::datey(2025), on_unreadable = "ignore")),
    "val_data"
  )
})

test_that("the error names the column and its type", {
  b <- datey::datey(1950)
  cols <- list(birth = c(b, b), when = as.Date(c("2020-01-01", "2021-01-01")))
  expect_error(val_data(cols, as_at = datey::datey(2025)), "when.*Date")
})

test_that("assigning an unreadable column errors", {
  b <- datey::datey(1950)
  v <- val_data(list(birth = c(b, b), x = c(1, 2)), as_at = datey::datey(2025))
  expect_error({ v$bad <- as.Date(c("2020-01-01", "2021-01-01")) }, "cannot read")
})
