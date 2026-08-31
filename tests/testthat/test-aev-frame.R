# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `as.data.frame.aev()` is the first layer of charting: everything drawn is read
# off this frame, so its column set and its behaviour on degenerate records are
# both load-bearing.
#
# THE SHAPE IS FIXED BY THE TYPE AND NOT BY THE VALUES. `name` and `group` are
# present whether or not the `aev` was labelled, which is what lets a frame from
# a breakdown and a frame from an ungrouped analysis be stacked. The tests below
# assert the whole column vector by identity rather than checking membership,
# because a column quietly appearing or vanishing is the failure being guarded
# against.

frame_columns <- c("name", "group", "A", "E", "V",
                   "A_over_E", "log_A_over_E_stddev", "deviance_residual")

labelled_aev <- function() {
  aev <- create_aev(A = c(1100, 40, 0), E = c(1000, 50, 20), V = c(2500, 125, 40))
  names(aev) <- c("65-70", "70-75", "75-80")
  group_names(aev) <- "age"
  aev
}

test_that("the frame has exactly the charting columns, in order", {

  frame <- as.data.frame(labelled_aev())

  expect_s3_class(frame, "data.frame")
  expect_identical(names(frame), frame_columns)
  expect_identical(nrow(frame), 3L)
})

test_that("one row per record, however many there are", {

  for (records in c(1L, 3L, 10L)) {
    aev <- create_aev(A = rep(1, records), E = rep(1, records), V = rep(1, records))
    expect_identical(nrow(as.data.frame(aev)), records)
    expect_identical(nrow(as.data.frame(aev)), length(aev))
  }
})

test_that("a zero-record aev gives a zero-row frame with the same columns", {

  frame <- as.data.frame(labelled_aev()[integer(0)])

  expect_identical(nrow(frame), 0L)
  expect_identical(names(frame), frame_columns)
})

test_that("the columns hold what the accessors hold", {

  aev <- labelled_aev()
  frame <- as.data.frame(aev)

  expect_identical(frame$A, aev$A)
  expect_identical(frame$E, aev$E)
  expect_identical(frame$V, aev$V)
  expect_identical(frame$A_over_E, aev$A_over_E)
  expect_identical(frame$log_A_over_E_stddev, aev$log_A_over_E_stddev)
  expect_identical(frame$deviance_residual, aev$deviance_residual)
})

test_that("labels arrive as character columns, not factors", {

  frame <- as.data.frame(labelled_aev())

  expect_identical(frame$name, c("65-70", "70-75", "75-80"))
  expect_identical(frame$group, rep("age", 3))
  expect_type(frame$name, "character")
  expect_type(frame$group, "character")
})

test_that("an unlabelled aev keeps the columns and fills them with NA", {

  frame <- as.data.frame(create_aev(A = 1100, E = 1000, V = 2500))

  expect_identical(names(frame), frame_columns)
  expect_identical(frame$name, NA_character_)
  expect_identical(frame$group, NA_character_)
})

test_that("frames from a labelled and an unlabelled aev stack", {

  # The reason the two label columns are unconditional. If either were dropped
  # when absent, this `rbind()` would fail on differing column sets.
  stacked <- rbind(as.data.frame(labelled_aev()),
                   as.data.frame(create_aev(A = 7, E = 10, V = 20)))

  expect_identical(nrow(stacked), 4L)
  expect_identical(names(stacked), frame_columns)
  expect_identical(stacked$name, c("65-70", "70-75", "75-80", NA))
})

test_that("names() alone, without group_names(), still gives both columns", {

  aev <- create_aev(A = c(1, 2), E = c(1, 2), V = c(1, 2))
  names(aev) <- c("a", "b")
  frame <- as.data.frame(aev)

  expect_identical(frame$name, c("a", "b"))
  expect_identical(frame$group, rep(NA_character_, 2))
})

test_that("row.names is honoured and defaults to the usual integers", {

  aev <- labelled_aev()

  expect_identical(rownames(as.data.frame(aev)), c("1", "2", "3"))
  expect_identical(rownames(as.data.frame(aev, row.names = c("x", "y", "z"))),
                   c("x", "y", "z"))
})

test_that("dispatch goes through the generic", {

  # `as.data.frame()` on an aev must find the method. Calling the method
  # directly would pass even with the S3 registration missing from NAMESPACE.
  expect_identical(as.data.frame(labelled_aev()),
                   as.data.frame.aev(labelled_aev()))
  expect_identical(class(as.data.frame(labelled_aev())), "data.frame")
})

test_that("empty and missing records are indistinguishable in the derived columns", {

  # THE REASON `autoplot()` HAS TO READ THE TRIPLE. An empty cell and a missing
  # one agree in every calculated column, so a chart deciding what to draw from
  # `A_over_E` or the residual alone cannot tell a true statement about nobody
  # from an absence of information.
  empty <- as.data.frame(create_aev(A = 0, E = 0, V = 0))
  missing <- as.data.frame(create_aev_unchecked(A = NaN, E = NaN, V = NaN))

  for (column in c("A_over_E", "log_A_over_E_stddev", "deviance_residual")) {
    expect_true(is.nan(empty[[column]]))
    expect_true(is.nan(missing[[column]]))
  }

  # They differ in the triple, which is where the distinction survives.
  expect_identical(empty$E, 0)
  expect_true(is.nan(missing$E))
})

test_that("a zero-death record is off the log scale but on the residual scale", {

  # A = 0 with exposure is a real answer, not a degenerate one: the residual
  # panel can draw it and the A/E panel cannot, since log(0) has no position.
  frame <- as.data.frame(create_aev(A = 0, E = 50, V = 125))

  expect_identical(frame$A_over_E, 0)
  expect_identical(log(frame$A_over_E), -Inf)
  expect_true(is.finite(frame$deviance_residual))
  expect_true(frame$deviance_residual < 0)
})

test_that("labels from a breakdown reach the frame", {

  data <- exp_data(
    list(
      birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
      pension   = c(5000, 12000, 30000, 8000, 15000),
      E2R_start = datey::datey(rep(2015, 5)),
      E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
      E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )

  result <- aev(data,
                mortality      = mortality_const(log_mu = -4),
                breakdown      = bands(.i$pension, thresholds = c(10000, 20000)),
                overdispersion = 1)

  frame <- as.data.frame(result)

  expect_identical(nrow(frame), length(result))
  expect_identical(frame$name, names(result))
  expect_identical(frame$group, group_names(result))
  expect_false(anyNA(frame$name))
})
