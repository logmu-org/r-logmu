# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("val_data can be populated by a data.frame", {

  b_1 <- datey::datey(1951)
  b_2 <- datey::datey(1952)
  b_3 <- datey::datey(1953)

  birth <- c(b_1, b_2, b_3, b_1, b_2, b_3)
  count <- c(10, 20, 30)

  count2 <- c(10, 20, 30, 10, 20, 30)

  df <- data.frame(birth, count)

  as_at <- datey::datey(2025)

  val_data <- val_data(df, as_at)

  expect_identical(val_data$birth, birth)
  expect_identical(val_data$count, count2)
})
