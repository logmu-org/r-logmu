# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("val_from_exp_data", {
  t_0 <- datey::datey(2000)
  t_2 <- datey::datey(2002)
  t_4 <- datey::datey(2004)
  t_6 <- datey::datey(2006)
  t_8 <- datey::datey(2008)

  b <- datey::datey(1950)
  b2 <- datey::datey(1960)

  birth     <- c(  b,  b2,   b,  b2,   b,  b2,   b,  b2,   b,  b2)
  count     <- c(  1,   2,   3,   4,   5,   6,   7,   8,   9,  10)
  birth4    <- c(      b2,       b2,   b,  b2,   b,  b2,   b     )
  count4    <- c(       2,        4,   5,   6,   7,   8,   9     )
  birth8     <- c(                                   b2,       b2)
  count8     <- c(                                    8,       10)
  E2R_start <- c(t_0, t_0, t_2, t_0, t_2, t_4, t_0, t_2, t_4, t_6)
  E2R_end   <- c(t_2, t_4, t_4, t_6, t_6, t_6, t_8, t_8, t_8, t_8)
  E2R_died  <- c(  T,   F,   T,   F,   T,   F,   T,   F,   T,   F)

  df <- data.frame(birth, count, E2R_start, E2R_end, E2R_died)
  exp_data <- exp_data(df, exp_start = t_0, exp_end = t_8)

  val_data4 <- val_from_exp_data(exp_data, as_at = t_4)
  expect_identical(val_data4$birth, birth4)
  expect_identical(val_data4$count, count4)

  val_data8 <- val_from_exp_data(exp_data) # Default as_at
  expect_identical(val_data8$birth, birth8)
  expect_identical(val_data8$count, count8)
})
