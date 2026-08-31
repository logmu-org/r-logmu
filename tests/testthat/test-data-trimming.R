# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("exp_data trims by E2R", {
  t_0 <- datey::datey(2000)
  t_2 <- datey::datey(2002)
  t_4 <- datey::datey(2004)
  t_6 <- datey::datey(2006)
  t_8 <- datey::datey(2008)

  b <- datey::datey(1950)


  birth <- c(b, b, b, b, b, b, b, b, b, b)
  count <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
  E2R_start <- c(t_0, t_0, t_2, t_0, t_2, t_4, t_0, t_2, t_4, t_6)
  E2R_end   <- c(t_2, t_4, t_4, t_6, t_6, t_6, t_8, t_8, t_8, t_8)
  E2R_died  <- c(  T,   F,   T,   F,   T,   F,   T,   F,   T,   F)

  df <- data.frame(birth, count, E2R_start, E2R_end, E2R_died)


  # All E2Rs
  exp_data <- exp_data(df, exp_start = t_0, exp_end = t_8)

  expect_identical(exp_data$birth, birth)
  expect_identical(exp_data$count, count)
  expect_identical(exp_data$E2R_start, E2R_start)
  expect_identical(exp_data$E2R_end, E2R_end)
  expect_identical(exp_data$E2R_died, E2R_died)


  # Trim left
  exp_data <- exp_data(df, exp_start = t_4, exp_end = t_8)

  birth <- c(b, b, b, b, b, b, b)
  count <- c(4, 5, 6, 7, 8, 9, 10)
  E2R_start <- c(t_4, t_4, t_4, t_4, t_4, t_4, t_6)
  E2R_end   <- c(t_6, t_6, t_6, t_8, t_8, t_8, t_8)
  E2R_died  <- c(  F,   T,   F,   T,   F,   T,   F)

  expect_identical(exp_data$birth, birth)
  expect_identical(exp_data$count, count)
  expect_identical(exp_data$E2R_start, E2R_start)
  expect_identical(exp_data$E2R_end, E2R_end)
  expect_identical(exp_data$E2R_died, E2R_died)


  # Trim right
  exp_data <- exp_data(df, exp_start = t_0, exp_end = t_4)

  birth <- c(b, b, b, b, b,    b, b)
  count <- c(1, 2, 3, 4, 5,    7, 8)
  E2R_start <- c(t_0, t_0, t_2, t_0, t_2,      t_0, t_2)
  E2R_end   <- c(t_2, t_4, t_4, t_4, t_4,      t_4, t_4)
  E2R_died  <- c(  T,   F,   T,   F,   F,        F,   F)

  expect_identical(exp_data$birth, birth)
  expect_identical(exp_data$count, count)
  expect_identical(exp_data$E2R_start, E2R_start)
  expect_identical(exp_data$E2R_end, E2R_end)
  expect_identical(exp_data$E2R_died, E2R_died)


  # Trim both ends
  exp_data <- exp_data(df, exp_start = t_2, exp_end = t_6)

  birth <- c(   b, b, b, b, b, b, b, b   )
  count <- c(   2, 3, 4, 5, 6, 7, 8, 9   )
  E2R_start <- c(     t_2, t_2, t_2, t_2, t_4, t_2, t_2, t_4     )
  E2R_end   <- c(     t_4, t_4, t_6, t_6, t_6, t_6, t_6, t_6     )
  E2R_died  <- c(       F,   T,   F,   T,   F,   F,   F,   F     )

  expect_identical(exp_data$birth, birth)
  expect_identical(exp_data$count, count)
  expect_identical(exp_data$E2R_start, E2R_start)
  expect_identical(exp_data$E2R_end, E2R_end)
  expect_identical(exp_data$E2R_died, E2R_died)
})

test_that("exp_data trims by count", {
  t_0 <- datey::datey(2000)
  t_2 <- datey::datey(2002)
  t_4 <- datey::datey(2004)
  t_6 <- datey::datey(2006)
  t_8 <- datey::datey(2008)

  b <- datey::datey(1950)


  birth <- c(b, b, b, b, b, b, b, b, b, b)
  count <- c(1, 2, 0, 4, 5, 6, 7, 0, 9, 10)
  E2R_start <- c(t_0, t_0, t_2, t_0, t_2, t_4, t_0, t_2, t_4, t_6)
  E2R_end   <- c(t_2, t_4, t_4, t_6, t_6, t_6, t_8, t_8, t_8, t_8)
  E2R_died  <- c(  T,   F,   T,   F,   T,   F,   T,   F,   T,   F)

  df <- data.frame(birth, count, E2R_start, E2R_end, E2R_died)


  exp_data <- exp_data(df, exp_start = t_2, exp_end = t_6)

  birth <- c(   b,    b, b, b, b,    b   )
  count <- c(   2,    4, 5, 6, 7,    9   )
  E2R_start <- c(     t_2,      t_2, t_2, t_4, t_2,      t_4     )
  E2R_end   <- c(     t_4,      t_6, t_6, t_6, t_6,      t_6     )
  E2R_died  <- c(       F,        F,   T,   F,   F,        F     )

  expect_identical(exp_data$birth, birth)
  expect_identical(exp_data$count, count)
  expect_identical(exp_data$E2R_start, E2R_start)
  expect_identical(exp_data$E2R_end, E2R_end)
  expect_identical(exp_data$E2R_died, E2R_died)
})

test_that("val_data trimms by count", {
  b <- datey::datey(1950)

  birth   <- c( b,  b,  b,  b,  b,  b,  b,  b,  b,  b)
  count   <- c( 1,  2,  3,  4,  5,  6,  7,  8,  9, 10)
  count_z <- c( 1,  2,  0,  4,  5,  6,  7,  0,  9, 10)
  x       <- c(11, 12, 13, 14, 15, 16, 17, 18, 19, 20)

  birth2 <- c( b,  b,      b,  b,  b,  b,      b,  b)
  count2 <- c( 1,  2,      4,  5,  6,  7,      9, 10)
  x2     <- c(11, 12,     14, 15, 16, 17,     19, 20)

  df <- data.frame(birth, count, x)
  val_data <- val_data(df, as_at = datey::datey(2025))
  val_data$count <- count_z
  expect_identical(val_data$birth, birth2)
  expect_identical(val_data$count, count2)
  expect_identical(val_data$x, x2)

  df <- data.frame(birth, count = count_z, x)
  val_data <- val_data(df, as_at = datey::datey(2025))
  expect_identical(val_data$birth, birth2)
  expect_identical(val_data$count, count2)
  expect_identical(val_data$x, x2)
})
