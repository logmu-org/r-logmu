# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("invalid exp_data and val_data", {

  t_1880 <- datey::datey(1880)
  t_1900 <- datey::datey(1900)
  t_1920 <- datey::datey(1920)
  t_1940 <- datey::datey(1940)
  t_1960 <- datey::datey(1960)
  t_1980 <- datey::datey(1980)
  t_2000 <- datey::datey(2000)
  t_2010 <- datey::datey(2010)
  t_2015 <- datey::datey(2015)
  t_2020 <- datey::datey(2020)
  t_2025 <- datey::datey(2025)
  t_2030 <- datey::datey(2030)

  birth_3 <- c(t_1880, t_1900, t_1920)
  birth_4 <- c(t_1940, t_1960, t_1980, t_2000)

  E2R_start_3 <- c(t_2010, t_2015, t_2020)
  E2R_end_3 <- c(t_2015, t_2020, t_2025)
  E2R_died_3 <- c(TRUE, FALSE, TRUE)

  E2R_start_4 <- c(t_2010, t_2015, t_2020, t_2010)
  E2R_end_4 <- c(t_2015, t_2020, t_2025, t_2015)
  E2R_died_4 <- c(TRUE, FALSE, TRUE, FALSE)

  double_len_3 <- c(1, 2, 3)
  double_len_4 <- c(4, 5, 6, 7)
  negative_3 <- c(1, -2, 3)
  NaN_3 <- c(1, NaN, 3)

  text_len_3 <- c("B", "C", "D")
  text_len_4 <- c("E", "F", "G", "H")

  exp_cols_3 <- list(
    birth = birth_3,
    pension = double_len_3,
    count = double_len_3,
    E2R_start = E2R_start_3, E2R_end = E2R_end_3, E2R_died = E2R_died_3,
    x = text_len_3
  )
  exp_cols_4 <- list(
    birth = birth_4,
    pension = double_len_4,
    E2R_start = E2R_start_4, E2R_end = E2R_end_4, E2R_died = E2R_died_4,
    x = text_len_4
  )

  val_cols_3 <- list(
    birth = birth_3,
    pension = double_len_3,
    count = double_len_3,
    x = text_len_3
  )
  val_cols_4 <- list(
    birth = birth_4,
    pension = double_len_4,
    x = text_len_4
  )

  # Check everything works first
  exp_data(exp_cols_3, exp_start = t_2000, exp_end = t_2025)
  val_data(val_cols_3, as_at= t_2025)
  exp_data(exp_cols_4, exp_start = t_2000, exp_end = t_2025)
  val_data(val_cols_4, as_at= t_2025)

  # No columns
  expect_error(exp_data(list(), exp_start = t_2000, exp_end = t_2025))
  expect_error(val_data(list(), as_at = t_2025))

  # NA attrs
  exp_cols_3_no_birth <- list(
    pension = double_len_3,
    E2R_start = E2R_start_3, E2R_end = E2R_end_3, E2R_died = E2R_died_3,
    x = text_len_3
  )
  exp_cols_3_no_E2R_start <- list(
    birth = birth_3,
    pension = double_len_3,
    E2R_end = E2R_end_3, E2R_died = E2R_died_3,
    x = text_len_3
  )
  exp_cols_3_E2R_end <- list(
    birth = birth_3,
    pension = double_len_3,
    E2R_start = E2R_start_3, E2R_died = E2R_died_3,
    x = text_len_3
  )
  exp_cols_3_no_E2R_died <- list(
    birth = birth_3,
    pension = double_len_3,
    E2R_start = E2R_start_3, E2R_end = E2R_end_3,
    x = text_len_3
  )
  val_cols_3_no_birth <- list(
    pension = double_len_3,
    x = text_len_3
  )
  val_cols_3_with_E2R_start <- list(
    birth = birth_3,
    pension = double_len_3,
    x = text_len_3,
    E2R_start_3
  )
  val_cols_3_with_E2R_end <- list(
    birth = birth_3,
    pension = double_len_3,
    x = text_len_3,
    E2R_end_3
  )
  val_cols_3_with_E2R_died <- list(
    birth = birth_3,
    pension = double_len_3,
    x = text_len_3,
    E2R_died_3
  )

  expect_error(exp_data(exp_cols_3, exp_start = NA, exp_end = t_2025))
  expect_error(exp_data(exp_cols_3, exp_start = "A", exp_end = t_2025))
  expect_error(exp_data(exp_cols_3, exp_start = datey::NA_datey_, exp_end = t_2025))
  expect_error(exp_data(exp_cols_3, exp_start = t_2000, exp_end = NA))
  expect_error(exp_data(exp_cols_3, exp_start = t_2000, exp_end = "A"))
  expect_error(exp_data(exp_cols_3, exp_start = t_2000, exp_end = datey::NA_datey_))
  expect_error(val_data(val_cols_3, as_at = NA))
  expect_error(val_data(val_cols_3, as_at = "A"))
  expect_error(val_data(val_cols_3, as_at = datey::NA_datey_))

  # Missing or extraneous columns
  expect_error(exp_data(exp_cols_3_no_birth, exp_start = t_2000, exp_end = t_2025))
  expect_error(exp_data(exp_cols_3_no_E2R_start, exp_start = t_2000, exp_end = t_2025))
  expect_error(exp_data(exp_cols_3_no_E2R_end, exp_start = t_2000, exp_end = t_2025))
  expect_error(exp_data(exp_cols_3_no_E2R_died, exp_start = t_2000, exp_end = t_2025))
  expect_error(val_data(val_cols_3_no_birth, as_at = t_2025))
  expect_error(val_data(val_cols_3_with_E2R_start, as_at = t_2025))
  expect_error(val_data(val_cols_3_with_E2R_end, as_at = t_2025))
  expect_error(val_data(val_cols_3_with_E2R_died, as_at = t_2025))

  # Incompatible column lengths
  exp_cols_3_and_4 <- list(
    birth = birth_3,
    pension = double_len_3,
    E2R_start = E2R_start_3, E2R_end = E2R_end_3, E2R_died = E2R_died_3,
    x = text_len_4
  )
  val_cols_3_and_4 <- list(
    birth = birth_4,
    pension = double_len_3,
    x = text_len_4
  )

  expect_error(exp_data(exp_cols_3_and_4, exp_start = t_2000, exp_end = t_2025))
  expect_error(val_data(val_cols_3_and_4, as_at = t_2025))

  # Invalid counts
  exp_cols_3_negative_count <- list(
    birth = birth_3,
    pension = double_len_3,
    E2R_start = E2R_start_3, E2R_end = E2R_end_3, E2R_died = E2R_died_3,
    x = text_len_4,
    count = negative_3
  )
  exp_cols_3_NaN_count <- list(
    birth = birth_3,
    pension = double_len_3,
    count = NaN_3,
    E2R_start = E2R_start_3, E2R_end = E2R_end_3, E2R_died = E2R_died_3,
    x = text_len_4
  )
  val_cols_3_negative_count <- list(
    birth = birth_4,
    count = negative_3,
    pension = double_len_3,
    x = text_len_4
  )
  val_cols_3_NaN_count <- list(
    birth = birth_4,
    pension = double_len_3,
    x = text_len_4,
    count = NaN_3
  )
  expect_error(exp_data(exp_cols_3_negative_count, exp_start = t_2000, exp_end = t_2025))
  expect_error(exp_data(exp_cols_3_NaN_count, exp_start = t_2000, exp_end = t_2025))
  expect_error(val_data(val_cols_3_negative_count, as_at = t_2025))
  expect_error(val_data(val_cols_3_NaN_count, as_at = t_2025))

  # Special names off by case
  exp_cols_special_name_off_case_A <- list(
    birth = birth_4,
    pension = double_len_4,
    E2r_start = E2R_start_4, E2R_end = E2R_end_4, E2R_died = E2R_died_4,
    x = text_len_4
  )
  exp_cols_special_name_off_case_B <- list(
    birth = birth_4,
    pension = double_len_4,
    UID = double_len_4,
    E2R_start = E2R_start_4, E2R_end = E2R_end_4, E2R_died = E2R_died_4,
    x = text_len_4
  )

  val_cols_special_name_off_case <- list(
    Birth = birth_4,
    pension = double_len_4,
    x = text_len_4
  )
  expect_error(exp_data(exp_cols_special_name_off_case_A, exp_start = t_2000, exp_end = t_2025))
  expect_error(exp_data(exp_cols_special_name_off_case_B, exp_start = t_2000, exp_end = t_2025))
  expect_error(val_data(val_cols_special_name_off_case, as_at = t_2025))

  # Illegal names (case insensitive)
  # ".data", ".birth", ".t`", ".age", "e2r"
  exp_cols_illegal_name.data <- list(
    birth = birth_4,
    .data = double_len_4,
    E2R_start = E2R_start_4, E2R_end = E2R_end_4, E2R_died = E2R_died_4,
    x = text_len_4
  )
  exp_cols_illegal_name_.BIRTH <- list(
    birth = birth_4,
    .BIRTH = double_len_4,
    uid = double_len_4,
    E2R_start = E2R_start_4, E2R_end = E2R_end_4, E2R_died = E2R_died_4,
    x = text_len_4
  )
  exp_cols_illegal_name_.E2R <- list(
    birth = birth_4,
    E2R = double_len_4,
    uid = double_len_4,
    E2R_start = E2R_start_4, E2R_end = E2R_end_4, E2R_died = E2R_died_4,
    x = text_len_4
  )

  val_cols_special_name_off_case.t <- list(
    birth = birth_4,
    .t = birth_4,
    x = text_len_4
  )
  val_cols_special_name_off_case.Age <- list(
    birth = birth_4,
    .Age = double_len_4,
    x = text_len_4
  )
  expect_error(exp_data(exp_cols_illegal_name.data,   exp_start = t_2000, exp_end = t_2025))
  expect_error(exp_data(exp_cols_illegal_name_.BIRTH, exp_start = t_2000, exp_end = t_2025))
  expect_error(exp_data(exp_cols_illegal_name_.E2R,   exp_start = t_2000, exp_end = t_2025))
  expect_error(val_data(val_cols_special_name_off_case.t,   as_at = t_2025))
  expect_error(val_data(val_cols_special_name_off_case.Age, as_at = t_2025))
})
