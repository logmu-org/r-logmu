# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the column scan: the reader materialises R columns into a ColumnSet, and the R-free scan
# reads off each column's NAs, whether its values are constant, and its per-type range. Text is not a
# column type -- a column of strings arrives as a factor and is read as a category.

test_that("a plain double column reports its range and integrality", {
  res <- cpp_veil_scan(list(age = c(30, 45, 60, 45)))
  i <- which(res$names == "age")

  expect_equal(res$types[i], "double")
  expect_false(res$has_nas[i])
  expect_false(res$constant[i])
  expect_equal(res$min[i], 30)
  expect_equal(res$max[i], 60)
  expect_true(res$all_integral[i])
})

test_that("a double column with an NA and a fraction is flagged", {
  res <- cpp_veil_scan(list(wt = c(1.5, 2.0, NA, 3.5)))
  i <- which(res$names == "wt")

  expect_true(res$has_nas[i])
  expect_false(res$all_integral[i])
  expect_equal(res$min[i], 1.5)
  expect_equal(res$max[i], 3.5)
})

test_that("a datey column is typed and scanned", {
  res <- cpp_veil_scan(list(birth = datey::datey(c(1960, 1970, 1960, 1980))))
  i <- which(res$names == "birth")

  expect_equal(res$types[i], "datey")
  expect_false(res$has_nas[i])
  expect_false(res$constant[i])
  expect_lt(res$min[i], res$max[i])
})

test_that("a mixed logical column spans 0 to 1", {
  res <- cpp_veil_scan(list(ind = c(TRUE, FALSE, TRUE, TRUE)))
  i <- which(res$names == "ind")

  expect_equal(res$types[i], "bool")
  expect_false(res$constant[i])
  expect_equal(res$min[i], 0)
  expect_equal(res$max[i], 1)
})

test_that("a constant column is recognised", {
  res <- cpp_veil_scan(list(one = c(5, 5, 5, 5)))
  i <- which(res$names == "one")

  expect_true(res$constant[i])
  expect_equal(res$min[i], res$max[i])
})

test_that("a missing logical value is refused", {
  expect_error(cpp_veil_scan(list(x = c(TRUE, NA, FALSE))), "missing logical")
})

test_that("a character column is refused, because text reaches veil as a factor", {
  expect_error(cpp_veil_scan(list(s = c("a", "b"))), "must be a factor")
})

# A factor is scanned as a category: its codes are remapped onto the crossing's shared numbering, so
# the range reported is in that numbering rather than in the factor's own.
test_that("a factor column scans as a category over its codes", {
  res <- cpp_veil_scan(list(sex = factor(c("male", "female", "male"))))
  i <- which(res$names == "sex")

  expect_equal(res$types[i], "category")
  expect_false(res$has_nas[i])
  expect_false(res$constant[i])
  expect_equal(res$min[i], 0)
  expect_equal(res$max[i], 1)
})

test_that("a constant factor column is seen to be constant", {
  res <- cpp_veil_scan(list(sex = factor(c("male", "male"))))
  i <- which(res$names == "sex")

  expect_true(res$constant[i])
  expect_false(res$has_nas[i])
})

# A category has no missing value: its storage is `int`, which carries no platform-independent NA, so
# an NA factor is refused where it is read rather than marked and carried. Only double, datey and
# durationy have a missing state. A category column therefore never reports `has_nas`.
test_that("a factor column carrying an NA is refused when it is read", {
  expect_error(cpp_veil_scan(list(sex = factor(c("male", NA, "female")))), "must not contain NA")
  expect_error(cpp_veil_scan(list(sex = factor(c(NA, NA), levels = c("male", "female")))),
               "must not contain NA")
})
