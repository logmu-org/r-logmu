# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the type annotation pass. Every node is given a base type by the datey rules (a number makes
# a value operator double, pure datey/durationy uses the exact integer algebra), and an ill-typed
# expression is rejected. Runs through cpp_veil_prepare, which reads each column's type.

cols <- list(
  birth  = datey::datey(1960),
  sex    = factor("male"),
  salary = 30000,
  smoker = TRUE
)

root_type <- function(res) res$types[[res$root + 1L]]

test_that("time is a datey and age is a durationy", {
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ .t), cols)), "datey")
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ .x), cols)), "durationy") # .t - .b
})

test_that("a number makes a value operator double", {
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ -10 + 0.1 * .x), cols)), "double")
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ clamp((105 - .x) / 40, 0, 1)), cols)), "double")
})

test_that("comparisons and logicals are bool", {
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ .i$sex == "male"), cols)), "bool")
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ .i$smoker & (.x > 65)), cols)), "bool")
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ is.na(.i$salary)), cols)), "bool")
})

# Text records WHICH value, not where it sits in a sequence, so ordering it is meaningless. Refused
# on type, in the annotation pass, rather than left to lowering -- which is why these run through
# cpp_veil_prepare and not through a calculation.
test_that("text can be compared for equality but not ordered", {
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ .i$sex == "male"), cols)), "bool")
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ .i$sex != "male"), cols)), "bool")

  expect_error(cpp_veil_prepare(it_ast(~ .i$sex <  "male"), cols), "cannot order")
  expect_error(cpp_veil_prepare(it_ast(~ .i$sex <= "male"), cols), "cannot order")
  expect_error(cpp_veil_prepare(it_ast(~ .i$sex >  "male"), cols), "cannot order")
  expect_error(cpp_veil_prepare(it_ast(~ .i$sex >= "male"), cols), "cannot order")

  # Reversed, so the refusal follows the operand type rather than which side it sits on.
  expect_error(cpp_veil_prepare(it_ast(~ "male" < .i$sex), cols), "cannot order")

  # Two text fields, so it is not the literal that is being refused.
  two <- c(cols, list(scheme = factor("A")))
  expect_error(cpp_veil_prepare(it_ast(~ .i$sex < .i$scheme), two), "cannot order")
  expect_equal(root_type(cpp_veil_prepare(it_ast(~ .i$sex == .i$scheme), two)), "bool")
})

test_that("a mortality table's vector_log_mu is a double vector", {
  log_mu <- matrix(log(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)), nrow = 3, ncol = 2)
  tbl <- mortality_table(x0 = 60, t0 = 2010, log_mu = log_mu)
  res <- cpp_veil_prepare(it_obj(tbl), cols)

  expect_equal(root_type(res), "double") # the vector_log_mu node
  expect_true("" %in% res$types)         # the table obj carries no value type
})

test_that("an include is a datey_interval", {
  expect_equal(root_type(cpp_veil_prepare(it_obj(age(65, 95)), cols)), "datey_interval")
})

test_that("datey plus datey is rejected", {
  expect_error(cpp_veil_prepare(it_ast(~ .b + .b), cols), "add")
})

test_that("a text field cannot enter a numeric operator", {
  expect_error(cpp_veil_prepare(it_ast(~ exp(.i$sex)), cols), "numeric")
})

test_that("select branches must have a common type", {
  expect_error(cpp_veil_prepare(it_ast(~ ifelse(.i$smoker, .x, .i$sex)), cols), "branches")
})

test_that("a bare integer column is rejected at prepare time", {
  wide <- c(cols, list(n = 1:2))
  expect_error(cpp_veil_prepare(it_ast(~ .i$n > 0), wide), "bare integer")
})
