# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the time-varying tagging pass. A node is time-varying (a vector) when it reaches `.t` or a
# vector source; otherwise it is a per-individual scalar. A select's condition must be time-invariant,
# a rule the R layer already enforces before an expression can cross.

cols <- list(
  birth  = datey::datey(1960),
  sex    = factor("male"),
  salary = 30000,
  smoker = TRUE
)

root_vector <- function(res) res$vectors[[res$root + 1L]]

test_that("time and age are time-varying", {
  expect_true(root_vector(cpp_veil_prepare(it_ast(~ .t), cols)))
  expect_true(root_vector(cpp_veil_prepare(it_ast(~ .x), cols)))
})

test_that("a field and a fact comparison are time-invariant", {
  expect_false(root_vector(cpp_veil_prepare(it_ast(~ .i$salary), cols)))
  expect_false(root_vector(cpp_veil_prepare(it_ast(~ .i$sex == "male"), cols)))
})

test_that("a mortality table is a time-varying vector", {
  log_mu <- matrix(log(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06)), nrow = 3, ncol = 2)
  tbl <- mortality_table(x0 = 60, t0 = 2010, log_mu = log_mu)
  expect_true(root_vector(cpp_veil_prepare(it_obj(tbl), cols)))
})

test_that("time-dependence propagates through arithmetic but a fact expression does not", {
  expect_true(root_vector(cpp_veil_prepare(it_ast(~ -10 + 0.1 * .x), cols)))
  expect_false(root_vector(cpp_veil_prepare(it_ast(~ .i$salary * 2 + 1), cols)))
})

test_that("a select is time-varying through a branch, invariant through its facts", {
  expect_true(root_vector(cpp_veil_prepare(it_ast(~ ifelse(.i$smoker, .x, 0.5)), cols)))
  expect_false(root_vector(cpp_veil_prepare(it_ast(~ ifelse(.i$smoker, .i$salary, 0)), cols)))
})

test_that("a time-dependent condition is refused before it can cross", {
  # it_check_conditional enforces this at the R level, so veil's own backstop is reached only by a
  # front end that does not do the check itself.
  expect_error(it_ast(~ if (.x < 65) 1 else 2), "depend on time")
})
