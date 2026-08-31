# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the coercion insertion pass. A value operator that widened to double gets its datey/durationy
# or bool operands wrapped in a ToDouble node; comparisons and pure datey/durationy arithmetic are
# left alone. A ToDouble node is a call whose moniker is `double`.

cols <- list(
  birth  = datey::datey(1960),
  sex    = "male",
  salary = 30000,
  smoker = TRUE
)

n_coerce <- function(res) sum(res$kinds == "call" & res$labels == "double")

test_that("a number coerces a durationy operand to double", {
  # 0.1 * .x : the durationy .x is widened for the multiply.
  expect_equal(n_coerce(cpp_veil_prepare(it_ast(~ 0.1 * .x), cols)), 1L)
})

test_that("a transcendental coerces its operand", {
  expect_equal(n_coerce(cpp_veil_prepare(it_ast(~ exp(.x)), cols)), 1L)
})

test_that("pure datey/durationy arithmetic is not coerced", {
  # .t - .b is datey - datey -> durationy, the integer algebra, so nothing widens.
  expect_equal(n_coerce(cpp_veil_prepare(it_ast(~ .x), cols)), 0L)
})

test_that("an all-double expression needs no coercion", {
  expect_equal(n_coerce(cpp_veil_prepare(it_ast(~ .i$salary * 2 + 1), cols)), 0L)
})

test_that("a select widens a branch but not the condition", {
  # ifelse(smoker, .x, 0.5): the .x branch is widened; the bool condition is left as it is.
  expect_equal(n_coerce(cpp_veil_prepare(it_ast(~ ifelse(.i$smoker, .x, 0.5)), cols)), 1L)
})

test_that("a comparison is left for the narrowing pass", {
  expect_equal(n_coerce(cpp_veil_prepare(it_ast(~ .x > 65), cols)), 0L)
})
