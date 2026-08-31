# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("scalar max/min/clamp evaluate element-wise", {
  expect_equal(value(variable(max(.i$a, .i$b)), list(a = 1, b = 5)), 5)
  expect_equal(value(variable(min(.i$a, .i$b)), list(a = 1, b = 5)), 1)
  expect_equal(value(variable(clamp(.i$a, 0, 10)), list(a = 15)), 10)
  expect_equal(value(variable(clamp(.i$a, 0, 10)), list(a = -3)), 0)
  expect_equal(value(variable(clamp(.i$a, 0, 10)), list(a = 5)), 5)
})

test_that("idiomatic aliases are accepted and normalised", {
  expect_equal(it_deparse(variable(pmax(.i$a, .i$b))$ast), "max(.i$a, .i$b)")
  expect_equal(value(variable(pmin(.i$a, .i$b)), list(a = 1, b = 5)), 1)
  expect_true(value(variable(.i$p && .i$q), list(p = TRUE, q = TRUE)))
  expect_false(value(variable(.i$p && .i$q), list(p = TRUE, q = FALSE)))
  expect_true(value(variable(.i$p || .i$q), list(p = FALSE, q = TRUE)))
})

test_that("%in% with a character vector works", {
  pen <- variable(.i$category %in% c("ret", "dep"))
  expect_true(value(pen, list(category = "ret")))
  expect_false(value(pen, list(category = "act")))
  expect_equal(it_deparse(pen$ast), '(.i$category %in% c("ret", "dep"))')
})

test_that("namespaced functions over pronouns work", {
  expect_equal(value(variable(base::abs(.i$x)), list(x = -5)), 5)
  expect_equal(value(variable(abs(.i$x)), list(x = -5)), 5)
  expect_equal(value(variable(log10(.i$x)), list(x = 1000)), 3)
  expect_equal(value(variable(sqrt(.i$x)), list(x = 9)), 3)
})

test_that("ifelse selects on a scalar condition", {
  v <- variable(ifelse(.i$smoker, 100, 200))
  expect_equal(value(v, list(smoker = TRUE)), 100)
  expect_equal(value(v, list(smoker = FALSE)), 200)
})

test_that("single-expression braces are transparent", {
  expect_equal(it_deparse(variable({ .i$x + 1 })$ast), "(.i$x + 1)")
})

test_that("conditionals must have a time-invariant condition", {
  expect_s3_class(variable(ifelse(.i$smoker, 1, 0)), "variable")
  expect_error(variable(ifelse(.x < 65, 1, 0)), "must not depend on time")
  expect_error(variable(if (.t > .b) 1 else 0), "must not depend on time")
})

test_that("forbidden constructs are rejected", {
  expect_error(variable(z <- .i$x), "not allowed")
  expect_error(variable({ q <- 1; q + .i$x }), "multiple statements")
  expect_error(variable(.i$x + (function(z) z)(1)), "not allowed")
})

test_that("non-whitelisted functions and two-arg log are rejected", {
  expect_error(variable(weirdfun(.i$x)), "not an allowed")
  expect_error(variable(isTRUE(.i$flag)), "not an allowed")
  expect_error(variable(log(.i$x, 10)), "not supported")
})
