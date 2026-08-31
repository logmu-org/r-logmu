# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("a referenced variable is spliced by its AST", {
  base <- variable(.i$pension)
  v <- variable(base * 2)
  expect_equal(it_deparse(v$ast), "(.i$pension * 2)")
  expect_equal(value(v, .i = list(pension = 100)), 200)
})

test_that("a referenced indicator is spliced into a variable", {
  male <- indicator(.i$sex == "male")
  v <- variable(.i$pension * male)
  expect_equal(value(v, .i = list(pension = 100, sex = "male")), 100)
  expect_equal(value(v, .i = list(pension = 100, sex = "female")), 0)
})

test_that("a referenced formula is spliced by its body", {
  f <- ~ .i$sex == "male"
  male <- indicator(f)
  expect_true(is_indicator(male))
  expect_true(logical_value(male, .i = list(sex = "male")))
})

test_that("an ordinary scalar reference still folds", {
  k <- 5
  v <- variable(.i$pension * k)
  expect_equal(it_deparse(v$ast), "(.i$pension * 5)")
})

test_that("a referenced mortality becomes an opaque obj leaf", {
  mort <- mortality_const(q = 0.01)
  node <- it_build_ast(quote(mort), environment())
  expect_equal(node$kind, "obj")
  expect_identical(node$value, mort)
  expect_equal(it_deparse(node), "<mortality_const>")
})

test_that("mortality selection splices conditions and obj leaves", {
  male   <- indicator(.i$sex == "male")
  mort_m <- mortality_const(q = 0.01)
  mort_f <- mortality_const(q = 0.02)

  sel <- it_build_ast(quote(if (male) mort_m else mort_f), environment())

  expect_equal(sel$fn, "if")
  expect_equal(sel$args[[2L]]$kind, "obj")
  expect_equal(it_deparse(sel), 'if ((.i$sex == "male")) <mortality_const> else <mortality_const>')

  # the reference evaluator selects the right mortality per individual
  expect_identical(it_eval(sel, .i = list(sex = "male")),   mort_m)
  expect_identical(it_eval(sel, .i = list(sex = "female")), mort_f)
})
