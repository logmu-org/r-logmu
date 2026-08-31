# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("pronouns parse to the expected node kinds", {
  expect_equal(it_build_ast(quote(.t))$kind, "time")
  expect_equal(it_build_ast(quote(.b)), it_field("birth"))
  expect_equal(it_build_ast(quote(.i$pension)), it_field("pension"))

  # .x is sugar for .t - .b
  expect_equal(it_build_ast(quote(.x)), it_build_ast(quote(.t - .b)))
  expect_equal(it_deparse(it_build_ast(quote(.x))), "(.t - .i$birth)")
})

test_that("redundant parentheses are dropped", {
  expect_equal(it_build_ast(quote((.t - .b))), it_build_ast(quote(.t - .b)))
  expect_equal(
    it_deparse(it_build_ast(quote(.i$pension * (.t - .b)))),
    "(.i$pension * (.t - .i$birth))"
  )
})

test_that("a comparison keeps its structure and folds the literal", {
  node <- it_build_ast(quote(.i$pension > 0))
  expect_equal(node$kind, "call")
  expect_equal(node$fn, ">")
  expect_equal(node$args[[1L]], it_field("pension"))
  expect_equal(node$args[[2L]], it_lit(0))
  expect_equal(it_deparse(node), "(.i$pension > 0)")
})

test_that("pronoun-free sub-expressions are folded once", {
  expect_equal(it_build_ast(quote(1 + 2 * 3)), it_lit(7))
  expect_equal(it_deparse(it_build_ast(quote(.i$x + 2 * 3))), "(.i$x + 6)")
})

test_that("free variables fold in the supplied environment", {
  threshold <- 1000
  node <- it_build_ast(quote(.i$pension > threshold), environment())
  expect_equal(it_deparse(node), "(.i$pension > 1000)")
})

test_that("fields and time use are reported", {
  node <- it_build_ast(quote(.i$category == "ret" | .i$category == "dep"))
  expect_equal(it_fields(node), "category")
  expect_false(it_uses_t(node))

  weighted <- it_build_ast(quote(.i$pension * (.t - .b)))
  expect_setequal(it_fields(weighted), c("pension", "birth"))
  expect_true(it_uses_t(weighted))
})

test_that("indicators are detected structurally", {
  expect_true(it_is_indicator(it_build_ast(quote(.i$sex == "male"))))
  expect_true(it_is_indicator(it_build_ast(quote(.i$pension > 0))))
  # depends on time -> not a (time-invariant) indicator
  expect_false(it_is_indicator(it_build_ast(quote(.t > .b))))
  # not logical -> a general variable, not an indicator
  expect_false(it_is_indicator(it_build_ast(quote(.i$pension))))
})

test_that("reserved fields and bad pronouns are rejected", {
  expect_error(it_build_ast(quote(.i$E2R_died)), "reserved field")
  expect_error(it_build_ast(quote(.i)), "must be used as")
  expect_error(it_build_ast(quote(.z)), "Unknown pronoun")
})

test_that("non-scalar constants are rejected", {
  expect_error(it_build_ast(quote(c(1, 2, 3))), "single constant")
})

test_that("datey literals fold and keep their type", {
  node <- it_build_ast(quote(.i$birth > datey::datey(1960)))
  expect_equal(node$fn, ">")
  expect_equal(node$args[[2L]]$kind, "lit")
  expect_true(datey::is_datey(node$args[[2L]]$value))
})

test_that("the public entry handles bare expressions and formulas", {
  expect_equal(
    it_deparse(pronoun_expressions(.i$pension > 0)),
    "(.i$pension > 0)"
  )
  expect_equal(
    it_deparse(pronoun_expressions(~ .t - .b)),
    "(.t - .i$birth)"
  )
})

# ---- missing constants ------------------------------------------------------
#
# A missing value WRITTEN INTO an expression is refused at construction; a missing
# value arriving from the DATA is legal and flows through to a NaN result. The two
# are tested in different places for that reason -- the data side lives in the
# engine tests, and nothing here should be read as constraining it.

test_that("a constant that folds to NaN is refused", {
  # Provably NaN before any data is read, so there is no dataset under which the
  # expression means anything. Caught here rather than answered at run time.
  expect_error(it_build_ast(quote(log(-1))), "is NaN whatever the data holds")
  expect_error(it_build_ast(quote(0 / 0)), "is NaN whatever the data holds")
  expect_error(it_build_ast(quote(NaN)), "is NaN whatever the data holds")

  # Reported against the offending sub-expression, not the whole thing.
  expect_error(it_build_ast(quote(.i$pension + NaN)), "`NaN` is NaN")
})

test_that("an NA constant of any type is refused", {
  # `double`, `datey` and `durationy` carry a portable missing value and the rest
  # do not, so an NA logical, integer or string has no representation to be given
  # rather than merely an unhelpful one.
  expect_error(it_build_ast(quote(NA)), "is missing")
  expect_error(it_build_ast(quote(NA_real_)), "is missing")
  expect_error(it_build_ast(quote(NA_integer_)), "is missing")
  expect_error(it_build_ast(quote(NA_character_)), "is missing")
  expect_error(it_build_ast(quote(.i$pension > NA)), "is missing")
})

test_that("an infinite constant is accepted", {
  # NOT AN OVERSIGHT. `-Inf` is clean IEEE and is the neutral value in log space:
  # a mortality of `-Inf` is a rate of zero, contributing no exposure, which is
  # what a user guarding a missing column should reach for. A check written
  # against non-finiteness rather than missingness would take it away.
  expect_equal(it_build_ast(quote(-Inf))$value, -Inf)
  expect_equal(it_build_ast(quote(Inf))$value, Inf)
})
