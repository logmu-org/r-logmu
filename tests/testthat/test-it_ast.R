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

test_that("`if` selects per individual, over a whole population at once", {
  # THIS IS WHAT `simulate()` LEANS ON: `.i` holding many individuals, and a
  # sex-switched mortality evaluated over all of them in one call. It used to
  # take `isTRUE()` of the condition, which picked the else branch for every
  # vector condition and was only ever right for one individual at a time.
  ast <- it_build_ast(quote(if (.i$male) 1 else 2))
  expect_identical(it_eval(ast, list(male = c(TRUE, FALSE, TRUE))), c(1, 2, 1))
  expect_identical(it_eval(ast, list(male = TRUE)), 1)
  expect_identical(it_eval(ast, list(male = FALSE)), 2)

  # An NA condition is an NA answer, not a silent else.
  expect_identical(it_eval(ast, list(male = c(TRUE, NA))), c(1, NA))

  # A branch that varies over time recycles against a per-individual condition.
  timed <- it_build_ast(quote(if (.i$male) .t - .b else 0))
  ages <- it_eval(timed, list(male = TRUE, birth = datey::datey(1950)),
                  datey::datey(c(2015, 2016)))
  expect_equal(as.double(ages), c(65, 66))

  # And a `datey`-valued branch keeps its class, which `ifelse()` would strip.
  dated <- it_build_ast(quote(if (.i$male) .b else .t))
  picked <- it_eval(dated, list(male = c(TRUE, FALSE), birth = datey::datey(c(1950, 1960))),
                    datey::datey(2015))
  expect_true(datey::is_datey(picked))
  expect_identical(picked, datey::datey(c(1950, 2015)))
})

test_that("reserved fields and bad pronouns are rejected", {
  expect_error(it_build_ast(quote(.i$E2R_died)), "reserved field")

  # THE PREFIX AND THE COUNT ARE WHAT IS RESERVED, case-insensitively, as
  # `exp_data()` recognises them: a column only `simulate()` knows about, one
  # that does not exist yet, and the record count in any casing.
  expect_error(it_build_ast(quote(.i$E2R_cap)), "names beginning `E2R_` are reserved (case insensitive)", fixed = TRUE)
  expect_error(it_build_ast(quote(.i$E2R_anything)), "reserved field")
  expect_error(it_build_ast(quote(.i$e2r_start)), "reserved field")
  expect_error(it_build_ast(quote(.i$count)), "name `count`", fixed = TRUE)
  expect_error(it_build_ast(quote(.i$COUNT)), "reserved field")

  # But only those: `birth` and `uid` are special to `exp_data()` and may be
  # read, and a near miss of the prefix is an ordinary field.
  expect_identical(it_build_ast(quote(.i$birth))$name, "birth")
  expect_identical(it_build_ast(quote(.i$uid))$name, "uid")
  expect_identical(it_build_ast(quote(.i$E2Rcap))$name, "E2Rcap")
  expect_identical(it_build_ast(quote(.i$counter))$name, "counter")
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
