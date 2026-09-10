# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The covariates algebra: the right-hand side of a fitted model, written as `covariates(...)` and
# combined with `*` and `c()`.
#
# Single quotes throughout for anything containing a double quote, so no backslash escape appears in
# this file. That is a shell trap rather than an R one, but the file is edited from a shell often
# enough to be worth avoiding.

deparsed <- function(x) vapply(unclass(x), function(item) it_deparse(item$ast), character(1L))

test_that("covariates captures pronoun expressions without evaluating them", {
  # `.i` IS NOT AN OBJECT and exists only inside non-standard evaluation, which is precisely why
  # this is a constructor rather than `list(...)`. Evaluating the argument would fail outright.
  expect_error(eval(quote(.i$sex)), "not found")

  cv <- covariates(.i$sex == "M", .x)
  expect_true(is_covariates(cv))
  expect_identical(length(cv), 2L)
  expect_identical(deparsed(cv), c('(.i$sex == "M")', "(.t - .i$birth)"))

  # The narrowest provable type is kept, exactly as `variable()` does: a time-invariant logical is
  # an indicator, and anything using `.t` is not static.
  expect_true(is_indicator(unclass(cv)[[1L]]))
  expect_false(is_static_variable(unclass(cv)[[2L]]))
})

test_that("names are optional and ride on the terms", {
  cv <- covariates(male = .i$sex == "M", .i$sex == "F")
  expect_identical(names(cv), c("male", ""))
  expect_null(names(covariates(.i$sex == "M")))
})

test_that("a formula and a prebuilt variable are both accepted", {
  shape <- variable(.x)

  # `it_capture()` already parses a formula and splices a symbol bound to a variable, so neither
  # needs a path of its own here.
  expect_identical(deparsed(covariates(~ .i$pension)), ".i$pension")
  expect_identical(deparsed(covariates(shape)), "(.t - .i$birth)")
})

test_that("a covariates object splices into another", {
  cv <- covariates(male = .i$sex == "M")
  both <- covariates(cv, .i$pension)

  expect_identical(length(both), 2L)
  expect_identical(names(both), c("male", ""))

  # A name cannot label a spliced GROUP, because a covariate name identifies one coefficient.
  expect_error(covariates(everything = cv), "names a group")
})

test_that("multiplying by a scalar function distributes over the list", {
  # THE COMMON SHAPE IS `X = I * phi`: a set of indicators times one age shape. That is the whole
  # reason the algebra exists rather than making the user write out every product.
  shape <- variable(.x)
  cv <- covariates(.i$sex == "M", .i$sex == "F")

  expect_identical(deparsed(cv * shape),
                   c('((.i$sex == "M") * (.t - .i$birth))',
                     '((.i$sex == "F") * (.t - .i$birth))'))

  # Either way round: `*` dispatches on whichever operand carries the method, and the list decides
  # the shape of the answer.
  expect_identical(deparsed(shape * cv), deparsed(cv * shape))

  # A plain number is a rescaling, which is meaningful.
  expect_identical(deparsed(cv * 2),
                   c('((.i$sex == "M") * 2)', '((.i$sex == "F") * 2)'))
})

test_that("multiplying two lists gives the Cartesian product", {
  region <- covariates(north = .i$region == "N", south = .i$region == "S")
  sex <- covariates(male = .i$sex == "M", female = .i$sex == "F")
  crossed <- region * sex

  expect_identical(length(crossed), 4L)
  expect_identical(names(crossed), c("north.male", "north.female", "south.male", "south.female"))

  # Row-major: the left side moves slowest.
  expect_identical(deparsed(crossed)[[2L]], '((.i$region == "N") * (.i$sex == "F"))')
})

test_that("a product is named only when both sides are", {
  # Joining a name to an empty one would label every product of that row identically, which is
  # worse than leaving them unnamed -- a name is supposed to identify one coefficient.
  named <- covariates(male = .i$sex == "M")
  anonymous <- covariates(.i$smoker, !.i$smoker)

  expect_null(names(named * anonymous))
  expect_identical(names(named * covariates(smoker = .i$smoker)), "male.smoker")
})

test_that("disjoint asserts exclusivity and refuses what R can see is wrong", {
  dj <- disjoint(.i$sex == "M", .i$sex == "F")

  expect_true(is_disjoint(dj))
  expect_false(is_disjoint(covariates(.i$sex == "M", .i$sex == "F")))

  # R CANNOT PROVE INDICATOR-NESS -- a column's type is unknown until the data arrives, so
  # `.i$is_male` reaches here as a static variable and the engine confirms it. What R CAN prove is
  # that a term uses `.t`, and a membership that moved during an individual's exposure would not
  # hoist out of the integral, which is where the whole saving comes from.
  expect_no_error(disjoint(.i$is_male))
  expect_error(disjoint(.x > 5), "must be an indicator that does not vary with time")
})

test_that("disjointness survives a shape and a crossing, but not concatenation", {
  dj <- disjoint(.i$sex == "M", .i$sex == "F")
  shape <- variable(.x)

  # A COMMON FACTOR CANNOT MAKE EXCLUSIVE TERMS OVERLAP: if I_j I_l = 0 then I_j phi I_l phi = 0.
  expect_true(is_disjoint(dj * shape))

  # And a crossing is exclusive if EITHER side is -- nobody in two regions is in two region-sex
  # cells, whatever the sexes do.
  expect_true(is_disjoint(dj * covariates(.i$smoker, !.i$smoker)))
  expect_true(is_disjoint(covariates(.i$smoker, !.i$smoker) * dj))

  # BUT IT DOES NOT COMPOSE. Two internally exclusive sets say nothing about each other, so a
  # concatenation asserts nothing -- exactly as `includes()` checks nothing across the includes it
  # composes. Claiming otherwise would assert something nobody stated and the engine would then go
  # on to omit integrals that are not zero.
  expect_false(is_disjoint(c(dj, dj)))
  expect_identical(length(c(dj, dj)), 4L)

  # Subsetting keeps it, since a subset of mutually exclusive terms is still mutually exclusive.
  expect_true(is_disjoint(dj[1L]))
})

test_that("the operator refuses what has no object-level meaning", {
  shape <- variable(.x)
  cv <- covariates(.i$sex == "M")

  # Arithmetic between two plain variables stays inside pronoun expressions. That is the existing
  # rule and none of this relaxes it.
  expect_error(shape * shape, "not defined for logmu function objects")
  expect_error(cv + shape, "not defined for logmu function objects")
  expect_error(cv * "text", "must be a `covariates`, a `variable` or a single number")
})

test_that("covariates objects are immutable", {
  cv <- covariates(.i$sex == "M")
  expect_error(cv[[1L]] <- 1, "immutable")
})
