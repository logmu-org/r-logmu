# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the constant-folding pass. Using the column scan's constraints, a field of a constant, NA-free
# column folds to a literal, and is.na folds when its answer is settled.

# A dataset's columns must all be the same length, as they would be in a data frame.
cols <- list(
  varying = c(10, 20, 30),   # not constant, no NAs
  fixed   = c(7, 7, 7),      # constant, no NAs
  gappy   = c(1, NA, 3),     # mixes NAs and values
  age     = c(30, 45, 60),   # ranges 30 to 60, no NAs -- for range-comparison folding
  birth   = datey::datey(c(1960, 1965, 1975)) # ranges 1960 to 1975 -- for comparison narrowing
)

root_kind <- function(res) res$kinds[[res$root + 1L]]
root_value <- function(res) res$values[[res$root + 1L]]

# Asserts that an expression folded all the way to a specific boolean, not merely that it folded to
# some literal -- a fold that reached the wrong answer would otherwise pass unnoticed.
expect_folds_to <- function(expr_res, expected) {
  expect_equal(root_kind(expr_res), "lit")
  expect_equal(expr_res$types[[expr_res$root + 1L]], "bool")
  expect_equal(root_value(expr_res), expected)
}

test_that("a constant column folds to a literal", {
  res <- cpp_veil_fold(it_ast(~ .i$fixed), cols)
  expect_equal(root_kind(res), "lit")
  expect_equal(res$types[[res$root + 1L]], "double")
})

test_that("a varying column stays a field", {
  res <- cpp_veil_fold(it_ast(~ .i$varying), cols)
  expect_equal(root_kind(res), "field")
})

test_that("is.na of a column with no NAs folds to false", {
  expect_folds_to(cpp_veil_fold(it_ast(~ is.na(.i$varying)), cols), "FALSE")
})

test_that("is.na of a column that mixes NAs and values does not fold", {
  res <- cpp_veil_fold(it_ast(~ is.na(.i$gappy)), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "is_na")
})

# Interval-based comparison folding. `age` ranges 30 to 60 across the data, so a comparison against a
# literal wholly above, below, or straddling that range settles or stays undecided accordingly. A
# bare field's interval is exactly its scanned range, which is why this covers what the earlier
# range-based fold did.

test_that("a comparison always false against the column's range folds to false", {
  expect_folds_to(cpp_veil_fold(it_ast(~ .i$age > 200), cols), "FALSE")
})

test_that("a comparison always true against the column's range folds to true", {
  expect_folds_to(cpp_veil_fold(it_ast(~ .i$age < 200), cols), "TRUE")
})

test_that("a comparison the range straddles is left undecided", {
  res <- cpp_veil_fold(it_ast(~ .i$age > 40), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "GT")
})

test_that("the literal-on-the-left ordering folds the same way", {
  # Every age exceeds 10, so this is true whichever side the literal is written on.
  expect_folds_to(cpp_veil_fold(it_ast(~ 10 < .i$age), cols), "TRUE")
})

test_that("a comparison against a column with NAs does not fold", {
  res <- cpp_veil_fold(it_ast(~ .i$gappy > 200), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "GT")
})

# Newly reachable now that folding reads intervals rather than a single column's scanned range.

test_that("a comparison over a derived expression folds", {
  # age is 30 to 60, so age - 100 is -70 to -40 and can never be positive.
  expect_folds_to(cpp_veil_fold(it_ast(~ .i$age - 100 > 0), cols), "FALSE")
})

test_that("a comparison of two columns folds when their ranges cannot overlap", {
  # varying is 10 to 30 and age is 30 to 60, so varying can never exceed age.
  expect_folds_to(cpp_veil_fold(it_ast(~ .i$varying > .i$age), cols), "FALSE")
})

test_that("a comparison of two columns whose ranges overlap is left undecided", {
  res <- cpp_veil_fold(it_ast(~ .i$varying < .i$age), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "LT")
})

test_that("an age expression in years folds against a datey column", {
  # `.x` is `.t - .b`, a durationy in clicks, while 65 is a plain number of YEARS. Narrowing puts
  # both sides in clicks first, so the fold compares like with like. `.t` is only bounded by the
  # representable calendar, so this one cannot settle -- it must stay a live comparison.
  res <- cpp_veil_fold(it_ast(~ .t - .b < 65), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "LE") # `< 65 years` narrowed to `<= threshold clicks`
  expect_true("34733399" %in% res$values)          # ceil(65 * 534360) - 1
})

# Comparison narrowing. A datey/durationy field compared to a plain number is, by the datey package's
# own rule, compared in years; narrowing rewrites the plain-number literal to an integer click
# threshold (rounded per operator direction) so it compares like-for-like with the field's own clicks
# -- which is what lets range folding, reading the field's scanned CLICK range, decide it correctly.
# `birth` ranges 1960 to 1975 across the data.

test_that("a datey comparison against a distant year folds FALSE once narrowed", {
  # Every birth year lies in 1960..1975, so none exceeds 1990. Without narrowing this compared the
  # column's CLICK range against the raw number 1990 and folded to TRUE -- the opposite answer.
  expect_folds_to(cpp_veil_fold(it_ast(~ .i$birth > 1990), cols), "FALSE")
})

test_that("a datey comparison the other way round folds TRUE", {
  expect_folds_to(cpp_veil_fold(it_ast(~ .i$birth < 1990), cols), "TRUE")
})

test_that("the literal-on-the-left ordering narrows and folds the same way", {
  expect_folds_to(cpp_veil_fold(it_ast(~ 1990 < .i$birth), cols), "FALSE")
})

test_that("narrowing rewrites the operator and the threshold to clicks", {
  # birth's minimum is exactly 1960, so `> 1960` becomes `>= 1960 in clicks, plus one click`. The
  # data reaches 1975, so the comparison stays live; the moniker and the threshold show the rewrite.
  res <- cpp_veil_fold(it_ast(~ .i$birth > 1960), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "GE")
  expect_true("1047345601" %in% res$values) # 1960 * 534360 + 1
})

test_that("equality against a year that lands between clicks is left alone", {
  # Narrowing declines a value that is not a whole click rather than settling it: no exact click
  # threshold exists, so the comparison survives to be evaluated at run time.
  res <- cpp_veil_fold(it_ast(~ .i$birth == (1960 + 1 / 13)), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "EQ")
})

test_that("an ordering against a year that lands between clicks is also left alone", {
  # The whole-click rule applies to every comparison, not just equality. Without it the year-space
  # and click-space answers could differ over the single click on the boundary.
  res <- cpp_veil_fold(it_ast(~ .i$birth > (1960 + 1 / 13)), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "GT") # unchanged: not rewritten to GE
})

test_that("equality against a whole click narrows to an exact click comparison", {
  res <- cpp_veil_fold(it_ast(~ .i$birth == 1965), cols)
  expect_equal(root_kind(res), "call")
  expect_equal(res$labels[[res$root + 1L]], "EQ")
  expect_true("1050017400" %in% res$values) # 1965 * 534360
})
