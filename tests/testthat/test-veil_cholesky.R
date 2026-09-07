# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# The small dense linear algebra a Newton step needs: a Cholesky, a triangular solve, an inverse, a
# quadratic form and a trace. `cpp_veil_cholesky()` is a window onto veil/Cholesky.hpp and nothing
# in the package calls it.
#
# THE PACKING IS ROW-MAJOR UPPER, the convention FitRecipe.hpp defines: (1,1), (1,2) ... (1,k),
# (2,2) ... The helpers below spell that out in a loop rather than borrowing `upper.tri()`, which
# walks the triangle in column order and would agree with the C++ only by an accident of symmetry.

pack_upper <- function(m)
{
  k <- nrow(m)
  out <- numeric(0)
  for (row in seq_len(k))
  {
    for (col in row:k) out <- c(out, m[row, col])
  }
  out
}

unpack_symmetric <- function(v, k)
{
  m <- matrix(0, k, k)
  i <- 1L
  for (row in seq_len(k))
  {
    for (col in row:k)
    {
      m[row, col] <- v[i]
      m[col, row] <- v[i]
      i <- i + 1L
    }
  }
  m
}

unpack_triangle <- function(v, k)
{
  m <- matrix(0, k, k)
  i <- 1L
  for (row in seq_len(k))
  {
    for (col in row:k)
    {
      m[row, col] <- v[i]
      i <- i + 1L
    }
  }
  m
}

# Everything comes from one factorisation, so one call answers everything.
factorise <- function(m, b = rep(0, nrow(m)), other = diag(nrow(m)))
{
  cpp_veil_cholesky(pack_upper(m), nrow(m), as.double(b), pack_upper(other))
}

test_that("the packing read back is row-major upper, not column-major", {
  # A = U^T U with U = [[2, 1, -1], [0, 3, 1], [0, 0, sqrt(3)]], worked out by hand:
  #   U11 = sqrt(4) = 2;  U12 = 2 / 2 = 1;  U13 = -2 / 2 = -1
  #   U22 = sqrt(10 - 1) = 3;  U23 = (2 - 1 * -1) / 3 = 1
  #   U33 = sqrt(5 - 1 - 1) = sqrt(3)
  a <- matrix(c(4, 2, -2, 2, 10, 2, -2, 2, 5), 3, 3)
  res <- factorise(a)

  expect_true(res$positive_definite)
  expect_identical(res$terms, 3L)

  # Row-major would be read as c(2, 1, -1, 3, 1, sqrt(3)); a column-major reading of the same
  # triangle would give c(2, 1, 3, -1, 1, sqrt(3)), so this pins the convention rather than just
  # the arithmetic.
  expect_equal(res$factor, c(2, 1, -1, 3, 1, sqrt(3)))
})

test_that("a hand-worked 2x2 factors, solves, inverts, and gives the quadratic form and trace", {
  # A = [[4, 2], [2, 5]]. U11 = 2, U12 = 1, U22 = sqrt(5 - 1) = 2.
  # det A = 16, so A^-1 = [[5, -2], [-2, 4]] / 16.
  a <- matrix(c(4, 2, 2, 5), 2, 2)
  b <- c(6, 3)
  other <- matrix(c(1, 2, 2, 3), 2, 2)
  res <- factorise(a, b, other)

  expect_true(res$positive_definite)
  expect_equal(res$factor, c(2, 1, 2))
  expect_equal(res$inverse, c(5 / 16, -2 / 16, 4 / 16))

  # A^-1 b = (5 * 6 - 2 * 3, -2 * 6 + 4 * 3) / 16 = (24, 0) / 16 = (1.5, 0).
  expect_equal(res$solve, c(1.5, 0))

  # b^T A^-1 b = 6 * 1.5 + 3 * 0 = 9.
  expect_equal(res$quadratic_form, 9)

  # tr(B A^-1), both symmetric, so the Frobenius inner product:
  #   1 * 5/16 + 3 * 4/16 + 2 * (2 * -2/16) = (5 + 12 - 8) / 16 = 9/16.
  expect_equal(res$trace, 9 / 16)
})

test_that("a hand-worked 3x3 solves, inverts and gives the quadratic form and trace", {
  a <- matrix(c(4, 2, -2, 2, 10, 2, -2, 2, 5), 3, 3)
  b <- c(2, 4, 6)
  other <- matrix(c(1, 0.5, 0.25, 0.5, 2, -1, 0.25, -1, 3), 3, 3)
  res <- factorise(a, b, other)

  # Forward U^T y = b: y = (1, 1, 6 / sqrt(3)). Back U x = y: x3 = 2, x2 = (1 - 2) / 3 = -1/3,
  # x1 = (1 + 1/3 + 2) / 2 = 5/3.
  expect_equal(res$solve, c(5 / 3, -1 / 3, 2))

  # det A = 108, and the cofactors give A^-1 * 108 = [[46, -14, 24], [-14, 16, -12], [24, -12, 36]].
  expect_equal(res$inverse, c(46, -14, 24, 16, -12, 36) / 108)

  # b^T A^-1 b = ||U^-T b||^2 = 1 + 1 + 12 = 14, which is also b . x = 10/3 - 4/3 + 12.
  expect_equal(res$quadratic_form, 14)
  expect_equal(res$quadratic_form, sum(b * res$solve))

  # The trace against an explicit product, summed independently of the C++.
  inverse <- unpack_symmetric(res$inverse, 3)
  expect_equal(res$trace, sum(diag(other %*% inverse)))
})

test_that("the identity factors to itself and leaves everything alone", {
  identity <- diag(4)
  b <- c(1, -2, 3, -4)
  other <- matrix(c(2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 5, 0, 0, 0, 0, 7), 4, 4)
  res <- factorise(identity, b, other)

  expect_true(res$positive_definite)
  expect_equal(res$factor, pack_upper(identity))
  expect_equal(res$inverse, pack_upper(identity))
  expect_equal(res$solve, b)
  expect_equal(res$quadratic_form, sum(b^2))
  expect_equal(res$trace, 2 + 3 + 5 + 7)
})

test_that("k = 0 means there is nothing to solve rather than a crash", {
  res <- cpp_veil_cholesky(numeric(0), 0L, numeric(0), numeric(0))

  expect_true(res$positive_definite)
  expect_identical(res$terms, 0L)
  expect_identical(res$failed_parameter, 0L) # `terms`, the past-the-end sentinel
  expect_identical(length(res$factor), 0L)
  expect_identical(length(res$solve), 0L)
  expect_identical(length(res$inverse), 0L)
  expect_identical(res$quadratic_form, 0)
  expect_identical(res$trace, 0)
})

test_that("k = 1 is the scalar case", {
  res <- cpp_veil_cholesky(4, 1L, 6, 8)

  expect_true(res$positive_definite)
  expect_equal(res$factor, 2)
  expect_equal(res$inverse, 0.25)
  expect_equal(res$solve, 1.5)
  expect_equal(res$quadratic_form, 9)
  expect_equal(res$trace, 2)
})

test_that("failed_parameter is the past-the-end sentinel when nothing failed", {
  res <- factorise(diag(3))
  expect_true(res$positive_definite)
  expect_identical(res$failed_parameter, 3L)
})

test_that("the factor multiplied back reproduces the input", {
  set.seed(20260904)
  for (k in 1:6)
  {
    x <- matrix(stats::rnorm(4 * k * k), 4 * k, k)
    a <- crossprod(x) + diag(k) # crossprod alone can be near-singular; the ridge keeps it honest
    res <- factorise(a)

    expect_true(res$positive_definite)
    upper <- unpack_triangle(res$factor, k)
    expect_equal(t(upper) %*% upper, a, tolerance = 1e-12)

    # The diagonal of U is positive by construction, which is what makes the factor unique.
    expect_true(all(diag(upper) > 0))
  }
})

test_that("A times its inverse is the identity", {
  set.seed(11223344)
  for (k in 1:6)
  {
    x <- matrix(stats::rnorm(4 * k * k), 4 * k, k)
    a <- crossprod(x) + diag(k)
    res <- factorise(a)

    inverse <- unpack_symmetric(res$inverse, k)
    expect_equal(a %*% inverse, diag(k), tolerance = 1e-9)
  }
})

test_that("the quadratic form agrees with the solve, and the trace with an explicit product", {
  set.seed(55667788)
  for (k in 1:6)
  {
    x <- matrix(stats::rnorm(4 * k * k), 4 * k, k)
    a <- crossprod(x) + diag(k)
    b <- stats::rnorm(k)
    y <- matrix(stats::rnorm(4 * k * k), 4 * k, k)
    other <- crossprod(y)
    res <- factorise(a, b, other)

    expect_equal(res$quadratic_form, sum(b * res$solve), tolerance = 1e-9)
    expect_true(res$quadratic_form > 0)

    inverse <- unpack_symmetric(res$inverse, k)
    expect_equal(res$trace, sum(diag(other %*% inverse)), tolerance = 1e-9)
  }
})

test_that("the triangular inverse is the inverse of the factor", {
  a <- matrix(c(4, 2, -2, 2, 10, 2, -2, 2, 5), 3, 3)
  res <- factorise(a)

  upper <- unpack_triangle(res$factor, 3)
  triangular <- unpack_triangle(res$triangular_inverse, 3)
  expect_equal(upper %*% triangular, diag(3), tolerance = 1e-12)
})

test_that("a zero row and column names the parameter it belongs to", {
  # Parameter 2 (zero-based 1) has no exposure behind it at all.
  a <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 1), 3, 3)
  res <- factorise(a)

  expect_false(res$positive_definite)
  expect_identical(res$failed_parameter, 1L)
})

test_that("two identical covariates are caught at the second of them", {
  # The realistic collinearity: the third covariate is a copy of the first. Chosen so that
  # A = t(X) %*% X has power-of-two entries and the failing pivot is exactly zero rather than a
  # rounding-sized negative, which keeps the test from depending on which way the last bit fell.
  x <- cbind(c(1, 1, 1, 1), c(1, -1, 1, -1), c(1, 1, 1, 1))
  a <- crossprod(x)
  expect_identical(a[1, 1], 4) # exactness of the construction, not of the C++

  res <- factorise(a)
  expect_false(res$positive_definite)
  expect_identical(res$failed_parameter, 2L)
})

test_that("a negative pivot is a failure, not a factorisation", {
  # Indefinite rather than singular: the leading minor is fine and the second pivot goes negative.
  a <- matrix(c(1, 2, 2, 1), 2, 2)
  res <- factorise(a)

  expect_false(res$positive_definite)
  expect_identical(res$failed_parameter, 1L)
})

test_that("a NaN on the diagonal does not pass as positive definite", {
  a <- matrix(c(4, 2, 2, NaN), 2, 2)
  res <- factorise(a)

  expect_false(res$positive_definite)
  expect_identical(res$failed_parameter, 1L)

  # Not a plausible-looking answer either.
  expect_true(all(is.nan(res$solve)))
  expect_true(all(is.nan(res$inverse)))
  expect_true(is.nan(res$quadratic_form))
  expect_true(is.nan(res$trace))
})

test_that("a NaN off the diagonal is caught at the column it feeds", {
  # The NaN sits at (1, 3) and enters only the third pivot, so nothing hides in a corner.
  a <- matrix(c(4, 2, NaN, 2, 10, 2, NaN, 2, 5), 3, 3)
  res <- factorise(a)

  expect_false(res$positive_definite)
  expect_identical(res$failed_parameter, 2L)
  expect_true(is.nan(res$quadratic_form))
})

test_that("an NA_real_ is a NaN as far as the factorisation is concerned", {
  a <- matrix(c(4, 2, 2, NA_real_), 2, 2)
  res <- factorise(a)

  expect_false(res$positive_definite)
  expect_identical(res$failed_parameter, 1L)
})

test_that("a failed factorisation answers NaN everywhere rather than throwing", {
  a <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 1), 3, 3)
  res <- factorise(a, b = c(1, 2, 3), other = diag(3))

  expect_false(res$positive_definite)
  expect_identical(length(res$solve), 3L)
  expect_true(all(is.nan(res$solve)))
  expect_true(all(is.nan(res$inverse)))
  expect_true(all(is.nan(res$triangular_inverse)))
  expect_true(is.nan(res$quadratic_form))
  expect_true(is.nan(res$trace))
})

test_that("the same input twice gives a bit-identical answer", {
  set.seed(998877)
  x <- matrix(stats::rnorm(60), 20, 3)
  a <- crossprod(x) + diag(3)
  b <- stats::rnorm(3)
  y <- matrix(stats::rnorm(60), 20, 3)
  other <- crossprod(y)

  first <- factorise(a, b, other)
  second <- factorise(a, b, other)
  expect_identical(first, second)
})

test_that("a shape mismatch is a caller error", {
  expect_error(cpp_veil_cholesky(c(1, 0, 1), 3L, c(0, 0, 0), c(1, 0, 0, 1, 0, 1)), "packed_upper")
  expect_error(cpp_veil_cholesky(c(1, 0, 0, 1, 0, 1), 3L, c(0, 0), c(1, 0, 0, 1, 0, 1)),
               "right_hand_side")
  expect_error(cpp_veil_cholesky(c(1, 0, 0, 1, 0, 1), 3L, c(0, 0, 0), c(1, 0)),
               "packed_upper_other")
  expect_error(cpp_veil_cholesky(numeric(0), -1L, numeric(0), numeric(0)), "terms")
})

# THE NEAR-SINGULARITY TEST AND THE FAILURE CODES.
#
# Exact positivity is too weak on its own. Two covariates agreeing to twelve digits give a small
# POSITIVE pivot, so the factorisation succeeds and the inverse comes back enormous -- a fit that
# looks like an answer and is not one. The pivot is therefore compared with the diagonal entry it
# came from, and that ratio is `1 - R^2`: the fraction of a covariate that is its own rather than a
# rehash of its predecessors.

test_that("a pivot that survives positivity but not the ratio test is refused", {

  # Two columns agreeing to five digits, as a Gram matrix. The third is independent.
  #
  # FIVE DIGITS RATHER THAN TWELVE, AND THAT IS NOT TIMIDITY. For columns differing by e the pivot
  # is of order e^2 against a diagonal of order 1, so e = 1e-12 puts the true pivot near 1e-24 --
  # four orders BELOW double precision relative to the diagonal. What the factorisation then computes
  # is rounding noise whose sign is arbitrary, and the matrix is caught by positivity rather than by
  # the ratio. At e = 1e-5 the pivot ratio is about 2e-11: far inside the 1e-7 band, and still five
  # digits clear of the noise floor, so it is genuinely a small POSITIVE pivot and the ratio test is
  # the only thing that can catch it.
  epsilon <- 1e-5
  x1 <- c(1, 1, 1)
  x2 <- x1 + c(epsilon, 0, 0)
  x3 <- c(1, -1, 0)
  design <- cbind(x1, x2, x3)
  gram <- t(design) %*% design

  result <- factorise(gram)

  expect_false(result$positive_definite)
  expect_identical(result$status, "nearly_singular")

  # Detected at column 2 -- zero-based 1 -- because columns 1 and 2 are the pair.
  expect_identical(result$failed_parameter, 1L)

  # And the diagnosis: column 2 is one times column 1, to within the epsilon that separates them.
  expect_equal(result$dependency, 1, tolerance = 1e-4)
})

test_that("the ratio test does not refuse a legitimately correlated model", {

  # Age and age squared over ages 60 to 90, which is about as collinear as a sane actuarial model
  # gets. Uncentred `1 - R^2` here is around 0.015, five orders clear of the 1e-7 threshold, so a
  # threshold that rejected this would be useless.
  ages <- seq(60, 90, by = 0.5)
  design <- cbind(ages, ages^2)
  result <- factorise(t(design) %*% design)

  expect_true(result$positive_definite)
  expect_identical(result$status, "ok")
})

test_that("the failure codes tell a NaN apart from a singularity", {

  # A NaN must be asked about FIRST. It fails every comparison, so a test written the other way
  # round would report it as an ordinary singularity and send the user hunting for a collinear
  # covariate that does not exist.
  nan_matrix <- diag(3)
  nan_matrix[2, 2] <- NaN
  nan_result <- factorise(nan_matrix)
  expect_identical(nan_result$status, "not_finite")
  expect_identical(nan_result$failed_parameter, 1L)

  # A NaN says nothing about collinearity, so there is no dependency to report.
  expect_identical(length(nan_result$dependency), 0L)

  # A negative pivot is not a near-singularity: the matrix is not a covariance at all.
  negative <- diag(c(1, -1, 1))
  negative_result <- factorise(negative)
  expect_identical(negative_result$status, "not_positive_definite")
  expect_identical(negative_result$failed_parameter, 1L)

  # An exactly singular matrix is caught by positivity rather than by the ratio.
  singular <- diag(c(1, 0, 1))
  expect_identical(factorise(singular)$status, "not_positive_definite")

  # And a healthy matrix says so.
  expect_identical(factorise(diag(3))$status, "ok")
})

test_that("the dependency names the earlier covariates rather than the pivot", {

  # THE PIVOT INDEX IS NOT THE CULPRIT. Columns 1 and 3 are identical, so the factorisation fails at
  # column 3 -- because 1 and 2 are still independent -- and reporting that index alone would send a
  # user to the wrong column. The back-substitution says which earlier columns account for it.
  x1 <- c(1, 0, 0, 1)
  x2 <- c(0, 1, 0, 1)
  x3 <- x1
  design <- cbind(x1, x2, x3)
  result <- factorise(t(design) %*% design)

  expect_false(result$positive_definite)
  expect_identical(result$failed_parameter, 2L)

  # Column 3 is 1.0 times column 1 plus 0.0 times column 2.
  expect_equal(result$dependency, c(1, 0), tolerance = 1e-8)

  # A failure at the very first parameter has no earlier column to blame, so the answer is empty --
  # that covariate is degenerate on its own rather than a duplicate of anything.
  first <- factorise(diag(c(0, 1, 1)))
  expect_identical(first$failed_parameter, 0L)
  expect_identical(length(first$dependency), 0L)
})
