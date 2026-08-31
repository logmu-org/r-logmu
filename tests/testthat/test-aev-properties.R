# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("aev properties", {

  aev_values <- matrix(c(
    # A, E, V
    0,  0,  0,
    0,  1,  1,
    1,  1,  1,
    1,  1,  3,
    1 + 0.00149,  1,  3,
    1 - 0.00149,  1,  3,
    1 + 0.00151,  1,  3,
    1 - 0.00151,  1,  3,
    40, 50, 60,
    10 + 10 + 100, 20 + 20 + 200, 20*20 + 20*20 + 200*200
  ), ncol = 3, byrow = TRUE)

  A <- aev_values[, 1]
  E <- aev_values[, 2]
  V <- aev_values[, 3]

  aev <- create_aev(A,E,V)

  expect_identical(aev$A, A)
  expect_identical(aev$E, E)
  expect_identical(aev$V, V)

  expect_equal(aev$A_minus_E, A - E)
  expect_equal(aev$A_minus_E_stddev, sqrt(V))
  expect_equal(aev$A_over_E, A / E)
  expect_equal(aev$log_A_over_E, log(A / E))
  expect_equal(aev$log_A_over_E_stddev, sqrt(V) / E)
  expect_equal(aev$log_A_over_E_95pc, qnorm(0.975) * sqrt(V) / E)
  expect_equal(aev$Pearson_residual, (A - E) / sqrt(V))

  deviance_residual <- sign(A-E) * sqrt(2*E/V*(A*log(A/E)-(A-E)))
  alpha <- A / E - 1
  deviance_residual_A_eq_E <- E / sqrt(V) * alpha *
    (1 + alpha * (-1/6 + alpha * (5/72 + alpha * (-83/2160))))

  deviance_residual <- ifelse(abs(alpha) > 0.0015, deviance_residual, deviance_residual_A_eq_E)
  deviance_residual <- ifelse(A == 0, -E * sqrt(2 / V), deviance_residual)

  expect_equal(aev$deviance_residual, deviance_residual)

  # Addition:
  # - None of `A`, `E` or `V` are negative
  #try(create_aev(A = -1, E = 2, V = 3))
  #> Error : A cannot be negative.
  # - `E` and `V` are either both zero or both non-zero
  #try(create_aev(A = 0, E = 0, V = 3))
  #> Error : V cannot be non-zero if E is zero.
  # - `A` cannot be non-zero if `E` and `V` are zero
  #try(create_aev(A = 1, E = 0, V = 0))
})

test_that("deviance residual is accurate across the series/general switch", {

  # AN INDEPENDENT REFERENCE, not a restatement of the implementation. The test
  # above mirrors `aev.cpp` and so can only catch transcription slips; this one
  # is derived from the deviance itself:
  #
  #   r_D * sqrt(V)/E = alpha * sqrt(f(alpha)),
  #   f(alpha) = sum_{n >= 0} 2*(-1)^n*alpha^n / ((n + 1)*(n + 2))
  #
  # `f` is dominated by its leading 1 and its terms fall like alpha^n, so
  # summing smallest-first suffers no cancellation and lands within a few ulps.
  reference <- function(alpha) {
    f <- 0
    for (n in 60:0) f <- f + 2 * (-1)^n * alpha^n / ((n + 1) * (n + 2))
    alpha * sqrt(f)
  }

  # `u` is built first and `alpha` read back from it, exactly as `aev.cpp` does,
  # so the reference and the implementation are talking about the same number
  # rather than about two roundings of it.
  magnitudes <- c(5e-5, 5e-4, 1e-3, 1.49e-3, 1.51e-3, 3e-3, 1e-2, 1e-1)
  u <- c(1 + magnitudes, 1 - magnitudes)
  alpha <- u - 1

  ones <- rep(1, length(u))
  actual <- create_aev(A = u, E = ones, V = ones)$deviance_residual
  expected <- vapply(alpha, reference, numeric(1))

  # THE TOLERANCE IS THE WHOLE POINT OF THIS TEST, so do not relax it to the
  # default. `testthat_tolerance()` is 1.49e-8, over a hundred times larger than
  # everything the cubic term of the series contributes -- a test at the default
  # cannot see whether that term is present at all, which was confirmed by
  # deleting it and watching the suite pass.
  #
  # Measured: this code errs by at most about 1.2e-13 across the sweep, and
  # dropping the cubic term takes that to 1.3e-10. 1e-11 sits between them with
  # room on both sides, and leaves slack for a different platform's log().
  expect_lt(max(abs((actual - expected) / expected)), 1e-11)
})

test_that("deviance residual treats an underflowed A/E as zero", {

  # A is positive but A/E falls below the smallest subnormal, so the quotient is
  # exactly zero. That is the same limit as A = 0 -- the bracket
  # A*log(A/E) - (A - E) tends to E in both cases -- and the general branch must
  # not be reached, because it would evaluate 0 * log(0) and return NaN.
  expect_equal(create_aev(A = 1e-200, E = 1e200, V = 1)$deviance_residual,
               -1e200 * sqrt(2))

  # The plain A = 0 case goes the same way, and is why the test is on the
  # quotient rather than on A.
  expect_equal(create_aev(A = 0, E = 100, V = 100)$deviance_residual,
               -100 * sqrt(2 / 100))
})

test_that("aev validity", {

  # (a) none of `A`, `E` or `V` are negative,
  expect_error(create_aev(-1, 2, 3))
  expect_error(create_aev(1, -2, 3))
  expect_error(create_aev(1, 2, -3))

  # (b) `E` and `V` are either both zero or both non-zero, and
  expect_error(create_aev(0, 0, 3))
  expect_error(create_aev(0, 2, 0))
  expect_error(create_aev(1, 0, 3))
  expect_error(create_aev(1, 2, 0))

  # (c) `A` cannot be non-zero if `E` and `V` are zero.
  expect_error(create_aev(1, 0, 0))
})

test_that("aev modification is illegal", {

  aev <- create_aev(1, 2, 3)

  expect_error(aev$A <- 1)
  expect_error(aev$B <- 1)
  expect_error(aev[1] <- 1)
  expect_error(aev[[1]] <- 1)
})

# An aev is a VECTOR OF RECORDS built on a list of three parallel fields, which
# is the shape POSIXlt has. Everything below pins the methods that make the two
# readings agree. They were wrong before and silently so: `[` picked fields and
# `length()` answered 3 whatever the aev held, neither of which bites until a
# breakdown returns a length-15 result and a user subsets it.

test_that("length counts records, not fields", {
  expect_equal(length(create_aev(1, 2, 3)), 1L)
  expect_equal(length(create_aev(c(1, 2, 3), c(1, 2, 3), c(1, 2, 3))), 3L)

  # The field count is still reachable, and `is_aev()` depends on asking for it
  # the long way round -- which is exactly what `length.aev` would have broken.
  expect_equal(length(unclass(create_aev(c(1, 2), c(1, 2), c(1, 2)))), 3L)
  expect_true(is_aev(create_aev(c(1, 2), c(1, 2), c(1, 2))))
})

test_that("subsetting takes records and keeps A, E and V together", {
  aev <- create_aev(A = c(10, 20, 30), E = c(1, 2, 3), V = c(4, 5, 6))

  one <- aev[2]
  expect_true(is_aev(one))
  expect_equal(length(one), 1L)
  expect_identical(one$A, 20)
  expect_identical(one$E, 2)
  expect_identical(one$V, 5)

  # A range, reversed order and repetition all behave as they do for a vector.
  expect_identical(aev[c(3, 1)]$A, c(30, 10))
  expect_identical(aev[c(1, 1)]$E, c(1, 1))
  expect_identical(aev[-1]$V, c(5, 6))

  # Calculated properties follow the records rather than being recomputed from
  # the whole, which is the point of subsetting before charting.
  expect_equal(aev[2]$A_over_E, 10)
})

test_that("labels are absent until set, and then subset with the records", {
  aev <- create_aev(A = c(10, 20, 30), E = c(1, 2, 3), V = c(4, 5, 6))

  # A hand-built aev has neither level of naming.
  expect_null(names(aev))
  expect_null(group_names(aev))

  names(aev) <- c("65-70", "70-75", "75-80")
  group_names(aev) <- "age"

  expect_identical(names(aev), c("65-70", "70-75", "75-80"))
  expect_identical(group_names(aev), rep("age", 3))

  # THE LABELS TRAVEL WITH THE RECORDS. A subset that kept all three labels, or
  # dropped them, would mislabel a chart rather than fail.
  expect_identical(names(aev[c(3, 1)]), c("75-80", "65-70"))
  expect_identical(group_names(aev[2]), "age")

  # And they can be taken off again.
  names(aev) <- NULL
  expect_null(names(aev))
})

test_that("labels must have one per record, and only group names recycle", {
  aev <- create_aev(A = c(10, 20, 30), E = c(1, 2, 3), V = c(4, 5, 6))

  # `names<-` does NOT recycle, deliberately: `names(x)[1] <- "y"` on an aev with
  # no names yet yields a length-1 vector, and recycling it would rename every
  # record instead of the first.
  expect_error(names(aev) <- "only one", "one label per record")
  expect_error(names(aev) <- c("a", "b", "c", "d"), "one label per record")
  expect_error(names(aev) <- c("a", "b", NA), "cannot be NA")

  # `group_names<-` DOES recycle a single value, matching `includes`: naming a
  # whole band set after one thing is the ordinary case.
  group_names(aev) <- "age"
  expect_identical(group_names(aev), rep("age", 3))
  expect_error(group_names(aev) <- c("age", "period"), "one label per record")
})

test_that("group_names is generic and refuses what it does not know", {
  # It serves both `includes` and `aev` now, so a type with no method must say so
  # rather than reaching for whichever implementation was written first.
  expect_error(group_names(1:3), "not implemented")
  expect_error(group_names("text"), "not implemented")
})

test_that("labels survive addition, and a mismatch refuses it", {
  # `Ops.aev` rebuilds from the unclassed members, so before this the sum of two
  # breakdowns came back with its rows unnamed -- silently, and only visible
  # once it reached a chart.
  male <- create_aev(A = c(10, 20), E = c(1, 2), V = c(4, 5))
  names(male) <- c("65-70", "70-75")
  group_names(male) <- "age"

  female <- create_aev(A = c(30, 40), E = c(3, 4), V = c(6, 7))
  names(female) <- c("65-70", "70-75")
  group_names(female) <- "age"

  both <- male + female
  expect_identical(names(both), c("65-70", "70-75"))
  expect_identical(group_names(both), rep("age", 2))
  expect_equal(both$A, c(40, 60))

  # Rows labelled differently are not the same rows, which is the misalignment
  # the length check already refuses one level up.
  other <- create_aev(A = c(30, 40), E = c(3, 4), V = c(6, 7))
  names(other) <- c("2000-2005", "2005-2010")
  expect_error(male + other, "`names` differ")

  same_names <- create_aev(A = c(30, 40), E = c(3, 4), V = c(6, 7))
  names(same_names) <- c("65-70", "70-75")
  group_names(same_names) <- "duration since entry"
  expect_error(male + same_names, "`group_names` differ")
})

test_that("an absent label set is not a disagreement", {
  # Adding a hand-built aev to a labelled one keeps the labels rather than
  # erasing them: the user asserted the rows correspond by adding at all.
  labelled <- create_aev(A = c(10, 20), E = c(1, 2), V = c(4, 5))
  names(labelled) <- c("65-70", "70-75")
  group_names(labelled) <- "age"

  bare <- create_aev(A = c(30, 40), E = c(3, 4), V = c(6, 7))

  expect_identical(names(labelled + bare), c("65-70", "70-75"))
  expect_identical(names(bare + labelled), c("65-70", "70-75"))
  expect_identical(group_names(bare + labelled), rep("age", 2))

  # Two bare aevs stay bare.
  expect_null(names(bare + bare))
  expect_null(group_names(bare + bare))
})

test_that("is_aev answers FALSE rather than throwing, whatever it is given", {
  # WRITTEN BECAUSE IT THREW. `class(x) == "aev"` is length 2 for anything with a
  # class vector of two, and R refuses a condition that is not length 1, so these
  # raised an error where a predicate has to answer. It went unnoticed while
  # `is_aev()` was internal and every caller already held an aev; exporting it
  # made the difference matter.
  expect_false(is_aev(matrix(1)))
  expect_false(is_aev(as.POSIXlt("2020-01-01", tz = "UTC")))
  expect_false(is_aev(data.frame(a = 1)))
  expect_false(is_aev(factor("a")))

  # Nothing with the right shape but the wrong class, and nothing with the right
  # class but the wrong shape.
  expect_false(is_aev(list(A = 1, E = 1, V = 1)))
  expect_false(is_aev(structure(list(A = 1, E = 1), class = "aev")))
  expect_false(is_aev(NULL))

  expect_true(is_aev(create_aev(1, 2, 3)))
})

test_that("a computed aev is not held to create_aev's invariants", {
  # THE RULE, AS A TEST. `create_aev()` still refuses a degenerate triple, because
  # somebody typed it and wants to be told. The same triple arrived at by
  # arithmetic is returned, because `exp()` underflowing to exactly zero is a true
  # statement about an impossible model and a job that has run for hours must not
  # be lost to it.
  expect_error(create_aev(A = 1, E = 0, V = 0), "A cannot be non-zero")

  degenerate <- create_aev_unchecked(A = 1, E = 0, V = 0)
  expect_true(is_aev(degenerate))
  expect_equal(degenerate$A, 1)

  # Adding is not a second chance to refuse it. Summing results is the last thing
  # a long run does, after the expensive part is already paid for.
  total <- degenerate + degenerate
  expect_equal(total$A, 2)
  expect_equal(total$E, 0)

  # The calculated properties take their IEEE values rather than being
  # special-cased -- pinned here so a well-meant guard cannot creep back in.
  expect_identical(degenerate$A_over_E, Inf)
  expect_identical(degenerate$log_A_over_E, Inf)
  expect_identical(degenerate$Pearson_residual, Inf)
  expect_identical(degenerate$A_minus_E, 1)
  expect_true(is.nan(degenerate$deviance_residual))
})
