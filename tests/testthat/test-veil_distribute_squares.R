# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `passDistributeSquares` and the rotation step in `passHoistFromIntegrate` are ONE OPTIMISATION IN
# TWO PLACES, and neither is worth anything alone -- measured: with the pass switched off, the
# rotation changes not one number.
#
#     (I phi) * (I phi)  --distribute-->  (I*I) * (phi*phi)
#                        --fold------->   I * (phi*phi)
#                        --rotate+hoist-> I * integrate(phi^2 mu)
#
# and that last integral is the SAME NODE for every term, so sharing turns k diagonal integrals into
# one. The shape it is for -- an indicator times an age curve -- is the commonest covariate there is.
#
# NOTHING NUMERIC CAN WITNESS IT. The rewrite is value-preserving, so the tests below watch
# `squares_distributed` for whether the pass fired and the instruction and slot counts for whether
# the rotation and the hoist then banked it.

clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

start_clicks <- c(2010, 2010, 2010) * clicks_per_year

squares_cols <- list(
  birth     = datey::datey(c(1940, 1945, 1950)),
  amount    = c(1000, 2500, 400),
  male      = c(TRUE, FALSE, TRUE),
  g1        = c(1, 0, 0),
  g2        = c(0, 1, 0),
  g3        = c(0, 0, 1),

  # A TIME-INVARIANT duration, which is where the coercion inside the square actually bites: a
  # duration on the time vector is already years, so `ToDouble` lowers to nothing, but a scalar one
  # is CLICKS until something reads it as years.
  service   = datey::durationy(c(2, 5, 1)),
  E2R_start = datey_clicks(start_clicks),
  E2R_end   = datey_clicks(start_clicks + c(12L, 8L, 4L) * quarter),
  E2R_died  = c(TRUE, FALSE, TRUE)
)

# A mortality that VARIES OVER TIME, which the hoist needs: against a constant one the whole
# integrand is invariant and the pass declines it, having nothing left to integrate.
squares_table <- mortality_table(
  x0 = 50, t0 = 2000,
  log_mu = matrix(seq(-6, -2, length.out = 60 * 30), nrow = 60, ncol = 30))

squares_aev <- function(weight, include = NULL) {
  cpp_veil_aev(it_obj(squares_table), weight, squares_cols, quarter, include, 1, 1L)
}

shaped_terms <- function(k) {
  lapply(seq_len(k), function(j) it_capture(str2lang(sprintf(".i$g%d * .x", j)), globalenv()))
}

squares_fit <- function(k) {
  cpp_veil_fit(it_obj(squares_table), shaped_terms(k), list(rep(0, k)), it_ast(~ .i$amount),
               NULL, NULL, squares_cols, quarter, NULL, FALSE, 1L)
}

test_that("an indicator times a shape has its square distributed", {
  # An AEV's V is `mu * (w * w)`, so a weight of `I * phi` puts `(I phi)^2` in front of the pass in
  # the simplest setting there is.
  res <- squares_aev(it_ast(~ .i$male * .x))
  expect_equal(res$squares_distributed, 1L)
})

test_that("the distributed square gives the identical answer, bit for bit", {
  # THE ORACLE IS A DIFFERENT CODE PATH ENTIRELY. `V` under a weight of `I * phi` is
  # `Omega integrate(mu I^2 phi^2)`, and since `I` is an indicator that is the same number as `V`
  # under a weight of `phi` with the individuals `I` selects taken by an `include` instead. One
  # reaches it by folding a square, the other by clipping exposure to nothing.
  #
  # AND IT IS EXACT. Distributing regroups the multiplications, which is normally enough to move the
  # last place -- but `I` is zero or one, so `I*I` folds to `I` and the surviving product is the
  # same one, in the same order.
  weighted <- squares_aev(it_ast(~ .i$male * .x))
  included <- squares_aev(it_ast(~ .x * 1), include = include(.i$male))

  expect_identical(weighted$contributions$V, included$contributions$V)
  expect_identical(weighted$contributions$E, included$contributions$E)

  # The woman is weighted to nothing, and that is a zero rather than a small number.
  expect_identical(weighted$contributions$V[[2]], 0)
})

test_that("it declines a bare indicator, which the fold already owns", {
  # `I * I` is a square of something ALREADY zero or one, so `passFoldIndicatorSquares` collapses it
  # outright. Collapsing beats distributing, and doing both would only add nodes.
  res <- squares_aev(it_ast(~ .i$male))
  expect_equal(res$squares_distributed, 0L)
  expect_equal(sum(res$monikers == "integrate"), 1L)
})

test_that("it declines a product with no indicator in it", {
  # `(amount * phi)^2` distributes to `(amount^2)(phi^2)`, and nothing folds: no square collapses
  # and no factor becomes time-invariant that was not already. That is a rounding change for
  # nothing, so it is refused.
  res <- squares_aev(it_ast(~ .i$amount * .x))
  expect_equal(res$squares_distributed, 0L)
})

test_that("it declines a literal factor, which is a constant rather than an indicator", {
  # 0 and 1 are the most zero-or-one values there are, so a naive test fires on `(.x * 1)^2` and
  # splits off a constant that folding and the hoist handle far more cheaply. The point of the pass
  # is an indicator that VARIES BETWEEN INDIVIDUALS.
  res <- squares_aev(it_ast(~ .x * 1))
  expect_equal(res$squares_distributed, 0L)
})

test_that("the rotation banks it: the indicator leaves the integral", {
  # WITHOUT THE ROTATION THE HOIST STOPS ONE LEVEL SHORT. The integrand is `((I*I)(phi*phi)) * mu`,
  # where both sides of the outer product vary, so a walk that follows the varying side never sees
  # the invariant `I*I` nested inside the left one. Rotating `(A*B)*C` to `A*(B*C)` puts it where
  # the next turn of the loop peels it.
  #
  # These are measured counts and they are the only witness there is. Switching off either half of
  # the optimisation moves them.
  res <- squares_aev(it_ast(~ .i$male * .x))
  expect_equal(res$instruction_count, 12L)
  expect_equal(res$slot_evaluations, 196L)
})

test_that("a fit distributes one square per diagonal, shared by both triangles", {
  # The recipe builds `X_j * X_j` ONCE and hands the same node to `ewXX` and to `ew2XX`, which
  # differ only in their weight factor. So one rewrite per diagonal moves two references.
  for (k in 1:3) {
    expect_equal(squares_fit(k)$squares_distributed, 2L * k)
  }
})

test_that("the diagonals then cost less than they did, and the answers do not move", {
  # Measured slot evaluations with the optimisation in place. What they are worth is the difference
  # from the numbers without it -- 4881, 9213 and 14625 for one, two and three shaped covariates in
  # the larger fixture this fit is a small version of -- which is a saving of 1080 slot evaluations
  # per term beyond the first, exactly one triangle entry's worth.
  expect_equal(squares_fit(1)$runs[[1]]$slot_evaluations, 230L)
  expect_equal(squares_fit(2)$runs[[1]]$slot_evaluations, 382L)
  expect_equal(squares_fit(3)$runs[[1]]$slot_evaluations, 582L)

  # And the arithmetic still agrees with an A/E, which is a different recipe: with one constant
  # covariate every family collapses onto its three answers.
  plain <- cpp_veil_fit(it_obj(squares_table), list(it_ast(~ 1)), list(0), it_ast(~ .i$amount),
                        NULL, NULL, squares_cols, quarter, NULL, FALSE, 1L)
  reference <- squares_aev(it_ast(~ .i$amount))
  expect_equal(plain$runs[[1]]$E, sum(reference$contributions$E), tolerance = 1e-12)
})

# The three tests below exist because the disable-and-recheck found their guards unwitnessed. Two
# are closed here; the third is recorded in the pass itself as defensive.

test_that("it declines a product of two indicators, which folds outright", {
  # `(I * J)` is itself zero or one, so `(I J)^2` is a square of something already an indicator and
  # `passFoldIndicatorSquares` collapses the whole thing to `I J` in one step. Distributing first
  # would reach the same value by way of five more nodes.
  res <- squares_aev(it_ast(~ .i$male * .i$g1))
  expect_equal(res$squares_distributed, 0L)
  expect_equal(res$instruction_count, 8L)
})

test_that("a scalar duration factor is read as years, not as clicks", {
  # THE CASE THE COERCION EXISTS FOR, and the only one that can see it. `.x` reaches the square as a
  # time vector, which already holds years, so `ToDouble` on it lowers to no instruction at all and
  # dropping it changes nothing. A duration COLUMN is a scalar of clicks -- 534,360 of them to the
  # year -- so squaring it without reading it as years is out by that factor squared.
  #
  # The oracle is the same V reached without any square being distributed: `.i$service * 1` has a
  # literal factor, which this pass declines, so the include route never goes near it.
  weighted <- squares_aev(it_ast(~ .i$male * .i$service))
  included <- squares_aev(it_ast(~ .i$service * 1), include = include(.i$male))

  expect_identical(weighted$contributions$V, included$contributions$V)
  expect_equal(weighted$squares_distributed, 1L)
  expect_equal(included$squares_distributed, 0L)
})

test_that("the click-backed factor is coerced whichever side it is on", {
  # `.x` is a `durationy` and `durationy * durationy` is undefined, so each factor is read as years
  # on its way into the square. Writing the product the other way round puts the age on the LEFT,
  # which is the operand the first version of this pass left uncoerced.
  swapped <- squares_aev(it_ast(~ .x * .i$male))
  written_out <- squares_aev(it_ast(~ .i$male * .x))

  expect_equal(swapped$squares_distributed, 1L)
  expect_identical(swapped$contributions$V, written_out$contributions$V)
  expect_identical(swapped$contributions$E, written_out$contributions$E)
})
