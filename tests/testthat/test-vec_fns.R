# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

create_unit_vector <- function(N, D = 0) ((1:N + D) * 1.6180339887) %% 1
create_signed_vector <- function(N, D = 0) create_unit_vector(N, D) * 2 - 1

# TEMPORARY, 2026-08-31: announce every operation, flushed.
#
# A Windows CI runner dies inside this file with no R error and no testthat
# output at all -- the process simply stops. These lines survive that, so the
# log names the operation and the vector length that killed it.
#
# Remove once the failure is located.
say <- function(what, N) {
  cat("  [vec]", what, "N =", N, "
")
  flush(stdout())
}

test_that("vector unary functions", {

  test_functions <- function(N)
  {
    x <- create_signed_vector(N)
    exp_x <- exp(x)

    say("neg", N)
    expect_identical(vec_neg(x), -x)
    say("exp", N)
    expect_equal(vec_exp(x), exp(x))
    say("expm1", N)
    expect_equal(vec_expm1(x), expm1(x))
    say("log", N)
    expect_equal(vec_log(exp_x), log(exp_x))
    say("log1p", N)
    expect_equal(vec_log1p(x), log1p(x))
  }

  test_functions(1)
  test_functions(3)
  test_functions(4)
  test_functions(5)
  test_functions(7)
  test_functions(8)
  test_functions(31)
  test_functions(32)
  test_functions(33)
  test_functions(1003)
})

test_that("vector binary functions", {

  test_functions <- function(N)
  {
    x <- create_signed_vector(N, 0)
    y <- create_signed_vector(N, 1)
    ux <- create_unit_vector(N, 0)

    say("add", N)
    expect_identical(vec_add(x, y), x + y)
    expect_identical(vec_add(x, 0.1234), x + 0.1234)
    expect_identical(vec_add(-0.2345, y), -0.2345 + y)

    say("sub", N)
    expect_identical(vec_sub(x, y), x - y)
    expect_identical(vec_sub(x, 0.1234), x - 0.1234)
    expect_identical(vec_sub(-0.2345, y), -0.2345 - y)

    say("mul", N)
    expect_identical(vec_mul(x, y), x * y)
    expect_identical(vec_mul(x, 0.1234), x * 0.1234)
    expect_identical(vec_mul(-0.2345, y), -0.2345 * y)

    say("div", N)
    expect_identical(vec_div(x, y), x / y)
    expect_identical(vec_div(x, 0.1234), x / 0.1234)
    expect_identical(vec_div(-0.2345, y), -0.2345 / y)

    # Temperamental
    say("pow", N)
    expect_equal(vec_pow(ux, y), ux ^ y)
    expect_equal(vec_pow(ux, 0.1234), ux ^ 0.1234)
    expect_equal(vec_pow(0.2345, y), 0.2345 ^ y)

    say("min", N)
    expect_identical(vec_min(x, y), pmin(x,y))
    expect_identical(vec_min(x, 0.1234), pmin(x, 0.1234))
    expect_identical(vec_min(-0.2345, y), pmin(-0.2345, y))

    say("max", N)
    expect_identical(vec_max(x, y), pmax(x,y))
    expect_identical(vec_max(x, 0.1234), pmax(x, 0.1234))
    expect_identical(vec_max(-0.2345, y), pmax(-0.2345, y))

    say("clamp", N)
    expect_identical(vec_clamp(x, -0.123, +0.234), pmin(pmax(x,-0.123),+0.234))
  }

  test_functions(1)
  test_functions(3)
  test_functions(4)
  test_functions(5)
  test_functions(7)
  test_functions(8)
  test_functions(31)
  test_functions(32)
  test_functions(33)
  test_functions(1003)
})

#### pow at the edges of the domain ####

# THE OLD TESTS COULD NOT HAVE CAUGHT A WRONG `pow`. They exercised bases in
# (0, 1) and exponents in (-1, 1): no zero, no negative base, no infinity, no
# NaN. When these were written on 2026-09-01 the then-current implementation
# failed SIXTEEN of them on the working AVX-512 tier -- `pow(NaN, Inf)` gave
# Inf, `pow(-2, 0.5)` was not NaN -- so `vec_pow` had been quietly wrong for
# every special value, quite apart from the AVX2 crash that started the hunt.
#
# `vec_pow` FOLLOWS IEEE 754, WHICH R'S `^` DOES NOT. `CLAUDE.md` requires IEEE,
# and the engine already evaluates `Op::Pow` with `std::pow`, so a pronoun
# expression behaves this way too. R's `^` returns NaN for a negative base with
# an infinite exponent where IEEE gives a finite answer or an infinity. The
# eleven disagreements are listed below rather than hidden behind a comparison.

# TEMPORARY, 2026-09-01: the crash moved into these blocks once the original
# tests started passing. Same trick as `say()` above -- flushed, so it
# survives an abrupt termination.
mark <- function(...) {
  cat("  [pow]", ..., "\n")
  flush(stdout())
}

pow_awkward_bases <- c(0, -0, 1, -1, 2, -2, 0.5, -0.5, Inf, -Inf, NaN)
pow_awkward_exponents <- c(0, 1, 2, 3, -1, -2, 0.5, -0.5, Inf, -Inf, NaN)

# Where IEEE and R's `^` part company. Everything else must agree exactly.
pow_ieee_departures <- data.frame(
  x        = c(  -0,  -Inf, -Inf,   -1,   -2, -0.5, -Inf,   -1,   -2, -0.5, -Inf),
  y        = c(  -1,   0.5, -0.5,  Inf,  Inf,  Inf,  Inf, -Inf, -Inf, -Inf, -Inf),
  expected = c(-Inf,   Inf,    0,    1,  Inf,    0,  Inf,    1,    0,  Inf,    0)
)

test_that("pow follows IEEE where R's ^ departs from it", {

  mark("block: ieee departures")

  for (i in seq_len(nrow(pow_ieee_departures))) {
    row <- pow_ieee_departures[i, ]
    label <- paste0("pow(", row$x, ", ", row$y, ")")
    mark("departure", row$x, row$y)

    expect_equal(vec_pow(row$x, row$y), row$expected, info = label)
    expect_equal(vec_pow(c(row$x, row$x), c(row$y, row$y)),
                 c(row$expected, row$expected), info = paste(label, "VV"))
    expect_equal(vec_pow(c(row$x, row$x), row$y),
                 c(row$expected, row$expected), info = paste(label, "VS"))
  }

  # And R really does disagree, so this test is not asserting a tautology.
  expect_true(is.nan((-2)^Inf))
})

# The grid minus the departures: everywhere else `vec_pow` and `^` must agree
# exactly, bit for bit.
pow_agreed_grid <- function() {
  grid <- expand.grid(x = pow_awkward_bases, y = pow_awkward_exponents)
  departure <- mapply(
    function(x, y) any(identical_pair(x, pow_ieee_departures$x) &
                         identical_pair(y, pow_ieee_departures$y)),
    grid$x, grid$y
  )
  grid[!departure, ]
}

# NA IS NOT FALSE, and here it has to be. `value == against` is NA whenever
# either side is NaN, and `NA | FALSE` stays NA -- which then selects NA
# elements out of the grid rather than dropping them.
identical_pair <- function(value, against) {
  same <- (value == against) | (is.nan(value) & is.nan(against))
  same[is.na(same)] <- FALSE
  same
}

test_that("pow matches R's ^ everywhere else, vector against vector", {

  mark("block: grid VV")

  grid <- pow_agreed_grid()
  expect_identical(vec_pow(grid$x, grid$y), grid$x^grid$y)
})

test_that("pow matches R's ^ everywhere else, with a scalar exponent", {

  mark("block: grid VS")

  grid <- pow_agreed_grid()
  for (y in unique(grid$y)) {
    mark("grid VS exponent", y)
    x <- grid$x[identical_pair(grid$y, y)]
    expect_identical(vec_pow(x, y), x^y, info = paste("exponent", y))
  }
})

test_that("the two rules everyone forgets", {

  mark("block: x^0 and 1^y")

  # `x^0` is 1 for every x, and `1^y` is 1 for every y. Both are exactly where
  # exp(y * log x) goes wrong: log(0) is -Inf and 0 * -Inf is NaN, while
  # Inf * log(1) is Inf * 0, also NaN. The scalar-base path guards for both.
  expect_identical(vec_pow(pow_awkward_bases, 0), rep(1, length(pow_awkward_bases)))
  expect_identical(vec_pow(1, pow_awkward_exponents), rep(1, length(pow_awkward_exponents)))
})

test_that("a negative base keeps its sign for integer exponents", {

  mark("block: negative base")

  # `(-2)^3` is -8, not NaN, which exp(y * log x) could never produce.
  expect_identical(vec_pow(c(-2, -2, -2), c(1, 2, 3)), c(-2, 4, -8))
  expect_true(is.nan(vec_pow(-2, 0.5)))
})

#### The one vectorised shape: a scalar base ####

test_that("a scalar base uses exp(y log x) and stays within a few ulp", {

  mark("block: scalar base fast path")

  # THE FAST PATH, and the only place `pow` is vectorised: `1.02 ^ t` over a
  # vector of durations, which is what a projection factor looks like.
  #
  # TOLERANCE, NOT EXACTNESS, and deliberately so. Relative error is about
  # |y log x| * epsilon because error in the exponent becomes relative error
  # after `exp`. Requiring bit-equality here would pass on this machine and
  # break on another.
  y <- c(-30, -3.5, -1, -0.25, 0, 0.25, 1, 3.5, 30)

  for (x in c(1.02, 0.5, 2, 1e-3, 1e3)) {
    mark("fast path base", x)
    expect_equal(vec_pow(x, y), x^y, tolerance = 1e-12, info = paste("base", x))
  }
})

test_that("a scalar base of one is exact for every exponent, including infinite", {

  mark("block: base one")

  # Needs its own branch: log(1) is 0, and Inf * 0 is NaN, where pow(1, Inf) is
  # 1. Filled directly rather than computed.
  expect_identical(vec_pow(1, c(0, 1, Inf, -Inf, NaN)), rep(1, 5))
})

test_that("the scalar-base fast path is exact for special exponents", {

  mark("block: fast path specials")

  # Inside the guard -- x finite, positive, not 1 -- every special exponent
  # comes out right through exp(y log x), with no special-casing needed.
  expect_identical(vec_pow(2, c(0, Inf, -Inf, NaN)), c(1, Inf, 0, NaN))
  expect_identical(vec_pow(0.5, c(0, Inf, -Inf, NaN)), c(1, 0, Inf, NaN))
})
