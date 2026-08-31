# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests lowering and the interpreter -- the scalar spine. An expression is compiled the whole way
# down to three-address code and run for every individual, and the answer is checked against R
# evaluating the same expression over the same columns.
#
# R IS THE ORACLE HERE. That is the point of doing the scalar spine first: for a time-invariant value
# per individual, R computes the expected answer itself, so these tests check veil against R rather
# than against a number worked out by hand.
#
# Two places where veil differs from R on purpose, and which the tests below therefore steer around:
#
#   - A comparison is never missing. veil follows IEEE 754 and draws no distinction between R's
#     NA_real_ and any other NaN, so `NaN == NaN` is false where R gives NA.
#   - Every individual gets a value, even from an expression that reads no column. R evaluating a
#     constant expression gives one value, so the oracle needs rep().

# A dataset's columns must all be the same length, as they would be in a data frame.
cols <- list(
  age    = c(30, 45, 60),
  weight = c(1.5, 2.5, 0.5),
  fixed  = c(7, 7, 7),                        # constant, so the tree phase folds it away
  gappy  = c(1, NA, 3),
  sign   = c(1, -1, 1),                       # multiplied by zero, gives +0 and -0
  birth  = datey::datey(c(1960, 1965, 1975))
)

records <- length(cols$age)

values_of <- function(expr) cpp_veil_eval(expr, cols)$values

# Asserts that veil's answer matches R's own for the same expression. `expected` is evaluated by R.
expect_matches_r <- function(expr, expected) {
  expect_equal(values_of(expr), expected)
}

# The same check, but with `it_eval` -- the package's own R evaluator of the same AST -- supplying
# the expected answer. It is the stronger of the two: it reads the very node tree that crossed into
# C++, so it cannot drift from the expression veil compiled the way a transcribed one can. A datey
# or durationy result comes back from R classed and from veil as bare clicks, which is what the
# unclass is for.
expect_matches_it_eval <- function(expr) {
  expect_equal(values_of(expr), unclass(it_eval(expr, .i = cols)))
}

test_that("a sweep of expressions agrees with R's own evaluator", {
  # Every expression here reads at least one column, so R returns one value per individual as veil
  # does. The selection is left out: it_eval gives `if` its scalar per-individual meaning, which
  # cannot be evaluated over a whole column at once.
  sweep <- list(
    it_ast(~ .i$age),
    it_ast(~ .i$age * 2 + 1),
    it_ast(~ .i$age - .i$weight),
    it_ast(~ .i$age / .i$weight),
    it_ast(~ -.i$age),
    it_ast(~ abs(.i$weight - 2)),
    it_ast(~ exp(.i$age / 100)),
    it_ast(~ log(.i$age)),
    it_ast(~ sqrt(.i$age)),
    it_ast(~ log10(.i$age)),
    it_ast(~ expm1(.i$weight)),
    it_ast(~ log1p(.i$weight)),
    it_ast(~ sin(.i$weight) + cos(.i$weight)),
    it_ast(~ .i$weight^3),
    it_ast(~ .i$age %/% 7),
    it_ast(~ .i$age %% 7),
    it_ast(~ floor(.i$age / 7)),
    it_ast(~ ceiling(.i$age / 7)),
    it_ast(~ trunc(.i$age / 7)),
    it_ast(~ sign(.i$weight - 2)),
    it_ast(~ min(.i$age, 40)),
    it_ast(~ max(.i$age, 40)),
    it_ast(~ clamp(.i$age, 40, 50)),
    it_ast(~ .i$age > 40),
    it_ast(~ .i$age <= 45),
    it_ast(~ .i$weight == 2.5),
    it_ast(~ .i$age > 40 & .i$weight > 1),
    it_ast(~ .i$age > 40 | .i$weight > 2),
    it_ast(~ !(.i$age > 40)),
    it_ast(~ is.na(.i$gappy)),
    it_ast(~ .b),
    it_ast(~ .b + datey::durationy(10)),
    it_ast(~ .b > 1968),
    it_ast(~ .b <= 1965)
  )
  for (expr in sweep) {
    expect_matches_it_eval(expr)
  }
})

# The sweep above trusts it_ast to have built the tree the source text describes. These write the
# expected answer out by hand instead, so that a fault in the AST builder itself would show up as a
# disagreement rather than being shared by both sides.

test_that("arithmetic over a column matches R written out by hand", {
  expect_matches_r(it_ast(~ .i$age * 2 + 1), with(cols, age * 2 + 1))
  expect_matches_r(it_ast(~ .i$age / .i$weight), with(cols, age / weight))
  expect_matches_r(it_ast(~ exp(.i$age / 100)), with(cols, exp(age / 100)))
  expect_matches_r(it_ast(~ .i$age %% 7), with(cols, age %% 7))
  expect_matches_r(it_ast(~ min(.i$age, 40)), with(cols, pmin(age, 40)))
  expect_matches_r(it_ast(~ clamp(.i$age, 40, 50)), with(cols, pmin(pmax(age, 40), 50)))
  expect_matches_r(it_ast(~ .i$age > 40 & .i$weight > 1), with(cols, age > 40 & weight > 1))
})

test_that("rounding is banker's rounding, as R's round is", {
  expect_equal(values_of(it_ast(~ round(2.5))), rep(2, records))
  expect_equal(values_of(it_ast(~ round(3.5))), rep(4, records))
  expect_equal(values_of(it_ast(~ round(-2.5))), rep(-2, records))
})

test_that("round leaves an infinity alone", {
  # IEEE 754-2019 roundToIntegralTiesToEven returns an infinity unchanged, and so does R. The
  # division is by a column rather than by a literal so that R cannot fold it away before it
  # crosses, which is what makes this reach the interpreter.
  expect_matches_r(it_ast(~ round(.i$weight / 0)), with(cols, round(weight / 0)))
  expect_matches_r(it_ast(~ round(-.i$weight / 0)), with(cols, round(-weight / 0)))
  expect_equal(is.na(values_of(it_ast(~ round(.i$gappy / 0)))), c(FALSE, TRUE, FALSE))
})

test_that("the interval round is given agrees with the values round produces", {
  # Interval propagation bounds `round` by rounding each end, which is only sound because the engine
  # rounds the same way. These pairs are the coupling: the second number is what the pass believes
  # the range to be, and the values are what the interpreter actually computes for the same
  # expression, checked against R. A change to either that broke the agreement would show up here
  # rather than as a comparison folded to the wrong answer.
  expect_equal(values_of(it_ast(~ round(.i$age * 0.09))), round(cols$age * 0.09))
  expect_equal(range(values_of(it_ast(~ round(.i$age * 0.09)))), c(3, 5))

  expect_equal(values_of(it_ast(~ round(.i$age * 0.125))), round(cols$age * 0.125))
  expect_equal(range(values_of(it_ast(~ round(.i$age * 0.125)))), c(4, 8))
})

test_that("a quotient carries an interval, so banding by decade folds", {
  # This is the spelling the round tightening could not reach until division propagated an interval:
  # age is 30 to 60, so age/10 is [3, 6] and rounding it stays [3, 6]. 45/10 is an exact half, and
  # rounds down to 4 on the even side -- R agrees, which is what the first assertion checks.
  expect_matches_r(it_ast(~ round(.i$age / 10)), with(cols, round(age / 10)))

  res <- cpp_veil_eval(it_ast(~ round(.i$age / 10) <= 6), cols)
  expect_equal(res$values, rep(TRUE, records))
  expect_equal(length(res$monikers), 0L)

  res <- cpp_veil_eval(it_ast(~ .i$age / 10 > 2), cols)
  expect_equal(res$values, rep(TRUE, records))
  expect_equal(length(res$monikers), 0L)
})

test_that("a quotient whose divisor can reach zero still runs", {
  # weight is 0.5 to 2.5, so it never reaches zero and the quotient is bounded. gappy has NAs, so it
  # gives nothing away and the division declines -- the comparison is left to run per individual.
  res <- cpp_veil_eval(it_ast(~ .i$age / .i$gappy > 2), cols)
  expect_true(length(res$monikers) > 0L)
  # The middle individual divides by a missing value, and a veil comparison against a NaN is FALSE
  # rather than missing -- the IEEE unordered rule -- which is what the `& !is.na` spells out.
  expect_equal(res$values, with(cols, (age / gappy > 2) & !is.na(gappy)))
})

test_that("the tightened round bound settles a comparison the looser one could not", {
  # round(age * 0.09) lies in [3, 5]. Under the floor-and-ceiling bracket this replaced the bound was
  # [2, 6], which straddles both thresholds below and left them to run per individual.
  res <- cpp_veil_eval(it_ast(~ round(.i$age * 0.09) < 6), cols)
  expect_equal(res$values, rep(TRUE, records))
  expect_equal(length(res$monikers), 0L)

  res <- cpp_veil_eval(it_ast(~ round(.i$age * 0.09) > 2), cols)
  expect_equal(res$values, rep(TRUE, records))
  expect_equal(length(res$monikers), 0L)
})

test_that("a selection matches ifelse", {
  # `if` in R is not vectorised, so the oracle is the vectorised spelling of the same choice.
  expect_matches_r(it_ast(~ if (.i$age > 40) .i$weight else 0),
                   with(cols, ifelse(age > 40, weight, 0)))
})

test_that("is.na matches R where the answer is not settled at compile time", {
  expect_matches_r(it_ast(~ is.na(.i$gappy)), with(cols, is.na(gappy)))
})

test_that("a datey comes back as clicks and matches R's own arithmetic", {
  # This entry point hands back bare clicks with no class attached, so the oracle is unclassed.
  expect_equal(values_of(it_ast(~ .b)), as.integer(unclass(cols$birth)))
  expect_equal(values_of(it_ast(~ .b + datey::durationy(10))),
               as.integer(unclass(cols$birth + datey::durationy(10))))
})

test_that("a duration read as a number is in years", {
  # Multiplying by a plain number widens the duration to a double, and a durationy reads as years.
  expect_matches_r(it_ast(~ (.b - datey::datey(1900)) * 1),
                   as.numeric(cols$birth - datey::datey(1900)))
})

test_that("comparing a datey against a plain number matches R", {
  # A plain number is a number of YEARS, and narrowing rewrites it into a click threshold long
  # before this runs. R compares by converting the date to years instead. This is the check that the
  # two routes reach the same answer, which is the whole reason narrowing has to be exact.
  expect_matches_r(it_ast(~ .b > 1968), cols$birth > 1968)
  expect_matches_r(it_ast(~ .b <= 1965), cols$birth <= 1965)
  expect_matches_r(it_ast(~ .b > datey::datey(1968)), cols$birth > datey::datey(1968))
})

# min, max and clamp follow IEEE 754-2019 `minimum` and `maximum`, which R's pmin and pmax do not.
# These assert against the standard, so R is NOT the oracle here -- it is the thing being diverged
# from, deliberately, and the comments say what R gives instead.

# A signed zero is only observable through something that reads its sign, and a reciprocal is the
# usual route: +0 gives Inf, -0 gives -Inf. The zeros are produced by arithmetic rather than written
# as literals, so that R's own constant folding cannot settle the answer before it ever crosses.
reciprocal_of <- function(expr) 1 / values_of(expr)

test_that("min and max order the signed zeros, whichever way round the arguments come", {
  # 754-2019: minimum(+0, -0) is -0 both ways round. R's pmin returns its FIRST argument on a tie,
  # so pmin(0, -0) gives +0 -- order-dependent, and wrong in one direction.
  expect_equal(reciprocal_of(it_ast(~ min(.i$sign * 0, 0))), c(Inf, -Inf, Inf))
  expect_equal(reciprocal_of(it_ast(~ min(0, .i$sign * 0))), c(Inf, -Inf, Inf))

  # 754-2019: maximum(+0, -0) is +0 both ways. R's pmax(-0, 0) gives -0.
  expect_equal(reciprocal_of(it_ast(~ max(.i$sign * 0, 0))), c(Inf, Inf, Inf))
  expect_equal(reciprocal_of(it_ast(~ max(0, .i$sign * 0))), c(Inf, Inf, Inf))
})

test_that("clamp is the composition of maximum and minimum", {
  # clamp is veil's own, so this is a choice rather than a standard: clamp(x, lo, hi) means
  # minimum(maximum(x, lo), hi). Defining it that way keeps it from becoming a third set of NaN and
  # signed-zero rules to hold in step with the other two.
  expect_equal(reciprocal_of(it_ast(~ clamp(.i$sign * 0, 0, 1))), c(Inf, Inf, Inf))
  expect_equal(values_of(it_ast(~ clamp(.i$age, 40, 50))), c(40, 45, 50))
})

test_that("min and max propagate a missing value from either argument", {
  # This is 754-2019 `minimum`, not a concession to R: the 2008 edition's minNum -- which is what
  # C's fmin implements -- discarded the NaN and returned 2. R happens to agree with 2019 here.
  expect_equal(is.na(values_of(it_ast(~ min(.i$gappy, 2)))), c(FALSE, TRUE, FALSE))
  expect_equal(is.na(values_of(it_ast(~ min(2, .i$gappy)))), c(FALSE, TRUE, FALSE))
  expect_equal(is.na(values_of(it_ast(~ max(.i$gappy, 2)))), c(FALSE, TRUE, FALSE))
  expect_equal(is.na(values_of(it_ast(~ clamp(.i$gappy, 0, 2)))), c(FALSE, TRUE, FALSE))
})

test_that("click min and max need none of that, and are plain comparisons", {
  # An integer has no NaN and no signed zero, so the comparison already IS 754-2019 minimum. This is
  # the case the include clip uses, which is why full compliance costs nothing where it matters.
  expect_equal(values_of(it_ast(~ min(.b, datey::datey(1965)))),
               as.integer(pmin(unclass(cols$birth), unclass(datey::datey(1965)))))
  expect_equal(values_of(it_ast(~ max(.b, datey::datey(1965)))),
               as.integer(pmax(unclass(cols$birth), unclass(datey::datey(1965)))))
})

# What the expression lowered to. A block is a table of operands plus a straight-line body, so these
# read the shape of the compilation rather than its answer.

test_that("a constant column is folded away, leaving nothing to run", {
  res <- cpp_veil_eval(it_ast(~ .i$fixed), cols)
  expect_equal(res$values, rep(7, records))
  expect_equal(length(res$monikers), 0L)   # nothing to compute
  expect_equal(res$column_count, 0L)       # and no column to read
  expect_equal(res$constant_count, 1L)
})

test_that("a comparison the data settles is folded away too", {
  res <- cpp_veil_eval(it_ast(~ .i$age > 200), cols)
  expect_equal(res$values, rep(FALSE, records))
  expect_equal(length(res$monikers), 0L)
})

test_that("a column read twice is read once", {
  res <- cpp_veil_eval(it_ast(~ .i$age + .i$age), cols)
  expect_equal(res$values, with(cols, age + age))
  expect_equal(res$column_count, 1L)
  expect_equal(res$monikers, "add")
})

test_that("the block reports the operations it holds, in order", {
  res <- cpp_veil_eval(it_ast(~ exp(.i$age / 100)), cols)
  expect_equal(res$monikers, c("rdiv", "exp"))
  expect_equal(res$type, "double")
})

test_that("a widening shows up as an explicit conversion", {
  # `.b * 1` widens the datey to a number, which lowering makes an instruction rather than leaving
  # the interpreter to guess how a click becomes a year.
  res <- cpp_veil_eval(it_ast(~ .b * 1), cols)
  expect_equal(res$monikers, c("double", "mul"))
})

# What the scalar spine does not do yet. These are refused with a message rather than half-handled.

test_that("a time-varying expression is refused", {
  # These columns carry no exposure, so there is nowhere for a time vector to come from and that is
  # what the message says. `test-veil_time.R` covers the other half: with an exposure present, the
  # complaint becomes the useful one, that such an expression has no single value per individual.
  expect_error(cpp_veil_eval(it_ast(~ .t), cols), "samples time")
  expect_error(cpp_veil_eval(it_ast(~ .x), cols), "samples time")
})

test_that("a mortality table needs an exposure to be read over", {
  # The table lowers perfectly well now -- see test-veil_log_mu.R -- so what stops it here is that
  # this data carries no exposure for the time vector to be built from. The complaint names that,
  # rather than the table.
  log_mu <- matrix(seq(-5, -3, length.out = 12), nrow = 3, ncol = 4)
  tbl <- mortality_table(x0 = 60, t0 = 2010, log_mu = log_mu)
  expect_error(cpp_veil_eval(it_obj(tbl), cols), "no exposure was supplied")
})

test_that("an include is refused", {
  expect_error(cpp_veil_eval(it_obj(age(65, 95)), cols), "cannot yet be lowered")
})

test_that("a constant mortality is just a number, and does run", {
  # It lowers to a literal rather than to anything needing the time vector, so nothing refuses it.
  m <- mortality_const(log_mu = -4.5)
  expect_equal(values_of(it_obj(m)), rep(-4.5, records))
})
