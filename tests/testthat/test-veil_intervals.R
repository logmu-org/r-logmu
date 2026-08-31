# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests interval propagation: the range of values each node can take, worked out from the literals and
# scanned column ranges upwards through the operators. The pass annotates only -- nothing is rewritten
# -- so these tests read the reported bounds rather than the shape of the tree.
#
# Intervals are in each node's OWN units, so a datey or durationy node reports CLICKS, not years.

clicks_per_year <- 534360

# A dataset's columns must all be the same length, as they would be in a data frame.
cols <- list(
  age     = c(30, 45, 60),
  gappy   = c(1, NA, 3),
  pos     = c(2, 4, 5),    # strictly positive, so it may be divided by
  neg     = c(-5, -4, -2), # strictly negative, likewise
  spans   = c(-2, 1, 3),   # straddles zero, so a quotient has a gap in it
  touches = c(0, 1, 2),    # reaches zero at one end, which is the trap
  birth   = datey::datey(c(1960, 1965, 1975))
)

root_lo <- function(res) res$lo[[res$root + 1L]]
root_hi <- function(res) res$hi[[res$root + 1L]]

test_that("a scanned column reports its own range", {
  res <- cpp_veil_intervals(it_ast(~ .i$age), cols)
  expect_equal(root_lo(res), 30)
  expect_equal(root_hi(res), 60)
})

test_that("a literal is a single point", {
  res <- cpp_veil_intervals(it_ast(~ 7), cols)
  expect_equal(root_lo(res), 7)
  expect_equal(root_hi(res), 7)
})

test_that("a column with NAs gives nothing away", {
  res <- cpp_veil_intervals(it_ast(~ .i$gappy), cols)
  expect_equal(root_lo(res), -Inf)
  expect_equal(root_hi(res), Inf)
})

test_that("arithmetic propagates through a product", {
  res <- cpp_veil_intervals(it_ast(~ .i$age * 2), cols)
  expect_equal(root_lo(res), 60)
  expect_equal(root_hi(res), 120)
})

test_that("a sum of two columns adds their ranges", {
  res <- cpp_veil_intervals(it_ast(~ .i$age + .i$age), cols)
  expect_equal(root_lo(res), 60)
  expect_equal(root_hi(res), 120)
})

test_that("a datey column reports clicks, not years", {
  res <- cpp_veil_intervals(it_ast(~ .i$birth), cols)
  expect_equal(root_lo(res), 1960 * clicks_per_year)
  expect_equal(root_hi(res), 1975 * clicks_per_year)
})

test_that("a duration between two dateys stays in clicks", {
  # `.t` is only known to lie in the representable calendar, 1000 to 3000, so the duration spans
  # that against birth's scanned range. Both ends are clicks.
  res <- cpp_veil_intervals(it_ast(~ .t - .b), cols)
  expect_equal(root_lo(res), 1000 * clicks_per_year - 1975 * clicks_per_year)
  expect_equal(root_hi(res), 3000 * clicks_per_year - 1960 * clicks_per_year)
})

test_that("converting a duration to a number rescales clicks to years", {
  # Multiplying by a plain number widens the duration to a double, which is the one place the pass
  # rescales. Halving the age in years is the clearest way to see it.
  res <- cpp_veil_intervals(it_ast(~ (.t - .b) * 0.5), cols)
  expect_equal(root_lo(res), 0.5 * (1000 - 1975))
  expect_equal(root_hi(res), 0.5 * (3000 - 1960))
})

test_that("a comparison is confined to 0 and 1", {
  res <- cpp_veil_intervals(it_ast(~ .i$age > 40), cols)
  expect_equal(root_lo(res), 0)
  expect_equal(root_hi(res), 1)
})

test_that("a selection spans both of its branches", {
  res <- cpp_veil_intervals(it_ast(~ if (.i$age > 40) 2 else 9), cols)
  expect_equal(root_lo(res), 2)
  expect_equal(root_hi(res), 9)
})

test_that("clamp is bounded by the limits it was given", {
  res <- cpp_veil_intervals(it_ast(~ clamp(.i$age, 0, 1)), cols)
  expect_equal(root_lo(res), 0)
  expect_equal(root_hi(res), 1)
})

test_that("exp of a bounded range stays bounded", {
  res <- cpp_veil_intervals(it_ast(~ exp(.i$age)), cols)
  expect_equal(root_lo(res), exp(30))
  expect_equal(root_hi(res), exp(60))
})

test_that("dividing by a literal scales the range", {
  res <- cpp_veil_intervals(it_ast(~ .i$age / 10), cols)
  expect_equal(root_lo(res), 3)
  expect_equal(root_hi(res), 6)

  # A negative divisor turns the interval over.
  res <- cpp_veil_intervals(it_ast(~ .i$age / -10), cols)
  expect_equal(root_lo(res), -6)
  expect_equal(root_hi(res), -3)
})

test_that("dividing by a column takes the extreme corners", {
  # age is 30 to 60 and pos is 2 to 5, so the quotient runs from 30/5 to 60/2.
  res <- cpp_veil_intervals(it_ast(~ .i$age / .i$pos), cols)
  expect_equal(root_lo(res), 6)
  expect_equal(root_hi(res), 30)

  # neg is -5 to -2, so the extremes are 60/-2 and 30/-5.
  res <- cpp_veil_intervals(it_ast(~ .i$age / .i$neg), cols)
  expect_equal(root_lo(res), -30)
  expect_equal(root_hi(res), -6)
})

test_that("a divisor that reaches zero gives nothing away", {
  # Straddling zero puts a gap in the quotient, which one interval cannot express.
  res <- cpp_veil_intervals(it_ast(~ .i$age / .i$spans), cols)
  expect_equal(root_lo(res), -Inf)
  expect_equal(root_hi(res), Inf)

  # Merely TOUCHING zero is refused too. Reaching it from below is the case that would be wrong
  # rather than loose -- an interval endpoint written 0 is a POSITIVE zero, so dividing by it gives
  # +Inf where the values just inside tend to -Inf, and the corner rule would claim a bound that
  # omits them. Reaching it from above, as `touches` does, would in fact be safe; the guard declines
  # it anyway, because an asymmetric test buys little and invites the sign mistake.
  res <- cpp_veil_intervals(it_ast(~ .i$age / .i$touches), cols)
  expect_equal(root_lo(res), -Inf)
  expect_equal(root_hi(res), Inf)
})

test_that("round rounds the bounds rather than bracketing them", {
  # `round` is banker's rounding, so the bound is the rounding applied to each end -- sound because
  # rounding is monotonic. The looser floor-and-ceiling bracket this replaced gave 2 and 6.
  # age runs 30 to 60, so age * 0.09 runs 2.7 to 5.4. (Multiplied rather than divided: division has
  # no interval rule yet, so a quotient is unknown and would tell us nothing about rounding.)
  res <- cpp_veil_intervals(it_ast(~ round(.i$age * 0.09)), cols)
  expect_equal(root_lo(res), 3)
  expect_equal(root_hi(res), 5)

  # An exact half at the top end lands on the even side: 60 * 0.125 is 7.5, which rounds to 8. The
  # bottom end, 3.75, is where the bracket would have said 3.
  res <- cpp_veil_intervals(it_ast(~ round(.i$age * 0.125)), cols)
  expect_equal(root_lo(res), 4)
  expect_equal(root_hi(res), 8)
})

test_that("round leaves an infinite bound alone", {
  # roundBankers is finite-only -- an infinity reaches std::fmod and comes back NaN -- so an infinite
  # end passes through untouched. Rounding it would collapse the whole interval to unknown, losing
  # the finite end as well. `max` against zero is the way to get a half-bounded interval here: gappy
  # has NAs, so it gives nothing away, and the maximum lifts the bottom to zero.
  res <- cpp_veil_intervals(it_ast(~ round(max(.i$gappy, 0))), cols)
  expect_equal(root_lo(res), 0)
  expect_equal(root_hi(res), Inf)
})

test_that("log is bounded only when its argument is known positive", {
  # age is known to start at 30, so the logarithm is bounded. An unknown argument could reach zero
  # or below, where the logarithm is undefined, so the pass declines to bound it.
  expect_equal(root_lo(cpp_veil_intervals(it_ast(~ log(.i$age)), cols)), log(30))
  expect_equal(root_lo(cpp_veil_intervals(it_ast(~ log(.i$gappy)), cols)), -Inf)
})
