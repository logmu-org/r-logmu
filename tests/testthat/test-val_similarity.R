# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests for `val_similarity` and `val_distance`, the second weighting factor.
#
# THE DEFINING PROPERTY, and what most of these pin, is that the similarity is
# LINEAR IN ALL THREE OUTPUTS while the weight is squared in V. Halving the
# similarity halves V; halving the weight quarters it. Anything that folded the
# similarity into the weight would pass a test of E alone and fail here.

data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955)),
    pension   = c(5000, 12000, 30000),
    male      = c(TRUE, FALSE, TRUE),
    fraction  = c(0.2, 0.5, 0.9),
    E2R_start = datey::datey(rep(2015, 3)),
    E2R_end   = datey::datey(c(2020, 2020, 2018)),
    E2R_died  = c(FALSE, FALSE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

flat <- mortality_const(log_mu = -4)

parts <- function(x) c(x$A, x$E, x$V)

test_that("similarity scales A, E and V alike", {
  plain <- aev(data, mortality = flat, overdispersion = 1)
  half  <- aev(data, mortality = flat, overdispersion = 1, val_similarity = 0.5)

  expect_equal(parts(half), parts(plain) / 2)
})

test_that("a distance is the same quantity written the other way round", {
  by_similarity <- aev(data, mortality = flat, overdispersion = 1, val_similarity = 0.5)
  by_distance  <- aev(data, mortality = flat, overdispersion = 1, val_distance = log(2))

  expect_equal(parts(by_distance), parts(by_similarity))
})

test_that("similarity is linear in V where the weight is squared", {
  # THE TEST THAT SEPARATES THE TWO FACTORS. Doubling the weight doubles E and
  # quadruples V. Halving the similarity halves both.
  plain    <- aev(data, mortality = flat, overdispersion = 1, weight = .i$pension)
  weighted <- aev(data, mortality = flat, overdispersion = 1, weight = 2 * .i$pension)
  similar  <- aev(data, mortality = flat, overdispersion = 1, weight = .i$pension,
                  val_similarity = 0.5)

  expect_equal(weighted$E, 2 * plain$E)
  expect_equal(weighted$V, 4 * plain$V)

  expect_equal(similar$E, plain$E / 2)
  expect_equal(similar$V, plain$V / 2)
})

test_that("an indicator weight still gives V = E when a similarity is applied", {
  # THE WITNESS THAT THE SIMILARITY IS NOT FOLDED INTO THE WEIGHT. V collapses
  # into E because `w * w` simplifies to `w` for an indicator. Multiply the
  # similarity in first and `0.5 * male` is no longer an indicator, so its square
  # would be a quarter and V would come back at half of E.
  indicator <- aev(data, mortality = flat, overdispersion = 1, weight = .i$male,
                   val_similarity = 0.5)

  expect_equal(indicator$V, indicator$E)
})

test_that("a similarity of one, or a distance of zero, changes nothing", {
  plain <- aev(data, mortality = flat, overdispersion = 1, weight = .i$pension)

  expect_equal(parts(aev(data, mortality = flat, overdispersion = 1,
                         weight = .i$pension, val_similarity = 1)),
               parts(plain))
  expect_equal(parts(aev(data, mortality = flat, overdispersion = 1,
                         weight = .i$pension, val_distance = 0)),
               parts(plain))
})

test_that("a distance may vary with time, and A takes its value at the death", {
  # The only death is at 2018, so its contribution to A is exp(-d) evaluated
  # there -- which is what pins that the factor reaches `died_value` rather than
  # being applied to the integral alone.
  decaying <- aev(data, mortality = flat, overdispersion = 1,
                  val_distance = (2020 - .t) / 10)

  expect_equal(decaying$A, exp(-0.2))

  # And a varying distance is not the same as a constant one, or the test above
  # would hold for the wrong reason.
  expect_false(isTRUE(all.equal(
    decaying$E,
    aev(data, mortality = flat, overdispersion = 1, val_distance = 0.2)$E
  )))
})

test_that("similarity reaches every element of a breakdown", {
  plain <- aev(data, mortality = flat, overdispersion = 1,
               breakdown = bands(.i$pension, thresholds = 10000))
  half  <- aev(data, mortality = flat, overdispersion = 1, val_similarity = 0.5,
               breakdown = bands(.i$pension, thresholds = 10000))

  expect_equal(half$E, plain$E / 2)
  expect_equal(half$V, plain$V / 2)
})

test_that("the two spellings are one quantity, so only one may be given", {
  expect_error(
    aev(data, mortality = flat, overdispersion = 1, val_similarity = 0.5, val_distance = 1),
    "not both"
  )
})

# ---- the bound, checked analytically ---------------------------------------

test_that("a similarity certainly outside [0, 1] is refused before any record is read", {
  expect_error(aev(data, mortality = flat, overdispersion = 1, val_similarity = 2),
               "always above 1")
  expect_error(aev(data, mortality = flat, overdispersion = 1, val_similarity = -0.5),
               "always below 0")

  # Said in the spelling the user wrote: a negative distance is a similarity
  # above one, but they did not write a similarity.
  expect_error(aev(data, mortality = flat, overdispersion = 1, val_distance = -1),
               "`val_distance` cannot be negative")
})

test_that("the bound is checked against the column's scanned range, not just literals", {
  # `pension` runs from 5000 to 30000, so no record could ever give a similarity
  # in range. That is settled from the scan, with nothing evaluated.
  expect_error(aev(data, mortality = flat, overdispersion = 1, val_similarity = .i$pension),
               "always above 1")

  # And a column that does lie in range passes, or the test above would hold for
  # the wrong reason.
  expect_silent(aev(data, mortality = flat, overdispersion = 1, val_similarity = .i$fraction))
})

test_that("a range that merely permits a violation is accepted", {
  # THE RESTRAINT THAT MAKES THE CHECK USABLE. The range of `.t` alone allows a
  # negative distance here, because it carries every representable date rather
  # than the span of this data. Refusing anything that might violate would
  # reject the ordinary decay kernel, and a check that fires on correct code is
  # worse than one that stays quiet on incorrect code.
  expect_silent(aev(data, mortality = flat, overdispersion = 1,
                    val_distance = (2025 - .t) / 10))
})

test_that("the range of `.t` is the data's exposure, so a bad kernel is caught", {
  # Over exposure running 2015 to 2020 this distance is negative at every
  # instant, and that is settled from the exposure columns before a record is
  # read. It is only knowable because `.t` is bounded by the DATA -- against the
  # representable calendar the expression could be positive somewhere.
  expect_error(aev(data, mortality = flat, overdispersion = 1,
                   val_distance = (2000 - .t) / 10),
               "`val_distance` cannot be negative")
})

test_that("the same kernel is accepted over data it suits", {
  # THE WITNESS THAT THE BOUND COMES FROM THE DATA rather than from anything
  # written in the expression: the identical distance is fine here, because this
  # exposure ends before 2000.
  earlier <- exp_data(
    list(
      birth     = datey::datey(1930),
      E2R_start = datey::datey(1990),
      E2R_end   = datey::datey(1995),
      E2R_died  = FALSE
    ),
    exp_start = datey::datey(1990),
    exp_end   = datey::datey(1995)
  )

  expect_silent(aev(earlier, mortality = flat, overdispersion = 1,
                    val_distance = (2000 - .t) / 10))
})

# ---- through a batch -------------------------------------------------------

test_that("similarity and distance default from a batch", {
  b <- batch(
    .exp_data = data, .overdispersion = 1, .val_distance = log(2),
    inherited = aev(mortality = flat)
  )

  expect_identical(b$inherited,
                   aev(data, mortality = flat, overdispersion = 1, val_distance = log(2)))
})

test_that("an analysis naming either spelling takes neither batch default", {
  # Injected argument by argument, a `.val_distance` default would join the
  # analysis's own `val_similarity` and make the pair that is refused -- a confusing
  # way to be told a default exists.
  b <- batch(
    .exp_data = data, .overdispersion = 1, .val_distance = log(2),
    own = aev(mortality = flat, val_similarity = 0.25)
  )

  expect_identical(b$own,
                   aev(data, mortality = flat, overdispersion = 1, val_similarity = 0.25))
})

test_that("a batched analysis giving both spellings is named in the complaint", {
  # The engine refuses the pair as well, and its message is the same words. What
  # only the batch can say is WHOSE fault it is, which is the whole reason for a
  # second check rather than leaving it to the boundary.
  expect_error(
    batch(.exp_data = data, .overdispersion = 1,
          fine    = aev(mortality = flat),
          muddled = aev(mortality = flat, val_similarity = 0.5, val_distance = 1)),
    "`muddled` gives both"
  )
})

test_that("a batch may not give both spellings either", {
  expect_error(
    batch(.exp_data = data, .overdispersion = 1, .val_similarity = 0.5, .val_distance = 1,
          a = aev(mortality = flat)),
    "not both"
  )
})
