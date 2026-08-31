# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the A/E/V recipe: one A, one E and one V out of one pass over the data.
#
#     A = died_value(w)        the weight at the moment of death, nothing for a survivor
#     E = integrate(mu * w)    expected deaths
#     V = integrate(mu * w^2)  the second moment
#
# `res$A` IS THE TOTAL, summed in the core. `res$contributions$A` is the per-individual vector, which
# is diagnostic and goes when the batch entry point lands. Nearly every assertion below is written
# against the contributions, because that is where the analytic oracles live: each individual has
# their own exposure and their own age, so each has their own answer worked out on paper. A total that
# disagrees says only that something somewhere is wrong.
#
# THE ORACLE IS WORKED OUT ON PAPER, not computed a second way. With a constant log mu the integrand
# is constant in time, and with a weight linear in time it is linear, and the midpoint rule is exact
# for both -- so E is `mu * w * exposure` or an integral of a straight line, and V follows. That keeps
# the expected answers independent of everything the engine does to reach them.
#
# Exposures are whole quarters, so there is no short final interval and none of the half-click
# rounding test-veil_time.R has to allow for.

clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

birth_years <- c(1940, 1945, 1950)
start_clicks <- c(2010, 2010, 2010) * clicks_per_year
end_clicks <- start_clicks + c(12L, 8L, 4L) * quarter

cols <- list(
  birth     = datey::datey(birth_years),
  amount    = c(1000, 2500, 400),
  male      = c(TRUE, FALSE, TRUE),
  E2R_start = datey_clicks(start_clicks),
  E2R_end   = datey_clicks(end_clicks),
  E2R_died  = c(TRUE, FALSE, TRUE)
)

exposure_years <- (end_clicks - start_clicks) / clicks_per_year
died <- cols$E2R_died

log_mu_value <- -3.2
mu_value <- exp(log_mu_value)

constant_mortality <- mortality_const(log_mu = log_mu_value)

aev <- function(mortality = it_obj(constant_mortality), weight = NULL, include = NULL,
                time_scale = quarter_scale, columns = cols,
                overdispersion = no_overdispersion, threads = 1L) {
  cpp_veil_aev(mortality, weight, columns, time_scale, include, overdispersion, threads)
}

test_that("a count of lives is exposure and deaths", {
  # No weight at all, so w is one: A counts deaths, E and V are both mu times the exposure.
  res <- aev()

  expect_equal(res$contributions$A, as.double(died), tolerance = 1e-12)
  expect_equal(res$contributions$E, mu_value * exposure_years, tolerance = 1e-12)
  expect_equal(res$contributions$V, mu_value * exposure_years, tolerance = 1e-12)
})

test_that("a concept object handed in as an obj leaf is lowered", {
  # THE R FRONT END SPLICES rather than passing these through: `it_classify_value`
  # replaces any concept object carrying an `ast` with that ast, so a `variable`
  # written into a pronoun expression never reaches the engine as a leaf. These
  # calls hand one over directly, which is the only way to reach that lowering
  # and the reason it is tested here rather than through `aev()`.
  amounts <- static_variable(.i$amount)

  by_object <- aev(weight = it_obj(amounts))
  by_expression <- aev(weight = it_ast(~ .i$amount))

  expect_equal(by_object$E, by_expression$E, tolerance = 1e-14)
  expect_equal(by_object$V, by_expression$V, tolerance = 1e-14)

  # An indicator is both a variable and an include; as a weight it is the value.
  flag <- indicator(.i$amount > 150)
  expect_equal(aev(weight = it_obj(flag))$E,
               aev(weight = it_ast(~ .i$amount > 150))$E, tolerance = 1e-14)
})

test_that("an indicator is accepted as an include, having no terms of its own", {
  # `include(...)` returns an indicator, which holds an `ast` where every other
  # include holds `terms`. It is the ordinary shape a population include arrives
  # in and it used to fail with "expected 'list' actual 'NULL'".
  # Amounts are 1000, 2500 and 400, so this excludes exactly one -- without which
  # the last assertion below could not tell a working gate from an absent one.
  flag <- indicator(.i$amount > 500)

  by_indicator <- aev(include = flag)
  by_gate <- aev(include = new_include(list(gate_term(it_ast(~ .i$amount > 500)))))

  expect_equal(by_indicator$E, by_gate$E, tolerance = 1e-14)
  expect_identical(by_indicator$records_included, by_gate$records_included)
  expect_lt(by_indicator$records_included, aev()$records_included)
})

test_that("overdispersion scales V and touches nothing else", {
  # THE SOLE WITNESS TO THE OVERDISPERSION MULTIPLY. Every other assertion in this suite runs at 1,
  # where the multiply is the identity, so without this a dropped scaling passes everything -- and it
  # is not a scaling anyone would notice by eye, since a V too small by a factor merely reads as a
  # tighter confidence interval.
  plain <- aev(weight = it_ast(~ .i$amount))
  dispersed <- aev(weight = it_ast(~ .i$amount), overdispersion = 2.5)

  # Overdispersion is a statement about VARIANCE alone: A and E must be bit-identical.
  expect_identical(dispersed$A, plain$A)
  expect_identical(dispersed$E, plain$E)
  expect_identical(dispersed$contributions$E, plain$contributions$E)

  expect_equal(dispersed$V, plain$V * 2.5, tolerance = 1e-14)
  expect_equal(dispersed$contributions$V, plain$contributions$V * 2.5, tolerance = 1e-14)
})

test_that("an indicator weight gives V = overdispersion x E at the output", {
  # V = E survives on the INTEGRALS, where w^2 = w folds away -- that is what the folding pass is
  # for. It is only the output that carries the overdispersion, so the identity a user meets is
  # V = Omega E. Getting this wrong in either direction is silent.
  res <- aev(overdispersion = 3)
  expect_equal(res$V, res$E * 3, tolerance = 1e-12)
  expect_equal(res$contributions$V, res$contributions$E * 3, tolerance = 1e-12)
})

test_that("a non-positive overdispersion is refused", {
  for (bad_value in c(0, -1, NaN)) {
    expect_error(aev(overdispersion = bad_value), "must be a positive number")
  }
})

test_that("the answer is one A, one E and one V, added up in the core", {
  # THE TOTALS ARE THE RESULT, and they are single numbers rather than vectors. Asserted two ways at
  # once: against the analytic answer summed in R, and against the sum of the engine's own
  # contributions. The second is what checks the accumulation rather than the arithmetic -- a chunk
  # folded twice or a partial never folded in would pass the first if the error were small enough.
  res <- aev(weight = it_ast(~ .i$amount))

  expect_length(res$A, 1L)
  expect_length(res$E, 1L)
  expect_length(res$V, 1L)

  expect_equal(res$A, sum(cols$amount * died), tolerance = 1e-12)
  expect_equal(res$E, sum(mu_value * cols$amount * exposure_years), tolerance = 1e-12)
  expect_equal(res$V, sum(mu_value * cols$amount^2 * exposure_years), tolerance = 1e-12)

  expect_equal(res$A, sum(res$contributions$A), tolerance = 1e-14)
  expect_equal(res$E, sum(res$contributions$E), tolerance = 1e-14)
  expect_equal(res$V, sum(res$contributions$V), tolerance = 1e-14)

  # Three individuals is one chunk, so there is one partial and nothing to fold.
  expect_equal(res$chunk_count, 1L)
  expect_equal(res$records_included, 3L)
})


test_that("a constant weight scales A, E and V as it should", {
  res <- aev(weight = it_ast(~ .i$amount))

  expect_equal(res$contributions$A, cols$amount * died, tolerance = 1e-12)
  expect_equal(res$contributions$E, mu_value * cols$amount * exposure_years, tolerance = 1e-12)

  # V takes the SQUARE of the weight, which is the whole reason it is computed separately.
  expect_equal(res$contributions$V, mu_value * cols$amount^2 * exposure_years, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(res$contributions$E, res$contributions$V)))
})

test_that("an indicator weight makes V equal E", {
  res <- aev(weight = it_ast(~ .i$male))

  expect_equal(res$contributions$A, as.double(cols$male & died), tolerance = 1e-12)
  expect_equal(res$contributions$E, mu_value * as.double(cols$male) * exposure_years, tolerance = 1e-12)
  expect_equal(res$contributions$V, res$contributions$E, tolerance = 1e-12)
})

test_that("V costs nothing at all for an indicator weight", {
  # `w * w` folds to `w`, which leaves V's expression identical to E's, and sharing then merges them.
  # So V is not merely cheap, it is the same operand -- one integral serves both.
  #
  # A varying mortality, so that the integrals are real rather than broadcasts of a constant.
  tbl <- mortality_table(
    x0 = 55, t0 = 2005,
    log_mu = matrix(log_mu_value, nrow = 40L, ncol = 20L)
  )
  res <- aev(mortality = it_obj(tbl), weight = it_ast(~ .i$male))

  expect_equal(sum(res$monikers == "integrate"), 1L)
  expect_equal(res$contributions$V, res$contributions$E, tolerance = 1e-12)
  expect_equal(res$contributions$E, mu_value * as.double(cols$male) * exposure_years, tolerance = 1e-12)

  # A non-indicator weight has a genuine square, so it keeps its own multiply and V differs from E.
  amounts <- aev(mortality = it_obj(tbl), weight = it_ast(~ .i$amount))
  expect_false(isTRUE(all.equal(amounts$contributions$E, amounts$contributions$V)))
})

test_that("a double column of zeros and ones is recognised, which R could not prove for itself", {
  # `flag` is a plain double, not a logical, so nothing structural says it is an indicator. The
  # column scan finds it holds only zeros and ones, and that is what makes the square foldable.
  scanned <- cols
  scanned$flag <- c(1, 0, 1)

  tbl <- mortality_table(
    x0 = 55, t0 = 2005,
    log_mu = matrix(log_mu_value, nrow = 40L, ncol = 20L)
  )
  res <- cpp_veil_aev(it_obj(tbl), it_ast(~ .i$flag), scanned, quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(sum(res$monikers == "integrate"), 1L)
  expect_equal(res$contributions$V, res$contributions$E, tolerance = 1e-12)
  expect_equal(res$contributions$E, mu_value * scanned$flag * exposure_years, tolerance = 1e-12)
})

test_that("a column of small whole numbers is not mistaken for an indicator", {
  # Integral, but a two is not zero or one -- and a half squares to a quarter, which is why an
  # interval of [0, 1] would not be enough either.
  counted <- cols
  counted$units <- c(1, 2, 1)

  tbl <- mortality_table(
    x0 = 55, t0 = 2005,
    log_mu = matrix(log_mu_value, nrow = 40L, ncol = 20L)
  )
  res <- cpp_veil_aev(it_obj(tbl), it_ast(~ .i$units), counted, quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(res$contributions$E, mu_value * counted$units * exposure_years, tolerance = 1e-12)
  expect_equal(res$contributions$V, mu_value * counted$units^2 * exposure_years, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(res$contributions$E, res$contributions$V)))
})

test_that("a comparison is an indicator wherever it came from", {
  tbl <- mortality_table(
    x0 = 55, t0 = 2005,
    log_mu = matrix(log_mu_value, nrow = 40L, ncol = 20L)
  )
  res <- aev(mortality = it_obj(tbl), weight = it_ast(~ .i$amount > 500))

  expect_equal(sum(res$monikers == "integrate"), 1L)
  expect_equal(res$contributions$V, res$contributions$E, tolerance = 1e-12)
  expect_equal(res$contributions$E, mu_value * as.double(cols$amount > 500) * exposure_years, tolerance = 1e-12)
})

test_that("a weight that varies over time is integrated, not sampled", {
  # w = age in years, so the integrand is mu * (t - b): linear in time, and the midpoint rule is
  # exact for it. A is the age at the end of the exposure, which is where a death happens.
  #
  # `.x * 1` rather than `.x`: a durationy times a number is a plain number, which is the datey
  # package's own rule and what a weight has to be. See the test below.
  res <- aev(weight = it_ast(~ .x * 1))

  from <- 2010 - birth_years
  to <- end_clicks / clicks_per_year - birth_years

  expect_equal(res$contributions$A, to * died, tolerance = 1e-9)
  expect_equal(res$contributions$E, mu_value * (to - from) * (to + from) / 2, tolerance = 1e-9)

  # V is not asserted here: its integrand is mu * (t - b) squared, and the midpoint rule is exact for
  # a linear integrand, not a quadratic one. There is no answer to write down without redoing the
  # rule, which would make the test a copy of the thing it checks.
})

test_that("a weight must be a number, as the datey rules make it", {
  # `.x` on its own is a durationy, and a duration squared is not a duration -- the datey package
  # refuses `durationy * durationy` too, so this agrees with it rather than inventing a rule.
  #
  # The message names `mul`, which is the recipe's squaring rather than anything the caller wrote.
  # Accurate but not friendly; worth improving when the recipe learns to check its weight up front.
  expect_error(aev(weight = it_ast(~ .x)), "not defined for durationy")
})

test_that("mortality is exponentiated once however many outputs read it", {
  res <- aev(weight = it_ast(~ .i$amount))

  # E and V share mu by construction -- both integrands point at the same node -- so the block holds
  # one exponential, not two. This is the property the recipe exists to get right.
  expect_equal(sum(res$monikers == "exp"), 1L)
  expect_equal(sum(res$monikers == "integrate"), 2L)
  expect_equal(sum(res$monikers == "died_value"), 1L)
})

test_that("an include clips what the AEV counts", {
  # Half the exposure of the first individual, and the clip cuts the end off, so the death goes too.
  res <- aev(include = period(2010, 2011.5))

  clipped <- pmin(end_clicks, as.integer(2011.5 * clicks_per_year)) - start_clicks
  clipped_years <- pmax(0, clipped) / clicks_per_year
  kept_death <- died & (end_clicks <= 2011.5 * clicks_per_year)

  expect_equal(res$contributions$E, mu_value * clipped_years, tolerance = 1e-12)
  expect_equal(res$contributions$A, as.double(kept_death), tolerance = 1e-12)
})

test_that("a constant weight is hoisted out of the integral, leaving one integral for E and V", {
  # This needs a mortality that VARIES, because hoisting only fires when something is left to
  # integrate. With a constant mortality both sides of `mu * w` are time-invariant and the pass
  # rightly declines -- the broadcast path already handles that better.
  tbl <- mortality_table(
    x0 = 55, t0 = 2005,
    log_mu = matrix(log_mu_value, nrow = 40L, ncol = 20L)
  )
  res <- aev(mortality = it_obj(tbl), weight = it_ast(~ .i$amount))

  # E becomes `w * integral(mu)` and V becomes `(w * w) * integral(mu)`. Those two integrals are the
  # same node once hoisting has run, so sharing merges them: ONE integral, not two, and the weight
  # never enters the time vector at all.
  expect_equal(sum(res$monikers == "integrate"), 1L)

  # And the answers are unchanged, because hoisting a constant out of an integral is exact.
  expect_equal(res$contributions$E, mu_value * cols$amount * exposure_years, tolerance = 1e-12)
  expect_equal(res$contributions$V, mu_value * cols$amount^2 * exposure_years, tolerance = 1e-12)
  expect_equal(res$contributions$A, cols$amount * died, tolerance = 1e-12)
})

test_that("a mortality table works in the recipe", {
  # Constant across the table, so mu is constant and the integral is exact. What is being tested is
  # that a table reaches the recipe at all, not the lookup, which test-veil_log_mu.R covers.
  tbl <- mortality_table(
    x0 = 55, t0 = 2005,
    log_mu = matrix(log_mu_value, nrow = 40L, ncol = 20L)
  )
  res <- aev(mortality = it_obj(tbl))

  expect_true("vector_log_mu" %in% res$monikers)
  expect_equal(res$contributions$E, mu_value * exposure_years, tolerance = 1e-12)
  expect_equal(res$contributions$A, as.double(died), tolerance = 1e-12)
})

test_that("a finer integration interval gives the same answer for a constant integrand", {
  quarterly <- aev(weight = it_ast(~ .i$amount))
  monthly <- aev(weight = it_ast(~ .i$amount), time_scale = month_scale)

  expect_equal(quarterly$E, monthly$E, tolerance = 1e-12)
  expect_equal(quarterly$V, monthly$V, tolerance = 1e-12)
})

test_that("an excluded individual is skipped rather than run, and contributes nothing", {
  res <- aev(include = period(2030, 2040))

  expect_equal(res$contributions$A, rep(0, 3), tolerance = 1e-12)
  expect_equal(res$contributions$E, rep(0, 3), tolerance = 1e-12)
  expect_equal(res$contributions$V, rep(0, 3), tolerance = 1e-12)
  expect_equal(res$A, 0)
  expect_equal(res$E, 0)
  expect_equal(res$V, 0)

  # SKIPPED, NOT MERELY ZEROED, and `records_included` is the only thing that says so. Running an
  # excluded individual over an empty grid answers zero as well, and fills no slots either, so
  # neither the values above nor the slot count below can tell the two apart. That is the point:
  # skipping is arithmetically identical and simply does less.
  expect_equal(res$records_included, 0L)
  expect_equal(res$slot_evaluations, 0L)
})

test_that("a NaN log mu gives a NaN E and V, leaving A alone", {
  # THE WRONG ANSWER THIS WAS WRITTEN FOR (found 2026-08-20): the sharing pass keyed a double
  # literal by its VALUE in a `std::map`, NaN compares unordered against everything, and so the
  # literal NaN merged with the recipe's own literal 1.0. E came back as exp(1) per unit of
  # exposure -- a perfectly ordinary-looking number for a mortality nobody could evaluate.
  #
  # Handed in as a raw literal node because the R front end refuses a written NaN at construction.
  # The path that still reaches here in normal use is a NaN in the DATA, which is legal and is
  # meant to arrive exactly like this.
  res <- aev(mortality = list(kind = "lit", value = NaN))

  expect_true(is.nan(res$E))
  expect_true(is.nan(res$V))
  expect_true(all(is.nan(res$contributions$E)))

  # A NEVER TOUCHES MU. It is the weighted count of deaths that actually happened, so a mortality
  # nobody can evaluate leaves it intact -- which is why the short circuit for a provably-NaN
  # output has to be per output rather than per block.
  expect_equal(res$A, sum(died), tolerance = 1e-12)
  expect_equal(res$contributions$A, as.double(died), tolerance = 1e-12)
})

test_that("a NaN in the data reaches E without being caught on the way", {
  # THE POLICY, STATED AS A TEST: NaN flows. No scan, no rejection, no warning -- a column may hold
  # a missing value and the result is simply NaN. This is the case the construction-time check must
  # never be extended to cover.
  dirty <- cols
  dirty$log_mu <- c(log_mu_value, NaN, log_mu_value)

  res <- aev(mortality = it_ast(~ .i$log_mu), columns = dirty)

  expect_true(is.nan(res$E))
  expect_true(is.nan(res$contributions$E[[2L]]))

  # The individuals either side are untouched: the NaN spreads through the arithmetic that uses
  # it and no further.
  expect_equal(res$contributions$E[[1L]], mu_value * exposure_years[[1L]], tolerance = 1e-12)
  expect_equal(res$contributions$E[[3L]], mu_value * exposure_years[[3L]], tolerance = 1e-12)
})
