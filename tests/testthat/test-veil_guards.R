# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests the guards that turn malformed input into a named error rather than a wrong number or a dead
# session. Every case here is one an ordinary user cannot reach through the R constructors -- an
# include is built by `period()` and friends, `indicator()` refuses `.t`, and `exp_data` validation
# will refuse missing dates. That is exactly why they are tested at this level: the engine is meant
# to be reachable from a Python or C# front end that has none of those constructors, so it has to
# hold its own preconditions rather than trust the caller to have held them.
#
# The includes below are therefore hand-built, deliberately bypassing the constructors.

clicks_per_year <- 534360L

datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}

raw_include <- function(terms) {
  structure(list(terms = terms), class = c("include", "logmu_function"))
}

cols <- list(
  birth     = datey::datey(c(1940, 1945, 1950)),
  amount    = c(1000, 2500, 400),
  male      = c(TRUE, FALSE, TRUE),
  E2R_start = datey_clicks(rep(2010, 3) * clicks_per_year),
  E2R_end   = datey_clicks(rep(2013, 3) * clicks_per_year),
  E2R_died  = c(TRUE, FALSE, TRUE)
)

flat_table <- function() {
  mortality_table(x0 = 55, t0 = 2005, log_mu = matrix(-3.2, nrow = 40L, ncol = 20L))
}

# ---- The include cannot sample time ------------------------------------------------------------

test_that("an include's gate that depends on time is refused, not run", {
  # THIS IS THE CASE THAT USED TO TAKE THE SESSION DOWN. Lowering a gate that samples time reaches
  # ensureTimeGrid from inside the include resolution that has not yet produced a grid, and recursed
  # until the stack ran out -- a crash, not a catchable error.
  #
  # The circularity is real rather than an implementation accident: the include is what decides where
  # the exposure starts and ends, so a gate needing the sample points is asking for a grid its own
  # result determines.
  timed_gate <- raw_include(list(list(kind = "gate", ast = it_ast(~ .t > datey::datey(2011)))))

  expect_error(
    cpp_veil_integrate(it_ast(~ 1), cols, quarter_scale, timed_gate),
    "gate must not depend on time"
  )
})

test_that("a missing bound in an include is refused at the crossing", {
  # NA would otherwise cross as INT_MIN and become an interval bound four millennia adrift.
  missing_bound <- raw_include(list(
    list(kind = "absolute", from = NA_integer_, to = as.integer(2012 * clicks_per_year))
  ))

  expect_error(cpp_veil_integrate(it_ast(~ 1), cols, quarter_scale, missing_bound), "missing")
})

# An offset says where a duration is measured FROM, so it has to be a date. Anything else would be
# adding years to a number and calling the result a date -- and a durationy is the case that would
# otherwise pass unnoticed, since it is clicks too and every check downstream would be satisfied.
test_that("an offset that is not a datey is refused at lowering", {
  typed <- cols
  typed$gap <- datey::durationy(c(1, 2, 3))

  offset_on_gap <- raw_include(list(
    list(kind = "offset", offset = it_ast(~ .i$gap),
         from = as.integer(0), to = as.integer(5 * clicks_per_year))
  ))

  expect_error(
    cpp_veil_integrate(it_ast(~ .i$gap * 1), typed, quarter_scale, offset_on_gap),
    "offset must be a datey"
  )
})

# The offset decides where the exposure starts, so one that sampled time would need the grid it is
# being used to build. Refused explicitly rather than left to recurse -- lowering a time-varying
# offset reaches ensureTimeGrid from inside the include resolution and runs the stack out.
test_that("an offset that depends on time is refused", {
  timed_offset <- raw_include(list(
    list(kind = "offset", offset = it_ast(~ .t),
         from = as.integer(0), to = as.integer(5 * clicks_per_year))
  ))

  expect_error(
    cpp_veil_integrate(it_ast(~ 1), cols, quarter_scale, timed_offset),
    "offset must not depend on time"
  )
})

# ---- Missing data is refused rather than silently included -------------------------------------

test_that("a missing exposure bound is refused rather than integrated", {
  # R's NA integer is the most negative int, so an NA start against a real end PASSES the ordering
  # test and describes an exposure running back four millennia. The individual would contribute an
  # enormous E rather than none, which is the failure that survives eyeballing the output.
  na_start <- cols
  na_start$E2R_start <- datey_clicks(c(2010 * clicks_per_year, NA, 2010 * clicks_per_year))

  expect_error(cpp_veil_integrate(it_ast(~ 1), na_start, quarter_scale, NULL), "missing")

  na_end <- cols
  na_end$E2R_end <- datey_clicks(c(2013 * clicks_per_year, NA, 2013 * clicks_per_year))

  expect_error(cpp_veil_integrate(it_ast(~ 1), na_end, quarter_scale, NULL), "missing")
})

test_that("a missing birth date is refused rather than read at the table edge", {
  # THE WORST OF THE MISSING-DATA CASES, because it does not look like a failure. The lattice index
  # goes hugely negative and the lookup CLAMPS to the edge of the table, so what comes back is a
  # perfectly plausible log mu for the wrong individual.
  na_birth <- cols
  na_birth$birth <- datey_clicks(c(1940 * clicks_per_year, NA, 1950 * clicks_per_year))

  expect_error(
    cpp_veil_aev(it_obj(flat_table()), NULL, na_birth, quarter_scale, NULL, no_overdispersion, 1L),
    "birth date is missing"
  )
})

test_that("a missing bound on an excluded individual does not stop the rest", {
  # The complement of the tests above, and the reason the check sits where it does. An individual an
  # include empties never reaches the grid at all, so their unusable dates are simply not consulted
  # -- exclusion is not an error, and only data that is actually going to be integrated is checked.
  na_birth <- cols
  na_birth$birth <- datey_clicks(c(1940 * clicks_per_year, NA, 1950 * clicks_per_year))
  na_birth$male <- c(TRUE, FALSE, TRUE)

  males_only <- raw_include(list(list(kind = "gate", ast = it_ast(~ .i$male))))
  res <- cpp_veil_aev(it_obj(flat_table()), NULL, na_birth, quarter_scale, males_only, no_overdispersion, 1L)

  expect_equal(res$contributions$E[2], 0)
  expect_true(all(res$contributions$E[c(1, 3)] > 0))

  # Two of the three contribute, which is what says the excluded one was skipped rather than run with
  # an empty grid. Skipping is why their missing birth date is never consulted.
  expect_equal(res$records_included, 2L)
})

# ---- Folding an indicator square follows through -----------------------------------------------

test_that("a square of a square folds the whole way down to the indicator", {
  # `(x*x)*(x*x)` is a square of `x*x`, which is itself a square of `x`, so the fold has to be
  # followed through rather than applied one step. Repointing one step leaves the parent reading
  # `x*x` -- the right VALUE, since a zero or a one is unchanged by squaring, so nothing asserted
  # about the numbers would notice.
  #
  # WHAT IT COSTS IS ONE MULTIPLY, AND THE ASSERTION HAS TO BE THAT. The single-step form was
  # measured at eight instructions against seven, with two `mul`s rather than one. Note what it does
  # NOT cost: the integral is still shared, because the sharing pass runs again after this fold and
  # merges the two structurally identical leftovers. So `integrate` count is the same either way and
  # would make a test that passes whether or not the fold followed through.
  nested <- it_ast(~ (.i$male * .i$male) * (.i$male * .i$male))
  res <- cpp_veil_aev(it_obj(flat_table()), nested, cols, quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(sum(res$monikers == "mul"), 1L)
  expect_equal(res$instruction_count, 7L)

  expect_equal(sum(res$monikers == "integrate"), 1L)
  expect_equal(res$contributions$V, res$contributions$E, tolerance = 1e-12)
  expect_equal(res$contributions$E, exp(-3.2) * as.double(cols$male) * 3, tolerance = 1e-12)
})
