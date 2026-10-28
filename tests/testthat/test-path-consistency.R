# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# One source per operation, and the paths that compile it agreeing to within two ulps.
#
# THE POLICY. Every path that computes `exp`, `log` or a table lookup uses the same source -- EVE's
# function where EVE has one, our single definition otherwise -- and two paths may differ by up to two
# ulps, because they are compiled in different places and nothing can stop a compiler fusing a
# multiply and add in one and not the other. The same binary on the same hardware still gives the
# same answer every time.

# ---- the ulp helper itself ------------------------------------------------------------------------

test_that("the ulp helper counts ulps", {
  eps <- .Machine$double.eps
  expect_identical(ulps_apart(1, 1 + 2 * eps), 2)
  expect_identical(ulps_apart(-4, -4 - 4 * eps), 1)
  expect_identical(ulps_apart(c(Inf, NaN, 0), c(Inf, NaN, 0)), c(0, 0, 0))
  expect_identical(ulps_apart(c(Inf, NaN, 1), c(1, 1, NaN)), c(Inf, Inf, Inf))
  expect_failure(expect_within_ulps(1 + 3 * eps, 1))
  expect_success(expect_within_ulps(1 + 2 * eps, 1))
})

# ---- exp and log: the engine's scalar path against the kernel -------------------------------------

# Log mortalities and mortalities, and both functions' edges. The extremes are where EVE and the C
# library part company grossly -- EVE's `exp` overflows from about 709.44 and flushes to zero below
# about -708.4, where `std::exp` does neither -- so they are what tells the two sources apart.
grid <- local({
  set.seed(20261028)
  c(seq(-15, 2, length.out = 2001),
    stats::runif(2000, -30, 5),
    exp(stats::runif(500, -40, 3)),
    -745.2, -708.5, -708.3, -1e-300, 0, 1e-300, 5e-324, 2.2e-308,
    709.4, 709.5, 709.8, 710, -Inf, Inf, NaN, -1)
})

test_that("the engine's scalar exp and log agree with the kernels", {
  # A per-record `exp(.i$x)` runs the engine's scalar path; `vec_exp()` runs the active tier's
  # kernel over a long vector, so nearly every value goes through full lanes.
  cols <- list(birth = datey::datey(rep(1950, length(grid))), x = grid)
  engine_exp <- cpp_veil_eval(it_ast(~ exp(.i$x)), cols)$values
  engine_log <- suppressWarnings(cpp_veil_eval(it_ast(~ log(.i$x)), cols)$values)
  expect_within_ulps(engine_exp, vec_exp(grid))
  expect_within_ulps(engine_log, suppressWarnings(vec_log(grid)))
})

# ---- the table lookup ------------------------------------------------------------------------------

clicks_per_year <- 534360L

# Irregular cells, so no interpolation can come out right by accident.
lookup_table <- local({
  set.seed(7)
  matrix(-9 + 0.09 * (0:59) + stats::rnorm(60 * 12, sd = 0.05), nrow = 60, ncol = 12)
})
x0_clicks <- 50L * clicks_per_year
t0_clicks <- 2010L * clicks_per_year

lookup <- function(b, t) {
  cpp_slow_lookup_log_mu(lookup_table, x0_clicks, t0_clicks, as.integer(b), as.integer(t))
}

test_that("the lookup recycles births and instants against each other", {
  # The table has 12 periods. Until 2026-10-28 the instants were indexed modulo the period count
  # rather than their own length, so these read the wrong instants, or past the end of `t`.
  set.seed(3)
  b <- as.integer(stats::runif(40, 1945, 1970) * clicks_per_year)
  t <- as.integer(stats::runif(40, 2010, 2021) * clicks_per_year)
  one_at_a_time <- vapply(seq_along(b), function(i) lookup(b[i], t[i]), numeric(1))

  expect_identical(lookup(b, t), one_at_a_time)                         # more instants than periods
  expect_identical(lookup(b, t[1]), vapply(b, lookup, numeric(1), t = t[1]))  # one instant, many births
  expect_identical(lookup(b[1], t), vapply(t, lookup, numeric(1), b = b[1]))  # one birth, many instants
})

test_that("the lookup is exact on the lattice and bilinear between", {
  # A whole cohort and a whole period make both weights zero, so the cell comes back exactly.
  # Age is period minus cohort.
  cohorts <- rep(c(1950L, 1955L, 1958L), each = 3)
  periods <- rep(c(2010L, 2013L, 2020L), times = 3)
  cells <- lookup_table[cbind(periods - cohorts - 50L + 1L, periods - 2010L + 1L)]
  expect_identical(lookup(cohorts * clicks_per_year, periods * clicks_per_year), cells)

  # Between lattice points, against an independent cohort-first interpolation written out in R.
  set.seed(11)
  b <- as.integer(stats::runif(5000, 1920, 1975) * clicks_per_year)
  t <- as.integer(stats::runif(5000, 2005, 2027) * clicks_per_year)
  b0 <- t0_clicks - x0_clicks
  cohort <- (b - b0) %/% clicks_per_year
  cohort_weight <- ((b - b0) %% clicks_per_year) / clicks_per_year
  period <- (t - t0_clicks) %/% clicks_per_year
  period_weight <- ((t - t0_clicks) %% clicks_per_year) / clicks_per_year
  cell <- function(p, x) {
    lookup_table[cbind(pmin(pmax(x, 0), nrow(lookup_table) - 1) + 1,
                       pmin(pmax(p, 0), ncol(lookup_table) - 1) + 1)]
  }
  older_now <- cell(period, period - cohort)
  younger_now <- cell(period, period - cohort - 1)
  older_next <- cell(period + 1, period + 1 - cohort)
  younger_next <- cell(period + 1, period - cohort)
  now <- older_now + cohort_weight * (younger_now - older_now)
  nxt <- older_next + cohort_weight * (younger_next - older_next)
  expect_within_ulps(lookup(b, t), now + period_weight * (nxt - now))
})

test_that("the lookup clamps age and period independently at every edge", {
  # Far outside on each side, at a whole cohort and period so the answer is one cell exactly.
  at <- function(birth_year, year) lookup(birth_year * clicks_per_year, year * clicks_per_year)
  last_age <- nrow(lookup_table)
  last_period <- ncol(lookup_table)

  expect_identical(at(1990, 2015), lookup_table[1, 6])                    # too young
  expect_identical(at(1900, 2015), lookup_table[last_age, 6])             # too old
  expect_identical(at(1950, 1990), lookup_table[1, 1])                    # before the table, so young
  expect_identical(at(1950, 2050), lookup_table[51, last_period])         # after the table
  expect_identical(at(1900, 2050), lookup_table[last_age, last_period])   # after it, and too old
})
