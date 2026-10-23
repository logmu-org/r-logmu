# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Scalar and vector forms of each SIMD tier, held to the bit.
#
# WHY TO THE BIT. Each tier compiles its scalar `exp`, `log` and table lookup from the same EVE code as
# its vector kernels, so the two must agree exactly; where they do, whether a value took the scalar or
# the vector path can never show in a result. A tolerance here would hide exactly the failure this is
# for -- a scalar form that quietly became somebody else's implementation, or a multiply and add that
# the compiler fused in one form and not the other.
#
# EVERY TIER THIS CPU CAN RUN IS CHECKED IN ONE PROCESS. The test hooks bypass selection, so this does
# not depend on `LOGMU_TIER`; a tier the CPU or the build lacks comes back NULL and is skipped.

tiers <- c("baseline", "avx2", "avx512")

# Values a log mortality and a mortality can take, plus the edges of both functions' ranges. The
# length is not a multiple of any lane count, so the vector kernels' own scalar tail is exercised too,
# but nearly everything goes through full lanes, which is where a disagreement would be.
grid <- local({
  set.seed(20261012)
  c(seq(-15, 2, length.out = 2001),
    stats::runif(2000, -30, 5),
    exp(stats::runif(500, -40, 3)),
    -745.2, -708.5, -708.3, -1e-300, 0, 1e-300, 5e-324, 2.2e-308,
    709.4, 709.5, 709.8, 710, -Inf, Inf, NaN, -1)
})

test_that("each tier's scalar and vector exp and log agree to the bit", {
  checked <- 0L
  for (tier in tiers) {
    forms <- cpp_tier_exp_log_forms(tier, grid)
    if (is.null(forms)) next
    checked <- checked + 1L
    expect_identical(forms$exp_scalar, forms$exp_vector, label = paste(tier, "exp scalar"))
    expect_identical(forms$log_scalar, forms$log_vector, label = paste(tier, "log scalar"))
  }
  # Baseline is built everywhere and runs on any CPU, so at least one tier was really checked.
  expect_gte(checked, 1L)
})

test_that("the engine's scalar exp and log are the active tier's", {
  # A per-record `exp(.i$x)` runs the engine's scalar path; `vec_exp()` runs the active tier's
  # kernel. Before the two were joined, about one value in ten differed in the last place.
  cols <- list(birth = datey::datey(rep(1950, length(grid))), x = grid)
  engine_exp <- cpp_veil_eval(it_ast(~ exp(.i$x)), cols)$values
  engine_log <- suppressWarnings(cpp_veil_eval(it_ast(~ log(.i$x)), cols)$values)
  expect_identical(engine_exp, vec_exp(grid))
  expect_identical(engine_log, suppressWarnings(vec_log(grid)))
})

# ---- the table lookup -----------------------------------------------------------------------------

clicks_per_year <- 534360L

# A table whose cells are irregular, so no interpolation can come out right by accident.
lookup_table <- local({
  set.seed(7)
  matrix(-9 + 0.09 * (0:59) + stats::rnorm(60 * 12, sd = 0.05), nrow = 60, ncol = 12)
})
x0_clicks <- 50L * clicks_per_year
t0_clicks <- 2010L * clicks_per_year

# Births and instants anywhere in the table and well beyond every edge.
lookup_points <- local({
  set.seed(11)
  n <- 5000L
  list(b = as.integer(stats::runif(n, 1920, 1975) * clicks_per_year),
       t = as.integer(stats::runif(n, 2005, 2027) * clicks_per_year))
})

tier_lookup <- function(tier, b, t) {
  cpp_tier_table_log_mu(tier, lookup_table, x0_clicks, t0_clicks, b, t)
}

test_that("R's lookup is exactly the active tier's", {
  r_side <- cpp_slow_lookup_log_mu(lookup_table, x0_clicks, t0_clicks,
                                   lookup_points$b, lookup_points$t)
  expect_identical(r_side, tier_lookup(vec_active_tier(), lookup_points$b, lookup_points$t))
})

test_that("the tiers' lookups differ by no more than a rounding or two", {
  # Fused on a tier with FMA, a separate multiply and add on one without: an ulp or two apart, as
  # agreed, and never more.
  reference <- tier_lookup("baseline", lookup_points$b, lookup_points$t)
  for (tier in tiers[-1]) {
    other <- tier_lookup(tier, lookup_points$b, lookup_points$t)
    if (is.null(other)) next
    expect_equal(other, reference, tolerance = 8 * .Machine$double.eps, label = tier)
  }
})

test_that("the lookup is cohort-first bilinear interpolation, exact on the lattice", {
  # Each lattice point is a whole cohort and a whole period, so both weights are zero and the cell
  # comes back exactly. Age is period minus cohort.
  cohorts <- c(1950L, 1955L, 1958L)
  periods <- c(2010L, 2013L, 2020L)
  b <- rep(cohorts, each = length(periods)) * clicks_per_year
  t <- rep(periods, times = length(cohorts)) * clicks_per_year
  cells <- lookup_table[cbind(periods - cohorts[rep(1:3, each = 3)] - 50L + 1L,
                              periods - 2010L + 1L)]
  for (tier in tiers) {
    got <- tier_lookup(tier, as.integer(b), as.integer(t))
    if (is.null(got)) next
    expect_identical(got, cells, label = tier)
  }

  # Between the lattice points, an independent bilinear interpolation written out in R.
  reference <- function(b, t) {
    cohort <- floor((b - (t0_clicks - x0_clicks)) / clicks_per_year)
    cohort_weight <- (b - (t0_clicks - x0_clicks)) / clicks_per_year - cohort
    period <- floor((t - t0_clicks) / clicks_per_year)
    period_weight <- (t - t0_clicks) / clicks_per_year - period
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
    now + period_weight * (nxt - now)
  }
  expected <- reference(lookup_points$b, lookup_points$t)
  expect_equal(tier_lookup("baseline", lookup_points$b, lookup_points$t), expected,
               tolerance = 16 * .Machine$double.eps)
})

test_that("the lookup clamps age and period independently at every edge", {
  # Far outside on each side, at a whole cohort and period so the answer is one cell exactly.
  corner <- function(age_index, period_index) lookup_table[age_index, period_index]
  last_age <- nrow(lookup_table)
  last_period <- ncol(lookup_table)

  at <- function(birth_year, year) {
    tier_lookup("baseline", as.integer(birth_year * clicks_per_year),
                as.integer(year * clicks_per_year))
  }
  expect_identical(at(1990, 2015), corner(1, 6))            # too young: youngest age
  expect_identical(at(1900, 2015), corner(last_age, 6))     # too old: oldest age
  expect_identical(at(1950, 1990), corner(1, 1))            # before the table: first period, and
                                                            # the age it implies is below the first
  expect_identical(at(1950, 2050), corner(51, last_period))        # after the table: last period
  expect_identical(at(1900, 2050), corner(last_age, last_period))  # and above the last age too
})
