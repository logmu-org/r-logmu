# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# A/E BENCHMARK HARNESS.
#
# WHAT THIS IS FOR. Every decision about putting the SIMD kernels into the
# engine turns on where an A/E's time actually goes. The point of this file is a
# cost breakdown, not a single headline number.
#
# THE BREAKDOWN COMES FROM DIFFERENCING VARIANTS, not from instrumenting the
# engine. Each pair below differs in exactly one thing, and the gap between them
# is the cost of that thing:
#
#   table vs expression mortality  ->  the table lattice and interpolation
#   expression vs constant         ->  the exponential and the vector arithmetic
#   amounts vs indicator vs none   ->  what carrying a weight costs
#   monthly vs annual time_scale   ->  slot-proportional vs per-individual cost
#   include all vs include few     ->  the scalar spine, which is the whole cost
#                                      of an excluded individual
#   1 vs N threads                 ->  what the pool is actually delivering
#
# A CONSTANT MORTALITY IS NOT A VALID CONTROL ON ITS OWN, and finding that out
# cost a re-run. `mortality_const` does not merely drop the table lookup: it
# makes mu time-invariant, so folding evaluates the exponential once and
# hoisting lifts the whole integrand out of the integral. Against that,
# `table - const` is the cost of having a time vector at all rather than the
# cost of the table. Gompertz written as a pronoun expression is the control
# that keeps the time vector and the exponential and removes only the lattice
# arithmetic. All three variants have to stay for the differencing to mean
# anything.
#
# SLOT COUNTS ARE REPORTED because they are the kernel call length, and whether
# a kernel call beats a scalar loop turns on them. A one-year window at monthly
# steps is twelve slots; a five-year investigation is sixty, which is a
# different proposition entirely.
#
# NOT PART OF THE PACKAGE. `bench/` is in `.Rbuildignore`, so this is version
# controlled but never ships, and it may depend on things the package itself
# does not.
#
# Run it from the package root or from `bench/`:
#
#   "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" bench/bench_aev.R [records] [reps]
#
# To compare SIMD tiers, set `LOGMU_TIER` to `baseline`, `avx2` or `avx512`.
# IT CAN ONLY EVER GO DOWN -- the CPU feature checks still gate every tier, so
# asking for one the machine lacks silently leaves you on the one it has, which
# is why this script asserts that the tier it got is the tier it asked for.
#
#   LOGMU_TIER=avx2 "/c/Program Files/R/R-4.6.1/bin/Rscript.exe" bench/bench_aev.R

package_root <- local({
  for (candidate in c(".", "..")) {
    if (file.exists(file.path(candidate, "DESCRIPTION"))) {
      return(normalizePath(candidate))
    }
  }
  stop("run this from the package root or from bench/", call. = FALSE)
})

# THE INSTALLED PACKAGE, NOT `devtools::load_all()`, AND THAT IS NOT A
# PREFERENCE. `load_all()` and `pkgbuild::compile_dll()` both compile with
# `-UNDEBUG -Wall -pedantic -g -O0`; `debug = FALSE` only drops the `-g`,
# because pkgbuild REPLACES R's own CXX20FLAGS rather than adding to them and so
# never passes an `-O` flag at all. Only `R CMD INSTALL` builds the package the
# way R's `Makeconf` says, at `-O2 -DNDEBUG`, which is what CRAN ships and what
# a user runs.
#
# So install before benchmarking, from the package root:
#
#   "/c/Program Files/R/R-4.6.1/bin/R.exe" CMD INSTALL --no-multiarch --no-docs .
#
# and remember that anything going through devtools -- the RStudio Test button
# included -- leaves `-O0` objects behind, so `rm -f src/*.o src/vec_ops/*.o`
# first if one has run since.
library(logmu)

args <- commandArgs(trailingOnly = TRUE)
records <- if (length(args) >= 1) as.integer(args[[1]]) else 100000L
reps <- if (length(args) >= 2) as.integer(args[[2]]) else 5L

# THE WITNESS FOR THE TIER OVERRIDE. It is a ceiling and can only lower the
# tier, so a request for one this CPU does not have leaves the selection where
# it was -- silently, because nothing can raise an error out of C++ static
# initialisation. Checking it here is what makes the environment variable safe
# to offer at all.
requested_tier <- Sys.getenv("LOGMU_TIER", unset = "")
if (nzchar(requested_tier) && !identical(requested_tier, vec_active_tier())) {
  stop(sprintf("asked for the `%s` tier and got `%s` -- this CPU cannot offer it",
               requested_tier, vec_active_tier()), call. = FALSE)
}

# THE WITNESS FOR THE BUILD, and it stops the run rather than warning about it.
# An unoptimised binary does not read slightly pessimistic; it inverts the
# conclusion. The veil `exp` kernel measured THIRTY TIMES SLOWER than the scalar
# loop it replaces when the package came from `load_all()`, and about twice as
# FAST when the same source was installed by `R CMD INSTALL`. A warning would
# have been scrolled past.
if (!logmu:::cpp_build_optimised() || !logmu:::cpp_build_asserts_disabled()) {
  stop(sprintf(
    "the installed logmu was not built for measurement (optimised: %s, asserts disabled: %s) -- see the note at the top of this file",
    logmu:::cpp_build_optimised(), logmu:::cpp_build_asserts_disabled()),
    call. = FALSE)
}

clicks_per_year <- 534360L

# ---------------------------------------------------------------- the portfolio

# A pension scheme experience investigation, shaped the way a real one is: most
# individuals exposed for the whole window, a minority joining or leaving part
# way through, and a small number dying.
#
# THE VARIETY MATTERS MORE THAN THE REALISM. What the engine's cost depends on
# is the spread of exposure lengths, because that is the spread of slot counts,
# and a portfolio where everybody is exposed for exactly the window would hide
# the short final interval and the ragged lengths entirely.
make_portfolio <- function(n, exp_years, seed = 20260903L) {
  set.seed(seed)

  exp_start_year <- 2015
  exp_end_year <- exp_start_year + exp_years

  # Ages 50 to 95 at the start of the window, which is a pensioner population.
  birth_year <- exp_start_year - stats::runif(n, 50, 95)

  # A fifth of the scheme is exposed for part of the window only. Their start is
  # anywhere in it and their exposure at least six months, so no exposure is
  # degenerate and a good number do not land on a whole number of months.
  partial <- stats::runif(n) < 0.2
  start_offset <- ifelse(partial, stats::runif(n) * (exp_years - 0.5), 0)
  full_length <- exp_years - start_offset
  length_years <- ifelse(partial, 0.5 + stats::runif(n) * (full_length - 0.5), full_length)

  # Two per cent a year, which is about right for this age range, and a death
  # ends the exposure where it happens. Scaled by the length already drawn
  # rather than set independently, so a death always falls inside the exposure.
  died <- stats::runif(n) < 0.02 * exp_years
  length_years <- ifelse(died, length_years * pmax(stats::runif(n), 0.05), length_years)

  data <- exp_data(
    list(
      birth     = datey::datey(birth_year),
      pension   = round(stats::rlnorm(n, log(12000), 0.8)),
      male      = stats::runif(n) < 0.55,
      E2R_start = datey::datey(exp_start_year + start_offset),
      E2R_end   = datey::datey(exp_start_year + start_offset + length_years),
      E2R_died  = died
    ),
    exp_start = datey::datey(exp_start_year),
    exp_end   = datey::datey(exp_end_year)
  )

  attr(data, "bench_length_years") <- length_years
  attr(data, "bench_died") <- died
  data
}

# The slot counts the engine will actually build, worked out the way
# `buildTimeGrid` does: whole intervals, plus a short final one where the
# exposure does not divide, plus the death slot for a death.
#
# COMPUTED HERE RATHER THAN READ BACK because the engine reports only the total
# `slotEvaluations`, which is slots times instructions and cannot be turned back
# into a distribution.
slot_counts <- function(data, time_scale) {
  span <- round(attr(data, "bench_length_years") * clicks_per_year)
  dt <- logmu:::time_scale_clicks(time_scale)
  whole <- span %/% dt
  whole + ifelse(span %% dt > 0, 1L, 0L) + ifelse(attr(data, "bench_died"), 1L, 0L)
}

# ------------------------------------------------------------------- the basis

# Gompertz in age with a mild improvement trend, on an annual age-period
# lattice. Wide enough that nobody in the portfolio reaches an edge.
make_table <- function() {
  ages <- 40:110
  years <- 2010:2030
  log_mu <- outer(ages, years, function(x, t) -10.5 + 0.09 * x - 0.015 * (t - 2015))
  mortality_table(
    x0 = datey::durationy(min(ages)),
    t0 = datey::datey(min(years)),
    log_mu = log_mu,
    name = "bench_gompertz")
}

# --------------------------------------------------------------------- timing

now <- function() as.numeric(Sys.time())

# Elapsed rather than CPU time, deliberately: the pool is being measured, and
# CPU time sums across its threads.
#
# THE MINIMUM IS THE HONEST FIGURE on a machine with other things running --
# noise only ever adds. The median is reported alongside so that a wide gap
# between them is visible rather than hidden.
time_call <- function(fn, reps) {
  fn()
  times <- vapply(seq_len(reps), function(i) {
    at <- now()
    fn()
    now() - at
  }, numeric(1))
  list(min = min(times), median = stats::median(times))
}

results <- list()

record <- function(label, fn, n, reps. = reps, note = "") {
  timing <- time_call(fn, reps.)
  value <- fn()
  results[[length(results) + 1L]] <<- data.frame(
    variant = label,
    records = n,
    seconds = timing$min,
    median = timing$median,
    ns_per_record = timing$min / n * 1e9,
    E = unclass(value)[["E"]][[1]],
    note = note,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  %-34s %8.3f s   %7.0f ns/record\n",
              label, timing$min, timing$min / n * 1e9))
  invisible(NULL)
}

# ----------------------------------------------------------------------- runs

cat("logmu A/E benchmark\n")
cat(sprintf("tier            : %s, %d lanes%s\n",
            vec_active_tier(), vec_active_lanes(),
            if (nzchar(requested_tier)) " (forced)" else ""))
cat(sprintf("default threads : %d\n", logmu:::cpp_veil_default_threads()))
# THE BINARY, NOT THE SOURCE. An install can be older than the working tree and
# nothing else here would notice, so say which file was measured and when it was
# built.
cat(sprintf("binary          : built %s\n",
            format(file.info(getLoadedDLLs()[["logmu"]][["path"]])$mtime,
                   "%Y-%m-%d %H:%M:%S")))
cat(sprintf("records         : %d\n", records))
cat(sprintf("reps            : %d\n\n", reps))

table_basis <- make_table()
flat_basis <- mortality_const(log_mu = -4)

# WARMED UP ONCE BEFORE ANYTHING IS TIMED. R byte-compiles a closure on its
# first few calls, and the whole NSE front end of `aev()` is closures, so
# without this the first variant measured carries the compilation of every
# variant after it. It showed up as an eighteen-millisecond first run against
# one millisecond for the same work later.
invisible(local({
  warm <- make_portfolio(1000L, 1)
  basis <- settings(overdispersion = 2, time_scale = 1 / 12)
  for (i in 1:3) {
    aev(warm, mortality = table_basis, weight = .i$pension, settings = basis)
    aev(warm, mortality = -10.5 + 0.09 * .x, weight = .i$pension, settings = basis)
    aev(warm, mortality = flat_basis, settings = basis)
    aev(warm, mortality = table_basis, weight = .i$pension,
        include = include(.i$pension > 40000), settings = basis)
  }
}))

for (exp_years in c(1, 5)) {

  data <- make_portfolio(records, exp_years)
  note <- sprintf("exp_years=%d", exp_years)

  cat(sprintf("=== %d-year experience period ===\n", exp_years))
  for (scale in list(list(name = "monthly", value = 1 / 12),
                     list(name = "annual", value = 1))) {
    counts <- slot_counts(data, scale$value)
    cat(sprintf("slot counts (%s): min %d, median %d, mean %.1f, max %d\n",
                scale$name, min(counts), stats::median(counts),
                mean(counts), max(counts)))
  }
  cat("\n")

  monthly <- settings(overdispersion = 2, time_scale = 1 / 12)
  annual <- settings(overdispersion = 2, time_scale = 1)

  # The reference run: a real table, a real weight, monthly steps.
  record("table, amounts, monthly",
         function() aev(data, mortality = table_basis, weight = .i$pension,
                        settings = monthly),
         records, note = note)

  # Keeps the time vector and the exponential, drops the table. So
  # `table - expression` is the lattice arithmetic and the interpolation.
  record("expr gompertz, amounts, monthly",
         function() aev(data, mortality = -10.5 + 0.09 * .x, weight = .i$pension,
                        settings = monthly),
         records, note = note)

  # Time-invariant, so it folds and hoists away entirely. `expression - const`
  # is the exponential plus the vector arithmetic around it.
  record("const, amounts, monthly",
         function() aev(data, mortality = flat_basis, weight = .i$pension,
                        settings = monthly),
         records, note = note)

  # No weight means w is one, so V's integrand collapses onto E's.
  record("table, no weight, monthly",
         function() aev(data, mortality = table_basis, settings = monthly),
         records, note = note)

  # An indicator weight, where `w * w` simplifies to `w` and V collapses onto E
  # through the simplification rather than through there being no weight.
  record("table, indicator, monthly",
         function() aev(data, mortality = table_basis, weight = .i$male,
                        settings = monthly),
         records, note = note)

  # THE ONE THAT DECIDES THE SLOT AXIS. Twelve times fewer slots per individual
  # for the same individuals. If this is much faster the cost is
  # slot-proportional and lanes over slots have something to bite on; if it is
  # barely faster the cost is per individual and they do not.
  record("table, amounts, annual",
         function() aev(data, mortality = table_basis, weight = .i$pension,
                        settings = annual),
         records, note = note)

  # The scalar spine. An include that rejects most of the portfolio leaves the
  # columns, the prologue and the rejection as the entire cost of nearly
  # everybody.
  record("table, amounts, few included",
         function() aev(data, mortality = table_basis, weight = .i$pension,
                        include = include(.i$pension > 40000),
                        settings = monthly),
         records, note = note)

  # THE THREAD SWEEP. Answers do not move with thread count -- the summation
  # order is fixed by the data -- so this is purely the pool's scaling.
  for (threads in unique(c(1L, 2L, 4L, logmu:::cpp_veil_default_threads()))) {
    local({
      th <- threads
      record(sprintf("table, amounts, %d thread(s)", th),
             function() aev(data, mortality = table_basis, weight = .i$pension,
                            settings = monthly, threads = th),
             records, note = note)
    })
  }

  cat("\n")
}

# ---------------------------------------------------------------- the summary

all <- do.call(rbind, results)
cat("=== all runs ===\n")
print(all[, c("variant", "note", "seconds", "ns_per_record", "E")], row.names = FALSE)

# ALWAYS WRITTEN INTO `bench/`, never into the working directory, and that is
# not a tidiness preference. `bench/` is covered by `^bench$` in `.Rbuildignore`
# and `bench/*.csv` in `.gitignore`; the package root is covered by neither. A
# first version wrote to the working directory, a run from the repo root left a
# CSV beside `DESCRIPTION`, and `R CMD build` put it straight into the tarball.
out_file <- file.path(package_root, "bench", "bench_aev_results.csv")
utils::write.csv(all, out_file, row.names = FALSE)
cat(sprintf("\nwritten to %s\n", normalizePath(out_file)))
