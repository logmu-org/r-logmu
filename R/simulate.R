# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# simulate -- experience data from a known mortality
# ============================================================================
#
# The one thing a real dataset can never offer is the truth behind it. This
# manufactures a dataset whose truth is known -- the mortality it was drawn
# from -- so that a fit can be checked against it, an A/E should sit at one,
# and a wrong answer is a wrong answer rather than a difference of opinion.

#' Simulate experience data from a mortality
#'
#' @description
#' Draws a death time for every life in a population under a known
#' `mortality`, and returns the result as an [exp_data()] with the same shape as
#' real data: `E2R_start`, `E2R_end` and `E2R_died`, and every other column the
#' population came with. A fit to the result can then be judged against the
#' mortality that made it.
#'
#' @section The population:
#' `lives` is a data frame (or list of columns) with a `birth` column of
#' `datey`s, whatever covariate columns the mortality and any later analysis
#' need, and two optional `datey` columns:
#'
#' * `E2R_start` -- when each life's exposure begins: retirement, joining, the
#'   first date the record is reliable. Absent, every life starts at
#'   `exp_start`. A life is alive at its start by construction, which is the
#'   left-truncation the E2R framework assumes.
#' * `E2R_cap` -- the date a life would leave observation for a reason other
#'   than death: withdrawal, transfer out, a record closed. `NA` means no such
#'   date. Absent, no life has one.
#'
#' Each life is observed from `max(E2R_start, exp_start)` to the earliest of
#' its death, its `E2R_cap` and `exp_end`. A death before the cap is a death and
#' the cap plays no part; a life still alive at the cap is censored there. Lives
#' with no observation at all are dropped, with a message saying how many.
#' `E2R_cap` is dropped from the result, so that simulated data has exactly
#' the columns real data has. Names beginning `E2R_` are reserved and cannot
#' be read by a pronoun expression, so it could not have leaked into a
#' calculation either way.
#'
#' @section How deaths are drawn:
#' Each life gets a target \eqn{-\log U} for uniform \eqn{U}, and dies where
#' its cumulative hazard \eqn{\int\mu\,dt} from its start first reaches it.
#' The hazard is accumulated on the same grid the engine integrates on: from
#' the life's own start in steps of `time_scale`, \eqn{\mu} read at each step's
#' midpoint, a shorter final step to the cap, and the death placed within the
#' crossing step at the constant hazard of that step. So a fit to simulated
#' data is unbiased against the mortality that made it, up to sampling noise.
#'
#' A death is recorded to the day it fell on, at [datey::mid_day()], exactly
#' as a real death is. Deaths are independent across lives with no extra
#' variation, so the data is Poisson-like: an overdispersion of 1 is the truth
#' for it.
#'
#' @section Reproducibility:
#' The targets are drawn with R's own generator, so `set.seed()` before the
#' call, or the `seed` argument, makes the result repeatable.
#'
#' @param object The `mortality` to simulate from: any mortality, including a
#'   fitted model's `$mortality` and a [mortality()] expression switching on
#'   the population's columns.
#' @param nsim How many datasets to draw. One returns an `exp_data`; more
#'   returns a list of them.
#' @param seed Passed to `set.seed()` first when not `NULL`.
#' @param lives The population. See *The population*.
#' @param exp_start,exp_end The experience period, as `datey`s or years.
#' @param time_scale The step of the hazard grid, as in [settings()].
#' @param ... Ignored.
#' @returns An [exp_data()], or a list of `nsim` of them.
#' @examples
#' set.seed(1)
#' lives <- data.frame(
#'   birth   = datey::datey(2015 - sample(60:85, 500, replace = TRUE)),
#'   sex     = sample(c("male", "female"), 500, replace = TRUE),
#'   pension = round(rlnorm(500, log(12000), 0.7))
#' )
#'
#' male <- gompertz_mortality(log_mu_0 = -3.6)
#' female <- gompertz_mortality(log_mu_0 = -4.0)
#' truth <- mortality(if (.i$sex == "male") male else female)
#'
#' data <- simulate(truth, lives = lives,
#'                  exp_start = datey::datey(2015), exp_end = datey::datey(2020))
#' data
#'
#' # The data is neutral against the mortality that made it, within noise.
#' aev(data, mortality = truth, overdispersion = 1)
#' @name simulate
NULL

#' @rdname simulate
#' @importFrom stats simulate
#' @export
simulate.mortality <- function(object,
                               nsim = 1,
                               seed = NULL,
                               lives,
                               exp_start,
                               exp_end,
                               time_scale = default_time_scale,
                               ...) {

  ensure_is_mortality(object)
  if (missing(lives)) stop("`lives` is required.", call. = FALSE)
  if (missing(exp_start) || missing(exp_end)) {
    stop("`exp_start` and `exp_end` are required.", call. = FALSE)
  }
  if (!is_single_pure_finite_numeric(nsim) || nsim < 1 || nsim != round(nsim)) {
    stop("`nsim` must be a single whole number of at least 1.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  population <- simulate_population(lives, exp_start, exp_end)
  clicks <- time_scale_clicks(time_scale)

  draws <- lapply(seq_len(nsim), function(draw) {
    simulate_one(object, population, clicks)
  })
  if (nsim == 1) draws[[1L]] else draws
}

# The population made ready: columns as a list, `birth` checked, each life's
# observation window `from` and `to` as `datey`s, and lives with no window
# dropped. Everything about the input is settled here, once, so `nsim` draws
# share it.
simulate_population <- function(lives, exp_start, exp_end) {
  if (is_exp_data(lives)) {
    stop("`lives` is already an `exp_data`. Give the population without its E2Rs.",
         call. = FALSE)
  }
  if (is.data.frame(lives)) lives <- as.list(lives)
  if (!is.list(lives) || is.null(names(lives)) || any(!nzchar(names(lives)))) {
    stop("`lives` must be a data frame or a named list of columns.", call. = FALSE)
  }
  for (name in c("E2R_end", "E2R_died")) {
    if (name %in% names(lives)) {
      stop("`lives` must not carry `", name, "`: that is what is simulated.", call. = FALSE)
    }
  }

  birth <- lives[["birth"]]
  if (is.null(birth) || !is_valid_datey(birth) || length(birth) == 0L) {
    stop("`lives` must have a `birth` column of valid `datey`s.", call. = FALSE)
  }
  n <- length(birth)
  lengths <- vapply(lives, length, integer(1L))
  if (any(lengths != n)) {
    stop("Every column of `lives` must have the same length as `birth`.", call. = FALSE)
  }

  exp_start <- get_single_valid_datey(exp_start)
  exp_end <- get_single_valid_datey(exp_end)
  if (!(exp_start < exp_end)) stop("`exp_start` must be before `exp_end`.", call. = FALSE)

  from <- simulate_bound(lives[["E2R_start"]], "E2R_start", n, exp_start, later = TRUE)
  to <- simulate_bound(lives[["E2R_cap"]], "E2R_cap", n, exp_end, later = FALSE)

  observed <- from < to
  if (!all(observed)) {
    message(sum(!observed), " of ", n, " lives had no exposure within the period and ",
            "were dropped.")
  }

  columns <- lives[setdiff(names(lives), c("E2R_start", "E2R_cap"))]
  columns <- lapply(columns, function(column) column[observed])

  list(columns = columns, from = from[observed], to = to[observed],
       exp_start = exp_start, exp_end = exp_end)
}

# A per-life bound clipped to the period, as a `datey`: the later of a start
# and `exp_start`, or the earlier of a cap and `exp_end`. An absent column, or
# an NA in it, is the period's own bound. `rep()`, `pmax()` and `pmin()` on a
# `datey` need datey 0.1.2, which DESCRIPTION pins -- 0.1.1 had no `rep`
# method, and base `pmax()` recycles with `rep()`, so it handed back the limit
# for every life.
simulate_bound <- function(column, name, n, limit, later) {
  if (is.null(column)) return(rep(limit, n))
  if (!datey::is_datey(column)) {
    stop("`lives$", name, "` must be a `datey` column.", call. = FALSE)
  }
  column[is.na(column)] <- limit
  if (later) pmax(column, limit) else pmin(column, limit)
}

# A `datey` from double clicks. `datey()` takes years and rounds to the click,
# so an integer-valued click count survives the round trip exactly.
simulate_datey <- function(clicks) {
  datey::datey(clicks / as.double(unclass(datey::durationy(1))))
}

# One dataset: the walk over the clicks and the E2R triple it leaves behind.
#
# EVERYTHING IS IN DOUBLE CLICKS AND MIRRORS `TimeGrid.hpp`: from each life's
# own `from`, full steps of `clicks` sampled at their exact midpoints, then
# whatever is left to `to` sampled at its midpoint with banker's halving where
# the sum is odd. Clicks rather than dateys because the midpoint rule is
# integer arithmetic on clicks, two click counts of a billion each overflow an
# R integer when added, and a double carries every integer to 2^53 exactly.
# The walk is vectorised over lives and loops over steps, so a population of
# any size costs one `log_mu()` call per step of the longest window.
simulate_one <- function(object, population, clicks) {
  from <- as.double(unclass(population$from))
  to <- as.double(unclass(population$to))
  n <- length(from)
  clicks <- as.double(clicks)
  clicks_per_year <- as.double(unclass(datey::durationy(1)))

  target <- -log(stats::runif(n))
  hazard <- numeric(n)
  death <- rep(NA_real_, n)

  step <- 0L
  repeat {
    begin <- from + step * clicks
    active <- which(is.na(death) & begin < to)
    if (length(active) == 0L) break

    begin <- begin[active]
    end <- pmin(begin + clicks, to[active])
    width <- end - begin
    mid <- simulate_midpoint(begin, end)

    .i <- lapply(population$columns, function(column) column[active])
    mu <- exp(log_mu(object, .i, simulate_datey(mid)))
    gained <- mu * width / clicks_per_year

    crossed <- hazard[active] + gained >= target[active]
    fraction <- (target[active] - hazard[active]) / gained
    death[active[crossed]] <- begin[crossed] + fraction[crossed] * width[crossed]
    hazard[active] <- hazard[active] + gained

    step <- step + 1L
  }

  died <- !is.na(death)

  # TO THE DAY, AT MIDDAY, as a real death is recorded. `mid_day()` on a datey
  # is the middle of the day the instant fell in, which can be later than the
  # instant by up to half a day and so is capped at the window's end -- or
  # earlier than it, and a life that entered in the afternoon and died that
  # evening would then have an empty E2R, so such a death keeps its instant.
  instant <- round(ifelse(died, death, from))
  midday <- as.double(unclass(datey::mid_day(simulate_datey(instant))))
  midday[midday <= from] <- instant[midday <= from]
  end <- ifelse(died, pmin(midday, to), to)

  columns <- population$columns
  columns$E2R_start <- population$from
  columns$E2R_end <- simulate_datey(end)
  columns$E2R_died <- died

  exp_data(columns, exp_start = population$exp_start, exp_end = population$exp_end)
}

# The midpoint of a step in clicks. A full step has an exact midpoint because
# half of every permitted step is a whole number of clicks; only the final,
# short step can have an odd sum, and it is halved to the nearer even click so
# the residual half-click errors cancel across a portfolio rather than
# accumulate -- the rule `TimeGrid.hpp` follows.
simulate_midpoint <- function(begin, end) {
  total <- begin + end
  half <- total %/% 2
  odd <- total %% 2 == 1
  half[odd & half %% 2 == 1] <- half[odd & half %% 2 == 1] + 1
  half
}
