# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Not `logmu_function`s, for the reason given at the top of `base_table.R`.
q_projection_table_class <- c("q_projection_table", "projection_table")
log_mu_projection_table_class <- c("log_mu_projection_table", "projection_table")

#' Create a mortality projection table
#'
#' @description
#' Creates an age-period table of log projection factors \eqn{P_{xt}}.
#' Adding it to a [base_table()] with `base + projection` gives a
#' [mortality_table()].
#'
#' Only differences in \eqn{P} between periods matter. For a base table with
#' \eqn{\log\mu^B_x}{log mu^B_x} as at \eqn{t_B}, a `log_mu_projection_table`
#' gives
#'
#' \deqn{\log\mu_{xt} = \log\mu^B_x + P_{xt} - P_{x t_B}}{log mu_xt = log mu^B_x + P_xt - P_x,t_B}
#'
#' and a `q_projection_table` gives
#'
#' \deqn{q_{xt} = q^B_x \exp(P_{xt} - P_{x t_B})}{q_xt = q^B_x exp(P_xt - P_x,t_B)}
#'
#' where \eqn{q^B_x = 1 - \exp(-\mu^B_x)}{q^B_x = 1 - exp(-mu^B_x)} and
#' \eqn{\mu_{xt} = -\log(1 - q_{xt})}{mu_xt = -log(1 - q_xt)}. There is no
#' convexity adjustment either way, so a projection with no change in \eqn{P}
#' gives back the base table. A projected \eqn{q} is capped at \eqn{1 - u},
#' where \eqn{u} is half the machine epsilon, so that a projection that would
#' push \eqn{q} to 1 or more, typically at old ages in late years, gives a very
#' large but finite \eqn{\mu}.
#'
#' \eqn{P} is interpolated between its grid points in cohort-period space,
#' exactly as [mortality_table()] interpolates \eqn{\log\mu}. Beyond its ages
#' and periods it is held at its edge values.
#'
#' The resulting mortality table has the base table's ages and periods at
#' whole years from its as-at date, extended in both directions to cover the
#' projection's ages and periods. Its terminal age is the base table's.
#' The as-at date must lie within the projection's periods. The two grids
#' need not line up.
#'
#' `as_q_projection_table()` and `as_log_mu_projection_table()` make a
#' projection from a [mortality_table()] with the same grid, setting
#' \eqn{P = \log q}{P = log q}, where \eqn{q = 1 - \exp(-\mu)}, and
#' \eqn{P = \log\mu}{P = log mu} respectively.
#'
#' @param x_start,t_start The youngest age and earliest time of the grid.
#' @param log_factor An age-period matrix of log projection factors
#' \eqn{P_{xt}}, all finite.
#' @param name An optional name for this projection, read back by [name()].
#' If both the base table and the projection are named, the resulting
#' mortality is named after both, e.g. `S4PFA + CMI_2025`.
#' @param x A [mortality_table()].
#' @returns A `projection_table`.
#' @examples
#' base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9, -4.8), name = "base")
#' improving <- log_mu_projection_table(
#'   x_start = 60,
#'   t_start = 2020,
#'   log_factor = outer(rep(1, 3), -0.02 * (0:10)),
#'   name = "two percent"
#' )
#' mortality <- base + improving
#' mortality
#' log_mu(mortality, list(birth = datey::datey(1950)), datey::datey(c(2020, 2025)))
#' @name projection_table
NULL

#' @rdname projection_table
#' @export
q_projection_table <- function(x_start, t_start, log_factor, name = NULL) {
  new_projection_table(x_start, t_start, log_factor, name, q_projection_table_class)
}

#' @rdname projection_table
#' @export
log_mu_projection_table <- function(x_start, t_start, log_factor, name = NULL) {
  new_projection_table(x_start, t_start, log_factor, name, log_mu_projection_table_class)
}

#' @rdname projection_table
#' @export
as_q_projection_table <- function(x) {
  ensure_is_mortality_table(x)
  log_factor <- log(q_from_log_mu(table_values(x)))
  if (!is_finite_matrix_of_double(log_factor)) {
    stop("The log q implied by `x` must all be finite.", call. = FALSE)
  }
  projection_like(x, log_factor, q_projection_table_class)
}

#' @rdname projection_table
#' @export
as_log_mu_projection_table <- function(x) {
  ensure_is_mortality_table(x)
  projection_like(x, table_values(x), log_mu_projection_table_class)
}

new_projection_table <- function(x_start, t_start, log_factor, name, class) {
  x_start <- get_single_valid_durationy(x_start)
  t_start <- get_single_valid_datey(t_start)

  if (!is.null(name)) ensure_is_valid_name(name)

  if (!is_finite_matrix_of_double(log_factor)) {
    stop("log_factor must be a matrix of finite values.", call. = FALSE)
  }

  period_covered <- datey::durationy(ncol(log_factor))
  ensure_is_single_valid_datey(t_start + period_covered)

  structure(
    log_factor,
    class = class,
    x_start = x_start,
    t_start = t_start,
    name = name
  )
}

ensure_is_mortality_table <- function(x) {
  if (!inherits(x, "mortality_table")) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `mortality_table`.", call. = FALSE)
  }
}

# A table's values as a bare matrix.
table_values <- function(x) {
  matrix(as.vector(unclass(x)), nrow = nrow(x), ncol = ncol(x))
}

projection_like <- function(table, log_factor, class) {
  structure(
    log_factor,
    class = class,
    x_start = attr(table, "x_start", exact = TRUE),
    t_start = attr(table, "t_start", exact = TRUE),
    name = name(table)
  )
}

#' @export
#' @noRd
print.projection_table <- function(x, ...) {
  cat(sprintf("<%s%s> %d ages from %s, %d years from %s\n",
              class(x)[[1L]],
              if (is.null(name(x))) "" else paste0(": ", name(x)),
              nrow(x), format(attr(x, "x_start", exact = TRUE)),
              ncol(x), format(attr(x, "t_start", exact = TRUE))))
  invisible(x)
}

#' @export
#' @noRd
`[<-.projection_table` <- function(x, i, value) stop_table_immutable("projection_table")

#' @export
#' @noRd
`[[<-.projection_table` <- function(x, i, value) stop_table_immutable("projection_table")

# ---- base + projection ------------------------------------------------------
#
# ONE FUNCTION FOR BOTH CLASSES. R's group generics compare the method found
# for each operand and use it only if the two are identical, so `base +
# projection` needs the same function registered as both methods.
#
# AND IT WINS AGAINST ANY OTHER METHOD. Where the other operand has an `Ops`
# method of its own (a `mortality_table`, a `datey`), R asks
# `chooseOpsMethod()`, and ours says yes, so the user gets the message below
# rather than R's `Incompatible methods` warning and plain arithmetic on the
# numbers underneath.

ops_base_projection <- function(e1, e2) {
  if (identical(.Generic, "+") && !missing(e2)) {
    if (inherits(e1, "base_table") && inherits(e2, "projection_table")) {
      return(project_base_table(e1, e2))
    }
    if (inherits(e1, "projection_table") && inherits(e2, "base_table")) {
      return(project_base_table(e2, e1))
    }
  }
  stop(
    "Operator `", .Generic, "` is not defined here. ",
    "A `base_table` and a `projection_table` can only be added to each other.",
    call. = FALSE
  )
}

#' @export
#' @noRd
Ops.base_table <- ops_base_projection

#' @export
#' @noRd
Ops.projection_table <- ops_base_projection

#' @export
#' @noRd
chooseOpsMethod.base_table <- function(x, y, mx, my, cl, reverse) TRUE

#' @export
#' @noRd
chooseOpsMethod.projection_table <- function(x, y, mx, my, cl, reverse) TRUE

# The mortality table that a base table and a projection make.
#
# THE GRID IS THE BASE TABLE'S: its ages, and whole years from its as-at date.
# It is extended in both directions to cover the projection's ages and
# periods, because beyond those both the base table and the projection are
# held at their edges, and so is the result: the table's own lookup then
# reproduces it outside its grid.
#
# P AT A GRID POINT IS LOOKED UP IN COHORT-PERIOD SPACE, by the same lookup
# that reads a mortality table, for the birth `t - x` at time `t`. Where the
# two grids line up this reads P exactly; where they do not, it interpolates
# P exactly as P would be read if it were a table of log mu.
project_base_table <- function(base, projection) {

  clicks_per_year <- as.double(unclass(datey::durationy(1)))
  clicks <- function(x) as.double(unclass(x))
  whole_years <- function(x) datey::durationy(x)

  x_base <- attr(base, "x_start", exact = TRUE)
  as_at <- attr(base, "as_at", exact = TRUE)
  log_mu_base <- as.vector(unclass(base))
  base_count <- length(log_mu_base)

  x_projection <- attr(projection, "x_start", exact = TRUE)
  t_projection <- attr(projection, "t_start", exact = TRUE)
  log_factor <- table_values(projection)
  projection_t_end <- t_projection + whole_years(ncol(log_factor) - 1L)

  if (as_at < t_projection || as_at > projection_t_end) {
    stop(
      "The base table's as-at date (", format(as_at), ") must lie within the ",
      "projection's periods (", format(t_projection), " to ",
      format(projection_t_end), ").",
      call. = FALSE
    )
  }

  # Grid indices, counted in whole years from the base table's first age and
  # its as-at date, rounded outwards.
  age_offset_start <- (clicks(x_projection) - clicks(x_base)) / clicks_per_year
  age_offset_end <- age_offset_start + nrow(log_factor) - 1
  age_index <- seq(
    min(0, floor(age_offset_start)),
    max(base_count - 1, ceiling(age_offset_end))
  )
  period_offset_start <- (clicks(t_projection) - clicks(as_at)) / clicks_per_year
  period_offset_end <- period_offset_start + ncol(log_factor) - 1
  period_index <- seq(floor(period_offset_start), ceiling(period_offset_end))

  age_count <- length(age_index)
  period_count <- length(period_index)
  ages <- x_base + whole_years(age_index)
  times <- as_at + whole_years(period_index)

  # Every grid point, ages varying fastest, as in the matrix.
  grid_ages <- rep(ages, times = period_count)
  grid_times <- rep(times, each = age_count)

  lookup <- function(at_ages, at_times) {
    cpp_slow_lookup_log_mu(
      log_factor,
      x0_clicks = x_projection, t0_clicks = t_projection,
      b_clicks = at_times - at_ages, t_clicks = at_times
    )
  }
  log_factor_change <-
    lookup(grid_ages, grid_times) -
    rep(lookup(ages, rep(as_at, age_count)), times = period_count)

  grid_log_mu_base <- rep(
    log_mu_base[pmin(pmax(age_index, 0L), base_count - 1L) + 1L],
    times = period_count
  )

  log_mu <- if (inherits(projection, "q_projection_table")) {
    q <- cap_q(q_from_log_mu(grid_log_mu_base) * exp(log_factor_change))
    log_mu_from_q(q)
  } else {
    grid_log_mu_base + log_factor_change
  }
  if (!all(is.finite(log_mu))) {
    stop("The projected log mu must all be finite.", call. = FALSE)
  }

  base_name <- name(base)
  projection_name <- name(projection)

  new_mortality_table(
    matrix(log_mu, nrow = age_count, ncol = period_count),
    x_start = ages[1L],
    t_start = times[1L],
    terminal_age = attr(base, "terminal_age", exact = TRUE),
    name = if (is.null(base_name) || is.null(projection_name)) NULL else {
      paste(base_name, "+", projection_name)
    }
  )
}
