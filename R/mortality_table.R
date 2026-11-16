# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' Create an `mortality` from a 2D age-period matrix of annual mortality rates
#'
#' @description
#' Creates an age-period mortality table from a 2D age-period matrix of
#' annual mortality rates that can be \eqn{q_{xt}}, \eqn{\mu_{xt}} or
#' \eqn{\log\mu_{xt}}.
#'
#' The first age and period are determined by the `x_start` and `t_start` parameters
#' respectively using \eqn{\mu}{mu}-timing. (See below for what this means for
#' \eqn{q} rates.)
#'
#' The terminal age of the table, the age by which everyone is assumed to have
#' died, is the start age `x_start` plus the number of age rows as years (plus
#' another half year for \eqn{q} rates -- see below). Read it back with
#' [terminal_age()].
#'
#' The annual mortality rates must be provided as a 2D age-period matrix in one
#' (and only one) of `q`, `mu` or `log_mu` as appropriate.
#'
#' It is an implicit assumption that these rates are 'smooth' at the annual
#' scale, i.e. the second differences of \eqn{\mu} and \eqn{\log\mu} are
#' 'small', e.g.
#' \eqn{\Delta^2\mu<10\\%\mu} and
#' \eqn{\Delta^2\log\mu<10\\%}.
#' If you want to allow for realistic, i.e. non-smooth, historical
#' annual and sub-annual noise then use a `variation`.
#'
#' Important notes for \eqn{q} rates:
#'
#' 1. All timing is \eqn{\mu}{mu}-timing. This means:
#'
#'     - The definition of \eqn{q_{xt}}{q_xt} is *centred* on \eqn{(x,t)}, i.e.
#'
#'         \deqn{q_{xt} = 1 - \exp\left(-\int_{-\frac{1}{2}}^{+\frac{1}{2}}\mu_{x+\varepsilon,\, t+\varepsilon}\,\mathrm{d}\varepsilon\right)}{q_xt = 1 - exp(-integral from -1/2 to +1/2 of mu_(x+e, t+e) de)}
#'
#'
#'     This differs from the normal convention whereby \eqn{q_{xt}}{q_xt} relates
#'     to the year from \eqn{(x,t)} to \eqn{(x+1,t+1)}.
#'
#'     - The `x_start` and `t_start` parameters relate to the *middle* of the year
#'     covered by the youngest and earliest \eqn{q} rate.
#'
#'     - The terminal age of the resulting mortality table is the start age `x_start`
#'     plus the number of age rows as years *plus an extra half year*.
#'
#' 1. The calculation of \eqn{\mu_{xt}} from \eqn{q_{xt}} includes an
#' allowance for estimated convexity determined by examining the two
#' neighbouring \eqn{q} rates
#' (by cohort in the interior of `q`,
#' by age in its first and last periods, and
#' not at all for its youngest and oldest ages,
#' i.e. simply \eqn{\mu = -\log(1 - q)}{mu = -log(1 - q)}).
#' This may produce artefacts if the rates are not smooth at an annual scale.
#'
#' 1. A common convention when specifying \eqn{q}-based mortality tables is to
#' include \eqn{q_\omega=1}, where \eqn{\omega} is the terminal age of the mortality
#' table.
#' *Do not include a \eqn{q=1} age row in the `q` argument.*
#' (If you are creating the mortality from a base table and a projection then
#' the \eqn{q_\omega=1} is likely in the base table.)
#'
#' A [base_table()] can be turned into a mortality by passing it as the only
#' argument, optionally with a `name`. The result is a table of one period at
#' the base table's as-at date, so its log mu does not vary with time.
#'
#' @param x_start,t_start The youngest age and earliest time respectively.
#' For \eqn{q} rates, these are the *middle* of the year covered by
#' the youngest and earliest \eqn{q} rate, i.e. \eqn{\mu}{mu}-timing.
#' Alternatively `x_start` can be a [base_table()], in which case `t_start`,
#' `q`, `mu` and `log_mu` must not be given.
#' @param q An age-period matrix of \eqn{q_{xt}}. It is required that
#' \eqn{0 < q < 1}.
#' @param mu An age-period matrix of \eqn{\mu_{xt}}. It is required that
#' \eqn{0 < \mu < +\infty}.
#' @param log_mu An age-period matrix of \eqn{\log\mu_{xt}}. It is required that
#' \eqn{-\infty < \log\mu < +\infty}.
#' @param name An optional name for this mortality, read back by [name()] and
#' printed wherever the mortality appears. From a [base_table()], the default
#' is the base table's own name.
#' @returns A `mortality`.
#' @export
mortality_table <- function(x_start, t_start, q = NULL, mu = NULL, log_mu = NULL, name = NULL) {

  if (inherits(x_start, "base_table")) {
    return(mortality_table_from_base(x_start, t_start, q, mu, log_mu, name))
  }

  x_start <- get_single_valid_durationy(x_start)
  t_start <- get_single_valid_datey(t_start)

  if (!is.null(name)) ensure_is_valid_name(name)

  log_mu <- rates_to_log_mu(q, mu, log_mu, shape = "matrix")

  new_mortality_table(
    log_mu,
    x_start = x_start,
    t_start = t_start,
    terminal_age = terminal_age_after(x_start, nrow(log_mu), from_q = !is.null(q)),
    name = name
  )
}

# A `base_table` given as `mortality_table()`'s first argument becomes a table
# of one period at its as-at date. The lookup holds a table at its edges, so
# the result is the base table's log mu at every time.
mortality_table_from_base <- function(base, t_start, q, mu, log_mu, name) {
  if (!missing(t_start) || !is.null(q) || !is.null(mu) || !is.null(log_mu)) {
    stop(
      "Give `mortality_table()` either a `base_table` alone or ",
      "`x_start`, `t_start` and one of `q`, `mu` or `log_mu`.",
      call. = FALSE
    )
  }
  if (is.null(name)) {
    name <- name(base)
  } else {
    ensure_is_valid_name(name)
  }
  new_mortality_table(
    matrix(as.vector(unclass(base)), ncol = 1L),
    x_start = attr(base, "x_start", exact = TRUE),
    t_start = attr(base, "as_at", exact = TRUE),
    terminal_age = attr(base, "terminal_age", exact = TRUE),
    name = name
  )
}

# Converts whichever one of `q`, `mu` or `log_mu` was supplied to log mu, as a
# matrix with the input's own attributes. A `q` is converted with the
# convexity adjustment in `src/mortality/convert_m_to_mu.hpp`. `shape` names
# what the caller was given, for the messages: `base_table()` hands its vector
# over as a one-column matrix.
rates_to_log_mu <- function(q, mu, log_mu, shape) {

  if (!is.null(log_mu)) {

    if (!is.null(q) || !is.null(mu)) {
      stop("Only one of q, mu or log_mu should be specified.", call. = FALSE)
    }
    if (!is_finite_matrix_of_double(log_mu)) {
      stop("log_mu must be a ", shape, " of finite values.", call. = FALSE)
    }

  } else if (!is.null(mu)) {

    if (!is.null(q)) {
      stop("Only one of q, mu or log_mu should be specified.", call. = FALSE)
    }
    if (!is_matrix_of_double(mu)) {
      stop("mu must be a ", shape, " of doubles.", call. = FALSE)
    }
    log_mu <- cpp_vec_log(mu) # 10x quicker on TJG machine
    attributes(log_mu) <- attributes(mu)
    if (!is_finite_matrix_of_double(log_mu)) {
      stop("log(mu) must all be finite.", call. = FALSE)
    }

  } else if (!is.null(q)) {

    if (!is_matrix_of_double(q)) {
      stop("q must be a ", shape, " of doubles.", call. = FALSE)
    }
    log_mu <- cpp_matrix_q_to_log_mu(q)
    attributes(log_mu) <- attributes(q)
    if (!is_finite_matrix_of_double(log_mu)) {
      stop("The log mu implied by q must all be finite.", call. = FALSE)
    }

  } else {

    stop("At least one of q, mu or log_mu must be specified.", call. = FALSE)

  }

  log_mu
}

# The row after the last is taken to be the one in which everyone dies. For
# mu that row is a single age; a q row covers a whole year centred on its age,
# so everyone has died half a year later.
terminal_age_after <- function(x_start, age_count, from_q) {
  terminal_age <- x_start + datey::durationy(age_count)
  if (from_q) {
    terminal_age <- terminal_age + datey::durationy(0.5)
  }
  ensure_is_single_valid_durationy(terminal_age)
  terminal_age
}

# Every `mortality_table` is made here, whether from rates, from a
# `base_table` or from a base table plus a projection.
new_mortality_table <- function(log_mu, x_start, t_start, terminal_age, name) {
  period_covered <- datey::durationy(ncol(log_mu))
  ensure_is_single_valid_datey(t_start + period_covered)

  structure(
    log_mu,
    class = mortality_table_class,
    x_start = x_start,
    t_start = t_start,
    terminal_age = terminal_age,
    name = name
  )
}

#' @export
#' @noRd
print.mortality_table <- function(x, ...) {
  cat(sprintf("<mortality_table%s> %d ages from %s, %d years from %s\n",
              if (is.null(name(x))) "" else paste0(": ", name(x)),
              nrow(x), format(attr(x, "x_start", exact = TRUE)),
              ncol(x), format(attr(x, "t_start", exact = TRUE))))
  invisible(x)
}

#' @rdname mortality
#' @export
log_mu.mortality_table <- function(x, .i, .t) {
  .b <- get_birth_dateys_from_.i(.i)
  ensure_is_valid_datey(.t)
  x_start <- attr(x, "x_start")
  t_start <- attr(x, "t_start")

  cpp_slow_lookup_log_mu(x,
    x0_clicks = x_start, t0_clicks = t_start,
    b_clicks = .b, t_clicks = .t)
}

#' @rdname mortality
#' @export
terminal.mortality_table <- function(x, .i) {
  .b <- get_birth_dateys_from_.i(.i)
  .b + attr(x, "terminal_age", exact = TRUE)
}

