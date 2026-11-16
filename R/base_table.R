# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# NOT A `logmu_function`. A base table is a function of age alone, as at a
# date, so it is not a function of `.i` and `.t`. Inheriting the class would
# also let a base table into a pronoun expression, where `it_ast.R` splices
# any `logmu_function` as an object leaf that the engine cannot read.
base_table_class <- "base_table"

#' Create a base table of annual mortality rates by age
#'
#' @description
#' Creates a base table, i.e. annual mortality rates by age as at a given
#' date. The rates can be \eqn{q_x}, \eqn{\mu_x} or \eqn{\log\mu_x}, and must
#' be provided as a vector in one (and only one) of `q`, `mu` or `log_mu`.
#'
#' A base table is not a `mortality` by itself. Add a projection to it,
#' e.g. `base + projection` (see [q_projection_table()]), or turn it into a
#' mortality that does not vary with time using [mortality_table()].
#'
#' Everything about ages, terminal age and \eqn{q} rates is as for
#' [mortality_table()], with the as-at date `as_at` in the role of `t_start`:
#'
#' - All timing is \eqn{\mu}{mu}-timing, so for \eqn{q} rates `x_start` and
#' `as_at` are the *middle* of the year covered by the youngest \eqn{q} rate.
#'
#' - The terminal age is `x_start` plus the number of ages as years, plus
#' another half year for \eqn{q} rates. Read it back with [terminal_age()].
#'
#' - \eqn{\mu_x} is obtained from \eqn{q_x} allowing for the convexity
#' estimated from the neighbouring ages, except at the youngest and oldest
#' ages. *Do not include a \eqn{q=1} age in `q`.*
#'
#' Only the resulting \eqn{\log\mu_x} is kept. Two base tables with the same
#' \eqn{\log\mu_x} are therefore treated identically, however they were
#' created.
#'
#' @param x_start The youngest age.
#' @param as_at The as-at date of the rates, read back by [as_at()].
#' @param q A vector of \eqn{q_x}. It is required that \eqn{0 < q < 1}.
#' @param mu A vector of \eqn{\mu_x}. It is required that
#' \eqn{0 < \mu < +\infty}.
#' @param log_mu A vector of \eqn{\log\mu_x}. It is required that
#' \eqn{-\infty < \log\mu < +\infty}.
#' @param name An optional name for this base table, read back by [name()].
#' @returns A `base_table`.
#' @examples
#' base <- base_table(x_start = 60, as_at = 2020, log_mu = c(-5, -4.9, -4.8), name = "base")
#' base
#' as_at(base)
#' terminal_age(base)
#' log_mu(mortality_table(base), list(birth = datey::datey(1950)), datey::datey(2030))
#' @export
base_table <- function(x_start, as_at, q = NULL, mu = NULL, log_mu = NULL, name = NULL) {

  x_start <- get_single_valid_durationy(x_start)
  as_at <- get_single_valid_datey(as_at)

  if (!is.null(name)) ensure_is_valid_name(name)

  # The conversion is `mortality_table()`'s own, applied to one period.
  as_column <- function(rates, what) {
    if (is.null(rates)) return(NULL)
    if (!is.double(rates) || !is.null(dim(rates)) || length(rates) == 0L) {
      stop(what, " must be a vector of doubles.", call. = FALSE)
    }
    matrix(rates, ncol = 1L)
  }
  log_mu <- rates_to_log_mu(
    as_column(q, "q"),
    as_column(mu, "mu"),
    as_column(log_mu, "log_mu"),
    shape = "vector"
  )

  new_base_table(
    as.vector(log_mu),
    x_start = x_start,
    as_at = as_at,
    terminal_age = terminal_age_after(x_start, length(log_mu), from_q = !is.null(q)),
    name = name
  )
}

new_base_table <- function(log_mu, x_start, as_at, terminal_age, name) {
  structure(
    log_mu,
    class = base_table_class,
    x_start = x_start,
    as_at = as_at,
    terminal_age = terminal_age,
    name = name
  )
}

ensure_is_base_table <- function(x) {
  if (!inherits(x, "base_table")) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `base_table`.", call. = FALSE)
  }
}

#' @export
#' @noRd
print.base_table <- function(x, ...) {
  cat(sprintf("<base_table%s> %d ages from %s, as at %s\n",
              if (is.null(name(x))) "" else paste0(": ", name(x)),
              length(x), format(attr(x, "x_start", exact = TRUE)),
              format(attr(x, "as_at", exact = TRUE))))
  invisible(x)
}

#' @export
#' @noRd
`[<-.base_table` <- function(x, i, value) stop_table_immutable("base_table")

#' @export
#' @noRd
`[[<-.base_table` <- function(x, i, value) stop_table_immutable("base_table")

stop_table_immutable <- function(class_name) {
  stop("A `", class_name, "` is immutable and cannot be modified in place.",
       call. = FALSE)
}

# ---- scaling ----------------------------------------------------------------

#' Scale a base table
#'
#' @description
#' `q_scaled()` multiplies the \eqn{q_x} of a base table by a scaling factor
#' and `mu_scaled()` multiplies its \eqn{\mu_x}.
#'
#' A base table holds only \eqn{\log\mu_x}, so `q_scaled()` reads each
#' \eqn{\mu_x} as the integral of \eqn{\mu} over the year, i.e. it works with
#' \eqn{q_x = 1 - \exp(-\mu_x)}, and converts back the same way. There is no
#' convexity adjustment either way, so a scaling factor of 1 gives back the
#' original table. A scaled \eqn{q} is capped at \eqn{1 - u}, where \eqn{u}
#' is half the machine epsilon, so that a scaling that would push \eqn{q} to 1
#' or more at the oldest ages gives a very large but finite \eqn{\mu}.
#'
#' The as-at date and terminal age are unchanged. A named table gives a named
#' result, e.g. scaling `S4PFA` by 1.1 gives `110%q S4PFA` or `110%mu S4PFA`
#' (written with the Greek letter mu).
#'
#' @param base_table A [base_table()].
#' @param scaling_factor A single positive number.
#' @returns A `base_table`.
#' @examples
#' base <- base_table(x_start = 60, as_at = 2020, q = c(0.01, 0.011, 0.012), name = "base")
#' q_scaled(base, 1.1)
#' mu_scaled(base, 0.9)
#' @name scaled
NULL

#' @rdname scaled
#' @export
q_scaled <- function(base_table, scaling_factor) {
  ensure_is_base_table(base_table)
  ensure_is_scaling_factor(scaling_factor)
  q <- cap_q(scaling_factor * q_from_log_mu(as.vector(unclass(base_table))))
  log_mu <- log_mu_from_q(q)
  if (!all(is.finite(log_mu))) {
    stop("The scaled log mu must all be finite.", call. = FALSE)
  }
  rescaled_base_table(base_table, log_mu, scaling_factor, "q")
}

#' @rdname scaled
#' @export
mu_scaled <- function(base_table, scaling_factor) {
  ensure_is_base_table(base_table)
  ensure_is_scaling_factor(scaling_factor)
  log_mu <- as.vector(unclass(base_table)) + log(scaling_factor)
  rescaled_base_table(base_table, log_mu, scaling_factor, "\u03bc")
}

ensure_is_scaling_factor <- function(x) {
  if (!is_single_pure_finite_numeric(x) || x <= 0) {
    stop("`scaling_factor` must be a single positive number.", call. = FALSE)
  }
}

rescaled_base_table <- function(base_table, log_mu, scaling_factor, rate_symbol) {
  base_name <- name(base_table)
  new_base_table(
    log_mu,
    x_start = attr(base_table, "x_start", exact = TRUE),
    as_at = attr(base_table, "as_at", exact = TRUE),
    terminal_age = attr(base_table, "terminal_age", exact = TRUE),
    name = if (is.null(base_name)) NULL else {
      paste0(format_percent(scaling_factor), "%", rate_symbol, " ", base_name)
    }
  )
}

# 1.1 gives "110" and 1.125 gives "112.5": four significant figures, so the
# representation error in 1.1 does not show.
format_percent <- function(x) {
  format(100 * x, digits = 4L, scientific = FALSE, trim = TRUE)
}

# ---- q and log mu without convexity -----------------------------------------
#
# Used where a log mu that logmu already holds is read as q and back again:
# scaling and projecting. The point mu is read as the integral of mu over the
# year, both ways, so a round trip with no change gives back the log mu.

q_from_log_mu <- function(log_mu) -expm1(-exp(log_mu))

log_mu_from_q <- function(q) log(-log1p(-q))

# A q computed by scaling or projecting may reach 1 or more, typically at
# old ages in late years that may never be used. Capped at 1 - u, where u is
# the unit round-off, it maps to a log mu of about 3.6.
#
# THERE IS NO FLOOR. A positive q times a positive factor stays positive
# unless it underflows, which needs a q below about 5e-324, and a floor would
# change a very small log mu even when the factor is 1. A q that did underflow
# gives a log mu of -Inf, which the callers refuse.
q_round_off <- .Machine$double.eps / 2

cap_q <- function(q) pmin(q, 1 - q_round_off)
