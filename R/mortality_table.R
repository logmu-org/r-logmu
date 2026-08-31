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
#' The first age and period are determined by the `x0` and `t0` parameters
#' respectively using \eqn{\mu}{mu}-timing. (See below for what this means for
#' \eqn{q} rates.)
#'
#' The end age of the table is determined as the start age `x0`
#' plus the number of age rows as years (plus another half year for
#' \eqn{q} rates -- see below).
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
#'         \deqn{q_{xt} = 1 - \exp\left(-\int_{t-\frac{1}{2}}^{t+\frac{1}{2}}\mu_{x+\varepsilon,\, t+\varepsilon}\,\mathrm{d}\varepsilon\right)}{q_xt = 1 - exp(-integral from t-1/2 to t+1/2 of mu_(x+e, t+e) de)}
#'
#'
#'     This differs from the normal convention whereby \eqn{q_{xt}}{q_xt} relates
#'     to the year from \eqn{(x,t)} to \eqn{(x+1,t+1)}.
#'
#'     - The `x0` and `t0` parameters relate to the *middle* of the year
#'     covered by the youngest and earliest \eqn{q} rate.
#'
#'     - The end age of the resulting mortality table is the start age `x0`
#'     plus the number of age rows as years *plus an extra half year*.
#'
#' 1. The calculation of \eqn{\mu_{xt}} from \eqn{q_{xt}} includes an
#' allowance for estimated convexity determined by examining the two
#' neighbouring \eqn{q} rates
#' (by cohort for the interior of the `annual_rates` and by age at its edges).
#' This may produce artefacts if the rates are not smooth at an annual scale.
#'
#' 1. A common convention when specifying \eqn{q}-based mortality tables is to
#' include \eqn{q_\omega=1}, where \eqn{\omega} is the end age of the mortality
#' table.
#' *Do not include a \eqn{q=1} age row in the `annual_rates` argument.*
#' (If you are creating the mortality from a base table and a projection then
#' the \eqn{q_\omega=1} is likely in the base table.)
#'
#' @param x0,t0 The youngest age and earliest time respectively.
#' For \eqn{q} rates, these are the *middle* of the year covered by
#' the youngest and earliest \eqn{q} rate, i.e. \eqn{\mu}{mu}-timing.
#' @param q An age-period matrix of \eqn{q_{xt}}. It is required that
#' \eqn{0 < q < 1}.
#' @param mu An age-period matrix of \eqn{\mu_{xt}}. It is required that
#' \eqn{0 < \mu < +\infty}.
#' @param log_mu An age-period matrix of \eqn{\log\mu_{xt}}. It is required that
#' \eqn{-\infty < \log\mu < +\infty}.
#' @param name An optional name for this mortality.
#' @returns A `mortality`.
#' @export
mortality_table <- function(x0, t0, q = NULL, mu = NULL, log_mu = NULL, name = NULL) {

  x0 <- get_single_valid_durationy(x0)
  t0 <- get_single_valid_datey(t0)

  if (!is.null(name)) ensure_is_valid_name(name)

  if (!is.null(log_mu)) {

    stopifnot(
      "log_mu must be a matrix of finite values." = is_finite_matrix_of_double(log_mu),
      "Only one of q, mu or log_mu should be specified." = is.null(q) && is.null(mu)
    )
  } else if (!is.null(mu)) {

    stopifnot(
      "mu must be a matrix of doubles." = is_matrix_of_double(mu),
      "Only one of q, mu or log_mu should be specified." = is.null(q) && is.null(log_mu)
    )
    log_mu <- cpp_vec_log(mu) # 10x quicker on TJG machine
    attributes(log_mu) <- attributes(mu)
    stopifnot(
      "log(mu) must all be finite." = is_finite_matrix_of_double(log_mu)
    )

  } else if (!is.null(q)) {

    stopifnot(
      "q must be a matrix of doubles." = is_matrix_of_double(q),
      "Only one of q, mu or log_mu should be specified." = is.null(log_mu) && is.null(mu)
    )
    log_mu <- cpp_matrix_q_to_log_mu(q)
    attributes(log_mu) <- attributes(q)
    stopifnot(
      "The log mu implied by q must all be finite." = is_finite_matrix_of_double(log_mu)
    )

  } else {

    stop("At least one of q, mu or log_mu must be specified.", call. = FALSE)

  }

  end_age <- x0 + datey::durationy(nrow(log_mu));
  if (!is.null(q)) {
    # q rates cover a whole year so we need to add half a year:
    end_age <- end_age + datey::durationy(0.5);
  }

  ensure_is_single_valid_durationy(end_age)

  period_covered <- datey::durationy(ncol(log_mu));
  ensure_is_single_valid_datey(t0 + period_covered)

  structure(log_mu, class = mortality_table_class, x0 = x0, t0 = t0, end_age = end_age, name = name)
}

#' @rdname mortality
#' @export
log_mu.mortality_table <- function(x, .i, .t) {
  .b <- get_birth_datey_from_.i(.i)
  ensure_is_valid_datey(.t)
  x0 <- attr(x, "x0")
  t0 <- attr(x, "t0")

  cpp_slow_lookup_log_mu(x,
    x0_clicks = x0, t0_clicks = t0,
    b_clicks = .b, t_clicks = .t)
}

