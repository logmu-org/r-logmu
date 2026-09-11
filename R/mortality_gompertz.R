# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' Create a Gompertz `mortality` with an age and a period slope
#'
#' @description
#' Creates a `mortality` whose \eqn{\log\mu} is linear in age and in time:
#'
#' \deqn{\log\mu_{xt} = \texttt{log\_mu\_0}
#'   + \texttt{slope\_x}\,(x - \texttt{x\_0})
#'   + \texttt{slope\_t}\,(t - \texttt{t\_0})}
#'
#' Both slopes are additive, so a period improvement is a *negative*
#' `slope_t`. That is the usual case: mortality rises with age and falls
#' with time.
#'
#' This is the natural choice of *test* mortality for the \eqn{Z} scaling of a
#' log-likelihood, where only `slope_x` materially matters, and it also serves
#' as a reference mortality or as a simple model in its own right.
#'
#' Every parameter has a default, so `gompertz_mortality()` on its own gives a
#' broadly reasonable pensioner basis:
#'
#' \deqn{\log\mu_{xt} = -3.8 + 0.1\,(x - 75) - 0.01\,(t - 2020)}
#'
#' @details
#' Internally the law is held in its cohort form, which is the shape that
#' integrates cheaply. Writing \eqn{b} for the date of birth and
#' \eqn{b_0 = \texttt{t\_0} - \texttt{x\_0}} for the cohort origin the two
#' offsets imply,
#'
#' \deqn{\log\mu_{bt} = \texttt{log\_mu\_0}
#'   - \texttt{slope\_x}\,(b - b_0)
#'   + (\texttt{slope\_x} + \texttt{slope\_t})\,(t - \texttt{t\_0})}
#'
#' which is identical to the age-period form above because
#' \eqn{x = t - b}. Birth is fixed for an individual, so the first two terms
#' are a per-individual constant that lifts out of the integral, leaving a
#' single term linear in \eqn{t} whose slope every individual shares.
#'
#' @param log_mu_0 The value of \eqn{\log\mu} at age `x_0` and time `t_0`.
#' @param slope_x The slope of \eqn{\log\mu} in age, per year. Positive for
#' mortality that rises with age.
#' @param x_0 The age at which \eqn{\log\mu} takes the value `log_mu_0`, as a
#' `durationy` or a number of years.
#' @param slope_t The slope of \eqn{\log\mu} in time, per year. Typically
#' negative, i.e. mortality improving.
#' @param t_0 The time at which \eqn{\log\mu} takes the value `log_mu_0`, as a
#' `datey` or a year.
#' @param name An optional name for this mortality.
#' @returns A `mortality`.
#' @examples
#' # The default basis.
#' log_mu(gompertz_mortality(), list(birth = datey::datey(1945)), datey::datey(2020))
#'
#' m <- gompertz_mortality(
#'   log_mu_0 = -4.6, slope_x = 0.10, x_0 = 70,
#'   slope_t = -0.02, t_0 = 2020
#' )
#' # log mu at age 70 in 2020, for someone born in 1950:
#' log_mu(m, list(birth = datey::datey(1950)), datey::datey(2020))
#' @export
gompertz_mortality <- function(
    log_mu_0 = -3.8,
    slope_x = 0.1,
    x_0 = 75,
    slope_t = -0.01,
    t_0 = 2020,
    name = NULL) {

  if (!is.null(name)) ensure_is_valid_name(name)

  stopifnot(
    "log_mu_0 must be a finite numeric scalar." = is_single_pure_finite_numeric(log_mu_0),
    "slope_x must be a finite numeric scalar." = is_single_pure_finite_numeric(slope_x),
    "slope_t must be a finite numeric scalar." = is_single_pure_finite_numeric(slope_t)
  )

  x_0 <- get_single_valid_durationy(x_0)
  t_0 <- get_single_valid_datey(t_0)

  # Fixing the age and period origins determines the cohort origin: the
  # age-period and cohort forms agree only when `b_0` is `t_0 - x_0`.
  b_0 <- t_0 - x_0
  ensure_is_single_valid_datey(b_0)

  # The one slope in time that everybody shares, once age has been resolved
  # into birth. Adding two finite doubles can only overflow to an infinity.
  cohort_slope <- slope_x + slope_t
  stopifnot(
    "slope_x + slope_t must be finite." = is_single_pure_finite_numeric(cohort_slope)
  )

  # Built through `it_capture()` so that this law takes exactly the same
  # coercion and folding path as a law a user writes out by hand.
  env <- new.env(parent = baseenv())
  env$log_mu_0 <- as.double(log_mu_0)
  env$slope_x <- as.double(slope_x)
  env$cohort_slope <- as.double(cohort_slope)
  env$b_0 <- b_0
  env$t_0 <- t_0

  ast <- it_capture(
    quote(log_mu_0 - slope_x * (.b - b_0) + cohort_slope * (.t - t_0)),
    env)

  structure(new_mortality_expr(ast), name = name)
}

# ============================================================================
# default_mortality -- the one logmu falls back on
# ============================================================================
#
# DEFINED AS `gompertz_mortality()` RATHER THAN BY REPEATING ITS NUMBERS, so the
# two cannot drift apart: there is one set of parameters and this names it. What
# a test has to pin is therefore not the agreement between them, which is
# structural, but the VALUES themselves against the specification, so that
# changing the Gompertz defaults trips something.

#' The default mortality
#'
#' @description
#' The general-purpose mortality **logmu** falls back on where one is needed and
#' none was given, and a reasonable thing to reach for whenever a broadly
#' sensible pensioner basis will do:
#'
#' \deqn{\log\mu_{xt} = -3.8 + 0.1\,(x - 75) - 0.01\,(t - 2020)}
#'
#' It sits about midway between S4PMA and S4PFA at age 75, with a shallower
#' slope to allow for higher mortality at younger ages and for plateauing at
#' higher ones.
#'
#' @details
#' Its main use is as the test mortality that sets the \eqn{Z} scale of a
#' log-likelihood, which is why [fit()] and [compare_models()] default
#' `test_mortality` to it. Only `slope_x` materially matters there, so a single
#' fixed choice is enough and a shared one is better than a good one: \eqn{Z} is
#' the yardstick that makes a difference of one in the penalised log-likelihood
#' mean one parameter's worth, so two analyses that used different test
#' mortalities are not on one scale.
#'
#' It takes no arguments deliberately. Vary it and it is no longer the default
#' mortality but a Gompertz law of your own, which is what
#' [gompertz_mortality()] is for.
#'
#' @returns A `mortality`.
#' @examples
#' default_mortality()
#'
#' # It is exactly the no-argument Gompertz.
#' identical(default_mortality(), gompertz_mortality())
#' @export
default_mortality <- function() gompertz_mortality()
