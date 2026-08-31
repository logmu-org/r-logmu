# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' Create a constant `mortality`
#'
#' @description
#' Creates a constant mortality from \eqn{q}, \eqn{\mu} or \eqn{\log\mu}.
#'
#' The annual mortality rate must be provided as one
#' (and only one) of `q`, `mu` or `log_mu` as appropriate.
#'
#' @param q The probability of dying over one years. It is required that
#' \eqn{0 < q < 1}.
#' @param mu The instantaneous annual rate ('force') of mortality. It is required that
#' \eqn{0 < \mu < +\infty}.
#' @param log_mu The natural logarithm of the instantaneous annual rate ('force') of mortality. It is required that
#' \eqn{-\infty < \log\mu < +\infty}.
#' @param name An optional name for this mortality.
#' @returns A `mortality`.
#' @export
mortality_const <- function(q = NULL, mu = NULL, log_mu = NULL, name = NULL) {

  if (!is.null(name)) ensure_is_valid_name(name)

  if (!is.null(log_mu)) {

    stopifnot(
      "log_mu must be a finite numeric scalar." = is_single_pure_finite_numeric(log_mu),
      "Only one of q, mu or log_mu should be specified." = is.null(q) && is.null(mu)
    )
  } else if (!is.null(mu)) {

    stopifnot(
      "mu must be a finite numeric scalar" = is_single_pure_finite_numeric(mu),
      "Only one of q, mu or log_mu should be specified." = is.null(q) && is.null(log_mu)
    )
    log_mu <- log(mu)
    stopifnot(
      "log(mu) must be finite." = is_single_pure_finite_numeric(log_mu)
    )

  } else if (!is.null(q)) {

    stopifnot(
      "q must be a finite numeric scalar." = is_single_pure_finite_numeric(q),
      "Only one of q, mu or log_mu should be specified." = is.null(log_mu) && is.null(mu)
    )
    log_mu <- log(-log1p(-q))
    stopifnot(
      "The log mu implied by q must be finite." = is_single_pure_finite_numeric(log_mu)
    )

  } else {

    stop("At least one of q, mu or log_mu must be specified.", call. = FALSE)

  }

  structure(log_mu, class = mortality_const_class, name = name)
}

#' @rdname mortality
#' @export
log_mu.mortality_const <- function(x, .i, .t) {
  .b <- get_birth_datey_from_.i(.i)
  ensure_is_valid_datey(.t)
  rep_len(x, length(.t))
}

