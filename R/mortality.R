# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

###### Typing ######

mortality_class <- c("mortality", "logmu_function")
mortality_const_class <- c("mortality_const", "mortality", "logmu_function")
mortality_table_class <- c("mortality_table", "mortality", "logmu_function")
mortality_expr_class <- c("mortality_expr", "mortality", "logmu_function")
mortality_prophazard_class <- c("mortality_prophazard", "mortality", "logmu_function")


#' Mortality-related functionality
#'
#' @description
#'
#' A **logmu** mortality is an object class that defines the annual mortality
#' rate \eqn{\mu_{it}} for individual \eqn{i} at time \eqn{t}, where
#' \eqn{i} contains all the available invariant information relating to the
#' relevant individual.
#'
#' All **logmu** actually works in terms of \eqn{\log\mu_{it}} for reasons set out
#' TBC.
#'
#' Access to \eqn{\log\mu_{it}} for a mortality object is implicit using
#' `.i` and `.t` pronouns.
#'
#' For testing and understanding, you can access the \eqn{\log\mu_{it}}
#' calculation by defining your own `.i` and `.t` variables and then calling
#' `log_mu()` on a mortality. Do *not* use this for performant scenarios.
#'
#' These are the core types of mortality object:
#'
#' |Type       |**datey** function |Notes|
#' |:----------|:------------------|:---
#' |Constant   |[mortality_const()]|
#' |Age-period |[mortality_table()]|Mortality rates that are smooth at an annual scale|
#' |Gompertz   |[gompertz_mortality()]|Log mu linear in age and in time|
#' |Model      |TBC|
#'
#' You can build on these:
#'
#' - Add a projection to a [base_table()] to make an age-period table, see
#' [projection_table]. Scale a base table with [q_scaled()] or [mu_scaled()].
#' - Select between tables using TBC.
#' - Scale mortality rates using TBC.
#' - Provide sub-annual mortality variation using TBC.
#'
#' LINK TO OTHER PACKAGES TO GET STANDARD TABLES
#'
#' Finally there are a couple of
#' mortality-related helper functions:
#'
#' - `is_mortality(x)` tests whether `x` is a `mortality`.
#'
#' - `terminal_age(x)` gets the terminal age of a [mortality_table()] or a
#' [base_table()], i.e. the age by which everyone is assumed to have died.
#' Any other `mortality` may have a terminal age that depends on the
#' individual, so use `terminal()`.
#'
#' - `terminal(x, .i)` gets the date by which each individual in `.i` is
#' assumed to have died. A `mortality` with no terminal age, such as a constant
#' or a Gompertz, gives `datey::datey(datey::valid_years_end)`. An expression
#' that refers to a [mortality_table()] is not yet supported.
#'
#' Terminal ages and dates are provided for valuation calculations and are
#' currently ignored for A/E calculations.
#'
#' @param x The `mortality` object.
#' @param .i A list of named scalar arguments.
#' @param .t A vector of `datey` representing time.
#' @returns
#' `is_mortality()` returns a scalar `logical`.
#'
#' `terminal_age()` returns a scalar `durationy`.
#'
#' `terminal()` returns a `datey` for each individual in `.i`.
#'
#' `log_mu()` returns a vector of \eqn{\log\mu}{log mu} values at `.t`.
#' @examples
#' mortality <- mortality_table(x_start = 70, t_start = 2020, q = matrix(0.01))
#' .t <- datey::datey(2020)
#' .i <- list(birth = datey::datey(1950))
#' log_mu(mortality, .i, .t) # log(-log(1 - 0.01)) = -4.600149
#' @name mortality
NULL

#' @rdname mortality
#' @export
is_mortality <- function(x) inherits(x, "mortality")

ensure_is_mortality <- function(x) {
  if (!is_mortality(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `mortality`.", call. = FALSE)
  }
}

#' @rdname mortality
#' @export
terminal_age <- function(x) {
  if (!inherits(x, c("mortality_table", "base_table"))) {
    stop(
      "`terminal_age()` needs a `mortality_table` or a `base_table`. ",
      "For a general `mortality`, use `terminal(x, .i)`.",
      call. = FALSE
    )
  }
  attr(x, "terminal_age", exact = TRUE)
}

#' @rdname mortality
#' @export
terminal <- function(x, .i) UseMethod("terminal")

#' @rdname mortality
#' @export
terminal.default <- function(x, .i) {
  arg_name <- deparse(substitute(x))
  stop("S3 function `terminal` is not implemented for `", arg_name, "`.", call. = FALSE)
}

# What `terminal()` gives for an individual under a mortality with no terminal
# age: the last date a `datey` can hold.
no_terminal <- function(.i) {
  .b <- get_birth_dateys_from_.i(.i)
  rep(datey::datey(datey::valid_years_end), length(.b))
}

#' @rdname mortality
#' @export
log_mu <- function(x, .i, .t) UseMethod("log_mu")

#' @rdname mortality
#' @export
log_mu.default <- function(x, .i, .t) {
  arg_name <- deparse(substitute(x))
  stop("S3 function `log_mu` is not implemented for `", arg_name, "`.", call. = FALSE)
}

# The list-like modification ops and immutability guard are inherited from the
# `logmu_function` base class (see logmu_function.R).
