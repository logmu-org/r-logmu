# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

###### Typing ######

# An indicator is both a (time-invariant, {0,1}) `static_variable` and an
# `include`: TRUE maps to all-of-time, FALSE to an empty interval.
indicator_class <- c("indicator", "static_variable", "variable", "include", "logmu_function")

#' Indicators: a 0/1 function of an individual
#'
#' @description
#' An `indicator` is a time-invariant function of an individual's facts that
#' takes the value 0 or 1. It is the one type that is both a `variable` (a
#' \eqn{\{0,1\}} weight, for which \eqn{E = V}) and an `include` (a subset:
#' TRUE selects all of time, FALSE selects nothing).
#'
#' Build one from a pronoun expression or a `~` formula. The expression must
#' not use time `.t`. A logical result is always a valid indicator; a numeric
#' result is accepted but is only checked to be \eqn{\{0,1\}} when evaluated.
#'
#' `logical_value()` evaluates the indicator against an individual's facts `.i`
#' and returns a `logical`. Because an indicator is also an `include`,
#' `period_included()` works too, returning all-of-time or the empty interval.
#'
#' @param expr A pronoun expression, optionally written as a `~` formula.
#' @param x An `indicator`.
#' @param .i A named list of scalar facts.
#' @returns
#' `indicator()` returns an `indicator`.
#'
#' `is_indicator()` returns a scalar `logical`.
#'
#' `logical_value()` returns a `logical`.
#' @examples
#' male <- indicator(.i$sex == "male")
#' logical_value(male, .i = list(sex = "male"))
#' period_included(male, .i = list(sex = "female"))   # empty interval
#' @name indicator
NULL

# `group_name` and `name` are the two levels of naming an `includes` carries --
# an indicator is an include, so it can be an element of a breakdown and needs
# to be nameable like any other.
new_indicator <- function(ast, group_name = NULL, name = NULL) {
  structure(list(ast = ast, group_name = group_name, name = name),
            class = indicator_class)
}

# Coerce an evaluated value to a logical {0,1}, or stop.
indicator_as_logical <- function(v) {
  if (is.logical(v)) return(v)
  if (is.numeric(v) && all(v %in% c(0, 1) | is.na(v))) return(as.logical(v))
  stop("An `indicator` must take only the values 0 or 1 (or TRUE / FALSE).", call. = FALSE)
}

#' @rdname indicator
#' @export
indicator <- function(expr) {
  ast <- it_capture(substitute(expr), parent.frame())
  if (it_uses_t(ast)) {
    stop("An `indicator` must not depend on time `.t`.", call. = FALSE)
  }
  # A constant expression can be checked for 0/1-ness immediately.
  if (identical(ast$kind, "lit")) indicator_as_logical(ast$value)
  new_indicator(ast)
}

#' @rdname indicator
#' @export
is_indicator <- function(x) inherits(x, "indicator")

ensure_is_indicator <- function(x) {
  if (!is_indicator(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be an `indicator`.", call. = FALSE)
  }
}

#' @rdname indicator
#' @export
logical_value <- function(x, .i) UseMethod("logical_value")

#' @rdname indicator
#' @export
logical_value.default <- function(x, .i) {
  arg_name <- deparse(substitute(x))
  stop("S3 function `logical_value` is not implemented for `", arg_name, "`.", call. = FALSE)
}

#' @rdname indicator
#' @export
logical_value.indicator <- function(x, .i) {
  indicator_as_logical(it_eval(x$ast, .i))
}

# An indicator resolves to all-of-time when TRUE, the empty interval otherwise.
#' @rdname include
#' @export
period_included.indicator <- function(x, .i) {
  if (isTRUE(logical_value(x, .i))) include_all_interval() else include_none_interval()
}

#' @export
#' @noRd
print.indicator <- function(x, ...) {
  cat(sprintf("<indicator> %s\n", it_deparse(x$ast)))
  invisible(x)
}
