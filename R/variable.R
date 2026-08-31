# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

###### Typing ######

variable_class        <- c("variable", "logmu_function")
static_variable_class <- c("static_variable", "variable", "logmu_function")

#' Variables: scalar functions of an individual and time
#'
#' @description
#' A `variable` is a scalar, real-valued function \eqn{f(i, t)} written with
#' the data pronouns `.i$field`, `.t`, `.b` and `.x`. A `static_variable` is a
#' variable that does not depend on time -- a function of the individual's
#' facts alone, \eqn{f(i)}.
#'
#' Construct one from a pronoun expression or a `~` formula. The constructors
#' return the **most specific** type they can prove: a time-invariant
#' expression becomes a `static_variable`, and one that is also structurally
#' logical becomes an `indicator`. `static_variable()` additionally asserts that
#' the expression does not use `.t` (and fails if it does); `indicator()`
#' asserts a \eqn{\{0,1\}} value. So `variable(.i$pension)` is a
#' `static_variable` and `variable(.i$pension > 0)` is an `indicator`.
#'
#' `value()` evaluates the variable against a single individual's facts `.i`
#' (a named list) and a time vector `.t`. It is the plain-R reference path for
#' testing and understanding -- not the performance path. For a
#' `static_variable`, `.t` may be omitted.
#'
#' @param expr A pronoun expression, optionally written as a `~` formula.
#' @param x A `variable`.
#' @param .i A named list of scalar facts.
#' @param .t A vector of `datey` (omit for a time-invariant variable).
#' @param ... Unused.
#' @returns
#' `variable()` and `static_variable()` return a `variable` (the latter also
#' classed `static_variable`).
#'
#' `is_variable()` and `is_static_variable()` return a scalar `logical`.
#'
#' `value()` returns a vector aligned to `.t`.
#' @examples
#' amounts <- static_variable(.i$pension)
#' value(amounts, .i = list(pension = 1000))
#'
#' v <- variable(.i$pension * 2)
#' value(v, .i = list(pension = 1000))
#' @name variable
NULL

# Wrap a parsed AST as a (static_)variable object.
new_variable <- function(ast, static) {
  structure(
    list(ast = ast),
    class = if (static) static_variable_class else variable_class
  )
}

# Return the narrowest concept type provable for `ast`:
#   uses .t                   -> variable
#   time-invariant, logical   -> indicator
#   time-invariant, otherwise -> static_variable
# (A time-invariant numeric expression that is really {0,1} cannot be proven so
# structurally; it stays a static_variable here and may still be optimised once
# the data is known.)
it_specialise <- function(ast) {
  if (it_uses_t(ast)) return(new_variable(ast, static = FALSE))
  if (it_is_logical(ast)) return(new_indicator(ast))
  new_variable(ast, static = TRUE)
}

#' @rdname variable
#' @export
variable <- function(expr) {
  it_specialise(it_capture(substitute(expr), parent.frame()))
}

#' @rdname variable
#' @export
static_variable <- function(expr) {
  ast <- it_capture(substitute(expr), parent.frame())
  if (it_uses_t(ast)) {
    stop("A `static_variable` must not depend on time `.t`; use `variable()` instead.",
         call. = FALSE)
  }
  it_specialise(ast)
}

#' @rdname variable
#' @export
is_variable <- function(x) inherits(x, "variable")

#' @rdname variable
#' @export
is_static_variable <- function(x) inherits(x, "static_variable")

ensure_is_variable <- function(x) {
  if (!is_variable(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `variable`.", call. = FALSE)
  }
}

#' @rdname variable
#' @export
value <- function(x, .i, .t = NULL, ...) UseMethod("value")

#' @rdname variable
#' @export
value.default <- function(x, .i, .t = NULL, ...) {
  arg_name <- deparse(substitute(x))
  stop("S3 function `value` is not implemented for `", arg_name, "`.", call. = FALSE)
}

#' @rdname variable
#' @export
value.variable <- function(x, .i, .t = NULL, ...) it_eval(x$ast, .i, .t)

#' @export
#' @noRd
print.variable <- function(x, ...) {
  kind <- if (is_static_variable(x)) "static_variable" else "variable"
  cat(sprintf("<%s> %s\n", kind, it_deparse(x$ast)))
  invisible(x)
}
