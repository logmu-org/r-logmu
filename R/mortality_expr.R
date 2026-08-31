# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' Define a `mortality` from a pronoun expression
#'
#' @description
#' Builds a `mortality` from a pronoun expression for \eqn{\log\mu}. Concrete
#' `mortality` objects may appear in the expression -- each contributes its own
#' \eqn{\log\mu} at the individual and time -- so you can
#'
#' - select between mortalities with `ifelse()` (the condition must be
#'   time-invariant, as for any conditional),
#' - scale or adjust a mortality (adding in \eqn{\log\mu} is multiplying
#'   \eqn{\mu}), or
#' - write a closed-form law directly from the pronouns.
#'
#' A bare reference to a concrete `mortality` is returned unchanged.
#'
#' @param expr A pronoun expression for \eqn{\log\mu}, optionally a `~` formula.
#' @returns A `mortality`.
#' @examples
#' base <- mortality_const(log_mu = -4)
#' mortality(base + 0.05)   # scale mu up by exp(0.05)
#' @export
mortality <- function(expr) {
  ast <- it_capture(substitute(expr), parent.frame())

  # A bare reference to a concrete mortality needs no wrapping.
  if (identical(ast$kind, "obj") && is_mortality(ast$value)) return(ast$value)

  # Every concept object referenced must itself be a mortality.
  for (leaf in it_obj_leaves(ast)) {
    if (!is_mortality(leaf)) {
      stop("A `mortality` expression may only reference `mortality` objects; got a `",
           class(leaf)[[1L]], "`.", call. = FALSE)
    }
  }

  new_mortality_expr(ast)
}

new_mortality_expr <- function(ast) {
  structure(list(ast = ast), class = mortality_expr_class)
}

#' @rdname mortality
#' @export
log_mu.mortality_expr <- function(x, .i, .t) {
  it_eval(x$ast, .i, .t, obj_fn = function(o, .i, .t) log_mu(o, .i, .t))
}
