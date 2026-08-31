# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' Parse a pronoun expression into an AST
#'
#' @description
#' Parses an R expression written with the data pronouns `.i$field`, `.t`,
#' `.b` and `.x` into a small abstract syntax tree. Either a bare expression
#' or a `~` formula may be supplied; constant sub-expressions are evaluated in
#' the calling environment (or the formula's environment) and folded.
#'
#' @param expression A pronoun expression, optionally written as a `~` formula.
#' @returns An `it_node`, the root of the parsed expression tree.
#' @examples
#' pronoun_expressions(.i$pension > 0)
#' pronoun_expressions(~ .t - .b)
#' @export
pronoun_expressions <- function(expression) {
  it_capture(substitute(expression), parent.frame())
}
