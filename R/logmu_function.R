# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# logmu_function -- the shared base of every pronoun-function object
# ============================================================================
#
# `logmu_function` is the S3 base class of every object that represents a
# function of an individual's facts (and possibly time): `mortality`,
# `variable`, `static_variable`, `include` and `indicator`. Behaviour common
# to all of them is defined here, once, and inherited through the shared class.
#
# The members follow a common design pattern:
#   * a class-vector constant ending in "logmu_function";
#   * a validating, immutable constructor;
#   * one S3 value/inspection generic with an erroring `.default` and the
#     class's own method (the testing/education path, not the engine path);
#   * an exported `is_XXX()` defined as `inherits(x, "XXX")`;
#   * a non-exported `ensure_is_XXX()` that names the offending argument;
#   * the immutability ops and (later) `Ops` inherited from here.

# ---- immutability ----------------------------------------------------------
#
# logmu function objects are values: the list-like modification ops are
# disabled so an object cannot be corrupted in place.

#' @export
#' @noRd
`$<-.logmu_function` <- function(x, name, value) stop_logmu_immutable()

#' @export
#' @noRd
`[[<-.logmu_function` <- function(x, i, value) stop_logmu_immutable()

#' @export
#' @noRd
`[<-.logmu_function` <- function(x, i, value) stop_logmu_immutable()

stop_logmu_immutable <- function() {
  stop("A `logmu_function` object is immutable and cannot be modified in place.",
       call. = FALSE)
}

# ---- operators -------------------------------------------------------------
#
# A single `Ops` method for the whole family, switching on the operator. The
# only object-level operator is `&`, which intersects two subsets (includes
# and/or indicators) into one. Logical `|`/`!` have no object-level meaning --
# they are ordinary R operators used *inside* a pronoun expression.

#' @export
#' @noRd
Ops.logmu_function <- function(e1, e2) {
  if (identical(.Generic, "&")) return(subset_intersect(e1, e2))
  stop(sprintf("Operator `%s` is not defined for logmu function objects.", .Generic),
       call. = FALSE)
}

# ---- type predicate --------------------------------------------------------

#' Test for a logmu function object
#'
#' @description
#' Tests whether `x` is any **logmu** function object -- a `mortality`,
#' `variable`, `static_variable`, `include` or `indicator`.
#'
#' @param x Object to test.
#' @returns A scalar `logical`.
#' @export
is_logmu_function <- function(x) inherits(x, "logmu_function")

ensure_is_logmu_function <- function(x) {
  if (!is_logmu_function(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `logmu_function`.", call. = FALSE)
  }
}
