# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# name -- an object's own name
# ============================================================================
#
# ONE RULE FOR STORAGE: an object's OWN name is its `name` attribute, and a
# collection's element names are `names()`. Two of the mortalities are atomic
# (a classed double and a matrix), so an attribute is the only place a name
# can live there, and making that the rule everywhere gives one default method
# instead of one per class. A `logmu_function` blocks `$<-` but not `attr<-`,
# so the setter works on an include or a variable without a special case.
#
# ONE READER, WRITTEN ONCE WITH `exact = TRUE`: `attr(x, "name")` partially
# matches, so on an object that has a `names` attribute and no `name` it
# silently returns the element names. Every call site going through here
# retires that trap.

#' The name of an object
#'
#' @description
#' Many **logmu** objects can carry a name of their own: a `mortality`, a
#' [model()], an [include] or [indicator] and each term of a [covariates()].
#' `name()` reads it and `name<-` sets or, with `NULL`, clears it. Where an
#' object has one, the name is what **logmu** prints for it, so a reference
#' mortality shows as `S4PMA` rather than as `<mortality_table>`, and a model
#' handed to [compare_models()] without a list name is labelled by its own.
#'
#' An include carries a second level, its `group_name`, which says which
#' breakdown it belongs to: `age` for `[60, 70)`. Over a whole [includes()]
#' the two levels are read as vectors by `names()` and [group_names()].
#'
#' # What a name may be
#'
#' A name is a single string that is not blank, has no leading, trailing or
#' consecutive spaces, and holds no control character such as a tab or a
#' newline. The rule is the same for every name in the package, whether you
#' wrote it or **logmu** generated it: `[60, 70)`, `>= 70` and `a+b` are all
#' names it makes for itself, and they follow the rule too.
#'
#' @param x An object.
#' @param value A name, or `NULL` to clear it.
#' @returns
#' `name()` and `group_name()` return a single string, or `NULL` where the
#' object has none.
#'
#' `name<-` and `group_name<-` return `x` with the name set.
#' @examples
#' base <- mortality_const(log_mu = -4, name = "flat")
#' name(base)
#'
#' m <- model(ref_mortality = base, covariates = covariates(level = 1))
#' name(m) <- "one level"
#' m
#'
#' std <- ages(65, 95, by = 15)
#' name(std[[1]])
#' group_name(std[[1]])
#' @name name
NULL

#' @rdname name
#' @export
name <- function(x) UseMethod("name")

#' @rdname name
#' @usage NULL
#' @export
name.default <- function(x) attr(x, "name", exact = TRUE)

#' @rdname name
#' @export
`name<-` <- function(x, value) UseMethod("name<-")

#' @rdname name
#' @usage NULL
#' @export
`name<-.default` <- function(x, value) {
  if (!is.null(value)) ensure_is_valid_name(value)
  attr(x, "name") <- value
  x
}

#' @rdname name
#' @export
group_name <- function(x) UseMethod("group_name")

#' @rdname name
#' @usage NULL
#' @export
group_name.default <- function(x) attr(x, "group_name", exact = TRUE)

#' @rdname name
#' @export
`group_name<-` <- function(x, value) UseMethod("group_name<-")

#' @rdname name
#' @usage NULL
#' @export
`group_name<-.default` <- function(x, value) {
  if (!is.null(value)) ensure_is_valid_name(value)
  attr(x, "group_name") <- value
  x
}
