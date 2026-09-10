# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# covariates -- the right-hand side of a fitted model
# ============================================================================
#
# A model is written
#
#     log mu = log mu^ref + beta_1 X^(1) + beta_2 X^(2) + ...
#
# and `covariates()` is how the X terms are given. The coefficients are the
# fitter's; nothing here names or supplies one.
#
# WHY A CONSTRUCTOR AND NOT `list(...)`. A bare `list(.i$is_male, ...)` cannot
# work: `.i` is not an exported object, so the pronoun exists only inside
# non-standard evaluation, and `Ops.logmu_function` deliberately refuses every
# operator but `&` because operators are meant to be written INSIDE a pronoun
# expression rather than between objects. Making the list form work would need
# `.i` exported and an object-level meaning for `!`, `&`, `==` and `*` -- which
# would leave the pronoun language existing in two forms, with two code paths
# that have to agree. `covariates()` captures its arguments unevaluated, so
# every operator means exactly what it means inside `aev()`: same parser, same
# path, no new semantics. It also matches `include()`, `band()`, `category()`
# and `includes()`, which are all built this way.
#
# NAMES ARE OPTIONAL AND ARE NOT REQUIRED ANYWHERE. Naming covariates is not
# usual practice and is a place for an error to hide; they are carried when
# given and ignored otherwise.

covariates_class <- c("covariates", "logmu_function")

# The items are `variable` objects and the object IS the list, which is what
# makes `length()`, `[[`, printing and `c()` work without writing any of them.
#
# `disjoint` records that the members were ASSERTED mutually exclusive. It is
# not checked here -- see `disjoint()` -- and it survives multiplication by a
# shape, because if `I^(j) I^(l) = 0` then `I^(j) phi * I^(l) phi = 0` too.
new_covariates <- function(items, disjoint = FALSE) {
  items <- unname(items)
  named <- vapply(items, function(x) attr(x, "covariate_name") %||% "", character(1L))
  if (any(nzchar(named))) names(items) <- named
  structure(items, disjoint = disjoint, class = covariates_class)
}

# A name rides on the item rather than on the list, so it survives the algebra:
# a Cartesian product has to build new names from the two it came from, and a
# names attribute on the outer list would be lost the moment anything rebuilt.
label_covariate <- function(x, name) {
  if (is.null(name) || !nzchar(name)) return(x)
  attr(x, "covariate_name") <- name
  x
}

covariate_name_of <- function(x) attr(x, "covariate_name") %||% ""

# Capture one unevaluated argument as a variable.
#
# A SYMBOL BOUND TO A `covariates` IS SPLICED. Everything else goes through
# `it_capture()`, which already splices a symbol bound to a `variable` and
# parses a `~` formula, so no separate path is needed for either. Only a symbol
# is tried, never a call, so nothing that could be a pronoun expression is
# evaluated by accident.
capture_covariate_argument <- function(argument, env) {
  if (is.symbol(argument)) {
    value <- tryCatch(eval(argument, env), error = function(e) NULL)
    if (is_covariates(value)) {
      return(list(items = unname(unclass(value)), spliced = TRUE))
    }
  }
  list(items = list(it_specialise(it_capture(argument, env))), spliced = FALSE)
}

build_covariates <- function(arguments, env) {
  given <- names(arguments)
  items <- list()

  for (i in seq_along(arguments)) {
    label <- if (!is.null(given) && nzchar(given[[i]])) given[[i]] else NULL
    captured <- capture_covariate_argument(arguments[[i]], env)

    # A label on a spliced group would have to be shared by all of it, which is
    # the two-level naming `includes` needs and a covariate list does not.
    # Judged by whether it WAS a splice rather than by how many it brought: a
    # group of one is still a group, and naming it is still the wrong thing.
    if (!is.null(label)) {
      if (captured$spliced) {
        stop(sprintf("`%s` names a group of %d covariates; name them individually.",
                     label, length(captured$items)),
             call. = FALSE)
      }
      captured$items[[1L]] <- label_covariate(captured$items[[1L]], label)
    }
    items <- c(items, captured$items)
  }

  items
}

#' Covariates: the fitted terms of a mortality model
#'
#' @description
#' `covariates()` collects the \eqn{X} terms of a proportional hazards model,
#' \eqn{\log\mu = \log\mu^{\mathrm{ref}} + \sum_j \beta_j X^{(j)}}. Each
#' argument is a pronoun expression, and the coefficients belong to the fit.
#'
#' `disjoint()` is the same, and additionally asserts that the terms are
#' mutually exclusive indicators -- no individual is in more than one. That is
#' worth stating: where it holds, every off-diagonal entry of the information
#' matrix is zero, and the fit does far less work.
#'
#' Multiplying by a variable distributes over the list, which is how the common
#' shape \eqn{X^{(k)} = I^{(k)} \varphi} is written: a set of indicators times a
#' shared age shape. Multiplying two covariate lists gives their Cartesian
#' product.
#'
#' @param ... Pronoun expressions, optionally named. A `covariates` object is
#'   spliced in.
#' @param x An object.
#' @returns
#' `covariates()` and `disjoint()` return a `covariates`.
#'
#' `is_covariates()` returns a scalar `logical`.
#'
#' `is_disjoint()` returns a scalar `logical` saying whether the terms were
#' asserted mutually exclusive.
#' @examples
#' age_shape <- variable(.x)
#' covariates(.i$sex == "M", .i$sex == "F") * age_shape
#' @name covariates
NULL

#' @rdname covariates
#' @export
covariates <- function(...) {
  new_covariates(build_covariates(match.call(expand.dots = FALSE)$..., parent.frame()))
}

#' @rdname covariates
#' @export
disjoint <- function(...) {
  items <- build_covariates(match.call(expand.dots = FALSE)$..., parent.frame())

  # WHAT R CAN CHECK, AND WHAT IT CANNOT. Whether a column holds a logical is
  # not known until the data arrives, so `.i$is_male` reaches here as a
  # `static_variable` and nothing here can prove it is an indicator -- the same
  # split as time-invariance being checked in R and category-ness in the engine.
  #
  # What R can prove is that a term uses `.t`, and a time-varying member is
  # refused. The whole point of the pattern is an indicator on the INDIVIDUAL
  # times a shape over time: a membership that moved during an individual's
  # exposure would not hoist out of the integral, and the saving would be lost
  # along with the meaning.
  for (i in seq_along(items)) {
    if (!is_static_variable(items[[i]])) {
      stop(sprintf(paste0("Every term in `disjoint()` must be an indicator that does not ",
                          "vary with time; term %d uses `.t`."), i),
           call. = FALSE)
    }
  }

  new_covariates(items, disjoint = TRUE)
}

#' @rdname covariates
#' @export
is_covariates <- function(x) inherits(x, "covariates")

#' @rdname covariates
#' @export
is_disjoint <- function(x) isTRUE(attr(x, "disjoint"))

ensure_is_covariates <- function(x) {
  if (!is_covariates(x)) {
    stop(sprintf("`%s` must be a `covariates` object.",
                 deparse(substitute(x))), call. = FALSE)
  }
  x
}

# ---- the algebra -----------------------------------------------------------
#
# THE PRODUCT OF TWO ASTs IS BUILT DIRECTLY rather than by re-parsing, because
# there is nothing to parse: both sides are already the trees the front end
# produces, and `it_specialise()` then re-derives the narrowest type for the
# result. A product of two indicators is an indicator; of a static variable and
# a time-varying one, time-varying.

multiply_asts <- function(left, right) {
  it_specialise(it_call("*", list(left$ast, right$ast)))
}

# One covariate times one scalar function, keeping the covariate's name.
scale_covariate <- function(item, factor) {
  label_covariate(multiply_asts(item, factor), covariate_name_of(item))
}

# Every covariate times the same scalar function. THE DISJOINTNESS SURVIVES:
# a common factor cannot make two mutually exclusive terms overlap.
distribute_over <- function(group, factor) {
  new_covariates(lapply(unclass(group), scale_covariate, factor = factor),
                 disjoint = is_disjoint(group))
}

# THE CARTESIAN PRODUCT, which is what crossing two breakdowns means: region by
# sex is every region with every sex.
#
# THE RESULT IS DISJOINT IF EITHER SIDE IS. If no individual is in two regions
# then no individual is in two region-sex cells, whatever the sexes do. Only one
# side needs to separate them.
cross_covariates <- function(left, right) {
  items <- list()
  for (a in unclass(left)) {
    for (b in unclass(right)) {
      # NAMED ONLY WHEN BOTH SIDES ARE. Joining a name to an empty one would give
      # every product of that row the same label -- `male` crossed with two
      # unnamed smoking states would be `male` twice -- which is worse than no
      # name at all, since a name is meant to identify one coefficient.
      leftName <- covariate_name_of(a)
      rightName <- covariate_name_of(b)
      name <- if (nzchar(leftName) && nzchar(rightName)) {
        paste(leftName, rightName, sep = ".")
      } else {
        ""
      }
      items <- c(items, list(label_covariate(multiply_asts(a, b), name)))
    }
  }
  new_covariates(items, disjoint = is_disjoint(left) || is_disjoint(right))
}

# A plain number is a rescaling of every term, which is meaningful and free.
as_scalar_function <- function(x, side) {
  if (is_covariates(x)) return(NULL)
  if (is_variable(x)) return(x)
  if (is.numeric(x) && length(x) == 1L && !is.na(x)) {
    return(it_specialise(it_lit(as.double(x))))
  }
  stop(sprintf(paste0("The %s of `*` must be a `covariates`, a `variable` or a single ",
                      "number."), side),
       call. = FALSE)
}

covariates_multiply <- function(e1, e2) {
  if (is_covariates(e1) && is_covariates(e2)) return(cross_covariates(e1, e2))
  if (is_covariates(e1)) return(distribute_over(e1, as_scalar_function(e2, "right side")))
  distribute_over(e2, as_scalar_function(e1, "left side"))
}

# ---- base methods ----------------------------------------------------------

#' @export
#' @noRd
c.covariates <- function(...) {
  parts <- list(...)
  for (part in parts) ensure_is_covariates(part)

  # DISJOINTNESS IS A PROPERTY OF A GROUP AND DOES NOT COMPOSE. Two internally
  # exclusive sets say nothing about each other, so a concatenation asserts
  # nothing -- exactly as `includes()` composes separately-built includes and
  # checks nothing across them. Claiming otherwise would be asserting something
  # nobody stated and the engine would then not check.
  new_covariates(unlist(lapply(parts, function(p) unname(unclass(p))), recursive = FALSE),
                 disjoint = FALSE)
}

#' @export
#' @noRd
`[.covariates` <- function(x, i) new_covariates(unclass(x)[i], disjoint = is_disjoint(x))

#' @export
#' @noRd
print.covariates <- function(x, ...) {
  cat(sprintf("<covariates: %d term%s%s>\n",
              length(x), if (length(x) == 1L) "" else "s",
              if (is_disjoint(x)) ", disjoint" else ""))
  labels <- names(x)
  for (i in seq_along(x)) {
    label <- if (!is.null(labels) && nzchar(labels[[i]])) paste0(labels[[i]], ": ") else ""
    cat(sprintf("  [%d] %s%s\n", i, label, it_deparse(unclass(x)[[i]]$ast)))
  }
  invisible(x)
}
