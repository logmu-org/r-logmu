# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

includes_class <- "includes"

#' Includes: a collection of includes
#'
#' @description
#' An `includes` is an optionally named list of [include]s. It is what the
#' `breakdown` argument of an analysis takes, and it is what every plural
#' constructor -- [bands()], [ages()], [durations()], [periods()] -- returns.
#'
#' It is a separate type rather than an [include] with a length, because giving
#' `include` a length would touch every working include path to buy generality
#' only collections need.
#'
#' Each element carries **two levels of naming**. Its own `names()` entry says
#' which band it is, and its [group_names()] entry says what was banded. A
#' fifteen-element breakdown of ages, periods and amounts can then tell a chart
#' which six are the ages, which a single flat name could not.
#'
#' Elements need not be disjoint. Examining A/E on intersecting subsets is
#' legitimate, so **logmu** does not check for overlaps and does not treat a
#' breakdown as a partition.
#'
#' @param ... `include` and `includes` objects to collect. An `includes`
#'   argument is flattened in. A name given to an `include` argument replaces
#'   its name; a name given to an `includes` argument replaces the group name of
#'   every element it contributes.
#' An `includes` is subsettable with `[` and `[[`, which take an index and keep
#' the labels of whatever they select.
#'
#' Both labels may be replaced. `names<-` takes exactly one name per include;
#' `group_names<-` takes one per include or a single value for all of them.
#'
#' @param x An `includes`.
#' @returns
#' `includes()` returns an `includes`.
#'
#' `is_includes()` returns a scalar `logical`.
#'
#' `group_names()` returns a `character` vector with one element per include.
#'
#' `names<-()` and `group_names<-()` return a new `includes`.
#' @examples
#' std <- includes(
#'   ages(65, 95, by = 5),
#'   periods(2000, 2020, by = 10)
#' )
#' std
#' names(std)
#' group_names(std)
#'
#' # Naming an includes argument renames the group of everything in it.
#' includes(cohort = ages(65, 95, by = 15))
#' @name includes
NULL

# THE ELEMENTS ARE THE SINGLE SOURCE OF TRUTH for both levels of naming. The
# list names are a cache filled from them here, which is what makes the base
# machinery -- printing, subsetting, data.frame interop -- work for free.
#
# Any name written to the cache alone would revert the moment anything rebuilt
# from the elements, so `names<-` below writes through rather than leaving the
# two to drift. Incoming names are discarded first, so the result never depends
# on what the caller's list happened to be carrying.
new_includes <- function(items) {
  items <- unname(items)
  named <- vapply(items, function(inc) inc$name %||% "", character(1L))
  if (any(nzchar(named))) names(items) <- named
  structure(items, class = includes_class)
}

# An include is immutable (`$<-` is blocked on every logmu_function), so a
# relabel builds a new one rather than editing in place. An indicator is an
# include but holds an AST rather than terms, and relabelling it must not
# demote it to a plain include -- that would lose the {0,1} guarantee V = E
# rests on.
relabel_include <- function(x, group_name = x$group_name, name = x$name) {
  if (is_indicator(x)) return(new_indicator(x$ast, group_name, name))
  new_include(x$terms, group_name, name)
}

#' @rdname includes
#' @export
includes <- function(...) {
  args <- list(...)
  given <- names(args)

  items <- list()
  for (i in seq_along(args)) {
    argument <- args[[i]]
    label <- if (!is.null(given) && nzchar(given[[i]])) given[[i]] else NULL

    if (is_includes(argument)) {
      contributed <- unname(unclass(argument))
      if (!is.null(label)) {
        contributed <- lapply(contributed, relabel_include, group_name = label)
      }
      items <- c(items, contributed)
    } else if (is_include(argument)) {
      if (!is.null(label)) argument <- relabel_include(argument, name = label)
      items <- c(items, list(argument))
    } else {
      stop("Every argument to `includes()` must be an `include` or an `includes`.",
           call. = FALSE)
    }
  }

  new_includes(items)
}

#' @rdname includes
#' @export
is_includes <- function(x) inherits(x, "includes")

ensure_is_includes <- function(x) {
  if (!is_includes(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be an `includes`.", call. = FALSE)
  }
}

# Coerce whatever a breakdown argument was given as into an `includes`. A single
# include is the one-element collection, which is what lets `breakdown =` take
# either without the caller checking.
as_includes <- function(x) {
  if (is_includes(x)) return(x)
  if (is_include(x)) return(new_includes(list(x)))
  stop("A breakdown must be an `include` or an `includes`.", call. = FALSE)
}

# GENERIC, because an `aev` carries the same two levels of naming once a
# breakdown has produced it -- the labels travel from the includes to the result.
#' @rdname includes
#' @export
group_names <- function(x) UseMethod("group_names")

#' @rdname includes
#' @usage NULL
#' @export
group_names.default <- function(x) {
  stop("S3 function `group_names` is not implemented for `",
       deparse(substitute(x)), "`.", call. = FALSE)
}

#' @rdname includes
#' @usage NULL
#' @export
# Unnamed deliberately: this is a vector PARALLEL to `names(x)`, and carrying the
# other level of naming as its own names would invite reading the wrong one.
group_names.includes <- function(x) {
  unname(vapply(unclass(x), function(inc) inc$group_name %||% NA_character_,
                character(1L)))
}

# Both setters rebuild rather than edit, because the elements are immutable --
# which is R's own copy-on-modify semantics, not a breach of the `$<-` block.
#' @rdname includes
#' @usage NULL
#' @export
`names<-.includes` <- function(x, value) {
  items <- unname(unclass(x))

  if (is.null(value)) {
    return(new_includes(lapply(items, relabel_include, name = NULL)))
  }

  value <- as.character(value)
  # EXACTLY ONE PER INCLUDE, with no recycling, deliberately. `names(x)[1] <- "y"`
  # on an includes that has no names yet produces a length-1 vector, and
  # recycling that would silently rename every element instead of the first.
  if (length(value) != length(items)) {
    stop("`names` must have one name per include (", length(items), ").", call. = FALSE)
  }
  if (anyNA(value)) stop("`names` cannot be NA.", call. = FALSE)

  new_includes(Map(function(inc, nm) relabel_include(inc, name = nm), items, value))
}

#' @rdname includes
#' @usage NULL
#' @export
`group_names<-` <- function(x, value) UseMethod("group_names<-")

#' @rdname includes
#' @usage NULL
#' @export
`group_names<-.default` <- function(x, value) {
  stop("S3 function `group_names<-` is not implemented for `",
       deparse(substitute(x)), "`.", call. = FALSE)
}

#' @rdname includes
#' @usage NULL
#' @export
`group_names<-.includes` <- function(x, value) {
  items <- unname(unclass(x))

  value <- as.character(value)
  # A SINGLE VALUE IS RECYCLED here, unlike `names<-`: naming a whole band set
  # after one thing is the ordinary case, and `group_names()` never returns NULL,
  # so the `[i] <-` idiom cannot produce a short vector by accident.
  if (length(value) == 1L) value <- rep(value, length(items))
  if (length(value) != length(items)) {
    stop("`group_names` must be a single value or one per include (",
         length(items), ").", call. = FALSE)
  }
  if (anyNA(value)) stop("`group_names` cannot be NA.", call. = FALSE)

  new_includes(Map(function(inc, nm) relabel_include(inc, group_name = nm), items, value))
}

#' @rdname includes
#' @usage NULL
#' @export
length.includes <- function(x) length(unclass(x))

#' @rdname includes
#' @usage NULL
#' @export
`[.includes` <- function(x, i) new_includes(unclass(x)[i])

#' @rdname includes
#' @usage NULL
#' @export
`[[.includes` <- function(x, i) unclass(x)[[i]]

#' @returns `x`, invisibly.
#' @export
print.includes <- function(x, ...) {
  items <- unclass(x)
  cat("<includes[", length(items), "]>\n", sep = "")
  for (inc in items) {
    cat("  ", inc$group_name %||% "?", ": ", inc$name %||% "?", "\n", sep = "")
  }
  invisible(x)
}
