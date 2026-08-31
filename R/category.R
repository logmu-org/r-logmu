# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' Categories: breaking a population down by a categorical field
#'
#' @description
#' `category()` and `categories()` build [include]s from the values of a
#' categorical field such as a scheme, a status or a sex. They complete the
#' constructor family: [bands()] cuts a continuous variable, `categories()`
#' divides a discrete one.
#'
#' Each argument in `...` is **one group**, so a character vector joins its
#' values into a single group rather than splitting into several. Names label
#' the groups.
#'
#' `.source` says where the field's values come from. With no `...` it supplies
#' the groups, one per value; with `...` it is what the written categories are
#' checked against.
#'
#' Every group is an [indicator], so \eqn{V = \Omega E} holds within each.
#'
#' @section Where the values come from:
#' `.source` may be:
#'
#' * a **factor**, giving its declared `levels()` in their declared order,
#'   including any level no record uses;
#' * a **character** vector or column, giving its distinct values in
#'   alphabetical order;
#' * a **dataset** (`exp_data`, `val_data` or a data frame), giving the column
#'   named by the field, which requires `variable` to be a plain field.
#'
#' Prefer a factor when the breakdown is to be reused across several datasets.
#' A factor's levels are a property of the column rather than of its contents,
#' so every dataset yields the same groups in the same order and the results
#' stay comparable. Two character columns need not.
#'
#' @section Checking:
#' Three checks guard the grouping, and each runs whenever it has what it
#' needs:
#'
#' * `.known` -- every category written in `...` is a value of `.source`. This
#'   catches a mistyped category.
#' * `.exhaustive` -- every value of `.source` falls in some group. This catches
#'   a value added to the data later, which would otherwise be dropped from the
#'   breakdown in silence.
#' * `.disjoint` -- no value falls in more than one group.
#'
#' The first two need `.source` and so default on when one is given; the third
#' needs only `...` and so is always on. Setting `.known` or `.exhaustive`
#' without a `.source` is an error, since there would be nothing to check
#' against.
#'
#' `.disjoint` is not the route to overlapping groups. Build those separately
#' and collect them with [includes()], which checks nothing across its
#' elements.
#'
#' @param variable A pronoun expression naming the field to divide, e.g.
#'   `.i$status`. It must not depend on time `.t`.
#' @param ... The groups. Each argument is one group; a character vector of
#'   several values is one group joining them. A name labels the group.
#' @param .source Where the field's values come from -- a factor, a character
#'   vector or column, or a dataset holding a column of the field's name. `NULL`
#'   for none. Named with a dot because every other name belongs to `...`.
#' @param .known Whether every category in `...` must be a value of `.source`.
#' @param .exhaustive Whether every value of `.source` must fall in some group.
#' @param .disjoint Whether no value may fall in more than one group.
#' @returns
#' `category()` returns an `indicator`, which is also an `include`.
#'
#' `categories()` returns an [includes()].
#' @examples
#' # One group per status, written out.
#' categories(.i$status, "act", "def", "pen")
#'
#' # A join, and named groups.
#' categories(.i$status, active = "act", other = c("def", "dep"))
#'
#' # Every level of a factor, in its declared order.
#' scheme <- factor(c("A", "B", "A"), levels = c("A", "B", "C"))
#' categories(.i$scheme, .source = scheme)
#'
#' # A single group of two joined values.
#' category(.i$status, c("def", "dep"))
#' @name category
NULL

# ---- helpers ---------------------------------------------------------------

# Categories are quoted in messages with `encodeString()` rather than
# `dQuote()`, which would emit directional quotes and so non-ASCII text.
category_quote <- function(values) {
  paste(encodeString(values, quote = "\""), collapse = ", ")
}

# One group's values. Each `...` argument is ONE group: a character vector
# joins its values rather than splitting into several groups.
category_group_values <- function(values, position) {
  where <- paste0("Argument ", position, " of `...`")
  if (!is_pure_text(values) || length(values) < 1L || anyNA(values)) {
    stop(where, " must be a character vector of one or more categories, with no `NA`s.",
         call. = FALSE)
  }
  repeated <- unique(values[duplicated(values)])
  if (length(repeated) > 0L) {
    # Across groups a repeat means overlapping groups, which `.disjoint`
    # governs because it is a thing a user can mean. Within one group it means
    # nothing at all, so it is always refused.
    stop(where, " repeats ", category_quote(repeated),
         ". A repeat within one group says nothing.", call. = FALSE)
  }
  values
}

# The field name, for the dataset form of `.source`. NULL for a computed
# variable, which has no column of its own.
category_field_name <- function(q) {
  if (it_is_field_access(q)) it_field_name(q) else NULL
}

# The field's value set. A factor gives its DECLARED levels, unused ones
# included: the level set is a property of the column rather than of its
# contents, so a breakdown shared across datasets keeps the same groups in the
# same order. A character column can only give what it happens to hold.
category_source_values <- function(source, field) {
  if (is.data.frame(source)) {
    if (is.null(field)) {
      stop("`.source` may be a dataset only when the variable is a plain field such as `.i$status`. ",
           "Pass the column itself for a computed variable.", call. = FALSE)
    }
    if (!field %in% names(source)) {
      stop(sprintf("`.source` has no column `%s`. Its columns are: %s.",
                   field, paste(names(source), collapse = ", ")),
           call. = FALSE)
    }
    source <- source[[field]]
  }

  values <- if (is.factor(source)) {
    levels(source)
  } else if (is_pure_text(source)) {
    sort(unique(source[!is.na(source)]))
  } else {
    stop("`.source` must be a factor, a character vector, or a dataset holding one.",
         call. = FALSE)
  }

  if (length(values) < 1L) {
    stop("`.source` holds no categories.", call. = FALSE)
  }
  values
}

# A group's own name defaults to what it holds: the value itself, or the
# joined values.
category_default_name <- function(values) paste(values, collapse = "+")

category_names <- function(groups) {
  defaults <- vapply(groups, category_default_name, character(1L))
  given <- names(groups)
  labels <- if (is.null(given)) defaults else ifelse(nzchar(given), given, defaults)

  repeated <- unique(labels[duplicated(labels)])
  if (length(repeated) > 0L) {
    # An `includes` carries its names through to charts and to `[[`, so two
    # groups sharing one label is a reporting trap rather than a nicety.
    stop("Two groups are named ", category_quote(repeated),
         ". Every group needs its own name.", call. = FALSE)
  }
  labels
}

# One group's membership test. A single value compares directly; a join reads
# and prints better as `%in%`, which the engine desugars to the same
# or-chain of equalities.
category_ast <- function(ast, values) {
  if (length(values) == 1L) {
    it_call("==", list(ast, it_lit(values)))
  } else {
    it_call("%in%", list(ast, it_lit(values)))
  }
}

category_check_flag <- function(value, name) {
  if (!is_single_valid_logical(value)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
}

# ---- the shared constructor ------------------------------------------------

# `.exhaustive` and `.disjoint` are FALSE from `category()`, which has one
# group: coverage would refuse the ordinary single-group call, and disjointness
# is a relation between groups.
build_categories <- function(ast, q, dots, source, known, exhaustive, disjoint) {
  if (it_uses_t(ast)) {
    stop("A category variable must not depend on time `.t`.", call. = FALSE)
  }
  category_check_flag(known, ".known")
  category_check_flag(exhaustive, ".exhaustive")
  category_check_flag(disjoint, ".disjoint")

  if (is.null(source) && (known || exhaustive)) {
    stop("`.known` and `.exhaustive` check the categories against `.source`, ",
         "so they need one to be given.", call. = FALSE)
  }

  field <- category_field_name(q)
  values <- if (is.null(source)) NULL else category_source_values(source, field)
  label <- band_group_name(q)

  if (length(dots) == 0L) {
    if (is.null(values)) {
      stop("Give one or more categories, or a `.source` to take them from.", call. = FALSE)
    }
    groups <- as.list(values)
  } else {
    groups <- lapply(seq_along(dots), function(i) category_group_values(dots[[i]], i))
    names(groups) <- names(dots)
  }

  written <- unlist(groups, use.names = FALSE)

  # `.known` is checked before `.exhaustive` because a typo trips both -- the
  # mistyped value is unknown AND the value it should have been is uncovered --
  # and only the first is actionable.
  if (known) {
    unknown <- setdiff(written, values)
    if (length(unknown) > 0L) {
      stop(sprintf("`%s` has no category %s. Set `.known = FALSE` if it is absent from this dataset only.",
                   label, category_quote(unknown)),
           call. = FALSE)
    }
  }

  if (exhaustive) {
    uncovered <- setdiff(values, written)
    if (length(uncovered) > 0L) {
      stop(sprintf("Categories of `%s` in no group: %s. Set `.exhaustive = FALSE` if they are meant to be left out.",
                   label, category_quote(uncovered)),
           call. = FALSE)
    }
  }

  if (disjoint) {
    overlapping <- unique(written[duplicated(written)])
    if (length(overlapping) > 0L) {
      stop(sprintf("Categories in more than one group: %s. Set `.disjoint = FALSE` if the groups are meant to overlap.",
                   category_quote(overlapping)),
           call. = FALSE)
    }
  }

  labels <- category_names(groups)
  lapply(seq_along(groups), function(i) {
    new_indicator(category_ast(ast, groups[[i]]), label, labels[[i]])
  })
}

# ---- the constructors ------------------------------------------------------

#' @rdname category
#' @export
category <- function(variable, ..., .source = NULL, .known = !is.null(.source)) {
  q <- substitute(variable)
  dots <- list(...)
  if (length(dots) != 1L) {
    stop(if (length(dots) == 0L) {
           "`category()` needs one category. Only `categories()` can take them from `.source`."
         } else {
           "`category()` takes one category. Use `categories()` for several, or `c()` to join them into one."
         },
         call. = FALSE)
  }
  build_categories(it_capture(q, parent.frame()), q, dots, .source, .known,
                   exhaustive = FALSE, disjoint = FALSE)[[1L]]
}

#' @rdname category
#' @export
categories <- function(variable, ..., .source = NULL,
                       .known = !is.null(.source),
                       .exhaustive = !is.null(.source),
                       .disjoint = TRUE) {
  q <- substitute(variable)
  new_includes(build_categories(it_capture(q, parent.frame()), q, list(...),
                                .source, .known, .exhaustive, .disjoint))
}
