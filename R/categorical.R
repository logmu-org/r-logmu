# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# categorical -- the domain of a categorical_encoded vector
# ============================================================================
#
# THE NAMES MIRROR NAXP. There the bare name is the thing declared (`naxp`) and
# the vector carries a suffix (`naxp_encoded`), so here the domain is a
# `categorical` and the vector a `categorical_encoded`. `encode()` and
# `decode()` serve both kinds of domain.
#
# CODE 0 IS ALWAYS INVALID, as it is for a naxp. The categories proper are
# codes 1..size. To treat unknowns as a category, a user declares one and maps
# to it; nothing is a category by default.

categorical_class <- "categorical"

# Codes 0..size, and positions 1..size + 1 of a dense map, must all fit an R
# integer.
max_categorical_size <- .Machine$integer.max - 1L

#' Categorical domains
#'
#' @description
#' A `categorical` is a domain: the fixed set of categories a
#' [categorical_encoded][encode()] vector can hold, with their labels and, if
#' it has one, their order. It is defined once and reused, so every dataset,
#' map and platform gives the same text the same code.
#'
#' The categories are codes `1..size`. Code 0 is always invalid, meaning no
#' valid value, as it is for a naxp. To treat unknown values as a category of
#' their own, declare one, such as `"unknown"`, and map to it.
#'
#' Labels are optional. Without them each code is written as its number in
#' plain ASCII digits, whatever the locale, so `"7"` for code 7.
#'
#' A domain is ordered if it has `sorted_codes`: every category, each exactly
#' once, from lowest to highest. The order is separate from the codes, so codes
#' taken from a classification scheme can keep their numbers while being ranked
#' by, say, deprivation.
#'
#' Two domains are the same if they have the same size, labels and
#' `sorted_codes`. Their names play no part, so two domains built separately,
#' even on different platforms, compare equal if their content does.
#'
#' A domain cannot be modified once made.
#'
#' @param size The number of categories. May be omitted if `labels` is given.
#' @param labels Optional labels for codes `1..size`, in code order. Each must
#'   be a valid [name()], and they must be unique.
#' @param sorted_codes For an ordered domain, every code `1..size` once, in
#'   order from lowest to highest. `NULL` for an unordered domain.
#' @param name An optional name for the domain, read back by [name()].
#' @param x An object.
#' @returns
#' `categorical()` returns a `categorical`.
#'
#' `is_categorical()` returns a scalar `logical`.
#' @examples
#' # Socio-economic types 1 to 60. Code 0 is a postcode that maps to none.
#' types <- categorical(60, name = "types")
#' types
#'
#' sex <- categorical(labels = c("female", "male"), name = "sex")
#' encode(c("male", NA, "female"), sex)
#'
#' # Deprivation bands, ranked, with unknown as a category of its own.
#' categorical(
#'   labels = c("low", "unknown", "mid", "high"),
#'   sorted_codes = c(1, 3, 4, 2)
#' )
#' @export
categorical <- function(
  size = NULL,
  labels = NULL,
  sorted_codes = NULL,
  name = NULL) {

  if (!is.null(labels)) {
    if (!is_pure_text(labels) || length(labels) == 0L) {
      stop("`labels` must be a character vector with at least one label.", call. = FALSE)
    }
    ensure_are_valid_names(labels, "labels")
    duplicated_labels <- unique(labels[duplicated(labels)])
    if (length(duplicated_labels) > 0L) {
      stop("`labels` must be unique; repeated: ",
           paste0("\"", duplicated_labels, "\"", collapse = ", "), ".", call. = FALSE)
    }
  }

  if (is.null(size)) {
    if (is.null(labels)) {
      stop("A `categorical` needs a `size` or `labels`.", call. = FALSE)
    }
    size <- length(labels)
  }
  ensure_is_categorical_size(size)
  size <- as.integer(size)
  if (!is.null(labels) && length(labels) != size) {
    stop("`labels` must have one label per code, i.e. ", size, ".", call. = FALSE)
  }

  sorted_codes <- checked_sorted_codes(sorted_codes, size)

  if (!is.null(name)) ensure_is_valid_name(name)

  new_categorical(size, labels, sorted_codes, name)
}

new_categorical <- function(size, labels, sorted_codes, name) {
  structure(
    list(
      size = size,
      labels = labels,
      sorted_codes = sorted_codes
    ),
    class = categorical_class,
    name = name
  )
}

ensure_is_categorical_size <- function(size) {
  if (!is_single_pure_finite_numeric(size) || size != round(size) ||
      size < 1 || size > max_categorical_size) {
    stop("`size` must be a whole number from 1 to ", max_categorical_size, ".",
         call. = FALSE)
  }
}

# Every category exactly once. Code 0 is not a category, so it has no place.
checked_sorted_codes <- function(sorted_codes, size) {
  if (is.null(sorted_codes)) return(NULL)
  if (!is_pure_numeric(sorted_codes) || anyNA(sorted_codes) ||
      any(sorted_codes != round(sorted_codes)) ||
      !identical(sort(as.integer(sorted_codes)), seq_len(size))) {
    stop("`sorted_codes` must hold every code from 1 to ", size, " exactly once.",
         call. = FALSE)
  }
  as.integer(sorted_codes)
}

#' @rdname categorical
#' @export
is_categorical <- function(x) inherits(x, categorical_class)

ensure_is_categorical <- function(x) {
  if (!is_categorical(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `categorical`.", call. = FALSE)
  }
}

is_ordered_domain <- function(domain) is_categorical(domain) && !is.null(domain$sorted_codes)

# The rank of each code in an ordered domain, indexed by code + 1. Code 0 is
# invalid and ranks below every category, as it does in naxp's `xtfrm()`.
category_ranks <- function(domain) {
  ranks <- integer(domain$size + 1L)
  ranks[domain$sorted_codes + 1L] <- seq_along(domain$sorted_codes)
  ranks
}

# UNLABELLED CODES ARE WRITTEN AS INTEGERS, which `as.character()` gives as
# plain ASCII digits with no grouping, no leading zeros and no scientific
# notation, whatever the locale. Neither direction builds the text of every
# code, because an unlabelled domain may be very large.

# The text of each code in `codes`, with NA for code 0, which has none.
category_text <- function(domain, codes) {
  codes <- as.integer(codes)
  text <- rep(NA_character_, length(codes))
  valid <- codes != 0L
  text[valid] <- if (is.null(domain$labels)) {
    as.character(codes[valid])
  } else {
    domain$labels[codes[valid]]
  }
  text
}

# The code 1..size an unlabelled domain writes as each string, or NA. Only the
# canonical form counts: "07" and "7.0" are not code 7.
unlabelled_codes <- function(text, size) {
  codes <- rep(NA_integer_, length(text))
  canonical <- !is.na(text) & grepl("^[1-9][0-9]{0,9}$", text)
  number <- as.numeric(text[canonical])
  codes[canonical] <- ifelse(number <= size, as.integer(pmin(number, size)), NA_integer_)
  codes
}

# The code of each string in `text`, with NA as code 0. Text that names no
# category is an error: a mistyped category must not quietly become invalid.
category_codes <- function(domain, text) {
  codes <- if (is.null(domain$labels)) {
    unlabelled_codes(text, domain$size)
  } else {
    match(text, domain$labels)
  }
  unknown <- unique(text[is.na(codes) & !is.na(text)])
  if (length(unknown) > 0L) {
    shown <- utils::head(unknown, 5L)
    stop("Not a category of this domain: ",
         paste0("\"", shown, "\"", collapse = ", "),
         if (length(unknown) > length(shown)) ", ..." else "", ".", call. = FALSE)
  }
  codes[is.na(codes)] <- 0L
  as.integer(codes)
}

# Whether two domains are the same: equal in content, whatever their names. A
# `naxp` is a domain too, and equals another with the same pattern.
same_domain <- function(a, b) {
  if (is_naxp(a) || is_naxp(b)) {
    return(is_naxp(a) && is_naxp(b) && same_naxp(a, b))
  }
  a$size == b$size &&
    identical(a$labels, b$labels) &&
    identical(a$sorted_codes, b$sorted_codes)
}

# How a domain is described in a single line: its name, else its pattern or
# size.
domain_label <- function(domain) {
  if (is_naxp(domain)) return(paste0("naxp ", naxp_pattern(domain)))
  name(domain) %||% paste0(domain$size, if (domain$size == 1L) " category" else " categories")
}

#' @export
#' @noRd
print.categorical <- function(x, ...) {
  cat(sprintf(
    "<categorical%s> %d %s, %s\n",
    if (is.null(name(x))) "" else paste0(": ", name(x)),
    x$size,
    if (x$size == 1L) "category" else "categories",
    if (is_ordered_domain(x)) "ordered" else "unordered"
  ))
  if (!is.null(x$labels)) {
    shown <- utils::head(x$labels, 6L)
    cat("Labels: ", paste(shown, collapse = ", "),
        if (x$size > length(shown)) ", ..." else "", "\n", sep = "")
  }
  if (is_ordered_domain(x)) {
    shown <- utils::head(x$sorted_codes, 6L)
    cat("Order: ", paste(category_text(x, shown), collapse = " < "),
        if (length(x$sorted_codes) > length(shown)) " < ..." else "", "\n", sep = "")
  }
  invisible(x)
}

#' @export
#' @noRd
`$<-.categorical` <- function(x, name, value) stop_table_immutable("categorical")

#' @export
#' @noRd
`[[<-.categorical` <- function(x, i, value) stop_table_immutable("categorical")

#' @export
#' @noRd
`[<-.categorical` <- function(x, i, value) stop_table_immutable("categorical")
