# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

category_domain_class <- "category_domain"

# Codes 0..size must all fit an R integer.
max_category_domain_size <- .Machine$integer.max - 1L

#' Category domains
#'
#' @description
#' A `category_domain` is the set of codes a [categorical][as_categorical()] can
#' take. Codes `1..size` are the categories proper. Code 0 means undefined, for
#' example a postcode that maps to no socio-economic type, and is itself a
#' valid category where `allows_undefined` is `TRUE`.
#'
#' Labels are optional. Without them each code is written as its number in
#' plain ASCII digits, whatever the locale, so `"7"` for code 7.
#'
#' Giving an `undefined_label` allows undefined, so `allows_undefined` need not
#' be set as well; setting it to `FALSE` alongside a label is an error. Where
#' undefined is allowed without a label, code 0 is written `"(undefined)"` in a
#' labelled domain and `"0"` in an unlabelled one.
#'
#' A domain is ordered if it has `sorted_codes`: every code it can hold, each
#' exactly once, from lowest to highest. That includes code 0 where undefined is
#' allowed, and code 0 may sit anywhere in the order. The order is separate from
#' the codes, so codes taken from a classification scheme can keep their numbers
#' while being ranked by, say, deprivation.
#'
#' Two domains are the same if they have the same size, labels,
#' `allows_undefined` and `sorted_codes`. Their names and undefined labels play
#' no part, so two domains built separately, even on different platforms,
#' compare equal if their content does.
#'
#' A domain cannot be modified once made.
#'
#' @param size The number of codes, not counting 0. May be omitted if `labels`
#'   is given.
#' @param labels Optional labels for codes `1..size`, in code order. Each must
#'   be a valid [name()], and they must be unique.
#' @param allows_undefined Whether code 0 is permitted. Defaults to whether
#'   `undefined_label` is given.
#' @param undefined_label The text for code 0. It must not be the text of
#'   another code.
#' @param sorted_codes For an ordered domain, every code of the domain once, in
#'   order from lowest to highest. `NULL` for an unordered domain.
#' @param name An optional name for the domain, read back by [name()].
#' @param x An object.
#' @returns
#' `category_domain()` returns a `category_domain`.
#'
#' `is_category_domain()` returns a scalar `logical`.
#' @examples
#' # Socio-economic types 1 to 60, with 0 for a postcode that maps to none.
#' types <- category_domain(60, allows_undefined = TRUE, name = "types")
#' types
#'
#' # Deprivation bands, with members whose band is unknown.
#' bands <- category_domain(labels = c("low", "mid", "high"), undefined_label = "unknown")
#' bands
#'
#' # The same bands ranked, with unknown placed between mid and high.
#' category_domain(
#'   labels = c("low", "mid", "high"),
#'   undefined_label = "unknown",
#'   sorted_codes = c(1, 2, 0, 3)
#' )
#' @export
category_domain <- function(
  size = NULL,
  labels = NULL,
  allows_undefined = !is.null(undefined_label),
  undefined_label = NULL,
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
      stop("A `category_domain` needs a `size` or `labels`.", call. = FALSE)
    }
    size <- length(labels)
  }
  ensure_is_category_domain_size(size)
  size <- as.integer(size)
  if (!is.null(labels) && length(labels) != size) {
    stop("`labels` must have one label per code, i.e. ", size, ".", call. = FALSE)
  }

  if (!is_single_valid_logical(allows_undefined)) {
    stop("`allows_undefined` must be TRUE or FALSE.", call. = FALSE)
  }
  undefined_label <- checked_undefined_label(undefined_label, allows_undefined, labels, size)
  sorted_codes <- checked_sorted_codes(sorted_codes, allows_undefined, size)

  if (!is.null(name)) ensure_is_valid_name(name)

  new_category_domain(size, labels, allows_undefined, undefined_label, sorted_codes, name)
}

new_category_domain <- function(size, labels, allows_undefined, undefined_label, sorted_codes, name) {
  structure(
    list(
      size = size,
      labels = labels,
      allows_undefined = allows_undefined,
      undefined_label = undefined_label,
      sorted_codes = sorted_codes
    ),
    class = category_domain_class,
    name = name
  )
}

ensure_is_category_domain_size <- function(size) {
  if (!is_single_pure_finite_numeric(size) || size != round(size) ||
      size < 1 || size > max_category_domain_size) {
    stop("`size` must be a whole number from 1 to ", max_category_domain_size, ".",
         call. = FALSE)
  }
}

# Code 0 has a label exactly when the domain allows it.
checked_undefined_label <- function(undefined_label, allows_undefined, labels, size) {
  if (!allows_undefined) {
    if (!is.null(undefined_label)) {
      stop("`undefined_label` is given but `allows_undefined` is FALSE.", call. = FALSE)
    }
    return(NULL)
  }
  if (is.null(undefined_label)) {
    return(if (is.null(labels)) "0" else "(undefined)")
  }
  ensure_is_valid_name(undefined_label)
  clash <- if (is.null(labels)) unlabelled_codes(undefined_label, size) else match(undefined_label, labels)
  if (!is.na(clash)) {
    stop("`undefined_label` \"", undefined_label, "\" is already the text of a code.",
         call. = FALSE)
  }
  undefined_label
}

# Every code of the domain exactly once, 0 included where it is allowed.
checked_sorted_codes <- function(sorted_codes, allows_undefined, size) {
  if (is.null(sorted_codes)) return(NULL)
  expected <- c(if (allows_undefined) 0L, seq_len(size))
  if (!is_pure_numeric(sorted_codes) || anyNA(sorted_codes) ||
      any(sorted_codes != round(sorted_codes)) ||
      !identical(sort(as.integer(sorted_codes)), expected)) {
    stop("`sorted_codes` must hold every code of the domain exactly once: 1 to ", size,
         if (allows_undefined) ", and 0" else "", ".", call. = FALSE)
  }
  as.integer(sorted_codes)
}

#' @rdname category_domain
#' @export
is_category_domain <- function(x) inherits(x, category_domain_class)

ensure_is_category_domain <- function(x) {
  if (!is_category_domain(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `category_domain`.", call. = FALSE)
  }
}

is_ordered_domain <- function(domain) is_category_domain(domain) && !is.null(domain$sorted_codes)

# The rank of each code in an ordered domain, indexed by code + 1.
category_domain_ranks <- function(domain) {
  ranks <- integer(domain$size + 1L)
  ranks[domain$sorted_codes + 1L] <- seq_along(domain$sorted_codes)
  ranks
}

# UNLABELLED CODES ARE WRITTEN AS INTEGERS, which `as.character()` gives as
# plain ASCII digits with no grouping, no leading zeros and no scientific
# notation, whatever the locale. Neither direction builds the text of every
# code, because an unlabelled domain may be very large.

# The text of each code in `codes`, with code 0 as the undefined label.
category_domain_text <- function(domain, codes) {
  codes <- as.integer(codes)
  text <- character(length(codes))
  defined <- codes != 0L
  text[!defined] <- domain$undefined_label %||% NA_character_
  text[defined] <- if (is.null(domain$labels)) {
    as.character(codes[defined])
  } else {
    domain$labels[codes[defined]]
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

# The code of each string in `text`. Text that names no code is an error: a
# mistyped category must not quietly become undefined.
category_domain_codes <- function(domain, text) {
  codes <- if (is.null(domain$labels)) {
    unlabelled_codes(text, domain$size)
  } else {
    match(text, domain$labels)
  }
  if (domain$allows_undefined) {
    codes[!is.na(text) & text == domain$undefined_label] <- 0L
  }
  unknown <- unique(text[is.na(codes)])
  if (length(unknown) > 0L) {
    shown <- utils::head(unknown, 5L)
    stop("Not a category of this domain: ",
         paste0("\"", shown, "\"", collapse = ", "),
         if (length(unknown) > length(shown)) ", ..." else "", ".", call. = FALSE)
  }
  as.integer(codes)
}

# Whether two domains are the same: equal in content, whatever their names or
# undefined labels. A `naxp` is a domain too, and equals another with the same
# pattern.
same_category_domain <- function(a, b) {
  if (is_naxp(a) || is_naxp(b)) {
    return(is_naxp(a) && is_naxp(b) && same_naxp(a, b))
  }
  a$size == b$size &&
    identical(a$labels, b$labels) &&
    a$allows_undefined == b$allows_undefined &&
    identical(a$sorted_codes, b$sorted_codes)
}

# How a domain is described in a single line: its name, else its pattern or
# size.
category_domain_label <- function(domain) {
  if (is_naxp(domain)) return(paste0("naxp ", naxp_pattern(domain)))
  name(domain) %||% paste0(domain$size, if (domain$size == 1L) " code" else " codes")
}

#' @export
#' @noRd
print.category_domain <- function(x, ...) {
  cat(sprintf(
    "<category_domain%s> %d %s, %s, %s\n",
    if (is.null(name(x))) "" else paste0(": ", name(x)),
    x$size,
    if (x$size == 1L) "code" else "codes",
    if (is_ordered_domain(x)) "ordered" else "unordered",
    if (x$allows_undefined) {
      paste0("undefined allowed as \"", x$undefined_label, "\"")
    } else {
      "undefined not allowed"
    }
  ))
  if (!is.null(x$labels)) {
    shown <- utils::head(x$labels, 6L)
    cat("Labels: ", paste(shown, collapse = ", "),
        if (x$size > length(shown)) ", ..." else "", "\n", sep = "")
  }
  if (is_ordered_domain(x)) {
    shown <- utils::head(x$sorted_codes, 6L)
    cat("Order: ", paste(category_domain_text(x, shown), collapse = " < "),
        if (length(x$sorted_codes) > length(shown)) " < ..." else "", "\n", sep = "")
  }
  invisible(x)
}

#' @export
#' @noRd
`$<-.category_domain` <- function(x, name, value) stop_table_immutable("category_domain")

#' @export
#' @noRd
`[[<-.category_domain` <- function(x, i, value) stop_table_immutable("category_domain")

#' @export
#' @noRd
`[<-.category_domain` <- function(x, i, value) stop_table_immutable("category_domain")
