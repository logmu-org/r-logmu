# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# categorical_encoded -- logmu's own categorical vector, and encode()/decode()
# ============================================================================
#
# NOT A FACTOR, AND NOT A SUBCLASS OF ONE. A factor's missing value is R's NA,
# which no other platform shares, and a subclass of factor is worse, because
# forcats keeps the subclass while reordering and merging the levels underneath
# it. So a categorical_encoded vector stores its own codes 0..size, exactly as
# the engine does, and converts to and from a factor at an explicit boundary.
# The pattern is naxp's `naxp_encoded`.
#
# CODE 0 IS INVALID AND NA IS NEVER STORED. NA is translated at the boundary:
# `encode()` turns it into 0, `decode()` turns 0 back into NA, and `is.na()` is
# TRUE for 0. Everything that could leave an NA in the codes -- reading past
# the end, lengthening, assigning NA -- stores 0 instead.
#
# NO NA LOGIC. logmu is cross-platform, so nothing here relies on R's NA
# semantics. Comparisons give plain TRUE or FALSE: code 0 compares as an
# ordinary code and ranks below every category. Comparing with NA is an error
# that points to `is.na()`.
#
# `encode(x, domain)` DISPATCHES ON THE DOMAIN, so one verb serves a
# `categorical` and a `naxp` alike. It takes text or a factor only. Codes go
# through `categorical_encoded()`, mirroring `naxp::naxp_encoded()`, because a
# number given to `encode()` could be a code or the text it stands for.

categorical_encoded_class <- "categorical_encoded"

# How code 0 is shown, the same as naxp shows it.
invalid_marker <- "<invalid (0)>"

#' Encoding text as categories
#'
#' @description
#' `encode()` turns text, or a factor, into codes in a domain: a
#' [categorical()], giving a `categorical_encoded` vector, or a naxp from the
#' naxp package, giving a `naxp_encoded` vector. `decode()` turns codes back
#' into text.
#'
#' Code 0 is invalid. `encode()` gives code 0 for `NA`, `decode()` gives `NA`
#' for code 0, and [is.na()] is `TRUE` for it. For a `categorical`, text that
#' names no category of the domain is an error rather than becoming invalid,
#' so that a mistyped category cannot vanish. For a naxp, text the pattern
#' does not accept encodes to 0, as [naxp::naxp_encode()] does.
#'
#' `encode()` always needs a domain. A factor never makes its own: give it one
#' with `encode(f, categorical(labels = levels(f)))`.
#'
#' `categorical_encoded()` makes a vector from codes that are already numbers,
#' such as socio-economic types supplied as integers. `encode()` refuses
#' numbers and points to it.
#'
#' A `categorical_encoded` vector prints as its labels, with code 0 shown as
#' `<invalid (0)>`. It compares by `==` and `!=` with another in the same
#' domain or with text. In an ordered domain `<`, `<=`, `>` and `>=` compare by
#' the domain's order, as do [sort()] and [order()]. Code 0 compares as an
#' ordinary code and ranks below every category, so comparisons always give
#' `TRUE` or `FALSE`, never `NA`; test for invalid values with [is.na()]. This
#' differs from a factor, deliberately, because logmu behaves the same on every
#' platform. Comparing with `NA` is an error. Arithmetic is an error. [c()]
#' combines vectors in the same domain.
#'
#' `%in%` and [match()] compare text, not codes, because R gives each side to
#' them separately. So `x %in% "act"` works, but a vector from another domain
#' is not refused as it is by `==`.
#'
#' `factor(x)` gives a factor with the levels that occur, as it does for any
#' vector, in the domain's order. To convert with every category of the domain
#' as a level, use `forcats::as_factor(x)`: it gives the levels in the domain's
#' order, an ordered factor for an ordered domain, and `NA` for code 0.
#'
#' `domain()` returns the domain of a `categorical_encoded` or a
#' `naxp_encoded` vector.
#'
#' @param x For `encode()`, text or a factor. For the others, an object.
#' @param domain A [categorical()], or a naxp.
#' @param codes Whole numbers from 0 to the domain's size.
#' @returns
#' `encode()` returns a `categorical_encoded` vector for a `categorical`
#' domain, and a `naxp_encoded` vector for a naxp. `categorical_encoded()`
#' returns a `categorical_encoded` vector.
#'
#' `decode()` returns a character vector.
#'
#' `is_categorical_encoded()` returns a scalar `logical`.
#'
#' `domain()` returns a `categorical`, or a naxp.
#' @examples
#' sex <- categorical(labels = c("female", "male"), name = "sex")
#' x <- encode(c("male", NA, "female"), sex)
#' x
#' is.na(x)
#' decode(x)
#' x == "male"
#'
#' # A factor needs a domain too.
#' status <- factor(c("act", "def", "act"))
#' encode(status, categorical(labels = levels(status)))
#'
#' # Codes supplied as numbers, with 0 for none.
#' types <- categorical(60, name = "types")
#' categorical_encoded(c(12, 0, 60), types)
#' @export
encode <- function(x, domain) {
  if (missing(domain)) {
    stop("`encode()` needs a `domain`: a `categorical` or a naxp.", call. = FALSE)
  }
  UseMethod("encode", domain)
}

#' @export
#' @noRd
encode.default <- function(x, domain) {
  stop("`domain` must be a `categorical` or a naxp.", call. = FALSE)
}

#' @export
#' @noRd
encode.categorical <- function(x, domain) {
  if (is_categorical_encoded(x)) {
    if (!same_domain(domain(x), domain)) {
      stop("These codes are in a different domain.", call. = FALSE)
    }
    return(x)
  }
  text <- encodable_text(x, "categorical_encoded()")
  new_categorical_encoded(category_codes(domain, text), domain)
}

#' @export
#' @noRd
encode.naxp <- function(x, domain) {
  if (is_naxp_encoded(x)) {
    if (!same_naxp(naxp_of(x), domain)) {
      stop("These codes were encoded with a different naxp.", call. = FALSE)
    }
    return(x)
  }
  text <- encodable_text(x, "naxp::naxp_encoded()")
  naxp_encode_text(text, domain)
}

# The text `encode()` works on. A factor is read through its labels, and a
# vector that is all NA is text that is all missing.
encodable_text <- function(x, codes_constructor) {
  if (is.factor(x)) return(as.character(x))
  if (is.character(x) && !is.object(x)) return(unname(x))
  if (is.logical(x) && !is.object(x) && all(is.na(x))) return(rep(NA_character_, length(x)))
  if (is.numeric(x)) {
    stop("`encode()` takes text or a factor. For codes, use ", codes_constructor, ".",
         call. = FALSE)
  }
  stop("`encode()` takes text or a factor, not an object of class `", class(x)[[1L]], "`.",
       call. = FALSE)
}

#' @rdname encode
#' @export
decode <- function(x) UseMethod("decode")

#' @export
#' @noRd
decode.default <- function(x) {
  stop("Only a `categorical_encoded` or a `naxp_encoded` vector can be decoded.", call. = FALSE)
}

#' @export
#' @noRd
decode.categorical_encoded <- function(x) category_text(domain(x), categorical_codes(x))

#' @export
#' @noRd
decode.naxp_encoded <- function(x) naxp_decode_codes(x)

#' @rdname encode
#' @export
categorical_encoded <- function(codes, domain) {
  ensure_is_categorical(domain)
  new_categorical_encoded(checked_category_codes(codes, domain), domain)
}

# Whole numbers in 0..size, with no NA: a number has no missing value to give.
checked_category_codes <- function(codes, domain) {
  if (!is_pure_numeric(codes) || anyNA(codes) || any(codes != round(codes)) ||
      any(codes < 0) || any(codes > domain$size)) {
    stop("Codes must be whole numbers from 0 to ", domain$size, ", with no NA.", call. = FALSE)
  }
  codes <- as.integer(codes)
  attributes(codes) <- NULL
  codes
}

new_categorical_encoded <- function(codes, domain) {
  structure(codes, domain = domain, class = categorical_encoded_class)
}

#' @rdname encode
#' @export
is_categorical_encoded <- function(x) inherits(x, categorical_encoded_class)

#' @rdname encode
#' @export
domain <- function(x) UseMethod("domain")

#' @export
#' @noRd
domain.default <- function(x) {
  stop("Only a `categorical_encoded` or a `naxp_encoded` vector has a domain.", call. = FALSE)
}

#' @export
#' @noRd
domain.categorical_encoded <- function(x) attr(x, "domain", exact = TRUE)

#' @export
#' @noRd
domain.naxp_encoded <- function(x) naxp_of(x)

# The codes without class, domain or names.
categorical_codes <- function(x) {
  codes <- unclass(x)
  attributes(codes) <- NULL
  codes
}

# Codes to compare with a vector in `domain`. NA is refused: comparisons give
# TRUE or FALSE only, and `is.na()` is how invalid values are found.
comparison_codes <- function(value, domain) {
  if (is_categorical_encoded(value)) {
    if (!same_domain(domain(value), domain)) {
      stop("These categorical_encoded vectors are in different domains.", call. = FALSE)
    }
    return(categorical_codes(value))
  }
  if (anyNA(value)) {
    stop("Comparing with NA is not defined: an invalid value is code 0, and `is.na()` ",
         "finds it.", call. = FALSE)
  }
  if (is.character(value)) {
    return(category_codes(domain, value))
  }
  stop("A categorical_encoded vector can be compared only with one in the same domain, ",
       "or with text.", call. = FALSE)
}

# Codes to assign into a vector in `domain`. NA stores code 0, as encoding NA
# does.
assignment_codes <- function(value, domain) {
  if (is_categorical_encoded(value)) {
    if (!same_domain(domain(value), domain)) {
      stop("These categorical_encoded vectors are in different domains.", call. = FALSE)
    }
    return(categorical_codes(value))
  }
  if (is.character(value)) {
    return(category_codes(domain, value))
  }
  if (is.logical(value) && all(is.na(value))) {
    return(integer(length(value)))
  }
  stop("Assign a categorical_encoded vector in the same domain, text to encode, or NA.",
       call. = FALSE)
}

# Positions gained by lengthening, or by indexing out of range or with NA, hold
# NA; codes never do, so they hold code 0.
fill_zero <- function(codes) {
  if (anyNA(codes)) codes[is.na(codes)] <- 0L
  codes
}

#' @export
#' @noRd
format.categorical_encoded <- function(x, ...) {
  text <- decode.categorical_encoded(x)
  text[is.na(text)] <- invalid_marker
  text
}

#' @export
#' @noRd
print.categorical_encoded <- function(x, ...) {
  cat("<categorical_encoded[", length(x), "]> ", domain_label(domain(x)), "\n", sep = "")
  if (length(x) > 0L) {
    print(format(x), quote = FALSE)
  }
  invisible(x)
}

#' @export
#' @noRd
as.character.categorical_encoded <- function(x, ...) decode.categorical_encoded(x)

#' @export
#' @noRd
is.na.categorical_encoded <- function(x) categorical_codes(x) == 0L

#' @export
#' @noRd
anyNA.categorical_encoded <- function(x, recursive = FALSE) any(categorical_codes(x) == 0L)

#' @export
#' @noRd
as.data.frame.categorical_encoded <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  nm = deparse1(substitute(x))) {
  as.data.frame.vector(x, row.names = row.names, optional = optional, ..., nm = nm)
}

#' @export
#' @noRd
`[.categorical_encoded` <- function(x, ...) {
  new_categorical_encoded(fill_zero(categorical_codes(NextMethod())), domain(x))
}

#' @export
#' @noRd
`[[.categorical_encoded` <- function(x, ...) {
  new_categorical_encoded(categorical_codes(NextMethod()), domain(x))
}

#' @export
#' @noRd
`[<-.categorical_encoded` <- function(x, ..., value) {
  domain <- domain(x)
  value <- assignment_codes(value, domain)
  codes <- categorical_codes(x)
  codes[...] <- value
  new_categorical_encoded(fill_zero(codes), domain)
}

#' @export
#' @noRd
`[[<-.categorical_encoded` <- function(x, ..., value) {
  domain <- domain(x)
  value <- assignment_codes(value, domain)
  codes <- categorical_codes(x)
  codes[[...]] <- value
  new_categorical_encoded(fill_zero(codes), domain)
}

#' @export
#' @noRd
`length<-.categorical_encoded` <- function(x, value) {
  codes <- categorical_codes(x)
  length(codes) <- value
  new_categorical_encoded(fill_zero(codes), domain(x))
}

#' @export
#' @noRd
c.categorical_encoded <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  domain <- domain(parts[[1L]])
  for (part in parts) {
    if (!is_categorical_encoded(part) || !same_domain(domain(part), domain)) {
      stop("Only categorical_encoded vectors in the same domain can be combined.", call. = FALSE)
    }
  }
  new_categorical_encoded(unlist(lapply(parts, categorical_codes)), domain)
}

#' @export
#' @noRd
rep.categorical_encoded <- function(x, ...) {
  new_categorical_encoded(categorical_codes(NextMethod()), domain(x))
}

#' @export
#' @noRd
unique.categorical_encoded <- function(x, incomparables = FALSE, ...) {
  new_categorical_encoded(unique(categorical_codes(x), incomparables = incomparables, ...), domain(x))
}

# NO `mtfrm()` METHOD, unlike `naxp_encoded`. Without one, match() and %in%
# compare text, so `x %in% "act"` works; with one they would compare codes
# against the text and quietly find nothing. Labels are unique within a domain,
# so matching by text is exact there.
# Sorting follows the domain's order where it has one, and the codes otherwise.
# Either way code 0 sorts first.
#' @export
#' @noRd
xtfrm.categorical_encoded <- function(x) {
  domain <- domain(x)
  codes <- categorical_codes(x)
  if (is_ordered_domain(domain)) category_ranks(domain)[codes + 1L] else codes
}

#' @export
#' @noRd
Ops.categorical_encoded <- function(e1, e2) {
  equality <- .Generic %in% c("==", "!=")
  ordering <- .Generic %in% c("<", "<=", ">", ">=")
  if (!equality && !ordering) {
    stop("`", .Generic, "` is not defined for categorical_encoded vectors, which support ",
         "only comparison.", call. = FALSE)
  }
  domain <- if (is_categorical_encoded(e1)) domain(e1) else domain(e2)
  if (ordering && !is_ordered_domain(domain)) {
    stop("`", .Generic, "` needs an ordered domain; this one is unordered.", call. = FALSE)
  }
  e1 <- comparison_codes(e1, domain)
  e2 <- comparison_codes(e2, domain)
  if (ordering) {
    # Compare ranks in the domain's order, in which code 0 is lowest.
    ranks <- category_ranks(domain)
    e1 <- ranks[e1 + 1L]
    e2 <- ranks[e2 + 1L]
  }
  get(.Generic, envir = baseenv())(e1, e2)
}

#' @export
#' @noRd
Math.categorical_encoded <- function(x, ...) {
  stop("`", .Generic, "()` is not defined for categorical_encoded vectors. Use as.integer() ",
       "for the codes.", call. = FALSE)
}

#' @export
#' @noRd
Summary.categorical_encoded <- function(..., na.rm = FALSE) {
  stop("`", .Generic, "()` is not defined for categorical_encoded vectors. Use as.integer() ",
       "for the codes.", call. = FALSE)
}

# Every category of the domain as a level, in the domain's order where it has
# one and in code order otherwise. Code 0 is not a category, so it becomes NA.
#' @exportS3Method forcats::as_factor
as_factor.categorical_encoded <- function(x, ...) {
  domain <- domain(x)
  ordered <- is_ordered_domain(domain)
  all_codes <- if (ordered) domain$sorted_codes else seq_len(domain$size)
  factor(
    decode.categorical_encoded(x),
    levels = category_text(domain, all_codes),
    ordered = ordered
  )
}
