# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# categorical -- logmu's own categorical vector
# ============================================================================
#
# NOT A FACTOR, AND NOT A SUBCLASS OF ONE. A factor cannot hold code 0: some
# base R and forcats functions refuse it as malformed and others shrink or
# miscount without warning. A subclass of factor is worse, because forcats
# keeps the subclass while reordering and merging the levels underneath it.
# So a categorical stores its own codes 0..size, exactly as the engine does, and
# converts to and from a factor at an explicit boundary where the rules are
# checked again. The pattern is naxp's `naxp_encoded`.
#
# NO NA, EVER. A categorical is int-backed and `int` has no portable missing value
# (see the NaN policy), so every operation that could produce an NA -- reading
# past the end, lengthening -- either fills code 0, where the domain allows it,
# or stops.
#
# The constructor is `as_categorical()`. `category()` is already the function that
# builds an indicator for one group of a categorical field.

categorical_class <- "categorical"

#' Categorical vectors
#'
#' @description
#' A `categorical` is a vector of codes in a [category_domain()]. Codes `1..size`
#' are the domain's categories and code 0 means undefined, which is allowed only
#' where the domain allows it. A categorical never holds `NA`.
#'
#' `as_categorical()` makes one:
#'
#' * from a **factor**, whose levels become the labels of a new domain, in
#'   their order, unless `domain` is given. Name the level that stands for
#'   code 0 with `undefined_level`. An ordered factor gives an ordered domain in
#'   the order of its levels, so the undefined level keeps its place in that
#'   order. `NA` is an error: give missing values a level of their own, such as
#'   `"unknown"`, and name it as `undefined_level`.
#' * from **text**, matched against the labels of `domain`, which is required.
#'   The domain says which categories exist; text cannot.
#' * from **whole numbers**, which are taken as codes in `domain`.
#'
#' Text that names no category of the domain is an error rather than becoming
#' undefined, so that a mistyped category cannot vanish.
#'
#' A categorical prints as its labels. It compares by equality with another in
#' the same domain, with text, or with numbers taken as codes. In an ordered
#' domain `<`, `<=`, `>` and `>=` compare by the domain's order, as do [sort()]
#' and [order()]. Arithmetic is an error. [c()] combines categorical vectors in
#' the same domain. Assigning into a categorical accepts another in the same
#' domain or text.
#'
#' `%in%` and [match()] compare labels, not codes, because R gives each side to
#' them separately. So `x %in% "act"` works, but numbers on the right are
#' compared with the labels, and a categorical from another domain is not
#' refused as it is by `==`.
#'
#' `factor(x)` gives a factor with the levels that occur, as it does for any
#' vector, in the domain's order. To convert with every category of the domain
#' as a level, use `forcats::as_factor(x)`: it gives the levels in the domain's
#' order, with the undefined label first where the domain is unordered and
#' allows undefined, and an ordered factor for an ordered domain.
#'
#' `domain()` returns the domain of a categorical, or the naxp of a
#' `naxp_encoded` vector, which logmu treats as a categorical whose domain is the
#' naxp.
#'
#' @param x An object.
#' @param domain A [category_domain()].
#' @param undefined_level The level of a factor that stands for code 0.
#' @param ... Unused.
#' @returns
#' `as_categorical()` returns a `categorical`.
#'
#' `is_categorical()` returns a scalar `logical`.
#'
#' `domain()` returns a `category_domain`, or a naxp.
#' @examples
#' status <- as_categorical(factor(c("act", "def", "act", "pen")))
#' status
#' domain(status)
#' status == "act"
#'
#' # Codes in a domain of socio-economic types, with 0 for none.
#' types <- category_domain(60, allows_undefined = TRUE)
#' as_categorical(c(12, 0, 60), types)
#'
#' # A factor level that stands for undefined.
#' as_categorical(factor(c("a", "none", "b")), undefined_level = "none")
#' @export
as_categorical <- function(x, ...) UseMethod("as_categorical")

#' @rdname as_categorical
#' @export
as_categorical.default <- function(x, ...) {
  stop("Cannot make a `categorical` from an object of class `", class(x)[[1L]], "`.",
       call. = FALSE)
}

#' @rdname as_categorical
#' @export
as_categorical.categorical <- function(x, domain = NULL, ...) {
  if (!is.null(domain)) {
    ensure_is_category_domain(domain)
    if (!same_category_domain(domain(x), domain)) {
      stop("This categorical is in a different domain.", call. = FALSE)
    }
  }
  x
}

#' @rdname as_categorical
#' @export
as_categorical.factor <- function(x, domain = NULL, undefined_level = NULL, ...) {
  if (anyNA(x)) {
    stop("A categorical cannot hold NA. Give the missing values a level of their own, ",
         "such as \"unknown\", and name it as `undefined_level`.", call. = FALSE)
  }
  if (!is.null(domain)) {
    if (!is.null(undefined_level)) {
      stop("Give `undefined_level` or `domain`, not both: a domain has its own ",
           "undefined label.", call. = FALSE)
    }
    return(as_categorical.character(as.character(x), domain))
  }
  if (!is.null(undefined_level)) ensure_is_single_valid_text(undefined_level)
  labels <- setdiff(levels(x), undefined_level)
  # An ordered factor's levels are its order, with the undefined level in its
  # own place. Its codes run 1, 2, ... over the other levels in the same order.
  level_codes <- match(levels(x), labels)
  level_codes[is.na(level_codes)] <- 0L
  if (is.ordered(x) && !is.null(undefined_level) && !undefined_level %in% levels(x)) {
    stop("The undefined level \"", undefined_level, "\" must be a level of an ordered ",
         "factor, so that its place in the order is known.", call. = FALSE)
  }
  domain <- category_domain(
    labels = labels,
    undefined_label = undefined_level,
    sorted_codes = if (is.ordered(x)) level_codes
  )
  new_categorical(level_codes[as.integer(x)], domain)
}

#' @rdname as_categorical
#' @export
as_categorical.character <- function(x, domain, ...) {
  if (missing(domain)) {
    stop("Text needs a `domain` to say which categories exist.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("A categorical cannot hold NA.", call. = FALSE)
  }
  ensure_is_category_domain(domain)
  new_categorical(category_domain_codes(domain, unname(x)), domain)
}

#' @rdname as_categorical
#' @export
as_categorical.numeric <- function(x, domain, ...) {
  if (missing(domain)) {
    stop("Codes need a `domain` to belong to.", call. = FALSE)
  }
  ensure_is_category_domain(domain)
  new_categorical(checked_category_codes(x, domain), domain)
}

# Whole numbers in 0..size, with 0 only where the domain allows it.
checked_category_codes <- function(codes, domain) {
  codes <- unclass(codes)
  attributes(codes) <- NULL
  if (anyNA(codes) || any(codes != round(codes)) ||
      any(codes < 0) || any(codes > domain$size)) {
    stop("Codes must be whole numbers from ", if (domain$allows_undefined) 0 else 1,
         " to ", domain$size, ", with no NA.", call. = FALSE)
  }
  codes <- as.integer(codes)
  ensure_undefined_allowed(codes, domain)
  codes
}

ensure_undefined_allowed <- function(codes, domain) {
  if (!domain$allows_undefined && any(codes == 0L)) {
    stop("This domain does not allow code 0 (undefined).", call. = FALSE)
  }
}

new_categorical <- function(codes, domain) {
  structure(codes, domain = domain, class = categorical_class)
}

#' @rdname as_categorical
#' @export
is_categorical <- function(x) inherits(x, categorical_class)

ensure_is_categorical <- function(x) {
  if (!is_categorical(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `categorical`.", call. = FALSE)
  }
}

#' @rdname as_categorical
#' @export
domain <- function(x) UseMethod("domain")

#' @rdname as_categorical
#' @export
domain.default <- function(x) {
  stop("Only a `categorical` or a `naxp_encoded` vector has a domain.", call. = FALSE)
}

#' @rdname as_categorical
#' @export
domain.categorical <- function(x) attr(x, "domain", exact = TRUE)

#' @rdname as_categorical
#' @export
domain.naxp_encoded <- function(x) naxp_of(x)

# The codes without class, domain or names.
categorical_codes <- function(x) {
  codes <- unclass(x)
  attributes(codes) <- NULL
  codes
}

# Codes to assign into, or compare with, a categorical in `domain`.
categorical_codes_under <- function(value, domain, numbers_allowed) {
  if (is_categorical(value)) {
    if (!same_category_domain(domain(value), domain)) {
      stop("These categorical vectors are in different domains.", call. = FALSE)
    }
    return(categorical_codes(value))
  }
  if (is.character(value)) {
    return(category_domain_codes(domain, value))
  }
  if (numbers_allowed && is.numeric(value) && !is.object(value)) {
    return(value)
  }
  stop(
    if (numbers_allowed) {
      "A categorical can be compared only with one in the same domain, with text, or with codes."
    } else {
      "Assign a categorical in the same domain, or text, which is matched to its labels."
    },
    call. = FALSE
  )
}

# Positions gained by lengthening hold NA. A categorical fills them with code 0
# where its domain allows it, and otherwise cannot be lengthened.
fill_undefined <- function(codes, domain) {
  if (anyNA(codes)) {
    if (!domain$allows_undefined) {
      stop("This would leave positions with no category, and the domain does not ",
           "allow code 0 (undefined).", call. = FALSE)
    }
    codes[is.na(codes)] <- 0L
  }
  codes
}

#' @export
#' @noRd
format.categorical <- function(x, ...) {
  category_domain_text(domain(x), categorical_codes(x))
}

#' @export
#' @noRd
print.categorical <- function(x, ...) {
  cat("<categorical[", length(x), "]> ", category_domain_label(domain(x)), "\n", sep = "")
  if (length(x) > 0L) {
    print(format(x), quote = FALSE)
  }
  invisible(x)
}

#' @export
#' @noRd
as.character.categorical <- function(x, ...) format.categorical(x)

#' @export
#' @noRd
as.data.frame.categorical <- function(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  nm = deparse1(substitute(x))) {
  as.data.frame.vector(x, row.names = row.names, optional = optional, ..., nm = nm)
}

#' @export
#' @noRd
`[.categorical` <- function(x, ...) {
  domain <- domain(x)
  new_categorical(fill_undefined(categorical_codes(NextMethod()), domain), domain)
}

#' @export
#' @noRd
`[[.categorical` <- function(x, ...) {
  new_categorical(categorical_codes(NextMethod()), domain(x))
}

#' @export
#' @noRd
`[<-.categorical` <- function(x, ..., value) {
  domain <- domain(x)
  value <- categorical_codes_under(value, domain, numbers_allowed = FALSE)
  codes <- categorical_codes(x)
  codes[...] <- value
  new_categorical(fill_undefined(codes, domain), domain)
}

#' @export
#' @noRd
`[[<-.categorical` <- function(x, ..., value) {
  domain <- domain(x)
  value <- categorical_codes_under(value, domain, numbers_allowed = FALSE)
  codes <- categorical_codes(x)
  codes[[...]] <- value
  new_categorical(fill_undefined(codes, domain), domain)
}

#' @export
#' @noRd
`length<-.categorical` <- function(x, value) {
  domain <- domain(x)
  codes <- categorical_codes(x)
  length(codes) <- value
  new_categorical(fill_undefined(codes, domain), domain)
}

#' @export
#' @noRd
c.categorical <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  domain <- domain(parts[[1L]])
  for (part in parts) {
    if (!is_categorical(part) || !same_category_domain(domain(part), domain)) {
      stop("Only categorical vectors in the same domain can be combined.", call. = FALSE)
    }
  }
  new_categorical(unlist(lapply(parts, categorical_codes)), domain)
}

#' @export
#' @noRd
rep.categorical <- function(x, ...) {
  new_categorical(categorical_codes(NextMethod()), domain(x))
}

#' @export
#' @noRd
unique.categorical <- function(x, incomparables = FALSE, ...) {
  new_categorical(unique(categorical_codes(x), incomparables = incomparables, ...), domain(x))
}

# NO `mtfrm()` METHOD, unlike `naxp_encoded`. Without one, match() and %in%
# compare labels, so `x %in% "act"` works; with one they would compare codes
# against the text and quietly find nothing. Labels are unique within a domain,
# so matching by label is exact there.
# Sorting follows the domain's order where it has one, and the codes otherwise.
#' @export
#' @noRd
xtfrm.categorical <- function(x) {
  domain <- domain(x)
  codes <- categorical_codes(x)
  if (is_ordered_domain(domain)) category_domain_ranks(domain)[codes + 1L] else codes
}

#' @export
#' @noRd
Ops.categorical <- function(e1, e2) {
  equality <- .Generic %in% c("==", "!=")
  ordering <- .Generic %in% c("<", "<=", ">", ">=")
  if (!equality && !ordering) {
    stop("`", .Generic, "` is not defined for categorical vectors, which support only comparison.",
         call. = FALSE)
  }
  domain <- if (is_categorical(e1)) domain(e1) else domain(e2)
  if (ordering && !is_ordered_domain(domain)) {
    stop("`", .Generic, "` needs an ordered domain; this one is unordered.", call. = FALSE)
  }
  e1 <- categorical_codes_under(e1, domain, numbers_allowed = TRUE)
  e2 <- categorical_codes_under(e2, domain, numbers_allowed = TRUE)
  if (ordering) {
    # Compare ranks in the domain's order, in which code 0 has a place like any
    # other code. Bare numbers must be codes of the domain to have a rank.
    ranks <- category_domain_ranks(domain)
    e1 <- ranks[checked_category_codes(e1, domain) + 1L]
    e2 <- ranks[checked_category_codes(e2, domain) + 1L]
  }
  get(.Generic, envir = baseenv())(e1, e2)
}

#' @export
#' @noRd
Math.categorical <- function(x, ...) {
  stop("`", .Generic, "()` is not defined for categorical vectors. Use as.integer() for the codes.",
       call. = FALSE)
}

#' @export
#' @noRd
Summary.categorical <- function(..., na.rm = FALSE) {
  stop("`", .Generic, "()` is not defined for categorical vectors. Use as.integer() for the codes.",
       call. = FALSE)
}

# Every category of the domain as a level: in the domain's order where it has
# one, and otherwise in code order with the undefined label first.
#' @exportS3Method forcats::as_factor
as_factor.categorical <- function(x, ...) {
  domain <- domain(x)
  ordered <- is_ordered_domain(domain)
  all_codes <- if (ordered) {
    domain$sorted_codes
  } else {
    c(if (domain$allows_undefined) 0L, seq_len(domain$size))
  }
  factor(
    format.categorical(x),
    levels = category_domain_text(domain, all_codes),
    ordered = ordered
  )
}
