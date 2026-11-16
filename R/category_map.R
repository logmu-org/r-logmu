# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# category_map -- maps from a naxp or a categorical to a value
# ============================================================================
#
# ONE CLASS. Its source and its value type are both attributes, because both
# are equally part of what a map is. The source decides only the storage,
# which is an internal matter:
#
#   naxp source         sparse: sorted naxp codes in `keys`, `values` parallel,
#                       and `undefined` for code 0 and every unlisted key
#   categorical source  dense: `values` indexed by code, one for every code of
#                       the source domain (size, or size + 1 where it allows
#                       code 0); unlisted codes and code 0 take `undefined` when
#                       the map is built, so nothing else needs storing
#
# The value type is logical, integer, double, naxp, categorical, or a flat
# struct of those (a data frame of columns in R).
#
# `undefined` DEFAULTS TO ALL-ZERO BITS: FALSE, 0L, +0.0, naxp code 0 or
# category code 0, field by field for a struct. It covers code 0 and every key
# the map does not list. Code 0 is therefore never a key.
#
# WHY `undefined` IS STORED. Without it, "unlisted keys get zero" is not closed
# under composition: if `f` sends code 0 to something other than zero, the
# single map `f` after `g` would need its unlisted keys to get that value, with
# nowhere to keep it. With it, composition always gives a single map.
#
# A MAP IS A FUNCTION, as `ecdf()` and `approxfun()` return functions: `oac(x)`
# looks `x` up. Its data are attributes, so the closure itself is one shared
# function whose environment is the namespace, and it finds its own attributes
# through `sys.function()`. Applying a map to another map is an error:
# composition is a deliberate act, through `compose_maps()`.
#
# Struct fields are rewritten at compile time in the engine (`coords(...)$easting`
# becomes a map onto the `easting` column alone), so a struct is never a value
# the engine carries. That is engine work and is not here yet.

category_map_class <- c("category_map", "function")

# Struct field names work directly as names in R, Python and C#. They are
# compared case-sensitively.
struct_field_pattern <- "^[A-Za-z][A-Za-z0-9_]*$"

#' Maps from naxps and categoricals
#'
#' @description
#' A `category_map` gives a value for every code of a naxp, such as a
#' postcode, or of a [categorical][as_categorical()]. It is a function: a map
#' `oac` applied to postcodes is `oac(postcodes)`.
#'
#' `keys` lists the codes the map gives values for, as a `naxp_encoded` vector
#' or a `categorical`; the map's source is its naxp or its domain. `values`
#' gives the value for each key, in the same order. It may be a logical,
#' integer or double vector, a `naxp_encoded` vector, a `categorical`, or a
#' data frame of those to give a struct of named fields. Logical and integer
#' values may not be `NA`. Field names must start with a letter and use only
#' ASCII letters, digits and underscores; they are case-sensitive.
#'
#' `undefined` is the value for code 0 and for every key the map does not
#' list. It defaults to zero: `FALSE`, `0`, naxp code 0 or category code 0,
#' and zero in every field of a struct. For a postcode, being listed is what
#' makes it valid, so code 0 is never itself a key. Where the default would be
#' category code 0 and the values' domain does not allow it, give `undefined`
#' explicitly.
#'
#' Note that a double map used as a weight gives a weight of 0 by default to a
#' key it does not list, which removes that record from the experience
#' entirely.
#'
#' Where both the source and a categorical value domain are ordered, the map
#' must keep the order: going up the source's order, the values may not go
#' down the values' order.
#'
#' Applying a map to a vector of its source gives a vector of its values, or a
#' data frame for a struct. A map keyed by a categorical also accepts text or a
#' factor, which is matched to its source's labels. `f(g(x))` applies `g` and
#' then `f`. To make a single map that does both, use `compose_maps(f, g)`,
#' which needs the values of `g` to be in the source of `f`. The result keeps
#' nothing of `f` or `g`.
#'
#' @param keys The keys: a `naxp_encoded` vector or a `categorical`. They must
#'   be unique and may not include code 0.
#' @param values The value for each key.
#' @param undefined The value for code 0 and every unlisted key, of the same
#'   type as `values`: a single value, or a one-row data frame for a struct.
#'   `NULL` for zero.
#' @param name An optional name for the map, read back by [name()].
#' @param f,g Maps; `compose_maps(f, g)` applies `g` first.
#' @param x An object.
#' @returns
#' `category_map()` and `compose_maps()` return a `category_map`, which is a
#' function.
#'
#' `is_category_map()` returns a scalar `logical`.
#' @examples
#' # Socio-economic types 1 to 3, remapped into two clusters.
#' types <- category_domain(3, allows_undefined = TRUE, name = "types")
#' clusters <- category_domain(labels = c("low", "high"), undefined_label = "none")
#' to_cluster <- category_map(
#'   keys = as_categorical(c(1, 2, 3), types),
#'   values = as_categorical(c("low", "low", "high"), clusters),
#'   name = "to_cluster"
#' )
#' to_cluster
#' to_cluster(as_categorical(c(3, 0, 1), types))
#'
#' # A weight for each cluster, with 1 where none applies.
#' weight <- category_map(
#'   keys = as_categorical(c("low", "high"), clusters),
#'   values = c(1.1, 0.9),
#'   undefined = 1
#' )
#' weight(to_cluster(as_categorical(c(3, 0, 1), types)))
#'
#' # The same as a single map.
#' type_weight <- compose_maps(weight, to_cluster)
#' type_weight(as_categorical(c(3, 0, 1), types))
#'
#' # A struct of two fields.
#' details <- category_map(
#'   keys = as_categorical(c(1, 3), types),
#'   values = data.frame(weight = c(0.9, 1.2), flagged = c(FALSE, TRUE))
#' )
#' details(as_categorical(c(1, 2, 3), types))
#' @export
category_map <- function(keys, values, undefined = NULL, name = NULL) {
  if (is_naxp_encoded(keys)) {
    source <- naxp_of(keys)
    key_codes <- naxp_codes(keys)
    key_text <- paste0("code ", key_codes)
  } else if (is_categorical(keys)) {
    source <- domain(keys)
    key_codes <- categorical_codes(keys)
    key_text <- format.categorical(keys)
  } else {
    stop("The keys of a `category_map` must be a `naxp_encoded` vector or a `categorical`.",
         call. = FALSE)
  }
  if (any(key_codes == 0)) {
    stop("Code 0 cannot be a key: give its value as `undefined`.", call. = FALSE)
  }
  repeated <- unique(key_text[duplicated(key_codes)])
  if (length(repeated) > 0L) {
    shown <- utils::head(repeated, 5L)
    stop("Keys must be unique; repeated: ", paste0("\"", shown, "\"", collapse = ", "),
         if (length(repeated) > length(shown)) ", ..." else "", ".", call. = FALSE)
  }
  if (!is.null(name)) ensure_is_valid_name(name)

  checked <- checked_map_values(values, length(key_codes))
  # A naxp has far too many codes to list, so its undefined is always reached.
  undefined_reached <- is_naxp(source) || source$allows_undefined || length(key_codes) < source$size
  undefined <- checked_undefined(undefined, checked, undefined_reached)

  if (is_naxp(source)) {
    order <- order(key_codes, method = "radix")
    return(new_category_map(
      source = source,
      keys = key_codes[order],
      values = take_map_values(checked$values, order),
      undefined = undefined,
      value_type = checked$value_type,
      name = name
    ))
  }

  # Dense storage: position code + 1 where the source allows code 0, code
  # otherwise. Unlisted codes and code 0 index nothing and so take `undefined`.
  offset <- as.integer(source$allows_undefined)
  index <- rep(NA_integer_, source$size + offset)
  index[key_codes + offset] <- seq_along(key_codes)
  dense_values <- take_map_values(checked$values, index, undefined)
  ensure_order_is_kept(source, dense_values)
  new_category_map(
    source = source,
    keys = NULL,
    values = dense_values,
    undefined = NULL,
    value_type = checked$value_type,
    name = name
  )
}

# The one closure every map shares. Its environment is the namespace, so a map
# carries nothing but its attributes.
category_map_function <- function(x) apply_category_map(sys.function(), x)

new_category_map <- function(source, keys, values, undefined, value_type, name) {
  structure(
    category_map_function,
    class = category_map_class,
    source = source,
    keys = keys,
    values = values,
    undefined = undefined,
    value_type = value_type,
    name = name
  )
}

#' @rdname category_map
#' @export
is_category_map <- function(x) inherits(x, "category_map")

ensure_is_category_map <- function(x) {
  if (!is_category_map(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `category_map`.", call. = FALSE)
  }
}

map_source <- function(map) attr(map, "source", exact = TRUE)
map_keys <- function(map) attr(map, "keys", exact = TRUE)
map_values <- function(map) attr(map, "values", exact = TRUE)
map_undefined <- function(map) attr(map, "undefined", exact = TRUE)
map_value_type <- function(map) attr(map, "value_type", exact = TRUE)

# ---- value types --------------------------------------------------------------

# The value type of `values`, and `values` themselves stripped of names: a
# list with `values` and `value_type`.
checked_map_values <- function(values, key_count, what = "`values`") {
  if (is.data.frame(values)) {
    ensure_struct_fields(values)
    if (nrow(values) != key_count) {
      stop(what, " must have ", key_count, if (key_count == 1L) " row." else " rows.",
           call. = FALSE)
    }
    fields <- lapply(names(values), function(field) {
      checked_column_values(values[[field]], paste0("Field `", field, "`"))
    })
    names(fields) <- names(values)
    return(list(values = new_struct(fields, key_count), value_type = "struct"))
  }
  value_type <- column_value_type(values, what)
  if (length(values) != key_count) {
    stop(what, " must have ", key_count, if (key_count == 1L) " value." else " values.",
         call. = FALSE)
  }
  list(values = checked_column_values(values, what), value_type = value_type)
}

ensure_struct_fields <- function(values) {
  fields <- names(values)
  if (length(fields) == 0L) {
    stop("A struct of values must have at least one field.", call. = FALSE)
  }
  bad <- fields[!grepl(struct_field_pattern, fields)]
  if (length(bad) > 0L) {
    stop("Field names must start with a letter and use only ASCII letters, digits and ",
         "underscores; got ", paste0("\"", bad, "\"", collapse = ", "), ".", call. = FALSE)
  }
  repeated <- unique(fields[duplicated(fields)])
  if (length(repeated) > 0L) {
    stop("Field names must be unique; repeated: ",
         paste0("\"", repeated, "\"", collapse = ", "), ".", call. = FALSE)
  }
}

# One of "logical", "integer", "double", "naxp" or "categorical". A struct may
# not nest, so a data frame or a list is not a column type.
column_value_type <- function(column, what) {
  if (is_categorical(column)) return("categorical")
  if (is_naxp_encoded(column)) return("naxp")
  if (is.factor(column)) {
    stop(what, " is a factor. Convert it with `as_categorical()` first.", call. = FALSE)
  }
  if (!is.object(column) && is.null(dim(column))) {
    if (is.logical(column)) return("logical")
    if (is.integer(column)) return("integer")
    if (is.double(column)) return("double")
  }
  stop(what, " must be logical, integer, double, a `naxp_encoded` vector or a ",
       "`categorical`. Structs may not nest.", call. = FALSE)
}

checked_column_values <- function(column, what) {
  value_type <- column_value_type(column, what)
  if (value_type %in% c("logical", "integer") && anyNA(column)) {
    stop(what, " may not be NA: a ", value_type, " has no portable missing value.",
         call. = FALSE)
  }
  if (value_type %in% c("logical", "integer", "double")) {
    column <- as.vector(column)
    names(column) <- NULL
  }
  column
}

new_struct <- function(fields, row_count) {
  structure(fields, class = "data.frame", row.names = c(NA_integer_, -row_count))
}

# Whether two columns have the same value type, in the same domain or naxp.
same_column_type <- function(a, b) {
  type <- column_value_type(a, "A value")
  if (type != column_value_type(b, "A value")) return(FALSE)
  switch(type,
    categorical = same_category_domain(domain(a), domain(b)),
    naxp = same_naxp(naxp_of(a), naxp_of(b)),
    TRUE
  )
}

# `undefined` as a single value of the map's value type: zero where it is not
# given. Zero must be a valid value where the map can reach it.
checked_undefined <- function(undefined, checked, undefined_reached) {
  if (is.null(undefined)) {
    zero <- zero_value_like(checked$values)
    if (undefined_reached) ensure_zero_is_valid(zero)
    return(zero)
  }
  given <- checked_map_values(undefined, 1L, "`undefined`")
  same <- given$value_type == checked$value_type && if (checked$value_type == "struct") {
    identical(names(given$values), names(checked$values)) &&
      all(mapply(same_column_type, given$values, checked$values))
  } else {
    same_column_type(given$values, checked$values)
  }
  if (!same) {
    stop("`undefined` must be of the same type as `values`, with the same fields, ",
         "domains and naxps.", call. = FALSE)
  }
  given$values
}

# A single all-zero value of the type of `values`. Category code 0 is made here
# even where the domain does not allow it; `ensure_zero_is_valid()` says so.
zero_value_like <- function(values) {
  if (is.data.frame(values)) {
    return(new_struct(lapply(values, zero_value_like), 1L))
  }
  if (is_categorical(values)) return(new_categorical(0L, domain(values)))
  if (is_naxp_encoded(values)) {
    return(new_naxp_encoded_codes(if (is.integer(naxp_codes(values))) 0L else 0, naxp_of(values)))
  }
  if (is.logical(values)) FALSE else if (is.integer(values)) 0L else 0
}

ensure_zero_is_valid <- function(zero) {
  columns <- if (is.data.frame(zero)) zero else list(zero)
  for (column in columns) {
    if (is_categorical(column) && !domain(column)$allows_undefined) {
      stop("Some keys are unlisted or code 0, and get `undefined`. It defaults to ",
           "category code 0, which the values' domain does not allow, so give ",
           "`undefined` explicitly.", call. = FALSE)
    }
  }
}

# ---- reading values by position -------------------------------------------------

# `values[index]`, where an NA index -- code 0 or an unlisted key -- gives `fill`.
take_map_values <- function(values, index, fill = NULL) {
  if (is.data.frame(values)) {
    fields <- lapply(names(values), function(field) {
      take_column_values(values[[field]], index, fill[[field]])
    })
    names(fields) <- names(values)
    return(new_struct(fields, length(index)))
  }
  take_column_values(values, index, fill)
}

take_column_values <- function(column, index, fill) {
  missing <- is.na(index)
  if (is_categorical(column)) {
    codes <- categorical_codes(column)[index]
    if (any(missing)) codes[missing] <- categorical_codes(fill)
    return(new_categorical(codes, domain(column)))
  }
  if (is_naxp_encoded(column)) {
    codes <- naxp_codes(column)[index]
    if (any(missing)) codes[missing] <- naxp_codes(fill)
    return(new_naxp_encoded_codes(codes, naxp_of(column)))
  }
  taken <- column[index]
  if (any(missing)) taken[missing] <- fill
  taken
}

# Where both domains are ordered, going up the source's order the values may
# not go down the values' order.
ensure_order_is_kept <- function(source, dense_values) {
  if (!is_ordered_domain(source) || !is_categorical(dense_values) ||
      !is_ordered_domain(domain(dense_values))) {
    return(invisible())
  }
  offset <- as.integer(source$allows_undefined)
  value_codes <- categorical_codes(dense_values)[source$sorted_codes + offset]
  value_ranks <- category_domain_ranks(domain(dense_values))[value_codes + 1L]
  if (is.unsorted(value_ranks)) {
    stop("Both domains are ordered, so the map must keep the order: going up the ",
         "source's order, the values may not go down.", call. = FALSE)
  }
}

# ---- applying and composing maps ----------------------------------------------------

apply_category_map <- function(map, x) {
  if (is_category_map(x)) {
    stop("Applying a map to a map does not compose them. Use `compose_maps()`.",
         call. = FALSE)
  }
  source <- map_source(map)
  if (is_naxp(source)) {
    if (!is_naxp_encoded(x) || !same_naxp(naxp_of(x), source)) {
      stop("This map looks up `naxp_encoded` codes of the naxp ", naxp_pattern(source), ".",
           call. = FALSE)
    }
    # Code 0 is never a key, so it finds nothing and gets `undefined` too.
    index <- match(naxp_codes(x), map_keys(map))
    return(take_map_values(map_values(map), index, map_undefined(map)))
  }
  if (is.character(x) || is.factor(x)) {
    x <- as_categorical(as.character(x), source)
  }
  if (!is_categorical(x) || !same_category_domain(domain(x), source)) {
    stop("This map looks up categories in the domain ", category_domain_label(source), ".",
         call. = FALSE)
  }
  take_map_values(map_values(map), categorical_codes(x) + as.integer(source$allows_undefined))
}

#' @rdname category_map
#' @export
compose_maps <- function(f, g, name = NULL) {
  ensure_is_category_map(f)
  ensure_is_category_map(g)
  g_values <- map_values(g)
  if (!map_value_type(g) %in% c("categorical", "naxp") ||
      !same_category_domain(domain(g_values), map_source(f))) {
    stop("`compose_maps(f, g)` needs the values of `g` to be in the source of `f`.",
         call. = FALSE)
  }
  if (!is.null(name)) ensure_is_valid_name(name)
  # `g`'s unlisted keys get its `undefined`, which `f` then maps like any other
  # value, so the single map's `undefined` is `f` of `g`'s.
  new_category_map(
    source = map_source(g),
    keys = map_keys(g),
    values = f(g_values),
    undefined = if (is_naxp(map_source(g))) f(map_undefined(g)),
    value_type = map_value_type(f),
    name = name
  )
}

# ---- printing -------------------------------------------------------------------

#' @export
#' @noRd
print.category_map <- function(x, ...) {
  source <- map_source(x)
  count <- if (is_naxp(source)) length(map_keys(x)) else NROW(map_values(x))
  cat(sprintf(
    "<category_map%s> %d %s from %s to %s\n",
    if (is.null(name(x))) "" else paste0(": ", name(x)),
    count,
    if (is_naxp(source)) {
      if (count == 1L) "key" else "keys"
    } else {
      if (count == 1L) "code" else "codes"
    },
    category_domain_label(source),
    describe_value_type(map_values(x), map_value_type(x))
  ))
  invisible(x)
}

describe_value_type <- function(values, value_type) {
  switch(value_type,
    categorical = paste0("categorical (", category_domain_label(domain(values)), ")"),
    naxp = paste0("naxp ", naxp_pattern(naxp_of(values))),
    struct = paste0("struct (", paste(names(values), collapse = ", "), ")"),
    value_type
  )
}
