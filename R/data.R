# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

###### Typing ######

exp_data_class <- c("exp_data", "logmu_data", "tbl", "data.frame")
val_data_class <- c("val_data", "logmu_data", "tbl", "data.frame")

#' Is `x` an `exp_data` or a `val_data`?
#'
#' @description
#'
#' `is_exp_data()` checks whether an object is a `exp_data`.
#'
#' `is_val_data()` checks whether an object is a `val_data`.
#'
#' @param x The object to test.
#' @returns A scalar `logical`.
#' @name is_exp_and_val_data_type
NULL

#' @rdname is_exp_and_val_data_type
#' @export
is_exp_data <- function(x) {
  typeof(x) == "list" && isa(x, exp_data_class)
}

#' @rdname is_exp_and_val_data_type
#' @export
is_val_data <- function(x) {
  typeof(x) == "list" && isa(x, val_data_class)
}

ensure_is_exp_data <- function(x) {
  if (!is_exp_data(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be `exp_data`.", call. = FALSE)
  }
}
ensure_is_val_data <- function(x) {
  if (!is_val_data(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be `val_data`.", call. = FALSE)
  }
}

###### Attributes ######

#' Get experience period of `exp_data`
#'
#' @description
#'
#' `exp_start` gets the start of the overall experience period (inclusive).
#'
#' `exp_end` gets the end of the overall experience period (exclusive).
#'
#' @param exp_data The mortality experience data.
#' @returns A scalar `datey`.
#' @name exp_data_info
NULL

#' @rdname exp_data_info
#' @export
exp_start <- function(exp_data) {
  ensure_is_exp_data(exp_data)
  attr(exp_data, "exp_start", exact = TRUE)
}

#' @rdname exp_data_info
#' @export
exp_end <- function(exp_data) {
  ensure_is_exp_data(exp_data)
  attr(exp_data, "exp_end", exact = TRUE)
}

#' Get the 'as at date' of a `val_data`
#'
#' @description
#'
#' Gets the 'as at' date of a valuation data.
#'
#' @param val_data The valuation data.
#' @returns A scalar `datey`.
#' @export
as_at <- function(val_data) {
  ensure_is_val_data(val_data)
  attr(val_data, "as_at", exact = TRUE)
}

###### Ctors ######

#' Create an `exp_data` or a `val_data`
#'
#' @description
#' Create a experience or valuation data from a list of columns, a `data.frame` or
#' `tibble` (or similar).
#'
#' The columns must be atomic, the same length and uniquely named.
#'
#' There *must* be a `birth` date column.
#'
#' A `count` column is optional. If present it must be numeric, non-negative
#' and finite. Rows with a nil count are dropped.
#'
#' Experience data must also contain (and valuation data must *not* contain)
#' - `E2R_start` and `E2R_end` date columns, and
#' - an `E2R_died` logical column.
#'
#' These special columns cannot contain NA and, if they are dates, they must be
#' `datey`.
#'
#' The experience data period and individual E2Rs must be *proper*:
#' - at the dataset level,  `exp_start` cannot be after `exp_end`, and
#' - at the record level, `E2R_start` cannot be after `E2R_end` and if they are
#' the same then `E2R_died` must be false.
#'
#' The experience data is trimmed to the experience period, i.e.
#' it is guaranteed that
#' `E2R_start` is not before `exp_start` and `E2R_end` is not after `exp_end`.
#' Records with empty E2Rs after trimming are dropped.
#'
#' @param columns The list of columns (which includes things like `data.frame`
#' and `tibble`) to convert to an `exp_data` or `val_data`.
#' @param exp_start The start of the overall experience period (inclusive).
#' Must be before `exp_end`.
#' @param exp_end The end of the overall experience period (exclusive).
#' Must be after `exp_start`.
#' @param as_at The 'as at' date of the valuation data.
#' @param on_unreadable What to do about columns whose type **logmu** cannot
#' read (anything other than `double`, `integer`, `logical`, `character`,
#' `factor`, `datey` or `durationy`): `"error"` (the default) rejects them,
#' `"warn"` keeps them with a warning, and `"ignore"` keeps them silently. The
#' special columns (`birth`, the `E2R_*` columns, `count`) are always checked.
#' @returns An `exp_data` or `val_data`.
#' @name exp_and_val_data
NULL

#' @rdname exp_and_val_data
#' @export
exp_data <- function(columns, exp_start, exp_end,
                     on_unreadable = c("error", "warn", "ignore")) {
  on_unreadable <- match.arg(on_unreadable)
  new_exp_data(columns, exp_start, exp_end, on_unreadable = on_unreadable)
}

#' @rdname exp_and_val_data
#' @export
val_data <- function(columns, as_at, on_unreadable = c("error", "warn", "ignore")) {
  on_unreadable <- match.arg(on_unreadable)
  new_val_data(columns, as_at, on_unreadable = on_unreadable)
}

#' Create a `val_data` from an `exp_data`
#'
#' @description
#'
#' Creates a `val_data` from those individuals  in an `exp_data` who
#' were alive at the specified date.
#'
#' @param exp_data The `exp_data`.
#' @param as_at The as at date (a `datey`). If omitted, the *end* of the
#' experience period is used.
#' @returns A `val_data`.
#' @export
val_from_exp_data <- function(exp_data, as_at = NULL) {

  # These ensure that exp_data really is an exp_data
  start <- exp_start(exp_data)
  end <- exp_end(exp_data)

  if (is.null(as_at)) as_at <- end

  ensure_is_single_valid_datey(as_at)

  columns <- unclass(exp_data)

  # Get E2R columns
  start <- columns[["E2R_start"]]
  end <- columns[["E2R_end"]]
  died <- columns[["E2R_died"]]

  # and drop them
  columns[["E2R_start"]] <- NULL
  columns[["E2R_end"]] <- NULL
  columns[["E2R_died"]] <- NULL

  # mask = t in [a,b) or t == b and alive
  keep <- (start <= as_at & as_at < end) | (as_at == end & !died)
  columns <- lapply(columns, function(x) {
    result <- x[keep]
    attributes(result) <- attributes(x)
    result
  })

  # NB Setting new_col_name to "" stops it re-checking count
  new_val_data(columns, as_at, new_col_name = "")
}

###### Column assignment ######

#' Add or assign a data column
#'
#' @description
#'
#' Add or assign a data column to `exp_data` or `val_data`.
#'
#' @param x The `exp_data` or `val_data`.
#' @param name The `name` of the column to assign.
#' @param value The column values.
#' @returns The updated `exp_data` or `val_data`.
#' @name assign_column
NULL

#' @rdname assign_column
#' @export
`$<-.exp_data` <- function(x, name, value) {
  columns <- unclass(x)
  columns[[name]] <- value
  new_exp_data(columns, exp_start = exp_start(x), exp_end = exp_end(x), new_col_name = name)
}

#' @rdname assign_column
#' @export
`$<-.val_data` <- function(x, name, value) {
  columns <- unclass(x)
  columns[[name]] <- value
  new_val_data(columns, as_at = as_at(x), new_col_name = name)
}

###### Internal ctors ######

new_exp_data <- function(columns, exp_start, exp_end, new_col_name = NULL,
                         on_unreadable = "error") {

  columns <- process_data_columns(
    columns,
    is_exp = TRUE,
    new_col_name = new_col_name,
    on_unreadable = on_unreadable
    )

  if (!is_single_valid_datey(exp_start))
    stop("`exp_start` must be a single valid `datey`.", call. = FALSE)

  if (!is_single_valid_datey(exp_end))
    stop("`exp_end` must be a single valid `datey`.", call. = FALSE)

  if (exp_start >= exp_end)
    stop("There must be a non-zero time period between `exp_start` and `exp_end`.", call. = FALSE)

  if (is.null(new_col_name) || new_col_name %in% c("E2R_start","E2R_end","E2R_died")) {

    E2R_start <- columns[["E2R_start"]]
    E2R_end <- columns[["E2R_end"]]
    E2R_died <- columns[["E2R_died"]]

    if (any(E2R_start > E2R_end))
      stop("Invalid E2Rs -- it is an error if `E2R_start` is after `E2R_end`.", call. = FALSE)
    if (any(E2R_start == E2R_end & E2R_died))
      stop("Invalid E2Rs -- it is an error if `E2R_start` is the same as `E2R_end` and `E2R_died` is true.", call. = FALSE)

    new_E2R_start <- ifelse(E2R_start >= exp_start, E2R_start, exp_start)
    keep_end <- E2R_end <= exp_end # <= not < because we want to retain deaths at end of E2R
    new_E2R_end <- ifelse(keep_end, E2R_end, exp_end)
    new_E2R_died <- ifelse(keep_end, E2R_died, FALSE)

    attributes(new_E2R_start) <- attributes(E2R_start)
    attributes(new_E2R_end) <- attributes(E2R_end)
    attributes(new_E2R_died) <- attributes(E2R_died)

    columns[["E2R_start"]] <- new_E2R_start
    columns[["E2R_end"]] <- new_E2R_end
    columns[["E2R_died"]] <- new_E2R_died

    # Remove E2Rs that have been made empty (or improper) by the above trimming
    has_E2R <- new_E2R_start < new_E2R_end
    if (!all(has_E2R)) {
      columns <- lapply(columns, function(x) {
        result <- x[has_E2R]
        attributes(result) <- attributes(x)
        result
      })
    }

  }

  structure(columns, class = exp_data_class, exp_start = exp_start, exp_end = exp_end, row.names = c(NA_integer_, -length(columns$birth)))
}

new_val_data <- function(columns, as_at, new_col_name = NULL,
                         on_unreadable = "error") {
  columns <- process_data_columns(
    columns,
    is_exp = FALSE,
    new_col_name = new_col_name,
    on_unreadable = on_unreadable
    )

  if (!is_single_valid_datey(as_at))
    stop("`as_at` must be a single valid `datey`.", call. = FALSE)

  structure(columns, class = val_data_class, as_at = as_at, row.names = c(NA_integer_, -length(columns$birth)))
}


###### Helper methods ######

process_data_columns <- function(columns, is_exp, new_col_name, on_unreadable = "error") {

  if (typeof(columns) != "list")
    stop("`columns` is not a list.", call. = FALSE)

  col_names <- names(columns)
  if (is.null(col_names))
    stop("The columns must be named.", call. = FALSE)

  # Clear all columns attributes (including class) other than names
  attributes(columns) <- NULL
  names(columns) <- col_names

  index_first_duplicate <- anyDuplicated(col_names)
  if(index_first_duplicate > 0L)
    stop("Duplicate column name `", names[index_first_duplicate], "`", call. = FALSE)

  seen_birth <- FALSE
  seen_E2R_start <- FALSE
  seen_E2R_end <- FALSE
  seen_E2R_died <- FALSE
  seen_count <- FALSE

  common_len <- NULL

  special_names_LC <- c("birth", "e2r_start", "e2r_end", "e2r_died", "count", "uid")
  illegal_names_LC <- "e2r"

  for (i in seq_along(columns)) {

    col_name <- col_names[i]
    if (is.null(col_name) || is.na(col_name) || col_name == "")
      stop("Missing name for column ", i, ".", call. = FALSE)
    if (grepl(".", col_name, fixed = TRUE))
      stop("`", col_name, "` cannot be used as a column name because it contains a full stop ('.').", call. = FALSE)

    if (col_name == "birth") seen_birth <- TRUE
    else if (col_name == "E2R_start") seen_E2R_start <- TRUE
    else if (col_name == "E2R_end") seen_E2R_end <- TRUE
    else if (col_name == "E2R_died") seen_E2R_died <- TRUE
    else if (col_name == "count") seen_count <- TRUE
    else if (col_name != "uid"){
      lower_col_name <- tolower(col_name)
      if (lower_col_name %in% special_names_LC)
        stop("`", col_name, "` is a special name but the case is wrong.", call. = FALSE)
      if (lower_col_name %in% illegal_names_LC)
        stop("`", col_name, "` cannot be used as a column name.", call. = FALSE)
    }

    column <- columns[[i]]

    if (!is.atomic(column))
      stop("Column name `", col_name, "` is not an atomic vector.", call. = FALSE)

    len <- length(column)
    if (is.null(common_len)) {
      common_len <- len
    } else if (len != common_len) {
      stop("Inconsistent column lengths -- column `", col_name,
        "` has ", len, " rows, but an earlier column has ",
        common_len, " rows,", call. = FALSE)
    }
  }

  if (!seen_birth) stop("Missing birth column.", call. = FALSE)

  if (is.null(new_col_name) || new_col_name == "birth")
    check_column_is_valid(columns, "birth")

  if (is_exp) {

    if (!seen_E2R_start || !seen_E2R_end || !seen_E2R_died)
      stop("Exp data must contain columns named `E2R_start`, `E2R_end` and `E2R_died`.", call. = FALSE)

    if (is.null(new_col_name) || new_col_name == "E2R_start")
      check_column_is_valid(columns, "E2R_start")

    if (is.null(new_col_name) || new_col_name == "E2R_end")
      check_column_is_valid(columns, "E2R_end")

    if (is.null(new_col_name) || new_col_name == "E2R_died")
      check_column_is_valid(columns, "E2R_died")

  } else {

    if (seen_E2R_start || seen_E2R_end || seen_E2R_died)
      stop("Val data cannot contain columns named `E2R_start`, `E2R_end` or `E2R_died`.", call. = FALSE)

  }

  # Check that ordinary columns are of a type logmu can read. The special
  # columns are validated separately above, so they are exempt here.
  check_columns_are_readable(columns, new_col_name, on_unreadable)

  # Drop records with count == 0
  if (seen_count && (is.null(new_col_name) || new_col_name == "count")) {
    check_column_is_valid(columns, "count")
    non_zero_count <- columns$count > 0
    if (!all(non_zero_count)) {
      columns <- lapply(columns, function(x) {
        result <- x[non_zero_count]
        attributes(result) <- attributes(x)
        result
      })
    }
  }

  columns
}

check_column_is_valid <- function(columns, col_name) {
  column <- columns[[col_name]]
  # All reserved columns are `datey` except `E2R_died` and `count`
  if (col_name == "E2R_died") {
    if (!is_pure_logical(column))
      stop("Column `E2R_died` must be a `logical`.", call. = FALSE)

  } else if (col_name == "count") {
    if (!is_pure_numeric(column))
      stop("Column `count` must be a numeric.", call. = FALSE)
    if (anyNA(column) || !all(is.finite(column) & column >= 0))
      stop("A `count` column must be a non-negative and finite.", call. = FALSE)
  } else {
    if (!datey::is_datey(column))
      stop("Column `", col_name, "` must be a `datey`.", call. = FALSE)
  }

  if (anyNA(column))
    stop("Column `", col_name, "` cannot contain NAs.", call. = FALSE)
}

# The column types logmu can read.
logmu_readable_types <- c("double", "integer", "logical", "character",
                          "factor", "datey", "durationy")

is_logmu_readable <- function(column) {
  is_pure_numeric(column) ||   # raw double or integer
    is_pure_logical(column) ||
    is_pure_text(column) ||
    is.factor(column) ||
    datey::is_datey(column) ||
    datey::is_durationy(column)
}

# Apply the `on_unreadable` policy to any ordinary column logmu cannot read.
# The special columns (birth, E2R_*, count) are validated elsewhere and exempt.
check_columns_are_readable <- function(columns, new_col_name, on_unreadable) {

  special <- c("birth", "E2R_start", "E2R_end", "E2R_died", "count")

  to_check <- if (is.null(new_col_name)) {
    setdiff(names(columns), special)
  } else if (new_col_name == "" || new_col_name %in% special) {
    character(0L)
  } else {
    new_col_name
  }

  if (length(to_check) == 0L) return(invisible())

  unreadable <- to_check[!vapply(columns[to_check], is_logmu_readable, logical(1L))]
  if (length(unreadable) == 0L) return(invisible())

  types <- vapply(unreadable, function(n) class(columns[[n]])[[1L]], character(1L))
  detail <- paste0("`", unreadable, "` (", types, ")", collapse = ", ")
  msg <- paste0(
    "logmu cannot read these column types: ", detail,
    ". Readable types are ", paste(logmu_readable_types, collapse = ", "), "."
  )

  switch(on_unreadable,
    error  = stop(msg, call. = FALSE),
    warn   = warning(msg, call. = FALSE),
    ignore = invisible()
  )
}
