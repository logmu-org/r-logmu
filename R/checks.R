# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

is_pure_logical <- function(x) is.logical(x) && is.null(attributes(x))
is_pure_numeric <- function(x) is.numeric(x) && is.null(attributes(x))
is_pure_double <- function(x) is.double(x) && is.null(attributes(x))
is_pure_text <- function(x) is.character(x) && is.null(attributes(x))
is_valid_datey <- function(x) datey::is_datey(x) && !anyNA(x)
is_single_valid_datey <- function(x) datey::is_datey(x) && length(x) == 1 && !is.na(x)
is_single_valid_durationy <- function(x) datey::is_durationy(x) && length(x) == 1 && !is.na(x)
is_single_pure_finite_numeric <- function(x) is_pure_numeric(x) && length(x) == 1 && is.finite(x)
is_single_valid_text <- function(x) is_pure_text(x) && length(x) == 1 && !is.na(x)
is_single_valid_logical <- function(x) is_pure_logical(x) && length(x) == 1 && !is.na(x)
is_matrix_of_double <- function(x) {
  if (!is.double(x) || length(x) == 0) return(FALSE)
  dims <- dim(x)
  if(is.null(dims) || length(dims) != 2) return(FALSE)
  if(dims[1] * dims[2] > 2147483647L) return(FALSE)
  TRUE
}
is_finite_matrix_of_double <- function(x) is_matrix_of_double(x) && cpp_all_finite(x)

get_single_valid_datey <- function(x) {
  if (is_single_valid_datey(x)) return(x)
  if (is_single_pure_finite_numeric(x)) {
    x2 <- datey::datey(x)
    if (is_single_valid_datey(x2)) return(x2)
  }
  stop("`", deparse(substitute(x)), "` must be a valid `datey` or year.", call. = FALSE)
}
get_single_valid_durationy <- function(x) {
  if (is_single_valid_durationy(x)) return(x)
  if (is_single_pure_finite_numeric(x)) {
    x2 <- datey::durationy(x)
    if (is_single_valid_durationy(x2)) return(x2)
  }
  stop("`", deparse(substitute(x)), "` must be a valid `durationy` or duration in years.", call. = FALSE)
}
get_birth_datey_from_.i <- function(.i) {
  .b <- .i[["birth"]]
  if (is.null(.b) || !is_single_valid_datey(.b)) {
      stop("`", deparse(substitute(.i)), "` must be a list with a valid `datey` named `birth`.", call. = FALSE)
  }
  .b
}


ensure_is_pure_numeric <- function(x) {
  if (!is_pure_numeric(x)) {
    arg_name <- deparse(substitute(x))
    stop("`", arg_name, "` must be a pure numeric.", call. = FALSE)
  }
}
ensure_is_pure_double <- function(x) {
  if (!is_pure_double(x)) {
    arg_name <- deparse(substitute(x))
    stop("`", arg_name, "` must be a pure double.", call. = FALSE)
  }
}
ensure_is_valid_datey <- function(x) {
  if (!is_valid_datey(x)) {
    arg_name <- deparse(substitute(x))
    stop("`", arg_name, "` must be valid `datey`(s).", call. = FALSE)
  }
}
ensure_is_single_valid_datey <- function(x) {
  if (!is_single_valid_datey(x)) {
    arg_name <- deparse(substitute(x))
    stop("`", arg_name, "` must be a (single) valid `datey`.", call. = FALSE)
  }
}
ensure_is_single_valid_durationy <- function(x) {
  if (!is_single_valid_durationy(x)) {
    arg_name <- deparse(substitute(x))
    stop("`", arg_name, "` must be a (single) valid `durationy`.", call. = FALSE)
  }
}
ensure_is_single_valid_text <- function(x) {
  if (!is_single_valid_text(x)) {
    arg_name <- deparse(substitute(x))
    stop("`", arg_name, "` must be a (single) text string.", call. = FALSE)
  }
}
ensure_is_valid_name <- function(x) {
  # TODO: Check not blank, no trailing, leading or consecutive spaces.
  if (!is_single_valid_text(x)) {
    arg_name <- deparse(substitute(x))
    stop("`", arg_name, "` must be a valid name, i.e. not blank and no trailing, leading or consecutive spaces.", call. = FALSE)
  }
}
ensure_is_finite_matrix_of_double <- function(x) {
  if (!is_finite_matrix_of_double(x)) {
    arg_name <- deparse(substitute(x))
    stop("`", arg_name, "` must be a matrix of `double`.", call. = FALSE)
  }
}
