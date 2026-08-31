# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

settings_class <- "logmu_settings"

# A quarter of a year: short enough that a mortality table's annual cells are
# sampled sensibly, long enough that a few years of exposure need only a dozen
# or so slots. Named once so `settings()` and `aev()` cannot drift apart.
default_time_scale <- 0.25

#' Settings for an analysis
#'
#' @description
#' A `settings` object carries the choices an analysis needs that are not part
#' of the question being asked: the overdispersion assumed, and the time scale
#' the integration uses.
#'
#' Build one once and pass it wherever it is wanted. Analytic functions take it
#' as a `settings` argument, and it is found among the arguments by its class,
#' so it needs no name and no position.
#'
#' @section Overdispersion:
#' Overdispersion, written \eqn{\Omega}, is defined by
#' \deqn{\mathrm{Var}(\mathrm{A}w - \mathrm{E}w) = \Omega\,\mathbb{E}\,\mathrm{E}w^2}
#' so it is the factor by which the variance of experience exceeds what
#' independent deaths under a deterministic mortality would give. It is
#' **required**, and has no default anywhere in **logmu**. Failing to allow for
#' it does not make results neutral -- it understates uncertainty by
#' \eqn{\sqrt\Omega} and selects overfitted models, so there is no safe value to
#' assume on a user's behalf.
#'
#' Values between 2 and 3 are usual for pensions longevity work, with higher
#' values making model selection more resistant to overfitting.
#'
#' @section Time scale:
#' The time scale is the width of one numerical integration interval. It may be
#' given as a `durationy` or as a number of years, and must be one of
#' 1, 1/4, 1/12 or 1/60 of a year.
#'
#' Those four are the intervals that, together with their halves, are a whole
#' number of clicks, so every sample point lands exactly on the click grid. They
#' also nest: refining from 1/4 to 1/12 to 1/60 keeps every sample already taken
#' and adds more between them.
#'
#' Smaller is more accurate and costs proportionally more. The default of a
#' quarter year is short enough to sample an annual mortality table sensibly.
#'
#' @param overdispersion A single positive number. Required.
#' @param time_scale The integration interval, as a `durationy` or a number of
#'   years. One of 1, 1/4, 1/12 or 1/60.
#' @param x An object.
#' @param ... Ignored.
#' @returns
#' `settings()` returns a `logmu_settings` object.
#'
#' `is_settings()` returns a scalar `logical`.
#' @examples
#' settings(overdispersion = 2)
#'
#' settings(overdispersion = 2.5, time_scale = 1 / 12)
#'
#' # A durationy says the same thing.
#' settings(overdispersion = 2.5, time_scale = datey::durationy(1 / 12))
#'
#' # Overdispersion is required.
#' try(settings(time_scale = 1))
#'
#' # And the time scale must be one of the four.
#' try(settings(overdispersion = 2, time_scale = 0.5))
#' @name settings
NULL

# The permitted scales come from the ENGINE, not from a list written out here.
# A second spelling in R is somewhere for the two to drift, and the symptom
# would be a scale that one accepts and the other refuses.
permitted_time_scale_clicks <- function() cpp_veil_time_scales()

# The permitted scales named the way a user wrote them, for an error message.
# Derived from the clicks so the message cannot describe a set that is no longer
# the set being enforced.
permitted_time_scales_text <- function() {
  clicks <- permitted_time_scale_clicks()
  per_year <- unclass(datey::durationy(1)) %/% clicks
  paste(ifelse(per_year == 1L, "1", paste0("1/", per_year)), collapse = ", ")
}

# A time scale in clicks, whatever form it was written in.
#
# THE CHECK IS EXACT INTEGER COMPARISON and needs no tolerance, because a
# `durationy` is stored as whole clicks: `durationy(1/12)` lands on exactly
# 44530 despite 1/12 being unrepresentable as a double. Routing a number through
# `durationy()` therefore forgives about half a click -- half a minute -- and
# nothing more.
time_scale_clicks <- function(time_scale) {
  if (!is_single_valid_durationy(time_scale) &&
      !is_single_pure_finite_numeric(time_scale)) {
    stop("`time_scale` must be a single `durationy` or a number of years.",
         call. = FALSE)
  }

  clicks <- unclass(get_single_valid_durationy(time_scale))

  if (!clicks %in% permitted_time_scale_clicks()) {
    stop("`time_scale` must be one of ", permitted_time_scales_text(),
         " of a year.", call. = FALSE)
  }

  clicks
}

# Overdispersion is REQUIRED and has no default -- see the note in the
# documentation above. `missing()` rather than a NULL default, so that the error
# a user meets names the argument rather than complaining about a NULL later on.
check_overdispersion <- function(overdispersion) {
  if (!is_single_pure_finite_numeric(overdispersion) || overdispersion <= 0) {
    stop("`overdispersion` must be a single positive number.", call. = FALSE)
  }
  as.double(overdispersion)
}

#' @rdname settings
#' @export
settings <- function(overdispersion, time_scale = default_time_scale) {
  if (missing(overdispersion)) {
    stop("`overdispersion` is required and has no default. ",
         "Values between 2 and 3 are usual for pensions longevity work.",
         call. = FALSE)
  }

  structure(
    list(
      overdispersion = check_overdispersion(overdispersion),
      time_scale_clicks = time_scale_clicks(time_scale)
    ),
    class = settings_class
  )
}

#' @rdname settings
#' @export
is_settings <- function(x) inherits(x, settings_class)

ensure_is_settings <- function(x) {
  if (!is_settings(x)) {
    arg_name <- deparse(substitute(x))
    stop("The argument `", arg_name, "` must be a `settings` object.", call. = FALSE)
  }
}

# The time scale as a durationy, for printing and for anyone who wants it back
# in the units it was given in.
settings_time_scale <- function(x) datey::durationy(x$time_scale_clicks / unclass(datey::durationy(1)))

#' @rdname settings
#' @usage NULL
#' @export
`$<-.logmu_settings` <- function(x, name, value) {
  stop("Modification of `settings` is invalid.", call. = FALSE)
}

#' @rdname settings
#' @usage NULL
#' @export
`[[<-.logmu_settings` <- function(x, i, value) {
  stop("Modification of `settings` is invalid.", call. = FALSE)
}

#' @rdname settings
#' @returns `x`, invisibly.
#' @export
print.logmu_settings <- function(x, ...) {
  u <- unclass(x)
  per_year <- unclass(datey::durationy(1)) %/% u$time_scale_clicks
  cat("<settings>\n")
  cat("  overdispersion: ", format(u$overdispersion), "\n", sep = "")
  cat("  time scale:     ",
      if (per_year == 1L) "1 year" else paste0("1/", per_year, " year"),
      "\n", sep = "")
  invisible(x)
}
