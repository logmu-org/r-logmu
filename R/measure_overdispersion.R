# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# measure_overdispersion -- a check on the overdispersion assumed a priori
# ============================================================================
#
# Overdispersion is the user's to choose, before any fitting, and it must be the
# same for every candidate being compared. This measures it from the data as a
# DIAGNOSTIC of that choice. Nothing here feeds anything else: the number it
# returns is for reading, and `overdispersion` stays a required argument
# everywhere it is used.

overdispersion_measure_class <- "overdispersion_measure"

# BELOW THIS MANY EFFECTIVE EXPECTED DEATHS a cell's Pearson term is far from
# normal, so the estimate is noisy and the chi-squared interval unreliable. Five
# is the usual rule of thumb for a Pearson statistic. A fixed constant rather
# than an argument, because it only decides whether to warn.
overdispersion_minimum_deaths <- 5

#' Measure overdispersion in experience data
#'
#' @description
#' Measures the overdispersion \eqn{\Omega} in experience data against a given
#' mortality, as a check on the value you have chosen to analyse and fit with.
#'
#' In **logmu** overdispersion is set *a priori*, and it must be the same for
#' every candidate model being compared. So this function is a diagnostic. It
#' tells you whether the data is broadly consistent with the \eqn{\Omega} you
#' chose, and so whether gaps in penalised log-likelihood between models are as
#' material as they look. Its result never feeds a fit.
#'
#' The arguments are those of [aev()], less `overdispersion` and plus
#' `parameter_count`.
#'
#' @section What it measures:
#' The data is divided into cells by `breakdown`, and A, E and V are worked out
#' in each cell as [aev()] would at an overdispersion of 1. The measure is the
#' Pearson statistic divided by its degrees of freedom:
#'
#' \deqn{\hat\Omega = \frac{1}{C - k}\sum_c \frac{(A_c - E_c)^2}{V_c}}
#'
#' where \eqn{C} is the number of cells with any exposure and \eqn{k} is
#' `parameter_count`. If the mortality is right and the overdispersion is
#' \eqn{\Omega}, each term has expectation about \eqn{\Omega}, whatever the
#' weight. A cell with no exposure, or a weight of zero throughout, has nothing
#' in it and is not counted.
#'
#' Because \eqn{\hat\Omega} rests on the residuals, it includes any lack of fit.
#' A mortality that misses a pattern the cells can see will show a larger
#' \eqn{\hat\Omega} than one that captures it, so the meaningful reading comes
#' from the richest model that fits adequately.
#'
#' @section Choosing the breakdown:
#' The cells are the elements of `breakdown`, and the choice matters:
#'
#' * \eqn{\hat\Omega} sees variation only along the dimensions the cells
#'   separate. Age cells cannot see a bad year, and calendar-year cells cannot
#'   see a poorly shaped age curve. Two breakdowns of the same data will give
#'   different measures.
#' * Cells that match a fitted model's own covariates carry no information. A
#'   model with a level for each sex makes A equal E in each sex, so cells by
#'   sex alone would have residuals of exactly zero. Useful cells cut across
#'   the covariates or are finer than them.
#' * The cells should not overlap. This is not checked.
#' * Each cell needs enough deaths. A cell with fewer than five *effective
#'   expected deaths*, \eqn{E_c^2 / V_c}, draws a warning. For lives this is
#'   simply \eqn{E_c}; with a weight it is the Kish effective count, the number
#'   of equally weighted expected deaths that would give the same precision.
#'   Small cells do not bias the measure, but they make it noisy and its
#'   interval unreliable, and a coarser breakdown is the remedy.
#'
#' With no `breakdown` there is a single cell, which leaves at most one degree
#' of freedom.
#'
#' @section The parameter count:
#' Fitting a parameter forces the residuals to satisfy one more equation, and
#' so removes about one cell's worth of free variation. `parameter_count` is
#' the number of parameters fitted *to this data* to produce `mortality`:
#'
#' * 0 for a standard table or any other mortality with nothing fitted;
#' * `length(coef(f))` for `f$mortality` from [fit()];
#' * in a layered fit, the parameters of the reference mortality as well, where
#'   it too was fitted to this data.
#'
#' The divisor \eqn{C - k} is exact for a count of lives and for any constant
#' weight. With a weight that varies it is an approximation, and so is the
#' interval.
#'
#' @section Relevance:
#' `val_relevance` and `val_distance` are used exactly as in [aev()]. A
#' relevance \eqn{r} multiplies the log-likelihood, which is to say that each
#' record is taken to have overdispersion \eqn{\Omega / r}, so V is already the
#' right yardstick. If the less relevant records in fact vary no more than the
#' rest, the measure will come out lower than it would without relevance. That
#' is the data saying the discount does not show in it.
#'
#' @section The interval:
#' If the terms were normal, the Pearson statistic would be \eqn{\Omega} times
#' a chi-squared variable on \eqn{C - k} degrees of freedom. The standard error
#' and the 95% interval for \eqn{\Omega} both rest on that, so they are a guide
#' rather than exact. Whether the value you chose lies inside the interval is
#' the question the interval answers.
#'
#' @param exp_data The experience data.
#' @param mortality A pronoun expression for \eqn{\log\mu}, or a `mortality`,
#'   as in [aev()]. Ordinarily the fitted mortality of the model you are
#'   checking.
#' @param include An `include` naming the population. `NULL` means everybody.
#' @param weight A pronoun expression for the weight \eqn{w}, or a `variable`.
#'   `NULL` means a weight of 1, i.e. a count of lives. The measure is the
#'   overdispersion for this weight.
#' @param val_relevance,val_distance Two spellings of the relevance \eqn{r},
#'   as in [aev()]. See *Relevance*.
#' @param breakdown An `include` or [includes()] dividing the population into
#'   cells. See *Choosing the breakdown*.
#' @param parameter_count The number of parameters fitted to this data to
#'   produce `mortality`. Required. See *The parameter count*.
#' @param settings A [settings()] object supplying `time_scale`. Any
#'   `overdispersion` it carries is not used.
#' @param time_scale Given directly, this overrides `settings`.
#' @param threads Worker threads to use. `0` asks for as many as the machine
#'   reports. Cannot change any answer.
#' @param x An `overdispersion_measure`.
#' @param ... Ignored.
#' @returns
#' `measure_overdispersion()` returns an `overdispersion_measure`, a list with
#' fields `overdispersion` (the measured \eqn{\hat\Omega}), `std_error` (one
#' standard deviation), `lower_95pc` and `upper_95pc` (the 95% interval),
#' `cells` (\eqn{C}), `parameter_count` (\eqn{k}) and `degrees_of_freedom`
#' (\eqn{C - k}). It warns if any cell has fewer than five effective expected
#' deaths, and stops if \eqn{C - k} is not positive.
#'
#' `print()` returns `x` invisibly.
#' @examples
#' set.seed(1)
#' lives <- data.frame(
#'   birth = datey::datey(2015 - sample(60:85, 10000, replace = TRUE))
#' )
#' truth <- gompertz_mortality(log_mu_0 = -3.6, slope_t = 0)
#'
#' data <- simulate(
#'   truth,
#'   lives = lives,
#'   exp_start = datey::datey(2015),
#'   exp_end = datey::datey(2020))
#'
#' # Simulated deaths are independent, so the true overdispersion is 1 and the
#' # interval should usually contain it.
#' measure_overdispersion(
#'   data,
#'   truth,
#'   breakdown = ages(60, 90, by = 2),
#'   parameter_count = 0)
#' @export
measure_overdispersion <- function(
    exp_data,
    mortality,
    include = NULL,
    weight = NULL,
    val_relevance = NULL,
    val_distance = NULL,
    breakdown = NULL,
    parameter_count,
    settings = NULL,
    time_scale = NULL,
    threads = cpp_veil_default_threads()) {

  ensure_is_exp_data(exp_data)

  if (missing(mortality)) {
    stop("`mortality` is required.", call. = FALSE)
  }
  if (missing(parameter_count)) {
    stop("`parameter_count` is required: the number of parameters fitted to this data ",
         "to produce `mortality`, or 0 if none were.", call. = FALSE)
  }
  parameter_count <- check_parameter_count(parameter_count)

  # AT AN OVERDISPERSION OF ONE, given in the call so that it wins over any a
  # `settings` carries. V is then E w^2 itself, and the Pearson sum over the
  # degrees of freedom is the measure, with no assumed value to divide out.
  resolved <- resolve_settings(settings, 1, time_scale)

  # The same capture and the same plan as `aev()`, so the cells here cannot
  # drift from the A/E a user would compute for themselves.
  caller <- parent.frame()
  plan <- aev_plan(
    it_capture(substitute(mortality), caller),
    aev_optional_ast(substitute(weight), missing(weight), caller),
    aev_optional_ast(substitute(val_relevance), missing(val_relevance), caller),
    aev_optional_ast(substitute(val_distance), missing(val_distance), caller),
    include,
    breakdown,
    resolved$overdispersion)

  run <- cpp_veil_run(
    plan$specs,
    exp_data_columns(exp_data),
    resolved$time_scale_clicks,
    FALSE,
    as.integer(threads))

  overdispersion_from_cells(plan$finalise(run$results), parameter_count)
}

# A single whole number, not negative, returned as an integer to match the `k`
# column of `compare_models()`.
check_parameter_count <- function(parameter_count) {
  if (!is.numeric(parameter_count) || length(parameter_count) != 1L ||
        !is.finite(parameter_count) || parameter_count < 0 ||
        parameter_count != round(parameter_count)) {
    stop("`parameter_count` must be a single whole number, 0 or more.", call. = FALSE)
  }
  as.integer(parameter_count)
}

# The measure from the cells' A, E and V, all at an overdispersion of one.
overdispersion_from_cells <- function(cells, parameter_count) {
  values <- unclass(cells)
  A <- values$A
  E <- values$E
  V <- values$V

  # A CELL WITH NOTHING IN IT says nothing. No exposure, or a weight of zero
  # throughout, gives A = E = V = 0 and a term of 0/0, so it is left out and not
  # counted. A cell with V of zero and deaths in it is KEPT: that is a mortality
  # that underflowed against deaths that happened, its term is infinite and so
  # is the measure, which is a true statement about an impossible model. NaN is
  # kept too, and flows through, as everywhere else.
  empty <- !is.na(A) & !is.na(V) & A == 0 & V == 0
  cell_count <- sum(!empty)

  degrees <- cell_count - parameter_count
  if (degrees <= 0L) {
    stop(sprintf(paste0(
      "There are no degrees of freedom left to measure overdispersion with: %d cell%s ",
      "with exposure, less %d fitted parameter%s. Give a `breakdown` with more cells ",
      "than parameters."),
      cell_count, if (cell_count == 1L) "" else "s",
      parameter_count, if (parameter_count == 1L) "" else "s"), call. = FALSE)
  }

  measure <- sum(((A - E)^2 / V)[!empty]) / degrees

  # EFFECTIVE EXPECTED DEATHS, E^2 / V: simply E for lives, and the Kish
  # effective count under a weight. NaN where both are zero, which the empty
  # mask has already set aside.
  effective <- E^2 / V
  small <- !empty & !is.na(effective) & effective < overdispersion_minimum_deaths
  if (any(small)) {
    labels <- overdispersion_cell_labels(cells)
    warning(sprintf(paste0(
      "%d cell%s %s fewer than %d effective expected deaths (%s), so the measure is ",
      "noisy and its interval unreliable. A coarser `breakdown` would help."),
      sum(small), if (sum(small) == 1L) "" else "s", if (sum(small) == 1L) "has" else "have",
      overdispersion_minimum_deaths, paste(labels[small], collapse = ", ")), call. = FALSE)
  }

  structure(
    list(
      overdispersion = measure,
      std_error = measure * sqrt(2 / degrees),
      lower_95pc = measure * degrees / stats::qchisq(0.975, degrees),
      upper_95pc = measure * degrees / stats::qchisq(0.025, degrees),
      cells = cell_count,
      parameter_count = parameter_count,
      degrees_of_freedom = degrees
    ),
    class = overdispersion_measure_class
  )
}

# The cells' names for a message, falling back to their positions where a
# hand-built `includes` carries none.
overdispersion_cell_labels <- function(cells) {
  count <- length(cells)
  given <- names(cells)
  if (is.null(given)) {
    given <- if (count == 1L) "all records" else rep(NA_character_, count)
  }
  ifelse(is.na(given) | !nzchar(given), paste("cell", seq_len(count)), given)
}

#' @rdname measure_overdispersion
#' @export
print.overdispersion_measure <- function(x, ...) {
  number <- function(value) formatC(value, format = "f", digits = 3L)
  cat("<overdispersion measure>\n")
  cat(sprintf("  overdispersion      %s\n", number(x$overdispersion)))
  cat(sprintf("  standard error      %s (1 sd)\n", number(x$std_error)))
  cat(sprintf("  95%% interval        %s to %s\n",
              number(x$lower_95pc), number(x$upper_95pc)))
  cat(sprintf("  degrees of freedom  %d (%d cell%s less %d parameter%s)\n",
              x$degrees_of_freedom,
              x$cells, if (x$cells == 1L) "" else "s",
              x$parameter_count, if (x$parameter_count == 1L) "" else "s"))
  invisible(x)
}
