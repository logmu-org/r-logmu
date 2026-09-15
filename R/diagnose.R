# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# diagnose -- the standard diagnostics of every fitted candidate
# ============================================================================
#
# A fit keeps no data, deliberately, so anything that needs the data again -- an
# A/E over a breakdown, above all -- has to be asked for with the data to hand.
# `summary()` gives what needs no data; this gives the rest, for a lone fit or
# for every candidate of a comparison at once.

diagnosis_class <- "diagnosis"

#' Standard diagnostics for a fit, or for every candidate of a comparison
#'
#' @description
#' Runs the A/E diagnostics of a fitted model over the experience data it was
#' fitted to, broken down however you ask, and returns them alongside the
#' model's coefficient table and correlation matrix. Given a
#' `model_comparison`, it does so for **every** converged candidate, not just
#' the winner, and the whole set crosses the data once.
#'
#' A `fit` keeps no data, which is why the data is passed again here rather
#' than remembered. What the fit *does* keep is what it was run with -- its
#' include, weight, relevance, overdispersion and time scale -- and those are
#' the defaults, so the diagnostics of a fit are by default over the same
#' population, weighted the same way, as the fit itself.
#'
#' @section Breakdowns:
#' `breakdown` may be a single `include` or [includes()], or a list of them for
#' several breakdowns side by side. Each gives one [aev()] per candidate, with
#' one record per element. A list's names label the breakdowns; unnamed ones are
#' labelled by what was written. `NULL` gives one ungrouped A/E per candidate.
#'
#' @section What it shows:
#' `print()` gives each candidate's coefficient table and correlation matrix,
#' then each of its A/E tables. `plot()` draws a grid of A/E charts, candidates
#' down the page and breakdowns across it, each drawn by [autoplot.aev()].
#'
#' @param exp_data The experience data, ordinarily the data the model was fitted
#'   to.
#' @param x A `fit` or a `model_comparison`.
#' @param breakdown An `include`, an [includes()], a list of either, or `NULL`.
#'   See *Breakdowns*.
#' @param include The population. Defaults to the one the fit was run over.
#' @param weight,val_relevance,val_distance Pronoun expressions, as in [aev()].
#'   Each defaults to what the fit was run with; give `NULL` to run without.
#' @param settings,overdispersion,time_scale As in [aev()]. Each defaults to the
#'   fit's own.
#' @param threads The number of worker threads.
#' @param ... For `plot()`, passed to [autoplot.aev()] for every chart, except
#'   `title` and `subtitle`, which name the candidate and the breakdown.
#' @returns
#' `diagnose()` returns a `diagnosis`, whose fields are `candidates` (one entry
#' per converged candidate, best first, each holding `summary` -- the
#' [summary()][fit] of its fit -- and `aev`, a list of one [aev()] per
#' breakdown), `breakdowns` (the breakdown labels), `labels` (the include and
#' weight the diagnostics were run with, as text), `overdispersion`, and
#' `pv_impact`, which is `NULL` until present-value impacts exist.
#'
#' `plot()` draws the chart grid and returns `x` invisibly.
#' @examples
#' data <- exp_data(
#'   list(
#'     birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948, 1952, 1938, 1960, 1943)),
#'     male      = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
#'     sex       = c("m", "f", "m", "f", "m", "m", "f", "m", "f"),
#'     E2R_start = datey::datey(rep(2015, 9)),
#'     E2R_end   = datey::datey(c(2020, 2018, 2020, 2019, 2020, 2020, 2020, 2017, 2020)),
#'     E2R_died  = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
#'   ),
#'   exp_start = datey::datey(2015),
#'   exp_end   = datey::datey(2020)
#' )
#'
#' reference <- gompertz_mortality(slope_t = 0)
#' comparison <- compare_models(
#'   data,
#'   list(level = model(reference, covariates(level = 1)),
#'        sex   = model(reference, covariates(male = .i$male, female = !.i$male))),
#'   overdispersion = 1
#' )
#'
#' d <- diagnose(data, comparison,
#'               breakdown = list(sex = categories(.i$sex, "m", "f"), age = ages(55, 85, by = 15)))
#' d
#' @name diagnose
NULL

#' @rdname diagnose
#' @export
diagnose <- function(exp_data,
                     x,
                     breakdown = NULL,
                     include,
                     weight,
                     val_relevance,
                     val_distance,
                     settings = NULL,
                     overdispersion = NULL,
                     time_scale = NULL,
                     threads = cpp_veil_default_threads()) {

  ensure_is_exp_data(exp_data)
  if (missing(x)) stop("`x` is required: a `fit` or a `model_comparison`.", call. = FALSE)
  fits <- diagnose_fits(x)
  defaults <- x$arguments

  # THE FIT'S OWN ARGUMENTS ARE THE DEFAULTS, and `missing()` is what tells "not
  # given, use the fit's" from "given as NULL, use none" -- a distinction
  # `aev()` does not need, because there nothing is remembered. Captured in the
  # caller's frame, as every pronoun argument is.
  caller <- parent.frame()
  given <- function(expr, absent, default) {
    if (absent) default else aev_optional_ast(expr, FALSE, caller)
  }
  weight_ast <- given(substitute(weight), missing(weight), defaults$weight)
  relevance_ast <- given(substitute(val_relevance), missing(val_relevance),
                          defaults$val_relevance)
  distance_ast <- given(substitute(val_distance), missing(val_distance), defaults$val_distance)

  include_absent <- missing(include)
  include_label <- if (include_absent) x$labels$include
                   else if (is.null(include)) "all records"
                   else it_short_deparse(substitute(include))
  if (include_absent) include <- defaults$include

  resolved <- diagnose_settings(settings, overdispersion, time_scale, x$overdispersion,
                                defaults$time_scale_clicks)

  breakdowns <- diagnose_breakdowns(breakdown, substitute(breakdown))

  # ONE PLAN PER CANDIDATE PER BREAKDOWN, AND ONE PASS FOR ALL OF THEM. This is
  # the seam `batch()` runs through: every plan's specifications are compiled and
  # scheduled together, the data is read once, and each plan reads back its own
  # slice of the answers. Ten candidates on three breakdowns is one crossing.
  plans <- lapply(fits, function(fit) {
    lapply(breakdowns, function(elements) {
      aev_plan(it_obj(fit$mortality), weight_ast, relevance_ast, distance_ast,
               include, elements, resolved$overdispersion)
    })
  })
  flat <- unlist(plans, recursive = FALSE)
  specs <- unlist(lapply(flat, function(plan) plan$specs), recursive = FALSE)

  run <- cpp_veil_run(specs, exp_data_columns(exp_data), resolved$time_scale_clicks,
                      FALSE, as.integer(threads))

  taken <- 0L
  results <- lapply(flat, function(plan) {
    wanted <- length(plan$specs)
    result <- plan$finalise(run$results[taken + seq_len(wanted)])
    taken <<- taken + wanted
    result
  })

  candidates <- lapply(seq_along(fits), function(at) {
    aevs <- results[(at - 1L) * length(breakdowns) + seq_along(breakdowns)]
    names(aevs) <- names(breakdowns)
    list(summary = summary(fits[[at]]), aev = aevs)
  })
  names(candidates) <- names(fits)

  structure(
    list(
      candidates = candidates,
      breakdowns = names(breakdowns),
      labels = list(
        include = include_label,
        weight = if (missing(weight)) x$labels$weight
                 else if (is.null(weight_ast)) "lives"
                 else it_deparse(weight_ast)
      ),
      overdispersion = resolved$overdispersion,

      # A SLOT FOR WHAT IS TO COME, so the shape of the object does not change
      # when present-value impacts land.
      pv_impact = NULL
    ),
    class = diagnosis_class
  )
}

# The converged fits to diagnose, named. A lone fit is a comparison of one.
diagnose_fits <- function(x) {
  if (is_fit(x)) return(list(fit = x))
  if (is_model_comparison(x)) {
    fits <- Filter(Negate(is.null), x$fits)
    if (length(fits) == 0L) {
      stop("No candidate of the comparison converged, so there is nothing to diagnose.",
           call. = FALSE)
    }
    return(fits)
  }
  stop("`x` must be a `fit` or a `model_comparison`.", call. = FALSE)
}

# `resolve_settings()` with the fit's own values as the last resort: the call
# wins, then a `settings`, then what the fit was run with. `overdispersion` is
# never absent here, because a fit cannot exist without one.
diagnose_settings <- function(settings, overdispersion, time_scale,
                              fit_overdispersion, fit_clicks) {
  if (!is.null(settings)) ensure_is_settings(settings)

  overdispersion <- overdispersion %||% settings$overdispersion %||% fit_overdispersion

  clicks <- if (!is.null(time_scale)) {
    time_scale_clicks(time_scale)
  } else if (!is.null(settings)) {
    settings$time_scale_clicks
  } else {
    fit_clicks
  }

  list(overdispersion = check_overdispersion(overdispersion), time_scale_clicks = clicks)
}

# The breakdowns as a named list, each element an `includes` or NULL.
#
# LABELLED BY WHAT WAS WRITTEN where no name was given, because the label heads
# a chart and `breakdown 2` heads nothing. A list's names win; failing those,
# the expression of each element; failing that, its position.
diagnose_breakdowns <- function(breakdown, expr) {
  if (is.null(breakdown)) return(list(all = NULL))

  if (is_include(breakdown) || is_includes(breakdown)) {
    result <- list(as_includes(breakdown))
    names(result) <- it_short_deparse(expr)
    return(result)
  }

  if (!is.list(breakdown)) {
    stop("`breakdown` must be an `include`, an `includes`, a list of them, or NULL.",
         call. = FALSE)
  }

  written <- if (is.call(expr) && identical(expr[[1L]], as.name("list")) &&
                 length(expr) == length(breakdown) + 1L) {
    vapply(as.list(expr)[-1L], it_short_deparse, character(1L))
  } else {
    paste("breakdown", seq_along(breakdown))
  }
  given <- names(breakdown) %||% rep("", length(breakdown))
  labels <- ifelse(nzchar(given), given, written)

  result <- lapply(breakdown, function(one) {
    if (is.null(one)) NULL else as_includes(one)
  })
  names(result) <- labels
  result
}

#' @rdname diagnose
#' @export
is_diagnosis <- function(x) inherits(x, diagnosis_class)

#' @rdname diagnose
#' @export
print.diagnosis <- function(x, ...) {
  cat(sprintf("<diagnosis: %d candidate%s, %d breakdown%s>\n",
              length(x$candidates), if (length(x$candidates) == 1L) "" else "s",
              length(x$breakdowns), if (length(x$breakdowns) == 1L) "" else "s"))
  cat(sprintf("  include         %s\n", x$labels$include))
  cat(sprintf("  weight          %s\n", x$labels$weight))
  cat(sprintf("  overdispersion  %s\n", format(x$overdispersion, digits = 5L)))

  for (label in names(x$candidates)) {
    candidate <- x$candidates[[label]]
    cat(sprintf("\n---- %s\n\n", label))
    print(candidate$summary$coefficients, digits = 4L)
    cat_correlation(candidate$summary$correlation)

    for (name in names(candidate$aev)) {
      cat(sprintf("\nA/E by %s\n", name))
      print(candidate$aev[[name]])
    }
  }
  invisible(x)
}

# The charts of a diagnosis, candidates down and breakdowns across, in row-major
# order. Each is an `autoplot.aev()` titled with the candidate and subtitled
# with the breakdown; everything else in `...` goes to every chart alike.
diagnosis_charts <- function(x, ...) {
  extra <- list(...)
  if (any(c("title", "subtitle") %in% names(extra))) {
    stop("`title` and `subtitle` are set by `plot()` on a diagnosis: the candidate ",
         "and the breakdown name each chart.", call. = FALSE)
  }
  unlist(lapply(names(x$candidates), function(label) {
    lapply(x$breakdowns, function(name) {
      do.call(ggplot2::autoplot,
              c(list(x$candidates[[label]]$aev[[name]], title = label, subtitle = name), extra))
    })
  }), recursive = FALSE)
}

#' @rdname diagnose
#' @export
plot.diagnosis <- function(x, ...) {
  charts <- diagnosis_charts(x, ...)
  rows <- length(x$candidates)
  columns <- length(x$breakdowns)

  # A GRID OF VIEWPORTS rather than a facet, because each chart is a complete
  # `ggplot` with its own two panels and its own type-driven sizing, and those
  # do not survive being faceted. Each is printed into its own cell.
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(rows, columns)))
  for (row in seq_len(rows)) {
    for (column in seq_len(columns)) {
      print(charts[[(row - 1L) * columns + column]],
            vp = grid::viewport(layout.pos.row = row, layout.pos.col = column))
    }
  }
  grid::upViewport()
  invisible(x)
}
