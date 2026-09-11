# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# compare_models -- ranking candidate models on penalised log-likelihood
# ============================================================================
#
# THERE IS NO NEW STATISTIC HERE. `fit()` already returns
# `penalised_log_likelihood`, which is `L(beta_hat) - p` with
# `p = Z^-1 tr(J I^-1)`, and that IS the selection criterion -- the AIC without
# the factor of -2. What this function adds is the machinery that makes several
# of those numbers comparable:
#
#   * ONE Z, computed once and handed to every candidate. Given as a test
#     mortality it would otherwise be re-measured per candidate over the same
#     data, which is wasted work and, once a candidate changes the include or
#     weight, a different scale for each -- so `resolve_Z()` runs here,
#     before the loop, and what reaches the engine is a number.
#   * ONE dataset, include, weight and overdispersion, because they are
#     arguments of this call rather than of a `model`. A `model` deliberately
#     holds nothing but `ref_mortality` and `covariates`, which is what makes
#     candidates comparable by construction rather than by a rule in the
#     documentation.
#
# A FAILING CANDIDATE LANDS ON THE CANDIDATE, NOT THE COMPARISON. `fit()` raises,
# because a fit returns a fit or it fails; here one non-converging candidate out
# of fifty must not destroy the other forty-nine. So the row carries the status
# and `NA`, the diagnosis is kept in `failures`, and the rest still rank.
#
# ONE ENGINE RUN PER CANDIDATE, deliberately for now. Fitting the candidates
# together in a single crossing needs a multi-block fit loop in the engine, with
# per-block convergence; the threading pool already parallelises within each fit,
# and sharing Z and the argument resolution is what this slice is for.

comparison_class <- "model_comparison"

#' Compare mortality models by penalised log-likelihood
#'
#' @description
#' Fits several candidate [model()]s against one set of experience data and
#' ranks them on penalised log-likelihood,
#' \deqn{L_{\mathrm{P}} = L(\hat\beta) - p, \qquad p = Z^{-1}\mathrm{tr}(\mathbf{J}\mathbf{I}^{-1}).}
#' Larger is better. This is Akaike's information criterion without the factor
#' of \eqn{-2}, generalised to an arbitrary log-likelihood weight.
#'
#' @details
#' # What makes the candidates comparable
#'
#' Every candidate is fitted over the same data, `include`, weight,
#' overdispersion and \eqn{Z}, because all of those are arguments of this call
#' and none of them is part of a `model`. \eqn{Z} in particular is resolved once
#' and shared: it is measured on `test_mortality` over this data one time, not
#' once per candidate.
#'
#' Candidates may differ in their reference mortality. The constant
#' \eqn{\mathrm{A}w\log\mu^{\mathrm{ref}}} is retained in \eqn{L}, so a
#' comparison across different base tables is meaningful.
#'
#' # Reading the shortfall
#'
#' Where the weight is absent or is an indicator, \eqn{w^2 = w} makes \eqn{p}
#' exactly `k`, \eqn{L_{\mathrm{P}}} exactly \eqn{-}AIC\eqn{/2}, and a shortfall
#' of 1 one parameter's worth -- which is the usual threshold for a difference
#' that matters. For an *ad hoc* weight the penalty is no longer a count and the
#' shortfall carries no such calibration; it still ranks the candidates, but how
#' large a shortfall is worth acting on is a question this function does not
#' answer.
#'
#' Overdispersion does real work in the ranking. The penalty is free of
#' \eqn{\Omega} and the fit term is not, so a larger \eqn{\Omega} favours
#' simpler models, which is the quasi-likelihood behaviour and is why
#' `overdispersion` is yours to supply.
#'
#' @param exp_data The experience data.
#' @param models A list of [model()]s, ideally named. Names appear in the
#'   ranking and index `fits`; unnamed candidates are labelled by position.
#' @param include An `include` restricting the population, or `NULL` for all of
#'   it.
#' @param weight A pronoun expression for the weight \eqn{w}, or `NULL` for
#'   lives. Weights must not be negative.
#' @param val_similarity,val_distance A pronoun expression for the second
#'   weighting factor, in either spelling. At most one may be given.
#' @param test_mortality The mortality on which the \eqn{Z} scale is measured,
#'   once, for every candidate. See [fit()].
#' @param Z The \eqn{Z} scale as a single positive number, given instead
#'   of a `test_mortality` rather than as well as one.
#' @param settings A [settings()] object carrying `overdispersion` and
#'   `time_scale`.
#' @param overdispersion,time_scale Given directly, these override the
#'   `settings`.
#' @param L_tolerance The convergence tolerance, as a gain in \eqn{L}.
#' @param armijo The fraction of the linearly predicted gain a damped step must
#'   actually deliver to be accepted.
#' @param max_iterations The most Newton steps to take.
#' @param max_halvings The most times one step may be halved before a candidate
#'   is abandoned.
#' @param threads The number of worker threads.
#' @param x A `model_comparison`.
#' @param ... Ignored.
#' @returns
#' `compare_models()` returns a `model_comparison`, whose fields are `table` (the
#' ranking, best first, with columns `k`, `log_likelihood`, `penalty`,
#' `penalised_log_likelihood` and `status`), `fits` (the `fit` of each candidate,
#' in the same order, `NULL` where one failed), `failures` (the diagnosis of each
#' failure, named, empty where there were none), `Z`, `overdispersion`,
#' `diagnostics` (the population and test-mortality numbers the header reports)
#' and `labels` (what was written for the include, the weight and the test
#' mortality).
#'
#' `table` holds the values as computed. `print()` takes the best candidate's
#' \eqn{L_{\mathrm{P}}} off both likelihood columns, because differences are what
#' matter and the absolute numbers are large and uninformative -- so the winner
#' prints as 0. Both columns shift by the same constant, which keeps
#' \eqn{L - L_{\mathrm{P}} = p} true row by row.
#'
#' `is_model_comparison()` returns a scalar `logical`.
#'
#' `best_fit()` returns the `fit` of the highest-ranked candidate.
#' @examples
#' data <- exp_data(
#'   list(
#'     birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
#'     male      = c(TRUE, FALSE, TRUE, FALSE, TRUE),
#'     E2R_start = datey::datey(rep(2015, 5)),
#'     E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
#'     E2R_died  = c(FALSE, TRUE, TRUE, FALSE, FALSE)
#'   ),
#'   exp_start = datey::datey(2015),
#'   exp_end   = datey::datey(2020)
#' )
#'
#' # One death of each sex, so the split is identifiable.
#' reference <- gompertz_mortality()
#' comparison <- compare_models(
#'   data,
#'   models = list(
#'     level = model(reference, covariates(level = 1)),
#'     sex   = model(reference, covariates(male = .i$male, female = !.i$male))
#'   ),
#'   settings = settings(overdispersion = 1)
#' )
#' comparison
#' coef(best_fit(comparison))
#' @name compare_models
NULL

#' @rdname compare_models
#' @export
compare_models <- function(exp_data,
                           models,
                           include = NULL,
                           weight = NULL,
                           val_similarity = NULL,
                           val_distance = NULL,
                           test_mortality = default_mortality(),
                           Z = NULL,
                           settings = NULL,
                           overdispersion = NULL,
                           time_scale = NULL,
                           L_tolerance = default_L_tolerance,
                           armijo = default_armijo,
                           max_iterations = default_max_iterations,
                           max_halvings = default_max_halvings,
                           threads = cpp_veil_default_threads()) {

  ensure_is_exp_data(exp_data)
  if (missing(models)) stop("`models` is required.", call. = FALSE)
  labels <- comparison_labels(models)

  resolved <- resolve_settings(settings, overdispersion, time_scale)

  # CAPTURED FROM THE CALLER'S FRAME, exactly as `fit()` does it and for the same
  # reason: `substitute()` on an unsupplied argument yields its default, so a
  # value test cannot tell `weight = NULL` from no weight.
  caller <- parent.frame()
  weight_ast <- aev_optional_ast(substitute(weight), missing(weight), caller)
  similarity_ast <- aev_optional_ast(substitute(val_similarity), missing(val_similarity), caller)
  distance_ast <- aev_optional_ast(substitute(val_distance), missing(val_distance), caller)
  test_ast <- fit_test_mortality_ast(substitute(test_mortality), missing(test_mortality),
                                     Z, caller)

  if (!is.null(include)) {
    if (is_includes(include)) {
      stop("`include` takes a single `include`. Every candidate is fitted over one population.",
           call. = FALSE)
    }
    ensure_is_include(include)
  }

  columns <- exp_data_columns(exp_data)

  # ONCE, FOR EVERY CANDIDATE. This is the whole mechanical reason this function
  # exists rather than a loop over `fit()` in user code, and it is also where the
  # header's population numbers come from.
  diagnostics <- fit_diagnostics(test_ast, weight_ast, include, columns,
                                 resolved$time_scale_clicks, threads)
  scale <- resolve_Z(Z, weight_ast, diagnostics)
  call_labels <- fit_call_labels(substitute(test_mortality), missing(test_mortality),
                                 substitute(include), include, weight_ast)

  # EVERY CANDIDATE STARTS AT `beta = 0`, AND THERE IS NO `start` ARGUMENT HERE.
  # `fit()` has one; a comparison deliberately does not. Newton-Raphson with
  # relaxation converges fast enough from the reference mortality that the
  # starting point does not matter in practice, and a per-candidate `start` can
  # only be a list, because candidates differ in how many terms they have --
  # which then wants naming, and a named list read positionally would produce a
  # slower fit rather than an error. If a candidate ever does need help, the
  # user's move is to translate `log mu^ref` by a weighted sum of the covariates
  # so that `beta = 0` is where the data already is.
  # `unname()` so the runs carry no names into the table's columns; the labels are
  # the row names and are applied once, in `new_comparison()`.
  runs <- lapply(unname(models), function(model) {
    fit_engine_run(model, include, weight_ast, similarity_ast, distance_ast,
                   columns, resolved, scale, NULL,
                   L_tolerance, armijo, max_iterations, max_halvings, threads)
  })

  new_comparison(runs, models, labels, scale, resolved$overdispersion, L_tolerance,
                 diagnostics, call_labels)
}

# ---- the candidates --------------------------------------------------------

# One label per candidate, and they must be UNIQUE, because they index `fits`.
# A candidate with no name is labelled by position rather than refused: a
# comparison is often built up programmatically, and `lapply()` over a set of
# covariate lists produces a list with no names at all.
comparison_labels <- function(models) {
  # A `model` IS A LIST, so this test has to come first or a single model handed
  # in here is silently taken apart and its `ref_mortality` reported as candidate
  # one.
  if (is_model(models)) {
    stop("`models` must be a list of `model` objects; one `model` was given. ",
         "Use `fit()` for a single model, or `list()` to compare it with others.",
         call. = FALSE)
  }
  if (!is.list(models)) {
    stop("`models` must be a list of `model` objects.", call. = FALSE)
  }
  if (length(models) < 1L) {
    stop("`models` is empty; there is nothing to compare.", call. = FALSE)
  }

  for (at in seq_along(models)) {
    if (!is_model(models[[at]])) {
      stop(sprintf("`models[[%d]]` is not a `model`, as built by `model()`.", at),
           call. = FALSE)
    }
  }

  labels <- names(models)
  if (is.null(labels)) labels <- rep("", length(models))
  anonymous <- is.na(labels) | !nzchar(labels)
  labels[anonymous] <- sprintf("model_%d", seq_along(models))[anonymous]

  duplicated_labels <- unique(labels[duplicated(labels)])
  if (length(duplicated_labels) > 0L) {
    stop("Every candidate needs its own name; ",
         paste(sprintf("`%s`", duplicated_labels), collapse = ", "),
         " is used more than once.", call. = FALSE)
  }
  labels
}

# ---- the result ------------------------------------------------------------

new_comparison <- function(runs, models, labels, Z, overdispersion, L_tolerance,
                           diagnostics, call_labels) {

  converged <- vapply(runs, function(run) identical(run$status, "converged"),
                      logical(1L))

  fits <- vector("list", length(runs))
  fits[converged] <- lapply(which(converged), function(at) {
    new_fit(runs[[at]], models[[at]], Z, overdispersion, diagnostics, call_labels)
  })
  names(fits) <- labels

  failures <- vapply(which(!converged), function(at) {
    fit_failure_message(runs[[at]], models[[at]], L_tolerance)
  }, character(1L))
  names(failures) <- labels[!converged]

  # WARNED ABOUT ONCE, NAMING THE CANDIDATES. The comparison is still returned,
  # but a candidate silently absent from the ranking is how somebody comes to
  # believe they compared fifty models when they compared forty-nine. The
  # diagnosis itself stays in `failures` rather than going into the warning,
  # which would be several paragraphs long at fifty candidates.
  if (length(failures) > 0L) {
    warning(sprintf("%d of %d candidate%s could not be fitted: %s. See `$failures` for why.",
                    length(failures), length(runs),
                    if (length(runs) == 1L) "" else "s",
                    paste(sprintf("`%s` (%s)", names(failures),
                                  vapply(runs[!converged], function(run) run$status,
                                         character(1L))),
                          collapse = ", ")),
            call. = FALSE)
  }

  # NA WHERE A CANDIDATE FAILED, AND THE FIELDS ARE NOT READ. They are all
  # present -- the engine fills them whatever the status -- but they hold
  # whatever the loop had reached when it gave up. A run that stops at
  # `not_identifiable` returns a log-likelihood at a point that was never
  # accepted and a penalty of NaN, so putting either in the table would be
  # reporting a number nobody should compare.
  from_run <- function(field) {
    vapply(seq_along(runs), function(at) {
      if (converged[[at]]) as.double(runs[[at]][[field]]) else NA_real_
    }, double(1L))
  }

  log_likelihood <- from_run("log_likelihood")
  penalty <- from_run("penalty")
  penalised <- log_likelihood - penalty

  # RAW, UNSHIFTED VALUES. `print()` deducts the maximum `L_P` from the two
  # likelihood columns, because differences are the whole point and the absolute
  # numbers are large and uninformative -- but the table keeps what was actually
  # computed, so nobody needs an argument to get at it.
  table <- data.frame(
    k = unname(vapply(models, function(model) length(model$covariates), integer(1L))),
    log_likelihood = log_likelihood,
    penalty = penalty,
    penalised_log_likelihood = penalised,
    status = vapply(runs, function(run) run$status, character(1L)),
    row.names = labels,
    stringsAsFactors = FALSE
  )

  # BEST FIRST, and a failure last whatever its status: `order()` puts `NA` at the
  # end by default, which is exactly the wanted behaviour and the reason the sort
  # is keyed on a likelihood rather than on `status`. `decreasing` because larger
  # `L_P` is better.
  ranking <- order(penalised, decreasing = TRUE, na.last = TRUE)

  structure(
    list(table = table[ranking, , drop = FALSE],
         fits = fits[ranking],
         failures = failures,
         Z = Z,
         overdispersion = overdispersion,
         diagnostics = diagnostics,
         labels = call_labels),
    class = comparison_class
  )
}

# The best `L_P`, or `NA` when nothing converged -- in which case there is no
# maximum to deduct and the shown columns are `NA` throughout, which is the
# honest answer rather than an arbitrary winner reading zero.
comparison_best_L_P <- function(x) {
  penalised <- x$table$penalised_log_likelihood
  if (all(is.na(penalised))) return(NA_real_)
  max(penalised, na.rm = TRUE)
}

# THE PRINTED VIEW, SEPARATED FROM PRINTING IT. Extracted so the shift can be
# tested: the arithmetic below is the load-bearing part and a `capture.output()`
# assertion on a formatted number is a poor way to hold it.
#
# BOTH LIKELIHOOD COLUMNS SHIFT BY THE SAME CONSTANT, the maximum `L_P`. That is
# what keeps `L - L_P = p` true row by row, so a reader can check the arithmetic
# on the face of the table -- which deducting each column's own maximum would
# break. `L_P` is then 0 at the top and negative below.
comparison_shown_table <- function(x) {
  best <- comparison_best_L_P(x)

  shown <- data.frame(
    L_P = x$table$penalised_log_likelihood - best,
    L = x$table$log_likelihood - best,
    p = x$table$penalty,
    k = x$table$k,
    row.names = rownames(x$table)
  )

  # THE STATUS COLUMN APPEARS ONLY WHEN IT SAYS SOMETHING. With everything
  # converged it reads `converged` all the way down and is wide enough to wrap
  # every row. It stays in `$table` either way.
  if (length(x$failures) > 0L) shown$status <- x$table$status
  shown
}

#' @rdname compare_models
#' @export
is_model_comparison <- function(x) inherits(x, comparison_class)

#' @rdname compare_models
#' @export
best_fit <- function(x) {
  if (!is_model_comparison(x)) {
    stop(sprintf("`%s` must be a `model_comparison`, as returned by `compare_models()`.",
                 deparse(substitute(x))), call. = FALSE)
  }
  best <- x$fits[[1L]]
  if (is.null(best)) {
    stop("No candidate converged, so there is no best fit. See `$failures` for why.",
         call. = FALSE)
  }
  best
}

# NO `logLik()` OR `AIC()` HERE EITHER, for the reason `fit()` gives: base R's
# `AIC()` penalises by the parameter count, and `p` equals that only where the
# weight is an indicator.

#' @rdname compare_models
#' @export
print.model_comparison <- function(x, ...) {
  candidates <- nrow(x$table)
  cat(sprintf("<model_comparison: %d candidate%s, %d converged>\n",
              candidates, if (candidates == 1L) "" else "s",
              candidates - length(x$failures)))

  cat_analysis_header(x)

  # BOTH LIKELIHOOD COLUMNS ARE SHIFTED BY THE SAME CONSTANT, the maximum `L_P`.
  # That is what keeps `L - L_P = p` true row by row, so a reader can check the
  # arithmetic on the face of the table -- which deducting each column's own
  # maximum would break. `L_P` is then 0 at the top and negative below.
  cat(sprintf("  max L_P         %s\n", format(comparison_best_L_P(x))))

  cat("\n")
  print(comparison_shown_table(x))

  # THE DIAGNOSIS IS SHOWN, NOT JUST THE STATUS. `fit()` puts real effort into
  # saying what went wrong and what to do about it, and a comparison that
  # reported only `not_identifiable` would throw all of it away.
  for (label in names(x$failures)) {
    cat(sprintf("\n%s: %s\n", label, x$failures[[label]]))
  }

  invisible(x)
}
