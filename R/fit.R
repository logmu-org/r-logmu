# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# fit -- maximum weighted log-likelihood for a proportional hazards model
# ============================================================================
#
# The whole iteration lives in the engine (`src/veil/FitLoop.hpp`), deliberately,
# so a front end in another language gets the loop rather than a
# reimplementation of it. Everything here turns a `model` and a dataset into the
# arguments `cpp_veil_fit_run` wants, works out Z, and reads the answer back.

fit_class <- "fit"

# The defaults of the numerical method, named once so the documentation and the
# signature cannot drift apart.
default_L_tolerance <- 1e-6
default_armijo <- 1e-4
default_max_iterations <- 25L
default_max_halvings <- 30L

#' Fit a mortality model by maximum weighted log-likelihood
#'
#' @description
#' Estimates the coefficients of a [model()] against experience data, by
#' Newton-Raphson on the weighted log-likelihood
#' \deqn{L^* = \Omega^{-1}(\mathrm{A}w\log\mu - \mathrm{E}w).}
#'
#' A fit returns a fit or it fails. A model that does not converge, or whose
#' information matrix cannot be factored, raises an error carrying the
#' diagnosis rather than returning a number nobody should use.
#'
#' @details
#' # The Z scale
#'
#' Results are reported on the \eqn{L = Z^{-1}L^*} scale, where
#' \eqn{Z = \mathrm{E}w^2/\mathrm{E}w} on a test mortality. \eqn{Z} is a change
#' of units, chosen so that one parameter costs about one unit of \eqn{L}, and
#' it must be the same for every model being compared -- otherwise the
#' penalised log-likelihoods are not on one scale and the comparison means
#' nothing.
#'
#' So `z` is given in one of two ways, and is never taken from the model's own
#' reference mortality:
#'
#' - as a number, which is what model selection needs, the caller having
#'   computed it once and passed it to every candidate; or
#' - as a test mortality, from which \eqn{Z} is computed over this same data,
#'   weight and include. [gompertz_mortality()] is a reasonable choice.
#'
#' `z` may be omitted only when the weight is absent or is an indicator, where
#' \eqn{w^2 = w} makes \eqn{Z} exactly one.
#'
#' # Convergence
#'
#' The iteration stops when the gain in \eqn{L} still available falls below
#' `L_tolerance`, which is a forward-looking test: it measures the distance left to
#' the optimum rather than the size of the last step, so it is not fooled by a
#' damped step or by a model whose maximum does not exist. Stopping within
#' \eqn{\epsilon} of the maximum leaves the coefficients within
#' \eqn{\sqrt{2\epsilon/I}} of theirs.
#'
#' @param exp_data The experience data.
#' @param model The [model()] to fit.
#' @param include An `include` restricting the population, or `NULL` for all of
#'   it.
#' @param weight A pronoun expression for the weight \eqn{w}, or `NULL` for
#'   lives. Weights must not be negative.
#' @param val_similarity,val_distance A pronoun expression for the second
#'   weighting factor, in either spelling. At most one may be given.
#' @param z The \eqn{Z} scale: a single positive number, or a test mortality
#'   from which to compute it. See Details.
#' @param settings A [settings()] object carrying `overdispersion` and
#'   `time_scale`.
#' @param overdispersion,time_scale Given directly, these override the
#'   `settings`.
#' @param start Starting coefficients, one per model term. Defaults to zeros.
#' @param L_tolerance The convergence tolerance, as a gain in \eqn{L}. One unit of
#'   \eqn{L} is about what one parameter costs.
#' @param armijo The fraction of the linearly predicted gain a damped step must
#'   actually deliver to be accepted. Must lie strictly between 0 and 0.5.
#' @param max_iterations The most Newton steps to take.
#' @param max_halvings The most times one step may be halved before the fit is
#'   abandoned.
#' @param threads The number of worker threads.
#' @param x,object A `fit`.
#' @param ... Ignored.
#' @returns
#' `fit()` returns a `fit`, whose fields are `beta`, `mortality` (the fitted
#' model as a `mortality`), `variance`, `log_likelihood`, `penalty`,
#' `penalised_log_likelihood`, `z`, `overdispersion`, `iterations`,
#' `evaluations`, `predicted_gain` and `model`.
#'
#' `is_fit()` returns a scalar `logical`.
#'
#' `coef()` returns the fitted coefficients as a named numeric vector, and
#' `vcov()` their variance-covariance matrix.
#' @examples
#' data <- exp_data(
#'   list(
#'     birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
#'     male      = c(TRUE, FALSE, TRUE, FALSE, TRUE),
#'     E2R_start = datey::datey(rep(2015, 5)),
#'     E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
#'     E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
#'   ),
#'   exp_start = datey::datey(2015),
#'   exp_end   = datey::datey(2020)
#' )
#'
#' m <- model(ref_mortality = gompertz_mortality(), covariates = covariates(level = 1))
#' fit(data, m, settings = settings(overdispersion = 1))
#' @name fit
NULL

#' @rdname fit
#' @export
fit <- function(exp_data,
                model,
                include = NULL,
                weight = NULL,
                val_similarity = NULL,
                val_distance = NULL,
                z = NULL,
                settings = NULL,
                overdispersion = NULL,
                time_scale = NULL,
                start = NULL,
                L_tolerance = default_L_tolerance,
                armijo = default_armijo,
                max_iterations = default_max_iterations,
                max_halvings = default_max_halvings,
                threads = cpp_veil_default_threads()) {

  ensure_is_exp_data(exp_data)
  if (missing(model)) stop("`model` is required.", call. = FALSE)
  ensure_is_model(model)

  resolved <- resolve_settings(settings, overdispersion, time_scale)

  # CAPTURED FROM THE CALLER'S FRAME. `missing()` rather than a value test, for
  # the reason `aev()` gives: `substitute()` on an unsupplied argument yields
  # its default, so a value test cannot tell `weight = NULL` from no weight.
  caller <- parent.frame()
  weight_ast <- aev_optional_ast(substitute(weight), missing(weight), caller)
  similarity_ast <- aev_optional_ast(substitute(val_similarity), missing(val_similarity), caller)
  distance_ast <- aev_optional_ast(substitute(val_distance), missing(val_distance), caller)
  z_ast <- aev_optional_ast(substitute(z), missing(z), caller)

  if (!is.null(include)) {
    if (is_includes(include)) {
      stop("`include` takes a single `include`. A fit is one model over one population.",
           call. = FALSE)
    }
    ensure_is_include(include)
  }

  terms <- model_term_asts(model)
  start <- fit_start(start, length(terms))

  columns <- exp_data_columns(exp_data)
  z_value <- fit_z_scale(z_ast, weight_ast, include, columns,
                         resolved$time_scale_clicks, threads)

  run <- cpp_veil_fit_run(
    model$ref_mortality,
    terms,
    weight_ast,
    similarity_ast,
    distance_ast,
    columns,
    resolved$time_scale_clicks,
    include,
    start,
    as.integer(max_iterations),
    as.double(L_tolerance),
    as.double(armijo),
    as.double(max_halvings),
    resolved$overdispersion,
    z_value,

    # THE ASSERTION TRAVELS, THE VERIFICATION HAPPENS IN THE ENGINE. R can prove a `disjoint()` term
    # does not vary with time and nothing more -- a column's type is unknown until the data arrives.
    # The engine walks the data once before compiling and refuses a claim that does not hold, which
    # is what earns the right to leave the off-diagonal integrals out of the block entirely.
    is_disjoint(model$covariates),
    as.integer(threads)
  )

  if (!identical(run$status, "converged")) {
    stop(fit_failure_message(run, model, L_tolerance), call. = FALSE)
  }

  new_fit(run, model, z_value, resolved$overdispersion)
}

# The starting coefficients. Zeros unless the user says otherwise -- an
# A/E-based warm start was considered and withdrawn, because the level is
# carried across age-based terms in a realistic model and a flat scaling is not
# what anybody fits.
fit_start <- function(start, terms) {
  if (is.null(start)) return(rep(0, terms))
  if (!is_pure_numeric(start) || length(start) != terms || !all(is.finite(start))) {
    stop(sprintf("`start` must be %d finite number%s, one for each model term.",
                 terms, if (terms == 1L) "" else "s"),
         call. = FALSE)
  }
  as.double(start)
}

# Z, by whichever of the two routes the caller asked for.
#
# NEVER FROM THE MODEL'S OWN REFERENCE MORTALITY. Candidates may differ in their
# reference, so a Z taken from it would move with the candidate and destroy the
# comparability it exists to provide.
fit_z_scale <- function(z_ast, weight_ast, include, columns, time_scale_clicks, threads) {

  if (is.null(z_ast)) {
    # `w^2 = w`, so `Z = Ew^2/Ew = 1` exactly and there is nothing to run. This
    # is the only case where omitting `z` is safe, and it covers lives-weighted
    # work entirely.
    if (is.null(weight_ast) || it_is_indicator(weight_ast)) return(1)

    stop("`z` is required when the weight is not an indicator. Give it as a single ",
         "number -- which is what comparing models needs, since they must share one Z ",
         "-- or as a test mortality such as `gompertz_mortality()`.",
         call. = FALSE)
  }

  # A number written into the call, which is the model-selection route.
  if (identical(z_ast$kind, "lit") && is.numeric(z_ast$value)) {
    value <- as.double(z_ast$value)
    if (!is_single_pure_finite_numeric(value) || value <= 0) {
      stop("`z` must be a single positive finite number.", call. = FALSE)
    }
    return(value)
  }

  # Otherwise a test mortality, and Z is `Ew^2/Ew` over the SAME data, weight
  # and include.
  #
  # PINNED AT `overdispersion = 1` AND READ AS `V / E`. The engine returns
  # `V = Omega Ew^2` for the Omega of the run that produced it, so at Omega = 1
  # that ratio is exactly `Ew^2/Ew`. Writing it as the general `V/(Omega E)`
  # invites substituting the fit's own Omega, which would scale L and p together
  # -- invisible in any ranking, and about 2.5x out in the calibration that says
  # a change of one is significant.
  ensure_mortality_leaves(z_ast)
  run <- cpp_veil_run(
    list(list(
      mortality      = z_ast,
      weight         = weight_ast,
      val_similarity = NULL,
      val_distance   = NULL,
      include        = include,
      overdispersion = 1
    )),
    columns,
    time_scale_clicks,
    FALSE,
    as.integer(threads)
  )

  result <- run$results[[1L]]
  value <- result$V / result$E
  if (!is_single_pure_finite_numeric(value) || value <= 0) {
    stop("The test mortality gave a Z of ", format(value), ", which cannot be used as a ",
         "scale. It is `Ew^2/Ew` over this data, so an empty or zero-weighted ",
         "population is the usual cause.", call. = FALSE)
  }
  value
}

# ---- failure ---------------------------------------------------------------
#
# FAILURE, NOT A WARNING. A warning is easily ignored, and a model whose
# information matrix is singular has no answer at all rather than one we did not
# reach: the coefficient along the null direction is not determined by the data,
# so any number returned would be arbitrary. Once refusing is the behaviour the
# message is the whole of the user's experience, so each one names what happened
# and what to do about it.

fit_failure_message <- function(run, model, L_tolerance) {
  switch(run$status,
    not_identifiable = fit_collinearity_message(run, model),
    not_finite = paste0(
      "The fit did not converge: the log-likelihood stopped being finite. ",
      "`Ew` overflows once a linear predictor reaches about 709, so a covariate on a ",
      "very large scale, or a `start` far from the data, is the usual cause."),
    did_not_converge = sprintf(
      paste0("The fit did not converge in %d iteration%s; the gain still available was %s ",
             "against an `L_tolerance` of %s. Raise `max_iterations`, or check that the model ",
             "is identifiable."),
      run$iterations, if (identical(run$iterations, 1L)) "" else "s",
      format(run$predicted_gain), format(L_tolerance)),
    step_collapsed = paste0(
      "The fit did not converge: a Newton step was halved to nothing without improving ",
      "the log-likelihood. That usually means the model is very nearly unidentifiable."),
    sprintf("The fit did not converge (%s).", run$status)
  )
}

# THE PIVOT INDEX NAMES WHERE THE DEPENDENCY WAS DETECTED, NOT THE CULPRIT --
# with terms 1 and 3 identical the failure is reported at 3, because 1 and 2 are
# still independent. So the message must never say "this term is bad"; it says
# which term could not be added, and, from the dependency the factorisation
# recovers, what it is a combination of.
fit_collinearity_message <- function(run, model) {
  labels <- model_term_labels(model)
  failed <- run$failed_parameter + 1L

  if (failed < 1L || failed > length(labels)) {
    return("The model is not identifiable: the information matrix is singular.")
  }

  head <- sprintf(paste0("The model is not identifiable. Term %d (%s) is a linear ",
                         "combination of the terms before it, so its coefficient is not ",
                         "determined by the data."),
                  failed, labels[[failed]])

  weights <- run$dependency
  earlier <- seq_len(min(length(weights), failed - 1L))
  if (length(earlier) < 1L) {
    return(paste0(head, " It is constant, or zero, over the included population."))
  }
  if (!all(is.finite(weights[earlier]))) return(head)

  parts <- sprintf("%s * [%d] %s", format(weights[earlier], digits = 4L),
                   earlier, labels[earlier])
  paste0(head, " It equals ", paste(parts, collapse = " + "), ".")
}

# ---- the result ------------------------------------------------------------

new_fit <- function(run, model, z, overdispersion) {
  labels <- model_term_labels(model)

  beta <- run$beta
  names(beta) <- labels

  # THE FITTED MODEL IS A `mortality`, and it is a FIELD rather than the object
  # itself, so the fit's convergence, variance and diagnostics stay separate
  # from the thing you hand on. Once the coefficients are numbers,
  # `log mu^ref + sum beta_j X_j` is a complete specification, which is what
  # gives layered fitting for nothing.
  fitted <- mortality_from_ast(fitted_mortality_ast(model, run$beta))

  structure(
    list(
      beta = beta,
      mortality = fitted,
      variance = unpack_triangle(run$variance, length(labels), labels),
      log_likelihood = run$log_likelihood,
      penalty = run$penalty,

      # L_P = L(beta_hat) - p, which is what models are compared on. Where the
      # weight is an indicator `p` is exactly the number of terms, and this is
      # then -AIC/2.
      penalised_log_likelihood = run$log_likelihood - run$penalty,

      # REPORTED, because it is the unit every one of these numbers is quoted
      # in. An unusual portfolio or a changed test mortality is then visible
      # rather than silent.
      z = z,
      overdispersion = overdispersion,

      iterations = run$iterations,
      evaluations = run$evaluations,
      predicted_gain = run$predicted_gain,
      model = model
    ),
    class = fit_class
  )
}

# `log mu^ref + sum_j beta_j X_j`, with the coefficients now numbers. Built as
# nodes rather than re-parsed, exactly as the covariate algebra builds its
# products: both sides are already the trees the front end produces.
fitted_mortality_ast <- function(model, beta) {
  ast <- model$ref_mortality
  terms <- model_term_asts(model)
  for (j in seq_along(terms)) {
    ast <- it_call("+", list(ast, it_call("*", list(it_lit(beta[[j]]), terms[[j]]))))
  }
  ast
}

# The engine packs symmetric matrices as an UPPER triangle, row-major, which is
# `packedTriangleIndex` and the only such convention in the package.
unpack_triangle <- function(packed, terms, labels) {
  out <- matrix(NA_real_, nrow = terms, ncol = terms, dimnames = list(labels, labels))
  if (length(packed) != terms * (terms + 1L) / 2L) return(out)

  at <- 1L
  for (row in seq_len(terms)) {
    for (column in row:terms) {
      out[row, column] <- packed[[at]]
      out[column, row] <- packed[[at]]
      at <- at + 1L
    }
  }
  out
}

#' @rdname fit
#' @export
is_fit <- function(x) inherits(x, fit_class)

#' @rdname fit
#' @export
coef.fit <- function(object, ...) object$beta

#' @rdname fit
#' @export
vcov.fit <- function(object, ...) object$variance

# NO `logLik()` METHOD, DELIBERATELY. It would make `AIC()` work and give the
# wrong answer: `AIC()` penalises by the number of parameters, where the right
# penalty here is `p = Z^-1 tr(J I^-1)`, which equals the term count only when
# the weight is an indicator. Compare models on `penalised_log_likelihood`.

#' @rdname fit
#' @export
print.fit <- function(x, ...) {
  cat(sprintf("<fit: %d term%s, %d iteration%s, %d walk%s of the data>\n",
              length(x$beta), if (length(x$beta) == 1L) "" else "s",
              x$iterations, if (identical(x$iterations, 1L)) "" else "s",
              x$evaluations, if (identical(x$evaluations, 1L)) "" else "s"))

  errors <- sqrt(diag(x$variance))
  table <- data.frame(estimate = unname(x$beta), std_error = unname(errors),
                      row.names = names(x$beta))
  print(table)

  cat(sprintf("\nlog-likelihood %s   penalty %s   penalised %s\n",
              format(x$log_likelihood), format(x$penalty),
              format(x$penalised_log_likelihood)))
  cat(sprintf("Z %s   overdispersion %s\n", format(x$z), format(x$overdispersion)))
  invisible(x)
}
