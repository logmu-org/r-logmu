# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# model -- the right-hand side of a proportional hazards mortality model
# ============================================================================
#
#     log mu = log mu^ref + beta_1 X^(1) + beta_2 X^(2) + ...
#
# A MODEL IS THE RIGHT-HAND SIDE AND NOTHING ELSE. The experience data, the
# include, the weight and the convergence tolerance are arguments of `fit()`,
# not parts of the model, because every candidate in a selection exercise shares
# them. Keeping them out is what makes two models comparable by construction
# rather than by a rule in the documentation.

model_class <- "model"

#' A proportional hazards mortality model
#'
#' @description
#' A `model` is the right-hand side of
#' \deqn{\log\mu = \log\mu^{\mathrm{ref}} + \sum_j \beta_j X^{(j)},}
#' that is a reference mortality and the covariates whose coefficients
#' [fit()] estimates. It carries nothing else: the experience data, the
#' `include`, the weight and the tolerance all belong to the fit, because every
#' candidate model in a selection exercise shares them.
#'
#' A fitted model is itself a `mortality`, so it can be handed back as the
#' `ref_mortality` of another model.
#'
#' @param ref_mortality The reference mortality: a `mortality`, or a pronoun
#'   expression for \eqn{\log\mu^{\mathrm{ref}}}. Anything held fixed goes here
#'   -- a reference mortality is already an offset, so fixed effects need no
#'   feature of their own.
#' @param covariates A [covariates()] object giving the \eqn{X} terms. The
#'   coefficients belong to the fit; nothing here names or supplies one.
#' @param x An object.
#' @param ... Ignored.
#' @returns
#' `model()` returns a `model`.
#'
#' `is_model()` returns a scalar `logical`.
#' @examples
#' age_shape <- variable(.x)
#' model(
#'   ref_mortality = gompertz_mortality(),
#'   covariates    = covariates(male = .i$male, female = !.i$male) * age_shape
#' )
#' @name model
NULL

#' @rdname model
#' @export
model <- function(ref_mortality, covariates) {

  if (missing(ref_mortality)) stop("`ref_mortality` is required.", call. = FALSE)
  if (missing(covariates)) stop("`covariates` is required.", call. = FALSE)

  # CAPTURED UNEVALUATED, in the caller's frame, exactly as `aev()` captures its
  # mortality: this may be a `mortality` object, a `~` formula or a bare pronoun
  # expression, and only the parser can tell which.
  ref_ast <- it_capture(substitute(ref_mortality), parent.frame())
  ensure_mortality_leaves(ref_ast)

  ensure_is_covariates(covariates)
  if (length(covariates) < 1L) {
    stop("A `model` needs at least one covariate; there is nothing to fit.", call. = FALSE)
  }

  structure(list(ref_mortality = ref_ast, covariates = covariates), class = model_class)
}

#' @rdname model
#' @export
is_model <- function(x) inherits(x, "model")

ensure_is_model <- function(x) {
  if (!is_model(x)) {
    stop(sprintf("`%s` must be a `model`, as built by `model()`.",
                 deparse(substitute(x))), call. = FALSE)
  }
  x
}

# The X terms as the engine wants them: the bare ASTs, in order.
model_term_asts <- function(x) {
  lapply(unclass(x$covariates), function(term) term$ast)
}

# One label per coefficient. A name where the user gave one, and the expression
# itself where they did not -- naming covariates is not usual practice, and the
# expression is what they wrote, so it identifies the term better than `X1`
# would.
model_term_labels <- function(x) {
  # UNNAMED: `covariates` names its own items, and a named result would put those
  # names inside the `dimnames` of the variance matrix.
  unname(vapply(unclass(x$covariates), function(term) {
    name <- covariate_name_of(term)
    if (nzchar(name)) name else it_deparse(term$ast)
  }, character(1L)))
}

#' @rdname model
#' @export
print.model <- function(x, ...) {
  cat(sprintf("<model: %d term%s%s>\n",
              length(x$covariates), if (length(x$covariates) == 1L) "" else "s",
              if (is_disjoint(x$covariates)) ", disjoint" else ""))
  cat(sprintf("  reference: %s\n", it_deparse(x$ref_mortality)))
  labels <- model_term_labels(x)
  for (i in seq_along(labels)) {
    cat(sprintf("  [%d] %s\n", i, labels[[i]]))
  }
  invisible(x)
}
