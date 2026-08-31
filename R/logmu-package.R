# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' @title Actuarial mortality experience analysis and modelling
#' @description
#'
#' TODO: CHANGE THIS TO A USER SUMMARY OF THE PACKAGE
#'
#' NOTE THAT FOR CONINTUOUS TIME **logmu** USES DATEY PACKAGE
#'
#' The **logmu** package provides
#' high-performance actuarial mortality experience analysis and
#' model fitting, built on a flexible mortality framework incorporating
#' time-based covariates and arbitrary proportional hazards models and
#' with particular support for postcode-based socio-economic mortality
#' models.
#'
#' A/E analysis includes confidence intervals and residuals with visualisation.
#' Model fitting features include covariate clustering,
#' maximum likelihood fitting and model selection using AIC.
#'
#' Weighting #' (e.g. amounts vs lives),
#' probabilistic similarity weighting (e.g. down-weight older data) and
#' optionally time-based inclusion criteria (e.g. sub-setting experience by age) are supported throughout.
#'
#' Calculations use SIMD vectorisation and multi-threading for performance.
#'
#' @seealso
#' - TODO
#' - vignette("XXX") for a worked introduction.
#' - The series of articles starting with `vignette("Measures matter")`
#'   for the **logmu** theoretical framework.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom lifecycle deprecated
#' @useDynLib logmu, .registration = TRUE
## usethis namespace: end
NULL
