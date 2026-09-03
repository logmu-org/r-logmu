# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' SIMD active lanes and tier
#'
#' @description
#'
#' `vec_active_tier()` names the SIMD kernel tier the package resolved when its
#' shared library was loaded, and `vec_active_lanes()` gives the number of
#' `double` lanes that tier works on.
#'
#' The tier is chosen from the instruction sets the CPU reports: `"avx512"`
#' (8 lanes), `"avx2"` (4 lanes), or `"baseline"` (2 lanes) where neither of
#' those is available.
#'
#' @section Forcing a lower tier:
#'
#' Setting the environment variable `LOGMU_TIER` to `"baseline"`, `"avx2"` or
#' `"avx512"` before the package is loaded places a ceiling on the tier that may
#' be selected. It is there for comparing tiers in a benchmark, and for
#' reproducing a fault that only one of them shows.
#'
#' The ceiling can only ever lower the tier, never raise it. A CPU that does not
#' support the tier asked for still falls back to the best one it does support,
#' because the processor feature checks apply exactly as they otherwise would.
#'
#' An unrecognised value places no ceiling at all, and nothing is reported, so
#' call `vec_active_tier()` to confirm which tier is in use rather than assuming
#' the request was honoured.
#'
#' The variable is read once, when the shared library is loaded. Changing it
#' later has no effect for the rest of the session.
#'
#' @name vec_active
NULL

#' @rdname vec_active
#' @returns A scalar `integer` representing the number of active SIMD lanes.
#' @export
vec_active_lanes <- function() cpp_vec_active_lanes()
#' @rdname vec_active
#' @returns A scalar `character` describing the active SIMD tier.
#' @export
vec_active_tier <- function() cpp_vec_active_tier()

#' Vector ops
#' @param x A vector of double.
#' @param y A vector of double.
#' @param min,max Vectors of double for the `clamp` operation.
#' @returns A vector of `double`.
#' @name vec_ops
NULL

#' @rdname vec_ops
#' @export
vec_neg <- function(x) cpp_vec_neg(x)
#' @rdname vec_ops
#' @export
vec_exp <- function(x) cpp_vec_exp(x)
#' @rdname vec_ops
#' @export
vec_expm1 <- function(x) cpp_vec_expm1(x)
#' @rdname vec_ops
#' @export
vec_log <- function(x) cpp_vec_log(x)
#' @rdname vec_ops
#' @export
vec_log1p <- function(x) cpp_vec_log1p(x)
#' @rdname vec_ops
#' @export
vec_m_from_q <- function(x) cpp_vec_m_from_q(x)

#' @rdname vec_ops
#' @export
vec_add <- function(x, y) cpp_vec_add(x, y)
#' @rdname vec_ops
#' @export
vec_sub <- function(x, y) cpp_vec_sub(x, y)
#' @rdname vec_ops
#' @export
vec_mul <- function(x, y) cpp_vec_mul(x, y)
#' @rdname vec_ops
#' @export
vec_div <- function(x, y) cpp_vec_div(x, y)
#' @rdname vec_ops
#' @export
vec_pow <- function(x, y) cpp_vec_pow(x, y)
#' @rdname vec_ops
#' @export
vec_min <- function(x, y) cpp_vec_min(x, y)
#' @rdname vec_ops
#' @export
vec_max <- function(x, y) cpp_vec_max(x, y)
#' @rdname vec_ops
#' @export
vec_clamp <- function(x, min, max) cpp_vec_clamp(x, min, max)
