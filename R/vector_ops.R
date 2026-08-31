# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

#' SIMD active lanes and tier
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
