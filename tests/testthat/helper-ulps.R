# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# THE AGREED TOLERANCE BETWEEN PATHS, for ELEMENTARY results: one exp, one log, one table lookup.
# Two paths that compile the same source may still differ by an ulp or two, through a multiply and
# add the compiler fused in one place and not the other, or a library choosing its algorithm by
# register width. Aggregates such as A, E and L accumulate those differences over many terms and are
# compared with a relative tolerance instead.

# The spacing of doubles at each |v|: 2^(exponent - 52), and the smallest subnormal at zero.
ulp_at <- function(v) {
  magnitude <- abs(v)
  pmax(2^(floor(log2(magnitude)) - 52), 2^-1074)
}

# How many ulps apart, measured at the larger of the two. Equal values, including equal infinities,
# are 0 apart, and two NaNs are treated as agreeing; anything else involving a NaN or an infinity is
# infinitely far apart.
ulps_apart <- function(x, y) {
  apart <- abs(x - y) / ulp_at(pmax(abs(x), abs(y)))
  apart[!is.na(x) & !is.na(y) & x == y] <- 0
  apart[is.nan(x) & is.nan(y)] <- 0
  apart[is.na(apart)] <- Inf
  apart
}

# The label comes from `substitute()` rather than rlang, which is not one of the package's declared
# dependencies, so `R CMD check` would flag it.
expect_within_ulps <- function(object, expected, ulps = 2) {
  label <- paste(deparse(substitute(object)), collapse = "")
  stopifnot(length(object) == length(expected))
  apart <- ulps_apart(object, expected)
  worst <- if (length(apart)) which.max(apart) else 0L
  testthat::expect(
    all(apart <= ulps),
    sprintf("`%s` differs by up to %s ulps (at element %d: %s against %s); allowed %d.",
            label, format(apart[worst]), worst,
            format(object[worst], digits = 17), format(expected[worst], digits = 17), ulps)
  )
  invisible(object)
}
