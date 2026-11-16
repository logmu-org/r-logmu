# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# Reading naxp objects without depending on the naxp package
# ============================================================================
#
# logmu accepts a `naxp` as a category domain and a `naxp_encoded` vector as a
# column of codes, but it never parses a pattern or encodes text, so it does
# not need the naxp package itself. These helpers read the two objects directly,
# and they are the ONLY place that knows their representation (naxp 0.1.0):
#
#   naxp          the pattern string, classed "naxp"
#   naxp_encoded  integer or double codes 0..max, classed "naxp_encoded",
#                 with the naxp in attribute "naxp"
#
# If naxp's representation changes, this file changes and nothing else does,
# and the test that reads real naxp objects in `test-categorical.R` says so.
#
# Codes from a `naxp_encoded` vector were checked by naxp when it was made, so
# logmu trusts their range. That is why a naxp map's keys must arrive already
# encoded rather than as bare numbers: logmu cannot check a bare number against
# a pattern's maximum without the naxp package.

is_naxp <- function(x) inherits(x, "naxp")
is_naxp_encoded <- function(x) inherits(x, "naxp_encoded")

naxp_pattern <- function(naxp) as.character(unclass(naxp))

same_naxp <- function(a, b) identical(naxp_pattern(a), naxp_pattern(b))

# The naxp a `naxp_encoded` vector was encoded with.
naxp_of <- function(x) attr(x, "naxp", exact = TRUE)

# The bare codes, integer or double, without class, naxp or names.
naxp_codes <- function(x) {
  codes <- unclass(x)
  attributes(codes) <- NULL
  codes
}

new_naxp_encoded_codes <- function(codes, naxp) {
  attr(codes, "naxp") <- naxp
  class(codes) <- "naxp_encoded"
  codes
}
