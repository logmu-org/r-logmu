# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# ============================================================================
# Reading naxp objects without depending on the naxp package
# ============================================================================
#
# logmu accepts a `naxp` as a domain and a `naxp_encoded` vector as a column of
# codes. It reads the two objects directly rather than importing naxp, and
# these helpers are the ONLY place that knows their representation (naxp 0.1.0
# and 0.2.0 alike):
#
#   naxp          the pattern string, classed "naxp"
#   naxp_encoded  integer or double codes 0..max, classed "naxp_encoded",
#                 with the naxp in attribute "naxp"
#
# If naxp's representation changes, this file changes and nothing else does,
# and the test that reads real naxp objects in `test-categorical_encoded.R`
# says so.
#
# Only `encode()` and `decode()` need naxp itself, to parse text against a
# pattern, and they hand the work to it here. naxp is suggested, so they say
# so plainly where it is not installed.
#
# Codes from a `naxp_encoded` vector were checked by naxp when it was made, so
# logmu trusts their range. That is why a naxp map's keys must arrive already
# encoded rather than as bare numbers: logmu cannot check a bare number against
# a pattern's maximum without the naxp package.
#
# `is.na()` on a `naxp_encoded` vector is naxp's own method, from naxp 0.2.0.
# logmu never registers methods for naxp's classes: a postcode vector must not
# behave differently depending on whether logmu happens to be loaded. logmu
# tests for code 0 directly.

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

ensure_naxp_is_installed <- function() {
  if (!requireNamespace("naxp", quietly = TRUE)) {
    stop("Encoding or decoding with a naxp needs the naxp package, which is not installed.",
         call. = FALSE)
  }
}

# Text, with NA where missing, as codes of `naxp`. NA and text the pattern
# does not accept both encode to 0.
naxp_encode_text <- function(text, naxp) {
  ensure_naxp_is_installed()
  naxp::naxp_encode(text, naxp)
}

# The text of each code, with NA for code 0.
naxp_decode_codes <- function(x) {
  ensure_naxp_is_installed()
  text <- naxp::naxp_decode(x)
  names(text) <- NULL
  text
}
