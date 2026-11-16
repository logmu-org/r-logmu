# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Real naxp objects, from the naxp package, which is suggested rather than
# imported. A test that uses these skips itself where naxp is not installed.
test_naxp <- function(pattern) naxp::naxp_parse(pattern)

test_naxp_encoded <- function(codes, naxp) naxp::naxp_encoded(codes, naxp)
