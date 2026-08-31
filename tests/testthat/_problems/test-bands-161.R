# Extracted from test-bands.R:161

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_equal(duration(.i$entry, 0, 5)$terms, band(.t - .i$entry, 0, 5)$terms)
