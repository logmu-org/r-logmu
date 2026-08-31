# Extracted from test-aev-properties.R:302

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(create_aev(A = 1, E = 0, V = 0), "A cannot be non-zero")
degenerate <- create_aev_unchecked(A = 1, E = 0, V = 0)
