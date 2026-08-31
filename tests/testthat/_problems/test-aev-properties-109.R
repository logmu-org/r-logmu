# Extracted from test-aev-properties.R:109

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_equal(create_aev(A = 1e-200, E = 1e200, V = 1)$deviance_residual,
               -1e200 * sqrt(2))
