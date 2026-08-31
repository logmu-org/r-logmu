# Extracted from test-bands.R:165

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_equal(duration(.i$birth, 65, 95)$terms, age(65, 95)$terms)
