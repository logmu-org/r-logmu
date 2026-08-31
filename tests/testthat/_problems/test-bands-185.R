# Extracted from test-bands.R:185

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(durations(.i$entry, 0, 11, by = 5), "divide")
