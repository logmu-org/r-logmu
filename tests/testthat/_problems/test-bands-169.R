# Extracted from test-bands.R:169

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
d <- duration(min(.i$entry, .i$retirement), 0, 5)
