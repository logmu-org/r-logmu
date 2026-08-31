# Extracted from test-it_ast.R:107

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(it_build_ast(quote(log(-1))), "is NaN whatever the data holds")
expect_error(it_build_ast(quote(0 / 0)), "is NaN whatever the data holds")
expect_error(it_build_ast(quote(NaN)), "is NaN whatever the data holds")
expect_error(it_build_ast(quote(.i$pension + NaN)), "`NaN` is NaN")
