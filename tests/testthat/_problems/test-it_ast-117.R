# Extracted from test-it_ast.R:117

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_error(it_build_ast(quote(NA)), "is missing")
expect_error(it_build_ast(quote(NA_real_)), "is missing")
expect_error(it_build_ast(quote(NA_integer_)), "is missing")
expect_error(it_build_ast(quote(NA_character_)), "is missing")
