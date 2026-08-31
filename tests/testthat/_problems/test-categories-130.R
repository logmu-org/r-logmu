# Extracted from test-categories.R:130

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
status_levels <- function() {
  factor(c("act", "def", "pen"), levels = c("act", "def", "pen", "dep"))
}

# test -------------------------------------------------------------------------
expect_error(categories(.i$status, "act", "def", "pen", "depp", .source = status_levels()),
               "no category")
