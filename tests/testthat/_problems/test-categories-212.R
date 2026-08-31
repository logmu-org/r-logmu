# Extracted from test-categories.R:212

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
status_levels <- function() {
  factor(c("act", "def", "pen"), levels = c("act", "def", "pen", "dep"))
}

# test -------------------------------------------------------------------------
both <- category(.i$status, "act") & age(65, 95)
