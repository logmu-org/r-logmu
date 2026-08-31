# Extracted from test-categories.R:206

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
status_levels <- function() {
  factor(c("act", "def", "pen"), levels = c("act", "def", "pen", "dep"))
}

# test -------------------------------------------------------------------------
breakdown <- includes(
    categories(.i$status, "act", "def", "pen", "dep"),
    ages(65, 85, by = 10)
  )
