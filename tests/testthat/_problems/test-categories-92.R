# Extracted from test-categories.R:92

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
status_levels <- function() {
  factor(c("act", "def", "pen"), levels = c("act", "def", "pen", "dep"))
}

# test -------------------------------------------------------------------------
data <- data.frame(status = status_levels(), scheme = c("A", "B", "A"))
expect_equal(names(categories(.i$status, .source = data)), c("act", "def", "pen", "dep"))
