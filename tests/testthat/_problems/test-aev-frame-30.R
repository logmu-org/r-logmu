# Extracted from test-aev-frame.R:30

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
frame_columns <- c("name", "group", "A", "E", "V",
                   "A_over_E", "log_A_over_E_stddev", "deviance_residual")
labelled_aev <- function() {
  aev <- create_aev(A = c(1100, 40, 0), E = c(1000, 50, 20), V = c(2500, 125, 40))
  names(aev) <- c("65-70", "70-75", "75-80")
  group_names(aev) <- "age"
  aev
}

# test -------------------------------------------------------------------------
frame <- as.data.frame(labelled_aev())
