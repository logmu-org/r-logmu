# Extracted from test-categories.R:242

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
status_levels <- function() {
  factor(c("act", "def", "pen"), levels = c("act", "def", "pen", "dep"))
}
exp_status <- function(status) {
  exp_data(
    list(
      birth     = datey::datey(c(1945, 1950, 1955)),
      status    = status,
      E2R_start = datey::datey(rep(2015, 3)),
      E2R_end   = datey::datey(rep(2020, 3)),
      E2R_died  = rep(FALSE, 3)
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )
}

# test -------------------------------------------------------------------------
data <- exp_status(factor(c("act", "def", "act"), levels = c("act", "def", "pen")))
flat <- mortality_const(log_mu = -4)
basis <- settings(overdispersion = 1)
total <- unclass(aev(data, mortality = flat, settings = basis))$E
groups <- unclass(aev(data, mortality = flat, settings = basis,
                        breakdown = categories(.i$status, .source = data)))$E
