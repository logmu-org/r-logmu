# Extracted from test-val_similarity.R:121

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955)),
    pension   = c(5000, 12000, 30000),
    male      = c(TRUE, FALSE, TRUE),
    fraction  = c(0.2, 0.5, 0.9),
    E2R_start = datey::datey(rep(2015, 3)),
    E2R_end   = datey::datey(c(2020, 2020, 2018)),
    E2R_died  = c(FALSE, FALSE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)
flat <- mortality_const(log_mu = -4)
parts <- function(x) c(x$A, x$E, x$V)

# test -------------------------------------------------------------------------
expect_error(aev(data, mortality = flat, overdispersion = 1, val_similarity = 2),
               "always above 1")
