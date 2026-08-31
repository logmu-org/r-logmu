# Extracted from test-batch.R:238

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940)),
    pension   = c(5000, 12000, 30000, 8000),
    male      = c(TRUE, FALSE, TRUE, FALSE),
    E2R_start = datey::datey(rep(2015, 4)),
    E2R_end   = datey::datey(c(2020, 2020, 2018, 2020)),
    E2R_died  = c(FALSE, FALSE, TRUE, FALSE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)
other_data <- exp_data(
  list(
    birth     = datey::datey(1950),
    E2R_start = datey::datey(2015),
    E2R_end   = datey::datey(2020),
    E2R_died  = FALSE
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)
light <- mortality_const(log_mu = -4.5)
heavy <- mortality_const(log_mu = -4.0)
ageing <- mortality(-10 + 0.1 * .x)

# test -------------------------------------------------------------------------
expect_error(
    batch(.overdispersion = 1, nowhere = aev(mortality = light)),
    "`nowhere` has no `exp_data`"
  )
