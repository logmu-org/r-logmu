# Extracted from test-aev.R:285

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
    pension   = c(5000, 12000, 30000, 8000, 15000),
    male      = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    E2R_start = datey::datey(rep(2015, 5)),
    E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
    E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)
pensions <- c(5000, 12000, 30000, 8000, 15000)
with_log_mu <- function(log_mu) {
  exp_data(
    list(
      birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
      pension   = pensions,
      log_mu    = log_mu,
      E2R_start = datey::datey(rep(2015, 5)),
      E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
      E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )
}
flat <- mortality_const(log_mu = -4)
basis <- settings(overdispersion = 2)

# test -------------------------------------------------------------------------
sentinel <- with_log_mu(rep(-800, length(pensions)))
b <- batch(.exp_data = sentinel, .settings = basis,
             ordinary  = aev(mortality = -4),
             underflow = aev(mortality = .i$log_mu),
             also_fine = aev(mortality = -5))
