# Extracted from test-aev-frame.R:176

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
data <- exp_data(
    list(
      birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
      pension   = c(5000, 12000, 30000, 8000, 15000),
      E2R_start = datey::datey(rep(2015, 5)),
      E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
      E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
    ),
    exp_start = datey::datey(2015),
    exp_end   = datey::datey(2020)
  )
result <- aev(data,
                mortality      = mortality_const(log_mu = -4),
                breakdown      = bands(.i$pension, thresholds = c(10000, 20000)),
                overdispersion = 1)
frame <- as.data.frame(result)
