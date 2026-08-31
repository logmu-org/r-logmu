# Extracted from test-veil_aev.R:372

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
clicks_per_year <- 534360L
quarter <- clicks_per_year %/% 4L
datey_clicks <- function(clicks) {
  structure(as.integer(clicks), class = class(datey::datey(2010)))
}
birth_years <- c(1940, 1945, 1950)
start_clicks <- c(2010, 2010, 2010) * clicks_per_year
end_clicks <- start_clicks + c(12L, 8L, 4L) * quarter
cols <- list(
  birth     = datey::datey(birth_years),
  amount    = c(1000, 2500, 400),
  male      = c(TRUE, FALSE, TRUE),
  E2R_start = datey_clicks(start_clicks),
  E2R_end   = datey_clicks(end_clicks),
  E2R_died  = c(TRUE, FALSE, TRUE)
)
exposure_years <- (end_clicks - start_clicks) / clicks_per_year
died <- cols$E2R_died
log_mu_value <- -3.2
mu_value <- exp(log_mu_value)
constant_mortality <- mortality_const(log_mu = log_mu_value)
aev <- function(mortality = it_obj(constant_mortality), weight = NULL, include = NULL,
                time_scale = quarter_scale, columns = cols,
                overdispersion = no_overdispersion, threads = 1L) {
  cpp_veil_aev(mortality, weight, columns, time_scale, include, overdispersion, threads)
}

# test -------------------------------------------------------------------------
res <- aev(mortality = list(kind = "lit", value = NaN))
expect_true(is.nan(res$E))
expect_true(is.nan(res$V))
expect_true(all(is.nan(res$contributions$E)))
