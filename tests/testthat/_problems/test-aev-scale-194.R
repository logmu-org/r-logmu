# Extracted from test-aev-scale.R:194

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
scale_aev <- function() {
  aev <- create_aev(A = c(1218, 1075, 1001), E = rep(1000, 3), V = c(876, 59, 2656))
  names(aev) <- c("a", "b", "c")
  aev
}
params <- function(p, panel = 1L) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplot_build(p)$layout$panel_params[[panel]]
}

# test -------------------------------------------------------------------------
aev <- create_aev(A = c(1100, 12), E = c(1000, 4), V = c(2500, 6))
expect_identical(as.character(aev_status(aev)), c("ok", "ok"))
