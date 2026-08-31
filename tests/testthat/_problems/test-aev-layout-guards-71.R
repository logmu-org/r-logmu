# Extracted from test-aev-layout-guards.R:71

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
guarded_aev <- function() {
  aev <- create_aev(A = c(1100, 970, 1050), E = rep(1000, 3), V = c(2500, 1600, 900))
  names(aev) <- c("a", "b", "c")
  aev
}
grouped_guarded_aev <- function() {
  aev <- guarded_aev()
  group_names(aev) <- c("one", "one", "two")
  aev
}
built <- function(aev = guarded_aev(), ...) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplotGrob(autoplot(aev, ...))
}
renaming <- function(table, from, to) {
  table$layout$name <- sub(paste0("^", from), to, table$layout$name)
  table
}

# test -------------------------------------------------------------------------
expect_no_error(aev_require_layout(theme_aev(), "the theme"))
expect_true(is.na(theme_aev()$legend.background$fill))
