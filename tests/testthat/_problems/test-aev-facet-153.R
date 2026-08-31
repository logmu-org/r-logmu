# Extracted from test-aev-facet.R:153

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
facet_aev_aev <- function() {
  aev <- create_aev(A = c(1100, 970, 1010), E = rep(1000, 3), V = c(2500, 1387, 900))
  names(aev) <- c("[65, 75)", "[75, 85)", "[85, 95)")
  aev
}
grouped_aev <- function() {
  aev <- create_aev(
    A = c(1100, 970, 1010, 1040, 990),
    E = rep(1000, 5),
    V = c(2500, 1387, 900, 1600, 1200)
  )
  names(aev) <- c(
    "[2010, 2015)", "[65, 75)", "[75, 85)", "[85, 95)", "[1940, 1960)"
  )
  group_names(aev) <- c("Period", "Age", "Age", "Age", "Cohort")
  aev
}
assembled <- function(p) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  ggplot2::ggplotGrob(p)
}
rows_of <- function(table, prefix) {
  table$layout$t[startsWith(table$layout$name, prefix)]
}
cell <- function(table, prefix, row) {
  table$grobs[[which(startsWith(table$layout$name, prefix) & table$layout$t == row)[[1]]]]
}

# test -------------------------------------------------------------------------
table <- assembled(autoplot(facet_aev_aev()))
