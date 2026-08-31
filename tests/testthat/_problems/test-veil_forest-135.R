# Extracted from test-veil_forest.R:135

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
cols <- list(
  age    = c(30, 45, 60),
  weight = c(1.5, 2.5, 0.5),
  height = c(1.6, 1.8, 1.7),
  birth  = datey::datey(c(1960, 1965, 1975))
)
eval_multi <- function(...) cpp_veil_eval_multi(list(...), cols)

# test -------------------------------------------------------------------------
res <- cpp_veil_eval_multi(list(list(kind = "lit", value = NaN),
                                  list(kind = "lit", value = 1)), cols)
expect_true(all(is.nan(res$values[[1]])))
