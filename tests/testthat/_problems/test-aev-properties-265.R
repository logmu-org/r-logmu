# Extracted from test-aev-properties.R:265

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
labelled <- create_aev(A = c(10, 20), E = c(1, 2), V = c(4, 5))
names(labelled) <- c("65-70", "70-75")
group_names(labelled) <- "age"
bare <- create_aev(A = c(30, 40), E = c(3, 4), V = c(6, 7))
expect_identical(names(labelled + bare), c("65-70", "70-75"))
