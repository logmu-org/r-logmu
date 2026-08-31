# Extracted from test-aev-properties.R:253

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "logmu", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
male <- create_aev(A = c(10, 20), E = c(1, 2), V = c(4, 5))
names(male) <- c("65-70", "70-75")
group_names(male) <- "age"
female <- create_aev(A = c(30, 40), E = c(3, 4), V = c(6, 7))
names(female) <- c("65-70", "70-75")
group_names(female) <- "age"
both <- male + female
expect_identical(names(both), c("65-70", "70-75"))
expect_identical(group_names(both), rep("age", 2))
expect_equal(both$A, c(40, 60))
other <- create_aev(A = c(30, 40), E = c(3, 4), V = c(6, 7))
names(other) <- c("2000-2005", "2005-2010")
expect_error(male + other, "`names` differ")
same_names <- create_aev(A = c(30, 40), E = c(3, 4), V = c(6, 7))
names(same_names) <- c("65-70", "70-75")
group_names(same_names) <- "duration since entry"
expect_error(male + same_names, "`group_names` differ")
