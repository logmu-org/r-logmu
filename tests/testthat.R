# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/testing-design.html#sec-tests-files-overview
# * https://testthat.r-lib.org/articles/special-files.html

library(testthat)
library(logmu)

# WHICH SIMD TIER THIS MACHINE CHOSE, printed before anything runs.
#
# If a kernel executes an instruction the CPU does not implement, the process
# dies on the spot: no R error, no testthat output, just a log that stops after
# `test_check()`. That happened on a CI runner on 2026-08-31 and there was
# nothing in the output to say which tier had been selected. One line here puts
# it in every test log, including a truncated one.
cat("logmu SIMD tier:", vec_active_tier(), "with", vec_active_lanes(), "lanes\n")

test_check("logmu")
