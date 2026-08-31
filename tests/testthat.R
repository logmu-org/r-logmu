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

# TEMPORARY, 2026-08-31: name each file as it starts, flushed.
#
# A Windows CI runner dies during the suite with no testthat output at all --
# stdout is buffered and the buffer goes with the process -- so there is nothing
# to say which file it reached. This subclass keeps the check reporter's
# behaviour and adds one flushed line per file, which survives a crash.
#
# Remove once the failure is located.
LoudCheck <- R6::R6Class(
  "LoudCheck",
  inherit = testthat::CheckReporter,
  public = list(
    start_file = function(filename) {
      cat("[file]", filename, "
")
      flush(stdout())
      super$start_file(filename)
    }
  )
)

# TEMPORARY, 2026-08-31: run ONLY the vector-function file.
#
# The Windows runner dies inside it, but at a DIFFERENT operation each run --
# `pow` at N = 4 one time, `exp` at N = 1 the next. A wandering crash point is
# memory corruption, not a bad instruction, and the damage may well be done by
# something that ran earlier: the engine and thread-pool tests come before this
# file alphabetically.
#
# If this file passes alone, it is a victim and the corruption is upstream.
# If it still dies, the fault really is here.
#
# Restore to `test_check("logmu", reporter = LoudCheck$new())` afterwards.
test_check("logmu", filter = "vec_fns", reporter = LoudCheck$new())
