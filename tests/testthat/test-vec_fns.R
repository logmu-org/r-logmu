# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

create_unit_vector <- function(N, D = 0) ((1:N + D) * 1.6180339887) %% 1
create_signed_vector <- function(N, D = 0) create_unit_vector(N, D) * 2 - 1

# TEMPORARY, 2026-08-31: announce every operation, flushed.
#
# A Windows CI runner dies inside this file with no R error and no testthat
# output at all -- the process simply stops. These lines survive that, so the
# log names the operation and the vector length that killed it.
#
# Remove once the failure is located.
say <- function(what, N) {
  cat("  [vec]", what, "N =", N, "
")
  flush(stdout())
}

# TEMPORARY, 2026-08-31: how does it die, not just where.
#
# A Windows CI runner terminates silently inside `vec_pow` at N = 4. Running the
# same call in a SUBPROCESS lets us read the exit code, and on Windows that code
# names the fault: 0xC000001D (3221225501) is an illegal instruction, meaning
# the binary used something the CPU does not have; 0xC0000005 (3221225477) is an
# access violation and 0xC00000FD (3221225725) a stack overflow, meaning a
# memory bug instead. The distinction decides where to look next.
#
# Remove once the failure is understood.
test_that("how vec_pow dies, if it dies", {

  script <- tempfile(fileext = ".R")
  writeLines(c(
    'library(logmu)',
    'ux <- ((1:4 + 0) * 1.6180339887) %% 1',
    'y  <- (((1:4 + 1) * 1.6180339887) %% 1) * 2 - 1',
    'invisible(vec_pow(ux, y))',
    'cat("PROBE SURVIVED VV\n")',
    'invisible(vec_pow(ux, 0.1234))',
    'cat("PROBE SURVIVED VS\n")',
    'invisible(vec_pow(0.2345, y))',
    'cat("PROBE SURVIVED SV\n")'
  ), script)

  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), shQuote(script),
    stdout = TRUE, stderr = TRUE
  ))

  status <- attr(out, "status")
  if (is.null(status)) status <- 0L

  cat("[probe] exit status:", status, "
")
  cat("[probe] output:", paste(out, collapse = " | "), "
")
  flush(stdout())

  succeed()
})

test_that("vector unary functions", {

  test_functions <- function(N)
  {
    x <- create_signed_vector(N)
    exp_x <- exp(x)

    say("neg", N)
    expect_identical(vec_neg(x), -x)
    say("exp", N)
    expect_equal(vec_exp(x), exp(x))
    say("expm1", N)
    expect_equal(vec_expm1(x), expm1(x))
    say("log", N)
    expect_equal(vec_log(exp_x), log(exp_x))
    say("log1p", N)
    expect_equal(vec_log1p(x), log1p(x))
  }

  test_functions(1)
  test_functions(3)
  test_functions(4)
  test_functions(5)
  test_functions(7)
  test_functions(8)
  test_functions(31)
  test_functions(32)
  test_functions(33)
  test_functions(1003)
})

test_that("vector binary functions", {

  test_functions <- function(N)
  {
    x <- create_signed_vector(N, 0)
    y <- create_signed_vector(N, 1)
    ux <- create_unit_vector(N, 0)

    say("add", N)
    expect_identical(vec_add(x, y), x + y)
    expect_identical(vec_add(x, 0.1234), x + 0.1234)
    expect_identical(vec_add(-0.2345, y), -0.2345 + y)

    say("sub", N)
    expect_identical(vec_sub(x, y), x - y)
    expect_identical(vec_sub(x, 0.1234), x - 0.1234)
    expect_identical(vec_sub(-0.2345, y), -0.2345 - y)

    say("mul", N)
    expect_identical(vec_mul(x, y), x * y)
    expect_identical(vec_mul(x, 0.1234), x * 0.1234)
    expect_identical(vec_mul(-0.2345, y), -0.2345 * y)

    say("div", N)
    expect_identical(vec_div(x, y), x / y)
    expect_identical(vec_div(x, 0.1234), x / 0.1234)
    expect_identical(vec_div(-0.2345, y), -0.2345 / y)

    # Temperamental
    say("pow", N)
    expect_equal(vec_pow(ux, y), ux ^ y)
    expect_equal(vec_pow(ux, 0.1234), ux ^ 0.1234)
    expect_equal(vec_pow(0.2345, y), 0.2345 ^ y)

    say("min", N)
    expect_identical(vec_min(x, y), pmin(x,y))
    expect_identical(vec_min(x, 0.1234), pmin(x, 0.1234))
    expect_identical(vec_min(-0.2345, y), pmin(-0.2345, y))

    say("max", N)
    expect_identical(vec_max(x, y), pmax(x,y))
    expect_identical(vec_max(x, 0.1234), pmax(x, 0.1234))
    expect_identical(vec_max(-0.2345, y), pmax(-0.2345, y))

    say("clamp", N)
    expect_identical(vec_clamp(x, -0.123, +0.234), pmin(pmax(x,-0.123),+0.234))
  }

  test_functions(1)
  test_functions(3)
  test_functions(4)
  test_functions(5)
  test_functions(7)
  test_functions(8)
  test_functions(31)
  test_functions(32)
  test_functions(33)
  test_functions(1003)
})
