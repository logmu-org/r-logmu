# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# Tests buffer allocation: how few time-vector stores a block needs at once, once liveness has said
# where each one stops being read.
#
# EVERY TEST ASSERTS THE VALUES AND THE COUNT TOGETHER, and neither alone would do. Allocation
# cannot change an answer, so a value-only test passes for a pass that does nothing at all; and a
# count-only test passes for a pass that hands two live vectors the same buffer, which is the one
# way this can go wrong and the way that produces plausible wrong numbers rather than a failure.
#
# THE ORACLE IS THE ANALYTIC INTEGRAL, as in test-veil_time.R. Exposures are whole quarters, so
# there is no short final interval and no half-click rounding to allow for, and the midpoint rule is
# exact for an integrand that is constant or linear in time.
#
# `log(exp(...))` IS THE CHAIN, and it is written that way on purpose: nothing in veil rewrites it,
# so it lowers to a real chain of vector instructions, while the value it computes is the age it
# started from to within a couple of units in the last place. That gives a chain of any length an
# answer that can be written down.

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
  E2R_start = datey_clicks(start_clicks),
  E2R_end   = datey_clicks(end_clicks),
  E2R_died  = c(TRUE, FALSE, TRUE)
)

# Ages at the two ends of the exposure, which is what a durationy integrand runs between.
age_from <- (start_clicks - unclass(cols$birth)) / clicks_per_year
age_to <- (end_clicks - unclass(cols$birth)) / clicks_per_year
exposure_years <- (end_clicks - start_clicks) / clicks_per_year

# The integral of age over the exposure, factored so the two squares are never subtracted from each
# other -- see the note in test-veil_time.R about what the obvious spelling costs.
integral_of_age <- (age_to - age_from) * (age_to + age_from) / 2

integrate_of <- function(expr) cpp_veil_integrate(expr, cols, quarter_scale, NULL)

test_that("a chain of vector operations needs two buffers however long it is", {
  # Each step is read by the next one and by nothing after it, so the whole chain ping-pongs between
  # a pair of buffers. That the count does not grow with the length is the point: the virtual form
  # names one vector per step.
  short <- integrate_of(it_ast(~ log(exp(log(exp(.x))))))
  long <- integrate_of(it_ast(~ log(exp(log(exp(log(exp(log(exp(.x))))))))))

  expect_equal(short$vector_operand_count, 5L)
  expect_equal(long$vector_operand_count, 9L)

  # TWO AND NOT ONE, stated as an equality so it pins the rule rather than merely the saving. A
  # buffer is released only once the instruction that last read it has finished, so a result never
  # lands on top of an argument it is still reading. Every vector op veil has today reads slot i to
  # write slot i and would survive the overlap; the next one to read a slot other than the one it
  # writes would not, and the values here would not notice.
  expect_equal(short$buffer_count, 2L)
  expect_equal(long$buffer_count, 2L)

  expect_equal(short$integral, integral_of_age, tolerance = 1e-12)
  expect_equal(long$integral, integral_of_age, tolerance = 1e-12)
})

test_that("two arms that are live at once do not share a buffer", {
  # `(.x + 1) + (.x + 2)` holds the age and both arms at the same time, so three stores are needed
  # where the chain above needed two. Sharing has already made the two mentions of `.x` one operand.
  res <- integrate_of(it_ast(~ (.x + 1) + (.x + 2)))

  expect_equal(res$vector_operand_count, 4L)
  expect_equal(res$buffer_count, 3L)

  # The integrand is 2 * age + 3, which is linear, so the midpoint rule gives it exactly.
  expect_equal(res$integral, 2 * integral_of_age + 3 * exposure_years, tolerance = 1e-12)
})

test_that("allocation leaves the answers alone", {
  # A sweep over shapes that allocate differently -- no vector at all, one, a chain, two live arms --
  # against answers worked out on paper. Nothing here should be sensitive to the layout, which is
  # exactly the claim being made.
  constant <- integrate_of(it_ast(~ 1))
  expect_equal(constant$buffer_count, 1L)
  expect_equal(constant$integral, exposure_years, tolerance = 1e-12)

  plain <- integrate_of(it_ast(~ .x))
  expect_equal(plain$buffer_count, 1L)
  expect_equal(plain$integral, integral_of_age, tolerance = 1e-12)

  # A per-individual scalar mixed into a time-varying expression stays a scalar, so it needs no
  # buffer of its own and the count is the same as the plain age.
  weighted <- integrate_of(it_ast(~ .x * .i$amount))
  expect_equal(weighted$buffer_count, 1L)
  expect_equal(weighted$integral, cols$amount * integral_of_age, tolerance = 1e-12)
})

test_that("lowering emits nothing that nothing reads", {
  # AN INVARIANT OF LOWERING RATHER THAN A LAW: it walks down from the roots and emits an
  # instruction only because something asked for its value, so a block as lowered holds no dead
  # code. Nothing removes dead code today because there has never been any, and this is the
  # assertion that will report the day that changes -- which the first TAC-level rewrite, fusing two
  # reductions over one grid or letting an if/else arm omit its instructions, will do on purpose.
  expect_equal(integrate_of(it_ast(~ .x))$dead_instruction_count, 0L)
  expect_equal(integrate_of(it_ast(~ log(exp(.x)) * .i$amount))$dead_instruction_count, 0L)

  res <- cpp_veil_aev(it_obj(mortality_const(log_mu = -3.2)), it_ast(~ .i$amount), cols, quarter_scale, NULL, no_overdispersion, 1L)
  expect_equal(res$dead_instruction_count, 0L)
})

test_that("an AEV reuses the buffer A is finished with", {
  # A's integrand is read by `died_value` and by nothing else, so the store it was broadcast into is
  # free by the time E's integrand needs one. Two vector operands, one buffer.
  log_mu_value <- -3.2
  res <- cpp_veil_aev(it_obj(mortality_const(log_mu = log_mu_value)), NULL, cols, quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(res$vector_operand_count, 2L)
  expect_equal(res$buffer_count, 1L)

  expect_equal(res$contributions$A, as.double(cols$E2R_died), tolerance = 1e-12)
  expect_equal(res$contributions$E, exp(log_mu_value) * exposure_years, tolerance = 1e-12)
  expect_equal(res$contributions$V, exp(log_mu_value) * exposure_years, tolerance = 1e-12)
})

test_that("a real mortality table with a weight that varies over time still fits three buffers", {
  # The heaviest shape the recipe currently produces: log mu over the grid, its exponential, the
  # weight over the grid, and the products for E and V. Seven vectors named, three needed at once.
  log_mu_value <- -3.2
  tbl <- mortality_table(
    x0 = 55, t0 = 2005,
    log_mu = matrix(log_mu_value, nrow = 40L, ncol = 20L)
  )

  # `.x * 1` rather than `.x`, because squaring a durationy for V is arithmetic the datey rules do
  # not admit; multiplying by one makes it a plain number of years.
  res <- cpp_veil_aev(it_obj(tbl), it_ast(~ .x * 1), cols, quarter_scale, NULL, no_overdispersion, 1L)

  expect_equal(res$vector_operand_count, 7L)
  expect_equal(res$buffer_count, 3L)

  # Constant log mu and a weight linear in time, so E's integrand is linear and exact. A is the age
  # at death, which is the exposure end for the individuals who died.
  expect_equal(res$contributions$A, ifelse(cols$E2R_died, age_to, 0), tolerance = 1e-12)
  expect_equal(res$contributions$E, exp(log_mu_value) * integral_of_age, tolerance = 1e-12)
})
