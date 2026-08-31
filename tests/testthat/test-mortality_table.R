# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

test_that("simple test of mortality_table", {

  get_test_q <- function(N_x, N_t) {
    noise <- function(x,t) {
      noise_x <- (x * 0.618033988749895) %% 1
      noise_xt <- ((noise_x + t) * 0.618033988749895) %% 1
      noise_xt - 0.5
    }
    gompertz <- function(x,t) -4 + 0.13 * x - 0.01 * t + 0.001 * noise(x,t)

    log_m <- outer(1:N_x, 1:N_t, gompertz)
    m <- exp(log_m)
    q <- -expm1(-m)
    q
  }

  N_x <- 4
  N_t <- 5

  q <- get_test_q(N_x = 4, N_t = 5)
  m <- -log1p(-q)
  mu_data <- c(
    # Left age col 1:
    m[1,1],
    m[2,1] - (m[3,1] - 2*m[2,1] + m[1,1])/24,
    m[3,1] - (m[4,1] - 2*m[3,1] + m[2,1])/24,
    m[4,1],

    # Col 2:
    m[1,2],
    m[2,2] - (m[3,3] - 2*m[2,2] + m[1,1])/24,
    m[3,2] - (m[4,3] - 2*m[3,2] + m[2,1])/24,
    m[4,2],

    # Col 3:
    m[1,3],
    m[2,3] - (m[3,4] - 2*m[2,3] + m[1,2])/24,
    m[3,3] - (m[4,4] - 2*m[3,3] + m[2,2])/24,
    m[4,3],

    # Col 4:
    m[1,4],
    m[2,4] - (m[3,5] - 2*m[2,4] + m[1,3])/24,
    m[3,4] - (m[4,5] - 2*m[3,4] + m[2,3])/24,
    m[4,4],

    # Right age col 5:
    m[1,5],
    m[2,5] - (m[3,5] - 2*m[2,5] + m[1,5])/24,
    m[3,5] - (m[4,5] - 2*m[3,5] + m[2,5])/24,
    m[4,5]
  )

  mu <- matrix(data = mu_data, nrow = N_x, ncol = N_t)
  log_mu <- log(mu)

  x0 <- datey::durationy(65.25)
  t0 <- datey::datey(2020.125)

  table_q <- mortality_table(x0, t0, q = q)
  table_mu <- mortality_table(x0, t0, mu = mu)
  table_log_mu <- mortality_table(x0, t0, log_mu = log_mu)

  strip <- function(x) {
    attr(x, "end_age") <- NULL
    x
  }

  expect_equal(strip(table_q), strip(table_log_mu))
  expect_equal(strip(table_mu), strip(table_log_mu))

  test_value_on_table <- function(mortality, .b, .t, expected) {
    .b <- datey::datey(.b)
    .t <- datey::datey(.t)
    .i <- list(birth = .b)
    actual <- log_mu(table_log_mu, .i, .t)
    expect_equal(actual, expected)
  }
  test_value <- function(.b, .t, expected) {
    test_value_on_table(table_q, .b, .t, expected)
    test_value_on_table(table_mu, .b, .t, expected)
    test_value_on_table(table_log_mu, .b, .t, expected)
  }

  b0 <- t0 - x0

  # Corners
  test_value(b0, t0, log_mu[1,1])
  test_value(b0 + 4, t0 + 4, log_mu[1,5])
  test_value(b0 - 3, t0, log_mu[4,1])
  test_value(b0 + 1, t0 + 4, log_mu[4,5])

  # Check x is correct when outside the time range
  test_value(b0 + (-10 - 2), t0 - 10, log_mu[3,1])
  test_value(b0 + (10 - 1), t0 + 10, log_mu[2,5])

  # Interpolation inside the range
  # x = 2.125, t = 3.25
  test_value(b0 + 0.875, t0 + 2.25,
    log_mu[3,3] * 0.125 * 0.75 +
    log_mu[4,4] * 0.125 * 0.25 +
    log_mu[2,3] * 0.875 * 0.75 +
    log_mu[3,4] * 0.875 * 0.25
  )

})
