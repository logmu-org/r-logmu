# LogMu mortality experience analysis and model fitting
#
# This file is licensed to you under the Apache Licence 2.0.
#
# Copyright (c) Tim Gordon

# `simulate()` draws experience data from a known mortality. THE ORACLES:
#
#   * with the uniform draw pinned, the death instant is a line of arithmetic
#     -- the constant-hazard crossing within the step -- so the walk is checked
#     to the click and not to a tolerance;
#   * with the draw free, the data must be neutral against the mortality that
#     made it: an A/E on the truth sits at one within its own interval, which
#     is the whole point of simulating.

years <- function(x) datey::datey(x)
clicks_of <- function(x) as.double(unclass(x))
click_year <- as.double(unclass(datey::durationy(1)))

population <- function(n, seed = 1) {
  set.seed(seed)
  data.frame(
    birth   = years(2015 - sample(60:85, n, replace = TRUE)),
    sex     = sample(c("male", "female"), n, replace = TRUE),
    pension = round(exp(rnorm(n, log(12000), 0.7)))
  )
}

male <- gompertz_mortality(log_mu_0 = -3.6, slope_t = 0)
female <- gompertz_mortality(log_mu_0 = -4.0, slope_t = 0)
by_sex <- mortality(if (.i$sex == "male") male else female)

test_that("simulate() returns experience data shaped exactly like real data", {
  lives <- population(200)
  data <- simulate(by_sex, seed = 7, lives = lives, exp_start = 2015, exp_end = 2020)

  expect_true(is_exp_data(data))
  expect_identical(names(unclass(data)),
                   c("birth", "sex", "pension", "E2R_start", "E2R_end", "E2R_died"))
  expect_identical(exp_start(data), years(2015))
  expect_identical(exp_end(data), years(2020))
  expect_identical(length(unclass(data)$birth), 200L)
  expect_identical(unclass(data)$sex, lives$sex)
  expect_true(all(unclass(data)$E2R_start == years(2015)))
  expect_true(is.logical(unclass(data)$E2R_died))

  # A survivor runs to the end of the period; a death is recorded at midday.
  survived <- !unclass(data)$E2R_died
  expect_true(all(unclass(data)$E2R_end[survived] == years(2020)))
  expect_true(all(datey::is_mid_day(unclass(data)$E2R_end[!survived])))
  expect_true(any(!survived))
})

test_that("the draw is reproducible and `nsim` gives a list", {
  lives <- population(100)
  one <- simulate(by_sex, seed = 3, lives = lives, exp_start = 2015, exp_end = 2020)
  two <- simulate(by_sex, seed = 3, lives = lives, exp_start = 2015, exp_end = 2020)
  expect_identical(one, two)

  set.seed(3)
  three <- simulate(by_sex, lives = lives, exp_start = 2015, exp_end = 2020)
  expect_identical(one, three)

  several <- simulate(by_sex, nsim = 2, seed = 3, lives = lives, exp_start = 2015, exp_end = 2020)
  expect_length(several, 2L)
  expect_identical(several[[1L]], one)
  expect_false(identical(several[[1L]], several[[2L]]))
})

test_that("with the draw pinned, a constant hazard dies exactly where the arithmetic says", {
  # Target -log(U) = 0.25 at mu = 0.1 is 2.5 years after the start, whatever
  # the step, because a constant hazard accumulates exactly on any grid. From
  # 2015-01-01.0 that is 2017-07-02.0, recorded at that day's midday.
  local_mocked_bindings(runif = function(n, ...) rep(exp(-0.25), n), .package = "stats")
  lives <- data.frame(birth = years(c(1950, 1940)))
  steady <- mortality_const(log_mu = log(0.1))

  for (scale in c(1, 1 / 4, 1 / 12)) {
    data <- simulate(steady, lives = lives, exp_start = 2015, exp_end = 2020, time_scale = scale)
    expect_identical(unclass(data)$E2R_died, c(TRUE, TRUE))
    expect_identical(unclass(data)$E2R_end, datey::mid_day(years(c(2017.5, 2017.5))))
  }

  # Not reached before the period ends: a survivor, censored at `exp_end`.
  short <- simulate(steady, lives = lives, exp_start = 2015, exp_end = 2017)
  expect_identical(unclass(short)$E2R_died, c(FALSE, FALSE))
  expect_true(all(unclass(short)$E2R_end == years(2017)))
})

test_that("with the draw pinned, a death whose midday would overrun its cap is recorded at the cap", {
  # The instant is 2017-07-02.2, its midday is .5, and the cap is .4: the
  # death is a death -- it came before the cap -- and it is recorded at the cap
  # rather than half a day past it.
  mu <- 0.1
  target <- mu * (2.5 - 0.3 / 365)
  local_mocked_bindings(runif = function(n, ...) rep(exp(-target), n), .package = "stats")
  lives <- data.frame(birth = years(c(1950, 1940)))
  # Subsetting rather than recycling, which would strip the `datey` class.
  lives$E2R_cap <- datey::datey(2017, 7, 2, day_fraction = 0.4)[c(1L, 1L)]

  data <- simulate(mortality_const(log_mu = log(mu)), lives = lives, exp_start = 2015, exp_end = 2020)
  expect_identical(unclass(data)$E2R_died, c(TRUE, TRUE))
  expect_identical(unclass(data)$E2R_end, lives$E2R_cap)

  # Without the cap the same death is at that day's midday, after the cap.
  free <- simulate(mortality_const(log_mu = log(mu)), lives = lives[, "birth", drop = FALSE],
                   exp_start = 2015, exp_end = 2020)
  expect_identical(unclass(free)$E2R_end, datey::mid_day(years(c(2017.5, 2017.5))))
  expect_true(all(unclass(free)$E2R_end > lives$E2R_cap))
})

test_that("with the draw pinned, a death on the evening of its start day keeps its instant", {
  # A life entering at 2015-01-01.75 and dying a tenth of a day later would be
  # recorded at that day's midday, BEFORE it entered, and have an empty E2R.
  # Such a death keeps its instant instead.
  mu <- 0.1
  days_in_2015 <- 365
  target <- mu * 0.1 / days_in_2015
  local_mocked_bindings(runif = function(n, ...) rep(exp(-target), n), .package = "stats")
  lives <- data.frame(birth = years(1950), E2R_start = datey::datey(2015, 1, 1, day_fraction = 0.75))
  data <- simulate(mortality_const(log_mu = log(mu)), lives = lives, exp_start = 2015, exp_end = 2020)

  instant <- clicks_of(lives$E2R_start) + (target / mu) * click_year
  expect_identical(unclass(data)$E2R_died, TRUE)
  expect_identical(unclass(data)$E2R_end, datey::datey(round(instant) / click_year))
  expect_false(datey::is_mid_day(unclass(data)$E2R_end))
  expect_true(unclass(data)$E2R_end > unclass(data)$E2R_start)
})

test_that("with the draw pinned, the hazard is read at each step's midpoint", {
  # A hazard that rises through the year, one-year steps, and a target reached
  # inside the first step. The rule is constant hazard at the step's MIDPOINT,
  # so the fraction of the step is target / (mu(mid) * 1 year), and reading mu
  # at the start of the step instead would put the death 5% later.
  local_mocked_bindings(runif = function(n, ...) rep(exp(-0.05), n), .package = "stats")
  rising <- mortality(-3 + 0.1 * (.t - 2015))
  lives <- data.frame(birth = years(1950))
  data <- simulate(rising, lives = lives, exp_start = 2015, exp_end = 2020, time_scale = 1)

  mu_mid <- exp(-3 + 0.1 * 0.5)
  instant <- clicks_of(years(2015)) + (0.05 / mu_mid) * click_year
  expected <- datey::mid_day(datey::datey(instant / click_year))
  expect_identical(unclass(data)$E2R_end, expected)

  # And it is the midpoint rather than the start: that answer differs.
  mu_start <- exp(-3)
  wrong <- datey::mid_day(datey::datey((clicks_of(years(2015)) + (0.05 / mu_start) * click_year) / click_year))
  expect_false(identical(unclass(data)$E2R_end, wrong))
})

test_that("the final short step and its midpoint follow the engine's rule", {
  # A full step's midpoint is exact; only a short final step can have an odd
  # sum, which is halved to the nearer EVEN click so the half-click residuals
  # cancel across a portfolio.
  expect_identical(simulate_midpoint(0, 133590), 66795)
  expect_identical(simulate_midpoint(c(1, 1, 2), c(4, 6, 6)), c(2, 4, 4))
})

test_that("`E2R_start` and `E2R_cap` bound each life's window", {
  lives <- population(300)
  lives$E2R_start <- years(2015 + c(rep(0, 100), rep(1.5, 100), rep(-1, 100)))
  lives$E2R_cap <- years(c(rep(NA, 100), rep(2019, 100), rep(2018.25, 100)))
  data <- simulate(by_sex, seed = 5, lives = lives, exp_start = 2015, exp_end = 2020)
  u <- unclass(data)

  expect_false("E2R_cap" %in% names(u))
  # A start before the period is the period's start; one inside it is kept.
  expect_true(all(u$E2R_start[1:100] == years(2015)))
  expect_true(all(u$E2R_start[101:200] == years(2016.5)))
  expect_true(all(u$E2R_start[201:300] == years(2015)))
  # No E2R runs past its cap, and a survivor is censored exactly there.
  expect_true(all(u$E2R_end[101:200] <= years(2019)))
  expect_true(all(u$E2R_end[101:200][!u$E2R_died[101:200]] == years(2019)))
  expect_true(all(u$E2R_end[201:300] <= years(2018.25)))
  expect_true(all(u$E2R_end[1:100][!u$E2R_died[1:100]] == years(2020)))
  # A death before the cap is a death, and it is recorded at midday -- unless
  # that midday would overrun a cap that is not itself a day boundary
  # (2018.25 is 2018-04-02.25), when the cap is kept.
  expect_true(any(u$E2R_died[101:200]))
  died_end <- u$E2R_end[u$E2R_died]
  expect_true(all(datey::is_mid_day(died_end) | died_end == years(2018.25)))
  expect_true(all(datey::is_mid_day(u$E2R_end[1:200][u$E2R_died[1:200]])))
})

test_that("lives with no window are dropped and counted", {
  lives <- population(50)
  lives$E2R_start <- years(c(rep(2015, 40), rep(2021, 10)))
  expect_message(
    data <- simulate(by_sex, seed = 1, lives = lives, exp_start = 2015, exp_end = 2020),
    "10 of 50 lives had no exposure"
  )
  expect_identical(length(unclass(data)$birth), 40L)
  expect_identical(unclass(data)$pension, lives$pension[1:40])
})

test_that("the data is neutral against the mortality that made it", {
  # THE POINT OF THE THING. An A/E on the truth sits at one within its own
  # 95% interval (four standard deviations here, so the test is not flaky), on
  # the whole population and by sex -- which also exercises the per-individual
  # `if` in the R-side evaluator over a whole population at once.
  lives <- population(6000)
  data <- simulate(by_sex, seed = 11, lives = lives, exp_start = 2015, exp_end = 2020,
                   time_scale = 1 / 12)

  whole <- aev(data, mortality = by_sex, overdispersion = 1)
  expect_lt(abs(log(whole$A / whole$E)), 4 * whole$log_A_over_E_stddev)

  split <- aev(data, mortality = by_sex, overdispersion = 1,
               breakdown = categories(.i$sex, "male", "female"))
  expect_true(all(abs(log(split$A / split$E)) < 4 * split$log_A_over_E_stddev))

  # And a fit of a level on the truth finds nothing to fit.
  f <- fit(data, model(by_sex, covariates(level = 1)), overdispersion = 1)
  expect_lt(abs(unname(f$beta)), 4 * sqrt(f$variance[1, 1]))
})

test_that("a mortality table simulates too", {
  # The table lookup is vectorised over births as well as times, which is what
  # a population walked click by click needs.
  ages <- 0:120
  log_mu <- outer(-9 + 0.09 * ages, 0:9 * -0.01, `+`)
  table <- mortality_table(x0 = 0, t0 = 2013, log_mu = log_mu)
  lives <- population(3000)
  data <- simulate(table, seed = 2, lives = lives, exp_start = 2015, exp_end = 2020)

  expect_true(is_exp_data(data))
  whole <- aev(data, mortality = table, overdispersion = 1)
  expect_lt(abs(log(whole$A / whole$E)), 4 * whole$log_A_over_E_stddev)
})

test_that("bad populations are refused with a reason", {
  lives <- population(10)
  expect_error(simulate(by_sex, exp_start = 2015, exp_end = 2020), "`lives` is required")
  expect_error(simulate(by_sex, lives = lives), "`exp_start` and `exp_end` are required")
  expect_error(simulate(by_sex, lives = lives, exp_start = 2020, exp_end = 2015),
               "`exp_start` must be before")
  expect_error(simulate(by_sex, nsim = 0, lives = lives, exp_start = 2015, exp_end = 2020),
               "`nsim`")
  expect_error(simulate(by_sex, lives = lives[, "sex", drop = FALSE], exp_start = 2015,
                        exp_end = 2020),
               "`birth` column")
  expect_error(simulate(by_sex, lives = list(years(1950), 1), exp_start = 2015, exp_end = 2020),
               "named list")

  done <- simulate(by_sex, seed = 1, lives = lives, exp_start = 2015, exp_end = 2020)
  expect_error(simulate(by_sex, lives = done, exp_start = 2015, exp_end = 2020),
               "already an `exp_data`")

  lives$E2R_died <- FALSE
  expect_error(simulate(by_sex, lives = lives, exp_start = 2015, exp_end = 2020),
               "must not carry `E2R_died`")
  lives$E2R_died <- NULL
  lives$E2R_cap <- 2017
  expect_error(simulate(by_sex, lives = lives, exp_start = 2015, exp_end = 2020),
               "must be a `datey` column")
  lives$E2R_cap <- NULL
  ragged <- as.list(lives)
  ragged$extra <- 1:3
  expect_error(simulate(by_sex, lives = ragged, exp_start = 2015, exp_end = 2020),
               "same length")

  expect_error(simulate(by_sex, lives = lives[1:10, 1:3], exp_start = 2015, exp_end = 2020,
                        time_scale = 0.3),
               "`time_scale`")
})
