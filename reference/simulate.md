# Simulate experience data from a mortality

Draws a death time for every life in a population under a known
`mortality`, and returns the result as an
[`exp_data()`](https://r-logmu.logmu.org/reference/exp_and_val_data.md)
with the same shape as real data: `E2R_start`, `E2R_end` and `E2R_died`,
and every other column the population came with. A fit to the result can
then be judged against the mortality that made it.

## Usage

``` r
# S3 method for class 'mortality'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  lives,
  exp_start,
  exp_end,
  time_scale = default_time_scale,
  ...
)
```

## Arguments

- object:

  The `mortality` to simulate from: any mortality, including a fitted
  model's `$mortality` and a
  [`mortality()`](https://r-logmu.logmu.org/reference/mortality.md)
  expression switching on the population's columns.

- nsim:

  How many datasets to draw. One returns an `exp_data`; more returns a
  list of them.

- seed:

  Passed to [`set.seed()`](https://rdrr.io/r/base/Random.html) first
  when not `NULL`.

- lives:

  The population. See *The population*.

- exp_start, exp_end:

  The experience period, as `datey`s or years.

- time_scale:

  The step of the hazard grid, as in
  [`settings()`](https://r-logmu.logmu.org/reference/settings.md).

- ...:

  Ignored.

## Value

An
[`exp_data()`](https://r-logmu.logmu.org/reference/exp_and_val_data.md),
or a list of `nsim` of them.

## The population

`lives` is a data frame (or list of columns) with a `birth` column of
`datey`s, whatever covariate columns the mortality and any later
analysis need, and two optional `datey` columns:

- `E2R_start` – when each life's exposure begins: retirement, joining,
  the first date the record is reliable. Absent, every life starts at
  `exp_start`. A life is alive at its start by construction, which is
  the left-truncation the E2R framework assumes.

- `E2R_cap` – the date a life would leave observation for a reason other
  than death: withdrawal, transfer out, a record closed. `NA` means no
  such date. Absent, no life has one.

Each life is observed from `max(E2R_start, exp_start)` to the earliest
of its death, its `E2R_cap` and `exp_end`. A death before the cap is a
death and the cap plays no part; a life still alive at the cap is
censored there. Lives with no observation at all are dropped, with a
message saying how many. `E2R_cap` is dropped from the result, so that
simulated data has exactly the columns real data has. Names beginning
`E2R_` are reserved and cannot be read by a pronoun expression, so it
could not have leaked into a calculation either way.

## How deaths are drawn

Each life gets a target \\-\log U\\ for uniform \\U\\, and dies where
its cumulative hazard \\\int\mu\\dt\\ from its start first reaches it.
The hazard is accumulated on the same grid the engine integrates on:
from the life's own start in steps of `time_scale`, \\\mu\\ read at each
step's midpoint, a shorter final step to the cap, and the death placed
within the crossing step at the constant hazard of that step. So a fit
to simulated data is unbiased against the mortality that made it, up to
sampling noise.

A death is recorded to the day it fell on, at
[`datey::mid_day()`](https://r-datey.logmu.org/reference/datey.html),
exactly as a real death is. Deaths are independent across lives with no
extra variation, so the data is Poisson-like: an overdispersion of 1 is
the truth for it.

## Reproducibility

The targets are drawn with R's own generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) before the call, or
the `seed` argument, makes the result repeatable.

## Examples

``` r
set.seed(1)
lives <- data.frame(
  birth   = datey::datey(2015 - sample(60:85, 500, replace = TRUE)),
  sex     = sample(c("male", "female"), 500, replace = TRUE),
  pension = round(rlnorm(500, log(12000), 0.7))
)

male <- gompertz_mortality(log_mu_0 = -3.6)
female <- gompertz_mortality(log_mu_0 = -4.0)
truth <- mortality(if (.i$sex == "male") male else female)

data <- simulate(truth, lives = lives,
                 exp_start = datey::datey(2015), exp_end = datey::datey(2020))
data
#> # Exp data: 500 × 6 (5 yr)
#> #   Start:  2015-01-01.0
#> #   End:    2020-01-01.0
#>           birth sex    pension    E2R_start      E2R_end E2R_died
#>         <datey> <chr>    <dbl>      <datey>      <datey> <lgl>   
#>  1 1931-01-01.0 female   15426 2015-01-01.0 2015-04-25.5 TRUE    
#>  2 1952-01-01.0 male     10470 2015-01-01.0 2015-11-28.5 TRUE    
#>  3 1949-01-01.0 female    9759 2015-01-01.0 2020-01-01.0 FALSE   
#>  4 1955-01-01.0 female   16989 2015-01-01.0 2020-01-01.0 FALSE   
#>  5 1954-01-01.0 female   16850 2015-01-01.0 2020-01-01.0 FALSE   
#>  6 1933-01-01.0 male     12159 2015-01-01.0 2020-01-01.0 FALSE   
#>  7 1945-01-01.0 male     18714 2015-01-01.0 2020-01-01.0 FALSE   
#>  8 1942-01-01.0 male     20349 2015-01-01.0 2019-03-08.5 TRUE    
#>  9 1938-01-01.0 male     21508 2015-01-01.0 2019-07-11.5 TRUE    
#> 10 1937-01-01.0 female   23593 2015-01-01.0 2020-01-01.0 FALSE   
#> # ℹ 490 more rows

# The data is neutral against the mortality that made it, within noise.
aev(data, mortality = truth, overdispersion = 1)
#> <aev[1]>
#>    A        E        V      A/E  95% conf dev resid
#> 1 73 66.88376 66.88376 1.091446 0.2396557 0.7368819
```
