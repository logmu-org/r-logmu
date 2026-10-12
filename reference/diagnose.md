# Standard diagnostics for a fit, or for every candidate of a comparison

Runs the A/E diagnostics of a fitted model over the experience data it
was fitted to, broken down however you ask, and returns them alongside
the model's coefficient table and correlation matrix. Given a
`model_comparison`, it does so for **every** converged candidate, not
just the winner, and the whole set crosses the data once.

A `fit` keeps no data, which is why the data is passed again here rather
than remembered. What the fit *does* keep is what it was run with – its
include, weight, relevance, overdispersion and time scale – and those
are the defaults, so the diagnostics of a fit are by default over the
same population, weighted the same way, as the fit itself.

## Usage

``` r
diagnose(
  exp_data,
  x,
  breakdown = NULL,
  include,
  weight,
  val_relevance,
  val_distance,
  settings = NULL,
  overdispersion = NULL,
  time_scale = NULL,
  threads = cpp_veil_default_threads()
)

is_diagnosis(x)

# S3 method for class 'diagnosis'
print(x, ...)

# S3 method for class 'diagnosis'
plot(x, ...)
```

## Arguments

- exp_data:

  The experience data, ordinarily the data the model was fitted to.

- x:

  A `fit` or a `model_comparison`. A lone fit is one candidate, labelled
  by its model's `name` where it has one and `fit` otherwise.

- breakdown:

  An `include`, an
  [`includes()`](https://r-logmu.logmu.org/reference/includes.md), a
  list of either, or `NULL`. See *Breakdowns*.

- include:

  The population. Defaults to the one the fit was run over.

- weight, val_relevance, val_distance:

  Pronoun expressions, as in
  [`aev()`](https://r-logmu.logmu.org/reference/aev.md). Each defaults
  to what the fit was run with; give `NULL` to run without.

- settings, overdispersion, time_scale:

  As in [`aev()`](https://r-logmu.logmu.org/reference/aev.md). Each
  defaults to the fit's own.

- threads:

  The number of worker threads.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html), passed
  to [`autoplot.aev()`](https://r-logmu.logmu.org/reference/aev_plot.md)
  for every chart, except `title` and `subtitle`, which name the
  candidate and the breakdown.

## Value

`diagnose()` returns a `diagnosis`, whose fields are `candidates` (one
entry per converged candidate, best first, each holding `summary` – the
[summary()](https://r-logmu.logmu.org/reference/fit.md) of its fit – and
`aev`, a list of one
[`aev()`](https://r-logmu.logmu.org/reference/aev.md) per breakdown),
`breakdowns` (the breakdown labels), `labels` (the include and weight
the diagnostics were run with, as text), `overdispersion`, and
`pv_impact`, which is `NULL` until present-value impacts exist.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the chart
grid and returns `x` invisibly.

## Breakdowns

`breakdown` may be a single `include` or
[`includes()`](https://r-logmu.logmu.org/reference/includes.md), or a
list of them for several breakdowns side by side. Each gives one
[`aev()`](https://r-logmu.logmu.org/reference/aev.md) per candidate,
with one record per element. A list's names label the breakdowns;
unnamed ones are labelled by what was written. `NULL` gives one
ungrouped A/E per candidate.

## What it shows

[`print()`](https://rdrr.io/r/base/print.html) gives each candidate's
coefficient table and correlation matrix, then each of its A/E tables.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws a grid of
A/E charts, candidates down the page and breakdowns across it, each
drawn by
[`autoplot.aev()`](https://r-logmu.logmu.org/reference/aev_plot.md).

## Examples

``` r
data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948, 1952, 1938, 1960, 1943)),
    male      = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE),
    sex       = c("m", "f", "m", "f", "m", "m", "f", "m", "f"),
    E2R_start = datey::datey(rep(2015, 9)),
    E2R_end   = datey::datey(c(2020, 2018, 2020, 2019, 2020, 2020, 2020, 2017, 2020)),
    E2R_died  = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

reference <- gompertz_mortality(slope_t = 0)
comparison <- compare_models(
  data,
  list(level = model(reference, covariates(level = 1)),
       sex   = model(reference, covariates(male = .i$male, female = !.i$male))),
  overdispersion = 1
)

d <- diagnose(data, comparison,
              breakdown = list(sex = categories(.i$sex, "m", "f"), age = ages(55, 85, by = 15)))
d
#> <diagnosis: 2 candidates, 2 breakdowns>
#>   include         all records
#>   weight          lives
#>   overdispersion  1
#> 
#> ---- level
#> 
#>       estimate std_error conf_95pc
#> level    2.026    0.4472    0.8765
#> 
#> A/E by sex
#> <aev[2]>
#>   group name A        E        V       A/E 95% conf  dev resid
#> 1   sex    m 3 1.789144 1.789144 1.6767791 1.465296  0.8243406
#> 2   sex    f 2 3.210883 3.210883 0.6228816 1.093795 -0.7267536
#> 
#> A/E by age
#> <aev[2]>
#>   group  name A        E        V       A/E  95% conf  dev resid
#> 1   age 55-70 2 1.112445 1.112445 1.7978412 1.8582698  0.7558023
#> 2   age 70-85 3 3.887582 3.887582 0.7716879 0.9940502 -0.4691626
#> 
#> ---- sex
#> 
#>        estimate std_error conf_95pc
#> male      2.542    0.5773     1.132
#> female    1.553    0.7070     1.386
#> 
#> correlation
#>        male
#> female 0.00
#> 
#> A/E by sex
#> <aev[2]>
#>   group name A        E        V       A/E 95% conf     dev resid
#> 1   sex    m 3 3.000007 3.000007 0.9999978 1.131584 -3.820109e-06
#> 2   sex    f 2 2.000708 2.000708 0.9996463 1.385659 -5.003728e-04
#> 
#> A/E by age
#> <aev[2]>
#>   group  name A        E        V       A/E 95% conf  dev resid
#> 1   age 55-70 2 1.635357 1.635357 1.2229743 1.532646  0.2754234
#> 2   age 70-85 3 3.365357 3.365357 0.8914359 1.068397 -0.2029370
```
