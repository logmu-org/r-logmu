# Compare mortality models by penalised log-likelihood

Fits several candidate
[`model()`](https://r-logmu.logmu.org/reference/model.md)s against one
set of experience data and ranks them on penalised log-likelihood,
\$\$L\_{\mathrm{P}} = L(\hat\beta) - p, \qquad p =
Z^{-1}\mathrm{tr}(\mathbf{J}\mathbf{I}^{-1}).\$\$ Larger is better. This
is Akaike's information criterion without the factor of \\-2\\,
generalised to an arbitrary log-likelihood weight.

## Usage

``` r
compare_models(
  exp_data,
  models,
  include = NULL,
  weight = NULL,
  val_relevance = NULL,
  val_distance = NULL,
  test_mortality = default_mortality(),
  Z = NULL,
  settings = NULL,
  overdispersion = NULL,
  time_scale = NULL,
  L_tolerance = default_L_tolerance,
  armijo = default_armijo,
  max_iterations = default_max_iterations,
  max_halvings = default_max_halvings,
  threads = cpp_veil_default_threads()
)

is_model_comparison(x)

best_fit(x)

# S3 method for class 'model_comparison'
summary(object, ...)

# S3 method for class 'summary.model_comparison'
print(x, ...)

# S3 method for class 'model_comparison'
print(x, ...)
```

## Arguments

- exp_data:

  The experience data.

- models:

  A list of [`model()`](https://r-logmu.logmu.org/reference/model.md)s,
  ideally named. Names appear in the ranking and index `fits`. A
  candidate the list does not name takes its model's own `name`, and one
  with neither is labelled by position.

- include:

  An `include` restricting the population, or `NULL` for all of it.

- weight:

  A pronoun expression for the weight \\w\\, or `NULL` for lives.
  Weights must not be negative.

- val_relevance, val_distance:

  A pronoun expression for the second weighting factor, in either
  spelling. At most one may be given.

- test_mortality:

  The mortality on which the \\Z\\ scale is measured, once, for every
  candidate. See [`fit()`](https://r-logmu.logmu.org/reference/fit.md).

- Z:

  The \\Z\\ scale as a single positive number, given instead of a
  `test_mortality` rather than as well as one.

- settings:

  A [`settings()`](https://r-logmu.logmu.org/reference/settings.md)
  object carrying `overdispersion` and `time_scale`.

- overdispersion, time_scale:

  Given directly, these override the `settings`.

- L_tolerance:

  The convergence tolerance, as a gain in \\L\\.

- armijo:

  The fraction of the linearly predicted gain a damped step must
  actually deliver to be accepted.

- max_iterations:

  The most Newton steps to take.

- max_halvings:

  The most times one step may be halved before a candidate is abandoned.

- threads:

  The number of worker threads.

- x, object:

  A `model_comparison`.

- ...:

  Ignored.

## Value

`compare_models()` returns a `model_comparison`, whose fields are
`table` (the ranking, best first, with columns `L_P`, `L`, `p`, `k` and
`status`), `fits` (the `fit` of each candidate, in the same order,
`NULL` where one failed), `failures` (the diagnosis of each failure,
named, empty where there were none), `Z`, `overdispersion`, `population`
(the population and test-mortality numbers the header reports),
`arguments` (what the comparison was run with, for
[`diagnose()`](https://r-logmu.logmu.org/reference/diagnose.md)) and
`labels` (what was written for the include, the weight and the test
mortality).

`table` holds the values as computed.
[`print()`](https://rdrr.io/r/base/print.html) takes the best
candidate's \\L\_{\mathrm{P}}\\ off both likelihood columns, because
differences are what matter and the absolute numbers are large and
uninformative – so the winner prints as 0. Both columns shift by the
same constant, which keeps \\L - L\_{\mathrm{P}} = p\\ true row by row.

`is_model_comparison()` returns a scalar `logical`.

`best_fit()` returns the `fit` of the highest-ranked candidate.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`summary.model_comparison`: the shifted `table` as printed, the header's
numbers, and in `fits` the
[summary()](https://r-logmu.logmu.org/reference/fit.md) of every
converged candidate, best first. Printing it shows the ranking and then
each candidate's coefficient table and correlation matrix.

## What makes the candidates comparable

Every candidate is fitted over the same data, `include`, weight,
overdispersion and \\Z\\, because all of those are arguments of this
call and none of them is part of a `model`. \\Z\\ in particular is
resolved once and shared: it is measured on `test_mortality` over this
data one time, not once per candidate.

Candidates may differ in their reference mortality. The constant
\\\mathrm{A}w\log\mu^{\mathrm{ref}}\\ is retained in \\L\\, so a
comparison across different base tables is meaningful.

## Reading the shortfall

Where the weight is absent or is an indicator, \\w^2 = w\\ makes \\p\\
exactly `k`, \\L\_{\mathrm{P}}\\ exactly \\-\\AIC\\/2\\, and a shortfall
of 1 one parameter's worth – which is the usual threshold for a
difference that matters. For an *ad hoc* weight the penalty is no longer
a count and the shortfall carries no such calibration; it still ranks
the candidates, but how large a shortfall is worth acting on is a
question this function does not answer.

Overdispersion does real work in the ranking. The penalty is free of
\\\Omega\\ and the fit term is not, so a larger \\\Omega\\ favours
simpler models, which is the quasi-likelihood behaviour and is why
`overdispersion` is yours to supply.

## Examples

``` r
data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
    male      = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    E2R_start = datey::datey(rep(2015, 5)),
    E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
    E2R_died  = c(FALSE, TRUE, TRUE, FALSE, FALSE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

# One death of each sex, so the split is identifiable.
reference <- gompertz_mortality()
comparison <- compare_models(
  data,
  models = list(
    level = model(reference, covariates(level = 1)),
    sex   = model(reference, covariates(male = .i$male, female = !.i$male))
  ),
  settings = settings(overdispersion = 1)
)
comparison
#> <model_comparison: 2 candidates, 2 converged>
#>   include         all records
#>   weight          lives
#>   data            2 deaths, 22 years of exposure
#>   test mortality  logmu default
#>                   Aw 2   Ew 0.36215   Ew^2 0.36215   Z 1
#>                   log A/E 1.7088 ± 1.6617 (1 sd) or ± 3.2569 (95%), at overdispersion 1
#>   overdispersion  1
#>   max L_P         -8.862328
#> 
#>          L_P     L     p k
#> level  0.000 1.000 1.000 1
#> sex   -0.985 1.015 2.000 2
coef(best_fit(comparison))
#>    level 
#> 1.708867 
```
