# Fit a mortality model by maximum weighted log-likelihood

Estimates the coefficients of a
[`model()`](https://r-logmu.logmu.org/reference/model.md) against
experience data, by Newton-Raphson on the weighted log-likelihood
\$\$L^\* = \Omega^{-1}(\mathrm{A}w\log\mu - \mathrm{E}w).\$\$

A fit returns a fit or it fails. A model that does not converge, or
whose information matrix cannot be factored, raises an error carrying
the diagnosis rather than returning a number nobody should use.

## Usage

``` r
fit(
  exp_data,
  model,
  include = NULL,
  weight = NULL,
  val_relevance = NULL,
  val_distance = NULL,
  test_mortality = default_mortality(),
  Z = NULL,
  settings = NULL,
  overdispersion = NULL,
  time_scale = NULL,
  start = NULL,
  L_tolerance = default_L_tolerance,
  armijo = default_armijo,
  max_iterations = default_max_iterations,
  max_halvings = default_max_halvings,
  threads = cpp_veil_default_threads()
)

is_fit(x)

# S3 method for class 'fit'
coef(object, ...)

# S3 method for class 'fit'
vcov(object, ...)

# S3 method for class 'fit'
residuals(object, ...)

# S3 method for class 'fit'
summary(object, ...)

# S3 method for class 'summary.fit'
print(x, ...)

# S3 method for class 'fit'
print(x, ...)
```

## Arguments

- exp_data:

  The experience data.

- model:

  The [`model()`](https://r-logmu.logmu.org/reference/model.md) to fit.

- include:

  An `include` restricting the population, or `NULL` for all of it.

- weight:

  A pronoun expression for the weight \\w\\, or `NULL` for lives.
  Weights must not be negative.

- val_relevance, val_distance:

  A pronoun expression for the second weighting factor, in either
  spelling. At most one may be given.

- test_mortality:

  The mortality on which the \\Z\\ scale is measured. See Details.

- Z:

  The \\Z\\ scale as a single positive number, given instead of a
  `test_mortality` rather than as well as one.

- settings:

  A [`settings()`](https://r-logmu.logmu.org/reference/settings.md)
  object carrying `overdispersion` and `time_scale`.

- overdispersion, time_scale:

  Given directly, these override the `settings`.

- start:

  Starting coefficients, one per model term. Defaults to zeros.

- L_tolerance:

  The convergence tolerance, as a gain in \\L\\. One unit of \\L\\ is
  about what one parameter costs.

- armijo:

  The fraction of the linearly predicted gain a damped step must
  actually deliver to be accepted. Must lie strictly between 0 and 0.5.

- max_iterations:

  The most Newton steps to take.

- max_halvings:

  The most times one step may be halved before the fit is abandoned.

- threads:

  The number of worker threads.

- x, object:

  A `fit`.

- ...:

  Ignored.

## Value

`fit()` returns a `fit`, whose fields are `beta`, `mortality` (the
fitted model as a `mortality`), `variance`, `L`, `p`, `L_P`, `Z`,
`overdispersion`, `population` (the population and test-mortality
numbers [`print()`](https://rdrr.io/r/base/print.html) reports),
`arguments` (the include, weight and relevance the fit was run with,
kept so that
[`diagnose()`](https://r-logmu.logmu.org/reference/diagnose.md) can
reuse them), `labels` (what was written for the include, the weight and
the test mortality), `iterations`, `evaluations`, `predicted_gain` and
`model`.

`is_fit()` returns a scalar `logical`.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the fitted
coefficients as a named numeric vector, and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) their variance-covariance
matrix.

[`residuals()`](https://rdrr.io/r/stats/residuals.html) stops with a
message: a `fit` keeps no data, so residuals come from an
[`aev()`](https://r-logmu.logmu.org/reference/aev.md) on the fitted
mortality.

[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`summary.fit`: the `coefficients` table (`estimate`, `std_error` for one
standard deviation and `conf_95pc` for the 95% half-width), the
`correlation` matrix of the estimates, and the header's numbers. A fit
prints as its summary.

## The Z scale

Results are reported on the \\L = Z^{-1}L^\*\\ scale, where \\Z =
\mathrm{E}w^2/\mathrm{E}w\\ on a test mortality. \\Z\\ is a change of
units, chosen so that one parameter costs about one unit of \\L\\, and
it must be the same for every model being compared – otherwise the
penalised log-likelihoods are not on one scale and the comparison means
nothing.

It is set in one of two ways, never from the model's own reference
mortality, and they are alternatives rather than a pair:

- `test_mortality` names the mortality \\Z\\ is measured on, over this
  same data, weight and include. It defaults to
  [`default_mortality()`](https://r-logmu.logmu.org/reference/default_mortality.md),
  and a fixed shared choice is worth more here than a good one, since
  \\Z\\'s job is to be one yardstick.

- `Z` gives \\Z\\ directly as a number, skipping the measurement. This
  is what model selection needs, the caller having obtained it once and
  passed it to every candidate, and
  [`compare_models()`](https://r-logmu.logmu.org/reference/compare_models.md)
  does exactly that.

Where the weight is absent, is the literal 1, or takes only the values 0
and 1 – a logical expression, or a column holding nothing but zeros and
ones – \\w^2 = w\\ makes \\Z\\ exactly one. The engine recognises such a
weight and accumulates \\\mathrm{E}w^2\\ as \\\mathrm{E}w\\ itself, so
the measured \\Z\\ is one to the last bit rather than to a tolerance.
The measurement is made regardless, since the same pass over the data
supplies the population figures that
[`print()`](https://rdrr.io/r/base/print.html) reports.

## Convergence

The iteration stops when the gain in \\L\\ still available falls below
`L_tolerance`, which is a forward-looking test: it measures the distance
left to the optimum rather than the size of the last step, so it is not
fooled by a damped step or by a model whose maximum does not exist.
Stopping within \\\epsilon\\ of the maximum leaves the coefficients
within \\\sqrt{2\epsilon/I}\\ of theirs.

## Examples

``` r
data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955, 1940, 1948)),
    male      = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    E2R_start = datey::datey(rep(2015, 5)),
    E2R_end   = datey::datey(c(2020, 2020, 2018, 2020, 2019)),
    E2R_died  = c(FALSE, FALSE, TRUE, FALSE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

m <- model(ref_mortality = gompertz_mortality(), covariates = covariates(level = 1))
fit(data, m, settings = settings(overdispersion = 1))
#> <fit: 1 term, 5 iterations, 6 walks of the data>
#>   include         all records
#>   weight          lives
#>   data            2 deaths, 22 years of exposure
#>   test mortality  logmu default
#>                   Aw 2   Ew 0.36215   Ew^2 0.36215   Z 1
#>                   log A/E 1.7088 ± 1.6617 (1 sd) or ± 3.2569 (95%), at overdispersion 1
#>   overdispersion  1
#> 
#>       estimate std_error conf_95pc
#> level    1.709    0.7071     1.386
#> 
#> L -7.752328   p 1   L_P -8.752328   k 1
```
