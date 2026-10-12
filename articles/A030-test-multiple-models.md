# Test multiple models

``` r

library(logmu)
```

[`fit()`](https://r-logmu.logmu.org/reference/fit.md) calibrates one
model. It cannot tell you whether that model was worth having, because a
model with more parameters always fits the data it was calibrated on
better than a model with fewer.
[`compare_models()`](https://r-logmu.logmu.org/reference/compare_models.md)
fits several candidates and ranks them on penalised log-likelihood,
which is Akaike’s information criterion with the clutter removed. The
reasoning behind it is in [Suddenly
AIC](https://r-logmu.logmu.org/articles/B050-suddenly-aic.md); this
article is about running it.

## A portfolio to work on

There is no data in the package, so here is some, drawn by
[`simulate()`](https://r-logmu.logmu.org/reference/simulate.md) from a
mortality we will pretend not to know. The truth the candidates have to
find is a level, an age slope steeper than the reference’s, and a
loading on males; how to build such a portfolio is in [Fit a
model](https://r-logmu.logmu.org/articles/A020-fit-a-model.md).

``` r

set.seed(20260910)

lives <- 4000
population <- data.frame(
  birth   = datey::datey(2015 - sample(55:89, lives, replace = TRUE)),
  male    = runif(lives) < 0.6,
  pension = round(exp(rnorm(lives, log(12000), 0.7)))
)

truth <- mortality(
  if (.i$male) {
    gompertz_mortality(log_mu_0 = -3.5, slope_x = 0.12, slope_t = 0)
  } else {
    gompertz_mortality(log_mu_0 = -3.9, slope_x = 0.12, slope_t = 0)
  }
)

portfolio <- simulate(truth, lives = population,
                      exp_start = datey::datey(2015), exp_end = datey::datey(2020))
sum(portfolio$E2R_died)
#> [1] 743
```

## Four candidates

Every candidate is a
[`model()`](https://r-logmu.logmu.org/reference/model.md): a reference
mortality, plus the covariates whose coefficients are to be estimated. A
`model` holds nothing else. The data, the population, the weight and the
overdispersion are arguments of the comparison, not properties of a
candidate, and that is what makes candidates comparable rather than a
rule you have to remember.

The reference here is a Gompertz curve with its time trend switched off,
so anything the candidates find is theirs.

``` r

reference <- gompertz_mortality(slope_t = 0)

candidates <- list(
  level    = model(reference, covariates(level = 1)),
  sex      = model(reference, covariates(male = .i$male, female = !.i$male)),
  age      = model(reference, covariates(level = 1, age = .x - 75)),
  sex_age  = model(reference, covariates(male = .i$male, female = !.i$male,
                                         age = .x - 75))
)
```

The `sex` and `sex_age` candidates split the population into two
indicators that cannot both be one. Saying so with
[`disjoint()`](https://r-logmu.logmu.org/reference/covariates.md)
instead of
[`covariates()`](https://r-logmu.logmu.org/reference/covariates.md) lets
the engine leave the off-diagonal integrals out altogether, which is
worth having once there are more than a handful of groups.

## The comparison

``` r

comparison <- compare_models(portfolio, candidates,
                             overdispersion = 1)
comparison
#> <model_comparison: 4 candidates, 4 converged>
#>   include         all records
#>   weight          lives
#>   data            743 deaths, 18183 years of exposure
#>   test mortality  logmu default
#>                   Aw 743   Ew 573   Ew^2 573   Z 1
#>                   log A/E 0.25981 ± 0.041776 (1 sd) or ± 0.081879 (95%), at overdispersion 1
#>   overdispersion  1
#>   max L_P         -2739.508
#> 
#>             L_P       L     p k
#> sex_age   0.000   3.000 3.000 3
#> sex      -2.862  -0.862 2.000 2
#> age     -19.217 -17.217 2.000 2
#> level   -21.702 -20.702 1.000 1
```

Larger L\_\text{P} is better, so the table is ordered best first. The
header describes the data and the scale everything below it is quoted
in.

The two likelihood columns have the best candidate’s L\_\text{P} taken
off them, so the winner reads 0 and everyone else is negative. Both are
shifted by the same constant, which is why L - L\_\text{P} = p still
holds row by row: the arithmetic can be checked on the face of the
table. `$table` keeps the unshifted numbers.

Note the `p` column against the `k` beside it: the penalty is exactly
the parameter count. That is not a coincidence and it is not hard-wired
– the penalty is p = Z^{-1}\operatorname{tr}(\mathbf{J}\mathbf{I}^{-1}),
and with no weight w^2 = w makes it collapse to \dim(\beta). So here,
and only here, L\_\text{P} is exactly -\text{AIC}/2 and a shortfall of 1
is one parameter’s worth. Anything under about 1 is not worth acting on.

The winner is a `fit` like any other:

``` r

best <- best_fit(comparison)
coef(best)
#>        male      female         age 
#>  0.35807146 -0.13867752  0.01309349
sqrt(diag(vcov(best)))
#>        male      female         age 
#> 0.057812446 0.078338077 0.004795356
```

The two sex coefficients differ by 0.50, against the 0.4 that went into
the data, about one standard error over. The age term comes out at 0.013
where 0.02 went in, about one and a half standard errors short – the
levels and the slope are estimated jointly, and are correlated, so they
should be read together rather than one at a time.

The fitted model is a `mortality`, so `best$mortality` can be handed
straight to [`aev()`](https://r-logmu.logmu.org/reference/aev.md), or
used as the `ref_mortality` of a further model.

## Weights change what the penalty means

Weight the log-likelihood by pension and the ranking still works, but
two things change.

First, Z is no longer 1, so it has to be measured. It is measured on
`test_mortality`, which defaults to
[`default_mortality()`](https://r-logmu.logmu.org/reference/default_mortality.md),
and it has to be the same for every candidate or the penalised
log-likelihoods are on different scales and the comparison means
nothing.
[`compare_models()`](https://r-logmu.logmu.org/reference/compare_models.md)
measures it a single time and hands the same value to every candidate,
whatever the reference mortalities of the candidates themselves happen
to be.

Nominate a different test mortality if you have a reason to. Only its
age slope materially matters, and a shared choice is worth more than a
good one, so the default is usually the right answer. Where you already
have the number, pass `Z` instead and no measurement happens at all.

Second, p stops being a count.

``` r

weighted <- compare_models(portfolio, candidates,
                           weight = .i$pension,
                           test_mortality = reference,
                           overdispersion = 1)
weighted$table[, c("L_P", "p", "k")]
#>               L_P        p k
#> sex_age -1661.860 3.068963 3
#> sex     -1664.756 2.000682 2
#> age     -1671.446 2.056258 2
#> level   -1674.190 1.000000 1
```

Mostly it stops being a count. `level` still costs exactly 1, and that
is not rounding: with a single constant covariate
\mathbf{J}\mathbf{I}^{-1} is the scalar \mathrm{E}w^2/\mathrm{E}w, which
is Z itself, so p = Z^{-1}Z = 1 however the data is weighted. The others
sit a little above their term counts, three parameters costing 3.07.

The spread can change too. Weighting concentrates the data on the larger
pensions, so the evidence for a term can grow or shrink depending on
where the money sits; here the gap between `sex_age` and `sex` happens
to be much the same as on lives, and the ranking survives. What you no
longer get for free is the calibration – a shortfall of 1 is no longer
one parameter’s worth, and [Suddenly
AIC](https://r-logmu.logmu.org/articles/B050-suddenly-aic.html#just-weight-a-moment)
sets out why.

## Overdispersion favours simpler models

Mortality data is overdispersed, always, and `overdispersion` has no
default for that reason. It is not a nuisance in this calculation: it
changes the answer, and in the right direction.

L carries \Omega^{-1} and the penalty does not, because the \Omega^{-1}
in \mathbf{I} and \mathbf{J} cancels in \mathbf{J}\mathbf{I}^{-1}. So
raising \Omega shrinks the advantage a richer model has while leaving
what its extra parameters cost untouched.

``` r

behind <- function(omega) {
  table <- compare_models(portfolio, candidates, overdispersion = omega)$table
  max(table$L_P) - table["level", "L_P"]
}

c(omega_1 = behind(1), omega_5 = behind(5), omega_25 = behind(25))
#>   omega_1   omega_5  omega_25 
#> 21.701713  2.968005  0.000000
```

At \Omega = 1 the level-only model is 22 behind. At 5 it is under 3. At
25 its shortfall is zero, meaning it has won outright – at that much
overdispersion the evidence for a sex split is not there, and a single
level is all the data supports. This is the quasi-AIC, and it is why
\Omega deserves an argument rather than a default.

## When a candidate cannot be fitted

[`fit()`](https://r-logmu.logmu.org/reference/fit.md) raises rather than
return a number nobody should use. A comparison cannot afford to: one
hopeless candidate out of fifty must not take the other forty-nine with
it. So a candidate that fails is reported in place – its row carries the
status and `NA`, it sorts last, and the diagnosis
[`fit()`](https://r-logmu.logmu.org/reference/fit.md) would have raised
is kept.

``` r

broken <- c(candidates,
            list(collinear = model(reference,
                                   covariates(level = 1, again = 1))))

comparison <- compare_models(portfolio, broken, overdispersion = 1)
#> Warning: 1 of 5 candidates could not be fitted: `collinear` (not_identifiable).
#> See `$failures` for why.
comparison$failures[["collinear"]]
#> [1] "The model is not identifiable. Term 2 (again) is a linear combination of the terms before it, so its coefficient is not determined by the data. It equals 1 * [1] level."
```

The warning names the failures as they happen, so a candidate cannot go
missing from a ranking quietly.
