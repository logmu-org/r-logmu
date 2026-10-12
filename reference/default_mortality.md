# The default mortality

The general-purpose mortality **logmu** falls back on where one is
needed and none was given, and a reasonable thing to reach for whenever
a broadly sensible pensioner basis will do:

\$\$\log\mu\_{xt} = -3.8 + 0.1\\(x - 75) - 0.01\\(t - 2020)\$\$

It sits about midway between S4PMA and S4PFA at age 75, with a shallower
slope to allow for higher mortality at younger ages and for plateauing
at higher ones.

## Usage

``` r
default_mortality()
```

## Value

A `mortality` named `logmu default`.

## Details

Its main use is as the test mortality that sets the \\Z\\ scale of a
log-likelihood, which is why
[`fit()`](https://r-logmu.logmu.org/reference/fit.md) and
[`compare_models()`](https://r-logmu.logmu.org/reference/compare_models.md)
default `test_mortality` to it. Only `slope_x` materially matters there,
so a single fixed choice is enough and a shared one is better than a
good one: \\Z\\ is the yardstick that makes a difference of one in the
penalised log-likelihood mean one parameter's worth, so two analyses
that used different test mortalities are not on one scale.

It takes no arguments deliberately. Vary it and it is no longer the
default mortality but a Gompertz law of your own, which is what
[`gompertz_mortality()`](https://r-logmu.logmu.org/reference/gompertz_mortality.md)
is for.

It is named `logmu default`, which is how it appears wherever a
mortality is printed, as in the header of a
[`fit()`](https://r-logmu.logmu.org/reference/fit.md) or a
[`compare_models()`](https://r-logmu.logmu.org/reference/compare_models.md).

## Examples

``` r
default_mortality()
#> <mortality_expr: logmu default> ((-3.8 - (0.1 * (.i$birth - 1945-01-01.0))) + (0.09 * (.t - 2020-01-01.0)))
name(default_mortality())
#> [1] "logmu default"

# It is exactly the no-argument Gompertz, given a name.
identical(default_mortality(), gompertz_mortality(name = "logmu default"))
#> [1] TRUE
```
