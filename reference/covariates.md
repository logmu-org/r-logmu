# Covariates: the fitted terms of a mortality model

`covariates()` collects the \\X\\ terms of a proportional hazards model,
\\\log\mu = \log\mu^{\mathrm{ref}} + \sum_j \beta_j X^{(j)}\\. Each
argument is a pronoun expression, and the coefficients belong to the
fit.

`disjoint()` is the same, and additionally asserts that the terms are
mutually exclusive indicators – no individual is in more than one. That
is worth stating: where it holds, every off-diagonal entry of the
information matrix is zero, and the fit does far less work.

Multiplying by a variable distributes over the list, which is how the
common shape \\X^{(k)} = I^{(k)} \varphi\\ is written: a set of
indicators times a shared age shape. Multiplying two covariate lists
gives their Cartesian product.

## Usage

``` r
covariates(...)

disjoint(...)

is_covariates(x)

is_disjoint(x)
```

## Arguments

- ...:

  Pronoun expressions, optionally named. A `covariates` object is
  spliced in.

- x:

  An object.

## Value

`covariates()` and `disjoint()` return a `covariates`.

`is_covariates()` returns a scalar `logical`.

`is_disjoint()` returns a scalar `logical` saying whether the terms were
asserted mutually exclusive.

## Examples

``` r
age_shape <- variable(.x)
covariates(.i$sex == "M", .i$sex == "F") * age_shape
#> <covariates: 2 terms>
#>   [1] ((.i$sex == "M") * (.t - .i$birth))
#>   [2] ((.i$sex == "F") * (.t - .i$birth))
```
