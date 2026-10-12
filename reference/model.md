# A proportional hazards mortality model

A `model` is the right-hand side of \$\$\log\mu =
\log\mu^{\mathrm{ref}} + \sum_j \beta_j X^{(j)},\$\$ that is a reference
mortality and the covariates whose coefficients
[`fit()`](https://r-logmu.logmu.org/reference/fit.md) estimates. It
carries nothing else: the experience data, the `include`, the weight and
the tolerance all belong to the fit, because every candidate model in a
selection exercise shares them.

A fitted model is itself a `mortality`, so it can be handed back as the
`ref_mortality` of another model.

A model may carry a `name`. It is how
[`compare_models()`](https://r-logmu.logmu.org/reference/compare_models.md)
labels the candidate when the list of models does not, how
[`diagnose()`](https://r-logmu.logmu.org/reference/diagnose.md) labels a
lone fit, and it is printed with the model and with its fit. A name
given in the list handed to
[`compare_models()`](https://r-logmu.logmu.org/reference/compare_models.md)
wins over the model's own.

## Usage

``` r
model(ref_mortality, covariates, name = NULL)

is_model(x)

# S3 method for class 'model'
print(x, ...)
```

## Arguments

- ref_mortality:

  The reference mortality: a `mortality`, or a pronoun expression for
  \\\log\mu^{\mathrm{ref}}\\. Anything held fixed goes here – a
  reference mortality is already an offset, so fixed effects need no
  feature of their own.

- covariates:

  A [`covariates()`](https://r-logmu.logmu.org/reference/covariates.md)
  object giving the \\X\\ terms. The coefficients belong to the fit;
  nothing here names or supplies one.

- name:

  An optional name for this model: a single string that is not blank,
  has no leading, trailing or consecutive spaces, and holds no control
  character such as a tab.

- x:

  An object.

- ...:

  Ignored.

## Value

`model()` returns a `model`.

`is_model()` returns a scalar `logical`.

## Examples

``` r
age_shape <- variable(.x)
model(
  ref_mortality = gompertz_mortality(),
  covariates    = covariates(male = .i$male, female = !.i$male) * age_shape
)
#> <model: 2 terms>
#>   reference: ((-3.8 - (0.1 * (.i$birth - 1945-01-01.0))) + (0.09 * (.t - 2020-01-01.0)))
#>   [1] male
#>   [2] female
```
