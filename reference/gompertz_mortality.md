# Create a Gompertz `mortality` with an age and a period slope

Creates a `mortality` whose \\\log\mu\\ is linear in age and in time:

\$\$\log\mu\_{xt} = \texttt{log\\mu\\0} + \texttt{slope\\x}\\(x -
\texttt{x\\0}) + \texttt{slope\\t}\\(t - \texttt{t\\0})\$\$

Both slopes are additive, so a period improvement is a *negative*
`slope_t`. That is the usual case: mortality rises with age and falls
with time.

This is the natural choice of *test* mortality for the \\Z\\ scaling of
a log-likelihood, where only `slope_x` materially matters, and it also
serves as a reference mortality or as a simple model in its own right.

Every parameter has a default, so `gompertz_mortality()` on its own
gives a broadly reasonable pensioner basis:

\$\$\log\mu\_{xt} = -3.8 + 0.1\\(x - 75) - 0.01\\(t - 2020)\$\$

## Usage

``` r
gompertz_mortality(
  log_mu_0 = -3.8,
  slope_x = 0.1,
  x_0 = 75,
  slope_t = -0.01,
  t_0 = 2020,
  name = NULL
)
```

## Arguments

- log_mu_0:

  The value of \\\log\mu\\ at age `x_0` and time `t_0`.

- slope_x:

  The slope of \\\log\mu\\ in age, per year. Positive for mortality that
  rises with age.

- x_0:

  The age at which \\\log\mu\\ takes the value `log_mu_0`, as a
  `durationy` or a number of years.

- slope_t:

  The slope of \\\log\mu\\ in time, per year. Typically negative, i.e.
  mortality improving.

- t_0:

  The time at which \\\log\mu\\ takes the value `log_mu_0`, as a `datey`
  or a year.

- name:

  An optional name for this mortality, read back by
  [`name()`](https://r-logmu.logmu.org/reference/name.md) and printed
  wherever the mortality appears.

## Value

A `mortality`.

## Details

Internally the law is held in its cohort form, which is the shape that
integrates cheaply. Writing \\b\\ for the date of birth and \\b_0 =
\texttt{t\\0} - \texttt{x\\0}\\ for the cohort origin the two offsets
imply,

\$\$\log\mu\_{bt} = \texttt{log\\mu\\0} - \texttt{slope\\x}\\(b - b_0) +
(\texttt{slope\\x} + \texttt{slope\\t})\\(t - \texttt{t\\0})\$\$

which is identical to the age-period form above because \\x = t - b\\.
Birth is fixed for an individual, so the first two terms are a
per-individual constant that lifts out of the integral, leaving a single
term linear in \\t\\ whose slope every individual shares.

## Examples

``` r
# The default basis.
log_mu(gompertz_mortality(), list(birth = datey::datey(1945)), datey::datey(2020))
#> [1] -3.8

m <- gompertz_mortality(
  log_mu_0 = -4.6, slope_x = 0.10, x_0 = 70,
  slope_t = -0.02, t_0 = 2020
)
# log mu at age 70 in 2020, for someone born in 1950:
log_mu(m, list(birth = datey::datey(1950)), datey::datey(2020))
#> [1] -4.6
```
