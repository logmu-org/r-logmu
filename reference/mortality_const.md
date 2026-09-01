# Create a constant `mortality`

Creates a constant mortality from \\q\\, \\\mu\\ or \\\log\mu\\.

The annual mortality rate must be provided as one (and only one) of `q`,
`mu` or `log_mu` as appropriate.

## Usage

``` r
mortality_const(q = NULL, mu = NULL, log_mu = NULL, name = NULL)
```

## Arguments

- q:

  The probability of dying over one years. It is required that \\0 \< q
  \< 1\\.

- mu:

  The instantaneous annual rate ('force') of mortality. It is required
  that \\0 \< \mu \< +\infty\\.

- log_mu:

  The natural logarithm of the instantaneous annual rate ('force') of
  mortality. It is required that \\-\infty \< \log\mu \< +\infty\\.

- name:

  An optional name for this mortality.

## Value

A `mortality`.
