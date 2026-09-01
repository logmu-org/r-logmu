# Create an `mortality` from a 2D age-period matrix of annual mortality rates

Creates an age-period mortality table from a 2D age-period matrix of
annual mortality rates that can be \\q\_{xt}\\, \\\mu\_{xt}\\ or
\\\log\mu\_{xt}\\.

The first age and period are determined by the `x0` and `t0` parameters
respectively using \\\mu\\-timing. (See below for what this means for
\\q\\ rates.)

The end age of the table is determined as the start age `x0` plus the
number of age rows as years (plus another half year for \\q\\ rates –
see below).

The annual mortality rates must be provided as a 2D age-period matrix in
one (and only one) of `q`, `mu` or `log_mu` as appropriate.

It is an implicit assumption that these rates are 'smooth' at the annual
scale, i.e. the second differences of \\\mu\\ and \\\log\mu\\ are
'small', e.g. \\\Delta^2\mu\<10\\\\\mu\\ and
\\\Delta^2\log\mu\<10\\\\\\. If you want to allow for realistic, i.e.
non-smooth, historical annual and sub-annual noise then use a
`variation`.

Important notes for \\q\\ rates:

1.  All timing is \\\mu\\-timing. This means:

    - The definition of \\q\_{xt}\\ is *centred* on \\(x,t)\\, i.e.

      \$\$q\_{xt} = 1 -
      \exp\left(-\int\_{t-\frac{1}{2}}^{t+\frac{1}{2}}\mu\_{x+\varepsilon,\\
      t+\varepsilon}\\\mathrm{d}\varepsilon\right)\$\$

    This differs from the normal convention whereby \\q\_{xt}\\ relates
    to the year from \\(x,t)\\ to \\(x+1,t+1)\\.

    - The `x0` and `t0` parameters relate to the *middle* of the year
      covered by the youngest and earliest \\q\\ rate.

    - The end age of the resulting mortality table is the start age `x0`
      plus the number of age rows as years *plus an extra half year*.

2.  The calculation of \\\mu\_{xt}\\ from \\q\_{xt}\\ includes an
    allowance for estimated convexity determined by examining the two
    neighbouring \\q\\ rates (by cohort for the interior of the
    `annual_rates` and by age at its edges). This may produce artefacts
    if the rates are not smooth at an annual scale.

3.  A common convention when specifying \\q\\-based mortality tables is
    to include \\q\_\omega=1\\, where \\\omega\\ is the end age of the
    mortality table. *Do not include a \\q=1\\ age row in the
    `annual_rates` argument.* (If you are creating the mortality from a
    base table and a projection then the \\q\_\omega=1\\ is likely in
    the base table.)

## Usage

``` r
mortality_table(x0, t0, q = NULL, mu = NULL, log_mu = NULL, name = NULL)
```

## Arguments

- x0, t0:

  The youngest age and earliest time respectively. For \\q\\ rates,
  these are the *middle* of the year covered by the youngest and
  earliest \\q\\ rate, i.e. \\\mu\\-timing.

- q:

  An age-period matrix of \\q\_{xt}\\. It is required that \\0 \< q \<
  1\\.

- mu:

  An age-period matrix of \\\mu\_{xt}\\. It is required that \\0 \< \mu
  \< +\infty\\.

- log_mu:

  An age-period matrix of \\\log\mu\_{xt}\\. It is required that
  \\-\infty \< \log\mu \< +\infty\\.

- name:

  An optional name for this mortality.

## Value

A `mortality`.
