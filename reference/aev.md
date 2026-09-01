# Actual, expected and variance over experience data

Calculates \\A\\, \\E\\ and \\V\\ over an `exp_data` under a given
mortality, optionally restricted to a sub-population and optionally
broken down into groups.

`mortality` and `weight` are **pronoun expressions**, so they may be
written out in place (`.i$pension`), or reference a `mortality`,
`variable` or `indicator` object by name, or combine the two.

## Usage

``` r
aev(
  exp_data,
  mortality,
  include = NULL,
  weight = NULL,
  val_similarity = NULL,
  val_distance = NULL,
  breakdown = NULL,
  settings = NULL,
  overdispersion = NULL,
  time_scale = NULL,
  threads = cpp_veil_default_threads()
)
```

## Arguments

- exp_data:

  The experience data.

- mortality:

  A pronoun expression for \\\log\mu\\, or a `mortality`.

- include:

  An `include` naming the population. `NULL` means everybody.

- weight:

  A pronoun expression for the weight \\w\\, or a `variable`. `NULL`
  means a weight of 1, i.e. a count of lives.

- val_similarity, val_distance:

  Two spellings of the second weighting factor \\s\\, related by \\d =
  -\log s\\. Give one or the other. A similarity belongs in \\\[0, 1\]\\
  and a distance must not be negative; see *Similarity and distance*,
  *Why the names begin with `val_`* and *The bound on a similarity*.

- breakdown:

  An `include` or
  [`includes()`](https://r-logmu.logmu.org/reference/includes.md)
  dividing the population into groups. `NULL` gives a single ungrouped
  result.

- settings:

  A [`settings()`](https://r-logmu.logmu.org/reference/settings.md)
  object supplying `overdispersion` and `time_scale`.

- overdispersion, time_scale:

  Given directly, these override `settings`.

- threads:

  Worker threads to use. `0` asks for as many as the machine reports.
  Cannot change any answer.

## Value

An `aev` with one record per breakdown element, or one record if there
is no breakdown.

## Breakdowns

With no `breakdown` the result is a length-1 `aev`. With one, the result
is a single `aev` with one record per element – not a list – carrying
the breakdown's [`names()`](https://rdrr.io/r/base/names.html) and
[`group_names()`](https://r-logmu.logmu.org/reference/includes.md) so a
chart can tell which records were ages and which were amounts.

Every element is intersected with `include`, so `include` says who is in
the population and `breakdown` says how to divide them. Elements need
not be disjoint and need not cover everybody: examining A/E on
intersecting subsets is legitimate, and a record outside every element
is simply absent.

## Overdispersion

`overdispersion` is required, either directly or through `settings`. It
scales `V`, and so scales every confidence interval and residual read
from the result. There is no default: see
[`settings()`](https://r-logmu.logmu.org/reference/settings.md).

## Similarity and distance

A second weighting factor \\s\\ says how much a record should count at
all, as opposed to how much of it counts. It differs from `weight` in
where it appears: the weight is squared in \\V\\ and the similarity is
not, so \\A\\ and \\E\\ take \\sw\\ while \\V\\ takes \\sw^2\\. The
effect is that halving a record's similarity doubles the variance it
implies, while leaving the number of parameters a model may support
unchanged.

It may be written either way round, and they are one quantity:

- `val_similarity` is \\s\\ itself, intended to lie in \\\[0, 1\]\\,
  where 1 counts a record fully and 0 not at all.

- `val_distance` is \\d = -\log s\\, so 0 counts a record fully,
  `log(2)` counts it at half, and larger counts it less.

Give one or the other, never both. `val_distance` is usually the easier
of the two: a decay kernel is `val_distance = (2025 - .t) / 10` and a
Gaussian is `val_distance = (x / h)^2`, where the similarity form would
need the exponential written out each time. Distances also add, so
several reasons to discount a record combine by addition.

A similarity of zero counts a record not at all, but
[`include()`](https://r-logmu.logmu.org/reference/include.md) is the way
to leave records out: it clips exposure rather than multiplying through
it.

## Why the names begin with `val_`

A similarity or a distance is ordinarily a relation between two things.
These arguments are functions of one record and one time, so the far end
of the relation is left implicit: it is the valuation the analysis is
aimed at, taken as a whole and at its as-at date. Leaving it implicit is
what allows them to be written as ordinary variables of the experience
record.

You write a `val_` argument yourself, and nothing derives it from
valuation data. The prefix keeps `similarity` and `distance` free for a
later, more general form that would take two records and two times.

## The bound on a similarity

A similarity is a proportion and must lie in \\\[0, 1\]\\; a distance,
being \\-\log s\\, must not be negative.

Neither is enforced value by value. Both are expressions of the record
and of time, so their values are not known until the walk, and testing
each one would put a comparison in the innermost loop to police
something already stated here.

What **logmu** does instead is analytic, and costs nothing. It works out
the range each expression can take – from the numbers written in it and
the range of each column in the data – and refuses the calculation
before reading a record if that range lies **wholly** outside the bound.
So `val_similarity = 2` is refused, and so is a similarity given as a
column whose values run from 5 to 200.

It refuses only what is certainly wrong. A range that merely permits a
violation is accepted, because the arithmetic that derives it is
conservative: `val_distance = (2025 - .t) / 10` would otherwise be
rejected, since the range of `.t` alone allows a negative distance that
no exposure in the data reaches. A check that fires on correct code
would be worse than one that occasionally stays quiet.

## Examples

``` r
data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955)),
    pension   = c(5000, 12000, 30000),
    E2R_start = datey::datey(c(2015, 2015, 2015)),
    E2R_end   = datey::datey(c(2020, 2020, 2018)),
    E2R_died  = c(FALSE, FALSE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

basis <- settings(overdispersion = 2)

aev(data, mortality = mortality_const(log_mu = -4), settings = basis)
#> <aev[1]>
#>   A         E         V      A/E 95% conf dev resid
#> 1 1 0.2381033 0.4762066 4.199858 5.680419 0.8204596

# Weighted by pension, broken down by amount.
aev(data,
    mortality = mortality_const(log_mu = -4),
    weight    = .i$pension,
    breakdown = bands(.i$pension, thresholds = c(10000, 20000)),
    settings  = basis)
#> <aev[3]>
#>     group        name     A        E        V      A/E  95% conf  dev resid
#> 1 pension     < 10000     0  457.891  4578910  0.00000  9.159401 -0.3026189
#> 2 pension 10000-20000     0 1098.938 26374520  0.00000  9.159401 -0.3026189
#> 3 pension    >= 20000 30000 1648.407 98904450 18.19938 11.824735  1.3986903
```
