# Includes: a single time interval per individual

An `include` maps each individual to a single clopen time interval
\\\[\nu, \tau)\\ of exposure to count. To guarantee the interval is
always a single, well-defined period, includes are built only through
these constructors – never from raw `.t` expressions.

The whole family is one construct with three shorthands:

- `band()` and `bands()` band any permitted variable.

- `age()` and `ages()` are shorthand for banding `.x`, i.e. `.t - .b`.

- `duration()` and `durations()` are shorthand for banding `.t` minus an
  origin, e.g. `.t - .i$entry`. `age()` is the special case whose origin
  is birth.

- `period()` and `periods()` are shorthand for banding `.t`.

Every singular constructor returns one `include`; every plural
constructor returns an
[`includes()`](https://r-logmu.logmu.org/reference/includes.md), even
when it holds a single band. So no constructor's return type depends on
which arguments were supplied.

Banding cuts a variable at edges, so it suits a continuous one. A
discrete field is divided by
[`category()`](https://r-logmu.logmu.org/reference/category.md) and
[`categories()`](https://r-logmu.logmu.org/reference/category.md)
instead, which follow the same singular and plural rule.

## Usage

``` r
band(variable, from = NULL, to = NULL)

bands(variable, from = NULL, to = NULL, by = NULL, thresholds = NULL)

age(from = NULL, to = NULL)

ages(from = NULL, to = NULL, by = NULL, thresholds = NULL)

duration(since, from = NULL, to = NULL)

durations(since, from = NULL, to = NULL, by = NULL, thresholds = NULL)

period(from = NULL, to = NULL)

periods(from = NULL, to = NULL, by = NULL, thresholds = NULL)

include(expr)

is_include(x)

period_included(x, .i)

# Default S3 method
period_included(x, .i)

# S3 method for class 'include'
period_included(x, .i)

# S3 method for class 'indicator'
period_included(x, .i)
```

## Arguments

- variable:

  A pronoun expression to band. See *What may be banded*.

- from, to:

  Outer band edges, or `NULL` for unbounded. Dates (`datey` or year
  numbers) when banding `.t`; durations (`durationy` or year counts)
  when banding `.t` minus an origin; ordinary numbers otherwise.

- by:

  The width of each band. Requires both `from` and `to`, and must divide
  `to - from` exactly.

- thresholds:

  Interior band edges, given explicitly instead of `by`.

- since:

  The origin a duration is measured from: a time-invariant `datey`
  expression, e.g. `.i$entry` or `min(.i$entry, .i$retirement)`.

- expr:

  A time-invariant pronoun expression.

- x:

  An `include`.

- .i:

  A named list of scalar facts.

## Value

`band()`, `age()`, `duration()` and `period()` return an `include`;
`include()` returns an `indicator`, which is also an `include`.

`bands()`, `ages()`, `durations()` and `periods()` return an `includes`.

`is_include()` returns a scalar `logical`.

`period_included()` returns a `datey_interval`.

## What may be banded

The banded variable must be **increasing in `.t` at unit slope**, which
permits exactly three shapes:

- a time-invariant variable, e.g. `.i$pension` – this *gates* exposure,
  since the individual is in one band throughout;

- `.t` itself – this *clips* exposure to calendar bounds;

- `.t` minus a time-invariant `datey` expression, e.g.
  `.t - .i$retirement`, `.t - min(.i$entry, .i$retirement)` (and `.x`,
  which is `.t - .b`) – this clips exposure to bounds measured from that
  origin.

A decreasing expression such as `2020 - .t` is refused: it would flip
the band to `(a, b]`, so adjacent bands would either double-count or
drop an instant.

## Band edges

The edges are `c(from, thresholds, to)`, and a `NULL` bound means
unbounded. `by` generates the interior thresholds from `from` and `to`,
and must divide `to - from` exactly. `by` and `thresholds` cannot both
be given.

So `ages(65, 95, by = 5)` gives six bands and excludes the outside,
because both bounds were supplied, while
`bands(.i$pension, thresholds = c(5000, 10000, 20000))` gives four bands
open at both ends, because neither was.

## Missing values

A `NaN` banded variable compares `FALSE` against every threshold,
following IEEE 754, so such a record falls in no band and is simply
absent. Note that base R differs here: `NaN < 5` is `NA` in R, which
needs a third logical state **logmu** does not have.

## Combining

Includes (and indicators) combine with `&`, which intersects: the result
selects the time within *both* operands. Internally an include is a
conjunction of terms – interval bounds plus indicator gates – resolved
together by `period_included()`.

`period_included()` resolves the interval for a single individual's
facts `.i`. It is the plain-R reference path for testing and
understanding, not the performance path. An offset that is `NA` yields
the empty interval.

## Examples

``` r
period_included(period(2010, 2020), .i = list())
#> [1] [2010-01-01.0, 2020-01-01.0)
period_included(age(65, 95), .i = list(birth = datey::datey(1950)))
#> [1] [2015-01-01.0, 2045-01-01.0)
period_included(age(65, 95) & period(2010, 2040),
                .i = list(birth = datey::datey(1950)))
#> [1] [2015-01-01.0, 2040-01-01.0)

# Six age bands, outside excluded.
ages(65, 95, by = 5)
#> <includes[6]>
#>   age: 65-70
#>   age: 70-75
#>   age: 75-80
#>   age: 80-85
#>   age: 85-90
#>   age: 90-95

# Four amount bands, open at both ends.
bands(.i$pension, thresholds = c(5000, 10000, 20000))
#> <includes[4]>
#>   pension: < 5000
#>   pension: 5000-10000
#>   pension: 10000-20000
#>   pension: >= 20000

# Duration since an origin other than birth.
durations(.i$retirement, 0, 10, by = 5)
#> <includes[2]>
#>   duration since retirement: 0-5
#>   duration since retirement: 5-10

# The origin may be computed, not only a bare field.
band(.t - min(.i$entry, .i$retirement), 0, 5)
#> <include: .t - min(.i$entry, .i$retirement) 0-5>
#>   since min(.i$entry, .i$retirement) [0 yr, 5 yr)
```
