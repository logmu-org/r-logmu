# Categories: breaking a population down by a categorical field

`category()` and `categories()` build
[include](https://r-logmu.logmu.org/reference/include.md)s from the
values of a categorical field such as a scheme, a status or a sex. They
complete the constructor family:
[`bands()`](https://r-logmu.logmu.org/reference/include.md) cuts a
continuous variable, `categories()` divides a discrete one.

Each argument in `...` is **one group**, so a character vector joins its
values into a single group rather than splitting into several. Names
label the groups.

`.source` says where the field's values come from. With no `...` it
supplies the groups, one per value; with `...` it is what the written
categories are checked against.

Every group is an
[indicator](https://r-logmu.logmu.org/reference/indicator.md), so \\V =
\Omega E\\ holds within each.

## Usage

``` r
category(variable, ..., .source = NULL, .known = !is.null(.source))

categories(
  variable,
  ...,
  .source = NULL,
  .known = !is.null(.source),
  .exhaustive = !is.null(.source),
  .disjoint = TRUE
)
```

## Arguments

- variable:

  A pronoun expression naming the field to divide, e.g. `.i$status`. It
  must not depend on time `.t`.

- ...:

  The groups. Each argument is one group; a character vector of several
  values is one group joining them. A name labels the group.

- .source:

  Where the field's values come from – a factor, a character vector or
  column, or a dataset holding a column of the field's name. `NULL` for
  none. Named with a dot because every other name belongs to `...`.

- .known:

  Whether every category in `...` must be a value of `.source`.

- .exhaustive:

  Whether every value of `.source` must fall in some group.

- .disjoint:

  Whether no value may fall in more than one group.

## Value

`category()` returns an `indicator`, which is also an `include`.

`categories()` returns an
[`includes()`](https://r-logmu.logmu.org/reference/includes.md).

## Where the values come from

`.source` may be:

- a **factor**, giving its declared
  [`levels()`](https://rdrr.io/r/base/levels.html) in their declared
  order, including any level no record uses;

- a **character** vector or column, giving its distinct values in
  alphabetical order;

- a **dataset** (`exp_data`, `val_data` or a data frame), giving the
  column named by the field, which requires `variable` to be a plain
  field.

Prefer a factor when the breakdown is to be reused across several
datasets. A factor's levels are a property of the column rather than of
its contents, so every dataset yields the same groups in the same order
and the results stay comparable. Two character columns need not.

## Checking

Three checks guard the grouping, and each runs whenever it has what it
needs:

- `.known` – every category written in `...` is a value of `.source`.
  This catches a mistyped category.

- `.exhaustive` – every value of `.source` falls in some group. This
  catches a value added to the data later, which would otherwise be
  dropped from the breakdown in silence.

- `.disjoint` – no value falls in more than one group.

The first two need `.source` and so default on when one is given; the
third needs only `...` and so is always on. Setting `.known` or
`.exhaustive` without a `.source` is an error, since there would be
nothing to check against.

`.disjoint` is not the route to overlapping groups. Build those
separately and collect them with
[`includes()`](https://r-logmu.logmu.org/reference/includes.md), which
checks nothing across its elements.

## Examples

``` r
# One group per status, written out.
categories(.i$status, "act", "def", "pen")
#> <includes[3]>
#>   status: act
#>   status: def
#>   status: pen

# A join, and named groups.
categories(.i$status, active = "act", other = c("def", "dep"))
#> <includes[2]>
#>   status: active
#>   status: other

# Every level of a factor, in its declared order.
scheme <- factor(c("A", "B", "A"), levels = c("A", "B", "C"))
categories(.i$scheme, .source = scheme)
#> <includes[3]>
#>   scheme: A
#>   scheme: B
#>   scheme: C

# A single group of two joined values.
category(.i$status, c("def", "dep"))
#> <indicator> (.i$status %in% c("def", "dep"))
```
