# Includes: a collection of includes

An `includes` is an optionally named list of
[include](https://r-logmu.logmu.org/reference/include.md)s. It is what
the `breakdown` argument of an analysis takes, and it is what every
plural constructor –
[`bands()`](https://r-logmu.logmu.org/reference/include.md),
[`ages()`](https://r-logmu.logmu.org/reference/include.md),
[`durations()`](https://r-logmu.logmu.org/reference/include.md),
[`periods()`](https://r-logmu.logmu.org/reference/include.md) – returns.

It is a separate type rather than an
[include](https://r-logmu.logmu.org/reference/include.md) with a length,
because giving `include` a length would touch every working include path
to buy generality only collections need.

Each element carries **two levels of naming**. Its own
[`names()`](https://rdrr.io/r/base/names.html) entry says which band it
is, and its `group_names()` entry says what was banded. A
fifteen-element breakdown of ages, periods and amounts can then tell a
chart which six are the ages, which a single flat name could not.

Elements need not be disjoint. Examining A/E on intersecting subsets is
legitimate, so **logmu** does not check for overlaps and does not treat
a breakdown as a partition.

## Usage

``` r
includes(...)

is_includes(x)

group_names(x)
```

## Arguments

- ...:

  `include` and `includes` objects to collect. An `includes` argument is
  flattened in. A name given to an `include` argument replaces its name;
  a name given to an `includes` argument replaces the group name of
  every element it contributes. An `includes` is subsettable with `[`
  and `[[`, which take an index and keep the labels of whatever they
  select.

  Both labels may be replaced. `names<-` takes exactly one name per
  include; `group_names<-` takes one per include or a single value for
  all of them.

- x:

  An `includes`.

## Value

`includes()` returns an `includes`.

`is_includes()` returns a scalar `logical`.

`group_names()` returns a `character` vector with one element per
include.

`names<-()` and `group_names<-()` return a new `includes`.

## Examples

``` r
std <- includes(
  ages(65, 95, by = 5),
  periods(2000, 2020, by = 10)
)
std
#> <includes[8]>
#>   age: 65-70
#>   age: 70-75
#>   age: 75-80
#>   age: 80-85
#>   age: 85-90
#>   age: 90-95
#>   period: 2000-2010
#>   period: 2010-2020
names(std)
#> [1] "65-70"     "70-75"     "75-80"     "80-85"     "85-90"     "90-95"    
#> [7] "2000-2010" "2010-2020"
group_names(std)
#> [1] "age"    "age"    "age"    "age"    "age"    "age"    "period" "period"

# Naming an includes argument renames the group of everything in it.
includes(cohort = ages(65, 95, by = 15))
#> <includes[2]>
#>   cohort: 65-80
#>   cohort: 80-95
```
