# Is an object an `aev`?

Tests whether `x` is an `aev` – the three-field record type
[`aev()`](https://r-logmu.logmu.org/reference/aev.md) returns and
[`create_aev()`](https://r-logmu.logmu.org/reference/aev_properties.md)
builds. See
[aev_properties](https://r-logmu.logmu.org/reference/aev_properties.md).

## Usage

``` r
is_aev(x)
```

## Arguments

- x:

  An object to test.

## Value

A single `TRUE` or `FALSE`. Never raises an error, whatever `x` is.

## Examples

``` r
is_aev(create_aev(A = 1, E = 2, V = 3))
#> [1] TRUE
is_aev(1:3)
#> [1] FALSE
```
