# Indicators: a 0/1 function of an individual

An `indicator` is a time-invariant function of an individual's facts
that takes the value 0 or 1. It is the one type that is both a
`variable` (a \\\\0,1\\\\ weight, for which \\E = V\\) and an `include`
(a subset: TRUE selects all of time, FALSE selects nothing).

Build one from a pronoun expression or a `~` formula. The expression
must not use time `.t`. A logical result is always a valid indicator; a
numeric result is accepted but is only checked to be \\\\0,1\\\\ when
evaluated.

`logical_value()` evaluates the indicator against an individual's facts
`.i` and returns a `logical`. Because an indicator is also an `include`,
[`period_included()`](https://r-logmu.logmu.org/reference/include.md)
works too, returning all-of-time or the empty interval.

## Usage

``` r
indicator(expr)

is_indicator(x)

logical_value(x, .i)

# Default S3 method
logical_value(x, .i)

# S3 method for class 'indicator'
logical_value(x, .i)
```

## Arguments

- expr:

  A pronoun expression, optionally written as a `~` formula.

- x:

  An `indicator`.

- .i:

  A named list of scalar facts.

## Value

`indicator()` returns an `indicator`.

`is_indicator()` returns a scalar `logical`.

`logical_value()` returns a `logical`.

## Examples

``` r
male <- indicator(.i$sex == "male")
logical_value(male, .i = list(sex = "male"))
#> [1] TRUE
period_included(male, .i = list(sex = "female"))   # empty interval
#> [1] <NA>
```
