# Variables: scalar functions of an individual and time

A `variable` is a scalar, real-valued function \\f(i, t)\\ written with
the data pronouns `.i$field`, `.t`, `.b` and `.x`. A `static_variable`
is a variable that does not depend on time – a function of the
individual's facts alone, \\f(i)\\.

Construct one from a pronoun expression or a `~` formula. The
constructors return the **most specific** type they can prove: a
time-invariant expression becomes a `static_variable`, and one that is
also structurally logical becomes an `indicator`. `static_variable()`
additionally asserts that the expression does not use `.t` (and fails if
it does);
[`indicator()`](https://r-logmu.logmu.org/reference/indicator.md)
asserts a \\\\0,1\\\\ value. So `variable(.i$pension)` is a
`static_variable` and `variable(.i$pension > 0)` is an `indicator`.

`value()` evaluates the variable against a single individual's facts
`.i` (a named list) and a time vector `.t`. It is the plain-R reference
path for testing and understanding – not the performance path. For a
`static_variable`, `.t` may be omitted.

## Usage

``` r
variable(expr)

static_variable(expr)

is_variable(x)

is_static_variable(x)

value(x, .i, .t = NULL, ...)

# Default S3 method
value(x, .i, .t = NULL, ...)

# S3 method for class 'variable'
value(x, .i, .t = NULL, ...)
```

## Arguments

- expr:

  A pronoun expression, optionally written as a `~` formula.

- x:

  A `variable`.

- .i:

  A named list of scalar facts.

- .t:

  A vector of `datey` (omit for a time-invariant variable).

- ...:

  Unused.

## Value

`variable()` and `static_variable()` return a `variable` (the latter
also classed `static_variable`).

`is_variable()` and `is_static_variable()` return a scalar `logical`.

`value()` returns a vector aligned to `.t`.

## Examples

``` r
amounts <- static_variable(.i$pension)
value(amounts, .i = list(pension = 1000))
#> [1] 1000

v <- variable(.i$pension * 2)
value(v, .i = list(pension = 1000))
#> [1] 2000
```
