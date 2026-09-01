# Mortality-related functionality

A **logmu** mortality is an object class that defines the annual
mortality rate \\\mu\_{it}\\ for individual \\i\\ at time \\t\\, where
\\i\\ contains all the available invariant information relating to the
relevant individual.

All **logmu** actually works in terms of \\\log\mu\_{it}\\ for reasons
set out TBC.

Access to \\\log\mu\_{it}\\ for a mortality object is implicit using
`.i` and `.t` pronouns.

For testing and understanding, you can access the \\\log\mu\_{it}\\
calculation by defining your own `.i` and `.t` variables and then
calling `log_mu()` on a mortality. Do *not* use this for performant
scenarios.

These are the core types of mortality object:

|  |  |  |
|----|----|----|
| Type | **datey** function | Notes |
| Constant | [`mortality_const()`](https://r-logmu.logmu.org/reference/mortality_const.md) |  |
| Age-period | [`mortality_table()`](https://r-logmu.logmu.org/reference/mortality_table.md) | Mortality rates that are smooth at an annual scale |
| Model | TBC |  |

You can build on these:

- Select between tables using TBC.

- Scale mortality rates using TBC.

- Provide sub-annual mortality variation using TBC.

LINK TO OTHER PACKAGES TO GET STANDARD TABLES

Finally there are a couple of mortality-related helper functions:

- `is_mortality(x)` tests whether `x` is a `mortality`.

- `end_age(x)` gets the end age of the `mortality`, i.e. the age after
  which everyone is assumed to be dead. (This is provided for valuation
  calculation and is currently ignored for AEV calculations.)

Builds a `mortality` from a pronoun expression for \\\log\mu\\. Concrete
`mortality` objects may appear in the expression – each contributes its
own \\\log\mu\\ at the individual and time – so you can

- select between mortalities with
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html) (the condition must
  be time-invariant, as for any conditional),

- scale or adjust a mortality (adding in \\\log\mu\\ is multiplying
  \\\mu\\), or

- write a closed-form law directly from the pronouns.

A bare reference to a concrete `mortality` is returned unchanged.

## Usage

``` r
is_mortality(x)

end_age(x)

log_mu(x, .i, .t)

# Default S3 method
log_mu(x, .i, .t)

# S3 method for class 'mortality_const'
log_mu(x, .i, .t)

mortality(expr)

# S3 method for class 'mortality_expr'
log_mu(x, .i, .t)

# S3 method for class 'mortality_table'
log_mu(x, .i, .t)
```

## Arguments

- x:

  The `mortality` object.

- .i:

  A list of named scalar arguments.

- .t:

  A vector of `datey` representing time.

- expr:

  A pronoun expression for \\\log\mu\\, optionally a `~` formula.

## Value

`is_mortality()` returns a scalar `logical`.

`end_age()` returns a scalar `durationy`.

A vector of \\\log\mu\\ values at `.t`.

A `mortality`.

## Examples

``` r
mortality <- mortality_table(x0 = 70, t0 = 2020, q = matrix(0.01))
.t <- datey::datey(2020)
.i <- list(birth = datey::datey(1950))
log_mu(mortality, .i, .t) # log(-log(1 - 0.01)) = -4.600149
#> [1] -4.600149
base <- mortality_const(log_mu = -4)
mortality(base + 0.05)   # scale mu up by exp(0.05)
#> $ast
#> call `+`
#>   obj <mortality_const>
#>   lit 0.05
#> 
#> attr(,"class")
#> [1] "mortality_expr" "mortality"      "logmu_function"
```
