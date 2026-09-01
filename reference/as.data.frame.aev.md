# Coerce an `aev` to a data frame

One row per record, carrying the labels, the raw triple and the three
calculated properties a chart reads. This is the frame **logmu**'s own
plotting is built on, and the way to take an `aev` into **ggplot2**,
**dplyr** or anything else that works on data frames.

## Usage

``` r
# S3 method for class 'aev'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `aev` object.

- row.names:

  Row names for the result, or `NULL` for the default.

- optional:

  Ignored. Present because the generic has it; every column name here is
  already syntactic.

- ...:

  Ignored.

## Value

A `data.frame` with one row per record of `x`.

## Columns

|  |  |
|----|----|
| Column | Contents |
| `name` | the element label, [`names()`](https://rdrr.io/r/base/names.html) on the `aev` |
| `group` | the group label, [`group_names()`](https://r-logmu.logmu.org/reference/includes.md) on the `aev` |
| `A`, `E`, `V` | the triple itself |
| `A_over_E` | \\A / E\\ |
| `log_A_over_E_stddev` | \\\sqrt{V} / E\\ |
| `deviance_residual` | see [aev_properties](https://r-logmu.logmu.org/reference/aev_properties.md) |

The calculated columns are named after the properties that produce them,
so `frame$A_over_E` and `aev$A_over_E` are the same word for the same
quantity. The five remaining properties are left out because each is a
line of arithmetic on `A`, `E` and `V`, and naming them here would fix
five more column names for no gain.

`name` and `group` are always present, and are `NA` on an `aev` that
carries no labels. The columns therefore depend on the input's type and
never on its values, which is what lets frames from a broken-down `aev`
and an ungrouped one be stacked with
[`rbind()`](https://rdrr.io/r/base/cbind.html).

## Examples

``` r
aev <- create_aev(A = c(1100, 40), E = c(1000, 50), V = c(2500, 125))
names(aev) <- c("65-70", "70-75")
group_names(aev) <- "age"

as.data.frame(aev)
#>    name group    A    E    V A_over_E log_A_over_E_stddev deviance_residual
#> 1 65-70   age 1100 1000 2500      1.1           0.0500000         1.9679833
#> 2 70-75   age   40   50  125      0.8           0.2236068        -0.9270417

# An unlabelled `aev` gives the same columns, with the labels NA.
as.data.frame(create_aev(A = 1100, E = 1000, V = 2500))
#>   name group    A    E    V A_over_E log_A_over_E_stddev deviance_residual
#> 1 <NA>  <NA> 1100 1000 2500      1.1                0.05          1.967983
```
