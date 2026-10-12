# The name of an object

Many **logmu** objects can carry a name of their own: a `mortality`, a
[`model()`](https://r-logmu.logmu.org/reference/model.md), an
[include](https://r-logmu.logmu.org/reference/include.md) or
[indicator](https://r-logmu.logmu.org/reference/indicator.md) and each
term of a
[`covariates()`](https://r-logmu.logmu.org/reference/covariates.md).
`name()` reads it and `name<-` sets or, with `NULL`, clears it. Where an
object has one, the name is what **logmu** prints for it, so a reference
mortality shows as `S4PMA` rather than as `<mortality_table>`, and a
model handed to
[`compare_models()`](https://r-logmu.logmu.org/reference/compare_models.md)
without a list name is labelled by its own.

An include carries a second level, its `group_name`, which says which
breakdown it belongs to: `age` for `[60, 70)`. Over a whole
[`includes()`](https://r-logmu.logmu.org/reference/includes.md) the two
levels are read as vectors by
[`names()`](https://rdrr.io/r/base/names.html) and
[`group_names()`](https://r-logmu.logmu.org/reference/includes.md).

## Usage

``` r
name(x)

name(x) <- value

group_name(x)

group_name(x) <- value
```

## Arguments

- x:

  An object.

- value:

  A name, or `NULL` to clear it.

## Value

`name()` and `group_name()` return a single string, or `NULL` where the
object has none.

`name<-` and `group_name<-` return `x` with the name set.

## What a name may be

A name is a single string that is not blank, has no leading, trailing or
consecutive spaces, and holds no control character such as a tab or a
newline. The rule is the same for every name in the package, whether you
wrote it or **logmu** generated it: `[60, 70)`, `>= 70` and `a+b` are
all names it makes for itself, and they follow the rule too.

## Examples

``` r
base <- mortality_const(log_mu = -4, name = "flat")
name(base)
#> [1] "flat"

m <- model(ref_mortality = base, covariates = covariates(level = 1))
name(m) <- "one level"
m
#> <model: 1 term>
#>   name:      one level
#>   reference: flat
#>   [1] level

std <- ages(65, 95, by = 15)
name(std[[1]])
#> [1] "65-80"
group_name(std[[1]])
#> [1] "age"
```
