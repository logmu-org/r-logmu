# Properties of an aev

An `aev` comprises the following properties:

|          |                                       |
|----------|---------------------------------------|
| Property | Formula                               |
| `A`      | Actual deaths weighted by \\w\\       |
| `E`      | 'Expected' deaths weighted by \\w\\   |
| `V`      | 'Expected' deaths weighted by \\w^2\\ |

where \\w\\ is an arbitrary (non-negative) weight used to calculate
them.

`create_aev()` refuses a triple that breaks any of the following, so an
`aev` you build by hand is guaranteed to satisfy them: (a) none of `A`,
`E` or `V` are negative, (b) `E` and `V` are either both zero or both
non-zero, and (c) `A` cannot be non-zero if `E` and `V` are zero.

A **computed** `aev` is not checked against them, and (b) and (c) can
fail when a mortality underflows:
[`exp()`](https://rdrr.io/r/base/Log.html) reaches exactly zero at a log
mortality of about -746, so `E` and `V` are then zero while `A` still
counts the deaths that happened. That is a legitimate statement about an
impossible model, so it is returned rather than raised – an analysis
running for hours, or a
[`batch()`](https://r-logmu.logmu.org/reference/batch.md) whose other
elements are sound, must not be lost to one underflowing cell. The
calculated properties below take their ordinary IEEE values there, so
`A_over_E` reads `Inf` and `deviance_residual` reads `NaN`.

These columns cannot be modified individually.

In typical use, `aev`s are created by **logmu** analytic functions. A
`create_aev` function is provided for testing and illustration.

`A`, `E` or `V` have the following statistical properties *if the
mortality used to calculate them is correct*:

\$\$\mathbb{E}\\(A - E) = 0\$\$ \$\$\mathrm{Var}(A-E) =
\mathbb{E}\\V\$\$

An `aev` provides the following *calculated* properties:

|  |  |
|----|----|
| Property | Formula |
| `A_minus_E` | \\A - E\\ |
| `A_minus_E_stddev` | \\\sqrt{V}\\ |
| `A_over_E` | \\A / E\\ |
| `log_A_over_E` | \\\log(A / E)\\ |
| `log_A_over_E_stddev` | \\\sqrt{V} / E\\ |
| `log_A_over_E_95pc` | \\k\sqrt{V} / E\\ where \\k \approx 1.96\\ – used by **logmu** to display 95% confidence intervals |
| `Pearson_residual` | \\(A - E) / \sqrt{V}\\ |
| `deviance_residual` | \\\mathrm{sign}(A-E)\sqrt{2E/V\cdot\Big\[A\cdot\log(A/E)-(A-E)\Big\]}\\ |

`aev`s can be added, which simply means adding their A, E and V
components.

The addition of `aev`s is statistically legitimate *provided they relate
to independent experience data*, i.e. data that does not intersect by
time *and* individual. This means that it *is* legitimate to add `aev`s
for the same individual provided they relate to non-overlapping time
periods, or for overlapping time periods provided they relate to
different individuals.

An `aev` is a **vector of records**:
[`length()`](https://rdrr.io/r/base/length.html) reports how many, and
`[` subsets records, keeping `A`, `E` and `V` together.

When an `aev` comes from a breakdown it also carries two levels of
labelling, [`names()`](https://rdrr.io/r/base/names.html) and
[`group_names()`](https://r-logmu.logmu.org/reference/includes.md),
exactly as the
[includes](https://r-logmu.logmu.org/reference/includes.md) it came from
did. Both are absent from a hand-built `aev`.

Two `aev`s may be added, which adds `A`, `E` and `V` record by record.
The labels travel with the sum. Adding record by record asserts that the
records correspond, so labels that disagree say they do not, and the
addition is refused. A missing set of labels is not a disagreement:
adding a hand-built `aev` to a labelled one keeps the labels.

Both labels may be replaced. `names<-` takes exactly one label per
record; `group_names<-` takes one per record or a single value for all
of them.

## Usage

``` r
create_aev(A, E, V)

# S3 method for class 'aev'
x[[i]]

# S3 method for class 'aev'
x$name

# S3 method for class 'aev'
x[i]
```

## Arguments

- A, E, V:

  Columns used to construct an `aev`.

- x:

  An `aev` object.

- i, name:

  The `aev` property being requested.

## Value

An `aev`.

## Examples

``` r
aev <- create_aev(A = c(1100, 0), E = c(1000, 1), V = c(2500, 1))
aev
#> <aev[2]>
#>      A    E    V A/E  95% conf dev resid
#> 1 1100 1000 2500 1.1 0.0979982  1.967983
#> 2    0    1    1 0.0 1.9599640 -1.414214

# Data properties:
aev$A
#> [1] 1100    0
aev$E
#> [1] 1000    1
aev$V
#> [1] 2500    1

# Calculated properties:
aev$A_minus_E
#> [1] 100  -1
aev$A_minus_E_stddev
#> [1] 50  1
aev$A_over_E
#> [1] 1.1 0.0
aev$log_A_over_E
#> [1] 0.09531018       -Inf
aev$log_A_over_E_stddev
#> [1] 0.05 1.00
aev$log_A_over_E_95pc
#> [1] 0.0979982 1.9599640
aev$Pearson_residual
#> [1]  2 -1
aev$deviance_residual
#> [1]  1.967983 -1.414214

# Addition:
aev + create_aev(A = c(100, 1), E = c(200, 10), V = c(1500, 10))
#> <aev[2]>
#>      A    E    V        A/E  95% conf dev resid
#> 1 1200 1200 4000 1.00000000 0.1032992  0.000000
#> 2    1   11   11 0.09090909 0.5909514 -3.899258

# A vector of records: `length()` counts records and `[` subsets them.
length(aev)
#> [1] 2
aev[1]
#> <aev[1]>
#>      A    E    V A/E  95% conf dev resid
#> 1 1100 1000 2500 1.1 0.0979982  1.967983

# Labels, as a breakdown would supply them:
names(aev) <- c("65-70", "70-75")
group_names(aev) <- "age"
aev
#> <aev[2]>
#>   group  name    A    E    V A/E  95% conf dev resid
#> 1   age 65-70 1100 1000 2500 1.1 0.0979982  1.967983
#> 2   age 70-75    0    1    1 0.0 1.9599640 -1.414214
names(aev[2])
#> [1] "70-75"

# Validity:
# (a) none of `A`, `E` or `V` are negative,
try(create_aev(A = -1, E = 2, V = 3))
#> Error : A cannot be negative.
# (b) `E` and `V` are either both zero or both non-zero, and
try(create_aev(A = 0, E = 0, V = 3))
#> Error : V cannot be non-zero if E is zero.
# (c) `A` cannot be non-zero if `E` and `V` are zero.
try(create_aev(A = 1, E = 0, V = 0))
#> Error : A cannot be non-zero if E and V are zero.
```
