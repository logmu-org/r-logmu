# Settings for an analysis

A `settings` object carries the choices an analysis needs that are not
part of the question being asked: the overdispersion assumed, and the
time scale the integration uses.

Build one once and pass it wherever it is wanted. Analytic functions
take it as a `settings` argument, and it is found among the arguments by
its class, so it needs no name and no position.

## Usage

``` r
settings(overdispersion, time_scale = default_time_scale)

is_settings(x)

# S3 method for class 'logmu_settings'
print(x, ...)
```

## Arguments

- overdispersion:

  A single positive number. Required.

- time_scale:

  The integration interval, as a `durationy` or a number of years. One
  of 1, 1/4, 1/12 or 1/60.

- x:

  An object.

- ...:

  Ignored.

## Value

`settings()` returns a `logmu_settings` object.

`is_settings()` returns a scalar `logical`.

`x`, invisibly.

## Overdispersion

Overdispersion, written \\\Omega\\, is defined by
\$\$\mathrm{Var}(\mathrm{A}w - \mathrm{E}w) =
\Omega\\\mathbb{E}\\\mathrm{E}w^2\$\$ so it is the factor by which the
variance of experience exceeds what independent deaths under a
deterministic mortality would give. It is **required**, and has no
default anywhere in **logmu**. Failing to allow for it does not make
results neutral – it understates uncertainty by \\\sqrt\Omega\\ and
selects overfitted models, so there is no safe value to assume on a
user's behalf.

Values between 2 and 3 are usual for pensions longevity work, with
higher values making model selection more resistant to overfitting.

## Time scale

The time scale is the width of one numerical integration interval. It
may be given as a `durationy` or as a number of years, and must be one
of 1, 1/4, 1/12 or 1/60 of a year.

Those four are the intervals that, together with their halves, are a
whole number of clicks, so every sample point lands exactly on the click
grid. They also nest: refining from 1/4 to 1/12 to 1/60 keeps every
sample already taken and adds more between them.

Smaller is more accurate and costs proportionally more. The default of a
quarter year is short enough to sample an annual mortality table
sensibly.

## Examples

``` r
settings(overdispersion = 2)
#> <settings>
#>   overdispersion: 2
#>   time scale:     1/4 year

settings(overdispersion = 2.5, time_scale = 1 / 12)
#> <settings>
#>   overdispersion: 2.5
#>   time scale:     1/12 year

# A durationy says the same thing.
settings(overdispersion = 2.5, time_scale = datey::durationy(1 / 12))
#> <settings>
#>   overdispersion: 2.5
#>   time scale:     1/12 year

# Overdispersion is required.
try(settings(time_scale = 1))
#> Error : `overdispersion` is required and has no default. Values between 2 and 3 are usual for pensions longevity work.

# And the time scale must be one of the four.
try(settings(overdispersion = 2, time_scale = 0.5))
#> Error : `time_scale` must be one of 1, 1/4, 1/12, 1/60 of a year.
```
