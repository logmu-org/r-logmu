# Run several analyses over one dataset in a single pass

`batch()` takes named calls to **logmu** analyses and runs them
together: their specifications are compiled and scheduled as one
crossing, the data is read once, and each analysis comes back as
whatever it would have returned on its own.

## Usage

``` r
batch(
  ...,
  .exp_data = NULL,
  .include = NULL,
  .weight = NULL,
  .val_similarity = NULL,
  .val_distance = NULL,
  .breakdown = NULL,
  .settings = NULL,
  .overdispersion = NULL,
  .time_scale = NULL,
  .threads = cpp_veil_default_threads()
)
```

## Arguments

- ...:

  Named calls to **logmu** analyses, currently
  [`aev()`](https://r-logmu.logmu.org/reference/aev.md). Each name
  becomes a name in the result and may not begin with a dot.

- .exp_data, .include, .weight, .val_similarity, .val_distance,
  .breakdown, .settings, .overdispersion, .time_scale:

  Defaults for the analyses, each standing in for the argument of the
  same name wherever an analysis does not supply its own.
  `.val_similarity` and `.val_distance` are two spellings of one
  quantity, so only one may be given and an analysis naming either of
  them takes neither default.

- .threads:

  Worker threads for the whole batch. `0` asks for as many as the
  machine reports.

## Value

A named list holding each analysis's own result, in the order written.
Nothing else – no class, no attributes.

## Settings are defaults, not overrides

Every dotted argument supplies a default for the analyses inside. An
analysis naming its own value wins; the batch's applies when it has
none; it is an error only when neither supplies one. So `overdispersion`
remains required without **logmu** ever assuming a value for it.

The dot marks a batch setting, and it is what stops a setting colliding
with an analysis you have named. **An element of a batch may not be
named with a leading dot**, which reserves the whole dotted namespace so
that settings added in future cannot break existing code.

`.threads` is the exception to the default rule: it belongs to the run
rather than to any one analysis, so the batch's value is used throughout
and a `threads` argument inside a batched call is ignored. It cannot
change an answer, only a duration.

## What a batch may not do

**No element may use another element's result.** Every specification is
compiled before a record is read and the results exist only once the
pass is over, so a dependency between elements could not be honoured.
`batch()` refuses a call that mentions a sibling's name rather than
letting it resolve silently against something of the same name in your
workspace.

A batch runs over **one experience dataset**. Analyses may differ in
`time_scale`, which is worth doing to see whether the integration
interval moves the answer; a batch then makes one pass per distinct
scale.

## Examples

``` r
data <- exp_data(
  list(
    birth     = datey::datey(c(1945, 1950, 1955)),
    pension   = c(5000, 12000, 30000),
    male      = c(TRUE, FALSE, TRUE),
    E2R_start = datey::datey(c(2015, 2015, 2015)),
    E2R_end   = datey::datey(c(2020, 2020, 2018)),
    E2R_died  = c(FALSE, FALSE, TRUE)
  ),
  exp_start = datey::datey(2015),
  exp_end   = datey::datey(2020)
)

b <- batch(
  .exp_data       = data,
  .overdispersion = 2,
  .weight         = .i$pension,
  light = aev(mortality = mortality_const(log_mu = -4.5)),
  heavy = aev(mortality = mortality_const(log_mu = -4.0), overdispersion = 1)
)

b$light
#> <aev[1]>
#>       A        E        V      A/E 95% conf dev resid
#> 1 30000 1944.074 78762785 15.43151 8.947377  1.633256
b$heavy
#> <aev[1]>
#>       A        E        V      A/E 95% conf dev resid
#> 1 30000 3205.237 64928940 9.359683 4.927279  1.994644
```
