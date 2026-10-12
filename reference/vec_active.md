# SIMD active lanes and tier

`vec_active_tier()` names the SIMD kernel tier the package resolved when
its shared library was loaded, and `vec_active_lanes()` gives the number
of `double` lanes that tier works on.

The tier is chosen from the instruction sets the CPU reports: `"avx512"`
(8 lanes), `"avx2"` (4 lanes), or `"baseline"` (2 lanes) where neither
of those is available.

## Usage

``` r
vec_active_lanes()

vec_active_tier()
```

## Value

A scalar `integer` representing the number of active SIMD lanes.

A scalar `character` describing the active SIMD tier.

## Forcing a lower tier

Setting the environment variable `LOGMU_TIER` to `"baseline"`, `"avx2"`
or `"avx512"` before the package is loaded places a ceiling on the tier
that may be selected. It is there for comparing tiers in a benchmark,
and for reproducing a fault that only one of them shows.

The ceiling can only ever lower the tier, never raise it. A CPU that
does not support the tier asked for still falls back to the best one it
does support, because the processor feature checks apply exactly as they
otherwise would.

An unrecognised value places no ceiling at all, and nothing is reported,
so call `vec_active_tier()` to confirm which tier is in use rather than
assuming the request was honoured.

The variable is read once, when the shared library is loaded. Changing
it later has no effect for the rest of the session.
