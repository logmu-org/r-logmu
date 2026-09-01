# Create a `val_data` from an `exp_data`

Creates a `val_data` from those individuals in an `exp_data` who were
alive at the specified date.

## Usage

``` r
val_from_exp_data(exp_data, as_at = NULL)
```

## Arguments

- exp_data:

  The `exp_data`.

- as_at:

  The as at date (a `datey`). If omitted, the *end* of the experience
  period is used.

## Value

A `val_data`.
