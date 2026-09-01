# Create an `exp_data` or a `val_data`

Create a experience or valuation data from a list of columns, a
`data.frame` or `tibble` (or similar).

The columns must be atomic, the same length and uniquely named.

There *must* be a `birth` date column.

A `count` column is optional. If present it must be numeric,
non-negative and finite. Rows with a nil count are dropped.

Experience data must also contain (and valuation data must *not*
contain)

- `E2R_start` and `E2R_end` date columns, and

- an `E2R_died` logical column.

These special columns cannot contain NA and, if they are dates, they
must be `datey`.

The experience data period and individual E2Rs must be *proper*:

- at the dataset level, `exp_start` cannot be after `exp_end`, and

- at the record level, `E2R_start` cannot be after `E2R_end` and if they
  are the same then `E2R_died` must be false.

The experience data is trimmed to the experience period, i.e. it is
guaranteed that `E2R_start` is not before `exp_start` and `E2R_end` is
not after `exp_end`. Records with empty E2Rs after trimming are dropped.

## Usage

``` r
exp_data(
  columns,
  exp_start,
  exp_end,
  on_unreadable = c("error", "warn", "ignore")
)

val_data(columns, as_at, on_unreadable = c("error", "warn", "ignore"))
```

## Arguments

- columns:

  The list of columns (which includes things like `data.frame` and
  `tibble`) to convert to an `exp_data` or `val_data`.

- exp_start:

  The start of the overall experience period (inclusive). Must be before
  `exp_end`.

- exp_end:

  The end of the overall experience period (exclusive). Must be after
  `exp_start`.

- on_unreadable:

  What to do about columns whose type **logmu** cannot read (anything
  other than `double`, `integer`, `logical`, `character`, `factor`,
  `datey` or `durationy`): `"error"` (the default) rejects them,
  `"warn"` keeps them with a warning, and `"ignore"` keeps them
  silently. The special columns (`birth`, the `E2R_*` columns, `count`)
  are always checked.

- as_at:

  The 'as at' date of the valuation data.

## Value

An `exp_data` or `val_data`.
