# Package index

## Package overview

- [`logmu`](https://r-logmu.logmu.org/reference/logmu-package.md)
  [`logmu-package`](https://r-logmu.logmu.org/reference/logmu-package.md)
  : Actuarial mortality experience analysis and modelling

## Experience and valuation data

- [`exp_data()`](https://r-logmu.logmu.org/reference/exp_and_val_data.md)
  [`val_data()`](https://r-logmu.logmu.org/reference/exp_and_val_data.md)
  :

  Create an `exp_data` or a `val_data`

- [`is_exp_data()`](https://r-logmu.logmu.org/reference/is_exp_and_val_data_type.md)
  [`is_val_data()`](https://r-logmu.logmu.org/reference/is_exp_and_val_data_type.md)
  :

  Is `x` an `exp_data` or a `val_data`?

- [`exp_start()`](https://r-logmu.logmu.org/reference/exp_data_info.md)
  [`exp_end()`](https://r-logmu.logmu.org/reference/exp_data_info.md) :

  Get experience period of `exp_data`

- [`as_at()`](https://r-logmu.logmu.org/reference/as_at.md) :

  Get the 'as at date' of a `val_data`

- [`val_from_exp_data()`](https://r-logmu.logmu.org/reference/val_from_exp_data.md)
  :

  Create a `val_data` from an `exp_data`

- [`` `$<-`( ``*`<exp_data>`*`)`](https://r-logmu.logmu.org/reference/assign_column.md)
  [`` `$<-`( ``*`<val_data>`*`)`](https://r-logmu.logmu.org/reference/assign_column.md)
  : Add or assign a data column

## Mortality

- [`is_mortality()`](https://r-logmu.logmu.org/reference/mortality.md)
  [`end_age()`](https://r-logmu.logmu.org/reference/mortality.md)
  [`log_mu()`](https://r-logmu.logmu.org/reference/mortality.md)
  [`mortality()`](https://r-logmu.logmu.org/reference/mortality.md) :
  Mortality-related functionality

- [`mortality_const()`](https://r-logmu.logmu.org/reference/mortality_const.md)
  :

  Create a constant `mortality`

- [`mortality_table()`](https://r-logmu.logmu.org/reference/mortality_table.md)
  :

  Create an `mortality` from a 2D age-period matrix of annual mortality
  rates

## Pronouns `.i` and `.t`

- [`pronoun_expressions()`](https://r-logmu.logmu.org/reference/pronoun_expressions.md)
  : Parse a pronoun expression into an AST
- [`variable()`](https://r-logmu.logmu.org/reference/variable.md)
  [`static_variable()`](https://r-logmu.logmu.org/reference/variable.md)
  [`is_variable()`](https://r-logmu.logmu.org/reference/variable.md)
  [`is_static_variable()`](https://r-logmu.logmu.org/reference/variable.md)
  [`value()`](https://r-logmu.logmu.org/reference/variable.md) :
  Variables: scalar functions of an individual and time
- [`indicator()`](https://r-logmu.logmu.org/reference/indicator.md)
  [`is_indicator()`](https://r-logmu.logmu.org/reference/indicator.md)
  [`logical_value()`](https://r-logmu.logmu.org/reference/indicator.md)
  : Indicators: a 0/1 function of an individual
- [`is_logmu_function()`](https://r-logmu.logmu.org/reference/is_logmu_function.md)
  : Test for a logmu function object

## Choosing which records to count

- [`band()`](https://r-logmu.logmu.org/reference/include.md)
  [`bands()`](https://r-logmu.logmu.org/reference/include.md)
  [`age()`](https://r-logmu.logmu.org/reference/include.md)
  [`ages()`](https://r-logmu.logmu.org/reference/include.md)
  [`duration()`](https://r-logmu.logmu.org/reference/include.md)
  [`durations()`](https://r-logmu.logmu.org/reference/include.md)
  [`period()`](https://r-logmu.logmu.org/reference/include.md)
  [`periods()`](https://r-logmu.logmu.org/reference/include.md)
  [`include()`](https://r-logmu.logmu.org/reference/include.md)
  [`is_include()`](https://r-logmu.logmu.org/reference/include.md)
  [`period_included()`](https://r-logmu.logmu.org/reference/include.md)
  : Includes: a single time interval per individual
- [`includes()`](https://r-logmu.logmu.org/reference/includes.md)
  [`is_includes()`](https://r-logmu.logmu.org/reference/includes.md)
  [`group_names()`](https://r-logmu.logmu.org/reference/includes.md) :
  Includes: a collection of includes
- [`category()`](https://r-logmu.logmu.org/reference/category.md)
  [`categories()`](https://r-logmu.logmu.org/reference/category.md) :
  Categories: breaking a population down by a categorical field

## A/E analysis

- [`aev()`](https://r-logmu.logmu.org/reference/aev.md) : Actual,
  expected and variance over experience data

- [`create_aev()`](https://r-logmu.logmu.org/reference/aev_properties.md)
  [`` `[[`( ``*`<aev>`*`)`](https://r-logmu.logmu.org/reference/aev_properties.md)
  [`` `$`( ``*`<aev>`*`)`](https://r-logmu.logmu.org/reference/aev_properties.md)
  [`` `[`( ``*`<aev>`*`)`](https://r-logmu.logmu.org/reference/aev_properties.md)
  : Properties of an aev

- [`is_aev()`](https://r-logmu.logmu.org/reference/is_aev.md) :

  Is an object an `aev`?

- [`as.data.frame(`*`<aev>`*`)`](https://r-logmu.logmu.org/reference/as.data.frame.aev.md)
  :

  Coerce an `aev` to a data frame

## Charting an A/E

- [`autoplot()`](https://r-logmu.logmu.org/reference/aev_plot.md)
  [`plot(`*`<aev>`*`)`](https://r-logmu.logmu.org/reference/aev_plot.md)
  :

  Chart an `aev`

- [`theme_aev()`](https://r-logmu.logmu.org/reference/theme_aev.md) :
  The chart theme, and how to resize a chart

## Running several analyses at once

- [`batch()`](https://r-logmu.logmu.org/reference/batch.md) : Run
  several analyses over one dataset in a single pass
- [`settings()`](https://r-logmu.logmu.org/reference/settings.md)
  [`is_settings()`](https://r-logmu.logmu.org/reference/settings.md)
  [`print(`*`<logmu_settings>`*`)`](https://r-logmu.logmu.org/reference/settings.md)
  : Settings for an analysis

## Vector operations

- [`vec_active_lanes()`](https://r-logmu.logmu.org/reference/vec_active.md)
  [`vec_active_tier()`](https://r-logmu.logmu.org/reference/vec_active.md)
  : SIMD active lanes and tier
- [`vec_neg()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_exp()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_expm1()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_log()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_log1p()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_m_from_q()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_add()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_sub()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_mul()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_div()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_pow()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_min()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_max()`](https://r-logmu.logmu.org/reference/vec_ops.md)
  [`vec_clamp()`](https://r-logmu.logmu.org/reference/vec_ops.md) :
  Vector ops
