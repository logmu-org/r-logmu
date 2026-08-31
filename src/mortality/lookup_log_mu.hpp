// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <cmath>
#include <algorithm>
#include "datey.h"
#include "math_extra.h"

namespace mortality
{

inline double get_value_from_tx_array(
  const double* p, int n_t, int n_x,
  int t_index, int x_index
)
{
  t_index = std::clamp(t_index, 0, n_t - 1);
  x_index = std::clamp(x_index, 0, n_x - 1);

  // We've checked that n_t * n_x <= 2^31 - 1 so this is safe.
  return p[t_index * n_x + x_index];
}


// Look up a single interpolated log mu(b,t) value from rectangular tx array.
// Parameters are `int`s because
// (a) R matrix dimensions are `int`s and
// (b) `datey` is `int` under the covers.
// Correctness, not performance is the primary issue here.
// (For performance we'll use vectorised operations.)
inline double get_interpolated_log_mu_bt_from_annual_tx_array(
  const double* p, int n_t, int n_x,
  int t0_clicks, int x0_clicks,
  int b_clicks, int t_clicks
)
{
  int b0_clicks = t0_clicks - x0_clicks;

  auto b_split = divmod(b_clicks - b0_clicks, ClicksPerYear);
  int b_index_0 = b_split.div;
  double b_epsilon = b_split.mod / (double)ClicksPerYear;

  auto t_split = divmod(t_clicks - t0_clicks, ClicksPerYear);
  int t_index_0 = t_split.div;
  double t_epsilon = t_split.mod / (double)ClicksPerYear;

  double log_mu_b0_t0 = get_value_from_tx_array(p, n_t, n_x, t_index_0, t_index_0 - b_index_0);
  double log_mu_b1_t0 = get_value_from_tx_array(p, n_t, n_x, t_index_0, t_index_0 - (b_index_0 + 1));
  double log_mu_b0_t1 = get_value_from_tx_array(p, n_t, n_x, t_index_0 + 1, t_index_0 + 1 - b_index_0);
  double log_mu_b1_t1 = get_value_from_tx_array(p, n_t, n_x, t_index_0 + 1, t_index_0 - b_index_0);

  double log_mu_b0_t = std::lerp(log_mu_b0_t0, log_mu_b0_t1, t_epsilon);
  double log_mu_b1_t = std::lerp(log_mu_b1_t0, log_mu_b1_t1, t_epsilon);

  return std::lerp(log_mu_b0_t, log_mu_b1_t, b_epsilon);
}

} // namespace mortality

