// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cmath>
#include <limits>

namespace veil
{

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
constexpr double Infinity = std::numeric_limits<double>::infinity();
inline bool isIntegral(double val) noexcept { return std::trunc(val) == val; }

}