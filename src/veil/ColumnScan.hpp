// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include "datey.h" // isValidDatey, isValidDurationy
#include "veil/ColumnView.hpp"
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"
#include "veil/TypeSpecificConstraint.hpp"
#include "veil/double_properties.hpp"

namespace veil
{

// What a column scan reads off one column's values: what the folding and interval passes need to
// specialise a calculation to its data. R-free -- it reads a ColumnView's raw values, never a SEXP.
struct TypeWithConstraints final
{
    TypeFull type;
    bool hasNAs;
    bool hasValues;         // At least one non-NA value.
    bool valuesAreConstant; // All non-NA values equal.
    TypeSpecificConstraint specific;
};

namespace detail
{

// Invalidity is by year range, from the datey framework (src/datey.h): a datey is a year in
// [1000, 3000] and a durationy a span in [-2000, +2000] years. This is the platform-independent test
// -- R's INT_MIN NA sentinel is an R-ism the core must not lean on -- and it catches an R integer NA
// anyway, since INT_MIN falls far outside both ranges.
//
// A CATEGORY IS NOT AMONG THEM, and needs no test here. Its storage is `int`, which has no
// platform-independent NA representation, so a category has no missing value: a factor column
// carrying one is refused when it is read, and every code that reaches this is a genuine index into
// the crossing's StringMapping.
inline bool isValidClick(int clicks, Type type) noexcept
{
    switch (type)
    {
        case Type::Datey: return isValidDatey(clicks);
        case Type::Durationy: return isValidDurationy(clicks);
        default: return true;
    }
}

inline TypeWithConstraints scanDoubles(const double* values, size_t count)
{
    bool hasNAs = false;
    double min = Infinity;
    double max = -Infinity;
    bool allIntegral = true;
    for (size_t i = 0; i < count; ++i)
    {
        const double value = values[i];
        if (std::isnan(value)) { hasNAs = true; continue; } // A NaN is an NA; veil draws no distinction.
        min = std::min(min, value);
        max = std::max(max, value);
        allIntegral &= isIntegral(value);
    }
    const bool hasValues = min <= max;
    return TypeWithConstraints{TypeFull::createDouble(), hasNAs, hasValues,
                             hasValues && min == max, DoubleConstraint(min, max, allIntegral)};
}

inline TypeWithConstraints scanInts(const int* values, size_t count, TypeFull type)
{
    bool hasNAs = false;
    int min = std::numeric_limits<int>::max();
    int max = std::numeric_limits<int>::min();
    for (size_t i = 0; i < count; ++i)
    {
        const int value = values[i];
        if (!isValidClick(value, type.type)) { hasNAs = true; continue; }
        min = std::min(min, value);
        max = std::max(max, value);
    }
    const bool hasValues = min <= max;
    return TypeWithConstraints{type, hasNAs, hasValues, hasValues && min == max, IntConstraint(min, max)};
}

inline TypeWithConstraints scanBools(const bool* values, size_t count)
{
    bool minValue = true;
    bool maxValue = false;
    for (size_t i = 0; i < count; ++i)
    {
        minValue = minValue && values[i];
        maxValue = maxValue || values[i];
    }
    const bool hasValues = count > 0;
    // A bool column carries no NA: the reader has already refused a missing logical.
    return TypeWithConstraints{TypeFull::createBool(), false, hasValues,
                             hasValues && minValue == maxValue, BoolConstraint(minValue, maxValue)};
}

} // namespace detail

// Reads one column's constraints. datey, durationy and category share the integer scan.
//
// Text is not a column type: a column of strings reaches the core already encoded to category codes,
// so Text survives only as the type of a literal a front end wrote, and only until passEncodeText.
inline TypeWithConstraints scanColumn(const ColumnView& column)
{
    switch (column.type.type)
    {
        case Type::Bool: return detail::scanBools(column.bools(), column.count);
        case Type::Double: return detail::scanDoubles(column.doubles(), column.count);
        case Type::Datey:
        case Type::Durationy:
        case Type::Category: return detail::scanInts(column.ints(), column.count, column.type);
        case Type::Text: throw std::runtime_error("veil: text is not a column type; a text column is held as category codes.");
        case Type::DateyInterval: throw std::runtime_error("veil: a datey_interval is not a column type.");
    }
    throw std::runtime_error("veil: unknown column type in scan.");
}

} // namespace veil
