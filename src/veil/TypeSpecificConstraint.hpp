// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <variant>

namespace veil
{

struct BoolConstraint final
{
    const bool min;
    const bool max;
    BoolConstraint(bool min, bool max) 
        : min(min), max(max) { }
};
struct DoubleConstraint final
{
    const double min;
    const double max;
    const bool allIntegral;
    DoubleConstraint(double min, double max, bool allIntegral) 
        : min(min), max(max), allIntegral(allIntegral) { }
};
struct IntConstraint final
{
    const int min;
    const int max;
    IntConstraint(int min, int max) 
        : min(min), max(max) { }
    bool includesZero() const { return this->min <= 0 && this->max >= 0; }
};
struct TextConstraint final
{
    const std::optional<std::string> constantText;
    const size_t minLength;
    const size_t maxLength;
    TextConstraint(std::optional<std::string> constantText_, size_t minLength, size_t maxLength) 
        : constantText(constantText_), minLength(minLength), maxLength(maxLength) { }
    bool hasValues() const { return this->minLength <= this->maxLength; }
    bool isConstant() const { return this->constantText.has_value(); }
    std::string constant() const { return this->constantText.value(); }
};

using TypeSpecificConstraint = std::variant<BoolConstraint, IntConstraint, DoubleConstraint, TextConstraint>;

}