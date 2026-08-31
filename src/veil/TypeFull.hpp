// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include "veil/Type.hpp"

namespace veil
{

class TypeFull final
{
public:
    const Type type;
    const int max;

private:
    constexpr TypeFull(Type type, int max = 0)
        : type(type), max(max) {}

public:
    static constexpr TypeFull createBool() { return TypeFull(Type::Bool); }
    static constexpr TypeFull createDouble() { return TypeFull(Type::Double); }
    static constexpr TypeFull createDatey() { return TypeFull(Type::Datey); }
    static constexpr TypeFull createDateyInterval() { return TypeFull(Type::DateyInterval); }
    static constexpr TypeFull createDurationy() { return TypeFull(Type::Durationy); }
    static constexpr TypeFull createText() { return TypeFull(Type::Text); }
    static constexpr TypeFull createCategory(int max) { return TypeFull(Type::Category, max); }

    constexpr bool operator==(const TypeFull& other) const = default;
};

}