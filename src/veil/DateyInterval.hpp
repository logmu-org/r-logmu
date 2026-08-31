// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once 

#include <cstdint>
#include <bit>

namespace veil
{

struct DateyInterval
{
    const int Start;
    const int End;
    DateyInterval(int start, int end) : Start(start), End(end) { }
    static DateyInterval FromUInt64(uint64_t bits)
    {
        auto start = static_cast<int32_t>(bits >> 32);
        auto end = static_cast<int32_t>(bits);
        return DateyInterval(start, end);
    }
    static uint64_t ToUInt64(DateyInterval dateyInterval)
    {
        // Avoid sign extension
        // std::bit_cast available from C++20
        auto start = static_cast<uint64_t>(std::bit_cast<uint32_t>(dateyInterval.Start));
        auto end = static_cast<uint64_t>(std::bit_cast<uint32_t>(dateyInterval.End));
        return start << 32 | end;
    }
};

}
