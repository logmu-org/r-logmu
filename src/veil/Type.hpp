// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstdint>

namespace veil
{

enum class Type : uint8_t
{
    Bool = 1,
    Double = 2,
    Datey = 3,
    DateyInterval = 4,
    Durationy = 5,
    Text = 6,
    Category = 7, // Needs a `max` parameter!
};

// A CATEGORY HAS NO MISSING VALUE, and there is no sentinel here for one.
//
// A category value is an index into the crossing's StringMapping, so its values are 0 to max - 1 and
// every string has one. Its storage is `int`, and `int` HAS NO PLATFORM-INDEPENDENT NA
// REPRESENTATION: R's `NA_INTEGER` is `INT_MIN`, a convention no other platform repeats, and veil
// inventing a sentinel of its own would be no more portable for being ours. Nor does a nullable
// wrapper count -- the test is whether the primitive itself carries a portable missing value.
//
// `double` does (NaN) and `datey`/`durationy` do (below their representable range), which is why
// those three types alone have a missing state. The rule is representational, NOT a judgement about
// which kinds of fact can be unrecorded.
//
// So a factor column carrying an NA is REFUSED WHEN IT IS READ, beside the identical refusal for a
// logical. A dataset that needs "unknown" says so with an explicit level, which is better data
// anyway: it appears in a breakdown and cannot silently vanish from a total.

}