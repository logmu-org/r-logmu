// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include "veil/Type.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// Identifies a column within ONE dataset. A batch may span several datasets, so a ColumnId is only
// meaningful against the ColumnSet of the dataset its calculation reads. String indices are different
// -- those are scoped to the whole crossing, and numbered by the StringMapping.
using ColumnId = uint32_t;

// A typed view of one column: a tag, a pointer to a block of `count` values, and nothing else.
//
// The pointer is `const void*` and stays private. Every read goes through a typed accessor which
// checks the tag first, so the cast happens in exactly one place and always back to the type the
// block genuinely holds. That is not type punning -- nothing reinterprets one representation as
// another -- and it aliases nothing.
//
// A `std::variant` of pointer types would be worse here rather than safer: Datey, Durationy and
// Category all map to `const int*`, so the variant index cannot encode TypeFull and both would have
// to be carried and kept in step. One tag, one invariant, one checkpoint.
//
// The block is either R's own memory (double and datey columns are read zero-copy; R has no
// compacting GC, so the address is stable while the SEXP is protected) or a buffer the ColumnSet
// owns (bools, which R holds as int32, and text, whose factor codes are remapped onto the crossing's
// shared numbering).
//
// LIFETIME: for a zero-copy column this views memory the binding does not own. The binding must
// keep the SEXP protected for as long as any view of it is live.
class ColumnView final
{
public:
    ColumnView(TypeFull type, const void* data, size_t count)
        : type(type), count(count), data(data) {}

    const TypeFull type;
    const size_t count;

    const double* doubles() const
    {
        this->requireType(Type::Double);
        return static_cast<const double*>(this->data);
    }

    // Datey, durationy and category are all 32-bit integer clicks or codes.
    const int* ints() const
    {
        if (this->type.type != Type::Datey && this->type.type != Type::Durationy
            && this->type.type != Type::Category)
        {
            throw std::logic_error("veil: column is not an integer-backed type.");
        }
        return static_cast<const int*>(this->data);
    }

    // R holds logicals as int32; the reader converts them so the core sees a standard 1-byte bool.
    const bool* bools() const
    {
        this->requireType(Type::Bool);
        return static_cast<const bool*>(this->data);
    }

private:
    void requireType(Type wanted) const
    {
        if (this->type.type != wanted)
        {
            throw std::logic_error("veil: column read at the wrong type.");
        }
    }

    const void* data;
};

} // namespace veil
