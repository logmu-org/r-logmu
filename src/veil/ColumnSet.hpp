// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include "veil/ColumnView.hpp"
#include "veil/TypeFull.hpp"

namespace veil
{

// The referenced columns of ONE dataset, ready to read.
//
// A batch may span several datasets, so a batch holds one of these per dataset it touches. What it
// does NOT hold is the StringMapping: string indices are scoped to the whole crossing, so that a
// literal means the same thing in every calculation and in every dataset.
//
// The set owns the columns it had to convert -- bools, which R holds as int32, and text, which is
// remapped onto the crossing's shared numbering -- and merely views the ones that needed no copy.
// Owned blocks are held by unique_ptr rather than inline, so a view keeps pointing at the same
// address when the set grows.
//
// LIFETIME: viewed columns point into memory the set does not own, so whatever supplied them (in R's
// case, the protected SEXPs) must outlive the set.
class ColumnSet final
{
public:
    // A column needing no copy: the block is somebody else's, and must outlive this set.
    ColumnId addView(std::string name, TypeFull type, const void* data, size_t count)
    {
        return this->add(std::move(name), ColumnView(type, data, count));
    }

    // R's int32 logicals, converted to standard 1-byte bools.
    ColumnId addBools(std::string name, std::unique_ptr<bool[]> values, size_t count)
    {
        const bool* data = values.get();
        this->ownedBools.push_back(std::move(values));
        return this->add(std::move(name), ColumnView(TypeFull::createBool(), data, count));
    }

    // A text column, encoded to indices into the crossing's StringMapping. Its type is Category
    // rather than Text: what the core holds IS a category, and the strings it was written as are the
    // mapping's business. `levels` is the mapping's size, which bounds every code the column can hold
    // and is what the scan tests a value against to tell a code from a missing one.
    ColumnId addCodes(std::string name, std::unique_ptr<int32_t[]> values, size_t count, int levels)
    {
        const int32_t* data = values.get();
        this->ownedCodes.push_back(std::move(values));
        return this->add(std::move(name), ColumnView(TypeFull::createCategory(levels), data, count));
    }

    std::optional<ColumnId> find(std::string_view name) const
    {
        for (size_t i = 0; i < this->names.size(); ++i)
        {
            if (this->names[i] == name) { return static_cast<ColumnId>(i); }
        }
        return std::nullopt;
    }

    const ColumnView& at(ColumnId id) const { return this->views[id]; }
    const std::string& name(ColumnId id) const { return this->names[id]; }

    size_t size() const noexcept { return this->views.size(); }

    // Every column of a dataset is the same length, so this is the individual count.
    size_t recordCount() const noexcept { return this->views.empty() ? 0 : this->views.front().count; }

private:
    ColumnId add(std::string name, ColumnView view)
    {
        if (!this->views.empty() && view.count != this->recordCount())
        {
            throw std::logic_error("veil: a dataset's columns must all be the same length.");
        }
        if (this->find(name).has_value())
        {
            throw std::logic_error("veil: a dataset cannot have two columns of the same name.");
        }

        this->names.push_back(std::move(name));
        this->views.push_back(view);
        return static_cast<ColumnId>(this->views.size() - 1);
    }

    std::vector<std::string> names;
    std::vector<ColumnView> views;
    std::vector<std::unique_ptr<bool[]>> ownedBools;
    std::vector<std::unique_ptr<int32_t[]>> ownedCodes;
};

} // namespace veil
