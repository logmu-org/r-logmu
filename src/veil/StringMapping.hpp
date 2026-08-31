// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace veil
{

// The strings of one crossing, addressed by index.
//
// Text values are held as int32 indices, never as strings: the only text operation is comparison for
// equality, so equality of indices is all the hot path needs, and that is an integer test rather
// than a string one. Ordering is not offered at all -- an index records WHICH string, not where it
// sits in a sequence -- and is refused in passAnnotateTypes.
//
// THE CORE NEVER LEARNS WHERE THE INDICES CAME FROM, which is the point of the type. The R binding
// merges the level tables of its factor columns into one of these and remaps each column's codes as
// it copies them, so what crosses is integers plus this table. Another front end supplies whatever
// numbering its own platform already has. Nothing here knows about R.
//
// THE WORK IS SPLIT BY ORDER OF COST, which is why this holds levels rather than rows. Turning a
// million strings into codes is O(rows) and R already does it, in C, as `factor()`; a second interner
// here would duplicate that in the wrong place. Merging a handful of level tables is O(levels), costs
// nothing wherever it runs, and has to happen somewhere the front ends can share.
//
// ONE MAPPING PER CROSSING, NOT PER COLUMN. A batch may span several datasets and several columns,
// and a shared numbering is what lets a literal mean the same thing in all of them -- so comparing
// two text columns needs no reconciliation, and no index can be read against the wrong table.
//
// It is filled once at the boundary and const thereafter. That is what keeps it safe to read from a
// worker: nothing mutates it while a calculation runs, and no R memory is touched during one.
//
// Cost is bounded by CARDINALITY, not row count: a million-row sex column contributes two strings.
class StringMapping final
{
public:
    // EVERY STRING HAS AN INDEX, AND EVERY INDEX A STRING. There is no sentinel here for either
    // absence, because neither absence reaches this table: a literal naming a string no column holds
    // is settled while compiling, where `indexOfString` gives nothing back and `passEncodeText` folds
    // the comparison to a constant; and a record holding no value at all is refused when its column
    // is read, since a category is `int`-backed and `int` carries no portable NA.
    StringMapping() = default;

    // The index of one string, adding it at the next index if it is not already here. The boundary
    // calls this while merging the front end's per-column tables into the one shared numbering.
    //
    // FIND-OR-ADD RATHER THAN A BARE APPEND, because the same string reaches it once per column that
    // has that level -- "male" in a sex column and again in a spouse-sex column -- and appending twice
    // would leave two indices for one string, which would then compare unequal. Idempotent by
    // construction, so no caller has to remember to look first.
    int32_t addString(std::string text)
    {
        if (const std::optional<int32_t> found = this->indexOfString(text)) { return *found; }

        const auto index = static_cast<int32_t>(this->strings.size());
        this->strings.push_back(text);
        this->indices.emplace(std::move(text), index);
        return index;
    }

    // Absent when no row in the crossing holds this string, so a comparison against it folds to a
    // constant false rather than failing. A typo and a level absent from one dataset both land here,
    // and both are ordinary rather than exceptional.
    std::optional<int32_t> indexOfString(std::string_view text) const
    {
        const auto found = this->indices.find(std::string(text));
        if (found == this->indices.end()) { return std::nullopt; }
        return found->second;
    }

    // Renders an index back to its string, for diagnostics. This is the only direction the core
    // uses, and it uses it only when reporting -- never while calculating.
    //
    // Checked rather than trusted, though every index it is asked for was handed out by this same
    // table: this runs when something is being explained to a user, so it is nowhere near hot enough
    // for a bounds test to matter, and reading off the end would be undefined behaviour.
    const std::string& stringAtIndex(int32_t index) const
    {
        if (index < 0 || static_cast<size_t>(index) >= this->strings.size())
        {
            throw std::out_of_range("veil: no string has this index.");
        }
        return this->strings[static_cast<size_t>(index)];
    }

    size_t size() const noexcept { return this->strings.size(); }

private:
    std::vector<std::string> strings;
    std::unordered_map<std::string, int32_t> indices;
};

} // namespace veil
