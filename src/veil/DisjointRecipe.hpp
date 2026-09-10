// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"

namespace veil
{

// The pre-flight walk that earns the right to omit a fit's off-diagonal integrals.
//
// WHY IT EXISTS AT ALL. `disjoint(...)` is an assertion the USER makes: no individual belongs to
// more than one term. R cannot check it -- a column's type is unknown until the data arrives, so it
// cannot even prove the terms are indicators -- and omitting an integral is a COMPILE-TIME decision,
// taken before a single record is read. That leaves three options and only one is safe: trust the
// claim and compute a Hessian for a model nobody wrote, ignore it unless provable and gain nothing,
// or verify it once against the data and then compile the reduced block. This is the third.
//
// WHAT IT COMPUTES. One root per unordered pair,
//
//     integrate(|X_j * X_l|)      for every j < l
//
// and the assertion holds exactly when every total is zero.
//
// WHY A SUM IS ENOUGH. Absolute values are non-negative, so a sum of them is zero if and only if
// every one is. That is what makes the check a handful of accumulators rather than a per-record
// interrogation: nothing has to be remembered, nothing has to short-circuit, and the record chunks
// still fold in chunk order like every other reduction.
//
// WHY THE ABSOLUTE VALUE, when indicators are non-negative anyway. Because non-negativity is
// PRECISELY THE PROPERTY NOBODY HAS PROVED -- it is the same unknown that makes this walk necessary.
// Without it a term that is -1 on one individual and +1 on another cancels in the total and a false
// assertion sails through into a silently wrong Hessian. `Abs` costs one scalar operation per pair
// and makes the check sound whatever the terms turn out to hold.
//
// WHY `Integrate` RATHER THAN A SUM OVER RECORDS. Two reasons, and the second is the real one. It
// needs no new reduction -- `Integrate` and `DiedValue` are the only two the engine has. And it
// checks the property over EXACTLY THE POPULATION THE FIT INTEGRATES: the same include, the same
// exposure clipping. An individual excluded, or clipped to no exposure at all, contributes zero to
// every off-diagonal integral whatever groups they belong to, so their membership genuinely does not
// matter -- and a record-level sum would refuse fits that are perfectly sound.
//
// WHAT IT COSTS. No mortality, no exponential and no beta: the products are time-invariant, so they
// are computed once per record and broadcast across the grid. Against one Newton iteration -- which
// exponentiates at every slot of every record and then integrates the whole triangle -- this is
// small, and it is paid once for a fit that walks the data many times.
//
// A REJECTED ALTERNATIVE, because it looks cheaper and is not sound. `integrate(max(0, sum_j X_j -
// 1))` is ONE root instead of k(k-1)/2, and for genuine 0/1 indicators "no individual sums above
// one" is the same statement as pairwise exclusivity. It leans on the terms being indicators, which
// is the assumption this walk exists to avoid making, and it cannot name the offending pair.

// The pairs, in the order `buildDisjointCheckRecipe` returns roots for them: `(0,1), (0,2), ...,
// (0,k-1), (1,2), ...`. One definition, so a non-zero total maps back to the terms a user wrote.
inline std::vector<std::pair<size_t, size_t>> disjointCheckPairs(size_t terms)
{
    std::vector<std::pair<size_t, size_t>> pairs;
    if (terms < 2) { return pairs; }

    pairs.reserve(terms * (terms - 1) / 2);
    for (size_t row = 0; row < terms; ++row)
    {
        for (size_t column = row + 1; column < terms; ++column)
        {
            pairs.emplace_back(row, column);
        }
    }
    return pairs;
}

// Roots for the check, in `disjointCheckPairs` order. Fewer than two terms has nothing to check and
// yields no roots, which the caller should read as "there is no block to run".
inline std::vector<NodeId> buildDisjointCheckRecipe(Tree& tree, const std::vector<NodeId>& terms)
{
    // COERCED THE SAME WAY THE FIT COERCES, and it must be. A logical term is what `Mul` would
    // otherwise refuse, and a `durationy` one squares to a product datey does not define -- the same
    // two reasons `buildFitRecipe` puts every term through `ToDouble` once. Checking a differently
    // coerced X from the one the fit integrates would be checking a different claim.
    std::vector<NodeId> covariates;
    covariates.reserve(terms.size());
    for (const NodeId term : terms)
    {
        covariates.push_back(tree.buildCall(Op::ToDouble, {term}));
    }

    std::vector<NodeId> roots;
    for (const auto& [row, column] : disjointCheckPairs(terms.size()))
    {
        const NodeId product = tree.buildCall(Op::Mul, {covariates[row], covariates[column]});
        const NodeId magnitude = tree.buildCall(Op::Abs, {product});
        roots.push_back(tree.buildCall(Op::Integrate, {magnitude}));
    }
    return roots;
}

// Which pair a set of totals convicts, or `pairs.size()` when they all came back zero.
//
// NaN MUST CONVICT, because an assertion that cannot be evaluated has not been verified. The trap is
// NOT `!=`, which is equally safe -- an inequality against a NaN is TRUE, so `total != 0.0` convicts
// one too. It is any ORDERED test: `total > 0.0` is false for a NaN and would wave it through, which
// is the same shape of mistake as the Cholesky pivot's `pivot <= 0.0`. The equality form is written
// negated so that the only comparison here is one NaN is guaranteed to fail.
inline size_t firstDisjointViolation(const std::vector<double>& totals, size_t pairs)
{
    if (totals.size() != pairs)
    {
        throw std::runtime_error("veil: a disjointness check produced the wrong number of outputs.");
    }
    for (size_t i = 0; i < pairs; ++i)
    {
        if (!(totals[i] == 0.0)) { return i; }
    }
    return pairs;
}

} // namespace veil
