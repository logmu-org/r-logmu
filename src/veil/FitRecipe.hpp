// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>
#include "veil/AevRecipe.hpp"
#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"

namespace veil
{

// The fit recipe: everything ONE Newton-Raphson iteration needs, accumulated in one walk of the
// data. It is the same machinery as the A/E recipe with more roots, exactly as AevRecipe.hpp
// anticipated -- an assembly that builds nodes and hands back roots, with no pass and no engine
// change behind it.
//
// The model is proportional hazards on a base table, and log mu is LINEAR in beta:
//
//     log mu(i,t) = log mu_base(i,t) + sum_j beta_j X_j(i,t)
//
// so with w the weight, and writing Aw for a sum over deaths and Ew for an integral over exposure:
//
//     L*   = Omega^-1 (Aw log mu - Ew)
//     L*'  = Omega^-1 (AwX - EwX)
//     -L*'' = I = Omega^-1 Ew XX^T          and       J = Omega^-1 Ew^2 XX^T
//
// which is six families of root:
//
//     a                = died_value(w log mu)              1
//     e                = integrate(mu w)                   1
//     scoreActual[j]   = died_value(w X_j)                 k
//     scoreExpected[j] = integrate(mu w X_j)               k
//     ewXX[j][l]       = integrate(mu w X_j X_l)           k(k+1)/2
//     ew2XX[j][l]      = integrate(mu w^2 X_j X_l)         k(k+1)/2
//
// BETA IS NOT HERE, AND THAT IS DELIBERATE. `logMu` arrives already carrying the current beta, the
// way `buildAevRecipe` takes a mortality it does not question. Whether the fitter varies beta by
// refolding literals or by some later mechanism is the loop's business, not the recipe's; keeping it
// out is what lets this stay an assembly.
//
// OMEGA IS NOT APPLIED HERE EITHER. These roots are the raw moments `Ew XX^T` and `Ew^2 XX^T`, named
// so, and `I` and `J` are `Omega^-1` times them. The scalar is left to the caller because it belongs
// in exactly two places -- it cancels in the penalty `tr(J I^-1)` and appears once in
// `Var(beta_hat)` -- and a recipe that had already applied it would make the second one easy to
// apply twice.
//
// WHY THE INTEGRANDS ARE SPELLED AS A CHAIN OF PLAIN BINARY PRODUCTS.
// `passHoistFromIntegrate` peels time-invariant factors off such a chain, testing BOTH sides of each
// `Mul` and descending into whichever side still varies. What it needs is that each invariant factor
// is a separate operand somewhere on that chain; it does not care about the association, which was
// checked by rewriting these products with `mu` outermost and finding every test unmoved. What would
// defeat it is burying an invariant factor inside a subtree that also carries `mu`.
//
// The payoff, for the common model of a time-invariant weight and covariates that do not vary over
// an individual's exposure: each of the k(k+1)/2 information integrals peels down to
// `(w X_j X_l) * integrate(mu)`, and that one `integrate(mu)` is shared by all of them. MEASURED --
// going from one term to three leaves the vector work at 202 slot evaluations and 18 instructions,
// where switching the hoist off takes the same model to 662 and 62.
//
// IT NEEDS `mu` TO VARY WITH TIME. Against a constant mortality the whole integrand is invariant and
// the pass declines it, because lifting the last factor out would need the exposure length, which is
// a per-individual value rather than a factor. That is not a defect and nothing here should try to
// work around it.
//
// X_j AND X_l ARE MULTIPLIED TOGETHER FIRST, adjacently, for the same reason V is written
// `mu * (w * w)`: where a term is an indicator, `X_j * X_j` on the diagonal simplifies to `X_j` and
// the fold that does it is looking for a bare product of a node with itself. Splitting them across
// the chain would hide it. The product node is also shared by construction between `ewXX` and
// `ew2XX`, which differ only in the weight factor.
//
// CONTRIBUTIONS MUST BE OFF FOR A FIT. The diagnostic per-individual array in `CalculationResult` is
// `records * outputCount` doubles, and `outputCount` here is `2 + 2k + k(k+1)`. At twenty terms and
// a million records that is 462 outputs and about 3.7 GB.

// A packed symmetric matrix, upper triangle, ROW-MAJOR: the order is (0,0), (0,1) ... (0,k-1),
// (1,1) ... (1,k-1), (2,2) ... and so on. One definition, used by the recipe that fills it and by
// everything that reads it back, so the two cannot drift apart.
inline size_t packedTriangleSize(size_t terms) noexcept
{
    return terms * (terms + 1) / 2;
}

// The slot holding entry (row, column) of that packed triangle. Symmetric, so the two are swapped
// into order rather than refused.
inline size_t packedTriangleIndex(size_t row, size_t column, size_t terms)
{
    if (row > column) { const size_t swap = row; row = column; column = swap; }
    if (column >= terms) { throw std::out_of_range("veil: packed triangle index past the end."); }

    // Rows above this one occupy `terms - r` slots each, for r = 0 .. row-1.
    return row * terms - row * (row - 1) / 2 + (column - row);
}

struct FitRoots final
{
    NodeId a = invalidNodeId; // died_value(w log mu)
    NodeId e = invalidNodeId; // integrate(mu w)

    std::vector<NodeId> scoreActual;   // k of them, died_value(w X_j)
    std::vector<NodeId> scoreExpected; // k of them, integrate(mu w X_j)

    // Both packed upper triangles, k(k+1)/2 each, indexed by packedTriangleIndex. An entry left
    // `invalidNodeId` is a KNOWN ZERO that no output carries -- see `offDiagonalsOmitted`.
    std::vector<NodeId> ewXX;  // integrate(mu w X_j X_l)
    std::vector<NodeId> ew2XX; // integrate(mu w^2 X_j X_l)

    // Set when the covariates were VERIFIED mutually exclusive, in which case every off-diagonal
    // integrand `X_j X_l` is identically zero over the included exposure and the integral is not
    // built at all. One flag rather than a mask, because the pattern is entirely determined by it,
    // and `fitRootOrder` and `unpackFit` both derive the same pattern from this one place rather
    // than each holding a copy that has to agree.
    bool offDiagonalsOmitted = false;
};

// How many outputs a fit block has, given the terms and whether the off-diagonals were omitted. The
// one definition, so the builder's check and the loop's unpacking cannot drift.
inline size_t fitOutputCount(size_t terms, bool offDiagonalsOmitted) noexcept
{
    const size_t triangle = offDiagonalsOmitted ? terms : packedTriangleSize(terms);
    return 2 + 2 * terms + 2 * triangle;
}

// THE RECIPE BUILDS THE LINEAR PREDICTOR ITSELF, and that is a correctness property rather than a
// convenience. An earlier draft took a `logMu` that already carried beta, alongside the terms. That
// left nothing forcing the X inside the mortality to be the X in the score: hand it a mismatched
// pair and it computes a gradient and a Hessian for one model while evaluating the likelihood of
// another, silently, converging to the maximum of nothing. Taking `logMuBase` and the terms and
// forming
//
//     log mu = log mu_base + sum_j coefficient_j * X_j
//
// here means the same coerced node provably feeds mu, the score and the information, and the
// invariant is the recipe's to keep rather than the caller's to remember.
//
// `coefficients` must be the same length as `terms`. They are ordinary nodes: this recipe does not
// care whether a coefficient is a literal, as it is when accumulating at a fixed beta, or something
// the fitter can set between iterations. That choice belongs to whoever compiles the block.
//
// `terms` may be empty, which degenerates to the log-likelihood alone -- the case that can be
// checked against `aev()`. `weight` and the second weighting factor follow `buildAevRecipe` exactly,
// including that a distance is built as its own `exp(-d)` node rather than fused into the mortality.
//
// Children are built before parents, so the arena still satisfies the ordering the type annotation
// pass relies on.
inline FitRoots buildFitRecipe(
    Tree& tree,
    NodeId logMuBase,
    const std::vector<NodeId>& terms,
    const std::vector<NodeId>& coefficients,
    NodeId weight,
    NodeId similarity = invalidNodeId,
    SimilarityForm form = SimilarityForm::Similarity,

    // ONLY EVER SET FROM A VERIFIED ASSERTION, never from one the user merely made. Omitting an
    // integral is a compile-time decision, so a false claim here does not produce a slow answer or a
    // noisy one -- it produces a Hessian for a model nobody wrote, and the fit converges to the
    // maximum of nothing. `buildDisjointCheckRecipe` is what earns the right to pass true.
    bool disjointCovariates = false)
{
    if (coefficients.size() != terms.size())
    {
        throw std::invalid_argument("veil: a fit needs one coefficient for each model term.");
    }

    // EVERY TERM IS COERCED TO A PLAIN NUMBER, ONCE, AND THAT IS NOT A CONVENIENCE.
    //
    // `.x` is an age, which is a `durationy`, and it is the most common covariate there is -- a
    // Gompertz model is `log mu = a + b * age`. But the information wants `X_j X_l`, and
    // `durationy * durationy` is a product datey does not define and veil rightly refuses: there is
    // no such type, and inventing one to serve this recipe would put a unit into the engine that the
    // rest of the package does not have. A model term is a COVARIATE, and a covariate is a number;
    // whatever units it carries live in beta, which the engine never sees.
    //
    // Coercing ONCE, here, is what keeps the score and the information consistent. The same node
    // feeds `died_value(w X_j)` and `integrate(mu w X_j X_l)`, so there is no way for beta to mean a
    // duration in one and a number in the other.
    //
    // It is free where it matters. `ToDouble` is source-aware -- a click-backed value reads as years
    // and a logical as zero or one -- and `passLowerToBlock` emits nothing at all for a `ToDouble`
    // applied to something already vectorised, because a time vector holds doubles regardless. So a
    // time-varying term pays no instruction, and a logical term such as `.i$male` becomes the proper
    // 0/1 indicator it was always meant to be rather than a bool that `Mul` would refuse.
    std::vector<NodeId> covariates;
    covariates.reserve(terms.size());
    for (const NodeId term : terms)
    {
        covariates.push_back(tree.buildCall(Op::ToDouble, {term}));
    }

    // log mu = log mu_base + (sum_j coefficient_j * X_j), and THE BRACKETS ARE THE POINT.
    //
    // The contributions are summed among THEMSELVES first and added to the base ONCE. Folding them
    // into the base one at a time -- `((base + c0) + c1) + c2` -- makes every Add touch the
    // time-varying accumulator, so each one is a vector operation and the vector work grows with the
    // number of terms. Summed separately, a model whose covariates do not vary over an individual's
    // exposure has a wholly scalar sum and pays exactly ONE vector add however many terms it has.
    // Measured on three time-invariant covariates against a Gompertz base: 188 slot evaluations
    // left-folded, 132 this way, against 132 for a single term.
    //
    // Where the covariates do vary with time the sum is time-varying too and the two spellings cost
    // the same, so this is never worse.
    //
    // With no terms there is no sum and no add: log mu is the base, and the recipe reduces to a
    // log-likelihood.
    NodeId logMu = logMuBase;
    if (!covariates.empty())
    {
        NodeId contributions = tree.buildCall(Op::Mul, {coefficients[0], covariates[0]});
        for (size_t index = 1; index < covariates.size(); ++index)
        {
            const NodeId term = tree.buildCall(Op::Mul, {coefficients[index], covariates[index]});
            contributions = tree.buildCall(Op::Add, {contributions, term});
        }
        logMu = tree.buildCall(Op::Add, {logMuBase, contributions});
    }

    const NodeId mu = tree.buildCall(Op::Exp, {logMu});

    const NodeId factor = similarity == invalidNodeId || form == SimilarityForm::Similarity
        ? similarity
        : tree.buildCall(Op::Exp, {tree.buildCall(Op::Neg, {similarity})});

    const auto weighted = [&tree, factor](NodeId term)
    {
        return factor == invalidNodeId ? term : tree.buildCall(Op::Mul, {factor, term});
    };

    // Linear in the weight for the likelihood, the score and the information; the second weighting
    // factor is linear in both, and never squared. Spelled `w * w` plainly so that an indicator
    // weight still collapses the two triangles onto each other.
    const NodeId weightTerm = weighted(weight);
    const NodeId weightSquaredTerm = weighted(tree.buildCall(Op::Mul, {weight, weight}));

    // `mu` innermost, so the hoist can peel outwards through the chain.
    const auto integrateWeighted = [&tree, mu](NodeId weightFactor, NodeId covariate)
    {
        const NodeId integrand = covariate == invalidNodeId
            ? tree.buildCall(Op::Mul, {weightFactor, mu})
            : tree.buildCall(Op::Mul, {weightFactor, tree.buildCall(Op::Mul, {covariate, mu})});
        return tree.buildCall(Op::Integrate, {integrand});
    };

    FitRoots roots;
    roots.a = tree.buildCall(Op::DiedValue, {tree.buildCall(Op::Mul, {weightTerm, logMu})});
    roots.e = integrateWeighted(weightTerm, invalidNodeId);

    const size_t count = covariates.size();
    roots.scoreActual.reserve(count);
    roots.scoreExpected.reserve(count);
    for (const NodeId covariate : covariates)
    {
        roots.scoreActual.push_back(
            tree.buildCall(Op::DiedValue, {tree.buildCall(Op::Mul, {weightTerm, covariate})}));
        roots.scoreExpected.push_back(integrateWeighted(weightTerm, covariate));
    }

    roots.offDiagonalsOmitted = disjointCovariates;
    roots.ewXX.resize(packedTriangleSize(count), invalidNodeId);
    roots.ew2XX.resize(packedTriangleSize(count), invalidNodeId);
    for (size_t row = 0; row < count; ++row)
    {
        for (size_t column = row; column < count; ++column)
        {
            // THE WHOLE POINT OF THE ASSERTION. `X_j X_l` is identically zero off the diagonal, so
            // the integral is too, and it is left unbuilt rather than built and discarded -- the
            // saving is the grid walk, which is where the cost of a fit actually is.
            if (disjointCovariates && row != column) { continue; }

            // Built once and shared by the two triangles, which differ only in their weight factor.
            const NodeId product = tree.buildCall(Op::Mul, {covariates[row], covariates[column]});
            const size_t slot = packedTriangleIndex(row, column, count);
            roots.ewXX[slot] = integrateWeighted(weightTerm, product);
            roots.ew2XX[slot] = integrateWeighted(weightSquaredTerm, product);
        }
    }

    return roots;
}

// The roots in the order the block's outputs will hold them, which is the order the binding unpacks.
// One definition rather than two loops that have to agree.
inline std::vector<NodeId> fitRootOrder(const FitRoots& roots)
{
    std::vector<NodeId> order;
    order.reserve(2 + roots.scoreActual.size() + roots.scoreExpected.size()
                  + roots.ewXX.size() + roots.ew2XX.size());

    order.push_back(roots.a);
    order.push_back(roots.e);
    for (const NodeId id : roots.scoreActual) { order.push_back(id); }
    for (const NodeId id : roots.scoreExpected) { order.push_back(id); }

    // An omitted off-diagonal has no root and therefore no output. The survivors keep their packed
    // order, so with the off-diagonals gone they are simply the k diagonals in term order, which is
    // exactly what `unpackFit` puts back.
    for (const NodeId id : roots.ewXX) { if (id != invalidNodeId) { order.push_back(id); } }
    for (const NodeId id : roots.ew2XX) { if (id != invalidNodeId) { order.push_back(id); } }

    return order;
}

} // namespace veil
