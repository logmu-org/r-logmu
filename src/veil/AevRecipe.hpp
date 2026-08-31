// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include "veil/Node.hpp"
#include "veil/Op.hpp"
#include "veil/Tree.hpp"

namespace veil
{

// The A/E/V recipe: three roots assembled from a mortality and a weight.
//
//     A = died_value(w)        the weight at the moment of death, nothing for a survivor
//     E = integrate(mu * w)    the expected deaths over the exposure
//     V = integrate(mu * w^2)  the second moment, which is what an A/E confidence interval needs
//
// A RECIPE IS AN ASSEMBLY, NOT A PASS. It builds nodes and hands back roots; every rewrite that
// follows is a general pass that knows nothing about AEV. That division is deliberate -- the
// log-likelihood and its first and second differentials are more recipes over the same machinery,
// with one root per output, and none of them should need the engine changed to accommodate it.
//
// MU IS ONE NODE, SHARED BY E AND V BY CONSTRUCTION. Both integrands point at the same NodeId, so
// lowering's memo computes the exponential once; that needs no help from the sharing pass, which is
// there for the sharing nobody arranged.
//
// THE SPELLING OF V MATTERS AND IS TIM'S (2026-07-28). It is `mu * (w * w)` -- the plain form -- and
// deliberately NOT `(mu * w) * w`, which would reuse E's product. Reusing it looks like a saving and
// is a trap: it forces `w` into the time vector, where it usually has no business being. A weight is
// normally constant, or a constant times a simple function of time, and the rewrite that pays is
// hoisting the time-invariant factor out of the integral altogether -- E = w * integral(mu),
// V = w^2 * integral(mu) -- which turns two vector multiplies and two reductions into one reduction
// and a few scalar multiplies. That is a general algebraic rewrite over `integrate`, and writing V to
// share E's product here would put it out of reach.
//
// V = E FALLS OUT, RATHER THAN BEING SPECIAL-CASED. Where `w` is provably zero or one, `w * w`
// simplifies to `w`, V's integrand becomes the same node as E's, and the sharing pass makes V free.
// Note an interval of [0, 1] is NOT enough to justify that -- 0.5 squared is not 0.5 -- so it needs
// the "only zero or one" fact the column scan and the always-bool rule produce.

struct AevRoots final
{
    NodeId a = invalidNodeId;
    NodeId e = invalidNodeId;
    NodeId v = invalidNodeId;
};

// Which way round the user wrote the second weighting factor. They are the same quantity -- a
// similarity s and a distance d are related by d = -log s -- but NEITHER IS CONVERTED INTO THE OTHER
// BEFORE IT GETS HERE, and that is deliberate. Wrapping one as the other in R would leave the other
// spelling's users needing a symbolic cancellation of log and exp to recover what they wrote, which
// no pass performs; and a log/exp round trip is not exact in IEEE 754, so a user writing
// `similarity = 0.5` would find their own number perturbed in the last bits for no purpose.
enum class SimilarityForm : unsigned char
{
    Similarity,  // the multiplier itself, used as written
    Distance,   // -log of the multiplier, so the multiplier is exp(-d)
};

// SIMILARITY MULTIPLIES THE LOG-LIKELIHOOD LINEARLY AND IS NEVER SQUARED. Where the weight w is
// linear in E and squared in V, the similarity s is linear in both:
//
//     A = died_value(s * w)
//     E = integrate(mu * (s * w))
//     V = integrate(mu * (s * (w * w)))
//
// S IS A SEPARATE FACTOR AND MUST NEVER BE FOLDED INTO THE WEIGHT. V collapses into E for free where
// w is provably zero or one, because `w * w` simplifies to `w` and the two integrands become the
// same node. Multiply s into w first and `s * w` is no longer an indicator for a fractional s, so
// that collapse is lost -- and with it the exactness of Z = 1 for an indicator weight. Spelled as
// above it survives untouched, since the simplification still sees a bare `w * w`.
//
// A DISTANCE IS BUILT AS ITS OWN `exp(-d)` NODE rather than fused into the mortality as
// `exp(log_mu - d)`. Fusing saves an exponential per time step, but only where d varies with time:
// where it does not -- the common case -- a separate factor hoists out of the integral altogether,
// and that hoist is worth far more than the exponential. Keeping d in the tree leaves the fusion
// available to a later pass, which is the other half of why nothing is converted in R.
//
// `logMu` is the mortality's tree, already ingested; `weight` likewise, and `similarity` may be
// `invalidNodeId` for none. Children are built before parents, so the arena still satisfies the
// ordering the type annotation pass relies on.
inline AevRoots buildAevRecipe(
    Tree& tree,
    NodeId logMu,
    NodeId weight,
    NodeId similarity = invalidNodeId,
    SimilarityForm form = SimilarityForm::Similarity)
{
    const NodeId mu = tree.buildCall(Op::Exp, {logMu});

    const NodeId factor = similarity == invalidNodeId || form == SimilarityForm::Similarity
        ? similarity
        : tree.buildCall(Op::Exp, {tree.buildCall(Op::Neg, {similarity})});

    const auto weighted = [&tree, factor](NodeId term)
    {
        return factor == invalidNodeId ? term : tree.buildCall(Op::Mul, {factor, term});
    };

    const NodeId weightTerm = weighted(weight);
    const NodeId a = tree.buildCall(Op::DiedValue, {weightTerm});
    const NodeId e = tree.buildCall(Op::Integrate, {tree.buildCall(Op::Mul, {mu, weightTerm})});

    const NodeId weightSquared = tree.buildCall(Op::Mul, {weight, weight});
    const NodeId v =
        tree.buildCall(Op::Integrate, {tree.buildCall(Op::Mul, {mu, weighted(weightSquared)})});

    return AevRoots{a, e, v};
}

} // namespace veil
