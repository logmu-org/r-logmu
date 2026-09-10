// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>
#include "veil/Block.hpp"
#include "veil/Cholesky.hpp"
#include "veil/ColumnView.hpp"
#include "veil/Engine.hpp"
#include "veil/FitRecipe.hpp"

namespace veil
{

// The Newton-Raphson loop: maximum weighted log-likelihood for a proportional hazards model.
//
// THE BLOCK IS COMPILED ONCE AND NEVER TOUCHED AGAIN. Between iterations the only thing that changes
// is k doubles, written with `Block::setParameter`. Recompiling would be wrong rather than merely
// slow: folding would bake each beta into constants afresh, so the block could change SHAPE from one
// iteration to the next -- different sharing, different hoisting -- and the objective the line
// search evaluates would not be structurally the one the Newton step was computed from.
//
// ONE WALK PER EVALUATION, AND ONE EVALUATION PER ITERATION when nothing is damped. The accepted
// trial's accumulation carries into the next iteration rather than being recomputed, which is why
// there is no second likelihood-only block: a trial evaluation and the next iteration's
// accumulation are the same walk.
//
// OMEGA AND Z NEVER ENTER THE LOOP. The Newton step is free of both -- neither moves the argmax --
// and they scale both sides of the acceptance test identically, so they cancel out of it. They
// appear only in the convergence test and in what is reported. Everything inside runs on the RAW
// accumulations, `A = Ew XX^T` and `B = Ew^2 XX^T`, with `I = Omega^-1 A` and `J = Omega^-1 B`.
//
// THE SCALAR ARITHMETIC IS O(k^3) AND THE WALK IS O(records * slots), so nothing here is worth
// optimising against the walk. At twenty parameters and a million records the factorisation is about
// one part in four million of an iteration.

// How a fit ended. **A FIT RETURNS A FIT OR IT FAILS** (Tim, 2026-09-07): not converged is also not
// a fit, and a warning is too easy to ignore. Separate codes because the message is the whole of
// what a user gets, and "your covariates are collinear" and "your data contains a NaN" send someone
// to entirely different places.
enum class FitStatus : unsigned char
{
    Converged,
    NotIdentifiable,   // the information matrix would not factorise: collinear, or no exposure
    NotFinite,         // a NaN or infinity reached the gradient or the likelihood
    DidNotConverge,    // the iteration budget ran out with the gain still above the tolerance
    StepCollapsed,     // backtracking halved past the smallest useful step without improving
    Interrupted,       // the host asked to stop
};

struct FitControl final
{
    // Beta's starting point. Empty means every coefficient starts at zero, which is the default.
    std::vector<double> start;

    size_t maxIterations = 25;

    // In units of L, where 1 is about one parameter's worth. 1e-6 rather than 1e-4 because the loss
    // in L from being off by d is half its square in standard-error units, so a tolerance eps leaves
    // beta within sqrt(2 eps) standard errors -- the square root, not eps. 1e-6 is 0.0014 of a
    // standard error, and Newton squares the error each late iteration so tightening it costs about
    // one more.
    double convergenceTolerance = 1e-6;

    // The Armijo constant: the fraction of the linearly predicted gain a trial step must actually
    // capture. Its CEILING IS ONE HALF, because one half is exactly what a full Newton step delivers
    // on a quadratic -- any larger and the full step would be rejected near the optimum and the
    // iteration would stall just as it should converge. 1e-4 leaves four orders of margin, so the
    // condition is close to inert until a step has overshot into territory where the function has
    // turned over.
    double armijo = 1e-4;

    // Backtracking gives up here. 2^-30 is far below any step that could still be making progress.
    size_t maxHalvings = 30;

    double overdispersion = 1.0;

    // Ew^2 / Ew on a test mortality over the same data, computed ONCE and held. It puts L on the
    // scale where one parameter costs about 1, which is what makes the tolerance mean what it says.
    double zScale = 1.0;
};

struct FitResult final
{
    FitStatus status = FitStatus::DidNotConverge;

    std::vector<double> beta;
    size_t iterations = 0;
    size_t evaluations = 0; // Walks of the data, so damping is visible in the cost.

    // On the L scale: `L = (Aw log mu - Ew) / (Omega Z)`.
    double logLikelihood = std::numeric_limits<double>::quiet_NaN();

    // What the convergence test last saw: the gain in L still available, `lambda^2 / (2 Omega Z)`.
    double predictedGain = std::numeric_limits<double>::quiet_NaN();

    // Packed upper triangles, `packedTriangleIndex`. Filled only on success.
    std::vector<double> variance; // Omega * A^-1 B A^-1
    double penalty = std::numeric_limits<double>::quiet_NaN(); // p = tr(B A^-1) / Z

    // Set when the status is NotIdentifiable. `dependency` expresses the offending covariate in
    // terms of the ones before it, so a message can name the actual collinearity rather than a
    // pivot index -- see `choleskyDependency`.
    size_t failedParameter = 0;
    CholeskyStatus factorStatus = CholeskyStatus::Ok;
    std::vector<double> dependency;
};

namespace detail
{

// Every accumulated family of one walk, unpacked from the block's outputs in `fitRootOrder`'s order.
struct FitAccumulation final
{
    double actual = 0.0;      // Aw log mu
    double expected = 0.0;    // Ew
    std::vector<double> gradient;    // AwX - EwX, the raw score
    std::vector<double> information; // Ew XX^T, packed
    std::vector<double> secondMoment; // Ew^2 XX^T, packed

    // The objective the line search compares, on the raw scale: Omega and Z cancel out of the
    // acceptance test, so they are left out of it entirely.
    double rawLogLikelihood() const noexcept { return this->actual - this->expected; }
};

// Spreads the k diagonal outputs of a disjoint block back over a full packed triangle, filling the
// off-diagonals with exact zeros -- which is what they are, having been verified so before the block
// was compiled. Everything downstream then sees the shape it always saw, so the Cholesky, the
// variance and the penalty need to know nothing about any of this.
inline std::vector<double> spreadDiagonal(const double* diagonal, size_t terms)
{
    std::vector<double> packed(packedTriangleSize(terms), 0.0);
    for (size_t j = 0; j < terms; ++j) { packed[packedTriangleIndex(j, j, terms)] = diagonal[j]; }
    return packed;
}

inline FitAccumulation unpackFit(const std::vector<double>& totals, size_t terms,
                                 bool offDiagonalsOmitted)
{
    // THE LENGTH CHECK IS THE GUARD AGAINST A MISMATCHED FLAG. A block compiled with the
    // off-diagonals omitted and unpacked without the flag -- or the reverse -- would otherwise read
    // the second moment as the information and be silently wrong, so the sizes are made to disagree
    // rather than left to line up by luck.
    const size_t triangle = packedTriangleSize(terms);
    const size_t carried = offDiagonalsOmitted ? terms : triangle;
    if (totals.size() != fitOutputCount(terms, offDiagonalsOmitted))
    {
        throw std::runtime_error("veil: a fit block produced the wrong number of outputs.");
    }

    FitAccumulation out;
    out.actual = totals[0];
    out.expected = totals[1];

    out.gradient.reserve(terms);
    for (size_t j = 0; j < terms; ++j)
    {
        out.gradient.push_back(totals[2 + j] - totals[2 + terms + j]);
    }

    const size_t informationAt = 2 + 2 * terms;
    if (offDiagonalsOmitted)
    {
        out.information = spreadDiagonal(totals.data() + informationAt, terms);
        out.secondMoment = spreadDiagonal(totals.data() + informationAt + carried, terms);
        return out;
    }

    out.information.assign(totals.begin() + static_cast<std::ptrdiff_t>(informationAt),
                           totals.begin() + static_cast<std::ptrdiff_t>(informationAt + carried));
    out.secondMoment.assign(totals.begin() + static_cast<std::ptrdiff_t>(informationAt + carried),
                            totals.begin() + static_cast<std::ptrdiff_t>(informationAt + 2 * carried));
    return out;
}

inline bool allFinite(const std::vector<double>& values)
{
    for (const double value : values)
    {
        if (!std::isfinite(value)) { return false; }
    }
    return true;
}

} // namespace detail

// Runs the loop. `block` must be a fit block with `terms` coefficients lowered as parameters, and it
// is taken by non-const reference for exactly one reason: the coefficients are set on it.
inline FitResult runFit(
    Block& block,
    const std::vector<const ColumnView*>& columns,
    size_t records,
    size_t terms,

    // How the BLOCK was built, not a choice the loop makes: true when the covariates were verified
    // mutually exclusive and the block therefore carries only the diagonal of each triangle.
    bool offDiagonalsOmitted,
    const FitControl& control,
    size_t threads,
    const std::function<bool()>& interrupted)
{
    if (!control.start.empty() && control.start.size() != terms)
    {
        throw std::runtime_error("veil: the starting beta has the wrong number of coefficients.");
    }
    if (!(control.overdispersion > 0.0) || !(control.zScale > 0.0))
    {
        throw std::runtime_error("veil: overdispersion and Z must both be positive.");
    }

    FitResult result;
    result.beta = control.start.empty() ? std::vector<double>(terms, 0.0) : control.start;

    // Both scalars enter here and nowhere else in the loop.
    const double scale = control.overdispersion * control.zScale;

    // Set the coefficients, walk the data, unpack. `stopped` is a flag rather than an exception
    // because an interruption is something the host ASKED for, not an error -- the same reason
    // `CalculationResult` reports it as a field.
    bool stopped = false;
    const auto walk = [&](const std::vector<double>& beta) -> detail::FitAccumulation
    {
        for (size_t j = 0; j < terms; ++j)
        {
            block.setParameter(static_cast<ParamId>(j), beta[j]);
        }
        const CalculationResult calculation =
            runCalculation(block, columns, records, false, threads, interrupted);
        if (calculation.interrupted)
        {
            // THE TOTALS ARE MEANINGLESS when a run is cut short -- some chunks contributed and some
            // did not -- so nothing is unpacked from them.
            stopped = true;
            return detail::FitAccumulation{};
        }
        ++result.evaluations;
        return detail::unpackFit(calculation.totals, terms, offDiagonalsOmitted);
    };

    detail::FitAccumulation current = walk(result.beta);
    if (stopped)
    {
        result.status = FitStatus::Interrupted;
        return result;
    }

    for (result.iterations = 1; result.iterations <= control.maxIterations; ++result.iterations)
    {
        result.logLikelihood = current.rawLogLikelihood() / scale;

        // ASKED BEFORE THE FACTORISATION, because a NaN pivot would otherwise be reported as an
        // ordinary singularity and send the user hunting for a collinear covariate that does not
        // exist. A NaN here is bad data, not a bad model.
        if (!std::isfinite(current.rawLogLikelihood())
            || !detail::allFinite(current.gradient)
            || !detail::allFinite(current.information))
        {
            result.status = FitStatus::NotFinite;
            return result;
        }

        const CholeskyFactor factor = choleskyFactorPacked(current.information, terms);
        if (!factor.positiveDefinite())
        {
            result.status = FitStatus::NotIdentifiable;
            result.factorStatus = factor.status;
            result.failedParameter = factor.failedParameter;
            result.dependency = choleskyDependency(factor);
            return result;
        }

        const std::vector<double> step = choleskySolve(factor, current.gradient);

        // FROM THE SQUARED FORM, NOT FROM `gradient . step`. The quadratic form is the squared norm
        // of the solve's own forward substitution, so it is a sum of squares and CANNOT come back
        // negative; the dot product can, by rounding, once the gradient is near zero. The two agree
        // to rounding everywhere else.
        //
        // The half is the second-order term: along the Newton direction the linear part promises
        // `grad . step` and the curvature gives half of it back, so half is the gain actually
        // available. `predictedGain` is therefore a change in L, and the tolerance is compared
        // against it directly -- which is what makes the tolerance mean "the largest insignificant
        // change in L" rather than half of one.
        const double decrement = choleskyQuadraticForm(factor, current.gradient);
        result.predictedGain = decrement / (2.0 * scale);

        if (result.predictedGain <= control.convergenceTolerance)
        {
            result.status = FitStatus::Converged;
            result.variance = choleskySandwich(factor, current.secondMoment);
            for (double& entry : result.variance) { entry *= control.overdispersion; }
            result.penalty = choleskyTraceOfProductWithInverse(factor, current.secondMoment)
                / control.zScale;
            return result;
        }

        // THE FULL STEP FIRST, because it is what is accepted in all but the earliest iterations,
        // and the accepted trial IS the next iteration's accumulation.
        std::vector<double> trialBeta(terms, 0.0);
        detail::FitAccumulation trial;
        bool accepted = false;
        double length = 1.0;

        for (size_t halving = 0; halving <= control.maxHalvings; ++halving)
        {
            for (size_t j = 0; j < terms; ++j)
            {
                trialBeta[j] = result.beta[j] + length * step[j];
            }

            trial = walk(trialBeta);
            if (stopped)
            {
                result.status = FitStatus::Interrupted;
                return result;
            }

            // A NON-FINITE TRIAL IS NOT AN IMPROVEMENT. `Ew` overflows once `beta^T X` reaches
            // about 709, which a first step from a badly levelled reference table can easily clear,
            // and the objective then comes back as an infinity or a NaN.
            //
            // THE `isfinite` TEST IS BELT AND BRACES AND IS KNOWN TO BE SO -- a disable-and-recheck
            // removed it and not one test moved. It is redundant because every comparison against a
            // NaN is false, so a NaN gain fails the Armijo inequality on its own and the step is
            // halved regardless; and an infinite `Ew` gives a gain of minus infinity, which fails it
            // too. It is kept because the redundancy depends on how the inequality below happens to
            // be written, and a later hand rearranging that should not silently lose the guard.
            // Do not claim a test covers this: none does, and none easily could.
            //
            // ARMIJO RATHER THAN ANY IMPROVEMENT. A step may improve the objective by an arbitrarily
            // small amount and still be accepted, which is how an iteration crawls to a halt away
            // from the optimum. Omega and Z cancel from both sides, so this runs on raw values.
            const double gain = trial.rawLogLikelihood() - current.rawLogLikelihood();
            if (std::isfinite(trial.rawLogLikelihood())
                && gain >= control.armijo * length * decrement)
            {
                accepted = true;
                break;
            }
            length *= 0.5;
        }

        if (!accepted)
        {
            result.status = FitStatus::StepCollapsed;
            return result;
        }

        result.beta = trialBeta;
        current = trial;
    }

    // The budget ran out. `iterations` is one past the last, so bring it back.
    result.iterations = control.maxIterations;
    result.logLikelihood = current.rawLogLikelihood() / scale;
    result.status = FitStatus::DidNotConverge;
    return result;
}

} // namespace veil
