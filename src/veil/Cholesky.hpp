// LogMu mortality experience analysis and model fitting
//
// This file is licensed to you under the Apache Licence 2.0.
//
// Copyright (c) Tim Gordon

#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>
#include "veil/FitRecipe.hpp" // packedTriangleSize, packedTriangleIndex

namespace veil
{

// The whole of the dense linear algebra a Newton-Raphson fit needs, and no more of it.
//
// One iteration wants, with k parameters:
//
//     I = Omega^-1 ewXX      (k x k, symmetric, positive definite when the model is identified)
//     J = Omega^-1 ew2XX     (k x k, symmetric)
//     grad                   (k-vector)
//
//     step               delta      = I^-1 grad
//     Newton decrement   lambda^2   = grad^T I^-1 grad = grad^T delta
//     penalty            p          proportional to tr(J I^-1)
//     variance           Var(beta)  = I^-1 J I^-1
//
// That is a Cholesky, a triangular solve, an inverse from the factor, a quadratic form and a trace.
// No eigenvalues, no SVD, no QR, nothing sparse and nothing iterative, so there is no library here
// and there is not going to be one. The reasons, written down so they do not have to be rediscovered:
//
//   - src/veil/ is R-free on purpose, which is what lets a Python or C# front end reuse it. R's
//     bundled LAPACK would end that; RcppEigen would drag in the Rcpp dependency the package has
//     stayed clear of.
//   - The CRAN check already carries one WARNING for vendored EVE's pragmas. A second vendored
//     library is a second thing a reviewer has to be talked through.
//   - The reproducibility standard is that the same binary on the same hardware agrees EXACTLY. An
//     unblocked hand-written factorisation is trivially deterministic. A tuned library's blocked or
//     threaded path is a promise nobody has made us.
//   - The failure path is the interesting part, and it is ours to shape -- see below.
//   - k is tens. Accumulating ewXX over the data is the expensive half and it is already built.
//
// PLAIN `double` THROUGHOUT, IEEE 754-2019. Not `long double`: on x86 that is the 80-bit x87 type,
// its width and its excess precision vary by platform and compiler, and that would break the
// reproducibility standard the moment someone built on a machine that spilled differently.
//
// THE FAILURE PATH IS THE POINT, NOT AN AFTERTHOUGHT. A singular or indefinite I means collinear
// covariates, or a level with no exposure behind it. A Cholesky knows exactly which pivot went
// non-positive, and because the factorisation runs parameter by parameter that pivot index IS a
// parameter index. So the failure is REPORTED, never thrown: `positiveDefinite` false and
// `failedParameter` set, and the fitter turns that into a message naming the covariate rather than a
// generic numerical failure.
//
// THE PIVOT IS TESTED AGAINST THE DIAGONAL ENTRY IT CAME FROM, and that ratio has a meaning rather
// than being an arbitrary epsilon. Factorising `A = EwXX^T`, the pivot at row j is what is left of
// `A[j][j]` after removing what the earlier covariates explain, so
//
//     pivot / A[j][j] = 1 - R^2
//
// the uncentred R-squared of covariate j on the covariates before it, weighted by `w mu`. The
// threshold therefore asks what fraction of a covariate has to be its own rather than a rehash of
// its predecessors. `NearlySingularTolerance` is 1e-7, matching the rank tolerance `lm` and
// `glm.fit` use in their QR.
//
// AN EARLIER VERSION TESTED EXACT POSITIVITY ONLY, and that was too weak. Two covariates agreeing to
// twelve digits give a small POSITIVE pivot, the factorisation succeeds, and `A^-1` comes back
// enormous -- a fit that looks like an answer and is not one. At 1e-7 a coefficient's standard error
// is already inflated about three thousandfold against the orthogonal case, so nothing that trips it
// is estimated in any useful sense. Nothing legitimate is near it either: age and age squared over
// ages 60 to 90, about as collinear as a sane model gets, sit at 1 - R^2 = 0.015.

// How a factorisation ended. Separate codes rather than one flag, because the fitter refuses on any
// of them and the message it writes is the whole of what the user gets -- and "this covariate is
// collinear with an earlier one" and "your data contains a NaN" send someone to entirely different
// places.
enum class CholeskyStatus : unsigned char
{
    Ok,
    NotFinite,           // a pivot was NaN or infinite, so the data carried one in
    NotPositiveDefinite, // a pivot was zero or negative: exactly singular, or not a covariance at all
    NearlySingular,      // a pivot survived positivity but not the ratio test below
};

// What fraction of a covariate must be independent of the covariates before it. See the note above
// for why this is `1 - R^2` and why 1e-7 is the number.
constexpr double NearlySingularTolerance = 1e-7;

// A factorisation of a symmetric positive-definite matrix A as
//
//     A = U^T U
//
// with U upper triangular and its diagonal positive. `factorPacked` holds U in the SAME layout as
// the input: a packed upper triangle, ROW-MAJOR, indexed by `packedTriangleIndex`. One packing
// convention for the whole fit, defined once in FitRecipe.hpp, so the recipe that fills a triangle
// and the algebra that reads it back cannot drift apart.
//
// WHEN `positiveDefinite` IS FALSE, `factorPacked` HOLDS A PARTIAL FACTORISATION AND MUST NOT BE
// USED. Every function below checks the flag and answers NaN rather than reading it.
struct CholeskyFactor final
{
    CholeskyStatus status = CholeskyStatus::Ok;

    // The one gate every function below gets past. Derived rather than stored, so it cannot fall out
    // of step with the status.
    bool positiveDefinite() const noexcept { return this->status == CholeskyStatus::Ok; }

    // The parameter whose pivot failed, zero-based. Equal to `terms` when nothing failed, which is
    // one past the last parameter and so cannot be mistaken for one.
    size_t failedParameter = 0;

    size_t terms = 0;

    // U, packed upper triangle, row-major. `packedTriangleSize(terms)` entries.
    std::vector<double> factorPacked;
};

// Factorise a symmetric matrix given as a packed upper triangle.
//
// k = 0 is a legitimate input and means there is nothing to solve: the result is positive definite,
// empty, and every function below answers the empty or zero thing for it. A fit with no parameters
// is a real case -- it is the one that must agree with `aev()` -- so it must not be a crash.
//
// Throws only on a shape error, which is a caller bug rather than a property of the data.
inline CholeskyFactor choleskyFactorPacked(const std::vector<double>& packedUpper, size_t terms)
{
    if (packedUpper.size() != packedTriangleSize(terms))
    {
        throw std::invalid_argument("veil: packed upper triangle of the wrong length to factorise.");
    }

    CholeskyFactor factor;
    factor.terms = terms;
    factor.failedParameter = terms;
    factor.factorPacked.assign(packedTriangleSize(terms), 0.0);

    for (size_t row = 0; row < terms; ++row)
    {
        double pivot = packedUpper[packedTriangleIndex(row, row, terms)];
        for (size_t above = 0; above < row; ++above)
        {
            const double entry = factor.factorPacked[packedTriangleIndex(above, row, terms)];
            pivot -= entry * entry;
        }

        // `!(pivot > 0.0)` RATHER THAN `pivot <= 0.0`, AND THE DIFFERENCE IS NOT COSMETIC. Every
        // comparison against NaN is false, so `pivot <= 0.0` would wave a NaN pivot through, and
        // `std::sqrt(NaN)` is NaN, and the whole factor would fill with NaN while still claiming to
        // be positive definite. The package's policy is that IEEE NaN is the only missing-value
        // marker and that it flows through the engine unchecked, so a NaN genuinely does arrive
        // here. Written this way the negation makes NaN fail, which is the honest answer.
        //
        // A NaN anywhere in the triangle reaches a diagonal: an off-diagonal at (row, column) enters
        // the pivot for `column`, so nothing can hide in a corner.
        const double originalDiagonal = packedUpper[packedTriangleIndex(row, row, terms)];
        if (!(pivot > 0.0))
        {
            // A NaN fails both of these, so it must be asked about FIRST or it would be reported as
            // an ordinary singularity and send the user hunting for a collinear covariate that does
            // not exist.
            factor.status = std::isfinite(pivot) ? CholeskyStatus::NotPositiveDefinite
                                                 : CholeskyStatus::NotFinite;
            factor.failedParameter = row;
            return factor;
        }
        if (pivot < NearlySingularTolerance * originalDiagonal)
        {
            // Positive, but so little of this covariate is its own that the inverse would be
            // meaningless. The partial factor is kept: `choleskyDependency` reads it to say WHICH
            // earlier covariates account for this one.
            factor.status = CholeskyStatus::NearlySingular;
            factor.failedParameter = row;
            return factor;
        }

        const double diagonal = std::sqrt(pivot);
        factor.factorPacked[packedTriangleIndex(row, row, terms)] = diagonal;

        for (size_t column = row + 1; column < terms; ++column)
        {
            double sum = packedUpper[packedTriangleIndex(row, column, terms)];
            for (size_t above = 0; above < row; ++above)
            {
                sum -= factor.factorPacked[packedTriangleIndex(above, row, terms)]
                    * factor.factorPacked[packedTriangleIndex(above, column, terms)];
            }
            factor.factorPacked[packedTriangleIndex(row, column, terms)] = sum / diagonal;
        }
    }

    return factor;
}

// Solve A x = b, in two triangular substitutions: U^T y = b forwards, then U x = y backwards.
//
// Answers a vector of NaN when the factorisation failed, rather than throwing or returning something
// plausible. Throws on a length mismatch, which is a caller bug.
inline std::vector<double> choleskySolve(
    const CholeskyFactor& factor,
    const std::vector<double>& rightHandSide)
{
    if (rightHandSide.size() != factor.terms)
    {
        throw std::invalid_argument("veil: right hand side of the wrong length to solve.");
    }

    const size_t terms = factor.terms;
    std::vector<double> solution(terms, std::numeric_limits<double>::quiet_NaN());
    if (!factor.positiveDefinite()) { return solution; }

    // Forward: U^T is lower triangular with (U^T)[row][column] = U[column][row].
    for (size_t row = 0; row < terms; ++row)
    {
        double sum = rightHandSide[row];
        for (size_t column = 0; column < row; ++column)
        {
            sum -= factor.factorPacked[packedTriangleIndex(column, row, terms)] * solution[column];
        }
        solution[row] = sum / factor.factorPacked[packedTriangleIndex(row, row, terms)];
    }

    // Backward, in place: the forward pass left y where x is going.
    for (size_t row = terms; row-- > 0;)
    {
        double sum = solution[row];
        for (size_t column = row + 1; column < terms; ++column)
        {
            sum -= factor.factorPacked[packedTriangleIndex(row, column, terms)] * solution[column];
        }
        solution[row] = sum / factor.factorPacked[packedTriangleIndex(row, row, terms)];
    }

    return solution;
}

// V = U^-1, upper triangular, in the same packed layout. A building block for the inverse rather
// than an answer anyone wants on its own, but exposed because it is worth being able to test.
//
// Answers NaN throughout when the factorisation failed.
inline std::vector<double> choleskyTriangularInversePacked(const CholeskyFactor& factor)
{
    const size_t terms = factor.terms;
    std::vector<double> triangular(
        packedTriangleSize(terms), std::numeric_limits<double>::quiet_NaN());
    if (!factor.positiveDefinite()) { return triangular; }

    for (size_t row = 0; row < terms; ++row)
    {
        triangular[packedTriangleIndex(row, row, terms)] =
            1.0 / factor.factorPacked[packedTriangleIndex(row, row, terms)];
    }

    // Column by column, upwards: V[row][column] needs the entries below it in the same column.
    for (size_t column = 0; column < terms; ++column)
    {
        for (size_t row = column; row-- > 0;)
        {
            double sum = 0.0;
            for (size_t between = row + 1; between <= column; ++between)
            {
                sum += factor.factorPacked[packedTriangleIndex(row, between, terms)]
                    * triangular[packedTriangleIndex(between, column, terms)];
            }
            triangular[packedTriangleIndex(row, column, terms)] =
                -sum / factor.factorPacked[packedTriangleIndex(row, row, terms)];
        }
    }

    return triangular;
}

// A^-1 as a packed upper triangle, the same layout in and out.
//
// A = U^T U, so A^-1 = U^-1 U^-T = V V^T with V = U^-1. Both operands of that product are upper
// triangular, so entry (row, column) sums only over p >= column.
//
// Answers NaN throughout when the factorisation failed.
inline std::vector<double> choleskyInversePacked(const CholeskyFactor& factor)
{
    const size_t terms = factor.terms;
    std::vector<double> inverse(packedTriangleSize(terms), std::numeric_limits<double>::quiet_NaN());
    if (!factor.positiveDefinite()) { return inverse; }

    const std::vector<double> triangular = choleskyTriangularInversePacked(factor);

    for (size_t row = 0; row < terms; ++row)
    {
        for (size_t column = row; column < terms; ++column)
        {
            double sum = 0.0;
            for (size_t between = column; between < terms; ++between)
            {
                sum += triangular[packedTriangleIndex(row, between, terms)]
                    * triangular[packedTriangleIndex(column, between, terms)];
            }
            inverse[packedTriangleIndex(row, column, terms)] = sum;
        }
    }

    return inverse;
}

// v^T A^-1 v, which is the Newton decrement lambda^2 when v is the score.
//
// A^-1 = U^-1 U^-T, so v^T A^-1 v = ||U^-T v||^2 and only the forward substitution is needed. That
// makes it half the work of `choleskySolve`, and it also makes it manifestly non-negative for a
// positive-definite A, where forming grad^T delta from the full solve leaves the sign at the mercy
// of cancellation. It agrees with grad^T delta to rounding, not to the bit.
//
// Zero for k = 0, which is the right answer: an empty quadratic form is an empty sum. NaN when the
// factorisation failed. Throws on a length mismatch.
inline double choleskyQuadraticForm(const CholeskyFactor& factor, const std::vector<double>& values)
{
    if (values.size() != factor.terms)
    {
        throw std::invalid_argument("veil: vector of the wrong length for the quadratic form.");
    }

    const size_t terms = factor.terms;
    if (!factor.positiveDefinite()) { return std::numeric_limits<double>::quiet_NaN(); }

    std::vector<double> forward(terms, 0.0);
    double total = 0.0;
    for (size_t row = 0; row < terms; ++row)
    {
        double sum = values[row];
        for (size_t column = 0; column < row; ++column)
        {
            sum -= factor.factorPacked[packedTriangleIndex(column, row, terms)] * forward[column];
        }
        forward[row] = sum / factor.factorPacked[packedTriangleIndex(row, row, terms)];
        total += forward[row] * forward[row];
    }

    return total;
}

// tr(B A^-1) for a symmetric B given as a packed upper triangle -- the penalty term tr(J I^-1), up
// to the Z that lives outside this file.
//
// Both matrices are symmetric, so the trace collapses to their Frobenius inner product: the
// diagonal once and every off-diagonal twice. No product matrix is formed.
//
// Zero for k = 0. NaN when the factorisation failed. Throws on a length mismatch.
inline double choleskyTraceOfProductWithInverse(
    const CholeskyFactor& factor,
    const std::vector<double>& packedUpperOther)
{
    if (packedUpperOther.size() != packedTriangleSize(factor.terms))
    {
        throw std::invalid_argument("veil: packed upper triangle of the wrong length for the trace.");
    }

    if (!factor.positiveDefinite()) { return std::numeric_limits<double>::quiet_NaN(); }

    const size_t terms = factor.terms;
    const std::vector<double> inverse = choleskyInversePacked(factor);

    double total = 0.0;
    for (size_t row = 0; row < terms; ++row)
    {
        const size_t diagonal = packedTriangleIndex(row, row, terms);
        total += packedUpperOther[diagonal] * inverse[diagonal];
        for (size_t column = row + 1; column < terms; ++column)
        {
            const size_t slot = packedTriangleIndex(row, column, terms);
            total += 2.0 * packedUpperOther[slot] * inverse[slot];
        }
    }

    return total;
}

// THE SANDWICH `A^-1 B A^-1`, packed upper, which is what `Var(beta_hat)` is.
//
// With `A = Ew XX^T` and `B = Ew^2 XX^T` the raw moments, `I = Omega^-1 A` and `J = Omega^-1 B`, so
//
//     Var(beta_hat) = I^-1 J I^-1 = Omega * A^-1 B A^-1
//
// and the Omega belongs to the caller, exactly as it does everywhere else in the fit. This function
// answers the bracket alone.
//
// Symmetric by construction, since B is symmetric and `A^-1` is too, so only the upper triangle is
// formed. NaN-filled when the factorisation failed, like every other reader of the factor.
inline std::vector<double> choleskySandwich(
    const CholeskyFactor& factor,
    const std::vector<double>& packedUpperOther)
{
    const size_t terms = factor.terms;
    if (packedUpperOther.size() != packedTriangleSize(terms))
    {
        throw std::invalid_argument("veil: the sandwich's middle matrix is the wrong length.");
    }

    std::vector<double> result(packedTriangleSize(terms), std::numeric_limits<double>::quiet_NaN());
    if (!factor.positiveDefinite()) { return result; }

    const std::vector<double> inverse = choleskyInversePacked(factor);
    const auto at = [terms](const std::vector<double>& packed, size_t row, size_t column)
    {
        return packed[packedTriangleIndex(row, column, terms)];
    };

    // `middle = A^-1 B`, which is NOT symmetric and so is held in full.
    std::vector<double> middle(terms * terms, 0.0);
    for (size_t row = 0; row < terms; ++row)
    {
        for (size_t column = 0; column < terms; ++column)
        {
            double total = 0.0;
            for (size_t inner = 0; inner < terms; ++inner)
            {
                total += at(inverse, row, inner) * at(packedUpperOther, inner, column);
            }
            middle[row * terms + column] = total;
        }
    }

    for (size_t row = 0; row < terms; ++row)
    {
        for (size_t column = row; column < terms; ++column)
        {
            double total = 0.0;
            for (size_t inner = 0; inner < terms; ++inner)
            {
                total += middle[row * terms + inner] * at(inverse, inner, column);
            }
            result[packedTriangleIndex(row, column, terms)] = total;
        }
    }

    return result;
}

// WHICH EARLIER COVARIATES ACCOUNT FOR THE ONE THAT FAILED.
//
// A pivot index says WHERE a dependency was detected, not WHAT depends on what, and the two are
// routinely different: with covariates 0 and 2 identical the factorisation fails at 2, because 0 and
// 1 are still independent. Reporting the index alone sends a user to the wrong column half the time.
//
// The coefficients are one back-substitution against the factor already computed, so this costs
// O(k^2) and no new decomposition. Failing at column j, the entries `U[0..j-1][j]` were written by
// the earlier rows and satisfy `U^T y = A[0..j-1][j]`, while the dependency z is defined by
// `A[0..j-1][0..j-1] z = A[0..j-1][j]`, which is `U^T U z = U^T y`, which is `U z = y`. So z is the
// back-substitution of the leading block against a column already in hand.
//
// The answer reads: column `failedParameter` is approximately the sum over i of `z[i]` times column
// i. Empty when the failure was at parameter 0, which means that covariate is degenerate on its own
// -- no exposure, or identically zero -- rather than a duplicate of anything.
//
// MEANINGLESS FOR `NotFinite`, and answers empty for it: a NaN says nothing about collinearity.
inline std::vector<double> choleskyDependency(const CholeskyFactor& factor)
{
    if (factor.status == CholeskyStatus::Ok || factor.status == CholeskyStatus::NotFinite)
    {
        return std::vector<double>();
    }

    const size_t terms = factor.terms;
    const size_t failed = factor.failedParameter;
    std::vector<double> coefficients(failed, 0.0);

    // Back-substitute U[0..failed-1][0..failed-1] z = U[0..failed-1][failed].
    for (size_t step = failed; step-- > 0;)
    {
        double sum = factor.factorPacked[packedTriangleIndex(step, failed, terms)];
        for (size_t column = step + 1; column < failed; ++column)
        {
            sum -= factor.factorPacked[packedTriangleIndex(step, column, terms)] * coefficients[column];
        }
        const double diagonal = factor.factorPacked[packedTriangleIndex(step, step, terms)];
        coefficients[step] = sum / diagonal;
    }

    return coefficients;
}

} // namespace veil
