// Group 3 (approximation): least-squares fitting.
//
// Every case draws a point set through one of the MakePoints helpers below,
// which alternate a uniform cloud with exactly representable lattice sets and
// with constructions aimed at a specific branch (all points coincident, a
// collinear or coplanar set, an exact cocircular/cospherical set, a set with
// a degenerate x-domain). Every helper records the same layout: the point
// count, the mode index, then the point coordinates.
//
// The fits are long arithmetic pipelines (covariance accumulation, Gaussian
// elimination, symmetric eigensolvers), but none of them reaches the C math
// library except through std::sqrt: SymmetricEigensolver2x2, the *iterative*
// SymmetricEigensolver3x3 (the non-iterative NISymmetricEigensolver3x3, which
// uses std::acos and std::cos, is a different class and no fitter here calls
// it), SymmetricEigensolver<Real> (Householder plus implicit-shift QL) and
// RootsPolynomial::Find (Cauchy bound plus bisection) all use only + - * /
// and sqrt. Every case is therefore declared exact on the TypeScript side.
#define ORACLE_FAMILY "v03-approximation"
#include "Oracle.h"

#include <Mathematics/ApprCircle2.h>
#include <Mathematics/ApprCurveByArcs.h>
#include <Mathematics/ApprGaussian2.h>
#include <Mathematics/ApprGaussian3.h>
#include <Mathematics/ApprHeightLine2.h>
#include <Mathematics/ApprHeightPlane3.h>
#include <Mathematics/ApprOrthogonalLine2.h>
#include <Mathematics/ApprOrthogonalLine3.h>
#include <Mathematics/ApprOrthogonalPlane3.h>
#include <Mathematics/ApprParallelLines2.h>
#include <Mathematics/ApprPolynomial2.h>
#include <Mathematics/ApprPolynomial3.h>
#include <Mathematics/ApprPolynomial4.h>
#include <Mathematics/ApprPolynomialSpecial2.h>
#include <Mathematics/ApprPolynomialSpecial3.h>
#include <Mathematics/ApprPolynomialSpecial4.h>
#include <Mathematics/ApprQuadratic2.h>
#include <Mathematics/ApprQuadratic3.h>
#include <Mathematics/ApprQuery.h>
#include <Mathematics/ApprSphere3.h>
#include <Mathematics/BezierCurve.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

using namespace gte;

namespace
{
    // The twelve lattice points of the circle of radius 5 about the origin.
    constexpr int32_t kCircleLattice[12][2] =
    {
        { 5, 0 }, { -5, 0 }, { 0, 5 }, { 0, -5 },
        { 3, 4 }, { 3, -4 }, { -3, 4 }, { -3, -4 },
        { 4, 3 }, { 4, -3 }, { -4, 3 }, { -4, -3 }
    };

    // Lattice points of the sphere of radius 5 about the origin.
    constexpr int32_t kSphereLattice[30][3] =
    {
        { 5, 0, 0 }, { -5, 0, 0 }, { 0, 5, 0 }, { 0, -5, 0 }, { 0, 0, 5 }, { 0, 0, -5 },
        { 3, 4, 0 }, { 3, -4, 0 }, { -3, 4, 0 }, { -3, -4, 0 },
        { 4, 3, 0 }, { 4, -3, 0 }, { -4, 3, 0 }, { -4, -3, 0 },
        { 3, 0, 4 }, { 3, 0, -4 }, { -3, 0, 4 }, { -3, 0, -4 },
        { 4, 0, 3 }, { 4, 0, -3 }, { -4, 0, 3 }, { -4, 0, -3 },
        { 0, 3, 4 }, { 0, 3, -4 }, { 0, -3, 4 }, { 0, -3, -4 },
        { 0, 4, 3 }, { 0, 4, -3 }, { 0, -4, 3 }, { 0, -4, -3 }
    };

    // A 2D point set. Records the count, the mode and then 2 doubles per
    // point. Modes:
    //   0  uniform cloud
    //   1  small integer lattice
    //   2  every point the same lattice point (zero covariance)
    //   3  collinear on an integer lattice line
    //   4  exactly cocircular: radius 5 about an integer center
    std::vector<Vector2<double>> MakePoints2(oracle::Ctx& io, int32_t minPoints,
        int32_t maxPoints, int32_t numModes = 5)
    {
        int32_t const n = io.integer(minPoints, maxPoints);
        int32_t const mode = io.integer(0, numModes - 1);
        std::vector<Vector2<double>> P(static_cast<size_t>(n));

        if (mode == 0)
        {
            for (int32_t i = 0; i < n; ++i) { P[i] = io.vec<2>(-10.0, 10.0); }
        }
        else if (mode == 1)
        {
            for (int32_t i = 0; i < n; ++i) { P[i] = io.latticeVec<2>(-6, 6); }
        }
        else if (mode == 2)
        {
            Vector2<double> q{ static_cast<double>(io.rawInteger(-6, 6)),
                static_cast<double>(io.rawInteger(-6, 6)) };
            for (int32_t i = 0; i < n; ++i) { P[i] = io.givenVec(q); }
        }
        else if (mode == 3)
        {
            Vector2<double> base{ static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) };
            Vector2<double> dir{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            for (int32_t i = 0; i < n; ++i)
            {
                double const t = static_cast<double>(io.rawInteger(-4, 4));
                P[i] = io.givenVec(Vector2<double>{ base[0] + t * dir[0], base[1] + t * dir[1] });
            }
        }
        else
        {
            Vector2<double> center{ static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) };
            for (int32_t i = 0; i < n; ++i)
            {
                int32_t const k = io.rawInteger(0, 11);
                P[i] = io.givenVec(Vector2<double>{
                    center[0] + static_cast<double>(kCircleLattice[k][0]),
                    center[1] + static_cast<double>(kCircleLattice[k][1]) });
            }
        }
        return P;
    }

    // A 3D point set. Modes:
    //   0  uniform cloud
    //   1  small integer lattice
    //   2  every point the same lattice point
    //   3  collinear on an integer lattice line
    //   4  coplanar on an integer lattice plane
    //   5  exactly cospherical: radius 5 about an integer center
    std::vector<Vector3<double>> MakePoints3(oracle::Ctx& io, int32_t minPoints,
        int32_t maxPoints, int32_t numModes = 6)
    {
        int32_t const n = io.integer(minPoints, maxPoints);
        int32_t const mode = io.integer(0, numModes - 1);
        std::vector<Vector3<double>> P(static_cast<size_t>(n));

        if (mode == 0)
        {
            for (int32_t i = 0; i < n; ++i) { P[i] = io.vec<3>(-10.0, 10.0); }
        }
        else if (mode == 1)
        {
            for (int32_t i = 0; i < n; ++i) { P[i] = io.latticeVec<3>(-6, 6); }
        }
        else if (mode == 2)
        {
            Vector3<double> q{ static_cast<double>(io.rawInteger(-6, 6)),
                static_cast<double>(io.rawInteger(-6, 6)),
                static_cast<double>(io.rawInteger(-6, 6)) };
            for (int32_t i = 0; i < n; ++i) { P[i] = io.givenVec(q); }
        }
        else if (mode == 3)
        {
            Vector3<double> base{ static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) };
            Vector3<double> dir{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            for (int32_t i = 0; i < n; ++i)
            {
                double const t = static_cast<double>(io.rawInteger(-4, 4));
                P[i] = io.givenVec(Vector3<double>{ base[0] + t * dir[0],
                    base[1] + t * dir[1], base[2] + t * dir[2] });
            }
        }
        else if (mode == 4)
        {
            Vector3<double> base{ static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) };
            Vector3<double> u{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            Vector3<double> v{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            for (int32_t i = 0; i < n; ++i)
            {
                double const s = static_cast<double>(io.rawInteger(-3, 3));
                double const t = static_cast<double>(io.rawInteger(-3, 3));
                P[i] = io.givenVec(Vector3<double>{
                    base[0] + s * u[0] + t * v[0],
                    base[1] + s * u[1] + t * v[1],
                    base[2] + s * u[2] + t * v[2] });
            }
        }
        else
        {
            Vector3<double> center{ static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) };
            for (int32_t i = 0; i < n; ++i)
            {
                int32_t const k = io.rawInteger(0, 29);
                P[i] = io.givenVec(Vector3<double>{
                    center[0] + static_cast<double>(kSphereLattice[k][0]),
                    center[1] + static_cast<double>(kSphereLattice[k][1]),
                    center[2] + static_cast<double>(kSphereLattice[k][2]) });
            }
        }
        return P;
    }

    // An index list into a set of n observations. Records the index count and
    // then the indices. Modes: identity, reversed, repeated, random.
    std::vector<int32_t> MakeIndices(oracle::Ctx& io, int32_t n, int32_t minIndices)
    {
        int32_t const count = io.integer(minIndices, n);
        int32_t const mode = io.integer(0, 3);
        std::vector<int32_t> indices(static_cast<size_t>(count));
        for (int32_t i = 0; i < count; ++i)
        {
            int32_t value;
            if (mode == 0) { value = i; }
            else if (mode == 1) { value = count - 1 - i; }
            else if (mode == 2) { value = i % 2; }
            else { value = io.rawInteger(0, n - 1); }
            indices[i] = static_cast<int32_t>(io.given(static_cast<double>(value)));
        }
        return indices;
    }

    // Observations (x,w). Records the count, the mode and then 2 doubles per
    // observation. Modes:
    //   0  uniform
    //   1  small integer lattice
    //   2  degenerate x-domain (every sample shares its x)
    //   3  exact polynomial samples w = sum c[i]*x^i on a lattice
    //   4  degenerate w-domain (every sample shares its w)
    std::vector<std::array<double, 2>> MakeObs2(oracle::Ctx& io, int32_t minN,
        int32_t maxN)
    {
        int32_t const n = io.integer(minN, maxN);
        int32_t const mode = io.integer(0, 4);
        std::vector<std::array<double, 2>> obs(static_cast<size_t>(n));
        double const fixedX = static_cast<double>(io.rawInteger(-4, 4));
        double const fixedW = static_cast<double>(io.rawInteger(-4, 4));
        std::array<double, 4> c{};
        for (int32_t i = 0; i < 4; ++i)
        {
            c[i] = static_cast<double>(io.rawInteger(-3, 3));
        }

        for (int32_t i = 0; i < n; ++i)
        {
            double x, w;
            if (mode == 0)
            {
                x = io.raw(-5.0, 5.0);
                w = io.raw(-5.0, 5.0);
            }
            else if (mode == 1)
            {
                x = static_cast<double>(io.rawInteger(-5, 5));
                w = static_cast<double>(io.rawInteger(-5, 5));
            }
            else if (mode == 2)
            {
                x = fixedX;
                w = static_cast<double>(io.rawInteger(-5, 5));
            }
            else if (mode == 3)
            {
                x = static_cast<double>(io.rawInteger(-3, 3));
                w = ((c[3] * x + c[2]) * x + c[1]) * x + c[0];
            }
            else
            {
                x = static_cast<double>(io.rawInteger(-5, 5));
                w = fixedW;
            }
            obs[i][0] = io.given(x);
            obs[i][1] = io.given(w);
        }
        return obs;
    }

    // Observations (x,y,w). Modes as for MakeObs2, with mode 3 sampling the
    // bilinear polynomial c0 + c1*x + c2*y + c3*x*y exactly.
    std::vector<std::array<double, 3>> MakeObs3(oracle::Ctx& io, int32_t minN,
        int32_t maxN)
    {
        int32_t const n = io.integer(minN, maxN);
        int32_t const mode = io.integer(0, 4);
        std::vector<std::array<double, 3>> obs(static_cast<size_t>(n));
        double const fixedX = static_cast<double>(io.rawInteger(-4, 4));
        double const fixedW = static_cast<double>(io.rawInteger(-4, 4));
        std::array<double, 4> c{};
        for (int32_t i = 0; i < 4; ++i)
        {
            c[i] = static_cast<double>(io.rawInteger(-3, 3));
        }

        for (int32_t i = 0; i < n; ++i)
        {
            double x, y, w;
            if (mode == 0)
            {
                x = io.raw(-5.0, 5.0);
                y = io.raw(-5.0, 5.0);
                w = io.raw(-5.0, 5.0);
            }
            else if (mode == 1)
            {
                x = static_cast<double>(io.rawInteger(-4, 4));
                y = static_cast<double>(io.rawInteger(-4, 4));
                w = static_cast<double>(io.rawInteger(-4, 4));
            }
            else if (mode == 2)
            {
                x = fixedX;
                y = static_cast<double>(io.rawInteger(-4, 4));
                w = static_cast<double>(io.rawInteger(-4, 4));
            }
            else if (mode == 3)
            {
                x = static_cast<double>(io.rawInteger(-3, 3));
                y = static_cast<double>(io.rawInteger(-3, 3));
                w = c[0] + c[1] * x + c[2] * y + c[3] * x * y;
            }
            else
            {
                x = static_cast<double>(io.rawInteger(-4, 4));
                y = static_cast<double>(io.rawInteger(-4, 4));
                w = fixedW;
            }
            obs[i][0] = io.given(x);
            obs[i][1] = io.given(y);
            obs[i][2] = io.given(w);
        }
        return obs;
    }

    // Observations (x,y,z,w). Modes as for MakeObs2, with mode 3 sampling the
    // affine polynomial c0 + c1*x + c2*y + c3*z exactly.
    std::vector<std::array<double, 4>> MakeObs4(oracle::Ctx& io, int32_t minN,
        int32_t maxN)
    {
        int32_t const n = io.integer(minN, maxN);
        int32_t const mode = io.integer(0, 4);
        std::vector<std::array<double, 4>> obs(static_cast<size_t>(n));
        double const fixedX = static_cast<double>(io.rawInteger(-3, 3));
        double const fixedW = static_cast<double>(io.rawInteger(-3, 3));
        std::array<double, 4> c{};
        for (int32_t i = 0; i < 4; ++i)
        {
            c[i] = static_cast<double>(io.rawInteger(-3, 3));
        }

        for (int32_t i = 0; i < n; ++i)
        {
            double x, y, z, w;
            if (mode == 0)
            {
                x = io.raw(-4.0, 4.0);
                y = io.raw(-4.0, 4.0);
                z = io.raw(-4.0, 4.0);
                w = io.raw(-4.0, 4.0);
            }
            else if (mode == 1)
            {
                x = static_cast<double>(io.rawInteger(-3, 3));
                y = static_cast<double>(io.rawInteger(-3, 3));
                z = static_cast<double>(io.rawInteger(-3, 3));
                w = static_cast<double>(io.rawInteger(-3, 3));
            }
            else if (mode == 2)
            {
                x = fixedX;
                y = static_cast<double>(io.rawInteger(-3, 3));
                z = static_cast<double>(io.rawInteger(-3, 3));
                w = static_cast<double>(io.rawInteger(-3, 3));
            }
            else if (mode == 3)
            {
                x = static_cast<double>(io.rawInteger(-2, 2));
                y = static_cast<double>(io.rawInteger(-2, 2));
                z = static_cast<double>(io.rawInteger(-2, 2));
                w = c[0] + c[1] * x + c[2] * y + c[3] * z;
            }
            else
            {
                x = static_cast<double>(io.rawInteger(-3, 3));
                y = static_cast<double>(io.rawInteger(-3, 3));
                z = static_cast<double>(io.rawInteger(-3, 3));
                w = fixedW;
            }
            obs[i][0] = io.given(x);
            obs[i][1] = io.given(y);
            obs[i][2] = io.given(z);
            obs[i][3] = io.given(w);
        }
        return obs;
    }

    // A strictly increasing nonnegative degree list with gaps. Records the
    // count and then the degrees.
    std::vector<int32_t> MakeDegrees(oracle::Ctx& io, int32_t maxCount)
    {
        int32_t const count = io.integer(1, maxCount);
        std::vector<int32_t> degrees(static_cast<size_t>(count));
        int32_t last = -1;
        for (int32_t i = 0; i < count; ++i)
        {
            last += io.rawInteger(1, 3);
            degrees[i] = static_cast<int32_t>(io.given(static_cast<double>(last)));
        }
        return degrees;
    }
}

// ---------------------------------------------------------------- ApprQuery

// The base-class Fit(observations) overload (the std::vector form; the
// (count, pointer) form builds the same identity index list).
ORACLE_CASE("ApprQuery.fit.all")
{
    auto P = MakePoints2(io, 1, 8);
    ApprHeightLine2<double> fitter;
    bool success = fitter.Fit(P);
    io.outBool(success);
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
}

// Fit(observations, imin, imax) over a contiguous subset, including the
// imin > imax rejection.
ORACLE_CASE("ApprQuery.fit.range")
{
    auto P = MakePoints2(io, 2, 8);
    int32_t imin = io.integer(0, static_cast<int32_t>(P.size()) - 1);
    int32_t imax = io.integer(0, static_cast<int32_t>(P.size()) - 1);
    ApprHeightLine2<double> fitter;
    bool success = fitter.Fit(P, static_cast<size_t>(imin), static_cast<size_t>(imax));
    io.outBool(success);
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
}

// Fit(observations, indices) with repeated and out-of-order indices.
ORACLE_CASE("ApprQuery.fit.indexed")
{
    auto P = MakePoints2(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(P.size()), 1);
    ApprHeightLine2<double> fitter;
    bool success = fitter.Fit(P, indices);
    io.outBool(success);
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
}

// Fit(observations, indices, numIndices): the prefix is clamped with
// std::min(numIndices, indices.size()), so a numIndices past the end is legal.
ORACLE_CASE("ApprQuery.fit.numIndices")
{
    auto P = MakePoints2(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(P.size()), 1);
    int32_t numIndices = io.integer(0, static_cast<int32_t>(indices.size()) + 2);
    ApprHeightLine2<double> fitter;
    bool success = fitter.Fit(P, indices, static_cast<size_t>(numIndices));
    io.outBool(success);
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
}

// FitIndexed called directly, plus GetMinimumRequired and Error.
ORACLE_CASE("ApprQuery.fitIndexed.direct")
{
    auto P = MakePoints2(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(P.size()), 1);
    auto probe = io.vec<2>(-10.0, 10.0);
    ApprHeightLine2<double> fitter;
    bool success = fitter.FitIndexed(P.size(), P.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
    io.outReal(fitter.Error(probe));
}

// RANSAC with too few observations: an immediate false, bestConsensus
// untouched.
ORACLE_CASE("ApprQuery.ransac.tooFew")
{
    auto P = MakePoints2(io, 1, 1);
    double maxError = io.real(0.0, 4.0);
    int32_t numIterations = io.integer(1, 4);
    ApprGaussian2<double> candidate, best;
    std::vector<int32_t> consensus;
    bool success = ApprQuery<double, Vector2<double>>::RANSAC(candidate, P, 1,
        maxError, static_cast<size_t>(numIterations), consensus, best);
    io.outBool(success);
    io.outInt(consensus.size());
    io.outVec(best.GetParameters().center);
    io.outVec(best.GetParameters().extent);
}

// RANSAC with exactly GetMinimumRequired() observations: the shuffle loop is
// skipped entirely, so this path is deterministic on both sides (see the
// group report for why the shuffling path is not).
ORACLE_CASE("ApprQuery.ransac.minimum")
{
    auto P = MakePoints2(io, 2, 2);
    double maxError = io.real(0.0, 4.0);
    int32_t numIterations = io.integer(1, 4);
    auto probe = io.vec<2>(-6.0, 6.0);
    ApprGaussian2<double> candidate, best;
    std::vector<int32_t> consensus;
    bool success = ApprQuery<double, Vector2<double>>::RANSAC(candidate, P, 2,
        maxError, static_cast<size_t>(numIterations), consensus, best);
    io.outBool(success);
    io.outInt(consensus.size());
    for (auto index : consensus) { io.outInt(index); }
    io.outVec(best.GetParameters().center);
    io.outVec(best.GetParameters().axis[0]);
    io.outVec(best.GetParameters().axis[1]);
    io.outVec(best.GetParameters().extent);
    io.outReal(best.Error(probe));
}

// -------------------------------------------------------------- ApprCircle2

ORACLE_CASE("ApprCircle2.fitUsingSquaredLengths")
{
    auto P = MakePoints2(io, 1, 10);
    Circle2<double> circle{};
    ApprCircle2<double> fitter;
    bool success = fitter.FitUsingSquaredLengths(static_cast<int32_t>(P.size()),
        P.data(), circle);
    io.outBool(success);
    io.outVec(circle.center);
    io.outReal(circle.radius);
}

// Points exactly on the radius-5 lattice circle: the fit is exact, and
// repeated points drive the determinant to zero.
ORACLE_CASE("ApprCircle2.fitUsingSquaredLengths.cocircular")
{
    int32_t const n = io.integer(3, 8);
    Vector2<double> center{ static_cast<double>(io.rawInteger(-4, 4)),
        static_cast<double>(io.rawInteger(-4, 4)) };
    std::vector<Vector2<double>> P(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i)
    {
        int32_t const k = io.rawInteger(0, 11);
        P[i] = io.givenVec(Vector2<double>{
            center[0] + static_cast<double>(kCircleLattice[k][0]),
            center[1] + static_cast<double>(kCircleLattice[k][1]) });
    }
    Circle2<double> circle{};
    ApprCircle2<double> fitter;
    bool success = fitter.FitUsingSquaredLengths(n, P.data(), circle);
    io.outBool(success);
    io.outVec(circle.center);
    io.outReal(circle.radius);
}

ORACLE_CASE("ApprCircle2.fitUsingLengths")
{
    auto P = MakePoints2(io, 1, 8);
    int32_t maxIterations = io.integer(0, 6);
    bool initialCenterIsAverage = io.boolean();
    auto initialCenter = io.vec<2>(-8.0, 8.0);
    double initialRadius = io.real(0.5, 4.0);
    Circle2<double> circle(initialCenter, initialRadius);
    ApprCircle2<double> fitter;
    uint32_t iterations = fitter.FitUsingLengths(static_cast<int32_t>(P.size()),
        P.data(), static_cast<uint32_t>(maxIterations), initialCenterIsAverage,
        circle);
    io.outInt(iterations);
    io.outVec(circle.center);
    io.outReal(circle.radius);
}

// A positive epsilon so the convergence break is reachable.
ORACLE_CASE("ApprCircle2.fitUsingLengths.epsilon")
{
    auto P = MakePoints2(io, 2, 8);
    int32_t maxIterations = io.integer(1, 12);
    bool initialCenterIsAverage = io.boolean();
    auto initialCenter = io.vec<2>(-8.0, 8.0);
    double initialRadius = io.real(0.5, 4.0);
    double epsilon = io.real(0.0, 0.25);
    Circle2<double> circle(initialCenter, initialRadius);
    ApprCircle2<double> fitter;
    uint32_t iterations = fitter.FitUsingLengths(static_cast<int32_t>(P.size()),
        P.data(), static_cast<uint32_t>(maxIterations), initialCenterIsAverage,
        circle, epsilon);
    io.outInt(iterations);
    io.outVec(circle.center);
    io.outReal(circle.radius);
}

// -------------------------------------------------------------- ApprSphere3

ORACLE_CASE("ApprSphere3.fitUsingSquaredLengths")
{
    auto P = MakePoints3(io, 1, 10);
    Sphere3<double> sphere{};
    ApprSphere3<double> fitter;
    bool success = fitter.FitUsingSquaredLengths(static_cast<int32_t>(P.size()),
        P.data(), sphere);
    io.outBool(success);
    io.outVec(sphere.center);
    io.outReal(sphere.radius);
}

ORACLE_CASE("ApprSphere3.fitUsingSquaredLengths.cospherical")
{
    int32_t const n = io.integer(4, 10);
    Vector3<double> center{ static_cast<double>(io.rawInteger(-4, 4)),
        static_cast<double>(io.rawInteger(-4, 4)),
        static_cast<double>(io.rawInteger(-4, 4)) };
    std::vector<Vector3<double>> P(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i)
    {
        int32_t const k = io.rawInteger(0, 29);
        P[i] = io.givenVec(Vector3<double>{
            center[0] + static_cast<double>(kSphereLattice[k][0]),
            center[1] + static_cast<double>(kSphereLattice[k][1]),
            center[2] + static_cast<double>(kSphereLattice[k][2]) });
    }
    Sphere3<double> sphere{};
    ApprSphere3<double> fitter;
    bool success = fitter.FitUsingSquaredLengths(n, P.data(), sphere);
    io.outBool(success);
    io.outVec(sphere.center);
    io.outReal(sphere.radius);
}

ORACLE_CASE("ApprSphere3.fitUsingLengths")
{
    auto P = MakePoints3(io, 1, 8);
    int32_t maxIterations = io.integer(0, 6);
    bool initialCenterIsAverage = io.boolean();
    auto initialCenter = io.vec<3>(-8.0, 8.0);
    double initialRadius = io.real(0.5, 4.0);
    Sphere3<double> sphere(initialCenter, initialRadius);
    ApprSphere3<double> fitter;
    uint32_t iterations = fitter.FitUsingLengths(static_cast<int32_t>(P.size()),
        P.data(), static_cast<uint32_t>(maxIterations), initialCenterIsAverage,
        sphere);
    io.outInt(iterations);
    io.outVec(sphere.center);
    io.outReal(sphere.radius);
}

ORACLE_CASE("ApprSphere3.fitUsingLengths.epsilon")
{
    auto P = MakePoints3(io, 2, 8);
    int32_t maxIterations = io.integer(1, 12);
    bool initialCenterIsAverage = io.boolean();
    auto initialCenter = io.vec<3>(-8.0, 8.0);
    double initialRadius = io.real(0.5, 4.0);
    double epsilon = io.real(0.0, 0.25);
    Sphere3<double> sphere(initialCenter, initialRadius);
    ApprSphere3<double> fitter;
    uint32_t iterations = fitter.FitUsingLengths(static_cast<int32_t>(P.size()),
        P.data(), static_cast<uint32_t>(maxIterations), initialCenterIsAverage,
        sphere, epsilon);
    io.outInt(iterations);
    io.outVec(sphere.center);
    io.outReal(sphere.radius);
}

// ------------------------------------------------------------ ApprGaussian2

ORACLE_CASE("ApprGaussian2.fit")
{
    auto P = MakePoints2(io, 1, 10);
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprGaussian2<double> fitter;
    bool success = fitter.Fit(P);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outVec(fitter.GetParameters().center);
    io.outVec(fitter.GetParameters().axis[0]);
    io.outVec(fitter.GetParameters().axis[1]);
    io.outVec(fitter.GetParameters().extent);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprGaussian2.fitIndexed")
{
    auto P = MakePoints2(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(P.size()), 1);
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprGaussian2<double> fitter;
    bool success = fitter.FitIndexed(P.size(), P.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outVec(fitter.GetParameters().center);
    io.outVec(fitter.GetParameters().axis[0]);
    io.outVec(fitter.GetParameters().axis[1]);
    io.outVec(fitter.GetParameters().extent);
    io.outReal(fitter.Error(probe));
}

// CopyParameters through the base-class pointer.
ORACLE_CASE("ApprGaussian2.copyParameters")
{
    auto P = MakePoints2(io, 2, 8);
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprGaussian2<double> source, target;
    source.Fit(P);
    target.CopyParameters(&source);
    io.outVec(target.GetParameters().center);
    io.outVec(target.GetParameters().axis[0]);
    io.outVec(target.GetParameters().axis[1]);
    io.outVec(target.GetParameters().extent);
    io.outReal(target.Error(probe));
}

// ------------------------------------------------------------ ApprGaussian3

ORACLE_CASE("ApprGaussian3.fit")
{
    auto P = MakePoints3(io, 1, 10);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprGaussian3<double> fitter;
    bool success = fitter.Fit(P);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outVec(fitter.GetParameters().center);
    io.outVec(fitter.GetParameters().axis[0]);
    io.outVec(fitter.GetParameters().axis[1]);
    io.outVec(fitter.GetParameters().axis[2]);
    io.outVec(fitter.GetParameters().extent);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprGaussian3.fitIndexed")
{
    auto P = MakePoints3(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(P.size()), 1);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprGaussian3<double> fitter;
    bool success = fitter.FitIndexed(P.size(), P.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outVec(fitter.GetParameters().center);
    io.outVec(fitter.GetParameters().axis[0]);
    io.outVec(fitter.GetParameters().axis[1]);
    io.outVec(fitter.GetParameters().axis[2]);
    io.outVec(fitter.GetParameters().extent);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprGaussian3.copyParameters")
{
    auto P = MakePoints3(io, 2, 8);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprGaussian3<double> source, target;
    source.Fit(P);
    target.CopyParameters(&source);
    io.outVec(target.GetParameters().center);
    io.outVec(target.GetParameters().axis[0]);
    io.outVec(target.GetParameters().axis[1]);
    io.outVec(target.GetParameters().axis[2]);
    io.outVec(target.GetParameters().extent);
    io.outReal(target.Error(probe));
}

// ----------------------------------------------------------- ApprHeightLine2

ORACLE_CASE("ApprHeightLine2.fit")
{
    auto P = MakePoints2(io, 1, 10);
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprHeightLine2<double> fitter;
    bool success = fitter.Fit(P);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
    io.outReal(fitter.Error(probe));
}

// A degenerate x-domain (every sample shares its x) makes covar00 exactly
// zero, which is the documented failure branch.
ORACLE_CASE("ApprHeightLine2.fit.verticalData")
{
    int32_t const n = io.integer(1, 8);
    double const x = static_cast<double>(io.rawInteger(-5, 5));
    std::vector<Vector2<double>> P(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i)
    {
        double const y = static_cast<double>(io.rawInteger(-5, 5));
        P[i] = io.givenVec(Vector2<double>{ x, y });
    }
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprHeightLine2<double> fitter;
    bool success = fitter.Fit(P);
    io.outBool(success);
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprHeightLine2.copyParameters")
{
    auto P = MakePoints2(io, 2, 8);
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprHeightLine2<double> source, target;
    source.Fit(P);
    target.CopyParameters(&source);
    io.outVec(target.GetParameters().first);
    io.outVec(target.GetParameters().second);
    io.outReal(target.Error(probe));
}

// ---------------------------------------------------------- ApprHeightPlane3

ORACLE_CASE("ApprHeightPlane3.fit")
{
    auto P = MakePoints3(io, 1, 10);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprHeightPlane3<double> fitter;
    bool success = fitter.Fit(P);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
    io.outReal(fitter.Error(probe));
}

// Samples whose (x,y) projections are collinear on a lattice line make the
// 2x2 determinant exactly zero, which is the documented failure branch.
ORACLE_CASE("ApprHeightPlane3.fit.collinearXY")
{
    int32_t const n = io.integer(1, 8);
    Vector2<double> base{ static_cast<double>(io.rawInteger(-4, 4)),
        static_cast<double>(io.rawInteger(-4, 4)) };
    Vector2<double> dir{ static_cast<double>(io.rawInteger(-3, 3)),
        static_cast<double>(io.rawInteger(-3, 3)) };
    std::vector<Vector3<double>> P(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i)
    {
        double const t = static_cast<double>(io.rawInteger(-4, 4));
        double const z = static_cast<double>(io.rawInteger(-5, 5));
        P[i] = io.givenVec(Vector3<double>{ base[0] + t * dir[0],
            base[1] + t * dir[1], z });
    }
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprHeightPlane3<double> fitter;
    bool success = fitter.Fit(P);
    io.outBool(success);
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprHeightPlane3.copyParameters")
{
    auto P = MakePoints3(io, 3, 8);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprHeightPlane3<double> source, target;
    source.Fit(P);
    target.CopyParameters(&source);
    io.outVec(target.GetParameters().first);
    io.outVec(target.GetParameters().second);
    io.outReal(target.Error(probe));
}

// ------------------------------------------------------- ApprOrthogonalLine2

// The return value is 'eval[0] < eval[1]', so the lattice and coincident
// modes (which produce exactly equal eigenvalues) reach the false branch.
ORACLE_CASE("ApprOrthogonalLine2.fit")
{
    auto P = MakePoints2(io, 1, 10);
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprOrthogonalLine2<double> fitter;
    bool unique = fitter.Fit(P);
    io.outBool(unique);
    io.outInt(fitter.GetMinimumRequired());
    io.outVec(fitter.GetParameters().origin);
    io.outVec(fitter.GetParameters().direction);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprOrthogonalLine2.fitIndexed")
{
    auto P = MakePoints2(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(P.size()), 1);
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprOrthogonalLine2<double> fitter;
    bool unique = fitter.FitIndexed(P.size(), P.data(), indices.size(), indices.data());
    io.outBool(unique);
    io.outVec(fitter.GetParameters().origin);
    io.outVec(fitter.GetParameters().direction);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprOrthogonalLine2.copyParameters")
{
    auto P = MakePoints2(io, 2, 8);
    auto probe = io.vec<2>(-12.0, 12.0);
    ApprOrthogonalLine2<double> source, target;
    source.Fit(P);
    target.CopyParameters(&source);
    io.outVec(target.GetParameters().origin);
    io.outVec(target.GetParameters().direction);
    io.outReal(target.Error(probe));
}

// ------------------------------------------------------- ApprOrthogonalLine3

ORACLE_CASE("ApprOrthogonalLine3.fit")
{
    auto P = MakePoints3(io, 1, 10);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprOrthogonalLine3<double> fitter;
    bool unique = fitter.Fit(P);
    io.outBool(unique);
    io.outInt(fitter.GetMinimumRequired());
    io.outVec(fitter.GetParameters().origin);
    io.outVec(fitter.GetParameters().direction);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprOrthogonalLine3.fitIndexed")
{
    auto P = MakePoints3(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(P.size()), 1);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprOrthogonalLine3<double> fitter;
    bool unique = fitter.FitIndexed(P.size(), P.data(), indices.size(), indices.data());
    io.outBool(unique);
    io.outVec(fitter.GetParameters().origin);
    io.outVec(fitter.GetParameters().direction);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprOrthogonalLine3.copyParameters")
{
    auto P = MakePoints3(io, 2, 8);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprOrthogonalLine3<double> source, target;
    source.Fit(P);
    target.CopyParameters(&source);
    io.outVec(target.GetParameters().origin);
    io.outVec(target.GetParameters().direction);
    io.outReal(target.Error(probe));
}

// ------------------------------------------------------ ApprOrthogonalPlane3

ORACLE_CASE("ApprOrthogonalPlane3.fit")
{
    auto P = MakePoints3(io, 1, 10);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprOrthogonalPlane3<double> fitter;
    bool unique = fitter.Fit(P);
    io.outBool(unique);
    io.outInt(fitter.GetMinimumRequired());
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprOrthogonalPlane3.fitIndexed")
{
    auto P = MakePoints3(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(P.size()), 1);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprOrthogonalPlane3<double> fitter;
    bool unique = fitter.FitIndexed(P.size(), P.data(), indices.size(), indices.data());
    io.outBool(unique);
    io.outVec(fitter.GetParameters().first);
    io.outVec(fitter.GetParameters().second);
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprOrthogonalPlane3.copyParameters")
{
    auto P = MakePoints3(io, 3, 8);
    auto probe = io.vec<3>(-12.0, 12.0);
    ApprOrthogonalPlane3<double> source, target;
    source.Fit(P);
    target.CopyParameters(&source);
    io.outVec(target.GetParameters().first);
    io.outVec(target.GetParameters().second);
    io.outReal(target.Error(probe));
}

// ---------------------------------------------------------- ApprPolynomial2

ORACLE_CASE("ApprPolynomial2.fit")
{
    auto obs = MakeObs2(io, 1, 8);
    int32_t degree = io.integer(1, 3);
    double probeX = io.real(-6.0, 6.0);
    std::array<double, 2> probe{ io.real(-6.0, 6.0), io.real(-6.0, 6.0) };
    ApprPolynomial2<double> fitter(degree);
    bool success = fitter.Fit(obs);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.Evaluate(probeX));
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprPolynomial2.fitIndexed")
{
    auto obs = MakeObs2(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(obs.size()), 1);
    int32_t degree = io.integer(1, 3);
    double probeX = io.real(-6.0, 6.0);
    ApprPolynomial2<double> fitter(degree);
    bool success = fitter.FitIndexed(obs.size(), obs.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.Evaluate(probeX));
}

// CopyParameters carries the degree, the domain and the coefficients.
ORACLE_CASE("ApprPolynomial2.copyParameters")
{
    auto obs = MakeObs2(io, 2, 8);
    int32_t degree = io.integer(1, 3);
    double probeX = io.real(-6.0, 6.0);
    ApprPolynomial2<double> source(degree), target(degree);
    source.Fit(obs);
    target.CopyParameters(&source);
    io.outInt(target.GetParameters().size());
    for (auto value : target.GetParameters()) { io.outReal(value); }
    io.outReal(target.GetXDomain()[0]);
    io.outReal(target.GetXDomain()[1]);
    io.outReal(target.Evaluate(probeX));
}

// ---------------------------------------------------------- ApprPolynomial3

ORACLE_CASE("ApprPolynomial3.fit")
{
    auto obs = MakeObs3(io, 1, 10);
    int32_t xDegree = io.integer(0, 2);
    int32_t yDegree = io.integer(0, 2);
    double probeX = io.real(-6.0, 6.0);
    double probeY = io.real(-6.0, 6.0);
    std::array<double, 3> probe{ io.real(-6.0, 6.0), io.real(-6.0, 6.0), io.real(-6.0, 6.0) };
    ApprPolynomial3<double> fitter(xDegree, yDegree);
    bool success = fitter.Fit(obs);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.GetYDomain()[0]);
    io.outReal(fitter.GetYDomain()[1]);
    io.outReal(fitter.Evaluate(probeX, probeY));
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprPolynomial3.fitIndexed")
{
    auto obs = MakeObs3(io, 2, 10);
    auto indices = MakeIndices(io, static_cast<int32_t>(obs.size()), 1);
    int32_t xDegree = io.integer(0, 2);
    int32_t yDegree = io.integer(0, 2);
    double probeX = io.real(-6.0, 6.0);
    double probeY = io.real(-6.0, 6.0);
    ApprPolynomial3<double> fitter(xDegree, yDegree);
    bool success = fitter.FitIndexed(obs.size(), obs.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.GetYDomain()[0]);
    io.outReal(fitter.GetYDomain()[1]);
    io.outReal(fitter.Evaluate(probeX, probeY));
}

ORACLE_CASE("ApprPolynomial3.copyParameters")
{
    auto obs = MakeObs3(io, 2, 10);
    int32_t xDegree = io.integer(0, 2);
    int32_t yDegree = io.integer(0, 2);
    double probeX = io.real(-6.0, 6.0);
    double probeY = io.real(-6.0, 6.0);
    ApprPolynomial3<double> source(xDegree, yDegree), target(xDegree, yDegree);
    source.Fit(obs);
    target.CopyParameters(&source);
    io.outInt(target.GetParameters().size());
    for (auto value : target.GetParameters()) { io.outReal(value); }
    io.outReal(target.GetXDomain()[0]);
    io.outReal(target.GetXDomain()[1]);
    io.outReal(target.GetYDomain()[0]);
    io.outReal(target.GetYDomain()[1]);
    io.outReal(target.Evaluate(probeX, probeY));
}

// ---------------------------------------------------------- ApprPolynomial4

ORACLE_CASE("ApprPolynomial4.fit")
{
    auto obs = MakeObs4(io, 1, 12);
    int32_t xDegree = io.integer(0, 2);
    int32_t yDegree = io.integer(0, 2);
    int32_t zDegree = io.integer(0, 2);
    double probeX = io.real(-5.0, 5.0);
    double probeY = io.real(-5.0, 5.0);
    double probeZ = io.real(-5.0, 5.0);
    std::array<double, 4> probe{ io.real(-5.0, 5.0), io.real(-5.0, 5.0),
        io.real(-5.0, 5.0), io.real(-5.0, 5.0) };
    ApprPolynomial4<double> fitter(xDegree, yDegree, zDegree);
    bool success = fitter.Fit(obs);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.GetYDomain()[0]);
    io.outReal(fitter.GetYDomain()[1]);
    io.outReal(fitter.GetZDomain()[0]);
    io.outReal(fitter.GetZDomain()[1]);
    io.outReal(fitter.Evaluate(probeX, probeY, probeZ));
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprPolynomial4.fitIndexed")
{
    auto obs = MakeObs4(io, 2, 12);
    auto indices = MakeIndices(io, static_cast<int32_t>(obs.size()), 1);
    int32_t xDegree = io.integer(0, 2);
    int32_t yDegree = io.integer(0, 1);
    int32_t zDegree = io.integer(0, 1);
    double probeX = io.real(-5.0, 5.0);
    double probeY = io.real(-5.0, 5.0);
    double probeZ = io.real(-5.0, 5.0);
    ApprPolynomial4<double> fitter(xDegree, yDegree, zDegree);
    bool success = fitter.FitIndexed(obs.size(), obs.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.GetYDomain()[0]);
    io.outReal(fitter.GetYDomain()[1]);
    io.outReal(fitter.GetZDomain()[0]);
    io.outReal(fitter.GetZDomain()[1]);
    io.outReal(fitter.Evaluate(probeX, probeY, probeZ));
}

ORACLE_CASE("ApprPolynomial4.copyParameters")
{
    auto obs = MakeObs4(io, 2, 12);
    int32_t xDegree = io.integer(0, 1);
    int32_t yDegree = io.integer(0, 1);
    int32_t zDegree = io.integer(0, 1);
    double probeX = io.real(-5.0, 5.0);
    double probeY = io.real(-5.0, 5.0);
    double probeZ = io.real(-5.0, 5.0);
    ApprPolynomial4<double> source(xDegree, yDegree, zDegree);
    ApprPolynomial4<double> target(xDegree, yDegree, zDegree);
    source.Fit(obs);
    target.CopyParameters(&source);
    io.outInt(target.GetParameters().size());
    for (auto value : target.GetParameters()) { io.outReal(value); }
    io.outReal(target.GetXDomain()[0]);
    io.outReal(target.GetXDomain()[1]);
    io.outReal(target.Evaluate(probeX, probeY, probeZ));
}

// --------------------------------------------------- ApprPolynomialSpecial2

ORACLE_CASE("ApprPolynomialSpecial2.fit")
{
    auto obs = MakeObs2(io, 1, 8);
    auto degrees = MakeDegrees(io, 3);
    double probeX = io.real(-6.0, 6.0);
    std::array<double, 2> probe{ io.real(-6.0, 6.0), io.real(-6.0, 6.0) };
    ApprPolynomialSpecial2<double> fitter(degrees);
    bool success = fitter.Fit(obs);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.Evaluate(probeX));
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprPolynomialSpecial2.fitIndexed")
{
    auto obs = MakeObs2(io, 2, 8);
    auto indices = MakeIndices(io, static_cast<int32_t>(obs.size()), 1);
    auto degrees = MakeDegrees(io, 3);
    double probeX = io.real(-6.0, 6.0);
    ApprPolynomialSpecial2<double> fitter(degrees);
    bool success = fitter.FitIndexed(obs.size(), obs.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.Evaluate(probeX));
}

ORACLE_CASE("ApprPolynomialSpecial2.copyParameters")
{
    auto obs = MakeObs2(io, 2, 8);
    auto degrees = MakeDegrees(io, 3);
    double probeX = io.real(-6.0, 6.0);
    ApprPolynomialSpecial2<double> source(degrees), target(degrees);
    source.Fit(obs);
    target.CopyParameters(&source);
    io.outInt(target.GetParameters().size());
    for (auto value : target.GetParameters()) { io.outReal(value); }
    io.outReal(target.GetXDomain()[0]);
    io.outReal(target.GetXDomain()[1]);
    io.outReal(target.Evaluate(probeX));
}

// The constructor's LogAssert on a degree list that is not strictly
// increasing (the port preserves both asserts).
ORACLE_CASE("ApprPolynomialSpecial2.constructor.assert")
{
    int32_t count = io.integer(1, 3);
    std::vector<int32_t> degrees(static_cast<size_t>(count));
    for (int32_t i = 0; i < count; ++i)
    {
        degrees[i] = static_cast<int32_t>(io.given(
            static_cast<double>(io.rawInteger(-1, 2))));
    }
    ApprPolynomialSpecial2<double> fitter(degrees);
    io.outInt(fitter.GetMinimumRequired());
}

// --------------------------------------------------- ApprPolynomialSpecial3

ORACLE_CASE("ApprPolynomialSpecial3.fit")
{
    auto obs = MakeObs3(io, 1, 10);
    auto xDegrees = MakeDegrees(io, 3);
    double probeX = io.real(-6.0, 6.0);
    double probeY = io.real(-6.0, 6.0);
    std::array<double, 3> probe{ io.real(-6.0, 6.0), io.real(-6.0, 6.0), io.real(-6.0, 6.0) };
    // The upstream constructor requires the two lists to have equal size and
    // each to be strictly increasing (the preserved defect below), so the y
    // list is the identity list of the same size.
    std::vector<int32_t> yDegrees(xDegrees.size());
    for (size_t i = 0; i < yDegrees.size(); ++i)
    {
        yDegrees[i] = static_cast<int32_t>(i);
    }
    ApprPolynomialSpecial3<double> fitter(xDegrees, yDegrees);
    bool success = fitter.Fit(obs);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.GetYDomain()[0]);
    io.outReal(fitter.GetYDomain()[1]);
    io.outReal(fitter.Evaluate(probeX, probeY));
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprPolynomialSpecial3.fitIndexed")
{
    auto obs = MakeObs3(io, 2, 10);
    auto indices = MakeIndices(io, static_cast<int32_t>(obs.size()), 1);
    auto xDegrees = MakeDegrees(io, 3);
    double probeX = io.real(-6.0, 6.0);
    double probeY = io.real(-6.0, 6.0);
    std::vector<int32_t> yDegrees(xDegrees.size());
    for (size_t i = 0; i < yDegrees.size(); ++i)
    {
        yDegrees[i] = static_cast<int32_t>(2 * i);
    }
    ApprPolynomialSpecial3<double> fitter(xDegrees, yDegrees);
    bool success = fitter.FitIndexed(obs.size(), obs.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.GetYDomain()[0]);
    io.outReal(fitter.GetYDomain()[1]);
    io.outReal(fitter.Evaluate(probeX, probeY));
}

ORACLE_CASE("ApprPolynomialSpecial3.copyParameters")
{
    auto obs = MakeObs3(io, 2, 10);
    auto xDegrees = MakeDegrees(io, 3);
    double probeX = io.real(-6.0, 6.0);
    double probeY = io.real(-6.0, 6.0);
    std::vector<int32_t> yDegrees(xDegrees.size());
    for (size_t i = 0; i < yDegrees.size(); ++i)
    {
        yDegrees[i] = static_cast<int32_t>(i);
    }
    ApprPolynomialSpecial3<double> source(xDegrees, yDegrees);
    ApprPolynomialSpecial3<double> target(xDegrees, yDegrees);
    source.Fit(obs);
    target.CopyParameters(&source);
    io.outInt(target.GetParameters().size());
    for (auto value : target.GetParameters()) { io.outReal(value); }
    io.outReal(target.GetXDomain()[0]);
    io.outReal(target.GetXDomain()[1]);
    io.outReal(target.GetYDomain()[0]);
    io.outReal(target.GetYDomain()[1]);
    io.outReal(target.Evaluate(probeX, probeY));
}

// The preserved upstream defect: the constructor asserts that each degree
// list is *separately* strictly increasing, which rejects the affine model
// {1, x, y} that the header documents as admissible.
ORACLE_CASE("ApprPolynomialSpecial3.constructor.affineAssert")
{
    int32_t count = io.integer(2, 3);
    std::vector<int32_t> xDegrees(static_cast<size_t>(count));
    std::vector<int32_t> yDegrees(static_cast<size_t>(count));
    for (int32_t i = 0; i < count; ++i)
    {
        xDegrees[i] = static_cast<int32_t>(io.given(static_cast<double>(i == 0 ? 0 : 1)));
        yDegrees[i] = static_cast<int32_t>(io.given(static_cast<double>(i)));
    }
    ApprPolynomialSpecial3<double> fitter(xDegrees, yDegrees);
    io.outInt(fitter.GetMinimumRequired());
}

// --------------------------------------------------- ApprPolynomialSpecial4

ORACLE_CASE("ApprPolynomialSpecial4.fit")
{
    auto obs = MakeObs4(io, 1, 12);
    auto xDegrees = MakeDegrees(io, 3);
    double probeX = io.real(-5.0, 5.0);
    double probeY = io.real(-5.0, 5.0);
    double probeZ = io.real(-5.0, 5.0);
    std::array<double, 4> probe{ io.real(-5.0, 5.0), io.real(-5.0, 5.0),
        io.real(-5.0, 5.0), io.real(-5.0, 5.0) };
    std::vector<int32_t> yDegrees(xDegrees.size()), zDegrees(xDegrees.size());
    for (size_t i = 0; i < xDegrees.size(); ++i)
    {
        yDegrees[i] = static_cast<int32_t>(i);
        zDegrees[i] = static_cast<int32_t>(2 * i);
    }
    ApprPolynomialSpecial4<double> fitter(xDegrees, yDegrees, zDegrees);
    bool success = fitter.Fit(obs);
    io.outBool(success);
    io.outInt(fitter.GetMinimumRequired());
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.GetYDomain()[0]);
    io.outReal(fitter.GetYDomain()[1]);
    io.outReal(fitter.GetZDomain()[0]);
    io.outReal(fitter.GetZDomain()[1]);
    io.outReal(fitter.Evaluate(probeX, probeY, probeZ));
    io.outReal(fitter.Error(probe));
}

ORACLE_CASE("ApprPolynomialSpecial4.fitIndexed")
{
    auto obs = MakeObs4(io, 2, 12);
    auto indices = MakeIndices(io, static_cast<int32_t>(obs.size()), 1);
    auto xDegrees = MakeDegrees(io, 3);
    double probeX = io.real(-5.0, 5.0);
    double probeY = io.real(-5.0, 5.0);
    double probeZ = io.real(-5.0, 5.0);
    std::vector<int32_t> yDegrees(xDegrees.size()), zDegrees(xDegrees.size());
    for (size_t i = 0; i < xDegrees.size(); ++i)
    {
        yDegrees[i] = static_cast<int32_t>(2 * i);
        zDegrees[i] = static_cast<int32_t>(i);
    }
    ApprPolynomialSpecial4<double> fitter(xDegrees, yDegrees, zDegrees);
    bool success = fitter.FitIndexed(obs.size(), obs.data(), indices.size(), indices.data());
    io.outBool(success);
    io.outInt(fitter.GetParameters().size());
    for (auto value : fitter.GetParameters()) { io.outReal(value); }
    io.outReal(fitter.GetXDomain()[0]);
    io.outReal(fitter.GetXDomain()[1]);
    io.outReal(fitter.GetYDomain()[0]);
    io.outReal(fitter.GetYDomain()[1]);
    io.outReal(fitter.GetZDomain()[0]);
    io.outReal(fitter.GetZDomain()[1]);
    io.outReal(fitter.Evaluate(probeX, probeY, probeZ));
}

ORACLE_CASE("ApprPolynomialSpecial4.copyParameters")
{
    auto obs = MakeObs4(io, 2, 12);
    auto xDegrees = MakeDegrees(io, 3);
    double probeX = io.real(-5.0, 5.0);
    double probeY = io.real(-5.0, 5.0);
    double probeZ = io.real(-5.0, 5.0);
    std::vector<int32_t> yDegrees(xDegrees.size()), zDegrees(xDegrees.size());
    for (size_t i = 0; i < xDegrees.size(); ++i)
    {
        yDegrees[i] = static_cast<int32_t>(i);
        zDegrees[i] = static_cast<int32_t>(i);
    }
    ApprPolynomialSpecial4<double> source(xDegrees, yDegrees, zDegrees);
    ApprPolynomialSpecial4<double> target(xDegrees, yDegrees, zDegrees);
    source.Fit(obs);
    target.CopyParameters(&source);
    io.outInt(target.GetParameters().size());
    for (auto value : target.GetParameters()) { io.outReal(value); }
    io.outReal(target.GetXDomain()[0]);
    io.outReal(target.GetXDomain()[1]);
    io.outReal(target.Evaluate(probeX, probeY, probeZ));
}

// The 'input arrays must have the same size' assert.
ORACLE_CASE("ApprPolynomialSpecial4.constructor.sizeAssert")
{
    int32_t count = io.integer(1, 3);
    int32_t extra = io.integer(0, 2);
    std::vector<int32_t> xDegrees(static_cast<size_t>(count));
    std::vector<int32_t> yDegrees(static_cast<size_t>(count));
    std::vector<int32_t> zDegrees(static_cast<size_t>(count + extra));
    for (int32_t i = 0; i < count; ++i) { xDegrees[i] = i; yDegrees[i] = i; }
    for (int32_t i = 0; i < count + extra; ++i) { zDegrees[i] = i; }
    ApprPolynomialSpecial4<double> fitter(xDegrees, yDegrees, zDegrees);
    io.outInt(fitter.GetMinimumRequired());
}
