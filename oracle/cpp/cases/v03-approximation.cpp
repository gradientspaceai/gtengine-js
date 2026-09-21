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
#include <Mathematics/SymmetricEigensolver.h>

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
    // Drawn but NOT recorded, so a rejection loop can redraw every quantity
    // its acceptance test depends on. RecordPoints2 writes the accepted draw.
    struct Points2
    {
        int32_t n = 0, mode = 0;
        std::vector<Vector2<double>> P;
    };

    Points2 DrawPoints2(oracle::Ctx& io, int32_t minPoints, int32_t maxPoints,
        int32_t numModes = 5, int32_t forcedMode = -1)
    {
        Points2 out;
        out.n = io.rawInteger(minPoints, maxPoints);
        out.mode = io.rawInteger(0, numModes - 1);
        if (forcedMode >= 0) { out.mode = forcedMode; }
        out.P.resize(static_cast<size_t>(out.n));
        int32_t const n = out.n;

        if (out.mode == 0)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                out.P[i][0] = io.raw(-10.0, 10.0);
                out.P[i][1] = io.raw(-10.0, 10.0);
            }
        }
        else if (out.mode == 1)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                out.P[i][0] = static_cast<double>(io.rawInteger(-6, 6));
                out.P[i][1] = static_cast<double>(io.rawInteger(-6, 6));
            }
        }
        else if (out.mode == 2)
        {
            Vector2<double> q{ static_cast<double>(io.rawInteger(-6, 6)),
                static_cast<double>(io.rawInteger(-6, 6)) };
            for (int32_t i = 0; i < n; ++i) { out.P[i] = q; }
        }
        else if (out.mode == 3)
        {
            Vector2<double> base{ static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) };
            Vector2<double> dir{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            for (int32_t i = 0; i < n; ++i)
            {
                double const t = static_cast<double>(io.rawInteger(-4, 4));
                out.P[i] = Vector2<double>{ base[0] + t * dir[0], base[1] + t * dir[1] };
            }
        }
        else
        {
            Vector2<double> center{ static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) };
            for (int32_t i = 0; i < n; ++i)
            {
                int32_t const k = io.rawInteger(0, 11);
                out.P[i] = Vector2<double>{
                    center[0] + static_cast<double>(kCircleLattice[k][0]),
                    center[1] + static_cast<double>(kCircleLattice[k][1]) };
            }
        }
        return out;
    }

    std::vector<Vector2<double>> RecordPoints2(oracle::Ctx& io, Points2 const& draw)
    {
        io.given(static_cast<double>(draw.n));
        io.given(static_cast<double>(draw.mode));
        for (auto const& p : draw.P) { io.givenVec(p); }
        return draw.P;
    }

    std::vector<Vector2<double>> MakePoints2(oracle::Ctx& io, int32_t minPoints,
        int32_t maxPoints, int32_t numModes = 5)
    {
        return RecordPoints2(io, DrawPoints2(io, minPoints, maxPoints, numModes));
    }

    // A 3D point set. Modes:
    //   0  uniform cloud
    //   1  small integer lattice
    //   2  every point the same lattice point
    //   3  collinear on an integer lattice line
    //   4  coplanar on an integer lattice plane
    //   5  exactly cospherical: radius 5 about an integer center
    struct Points3
    {
        int32_t n = 0, mode = 0;
        std::vector<Vector3<double>> P;
    };

    Points3 DrawPoints3(oracle::Ctx& io, int32_t minPoints, int32_t maxPoints,
        int32_t numModes = 6, int32_t forcedMode = -1)
    {
        Points3 out;
        out.n = io.rawInteger(minPoints, maxPoints);
        out.mode = io.rawInteger(0, numModes - 1);
        if (forcedMode >= 0) { out.mode = forcedMode; }
        out.P.resize(static_cast<size_t>(out.n));
        int32_t const n = out.n, mode = out.mode;
        std::vector<Vector3<double>>& P = out.P;

        if (mode == 0)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                P[i][0] = io.raw(-10.0, 10.0);
                P[i][1] = io.raw(-10.0, 10.0);
                P[i][2] = io.raw(-10.0, 10.0);
            }
        }
        else if (mode == 1)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                P[i][0] = static_cast<double>(io.rawInteger(-6, 6));
                P[i][1] = static_cast<double>(io.rawInteger(-6, 6));
                P[i][2] = static_cast<double>(io.rawInteger(-6, 6));
            }
        }
        else if (mode == 2)
        {
            Vector3<double> q{ static_cast<double>(io.rawInteger(-6, 6)),
                static_cast<double>(io.rawInteger(-6, 6)),
                static_cast<double>(io.rawInteger(-6, 6)) };
            for (int32_t i = 0; i < n; ++i) { P[i] = q; }
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
                P[i] = Vector3<double>{ base[0] + t * dir[0],
                    base[1] + t * dir[1], base[2] + t * dir[2] };
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
                P[i] = Vector3<double>{
                    base[0] + s * u[0] + t * v[0],
                    base[1] + s * u[1] + t * v[1],
                    base[2] + s * u[2] + t * v[2] };
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
                P[i] = Vector3<double>{
                    center[0] + static_cast<double>(kSphereLattice[k][0]),
                    center[1] + static_cast<double>(kSphereLattice[k][1]),
                    center[2] + static_cast<double>(kSphereLattice[k][2]) };
            }
        }
        return out;
    }

    std::vector<Vector3<double>> RecordPoints3(oracle::Ctx& io, Points3 const& draw)
    {
        io.given(static_cast<double>(draw.n));
        io.given(static_cast<double>(draw.mode));
        for (auto const& p : draw.P) { io.givenVec(p); }
        return draw.P;
    }

    std::vector<Vector3<double>> MakePoints3(oracle::Ctx& io, int32_t minPoints,
        int32_t maxPoints, int32_t numModes = 6)
    {
        return RecordPoints3(io, DrawPoints3(io, minPoints, maxPoints, numModes));
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

// ----------------------------------------------------------- ApprQuadratic2
//
// The four ApprQuadratic entry points minimize C^T M C over unit C by taking
// the eigenvector of the smallest eigenvalue of M from SymmetricEigensolver
// (the NxN Householder + implicit-shift QL solver). The port deliberately
// fixes a result-corrupting defect of that solver: when a Householder step is
// degenerate (the subcolumn below the subdiagonal is already zero, so
// 'length == 0' and the reflection actually applied is the identity),
// upstream still stores the reflection parameter 2/Dot(v,v) == 2, and
// GetEigenvector rebuilds H = I - 2*e*e^T instead of the identity. The
// eigenvalues are unaffected (the stored value lands in the lower triangle,
// which the diagonal copy never reads) but Q is corrupted.
// docs/UPSTREAM-FINDINGS.md, issue #80.
//
// EigenDecouples below replicates upstream's Tridiagonalize verbatim and
// reports whether any step is degenerate. That is an exact characterization
// of the inputs on which the two implementations differ: the fix changes only
// the value stored for a degenerate step. The main cases reject those inputs
// by rejection sampling; the .decoupledDeviation cases keep only them.

namespace
{
    // A verbatim copy of SymmetricEigensolver<double>::Tridiagonalize that
    // returns true as soon as a degenerate Householder step occurs. 'matrix'
    // is a copy, since the routine overwrites it.
    bool EigenDecouples(std::vector<double> matrix, int32_t size)
    {
        std::vector<double> v(size, 0.0), p(size, 0.0), w(size, 0.0);
        int32_t r, c;
        for (int32_t i = 0, ip1 = 1; i < size - 2; ++i, ++ip1)
        {
            double length = 0.0;
            for (r = 0; r < ip1; ++r) { v[r] = 0.0; }
            for (r = ip1; r < size; ++r)
            {
                double vr = matrix[r + static_cast<size_t>(size) * i];
                v[r] = vr;
                length += vr * vr;
            }
            double vdv = 1.0;
            length = std::sqrt(length);
            if (length > 0.0)
            {
                double v1 = v[ip1];
                double sgn = (v1 >= 0.0 ? 1.0 : -1.0);
                double invDenom = 1.0 / (v1 + sgn * length);
                v[ip1] = 1.0;
                for (r = ip1 + 1; r < size; ++r)
                {
                    v[r] *= invDenom;
                    vdv += v[r] * v[r];
                }
            }
            else
            {
                // The degenerate step: upstream stores 2 for a reflection
                // that is the identity.
                return true;
            }

            double invvdv = 1.0 / vdv;
            double twoinvvdv = invvdv * 2.0;
            double pdvtvdv = 0.0;
            for (r = i; r < size; ++r)
            {
                p[r] = 0.0;
                for (c = i; c < r; ++c)
                {
                    p[r] += matrix[r + static_cast<size_t>(size) * c] * v[c];
                }
                for (/**/; c < size; ++c)
                {
                    p[r] += matrix[c + static_cast<size_t>(size) * r] * v[c];
                }
                p[r] *= twoinvvdv;
                pdvtvdv += p[r] * v[r];
            }

            pdvtvdv *= invvdv;
            for (r = i; r < size; ++r) { w[r] = p[r] - pdvtvdv * v[r]; }

            for (r = i; r < size; ++r)
            {
                double vr = v[r], wr = w[r];
                double offset = vr * wr * 2.0;
                matrix[r + static_cast<size_t>(size) * r] -= offset;
                for (c = r + 1; c < size; ++c)
                {
                    offset = vr * w[c] + wr * v[c];
                    matrix[c + static_cast<size_t>(size) * r] -= offset;
                }
            }

            matrix[i + static_cast<size_t>(size) * ip1] = twoinvvdv;
            for (r = ip1 + 1; r < size; ++r)
            {
                matrix[i + static_cast<size_t>(size) * r] = v[r];
            }
        }
        return false;
    }

    // A verbatim copy of ApprQuadratic2::operator()'s matrix assembly.
    std::vector<double> BuildQuadratic2M(std::vector<Vector2<double>> const& P)
    {
        int32_t const numPoints = static_cast<int32_t>(P.size());
        Matrix<6, 6, double> M{};
        for (int32_t i = 0; i < numPoints; ++i)
        {
            double x = P[i][0];
            double y = P[i][1];
            double x2 = x * x;
            double y2 = y * y;
            double xy = x * y;
            double x3 = x * x2;
            double xy2 = x * y2;
            double x2y = x * xy;
            double y3 = y * y2;
            double x4 = x * x3;
            double x2y2 = x * xy2;
            double x3y = x * x2y;
            double y4 = y * y3;
            double xy3 = x * y3;

            M(0, 1) += x;
            M(0, 2) += y;
            M(0, 3) += x2;
            M(0, 4) += xy;
            M(0, 5) += y2;
            M(1, 3) += x3;
            M(1, 4) += x2y;
            M(1, 5) += xy2;
            M(2, 5) += y3;
            M(3, 3) += x4;
            M(3, 4) += x3y;
            M(3, 5) += x2y2;
            M(4, 5) += xy3;
            M(5, 5) += y4;
        }

        double const rNumPoints = static_cast<double>(numPoints);
        M(0, 0) = rNumPoints;
        M(1, 1) = M(0, 3);
        M(1, 2) = M(0, 4);
        M(2, 2) = M(0, 5);
        M(2, 3) = M(1, 4);
        M(2, 4) = M(1, 5);
        M(4, 4) = M(3, 5);
        for (int32_t row = 0; row < 6; ++row)
        {
            for (int32_t col = 0; col < row; ++col) { M(row, col) = M(col, row); }
        }
        for (int32_t row = 0; row < 6; ++row)
        {
            for (int32_t col = 0; col < 6; ++col) { M(row, col) /= rNumPoints; }
        }
        M(0, 0) = 1.0;
        return std::vector<double>(&M[0], &M[0] + 36);
    }

    // A verbatim copy of ApprQuadraticCircle2::operator()'s assembly.
    std::vector<double> BuildQuadraticCircle2M(std::vector<Vector2<double>> const& P)
    {
        int32_t const numPoints = static_cast<int32_t>(P.size());
        Matrix<4, 4, double> M{};
        for (int32_t i = 0; i < numPoints; ++i)
        {
            double x = P[i][0];
            double y = P[i][1];
            double x2 = x * x;
            double y2 = y * y;
            double xy = x * y;
            double r2 = x2 + y2;
            double xr2 = x * r2;
            double yr2 = y * r2;
            double r4 = r2 * r2;

            M(0, 1) += x;
            M(0, 2) += y;
            M(0, 3) += r2;
            M(1, 1) += x2;
            M(1, 2) += xy;
            M(1, 3) += xr2;
            M(2, 2) += y2;
            M(2, 3) += yr2;
            M(3, 3) += r4;
        }

        double const rNumPoints = static_cast<double>(numPoints);
        M(0, 0) = rNumPoints;
        for (int32_t row = 0; row < 4; ++row)
        {
            for (int32_t col = 0; col < row; ++col) { M(row, col) = M(col, row); }
        }
        for (int32_t row = 0; row < 4; ++row)
        {
            for (int32_t col = 0; col < 4; ++col) { M(row, col) /= rNumPoints; }
        }
        M(0, 0) = 1.0;
        return std::vector<double>(&M[0], &M[0] + 16);
    }

    // A verbatim copy of ApprQuadratic3::operator()'s assembly. Note that
    // upstream does NOT set M(0,0) = 1 here (the documented omission).
    std::vector<double> BuildQuadratic3M(std::vector<Vector3<double>> const& P)
    {
        int32_t const numPoints = static_cast<int32_t>(P.size());
        Matrix<10, 10, double> M{};
        for (int32_t i = 0; i < numPoints; ++i)
        {
            double x = P[i][0], y = P[i][1], z = P[i][2];
            double x2 = x * x, y2 = y * y, z2 = z * z;
            double xy = x * y, xz = x * z, yz = y * z;
            double x3 = x * x2, xy2 = x * y2, xz2 = x * z2;
            double x2y = x * xy, x2z = x * xz, xyz = x * yz;
            double y3 = y * y2, yz2 = y * z2, y2z = y * yz, z3 = z * z2;
            double x4 = x * x3, x2y2 = x * xy2, x2z2 = x * xz2;
            double x3y = x * x2y, x3z = x * x2z, x2yz = x * xyz;
            double y4 = y * y3, y2z2 = y * yz2, xy3 = x * y3;
            double xy2z = x * y2z, y3z = y * y2z, z4 = z * z3;
            double xyz2 = x * yz2, xz3 = x * z3, yz3 = y * z3;

            M(0, 1) += x;
            M(0, 2) += y;
            M(0, 3) += z;
            M(0, 4) += x2;
            M(0, 5) += xy;
            M(0, 6) += xz;
            M(0, 7) += y2;
            M(0, 8) += yz;
            M(0, 9) += z2;
            M(1, 4) += x3;
            M(1, 5) += x2y;
            M(1, 6) += x2z;
            M(1, 7) += xy2;
            M(1, 8) += xyz;
            M(1, 9) += xz2;
            M(2, 5) += xy2;
            M(2, 7) += y3;
            M(2, 8) += y2z;
            M(2, 9) += yz2;
            M(3, 9) += z3;
            M(4, 4) += x4;
            M(4, 5) += x3y;
            M(4, 6) += x3z;
            M(4, 7) += x2y2;
            M(4, 8) += x2yz;
            M(4, 9) += x2z2;
            M(5, 7) += xy3;
            M(5, 8) += xy2z;
            M(5, 9) += xyz2;
            M(6, 9) += xz3;
            M(7, 7) += y4;
            M(7, 8) += y3z;
            M(7, 9) += y2z2;
            M(8, 9) += yz3;
            M(9, 9) += z4;
        }

        double const rNumPoints = static_cast<double>(numPoints);
        M(0, 0) = rNumPoints;
        M(1, 1) = M(0, 4);
        M(1, 2) = M(0, 5);
        M(1, 3) = M(0, 6);
        M(2, 2) = M(0, 7);
        M(2, 3) = M(0, 8);
        M(2, 4) = M(1, 5);
        M(2, 6) = M(1, 8);
        M(3, 3) = M(0, 9);
        M(3, 4) = M(1, 6);
        M(3, 5) = M(1, 8);
        M(3, 6) = M(1, 9);
        M(3, 7) = M(2, 8);
        M(3, 8) = M(2, 9);
        M(5, 5) = M(4, 7);
        M(5, 6) = M(4, 8);
        M(6, 6) = M(4, 9);
        M(6, 7) = M(5, 8);
        M(6, 8) = M(5, 9);
        M(8, 8) = M(7, 9);
        for (int32_t row = 0; row < 10; ++row)
        {
            for (int32_t col = 0; col < row; ++col) { M(row, col) = M(col, row); }
        }
        for (int32_t row = 0; row < 10; ++row)
        {
            for (int32_t col = 0; col < 10; ++col) { M(row, col) /= rNumPoints; }
        }
        return std::vector<double>(&M[0], &M[0] + 100);
    }

    // A verbatim copy of ApprQuadraticSphere3::operator()'s assembly.
    std::vector<double> BuildQuadraticSphere3M(std::vector<Vector3<double>> const& P)
    {
        int32_t const numPoints = static_cast<int32_t>(P.size());
        Matrix<5, 5, double> M{};
        for (int32_t i = 0; i < numPoints; ++i)
        {
            double x = P[i][0], y = P[i][1], z = P[i][2];
            double x2 = x * x, y2 = y * y, z2 = z * z;
            double xy = x * y, xz = x * z, yz = y * z;
            double r2 = x2 + y2 + z2;
            double xr2 = x * r2, yr2 = y * r2, zr2 = z * r2;
            double r4 = r2 * r2;

            M(0, 1) += x;
            M(0, 2) += y;
            M(0, 3) += z;
            M(0, 4) += r2;
            M(1, 1) += x2;
            M(1, 2) += xy;
            M(1, 3) += xz;
            M(1, 4) += xr2;
            M(2, 2) += y2;
            M(2, 3) += yz;
            M(2, 4) += yr2;
            M(3, 3) += z2;
            M(3, 4) += zr2;
            M(4, 4) += r4;
        }

        double const rNumPoints = static_cast<double>(numPoints);
        M(0, 0) = rNumPoints;
        for (int32_t row = 0; row < 5; ++row)
        {
            for (int32_t col = 0; col < row; ++col) { M(row, col) = M(col, row); }
        }
        for (int32_t row = 0; row < 5; ++row)
        {
            for (int32_t col = 0; col < 5; ++col) { M(row, col) /= rNumPoints; }
        }
        M(0, 0) = 1.0;
        return std::vector<double>(&M[0], &M[0] + 25);
    }

    // Rejection sampling for the decoupling predicate. The loop is capped and
    // redraws the whole point set (count, mode and every coordinate) on every
    // attempt, which is everything the predicate depends on; on exhaustion it
    // falls back to the last candidate.
    using Build2 = std::vector<double>(*)(std::vector<Vector2<double>> const&);
    using Build3 = std::vector<double>(*)(std::vector<Vector3<double>> const&);
    constexpr int32_t kMaxAttempts = 24;

    Points2 DrawPoints2Decoupling(oracle::Ctx& io, int32_t minPoints,
        int32_t maxPoints, int32_t numModes, int32_t forcedMode,
        Build2 build, int32_t size, bool want)
    {
        Points2 draw;
        for (int32_t attempt = 0; attempt < kMaxAttempts; ++attempt)
        {
            draw = DrawPoints2(io, minPoints, maxPoints, numModes, forcedMode);
            if (EigenDecouples(build(draw.P), size) == want) { break; }
        }
        return draw;
    }

    Points3 DrawPoints3Decoupling(oracle::Ctx& io, int32_t minPoints,
        int32_t maxPoints, int32_t numModes, int32_t forcedMode,
        Build3 build, int32_t size, bool want)
    {
        Points3 draw;
        for (int32_t attempt = 0; attempt < kMaxAttempts; ++attempt)
        {
            draw = DrawPoints3(io, minPoints, maxPoints, numModes, forcedMode);
            if (EigenDecouples(build(draw.P), size) == want) { break; }
        }
        return draw;
    }
}

ORACLE_CASE("ApprQuadratic2.fit")
{
    auto P = RecordPoints2(io,
        DrawPoints2Decoupling(io, 1, 10, 5, -1, BuildQuadratic2M, 6, false));
    std::array<double, 6> coefficients{};
    ApprQuadratic2<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), coefficients);
    io.outReal(measure);
    for (auto value : coefficients) { io.outReal(value); }
}

// Points exactly on a lattice conic (the radius-5 circle), where the minimum
// eigenvalue is at the round-off floor and the clamp to zero is reachable.
ORACLE_CASE("ApprQuadratic2.fit.cocircular")
{
    auto P = RecordPoints2(io,
        DrawPoints2Decoupling(io, 5, 10, 5, 4, BuildQuadratic2M, 6, false));
    std::array<double, 6> coefficients{};
    ApprQuadratic2<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), coefficients);
    io.outReal(measure);
    for (auto value : coefficients) { io.outReal(value); }
}

// The degenerate Householder step of SymmetricEigensolver::Tridiagonalize
// (issue #80): upstream's Q is rebuilt with a spurious reflection, so its
// eigenvector differs from the port's.
ORACLE_CASE("ApprQuadratic2.fit.decoupledDeviation")
{
    auto P = RecordPoints2(io,
        DrawPoints2Decoupling(io, 1, 5, 5, -1, BuildQuadratic2M, 6, true));
    std::array<double, 6> coefficients{};
    ApprQuadratic2<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), coefficients);
    io.outReal(measure);
    for (auto value : coefficients) { io.outReal(value); }
}

ORACLE_CASE("ApprQuadraticCircle2.fit")
{
    auto P = RecordPoints2(io,
        DrawPoints2Decoupling(io, 1, 10, 5, -1, BuildQuadraticCircle2M, 4, false));
    Circle2<double> circle{};
    ApprQuadraticCircle2<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), circle);
    io.outReal(measure);
    io.outVec(circle.center);
    io.outReal(circle.radius);
}

ORACLE_CASE("ApprQuadraticCircle2.fit.cocircular")
{
    auto P = RecordPoints2(io,
        DrawPoints2Decoupling(io, 4, 10, 5, 4, BuildQuadraticCircle2M, 4, false));
    Circle2<double> circle{};
    ApprQuadraticCircle2<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), circle);
    io.outReal(measure);
    io.outVec(circle.center);
    io.outReal(circle.radius);
}

ORACLE_CASE("ApprQuadraticCircle2.fit.decoupledDeviation")
{
    auto P = RecordPoints2(io,
        DrawPoints2Decoupling(io, 1, 5, 5, -1, BuildQuadraticCircle2M, 4, true));
    Circle2<double> circle{};
    ApprQuadraticCircle2<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), circle);
    io.outReal(measure);
    io.outVec(circle.center);
    io.outReal(circle.radius);
}

// ----------------------------------------------------------- ApprQuadratic3

ORACLE_CASE("ApprQuadratic3.fit")
{
    auto P = RecordPoints3(io,
        DrawPoints3Decoupling(io, 1, 12, 6, -1, BuildQuadratic3M, 10, false));
    std::array<double, 10> coefficients{};
    ApprQuadratic3<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), coefficients);
    io.outReal(measure);
    for (auto value : coefficients) { io.outReal(value); }
}

ORACLE_CASE("ApprQuadratic3.fit.cospherical")
{
    auto P = RecordPoints3(io,
        DrawPoints3Decoupling(io, 8, 14, 6, 5, BuildQuadratic3M, 10, false));
    std::array<double, 10> coefficients{};
    ApprQuadratic3<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), coefficients);
    io.outReal(measure);
    for (auto value : coefficients) { io.outReal(value); }
}

ORACLE_CASE("ApprQuadratic3.fit.decoupledDeviation")
{
    auto P = RecordPoints3(io,
        DrawPoints3Decoupling(io, 1, 6, 6, -1, BuildQuadratic3M, 10, true));
    std::array<double, 10> coefficients{};
    ApprQuadratic3<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), coefficients);
    io.outReal(measure);
    for (auto value : coefficients) { io.outReal(value); }
}

ORACLE_CASE("ApprQuadraticSphere3.fit")
{
    auto P = RecordPoints3(io,
        DrawPoints3Decoupling(io, 1, 12, 6, -1, BuildQuadraticSphere3M, 5, false));
    Sphere3<double> sphere{};
    ApprQuadraticSphere3<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), sphere);
    io.outReal(measure);
    io.outVec(sphere.center);
    io.outReal(sphere.radius);
}

ORACLE_CASE("ApprQuadraticSphere3.fit.cospherical")
{
    auto P = RecordPoints3(io,
        DrawPoints3Decoupling(io, 5, 12, 6, 5, BuildQuadraticSphere3M, 5, false));
    Sphere3<double> sphere{};
    ApprQuadraticSphere3<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), sphere);
    io.outReal(measure);
    io.outVec(sphere.center);
    io.outReal(sphere.radius);
}

ORACLE_CASE("ApprQuadraticSphere3.fit.decoupledDeviation")
{
    auto P = RecordPoints3(io,
        DrawPoints3Decoupling(io, 1, 6, 6, -1, BuildQuadraticSphere3M, 5, true));
    Sphere3<double> sphere{};
    ApprQuadraticSphere3<double> fitter;
    double measure = fitter(static_cast<int32_t>(P.size()), P.data(), sphere);
    io.outReal(measure);
    io.outVec(sphere.center);
    io.outReal(sphere.radius);
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

// ------------------------------------------------------- ApprParallelLines2
//
// The port deliberately fixes four defects of upstream's Fit (issue #91,
// docs/UPSTREAM-FINDINGS.md):
//   1. ComputeF writes a30[1] = -3 instead of -3*Z12 (result-corrupting).
//   2. Fit accepts roots with sigma^2 > 1, which cannot come from a unit
//      direction (gamma^2 = 1 - sigma^2 < 0). The port accepts a root up to 1e-4
//      above 1, the accuracy of the bisected roots; a strict bound loses the
//      minimizer of nearly vertical lines.
//   3. The f1 == 0 branch uses gamma = sqrt(sigma) instead of
//      sqrt(1 - sigma^2).
//   4. Fit reads Polynomial1::operator[] past the end of its coefficient
//      vector when the arithmetic operators eliminate leading zeros. That is
//      undefined behaviour in C++ (MSVC reads heap memory in release), so
//      *every* case below rejects the draws that reach it.
// The namespace below carries a copy of upstream's Fit with the first three
// fixes switchable, plus the degree check for the fourth. The main case keeps
// the draws on which the two agree bit for bit, the deviation case those on
// which they do not.

namespace parallel
{
    struct ZValues
    {
        double Z20, Z11, Z02, Z30, Z21, Z12, Z03, Z40, Z31, Z22, Z13, Z04;

        explicit ZValues(std::vector<Vector2<double>> const& P)
            :
            Z20(0), Z11(0), Z02(0), Z30(0), Z21(0), Z12(0), Z03(0),
            Z40(0), Z31(0), Z22(0), Z13(0), Z04(0)
        {
            double const invN = 1.0 / static_cast<double>(P.size());
            for (auto const& sample : P)
            {
                double xx = sample[0] * sample[0];
                double xy = sample[0] * sample[1];
                double yy = sample[1] * sample[1];
                double xxx = xx * sample[0];
                double xxy = xy * sample[0];
                double xyy = xy * sample[1];
                double yyy = yy * sample[1];
                double xxxx = xxx * sample[0];
                double xxxy = xxx * sample[1];
                double xxyy = xx * yy;
                double xyyy = yyy * sample[0];
                double yyyy = yyy * sample[1];
                Z20 += xx; Z11 += xy; Z02 += yy;
                Z30 += xxx; Z21 += xxy; Z12 += xyy; Z03 += yyy;
                Z40 += xxxx; Z31 += xxxy; Z22 += xxyy; Z13 += xyyy; Z04 += yyyy;
            }
            Z20 *= invN; Z11 *= invN; Z02 *= invN;
            Z30 *= invN; Z21 *= invN; Z12 *= invN; Z03 *= invN;
            Z40 *= invN; Z31 *= invN; Z22 *= invN; Z13 *= invN; Z04 *= invN;
        }
    };

    void ComputeProduct(Polynomial1<double> const& A0, Polynomial1<double> const& B0,
        Polynomial1<double> const& A1, Polynomial1<double> const& B1,
        Polynomial1<double>& A2, Polynomial1<double>& B2)
    {
        Polynomial1<double> gammaSqr{ 1.0, 0.0, -1.0 };
        A2 = A0 * A1 + gammaSqr * B0 * B1;
        B2 = A0 * B1 + B0 * A1;
    }

    // 'fixA30' selects the port's coefficient (-3*Z12) over upstream's (-3).
    void ComputeF(ZValues const& data, Polynomial1<double>& f0,
        Polynomial1<double>& f1, bool fixA30)
    {
        Polynomial1<double> a11(2);
        a11[0] = data.Z11;
        a11[2] = -2.0 * data.Z11;
        Polynomial1<double> b11(1);
        b11[1] = data.Z02 - data.Z20;
        Polynomial1<double> a20(2);
        a20[0] = data.Z02;
        a20[2] = data.Z20 - data.Z02;
        Polynomial1<double> b20(1);
        b20[1] = -2.0 * data.Z11;
        Polynomial1<double> a30(3);
        a30[1] = (fixA30 ? -3.0 * data.Z12 : -3.0);
        a30[3] = 3.0 * data.Z12 - data.Z30;
        Polynomial1<double> b30(2);
        b30[0] = data.Z03;
        b30[2] = 3.0 * data.Z21 - data.Z03;
        Polynomial1<double> a21(3);
        a21[1] = data.Z03 - 2.0 * data.Z21;
        a21[3] = 3.0 * data.Z21 - data.Z03;
        Polynomial1<double> b21(2);
        b21[0] = data.Z12;
        b21[2] = data.Z30 - 3.0 * data.Z12;
        Polynomial1<double> a40(4);
        a40[0] = data.Z04;
        a40[2] = 6.0 * data.Z22 - 2.0 * data.Z04;
        a40[4] = data.Z40 - 6.0 * data.Z22 + data.Z04;
        Polynomial1<double> b40(3);
        b40[1] = -4.0 * data.Z13;
        b40[3] = 4.0 * (data.Z13 - data.Z31);
        Polynomial1<double> a31(4);
        a31[0] = data.Z13;
        a31[2] = 3.0 * data.Z31 - 5.0 * data.Z13;
        a31[4] = 4.0 * (data.Z13 - data.Z31);
        Polynomial1<double> b31(3);
        b31[1] = data.Z04 - 3.0 * data.Z22;
        b31[3] = 6.0 * data.Z22 - data.Z40 - data.Z04;

        Polynomial1<double> c0, d0, c1, d1, c2, d2, c3, d3;
        Polynomial1<double> c4, d4, c5, d5, c6, d6, c7, d7;
        ComputeProduct(a20, b20, a20, b20, c0, d0);
        ComputeProduct(a31, b31, c0, d0, c1, d1);
        ComputeProduct(a21, b21, a20, b20, c2, d2);
        ComputeProduct(a30, b30, c2, d2, c3, d3);
        ComputeProduct(a30, b30, a11, b11, c4, d4);
        ComputeProduct(a30, b30, c4, d4, c5, d5);
        ComputeProduct(c0, d0, a11, b11, c6, d6);
        ComputeProduct(a20, b20, c6, d6, c7, d7);
        f0 = 2.0 * (c1 - c7) - 3.0 * c3 + c5;
        f1 = 2.0 * (d1 - d7) - 3.0 * d3 + d5;
    }

    void UpdateParameters(ZValues const& data, double sigma, double sigmaSqr,
        double gamma, double& minSigma, double& minGamma, double& minK,
        double& minRSqr, double& minError)
    {
        double A20 = data.Z02 + (data.Z20 - data.Z02) * sigmaSqr;
        double B20 = -2.0 * data.Z11 * sigma;
        double S20 = A20 + gamma * B20;
        double A30 = -sigma * (3.0 * data.Z12 + (data.Z30 - 3.0 * data.Z12) * sigmaSqr);
        double B30 = data.Z03 + (3.0 * data.Z21 - data.Z03) * sigmaSqr;
        double S30 = A30 + gamma * B30;
        double A40 = data.Z04 + ((6.0 * data.Z22 - 2.0 * data.Z04)
            + (data.Z40 - 6.0 * data.Z22 + data.Z04) * sigmaSqr) * sigmaSqr;
        double B40 = -4.0 * sigma * (data.Z13 + (data.Z31 - data.Z13) * sigmaSqr);
        double S40 = A40 + gamma * B40;
        double k = S30 / (2.0 * S20);
        double ksqr = k * k;
        double rsqr = ksqr + S20;
        double error = S40 - 4.0 * k * S30 + (4.0 * ksqr - S20) * S20;
        if (error < minError)
        {
            minSigma = sigma;
            minGamma = gamma;
            minK = k;
            minRSqr = rsqr;
            minError = error;
        }
    }

    double Coefficient(Polynomial1<double> const& p, uint32_t i)
    {
        return i <= p.GetDegree() ? p[i] : 0.0;
    }

    // Upstream's Fit with the port's fixes 1-3 applied. Also reports through
    // 'inRange' whether upstream's own reads of f0, f1 and h stay inside
    // their coefficient vectors (fix 4).
    void FitFixed(std::vector<Vector2<double>> const& P, uint32_t maxIterations,
        Vector2<double>& C, Vector2<double>& V, double& radius, bool& inRange)
    {
        size_t const n = P.size();
        double const invN = 1.0 / static_cast<double>(n);
        std::vector<Vector2<double>> PAdjust = P;
        Vector2<double> A{ 0.0, 0.0 };
        for (auto const& sample : PAdjust) { A += sample; }
        A *= invN;
        for (auto& sample : PAdjust) { sample -= A; }

        ZValues data(PAdjust);
        Polynomial1<double> f0, f1;
        ComputeF(data, f0, f1, true);

        // The degree check uses upstream's own f0, f1 and h.
        {
            Polynomial1<double> u0, u1;
            ComputeF(data, u0, u1, false);
            Polynomial1<double> sigmaSqrPoly{ 0.0, 0.0, 1.0 };
            Polynomial1<double> u0Sqr = u0 * u0, u1Sqr = u1 * u1;
            Polynomial1<double> h = sigmaSqrPoly * u1Sqr + (u0Sqr - u1Sqr);
            bool const f1IsZero = !(u1 != Polynomial1<double>{ 0.0 });
            inRange = (u0.GetDegree() >= 8u)
                && (f1IsZero || (u1.GetDegree() >= 7u && h.GetDegree() >= 16u));
        }

        Polynomial1<double> freduced0(4), freduced1(3);
        for (uint32_t i = 0; i <= 4; ++i) { freduced0[i] = Coefficient(f0, 2 * i); }
        for (uint32_t i = 0; i <= 3; ++i) { freduced1[i] = Coefficient(f1, 2 * i + 1); }

        double minSigma = 0.0, minGamma = 1.0;
        double minK = data.Z03 / (2.0 * data.Z02);
        double minKSqr = minK * minK;
        double minRSqr = minKSqr + data.Z02;
        double minError = data.Z04 - 4.0 * minK * data.Z03
            + (4.0 * minKSqr - data.Z02) * data.Z02;

        if (f1 != Polynomial1<double>{ 0.0 })
        {
            Polynomial1<double> sigmaSqrPoly{ 0.0, 0.0, 1.0 };
            Polynomial1<double> f0Sqr = f0 * f0, f1Sqr = f1 * f1;
            Polynomial1<double> h = sigmaSqrPoly * f1Sqr + (f0Sqr - f1Sqr);
            Polynomial1<double> hreduced(8);
            for (uint32_t i = 0; i <= 8; ++i) { hreduced[i] = Coefficient(h, 2 * i); }

            std::array<double, 8> roots{};
            int32_t numRoots = RootsPolynomial<double>::Find(8, &hreduced[0],
                maxIterations, roots.data());
            for (int32_t i = 0; i < numRoots; ++i)
            {
                double sigmaSqr = roots[i];
                if (sigmaSqr > 0.0 && sigmaSqr <= 1.0 + 1e-4)
                {
                    double sigma = std::sqrt(sigmaSqr);
                    double gamma = -freduced0(sigmaSqr) / (sigma * freduced1(sigmaSqr));
                    UpdateParameters(data, sigma, sigmaSqr, gamma,
                        minSigma, minGamma, minK, minRSqr, minError);
                }
            }
        }
        else
        {
            Polynomial1<double> hreduced(4);
            for (uint32_t i = 0; i <= 4; ++i) { hreduced[i] = Coefficient(f0, 2 * i); }

            std::array<double, 4> roots{};
            int32_t numRoots = RootsPolynomial<double>::Find(4, &hreduced[0],
                maxIterations, roots.data());
            for (int32_t i = 0; i < numRoots; ++i)
            {
                double sigmaSqr = roots[i];
                if (sigmaSqr > 0.0 && sigmaSqr <= 1.0 + 1e-4)
                {
                    double sigma = std::sqrt(sigmaSqr);
                    double gamma = std::sqrt(std::max(1.0 - sigmaSqr, 0.0));
                    UpdateParameters(data, sigma, sigmaSqr, gamma,
                        minSigma, minGamma, minK, minRSqr, minError);
                    gamma = -gamma;
                    UpdateParameters(data, sigma, sigmaSqr, gamma,
                        minSigma, minGamma, minK, minRSqr, minError);
                }
            }
        }

        V = Vector2<double>{ minGamma, minSigma };
        C = A + minK * Vector2<double>{ -minSigma, minGamma };
        C -= Dot(C, V) * V;
        radius = std::sqrt(minRSqr);
    }
}

namespace parallel
{
    // Points for the parallel-lines fit, drawn but not recorded. Modes:
    //   0  uniform cloud
    //   1  exact lattice points on two parallel lattice lines
    //   2  small lattice cloud
    Points2 Draw(oracle::Ctx& io, int32_t minPoints, int32_t maxPoints)
    {
        Points2 out;
        out.n = io.rawInteger(minPoints, maxPoints);
        out.mode = io.rawInteger(0, 2);
        out.P.resize(static_cast<size_t>(out.n));
        int32_t const n = out.n;

        if (out.mode == 0)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                out.P[i][0] = io.raw(-10.0, 10.0);
                out.P[i][1] = io.raw(-10.0, 10.0);
            }
        }
        else if (out.mode == 1)
        {
            Vector2<double> center{ static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) };
            Vector2<double> dir{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            if (dir[0] == 0.0 && dir[1] == 0.0) { dir[0] = 1.0; }
            Vector2<double> normal{ -dir[1], dir[0] };
            for (int32_t i = 0; i < n; ++i)
            {
                double const s = static_cast<double>(io.rawInteger(0, 1) * 2 - 1);
                double const t = static_cast<double>(io.rawInteger(-4, 4));
                out.P[i] = Vector2<double>{ center[0] + s * normal[0] + t * dir[0],
                    center[1] + s * normal[1] + t * dir[1] };
            }
        }
        else
        {
            for (int32_t i = 0; i < n; ++i)
            {
                out.P[i][0] = static_cast<double>(io.rawInteger(-6, 6));
                out.P[i][1] = static_cast<double>(io.rawInteger(-6, 6));
            }
        }
        return out;
    }

    // Draws until upstream's own reads stay in range and the corrected fit
    // either agrees (want == false) or differs (want == true). The loop is
    // capped and redraws the point set and maxIterations, which is everything
    // the test depends on; on exhaustion it falls back to the last in-range
    // candidate, or to the last candidate at all.
    struct Draws
    {
        Points2 points;
        uint32_t maxIterations = 16;
    };

    // Upstream's a30[1] coefficient is wrong, so its h(sigma^2) is a
    // different polynomial from the corrected one on essentially every
    // input; the two fits then agree only when neither polynomial's roots
    // improve the seed (sigma,gamma) = (0,1). That happens on a few percent
    // of the draws, so the agreeing case needs a generous attempt cap.
    constexpr int32_t kMaxParallelAttempts = 400;

    Draws DrawUntil(oracle::Ctx& io, bool wantDeviation)
    {
        Draws best, fallback;
        bool haveBest = false, haveFallback = false;
        int32_t const cap = (wantDeviation ? kMaxAttempts : kMaxParallelAttempts);
        for (int32_t attempt = 0; attempt < cap; ++attempt)
        {
            Draws candidate;
            candidate.points = Draw(io, 6, 14);
            candidate.maxIterations =
                static_cast<uint32_t>(io.rawInteger(8, 32));

            Vector2<double> fixedC{}, fixedV{};
            double fixedRadius = 0.0;
            bool inRange = false;
            FitFixed(candidate.points.P, candidate.maxIterations,
                fixedC, fixedV, fixedRadius, inRange);
            if (!haveFallback) { fallback = candidate; haveFallback = true; }
            if (!inRange) { continue; }

            Vector2<double> C{}, V{};
            double radius = 0.0;
            ApprParallelLines2<double> fitter;
            fitter.Fit(candidate.points.P, candidate.maxIterations, C, V, radius);
            bool same = (C[0] == fixedC[0] && C[1] == fixedC[1]
                && V[0] == fixedV[0] && V[1] == fixedV[1]
                && radius == fixedRadius);
            if (!haveBest) { best = candidate; haveBest = true; }
            if (same != wantDeviation) { return candidate; }
        }
        return haveBest ? best : fallback;
    }
}

// The main case keeps the draws on which upstream's four defects have no
// effect on the returned (C, V, radius).
ORACLE_CASE("ApprParallelLines2.fit")
{
    auto draws = parallel::DrawUntil(io, false);
    auto P = RecordPoints2(io, draws.points);
    uint32_t maxIterations =
        static_cast<uint32_t>(io.given(static_cast<double>(draws.maxIterations)));
    Vector2<double> C{}, V{};
    double radius = 0.0;
    ApprParallelLines2<double> fitter;
    fitter.Fit(P, maxIterations, C, V, radius);
    io.outVec(C);
    io.outVec(V);
    io.outReal(radius);
}

// The deviation case keeps the draws on which they differ (issue #91).
ORACLE_CASE("ApprParallelLines2.fit.deviation")
{
    auto draws = parallel::DrawUntil(io, true);
    auto P = RecordPoints2(io, draws.points);
    uint32_t maxIterations =
        static_cast<uint32_t>(io.given(static_cast<double>(draws.maxIterations)));
    Vector2<double> C{}, V{};
    double radius = 0.0;
    ApprParallelLines2<double> fitter;
    fitter.Fit(P, maxIterations, C, V, radius);
    io.outVec(C);
    io.outVec(V);
    io.outReal(radius);
}

// -------------------------------------------------------- ApprCurveByArcs

namespace
{
    void EmitArcs(oracle::Ctx& io, std::vector<double> const& times,
        std::vector<Vector2<double>> const& points,
        std::vector<Arc2<double>> const& arcs)
    {
        io.outInt(times.size());
        for (auto t : times) { io.outReal(t); }
        for (auto const& p : points) { io.outVec(p); }
        io.outInt(arcs.size());
        for (auto const& arc : arcs)
        {
            io.outVec(arc.center);
            io.outReal(arc.radius);
            io.outVec(arc.end[0]);
            io.outVec(arc.end[1]);
        }
    }
}

// A Bezier curve with lattice control points. The arc-length subdivision
// (Romberg integration of Length(X'(t)) and a bisection for the time of a
// given length) and the circumscribed-circle solve use only + - * / and sqrt.
ORACLE_CASE("ApprCurveByArcs.compute")
{
    int32_t degree = io.integer(2, 4);
    std::vector<Vector2<double>> controls(static_cast<size_t>(degree) + 1);
    for (int32_t i = 0; i <= degree; ++i)
    {
        controls[i] = io.latticeVec<2>(-6, 6);
    }
    int32_t numArcs = io.integer(1, 3);
    auto curve = std::make_shared<BezierCurve<2, double>>(degree, controls.data());
    std::vector<double> times;
    std::vector<Vector2<double>> points;
    std::vector<Arc2<double>> arcs;
    ApproximateCurveByArcs(std::static_pointer_cast<ParametricCurve<2, double>>(curve),
        static_cast<size_t>(numArcs), times, points, arcs);
    EmitArcs(io, times, points, arcs);
}

// Collinear lattice control points, so the curve is a straight segment and
// every {P0,M,P1} triple is (numerically) collinear. With epsilon = 0 the
// test 'fabs(det) >= epsilon' is still true and upstream divides by a
// near-zero or exactly zero determinant (the preserved defect of issue
// #163); with epsilon > 0 the numeric_limits::max() line-segment sentinel is
// taken.
ORACLE_CASE("ApprCurveByArcs.compute.collinear")
{
    int32_t degree = io.integer(2, 4);
    Vector2<double> base{ static_cast<double>(io.rawInteger(-4, 4)),
        static_cast<double>(io.rawInteger(-4, 4)) };
    Vector2<double> dir{ static_cast<double>(io.rawInteger(-3, 3)),
        static_cast<double>(io.rawInteger(-3, 3)) };
    if (dir[0] == 0.0 && dir[1] == 0.0) { dir[0] = 1.0; }
    std::vector<Vector2<double>> controls(static_cast<size_t>(degree) + 1);
    for (int32_t i = 0; i <= degree; ++i)
    {
        double const t = static_cast<double>(io.rawInteger(-4, 4));
        controls[i] = io.givenVec(Vector2<double>{ base[0] + t * dir[0],
            base[1] + t * dir[1] });
    }
    int32_t numArcs = io.integer(1, 3);
    double epsilon = io.real(0.0, 0.5);
    bool useEpsilon = io.boolean();
    auto curve = std::make_shared<BezierCurve<2, double>>(degree, controls.data());
    std::vector<double> times;
    std::vector<Vector2<double>> points;
    std::vector<Arc2<double>> arcs;
    ApproximateCurveByArcs(std::static_pointer_cast<ParametricCurve<2, double>>(curve),
        static_cast<size_t>(numArcs), times, points, arcs,
        useEpsilon ? epsilon : 0.0);
    EmitArcs(io, times, points, arcs);
}

// A positive epsilon on an ordinary curve: the sentinel branch is taken
// wherever the circumscribing determinant is small.
ORACLE_CASE("ApprCurveByArcs.compute.epsilon")
{
    int32_t degree = io.integer(2, 4);
    std::vector<Vector2<double>> controls(static_cast<size_t>(degree) + 1);
    for (int32_t i = 0; i <= degree; ++i)
    {
        controls[i] = io.latticeVec<2>(-6, 6);
    }
    int32_t numArcs = io.integer(1, 4);
    double epsilon = io.real(0.0, 20.0);
    auto curve = std::make_shared<BezierCurve<2, double>>(degree, controls.data());
    std::vector<double> times;
    std::vector<Vector2<double>> points;
    std::vector<Arc2<double>> arcs;
    ApproximateCurveByArcs(std::static_pointer_cast<ParametricCurve<2, double>>(curve),
        static_cast<size_t>(numArcs), times, points, arcs, epsilon);
    EmitArcs(io, times, points, arcs);
}

// The 'Invalid input' LogAssert for numArcs == 0.
ORACLE_CASE("ApprCurveByArcs.compute.assert")
{
    int32_t degree = io.integer(2, 3);
    std::vector<Vector2<double>> controls(static_cast<size_t>(degree) + 1);
    for (int32_t i = 0; i <= degree; ++i)
    {
        controls[i] = io.latticeVec<2>(-4, 4);
    }
    int32_t numArcs = io.integer(0, 1);
    auto curve = std::make_shared<BezierCurve<2, double>>(degree, controls.data());
    std::vector<double> times;
    std::vector<Vector2<double>> points;
    std::vector<Arc2<double>> arcs;
    ApproximateCurveByArcs(std::static_pointer_cast<ParametricCurve<2, double>>(curve),
        static_cast<size_t>(numArcs), times, points, arcs);
    EmitArcs(io, times, points, arcs);
}
