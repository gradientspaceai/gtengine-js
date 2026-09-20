// Verify group 4 (approximation): differential cases for the Appr* headers
// listed in plan/verify-groups.json group 4.
//
// These are least-squares fitters and iterative minimizers. Point sets are
// kept small (6-12 points) and iteration counts tiny so that one record is
// tens of doubles and a deep run stays fast.
//
// Comparison policy per case (see the TypeScript replay for the declarations):
//
//  * ApprParabola2, ApprParaboloid3 and ApprGreatCircle3 use only + - * /
//    and sqrt (the 3x3 closed-form inverse, an LDL^T decomposition and the
//    *iterative* SymmetricEigensolver3x3, which is the sqrt-only solver, not
//    the acos/cos NISymmetricEigensolver3x3), so they are compared bit for
//    bit.
//  * ApprGreatArc3 (atan2), ApprEllipse2 / ApprEllipsoid3 (RootsPolynomial::
//    SolveCubic uses pow, atan2, sin and cos), ApprEllipseByArcs (pow),
//    ApprTorus3 (sin, cos and SolveCubic), ApprCylinder3's hemisphere and
//    mesh searches (sin, cos) and ApprCone3EllipseAndPoints (sin, cos through
//    Minimize1) reach the C math library and are compared with a tolerance.
//  * ApprCone3's minimizer path is arithmetic-only once the initial guess is
//    fixed, so the iterated cases feed it an initial cone whose cosine is
//    bit-identical in the MSVC runtime and in V8 (see AgreedAngle below) and
//    compare every output but the final acos-derived cone angle exactly.
#define ORACLE_FAMILY "v04-approximation"
#include "Oracle.h"

#include <Mathematics/ApprCone3.h>
#include <Mathematics/ApprCone3EllipseAndPoints.h>
#include <Mathematics/ApprCylinder3.h>
#include <Mathematics/ApprEllipse2.h>
#include <Mathematics/ApprEllipseByArcs.h>
#include <Mathematics/ApprEllipsoid3.h>
#include <Mathematics/ApprGreatCircle3.h>
#include <Mathematics/ApprParabola2.h>
#include <Mathematics/ApprParaboloid3.h>
#include <Mathematics/ApprTorus3.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

using namespace gte;

namespace
{
    // ---- shared generators -------------------------------------------------
    // Every helper records a fixed number of doubles per point in every mode,
    // and the point count is itself recorded first, so the TypeScript replay
    // reads a self-describing layout.

    // A dyadic rational in [-2, 2] with denominator 4: exactly representable
    // and small enough that products of lattice coordinates stay exact.
    double RawDyadic(oracle::Ctx& io)
    {
        return static_cast<double>(io.rawInteger(-8, 8)) * 0.25;
    }

    // 2D point sets for the parabola fitters.
    //   mode 0: integer lattice
    //   mode 1: uniform cloud
    //   mode 2: exactly on y = u0*x^2 + u1*x + u2 (dyadic u, lattice x)
    //   mode 3: degenerate, every point on one vertical line x = c, which
    //           makes the 3x3 system singular
    std::vector<Vector2<double>> Points2(oracle::Ctx& io, int mode, int minCount)
    {
        int n = io.integer(minCount, minCount + 4);
        std::vector<Vector2<double>> points(static_cast<size_t>(n));
        if (mode == 0)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<2>(-4, 4); }
        }
        else if (mode == 1)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.vec<2>(-5.0, 5.0); }
        }
        else if (mode == 2)
        {
            double u0 = RawDyadic(io);
            double u1 = RawDyadic(io);
            double u2 = RawDyadic(io);
            for (int i = 0; i < n; ++i)
            {
                double x = static_cast<double>(io.rawInteger(-4, 4));
                Vector2<double> p{ x, u0 * x * x + u1 * x + u2 };
                points[i] = io.givenVec(p);
            }
        }
        else
        {
            double x = static_cast<double>(io.rawInteger(-3, 3));
            for (int i = 0; i < n; ++i)
            {
                Vector2<double> p{ x, static_cast<double>(io.rawInteger(-3, 3)) };
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }

    // 3D point sets for the paraboloid fitters.
    //   mode 0: integer lattice
    //   mode 1: uniform cloud
    //   mode 2: exactly on z = u0*x^2 + u1*x*y + u2*y^2 + u3*x + u4*y + u5
    //   mode 3: degenerate, every point on one vertical line (x,y) = c
    std::vector<Vector3<double>> Points3(oracle::Ctx& io, int mode, int minCount)
    {
        int n = io.integer(minCount, minCount + 4);
        std::vector<Vector3<double>> points(static_cast<size_t>(n));
        if (mode == 0)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<3>(-3, 3); }
        }
        else if (mode == 1)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.vec<3>(-5.0, 5.0); }
        }
        else if (mode == 2)
        {
            double u0 = RawDyadic(io);
            double u1 = RawDyadic(io);
            double u2 = RawDyadic(io);
            double u3 = RawDyadic(io);
            double u4 = RawDyadic(io);
            double u5 = RawDyadic(io);
            for (int i = 0; i < n; ++i)
            {
                double x = static_cast<double>(io.rawInteger(-3, 3));
                double y = static_cast<double>(io.rawInteger(-3, 3));
                double z = u0 * x * x + u1 * x * y + u2 * y * y + u3 * x + u4 * y + u5;
                Vector3<double> p{ x, y, z };
                points[i] = io.givenVec(p);
            }
        }
        else
        {
            double x = static_cast<double>(io.rawInteger(-3, 3));
            double y = static_cast<double>(io.rawInteger(-3, 3));
            for (int i = 0; i < n; ++i)
            {
                Vector3<double> p{ x, y, static_cast<double>(io.rawInteger(-3, 3)) };
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }

    // An orthonormal right-handed frame { N, U, V } built from unrecorded
    // draws. Only vectors derived from it are recorded, so the sin/cos used
    // here never enter the compared computation.
    std::array<Vector3<double>, 3> RawFrame3(oracle::Ctx& io)
    {
        std::array<Vector3<double>, 3> basis{};
        double len = 0.0;
        do
        {
            for (int i = 0; i < 3; ++i) { basis[0][i] = io.raw(-1.0, 1.0); }
            len = Length(basis[0]);
        } while (len < 0.1 || len > 1.0);
        Normalize(basis[0]);
        ComputeOrthogonalComplement(1, basis.data());
        return basis;
    }

    // Unit-length 3D samples for the great-circle fitters.
    //   mode 0: signed coordinate axes, which give a diagonal covariance
    //           matrix with repeated eigenvalues (the degenerate branches of
    //           SymmetricEigensolver3x3)
    //   mode 1: uniform unit vectors
    //   mode 2: unit vectors near one great circle
    std::vector<Vector3<double>> UnitPoints3(oracle::Ctx& io, int mode)
    {
        int n = io.integer(3, 7);
        std::vector<Vector3<double>> points(static_cast<size_t>(n));
        if (mode == 0)
        {
            for (int i = 0; i < n; ++i)
            {
                int k = io.rawInteger(0, 5);
                Vector3<double> p{};
                p.MakeZero();
                p[k % 3] = (k < 3 ? 1.0 : -1.0);
                points[i] = io.givenVec(p);
            }
        }
        else if (mode == 1)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.unit<3>(); }
        }
        else
        {
            auto basis = RawFrame3(io);
            for (int i = 0; i < n; ++i)
            {
                double angle = io.raw(-3.1, 3.1);
                double w = io.raw(-0.2, 0.2);
                Vector3<double> p = std::cos(angle) * basis[1]
                    + std::sin(angle) * basis[2] + w * basis[0];
                Normalize(p);
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }
}

// ---- ApprParabola2 ---------------------------------------------------------

ORACLE_CASE("ApprParabola2.fit")
{
    int mode = io.index() % 4;
    auto points = Points2(io, mode, 3);
    std::array<double, 3> u{};
    double meanSquareError = 0.0;
    bool success = ApprParabola2<double>::Fit(points, u, &meanSquareError);
    io.outBool(success);
    for (int i = 0; i < 3; ++i) { io.outReal(u[i]); }
    io.outReal(meanSquareError);
}

ORACLE_CASE("ApprParabola2.fitRobust")
{
    int mode = io.index() % 4;
    auto points = Points2(io, mode, 3);
    Vector2<double> average{};
    std::array<double, 3> v{};
    double meanSquareError = 0.0;
    bool success = ApprParabola2<double>::FitRobust(points, average, v, &meanSquareError);
    io.outBool(success);
    io.outVec(average);
    for (int i = 0; i < 3; ++i) { io.outReal(v[i]); }
    io.outReal(meanSquareError);
}

// Throw parity for the "Insufficient points" assertion. The count straddles
// the LogAssert bound so that both branches occur.
ORACLE_CASE("ApprParabola2.fit.throw")
{
    int n = io.integer(0, 4);
    std::vector<Vector2<double>> points(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<2>(-3, 3); }
    std::array<double, 3> u{};
    double meanSquareError = 0.0;
    bool success = ApprParabola2<double>::Fit(points, u, &meanSquareError);
    io.outBool(success);
    for (int i = 0; i < 3; ++i) { io.outReal(u[i]); }
    io.outReal(meanSquareError);
}

ORACLE_CASE("ApprParabola2.fitRobust.throw")
{
    int n = io.integer(0, 4);
    std::vector<Vector2<double>> points(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<2>(-3, 3); }
    Vector2<double> average{};
    std::array<double, 3> v{};
    double meanSquareError = 0.0;
    bool success = ApprParabola2<double>::FitRobust(points, average, v, &meanSquareError);
    io.outBool(success);
    io.outVec(average);
    for (int i = 0; i < 3; ++i) { io.outReal(v[i]); }
    io.outReal(meanSquareError);
}

// ---- ApprParaboloid3 -------------------------------------------------------

ORACLE_CASE("ApprParaboloid3.fit")
{
    int mode = io.index() % 4;
    auto points = Points3(io, mode, 6);
    std::array<double, 6> u{};
    double meanSquareError = 0.0;
    bool success = ApprParaboloid3<double>::Fit(points, u, &meanSquareError);
    io.outBool(success);
    for (int i = 0; i < 6; ++i) { io.outReal(u[i]); }
    io.outReal(meanSquareError);
}

ORACLE_CASE("ApprParaboloid3.fitRobust")
{
    int mode = io.index() % 4;
    auto points = Points3(io, mode, 6);
    Vector3<double> average{};
    std::array<double, 6> v{};
    double meanSquareError = 0.0;
    bool success = ApprParaboloid3<double>::FitRobust(points, average, v, &meanSquareError);
    io.outBool(success);
    io.outVec(average);
    for (int i = 0; i < 6; ++i) { io.outReal(v[i]); }
    io.outReal(meanSquareError);
}

ORACLE_CASE("ApprParaboloid3.fit.throw")
{
    int n = io.integer(3, 7);
    std::vector<Vector3<double>> points(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<3>(-3, 3); }
    std::array<double, 6> u{};
    double meanSquareError = 0.0;
    bool success = ApprParaboloid3<double>::Fit(points, u, &meanSquareError);
    io.outBool(success);
    for (int i = 0; i < 6; ++i) { io.outReal(u[i]); }
    io.outReal(meanSquareError);
}

ORACLE_CASE("ApprParaboloid3.fitRobust.throw")
{
    int n = io.integer(3, 7);
    std::vector<Vector3<double>> points(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<3>(-3, 3); }
    Vector3<double> average{};
    std::array<double, 6> v{};
    double meanSquareError = 0.0;
    bool success = ApprParaboloid3<double>::FitRobust(points, average, v, &meanSquareError);
    io.outBool(success);
    io.outVec(average);
    for (int i = 0; i < 6; ++i) { io.outReal(v[i]); }
    io.outReal(meanSquareError);
}

// ---- ApprGreatCircle3 ------------------------------------------------------
// The covariance sums and the iterative SymmetricEigensolver3x3 (the QL
// solver, whose only libm call is sqrt) use arithmetic only.

ORACLE_CASE("ApprGreatCircle3.compute")
{
    int mode = io.index() % 3;
    auto points = UnitPoints3(io, mode);
    Vector3<double> normal{};
    ApprGreatCircle3<double>()(static_cast<int32_t>(points.size()),
        points.data(), normal);
    io.outVec(normal);
}

// ApprGreatArc3 calls atan2 on every projected sample and then picks the
// largest gap between consecutive angles, so its control flow depends on
// libm. The two generator modes place the samples so that the winning gap is
// separated from every other gap by more than a radian:
//   mode 0: one narrow cluster, so the wrap-around gap wins (end0 = 0);
//   mode 1: two clusters with a 2-radian hole between them, so an interior
//           gap wins and the loop's update branch is exercised.
// std::sort is not stable and Array.prototype.sort is, so the generator
// avoids exactly equal angles: the samples are distinct unit vectors drawn
// from a continuum.
ORACLE_CASE("ApprGreatArc3.compute")
{
    int mode = io.index() % 2;
    int n = io.integer(3, 7);
    auto basis = RawFrame3(io);
    std::vector<Vector3<double>> points(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
    {
        double angle = io.raw(-0.6, 0.6);
        if (mode != 0)
        {
            double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
            angle = sign * io.raw(1.0, 2.8);
        }
        double w = io.raw(-0.15, 0.15);
        Vector3<double> p = std::cos(angle) * basis[1]
            + std::sin(angle) * basis[2] + w * basis[0];
        Normalize(p);
        points[i] = io.givenVec(p);
    }
    Vector3<double> normal{}, arcEnd0{}, arcEnd1{};
    ApprGreatArc3<double>()(static_cast<int32_t>(points.size()),
        points.data(), normal, arcEnd0, arcEnd1);
    io.outVec(normal);
    io.outVec(arcEnd0);
    io.outVec(arcEnd1);
}

// ---- ApprEllipseByArcs -----------------------------------------------------

namespace
{
    void EmitArcs(oracle::Ctx& io, bool success,
        std::vector<Vector2<double>> const& points,
        std::vector<Vector2<double>> const& centers,
        std::vector<double> const& radii)
    {
        io.outBool(success);
        io.outInt(points.size());
        for (auto const& p : points) { io.outVec(p); }
        io.outInt(centers.size());
        for (auto const& c : centers) { io.outVec(c); }
        for (auto const& r : radii) { io.outReal(r); }
    }
}

// std::pow is on the path for every intermediate point, so this case is
// compared with a tolerance.
//   mode 0: small integer semi-axes; numArcs straddles the "at least 2 arcs"
//           test and a == b occurs, so both early-return branches are hit
//   mode 1: uniform semi-axes
//   mode 2: extreme aspect ratios, where (tmp - a2) * invB2mA2 cancels
ORACLE_CASE("ApprEllipseByArcs.approximate")
{
    int mode = io.index() % 3;
    double a = 0.0, b = 0.0;
    int numArcs = 0;
    if (mode == 0)
    {
        a = io.lattice(1, 5);
        b = io.lattice(1, 5);
        numArcs = io.integer(0, 6);
    }
    else if (mode == 1)
    {
        a = io.real(0.25, 8.0);
        b = io.real(0.25, 8.0);
        numArcs = io.integer(2, 6);
    }
    else
    {
        a = 1.0;
        a = io.given(a);
        b = io.real(1e-5, 1e5);
        numArcs = io.integer(2, 6);
    }
    std::vector<Vector2<double>> points{}, centers{};
    std::vector<double> radii{};
    bool success = ApproximateEllipseByArcs(a, b, numArcs, points, centers, radii);
    EmitArcs(io, success, points, centers, radii);
}

// Deliberate port fix, issue #322: the intermediate-arc loop ignores
// Circumscribe's return value, so a degenerate point triple leaves 'circle'
// holding the previous arc's centre and radius and the function still
// reports success. The port propagates the failure.
//
// The defect is reached when two consecutive interior points coincide. That
// happens when a and b are adjacent doubles: the curvature range is then
// about 3 * 2^-52 relative, so pow(ab/curv, 2/3) rounds to the same double
// for consecutive i. Ordinary semi-axes never reach it (50000 uniform and
// 50000 extreme-aspect-ratio probes of the same expression chain produced no
// degenerate triple), so the main case above stays broad.
ORACLE_CASE("ApprEllipseByArcs.approximate.deviation")
{
    int k = io.integer(1, 4);
    double a = 1.0;
    a = io.given(a);
    double b = 1.0 + static_cast<double>(k) * 0x1.0p-52;
    b = io.given(b);
    int numArcs = io.integer(3, 6);
    std::vector<Vector2<double>> points{}, centers{};
    std::vector<double> radii{};
    bool success = ApproximateEllipseByArcs(a, b, numArcs, points, centers, radii);
    EmitArcs(io, success, points, centers, radii);
}

// ---- ApprEllipse2 ----------------------------------------------------------
// The two-step gradient descent solves a cubic for the centre update, and
// RootsPolynomial::SolveCubic calls pow, atan2, sin and cos, so these cases
// are compared with a tolerance. The iteration count is fixed by the caller
// (there is no convergence test), so it cannot drift.

namespace
{
    // 2D samples for the ellipse fitters.
    //   mode 0: integer lattice, with counts small enough that the initial
    //           oriented box is degenerate (a single point or a segment gives
    //           a zero extent, and Matrix.h's operator/ by zero yields the
    //           zero matrix)
    //   mode 1: uniform cloud
    //   mode 2: samples on an ellipse, so the fit has a genuine minimum
    std::vector<Vector2<double>> EllipsePoints2(oracle::Ctx& io, int mode)
    {
        int n = io.integer(1, 8);
        std::vector<Vector2<double>> points(static_cast<size_t>(n));
        if (mode == 0)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<2>(-3, 3); }
        }
        else if (mode == 1)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.vec<2>(-4.0, 4.0); }
        }
        else
        {
            double angle0 = io.raw(-3.14, 3.14);
            double cs = std::cos(angle0);
            double sn = std::sin(angle0);
            Vector2<double> u{ cs, sn };
            Vector2<double> v{ -sn, cs };
            Vector2<double> center{ io.raw(-2.0, 2.0), io.raw(-2.0, 2.0) };
            double e0 = io.raw(0.5, 3.0);
            double e1 = io.raw(0.5, 3.0);
            for (int i = 0; i < n; ++i)
            {
                double t = io.raw(-3.14, 3.14);
                double r = io.raw(0.95, 1.05);
                Vector2<double> p = center + (r * e0 * std::cos(t)) * u
                    + (r * e1 * std::sin(t)) * v;
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }

    // A right-handed 2D frame with positive extents, used as the caller's
    // initial ellipse. GTE's Perp(x,y) = (y,-x) is the clockwise rotation, so
    // axis[1] = -Perp(axis[0]) is the right-handed choice.
    Ellipse2<double> RawEllipse2(oracle::Ctx& io)
    {
        double angle = io.raw(-3.14, 3.14);
        double cs = std::cos(angle);
        double sn = std::sin(angle);
        Ellipse2<double> ellipse{};
        Vector2<double> center{ io.raw(-2.0, 2.0), io.raw(-2.0, 2.0) };
        ellipse.center = io.givenVec(center);
        Vector2<double> axis0{ cs, sn };
        ellipse.axis[0] = io.givenVec(axis0);
        Vector2<double> axis1{ -sn, cs };
        ellipse.axis[1] = io.givenVec(axis1);
        Vector2<double> extent{ io.raw(0.5, 3.0), io.raw(0.5, 3.0) };
        ellipse.extent = io.givenVec(extent);
        return ellipse;
    }

    void EmitEllipse2(oracle::Ctx& io, double error, Ellipse2<double> const& e)
    {
        io.outReal(error);
        io.outVec(e.center);
        io.outVec(e.axis[0]);
        io.outVec(e.axis[1]);
        io.outVec(e.extent);
    }
}

ORACLE_CASE("ApprEllipse2.compute.box")
{
    int mode = io.index() % 3;
    auto points = EllipsePoints2(io, mode);
    int numIterations = io.integer(0, 2);
    Ellipse2<double> ellipse{};
    ApprEllipse2<double> fitter{};
    double error = fitter(points, static_cast<size_t>(numIterations), false, ellipse);
    EmitEllipse2(io, error, ellipse);
}

ORACLE_CASE("ApprEllipse2.compute.ellipse")
{
    int mode = io.index() % 3;
    auto points = EllipsePoints2(io, mode);
    int numIterations = io.integer(0, 2);
    Ellipse2<double> ellipse = RawEllipse2(io);
    ApprEllipse2<double> fitter{};
    double error = fitter(points, static_cast<size_t>(numIterations), true, ellipse);
    EmitEllipse2(io, error, ellipse);
}

// No deviation case exists for the deliberate port fix of issues #224/#322
// (ApprEllipse2 and ApprEllipsoid3 ignore GetContainer's failure flag). In
// the default upstream configuration GTE_APPR_QUERY_VALIDATE_INDICES is not
// defined, so ApprQuery::ValidIndices returns true unconditionally,
// ApprGaussian2/3::Fit never reports failure and GetContainer never returns
// false. The one input that would exercise the discarded flag, the empty
// point set, instead makes GetContainer dereference points[0]: it terminates
// the process with an access violation, so it cannot be recorded. See the
// group report.

// ---- ApprEllipsoid3 --------------------------------------------------------
// The 3D counterpart of ApprEllipse2, with the same SolveCubic dependency.

namespace
{
    std::vector<Vector3<double>> EllipsoidPoints3(oracle::Ctx& io, int mode)
    {
        int n = io.integer(1, 8);
        std::vector<Vector3<double>> points(static_cast<size_t>(n));
        if (mode == 0)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<3>(-3, 3); }
        }
        else if (mode == 1)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.vec<3>(-4.0, 4.0); }
        }
        else
        {
            auto basis = RawFrame3(io);
            Vector3<double> center{};
            for (int i = 0; i < 3; ++i) { center[i] = io.raw(-2.0, 2.0); }
            double e0 = io.raw(0.5, 3.0);
            double e1 = io.raw(0.5, 3.0);
            double e2 = io.raw(0.5, 3.0);
            for (int i = 0; i < n; ++i)
            {
                double theta = io.raw(-3.14, 3.14);
                double phi = io.raw(0.0, 3.14);
                double r = io.raw(0.95, 1.05);
                Vector3<double> p = center
                    + (r * e0 * std::cos(theta) * std::sin(phi)) * basis[1]
                    + (r * e1 * std::sin(theta) * std::sin(phi)) * basis[2]
                    + (r * e2 * std::cos(phi)) * basis[0];
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }

    // A right-handed 3D frame with positive extents. RawFrame3 returns the
    // right-handed { N, U, V }, and the cyclic permutation { U, V, N } is
    // right-handed too.
    Ellipsoid3<double> RawEllipsoid3(oracle::Ctx& io)
    {
        auto basis = RawFrame3(io);
        Ellipsoid3<double> ellipsoid{};
        Vector3<double> center{};
        for (int i = 0; i < 3; ++i) { center[i] = io.raw(-2.0, 2.0); }
        ellipsoid.center = io.givenVec(center);
        ellipsoid.axis[0] = io.givenVec(basis[1]);
        ellipsoid.axis[1] = io.givenVec(basis[2]);
        ellipsoid.axis[2] = io.givenVec(basis[0]);
        Vector3<double> extent{};
        for (int i = 0; i < 3; ++i) { extent[i] = io.raw(0.5, 3.0); }
        ellipsoid.extent = io.givenVec(extent);
        return ellipsoid;
    }

    void EmitEllipsoid3(oracle::Ctx& io, double error, Ellipsoid3<double> const& e)
    {
        io.outReal(error);
        io.outVec(e.center);
        io.outVec(e.axis[0]);
        io.outVec(e.axis[1]);
        io.outVec(e.axis[2]);
        io.outVec(e.extent);
    }
}

ORACLE_CASE("ApprEllipsoid3.compute.box")
{
    int mode = io.index() % 3;
    auto points = EllipsoidPoints3(io, mode);
    int numIterations = io.integer(0, 2);
    Ellipsoid3<double> ellipsoid{};
    ApprEllipsoid3<double> fitter{};
    double error = fitter(points, static_cast<size_t>(numIterations), false, ellipsoid);
    EmitEllipsoid3(io, error, ellipsoid);
}

ORACLE_CASE("ApprEllipsoid3.compute.ellipsoid")
{
    int mode = io.index() % 3;
    auto points = EllipsoidPoints3(io, mode);
    int numIterations = io.integer(0, 2);
    Ellipsoid3<double> ellipsoid = RawEllipsoid3(io);
    ApprEllipsoid3<double> fitter{};
    double error = fitter(points, static_cast<size_t>(numIterations), true, ellipsoid);
    EmitEllipsoid3(io, error, ellipsoid);
}

// ApprEllipsoid3 carries the same discarded GetContainer flag as
// ApprEllipse2 and is unreachable for the same reason; see the comment
// above.
