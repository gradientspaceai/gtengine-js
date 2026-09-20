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
    // At least four points, so that the initial oriented box generically has
    // three positive extents and the matrix M has three distinct eigenvalues.
    // Fewer points give M a repeated eigenvalue, whose eigenvectors are not
    // determined by the data; the degenerate case below covers those.
    std::vector<Vector3<double>> EllipsoidPoints3(oracle::Ctx& io, int mode)
    {
        int n = io.integer(4, 8);
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

// One, two or three points make the initial oriented box degenerate: a zero
// extent turns Matrix.h's operator/ into the zero matrix, so M has a repeated
// eigenvalue and its eigenvectors are whatever the eigensolver happens to
// produce for that eigenspace. The axes are therefore not a function of the
// data and are not emitted; the error, the centre and the extents are.
ORACLE_CASE("ApprEllipsoid3.compute.degenerateBox")
{
    int n = io.integer(1, 3);
    std::vector<Vector3<double>> points(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<3>(-3, 3); }
    int numIterations = io.integer(0, 2);
    Ellipsoid3<double> ellipsoid{};
    ApprEllipsoid3<double> fitter{};
    double error = fitter(points, static_cast<size_t>(numIterations), false, ellipsoid);
    io.outReal(error);
    io.outVec(ellipsoid.center);
    io.outVec(ellipsoid.extent);
}

// ApprEllipsoid3 carries the same discarded GetContainer flag as
// ApprEllipse2 and is unreachable for the same reason; see the comment
// above.

// ---- ApprCylinder3 ---------------------------------------------------------

namespace
{
    // 3D samples for the cylinder fitters.
    //   mode 0: integer lattice (many coplanar and coincident configurations)
    //   mode 1: uniform cloud
    //   mode 2: on a cylinder of radius r about a random axis, with a small
    //           radial perturbation so the least-squares problem is not
    //           exactly degenerate
    std::vector<Vector3<double>> CylinderPoints3(oracle::Ctx& io, int mode,
        Vector3<double>* generatingAxis = nullptr)
    {
        int n = io.integer(6, 10);
        std::vector<Vector3<double>> points(static_cast<size_t>(n));
        if (generatingAxis != nullptr) { generatingAxis->MakeZero(); }
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
            if (generatingAxis != nullptr) { *generatingAxis = basis[0]; }
            Vector3<double> center{};
            for (int i = 0; i < 3; ++i) { center[i] = io.raw(-2.0, 2.0); }
            double radius = io.raw(0.5, 3.0);
            for (int i = 0; i < n; ++i)
            {
                double t = io.raw(-3.14, 3.14);
                double h = io.raw(-3.0, 3.0);
                double s = io.raw(0.97, 1.03);
                Vector3<double> p = center
                    + (s * radius * std::cos(t)) * basis[1]
                    + (s * radius * std::sin(t)) * basis[2]
                    + h * basis[0];
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }

    void EmitCylinder(oracle::Ctx& io, Cylinder3<double> const& cylinder)
    {
        io.outVec(cylinder.axis.origin);
        io.outVec(cylinder.axis.direction);
        io.outReal(cylinder.radius);
        io.outReal(cylinder.height);
    }
}

// The covariance eigenvector constructor. Preprocess, G and the iterative
// SymmetricEigensolver3x3 use arithmetic and sqrt only.
ORACLE_CASE("ApprCylinder3.compute.eigenIndex")
{
    int mode = io.index() % 3;
    auto points = CylinderPoints3(io, mode);
    int eigenIndex = io.integer(0, 2);
    ApprCylinder3<double> fitter(static_cast<size_t>(eigenIndex));
    Cylinder3<double> cylinder{};
    double error = fitter(points.size(), points.data(), cylinder);
    io.outReal(error);
    EmitCylinder(io, cylinder);
}

// The specified-axis constructor. The axis is a recorded unit vector, so no
// libm enters the compared computation.
ORACLE_CASE("ApprCylinder3.compute.specifiedAxis")
{
    int mode = io.index() % 3;
    Vector3<double> generatingAxis{};
    auto points = CylinderPoints3(io, mode, &generatingAxis);
    // On the on-cylinder mode the specified axis is the generating axis, so
    // the fit has to recover the generating cylinder; that is the
    // independent-reference check reported for this group.
    auto axis = (mode == 2 ? io.givenVec(generatingAxis) : io.unit<3>());
    ApprCylinder3<double> fitter(axis);
    Cylinder3<double> cylinder{};
    double error = fitter(points.size(), points.data(), cylinder);
    io.outReal(error);
    EmitCylinder(io, cylinder);
}

// Throw parity for the two LogAssert guards of the point-fitting operator:
// a zero cylinder axis and fewer than 6 points. Normalize(axis, true) leaves
// a zero vector zero, so the assertion in operator() fires.
ORACLE_CASE("ApprCylinder3.compute.throw")
{
    int n = io.integer(4, 7);
    std::vector<Vector3<double>> points(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { points[i] = io.latticeVec<3>(-3, 3); }
    Vector3<double> axis{};
    axis.MakeZero();
    if (io.index() % 3 != 0)
    {
        for (int attempt = 0; attempt < 8; ++attempt)
        {
            for (int i = 0; i < 3; ++i)
            {
                axis[i] = static_cast<double>(io.rawInteger(-1, 1));
            }
            if (axis != Vector3<double>::Zero()) { break; }
        }
    }
    axis = io.givenVec(axis);
    ApprCylinder3<double> fitter(axis);
    Cylinder3<double> cylinder{};
    double error = fitter(points.size(), points.data(), cylinder);
    io.outReal(error);
    EmitCylinder(io, cylinder);
}

// The single-threaded hemisphere search. cos and sin generate the candidate
// directions and the winner is chosen by comparing G values, so this case is
// compared with a tolerance; the sample grid is tiny (at most 13 candidate
// directions) to keep the argmin well separated and the record small.
ORACLE_CASE("ApprCylinder3.compute.hemisphere")
{
    int mode = io.index() % 3;
    auto points = CylinderPoints3(io, mode);
    // An even number of theta samples puts two candidate directions on
    // perpendicular coordinate planes, whose projected measures are exactly
    // equal for a symmetric lattice point set; the winner is then decided by
    // the last bit of cos and sin. Three to five samples avoid that tie.
    double rawTheta = 3.0 + 2.0 * static_cast<double>(io.rawInteger(0, 1));
    int numThetaSamples = static_cast<int>(io.given(rawTheta));
    int numPhiSamples = io.integer(1, 3);
    ApprCylinder3<double> fitter(0, static_cast<size_t>(numThetaSamples),
        static_cast<size_t>(numPhiSamples), true);
    Cylinder3<double> cylinder{};
    double error = fitter(points.size(), points.data(), cylinder);
    io.outReal(error);
    EmitCylinder(io, cylinder);
}

// The single-threaded mesh fit. The triangles are a fan over the recorded
// points, which the replay rebuilds from the point count.
ORACLE_CASE("ApprCylinder3.computeMesh")
{
    int mode = io.index() % 3;
    auto points = CylinderPoints3(io, mode);
    // An even number of theta samples puts two candidate directions on
    // perpendicular coordinate planes, whose projected measures are exactly
    // equal for a symmetric lattice point set; the winner is then decided by
    // the last bit of cos and sin. Three to five samples avoid that tie.
    double rawTheta = 3.0 + 2.0 * static_cast<double>(io.rawInteger(0, 1));
    int numThetaSamples = static_cast<int>(io.given(rawTheta));
    int numPhiSamples = io.integer(1, 3);
    int numTriangles = static_cast<int>(points.size()) - 2;
    std::vector<int32_t> indices(static_cast<size_t>(3 * numTriangles));
    for (int t = 0; t < numTriangles; ++t)
    {
        indices[3 * t + 0] = 0;
        indices[3 * t + 1] = t + 1;
        indices[3 * t + 2] = t + 2;
    }
    ApprCylinder3<double> fitter(0, static_cast<size_t>(numThetaSamples),
        static_cast<size_t>(numPhiSamples), false);
    Cylinder3<double> cylinder{};
    fitter(points.size(), points.data(), static_cast<size_t>(numTriangles),
        indices.data(), cylinder);
    EmitCylinder(io, cylinder);
}

// ---- ApprCone3 -------------------------------------------------------------

namespace
{
    // A cone angle whose cosine is bit-identical in the MSVC runtime and in
    // V8. The candidates are the 24 dyadic angles j/16 in [0.0625, 1.5],
    // which is inside the documented range (0, pi/2). Every k/64 for
    // k = 1..100 was compared between std::cos and Math.cos; only k = 25 and
    // k = 54 disagree, and neither is a multiple of 4, so j/16 = 4j/64 always
    // agrees. Fixing the cosine removes the only libm call that precedes the
    // Gauss-Newton and Levenberg-Marquardt iterations, which are otherwise
    // pure arithmetic, linear solves and sqrt.
    double AgreedAngle(oracle::Ctx& io)
    {
        int j = io.integer(1, 24);
        return io.given(static_cast<double>(j) / 16.0);
    }

    // 3D samples for the cone fitters.
    //   mode 0: integer lattice
    //   mode 1: uniform cloud
    //   mode 2: on a circular cone about a random axis, with a small radial
    //           perturbation
    std::vector<Vector3<double>> ConePoints3(oracle::Ctx& io, int mode)
    {
        int n = io.integer(6, 10);
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
            Vector3<double> vertex{};
            for (int i = 0; i < 3; ++i) { vertex[i] = io.raw(-2.0, 2.0); }
            double tanAngle = io.raw(0.2, 1.2);
            for (int i = 0; i < n; ++i)
            {
                double h = io.raw(0.5, 3.0);
                double t = io.raw(-3.14, 3.14);
                double s = io.raw(0.97, 1.03);
                Vector3<double> p = vertex + h * basis[0]
                    + (s * h * tanAngle * std::cos(t)) * basis[1]
                    + (s * h * tanAngle * std::sin(t)) * basis[2];
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }

    void EmitConeAndResult(oracle::Ctx& io, Vector3<double> const& vertex,
        Vector3<double> const& axis, double angle)
    {
        io.outVec(vertex);
        io.outVec(axis);
        io.outReal(angle);
    }
}

// Gauss-Newton with a caller-supplied initial cone. The tolerances are zero,
// so the convergence test can only fire on an exactly zero update and the
// iteration count is fixed by maxIterations. Every output but the final cone
// angle, which acos produces, is arithmetic-only and is compared bit for bit
// (io.outRealExact on the replay side).
ORACLE_CASE("ApprCone3.gaussNewton.initialGuess")
{
    int mode = io.index() % 3;
    auto points = ConePoints3(io, mode);
    int maxIterations = io.integer(1, 4);
    Vector3<double> coneVertex = io.vec<3>(-2.0, 2.0);
    Vector3<double> coneAxis = io.unit<3>();
    double coneAngle = AgreedAngle(io);
    ApprCone3<double> fitter{};
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        static_cast<size_t>(maxIterations), 0.0, 0.0, true,
        coneVertex, coneAxis, coneAngle);
    io.outInt(result.numIterations);
    io.outBool(result.converged);
    io.outReal(result.minError);
    io.outReal(result.minErrorDifference);
    io.outReal(result.minUpdateLength);
    for (int i = 0; i < 6; ++i) { io.outReal(result.minLocation[i]); }
    EmitConeAndResult(io, coneVertex, coneAxis, coneAngle);
}

// Levenberg-Marquardt with a caller-supplied initial cone. Every fourth
// record passes a nonpositive lambdaFactor, which upstream turns into a
// single Gauss-Newton adjustment.
//
// LevenbergMarquardtMinimizer::DoIteration builds -J^T*F from the member mF,
// which holds F at the previously *rejected* candidate whenever the inner
// lambda-adjustment loop runs DoIteration more than once for the same
// pCurrent. The port re-evaluates F at pCurrent (issue #261, "fixed"), so
// upstream and the port agree exactly on the records where no outer
// iteration repeats a DoIteration, and that is exactly the records where
// result.numAdjustments is 0 for every prefix of the iteration: DoIteration
// is repeated iff the inner loop increments numAdjustments. The probe below
// runs upstream itself with maxIterations = 1, 2, ... and keeps the longest
// prefix that stays sound, so the recorded maxIterations is an aimed
// construction rather than a rejection loop.
namespace
{
    // The number of Levenberg-Marquardt iterations of 'fitter' on these
    // inputs for which upstream never repeats a DoIteration, capped at
    // maxCap. Zero means even the first iteration repeats one.
    size_t SoundLMIterations(ApprCone3<double>& fitter,
        std::vector<Vector3<double>> const& points,
        Vector3<double> const& vertex, Vector3<double> const& axis,
        double angle, double lambdaFactor, double lambdaAdjust,
        size_t maxAdjustments, size_t maxCap, bool wantSound)
    {
        size_t chosen = (wantSound ? 0 : maxCap);
        for (size_t m = 1; m <= maxCap; ++m)
        {
            Vector3<double> v = vertex, a = axis;
            double g = angle;
            auto probe = fitter(static_cast<int32_t>(points.size()),
                points.data(), m, 0.0, 0.0, lambdaFactor, lambdaAdjust,
                maxAdjustments, true, v, a, g);
            if (wantSound)
            {
                if (probe.numAdjustments != 0) { break; }
                chosen = m;
            }
            else if (probe.numAdjustments != 0)
            {
                chosen = m;
                break;
            }
        }
        return chosen;
    }
}

ORACLE_CASE("ApprCone3.levenbergMarquardt.initialGuess")
{
    int mode = io.index() % 3;
    auto points = ConePoints3(io, mode);
    int maxAdjustments = io.integer(1, 3);
    double lambdaFactor = (io.index() % 4 == 0 ? io.real(-1.0, 0.0)
        : io.real(1e-4, 1e-2));
    double lambdaAdjust = io.real(2.0, 10.0);
    Vector3<double> coneVertex = io.vec<3>(-2.0, 2.0);
    Vector3<double> coneAxis = io.unit<3>();
    double coneAngle = AgreedAngle(io);
    ApprCone3<double> fitter{};
    size_t maxIterations = SoundLMIterations(fitter, points, coneVertex,
        coneAxis, coneAngle, lambdaFactor, lambdaAdjust,
        static_cast<size_t>(maxAdjustments), 4, true);
    io.given(static_cast<double>(maxIterations));
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        maxIterations, 0.0, 0.0,
        lambdaFactor, lambdaAdjust, static_cast<size_t>(maxAdjustments), true,
        coneVertex, coneAxis, coneAngle);
    io.outInt(result.numIterations);
    io.outInt(result.numAdjustments);
    io.outBool(result.converged);
    io.outReal(result.minError);
    io.outReal(result.minErrorDifference);
    io.outReal(result.minUpdateLength);
    for (int i = 0; i < 6; ++i) { io.outReal(result.minLocation[i]); }
    EmitConeAndResult(io, coneVertex, coneAxis, coneAngle);
}

// The deliberate port fix of issue #261, demonstrated through ApprCone3:
// maxIterations is the first prefix length at which upstream repeats a
// DoIteration and therefore builds the step from a stale residual.
ORACLE_CASE("ApprCone3.levenbergMarquardt.staleResidual.deviation")
{
    int mode = io.index() % 3;
    auto points = ConePoints3(io, mode);
    int maxAdjustments = io.integer(1, 3);
    double lambdaFactor = io.real(1e-4, 1e-2);
    double lambdaAdjust = io.real(2.0, 10.0);
    Vector3<double> coneVertex = io.vec<3>(-2.0, 2.0);
    Vector3<double> coneAxis = io.unit<3>();
    double coneAngle = AgreedAngle(io);
    ApprCone3<double> fitter{};
    size_t maxIterations = SoundLMIterations(fitter, points, coneVertex,
        coneAxis, coneAngle, lambdaFactor, lambdaAdjust,
        static_cast<size_t>(maxAdjustments), 4, false);
    io.given(static_cast<double>(maxIterations));
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        maxIterations, 0.0, 0.0,
        lambdaFactor, lambdaAdjust, static_cast<size_t>(maxAdjustments), true,
        coneVertex, coneAxis, coneAngle);
    io.outInt(result.numIterations);
    io.outInt(result.numAdjustments);
    io.outBool(result.converged);
    io.outReal(result.minError);
    io.outReal(result.minErrorDifference);
    io.outReal(result.minUpdateLength);
    for (int i = 0; i < 6; ++i) { io.outReal(result.minLocation[i]); }
    EmitConeAndResult(io, coneVertex, coneAxis, coneAngle);
}

// ComputeInitialCone, the private initial-guess helper, reached with
// maxIterations = 0 so that the minimizer returns its input untouched. Its
// atan2 and the cos/acos round trip around it are the only libm calls, and
// none of them decides a branch (the only test, hrSlope < 0, comes from
// ApprHeightLine2::Fit, which is arithmetic-only), so a tolerance is enough.
ORACLE_CASE("ApprCone3.gaussNewton.computeInitialCone")
{
    int mode = io.index() % 3;
    auto points = ConePoints3(io, mode);
    Vector3<double> coneVertex{};
    Vector3<double> coneAxis{};
    double coneAngle = 0.0;
    ApprCone3<double> fitter{};
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        0, 0.0, 0.0, false, coneVertex, coneAxis, coneAngle);
    io.outInt(result.numIterations);
    io.outBool(result.converged);
    io.outReal(result.minError);
    for (int i = 0; i < 6; ++i) { io.outReal(result.minLocation[i]); }
    EmitConeAndResult(io, coneVertex, coneAxis, coneAngle);
}

// ---- ApprTorus3 ------------------------------------------------------------
// Every evaluation of F and J calls sin and cos, the non-iterative fit goes
// through RootsPolynomial::SolveCubic, and the spherical angles of the
// initial guess come from atan2 and acos, so all of these cases are compared
// with a tolerance.

namespace
{
    // 3D samples for the torus fitters.
    //   mode 0: integer lattice
    //   mode 1: uniform cloud
    //   mode 2: on a torus about a random axis, sampled in all octants as the
    //           header requires
    std::vector<Vector3<double>> TorusPoints3(oracle::Ctx& io, int mode)
    {
        int n = io.integer(8, 12);
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
            double r0 = io.raw(1.5, 3.0);
            double r1 = io.raw(0.3, 1.0);
            for (int i = 0; i < n; ++i)
            {
                double t = io.raw(-3.14, 3.14);
                double psi = io.raw(-3.14, 3.14);
                double radial = r0 + r1 * std::cos(psi);
                Vector3<double> p = center
                    + (radial * std::cos(t)) * basis[1]
                    + (radial * std::sin(t)) * basis[2]
                    + (r1 * std::sin(psi)) * basis[0];
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }

    // The caller-supplied initial torus. Every third record uses the exact
    // axis (0,0,1) so that the |N[2]| < 1 test of the initial guess takes its
    // other branch.
    void RawTorus(oracle::Ctx& io, Vector3<double>& C, Vector3<double>& N,
        double& r0, double& r1)
    {
        Vector3<double> center{};
        for (int i = 0; i < 3; ++i) { center[i] = io.raw(-2.0, 2.0); }
        C = io.givenVec(center);
        if (io.index() % 3 == 0)
        {
            Vector3<double> axis{ 0.0, 0.0, 1.0 };
            N = io.givenVec(axis);
        }
        else
        {
            N = io.unit<3>();
        }
        r0 = io.real(1.0, 3.0);
        r1 = io.real(0.2, 0.9);
    }

    void EmitTorus(oracle::Ctx& io, Vector3<double> const& C,
        Vector3<double> const& N, double r0, double r1)
    {
        io.outVec(C);
        io.outVec(N);
        io.outReal(r0);
        io.outReal(r1);
    }

    // The torus counterpart of SoundLMIterations; see that comment.
    size_t SoundLMIterationsTorus(ApprTorus3<double> const& fitter,
        std::vector<Vector3<double>> const& points, Vector3<double> const& C,
        Vector3<double> const& N, double r0, double r1, double lambdaFactor,
        double lambdaAdjust, size_t maxAdjustments, size_t maxCap,
        bool wantSound)
    {
        size_t chosen = (wantSound ? 0 : maxCap);
        for (size_t m = 1; m <= maxCap; ++m)
        {
            Vector3<double> c = C, n = N;
            double a = r0, b = r1;
            auto probe = fitter(static_cast<int32_t>(points.size()),
                points.data(), m, 0.0, 0.0, lambdaFactor, lambdaAdjust,
                maxAdjustments, true, c, n, a, b);
            if (wantSound)
            {
                if (probe.numAdjustments != 0) { break; }
                chosen = m;
            }
            else if (probe.numAdjustments != 0)
            {
                chosen = m;
                break;
            }
        }
        return chosen;
    }
}

// The non-iterative fit: an orthogonal-plane fit followed by SolveCubic and
// an argmin of H over the positive roots.
ORACLE_CASE("ApprTorus3.compute")
{
    int mode = io.index() % 3;
    auto points = TorusPoints3(io, mode);
    Vector3<double> C{}, N{};
    double r0 = 0.0, r1 = 0.0;
    ApprTorus3<double> fitter{};
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        C, N, r0, r1);
    io.outBool(result.first);
    io.outReal(result.second);
    EmitTorus(io, C, N, r0, r1);
}

ORACLE_CASE("ApprTorus3.gaussNewton.initialGuess")
{
    int mode = io.index() % 3;
    auto points = TorusPoints3(io, mode);
    int maxIterations = io.integer(1, 2);
    Vector3<double> C{}, N{};
    double r0 = 0.0, r1 = 0.0;
    RawTorus(io, C, N, r0, r1);
    ApprTorus3<double> fitter{};
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        static_cast<size_t>(maxIterations), 0.0, 0.0, true, C, N, r0, r1);
    io.outInt(result.numIterations);
    io.outBool(result.converged);
    io.outReal(result.minError);
    io.outReal(result.minErrorDifference);
    io.outReal(result.minUpdateLength);
    for (int i = 0; i < 7; ++i) { io.outReal(result.minLocation[i]); }
    EmitTorus(io, C, N, r0, r1);
}

// Restricted to the iteration prefixes on which upstream never repeats a
// DoIteration; see the ApprCone3 Levenberg-Marquardt comment and issue #261.
ORACLE_CASE("ApprTorus3.levenbergMarquardt.initialGuess")
{
    int mode = io.index() % 3;
    auto points = TorusPoints3(io, mode);
    int maxAdjustments = io.integer(1, 3);
    double lambdaFactor = io.real(1e-4, 1e-2);
    double lambdaAdjust = io.real(2.0, 10.0);
    Vector3<double> C{}, N{};
    double r0 = 0.0, r1 = 0.0;
    RawTorus(io, C, N, r0, r1);
    ApprTorus3<double> fitter{};
    size_t maxIterations = SoundLMIterationsTorus(fitter, points, C, N, r0, r1,
        lambdaFactor, lambdaAdjust, static_cast<size_t>(maxAdjustments), 2, true);
    io.given(static_cast<double>(maxIterations));
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        maxIterations, 0.0, 0.0, lambdaFactor, lambdaAdjust,
        static_cast<size_t>(maxAdjustments), true, C, N, r0, r1);
    io.outInt(result.numIterations);
    io.outInt(result.numAdjustments);
    io.outBool(result.converged);
    io.outReal(result.minError);
    io.outReal(result.minErrorDifference);
    io.outReal(result.minUpdateLength);
    for (int i = 0; i < 7; ++i) { io.outReal(result.minLocation[i]); }
    EmitTorus(io, C, N, r0, r1);
}

// The initial-guess path of the two minimizer overloads, reached with
// maxIterations = 0 so that the minimizer returns its input untouched.
ORACLE_CASE("ApprTorus3.gaussNewton.computeInitialTorus")
{
    int mode = io.index() % 3;
    auto points = TorusPoints3(io, mode);
    Vector3<double> C{}, N{};
    double r0 = 0.0, r1 = 0.0;
    ApprTorus3<double> fitter{};
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        0, 0.0, 0.0, false, C, N, r0, r1);
    io.outInt(result.numIterations);
    io.outBool(result.converged);
    io.outReal(result.minError);
    for (int i = 0; i < 7; ++i) { io.outReal(result.minLocation[i]); }
    EmitTorus(io, C, N, r0, r1);
}

ORACLE_CASE("ApprTorus3.levenbergMarquardt.computeInitialTorus")
{
    int mode = io.index() % 3;
    auto points = TorusPoints3(io, mode);
    Vector3<double> C{}, N{};
    double r0 = 0.0, r1 = 0.0;
    ApprTorus3<double> fitter{};
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        0, 0.0, 0.0, 1e-3, 10.0, 2, false, C, N, r0, r1);
    io.outInt(result.numIterations);
    io.outInt(result.numAdjustments);
    io.outBool(result.converged);
    io.outReal(result.minError);
    for (int i = 0; i < 7; ++i) { io.outReal(result.minLocation[i]); }
    EmitTorus(io, C, N, r0, r1);
}

ORACLE_CASE("ApprCone3.levenbergMarquardt.computeInitialCone")
{
    int mode = io.index() % 3;
    auto points = ConePoints3(io, mode);
    Vector3<double> coneVertex{};
    Vector3<double> coneAxis{};
    double coneAngle = 0.0;
    ApprCone3<double> fitter{};
    auto result = fitter(static_cast<int32_t>(points.size()), points.data(),
        0, 0.0, 0.0, 1e-3, 10.0, 2, false, coneVertex, coneAxis, coneAngle);
    io.outInt(result.numIterations);
    io.outInt(result.numAdjustments);
    io.outBool(result.converged);
    io.outReal(result.minError);
    for (int i = 0; i < 6; ++i) { io.outReal(result.minLocation[i]); }
    EmitConeAndResult(io, coneVertex, coneAxis, coneAngle);
}

// ---- ApprCone3EllipseAndPoints ---------------------------------------------
// ApprCone3EllipseAndPoints::Fit minimizes a least-squares error over the
// cone angle with Minimize1, and the error function calls sin, cos and sqrt,
// so the case is compared with a tolerance. Minimize1's bracket search
// compares those libm-derived values, which is why mode 1 places the samples
// exactly on the cone that the (ellipse, theta) pair determines: the error
// function then has one deep minimum instead of a flat landscape of
// near-ties.

namespace
{
    // A right-handed Ellipse3: RawFrame3 returns { N, U, V } with N = U x V,
    // so normal = basis[0], axis[0] = basis[1] and axis[1] = basis[2].
    Ellipse3<double> RawEllipse3(oracle::Ctx& io, bool degenerateExtents)
    {
        auto basis = RawFrame3(io);
        Ellipse3<double> ellipse{};
        Vector3<double> center{};
        for (int i = 0; i < 3; ++i) { center[i] = io.raw(-2.0, 2.0); }
        ellipse.center = io.givenVec(center);
        ellipse.normal = io.givenVec(basis[0]);
        ellipse.axis[0] = io.givenVec(basis[1]);
        ellipse.axis[1] = io.givenVec(basis[2]);
        Vector2<double> extent{};
        if (degenerateExtents)
        {
            // Zero or negative extents, the inputs on which upstream's
            // ComputeCone divides by a without validating it.
            extent[0] = io.rawInteger(0, 1) == 0 ? 0.0 : io.raw(-3.0, -0.5);
            extent[1] = io.raw(0.4, 2.0);
            if (io.rawInteger(0, 1) == 0) { std::swap(extent[0], extent[1]); }
        }
        else
        {
            // a >= b > 0, the convention Ellipse3 uses for its extents.
            extent[0] = io.raw(1.0, 3.0);
            extent[1] = io.raw(0.4, 1.0);
        }
        ellipse.extent = io.givenVec(extent);
        return ellipse;
    }

    // A copy of ApprCone3EllipseAndPoints::ComputeCone, used only to place
    // the sample points on the cone whose cross-section is 'ellipse'. It is a
    // generator, not part of the compared computation.
    void ReferenceCone(double theta, double sigma0, double sigma1,
        Ellipse3<double> const& ellipse, Vector3<double>& K,
        Vector3<double>& D)
    {
        double const zero = 0.0, one = 1.0;
        auto const& C = ellipse.center;
        auto const& N = ellipse.normal;
        auto const& U = ellipse.axis[0];
        double a = ellipse.extent[0];
        double b = ellipse.extent[1];
        double bDivA = b / a;
        double eSqr = std::max(zero, one - bDivA * bDivA);
        double omesqr = one - eSqr;
        double e = std::sqrt(eSqr);
        double snTheta = std::sin(theta);
        double csTheta = std::cos(theta);
        double snPhi = sigma0 * e * csTheta;
        double snPhiSqr = snPhi * snPhi;
        double csPhi = sigma1 * std::sqrt(std::max(zero, one - snPhiSqr));
        double h = a * omesqr * csTheta / (snTheta * std::fabs(csPhi));
        D = csPhi * N + snPhi * U;
        double snThetaSqr = snTheta * snTheta;
        double csThetaSqr = csTheta * csTheta;
        Vector3<double> Q = C - ((h * snPhi * snThetaSqr) / (csThetaSqr - snPhiSqr)) * U;
        K = Q - h * D;
    }

    // Samples for the cone fit.
    //   mode 0: a uniform cloud, so the minimizer sees a generic landscape
    //   mode 1: exactly on the cone determined by (ellipse, theta), which is
    //           the configuration the header is written for
    std::vector<Vector3<double>> ConeFitPoints(oracle::Ctx& io, int mode,
        Ellipse3<double> const& ellipse)
    {
        int n = io.integer(4, 8);
        std::vector<Vector3<double>> points(static_cast<size_t>(n));
        if (mode == 0)
        {
            for (int i = 0; i < n; ++i) { points[i] = io.vec<3>(-4.0, 4.0); }
        }
        else
        {
            double theta = io.raw(0.3, 1.2);
            Vector3<double> K{}, D{};
            ReferenceCone(theta, 1.0, 1.0, ellipse, K, D);
            std::array<Vector3<double>, 3> frame{};
            frame[0] = D;
            ComputeOrthogonalComplement(1, frame.data());
            double csTheta = std::cos(theta);
            double snTheta = std::sin(theta);
            for (int i = 0; i < n; ++i)
            {
                double s = io.raw(0.5, 4.0);
                double alpha = io.raw(-3.14, 3.14);
                Vector3<double> dir = csTheta * frame[0]
                    + (snTheta * std::cos(alpha)) * frame[1]
                    + (snTheta * std::sin(alpha)) * frame[2];
                Vector3<double> p = K + s * dir;
                points[i] = io.givenVec(p);
            }
        }
        return points;
    }

    // A copy of ApprCone3ExtractEllipses::ProcessPlane.
    void ProcessPlaneRef(std::vector<Plane3<double>>& planes,
        Plane3<double> const& plane, double epsilon)
    {
        double const zero = 0.0, one = 1.0;
        double const oneMinusEpsilon = one - epsilon;
        for (size_t i = 0; i < planes.size(); ++i)
        {
            double cosAngle = Dot(plane.normal, planes[i].normal);
            double absDiff{};
            if (cosAngle > zero)
            {
                absDiff = std::fabs(plane.constant - planes[i].constant);
                if (cosAngle >= oneMinusEpsilon && absDiff <= epsilon) { return; }
            }
            else
            {
                cosAngle = -cosAngle;
                absDiff = std::fabs(plane.constant + planes[i].constant);
                if (cosAngle >= oneMinusEpsilon && absDiff <= epsilon) { return; }
            }
        }
        planes.push_back(plane);
    }

    // A copy of ApprCone3ExtractEllipses::LocatePlanes.
    void LocatePlanesRef(
        std::vector<ApprCone3ExtractEllipses<double>::OBBNode> const& nodes,
        size_t nodeIndex, double boxExtentEpsilon, double cosAngleEpsilon,
        std::vector<Plane3<double>>& planes)
    {
        auto const& node = nodes[nodeIndex];
        if (node.maxIndex >= node.minIndex + 2)
        {
            auto const& box = node.boundingVolume.box;
            for (int j = 0; j < 3; ++j)
            {
                if (box.extent[j] <= boxExtentEpsilon)
                {
                    Plane3<double> plane(box.axis[j], box.center);
                    ProcessPlaneRef(planes, plane, cosAngleEpsilon);
                    return;
                }
            }
        }
        if (node.leftChild != std::numeric_limits<size_t>::max())
        {
            LocatePlanesRef(nodes, node.leftChild, boxExtentEpsilon,
                cosAngleEpsilon, planes);
        }
        if (node.rightChild != std::numeric_limits<size_t>::max())
        {
            LocatePlanesRef(nodes, node.rightChild, boxExtentEpsilon,
                cosAngleEpsilon, planes);
        }
    }

    // The number of points each plane of ApprCone3ExtractEllipses would
    // receive, computed from upstream's own LocatePlanes, ProcessPlane and
    // AssociatePointsWithPlanes without calling Extract. Extract itself must
    // not be called unless every count is at least 3: an empty plane sends
    // ApprEllipse2 over an empty point set, whose GetContainer dereferences
    // points[0] and kills the process (issue #349).
    std::vector<size_t> PlanePointCounts(
        std::vector<ApprCone3ExtractEllipses<double>::OBBNode> const& nodes,
        std::vector<Vector3<double>> const& points, double boxExtentEpsilon,
        double cosAngleEpsilon)
    {
        std::vector<Plane3<double>> planes{};
        LocatePlanesRef(nodes, 0, std::max(boxExtentEpsilon, 0.0),
            std::max(cosAngleEpsilon, 0.0), planes);
        std::vector<size_t> counts(planes.size(), 0);
        for (auto const& point : points)
        {
            double minDistance = std::numeric_limits<double>::max();
            size_t minJ = std::numeric_limits<size_t>::max();
            for (size_t j = 0; j < planes.size(); ++j)
            {
                auto diff = point - planes[j].origin;
                double distance = std::fabs(Dot(planes[j].normal, diff));
                if (distance < minDistance)
                {
                    minDistance = distance;
                    minJ = j;
                }
            }
            if (minJ == std::numeric_limits<size_t>::max()) { return {}; }
            ++counts[minJ];
        }
        return counts;
    }

    void EmitCone3(oracle::Ctx& io, Cone3<double> const& cone)
    {
        io.outVec(cone.ray.origin);
        io.outVec(cone.ray.direction);
        io.outReal(cone.angle);
        io.outReal(cone.cosAngle);
        io.outReal(cone.sinAngle);
        io.outReal(cone.tanAngle);
        io.outReal(cone.cosAngleSqr);
        io.outReal(cone.sinAngleSqr);
        io.outReal(cone.invSinAngle);
    }
}

ORACLE_CASE("ApprCone3EllipseAndPoints.fit")
{
    int mode = io.index() % 2;
    auto ellipse = RawEllipse3(io, false);
    auto points = ConeFitPoints(io, mode, ellipse);
    ApprCone3EllipseAndPoints<double>::Control control{};
    control.maxSubdivisions = io.integer(1, 4);
    control.maxBisections = io.integer(1, 2);
    control.epsilon = io.real(1e-9, 1e-7);
    control.tolerance = io.real(1e-5, 1e-3);
    control.padding = io.real(1e-4, 1e-2);
    control.penalty = io.real(0.5, 2.0);
    auto cone = ApprCone3EllipseAndPoints<double>::Fit(ellipse, points, control);
    EmitCone3(io, cone);
}

// Deliberate port fix, issue #349: ComputeCone divides by the ellipse extent
// a without validating it and Fit divides the accumulated error by
// points.size() with no empty-set guard. The port asserts on both. Only the
// extent guard is observable: with an empty point set the error function is
// NaN, upstream's own "Failed to find fitted cone" assertion fires and the
// two sides agree on throwing. A zero or negative extent instead lets
// upstream return a degenerate cone, which the port refuses.
ORACLE_CASE("ApprCone3EllipseAndPoints.fit.deviation")
{
    int mode = io.index() % 2;
    auto ellipse = RawEllipse3(io, true);
    auto points = ConeFitPoints(io, mode, ellipse);
    ApprCone3EllipseAndPoints<double>::Control control{};
    control.maxSubdivisions = io.integer(1, 4);
    control.maxBisections = io.integer(1, 2);
    control.epsilon = io.real(1e-9, 1e-7);
    control.tolerance = io.real(1e-5, 1e-3);
    control.padding = io.real(1e-4, 1e-2);
    control.penalty = io.real(0.5, 2.0);
    auto cone = ApprCone3EllipseAndPoints<double>::Fit(ellipse, points, control);
    EmitCone3(io, cone);
}

// ApprCone3ExtractEllipses. The samples are two exact circular cross-sections
// of one cone, which is what the class is written for. Every plane must end
// up with at least three supporting points: upstream runs ApprGaussian3 and
// ApprEllipse2 over the empty index list of a point-less plane (issue #349,
// fixed in the port by discarding such planes), and a plane with one or two
// points gives an ellipse with infinite extents, whose all-NaN outputs would
// test nothing. The generator therefore picks the boxExtentEpsilon from a
// fixed ladder, keeping the first value for which upstream's own GetIndices
// reports at least three points for every plane; the choice is recorded.
ORACLE_CASE("ApprCone3ExtractEllipses.extract")
{
    // Two planar elliptical cross-sections in two different planes. The
    // semi-axes are deliberately unequal: on a circular section the
    // least-squares ellipse is rotationally symmetric, its axes are not
    // determined by the data, and the 1024 gradient-descent iterations of
    // ApprEllipse2 wander inside that symmetry.
    int perSection = io.integer(6, 8);
    std::vector<Vector3<double>> points{};
    for (int s = 0; s < 2; ++s)
    {
        auto basis = RawFrame3(io);
        Vector3<double> center{};
        for (int i = 0; i < 3; ++i) { center[i] = io.raw(-3.0, 3.0); }
        double e0 = io.raw(2.0, 4.0);
        double e1 = io.raw(0.5, 1.0);
        double phase = io.raw(-3.14, 3.14);
        for (int i = 0; i < perSection; ++i)
        {
            double a = phase + 6.283185307179586 * static_cast<double>(i)
                / static_cast<double>(perSection);
            Vector3<double> p = center
                + (e0 * std::cos(a)) * basis[1]
                + (e1 * std::sin(a)) * basis[2];
            points.push_back(io.givenVec(p));
        }
    }
    double cosAngleEpsilon = io.real(1e-4, 1e-2);
    // Extract must not be called with a boxExtentEpsilon that leaves a plane
    // without points: AssociatePointsWithPlanes then indexes
    // mIndices[size_t(-1)] when the plane set is empty, and a point-less
    // plane sends ApprEllipse2 over an empty vector whose GetContainer
    // dereferences points[0]. Both kill the process (issue #349; the port
    // returns early and discards point-less planes). PlanePointCounts
    // replays upstream's own plane location and association to pick a
    // boxExtentEpsilon that avoids those states; the fallback 1e9 makes the
    // root node flat, which always yields exactly one plane holding every
    // point.
    ApprCone3ExtractEllipses<double>::OBBTree tree{};
    tree.Create(points);
    auto const& nodes = tree.GetNodes();
    double const ladder[6] = { 1e-8, 1e-6, 1e-4, 1e-2, 1e-1, 1e9 };
    double boxExtentEpsilon = ladder[5];
    for (int k = 0; k < 6; ++k)
    {
        auto counts = PlanePointCounts(nodes, points, ladder[k], cosAngleEpsilon);
        bool ok = !counts.empty();
        for (auto count : counts)
        {
            if (count < 3) { ok = false; }
        }
        if (ok) { boxExtentEpsilon = ladder[k]; break; }
    }
    boxExtentEpsilon = io.given(boxExtentEpsilon);
    ApprCone3ExtractEllipses<double> extractor{};
    std::vector<Ellipse3<double>> ellipses{};
    extractor.Extract(points, boxExtentEpsilon, cosAngleEpsilon, ellipses);
    io.outInt(extractor.GetPlanes().size());
    for (auto const& plane : extractor.GetPlanes())
    {
        io.outVec(plane.normal);
        io.outReal(plane.constant);
    }
    io.outInt(extractor.GetIndices().size());
    for (auto const& list : extractor.GetIndices())
    {
        io.outInt(list.size());
        for (auto index : list) { io.outInt(index); }
    }
    // Only the plane normal of each extracted ellipse is emitted. The centre,
    // axes and extents come from ApprEllipse2 run for 1024 iterations of a
    // two-step gradient descent, each step solving a cubic with pow, atan2,
    // sin and cos. On 6-8 samples that iterate does not converge (it returns
    // extents an order of magnitude larger than the generating ones), so the
    // two sides wander apart: 247 of 2000 records differ, by up to 0.4
    // relative. There is no tolerance that makes those numbers a test.
    // ApprEllipse2 itself is compared directly by its own cases above, at 0
    // and 1e-15 relative error for 0 to 2 iterations.
    io.outInt(ellipses.size());
    for (auto const& e : ellipses)
    {
        io.outVec(e.normal);
    }
}

// Throw parity for the two LogAssert guards upstream keeps: an invalid
// Control and the "Failed to find fitted cone" assertion, which fires on an
// empty point set because the error function is NaN for every theta.
ORACLE_CASE("ApprCone3EllipseAndPoints.fit.throw")
{
    auto ellipse = RawEllipse3(io, false);
    int n = io.integer(0, 3);
    std::vector<Vector3<double>> points(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) { points[i] = io.vec<3>(-4.0, 4.0); }
    ApprCone3EllipseAndPoints<double>::Control control{};
    control.maxSubdivisions = io.integer(0, 3);
    control.maxBisections = io.integer(0, 3);
    control.epsilon = io.real(-1e-8, 1e-7);
    control.tolerance = io.real(1e-5, 1e-3);
    control.padding = io.real(1e-4, 1e-2);
    control.penalty = io.real(0.5, 2.0);
    auto cone = ApprCone3EllipseAndPoints<double>::Fit(ellipse, points, control);
    EmitCone3(io, cone);
}
