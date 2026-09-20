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
