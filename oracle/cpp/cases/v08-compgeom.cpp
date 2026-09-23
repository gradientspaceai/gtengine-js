// Verify group 8 (computational geometry): differential cases for the headers
// listed in plan/verify-groups.json group 8.
//
// Comparison policy for this family. Almost everything here is + - * / sqrt
// and comparisons, so every ordinary case is declared exact on the TypeScript
// side. The three places where the C math library enters are handled as
// follows.
//
//   * SortPointsOnCircle::ByAngle compares std::atan2 values. A one-ulp
//     disagreement between the MSVC runtime and V8 would change the sort
//     order, which is a discrete output that no tolerance can rescue. The
//     generator therefore rejects point pairs whose directions from the sort
//     centre are distinct but (nearly) parallel: with |DotPerp(W0,W1)| >=
//     1e-6 * |W0| * |W1| the angular gap is more than ten orders of magnitude
//     above an ulp of atan2, so the comparison is decided by the geometry and
//     not by the library. Exactly coincident W vectors are allowed, because
//     atan2 of identical arguments is identical in any one library and the
//     comparator then falls through to the arithmetic tie-break.
//
//   * InscribedFixedAspectRectInQuad::Execute uses std::atan2 only to pick
//     the quadrant index j = floor((2/pi) * angle) of each edge normal. The
//     angle value itself never reaches the result. The generator uses lattice
//     quads whose four inner normals all have both components nonzero, so
//     every angle is at least atan(1/12) away from a quadrant boundary and j
//     is decided identically by both libraries. Everything after that is
//     arithmetic, so the case is exact.
//
//   * MinimumVolumeBox3 point clouds are built by rotating lattice clouds
//     with sin/cos, but only the final rotated coordinates are recorded as
//     inputs, so no libm value enters a compared computation.
//
// Container orders that C++ leaves unspecified are canonicalized identically
// on both sides and the case comment says so (BoxManager's overlap set is a
// std::set and therefore already ordered; NearestNeighborQuery's leaf site
// lists are sorted).
#define ORACLE_FAMILY "v08-compgeom"
#include "Oracle.h"

#include <Mathematics/DisjointIntervals.h>
#include <Mathematics/DisjointRectangles.h>
#include <Mathematics/SortPointsOnCircle.h>
#include <Mathematics/CircleThroughPointSpecifiedTangentAndRadius.h>
#include <Mathematics/CircleThroughTwoPointsSpecifiedRadius.h>
#include <Mathematics/ConvexHullSimplePolygon.h>
#include <Mathematics/ConvexPolyhedron3.h>
#include <Mathematics/PrimalQuery2.h>
#include <Mathematics/PrimalQuery3.h>
#include <Mathematics/ExtremalQuery3PRJ.h>
#include <Mathematics/NearestNeighborQuery.h>
#include <Mathematics/BoxManager.h>
// InscribedFixedAspectRectInQuad.h uses GTE_C_TWO_PI and GTE_C_INV_HALF_PI
// but does not include Constants.h, so it does not compile on its own. The
// include below is the workaround; the missing include is reported as an
// upstream suspect in oracle/reports/v08-compgeom.md.
#include <Mathematics/Constants.h>
#include <Mathematics/InscribedFixedAspectRectInQuad.h>
#include <Mathematics/TriangulateEC.h>
#include <Mathematics/MinimumVolumeBox3FloatingPoint.h>
#include <Mathematics/MinimumVolumeBox3Rational.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <exception>
#include <functional>
#include <limits>
#include <set>
#include <vector>

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

using namespace gte;

namespace
{
    // Both MinimumVolumeBox3 specializations put large exact-arithmetic
    // objects on the stack (the rational one uses UIntegerFP32<2561>, about
    // 10 KiB per number, and its own header asks for a 1 GiB stack reserve;
    // the floating-point one reaches the same exact types through
    // ConvexHull3). The default 1 MiB main-thread stack overflows, so those
    // cases run on a thread with a large stack reservation. Nothing else
    // changes: the callable runs to completion and its exception, if any, is
    // rethrown on the calling thread so the harness records the throw.
    // A single worker thread with a large stack, created once and reused for
    // every call. Creating one thread per record exhausts the user-mode
    // address space: a 1 GiB stack reservation multiplied by the thousands
    // of records of a deep run crashes the generator partway through.
    struct BigStackCall
    {
        std::function<void()> fn;
        std::exception_ptr error;
    };

    BigStackCall gBigStackCall{};
    HANDLE gBigStackThread = nullptr;
    HANDLE gBigStackRequest = nullptr;
    HANDLE gBigStackDone = nullptr;

    DWORD WINAPI BigStackWorker(LPVOID)
    {
        for (;;)
        {
            WaitForSingleObject(gBigStackRequest, INFINITE);
            try
            {
                gBigStackCall.fn();
            }
            catch (...)
            {
                gBigStackCall.error = std::current_exception();
            }
            SetEvent(gBigStackDone);
        }
    }

    void RunWithBigStack(std::function<void()> fn)
    {
        if (gBigStackThread == nullptr)
        {
            gBigStackRequest = CreateEventW(nullptr, FALSE, FALSE, nullptr);
            gBigStackDone = CreateEventW(nullptr, FALSE, FALSE, nullptr);
            gBigStackThread = CreateThread(nullptr, 512ull * 1024ull * 1024ull,
                BigStackWorker, nullptr, STACK_SIZE_PARAM_IS_A_RESERVATION, nullptr);
            if (gBigStackThread == nullptr)
            {
                throw std::runtime_error("CreateThread failed");
            }
        }

        gBigStackCall.fn = std::move(fn);
        gBigStackCall.error = nullptr;
        SetEvent(gBigStackRequest);
        WaitForSingleObject(gBigStackDone, INFINITE);
        if (gBigStackCall.error)
        {
            std::exception_ptr error = gBigStackCall.error;
            gBigStackCall.error = nullptr;
            std::rethrow_exception(error);
        }
    }
}

namespace
{
    // ---- DisjointIntervals / DisjointRectangles -------------------------

    // Membership in a set of half-open intervals, computed by a linear scan
    // of the reported intervals. This is independent of the merge algorithms
    // under test and is used as the reference oracle for the Boolean
    // operations.
    bool InIntervalSet(DisjointIntervals<double> const& s, double t)
    {
        int32_t const n = s.GetNumIntervals();
        for (int32_t i = 0; i < n; ++i)
        {
            double a = 0.0, b = 0.0;
            s.GetInterval(i, a, b);
            if (a <= t && t < b)
            {
                return true;
            }
        }
        return false;
    }

    // The reported intervals must be nonempty, ascending and separated.
    bool IntervalSetIsCanonical(DisjointIntervals<double> const& s)
    {
        int32_t const n = s.GetNumIntervals();
        double prevMax = -std::numeric_limits<double>::max();
        for (int32_t i = 0; i < n; ++i)
        {
            double a = 0.0, b = 0.0;
            if (!s.GetInterval(i, a, b))
            {
                return false;
            }
            if (!(a < b) || !(prevMax < a))
            {
                return false;
            }
            prevMax = b;
        }
        return true;
    }

    // Record a set: the interval count, every interval, and the out-of-range
    // probes at i = -1 and i = count (upstream returns false and writes zeros
    // there; the port returns null and the replay emits the same zeros).
    void OutIntervalSet(oracle::Ctx& io, DisjointIntervals<double> const& s)
    {
        int32_t const n = s.GetNumIntervals();
        io.outInt(n);
        for (int32_t i = 0; i < n; ++i)
        {
            double a = 0.0, b = 0.0;
            bool ok = s.GetInterval(i, a, b);
            io.outBool(ok);
            io.outReal(a);
            io.outReal(b);
        }

        double a = 0.0, b = 0.0;
        bool okHigh = s.GetInterval(n, a, b);
        io.outBool(okHigh);
        io.outReal(a);
        io.outReal(b);

        double c = 0.0, d = 0.0;
        bool okLow = s.GetInterval(-1, c, d);
        io.outBool(okLow);
        io.outReal(c);
        io.outReal(d);

        io.outBool(IntervalSetIsCanonical(s));
    }

    // The reference check for a Boolean operation on interval sets: sample
    // the real line on a half-integer lattice and compare membership in the
    // computed result with the Boolean combination of the memberships in the
    // operands. 'op' is 0 = union, 1 = intersection, 2 = difference,
    // 3 = exclusive or.
    bool IntervalOpAgreesWithReference(DisjointIntervals<double> const& s0,
        DisjointIntervals<double> const& s1,
        DisjointIntervals<double> const& result, int32_t op)
    {
        for (int32_t k = -40; k <= 40; ++k)
        {
            double t = 0.25 * static_cast<double>(k);
            bool in0 = InIntervalSet(s0, t);
            bool in1 = InIntervalSet(s1, t);
            bool expected =
                (op == 0 ? (in0 || in1) :
                (op == 1 ? (in0 && in1) :
                (op == 2 ? (in0 && !in1) : (in0 != in1))));
            if (InIntervalSet(result, t) != expected)
            {
                return false;
            }
        }
        return true;
    }

    // Draw the endpoints of one interval. Mode 0 uses a very tight lattice so
    // that touching and nested intervals and empty (xmin >= xmax) requests
    // occur constantly; mode 1 a wider lattice; mode 2 uniform reals. Every
    // mode records exactly two doubles.
    void DrawInterval(oracle::Ctx& io, int32_t mode, double& xmin, double& xmax)
    {
        if (mode == 0)
        {
            xmin = io.lattice(-3, 3);
            xmax = io.lattice(-3, 3);
        }
        else if (mode == 1)
        {
            xmin = io.lattice(-8, 8);
            xmax = io.lattice(-8, 8);
        }
        else
        {
            xmin = io.real(-8.0, 8.0);
            xmax = io.real(-8.0, 8.0);
        }
    }

    // Build an interval set from 'numOps' insertions.
    DisjointIntervals<double> MakeIntervalSet(oracle::Ctx& io, int32_t mode, int32_t numOps)
    {
        DisjointIntervals<double> s{};
        for (int32_t k = 0; k < numOps; ++k)
        {
            double xmin = 0.0, xmax = 0.0;
            DrawInterval(io, mode, xmin, xmax);
            s.Insert(xmin, xmax);
        }
        return s;
    }
}

ORACLE_CASE("DisjointIntervals.insertRemove")
{
    int32_t mode = io.integer(0, 2);
    int32_t numOps = io.integer(1, 7);
    DisjointIntervals<double> s{};
    for (int32_t k = 0; k < numOps; ++k)
    {
        int32_t op = io.integer(0, 1);
        double xmin = 0.0, xmax = 0.0;
        DrawInterval(io, mode, xmin, xmax);
        bool success = (op == 0 ? s.Insert(xmin, xmax) : s.Remove(xmin, xmax));
        io.outBool(success);
        io.outInt(s.GetNumIntervals());
    }
    OutIntervalSet(io, s);
}

ORACLE_CASE("DisjointIntervals.operators")
{
    int32_t mode = io.integer(0, 2);
    int32_t numOps0 = io.integer(1, 5);
    int32_t numOps1 = io.integer(1, 5);
    DisjointIntervals<double> s0 = MakeIntervalSet(io, mode, numOps0);
    DisjointIntervals<double> s1 = MakeIntervalSet(io, mode, numOps1);

    OutIntervalSet(io, s0);
    OutIntervalSet(io, s1);

    DisjointIntervals<double> rUnion = s0 | s1;
    OutIntervalSet(io, rUnion);
    io.outBool(IntervalOpAgreesWithReference(s0, s1, rUnion, 0));

    DisjointIntervals<double> rIntersect = s0 & s1;
    OutIntervalSet(io, rIntersect);
    io.outBool(IntervalOpAgreesWithReference(s0, s1, rIntersect, 1));

    DisjointIntervals<double> rDifference = s0 - s1;
    OutIntervalSet(io, rDifference);
    io.outBool(IntervalOpAgreesWithReference(s0, s1, rDifference, 2));

    DisjointIntervals<double> rXor = s0 ^ s1;
    OutIntervalSet(io, rXor);
    io.outBool(IntervalOpAgreesWithReference(s0, s1, rXor, 3));

    // The empty set is the identity for union / xor and the annihilator for
    // intersection; the operators are exercised with an empty operand too.
    DisjointIntervals<double> empty{};
    DisjointIntervals<double> rEmptyUnion = empty | s0;
    OutIntervalSet(io, rEmptyUnion);
    DisjointIntervals<double> rEmptyIntersect = s0 & empty;
    OutIntervalSet(io, rEmptyIntersect);
    DisjointIntervals<double> rEmptyDifference = empty - s0;
    OutIntervalSet(io, rEmptyDifference);
    DisjointIntervals<double> rEmptyXor = s0 ^ empty;
    OutIntervalSet(io, rEmptyXor);
}

namespace
{
    // ---- DisjointRectangles ---------------------------------------------

    bool InRectangleSet(DisjointRectangles<double> const& s, double x, double y)
    {
        int32_t const n = s.GetNumRectangles();
        for (int32_t i = 0; i < n; ++i)
        {
            double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
            s.GetRectangle(i, xmin, xmax, ymin, ymax);
            if (xmin <= x && x < xmax && ymin <= y && y < ymax)
            {
                return true;
            }
        }
        return false;
    }

    // The strips must be nonempty and ascending, and each strip's interval
    // set must itself be canonical.
    bool RectangleSetIsCanonical(DisjointRectangles<double> const& s)
    {
        int32_t const n = s.GetNumStrips();
        double prevMax = -std::numeric_limits<double>::max();
        int32_t total = 0;
        for (int32_t i = 0; i < n; ++i)
        {
            double ymin = 0.0, ymax = 0.0;
            DisjointIntervals<double> xset{};
            if (!s.GetStrip(i, ymin, ymax, xset))
            {
                return false;
            }
            if (!(ymin < ymax) || !(prevMax <= ymin))
            {
                return false;
            }
            if (!IntervalSetIsCanonical(xset))
            {
                return false;
            }
            prevMax = ymax;
            total += xset.GetNumIntervals();
        }
        return total == s.GetNumRectangles();
    }

    void OutRectangleSet(oracle::Ctx& io, DisjointRectangles<double> const& s)
    {
        int32_t const numRectangles = s.GetNumRectangles();
        io.outInt(numRectangles);
        for (int32_t i = 0; i < numRectangles; ++i)
        {
            double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
            bool ok = s.GetRectangle(i, xmin, xmax, ymin, ymax);
            io.outBool(ok);
            io.outReal(xmin);
            io.outReal(xmax);
            io.outReal(ymin);
            io.outReal(ymax);
        }

        // Out-of-range rectangle probe.
        double x0 = 0.0, x1 = 0.0, y0 = 0.0, y1 = 0.0;
        bool okHigh = s.GetRectangle(numRectangles, x0, x1, y0, y1);
        io.outBool(okHigh);
        io.outReal(x0);
        io.outReal(x1);
        io.outReal(y0);
        io.outReal(y1);

        int32_t const numStrips = s.GetNumStrips();
        io.outInt(numStrips);
        for (int32_t i = 0; i < numStrips; ++i)
        {
            double ymin = 0.0, ymax = 0.0;
            DisjointIntervals<double> xset{};
            bool ok = s.GetStrip(i, ymin, ymax, xset);
            io.outBool(ok);
            io.outReal(ymin);
            io.outReal(ymax);
            OutIntervalSet(io, xset);
        }

        // Out-of-range strip probes. Upstream leaves ymin, ymax and the
        // interval set untouched when it returns false, so only the bool is
        // recorded here.
        double sy0 = 0.0, sy1 = 0.0;
        DisjointIntervals<double> sxset{};
        bool okStripHigh = s.GetStrip(numStrips, sy0, sy1, sxset);
        io.outBool(okStripHigh);
        bool okStripLow = s.GetStrip(-1, sy0, sy1, sxset);
        io.outBool(okStripLow);

        io.outBool(RectangleSetIsCanonical(s));
    }

    bool RectangleOpAgreesWithReference(DisjointRectangles<double> const& s0,
        DisjointRectangles<double> const& s1,
        DisjointRectangles<double> const& result, int32_t op)
    {
        for (int32_t kx = -14; kx <= 14; ++kx)
        {
            double x = 0.5 * static_cast<double>(kx);
            for (int32_t ky = -14; ky <= 14; ++ky)
            {
                double y = 0.5 * static_cast<double>(ky);
                bool in0 = InRectangleSet(s0, x, y);
                bool in1 = InRectangleSet(s1, x, y);
                bool expected =
                    (op == 0 ? (in0 || in1) :
                    (op == 1 ? (in0 && in1) :
                    (op == 2 ? (in0 && !in1) : (in0 != in1))));
                if (InRectangleSet(result, x, y) != expected)
                {
                    return false;
                }
            }
        }
        return true;
    }

    void DrawRectangle(oracle::Ctx& io, int32_t mode, double& xmin, double& xmax,
        double& ymin, double& ymax)
    {
        DrawInterval(io, mode, xmin, xmax);
        DrawInterval(io, mode, ymin, ymax);
    }

    DisjointRectangles<double> MakeRectangleSet(oracle::Ctx& io, int32_t mode, int32_t numOps)
    {
        DisjointRectangles<double> s{};
        for (int32_t k = 0; k < numOps; ++k)
        {
            double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
            DrawRectangle(io, mode, xmin, xmax, ymin, ymax);
            s.Insert(xmin, xmax, ymin, ymax);
        }
        return s;
    }
}

ORACLE_CASE("DisjointRectangles.insertRemove")
{
    int32_t mode = io.integer(0, 2);
    int32_t numOps = io.integer(1, 6);
    DisjointRectangles<double> s{};
    for (int32_t k = 0; k < numOps; ++k)
    {
        int32_t op = io.integer(0, 1);
        double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
        DrawRectangle(io, mode, xmin, xmax, ymin, ymax);
        bool success = (op == 0 ? s.Insert(xmin, xmax, ymin, ymax)
            : s.Remove(xmin, xmax, ymin, ymax));
        io.outBool(success);
        io.outInt(s.GetNumRectangles());
    }
    OutRectangleSet(io, s);
}

ORACLE_CASE("DisjointRectangles.operators")
{
    int32_t mode = io.integer(0, 2);
    int32_t numOps0 = io.integer(1, 4);
    int32_t numOps1 = io.integer(1, 4);
    DisjointRectangles<double> s0 = MakeRectangleSet(io, mode, numOps0);
    DisjointRectangles<double> s1 = MakeRectangleSet(io, mode, numOps1);

    OutRectangleSet(io, s0);
    OutRectangleSet(io, s1);

    DisjointRectangles<double> rUnion = s0 | s1;
    OutRectangleSet(io, rUnion);
    io.outBool(RectangleOpAgreesWithReference(s0, s1, rUnion, 0));

    DisjointRectangles<double> rIntersect = s0 & s1;
    OutRectangleSet(io, rIntersect);
    io.outBool(RectangleOpAgreesWithReference(s0, s1, rIntersect, 1));

    DisjointRectangles<double> rDifference = s0 - s1;
    OutRectangleSet(io, rDifference);
    io.outBool(RectangleOpAgreesWithReference(s0, s1, rDifference, 2));

    DisjointRectangles<double> rXor = s0 ^ s1;
    OutRectangleSet(io, rXor);
    io.outBool(RectangleOpAgreesWithReference(s0, s1, rXor, 3));

    DisjointRectangles<double> empty{};
    DisjointRectangles<double> rEmptyUnion = empty | s0;
    OutRectangleSet(io, rEmptyUnion);
    DisjointRectangles<double> rEmptyIntersect = s0 & empty;
    OutRectangleSet(io, rEmptyIntersect);
    DisjointRectangles<double> rEmptyDifference = empty - s0;
    OutRectangleSet(io, rEmptyDifference);
    DisjointRectangles<double> rEmptyXor = s0 ^ empty;
    OutRectangleSet(io, rEmptyXor);
}

namespace
{
    // ---- PrimalQuery2 / PrimalQuery3 ------------------------------------
    //
    // Upstream is templated on Real and intended for an exact type; the port
    // is number-only, so the C++ side is instantiated with double. All four
    // generator modes record exactly two doubles per 2D point and three per
    // 3D point, so the replay reads the mode, the point count and then the
    // coordinates.
    //
    // Modes 1 and 2 construct exact degeneracies on the integer lattice: the
    // 12 lattice points of the circle x^2+y^2 = 25 (every 4 of them are
    // exactly cocircular and every 3 of them give a nonzero determinant that
    // is an exact integer), the 30 lattice points of the sphere
    // x^2+y^2+z^2 = 9, and collinear / coplanar runs V0 + k*dir. With
    // coordinates bounded by 10 every product and sum in the determinants is
    // an integer below 2^53, so the predicates are evaluated at exact zero
    // where the configuration is degenerate.

    int32_t const circle5[12][2] =
    {
        { 5, 0 }, { 4, 3 }, { 3, 4 }, { 0, 5 }, { -3, 4 }, { -4, 3 },
        { -5, 0 }, { -4, -3 }, { -3, -4 }, { 0, -5 }, { 3, -4 }, { 4, -3 }
    };

    int32_t const sphere3[30][3] =
    {
        { 3, 0, 0 }, { -3, 0, 0 }, { 0, 3, 0 }, { 0, -3, 0 }, { 0, 0, 3 }, { 0, 0, -3 },
        { 1, 2, 2 }, { 1, 2, -2 }, { 1, -2, 2 }, { 1, -2, -2 },
        { -1, 2, 2 }, { -1, 2, -2 }, { -1, -2, 2 }, { -1, -2, -2 },
        { 2, 1, 2 }, { 2, 1, -2 }, { 2, -1, 2 }, { 2, -1, -2 },
        { -2, 1, 2 }, { -2, 1, -2 }, { -2, -1, 2 }, { -2, -1, -2 },
        { 2, 2, 1 }, { 2, 2, -1 }, { 2, -2, 1 }, { 2, -2, -1 },
        { -2, 2, 1 }, { -2, 2, -1 }, { -2, -2, 1 }, { -2, -2, -1 }
    };

    // One 2D point in the given mode. Exactly two doubles are recorded.
    Vector2<double> MakePoint2(oracle::Ctx& io, int32_t mode,
        Vector2<double> const& base, Vector2<double> const& dir)
    {
        if (mode == 0)
        {
            return io.latticeVec<2>(-4, 4);
        }
        if (mode == 1)
        {
            int32_t k = io.rawInteger(0, 11);
            Vector2<double> p
            {
                base[0] + static_cast<double>(circle5[k][0]),
                base[1] + static_cast<double>(circle5[k][1])
            };
            return io.givenVec(p);
        }
        if (mode == 2)
        {
            if (io.rawInteger(0, 3) != 0)
            {
                double k = static_cast<double>(io.rawInteger(-3, 3));
                Vector2<double> p{ base[0] + k * dir[0], base[1] + k * dir[1] };
                return io.givenVec(p);
            }
            Vector2<double> p
            {
                static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4))
            };
            return io.givenVec(p);
        }
        return io.vec<2>(-5.0, 5.0);
    }

    std::vector<Vector2<double>> MakePoints2(oracle::Ctx& io, int32_t mode, int32_t n,
        Vector2<double>& base, Vector2<double>& dir)
    {
        base[0] = static_cast<double>(io.rawInteger(-3, 3));
        base[1] = static_cast<double>(io.rawInteger(-3, 3));
        dir = Vector2<double>{ 0.0, 0.0 };
        do
        {
            dir[0] = static_cast<double>(io.rawInteger(-2, 2));
            dir[1] = static_cast<double>(io.rawInteger(-2, 2));
        }
        while (dir[0] == 0.0 && dir[1] == 0.0);

        std::vector<Vector2<double>> P(static_cast<size_t>(n));
        for (int32_t i = 0; i < n; ++i)
        {
            P[static_cast<size_t>(i)] = MakePoint2(io, mode, base, dir);
        }
        return P;
    }

    // One 3D point in the given mode. Exactly three doubles are recorded.
    Vector3<double> MakePoint3(oracle::Ctx& io, int32_t mode,
        Vector3<double> const& base, Vector3<double> const& dir0,
        Vector3<double> const& dir1)
    {
        if (mode == 0)
        {
            return io.latticeVec<3>(-3, 3);
        }
        if (mode == 1)
        {
            int32_t k = io.rawInteger(0, 29);
            Vector3<double> p
            {
                base[0] + static_cast<double>(sphere3[k][0]),
                base[1] + static_cast<double>(sphere3[k][1]),
                base[2] + static_cast<double>(sphere3[k][2])
            };
            return io.givenVec(p);
        }
        if (mode == 2)
        {
            if (io.rawInteger(0, 3) != 0)
            {
                // Coplanar (and, when the second coefficient is zero,
                // collinear) with the frame <base, dir0, dir1>.
                double a = static_cast<double>(io.rawInteger(-2, 2));
                double b = static_cast<double>(io.rawInteger(-2, 2));
                Vector3<double> p = base + a * dir0 + b * dir1;
                return io.givenVec(p);
            }
            Vector3<double> p
            {
                static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3))
            };
            return io.givenVec(p);
        }
        return io.vec<3>(-4.0, 4.0);
    }

    std::vector<Vector3<double>> MakePoints3(oracle::Ctx& io, int32_t mode, int32_t n,
        Vector3<double>& base, Vector3<double>& dir0, Vector3<double>& dir1)
    {
        for (int32_t j = 0; j < 3; ++j)
        {
            base[j] = static_cast<double>(io.rawInteger(-2, 2));
        }
        dir0 = Vector3<double>{ 0.0, 0.0, 0.0 };
        dir1 = Vector3<double>{ 0.0, 0.0, 0.0 };
        do
        {
            for (int32_t j = 0; j < 3; ++j)
            {
                dir0[j] = static_cast<double>(io.rawInteger(-2, 2));
                dir1[j] = static_cast<double>(io.rawInteger(-2, 2));
            }
        }
        while (Length(Cross(dir0, dir1)) == 0.0);

        std::vector<Vector3<double>> P(static_cast<size_t>(n));
        for (int32_t i = 0; i < n; ++i)
        {
            P[static_cast<size_t>(i)] = MakePoint3(io, mode, base, dir0, dir1);
        }
        return P;
    }
}

ORACLE_CASE("PrimalQuery2.toLine")
{
    int32_t mode = io.integer(0, 3);
    int32_t n = io.integer(4, 9);
    Vector2<double> base{ 0.0, 0.0 }, dir{ 0.0, 0.0 };
    std::vector<Vector2<double>> P = MakePoints2(io, mode, n, base, dir);
    PrimalQuery2<double> query(n, P.data());
    io.outInt(query.GetNumVertices());

    for (int32_t q = 0; q < 4; ++q)
    {
        int32_t i = io.integer(0, n - 1);
        int32_t v0 = io.integer(0, n - 1);
        int32_t v1 = io.integer(0, n - 1);
        io.outInt(query.ToLine(i, v0, v1));

        Vector2<double> test = MakePoint2(io, mode, base, dir);
        io.outInt(query.ToLine(test, v0, v1));
    }
}

ORACLE_CASE("PrimalQuery2.toLineWithOrder")
{
    // The four-argument overload's collinear 'order' has a known upstream
    // defect (issue #100: it squares P-V0 instead of V1-V0) which the port
    // preserves verbatim, so the values are compared bit for bit. Modes 1 and
    // 2 make the collinear branch common.
    int32_t mode = io.integer(0, 3);
    int32_t n = io.integer(4, 9);
    Vector2<double> base{ 0.0, 0.0 }, dir{ 0.0, 0.0 };
    std::vector<Vector2<double>> P = MakePoints2(io, mode, n, base, dir);
    PrimalQuery2<double> query(n, P.data());

    for (int32_t q = 0; q < 4; ++q)
    {
        int32_t i = io.integer(0, n - 1);
        int32_t v0 = io.integer(0, n - 1);
        int32_t v1 = io.integer(0, n - 1);

        int32_t orderI = 0;
        int32_t signI = query.ToLine(i, v0, v1, orderI);
        io.outInt(signI);
        io.outInt(orderI);

        Vector2<double> test = MakePoint2(io, mode, base, dir);
        int32_t orderT = 0;
        int32_t signT = query.ToLine(test, v0, v1, orderT);
        io.outInt(signT);
        io.outInt(orderT);
    }
}

ORACLE_CASE("PrimalQuery2.toTriangle")
{
    int32_t mode = io.integer(0, 3);
    int32_t n = io.integer(4, 9);
    Vector2<double> base{ 0.0, 0.0 }, dir{ 0.0, 0.0 };
    std::vector<Vector2<double>> P = MakePoints2(io, mode, n, base, dir);
    PrimalQuery2<double> query(n, P.data());

    for (int32_t q = 0; q < 4; ++q)
    {
        int32_t i = io.integer(0, n - 1);
        int32_t v0 = io.integer(0, n - 1);
        int32_t v1 = io.integer(0, n - 1);
        int32_t v2 = io.integer(0, n - 1);
        io.outInt(query.ToTriangle(i, v0, v1, v2));

        Vector2<double> test = MakePoint2(io, mode, base, dir);
        io.outInt(query.ToTriangle(test, v0, v1, v2));
    }
}

ORACLE_CASE("PrimalQuery2.toCircumcircle")
{
    int32_t mode = io.integer(0, 3);
    int32_t n = io.integer(4, 9);
    Vector2<double> base{ 0.0, 0.0 }, dir{ 0.0, 0.0 };
    std::vector<Vector2<double>> P = MakePoints2(io, mode, n, base, dir);
    PrimalQuery2<double> query(n, P.data());

    for (int32_t q = 0; q < 4; ++q)
    {
        int32_t i = io.integer(0, n - 1);
        int32_t v0 = io.integer(0, n - 1);
        int32_t v1 = io.integer(0, n - 1);
        int32_t v2 = io.integer(0, n - 1);
        io.outInt(query.ToCircumcircle(i, v0, v1, v2));

        Vector2<double> test = MakePoint2(io, mode, base, dir);
        io.outInt(query.ToCircumcircle(test, v0, v1, v2));
    }
}

ORACLE_CASE("PrimalQuery2.toLineExtended")
{
    int32_t mode = io.integer(0, 3);
    int32_t n = io.integer(4, 9);
    Vector2<double> base{ 0.0, 0.0 }, dir{ 0.0, 0.0 };
    std::vector<Vector2<double>> P = MakePoints2(io, mode, n, base, dir);
    PrimalQuery2<double> query(n, P.data());

    for (int32_t q = 0; q < 5; ++q)
    {
        int32_t i = io.integer(0, n - 1);
        int32_t v0 = io.integer(0, n - 1);
        int32_t v1 = io.integer(0, n - 1);
        auto order = query.ToLineExtended(P[static_cast<size_t>(i)],
            P[static_cast<size_t>(v0)], P[static_cast<size_t>(v1)]);
        io.outInt(static_cast<int32_t>(order));
    }
}

ORACLE_CASE("PrimalQuery3.toPlane")
{
    int32_t mode = io.integer(0, 3);
    int32_t n = io.integer(5, 10);
    Vector3<double> base{ 0.0, 0.0, 0.0 }, dir0{ 0.0, 0.0, 0.0 }, dir1{ 0.0, 0.0, 0.0 };
    std::vector<Vector3<double>> P = MakePoints3(io, mode, n, base, dir0, dir1);
    PrimalQuery3<double> query(n, P.data());
    io.outInt(query.GetNumVertices());

    for (int32_t q = 0; q < 4; ++q)
    {
        int32_t i = io.integer(0, n - 1);
        int32_t v0 = io.integer(0, n - 1);
        int32_t v1 = io.integer(0, n - 1);
        int32_t v2 = io.integer(0, n - 1);
        io.outInt(query.ToPlane(i, v0, v1, v2));

        Vector3<double> test = MakePoint3(io, mode, base, dir0, dir1);
        io.outInt(query.ToPlane(test, v0, v1, v2));
    }
}

ORACLE_CASE("PrimalQuery3.toTetrahedron")
{
    int32_t mode = io.integer(0, 3);
    int32_t n = io.integer(5, 10);
    Vector3<double> base{ 0.0, 0.0, 0.0 }, dir0{ 0.0, 0.0, 0.0 }, dir1{ 0.0, 0.0, 0.0 };
    std::vector<Vector3<double>> P = MakePoints3(io, mode, n, base, dir0, dir1);
    PrimalQuery3<double> query(n, P.data());

    for (int32_t q = 0; q < 4; ++q)
    {
        int32_t i = io.integer(0, n - 1);
        int32_t v0 = io.integer(0, n - 1);
        int32_t v1 = io.integer(0, n - 1);
        int32_t v2 = io.integer(0, n - 1);
        int32_t v3 = io.integer(0, n - 1);
        io.outInt(query.ToTetrahedron(i, v0, v1, v2, v3));

        Vector3<double> test = MakePoint3(io, mode, base, dir0, dir1);
        io.outInt(query.ToTetrahedron(test, v0, v1, v2, v3));
    }
}

ORACLE_CASE("PrimalQuery3.toCircumsphere")
{
    int32_t mode = io.integer(0, 3);
    int32_t n = io.integer(5, 10);
    Vector3<double> base{ 0.0, 0.0, 0.0 }, dir0{ 0.0, 0.0, 0.0 }, dir1{ 0.0, 0.0, 0.0 };
    std::vector<Vector3<double>> P = MakePoints3(io, mode, n, base, dir0, dir1);
    PrimalQuery3<double> query(n, P.data());

    for (int32_t q = 0; q < 4; ++q)
    {
        int32_t i = io.integer(0, n - 1);
        int32_t v0 = io.integer(0, n - 1);
        int32_t v1 = io.integer(0, n - 1);
        int32_t v2 = io.integer(0, n - 1);
        int32_t v3 = io.integer(0, n - 1);
        io.outInt(query.ToCircumsphere(i, v0, v1, v2, v3));

        Vector3<double> test = MakePoint3(io, mode, base, dir0, dir1);
        io.outInt(query.ToCircumsphere(test, v0, v1, v2, v3));
    }
}

namespace
{
    // ---- SortPointsOnCircle ---------------------------------------------

    struct SortInput
    {
        std::vector<std::array<double, 2>> P;
        std::array<double, 2> C;
        std::array<double, 2> D;
        bool sortCCW;
    };

    // The W vectors the sort operates on, computed exactly as upstream does.
    std::vector<std::array<double, 2>> SortW(SortInput const& in)
    {
        std::array<double, 2> Dperp = (in.sortCCW
            ? std::array<double, 2>{ -in.D[1], in.D[0] }
            : std::array<double, 2>{ in.D[1], -in.D[0] });
        std::vector<std::array<double, 2>> W(in.P.size());
        for (size_t i = 0; i < in.P.size(); ++i)
        {
            std::array<double, 2> V = { in.P[i][0] - in.C[0], in.P[i][1] - in.C[1] };
            W[i] = { in.D[0] * V[0] + in.D[1] * V[1], Dperp[0] * V[0] + Dperp[1] * V[1] };
        }
        return W;
    }

    // No point may coincide with the sort centre: upstream's
    // LessThanByGeometry is then not a strict weak ordering, which is
    // undefined behaviour for std::sort (finding #394). Such inputs are kept
    // out of every generator here.
    bool NoZeroW(std::vector<std::array<double, 2>> const& W)
    {
        for (auto const& w : W)
        {
            if (w[0] == 0.0 && w[1] == 0.0)
            {
                return false;
            }
        }
        return true;
    }

    // For the ByAngle comparison, two W vectors must either be identical (so
    // that std::atan2 returns the same double for both in any one library) or
    // have an angular gap far above an ulp of atan2.
    bool AnglesAreWellSeparated(std::vector<std::array<double, 2>> const& W)
    {
        for (size_t i = 0; i + 1 < W.size(); ++i)
        {
            for (size_t j = i + 1; j < W.size(); ++j)
            {
                if (W[i][0] == W[j][0] && W[i][1] == W[j][1])
                {
                    continue;
                }
                double cross = W[i][0] * W[j][1] - W[j][0] * W[i][1];
                double li = std::sqrt(W[i][0] * W[i][0] + W[i][1] * W[i][1]);
                double lj = std::sqrt(W[j][0] * W[j][0] + W[j][1] * W[j][1]);
                if (std::fabs(cross) < 1e-6 * li * lj)
                {
                    return false;
                }
            }
        }
        return true;
    }

    // Draw a candidate sort input without recording anything. Mode 0 is a
    // small lattice, mode 1 is uniform and mode 2 repeats a few of the drawn
    // points verbatim so that the comparators see exactly equivalent
    // elements.
    SortInput DrawSortInput(oracle::Ctx& io, int32_t mode, int32_t n)
    {
        SortInput in{};
        in.P.resize(static_cast<size_t>(n));
        if (mode == 1)
        {
            in.C = { io.raw(-3.0, 3.0), io.raw(-3.0, 3.0) };
            do
            {
                in.D = { io.raw(-2.0, 2.0), io.raw(-2.0, 2.0) };
            }
            while (in.D[0] == 0.0 && in.D[1] == 0.0);
            for (int32_t i = 0; i < n; ++i)
            {
                in.P[static_cast<size_t>(i)] = { io.raw(-6.0, 6.0), io.raw(-6.0, 6.0) };
            }
        }
        else
        {
            in.C = { static_cast<double>(io.rawInteger(-2, 2)),
                static_cast<double>(io.rawInteger(-2, 2)) };
            do
            {
                in.D = { static_cast<double>(io.rawInteger(-2, 2)),
                    static_cast<double>(io.rawInteger(-2, 2)) };
            }
            while (in.D[0] == 0.0 && in.D[1] == 0.0);
            for (int32_t i = 0; i < n; ++i)
            {
                if (mode == 2 && i > 0 && io.rawInteger(0, 2) == 0)
                {
                    in.P[static_cast<size_t>(i)] =
                        in.P[static_cast<size_t>(io.rawInteger(0, i - 1))];
                }
                else
                {
                    in.P[static_cast<size_t>(i)] =
                        { static_cast<double>(io.rawInteger(-5, 5)),
                          static_cast<double>(io.rawInteger(-5, 5)) };
                }
            }
        }
        in.sortCCW = (io.rawInteger(0, 1) != 0);
        return in;
    }

    // A construction that always satisfies both acceptance tests: distinct
    // directions taken from the 12 lattice points of the circle of radius 5.
    SortInput FallbackSortInput(oracle::Ctx& io, int32_t n)
    {
        SortInput in{};
        in.C = { 0.0, 0.0 };
        in.D = { 1.0, 0.0 };
        in.sortCCW = (io.rawInteger(0, 1) != 0);
        in.P.resize(static_cast<size_t>(n));
        for (int32_t i = 0; i < n; ++i)
        {
            int32_t k = i % 12;
            double scale = static_cast<double>(1 + (i / 12));
            in.P[static_cast<size_t>(i)] =
                { scale * static_cast<double>(circle5[k][0]),
                  scale * static_cast<double>(circle5[k][1]) };
        }
        return in;
    }

    // Record a sort input: the points, the centre, the direction and the
    // ordering flag.
    void GiveSortInput(oracle::Ctx& io, SortInput const& in)
    {
        for (auto const& p : in.P)
        {
            io.given(p[0]);
            io.given(p[1]);
        }
        io.given(in.C[0]);
        io.given(in.C[1]);
        io.given(in.D[0]);
        io.given(in.D[1]);
        io.given(in.sortCCW ? 1.0 : 0.0);
    }
}

ORACLE_CASE("SortPointsOnCircle.byAngleAndByGeometry")
{
    // Both sorts of the same input. The generator rejects points coinciding
    // with the centre (upstream's geometric comparator is then not a strict
    // weak ordering, finding #394) and directions that are distinct but
    // nearly parallel (an ulp of std::atan2 would then decide the order).
    // Coincident points are allowed and are the interesting ties: MSVC's
    // std::sort is an insertion sort, hence stable, for at most 32 elements,
    // and Array.prototype.sort is stable, so the two agree on equivalent
    // elements for these sizes. Anything larger than 32 points would not be
    // comparable.
    int32_t mode = io.integer(0, 2);
    int32_t n = io.integer(2, 10);

    SortInput in{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        in = DrawSortInput(io, mode, n);
        std::vector<std::array<double, 2>> W = SortW(in);
        accepted = NoZeroW(W) && AnglesAreWellSeparated(W);
    }
    if (!accepted)
    {
        in = FallbackSortInput(io, n);
    }
    GiveSortInput(io, in);

    std::vector<size_t> byAngle{}, byGeometry{};
    SortPointsOnCircle<double>::ByAngle(in.P, in.C, in.D, in.sortCCW, byAngle);
    SortPointsOnCircle<double>::ByGeometry(in.P, in.C, in.D, in.sortCCW, byGeometry);

    io.outInt(static_cast<int32_t>(byAngle.size()));
    for (size_t i = 0; i < byAngle.size(); ++i)
    {
        io.outInt(static_cast<int32_t>(byAngle[i]));
    }
    for (size_t i = 0; i < byGeometry.size(); ++i)
    {
        io.outInt(static_cast<int32_t>(byGeometry[i]));
    }

    // Reference check: with a strict weak ordering and no zero W, the two
    // independent algorithms must produce the same permutation.
    io.outBool(byAngle == byGeometry);
}

ORACLE_CASE("SortPointsOnCircle.byGeometry.ties")
{
    // ByGeometry uses arithmetic only, so exact ties (coincident points and
    // points on a common ray from the centre) are compared here without the
    // angular separation requirement of the combined case. Points coinciding
    // with the centre remain excluded (finding #394).
    int32_t n = io.integer(2, 10);

    SortInput in{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        in = DrawSortInput(io, 2, n);
        accepted = NoZeroW(SortW(in));
    }
    if (!accepted)
    {
        in = FallbackSortInput(io, n);
    }
    GiveSortInput(io, in);

    std::vector<size_t> byGeometry{};
    SortPointsOnCircle<double>::ByGeometry(in.P, in.C, in.D, in.sortCCW, byGeometry);
    io.outInt(static_cast<int32_t>(byGeometry.size()));
    for (size_t i = 0; i < byGeometry.size(); ++i)
    {
        io.outInt(static_cast<int32_t>(byGeometry[i]));
    }
}

namespace
{
    // ---- the two circle constructions -----------------------------------

    void OutCircles(oracle::Ctx& io, size_t numCircles,
        std::array<Circle2<double>, 2> const& circle)
    {
        io.outInt(static_cast<int32_t>(numCircles));
        for (size_t i = 0; i < 2; ++i)
        {
            io.outVec(circle[i].center);
            io.outReal(circle[i].radius);
        }
    }

    // The 10 lattice vectors used to make |P-Q| exactly representable, so
    // that r = |P-Q|/2 gives an exactly zero 'argument' (the single-circle
    // branch).
    int32_t const pythag[10][2] =
    {
        { 3, 4 }, { 4, 3 }, { -3, 4 }, { 5, 0 }, { 0, 5 },
        { 6, 8 }, { 8, -6 }, { 1, 0 }, { 0, 2 }, { -2, 0 }
    };
}

ORACLE_CASE("CircleThroughTwoPointsSpecifiedRadius.compute")
{
    // Mode 0 constructs r = |P-Q|/2 exactly (argument == 0, one circle),
    // mode 1 makes r slightly smaller (argument < 0, no circle), mode 2 sets
    // P == Q (sqrLengthPmQ == 0, no circle), mode 3 is a lattice pair with a
    // lattice radius and mode 4 is uniform. Every mode records five doubles.
    int32_t mode = io.integer(0, 4);
    Vector2<double> P{ 0.0, 0.0 }, Q{ 0.0, 0.0 };
    double r = 0.0;

    if (mode == 0 || mode == 1)
    {
        int32_t k = io.rawInteger(0, 9);
        Vector2<double> base{ static_cast<double>(io.rawInteger(-4, 4)),
            static_cast<double>(io.rawInteger(-4, 4)) };
        Vector2<double> d{ static_cast<double>(pythag[k][0]),
            static_cast<double>(pythag[k][1]) };
        double length = std::sqrt(d[0] * d[0] + d[1] * d[1]);
        P = io.givenVec(base);
        Q = io.givenVec(Vector2<double>{ base[0] - d[0], base[1] - d[1] });
        r = io.given(mode == 0 ? 0.5 * length : 0.25 * length);
    }
    else if (mode == 2)
    {
        Vector2<double> base{ static_cast<double>(io.rawInteger(-4, 4)),
            static_cast<double>(io.rawInteger(-4, 4)) };
        P = io.givenVec(base);
        Q = io.givenVec(base);
        r = io.given(static_cast<double>(io.rawInteger(1, 5)));
    }
    else if (mode == 3)
    {
        P = io.latticeVec<2>(-5, 5);
        Q = io.latticeVec<2>(-5, 5);
        r = io.lattice(0, 6);
    }
    else
    {
        P = io.vec<2>(-5.0, 5.0);
        Q = io.vec<2>(-5.0, 5.0);
        r = io.real(0.0, 8.0);
    }

    std::array<Circle2<double>, 2> circle{};
    size_t numCircles = CircleThroughTwoPointsSpecifiedRadius(P, Q, r, circle);
    OutCircles(io, numCircles, circle);

    // Reference check: each reported circle must contain P and Q, i.e. the
    // squared distances from the centre must both equal r^2 to within the
    // conditioning of the construction.
    bool referenceOk = true;
    for (size_t i = 0; i < numCircles; ++i)
    {
        Vector2<double> dp = P - circle[i].center;
        Vector2<double> dq = Q - circle[i].center;
        double scale = std::max(1.0, r * r);
        referenceOk = referenceOk
            && std::fabs(Dot(dp, dp) - r * r) <= 1e-9 * scale
            && std::fabs(Dot(dq, dq) - r * r) <= 1e-9 * scale;
    }
    io.outBool(referenceOk);
}

ORACLE_CASE("CircleThroughPointSpecifiedTangentAndRadius.compute")
{
    // The normal must be unit length. Mode 0 uses an axis-aligned normal
    // (exactly unit) and a lattice point pair, so the signed distance s is an
    // exact integer; the radius is then chosen to hit s == 0, s == r,
    // s == 2*r, s > 2*r and 0 < s < r / r < s < 2*r exactly. Mode 1 is a
    // uniform normal and radius. Every mode records seven doubles.
    int32_t mode = io.integer(0, 1);
    int32_t branch = io.integer(0, 5);
    Vector2<double> P{ 0.0, 0.0 }, A{ 0.0, 0.0 }, N{ 0.0, 0.0 };
    double r = 0.0;

    if (mode == 0)
    {
        int32_t axis = io.rawInteger(0, 3);
        Vector2<double> n{ 0.0, 0.0 };
        n[axis % 2] = (axis < 2 ? 1.0 : -1.0);
        Vector2<double> a{ static_cast<double>(io.rawInteger(-4, 4)),
            static_cast<double>(io.rawInteger(-4, 4)) };
        // s = Dot(n, p - a); choose p so that s takes the wanted value. The
        // sign of s is drawn independently so that upstream's "negate N and
        // s" branch is reached as often as not.
        double magnitude = static_cast<double>(io.rawInteger(1, 4)) * 2.0;
        double sign = (io.rawInteger(0, 1) != 0 ? 1.0 : -1.0);
        double wanted = (branch == 0 ? 0.0 : sign * magnitude);
        Vector2<double> tangent{ n[1], -n[0] };
        double along = static_cast<double>(io.rawInteger(-4, 4));
        Vector2<double> p{ a[0] + wanted * n[0] + along * tangent[0],
            a[1] + wanted * n[1] + along * tangent[1] };
        P = io.givenVec(p);
        A = io.givenVec(a);
        N = io.givenVec(n);

        double radius = 1.0;
        if (branch == 0) { radius = static_cast<double>(io.rawInteger(1, 4)); }
        else if (branch == 1) { radius = magnitude; }          // s == r
        else if (branch == 2) { radius = 0.5 * magnitude; }    // s == 2*r
        else if (branch == 3) { radius = 0.25 * magnitude; }   // s > 2*r
        else if (branch == 4) { radius = 0.75 * magnitude; }   // r < s < 2*r
        else { radius = 2.0 * magnitude; }                     // 0 < s < r
        r = io.given(radius);
    }
    else
    {
        P = io.vec<2>(-5.0, 5.0);
        A = io.vec<2>(-5.0, 5.0);
        N = io.unit<2>();
        r = io.real(0.25, 6.0);
    }

    std::array<Circle2<double>, 2> circle{};
    size_t numCircles = CircleThroughPointSpecifiedTangentAndRadius(P, A, N, r, circle);
    OutCircles(io, numCircles, circle);

    // Reference check: each reported circle must contain P and be tangent to
    // the line Dot(N, X - A) = 0, i.e. |Dot(N, C - A)| == r.
    bool referenceOk = true;
    for (size_t i = 0; i < numCircles; ++i)
    {
        Vector2<double> dp = P - circle[i].center;
        double scale = std::max(1.0, r * r);
        referenceOk = referenceOk
            && std::fabs(Dot(dp, dp) - r * r) <= 1e-9 * scale
            && std::fabs(std::fabs(Dot(N, circle[i].center - A)) - r) <= 1e-9 * std::max(1.0, r);
    }
    io.outBool(referenceOk);
}

namespace
{
    // ---- ConvexHullSimplePolygon ----------------------------------------
    //
    // The input must be a simple counterclockwise polygon. Mode 0 walks the
    // boundary lattice points of an axis-aligned rectangle, which puts long
    // collinear runs on every edge so that WhichSide returns 0 and the
    // >= 0 / <= 0 branches are both exercised. Modes 1 and 2 are star-shaped
    // polygons about the origin (strictly increasing vertex angles with
    // arbitrary positive radii is enough for simplicity and for
    // counterclockwise order), lattice and uniform respectively. Only the
    // final vertex coordinates are recorded, so the sin/cos of mode 2 never
    // enters a compared computation.

    std::vector<Vector2<double>> MakeSimplePolygon(oracle::Ctx& io, int32_t mode)
    {
        std::vector<Vector2<double>> polygon{};
        if (mode == 0)
        {
            int32_t w = io.rawInteger(2, 4);
            int32_t h = io.rawInteger(2, 4);
            double ox = static_cast<double>(io.rawInteger(-3, 3));
            double oy = static_cast<double>(io.rawInteger(-3, 3));
            for (int32_t x = 0; x < w; ++x)
            {
                polygon.push_back(Vector2<double>{ ox + x, oy });
            }
            for (int32_t y = 0; y < h; ++y)
            {
                polygon.push_back(Vector2<double>{ ox + w, oy + y });
            }
            for (int32_t x = w; x > 0; --x)
            {
                polygon.push_back(Vector2<double>{ ox + x, oy + h });
            }
            for (int32_t y = h; y > 0; --y)
            {
                polygon.push_back(Vector2<double>{ ox, oy + y });
            }
        }
        else if (mode == 1)
        {
            int32_t n = io.rawInteger(3, 12);
            int32_t start = io.rawInteger(0, 11);
            std::vector<int32_t> pick{};
            for (int32_t k = 0; k < 12 && static_cast<int32_t>(pick.size()) < n; ++k)
            {
                if (io.rawInteger(0, 11) < n)
                {
                    pick.push_back((start + k) % 12);
                }
            }
            if (pick.size() < 3)
            {
                pick = { start % 12, (start + 4) % 12, (start + 8) % 12 };
            }
            for (int32_t k : pick)
            {
                double scale = static_cast<double>(io.rawInteger(1, 3));
                polygon.push_back(Vector2<double>{
                    scale * static_cast<double>(circle5[k][0]),
                    scale * static_cast<double>(circle5[k][1]) });
            }
        }
        else
        {
            int32_t n = io.rawInteger(3, 10);
            std::vector<double> angles(static_cast<size_t>(n));
            for (int32_t k = 0; k < n; ++k)
            {
                angles[static_cast<size_t>(k)] = io.raw(0.0, 6.28318530717958647692);
            }
            std::sort(angles.begin(), angles.end());
            for (int32_t k = 0; k < n; ++k)
            {
                double radius = io.raw(0.5, 5.0);
                polygon.push_back(Vector2<double>{
                    radius * std::cos(angles[static_cast<size_t>(k)]),
                    radius * std::sin(angles[static_cast<size_t>(k)]) });
            }
        }
        return polygon;
    }

    // The hull must be convex, counterclockwise and contain every polygon
    // vertex. This is an independent check of the Melkman construction.
    bool HullIsValid(std::vector<Vector2<double>> const& polygon,
        std::vector<size_t> const& hull)
    {
        size_t const m = hull.size();
        if (m < 3)
        {
            return false;
        }
        double area = 0.0;
        for (size_t i = 0; i < m; ++i)
        {
            Vector2<double> const& a = polygon[hull[i]];
            Vector2<double> const& b = polygon[hull[(i + 1) % m]];
            area += a[0] * b[1] - b[0] * a[1];
            Vector2<double> const& c = polygon[hull[(i + 2) % m]];
            double turn = (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]);
            if (turn < 0.0)
            {
                return false;
            }
        }
        if (!(area > 0.0))
        {
            return false;
        }
        for (auto const& p : polygon)
        {
            for (size_t i = 0; i < m; ++i)
            {
                Vector2<double> const& a = polygon[hull[i]];
                Vector2<double> const& b = polygon[hull[(i + 1) % m]];
                double side = (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
                if (side < 0.0)
                {
                    return false;
                }
            }
        }
        return true;
    }
}

ORACLE_CASE("ConvexHullSimplePolygon.compute")
{
    int32_t mode = io.integer(0, 2);
    std::vector<Vector2<double>> polygon = MakeSimplePolygon(io, mode);
    int32_t n = io.integer(static_cast<int32_t>(polygon.size()),
        static_cast<int32_t>(polygon.size()));
    for (int32_t i = 0; i < n; ++i)
    {
        io.givenVec(polygon[static_cast<size_t>(i)]);
    }

    ConvexHullSimplePolygon<double> query{};
    std::vector<size_t> hull{};
    query(polygon, hull);

    io.outInt(static_cast<int32_t>(hull.size()));
    for (size_t i = 0; i < hull.size(); ++i)
    {
        io.outInt(static_cast<int32_t>(hull[i]));
    }
    io.outBool(HullIsValid(polygon, hull));
}

namespace
{
    // ---- convex lattice polytopes ---------------------------------------
    //
    // Three convex polytopes with counterclockwise-when-seen-from-outside
    // triangle faces, used by ConvexPolyhedron3, ExtremalQuery3 and
    // ExtremalQuery3PRJ. All vertices are on the integer lattice, so the face
    // normals and the support values are exact integers and the extremal
    // queries have exact ties.

    int32_t const tetraIndices[12] =
    {
        0, 2, 1,  0, 1, 3,  0, 3, 2,  1, 2, 3
    };

    int32_t const cubeIndices[36] =
    {
        0, 3, 2,  0, 2, 1,  4, 5, 6,  4, 6, 7,
        0, 1, 5,  0, 5, 4,  1, 2, 6,  1, 6, 5,
        2, 3, 7,  2, 7, 6,  3, 0, 4,  3, 4, 7
    };

    int32_t const octaIndices[24] =
    {
        0, 2, 4,  2, 1, 4,  1, 3, 4,  3, 0, 4,
        2, 0, 5,  1, 2, 5,  3, 1, 5,  0, 3, 5
    };

    // Build one of the three polytopes, scaled and translated on the lattice.
    // Nothing is recorded here; the caller records the final vertices and
    // indices.
    void MakePolytope(oracle::Ctx& io, int32_t shape,
        std::vector<Vector3<double>>& vertices, std::vector<int32_t>& indices)
    {
        double sx = static_cast<double>(io.rawInteger(1, 4));
        double sy = static_cast<double>(io.rawInteger(1, 4));
        double sz = static_cast<double>(io.rawInteger(1, 4));
        double ox = static_cast<double>(io.rawInteger(-3, 3));
        double oy = static_cast<double>(io.rawInteger(-3, 3));
        double oz = static_cast<double>(io.rawInteger(-3, 3));

        vertices.clear();
        indices.clear();
        if (shape == 0)
        {
            vertices =
            {
                Vector3<double>{ ox, oy, oz },
                Vector3<double>{ ox + sx, oy, oz },
                Vector3<double>{ ox, oy + sy, oz },
                Vector3<double>{ ox, oy, oz + sz }
            };
            indices.assign(tetraIndices, tetraIndices + 12);
        }
        else if (shape == 1)
        {
            for (int32_t k = 0; k < 8; ++k)
            {
                double x = ((k == 1 || k == 2 || k == 5 || k == 6) ? sx : 0.0);
                double y = ((k == 2 || k == 3 || k == 6 || k == 7) ? sy : 0.0);
                double z = (k >= 4 ? sz : 0.0);
                vertices.push_back(Vector3<double>{ ox + x, oy + y, oz + z });
            }
            indices.assign(cubeIndices, cubeIndices + 36);
        }
        else
        {
            vertices =
            {
                Vector3<double>{ ox + sx, oy, oz },
                Vector3<double>{ ox - sx, oy, oz },
                Vector3<double>{ ox, oy + sy, oz },
                Vector3<double>{ ox, oy - sy, oz },
                Vector3<double>{ ox, oy, oz + sz },
                Vector3<double>{ ox, oy, oz - sz }
            };
            indices.assign(octaIndices, octaIndices + 24);
        }
    }

    // Record a mesh: the vertex count, the vertices, the index count, the
    // indices.
    void GiveMesh(oracle::Ctx& io, std::vector<Vector3<double>> const& vertices,
        std::vector<int32_t> const& indices)
    {
        int32_t numVertices = static_cast<int32_t>(vertices.size());
        io.integer(numVertices, numVertices);
        for (auto const& v : vertices)
        {
            io.givenVec(v);
        }
        int32_t numIndices = static_cast<int32_t>(indices.size());
        io.integer(numIndices, numIndices);
        for (int32_t i : indices)
        {
            io.integer(i, i);
        }
    }
}

ORACLE_CASE("ConvexPolyhedron3.construct")
{
    // Mode 0-2 are the three convex lattice polytopes; mode 3 is an
    // unvalidated mesh of random lattice vertices and random triangles (the
    // class does not check convexity); mode 4 appends one or two stray
    // indices, which upstream silently drops because it validates
    // indices.size() >= 12 but not indices.size() % 3 == 0 (finding #175,
    // preserved); mode 5 is a too-small input for which the constructor
    // leaves every member array empty.
    int32_t mode = io.integer(0, 5);
    std::vector<Vector3<double>> vertices{};
    std::vector<int32_t> indices{};

    if (mode <= 2)
    {
        MakePolytope(io, mode, vertices, indices);
    }
    else if (mode == 3 || mode == 4)
    {
        int32_t numVertices = io.rawInteger(4, 7);
        for (int32_t i = 0; i < numVertices; ++i)
        {
            vertices.push_back(Vector3<double>{
                static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) });
        }
        int32_t numTriangles = io.rawInteger(4, 6);
        for (int32_t t = 0; t < 3 * numTriangles; ++t)
        {
            indices.push_back(io.rawInteger(0, numVertices - 1));
        }
        if (mode == 4)
        {
            int32_t extra = io.rawInteger(1, 2);
            for (int32_t k = 0; k < extra; ++k)
            {
                indices.push_back(io.rawInteger(0, numVertices - 1));
            }
        }
    }
    else
    {
        int32_t numVertices = io.rawInteger(1, 4);
        for (int32_t i = 0; i < numVertices; ++i)
        {
            vertices.push_back(Vector3<double>{
                static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)),
                static_cast<double>(io.rawInteger(-4, 4)) });
        }
        int32_t numIndices = io.rawInteger(0, 11);
        for (int32_t k = 0; k < numIndices; ++k)
        {
            indices.push_back(io.rawInteger(0, numVertices - 1));
        }
    }

    GiveMesh(io, vertices, indices);
    bool wantPlanes = io.boolean();
    bool wantAlignedBox = io.boolean();

    std::vector<Vector3<double>> movedVertices = vertices;
    std::vector<int32_t> movedIndices = indices;
    ConvexPolyhedron3<double> polyhedron(std::move(movedVertices),
        std::move(movedIndices), wantPlanes, wantAlignedBox);

    io.outInt(static_cast<int32_t>(polyhedron.vertices.size()));
    io.outInt(static_cast<int32_t>(polyhedron.indices.size()));
    io.outInt(static_cast<int32_t>(polyhedron.planes.size()));
    for (auto const& plane : polyhedron.planes)
    {
        io.outVec(plane);
    }
    io.outVec(polyhedron.alignedBox.min);
    io.outVec(polyhedron.alignedBox.max);

    // Reference check for the convex modes: every generated plane must have
    // an outward normal, i.e. Dot(N, V) + d <= 0 for every polytope vertex V.
    bool outward = true;
    if (mode <= 2 && wantPlanes)
    {
        for (auto const& plane : polyhedron.planes)
        {
            for (auto const& v : polyhedron.vertices)
            {
                double value = plane[0] * v[0] + plane[1] * v[1] + plane[2] * v[2] + plane[3];
                outward = outward && (value <= 0.0);
            }
        }
    }
    io.outBool(outward);

    // Regenerating after construction must reproduce the same values.
    polyhedron.GeneratePlanes();
    polyhedron.GenerateAlignedBox();
    io.outInt(static_cast<int32_t>(polyhedron.planes.size()));
    for (auto const& plane : polyhedron.planes)
    {
        io.outVec(plane);
    }
    io.outVec(polyhedron.alignedBox.min);
    io.outVec(polyhedron.alignedBox.max);
}

ORACLE_CASE("ExtremalQuery3PRJ.getExtremeVertices")
{
    // The base class ExtremalQuery3 computes the face normals with UnitCross
    // (a cross product and a sqrt) and ExtremalQuery3PRJ projects the
    // vertices onto the direction after subtracting the vertex average, so
    // the whole query is arithmetic. Directions drawn on the lattice are
    // frequently perpendicular to a face or parallel to an edge of these
    // polytopes, which is where the strict comparisons produce ties; the
    // first vertex visited in ascending std::set order wins.
    int32_t shape = io.integer(0, 2);
    std::vector<Vector3<double>> vertices{};
    std::vector<int32_t> indices{};
    MakePolytope(io, shape, vertices, indices);
    GiveMesh(io, vertices, indices);

    auto vertexPool = std::make_shared<std::vector<Vector3<double>>>(vertices);
    Polyhedron3<double> polytope(vertexPool, static_cast<int32_t>(indices.size()),
        indices.data(), true);
    ExtremalQuery3PRJ<double> query(polytope);

    auto const& normals = query.GetFaceNormals();
    io.outInt(static_cast<int32_t>(normals.size()));
    for (auto const& n : normals)
    {
        io.outVec(n);
    }
    io.outInt(static_cast<int32_t>(query.GetPolytope().GetUniqueIndices().size()));

    for (int32_t q = 0; q < 4; ++q)
    {
        int32_t dirMode = io.integer(0, 1);
        Vector3<double> direction = (dirMode == 0
            ? io.latticeDir<3>(-2, 2) : io.unit<3>());
        int32_t positiveDirection = 0, negativeDirection = 0;
        query.GetExtremeVertices(direction, positiveDirection, negativeDirection);
        io.outInt(positiveDirection);
        io.outInt(negativeDirection);

        // Reference check: brute force over the unrotated, unshifted
        // vertices. The reported extreme vertices must realize the maximum
        // and the minimum of Dot(direction, V) over the polytope vertices.
        double best = -std::numeric_limits<double>::max();
        double worst = std::numeric_limits<double>::max();
        for (int32_t i : polytope.GetUniqueIndices())
        {
            double d = Dot(direction, vertices[static_cast<size_t>(i)]);
            best = std::max(best, d);
            worst = std::min(worst, d);
        }
        double dPos = Dot(direction, vertices[static_cast<size_t>(positiveDirection)]);
        double dNeg = Dot(direction, vertices[static_cast<size_t>(negativeDirection)]);
        io.outBool(dPos == best && dNeg == worst);
    }
}

namespace
{
    // ---- NearestNeighborQuery -------------------------------------------
    //
    // std::nth_element only guarantees the partition postconditions, not a
    // particular permutation, and the port replaces it with its own
    // quickselect. The permutation is therefore not comparable and is never
    // emitted; what is emitted is the tree shape, the split values and, for
    // every leaf, the SORTED list of the original site indices it holds.
    // Those are all determined by the postconditions as long as the
    // coordinates along each split axis are pairwise distinct, which the
    // generator guarantees by drawing a permutation of distinct integers for
    // each axis.

    using Site3 = PositionSite<3, double>;
    using NNQuery3 = NearestNeighborQuery<3, double, Site3>;

    // n sites whose coordinates along each axis are pairwise distinct.
    std::vector<Vector3<double>> MakeDistinctSites(oracle::Ctx& io, int32_t n)
    {
        std::vector<Vector3<double>> points(static_cast<size_t>(n));
        for (int32_t d = 0; d < 3; ++d)
        {
            std::vector<int32_t> values(static_cast<size_t>(n));
            for (int32_t i = 0; i < n; ++i)
            {
                values[static_cast<size_t>(i)] = i;
            }
            // Fisher-Yates with the oracle's own generator (std::shuffle with
            // a standard engine would not be reproducible across libraries,
            // but nothing here is recorded: only the final coordinates are).
            for (int32_t i = n - 1; i > 0; --i)
            {
                int32_t j = io.rawInteger(0, i);
                std::swap(values[static_cast<size_t>(i)], values[static_cast<size_t>(j)]);
            }
            double offset = static_cast<double>(io.rawInteger(-4, 4));
            for (int32_t i = 0; i < n; ++i)
            {
                points[static_cast<size_t>(i)][d] =
                    offset + static_cast<double>(values[static_cast<size_t>(i)]);
            }
        }
        return points;
    }

    // n sites drawn on a coarse lattice, so that coordinates repeat along
    // every axis.
    std::vector<Vector3<double>> MakeTiedSites(oracle::Ctx& io, int32_t n)
    {
        std::vector<Vector3<double>> points(static_cast<size_t>(n));
        for (int32_t i = 0; i < n; ++i)
        {
            for (int32_t d = 0; d < 3; ++d)
            {
                points[static_cast<size_t>(i)][d] =
                    static_cast<double>(io.rawInteger(-2, 2));
            }
        }
        return points;
    }

    void OutTree(oracle::Ctx& io, NNQuery3 const& query, int32_t n)
    {
        io.outInt(query.GetMaxLeafSize());
        io.outInt(query.GetMaxLevel());
        io.outInt(query.GetDepth());
        io.outInt(query.GetLargestNodeSize());
        io.outInt(query.GetNumNodes());

        auto const& sortedPoints = query.GetSortedPoints();
        io.outInt(static_cast<int32_t>(sortedPoints.size()));

        for (auto const& node : query.GetNodes())
        {
            io.outReal(node.split);
            io.outInt(node.axis);
            io.outInt(node.numSites);
            io.outInt(node.siteOffset);
            io.outInt(node.left);
            io.outInt(node.right);
            if (node.siteOffset != -1)
            {
                // Canonicalize: the site indices of a leaf, sorted. The order
                // inside the leaf is whatever std::nth_element left there and
                // is not comparable.
                std::vector<int32_t> leaf{};
                for (int32_t k = 0; k < node.numSites; ++k)
                {
                    leaf.push_back(sortedPoints[static_cast<size_t>(node.siteOffset + k)].second);
                }
                std::sort(leaf.begin(), leaf.end());
                for (int32_t v : leaf)
                {
                    io.outInt(v);
                }
            }
        }
        io.outInt(n);
    }
}

ORACLE_CASE("NearestNeighborQuery.build")
{
    int32_t n = io.integer(1, 12);
    int32_t maxLeafSize = io.integer(1, 4);
    int32_t maxLevel = io.integer(1, 5);
    std::vector<Vector3<double>> points = MakeDistinctSites(io, n);
    std::vector<Site3> sites{};
    for (int32_t i = 0; i < n; ++i)
    {
        sites.push_back(Site3(io.givenVec(points[static_cast<size_t>(i)])));
    }

    NNQuery3 query(sites, maxLeafSize, maxLevel);
    OutTree(io, query, n);
}

ORACLE_CASE("NearestNeighborQuery.build.tiedCoordinates")
{
    // Coordinates repeat along every axis, so the median value is attained by
    // several sites and which of them std::nth_element leaves on each side of
    // the median is unspecified. Measured: emitting the leaf contents makes
    // 11 of 20 records disagree, and so does the split value of any internal
    // node below the root, because that node's subrange holds a different
    // multiset. Only what follows from the site counts alone is emitted here:
    // the depth, the largest leaf size, the node count, the node ranges and
    // the root split (an order statistic over the whole array, hence
    // determined by value). See the report's "Not covered" section.
    int32_t n = io.integer(4, 12);
    int32_t maxLeafSize = io.integer(1, 3);
    int32_t maxLevel = io.integer(1, 5);
    std::vector<Vector3<double>> points = MakeTiedSites(io, n);
    std::vector<Site3> sites{};
    for (int32_t i = 0; i < n; ++i)
    {
        sites.push_back(Site3(io.givenVec(points[static_cast<size_t>(i)])));
    }

    NNQuery3 query(sites, maxLeafSize, maxLevel);
    io.outInt(query.GetDepth());
    io.outInt(query.GetLargestNodeSize());
    io.outInt(query.GetNumNodes());
    io.outReal(query.GetNodes()[0].split);

    for (auto const& node : query.GetNodes())
    {
        io.outInt(node.axis);
        io.outInt(node.numSites);
        io.outInt(node.siteOffset);
        io.outInt(node.left);
        io.outInt(node.right);
    }
}

ORACLE_CASE("NearestNeighborQuery.build.maxLevelAssert")
{
    // Throw-parity case: the constructor asserts 0 < maxLevel <= 32. The
    // sampled values straddle both ends of the valid range.
    int32_t const levels[6] = { -1, 0, 1, 5, 32, 33 };
    int32_t n = io.integer(1, 6);
    int32_t maxLeafSize = io.integer(1, 3);
    int32_t chosen = levels[io.rawInteger(0, 5)];
    int32_t maxLevel = io.integer(chosen, chosen);
    std::vector<Vector3<double>> points = MakeDistinctSites(io, n);
    std::vector<Site3> sites{};
    for (int32_t i = 0; i < n; ++i)
    {
        sites.push_back(Site3(io.givenVec(points[static_cast<size_t>(i)])));
    }

    NNQuery3 query(sites, maxLeafSize, maxLevel);
    OutTree(io, query, n);
}

ORACLE_CASE("NearestNeighborQuery.findNeighbors.all")
{
    // maxNeighbors is at least the number of sites, so every site inside the
    // radius is retained and the result does not depend on the order in
    // which the leaves were visited. Distance ties are therefore allowed and
    // are the point of this case: the heap pops them in decreasing site
    // index order.
    int32_t n = io.integer(1, 12);
    int32_t maxLeafSize = io.integer(1, 4);
    int32_t maxLevel = io.integer(1, 5);
    std::vector<Vector3<double>> points = MakeDistinctSites(io, n);
    std::vector<Site3> sites{};
    for (int32_t i = 0; i < n; ++i)
    {
        sites.push_back(Site3(io.givenVec(points[static_cast<size_t>(i)])));
    }

    NNQuery3 query(sites, maxLeafSize, maxLevel);
    OutTree(io, query, n);

    for (int32_t q = 0; q < 4; ++q)
    {
        Vector3<double> point = io.latticeVec<3>(-6, 6);
        double radius = io.lattice(0, 8);
        std::array<int32_t, 12> neighbors{};
        int32_t numNeighbors = query.FindNeighbors<12>(point, radius, neighbors);
        io.outInt(numNeighbors);
        for (int32_t k = 0; k < numNeighbors; ++k)
        {
            io.outInt(neighbors[static_cast<size_t>(k)]);
        }

        // Reference check: brute force over all sites.
        int32_t expected = 0;
        for (int32_t i = 0; i < n; ++i)
        {
            Vector3<double> diff = points[static_cast<size_t>(i)] - point;
            if (Dot(diff, diff) <= radius * radius)
            {
                ++expected;
            }
        }
        io.outBool(numNeighbors == expected);
    }
}

ORACLE_CASE("NearestNeighborQuery.findNeighbors.limited")
{
    // maxNeighbors is smaller than the number of sites. The retained set then
    // depends on the visiting order whenever two candidates are equidistant
    // from the query point (upstream replaces only on a strict improvement),
    // and the visiting order depends on the std::nth_element permutation,
    // which is not comparable. The generator therefore rejects query points
    // with equidistant sites; it is capped and falls back to a point with the
    // fewest ties seen.
    int32_t n = io.integer(6, 12);
    int32_t maxLeafSize = io.integer(1, 3);
    int32_t maxLevel = io.integer(1, 5);
    int32_t maxNeighbors = io.integer(1, 4);
    std::vector<Vector3<double>> points = MakeDistinctSites(io, n);
    std::vector<Site3> sites{};
    for (int32_t i = 0; i < n; ++i)
    {
        sites.push_back(Site3(io.givenVec(points[static_cast<size_t>(i)])));
    }

    NNQuery3 query(sites, maxLeafSize, maxLevel);

    for (int32_t q = 0; q < 3; ++q)
    {
        Vector3<double> point{ 0.0, 0.0, 0.0 };
        Vector3<double> best{ 0.0, 0.0, 0.0 };
        int32_t bestTies = -1;
        double radius = 0.0, bestRadius = 0.0;
        for (int32_t attempt = 0; attempt < 32; ++attempt)
        {
            for (int32_t d = 0; d < 3; ++d)
            {
                point[d] = 0.5 * static_cast<double>(io.rawInteger(-12, 12));
            }
            radius = static_cast<double>(io.rawInteger(1, 8));
            std::vector<double> sqr{};
            for (int32_t i = 0; i < n; ++i)
            {
                Vector3<double> diff = points[static_cast<size_t>(i)] - point;
                double d2 = Dot(diff, diff);
                if (d2 <= radius * radius)
                {
                    sqr.push_back(d2);
                }
            }
            std::sort(sqr.begin(), sqr.end());
            int32_t ties = 0;
            for (size_t k = 0; k + 1 < sqr.size(); ++k)
            {
                if (sqr[k] == sqr[k + 1])
                {
                    ++ties;
                }
            }
            if (bestTies < 0 || ties < bestTies)
            {
                bestTies = ties;
                best = point;
                bestRadius = radius;
            }
            if (ties == 0)
            {
                break;
            }
        }
        Vector3<double> usedPoint = io.givenVec(best);
        double usedRadius = io.given(bestRadius);

        std::array<int32_t, 4> neighbors{};
        int32_t numNeighbors = 0;
        if (maxNeighbors == 1)
        {
            std::array<int32_t, 1> nb{};
            numNeighbors = query.FindNeighbors<1>(usedPoint, usedRadius, nb);
            for (int32_t k = 0; k < numNeighbors; ++k) { neighbors[static_cast<size_t>(k)] = nb[static_cast<size_t>(k)]; }
        }
        else if (maxNeighbors == 2)
        {
            std::array<int32_t, 2> nb{};
            numNeighbors = query.FindNeighbors<2>(usedPoint, usedRadius, nb);
            for (int32_t k = 0; k < numNeighbors; ++k) { neighbors[static_cast<size_t>(k)] = nb[static_cast<size_t>(k)]; }
        }
        else if (maxNeighbors == 3)
        {
            std::array<int32_t, 3> nb{};
            numNeighbors = query.FindNeighbors<3>(usedPoint, usedRadius, nb);
            for (int32_t k = 0; k < numNeighbors; ++k) { neighbors[static_cast<size_t>(k)] = nb[static_cast<size_t>(k)]; }
        }
        else
        {
            std::array<int32_t, 4> nb{};
            numNeighbors = query.FindNeighbors<4>(usedPoint, usedRadius, nb);
            for (int32_t k = 0; k < numNeighbors; ++k) { neighbors[static_cast<size_t>(k)] = nb[static_cast<size_t>(k)]; }
        }

        io.outInt(numNeighbors);
        for (int32_t k = 0; k < numNeighbors; ++k)
        {
            io.outInt(neighbors[static_cast<size_t>(k)]);
        }

        // Reference check: the returned indices must be the maxNeighbors
        // closest sites within the radius, computed by brute force.
        std::vector<std::pair<double, int32_t>> all{};
        for (int32_t i = 0; i < n; ++i)
        {
            Vector3<double> diff = points[static_cast<size_t>(i)] - usedPoint;
            double d2 = Dot(diff, diff);
            if (d2 <= usedRadius * usedRadius)
            {
                all.push_back({ d2, i });
            }
        }
        std::sort(all.begin(), all.end());
        int32_t expected = std::min(static_cast<int32_t>(all.size()), maxNeighbors);
        bool referenceOk = (numNeighbors == expected);
        if (referenceOk)
        {
            std::vector<int32_t> got(neighbors.begin(), neighbors.begin() + numNeighbors);
            std::sort(got.begin(), got.end());
            std::vector<int32_t> want{};
            for (int32_t k = 0; k < expected; ++k)
            {
                want.push_back(all[static_cast<size_t>(k)].second);
            }
            std::sort(want.begin(), want.end());
            referenceOk = (got == want);
        }
        io.outBool(referenceOk);
    }
}

namespace
{
    // ---- BoxManager -----------------------------------------------------
    //
    // The overlap container is a std::set<EdgeKey<false>>, so its iteration
    // order is already the lexicographic order of (V[0], V[1]); the port
    // sorts its Map values the same way and the case comment in the replay
    // says so. Initialize sorts the endpoints with std::sort on a comparator
    // that ties on (value, type); MSVC's std::sort is an insertion sort,
    // hence stable, for at most 32 elements, and Array.prototype.sort is
    // stable, so for at most 16 boxes the two builds produce the same
    // endpoint order. Larger inputs would not be comparable.

    AlignedBox3<double> DrawBox(oracle::Ctx& io)
    {
        AlignedBox3<double> box{};
        for (int32_t d = 0; d < 3; ++d)
        {
            double a = io.lattice(-4, 4);
            double b = io.lattice(-4, 4);
            box.min[d] = std::min(a, b);
            box.max[d] = std::max(a, b);
        }
        return box;
    }

    bool BoxesOverlap(AlignedBox3<double> const& b0, AlignedBox3<double> const& b1)
    {
        for (int32_t d = 0; d < 3; ++d)
        {
            if (b0.max[d] < b1.min[d] || b0.min[d] > b1.max[d])
            {
                return false;
            }
        }
        return true;
    }

    void OutOverlap(oracle::Ctx& io, BoxManager<double> const& manager,
        std::vector<AlignedBox3<double>> const& boxes)
    {
        auto const& overlap = manager.GetOverlap();
        io.outInt(static_cast<int32_t>(overlap.size()));
        for (auto const& key : overlap)
        {
            io.outInt(key.V[0]);
            io.outInt(key.V[1]);
        }

        // Reference check: brute force over all pairs.
        std::set<std::pair<int32_t, int32_t>> expected{};
        int32_t const n = static_cast<int32_t>(boxes.size());
        for (int32_t i = 0; i < n; ++i)
        {
            for (int32_t j = i + 1; j < n; ++j)
            {
                if (BoxesOverlap(boxes[static_cast<size_t>(i)], boxes[static_cast<size_t>(j)]))
                {
                    expected.insert({ i, j });
                }
            }
        }
        std::set<std::pair<int32_t, int32_t>> reported{};
        for (auto const& key : overlap)
        {
            reported.insert({ key.V[0], key.V[1] });
        }
        io.outBool(reported == expected);
    }
}

ORACLE_CASE("BoxManager.initializeAndUpdate")
{
    int32_t n = io.integer(2, 8);
    std::vector<AlignedBox3<double>> boxes(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i)
    {
        boxes[static_cast<size_t>(i)] = DrawBox(io);
    }

    BoxManager<double> manager(boxes);
    OutOverlap(io, manager, boxes);

    int32_t numRounds = io.integer(1, 3);
    for (int32_t round = 0; round < numRounds; ++round)
    {
        int32_t numMoves = io.integer(1, 3);
        for (int32_t m = 0; m < numMoves; ++m)
        {
            int32_t index = io.integer(0, n - 1);
            AlignedBox3<double> box = DrawBox(io);
            manager.SetBox(index, box);
        }
        manager.Update();
        OutOverlap(io, manager, boxes);

        int32_t probe = io.integer(0, n - 1);
        AlignedBox3<double> got{};
        manager.GetBox(probe, got);
        io.outVec(got.min);
        io.outVec(got.max);
    }

    // A full re-initialization from the moved boxes must agree with the
    // incrementally maintained set.
    manager.Initialize();
    OutOverlap(io, manager, boxes);
}

namespace
{
    // ---- InscribedFixedAspectRectInQuad ---------------------------------
    //
    // Execute uses std::atan2 only to compute the quadrant index
    // j = floor((2/pi) * angle) of each inner edge normal; the angle itself
    // never reaches the result. The generator therefore requires every
    // normal to have both components nonzero, i.e. no edge of the quad is
    // axis parallel. With lattice coordinates bounded by 6 the smallest
    // possible angle to a quadrant boundary is atan(1/12) = 0.083 radians,
    // twelve orders of magnitude above an ulp of atan2, so j is decided by
    // the geometry. Everything after that is arithmetic.
    //
    // Upstream asserts that the line of constraints 0 and 2 is not parallel
    // to constraints 1 and 3, which is false for legitimate convex quads
    // (finding #395 (a)), and its untoleranced interval test can report a
    // degenerate feasible interval empty (finding #395 (b)). The port
    // preserves both, so such quads are kept in the generator and compared
    // for throw parity.

    bool QuadIsConvexCCW(std::array<Vector2<double>, 4> const& quad)
    {
        for (int32_t i = 0; i < 4; ++i)
        {
            Vector2<double> e0 = quad[(i + 1) % 4] - quad[i];
            Vector2<double> e1 = quad[(i + 2) % 4] - quad[(i + 1) % 4];
            if (e0[0] * e1[1] - e0[1] * e1[0] <= 0.0)
            {
                return false;
            }
        }
        return true;
    }

    bool QuadHasNoAxisParallelEdge(std::array<Vector2<double>, 4> const& quad)
    {
        for (int32_t i = 0; i < 4; ++i)
        {
            Vector2<double> e = quad[(i + 1) % 4] - quad[i];
            if (e[0] == 0.0 || e[1] == 0.0)
            {
                return false;
            }
        }
        return true;
    }
}

ORACLE_CASE("InscribedFixedAspectRectInQuad.execute")
{
    int32_t mode = io.integer(0, 1);
    std::array<Vector2<double>, 4> quad{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 128 && !accepted; ++attempt)
    {
        for (int32_t i = 0; i < 4; ++i)
        {
            if (mode == 0)
            {
                quad[static_cast<size_t>(i)] =
                    Vector2<double>{ static_cast<double>(io.rawInteger(-6, 6)),
                        static_cast<double>(io.rawInteger(-6, 6)) };
            }
            else
            {
                quad[static_cast<size_t>(i)] =
                    Vector2<double>{ io.raw(-6.0, 6.0), io.raw(-6.0, 6.0) };
            }
        }
        accepted = QuadIsConvexCCW(quad) && QuadHasNoAxisParallelEdge(quad);
    }
    if (!accepted)
    {
        // A convex counterclockwise quad with no axis-parallel edge.
        quad = { Vector2<double>{ 3.0, 1.0 }, Vector2<double>{ 1.0, 4.0 },
            Vector2<double>{ -3.0, 1.0 }, Vector2<double>{ -1.0, -4.0 } };
    }
    for (int32_t i = 0; i < 4; ++i)
    {
        io.givenVec(quad[static_cast<size_t>(i)]);
    }

    double aspectRatio = 1.0;
    if (mode == 0)
    {
        double numerator = static_cast<double>(io.rawInteger(1, 4));
        double denominator = static_cast<double>(io.rawInteger(1, 4));
        aspectRatio = io.given(numerator / denominator);
    }
    else
    {
        aspectRatio = io.real(0.25, 4.0);
    }

    Vector2<double> rectOrigin{ 0.0, 0.0 };
    double rectWidth = 0.0, rectHeight = 0.0;
    bool isUnique = InscribedFixedAspectRectInQuad<double>::Execute(quad, aspectRatio,
        rectOrigin, rectWidth, rectHeight);
    io.outBool(isUnique);
    io.outVec(rectOrigin);
    io.outReal(rectWidth);
    io.outReal(rectHeight);

    // Reference check: the rectangle must lie inside the quad (each of its
    // four corners on the inner side of each of the four edges) and have the
    // requested aspect ratio.
    bool referenceOk = (rectWidth >= 0.0);
    std::array<Vector2<double>, 4> rect =
    {
        rectOrigin,
        Vector2<double>{ rectOrigin[0] + rectWidth, rectOrigin[1] },
        Vector2<double>{ rectOrigin[0] + rectWidth, rectOrigin[1] + rectHeight },
        Vector2<double>{ rectOrigin[0], rectOrigin[1] + rectHeight }
    };
    double scale = 0.0;
    for (int32_t i = 0; i < 4; ++i)
    {
        scale = std::max(scale, std::fabs(quad[static_cast<size_t>(i)][0]));
        scale = std::max(scale, std::fabs(quad[static_cast<size_t>(i)][1]));
    }
    for (int32_t i = 0; i < 4; ++i)
    {
        Vector2<double> e = quad[static_cast<size_t>((i + 1) % 4)] - quad[static_cast<size_t>(i)];
        for (int32_t k = 0; k < 4; ++k)
        {
            Vector2<double> d = rect[static_cast<size_t>(k)] - quad[static_cast<size_t>(i)];
            referenceOk = referenceOk && (e[0] * d[1] - e[1] * d[0] >= -1e-9 * scale * scale);
        }
    }
    io.outBool(referenceOk);
}

namespace
{
    // ---- TriangulateEC ---------------------------------------------------
    //
    // Upstream is TriangulateEC<InputType, ComputeType>; the port is
    // number-only, so the C++ side is instantiated as
    // TriangulateEC<double, double>. The triangle list is produced by a
    // deterministic ear-clipping order with no container whose order C++
    // leaves unspecified, so it is compared in the order the algorithm
    // produced it, triple by triple.
    //
    // The outer polygon of the hole cases is always the 12 lattice
    // directions of the circle of radius 5 scaled by 2 or 3, so every outer
    // vertex is at distance 10 or 15 from the origin and every outer edge is
    // at distance at least 10*cos(15 degrees) = 9.66 from it. Holes are built
    // inside the disk of radius 6 about a centre of length at most 4, hence
    // strictly inside the outer polygon, and the two holes of the
    // multiple-hole case, centred at (-4,0) and (4,0) with vertex offsets of
    // at most 2, are disjoint.

    int32_t const smallDirs[8][2] =
    {
        { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 },
        { -1, 0 }, { -1, -1 }, { 0, -1 }, { 1, -1 }
    };

    std::vector<Vector2<double>> MakeOuterPolygon(oracle::Ctx& io)
    {
        std::vector<Vector2<double>> P{};
        for (int32_t k = 0; k < 12; ++k)
        {
            double scale = static_cast<double>(io.rawInteger(2, 3));
            P.push_back(Vector2<double>{
                scale * static_cast<double>(circle5[k][0]),
                scale * static_cast<double>(circle5[k][1]) });
        }
        return P;
    }

    // A counterclockwise star polygon of 3 to 8 vertices about 'centre'.
    std::vector<Vector2<double>> MakeSmallStar(oracle::Ctx& io, double cx, double cy)
    {
        int32_t numVertices = io.rawInteger(3, 8);
        int32_t start = io.rawInteger(0, 7);
        std::vector<Vector2<double>> ccw{};
        for (int32_t k = 0; k < numVertices; ++k)
        {
            int32_t d = (start + (k * 8) / numVertices) % 8;
            double scale = static_cast<double>(io.rawInteger(1, 2));
            ccw.push_back(Vector2<double>{
                cx + scale * static_cast<double>(smallDirs[d][0]),
                cy + scale * static_cast<double>(smallDirs[d][1]) });
        }
        return ccw;
    }

    // The clockwise version, for a hole.
    std::vector<Vector2<double>> MakeHolePolygon(oracle::Ctx& io, double cx, double cy)
    {
        std::vector<Vector2<double>> ccw = MakeSmallStar(io, cx, cy);
        std::reverse(ccw.begin(), ccw.end());
        return ccw;
    }

    double ShoelaceArea(std::vector<Vector2<double>> const& points,
        std::vector<int32_t> const& polygon)
    {
        double twiceArea = 0.0;
        size_t const m = polygon.size();
        for (size_t i = 0; i < m; ++i)
        {
            Vector2<double> const& a = points[static_cast<size_t>(polygon[i])];
            Vector2<double> const& b = points[static_cast<size_t>(polygon[(i + 1) % m])];
            twiceArea += a[0] * b[1] - b[0] * a[1];
        }
        return twiceArea;
    }

    // The triangles must be counterclockwise (or degenerate) and their total
    // signed area must equal the signed area of the region. With lattice
    // coordinates every shoelace term is an exact integer, so the comparison
    // is exact.
    bool TriangulationIsValid(std::vector<Vector2<double>> const& points,
        std::vector<std::array<int32_t, 3>> const& triangles, double expectedTwiceArea)
    {
        double twiceArea = 0.0;
        for (auto const& t : triangles)
        {
            Vector2<double> const& a = points[static_cast<size_t>(t[0])];
            Vector2<double> const& b = points[static_cast<size_t>(t[1])];
            Vector2<double> const& c = points[static_cast<size_t>(t[2])];
            double d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            if (d < 0.0)
            {
                return false;
            }
            twiceArea += d;
        }
        return twiceArea == expectedTwiceArea;
    }

    void GivePoints(oracle::Ctx& io, std::vector<Vector2<double>> const& points)
    {
        int32_t n = static_cast<int32_t>(points.size());
        io.integer(n, n);
        for (auto const& p : points)
        {
            io.givenVec(p);
        }
    }

    void GivePolygon(oracle::Ctx& io, std::vector<int32_t> const& polygon)
    {
        int32_t n = static_cast<int32_t>(polygon.size());
        io.integer(n, n);
        for (int32_t i : polygon)
        {
            io.integer(i, i);
        }
    }

    void OutTriangles(oracle::Ctx& io,
        std::vector<std::array<int32_t, 3>> const& triangles)
    {
        io.outInt(static_cast<int32_t>(triangles.size()));
        for (auto const& t : triangles)
        {
            io.outInt(t[0]);
            io.outInt(t[1]);
            io.outInt(t[2]);
        }
    }
}

ORACLE_CASE("TriangulateEC.triangulate")
{
    // The whole point array is one simple counterclockwise polygon. Mode 0
    // puts long collinear runs on the boundary of an axis-aligned rectangle,
    // mode 1 is a star-shaped lattice polygon and mode 2 a star-shaped
    // uniform one.
    int32_t mode = io.integer(0, 2);
    std::vector<Vector2<double>> points = MakeSimplePolygon(io, mode);
    GivePoints(io, points);

    TriangulateEC<double, double> triangulator(points);
    triangulator();
    OutTriangles(io, triangulator.GetTriangles());

    std::vector<int32_t> polygon(points.size());
    for (size_t i = 0; i < points.size(); ++i)
    {
        polygon[i] = static_cast<int32_t>(i);
    }
    bool referenceOk = (triangulator.GetTriangles().size() == points.size() - 2)
        && TriangulationIsValid(points, triangulator.GetTriangles(),
            ShoelaceArea(points, polygon));
    io.outBool(referenceOk);
}

ORACLE_CASE("TriangulateEC.triangulatePolygon")
{
    // The point pool holds unused points as well; the polygon is given by an
    // index list into the pool, in counterclockwise order.
    int32_t mode = io.integer(0, 2);
    std::vector<Vector2<double>> polygonPoints = MakeSimplePolygon(io, mode);
    int32_t numExtra = io.rawInteger(1, 4);
    std::vector<Vector2<double>> points{};
    std::vector<int32_t> polygon{};
    for (int32_t k = 0; k < numExtra; ++k)
    {
        points.push_back(Vector2<double>{
            static_cast<double>(io.rawInteger(-9, 9)),
            static_cast<double>(io.rawInteger(-9, 9)) });
    }
    for (auto const& p : polygonPoints)
    {
        polygon.push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }
    GivePoints(io, points);
    GivePolygon(io, polygon);

    TriangulateEC<double, double> triangulator(points);
    triangulator(polygon);
    OutTriangles(io, triangulator.GetTriangles());

    bool referenceOk = (triangulator.GetTriangles().size() == polygon.size() - 2)
        && TriangulationIsValid(points, triangulator.GetTriangles(),
            ShoelaceArea(points, polygon));
    io.outBool(referenceOk);
}

ORACLE_CASE("TriangulateEC.triangulateWithHole")
{
    std::vector<Vector2<double>> outerPoints = MakeOuterPolygon(io);
    double cx = static_cast<double>(io.rawInteger(-4, 4));
    double cy = static_cast<double>(io.rawInteger(-4, 4));
    std::vector<Vector2<double>> innerPoints = MakeHolePolygon(io, cx, cy);

    std::vector<Vector2<double>> points{};
    std::vector<int32_t> outer{}, inner{};
    for (auto const& p : outerPoints)
    {
        outer.push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }
    for (auto const& p : innerPoints)
    {
        inner.push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }
    GivePoints(io, points);
    GivePolygon(io, outer);
    GivePolygon(io, inner);

    TriangulateEC<double, double> triangulator(points);
    triangulator(outer, inner);
    OutTriangles(io, triangulator.GetTriangles());

    double expected = ShoelaceArea(points, outer) + ShoelaceArea(points, inner);
    io.outBool(TriangulationIsValid(points, triangulator.GetTriangles(), expected));
}

ORACLE_CASE("TriangulateEC.triangulateWithHoles")
{
    std::vector<Vector2<double>> outerPoints = MakeOuterPolygon(io);
    std::vector<Vector2<double>> hole0 = MakeHolePolygon(io, -4.0, 0.0);
    std::vector<Vector2<double>> hole1 = MakeHolePolygon(io, 4.0, 0.0);

    std::vector<Vector2<double>> points{};
    std::vector<int32_t> outer{};
    std::vector<std::vector<int32_t>> inners(2);
    for (auto const& p : outerPoints)
    {
        outer.push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }
    for (auto const& p : hole0)
    {
        inners[0].push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }
    for (auto const& p : hole1)
    {
        inners[1].push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }
    GivePoints(io, points);
    GivePolygon(io, outer);
    GivePolygon(io, inners[0]);
    GivePolygon(io, inners[1]);

    TriangulateEC<double, double> triangulator(points);
    triangulator(outer, inners);
    OutTriangles(io, triangulator.GetTriangles());

    double expected = ShoelaceArea(points, outer)
        + ShoelaceArea(points, inners[0]) + ShoelaceArea(points, inners[1]);
    io.outBool(TriangulationIsValid(points, triangulator.GetTriangles(), expected));
}

ORACLE_CASE("TriangulateEC.triangulateTree")
{
    // A three-level tree: the 12-direction outer polygon, a clockwise hole
    // about the origin at radius 5, and a counterclockwise polygon inside the
    // hole at radius at most 2*sqrt(2), which is inside the hole's inradius
    // of 5*cos(15 degrees) = 4.83.
    std::vector<Vector2<double>> outerPoints = MakeOuterPolygon(io);
    std::vector<Vector2<double>> holePoints{};
    for (int32_t k = 11; k >= 0; --k)
    {
        holePoints.push_back(Vector2<double>{
            static_cast<double>(circle5[k][0]),
            static_cast<double>(circle5[k][1]) });
    }
    std::vector<Vector2<double>> innerOuterPoints = MakeSmallStar(io, 0.0, 0.0);

    std::vector<Vector2<double>> points{};
    std::vector<int32_t> outer{}, hole{}, innerOuter{};
    for (auto const& p : outerPoints)
    {
        outer.push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }
    for (auto const& p : holePoints)
    {
        hole.push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }
    for (auto const& p : innerOuterPoints)
    {
        innerOuter.push_back(static_cast<int32_t>(points.size()));
        points.push_back(p);
    }

    GivePoints(io, points);
    GivePolygon(io, outer);
    GivePolygon(io, hole);
    GivePolygon(io, innerOuter);

    auto root = std::make_shared<PolygonTree>();
    root->polygon = outer;
    auto holeNode = std::make_shared<PolygonTree>();
    holeNode->polygon = hole;
    auto innerNode = std::make_shared<PolygonTree>();
    innerNode->polygon = innerOuter;
    holeNode->child.push_back(innerNode);
    root->child.push_back(holeNode);

    TriangulateEC<double, double> triangulator(points);
    triangulator(root);
    OutTriangles(io, triangulator.GetTriangles());

    double expected = ShoelaceArea(points, outer) + ShoelaceArea(points, hole)
        + ShoelaceArea(points, innerOuter);
    io.outBool(TriangulationIsValid(points, triangulator.GetTriangles(), expected));
}

namespace
{
    // ---- MinimumVolumeBox3 ----------------------------------------------
    //
    // The port implements the MVB3FloatingPoint and MVB3Rational compute
    // types; MVB3GPU is not ported. Both are instantiated here with
    // <double, int32_t, ...> and run single threaded (numThreads = 0), and
    // lgMaxSample stays small so that the rational pipeline finishes.
    //
    // The port deliberately fixes four upstream defects in these two files,
    // so the generators of the main cases have to stay on the sound side of
    // each of them and each gets its own deviation case:
    //   #352 / #355  the dimension-2 Newell normal loop drops the
    //                wrap-around term in both files, so every coplanar input
    //                is wrong upstream. The main cases reject coplanar point
    //                sets (hull dimension < 3) and the deviation cases aim
    //                straight at them.
    //   #405         ComputeVolume assumes a hull-edge vertex realizes the
    //                axis[0]/axis[1] minima.
    //   #426         GetExtreme's strict-improvement hill climb stalls on a
    //                floating-point plateau.
    // Both #405 and #426 show up as the same observable: the returned box
    // does not contain the input points. The main floating-point case
    // rejects such records with a capped loop that redraws the whole cloud,
    // and the deviation case keeps only the worst offender it finds.
    //
    // ORDER DEPENDENCE (new upstream finding, see the group report).
    // ExtractMeshTopology numbers the edges and the triangles by iterating
    // ETManifoldMesh's mEMap and mTMap, which are std::unordered_map. The
    // candidate list mEdgeIndices is then every ordered pair (e0, e1) with
    // e0 < e1 IN THAT NUMBERING, and ProcessEdgePair is not symmetric in its
    // two edges (the first supplies N and the second M, and the level-curve
    // branch and the sampled axes differ between the two roles). So both the
    // set of candidates examined and the tie-breaking "first candidate with
    // a strictly smaller volume wins" depend on a hash-container iteration
    // order, which is not part of any specification. Measured with the probe
    // below: reversing the triangle order of the input mesh changes the
    // reported minimum volume on a noticeable fraction of small clouds, by
    // up to 0.8 percent. The oracle therefore
    //   * emits the box in a form invariant under the axis signs and the
    //     axis order, the two freedoms the algorithm genuinely leaves open,
    //     and
    //   * restricts the main cases to inputs on which upstream's own answer
    //     does not change when the mesh is presented in a different order.

    // Maximum amount by which a point sticks out of the box, measured along
    // the box axes. Zero when the box contains every point.
    double ContainmentViolation(OrientedBox3<double> const& box,
        std::vector<Vector3<double>> const& points)
    {
        double worst = 0.0;
        for (auto const& p : points)
        {
            Vector3<double> d = p - box.center;
            for (int32_t i = 0; i < 3; ++i)
            {
                double t = std::fabs(Dot(d, box.axis[i])) - box.extent[i];
                worst = std::max(worst, t);
            }
        }
        return worst;
    }

    double PointScale(std::vector<Vector3<double>> const& points)
    {
        double scale = 1.0;
        for (auto const& p : points)
        {
            for (int32_t i = 0; i < 3; ++i)
            {
                scale = std::max(scale, std::fabs(p[i]));
            }
        }
        return scale;
    }

    // Draw a point cloud. Mode 0 is a small lattice; mode 1 rotates a lattice
    // cloud by a fraction of a degree, the configuration in which upstream's
    // hill climb stalls (#426); mode 2 is uniform. The sin/cos of mode 1 is
    // applied here, in the generator; only the final coordinates are
    // recorded, so no libm value enters a compared computation.
    std::vector<Vector3<double>> DrawCloud(oracle::Ctx& io, int32_t mode, int32_t n)
    {
        std::vector<Vector3<double>> points(static_cast<size_t>(n));
        for (int32_t i = 0; i < n; ++i)
        {
            for (int32_t d = 0; d < 3; ++d)
            {
                points[static_cast<size_t>(i)][d] = (mode == 2
                    ? io.raw(-5.0, 5.0)
                    : static_cast<double>(io.rawInteger(-4, 4)));
            }
        }
        if (mode == 1)
        {
            double angle = io.raw(0.0005, 0.005);
            double c = std::cos(angle), s = std::sin(angle);
            int32_t axis = io.rawInteger(0, 2);
            int32_t a0 = (axis + 1) % 3, a1 = (axis + 2) % 3;
            for (auto& p : points)
            {
                double u = p[a0], v = p[a1];
                p[a0] = c * u - s * v;
                p[a1] = s * u + c * v;
            }
        }
        return points;
    }

    // A coplanar cloud: 'n' lattice points on a lattice plane through
    // 'origin' spanned by 'dir0' and 'dir1'. Sub-mode 0 forces a triangular
    // hull (three corners plus interior points), sub-mode 1 allows any hull.
    std::vector<Vector3<double>> DrawCoplanarCloud(oracle::Ctx& io, int32_t subMode,
        int32_t n)
    {
        Vector3<double> origin{ 0.0, 0.0, 0.0 }, dir0{ 0.0, 0.0, 0.0 }, dir1{ 0.0, 0.0, 0.0 };
        for (int32_t d = 0; d < 3; ++d)
        {
            origin[d] = static_cast<double>(io.rawInteger(-3, 3));
        }
        do
        {
            for (int32_t d = 0; d < 3; ++d)
            {
                dir0[d] = static_cast<double>(io.rawInteger(-2, 2));
                dir1[d] = static_cast<double>(io.rawInteger(-2, 2));
            }
        }
        while (Length(Cross(dir0, dir1)) == 0.0);

        std::vector<std::array<double, 2>> uv{};
        if (subMode == 0)
        {
            double scale = static_cast<double>(io.rawInteger(2, 4));
            uv.push_back({ 0.0, 0.0 });
            uv.push_back({ scale, 0.0 });
            uv.push_back({ 0.0, scale });
            for (int32_t i = 3; i < n; ++i)
            {
                // Strictly inside the triangle with corners (0,0), (s,0) and
                // (0,s), so the hull stays a triangle.
                double a = 0.25 * static_cast<double>(io.rawInteger(1, 2));
                double b = 0.25 * static_cast<double>(io.rawInteger(1, 2));
                uv.push_back({ a * scale, b * scale });
            }
        }
        else
        {
            for (int32_t i = 0; i < n; ++i)
            {
                uv.push_back({ static_cast<double>(io.rawInteger(-4, 4)),
                    static_cast<double>(io.rawInteger(-4, 4)) });
            }
        }

        std::vector<Vector3<double>> points{};
        for (auto const& p : uv)
        {
            points.push_back(origin + p[0] * dir0 + p[1] * dir1);
        }
        return points;
    }

    void GiveCloud(oracle::Ctx& io, std::vector<Vector3<double>> const& points)
    {
        int32_t n = static_cast<int32_t>(points.size());
        io.integer(n, n);
        for (auto const& p : points)
        {
            io.givenVec(p);
        }
    }

    // The box is emitted in a form that is invariant under the two freedoms
    // the algorithm leaves unspecified: the sign of each axis and the order
    // of the three axes. See the comment on the order dependence above. The
    // invariants are the centre, the sorted extents, the six distinct
    // entries of M = sum_i extent[i]^2 * axis[i] * axis[i]^T, and the volume;
    // together they determine the box up to axis signs and permutation.
    void OutBox(oracle::Ctx& io, size_t dimension, OrientedBox3<double> const& box,
        double volume)
    {
        io.outInt(static_cast<int32_t>(dimension));
        io.outVec(box.center);

        std::array<double, 3> extent{ box.extent[0], box.extent[1], box.extent[2] };
        std::sort(extent.begin(), extent.end());
        for (int32_t i = 0; i < 3; ++i)
        {
            io.outReal(extent[static_cast<size_t>(i)]);
        }

        std::array<double, 6> m{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        for (int32_t i = 0; i < 3; ++i)
        {
            double w = box.extent[i] * box.extent[i];
            Vector3<double> const& a = box.axis[i];
            m[0] += w * a[0] * a[0];
            m[1] += w * a[0] * a[1];
            m[2] += w * a[0] * a[2];
            m[3] += w * a[1] * a[1];
            m[4] += w * a[1] * a[2];
            m[5] += w * a[2] * a[2];
        }
        for (int32_t i = 0; i < 6; ++i)
        {
            io.outReal(m[static_cast<size_t>(i)]);
        }

        io.outReal(volume);
    }

    // The invariants of a box, for the order-stability probe.
    std::array<double, 13> BoxInvariants(OrientedBox3<double> const& box, double volume)
    {
        std::array<double, 13> out{};
        out[0] = box.center[0];
        out[1] = box.center[1];
        out[2] = box.center[2];
        std::array<double, 3> extent{ box.extent[0], box.extent[1], box.extent[2] };
        std::sort(extent.begin(), extent.end());
        out[3] = extent[0];
        out[4] = extent[1];
        out[5] = extent[2];
        for (int32_t i = 0; i < 3; ++i)
        {
            double w = box.extent[i] * box.extent[i];
            Vector3<double> const& a = box.axis[i];
            out[6] += w * a[0] * a[0];
            out[7] += w * a[0] * a[1];
            out[8] += w * a[0] * a[2];
            out[9] += w * a[1] * a[1];
            out[10] += w * a[1] * a[2];
            out[11] += w * a[2] * a[2];
        }
        out[12] = volume;
        return out;
    }

    bool InvariantsAgree(std::array<double, 13> const& a, std::array<double, 13> const& b,
        double tolerance)
    {
        for (size_t i = 0; i < a.size(); ++i)
        {
            double scale = std::max(1.0, std::max(std::fabs(a[i]), std::fabs(b[i])));
            if (std::fabs(a[i] - b[i]) > tolerance * scale)
            {
                return false;
            }
        }
        return true;
    }

    using MVB3FP = MinimumVolumeBox3<double, int32_t, MVB3FloatingPoint>;
    using MVB3R = MinimumVolumeBox3<double, int32_t, MVB3Rational>;

    // Upstream's MinimizerVariableT declares its three parameters as
    // "T const&" where the sibling MinimizerVariableS correctly uses
    // "Number const&" (finding #355). Every caller passes BSNumber
    // expressions, so each one is silently rounded to double and the whole
    // t-variable level curve is sampled at rounded parameters; the port
    // samples them exactly. This subclass changes nothing but records whether
    // the defective routine was reached, which is the exact separator between
    // the inputs on which upstream's rational pipeline is exact and the ones
    // on which it is not.
    class MVB3RProbe : public MVB3R
    {
    public:
        MVB3RProbe() : MVB3R(0), variableTCalls(0) {}

        std::size_t variableTCalls;

    protected:
        void MinimizerVariableT(double const& tminNumer, double const& tmaxNumer,
            double const& tDenom, Candidate& c, Candidate& mvc) override
        {
            ++variableTCalls;
            MVB3R::MinimizerVariableT(tminNumer, tmaxNumer, tDenom, c, mvc);
        }
    };

    void MVB3FPQuery(std::vector<Vector3<double>> const& vertices,
        std::vector<int32_t> const& indices, int32_t lgMaxSample,
        OrientedBox3<double>& box, double& volume)
    {
        MVB3FP query(0);
        query(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);
    }

    // Run the query twice, the second time with the triangles of the hull
    // presented in reverse order, and report whether the two boxes agree.
    // This is upstream's own control flow used as the probe.
    bool HullOrderStable(std::vector<Vector3<double>> const& vertices,
        std::vector<int32_t> const& indices, int32_t lgMaxSample)
    {
        size_t const numTriangles = indices.size() / 3;
        std::vector<std::vector<int32_t>> orders{};
        orders.push_back(indices);
        std::vector<int32_t> reversed{};
        for (size_t t = numTriangles; t > 0; --t)
        {
            reversed.push_back(indices[3 * (t - 1) + 0]);
            reversed.push_back(indices[3 * (t - 1) + 1]);
            reversed.push_back(indices[3 * (t - 1) + 2]);
        }
        orders.push_back(reversed);
        for (size_t shift = 1; shift < numTriangles; ++shift)
        {
            std::vector<int32_t> rotated{};
            for (size_t t = 0; t < numTriangles; ++t)
            {
                size_t u = (t + shift) % numTriangles;
                rotated.push_back(indices[3 * u + 0]);
                rotated.push_back(indices[3 * u + 1]);
                rotated.push_back(indices[3 * u + 2]);
            }
            orders.push_back(rotated);
        }

        std::array<double, 13> reference{};
        for (size_t k = 0; k < orders.size(); ++k)
        {
            OrientedBox3<double> box{};
            double volume = 0.0;
            try
            {
                MVB3FPQuery(vertices, orders[k], lgMaxSample, box, volume);
            }
            catch (std::exception const&)
            {
                return false;
            }
            std::array<double, 13> invariants = BoxInvariants(box, volume);
            if (k == 0)
            {
                reference = invariants;
            }
            else if (!InvariantsAgree(reference, invariants, 0.0))
            {
                return false;
            }
        }
        return true;
    }

    // The same probe for the arbitrary-point query: run it on the points as
    // drawn and on the reversed point list. Reversing changes both the
    // ConvexHull3 output ordering and the mesh insertion order, and it moves
    // the translation origin mTOrigin, so bit-identical invariants are a
    // strong statement that the answer does not depend on either.
    bool CloudOrderStable(std::vector<Vector3<double>> const& points, int32_t lgMaxSample)
    {
        std::vector<std::vector<Vector3<double>>> orders{};
        orders.push_back(points);
        orders.push_back(std::vector<Vector3<double>>(points.rbegin(), points.rend()));
        for (std::size_t shift = 1; shift < points.size(); ++shift)
        {
            std::vector<Vector3<double>> rotated{};
            for (std::size_t i = 0; i < points.size(); ++i)
            {
                rotated.push_back(points[(i + shift) % points.size()]);
            }
            orders.push_back(rotated);
        }

        std::array<double, 13> reference{};
        for (std::size_t k = 0; k < orders.size(); ++k)
        {
            OrientedBox3<double> box{};
            double volume = 0.0;
            std::size_t dimension = 0;
            try
            {
                MVB3FP query(0);
                dimension = query(orders[k], static_cast<std::size_t>(lgMaxSample),
                    box, volume);
            }
            catch (std::exception const&)
            {
                return false;
            }
            if (dimension != 3)
            {
                return false;
            }
            std::array<double, 13> invariants = BoxInvariants(box, volume);
            if (k == 0)
            {
                reference = invariants;
            }
            else if (!InvariantsAgree(reference, invariants, 0.0))
            {
                return false;
            }
        }
        return true;
    }

}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.compute")
{
    RunWithBigStack([&io]() {

    // Restricted to clouds whose hull is 3-dimensional (the dimension-0, -1
    // and -2 paths get their own cases) and on which upstream's box really
    // contains the input points, which is the observable symptom of both
    // #405 and #426. The loop is capped at 24 attempts and falls back to the
    // least defective cloud it saw.
    int32_t mode = io.integer(0, 2);
    int32_t lgMaxSample = io.integer(2, 4);
    int32_t n = io.integer(5, 10);

    std::vector<Vector3<double>> best{};
    double bestViolation = std::numeric_limits<double>::max();
    for (int32_t attempt = 0; attempt < 24; ++attempt)
    {
        std::vector<Vector3<double>> points = DrawCloud(io, mode, n);
        OrientedBox3<double> box{};
        double volume = 0.0;
        MVB3FP query(0);
        std::size_t dimension = 3;
        try
        {
            dimension = query(points, static_cast<std::size_t>(lgMaxSample), box, volume);
        }
        catch (std::exception const&)
        {
            continue;
        }
        if (dimension != 3)
        {
            continue;
        }
        double violation = ContainmentViolation(box, points) / PointScale(points);
        if (!CloudOrderStable(points, lgMaxSample))
        {
            // The reported box depends on the order in which the mesh edges
            // and triangles happen to be enumerated, which comes out of a
            // std::unordered_map; such a record is not comparable.
            continue;
        }
        if (violation < bestViolation)
        {
            bestViolation = violation;
            best = points;
        }
        if (violation <= 1e-9)
        {
            break;
        }
    }
    if (best.empty())
    {
        best = { Vector3<double>{ 0.0, 0.0, 0.0 }, Vector3<double>{ 4.0, 0.0, 0.0 },
            Vector3<double>{ 0.0, 3.0, 0.0 }, Vector3<double>{ 0.0, 0.0, 2.0 },
            Vector3<double>{ 1.0, 1.0, 1.0 } };
    }
    GiveCloud(io, best);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3FP query(0);
    std::size_t dimension = query(best, static_cast<std::size_t>(lgMaxSample), box, volume);

    // Only the hull dimension, the minimum volume and the reference checks
    // are compared. WHICH candidate attains the minimum depends on the
    // std::unordered_map enumeration order of the mesh edges and triangles
    // (see the section comment), and for a point cloud the hull itself is
    // produced by ConvexHull3, whose vertex ordering is a second
    // container-dependent input to the same search. The sibling computeHull
    // case, which fixes the mesh, compares the full box bit for bit.
    io.outInt(static_cast<int32_t>(dimension));
    io.outReal(volume);

    // Reference check: the box contains every input point and its volume is
    // the product of twice the extents.
    double scale = PointScale(best);
    io.outBool(ContainmentViolation(box, best) <= 1e-9 * scale);
    double product = 8.0 * box.extent[0] * box.extent[1] * box.extent[2];
    io.outBool(std::fabs(product - volume) <= 1e-9 * std::max(1.0, std::fabs(volume)));
    });
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.compute.lowDimension")
{
    RunWithBigStack([&io]() {

    // Hull dimension 0 (all points equal) and 1 (all points collinear). The
    // Newell defect of #352 is confined to the dimension-2 branch, so these
    // two agree with the port.
    int32_t dimensionWanted = io.integer(0, 1);
    int32_t lgMaxSample = io.integer(2, 4);
    int32_t n = io.integer(4, 8);

    std::vector<Vector3<double>> points{};
    Vector3<double> origin{ 0.0, 0.0, 0.0 }, dir{ 0.0, 0.0, 0.0 };
    for (int32_t d = 0; d < 3; ++d)
    {
        origin[d] = static_cast<double>(io.rawInteger(-4, 4));
    }
    if (dimensionWanted == 0)
    {
        for (int32_t i = 0; i < n; ++i)
        {
            points.push_back(origin);
        }
    }
    else
    {
        do
        {
            for (int32_t d = 0; d < 3; ++d)
            {
                dir[d] = static_cast<double>(io.rawInteger(-3, 3));
            }
        }
        while (dir[0] == 0.0 && dir[1] == 0.0 && dir[2] == 0.0);
        for (int32_t i = 0; i < n; ++i)
        {
            double t = static_cast<double>(io.rawInteger(-4, 4));
            points.push_back(origin + t * dir);
        }
    }
    GiveCloud(io, points);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3FP query(0);
    std::size_t dimension = query(points, static_cast<std::size_t>(lgMaxSample), box, volume);
    OutBox(io, dimension, box, volume);
    });
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.computeHull")
{
    RunWithBigStack([&io]() {

    // The vertices-and-indices query, fed the three convex lattice polytopes
    // (the cube's coplanar face triangles exercise
    // RemoveCoplanarTriangleAdjacencies).
    int32_t shape = io.integer(0, 2);
    int32_t lgMaxSample = io.integer(2, 4);
    std::vector<Vector3<double>> vertices{};
    std::vector<int32_t> indices{};
    MakePolytope(io, shape, vertices, indices);
    GiveMesh(io, vertices, indices);

    // The three polytope shapes are fixed, so a rejection loop could not
    // redraw its way out of a defect. Instead the two acceptance tests are
    // evaluated here and RECORDED AS AN INPUT, so the replay takes the same
    // branch: the box is compared only when upstream's answer does not
    // depend on the mesh order and upstream's box contains the polytope.
    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3FP query(0);
    query(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);
    double scale = PointScale(vertices);
    bool contains = (ContainmentViolation(box, vertices) <= 1e-9 * scale);
    bool stable = HullOrderStable(vertices, indices, lgMaxSample);
    int32_t usable = ((stable && contains) ? 1 : 0);
    io.integer(usable, usable);
    if (usable != 0)
    {
        OutBox(io, 3, box, volume);
        io.outBool(ContainmentViolation(box, vertices) <= 1e-9 * scale);
    }
    });
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.compute.coplanar")
{
    RunWithBigStack([&io]() {

    // Deviation: the dimension-2 Newell normal loop
    // "for (i0 = numHull - 1, i1 = 1; i1 < numHull; i0 = i1++)" drops the
    // wrap-around cross-product term, which makes the plane normal wrong for
    // every coplanar point set and exactly zero when the hull is a triangle.
    // The port starts the loop at i1 = 0. Sub-mode 0 forces a triangular
    // hull; sub-mode 1 allows any coplanar hull.
    int32_t subMode = io.integer(0, 1);
    int32_t lgMaxSample = io.integer(2, 4);
    int32_t n = io.integer(4, 8);
    std::vector<Vector3<double>> points = DrawCoplanarCloud(io, subMode, n);
    GiveCloud(io, points);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3FP query(0);
    std::size_t dimension = query(points, static_cast<std::size_t>(lgMaxSample), box, volume);
    OutBox(io, dimension, box, volume);
    });
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.compute.nonContaining")
{
    RunWithBigStack([&io]() {

    // Deviation: clouds on which upstream's own box fails to contain the
    // input points, the observable symptom of #405 (ComputeVolume assumes a
    // hull-edge vertex realizes the axis[0]/axis[1] minima) and #426
    // (GetExtreme's strict-improvement hill climb stalls on a plateau). The
    // loop is capped at 48 attempts and keeps the worst offender it saw; the
    // port's box contains the points, so the two disagree.
    int32_t mode = io.integer(0, 1);
    int32_t lgMaxSample = io.integer(2, 4);
    int32_t n = io.integer(6, 10);

    std::vector<Vector3<double>> best{};
    double bestViolation = -1.0;
    for (int32_t attempt = 0; attempt < 48; ++attempt)
    {
        std::vector<Vector3<double>> points = DrawCloud(io, mode, n);
        OrientedBox3<double> box{};
        double volume = 0.0;
        MVB3FP query(0);
        std::size_t dimension = 3;
        try
        {
            dimension = query(points, static_cast<std::size_t>(lgMaxSample), box, volume);
        }
        catch (std::exception const&)
        {
            continue;
        }
        if (dimension != 3)
        {
            continue;
        }
        double violation = ContainmentViolation(box, points) / PointScale(points);
        if (violation > bestViolation)
        {
            bestViolation = violation;
            best = points;
        }
        if (violation > 1e-3)
        {
            break;
        }
    }
    if (best.empty())
    {
        best = { Vector3<double>{ 0.0, 0.0, 0.0 }, Vector3<double>{ 4.0, 0.0, 0.0 },
            Vector3<double>{ 0.0, 3.0, 0.0 }, Vector3<double>{ 0.0, 0.0, 2.0 },
            Vector3<double>{ 1.0, 1.0, 1.0 } };
    }
    GiveCloud(io, best);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3FP query(0);
    std::size_t dimension = query(best, static_cast<std::size_t>(lgMaxSample), box, volume);
    OutBox(io, dimension, box, volume);
    io.outBool(ContainmentViolation(box, best) <= 1e-9 * PointScale(best));
    });
}


namespace
{
    // The rational pipeline. Its arithmetic is exact, so the minimum over the
    // sampled candidates is a well-defined number and does not depend on the
    // enumeration order; WHICH candidate attains it still does, so the same
    // invariant output form and the same order-stability probe are used. The
    // rational query is much slower, so lgMaxSample stays at 2, the clouds
    // stay tiny and the probe uses two alternative orders instead of all of
    // them.
    bool RationalCloudOrderStable(std::vector<Vector3<double>> const& points,
        int32_t lgMaxSample)
    {
        std::vector<std::vector<Vector3<double>>> orders{};
        orders.push_back(points);
        orders.push_back(std::vector<Vector3<double>>(points.rbegin(), points.rend()));
        std::vector<Vector3<double>> rotated{};
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            rotated.push_back(points[(i + 1) % points.size()]);
        }
        orders.push_back(rotated);

        std::array<double, 13> reference{};
        for (std::size_t k = 0; k < orders.size(); ++k)
        {
            OrientedBox3<double> box{};
            double volume = 0.0;
            std::size_t dimension = 0;
            try
            {
                MVB3R query(0);
                dimension = query(orders[k], static_cast<std::size_t>(lgMaxSample),
                    box, volume);
            }
            catch (std::exception const&)
            {
                return false;
            }
            if (dimension != 3)
            {
                return false;
            }
            std::array<double, 13> invariants = BoxInvariants(box, volume);
            if (k == 0)
            {
                reference = invariants;
            }
            else if (!InvariantsAgree(reference, invariants, 0.0))
            {
                return false;
            }
        }
        return true;
    }
}

ORACLE_CASE("MinimumVolumeBox3Rational.compute")
{
    RunWithBigStack([&io]() {
    // Uniform clouds only. The rational pipeline compares candidate volumes
    // EXACTLY, so on a lattice cloud many candidates tie exactly and which
    // one wins is decided by the std::unordered_map enumeration order; on a
    // uniform cloud an exact tie has probability zero and the minimizing
    // candidate is unique. Lattice clouds are not comparable here and are
    // reported under "Not covered".
    int32_t mode = io.integer(2, 2);
    int32_t lgMaxSample = io.integer(2, 2);
    int32_t n = io.integer(5, 7);

    // Accept only clouds on which upstream's rational pipeline is exact,
    // i.e. on which the defective MinimizerVariableT is never reached, AND on
    // which the answer does not change with the input order.
    std::vector<Vector3<double>> best{};
    for (int32_t attempt = 0; attempt < 6 && best.empty(); ++attempt)
    {
        std::vector<Vector3<double>> points = DrawCloud(io, mode, n);
        OrientedBox3<double> probeBox{};
        double probeVolume = 0.0;
        MVB3RProbe probe{};
        try
        {
            probe(points, static_cast<std::size_t>(lgMaxSample), probeBox, probeVolume);
        }
        catch (std::exception const&)
        {
            continue;
        }
        if (probe.variableTCalls == 0 && RationalCloudOrderStable(points, lgMaxSample))
        {
            best = points;
        }
    }
    if (best.empty())
    {
        best = { Vector3<double>{ 0.0, 0.0, 0.0 }, Vector3<double>{ 4.0, 0.0, 0.0 },
            Vector3<double>{ 0.0, 3.0, 0.0 }, Vector3<double>{ 0.0, 0.0, 2.0 },
            Vector3<double>{ 1.0, 1.0, 1.0 } };
    }
    GiveCloud(io, best);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3R query(0);
    std::size_t dimension = query(best, static_cast<std::size_t>(lgMaxSample), box, volume);

    // Only the hull dimension, the minimum volume and the containment of the
    // input points are compared here. WHICH candidate attains the minimum is
    // not comparable: it depends on the std::unordered_map enumeration order
    // (see the section comment) and, additionally, upstream's
    // MinimizerVariableT declares its parameters as T rather than Number, so
    // every t-variable level curve is sampled at parameters rounded to
    // double while the port samples them exactly (finding #355). The
    // sibling computeHull case compares the full box on a fixed mesh.
    io.outInt(static_cast<int32_t>(dimension));
    io.outReal(volume);
    double scale = PointScale(best);
    io.outBool(ContainmentViolation(box, best) <= 1e-9 * scale);
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.compute.lowDimension")
{
    RunWithBigStack([&io]() {
    int32_t dimensionWanted = io.integer(0, 1);
    int32_t lgMaxSample = io.integer(2, 2);
    int32_t n = io.integer(4, 6);

    std::vector<Vector3<double>> points{};
    Vector3<double> origin{ 0.0, 0.0, 0.0 }, dir{ 0.0, 0.0, 0.0 };
    for (int32_t d = 0; d < 3; ++d)
    {
        origin[d] = static_cast<double>(io.rawInteger(-4, 4));
    }
    if (dimensionWanted == 0)
    {
        for (int32_t i = 0; i < n; ++i)
        {
            points.push_back(origin);
        }
    }
    else
    {
        do
        {
            for (int32_t d = 0; d < 3; ++d)
            {
                dir[d] = static_cast<double>(io.rawInteger(-3, 3));
            }
        }
        while (dir[0] == 0.0 && dir[1] == 0.0 && dir[2] == 0.0);
        for (int32_t i = 0; i < n; ++i)
        {
            double t = static_cast<double>(io.rawInteger(-4, 4));
            points.push_back(origin + t * dir);
        }
    }
    GiveCloud(io, points);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3R query(0);
    std::size_t dimension = query(points, static_cast<std::size_t>(lgMaxSample), box, volume);
    OutBox(io, dimension, box, volume);
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.computeHull")
{
    RunWithBigStack([&io]() {
    int32_t shape = io.integer(0, 2);
    int32_t lgMaxSample = io.integer(2, 2);
    std::vector<Vector3<double>> vertices{};
    std::vector<int32_t> indices{};
    MakePolytope(io, shape, vertices, indices);
    GiveMesh(io, vertices, indices);

    // The order-stability probe, recorded as an input so that the replay
    // takes the same branch. See the floating-point sibling.
    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3R query(0);
    query(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);

    size_t const numTriangles = indices.size() / 3;
    bool stable = true;
    std::vector<size_t> shifts{ 1, numTriangles / 2, numTriangles - 1 };
    for (size_t si = 0; si < shifts.size() && stable; ++si)
    {
        size_t shift = shifts[si];
        if (shift == 0)
        {
            continue;
        }
        std::vector<int32_t> rotated{};
        for (size_t t = 0; t < numTriangles; ++t)
        {
            size_t u = (t + shift) % numTriangles;
            rotated.push_back(indices[3 * u + 0]);
            rotated.push_back(indices[3 * u + 1]);
            rotated.push_back(indices[3 * u + 2]);
        }
        OrientedBox3<double> boxOther{};
        double volumeOther = 0.0;
        MVB3R queryOther(0);
        queryOther(vertices, rotated, static_cast<std::size_t>(lgMaxSample),
            boxOther, volumeOther);
        stable = InvariantsAgree(BoxInvariants(box, volume),
            BoxInvariants(boxOther, volumeOther), 0.0);
    }

    double scale = PointScale(vertices);
    bool contains = (ContainmentViolation(box, vertices) <= 1e-9 * scale);
    int32_t usable = ((stable && contains) ? 1 : 0);
    io.integer(usable, usable);
    if (usable != 0)
    {
        OutBox(io, 3, box, volume);
        io.outBool(contains);
    }
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.compute.variableT")
{
    RunWithBigStack([&io]() {
    // Deviation: clouds on which upstream reaches MinimizerVariableT, whose
    // tminNumer, tmaxNumer and tDenom parameters are declared "T const&"
    // instead of "Number const&", so the exact BSNumber expressions the
    // callers pass are rounded to double and the whole t-variable level curve
    // is sampled at rounded parameters (finding #355). The port samples them
    // exactly, so a different candidate can win and the reported minimum
    // volume differs. The loop is capped and falls back to the cloud with the
    // most calls seen.
    int32_t lgMaxSample = io.integer(2, 2);
    int32_t n = io.integer(5, 7);

    std::vector<Vector3<double>> best{};
    std::size_t bestCalls = 0;
    for (int32_t attempt = 0; attempt < 12; ++attempt)
    {
        std::vector<Vector3<double>> points = DrawCloud(io, 2, n);
        OrientedBox3<double> probeBox{};
        double probeVolume = 0.0;
        MVB3RProbe probe{};
        std::size_t dimension = 0;
        try
        {
            dimension = probe(points, static_cast<std::size_t>(lgMaxSample),
                probeBox, probeVolume);
        }
        catch (std::exception const&)
        {
            continue;
        }
        if (dimension == 3 && probe.variableTCalls > bestCalls)
        {
            bestCalls = probe.variableTCalls;
            best = points;
        }
        if (bestCalls > 0)
        {
            break;
        }
    }
    if (best.empty())
    {
        best = { Vector3<double>{ 0.0, 0.0, 0.0 }, Vector3<double>{ 4.0, 0.0, 0.0 },
            Vector3<double>{ 0.0, 3.0, 0.0 }, Vector3<double>{ 0.0, 0.0, 2.0 },
            Vector3<double>{ 1.0, 1.0, 1.0 } };
    }
    GiveCloud(io, best);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3R query(0);
    std::size_t dimension = query(best, static_cast<std::size_t>(lgMaxSample), box, volume);
    io.outInt(static_cast<int32_t>(dimension));
    io.outReal(volume);
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.compute.coplanar")
{
    RunWithBigStack([&io]() {
    // Deviation: the same dimension-2 Newell wrap-around omission as the
    // floating-point sibling, in the rational file (finding #355).
    int32_t subMode = io.integer(0, 1);
    int32_t lgMaxSample = io.integer(2, 2);
    int32_t n = io.integer(4, 7);
    std::vector<Vector3<double>> points = DrawCoplanarCloud(io, subMode, n);
    GiveCloud(io, points);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MVB3R query(0);
    std::size_t dimension = query(points, static_cast<std::size_t>(lgMaxSample), box, volume);
    OutBox(io, dimension, box, volume);
    });
}
