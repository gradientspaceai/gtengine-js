// Verify group 30 (intersection): differential cases for the Intr* headers
// listed in plan/verify-groups.json group 30.
//
// Every generator mixes an exact small-lattice mode with a uniform mode (and,
// where a branch needs it, a constructed mode) so that the touching, parallel,
// coincident and contained branches are reached. All of these queries use only
// + - * / sqrt fabs min max and comparisons, so every case is declared exact
// on the TypeScript side; the sin/cos used to build orthonormal frames are
// applied to *unrecorded* draws and only the resulting frame is recorded as an
// input, so libm never enters the compared computation.
#define ORACLE_FAMILY "v30-intersection"
#include "Oracle.h"

#include <Mathematics/IntrIntervals.h>
#include <Mathematics/IntrAlignedBox2AlignedBox2.h>
#include <Mathematics/IntrAlignedBox2OrientedBox2.h>
#include <Mathematics/IntrAlignedBox3AlignedBox3.h>
#include <Mathematics/IntrAlignedBox3OrientedBox3.h>
#include <Mathematics/IntrCircle2Circle2.h>
#include <Mathematics/IntrHalfspace2Polygon2.h>
#include <Mathematics/IntrHalfspace3OrientedBox3.h>
#include <Mathematics/IntrHalfspace3Segment3.h>
#include <Mathematics/IntrHalfspace3Sphere3.h>
#include <Mathematics/IntrHalfspace3Triangle3.h>
#include <Mathematics/IntrLine2AlignedBox2.h>
#include <Mathematics/IntrLine2Line2.h>
#include <Mathematics/IntrLine2Triangle2.h>
#include <Mathematics/IntrLine3AlignedBox3.h>
#include <Mathematics/IntrLine3Rectangle3.h>
#include <Mathematics/IntrLine3Sphere3.h>
#include <Mathematics/IntrLine3Triangle3.h>
#include <Mathematics/IntrOrientedBox2OrientedBox2.h>
#include <Mathematics/IntrOrientedBox3OrientedBox3.h>
#include <Mathematics/IntrRay3Triangle3.h>
#include <Mathematics/IntrSegment3Triangle3.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace gte;

namespace
{
    // ---- shared generators -------------------------------------------------
    // Each helper records the same number of inputs in every mode, so the
    // TypeScript replay reads a fixed layout.

    template <int N>
    Vector<N, double> Pt(oracle::Ctx& io, int mode, int lat, double range)
    {
        return (mode == 0 ? io.latticeVec<N>(-lat, lat) : io.vec<N>(-range, range));
    }

    // A positive extent vector: lattice mode gives exactly representable
    // extents so that touching configurations occur.
    template <int N>
    Vector<N, double> Extent(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> e{};
        for (int i = 0; i < N; ++i)
        {
            e[i] = (mode == 0 ? static_cast<double>(io.integer(1, lat))
                : io.real(0.25, range));
        }
        return e;
    }

    // Axis-aligned box with min <= max, componentwise.
    template <int N>
    AlignedBox<N, double> ABox(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> lo{}, hi{};
        for (int i = 0; i < N; ++i)
        {
            double a{}, b{};
            if (mode == 0)
            {
                a = static_cast<double>(io.rawInteger(-lat, lat));
                b = static_cast<double>(io.rawInteger(-lat, lat));
            }
            else
            {
                a = io.raw(-range, range);
                b = io.raw(-range, range);
            }
            lo[i] = std::min(a, b);
            hi[i] = std::max(a, b);
        }
        AlignedBox<N, double> box{};
        box.min = io.givenVec(lo);
        box.max = io.givenVec(hi);
        return box;
    }

    // Orthonormal 2D frame. Mode 0 produces exact coordinate axes (with signs
    // and a swap) so that the parallel-axis branches of the separating-axis
    // tests are reached with exact arithmetic.
    std::array<Vector2<double>, 2> Frame2(oracle::Ctx& io, int mode)
    {
        double c{}, s{};
        if (mode == 0)
        {
            int k = io.rawInteger(0, 3);
            c = (k == 0 ? 1.0 : (k == 1 ? 0.0 : (k == 2 ? -1.0 : 0.0)));
            s = (k == 0 ? 0.0 : (k == 1 ? 1.0 : (k == 2 ? 0.0 : -1.0)));
        }
        else
        {
            double angle = io.raw(-3.141592653589793, 3.141592653589793);
            c = std::cos(angle);
            s = std::sin(angle);
        }
        std::array<Vector2<double>, 2> axis{};
        axis[0] = io.givenVec(Vector2<double>{ c, s });
        axis[1] = io.givenVec(Vector2<double>{ -s, c });
        return axis;
    }

    // Orthonormal 3D frame. Mode 0 produces a signed permutation of the
    // coordinate axes (exact), which makes box axes parallel.
    std::array<Vector3<double>, 3> Frame3(oracle::Ctx& io, int mode)
    {
        std::array<Vector3<double>, 3> v{};
        if (mode == 0)
        {
            int const perm[6][3] =
            {
                { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 },
                { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 }
            };
            int p = io.rawInteger(0, 5);
            for (int i = 0; i < 3; ++i)
            {
                double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
                v[i].MakeZero();
                v[i][perm[p][i]] = sign;
            }
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { v[0][i] = io.raw(-1.0, 1.0); }
                len = Length(v[0]);
            } while (len < 0.1 || len > 1.0);
            Normalize(v[0]);
            ComputeOrthogonalComplement(1, v.data());
        }
        std::array<Vector3<double>, 3> axis{};
        axis[0] = io.givenVec(v[0]);
        axis[1] = io.givenVec(v[1]);
        axis[2] = io.givenVec(v[2]);
        return axis;
    }

    template <int N>
    OrientedBox<N, double> OBox(oracle::Ctx& io, int mode, int lat, double range);

    template <>
    OrientedBox<2, double> OBox<2>(oracle::Ctx& io, int mode, int lat, double range)
    {
        OrientedBox2<double> box{};
        box.center = Pt<2>(io, mode, lat, range);
        auto axis = Frame2(io, mode);
        box.axis[0] = axis[0];
        box.axis[1] = axis[1];
        box.extent = Extent<2>(io, mode, lat, range);
        return box;
    }

    template <>
    OrientedBox<3, double> OBox<3>(oracle::Ctx& io, int mode, int lat, double range)
    {
        OrientedBox3<double> box{};
        box.center = Pt<3>(io, mode, lat, range);
        auto axis = Frame3(io, mode);
        box.axis[0] = axis[0];
        box.axis[1] = axis[1];
        box.axis[2] = axis[2];
        box.extent = Extent<3>(io, mode, lat, range);
        return box;
    }

    // A halfspace whose boundary sometimes passes exactly through 'onPoint',
    // which reaches the "distance is exactly zero" branches.
    template <int N>
    Halfspace<N, double> Space(oracle::Ctx& io, int mode, Vector<N, double> const& onPoint)
    {
        Vector<N, double> normal{};
        if (mode == 0)
        {
            normal = io.latticeDir<N>(-2, 2);
        }
        else
        {
            normal = io.unit<N>();
        }
        double constant = (mode == 0 ? io.given(Dot(normal, onPoint))
            : io.real(-3.0, 3.0));
        Halfspace<N, double> halfspace{};
        halfspace.normal = normal;
        halfspace.constant = constant;
        return halfspace;
    }

    // A triangle; lattice mode makes vertices land exactly on lines/planes.
    template <int N>
    Triangle<N, double> Tri(oracle::Ctx& io, int mode, int lat, double range)
    {
        Triangle<N, double> triangle{};
        triangle.v[0] = Pt<N>(io, mode, lat, range);
        triangle.v[1] = Pt<N>(io, mode, lat, range);
        triangle.v[2] = Pt<N>(io, mode, lat, range);
        return triangle;
    }
}

// ============================ IntrIntervals ==============================

namespace
{
    // Two finite intervals with t0 <= t1. Lattice mode produces touching and
    // coincident intervals.
    std::array<double, 2> Interval(oracle::Ctx& io, int mode)
    {
        double a{}, b{};
        if (mode == 0)
        {
            a = static_cast<double>(io.rawInteger(-4, 4));
            b = static_cast<double>(io.rawInteger(-4, 4));
        }
        else
        {
            a = io.raw(-5.0, 5.0);
            b = io.raw(-5.0, 5.0);
        }
        std::array<double, 2> interval{ std::min(a, b), std::max(a, b) };
        io.given(interval[0]);
        io.given(interval[1]);
        return interval;
    }
}

ORACLE_CASE("IntrIntervals.test.finiteFinite")
{
    int mode = io.index() % 2;
    auto interval0 = Interval(io, mode);
    auto interval1 = Interval(io, mode);
    TIIntervalInterval<double> query;
    auto r = query(interval0, interval1);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrIntervals.test.finiteSemiInfinite")
{
    int mode = io.index() % 2;
    auto finite = Interval(io, mode);
    double a = (mode == 0 ? io.lattice(-4, 4) : io.real(-5.0, 5.0));
    bool isPositiveInfinite = io.boolean();
    TIIntervalInterval<double> query;
    auto r = query(finite, a, isPositiveInfinite);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrIntervals.test.semiInfiniteSemiInfinite")
{
    int mode = io.index() % 2;
    double a0 = (mode == 0 ? io.lattice(-4, 4) : io.real(-5.0, 5.0));
    bool p0 = io.boolean();
    double a1 = (mode == 0 ? io.lattice(-4, 4) : io.real(-5.0, 5.0));
    bool p1 = io.boolean();
    TIIntervalInterval<double> query;
    auto r = query(a0, p0, a1, p1);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrIntervals.test.dynamic")
{
    int mode = io.index() % 2;
    double maxTime = io.real(0.5, 4.0);
    auto interval0 = Interval(io, mode);
    double speed0 = (mode == 0 ? io.lattice(-3, 3) : io.real(-4.0, 4.0));
    auto interval1 = Interval(io, mode);
    double speed1 = (mode == 0 ? io.lattice(-3, 3) : io.real(-4.0, 4.0));
    TIIntervalInterval<double> query;
    auto r = query(maxTime, interval0, speed0, interval1, speed1);
    io.outBool(r.intersect);
    io.outReal(r.firstTime);
    io.outReal(r.lastTime);
}

ORACLE_CASE("IntrIntervals.find.finiteFinite")
{
    int mode = io.index() % 2;
    auto interval0 = Interval(io, mode);
    auto interval1 = Interval(io, mode);
    FIIntervalInterval<double> query;
    auto r = query(interval0, interval1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outInt(r.type);
    io.outReal(r.overlap[0]);
    io.outReal(r.overlap[1]);
}

ORACLE_CASE("IntrIntervals.find.finiteSemiInfinite")
{
    int mode = io.index() % 2;
    auto finite = Interval(io, mode);
    double a = (mode == 0 ? io.lattice(-4, 4) : io.real(-5.0, 5.0));
    bool isPositiveInfinite = io.boolean();
    FIIntervalInterval<double> query;
    auto r = query(finite, a, isPositiveInfinite);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outInt(r.type);
    io.outReal(r.overlap[0]);
    io.outReal(r.overlap[1]);
}

ORACLE_CASE("IntrIntervals.find.semiInfiniteSemiInfinite")
{
    int mode = io.index() % 2;
    double a0 = (mode == 0 ? io.lattice(-4, 4) : io.real(-5.0, 5.0));
    bool p0 = io.boolean();
    double a1 = (mode == 0 ? io.lattice(-4, 4) : io.real(-5.0, 5.0));
    bool p1 = io.boolean();
    FIIntervalInterval<double> query;
    auto r = query(a0, p0, a1, p1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outInt(r.type);
    io.outReal(r.overlap[0]);
    io.outReal(r.overlap[1]);
}

// The dynamic FIQuery's "interval0 initially to the left of interval1" branch
// computes the contact point from the moved *left* endpoint of interval0, which
// the port fixes (UPSTREAM-FINDINGS, issue #62). The generator here excludes
// that branch by swapping the intervals when interval0 lies strictly to the
// left of interval1; the defect itself is covered by the deviation case below.
ORACLE_CASE("IntrIntervals.findDynamic.soundBranches")
{
    int mode = io.index() % 2;
    double maxTime = io.real(0.5, 4.0);
    double a0{}, b0{}, a1{}, b1{};
    if (mode == 0)
    {
        a0 = static_cast<double>(io.rawInteger(-4, 4));
        b0 = static_cast<double>(io.rawInteger(-4, 4));
        a1 = static_cast<double>(io.rawInteger(-4, 4));
        b1 = static_cast<double>(io.rawInteger(-4, 4));
    }
    else
    {
        a0 = io.raw(-5.0, 5.0);
        b0 = io.raw(-5.0, 5.0);
        a1 = io.raw(-5.0, 5.0);
        b1 = io.raw(-5.0, 5.0);
    }
    std::array<double, 2> interval0{ std::min(a0, b0), std::max(a0, b0) };
    std::array<double, 2> interval1{ std::min(a1, b1), std::max(a1, b1) };
    if (interval0[1] < interval1[0])
    {
        std::swap(interval0, interval1);
    }
    io.given(interval0[0]);
    io.given(interval0[1]);
    double speed0 = (mode == 0 ? io.lattice(-3, 3) : io.real(-4.0, 4.0));
    io.given(interval1[0]);
    io.given(interval1[1]);
    double speed1 = (mode == 0 ? io.lattice(-3, 3) : io.real(-4.0, 4.0));
    FIIntervalInterval<double> query;
    auto r = query(maxTime, interval0, speed0, interval1, speed1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outInt(r.type);
    io.outReal(r.firstTime);
    io.outReal(r.lastTime);
    io.outReal(r.overlap[0]);
    io.outReal(r.overlap[1]);
}

// Deliberate deviation: interval0 is strictly to the left of interval1 and the
// intervals approach each other, which is the branch where upstream reports
// interval0[0] + firstTime * speed0 instead of the contact point
// interval0[1] + firstTime * speed0. interval0 is non-degenerate so the two
// differ on every record.
ORACLE_CASE("IntrIntervals.findDynamic.leftApproachDeviation")
{
    double maxTime = io.real(0.5, 4.0);
    double left = io.real(-6.0, -2.0);
    double width = io.real(0.5, 2.0);
    std::array<double, 2> interval0{ left, left + width };
    io.given(interval0[0]);
    io.given(interval0[1]);
    double speed1 = io.real(-2.0, 2.0);
    double approach = io.real(0.25, 3.0);
    double speed0 = io.given(speed1 + approach);
    double gap = io.real(0.25, 4.0);
    double right = interval0[1] + gap;
    double width1 = io.real(0.5, 2.0);
    std::array<double, 2> interval1{ right, right + width1 };
    io.given(interval1[0]);
    io.given(interval1[1]);
    io.given(speed1);
    FIIntervalInterval<double> query;
    auto r = query(maxTime, interval0, speed0, interval1, speed1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outInt(r.type);
    io.outReal(r.firstTime);
    io.outReal(r.lastTime);
    io.outReal(r.overlap[0]);
    io.outReal(r.overlap[1]);
}

// ======================= aligned box / aligned box ========================

ORACLE_CASE("IntrAlignedBox2AlignedBox2.test")
{
    int mode = io.index() % 2;
    auto box0 = ABox<2>(io, mode, 3, 4.0);
    auto box1 = ABox<2>(io, mode, 3, 4.0);
    TIQuery<double, AlignedBox2<double>, AlignedBox2<double>> query;
    auto r = query(box0, box1);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrAlignedBox2AlignedBox2.find")
{
    int mode = io.index() % 2;
    auto box0 = ABox<2>(io, mode, 3, 4.0);
    auto box1 = ABox<2>(io, mode, 3, 4.0);
    FIQuery<double, AlignedBox2<double>, AlignedBox2<double>> query;
    auto r = query(box0, box1);
    io.outBool(r.intersect);
    // The result box is default constructed (min = -1, max = +1) when the
    // boxes are separated, so it is well defined on both paths.
    io.outVec(r.box.min);
    io.outVec(r.box.max);
}

ORACLE_CASE("IntrAlignedBox3AlignedBox3.test")
{
    int mode = io.index() % 2;
    auto box0 = ABox<3>(io, mode, 3, 4.0);
    auto box1 = ABox<3>(io, mode, 3, 4.0);
    TIQuery<double, AlignedBox3<double>, AlignedBox3<double>> query;
    auto r = query(box0, box1);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrAlignedBox3AlignedBox3.find")
{
    int mode = io.index() % 2;
    auto box0 = ABox<3>(io, mode, 3, 4.0);
    auto box1 = ABox<3>(io, mode, 3, 4.0);
    FIQuery<double, AlignedBox3<double>, AlignedBox3<double>> query;
    auto r = query(box0, box1);
    io.outBool(r.intersect);
    io.outVec(r.box.min);
    io.outVec(r.box.max);
}

// ======================= aligned box / oriented box =======================

ORACLE_CASE("IntrAlignedBox2OrientedBox2.test")
{
    int mode = io.index() % 2;
    auto box0 = ABox<2>(io, mode, 3, 4.0);
    auto box1 = OBox<2>(io, mode, 3, 3.0);
    TIQuery<double, AlignedBox2<double>, OrientedBox2<double>> query;
    auto r = query(box0, box1);
    io.outBool(r.intersect);
    io.outInt(r.separating);
}

ORACLE_CASE("IntrAlignedBox3OrientedBox3.test")
{
    int mode = io.index() % 2;
    auto box0 = ABox<3>(io, mode, 3, 4.0);
    auto box1 = OBox<3>(io, mode, 3, 3.0);
    // A negative epsilon is clamped to zero by the query; a large epsilon
    // reaches the existsParallelPair short circuit for non-parallel frames.
    double epsilon = io.real(-0.25, 0.6);
    TIQuery<double, AlignedBox3<double>, OrientedBox3<double>> query;
    auto r = query(box0, box1, epsilon);
    io.outBool(r.intersect);
    io.outInt(r.separating[0]);
    io.outInt(r.separating[1]);
}

// ======================= oriented box / oriented box ======================

ORACLE_CASE("IntrOrientedBox2OrientedBox2.test")
{
    int mode = io.index() % 2;
    auto box0 = OBox<2>(io, mode, 3, 3.0);
    auto box1 = OBox<2>(io, mode, 3, 3.0);
    TIQuery<double, OrientedBox2<double>, OrientedBox2<double>> query;
    auto r = query(box0, box1);
    io.outBool(r.intersect);
    io.outInt(r.separating);
}

ORACLE_CASE("IntrOrientedBox2OrientedBox2.find")
{
    int mode = io.index() % 2;
    auto box0 = OBox<2>(io, mode, 2, 2.5);
    auto box1 = OBox<2>(io, mode, 2, 2.5);
    FIQuery<double, OrientedBox2<double>, OrientedBox2<double>> query;
    auto r = query(box0, box1);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.polygon.size()));
    for (auto const& p : r.polygon) { io.outVec(p); }
}

ORACLE_CASE("IntrOrientedBox3OrientedBox3.test")
{
    int mode = io.index() % 2;
    auto box0 = OBox<3>(io, mode, 3, 3.0);
    auto box1 = OBox<3>(io, mode, 3, 3.0);
    double epsilon = io.real(-0.25, 0.6);
    TIQuery<double, OrientedBox3<double>, OrientedBox3<double>> query;
    auto r = query(box0, box1, epsilon);
    io.outBool(r.intersect);
    io.outInt(r.separating[0]);
    io.outInt(r.separating[1]);
}

// ============================ circle / circle =============================

namespace
{
    // Two circles; mode 2 forces the tangency and coincidence branches.
    void Circles(oracle::Ctx& io, int mode, Circle2<double>& circle0,
        Circle2<double>& circle1)
    {
        if (mode == 0)
        {
            circle0.center = io.latticeVec<2>(-3, 3);
            circle0.radius = io.lattice(1, 4);
            circle1.center = io.latticeVec<2>(-3, 3);
            circle1.radius = io.lattice(1, 4);
        }
        else if (mode == 1)
        {
            circle0.center = io.vec<2>(-3.0, 3.0);
            circle0.radius = io.real(0.25, 3.0);
            circle1.center = io.vec<2>(-3.0, 3.0);
            circle1.radius = io.real(0.25, 3.0);
        }
        else
        {
            // Externally tangent, internally tangent or identical circles,
            // built on a lattice so the tangency tests are exact.
            Vector2<double> c0{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            double r0 = static_cast<double>(io.rawInteger(1, 4));
            int kind = io.rawInteger(0, 2);
            Vector2<double> dir{ 1.0, 0.0 };
            if (io.rawInteger(0, 1) == 1) { dir = Vector2<double>{ 0.0, 1.0 }; }
            double r1 = static_cast<double>(io.rawInteger(1, 4));
            Vector2<double> c1{};
            if (kind == 0)
            {
                c1 = c0 + (r0 + r1) * dir;      // externally tangent
            }
            else if (kind == 1)
            {
                c1 = c0 + std::fabs(r0 - r1) * dir;  // internally tangent
            }
            else
            {
                c1 = c0;
                r1 = r0;                         // identical
            }
            circle0.center = io.givenVec(c0);
            circle0.radius = io.given(r0);
            circle1.center = io.givenVec(c1);
            circle1.radius = io.given(r1);
        }
    }
}

ORACLE_CASE("IntrCircle2Circle2.test")
{
    int mode = io.index() % 3;
    Circle2<double> circle0{}, circle1{};
    Circles(io, mode, circle0, circle1);
    TIQuery<double, Circle2<double>, Circle2<double>> query;
    auto r = query(circle0, circle1);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrCircle2Circle2.find")
{
    int mode = io.index() % 3;
    Circle2<double> circle0{}, circle1{};
    Circles(io, mode, circle0, circle1);
    FIQuery<double, Circle2<double>, Circle2<double>> query;
    auto r = query(circle0, circle1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    if (r.numIntersections == std::numeric_limits<int32_t>::max())
    {
        io.outVec(r.circle.center);
        io.outReal(r.circle.radius);
    }
    else if (r.numIntersections > 0)
    {
        io.outVec(r.point[0]);
        if (r.numIntersections == 2) { io.outVec(r.point[1]); }
    }
}

// ============================ halfspace queries ===========================

namespace
{
    // A convex polygon listed counterclockwise. Mode 0 is an axis-aligned
    // lattice rectangle (exact), the other modes inscribe the polygon in a
    // circle with increasing angles.
    std::vector<Vector2<double>> ConvexPolygon(oracle::Ctx& io, int mode)
    {
        std::vector<Vector2<double>> polygon;
        if (mode == 0)
        {
            io.integer(4, 4);
            double x0 = static_cast<double>(io.rawInteger(-3, 0));
            double x1 = x0 + static_cast<double>(io.rawInteger(1, 3));
            double y0 = static_cast<double>(io.rawInteger(-3, 0));
            double y1 = y0 + static_cast<double>(io.rawInteger(1, 3));
            polygon.push_back(io.givenVec(Vector2<double>{ x0, y0 }));
            polygon.push_back(io.givenVec(Vector2<double>{ x1, y0 }));
            polygon.push_back(io.givenVec(Vector2<double>{ x1, y1 }));
            polygon.push_back(io.givenVec(Vector2<double>{ x0, y1 }));
        }
        else
        {
            int n = io.integer(3, 6);
            double radius = io.raw(1.0, 3.0);
            Vector2<double> center{ io.raw(-2.0, 2.0), io.raw(-2.0, 2.0) };
            std::vector<double> angle(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i)
            {
                angle[static_cast<size_t>(i)] = io.raw(0.0, 6.283185307179586);
            }
            std::sort(angle.begin(), angle.end());
            for (int i = 0; i < n; ++i)
            {
                double a = angle[static_cast<size_t>(i)];
                Vector2<double> v{ center[0] + radius * std::cos(a),
                    center[1] + radius * std::sin(a) };
                polygon.push_back(io.givenVec(v));
            }
        }
        return polygon;
    }
}

ORACLE_CASE("IntrHalfspace2Polygon2.find")
{
    int mode = io.index() % 3;
    auto polygon = ConvexPolygon(io, mode);
    // Mode 2 puts the halfspace boundary exactly through one polygon vertex,
    // which reaches the distance == 0 branches of the clipper.
    Vector2<double> normal{};
    if (mode == 0)
    {
        normal = io.latticeDir<2>(-2, 2);
    }
    else
    {
        normal = io.unit<2>();
    }
    double constant{};
    if (mode == 1)
    {
        constant = io.real(-3.0, 3.0);
    }
    else
    {
        size_t k = static_cast<size_t>(io.rawInteger(0,
            static_cast<int>(polygon.size()) - 1));
        constant = io.given(Dot(normal, polygon[k]));
    }
    Halfspace<2, double> halfspace{};
    halfspace.normal = normal;
    halfspace.constant = constant;
    FIQuery<double, Halfspace<2, double>, std::vector<Vector2<double>>> query;
    auto r = query(halfspace, polygon);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.polygon.size()));
    for (auto const& p : r.polygon) { io.outVec(p); }
}

ORACLE_CASE("IntrHalfspace3OrientedBox3.test")
{
    int mode = io.index() % 2;
    auto box = OBox<3>(io, mode, 3, 3.0);
    auto halfspace = Space<3>(io, mode, box.center);
    TIQuery<double, Halfspace3<double>, OrientedBox3<double>> query;
    auto r = query(halfspace, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrHalfspace3Sphere3.test")
{
    int mode = io.index() % 2;
    Sphere3<double> sphere{};
    sphere.center = Pt<3>(io, mode, 3, 3.0);
    sphere.radius = (mode == 0 ? io.lattice(1, 3) : io.real(0.25, 3.0));
    auto halfspace = Space<3>(io, mode, sphere.center);
    TIQuery<double, Halfspace3<double>, Sphere3<double>> query;
    auto r = query(halfspace, sphere);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrHalfspace3Segment3.test")
{
    int mode = io.index() % 2;
    auto p0 = Pt<3>(io, mode, 3, 3.0);
    auto p1 = Pt<3>(io, mode, 3, 3.0);
    auto halfspace = Space<3>(io, mode, p0);
    Segment3<double> segment(p0, p1);
    TIQuery<double, Halfspace3<double>, Segment3<double>> query;
    auto r = query(halfspace, segment);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrHalfspace3Segment3.find")
{
    int mode = io.index() % 2;
    auto p0 = Pt<3>(io, mode, 3, 3.0);
    auto p1 = Pt<3>(io, mode, 3, 3.0);
    auto halfspace = Space<3>(io, mode, p0);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Halfspace3<double>, Segment3<double>> query;
    auto r = query(halfspace, segment);
    io.outBool(r.intersect);
    io.outInt(r.numPoints);
    for (int32_t i = 0; i < r.numPoints; ++i) { io.outVec(r.point[i]); }
}

ORACLE_CASE("IntrHalfspace3Triangle3.test")
{
    int mode = io.index() % 2;
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    auto halfspace = Space<3>(io, mode, triangle.v[0]);
    TIQuery<double, Halfspace3<double>, Triangle3<double>> query;
    auto r = query(halfspace, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrHalfspace3Triangle3.find")
{
    int mode = io.index() % 2;
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    auto halfspace = Space<3>(io, mode, triangle.v[0]);
    FIQuery<double, Halfspace3<double>, Triangle3<double>> query;
    auto r = query(halfspace, triangle);
    io.outBool(r.intersect);
    io.outInt(r.numPoints);
    for (int32_t i = 0; i < r.numPoints; ++i) { io.outVec(r.point[i]); }
}

// ============================== line queries ==============================

namespace
{
    // A line whose direction is not required to be unit length. Mode 0 uses
    // a small lattice for origin and direction so parallel and collinear
    // configurations occur exactly.
    template <int N>
    Line<N, double> Ln(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> origin = Pt<N>(io, mode, lat, range);
        Vector<N, double> direction = (mode == 0 ? io.latticeDir<N>(-2, 2) : io.unit<N>());
        return Line<N, double>(origin, direction);
    }

    // Exposes the protected DoQuery helpers, which the port exports as public
    // free functions.
    template <typename Q>
    struct Expose : public Q
    {
        using Q::DoQuery;
    };
}

ORACLE_CASE("IntrLine2AlignedBox2.test")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 4.0);
    auto box = ABox<2>(io, mode, 3, 4.0);
    TIQuery<double, Line2<double>, AlignedBox2<double>> query;
    auto r = query(line, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine2AlignedBox2.find")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 4.0);
    auto box = ABox<2>(io, mode, 3, 4.0);
    FIQuery<double, Line2<double>, AlignedBox2<double>> query;
    auto r = query(line, box);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

ORACLE_CASE("IntrLine2AlignedBox2.doQuery.ti")
{
    int mode = io.index() % 2;
    auto lineOrigin = Pt<2>(io, mode, 3, 4.0);
    auto lineDirection = (mode == 0 ? io.latticeDir<2>(-2, 2) : io.unit<2>());
    auto boxExtent = Extent<2>(io, mode, 3, 3.0);
    Expose<TIQuery<double, Line2<double>, AlignedBox2<double>>> query;
    TIQuery<double, Line2<double>, AlignedBox2<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, boxExtent, r);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine2AlignedBox2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto lineOrigin = Pt<2>(io, mode, 3, 4.0);
    auto lineDirection = (mode == 0 ? io.latticeDir<2>(-2, 2) : io.unit<2>());
    auto boxExtent = Extent<2>(io, mode, 3, 3.0);
    Expose<FIQuery<double, Line2<double>, AlignedBox2<double>>> query;
    FIQuery<double, Line2<double>, AlignedBox2<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, boxExtent, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

ORACLE_CASE("IntrLine3AlignedBox3.test")
{
    int mode = io.index() % 2;
    auto line = Ln<3>(io, mode, 3, 4.0);
    auto box = ABox<3>(io, mode, 3, 4.0);
    TIQuery<double, Line3<double>, AlignedBox3<double>> query;
    auto r = query(line, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3AlignedBox3.find")
{
    int mode = io.index() % 2;
    auto line = Ln<3>(io, mode, 3, 4.0);
    auto box = ABox<3>(io, mode, 3, 4.0);
    FIQuery<double, Line3<double>, AlignedBox3<double>> query;
    auto r = query(line, box);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    if (r.intersect)
    {
        // The query fills both parameters and both points whenever it
        // reports an intersection (parameter[1] = parameter[0] when the
        // intersection is a single point).
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrLine3AlignedBox3.doQuery.ti")
{
    int mode = io.index() % 2;
    auto lineOrigin = Pt<3>(io, mode, 3, 4.0);
    auto lineDirection = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    auto boxExtent = Extent<3>(io, mode, 3, 3.0);
    Expose<TIQuery<double, Line3<double>, AlignedBox3<double>>> query;
    TIQuery<double, Line3<double>, AlignedBox3<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, boxExtent, r);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3AlignedBox3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto lineOrigin = Pt<3>(io, mode, 3, 4.0);
    auto lineDirection = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    auto boxExtent = Extent<3>(io, mode, 3, 3.0);
    Expose<FIQuery<double, Line3<double>, AlignedBox3<double>>> query;
    FIQuery<double, Line3<double>, AlignedBox3<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, boxExtent, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

namespace
{
    // Two 2D lines. Modes 2 and 3 force the parallel branches: the same line
    // and parallel-but-distinct lines.
    void Lines2(oracle::Ctx& io, int mode, Line2<double>& line0, Line2<double>& line1)
    {
        Vector2<double> origin0{}, direction0{}, origin1{}, direction1{};
        if (mode == 0)
        {
            for (int i = 0; i < 2; ++i)
            {
                origin0[i] = static_cast<double>(io.rawInteger(-3, 3));
            }
            do
            {
                for (int i = 0; i < 2; ++i)
                {
                    direction0[i] = static_cast<double>(io.rawInteger(-2, 2));
                }
            } while (direction0[0] == 0.0 && direction0[1] == 0.0);
            for (int i = 0; i < 2; ++i)
            {
                origin1[i] = static_cast<double>(io.rawInteger(-3, 3));
            }
            do
            {
                for (int i = 0; i < 2; ++i)
                {
                    direction1[i] = static_cast<double>(io.rawInteger(-2, 2));
                }
            } while (direction1[0] == 0.0 && direction1[1] == 0.0);
        }
        else if (mode == 1)
        {
            for (int i = 0; i < 2; ++i) { origin0[i] = io.raw(-3.0, 3.0); }
            double len = 0.0;
            do
            {
                for (int i = 0; i < 2; ++i) { direction0[i] = io.raw(-1.0, 1.0); }
                len = Length(direction0);
            } while (len < 0.1 || len > 1.0);
            Normalize(direction0);
            for (int i = 0; i < 2; ++i) { origin1[i] = io.raw(-3.0, 3.0); }
            do
            {
                for (int i = 0; i < 2; ++i) { direction1[i] = io.raw(-1.0, 1.0); }
                len = Length(direction1);
            } while (len < 0.1 || len > 1.0);
            Normalize(direction1);
        }
        else
        {
            Vector2<double> p0{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            Vector2<double> d0{ static_cast<double>(io.rawInteger(-2, 2)),
                static_cast<double>(io.rawInteger(-2, 2)) };
            if (d0[0] == 0.0 && d0[1] == 0.0) { d0[0] = 1.0; }
            double scale = static_cast<double>(io.rawInteger(1, 3));
            if (io.rawInteger(0, 1) == 1) { scale = -scale; }
            double t = static_cast<double>(io.rawInteger(-3, 3));
            Vector2<double> p1 = p0 + t * d0;
            if (mode == 3)
            {
                // Offset perpendicular to the common direction, so the lines
                // are parallel but distinct.
                double offset = static_cast<double>(io.rawInteger(1, 3));
                p1 += offset * Vector2<double>{ -d0[1], d0[0] };
            }
            origin0 = p0;
            direction0 = d0;
            origin1 = p1;
            direction1 = scale * d0;
        }
        line0.origin = io.givenVec(origin0);
        line0.direction = io.givenVec(direction0);
        line1.origin = io.givenVec(origin1);
        line1.direction = io.givenVec(direction1);
    }
}

ORACLE_CASE("IntrLine2Line2.test")
{
    int mode = io.index() % 4;
    Line2<double> line0{}, line1{};
    Lines2(io, mode, line0, line1);
    TIQuery<double, Line2<double>, Line2<double>> query;
    auto r = query(line0, line1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

ORACLE_CASE("IntrLine2Line2.find")
{
    int mode = io.index() % 4;
    Line2<double> line0{}, line1{};
    Lines2(io, mode, line0, line1);
    FIQuery<double, Line2<double>, Line2<double>> query;
    auto r = query(line0, line1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.line0Parameter[0]);
    io.outReal(r.line0Parameter[1]);
    io.outReal(r.line1Parameter[0]);
    io.outReal(r.line1Parameter[1]);
    io.outVec(r.point);
}

ORACLE_CASE("IntrLine2Triangle2.test")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 3.0);
    auto triangle = Tri<2>(io, mode, 3, 3.0);
    TIQuery<double, Line2<double>, Triangle2<double>> query;
    auto r = query(line, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine2Triangle2.find")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 3.0);
    auto triangle = Tri<2>(io, mode, 3, 3.0);
    FIQuery<double, Line2<double>, Triangle2<double>> query;
    auto r = query(line, triangle);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    if (r.intersect)
    {
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrLine2Triangle2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto origin = Pt<2>(io, mode, 3, 3.0);
    auto direction = (mode == 0 ? io.latticeDir<2>(-2, 2) : io.unit<2>());
    auto triangle = Tri<2>(io, mode, 3, 3.0);
    Expose<FIQuery<double, Line2<double>, Triangle2<double>>> query;
    FIQuery<double, Line2<double>, Triangle2<double>>::Result r{};
    query.DoQuery(origin, direction, triangle, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

namespace
{
    // A rectangle with unit-length orthogonal axes.
    Rectangle3<double> Rect3(oracle::Ctx& io, int mode, int lat, double range)
    {
        Rectangle3<double> rectangle{};
        rectangle.center = Pt<3>(io, mode, lat, range);
        auto axis = Frame3(io, mode);
        rectangle.axis[0] = axis[0];
        rectangle.axis[1] = axis[1];
        rectangle.extent = Extent<2>(io, mode, lat, range);
        return rectangle;
    }
}

ORACLE_CASE("IntrLine3Rectangle3.test")
{
    int mode = io.index() % 2;
    auto line = Ln<3>(io, mode, 3, 3.0);
    auto rectangle = Rect3(io, mode, 3, 3.0);
    TIQuery<double, Line3<double>, Rectangle3<double>> query;
    auto r = query(line, rectangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3Rectangle3.find")
{
    int mode = io.index() % 2;
    auto line = Ln<3>(io, mode, 3, 3.0);
    auto rectangle = Rect3(io, mode, 3, 3.0);
    FIQuery<double, Line3<double>, Rectangle3<double>> query;
    auto r = query(line, rectangle);
    io.outBool(r.intersect);
    if (r.intersect)
    {
        io.outReal(r.parameter);
        io.outReal(r.rectCoord[0]);
        io.outReal(r.rectCoord[1]);
        io.outReal(r.rectCoord[2]);
        io.outVec(r.point);
    }
}

ORACLE_CASE("IntrLine3Sphere3.test")
{
    int mode = io.index() % 2;
    auto origin = Pt<3>(io, mode, 3, 3.0);
    auto direction = io.unit<3>();
    Sphere3<double> sphere{};
    sphere.center = Pt<3>(io, mode, 3, 3.0);
    sphere.radius = (mode == 0 ? io.lattice(1, 3) : io.real(0.25, 3.0));
    Line3<double> line(origin, direction);
    TIQuery<double, Line3<double>, Sphere3<double>> query;
    auto r = query(line, sphere);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3Sphere3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto lineOrigin = Pt<3>(io, mode, 3, 3.0);
    auto lineDirection = io.unit<3>();
    Sphere3<double> sphere{};
    sphere.center = Pt<3>(io, mode, 3, 3.0);
    sphere.radius = (mode == 0 ? io.lattice(1, 3) : io.real(0.25, 3.0));
    Expose<FIQuery<double, Line3<double>, Sphere3<double>>> query;
    FIQuery<double, Line3<double>, Sphere3<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, sphere, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

ORACLE_CASE("IntrLine3Sphere3.find")
{
    int mode = io.index() % 2;
    auto origin = Pt<3>(io, mode, 3, 3.0);
    auto direction = io.unit<3>();
    Sphere3<double> sphere{};
    sphere.center = Pt<3>(io, mode, 3, 3.0);
    sphere.radius = (mode == 0 ? io.lattice(1, 3) : io.real(0.25, 3.0));
    Line3<double> line(origin, direction);
    FIQuery<double, Line3<double>, Sphere3<double>> query;
    auto r = query(line, sphere);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    if (r.intersect)
    {
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrLine3Triangle3.test")
{
    int mode = io.index() % 2;
    auto line = Ln<3>(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    TIQuery<double, Line3<double>, Triangle3<double>> query;
    auto r = query(line, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3Triangle3.find")
{
    int mode = io.index() % 2;
    auto line = Ln<3>(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    FIQuery<double, Line3<double>, Triangle3<double>> query;
    auto r = query(line, triangle);
    io.outBool(r.intersect);
    if (r.intersect)
    {
        io.outReal(r.parameter);
        io.outReal(r.triangleBary[0]);
        io.outReal(r.triangleBary[1]);
        io.outReal(r.triangleBary[2]);
        io.outVec(r.point);
    }
}

ORACLE_CASE("IntrRay3Triangle3.test")
{
    int mode = io.index() % 2;
    auto origin = Pt<3>(io, mode, 3, 3.0);
    auto direction = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    Ray3<double> ray(origin, direction);
    TIQuery<double, Ray3<double>, Triangle3<double>> query;
    auto r = query(ray, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3Triangle3.find")
{
    int mode = io.index() % 2;
    auto origin = Pt<3>(io, mode, 3, 3.0);
    auto direction = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, Triangle3<double>> query;
    auto r = query(ray, triangle);
    io.outBool(r.intersect);
    if (r.intersect)
    {
        io.outReal(r.parameter);
        io.outReal(r.triangleBary[0]);
        io.outReal(r.triangleBary[1]);
        io.outReal(r.triangleBary[2]);
        io.outVec(r.point);
    }
}

ORACLE_CASE("IntrSegment3Triangle3.test")
{
    int mode = io.index() % 2;
    auto p0 = Pt<3>(io, mode, 3, 3.0);
    auto p1 = Pt<3>(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    Segment3<double> segment(p0, p1);
    TIQuery<double, Segment3<double>, Triangle3<double>> query;
    auto r = query(segment, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3Triangle3.find")
{
    int mode = io.index() % 2;
    auto p0 = Pt<3>(io, mode, 3, 3.0);
    auto p1 = Pt<3>(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Triangle3<double>> query;
    auto r = query(segment, triangle);
    io.outBool(r.intersect);
    if (r.intersect)
    {
        io.outReal(r.parameter);
        io.outReal(r.triangleBary[0]);
        io.outReal(r.triangleBary[1]);
        io.outReal(r.triangleBary[2]);
        io.outVec(r.point);
    }
}

// ===================== constructed hits and tangencies =====================
//
// Uniform inputs hit a triangle, a rectangle or a sphere only a few percent of
// the time and never land exactly on a boundary. These cases aim the line, ray
// or segment at a target point computed from the primitive itself, on a small
// lattice so that every product and sum stays exactly representable. The
// boundary comparisons of the upstream code (b1 == 0, b2 == 0, b1 + b2 == 1,
// |s| == extent, discriminant == 0, ray parameter == 0) are then reached with
// exact equality, which is where a reassociated port formula would show up.

namespace
{
    // A point of the plane of 'triangle' with barycentric coordinates
    // b1 = k1/4, b2 = k2/4, k1, k2 in [0, 5]. k = 0 lands on an edge through
    // v0, k1 + k2 == 4 on the opposite edge, k1 + k2 > 4 outside. The
    // triangle is integer valued, so the target is exact.
    Vector3<double> TriTarget(oracle::Ctx& io, Triangle3<double> const& triangle)
    {
        double k1 = static_cast<double>(io.rawInteger(0, 5));
        double k2 = static_cast<double>(io.rawInteger(0, 5));
        Vector3<double> edge1 = triangle.v[1] - triangle.v[0];
        Vector3<double> edge2 = triangle.v[2] - triangle.v[0];
        return triangle.v[0] + (k1 * edge1 + k2 * edge2) / 4.0;
    }

    // A lattice triangle, a lattice direction and a line origin placed so that
    // the line meets the triangle's plane at the target point at parameter
    // 't' (returned through 'parameter'). Records: 9 triangle, 3 direction,
    // 3 origin.
    void LineAtTriangle(oracle::Ctx& io, Triangle3<double>& triangle,
        Vector3<double>& origin, Vector3<double>& direction, double& parameter)
    {
        triangle = Tri<3>(io, 0, 3, 3.0);
        Vector3<double> target = TriTarget(io, triangle);
        direction = io.latticeDir<3>(-2, 2);
        parameter = static_cast<double>(io.rawInteger(-3, 3));
        origin = io.givenVec(target - parameter * direction);
    }
}

ORACLE_CASE("IntrLine3Triangle3.test.throughPoint")
{
    Triangle3<double> triangle{};
    Vector3<double> origin{}, direction{};
    double parameter = 0.0;
    LineAtTriangle(io, triangle, origin, direction, parameter);
    Line3<double> line(origin, direction);
    TIQuery<double, Line3<double>, Triangle3<double>> query;
    auto r = query(line, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3Triangle3.find.throughPoint")
{
    Triangle3<double> triangle{};
    Vector3<double> origin{}, direction{};
    double parameter = 0.0;
    LineAtTriangle(io, triangle, origin, direction, parameter);
    Line3<double> line(origin, direction);
    FIQuery<double, Line3<double>, Triangle3<double>> query;
    auto r = query(line, triangle);
    io.outBool(r.intersect);
    if (r.intersect)
    {
        io.outReal(r.parameter);
        io.outReal(r.triangleBary[0]);
        io.outReal(r.triangleBary[1]);
        io.outReal(r.triangleBary[2]);
        io.outVec(r.point);
    }
}

namespace
{
    // As LineAtTriangle, but the parameter of the target point ranges over
    // negative, zero (the ray origin lies on the triangle, the QdN >= 0
    // boundary) and positive values.
    void RayAtTriangle(oracle::Ctx& io, Triangle3<double>& triangle,
        Vector3<double>& origin, Vector3<double>& direction)
    {
        triangle = Tri<3>(io, 0, 3, 3.0);
        Vector3<double> target = TriTarget(io, triangle);
        direction = io.latticeDir<3>(-2, 2);
        double parameter = static_cast<double>(io.rawInteger(-2, 3));
        origin = io.givenVec(target - parameter * direction);
    }
}

ORACLE_CASE("IntrRay3Triangle3.test.throughPoint")
{
    Triangle3<double> triangle{};
    Vector3<double> origin{}, direction{};
    RayAtTriangle(io, triangle, origin, direction);
    Ray3<double> ray(origin, direction);
    TIQuery<double, Ray3<double>, Triangle3<double>> query;
    auto r = query(ray, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3Triangle3.find.throughPoint")
{
    Triangle3<double> triangle{};
    Vector3<double> origin{}, direction{};
    RayAtTriangle(io, triangle, origin, direction);
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, Triangle3<double>> query;
    auto r = query(ray, triangle);
    io.outBool(r.intersect);
    if (r.intersect)
    {
        io.outReal(r.parameter);
        io.outReal(r.triangleBary[0]);
        io.outReal(r.triangleBary[1]);
        io.outReal(r.triangleBary[2]);
        io.outVec(r.point);
    }
}

namespace
{
    // A segment whose supporting line meets the triangle at the target point.
    // t0 or t1 equal to zero puts an endpoint exactly on the triangle, which
    // is the |t| <= extent boundary; t0 == t1 == 0 makes the segment
    // degenerate, which upstream treats as parallel (DdN == 0).
    void SegmentAtTriangle(oracle::Ctx& io, Triangle3<double>& triangle,
        Vector3<double>& p0, Vector3<double>& p1)
    {
        triangle = Tri<3>(io, 0, 3, 3.0);
        Vector3<double> target = TriTarget(io, triangle);
        Vector3<double> direction = io.latticeDir<3>(-2, 2);
        double t0 = static_cast<double>(io.rawInteger(0, 3));
        p0 = io.givenVec(target - t0 * direction);
        double t1 = static_cast<double>(io.rawInteger(0, 3));
        p1 = io.givenVec(target + t1 * direction);
    }
}

ORACLE_CASE("IntrSegment3Triangle3.test.throughPoint")
{
    Triangle3<double> triangle{};
    Vector3<double> p0{}, p1{};
    SegmentAtTriangle(io, triangle, p0, p1);
    Segment3<double> segment(p0, p1);
    TIQuery<double, Segment3<double>, Triangle3<double>> query;
    auto r = query(segment, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3Triangle3.find.throughPoint")
{
    Triangle3<double> triangle{};
    Vector3<double> p0{}, p1{};
    SegmentAtTriangle(io, triangle, p0, p1);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Triangle3<double>> query;
    auto r = query(segment, triangle);
    io.outBool(r.intersect);
    if (r.intersect)
    {
        io.outReal(r.parameter);
        io.outReal(r.triangleBary[0]);
        io.outReal(r.triangleBary[1]);
        io.outReal(r.triangleBary[2]);
        io.outVec(r.point);
    }
}

namespace
{
    // A rectangle and a line through the point
    //   center + (a/2) * extent[0] * axis[0] + (b/2) * extent[1] * axis[1],
    // a, b in [-3, 3]. |a| == 2 is the |W1dDxQ| == extent[0] * |DdN| boundary
    // (a corner when |a| == |b| == 2), |a| == 3 is outside. Mode 0 uses the
    // exact signed-permutation frame and integer extents, so the boundary
    // comparisons are exact.
    void LineAtRectangle(oracle::Ctx& io, int mode, Rectangle3<double>& rectangle,
        Vector3<double>& origin, Vector3<double>& direction)
    {
        rectangle = Rect3(io, mode, 3, 3.0);
        double a = static_cast<double>(io.rawInteger(-3, 3)) / 2.0;
        double b = static_cast<double>(io.rawInteger(-3, 3)) / 2.0;
        Vector3<double> target = rectangle.center
            + (a * rectangle.extent[0]) * rectangle.axis[0]
            + (b * rectangle.extent[1]) * rectangle.axis[1];
        direction = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
        double parameter = static_cast<double>(io.rawInteger(-3, 3));
        origin = io.givenVec(target - parameter * direction);
    }
}

ORACLE_CASE("IntrLine3Rectangle3.test.throughPoint")
{
    int mode = io.index() % 2;
    Rectangle3<double> rectangle{};
    Vector3<double> origin{}, direction{};
    LineAtRectangle(io, mode, rectangle, origin, direction);
    Line3<double> line(origin, direction);
    TIQuery<double, Line3<double>, Rectangle3<double>> query;
    auto r = query(line, rectangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3Rectangle3.find.throughPoint")
{
    int mode = io.index() % 2;
    Rectangle3<double> rectangle{};
    Vector3<double> origin{}, direction{};
    LineAtRectangle(io, mode, rectangle, origin, direction);
    Line3<double> line(origin, direction);
    FIQuery<double, Line3<double>, Rectangle3<double>> query;
    auto r = query(line, rectangle);
    io.outBool(r.intersect);
    if (r.intersect)
    {
        io.outReal(r.parameter);
        io.outReal(r.rectCoord[0]);
        io.outReal(r.rectCoord[1]);
        io.outReal(r.rectCoord[2]);
        io.outVec(r.point);
    }
}

namespace
{
    // Sphere of integer radius r centred on a lattice point, with a line whose
    // direction is a signed coordinate axis e_j and whose origin is
    //   center + (r + d) * e_i + s * e_j,   i != j,  d in {-1, 0, 1}.
    // Then a1 = +-s, a0 = (r + d)^2 + s^2 - r^2 and the discriminant is
    // exactly r^2 - (r + d)^2: negative, zero (a tangent line) or positive
    // with no rounding at all.
    void TangentLine(oracle::Ctx& io, Sphere3<double>& sphere,
        Vector3<double>& origin, Vector3<double>& direction)
    {
        sphere.center = io.latticeVec<3>(-3, 3);
        sphere.radius = io.lattice(1, 4);
        int i = io.rawInteger(0, 2);
        int j = (i + 1 + io.rawInteger(0, 1)) % 3;
        double signJ = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        Vector3<double> d3{};
        d3.MakeZero();
        d3[j] = signJ;
        direction = io.givenVec(d3);
        double signI = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        double offset = static_cast<double>(io.rawInteger(-1, 1));
        double along = static_cast<double>(io.rawInteger(-3, 3));
        Vector3<double> p = sphere.center;
        p[i] += signI * (sphere.radius + offset);
        p[j] += along;
        origin = io.givenVec(p);
    }
}

ORACLE_CASE("IntrLine3Sphere3.test.tangent")
{
    Sphere3<double> sphere{};
    Vector3<double> origin{}, direction{};
    TangentLine(io, sphere, origin, direction);
    Line3<double> line(origin, direction);
    TIQuery<double, Line3<double>, Sphere3<double>> query;
    auto r = query(line, sphere);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3Sphere3.find.tangent")
{
    Sphere3<double> sphere{};
    Vector3<double> origin{}, direction{};
    TangentLine(io, sphere, origin, direction);
    Line3<double> line(origin, direction);
    FIQuery<double, Line3<double>, Sphere3<double>> query;
    auto r = query(line, sphere);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    if (r.intersect)
    {
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}
