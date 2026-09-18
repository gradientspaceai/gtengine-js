// Verify group 19 (distance): differential cases for the Dist* headers that
// group 19 covers. See ORACLE.md.
//
// Generator modes are selected from io.index() so that every case mixes
// uniform random inputs, small-lattice inputs (exact integer arithmetic, which
// produces parallel, touching, tangent and coincident configurations) and
// constructions aimed at a specific branch. Every generator records exactly
// the same number of doubles in every mode, so the TypeScript replay reads the
// inputs without knowing which mode produced them.
#define ORACLE_FAMILY "v19-distance"
#include "Oracle.h"

#include <Mathematics/DistAlignedBoxAlignedBox.h>
#include <Mathematics/DistCircle2Circle2.h>
#include <Mathematics/DistLine2AlignedBox2.h>
#include <Mathematics/DistLine2Circle2.h>
#include <Mathematics/DistLine2Triangle2.h>
#include <Mathematics/DistLine3CanonicalBox3.h>
#include <Mathematics/DistLineLine.h>
#include <Mathematics/DistLineRay.h>
#include <Mathematics/DistLineSegment.h>
#include <Mathematics/DistPoint2Circle2.h>
#include <Mathematics/DistPointCanonicalBox.h>
#include <Mathematics/DistPointLine.h>
#include <Mathematics/DistPointRay.h>
#include <Mathematics/DistPointRectangle.h>
#include <Mathematics/DistPointSegment.h>
#include <Mathematics/DistPointTriangle.h>
#include <Mathematics/DistRayRay.h>
#include <Mathematics/DistRaySegment.h>
#include <Mathematics/DistSegmentSegment.h>

#include <cmath>

using namespace gte;

namespace
{
    // ---- generators ------------------------------------------------------
    //
    // mode 0,1: uniform random
    // mode 2,3: small lattice (integers), so degenerate configurations occur
    // mode 4:   constructed (parallel directions, coincident points, ...)
    int Mode(oracle::Ctx& io)
    {
        return io.index() % 5;
    }

    template <int32_t N>
    Vector<N, double> Point(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 2 || mode == 3)
        {
            return io.latticeVec<N>(-3, 3);
        }
        return io.vec<N>(-5.0, 5.0);
    }

    template <int32_t N>
    Vector<N, double> Direction(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 2 || mode == 3)
        {
            return io.latticeDir<N>(-2, 2);
        }
        return io.unit<N>();
    }

    // A direction that is a scalar multiple of d0 in mode 4, so that the
    // parallel branches (det == 0, or det one ulp above zero) are reached.
    template <int32_t N>
    Vector<N, double> SecondDirection(oracle::Ctx& io, int32_t mode,
        Vector<N, double> const& d0)
    {
        if (mode == 4)
        {
            double s = io.raw(-2.0, 2.0);
            double t = (s >= 0.0 ? s + 0.25 : s - 0.25);
            return io.givenVec(t * d0);
        }
        return Direction<N>(io, mode);
    }

    // A point that coincides with p0 in mode 4.
    template <int32_t N>
    Vector<N, double> SecondPoint(oracle::Ctx& io, int32_t mode,
        Vector<N, double> const& p0)
    {
        if (mode == 4)
        {
            return io.givenVec(p0);
        }
        return Point<N>(io, mode);
    }

    // Unrecorded draws, for building derived inputs.
    template <int32_t N>
    Vector<N, double> RawVec(oracle::Ctx& io, double lo, double hi)
    {
        Vector<N, double> v{};
        for (int32_t i = 0; i < N; ++i) { v[i] = io.raw(lo, hi); }
        return v;
    }

    template <int32_t N>
    Vector<N, double> RawUnit(oracle::Ctx& io)
    {
        Vector<N, double> v{};
        double len = 0.0;
        do
        {
            for (int32_t i = 0; i < N; ++i) { v[i] = io.raw(-1.0, 1.0); }
            len = Length(v);
        } while (len < 0.1 || len > 1.0);
        Normalize(v);
        return v;
    }

    // ---- output helpers --------------------------------------------------
    template <typename Result, int32_t N>
    void EmitTwoClosest(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    template <typename Result, int32_t N>
    void EmitOneParameter(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    template <typename Result, int32_t N>
    void EmitTwoParameters(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }
}

// =====================================================================
// DistPointLine.h
// =====================================================================
namespace
{
    template <int32_t N>
    void PointLine(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto point = Point<N>(io, mode);
        auto origin = Point<N>(io, mode);
        auto direction = Direction<N>(io, mode);
        Line<N, double> line(origin, direction);
        DCPQuery<double, Vector<N, double>, Line<N, double>> query{};
        auto r = query(point, line);
        EmitOneParameter<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistPointLine.compute.2d") { PointLine<2>(io); }
ORACLE_CASE("DistPointLine.compute.3d") { PointLine<3>(io); }

// =====================================================================
// DistPointRay.h
// =====================================================================
namespace
{
    template <int32_t N>
    void PointRay(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto point = Point<N>(io, mode);
        auto origin = Point<N>(io, mode);
        auto direction = Direction<N>(io, mode);
        Ray<N, double> ray(origin, direction);
        DCPQuery<double, Vector<N, double>, Ray<N, double>> query{};
        auto r = query(point, ray);
        EmitOneParameter<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistPointRay.compute.2d") { PointRay<2>(io); }
ORACLE_CASE("DistPointRay.compute.3d") { PointRay<3>(io); }

// =====================================================================
// DistPointSegment.h
// =====================================================================
namespace
{
    template <int32_t N>
    void PointSegment(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto point = Point<N>(io, mode);
        auto p0 = Point<N>(io, mode);
        auto p1 = SecondPoint<N>(io, mode, p0);
        Segment<N, double> segment(p0, p1);
        DCPQuery<double, Vector<N, double>, Segment<N, double>> query{};
        auto r = query(point, segment);
        EmitOneParameter<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistPointSegment.compute.2d") { PointSegment<2>(io); }
ORACLE_CASE("DistPointSegment.compute.3d") { PointSegment<3>(io); }

// =====================================================================
// DistPointCanonicalBox.h
// =====================================================================
namespace
{
    // Extents are nonnegative; a zero extent (a degenerate box) is allowed by
    // the header and is generated in the lattice modes.
    template <int32_t N>
    Vector<N, double> Extent(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 2 || mode == 3)
        {
            return io.latticeVec<N>(0, 2);
        }
        return io.vec<N>(0.25, 3.0);
    }

    template <int32_t N>
    void PointCanonicalBox(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto point = Point<N>(io, mode);
        auto extent = Extent<N>(io, mode);
        CanonicalBox<N, double> box(extent);
        DCPQuery<double, Vector<N, double>, CanonicalBox<N, double>> query{};
        auto r = query(point, box);
        EmitTwoClosest<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistPointCanonicalBox.compute.2d") { PointCanonicalBox<2>(io); }
ORACLE_CASE("DistPointCanonicalBox.compute.3d") { PointCanonicalBox<3>(io); }

// =====================================================================
// DistPointRectangle.h
// =====================================================================
namespace
{
    // The rectangle axes must be unit length and perpendicular; they are built
    // from unrecorded draws and recorded with givenVec.
    template <int32_t N>
    void PointRectangle(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto point = Point<N>(io, mode);
        auto center = Point<N>(io, mode);
        Vector<N, double> a0 = RawUnit<N>(io);
        Vector<N, double> a1{};
        double len = 0.0;
        do
        {
            a1 = RawVec<N>(io, -1.0, 1.0);
            a1 = a1 - Dot(a1, a0) * a0;
            len = Length(a1);
        } while (len < 0.25);
        Normalize(a1);
        std::array<Vector<N, double>, 2> axis{};
        axis[0] = io.givenVec(a0);
        axis[1] = io.givenVec(a1);
        auto extent = io.vec<2>(0.25, 3.0);
        Rectangle<N, double> rectangle(center, axis, extent);
        DCPQuery<double, Vector<N, double>, Rectangle<N, double>> query{};
        auto r = query(point, rectangle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.cartesian[0]);
        io.outReal(r.cartesian[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }
}

ORACLE_CASE("DistPointRectangle.compute.2d") { PointRectangle<2>(io); }
ORACLE_CASE("DistPointRectangle.compute.3d") { PointRectangle<3>(io); }

// =====================================================================
// DistPoint2Circle2.h
// =====================================================================
ORACLE_CASE("DistPoint2Circle2.compute")
{
    // Mode 4 puts the point exactly at the circle center, which is the
    // 'equidistant' branch.
    int32_t mode = Mode(io);
    auto center = Point<2>(io, mode);
    auto point = SecondPoint<2>(io, mode, center);
    double radius = (mode == 2 || mode == 3)
        ? static_cast<double>(io.integer(1, 3)) : io.real(0.25, 4.0);
    Circle2<double> circle(center, radius);
    DCPQuery<double, Vector2<double>, Circle2<double>> query{};
    auto r = query(point, circle);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
    io.outBool(r.equidistant);
}

// =====================================================================
// DistLineLine.h
// =====================================================================
namespace
{
    template <int32_t N>
    void LineLine(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto origin0 = Point<N>(io, mode);
        auto direction0 = Direction<N>(io, mode);
        auto origin1 = Point<N>(io, mode);
        auto direction1 = SecondDirection<N>(io, mode, direction0);
        Line<N, double> line0(origin0, direction0);
        Line<N, double> line1(origin1, direction1);
        DCPQuery<double, Line<N, double>, Line<N, double>> query{};
        auto r = query(line0, line1);
        EmitTwoParameters<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistLineLine.compute.2d") { LineLine<2>(io); }
ORACLE_CASE("DistLineLine.compute.3d") { LineLine<3>(io); }

// =====================================================================
// DistLineRay.h
// =====================================================================
namespace
{
    template <int32_t N>
    void LineRay(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto lineOrigin = Point<N>(io, mode);
        auto lineDirection = Direction<N>(io, mode);
        auto rayOrigin = Point<N>(io, mode);
        auto rayDirection = SecondDirection<N>(io, mode, lineDirection);
        Line<N, double> line(lineOrigin, lineDirection);
        Ray<N, double> ray(rayOrigin, rayDirection);
        DCPQuery<double, Line<N, double>, Ray<N, double>> query{};
        auto r = query(line, ray);
        EmitTwoParameters<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistLineRay.compute.2d") { LineRay<2>(io); }
ORACLE_CASE("DistLineRay.compute.3d") { LineRay<3>(io); }

// =====================================================================
// DistLineSegment.h
// =====================================================================
namespace
{
    template <int32_t N>
    void LineSegment(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto lineOrigin = Point<N>(io, mode);
        auto lineDirection = Direction<N>(io, mode);
        auto p0 = Point<N>(io, mode);
        Vector<N, double> p1{};
        if (mode == 4)
        {
            // A segment parallel to the line.
            double s = io.raw(0.25, 2.0);
            p1 = io.givenVec(p0 + s * lineDirection);
        }
        else
        {
            p1 = Point<N>(io, mode);
        }
        Line<N, double> line(lineOrigin, lineDirection);
        Segment<N, double> segment(p0, p1);
        DCPQuery<double, Line<N, double>, Segment<N, double>> query{};
        auto r = query(line, segment);
        EmitTwoParameters<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistLineSegment.compute.2d") { LineSegment<2>(io); }
ORACLE_CASE("DistLineSegment.compute.3d") { LineSegment<3>(io); }

// =====================================================================
// DistRayRay.h
// =====================================================================
namespace
{
    template <int32_t N>
    void RayRay(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto origin0 = Point<N>(io, mode);
        auto direction0 = Direction<N>(io, mode);
        auto origin1 = Point<N>(io, mode);
        auto direction1 = SecondDirection<N>(io, mode, direction0);
        Ray<N, double> ray0(origin0, direction0);
        Ray<N, double> ray1(origin1, direction1);
        DCPQuery<double, Ray<N, double>, Ray<N, double>> query{};
        auto r = query(ray0, ray1);
        EmitTwoParameters<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistRayRay.compute.2d") { RayRay<2>(io); }
ORACLE_CASE("DistRayRay.compute.3d") { RayRay<3>(io); }

// =====================================================================
// DistRaySegment.h
// =====================================================================
namespace
{
    // Unrecorded lattice draws, so a rejection loop can retry without
    // perturbing the recorded input count.
    template <int32_t N>
    Vector<N, double> RawLatticeVec(oracle::Ctx& io, int32_t lo, int32_t hi)
    {
        Vector<N, double> v{};
        for (int32_t i = 0; i < N; ++i)
        {
            v[i] = static_cast<double>(io.rawInteger(lo, hi));
        }
        return v;
    }

    template <int32_t N>
    Vector<N, double> RawLatticeDir(oracle::Ctx& io, int32_t lo, int32_t hi)
    {
        Vector<N, double> v{};
        bool zero = true;
        do
        {
            v = RawLatticeVec<N>(io, lo, hi);
            zero = true;
            for (int32_t i = 0; i < N; ++i) { zero = zero && v[i] == 0.0; }
        } while (zero);
        return v;
    }

    // Upstream never clamps the ray parameter to s0 >= 0 in its parallel
    // branch and in regions 1 and 5 (UPSTREAM-FINDINGS, issue #126); the port
    // clamps. The main case therefore keeps only configurations for which
    // upstream reports s0 >= 0, where the port is identical to upstream. The
    // companion '.deviation' case keeps only the others.
    //
    // The main case draws its candidates through a rejection loop, so the
    // generator modes are built from unrecorded draws and only the accepted
    // configuration is recorded. The lattice modes reach det == 0 exactly (the
    // parallel branch) and the exact region boundaries; mode 4 makes the
    // segment a scalar multiple of the ray direction.
    template <int32_t N>
    void RaySegment(oracle::Ctx& io, bool wantNegative)
    {
        DCPQuery<double, Ray<N, double>, Segment<N, double>> query{};
        int32_t mode = Mode(io);
        Vector<N, double> origin{}, direction{}, p0{}, p1{};
        bool found = false;
        for (int32_t attempt = 0; attempt < 20000 && !found; ++attempt)
        {
            if (wantNegative)
            {
                origin = RawVec<N>(io, -4.0, 4.0);
                direction = RawUnit<N>(io);
                // Bias towards a segment behind the ray origin, where the
                // unconstrained minimizer has s0 < 0.
                double back = io.raw(0.25, 4.0);
                Vector<N, double> base = origin - back * direction;
                p0 = base + RawVec<N>(io, -1.5, 1.5);
                p1 = base + RawVec<N>(io, -1.5, 1.5);
            }
            else if (mode == 2 || mode == 3)
            {
                origin = RawLatticeVec<N>(io, -3, 3);
                direction = RawLatticeDir<N>(io, -2, 2);
                p0 = RawLatticeVec<N>(io, -3, 3);
                p1 = RawLatticeVec<N>(io, -3, 3);
            }
            else if (mode == 4)
            {
                origin = RawVec<N>(io, -4.0, 4.0);
                direction = RawUnit<N>(io);
                p0 = RawVec<N>(io, -4.0, 4.0);
                double s = io.raw(-2.0, 2.0);
                p1 = p0 + (s >= 0.0 ? s + 0.25 : s - 0.25) * direction;
            }
            else
            {
                origin = RawVec<N>(io, -4.0, 4.0);
                direction = RawUnit<N>(io);
                p0 = RawVec<N>(io, -4.0, 4.0);
                p1 = RawVec<N>(io, -4.0, 4.0);
            }
            auto probe = query(Ray<N, double>(origin, direction),
                Segment<N, double>(p0, p1));
            found = (wantNegative ? probe.parameter[0] < 0.0
                : probe.parameter[0] >= 0.0);
        }

        io.givenVec(origin);
        io.givenVec(direction);
        io.givenVec(p0);
        io.givenVec(p1);
        Ray<N, double> ray(origin, direction);
        Segment<N, double> segment(p0, p1);
        auto r = query(ray, segment);
        EmitTwoParameters<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistRaySegment.compute.2d") { RaySegment<2>(io, false); }
ORACLE_CASE("DistRaySegment.compute.3d") { RaySegment<3>(io, false); }
ORACLE_CASE("DistRaySegment.compute.3d.deviation") { RaySegment<3>(io, true); }

// =====================================================================
// DistSegmentSegment.h
// =====================================================================
namespace
{
    template <int32_t N>
    void SegmentSegmentInputs(oracle::Ctx& io, int32_t mode,
        Vector<N, double>& p0, Vector<N, double>& p1,
        Vector<N, double>& q0, Vector<N, double>& q1)
    {
        p0 = Point<N>(io, mode);
        p1 = SecondPoint<N>(io, mode, p0);
        q0 = Point<N>(io, mode);
        if (mode == 4)
        {
            // A segment parallel to segment0 (which is degenerate in mode 4,
            // so this is the degenerate-segment branch of ComputeRobust).
            double s = io.raw(0.25, 2.0);
            Vector<N, double> dir = RawUnit<N>(io);
            q1 = io.givenVec(q0 + s * dir);
        }
        else
        {
            q1 = Point<N>(io, mode);
        }
    }

    template <int32_t N>
    void SegmentSegment(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        Vector<N, double> p0{}, p1{}, q0{}, q1{};
        SegmentSegmentInputs<N>(io, mode, p0, p1, q0, q1);
        Segment<N, double> segment0(p0, p1);
        Segment<N, double> segment1(q0, q1);
        DCPQuery<double, Segment<N, double>, Segment<N, double>> query{};
        auto r = query(segment0, segment1);
        EmitTwoParameters<decltype(r), N>(io, r);
    }

    template <int32_t N>
    void SegmentSegmentRobust(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        Vector<N, double> p0{}, p1{}, q0{}, q1{};
        SegmentSegmentInputs<N>(io, mode, p0, p1, q0, q1);
        Segment<N, double> segment0(p0, p1);
        Segment<N, double> segment1(q0, q1);
        DCPQuery<double, Segment<N, double>, Segment<N, double>> query{};
        auto r = query.ComputeRobust(segment0, segment1);
        EmitTwoParameters<decltype(r), N>(io, r);
    }

    // The upstream GetClampedRoot, replicated so the generator can tell
    // whether ComputeRobust reaches the defect described below.
    double ClampedRoot(double slope, double h0, double h1)
    {
        if (h0 < 0.0)
        {
            if (h1 > 0.0)
            {
                double r = -h0 / slope;
                return (r > 1.0 ? 0.5 : r);
            }
            return 1.0;
        }
        return 0.0;
    }

    // True when ComputeIntersection computes an endpoint t-coordinate f/b
    // outside [0,1] and replaces it by 1/2 (UPSTREAM-FINDINGS, issue #418).
    // The port clamps to the nearest endpoint of [0,1] instead.
    template <int32_t N>
    bool UsesOutOfRangeRatio(Vector<N, double> const& P0,
        Vector<N, double> const& P1, Vector<N, double> const& Q0,
        Vector<N, double> const& Q1)
    {
        Vector<N, double> P1mP0 = P1 - P0;
        Vector<N, double> Q1mQ0 = Q1 - Q0;
        Vector<N, double> P0mQ0 = P0 - Q0;
        double a = Dot(P1mP0, P1mP0);
        double b = Dot(P1mP0, Q1mQ0);
        double c = Dot(Q1mQ0, Q1mQ0);
        double d = Dot(P1mP0, P0mQ0);
        if (!(a > 0.0 && c > 0.0))
        {
            return false;
        }
        double f00 = d;
        double f10 = f00 + a;
        double f01 = f00 - b;
        double f11 = f10 - b;
        double s0 = ClampedRoot(a, f00, f10);
        double s1 = ClampedRoot(a, f01, f11);
        int32_t c0 = (s0 <= 0.0 ? -1 : (s0 >= 1.0 ? +1 : 0));
        int32_t c1 = (s1 <= 0.0 ? -1 : (s1 >= 1.0 ? +1 : 0));
        if ((c0 == -1 && c1 == -1) || (c0 == +1 && c1 == +1))
        {
            return false;
        }
        auto bad = [](double v) { return !(v >= 0.0 && v <= 1.0); };
        if (c0 < 0)
        {
            if (bad(f00 / b)) { return true; }
            if (c1 != 0 && bad(f10 / b)) { return true; }
        }
        else if (c0 == 0)
        {
            if (c1 < 0 && bad(f00 / b)) { return true; }
            if (c1 > 0 && bad(f10 / b)) { return true; }
        }
        else
        {
            if (bad(f10 / b)) { return true; }
            if (c1 != 0 && bad(f00 / b)) { return true; }
        }
        return false;
    }
}

ORACLE_CASE("DistSegmentSegment.compute.2d") { SegmentSegment<2>(io); }
ORACLE_CASE("DistSegmentSegment.compute.3d") { SegmentSegment<3>(io); }
ORACLE_CASE("DistSegmentSegment.computeRobust.2d") { SegmentSegmentRobust<2>(io); }
ORACLE_CASE("DistSegmentSegment.computeRobust.3d") { SegmentSegmentRobust<3>(io); }

ORACLE_CASE("DistSegmentSegment.computeEndpoints.3d")
{
    // The four-endpoint overload; the same code as the two-segment overload,
    // exercised through the other public entry point.
    int32_t mode = Mode(io);
    Vector3<double> p0{}, p1{}, q0{}, q1{};
    SegmentSegmentInputs<3>(io, mode, p0, p1, q0, q1);
    DCPQuery<double, Segment3<double>, Segment3<double>> query{};
    auto r = query(p0, p1, q0, q1);
    EmitTwoParameters<decltype(r), 3>(io, r);
}

ORACLE_CASE("DistSegmentSegment.computeRobustEndpoints.3d")
{
    int32_t mode = Mode(io);
    Vector3<double> p0{}, p1{}, q0{}, q1{};
    SegmentSegmentInputs<3>(io, mode, p0, p1, q0, q1);
    DCPQuery<double, Segment3<double>, Segment3<double>> query{};
    auto r = query.ComputeRobust(p0, p1, q0, q1);
    EmitTwoParameters<decltype(r), 3>(io, r);
}

ORACLE_CASE("DistSegmentSegment.computeRobust.3d.deviation")
{
    // f10 = b exactly when P1-P0 is perpendicular to P1-Q1, so the line
    // dR/ds = 0 passes through the domain corner (1,1) and the ratio f10/b
    // rounds to one ulp outside [0,1]. Upstream then replaces it by 1/2.
    DCPQuery<double, Segment3<double>, Segment3<double>> query{};
    Vector3<double> p0{}, p1{}, q0{}, q1{};
    bool found = false;
    for (int32_t attempt = 0; attempt < 20000 && !found; ++attempt)
    {
        p0 = RawVec<3>(io, -4.0, 4.0);
        p1 = RawVec<3>(io, -4.0, 4.0);
        q0 = RawVec<3>(io, -4.0, 4.0);
        Vector3<double> u = p1 - p0;
        Vector3<double> w = RawVec<3>(io, -4.0, 4.0);
        double uu = Dot(u, u);
        if (uu <= 0.0) { continue; }
        w = w - (Dot(w, u) / uu) * u;
        q1 = p1 + w;
        found = UsesOutOfRangeRatio<3>(p0, p1, q0, q1);
    }

    io.givenVec(p0);
    io.givenVec(p1);
    io.givenVec(q0);
    io.givenVec(q1);
    Segment3<double> segment0(p0, p1);
    Segment3<double> segment1(q0, q1);
    auto r = query.ComputeRobust(segment0, segment1);
    EmitTwoParameters<decltype(r), 3>(io, r);
}

// =====================================================================
// DistPointTriangle.h
// =====================================================================
namespace
{
    template <int32_t N>
    void TriangleInputs(oracle::Ctx& io, int32_t mode, Vector<N, double>& point,
        std::array<Vector<N, double>, 3>& v)
    {
        point = Point<N>(io, mode);
        v[0] = Point<N>(io, mode);
        v[1] = Point<N>(io, mode);
        if (mode == 4)
        {
            // A degenerate (collinear) triangle.
            double s = io.raw(-2.0, 2.0);
            v[2] = io.givenVec(v[0] + s * (v[1] - v[0]));
        }
        else
        {
            v[2] = Point<N>(io, mode);
        }
    }

    template <int32_t N>
    void PointTriangle(oracle::Ctx& io, bool conjugateGradient)
    {
        int32_t mode = Mode(io);
        Vector<N, double> point{};
        std::array<Vector<N, double>, 3> v{};
        TriangleInputs<N>(io, mode, point, v);
        Triangle<N, double> triangle(v[0], v[1], v[2]);
        DCPQuery<double, Vector<N, double>, Triangle<N, double>> query{};
        auto r = conjugateGradient ? query.UseConjugateGradient(point, triangle)
            : query(point, triangle);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.barycentric[0]);
        io.outReal(r.barycentric[1]);
        io.outReal(r.barycentric[2]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }
}

ORACLE_CASE("DistPointTriangle.compute.2d") { PointTriangle<2>(io, false); }
ORACLE_CASE("DistPointTriangle.compute.3d") { PointTriangle<3>(io, false); }
ORACLE_CASE("DistPointTriangle.useConjugateGradient.2d") { PointTriangle<2>(io, true); }
ORACLE_CASE("DistPointTriangle.useConjugateGradient.3d") { PointTriangle<3>(io, true); }

// =====================================================================
// DistAlignedBoxAlignedBox.h
// =====================================================================
namespace
{
    template <int32_t N>
    void AlignedBoxInputs(oracle::Ctx& io, int32_t mode,
        Vector<N, double>& lo, Vector<N, double>& hi)
    {
        lo = Point<N>(io, mode);
        Vector<N, double> size{};
        if (mode == 2 || mode == 3)
        {
            size = io.latticeVec<N>(0, 3);
        }
        else
        {
            size = io.vec<N>(0.0, 4.0);
        }
        hi = lo + size;
    }

    template <int32_t N>
    void AlignedBoxAlignedBox(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        Vector<N, double> lo0{}, hi0{}, lo1{}, hi1{};
        AlignedBoxInputs<N>(io, mode, lo0, hi0);
        AlignedBoxInputs<N>(io, mode, lo1, hi1);
        AlignedBox<N, double> box0(lo0, hi0);
        AlignedBox<N, double> box1(lo1, hi1);
        DCPQuery<double, AlignedBox<N, double>, AlignedBox<N, double>> query{};
        auto r = query(box0, box1);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0].min);
        io.outVec(r.closest[0].max);
        io.outVec(r.closest[1].min);
        io.outVec(r.closest[1].max);
    }
}

ORACLE_CASE("DistAlignedBoxAlignedBox.compute.2d") { AlignedBoxAlignedBox<2>(io); }
ORACLE_CASE("DistAlignedBoxAlignedBox.compute.3d") { AlignedBoxAlignedBox<3>(io); }

// =====================================================================
// DistCircle2Circle2.h
// =====================================================================
ORACLE_CASE("DistCircle2Circle2.compute")
{
    // mode 0,1: uniform (separated, nested or intersecting)
    // mode 2,3: lattice centers and integer radii, which produce exact
    //           tangency (distance exactly 0) and exact containment
    // mode 4:   concentric, and cocircular when the radii agree
    int32_t mode = Mode(io);
    auto center0 = Point<2>(io, mode);
    double radius0 = (mode == 2 || mode == 3)
        ? static_cast<double>(io.integer(1, 4)) : io.real(0.25, 4.0);
    auto center1 = SecondPoint<2>(io, mode, center0);
    double radius1 = 0.0;
    if (mode == 2 || mode == 3)
    {
        radius1 = static_cast<double>(io.integer(1, 4));
    }
    else if (mode == 4)
    {
        // Concentric; half of these records are also cocircular. The choice is
        // made from unrecorded draws so that the number of recorded doubles
        // does not depend on the mode.
        double alternate = io.raw(0.25, 4.0);
        double pick = io.raw(0.0, 1.0);
        radius1 = io.given(pick < 0.5 ? radius0 : alternate);
    }
    else
    {
        radius1 = io.real(0.25, 4.0);
    }
    Circle2<double> circle0(center0, radius0);
    Circle2<double> circle1(center1, radius1);
    DCPQuery<double, Circle2<double>, Circle2<double>> query{};
    auto r = query(circle0, circle1);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outInt(r.numClosestPairs);
    io.outBool(r.concentric);
    io.outBool(r.cocircular);
    for (size_t j = 0; j < r.numClosestPairs; ++j)
    {
        io.outVec(r.closest[j][0]);
        io.outVec(r.closest[j][1]);
    }
}

// =====================================================================
// DistLine2AlignedBox2.h
// =====================================================================
ORACLE_CASE("DistLine2AlignedBox2.compute")
{
    // mode 2,3 use lattice directions, which have zero components and so
    // reach DoQuery1D; mode 4 uses the zero direction, which reaches
    // DoQuery0D.
    int32_t mode = Mode(io);
    Vector2<double> lo{}, hi{};
    AlignedBoxInputs<2>(io, mode, lo, hi);
    auto origin = Point<2>(io, mode);
    Vector2<double> direction{};
    if (mode == 4)
    {
        direction = io.givenVec(Vector2<double>{ 0.0, 0.0 });
    }
    else
    {
        direction = Direction<2>(io, mode);
    }
    Line2<double> line(origin, direction);
    AlignedBox2<double> box(lo, hi);
    DCPQuery<double, Line2<double>, AlignedBox2<double>> query{};
    auto r = query(line, box);
    EmitOneParameter<decltype(r), 2>(io, r);
}

namespace
{
    // Upstream declares DoQuery protected and grants friendship to the
    // line/oriented-box query; the port exports it as a free function
    // (distLine2AlignedBox2DoQuery), which DistLine2OrientedBox2 calls. A
    // derived class re-publishes it so this case can call it directly.
    struct Line2Box2Access
        : public DCPQuery<double, Line2<double>, AlignedBox2<double>>
    {
        using Base = DCPQuery<double, Line2<double>, AlignedBox2<double>>;
        using Base::DoQuery;
    };
}

ORACLE_CASE("DistLine2AlignedBox2.doQuery")
{
    // DoQuery works in the box frame (the box is centered at the origin) and
    // sets only parameter and closest; distance and sqrDistance are left at
    // their Result{} defaults and are not emitted. The origin and direction
    // are in/out parameters: DoQuery reflects them in place and does not undo
    // the reflection, so the reflected values are emitted as outputs too.
    int32_t mode = Mode(io);
    auto extent = Extent<2>(io, mode);
    auto origin = Point<2>(io, mode);
    Vector2<double> direction{};
    if (mode == 4)
    {
        direction = io.givenVec(Vector2<double>{ 0.0, 0.0 });
    }
    else
    {
        direction = Direction<2>(io, mode);
    }
    Line2Box2Access::Result r{};
    Line2Box2Access::DoQuery(origin, direction, extent, r);
    io.outReal(r.parameter);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
    io.outVec(origin);
    io.outVec(direction);
}

// =====================================================================
// DistLine2Circle2.h
// =====================================================================
ORACLE_CASE("DistLine2Circle2.compute")
{
    // The lattice modes reach exact tangency (test == 0) and the
    // two-intersection branch with exactly representable arithmetic.
    int32_t mode = Mode(io);
    auto center = Point<2>(io, mode);
    double radius = (mode == 2 || mode == 3)
        ? static_cast<double>(io.integer(1, 3)) : io.real(0.25, 4.0);
    auto origin = Point<2>(io, mode);
    auto direction = Direction<2>(io, mode);
    Line2<double> line(origin, direction);
    Circle2<double> circle(center, radius);
    DCPQuery<double, Line2<double>, Circle2<double>> query{};
    auto r = query(line, circle);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outInt(r.numClosestPairs);
    for (size_t j = 0; j < r.numClosestPairs; ++j)
    {
        io.outReal(r.parameter[j]);
        io.outVec(r.closest[j][0]);
        io.outVec(r.closest[j][1]);
    }
}

// =====================================================================
// DistLine2Triangle2.h
// =====================================================================
namespace
{
    void EmitLine2Triangle2(oracle::Ctx& io,
        DCPQuery<double, Line2<double>, Triangle2<double>>::Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outReal(r.barycentric[0]);
        io.outReal(r.barycentric[1]);
        io.outReal(r.barycentric[2]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }
}

ORACLE_CASE("DistLine2Triangle2.compute")
{
    // The lattice modes place vertices exactly on the line (the sign-zero
    // rows of the classification table) and make the triangle degenerate.
    int32_t mode = Mode(io);
    auto origin = Point<2>(io, mode);
    auto direction = Direction<2>(io, mode);
    Vector2<double> v0{}, v1{}, v2{};
    v0 = Point<2>(io, mode);
    v1 = Point<2>(io, mode);
    if (mode == 4)
    {
        // A vertex exactly on the line, so a sign is exactly zero.
        double s = io.raw(-3.0, 3.0);
        v2 = io.givenVec(origin + s * direction);
    }
    else
    {
        v2 = Point<2>(io, mode);
    }
    Line2<double> line(origin, direction);
    Triangle2<double> triangle(v0, v1, v2);
    DCPQuery<double, Line2<double>, Triangle2<double>> query{};
    auto r = query(line, triangle);
    EmitLine2Triangle2(io, r);
}

ORACLE_CASE("DistLine2Triangle2.compute.deviation")
{
    // Two vertices are placed on the line, so their normal components are
    // pure rounding noise. When the signs come out opposite but
    // DotPerp(D, V1-V0) rounds to exactly zero, upstream divides by zero and
    // returns NaN (UPSTREAM-FINDINGS, issue #441). The port falls back to the
    // normal-component difference, which the sign condition guarantees to be
    // nonzero.
    DCPQuery<double, Line2<double>, Triangle2<double>> query{};
    Vector2<double> origin{}, direction{}, v0{}, v1{}, v2{};
    bool found = false;
    for (int32_t attempt = 0; attempt < 20000 && !found; ++attempt)
    {
        origin = RawVec<2>(io, -4.0, 4.0);
        direction = RawUnit<2>(io);
        double a0 = io.raw(-3.0, 3.0);
        double a1 = io.raw(-3.0, 3.0);
        double a2 = io.raw(-3.0, 3.0);
        double off = io.raw(1.0, 4.0);
        v0 = origin + a0 * direction;
        v1 = origin + a1 * direction;
        v2 = origin + a2 * direction + off * Perp(direction);
        auto probe = query(Line2<double>(origin, direction),
            Triangle2<double>(v0, v1, v2));
        found = std::isnan(probe.distance);
    }

    io.givenVec(origin);
    io.givenVec(direction);
    io.givenVec(v0);
    io.givenVec(v1);
    io.givenVec(v2);
    Line2<double> line(origin, direction);
    Triangle2<double> triangle(v0, v1, v2);
    auto r = query(line, triangle);
    EmitLine2Triangle2(io, r);
}

// =====================================================================
// DistLine3CanonicalBox3.h
// =====================================================================
ORACLE_CASE("DistLine3CanonicalBox3.compute")
{
    // The lattice modes give directions with zero components (DoQuery2D and
    // DoQuery1D) and zero extents (a flat box); mode 4 gives the zero
    // direction (DoQuery0D).
    int32_t mode = Mode(io);
    auto extent = Extent<3>(io, mode);
    auto origin = Point<3>(io, mode);
    Vector3<double> direction{};
    if (mode == 4)
    {
        direction = io.givenVec(Vector3<double>{ 0.0, 0.0, 0.0 });
    }
    else
    {
        direction = Direction<3>(io, mode);
    }
    Line3<double> line(origin, direction);
    CanonicalBox3<double> box(extent);
    DCPQuery<double, Line3<double>, CanonicalBox3<double>> query{};
    auto r = query(line, box);
    EmitOneParameter<decltype(r), 3>(io, r);
}

ORACLE_CASE("DistLine3CanonicalBox3.compute.deviation")
{
    // A line through a point of a flat box. The squared distance is
    // accumulated as pme^2 + tmp^2 + PpE^2 + delta*parameter, whose terms
    // cancel; the rounding can make it slightly negative and upstream then
    // returns distance = NaN (UPSTREAM-FINDINGS, issue #421). The port clamps
    // the squared distance at zero.
    DCPQuery<double, Line3<double>, CanonicalBox3<double>> query{};
    Vector3<double> extent{}, origin{}, direction{};
    bool found = false;
    for (int32_t attempt = 0; attempt < 20000 && !found; ++attempt)
    {
        extent = RawVec<3>(io, 0.25, 3.0);
        extent[io.rawInteger(0, 2)] = 0.0;
        Vector3<double> surface{};
        for (int32_t i = 0; i < 3; ++i)
        {
            surface[i] = io.raw(-extent[i], extent[i]);
        }
        int32_t face = io.rawInteger(0, 2);
        surface[face] = (io.rawInteger(0, 1) == 0 ? -extent[face] : extent[face]);
        direction = RawUnit<3>(io);
        double t = io.raw(0.5, 4.0);
        origin = surface - t * direction;
        auto probe = query(Line3<double>(origin, direction),
            CanonicalBox3<double>(extent));
        found = (probe.sqrDistance < 0.0);
    }

    io.givenVec(extent);
    io.givenVec(origin);
    io.givenVec(direction);
    Line3<double> line(origin, direction);
    CanonicalBox3<double> box(extent);
    auto r = query(line, box);
    EmitOneParameter<decltype(r), 3>(io, r);
}
