// Verify group 33 (intersection): differential cases for the Intr* headers
// listed in plan/verify-groups.json group 33.
//
// Every generator mixes an exactly representable small-lattice mode with a
// uniform mode and, where a branch needs it, a constructed mode, so that the
// touching, tangent, parallel, coplanar and degenerate branches are reached.
// Almost every query in this group uses only + - * / sqrt fabs min max and
// comparisons, so almost every case is declared exact on the TypeScript side;
// the sin/cos used to build orthonormal frames and unit directions are applied
// to *unrecorded* draws and only the resulting frame or unit vector is
// recorded as an input, so libm never enters the compared computation. The one
// exception is IntrEllipsoid3Ellipsoid3, whose eigensolver and bisection
// brackets call the C math library; it is marked in its case comments.
#define ORACLE_FAMILY "v33-intersection"
#include "Oracle.h"

#include <Mathematics/IntrEllipsoid3Ellipsoid3.h>
#include <Mathematics/IntrHalfspace3Ellipsoid3.h>
#include <Mathematics/IntrLine2Arc2.h>
#include <Mathematics/IntrLine2SegmentMesh2.h>
#include <Mathematics/IntrLine3Ellipsoid3.h>
#include <Mathematics/IntrOrientedBox3Cylinder3.h>
#include <Mathematics/IntrPlane3Plane3.h>
#include <Mathematics/IntrPlane3Triangle3.h>
#include <Mathematics/IntrRay2Circle2.h>
#include <Mathematics/IntrRay2OrientedBox2.h>
#include <Mathematics/IntrRay3Capsule3.h>
#include <Mathematics/IntrRay3Cylinder3.h>
#include <Mathematics/IntrRay3OrientedBox3.h>
#include <Mathematics/IntrSegment2Circle2.h>
#include <Mathematics/IntrSegment2OrientedBox2.h>
#include <Mathematics/IntrSegment3Capsule3.h>
#include <Mathematics/IntrSegment3Cylinder3.h>
#include <Mathematics/IntrSegment3OrientedBox3.h>
#include <Mathematics/IntrSphere3Frustum3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
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

    // A positive extent vector; lattice mode gives exactly representable
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

    // A signed coordinate axis (exact unit vector).
    template <int N>
    Vector<N, double> AxisUnit(oracle::Ctx& io)
    {
        int k = io.rawInteger(0, 2 * N - 1);
        Vector<N, double> v{};
        v.MakeZero();
        v[k % N] = (k < N ? 1.0 : -1.0);
        return io.givenVec(v);
    }

    // Orthonormal 2D frame. Mode 0 produces exact coordinate axes (with signs
    // and a swap) so that the parallel-axis branches are reached exactly.
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

    OrientedBox2<double> OBox2(oracle::Ctx& io, int mode, int lat, double range)
    {
        OrientedBox2<double> box{};
        box.center = Pt<2>(io, mode, lat, range);
        auto axis = Frame2(io, mode);
        box.axis[0] = axis[0];
        box.axis[1] = axis[1];
        box.extent = Extent<2>(io, mode, lat, range);
        return box;
    }

    OrientedBox3<double> OBox3(oracle::Ctx& io, int mode, int lat, double range)
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

    template <int N>
    Triangle<N, double> Tri(oracle::Ctx& io, int mode, int lat, double range)
    {
        Triangle<N, double> triangle{};
        triangle.v[0] = Pt<N>(io, mode, lat, range);
        triangle.v[1] = Pt<N>(io, mode, lat, range);
        triangle.v[2] = Pt<N>(io, mode, lat, range);
        return triangle;
    }

    // A line with a unit-length direction, for the queries that document that
    // precondition (circle, ellipsoid, capsule, cylinder).
    template <int N>
    Line<N, double> UnitLn(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> origin = Pt<N>(io, mode, lat, range);
        Vector<N, double> direction = (mode == 0 ? AxisUnit<N>(io) : io.unit<N>());
        return Line<N, double>(origin, direction);
    }

    // A ray whose direction is not required to be unit length (the box
    // queries clip parametrically and do not need a unit direction).
    template <int N>
    Ray<N, double> Ry(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> origin = Pt<N>(io, mode, lat, range);
        Vector<N, double> direction = (mode == 0 ? io.latticeDir<N>(-2, 2) : io.unit<N>());
        return Ray<N, double>(origin, direction);
    }

    template <int N>
    Ray<N, double> UnitRy(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> origin = Pt<N>(io, mode, lat, range);
        Vector<N, double> direction = (mode == 0 ? AxisUnit<N>(io) : io.unit<N>());
        return Ray<N, double>(origin, direction);
    }

    template <int N>
    Segment<N, double> Sg(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> p0 = Pt<N>(io, mode, lat, range);
        Vector<N, double> p1 = Pt<N>(io, mode, lat, range);
        return Segment<N, double>(p0, p1);
    }

    Circle2<double> Circ(oracle::Ctx& io, int mode, int lat, double range)
    {
        Circle2<double> circle{};
        circle.center = Pt<2>(io, mode, lat, range);
        circle.radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        return circle;
    }

    Cylinder3<double> Cyl(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector3<double> origin = Pt<3>(io, mode, lat, range);
        Vector3<double> direction = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
        double radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        double height = (mode == 0 ? io.lattice(1, 2 * lat) : io.real(0.5, 2.0 * range));
        Line3<double> axis(origin, direction);
        return Cylinder3<double>(axis, radius, height);
    }

    // Extents that are powers of two in lattice mode, so that 1/e and e*(1/e)
    // are exact and the ellipsoid matrices are exactly representable.
    template <int N>
    Vector<N, double> ExtentPow2(oracle::Ctx& io, int mode, double range)
    {
        Vector<N, double> e{};
        for (int i = 0; i < N; ++i)
        {
            e[i] = (mode == 0 ? io.given(std::ldexp(1.0, io.rawInteger(-1, 2)))
                : io.real(0.5, range));
        }
        return e;
    }

    Ellipsoid3<double> Ellip(oracle::Ctx& io, int mode, int lat, double range)
    {
        Ellipsoid3<double> ellipsoid{};
        ellipsoid.center = Pt<3>(io, mode, lat, range);
        auto axis = Frame3(io, mode);
        ellipsoid.axis[0] = axis[0];
        ellipsoid.axis[1] = axis[1];
        ellipsoid.axis[2] = axis[2];
        ellipsoid.extent = ExtentPow2<3>(io, mode, range);
        return ellipsoid;
    }

    // Exposes the protected DoQuery helpers, which the port exports as public
    // free functions.
    template <typename Q>
    struct Expose : public Q
    {
        using Q::DoQuery;
    };

    // The twelve lattice points on the circle of radius 5 centred at the
    // origin, in counterclockwise order. Arc endpoints drawn from this table
    // lie exactly on the circle, which is what Arc2::Contains assumes, and the
    // line/circle roots through them are exact integers.
    Vector2<double> CirclePoint5(int k)
    {
        static int const p[12][2] =
        {
            { 5, 0 }, { 4, 3 }, { 3, 4 }, { 0, 5 }, { -3, 4 }, { -4, 3 },
            { -5, 0 }, { -4, -3 }, { -3, -4 }, { 0, -5 }, { 3, -4 }, { 4, -3 }
        };
        return Vector2<double>{ static_cast<double>(p[k][0]),
            static_cast<double>(p[k][1]) };
    }

    // An arc whose endpoints lie on its circle (the precondition of
    // Arc2::Contains). Mode 0 puts the arc on the lattice circle of radius 5,
    // where the endpoints and the side-of-line test are exact.
    Arc2<double> Arc(oracle::Ctx& io, int mode)
    {
        Arc2<double> arc{};
        if (mode == 0)
        {
            Vector2<double> center{
                static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            arc.center = io.givenVec(center);
            arc.radius = io.given(5.0);
            arc.end[0] = io.givenVec(center + CirclePoint5(io.rawInteger(0, 11)));
            arc.end[1] = io.givenVec(center + CirclePoint5(io.rawInteger(0, 11)));
        }
        else
        {
            Vector2<double> center{ io.raw(-3.0, 3.0), io.raw(-3.0, 3.0) };
            double radius = io.raw(1.0, 4.0);
            double a0 = io.raw(-3.141592653589793, 3.141592653589793);
            double a1 = io.raw(-3.141592653589793, 3.141592653589793);
            arc.center = io.givenVec(center);
            arc.radius = io.given(radius);
            arc.end[0] = io.givenVec(center
                + radius * Vector2<double>{ std::cos(a0), std::sin(a0) });
            arc.end[1] = io.givenVec(center
                + radius * Vector2<double>{ std::cos(a1), std::sin(a1) });
        }
        return arc;
    }
}

// =============================== IntrLine2Arc2 ===========================

ORACLE_CASE("IntrLine2Arc2.test")
{
    int mode = io.index() % 2;
    auto arc = Arc(io, mode);
    auto line = UnitLn<2>(io, mode, 5, 6.0);
    TIQuery<double, Line2<double>, Arc2<double>> query;
    auto r = query(line, arc);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine2Arc2.find")
{
    int mode = io.index() % 2;
    auto arc = Arc(io, mode);
    auto line = UnitLn<2>(io, mode, 5, 6.0);
    FIQuery<double, Line2<double>, Arc2<double>> query;
    auto r = query(line, arc);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

// The line passes exactly through lattice points of the arc's circle. With the
// circle of radius 5 centred at an integer point, an axis-parallel unit
// direction and a lattice offset 'a' in {0,+-3,+-4,+-5}, the discriminant is
// exactly 25 - a^2 and the two roots are integers, so the intersection points
// are exactly the lattice circle points that the arc endpoints are drawn from
// and Arc2::Contains evaluates DotPerp of integer vectors (exactly zero when
// the hit is an endpoint). a = +-5 gives the tangent branch.
ORACLE_CASE("IntrLine2Arc2.find.throughLatticePoint")
{
    Vector2<double> center{
        static_cast<double>(io.rawInteger(-3, 3)),
        static_cast<double>(io.rawInteger(-3, 3)) };
    Arc2<double> arc{};
    arc.center = io.givenVec(center);
    arc.radius = io.given(5.0);
    arc.end[0] = io.givenVec(center + CirclePoint5(io.rawInteger(0, 11)));
    arc.end[1] = io.givenVec(center + CirclePoint5(io.rawInteger(0, 11)));

    static int const offsets[7] = { 0, 3, -3, 4, -4, 5, -5 };
    double a = static_cast<double>(offsets[io.rawInteger(0, 6)]);
    double t = static_cast<double>(io.rawInteger(-5, 5));
    int axis = io.rawInteger(0, 1);
    Vector2<double> direction{}, origin{};
    direction.MakeZero();
    direction[axis] = 1.0;
    origin[axis] = center[axis] - t;
    origin[1 - axis] = center[1 - axis] + a;
    io.givenVec(origin);
    io.givenVec(direction);
    Line2<double> line(origin, direction);
    FIQuery<double, Line2<double>, Arc2<double>> query;
    auto r = query(line, arc);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

// ============================== IntrRay2Circle2 ==========================

ORACLE_CASE("IntrRay2Circle2.test")
{
    int mode = io.index() % 2;
    auto ray = UnitRy<2>(io, mode, 4, 5.0);
    auto circle = Circ(io, mode, 3, 4.0);
    TIQuery<double, Ray2<double>, Circle2<double>> query;
    auto r = query(ray, circle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay2Circle2.find")
{
    int mode = io.index() % 2;
    auto ray = UnitRy<2>(io, mode, 4, 5.0);
    auto circle = Circ(io, mode, 3, 4.0);
    FIQuery<double, Ray2<double>, Circle2<double>> query;
    auto r = query(ray, circle);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

ORACLE_CASE("IntrRay2Circle2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto rayOrigin = Pt<2>(io, mode, 4, 5.0);
    auto rayDirection = (mode == 0 ? AxisUnit<2>(io) : io.unit<2>());
    auto circle = Circ(io, mode, 3, 4.0);
    Expose<FIQuery<double, Ray2<double>, Circle2<double>>> query;
    FIQuery<double, Ray2<double>, Circle2<double>>::Result r{};
    query.DoQuery(rayOrigin, rayDirection, circle, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// The exact construction: the circle of radius 5 centred at an integer point,
// an axis-parallel unit direction and a lattice offset 'a' perpendicular to
// it. The discriminant is exactly 25 - a^2, so the roots are the integers
// t -+ sqrt(25-a^2); a = +-5 makes them coincide (tangency) and t = +-root
// puts a root exactly at the ray origin, where the [0,+infinity) clip is
// evaluated at exact equality.
ORACLE_CASE("IntrRay2Circle2.find.latticeTangent")
{
    Vector2<double> center{
        static_cast<double>(io.rawInteger(-3, 3)),
        static_cast<double>(io.rawInteger(-3, 3)) };
    Circle2<double> circle{};
    circle.center = io.givenVec(center);
    circle.radius = io.given(5.0);
    static int const offsets[7] = { 0, 3, -3, 4, -4, 5, -5 };
    double a = static_cast<double>(offsets[io.rawInteger(0, 6)]);
    double t = static_cast<double>(io.rawInteger(-5, 5));
    int axis = io.rawInteger(0, 1);
    Vector2<double> direction{}, origin{};
    direction.MakeZero();
    direction[axis] = 1.0;
    origin[axis] = center[axis] - t;
    origin[1 - axis] = center[1 - axis] + a;
    io.givenVec(origin);
    io.givenVec(direction);
    Ray2<double> ray(origin, direction);
    FIQuery<double, Ray2<double>, Circle2<double>> query;
    auto r = query(ray, circle);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

// =========================== IntrSegment2Circle2 =========================

ORACLE_CASE("IntrSegment2Circle2.test")
{
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 4, 5.0);
    auto circle = Circ(io, mode, 3, 4.0);
    TIQuery<double, Segment2<double>, Circle2<double>> query;
    auto r = query(segment, circle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment2Circle2.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 4, 5.0);
    auto circle = Circ(io, mode, 3, 4.0);
    FIQuery<double, Segment2<double>, Circle2<double>> query;
    auto r = query(segment, circle);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

ORACLE_CASE("IntrSegment2Circle2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto segOrigin = Pt<2>(io, mode, 4, 5.0);
    auto segDirection = (mode == 0 ? AxisUnit<2>(io) : io.unit<2>());
    double segExtent = (mode == 0 ? io.lattice(1, 5) : io.real(0.25, 5.0));
    auto circle = Circ(io, mode, 3, 4.0);
    Expose<FIQuery<double, Segment2<double>, Circle2<double>>> query;
    FIQuery<double, Segment2<double>, Circle2<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, circle, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// The same exact construction as IntrRay2Circle2.find.latticeTangent, with the
// segment endpoints placed symmetrically about the origin so that the centred
// form (centre, unit direction, integer extent) is exact and the [-e,e] clip
// is evaluated at exact equality.
ORACLE_CASE("IntrSegment2Circle2.find.latticeTangent")
{
    Vector2<double> center{
        static_cast<double>(io.rawInteger(-3, 3)),
        static_cast<double>(io.rawInteger(-3, 3)) };
    Circle2<double> circle{};
    circle.center = io.givenVec(center);
    circle.radius = io.given(5.0);
    static int const offsets[7] = { 0, 3, -3, 4, -4, 5, -5 };
    double a = static_cast<double>(offsets[io.rawInteger(0, 6)]);
    double t = static_cast<double>(io.rawInteger(-5, 5));
    double e = static_cast<double>(io.rawInteger(1, 6));
    int axis = io.rawInteger(0, 1);
    Vector2<double> direction{}, origin{};
    direction.MakeZero();
    direction[axis] = 1.0;
    origin[axis] = center[axis] - t;
    origin[1 - axis] = center[1 - axis] + a;
    Vector2<double> p0 = origin - e * direction;
    Vector2<double> p1 = origin + e * direction;
    io.givenVec(p0);
    io.givenVec(p1);
    Segment2<double> segment(p0, p1);
    FIQuery<double, Segment2<double>, Circle2<double>> query;
    auto r = query(segment, circle);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

namespace
{
    // A lattice oriented box: integer centre and extents, axes a signed
    // permutation of the coordinate axes. Every world/box change of basis is
    // then exact, so a dyadic point of the box boundary stays dyadic in world
    // coordinates.
    template <int N>
    struct LatticeBox
    {
        Vector<N, double> center{}, extent{};
        std::array<Vector<N, double>, N> axis{};
    };

    LatticeBox<2> DrawLatticeBox2(oracle::Ctx& io)
    {
        LatticeBox<2> box{};
        box.center = io.latticeVec<2>(-3, 3);
        auto axis = Frame2(io, 0);
        box.axis[0] = axis[0];
        box.axis[1] = axis[1];
        box.extent = Extent<2>(io, 0, 3, 3.0);
        return box;
    }

    LatticeBox<3> DrawLatticeBox3(oracle::Ctx& io)
    {
        LatticeBox<3> box{};
        box.center = io.latticeVec<3>(-3, 3);
        auto axis = Frame3(io, 0);
        box.axis[0] = axis[0];
        box.axis[1] = axis[1];
        box.axis[2] = axis[2];
        box.extent = Extent<3>(io, 0, 3, 3.0);
        return box;
    }

    // A point of the boundary of 'box': one coordinate is +-E[f], the others
    // are (k/2)*E[i] with integer k in [-2,2]. The world point is exact.
    template <int N>
    Vector<N, double> BoundaryPoint(oracle::Ctx& io, LatticeBox<N> const& box)
    {
        int f = io.rawInteger(0, N - 1);
        Vector<N, double> y{};
        for (int i = 0; i < N; ++i)
        {
            y[i] = (i == f ? (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0) * box.extent[i]
                : 0.5 * static_cast<double>(io.rawInteger(-2, 2)) * box.extent[i]);
        }
        Vector<N, double> point = box.center;
        for (int i = 0; i < N; ++i) { point = point + y[i] * box.axis[i]; }
        return point;
    }
}

// ========================== IntrRay2OrientedBox2 =========================

ORACLE_CASE("IntrRay2OrientedBox2.test")
{
    int mode = io.index() % 2;
    auto ray = Ry<2>(io, mode, 4, 5.0);
    auto box = OBox2(io, mode, 3, 4.0);
    TIQuery<double, Ray2<double>, OrientedBox2<double>> query;
    auto r = query(ray, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay2OrientedBox2.find")
{
    int mode = io.index() % 2;
    auto ray = Ry<2>(io, mode, 4, 5.0);
    auto box = OBox2(io, mode, 3, 4.0);
    FIQuery<double, Ray2<double>, OrientedBox2<double>> query;
    auto r = query(ray, box);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

// The ray is aimed at a dyadic point of the box boundary with an integer
// direction, so the Liang-Barsky clip evaluates its '<=' tests at exact
// equality and edge and corner contacts occur.
ORACLE_CASE("IntrRay2OrientedBox2.find.throughPoint")
{
    auto box = DrawLatticeBox2(io);
    Vector2<double> target = BoundaryPoint<2>(io, box);
    Vector2<double> rawDirection{
        static_cast<double>(io.rawInteger(-2, 2)),
        static_cast<double>(io.rawInteger(-2, 2)) };
    if (rawDirection[0] == 0.0 && rawDirection[1] == 0.0) { rawDirection[0] = 1.0; }
    double t = 0.5 * static_cast<double>(io.rawInteger(-4, 4));
    Vector2<double> origin = io.givenVec(target - t * rawDirection);
    Vector2<double> direction = io.givenVec(rawDirection);
    OrientedBox2<double> obox(box.center, box.axis, box.extent);
    Ray2<double> ray(origin, direction);
    FIQuery<double, Ray2<double>, OrientedBox2<double>> query;
    auto r = query(ray, obox);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

// ========================== IntrRay3OrientedBox3 =========================

ORACLE_CASE("IntrRay3OrientedBox3.test")
{
    int mode = io.index() % 2;
    auto ray = Ry<3>(io, mode, 4, 5.0);
    auto box = OBox3(io, mode, 3, 4.0);
    TIQuery<double, Ray3<double>, OrientedBox3<double>> query;
    auto r = query(ray, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3OrientedBox3.find")
{
    int mode = io.index() % 2;
    auto ray = Ry<3>(io, mode, 4, 5.0);
    auto box = OBox3(io, mode, 3, 4.0);
    FIQuery<double, Ray3<double>, OrientedBox3<double>> query;
    auto r = query(ray, box);
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

ORACLE_CASE("IntrRay3OrientedBox3.find.throughPoint")
{
    auto box = DrawLatticeBox3(io);
    Vector3<double> target = BoundaryPoint<3>(io, box);
    Vector3<double> rawDirection{};
    bool zero = true;
    do
    {
        zero = true;
        for (int i = 0; i < 3; ++i)
        {
            rawDirection[i] = static_cast<double>(io.rawInteger(-2, 2));
            zero = zero && rawDirection[i] == 0.0;
        }
    } while (zero);
    double t = 0.5 * static_cast<double>(io.rawInteger(-4, 4));
    Vector3<double> origin = io.givenVec(target - t * rawDirection);
    Vector3<double> direction = io.givenVec(rawDirection);
    OrientedBox3<double> obox(box.center, box.axis, box.extent);
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, OrientedBox3<double>> query;
    auto r = query(ray, obox);
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

// ======================== IntrSegment3OrientedBox3 =======================

ORACLE_CASE("IntrSegment3OrientedBox3.test")
{
    int mode = io.index() % 2;
    auto segment = Sg<3>(io, mode, 4, 5.0);
    auto box = OBox3(io, mode, 3, 4.0);
    TIQuery<double, Segment3<double>, OrientedBox3<double>> query;
    auto r = query(segment, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3OrientedBox3.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<3>(io, mode, 4, 5.0);
    auto box = OBox3(io, mode, 3, 4.0);
    FIQuery<double, Segment3<double>, OrientedBox3<double>> query;
    auto r = query(segment, box);
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

// The segment is aimed at a dyadic point of the box boundary along a signed
// coordinate axis, so the centred form (centre, unit direction, extent) and
// the box-frame coordinates are exact and the Liang-Barsky '<=' tests are
// evaluated at exact equality.
ORACLE_CASE("IntrSegment3OrientedBox3.find.throughPoint")
{
    auto box = DrawLatticeBox3(io);
    Vector3<double> target = BoundaryPoint<3>(io, box);
    int k = io.rawInteger(0, 5);
    Vector3<double> direction{};
    direction.MakeZero();
    direction[k % 3] = (k < 3 ? 1.0 : -1.0);
    // a = 0 puts an endpoint exactly on the boundary point, which is where the
    // clip collapses to a single contact (numIntersections == 1).
    double a = static_cast<double>(io.rawInteger(0, 4));
    double b = static_cast<double>(io.rawInteger(1, 4));
    Vector3<double> p0 = io.givenVec(target - a * direction);
    Vector3<double> p1 = io.givenVec(target + b * direction);
    OrientedBox3<double> obox(box.center, box.axis, box.extent);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, OrientedBox3<double>> query;
    auto r = query(segment, obox);
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

// ======================== IntrSegment2OrientedBox2 =======================
//
// Upstream's FIQuery builds the intersection point as
//   result.point[i] = box.center + (segOrigin + parameter[i] * segDirection)
// where segOrigin and segDirection are box-frame vectors and box.center is
// world-space. That is wrong for any box whose axes are not the coordinate
// axes; the port evaluates the world-space centred form instead
// (docs/UPSTREAM-FINDINGS.md IntrSegment2OrientedBox2.h FIQuery, issue #255).
// Even for an axis-aligned box the two expressions are only *mathematically*
// equal: upstream adds and subtracts box.center, which rounds. The set on
// which upstream's expression is bit-identical to the port's is therefore
// box.center == 0 with the identity axes, which is what
// 'find.originBox' uses. 'find' covers every other field of the result on a
// general rotated box, 'find.degenerate' covers the zero-extent branch (whose
// points are the segment endpoints, with no change of basis), and
// 'find.frameDeviation' demonstrates the fix.

namespace
{
    void OutSeg2Box2FI(oracle::Ctx& io,
        FIQuery<double, Segment2<double>, OrientedBox2<double>>::Result const& r,
        bool withPoints)
    {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outReal(r.cdeParameter[0]);
        io.outReal(r.cdeParameter[1]);
        if (withPoints)
        {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }
}

ORACLE_CASE("IntrSegment2OrientedBox2.test")
{
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 4, 5.0);
    auto box = OBox2(io, mode, 3, 4.0);
    TIQuery<double, Segment2<double>, OrientedBox2<double>> query;
    auto r = query(segment, box);
    io.outBool(r.intersect);
}

// Every field but 'point'; see the header comment above.
ORACLE_CASE("IntrSegment2OrientedBox2.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 4, 5.0);
    auto box = OBox2(io, mode, 3, 4.0);
    FIQuery<double, Segment2<double>, OrientedBox2<double>> query;
    auto r = query(segment, box);
    OutSeg2Box2FI(io, r, false);
}

// The segment is aimed at a dyadic point of a lattice box's boundary along a
// signed coordinate axis, so the centred form and the box-frame coordinates
// are exact; a = 0 puts an endpoint exactly on the boundary point, which is
// where the clip collapses to a single contact. Every field but 'point'.
ORACLE_CASE("IntrSegment2OrientedBox2.find.throughPoint")
{
    auto box = DrawLatticeBox2(io);
    Vector2<double> target = BoundaryPoint<2>(io, box);
    int k = io.rawInteger(0, 3);
    Vector2<double> direction{};
    direction.MakeZero();
    direction[k % 2] = (k < 2 ? 1.0 : -1.0);
    double a = static_cast<double>(io.rawInteger(0, 4));
    double b = static_cast<double>(io.rawInteger(1, 4));
    Vector2<double> p0 = io.givenVec(target - a * direction);
    Vector2<double> p1 = io.givenVec(target + b * direction);
    OrientedBox2<double> obox(box.center, box.axis, box.extent);
    Segment2<double> segment(p0, p1);
    FIQuery<double, Segment2<double>, OrientedBox2<double>> query;
    auto r = query(segment, obox);
    OutSeg2Box2FI(io, r, false);
}

// A box centred at the origin with the identity axes: upstream's expression
// box.center + (segOrigin + t * segDirection) is then bit-identical to the
// port's tmpOrigin + t * tmpDirection, so every field including 'point' is
// compared.
ORACLE_CASE("IntrSegment2OrientedBox2.find.originBox")
{
    OrientedBox2<double> box{};
    box.center = io.givenVec(Vector2<double>{ 0.0, 0.0 });
    box.axis[0] = io.givenVec(Vector2<double>{ 1.0, 0.0 });
    box.axis[1] = io.givenVec(Vector2<double>{ 0.0, 1.0 });
    box.extent = Extent<2>(io, 0, 3, 3.0);
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 4, 5.0);
    FIQuery<double, Segment2<double>, OrientedBox2<double>> query;
    auto r = query(segment, box);
    OutSeg2Box2FI(io, r, true);
}

// The zero-extent branch: a degenerate (point) segment, reported as contained
// or not. The points upstream stores are the segment endpoints themselves, so
// no change of basis is involved and the branch is compared in full.
ORACLE_CASE("IntrSegment2OrientedBox2.find.degenerate")
{
    int mode = io.index() % 2;
    auto box = OBox2(io, mode, 3, 4.0);
    // The point is placed in box coordinates so that inside and outside occur
    // about equally often and the |coordinate| <= extent test is reached at
    // exact equality in the lattice mode.
    double c0 = 0.5 * static_cast<double>(io.rawInteger(-3, 3));
    double c1 = 0.5 * static_cast<double>(io.rawInteger(-3, 3));
    Vector2<double> p = io.givenVec(box.center + (c0 * box.extent[0]) * box.axis[0]
        + (c1 * box.extent[1]) * box.axis[1]);
    Segment2<double> segment(p, p);
    FIQuery<double, Segment2<double>, OrientedBox2<double>> query;
    auto r = query(segment, box);
    OutSeg2Box2FI(io, r, true);
}

// Deliberate deviation: a rotated box, where upstream adds a box-frame vector
// to the world-space box centre. docs/UPSTREAM-FINDINGS.md
// IntrSegment2OrientedBox2.h FIQuery, issue #255.
ORACLE_CASE("IntrSegment2OrientedBox2.find.frameDeviation")
{
    auto box = OBox2(io, 1, 3, 4.0);
    // A segment through the box centre, so that every record has an
    // intersection and therefore a reported point.
    double angle = io.raw(-3.141592653589793, 3.141592653589793);
    Vector2<double> u{ std::cos(angle), std::sin(angle) };
    double halfLength = io.raw(6.0, 10.0);
    Vector2<double> p0 = io.givenVec(box.center + halfLength * u);
    Vector2<double> p1 = io.givenVec(box.center - halfLength * u);
    Segment2<double> segment(p0, p1);
    FIQuery<double, Segment2<double>, OrientedBox2<double>> query;
    auto r = query(segment, box);
    OutSeg2Box2FI(io, r, true);
}

// ============================= IntrPlane3Plane3 ==========================

namespace
{
    // A plane with a unit-length normal, the documented precondition of both
    // queries. Mode 0 uses a signed coordinate axis and an integer constant,
    // which makes |Dot(N0,N1)| exactly 0 or 1 and the coplanarity test exact.
    Plane3<double> Pln(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector3<double> normal = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
        double constant = (mode == 0 ? io.lattice(-lat, lat) : io.real(-range, range));
        return Plane3<double>(normal, constant);
    }

    void OutPlane3Plane3FI(oracle::Ctx& io,
        FIQuery<double, Plane3<double>, Plane3<double>>::Result const& r)
    {
        io.outBool(r.intersect);
        io.outBool(r.isLine);
        if (r.intersect)
        {
            if (r.isLine)
            {
                io.outVec(r.line.origin);
                io.outVec(r.line.direction);
            }
            else
            {
                io.outVec(r.plane.normal);
                io.outReal(r.plane.constant);
            }
        }
    }
}

ORACLE_CASE("IntrPlane3Plane3.test")
{
    int mode = io.index() % 2;
    auto plane0 = Pln(io, mode, 3, 4.0);
    auto plane1 = Pln(io, mode, 3, 4.0);
    TIQuery<double, Plane3<double>, Plane3<double>> query;
    auto r = query(plane0, plane1);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrPlane3Plane3.find")
{
    int mode = io.index() % 2;
    auto plane0 = Pln(io, mode, 3, 4.0);
    auto plane1 = Pln(io, mode, 3, 4.0);
    FIQuery<double, Plane3<double>, Plane3<double>> query;
    auto r = query(plane0, plane1);
    OutPlane3Plane3FI(io, r);
}

// The second normal is exactly +-the first, so |Dot(N0,N1)| is exactly 1 and
// the coplanar / parallel-but-distinct branches are both reached (the constant
// is either the matching one or a perturbation of it). A uniform unit normal
// is used, which is the regime of the documented upstream quirk that
// |Dot(N,N)| is only about half the time exactly 1 (issue #465): here both
// planes share the same normal vector, so the dot product is Dot(N,N) itself.
ORACLE_CASE("IntrPlane3Plane3.find.parallel")
{
    Vector3<double> normal0 = io.unit<3>();
    double constant0 = io.real(-4.0, 4.0);
    double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    int kind = io.rawInteger(0, 2);
    double delta = (kind == 0 ? 0.0 : (kind == 1 ? io.raw(-1.0, 1.0) : 1.0));
    Vector3<double> normal1 = io.givenVec(sign * normal0);
    double constant1 = io.given(sign * constant0 + delta);
    Plane3<double> plane0(normal0, constant0);
    Plane3<double> plane1(normal1, constant1);
    FIQuery<double, Plane3<double>, Plane3<double>> query;
    auto r = query(plane0, plane1);
    OutPlane3Plane3FI(io, r);
    TIQuery<double, Plane3<double>, Plane3<double>> tiQuery;
    io.outBool(tiQuery(plane0, plane1).intersect);
}

// =========================== IntrPlane3Triangle3 =========================

namespace
{
    void OutPlane3Triangle3FI(oracle::Ctx& io,
        FIQuery<double, Plane3<double>, Triangle3<double>>::Result const& r)
    {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outBool(r.isInterior);
        // Every point is emitted on every path; the entries upstream does not
        // write keep their default-constructed zero value, which the port's
        // default result also has.
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
        io.outVec(r.point[2]);
    }
}

// The lattice mode puts the plane normal on a coordinate axis with an integer
// constant and the triangle on the integer lattice, so every signed distance
// s[i] = Dot(N,V[i]) - c is an exact integer and the numZero = 1, 2 and 3
// branches all occur.
ORACLE_CASE("IntrPlane3Triangle3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln(io, mode, 3, 4.0);
    auto triangle = Tri<3>(io, mode, 3, 4.0);
    TIQuery<double, Plane3<double>, Triangle3<double>> query;
    auto r = query(plane, triangle);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outBool(r.isInterior);
}

ORACLE_CASE("IntrPlane3Triangle3.find")
{
    int mode = io.index() % 2;
    auto plane = Pln(io, mode, 3, 4.0);
    auto triangle = Tri<3>(io, mode, 3, 4.0);
    FIQuery<double, Plane3<double>, Triangle3<double>> query;
    auto r = query(plane, triangle);
    OutPlane3Triangle3FI(io, r);
}

// Vertices constructed to make a chosen number of signed distances exactly
// zero: the plane is the coordinate plane x_k = c and each vertex either has
// x_k = c exactly or is offset by a nonzero integer. All ten (n, p, z) rows of
// upstream's table are reachable.
ORACLE_CASE("IntrPlane3Triangle3.find.onPlane")
{
    int k = io.rawInteger(0, 2);
    Vector3<double> normal{};
    normal.MakeZero();
    normal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(normal);
    double constant = io.lattice(-3, 3);
    Triangle3<double> triangle{};
    for (int i = 0; i < 3; ++i)
    {
        Vector3<double> v{};
        for (int j = 0; j < 3; ++j)
        {
            v[j] = static_cast<double>(io.rawInteger(-3, 3));
        }
        // Offset 0 puts the vertex exactly on the plane.
        double offset = static_cast<double>(io.rawInteger(-1, 1));
        v[k] = normal[k] * (constant + offset);
        triangle.v[i] = io.givenVec(v);
    }
    Plane3<double> plane(normal, constant);
    FIQuery<double, Plane3<double>, Triangle3<double>> query;
    auto r = query(plane, triangle);
    OutPlane3Triangle3FI(io, r);
    TIQuery<double, Plane3<double>, Triangle3<double>> tiQuery;
    auto t = tiQuery(plane, triangle);
    io.outBool(t.intersect);
    io.outInt(t.numIntersections);
    io.outBool(t.isInterior);
}

// =========================== IntrLine3Ellipsoid3 =========================

ORACLE_CASE("IntrLine3Ellipsoid3.test")
{
    int mode = io.index() % 2;
    auto line = UnitLn<3>(io, mode, 3, 4.0);
    auto ellipsoid = Ellip(io, mode, 2, 3.0);
    TIQuery<double, Line3<double>, Ellipsoid3<double>> query;
    auto r = query(line, ellipsoid);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3Ellipsoid3.find")
{
    int mode = io.index() % 2;
    auto line = UnitLn<3>(io, mode, 3, 4.0);
    auto ellipsoid = Ellip(io, mode, 2, 3.0);
    FIQuery<double, Line3<double>, Ellipsoid3<double>> query;
    auto r = query(line, ellipsoid);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrLine3Ellipsoid3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto lineOrigin = Pt<3>(io, mode, 3, 4.0);
    auto lineDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    auto ellipsoid = Ellip(io, mode, 2, 3.0);
    Expose<FIQuery<double, Line3<double>, Ellipsoid3<double>>> query;
    FIQuery<double, Line3<double>, Ellipsoid3<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, ellipsoid, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// With a signed-permutation frame, an integer centre and extents that are
// powers of two, M = R*diag(1/e^2)*R^T is exact. Put the line origin at
// C + (k/2)*e_j*U_j + t*U_i with the unit direction U_i: then
// a2 = 1/e_i^2, a1 = t/e_i^2 and a0 = k^2/4 - 1 + t^2/e_i^2, all exact, and
// the discriminant is exactly (1 - k^2/4)/e_i^2. k = +-2 makes it exactly zero
// (the tangent branch), k = 0 gives the integer parameters -t -+ e_i and
// t = 0 puts a root exactly at the line origin.
ORACLE_CASE("IntrLine3Ellipsoid3.find.latticeTangent")
{
    Ellipsoid3<double> ellipsoid{};
    ellipsoid.center = io.latticeVec<3>(-3, 3);
    auto axis = Frame3(io, 0);
    ellipsoid.axis[0] = axis[0];
    ellipsoid.axis[1] = axis[1];
    ellipsoid.axis[2] = axis[2];
    ellipsoid.extent = ExtentPow2<3>(io, 0, 3.0);
    int i = io.rawInteger(0, 2);
    int j = (i + 1 + io.rawInteger(0, 1)) % 3;
    double k = static_cast<double>(io.rawInteger(-2, 2));
    double t = static_cast<double>(io.rawInteger(-3, 3));
    Vector3<double> origin = ellipsoid.center
        + (0.5 * k * ellipsoid.extent[j]) * ellipsoid.axis[j]
        + t * ellipsoid.axis[i];
    io.givenVec(origin);
    io.givenVec(ellipsoid.axis[i]);
    Line3<double> line(origin, ellipsoid.axis[i]);
    FIQuery<double, Line3<double>, Ellipsoid3<double>> query;
    auto r = query(line, ellipsoid);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

// ========================= IntrHalfspace3Ellipsoid3 ======================

ORACLE_CASE("IntrHalfspace3Ellipsoid3.test")
{
    int mode = io.index() % 2;
    Halfspace3<double> halfspace{};
    halfspace.normal = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    halfspace.constant = (mode == 0 ? io.lattice(-4, 4) : io.real(-4.0, 4.0));
    auto ellipsoid = Ellip(io, mode, 2, 3.0);
    TIQuery<double, Halfspace3<double>, Ellipsoid3<double>> query;
    auto r = query(halfspace, ellipsoid);
    io.outBool(r.intersect);
}

// The halfspace plane is tangent to the ellipsoid: with a signed-permutation
// frame and power-of-two extents, MInverse = R*diag(e^2)*R^T is exact, the
// normal is the principal axis U_j, so discr is exactly e_j^2 and the
// projection extent is exactly e_j. Choosing the constant as
// Dot(N,C) + e_j + d makes tmax exactly -d, which is zero for d = 0.
ORACLE_CASE("IntrHalfspace3Ellipsoid3.test.tangent")
{
    Ellipsoid3<double> ellipsoid{};
    ellipsoid.center = io.latticeVec<3>(-3, 3);
    auto axis = Frame3(io, 0);
    ellipsoid.axis[0] = axis[0];
    ellipsoid.axis[1] = axis[1];
    ellipsoid.axis[2] = axis[2];
    ellipsoid.extent = ExtentPow2<3>(io, 0, 3.0);
    int j = io.rawInteger(0, 2);
    double d = static_cast<double>(io.rawInteger(-1, 1));
    Halfspace3<double> halfspace{};
    halfspace.normal = io.givenVec(ellipsoid.axis[j]);
    halfspace.constant = io.given(Dot(halfspace.normal, ellipsoid.center)
        + ellipsoid.extent[j] + d);
    TIQuery<double, Halfspace3<double>, Ellipsoid3<double>> query;
    auto r = query(halfspace, ellipsoid);
    io.outBool(r.intersect);
}

// ======================== IntrEllipsoid3Ellipsoid3 =======================
//
// This is the only query of the group that reaches the C math library on the
// compared path: SymmetricEigensolver3x3 (atan2, cos, sqrt) factors M2 and
// RootsBisection isolates the roots of f(s). Both outputs are discrete, so the
// comparison is still exact; the generators reject the inputs on which
// upstream's ad-hoc bracketing asserts fire (docs/UPSTREAM-FINDINGS.md
// IntrEllipsoid3Ellipsoid3.h GetRoots, issues #255 and #461), which the
// dedicated throw-parity case covers instead.

namespace
{
    struct EllipsoidData
    {
        Vector3<double> center{}, extent{};
        std::array<Vector3<double>, 3> axis{};
    };

    void RawEllipsoidAxes(oracle::Ctx& io, int mode, EllipsoidData& e)
    {
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
                e.axis[i].MakeZero();
                e.axis[i][perm[p][i]] = sign;
            }
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { e.axis[0][i] = io.raw(-1.0, 1.0); }
                len = Length(e.axis[0]);
            } while (len < 0.1 || len > 1.0);
            Normalize(e.axis[0]);
            ComputeOrthogonalComplement(1, e.axis.data());
        }
    }

    // The layout of Ellip(): centre, three axes, extents.
    void EmitEllipsoid(oracle::Ctx& io, EllipsoidData const& e)
    {
        io.givenVec(e.center);
        io.givenVec(e.axis[0]);
        io.givenVec(e.axis[1]);
        io.givenVec(e.axis[2]);
        io.givenVec(e.extent);
    }

    Ellipsoid3<double> MakeEllipsoid(EllipsoidData const& e)
    {
        Ellipsoid3<double> ellipsoid{};
        ellipsoid.center = e.center;
        ellipsoid.axis[0] = e.axis[0];
        ellipsoid.axis[1] = e.axis[1];
        ellipsoid.axis[2] = e.axis[2];
        ellipsoid.extent = e.extent;
        return ellipsoid;
    }

    // True when upstream classifies the pair without throwing.
    bool Ellipsoid3Ellipsoid3Succeeds(EllipsoidData const& a,
        EllipsoidData const& b)
    {
        try
        {
            TIQuery<double, Ellipsoid3<double>, Ellipsoid3<double>> query;
            (void)query(MakeEllipsoid(a), MakeEllipsoid(b));
            return true;
        }
        catch (std::exception const&)
        {
            return false;
        }
    }
}

// The generator rejects the configurations on which upstream's bracketing
// asserts fire; the throw-parity case below covers one of those.
ORACLE_CASE("IntrEllipsoid3Ellipsoid3.test")
{
    int mode = io.index() % 2;
    EllipsoidData a{}, b{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        for (int i = 0; i < 3; ++i)
        {
            a.center[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                : io.raw(-2.0, 2.0));
        }
        RawEllipsoidAxes(io, mode, a);
        for (int i = 0; i < 3; ++i)
        {
            a.extent[i] = (mode == 0 ? std::ldexp(1.0, io.rawInteger(-1, 1))
                : io.raw(0.5, 2.5));
        }
        for (int i = 0; i < 3; ++i)
        {
            b.center[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                : io.raw(-2.0, 2.0));
        }
        RawEllipsoidAxes(io, mode, b);
        for (int i = 0; i < 3; ++i)
        {
            b.extent[i] = (mode == 0 ? std::ldexp(1.0, io.rawInteger(-1, 1))
                : io.raw(0.5, 2.5));
        }
        if (Ellipsoid3Ellipsoid3Succeeds(a, b)) { break; }
    }
    EmitEllipsoid(io, a);
    EmitEllipsoid(io, b);
    TIQuery<double, Ellipsoid3<double>, Ellipsoid3<double>> query;
    auto r = query(MakeEllipsoid(a), MakeEllipsoid(b));
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.classification));
}

// The common-centre branch: K2 and therefore K are exactly zero, so the
// closed-form 'K == 0' path is taken and no bisection is performed. All three
// classifications of that path occur.
ORACLE_CASE("IntrEllipsoid3Ellipsoid3.test.concentric")
{
    int mode = io.index() % 2;
    Vector3<double> center = Pt<3>(io, mode, 2, 2.0);
    EllipsoidData a{}, b{};
    a.center = center;
    b.center = center;
    RawEllipsoidAxes(io, mode, a);
    io.givenVec(a.axis[0]);
    io.givenVec(a.axis[1]);
    io.givenVec(a.axis[2]);
    a.extent = ExtentPow2<3>(io, mode, 2.5);
    RawEllipsoidAxes(io, mode, b);
    io.givenVec(b.axis[0]);
    io.givenVec(b.axis[1]);
    io.givenVec(b.axis[2]);
    b.extent = ExtentPow2<3>(io, mode, 2.5);
    TIQuery<double, Ellipsoid3<double>, Ellipsoid3<double>> query;
    auto r = query(MakeEllipsoid(a), MakeEllipsoid(b));
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.classification));
}

// Both frames are the identity and the extents are powers of two, so M2 is
// exactly diagonal and its eigenvalues are exactly the diagonal entries.
// Repeated extents then make two (or three) of the d-values exactly equal,
// which is what reaches the 'd0 > d1 = d2', 'd0 = d1 > d2' and 'd0 = d1 = d2'
// branches of the valid-pair analysis and with them the one- and two-argument
// GetRoots overloads.
ORACLE_CASE("IntrEllipsoid3Ellipsoid3.test.equalEigenvalues")
{
    EllipsoidData a{}, b{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        for (int i = 0; i < 3; ++i)
        {
            a.center[i] = static_cast<double>(io.rawInteger(-2, 2));
            a.axis[i].MakeZero();
            a.axis[i][i] = 1.0;
            b.axis[i].MakeZero();
            b.axis[i][i] = 1.0;
        }
        double ra = std::ldexp(1.0, io.rawInteger(-1, 1));
        a.extent = Vector3<double>{ ra, ra, ra };
        for (int i = 0; i < 3; ++i)
        {
            b.center[i] = static_cast<double>(io.rawInteger(-2, 2));
        }
        double s0 = std::ldexp(1.0, io.rawInteger(-1, 1));
        double s1 = std::ldexp(1.0, io.rawInteger(-1, 1));
        b.extent[0] = s0;
        b.extent[1] = (io.rawInteger(0, 1) == 0 ? s0 : s1);
        b.extent[2] = (io.rawInteger(0, 1) == 0 ? s0 : s1);
        if (Ellipsoid3Ellipsoid3Succeeds(a, b)) { break; }
    }
    EmitEllipsoid(io, a);
    EmitEllipsoid(io, b);
    TIQuery<double, Ellipsoid3<double>, Ellipsoid3<double>> query;
    auto r = query(MakeEllipsoid(a), MakeEllipsoid(b));
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.classification));
}

// Throw parity for the reachable bracketing assert: a ball of radius 0.4 at
// the origin against an ellipsoid whose centre is 0.001 away and whose frame
// is a near-180-degree rotation about x. Upstream throws "Unexpected
// condition." instead of classifying, and the port preserves that (the same
// configuration is pinned in test/IntrEllipsoid3Ellipsoid3.test.ts).
ORACLE_CASE("IntrEllipsoid3Ellipsoid3.test.bracketAssert")
{
    EllipsoidData a{}, b{};
    for (int i = 0; i < 3; ++i) { a.axis[i].MakeZero(); a.axis[i][i] = 1.0; }
    a.extent = Vector3<double>{ 0.4, 0.4, 0.4 };
    b.center = Vector3<double>{ 0.001, 0.0, 0.0 };
    b.axis[0] = Vector3<double>{ 0.9999995000000417, 0.0, -0.0009999998333333417 };
    b.axis[1] = Vector3<double>{ 8.586571855778066e-11, -0.9999999999999963,
        8.586568993587257e-8 };
    b.axis[2] = Vector3<double>{ -0.000999999833333338, -8.586573286873543e-8,
        -0.999999500000038 };
    b.extent = Vector3<double>{ 0.4904589041303176, 0.41139783664312857,
        0.9127878433790204 };
    EmitEllipsoid(io, a);
    EmitEllipsoid(io, b);
    TIQuery<double, Ellipsoid3<double>, Ellipsoid3<double>> query;
    auto r = query(MakeEllipsoid(a), MakeEllipsoid(b));
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.classification));
}

// ======================== IntrOrientedBox3Cylinder3 ======================
//
// The query is a thin wrapper around IntrCanonicalBox3Cylinder3, so it
// inherits the (U1,-D) sign typo of DoQueryNoZeros and the port's correction
// of it (docs/UPSTREAM-FINDINGS.md IntrCanonicalBox3Cylinder3.h, issue #197;
// first demonstrated by group 31, reused for the aligned-box wrapper by group
// 32). The code below is upstream's DoQueryNoZeros with the single corrected
// sign, used only by the generator to tell the sound inputs from the
// defective ones.

namespace
{
    double CbcSqrDistance(Vector3<double> const& P0, Vector3<double> const& P1,
        Vector3<double> const& C, Vector3<double> const& W0,
        Vector3<double> const& W1)
    {
        Vector3<double> P0mC = P0 - C;
        Vector3<double> P1mC = P1 - C;
        Vector2<double> Q0{ Dot(W0, P0mC), Dot(W1, P0mC) };
        Vector2<double> Q1{ Dot(W0, P1mC), Dot(W1, P1mC) };

        Vector2<double> direction = Q1 - Q0;
        double s = Dot(direction, Q1);
        if (s <= 0.0)
        {
            return Dot(Q1, Q1);
        }
        s = Dot(direction, Q0);
        if (s >= 0.0)
        {
            return Dot(Q0, Q0);
        }
        s /= Dot(direction, direction);
        Vector2<double> closest = Q0 - s * direction;
        return Dot(closest, closest);
    }

    bool CbcFixedNoZeros(Vector3<double> const& C, Vector3<double> const& D,
        double r, double hDiv2, Vector3<double> const& E)
    {
        std::array<double, 3> negEmCDivD
        {
            (-E[0] - C[0]) / D[0], (-E[1] - C[1]) / D[1], (-E[2] - C[2]) / D[2]
        };
        std::array<double, 3> posEmCDivD
        {
            (E[0] - C[0]) / D[0], (E[1] - C[1]) / D[1], (E[2] - C[2]) / D[2]
        };
        double max01 = std::max(negEmCDivD[0], negEmCDivD[1]);
        double max23 = std::max(negEmCDivD[2], -hDiv2);
        double lower = std::max(max01, max23);
        double min01 = std::min(posEmCDivD[0], posEmCDivD[1]);
        double min23 = std::min(posEmCDivD[2], hDiv2);
        double upper = std::min(min01, min23);
        if (lower <= upper)
        {
            return true;
        }

        double dotDC = Dot(D, C);
        double d0e0 = D[0] * E[0], d1e1 = D[1] * E[1], d2e2 = D[2] * E[2];
        double t1 = +d0e0 - d1e1 - d2e2 - dotDC, s1p = t1 + hDiv2, s1n = t1 - hDiv2;
        double t2 = -d0e0 + d1e1 - d2e2 - dotDC, s2p = t2 + hDiv2, s2n = t2 - hDiv2;
        double t3 = +d0e0 + d1e1 - d2e2 - dotDC, s3p = t3 + hDiv2, s3n = t3 - hDiv2;
        double t4 = -d0e0 - d1e1 + d2e2 - dotDC, s4p = t4 + hDiv2, s4n = t4 - hDiv2;
        double t5 = +d0e0 - d1e1 + d2e2 - dotDC, s5p = t5 + hDiv2, s5n = t5 - hDiv2;
        double t6 = -d0e0 + d1e1 + d2e2 - dotDC, s6p = t6 + hDiv2, s6n = t6 - hDiv2;

        std::array<Vector3<double>, 3> basis{};
        basis[0] = D;
        ComputeOrthogonalComplement(1, basis.data());
        auto const& W0 = basis[1];
        auto const& W1 = basis[2];

        double sqrRadius = r * r;
        double sqrDistance{};
        Vector3<double> P0{}, P1{};

        // (U0, -U1)
        lower = (s1p >= 0.0 ? -E[2] : -E[2] - s1p / D[2]);
        upper = (s5n <= 0.0 ? +E[2] : +E[2] - s5n / D[2]);
        if (lower <= upper)
        {
            P0 = { +E[0], -E[1], lower };
            P1 = { +E[0], -E[1], upper };
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (U1, -U0)
        lower = (s2p >= 0.0 ? -E[2] : -E[2] - s2p / D[2]);
        upper = (s6n <= 0.0 ? +E[2] : +E[2] - s6n / D[2]);
        if (lower <= upper)
        {
            P0 = { -E[0], +E[1], lower };
            P1 = { -E[0], +E[1], upper };
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (U0, -U2)
        lower = (s1p >= 0.0 ? -E[1] : -E[1] - s1p / D[1]);
        upper = (s3n <= 0.0 ? +E[1] : +E[1] - s3n / D[1]);
        if (lower <= upper)
        {
            P0 = { +E[0], lower, -E[2] };
            P1 = { +E[0], upper, -E[2] };
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (U2, -U0)
        lower = (s4p >= 0.0 ? -E[1] : -E[1] - s4p / D[1]);
        upper = (s6n <= 0.0 ? +E[1] : +E[1] - s6n / D[1]);
        if (lower <= upper)
        {
            P0 = { -E[0], lower, +E[2] };
            P1 = { -E[0], upper, +E[2] };
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (U1, -U2)
        lower = (s2p >= 0.0 ? -E[0] : -E[0] - s2p / D[0]);
        upper = (s3n <= 0.0 ? +E[0] : +E[0] - s3n / D[0]);
        if (lower <= upper)
        {
            P0 = { lower, +E[1], -E[2] };
            P1 = { upper, +E[1], -E[2] };
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (U2, -U1)
        lower = (s4p >= 0.0 ? -E[0] : -E[0] - s4p / D[0]);
        upper = (s5n <= 0.0 ? +E[0] : +E[0] - s5n / D[0]);
        if (lower <= upper)
        {
            P0 = { lower, -E[1], +E[2] };
            P1 = { upper, -E[1], +E[2] };
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (U0, -D)
        lower = (s3p >= 0.0 ? -E[2] : -E[2] - s3p / D[2]);
        upper = (s5p <= 0.0 ? +E[2] : +E[2] - s5p / D[2]);
        if (lower <= upper)
        {
            P0 = (s3p >= 0.0 ? Vector3<double>{ +E[0], +E[1] - s3p / D[1], -E[2] }
                : Vector3<double>{ +E[0], +E[1], -E[2] - s3p / D[2] });
            P1 = (s5p <= 0.0 ? Vector3<double>{ +E[0], -E[1] - s5p / D[1], +E[2] }
                : Vector3<double>{ +E[0], -E[1], +E[2] - s5p / D[2] });
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (D, -U0)
        lower = (s2n >= 0.0 ? -E[2] : -E[2] - s2n / D[2]);
        upper = (s4n <= 0.0 ? +E[2] : +E[2] - s4n / D[2]);
        if (lower <= upper)
        {
            P0 = (s2n >= 0.0 ? Vector3<double>{ -E[0], +E[1] - s2n / D[1], -E[2] }
                : Vector3<double>{ -E[0], +E[1], -E[2] - s2n / D[2] });
            P1 = (s4n <= 0.0 ? Vector3<double>{ -E[0], -E[1] - s4n / D[1], +E[2] }
                : Vector3<double>{ -E[0], -E[1], +E[2] - s4n / D[2] });
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (U1, -D); the corrected sign of E[1] in the second P1 is the only
        // difference from the header.
        lower = (s6p >= 0.0 ? -E[0] : -E[0] - s6p / D[0]);
        upper = (s3p <= 0.0 ? +E[0] : +E[0] - s3p / D[0]);
        if (lower <= upper)
        {
            P0 = (s6p >= 0.0 ? Vector3<double>{ -E[0], +E[1], +E[2] - s6p / D[2] }
                : Vector3<double>{ -E[0] - s6p / D[0], +E[1], +E[2] });
            P1 = (s3p <= 0.0 ? Vector3<double>{ +E[0], +E[1], -E[2] - s3p / D[2] }
                : Vector3<double>{ +E[0] - s3p / D[0], +E[1], -E[2] });
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (D, -U1)
        lower = (s4n >= 0.0 ? -E[0] : -E[0] - s4n / D[0]);
        upper = (s1n <= 0.0 ? +E[0] : +E[0] - s1n / D[0]);
        if (lower <= upper)
        {
            P0 = (s4n >= 0.0 ? Vector3<double>{ -E[0], -E[1], +E[2] - s4n / D[2] }
                : Vector3<double>{ -E[0] - s4n / D[0], -E[1], +E[2] });
            P1 = (s1n <= 0.0 ? Vector3<double>{ +E[0], -E[1], -E[2] - s1n / D[2] }
                : Vector3<double>{ +E[0] - s1n / D[0], -E[1], -E[2] });
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (U2, -D)
        lower = (s5p >= 0.0 ? -E[1] : -E[1] - s5p / D[1]);
        upper = (s6p <= 0.0 ? +E[1] : +E[1] - s6p / D[1]);
        if (lower <= upper)
        {
            P0 = (s5p >= 0.0 ? Vector3<double>{ +E[0] - s5p / D[0], -E[1], +E[2] }
                : Vector3<double>{ +E[0], -E[1] - s5p / D[1], +E[2] });
            P1 = (s6p <= 0.0 ? Vector3<double>{ -E[0] - s6p / D[0], +E[1], +E[2] }
                : Vector3<double>{ -E[0], E[1] - s6p / D[1], +E[2] });
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        // (D, -U2)
        lower = (s1n >= 0.0 ? -E[1] : -E[1] - s1n / D[1]);
        upper = (s2n <= 0.0 ? +E[1] : +E[1] - s2n / D[1]);
        if (lower <= upper)
        {
            P0 = (s1n >= 0.0 ? Vector3<double>{ +E[0] - s1n / D[0], -E[1], -E[2] }
                : Vector3<double>{ +E[0], -E[1] - s1n / D[1], -E[2] });
            P1 = (s2n <= 0.0 ? Vector3<double>{ -E[0] - s2n / D[0], +E[1], -E[2] }
                : Vector3<double>{ -E[0], E[1] - s2n / D[1], -E[2] });
            sqrDistance = CbcSqrDistance(P0, P1, C, W0, W1);
            if (sqrDistance <= sqrRadius) { return true; }
        }

        return false;
    }

    // The answer upstream would report if the (U1,-D) typo were corrected.
    // Only the all-nonzero-component path differs, so the other paths return
    // upstream's own answer.
    bool CbcFixed(CanonicalBox3<double> const& box,
        Cylinder3<double> const& cylinder, bool upstreamAnswer)
    {
        Vector3<double> const& C0 = cylinder.axis.origin;
        Vector3<double> const& D0 = cylinder.axis.direction;
        Vector3<double> absD{ std::fabs(D0[0]), std::fabs(D0[1]), std::fabs(D0[2]) };
        double hDiv2 = 0.5 * cylinder.height;
        Vector3<double> const& E = box.extent;
        double intervalCenter = -Dot(D0, C0);
        double intervalRadius = Dot(E, absD);
        if (std::fabs(intervalCenter) > intervalRadius + hDiv2)
        {
            return false;
        }

        Vector3<double> C = C0;
        Vector3<double> D = D0;
        for (int32_t i = 0; i < 3; ++i)
        {
            if (D[i] < 0.0)
            {
                C[i] = -C[i];
                D[i] = -D[i];
            }
        }
        if (D[0] > 0.0 && D[1] > 0.0 && D[2] > 0.0)
        {
            return CbcFixedNoZeros(C, D, cylinder.radius, hDiv2, E);
        }
        return upstreamAnswer;
    }

    struct BoxCylinderData
    {
        Vector3<double> boxCenter{}, boxExtent{};
        std::array<Vector3<double>, 3> boxAxis{};
        Vector3<double> origin{}, direction{};
        double radius{}, height{};
    };

    void RawBoxCylinder(oracle::Ctx& io, int mode, BoxCylinderData& d)
    {
        for (int i = 0; i < 3; ++i)
        {
            d.boxCenter[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                : io.raw(-2.0, 2.0));
        }
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
                d.boxAxis[i].MakeZero();
                d.boxAxis[i][perm[p][i]] = sign;
            }
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { d.boxAxis[0][i] = io.raw(-1.0, 1.0); }
                len = Length(d.boxAxis[0]);
            } while (len < 0.1 || len > 1.0);
            Normalize(d.boxAxis[0]);
            ComputeOrthogonalComplement(1, d.boxAxis.data());
        }
        for (int i = 0; i < 3; ++i)
        {
            d.boxExtent[i] = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
                : io.raw(0.25, 3.0));
        }
        for (int i = 0; i < 3; ++i)
        {
            d.origin[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                : io.raw(-2.0, 2.0));
        }
        if (mode == 0)
        {
            int k = io.rawInteger(0, 5);
            d.direction.MakeZero();
            d.direction[k % 3] = (k < 3 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { d.direction[i] = io.raw(-1.0, 1.0); }
                len = Length(d.direction);
            } while (len < 0.1 || len > 1.0);
            Normalize(d.direction);
        }
        d.radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 2))
            : io.raw(0.25, 2.0));
        d.height = (mode == 0 ? static_cast<double>(io.rawInteger(1, 4))
            : io.raw(0.5, 4.0));
    }

    void EmitBoxCylinder(oracle::Ctx& io, BoxCylinderData const& d)
    {
        io.givenVec(d.boxCenter);
        io.givenVec(d.boxAxis[0]);
        io.givenVec(d.boxAxis[1]);
        io.givenVec(d.boxAxis[2]);
        io.givenVec(d.boxExtent);
        io.givenVec(d.origin);
        io.givenVec(d.direction);
        io.given(d.radius);
        io.given(d.height);
    }

    // 'requireDeviation' selects the inputs on which the corrected
    // IntrCanonicalBox3Cylinder3 answer differs from upstream's (the deviation
    // case) or agrees with it (the main case).
    void OrientedBoxCylinderCase(oracle::Ctx& io, int mode, bool requireDeviation)
    {
        BoxCylinderData d{};
        for (int attempt = 0; attempt < 20000; ++attempt)
        {
            RawBoxCylinder(io, mode, d);
            OrientedBox3<double> probeBox(d.boxCenter, d.boxAxis, d.boxExtent);
            Line3<double> probeAxis(d.origin, d.direction);
            Cylinder3<double> probeCylinder(probeAxis, d.radius, d.height);
            TIQuery<double, OrientedBox3<double>, Cylinder3<double>> probeQuery;
            bool upstream = probeQuery(probeBox, probeCylinder).intersect;
            // The canonical-box query is applied to the transformed cylinder
            // the wrapper builds.
            CanonicalBox3<double> cbox(d.boxExtent);
            Vector3<double> diff = d.origin - d.boxCenter;
            Cylinder3<double> transformed{};
            transformed.radius = d.radius;
            transformed.height = d.height;
            for (int32_t i = 0; i < 3; ++i)
            {
                transformed.axis.origin[i] = Dot(d.boxAxis[i], diff);
                transformed.axis.direction[i] = Dot(d.boxAxis[i], d.direction);
            }
            if ((upstream != CbcFixed(cbox, transformed, upstream)) == requireDeviation)
            {
                break;
            }
        }
        EmitBoxCylinder(io, d);
        OrientedBox3<double> box(d.boxCenter, d.boxAxis, d.boxExtent);
        Line3<double> axis(d.origin, d.direction);
        Cylinder3<double> cylinder(axis, d.radius, d.height);
        TIQuery<double, OrientedBox3<double>, Cylinder3<double>> query;
        auto r = query(box, cylinder);
        io.outBool(r.intersect);
    }
}

// The generator rejects the inputs on which the (U1,-D) sign typo of
// IntrCanonicalBox3Cylinder3::DoQueryNoZeros changes the answer, which is
// exactly the set on which the port deliberately deviates.
ORACLE_CASE("IntrOrientedBox3Cylinder3.test")
{
    OrientedBoxCylinderCase(io, io.index() % 2, false);
}

// Deliberate deviation inherited from the delegate. docs/UPSTREAM-FINDINGS.md
// IntrCanonicalBox3Cylinder3.h DoQueryNoZeros, issue #197.
ORACLE_CASE("IntrOrientedBox3Cylinder3.test.edgeTypoDeviation")
{
    OrientedBoxCylinderCase(io, 1, true);
}

// Upstream asserts that the cylinder is finite, so both sides throw. The
// record has no outputs.
ORACLE_CASE("IntrOrientedBox3Cylinder3.test.infinite")
{
    auto box = OBox3(io, 1, 3, 3.0);
    Vector3<double> origin = io.vec<3>(-2.0, 2.0);
    Vector3<double> direction = io.unit<3>();
    double radius = io.real(0.25, 2.0);
    double height = io.given(-1.0);
    Line3<double> axis(origin, direction);
    Cylinder3<double> cylinder(axis, radius, height);
    TIQuery<double, OrientedBox3<double>, Cylinder3<double>> query;
    auto r = query(box, cylinder);
    io.outBool(r.intersect);
}

// ============================ IntrRay3Cylinder3 ==========================
//
// Upstream provides only an FIQuery for this pair. The delegate
// IntrLine3Cylinder3 has no defect the port corrects, so the main cases need
// no rejection; the one deviation is the missing IsFinite guard, which the
// port turns into an assert (docs/UPSTREAM-FINDINGS.md cylinder headers,
// issues #197, #206, #255).

ORACLE_CASE("IntrRay3Cylinder3.find")
{
    int mode = io.index() % 2;
    auto ray = UnitRy<3>(io, mode, 3, 3.0);
    auto cylinder = Cyl(io, mode, 3, 3.0);
    FIQuery<double, Ray3<double>, Cylinder3<double>> query;
    auto r = query(ray, cylinder);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrRay3Cylinder3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto rayOrigin = Pt<3>(io, mode, 3, 3.0);
    auto rayDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    auto cylinder = Cyl(io, mode, 3, 3.0);
    Expose<FIQuery<double, Ray3<double>, Cylinder3<double>>> query;
    FIQuery<double, Ray3<double>, Cylinder3<double>>::Result r{};
    query.DoQuery(rayOrigin, rayDirection, cylinder, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// A cylinder whose axis is a coordinate axis with an even integer height, an
// integer radius from a Pythagorean pair and a ray aimed at a lattice point of
// the wall or of an end disk, so that the wall roots and the |z| <= h/2 clip
// are evaluated at exact equality.
namespace
{
    // A cylinder whose axis is a coordinate axis with radius 5 and an even
    // integer height, together with a point aimed at a lattice point of the
    // wall circle at a dyadic height. Every coordinate is exact, so the wall
    // quadratic and the |z| <= h/2 clip are evaluated exactly, and
    // z = +-h/2 puts the hit exactly on an end disk.
    struct LatticeCylinderAim
    {
        Vector3<double> center{}, w{}, u{}, origin{};
        double height{};
    };

    LatticeCylinderAim DrawLatticeCylinderAim(oracle::Ctx& io)
    {
        LatticeCylinderAim aim{};
        for (int i = 0; i < 3; ++i)
        {
            aim.center[i] = static_cast<double>(io.rawInteger(-3, 3));
        }
        int k = io.rawInteger(0, 2);
        Vector3<double> v{};
        aim.w.MakeZero();
        aim.w[k] = 1.0;
        aim.u.MakeZero();
        aim.u[(k + 1) % 3] = 1.0;
        v.MakeZero();
        v[(k + 2) % 3] = 1.0;
        aim.height = 2.0 * static_cast<double>(io.rawInteger(1, 4));
        Vector2<double> cp = CirclePoint5(io.rawInteger(0, 11));
        double z = 0.5 * static_cast<double>(io.rawInteger(-2, 2)) * aim.height;
        double a = 0.5 * static_cast<double>(io.rawInteger(-8, 8));
        aim.origin = aim.center + (cp[0] - a) * aim.u + cp[1] * v + z * aim.w;
        return aim;
    }
}

ORACLE_CASE("IntrRay3Cylinder3.find.latticeAimed")
{
    auto aim = DrawLatticeCylinderAim(io);
    io.givenVec(aim.center);
    io.givenVec(aim.w);
    double radius = io.given(5.0);
    double height = io.given(aim.height);
    io.givenVec(aim.origin);
    io.givenVec(aim.u);
    Line3<double> axis(aim.center, aim.w);
    Cylinder3<double> cylinder(axis, radius, height);
    Ray3<double> ray(aim.origin, aim.u);
    FIQuery<double, Ray3<double>, Cylinder3<double>> query;
    auto r = query(ray, cylinder);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

// Deliberate deviation: an infinite cylinder (the height = -1 sentinel).
// Upstream feeds the sentinel into the finite-cylinder formulas, where the
// half-height is negative and the slab is empty; the port asserts.
ORACLE_CASE("IntrRay3Cylinder3.find.infiniteDeviation")
{
    auto ray = UnitRy<3>(io, 1, 3, 3.0);
    Vector3<double> origin = io.vec<3>(-3.0, 3.0);
    Vector3<double> direction = io.unit<3>();
    double radius = io.real(0.5, 3.0);
    double height = io.given(-1.0);
    Line3<double> axis(origin, direction);
    Cylinder3<double> cylinder(axis, radius, height);
    FIQuery<double, Ray3<double>, Cylinder3<double>> query;
    auto r = query(ray, cylinder);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// ========================== IntrSegment3Cylinder3 ========================

ORACLE_CASE("IntrSegment3Cylinder3.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<3>(io, mode, 3, 3.0);
    auto cylinder = Cyl(io, mode, 3, 3.0);
    FIQuery<double, Segment3<double>, Cylinder3<double>> query;
    auto r = query(segment, cylinder);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrSegment3Cylinder3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto segOrigin = Pt<3>(io, mode, 3, 3.0);
    auto segDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    double segExtent = (mode == 0 ? io.lattice(1, 4) : io.real(0.25, 4.0));
    auto cylinder = Cyl(io, mode, 3, 3.0);
    Expose<FIQuery<double, Segment3<double>, Cylinder3<double>>> query;
    FIQuery<double, Segment3<double>, Cylinder3<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, cylinder, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

ORACLE_CASE("IntrSegment3Cylinder3.find.latticeAimed")
{
    auto aim = DrawLatticeCylinderAim(io);
    double e = static_cast<double>(io.rawInteger(1, 6));
    io.givenVec(aim.center);
    io.givenVec(aim.w);
    double radius = io.given(5.0);
    double height = io.given(aim.height);
    Vector3<double> p0 = io.givenVec(aim.origin - e * aim.u);
    Vector3<double> p1 = io.givenVec(aim.origin + e * aim.u);
    Line3<double> axis(aim.center, aim.w);
    Cylinder3<double> cylinder(axis, radius, height);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Cylinder3<double>> query;
    auto r = query(segment, cylinder);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

// The documented upstream limitation, preserved by the port: a zero-length
// segment has no centred-form direction, so the line-cylinder quadratic has a
// vanishing leading coefficient and the reported parameters are NaN or
// infinite. Both sides must produce the same bits.
// docs/UPSTREAM-FINDINGS.md IntrSegment3Cylinder3.h GetCenteredForm, #197.
ORACLE_CASE("IntrSegment3Cylinder3.find.degenerateSegment")
{
    int mode = io.index() % 2;
    auto cylinder = Cyl(io, mode, 3, 3.0);
    Vector3<double> p = Pt<3>(io, mode, 3, 3.0);
    Segment3<double> segment(p, p);
    FIQuery<double, Segment3<double>, Cylinder3<double>> query;
    auto r = query(segment, cylinder);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

// Deliberate deviation: an infinite cylinder, as for the ray query.
ORACLE_CASE("IntrSegment3Cylinder3.find.infiniteDeviation")
{
    auto segment = Sg<3>(io, 1, 3, 3.0);
    Vector3<double> origin = io.vec<3>(-3.0, 3.0);
    Vector3<double> direction = io.unit<3>();
    double radius = io.real(0.5, 3.0);
    double height = io.given(-1.0);
    Line3<double> axis(origin, direction);
    Cylinder3<double> cylinder(axis, radius, height);
    FIQuery<double, Segment3<double>, Cylinder3<double>> query;
    auto r = query(segment, cylinder);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// ===================== IntrRay3Capsule3, IntrSegment3Capsule3 ============
//
// The FI queries delegate to IntrLine3Capsule3, whose two result-corrupting
// defects the port fixes (docs/UPSTREAM-FINDINGS.md IntrLine3Capsule3.h,
// issue #461), so they inherit that deviation; the reference interval and the
// probe below are the ones group 31 uses for the delegate itself. The ray TI
// query delegates to DistRaySegment, whose missing 's0 >= 0' clamp the port
// fixes (issue #126): the port's distance is then larger than upstream's, so
// the main case keeps only the inputs with a nonnegative ray parameter, where
// the clamp is a no-op. The segment TI query delegates to
// DistSegmentSegment::operator(), which the port does not change (only
// ComputeRobust deviates), so it needs no rejection.

namespace
{
    struct CapsuleInterval
    {
        bool intersect;
        int32_t numIntersections;
        std::array<double, 2> parameter;
    };

    // An independent reference for the line/capsule intersection interval. The
    // solid capsule is the union of the solid finite cylinder and the two
    // solid end balls, and a line meets a convex set in an interval, so the
    // interval is [min,max] over the candidate roots on the capsule boundary.
    CapsuleInterval ReferenceLineCapsule(Vector3<double> const& lineOrigin,
        Vector3<double> const& lineDirection, Capsule3<double> const& capsule)
    {
        CapsuleInterval out{ false, 0, { 0.0, 0.0 } };

        Vector3<double> segOrigin{}, segDirection{};
        double segExtent{};
        capsule.segment.GetCenteredForm(segOrigin, segDirection, segExtent);
        std::array<Vector3<double>, 3> basis{};
        basis[0] = segDirection;
        ComputeOrthogonalComplement(1, basis.data());
        double rSqr = capsule.radius * capsule.radius;
        Vector3<double> const& W = basis[0];
        Vector3<double> const& U = basis[1];
        Vector3<double> const& V = basis[2];

        Vector3<double> diff = lineOrigin - segOrigin;
        Vector3<double> P{ Dot(U, diff), Dot(V, diff), Dot(W, diff) };
        double dz = Dot(W, lineDirection);
        if (std::fabs(dz) == 1.0)
        {
            double radialSqrDist = rSqr - P[0] * P[0] - P[1] * P[1];
            if (radialSqrDist >= 0.0)
            {
                out.intersect = true;
                out.numIntersections = 2;
                double zOffset = std::sqrt(radialSqrDist) + segExtent;
                if (dz > 0.0)
                {
                    out.parameter[0] = -P[2] - zOffset;
                    out.parameter[1] = -P[2] + zOffset;
                }
                else
                {
                    out.parameter[0] = P[2] - zOffset;
                    out.parameter[1] = P[2] + zOffset;
                }
            }
            return out;
        }

        Vector3<double> D{ Dot(U, lineDirection), Dot(V, lineDirection), dz };
        double a0 = P[0] * P[0] + P[1] * P[1] - rSqr;
        double a1 = P[0] * D[0] + P[1] * D[1];
        double a2 = D[0] * D[0] + D[1] * D[1];
        double discr = a1 * a1 - a0 * a2;
        if (discr < 0.0)
        {
            return out;
        }

        int32_t count = 0;
        double tMin = 0.0, tMax = 0.0;
        auto accept = [&count, &tMin, &tMax](double t)
        {
            if (count == 0) { tMin = t; tMax = t; }
            else
            {
                if (t < tMin) { tMin = t; }
                if (t > tMax) { tMax = t; }
            }
            ++count;
        };
        auto finish = [&out, &count, &tMin, &tMax]()
        {
            if (count > 0)
            {
                out.intersect = true;
                out.numIntersections = (tMin == tMax ? 1 : 2);
                out.parameter[0] = tMin;
                out.parameter[1] = tMax;
            }
        };

        double root{}, tValue{}, zValue{};
        if (discr > 0.0)
        {
            root = std::sqrt(discr);
            tValue = (-a1 - root) / a2;
            zValue = P[2] + tValue * D[2];
            if (std::fabs(zValue) <= segExtent) { accept(tValue); }

            tValue = (-a1 + root) / a2;
            zValue = P[2] + tValue * D[2];
            if (std::fabs(zValue) <= segExtent) { accept(tValue); }

            if (count == 2)
            {
                finish();
                return out;
            }
        }
        else
        {
            tValue = -a1 / a2;
            zValue = P[2] + tValue * D[2];
            if (std::fabs(zValue) <= segExtent)
            {
                accept(tValue);
                finish();
                return out;
            }
        }

        double PZpE = P[2] + segExtent;
        a1 += PZpE * D[2];
        a0 += PZpE * PZpE;
        discr = a1 * a1 - a0;
        if (discr > 0.0)
        {
            root = std::sqrt(discr);
            accept(-a1 - root);
            accept(-a1 + root);
        }
        else if (discr == 0.0)
        {
            accept(-a1);
        }

        a1 -= 2.0 * segExtent * D[2];
        a0 -= 4.0 * segExtent * P[2];
        discr = a1 * a1 - a0;
        if (discr > 0.0)
        {
            root = std::sqrt(discr);
            accept(-a1 - root);
            accept(-a1 + root);
        }
        else if (discr == 0.0)
        {
            accept(-a1);
        }

        finish();
        return out;
    }

    bool Line3Capsule3IsSound(Vector3<double> const& lineOrigin,
        Vector3<double> const& lineDirection, Capsule3<double> const& capsule)
    {
        Expose<FIQuery<double, Line3<double>, Capsule3<double>>> query;
        FIQuery<double, Line3<double>, Capsule3<double>>::Result r{};
        query.DoQuery(lineOrigin, lineDirection, capsule, r);
        auto reference = ReferenceLineCapsule(lineOrigin, lineDirection, capsule);
        return r.intersect == reference.intersect
            && static_cast<int32_t>(r.numIntersections) == reference.numIntersections
            && r.parameter[0] == reference.parameter[0]
            && r.parameter[1] == reference.parameter[1];
    }

    // Unrecorded draws for a unit-direction ray (or line) and a capsule.
    void RawLineCapsule(oracle::Ctx& io, int mode, Vector3<double>& lineOrigin,
        Vector3<double>& lineDirection, Vector3<double>& p0,
        Vector3<double>& p1, double& radius)
    {
        for (int i = 0; i < 3; ++i)
        {
            lineOrigin[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-3.0, 3.0));
        }
        if (mode == 0)
        {
            int k = io.rawInteger(0, 5);
            lineDirection.MakeZero();
            lineDirection[k % 3] = (k < 3 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { lineDirection[i] = io.raw(-1.0, 1.0); }
                len = Length(lineDirection);
            } while (len < 0.1 || len > 1.0);
            Normalize(lineDirection);
        }
        for (int i = 0; i < 3; ++i)
        {
            p0[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-3.0, 3.0));
        }
        for (int i = 0; i < 3; ++i)
        {
            p1[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-3.0, 3.0));
        }
        radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.25, 3.0));
    }

    // A capsule whose segment is parallel to a coordinate axis, in capsule
    // coordinates, with integer half-length e and radius r.
    Capsule3<double> AxisCapsule(oracle::Ctx& io, int k, double e, double r,
        Vector3<double> const& center)
    {
        Vector3<double> ek{};
        ek.MakeZero();
        ek[k] = 1.0;
        Vector3<double> p0 = center - e * ek;
        Vector3<double> p1 = center + e * ek;
        // One io draw per statement: MSVC evaluates constructor arguments
        // right to left.
        io.givenVec(p0);
        io.givenVec(p1);
        double radius = io.given(r);
        Segment3<double> segment(p0, p1);
        return Capsule3<double>(segment, radius);
    }
}

// The generator keeps only the configurations for which upstream's ray/segment
// distance query reports a nonnegative ray parameter, where the port's clamp
// is a no-op. The capsule radius is drawn as a multiple of that distance, so
// both answers occur about equally often.
ORACLE_CASE("IntrRay3Capsule3.test")
{
    int mode = io.index() % 2;
    DCPQuery<double, Ray3<double>, Segment3<double>> rsQuery{};
    Vector3<double> origin{}, direction{}, p0{}, p1{};
    double unusedRadius{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawLineCapsule(io, mode, origin, direction, p0, p1, unusedRadius);
        auto probe = rsQuery(Ray3<double>(origin, direction),
            Segment3<double>(p0, p1));
        if (probe.parameter[0] >= 0.0) { break; }
    }
    io.givenVec(origin);
    io.givenVec(direction);
    io.givenVec(p0);
    io.givenVec(p1);
    Ray3<double> ray(origin, direction);
    Segment3<double> segment(p0, p1);
    double distance = rsQuery(ray, segment).distance;
    double radius = io.given(distance * io.raw(0.5, 1.5));
    Capsule3<double> capsule(segment, radius);
    TIQuery<double, Ray3<double>, Capsule3<double>> query;
    auto r = query(ray, capsule);
    io.outBool(r.intersect);
}

// Deliberate deviation: upstream's ray/segment distance query never clamps the
// ray parameter to s0 >= 0 in its parallel branch and in regions 1 and 5, so
// it reports a distance measured to a point behind the ray origin. The port
// clamps, which can only increase the distance. Setting the capsule radius to
// exactly upstream's distance makes upstream report an intersection and the
// port report none. docs/UPSTREAM-FINDINGS.md DistRaySegment.h, issue #126.
ORACLE_CASE("IntrRay3Capsule3.test.clampDeviation")
{
    DCPQuery<double, Ray3<double>, Segment3<double>> rsQuery{};
    Vector3<double> origin{}, direction{}, p0{}, p1{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        for (int i = 0; i < 3; ++i) { origin[i] = io.raw(-4.0, 4.0); }
        double len = 0.0;
        do
        {
            for (int i = 0; i < 3; ++i) { direction[i] = io.raw(-1.0, 1.0); }
            len = Length(direction);
        } while (len < 0.1 || len > 1.0);
        Normalize(direction);
        // Bias towards a segment behind the ray origin.
        double back = io.raw(0.5, 4.0);
        Vector3<double> base = origin - back * direction;
        for (int i = 0; i < 3; ++i) { p0[i] = base[i] + io.raw(-1.5, 1.5); }
        for (int i = 0; i < 3; ++i) { p1[i] = base[i] + io.raw(-1.5, 1.5); }
        auto probe = rsQuery(Ray3<double>(origin, direction),
            Segment3<double>(p0, p1));
        if (probe.parameter[0] < 0.0 && probe.distance > 0.25) { break; }
    }
    io.givenVec(origin);
    io.givenVec(direction);
    io.givenVec(p0);
    io.givenVec(p1);
    Ray3<double> ray(origin, direction);
    Segment3<double> segment(p0, p1);
    double radius = io.given(rsQuery(ray, segment).distance);
    Capsule3<double> capsule(segment, radius);
    TIQuery<double, Ray3<double>, Capsule3<double>> query;
    auto r = query(ray, capsule);
    io.outBool(r.intersect);
}

// The generator rejects the inputs on which the delegate IntrLine3Capsule3
// disagrees with the reference interval, which is exactly the set on which the
// port deliberately deviates (issue #461), and the degenerate zero-length
// capsule segment whose frame is undefined.
ORACLE_CASE("IntrRay3Capsule3.find")
{
    int mode = io.index() % 2;
    Vector3<double> origin{}, direction{}, p0{}, p1{};
    double radius{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawLineCapsule(io, mode, origin, direction, p0, p1, radius);
        Segment3<double> probeSegment(p0, p1);
        Capsule3<double> probeCapsule(probeSegment, radius);
        if (Line3Capsule3IsSound(origin, direction, probeCapsule)) { break; }
    }
    io.givenVec(origin);
    io.givenVec(direction);
    io.givenVec(p0);
    io.givenVec(p1);
    io.given(radius);
    Segment3<double> segment(p0, p1);
    Capsule3<double> capsule(segment, radius);
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, Capsule3<double>> query;
    auto r = query(ray, capsule);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrRay3Capsule3.doQuery.fi")
{
    int mode = io.index() % 2;
    Vector3<double> origin{}, direction{}, p0{}, p1{};
    double radius{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawLineCapsule(io, mode, origin, direction, p0, p1, radius);
        Segment3<double> probeSegment(p0, p1);
        Capsule3<double> probeCapsule(probeSegment, radius);
        if (Line3Capsule3IsSound(origin, direction, probeCapsule)) { break; }
    }
    io.givenVec(origin);
    io.givenVec(direction);
    io.givenVec(p0);
    io.givenVec(p1);
    io.given(radius);
    Segment3<double> segment(p0, p1);
    Capsule3<double> capsule(segment, radius);
    Expose<FIQuery<double, Ray3<double>, Capsule3<double>>> query;
    FIQuery<double, Ray3<double>, Capsule3<double>>::Result r{};
    query.DoQuery(origin, direction, capsule, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// Deliberate deviation inherited from the delegate: the ray is tangent to the
// bottom end sphere at its pole, so the single accepted root comes from a
// hemisphere and upstream never sets 'intersect'; the ray clip is guarded on
// that flag, so upstream drops the contact. issue #461 item 1.
ORACLE_CASE("IntrRay3Capsule3.doQuery.capTangentDeviation")
{
    Vector3<double> center = io.latticeVec<3>(-3, 3);
    int k = io.rawInteger(0, 2);
    double e = static_cast<double>(io.rawInteger(1, 3));
    double r = static_cast<double>(io.rawInteger(1, 3));
    auto capsule = AxisCapsule(io, k, e, r, center);
    Vector3<double> ek{}, em{};
    ek.MakeZero();
    ek[k] = 1.0;
    em.MakeZero();
    em[(k + 1) % 3] = 1.0;
    Vector3<double> origin = center - (e + r) * ek
        - static_cast<double>(io.rawInteger(1, 3)) * em;
    io.givenVec(origin);
    io.givenVec(em);
    Expose<FIQuery<double, Ray3<double>, Capsule3<double>>> query;
    FIQuery<double, Ray3<double>, Capsule3<double>>::Result result{};
    query.DoQuery(origin, em, capsule, result);
    io.outBool(result.intersect);
    io.outInt(static_cast<int32_t>(result.numIntersections));
    io.outReal(result.parameter[0]);
    io.outReal(result.parameter[1]);
}

// ========================== IntrSegment3Capsule3 =========================

// DistSegmentSegment::operator() is ported unchanged (only its ComputeRobust
// deviates), so this case needs no rejection.
ORACLE_CASE("IntrSegment3Capsule3.test")
{
    int mode = io.index() % 2;
    Vector3<double> q0{}, q1{}, p0{}, p1{};
    double unusedRadius{};
    RawLineCapsule(io, mode, q0, q1, p0, p1, unusedRadius);
    // RawLineCapsule's second vector is a unit direction; turn it into the
    // second endpoint of the query segment.
    Vector3<double> e0 = io.givenVec(q0);
    Vector3<double> e1 = io.givenVec(q0 + (mode == 0 ? 2.0 : 3.0) * q1);
    io.givenVec(p0);
    io.givenVec(p1);
    Segment3<double> segment(e0, e1);
    Segment3<double> capsuleSegment(p0, p1);
    DCPQuery<double, Segment3<double>, Segment3<double>> ssQuery{};
    double distance = ssQuery(segment, capsuleSegment).distance;
    double radius = io.given(distance * io.raw(0.5, 1.5));
    Capsule3<double> capsule(capsuleSegment, radius);
    TIQuery<double, Segment3<double>, Capsule3<double>> query;
    auto r = query(segment, capsule);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3Capsule3.find")
{
    int mode = io.index() % 2;
    Vector3<double> segOrigin{}, segDirection{}, p0{}, p1{};
    double radius{}, segExtent{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawLineCapsule(io, mode, segOrigin, segDirection, p0, p1, radius);
        segExtent = (mode == 0 ? static_cast<double>(io.rawInteger(1, 4))
            : io.raw(0.25, 4.0));
        Segment3<double> probeSegment(p0, p1);
        Capsule3<double> probeCapsule(probeSegment, radius);
        // Probe with the centred form the query itself will compute from the
        // recorded endpoints, not with the draws they were built from.
        Segment3<double> probeQuerySegment(segOrigin - segExtent * segDirection,
            segOrigin + segExtent * segDirection);
        Vector3<double> co{}, cd{};
        double ce{};
        probeQuerySegment.GetCenteredForm(co, cd, ce);
        if (Line3Capsule3IsSound(co, cd, probeCapsule)) { break; }
    }
    Vector3<double> e0 = io.givenVec(segOrigin - segExtent * segDirection);
    Vector3<double> e1 = io.givenVec(segOrigin + segExtent * segDirection);
    io.givenVec(p0);
    io.givenVec(p1);
    io.given(radius);
    Segment3<double> segment(e0, e1);
    Segment3<double> capsuleSegment(p0, p1);
    Capsule3<double> capsule(capsuleSegment, radius);
    FIQuery<double, Segment3<double>, Capsule3<double>> query;
    auto r = query(segment, capsule);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    if (r.intersect)
    {
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrSegment3Capsule3.doQuery.fi")
{
    int mode = io.index() % 2;
    Vector3<double> segOrigin{}, segDirection{}, p0{}, p1{};
    double radius{}, segExtent{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawLineCapsule(io, mode, segOrigin, segDirection, p0, p1, radius);
        segExtent = (mode == 0 ? static_cast<double>(io.rawInteger(1, 4))
            : io.raw(0.25, 4.0));
        Segment3<double> probeSegment(p0, p1);
        Capsule3<double> probeCapsule(probeSegment, radius);
        if (Line3Capsule3IsSound(segOrigin, segDirection, probeCapsule)) { break; }
    }
    io.givenVec(segOrigin);
    io.givenVec(segDirection);
    io.given(segExtent);
    io.givenVec(p0);
    io.givenVec(p1);
    io.given(radius);
    Segment3<double> capsuleSegment(p0, p1);
    Capsule3<double> capsule(capsuleSegment, radius);
    Expose<FIQuery<double, Segment3<double>, Capsule3<double>>> query;
    FIQuery<double, Segment3<double>, Capsule3<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, capsule, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// Deliberate deviation inherited from the delegate: the segment meets the
// capsule wall exactly on the cap-junction circle z = -e. Upstream's wall test
// and bottom-hemisphere test both accept that one root, the count reaches two
// and the early return collapses the intersection to a single point,
// discarding the true far endpoint. issue #461 item 2.
ORACLE_CASE("IntrSegment3Capsule3.doQuery.junctionDeviation")
{
    Vector3<double> center = io.latticeVec<3>(-3, 3);
    double e = static_cast<double>(io.rawInteger(1, 3));
    double r = static_cast<double>(io.rawInteger(1, 3));
    auto capsule = AxisCapsule(io, 2, e, r, center);
    // Capsule coordinates are U = (0,1,0), V = (-1,0,0), W = (0,0,1), so the
    // world point center + (0,r,-e) has capsule coordinates (r,0,-e) and the
    // world direction (0,1/2,-sqrt(3)/2) has capsule coordinates
    // (1/2,0,-sqrt(3)/2). The wall root is then exactly 0 and its z-value is
    // exactly -e.
    Vector3<double> segOrigin{ center[0], center[1] + r, center[2] - e };
    Vector3<double> segDirection{ 0.0, 0.5, -std::sqrt(0.75) };
    io.givenVec(segOrigin);
    io.givenVec(segDirection);
    double segExtent = io.given(static_cast<double>(io.rawInteger(4, 8)));
    Expose<FIQuery<double, Segment3<double>, Capsule3<double>>> query;
    FIQuery<double, Segment3<double>, Capsule3<double>>::Result result{};
    query.DoQuery(segOrigin, segDirection, segExtent, capsule, result);
    io.outBool(result.intersect);
    io.outInt(static_cast<int32_t>(result.numIntersections));
    io.outReal(result.parameter[0]);
    io.outReal(result.parameter[1]);
}

// ========================== IntrLine2SegmentMesh2 ========================
//
// Upstream sorts the intersection records with std::sort, which is not stable,
// so records with equal line parameters (two mesh segments meeting the line at
// a shared vertex, or the two endpoints of a coincident segment) have an
// unspecified relative order. Both sides therefore re-sort the records by a
// total order on (lineParameter, index pair, mesh parameter, point) before
// emitting them; records with equal keys are equal in every field.

namespace
{
    struct MeshHit
    {
        double lineParameter, meshParameter, px, py;
        int32_t i0, i1;
    };

    bool MeshHitLess(MeshHit const& a, MeshHit const& b)
    {
        if (a.lineParameter != b.lineParameter) { return a.lineParameter < b.lineParameter; }
        if (a.i0 != b.i0) { return a.i0 < b.i0; }
        if (a.i1 != b.i1) { return a.i1 < b.i1; }
        if (a.meshParameter != b.meshParameter) { return a.meshParameter < b.meshParameter; }
        if (a.px != b.px) { return a.px < b.px; }
        return a.py < b.py;
    }

    void OutMeshResult(oracle::Ctx& io,
        FIQuery<double, Line2<double>, SegmentMesh2<double>>::Result const& r)
    {
        std::vector<MeshHit> hits;
        hits.reserve(r.intersections.size());
        for (auto const& x : r.intersections)
        {
            hits.push_back(MeshHit{ x.lineParameter, x.meshSegmentParameter,
                x.point[0], x.point[1], static_cast<int32_t>(x.indexPair[0]),
                static_cast<int32_t>(x.indexPair[1]) });
        }
        std::sort(hits.begin(), hits.end(), MeshHitLess);
        io.outInt(static_cast<int32_t>(hits.size()));
        for (auto const& h : hits)
        {
            io.outInt(h.i0);
            io.outInt(h.i1);
            io.outReal(h.lineParameter);
            io.outReal(h.meshParameter);
            io.outReal(h.px);
            io.outReal(h.py);
        }
    }

    // Six vertices, four index pairs and a topology selector are recorded on
    // every record, whichever constructor is used.
    SegmentMesh2<double> DrawMesh(oracle::Ctx& io, int mode, int topology,
        std::vector<Vector2<double>> const& vertices,
        std::vector<std::array<size_t, 2>> const& indices)
    {
        (void)mode;
        if (topology == 0) { return SegmentMesh2<double>(vertices); }
        if (topology == 1) { return SegmentMesh2<double>(vertices, true); }
        if (topology == 2) { return SegmentMesh2<double>(vertices, false); }
        return SegmentMesh2<double>(vertices, indices, true);
    }
}

ORACLE_CASE("IntrLine2SegmentMesh2.find")
{
    int mode = io.index() % 2;
    int topology = io.integer(0, 3);
    std::vector<Vector2<double>> vertices(6);
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i] = Pt<2>(io, mode, 4, 5.0);
    }
    std::vector<std::array<size_t, 2>> indices(4);
    for (size_t i = 0; i < indices.size(); ++i)
    {
        int a = io.integer(0, 5);
        int b = io.integer(0, 5);
        indices[i] = { static_cast<size_t>(a), static_cast<size_t>(b) };
    }
    auto line = UnitLn<2>(io, mode, 4, 5.0);
    auto mesh = DrawMesh(io, mode, topology, vertices, indices);
    FIQuery<double, Line2<double>, SegmentMesh2<double>> query;
    auto r = query(line, mesh);
    OutMeshResult(io, r);
}

// Four of the six vertices lie on one lattice line, which is also the query
// line, so the collinear branch of IntrLine2Segment2 is reached exactly and
// upstream's +-max() line-parameter sentinel is recorded. Two mesh segments
// also share a vertex on the line, which is where the unstable sort produces
// records with equal line parameters.
ORACLE_CASE("IntrLine2SegmentMesh2.find.collinear")
{
    int topology = io.integer(0, 3);
    Vector2<double> base{
        static_cast<double>(io.rawInteger(-3, 3)),
        static_cast<double>(io.rawInteger(-3, 3)) };
    Vector2<double> step{};
    do
    {
        step[0] = static_cast<double>(io.rawInteger(-2, 2));
        step[1] = static_cast<double>(io.rawInteger(-2, 2));
    } while (step[0] == 0.0 && step[1] == 0.0);
    std::vector<Vector2<double>> vertices(6);
    for (size_t i = 0; i < 4; ++i)
    {
        double k = static_cast<double>(io.rawInteger(-3, 3));
        vertices[i] = io.givenVec(base + k * step);
    }
    vertices[4] = io.latticeVec<2>(-4, 4);
    vertices[5] = io.latticeVec<2>(-4, 4);
    std::vector<std::array<size_t, 2>> indices(4);
    for (size_t i = 0; i < indices.size(); ++i)
    {
        int a = io.integer(0, 5);
        int b = io.integer(0, 5);
        indices[i] = { static_cast<size_t>(a), static_cast<size_t>(b) };
    }
    Vector2<double> origin = io.givenVec(base);
    Vector2<double> direction = io.givenVec(step);
    Line2<double> line(origin, direction);
    auto mesh = DrawMesh(io, 0, topology, vertices, indices);
    FIQuery<double, Line2<double>, SegmentMesh2<double>> query;
    auto r = query(line, mesh);
    OutMeshResult(io, r);
}

// =========================== IntrSphere3Frustum3 =========================
//
// The query is 'distance from the sphere centre to the frustum <= radius', so
// it inherits the deviation of DistPoint3Frustum3, two of whose ten far-edge
// assignments do not clamp the free coordinate: upstream's closest point can
// lie off the end of the edge, outside the frustum, and its distance is then
// smaller than the port's. The predicate below is the one group 20 uses for
// that query; the main case rejects it (where the port's clamp is a no-op)
// and the deviation case keeps only it.
// docs/UPSTREAM-FINDINGS.md DistPoint3Frustum3.h, issue #421.

namespace
{
    struct FrustumData
    {
        Vector3<double> origin{}, rVector{}, uVector{}, dVector{};
        double dMin{}, dMax{}, uBound{}, rBound{};
        Vector3<double> point{};
    };

    bool FrustumClampChanges(FrustumData const& f)
    {
        Vector3<double> diff = f.point - f.origin;
        Vector3<double> test{ Dot(diff, f.rVector), Dot(diff, f.uVector),
            Dot(diff, f.dVector) };
        if (test[0] < 0.0) { test[0] = -test[0]; }
        if (test[1] < 0.0) { test[1] = -test[1]; }

        double dRatio = f.dMax / f.dMin;
        double rmin = f.rBound;
        double rmax = dRatio * rmin;
        double umin = f.uBound;
        double umax = dRatio * umin;
        double dmin = f.dMin;
        double dmax = f.dMax;
        double rminSqr = rmin * rmin;
        double uminSqr = umin * umin;
        double dminSqr = dmin * dmin;
        double minRDDot = rminSqr + dminSqr;
        double minUDDot = uminSqr + dminSqr;
        double minRUDDot = rminSqr + minUDDot;
        double maxRDDot = dRatio * minRDDot;
        double maxUDDot = dRatio * minUDDot;

        if (test[2] >= dmax)
        {
            return false;
        }

        if (test[2] <= dmin)
        {
            if (test[0] <= rmin)
            {
                return false;
            }
            if (test[1] <= umin)
            {
                return false;
            }
            double rudDot = rmin * test[0] + umin * test[1] + dmin * test[2];
            double rEdgeDot = umin * rudDot - minRUDDot * test[1];
            if (rEdgeDot >= 0.0)
            {
                double rdDot = rmin * test[0] + dmin * test[2];
                return rdDot >= maxRDDot && test[1] > umax;
            }
            double uEdgeDot = rmin * rudDot - minRUDDot * test[0];
            if (uEdgeDot >= 0.0)
            {
                double udDot = umin * test[1] + dmin * test[2];
                return udDot >= maxUDDot && test[0] > rmax;
            }
            return false;
        }

        double rDot = dmin * test[0] - rmin * test[2];
        double uDot = dmin * test[1] - umin * test[2];
        if (rDot <= 0.0 || uDot <= 0.0)
        {
            return false;
        }
        double rudDot = rmin * test[0] + umin * test[1] + dmin * test[2];
        double rEdgeDot = umin * rudDot - minRUDDot * test[1];
        if (rEdgeDot >= 0.0)
        {
            double rdDot = rmin * test[0] + dmin * test[2];
            return rdDot >= maxRDDot && test[1] > umax;
        }
        double uEdgeDot = rmin * rudDot - minRUDDot * test[0];
        if (uEdgeDot >= 0.0)
        {
            double udDot = umin * test[1] + dmin * test[2];
            return udDot >= maxUDDot && test[0] > rmax;
        }
        return false;
    }

    // Unrecorded draws. The test point is built in frustum coordinates so that
    // every Voronoi region is reachable.
    FrustumData RawFrustum(oracle::Ctx& io)
    {
        FrustumData f{};
        for (int i = 0; i < 3; ++i) { f.origin[i] = io.raw(-4.0, 4.0); }
        double len = 0.0;
        do
        {
            for (int i = 0; i < 3; ++i) { f.rVector[i] = io.raw(-1.0, 1.0); }
            len = Length(f.rVector);
        } while (len < 0.1 || len > 1.0);
        Normalize(f.rVector);
        Vector3<double> u{};
        do
        {
            for (int i = 0; i < 3; ++i) { u[i] = io.raw(-1.0, 1.0); }
            u = u - Dot(u, f.rVector) * f.rVector;
            len = Length(u);
        } while (len < 0.25);
        Normalize(u);
        f.uVector = u;
        f.dVector = Cross(f.rVector, f.uVector);
        f.dMin = io.raw(0.5, 2.0);
        f.dMax = f.dMin + io.raw(0.5, 4.0);
        f.uBound = io.raw(0.25, 2.0);
        f.rBound = io.raw(0.25, 2.0);
        double a = io.raw(-8.0, 8.0);
        double b = io.raw(-8.0, 8.0);
        double c = io.raw(-2.0, 8.0);
        f.point = f.origin + a * f.rVector + b * f.uVector + c * f.dVector;
        return f;
    }

    void EmitFrustum(oracle::Ctx& io, FrustumData const& f)
    {
        io.givenVec(f.origin);
        io.givenVec(f.dVector);
        io.givenVec(f.uVector);
        io.givenVec(f.rVector);
        io.given(f.dMin);
        io.given(f.dMax);
        io.given(f.uBound);
        io.given(f.rBound);
        io.givenVec(f.point);
    }
}

ORACLE_CASE("IntrSphere3Frustum3.test")
{
    FrustumData f{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        f = RawFrustum(io);
        if (!FrustumClampChanges(f)) { break; }
    }
    EmitFrustum(io, f);
    Frustum3<double> frustum(f.origin, f.dVector, f.uVector, f.rVector,
        f.dMin, f.dMax, f.uBound, f.rBound);
    DCPQuery<double, Vector3<double>, Frustum3<double>> pfQuery{};
    double distance = pfQuery(f.point, frustum).distance;
    Sphere3<double> sphere{};
    sphere.center = f.point;
    double factor = io.raw(0.5, 1.5);
    double slack = io.raw(0.0, 0.25);
    sphere.radius = io.given(distance * factor + slack);
    TIQuery<double, Sphere3<double>, Frustum3<double>> query;
    auto r = query(sphere, frustum);
    io.outBool(r.intersect);
}

// Deliberate deviation: the radius is exactly upstream's distance, so upstream
// reports an intersection; the port's clamped closest point is farther away,
// so the port reports none.
ORACLE_CASE("IntrSphere3Frustum3.test.clampDeviation")
{
    FrustumData f{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        f = RawFrustum(io);
        if (FrustumClampChanges(f)) { break; }
    }
    EmitFrustum(io, f);
    Frustum3<double> frustum(f.origin, f.dVector, f.uVector, f.rVector,
        f.dMin, f.dMax, f.uBound, f.rBound);
    DCPQuery<double, Vector3<double>, Frustum3<double>> pfQuery{};
    Sphere3<double> sphere{};
    sphere.center = f.point;
    sphere.radius = io.given(pfQuery(f.point, frustum).distance);
    TIQuery<double, Sphere3<double>, Frustum3<double>> query;
    auto r = query(sphere, frustum);
    io.outBool(r.intersect);
}
