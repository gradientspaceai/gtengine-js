// Verify group 31 (intersection): differential cases for the Intr* headers
// listed in plan/verify-groups.json group 31.
//
// Every generator mixes an exact small-lattice mode with a uniform mode (and,
// where a branch needs it, a constructed mode) so that the touching, parallel,
// collinear, tangent and degenerate branches are reached. Almost all of these
// queries use only + - * / sqrt fabs min max and comparisons, so almost every
// case is declared exact on the TypeScript side; the sin/cos used to build
// orthonormal frames, sector angles and unit directions are applied to
// *unrecorded* draws and only the resulting frame, cosine/sine pair or unit
// vector is recorded as an input, so libm never enters the compared
// computation. The two exceptions are marked in the case comments:
// IntrLine3Torus3 (quartic root finder and atan2) and IntrCylinder3Cylinder3
// (cos/sin of the sampled hemisphere directions).
#define ORACLE_FAMILY "v31-intersection"
#include "Oracle.h"

#include <Mathematics/IntrArc2Arc2.h>
#include <Mathematics/IntrCanonicalBox3Cylinder3.h>
#include <Mathematics/IntrCapsule3Capsule3.h>
#include <Mathematics/IntrCircle2Arc2.h>
#include <Mathematics/IntrCylinder3Cylinder3.h>
#include <Mathematics/IntrDisk2Sector2.h>
#include <Mathematics/IntrHalfspace3Capsule3.h>
#include <Mathematics/IntrHalfspace3Cylinder3.h>
#include <Mathematics/IntrLine2Circle2.h>
#include <Mathematics/IntrLine2OrientedBox2.h>
#include <Mathematics/IntrLine2Ray2.h>
#include <Mathematics/IntrLine2Segment2.h>
#include <Mathematics/IntrLine3Capsule3.h>
#include <Mathematics/IntrLine3Cylinder3.h>
#include <Mathematics/IntrLine3OrientedBox3.h>
#include <Mathematics/IntrLine3Torus3.h>
#include <Mathematics/IntrOrientedBox2Sector2.h>
#include <Mathematics/IntrOrientedBox3Frustum3.h>
#include <Mathematics/IntrRay2AlignedBox2.h>
#include <Mathematics/IntrRay2Ray2.h>
#include <Mathematics/IntrRay2Segment2.h>
#include <Mathematics/IntrRay2Triangle2.h>
#include <Mathematics/IntrRay3AlignedBox3.h>
#include <Mathematics/IntrRay3Rectangle3.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
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

    // A rectangle in 3D. The generator records a full orthonormal frame; the
    // rectangle uses only the first two axes.
    Rectangle3<double> Rect(oracle::Ctx& io, int mode, int lat, double range)
    {
        Rectangle3<double> rectangle{};
        rectangle.center = Pt<3>(io, mode, lat, range);
        auto axis = Frame3(io, mode);
        rectangle.axis[0] = axis[0];
        rectangle.axis[1] = axis[1];
        rectangle.extent = Extent<2>(io, mode, lat, range);
        return rectangle;
    }

    // A triangle; lattice mode makes vertices land exactly on lines.
    template <int N>
    Triangle<N, double> Tri(oracle::Ctx& io, int mode, int lat, double range)
    {
        Triangle<N, double> triangle{};
        triangle.v[0] = Pt<N>(io, mode, lat, range);
        triangle.v[1] = Pt<N>(io, mode, lat, range);
        triangle.v[2] = Pt<N>(io, mode, lat, range);
        return triangle;
    }

    // A line whose direction is not required to be unit length. Mode 0 uses a
    // small lattice for origin and direction so parallel and collinear
    // configurations occur exactly.
    template <int N>
    Line<N, double> Ln(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> origin = Pt<N>(io, mode, lat, range);
        Vector<N, double> direction = (mode == 0 ? io.latticeDir<N>(-2, 2) : io.unit<N>());
        return Line<N, double>(origin, direction);
    }

    // A line with a unit-length direction, for the queries that document that
    // precondition (circle, sphere, capsule, cylinder).
    template <int N>
    Line<N, double> UnitLn(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> origin = Pt<N>(io, mode, lat, range);
        Vector<N, double> direction = (mode == 0 ? AxisUnit<N>(io) : io.unit<N>());
        return Line<N, double>(origin, direction);
    }

    template <int N>
    Ray<N, double> Ry(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> origin = Pt<N>(io, mode, lat, range);
        Vector<N, double> direction = (mode == 0 ? io.latticeDir<N>(-2, 2) : io.unit<N>());
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

    Capsule3<double> Caps(oracle::Ctx& io, int mode, int lat, double range)
    {
        Segment3<double> segment = Sg<3>(io, mode, lat, range);
        double radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        return Capsule3<double>(segment, radius);
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

    // A sector whose half-angle is in (0, pi/2], the convexity precondition of
    // both sector queries. The cosine and sine are recorded, so libm does not
    // enter the compared computation.
    Sector2<double> Sect(oracle::Ctx& io, int mode, int lat, double range)
    {
        Sector2<double> sector{};
        sector.vertex = Pt<2>(io, mode, lat, range);
        sector.radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.5, range));
        sector.direction = (mode == 0 ? AxisUnit<2>(io) : io.unit<2>());
        double c{}, s{};
        if (mode == 0)
        {
            int k = io.rawInteger(0, 2);
            c = (k == 0 ? 0.0 : (k == 1 ? 0.6 : 0.8));
            s = (k == 0 ? 1.0 : (k == 1 ? 0.8 : 0.6));
        }
        else
        {
            double angle = io.raw(0.05, 1.5707963267948966);
            c = std::cos(angle);
            s = std::sin(angle);
        }
        sector.cosAngle = io.given(c);
        sector.sinAngle = io.given(s);
        return sector;
    }

    // Exposes the protected DoQuery helpers, which the port exports as public
    // free functions.
    template <typename Q>
    struct Expose : public Q
    {
        using Q::DoQuery;
    };
}

// ============================== IntrLine2Circle2 =========================

ORACLE_CASE("IntrLine2Circle2.test")
{
    int mode = io.index() % 2;
    auto line = UnitLn<2>(io, mode, 3, 4.0);
    auto circle = Circ(io, mode, 3, 4.0);
    TIQuery<double, Line2<double>, Circle2<double>> query;
    auto r = query(line, circle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine2Circle2.find")
{
    int mode = io.index() % 2;
    auto line = UnitLn<2>(io, mode, 3, 4.0);
    auto circle = Circ(io, mode, 3, 4.0);
    FIQuery<double, Line2<double>, Circle2<double>> query;
    auto r = query(line, circle);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

ORACLE_CASE("IntrLine2Circle2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto lineOrigin = Pt<2>(io, mode, 3, 4.0);
    auto lineDirection = (mode == 0 ? AxisUnit<2>(io) : io.unit<2>());
    auto circle = Circ(io, mode, 3, 4.0);
    Expose<FIQuery<double, Line2<double>, Circle2<double>>> query;
    FIQuery<double, Line2<double>, Circle2<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, circle, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// The tangent branch (discr exactly zero) is unreachable with uniform inputs.
// With circle center C and radius R on the lattice, line origin
// P = C + (R+d)*e_i + s*e_j and unit direction e_j, the discriminant is
// exactly R^2 - (R+d)^2: zero for d = 0, positive for d = -1, negative for
// d = +1.
ORACLE_CASE("IntrLine2Circle2.find.tangent")
{
    Circle2<double> circle{};
    circle.center = io.latticeVec<2>(-3, 3);
    circle.radius = io.lattice(1, 3);
    int i = io.rawInteger(0, 1);
    double d = static_cast<double>(io.rawInteger(-1, 1));
    double s = static_cast<double>(io.rawInteger(-3, 3));
    Vector2<double> ei{}, ej{};
    ei.MakeUnit(i);
    ej.MakeUnit(1 - i);
    Vector2<double> origin = circle.center + (circle.radius + d) * ei + s * ej;
    io.givenVec(origin);
    io.givenVec(ej);
    Line2<double> line(origin, ej);
    FIQuery<double, Line2<double>, Circle2<double>> query;
    auto r = query(line, circle);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t k = 0; k < r.numIntersections; ++k)
    {
        io.outReal(r.parameter[k]);
        io.outVec(r.point[k]);
    }
}

// ============================ IntrLine2OrientedBox2 ======================

ORACLE_CASE("IntrLine2OrientedBox2.test")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 4.0);
    auto box = OBox<2>(io, mode, 3, 3.0);
    TIQuery<double, Line2<double>, OrientedBox2<double>> query;
    auto r = query(line, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine2OrientedBox2.find")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 4.0);
    auto box = OBox<2>(io, mode, 3, 3.0);
    FIQuery<double, Line2<double>, OrientedBox2<double>> query;
    auto r = query(line, box);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

// ============================ IntrLine3OrientedBox3 ======================

ORACLE_CASE("IntrLine3OrientedBox3.test")
{
    int mode = io.index() % 2;
    auto line = Ln<3>(io, mode, 3, 4.0);
    auto box = OBox<3>(io, mode, 3, 3.0);
    TIQuery<double, Line3<double>, OrientedBox3<double>> query;
    auto r = query(line, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3OrientedBox3.find")
{
    int mode = io.index() % 2;
    auto line = Ln<3>(io, mode, 3, 4.0);
    auto box = OBox<3>(io, mode, 3, 3.0);
    FIQuery<double, Line3<double>, OrientedBox3<double>> query;
    auto r = query(line, box);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    if (r.intersect)
    {
        // The query fills both parameters and both points whenever it reports
        // an intersection (parameter[1] = parameter[0] for a single point).
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

// ============================== IntrLine2Ray2 ============================

ORACLE_CASE("IntrLine2Ray2.test")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 4.0);
    auto ray = Ry<2>(io, mode, 3, 4.0);
    TIQuery<double, Line2<double>, Ray2<double>> query;
    auto r = query(line, ray);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

ORACLE_CASE("IntrLine2Ray2.find")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 4.0);
    auto ray = Ry<2>(io, mode, 3, 4.0);
    FIQuery<double, Line2<double>, Ray2<double>> query;
    auto r = query(line, ray);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.lineParameter[0]);
    io.outReal(r.lineParameter[1]);
    io.outReal(r.rayParameter[0]);
    io.outReal(r.rayParameter[1]);
    io.outVec(r.point);
}

// ============================= IntrLine2Segment2 =========================

ORACLE_CASE("IntrLine2Segment2.test")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 4.0);
    auto segment = Sg<2>(io, mode, 3, 4.0);
    TIQuery<double, Line2<double>, Segment2<double>> query;
    auto r = query(line, segment);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

ORACLE_CASE("IntrLine2Segment2.find")
{
    int mode = io.index() % 2;
    auto line = Ln<2>(io, mode, 3, 4.0);
    auto segment = Sg<2>(io, mode, 3, 4.0);
    FIQuery<double, Line2<double>, Segment2<double>> query;
    auto r = query(line, segment);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.lineParameter[0]);
    io.outReal(r.lineParameter[1]);
    io.outReal(r.segmentParameter[0]);
    io.outReal(r.segmentParameter[1]);
    io.outVec(r.point);
}

// =============================== IntrRay2Ray2 ============================

ORACLE_CASE("IntrRay2Ray2.test")
{
    int mode = io.index() % 2;
    auto ray0 = Ry<2>(io, mode, 3, 4.0);
    auto ray1 = Ry<2>(io, mode, 3, 4.0);
    TIQuery<double, Ray2<double>, Ray2<double>> query;
    auto r = query(ray0, ray1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

ORACLE_CASE("IntrRay2Ray2.find")
{
    int mode = io.index() % 2;
    auto ray0 = Ry<2>(io, mode, 3, 4.0);
    auto ray1 = Ry<2>(io, mode, 3, 4.0);
    FIQuery<double, Ray2<double>, Ray2<double>> query;
    auto r = query(ray0, ray1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.ray0Parameter[0]);
    io.outReal(r.ray0Parameter[1]);
    io.outReal(r.ray1Parameter[0]);
    io.outReal(r.ray1Parameter[1]);
    io.outVec(r.point[0]);
    io.outVec(r.point[1]);
}

// Collinear rays: both rays lie on one lattice line, so the "same line" branch
// of IntrLine2Line2 is reached exactly and both the same-direction
// (numIntersections = maxInt) and opposite-direction (0, 1 or 2) branches are
// exercised, including the t == 0 case where the origins coincide.
ORACLE_CASE("IntrRay2Ray2.find.collinear")
{
    Vector2<double> base = io.latticeVec<2>(-3, 3);
    Vector2<double> direction = io.latticeDir<2>(-2, 2);
    double t0 = io.lattice(-3, 3);
    double t1 = io.lattice(-3, 3);
    bool flip = io.boolean();
    Vector2<double> origin0 = base + t0 * direction;
    Vector2<double> origin1 = base + t1 * direction;
    Vector2<double> direction1 = (flip ? -direction : direction);
    Ray2<double> ray0(origin0, direction);
    Ray2<double> ray1(origin1, direction1);
    FIQuery<double, Ray2<double>, Ray2<double>> query;
    auto r = query(ray0, ray1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.ray0Parameter[0]);
    io.outReal(r.ray0Parameter[1]);
    io.outReal(r.ray1Parameter[0]);
    io.outReal(r.ray1Parameter[1]);
    io.outVec(r.point[0]);
    io.outVec(r.point[1]);
}

// The collinear branch of the test query has three sub-branches (the second
// ray's origin ahead of, behind, or exactly at the first ray's origin) that
// uniform inputs reach only a handful of times in 2000 records.
ORACLE_CASE("IntrRay2Ray2.test.collinear")
{
    Vector2<double> base = io.latticeVec<2>(-3, 3);
    Vector2<double> direction = io.latticeDir<2>(-2, 2);
    double t0 = io.lattice(-3, 3);
    double t1 = io.lattice(-3, 3);
    bool flip = io.boolean();
    Vector2<double> origin0 = base + t0 * direction;
    Vector2<double> origin1 = base + t1 * direction;
    Vector2<double> direction1 = (flip ? -direction : direction);
    Ray2<double> ray0(origin0, direction);
    Ray2<double> ray1(origin1, direction1);
    TIQuery<double, Ray2<double>, Ray2<double>> query;
    auto r = query(ray0, ray1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

// ============================= IntrRay2Segment2 ==========================

ORACLE_CASE("IntrRay2Segment2.test")
{
    int mode = io.index() % 2;
    auto ray = Ry<2>(io, mode, 3, 4.0);
    auto segment = Sg<2>(io, mode, 3, 4.0);
    TIQuery<double, Ray2<double>, Segment2<double>> query;
    auto r = query(ray, segment);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

ORACLE_CASE("IntrRay2Segment2.find")
{
    int mode = io.index() % 2;
    auto ray = Ry<2>(io, mode, 3, 4.0);
    auto segment = Sg<2>(io, mode, 3, 4.0);
    FIQuery<double, Ray2<double>, Segment2<double>> query;
    auto r = query(ray, segment);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.rayParameter[0]);
    io.outReal(r.rayParameter[1]);
    io.outReal(r.segmentParameter[0]);
    io.outReal(r.segmentParameter[1]);
    io.outVec(r.point[0]);
    io.outVec(r.point[1]);
}

// A ray and a segment on the same lattice line, which reaches the collinear
// branch (and its interval-interval clip) exactly.
ORACLE_CASE("IntrRay2Segment2.find.collinear")
{
    Vector2<double> base = io.latticeVec<2>(-3, 3);
    Vector2<double> direction = io.latticeDir<2>(-2, 2);
    double t0 = io.lattice(-3, 3);
    double t1 = io.lattice(-3, 3);
    double t2 = io.lattice(-3, 3);
    Vector2<double> origin = base + t0 * direction;
    Vector2<double> p0 = base + t1 * direction;
    Vector2<double> p1 = base + t2 * direction;
    Ray2<double> ray(origin, direction);
    Segment2<double> segment(p0, p1);
    FIQuery<double, Ray2<double>, Segment2<double>> query;
    auto r = query(ray, segment);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.rayParameter[0]);
    io.outReal(r.rayParameter[1]);
    io.outReal(r.segmentParameter[0]);
    io.outReal(r.segmentParameter[1]);
    io.outVec(r.point[0]);
    io.outVec(r.point[1]);
}

// The collinear branch of the test query, whose t == 0 sub-branch (the
// right-most segment point exactly at the ray origin) uniform inputs never
// reach.
ORACLE_CASE("IntrRay2Segment2.test.collinear")
{
    Vector2<double> base = io.latticeVec<2>(-3, 3);
    Vector2<double> direction = io.latticeDir<2>(-2, 2);
    double t0 = io.lattice(-3, 3);
    double t1 = io.lattice(-3, 3);
    double t2 = io.lattice(-3, 3);
    Vector2<double> origin = base + t0 * direction;
    Vector2<double> p0 = base + t1 * direction;
    Vector2<double> p1 = base + t2 * direction;
    Ray2<double> ray(origin, direction);
    Segment2<double> segment(p0, p1);
    TIQuery<double, Ray2<double>, Segment2<double>> query;
    auto r = query(ray, segment);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

// ============================ IntrRay2AlignedBox2 ========================

ORACLE_CASE("IntrRay2AlignedBox2.test")
{
    int mode = io.index() % 2;
    auto ray = Ry<2>(io, mode, 3, 4.0);
    auto box = ABox<2>(io, mode, 3, 4.0);
    TIQuery<double, Ray2<double>, AlignedBox2<double>> query;
    auto r = query(ray, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay2AlignedBox2.find")
{
    int mode = io.index() % 2;
    auto ray = Ry<2>(io, mode, 3, 4.0);
    auto box = ABox<2>(io, mode, 3, 4.0);
    FIQuery<double, Ray2<double>, AlignedBox2<double>> query;
    auto r = query(ray, box);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int32_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

ORACLE_CASE("IntrRay2AlignedBox2.doQuery.ti")
{
    int mode = io.index() % 2;
    auto rayOrigin = Pt<2>(io, mode, 3, 4.0);
    auto rayDirection = (mode == 0 ? io.latticeDir<2>(-2, 2) : io.unit<2>());
    auto boxExtent = Extent<2>(io, mode, 3, 3.0);
    Expose<TIQuery<double, Ray2<double>, AlignedBox2<double>>> query;
    TIQuery<double, Ray2<double>, AlignedBox2<double>>::Result r{};
    query.DoQuery(rayOrigin, rayDirection, boxExtent, r);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay2AlignedBox2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto rayOrigin = Pt<2>(io, mode, 3, 4.0);
    auto rayDirection = (mode == 0 ? io.latticeDir<2>(-2, 2) : io.unit<2>());
    auto boxExtent = Extent<2>(io, mode, 3, 3.0);
    Expose<FIQuery<double, Ray2<double>, AlignedBox2<double>>> query;
    FIQuery<double, Ray2<double>, AlignedBox2<double>>::Result r{};
    query.DoQuery(rayOrigin, rayDirection, boxExtent, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// ============================ IntrRay3AlignedBox3 ========================

ORACLE_CASE("IntrRay3AlignedBox3.test")
{
    int mode = io.index() % 2;
    auto ray = Ry<3>(io, mode, 3, 4.0);
    auto box = ABox<3>(io, mode, 3, 4.0);
    TIQuery<double, Ray3<double>, AlignedBox3<double>> query;
    auto r = query(ray, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3AlignedBox3.find")
{
    int mode = io.index() % 2;
    auto ray = Ry<3>(io, mode, 3, 4.0);
    auto box = ABox<3>(io, mode, 3, 4.0);
    FIQuery<double, Ray3<double>, AlignedBox3<double>> query;
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

ORACLE_CASE("IntrRay3AlignedBox3.doQuery.ti")
{
    int mode = io.index() % 2;
    auto rayOrigin = Pt<3>(io, mode, 3, 4.0);
    auto rayDirection = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    auto boxExtent = Extent<3>(io, mode, 3, 3.0);
    Expose<TIQuery<double, Ray3<double>, AlignedBox3<double>>> query;
    TIQuery<double, Ray3<double>, AlignedBox3<double>>::Result r{};
    query.DoQuery(rayOrigin, rayDirection, boxExtent, r);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3AlignedBox3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto rayOrigin = Pt<3>(io, mode, 3, 4.0);
    auto rayDirection = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    auto boxExtent = Extent<3>(io, mode, 3, 3.0);
    Expose<FIQuery<double, Ray3<double>, AlignedBox3<double>>> query;
    FIQuery<double, Ray3<double>, AlignedBox3<double>>::Result r{};
    query.DoQuery(rayOrigin, rayDirection, boxExtent, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// ============================= IntrRay2Triangle2 =========================

ORACLE_CASE("IntrRay2Triangle2.test")
{
    int mode = io.index() % 2;
    auto ray = Ry<2>(io, mode, 3, 4.0);
    auto triangle = Tri<2>(io, mode, 3, 4.0);
    TIQuery<double, Ray2<double>, Triangle2<double>> query;
    auto r = query(ray, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay2Triangle2.find")
{
    int mode = io.index() % 2;
    auto ray = Ry<2>(io, mode, 3, 4.0);
    auto triangle = Tri<2>(io, mode, 3, 4.0);
    FIQuery<double, Ray2<double>, Triangle2<double>> query;
    auto r = query(ray, triangle);
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

ORACLE_CASE("IntrRay2Triangle2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto origin = Pt<2>(io, mode, 3, 4.0);
    auto direction = (mode == 0 ? io.latticeDir<2>(-2, 2) : io.unit<2>());
    auto triangle = Tri<2>(io, mode, 3, 4.0);
    Expose<FIQuery<double, Ray2<double>, Triangle2<double>>> query;
    FIQuery<double, Ray2<double>, Triangle2<double>>::Result r{};
    query.DoQuery(origin, direction, triangle, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// Uniform rays hit a triangle a few percent of the time and never on a
// boundary. Here the ray is aimed at v0 + (k1*E1 + k2*E2)/4 on an integer
// triangle with an integer direction and origin = target - t*direction, so the
// clipping comparisons of IntrLine2Triangle2::DoQuery and the ray's
// [0,+infinity) interval clip are evaluated at exact equality.
ORACLE_CASE("IntrRay2Triangle2.find.throughPoint")
{
    Triangle2<double> triangle{};
    triangle.v[0] = io.latticeVec<2>(-4, 4);
    triangle.v[1] = io.latticeVec<2>(-4, 4);
    triangle.v[2] = io.latticeVec<2>(-4, 4);
    double k1 = 0.25 * static_cast<double>(io.rawInteger(0, 5));
    double k2 = 0.25 * static_cast<double>(io.rawInteger(0, 5));
    double t = 0.5 * static_cast<double>(io.rawInteger(-4, 6));
    Vector2<double> target = triangle.v[0] + k1 * (triangle.v[1] - triangle.v[0])
        + k2 * (triangle.v[2] - triangle.v[0]);
    Vector2<double> direction = io.latticeDir<2>(-2, 2);
    Vector2<double> origin = target - t * direction;
    io.givenVec(origin);
    Ray2<double> ray(origin, direction);
    FIQuery<double, Ray2<double>, Triangle2<double>> query;
    auto r = query(ray, triangle);
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

// ============================= IntrRay3Rectangle3 ========================

ORACLE_CASE("IntrRay3Rectangle3.test")
{
    int mode = io.index() % 2;
    auto ray = Ry<3>(io, mode, 3, 4.0);
    auto rectangle = Rect(io, mode, 3, 3.0);
    TIQuery<double, Ray3<double>, Rectangle3<double>> query;
    auto r = query(ray, rectangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3Rectangle3.find")
{
    int mode = io.index() % 2;
    auto ray = Ry<3>(io, mode, 3, 4.0);
    auto rectangle = Rect(io, mode, 3, 3.0);
    FIQuery<double, Ray3<double>, Rectangle3<double>> query;
    auto r = query(ray, rectangle);
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

// The rectangle is axis aligned on the lattice and the ray is aimed at
// C + (a/2)*e0*W0 + (b/2)*e1*W1 with a, b in [-3,3], so the |rectCoord| <=
// extent comparisons are evaluated at exact equality (including the corner
// and edge cases), and t = 0 puts the ray origin exactly on the rectangle.
ORACLE_CASE("IntrRay3Rectangle3.find.throughPoint")
{
    Rectangle3<double> rectangle{};
    rectangle.center = io.latticeVec<3>(-3, 3);
    auto axis = Frame3(io, 0);
    rectangle.axis[0] = axis[0];
    rectangle.axis[1] = axis[1];
    rectangle.extent = Extent<2>(io, 0, 3, 3.0);
    double a = 0.5 * static_cast<double>(io.rawInteger(-3, 3));
    double b = 0.5 * static_cast<double>(io.rawInteger(-3, 3));
    double t = 0.5 * static_cast<double>(io.rawInteger(-4, 6));
    Vector3<double> target = rectangle.center
        + (a * rectangle.extent[0]) * rectangle.axis[0]
        + (b * rectangle.extent[1]) * rectangle.axis[1];
    Vector3<double> direction = io.latticeDir<3>(-2, 2);
    Vector3<double> origin = target - t * direction;
    io.givenVec(origin);
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, Rectangle3<double>> query;
    auto r = query(ray, rectangle);
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

// =============================== arcs ====================================

namespace
{
    // The twelve lattice points on the circle of radius 5 centered at the
    // origin. Arc endpoints drawn from this table lie exactly on the circle,
    // which is what the cocircular branches of the arc queries require.
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

    // An arc. Mode 0 puts the arc on the lattice circle of radius 5 about
    // 'center' so that two such arcs are cocircular and their endpoints
    // compare exactly.
    Arc2<double> Arc(oracle::Ctx& io, int mode, Vector2<double> const& center)
    {
        Arc2<double> arc{};
        if (mode == 0)
        {
            arc.center = io.givenVec(center);
            arc.radius = io.given(5.0);
            arc.end[0] = io.givenVec(center + CirclePoint5(io.rawInteger(0, 11)));
            arc.end[1] = io.givenVec(center + CirclePoint5(io.rawInteger(0, 11)));
        }
        else
        {
            arc.center = io.vec<2>(-3.0, 3.0);
            arc.radius = io.real(1.0, 4.0);
            arc.end[0] = io.vec<2>(-6.0, 6.0);
            arc.end[1] = io.vec<2>(-6.0, 6.0);
        }
        return arc;
    }

    void OutArc(oracle::Ctx& io, Arc2<double> const& arc)
    {
        io.outVec(arc.center);
        io.outReal(arc.radius);
        io.outVec(arc.end[0]);
        io.outVec(arc.end[1]);
    }
}

// Mode 0 gives two cocircular arcs (the eight COCIRCULAR_* configurations);
// mode 1 gives two unrelated circles (the NONCOCIRCULAR_* configurations and
// NO_INTERSECTION).
ORACLE_CASE("IntrArc2Arc2.find")
{
    int mode = io.index() % 2;
    Vector2<double> center = io.latticeVec<2>(-3, 3);
    auto arc0 = Arc(io, mode, center);
    auto arc1 = Arc(io, mode, center);
    FIQuery<double, Arc2<double>, Arc2<double>> query;
    auto r = query(arc0, arc1);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.configuration));
    io.outVec(r.point[0]);
    io.outVec(r.point[1]);
    OutArc(io, r.arc[0]);
    OutArc(io, r.arc[1]);
}

ORACLE_CASE("IntrCircle2Arc2.find")
{
    int mode = io.index() % 2;
    Vector2<double> center = io.latticeVec<2>(-3, 3);
    auto arc = Arc(io, mode, center);
    Circle2<double> circle{};
    if (mode == 0)
    {
        // The circle of the arc, so the "arc is on the circle" branch is
        // reached with exact equality.
        circle.center = io.givenVec(center);
        circle.radius = io.given(5.0);
    }
    else
    {
        circle.center = io.vec<2>(-3.0, 3.0);
        circle.radius = io.real(1.0, 4.0);
    }
    FIQuery<double, Circle2<double>, Arc2<double>> query;
    auto r = query(circle, arc);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outVec(r.point[0]);
    io.outVec(r.point[1]);
    OutArc(io, r.arc);
}

// ============================== IntrDisk2Sector2 =========================

ORACLE_CASE("IntrDisk2Sector2.test")
{
    int mode = io.index() % 2;
    auto disk = Circ(io, mode, 3, 3.0);
    auto sector = Sect(io, mode, 3, 3.0);
    TIQuery<double, Circle2<double>, Sector2<double>> query;
    auto r = query(disk, sector);
    io.outBool(r.intersect);
}

// ========================== IntrOrientedBox2Sector2 ======================

namespace
{
    // Upstream's IntrOrientedBox2Sector2 assigns the clip result's polygon
    // unconditionally. IntrHalfspace2Polygon2 reports intersect = true with an
    // EMPTY polygon when no clipping was necessary, so the assignment discards
    // the working polygon and the query then reports "no intersection"
    // (gtengine-js issue #200; the port keeps the polygon). This probe
    // replicates upstream's control flow, reports whether either clip is such
    // a no-op ('sound' is false when one is), and reports whether keeping the
    // polygon changes the answer ('differs').
    struct Ob2S2Probe
    {
        bool sound;
        bool differs;
    };

    std::vector<Vector2<double>> Ob2S2Clip(
        std::vector<Vector2<double>> const& polygon,
        Vector2<double> const& normal, Vector2<double> const& vertex,
        bool& wasNoOp)
    {
        Halfspace<2, double> halfspace{};
        halfspace.normal = normal;
        halfspace.constant = Dot(normal, vertex);
        FIQuery<double, Halfspace<2, double>, std::vector<Vector2<double>>> hpQuery;
        auto hpResult = hpQuery(halfspace, polygon);
        wasNoOp = (hpResult.intersect && hpResult.polygon.empty());
        return hpResult.polygon;
    }

    bool Ob2S2Reaches(std::vector<Vector2<double>> const& polygon,
        Sector2<double> const& sector)
    {
        DCPQuery<double, Vector2<double>, Segment2<double>> psQuery;
        int32_t const numVertices = static_cast<int32_t>(polygon.size());
        if (numVertices >= 2)
        {
            for (int32_t i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++)
            {
                Segment2<double> segment(polygon[i0], polygon[i1]);
                auto psResult = psQuery(sector.vertex, segment);
                if (psResult.distance <= sector.radius)
                {
                    return true;
                }
            }
        }
        return false;
    }

    Ob2S2Probe ProbeOb2S2(OrientedBox2<double> const& box,
        Sector2<double> const& sector)
    {
        Vector2<double> CmV = box.center - sector.vertex;
        Vector2<double> P{ Dot(box.axis[0], CmV), Dot(box.axis[1], CmV) };
        if (std::fabs(P[0]) <= box.extent[0] && std::fabs(P[1]) <= box.extent[1])
        {
            return { true, false };
        }

        Vector2<double> U0
        {
            +sector.cosAngle * sector.direction[0] + sector.sinAngle * sector.direction[1],
            -sector.sinAngle * sector.direction[0] + sector.cosAngle * sector.direction[1]
        };
        Vector2<double> N0 = Perp(U0);
        double prjcen0 = Dot(N0, CmV);
        double radius0 = box.extent[0] * std::fabs(Dot(N0, box.axis[0]))
            + box.extent[1] * std::fabs(Dot(N0, box.axis[1]));
        if (prjcen0 > radius0)
        {
            return { true, false };
        }

        Vector2<double> U1
        {
            +sector.cosAngle * sector.direction[0] - sector.sinAngle * sector.direction[1],
            +sector.sinAngle * sector.direction[0] + sector.cosAngle * sector.direction[1]
        };
        Vector2<double> N1 = -Perp(U1);
        double prjcen1 = Dot(N1, CmV);
        double radius1 = box.extent[0] * std::fabs(Dot(N1, box.axis[0]))
            + box.extent[1] * std::fabs(Dot(N1, box.axis[1]));
        if (prjcen1 > radius1)
        {
            return { true, false };
        }

        Vector2<double> e0U0 = box.extent[0] * box.axis[0];
        Vector2<double> e1U1 = box.extent[1] * box.axis[1];
        std::vector<Vector2<double>> boxPolygon;
        boxPolygon.push_back(box.center - e0U0 - e1U1);
        boxPolygon.push_back(box.center + e0U0 - e1U1);
        boxPolygon.push_back(box.center + e0U0 + e1U1);
        boxPolygon.push_back(box.center - e0U0 + e1U1);

        std::vector<Vector2<double>> upstream = boxPolygon, fixed = boxPolygon;
        bool noOp = false, sound = true;
        if (prjcen0 >= -radius0)
        {
            upstream = Ob2S2Clip(upstream, -N0, sector.vertex, noOp);
            sound = sound && !noOp;
            auto clipped = Ob2S2Clip(fixed, -N0, sector.vertex, noOp);
            if (!noOp) { fixed = clipped; }
        }
        if (prjcen1 >= -radius1)
        {
            upstream = Ob2S2Clip(upstream, -N1, sector.vertex, noOp);
            sound = sound && !noOp;
            auto clipped = Ob2S2Clip(fixed, -N1, sector.vertex, noOp);
            if (!noOp) { fixed = clipped; }
        }

        bool upstreamIntersect = Ob2S2Reaches(upstream, sector);
        bool fixedIntersect = Ob2S2Reaches(fixed, sector);
        return { sound, upstreamIntersect != fixedIntersect };
    }

    // Unrecorded draws for a box and a sector; the caller records the values.
    void RawBoxSector(oracle::Ctx& io, int mode, OrientedBox2<double>& box,
        Sector2<double>& sector)
    {
        for (int i = 0; i < 2; ++i)
        {
            box.center[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-3.0, 3.0));
            box.extent[i] = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
                : io.raw(0.25, 3.0));
        }
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
        box.axis[0] = Vector2<double>{ c, s };
        box.axis[1] = Vector2<double>{ -s, c };

        for (int i = 0; i < 2; ++i)
        {
            sector.vertex[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-3.0, 3.0));
        }
        sector.radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.5, 3.0));
        double dc{}, ds{};
        if (mode == 0)
        {
            int k = io.rawInteger(0, 3);
            dc = (k == 0 ? 1.0 : (k == 1 ? 0.0 : (k == 2 ? -1.0 : 0.0)));
            ds = (k == 0 ? 0.0 : (k == 1 ? 1.0 : (k == 2 ? 0.0 : -1.0)));
        }
        else
        {
            double angle = io.raw(-3.141592653589793, 3.141592653589793);
            dc = std::cos(angle);
            ds = std::sin(angle);
        }
        sector.direction = Vector2<double>{ dc, ds };
        if (mode == 0)
        {
            int k = io.rawInteger(0, 2);
            sector.cosAngle = (k == 0 ? 0.0 : (k == 1 ? 0.6 : 0.8));
            sector.sinAngle = (k == 0 ? 1.0 : (k == 1 ? 0.8 : 0.6));
        }
        else
        {
            double angle = io.raw(0.05, 1.5707963267948966);
            sector.cosAngle = std::cos(angle);
            sector.sinAngle = std::sin(angle);
        }
    }

    void RecordBoxSector(oracle::Ctx& io, OrientedBox2<double> const& box,
        Sector2<double> const& sector)
    {
        io.givenVec(box.center);
        io.givenVec(box.axis[0]);
        io.givenVec(box.axis[1]);
        io.givenVec(box.extent);
        io.givenVec(sector.vertex);
        io.given(sector.radius);
        io.givenVec(sector.direction);
        io.given(sector.cosAngle);
        io.given(sector.sinAngle);
    }
}

// The generator rejects the configurations on which upstream discards the
// working polygon, i.e. exactly the inputs on which the port deliberately
// deviates (see the deviation case below).
ORACLE_CASE("IntrOrientedBox2Sector2.test")
{
    int mode = io.index() % 2;
    OrientedBox2<double> box{};
    Sector2<double> sector{};
    for (;;)
    {
        RawBoxSector(io, mode, box, sector);
        if (ProbeOb2S2(box, sector).sound) { break; }
    }
    RecordBoxSector(io, box, sector);
    TIQuery<double, OrientedBox2<double>, Sector2<double>> query;
    auto r = query(box, sector);
    io.outBool(r.intersect);
}

// Deliberate deviation: configurations on which upstream's unconditional
// 'polygon = std::move(hpResult.polygon)' throws the clipped box away and the
// reported answer changes. gtengine-js issue #200.
ORACLE_CASE("IntrOrientedBox2Sector2.test.clipDeviation")
{
    OrientedBox2<double> box{};
    Sector2<double> sector{};
    for (;;)
    {
        RawBoxSector(io, 1, box, sector);
        auto probe = ProbeOb2S2(box, sector);
        if (!probe.sound && probe.differs) { break; }
    }
    RecordBoxSector(io, box, sector);
    TIQuery<double, OrientedBox2<double>, Sector2<double>> query;
    auto r = query(box, sector);
    io.outBool(r.intersect);
}

// ============================ IntrCapsule3Capsule3 =======================

ORACLE_CASE("IntrCapsule3Capsule3.test")
{
    int mode = io.index() % 2;
    auto capsule0 = Caps(io, mode, 3, 3.0);
    auto capsule1 = Caps(io, mode, 3, 3.0);
    TIQuery<double, Capsule3<double>, Capsule3<double>> query;
    auto r = query(capsule0, capsule1);
    io.outBool(r.intersect);
}

// =========================== IntrHalfspace3Capsule3 ======================

ORACLE_CASE("IntrHalfspace3Capsule3.test")
{
    int mode = io.index() % 2;
    auto capsule = Caps(io, mode, 3, 3.0);
    Halfspace3<double> halfspace{};
    halfspace.normal = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    // In lattice mode the boundary passes exactly through an endpoint, which
    // reaches the max(e0,e1) + radius == 0 boundary.
    halfspace.constant = (mode == 0
        ? io.given(Dot(halfspace.normal, capsule.segment.p[0]) + capsule.radius)
        : io.real(-3.0, 3.0));
    TIQuery<double, Halfspace3<double>, Capsule3<double>> query;
    auto r = query(halfspace, capsule);
    io.outBool(r.intersect);
}

// =========================== IntrHalfspace3Cylinder3 =====================

// Upstream computes root = sqrt(max(1, 1 - Dot(N,W)^2)), and since
// 1 - x^2 <= 1 the max always selects 1, so the documented
// r*sqrt(1-Dot(N,W)^2) term is never computed. The port evaluates
// sqrt(max(0, 1 - Dot(N,W)^2)) (gtengine-js issue #197). Upstream is sound
// exactly when Dot(N,W) is zero, so the main case makes the halfspace normal
// orthogonal to the cylinder axis: in lattice mode by construction from
// coordinate axes (Dot exactly 0), in uniform mode from an orthonormal frame
// (Dot within one ulp of 0, for which 1 - Dot^2 rounds to exactly 1).
ORACLE_CASE("IntrHalfspace3Cylinder3.test.perpendicular")
{
    int mode = io.index() % 2;
    Vector3<double> axisDirection{};
    Vector3<double> normal{};
    if (mode == 0)
    {
        int k = io.rawInteger(0, 2);
        int m = (k + 1 + io.rawInteger(0, 1)) % 3;
        axisDirection.MakeZero();
        axisDirection[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        normal.MakeZero();
        normal[m] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    }
    else
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
        double angle = io.raw(-3.141592653589793, 3.141592653589793);
        axisDirection = basis[0];
        normal = std::cos(angle) * basis[1] + std::sin(angle) * basis[2];
        Normalize(normal);
    }
    Vector3<double> origin = Pt<3>(io, mode, 3, 3.0);
    io.givenVec(axisDirection);
    double radius = (mode == 0 ? io.lattice(1, 3) : io.real(0.25, 3.0));
    double height = (mode == 0 ? io.lattice(1, 6) : io.real(0.5, 6.0));
    Line3<double> axis(origin, axisDirection);
    Cylinder3<double> cylinder(axis, radius, height);
    Halfspace3<double> halfspace{};
    halfspace.normal = io.givenVec(normal);
    halfspace.constant = io.real(-6.0, 6.0);
    TIQuery<double, Halfspace3<double>, Cylinder3<double>> query;
    auto r = query(halfspace, cylinder);
    io.outBool(r.intersect);
}

// Deliberate deviation: the halfspace constant is placed so that upstream's
// root = 1 reports an intersection and the corrected root = sqrt(1-Dot(N,W)^2)
// does not. gtengine-js issue #197.
ORACLE_CASE("IntrHalfspace3Cylinder3.test.rootDeviation")
{
    Vector3<double> origin = io.vec<3>(-3.0, 3.0);
    Vector3<double> axisDirection = io.unit<3>();
    double radius = io.real(0.5, 3.0);
    double height = io.real(0.5, 6.0);
    Vector3<double> normal = io.unit<3>();
    double absNdW = std::fabs(Dot(normal, axisDirection));
    double correct = std::sqrt(std::max(0.0, 1.0 - absNdW * absNdW));
    // center + radius*1 + (h/2)*absNdW = gap > 0 but
    // center + radius*correct + (h/2)*absNdW < 0.
    double gap = 0.5 * radius * (1.0 - correct);
    double center = gap - radius - 0.5 * height * absNdW;
    double constant = io.given(Dot(normal, origin) - center);
    Line3<double> axis(origin, axisDirection);
    Cylinder3<double> cylinder(axis, radius, height);
    Halfspace3<double> halfspace{};
    halfspace.normal = normal;
    halfspace.constant = constant;
    TIQuery<double, Halfspace3<double>, Cylinder3<double>> query;
    auto r = query(halfspace, cylinder);
    io.outBool(r.intersect);
}

// ============================= IntrLine3Capsule3 =========================

namespace
{
    // An independent reference for the line/capsule intersection interval.
    // The solid capsule is the union of the solid finite cylinder and the two
    // solid end balls, and a line meets a convex set in an interval, so the
    // interval is [min,max] over the candidate roots that lie on the capsule
    // boundary: the roots of the infinite-cylinder quadratic with |z| <= e and
    // the roots of the two end-sphere quadratics. Each candidate is a point of
    // the capsule and each endpoint of the true interval is one of the
    // candidates, so the extremes are the true endpoints. Upstream instead
    // accumulates roots one at a time, returns as soon as two are accepted and
    // sets 'intersect' only in the cylinder-wall branches, which makes it
    // disagree with this reference on two families of inputs (gtengine-js
    // issue #461, both fixed in the port). The main cases restrict their
    // generators to inputs where upstream agrees with the reference; the two
    // deviation cases construct inputs where it does not.
    struct CapsuleInterval
    {
        bool intersect;
        int32_t numIntersections;
        std::array<double, 2> parameter;
    };

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

    // Unrecorded draws for a unit-direction line and a capsule; the caller
    // records the values.
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
}

ORACLE_CASE("IntrLine3Capsule3.test")
{
    int mode = io.index() % 2;
    auto line = UnitLn<3>(io, mode, 3, 3.0);
    auto capsule = Caps(io, mode, 3, 3.0);
    TIQuery<double, Line3<double>, Capsule3<double>> query;
    auto r = query(line, capsule);
    io.outBool(r.intersect);
}

// The generator rejects the inputs on which upstream disagrees with the
// reference interval above, that is exactly the inputs on which the port
// deliberately deviates (issue #461), and the degenerate zero-length capsule
// segment whose frame is undefined (a documented upstream limitation the port
// preserves; both sides then produce NaN, which the probe also rejects).
ORACLE_CASE("IntrLine3Capsule3.find")
{
    int mode = io.index() % 2;
    Vector3<double> lineOrigin{}, lineDirection{}, p0{}, p1{};
    double radius{};
    for (;;)
    {
        RawLineCapsule(io, mode, lineOrigin, lineDirection, p0, p1, radius);
        Segment3<double> probeSegment(p0, p1);
        Capsule3<double> probeCapsule(probeSegment, radius);
        if (Line3Capsule3IsSound(lineOrigin, lineDirection, probeCapsule)) { break; }
    }
    io.givenVec(lineOrigin);
    io.givenVec(lineDirection);
    io.givenVec(p0);
    io.givenVec(p1);
    io.given(radius);
    Segment3<double> segment(p0, p1);
    Capsule3<double> capsule(segment, radius);
    Line3<double> line(lineOrigin, lineDirection);
    FIQuery<double, Line3<double>, Capsule3<double>> query;
    auto r = query(line, capsule);
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

ORACLE_CASE("IntrLine3Capsule3.doQuery.fi")
{
    int mode = io.index() % 2;
    Vector3<double> lineOrigin{}, lineDirection{}, p0{}, p1{};
    double radius{};
    for (;;)
    {
        RawLineCapsule(io, mode, lineOrigin, lineDirection, p0, p1, radius);
        Segment3<double> probeSegment(p0, p1);
        Capsule3<double> probeCapsule(probeSegment, radius);
        if (Line3Capsule3IsSound(lineOrigin, lineDirection, probeCapsule)) { break; }
    }
    io.givenVec(lineOrigin);
    io.givenVec(lineDirection);
    io.givenVec(p0);
    io.givenVec(p1);
    io.given(radius);
    Segment3<double> segment(p0, p1);
    Capsule3<double> capsule(segment, radius);
    Expose<FIQuery<double, Line3<double>, Capsule3<double>>> query;
    FIQuery<double, Line3<double>, Capsule3<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, capsule, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

namespace
{
    // A capsule whose segment is parallel to a coordinate axis, in capsule
    // coordinates {U,V,W} = {(0,1,0), (-1,0,0), (0,0,1)} after the axis
    // permutation, with integer half-length e and radius r.
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

// Deliberate deviation 1: the line is tangent to the bottom end sphere at its
// pole, so the single accepted root comes from a hemisphere. Upstream copies
// the parameter but never sets 'intersect' (gtengine-js issue #461 item 2).
ORACLE_CASE("IntrLine3Capsule3.doQuery.capTangentDeviation")
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
    Vector3<double> lineOrigin = center - (e + r) * ek;
    io.givenVec(lineOrigin);
    io.givenVec(em);
    Expose<FIQuery<double, Line3<double>, Capsule3<double>>> query;
    FIQuery<double, Line3<double>, Capsule3<double>>::Result result{};
    query.DoQuery(lineOrigin, em, capsule, result);
    io.outBool(result.intersect);
    io.outInt(static_cast<int32_t>(result.numIntersections));
    io.outReal(result.parameter[0]);
    io.outReal(result.parameter[1]);
}

// Deliberate deviation 2: the line meets the capsule wall exactly on the
// cap-junction circle z = -e. Upstream's wall test and bottom-hemisphere test
// both accept that one root, the count reaches two and the early return
// collapses the intersection to a single point, discarding the true far
// endpoint (gtengine-js issue #461 item 3).
ORACLE_CASE("IntrLine3Capsule3.doQuery.junctionDeviation")
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
    Vector3<double> lineOrigin{ center[0], center[1] + r, center[2] - e };
    Vector3<double> lineDirection{ 0.0, 0.5, -std::sqrt(0.75) };
    io.givenVec(lineOrigin);
    io.givenVec(lineDirection);
    Expose<FIQuery<double, Line3<double>, Capsule3<double>>> query;
    FIQuery<double, Line3<double>, Capsule3<double>>::Result result{};
    query.DoQuery(lineOrigin, lineDirection, capsule, result);
    io.outBool(result.intersect);
    io.outInt(static_cast<int32_t>(result.numIntersections));
    io.outReal(result.parameter[0]);
    io.outReal(result.parameter[1]);
}

// ============================ IntrLine3Cylinder3 =========================

ORACLE_CASE("IntrLine3Cylinder3.find")
{
    int mode = io.index() % 2;
    auto line = UnitLn<3>(io, mode, 3, 3.0);
    auto cylinder = Cyl(io, mode, 3, 3.0);
    FIQuery<double, Line3<double>, Cylinder3<double>> query;
    auto r = query(line, cylinder);
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

ORACLE_CASE("IntrLine3Cylinder3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto lineOrigin = Pt<3>(io, mode, 3, 3.0);
    auto lineDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    auto cylinder = Cyl(io, mode, 3, 3.0);
    Expose<FIQuery<double, Line3<double>, Cylinder3<double>>> query;
    FIQuery<double, Line3<double>, Cylinder3<double>>::Result r{};
    query.DoQuery(lineOrigin, lineDirection, cylinder, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// The axis-parallel branch (|Dot(W,D)| == 1) and the axis-perpendicular branch
// (D[2] == 0) are unreachable with uniform directions. Here the cylinder axis
// and the line direction are coordinate axes, so both are reached exactly, and
// the line origin is placed at the radial distance r + d from the axis so the
// tangency (discriminant exactly zero) also occurs.
ORACLE_CASE("IntrLine3Cylinder3.find.axisAligned")
{
    Vector3<double> axisOrigin{};
    for (int i = 0; i < 3; ++i)
    {
        axisOrigin[i] = static_cast<double>(io.rawInteger(-3, 3));
    }
    int k = io.rawInteger(0, 2);
    int m = io.rawInteger(0, 2);
    double radius = static_cast<double>(io.rawInteger(1, 3));
    double height = static_cast<double>(io.rawInteger(1, 6));
    double d = static_cast<double>(io.rawInteger(-1, 1));
    double s = static_cast<double>(io.rawInteger(-4, 4));
    Vector3<double> ek{}, em{}, en{};
    ek.MakeZero();
    ek[k] = 1.0;
    em.MakeZero();
    em[m] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    en.MakeZero();
    en[(k + 1) % 3] = 1.0;
    Vector3<double> lineOrigin = axisOrigin + (radius + d) * en + s * ek;
    io.givenVec(axisOrigin);
    io.givenVec(ek);
    io.given(radius);
    io.given(height);
    io.givenVec(lineOrigin);
    io.givenVec(em);
    Line3<double> axis(axisOrigin, ek);
    Cylinder3<double> cylinder(axis, radius, height);
    Line3<double> line(lineOrigin, em);
    FIQuery<double, Line3<double>, Cylinder3<double>> query;
    auto r = query(line, cylinder);
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

// ============================== IntrLine3Torus3 ==========================

// This is the only find-intersection query in the group whose compared path
// calls the C math library: RootsPolynomial::SolveQuartic uses sqrt, cbrt and
// trigonometric functions, and Torus3::GetParameters uses atan2. The case is
// compared with the default scaled tolerance.
ORACLE_CASE("IntrLine3Torus3.find")
{
    Torus3<double> torus{};
    torus.center = io.latticeVec<3>(-2, 2);
    auto frame = Frame3(io, io.index() % 2);
    torus.direction0 = frame[1];
    torus.direction1 = frame[2];
    torus.normal = frame[0];
    torus.radius1 = io.real(0.5, 1.5);
    torus.radius0 = torus.radius1 + io.real(0.5, 2.0);
    Vector3<double> lineOrigin = io.vec<3>(-4.0, 4.0);
    Vector3<double> lineDirection = io.unit<3>();
    Line3<double> line(lineOrigin, lineDirection);
    FIQuery<double, Line3<double>, Torus3<double>> query;
    auto r = query(line, torus);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    for (size_t i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.lineParameter[i]);
        io.outVec(r.point[i]);
        io.outReal(r.torusParameter[i][0]);
        io.outReal(r.torusParameter[i][1]);
    }
}

// ========================== IntrCylinder3Cylinder3 =======================

// The four analytic separating-axis tests use only arithmetic, but the
// hemisphere sampling loop calls cos and sin, so this case is compared with
// the default scaled tolerance. numThreads is always 1: the upstream
// multithreaded path reports whichever thread found a separating direction
// first, which is not deterministic, and the port always runs the
// single-threaded loop.
ORACLE_CASE("IntrCylinder3Cylinder3.test")
{
    int mode = io.index() % 2;
    int numTheta = io.integer(2, 6);
    int numPhi = io.integer(2, 5);
    auto cylinder0 = Cyl(io, mode, 2, 2.0);
    auto cylinder1 = Cyl(io, mode, 2, 2.0);
    TIQuery<double, Cylinder3<double>, Cylinder3<double>> query(1,
        static_cast<size_t>(numTheta), static_cast<size_t>(numPhi));
    auto r = query(cylinder0, cylinder1);
    io.outBool(r.separated);
    io.outVec(r.separatingDirection);
}

// Parallel axes take the branch that never calls cos or sin, so this variant
// is exact. The coincident-origin early return (Length(Delta) == 0) is also
// reached.
ORACLE_CASE("IntrCylinder3Cylinder3.test.parallel")
{
    Vector3<double> origin0 = io.latticeVec<3>(-3, 3);
    Vector3<double> direction = AxisUnit<3>(io);
    double radius0 = io.lattice(1, 3);
    double height0 = io.lattice(1, 6);
    Vector3<double> origin1 = io.latticeVec<3>(-3, 3);
    bool flip = io.boolean();
    double radius1 = io.lattice(1, 3);
    double height1 = io.lattice(1, 6);
    Vector3<double> direction1 = (flip ? -direction : direction);
    Line3<double> axis0(origin0, direction);
    Line3<double> axis1(origin1, direction1);
    Cylinder3<double> cylinder0(axis0, radius0, height0);
    Cylinder3<double> cylinder1(axis1, radius1, height1);
    TIQuery<double, Cylinder3<double>, Cylinder3<double>> query(1, 4, 3);
    auto r = query(cylinder0, cylinder1);
    io.outBool(r.separated);
    io.outVec(r.separatingDirection);
}

// Deliberate deviation: Cylinder3 marks an infinite cylinder with the
// sentinel height = -1, which upstream feeds into the finite-cylinder
// formulas; the port asserts IsFinite instead (gtengine-js issue #197).
ORACLE_CASE("IntrCylinder3Cylinder3.test.infiniteDeviation")
{
    auto cylinder0 = Cyl(io, 1, 2, 2.0);
    Vector3<double> origin1 = io.vec<3>(-2.0, 2.0);
    Vector3<double> direction1 = io.unit<3>();
    double radius1 = io.real(0.25, 2.0);
    double height1 = io.given(-1.0);
    Line3<double> axis1(origin1, direction1);
    Cylinder3<double> cylinder1(axis1, radius1, height1);
    TIQuery<double, Cylinder3<double>, Cylinder3<double>> query(1, 4, 3);
    auto r = query(cylinder0, cylinder1);
    io.outBool(r.separated);
    io.outVec(r.separatingDirection);
}

// ========================= IntrCanonicalBox3Cylinder3 ====================

namespace
{
    // A cylinder axis direction with a controlled number of zero components:
    // zeros = 2 reaches DoQueryTwoZeros, zeros = 1 reaches DoQueryOneZero and
    // zeros = 0 reaches DoQueryNoZeros. The draws are unrecorded; the caller
    // records the resulting direction.
    void RawCylinderDir(oracle::Ctx& io, int zeros, Vector3<double>& v)
    {
        v.MakeZero();
        if (zeros == 2)
        {
            int k = io.rawInteger(0, 2);
            v[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        }
        else if (zeros == 1)
        {
            int k = io.rawInteger(0, 2);
            double len = 0.0;
            do
            {
                v[(k + 1) % 3] = io.raw(-1.0, 1.0);
                v[(k + 2) % 3] = io.raw(-1.0, 1.0);
                len = Length(v);
            } while (len < 0.25 || len > 1.0);
            Normalize(v);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { v[i] = io.raw(-1.0, 1.0); }
                len = Length(v);
            } while (len < 0.25 || len > 1.0 || v[0] == 0.0 || v[1] == 0.0 || v[2] == 0.0);
            Normalize(v);
        }
    }

    // ---- upstream's DoQueryNoZeros with the (U1,-D) sign typo corrected ----
    // The typo is documented in docs/UPSTREAM-FINDINGS.md (gtengine-js issue
    // #197): the block builds a segment endpoint from the box vertex
    // V3 = (+E[0],+E[1],-E[2]) but writes -E[1]. Every other branch of every
    // other edge builds the point from the vertex that generated its s-value,
    // replacing exactly one coordinate. The port corrects it, so the generator
    // needs the corrected answer to tell the sound inputs from the defective
    // ones. Everything else here is copied verbatim from the header.

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

    // 'requireDeviation' selects the inputs on which the corrected answer
    // differs from upstream's (the deviation case) or agrees with it (the
    // main cases).
    void CanonicalBoxCylinderCase(oracle::Ctx& io, int zeros, bool requireDeviation)
    {
        Vector3<double> extent{}, origin{}, direction{};
        double radius{}, height{};
        for (;;)
        {
            for (int i = 0; i < 3; ++i) { extent[i] = io.raw(0.5, 3.0); }
            for (int i = 0; i < 3; ++i) { origin[i] = io.raw(-4.0, 4.0); }
            RawCylinderDir(io, zeros, direction);
            radius = io.raw(0.25, 3.0);
            height = io.raw(0.5, 6.0);
            CanonicalBox3<double> probeBox(extent);
            Line3<double> probeAxis(origin, direction);
            Cylinder3<double> probeCylinder(probeAxis, radius, height);
            TIQuery<double, CanonicalBox3<double>, Cylinder3<double>> probeQuery;
            bool upstream = probeQuery(probeBox, probeCylinder).intersect;
            if ((upstream != CbcFixed(probeBox, probeCylinder, upstream)) == requireDeviation)
            {
                break;
            }
        }
        io.givenVec(extent);
        io.givenVec(origin);
        io.givenVec(direction);
        io.given(radius);
        io.given(height);
        CanonicalBox3<double> box(extent);
        Line3<double> axis(origin, direction);
        Cylinder3<double> cylinder(axis, radius, height);
        TIQuery<double, CanonicalBox3<double>, Cylinder3<double>> query;
        auto r = query(box, cylinder);
        io.outBool(r.intersect);
    }
}

ORACLE_CASE("IntrCanonicalBox3Cylinder3.test")
{
    CanonicalBoxCylinderCase(io, 0, false);
}

ORACLE_CASE("IntrCanonicalBox3Cylinder3.test.oneZero")
{
    CanonicalBoxCylinderCase(io, 1, false);
}

ORACLE_CASE("IntrCanonicalBox3Cylinder3.test.twoZeros")
{
    CanonicalBoxCylinderCase(io, 2, false);
}

// Deliberate deviation: the sign typo in the (U1,-D) block of DoQueryNoZeros
// puts a segment endpoint on the wrong box edge, so upstream reports no
// intersection for a cylinder that does meet the box. gtengine-js issue #197.
ORACLE_CASE("IntrCanonicalBox3Cylinder3.test.edgeTypoDeviation")
{
    CanonicalBoxCylinderCase(io, 0, true);
}

// ========================== IntrOrientedBox3Frustum3 =====================

ORACLE_CASE("IntrOrientedBox3Frustum3.test")
{
    int mode = io.index() % 2;
    Vector3<double> frustumOrigin = Pt<3>(io, mode, 2, 2.0);
    auto frame = Frame3(io, mode);
    double dMin = (mode == 0 ? io.lattice(1, 2) : io.real(0.5, 2.0));
    double dMax = dMin + (mode == 0 ? io.lattice(1, 4) : io.real(0.5, 4.0));
    double uBound = (mode == 0 ? io.lattice(1, 2) : io.real(0.25, 2.0));
    double rBound = (mode == 0 ? io.lattice(1, 2) : io.real(0.25, 2.0));
    Frustum3<double> frustum(frustumOrigin, frame[0], frame[1], frame[2],
        dMin, dMax, uBound, rBound);
    // Place the box near the frustum's bounding region: inside it often
    // enough that the later separating axes are reached, outside it often
    // enough that each of the early rejections fires.
    double t = io.raw(-2.0, 3.0);
    double a = io.raw(-6.0, 6.0);
    double b = io.raw(-6.0, 6.0);
    OrientedBox3<double> box{};
    Vector3<double> center = frustumOrigin + (dMin + t * (dMax - dMin)) * frame[0]
        + (a * uBound) * frame[1] + (b * rBound) * frame[2];
    box.center = io.givenVec(center);
    auto boxAxis = Frame3(io, mode);
    box.axis[0] = boxAxis[0];
    box.axis[1] = boxAxis[1];
    box.axis[2] = boxAxis[2];
    box.extent = Extent<3>(io, mode, 2, 2.0);
    TIQuery<double, OrientedBox3<double>, Frustum3<double>> query;
    auto r = query(box, frustum);
    io.outBool(r.intersect);
}
