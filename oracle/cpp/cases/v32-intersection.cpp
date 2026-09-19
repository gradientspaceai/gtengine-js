// Verify group 32 (intersection): differential cases for the Intr* headers
// listed in plan/verify-groups.json group 32.
//
// Every generator mixes an exactly representable small-lattice mode with a
// uniform mode and, where a branch needs it, a constructed mode, so that the
// touching, tangent, collinear, contained and degenerate branches are reached.
// Almost every query in this group uses only + - * / sqrt fabs and
// comparisons, so almost every case is declared exact on the TypeScript side;
// the sin/cos used to build orthonormal frames and unit directions are applied
// to *unrecorded* draws and only the resulting frame or unit vector is
// recorded as an input, so libm never enters the compared computation. The
// exceptions are the two IntrEllipse2Ellipse2 queries (bisection with std::pow,
// and RootsPolynomial::SolveQuartic), marked in their case comments.
#define ORACLE_FAMILY "v32-intersection"
#include "Oracle.h"

#include <Mathematics/IntrAlignedBox2Circle2.h>
#include <Mathematics/IntrAlignedBox3Cylinder3.h>
#include <Mathematics/IntrAlignedBox3Sphere3.h>
#include <Mathematics/IntrConvexPolygonHyperplane.h>
#include <Mathematics/IntrEllipse2Ellipse2.h>
#include <Mathematics/IntrRay3Sphere3.h>
#include <Mathematics/IntrSegment2AlignedBox2.h>
#include <Mathematics/IntrSegment2Segment2.h>
#include <Mathematics/IntrSegment2Triangle2.h>
#include <Mathematics/IntrSegment3AlignedBox3.h>
#include <Mathematics/IntrSegment3Rectangle3.h>
#include <Mathematics/IntrSegment3Sphere3.h>
#include <Mathematics/IntrSphere3Sphere3.h>
#include <Mathematics/IntrSphere3Triangle3.h>
#include <Mathematics/IntrTriangle3Cylinder3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
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

    // Orthonormal 3D frame. Mode 0 produces a signed permutation of the
    // coordinate axes (exact).
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

    template <int N>
    Segment<N, double> Sg(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> p0 = Pt<N>(io, mode, lat, range);
        Vector<N, double> p1 = Pt<N>(io, mode, lat, range);
        return Segment<N, double>(p0, p1);
    }

    // A ray with a unit-length direction, for the queries that document that
    // precondition (sphere).
    template <int N>
    Ray<N, double> UnitRy(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> origin = Pt<N>(io, mode, lat, range);
        Vector<N, double> direction = (mode == 0 ? AxisUnit<N>(io) : io.unit<N>());
        return Ray<N, double>(origin, direction);
    }

    template <int N>
    Hypersphere<N, double> Sph(oracle::Ctx& io, int mode, int lat, double range)
    {
        Hypersphere<N, double> sphere{};
        sphere.center = Pt<N>(io, mode, lat, range);
        sphere.radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        return sphere;
    }

    // Exposes the protected DoQuery helpers, which the port exports as public
    // free functions.
    template <typename Q>
    struct Expose : public Q
    {
        using Q::DoQuery;
    };
}

// ============================== IntrRay3Sphere3 ==========================

ORACLE_CASE("IntrRay3Sphere3.test")
{
    int mode = io.index() % 2;
    auto ray = UnitRy<3>(io, mode, 3, 4.0);
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    TIQuery<double, Ray3<double>, Sphere3<double>> query;
    auto r = query(ray, sphere);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3Sphere3.find")
{
    int mode = io.index() % 2;
    auto ray = UnitRy<3>(io, mode, 3, 4.0);
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    FIQuery<double, Ray3<double>, Sphere3<double>> query;
    auto r = query(ray, sphere);
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

ORACLE_CASE("IntrRay3Sphere3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto rayOrigin = Pt<3>(io, mode, 3, 4.0);
    auto rayDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    Expose<FIQuery<double, Ray3<double>, Sphere3<double>>> query;
    FIQuery<double, Ray3<double>, Sphere3<double>>::Result r{};
    query.DoQuery(rayOrigin, rayDirection, sphere, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// The tangent branch (discr exactly zero) and the interval clip at exactly
// t = 0 are unreachable with uniform inputs. With sphere center C and radius R
// on the lattice, ray origin P = C + (R+d)*e_i + s*e_j and unit direction e_j,
// the discriminant is exactly R^2 - (R+d)^2 (zero for d = 0) and the roots are
// -s +- root, so s = 0 places a root exactly at the ray origin.
ORACLE_CASE("IntrRay3Sphere3.find.tangent")
{
    Sphere3<double> sphere{};
    sphere.center = io.latticeVec<3>(-3, 3);
    sphere.radius = io.lattice(1, 3);
    int i = io.rawInteger(0, 2);
    int j = (i + 1 + io.rawInteger(0, 1)) % 3;
    double d = static_cast<double>(io.rawInteger(-1, 1));
    double s = static_cast<double>(io.rawInteger(-3, 3));
    Vector3<double> ei{}, ej{};
    ei.MakeUnit(i);
    ej.MakeUnit(j);
    Vector3<double> origin = sphere.center + (sphere.radius + d) * ei + s * ej;
    io.givenVec(origin);
    io.givenVec(ej);
    Ray3<double> ray(origin, ej);
    FIQuery<double, Ray3<double>, Sphere3<double>> query;
    auto r = query(ray, sphere);
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

// ============================= IntrSegment3Sphere3 =======================

namespace
{
    // Upstream's Q(-e) and Q(e) for the segment/sphere TI query. Upstream
    // reports "no intersection" when both are negative, that is when the whole
    // segment lies strictly inside the solid sphere; the port reports an
    // intersection (docs/UPSTREAM-FINDINGS.md IntrSegment3Sphere3.h TIQuery,
    // issue #203). Upstream's other branches agree with the port, so this
    // predicate is exactly the set on which the port deviates.
    bool Seg3Sph3StrictlyInside(Vector3<double> const& p0,
        Vector3<double> const& p1, Sphere3<double> const& sphere)
    {
        Segment3<double> segment(p0, p1);
        Vector3<double> segOrigin{}, segDirection{};
        double segExtent{};
        segment.GetCenteredForm(segOrigin, segDirection, segExtent);
        Vector3<double> diff = segOrigin - sphere.center;
        double a0 = Dot(diff, diff) - sphere.radius * sphere.radius;
        double a1 = Dot(segDirection, diff);
        double discr = a1 * a1 - a0;
        if (discr < 0.0)
        {
            return false;
        }
        double tmp0 = segExtent * segExtent + a0;
        double tmp1 = 2.0 * a1 * segExtent;
        return (tmp0 - tmp1) < 0.0 && (tmp0 + tmp1) < 0.0;
    }

    // Unrecorded draws for a segment and a sphere; the caller records them.
    void RawSegSphere(oracle::Ctx& io, int mode, Vector3<double>& p0,
        Vector3<double>& p1, Vector3<double>& center, double& radius)
    {
        for (int i = 0; i < 3; ++i)
        {
            p0[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        for (int i = 0; i < 3; ++i)
        {
            p1[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        for (int i = 0; i < 3; ++i)
        {
            center[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.25, 4.0));
    }
}

// The generator rejects the segments that lie strictly inside the sphere,
// which is exactly the set on which the port deliberately deviates (see the
// deviation case below).
ORACLE_CASE("IntrSegment3Sphere3.test")
{
    int mode = io.index() % 2;
    Vector3<double> p0{}, p1{}, center{};
    double radius{};
    for (;;)
    {
        RawSegSphere(io, mode, p0, p1, center, radius);
        Sphere3<double> probe{};
        probe.center = center;
        probe.radius = radius;
        if (!Seg3Sph3StrictlyInside(p0, p1, probe)) { break; }
    }
    io.givenVec(p0);
    io.givenVec(p1);
    io.givenVec(center);
    io.given(radius);
    Segment3<double> segment(p0, p1);
    Sphere3<double> sphere{};
    sphere.center = center;
    sphere.radius = radius;
    TIQuery<double, Segment3<double>, Sphere3<double>> query;
    auto r = query(segment, sphere);
    io.outBool(r.intersect);
}

// Deliberate deviation: a segment strictly inside the solid sphere. Both
// endpoint values Q(-e) and Q(e) are negative, so upstream's 'qm * qp <= 0'
// test fails, it falls through to 'qm > 0 && |a1| < e' and reports no
// intersection, disagreeing with its own FI query.
// docs/UPSTREAM-FINDINGS.md IntrSegment3Sphere3.h TIQuery, issue #203.
ORACLE_CASE("IntrSegment3Sphere3.test.containedDeviation")
{
    Sphere3<double> sphere{};
    sphere.center = io.latticeVec<3>(-3, 3);
    sphere.radius = io.lattice(2, 4);
    Vector3<double> p0{}, p1{};
    // Endpoints strictly inside the sphere: |p - C| <= radius/2 < radius.
    for (;;)
    {
        for (int i = 0; i < 3; ++i) { p0[i] = io.raw(-1.0, 1.0); }
        if (Length(p0) <= 1.0 && Length(p0) > 0.0) { break; }
    }
    for (;;)
    {
        for (int i = 0; i < 3; ++i) { p1[i] = io.raw(-1.0, 1.0); }
        if (Length(p1) <= 1.0 && Length(p1) > 0.0) { break; }
    }
    p0 = sphere.center + (0.5 * sphere.radius) * p0;
    p1 = sphere.center + (0.5 * sphere.radius) * p1;
    io.givenVec(p0);
    io.givenVec(p1);
    Segment3<double> segment(p0, p1);
    TIQuery<double, Segment3<double>, Sphere3<double>> query;
    auto r = query(segment, sphere);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3Sphere3.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<3>(io, mode, 3, 4.0);
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    FIQuery<double, Segment3<double>, Sphere3<double>> query;
    auto r = query(segment, sphere);
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

ORACLE_CASE("IntrSegment3Sphere3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto segOrigin = Pt<3>(io, mode, 3, 4.0);
    auto segDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    double segExtent = (mode == 0 ? io.lattice(1, 4) : io.real(0.25, 4.0));
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    Expose<FIQuery<double, Segment3<double>, Sphere3<double>>> query;
    FIQuery<double, Segment3<double>, Sphere3<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, sphere, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// The tangent branch and the clips at exactly +-segExtent: the segment centre
// is C + (R+d)*e_i + s*e_j with unit direction e_j, so the discriminant is
// exactly R^2 - (R+d)^2 and the roots are -s +- root, which the lattice extent
// clips at exact equality.
ORACLE_CASE("IntrSegment3Sphere3.doQuery.fi.tangent")
{
    Sphere3<double> sphere{};
    sphere.center = io.latticeVec<3>(-3, 3);
    sphere.radius = io.lattice(1, 3);
    int i = io.rawInteger(0, 2);
    int j = (i + 1 + io.rawInteger(0, 1)) % 3;
    double d = static_cast<double>(io.rawInteger(-1, 1));
    double s = static_cast<double>(io.rawInteger(-3, 3));
    Vector3<double> ei{}, ej{};
    ei.MakeUnit(i);
    ej.MakeUnit(j);
    Vector3<double> segOrigin = sphere.center + (sphere.radius + d) * ei + s * ej;
    io.givenVec(segOrigin);
    io.givenVec(ej);
    double segExtent = io.lattice(1, 4);
    Expose<FIQuery<double, Segment3<double>, Sphere3<double>>> query;
    FIQuery<double, Segment3<double>, Sphere3<double>>::Result r{};
    query.DoQuery(segOrigin, ej, segExtent, sphere, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// ============================== IntrSphere3Sphere3 =======================

namespace
{
    void OutSphereFI(oracle::Ctx& io,
        FIQuery<double, Sphere3<double>, Sphere3<double>>::Result const& r)
    {
        io.outBool(r.intersect);
        io.outInt(r.type);
        // Upstream leaves 'point' at zero for type 2 and 'circle' at zero for
        // the other types, so both are emitted on every path: the default
        // values are as much a part of the result as the computed ones.
        io.outVec(r.point);
        io.outVec(r.circle.center);
        io.outVec(r.circle.normal);
        io.outReal(r.circle.radius);
    }
}

ORACLE_CASE("IntrSphere3Sphere3.test")
{
    int mode = io.index() % 2;
    auto sphere0 = Sph<3>(io, mode, 3, 4.0);
    auto sphere1 = Sph<3>(io, mode, 3, 4.0);
    TIQuery<double, Sphere3<double>, Sphere3<double>> query;
    auto r = query(sphere0, sphere1);
    io.outBool(r.intersect);
}

// The generator rejects exact internal tangency with r0 <= r1 (upstream's
// type 4), which is exactly the set on which the port deliberately deviates
// (see the deviation case below). Lattice mode reaches it often enough that
// the rejection matters.
ORACLE_CASE("IntrSphere3Sphere3.find")
{
    int mode = io.index() % 2;
    Vector3<double> c0{}, c1{};
    double r0{}, r1{};
    for (;;)
    {
        for (int i = 0; i < 3; ++i)
        {
            c0[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        r0 = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.25, 4.0));
        for (int i = 0; i < 3; ++i)
        {
            c1[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        r1 = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.25, 4.0));
        Vector3<double> diff = c1 - c0;
        double sqrLen = Dot(diff, diff);
        double rDif = r0 - r1;
        if (!(sqrLen == rDif * rDif && rDif <= 0.0)) { break; }
    }
    io.givenVec(c0);
    io.given(r0);
    io.givenVec(c1);
    io.given(r1);
    Sphere3<double> sphere0{}, sphere1{};
    sphere0.center = c0;
    sphere0.radius = r0;
    sphere1.center = c1;
    sphere1.radius = r1;
    FIQuery<double, Sphere3<double>, Sphere3<double>> query;
    auto r = query(sphere0, sphere1);
    OutSphereFI(io, r);
}

// Exact tangency on the lattice. With C1 = C0 + d*e_k the squared distance is
// exactly d^2, so d = r0 + r1 gives external tangency (type 1) and
// d = r0 - r1 > 0 gives internal tangency with sphere1 inside sphere0
// (type 6); the type-4 mirror image is the deviation case below. Mode 2 makes
// the spheres exactly equal, which upstream reports as type 2 with a
// zero-radius circle.
ORACLE_CASE("IntrSphere3Sphere3.find.tangent")
{
    Vector3<double> c0 = io.latticeVec<3>(-3, 3);
    int k = io.rawInteger(0, 2);
    int kind = io.rawInteger(0, 2);
    double r0 = static_cast<double>(io.rawInteger(2, 4));
    double r1 = static_cast<double>(io.rawInteger(1, static_cast<int>(r0) - 1));
    if (kind == 2) { r1 = r0; }
    double d = (kind == 0 ? r0 + r1 : (kind == 1 ? r0 - r1 : 0.0));
    Vector3<double> ek{};
    ek.MakeUnit(k);
    Vector3<double> c1 = c0 + d * ek;
    io.given(r0);
    io.givenVec(c1);
    io.given(r1);
    Sphere3<double> sphere0{}, sphere1{};
    sphere0.center = c0;
    sphere0.radius = r0;
    sphere1.center = c1;
    sphere1.radius = r1;
    FIQuery<double, Sphere3<double>, Sphere3<double>> query;
    auto r = query(sphere0, sphere1);
    OutSphereFI(io, r);
}

// Deliberate deviation: exact internal tangency with sphere0 inside sphere1
// (upstream's type 4). Upstream reports C1 + r1*(C1-C0)/|C1-C0|, the ANTIPODE
// of the contact point; the port reports C1 - r1*(C1-C0)/|C1-C0|.
// docs/UPSTREAM-FINDINGS.md IntrSphere3Sphere3.h FIQuery type-4 branch, issue
// #203.
ORACLE_CASE("IntrSphere3Sphere3.find.internalTangentDeviation")
{
    Vector3<double> c0 = io.latticeVec<3>(-3, 3);
    int k = io.rawInteger(0, 2);
    double r0 = static_cast<double>(io.rawInteger(1, 3));
    double d = static_cast<double>(io.rawInteger(1, 3));
    double r1 = r0 + d;
    Vector3<double> ek{};
    ek.MakeUnit(k);
    Vector3<double> c1 = c0 + d * ek;
    io.given(r0);
    io.givenVec(c1);
    io.given(r1);
    Sphere3<double> sphere0{}, sphere1{};
    sphere0.center = c0;
    sphere0.radius = r0;
    sphere1.center = c1;
    sphere1.radius = r1;
    FIQuery<double, Sphere3<double>, Sphere3<double>> query;
    auto r = query(sphere0, sphere1);
    OutSphereFI(io, r);
}

// =========================== IntrSegment2AlignedBox2 =====================

namespace
{
    void OutSeg2Box2FI(oracle::Ctx& io,
        FIQuery<double, Segment2<double>, AlignedBox2<double>>::Result const& r)
    {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        // Every field is emitted on every path: upstream's degenerate branch
        // fills point[1] without counting it, and the unset values are the
        // default-constructed ones, which are as much a part of the result.
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outReal(r.cdeParameter[0]);
        io.outReal(r.cdeParameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }
}

ORACLE_CASE("IntrSegment2AlignedBox2.test")
{
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 3, 4.0);
    auto box = ABox<2>(io, mode, 3, 4.0);
    TIQuery<double, Segment2<double>, AlignedBox2<double>> query;
    auto r = query(segment, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment2AlignedBox2.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 3, 4.0);
    auto box = ABox<2>(io, mode, 3, 4.0);
    FIQuery<double, Segment2<double>, AlignedBox2<double>> query;
    auto r = query(segment, box);
    OutSeg2Box2FI(io, r);
}

// The degenerate branch: a zero-length segment, which upstream handles with a
// containment test instead of the clipping code. The lattice box makes the
// point land on the boundary often.
ORACLE_CASE("IntrSegment2AlignedBox2.find.degenerate")
{
    Vector2<double> p = io.latticeVec<2>(-3, 3);
    auto box = ABox<2>(io, 0, 3, 4.0);
    Segment2<double> segment(p, p);
    FIQuery<double, Segment2<double>, AlignedBox2<double>> query;
    auto r = query(segment, box);
    OutSeg2Box2FI(io, r);
}

ORACLE_CASE("IntrSegment2AlignedBox2.doQuery.ti")
{
    int mode = io.index() % 2;
    auto segOrigin = Pt<2>(io, mode, 3, 4.0);
    auto segDirection = (mode == 0 ? AxisUnit<2>(io) : io.unit<2>());
    double segExtent = (mode == 0 ? io.lattice(1, 4) : io.real(0.25, 4.0));
    auto boxExtent = Extent<2>(io, mode, 3, 3.0);
    Expose<TIQuery<double, Segment2<double>, AlignedBox2<double>>> query;
    TIQuery<double, Segment2<double>, AlignedBox2<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, boxExtent, r);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment2AlignedBox2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto segOrigin = Pt<2>(io, mode, 3, 4.0);
    auto segDirection = (mode == 0 ? AxisUnit<2>(io) : io.unit<2>());
    double segExtent = (mode == 0 ? io.lattice(1, 4) : io.real(0.25, 4.0));
    auto boxExtent = Extent<2>(io, mode, 3, 3.0);
    Expose<FIQuery<double, Segment2<double>, AlignedBox2<double>>> query;
    FIQuery<double, Segment2<double>, AlignedBox2<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, boxExtent, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// =========================== IntrSegment3AlignedBox3 =====================

ORACLE_CASE("IntrSegment3AlignedBox3.test")
{
    int mode = io.index() % 2;
    auto segment = Sg<3>(io, mode, 3, 4.0);
    auto box = ABox<3>(io, mode, 3, 4.0);
    TIQuery<double, Segment3<double>, AlignedBox3<double>> query;
    auto r = query(segment, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3AlignedBox3.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<3>(io, mode, 3, 4.0);
    auto box = ABox<3>(io, mode, 3, 4.0);
    FIQuery<double, Segment3<double>, AlignedBox3<double>> query;
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

ORACLE_CASE("IntrSegment3AlignedBox3.doQuery.ti")
{
    int mode = io.index() % 2;
    auto segOrigin = Pt<3>(io, mode, 3, 4.0);
    auto segDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    double segExtent = (mode == 0 ? io.lattice(1, 4) : io.real(0.25, 4.0));
    auto boxExtent = Extent<3>(io, mode, 3, 3.0);
    Expose<TIQuery<double, Segment3<double>, AlignedBox3<double>>> query;
    TIQuery<double, Segment3<double>, AlignedBox3<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, boxExtent, r);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3AlignedBox3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto segOrigin = Pt<3>(io, mode, 3, 4.0);
    auto segDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    double segExtent = (mode == 0 ? io.lattice(1, 4) : io.real(0.25, 4.0));
    auto boxExtent = Extent<3>(io, mode, 3, 3.0);
    Expose<FIQuery<double, Segment3<double>, AlignedBox3<double>>> query;
    FIQuery<double, Segment3<double>, AlignedBox3<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, boxExtent, r);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.numIntersections));
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// A segment aimed at an exact face point of a lattice box: the endpoints are
// target + t * direction with integer direction and dyadic t, so the
// Liang-Barsky clips and the [-e,e] interval clip are evaluated at exact
// equality (including the corner and edge cases).
ORACLE_CASE("IntrSegment3AlignedBox3.find.throughPoint")
{
    auto box = ABox<3>(io, 0, 3, 4.0);
    Vector3<double> center{}, extent{};
    box.GetCenteredForm(center, extent);
    Vector3<double> target{};
    for (int i = 0; i < 3; ++i)
    {
        target[i] = center[i]
            + 0.5 * static_cast<double>(io.rawInteger(-2, 2)) * extent[i];
    }
    Vector3<double> direction{};
    for (int i = 0; i < 3; ++i)
    {
        direction[i] = static_cast<double>(io.rawInteger(-2, 2));
    }
    double t0 = 0.5 * static_cast<double>(io.rawInteger(-4, 0));
    double t1 = 0.5 * static_cast<double>(io.rawInteger(0, 4));
    Vector3<double> p0 = target + t0 * direction;
    Vector3<double> p1 = target + t1 * direction;
    io.givenVec(p0);
    io.givenVec(p1);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, AlignedBox3<double>> query;
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

// ============================ IntrSegment2Triangle2 ======================

ORACLE_CASE("IntrSegment2Triangle2.test")
{
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 3, 4.0);
    auto triangle = Tri<2>(io, mode, 3, 4.0);
    TIQuery<double, Segment2<double>, Triangle2<double>> query;
    auto r = query(segment, triangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment2Triangle2.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<2>(io, mode, 3, 4.0);
    auto triangle = Tri<2>(io, mode, 3, 4.0);
    FIQuery<double, Segment2<double>, Triangle2<double>> query;
    auto r = query(segment, triangle);
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

ORACLE_CASE("IntrSegment2Triangle2.doQuery.fi")
{
    int mode = io.index() % 2;
    auto origin = Pt<2>(io, mode, 3, 4.0);
    auto direction = (mode == 0 ? io.latticeDir<2>(-2, 2) : io.unit<2>());
    auto triangle = Tri<2>(io, mode, 3, 4.0);
    Expose<FIQuery<double, Segment2<double>, Triangle2<double>>> query;
    FIQuery<double, Segment2<double>, Triangle2<double>>::Result r{};
    query.DoQuery(origin, direction, triangle, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
}

// Uniform segments hit a triangle a few percent of the time and never on a
// boundary. Here the segment passes through v0 + (k1*E1 + k2*E2)/4 on an
// integer triangle, with p0 = target + t0*direction and p1 = target +
// t1*direction for an integer direction and dyadic t-values, so the clipping
// comparisons of IntrLine2Triangle2::DoQuery and the [0,1] interval clip are
// evaluated at exact equality.
ORACLE_CASE("IntrSegment2Triangle2.find.throughPoint")
{
    Triangle2<double> triangle{};
    triangle.v[0] = io.latticeVec<2>(-4, 4);
    triangle.v[1] = io.latticeVec<2>(-4, 4);
    triangle.v[2] = io.latticeVec<2>(-4, 4);
    double k1 = 0.25 * static_cast<double>(io.rawInteger(0, 5));
    double k2 = 0.25 * static_cast<double>(io.rawInteger(0, 5));
    double t0 = 0.5 * static_cast<double>(io.rawInteger(-4, 1));
    double t1 = 0.5 * static_cast<double>(io.rawInteger(-1, 4));
    Vector2<double> target = triangle.v[0] + k1 * (triangle.v[1] - triangle.v[0])
        + k2 * (triangle.v[2] - triangle.v[0]);
    Vector2<double> direction = io.latticeDir<2>(-2, 2);
    Vector2<double> p0 = target + t0 * direction;
    Vector2<double> p1 = target + t1 * direction;
    io.givenVec(p0);
    io.givenVec(p1);
    Segment2<double> segment(p0, p1);
    FIQuery<double, Segment2<double>, Triangle2<double>> query;
    auto r = query(segment, triangle);
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

// =========================== IntrSegment3Rectangle3 ======================

ORACLE_CASE("IntrSegment3Rectangle3.test")
{
    int mode = io.index() % 2;
    auto segment = Sg<3>(io, mode, 3, 4.0);
    auto rectangle = Rect(io, mode, 3, 3.0);
    TIQuery<double, Segment3<double>, Rectangle3<double>> query;
    auto r = query(segment, rectangle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3Rectangle3.find")
{
    int mode = io.index() % 2;
    auto segment = Sg<3>(io, mode, 3, 4.0);
    auto rectangle = Rect(io, mode, 3, 3.0);
    FIQuery<double, Segment3<double>, Rectangle3<double>> query;
    auto r = query(segment, rectangle);
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

// The rectangle is axis aligned on the lattice and the segment passes through
// C + (a/2)*e0*W0 + (b/2)*e1*W1 with a, b in [-3,3], so the
// |rectCoord| <= extent comparisons are evaluated at exact equality (corners
// and edges included) and the 0 <= t <= 1 test is evaluated at its endpoints.
ORACLE_CASE("IntrSegment3Rectangle3.find.throughPoint")
{
    Rectangle3<double> rectangle{};
    rectangle.center = io.latticeVec<3>(-3, 3);
    auto axis = Frame3(io, 0);
    rectangle.axis[0] = axis[0];
    rectangle.axis[1] = axis[1];
    rectangle.extent = Extent<2>(io, 0, 3, 3.0);
    double a = 0.5 * static_cast<double>(io.rawInteger(-3, 3));
    double b = 0.5 * static_cast<double>(io.rawInteger(-3, 3));
    double t0 = 0.5 * static_cast<double>(io.rawInteger(-4, 1));
    double t1 = 0.5 * static_cast<double>(io.rawInteger(-1, 4));
    Vector3<double> target = rectangle.center
        + (a * rectangle.extent[0]) * rectangle.axis[0]
        + (b * rectangle.extent[1]) * rectangle.axis[1];
    Vector3<double> direction = io.latticeDir<3>(-2, 2);
    Vector3<double> p0 = target + t0 * direction;
    Vector3<double> p1 = target + t1 * direction;
    io.givenVec(p0);
    io.givenVec(p1);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Rectangle3<double>> query;
    auto r = query(segment, rectangle);
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

// ============================ IntrSegment2Segment2 =======================

namespace
{
    int32_t const kInt32Max = 2147483647;

    void OutSeg2Seg2FI(oracle::Ctx& io,
        FIQuery<double, Segment2<double>, Segment2<double>>::Result const& r)
    {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.segment0Parameter[0]);
        io.outReal(r.segment0Parameter[1]);
        io.outReal(r.segment1Parameter[0]);
        io.outReal(r.segment1Parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }

    // The collinear branch of the centered-form FI query with the two centered
    // directions antiparallel, which is exactly the set on which the port
    // deviates: upstream reports segment1Parameter[i] = overlap[i] - t, a
    // parameter measured along segment0's direction, so
    // C1 + segment1Parameter[i] * D1 is not the intersection point. The port
    // negates it there (docs/UPSTREAM-FINDINGS.md IntrSegment2Segment2.h
    // FIQuery, issue #458).
    bool Seg2Seg2AntiparallelCollinear(Segment2<double> const& segment0,
        Segment2<double> const& segment1)
    {
        Vector2<double> seg0Origin{}, seg0Direction{}, seg1Origin{}, seg1Direction{};
        double seg0Extent{}, seg1Extent{};
        segment0.GetCenteredForm(seg0Origin, seg0Direction, seg0Extent);
        segment1.GetCenteredForm(seg1Origin, seg1Direction, seg1Extent);
        FIQuery<double, Line2<double>, Line2<double>> llQuery{};
        Line2<double> line0(seg0Origin, seg0Direction);
        Line2<double> line1(seg1Origin, seg1Direction);
        auto llResult = llQuery(line0, line1);
        if (llResult.numIntersections != kInt32Max)
        {
            return false;
        }
        if (Dot(seg0Direction, seg1Direction) >= 0.0)
        {
            return false;
        }
        Vector2<double> diff = seg1Origin - seg0Origin;
        double t = Dot(seg0Direction, diff);
        std::array<double, 2> interval0 = { -seg0Extent, seg0Extent };
        std::array<double, 2> interval1 = { t - seg1Extent, t + seg1Extent };
        FIQuery<double, std::array<double, 2>, std::array<double, 2>> iiQuery{};
        return iiQuery(interval0, interval1).intersect;
    }
}

ORACLE_CASE("IntrSegment2Segment2.test")
{
    int mode = io.index() % 2;
    auto segment0 = Sg<2>(io, mode, 3, 4.0);
    auto segment1 = Sg<2>(io, mode, 3, 4.0);
    TIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query(segment0, segment1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

ORACLE_CASE("IntrSegment2Segment2.testExact")
{
    int mode = io.index() % 2;
    auto segment0 = Sg<2>(io, mode, 3, 4.0);
    auto segment1 = Sg<2>(io, mode, 3, 4.0);
    TIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query.Exact(segment0, segment1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

// The generator rejects the antiparallel collinear configurations, which are
// exactly the inputs on which the port deliberately deviates (see the
// deviation case below).
ORACLE_CASE("IntrSegment2Segment2.find")
{
    int mode = io.index() % 2;
    Vector2<double> p0{}, p1{}, q0{}, q1{};
    for (;;)
    {
        for (int i = 0; i < 2; ++i)
        {
            p0[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        for (int i = 0; i < 2; ++i)
        {
            p1[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        for (int i = 0; i < 2; ++i)
        {
            q0[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        for (int i = 0; i < 2; ++i)
        {
            q1[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        Segment2<double> probe0(p0, p1);
        Segment2<double> probe1(q0, q1);
        if (!Seg2Seg2AntiparallelCollinear(probe0, probe1)) { break; }
    }
    io.givenVec(p0);
    io.givenVec(p1);
    io.givenVec(q0);
    io.givenVec(q1);
    Segment2<double> segment0(p0, p1);
    Segment2<double> segment1(q0, q1);
    FIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query(segment0, segment1);
    OutSeg2Seg2FI(io, r);
}

ORACLE_CASE("IntrSegment2Segment2.findExact")
{
    int mode = io.index() % 2;
    auto segment0 = Sg<2>(io, mode, 3, 4.0);
    auto segment1 = Sg<2>(io, mode, 3, 4.0);
    FIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query.Exact(segment0, segment1);
    OutSeg2Seg2FI(io, r);
}

namespace
{
    // Two segments on one lattice line. 'flip' reverses the endpoint order of
    // the second segment, which makes the two centered directions
    // antiparallel. The four parameters are drawn so that the segments have
    // positive length; 'requireOverlap' additionally forces the two parameter
    // intervals to overlap.
    void CollinearSegments(oracle::Ctx& io, bool flip, bool requireOverlap,
        Segment2<double>& segment0, Segment2<double>& segment1)
    {
        Vector2<double> base{}, direction{};
        int a{}, b{}, u0{}, u1{};
        for (;;)
        {
            for (int i = 0; i < 2; ++i)
            {
                base[i] = static_cast<double>(io.rawInteger(-3, 3));
            }
            bool zero = true;
            do
            {
                zero = true;
                for (int i = 0; i < 2; ++i)
                {
                    direction[i] = static_cast<double>(io.rawInteger(-2, 2));
                    zero = zero && direction[i] == 0.0;
                }
            } while (zero);
            a = io.rawInteger(-3, 2);
            b = a + io.rawInteger(1, 3);
            u0 = io.rawInteger(-3, 2);
            u1 = u0 + io.rawInteger(1, 3);
            if (!requireOverlap) { break; }
            if (std::max(a, u0) <= std::min(b, u1)) { break; }
        }
        Vector2<double> p0 = base + static_cast<double>(a) * direction;
        Vector2<double> p1 = base + static_cast<double>(b) * direction;
        Vector2<double> q0 = base + static_cast<double>(u0) * direction;
        Vector2<double> q1 = base + static_cast<double>(u1) * direction;
        io.givenVec(p0);
        io.givenVec(p1);
        io.givenVec(flip ? q1 : q0);
        io.givenVec(flip ? q0 : q1);
        segment0 = Segment2<double>(p0, p1);
        segment1 = (flip ? Segment2<double>(q1, q0) : Segment2<double>(q0, q1));
    }
}

// Both segments on one lattice line with the same orientation, which reaches
// the "same line" branch of IntrLine2Line2 exactly and clips two overlapping
// intervals at exact equality (including the single-point touch).
ORACLE_CASE("IntrSegment2Segment2.find.collinear")
{
    Segment2<double> segment0{}, segment1{};
    CollinearSegments(io, false, false, segment0, segment1);
    FIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query(segment0, segment1);
    OutSeg2Seg2FI(io, r);
}

ORACLE_CASE("IntrSegment2Segment2.findExact.collinear")
{
    int flip = io.integer(0, 1);
    Segment2<double> segment0{}, segment1{};
    CollinearSegments(io, flip != 0, false, segment0, segment1);
    FIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query.Exact(segment0, segment1);
    OutSeg2Seg2FI(io, r);
}

ORACLE_CASE("IntrSegment2Segment2.test.collinear")
{
    int flip = io.integer(0, 1);
    Segment2<double> segment0{}, segment1{};
    CollinearSegments(io, flip != 0, false, segment0, segment1);
    TIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query(segment0, segment1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

ORACLE_CASE("IntrSegment2Segment2.testExact.collinear")
{
    int flip = io.integer(0, 1);
    Segment2<double> segment0{}, segment1{};
    CollinearSegments(io, flip != 0, false, segment0, segment1);
    TIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query.Exact(segment0, segment1);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
}

// Deliberate deviation: two overlapping collinear segments whose centered
// directions are antiparallel. Upstream's segment1Parameter has the wrong
// sign, so C1 + segment1Parameter[i] * D1 is not the reported point.
// docs/UPSTREAM-FINDINGS.md IntrSegment2Segment2.h FIQuery, issue #458.
ORACLE_CASE("IntrSegment2Segment2.find.antiparallelDeviation")
{
    Segment2<double> segment0{}, segment1{};
    CollinearSegments(io, true, true, segment0, segment1);
    FIQuery<double, Segment2<double>, Segment2<double>> query;
    auto r = query(segment0, segment1);
    OutSeg2Seg2FI(io, r);
}

// ============================ IntrSphere3Triangle3 =======================

namespace
{
    void OutSphereTriangle(oracle::Ctx& io,
        FIQuery<double, Sphere3<double>, Triangle3<double>>::Result const& r)
    {
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }
}

// The dynamic sphere-versus-triangle query. Uniform velocities reach the
// initial-overlap and no-contact branches; the aimed variant below reaches the
// face, half-cylinder and sphere-wedge branches of the sphere-swept volume.
ORACLE_CASE("IntrSphere3Triangle3.find")
{
    int mode = io.index() % 2;
    auto sphere = Sph<3>(io, mode, 3, 3.0);
    auto sphereVelocity = (mode == 0 ? io.latticeVec<3>(-2, 2) : io.vec<3>(-2.0, 2.0));
    auto triangle = Tri<3>(io, mode, 3, 3.0);
    auto triangleVelocity = (mode == 0 ? io.latticeVec<3>(-2, 2) : io.vec<3>(-2.0, 2.0));
    FIQuery<double, Sphere3<double>, Triangle3<double>> query;
    auto r = query(sphere, sphereVelocity, triangle, triangleVelocity);
    OutSphereTriangle(io, r);
}

// The sphere starts away from the triangle and moves towards the point
// v0 + k1*E1 + k2*E2 with dyadic k1, k2 in [-1/2, 3/2], so the target is
// inside the triangle, on an edge, on a vertex or just outside it; the
// relative velocity then meets the triangular face, a half cylinder or a
// sphere wedge of the sphere-swept volume. The triangle velocity is nonzero on
// half the records, so the relative-velocity path is exercised too.
ORACLE_CASE("IntrSphere3Triangle3.find.aimed")
{
    Triangle3<double> triangle{};
    triangle.v[0] = io.latticeVec<3>(-3, 3);
    triangle.v[1] = io.latticeVec<3>(-3, 3);
    triangle.v[2] = io.latticeVec<3>(-3, 3);
    double k1 = 0.25 * static_cast<double>(io.rawInteger(-2, 6));
    double k2 = 0.25 * static_cast<double>(io.rawInteger(-2, 6));
    Vector3<double> target = triangle.v[0] + k1 * (triangle.v[1] - triangle.v[0])
        + k2 * (triangle.v[2] - triangle.v[0]);
    Vector3<double> offset{};
    for (int i = 0; i < 3; ++i) { offset[i] = static_cast<double>(io.rawInteger(-3, 3)); }
    double radius = io.given(static_cast<double>(io.rawInteger(1, 3)));
    Vector3<double> center = io.givenVec(target + offset);
    Vector3<double> triangleVelocity = io.latticeVec<3>(-1, 1);
    // The sphere moves towards the target with the given relative speed.
    double speed = 0.25 * static_cast<double>(io.rawInteger(1, 8));
    Vector3<double> sphereVelocity =
        io.givenVec(triangleVelocity - speed * offset);
    Sphere3<double> sphere{};
    sphere.center = center;
    sphere.radius = radius;
    FIQuery<double, Sphere3<double>, Triangle3<double>> query;
    auto r = query(sphere, sphereVelocity, triangle, triangleVelocity);
    OutSphereTriangle(io, r);
}

// =========================== IntrAlignedBox2Circle2 ======================

ORACLE_CASE("IntrAlignedBox2Circle2.test")
{
    int mode = io.index() % 2;
    auto box = ABox<2>(io, mode, 3, 4.0);
    auto circle = Sph<2>(io, mode, 3, 4.0);
    TIQuery<double, AlignedBox2<double>, Circle2<double>> query;
    auto r = query(box, circle);
    io.outBool(r.intersect);
}

namespace
{
    void OutBox2Circle2(oracle::Ctx& io,
        FIQuery<double, AlignedBox2<double>, Circle2<double>>::Result const& r)
    {
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }
}

// The dynamic box-versus-circle query. The lattice mode makes the circle
// centre land exactly on the box edges and corners, which reaches the
// overlap branches (InteriorOverlap, EdgeOverlap, VertexOverlap) at exact
// equality; the uniform mode reaches the separated branches.
ORACLE_CASE("IntrAlignedBox2Circle2.find")
{
    int mode = io.index() % 2;
    auto box = ABox<2>(io, mode, 3, 4.0);
    auto boxVelocity = (mode == 0 ? io.latticeVec<2>(-2, 2) : io.vec<2>(-2.0, 2.0));
    auto circle = Sph<2>(io, mode, 3, 4.0);
    auto circleVelocity = (mode == 0 ? io.latticeVec<2>(-2, 2) : io.vec<2>(-2.0, 2.0));
    FIQuery<double, AlignedBox2<double>, Circle2<double>> query;
    auto r = query(box, boxVelocity, circle, circleVelocity);
    OutBox2Circle2(io, r);
}

// The circle starts outside the box and moves towards a point of the box
// boundary, which reaches the IntersectsEdge and IntersectsVertex branches of
// the unbounded cases (and the rounded-corner quadratic) far more often than
// uniform velocities do.
ORACLE_CASE("IntrAlignedBox2Circle2.find.aimed")
{
    auto box = ABox<2>(io, 0, 3, 4.0);
    Vector2<double> boxCenter{}, extent{};
    box.GetCenteredForm(boxCenter, extent);
    Vector2<double> target{};
    for (int i = 0; i < 2; ++i)
    {
        target[i] = boxCenter[i]
            + 0.5 * static_cast<double>(io.rawInteger(-3, 3)) * extent[i];
    }
    Vector2<double> offset{};
    for (int i = 0; i < 2; ++i)
    {
        offset[i] = static_cast<double>(io.rawInteger(-4, 4));
    }
    double radius = io.given(static_cast<double>(io.rawInteger(1, 3)));
    Vector2<double> center = io.givenVec(target + offset);
    Vector2<double> boxVelocity = io.latticeVec<2>(-1, 1);
    double speed = 0.25 * static_cast<double>(io.rawInteger(1, 8));
    Vector2<double> circleVelocity = io.givenVec(boxVelocity - speed * offset);
    Circle2<double> circle{};
    circle.center = center;
    circle.radius = radius;
    FIQuery<double, AlignedBox2<double>, Circle2<double>> query;
    auto r = query(box, boxVelocity, circle, circleVelocity);
    OutBox2Circle2(io, r);
}

namespace
{
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
}

// ========================== IntrAlignedBox3Cylinder3 =====================

namespace
{
    // Unrecorded draws for a box and a finite cylinder, in the layout of
    // ABox<3> followed by Cyl.
    void RawBoxCylinder(oracle::Ctx& io, int mode, Vector3<double>& boxMin,
        Vector3<double>& boxMax, Vector3<double>& origin,
        Vector3<double>& direction, double& radius, double& height)
    {
        for (int i = 0; i < 3; ++i)
        {
            double a{}, b{};
            if (mode == 0)
            {
                a = static_cast<double>(io.rawInteger(-3, 3));
                b = static_cast<double>(io.rawInteger(-3, 3));
            }
            else
            {
                a = io.raw(-3.0, 3.0);
                b = io.raw(-3.0, 3.0);
            }
            boxMin[i] = std::min(a, b);
            boxMax[i] = std::max(a, b);
        }
        for (int i = 0; i < 3; ++i)
        {
            origin[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                : io.raw(-2.0, 2.0));
        }
        if (mode == 0)
        {
            int k = io.rawInteger(0, 5);
            direction.MakeZero();
            direction[k % 3] = (k < 3 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { direction[i] = io.raw(-1.0, 1.0); }
                len = Length(direction);
            } while (len < 0.1 || len > 1.0);
            Normalize(direction);
        }
        radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 2))
            : io.raw(0.25, 2.0));
        height = (mode == 0 ? static_cast<double>(io.rawInteger(1, 4))
            : io.raw(0.5, 4.0));
    }

    // 'requireDeviation' selects the inputs on which the corrected
    // IntrCanonicalBox3Cylinder3 answer differs from upstream's (the deviation
    // case) or agrees with it (the main case). The aligned-box query is a thin
    // wrapper around the canonical-box query, so it inherits the (U1,-D) sign
    // typo of DoQueryNoZeros and the port's correction of it.
    void AlignedBoxCylinderCase(oracle::Ctx& io, int mode, bool requireDeviation)
    {
        Vector3<double> boxMin{}, boxMax{}, origin{}, direction{};
        double radius{}, height{};
        for (;;)
        {
            RawBoxCylinder(io, mode, boxMin, boxMax, origin, direction, radius, height);
            AlignedBox3<double> probeBox(boxMin, boxMax);
            Line3<double> probeAxis(origin, direction);
            Cylinder3<double> probeCylinder(probeAxis, radius, height);
            TIQuery<double, AlignedBox3<double>, Cylinder3<double>> probeQuery;
            bool upstream = probeQuery(probeBox, probeCylinder).intersect;
            // The canonical-box query is applied to the translated cylinder.
            Vector3<double> boxCenter{}, boxExtent{};
            probeBox.GetCenteredForm(boxCenter, boxExtent);
            CanonicalBox3<double> cbox(boxExtent);
            Cylinder3<double> translated = probeCylinder;
            translated.axis.origin -= boxCenter;
            if ((upstream != CbcFixed(cbox, translated, upstream)) == requireDeviation)
            {
                break;
            }
        }
        io.givenVec(boxMin);
        io.givenVec(boxMax);
        io.givenVec(origin);
        io.givenVec(direction);
        io.given(radius);
        io.given(height);
        AlignedBox3<double> box(boxMin, boxMax);
        Line3<double> axis(origin, direction);
        Cylinder3<double> cylinder(axis, radius, height);
        TIQuery<double, AlignedBox3<double>, Cylinder3<double>> query;
        auto r = query(box, cylinder);
        io.outBool(r.intersect);
    }
}

// The generator rejects the inputs on which the (U1,-D) sign typo of
// IntrCanonicalBox3Cylinder3::DoQueryNoZeros changes the answer, which is
// exactly the set on which the port deliberately deviates (see the deviation
// case below).
ORACLE_CASE("IntrAlignedBox3Cylinder3.test")
{
    AlignedBoxCylinderCase(io, io.index() % 2, false);
}

// Deliberate deviation: the aligned-box query inherits the (U1,-D) sign typo
// of IntrCanonicalBox3Cylinder3::DoQueryNoZeros, which puts a segment endpoint
// on the wrong box edge. docs/UPSTREAM-FINDINGS.md
// IntrCanonicalBox3Cylinder3.h, issue #197.
ORACLE_CASE("IntrAlignedBox3Cylinder3.test.edgeTypoDeviation")
{
    AlignedBoxCylinderCase(io, 1, true);
}

// Upstream asserts that the cylinder is finite, so both sides throw. The
// record has no outputs.
ORACLE_CASE("IntrAlignedBox3Cylinder3.test.infinite")
{
    auto box = ABox<3>(io, 1, 3, 3.0);
    Vector3<double> origin = io.vec<3>(-2.0, 2.0);
    Vector3<double> direction = io.unit<3>();
    double radius = io.real(0.25, 2.0);
    double height = io.given(-1.0);
    Line3<double> axis(origin, direction);
    Cylinder3<double> cylinder(axis, radius, height);
    TIQuery<double, AlignedBox3<double>, Cylinder3<double>> query;
    auto r = query(box, cylinder);
    io.outBool(r.intersect);
}


// ======================= IntrConvexPolygonHyperplane =====================

namespace
{
    // The twelve lattice points on the circle of radius 5 centred at the
    // origin, in counterclockwise order. Any cyclically increasing subset of
    // them is a convex polygon with exactly representable vertices.
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

    // A convex polygon in 2D with 'n' vertices in counterclockwise order.
    // Mode 0 uses the lattice circle (exact); mode 1 uses increasing angles,
    // with the cos/sin applied to unrecorded draws so that only the resulting
    // vertices are recorded.
    std::vector<Vector2<double>> Poly2(oracle::Ctx& io, int mode, int n)
    {
        std::vector<Vector2<double>> polygon{};
        if (mode == 0)
        {
            Vector2<double> center{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            int k = io.rawInteger(0, 11);
            for (int i = 0; i < n; ++i)
            {
                polygon.push_back(center + CirclePoint5(k % 12));
                k += io.rawInteger(1, 2);
            }
        }
        else
        {
            Vector2<double> center{ io.raw(-3.0, 3.0), io.raw(-3.0, 3.0) };
            double radius = io.raw(1.0, 4.0);
            double angle = io.raw(0.0, 6.283185307179586);
            for (int i = 0; i < n; ++i)
            {
                polygon.push_back(center + radius * Vector2<double>{
                    std::cos(angle), std::sin(angle) });
                angle += io.raw(0.3, 1.0);
            }
        }
        for (auto const& p : polygon) { io.givenVec(p); }
        return polygon;
    }

    // A convex polygon in 3D. Mode 0 puts it in a coordinate plane on the
    // lattice circle, so the vertices and the plane normal are exact; mode 1
    // uses a general orthonormal frame.
    std::vector<Vector3<double>> Poly3(oracle::Ctx& io, int mode, int n,
        Vector3<double>& planeNormal)
    {
        std::vector<Vector3<double>> polygon{};
        if (mode == 0)
        {
            int axis = io.rawInteger(0, 2);
            int i0 = (axis + 1) % 3, i1 = (axis + 2) % 3;
            Vector3<double> center{};
            for (int i = 0; i < 3; ++i)
            {
                center[i] = static_cast<double>(io.rawInteger(-3, 3));
            }
            int k = io.rawInteger(0, 11);
            for (int i = 0; i < n; ++i)
            {
                Vector2<double> q = CirclePoint5(k % 12);
                Vector3<double> v = center;
                v[i0] += q[0];
                v[i1] += q[1];
                polygon.push_back(v);
                k += io.rawInteger(1, 2);
            }
            planeNormal.MakeZero();
            planeNormal[axis] = 1.0;
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
            Vector3<double> center{ io.raw(-3.0, 3.0), io.raw(-3.0, 3.0),
                io.raw(-3.0, 3.0) };
            double radius = io.raw(1.0, 4.0);
            double angle = io.raw(0.0, 6.283185307179586);
            for (int i = 0; i < n; ++i)
            {
                polygon.push_back(center
                    + (radius * std::cos(angle)) * basis[1]
                    + (radius * std::sin(angle)) * basis[2]);
                angle += io.raw(0.3, 1.0);
            }
            planeNormal = basis[0];
        }
        for (auto const& p : polygon) { io.givenVec(p); }
        return polygon;
    }

    // Unrecorded unit normals: the caller records the normal exactly once,
    // after the kind-2 branch has possibly replaced it by an edge normal.
    Vector2<double> RawUnit2(oracle::Ctx& io, int mode)
    {
        if (mode == 0)
        {
            int k = io.rawInteger(0, 3);
            Vector2<double> v{};
            v.MakeZero();
            v[k % 2] = (k < 2 ? 1.0 : -1.0);
            return v;
        }
        Vector2<double> v{};
        double len = 0.0;
        do
        {
            for (int i = 0; i < 2; ++i) { v[i] = io.raw(-1.0, 1.0); }
            len = Length(v);
        } while (len < 0.1 || len > 1.0);
        Normalize(v);
        return v;
    }

    Vector3<double> RawUnit3(oracle::Ctx& io, int mode)
    {
        if (mode == 0)
        {
            int k = io.rawInteger(0, 5);
            Vector3<double> v{};
            v.MakeZero();
            v[k % 3] = (k < 3 ? 1.0 : -1.0);
            return v;
        }
        Vector3<double> v{};
        double len = 0.0;
        do
        {
            for (int i = 0; i < 3; ++i) { v[i] = io.raw(-1.0, 1.0); }
            len = Length(v);
        } while (len < 0.1 || len > 1.0);
        Normalize(v);
        return v;
    }

    // The hyperplane of each case is built inline: 'kind' selects the
    // configuration the query reaches. Kinds 1, 2 and 3 build the constant
    // from the polygon itself, so the heights of the touched vertices are
    // exactly zero (in lattice mode also for the whole edge of kind 2, whose
    // normal is the integer perpendicular of the edge rather than a unit
    // vector: this query compares only the signs of Dot(N,X) - c, which a
    // positive scaling of (N,c) leaves unchanged).
}

ORACLE_CASE("IntrConvexPolygonHyperplane.test.2d")
{
    int mode = io.index() % 2;
    int n = io.integer(2, 6);
    auto polygon = Poly2(io, mode, n);
    int kind = io.integer(0, 2);
    Vector2<double> normal{};
    if (kind == 2 && n >= 2)
    {
        // The perpendicular of one polygon edge: the edge lies on the plane.
        size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
        size_t jp1 = (j + 1) % static_cast<size_t>(n);
        normal = Perp(polygon[jp1] - polygon[j]);
    }
    else
    {
        normal = RawUnit2(io, mode);
    }
    Hyperplane<2, double> hyperplane{};
    hyperplane.normal = io.givenVec(normal);
    if (kind == 0)
    {
        hyperplane.constant = io.real(-6.0, 6.0);
    }
    else
    {
        size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
        hyperplane.constant = io.given(Dot(normal, polygon[j]));
    }
    TIQuery<double, std::vector<Vector2<double>>, Hyperplane<2, double>> query;
    auto r = query(polygon, hyperplane);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.configuration));
}

namespace
{
    template <int N>
    void OutPolyPlaneFI(oracle::Ctx& io,
        typename FIQuery<double, std::vector<Vector<N, double>>,
            Hyperplane<N, double>>::Result const& r)
    {
        io.outInt(static_cast<int32_t>(r.configuration));
        io.outInt(static_cast<int32_t>(r.intersection.size()));
        for (auto const& p : r.intersection) { io.outVec(p); }
        io.outInt(static_cast<int32_t>(r.positivePolygon.size()));
        for (auto const& p : r.positivePolygon) { io.outVec(p); }
        io.outInt(static_cast<int32_t>(r.negativePolygon.size()));
        for (auto const& p : r.negativePolygon) { io.outVec(p); }
    }
}

ORACLE_CASE("IntrConvexPolygonHyperplane.find.2d")
{
    int mode = io.index() % 2;
    int n = io.integer(2, 6);
    auto polygon = Poly2(io, mode, n);
    int kind = io.integer(0, 2);
    Vector2<double> normal{};
    if (kind == 2 && n >= 2)
    {
        size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
        size_t jp1 = (j + 1) % static_cast<size_t>(n);
        normal = Perp(polygon[jp1] - polygon[j]);
    }
    else
    {
        normal = RawUnit2(io, mode);
    }
    Hyperplane<2, double> hyperplane{};
    hyperplane.normal = io.givenVec(normal);
    if (kind == 0)
    {
        hyperplane.constant = io.real(-6.0, 6.0);
    }
    else
    {
        size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
        hyperplane.constant = io.given(Dot(normal, polygon[j]));
    }
    FIQuery<double, std::vector<Vector2<double>>, Hyperplane<2, double>> query;
    auto r = query(polygon, hyperplane);
    OutPolyPlaneFI<2>(io, r);
}

ORACLE_CASE("IntrConvexPolygonHyperplane.test.3d")
{
    int mode = io.index() % 2;
    int n = io.integer(2, 6);
    Vector3<double> planeNormal{};
    auto polygon = Poly3(io, mode, n, planeNormal);
    int kind = io.integer(0, 3);
    Vector3<double> normal{};
    if (kind == 3)
    {
        // The plane of the polygon: every height is zero (CONTAINED).
        normal = planeNormal;
    }
    else if (kind == 2 && n >= 2)
    {
        // A plane containing one polygon edge and the polygon normal.
        size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
        size_t jp1 = (j + 1) % static_cast<size_t>(n);
        normal = Cross(polygon[jp1] - polygon[j], planeNormal);
    }
    else
    {
        normal = RawUnit3(io, mode);
    }
    Hyperplane<3, double> hyperplane{};
    hyperplane.normal = io.givenVec(normal);
    if (kind == 0)
    {
        hyperplane.constant = io.real(-6.0, 6.0);
    }
    else
    {
        size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
        hyperplane.constant = io.given(Dot(normal, polygon[j]));
    }
    TIQuery<double, std::vector<Vector3<double>>, Hyperplane<3, double>> query;
    auto r = query(polygon, hyperplane);
    io.outBool(r.intersect);
    io.outInt(static_cast<int32_t>(r.configuration));
}

ORACLE_CASE("IntrConvexPolygonHyperplane.find.3d")
{
    int mode = io.index() % 2;
    int n = io.integer(2, 6);
    Vector3<double> planeNormal{};
    auto polygon = Poly3(io, mode, n, planeNormal);
    int kind = io.integer(0, 3);
    Vector3<double> normal{};
    if (kind == 3)
    {
        normal = planeNormal;
    }
    else if (kind == 2 && n >= 2)
    {
        size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
        size_t jp1 = (j + 1) % static_cast<size_t>(n);
        normal = Cross(polygon[jp1] - polygon[j], planeNormal);
    }
    else
    {
        normal = RawUnit3(io, mode);
    }
    Hyperplane<3, double> hyperplane{};
    hyperplane.normal = io.givenVec(normal);
    if (kind == 0)
    {
        hyperplane.constant = io.real(-6.0, 6.0);
    }
    else
    {
        size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
        hyperplane.constant = io.given(Dot(normal, polygon[j]));
    }
    FIQuery<double, std::vector<Vector3<double>>, Hyperplane<3, double>> query;
    auto r = query(polygon, hyperplane);
    OutPolyPlaneFI<3>(io, r);
}

// =========================== IntrAlignedBox3Sphere3 ======================

namespace
{
    using Ab3S3Result = FIQuery<double, AlignedBox3<double>, Sphere3<double>>::Result;

    // ---- the port's corrected DoQuery ---------------------------------------
    // Upstream probes at most one rounded piece of the Minkowski sum per case
    // and accepts the first probe that reports a contact, so it reports a
    // contact later than the true first contact, or misses one entirely
    // (docs/UPSTREAM-FINDINGS.md IntrAlignedBox3Sphere3.h DoQueryRayRoundedFace
    // / DoQuery, issues #458 and #465). The port probes every candidate piece
    // and keeps the earliest contact. Everything below is upstream's code
    // except the three marked places, which follow the port; the generator
    // needs the corrected answer to tell the sound inputs from the defective
    // ones.

    void FixRoundedVertex(Vector3<double> const& K, double radius,
        Vector3<double> const& delta, Vector3<double> const& V,
        Ab3S3Result& result)
    {
        double a1 = Dot(V, delta);
        if (a1 < 0.0)
        {
            double a0 = Dot(delta, delta) - radius * radius;
            double a2 = Dot(V, V);
            double adiscr = a1 * a1 - a2 * a0;
            if (adiscr >= 0.0)
            {
                result.intersectionType = 1;
                result.contactTime = -(a1 + std::sqrt(adiscr)) / a2;
                result.contactPoint = K;
            }
        }
    }

    void FixRoundedEdge(int32_t i0, int32_t i1, int32_t i2,
        Vector3<double> const& K, Vector3<double> const& C, double radius,
        Vector3<double> const& delta, Vector3<double> const& V,
        Ab3S3Result& result)
    {
        double b1 = V[i0] * delta[i0] + V[i1] * delta[i1];
        if (b1 < 0.0)
        {
            double b0 = delta[i0] * delta[i0] + delta[i1] * delta[i1] - radius * radius;
            double b2 = V[i0] * V[i0] + V[i1] * V[i1];
            double bdiscr = b1 * b1 - b2 * b0;
            if (bdiscr >= 0.0)
            {
                double tmax = -(b1 + std::sqrt(bdiscr)) / b2;
                double p2 = C[i2] + tmax * V[i2];
                if (-K[i2] <= p2)
                {
                    if (p2 <= K[i2])
                    {
                        result.intersectionType = 1;
                        result.contactTime = tmax;
                        result.contactPoint[i0] = K[i0];
                        result.contactPoint[i1] = K[i1];
                        result.contactPoint[i2] = p2;
                    }
                    else
                    {
                        FixRoundedVertex(K, radius, delta, V, result);
                    }
                }
                else
                {
                    Vector3<double> otherK{}, otherDelta{};
                    otherK[i0] = K[i0];
                    otherK[i1] = K[i1];
                    otherK[i2] = -K[i2];
                    otherDelta[i0] = C[i0] - otherK[i0];
                    otherDelta[i1] = C[i1] - otherK[i1];
                    otherDelta[i2] = C[i2] - otherK[i2];
                    FixRoundedVertex(otherK, radius, otherDelta, V, result);
                }
            }
        }
    }

    // Port addition: keep the earliest of the candidate contacts.
    void FixKeepEarliest(Ab3S3Result& result, Ab3S3Result const& candidate)
    {
        if (candidate.intersectionType == 0)
        {
            return;
        }
        if (result.intersectionType == 0 || candidate.contactTime < result.contactTime)
        {
            result.intersectionType = candidate.intersectionType;
            result.contactTime = candidate.contactTime;
            result.contactPoint = candidate.contactPoint;
        }
    }

    // Port addition: probe the rounded edge parallel to axis i2 whose cylinder
    // axis passes through (sign0 * K[i0], sign1 * K[i1]).
    void FixProbeRoundedEdge(int32_t i0, int32_t i1, int32_t i2,
        Vector3<double> const& K, Vector3<double> const& C, double radius,
        Vector3<double> const& V, double sign0, double sign1,
        Ab3S3Result& result)
    {
        Vector3<double> edgeK{}, edgeDelta{};
        edgeK[i0] = sign0 * K[i0];
        edgeK[i1] = sign1 * K[i1];
        edgeK[i2] = K[i2];
        edgeDelta[i0] = C[i0] - edgeK[i0];
        edgeDelta[i1] = C[i1] - edgeK[i1];
        edgeDelta[i2] = C[i2] - edgeK[i2];
        double b0 = edgeDelta[i0] * edgeDelta[i0] + edgeDelta[i1] * edgeDelta[i1]
            - radius * radius;
        if (b0 <= 0.0)
        {
            return;
        }
        Ab3S3Result candidate{};
        FixRoundedEdge(i0, i1, i2, edgeK, C, radius, edgeDelta, V, candidate);
        FixKeepEarliest(result, candidate);
    }

    void FixRoundedFace(int32_t i0, int32_t i1, int32_t i2,
        Vector3<double> const& K, Vector3<double> const& C, double radius,
        Vector3<double> const& delta, Vector3<double> const& V,
        Ab3S3Result& result)
    {
        double tmax = (radius - delta[i0]) / V[i0];
        double p1 = C[i1] + tmax * V[i1];
        double p2 = C[i2] + tmax * V[i2];
        double s1 = (p1 < -K[i1] ? -1.0 : (p1 > K[i1] ? 1.0 : 0.0));
        double s2 = (p2 < -K[i2] ? -1.0 : (p2 > K[i2] ? 1.0 : 0.0));
        if (s1 == 0.0 && s2 == 0.0)
        {
            result.intersectionType = 1;
            result.contactTime = tmax;
            result.contactPoint[i0] = K[i0];
            result.contactPoint[i1] = p1;
            result.contactPoint[i2] = p2;
            return;
        }
        // Port change: probe every candidate piece, keep the earliest.
        if (s1 != 0.0)
        {
            FixProbeRoundedEdge(i0, i1, i2, K, C, radius, V, 1.0, s1, result);
        }
        if (s2 != 0.0)
        {
            FixProbeRoundedEdge(i2, i0, i1, K, C, radius, V, s2, 1.0, result);
        }
        if (s1 != 0.0 && s2 != 0.0)
        {
            FixProbeRoundedEdge(i1, i2, i0, K, C, radius, V, s1, s2, result);
        }
    }

    void FixDoQuery(Vector3<double> const& K, Vector3<double> const& inC,
        double radius, Vector3<double> const& inV, Ab3S3Result& result)
    {
        Vector3<double> C = inC, V = inV;
        std::array<double, 3> sign = { 0.0, 0.0, 0.0 };
        for (int32_t i = 0; i < 3; ++i)
        {
            if (C[i] >= 0.0)
            {
                sign[i] = 1.0;
            }
            else
            {
                C[i] = -C[i];
                V[i] = -V[i];
                sign[i] = -1.0;
            }
        }

        Vector3<double> delta = C - K;
        auto interiorOverlap = [&result](Vector3<double> const& c)
        {
            result.intersectionType = -1;
            result.contactTime = 0.0;
            result.contactPoint = c;
        };
        auto vertexOverlap = [&result](Vector3<double> const& k, double r,
            Vector3<double> const& d)
        {
            result.intersectionType = (Dot(d, d) < r * r ? -1 : 1);
            result.contactTime = 0.0;
            result.contactPoint = k;
        };
        auto edgeOverlap = [&result](int32_t i0, int32_t i1, int32_t i2,
            Vector3<double> const& k, Vector3<double> const& c, double r,
            Vector3<double> const& d)
        {
            result.intersectionType =
                (d[i0] * d[i0] + d[i1] * d[i1] < r * r ? -1 : 1);
            result.contactTime = 0.0;
            result.contactPoint[i0] = k[i0];
            result.contactPoint[i1] = k[i1];
            result.contactPoint[i2] = c[i2];
        };
        auto faceOverlap = [&result](int32_t i0, int32_t i1, int32_t i2,
            Vector3<double> const& k, Vector3<double> const& c, double r,
            Vector3<double> const& d)
        {
            result.intersectionType = (d[i0] < r ? -1 : 1);
            result.contactTime = 0.0;
            result.contactPoint[i0] = k[i0];
            result.contactPoint[i1] = c[i1];
            result.contactPoint[i2] = c[i2];
        };
        auto vertexSeparated = [&](Vector3<double> const& k, double r,
            Vector3<double> const& d, Vector3<double> const& v)
        {
            if (v[0] < 0.0 || v[1] < 0.0 || v[2] < 0.0)
            {
                // Port change: the three rounded edges that meet the vertex
                // are probed as well and the earliest contact is kept.
                FixRoundedVertex(k, r, d, v, result);
                FixProbeRoundedEdge(0, 1, 2, k, C, r, v, 1.0, 1.0, result);
                FixProbeRoundedEdge(1, 2, 0, k, C, r, v, 1.0, 1.0, result);
                FixProbeRoundedEdge(2, 0, 1, k, C, r, v, 1.0, 1.0, result);
            }
        };
        auto edgeSeparated = [&](int32_t i0, int32_t i1, int32_t i2,
            Vector3<double> const& k, Vector3<double> const& c, double r,
            Vector3<double> const& d, Vector3<double> const& v)
        {
            if (v[i0] < 0.0 || v[i1] < 0.0)
            {
                FixRoundedEdge(i0, i1, i2, k, c, r, d, v, result);
            }
        };
        auto vertexUnbounded = [&](Vector3<double> const& k,
            Vector3<double> const& c, double r, Vector3<double> const& d,
            Vector3<double> const& v)
        {
            if (v[0] < 0.0 && v[1] < 0.0 && v[2] < 0.0)
            {
                double tmax = (r - d[0]) / v[0];
                int32_t j0 = 0;
                double temp = (r - d[1]) / v[1];
                if (temp > tmax) { tmax = temp; j0 = 1; }
                temp = (r - d[2]) / v[2];
                if (temp > tmax) { tmax = temp; j0 = 2; }
                int32_t j1 = (j0 + 1) % 3;
                int32_t j2 = (j1 + 1) % 3;
                FixRoundedFace(j0, j1, j2, k, c, r, d, v, result);
            }
        };
        auto edgeUnbounded = [&](int32_t i0, int32_t i1, Vector3<double> const& k,
            Vector3<double> const& c, double r, Vector3<double> const& d,
            Vector3<double> const& v)
        {
            if (v[i0] < 0.0 && v[i1] < 0.0)
            {
                double tmax = (r - d[i0]) / v[i0];
                int32_t j0 = i0;
                double temp = (r - d[i1]) / v[i1];
                if (temp > tmax) { tmax = temp; j0 = i1; }
                int32_t j1 = (j0 + 1) % 3;
                int32_t j2 = (j1 + 1) % 3;
                FixRoundedFace(j0, j1, j2, k, c, r, d, v, result);
            }
        };
        auto faceUnbounded = [&](int32_t i0, int32_t i1, int32_t i2,
            Vector3<double> const& k, Vector3<double> const& c, double r,
            Vector3<double> const& d, Vector3<double> const& v)
        {
            if (v[i0] < 0.0)
            {
                FixRoundedFace(i0, i1, i2, k, c, r, d, v, result);
            }
        };

        if (delta[2] <= radius)
        {
            if (delta[1] <= radius)
            {
                if (delta[0] <= radius)
                {
                    if (delta[2] <= 0.0)
                    {
                        if (delta[1] <= 0.0)
                        {
                            if (delta[0] <= 0.0) { interiorOverlap(C); }
                            else { faceOverlap(0, 1, 2, K, C, radius, delta); }
                        }
                        else
                        {
                            if (delta[0] <= 0.0)
                            {
                                faceOverlap(1, 2, 0, K, C, radius, delta);
                            }
                            else if (delta[0] * delta[0] + delta[1] * delta[1] <= radius * radius)
                            {
                                edgeOverlap(0, 1, 2, K, C, radius, delta);
                            }
                            else
                            {
                                edgeSeparated(0, 1, 2, K, C, radius, delta, V);
                            }
                        }
                    }
                    else
                    {
                        if (delta[1] <= 0.0)
                        {
                            if (delta[0] <= 0.0)
                            {
                                faceOverlap(2, 0, 1, K, C, radius, delta);
                            }
                            else if (delta[0] * delta[0] + delta[2] * delta[2] <= radius * radius)
                            {
                                edgeOverlap(2, 0, 1, K, C, radius, delta);
                            }
                            else
                            {
                                edgeSeparated(2, 0, 1, K, C, radius, delta, V);
                            }
                        }
                        else
                        {
                            if (delta[0] <= 0.0)
                            {
                                if (delta[1] * delta[1] + delta[2] * delta[2] <= radius * radius)
                                {
                                    edgeOverlap(1, 2, 0, K, C, radius, delta);
                                }
                                else
                                {
                                    edgeSeparated(1, 2, 0, K, C, radius, delta, V);
                                }
                            }
                            else
                            {
                                if (Dot(delta, delta) <= radius * radius)
                                {
                                    vertexOverlap(K, radius, delta);
                                }
                                else
                                {
                                    vertexSeparated(K, radius, delta, V);
                                }
                            }
                        }
                    }
                }
                else
                {
                    faceUnbounded(0, 1, 2, K, C, radius, delta, V);
                }
            }
            else
            {
                if (delta[0] <= radius)
                {
                    faceUnbounded(1, 2, 0, K, C, radius, delta, V);
                }
                else
                {
                    edgeUnbounded(0, 1, K, C, radius, delta, V);
                }
            }
        }
        else
        {
            if (delta[1] <= radius)
            {
                if (delta[0] <= radius)
                {
                    faceUnbounded(2, 0, 1, K, C, radius, delta, V);
                }
                else
                {
                    edgeUnbounded(2, 0, K, C, radius, delta, V);
                }
            }
            else
            {
                if (delta[0] <= radius)
                {
                    edgeUnbounded(1, 2, K, C, radius, delta, V);
                }
                else
                {
                    vertexUnbounded(K, C, radius, delta, V);
                }
            }
        }

        if (result.intersectionType != 0)
        {
            for (int32_t i = 0; i < 3; ++i)
            {
                if (sign[i] < 0.0)
                {
                    result.contactPoint[i] = -result.contactPoint[i];
                }
            }
        }
    }

    // The answer the port reports: upstream's operator() with the corrected
    // DoQuery.
    Ab3S3Result Ab3S3Fixed(AlignedBox3<double> const& box,
        Vector3<double> const& boxVelocity, Sphere3<double> const& sphere,
        Vector3<double> const& sphereVelocity)
    {
        Ab3S3Result result{};
        Vector3<double> boxCenter = (box.max + box.min) * 0.5;
        Vector3<double> extent = (box.max - box.min) * 0.5;
        Vector3<double> C = sphere.center - boxCenter;
        Vector3<double> V = sphereVelocity - boxVelocity;
        AlignedBox3<double> superBox{};
        for (int32_t i = 0; i < 3; ++i)
        {
            superBox.max[i] = extent[i] + sphere.radius;
            superBox.min[i] = -superBox.max[i];
        }
        TIQuery<double, Ray3<double>, AlignedBox3<double>> rbQuery;
        auto rbResult = rbQuery(Ray3<double>(C, V), superBox);
        if (rbResult.intersect)
        {
            FixDoQuery(extent, C, sphere.radius, V, result);
            result.contactPoint += boxCenter;
        }
        return result;
    }

    bool Ab3S3Same(Ab3S3Result const& a, Ab3S3Result const& b)
    {
        if (a.intersectionType != b.intersectionType) { return false; }
        if (std::memcmp(&a.contactTime, &b.contactTime, sizeof(double)) != 0)
        {
            return false;
        }
        for (int32_t i = 0; i < 3; ++i)
        {
            if (std::memcmp(&a.contactPoint[i], &b.contactPoint[i], sizeof(double)) != 0)
            {
                return false;
            }
        }
        return true;
    }

    // Unrecorded draws for the dynamic box/sphere query; the caller records
    // them. 'aimed' sends the sphere towards a point of the box boundary,
    // which reaches the rounded face, edge and vertex branches far more often
    // than uniform velocities do.
    void RawBoxSphere(oracle::Ctx& io, int mode, bool aimed,
        Vector3<double>& boxMin, Vector3<double>& boxMax,
        Vector3<double>& boxVelocity, Vector3<double>& center, double& radius,
        Vector3<double>& sphereVelocity)
    {
        Vector3<double> extent{}, boxCenter{};
        for (int i = 0; i < 3; ++i)
        {
            double a{}, b{};
            if (mode == 0)
            {
                a = static_cast<double>(io.rawInteger(-3, 3));
                b = static_cast<double>(io.rawInteger(-3, 3));
            }
            else
            {
                a = io.raw(-3.0, 3.0);
                b = io.raw(-3.0, 3.0);
            }
            boxMin[i] = std::min(a, b);
            boxMax[i] = std::max(a, b);
            boxCenter[i] = 0.5 * (boxMax[i] + boxMin[i]);
            extent[i] = 0.5 * (boxMax[i] - boxMin[i]);
        }
        radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.25, 3.0));
        for (int i = 0; i < 3; ++i)
        {
            boxVelocity[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-1, 1))
                : io.raw(-2.0, 2.0));
        }
        if (aimed)
        {
            Vector3<double> target{}, offset{};
            for (int i = 0; i < 3; ++i)
            {
                target[i] = boxCenter[i]
                    + 0.5 * static_cast<double>(io.rawInteger(-3, 3)) * extent[i];
                offset[i] = static_cast<double>(io.rawInteger(-4, 4));
            }
            center = target + offset;
            double speed = 0.25 * static_cast<double>(io.rawInteger(1, 8));
            sphereVelocity = boxVelocity - speed * offset;
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                center[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-4, 4))
                    : io.raw(-4.0, 4.0));
            }
            for (int i = 0; i < 3; ++i)
            {
                sphereVelocity[i] = (mode == 0
                    ? static_cast<double>(io.rawInteger(-2, 2)) : io.raw(-2.0, 2.0));
            }
        }
    }

    void RecordBoxSphere(oracle::Ctx& io, Vector3<double> const& boxMin,
        Vector3<double> const& boxMax, Vector3<double> const& boxVelocity,
        Vector3<double> const& center, double radius,
        Vector3<double> const& sphereVelocity)
    {
        io.givenVec(boxMin);
        io.givenVec(boxMax);
        io.givenVec(boxVelocity);
        io.givenVec(center);
        io.given(radius);
        io.givenVec(sphereVelocity);
    }

    void Box3Sphere3Case(oracle::Ctx& io, int mode, bool aimed, bool requireDeviation)
    {
        Vector3<double> boxMin{}, boxMax{}, boxVelocity{}, center{}, sphereVelocity{};
        double radius{};
        for (;;)
        {
            RawBoxSphere(io, mode, aimed, boxMin, boxMax, boxVelocity, center,
                radius, sphereVelocity);
            AlignedBox3<double> probeBox(boxMin, boxMax);
            Sphere3<double> probeSphere{};
            probeSphere.center = center;
            probeSphere.radius = radius;
            FIQuery<double, AlignedBox3<double>, Sphere3<double>> probeQuery;
            auto upstream = probeQuery(probeBox, boxVelocity, probeSphere, sphereVelocity);
            auto fixed = Ab3S3Fixed(probeBox, boxVelocity, probeSphere, sphereVelocity);
            if (Ab3S3Same(upstream, fixed) != requireDeviation) { break; }
        }
        RecordBoxSphere(io, boxMin, boxMax, boxVelocity, center, radius,
            sphereVelocity);
        AlignedBox3<double> box(boxMin, boxMax);
        Sphere3<double> sphere{};
        sphere.center = center;
        sphere.radius = radius;
        FIQuery<double, AlignedBox3<double>, Sphere3<double>> query;
        auto r = query(box, boxVelocity, sphere, sphereVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }
}

ORACLE_CASE("IntrAlignedBox3Sphere3.test")
{
    int mode = io.index() % 2;
    auto box = ABox<3>(io, mode, 3, 4.0);
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    TIQuery<double, AlignedBox3<double>, Sphere3<double>> query;
    auto r = query(box, sphere);
    io.outBool(r.intersect);
}

// The generator rejects the inputs on which upstream's first-wins probe
// selection disagrees with the corrected one, that is exactly the inputs on
// which the port deliberately deviates (see the deviation case below).
ORACLE_CASE("IntrAlignedBox3Sphere3.find")
{
    Box3Sphere3Case(io, io.index() % 2, false, false);
}

ORACLE_CASE("IntrAlignedBox3Sphere3.find.aimed")
{
    Box3Sphere3Case(io, io.index() % 2, true, false);
}

// Deliberate deviation: the sphere enters the Minkowski sum through a rounded
// piece that upstream does not probe, or earlier than the piece upstream
// accepts. docs/UPSTREAM-FINDINGS.md IntrAlignedBox3Sphere3.h, issues #458
// and #465.
ORACLE_CASE("IntrAlignedBox3Sphere3.find.probeDeviation")
{
    Box3Sphere3Case(io, 1, true, true);
}

// =========================== IntrTriangle3Cylinder3 ======================

namespace
{
    // Upstream's DiskOverlapsPolygon reports containment of the origin when
    // every DotPerp(Q[i0], Q[i0]-Q[i1]) is zero, which is not containment but
    // degeneracy: all polygon vertices lie on one line through the origin. The
    // port falls through to the edge tests instead
    // (docs/UPSTREAM-FINDINGS.md IntrTriangle3Cylinder3.h DiskOverlapsPolygon,
    // issues #206 and #458).
    //
    // DotPerp(A, A-B) = -DotPerp(A,B), so "every term is zero" means that
    // consecutive vertices are parallel, which for a polygon with at least
    // three vertices means zero area. Every polygon the query builds is a
    // convex combination of the three projected triangle vertices, so a
    // nonzero projected triangle area implies a nonzero polygon area (the
    // clipped polygons collapse only when the triangle projection does).
    // Rejecting a zero-area projection is therefore a superset of the
    // deviating inputs, and it is exact.
    bool T3C3ProjectedAreaZero(Triangle3<double> const& triangle,
        Cylinder3<double> const& cylinder)
    {
        std::array<Vector3<double>, 3> basis{};
        basis[0] = cylinder.axis.direction;
        ComputeOrthogonalComplement(1, basis.data());
        std::array<Vector2<double>, 3> Q{};
        for (size_t i = 0; i < 3; ++i)
        {
            Vector3<double> delta = triangle.v[i] - cylinder.axis.origin;
            Q[i] = Vector2<double>{ Dot(basis[1], delta), Dot(basis[2], delta) };
        }
        return DotPerp(Q[1] - Q[0], Q[2] - Q[0]) == 0.0;
    }

    void RawTriangleCylinder(oracle::Ctx& io, int mode, Vector3<double>& v0,
        Vector3<double>& v1, Vector3<double>& v2, Vector3<double>& origin,
        Vector3<double>& direction, double& radius, double& height)
    {
        Vector3<double>* v[3] = { &v0, &v1, &v2 };
        for (int k = 0; k < 3; ++k)
        {
            for (int i = 0; i < 3; ++i)
            {
                (*v[k])[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                    : io.raw(-3.0, 3.0));
            }
        }
        for (int i = 0; i < 3; ++i)
        {
            origin[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                : io.raw(-2.0, 2.0));
        }
        if (mode == 0)
        {
            int k = io.rawInteger(0, 5);
            direction.MakeZero();
            direction[k % 3] = (k < 3 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { direction[i] = io.raw(-1.0, 1.0); }
                len = Length(direction);
            } while (len < 0.1 || len > 1.0);
            Normalize(direction);
        }
        radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.25, 3.0));
        height = (mode == 0 ? static_cast<double>(2 * io.rawInteger(1, 3))
            : io.raw(0.5, 6.0));
    }
}

// The generator rejects the triangles whose projection onto the plane
// perpendicular to the cylinder axis has zero area, a superset of the inputs
// on which the port deliberately deviates (see the deviation case below). The
// lattice mode uses an even height so that h/2 is an integer and the vertices
// land exactly on the bottom and top planes of the slab, which reaches the
// single-point and single-segment clip cases (1a, 1b, 2a, 2b of the PDF) and
// every polygon clip case.
ORACLE_CASE("IntrTriangle3Cylinder3.test")
{
    int mode = io.index() % 2;
    Vector3<double> v0{}, v1{}, v2{}, origin{}, direction{};
    double radius{}, height{};
    for (;;)
    {
        RawTriangleCylinder(io, mode, v0, v1, v2, origin, direction, radius, height);
        Triangle3<double> probeTriangle(v0, v1, v2);
        Line3<double> probeAxis(origin, direction);
        Cylinder3<double> probeCylinder(probeAxis, radius, height);
        if (!T3C3ProjectedAreaZero(probeTriangle, probeCylinder)) { break; }
    }
    io.givenVec(v0);
    io.givenVec(v1);
    io.givenVec(v2);
    io.givenVec(origin);
    io.givenVec(direction);
    io.given(radius);
    io.given(height);
    Triangle3<double> triangle(v0, v1, v2);
    Line3<double> axis(origin, direction);
    Cylinder3<double> cylinder(axis, radius, height);
    TIQuery<double, Triangle3<double>, Cylinder3<double>> query;
    auto r = query(triangle, cylinder);
    io.outBool(r.intersect);
}

// Deliberate deviation: the triangle lies in a plane containing the cylinder
// axis and well outside the cylinder, so its projection is a segment on a line
// through the origin. Every DotPerp is zero and upstream reports containment
// of the origin, hence an intersection, however far the triangle is from the
// cylinder. docs/UPSTREAM-FINDINGS.md IntrTriangle3Cylinder3.h
// DiskOverlapsPolygon, issues #206 and #458.
ORACLE_CASE("IntrTriangle3Cylinder3.test.degeneratePolygonDeviation")
{
    Vector3<double> origin = io.latticeVec<3>(-2, 2);
    int k = io.rawInteger(0, 2);
    int m = (k + 1 + io.rawInteger(0, 1)) % 3;
    Vector3<double> ek{}, em{};
    ek.MakeZero();
    ek[k] = 1.0;
    em.MakeZero();
    em[m] = 1.0;
    io.givenVec(ek);
    double radius = io.lattice(1, 3);
    int half = io.rawInteger(1, 3);
    double height = io.given(static_cast<double>(2 * half));
    // The three vertices are origin + a*W + b*U with |a| <= h/2 (inside the
    // slab) and b >= r + 1 (outside the disk of the cylinder).
    Vector3<double> v[3]{};
    for (int i = 0; i < 3; ++i)
    {
        double a = static_cast<double>(io.rawInteger(-half, half));
        double b = radius + static_cast<double>(io.rawInteger(1, 3));
        v[i] = origin + a * ek + b * em;
    }
    io.givenVec(v[0]);
    io.givenVec(v[1]);
    io.givenVec(v[2]);
    Triangle3<double> triangle(v[0], v[1], v[2]);
    Line3<double> axis(origin, ek);
    Cylinder3<double> cylinder(axis, radius, height);
    TIQuery<double, Triangle3<double>, Cylinder3<double>> query;
    auto r = query(triangle, cylinder);
    io.outBool(r.intersect);
}

// Deliberate deviation: Cylinder3 marks an infinite cylinder with the sentinel
// height = -1, which upstream feeds into the finite-cylinder formulas (an
// empty slab); the port asserts IsFinite instead.
// docs/UPSTREAM-FINDINGS.md, issues #197, #206 and #255.
ORACLE_CASE("IntrTriangle3Cylinder3.test.infiniteDeviation")
{
    auto triangle = Tri<3>(io, 1, 3, 3.0);
    Vector3<double> origin = io.vec<3>(-2.0, 2.0);
    Vector3<double> direction = io.unit<3>();
    double radius = io.real(0.25, 3.0);
    double height = io.given(-1.0);
    Line3<double> axis(origin, direction);
    Cylinder3<double> cylinder(axis, radius, height);
    TIQuery<double, Triangle3<double>, Cylinder3<double>> query;
    auto r = query(triangle, cylinder);
    io.outBool(r.intersect);
}

// ============================ IntrEllipse2Ellipse2 =======================

namespace
{
    // An ellipse. Mode 0 uses coordinate axes, lattice centres and extents
    // that are powers of two, so the standard-form matrices are exact; mode 1
    // uses a rotated frame and uniform values. The draws are unrecorded, so
    // that a case can reject a candidate; RecordEllipse writes the accepted
    // one and every case records the same layout.
    void RawEllipse(oracle::Ctx& io, int mode, Ellipse2<double>& ellipse)
    {
        for (int i = 0; i < 2; ++i)
        {
            ellipse.center[i] = (mode == 0
                ? static_cast<double>(io.rawInteger(-3, 3)) : io.raw(-3.0, 3.0));
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
        ellipse.axis[0] = Vector2<double>{ c, s };
        ellipse.axis[1] = Vector2<double>{ -s, c };
        for (int i = 0; i < 2; ++i)
        {
            ellipse.extent[i] = (mode == 0
                ? std::ldexp(1.0, io.rawInteger(-1, 1)) : io.raw(0.25, 3.0));
        }
    }

    void RecordEllipse(oracle::Ctx& io, Ellipse2<double> const& ellipse)
    {
        io.givenVec(ellipse.center);
        io.givenVec(ellipse.axis[0]);
        io.givenVec(ellipse.axis[1]);
        io.givenVec(ellipse.extent);
    }

    Ellipse2<double> Ellipse(oracle::Ctx& io, int mode)
    {
        Ellipse2<double> ellipse{};
        RawEllipse(io, mode, ellipse);
        RecordEllipse(io, ellipse);
        return ellipse;
    }

    // Upstream's TI query drops the c_i = 0 term of f(s) when exactly one of
    // c0, c1 is zero, which loses the two critical points at the pole
    // s = 1/d_i; the port adds them back (docs/UPSTREAM-FINDINGS.md
    // IntrEllipse2Ellipse2.h TIQuery, issue #458). This replicates upstream's
    // control flow down to the 'valid' list and reports whether a term is
    // dropped, which is exactly the set on which the port deviates.
    bool E2E2OneTermDropped(Ellipse2<double> const& ellipse0,
        Ellipse2<double> const& ellipse1)
    {
        double const zero = 0.0, one = 1.0;
        Matrix2x2<double> R0{};
        R0.SetCol(0, ellipse0.axis[0]);
        R0.SetCol(1, ellipse0.axis[1]);
        Matrix2x2<double> R1{};
        R1.SetCol(0, ellipse1.axis[0]);
        R1.SetCol(1, ellipse1.axis[1]);
        Matrix2x2<double> D1
        {
            one / (ellipse1.extent[0] * ellipse1.extent[0]), zero,
            zero, one / (ellipse1.extent[1] * ellipse1.extent[1])
        };
        Matrix2x2<double> D0NegHalf{ ellipse0.extent[0], zero, zero, ellipse0.extent[1] };
        Matrix2x2<double> D0Half
        {
            one / ellipse0.extent[0], zero, zero, one / ellipse0.extent[1]
        };
        Vector2<double> K2 = D0Half * ((ellipse1.center - ellipse0.center) * R0);
        Matrix2x2<double> R1TR0D0NegHalf = MultiplyATB(R1, R0 * D0NegHalf);
        Matrix2x2<double> M2 = MultiplyATB(R1TR0D0NegHalf, D1) * R1TR0D0NegHalf;
        SymmetricEigensolver2x2<double> es{};
        std::array<double, 2> D{};
        std::array<std::array<double, 2>, 2> evec{};
        es(M2(0, 0), M2(0, 1), M2(1, 1), +1, D, evec);
        Matrix2x2<double> R{};
        R.SetCol(0, evec[0]);
        R.SetCol(1, evec[1]);
        Vector2<double> K = K2 * R;
        if (K == Vector2<double>::Zero())
        {
            return false;
        }
        double d0 = D[0], d1 = D[1];
        double c0 = K[0] * K[0], c1 = K[1] * K[1];
        double first0{}, second0{}, first1{}, second1{};
        if (d0 >= d1)
        {
            first0 = d0; second0 = c0; first1 = d1; second1 = c1;
        }
        else
        {
            first0 = d1; second0 = c1; first1 = d0; second1 = c0;
        }
        if (first0 > first1)
        {
            int numValid = (second0 > zero ? 1 : 0) + (second1 > zero ? 1 : 0);
            return numValid == 1;
        }
        return false;
    }

    // Upstream reports "the ellipses are the same" with
    // numPoints = SIZE_MAX, which is not representable as a double below
    // 2^53; the sentinel is emitted as -1 on both sides.
    int32_t E2E2NumPoints(size_t numPoints)
    {
        return (numPoints == std::numeric_limits<size_t>::max()
            ? -1 : static_cast<int32_t>(numPoints));
    }

    void OutEllipseFI(oracle::Ctx& io,
        FIQuery<double, Ellipse2<double>, Ellipse2<double>>::Result const& r)
    {
        io.outBool(r.intersect);
        io.outInt(E2E2NumPoints(r.numPoints));
        for (size_t i = 0; i < 4; ++i)
        {
            io.outBool(r.isTransverse[i]);
        }
        for (size_t i = 0; i < 4; ++i)
        {
            io.outVec(r.points[i]);
        }
    }
}

// The classification is a discrete output, so it is compared exactly even
// though the query calls the C math library (std::pow for the bracket of the
// middle roots, and the RootsBisection iteration). The generator rejects the
// inputs on which upstream drops one term of f(s), which is exactly the set on
// which the port deliberately deviates (see the deviation case below).
ORACLE_CASE("IntrEllipse2Ellipse2.test")
{
    int mode = io.index() % 2;
    Ellipse2<double> ellipse0{}, ellipse1{};
    for (;;)
    {
        RawEllipse(io, mode, ellipse0);
        RawEllipse(io, mode, ellipse1);
        if (!E2E2OneTermDropped(ellipse0, ellipse1)) { break; }
    }
    RecordEllipse(io, ellipse0);
    RecordEllipse(io, ellipse1);
    TIQuery<double, Ellipse2<double>, Ellipse2<double>> query;
    auto classification = query(ellipse0, ellipse1);
    io.outInt(static_cast<int32_t>(classification));
}

// Deliberate deviation: ellipse0 is the unit circle and ellipse1 is
// axis-aligned with its centre offset along one coordinate axis, so exactly
// one of c0, c1 is zero. Upstream drops that term of f(s) and loses the two
// critical points at the pole s = 1/d_j, which makes it report containment for
// ellipses that overlap. docs/UPSTREAM-FINDINGS.md IntrEllipse2Ellipse2.h
// TIQuery, issue #458.
ORACLE_CASE("IntrEllipse2Ellipse2.test.poleDeviation")
{
    Vector2<double> center = io.latticeVec<2>(-3, 3);
    Vector2<double> e0{ 1.0, 0.0 }, e1{ 0.0, 1.0 };
    io.givenVec(e0);
    io.givenVec(e1);
    Vector2<double> unitExtent{ 1.0, 1.0 };
    io.givenVec(unitExtent);
    Ellipse2<double> ellipse0{};
    ellipse0.center = center;
    ellipse0.axis[0] = e0;
    ellipse0.axis[1] = e1;
    ellipse0.extent = unitExtent;

    int k = io.rawInteger(0, 1);
    double t = 0.25 * static_cast<double>(io.rawInteger(1, 4));
    Vector2<double> offset{};
    offset.MakeZero();
    offset[k] = t;
    Vector2<double> center1 = io.givenVec(center + offset);
    io.givenVec(e0);
    io.givenVec(e1);
    // Narrow across the offset direction and tall along it, the configuration
    // of the documented example.
    Vector2<double> extent1{};
    extent1[k] = std::ldexp(1.0, io.rawInteger(1, 2));
    extent1[1 - k] = std::ldexp(1.0, io.rawInteger(-2, -1));
    io.givenVec(extent1);
    Ellipse2<double> ellipse1{};
    ellipse1.center = center1;
    ellipse1.axis[0] = e0;
    ellipse1.axis[1] = e1;
    ellipse1.extent = extent1;

    TIQuery<double, Ellipse2<double>, Ellipse2<double>> query;
    auto classification = query(ellipse0, ellipse1);
    io.outInt(static_cast<int32_t>(classification));
}

// The find-intersection query. Its compared path calls the C math library
// through RootsPolynomial::SolveQuartic (cbrt and the trigonometric resolvent),
// so the points are compared with a tolerance; 'numPoints' and 'isTransverse'
// are discrete and compared exactly.
ORACLE_CASE("IntrEllipse2Ellipse2.find")
{
    int mode = io.index() % 2;
    auto ellipse0 = Ellipse(io, mode);
    auto ellipse1 = Ellipse(io, mode);
    bool earlyExit = io.boolean();
    FIQuery<double, Ellipse2<double>, Ellipse2<double>> query;
    auto r = query(ellipse0, ellipse1, earlyExit);
    OutEllipseFI(io, r);
}

// The same query through the standard-form overload. The matrices are built
// from ellipses (so they are symmetric and positive definite) and recorded
// entry by entry.
ORACLE_CASE("IntrEllipse2Ellipse2.findStandardForm")
{
    int mode = io.index() % 2;
    Ellipse2<double> ellipse0{}, ellipse1{};
    RawEllipse(io, mode, ellipse0);
    RawEllipse(io, mode, ellipse1);
    FIQuery<double, Ellipse2<double>, Ellipse2<double>> query;
    Vector2<double> C0{}, C1{};
    Matrix2x2<double> M0{}, M1{};
    query.GetStandardForm(ellipse0, C0, M0);
    query.GetStandardForm(ellipse1, C1, M1);
    io.givenVec(C0);
    io.given(M0(0, 0));
    io.given(M0(0, 1));
    io.given(M0(1, 1));
    io.givenVec(C1);
    io.given(M1(0, 0));
    io.given(M1(0, 1));
    io.given(M1(1, 1));
    bool earlyExit = io.boolean();
    auto r = query(C0, M0, C1, M1, earlyExit);
    OutEllipseFI(io, r);
}

// GetStandardForm and both ComputeAlignedBox overloads.
ORACLE_CASE("IntrEllipse2Ellipse2.getStandardForm")
{
    int mode = io.index() % 2;
    auto ellipse = Ellipse(io, mode);
    FIQuery<double, Ellipse2<double>, Ellipse2<double>> query;
    Vector2<double> C{};
    Matrix2x2<double> M{};
    query.GetStandardForm(ellipse, C, M);
    io.outVec(C);
    io.outReal(M(0, 0));
    io.outReal(M(0, 1));
    io.outReal(M(1, 0));
    io.outReal(M(1, 1));
    AlignedBox2<double> box0{}, box1{};
    query.ComputeAlignedBox(ellipse, box0);
    query.ComputeAlignedBox(C, M, box1);
    io.outVec(box0.min);
    io.outVec(box0.max);
    io.outVec(box1.min);
    io.outVec(box1.max);
}

// Deliberate deviation: ellipse0 is a circle of radius 2^m centred at K and
// ellipse1 is concentric with it, with the exact (non-unit) axes (p,q) and
// (-q,p) and extents chosen so that M1(1,1) = M0(1,1) exactly. Then e0, e1 and
// e2 are all zero and y0 = 0 is a double root of the quartic, so
// 'divisor = e2 + e4*y0' is exactly zero. In that branch upstream writes both
// symmetric intersection points to result.points[numPoints] without
// incrementing numPoints between the two writes, so the first point is lost;
// the port stores both. docs/UPSTREAM-FINDINGS.md IntrEllipse2Ellipse2.h
// CaseE4NotZero, issue #250.
//
// The query divides by |axis[0]|^2 (GetStandardForm's USqrLen), which is what
// makes the non-unit axes legitimate here and keeps the matrices exact.
ORACLE_CASE("IntrEllipse2Ellipse2.find.divisorDeviation")
{
    Vector2<double> center = io.latticeVec<2>(-3, 3);
    int m = io.rawInteger(-1, 1);
    double radius = std::ldexp(1.0, m);
    Vector2<double> e0{ 1.0, 0.0 }, e1{ 0.0, 1.0 };
    io.givenVec(e0);
    io.givenVec(e1);
    Vector2<double> extent0{ radius, radius };
    io.givenVec(extent0);
    Ellipse2<double> ellipse0{};
    ellipse0.center = center;
    ellipse0.axis[0] = e0;
    ellipse0.axis[1] = e1;
    ellipse0.extent = extent0;

    // (p,q) with p^2 + q^2 = 5 and the extents (2^(m+1), 2^(m-1)) or their
    // swap, so that (q^2/a^2 + p^2/b^2)/(p^2+q^2) = 1/radius^2 exactly.
    int swap = io.rawInteger(0, 1);
    double sign0 = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    double sign1 = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    double p = (swap == 0 ? 1.0 : 2.0) * sign0;
    double q = (swap == 0 ? 2.0 : 1.0) * sign1;
    Vector2<double> u{ p, q };
    Vector2<double> v{ -q, p };
    io.givenVec(center);
    io.givenVec(u);
    io.givenVec(v);
    Vector2<double> extent1{};
    extent1[0] = std::ldexp(1.0, (swap == 0 ? m + 1 : m - 1));
    extent1[1] = std::ldexp(1.0, (swap == 0 ? m - 1 : m + 1));
    io.givenVec(extent1);
    Ellipse2<double> ellipse1{};
    ellipse1.center = center;
    ellipse1.axis[0] = u;
    ellipse1.axis[1] = v;
    ellipse1.extent = extent1;

    FIQuery<double, Ellipse2<double>, Ellipse2<double>> query;
    auto r = query(ellipse0, ellipse1, false);
    OutEllipseFI(io, r);
}
