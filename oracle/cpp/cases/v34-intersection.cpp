// Verify group 34 (intersection): differential cases for the Intr* headers
// listed in plan/verify-groups.json group 34.
//
// Every generator mixes an exactly representable small-lattice mode with a
// uniform mode and, where a branch needs it, a constructed mode, so that the
// touching, tangent, coplanar, parallel and degenerate branches are reached.
// Most of the queries in this group use only + - * / sqrt fabs and
// comparisons, so most cases are declared exact on the TypeScript side; the
// sin/cos used to build orthonormal frames, unit directions and arc endpoints
// are applied to *unrecorded* draws and only the resulting frame, unit vector
// or endpoint is recorded as an input, so libm never enters the compared
// computation. The trigonometry that Cone::SetAngle derives from the angle is
// likewise recorded as inputs and assigned by the replay. The only case whose
// compared path reaches the C math library is IntrAreaEllipse2Ellipse2
// (atan2, atan, sin, cos, and RootsPolynomial::SolveQuartic underneath), which
// is the one case compared with a tolerance; the symmetric eigensolve that
// Ellipse2::FromCoefficients runs for IntrPlane3Cylinder3 uses only arithmetic
// and sqrt and is compared exactly.
#define ORACLE_FAMILY "v34-intersection"
#include "Oracle.h"

#include <Mathematics/IntrAlignedBox3Cone3.h>
#include <Mathematics/IntrAreaEllipse2Ellipse2.h>
#include <Mathematics/IntrConvexMesh3Plane3.h>
#include <Mathematics/IntrLine3Cone3.h>
#include <Mathematics/IntrLine3Plane3.h>
#include <Mathematics/IntrOrientedBox2Circle2.h>
#include <Mathematics/IntrOrientedBox2Cone2.h>
#include <Mathematics/IntrOrientedBox3Sphere3.h>
#include <Mathematics/IntrPlane3Capsule3.h>
#include <Mathematics/IntrPlane3Circle3.h>
#include <Mathematics/IntrPlane3Cylinder3.h>
#include <Mathematics/IntrPlane3Ellipsoid3.h>
#include <Mathematics/IntrPlane3OrientedBox3.h>
#include <Mathematics/IntrPlane3Sphere3.h>
#include <Mathematics/IntrRay2Arc2.h>
#include <Mathematics/IntrRay2SegmentMesh2.h>
#include <Mathematics/IntrRay3Ellipsoid3.h>
#include <Mathematics/IntrSegment2Arc2.h>
#include <Mathematics/IntrSegment2SegmentMesh2.h>
#include <Mathematics/IntrSegment3Ellipsoid3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
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

    // A unit-length direction: an exact coordinate axis in lattice mode, a
    // uniformly distributed unit vector otherwise.
    template <int N>
    Vector<N, double> Dir(oracle::Ctx& io, int mode)
    {
        return (mode == 0 ? AxisUnit<N>(io) : io.unit<N>());
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

    // Orthonormal 2D frame; mode 0 gives a signed coordinate pair (exact).
    std::array<Vector2<double>, 2> Frame2(oracle::Ctx& io, int mode)
    {
        std::array<Vector2<double>, 2> v{};
        if (mode == 0)
        {
            int k = io.rawInteger(0, 3);
            v[0].MakeZero();
            v[0][k % 2] = (k < 2 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 2; ++i) { v[0][i] = io.raw(-1.0, 1.0); }
                len = Length(v[0]);
            } while (len < 0.1 || len > 1.0);
            Normalize(v[0]);
        }
        v[1] = Perp(v[0]);
        std::array<Vector2<double>, 2> axis{};
        axis[0] = io.givenVec(v[0]);
        axis[1] = io.givenVec(v[1]);
        return axis;
    }

    template <int N>
    Segment<N, double> Sg(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> p0 = Pt<N>(io, mode, lat, range);
        Vector<N, double> p1 = Pt<N>(io, mode, lat, range);
        return Segment<N, double>(p0, p1);
    }

    template <int N>
    Hypersphere<N, double> Sph(oracle::Ctx& io, int mode, int lat, double range)
    {
        Hypersphere<N, double> sphere{};
        sphere.center = Pt<N>(io, mode, lat, range);
        sphere.radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        return sphere;
    }

    // A plane through a point with a unit-length (or exact coordinate-axis)
    // normal. Records the normal and then the point; the plane constant is
    // Dot(normal, point), computed identically on both sides.
    Plane3<double> Pln3(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector3<double> normal = Dir<3>(io, mode);
        Vector3<double> origin = Pt<3>(io, mode, lat, range);
        return Plane3<double>(normal, origin);
    }

    // Oriented boxes in 2D and 3D.
    OrientedBox2<double> OBox2(oracle::Ctx& io, int mode, int lat, double range)
    {
        OrientedBox2<double> box{};
        box.center = Pt<2>(io, mode, lat, range);
        box.axis = Frame2(io, mode);
        box.extent = Extent<2>(io, mode, lat, range);
        return box;
    }

    OrientedBox3<double> OBox3(oracle::Ctx& io, int mode, int lat, double range)
    {
        OrientedBox3<double> box{};
        box.center = Pt<3>(io, mode, lat, range);
        box.axis = Frame3(io, mode);
        box.extent = Extent<3>(io, mode, lat, range);
        return box;
    }

    // An ellipsoid (or ellipse) with an orthonormal axis frame and positive
    // extents.
    Ellipsoid3<double> Ellipsoid(oracle::Ctx& io, int mode, int lat, double range)
    {
        Ellipsoid3<double> ellipsoid{};
        ellipsoid.center = Pt<3>(io, mode, lat, range);
        auto axis = Frame3(io, mode);
        ellipsoid.axis = axis;
        ellipsoid.extent = Extent<3>(io, mode, lat, range);
        return ellipsoid;
    }

    Ellipse2<double> Ellipse(oracle::Ctx& io, int mode, int lat, double range)
    {
        Ellipse2<double> ellipse{};
        ellipse.center = Pt<2>(io, mode, lat, range);
        auto axis = Frame2(io, mode);
        ellipse.axis = axis;
        ellipse.extent = Extent<2>(io, mode, lat, range);
        return ellipse;
    }

    // A cone. The six values derived from the angle by Cone::SetAngle are
    // computed with cos/sin/tan and RECORDED as inputs, so the C math library
    // never enters the compared computation: the replay assigns them.
    // 'kind' selects infinite (0), infinite truncated (1), finite (2) or
    // frustum (3). The height pair is recorded last, with maxHeight = -1 for
    // the infinite kinds, which is upstream's sentinel.
    template <int N>
    Cone<N, double> Cn(oracle::Ctx& io, int mode, int lat, double range, int kind)
    {
        Cone<N, double> cone{};
        cone.ray.origin = Pt<N>(io, mode, lat, range);
        cone.ray.direction = Dir<N>(io, mode);
        double angle = io.raw(0.15, 1.35);
        cone.SetAngle(angle);
        io.given(cone.angle);
        io.given(cone.cosAngle);
        io.given(cone.sinAngle);
        io.given(cone.tanAngle);
        io.given(cone.cosAngleSqr);
        io.given(cone.sinAngleSqr);
        io.given(cone.invSinAngle);
        double minHeight = 0.0, maxHeight = -1.0;
        if (kind == 1 || kind == 3)
        {
            minHeight = (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
                : io.raw(0.25, range));
        }
        if (kind == 2 || kind == 3)
        {
            maxHeight = minHeight + (mode == 0
                ? static_cast<double>(io.rawInteger(1, lat)) : io.raw(0.25, range));
        }
        io.given(minHeight);
        io.given(maxHeight);
        if (maxHeight < 0.0)
        {
            cone.MakeInfiniteTruncatedCone(minHeight);
        }
        else
        {
            cone.MakeConeFrustum(minHeight, maxHeight);
        }
        return cone;
    }

    // Exposes the protected DoQuery helpers, which the port exports as public
    // free functions.
    template <typename Q>
    struct Expose : public Q
    {
        using Q::DoQuery;
    };
}

// ============================== IntrLine3Plane3 ==========================

ORACLE_CASE("IntrLine3Plane3.test")
{
    int mode = io.index() % 2;
    auto origin = Pt<3>(io, mode, 3, 4.0);
    auto direction = Dir<3>(io, mode);
    auto plane = Pln3(io, mode, 3, 4.0);
    Line3<double> line(origin, direction);
    TIQuery<double, Line3<double>, Plane3<double>> query;
    auto r = query(line, plane);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrLine3Plane3.find")
{
    int mode = io.index() % 2;
    auto origin = Pt<3>(io, mode, 3, 4.0);
    auto direction = Dir<3>(io, mode);
    auto plane = Pln3(io, mode, 3, 4.0);
    Line3<double> line(origin, direction);
    FIQuery<double, Line3<double>, Plane3<double>> query;
    auto r = query(line, plane);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter);
    io.outVec(r.point);
}

ORACLE_CASE("IntrLine3Plane3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto origin = Pt<3>(io, mode, 3, 4.0);
    auto direction = Dir<3>(io, mode);
    auto plane = Pln3(io, mode, 3, 4.0);
    Expose<FIQuery<double, Line3<double>, Plane3<double>>> query;
    FIQuery<double, Line3<double>, Plane3<double>>::Result r{};
    query.DoQuery(origin, direction, plane, r);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter);
    io.outVec(r.point);
}

// A uniform line is never parallel to a uniform plane, so neither the
// "parallel and distinct" branch nor the "line in the plane" branch
// (numIntersections = max int32) is reached by the cases above. Here the
// plane normal is an exact coordinate axis, the line direction is an integer
// vector perpendicular to it (so Dot(D,N) is exactly zero) and the line
// origin is on the plane on half the records.
ORACLE_CASE("IntrLine3Plane3.find.parallel")
{
    Vector3<double> normal{};
    int k = io.rawInteger(0, 2);
    normal.MakeZero();
    normal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(normal);
    Vector3<double> planeOrigin = io.latticeVec<3>(-3, 3);
    Vector3<double> direction{};
    for (int i = 0; i < 3; ++i)
    {
        direction[i] = (i == k ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
    }
    if (direction == Vector3<double>::Zero()) { direction[(k + 1) % 3] = 1.0; }
    io.givenVec(direction);
    // On half the records the line origin is on the plane exactly.
    Vector3<double> lineOrigin = planeOrigin;
    bool offPlane = (io.rawInteger(0, 1) == 1);
    for (int i = 0; i < 3; ++i)
    {
        lineOrigin[i] += static_cast<double>(io.rawInteger(-3, 3))
            * (i == k ? 0.0 : 1.0);
    }
    if (offPlane) { lineOrigin[k] += 1.0; }
    io.givenVec(lineOrigin);
    Plane3<double> plane(normal, planeOrigin);
    Line3<double> line(lineOrigin, direction);
    FIQuery<double, Line3<double>, Plane3<double>> query;
    auto r = query(line, plane);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    io.outReal(r.parameter);
    io.outVec(r.point);
}

// ============================= IntrPlane3Capsule3 ========================

namespace
{
    Capsule3<double> Cap3(oracle::Ctx& io, int mode, int lat, double range)
    {
        Capsule3<double> capsule{};
        capsule.segment = Sg<3>(io, mode, lat, range);
        capsule.radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        return capsule;
    }
}

ORACLE_CASE("IntrPlane3Capsule3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto capsule = Cap3(io, mode, 3, 4.0);
    TIQuery<double, Plane3<double>, Capsule3<double>> query;
    auto r = query(plane, capsule);
    io.outBool(r.intersect);
}

// The "endpoint exactly on the plane" and "|signed distance| exactly equal to
// the radius" boundaries. The plane normal is a coordinate axis, so the
// signed distance of a lattice point is an exact integer.
ORACLE_CASE("IntrPlane3Capsule3.test.touching")
{
    Vector3<double> normal{};
    int k = io.rawInteger(0, 2);
    normal.MakeZero();
    normal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(normal);
    Vector3<double> planeOrigin = io.latticeVec<3>(-3, 3);
    double radius = io.lattice(1, 3);
    // Place both endpoints at integer heights relative to the plane, in
    // [-radius-1, radius+1], so the tests sdistance0*sdistance1 <= 0 and
    // |sdistance| <= radius are evaluated at exact equality often.
    int span = static_cast<int>(radius) + 1;
    Vector3<double> p0 = planeOrigin, p1 = planeOrigin;
    for (int i = 0; i < 3; ++i)
    {
        p0[i] += static_cast<double>(io.rawInteger(-3, 3));
        p1[i] += static_cast<double>(io.rawInteger(-3, 3));
    }
    p0[k] = planeOrigin[k] + normal[k] * static_cast<double>(io.rawInteger(-span, span));
    p1[k] = planeOrigin[k] + normal[k] * static_cast<double>(io.rawInteger(-span, span));
    io.givenVec(p0);
    io.givenVec(p1);
    Plane3<double> plane(normal, planeOrigin);
    Capsule3<double> capsule{};
    capsule.segment = Segment3<double>(p0, p1);
    capsule.radius = radius;
    TIQuery<double, Plane3<double>, Capsule3<double>> query;
    auto r = query(plane, capsule);
    io.outBool(r.intersect);
}

// ============================ IntrPlane3Ellipsoid3 =======================

ORACLE_CASE("IntrPlane3Ellipsoid3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto ellipsoid = Ellipsoid(io, mode, 3, 4.0);
    TIQuery<double, Plane3<double>, Ellipsoid3<double>> query;
    auto r = query(plane, ellipsoid);
    io.outBool(r.intersect);
}

// A plane tangent to the ellipsoid: with the ellipsoid axis-aligned on the
// lattice and the plane normal a coordinate axis, sqrt(N^T M^{-1} N) is the
// corresponding extent exactly, so 'distance <= root' is evaluated at exact
// equality when the plane is placed at that distance from the centre.
ORACLE_CASE("IntrPlane3Ellipsoid3.test.tangent")
{
    Ellipsoid3<double> ellipsoid{};
    ellipsoid.center = io.latticeVec<3>(-3, 3);
    std::array<Vector3<double>, 3> axis{};
    for (int i = 0; i < 3; ++i) { axis[i].MakeZero(); axis[i][i] = 1.0; }
    io.givenVec(axis[0]);
    io.givenVec(axis[1]);
    io.givenVec(axis[2]);
    ellipsoid.axis = axis;
    ellipsoid.extent = io.latticeVec<3>(1, 4);
    int k = io.rawInteger(0, 2);
    Vector3<double> normal{};
    normal.MakeZero();
    normal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(normal);
    double offset = static_cast<double>(io.rawInteger(-1, 1));
    Vector3<double> planeOrigin = ellipsoid.center;
    planeOrigin[k] += normal[k] * (ellipsoid.extent[k] + offset);
    io.givenVec(planeOrigin);
    Plane3<double> plane(normal, planeOrigin);
    TIQuery<double, Plane3<double>, Ellipsoid3<double>> query;
    auto r = query(plane, ellipsoid);
    io.outBool(r.intersect);
}

// =========================== IntrPlane3OrientedBox3 ======================

ORACLE_CASE("IntrPlane3OrientedBox3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto box = OBox3(io, mode, 3, 4.0);
    TIQuery<double, Plane3<double>, OrientedBox3<double>> query;
    auto r = query(plane, box);
    io.outBool(r.intersect);
}

// An axis-aligned lattice box and a coordinate-axis plane normal, so that the
// projection radius and the centre distance are exact integers and
// 'distance <= radius' is evaluated at exact equality.
ORACLE_CASE("IntrPlane3OrientedBox3.test.touching")
{
    OrientedBox3<double> box{};
    box.center = io.latticeVec<3>(-3, 3);
    std::array<Vector3<double>, 3> axis{};
    for (int i = 0; i < 3; ++i) { axis[i].MakeZero(); axis[i][i] = 1.0; }
    io.givenVec(axis[0]);
    io.givenVec(axis[1]);
    io.givenVec(axis[2]);
    box.axis = axis;
    box.extent = io.latticeVec<3>(1, 3);
    int k = io.rawInteger(0, 2);
    Vector3<double> normal{};
    normal.MakeZero();
    normal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(normal);
    Vector3<double> planeOrigin = box.center;
    planeOrigin[k] += normal[k] * (box.extent[k]
        + static_cast<double>(io.rawInteger(-1, 1)));
    io.givenVec(planeOrigin);
    Plane3<double> plane(normal, planeOrigin);
    TIQuery<double, Plane3<double>, OrientedBox3<double>> query;
    auto r = query(plane, box);
    io.outBool(r.intersect);
}

// ============================= IntrPlane3Sphere3 =========================

ORACLE_CASE("IntrPlane3Sphere3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    TIQuery<double, Plane3<double>, Sphere3<double>> query;
    auto r = query(plane, sphere);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrPlane3Sphere3.find")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    FIQuery<double, Plane3<double>, Sphere3<double>> query;
    auto r = query(plane, sphere);
    io.outBool(r.intersect);
    io.outBool(r.isCircle);
    io.outVec(r.circle.center);
    io.outVec(r.circle.normal);
    io.outReal(r.circle.radius);
    io.outVec(r.point);
}

// The tangency branch (distance exactly equal to the radius) is unreachable
// with uniform inputs. With a coordinate-axis normal and a lattice sphere the
// distance is an exact integer.
ORACLE_CASE("IntrPlane3Sphere3.find.tangent")
{
    Sphere3<double> sphere{};
    sphere.center = io.latticeVec<3>(-3, 3);
    sphere.radius = io.lattice(1, 3);
    int k = io.rawInteger(0, 2);
    Vector3<double> normal{};
    normal.MakeZero();
    normal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(normal);
    Vector3<double> planeOrigin = sphere.center;
    planeOrigin[k] += normal[k] * (sphere.radius
        + static_cast<double>(io.rawInteger(-2, 1)));
    io.givenVec(planeOrigin);
    Plane3<double> plane(normal, planeOrigin);
    FIQuery<double, Plane3<double>, Sphere3<double>> query;
    auto r = query(plane, sphere);
    io.outBool(r.intersect);
    io.outBool(r.isCircle);
    io.outVec(r.circle.center);
    io.outVec(r.circle.normal);
    io.outReal(r.circle.radius);
    io.outVec(r.point);
}

// ========================== IntrOrientedBox2Circle2 ======================

namespace
{
    using Ab2C2Result = FIQuery<double, AlignedBox2<double>, Circle2<double>>::Result;

    // Unrecorded draws for the dynamic oriented-box/circle query. 'aimed'
    // sends the circle towards a point of the box boundary, which reaches the
    // rounded face, edge and vertex branches of the Minkowski sum far more
    // often than uniform velocities do.
    void RawOBox2Circle2(oracle::Ctx& io, int mode, bool aimed,
        Vector2<double>& center, std::array<Vector2<double>, 2>& axis,
        Vector2<double>& extent, Vector2<double>& boxVelocity,
        Vector2<double>& circleCenter, double& radius,
        Vector2<double>& circleVelocity)
    {
        for (int i = 0; i < 2; ++i)
        {
            center[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
        }
        if (mode == 0)
        {
            int k = io.rawInteger(0, 3);
            axis[0].MakeZero();
            axis[0][k % 2] = (k < 2 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 2; ++i) { axis[0][i] = io.raw(-1.0, 1.0); }
                len = Length(axis[0]);
            } while (len < 0.1 || len > 1.0);
            Normalize(axis[0]);
        }
        axis[1] = Perp(axis[0]);
        for (int i = 0; i < 2; ++i)
        {
            extent[i] = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
                : io.raw(0.25, 3.0));
        }
        for (int i = 0; i < 2; ++i)
        {
            boxVelocity[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                : io.raw(-2.0, 2.0));
        }
        radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.25, 3.0));
        if (aimed)
        {
            // A dyadic point of the box boundary, in world coordinates.
            int side = io.rawInteger(0, 3);
            double t = 0.25 * static_cast<double>(io.rawInteger(-4, 4));
            Vector2<double> local{};
            local[side % 2] = (side < 2 ? 1.0 : -1.0) * extent[side % 2];
            local[1 - side % 2] = t * extent[1 - side % 2];
            Vector2<double> target = center + local[0] * axis[0] + local[1] * axis[1];
            Vector2<double> offset{};
            for (int i = 0; i < 2; ++i)
            {
                offset[i] = static_cast<double>(io.rawInteger(-4, 4));
            }
            circleCenter = target + offset;
            double speed = 0.25 * static_cast<double>(io.rawInteger(1, 8));
            circleVelocity = boxVelocity - speed * offset;
        }
        else
        {
            for (int i = 0; i < 2; ++i)
            {
                circleCenter[i] = (mode == 0
                    ? static_cast<double>(io.rawInteger(-4, 4)) : io.raw(-4.0, 4.0));
            }
            for (int i = 0; i < 2; ++i)
            {
                circleVelocity[i] = (mode == 0
                    ? static_cast<double>(io.rawInteger(-2, 2)) : io.raw(-2.0, 2.0));
            }
        }
    }

    void OBox2Circle2Case(oracle::Ctx& io, int mode, bool aimed)
    {
        Vector2<double> center{}, extent{}, boxVelocity{}, circleCenter{},
            circleVelocity{};
        std::array<Vector2<double>, 2> axis{};
        double radius{};
        RawOBox2Circle2(io, mode, aimed, center, axis, extent, boxVelocity,
            circleCenter, radius, circleVelocity);
        OrientedBox2<double> box{};
        box.center = io.givenVec(center);
        box.axis[0] = io.givenVec(axis[0]);
        box.axis[1] = io.givenVec(axis[1]);
        box.extent = io.givenVec(extent);
        io.givenVec(boxVelocity);
        Circle2<double> circle{};
        circle.center = io.givenVec(circleCenter);
        circle.radius = io.given(radius);
        io.givenVec(circleVelocity);
        FIQuery<double, OrientedBox2<double>, Circle2<double>> query;
        auto r = query(box, boxVelocity, circle, circleVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }
}

ORACLE_CASE("IntrOrientedBox2Circle2.test")
{
    int mode = io.index() % 2;
    auto box = OBox2(io, mode, 3, 4.0);
    auto circle = Sph<2>(io, mode, 3, 4.0);
    TIQuery<double, OrientedBox2<double>, Circle2<double>> query;
    auto r = query(box, circle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrOrientedBox2Circle2.find")
{
    OBox2Circle2Case(io, io.index() % 2, false);
}

ORACLE_CASE("IntrOrientedBox2Circle2.find.aimed")
{
    OBox2Circle2Case(io, io.index() % 2, true);
}

// =========================== IntrOrientedBox2Cone2 =======================

ORACLE_CASE("IntrOrientedBox2Cone2.test")
{
    int mode = io.index() % 2;
    auto box = OBox2(io, mode, 3, 4.0);
    auto cone = Cn<2>(io, mode, 3, 4.0, io.index() % 4);
    TIQuery<double, OrientedBox2<double>, Cone<2, double>> query;
    auto r = query(box, cone);
    io.outBool(r.intersect);
}

// A uniform cone rarely contains a box corner, so the quick-acceptance branch
// (the cone axis hits the box) and the corner branch are both rare. Here the
// cone vertex is placed so that the axis points at the box centre, with an
// integer perturbation, which makes the two branches roughly equally likely.
ORACLE_CASE("IntrOrientedBox2Cone2.test.aimed")
{
    OrientedBox2<double> box{};
    box.center = io.latticeVec<2>(-3, 3);
    box.axis = Frame2(io, 0);
    box.extent = Extent<2>(io, 0, 3, 4.0);
    Vector2<double> direction = AxisUnit<2>(io);
    // A negative distance puts the box behind the cone vertex, and the offset
    // is comparable to the box size, so both answers occur.
    double distance = static_cast<double>(io.rawInteger(-4, 8));
    Vector2<double> offset{};
    offset[0] = static_cast<double>(io.rawInteger(-5, 5));
    offset[1] = static_cast<double>(io.rawInteger(-5, 5));
    Vector2<double> origin = box.center - distance * direction + offset;
    io.givenVec(origin);
    Cone<2, double> cone{};
    cone.ray.origin = origin;
    cone.ray.direction = direction;
    double angle = io.raw(0.15, 1.35);
    cone.SetAngle(angle);
    io.given(cone.angle);
    io.given(cone.cosAngle);
    io.given(cone.sinAngle);
    io.given(cone.tanAngle);
    io.given(cone.cosAngleSqr);
    io.given(cone.sinAngleSqr);
    io.given(cone.invSinAngle);
    TIQuery<double, OrientedBox2<double>, Cone<2, double>> query;
    auto r = query(box, cone);
    io.outBool(r.intersect);
}

// ========================== IntrOrientedBox3Sphere3 ======================
//
// The find-intersection query is a thin wrapper around the protected DoQuery
// of FIQuery<AlignedBox3,Sphere3>, so it inherits that query's defect and the
// port's correction of it, first demonstrated by group 32
// (oracle/cpp/cases/v32-intersection.cpp). The corrected DoQuery below is
// copied verbatim from that file.

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

    // The answer the port reports: upstream's oriented-box wrapper with the
    // corrected DoQuery.
    Ab3S3Result Ob3S3Fixed(OrientedBox3<double> const& box,
        Vector3<double> const& boxVelocity, Sphere3<double> const& sphere,
        Vector3<double> const& sphereVelocity)
    {
        Ab3S3Result result{};
        result.intersectionType = 0;
        result.contactTime = 0.0;
        result.contactPoint = { 0.0, 0.0, 0.0 };
        Vector3<double> temp = sphere.center - box.center;
        Vector3<double> C
        {
            Dot(temp, box.axis[0]),
            Dot(temp, box.axis[1]),
            Dot(temp, box.axis[2])
        };
        temp = sphereVelocity - boxVelocity;
        Vector3<double> V
        {
            Dot(temp, box.axis[0]),
            Dot(temp, box.axis[1]),
            Dot(temp, box.axis[2])
        };
        FixDoQuery(box.extent, C, sphere.radius, V, result);
        if (result.intersectionType != 0)
        {
            auto& P = result.contactPoint;
            P = box.center + P[0] * box.axis[0] + P[1] * box.axis[1]
                + P[2] * box.axis[2];
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

    // Unrecorded draws for the dynamic oriented-box/sphere query.
    void RawOBox3Sphere3(oracle::Ctx& io, int mode, bool aimed,
        Vector3<double>& center, std::array<Vector3<double>, 3>& axis,
        Vector3<double>& extent, Vector3<double>& boxVelocity,
        Vector3<double>& sphereCenter, double& radius,
        Vector3<double>& sphereVelocity)
    {
        for (int i = 0; i < 3; ++i)
        {
            center[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-4.0, 4.0));
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
                axis[i].MakeZero();
                axis[i][perm[p][i]] = sign;
            }
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { axis[0][i] = io.raw(-1.0, 1.0); }
                len = Length(axis[0]);
            } while (len < 0.1 || len > 1.0);
            Normalize(axis[0]);
            ComputeOrthogonalComplement(1, axis.data());
        }
        for (int i = 0; i < 3; ++i)
        {
            extent[i] = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
                : io.raw(0.25, 3.0));
        }
        for (int i = 0; i < 3; ++i)
        {
            boxVelocity[i] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                : io.raw(-2.0, 2.0));
        }
        radius = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
            : io.raw(0.25, 3.0));
        if (aimed)
        {
            int face = io.rawInteger(0, 5);
            int i0 = face % 3, i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
            double s1 = 0.25 * static_cast<double>(io.rawInteger(-4, 4));
            double s2 = 0.25 * static_cast<double>(io.rawInteger(-4, 4));
            Vector3<double> local{};
            local[i0] = (face < 3 ? 1.0 : -1.0) * extent[i0];
            local[i1] = s1 * extent[i1];
            local[i2] = s2 * extent[i2];
            Vector3<double> target = center + local[0] * axis[0]
                + local[1] * axis[1] + local[2] * axis[2];
            Vector3<double> offset{};
            for (int i = 0; i < 3; ++i)
            {
                offset[i] = static_cast<double>(io.rawInteger(-4, 4));
            }
            sphereCenter = target + offset;
            double speed = 0.25 * static_cast<double>(io.rawInteger(1, 8));
            sphereVelocity = boxVelocity - speed * offset;
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                sphereCenter[i] = (mode == 0
                    ? static_cast<double>(io.rawInteger(-4, 4)) : io.raw(-4.0, 4.0));
            }
            for (int i = 0; i < 3; ++i)
            {
                sphereVelocity[i] = (mode == 0
                    ? static_cast<double>(io.rawInteger(-2, 2)) : io.raw(-2.0, 2.0));
            }
        }
    }

    void OBox3Sphere3Case(oracle::Ctx& io, int mode, bool aimed, bool requireDeviation)
    {
        Vector3<double> center{}, extent{}, boxVelocity{}, sphereCenter{},
            sphereVelocity{};
        std::array<Vector3<double>, 3> axis{};
        double radius{};
        for (;;)
        {
            RawOBox3Sphere3(io, mode, aimed, center, axis, extent, boxVelocity,
                sphereCenter, radius, sphereVelocity);
            OrientedBox3<double> probeBox{};
            probeBox.center = center;
            probeBox.axis = axis;
            probeBox.extent = extent;
            Sphere3<double> probeSphere{};
            probeSphere.center = sphereCenter;
            probeSphere.radius = radius;
            FIQuery<double, OrientedBox3<double>, Sphere3<double>> probeQuery;
            auto upstream = probeQuery(probeBox, boxVelocity, probeSphere, sphereVelocity);
            auto fixed = Ob3S3Fixed(probeBox, boxVelocity, probeSphere, sphereVelocity);
            Ab3S3Result up{};
            up.intersectionType = upstream.intersectionType;
            up.contactTime = upstream.contactTime;
            up.contactPoint = upstream.contactPoint;
            if (Ab3S3Same(up, fixed) != requireDeviation) { break; }
        }
        OrientedBox3<double> box{};
        box.center = io.givenVec(center);
        box.axis[0] = io.givenVec(axis[0]);
        box.axis[1] = io.givenVec(axis[1]);
        box.axis[2] = io.givenVec(axis[2]);
        box.extent = io.givenVec(extent);
        io.givenVec(boxVelocity);
        Sphere3<double> sphere{};
        sphere.center = io.givenVec(sphereCenter);
        sphere.radius = io.given(radius);
        io.givenVec(sphereVelocity);
        FIQuery<double, OrientedBox3<double>, Sphere3<double>> query;
        auto r = query(box, boxVelocity, sphere, sphereVelocity);
        io.outInt(r.intersectionType);
        io.outReal(r.contactTime);
        io.outVec(r.contactPoint);
    }
}

ORACLE_CASE("IntrOrientedBox3Sphere3.test")
{
    int mode = io.index() % 2;
    auto box = OBox3(io, mode, 3, 4.0);
    auto sphere = Sph<3>(io, mode, 3, 4.0);
    TIQuery<double, OrientedBox3<double>, Sphere3<double>> query;
    auto r = query(box, sphere);
    io.outBool(r.intersect);
}

// The generator rejects the inputs on which upstream's first-wins probe
// selection disagrees with the corrected one, that is exactly the inputs on
// which the port deliberately deviates (see the deviation case below).
ORACLE_CASE("IntrOrientedBox3Sphere3.find")
{
    OBox3Sphere3Case(io, io.index() % 2, false, false);
}

ORACLE_CASE("IntrOrientedBox3Sphere3.find.aimed")
{
    OBox3Sphere3Case(io, io.index() % 2, true, false);
}

// Deliberate deviation inherited from FIQuery<AlignedBox3,Sphere3>::DoQuery:
// the sphere enters the Minkowski sum through a rounded piece that upstream
// does not probe, or earlier than the piece upstream accepts.
// docs/UPSTREAM-FINDINGS.md IntrAlignedBox3Sphere3.h, issues #458 and #465.
ORACLE_CASE("IntrOrientedBox3Sphere3.find.probeDeviation")
{
    OBox3Sphere3Case(io, 1, true, true);
}

// ============================= IntrPlane3Circle3 =========================

namespace
{
    Circle3<double> Cir3(oracle::Ctx& io, int mode, int lat, double range)
    {
        Circle3<double> circle{};
        circle.center = Pt<3>(io, mode, lat, range);
        circle.normal = Dir<3>(io, mode);
        circle.radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        return circle;
    }

    // FIQuery<Plane3,Circle3>::Result::numIntersections is a size_t whose
    // "the whole circle" sentinel is SIZE_MAX, which is not representable as
    // a double; the port uses Number.MAX_SAFE_INTEGER. Both sides emit the
    // code 3 for that sentinel instead of the raw value.
    void OutNumIntersections(oracle::Ctx& io, size_t n)
    {
        io.outInt(n == std::numeric_limits<size_t>::max()
            ? 3 : static_cast<int32_t>(n));
    }

    void Plane3Circle3Find(oracle::Ctx& io, Plane3<double> const& plane,
        Circle3<double> const& circle)
    {
        FIQuery<double, Plane3<double>, Circle3<double>> query;
        auto r = query(plane, circle);
        io.outBool(r.intersect);
        OutNumIntersections(io, r.numIntersections);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
        io.outVec(r.circle.center);
        io.outVec(r.circle.normal);
        io.outReal(r.circle.radius);
    }
}

ORACLE_CASE("IntrPlane3Circle3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto circle = Cir3(io, mode, 3, 4.0);
    TIQuery<double, Plane3<double>, Circle3<double>> query;
    auto r = query(plane, circle);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrPlane3Circle3.find")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto circle = Cir3(io, mode, 3, 4.0);
    Plane3Circle3Find(io, plane, circle);
}

// The parallel branches of the delegated IntrPlane3Plane3 query (coplanar,
// hence the whole circle, and parallel but distinct, hence empty) need
// |Dot(N0,N1)| >= 1 exactly, which uniform normals essentially never satisfy.
// Both normals are the same signed coordinate axis here, so the dot product is
// exactly +-1, and the circle centre is on the plane on half the records.
ORACLE_CASE("IntrPlane3Circle3.find.parallel")
{
    int k = io.rawInteger(0, 2);
    Vector3<double> planeNormal{}, circleNormal{};
    planeNormal.MakeZero();
    planeNormal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    circleNormal.MakeZero();
    circleNormal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(planeNormal);
    Vector3<double> planeOrigin = io.latticeVec<3>(-3, 3);
    io.givenVec(circleNormal);
    Vector3<double> circleCenter = planeOrigin;
    for (int i = 0; i < 3; ++i)
    {
        circleCenter[i] += static_cast<double>(io.rawInteger(-3, 3))
            * (i == k ? 0.0 : 1.0);
    }
    if (io.rawInteger(0, 1) == 1) { circleCenter[k] += 1.0; }
    io.givenVec(circleCenter);
    double radius = io.lattice(1, 3);
    Plane3<double> plane(planeNormal, planeOrigin);
    Circle3<double> circle{};
    circle.center = circleCenter;
    circle.normal = circleNormal;
    circle.radius = radius;
    Plane3Circle3Find(io, plane, circle);
}

// The single-point branch (discriminant exactly zero, the plane tangent to the
// circle). The circle lies in a coordinate plane on the lattice and the input
// plane is perpendicular to it at signed distance 'radius + offset' from the
// circle centre along a coordinate axis, so the discriminant is an exact
// integer expression that is zero for offset = 0.
ORACLE_CASE("IntrPlane3Circle3.find.tangent")
{
    int k = io.rawInteger(0, 2);
    int j = (k + 1 + io.rawInteger(0, 1)) % 3;
    Vector3<double> circleNormal{}, planeNormal{};
    circleNormal.MakeZero();
    circleNormal[k] = 1.0;
    planeNormal.MakeZero();
    planeNormal[j] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(planeNormal);
    Vector3<double> circleCenter = io.latticeVec<3>(-3, 3);
    io.givenVec(circleNormal);
    double radius = io.lattice(1, 3);
    Vector3<double> planeOrigin = circleCenter;
    planeOrigin[j] += planeNormal[j] * (radius
        + static_cast<double>(io.rawInteger(-2, 1)));
    io.givenVec(planeOrigin);
    Plane3<double> plane(planeNormal, planeOrigin);
    Circle3<double> circle{};
    circle.center = circleCenter;
    circle.normal = circleNormal;
    circle.radius = radius;
    Plane3Circle3Find(io, plane, circle);
}

// ============================ IntrPlane3Cylinder3 ========================

namespace
{
    // A cylinder. An infinite cylinder is marked with height = -1, which is
    // upstream's sentinel.
    Cylinder3<double> Cyl3(oracle::Ctx& io, int mode, int lat, double range,
        bool infinite)
    {
        Cylinder3<double> cylinder{};
        cylinder.axis.origin = Pt<3>(io, mode, lat, range);
        cylinder.axis.direction = Dir<3>(io, mode);
        cylinder.radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        double height = (infinite ? -1.0
            : (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
                : io.raw(0.25, range)));
        cylinder.height = io.given(height);
        return cylinder;
    }

    void Plane3Cylinder3Find(oracle::Ctx& io, Plane3<double> const& plane,
        Cylinder3<double> const& cylinder)
    {
        FIQuery<double, Plane3<double>, Cylinder3<double>> query;
        auto r = query(plane, cylinder);
        io.outBool(r.intersect);
        io.outInt(static_cast<int32_t>(r.type));
        for (size_t i = 0; i < 2; ++i)
        {
            io.outVec(r.line[i].origin);
            io.outVec(r.line[i].direction);
        }
        io.outVec(r.ellipse.center);
        io.outVec(r.ellipse.normal);
        io.outVec(r.ellipse.axis[0]);
        io.outVec(r.ellipse.axis[1]);
        io.outReal(r.ellipse.extent[0]);
        io.outReal(r.ellipse.extent[1]);
        for (size_t i = 0; i < 2; ++i)
        {
            io.outVec(r.trimLine[i].origin);
            io.outVec(r.trimLine[i].direction);
        }
    }
}

ORACLE_CASE("IntrPlane3Cylinder3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto cylinder = Cyl3(io, mode, 3, 4.0, false);
    TIQuery<double, Plane3<double>, Cylinder3<double>> query;
    auto r = query(plane, cylinder);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrPlane3Cylinder3.test.infinite")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto cylinder = Cyl3(io, mode, 3, 4.0, true);
    TIQuery<double, Plane3<double>, Cylinder3<double>> query;
    auto r = query(plane, cylinder);
    io.outBool(r.intersect);
}

// The ELLIPSE and CIRCLE branches, which run Ellipse2::FromCoefficients (a
// symmetric eigensolve; arithmetic and sqrt only, no libm transcendental)
// and the two trim lines.
ORACLE_CASE("IntrPlane3Cylinder3.find")
{
    int mode = io.index() % 2;
    auto plane = Pln3(io, mode, 3, 4.0);
    auto cylinder = Cyl3(io, mode, 3, 4.0, false);
    Plane3Cylinder3Find(io, plane, cylinder);
}

// Dot(N,W) exactly zero reaches GetLinesOfIntersection, whose three branches
// (two lines, tangent single line, no intersection) are selected by the sign
// of r^2 - Dot(N,C-P)^2. Both the plane normal and the cylinder axis are
// coordinate axes here and the offset is an integer, so all three occur with
// exact arithmetic. Half the records use an infinite cylinder.
ORACLE_CASE("IntrPlane3Cylinder3.find.parallel")
{
    int k = io.rawInteger(0, 2);
    int j = (k + 1 + io.rawInteger(0, 1)) % 3;
    Vector3<double> normal{}, direction{};
    normal.MakeZero();
    normal[k] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    direction.MakeZero();
    direction[j] = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    io.givenVec(normal);
    Vector3<double> planeOrigin = io.latticeVec<3>(-3, 3);
    Vector3<double> center = planeOrigin;
    double radius = static_cast<double>(io.rawInteger(1, 3));
    for (int i = 0; i < 3; ++i)
    {
        center[i] += static_cast<double>(io.rawInteger(-3, 3))
            * (i == k ? 0.0 : 1.0);
    }
    center[k] += normal[k] * (radius + static_cast<double>(io.rawInteger(-3, 1)));
    io.givenVec(center);
    io.givenVec(direction);
    io.given(radius);
    double height = (io.rawInteger(0, 1) == 0 ? -1.0
        : static_cast<double>(io.rawInteger(1, 4)));
    io.given(height);
    Plane3<double> plane(normal, planeOrigin);
    Cylinder3<double> cylinder{};
    cylinder.axis.origin = center;
    cylinder.axis.direction = direction;
    cylinder.radius = radius;
    cylinder.height = height;
    Plane3Cylinder3Find(io, plane, cylinder);
}

// ============================ IntrRay2Arc2, IntrSegment2Arc2 =============

namespace
{
    // The twelve lattice points on the circle of radius 5, in counterclockwise
    // order starting at (5,0). An arc whose endpoints are two of these is
    // exactly on the circle, so the side-of-line test of Arc2::Contains is
    // evaluated in exact arithmetic.
    Vector2<double> const gCircle5[12] =
    {
        {  5.0,  0.0 }, {  4.0,  3.0 }, {  3.0,  4.0 }, {  0.0,  5.0 },
        { -3.0,  4.0 }, { -4.0,  3.0 }, { -5.0,  0.0 }, { -4.0, -3.0 },
        { -3.0, -4.0 }, {  0.0, -5.0 }, {  3.0, -4.0 }, {  4.0, -3.0 }
    };

    // An arc. Mode 0 puts the centre on the lattice, the radius at 5 and both
    // endpoints on gCircle5; mode 1 draws the endpoints from unrecorded
    // angles. Either way the recorded layout is centre, radius, end0, end1.
    Arc2<double> Arc(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector2<double> center{};
        double radius{};
        Vector2<double> e0{}, e1{};
        if (mode == 0)
        {
            for (int i = 0; i < 2; ++i)
            {
                center[i] = static_cast<double>(io.rawInteger(-lat, lat));
            }
            radius = 5.0;
            int i0 = io.rawInteger(0, 11);
            int i1 = (i0 + 1 + io.rawInteger(0, 10)) % 12;
            e0 = center + gCircle5[i0];
            e1 = center + gCircle5[i1];
        }
        else
        {
            for (int i = 0; i < 2; ++i) { center[i] = io.raw(-range, range); }
            radius = io.raw(0.5, range);
            double a0 = io.raw(0.0, 6.28318530717958647692);
            double a1 = a0 + io.raw(0.1, 6.0);
            e0 = center + radius * Vector2<double>{ std::cos(a0), std::sin(a0) };
            e1 = center + radius * Vector2<double>{ std::cos(a1), std::sin(a1) };
        }
        Arc2<double> arc{};
        arc.center = io.givenVec(center);
        arc.radius = io.given(radius);
        arc.end[0] = io.givenVec(e0);
        arc.end[1] = io.givenVec(e1);
        return arc;
    }

    using Ray2Arc2Result = FIQuery<double, Ray2<double>, Arc2<double>>::Result;
    using Seg2Arc2Result = FIQuery<double, Segment2<double>, Arc2<double>>::Result;

    // The answer the port reports: the line/circle CURVE intersection filtered
    // by the ray (or segment) parameter range, instead of upstream's
    // solid-disk clip, whose clipped endpoint is not on the circle at all
    // (docs/UPSTREAM-FINDINGS.md IntrRay2Arc2.h / IntrSegment2Arc2.h, issue
    // #304).
    template <typename R>
    R ArcFixed(Vector2<double> const& origin, Vector2<double> const& direction,
        Arc2<double> const& arc, double lo, double hi, bool isRay)
    {
        R result{};
        FIQuery<double, Line2<double>, Circle2<double>> lcQuery;
        Circle2<double> circle(arc.center, arc.radius);
        Line2<double> line(origin, direction);
        auto lc = lcQuery(line, circle);
        if (lc.intersect)
        {
            result.numIntersections = 0;
            for (int32_t i = 0; i < lc.numIntersections; ++i)
            {
                bool onComponent = (isRay ? lc.parameter[i] >= lo
                    : std::fabs(lc.parameter[i]) <= hi);
                if (onComponent && arc.Contains(lc.point[i]))
                {
                    result.intersect = true;
                    result.parameter[result.numIntersections] = lc.parameter[i];
                    result.point[result.numIntersections] = lc.point[i];
                    ++result.numIntersections;
                }
            }
        }
        return result;
    }

    template <typename R>
    bool SameArcResult(R const& a, R const& b)
    {
        if (a.intersect != b.intersect) { return false; }
        if (a.numIntersections != b.numIntersections) { return false; }
        for (int32_t i = 0; i < 2; ++i)
        {
            if (std::memcmp(&a.parameter[i], &b.parameter[i], sizeof(double)) != 0)
            {
                return false;
            }
            for (int32_t d = 0; d < 2; ++d)
            {
                if (std::memcmp(&a.point[i][d], &b.point[i][d], sizeof(double)) != 0)
                {
                    return false;
                }
            }
        }
        return true;
    }

    template <typename R>
    void OutArcResult(oracle::Ctx& io, R const& r)
    {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        io.outVec(r.point[0]);
        io.outVec(r.point[1]);
    }

    // Unrecorded draws of a ray aimed near the arc, so that the query reaches
    // the circle often. 'inside' places the origin strictly inside the disk,
    // which is where upstream's solid-disk clip corrupts the reported point.
    void RawRayForArc(oracle::Ctx& io, int mode, bool inside,
        Vector2<double> const& center, double radius, Vector2<double>& origin,
        Vector2<double>& direction)
    {
        if (mode == 0)
        {
            do
            {
                direction.MakeZero();
                int k = io.rawInteger(0, 3);
                direction[k % 2] = (k < 2 ? 1.0 : -1.0);
            } while (false);
            Vector2<double> target = center + gCircle5[io.rawInteger(0, 11)];
            double t = static_cast<double>(io.rawInteger(inside ? 1 : -6, 6));
            if (inside)
            {
                // Put the origin strictly inside the disk on the lattice.
                Vector2<double> offset{};
                offset[0] = static_cast<double>(io.rawInteger(-2, 2));
                offset[1] = static_cast<double>(io.rawInteger(-2, 2));
                origin = center + offset;
            }
            else
            {
                origin = target - t * direction;
            }
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 2; ++i) { direction[i] = io.raw(-1.0, 1.0); }
                len = Length(direction);
            } while (len < 0.1 || len > 1.0);
            Normalize(direction);
            if (inside)
            {
                double r = io.raw(0.0, 0.9) * radius;
                double a = io.raw(0.0, 6.28318530717958647692);
                origin = center + r * Vector2<double>{ std::cos(a), std::sin(a) };
            }
            else
            {
                double t = io.raw(-6.0, 6.0);
                double a = io.raw(0.0, 6.28318530717958647692);
                Vector2<double> target = center
                    + radius * Vector2<double>{ std::cos(a), std::sin(a) };
                origin = target - t * direction;
            }
        }
    }
}

ORACLE_CASE("IntrRay2Arc2.test")
{
    int mode = io.index() % 2;
    auto arc = Arc(io, mode, 3, 4.0);
    Vector2<double> origin{}, direction{};
    for (;;)
    {
        RawRayForArc(io, mode, false, arc.center, arc.radius, origin, direction);
        auto upstream = FIQuery<double, Ray2<double>, Arc2<double>>()(
            Ray2<double>(origin, direction), arc);
        auto fixed = ArcFixed<Ray2Arc2Result>(origin, direction, arc, 0.0, 0.0, true);
        if (SameArcResult(upstream, fixed)) { break; }
    }
    io.givenVec(origin);
    io.givenVec(direction);
    TIQuery<double, Ray2<double>, Arc2<double>> query;
    auto r = query(Ray2<double>(origin, direction), arc);
    io.outBool(r.intersect);
}

// The generator rejects the rays on which upstream's solid-disk clip replaces
// a line-circle root by the ray origin, which is exactly where the port
// deviates (see the deviation case below).
ORACLE_CASE("IntrRay2Arc2.find")
{
    int mode = io.index() % 2;
    auto arc = Arc(io, mode, 3, 4.0);
    Vector2<double> origin{}, direction{};
    for (;;)
    {
        RawRayForArc(io, mode, false, arc.center, arc.radius, origin, direction);
        auto upstream = FIQuery<double, Ray2<double>, Arc2<double>>()(
            Ray2<double>(origin, direction), arc);
        auto fixed = ArcFixed<Ray2Arc2Result>(origin, direction, arc, 0.0, 0.0, true);
        if (SameArcResult(upstream, fixed)) { break; }
    }
    io.givenVec(origin);
    io.givenVec(direction);
    FIQuery<double, Ray2<double>, Arc2<double>> query;
    auto r = query(Ray2<double>(origin, direction), arc);
    OutArcResult(io, r);
}

// Deliberate deviation: the ray origin is inside the disk, so upstream's
// ray/circle query clips the near root to the ray origin and hands that
// interior point to Arc2::Contains, reporting a "hit" at a point that is not
// on the circle. docs/UPSTREAM-FINDINGS.md IntrRay2Arc2.h, issue #304.
ORACLE_CASE("IntrRay2Arc2.find.insideDeviation")
{
    int mode = io.index() % 2;
    auto arc = Arc(io, mode, 3, 4.0);
    Vector2<double> origin{}, direction{};
    RawRayForArc(io, mode, true, arc.center, arc.radius, origin, direction);
    io.givenVec(origin);
    io.givenVec(direction);
    FIQuery<double, Ray2<double>, Arc2<double>> query;
    auto r = query(Ray2<double>(origin, direction), arc);
    OutArcResult(io, r);
}

ORACLE_CASE("IntrSegment2Arc2.test")
{
    int mode = io.index() % 2;
    auto arc = Arc(io, mode, 3, 4.0);
    Vector2<double> p0{}, p1{};
    Vector2<double> segOrigin{}, segDirection{};
    double segExtent{};
    for (;;)
    {
        Vector2<double> origin{}, direction{};
        RawRayForArc(io, mode, false, arc.center, arc.radius, origin, direction);
        double half = (mode == 0 ? static_cast<double>(io.rawInteger(1, 6))
            : io.raw(0.5, 6.0));
        p0 = origin - half * direction;
        p1 = origin + half * direction;
        Segment2<double> probe(p0, p1);
        probe.GetCenteredForm(segOrigin, segDirection, segExtent);
        auto upstream = FIQuery<double, Segment2<double>, Arc2<double>>()(probe, arc);
        auto fixed = ArcFixed<Seg2Arc2Result>(segOrigin, segDirection, arc,
            0.0, segExtent, false);
        if (SameArcResult(upstream, fixed)) { break; }
    }
    io.givenVec(p0);
    io.givenVec(p1);
    TIQuery<double, Segment2<double>, Arc2<double>> query;
    auto r = query(Segment2<double>(p0, p1), arc);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment2Arc2.find")
{
    int mode = io.index() % 2;
    auto arc = Arc(io, mode, 3, 4.0);
    Vector2<double> p0{}, p1{};
    Vector2<double> segOrigin{}, segDirection{};
    double segExtent{};
    for (;;)
    {
        Vector2<double> origin{}, direction{};
        RawRayForArc(io, mode, false, arc.center, arc.radius, origin, direction);
        double half = (mode == 0 ? static_cast<double>(io.rawInteger(1, 6))
            : io.raw(0.5, 6.0));
        p0 = origin - half * direction;
        p1 = origin + half * direction;
        Segment2<double> probe(p0, p1);
        probe.GetCenteredForm(segOrigin, segDirection, segExtent);
        auto upstream = FIQuery<double, Segment2<double>, Arc2<double>>()(probe, arc);
        auto fixed = ArcFixed<Seg2Arc2Result>(segOrigin, segDirection, arc,
            0.0, segExtent, false);
        if (SameArcResult(upstream, fixed)) { break; }
    }
    io.givenVec(p0);
    io.givenVec(p1);
    FIQuery<double, Segment2<double>, Arc2<double>> query;
    auto r = query(Segment2<double>(p0, p1), arc);
    OutArcResult(io, r);
}

// Deliberate deviation: a segment endpoint is inside the disk, so upstream's
// segment/circle query clips that root to the endpoint and hands the interior
// point to Arc2::Contains. docs/UPSTREAM-FINDINGS.md IntrSegment2Arc2.h,
// issue #304.
ORACLE_CASE("IntrSegment2Arc2.find.insideDeviation")
{
    int mode = io.index() % 2;
    auto arc = Arc(io, mode, 3, 4.0);
    Vector2<double> origin{}, direction{};
    RawRayForArc(io, mode, true, arc.center, arc.radius, origin, direction);
    double half = (mode == 0 ? static_cast<double>(io.rawInteger(4, 9))
        : io.raw(3.0, 9.0));
    Vector2<double> p0 = origin;
    Vector2<double> p1 = origin + half * direction;
    io.givenVec(p0);
    io.givenVec(p1);
    FIQuery<double, Segment2<double>, Arc2<double>> query;
    auto r = query(Segment2<double>(p0, p1), arc);
    OutArcResult(io, r);
}

// ================= IntrRay2SegmentMesh2, IntrSegment2SegmentMesh2 ========
//
// Both queries delegate to FIQuery<Line2,SegmentMesh2>, which orders its
// output with std::sort on the line parameter alone. std::sort is not stable
// and Array.prototype.sort is, so records with equal line parameters (a line
// through a shared vertex of a contiguous polysegment, which the aimed
// generator produces on purpose) may come out in either order. Both sides
// therefore re-sort the reported intersections by the TOTAL key
// (parameter, indexPair[0], indexPair[1], meshSegmentParameter, point) before
// emitting them.

namespace
{
    // A polysegment mesh. Mode 0 puts the vertices on a small lattice; the
    // topology flag chooses the open (numVertices-1 segments) or closed
    // (numVertices segments) contiguous constructor. The recorded layout is
    // the vertex count, the vertices, then the topology flag.
    SegmentMesh2<double> Mesh2(oracle::Ctx& io, int mode,
        std::vector<Vector2<double>>& vertices)
    {
        int numVertices = io.integer(3, 6);
        vertices.resize(static_cast<size_t>(numVertices));
        for (int i = 0; i < numVertices; ++i)
        {
            for (int d = 0; d < 2; ++d)
            {
                vertices[static_cast<size_t>(i)][d] = (mode == 0
                    ? io.lattice(-4, 4) : io.real(-5.0, 5.0));
            }
        }
        bool isOpen = io.boolean();
        return SegmentMesh2<double>(vertices, isOpen);
    }

    struct MeshHit
    {
        double parameter;
        double meshSegmentParameter;
        std::array<size_t, 2> indexPair;
        Vector2<double> point;
    };

    bool MeshHitLess(MeshHit const& a, MeshHit const& b)
    {
        if (a.parameter != b.parameter) { return a.parameter < b.parameter; }
        if (a.indexPair[0] != b.indexPair[0]) { return a.indexPair[0] < b.indexPair[0]; }
        if (a.indexPair[1] != b.indexPair[1]) { return a.indexPair[1] < b.indexPair[1]; }
        if (a.meshSegmentParameter != b.meshSegmentParameter)
        {
            return a.meshSegmentParameter < b.meshSegmentParameter;
        }
        if (a.point[0] != b.point[0]) { return a.point[0] < b.point[0]; }
        return a.point[1] < b.point[1];
    }

    void OutMeshHits(oracle::Ctx& io, std::vector<MeshHit>& hits)
    {
        std::sort(hits.begin(), hits.end(), MeshHitLess);
        io.outInt(static_cast<int32_t>(hits.size()));
        for (auto const& hit : hits)
        {
            io.outInt(static_cast<int32_t>(hit.indexPair[0]));
            io.outInt(static_cast<int32_t>(hit.indexPair[1]));
            io.outReal(hit.parameter);
            io.outReal(hit.meshSegmentParameter);
            io.outVec(hit.point);
        }
    }

    // Unrecorded draws of a line aimed at a mesh vertex (mode 0, integer
    // direction, so the line passes through the vertex exactly and the two
    // incident segments report the same line parameter) or at a uniform
    // point.
    void RawLineForMesh(oracle::Ctx& io, int mode,
        std::vector<Vector2<double>> const& vertices, bool aimed,
        Vector2<double>& origin, Vector2<double>& direction)
    {
        if (mode == 0 && aimed)
        {
            do
            {
                for (int d = 0; d < 2; ++d)
                {
                    direction[d] = static_cast<double>(io.rawInteger(-2, 2));
                }
            } while (direction == Vector2<double>::Zero());
            Vector2<double> target =
                vertices[static_cast<size_t>(io.rawInteger(0,
                    static_cast<int>(vertices.size()) - 1))];
            double t = 0.5 * static_cast<double>(io.rawInteger(-4, 4));
            origin = target - t * direction;
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int d = 0; d < 2; ++d) { direction[d] = io.raw(-1.0, 1.0); }
                len = Length(direction);
            } while (len < 0.1 || len > 1.0);
            Normalize(direction);
            for (int d = 0; d < 2; ++d) { origin[d] = io.raw(-5.0, 5.0); }
        }
    }

    void Ray2MeshCase(oracle::Ctx& io, int mode, bool aimed)
    {
        std::vector<Vector2<double>> vertices;
        auto mesh = Mesh2(io, mode, vertices);
        Vector2<double> origin{}, direction{};
        RawLineForMesh(io, mode, vertices, aimed, origin, direction);
        io.givenVec(origin);
        io.givenVec(direction);
        FIQuery<double, Ray2<double>, SegmentMesh2<double>> query;
        auto r = query(Ray2<double>(origin, direction), mesh);
        std::vector<MeshHit> hits;
        for (auto const& object : r.intersections)
        {
            hits.push_back({ object.rayParameter, object.meshSegmentParameter,
                object.indexPair, object.point });
        }
        OutMeshHits(io, hits);
    }

    void Segment2MeshCase(oracle::Ctx& io, int mode, bool aimed)
    {
        std::vector<Vector2<double>> vertices;
        auto mesh = Mesh2(io, mode, vertices);
        Vector2<double> origin{}, direction{};
        RawLineForMesh(io, mode, vertices, aimed, origin, direction);
        double half = (mode == 0 ? static_cast<double>(io.rawInteger(1, 5))
            : io.raw(1.0, 6.0));
        Vector2<double> p0 = origin;
        Vector2<double> p1 = origin + half * direction;
        io.givenVec(p0);
        io.givenVec(p1);
        FIQuery<double, Segment2<double>, SegmentMesh2<double>> query;
        auto r = query(Segment2<double>(p0, p1), mesh);
        std::vector<MeshHit> hits;
        for (auto const& object : r.intersections)
        {
            hits.push_back({ object.segmentParameter, object.meshSegmentParameter,
                object.indexPair, object.point });
        }
        OutMeshHits(io, hits);
    }
}

ORACLE_CASE("IntrRay2SegmentMesh2.find")
{
    Ray2MeshCase(io, io.index() % 2, false);
}

ORACLE_CASE("IntrRay2SegmentMesh2.find.throughVertex")
{
    Ray2MeshCase(io, 0, true);
}

ORACLE_CASE("IntrSegment2SegmentMesh2.find")
{
    Segment2MeshCase(io, io.index() % 2, false);
}

ORACLE_CASE("IntrSegment2SegmentMesh2.find.throughVertex")
{
    Segment2MeshCase(io, 0, true);
}

// ================= IntrRay3Ellipsoid3, IntrSegment3Ellipsoid3 ============

namespace
{
    using Line3Ell3Result = FIQuery<double, Line3<double>, Ellipsoid3<double>>::Result;

    void OutEllipsoidResult(oracle::Ctx& io, Line3Ell3Result const& r,
        bool withPoints)
    {
        io.outBool(r.intersect);
        io.outInt(static_cast<int32_t>(r.numIntersections));
        io.outReal(r.parameter[0]);
        io.outReal(r.parameter[1]);
        if (withPoints)
        {
            io.outVec(r.point[0]);
            io.outVec(r.point[1]);
        }
    }

    // An axis-aligned ellipsoid on the lattice whose extents are powers of
    // two, so that M = diag(1/e_i^2) is exact. Used by the aimed cases, where
    // the discriminant must be exactly zero for the tangency branch.
    Ellipsoid3<double> DyadicEllipsoid(oracle::Ctx& io)
    {
        Ellipsoid3<double> ellipsoid{};
        ellipsoid.center = io.latticeVec<3>(-3, 3);
        std::array<Vector3<double>, 3> axis{};
        for (int i = 0; i < 3; ++i) { axis[i].MakeZero(); axis[i][i] = 1.0; }
        io.givenVec(axis[0]);
        io.givenVec(axis[1]);
        io.givenVec(axis[2]);
        ellipsoid.axis = axis;
        Vector3<double> extent{};
        for (int i = 0; i < 3; ++i)
        {
            extent[i] = static_cast<double>(1 << io.rawInteger(0, 2));
        }
        ellipsoid.extent = io.givenVec(extent);
        return ellipsoid;
    }

    // A line through the point C + (e_i + d) * e_i with direction e_j, whose
    // discriminant is exactly zero for d = 0 (tangency) and whose near root is
    // at the recorded parameter offset, so the [0,+infinity) and [-e,e] clips
    // are evaluated at exact equality.
    void RawTangentLine(oracle::Ctx& io, Ellipsoid3<double> const& ellipsoid,
        Vector3<double>& origin, Vector3<double>& direction, int& i, int& j)
    {
        i = io.rawInteger(0, 2);
        j = (i + 1 + io.rawInteger(0, 1)) % 3;
        double d = static_cast<double>(io.rawInteger(-1, 1));
        double s = static_cast<double>(io.rawInteger(-3, 3));
        Vector3<double> ei{}, ej{};
        ei.MakeUnit(i);
        ej.MakeUnit(j);
        origin = ellipsoid.center + (ellipsoid.extent[i] + d) * ei + s * ej;
        direction = ej;
    }
}

namespace
{
    // Unrecorded draws of a linear component aimed at the neighbourhood of
    // the ellipsoid. A uniform origin and direction hit a small ellipsoid on
    // about 5 % of records, which leaves the two-root and one-root branches
    // nearly untested; here the component is aimed at a point at relative
    // radius 0.5 to 1.5 in the ellipsoid's own frame, so about half the
    // records hit.
    void RawLineForEllipsoid(oracle::Ctx& io, int mode,
        Ellipsoid3<double> const& ellipsoid, Vector3<double>& origin,
        Vector3<double>& direction)
    {
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
                for (int d = 0; d < 3; ++d) { direction[d] = io.raw(-1.0, 1.0); }
                len = Length(direction);
            } while (len < 0.1 || len > 1.0);
            Normalize(direction);
        }
        Vector3<double> u{};
        double ulen = 0.0;
        do
        {
            for (int d = 0; d < 3; ++d) { u[d] = io.raw(-1.0, 1.0); }
            ulen = Length(u);
        } while (ulen < 0.1 || ulen > 1.0);
        Normalize(u);
        double s = io.raw(0.5, 1.5);
        Vector3<double> target = ellipsoid.center;
        for (int d = 0; d < 3; ++d)
        {
            target += (s * ellipsoid.extent[d] * u[d]) * ellipsoid.axis[d];
        }
        double t = (mode == 0 ? static_cast<double>(io.rawInteger(-6, 6))
            : io.raw(-6.0, 6.0));
        origin = target - t * direction;
    }
}

ORACLE_CASE("IntrRay3Ellipsoid3.test")
{
    int mode = io.index() % 2;
    auto ellipsoid = Ellipsoid(io, mode, 3, 4.0);
    Vector3<double> origin{}, direction{};
    RawLineForEllipsoid(io, mode, ellipsoid, origin, direction);
    io.givenVec(origin);
    io.givenVec(direction);
    TIQuery<double, Ray3<double>, Ellipsoid3<double>> query;
    auto r = query(Ray3<double>(origin, direction), ellipsoid);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3Ellipsoid3.find")
{
    int mode = io.index() % 2;
    auto ellipsoid = Ellipsoid(io, mode, 3, 4.0);
    Vector3<double> origin{}, direction{};
    RawLineForEllipsoid(io, mode, ellipsoid, origin, direction);
    io.givenVec(origin);
    io.givenVec(direction);
    FIQuery<double, Ray3<double>, Ellipsoid3<double>> query;
    auto r = query(Ray3<double>(origin, direction), ellipsoid);
    OutEllipsoidResult(io, r, true);
}

ORACLE_CASE("IntrRay3Ellipsoid3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto ellipsoid = Ellipsoid(io, mode, 3, 4.0);
    Vector3<double> origin{}, direction{};
    RawLineForEllipsoid(io, mode, ellipsoid, origin, direction);
    io.givenVec(origin);
    io.givenVec(direction);
    Expose<FIQuery<double, Ray3<double>, Ellipsoid3<double>>> query;
    FIQuery<double, Ray3<double>, Ellipsoid3<double>>::Result r{};
    query.DoQuery(origin, direction, ellipsoid, r);
    OutEllipsoidResult(io, r, false);
}

ORACLE_CASE("IntrRay3Ellipsoid3.find.tangent")
{
    auto ellipsoid = DyadicEllipsoid(io);
    Vector3<double> origin{}, direction{};
    int i = 0, j = 0;
    RawTangentLine(io, ellipsoid, origin, direction, i, j);
    io.givenVec(origin);
    io.givenVec(direction);
    FIQuery<double, Ray3<double>, Ellipsoid3<double>> query;
    auto r = query(Ray3<double>(origin, direction), ellipsoid);
    OutEllipsoidResult(io, r, true);
}

namespace
{
    // Upstream's IntrSegment3Ellipsoid3 TIQuery reports "no intersection"
    // when Q(-e) and Q(e) are both negative, that is when the whole segment
    // lies strictly inside the solid ellipsoid; the port reports an
    // intersection (docs/UPSTREAM-FINDINGS.md IntrSegment3Ellipsoid3.h
    // TIQuery, issue #304). Upstream's other branches agree with the port, so
    // this predicate is exactly the set on which the port deviates. It is
    // evaluated on upstream's own expressions, after upstream's discr < 0
    // early exit and on the centred form the query itself uses.
    bool Seg3Ell3StrictlyInside(Vector3<double> const& p0,
        Vector3<double> const& p1, Ellipsoid3<double> const& ellipsoid)
    {
        Segment3<double> segment(p0, p1);
        Vector3<double> segOrigin{}, segDirection{};
        double segExtent{};
        segment.GetCenteredForm(segOrigin, segDirection, segExtent);
        Matrix3x3<double> M{};
        ellipsoid.GetM(M);
        Vector3<double> diff = segOrigin - ellipsoid.center;
        Vector3<double> matDiff = M * diff;
        Vector3<double> matDir = M * segDirection;
        double a0 = Dot(diff, matDiff) - 1.0;
        double a1 = Dot(segDirection, matDiff);
        double a2 = Dot(segDirection, matDir);
        double discr = a1 * a1 - a0 * a2;
        if (discr < 0.0) { return false; }
        double a2e = a2 * segExtent;
        double tmp0 = a2e * segExtent + a0;
        double tmp1 = 2.0 * a1 * segExtent;
        return (tmp0 - tmp1) < 0.0 && (tmp0 + tmp1) < 0.0;
    }

    void RawSegForEllipsoid(oracle::Ctx& io, int mode, bool inside,
        Ellipsoid3<double> const& ellipsoid, Vector3<double>& p0,
        Vector3<double>& p1)
    {
        if (inside)
        {
            // Both endpoints strictly inside the solid ellipsoid: draw a
            // direction and a radius under 1/2 in the ellipsoid's own frame.
            for (int k = 0; k < 2; ++k)
            {
                Vector3<double> u{};
                double len = 0.0;
                do
                {
                    for (int d = 0; d < 3; ++d) { u[d] = io.raw(-1.0, 1.0); }
                    len = Length(u);
                } while (len < 0.1 || len > 1.0);
                Normalize(u);
                double s = io.raw(0.0, 0.5);
                Vector3<double>& p = (k == 0 ? p0 : p1);
                p = ellipsoid.center;
                for (int d = 0; d < 3; ++d)
                {
                    p += (s * ellipsoid.extent[d] * u[d]) * ellipsoid.axis[d];
                }
            }
        }
        else
        {
            for (int d = 0; d < 3; ++d)
            {
                p0[d] = (mode == 0 ? static_cast<double>(io.rawInteger(-4, 4))
                    : io.raw(-5.0, 5.0));
                p1[d] = (mode == 0 ? static_cast<double>(io.rawInteger(-4, 4))
                    : io.raw(-5.0, 5.0));
            }
        }
    }
}

// The generator rejects the segments that lie strictly inside the ellipsoid,
// which is exactly where the port deviates (see the deviation case below).
ORACLE_CASE("IntrSegment3Ellipsoid3.test")
{
    int mode = io.index() % 2;
    auto ellipsoid = Ellipsoid(io, mode, 3, 4.0);
    Vector3<double> p0{}, p1{};
    for (;;)
    {
        RawSegForEllipsoid(io, mode, false, ellipsoid, p0, p1);
        if (!Seg3Ell3StrictlyInside(p0, p1, ellipsoid)) { break; }
    }
    io.givenVec(p0);
    io.givenVec(p1);
    TIQuery<double, Segment3<double>, Ellipsoid3<double>> query;
    auto r = query(Segment3<double>(p0, p1), ellipsoid);
    io.outBool(r.intersect);
}

// Deliberate deviation: a segment strictly inside the solid ellipsoid, which
// upstream reports as not intersecting, contradicting its own FI query.
// docs/UPSTREAM-FINDINGS.md IntrSegment3Ellipsoid3.h TIQuery, issue #304.
ORACLE_CASE("IntrSegment3Ellipsoid3.test.containedDeviation")
{
    int mode = io.index() % 2;
    auto ellipsoid = Ellipsoid(io, mode, 3, 4.0);
    Vector3<double> p0{}, p1{};
    RawSegForEllipsoid(io, mode, true, ellipsoid, p0, p1);
    io.givenVec(p0);
    io.givenVec(p1);
    TIQuery<double, Segment3<double>, Ellipsoid3<double>> query;
    auto r = query(Segment3<double>(p0, p1), ellipsoid);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3Ellipsoid3.find")
{
    int mode = io.index() % 2;
    auto ellipsoid = Ellipsoid(io, mode, 3, 4.0);
    Vector3<double> origin{}, direction{};
    RawLineForEllipsoid(io, mode, ellipsoid, origin, direction);
    double half = (mode == 0 ? static_cast<double>(io.rawInteger(1, 6))
        : io.raw(0.5, 6.0));
    Vector3<double> p0 = origin - half * direction;
    Vector3<double> p1 = origin + half * direction;
    io.givenVec(p0);
    io.givenVec(p1);
    FIQuery<double, Segment3<double>, Ellipsoid3<double>> query;
    auto r = query(Segment3<double>(p0, p1), ellipsoid);
    OutEllipsoidResult(io, r, true);
}

ORACLE_CASE("IntrSegment3Ellipsoid3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto ellipsoid = Ellipsoid(io, mode, 3, 4.0);
    auto segOrigin = Pt<3>(io, mode, 3, 4.0);
    auto segDirection = Dir<3>(io, mode);
    double segExtent = (mode == 0 ? io.lattice(1, 4) : io.real(0.5, 5.0));
    Expose<FIQuery<double, Segment3<double>, Ellipsoid3<double>>> query;
    FIQuery<double, Segment3<double>, Ellipsoid3<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, ellipsoid, r);
    OutEllipsoidResult(io, r, false);
}

ORACLE_CASE("IntrSegment3Ellipsoid3.find.tangent")
{
    auto ellipsoid = DyadicEllipsoid(io);
    Vector3<double> origin{}, direction{};
    int i = 0, j = 0;
    RawTangentLine(io, ellipsoid, origin, direction, i, j);
    double half = static_cast<double>(io.rawInteger(1, 5));
    Vector3<double> p0 = origin - half * direction;
    Vector3<double> p1 = origin + half * direction;
    io.givenVec(p0);
    io.givenVec(p1);
    FIQuery<double, Segment3<double>, Ellipsoid3<double>> query;
    auto r = query(Segment3<double>(p0, p1), ellipsoid);
    OutEllipsoidResult(io, r, true);
}

// ============================== IntrLine3Cone3 ===========================

namespace
{
    using L3C3QFN = QFNumber<double, 1>;
    using L3C3Result = FIQuery<double, Line3<double>, Cone3<double>>::Result;

    void OutQFN(oracle::Ctx& io, L3C3QFN const& q)
    {
        io.outReal(q.x[0]);
        io.outReal(q.x[1]);
        io.outReal(q.d);
    }

    void OutL3C3(oracle::Ctx& io, L3C3Result const& r, bool withPoints)
    {
        io.outBool(r.intersect);
        io.outInt(r.type);
        OutQFN(io, r.t[0]);
        OutQFN(io, r.t[1]);
        if (withPoints)
        {
            for (int32_t i = 0; i < 2; ++i)
            {
                for (int32_t k = 0; k < 3; ++k) { OutQFN(io, r.P[i][k]); }
            }
        }
    }

    // True when the port's correction of upstream's through-vertex handling
    // can change the answer (docs/UPSTREAM-FINDINGS.md IntrLine3Cone3.h,
    // issues #304 and #465; the port's fix is in four places, all of them
    // covered here). The predicate follows upstream's own control flow, on
    // the direction DoQuery actually uses:
    //   1. the line contains the cone vertex exactly, where the port answers
    //      from the sign of c2 alone instead of from the rounded
    //      discriminant (a superset of the deviating inputs: upstream
    //      sometimes agrees);
    //   2. discr < 0 with c2 > 0, which cannot happen in exact arithmetic and
    //      is the rounding of a double root;
    //   3. discr > 0 with h[0] < 0 < h[1] and c2 < 0, where upstream's
    //      SetRayClamp assumes c2 > 0;
    //   4. discr == 0 with c2 > 0 away from the vertex, which is again a
    //      rounded double root.
    bool Line3Cone3PortDeviates(Vector3<double> const& lineOrigin,
        Vector3<double> const& inU, Cone3<double> const& cone)
    {
        Vector3<double> U =
            (Dot(inU, cone.ray.direction) >= 0.0 ? inU : -inU);
        Vector3<double> PmV = lineOrigin - cone.ray.origin;
        double UdU = Dot(U, U);
        double DdU = Dot(cone.ray.direction, U);
        double DdPmV = Dot(cone.ray.direction, PmV);
        double UdPmV = Dot(U, PmV);
        double PmVdPmV = Dot(PmV, PmV);
        double c2 = DdU * DdU - cone.cosAngleSqr * UdU;
        double c1 = DdU * DdPmV - cone.cosAngleSqr * UdPmV;
        double c0 = DdPmV * DdPmV - cone.cosAngleSqr * PmVdPmV;

        double tv = -UdPmV / UdU;
        bool onVertex = true;
        for (int32_t i = 0; i < 3; ++i)
        {
            if (PmV[i] + tv * U[i] != 0.0) { onVertex = false; break; }
        }
        if (onVertex) { return true; }

        if (c2 != 0.0)
        {
            double discr = c1 * c1 - c0 * c2;
            if (discr < 0.0) { return c2 > 0.0; }
            if (discr > 0.0)
            {
                double x = -c1 / c2;
                double y = (c2 > 0.0 ? 1.0 / c2 : -1.0 / c2);
                std::array<L3C3QFN, 2> t = { L3C3QFN(x, -y, discr), L3C3QFN(x, y, discr) };
                std::array<L3C3QFN, 2> h = { t[0] * DdU + DdPmV, t[1] * DdU + DdPmV };
                L3C3QFN zero(0, 0, discr);
                if (h[0] >= zero) { return false; }
                if (h[1] <= zero) { return false; }
                return c2 < 0.0;
            }
            return c2 > 0.0;
        }
        return false;
    }

    // Unrecorded draws of a line for the cone query. 'throughVertex' places
    // the line origin at V - t*U with an integer direction U, so that the
    // vertex is on the line in exact arithmetic.
    void RawLineForCone(oracle::Ctx& io, int mode, bool throughVertex,
        Cone3<double> const& cone, Vector3<double>& origin,
        Vector3<double>& direction)
    {
        if (mode == 0 || throughVertex)
        {
            do
            {
                for (int d = 0; d < 3; ++d)
                {
                    direction[d] = static_cast<double>(io.rawInteger(-3, 3));
                }
            } while (direction == Vector3<double>::Zero());
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int d = 0; d < 3; ++d) { direction[d] = io.raw(-1.0, 1.0); }
                len = Length(direction);
            } while (len < 0.1 || len > 1.0);
            Normalize(direction);
        }
        if (throughVertex)
        {
            double t = static_cast<double>(io.rawInteger(-4, 4));
            origin = cone.ray.origin - t * direction;
        }
        else
        {
            for (int d = 0; d < 3; ++d)
            {
                origin[d] = (mode == 0 ? static_cast<double>(io.rawInteger(-4, 4))
                    : io.raw(-5.0, 5.0));
            }
        }
    }

    void Line3Cone3Case(oracle::Ctx& io, int mode, int kind, bool throughVertex,
        bool doQueryOnly)
    {
        auto cone = Cn<3>(io, mode, 3, 4.0, kind);
        Vector3<double> origin{}, direction{};
        for (;;)
        {
            RawLineForCone(io, mode, throughVertex, cone, origin, direction);
            if (Line3Cone3PortDeviates(origin, direction, cone) == throughVertex)
            {
                break;
            }
        }
        io.givenVec(origin);
        io.givenVec(direction);
        Line3<double> line(origin, direction);
        if (doQueryOnly)
        {
            Expose<FIQuery<double, Line3<double>, Cone3<double>>> query;
            L3C3Result r{};
            query.DoQuery(line.origin, line.direction, cone, r);
            OutL3C3(io, r, false);
        }
        else
        {
            FIQuery<double, Line3<double>, Cone3<double>> query;
            auto r = query(line, cone);
            OutL3C3(io, r, true);
        }
    }
}

// The generator rejects the lines on which the port's through-vertex
// correction can change the answer (see the deviation case below). The four
// cone kinds are visited in turn: infinite, infinite truncated, finite and
// frustum.
ORACLE_CASE("IntrLine3Cone3.find")
{
    Line3Cone3Case(io, io.index() % 2, io.index() % 4, false, false);
}

ORACLE_CASE("IntrLine3Cone3.doQuery.fi")
{
    Line3Cone3Case(io, io.index() % 2, io.index() % 4, false, true);
}

// Deliberate deviation: the line contains the cone vertex exactly, so the
// quadratic has a double root and upstream's case analysis is decided by the
// sign of a discriminant that has no significant digits. Upstream reports a
// point, the empty set, or a segment reaching the cone's maximum height where
// the answer is a ray or the vertex alone. docs/UPSTREAM-FINDINGS.md
// IntrLine3Cone3.h, issues #304 and #465.
ORACLE_CASE("IntrLine3Cone3.find.vertexDeviation")
{
    Line3Cone3Case(io, 0, io.index() % 4, true, false);
}

// ============================ IntrAlignedBox3Cone3 =======================

namespace
{
    using Box3Cone3Query = TIQuery<double, AlignedBox3<double>, Cone3<double>>;

    // The port exports the three protected static helpers as module
    // functions, so they are called directly on both sides.
    struct ExposeBox3Cone3 : public Box3Cone3Query
    {
        using Box3Cone3Query::ComputeBoxHeightInterval;
        using Box3Cone3Query::ConeAxisIntersectsBox;
        using Box3Cone3Query::HasPointInsideCone;
    };

    // An axis-aligned box with min <= max, componentwise.
    AlignedBox3<double> ABox3(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector3<double> lo{}, hi{};
        for (int i = 0; i < 3; ++i)
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
        AlignedBox3<double> box{};
        box.min = io.givenVec(lo);
        box.max = io.givenVec(hi);
        return box;
    }

    // A box that straddles a cone height plane, so that the query takes the
    // clipping path (ComputeCandidatesOnBoxEdges and
    // ComputeCandidatesOnBoxFaces) rather than the quick-reject, the
    // quick-accept or the fully-in-slab path.
    AlignedBox3<double> RawStraddlingBox(oracle::Ctx& io,
        Cone3<double> const& cone, Vector3<double>& lo, Vector3<double>& hi)
    {
        double h = (io.rawInteger(0, 1) == 0 ? cone.GetMinHeight()
            : (cone.IsFinite() ? cone.GetMaxHeight() : cone.GetMinHeight()));
        Vector3<double> center = cone.ray.origin + h * cone.ray.direction;
        for (int i = 0; i < 3; ++i)
        {
            center[i] += static_cast<double>(io.rawInteger(-2, 2));
        }
        for (int i = 0; i < 3; ++i)
        {
            double e = static_cast<double>(io.rawInteger(1, 3));
            lo[i] = center[i] - e;
            hi[i] = center[i] + e;
        }
        AlignedBox3<double> box{};
        box.min = lo;
        box.max = hi;
        return box;
    }

    // A box entirely between the two height planes, which is the path that
    // leaves stale adjacency bits behind in upstream's query object.
    AlignedBox3<double> RawInSlabBox(oracle::Ctx& io, Cone3<double> const& cone,
        Vector3<double>& lo, Vector3<double>& hi)
    {
        double hmin = cone.GetMinHeight();
        double hmax = (cone.IsFinite() ? cone.GetMaxHeight() : hmin + 8.0);
        double mid = 0.5 * (hmin + hmax);
        Vector3<double> center = cone.ray.origin + mid * cone.ray.direction;
        double e = 0.25 * (hmax - hmin);
        for (int i = 0; i < 3; ++i)
        {
            center[i] += static_cast<double>(io.rawInteger(-1, 1));
            lo[i] = center[i] - e;
            hi[i] = center[i] + e;
        }
        AlignedBox3<double> box{};
        box.min = lo;
        box.max = hi;
        return box;
    }
}

namespace
{
    // Unrecorded draws of a box near the cone. A uniform box in [-4,4]^3
    // intersects a uniform cone on about 10 % of records, which leaves the
    // acceptance branches nearly untested; here the box centre is a point of
    // the cone axis perturbed by an offset comparable to the box size, so
    // both answers occur about equally often.
    AlignedBox3<double> RawBoxNearCone(oracle::Ctx& io, int mode,
        Cone3<double> const& cone, Vector3<double>& lo, Vector3<double>& hi)
    {
        double hmin = cone.GetMinHeight();
        double hmax = (cone.IsFinite() ? cone.GetMaxHeight() : hmin + 6.0);
        double h = hmin + (hmax - hmin) * io.raw(0.0, 1.0);
        Vector3<double> center = cone.ray.origin + h * cone.ray.direction;
        for (int i = 0; i < 3; ++i)
        {
            double e = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
                : io.raw(0.25, 3.0));
            double offset = (mode == 0
                ? static_cast<double>(io.rawInteger(-3, 3)) : io.raw(-3.0, 3.0));
            if (mode == 0) { center[i] = std::floor(center[i]); }
            lo[i] = center[i] + offset - e;
            hi[i] = center[i] + offset + e;
        }
        AlignedBox3<double> box{};
        box.min = lo;
        box.max = hi;
        return box;
    }
}

ORACLE_CASE("IntrAlignedBox3Cone3.test")
{
    int mode = io.index() % 2;
    auto cone = Cn<3>(io, mode, 3, 4.0, io.index() % 4);
    Vector3<double> lo{}, hi{};
    auto box = RawBoxNearCone(io, mode, cone, lo, hi);
    io.givenVec(lo);
    io.givenVec(hi);
    Box3Cone3Query query;
    auto r = query(box, cone);
    io.outBool(r.intersect);
}

// A uniform box rarely straddles a cone height plane, so the clipping path
// and its 81 face configurations are rarely reached. Here the box is centred
// on one of the two height planes on the lattice, with a coordinate-axis cone
// axis and integer heights, so the projections are exact integers and the
// '0' (Z) face configurations occur.
ORACLE_CASE("IntrAlignedBox3Cone3.test.straddle")
{
    auto cone = Cn<3>(io, 0, 3, 4.0, io.index() % 4);
    Vector3<double> lo{}, hi{};
    auto box = RawStraddlingBox(io, cone, lo, hi);
    io.givenVec(lo);
    io.givenVec(hi);
    Box3Cone3Query query;
    auto r = query(box, cone);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrAlignedBox3Cone3.computeBoxHeightInterval")
{
    int mode = io.index() % 2;
    auto box = ABox3(io, mode, 3, 4.0);
    auto cone = Cn<3>(io, mode, 3, 4.0, io.index() % 4);
    double boxMinHeight = 0.0, boxMaxHeight = 0.0;
    ExposeBox3Cone3::ComputeBoxHeightInterval(box, cone, boxMinHeight, boxMaxHeight);
    io.outReal(boxMinHeight);
    io.outReal(boxMaxHeight);
}

ORACLE_CASE("IntrAlignedBox3Cone3.coneAxisIntersectsBox")
{
    int mode = io.index() % 2;
    auto cone = Cn<3>(io, mode, 3, 4.0, io.index() % 4);
    Vector3<double> lo{}, hi{};
    auto box = RawBoxNearCone(io, mode, cone, lo, hi);
    io.givenVec(lo);
    io.givenVec(hi);
    io.outBool(ExposeBox3Cone3::ConeAxisIntersectsBox(box, cone));
}

ORACLE_CASE("IntrAlignedBox3Cone3.hasPointInsideCone")
{
    int mode = io.index() % 2;
    auto cone = Cn<3>(io, mode, 3, 4.0, io.index() % 4);
    // The two points are relative to the cone vertex, as the helper expects.
    auto P0 = Pt<3>(io, mode, 4, 5.0);
    auto P1 = Pt<3>(io, mode, 4, 5.0);
    io.outBool(ExposeBox3Cone3::HasPointInsideCone(P0, P1, cone));
}

// Deliberate deviation: upstream's BoxFullyInConeSlab copies the twelve box
// edges into mCandidateEdges and sets mNumCandidateEdges to 12 without
// clearing mAdjacencyMatrix, so the adjacency bits that an earlier clipping
// query set for the clipped vertices (indices >= 8) survive; the next
// ClearCandidates visits only the twelve box edges, and the next clipping
// query then silently drops those candidate edges and reports a false
// negative. The case runs three queries on ONE query object -- a clipping
// configuration, a fully-in-slab configuration and the first configuration
// again -- and emits all three answers. docs/UPSTREAM-FINDINGS.md
// IntrAlignedBox3Cone3.h BoxFullyInConeSlab, issue #301.
ORACLE_CASE("IntrAlignedBox3Cone3.test.staleAdjacencyDeviation")
{
    auto cone = Cn<3>(io, 0, 3, 4.0, 3);
    Vector3<double> loA{}, hiA{}, loB{}, hiB{};
    AlignedBox3<double> boxA{}, boxB{};
    for (;;)
    {
        boxA = RawStraddlingBox(io, cone, loA, hiA);
        boxB = RawInSlabBox(io, cone, loB, hiB);
        // The reused-object answer for the third query differs from the
        // fresh-object answer exactly when the stale bits matter.
        Box3Cone3Query reused;
        (void)reused(boxA, cone);
        (void)reused(boxB, cone);
        auto stale = reused(boxA, cone);
        Box3Cone3Query fresh;
        auto clean = fresh(boxA, cone);
        if (stale.intersect != clean.intersect) { break; }
    }
    io.givenVec(loA);
    io.givenVec(hiA);
    io.givenVec(loB);
    io.givenVec(hiB);
    Box3Cone3Query query;
    auto rA = query(boxA, cone);
    auto rB = query(boxB, cone);
    auto rC = query(boxA, cone);
    io.outBool(rA.intersect);
    io.outBool(rB.intersect);
    io.outBool(rC.intersect);
}

// ========================== IntrAreaEllipse2Ellipse2 =====================
//
// Upstream's AreaEllipse2Ellipse2 declares "T mZero, mOne, mTwo, mPi, mTwoPi;"
// with no constructor and never assigns them, so every area, every pi scaling
// and the 'dtheta <= mPi' branch read indeterminate values
// (docs/UPSTREAM-FINDINGS.md IntrAreaEllipse2Ellipse2.h, issue #301). The port
// fixes that with literals and the pi constants, and additionally normalizes
// private copies of the ellipse axes, because the polar angles and the sector
// integral are correct only for unit-length axes although the class comment
// claims otherwise (issue #465). Upstream's own code is therefore not a usable
// reference for the main case: the class below is upstream's code verbatim
// with exactly those two corrections applied, and it is what the main case
// compares the port against. The unmodified upstream class is exercised by the
// deviation case, which value-initializes it so that the five uninitialized
// members are deterministically zero.
//
// The compared computation calls atan2, atan, sin and cos, and the underlying
// FIQuery<Ellipse2,Ellipse2> calls RootsPolynomial::SolveQuartic (pow, cos,
// atan2), so the case is compared with a tolerance.

namespace
{
    class FixedAreaEllipse2Ellipse2
    {
    public:
        using EIQuery = FIQuery<double, Ellipse2<double>, Ellipse2<double>>;

        struct Result
        {
            enum class Configuration
            {
                ELLIPSES_ARE_EQUAL,
                ELLIPSES_ARE_SEPARATED,
                E0_CONTAINS_E1,
                E1_CONTAINS_E0,
                ONE_CHORD_REGION,
                FOUR_CHORD_REGION,
                INVALID
            };

            Result()
                :
                configuration(Configuration::INVALID),
                findResult{},
                area(0.0)
            {
            }

            Configuration configuration;
            typename EIQuery::Result findResult;
            double area;
        };

        Result operator()(Ellipse2<double> const& ellipse0, Ellipse2<double> const& ellipse1)
        {
            EllipseInfo E0{};
            E0.center = ellipse0.center;
            E0.axis = ellipse0.axis;
            E0.extent = ellipse0.extent;
            FinishEllipseInfo(E0);

            EllipseInfo E1{};
            E1.center = ellipse1.center;
            E1.axis = ellipse1.axis;
            E1.extent = ellipse1.extent;
            FinishEllipseInfo(E1);

            Result ar{};
            ar.configuration = Result::Configuration::INVALID;
            ar.findResult = EIQuery{}(ellipse0, ellipse1);
            ar.area = mZero;
            AreaDispatch(E0, E1, ar);
            return ar;
        }

    private:
        struct EllipseInfo
        {
            EllipseInfo()
                :
                center(Vector2<double>::Zero()),
                axis{ Vector2<double>::Zero() , Vector2<double>::Zero() },
                extent(Vector2<double>::Zero()),
                M{},
                AB(0.0),
                halfAB(0.0),
                BpA(0.0),
                BmA(0.0)
            {
            }

            Vector2<double> center;
            std::array<Vector2<double>, 2> axis;
            Vector2<double> extent;
            Matrix2x2<double> M;
            double AB, halfAB, BpA, BmA;
        };

        void FinishEllipseInfo(EllipseInfo& E)
        {
            E.M = OuterProduct(E.axis[0], E.axis[0]) /
                (E.extent[0] * E.extent[0] * Dot(E.axis[0], E.axis[0]));
            E.M += OuterProduct(E.axis[1], E.axis[1]) /
                (E.extent[1] * E.extent[1] * Dot(E.axis[1], E.axis[1]));
            // Port correction: normalize the private axis copies, which the
            // polar angles and the sector integral require. Componentwise
            // division, as the port does.
            double length0 = std::sqrt(Dot(E.axis[0], E.axis[0]));
            double length1 = std::sqrt(Dot(E.axis[1], E.axis[1]));
            for (int32_t d = 0; d < 2; ++d)
            {
                E.axis[0][d] /= length0;
                E.axis[1][d] /= length1;
            }
            E.AB = E.extent[0] * E.extent[1];
            E.halfAB = E.AB / mTwo;
            E.BpA = E.extent[1] + E.extent[0];
            E.BmA = E.extent[1] - E.extent[0];
        }

        void AreaDispatch(EllipseInfo const& E0, EllipseInfo const& E1, Result& ar)
        {
            if (ar.findResult.intersect)
            {
                if (ar.findResult.numPoints == 1)
                {
                    AreaCS(E0, E1, ar);
                }
                else if (ar.findResult.numPoints == 2)
                {
                    if (ar.findResult.isTransverse[0])
                    {
                        Area2(E0, E1, 0, 1, ar);
                    }
                    else
                    {
                        AreaCS(E0, E1, ar);
                    }
                }
                else if (ar.findResult.numPoints == 3)
                {
                    if (!ar.findResult.isTransverse[0])
                    {
                        Area2(E0, E1, 1, 2, ar);
                    }
                    else if (!ar.findResult.isTransverse[1])
                    {
                        Area2(E0, E1, 2, 0, ar);
                    }
                    else
                    {
                        Area2(E0, E1, 0, 1, ar);
                    }
                }
                else
                {
                    Area4(E0, E1, ar);
                }
            }
            else
            {
                AreaCS(E0, E1, ar);
            }
        }

        void AreaCS(EllipseInfo const& E0, EllipseInfo const& E1, Result& ar)
        {
            if (ar.findResult.numPoints <= 1)
            {
                Vector2<double> diff = E0.center - E1.center;
                double qform0 = Dot(diff, E0.M * diff);
                double qform1 = Dot(diff, E1.M * diff);
                if (qform0 > mOne && qform1 > mOne)
                {
                    ar.configuration = Result::Configuration::ELLIPSES_ARE_SEPARATED;
                    ar.area = mZero;
                }
                else
                {
                    if (E0.AB < E1.AB)
                    {
                        ar.configuration = Result::Configuration::E1_CONTAINS_E0;
                        ar.area = mPi * E0.AB;
                    }
                    else
                    {
                        ar.configuration = Result::Configuration::E0_CONTAINS_E1;
                        ar.area = mPi * E1.AB;
                    }
                }
            }
            else
            {
                ar.configuration = Result::Configuration::ELLIPSES_ARE_EQUAL;
                ar.area = mPi * E0.AB;
            }
        }

        void Area2(EllipseInfo const& E0, EllipseInfo const& E1, int32_t i0,
            int32_t i1, Result& ar)
        {
            ar.configuration = Result::Configuration::ONE_CHORD_REGION;
            Vector2<double> const& P0 = ar.findResult.points[i0];
            Vector2<double> const& P1 = ar.findResult.points[i1];
            Vector2<double> P0mC0 = P0 - E0.center, P0mC1 = P0 - E1.center;
            Vector2<double> P1mC0 = P1 - E0.center, P1mC1 = P1 - E1.center;
            Vector2<double> N0 = E0.M * P0mC0, N1 = E1.M * P0mC1;
            double dotperp = DotPerp(N1, N0);
            if (dotperp > mZero)
            {
                ar.area =
                    ComputeAreaChordRegion(E0, P0mC0, P1mC0) +
                    ComputeAreaChordRegion(E1, P1mC1, P0mC1);
            }
            else
            {
                ar.area =
                    ComputeAreaChordRegion(E0, P1mC0, P0mC0) +
                    ComputeAreaChordRegion(E1, P0mC1, P1mC1);
            }
        }

        void Area4(EllipseInfo const& E0, EllipseInfo const& E1, Result& ar)
        {
            ar.configuration = Result::Configuration::FOUR_CHORD_REGION;
            std::multimap<double, int32_t> ordering;
            int32_t i;
            for (i = 0; i < 4; ++i)
            {
                Vector2<double> PmC = ar.findResult.points[i] - E0.center;
                double x = Dot(E0.axis[0], PmC);
                double y = Dot(E0.axis[1], PmC);
                double theta = std::atan2(y, x);
                ordering.insert(std::make_pair(theta, i));
            }

            std::array<int32_t, 4> permute{ 0, 0, 0, 0 };
            i = 0;
            for (auto const& element : ordering)
            {
                permute[i++] = element.second;
            }

            Vector2<double> diag20 =
                ar.findResult.points[permute[2]] - ar.findResult.points[permute[0]];
            Vector2<double> diag31 =
                ar.findResult.points[permute[3]] - ar.findResult.points[permute[1]];
            ar.area = std::fabs(DotPerp(diag20, diag31)) / mTwo;

            for (int32_t i0 = 3, i1 = 0; i1 < 4; i0 = i1++)
            {
                Vector2<double> const& P0 = ar.findResult.points[permute[i0]];
                Vector2<double> const& P1 = ar.findResult.points[permute[i1]];
                Vector2<double> P0mC0 = P0 - E0.center, P0mC1 = P0 - E1.center;
                Vector2<double> P1mC0 = P1 - E0.center, P1mC1 = P1 - E1.center;
                Vector2<double> N0 = E0.M * P0mC0, N1 = E1.M * P0mC1;
                double dotperp = DotPerp(N1, N0);
                if (dotperp > mZero)
                {
                    ar.area += ComputeAreaChordRegion(E0, P0mC0, P1mC0);
                }
                else
                {
                    ar.area += ComputeAreaChordRegion(E1, P0mC1, P1mC1);
                }
            }
        }

        double ComputeAreaChordRegion(EllipseInfo const& E,
            Vector2<double> const& P0mC, Vector2<double> const& P1mC)
        {
            double x0 = Dot(E.axis[0], P0mC);
            double y0 = Dot(E.axis[1], P0mC);
            double theta0 = std::atan2(y0, x0);
            double x1 = Dot(E.axis[0], P1mC);
            double y1 = Dot(E.axis[1], P1mC);
            double theta1 = std::atan2(y1, x1);
            if (theta1 < theta0)
            {
                theta1 += mTwoPi;
            }
            double triArea = std::fabs(DotPerp(P0mC, P1mC)) / mTwo;
            double dtheta = theta1 - theta0;
            double F0, F1, sectorArea;
            if (dtheta <= mPi)
            {
                F0 = ComputeIntegral(E, theta0);
                F1 = ComputeIntegral(E, theta1);
                sectorArea = F1 - F0;
                return sectorArea - triArea;
            }
            else
            {
                theta0 += mTwoPi;
                F0 = ComputeIntegral(E, theta0);
                F1 = ComputeIntegral(E, theta1);
                sectorArea = F0 - F1;
                return mPi * E.AB - (sectorArea - triArea);
            }
        }

        double ComputeIntegral(EllipseInfo const& E, double const& theta)
        {
            double twoTheta = mTwo * theta;
            double sn = std::sin(twoTheta);
            double cs = std::cos(twoTheta);
            double arg = E.BmA * sn / (E.BpA + E.BmA * cs);
            return E.halfAB * (theta - std::atan(arg));
        }

        // Port correction: the five constants upstream never assigns.
        double mZero = 0.0;
        double mOne = 1.0;
        double mTwo = 2.0;
        double mPi = GTE_C_PI;
        double mTwoPi = GTE_C_TWO_PI;
    };

    // Two ellipses that overlap often enough for the chord-region branches to
    // be reached: the second centre is a short offset from the first.
    void RawAreaEllipses(oracle::Ctx& io, int mode, Ellipse2<double>& e0,
        Ellipse2<double>& e1)
    {
        for (int k = 0; k < 2; ++k)
        {
            Ellipse2<double>& e = (k == 0 ? e0 : e1);
            for (int d = 0; d < 2; ++d)
            {
                e.center[d] = (mode == 0 ? static_cast<double>(io.rawInteger(-2, 2))
                    : io.raw(-2.0, 2.0));
            }
            if (mode == 0)
            {
                int j = io.rawInteger(0, 3);
                e.axis[0].MakeZero();
                e.axis[0][j % 2] = (j < 2 ? 1.0 : -1.0);
            }
            else
            {
                double len = 0.0;
                do
                {
                    for (int d = 0; d < 2; ++d) { e.axis[0][d] = io.raw(-1.0, 1.0); }
                    len = Length(e.axis[0]);
                } while (len < 0.1 || len > 1.0);
                Normalize(e.axis[0]);
            }
            e.axis[1] = Perp(e.axis[0]);
            for (int d = 0; d < 2; ++d)
            {
                e.extent[d] = (mode == 0 ? static_cast<double>(io.rawInteger(1, 3))
                    : io.raw(0.5, 3.0));
            }
        }
    }

    void RecordEllipse(oracle::Ctx& io, Ellipse2<double> const& e)
    {
        io.givenVec(e.center);
        io.givenVec(e.axis[0]);
        io.givenVec(e.axis[1]);
        io.givenVec(e.extent);
    }

    // FIQuery<Ellipse2,Ellipse2>::Result::numPoints uses SIZE_MAX for
    // "the ellipses are the same", which is not representable as a double;
    // the port uses Number.MAX_SAFE_INTEGER. Both sides emit -1 for it, as
    // group 32 does.
    void OutNumPoints(oracle::Ctx& io, size_t n)
    {
        io.outInt(n == std::numeric_limits<size_t>::max()
            ? -1 : static_cast<int32_t>(n));
    }
}

// The main case compares the port against upstream's algorithm with the two
// documented corrections applied (the class above).
ORACLE_CASE("IntrAreaEllipse2Ellipse2.compute")
{
    int mode = io.index() % 2;
    Ellipse2<double> e0{}, e1{};
    RawAreaEllipses(io, mode, e0, e1);
    RecordEllipse(io, e0);
    RecordEllipse(io, e1);
    FixedAreaEllipse2Ellipse2 query;
    auto r = query(e0, e1);
    io.outInt(static_cast<int32_t>(r.configuration));
    OutNumPoints(io, r.findResult.numPoints);
    io.outReal(r.area);
}

// Deliberate deviation: the unmodified upstream class. Value-initializing it
// makes the five members that upstream never assigns deterministically zero
// (the class has no user-provided constructor, so 'query{}' zero-initializes
// it), which is the only way to get a reproducible golden out of code that
// otherwise reads indeterminate stack values. With mTwo = 0 the half-area is
// infinite, with mPi = 0 every containment area is zero and with mTwoPi = 0
// the atan2 wrap never happens, so the reported area is garbage on every
// record. docs/UPSTREAM-FINDINGS.md IntrAreaEllipse2Ellipse2.h, issue #301.
ORACLE_CASE("IntrAreaEllipse2Ellipse2.compute.uninitializedDeviation")
{
    int mode = io.index() % 2;
    Ellipse2<double> e0{}, e1{};
    RawAreaEllipses(io, mode, e0, e1);
    RecordEllipse(io, e0);
    RecordEllipse(io, e1);
    AreaEllipse2Ellipse2<double> query{};
    auto r = query(e0, e1);
    io.outInt(static_cast<int32_t>(r.configuration));
    OutNumPoints(io, r.findResult.numPoints);
    io.outReal(r.area);
}

// ========================== IntrConvexMesh3Plane3 ========================
//
// Upstream static_asserts that Real is an arbitrary-precision type with a
// division operator, because only then are the sign tests Dot(N,X) - c exact.
// The port instantiates the query for 'number' like every other query in the
// library, and its header comment says the floating-point behaviour is
// exactly what a C++ instantiation with double would do if the static_assert
// were removed. APD below is a transparent wrapper around double -- every
// operation is the corresponding IEEE double operation -- that is *declared*
// to be an arbitrary-precision division type, which instantiates upstream's
// algorithm on plain doubles without touching the trait for 'double' itself
// and without editing the upstream header. The comparison is therefore
// bit-for-bit against upstream's own code.

namespace v34ap
{
    struct APD
    {
        APD() : v(0.0) {}
        APD(double x) : v(x) {}
        explicit operator double() const { return v; }
        double v;
    };

    inline APD operator+(APD a) { return a; }
    inline APD operator-(APD a) { return APD(-a.v); }
    inline APD operator+(APD a, APD b) { return APD(a.v + b.v); }
    inline APD operator-(APD a, APD b) { return APD(a.v - b.v); }
    inline APD operator*(APD a, APD b) { return APD(a.v * b.v); }
    inline APD operator/(APD a, APD b) { return APD(a.v / b.v); }
    inline APD& operator+=(APD& a, APD b) { a.v += b.v; return a; }
    inline APD& operator-=(APD& a, APD b) { a.v -= b.v; return a; }
    inline APD& operator*=(APD& a, APD b) { a.v *= b.v; return a; }
    inline APD& operator/=(APD& a, APD b) { a.v /= b.v; return a; }
    inline bool operator==(APD a, APD b) { return a.v == b.v; }
    inline bool operator!=(APD a, APD b) { return a.v != b.v; }
    inline bool operator< (APD a, APD b) { return a.v < b.v; }
    inline bool operator<=(APD a, APD b) { return a.v <= b.v; }
    inline bool operator> (APD a, APD b) { return a.v > b.v; }
    inline bool operator>=(APD a, APD b) { return a.v >= b.v; }
}

namespace gte
{
    template <>
    struct _is_arbitrary_precision_internal<v34ap::APD> : std::true_type {};

    template <>
    struct _has_division_operator_internal<v34ap::APD> : std::true_type {};
}

namespace
{
    using APD = v34ap::APD;
    using APMesh = ConvexMesh3<APD>;
    using APQuery = FIQuery<APD, ConvexMesh3<APD>, Plane3<APD>>;

    Vector3<APD> ToAP(Vector3<double> const& v)
    {
        return Vector3<APD>{ v[0], v[1], v[2] };
    }

    void OutAPVec(oracle::Ctx& io, Vector3<APD> const& v)
    {
        for (int32_t d = 0; d < 3; ++d) { io.outReal(static_cast<double>(v[d])); }
    }

    void OutAPMesh(oracle::Ctx& io, APMesh const& mesh)
    {
        io.outInt(mesh.configuration);
        io.outInt(static_cast<int32_t>(mesh.vertices.size()));
        for (auto const& vertex : mesh.vertices) { OutAPVec(io, vertex); }
        io.outInt(static_cast<int32_t>(mesh.triangles.size()));
        for (auto const& triangle : mesh.triangles)
        {
            io.outInt(triangle[0]);
            io.outInt(triangle[1]);
            io.outInt(triangle[2]);
        }
    }

    void OutAPResult(oracle::Ctx& io, APQuery::Result const& r, bool withPolyhedra)
    {
        io.outInt(r.configuration);
        io.outInt(r.requested);
        OutAPMesh(io, r.intersectionMesh);
        io.outInt(static_cast<int32_t>(r.intersectionPolygon.size()));
        for (auto const& vertex : r.intersectionPolygon) { OutAPVec(io, vertex); }
        if (withPolyhedra)
        {
            OutAPMesh(io, r.positivePolyhedron);
            OutAPMesh(io, r.negativePolyhedron);
        }
    }

    // Unrecorded unit-length (or exact coordinate-axis) direction.
    void RawDir3(oracle::Ctx& io, int mode, Vector3<double>& v)
    {
        if (mode == 0)
        {
            int k = io.rawInteger(0, 5);
            v.MakeZero();
            v[k % 3] = (k < 3 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int d = 0; d < 3; ++d) { v[d] = io.raw(-1.0, 1.0); }
                len = Length(v);
            } while (len < 0.1 || len > 1.0);
            Normalize(v);
        }
    }

    // A tetrahedron with vertices T, T+(a,0,0), T+(0,b,0), T+(0,0,c) and the
    // four faces in outward-counterclockwise order. Four vertices keep the
    // split output small. The draws are unrecorded; the caller records the
    // vertices.
    APMesh Tetra(oracle::Ctx& io, int mode, std::array<Vector3<double>, 4>& v)
    {
        Vector3<double> T{}, e{};
        for (int32_t d = 0; d < 3; ++d)
        {
            T[d] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-3.0, 3.0));
            e[d] = (mode == 0 ? static_cast<double>(io.rawInteger(1, 4))
                : io.raw(0.5, 4.0));
        }
        v[0] = T;
        v[1] = T; v[1][0] += e[0];
        v[2] = T; v[2][1] += e[1];
        v[3] = T; v[3][2] += e[2];
        APMesh mesh{};
        mesh.configuration = APMesh::CFG_POLYHEDRON;
        mesh.vertices = { ToAP(v[0]), ToAP(v[1]), ToAP(v[2]), ToAP(v[3]) };
        mesh.triangles = { { 0, 2, 1 }, { 0, 1, 3 }, { 0, 3, 2 }, { 1, 2, 3 } };
        return mesh;
    }

    // An axis-aligned box as a 12-triangle convex mesh with outward
    // counterclockwise faces. Each box face is two triangles, so a plane
    // containing a face gives a four-vertex coplanar polygon.
    APMesh Box(oracle::Ctx& io, int mode, std::array<Vector3<double>, 8>& v)
    {
        Vector3<double> lo{}, hi{};
        for (int32_t d = 0; d < 3; ++d)
        {
            double a = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
                : io.raw(-3.0, 3.0));
            double e = (mode == 0 ? static_cast<double>(io.rawInteger(1, 4))
                : io.raw(0.5, 4.0));
            lo[d] = a;
            hi[d] = a + e;
        }
        double const x[2] = { lo[0], hi[0] };
        double const y[2] = { lo[1], hi[1] };
        double const z[2] = { lo[2], hi[2] };
        int32_t k = 0;
        for (int32_t iz = 0; iz < 2; ++iz)
        {
            for (int32_t iy = 0; iy < 2; ++iy)
            {
                for (int32_t ix = 0; ix < 2; ++ix)
                {
                    v[static_cast<size_t>(k++)] = Vector3<double>{ x[ix], y[iy], z[iz] };
                }
            }
        }
        APMesh mesh{};
        mesh.configuration = APMesh::CFG_POLYHEDRON;
        mesh.vertices.resize(8);
        for (size_t i = 0; i < 8; ++i) { mesh.vertices[i] = ToAP(v[i]); }
        mesh.triangles =
        {
            { 0, 2, 3 }, { 0, 3, 1 },   // -z
            { 4, 5, 7 }, { 4, 7, 6 },   // +z
            { 0, 1, 5 }, { 0, 5, 4 },   // -y
            { 2, 6, 7 }, { 2, 7, 3 },   // +y
            { 0, 4, 6 }, { 0, 6, 2 },   // -x
            { 1, 3, 7 }, { 1, 7, 5 }    // +x
        };
        return mesh;
    }

    // True when the plane contains three or more mesh vertices, that is when
    // the query reaches GetIntersectionPolygon for a coplanar face, which is
    // exactly the configuration on which the port deviates (see the deviation
    // case). The bit 4 of the configuration is upstream's own "numZero >= 3"
    // flag.
    bool CoplanarFace(APMesh const& mesh, Plane3<APD> const& plane)
    {
        APQuery query;
        auto r = query(mesh, plane, APQuery::REQ_CONFIGURATION_ONLY);
        return (r.configuration & 4) != 0;
    }
}

// The transverse split of a tetrahedron, with every feature requested: the
// polygon of intersection, its triangulation and both split polyhedra. The
// generator rejects the planes that contain a face of the mesh.
ORACLE_CASE("IntrConvexMesh3Plane3.find.tetrahedron")
{
    int mode = io.index() % 2;
    std::array<Vector3<double>, 4> v{};
    APMesh mesh{};
    Vector3<double> normal{}, planeOrigin{};
    for (;;)
    {
        mesh = Tetra(io, mode, v);
        RawDir3(io, mode, normal);
        for (int32_t d = 0; d < 3; ++d)
        {
            planeOrigin[d] = 0.25 * (v[0][d] + v[1][d] + v[2][d] + v[3][d]);
            planeOrigin[d] += (mode == 0
                ? 0.5 * static_cast<double>(io.rawInteger(-2, 2)) : io.raw(-1.0, 1.0));
        }
        Plane3<APD> probe(ToAP(normal), ToAP(planeOrigin));
        if (!CoplanarFace(mesh, probe)) { break; }
    }
    for (size_t i = 0; i < 4; ++i) { io.givenVec(v[i]); }
    io.givenVec(normal);
    io.givenVec(planeOrigin);
    Plane3<APD> plane(ToAP(normal), ToAP(planeOrigin));
    APQuery query;
    auto r = query(mesh, plane, APQuery::REQ_ALL);
    OutAPResult(io, r, true);
}

// Only the configuration is requested, so none of the intersection code runs.
// This case does not reject the coplanar-face planes, which is how the
// CFG_*_SIDE_POLYGON configuration values are covered.
ORACLE_CASE("IntrConvexMesh3Plane3.find.configurationOnly")
{
    int mode = io.index() % 2;
    std::array<Vector3<double>, 4> v{};
    auto mesh = Tetra(io, mode, v);
    for (size_t i = 0; i < 4; ++i) { io.givenVec(v[i]); }
    Vector3<double> normal{};
    RawDir3(io, mode, normal);
    io.givenVec(normal);
    Vector3<double> planeOrigin{};
    for (int32_t d = 0; d < 3; ++d)
    {
        planeOrigin[d] = (mode == 0 ? static_cast<double>(io.rawInteger(-3, 3))
            : io.raw(-4.0, 4.0));
    }
    io.givenVec(planeOrigin);
    Plane3<APD> plane(ToAP(normal), ToAP(planeOrigin));
    APQuery query;
    auto r = query(mesh, plane, APQuery::REQ_CONFIGURATION_ONLY);
    OutAPResult(io, r, true);
}

// The tangential configurations that uniform planes never reach: the plane
// contains exactly one vertex (CFG_*_SIDE_VERTEX) or exactly one edge
// (CFG_*_SIDE_EDGE) of the lattice tetrahedron, with all-integer plane
// coefficients so that the sign tests are exact. The face-containing plane is
// deliberately absent: upstream is defective for it and it is covered by the
// deviation case below.
ORACLE_CASE("IntrConvexMesh3Plane3.find.tangent")
{
    std::array<Vector3<double>, 4> v{};
    auto mesh = Tetra(io, 0, v);
    for (size_t i = 0; i < 4; ++i) { io.givenVec(v[i]); }
    int32_t which = io.rawInteger(0, 2);
    double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    Vector3<double> normal{};
    Vector3<double> through = v[0];
    if (which == 0)
    {
        // The apex v3 alone.
        normal = { 0.0, 0.0, sign };
        through = v[3];
    }
    else if (which == 1)
    {
        // The edge <v0,v1>.
        normal = { 0.0, sign, sign };
    }
    else
    {
        // The vertex v2 alone.
        normal = { -sign, sign, -sign };
        through = v[2];
    }
    io.givenVec(normal);
    io.givenVec(through);
    Plane3<APD> plane(ToAP(normal), ToAP(through));
    APQuery query;
    auto r = query(mesh, plane, APQuery::REQ_ALL);
    OutAPResult(io, r, true);
}

// A box split transversely, which gives polygons of intersection with four,
// five or six vertices. Only the intersection mesh and polygon are requested,
// to keep the record small. The generator rejects the planes that contain a
// face of the box.
ORACLE_CASE("IntrConvexMesh3Plane3.find.box")
{
    int mode = io.index() % 2;
    std::array<Vector3<double>, 8> v{};
    APMesh mesh{};
    Vector3<double> normal{}, planeOrigin{};
    for (;;)
    {
        mesh = Box(io, mode, v);
        RawDir3(io, mode, normal);
        for (int32_t d = 0; d < 3; ++d)
        {
            planeOrigin[d] = 0.5 * (v[0][d] + v[7][d]);
            planeOrigin[d] += (mode == 0
                ? 0.5 * static_cast<double>(io.rawInteger(-2, 2)) : io.raw(-1.0, 1.0));
        }
        Plane3<APD> probe(ToAP(normal), ToAP(planeOrigin));
        if (!CoplanarFace(mesh, probe)) { break; }
    }
    for (size_t i = 0; i < 8; ++i) { io.givenVec(v[i]); }
    io.givenVec(normal);
    io.givenVec(planeOrigin);
    Plane3<APD> plane(ToAP(normal), ToAP(planeOrigin));
    APQuery query;
    auto r = query(mesh, plane, APQuery::REQ_INTR_BOTH);
    OutAPResult(io, r, false);
}

// Deliberate deviation: the plane contains a whole box face, so
// GetIntersectionPolygon builds the predecessor map 'polygonIndices' and then
// reads it back in index order instead of traversing it as a cycle, which
// reports an arbitrary permutation of the boundary instead of the polygon.
// For the four-vertex box face the boundary cycle is 0 -> 2 -> 3 -> 1 after
// the vertices are repacked, and upstream's readback is (1,3,0,2), which is
// not a rotation of it. docs/UPSTREAM-FINDINGS.md IntrConvexMesh3Plane3.h
// GetIntersectionPolygon, issue #301.
ORACLE_CASE("IntrConvexMesh3Plane3.find.coplanarFaceDeviation")
{
    int mode = io.index() % 2;
    std::array<Vector3<double>, 8> v{};
    auto mesh = Box(io, mode, v);
    for (size_t i = 0; i < 8; ++i) { io.givenVec(v[i]); }
    int32_t k = io.rawInteger(0, 2);
    double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    Vector3<double> normal{};
    normal.MakeZero();
    normal[k] = sign;
    io.givenVec(normal);
    // The face plane: the low face for sign > 0, the high face otherwise.
    Vector3<double> through = (sign > 0.0 ? v[0] : v[7]);
    io.givenVec(through);
    Plane3<APD> plane(ToAP(normal), ToAP(through));
    APQuery query;
    auto r = query(mesh, plane, APQuery::REQ_INTR_BOTH);
    OutAPResult(io, r, false);
}
