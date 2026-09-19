// Verify group 22 (distance): differential cases for the Dist* headers that
// group 22 covers. See ORACLE.md.
//
// Generator modes are selected from io.index() so that every case mixes
// uniform random inputs, small-lattice inputs (exact integer arithmetic, which
// produces parallel, touching, coplanar and coincident configurations) and a
// construction aimed at a specific branch. Every generator records exactly the
// same number of doubles in every mode, so the TypeScript replay reads the
// inputs without knowing which mode produced them.
#define ORACLE_FAMILY "v22-distance"
#include "Oracle.h"

#include <Mathematics/DistAlignedBox3OrientedBox3.h>
#include <Mathematics/DistCircle3Circle3.h>
#include <Mathematics/DistOrientedBox3Cone3.h>
#include <Mathematics/DistOrientedBox3OrientedBox3.h>
#include <Mathematics/DistPlane3AlignedBox3.h>
#include <Mathematics/DistPlane3OrientedBox3.h>
#include <Mathematics/DistPoint3Parallelepiped3.h>
#include <Mathematics/DistPoint3Tetrahedron3.h>
#include <Mathematics/DistRectangle3AlignedBox3.h>
#include <Mathematics/DistRectangle3CanonicalBox3.h>
#include <Mathematics/DistRectangle3OrientedBox3.h>
#include <Mathematics/DistRectangle3Rectangle3.h>
#include <Mathematics/DistTetrahedron3Tetrahedron3.h>
#include <Mathematics/DistTriangle3AlignedBox3.h>
#include <Mathematics/DistTriangle3CanonicalBox3.h>
#include <Mathematics/DistTriangle3OrientedBox3.h>
#include <Mathematics/DistTriangle3Rectangle3.h>
#include <Mathematics/DistTriangle3Triangle3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>

using namespace gte;

namespace
{
    // ---- generators ------------------------------------------------------
    //
    // mode 0,1: uniform random
    // mode 2,3: small lattice (integers), so degenerate configurations occur
    // mode 4:   constructed (contact, coplanar, coincident, ...)
    int32_t Mode(oracle::Ctx& io)
    {
        return io.index() % 5;
    }

    bool Lattice(int32_t mode)
    {
        return mode == 2 || mode == 3;
    }

    template <int32_t N>
    Vector<N, double> Point(oracle::Ctx& io, int32_t mode)
    {
        if (Lattice(mode))
        {
            return io.latticeVec<N>(-3, 3);
        }
        return io.vec<N>(-5.0, 5.0);
    }

    // Extents are nonnegative; a zero extent (a degenerate box or rectangle)
    // is allowed by the headers and is generated in the lattice modes.
    template <int32_t N>
    Vector<N, double> Extent(oracle::Ctx& io, int32_t mode)
    {
        if (Lattice(mode))
        {
            return io.latticeVec<N>(0, 2);
        }
        return io.vec<N>(0.25, 3.0);
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

    // An orthonormal triple of 3D axes; the coordinate axes (in a random
    // order and with random signs, so the frame stays right-handed) in the
    // lattice modes. A box whose axes are the coordinate axes keeps the
    // exactly zero components of lattice inputs in the box frame, which is
    // where the canonical-box queries take their degenerate branches.
    std::array<Vector3<double>, 3> Axes3(oracle::Ctx& io, int32_t mode)
    {
        Vector3<double> a0{}, a1{}, a2{};
        if (Lattice(mode))
        {
            int32_t k = io.rawInteger(0, 2);
            int32_t s = io.rawInteger(0, 1);
            a0.MakeUnit(k);
            a1.MakeUnit((k + 1) % 3);
            if (s == 1) { a0 = -a0; a1 = -a1; }
            a2 = Cross(a0, a1);
        }
        else
        {
            a0 = RawUnit<3>(io);
            double len = 0.0;
            do
            {
                a1 = RawVec<3>(io, -1.0, 1.0);
                a1 = a1 - Dot(a1, a0) * a0;
                len = Length(a1);
            } while (len < 0.25);
            Normalize(a1);
            a2 = Cross(a0, a1);
        }
        std::array<Vector3<double>, 3> axis{};
        axis[0] = io.givenVec(a0);
        axis[1] = io.givenVec(a1);
        axis[2] = io.givenVec(a2);
        return axis;
    }

    // An orthonormal pair of 3D directions, for a rectangle. Records 6.
    std::array<Vector3<double>, 2> RectAxes(oracle::Ctx& io, int32_t mode)
    {
        Vector3<double> a0{}, a1{};
        if (Lattice(mode))
        {
            int32_t k = io.rawInteger(0, 2);
            int32_t s = io.rawInteger(0, 1);
            a0.MakeUnit(k);
            a1.MakeUnit((k + 1) % 3);
            if (s == 1) { a1 = -a1; }
        }
        else
        {
            a0 = RawUnit<3>(io);
            double len = 0.0;
            do
            {
                a1 = RawVec<3>(io, -1.0, 1.0);
                a1 = a1 - Dot(a1, a0) * a0;
                len = Length(a1);
            } while (len < 0.25);
            Normalize(a1);
        }
        std::array<Vector3<double>, 2> axis{};
        axis[0] = io.givenVec(a0);
        axis[1] = io.givenVec(a1);
        return axis;
    }

    // A unit-length plane normal. The lattice modes use a signed coordinate
    // axis so that Dot(N, X) is exact on lattice points and the plane-box
    // query sees exact zeros.
    Vector3<double> UnitNormal(oracle::Ctx& io, int32_t mode)
    {
        if (Lattice(mode))
        {
            int32_t k = io.rawInteger(0, 2);
            int32_t s = io.rawInteger(0, 1);
            Vector3<double> n{ 0.0, 0.0, 0.0 };
            n[k] = (s == 0 ? -1.0 : 1.0);
            return io.givenVec(n);
        }
        return io.unit<3>();
    }

    // A box in min/size form; the maximum corner is derived so that it is
    // never below the minimum corner. Records 2N doubles.
    template <int32_t N>
    void AlignedBoxInputs(oracle::Ctx& io, int32_t mode,
        Vector<N, double>& lo, Vector<N, double>& hi)
    {
        lo = Point<N>(io, mode);
        Vector<N, double> size{};
        if (Lattice(mode))
        {
            size = io.latticeVec<N>(0, 3);
        }
        else
        {
            size = io.vec<N>(0.0, 4.0);
        }
        hi = lo + size;
    }

    // ---- output helpers --------------------------------------------------
    template <typename Result>
    void EmitTwoClosest(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    template <typename Result>
    void EmitCartesian(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.cartesian[0]);
        io.outReal(r.cartesian[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    template <typename Result>
    void EmitBarycentric3(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.barycentric[0]);
        io.outReal(r.barycentric[1]);
        io.outReal(r.barycentric[2]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    // A unit-length normal with exactly one zero component, which is the
    // DoQuery2D dispatch arm of the plane-canonical-box query when the box
    // axes are the coordinate axes.
    Vector3<double> PlanarNormal(oracle::Ctx& io)
    {
        Vector2<double> u = RawUnit<2>(io);
        int32_t k = io.rawInteger(0, 2);
        Vector3<double> n{ 0.0, 0.0, 0.0 };
        n[(k + 1) % 3] = u[0];
        n[(k + 2) % 3] = u[1];
        return io.givenVec(n);
    }
}

// =====================================================================
// DistPlane3AlignedBox3.h, DistPlane3OrientedBox3.h
// =====================================================================
ORACLE_CASE("DistPlane3AlignedBox3.compute")
{
    // The lattice modes give a coordinate-axis normal (DoQuery1D) and mode 4
    // a normal with one zero component (DoQuery2D) whose plane passes exactly
    // through the center of the box, so the distance is exactly zero.
    int32_t mode = Mode(io);
    Vector3<double> lo{}, hi{};
    AlignedBoxInputs<3>(io, mode, lo, hi);
    Vector3<double> normal{}, origin{};
    if (mode == 4)
    {
        normal = PlanarNormal(io);
        origin = io.givenVec(0.5 * (lo + hi));
    }
    else
    {
        normal = UnitNormal(io, mode);
        origin = Point<3>(io, mode);
    }
    Plane3<double> plane(normal, origin);
    AlignedBox3<double> box(lo, hi);
    DCPQuery<double, Plane3<double>, AlignedBox3<double>> query{};
    auto r = query(plane, box);
    EmitTwoClosest(io, r);
}

ORACLE_CASE("DistPlane3OrientedBox3.compute")
{
    // Mode 4 sends the plane exactly through the center of the box.
    int32_t mode = Mode(io);
    auto center = Point<3>(io, mode);
    auto axis = Axes3(io, mode);
    auto extent = Extent<3>(io, mode);
    Vector3<double> normal{}, origin{};
    if (mode == 4)
    {
        normal = PlanarNormal(io);
        origin = io.givenVec(center);
    }
    else
    {
        normal = UnitNormal(io, mode);
        origin = Point<3>(io, mode);
    }
    Plane3<double> plane(normal, origin);
    OrientedBox3<double> box(center, axis, extent);
    DCPQuery<double, Plane3<double>, OrientedBox3<double>> query{};
    auto r = query(plane, box);
    EmitTwoClosest(io, r);
}

// =====================================================================
// DistRectangle3CanonicalBox3.h, DistRectangle3AlignedBox3.h,
// DistRectangle3OrientedBox3.h
// =====================================================================
namespace
{
    // center (3), the two axes (6) and the extent (2).
    Rectangle3<double> RectangleInputs(oracle::Ctx& io, int32_t mode,
        Vector3<double> const& contact, bool useContact)
    {
        Vector3<double> center{};
        if (useContact)
        {
            center = io.givenVec(contact);
        }
        else
        {
            center = Point<3>(io, mode);
        }
        std::array<Vector3<double>, 2> axis = RectAxes(io, mode);
        Vector2<double> extent = Extent<2>(io, mode);
        return Rectangle3<double>(center, axis, extent);
    }

    // A point strictly inside the canonical box of the given extent.
    Vector3<double> RawInsideBox(oracle::Ctx& io, Vector3<double> const& extent)
    {
        Vector3<double> p{};
        for (int32_t i = 0; i < 3; ++i)
        {
            p[i] = io.raw(-1.0, 1.0) * extent[i];
        }
        return p;
    }
}

ORACLE_CASE("DistRectangle3CanonicalBox3.compute")
{
    // Mode 4 centers the rectangle inside the box, so the plane-box closest
    // point lands inside the rectangle and the early (plane) branch is taken
    // with distance exactly zero; the lattice modes make the rectangle
    // coplanar with a face of the box and reach the four-edge branch.
    int32_t mode = Mode(io);
    auto extent = Extent<3>(io, mode);
    Vector3<double> contact{};
    if (mode == 4) { contact = RawInsideBox(io, extent); }
    Rectangle3<double> rectangle = RectangleInputs(io, mode, contact, mode == 4);
    CanonicalBox3<double> box(extent);
    DCPQuery<double, Rectangle3<double>, CanonicalBox3<double>> query{};
    auto r = query(rectangle, box);
    EmitCartesian(io, r);
}

ORACLE_CASE("DistRectangle3AlignedBox3.compute")
{
    // Mode 4 centers the rectangle at the center of the box.
    int32_t mode = Mode(io);
    Vector3<double> lo{}, hi{};
    AlignedBoxInputs<3>(io, mode, lo, hi);
    Vector3<double> contact = 0.5 * (lo + hi);
    Rectangle3<double> rectangle = RectangleInputs(io, mode, contact, mode == 4);
    AlignedBox3<double> box(lo, hi);
    DCPQuery<double, Rectangle3<double>, AlignedBox3<double>> query{};
    auto r = query(rectangle, box);
    EmitCartesian(io, r);
}

ORACLE_CASE("DistRectangle3OrientedBox3.compute")
{
    // Mode 4 centers the rectangle at the center of the box.
    int32_t mode = Mode(io);
    auto center = Point<3>(io, mode);
    auto axis = Axes3(io, mode);
    auto extent = Extent<3>(io, mode);
    Rectangle3<double> rectangle = RectangleInputs(io, mode, center, mode == 4);
    OrientedBox3<double> box(center, axis, extent);
    DCPQuery<double, Rectangle3<double>, OrientedBox3<double>> query{};
    auto r = query(rectangle, box);
    EmitCartesian(io, r);
}

// =====================================================================
// DistTriangle3CanonicalBox3.h, DistTriangle3AlignedBox3.h,
// DistTriangle3OrientedBox3.h
// =====================================================================
namespace
{
    // Three vertices (9). When 'useContact' is true the first vertex is the
    // given point, which puts the triangle in contact with the solid.
    Triangle3<double> TriangleInputs(oracle::Ctx& io, int32_t mode,
        Vector3<double> const& contact, bool useContact)
    {
        Vector3<double> v0{};
        if (useContact)
        {
            v0 = io.givenVec(contact);
        }
        else
        {
            v0 = Point<3>(io, mode);
        }
        Vector3<double> v1 = Point<3>(io, mode);
        Vector3<double> v2 = Point<3>(io, mode);
        return Triangle3<double>(v0, v1, v2);
    }
}

ORACLE_CASE("DistTriangle3CanonicalBox3.compute")
{
    // Mode 4 puts the first vertex strictly inside the box (distance exactly
    // zero); the lattice modes produce degenerate (collinear) triangles, for
    // which upstream divides by a zero squared normal length and the
    // barycentric coordinates are NaN. That quirk is preserved by the port.
    int32_t mode = Mode(io);
    auto extent = Extent<3>(io, mode);
    Vector3<double> contact{};
    if (mode == 4) { contact = RawInsideBox(io, extent); }
    Triangle3<double> triangle = TriangleInputs(io, mode, contact, mode == 4);
    CanonicalBox3<double> box(extent);
    DCPQuery<double, Triangle3<double>, CanonicalBox3<double>> query{};
    auto r = query(triangle, box);
    EmitBarycentric3(io, r);
}

ORACLE_CASE("DistTriangle3AlignedBox3.compute")
{
    int32_t mode = Mode(io);
    Vector3<double> lo{}, hi{};
    AlignedBoxInputs<3>(io, mode, lo, hi);
    Vector3<double> contact = 0.5 * (lo + hi);
    Triangle3<double> triangle = TriangleInputs(io, mode, contact, mode == 4);
    AlignedBox3<double> box(lo, hi);
    DCPQuery<double, Triangle3<double>, AlignedBox3<double>> query{};
    auto r = query(triangle, box);
    EmitBarycentric3(io, r);
}

ORACLE_CASE("DistTriangle3OrientedBox3.compute")
{
    int32_t mode = Mode(io);
    auto center = Point<3>(io, mode);
    auto axis = Axes3(io, mode);
    auto extent = Extent<3>(io, mode);
    Triangle3<double> triangle = TriangleInputs(io, mode, center, mode == 4);
    OrientedBox3<double> box(center, axis, extent);
    DCPQuery<double, Triangle3<double>, OrientedBox3<double>> query{};
    auto r = query(triangle, box);
    EmitBarycentric3(io, r);
}

// =====================================================================
// DistRectangle3Rectangle3.h, DistTriangle3Rectangle3.h,
// DistTriangle3Triangle3.h
// =====================================================================
ORACLE_CASE("DistRectangle3Rectangle3.compute")
{
    // Mode 4 gives the two rectangles a common center, so they intersect and
    // the distance is exactly zero; the lattice modes make them coplanar or
    // parallel.
    int32_t mode = Mode(io);
    Vector3<double> unused{};
    Rectangle3<double> rectangle0 = RectangleInputs(io, mode, unused, false);
    Rectangle3<double> rectangle1 =
        RectangleInputs(io, mode, rectangle0.center, mode == 4);
    DCPQuery<double, Rectangle3<double>, Rectangle3<double>> query{};
    auto r = query(rectangle0, rectangle1);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.cartesian0[0]);
    io.outReal(r.cartesian0[1]);
    io.outReal(r.cartesian1[0]);
    io.outReal(r.cartesian1[1]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

ORACLE_CASE("DistTriangle3Rectangle3.compute")
{
    // Mode 4 centers the rectangle on the first triangle vertex.
    int32_t mode = Mode(io);
    Vector3<double> unused{};
    Triangle3<double> triangle = TriangleInputs(io, mode, unused, false);
    Rectangle3<double> rectangle =
        RectangleInputs(io, mode, triangle.v[0], mode == 4);
    DCPQuery<double, Triangle3<double>, Rectangle3<double>> query{};
    auto r = query(triangle, rectangle);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.barycentric[0]);
    io.outReal(r.barycentric[1]);
    io.outReal(r.barycentric[2]);
    io.outReal(r.cartesian[0]);
    io.outReal(r.cartesian[1]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

ORACLE_CASE("DistTriangle3Triangle3.compute")
{
    // Mode 4 gives the two triangles a common vertex.
    int32_t mode = Mode(io);
    Vector3<double> unused{};
    Triangle3<double> triangle0 = TriangleInputs(io, mode, unused, false);
    Triangle3<double> triangle1 =
        TriangleInputs(io, mode, triangle0.v[0], mode == 4);
    DCPQuery<double, Triangle3<double>, Triangle3<double>> query{};
    auto r = query(triangle0, triangle1);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.barycentric0[0]);
    io.outReal(r.barycentric0[1]);
    io.outReal(r.barycentric0[2]);
    io.outReal(r.barycentric1[0]);
    io.outReal(r.barycentric1[1]);
    io.outReal(r.barycentric1[2]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

// =====================================================================
// DistOrientedBox3OrientedBox3.h, DistAlignedBox3OrientedBox3.h
// =====================================================================
namespace
{
    // center (3), the three axes (9) and the extent (3).
    OrientedBox3<double> OrientedBoxInputs(oracle::Ctx& io, int32_t mode,
        Vector3<double> const& contact, bool useContact)
    {
        Vector3<double> center{};
        if (useContact)
        {
            center = io.givenVec(contact);
        }
        else
        {
            center = Point<3>(io, mode);
        }
        std::array<Vector3<double>, 3> axis = Axes3(io, mode);
        Vector3<double> extent = Extent<3>(io, mode);
        return OrientedBox3<double>(center, axis, extent);
    }
}

ORACLE_CASE("DistOrientedBox3OrientedBox3.compute")
{
    // Mode 4 gives the two boxes a common center, so they overlap and the
    // distance is exactly zero.
    int32_t mode = Mode(io);
    Vector3<double> unused{};
    OrientedBox3<double> box0 = OrientedBoxInputs(io, mode, unused, false);
    OrientedBox3<double> box1 =
        OrientedBoxInputs(io, mode, box0.center, mode == 4);
    DCPQuery<double, OrientedBox3<double>, OrientedBox3<double>> query{};
    auto r = query(box0, box1);
    EmitTwoClosest(io, r);
}

ORACLE_CASE("DistAlignedBox3OrientedBox3.compute")
{
    // Mode 4 centers the oriented box at the center of the aligned box.
    int32_t mode = Mode(io);
    Vector3<double> lo{}, hi{};
    AlignedBoxInputs<3>(io, mode, lo, hi);
    Vector3<double> contact = 0.5 * (lo + hi);
    OrientedBox3<double> box1 = OrientedBoxInputs(io, mode, contact, mode == 4);
    AlignedBox3<double> box0(lo, hi);
    DCPQuery<double, AlignedBox3<double>, OrientedBox3<double>> query{};
    auto r = query(box0, box1);
    EmitTwoClosest(io, r);
}

// =====================================================================
// DistPoint3Tetrahedron3.h, DistTetrahedron3Tetrahedron3.h
// =====================================================================
namespace
{
    // Four vertices (12). When 'useContact' is true the first vertex is the
    // given point.
    Tetrahedron3<double> TetrahedronInputs(oracle::Ctx& io, int32_t mode,
        Vector3<double> const& contact, bool useContact)
    {
        Vector3<double> v0{};
        if (useContact)
        {
            v0 = io.givenVec(contact);
        }
        else
        {
            v0 = Point<3>(io, mode);
        }
        Vector3<double> v1 = Point<3>(io, mode);
        Vector3<double> v2 = Point<3>(io, mode);
        Vector3<double> v3 = Point<3>(io, mode);
        return Tetrahedron3<double>(v0, v1, v2, v3);
    }

    template <typename Result>
    void EmitBarycentric4Pair(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        for (size_t i = 0; i < 4; ++i) { io.outReal(r.barycentric0[i]); }
        for (size_t i = 0; i < 4; ++i) { io.outReal(r.barycentric1[i]); }
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }
}

ORACLE_CASE("DistPoint3Tetrahedron3.compute")
{
    // Mode 4 puts the query point at the centroid, which is the 'inside the
    // solid tetrahedron' branch where no face is visible. The lattice modes
    // produce degenerate (coplanar) tetrahedra, for which GetPlanes
    // normalizes a zero cross product and ComputeBarycentrics reports
    // failure and zeroes the barycentric coordinates.
    int32_t mode = Mode(io);
    Vector3<double> unused{};
    Tetrahedron3<double> tetra = TetrahedronInputs(io, mode, unused, false);
    Vector3<double> point{};
    if (mode == 4)
    {
        point = io.givenVec(tetra.ComputeCentroid());
    }
    else
    {
        point = Point<3>(io, mode);
    }
    DCPQuery<double, Vector3<double>, Tetrahedron3<double>> query{};
    auto r = query(point, tetra);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    for (size_t i = 0; i < 4; ++i) { io.outReal(r.barycentric[i]); }
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

ORACLE_CASE("DistTetrahedron3Tetrahedron3.compute")
{
    // Mode 4 places the first vertex of the second tetrahedron at the
    // centroid of the first, which reaches the 'foundZeroDistance' early exit
    // of the face-pair loop.
    int32_t mode = Mode(io);
    Vector3<double> unused{};
    Tetrahedron3<double> tetra0 = TetrahedronInputs(io, mode, unused, false);
    Vector3<double> contact = tetra0.ComputeCentroid();
    Tetrahedron3<double> tetra1 =
        TetrahedronInputs(io, mode, contact, mode == 4);
    DCPQuery<double, Tetrahedron3<double>, Tetrahedron3<double>> query{};
    auto r = query(tetra0, tetra1);
    EmitBarycentric4Pair(io, r);
}

ORACLE_CASE("DistTetrahedron3Tetrahedron3.compute.nested")
{
    // The second tetrahedron is the first shrunk toward its centroid, so no
    // pair of faces meets and the query falls through to the containment test
    // of the centroids, the only branch that reports a zero distance without
    // a zero face-pair distance.
    Vector3<double> v[4]{};
    for (int32_t i = 0; i < 4; ++i)
    {
        v[i] = io.givenVec(RawVec<3>(io, -4.0, 4.0));
    }
    Tetrahedron3<double> tetra0(v[0], v[1], v[2], v[3]);
    Vector3<double> centroid = tetra0.ComputeCentroid();
    double s = io.raw(0.05, 0.3);
    Vector3<double> w[4]{};
    for (int32_t i = 0; i < 4; ++i)
    {
        w[i] = io.givenVec(centroid + s * (v[i] - centroid));
    }
    Tetrahedron3<double> tetra1(w[0], w[1], w[2], w[3]);
    DCPQuery<double, Tetrahedron3<double>, Tetrahedron3<double>> query{};
    auto r = query(tetra0, tetra1);
    EmitBarycentric4Pair(io, r);
}

// =====================================================================
// DistPoint3Parallelepiped3.h
// =====================================================================
namespace
{
    // center (3) and three right-handed axes (9). The lattice mode rejects
    // integer triples until the basis is right-handed; the fallback is the
    // coordinate frame.
    Parallelepiped3<double> ParallelepipedInputs(oracle::Ctx& io, int32_t mode)
    {
        Vector3<double> center{};
        std::array<Vector3<double>, 3> axis{};
        if (Lattice(mode))
        {
            center = Point<3>(io, mode);
            axis[0] = Vector3<double>{ 1.0, 0.0, 0.0 };
            axis[1] = Vector3<double>{ 0.0, 1.0, 0.0 };
            axis[2] = Vector3<double>{ 0.0, 0.0, 1.0 };
            for (int32_t attempt = 0; attempt < 200; ++attempt)
            {
                Vector3<double> a0 = RawLatticeVec<3>(io, -2, 2);
                Vector3<double> a1 = RawLatticeVec<3>(io, -2, 2);
                Vector3<double> a2 = RawLatticeVec<3>(io, -2, 2);
                if (DotCross(a0, a1, a2) > 0.0)
                {
                    axis[0] = a0;
                    axis[1] = a1;
                    axis[2] = a2;
                    break;
                }
            }
        }
        else
        {
            center = Point<3>(io, mode);
            Vector3<double> a0{}, a1{}, a2{};
            for (int32_t attempt = 0; attempt < 200; ++attempt)
            {
                a0 = RawVec<3>(io, -2.0, 2.0);
                a1 = RawVec<3>(io, -2.0, 2.0);
                a2 = RawVec<3>(io, -2.0, 2.0);
                if (DotCross(a0, a1, a2) > 0.0) { break; }
                std::swap(a0, a1);
                if (DotCross(a0, a1, a2) > 0.0) { break; }
            }
            axis[0] = a0;
            axis[1] = a1;
            axis[2] = a2;
        }
        io.givenVec(axis[0]);
        io.givenVec(axis[1]);
        io.givenVec(axis[2]);
        return Parallelepiped3<double>(center, axis);
    }
}

ORACLE_CASE("DistPoint3Parallelepiped3.compute")
{
    // Mode 4 puts the query point at the center of the parallelepiped, the
    // Rzzz region where the minimizer is Z itself.
    int32_t mode = Mode(io);
    Parallelepiped3<double> ppd = ParallelepipedInputs(io, mode);
    Vector3<double> point{};
    if (mode == 4)
    {
        point = io.givenVec(ppd.center);
    }
    else
    {
        point = Point<3>(io, mode);
    }
    DCPQuery<double, Vector3<double>, Parallelepiped3<double>> query{};
    auto r = query(point, ppd);
    EmitTwoClosest(io, r);
}

ORACLE_CASE("DistPoint3Parallelepiped3.ctor.leftHanded")
{
    // Throw parity: Parallelepiped3 asserts a right-handed basis, so both
    // implementations must reject a left-handed one. Swapping two axes of a
    // right-handed triple negates the triple product.
    int32_t mode = Mode(io);
    Vector3<double> center = Point<3>(io, mode);
    Vector3<double> a0{}, a1{}, a2{};
    for (int32_t attempt = 0; attempt < 200; ++attempt)
    {
        a0 = RawVec<3>(io, -2.0, 2.0);
        a1 = RawVec<3>(io, -2.0, 2.0);
        a2 = RawVec<3>(io, -2.0, 2.0);
        if (DotCross(a0, a1, a2) < 0.0) { break; }
        std::swap(a0, a1);
        if (DotCross(a0, a1, a2) < 0.0) { break; }
    }
    std::array<Vector3<double>, 3> axis{};
    axis[0] = io.givenVec(a0);
    axis[1] = io.givenVec(a1);
    axis[2] = io.givenVec(a2);
    Vector3<double> point = Point<3>(io, mode);
    Parallelepiped3<double> ppd(center, axis);
    DCPQuery<double, Vector3<double>, Parallelepiped3<double>> query{};
    auto r = query(point, ppd);
    EmitTwoClosest(io, r);
}

ORACLE_CASE("DistPoint3Parallelepiped3.getMinimizer")
{
    // GetMinimizer is public upstream and in the port. Z cycles through the
    // 27 regions of the dispatch tree; within a region the coordinate is
    // either strictly inside, exactly on a boundary (-1 or +1, where the
    // '<= PosOne()' comparisons are evaluated at equality) or well outside.
    // A is the Gram matrix of a basis, exactly as operator() builds it.
    int32_t region = io.index() % 27;
    Vector3<double> a0 = RawVec<3>(io, -2.0, 2.0);
    Vector3<double> a1 = RawVec<3>(io, -2.0, 2.0);
    Vector3<double> a2 = RawVec<3>(io, -2.0, 2.0);
    Matrix3x3<double> B{};
    B.SetCol(0, a0);
    B.SetCol(1, a1);
    B.SetCol(2, a2);
    Matrix3x3<double> A = MultiplyATB(B, B);
    for (int32_t r = 0; r < 3; ++r)
    {
        for (int32_t c = 0; c < 3; ++c) { io.given(A(r, c)); }
    }

    Vector3<double> Z{};
    int32_t code = region;
    for (int32_t i = 0; i < 3; ++i)
    {
        int32_t digit = code % 3;
        code /= 3;
        if (digit == 0)
        {
            Z[i] = -1.0 - io.raw(0.0, 2.0);
        }
        else if (digit == 1)
        {
            int32_t kind = io.rawInteger(0, 3);
            Z[i] = (kind == 0 ? -1.0 : (kind == 1 ? 1.0 : io.raw(-1.0, 1.0)));
        }
        else
        {
            Z[i] = 1.0 + io.raw(0.0, 2.0);
        }
    }
    io.givenVec(Z);

    DCPQuery<double, Vector3<double>, Parallelepiped3<double>> query{};
    Vector3<double> K = query.GetMinimizer(A, Z);
    io.outVec(K);
}

// =====================================================================
// DistOrientedBox3Cone3.h
// =====================================================================
namespace
{
    // The box and the cone frustum, drawn unrecorded so that a rejection
    // loop can probe them before they are recorded.
    struct BoxConeInputs
    {
        Vector3<double> origin{}, direction{};
        double angle{}, tanAngle{}, minHeight{}, maxInc{};
        Vector3<double> center{}, extent{};
        std::array<Vector3<double>, 3> axis{};
    };

    std::array<Vector3<double>, 3> RawAxes3(oracle::Ctx& io, int32_t mode)
    {
        std::array<Vector3<double>, 3> axis{};
        if (Lattice(mode))
        {
            int32_t k = io.rawInteger(0, 2);
            int32_t s = io.rawInteger(0, 1);
            axis[0].MakeUnit(k);
            axis[1].MakeUnit((k + 1) % 3);
            if (s == 1) { axis[0] = -axis[0]; axis[1] = -axis[1]; }
            axis[2] = Cross(axis[0], axis[1]);
        }
        else
        {
            axis[0] = RawUnit<3>(io);
            double len = 0.0;
            do
            {
                axis[1] = RawVec<3>(io, -1.0, 1.0);
                axis[1] = axis[1] - Dot(axis[1], axis[0]) * axis[0];
                len = Length(axis[1]);
            } while (len < 0.25);
            Normalize(axis[1]);
            axis[2] = Cross(axis[0], axis[1]);
        }
        return axis;
    }

    BoxConeInputs RawBoxCone(oracle::Ctx& io, int32_t mode)
    {
        BoxConeInputs b{};
        b.origin = RawVec<3>(io, -3.0, 3.0);
        b.direction = RawUnit<3>(io);
        b.angle = io.raw(0.15, 1.25);
        Ray3<double> ray(b.origin, b.direction);
        Cone3<double> probe(ray, b.angle);
        b.tanAngle = probe.tanAngle;
        b.minHeight = io.raw(0.0, 2.0);
        b.maxInc = io.raw(0.5, 4.0);
        if (mode == 4)
        {
            double mid = 0.5 * (b.minHeight + (b.minHeight + b.maxInc));
            b.center = b.origin + mid * b.direction;
        }
        else if (Lattice(mode))
        {
            b.center = RawLatticeVec<3>(io, -3, 3);
        }
        else
        {
            b.center = RawVec<3>(io, -5.0, 5.0);
        }
        b.axis = RawAxes3(io, mode);
        b.extent = Lattice(mode) ? RawLatticeVec<3>(io, 0, 2)
            : RawVec<3>(io, 0.25, 3.0);
        return b;
    }

    // origin (3), direction (3), angle, Cone::SetAngle's tanAngle, the
    // minimum height and the increment to the maximum height, then the box
    // center (3), axes (9) and extent (3). tanAngle is recorded because
    // Cone::SetAngle computes it with std::tan, which belongs to Cone.h
    // rather than to this group; the replay puts the recorded value into the
    // port's cone so that a libm difference of another group's header stays
    // out of this comparison. Everything else the query reads from the cone
    // is recorded directly.
    void EmitBoxConeInputs(oracle::Ctx& io, BoxConeInputs const& b)
    {
        io.givenVec(b.origin);
        io.givenVec(b.direction);
        io.given(b.angle);
        io.given(b.tanAngle);
        io.given(b.minHeight);
        io.given(b.maxInc);
        io.givenVec(b.center);
        io.givenVec(b.axis[0]);
        io.givenVec(b.axis[1]);
        io.givenVec(b.axis[2]);
        io.givenVec(b.extent);
    }

    Cone3<double> BuildCone(BoxConeInputs const& b, bool finite)
    {
        Ray3<double> ray(b.origin, b.direction);
        if (finite)
        {
            return Cone3<double>(ray, b.angle, b.minHeight,
                b.minHeight + b.maxInc);
        }
        return Cone3<double>(ray, b.angle, b.minHeight);
    }

    OrientedBox3<double> BuildBox(BoxConeInputs const& b)
    {
        return OrientedBox3<double>(b.center, b.axis, b.extent);
    }
}

namespace
{
    // ---- probe: which configurations can be compared bit for bit? --------
    //
    // DistOrientedBox3Cone3 minimizes F(angle), the distance between the box
    // and the planar quadrilateral cut from the frustum at that angle, with
    // Minimize1. Two things make the default search incomparable:
    //
    //  * The port fixes two defects of Minimize1::GetBracketedMinimum
    //    (UPSTREAM-FINDINGS Minimize1.h item 1, issue #298), so the two
    //    searches visit different angles. This is the '.deviation' case.
    //  * F calls std::cos and std::sin, and MSVC's and V8's differ by one
    //    ulp on a couple of percent of their arguments -- including
    //    sin(+/-pi/4), one of the angles the subdivision visits. The search
    //    then compares perturbed values and can settle on a different local
    //    minimum, so the difference is not a rounding difference that a
    //    tolerance could cover.
    //
    // The '.control' case removes both by exercising the Control overload
    // with maxSubdivisions = 1, maxBisections = 1, epsilon = 10 and
    // tolerance = 0, and by keeping only configurations whose initial
    // polyline {(-pi/2,f0),(0,fm),(pi/2,f1)} is V-shaped. GetMinimum then
    // evaluates F at exactly -pi/2, 0 and +pi/2 and calls
    // GetBracketedMinimum, whose first act is the convergence test
    // 'pi <= 2*0*|tm| + 10', so it returns without evaluating F again. The
    // answer is the smallest of three box-quadrilateral distances computed
    // at three fixed angles at which the two libm implementations agree bit
    // for bit (cos(0) = 1, sin(0) = 0, sin(+/-pi/2) = +/-1 and
    // cos(+/-pi/2) = 6.123233995736766e-17 in both). Everything the case
    // measures is then DistOrientedBox3Cone3's own arithmetic: the 5x5
    // quadratic form, the 10-D LCP and the mapping of the LCP solution back
    // to the two closest points.
    //
    // The probe below is upstream's DoBoxQuadQuery, copied verbatim, used to
    // evaluate the V-shape condition before the inputs are recorded.
    double ProbeF(LCPSolver<double, 10>& lcp, OrientedBox3<double> const& box,
        Cone3<double> const& cone, Vector3<double> const& coneW0,
        Vector3<double> const& coneW1, double quadAngle)
    {
        double const zero = 0.0, one = 1.0, two = 2.0;

        Vector3<double> K = box.center, ell{};
        for (int32_t i = 0; i < 3; ++i)
        {
            K -= box.extent[i] * box.axis[i];
            ell[i] = two * box.extent[i];
        }

        double cs = std::cos(quadAngle), sn = std::sin(quadAngle);
        Vector3<double> term = cone.tanAngle * (cs * coneW0 + sn * coneW1);
        std::array<Vector3<double>, 2> G{};
        G[0] = cone.ray.direction - term;
        G[1] = cone.ray.direction + term;

        Matrix<5, 5, double> A{};
        A(0, 0) = one;
        A(0, 1) = zero;
        A(0, 2) = zero;
        A(0, 3) = -Dot(box.axis[0], G[0]);
        A(0, 4) = -Dot(box.axis[0], G[1]);
        A(1, 0) = A(0, 1);
        A(1, 1) = one;
        A(1, 2) = zero;
        A(1, 3) = -Dot(box.axis[1], G[0]);
        A(1, 4) = -Dot(box.axis[1], G[1]);
        A(2, 0) = A(0, 2);
        A(2, 1) = A(1, 2);
        A(2, 2) = one;
        A(2, 3) = -Dot(box.axis[2], G[0]);
        A(2, 4) = -Dot(box.axis[2], G[1]);
        A(3, 0) = A(0, 3);
        A(3, 1) = A(1, 3);
        A(3, 2) = A(2, 3);
        A(3, 3) = Dot(G[0], G[0]);
        A(3, 4) = Dot(G[0], G[1]);
        A(4, 0) = A(0, 4);
        A(4, 1) = A(1, 4);
        A(4, 2) = A(2, 4);
        A(4, 3) = A(3, 4);
        A(4, 4) = Dot(G[1], G[1]);

        Vector3<double> KmV = K - cone.ray.origin;
        Vector<5, double> b{};
        b[0] = Dot(box.axis[0], KmV);
        b[1] = Dot(box.axis[1], KmV);
        b[2] = Dot(box.axis[2], KmV);
        b[3] = -Dot(G[0], KmV);
        b[4] = -Dot(G[1], KmV);

        Matrix<5, 5, double> D{};
        D(0, 0) = -one;
        D(1, 1) = -one;
        D(2, 2) = -one;
        D(3, 3) = one;
        D(3, 4) = one;
        D(4, 3) = -one;
        D(4, 4) = -one;

        Vector<5, double> e{};
        e[0] = -ell[0];
        e[1] = -ell[1];
        e[2] = -ell[2];
        e[3] = cone.GetMinHeight();
        e[4] = -cone.GetMaxHeight();

        std::array<double, 10> q;
        for (int32_t i = 0, ip5 = 5; i < 5; ++i, ++ip5)
        {
            q[i] = b[i];
            q[ip5] = -e[i];
        }

        std::array<std::array<double, 10>, 10> M;
        for (int32_t r = 0, rp5 = 5; r < 5; ++r, ++rp5)
        {
            for (int32_t c = 0, cp5 = 5; c < 5; ++c, ++cp5)
            {
                M[r][c] = A(r, c);
                M[rp5][c] = D(r, c);
                M[r][cp5] = -D(c, r);
                M[rp5][cp5] = zero;
            }
        }

        std::array<double, 10> w, z;
        if (lcp.Solve(q, M, w, z))
        {
            Vector3<double> boxClosestPoint = K;
            for (int32_t i = 0; i < 3; ++i)
            {
                boxClosestPoint += z[i] * box.axis[i];
            }

            Vector3<double> quadClosestPoint = cone.ray.origin;
            for (int32_t i = 0, ip3 = 3; i < 2; ++i, ++ip3)
            {
                quadClosestPoint += z[ip3] * G[i];
            }

            return Length(boxClosestPoint - quadClosestPoint);
        }
        return std::numeric_limits<double>::max();
    }

    // True when Minimize1::GetMinimum's initial polyline is V-shaped, which
    // is upstream's own test, evaluated on upstream's own F.
    bool ProbeVShaped(OrientedBox3<double> const& box, Cone3<double> const& cone)
    {
        LCPSolver<double, 10> lcp{};
        std::array<Vector3<double>, 3> basis{};
        basis[0] = cone.ray.direction;
        ComputeOrthogonalComplement(1, basis.data());
        Vector3<double> coneW0 = basis[1];
        Vector3<double> coneW1 = basis[2];

        double t0 = static_cast<double>(-GTE_C_HALF_PI);
        double t1 = static_cast<double>(+GTE_C_HALF_PI);
        double tInitial = 0.5 * (t0 + t1);
        double f0 = ProbeF(lcp, box, cone, coneW0, coneW1, t0);
        double fInitial = ProbeF(lcp, box, cone, coneW0, coneW1, tInitial);
        double f1 = ProbeF(lcp, box, cone, coneW0, coneW1, t1);
        return ((fInitial < f0) && (f1 >= fInitial)) ||
            ((f1 > fInitial) && (f0 >= fInitial));
    }

    // Mode 4 centers the box on the axis of the frustum halfway between the
    // two heights, so the box meets the frustum and the minimum distance is
    // attained in the interior of the angle sweep. When 'control' is not
    // null the generator keeps only the configurations whose initial polyline
    // is V-shaped, which is what confines the search to the three angles at
    // which the two libm implementations agree (see the probe above).
    void BoxCone(oracle::Ctx& io,
        DCPQuery<double, OrientedBox3<double>, Cone3<double>>::Control const*
        control)
    {
        int32_t mode = Mode(io);
        BoxConeInputs b{};
        bool haveFallback = false, accepted = (control == nullptr);
        BoxConeInputs fallback{};
        for (int32_t attempt = 0; attempt < 2000; ++attempt)
        {
            b = RawBoxCone(io, mode);
            if (accepted) { break; }
            fallback = b;
            haveFallback = true;
            if (ProbeVShaped(BuildBox(b), BuildCone(b, true)))
            {
                accepted = true;
                break;
            }
        }
        if (!accepted && haveFallback) { b = fallback; }

        EmitBoxConeInputs(io, b);
        OrientedBox3<double> box = BuildBox(b);
        Cone3<double> cone = BuildCone(b, true);
        DCPQuery<double, OrientedBox3<double>, Cone3<double>> query{};
        auto r = query(box, cone, control);
        io.outReal(r.distance);
        io.outVec(r.boxClosestPoint);
        io.outVec(r.coneClosestPoint);
    }
}

ORACLE_CASE("DistOrientedBox3Cone3.compute.control")
{
    // The Control overload of operator(), with controls that stop the search
    // at its three initial samples. See the probe above for why that is what
    // makes this query comparable at all.
    DCPQuery<double, OrientedBox3<double>, Cone3<double>>::Control control(
        1, 1, 10.0, 0.0);
    BoxCone(io, &control);
}

ORACLE_CASE("DistOrientedBox3Cone3.compute.deviation")
{
    // The default controls. Minimize1::GetBracketedMinimum collapses its
    // bracket whenever the parabola vertex lands within round-off of the
    // middle sample, which for this query is the normal situation: the two
    // endpoint slices of the angle sweep are the same point set, so
    // F(-pi/2) and F(+pi/2) agree to a few ulps and the exact 'tv == tm'
    // test is missed. The port compares them to within the resolution of the
    // bracket and continues the search where upstream returns
    // (UPSTREAM-FINDINGS DistOrientedBox3Cone3.h / Minimize1.h item 1, issue
    // #298).
    BoxCone(io, nullptr);
}

ORACLE_CASE("DistOrientedBox3Cone3.compute.infinite")
{
    // The header documents the precondition hmax < infinity but upstream
    // never checks it. Cone encodes an infinite cone as maxHeight = -1, which
    // upstream feeds straight into the LCP as the constraint z3 + z4 <= -1,
    // an infeasible program. The port asserts the precondition instead
    // (UPSTREAM-FINDINGS DistOrientedBox3Cone3.h item 2, issue #298), so it
    // throws on every record of this case.
    int32_t mode = Mode(io);
    BoxConeInputs b = RawBoxCone(io, mode);
    EmitBoxConeInputs(io, b);
    OrientedBox3<double> box = BuildBox(b);
    Cone3<double> cone = BuildCone(b, false);
    DCPQuery<double, OrientedBox3<double>, Cone3<double>> query{};
    auto r = query(box, cone);
    io.outReal(r.distance);
    io.outVec(r.boxClosestPoint);
    io.outVec(r.coneClosestPoint);
}

// =====================================================================
// DistCircle3Circle3.h
// =====================================================================
namespace
{
    using Circle3Query = DCPQuery<double, Circle3<double>, Circle3<double>>;

    // The closed form of DistPoint3Circle3, used only by the independent
    // reference below.
    double PointCircleDistance(Vector3<double> const& P, Circle3<double> const& c)
    {
        Vector3<double> delta = P - c.center;
        double h = Dot(c.normal, delta);
        Vector3<double> inPlane = delta - h * c.normal;
        double len = Length(inPlane);
        double radial = len - c.radius;
        return std::sqrt(h * h + radial * radial);
    }

    // An independent reference for the circle-circle distance: the minimum
    // over circle1 of the point-to-circle0 distance. The objective is smooth
    // and bounded in the circle angle, so a dense scan plus a golden-section
    // refinement resolves it to near machine precision. It is used only to
    // decide which records a case keeps.
    double MinDistanceBetweenCircles(Circle3<double> const& c0,
        Circle3<double> const& c1)
    {
        double const twoPi = 6.283185307179586476925286766559;
        Vector3<double> const& n = c1.normal;
        Vector3<double> seed = (std::fabs(n[0]) < 0.9
            ? Vector3<double>{ 1.0, 0.0, 0.0 }
            : Vector3<double>{ 0.0, 1.0, 0.0 });
        Vector3<double> u = seed - Dot(seed, n) * n;
        Normalize(u);
        Vector3<double> w = Cross(n, u);
        auto f = [&](double a)
        {
            Vector3<double> q = c1.center + (c1.radius * std::cos(a)) * u
                + (c1.radius * std::sin(a)) * w;
            return PointCircleDistance(q, c0);
        };

        int32_t const samples = 4096;
        double best = f(0.0);
        int32_t bestIndex = 0;
        for (int32_t i = 1; i < samples; ++i)
        {
            double y = f(twoPi * i / samples);
            if (y < best) { best = y; bestIndex = i; }
        }
        double h = twoPi / samples;
        double lo = twoPi * bestIndex / samples - h;
        double hi = lo + 2.0 * h;
        double const phi = 0.61803398874989484820458683436564;
        for (int32_t i = 0; i < 200; ++i)
        {
            double m0 = hi - phi * (hi - lo);
            double m1 = lo + phi * (hi - lo);
            if (f(m0) <= f(m1)) { hi = m1; } else { lo = m0; }
        }
        double mid = f(0.5 * (lo + hi));
        return (mid < best ? mid : best);
    }

    // The relative residual between upstream's reported distance and the
    // independent reference, which resolves the minimum to about 1e-13
    // relative (4096 samples followed by 200 golden-section steps on a
    // smooth objective).
    //
    // A case that wants sound inputs keeps the records below 1e-12. That
    // bound is what makes the 1e-12 comparison tolerance of the main cases
    // meaningful: on an accepted record upstream's own answer is right to
    // 1e-12, the port's is at least as good, so a difference bigger than
    // that is a finding rather than the conditioning of the substituted
    // expression. A deviation case keeps only the records above 1e-8, four
    // orders of magnitude above the reference's own resolution, so that
    // resolution cannot be mistaken for an upstream defect.
    double CirclesUpstreamResidual(Circle3<double> const& c0,
        Circle3<double> const& c1, double distance)
    {
        double reference = MinDistanceBetweenCircles(c0, c1);
        double scale = std::fmax(1.0, std::fmax(std::fabs(distance),
            std::fabs(reference)));
        return std::fabs(distance - reference) / scale;
    }

    struct Circle3Pair
    {
        Vector3<double> center0{}, normal0{}, center1{}, normal1{};
        double radius0{}, radius1{};
    };

    void EmitCirclePair(oracle::Ctx& io, Circle3Pair const& p)
    {
        io.givenVec(p.center0);
        io.givenVec(p.normal0);
        io.given(p.radius0);
        io.givenVec(p.center1);
        io.givenVec(p.normal1);
        io.given(p.radius1);
    }

    Circle3<double> Circle0(Circle3Pair const& p)
    {
        return Circle3<double>(p.center0, p.normal0, p.radius0);
    }

    Circle3<double> Circle1(Circle3Pair const& p)
    {
        return Circle3<double>(p.center1, p.normal1, p.radius1);
    }

    // numClosestPairs, distance, sqrDistance, equidistant and then the pairs.
    // 'firstPairOnly' emits pair 0 alone and omits the count, for the cases
    // where the port's extra closest-point candidates can make a tie visible
    // that upstream's candidate set does not contain (see the report).
    void EmitCircle3Circle3(oracle::Ctx& io, Circle3Query::Result const& r,
        bool firstPairOnly)
    {
        if (!firstPairOnly) { io.outInt(r.numClosestPairs); }
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outBool(r.equidistant);
        size_t n = firstPairOnly ? 1 : r.numClosestPairs;
        for (size_t j = 0; j < n; ++j)
        {
            io.outVec(r.circle0Closest[j]);
            io.outVec(r.circle1Closest[j]);
        }
    }

    enum class CircleFamily
    {
        // Generic normals, upstream's answer agrees with the reference.
        generic,
        // A mirror-symmetric configuration, built in the frame PrepareCircles
        // produces so that the symmetry is exact: circle1 is the unit circle
        // in the plane z = 0, circle0's center is on the plane x = 0 and its
        // normal has an exactly zero x-component. The reflection x -> -x maps
        // each circle to itself, so every critical angle comes in a pair with
        // the same distance and phi has a double root there. A floating-point
        // bisection cannot separate a double root, the sign-change search
        // reports neither half, and upstream's sn = -p6(cs)/p7(cs) evaluated
        // at whatever other root it did find lands off the unit circle
        // (issue #331). Random sampling never reaches this: the symmetry has
        // to be constructed.
        mirrorSymmetric,
        // Both normals are exactly +/-(0,0,1), so PrepareCircles' rotation is
        // the identity and the port's exact common normal is upstream's
        // rounded one. The five branches of DoQueryParallelPlanes are then
        // compared bit for bit.
        parallelPlanes,
        // Anti-parallel normals orthogonal to the z-axis.
        antiParallel,
        // Coaxial circles in parallel planes.
        coaxial,
        // The larger circle's normal is exactly (0,0,1) and circle0's normal
        // is generic, so PrepareCircles' acos is acos(1) = 0 and Rotation of
        // the zero axis-angle is exactly the identity matrix. No libm
        // function other than sqrt is then evaluated anywhere on the path, so
        // the residual of this family is the port's sn substitution alone.
        identityRotation
    };

    Circle3Pair RawCirclePair(oracle::Ctx& io, CircleFamily family,
        int32_t mode)
    {
        Circle3Pair p{};
        if (family == CircleFamily::parallelPlanes)
        {
            Vector3<double> n{ 0.0, 0.0, 1.0 };
            p.normal0 = (io.rawInteger(0, 1) == 0 ? -n : n);
            p.normal1 = (io.rawInteger(0, 1) == 0 ? -n : n);
            if (Lattice(mode))
            {
                p.center0 = RawLatticeVec<3>(io, -3, 3);
                p.center1 = RawLatticeVec<3>(io, -3, 3);
                p.radius0 = static_cast<double>(io.rawInteger(1, 4));
                p.radius1 = static_cast<double>(io.rawInteger(1, 4));
            }
            else
            {
                p.center0 = RawVec<3>(io, -4.0, 4.0);
                p.center1 = RawVec<3>(io, -4.0, 4.0);
                p.radius0 = io.raw(0.25, 4.0);
                p.radius1 = io.raw(0.25, 4.0);
            }
            if (mode == 4)
            {
                // Concentric circles in the same plane: the 'equidistant'
                // branches of DoQueryParallelPlanes.
                p.center1 = p.center0;
            }
        }
        else if (family == CircleFamily::mirrorSymmetric)
        {
            p.center1 = Vector3<double>{ 0.0, 0.0, 0.0 };
            p.normal1 = Vector3<double>{ 0.0, 0.0, 1.0 };
            p.radius1 = io.raw(1.0, 4.0);
            // circle0's center is offset in +y and its plane is close to
            // perpendicular to circle1's, which puts the global minimum on
            // the mirror axis itself: the closest point of circle1 is the one
            // with the largest y, cs = 0. That is exactly the critical angle
            // whose root of phi has even multiplicity, so upstream cannot see
            // it and the two implementations differ on nearly every draw.
            p.center0 = Vector3<double>{ 0.0, io.raw(1.0, 3.0),
                io.raw(-1.0, 1.0) };
            double t = io.raw(1.2, 1.94);
            p.normal0 = Vector3<double>{ 0.0, std::sin(t), std::cos(t) };
            p.radius0 = io.raw(0.25, p.radius1);
        }
        else if (family == CircleFamily::identityRotation)
        {
            p.normal1 = Vector3<double>{ 0.0, 0.0, 1.0 };
            p.normal0 = RawUnit<3>(io);
            if (Lattice(mode))
            {
                p.center0 = RawLatticeVec<3>(io, -3, 3);
                p.center1 = RawLatticeVec<3>(io, -3, 3);
                p.radius1 = static_cast<double>(io.rawInteger(2, 4));
                p.radius0 = static_cast<double>(io.rawInteger(1, 2));
            }
            else
            {
                p.center0 = RawVec<3>(io, -4.0, 4.0);
                p.center1 = RawVec<3>(io, -4.0, 4.0);
                p.radius1 = io.raw(1.0, 4.0);
                p.radius0 = io.raw(0.25, p.radius1);
            }
            if (mode == 4)
            {
                // The x-components of the two centers agree exactly, so the
                // transformed circle0 center already has a zero x-component
                // and the second rotation of PrepareCircles is skipped too.
                p.center0 = Vector3<double>{ p.center1[0],
                    p.center1[1] + io.raw(0.5, 3.0),
                    p.center1[2] + io.raw(-3.0, 3.0) };
            }
        }
        else if (family == CircleFamily::antiParallel)
        {
            int32_t k = io.rawInteger(0, 1);
            Vector3<double> n{ 0.0, 0.0, 0.0 };
            n[k] = 1.0;
            p.normal0 = n;
            p.normal1 = -n;
            p.center0 = RawVec<3>(io, -4.0, 4.0);
            p.center1 = RawVec<3>(io, -4.0, 4.0);
            p.radius0 = io.raw(0.25, 4.0);
            p.radius1 = io.raw(0.25, 4.0);
        }
        else if (family == CircleFamily::coaxial)
        {
            Vector3<double> n = RawUnit<3>(io);
            p.normal0 = n;
            p.normal1 = n;
            p.center0 = RawVec<3>(io, -4.0, 4.0);
            p.center1 = p.center0 + io.raw(-4.0, 4.0) * n;
            p.radius0 = io.raw(0.25, 4.0);
            p.radius1 = io.raw(0.25, 4.0);
        }
        else
        {
            p.normal0 = RawUnit<3>(io);
            p.normal1 = RawUnit<3>(io);
            p.center0 = RawVec<3>(io, -4.0, 4.0);
            p.center1 = RawVec<3>(io, -4.0, 4.0);
            p.radius0 = io.raw(0.25, 4.0);
            p.radius1 = io.raw(0.25, 4.0);
        }
        return p;
    }

    // 'selector' says how a candidate configuration is accepted:
    //   sound        keep it when upstream's own distance agrees with the
    //                independent reference to 1e-12 relative;
    //   wrongAnswer  keep it when upstream's own distance is wrong by more
    //                than 1e-8 relative, or upstream throws;
    //   byFamily     keep the first candidate. The construction itself is
    //                the defect condition: an anti-parallel pair of normals
    //                orthogonal to the z-axis, or a coaxial pair, always
    //                reaches the misclassified branch. The distance is not a
    //                usable selector for those two, because upstream's
    //                *distance* is often right there while its closest
    //                points are not.
    enum class CircleSelector { sound, wrongAnswer, byFamily };

    void Circle3Circle3(oracle::Ctx& io, CircleFamily family,
        CircleSelector selector, bool firstPairOnly)
    {
        int32_t mode = Mode(io);
        bool wantSound = (selector == CircleSelector::sound);
        Circle3Query query{};
        Circle3Pair p{}, fallback{};
        bool haveFallback = false, accepted = false;
        double bestResidual = -1.0;
        // A deviation case needs few attempts because its fallback is the
        // worst candidate the search saw; the cap keeps the deep run
        // tractable, since the rational root finder is slow on the
        // even-multiplicity roots these families produce.
        int32_t maxAttempts = (selector == CircleSelector::byFamily ? 1
            : (wantSound ? 200 : 24));
        for (int32_t attempt = 0; attempt < maxAttempts; ++attempt)
        {
            if (selector == CircleSelector::byFamily)
            {
                p = RawCirclePair(io, family, mode);
                accepted = true;
                if (family == CircleFamily::coaxial && mode == 4)
                {
                    // Mode 4 of the coaxial family aims at the other half of
                    // the defect: p6 and p7 vanish identically and upstream
                    // throws "Unexpected degree for p6" instead of answering.
                    // Whether the transformed normal lands exactly on (0,0,1)
                    // decides between the two, and only sampling finds out.
                    accepted = false;
                    for (int32_t k = 0; k < 40; ++k)
                    {
                        try
                        {
                            auto probe = query(Circle0(p), Circle1(p));
                            (void)probe;
                        }
                        catch (std::exception const&)
                        {
                            accepted = true;
                            break;
                        }
                        p = RawCirclePair(io, family, mode);
                    }
                }
                break;
            }
            p = RawCirclePair(io, family, mode);
            bool ok = false;
            try
            {
                auto probe = query(Circle0(p), Circle1(p));
                double residual = CirclesUpstreamResidual(Circle0(p),
                    Circle1(p), probe.distance);
                ok = wantSound ? (residual <= 1e-12) : (residual > 1e-8);
                // The fallback is the best candidate the search saw rather
                // than the last one: the least defective one for a case that
                // wants a sound upstream answer, the most defective one for a
                // deviation case.
                bool better = wantSound
                    ? (!haveFallback || residual < bestResidual)
                    : (residual > bestResidual);
                if (better)
                {
                    bestResidual = residual;
                    fallback = p;
                    haveFallback = true;
                }
            }
            catch (std::exception const&)
            {
                // Upstream's "Unexpected degree" assertion. That is itself a
                // defect the port fixes, so it is kept by a deviation case
                // and never by a case that wants a sound upstream answer.
                ok = !wantSound;
                if (!wantSound)
                {
                    fallback = p;
                    haveFallback = true;
                }
            }
            if (ok) { accepted = true; break; }
        }
        if (!accepted && haveFallback) { p = fallback; }

        EmitCirclePair(io, p);
        auto r = query(Circle0(p), Circle1(p));
        EmitCircle3Circle3(io, r, firstPairOnly);
    }
}

ORACLE_CASE("DistCircle3Circle3.compute")
{
    // Generic circles on which upstream's own answer agrees with the
    // independent reference. The port replaces upstream's
    // sn = -p6(cs)/p7(cs) with both signs of sqrt(1 - cs^2) and adds further
    // closest-point candidates (issues #331 and #431), so the reported
    // closest points are not bit-identical even here; see the report for the
    // measured residual.
    Circle3Circle3(io, CircleFamily::generic, CircleSelector::sound, false);
}

ORACLE_CASE("DistCircle3Circle3.compute.identityRotation")
{
    // The general polynomial path with PrepareCircles' rotation forced to the
    // identity: the larger circle's normal is exactly (0,0,1), so
    // std::acos(1) is 0 and Rotation of the zero axis-angle is the identity
    // matrix in both implementations. No libm function other than sqrt is
    // evaluated on the path, so the residual of this case isolates the port's
    // substitution of sn = +/-sqrt(1 - cs^2) for upstream's
    // sn = -p6(cs)/p7(cs) (issue #331). Putting upstream's expression and
    // upstream's candidate set back into src/DistCircle3Circle3.ts makes this
    // case and '.compute' bit-identical on all 20 records, which is the proof
    // that nothing else on the path deviates.
    Circle3Circle3(io, CircleFamily::identityRotation, CircleSelector::sound,
        false);
}

ORACLE_CASE("DistCircle3Circle3.compute.doubleRoot.deviation")
{
    // A constructed mirror-symmetric configuration, where phi has a double
    // root that a floating-point sign-change bisection cannot find and
    // upstream's sn = -p6(cs)/p7(cs) lands off the unit circle (issue #331).
    // Kept only when upstream's own reported distance is wrong by more than
    // 1e-8 relative to the independent reference, so that the resolution of
    // the reference cannot be mistaken for the defect; a record that never
    // reaches the threshold keeps the worst candidate the search saw.
    Circle3Circle3(io, CircleFamily::mirrorSymmetric,
        CircleSelector::wrongAnswer, true);
}

ORACLE_CASE("DistCircle3Circle3.compute.parallelPlanes")
{
    // Both normals are exactly +/-(0,0,1), so PrepareCircles' rotation is the
    // identity, the transformed circle0.normal is exactly (0,0,1) and the
    // port's withCommonNormal is a no-op: the two implementations run the
    // same DoQueryParallelPlanes on the same values. Mode 4 makes the circles
    // concentric, which is the 'equidistant' branch; the lattice modes give
    // exact tangency and exact containment.
    Circle3Circle3(io, CircleFamily::parallelPlanes, CircleSelector::sound,
        false);
}

ORACLE_CASE("DistCircle3Circle3.compute.antiParallel.deviation")
{
    // Anti-parallel normals orthogonal to the z-axis. PrepareCircles aligns
    // normals by negating any whose z-component is negative, which does
    // nothing here, so the transformed circle0.normal is (0,0,-1), the test
    // 'circle0.normal[2] < 1' misclassifies two circles in parallel planes as
    // skew and the general polynomial path answers the wrong question
    // (issue #431). The port negates the transformed normal.
    Circle3Circle3(io, CircleFamily::antiParallel, CircleSelector::byFamily,
        true);
}

ORACLE_CASE("DistCircle3Circle3.compute.coaxial.deviation")
{
    // Coaxial circles in parallel planes. Whenever the rotation leaves the
    // transformed circle0.normal a few ulps below (0,0,1), upstream takes the
    // polynomial path, where a1 = 0 makes phi and p6 perfect squares: the
    // sign-change bisection finds no root, upstream reads a
    // default-constructed candidate and reports distance 0 at the origin, or
    // p6 and p7 vanish identically and it reports "Unexpected degree for p6"
    // (issue #431). The port answers with the parallel-planes code.
    Circle3Circle3(io, CircleFamily::coaxial, CircleSelector::byFamily, true);
}
