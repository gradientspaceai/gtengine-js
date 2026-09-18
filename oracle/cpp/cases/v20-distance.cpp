// Verify group 20 (distance): differential cases for the Dist* headers that
// group 20 covers. See ORACLE.md.
//
// Generator modes are selected from io.index() so that every case mixes
// uniform random inputs, small-lattice inputs (exact integer arithmetic, which
// produces parallel, touching, tangent and coincident configurations) and
// constructions aimed at a specific branch. Every generator records exactly
// the same number of doubles in every mode, so the TypeScript replay reads the
// inputs without knowing which mode produced them.
#define ORACLE_FAMILY "v20-distance"
#include "Oracle.h"

#include <Mathematics/DistLine2Arc2.h>
#include <Mathematics/DistLine2OrientedBox2.h>
#include <Mathematics/DistLine3AlignedBox3.h>
#include <Mathematics/DistLine3Circle3.h>
#include <Mathematics/DistLine3OrientedBox3.h>
#include <Mathematics/DistLine3Rectangle3.h>
#include <Mathematics/DistLine3Triangle3.h>
#include <Mathematics/DistPoint2Arc2.h>
#include <Mathematics/DistPoint3Circle3.h>
#include <Mathematics/DistPoint3Cylinder3.h>
#include <Mathematics/DistPoint3Frustum3.h>
#include <Mathematics/DistPointAlignedBox.h>
#include <Mathematics/DistPointOrientedBox.h>
#include <Mathematics/DistRay2Circle2.h>
#include <Mathematics/DistRay2Triangle2.h>
#include <Mathematics/DistRay3CanonicalBox3.h>
#include <Mathematics/DistSegment2Circle2.h>
#include <Mathematics/DistSegment2Triangle2.h>
#include <Mathematics/DistSegment3CanonicalBox3.h>

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
    // mode 4:   constructed (parallel directions, coincident points, ...)
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

    template <int32_t N>
    Vector<N, double> Direction(oracle::Ctx& io, int32_t mode)
    {
        if (Lattice(mode))
        {
            return io.latticeDir<N>(-2, 2);
        }
        return io.unit<N>();
    }

    // Extents are nonnegative; a zero extent (a degenerate box) is allowed by
    // the headers and is generated in the lattice modes.
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

    // An orthonormal pair of 2D box/rectangle axes. The lattice modes use the
    // coordinate axes so that the box frame is the world frame and lattice
    // line directions keep their exactly zero components (the DoQuery1D and
    // DoQuery0D dispatch arms).
    std::array<Vector2<double>, 2> Axes2(oracle::Ctx& io, int32_t mode)
    {
        Vector2<double> a0{};
        if (Lattice(mode))
        {
            a0 = Vector2<double>{ 1.0, 0.0 };
        }
        else
        {
            a0 = RawUnit<2>(io);
        }
        Vector2<double> a1{ -a0[1], a0[0] };
        std::array<Vector2<double>, 2> axis{};
        axis[0] = io.givenVec(a0);
        axis[1] = io.givenVec(a1);
        return axis;
    }

    // An orthonormal triple of 3D axes, coordinate axes in the lattice modes.
    std::array<Vector3<double>, 3> Axes3(oracle::Ctx& io, int32_t mode)
    {
        Vector3<double> a0{}, a1{}, a2{};
        if (Lattice(mode))
        {
            a0 = Vector3<double>{ 1.0, 0.0, 0.0 };
            a1 = Vector3<double>{ 0.0, 1.0, 0.0 };
            a2 = Vector3<double>{ 0.0, 0.0, 1.0 };
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

    // A box in min/size form; the maximum corner is derived so that it is
    // never below the minimum corner.
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
    void EmitTwoClosest3(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    // distance, sqrDistance, parameter and the two closest points.
    template <typename Result, int32_t N>
    void EmitOneParameter(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    // The shared line/ray/segment versus circle result in 2D.
    template <typename Result>
    void EmitCircle2Pairs(oracle::Ctx& io, Result const& r)
    {
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

    template <typename Result>
    void EmitTriangle2(oracle::Ctx& io, Result const& r)
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

// =====================================================================
// DistPointAlignedBox.h
// =====================================================================
namespace
{
    template <int32_t N>
    void PointAlignedBox(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto point = Point<N>(io, mode);
        Vector<N, double> lo{}, hi{};
        AlignedBoxInputs<N>(io, mode, lo, hi);
        AlignedBox<N, double> box(lo, hi);
        DCPQuery<double, Vector<N, double>, AlignedBox<N, double>> query{};
        auto r = query(point, box);
        EmitTwoClosest3(io, r);
    }
}

ORACLE_CASE("DistPointAlignedBox.compute.2d") { PointAlignedBox<2>(io); }
ORACLE_CASE("DistPointAlignedBox.compute.3d") { PointAlignedBox<3>(io); }

// =====================================================================
// DistPointOrientedBox.h
// =====================================================================
namespace
{
    void OrientedBoxInputs2(oracle::Ctx& io, int32_t mode,
        Vector2<double>& center, std::array<Vector2<double>, 2>& axis,
        Vector2<double>& extent)
    {
        center = Point<2>(io, mode);
        axis = Axes2(io, mode);
        extent = Extent<2>(io, mode);
    }

    void OrientedBoxInputs3(oracle::Ctx& io, int32_t mode,
        Vector3<double>& center, std::array<Vector3<double>, 3>& axis,
        Vector3<double>& extent)
    {
        center = Point<3>(io, mode);
        axis = Axes3(io, mode);
        extent = Extent<3>(io, mode);
    }
}

ORACLE_CASE("DistPointOrientedBox.compute.2d")
{
    int32_t mode = Mode(io);
    auto point = Point<2>(io, mode);
    Vector2<double> center{}, extent{};
    std::array<Vector2<double>, 2> axis{};
    OrientedBoxInputs2(io, mode, center, axis, extent);
    OrientedBox2<double> box(center, axis, extent);
    DCPQuery<double, Vector2<double>, OrientedBox2<double>> query{};
    auto r = query(point, box);
    EmitTwoClosest3(io, r);
}

ORACLE_CASE("DistPointOrientedBox.compute.3d")
{
    int32_t mode = Mode(io);
    auto point = Point<3>(io, mode);
    Vector3<double> center{}, extent{};
    std::array<Vector3<double>, 3> axis{};
    OrientedBoxInputs3(io, mode, center, axis, extent);
    OrientedBox3<double> box(center, axis, extent);
    DCPQuery<double, Vector3<double>, OrientedBox3<double>> query{};
    auto r = query(point, box);
    EmitTwoClosest3(io, r);
}

// =====================================================================
// DistPoint2Arc2.h, DistLine2Arc2.h
// =====================================================================
namespace
{
    // The twelve lattice points of the circle of radius 5 (the Pythagorean
    // triple 3-4-5), so that the lattice modes put the arc endpoints exactly
    // on the circle and Arc2::Contains is evaluated with exact arithmetic.
    Vector2<double> const kCirclePoints[12] =
    {
        { 5.0, 0.0 }, { 4.0, 3.0 }, { 3.0, 4.0 }, { 0.0, 5.0 },
        { -3.0, 4.0 }, { -4.0, 3.0 }, { -5.0, 0.0 }, { -4.0, -3.0 },
        { -3.0, -4.0 }, { 0.0, -5.0 }, { 3.0, -4.0 }, { 4.0, -3.0 }
    };

    // An arc whose endpoints lie on its circle. In the lattice modes the
    // endpoints are exact; otherwise they are built from unrecorded angle
    // draws and only the resulting points are recorded, so no libm function
    // appears on the replayed path.
    void ArcInputs(oracle::Ctx& io, int32_t mode, Vector2<double>& center,
        double& radius, std::array<Vector2<double>, 2>& end)
    {
        double const twoPi = 6.283185307179586476925286766559;
        if (Lattice(mode))
        {
            center = io.latticeVec<2>(-3, 3);
            radius = io.given(5.0);
            int32_t i0 = io.rawInteger(0, 11);
            int32_t i1 = io.rawInteger(0, 11);
            end[0] = io.givenVec(center + kCirclePoints[i0]);
            end[1] = io.givenVec(center + kCirclePoints[i1]);
        }
        else
        {
            center = io.vec<2>(-5.0, 5.0);
            radius = io.real(0.25, 4.0);
            double a0 = io.raw(0.0, twoPi);
            double a1 = io.raw(0.0, twoPi);
            Vector2<double> e0{ center[0] + radius * std::cos(a0),
                center[1] + radius * std::sin(a0) };
            Vector2<double> e1{ center[0] + radius * std::cos(a1),
                center[1] + radius * std::sin(a1) };
            end[0] = io.givenVec(e0);
            end[1] = io.givenVec(e1);
        }
    }
}

ORACLE_CASE("DistPoint2Arc2.compute")
{
    // Mode 4 puts the point exactly at the arc center, which is the
    // 'equidistant' branch of the point-circle query.
    int32_t mode = Mode(io);
    Vector2<double> center{};
    double radius = 0.0;
    std::array<Vector2<double>, 2> end{};
    ArcInputs(io, mode, center, radius, end);
    Vector2<double> point{};
    if (mode == 4)
    {
        point = io.givenVec(center);
    }
    else
    {
        point = Point<2>(io, mode);
    }
    Arc2<double> arc(center, radius, end[0], end[1]);
    DCPQuery<double, Vector2<double>, Arc2<double>> query{};
    auto r = query(point, arc);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
    io.outBool(r.equidistant);
}

ORACLE_CASE("DistLine2Arc2.compute")
{
    // Mode 4 sends the line through the arc center, which makes the
    // line-circle query report two closest pairs; whether they are on the arc
    // decides between the loop result and the endpoint comparison.
    int32_t mode = Mode(io);
    Vector2<double> center{};
    double radius = 0.0;
    std::array<Vector2<double>, 2> end{};
    ArcInputs(io, mode, center, radius, end);
    Vector2<double> origin{}, direction{};
    if (mode == 4)
    {
        Vector2<double> d = RawUnit<2>(io);
        double t = io.raw(-3.0, 3.0);
        origin = io.givenVec(center - t * d);
        direction = io.givenVec(d);
    }
    else
    {
        origin = Point<2>(io, mode);
        direction = Direction<2>(io, mode);
    }
    Line2<double> line(origin, direction);
    Arc2<double> arc(center, radius, end[0], end[1]);
    DCPQuery<double, Line2<double>, Arc2<double>> query{};
    auto r = query(line, arc);
    EmitCircle2Pairs(io, r);
}

// =====================================================================
// DistLine2OrientedBox2.h
// =====================================================================
ORACLE_CASE("DistLine2OrientedBox2.compute")
{
    // Mode 4 uses the zero direction, which reaches DoQuery0D; the lattice
    // modes give directions with a zero component (DoQuery1D) in the box
    // frame because the box axes are the coordinate axes there.
    int32_t mode = Mode(io);
    Vector2<double> center{}, extent{};
    std::array<Vector2<double>, 2> axis{};
    OrientedBoxInputs2(io, mode, center, axis, extent);
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
    OrientedBox2<double> box(center, axis, extent);
    DCPQuery<double, Line2<double>, OrientedBox2<double>> query{};
    auto r = query(line, box);
    EmitOneParameter<decltype(r), 2>(io, r);
}

// =====================================================================
// DistLine3AlignedBox3.h
// =====================================================================
ORACLE_CASE("DistLine3AlignedBox3.compute")
{
    int32_t mode = Mode(io);
    Vector3<double> lo{}, hi{};
    AlignedBoxInputs<3>(io, mode, lo, hi);
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
    AlignedBox3<double> box(lo, hi);
    DCPQuery<double, Line3<double>, AlignedBox3<double>> query{};
    auto r = query(line, box);
    EmitOneParameter<decltype(r), 3>(io, r);
}

// =====================================================================
// DistLine3OrientedBox3.h
// =====================================================================
namespace
{
    void Line3OrientedBox3Inputs(oracle::Ctx& io, int32_t mode,
        Vector3<double>& center, std::array<Vector3<double>, 3>& axis,
        Vector3<double>& extent, Vector3<double>& origin,
        Vector3<double>& direction)
    {
        OrientedBoxInputs3(io, mode, center, axis, extent);
        origin = Point<3>(io, mode);
        if (mode == 4)
        {
            direction = io.givenVec(Vector3<double>{ 0.0, 0.0, 0.0 });
        }
        else
        {
            direction = Direction<3>(io, mode);
        }
    }
}

ORACLE_CASE("DistLine3OrientedBox3.compute")
{
    // closest[0] is NOT emitted here: upstream corrupts it (see the
    // '.deviation' case below). Everything else on this path is shared with
    // the canonical-box query and is unaffected by that defect.
    int32_t mode = Mode(io);
    Vector3<double> center{}, extent{}, origin{}, direction{};
    std::array<Vector3<double>, 3> axis{};
    Line3OrientedBox3Inputs(io, mode, center, axis, extent, origin, direction);
    Line3<double> line(origin, direction);
    OrientedBox3<double> box(center, axis, extent);
    DCPQuery<double, Line3<double>, OrientedBox3<double>> query{};
    auto r = query(line, box);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter);
    io.outVec(r.closest[1]);
}

ORACLE_CASE("DistLine3OrientedBox3.compute.deviation")
{
    // Upstream assigns the world-space line point to result.closest[0] and
    // then transforms it again as if it were in the box frame, so closest[0]
    // is box.center + sum_j (world point)[j] * box.axis[j]. The port drops the
    // premature assignment, which makes closest[0] the line point.
    int32_t mode = Mode(io);
    Vector3<double> center{}, extent{}, origin{}, direction{};
    std::array<Vector3<double>, 3> axis{};
    Line3OrientedBox3Inputs(io, mode, center, axis, extent, origin, direction);
    Line3<double> line(origin, direction);
    OrientedBox3<double> box(center, axis, extent);
    DCPQuery<double, Line3<double>, OrientedBox3<double>> query{};
    auto r = query(line, box);
    io.outVec(r.closest[0]);
}

// =====================================================================
// DistLine3Rectangle3.h
// =====================================================================
ORACLE_CASE("DistLine3Rectangle3.compute")
{
    // Mode 4 aims the line at a point of the rectangle, which reaches the
    // 'line intersects the rectangle' early return; the lattice modes make
    // the line parallel to the plane of the rectangle (NdD exactly zero).
    int32_t mode = Mode(io);
    auto center = Point<3>(io, mode);
    std::array<Vector3<double>, 3> frame = Axes3(io, mode);
    std::array<Vector3<double>, 2> axis{ frame[0], frame[1] };
    auto extent = Extent<2>(io, mode);
    Vector3<double> origin{}, direction{};
    if (mode == 4)
    {
        double s0 = io.raw(-1.0, 1.0);
        double s1 = io.raw(-1.0, 1.0);
        Vector3<double> target = center + (s0 * extent[0]) * axis[0]
            + (s1 * extent[1]) * axis[1];
        Vector3<double> d = RawUnit<3>(io);
        double t = io.raw(0.5, 4.0);
        origin = io.givenVec(target - t * d);
        direction = io.givenVec(d);
    }
    else
    {
        origin = Point<3>(io, mode);
        direction = Direction<3>(io, mode);
    }
    Line3<double> line(origin, direction);
    Rectangle3<double> rectangle(center, axis, extent);
    DCPQuery<double, Line3<double>, Rectangle3<double>> query{};
    auto r = query(line, rectangle);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter);
    io.outReal(r.cartesian[0]);
    io.outReal(r.cartesian[1]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

// =====================================================================
// DistLine3Triangle3.h
// =====================================================================
ORACLE_CASE("DistLine3Triangle3.compute")
{
    // Mode 4 aims the line at an interior point of the triangle, which
    // reaches the 'line intersects the triangle' early return; the lattice
    // modes produce degenerate triangles and lines exactly in the plane of
    // the triangle (NdD exactly zero).
    int32_t mode = Mode(io);
    auto v0 = Point<3>(io, mode);
    auto v1 = Point<3>(io, mode);
    auto v2 = Point<3>(io, mode);
    Vector3<double> origin{}, direction{};
    if (mode == 4)
    {
        double b1 = io.raw(0.0, 0.5);
        double b2 = io.raw(0.0, 0.5);
        Vector3<double> target = v0 + b1 * (v1 - v0) + b2 * (v2 - v0);
        Vector3<double> d = RawUnit<3>(io);
        double t = io.raw(0.5, 4.0);
        origin = io.givenVec(target - t * d);
        direction = io.givenVec(d);
    }
    else
    {
        origin = Point<3>(io, mode);
        direction = Direction<3>(io, mode);
    }
    Line3<double> line(origin, direction);
    Triangle3<double> triangle(v0, v1, v2);
    DCPQuery<double, Line3<double>, Triangle3<double>> query{};
    auto r = query(line, triangle);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter);
    io.outReal(r.barycentric[0]);
    io.outReal(r.barycentric[1]);
    io.outReal(r.barycentric[2]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

// =====================================================================
// DistPoint3Circle3.h
// =====================================================================
namespace
{
    // A circle with a unit-length normal. The lattice modes use a coordinate
    // axis as the normal and an integer radius so that the projection
    // arithmetic is exact.
    void Circle3Inputs(oracle::Ctx& io, int32_t mode, Vector3<double>& center,
        Vector3<double>& normal, double& radius)
    {
        center = Point<3>(io, mode);
        if (Lattice(mode))
        {
            int32_t k = io.rawInteger(0, 2);
            int32_t s = io.rawInteger(0, 1);
            Vector3<double> n{ 0.0, 0.0, 0.0 };
            n[k] = (s == 0 ? -1.0 : 1.0);
            normal = io.givenVec(n);
            radius = static_cast<double>(io.integer(1, 3));
        }
        else
        {
            normal = io.unit<3>();
            radius = io.real(0.25, 4.0);
        }
    }
}

ORACLE_CASE("DistPoint3Circle3.compute")
{
    // Mode 4 places the point exactly on the normal line through the circle
    // center, which is the 'equidistant' branch. The normal is a coordinate
    // axis there so that the scaled projection is exactly the zero vector.
    int32_t mode = Mode(io);
    Vector3<double> center{}, normal{};
    double radius = 0.0;
    if (mode == 4)
    {
        center = Point<3>(io, mode);
        int32_t k = io.rawInteger(0, 2);
        Vector3<double> n{ 0.0, 0.0, 0.0 };
        n[k] = 1.0;
        normal = io.givenVec(n);
        radius = io.real(0.25, 4.0);
    }
    else
    {
        Circle3Inputs(io, mode, center, normal, radius);
    }
    Vector3<double> point{};
    if (mode == 4)
    {
        double t = io.raw(-4.0, 4.0);
        point = io.givenVec(center + t * normal);
    }
    else
    {
        point = Point<3>(io, mode);
    }
    Circle3<double> circle(center, normal, radius);
    DCPQuery<double, Vector3<double>, Circle3<double>> query{};
    auto r = query(point, circle);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
    io.outBool(r.equidistant);
}

// =====================================================================
// DistPoint3Cylinder3.h
// =====================================================================
namespace
{
    void Cylinder3Inputs(oracle::Ctx& io, int32_t mode,
        Vector3<double>& axisOrigin, Vector3<double>& axisDirection,
        double& radius, double& height)
    {
        axisOrigin = Point<3>(io, mode);
        if (Lattice(mode))
        {
            int32_t k = io.rawInteger(0, 2);
            int32_t s = io.rawInteger(0, 1);
            Vector3<double> d{ 0.0, 0.0, 0.0 };
            d[k] = (s == 0 ? -1.0 : 1.0);
            axisDirection = io.givenVec(d);
            radius = static_cast<double>(io.integer(1, 3));
            height = static_cast<double>(io.integer(1, 4));
        }
        else
        {
            axisDirection = io.unit<3>();
            radius = io.real(0.25, 3.0);
            height = io.real(0.5, 5.0);
        }
    }
}

ORACLE_CASE("DistPoint3Cylinder3.compute")
{
    // A finite cylinder. Mode 4 places the point exactly on the cylinder
    // axis, which is the 'inside' branch of the infinite-cylinder query.
    int32_t mode = Mode(io);
    Vector3<double> axisOrigin{}, axisDirection{};
    double radius = 0.0, height = 0.0;
    Cylinder3Inputs(io, mode, axisOrigin, axisDirection, radius, height);
    Vector3<double> point{};
    if (mode == 4)
    {
        double t = io.raw(-4.0, 4.0);
        point = io.givenVec(axisOrigin + t * axisDirection);
    }
    else
    {
        point = Point<3>(io, mode);
    }
    Line3<double> axis(axisOrigin, axisDirection);
    Cylinder3<double> cylinder(axis, radius, height);
    DCPQuery<double, Vector3<double>, Cylinder3<double>> query{};
    auto r = query(point, cylinder);
    EmitTwoClosest3(io, r);
}

ORACLE_CASE("DistPoint3Cylinder3.compute.deviation")
{
    // Cylinder3::MakeInfiniteCylinder sets height = -1, but the query detects
    // an infinite cylinder with 'height == std::numeric_limits<T>::max()'. The
    // sentinel therefore falls into the finite branch and trips
    // LogAssert(height > 0), so upstream throws on every record here; the port
    // uses Cylinder3::IsInfinite and answers the infinite-cylinder query.
    int32_t mode = Mode(io);
    Vector3<double> axisOrigin{}, axisDirection{};
    double radius = 0.0, height = 0.0;
    Cylinder3Inputs(io, mode, axisOrigin, axisDirection, radius, height);
    double infinite = io.given(-1.0);
    auto point = Point<3>(io, mode);
    Line3<double> axis(axisOrigin, axisDirection);
    Cylinder3<double> cylinder(axis, radius, infinite);
    DCPQuery<double, Vector3<double>, Cylinder3<double>> query{};
    auto r = query(point, cylinder);
    EmitTwoClosest3(io, r);
}

// =====================================================================
// DistPoint3Frustum3.h
// =====================================================================
namespace
{
    struct FrustumData
    {
        Vector3<double> origin{}, rVector{}, uVector{}, dVector{};
        double dMin{}, dMax{}, uBound{}, rBound{};
        Vector3<double> point{};
    };

    // True when the query reaches one of the two far-edge assignments whose
    // free coordinate is not bounded by the surrounding case analysis, and
    // the coordinate actually exceeds the far half-extent. Upstream copies
    // the test coordinate straight through there, so the reported closest
    // point is off the end of the edge and outside the frustum; the port
    // clamps it. The other eight far-edge assignments are already bounded by
    // their branch conditions, so the clamp is a no-op for them.
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
            // Every far-edge assignment in this group is guarded.
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

    // Draw the frustum and the test point from unrecorded draws. The point is
    // built in frustum coordinates so that every Voronoi region is reachable.
    FrustumData RawFrustum(oracle::Ctx& io)
    {
        FrustumData f{};
        f.origin = RawVec<3>(io, -4.0, 4.0);
        f.rVector = RawUnit<3>(io);
        Vector3<double> u{};
        double len = 0.0;
        do
        {
            u = RawVec<3>(io, -1.0, 1.0);
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
        io.givenVec(f.rVector);
        io.givenVec(f.uVector);
        io.givenVec(f.dVector);
        io.given(f.dMin);
        io.given(f.dMax);
        io.given(f.uBound);
        io.given(f.rBound);
        io.givenVec(f.point);
        Frustum3<double> frustum(f.origin, f.dVector, f.uVector, f.rVector,
            f.dMin, f.dMax, f.uBound, f.rBound);
        DCPQuery<double, Vector3<double>, Frustum3<double>> query{};
        auto r = query(f.point, frustum);
        EmitTwoClosest3(io, r);
    }
}

ORACLE_CASE("DistPoint3Frustum3.compute")
{
    // Keep only configurations on which upstream's far-edge assignments are
    // bounded by their branch conditions, where the port's clamp is a no-op.
    FrustumData f{};
    for (int32_t attempt = 0; attempt < 20000; ++attempt)
    {
        f = RawFrustum(io);
        if (!FrustumClampChanges(f)) { break; }
    }
    EmitFrustum(io, f);
}

ORACLE_CASE("DistPoint3Frustum3.compute.deviation")
{
    // Keep only configurations that reach one of the two unbounded far-edge
    // assignments with a coordinate past the end of the edge.
    FrustumData f{};
    for (int32_t attempt = 0; attempt < 20000; ++attempt)
    {
        f = RawFrustum(io);
        if (FrustumClampChanges(f)) { break; }
    }
    EmitFrustum(io, f);
}

// =====================================================================
// DistRay2Circle2.h, DistSegment2Circle2.h
// =====================================================================
namespace
{
    // A circle with an integer radius and lattice center in the lattice
    // modes, which produces exact tangency and exact intersection parameters.
    void Circle2Inputs(oracle::Ctx& io, int32_t mode, Vector2<double>& center,
        double& radius)
    {
        center = Point<2>(io, mode);
        radius = Lattice(mode) ? static_cast<double>(io.integer(1, 3))
            : io.real(0.25, 4.0);
    }
}

ORACLE_CASE("DistRay2Circle2.compute")
{
    // Mode 4 puts the ray origin strictly inside the circle, which reaches
    // the 't0 < 0 <= t1' branch that drops the first intersection point.
    int32_t mode = Mode(io);
    Vector2<double> center{};
    double radius = 0.0;
    Circle2Inputs(io, mode, center, radius);
    Vector2<double> origin{}, direction{};
    if (mode == 4)
    {
        Vector2<double> inside = RawUnit<2>(io);
        double t = io.raw(0.0, 0.75);
        origin = io.givenVec(center + (t * radius) * inside);
        direction = io.unit<2>();
    }
    else
    {
        origin = Point<2>(io, mode);
        direction = Direction<2>(io, mode);
    }
    Ray2<double> ray(origin, direction);
    Circle2<double> circle(center, radius);
    DCPQuery<double, Ray2<double>, Circle2<double>> query{};
    auto r = query(ray, circle);
    EmitCircle2Pairs(io, r);
}

ORACLE_CASE("DistSegment2Circle2.compute")
{
    // Mode 4 puts both endpoints strictly inside the circle, which reaches
    // the 't0 < 0 and t1 > 1' branch that resets the whole result.
    int32_t mode = Mode(io);
    Vector2<double> center{};
    double radius = 0.0;
    Circle2Inputs(io, mode, center, radius);
    Vector2<double> p0{}, p1{};
    if (mode == 4)
    {
        Vector2<double> d0 = RawUnit<2>(io);
        double t0 = io.raw(0.0, 0.6);
        p0 = io.givenVec(center + (t0 * radius) * d0);
        Vector2<double> d1 = RawUnit<2>(io);
        double t1 = io.raw(0.0, 0.6);
        p1 = io.givenVec(center + (t1 * radius) * d1);
    }
    else
    {
        p0 = Point<2>(io, mode);
        p1 = Point<2>(io, mode);
    }
    Segment2<double> segment(p0, p1);
    Circle2<double> circle(center, radius);
    DCPQuery<double, Segment2<double>, Circle2<double>> query{};
    auto r = query(segment, circle);
    EmitCircle2Pairs(io, r);
}

// =====================================================================
// DistRay2Triangle2.h, DistSegment2Triangle2.h
// =====================================================================
namespace
{
    void Triangle2Inputs(oracle::Ctx& io, int32_t mode, Vector2<double>& v0,
        Vector2<double>& v1, Vector2<double>& v2)
    {
        v0 = Point<2>(io, mode);
        v1 = Point<2>(io, mode);
        if (mode == 4)
        {
            // A degenerate (collinear) triangle.
            double s = io.raw(-2.0, 2.0);
            v2 = io.givenVec(v0 + s * (v1 - v0));
        }
        else
        {
            v2 = Point<2>(io, mode);
        }
    }
}

ORACLE_CASE("DistRay2Triangle2.compute")
{
    // Mode 4 aims the ray away from the triangle, which makes the line-
    // triangle parameter negative and reaches the point-triangle branch.
    int32_t mode = Mode(io);
    Vector2<double> v0{}, v1{}, v2{};
    Triangle2Inputs(io, mode, v0, v1, v2);
    Vector2<double> origin{}, direction{};
    if (mode == 4)
    {
        Vector2<double> away = RawUnit<2>(io);
        double t = io.raw(1.0, 5.0);
        origin = io.givenVec(v0 + t * away);
        direction = io.givenVec(away);
    }
    else
    {
        origin = Point<2>(io, mode);
        direction = Direction<2>(io, mode);
    }
    Ray2<double> ray(origin, direction);
    Triangle2<double> triangle(v0, v1, v2);
    DCPQuery<double, Ray2<double>, Triangle2<double>> query{};
    auto r = query(ray, triangle);
    EmitTriangle2(io, r);
}

ORACLE_CASE("DistSegment2Triangle2.compute")
{
    // Mode 4 puts the whole segment past the triangle along its own
    // direction, which reaches the 'parameter > 1' endpoint branch.
    int32_t mode = Mode(io);
    Vector2<double> v0{}, v1{}, v2{};
    Triangle2Inputs(io, mode, v0, v1, v2);
    Vector2<double> p0{}, p1{};
    if (mode == 4)
    {
        Vector2<double> away = RawUnit<2>(io);
        double t = io.raw(1.0, 5.0);
        p0 = io.givenVec(v0 + t * away);
        double s = io.raw(0.25, 2.0);
        p1 = io.givenVec(p0 + s * away);
    }
    else
    {
        p0 = Point<2>(io, mode);
        p1 = Point<2>(io, mode);
    }
    Segment2<double> segment(p0, p1);
    Triangle2<double> triangle(v0, v1, v2);
    DCPQuery<double, Segment2<double>, Triangle2<double>> query{};
    auto r = query(segment, triangle);
    EmitTriangle2(io, r);
}

// =====================================================================
// DistRay3CanonicalBox3.h, DistSegment3CanonicalBox3.h
// =====================================================================
ORACLE_CASE("DistRay3CanonicalBox3.compute")
{
    // Mode 4 aims the ray away from the box, which makes the line parameter
    // negative and reaches the point-box branch.
    int32_t mode = Mode(io);
    auto extent = Extent<3>(io, mode);
    Vector3<double> origin{}, direction{};
    if (mode == 4)
    {
        Vector3<double> away = RawUnit<3>(io);
        double t = io.raw(4.0, 10.0);
        origin = io.givenVec(t * away);
        direction = io.givenVec(away);
    }
    else
    {
        origin = Point<3>(io, mode);
        direction = Direction<3>(io, mode);
    }
    Ray3<double> ray(origin, direction);
    CanonicalBox3<double> box(extent);
    DCPQuery<double, Ray3<double>, CanonicalBox3<double>> query{};
    auto r = query(ray, box);
    EmitOneParameter<decltype(r), 3>(io, r);
}

ORACLE_CASE("DistSegment3CanonicalBox3.compute")
{
    // Mode 4 puts the whole segment past the box along its own direction,
    // which reaches the 'parameter > 1' endpoint branch.
    int32_t mode = Mode(io);
    auto extent = Extent<3>(io, mode);
    Vector3<double> p0{}, p1{};
    if (mode == 4)
    {
        Vector3<double> away = RawUnit<3>(io);
        double t = io.raw(4.0, 10.0);
        p0 = io.givenVec(t * away);
        double s = io.raw(0.25, 2.0);
        p1 = io.givenVec(p0 + s * away);
    }
    else
    {
        p0 = Point<3>(io, mode);
        p1 = Point<3>(io, mode);
    }
    Segment3<double> segment(p0, p1);
    CanonicalBox3<double> box(extent);
    DCPQuery<double, Segment3<double>, CanonicalBox3<double>> query{};
    auto r = query(segment, box);
    EmitOneParameter<decltype(r), 3>(io, r);
}

// =====================================================================
// DistLine3Circle3.h
// =====================================================================
namespace
{
    void EmitLine3Circle3(oracle::Ctx& io,
        DCPQuery<double, Line3<double>, Circle3<double>>::Result const& r)
    {
        io.outInt(r.numClosestPairs);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outBool(r.equidistant);
        for (size_t j = 0; j < r.numClosestPairs; ++j)
        {
            io.outVec(r.linearClosest[j]);
            io.outVec(r.circularClosest[j]);
        }
    }
}

namespace
{
    // Finalize projects each critical line point onto the plane of the circle
    // and normalizes the projection. Upstream does not check that the
    // projection is nonzero; Vector::Normalize leaves a zero vector at zero,
    // so the reported circle point becomes the circle center, which is not on
    // the circle, and the reported distance is the distance to the center
    // (UPSTREAM-FINDINGS DistLine3Circle3.h item 3, issue #421). PDFSection421
    // reaches it exactly whenever a critical parameter makes the in-plane
    // component of the line point vanish, which the integer lattice produces
    // readily. This replicates the condition so the main case can keep only
    // the records upstream gets right.
    //
    // A zero projection means the critical line point is on the axis of the
    // circle, so every circle point is equidistant from it and the port's
    // corrected distance sqrt(h^2 + r^2) always exceeds upstream's |h|. The
    // port therefore keeps upstream's choice of critical point whenever
    // upstream did not choose the defective one, and the two disagree exactly
    // when upstream reports the circle center as a closest circle point.
    // That symptom is the selector, so it is exact and needs no replica of
    // the branch analysis.
    bool ReportsCircleCenter(
        DCPQuery<double, Line3<double>, Circle3<double>>::Result const& r,
        Circle3<double> const& circle)
    {
        for (size_t j = 0; j < r.numClosestPairs; ++j)
        {
            if (r.circularClosest[j] == circle.center) { return true; }
        }
        return false;
    }

    // The three closed-form branches. The circle normal is a coordinate axis
    // and the line direction and the offset of the line origin are integer
    // multiples of coordinate axes, so Cross(N, M) and Cross(N, D) are
    // exactly zero where the branch conditions require it:
    //   mode 0,1: M parallel to N and the origin off the axis  -> PDFSection412
    //   mode 2,3: the origin on the axis, M not parallel to N  -> PDFSection421
    //   mode 4:   M parallel to N and the origin on the axis   -> PDFSection411
    // No libm call other than sqrt appears on any of these paths.
    void Line3Circle3Axis(oracle::Ctx& io, bool wantZeroProjection)
    {
        // The zero projection is reachable only from PDFSection421, which the
        // lattice mode produces, so the deviation case uses that mode only.
        int32_t mode = wantZeroProjection ? 2 : Mode(io);
        DCPQuery<double, Line3<double>, Circle3<double>> query{};
        Vector3<double> center{}, normal{}, origin{}, direction{};
        double radius = 0.0;
        for (int32_t attempt = 0; attempt < 4000; ++attempt)
        {
            center = RawLatticeVec<3>(io, -3, 3);
            int32_t k = io.rawInteger(0, 2);
            int32_t sign = io.rawInteger(0, 1);
            Vector3<double> n{ 0.0, 0.0, 0.0 };
            n[k] = (sign == 0 ? -1.0 : 1.0);
            normal = n;
            radius = static_cast<double>(io.rawInteger(1, 3));
            if (mode == 2 || mode == 3)
            {
                double t = static_cast<double>(io.rawInteger(-4, 4));
                origin = center + t * normal;
                direction = RawLatticeDir<3>(io, -2, 2);
            }
            else
            {
                double m = static_cast<double>(io.rawInteger(1, 3));
                double s = (io.rawInteger(0, 1) == 0 ? -m : m);
                Vector3<double> offAxis = RawLatticeVec<3>(io, -3, 3);
                offAxis[k] = 0.0;
                if (mode == 4) { offAxis.MakeZero(); }
                double t = static_cast<double>(io.rawInteger(-4, 4));
                origin = center + t * normal + offAxis;
                direction = s * normal;
            }
            Circle3<double> candidateCircle(center, normal, radius);
            auto probe = query(Line3<double>(origin, direction),
                candidateCircle);
            if (ReportsCircleCenter(probe, candidateCircle)
                == wantZeroProjection)
            {
                break;
            }
        }

        io.givenVec(center);
        io.givenVec(normal);
        io.given(radius);
        io.givenVec(origin);
        io.givenVec(direction);
        Line3<double> line(origin, direction);
        Circle3<double> circle(center, normal, radius);
        auto r = query(line, circle);
        EmitLine3Circle3(io, r);
    }
}

ORACLE_CASE("DistLine3Circle3.compute.axis") { Line3Circle3Axis(io, false); }
ORACLE_CASE("DistLine3Circle3.compute.axis.deviation")
{
    Line3Circle3Axis(io, true);
}

namespace
{
    // tauHat solves G'(tau) = 1 and is used by PDFSection422 both as a
    // bisection bracket endpoint and, through 'intercept', as the branch
    // selector. Upstream omits the division by a2 (UPSTREAM-FINDINGS, issue
    // #247), so the port and upstream disagree grossly wherever tauHat can
    // influence the answer. That is exactly when a1 > sqrt(a3) (tauHat is
    // formed at all) and |a0| does not exceed both intercepts: when |a0| is
    // above both, upstream and the port take the same single-critical-point
    // branch with the same bracket, which involves no tauHat.
    bool TauHatMatters(Line3<double> const& line, Circle3<double> const& circle)
    {
        Vector3<double> const& N = circle.normal;
        Vector3<double> const& M = line.direction;
        Vector3<double> D = line.origin - circle.center;
        Vector3<double> NxM = Cross(N, M);
        Vector3<double> NxD = Cross(N, D);
        Vector3<double> const vzero{};
        if (NxM == vzero || NxD == vzero)
        {
            return false;
        }
        double NxMdNxM = Dot(NxM, NxM);
        if (!(NxMdNxM > 0.0))
        {
            return false;
        }
        double s = -Dot(NxM, NxD) / NxMdNxM;
        Vector3<double> E = s * M + D;
        double MdM = Dot(M, M);
        Vector3<double> NxE = Cross(N, E);
        double a0 = Dot(M, E) / MdM;
        double a1 = circle.radius * NxMdNxM / MdM;
        double a2 = NxMdNxM;
        double a3 = Dot(NxE, NxE);
        if (a3 == 0.0 || !(a1 > std::sqrt(a3)))
        {
            return false;
        }
        double const twoThirds = 2.0 / 3.0;
        double core = std::fabs(std::pow(a1 * a3, twoThirds) - a3);
        double tauHatUpstream = std::sqrt(core);
        double tauHatFixed = std::sqrt(core / a2);
        auto G = [&](double t)
        {
            return a1 * t / std::sqrt(a2 * t * t + a3);
        };
        double interceptUpstream = G(tauHatUpstream) - tauHatUpstream;
        double interceptFixed = G(tauHatFixed) - tauHatFixed;
        return std::fabs(a0)
            <= std::max(interceptUpstream, interceptFixed);
    }

    // An independent reference for the line-circle distance: the minimum over
    // the circle of the point-to-line distance. The objective is smooth and
    // bounded in the circle angle, so a dense scan plus a golden-section
    // refinement resolves it to near machine precision however badly the line
    // parameterization is conditioned. It is used only to decide which
    // records a case keeps.
    double MinDistanceOverCircle(Line3<double> const& line,
        Circle3<double> const& circle)
    {
        double const twoPi = 6.283185307179586476925286766559;
        Vector3<double> const& n = circle.normal;
        Vector3<double> seed = (std::fabs(n[0]) < 0.9
            ? Vector3<double>{ 1.0, 0.0, 0.0 }
            : Vector3<double>{ 0.0, 1.0, 0.0 });
        Vector3<double> u = seed - Dot(seed, n) * n;
        Normalize(u);
        Vector3<double> w = Cross(n, u);
        double MdM = Dot(line.direction, line.direction);
        auto f = [&](double a)
        {
            Vector3<double> q = circle.center
                + (circle.radius * std::cos(a)) * u
                + (circle.radius * std::sin(a)) * w;
            Vector3<double> diff = q - line.origin;
            double t = Dot(diff, line.direction) / MdM;
            return Length(diff - t * line.direction);
        };

        int32_t const samples = 2048;
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

    // True when upstream's own answer is correct: the reported distance
    // agrees with the independent reference to 1e-9 relative.
    bool UpstreamIsSound(Line3<double> const& line, Circle3<double> const& circle,
        double distance)
    {
        double reference = MinDistanceOverCircle(line, circle);
        double scale = std::fmax(1.0, std::fmax(std::fabs(distance),
            std::fabs(reference)));
        return std::fabs(distance - reference) <= 1e-9 * scale;
    }

    enum class Line3Circle3Family
    {
        // Upstream is sound and tauHat cannot influence the answer.
        sound,
        // tauHat influences the branch choice or a bisection bracket.
        tauHat,
        // A line nearly perpendicular to the plane of the circle, where
        // upstream's back-substitution t = tau + s cancels away the
        // significant digits of the critical parameter.
        nearPerpendicular
    };

    void Line3Circle3Generic(oracle::Ctx& io, Line3Circle3Family family)
    {
        DCPQuery<double, Line3<double>, Circle3<double>> query{};
        Vector3<double> center{}, normal{}, origin{}, direction{};
        double radius = 0.0;
        // The last candidate of the family that upstream answered without
        // throwing, used when no attempt satisfies the acceptance test. It
        // keeps a record of the intended family in the golden file instead of
        // a C++ exception.
        Vector3<double> fbCenter{}, fbNormal{}, fbOrigin{}, fbDirection{};
        double fbRadius = 0.0;
        bool haveFallback = false, accepted = false;
        for (int32_t attempt = 0; attempt < 4000; ++attempt)
        {
            center = RawVec<3>(io, -4.0, 4.0);
            normal = RawUnit<3>(io);
            radius = io.raw(0.25, 4.0);
            Vector3<double> perp = RawUnit<3>(io);
            perp = perp - Dot(perp, normal) * normal;
            Normalize(perp);
            if (family == Line3Circle3Family::nearPerpendicular)
            {
                // The direction is the circle normal perturbed in plane by
                // eps, drawn log-uniformly, and scaled; the origin sits at a
                // height on the axis plus an in-plane offset. |NxM| ~ eps, so
                // s ~ 1/eps^2 while the bisection bracket has width
                // r*|NxM|/Dot(M,M) ~ eps.
                double eps = std::pow(10.0, io.raw(-8.5, -6.0));
                double scale = std::pow(2.0, io.raw(-2.0, 2.0));
                direction = scale * (normal + eps * perp);
                double height = io.raw(-8.0, 8.0);
                double offset = io.raw(0.25, 4.0);
                Vector3<double> perp2 = RawUnit<3>(io);
                perp2 = perp2 - Dot(perp2, normal) * normal;
                Normalize(perp2);
                origin = center + height * normal + offset * perp2;
            }
            else if (family == Line3Circle3Family::tauHat)
            {
                // Aim the line near the axis of the circle, where the squared
                // distance can have three critical points. The direction is
                // long, so a2 = |NxM|^2 can exceed 1, which is the regime in
                // which upstream's missing division by a2 misplaces the
                // bisection bracket. That is how the segment-circle query
                // calls this file (M = P1 - P0).
                double dirScale = std::pow(2.0, io.raw(1.0, 4.0));
                Vector3<double> dirUnit = RawUnit<3>(io);
                direction = dirScale * dirUnit;
                double h = io.raw(-4.0, 4.0);
                Vector3<double> onAxis = center + h * normal;
                double off = io.raw(0.0, 0.25) * radius;
                double t = io.raw(0.5, 4.0);
                origin = onAxis + off * perp - t * direction;
            }
            else
            {
                direction = RawUnit<3>(io);
                origin = RawVec<3>(io, -4.0, 4.0);
            }

            Line3<double> candidateLine(origin, direction);
            Circle3<double> candidateCircle(center, normal, radius);
            bool matters = false;
            bool sound = false;
            try
            {
                matters = TauHatMatters(candidateLine, candidateCircle);
                auto probe = query(candidateLine, candidateCircle);
                sound = UpstreamIsSound(candidateLine, candidateCircle,
                    probe.distance);
            }
            catch (std::exception const&)
            {
                // A collapsed bisection bracket; upstream throws. Never kept.
                continue;
            }

            fbCenter = center;
            fbNormal = normal;
            fbRadius = radius;
            fbOrigin = origin;
            fbDirection = direction;
            haveFallback = true;

            if (family == Line3Circle3Family::sound)
            {
                accepted = (!matters && sound);
            }
            else if (family == Line3Circle3Family::tauHat)
            {
                // Keep only the records on which tauHat can influence the
                // answer AND upstream's reported distance is actually wrong.
                accepted = (matters && !sound);
            }
            else
            {
                accepted = !sound;
            }
            if (accepted) { break; }
        }

        if (!accepted && haveFallback)
        {
            center = fbCenter;
            normal = fbNormal;
            radius = fbRadius;
            origin = fbOrigin;
            direction = fbDirection;
        }

        io.givenVec(center);
        io.givenVec(normal);
        io.given(radius);
        io.givenVec(origin);
        io.givenVec(direction);
        Line3<double> line(origin, direction);
        Circle3<double> circle(center, normal, radius);
        auto r = query(line, circle);
        EmitLine3Circle3(io, r);
    }
}

ORACLE_CASE("DistLine3Circle3.compute")
{
    Line3Circle3Generic(io, Line3Circle3Family::sound);
}

ORACLE_CASE("DistLine3Circle3.compute.tauHat")
{
    Line3Circle3Generic(io, Line3Circle3Family::tauHat);
}

ORACLE_CASE("DistLine3Circle3.compute.nearPerpendicular")
{
    Line3Circle3Generic(io, Line3Circle3Family::nearPerpendicular);
}
