// Verify group 21 (distance): differential cases for the Dist* headers that
// group 21 covers. See ORACLE.md.
//
// Generator modes are selected from io.index() so that every case mixes
// uniform random inputs, small-lattice inputs (exact integer arithmetic, which
// produces parallel, touching, tangent and coincident configurations) and
// constructions aimed at a specific branch. Every generator records exactly
// the same number of doubles in every mode, so the TypeScript replay reads the
// inputs without knowing which mode produced them.
#define ORACLE_FAMILY "v21-distance"
#include "Oracle.h"

#include <Mathematics/DistPlane3CanonicalBox3.h>
#include <Mathematics/DistPoint2Parallelogram2.h>
#include <Mathematics/DistPoint3ConvexPolyhedron3.h>
#include <Mathematics/DistPointHyperellipsoid.h>
#include <Mathematics/DistPointHyperplane.h>
#include <Mathematics/DistRay2AlignedBox2.h>
#include <Mathematics/DistRay2Arc2.h>
#include <Mathematics/DistRay2OrientedBox2.h>
#include <Mathematics/DistRay3AlignedBox3.h>
#include <Mathematics/DistRay3Circle3.h>
#include <Mathematics/DistRay3OrientedBox3.h>
#include <Mathematics/DistRay3Rectangle3.h>
#include <Mathematics/DistRay3Triangle3.h>
#include <Mathematics/DistSegment2AlignedBox2.h>
#include <Mathematics/DistSegment2Arc2.h>
#include <Mathematics/DistSegment2OrientedBox2.h>
#include <Mathematics/DistSegment3AlignedBox3.h>
#include <Mathematics/DistSegment3Circle3.h>
#include <Mathematics/DistSegment3OrientedBox3.h>
#include <Mathematics/DistSegment3Rectangle3.h>
#include <Mathematics/DistSegment3Triangle3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <utility>
#include <vector>

using namespace gte;

namespace
{
    // ---- generators ------------------------------------------------------
    //
    // mode 0,1: uniform random
    // mode 2,3: small lattice (integers), so degenerate configurations occur
    // mode 4:   constructed (degenerate primitives, aimed configurations, ...)
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
    // the box headers and is generated in the lattice modes.
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

    // An orthonormal pair of 2D box axes. The lattice modes use the coordinate
    // axes so that the box frame is the world frame and lattice directions
    // keep their exactly zero components (the DoQuery1D and DoQuery0D dispatch
    // arms of the underlying line-box query).
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

    // The result shape shared by the line/ray/segment versus box queries.
    template <typename Result, int32_t N>
    void EmitOneParameter(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    // The oriented-box queries in 3D inherit the closest[0] defect of
    // DistLine3OrientedBox3 (UPSTREAM-FINDINGS, issue #421): the world-space
    // line point is written into closest[0] before the loop that maps the
    // box-frame closest points back to the world, so it is transformed twice.
    // distance, sqrDistance, parameter and closest[1] come from the canonical
    // box query and are unaffected, so the ordinary cases emit those and the
    // deviation cases emit closest[0] alone.
    template <typename Result>
    void EmitNoLinearClosest(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outVec(r.closest[1]);
    }

    template <typename Result, int32_t N>
    void EmitTwoClosest(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }
}

// =====================================================================
// DistRay2AlignedBox2.h, DistRay3AlignedBox3.h
// DistSegment2AlignedBox2.h, DistSegment3AlignedBox3.h
// =====================================================================
namespace
{
    // A ray aimed at 'target', recording the origin and the direction (2N
    // doubles, as the unconstrained modes do). Half of the records place the
    // target behind the ray origin, which is the point-primitive clamp.
    template <int32_t N>
    void AimedRayN(oracle::Ctx& io, Vector<N, double> const& target,
        Vector<N, double>& origin, Vector<N, double>& direction)
    {
        Vector<N, double> unit = RawUnit<N>(io);
        double t = io.raw(-3.0, 3.0);
        origin = io.givenVec(target - t * unit);
        direction = io.givenVec(unit);
    }

    // A segment through 'target', recording both endpoints (2N doubles). The
    // draws reach both the t < 0 and the t > 1 endpoint clamps.
    template <int32_t N>
    void AimedSegmentN(oracle::Ctx& io, Vector<N, double> const& target,
        Vector<N, double>& p0, Vector<N, double>& p1)
    {
        Vector<N, double> unit = RawUnit<N>(io);
        double t = io.raw(-3.0, 3.0);
        double len = io.raw(0.25, 4.0);
        Vector<N, double> start = target - t * unit;
        p0 = io.givenVec(start);
        p1 = io.givenVec(start + len * unit);
    }

    template <int32_t N>
    Vector<N, double> RawAlignedBoxPoint(oracle::Ctx& io,
        Vector<N, double> const& lo, Vector<N, double> const& hi)
    {
        Vector<N, double> target{};
        for (int32_t i = 0; i < N; ++i)
        {
            target[i] = io.raw(lo[i], hi[i]);
        }
        return target;
    }

    template <int32_t N>
    void RayAlignedBox(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        Vector<N, double> lo{}, hi{}, origin{}, direction{};
        AlignedBoxInputs<N>(io, mode, lo, hi);
        if (mode == 4)
        {
            AimedRayN<N>(io, RawAlignedBoxPoint<N>(io, lo, hi), origin,
                direction);
        }
        else
        {
            origin = Point<N>(io, mode);
            direction = Direction<N>(io, mode);
        }
        Ray<N, double> ray(origin, direction);
        AlignedBox<N, double> box(lo, hi);
        DCPQuery<double, Ray<N, double>, AlignedBox<N, double>> query{};
        auto r = query(ray, box);
        EmitOneParameter<decltype(r), N>(io, r);
    }

    template <int32_t N>
    void SegmentAlignedBox(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        Vector<N, double> lo{}, hi{}, p0{}, p1{};
        AlignedBoxInputs<N>(io, mode, lo, hi);
        if (mode == 4)
        {
            AimedSegmentN<N>(io, RawAlignedBoxPoint<N>(io, lo, hi), p0, p1);
        }
        else
        {
            p0 = Point<N>(io, mode);
            p1 = Point<N>(io, mode);
        }
        Segment<N, double> segment(p0, p1);
        AlignedBox<N, double> box(lo, hi);
        DCPQuery<double, Segment<N, double>, AlignedBox<N, double>> query{};
        auto r = query(segment, box);
        EmitOneParameter<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistRay2AlignedBox2.compute") { RayAlignedBox<2>(io); }
ORACLE_CASE("DistRay3AlignedBox3.compute") { RayAlignedBox<3>(io); }
ORACLE_CASE("DistSegment2AlignedBox2.compute") { SegmentAlignedBox<2>(io); }
ORACLE_CASE("DistSegment3AlignedBox3.compute") { SegmentAlignedBox<3>(io); }

// =====================================================================
// DistRay2OrientedBox2.h, DistRay3OrientedBox3.h
// DistSegment2OrientedBox2.h, DistSegment3OrientedBox3.h
// =====================================================================
namespace
{
    // A point of the box, in world coordinates, from unrecorded draws.
    template <int32_t N>
    Vector<N, double> RawBoxPoint(oracle::Ctx& io,
        Vector<N, double> const& center,
        std::array<Vector<N, double>, N> const& axis,
        Vector<N, double> const& extent)
    {
        Vector<N, double> p = center;
        for (int32_t i = 0; i < N; ++i)
        {
            p = p + io.raw(-extent[i], extent[i]) * axis[i];
        }
        return p;
    }

    // The recorded inputs of an oriented box, in the order
    // center, axis[0..N-1], extent.
    std::array<Vector2<double>, 2> BoxInputs2(oracle::Ctx& io, int32_t mode,
        Vector2<double>& center, Vector2<double>& extent)
    {
        center = Point<2>(io, mode);
        auto axis = Axes2(io, mode);
        extent = Extent<2>(io, mode);
        return axis;
    }

    std::array<Vector3<double>, 3> BoxInputs3(oracle::Ctx& io, int32_t mode,
        Vector3<double>& center, Vector3<double>& extent)
    {
        center = Point<3>(io, mode);
        auto axis = Axes3(io, mode);
        extent = Extent<3>(io, mode);
        return axis;
    }

    // A ray origin and direction: aimed at the box in mode 4 (half of those
    // records place the box behind the ray origin, which is the point-box
    // clamp), otherwise unconstrained.
    template <int32_t N>
    void RayInputsForBox(oracle::Ctx& io, int32_t mode,
        Vector<N, double> const& center,
        std::array<Vector<N, double>, N> const& axis,
        Vector<N, double> const& extent,
        Vector<N, double>& origin, Vector<N, double>& direction)
    {
        if (mode == 4)
        {
            AimedRayN<N>(io, RawBoxPoint<N>(io, center, axis, extent), origin,
                direction);
        }
        else
        {
            origin = Point<N>(io, mode);
            direction = Direction<N>(io, mode);
        }
    }

    // Segment endpoints: in mode 4 the segment starts before the box and is
    // long enough to reach past it about half the time, so both the t < 0 and
    // the t > 1 endpoint clamps are exercised.
    template <int32_t N>
    void SegmentInputsForBox(oracle::Ctx& io, int32_t mode,
        Vector<N, double> const& center,
        std::array<Vector<N, double>, N> const& axis,
        Vector<N, double> const& extent,
        Vector<N, double>& p0, Vector<N, double>& p1)
    {
        if (mode == 4)
        {
            AimedSegmentN<N>(io, RawBoxPoint<N>(io, center, axis, extent), p0,
                p1);
        }
        else
        {
            p0 = Point<N>(io, mode);
            p1 = Point<N>(io, mode);
        }
    }
}

ORACLE_CASE("DistRay2OrientedBox2.compute")
{
    int32_t mode = Mode(io);
    Vector2<double> center{}, extent{}, origin{}, direction{};
    auto axis = BoxInputs2(io, mode, center, extent);
    RayInputsForBox<2>(io, mode, center, axis, extent, origin, direction);
    Ray2<double> ray(origin, direction);
    OrientedBox2<double> box(center, axis, extent);
    DCPQuery<double, Ray2<double>, OrientedBox2<double>> query{};
    auto r = query(ray, box);
    EmitOneParameter<decltype(r), 2>(io, r);
}

ORACLE_CASE("DistSegment2OrientedBox2.compute")
{
    int32_t mode = Mode(io);
    Vector2<double> center{}, extent{}, p0{}, p1{};
    auto axis = BoxInputs2(io, mode, center, extent);
    SegmentInputsForBox<2>(io, mode, center, axis, extent, p0, p1);
    Segment2<double> segment(p0, p1);
    OrientedBox2<double> box(center, axis, extent);
    DCPQuery<double, Segment2<double>, OrientedBox2<double>> query{};
    auto r = query(segment, box);
    EmitOneParameter<decltype(r), 2>(io, r);
}

ORACLE_CASE("DistRay3OrientedBox3.compute")
{
    // closest[0] is not emitted; see EmitNoLinearClosest and the
    // '.deviation' case below.
    int32_t mode = Mode(io);
    Vector3<double> center{}, extent{}, origin{}, direction{};
    auto axis = BoxInputs3(io, mode, center, extent);
    RayInputsForBox<3>(io, mode, center, axis, extent, origin, direction);
    Ray3<double> ray(origin, direction);
    OrientedBox3<double> box(center, axis, extent);
    DCPQuery<double, Ray3<double>, OrientedBox3<double>> query{};
    auto r = query(ray, box);
    EmitNoLinearClosest(io, r);
}

ORACLE_CASE("DistRay3OrientedBox3.compute.deviation")
{
    // Same generator as the ordinary case, restricted to the branch that
    // returns the line-box output (parameter >= 0), which is where upstream's
    // doubly-transformed closest[0] appears. When the ray parameter is clamped
    // to 0 the query writes ray.origin into closest[0] itself and the port
    // agrees, so those records are rejected: they would make the deviation
    // look narrower than it is.
    DCPQuery<double, Ray3<double>, OrientedBox3<double>> query{};
    Vector3<double> center{}, extent{}, origin{}, direction{};
    std::array<Vector3<double>, 3> axis{};
    bool found = false;
    for (int32_t attempt = 0; attempt < 4000 && !found; ++attempt)
    {
        center = RawVec<3>(io, -5.0, 5.0);
        axis[0] = RawUnit<3>(io);
        double len = 0.0;
        Vector3<double> a1{};
        do
        {
            a1 = RawVec<3>(io, -1.0, 1.0);
            a1 = a1 - Dot(a1, axis[0]) * axis[0];
            len = Length(a1);
        } while (len < 0.25);
        Normalize(a1);
        axis[1] = a1;
        axis[2] = Cross(axis[0], axis[1]);
        extent = RawVec<3>(io, 0.25, 3.0);
        origin = RawVec<3>(io, -5.0, 5.0);
        direction = RawUnit<3>(io);
        auto probe = query(Ray3<double>(origin, direction),
            OrientedBox3<double>(center, axis, extent));
        found = (probe.parameter > 0.0);
    }

    io.givenVec(center);
    io.givenVec(axis[0]);
    io.givenVec(axis[1]);
    io.givenVec(axis[2]);
    io.givenVec(extent);
    io.givenVec(origin);
    io.givenVec(direction);
    Ray3<double> ray(origin, direction);
    OrientedBox3<double> box(center, axis, extent);
    auto r = query(ray, box);
    io.outVec(r.closest[0]);
}

ORACLE_CASE("DistSegment3OrientedBox3.compute")
{
    // closest[0] is not emitted; see EmitNoLinearClosest and the
    // '.deviation' case below.
    int32_t mode = Mode(io);
    Vector3<double> center{}, extent{}, p0{}, p1{};
    auto axis = BoxInputs3(io, mode, center, extent);
    SegmentInputsForBox<3>(io, mode, center, axis, extent, p0, p1);
    Segment3<double> segment(p0, p1);
    OrientedBox3<double> box(center, axis, extent);
    DCPQuery<double, Segment3<double>, OrientedBox3<double>> query{};
    auto r = query(segment, box);
    EmitNoLinearClosest(io, r);
}

ORACLE_CASE("DistSegment3OrientedBox3.compute.deviation")
{
    // As for the ray: keep only the records that return the line-box output,
    // which is where upstream's doubly-transformed closest[0] appears.
    DCPQuery<double, Segment3<double>, OrientedBox3<double>> query{};
    Vector3<double> center{}, extent{}, p0{}, p1{};
    std::array<Vector3<double>, 3> axis{};
    bool found = false;
    for (int32_t attempt = 0; attempt < 4000 && !found; ++attempt)
    {
        center = RawVec<3>(io, -5.0, 5.0);
        axis[0] = RawUnit<3>(io);
        double len = 0.0;
        Vector3<double> a1{};
        do
        {
            a1 = RawVec<3>(io, -1.0, 1.0);
            a1 = a1 - Dot(a1, axis[0]) * axis[0];
            len = Length(a1);
        } while (len < 0.25);
        Normalize(a1);
        axis[1] = a1;
        axis[2] = Cross(axis[0], axis[1]);
        extent = RawVec<3>(io, 0.25, 3.0);
        Vector3<double> target = RawBoxPoint<3>(io, center, axis, extent);
        Vector3<double> unit = RawUnit<3>(io);
        double t = io.raw(0.25, 3.0);
        double seg = io.raw(1.5, 6.0);
        p0 = target - t * unit;
        p1 = p0 + seg * unit;
        auto probe = query(Segment3<double>(p0, p1),
            OrientedBox3<double>(center, axis, extent));
        found = (probe.parameter > 0.0 && probe.parameter < 1.0);
    }

    io.givenVec(center);
    io.givenVec(axis[0]);
    io.givenVec(axis[1]);
    io.givenVec(axis[2]);
    io.givenVec(extent);
    io.givenVec(p0);
    io.givenVec(p1);
    Segment3<double> segment(p0, p1);
    OrientedBox3<double> box(center, axis, extent);
    auto r = query(segment, box);
    io.outVec(r.closest[0]);
}

// =====================================================================
// DistRay3Rectangle3.h, DistSegment3Rectangle3.h
// DistRay3Triangle3.h, DistSegment3Triangle3.h
// =====================================================================
namespace
{
    // A rectangle: center, two orthonormal 3D axes, two positive extents. The
    // lattice modes use two coordinate axes and integer extents. Axes3 records
    // three axes; the third is recorded too (and ignored) so that the input
    // layout does not depend on the mode.
    std::array<Vector3<double>, 2> RectangleInputs(oracle::Ctx& io,
        int32_t mode, Vector3<double>& center, Vector2<double>& extent)
    {
        center = Point<3>(io, mode);
        auto frame = Axes3(io, mode);
        std::array<Vector3<double>, 2> axis{ frame[0], frame[1] };
        if (Lattice(mode))
        {
            extent = io.latticeVec<2>(1, 3);
        }
        else
        {
            extent = io.vec<2>(0.25, 3.0);
        }
        return axis;
    }

    Vector3<double> RawRectanglePoint(oracle::Ctx& io,
        Vector3<double> const& center,
        std::array<Vector3<double>, 2> const& axis,
        Vector2<double> const& extent)
    {
        double c0 = io.raw(-extent[0], extent[0]);
        double c1 = io.raw(-extent[1], extent[1]);
        return center + c0 * axis[0] + c1 * axis[1];
    }

    Vector3<double> RawTrianglePoint(oracle::Ctx& io,
        std::array<Vector3<double>, 3> const& v)
    {
        double b0 = io.raw(0.0, 1.0);
        double b1 = io.raw(0.0, 1.0 - b0);
        double b2 = 1.0 - b0 - b1;
        return b0 * v[0] + b1 * v[1] + b2 * v[2];
    }

    template <typename Result>
    void EmitCartesian(oracle::Ctx& io, Result const& r)
    {
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outReal(r.parameter);
        io.outReal(r.cartesian[0]);
        io.outReal(r.cartesian[1]);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }

    template <typename Result>
    void EmitBarycentric(oracle::Ctx& io, Result const& r)
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

ORACLE_CASE("DistRay3Rectangle3.compute")
{
    int32_t mode = Mode(io);
    Vector3<double> center{}, origin{}, direction{};
    Vector2<double> extent{};
    auto axis = RectangleInputs(io, mode, center, extent);
    if (mode == 4)
    {
        Vector3<double> target = RawRectanglePoint(io, center, axis, extent);
        AimedRayN<3>(io, target, origin, direction);
    }
    else
    {
        origin = Point<3>(io, mode);
        direction = Direction<3>(io, mode);
    }
    Ray3<double> ray(origin, direction);
    Rectangle3<double> rectangle(center, axis, extent);
    DCPQuery<double, Ray3<double>, Rectangle3<double>> query{};
    auto r = query(ray, rectangle);
    EmitCartesian(io, r);
}

ORACLE_CASE("DistSegment3Rectangle3.compute")
{
    int32_t mode = Mode(io);
    Vector3<double> center{}, p0{}, p1{};
    Vector2<double> extent{};
    auto axis = RectangleInputs(io, mode, center, extent);
    if (mode == 4)
    {
        Vector3<double> target = RawRectanglePoint(io, center, axis, extent);
        AimedSegmentN<3>(io, target, p0, p1);
    }
    else
    {
        p0 = Point<3>(io, mode);
        p1 = Point<3>(io, mode);
    }
    Segment3<double> segment(p0, p1);
    Rectangle3<double> rectangle(center, axis, extent);
    DCPQuery<double, Segment3<double>, Rectangle3<double>> query{};
    auto r = query(segment, rectangle);
    EmitCartesian(io, r);
}

ORACLE_CASE("DistRay3Triangle3.compute")
{
    int32_t mode = Mode(io);
    std::array<Vector3<double>, 3> v{};
    v[0] = Point<3>(io, mode);
    v[1] = Point<3>(io, mode);
    v[2] = Point<3>(io, mode);
    Vector3<double> origin{}, direction{};
    if (mode == 4)
    {
        Vector3<double> target = RawTrianglePoint(io, v);
        AimedRayN<3>(io, target, origin, direction);
    }
    else
    {
        origin = Point<3>(io, mode);
        direction = Direction<3>(io, mode);
    }
    Ray3<double> ray(origin, direction);
    Triangle3<double> triangle(v);
    DCPQuery<double, Ray3<double>, Triangle3<double>> query{};
    auto r = query(ray, triangle);
    EmitBarycentric(io, r);
}

ORACLE_CASE("DistSegment3Triangle3.compute")
{
    int32_t mode = Mode(io);
    std::array<Vector3<double>, 3> v{};
    v[0] = Point<3>(io, mode);
    v[1] = Point<3>(io, mode);
    v[2] = Point<3>(io, mode);
    Vector3<double> p0{}, p1{};
    if (mode == 4)
    {
        Vector3<double> target = RawTrianglePoint(io, v);
        AimedSegmentN<3>(io, target, p0, p1);
    }
    else
    {
        p0 = Point<3>(io, mode);
        p1 = Point<3>(io, mode);
    }
    Segment3<double> segment(p0, p1);
    Triangle3<double> triangle(v);
    DCPQuery<double, Segment3<double>, Triangle3<double>> query{};
    auto r = query(segment, triangle);
    EmitBarycentric(io, r);
}

// =====================================================================
// DistPointHyperplane.h
// =====================================================================
namespace
{
    // A hyperplane whose normal need not be a coordinate direction. The
    // lattice modes use an integer normal scaled to unit length only when it
    // is already unit length (a signed coordinate axis), which keeps the
    // signed-distance arithmetic exact; the query itself requires a
    // unit-length normal.
    template <int32_t N>
    Vector<N, double> PlaneNormal(oracle::Ctx& io, int32_t mode)
    {
        if (Lattice(mode))
        {
            int32_t k = io.rawInteger(0, N - 1);
            int32_t sign = io.rawInteger(0, 1);
            Vector<N, double> n{};
            n.MakeZero();
            n[k] = (sign == 0 ? -1.0 : 1.0);
            return io.givenVec(n);
        }
        return io.unit<N>();
    }

    template <int32_t N>
    void PointHyperplane(oracle::Ctx& io)
    {
        // Mode 4 places the query point exactly on the hyperplane, where the
        // signed distance is the rounding of a cancelling dot product.
        int32_t mode = Mode(io);
        auto normal = PlaneNormal<N>(io, mode);
        auto planeOrigin = Point<N>(io, mode);
        Vector<N, double> point{};
        if (mode == 4)
        {
            Vector<N, double> tangent = RawVec<N>(io, -4.0, 4.0);
            tangent = tangent - Dot(tangent, normal) * normal;
            point = io.givenVec(planeOrigin + tangent);
        }
        else
        {
            point = Point<N>(io, mode);
        }
        Hyperplane<N, double> plane(normal, planeOrigin);
        DCPQuery<double, Vector<N, double>, Hyperplane<N, double>> query{};
        auto r = query(point, plane);
        io.outReal(r.distance);
        io.outReal(r.signedDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }
}

ORACLE_CASE("DistPointHyperplane.compute.2d") { PointHyperplane<2>(io); }
ORACLE_CASE("DistPointHyperplane.compute.3d") { PointHyperplane<3>(io); }
ORACLE_CASE("DistPointHyperplane.compute.4d") { PointHyperplane<4>(io); }

// =====================================================================
// DistPlane3CanonicalBox3.h
// =====================================================================
ORACLE_CASE("DistPlane3CanonicalBox3.compute")
{
    // The query reflects the plane normal into the first octant and then
    // dispatches on how many components are strictly positive:
    //   3 -> DoQuery3D, 2 -> DoQuery2D, 1 -> DoQuery1D, 0 -> DoQuery0D.
    // The uniform modes give three nonzero components. The lattice modes draw
    // an integer normal on a random support and normalize it, which leaves the
    // unsupported components exactly zero (DoQuery2D and DoQuery1D), and give
    // integer extents including zero (a flat box). Mode 4 passes the zero
    // normal, the "low-probability event" branch DoQuery0D, where the query
    // consumes plane.origin rather than the reflected origin.
    int32_t mode = Mode(io);
    Vector3<double> normal{};
    if (mode == 4)
    {
        normal = io.givenVec(Vector3<double>{ 0.0, 0.0, 0.0 });
    }
    else if (Lattice(mode))
    {
        Vector3<double> n{};
        do
        {
            int32_t support = io.rawInteger(1, 7);
            n.MakeZero();
            for (int32_t i = 0; i < 3; ++i)
            {
                if ((support & (1 << i)) != 0)
                {
                    n[i] = static_cast<double>(io.rawInteger(-3, 3));
                }
            }
        } while (Length(n) == 0.0);
        Normalize(n);
        normal = io.givenVec(n);
    }
    else
    {
        normal = io.unit<3>();
    }
    auto planeOrigin = Point<3>(io, mode);
    auto extent = Extent<3>(io, mode);
    Plane3<double> plane(normal, planeOrigin);
    CanonicalBox3<double> box(extent);
    DCPQuery<double, Plane3<double>, CanonicalBox3<double>> query{};
    auto r = query(plane, box);
    EmitTwoClosest<decltype(r), 3>(io, r);
}

// =====================================================================
// DistPoint2Parallelogram2.h
// =====================================================================
namespace
{
    // Parallelogram axes. They are not required to be orthogonal or unit
    // length, but Parallelogram2's constructor asserts DotPerp(a0, a1) > 0, so
    // the pair is swapped when it is left-handed and redrawn when it is
    // degenerate. The lattice modes make them integer vectors, which gives
    // exact edge parameters and nearly degenerate Gram matrices.
    std::array<Vector2<double>, 2> ParallelogramAxes(oracle::Ctx& io,
        int32_t mode)
    {
        std::array<Vector2<double>, 2> axis{};
        for (int32_t attempt = 0; attempt < 1000; ++attempt)
        {
            if (Lattice(mode))
            {
                axis[0] = RawLatticeVec<2>(io, -2, 2);
                axis[1] = RawLatticeVec<2>(io, -2, 2);
            }
            else
            {
                axis[0] = RawVec<2>(io, -2.0, 2.0);
                axis[1] = RawVec<2>(io, -2.0, 2.0);
            }
            double dp = DotPerp(axis[0], axis[1]);
            if (dp < 0.0)
            {
                std::swap(axis[0], axis[1]);
                dp = -dp;
            }
            if (dp > 0.0) { break; }
        }
        io.givenVec(axis[0]);
        io.givenVec(axis[1]);
        return axis;
    }
}

ORACLE_CASE("DistPoint2Parallelogram2.compute")
{
    // Mode 4 places the query point exactly on an edge of the parallelogram
    // (a corner of the [-1,1]^2 domain of GetMinimizer).
    int32_t mode = Mode(io);
    auto center = Point<2>(io, mode);
    auto axis = ParallelogramAxes(io, mode);
    Vector2<double> point{};
    if (mode == 4)
    {
        double k = static_cast<double>(io.rawInteger(-4, 4)) / 4.0;
        int32_t edge = io.rawInteger(0, 3);
        double y0 = (edge == 0 ? -1.0 : (edge == 1 ? 1.0 : k));
        double y1 = (edge == 2 ? -1.0 : (edge == 3 ? 1.0 : k));
        point = io.givenVec(center + y0 * axis[0] + y1 * axis[1]);
    }
    else
    {
        point = Point<2>(io, mode);
    }
    Parallelogram2<double> pgm(center, axis);
    DCPQuery<double, Vector2<double>, Parallelogram2<double>> query{};
    auto r = query(point, pgm);
    EmitTwoClosest<decltype(r), 2>(io, r);
}

ORACLE_CASE("DistPoint2Parallelogram2.getMinimizer")
{
    // GetMinimizer is public upstream and is exported by the port. A is the
    // Gram matrix B^T*B of the two axes, exactly as operator() builds it; Z
    // ranges over the nine regions of the [-1,1]^2 domain and lands exactly on
    // the boundary in the lattice modes (quarter-integer coordinates) and
    // exactly on a corner in mode 4.
    int32_t mode = Mode(io);
    auto axis = ParallelogramAxes(io, mode);
    Vector2<double> Z{};
    if (mode == 4)
    {
        int32_t k0 = io.rawInteger(-1, 1);
        int32_t k1 = io.rawInteger(-1, 1);
        Z = io.givenVec(Vector2<double>{
            static_cast<double>(k0), static_cast<double>(k1) });
    }
    else if (Lattice(mode))
    {
        double z0 = static_cast<double>(io.rawInteger(-8, 8)) / 4.0;
        double z1 = static_cast<double>(io.rawInteger(-8, 8)) / 4.0;
        Z = io.givenVec(Vector2<double>{ z0, z1 });
    }
    else
    {
        Z = io.vec<2>(-3.0, 3.0);
    }
    Matrix2x2<double> B{};
    B.SetCol(0, axis[0]);
    B.SetCol(1, axis[1]);
    Matrix2x2<double> A = MultiplyATB(B, B);
    DCPQuery<double, Vector2<double>, Parallelogram2<double>> query{};
    auto K = query.GetMinimizer(A, Z);
    io.outVec(K);
}

ORACLE_CASE("DistPoint2Parallelogram2.compute.degenerate")
{
    // Parallelogram2's constructor asserts DotPerp(axis[0], axis[1]) > 0. The
    // axes here are left-handed or collinear, so the C++ build throws on every
    // record and the port must throw on the same records.
    int32_t mode = Mode(io);
    auto center = Point<2>(io, mode);
    std::array<Vector2<double>, 2> axis{};
    for (int32_t attempt = 0; attempt < 1000; ++attempt)
    {
        axis[0] = RawLatticeVec<2>(io, -2, 2);
        if (Lattice(mode))
        {
            // Collinear: DotPerp is exactly zero.
            double s = static_cast<double>(io.rawInteger(-2, 2));
            axis[1] = s * axis[0];
        }
        else
        {
            axis[1] = RawVec<2>(io, -2.0, 2.0);
        }
        if (DotPerp(axis[0], axis[1]) <= 0.0) { break; }
    }
    io.givenVec(axis[0]);
    io.givenVec(axis[1]);
    auto point = Point<2>(io, mode);
    Parallelogram2<double> pgm(center, axis);
    DCPQuery<double, Vector2<double>, Parallelogram2<double>> query{};
    auto r = query(point, pgm);
    EmitTwoClosest<decltype(r), 2>(io, r);
}

// =====================================================================
// DistPointHyperellipsoid.h
// =====================================================================
namespace
{
    // Extents must be positive (the algorithm divides by the smallest one).
    // The lattice modes use integer extents, which makes equal extents (and
    // therefore the tie-breaking of the std::pair sort) and exact
    // sumZSqr == 1 reachable.
    template <int32_t N>
    Vector<N, double> EllipsoidExtent(oracle::Ctx& io, int32_t mode)
    {
        if (Lattice(mode))
        {
            return io.latticeVec<N>(1, 3);
        }
        return io.vec<N>(0.25, 3.0);
    }

    // The query point in hyperellipsoid coordinates, from unrecorded draws;
    // the caller records the world-space point. Mode 4 zeroes one coordinate
    // exactly, which is the y[i] == 0 branch of SqrDistanceSpecial (and, when
    // it is the smallest-extent coordinate, the subhyperellipsoid branch), and
    // places half of the records exactly on the hyperellipsoid, which is the
    // sumZSqr == 1 branch of Bisector. The lattice modes reach that branch
    // exactly whenever the point is an axis vertex.
    template <int32_t N>
    Vector<N, double> RawEllipsoidLocalPoint(oracle::Ctx& io, int32_t mode,
        Vector<N, double> const& extent)
    {
        if (mode == 4)
        {
            Vector<N, double> y{};
            for (int32_t i = 0; i < N; ++i)
            {
                y[i] = io.raw(-4.0, 4.0);
            }
            int32_t k = io.rawInteger(0, N - 1);
            y[k] = 0.0;
            if (io.rawInteger(0, 1) == 0)
            {
                double sum = 0.0;
                for (int32_t i = 0; i < N; ++i)
                {
                    double z = y[i] / extent[i];
                    sum += z * z;
                }
                if (sum > 0.0)
                {
                    y = y / std::sqrt(sum);
                    y[k] = 0.0;
                }
            }
            return y;
        }
        if (Lattice(mode))
        {
            return RawLatticeVec<N>(io, -3, 3);
        }
        return RawVec<N>(io, -4.0, 4.0);
    }

    // The query point is built in hyperellipsoid coordinates and recorded in
    // world coordinates, so the branch-reaching coordinates survive the change
    // of basis exactly in the lattice modes (where the axes are the coordinate
    // axes).
    template <int32_t N>
    void PointHyperellipsoid(oracle::Ctx& io,
        std::array<Vector<N, double>, N> const& axis,
        Vector<N, double> const& center, int32_t mode)
    {
        auto extent = EllipsoidExtent<N>(io, mode);
        auto local = RawEllipsoidLocalPoint<N>(io, mode, extent);
        Vector<N, double> point = center;
        for (int32_t i = 0; i < N; ++i) { point = point + local[i] * axis[i]; }
        io.givenVec(point);
        Hyperellipsoid<N, double> hyperellipsoid(center, axis, extent);
        DCPQuery<double, Vector<N, double>, Hyperellipsoid<N, double>> query{};
        auto r = query(point, hyperellipsoid);
        EmitTwoClosest<decltype(r), N>(io, r);
    }

    // The second public overload: the hyperellipsoid is axis-aligned and
    // centered at the origin, so only the extents are used.
    template <int32_t N>
    void PointHyperellipsoidAxisAligned(oracle::Ctx& io)
    {
        int32_t mode = Mode(io);
        auto extent = EllipsoidExtent<N>(io, mode);
        auto point = io.givenVec(RawEllipsoidLocalPoint<N>(io, mode, extent));
        DCPQuery<double, Vector<N, double>, Hyperellipsoid<N, double>> query{};
        auto r = query(point, extent);
        EmitTwoClosest<decltype(r), N>(io, r);
    }
}

ORACLE_CASE("DistPointHyperellipsoid.compute.2d")
{
    int32_t mode = Mode(io);
    auto center = Point<2>(io, mode);
    auto axis = Axes2(io, mode);
    PointHyperellipsoid<2>(io, axis, center, mode);
}

ORACLE_CASE("DistPointHyperellipsoid.compute.3d")
{
    int32_t mode = Mode(io);
    auto center = Point<3>(io, mode);
    auto axis = Axes3(io, mode);
    PointHyperellipsoid<3>(io, axis, center, mode);
}

ORACLE_CASE("DistPointHyperellipsoid.computeAxisAligned.2d")
{
    PointHyperellipsoidAxisAligned<2>(io);
}

ORACLE_CASE("DistPointHyperellipsoid.computeAxisAligned.3d")
{
    PointHyperellipsoidAxisAligned<3>(io);
}

// =====================================================================
// DistPoint3ConvexPolyhedron3.h
// =====================================================================
namespace
{
    // The face triangles of a tetrahedron whose vertices are ordered so that
    // the signed volume Dot(Cross(V1-V0, V2-V0), V3-V0) is positive; the
    // windings are then counterclockwise seen from outside.
    std::vector<int32_t> const kTetraIndices{
        0, 2, 1,
        0, 3, 2,
        0, 1, 3,
        1, 2, 3 };

    // The face triangles of a box whose vertex i has the sign pattern
    // (bit0 -> x, bit1 -> y, bit2 -> z), with a right-handed axis frame.
    std::vector<int32_t> const kBoxIndices{
        0, 4, 6,  0, 6, 2,
        1, 3, 7,  1, 7, 5,
        0, 1, 5,  0, 5, 4,
        2, 7, 3,  2, 6, 7,
        0, 2, 3,  0, 3, 1,
        4, 5, 7,  4, 7, 6 };

    // The face topology is a constant of the case, but it is recorded so that
    // the replay reads it from the golden record rather than duplicating it.
    void RecordIndices(oracle::Ctx& io, std::vector<int32_t> const& indices)
    {
        io.given(static_cast<double>(indices.size()));
        for (int32_t i : indices) { io.given(static_cast<double>(i)); }
    }

    // A tetrahedron with positive signed volume. The lattice mode uses integer
    // vertices, which makes coplanar (degenerate) candidates common; those are
    // rejected because upstream's LCP then has no meaningful solution.
    std::array<Vector3<double>, 4> RawTetrahedron(oracle::Ctx& io, bool lattice)
    {
        std::array<Vector3<double>, 4> v{};
        for (int32_t attempt = 0; attempt < 1000; ++attempt)
        {
            for (int32_t i = 0; i < 4; ++i)
            {
                v[i] = lattice ? RawLatticeVec<3>(io, -3, 3)
                    : RawVec<3>(io, -4.0, 4.0);
            }
            double volume = Dot(Cross(v[1] - v[0], v[2] - v[0]), v[3] - v[0]);
            if (volume < 0.0)
            {
                std::swap(v[1], v[2]);
                volume = -volume;
            }
            if (volume > 0.5) { break; }
        }
        return v;
    }

    void EmitConvexPolyhedron(oracle::Ctx& io,
        DCPQuery<double, Vector3<double>, ConvexPolyhedron3<double>>::Result const& r)
    {
        io.outBool(r.queryIsSuccessful);
        io.outInt(r.numLCPIterations);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        io.outVec(r.closest[0]);
        io.outVec(r.closest[1]);
    }
}

ORACLE_CASE("DistPoint3ConvexPolyhedron3.compute.tetrahedron")
{
    // The per-query LCP solver (numTriangles == 0 construction). Mode 4 places
    // the query point inside the tetrahedron, where the distance is 0.
    int32_t mode = Mode(io);
    auto v = RawTetrahedron(io, Lattice(mode));
    std::vector<Vector3<double>> vertices(v.begin(), v.end());
    for (auto const& vertex : vertices) { io.givenVec(vertex); }
    RecordIndices(io, kTetraIndices);
    Vector3<double> point{};
    if (mode == 4)
    {
        double b0 = io.raw(0.0, 1.0);
        double b1 = io.raw(0.0, 1.0 - b0);
        double b2 = io.raw(0.0, 1.0 - b0 - b1);
        double b3 = 1.0 - b0 - b1 - b2;
        point = io.givenVec(b0 * v[0] + b1 * v[1] + b2 * v[2] + b3 * v[3]);
    }
    else
    {
        point = Point<3>(io, mode);
    }
    std::vector<int32_t> indices(kTetraIndices);
    ConvexPolyhedron3<double> polyhedron(std::move(vertices), std::move(indices),
        true, true);
    DCPQuery<double, Vector3<double>, ConvexPolyhedron3<double>> query{};
    auto r = query(point, polyhedron);
    EmitConvexPolyhedron(io, r);
}

ORACLE_CASE("DistPoint3ConvexPolyhedron3.compute.box")
{
    // The cached LCP solver (numTriangles > 0 construction, 12 triangles).
    // The box is built from unrecorded draws and its eight vertices are the
    // recorded input, so the replay does not have to reproduce the change of
    // basis. The lattice mode gives an axis-aligned box with integer corners.
    int32_t mode = Mode(io);
    bool lattice = Lattice(mode);
    Vector3<double> center = lattice ? RawLatticeVec<3>(io, -3, 3)
        : RawVec<3>(io, -5.0, 5.0);
    std::array<Vector3<double>, 3> axis{};
    if (lattice)
    {
        axis[0] = Vector3<double>{ 1.0, 0.0, 0.0 };
        axis[1] = Vector3<double>{ 0.0, 1.0, 0.0 };
        axis[2] = Vector3<double>{ 0.0, 0.0, 1.0 };
    }
    else
    {
        axis[0] = RawUnit<3>(io);
        double len = 0.0;
        Vector3<double> a1{};
        do
        {
            a1 = RawVec<3>(io, -1.0, 1.0);
            a1 = a1 - Dot(a1, axis[0]) * axis[0];
            len = Length(a1);
        } while (len < 0.25);
        Normalize(a1);
        axis[1] = a1;
        axis[2] = Cross(axis[0], axis[1]);
    }
    Vector3<double> extent = lattice ? RawLatticeVec<3>(io, 1, 3)
        : RawVec<3>(io, 0.25, 3.0);
    std::vector<Vector3<double>> vertices(8);
    for (int32_t i = 0; i < 8; ++i)
    {
        double s0 = ((i & 1) != 0 ? extent[0] : -extent[0]);
        double s1 = ((i & 2) != 0 ? extent[1] : -extent[1]);
        double s2 = ((i & 4) != 0 ? extent[2] : -extent[2]);
        vertices[i] = io.givenVec(
            center + s0 * axis[0] + s1 * axis[1] + s2 * axis[2]);
    }
    RecordIndices(io, kBoxIndices);
    auto point = Point<3>(io, mode);
    std::vector<int32_t> indices(kBoxIndices);
    ConvexPolyhedron3<double> polyhedron(std::move(vertices), std::move(indices),
        true, true);
    DCPQuery<double, Vector3<double>, ConvexPolyhedron3<double>> query(12);
    auto r = query(point, polyhedron);
    EmitConvexPolyhedron(io, r);
}

ORACLE_CASE("DistPoint3ConvexPolyhedron3.compute.noPlanes")
{
    // A polyhedron with no planes: the query returns the zeroed result with
    // queryIsSuccessful false. The constructor rejects fewer than 4 vertices
    // or fewer than 12 indices, so the member arrays stay empty.
    int32_t mode = Mode(io);
    auto point = Point<3>(io, mode);
    std::vector<Vector3<double>> vertices{};
    for (int32_t i = 0; i < 3; ++i)
    {
        vertices.push_back(Point<3>(io, mode));
    }
    std::vector<int32_t> indices{ 0, 1, 2 };
    RecordIndices(io, indices);
    ConvexPolyhedron3<double> polyhedron(std::move(vertices), std::move(indices),
        true, true);
    DCPQuery<double, Vector3<double>, ConvexPolyhedron3<double>> query{};
    auto r = query(point, polyhedron);
    EmitConvexPolyhedron(io, r);
}

ORACLE_CASE("DistPoint3ConvexPolyhedron3.setMaxLCPIterations")
{
    // SetMaxLCPIterations is the third public entry point. One iteration is
    // never enough for a tetrahedron, so the solver fails to converge and the
    // query reports queryIsSuccessful == false with the distance members left
    // at zero.
    int32_t mode = Mode(io);
    auto v = RawTetrahedron(io, Lattice(mode));
    std::vector<Vector3<double>> vertices(v.begin(), v.end());
    for (auto const& vertex : vertices) { io.givenVec(vertex); }
    RecordIndices(io, kTetraIndices);
    auto point = Point<3>(io, mode);
    int32_t maxIterations = io.integer(1, 3);
    std::vector<int32_t> indices(kTetraIndices);
    ConvexPolyhedron3<double> polyhedron(std::move(vertices), std::move(indices),
        true, true);
    DCPQuery<double, Vector3<double>, ConvexPolyhedron3<double>> query{};
    query.SetMaxLCPIterations(maxIterations);
    auto r = query(point, polyhedron);
    EmitConvexPolyhedron(io, r);
}

// =====================================================================
// DistRay2Arc2.h, DistSegment2Arc2.h
// =====================================================================
namespace
{
    // The 12 points of the integer lattice on the circle of radius 5.
    std::array<Vector2<double>, 12> const kCircle5{
        Vector2<double>{ 5.0, 0.0 }, Vector2<double>{ 4.0, 3.0 },
        Vector2<double>{ 3.0, 4.0 }, Vector2<double>{ 0.0, 5.0 },
        Vector2<double>{ -3.0, 4.0 }, Vector2<double>{ -4.0, 3.0 },
        Vector2<double>{ -5.0, 0.0 }, Vector2<double>{ -4.0, -3.0 },
        Vector2<double>{ -3.0, -4.0 }, Vector2<double>{ 0.0, -5.0 },
        Vector2<double>{ 3.0, -4.0 }, Vector2<double>{ 4.0, -3.0 } };

    // An arc: center, radius and two endpoints on the circle, traversed
    // counterclockwise from end[0] to end[1]. The lattice modes place both
    // endpoints exactly on the circle of radius 5 with integer coordinates, so
    // Arc2::Contains and the circle queries evaluate their tests exactly.
    void ArcInputs(oracle::Ctx& io, int32_t mode, Vector2<double>& center,
        double& radius, std::array<Vector2<double>, 2>& end)
    {
        double const twoPi = 6.283185307179586476925286766559;
        center = Point<2>(io, mode);
        if (Lattice(mode))
        {
            radius = io.given(5.0);
            int32_t i0 = io.rawInteger(0, 11);
            int32_t span = io.rawInteger(1, 11);
            end[0] = io.givenVec(center + kCircle5[i0]);
            end[1] = io.givenVec(center + kCircle5[(i0 + span) % 12]);
        }
        else
        {
            radius = io.real(0.5, 4.0);
            double a0 = io.raw(0.0, twoPi);
            double a1 = a0 + io.raw(0.1, twoPi - 0.1);
            end[0] = io.givenVec(center + radius
                * Vector2<double>{ std::cos(a0), std::sin(a0) });
            end[1] = io.givenVec(center + radius
                * Vector2<double>{ std::cos(a1), std::sin(a1) });
        }
    }

    // A point of the arc, from unrecorded draws (used to aim a ray or a
    // segment at the arc). The angles come from the recorded endpoints, so the
    // aim point is on the circle of the arc up to rounding.
    Vector2<double> RawArcPoint(oracle::Ctx& io, Vector2<double> const& center,
        double radius, std::array<Vector2<double>, 2> const& end)
    {
        Vector2<double> d0 = end[0] - center;
        Vector2<double> d1 = end[1] - center;
        double a0 = std::atan2(d0[1], d0[0]);
        double a1 = std::atan2(d1[1], d1[0]);
        double const twoPi = 6.283185307179586476925286766559;
        while (a1 <= a0) { a1 += twoPi; }
        double a = a0 + io.raw(0.0, 1.0) * (a1 - a0);
        return center + radius * Vector2<double>{ std::cos(a), std::sin(a) };
    }

    void EmitArc(oracle::Ctx& io,
        DCPQuery<double, Line2<double>, Circle2<double>>::Result const& r)
    {
        io.outInt(r.numClosestPairs);
        io.outReal(r.distance);
        io.outReal(r.sqrDistance);
        for (size_t j = 0; j < r.numClosestPairs; ++j)
        {
            io.outReal(r.parameter[j]);
            io.outVec(r.closest[j][0]);
            io.outVec(r.closest[j][1]);
        }
    }
}

ORACLE_CASE("DistRay2Arc2.compute")
{
    // Mode 4 aims the ray at a point of the arc; half of those records put the
    // arc behind the ray origin, where the endpoint/arc sort decides.
    int32_t mode = Mode(io);
    Vector2<double> center{}, origin{}, direction{};
    double radius = 0.0;
    std::array<Vector2<double>, 2> end{};
    ArcInputs(io, mode, center, radius, end);
    if (mode == 4)
    {
        Vector2<double> target = RawArcPoint(io, center, radius, end);
        Vector2<double> unit = RawUnit<2>(io);
        double t = io.raw(-3.0, 3.0);
        origin = io.givenVec(target - t * unit);
        direction = io.givenVec(unit);
    }
    else
    {
        origin = Point<2>(io, mode);
        direction = Direction<2>(io, mode);
    }
    Ray2<double> ray(origin, direction);
    Arc2<double> arc(center, radius, end[0], end[1]);
    DCPQuery<double, Ray2<double>, Arc2<double>> query{};
    auto r = query(ray, arc);
    EmitArc(io, r);
}

ORACLE_CASE("DistSegment2Arc2.compute")
{
    int32_t mode = Mode(io);
    Vector2<double> center{}, p0{}, p1{};
    double radius = 0.0;
    std::array<Vector2<double>, 2> end{};
    ArcInputs(io, mode, center, radius, end);
    if (mode == 4)
    {
        Vector2<double> target = RawArcPoint(io, center, radius, end);
        Vector2<double> unit = RawUnit<2>(io);
        double t = io.raw(-3.0, 3.0);
        double len = io.raw(0.25, 4.0);
        Vector2<double> start = target - t * unit;
        p0 = io.givenVec(start);
        p1 = io.givenVec(start + len * unit);
    }
    else
    {
        p0 = Point<2>(io, mode);
        p1 = Point<2>(io, mode);
    }
    Segment2<double> segment(p0, p1);
    Arc2<double> arc(center, radius, end[0], end[1]);
    DCPQuery<double, Segment2<double>, Arc2<double>> query{};
    auto r = query(segment, arc);
    EmitArc(io, r);
}

// =====================================================================
// DistRay3Circle3.h, DistSegment3Circle3.h
// =====================================================================
namespace
{
    using LineCircleResult =
        DCPQuery<double, Line3<double>, Circle3<double>>::Result;

    void EmitCircle3(oracle::Ctx& io, LineCircleResult const& r)
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

    // DistLine3Circle3::Finalize normalizes the in-plane component of a
    // critical line point without checking that it is nonzero, so upstream can
    // report the circle center as a closest circle point (UPSTREAM-FINDINGS,
    // issue #421); the port reports the true closest point. The ray and
    // segment queries inherit it through the critical points they keep. The
    // symptom is exact, so it is the selector the ordinary cases reject on.
    bool ReportsCircleCenter(LineCircleResult const& r,
        Circle3<double> const& circle)
    {
        for (size_t j = 0; j < r.numClosestPairs; ++j)
        {
            if (r.circularClosest[j] == circle.center) { return true; }
        }
        return false;
    }

    // tauHat solves G'(tau) = 1 in PDFSection422; upstream omits the division
    // by a2 (UPSTREAM-FINDINGS, issue #247). It can influence the answer only
    // when a1 > sqrt(a3) and |a0| does not exceed both intercepts. This is the
    // same predicate group 20 uses for DistLine3Circle3, applied to the line
    // that carries the ray or the segment.
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
        return std::fabs(a0) <= std::max(interceptUpstream, interceptFixed);
    }

    // An independent reference for the distance between a circle and a linear
    // component restricted to t in [tmin, tmax]: the minimum over the circle
    // of the point-to-clamped-line distance, resolved by a dense scan plus a
    // golden-section refinement. Used only to decide which records to keep.
    double MinDistanceOverCircle(Line3<double> const& line,
        Circle3<double> const& circle, double tmin, double tmax)
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
            if (t < tmin) { t = tmin; }
            if (t > tmax) { t = tmax; }
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

    bool UpstreamIsSound(Line3<double> const& line, Circle3<double> const& circle,
        double tmin, double tmax, double distance)
    {
        double reference = MinDistanceOverCircle(line, circle, tmin, tmax);
        double scale = std::fmax(1.0, std::fmax(std::fabs(distance),
            std::fabs(reference)));
        return std::fabs(distance - reference) <= 1e-9 * scale;
    }

    // The closed-form branches of DistLine3Circle3, reached through the ray or
    // segment wrapper. The circle normal is a coordinate axis and the
    // direction and the offset of the origin are integer multiples of
    // coordinate axes, so Cross(N, M) and Cross(N, D) are exactly zero where
    // the branch conditions require it. No libm call other than sqrt appears
    // on these paths.
    //   mode 0,1: M parallel to N and the origin off the axis  -> PDFSection412
    //   mode 2,3: the origin on the axis, M not parallel to N  -> PDFSection421
    //   mode 4:   M parallel to N and the origin on the axis   -> PDFSection411
    void Circle3AxisInputs(oracle::Ctx& io, int32_t mode,
        Vector3<double>& center, Vector3<double>& normal, double& radius,
        Vector3<double>& origin, Vector3<double>& direction)
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
            Vector3<double> d{};
            do
            {
                d = RawLatticeVec<3>(io, -2, 2);
            } while (d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.0);
            direction = d;
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
    }
}

ORACLE_CASE("DistRay3Circle3.compute.axis")
{
    // The closed-form branches. Records on which upstream reports the circle
    // center (the Finalize defect the port fixes) are rejected; the
    // '.axis.deviation' case keeps exactly those.
    DCPQuery<double, Ray3<double>, Circle3<double>> query{};
    Vector3<double> center{}, normal{}, origin{}, direction{};
    double radius = 0.0;
    int32_t mode = Mode(io);
    for (int32_t attempt = 0; attempt < 4000; ++attempt)
    {
        Circle3AxisInputs(io, mode, center, normal, radius, origin, direction);
        Circle3<double> candidate(center, normal, radius);
        auto probe = query(Ray3<double>(origin, direction), candidate);
        if (!ReportsCircleCenter(probe, candidate)) { break; }
    }

    io.givenVec(center);
    io.givenVec(normal);
    io.given(radius);
    io.givenVec(origin);
    io.givenVec(direction);
    Ray3<double> ray(origin, direction);
    Circle3<double> circle(center, normal, radius);
    auto r = query(ray, circle);
    EmitCircle3(io, r);
}

ORACLE_CASE("DistRay3Circle3.compute.axis.deviation")
{
    // Keeps only the records on which upstream reports the circle center.
    DCPQuery<double, Ray3<double>, Circle3<double>> query{};
    Vector3<double> center{}, normal{}, origin{}, direction{};
    double radius = 0.0;
    for (int32_t attempt = 0; attempt < 8000; ++attempt)
    {
        Circle3AxisInputs(io, 2, center, normal, radius, origin, direction);
        Circle3<double> candidate(center, normal, radius);
        auto probe = query(Ray3<double>(origin, direction), candidate);
        if (ReportsCircleCenter(probe, candidate)) { break; }
    }

    io.givenVec(center);
    io.givenVec(normal);
    io.given(radius);
    io.givenVec(origin);
    io.givenVec(direction);
    Ray3<double> ray(origin, direction);
    Circle3<double> circle(center, normal, radius);
    auto r = query(ray, circle);
    EmitCircle3(io, r);
}

ORACLE_CASE("DistSegment3Circle3.compute.axis")
{
    DCPQuery<double, Segment3<double>, Circle3<double>> query{};
    Vector3<double> center{}, normal{}, origin{}, direction{};
    double radius = 0.0;
    int32_t mode = Mode(io);
    for (int32_t attempt = 0; attempt < 4000; ++attempt)
    {
        Circle3AxisInputs(io, mode, center, normal, radius, origin, direction);
        Circle3<double> candidate(center, normal, radius);
        auto probe = query(Segment3<double>(origin, origin + direction),
            candidate);
        if (!ReportsCircleCenter(probe, candidate)) { break; }
    }

    io.givenVec(center);
    io.givenVec(normal);
    io.given(radius);
    Vector3<double> p0 = io.givenVec(origin);
    Vector3<double> p1 = io.givenVec(origin + direction);
    Segment3<double> segment(p0, p1);
    Circle3<double> circle(center, normal, radius);
    auto r = query(segment, circle);
    EmitCircle3(io, r);
}

ORACLE_CASE("DistSegment3Circle3.compute.axis.deviation")
{
    DCPQuery<double, Segment3<double>, Circle3<double>> query{};
    Vector3<double> center{}, normal{}, origin{}, direction{};
    double radius = 0.0;
    for (int32_t attempt = 0; attempt < 8000; ++attempt)
    {
        Circle3AxisInputs(io, 2, center, normal, radius, origin, direction);
        Circle3<double> candidate(center, normal, radius);
        auto probe = query(Segment3<double>(origin, origin + direction),
            candidate);
        if (ReportsCircleCenter(probe, candidate)) { break; }
    }

    io.givenVec(center);
    io.givenVec(normal);
    io.given(radius);
    Vector3<double> p0 = io.givenVec(origin);
    Vector3<double> p1 = io.givenVec(origin + direction);
    Segment3<double> segment(p0, p1);
    Circle3<double> circle(center, normal, radius);
    auto r = query(segment, circle);
    EmitCircle3(io, r);
}

namespace
{
    // The generic (bisection) path. The generator keeps only configurations on
    // which upstream's own reported distance agrees with the independent
    // reference and on which the tauHat defect cannot influence the answer, so
    // that the case measures the two implementations where upstream is sound.
    // The remaining difference is the documented conditioning fix of the
    // back-substitution t = tau + s (UPSTREAM-FINDINGS, issue #421 item 3),
    // which is why these two cases are compared with a tolerance rather than
    // bit for bit; see oracle/reports/v21-distance.md.
    void Circle3Generic(oracle::Ctx& io, bool isSegment)
    {
        DCPQuery<double, Ray3<double>, Circle3<double>> rayQuery{};
        DCPQuery<double, Segment3<double>, Circle3<double>> segQuery{};
        Vector3<double> center{}, normal{}, origin{}, direction{};
        double radius = 0.0;
        Vector3<double> fbCenter{}, fbNormal{}, fbOrigin{}, fbDirection{};
        double fbRadius = 0.0;
        bool haveFallback = false, accepted = false;
        for (int32_t attempt = 0; attempt < 4000; ++attempt)
        {
            center = RawVec<3>(io, -4.0, 4.0);
            normal = RawUnit<3>(io);
            radius = io.raw(0.25, 4.0);
            direction = RawUnit<3>(io);
            origin = RawVec<3>(io, -4.0, 4.0);
            // Half of the candidates are aimed near the circle so that the
            // interesting clamping branches of the wrappers are reached.
            if (io.rawInteger(0, 1) == 0)
            {
                Vector3<double> perp = RawUnit<3>(io);
                perp = perp - Dot(perp, normal) * normal;
                Normalize(perp);
                double h = io.raw(-1.0, 1.0);
                Vector3<double> target = center + radius * perp + h * normal;
                double t = io.raw(-1.0, 2.0);
                origin = target - t * direction;
            }

            Line3<double> candidateLine(origin, direction);
            Circle3<double> candidateCircle(center, normal, radius);
            double tmax = (isSegment ? 1.0
                : std::numeric_limits<double>::max());
            bool matters = false, sound = false;
            try
            {
                matters = TauHatMatters(candidateLine, candidateCircle);
                double distance = 0.0;
                if (isSegment)
                {
                    auto probe = segQuery(
                        Segment3<double>(origin, origin + direction),
                        candidateCircle);
                    distance = probe.distance;
                }
                else
                {
                    auto probe = rayQuery(Ray3<double>(origin, direction),
                        candidateCircle);
                    distance = probe.distance;
                }
                sound = UpstreamIsSound(candidateLine, candidateCircle,
                    0.0, tmax, distance);
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
            accepted = (!matters && sound);
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
        Circle3<double> circle(center, normal, radius);
        if (isSegment)
        {
            Vector3<double> p0 = io.givenVec(origin);
            Vector3<double> p1 = io.givenVec(origin + direction);
            auto r = segQuery(Segment3<double>(p0, p1), circle);
            EmitCircle3(io, r);
        }
        else
        {
            io.givenVec(origin);
            io.givenVec(direction);
            auto r = rayQuery(Ray3<double>(origin, direction), circle);
            EmitCircle3(io, r);
        }
    }
}

ORACLE_CASE("DistRay3Circle3.compute") { Circle3Generic(io, false); }
ORACLE_CASE("DistSegment3Circle3.compute") { Circle3Generic(io, true); }

namespace
{
    // numClosestPairs == 2 is the rarest output shape of DistSegment3Circle3:
    // it is the tie arm of SelectClosestPoint, reached from the branch
    // 't0 < 0 and t1 >= 1' where the two candidates are the segment endpoints.
    // The construction is symmetric about the circle center, X -> 2C - X,
    // which maps the circle to itself, so the two endpoints are exactly
    // equidistant from it; everything is on the integer lattice, so the two
    // point-circle distances come out bit-identical. The segment is short and
    // the circle is large, which puts both critical parameters well outside
    // [0,1] and makes the branch choice insensitive to the last bits of the
    // bisected critical parameters. A rejection loop keeps only the
    // configurations for which upstream reports two pairs at the two
    // endpoints, so the emitted output comes from the two point-circle
    // queries alone and is exact.
    //
    // The ray query has no such case: its two-pair results come from the
    // 't0 >= 0' pass-through of a two-pair line result, whose points are
    // produced by the bisection that the port deliberately back-substitutes
    // differently (issue #421 item 3). See oracle/reports/v21-distance.md.
    void Circle3TwoPairs(oracle::Ctx& io)
    {
        DCPQuery<double, Segment3<double>, Circle3<double>> segQuery{};
        Vector3<double> center{}, normal{}, p0{}, p1{};
        double radius = 0.0;
        for (int32_t attempt = 0; attempt < 20000; ++attempt)
        {
            center = RawLatticeVec<3>(io, -3, 3);
            int32_t k = io.rawInteger(0, 2);
            int32_t sign = io.rawInteger(0, 1);
            Vector3<double> n{ 0.0, 0.0, 0.0 };
            n[k] = (sign == 0 ? -1.0 : 1.0);
            normal = n;
            radius = static_cast<double>(io.rawInteger(4, 5));
            Vector3<double> u{};
            do
            {
                u = RawLatticeVec<3>(io, -1, 1);
                u[k] = 0.0;
            } while (u[0] == 0.0 && u[1] == 0.0 && u[2] == 0.0);
            double h = static_cast<double>(io.rawInteger(1, 2));
            Vector3<double> offset = u + h * normal;
            p0 = center + offset;
            p1 = center - offset;
            Circle3<double> candidate(center, normal, radius);
            auto probe = segQuery(Segment3<double>(p0, p1), candidate);
            if (probe.numClosestPairs == 2
                && probe.linearClosest[0] == p0
                && probe.linearClosest[1] == p1)
            {
                break;
            }
        }

        io.givenVec(center);
        io.givenVec(normal);
        io.given(radius);
        io.givenVec(p0);
        io.givenVec(p1);
        Circle3<double> circle(center, normal, radius);
        auto r = segQuery(Segment3<double>(p0, p1), circle);
        EmitCircle3(io, r);
    }
}

ORACLE_CASE("DistSegment3Circle3.compute.twoPairs") { Circle3TwoPairs(io); }
