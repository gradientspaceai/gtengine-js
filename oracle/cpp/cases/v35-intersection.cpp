// Verify group 35 (intersection): differential cases for the Intr* headers
// listed in plan/verify-groups.json group 35.
//
// Every generator mixes an exactly representable small-lattice mode with a
// uniform mode and, where a branch needs it, a constructed mode, so that the
// touching, tangent, parallel, coplanar, through-vertex and degenerate
// branches are reached. Every query in this group uses only + - * / sqrt fabs
// min max and comparisons, so every case is declared exact on the TypeScript
// side. The sin/cos/tan used to build orthonormal frames, unit directions and
// cone angles are applied to *unrecorded* draws and only the resulting frame,
// unit vector or trigonometric value is recorded as an input, so the C math
// library never enters the compared computation.
//
// Five headers of this group inherit deliberate port fixes of upstream
// defects; each has a 'Deviation' case, and the corresponding main case
// rejects, by an explicit predicate evaluated on upstream's own control flow,
// exactly the inputs on which the port deviates. See oracle/reports for the
// details.
#define ORACLE_FAMILY "v35-intersection"
#include "Oracle.h"

#include <Mathematics/IntrOrientedBox3Cone3.h>
#include <Mathematics/IntrRay3Cone3.h>
#include <Mathematics/IntrRay3Plane3.h>
#include <Mathematics/IntrSegment3Cone3.h>
#include <Mathematics/IntrSegment3Plane3.h>
#include <Mathematics/IntrSphere3Cone3.h>
#include <Mathematics/IntrTetrahedron3Tetrahedron3.h>
#include <Mathematics/IntrTriangle2Triangle2.h>
#include <Mathematics/IntrTriangle3AlignedBox3.h>
#include <Mathematics/IntrTriangle3CanonicalBox3.h>
#include <Mathematics/IntrTriangle3OrientedBox3.h>
#include <Mathematics/IntrTriangle3Triangle3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

using namespace gte;

namespace
{
    double const kPi = 3.141592653589793;

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

    // Orthonormal 3D frame. Mode 0 produces a signed permutation of the
    // coordinate axes (exact), which makes box axes parallel to the cone axis
    // and to the coordinate planes.
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
    Triangle<N, double> Tri(oracle::Ctx& io, int mode, int lat, double range)
    {
        Triangle<N, double> triangle{};
        triangle.v[0] = Pt<N>(io, mode, lat, range);
        triangle.v[1] = Pt<N>(io, mode, lat, range);
        triangle.v[2] = Pt<N>(io, mode, lat, range);
        return triangle;
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

    AlignedBox3<double> ABox3(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector3<double> center = Pt<3>(io, mode, lat, range);
        Vector3<double> extent = Extent<3>(io, mode, lat, range);
        AlignedBox3<double> box{};
        box.min = center - extent;
        box.max = center + extent;
        return box;
    }

    // A plane whose normal is unit length, as DistPointHyperplane assumes.
    // Mode 0 uses a signed coordinate axis and an integer constant, so every
    // signed distance is an exact integer for an integer point.
    Plane3<double> Pln(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector3<double> normal = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
        double constant = (mode == 0 ? io.lattice(-lat, lat) : io.real(-range, range));
        return Plane3<double>(normal, constant);
    }

    // Exposes the protected DoQuery helpers, which the port exports as public
    // free functions.
    template <typename Q>
    struct Expose : public Q
    {
        using Q::DoQuery;
    };

    void EmitLinePlaneFI(oracle::Ctx& io,
        FIQuery<double, Line3<double>, Plane3<double>>::Result const& r)
    {
        io.outBool(r.intersect);
        io.outInt(r.numIntersections);
        io.outReal(r.parameter);
        io.outVec(r.point);
    }

    // An integer point exactly on the plane 'Dot(N,X) = c' when N is a signed
    // coordinate axis and c is an integer: set the component along the axis to
    // c * N[j] and leave the other two components free. The point is built
    // from unrecorded draws; the case records the ray origin or the segment
    // endpoints derived from it.
    Vector3<double> LatticePointOnAxisPlane(oracle::Ctx& io,
        Plane3<double> const& plane)
    {
        int j = 0;
        for (int i = 0; i < 3; ++i)
        {
            if (plane.normal[i] != 0.0) { j = i; }
        }
        Vector3<double> Q{};
        for (int i = 0; i < 3; ++i)
        {
            Q[i] = static_cast<double>(io.rawInteger(-4, 4));
        }
        Q[j] = plane.constant * plane.normal[j];
        return Q;
    }
}

// ============================== IntrRay3Plane3 ===========================
//
// The TIQuery is self contained; the FIQuery derives from the line/plane
// FIQuery and adds a clamp to t >= 0 in its protected DoQuery, which the port
// exports as a module function. Nothing in this chain is a deliberate port
// deviation, so no rejection is needed.
//
// 'result.point' is assigned only when the query reports an intersection;
// otherwise it keeps the zeros of the Result constructor, which the port's
// default result also has, so every field is emitted on every record.

ORACLE_CASE("IntrRay3Plane3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln(io, mode, 4, 5.0);
    auto origin = Pt<3>(io, mode, 4, 5.0);
    auto direction = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    Ray3<double> ray(origin, direction);
    TIQuery<double, Ray3<double>, Plane3<double>> query;
    auto r = query(ray, plane);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrRay3Plane3.find")
{
    int mode = io.index() % 2;
    auto plane = Pln(io, mode, 4, 5.0);
    auto origin = Pt<3>(io, mode, 4, 5.0);
    auto direction = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, Plane3<double>> query;
    auto r = query(ray, plane);
    EmitLinePlaneFI(io, r);
}

ORACLE_CASE("IntrRay3Plane3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto plane = Pln(io, mode, 4, 5.0);
    auto origin = Pt<3>(io, mode, 4, 5.0);
    auto direction = (mode == 0 ? io.latticeDir<3>(-2, 2) : io.unit<3>());
    Expose<FIQuery<double, Ray3<double>, Plane3<double>>> query;
    FIQuery<double, Ray3<double>, Plane3<double>>::Result r{};
    query.DoQuery(origin, direction, plane, r);
    EmitLinePlaneFI(io, r);
}

// The ray origin is placed at 'Q - k * direction' for an integer point Q
// exactly on an axis plane, so the intersection parameter is exactly the
// integer k and the 'parameter < 0' clamp is evaluated at exact equality
// (k = 0). Integer directions perpendicular to the normal give the
// line-on-plane branch, where numIntersections is int32 max.
ORACLE_CASE("IntrRay3Plane3.find.throughLatticePoint")
{
    auto plane = Pln(io, 0, 4, 5.0);
    auto direction = io.latticeDir<3>(-2, 2);
    Vector3<double> Q = LatticePointOnAxisPlane(io, plane);
    double k = static_cast<double>(io.rawInteger(-2, 2));
    Vector3<double> origin = io.givenVec(Q - k * direction);
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, Plane3<double>> query;
    auto r = query(ray, plane);
    EmitLinePlaneFI(io, r);
}

ORACLE_CASE("IntrRay3Plane3.test.throughLatticePoint")
{
    auto plane = Pln(io, 0, 4, 5.0);
    auto direction = io.latticeDir<3>(-2, 2);
    Vector3<double> Q = LatticePointOnAxisPlane(io, plane);
    double k = static_cast<double>(io.rawInteger(-2, 2));
    Vector3<double> origin = io.givenVec(Q - k * direction);
    Ray3<double> ray(origin, direction);
    TIQuery<double, Ray3<double>, Plane3<double>> query;
    auto r = query(ray, plane);
    io.outBool(r.intersect);
}

// ============================ IntrSegment3Plane3 =========================
//
// The TIQuery tests the endpoint signed distances; the FIQuery reduces the
// segment to its centered form and clamps |t| <= extent in its protected
// DoQuery, which the port exports as a module function.

ORACLE_CASE("IntrSegment3Plane3.test")
{
    int mode = io.index() % 2;
    auto plane = Pln(io, mode, 4, 5.0);
    auto p0 = Pt<3>(io, mode, 4, 5.0);
    auto p1 = Pt<3>(io, mode, 4, 5.0);
    Segment3<double> segment(p0, p1);
    TIQuery<double, Segment3<double>, Plane3<double>> query;
    auto r = query(segment, plane);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSegment3Plane3.find")
{
    int mode = io.index() % 2;
    auto plane = Pln(io, mode, 4, 5.0);
    auto p0 = Pt<3>(io, mode, 4, 5.0);
    auto p1 = Pt<3>(io, mode, 4, 5.0);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Plane3<double>> query;
    auto r = query(segment, plane);
    EmitLinePlaneFI(io, r);
}

ORACLE_CASE("IntrSegment3Plane3.doQuery.fi")
{
    int mode = io.index() % 2;
    auto plane = Pln(io, mode, 4, 5.0);
    auto segOrigin = Pt<3>(io, mode, 4, 5.0);
    auto segDirection = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
    double segExtent = (mode == 0 ? io.lattice(1, 4) : io.real(0.25, 4.0));
    Expose<FIQuery<double, Segment3<double>, Plane3<double>>> query;
    FIQuery<double, Segment3<double>, Plane3<double>>::Result r{};
    query.DoQuery(segOrigin, segDirection, segExtent, plane, r);
    EmitLinePlaneFI(io, r);
}

// The segment runs along a signed coordinate axis with a power-of-two length
// L, so GetCenteredForm is exact (the reciprocal 1/L is exact and the unit
// direction is a coordinate axis) and the extent is exactly L/2. The endpoint
// p0 sits k units before a lattice point Q of the plane, so the intersection
// parameter is exactly k - L/2 and the clamp |t| <= extent is evaluated at
// exact equality for k = 0 and k = L.
ORACLE_CASE("IntrSegment3Plane3.find.throughLatticePoint")
{
    auto plane = Pln(io, 0, 4, 5.0);
    auto direction = AxisUnit<3>(io);
    double L = std::ldexp(1.0, io.rawInteger(1, 3));
    Vector3<double> Q = LatticePointOnAxisPlane(io, plane);
    double k = static_cast<double>(io.rawInteger(-1, 9));
    Vector3<double> p0 = io.givenVec(Q - k * direction);
    Vector3<double> p1 = io.givenVec(p0 + L * direction);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Plane3<double>> query;
    auto r = query(segment, plane);
    EmitLinePlaneFI(io, r);
}

ORACLE_CASE("IntrSegment3Plane3.test.throughLatticePoint")
{
    auto plane = Pln(io, 0, 4, 5.0);
    auto direction = io.latticeDir<3>(-2, 2);
    Vector3<double> Q = LatticePointOnAxisPlane(io, plane);
    double k = static_cast<double>(io.rawInteger(-1, 5));
    double m = static_cast<double>(io.rawInteger(-1, 5));
    Vector3<double> p0 = io.givenVec(Q - k * direction);
    Vector3<double> p1 = io.givenVec(Q + m * direction);
    Segment3<double> segment(p0, p1);
    TIQuery<double, Segment3<double>, Plane3<double>> query;
    auto r = query(segment, plane);
    io.outBool(r.intersect);
}

// =========================== IntrTriangle2Triangle2 ======================
//
// Both queries document that the input vertices are counterclockwise
// ordered. A zero-area triangle has no orientation, so the generator redraws
// it; that is a precondition of the query, not a narrowing of the input set.

namespace
{
    Triangle2<double> CcwTri2(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector2<double> v[3]{};
        for (int attempt = 0; attempt < 1000; ++attempt)
        {
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 2; ++j)
                {
                    v[i][j] = (mode == 0
                        ? static_cast<double>(io.rawInteger(-lat, lat))
                        : io.raw(-range, range));
                }
            }
            double area = DotPerp(v[1] - v[0], v[2] - v[0]);
            if (area > 0.0) { break; }
            if (area < 0.0) { std::swap(v[1], v[2]); break; }
        }
        Triangle2<double> t{};
        t.v[0] = io.givenVec(v[0]);
        t.v[1] = io.givenVec(v[1]);
        t.v[2] = io.givenVec(v[2]);
        return t;
    }

    // A counterclockwise triangle that shares the vertex Q with another
    // triangle, so that the separating-axis tests and the clipping operations
    // are evaluated at exact equality on at least one vertex.
    Triangle2<double> CcwTri2AtVertex(oracle::Ctx& io, Vector2<double> const& Q,
        int lat)
    {
        Vector2<double> v[3]{};
        for (int attempt = 0; attempt < 1000; ++attempt)
        {
            v[0] = Q;
            for (int i = 1; i < 3; ++i)
            {
                for (int j = 0; j < 2; ++j)
                {
                    v[i][j] = Q[j] + static_cast<double>(io.rawInteger(-lat, lat));
                }
            }
            double area = DotPerp(v[1] - v[0], v[2] - v[0]);
            if (area > 0.0) { break; }
            if (area < 0.0) { std::swap(v[1], v[2]); break; }
        }
        Triangle2<double> t{};
        t.v[0] = io.givenVec(v[0]);
        t.v[1] = io.givenVec(v[1]);
        t.v[2] = io.givenVec(v[2]);
        return t;
    }

    void EmitPolygon2(oracle::Ctx& io, std::vector<Vector2<double>> const& poly)
    {
        io.outInt(static_cast<int32_t>(poly.size()));
        for (auto const& p : poly) { io.outVec(p); }
    }
}

ORACLE_CASE("IntrTriangle2Triangle2.test")
{
    int mode = io.index() % 2;
    auto triangle0 = CcwTri2(io, mode, 3, 4.0);
    auto triangle1 = CcwTri2(io, mode, 3, 4.0);
    TIQuery<double, Triangle2<double>, Triangle2<double>> query;
    auto r = query(triangle0, triangle1);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrTriangle2Triangle2.find")
{
    int mode = io.index() % 2;
    auto triangle0 = CcwTri2(io, mode, 3, 4.0);
    auto triangle1 = CcwTri2(io, mode, 3, 4.0);
    FIQuery<double, Triangle2<double>, Triangle2<double>> query;
    auto r = query(triangle0, triangle1);
    EmitPolygon2(io, r.intersection);
}

// The two triangles share a vertex exactly, so 'WhichSide' evaluates a
// projection of exactly zero and the clipping query meets its vertex-on-line
// configurations.
ORACLE_CASE("IntrTriangle2Triangle2.test.sharedVertex")
{
    auto triangle0 = CcwTri2(io, 0, 3, 4.0);
    int which = io.rawInteger(0, 2);
    auto triangle1 = CcwTri2AtVertex(io, triangle0.v[which], 3);
    TIQuery<double, Triangle2<double>, Triangle2<double>> query;
    auto r = query(triangle0, triangle1);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrTriangle2Triangle2.find.sharedVertex")
{
    auto triangle0 = CcwTri2(io, 0, 3, 4.0);
    int which = io.rawInteger(0, 2);
    auto triangle1 = CcwTri2AtVertex(io, triangle0.v[which], 3);
    FIQuery<double, Triangle2<double>, Triangle2<double>> query;
    auto r = query(triangle0, triangle1);
    EmitPolygon2(io, r.intersection);
}

// ====================== IntrTetrahedron3Tetrahedron3 =====================
//
// The port corrects two result-corrupting defects of the edge-edge phase
// (docs/UPSTREAM-FINDINGS.md, issue #307): upstream compares |Dot(E0,E1)| of
// *unnormalised* edge vectors against the cosine cutoff 1 - epsilon, and its
// separation test is a plane-side test through an edge endpoint instead of
// projection-interval disjointness. 'TetraFixed' below is the port's
// algorithm; the main case keeps the inputs where it agrees with upstream and
// the two deviation cases keep the inputs where it does not, one per defect.
//
// The face-normal phase needs outward normals, that is, positively oriented
// tetrahedra (Dot(Cross(v1-v0,v2-v0), v3-v0) > 0); the generator swaps two
// vertices when the orientation is negative, which is the documented
// precondition rather than a narrowing.
//
// 'separating' is std::numeric_limits<size_t>::max() when an entry does not
// participate; that value is not representable as a double, and the port uses
// Number.MAX_SAFE_INTEGER for it, so both sides emit -1 for the sentinel.

namespace
{
    int32_t TetraWhichSide(Tetrahedron3<double> const& tetra,
        Vector3<double> const& P, Vector3<double> const& N)
    {
        size_t positive = 0, negative = 0;
        for (size_t i = 0; i < 4; ++i)
        {
            double t = Dot(N, tetra.v[i] - P);
            if (t > 0.0) { ++positive; }
            else if (t < 0.0) { ++negative; }
            if (positive > 0 && negative > 0) { return 0; }
        }
        return (positive > 0 ? +1 : -1);
    }

    void TetraProjection(Tetrahedron3<double> const& tetra,
        Vector3<double> const& N, double& mn, double& mx)
    {
        mn = Dot(N, tetra.v[0]);
        mx = mn;
        for (size_t i = 1; i < 4; ++i)
        {
            double d = Dot(N, tetra.v[i]);
            if (d < mn) { mn = d; }
            else if (d > mx) { mx = d; }
        }
    }

    // The port's corrected query, transcribed from
    // src/IntrTetrahedron3Tetrahedron3.ts.
    bool TetraFixed(Tetrahedron3<double> const& tetra0,
        Tetrahedron3<double> const& tetra1, double epsilon,
        std::array<size_t, 2>& separating)
    {
        size_t const invalid = std::numeric_limits<size_t>::max();
        separating = { invalid, invalid };

        for (size_t i = 0; i < 4; ++i)
        {
            auto const& faceIndices = Tetrahedron3<double>::GetFaceIndices(i);
            Vector3<double> P = tetra0.v[faceIndices[0]];
            Vector3<double> N = tetra0.ComputeFaceNormal(i);
            if (TetraWhichSide(tetra1, P, N) > 0)
            {
                separating = { i, invalid };
                return false;
            }
        }

        for (size_t i = 0; i < 4; ++i)
        {
            auto const& faceIndices = Tetrahedron3<double>::GetFaceIndices(i);
            Vector3<double> P = tetra1.v[faceIndices[0]];
            Vector3<double> N = tetra1.ComputeFaceNormal(i);
            if (TetraWhichSide(tetra0, P, N) > 0)
            {
                separating = { invalid, i };
                return false;
            }
        }

        double const cutoff = std::min(std::max(1.0 - epsilon, 0.0), 1.0);
        for (size_t i0 = 0; i0 < 6; ++i0)
        {
            auto const& edge0Indices = Tetrahedron3<double>::GetEdgeIndices(i0);
            Vector3<double> P0 = tetra0.v[edge0Indices[0]];
            Vector3<double> E0 = tetra0.v[edge0Indices[1]] - P0;
            double lengthE0 = std::sqrt(Dot(E0, E0));
            for (size_t i1 = 0; i1 < 6; ++i1)
            {
                auto const& edge1Indices = Tetrahedron3<double>::GetEdgeIndices(i1);
                Vector3<double> E1 = tetra1.v[edge1Indices[1]] - tetra1.v[edge1Indices[0]];
                double lengthE1 = std::sqrt(Dot(E1, E1));
                double denom = lengthE0 * lengthE1;
                if (denom == 0.0) { continue; }
                double cosAngle = std::fabs(Dot(E0, E1)) / denom;
                if (cosAngle < cutoff)
                {
                    Vector3<double> N = Cross(E0, E1);
                    if (Dot(N, N) == 0.0) { continue; }
                    double min0{}, max0{}, min1{}, max1{};
                    TetraProjection(tetra0, N, min0, max0);
                    TetraProjection(tetra1, N, min1, max1);
                    if (max0 <= min1 || max1 <= min0)
                    {
                        separating = { i0, i1 };
                        return false;
                    }
                }
            }
        }
        return true;
    }

    // True when upstream's unnormalised cutoff test decides every one of the
    // 36 edge pairs the same way the port's cosine test does, so that a
    // disagreement between the two queries can only come from the separation
    // test itself.
    bool TetraCutoffAgrees(Tetrahedron3<double> const& tetra0,
        Tetrahedron3<double> const& tetra1, double epsilon)
    {
        double const cutoff = std::min(std::max(1.0 - epsilon, 0.0), 1.0);
        for (size_t i0 = 0; i0 < 6; ++i0)
        {
            auto const& e0i = Tetrahedron3<double>::GetEdgeIndices(i0);
            Vector3<double> E0 = tetra0.v[e0i[1]] - tetra0.v[e0i[0]];
            double lengthE0 = std::sqrt(Dot(E0, E0));
            for (size_t i1 = 0; i1 < 6; ++i1)
            {
                auto const& e1i = Tetrahedron3<double>::GetEdgeIndices(i1);
                Vector3<double> E1 = tetra1.v[e1i[1]] - tetra1.v[e1i[0]];
                double lengthE1 = std::sqrt(Dot(E1, E1));
                double denom = lengthE0 * lengthE1;
                bool upstreamUses = std::fabs(Dot(E0, E1)) < cutoff;
                bool portUses = (denom != 0.0) &&
                    (std::fabs(Dot(E0, E1)) / denom < cutoff);
                if (upstreamUses != portUses) { return false; }
            }
        }
        return true;
    }

    void RawTetra(oracle::Ctx& io, int mode, int lat, double range,
        Vector3<double> center, std::array<Vector3<double>, 4>& v)
    {
        for (int attempt = 0; attempt < 1000; ++attempt)
        {
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    v[i][j] = center[j] + (mode == 0
                        ? static_cast<double>(io.rawInteger(-lat, lat))
                        : io.raw(-range, range));
                }
            }
            double det = Dot(Cross(v[1] - v[0], v[2] - v[0]), v[3] - v[0]);
            if (det > 0.0) { break; }
            if (det < 0.0) { std::swap(v[1], v[2]); break; }
        }
    }

    Tetrahedron3<double> EmitTetra(oracle::Ctx& io,
        std::array<Vector3<double>, 4> const& v)
    {
        Tetrahedron3<double> tetra{};
        for (int i = 0; i < 4; ++i) { tetra.v[i] = io.givenVec(v[i]); }
        return tetra;
    }

    void EmitSeparating(oracle::Ctx& io, std::array<size_t, 2> const& separating)
    {
        for (int i = 0; i < 2; ++i)
        {
            io.outInt(separating[i] == std::numeric_limits<size_t>::max()
                ? -1 : static_cast<int32_t>(separating[i]));
        }
    }

    // 'requirement' is 0 for the main case (the port's corrected answer equals
    // upstream's), 1 for the cutoff deviation (the answers differ) and 2 for
    // the separation deviation (the answers differ although the cutoff test
    // decides every edge pair the same way).
    void TetraCase(oracle::Ctx& io, int mode, int lat, double range,
        int requirement)
    {
        std::array<Vector3<double>, 4> v0{}, v1{};
        double epsilon = 0.0;
        Vector3<double> c0{ 0.0, 0.0, 0.0 };
        for (int attempt = 0; attempt < 20000; ++attempt)
        {
            epsilon = (io.rawInteger(0, 1) == 0 ? 0.0 : io.raw(0.0, 0.25));
            Vector3<double> c1{};
            for (int j = 0; j < 3; ++j)
            {
                c1[j] = (mode == 0
                    ? static_cast<double>(io.rawInteger(-lat, lat))
                    : io.raw(-range, range));
            }
            RawTetra(io, mode, lat, range, c0, v0);
            RawTetra(io, mode, lat, range, c1, v1);
            Tetrahedron3<double> p0{}, p1{};
            for (int i = 0; i < 4; ++i) { p0.v[i] = v0[i]; p1.v[i] = v1[i]; }
            TIQuery<double, Tetrahedron3<double>, Tetrahedron3<double>> probe;
            auto up = probe(p0, p1, epsilon);
            std::array<size_t, 2> fixedSeparating{};
            bool fixedIntersect = TetraFixed(p0, p1, epsilon, fixedSeparating);
            bool differs = (up.intersect != fixedIntersect)
                || (up.separating != fixedSeparating);
            if (requirement == 0) { if (!differs) { break; } }
            else if (requirement == 1) { if (differs) { break; } }
            else { if (differs && TetraCutoffAgrees(p0, p1, epsilon)) { break; } }
        }
        io.given(epsilon);
        auto tetra0 = EmitTetra(io, v0);
        auto tetra1 = EmitTetra(io, v1);
        TIQuery<double, Tetrahedron3<double>, Tetrahedron3<double>> query;
        auto r = query(tetra0, tetra1, epsilon);
        io.outBool(r.intersect);
        EmitSeparating(io, r.separating);
    }
}

ORACLE_CASE("IntrTetrahedron3Tetrahedron3.test")
{
    TetraCase(io, io.index() % 2, 3, 2.0, 0);
}

// Tetrahedra small enough that every |Dot(E0,E1)| is below the cutoff, so
// upstream runs its edge-edge phase too; this is the only regime in which
// upstream's edge-edge separations are reachable at all, and the case keeps
// the records on which the two separation tests agree.
ORACLE_CASE("IntrTetrahedron3Tetrahedron3.test.smallEdges")
{
    TetraCase(io, 1, 1, 0.25, 0);
}

// Deliberate deviation: the unnormalised cutoff test. The tetrahedra are
// large, so almost every |Dot(E0,E1)| exceeds the cutoff and upstream skips
// the edge-edge phase entirely.
ORACLE_CASE("IntrTetrahedron3Tetrahedron3.test.edgeCutoffDeviation")
{
    TetraCase(io, io.index() % 2, 4, 4.0, 1);
}

// Deliberate deviation: the plane-side separation test. The tetrahedra are
// small, so the cutoff test decides every edge pair the same way for both
// queries (checked by TetraCutoffAgrees) and the only remaining difference is
// the separation test.
ORACLE_CASE("IntrTetrahedron3Tetrahedron3.test.edgeSeparationDeviation")
{
    TetraCase(io, 1, 1, 0.25, 2);
}

// ====================== IntrTriangle3*Box3 (3 headers) ===================
//
// The canonical-box TIQuery is the separating-axis test; the aligned-box and
// oriented-box TIQueries transform the triangle and delegate to it. All three
// FIQueries clip the triangle against the six face planes with
// IntrConvexPolygonHyperplane (group 32) and report the inside polygon plus
// one outside polygon per clipping face, in face order, so no canonicalisation
// is needed.
//
// The lattice mode uses integer extents and integer triangle vertices, so
// vertices land exactly on faces, edges and corners of the box and every
// separating-axis and clipping comparison is evaluated at exact equality. For
// the oriented box the lattice axes are a signed permutation of the coordinate
// axes, which makes the change of basis exact.

namespace
{
    void EmitPolygon3(oracle::Ctx& io, std::vector<Vector3<double>> const& poly)
    {
        io.outInt(static_cast<int32_t>(poly.size()));
        for (auto const& p : poly) { io.outVec(p); }
    }

    template <typename R>
    void EmitTriangleBoxFI(oracle::Ctx& io, R const& r)
    {
        EmitPolygon3(io, r.insidePolygon);
        io.outInt(static_cast<int32_t>(r.outsidePolygons.size()));
        for (auto const& poly : r.outsidePolygons) { EmitPolygon3(io, poly); }
    }

    // A triangle whose vertices are drawn on the integer lattice around the
    // box centre, scaled by the box extents, so that a vertex on a face, an
    // edge or a corner of the box is common.
    Triangle3<double> LatticeTriangleAtBox(oracle::Ctx& io,
        Vector3<double> const& center, Vector3<double> const& extent,
        std::array<Vector3<double>, 3> const& axis)
    {
        Triangle3<double> triangle{};
        for (int i = 0; i < 3; ++i)
        {
            Vector3<double> p = center;
            for (int j = 0; j < 3; ++j)
            {
                double c = static_cast<double>(io.rawInteger(-2, 2)) * extent[j];
                p = p + c * axis[j];
            }
            triangle.v[i] = io.givenVec(p);
        }
        return triangle;
    }

    std::array<Vector3<double>, 3> CoordinateAxes()
    {
        std::array<Vector3<double>, 3> axis{};
        axis[0] = Vector3<double>{ 1.0, 0.0, 0.0 };
        axis[1] = Vector3<double>{ 0.0, 1.0, 0.0 };
        axis[2] = Vector3<double>{ 0.0, 0.0, 1.0 };
        return axis;
    }
}

ORACLE_CASE("IntrTriangle3CanonicalBox3.test")
{
    int mode = io.index() % 2;
    auto extent = Extent<3>(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 4, 5.0);
    CanonicalBox3<double> box(extent);
    TIQuery<double, Triangle3<double>, CanonicalBox3<double>> query;
    auto r = query(triangle, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrTriangle3CanonicalBox3.find")
{
    int mode = io.index() % 2;
    auto extent = Extent<3>(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 4, 5.0);
    CanonicalBox3<double> box(extent);
    FIQuery<double, Triangle3<double>, CanonicalBox3<double>> query;
    auto r = query(triangle, box);
    EmitTriangleBoxFI(io, r);
}

ORACLE_CASE("IntrTriangle3CanonicalBox3.test.onFaces")
{
    auto extent = Extent<3>(io, 0, 3, 3.0);
    Vector3<double> center{ 0.0, 0.0, 0.0 };
    auto triangle = LatticeTriangleAtBox(io, center, extent, CoordinateAxes());
    CanonicalBox3<double> box(extent);
    TIQuery<double, Triangle3<double>, CanonicalBox3<double>> query;
    auto r = query(triangle, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrTriangle3CanonicalBox3.find.onFaces")
{
    auto extent = Extent<3>(io, 0, 3, 3.0);
    Vector3<double> center{ 0.0, 0.0, 0.0 };
    auto triangle = LatticeTriangleAtBox(io, center, extent, CoordinateAxes());
    CanonicalBox3<double> box(extent);
    FIQuery<double, Triangle3<double>, CanonicalBox3<double>> query;
    auto r = query(triangle, box);
    EmitTriangleBoxFI(io, r);
}

ORACLE_CASE("IntrTriangle3AlignedBox3.test")
{
    int mode = io.index() % 2;
    auto box = ABox3(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 4, 5.0);
    TIQuery<double, Triangle3<double>, AlignedBox3<double>> query;
    auto r = query(triangle, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrTriangle3AlignedBox3.find")
{
    int mode = io.index() % 2;
    auto box = ABox3(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 4, 5.0);
    FIQuery<double, Triangle3<double>, AlignedBox3<double>> query;
    auto r = query(triangle, box);
    EmitTriangleBoxFI(io, r);
}

ORACLE_CASE("IntrTriangle3AlignedBox3.test.onFaces")
{
    auto box = ABox3(io, 0, 3, 3.0);
    Vector3<double> center = 0.5 * (box.max + box.min);
    Vector3<double> extent = 0.5 * (box.max - box.min);
    auto triangle = LatticeTriangleAtBox(io, center, extent, CoordinateAxes());
    TIQuery<double, Triangle3<double>, AlignedBox3<double>> query;
    auto r = query(triangle, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrTriangle3AlignedBox3.find.onFaces")
{
    auto box = ABox3(io, 0, 3, 3.0);
    Vector3<double> center = 0.5 * (box.max + box.min);
    Vector3<double> extent = 0.5 * (box.max - box.min);
    auto triangle = LatticeTriangleAtBox(io, center, extent, CoordinateAxes());
    FIQuery<double, Triangle3<double>, AlignedBox3<double>> query;
    auto r = query(triangle, box);
    EmitTriangleBoxFI(io, r);
}

ORACLE_CASE("IntrTriangle3OrientedBox3.test")
{
    int mode = io.index() % 2;
    auto box = OBox3(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 4, 5.0);
    TIQuery<double, Triangle3<double>, OrientedBox3<double>> query;
    auto r = query(triangle, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrTriangle3OrientedBox3.find")
{
    int mode = io.index() % 2;
    auto box = OBox3(io, mode, 3, 3.0);
    auto triangle = Tri<3>(io, mode, 4, 5.0);
    FIQuery<double, Triangle3<double>, OrientedBox3<double>> query;
    auto r = query(triangle, box);
    EmitTriangleBoxFI(io, r);
}

ORACLE_CASE("IntrTriangle3OrientedBox3.test.onFaces")
{
    auto box = OBox3(io, 0, 3, 3.0);
    std::array<Vector3<double>, 3> axis{ box.axis[0], box.axis[1], box.axis[2] };
    auto triangle = LatticeTriangleAtBox(io, box.center, box.extent, axis);
    TIQuery<double, Triangle3<double>, OrientedBox3<double>> query;
    auto r = query(triangle, box);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrTriangle3OrientedBox3.find.onFaces")
{
    auto box = OBox3(io, 0, 3, 3.0);
    std::array<Vector3<double>, 3> axis{ box.axis[0], box.axis[1], box.axis[2] };
    auto triangle = LatticeTriangleAtBox(io, box.center, box.extent, axis);
    FIQuery<double, Triangle3<double>, OrientedBox3<double>> query;
    auto r = query(triangle, box);
    EmitTriangleBoxFI(io, r);
}

// ========================= IntrTriangle3Triangle3 ========================
//
// Four entry points: the stationary and moving overloads of TIQuery and of
// FIQuery. The stationary queries are transcribed unchanged. The two moving
// overloads carry the one deliberate port fix of this header
// (docs/UPSTREAM-FINDINGS.md, issue #334): upstream decides "the triangles are
// parallel" with 'fabs(Dot(N0,N1)) < 1' on the *unnormalised* edge cross
// products, so the comparison against 1 scales with the product of the
// triangle areas; the port uses the stationary query's own criterion,
// |Cross(N0,N1)|^2 > 0. The two criteria agree exactly when
//   (|Dot(N0,N1)| < 1) == (Dot(Cross(N0,N1),Cross(N0,N1)) > 0),
// which is the predicate the main moving cases require and the deviation
// cases forbid. It is evaluated on the translated triangles, as the query
// itself does.

namespace
{
    // Two triangles that share an axis-aligned plane: every vertex of
    // triangle0 has the same k-th coordinate c, and each vertex of triangle1
    // has k-th coordinate c + d with a small integer d, so the signed
    // distances to the plane of triangle0 are exact integers and every
    // (numNegative, numPositive, numZero) row of the header's table is
    // reachable.
    void AxisPlaneTrianglePair(oracle::Ctx& io, int maxOffset,
        Triangle3<double>& triangle0, Triangle3<double>& triangle1)
    {
        int k = io.rawInteger(0, 2);
        double c = static_cast<double>(io.rawInteger(-3, 3));
        for (int i = 0; i < 3; ++i)
        {
            Vector3<double> p{};
            for (int j = 0; j < 3; ++j)
            {
                p[j] = static_cast<double>(io.rawInteger(-3, 3));
            }
            p[k] = c;
            triangle0.v[i] = io.givenVec(p);
        }
        for (int i = 0; i < 3; ++i)
        {
            Vector3<double> p{};
            for (int j = 0; j < 3; ++j)
            {
                p[j] = static_cast<double>(io.rawInteger(-3, 3));
            }
            p[k] = c + static_cast<double>(io.rawInteger(-maxOffset, maxOffset));
            triangle1.v[i] = io.givenVec(p);
        }
    }

    void EmitTriangleTriangleFI(oracle::Ctx& io,
        FIQuery<double, Triangle3<double>, Triangle3<double>>::Result const& r)
    {
        io.outBool(r.intersect);
        io.outReal(r.contactTime);
        EmitPolygon3(io, r.intersection);
    }

    // ---- a parameterised copy of the moving separating-axis loop ----
    // Upstream's helpers are private, so they are transcribed here. The only
    // difference from upstream is the 'fixedParallelTest' switch, which
    // selects the port's parallel criterion. The copy is used only to select
    // records; the compared outputs always come from the real upstream query.

    std::array<double, 2> MtProject(Triangle3<double> const& tri,
        Vector3<double> const& dir)
    {
        double t = Dot(dir, tri.v[0]);
        std::array<double, 2> e{ t, t };
        for (size_t i = 1; i < 3; ++i)
        {
            t = Dot(dir, tri.v[i]);
            if (t < e[0]) { e[0] = t; }
            else if (t > e[1]) { e[1] = t; }
        }
        return e;
    }

    bool MtOverlap(double tMax, double speed, std::array<double, 2> const& e0,
        std::array<double, 2> const& e1, double& tFirst, double& tLast)
    {
        if (e1[1] < e0[0])
        {
            if (speed <= 0.0) { return false; }
            double t = (e0[0] - e1[1]) / speed;
            if (t > tFirst) { tFirst = t; }
            if (tFirst > tMax) { return false; }
            t = (e0[1] - e1[0]) / speed;
            if (t < tLast) { tLast = t; }
            if (tFirst > tLast) { return false; }
        }
        else if (e0[1] < e1[0])
        {
            if (speed >= 0.0) { return false; }
            double t = (e0[1] - e1[0]) / speed;
            if (t > tFirst) { tFirst = t; }
            if (tFirst > tMax) { return false; }
            t = (e0[0] - e1[1]) / speed;
            if (t < tLast) { tLast = t; }
            if (tFirst > tLast) { return false; }
        }
        else
        {
            if (speed > 0.0)
            {
                double t = (e0[1] - e1[0]) / speed;
                if (t < tLast) { tLast = t; }
                if (tFirst > tLast) { return false; }
            }
            else if (speed < 0.0)
            {
                double t = (e0[0] - e1[1]) / speed;
                if (t < tLast) { tLast = t; }
                if (tFirst > tLast) { return false; }
            }
        }
        return true;
    }

    bool MtOverlapDir(Triangle3<double> const& t0, Triangle3<double> const& t1,
        Vector3<double> const& dir, double tMax, Vector3<double> const& vel,
        double& tFirst, double& tLast)
    {
        auto e0 = MtProject(t0, dir);
        auto e1 = MtProject(t1, dir);
        double speed = Dot(dir, vel);
        return MtOverlap(tMax, speed, e0, e1, tFirst, tLast);
    }

    bool MovingTri(double tMax, Triangle3<double> const& inT0,
        Vector3<double> const& vel0, Triangle3<double> const& inT1,
        Vector3<double> const& vel1, bool fixedParallelTest, double& contactTime)
    {
        contactTime = 0.0;
        double tFirst = 0.0, tLast = std::numeric_limits<double>::max();
        Vector3<double> relVelocity = vel1 - vel0;
        Vector3<double> const& origin = inT0.v[0];
        Triangle3<double> t0(Vector3<double>{ 0.0, 0.0, 0.0 },
            inT0.v[1] - origin, inT0.v[2] - origin);
        Triangle3<double> t1(inT1.v[0] - origin, inT1.v[1] - origin,
            inT1.v[2] - origin);

        std::array<Vector3<double>, 3> E0{
            t0.v[1] - t0.v[0], t0.v[2] - t0.v[1], t0.v[0] - t0.v[2] };
        Vector3<double> N0 = Cross(E0[0], E0[1]);
        if (!MtOverlapDir(t0, t1, N0, tMax, relVelocity, tFirst, tLast))
        {
            return false;
        }

        std::array<Vector3<double>, 3> E1{
            t1.v[1] - t1.v[0], t1.v[2] - t1.v[1], t1.v[0] - t1.v[2] };
        Vector3<double> N1 = Cross(E1[0], E1[1]);

        bool notParallel{};
        if (fixedParallelTest)
        {
            Vector3<double> N0xN1 = Cross(N0, N1);
            notParallel = (Dot(N0xN1, N0xN1) > 0.0);
        }
        else
        {
            notParallel = (std::fabs(Dot(N0, N1)) < 1.0);
        }

        if (notParallel)
        {
            if (!MtOverlapDir(t0, t1, N1, tMax, relVelocity, tFirst, tLast))
            {
                return false;
            }
            for (size_t i1 = 0; i1 < 3; ++i1)
            {
                for (size_t i0 = 0; i0 < 3; ++i0)
                {
                    Vector3<double> dir = UnitCross(E0[i0], E1[i1]);
                    if (!MtOverlapDir(t0, t1, dir, tMax, relVelocity, tFirst, tLast))
                    {
                        return false;
                    }
                }
            }
        }
        else
        {
            for (size_t i0 = 0; i0 < 3; ++i0)
            {
                Vector3<double> dir = UnitCross(N0, E0[i0]);
                if (!MtOverlapDir(t0, t1, dir, tMax, relVelocity, tFirst, tLast))
                {
                    return false;
                }
            }
            for (size_t i1 = 0; i1 < 3; ++i1)
            {
                Vector3<double> dir = UnitCross(N1, E1[i1]);
                if (!MtOverlapDir(t0, t1, dir, tMax, relVelocity, tFirst, tLast))
                {
                    return false;
                }
            }
        }
        contactTime = tFirst;
        return true;
    }

    // The exact predicate: upstream's parallel test and the port's agree.
    bool MovingParallelTestAgrees(Triangle3<double> const& inT0,
        Triangle3<double> const& inT1)
    {
        Vector3<double> const& origin = inT0.v[0];
        Triangle3<double> t0(Vector3<double>{ 0.0, 0.0, 0.0 },
            inT0.v[1] - origin, inT0.v[2] - origin);
        Triangle3<double> t1(inT1.v[0] - origin, inT1.v[1] - origin,
            inT1.v[2] - origin);
        Vector3<double> N0 = Cross(t0.v[1] - t0.v[0], t0.v[2] - t0.v[1]);
        Vector3<double> N1 = Cross(t1.v[1] - t1.v[0], t1.v[2] - t1.v[1]);
        Vector3<double> N0xN1 = Cross(N0, N1);
        return (std::fabs(Dot(N0, N1)) < 1.0) == (Dot(N0xN1, N0xN1) > 0.0);
    }
}

ORACLE_CASE("IntrTriangle3Triangle3.test")
{
    int mode = io.index() % 2;
    auto triangle0 = Tri<3>(io, mode, 3, 4.0);
    auto triangle1 = Tri<3>(io, mode, 3, 4.0);
    TIQuery<double, Triangle3<double>, Triangle3<double>> query;
    auto r = query(triangle0, triangle1);
    io.outBool(r.intersect);
    io.outReal(r.contactTime);
}

ORACLE_CASE("IntrTriangle3Triangle3.find")
{
    int mode = io.index() % 2;
    auto triangle0 = Tri<3>(io, mode, 3, 4.0);
    auto triangle1 = Tri<3>(io, mode, 3, 4.0);
    FIQuery<double, Triangle3<double>, Triangle3<double>> query;
    auto r = query(triangle0, triangle1);
    EmitTriangleTriangleFI(io, r);
}

// Triangle1 straddles the plane of triangle0 with exactly integer signed
// distances, which reaches every row of the header's (n,p,z) table: the
// coplanar branch (0,0,3), the two-zero edge branch, the (1,1,1) branch and
// the transverse (1,2,0) and (2,1,0) branches.
ORACLE_CASE("IntrTriangle3Triangle3.find.onPlane")
{
    Triangle3<double> triangle0{}, triangle1{};
    AxisPlaneTrianglePair(io, 1, triangle0, triangle1);
    FIQuery<double, Triangle3<double>, Triangle3<double>> query;
    auto r = query(triangle0, triangle1);
    EmitTriangleTriangleFI(io, r);
}

ORACLE_CASE("IntrTriangle3Triangle3.test.onPlane")
{
    Triangle3<double> triangle0{}, triangle1{};
    AxisPlaneTrianglePair(io, 1, triangle0, triangle1);
    TIQuery<double, Triangle3<double>, Triangle3<double>> query;
    auto r = query(triangle0, triangle1);
    io.outBool(r.intersect);
    io.outReal(r.contactTime);
}

// Both triangles are in the same axis plane, so Cross(N0,N1) is exactly zero
// and the coplanar branch of the stationary separating-axis test is taken.
ORACLE_CASE("IntrTriangle3Triangle3.test.coplanar")
{
    Triangle3<double> triangle0{}, triangle1{};
    AxisPlaneTrianglePair(io, 0, triangle0, triangle1);
    TIQuery<double, Triangle3<double>, Triangle3<double>> query;
    auto r = query(triangle0, triangle1);
    io.outBool(r.intersect);
    io.outReal(r.contactTime);
}

namespace
{
    struct MovingData
    {
        std::array<Vector3<double>, 3> v0, v1;
        Vector3<double> vel0, vel1;
        double tMax;
    };

    // mode 0: both triangles in a common axis plane on the integer lattice,
    // so Cross(N0,N1) is exactly zero and |Dot(N0,N1)| >= 1 -- the one regime
    // in which upstream's parallel test is right about parallel triangles.
    // mode 1: uniform vertices in [-range,range]; with range below about 1 the
    // edge cross products are short enough that |Dot(N0,N1)| < 1, which is the
    // regime in which upstream's test is right about non-parallel triangles.
    void RawMoving(oracle::Ctx& io, int mode, double range, MovingData& d)
    {
        if (mode == 0)
        {
            int k = io.rawInteger(0, 2);
            double c = static_cast<double>(io.rawInteger(-3, 3));
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    d.v0[i][j] = static_cast<double>(io.rawInteger(-3, 3));
                    d.v1[i][j] = static_cast<double>(io.rawInteger(-3, 3));
                }
                d.v0[i][k] = c;
                d.v1[i][k] = c;
            }
        }
        else
        {
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    d.v0[i][j] = io.raw(-range, range);
                    d.v1[i][j] = io.raw(-range, range) + range;
                }
            }
        }
        // Aim the relative velocity from the centroid of triangle1 at the
        // centroid of triangle0 with noise, so that a contact occurs within
        // [0,tMax] on a good fraction of the records; a purely random relative
        // velocity almost never brings two small triangles together.
        Vector3<double> c0 = (d.v0[0] + d.v0[1] + d.v0[2]) / 3.0;
        Vector3<double> c1 = (d.v1[0] + d.v1[1] + d.v1[2]) / 3.0;
        Vector3<double> aim = c0 - c1;
        double aimLength = Length(aim);
        if (aimLength > 0.0) { Normalize(aim); }
        double speed = io.raw(0.2, 2.0);
        for (int j = 0; j < 3; ++j) { d.vel0[j] = io.raw(-0.3, 0.3); }
        for (int j = 0; j < 3; ++j)
        {
            d.vel1[j] = d.vel0[j] + speed * aim[j] + io.raw(-0.3, 0.3);
        }
        d.tMax = io.raw(0.25, 4.0);
    }

    void EmitMoving(oracle::Ctx& io, MovingData const& d,
        Triangle3<double>& triangle0, Triangle3<double>& triangle1)
    {
        for (int i = 0; i < 3; ++i) { triangle0.v[i] = io.givenVec(d.v0[i]); }
        for (int i = 0; i < 3; ++i) { triangle1.v[i] = io.givenVec(d.v1[i]); }
        io.givenVec(d.vel0);
        io.givenVec(d.vel1);
        io.given(d.tMax);
    }

    // 'requireDeviation' false keeps the records on which upstream's parallel
    // test and the port's take the same branch (so the whole computation is
    // identical); true keeps the records on which upstream's answer differs
    // from the port's.
    void MovingCase(oracle::Ctx& io, int mode, double range,
        bool requireDeviation, bool findQuery)
    {
        MovingData d{};
        for (int attempt = 0; attempt < 20000; ++attempt)
        {
            RawMoving(io, mode, range, d);
            Triangle3<double> p0{}, p1{};
            for (int i = 0; i < 3; ++i) { p0.v[i] = d.v0[i]; p1.v[i] = d.v1[i]; }
            if (!requireDeviation)
            {
                if (MovingParallelTestAgrees(p0, p1)) { break; }
                continue;
            }
            TIQuery<double, Triangle3<double>, Triangle3<double>> probe;
            auto up = probe(d.tMax, p0, d.vel0, p1, d.vel1);
            double fixedTime = 0.0;
            bool fixedIntersect = MovingTri(d.tMax, p0, d.vel0, p1, d.vel1,
                true, fixedTime);
            if (up.intersect != fixedIntersect
                || (up.intersect && up.contactTime != fixedTime))
            {
                break;
            }
        }
        Triangle3<double> triangle0{}, triangle1{};
        EmitMoving(io, d, triangle0, triangle1);
        if (findQuery)
        {
            FIQuery<double, Triangle3<double>, Triangle3<double>> query;
            auto r = query(d.tMax, triangle0, d.vel0, triangle1, d.vel1);
            EmitTriangleTriangleFI(io, r);
        }
        else
        {
            TIQuery<double, Triangle3<double>, Triangle3<double>> query;
            auto r = query(d.tMax, triangle0, d.vel0, triangle1, d.vel1);
            io.outBool(r.intersect);
            io.outReal(r.contactTime);
        }
    }
}

ORACLE_CASE("IntrTriangle3Triangle3.test.moving")
{
    MovingCase(io, 1, 0.6, false, false);
}

ORACLE_CASE("IntrTriangle3Triangle3.find.moving")
{
    MovingCase(io, 1, 0.6, false, true);
}

ORACLE_CASE("IntrTriangle3Triangle3.test.moving.coplanar")
{
    MovingCase(io, 0, 0.0, false, false);
}

ORACLE_CASE("IntrTriangle3Triangle3.find.moving.coplanar")
{
    MovingCase(io, 0, 0.0, false, true);
}

// Deliberate deviation: triangles large enough that |Dot(N0,N1)| >= 1 although
// Cross(N0,N1) != 0. Upstream calls them parallel and skips the N1 axis and
// all nine E0[i0] x E1[i1] axes; the port runs them.
ORACLE_CASE("IntrTriangle3Triangle3.test.moving.parallelDeviation")
{
    MovingCase(io, 1, 3.0, true, false);
}

ORACLE_CASE("IntrTriangle3Triangle3.find.moving.parallelDeviation")
{
    MovingCase(io, 1, 3.0, true, true);
}

// ================================ cones ==================================
//
// Cone3::SetAngle calls cos, sin and tan, which MSVC and V8 may round
// differently, so the generator draws the angle *unrecorded* and records
// cosAngle, sinAngle and tanAngle. The remaining members are exact functions
// of those three (cosAngleSqr = cosAngle*cosAngle, sinAngleSqr =
// sinAngle*sinAngle, invSinAngle = 1/sinAngle), which is precisely what
// SetAngle computes, so both sides derive them identically and the C math
// library never enters the compared path. Cone3::angle itself is not read by
// any query in this group and is set to zero on both sides.
//
// The cone kind is 0 (infinite), 1 (infinite truncated), 2 (finite) or
// 3 (frustum); minHeight and maxHeight are always recorded, with the
// maxHeight = -1 sentinel of an infinite cone.

namespace
{
    void SetConeTrig(Cone3<double>& cone, double cosA, double sinA, double tanA)
    {
        cone.angle = 0.0;
        cone.cosAngle = cosA;
        cone.sinAngle = sinA;
        cone.tanAngle = tanA;
        cone.cosAngleSqr = cosA * cosA;
        cone.sinAngleSqr = sinA * sinA;
        cone.invSinAngle = 1.0 / sinA;
    }

    Cone3<double> Cn(oracle::Ctx& io, int mode, int kind, int lat, double range)
    {
        Vector3<double> origin = Pt<3>(io, mode, lat, range);
        Vector3<double> direction = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
        double angle = io.raw(0.3, 1.2);
        double cosA = io.given(std::cos(angle));
        double sinA = io.given(std::sin(angle));
        double tanA = io.given(std::tan(angle));
        double rawMin = (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
            : io.raw(0.25, range));
        double rawSpan = (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
            : io.raw(0.5, range));
        double hmin = io.given((kind == 1 || kind == 3) ? rawMin : 0.0);
        double hmax = io.given((kind == 2 || kind == 3) ? hmin + rawSpan : -1.0);
        Cone3<double> cone{};
        cone.ray.origin = origin;
        cone.ray.direction = direction;
        SetConeTrig(cone, cosA, sinA, tanA);
        if (hmax == -1.0) { cone.MakeInfiniteTruncatedCone(hmin); }
        else { cone.MakeConeFrustum(hmin, hmax); }
        return cone;
    }

    // The same cone with the vertex pinned to the origin, which is the one
    // configuration in which upstream's FI 'point' expression (which omits the
    // cone vertex) is bit-identical to the port's corrected one.
    Cone3<double> CnAtOrigin(oracle::Ctx& io, int mode, int kind, double range)
    {
        Vector3<double> origin = io.givenVec(Vector3<double>{ 0.0, 0.0, 0.0 });
        Vector3<double> direction = (mode == 0 ? AxisUnit<3>(io) : io.unit<3>());
        double angle = io.raw(0.3, 1.2);
        double cosA = io.given(std::cos(angle));
        double sinA = io.given(std::sin(angle));
        double tanA = io.given(std::tan(angle));
        double rawMin = io.raw(0.25, range);
        double rawSpan = io.raw(0.5, range);
        double hmin = io.given((kind == 1 || kind == 3) ? rawMin : 0.0);
        double hmax = io.given((kind == 2 || kind == 3) ? hmin + rawSpan : -1.0);
        Cone3<double> cone{};
        cone.ray.origin = origin;
        cone.ray.direction = direction;
        SetConeTrig(cone, cosA, sinA, tanA);
        if (hmax == -1.0) { cone.MakeInfiniteTruncatedCone(hmin); }
        else { cone.MakeConeFrustum(hmin, hmax); }
        return cone;
    }
}

// ============================= IntrSphere3Cone3 ==========================
//
// The TIQuery dispatches on the cone kind to one of four private helpers; the
// four cases below select the kind explicitly so that each helper is covered.
// The sphere radius is drawn as a multiple of the distance from the sphere
// centre to the cone vertex, which splits the answer between true and false.

namespace
{
    Sphere3<double> Sph(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector3<double> center = Pt<3>(io, mode, lat, range);
        double radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        Sphere3<double> sphere{};
        sphere.center = center;
        sphere.radius = radius;
        return sphere;
    }

    void SphereConeTest(oracle::Ctx& io, int mode, int kind)
    {
        auto cone = Cn(io, mode, kind, 3, 3.0);
        auto sphere = Sph(io, mode, 3, 3.0);
        TIQuery<double, Sphere3<double>, Cone3<double>> query;
        auto r = query(sphere, cone);
        io.outBool(r.intersect);
    }
}

ORACLE_CASE("IntrSphere3Cone3.test.infiniteCone")
{
    SphereConeTest(io, io.index() % 2, 0);
}

ORACLE_CASE("IntrSphere3Cone3.test.infiniteTruncatedCone")
{
    SphereConeTest(io, io.index() % 2, 1);
}

ORACLE_CASE("IntrSphere3Cone3.test.finiteCone")
{
    SphereConeTest(io, io.index() % 2, 2);
}

ORACLE_CASE("IntrSphere3Cone3.test.coneFrustum")
{
    SphereConeTest(io, io.index() % 2, 3);
}

// The FI query ignores the cone height range. Its 'point' is the field the
// port corrects (docs/UPSTREAM-FINDINGS.md, issue #307: upstream's
// 'point = t * (cosAngle*D + tmp*B)' omits the cone vertex V), so the general
// case compares 'intersect' only; 'find.originVertex' compares the point where
// the two expressions are bit-identical and 'find.vertexDeviation'
// demonstrates the fix.
ORACLE_CASE("IntrSphere3Cone3.find")
{
    int mode = io.index() % 2;
    int kind = io.rawInteger(0, 3);
    auto cone = Cn(io, mode, kind, 3, 3.0);
    auto sphere = Sph(io, mode, 3, 3.0);
    FIQuery<double, Sphere3<double>, Cone3<double>> query;
    auto r = query(sphere, cone);
    io.outBool(r.intersect);
}

ORACLE_CASE("IntrSphere3Cone3.find.originVertex")
{
    int mode = io.index() % 2;
    int kind = io.rawInteger(0, 3);
    auto cone = CnAtOrigin(io, mode, kind, 3.0);
    auto sphere = Sph(io, mode, 3, 3.0);
    FIQuery<double, Sphere3<double>, Cone3<double>> query;
    auto r = query(sphere, cone);
    io.outBool(r.intersect);
    io.outVec(r.point);
}

// Deliberate deviation: a cone whose vertex is not the origin, on inputs that
// reach the third branch of the query (the one that computes a point from the
// reduced ray/circle problem). Upstream reports t*D and the port V + t*D.
ORACLE_CASE("IntrSphere3Cone3.find.vertexDeviation")
{
    Cone3<double> cone{};
    Sphere3<double> sphere{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        Vector3<double> coneOrigin{
            io.raw(-3.0, 3.0), io.raw(-3.0, 3.0), io.raw(-3.0, 3.0) };
        Vector3<double> coneDir{};
        double len = 0.0;
        do
        {
            for (int i = 0; i < 3; ++i) { coneDir[i] = io.raw(-1.0, 1.0); }
            len = Length(coneDir);
        } while (len < 0.1 || len > 1.0);
        Normalize(coneDir);
        double angle = io.raw(0.3, 1.2);
        Vector3<double> sphereCenter{
            io.raw(-3.0, 3.0), io.raw(-3.0, 3.0), io.raw(-3.0, 3.0) };
        double radius = io.raw(0.25, 3.0);
        cone.ray.origin = coneOrigin;
        cone.ray.direction = coneDir;
        SetConeTrig(cone, std::cos(angle), std::sin(angle), std::tan(angle));
        cone.MakeInfiniteCone();
        sphere.center = sphereCenter;
        sphere.radius = radius;

        // The third branch is reached when the cone vertex is outside the
        // sphere and the sphere centre is outside the (double-sided) cone;
        // it reports a point only when the reduced quadratic has a root.
        Vector3<double> diff = sphere.center - cone.ray.origin;
        double rSqr = radius * radius;
        double lenSqr = Dot(diff, diff);
        if (lenSqr <= rSqr) { continue; }
        double dotAD = Dot(diff, cone.ray.direction);
        double dotSqr = dotAD * dotAD;
        if (dotSqr >= lenSqr * cone.cosAngleSqr && dotAD > 0.0) { continue; }
        double uLen = std::sqrt(std::max(lenSqr - dotSqr, 0.0));
        double test = cone.cosAngle * dotAD + cone.sinAngle * uLen;
        double discr = test * test - lenSqr + rSqr;
        if (discr >= 0.0 && test >= 0.0) { break; }
    }
    io.givenVec(cone.ray.origin);
    io.givenVec(cone.ray.direction);
    io.given(cone.cosAngle);
    io.given(cone.sinAngle);
    io.given(cone.tanAngle);
    io.givenVec(sphere.center);
    io.given(sphere.radius);
    FIQuery<double, Sphere3<double>, Cone3<double>> query;
    auto r = query(sphere, cone);
    io.outBool(r.intersect);
    io.outVec(r.point);
}

// ==================== IntrRay3Cone3 / IntrSegment3Cone3 ==================
//
// Both derive from the line/cone FIQuery and clip its t-interval, and both
// therefore inherit the port's corrections of IntrLine3Cone3
// (docs/UPSTREAM-FINDINGS.md, issues #304 and #465): the through-vertex test
// is on the full vector equation rather than its U-component, and the
// c2 > 0 branches of the discriminant analysis are re-derived, because a line
// through the cone vertex has a double root whose discriminant is a cancelling
// difference. 'LineConePortDeviates' below evaluates, on upstream's own
// quantities and after the direction flip DoQuery performs, exactly the
// condition under which the port takes a different path; the main cases reject
// it and the deviation cases construct it.

namespace
{
    using QFN1 = QFNumber<double, 1>;

    template <typename R>
    void EmitConeFI(oracle::Ctx& io, R const& r)
    {
        io.outBool(r.intersect);
        io.outInt(r.type);
        for (int i = 0; i < 2; ++i)
        {
            io.outReal(r.t[i].x[0]);
            io.outReal(r.t[i].x[1]);
            io.outReal(r.t[i].d);
        }
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                io.outReal(r.P[i][j].x[0]);
                io.outReal(r.P[i][j].x[1]);
                io.outReal(r.P[i][j].d);
            }
        }
    }

    bool LineConeOnVertex(Vector3<double> const& PmV, Vector3<double> const& U,
        double t)
    {
        for (int i = 0; i < 3; ++i)
        {
            if (PmV[i] + t * U[i] != 0.0) { return false; }
        }
        return true;
    }

    bool LineConePortDeviates(Vector3<double> const& lineOrigin,
        Vector3<double> const& inDirection, Cone3<double> const& cone)
    {
        Vector3<double> U = (Dot(inDirection, cone.ray.direction) >= 0.0
            ? inDirection : -inDirection);
        Vector3<double> PmV = lineOrigin - cone.ray.origin;
        double UdU = Dot(U, U);
        double DdU = Dot(cone.ray.direction, U);
        double DdPmV = Dot(cone.ray.direction, PmV);
        double UdPmV = Dot(U, PmV);
        double PmVdPmV = Dot(PmV, PmV);
        double c2 = DdU * DdU - cone.cosAngleSqr * UdU;
        double c1 = DdU * DdPmV - cone.cosAngleSqr * UdPmV;
        double c0 = DdPmV * DdPmV - cone.cosAngleSqr * PmVdPmV;

        // The port hoists a through-vertex test to the top of doQuerySpecial:
        // tv is the parameter of the point of the line closest to the cone
        // vertex, and the line contains the vertex exactly when
        // (P - V) + tv * U is the zero vector.
        double tv = -UdPmV / UdU;
        if (LineConeOnVertex(PmV, U, tv)) { return true; }

        if (c2 != 0.0)
        {
            double discr = c1 * c1 - c0 * c2;
            if (discr < 0.0)
            {
                // The port reports a ray when c2 > 0 (a negative discriminant
                // there is the rounding of a double root); upstream reports
                // the empty set.
                return c2 > 0.0;
            }
            if (discr > 0.0)
            {
                if (c2 > 0.0) { return false; }
                // Block 3 (h[0] < 0 < h[1]) is where the port re-derives the
                // c2 < 0 case.
                double x = -c1 / c2;
                double y = -1.0 / c2;
                QFN1 t0(x, -y, discr), t1(x, y, discr);
                QFN1 h0 = t0 * DdU + DdPmV;
                QFN1 h1 = t1 * DdU + DdPmV;
                QFN1 zero(0.0, 0.0, discr);
                return !(h0 >= zero) && !(h1 <= zero);
            }
            // discr == 0: the port uses the full vector vertex test where
            // upstream uses only its U-component, and takes the ray branch
            // when c2 > 0 and the tangency point is not the vertex.
            double t = -c1 / c2;
            bool onVertexT = LineConeOnVertex(PmV, U, t);
            bool upstreamVertex = (t * UdU + UdPmV == 0.0);
            if (onVertexT != upstreamVertex) { return true; }
            return (!onVertexT && c2 > 0.0);
        }
        // The c2 == 0 branches are transcribed unchanged.
        return false;
    }

    struct ConeData
    {
        Vector3<double> origin, direction;
        double cosA, sinA, tanA, hmin, hmax;
    };

    void RawCone(oracle::Ctx& io, int mode, int kind, int lat, double range,
        ConeData& c)
    {
        for (int i = 0; i < 3; ++i)
        {
            c.origin[i] = (mode == 0
                ? static_cast<double>(io.rawInteger(-lat, lat))
                : io.raw(-range, range));
        }
        if (mode == 0)
        {
            int k = io.rawInteger(0, 5);
            c.direction.MakeZero();
            c.direction[k % 3] = (k < 3 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { c.direction[i] = io.raw(-1.0, 1.0); }
                len = Length(c.direction);
            } while (len < 0.1 || len > 1.0);
            Normalize(c.direction);
        }
        double angle = io.raw(0.3, 1.2);
        c.cosA = std::cos(angle);
        c.sinA = std::sin(angle);
        c.tanA = std::tan(angle);
        double rawMin = (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
            : io.raw(0.25, range));
        double rawSpan = (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
            : io.raw(0.5, range));
        c.hmin = ((kind == 1 || kind == 3) ? rawMin : 0.0);
        c.hmax = ((kind == 2 || kind == 3) ? c.hmin + rawSpan : -1.0);
    }

    Cone3<double> MakeCone(ConeData const& c)
    {
        Cone3<double> cone{};
        cone.ray.origin = c.origin;
        cone.ray.direction = c.direction;
        SetConeTrig(cone, c.cosA, c.sinA, c.tanA);
        if (c.hmax == -1.0) { cone.MakeInfiniteTruncatedCone(c.hmin); }
        else { cone.MakeConeFrustum(c.hmin, c.hmax); }
        return cone;
    }

    void EmitCone(oracle::Ctx& io, ConeData const& c)
    {
        io.givenVec(c.origin);
        io.givenVec(c.direction);
        io.given(c.cosA);
        io.given(c.sinA);
        io.given(c.tanA);
        io.given(c.hmin);
        io.given(c.hmax);
    }

    void RawIntegerDirection(oracle::Ctx& io, Vector3<double>& direction)
    {
        bool zero = true;
        do
        {
            zero = true;
            for (int i = 0; i < 3; ++i)
            {
                direction[i] = static_cast<double>(io.rawInteger(-2, 2));
                zero = zero && direction[i] == 0.0;
            }
        } while (zero);
    }

    void RawOriginDirection(oracle::Ctx& io, int mode, int lat, double range,
        Vector3<double>& origin, Vector3<double>& direction)
    {
        for (int i = 0; i < 3; ++i)
        {
            origin[i] = (mode == 0
                ? static_cast<double>(io.rawInteger(-lat, lat))
                : io.raw(-range, range));
        }
        if (mode == 0)
        {
            RawIntegerDirection(io, direction);
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
    }
}

// In lattice mode the cone axis is a signed coordinate axis with an integer
// vertex and integer heights, and the ray origin and direction are integers,
// so every height Dot(D, P - V) + t * Dot(D, U) is an exact integer and the
// height-range clamps are evaluated at exact equality.
ORACLE_CASE("IntrRay3Cone3.find")
{
    int mode = io.index() % 2;
    int kind = io.rawInteger(0, 3);
    ConeData c{};
    Vector3<double> origin{}, direction{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawCone(io, mode, kind, 3, 3.0, c);
        RawOriginDirection(io, mode, 3, 3.0, origin, direction);
        Cone3<double> probe = MakeCone(c);
        if (!LineConePortDeviates(origin, direction, probe)) { break; }
    }
    EmitCone(io, c);
    io.givenVec(origin);
    io.givenVec(direction);
    Cone3<double> cone = MakeCone(c);
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, Cone3<double>> query;
    auto r = query(ray, cone);
    EmitConeFI(io, r);
}

ORACLE_CASE("IntrSegment3Cone3.find")
{
    int mode = io.index() % 2;
    int kind = io.rawInteger(0, 3);
    ConeData c{};
    Vector3<double> p0{}, p1{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawCone(io, mode, kind, 3, 3.0, c);
        Vector3<double> direction{};
        RawOriginDirection(io, mode, 3, 3.0, p0, direction);
        double s = static_cast<double>(io.rawInteger(1, 4));
        p1 = p0 + s * direction;
        Cone3<double> probe = MakeCone(c);
        if (!LineConePortDeviates(p0, p1 - p0, probe)) { break; }
    }
    EmitCone(io, c);
    io.givenVec(p0);
    io.givenVec(p1);
    Cone3<double> cone = MakeCone(c);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Cone3<double>> query;
    auto r = query(segment, cone);
    EmitConeFI(io, r);
}

namespace
{
    // A point in the solid cone: V + h*D + r*(cos p * U1 + sin p * U2) with h
    // in the height range and 0 <= r <= h*tan(angle). Built from unrecorded
    // draws; the case records the segment endpoints derived from it.
    Vector3<double> RawPointInCone(oracle::Ctx& io, ConeData const& c)
    {
        std::array<Vector3<double>, 3> f{};
        f[0] = c.direction;
        ComputeOrthogonalComplement(1, f.data());
        double hmaxUse = (c.hmax == -1.0 ? c.hmin + 3.0 : c.hmax);
        double h = c.hmin + io.raw(0.0, 1.0) * (hmaxUse - c.hmin);
        double r = io.raw(0.0, 1.0) * (h * c.tanA);
        double phi = io.raw(-kPi, kPi);
        Vector3<double> p = c.origin + h * c.direction;
        p = p + (r * std::cos(phi)) * f[1];
        p = p + (r * std::sin(phi)) * f[2];
        return p;
    }
}

// A segment through a point of the solid cone, so that the query reports an
// intersection on essentially every record and the clipping blocks 21 to 31
// are exercised; the uniform 'find' case above mostly reports the empty set.
ORACLE_CASE("IntrSegment3Cone3.find.throughInterior")
{
    int kind = io.rawInteger(0, 3);
    ConeData c{};
    Vector3<double> p0{}, p1{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawCone(io, 1, kind, 3, 3.0, c);
        Vector3<double> target = RawPointInCone(io, c);
        Vector3<double> u{};
        double len = 0.0;
        do
        {
            for (int i = 0; i < 3; ++i) { u[i] = io.raw(-1.0, 1.0); }
            len = Length(u);
        } while (len < 0.1 || len > 1.0);
        Normalize(u);
        double a = io.raw(0.25, 3.0);
        double b = io.raw(0.25, 3.0);
        p0 = target - a * u;
        p1 = target + b * u;
        Cone3<double> probe = MakeCone(c);
        if (!LineConePortDeviates(p0, p1 - p0, probe)) { break; }
    }
    EmitCone(io, c);
    io.givenVec(p0);
    io.givenVec(p1);
    Cone3<double> cone = MakeCone(c);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Cone3<double>> query;
    auto r = query(segment, cone);
    EmitConeFI(io, r);
}

// Deliberate deviation: the ray's line contains the cone vertex exactly. With
// an integer vertex V, an integer direction U and the origin at V + k * U for
// a nonzero integer k, the vector (P - V) + tv * U is exactly zero, so the
// quadratic has a double root; upstream's discriminant is a cancelling
// difference that rounds to a nonzero value, and its case analysis then
// reports a point, the empty set or a segment reaching the cone's maximum
// height where the intersection is a ray or the vertex alone.
namespace
{
    // The configuration in which the port's correction is visible after the
    // ray or segment clip: the line contains the cone vertex, c2 > 0 (the
    // direction is interior to the cone's angular region, so the true
    // intersection with the positive cone is the ray of nonnegative heights
    // from the vertex), the vertex parameter lies in the clipped range
    // [lo,hi], and the cone's minimum height is zero (otherwise both sides
    // clamp the ray to the same hmin plane and agree).
    bool ThroughVertexRayConfiguration(Vector3<double> const& origin,
        Vector3<double> const& direction, Cone3<double> const& cone,
        double lo, double hi)
    {
        if (cone.GetMinHeight() != 0.0) { return false; }
        Vector3<double> U = (Dot(direction, cone.ray.direction) >= 0.0
            ? direction : -direction);
        double UdU = Dot(U, U);
        double DdU = Dot(cone.ray.direction, U);
        double c2 = DdU * DdU - cone.cosAngleSqr * UdU;
        if (c2 <= 0.0) { return false; }
        Vector3<double> PmV = origin - cone.ray.origin;
        double tv = -Dot(direction, PmV) / Dot(direction, direction);
        return lo <= tv && tv <= hi;
    }
}

ORACLE_CASE("IntrRay3Cone3.find.throughVertexDeviation")
{
    ConeData c{};
    Vector3<double> direction{}, origin{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        int kind = (io.rawInteger(0, 1) == 0 ? 0 : 2);
        RawCone(io, 0, kind, 3, 3.0, c);
        RawIntegerDirection(io, direction);
        int k = io.rawInteger(1, 3);
        int sign = io.rawInteger(0, 1);
        double kk = static_cast<double>(sign == 0 ? -k : k);
        origin = c.origin + kk * direction;
        Cone3<double> probe = MakeCone(c);
        if (ThroughVertexRayConfiguration(origin, direction, probe, 0.0,
            std::numeric_limits<double>::max()))
        {
            break;
        }
    }
    EmitCone(io, c);
    io.givenVec(origin);
    io.givenVec(direction);
    Cone3<double> cone = MakeCone(c);
    Ray3<double> ray(origin, direction);
    FIQuery<double, Ray3<double>, Cone3<double>> query;
    auto r = query(ray, cone);
    EmitConeFI(io, r);
}

// The same construction for the segment. The segment length along the integer
// direction is a power of two, so the centred quantity
// (P - V) + tv * (p1 - p0) is exactly zero as well.
ORACLE_CASE("IntrSegment3Cone3.find.throughVertexDeviation")
{
    ConeData c{};
    Vector3<double> p0{}, p1{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        int kind = (io.rawInteger(0, 1) == 0 ? 0 : 2);
        RawCone(io, 0, kind, 3, 3.0, c);
        Vector3<double> direction{};
        RawIntegerDirection(io, direction);
        int k = io.rawInteger(1, 3);
        int sign = io.rawInteger(0, 1);
        double kk = static_cast<double>(sign == 0 ? -k : k);
        double L = std::ldexp(1.0, io.rawInteger(0, 2));
        p0 = c.origin + kk * direction;
        p1 = p0 + L * direction;
        Cone3<double> probe = MakeCone(c);
        if (ThroughVertexRayConfiguration(p0, p1 - p0, probe, 0.0, 1.0))
        {
            break;
        }
    }
    EmitCone(io, c);
    io.givenVec(p0);
    io.givenVec(p1);
    Cone3<double> cone = MakeCone(c);
    Segment3<double> segment(p0, p1);
    FIQuery<double, Segment3<double>, Cone3<double>> query;
    auto r = query(segment, cone);
    EmitConeFI(io, r);
}

// ========================== IntrOrientedBox3Cone3 ========================
//
// A thin wrapper that transforms the cone into the box frame and calls the
// aligned-box/cone query, from which it derives; it therefore inherits that
// query's one deliberate port fix (docs/UPSTREAM-FINDINGS.md, issue #301):
// upstream's 'BoxFullyInConeSlab' copies the twelve box edges over
// mCandidateEdges and sets mNumCandidateEdges to 12 *without* touching
// mAdjacencyMatrix, so adjacency bits that an earlier clipping query set for
// the clipped vertices (indices 8..31) are never cleared. A later
// ClearCandidates only visits the twelve box edges, so those bits survive,
// InsertEdge then treats the corresponding candidate edges as already present
// and drops them, and the query reports a false negative.
//
// The defect needs a query object that is reused, so the main case uses a
// fresh object per record (which is what every other case in the family does)
// and the deviation case runs three calls on one object: a clipping
// configuration, a fully-inside-the-slab configuration, and the first
// configuration again.

ORACLE_CASE("IntrOrientedBox3Cone3.test")
{
    int mode = io.index() % 2;
    int kind = io.rawInteger(0, 3);
    auto box = OBox3(io, mode, 3, 3.0);
    auto cone = Cn(io, mode, kind, 3, 3.0);
    TIQuery<double, OrientedBox3<double>, Cone3<double>> query;
    auto r = query(box, cone);
    io.outBool(r.intersect);
}

namespace
{
    struct StaleData
    {
        Vector3<double> boxCenter, boxExtent;
        std::array<Vector3<double>, 3> boxAxis;
        ConeData coneA, coneB;
    };

    void RawStale(oracle::Ctx& io, StaleData& d)
    {
        std::array<Vector3<double>, 3> frame{};
        double len = 0.0;
        do
        {
            for (int i = 0; i < 3; ++i) { frame[0][i] = io.raw(-1.0, 1.0); }
            len = Length(frame[0]);
        } while (len < 0.1 || len > 1.0);
        Normalize(frame[0]);
        ComputeOrthogonalComplement(1, frame.data());
        d.boxAxis = frame;
        for (int i = 0; i < 3; ++i)
        {
            d.boxCenter[i] = io.raw(-2.0, 2.0);
            d.boxExtent[i] = io.raw(0.5, 1.5);
        }

        std::array<Vector3<double>, 3> coneFrame{};
        do
        {
            for (int i = 0; i < 3; ++i) { coneFrame[0][i] = io.raw(-1.0, 1.0); }
            len = Length(coneFrame[0]);
        } while (len < 0.1 || len > 1.0);
        Normalize(coneFrame[0]);
        ComputeOrthogonalComplement(1, coneFrame.data());
        Vector3<double> D = coneFrame[0];
        Vector3<double> perp = coneFrame[1];

        double angleA = io.raw(0.4, 1.2);
        double angleB = io.raw(0.4, 1.2);
        double a = io.raw(0.5, 2.5);
        double offsetA = io.raw(0.5, 2.5);
        double extra = io.raw(1.0, 3.0);
        double offsetB = io.raw(0.5, 2.5);
        double sumExtent = d.boxExtent[0] + d.boxExtent[1] + d.boxExtent[2];

        // Cone A is finite with its maximum-height plane through the box
        // centre, so the box straddles that plane and the clipping path runs.
        d.coneA.origin = d.boxCenter - a * D + offsetA * perp;
        d.coneA.direction = D;
        d.coneA.cosA = std::cos(angleA);
        d.coneA.sinA = std::sin(angleA);
        d.coneA.tanA = std::tan(angleA);
        d.coneA.hmin = 0.0;
        d.coneA.hmax = a;

        // Cone B is infinite with its vertex far enough below the box that
        // every box height is positive, so the box is fully inside the slab;
        // the perpendicular offset keeps the cone axis out of the box, which
        // would otherwise be a quick acceptance.
        d.coneB.origin = d.boxCenter - (sumExtent + extra) * D + offsetB * perp;
        d.coneB.direction = D;
        d.coneB.cosA = std::cos(angleB);
        d.coneB.sinA = std::sin(angleB);
        d.coneB.tanA = std::tan(angleB);
        d.coneB.hmin = 0.0;
        d.coneB.hmax = -1.0;
    }

    void EmitStale(oracle::Ctx& io, StaleData const& d)
    {
        io.givenVec(d.boxCenter);
        io.givenVec(d.boxAxis[0]);
        io.givenVec(d.boxAxis[1]);
        io.givenVec(d.boxAxis[2]);
        io.givenVec(d.boxExtent);
        EmitCone(io, d.coneA);
        EmitCone(io, d.coneB);
    }
}

// Deliberate deviation: the stale adjacency bits above. The record keeps the
// inputs on which the third call's answer differs from the answer a fresh
// query object gives for the same inputs.
ORACLE_CASE("IntrOrientedBox3Cone3.test.staleAdjacencyDeviation")
{
    StaleData d{};
    for (int attempt = 0; attempt < 20000; ++attempt)
    {
        RawStale(io, d);
        OrientedBox3<double> probeBox(d.boxCenter, d.boxAxis, d.boxExtent);
        Cone3<double> probeA = MakeCone(d.coneA);
        Cone3<double> probeB = MakeCone(d.coneB);
        TIQuery<double, OrientedBox3<double>, Cone3<double>> reused;
        reused(probeBox, probeA);
        reused(probeBox, probeB);
        bool third = reused(probeBox, probeA).intersect;
        TIQuery<double, OrientedBox3<double>, Cone3<double>> fresh;
        bool reference = fresh(probeBox, probeA).intersect;
        if (third != reference) { break; }
    }
    EmitStale(io, d);
    OrientedBox3<double> box(d.boxCenter, d.boxAxis, d.boxExtent);
    Cone3<double> coneA = MakeCone(d.coneA);
    Cone3<double> coneB = MakeCone(d.coneB);
    TIQuery<double, OrientedBox3<double>, Cone3<double>> query;
    io.outBool(query(box, coneA).intersect);
    io.outBool(query(box, coneB).intersect);
    io.outBool(query(box, coneA).intersect);
}
