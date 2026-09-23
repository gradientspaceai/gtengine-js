// Verify group 10 (computational geometry): differential cases for
// ExtremalQuery3BSP.h, IncrementalDelaunay2.h, MinimumAreaCircle2.h,
// MinimumVolumeSphere3.h and RotatingCalipers.h.
//
// WHAT IS COMPARABLE HERE. Three of the five headers route through standard
// library facilities whose order is not fixed by the standard:
//
//  * MinimumAreaCircle2 / MinimumVolumeSphere3 call
//    std::shuffle(permuted, mDRE) with a default-constructed
//    std::default_random_engine. On MSVC that type is std::mt19937 (checked
//    with typeid); the port uses a minstd_rand0-style Lehmer generator, so
//    the two permutations differ. The permutation decides the *order* in
//    which the support circles/spheres are built, and exactCircle3 /
//    exactSphere4 are not symmetric in their arguments, so a general input
//    cannot be compared bit for bit. The main cases therefore accept only
//    point sets whose result is invariant under the permutation, which is
//    tested by calling the real upstream query 48 times on one object: the
//    engine state advances between calls, so each call uses a different
//    shuffle. Measured acceptance: 63% (uniform doubles) to 78% (dense
//    lattice) of draws.
//
//  * ExtremalQuery3BSP builds its BSP tree from VETManifoldMesh, whose
//    GetEdges(), GetVertices() and Vertex::TAdjacent are std::unordered_*
//    containers (TAdjacent is even keyed by Triangle*). The port reads them
//    through sorted accessors, so the two BSP trees are different (measured:
//    the octahedron gives 12 nodes upstream and 15 in the port). Only the
//    query answer is compared, and only on directions whose extreme vertex
//    is unique - a tie is resolved by the leaf region the direction lands
//    in, which is tree shape. See "Not covered" in the group report.
//
//  * IncrementalDelaunay2's mGraph is a VETManifoldMesh as well, so the
//    triangle numbering follows hash-table order. The triangle *set* and
//    each triangle's stored vertex rotation are order independent (the
//    Delaunay triangulation is unique for points in general position and
//    VETManifoldMesh::Insert keeps the rotation the caller passed), so the
//    cases sort the triangles by their stored vertex tuple and remap the
//    adjacency indices through the same permutation, exactly as v09 does for
//    Delaunay2/Delaunay3. GetHull is unaffected: it walks a std::map.
//
// Everything compared is either combinatorial or computed with + - * / sqrt
// only (the circle/sphere centers and radii, the face normals through
// UnitCross), so every non-deviation case is declared { exact: true }.
//
// One io draw per C++ statement: MSVC evaluates function arguments right to
// left, so every generated value goes into a named local before it is used.
#define ORACLE_FAMILY "v10-compgeom"
#include "Oracle.h"

#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/ConvexHull2.h>
#include <Mathematics/ExtremalQuery3BSP.h>
#include <Mathematics/IncrementalDelaunay2.h>
#include <Mathematics/MinimumAreaCircle2.h>
#include <Mathematics/MinimumVolumeSphere3.h>
#include <Mathematics/RotatingCalipers.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <exception>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

using namespace gte;

namespace
{
    // BSNumber<UIntegerAP32> grows as needed and the expressions below use
    // only + - *, so these signs are exact for any double input. They are
    // used by the generator probes and by the independent reference checks;
    // nothing they compute reaches an output.
    using Exact = BSNumber<UIntegerAP32>;

    // Sign of DotPerp(B - A, C - A): +1 when <A,B,C> is counterclockwise.
    int32_t ExactOrient2(Vector2<double> const& A, Vector2<double> const& B,
        Vector2<double> const& C)
    {
        Exact ax(A[0]), ay(A[1]), bx(B[0]), by(B[1]), cx(C[0]), cy(C[1]);
        Exact x0 = bx - ax, y0 = by - ay, x1 = cx - ax, y1 = cy - ay;
        Exact det = x0 * y1 - x1 * y0;
        return det.GetSign();
    }

    // ---- lattice point generators ---------------------------------------
    // The 12 lattice points of the circle x^2 + y^2 = 25; a subset of them
    // is in convex position and cocircular, which is the degenerate input
    // the circle and Delaunay queries care about.
    std::array<std::array<int32_t, 2>, 12> const gCircle25
    { {
        { 5, 0 }, { 4, 3 }, { 3, 4 }, { 0, 5 }, { -3, 4 }, { -4, 3 },
        { -5, 0 }, { -4, -3 }, { -3, -4 }, { 0, -5 }, { 3, -4 }, { 4, -3 }
    } };

    // The counterclockwise order of gCircle25 (it is listed by decreasing
    // polar angle start, i.e. already counterclockwise from (5,0)).
    std::vector<Vector2<double>> Circle25Subset(oracle::Ctx& io, size_t k,
        double cx, double cy, double scale)
    {
        // Choose k of the 12 cocircular points, keeping their cyclic order.
        std::array<bool, 12> take{};
        size_t taken = 0;
        for (int32_t attempt = 0; attempt < 200 && taken < k; ++attempt)
        {
            int32_t j = io.rawInteger(0, 11);
            if (!take[static_cast<size_t>(j)])
            {
                take[static_cast<size_t>(j)] = true;
                ++taken;
            }
        }
        std::vector<Vector2<double>> pts{};
        for (size_t j = 0; j < 12; ++j)
        {
            if (take[j])
            {
                pts.push_back(Vector2<double>{
                    cx + scale * static_cast<double>(gCircle25[j][0]),
                    cy + scale * static_cast<double>(gCircle25[j][1]) });
            }
        }
        return pts;
    }

    // A counterclockwise convex polygon with distinct vertices, no three of
    // them collinear.
    //   mode 0: k cocircular lattice points of a scaled circle of radius 5;
    //           every caliper angle comparison is then an exact tie
    //           candidate
    //   mode 1: the convex hull of random lattice points
    //   mode 2: the convex hull of random real points
    // Returns an empty vector when the draw did not produce a polygon with
    // at least three vertices; the caller redraws.
    std::vector<Vector2<double>> RawConvexPolygon(oracle::Ctx& io, int32_t mode)
    {
        std::vector<Vector2<double>> poly{};
        if (mode == 0)
        {
            size_t k = static_cast<size_t>(io.rawInteger(3, 8));
            double cx = static_cast<double>(io.rawInteger(-3, 3));
            double cy = static_cast<double>(io.rawInteger(-3, 3));
            double scale = static_cast<double>(io.rawInteger(1, 3));
            poly = Circle25Subset(io, k, cx, cy, scale);
            if (poly.size() < 3)
            {
                poly.clear();
            }
            return poly;
        }

        size_t n = static_cast<size_t>(io.rawInteger(4, 10));
        std::vector<Vector2<double>> pts(n);
        for (size_t i = 0; i < n; ++i)
        {
            if (mode == 1)
            {
                pts[i][0] = static_cast<double>(io.rawInteger(-6, 6));
                pts[i][1] = static_cast<double>(io.rawInteger(-6, 6));
            }
            else
            {
                pts[i][0] = io.raw(-6.0, 6.0);
                pts[i][1] = io.raw(-6.0, 6.0);
            }
        }

        ConvexHull2<double> hull{};
        if (!hull(pts) || hull.GetDimension() != 2)
        {
            return poly;
        }
        auto const& indices = hull.GetHull();
        if (indices.size() < 3)
        {
            return poly;
        }
        for (auto index : indices)
        {
            poly.push_back(pts[static_cast<size_t>(index)]);
        }

        // ConvexHull2 can leave collinear hull points in the list and the
        // rotating calipers algorithm is specified for a polygon with no
        // duplicate and no collinear vertices, so reject those draws: the
        // main case is about the generic path. (The collinear-removal path of
        // CreatePolygon is covered by the .collinear case below.)
        size_t const m = poly.size();
        for (size_t i1 = 0; i1 < m; ++i1)
        {
            size_t i0 = (i1 + m - 1) % m;
            size_t i2 = (i1 + 1) % m;
            if (ExactOrient2(poly[i0], poly[i1], poly[i2]) <= 0)
            {
                poly.clear();
                return poly;
            }
        }
        return poly;
    }

    // Record a polygon and emit the antipode list of upstream's query.
    void EmitAntipodes(oracle::Ctx& io, std::vector<Vector2<double>> const& poly)
    {
        std::vector<RotatingCalipers<double>::Antipode> antipodes{};
        RotatingCalipers<double>::ComputeAntipodes(poly, antipodes);
        io.outInt(antipodes.size());
        for (auto const& a : antipodes)
        {
            io.outInt(a.vertex);
            io.outInt(a.edge[0]);
            io.outInt(a.edge[1]);
        }
    }
}

// ---- RotatingCalipers ----------------------------------------------------

// ComputeAntipodes is the only public entry point. Everything on the path is
// exact rational arithmetic (BSNumber) over the input coordinates, so the
// antipode list is a pure combinatorial function of the input and is
// compared exactly, in upstream's own order.
//
// The cocircular mode (index % 3 == 0) puts every vertex on one circle, which
// makes the caliper angle comparisons of ComputeNextAntipode tie as often as
// possible; the tie-duplication quirk of issue #286 (an edge emitted twice
// and another omitted) is preserved by the port and is compared here rather
// than excluded.
ORACLE_CASE("RotatingCalipers.computeAntipodes")
{
    int32_t mode = io.index() % 3;
    std::vector<Vector2<double>> poly{};
    for (int32_t attempt = 0; attempt < 32 && poly.empty(); ++attempt)
    {
        poly = RawConvexPolygon(io, mode);
    }
    if (poly.empty())
    {
        // The fallback is a fixed lattice triangle; every mode can reach the
        // cap, and the record must still contain a polygon.
        poly = { Vector2<double>{ 0.0, 0.0 }, Vector2<double>{ 4.0, 0.0 },
            Vector2<double>{ 0.0, 3.0 } };
    }
    io.given(static_cast<double>(poly.size()));
    for (auto const& v : poly)
    {
        io.givenVec<2>(v);
    }
    EmitAntipodes(io, poly);
}

// CreatePolygon's collinear removal: extra vertices are placed on the edges
// of a convex polygon at the exactly representable midpoint, so the arriving
// and leaving edges are exactly parallel and DotPerp is exactly zero. No
// duplicate point is inserted, so the port's most-recent-nonzero-edge fix of
// issue #286 is inactive and the two implementations must agree.
ORACLE_CASE("RotatingCalipers.computeAntipodes.collinear")
{
    int32_t mode = io.index() % 3;
    std::vector<Vector2<double>> poly{};
    for (int32_t attempt = 0; attempt < 32 && poly.empty(); ++attempt)
    {
        poly = RawConvexPolygon(io, mode);
    }
    if (poly.empty())
    {
        poly = { Vector2<double>{ 0.0, 0.0 }, Vector2<double>{ 4.0, 0.0 },
            Vector2<double>{ 0.0, 3.0 } };
    }

    // Subdivide a random subset of the edges at their midpoints. Halving is
    // exact in binary floating point, so the inserted vertex is exactly on
    // the segment and DotPerp of the two edges it creates is exactly zero.
    std::vector<Vector2<double>> refined{};
    size_t const m = poly.size();
    for (size_t i0 = 0; i0 < m; ++i0)
    {
        size_t i1 = (i0 + 1) % m;
        refined.push_back(poly[i0]);
        if (io.rawInteger(0, 1) == 1)
        {
            Vector2<double> mid
            {
                0.5 * (poly[i0][0] + poly[i1][0]),
                0.5 * (poly[i0][1] + poly[i1][1])
            };
            refined.push_back(mid);
        }
    }

    io.given(static_cast<double>(refined.size()));
    for (auto const& v : refined)
    {
        io.givenVec<2>(v);
    }
    EmitAntipodes(io, refined);
}

// Deliberate port fix of issue #286: CreatePolygon tests collinearity against
// the immediately preceding edge of the input array, so a duplicated vertex
// produces a zero-length edge whose DotPerp is zero against the next edge and
// the *next real corner* is discarded with the duplicate. The port compares
// against the most recent nonzero edge instead. This case duplicates one
// vertex of a convex polygon, which upstream answers with one corner fewer
// (and, when that drops the retained count below three, with a LogAssert).
ORACLE_CASE("RotatingCalipers.computeAntipodes.deviation.duplicate")
{
    int32_t mode = io.index() % 3;
    std::vector<Vector2<double>> poly{};
    for (int32_t attempt = 0; attempt < 32 && poly.empty(); ++attempt)
    {
        poly = RawConvexPolygon(io, mode);
    }
    if (poly.empty())
    {
        poly = { Vector2<double>{ 0.0, 0.0 }, Vector2<double>{ 4.0, 0.0 },
            Vector2<double>{ 0.0, 3.0 } };
    }

    size_t const m = poly.size();
    size_t dup = static_cast<size_t>(io.rawInteger(0, static_cast<int32_t>(m) - 1));
    std::vector<Vector2<double>> withDup{};
    for (size_t i = 0; i < m; ++i)
    {
        withDup.push_back(poly[i]);
        if (i == dup)
        {
            withDup.push_back(poly[i]);
        }
    }

    io.given(static_cast<double>(withDup.size()));
    for (auto const& v : withDup)
    {
        io.givenVec<2>(v);
    }
    EmitAntipodes(io, withDup);
}

// ---- MinimumAreaCircle2 / MinimumVolumeSphere3 ---------------------------

namespace
{
    bool SameBits(double a, double b)
    {
        uint64_t ua, ub;
        std::memcpy(&ua, &a, 8);
        std::memcpy(&ub, &b, 8);
        return ua == ub;
    }

    struct Result2
    {
        bool ok;
        double cx, cy, radius;
        int32_t numSupport;
        std::array<int32_t, 3> support;

        bool Same(Result2 const& o) const
        {
            return ok == o.ok && SameBits(cx, o.cx) && SameBits(cy, o.cy)
                && SameBits(radius, o.radius) && numSupport == o.numSupport
                && support == o.support;
        }
    };

    Result2 RunMAC(MinimumAreaCircle2<double, double>& q,
        std::vector<Vector2<double>> const& pts)
    {
        Circle2<double> circle{};
        Result2 r{};
        r.ok = q(static_cast<int32_t>(pts.size()), pts.data(), circle);
        r.cx = circle.center[0];
        r.cy = circle.center[1];
        r.radius = circle.radius;
        r.numSupport = q.GetNumSupport();
        r.support = q.GetSupport();
        std::sort(r.support.begin(), r.support.begin() + r.numSupport);
        for (int32_t i = r.numSupport; i < 3; ++i)
        {
            r.support[i] = -1;
        }
        return r;
    }

    struct Result3
    {
        bool ok;
        double cx, cy, cz, radius;
        int32_t numSupport;
        std::array<int32_t, 4> support;

        bool Same(Result3 const& o) const
        {
            return ok == o.ok && SameBits(cx, o.cx) && SameBits(cy, o.cy)
                && SameBits(cz, o.cz) && SameBits(radius, o.radius)
                && numSupport == o.numSupport && support == o.support;
        }
    };

    Result3 RunMVS(MinimumVolumeSphere3<double, double>& q,
        std::vector<Vector3<double>> const& pts)
    {
        Sphere3<double> sphere{};
        Result3 r{};
        r.ok = q(static_cast<int32_t>(pts.size()), pts.data(), sphere);
        r.cx = sphere.center[0];
        r.cy = sphere.center[1];
        r.cz = sphere.center[2];
        r.radius = sphere.radius;
        r.numSupport = q.GetNumSupport();
        r.support = q.GetSupport();
        std::sort(r.support.begin(), r.support.begin() + r.numSupport);
        for (int32_t i = r.numSupport; i < 4; ++i)
        {
            r.support[i] = -1;
        }
        return r;
    }

    // The engine member advances between calls, so calling the same object
    // repeatedly samples different shuffles of the same point set. A point
    // set whose result is the same for all of them is one on which the
    // implementation-defined permutation cannot be observed, and it is the
    // only kind this case can compare bit for bit.
    int32_t const gNumShuffles = 48;

    bool Invariant2(std::vector<Vector2<double>> const& pts, Result2& first)
    {
        MinimumAreaCircle2<double, double> q{};
        first = RunMAC(q, pts);
        if (!first.ok)
        {
            // The trapped-failure fallback is a deliberate port fix of
            // issue #286; the .deviation case owns it.
            return false;
        }
        for (int32_t k = 1; k < gNumShuffles; ++k)
        {
            Result2 r = RunMAC(q, pts);
            if (!r.Same(first))
            {
                return false;
            }
        }
        return true;
    }

    bool Invariant3(std::vector<Vector3<double>> const& pts, Result3& first)
    {
        MinimumVolumeSphere3<double, double> q{};
        first = RunMVS(q, pts);
        if (!first.ok)
        {
            return false;
        }
        for (int32_t k = 1; k < gNumShuffles; ++k)
        {
            Result3 r = RunMVS(q, pts);
            if (!r.Same(first))
            {
                return false;
            }
        }
        return true;
    }

    // mode 0: uniform doubles; the generic path
    // mode 1: lattice [-4,4]; exact arithmetic, cocircular triples, ties
    // mode 2: cocircular lattice about a lattice center, center and nearby
    //         lattice points mixed in
    // mode 3: dense lattice [-2,2]; duplicates and degeneracies are common
    void RawPoints2(oracle::Ctx& io, int32_t mode, size_t n,
        std::vector<Vector2<double>>& pts)
    {
        pts.resize(n);
        if (mode == 2)
        {
            double cx = static_cast<double>(io.rawInteger(-2, 2));
            double cy = static_cast<double>(io.rawInteger(-2, 2));
            for (size_t i = 0; i < n; ++i)
            {
                int32_t k = io.rawInteger(0, 13);
                if (k < 12)
                {
                    pts[i][0] = cx + static_cast<double>(gCircle25[k][0]);
                    pts[i][1] = cy + static_cast<double>(gCircle25[k][1]);
                }
                else if (k == 12)
                {
                    pts[i][0] = cx;
                    pts[i][1] = cy;
                }
                else
                {
                    pts[i][0] = cx + static_cast<double>(io.rawInteger(-2, 2));
                    pts[i][1] = cy + static_cast<double>(io.rawInteger(-2, 2));
                }
            }
            return;
        }

        for (size_t i = 0; i < n; ++i)
        {
            if (mode == 0)
            {
                pts[i][0] = io.raw(-4.0, 4.0);
                pts[i][1] = io.raw(-4.0, 4.0);
            }
            else if (mode == 1)
            {
                pts[i][0] = static_cast<double>(io.rawInteger(-4, 4));
                pts[i][1] = static_cast<double>(io.rawInteger(-4, 4));
            }
            else
            {
                pts[i][0] = static_cast<double>(io.rawInteger(-2, 2));
                pts[i][1] = static_cast<double>(io.rawInteger(-2, 2));
            }
        }
    }

    void RawPoints3(oracle::Ctx& io, int32_t mode, size_t n,
        std::vector<Vector3<double>>& pts)
    {
        pts.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            for (int32_t j = 0; j < 3; ++j)
            {
                if (mode == 0)
                {
                    pts[i][j] = io.raw(-4.0, 4.0);
                }
                else if (mode == 1)
                {
                    pts[i][j] = static_cast<double>(io.rawInteger(-4, 4));
                }
                else
                {
                    pts[i][j] = static_cast<double>(io.rawInteger(-2, 2));
                }
            }
        }
    }
}

// The full public surface of MinimumAreaCircle2<double,double>: operator(),
// GetNumSupport() and GetSupport(). The support indices are sorted before
// they are emitted: upstream writes them in the order the update functions
// happened to fill mSupport, which depends on the shuffle even when the
// support *set* does not.
//
// The generator keeps at most 8 points, which matters twice: the invariance
// probe costs 48 queries per draw, and MSVC's std::sort is an insertion sort
// (hence stable) for at most 32 elements, which is what makes upstream's
// choice of representative among equal points - the first in sorted order -
// agree with the port's explicit index tie-break.
ORACLE_CASE("MinimumAreaCircle2.compute")
{
    int32_t mode = io.index() % 4;
    int32_t n = io.integer(1, 8);
    std::vector<Vector2<double>> pts{};
    Result2 res{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 24 && !accepted; ++attempt)
    {
        RawPoints2(io, mode, static_cast<size_t>(n), pts);
        accepted = Invariant2(pts, res);
    }
    if (!accepted)
    {
        // Points on a parabola are in general position; the minimum-area
        // circle of such a set is supported by two or three of them with a
        // wide margin, so the result does not depend on the permutation.
        pts.resize(static_cast<size_t>(n));
        for (int32_t i = 0; i < n; ++i)
        {
            double t = static_cast<double>(i) - 2.0;
            pts[static_cast<size_t>(i)][0] = t;
            pts[static_cast<size_t>(i)][1] = t * t;
        }
        MinimumAreaCircle2<double, double> q{};
        res = RunMAC(q, pts);
    }

    for (int32_t i = 0; i < n; ++i)
    {
        io.givenVec<2>(pts[static_cast<size_t>(i)]);
    }

    io.outBool(res.ok);
    io.outReal(res.cx);
    io.outReal(res.cy);
    io.outReal(res.radius);
    io.outInt(res.numSupport);
    for (int32_t i = 0; i < res.numSupport; ++i)
    {
        io.outInt(res.support[i]);
    }
}

// Deliberate port fix of issue #286: in the trapped-failure branch upstream
// calls GetContainer(numPoints, points, minimal) after numPoints has been
// overwritten with the *unique* point count, so the fallback circle bounds
// only a prefix of the input array. The port passes the whole array.
//
// The trapped branch is reached about 3 times in 10000 random lattice draws,
// far too rarely for a rejection loop, so the case draws from a catalogue of
// four-point lattice sets that were found to reach it (collected with the
// real upstream query; see the group report). The points are permuted and a
// random number of duplicates of one of them is prefixed, which leaves the
// set of distinct points - and therefore the control flow - unchanged while
// making the unique count smaller than the array length, so upstream's
// prefix misses at least one genuine point.
ORACLE_CASE("MinimumAreaCircle2.compute.deviation.trappedFallback")
{
    static std::array<std::array<double, 8>, 16> const gTrapped2
    { {
        { -1, 4, -3, 0, 2, 4, 4, 0 },
        { 4, 2, 0, 4, 0, -3, 4, -1 },
        { -3, -4, 2, 3, -3, 4, 2, -3 },
        { -2, -1, 0, 2, 3, -4, 4, 1 },
        { -2, 3, 2, 1, -2, -4, 2, -2 },
        { 2, 2, -3, -2, -1, 2, 4, -2 },
        { 2, -3, 2, 3, -3, 4, -3, -4 },
        { 1, -3, -2, -3, 3, 1, -4, 1 },
        { -4, -4, 3, 0, -1, 2, -4, 0 },
        { -4, 2, 4, 2, 3, -3, -3, -3 },
        { -1, -4, 4, 0, -3, 0, 2, -4 },
        { -3, -4, -4, 2, 0, 2, -1, -4 },
        { -2, -1, 4, 1, 3, -4, 0, 2 },
        { 3, -3, 4, 2, -4, 2, -3, -3 },
        { 3, -1, 1, 3, -2, 3, -4, -1 },
        { -2, 1, 0, -3, 4, 4, 0, 4 }
    } };

    int32_t which = io.rawInteger(0, 15);
    auto const& raw = gTrapped2[static_cast<size_t>(which)];
    std::vector<Vector2<double>> base(4);
    for (size_t i = 0; i < 4; ++i)
    {
        base[i][0] = raw[2 * i];
        base[i][1] = raw[2 * i + 1];
    }

    std::array<int32_t, 4> order{ 0, 1, 2, 3 };
    for (int32_t i = 3; i > 0; --i)
    {
        int32_t j = io.rawInteger(0, i);
        std::swap(order[static_cast<size_t>(i)], order[static_cast<size_t>(j)]);
    }
    int32_t numDup = io.rawInteger(1, 3);
    int32_t dupOf = io.rawInteger(0, 3);

    std::vector<Vector2<double>> pts{};
    for (int32_t k = 0; k < numDup; ++k)
    {
        pts.push_back(base[static_cast<size_t>(dupOf)]);
    }
    for (int32_t i = 0; i < 4; ++i)
    {
        pts.push_back(base[static_cast<size_t>(order[static_cast<size_t>(i)])]);
    }

    io.given(static_cast<double>(pts.size()));
    for (auto const& p : pts)
    {
        io.givenVec<2>(p);
    }

    MinimumAreaCircle2<double, double> q{};
    Result2 res = RunMAC(q, pts);
    io.outBool(res.ok);
    io.outReal(res.cx);
    io.outReal(res.cy);
    io.outReal(res.radius);
    io.outInt(res.numSupport);
    for (int32_t i = 0; i < res.numSupport; ++i)
    {
        io.outInt(res.support[i]);
    }
}

// Throw parity: operator() calls LogError when there are no points.
ORACLE_CASE("MinimumAreaCircle2.compute.empty")
{
    Vector2<double> unused = io.latticeVec<2>(-2, 2);
    (void)unused;
    MinimumAreaCircle2<double, double> q{};
    Circle2<double> circle{};
    bool ok = q(0, static_cast<Vector2<double> const*>(nullptr), circle);
    io.outBool(ok);
}

// The full public surface of MinimumVolumeSphere3<double,double>. Same
// comparison policy as the 2D query above.
ORACLE_CASE("MinimumVolumeSphere3.compute")
{
    int32_t mode = io.index() % 3;
    int32_t n = io.integer(1, 8);
    std::vector<Vector3<double>> pts{};
    Result3 res{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 24 && !accepted; ++attempt)
    {
        RawPoints3(io, mode, static_cast<size_t>(n), pts);
        accepted = Invariant3(pts, res);
    }
    if (!accepted)
    {
        // Points on the moment curve are in general position.
        pts.resize(static_cast<size_t>(n));
        for (int32_t i = 0; i < n; ++i)
        {
            double t = static_cast<double>(i) - 2.0;
            pts[static_cast<size_t>(i)][0] = t;
            pts[static_cast<size_t>(i)][1] = t * t;
            pts[static_cast<size_t>(i)][2] = t * t * t;
        }
        MinimumVolumeSphere3<double, double> q{};
        res = RunMVS(q, pts);
    }

    for (int32_t i = 0; i < n; ++i)
    {
        io.givenVec<3>(pts[static_cast<size_t>(i)]);
    }

    io.outBool(res.ok);
    io.outReal(res.cx);
    io.outReal(res.cy);
    io.outReal(res.cz);
    io.outReal(res.radius);
    io.outInt(res.numSupport);
    for (int32_t i = 0; i < res.numSupport; ++i)
    {
        io.outInt(res.support[i]);
    }
}

// The 3D twin of the trapped-fallback deviation; see the 2D case.
ORACLE_CASE("MinimumVolumeSphere3.compute.deviation.trappedFallback")
{
    static std::array<std::array<double, 15>, 16> const gTrapped3
    { {
        { 2, 0, 1, 1, -1, 1, 0, 1, 2, -1, 2, -1, 1, 1, -2 },
        { 0, -1, 1, 0, -2, -1, -1, 1, 0, 2, 0, -1, 0, 1, -2 },
        { 0, 0, 2, 0, 1, 0, 2, -1, 0, 0, -2, -1, -1, -2, 1 },
        { 1, -1, 0, -2, -2, 2, -1, -2, -2, 1, 1, 2, 1, 0, -2 },
        { 2, -2, 2, 1, 2, 2, -1, -2, -1, -1, -2, 0, -1, 2, 0 },
        { 2, -1, -1, -1, -1, -2, 1, 0, 2, -2, 0, 1, 1, 2, -1 },
        { 1, 2, 2, -1, -2, 2, 2, 0, -1, 0, 2, 2, -1, 0, -1 },
        { 1, -2, 0, 1, -2, -1, -1, 1, 2, 2, 1, -1, -1, -2, 1 },
        { -2, 1, 0, 1, -1, 0, -1, 0, -2, -1, -2, 0, 0, 1, 1 },
        { 0, 1, 0, -1, 0, -1, 1, 2, -2, -2, 2, -2, -2, 1, -2 },
        { -2, 2, -1, -1, 2, 2, 2, 1, 0, -2, -2, -1, 2, -1, 0 },
        { 1, -1, 1, 0, 0, -2, 1, 1, -2, 1, -2, -1, 2, 0, 1 },
        { 2, 0, 1, -2, 2, -1, 1, 2, 1, -2, -1, -2, -1, -1, -1 },
        { -1, 1, 2, -1, -2, -1, 2, 2, 1, 2, -2, 2, 0, 2, -1 },
        { 1, 2, 2, 0, 2, 0, 1, 1, 0, 0, 1, 1, 2, 2, 1 },
        { -2, -1, 1, -1, -1, 2, 0, 0, 1, -1, 0, 0, -2, 0, 2 }
    } };

    int32_t which = io.rawInteger(0, 15);
    auto const& raw = gTrapped3[static_cast<size_t>(which)];
    std::vector<Vector3<double>> base(5);
    for (size_t i = 0; i < 5; ++i)
    {
        base[i][0] = raw[3 * i];
        base[i][1] = raw[3 * i + 1];
        base[i][2] = raw[3 * i + 2];
    }

    std::array<int32_t, 5> order{ 0, 1, 2, 3, 4 };
    for (int32_t i = 4; i > 0; --i)
    {
        int32_t j = io.rawInteger(0, i);
        std::swap(order[static_cast<size_t>(i)], order[static_cast<size_t>(j)]);
    }
    int32_t numDup = io.rawInteger(1, 3);
    int32_t dupOf = io.rawInteger(0, 4);

    std::vector<Vector3<double>> pts{};
    for (int32_t k = 0; k < numDup; ++k)
    {
        pts.push_back(base[static_cast<size_t>(dupOf)]);
    }
    for (int32_t i = 0; i < 5; ++i)
    {
        pts.push_back(base[static_cast<size_t>(order[static_cast<size_t>(i)])]);
    }

    io.given(static_cast<double>(pts.size()));
    for (auto const& p : pts)
    {
        io.givenVec<3>(p);
    }

    MinimumVolumeSphere3<double, double> q{};
    Result3 res = RunMVS(q, pts);
    io.outBool(res.ok);
    io.outReal(res.cx);
    io.outReal(res.cy);
    io.outReal(res.cz);
    io.outReal(res.radius);
    io.outInt(res.numSupport);
    for (int32_t i = 0; i < res.numSupport; ++i)
    {
        io.outInt(res.support[i]);
    }
}

// Throw parity: operator() calls LogError when there are no points.
ORACLE_CASE("MinimumVolumeSphere3.compute.empty")
{
    Vector3<double> unused = io.latticeVec<3>(-2, 2);
    (void)unused;
    MinimumVolumeSphere3<double, double> q{};
    Sphere3<double> sphere{};
    bool ok = q(0, static_cast<Vector3<double> const*>(nullptr), sphere);
    io.outBool(ok);
}

// ---- ExtremalQuery3BSP ---------------------------------------------------

namespace
{
    struct BasePolytope
    {
        std::vector<Vector3<double>> vertices;
        std::vector<int32_t> indices;
    };

    // Three strictly convex simplicial polytopes. Strict convexity matters:
    // when two triangles of a face are coplanar their face normals are equal
    // and the arc normal Cross(N0, N1) is the zero vector, after which every
    // isign(Dot(D, 0)) is zero and the BSP degenerates (a box triangulated
    // with two triangles per face answers even unique-argmax directions
    // wrongly, on both sides and differently). Such polytopes are outside
    // the algorithm's precondition and are not generated.
    BasePolytope const& BasePolytopeOf(int32_t which)
    {
        static std::array<BasePolytope, 3> const gBases
        { {
            // Regular tetrahedron on the lattice.
            {
                { { 1.0, 1.0, 1.0 }, { 1.0, -1.0, -1.0 },
                  { -1.0, 1.0, -1.0 }, { -1.0, -1.0, 1.0 } },
                { 0, 1, 2,  0, 2, 3,  0, 3, 1,  1, 3, 2 }
            },
            // Octahedron.
            {
                { { 1.0, 0.0, 0.0 }, { -1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 },
                  { 0.0, -1.0, 0.0 }, { 0.0, 0.0, 1.0 }, { 0.0, 0.0, -1.0 } },
                { 0, 2, 4,  2, 1, 4,  1, 3, 4,  3, 0, 4,
                  2, 0, 5,  1, 2, 5,  3, 1, 5,  0, 3, 5 }
            },
            // Triangular bipyramid (five vertices, six faces).
            {
                { { 0.0, 0.0, 3.0 }, { 0.0, 0.0, -2.0 }, { 2.0, 0.0, 0.0 },
                  { -1.0, 2.0, 0.0 }, { -1.0, -2.0, 0.0 } },
                { 0, 2, 3,  0, 3, 4,  0, 4, 2,  1, 3, 2,  1, 4, 3,  1, 2, 4 }
            }
        } };
        return gBases[static_cast<size_t>(which)];
    }

    // Exact sign of Dot(D, P) - Dot(D, Q).
    int32_t ExactDotCompare(Vector3<double> const& D, Vector3<double> const& P,
        Vector3<double> const& Q)
    {
        Exact s(0.0);
        for (int32_t j = 0; j < 3; ++j)
        {
            s = s + Exact(D[j]) * (Exact(P[j]) - Exact(Q[j]));
        }
        return s.GetSign();
    }

    // The brute-force extreme vertex: the unique argmax of Dot(D, V) over the
    // polytope vertices, or -1 when the maximum is attained more than once.
    // Ties are decided by the leaf region the direction falls into, which is
    // BSP tree shape, so they are not comparable and are rejected.
    int32_t UniqueArgMax(std::vector<Vector3<double>> const& verts,
        Vector3<double> const& D)
    {
        int32_t best = 0;
        bool tied = false;
        for (size_t i = 1; i < verts.size(); ++i)
        {
            int32_t s = ExactDotCompare(D, verts[i], verts[static_cast<size_t>(best)]);
            if (s > 0)
            {
                best = static_cast<int32_t>(i);
                tied = false;
            }
            else if (s == 0)
            {
                tied = true;
            }
        }
        return tied ? -1 : best;
    }

    // A direction is comparable when
    //  (a) the extreme vertex is unique for both D and -D (a tie is decided
    //      by the leaf region the direction lands in, which is tree shape),
    //      and
    //  (b) upstream's BSP answers both with that unique vertex.
    // Condition (b) excludes the directions on which the preserved
    // construction defect of issue #290 fires. It is not confined to the
    // icosahedra of the original measurement: a sheared octahedron (the
    // octahedron below under an integer linear map) answers about 2% of
    // random directions with a vertex that is not extreme. The port's BSP,
    // built from the same arcs in a different order, is right on those
    // directions, so nothing can be compared there. See the group report.
    bool SoundDirection(ExtremalQuery3BSP<double>& query,
        std::vector<Vector3<double>> const& verts, Vector3<double> const& D)
    {
        if (D[0] == 0.0 && D[1] == 0.0 && D[2] == 0.0)
        {
            return false;
        }
        Vector3<double> negD{ -D[0], -D[1], -D[2] };
        int32_t pos = UniqueArgMax(verts, D);
        int32_t neg = UniqueArgMax(verts, negD);
        if (pos < 0 || neg < 0)
        {
            return false;
        }
        int32_t bspPos = -1, bspNeg = -1;
        query.GetExtremeVertices(D, bspPos, bspNeg);
        return bspPos == pos && bspNeg == neg;
    }

    // Apply an integer linear map with positive determinant, which preserves
    // convexity, the face lattice and the counterclockwise orientation while
    // producing a new set of face normals.
    void TransformPolytope(oracle::Ctx& io, BasePolytope const& base,
        std::vector<Vector3<double>>& verts)
    {
        std::array<std::array<double, 3>, 3> M{};
        double det = 0.0;
        for (int32_t attempt = 0; attempt < 32 && det <= 0.0; ++attempt)
        {
            for (int32_t r = 0; r < 3; ++r)
            {
                for (int32_t c = 0; c < 3; ++c)
                {
                    M[static_cast<size_t>(r)][static_cast<size_t>(c)] =
                        static_cast<double>(io.rawInteger(r == c ? 1 : -1,
                            r == c ? 2 : 1));
                }
            }
            det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
        }
        if (det <= 0.0)
        {
            M = { { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } } };
        }

        verts.resize(base.vertices.size());
        for (size_t i = 0; i < base.vertices.size(); ++i)
        {
            for (int32_t r = 0; r < 3; ++r)
            {
                verts[i][r] =
                    M[static_cast<size_t>(r)][0] * base.vertices[i][0]
                    + M[static_cast<size_t>(r)][1] * base.vertices[i][1]
                    + M[static_cast<size_t>(r)][2] * base.vertices[i][2];
            }
        }
    }
}

// GetExtremeVertices, plus GetFaceNormals() of the ExtremalQuery3 base (pure
// UnitCross arithmetic, compared bit for bit).
//
// NOT compared: GetNumNodes() and GetTreeDepth(). The BSP tree is built from
// VETManifoldMesh's std::unordered_map / std::unordered_set<Triangle*>
// containers, whose iteration order MSVC decides from pointer hashes, while
// the port reads them through sorted accessors. The two trees are different
// partitions of the same Gauss map: measured on the octahedron above,
// upstream builds 12 nodes and the port 15.
//
// Directions whose extreme vertex is not unique are rejected for the same
// reason: a tie is resolved by which leaf region the direction lands in.
ORACLE_CASE("ExtremalQuery3BSP.getExtremeVertices")
{
    int32_t which = io.integer(0, 2);
    BasePolytope const& base = BasePolytopeOf(which);
    std::vector<Vector3<double>> verts{};
    TransformPolytope(io, base, verts);
    for (auto const& v : verts)
    {
        io.givenVec<3>(v);
    }

    auto pool = std::make_shared<std::vector<Vector3<double>>>(verts);
    std::vector<int32_t> indices = base.indices;
    Polyhedron3<double> polytope(pool, static_cast<int32_t>(indices.size()),
        indices.data(), true);
    ExtremalQuery3BSP<double> query(polytope);

    auto const& normals = query.GetFaceNormals();
    io.outInt(normals.size());
    for (auto const& n : normals)
    {
        io.outVec(n);
    }

    int32_t const numDirections = 6;
    for (int32_t k = 0; k < numDirections; ++k)
    {
        Vector3<double> D{ 0.0, 0.0, 0.0 };
        bool accepted = false;
        for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
        {
            if (k % 2 == 0)
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    D[j] = static_cast<double>(io.rawInteger(-4, 4));
                }
            }
            else
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    D[j] = io.raw(-1.0, 1.0);
                }
            }
            accepted = SoundDirection(query, verts, D);
        }
        if (!accepted)
        {
            // Walk a fixed list of small lattice directions for a sound one.
            for (int32_t a = -1; a <= 1 && !accepted; ++a)
            {
                for (int32_t b = -1; b <= 1 && !accepted; ++b)
                {
                    for (int32_t c = -1; c <= 1 && !accepted; ++c)
                    {
                        D = { static_cast<double>(3 * a + 1),
                            static_cast<double>(5 * b + 2),
                            static_cast<double>(7 * c + 4) };
                        accepted = SoundDirection(query, verts, D);
                    }
                }
            }
        }
        io.givenVec<3>(D);

        int32_t posVertex = -1, negVertex = -1;
        query.GetExtremeVertices(D, posVertex, negVertex);
        io.outInt(posVertex);
        io.outInt(negVertex);
    }
}

// ---- IncrementalDelaunay2 ------------------------------------------------

namespace
{
    using ID2 = IncrementalDelaunay2<double>;
    size_t const gInvalid = ID2::invalid;

    // size_t(-1) does not survive a double, and the port uses -1 for the
    // same sentinel, so the invalid index is emitted as -1.
    double AsIndex(size_t v)
    {
        return v == gInvalid ? -1.0 : static_cast<double>(v);
    }

    // The triangle numbering of GetTriangles()/GetAdjacencies() follows the
    // iteration order of VETManifoldMesh's std::unordered_map, which MSVC
    // decides from the TriangleKey hash while the port reads a sorted list.
    // The triangle *set* and each triangle's stored vertex rotation are
    // order independent (VETManifoldMesh::Insert keeps the rotation the
    // caller passed and each triangle is inserted once), so the features are
    // emitted sorted by their stored vertex tuple with the adjacency indices
    // remapped through the same permutation. This is exactly what the v09
    // cases do for Delaunay2/Delaunay3.
    void EmitTriangles(oracle::Ctx& io, ID2& del)
    {
        auto const& tris = del.GetTriangles();
        auto const& adjs = del.GetAdjacencies();
        size_t const count = tris.size();
        std::vector<size_t> order(count);
        std::iota(order.begin(), order.end(), size_t(0));
        std::sort(order.begin(), order.end(), [&tris](size_t a, size_t b)
        {
            return tris[a] < tris[b];
        });
        std::vector<size_t> rank(count);
        for (size_t r = 0; r < count; ++r)
        {
            rank[order[r]] = r;
        }

        io.outInt(count);
        for (size_t r = 0; r < count; ++r)
        {
            size_t t = order[r];
            for (size_t j = 0; j < 3; ++j)
            {
                io.outInt(tris[t][j]);
            }
            for (size_t j = 0; j < 3; ++j)
            {
                size_t a = adjs[t][j];
                io.outReal(a == gInvalid ? -1.0 : static_cast<double>(rank[a]));
            }
        }
    }

    void EmitHull(oracle::Ctx& io, ID2& del)
    {
        std::vector<size_t> hull{};
        del.GetHull(hull);
        io.outInt(hull.size());
        for (auto v : hull)
        {
            io.outInt(v);
        }
    }

    // Points strictly inside [-5,5]^2, which is strictly inside every domain
    // rectangle this family draws.
    //   mode 0: uniform doubles
    //   mode 1: lattice [-5,5]
    //   mode 2: cocircular lattice points of the circle of radius 5 about a
    //           lattice center, plus the center itself
    //   mode 3: dense lattice [-2,2]; duplicate insertions are common
    Vector2<double> RawInsidePoint(oracle::Ctx& io, int32_t mode)
    {
        Vector2<double> p{ 0.0, 0.0 };
        if (mode == 0)
        {
            p[0] = io.raw(-5.0, 5.0);
            p[1] = io.raw(-5.0, 5.0);
        }
        else if (mode == 1)
        {
            p[0] = static_cast<double>(io.rawInteger(-5, 5));
            p[1] = static_cast<double>(io.rawInteger(-5, 5));
        }
        else if (mode == 2)
        {
            int32_t k = io.rawInteger(0, 12);
            if (k < 12)
            {
                p[0] = static_cast<double>(gCircle25[k][0]);
                p[1] = static_cast<double>(gCircle25[k][1]);
            }
        }
        else
        {
            p[0] = static_cast<double>(io.rawInteger(-2, 2));
            p[1] = static_cast<double>(io.rawInteger(-2, 2));
        }
        return p;
    }

    // Draw and record the domain rectangle. Every generated point lies in
    // [-5,5]^2, hence strictly inside.
    std::array<double, 4> DomainRectangle(oracle::Ctx& io)
    {
        double xMin = io.lattice(-10, -6);
        double yMin = io.lattice(-10, -6);
        double xMax = io.lattice(6, 10);
        double yMax = io.lattice(6, 10);
        return { xMin, yMin, xMax, yMax };
    }
}

// Insert: the return index (an existing vertex index for a repeated
// position), then the whole triangulation state - GetNumVertices,
// GetNumTriangles, GetTriangles, GetAdjacencies and GetHull.
ORACLE_CASE("IncrementalDelaunay2.insert")
{
    int32_t mode = io.index() % 4;
    auto rect = DomainRectangle(io);
    int32_t n = io.integer(1, 8);
    ID2 del(rect[0], rect[1], rect[2], rect[3]);
    for (int32_t i = 0; i < n; ++i)
    {
        Vector2<double> p = RawInsidePoint(io, mode);
        io.givenVec<2>(p);
        size_t index = del.Insert(p);
        io.outReal(AsIndex(index));
    }
    io.outInt(del.GetNumVertices());
    io.outInt(del.GetNumTriangles());
    EmitTriangles(io, del);
    EmitHull(io, del);
}

// Insert, then Remove: positions that are vertices, positions that are not,
// and positions of the enclosing rectangle are all exercised. The return
// value of Remove is the vertex index, or invalid when the position is not a
// vertex of the triangulation.
ORACLE_CASE("IncrementalDelaunay2.remove")
{
    int32_t mode = io.index() % 4;
    auto rect = DomainRectangle(io);
    int32_t n = io.integer(3, 8);
    int32_t numRemove = io.integer(1, 4);
    ID2 del(rect[0], rect[1], rect[2], rect[3]);
    std::vector<Vector2<double>> inserted{};
    for (int32_t i = 0; i < n; ++i)
    {
        Vector2<double> p = RawInsidePoint(io, mode);
        io.givenVec<2>(p);
        del.Insert(p);
        inserted.push_back(p);
    }
    for (int32_t k = 0; k < numRemove; ++k)
    {
        // Two thirds of the removals hit an inserted position, one third is
        // a fresh draw that is usually not a vertex.
        Vector2<double> p{ 0.0, 0.0 };
        if (io.rawInteger(0, 2) != 0)
        {
            size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
            p = inserted[j];
        }
        else
        {
            p = RawInsidePoint(io, mode);
        }
        io.givenVec<2>(p);
        size_t index = del.Remove(p);
        io.outReal(AsIndex(index));
    }
    io.outInt(del.GetNumVertices());
    io.outInt(del.GetNumTriangles());
    EmitTriangles(io, del);
    EmitHull(io, del);
}

// GetContainingTriangle, GetTriangle, GetAdjacent and GetTriangulation.
//
// The SearchInfo path and finalTriangle index are NOT emitted: they are
// triangle indices in hash-table numbering and the walk starts at triangle 0
// of that numbering, so both the path and, for a query point outside the
// hull, the exit edge depend on the order. What is emitted is the outcome
// (the stored vertex triple of the containing triangle, or -1), which is
// order independent as long as the containing triangle is unique. Query
// points are accepted only when they are strictly inside exactly one
// triangle (a point on a shared edge stops the walk at whichever of the two
// triangles it reaches first) or strictly outside every triangle.
ORACLE_CASE("IncrementalDelaunay2.getContainingTriangle")
{
    int32_t mode = io.index() % 4;
    auto rect = DomainRectangle(io);
    int32_t n = io.integer(1, 8);
    ID2 del(rect[0], rect[1], rect[2], rect[3]);
    for (int32_t i = 0; i < n; ++i)
    {
        Vector2<double> p = RawInsidePoint(io, mode);
        io.givenVec<2>(p);
        del.Insert(p);
    }

    auto const& tris = del.GetTriangles();
    auto const& verts = del.GetVertices();

    // GetTriangulation: all vertices and all triangles of the graph, the
    // supervertex triangles included. The triangle keys are sorted tuples
    // (TriangleKey<true> orders its indices) and the container is a hash map,
    // so the list is sorted lexicographically on both sides.
    std::vector<Vector2<double>> tvertices{};
    std::vector<std::array<size_t, 3>> ttriangles{};
    del.GetTriangulation(tvertices, ttriangles);
    io.outInt(tvertices.size());
    for (auto const& v : tvertices)
    {
        io.outVec(v);
    }
    std::sort(ttriangles.begin(), ttriangles.end());
    io.outInt(ttriangles.size());
    for (auto const& t : ttriangles)
    {
        io.outInt(t[0]);
        io.outInt(t[1]);
        io.outInt(t[2]);
    }

    int32_t const numQueries = 4;
    for (int32_t k = 0; k < numQueries; ++k)
    {
        Vector2<double> q{ 0.0, 0.0 };
        bool accepted = false;
        for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
        {
            q[0] = io.raw(-7.0, 7.0);
            q[1] = io.raw(-7.0, 7.0);
            int32_t inside = 0;
            for (auto const& t : tris)
            {
                int32_t s0 = ExactOrient2(verts[t[0]], verts[t[1]], q);
                int32_t s1 = ExactOrient2(verts[t[1]], verts[t[2]], q);
                int32_t s2 = ExactOrient2(verts[t[2]], verts[t[0]], q);
                if (s0 > 0 && s1 > 0 && s2 > 0)
                {
                    ++inside;
                }
                else if (s0 >= 0 && s1 >= 0 && s2 >= 0)
                {
                    // On the boundary of a triangle: the walk can stop at
                    // either side of a shared edge.
                    inside = 2;
                    break;
                }
            }
            accepted = (inside <= 1);
        }
        io.givenVec<2>(q);

        ID2::SearchInfo info{};
        size_t t = del.GetContainingTriangle(q, info);
        // The index itself is hash-table numbering; only whether a
        // containing triangle was found, and which triangle it is, are
        // order independent.
        io.outBool(t != gInvalid);
        if (t != gInvalid)
        {
            std::array<size_t, 3> triangle{};
            bool valid = del.GetTriangle(t, triangle);
            io.outBool(valid);
            io.outInt(triangle[0]);
            io.outInt(triangle[1]);
            io.outInt(triangle[2]);
        }
    }

    // Out-of-range accessors return false.
    std::array<size_t, 3> triangle{}, adjacent{};
    bool validTriangle = del.GetTriangle(tris.size(), triangle);
    io.outBool(validTriangle);
    bool validAdjacent = del.GetAdjacent(tris.size(), adjacent);
    io.outBool(validAdjacent);
}

// FinalizeTriangulation removes the four rectangle vertices; the remaining
// Delaunay triangles are those of the inserted points alone. Insert and
// Remove then return 'invalid' and a second FinalizeTriangulation returns
// false. The inserted points must not be collinear, otherwise there are no
// Delaunay triangles left and GetHull is the undefined-behaviour path of
// issue #290 (the .deviation case below).
ORACLE_CASE("IncrementalDelaunay2.finalizeTriangulation")
{
    int32_t mode = io.index() % 4;
    auto rect = DomainRectangle(io);
    int32_t n = io.integer(3, 8);
    std::vector<Vector2<double>> pts{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        pts.clear();
        for (int32_t i = 0; i < n; ++i)
        {
            pts.push_back(RawInsidePoint(io, mode));
        }
        // Accept only when three of the points are not collinear, which is
        // what makes the finalized triangulation nonempty.
        for (size_t i = 0; i < pts.size() && !accepted; ++i)
        {
            for (size_t j = i + 1; j < pts.size() && !accepted; ++j)
            {
                for (size_t k = j + 1; k < pts.size() && !accepted; ++k)
                {
                    accepted = ExactOrient2(pts[i], pts[j], pts[k]) != 0;
                }
            }
        }
    }
    if (!accepted)
    {
        pts.clear();
        for (int32_t i = 0; i < n; ++i)
        {
            double t = static_cast<double>(i) - 2.0;
            pts.push_back(Vector2<double>{ t, t * t });
        }
    }

    ID2 del(rect[0], rect[1], rect[2], rect[3]);
    for (int32_t i = 0; i < n; ++i)
    {
        io.givenVec<2>(pts[static_cast<size_t>(i)]);
        del.Insert(pts[static_cast<size_t>(i)]);
    }

    bool finalized = del.FinalizeTriangulation();
    io.outBool(finalized);
    io.outInt(del.GetNumVertices());
    io.outInt(del.GetNumTriangles());
    EmitTriangles(io, del);
    EmitHull(io, del);

    bool again = del.FinalizeTriangulation();
    io.outBool(again);
    size_t insertAfter = del.Insert(pts[0]);
    io.outReal(AsIndex(insertAfter));
    size_t removeAfter = del.Remove(pts[0]);
    io.outReal(AsIndex(removeAfter));
}

// Deliberate port fix of issue #290: GetHull walks the edge map with an
// unbounded while (vNext != vStart) loop that writes into hull[], which was
// sized to the number of edges. For a finalized triangulation whose input
// points are collinear there are no Delaunay triangles at all, yet the
// triangles sharing a supervertex still contribute edges, and those edges
// form a path ending in a 2-cycle rather than a closed polygon. Upstream
// then loops forever and writes past the end of the vector - undefined
// behaviour that cannot be executed inside the generator.
//
// The case therefore runs a VERBATIM COPY of upstream's edge collection
// (GetGraph() exposes everything it reads) followed by upstream's walk with
// a step cap, and emits the number of edges and the prefix of hull[] that
// upstream would write before it runs off the end. The port's getHull()
// detects the open walk and throws, so every record of this case is a
// disagreement: that is the demonstration of the fix.
ORACLE_CASE("IncrementalDelaunay2.getHull.deviation.collinear")
{
    auto rect = DomainRectangle(io);
    int32_t n = io.integer(2, 6);

    // A lattice line through a lattice point, with lattice multiples of a
    // lattice direction: every inserted point is exactly collinear.
    int32_t bx = 0, by = 0, dx = 0, dy = 0;
    for (int32_t attempt = 0; attempt < 32 && dx == 0 && dy == 0; ++attempt)
    {
        bx = io.rawInteger(-2, 2);
        by = io.rawInteger(-2, 2);
        dx = io.rawInteger(-2, 2);
        dy = io.rawInteger(-2, 2);
    }
    if (dx == 0 && dy == 0)
    {
        dx = 1;
        dy = 2;
    }

    ID2 del(rect[0], rect[1], rect[2], rect[3]);
    for (int32_t i = 0; i < n; ++i)
    {
        double k = static_cast<double>(io.rawInteger(-2, 2));
        Vector2<double> p{ static_cast<double>(bx) + k * static_cast<double>(dx),
            static_cast<double>(by) + k * static_cast<double>(dy) };
        io.givenVec<2>(p);
        del.Insert(p);
    }
    bool finalized = del.FinalizeTriangulation();
    io.outBool(finalized);
    io.outInt(del.GetNumTriangles());

    // Verbatim copy of GetHull's edge collection.
    std::map<size_t, size_t> edges{};
    auto const& vmap = del.GetGraph().GetVertices();
    for (int32_t v = 0; v < 3; ++v)
    {
        auto vIter = vmap.find(v);
        if (vIter == vmap.end())
        {
            continue;
        }
        for (auto const& adj : vIter->second->TAdjacent)
        {
            for (size_t i0 = 1, i1 = 2, i2 = 0; i2 < 3; i0 = i1, i1 = i2, ++i2)
            {
                if (adj->V[i0] == v)
                {
                    if (adj->V[i1] >= 3 && adj->V[i2] >= 3)
                    {
                        edges.insert(std::make_pair(
                            static_cast<size_t>(adj->V[i2]),
                            static_cast<size_t>(adj->V[i1])));
                        break;
                    }
                }
            }
        }
    }

    if (edges.empty())
    {
        // Upstream dereferences edges.begin() on an empty map, which is
        // undefined behaviour and cannot be executed. The port returns an
        // empty hull, and so does this record: those records agree and the
        // defect is only described in the group report.
        io.outInt(0);
        return;
    }

    // Upstream's walk, with a step cap in place of the unbounded loop. The
    // LogAssert below is upstream's own; the cap replaces the out-of-bounds
    // write that follows when the edges do not form a closed cycle through
    // the smallest key, and the emitted prefix is what upstream writes
    // before it overruns hull[].
    auto eIter = edges.begin();
    size_t vStart = eIter->first;
    size_t vNext = eIter->second;
    std::vector<size_t> hull{ vStart };
    size_t const cap = edges.size() + 4;
    while (vNext != vStart && hull.size() < cap)
    {
        hull.push_back(vNext);
        auto it = edges.find(vNext);
        LogAssert(it != edges.end(), "Expecting to find a hull edge.");
        vNext = it->second;
    }
    io.outInt(hull.size());
    for (auto v : hull)
    {
        io.outInt(v);
    }
}

// ---- a 2D case aimed at the three-point support --------------------------

namespace
{
    // Sign of the in-circle determinant for a counterclockwise triangle
    // <A,B,C>: +1 when D is strictly inside the circumcircle.
    int32_t ExactInCircle(Vector2<double> const& A, Vector2<double> const& B,
        Vector2<double> const& C, Vector2<double> const& D)
    {
        std::array<Exact, 3> x{ Exact(A[0]) - Exact(D[0]), Exact(B[0]) - Exact(D[0]),
            Exact(C[0]) - Exact(D[0]) };
        std::array<Exact, 3> y{ Exact(A[1]) - Exact(D[1]), Exact(B[1]) - Exact(D[1]),
            Exact(C[1]) - Exact(D[1]) };
        std::array<Exact, 3> w{ x[0] * x[0] + y[0] * y[0], x[1] * x[1] + y[1] * y[1],
            x[2] * x[2] + y[2] * y[2] };
        Exact det = x[0] * (y[1] * w[2] - y[2] * w[1])
            - y[0] * (x[1] * w[2] - x[2] * w[1])
            + w[0] * (x[1] * y[2] - x[2] * y[1]);
        return det.GetSign();
    }

    // Exact test that the triangle is acute: the dot product of the two edges
    // at each vertex is positive. The minimum-area circle of an acute
    // triangle is its circumcircle, so the support set has three points.
    bool ExactAcute(Vector2<double> const& A, Vector2<double> const& B,
        Vector2<double> const& C)
    {
        auto dot = [](Vector2<double> const& P, Vector2<double> const& Q,
            Vector2<double> const& R)
        {
            Exact ux = Exact(Q[0]) - Exact(P[0]), uy = Exact(Q[1]) - Exact(P[1]);
            Exact vx = Exact(R[0]) - Exact(P[0]), vy = Exact(R[1]) - Exact(P[1]);
            return (ux * vx + uy * vy).GetSign();
        };
        return dot(A, B, C) > 0 && dot(B, A, C) > 0 && dot(C, A, B) > 0;
    }
}

// The general case above ends with a three-point support in about 4% of its
// records, because the permutation-invariance filter prefers the symmetric
// two-point circles. This case aims at the three-point support directly: an
// acute lattice triangle (its minimum-area circle is its circumcircle) plus
// extra points that the exact in-circle predicate places strictly inside.
// ExactCircle3's center and squared radius are then the emitted result.
ORACLE_CASE("MinimumAreaCircle2.compute.circumcircle")
{
    int32_t extra = io.integer(0, 4);
    std::vector<Vector2<double>> pts{};
    Result2 res{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        pts.clear();
        Vector2<double> A{ 0.0, 0.0 }, B{ 0.0, 0.0 }, C{ 0.0, 0.0 };
        bool acute = false;
        for (int32_t k = 0; k < 32 && !acute; ++k)
        {
            A[0] = static_cast<double>(io.rawInteger(-6, 6));
            A[1] = static_cast<double>(io.rawInteger(-6, 6));
            B[0] = static_cast<double>(io.rawInteger(-6, 6));
            B[1] = static_cast<double>(io.rawInteger(-6, 6));
            C[0] = static_cast<double>(io.rawInteger(-6, 6));
            C[1] = static_cast<double>(io.rawInteger(-6, 6));
            if (ExactOrient2(A, B, C) < 0)
            {
                std::swap(B, C);
            }
            acute = ExactOrient2(A, B, C) > 0 && ExactAcute(A, B, C);
        }
        if (!acute)
        {
            continue;
        }
        pts = { A, B, C };
        for (int32_t e = 0; e < extra; ++e)
        {
            for (int32_t k = 0; k < 32; ++k)
            {
                Vector2<double> D{ static_cast<double>(io.rawInteger(-6, 6)),
                    static_cast<double>(io.rawInteger(-6, 6)) };
                if (ExactInCircle(A, B, C, D) > 0)
                {
                    pts.push_back(D);
                    break;
                }
            }
        }
        accepted = Invariant2(pts, res);
    }
    if (!accepted)
    {
        pts = { Vector2<double>{ -4.0, 0.0 }, Vector2<double>{ 4.0, 0.0 },
            Vector2<double>{ 0.0, 3.0 } };
        MinimumAreaCircle2<double, double> q{};
        res = RunMAC(q, pts);
    }

    io.given(static_cast<double>(pts.size()));
    for (auto const& p : pts)
    {
        io.givenVec<2>(p);
    }

    io.outBool(res.ok);
    io.outReal(res.cx);
    io.outReal(res.cy);
    io.outReal(res.radius);
    io.outInt(res.numSupport);
    for (int32_t i = 0; i < res.numSupport; ++i)
    {
        io.outInt(res.support[i]);
    }
}
