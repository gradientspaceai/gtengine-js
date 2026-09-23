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
