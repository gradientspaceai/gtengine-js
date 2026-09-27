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
//    which the support circles/spheres are built, and ExactCircle3 /
//    ExactSphere4 are not symmetric in their arguments, so the result of a
//    general input depends on it. The cases make upstream run the port's
//    permutation: a private engine runs in lockstep with the query's mDRE
//    and the query is called until its shuffle equals the permutation the
//    port's first compute() applies (see PortOrderResult). Every input with
//    at most 6 unique points is then comparable bit for bit.
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
#include <cstdlib>
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
    //   mode 3: the convex hull of random dyadic points (multiples of 2^-10)
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
            else if (mode == 3)
            {
                // Dyadic reals, multiples of 2^-10: sums and halves of two
                // coordinates are exact, so edge midpoints lie exactly on
                // their edge (the .collinear case).
                pts[i][0] = static_cast<double>(io.rawInteger(-6144, 6144)) / 1024.0;
                pts[i][1] = static_cast<double>(io.rawInteger(-6144, 6144)) / 1024.0;
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
    // Modes 0, 1 and 3: every coordinate is an integer or a multiple of
    // 2^-10, so the midpoints below are exact. (Mode 2's arbitrary reals
    // would round a + b, leaving the "midpoint" a hair off the edge, and
    // when it falls inside, the polygon is not convex any more.)
    int32_t const modes[3] = { 0, 1, 3 };
    int32_t mode = modes[io.index() % 3];
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

    // Subdivide a random subset of the edges at their midpoints. With these
    // coordinates the sum and the halving are exact, so the inserted vertex
    // is exactly on the segment and DotPerp of the two edges it creates is
    // exactly zero; the exact orientation test below confirms it.
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
            LogAssert(ExactOrient2(poly[i0], mid, poly[i1]) == 0,
                "The midpoint must be exactly on the edge.");
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

// Throw parity: CreatePolygon keeps only the vertices at which the polygon
// turns, so 2 to 6 exactly collinear points (in any order, repeats allowed)
// leave fewer than 3 and ComputeAntipodes' LogAssert fires. One point, or
// none, is excluded: upstream then reads vertices[1] (or vertices.back())
// out of bounds before the assert, which the port replaces by an up-front
// assert (issue #286, "reads vertices.back() and vertices[1] before the size
// assert").
ORACLE_CASE("RotatingCalipers.computeAntipodes.collinearThrows")
{
    int32_t n = io.integer(2, 6);
    int32_t bx = io.rawInteger(-3, 3);
    int32_t by = io.rawInteger(-3, 3);
    int32_t dx = io.rawInteger(1, 3);
    int32_t dy = io.rawInteger(-3, 3);
    std::vector<Vector2<double>> pts{};
    for (int32_t i = 0; i < n; ++i)
    {
        double k = static_cast<double>(io.rawInteger(-2, 2));
        pts.push_back(Vector2<double>{ static_cast<double>(bx) + k * static_cast<double>(dx),
            static_cast<double>(by) + k * static_cast<double>(dy) });
    }
    for (auto const& p : pts)
    {
        io.givenVec<2>(p);
    }
    EmitAntipodes(io, pts);
}

// ---- MinimumAreaCircle2 / MinimumVolumeSphere3 ---------------------------

namespace
{
    struct Result2
    {
        bool ok;
        double cx, cy, radius;
        int32_t numSupport;
        std::array<int32_t, 3> support;
    };

    // The support indices are emitted in upstream's own order: the case
    // reproduces the port's permutation exactly (see below), so the update
    // functions fill mSupport in the same order on both sides.
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
        return r;
    }

    struct Result3
    {
        bool ok;
        double cx, cy, cz, radius;
        int32_t numSupport;
        std::array<int32_t, 4> support;
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
        return r;
    }

    // THE SHUFFLE, EXACTLY. operator() permutes the unique points with
    // std::shuffle(permuted, mDRE), and mDRE is a default-constructed
    // std::default_random_engine - std::mt19937 on MSVC. The port uses a
    // minstd_rand0-style Lehmer generator with its own Fisher-Yates loop, so
    // the two permutations differ, and the result depends on the permutation
    // (ExactCircle3/ExactSphere4 are not symmetric in their arguments, and
    // the support set of a cocircular input depends on the visiting order).
    // Instead of restricting the cases to permutation-independent inputs,
    // the C++ side makes upstream run THE PORT'S permutation:
    //
    //  * PortPermutation reproduces the permutation the port's first
    //    compute() applies to m unique points (the port's shuffle with its
    //    own generator, seed 1; the replay constructs a fresh query per
    //    record).
    //  * A private std::default_random_engine runs in lockstep with the
    //    query's mDRE: both start in the default state, mDRE is used for
    //    nothing but the one std::shuffle of m elements per operator() call,
    //    and std::shuffle's engine consumption depends only on the length,
    //    so shuffling a local index array once per call yields exactly the
    //    permutation that call applied.
    //
    // The query is called until the tracked permutation equals the port's,
    // and THAT call's result is the one recorded. Before the shuffle both
    // sides hold the same array (sorted unique indices; MSVC's std::sort is
    // an insertion sort, hence stable, below 32 elements, and the port sorts
    // with an explicit index tie-break), so from the shuffle on both run the
    // same computation on the same point order. The expected number of calls
    // is m!, so the unique count is capped at 6 (720) and the loop at
    // gPermutationCap calls (the probability of not meeting a given
    // permutation of 6 is below e^-27).
    int32_t const gMaxUnique = 6;
    int32_t const gPermutationCap = 20000;

    std::vector<int32_t> PortPermutation(int32_t m)
    {
        std::vector<int32_t> p(static_cast<size_t>(m));
        std::iota(p.begin(), p.end(), 0);
        uint64_t state = 1;
        for (int32_t i = m - 1; i > 0; --i)
        {
            state = (16807ull * state) % 2147483647ull;
            int32_t j = static_cast<int32_t>(state % static_cast<uint64_t>(i + 1));
            std::swap(p[static_cast<size_t>(i)], p[static_cast<size_t>(j)]);
        }
        return p;
    }

    // The number of unique points, which is the length of the array
    // operator() shuffles (upstream's sort/unique with the Vector operators).
    template <typename V>
    int32_t UniqueCount(std::vector<V> const& pts)
    {
        std::vector<int32_t> idx(pts.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&pts](int32_t a, int32_t b)
        {
            return pts[static_cast<size_t>(a)] < pts[static_cast<size_t>(b)];
        });
        auto end = std::unique(idx.begin(), idx.end(), [&pts](int32_t a, int32_t b)
        {
            return pts[static_cast<size_t>(a)] == pts[static_cast<size_t>(b)];
        });
        return static_cast<int32_t>(end - idx.begin());
    }

    // Upstream's result under the port's permutation. Returns false when the
    // unique count exceeds gMaxUnique or the permutation did not turn up.
    template <typename Query, typename Result, typename V, typename Run>
    bool PortOrderResult(std::vector<V> const& pts, Result& res, Run run)
    {
        int32_t m = UniqueCount(pts);
        if (m > gMaxUnique)
        {
            return false;
        }
        std::vector<int32_t> const target = PortPermutation(m);
        Query q{};
        std::default_random_engine dre{};
        for (int32_t k = 0; k < gPermutationCap; ++k)
        {
            std::vector<int32_t> p(static_cast<size_t>(m));
            std::iota(p.begin(), p.end(), 0);
            std::shuffle(p.begin(), p.end(), dre);
            Result r = run(q, pts);
            if (p == target)
            {
                res = r;
                return true;
            }
        }
        return false;
    }

    // A result the main cases can compare. The trapped-failure branch calls
    // GetContainer(numPoints, points) with the unique count and the full
    // array (issue #286, fixed in the port); it is sound exactly when there
    // is no duplicate, because then the prefix is the whole array.
    bool PortOrderMAC(std::vector<Vector2<double>> const& pts, Result2& res)
    {
        bool found = PortOrderResult<MinimumAreaCircle2<double, double>>(pts, res, RunMAC);
        return found && (res.ok || UniqueCount(pts) == static_cast<int32_t>(pts.size()));
    }

    bool PortOrderMVS(std::vector<Vector3<double>> const& pts, Result3& res)
    {
        bool found = PortOrderResult<MinimumVolumeSphere3<double, double>>(pts, res, RunMVS);
        return found && (res.ok || UniqueCount(pts) == static_cast<int32_t>(pts.size()));
    }

    void EmitMAC(oracle::Ctx& io, Result2 const& res)
    {
        io.outBool(res.ok);
        io.outReal(res.cx);
        io.outReal(res.cy);
        io.outReal(res.radius);
        io.outInt(res.numSupport);
        for (int32_t i = 0; i < res.numSupport; ++i)
        {
            io.outInt(res.support[static_cast<size_t>(i)]);
        }
    }

    void EmitMVS(oracle::Ctx& io, Result3 const& res)
    {
        io.outBool(res.ok);
        io.outReal(res.cx);
        io.outReal(res.cy);
        io.outReal(res.cz);
        io.outReal(res.radius);
        io.outInt(res.numSupport);
        for (int32_t i = 0; i < res.numSupport; ++i)
        {
            io.outInt(res.support[static_cast<size_t>(i)]);
        }
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

    // mode 0: uniform doubles; mode 1: lattice [-4,4]; mode 2: dense lattice
    // [-2,2]; mode 3: cospherical lattice points about a lattice center.
    void RawPoints3(oracle::Ctx& io, int32_t mode, size_t n,
        std::vector<Vector3<double>>& pts)
    {
        pts.resize(n);
        if (mode == 3)
        {
            // The 30 lattice points of the sphere x^2 + y^2 + z^2 = 9 about
            // a lattice center, plus the center: cospherical quadruples and
            // quintuples, where Contains ties and the four-point support is
            // not unique.
            static std::vector<Vector3<double>> const gSphere9 = []()
            {
                std::vector<Vector3<double>> s{};
                for (int32_t x = -3; x <= 3; ++x)
                {
                    for (int32_t y = -3; y <= 3; ++y)
                    {
                        for (int32_t z = -3; z <= 3; ++z)
                        {
                            if (x * x + y * y + z * z == 9)
                            {
                                s.push_back(Vector3<double>{ static_cast<double>(x),
                                    static_cast<double>(y), static_cast<double>(z) });
                            }
                        }
                    }
                }
                return s;
            }();
            double cx = static_cast<double>(io.rawInteger(-2, 2));
            double cy = static_cast<double>(io.rawInteger(-2, 2));
            double cz = static_cast<double>(io.rawInteger(-2, 2));
            int32_t const count = static_cast<int32_t>(gSphere9.size());
            for (size_t i = 0; i < n; ++i)
            {
                int32_t k = io.rawInteger(0, count);
                Vector3<double> offset{ 0.0, 0.0, 0.0 };
                if (k < count)
                {
                    offset = gSphere9[static_cast<size_t>(k)];
                }
                pts[i] = { cx + offset[0], cy + offset[1], cz + offset[2] };
            }
            return;
        }

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
// GetNumSupport() and GetSupport(), under the port's permutation (see
// PortOrderResult). Every input is comparable, including cocircular and
// tied ones whose result depends on the permutation; the only rejections
// are more than 6 unique points (the lockstep cost) and a trapped failure
// with a duplicate point (the #286 prefix defect, owned by the .deviation
// case). The support indices are emitted in upstream's order.
//
// The points are at most 8 (6 unique), well below the 32 elements at which
// MSVC's std::sort stops being an insertion sort, so upstream's choice of
// representative among equal points - the first in index order - agrees
// with the port's explicit index tie-break.
ORACLE_CASE("MinimumAreaCircle2.compute")
{
    int32_t mode = io.index() % 4;
    std::vector<Vector2<double>> pts{};
    Result2 res{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        // Uniform points are unique, so they are capped at 6; the lattice
        // modes may draw up to 8 and rely on duplicates.
        int32_t n = io.rawInteger(1, mode == 0 ? gMaxUnique : 8);
        RawPoints2(io, mode, static_cast<size_t>(n), pts);
        accepted = PortOrderMAC(pts, res);
    }
    if (!accepted)
    {
        pts = { Vector2<double>{ -2.0, 4.0 }, Vector2<double>{ 0.0, 0.0 },
            Vector2<double>{ 1.0, 1.0 } };
        accepted = PortOrderMAC(pts, res);
        LogAssert(accepted, "The fallback must be comparable.");
    }

    io.given(static_cast<double>(pts.size()));
    for (auto const& p : pts)
    {
        io.givenVec<2>(p);
    }
    EmitMAC(io, res);
}

namespace
{
    // Deviation input for issue #286: a catalogued point set whose query
    // traps under the port's permutation, presented in a random order with
    // 1..3 copies of one of its points prefixed. The set of distinct points,
    // hence the sorted unique array, the shuffle and the control flow, do not
    // depend on the input order or on the duplicates, so both sides reach
    // the trapped branch at the same point; the duplicates make the unique
    // count smaller than the array length, so upstream's GetContainer prefix
    // leaves out genuine points while the port bounds them all.
    template <int32_t N>
    std::vector<Vector<N, double>> TrappedInput(oracle::Ctx& io,
        std::vector<double> const& raw)
    {
        size_t const m = raw.size() / static_cast<size_t>(N);
        std::vector<Vector<N, double>> base(m);
        for (size_t i = 0; i < m; ++i)
        {
            for (int32_t j = 0; j < N; ++j)
            {
                base[i][j] = raw[static_cast<size_t>(N) * i + static_cast<size_t>(j)];
            }
        }
        for (size_t i = m - 1; i > 0; --i)
        {
            size_t j = static_cast<size_t>(io.rawInteger(0, static_cast<int32_t>(i)));
            std::swap(base[i], base[j]);
        }
        int32_t numDup = io.rawInteger(1, 3);
        size_t dupOf = static_cast<size_t>(io.rawInteger(0, static_cast<int32_t>(m) - 1));
        std::vector<Vector<N, double>> pts(static_cast<size_t>(numDup), base[dupOf]);
        pts.insert(pts.end(), base.begin(), base.end());
        return pts;
    }
}

// Deliberate port fix of issue #286: in the trapped-failure branch upstream
// calls GetContainer(numPoints, points, minimal) after numPoints has been
// overwritten with the *unique* point count, so the fallback circle bounds
// only a prefix of the input array. The port passes the whole array.
//
// The trapped branch is reached by about 1 in 13000 random lattice sets of
// 4-6 distinct points in [-4,4]^2 (1 in 7000 sets of 5-6 points in
// [-2,2]^3 for the sphere), far too rarely for a rejection loop, so the case draws from
// a catalogue of lattice sets on which the port's compute() traps (found
// with the port; the C++ side confirms that upstream traps as well under the
// same permutation: every record's first output is false).
ORACLE_CASE("MinimumAreaCircle2.compute.deviation.trappedFallback")
{
    static std::vector<std::vector<double>> const gTrapped2
    {
        { -1, 2, 0, -4, 4, -2, 1, -2, 2, 2, -3, -2 },
        { -4, 1, -2, 4, 2, 3, 1, -2, -1, 4 },
        { -1, -4, 3, 1, 3, -2, -1, 3, -1, 1 },
        { 4, 1, -1, -3, -3, 1, 2, -3, 0, 2 },
        { 3, 3, -3, 3, -4, -2, 0, -4, 4, -2 },
        { -4, 0, -2, 4, 3, 0, -1, 4, -3, 1, 1, 4 },
        { 2, -4, -3, -1, 1, 2, -1, 2, 3, 1 },
        { -4, -1, -1, -1, 1, 2, 2, 0, -2, -4, 2, -3 },
        { 3, 4, -1, 2, -1, -3, -1, 4, -1, 0, -3, 1 },
        { -1, 4, 3, 2, -4, 2, -4, 3, -4, -2 },
        { -3, 3, 3, 3, 4, -2, -3, 1, -4, -2 },
        { 1, -1, -4, 3, 1, 3, 3, 3, -2, -1, -1, 4 },
        { -3, 2, -4, -3, -4, -1, 0, -1, 0, 2, -4, 4 },
        { 0, 4, 2, 0, -2, 1, 4, 3, 3, -2 },
        { -1, 3, 3, 1, 4, -1, 3, -2, 2, 3, -3, -1 },
        { 2, -3, 4, 1, -1, -3, -3, 1 }
    };

    int32_t which = io.rawInteger(0, static_cast<int32_t>(gTrapped2.size()) - 1);
    std::vector<Vector2<double>> pts = TrappedInput<2>(io, gTrapped2[static_cast<size_t>(which)]);
    io.given(static_cast<double>(pts.size()));
    for (auto const& p : pts)
    {
        io.givenVec<2>(p);
    }

    Result2 res{};
    bool found = PortOrderResult<MinimumAreaCircle2<double, double>>(pts, res, RunMAC);
    LogAssert(found, "The catalogue sets have at most 6 unique points.");
    EmitMAC(io, res);
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
// comparison policy as the 2D query above: upstream runs the port's
// permutation, and every input with at most 6 unique points is comparable
// except a trapped failure with a duplicate point (#286, the .deviation
// case).
ORACLE_CASE("MinimumVolumeSphere3.compute")
{
    int32_t mode = io.index() % 4;
    std::vector<Vector3<double>> pts{};
    Result3 res{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        int32_t n = io.rawInteger(1, mode == 0 ? gMaxUnique : 8);
        RawPoints3(io, mode, static_cast<size_t>(n), pts);
        accepted = PortOrderMVS(pts, res);
    }
    if (!accepted)
    {
        pts = { Vector3<double>{ -1.0, 1.0, -1.0 }, Vector3<double>{ 0.0, 0.0, 0.0 },
            Vector3<double>{ 1.0, 1.0, 1.0 }, Vector3<double>{ 2.0, 4.0, 8.0 } };
        accepted = PortOrderMVS(pts, res);
        LogAssert(accepted, "The fallback must be comparable.");
    }

    io.given(static_cast<double>(pts.size()));
    for (auto const& p : pts)
    {
        io.givenVec<3>(p);
    }
    EmitMVS(io, res);
}

// The 3D twin of the trapped-fallback deviation; see the 2D case.
ORACLE_CASE("MinimumVolumeSphere3.compute.deviation.trappedFallback")
{
    static std::vector<std::vector<double>> const gTrapped3
    {
        { -1, 1, -1, 0, -1, -1, -1, 0, 1, -2, -2, 0, 0, -1, 1, -1, -2, -2 },
        { 0, -1, 1, 0, -2, -1, 2, -1, -1, -1, 2, 2, 2, 2, 1 },
        { -1, 1, -1, -1, -1, -2, 0, 1, -2, 0, -1, -1, 1, -2, 0, 2, 0, 0 },
        { -1, 2, 2, -1, 0, 1, 0, 2, 0, 1, -2, 0, 1, 1, 2, -2, -2, -1 },
        { -1, 0, 1, 1, -2, 0, -1, -1, 0, 2, 1, 0, 2, 0, -1 },
        { -2, 0, -1, 0, -1, 2, 1, -2, 0, 1, 0, 2, -1, 2, -1, -1, 1, -2 },
        { 2, 0, -1, -1, -1, -2, 2, 0, -2, 2, 2, 0, 0, -1, 1, -1, 2, 1 },
        { -1, 1, -1, 1, 2, 2, -1, 0, -1, -2, 2, 2, -2, 2, 1, 1, -1, 2 },
        { -2, 1, 1, 0, 0, -2, 2, 1, -1, 1, 1, 1, 2, -2, 0, 0, 2, 2 },
        { -2, 1, 1, 1, -2, 1, -1, 2, -2, -2, -2, -2, 1, 2, 0, -2, -2, -1 },
        { 0, -2, 2, -1, -2, -1, 0, 2, 0, -1, 1, -1, 0, 0, -2, -1, 1, 2 },
        { -2, -1, 2, -2, 1, 2, 1, -1, 2, 0, 0, -2, 1, 2, 0, -1, 0, -2 },
        { -2, 1, 1, 2, 0, 1, 2, -2, -1, -2, -2, -2, -1, 0, 0, -2, 0, -1 },
        { -2, -2, -1, 1, 0, 2, 2, 1, -1, 2, 0, -1, 1, 2, 0, -2, -1, -2 },
        { 1, -1, -2, 0, -2, -1, -2, 2, 1, -1, -2, 0, 1, 0, 0, 0, 2, 2 },
        { 2, 2, 1, 2, 2, -1, -1, -1, -1, -1, 2, 2, 2, 0, 1, 2, 0, -1 }
    };

    int32_t which = io.rawInteger(0, static_cast<int32_t>(gTrapped3.size()) - 1);
    std::vector<Vector3<double>> pts = TrappedInput<3>(io, gTrapped3[static_cast<size_t>(which)]);
    io.given(static_cast<double>(pts.size()));
    for (auto const& p : pts)
    {
        io.givenVec<3>(p);
    }

    Result3 res{};
    bool found = PortOrderResult<MinimumVolumeSphere3<double, double>>(pts, res, RunMVS);
    LogAssert(found, "The catalogue sets have at most 6 unique points.");
    EmitMVS(io, res);
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
    // Condition (b) is a guard for the preserved construction defect of
    // issue #290; on the polytopes this case draws it never fired in the
    // measurements of the group report. Sheared octahedra are NOT drawn:
    // there the MSVC build answers about 0.5% of random directions with a
    // vertex that is not extreme, and which directions those are changes
    // from one construction to the next (the tree depends on the pointer
    // hash order of unordered_set<Triangle*>, i.e. on heap addresses), so
    // the accepted directions - and the goldens - would not be reproducible.
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
// partitions of the same Gauss map (the port builds 15 nodes for the
// octahedron; upstream's count varies from one construction to the next,
// between 11 and 28 nodes for sheared octahedra, because it depends on heap
// addresses).
//
// Directions whose extreme vertex is not unique are rejected for the same
// reason: a tie is resolved by which leaf region the direction lands in.
ORACLE_CASE("ExtremalQuery3BSP.getExtremeVertices")
{
    int32_t which = io.integer(0, 2);
    BasePolytope const& base = BasePolytopeOf(which);
    std::vector<Vector3<double>> verts{};
    if (which == 1)
    {
        // The regular octahedron only; see SoundDirection.
        verts = base.vertices;
    }
    else
    {
        TransformPolytope(io, base, verts);
    }
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
    // [-5,5]^2, hence strictly inside. Every other record uses real bounds,
    // so the supervertex arithmetic of the constructor (xMin - dx,
    // xMin + 5 * dx, ...) is not exact and GetTriangulation's vertex list
    // discriminates its evaluation order.
    std::array<double, 4> DomainRectangle(oracle::Ctx& io)
    {
        if (io.index() % 2 == 0)
        {
            double xMin = io.lattice(-10, -6);
            double yMin = io.lattice(-10, -6);
            double xMax = io.lattice(6, 10);
            double yMax = io.lattice(6, 10);
            return { xMin, yMin, xMax, yMax };
        }
        double xMin = io.real(-10.0, -6.0);
        double yMin = io.real(-10.0, -6.0);
        double xMax = io.real(6.0, 10.0);
        double yMax = io.real(6.0, 10.0);
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

// Insert, Remove, then Insert again. Removals hit positions that are
// vertices, positions that are not vertices and positions already removed.
// An inserted point is strictly inside the rectangle, so it is never
// adjacent to a supervertex and its removal always takes
// RetriangulateInteriorRemovalPolygon; the boundary branch is reached by
// FinalizeTriangulation's removal of the rectangle corners (the
// .finalizeTriangulation case). The return value of Remove is the vertex index,
// or invalid when the position is not a vertex. The re-insertions put half
// of their points back at a previously inserted (possibly removed) position,
// which receives a fresh index because mVertices never shrinks.
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
    int32_t numAfter = io.integer(0, 3);
    for (int32_t k = 0; k < numAfter; ++k)
    {
        Vector2<double> p{ 0.0, 0.0 };
        if (io.rawInteger(0, 1) != 0)
        {
            size_t j = static_cast<size_t>(io.rawInteger(0, n - 1));
            p = inserted[j];
        }
        else
        {
            p = RawInsidePoint(io, mode);
        }
        io.givenVec<2>(p);
        size_t index = del.Insert(p);
        io.outReal(AsIndex(index));
    }
    io.outInt(del.GetNumVertices());
    io.outInt(del.GetNumTriangles());
    EmitTriangles(io, del);
    EmitHull(io, del);
}

// GetContainingTriangle (with its whole SearchInfo), GetTriangle,
// GetAdjacent and GetTriangulation.
//
// The walk starts at info.initialTriangle, and the default start (triangle 0)
// is hash-table numbering. As in v09's Delaunay2 cases, the start is chosen
// in the canonical numbering instead (the triangles sorted by their stored
// vertex tuple): the record holds the canonical rank of the start, each side
// maps it to its own index, and the path, finalTriangle and initialTriangle
// are emitted translated back to canonical ranks. The walk only follows
// adjacency, which is order independent, so the whole search is comparable,
// including query points on a shared edge or at a vertex (the walk stops at
// whichever triangle it reaches first, which is now the same triangle on
// both sides). Query points: uniform, lattice, an existing vertex, and the
// exact midpoint of a triangle edge.
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

    // GetTriangulation: all vertices (the three supervertices and the four
    // rectangle corners first) and all triangles of the graph, the
    // supervertex triangles included. The triangle keys are TriangleKey<true>
    // tuples and the container is a hash map, so the list is sorted
    // lexicographically on both sides.
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
    auto asRank = [&rank](size_t t)
    {
        return t == gInvalid ? -1.0 : static_cast<double>(rank[t]);
    };

    int32_t const numQueries = 4;
    for (int32_t k = 0; k < numQueries; ++k)
    {
        Vector2<double> q{ 0.0, 0.0 };
        int32_t qmode = io.rawInteger(0, 3);
        if (qmode == 0)
        {
            // Often outside the rectangle, hence outside the hull.
            q[0] = io.raw(-11.0, 11.0);
            q[1] = io.raw(-11.0, 11.0);
        }
        else if (qmode == 1)
        {
            q[0] = static_cast<double>(io.rawInteger(-6, 6));
            q[1] = static_cast<double>(io.rawInteger(-6, 6));
        }
        else
        {
            auto const& t = tris[static_cast<size_t>(io.rawInteger(0,
                static_cast<int32_t>(count) - 1))];
            size_t j0 = static_cast<size_t>(io.rawInteger(0, 2));
            Vector2<double> const& a = verts[t[j0]];
            Vector2<double> const& b = verts[t[(j0 + 1) % 3]];
            q = (qmode == 2 ? a : Vector2<double>{ 0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]) });
        }
        io.givenVec<2>(q);
        int32_t startRank = io.integer(0, static_cast<int32_t>(count) - 1);

        ID2::SearchInfo info{};
        info.initialTriangle = order[static_cast<size_t>(startRank)];
        size_t t = del.GetContainingTriangle(q, info);
        io.outReal(asRank(t));
        io.outReal(asRank(info.initialTriangle));
        io.outReal(asRank(info.finalTriangle));
        io.outInt(info.finalV[0]);
        io.outInt(info.finalV[1]);
        io.outInt(info.finalV[2]);
        io.outInt(info.numPath);
        for (size_t i = 0; i < info.numPath; ++i)
        {
            io.outReal(asRank(info.path[i]));
        }
        if (t != gInvalid)
        {
            std::array<size_t, 3> triangle{};
            bool valid = del.GetTriangle(t, triangle);
            io.outBool(valid);
            io.outInt(triangle[0]);
            io.outInt(triangle[1]);
            io.outInt(triangle[2]);
            std::array<size_t, 3> adjacent{};
            bool validAdj = del.GetAdjacent(t, adjacent);
            io.outBool(validAdj);
            for (size_t j = 0; j < 3; ++j)
            {
                io.outReal(asRank(adjacent[j]));
            }
        }
    }

    // Out-of-range accessors return false.
    std::array<size_t, 3> triangle{}, adjacent{};
    bool validTriangle = del.GetTriangle(count, triangle);
    io.outBool(validTriangle);
    bool validAdjacent = del.GetAdjacent(count, adjacent);
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

// Throw parity for the domain preconditions: the constructor asserts
// xMin < xMax and yMin < yMax, Insert asserts that the position is strictly
// inside the rectangle, and Remove asserts the same before
// FinalizeTriangulation (and returns invalid afterwards instead).
//   mode 0: an empty or inverted rectangle (the constructor throws)
//   mode 1: Insert of a position on or outside the rectangle (throws)
//   mode 2: Remove of a position on or outside the rectangle (throws)
//   mode 3: Remove of the same kind of position after
//           FinalizeTriangulation (returns invalid, no throw)
ORACLE_CASE("IncrementalDelaunay2.domainAsserts")
{
    int32_t mode = io.index() % 4;
    double xMin = io.lattice(-4, 0);
    double yMin = io.lattice(-4, 0);
    double xMax = 0.0, yMax = 0.0;
    if (mode == 0)
    {
        // One of the two extents is empty or negative.
        bool flatX = io.boolean();
        double lo = io.lattice(-4, 0);
        xMax = flatX ? xMin + lo : xMin + 4.0;
        yMax = flatX ? yMin + 4.0 : yMin + lo;
    }
    else
    {
        xMax = xMin + 8.0;
        yMax = yMin + 8.0;
    }
    ID2 del(xMin, yMin, xMax, yMax);

    Vector2<double> inside{ xMin + 3.0, yMin + 5.0 };
    size_t first = del.Insert(inside);

    // A position on the boundary (a coordinate equal to a bound) or just
    // outside it.
    double t = io.lattice(0, 8);
    int32_t side = io.integer(0, 3);
    double off = io.lattice(0, 1);
    Vector2<double> p{ 0.0, 0.0 };
    if (side == 0) { p = { xMin - off, yMin + t }; }
    else if (side == 1) { p = { xMax + off, yMin + t }; }
    else if (side == 2) { p = { xMin + t, yMin - off }; }
    else { p = { xMin + t, yMax + off }; }

    // A throw record has no outputs at all, so nothing is emitted before
    // the call that may throw.
    if (mode == 1)
    {
        size_t index = del.Insert(p);
        io.outReal(AsIndex(first));
        io.outReal(AsIndex(index));
    }
    else if (mode == 2)
    {
        size_t index = del.Remove(p);
        io.outReal(AsIndex(first));
        io.outReal(AsIndex(index));
    }
    else if (mode == 3)
    {
        io.outReal(AsIndex(first));
        bool finalized = del.FinalizeTriangulation();
        io.outBool(finalized);
        size_t index = del.Remove(p);
        io.outReal(AsIndex(index));
    }
}

// Deliberate port fix of issue #290: GetHull walks the edge map with an
// unbounded while (vNext != vStart) loop that writes into hull[], which is
// sized to the number of edges. For a finalized triangulation whose input
// points are collinear there are no Delaunay triangles at all, yet the
// triangles sharing a supervertex still contribute edges, and those edges
// do not form one cycle through all of them. Two outcomes:
//  * the walk returns to its start early (three or more collinear points
//    with the start at one end): upstream returns hull[] with numEdges
//    entries, the unwritten tail left at the zeros of the resize, i.e. a
//    "polygon" containing supervertex 0;
//  * the walk enters a 2-cycle that avoids the start: upstream loops forever
//    and writes past the end of hull[] - undefined behaviour that cannot be
//    executed inside the generator.
// The port throws in both situations (the degenerate triangulation has no
// convex hull polygon). The case runs a verbatim copy of upstream's edge
// collection (GetGraph() exposes everything it reads) and upstream's walk as
// a probe, calls the real GetHull when that is safe and emits the replica's
// prefix when it is not; see below. Records with one or two distinct points
// are the residue on which both sides agree: no edge at all (upstream would
// dereference edges.begin() of an empty map; the replica and the port both
// report an empty hull) or a 2-cycle through both edges (upstream is sound).
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

    // Keep every point in [-5,5]^2, strictly inside every domain rectangle
    // (a point on the rectangle would make Insert throw and the record would
    // test nothing).
    int32_t const reach = std::max(std::abs(bx), std::abs(by));
    int32_t const step = std::max(std::abs(dx), std::abs(dy));
    int32_t const kMax = (5 - reach) / step;

    ID2 del(rect[0], rect[1], rect[2], rect[3]);
    for (int32_t i = 0; i < n; ++i)
    {
        double k = static_cast<double>(io.rawInteger(-kMax, kMax));
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
        io.outBool(false);
        io.outInt(0);
        return;
    }

    // Upstream's walk as a probe, bounded by the size of hull[]. Upstream
    // writes hull[++i] for every step, and hull[] has numEdges entries, so
    // the call is safe exactly when the walk returns to vStart within
    // numEdges - 1 steps (or its LogAssert fires first, which is an ordinary
    // exception). Then the REAL upstream GetHull is called and its output is
    // emitted verbatim: the numEdges entries of hull[], where the entries
    // after an early return to vStart keep the zeros of the resize (a
    // supervertex index, i.e. garbage). Otherwise upstream overruns hull[]
    // and never terminates, which cannot be executed, and the record emits
    // the replica's prefix: what upstream writes up to the end of hull[].
    size_t const numEdges = edges.size();
    auto eIter = edges.begin();
    size_t vStart = eIter->first;
    size_t vNext = eIter->second;
    std::vector<size_t> prefix{ vStart };
    bool overrun = false;
    while (vNext != vStart)
    {
        if (prefix.size() == numEdges)
        {
            overrun = true;
            break;
        }
        prefix.push_back(vNext);
        auto it = edges.find(vNext);
        LogAssert(it != edges.end(), "Expecting to find a hull edge.");
        vNext = it->second;
    }

    io.outBool(overrun);
    if (overrun)
    {
        io.outInt(prefix.size());
        for (auto v : prefix)
        {
            io.outInt(v);
        }
    }
    else
    {
        std::vector<size_t> hull{};
        del.GetHull(hull);
        io.outInt(hull.size());
        for (auto v : hull)
        {
            io.outInt(v);
        }
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

// Aimed at the three-point support: an acute lattice triangle (its
// minimum-area circle is its circumcircle) plus extra points that the exact
// in-circle predicate places strictly inside. ExactCircle3's center and
// squared radius are then the emitted result. Upstream runs the port's
// permutation, as in the general case.
ORACLE_CASE("MinimumAreaCircle2.compute.circumcircle")
{
    // At most 3 extra points: 6 unique points is the lockstep cap.
    int32_t extra = io.integer(0, 3);
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
        accepted = PortOrderMAC(pts, res);
    }
    if (!accepted)
    {
        // An acute lattice triangle.
        pts = { Vector2<double>{ -4.0, 0.0 }, Vector2<double>{ 4.0, 0.0 },
            Vector2<double>{ 0.0, 5.0 } };
        accepted = PortOrderMAC(pts, res);
        LogAssert(accepted, "The fallback must be comparable.");
    }

    io.given(static_cast<double>(pts.size()));
    for (auto const& p : pts)
    {
        io.givenVec<2>(p);
    }
    EmitMAC(io, res);
}
