// Verify group 9 (computational geometry): differential cases for
// ConvexHull2.h, Delaunay2.h and Delaunay3.h.
//
// The three headers the port implements are the *replacement* class
// templates ConvexHull2<T>, Delaunay2<T> and Delaunay3<T> (the deprecated
// Delaunay2<InputType, ComputeType> / Delaunay3<InputType, ComputeType>
// specializations that route through PrimalQuery2/3 are not ported; see
// porting-status.json and the header comments of src/Delaunay2.ts). The C++
// side therefore instantiates Delaunay2<double> / Delaunay3<double> /
// ConvexHull2<double> with their own default rational types, which is the
// SWInterval fast path with a BSNumber exact fallback the port reproduces.
//
// Everything compared here is either a combinatorial output (indices,
// adjacencies, dimensions, duplicate maps, predicate signs) or a value
// computed with + - * / sqrt only (the dimension-1 line direction and the
// dimension-2 plane normal, through Normalize / UnitCross). No libm call is
// on any compared path, so every case is declared { exact: true }.
//
// CANONICALIZATION. ETManifoldMesh and TSManifoldMesh store their features in
// std::unordered_map, so upstream enumerates triangles / tetrahedra in
// hash-table order while the port enumerates them in sorted-key order. The
// triangle (tetrahedron) *set* and each feature's stored vertex tuple are
// order independent; only the numbering is not. Every case therefore sorts
// the features by their stored vertex tuple (which ETManifoldMesh::Insert
// keeps in the rotation the caller passed, so the tuple is a total key) and
// emits them in that order, remapping the adjacency indices through the same
// permutation. Hull edge lists (unordered by upstream's own comment) are
// sorted lexicographically. Everything else - dimension, counts, duplicate
// maps, search paths, per-vertex data - is emitted verbatim.
//
// Generator modes are selected by io.index() so that the record layout is
// fixed and the replay never has to know the mode.
//
// One io draw per C++ statement: MSVC evaluates function arguments right to
// left, so every generated value goes into a named local before it is used.
#define ORACLE_FAMILY "v09-compgeom"
#include "Oracle.h"

#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/ConvexHull2.h>
#include <Mathematics/Delaunay2.h>
#include <Mathematics/Delaunay3.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <exception>
#include <numeric>
#include <vector>

using namespace gte;

namespace
{
    // ---- exact predicates ------------------------------------------------
    // BSNumber<UIntegerAP32> grows as needed and the expressions below use
    // only + - *, so these signs are exact for any double input. They are
    // used only by the generator probes and by the independent reference
    // checks; nothing they compute reaches an output.
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

    // Sign of Dot(Cross(B - A, C - A), D - A): +1 when <A,B,C,D> is
    // positively oriented.
    int32_t ExactOrient3(Vector3<double> const& A, Vector3<double> const& B,
        Vector3<double> const& C, Vector3<double> const& D)
    {
        std::array<Exact, 3> u{ Exact(B[0]) - Exact(A[0]), Exact(B[1]) - Exact(A[1]),
            Exact(B[2]) - Exact(A[2]) };
        std::array<Exact, 3> v{ Exact(C[0]) - Exact(A[0]), Exact(C[1]) - Exact(A[1]),
            Exact(C[2]) - Exact(A[2]) };
        std::array<Exact, 3> w{ Exact(D[0]) - Exact(A[0]), Exact(D[1]) - Exact(A[1]),
            Exact(D[2]) - Exact(A[2]) };
        Exact c0 = u[1] * v[2] - u[2] * v[1];
        Exact c1 = u[2] * v[0] - u[0] * v[2];
        Exact c2 = u[0] * v[1] - u[1] * v[0];
        Exact det = c0 * w[0] + c1 * w[1] + c2 * w[2];
        return det.GetSign();
    }

    // The exact intrinsic dimension of a point set, which is what the port
    // computes with its own ToLine / ToPlane predicates.
    size_t ExactDimension2(std::vector<Vector2<double>> const& pts)
    {
        size_t n = pts.size();
        size_t i1 = n;
        for (size_t i = 1; i < n; ++i)
        {
            if (!(pts[i] == pts[0]))
            {
                i1 = i;
                break;
            }
        }
        if (i1 == n)
        {
            return 0;
        }
        for (size_t i = 1; i < n; ++i)
        {
            if (i != i1 && ExactOrient2(pts[0], pts[i1], pts[i]) != 0)
            {
                return 2;
            }
        }
        return 1;
    }

    size_t ExactDimension3(std::vector<Vector3<double>> const& pts)
    {
        size_t n = pts.size();
        size_t i1 = n;
        for (size_t i = 1; i < n; ++i)
        {
            if (!(pts[i] == pts[0]))
            {
                i1 = i;
                break;
            }
        }
        if (i1 == n)
        {
            return 0;
        }
        size_t i2 = n;
        for (size_t i = 1; i < n; ++i)
        {
            if (i == i1)
            {
                continue;
            }
            Vector3<double> u = pts[i1] - pts[0];
            Vector3<double> v = pts[i] - pts[0];
            // Exact test for Cross(u, v) != 0 using the 2x2 minors.
            Exact ux(u[0]), uy(u[1]), uz(u[2]), vx(v[0]), vy(v[1]), vz(v[2]);
            Exact m0 = uy * vz - uz * vy;
            Exact m1 = uz * vx - ux * vz;
            Exact m2 = ux * vy - uy * vx;
            if (m0.GetSign() != 0 || m1.GetSign() != 0 || m2.GetSign() != 0)
            {
                i2 = i;
                break;
            }
        }
        if (i2 == n)
        {
            return 1;
        }
        for (size_t i = 1; i < n; ++i)
        {
            if (i != i1 && i != i2
                && ExactOrient3(pts[0], pts[i1], pts[i2], pts[i]) != 0)
            {
                return 3;
            }
        }
        return 2;
    }

    // ---- generator probes (upstream's own control flow) ------------------
    // Delaunay2<T>::operator() builds IntrinsicsVector2<T> with a hardcoded
    // epsilon of 0 and trusts info.dimension, info.extreme[] and
    // info.extremeCCW, all of which are decided by floating-point distances
    // to a *normalized* frame. The port replaces that classification with its
    // own exact ToLine predicate (upstream finding #391). This probe is the
    // exact separator: it accepts exactly those inputs on which upstream's
    // floating-point classification agrees with exact arithmetic, so on an
    // accepted input the port's seed triangle is upstream's seed triangle and
    // the two algorithms are the same algorithm.
    bool Sound2(std::vector<Vector2<double>> const& pts)
    {
        IntrinsicsVector2<double> info(static_cast<int32_t>(pts.size()), pts.data(), 0.0);
        if (info.dimension != 2)
        {
            return false;
        }
        int32_t sign = ExactOrient2(pts[info.extreme[0]], pts[info.extreme[1]],
            pts[info.extreme[2]]);
        if (sign == 0)
        {
            return false;
        }
        return info.extremeCCW == (sign > 0);
    }

    // The 3D analogue. Delaunay3<T> additionally has a duplicate-detection
    // defect (#283): ProcessedVertex hashes and compares the 'location'
    // member, so a repeated coordinate triple is never recognized and Update
    // runs a second time for a point already in the mesh. The probe therefore
    // also rejects point sets with repeated coordinates; the deviation case
    // aims at exactly those.
    bool HasDuplicate3(std::vector<Vector3<double>> const& pts)
    {
        for (size_t i = 1; i < pts.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                if (pts[i] == pts[j])
                {
                    return true;
                }
            }
        }
        return false;
    }

    // The seed half of the probe: upstream's intrinsics report dimension 3
    // and its extremeCCW agrees with the exact orientation, so the seed
    // tetrahedron the port inserts is upstream's seed tetrahedron.
    bool SoundSeed3(std::vector<Vector3<double>> const& pts)
    {
        IntrinsicsVector3<double> info(static_cast<int32_t>(pts.size()), pts.data(), 0.0);
        if (info.dimension != 3)
        {
            return false;
        }
        int32_t sign = ExactOrient3(pts[info.extreme[0]], pts[info.extreme[1]],
            pts[info.extreme[2]], pts[info.extreme[3]]);
        if (sign == 0)
        {
            return false;
        }
        return info.extremeCCW == (sign > 0);
    }

    bool Sound3(std::vector<Vector3<double>> const& pts)
    {
        return !HasDuplicate3(pts) && SoundSeed3(pts);
    }

    // ---- point-set generators (raw, unrecorded) --------------------------
    // The 12 lattice points on the circle x^2 + y^2 = 25 and the 30 lattice
    // points on the sphere x^2 + y^2 + z^2 = 9. Subsets of these are exactly
    // cocircular / cospherical, so ToCircumcircle and ToCircumsphere see an
    // interval that straddles zero and the exact rational fallback decides
    // the sign. Adding the center of the circle (sphere) to the subset makes
    // the configuration non-degenerate again, so the mode reaches both.
    std::array<std::array<int32_t, 2>, 12> const gCircle25
    { {
        { 5, 0 }, { -5, 0 }, { 0, 5 }, { 0, -5 },
        { 3, 4 }, { 3, -4 }, { -3, 4 }, { -3, -4 },
        { 4, 3 }, { 4, -3 }, { -4, 3 }, { -4, -3 }
    } };

    std::array<std::array<int32_t, 3>, 30> const gSphere9
    { {
        { 3, 0, 0 }, { -3, 0, 0 }, { 0, 3, 0 }, { 0, -3, 0 }, { 0, 0, 3 }, { 0, 0, -3 },
        { 1, 2, 2 }, { 1, 2, -2 }, { 1, -2, 2 }, { 1, -2, -2 },
        { -1, 2, 2 }, { -1, 2, -2 }, { -1, -2, 2 }, { -1, -2, -2 },
        { 2, 1, 2 }, { 2, 1, -2 }, { 2, -1, 2 }, { 2, -1, -2 },
        { -2, 1, 2 }, { -2, 1, -2 }, { -2, -1, 2 }, { -2, -1, -2 },
        { 2, 2, 1 }, { 2, 2, -1 }, { 2, -2, 1 }, { 2, -2, -1 },
        { -2, 2, 1 }, { -2, 2, -1 }, { -2, -2, 1 }, { -2, -2, -1 }
    } };

    // mode 0: uniform doubles; the generic interval branch.
    // mode 1: lattice [-3,3]; exact arithmetic, collinear triples, ties.
    // mode 2: cocircular lattice about a lattice center, with the center and
    //         nearby lattice points mixed in.
    // mode 3: dense lattice [-2,2]; duplicates and degeneracies are common.
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
                pts[i][0] = static_cast<double>(io.rawInteger(-3, 3));
                pts[i][1] = static_cast<double>(io.rawInteger(-3, 3));
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
        if (mode == 2)
        {
            double cx = static_cast<double>(io.rawInteger(-2, 2));
            double cy = static_cast<double>(io.rawInteger(-2, 2));
            double cz = static_cast<double>(io.rawInteger(-2, 2));
            for (size_t i = 0; i < n; ++i)
            {
                int32_t k = io.rawInteger(0, 32);
                if (k < 30)
                {
                    pts[i][0] = cx + static_cast<double>(gSphere9[k][0]);
                    pts[i][1] = cy + static_cast<double>(gSphere9[k][1]);
                    pts[i][2] = cz + static_cast<double>(gSphere9[k][2]);
                }
                else if (k == 30)
                {
                    pts[i][0] = cx;
                    pts[i][1] = cy;
                    pts[i][2] = cz;
                }
                else
                {
                    pts[i][0] = cx + static_cast<double>(io.rawInteger(-2, 2));
                    pts[i][1] = cy + static_cast<double>(io.rawInteger(-2, 2));
                    pts[i][2] = cz + static_cast<double>(io.rawInteger(-2, 2));
                }
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
                    pts[i][j] = static_cast<double>(io.rawInteger(-3, 3));
                }
                else
                {
                    pts[i][j] = static_cast<double>(io.rawInteger(-2, 2));
                }
            }
        }
    }

    // The fallbacks used when a capped rejection loop expires. Points on a
    // parabola / on the moment curve are in general position, so upstream's
    // floating-point classification is never in doubt there.
    void Fallback2(size_t n, std::vector<Vector2<double>>& pts)
    {
        pts.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            double t = static_cast<double>(i) - 2.0;
            pts[i][0] = t;
            pts[i][1] = t * t;
        }
    }

    void Fallback3(size_t n, std::vector<Vector3<double>>& pts)
    {
        pts.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            double t = static_cast<double>(i) - 2.0;
            pts[i][0] = t;
            pts[i][1] = t * t;
            pts[i][2] = t * t * t;
        }
    }

    // Draw a point set on which upstream's own classification is sound. The
    // loop is capped and every coordinate is redrawn on each attempt; the
    // fallback is used when the cap expires (only the degenerate modes can
    // get there, and only rarely). The accepted set is recorded here, once.
    std::vector<Vector2<double>> SoundPoints2(oracle::Ctx& io, int32_t mode, size_t n)
    {
        std::vector<Vector2<double>> pts{};
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            RawPoints2(io, mode, n, pts);
            if (Sound2(pts))
            {
                break;
            }
            pts.clear();
        }
        if (pts.empty())
        {
            Fallback2(n, pts);
        }
        for (size_t i = 0; i < n; ++i)
        {
            io.givenVec<2>(pts[i]);
        }
        return pts;
    }

    std::vector<Vector3<double>> SoundPoints3(oracle::Ctx& io, int32_t mode, size_t n)
    {
        std::vector<Vector3<double>> pts{};
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            RawPoints3(io, mode, n, pts);
            if (Sound3(pts))
            {
                break;
            }
            pts.clear();
        }
        if (pts.empty())
        {
            Fallback3(n, pts);
        }
        for (size_t i = 0; i < n; ++i)
        {
            io.givenVec<3>(pts[i]);
        }
        return pts;
    }

    // ---- canonicalization ------------------------------------------------
    // order[r] is the index (in upstream's hash-table numbering) of the
    // triangle / tetrahedron whose stored vertex tuple has rank r in
    // lexicographic order. ETManifoldMesh::Insert and TSManifoldMesh::Insert
    // store the tuple in the rotation the caller passed and each feature is
    // inserted exactly once, so the tuple is a total key.
    template <size_t N>
    std::vector<size_t> CanonicalOrder(std::vector<int32_t> const& indices, size_t count)
    {
        std::vector<size_t> order(count);
        std::iota(order.begin(), order.end(), size_t(0));
        std::sort(order.begin(), order.end(), [&indices](size_t a, size_t b)
        {
            for (size_t j = 0; j < N; ++j)
            {
                if (indices[N * a + j] != indices[N * b + j])
                {
                    return indices[N * a + j] < indices[N * b + j];
                }
            }
            return false;
        });
        return order;
    }

    template <size_t N>
    std::vector<size_t> RankOf(std::vector<size_t> const& order)
    {
        std::vector<size_t> rank(order.size());
        for (size_t r = 0; r < order.size(); ++r)
        {
            rank[order[r]] = r;
        }
        return rank;
    }

    // Emit the feature list in canonical order with the adjacency indices
    // remapped through the same permutation.
    template <size_t N>
    void EmitFeatures(oracle::Ctx& io, std::vector<int32_t> const& indices,
        std::vector<int32_t> const& adjacencies, std::vector<size_t> const& order,
        std::vector<size_t> const& rank)
    {
        io.outInt(order.size());
        for (size_t r = 0; r < order.size(); ++r)
        {
            size_t t = order[r];
            for (size_t j = 0; j < N; ++j)
            {
                io.outInt(indices[N * t + j]);
            }
            for (size_t j = 0; j < N; ++j)
            {
                int32_t a = adjacencies[N * t + j];
                io.outInt(a < 0 ? -1 : static_cast<int32_t>(rank[static_cast<size_t>(a)]));
            }
        }
    }

    // Emit a hull index list whose element order upstream leaves unspecified
    // (it follows the hash-table order of the features). The tuples are
    // sorted lexicographically on both sides.
    template <size_t N>
    void EmitSortedTuples(oracle::Ctx& io, std::vector<size_t> const& flat)
    {
        size_t count = flat.size() / N;
        std::vector<std::array<size_t, N>> tuples(count);
        for (size_t i = 0; i < count; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                tuples[i][j] = flat[N * i + j];
            }
        }
        std::sort(tuples.begin(), tuples.end());
        io.outInt(count);
        for (auto const& tuple : tuples)
        {
            for (size_t j = 0; j < N; ++j)
            {
                io.outInt(tuple[j]);
            }
        }
    }

    // ConvexHull2 has no dimension precondition, so its generator adds the
    // two configurations Delaunay2 cannot use: exactly collinear sets (hull
    // dimension 1) and a single repeated point (hull dimension 0).
    //   mode 4: lattice points on a line through a lattice point with a
    //           lattice direction, repeats included
    //   mode 5: one lattice point repeated n times
    void RawHullPoints2(oracle::Ctx& io, int32_t mode, size_t n,
        std::vector<Vector2<double>>& pts)
    {
        if (mode < 4)
        {
            RawPoints2(io, mode, n, pts);
            return;
        }

        pts.resize(n);
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        if (mode == 5)
        {
            for (size_t i = 0; i < n; ++i)
            {
                pts[i][0] = bx;
                pts[i][1] = by;
            }
            return;
        }

        int32_t dx = 0, dy = 0;
        for (int32_t attempt = 0; attempt < 32 && dx == 0 && dy == 0; ++attempt)
        {
            dx = io.rawInteger(-2, 2);
            dy = io.rawInteger(-2, 2);
        }
        if (dx == 0 && dy == 0)
        {
            dx = 1;
            dy = 2;
        }
        for (size_t i = 0; i < n; ++i)
        {
            double k = static_cast<double>(io.rawInteger(-3, 3));
            pts[i][0] = bx + k * static_cast<double>(dx);
            pts[i][1] = by + k * static_cast<double>(dy);
        }
    }
}

// ---- ConvexHull2 ---------------------------------------------------------

// The full public surface of ConvexHull2<T>: operator(), GetDimension(),
// GetLine(), GetNumPoints(), GetNumUniquePoints() and GetHull(). The hull
// index list is emitted in upstream's own order (the divide-and-conquer
// merge order), not canonicalized: the order is part of the result.
//
// The two std facilities on this path are std::sort and std::unique over
// at most 10 indices. MSVC's std::sort is an insertion sort for at most 32
// elements and therefore stable there, which is what makes the port's stable
// Array.prototype.sort agree on inputs with repeated points; see the group
// report.
ORACLE_CASE("ConvexHull2.compute")
{
    int32_t mode = io.index() % 6;
    int32_t n = io.integer(1, 10);
    std::vector<Vector2<double>> pts{};
    RawHullPoints2(io, mode, static_cast<size_t>(n), pts);
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    ConvexHull2<double> hull{};
    bool result = hull(pts);
    io.outBool(result);
    io.outInt(hull.GetDimension());
    io.outInt(hull.GetNumPoints());
    io.outInt(hull.GetNumUniquePoints());
    auto const& indices = hull.GetHull();
    io.outInt(indices.size());
    for (auto index : indices)
    {
        io.outInt(index);
    }
    io.outVec(hull.GetLine().origin);
    io.outVec(hull.GetLine().direction);
}

// ---- Delaunay2 -----------------------------------------------------------

namespace
{
    // The whole result of a 2D triangulation, canonicalized. Also covers the
    // GetIndices(t, array) / GetAdjacencies(t, array) accessor overloads and
    // their out-of-range 'false' return.
    void EmitDelaunay2(oracle::Ctx& io, Delaunay2<double>& del, size_t n)
    {
        io.outInt(del.GetDimension());
        io.outInt(del.GetNumVertices());
        io.outInt(del.GetNumUniqueVertices());
        io.outInt(del.GetNumTriangles());
        auto const& duplicates = del.GetDuplicates();
        io.outInt(duplicates.size());
        for (auto duplicate : duplicates)
        {
            io.outInt(duplicate);
        }
        io.outInt(n);

        size_t numTriangles = del.GetNumTriangles();
        auto const& indices = del.GetIndices();
        auto const& adjacencies = del.GetAdjacencies();
        auto order = CanonicalOrder<3>(indices, numTriangles);
        auto rank = RankOf<3>(order);
        EmitFeatures<3>(io, indices, adjacencies, order, rank);

        std::vector<size_t> hull{};
        bool hullOk = del.GetHull(hull);
        io.outBool(hullOk);
        EmitSortedTuples<2>(io, hull);

        io.outVec(del.GetLine().origin);
        io.outVec(del.GetLine().direction);
        io.outInt(del.GetGraph().GetTriangles().size());
        io.outInt(del.GetGraph().GetEdges().size());

        std::array<int32_t, 3> tri{ 0, 0, 0 }, adj{ 0, 0, 0 };
        bool triOk = del.GetIndices(order[0], tri);
        io.outBool(triOk);
        for (size_t j = 0; j < 3; ++j)
        {
            io.outInt(tri[j]);
        }
        bool adjOk = del.GetAdjacencies(order[0], adj);
        io.outBool(adjOk);
        for (size_t j = 0; j < 3; ++j)
        {
            io.outInt(adj[j] < 0 ? -1 : static_cast<int32_t>(rank[static_cast<size_t>(adj[j])]));
        }
        std::array<int32_t, 3> outOfRange{ 0, 0, 0 };
        bool badOk = del.GetIndices(numTriangles, outOfRange);
        io.outBool(badOk);
        bool badAdjOk = del.GetAdjacencies(numTriangles, outOfRange);
        io.outBool(badAdjOk);
    }
}

// Delaunay2<T>::operator() on inputs where upstream's own floating-point
// dimension classification agrees with exact arithmetic (see Sound2). Modes
// 0-3 are uniform, lattice, cocircular lattice and dense lattice; duplicate
// input vertices are allowed here because the 2D ProcessedVertex hashes and
// compares the vertex only, so upstream's duplicate handling is correct.
ORACLE_CASE("Delaunay2.compute")
{
    int32_t mode = io.index() % 4;
    int32_t n = io.integer(3, 10);
    auto pts = SoundPoints2(io, mode, static_cast<size_t>(n));
    Delaunay2<double> del{};
    bool result = del(pts);
    io.outBool(result);
    EmitDelaunay2(io, del, pts.size());
}

// Degenerate inputs that upstream classifies correctly even with epsilon = 0:
// a single repeated point (dimension 0, the test is 'maxRange == 0') and an
// axis-aligned collinear set (dimension 1, because the normalized direction
// is exactly (+-1, 0) or (0, +-1) and every perpendicular distance is exactly
// zero). GetNumVertices() is NOT emitted here; it is the subject of the
// Delaunay2.compute.deviation.numVertices case.
ORACLE_CASE("Delaunay2.compute.lowDimension")
{
    int32_t mode = io.index() % 3;
    int32_t n = io.integer(1, 8);
    std::vector<Vector2<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double k = (mode == 0 ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
            pts[i][0] = bx + (mode == 1 ? k : 0.0);
            pts[i][1] = by + (mode == 2 ? k : 0.0);
        }
        IntrinsicsVector2<double> info(n, pts.data(), 0.0);
        size_t exact = ExactDimension2(pts);
        accepted = (exact <= 1 && static_cast<size_t>(info.dimension) == exact);
    }
    if (!accepted)
    {
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i] = { 0.0, 0.0 };
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    Delaunay2<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetDimension());
    io.outInt(del.GetNumUniqueVertices());
    io.outInt(del.GetNumTriangles());
    io.outInt(del.GetDuplicates().size());
    io.outVec(del.GetLine().origin);
    io.outVec(del.GetLine().direction);
}

// DELIBERATE DEVIATION (#277). GetNumVertices() returns mIRVertices.size(),
// and operator() clears mIRVertices and returns before resizing it whenever
// the intrinsic dimension is 0 or 1, so upstream reports 0 input vertices
// while GetVertices() still hands back the caller's non-empty array. The port
// returns the stored vertex count on every path. The generator is the
// lowDimension generator, restricted to the inputs upstream classifies
// correctly, so the only difference is the accessor.
ORACLE_CASE("Delaunay2.compute.deviation.numVertices")
{
    int32_t mode = io.index() % 2;
    int32_t n = io.integer(1, 8);
    std::vector<Vector2<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double k = (mode == 0 ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
            pts[i][0] = bx + k;
            pts[i][1] = by;
        }
        IntrinsicsVector2<double> info(n, pts.data(), 0.0);
        size_t exact = ExactDimension2(pts);
        accepted = (exact <= 1 && static_cast<size_t>(info.dimension) == exact);
    }
    if (!accepted)
    {
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i] = { 0.0, 0.0 };
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    Delaunay2<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetDimension());
    io.outInt(del.GetNumVertices());
}

// Throw parity for the two accessors that reject a degenerate result:
// GetHull() calls LogError("The dimension must be 2.") and
// GetContainingTriangle() calls LogAssert(mDimension == 2, ...). Every record
// of this case is a throw record; the generator alternates between the two
// entry points so that both diagnostics are exercised. The inputs are the
// degenerate sets upstream classifies correctly, so the dimension itself is
// not in question.
ORACLE_CASE("Delaunay2.degenerateAccessorsThrow")
{
    int32_t mode = (io.index() / 2) % 2;
    int32_t n = io.integer(1, 8);
    std::vector<Vector2<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double k = (mode == 0 ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
            pts[i][0] = bx + k;
            pts[i][1] = by;
        }
        IntrinsicsVector2<double> info(n, pts.data(), 0.0);
        size_t exact = ExactDimension2(pts);
        accepted = (exact <= 1 && static_cast<size_t>(info.dimension) == exact);
    }
    if (!accepted)
    {
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i] = { 0.0, 0.0 };
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    // Nothing is emitted before the throwing call: the harness discards the
    // outputs of a throw record, so the replay must not produce any either.
    Delaunay2<double> del{};
    bool result = del(pts);
    if (io.index() % 2 == 0)
    {
        std::vector<size_t> hull{};
        bool hullOk = del.GetHull(hull);
        io.outBool(hullOk);
    }
    else
    {
        Vector2<double> q{ 0.5, 0.25 };
        Delaunay2<double>::SearchInfo info{};
        size_t found = del.GetContainingTriangle(q, info);
        io.outInt(found == Delaunay2<double>::negOne ? -1 : static_cast<int32_t>(found));
    }
    io.outBool(result);
    io.outInt(del.GetDimension());
}

// DELIBERATE DEVIATION (#391). Delaunay2<T> constructs IntrinsicsVector2<T>
// with a hardcoded epsilon of 0 and the intrinsics measure distances against
// a *normalized* frame, so an exactly collinear set whose direction is not
// axis aligned has perpendicular distances that are nonzero by an ulp: the
// strict 'maxDistance <= 0' test fails and the set is reported as
// 2-dimensional with an exactly degenerate extreme triangle. Upstream then
// either returns true with a degenerate triangle list or throws from the
// incremental update. The port classifies the dimension with its own exact
// ToLine predicate and reports dimension 1.
//
// The generator draws exactly collinear lattice sets and accepts only those
// on which upstream's own IntrinsicsVector2 reports dimension 2, which is the
// exact separator: when it reports 1 the port agrees with it bit for bit.
ORACLE_CASE("Delaunay2.compute.deviation.epsilon")
{
    int32_t n = io.integer(2, 7);
    std::vector<Vector2<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        int32_t dx = 0, dy = 0;
        for (int32_t k = 0; k < 32 && (dx == 0 || dy == 0); ++k)
        {
            dx = io.rawInteger(-3, 3);
            dy = io.rawInteger(-3, 3);
        }
        if (dx == 0 || dy == 0)
        {
            dx = -3;
            dy = -9;
        }
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double k = static_cast<double>(io.rawInteger(-3, 3));
            pts[i][0] = bx + k * static_cast<double>(dx);
            pts[i][1] = by + k * static_cast<double>(dy);
        }
        IntrinsicsVector2<double> info(n, pts.data(), 0.0);
        accepted = (info.dimension == 2);
    }
    if (!accepted)
    {
        // The reproduction recorded in docs/UPSTREAM-FINDINGS.md, extended to
        // n points along the same line.
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i][0] = -3.0 * static_cast<double>(i);
            pts[i][1] = -9.0 * static_cast<double>(i);
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    Delaunay2<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetDimension());
    io.outInt(del.GetNumTriangles());
}

namespace
{
    // A query point aimed at a specific configuration relative to the
    // triangulation. Raw draws only; the caller records the result. The
    // triangle is selected by its canonical rank so that the C++ side and
    // the replay would agree, but the replay does not need to reproduce the
    // selection: it reads the recorded point.
    Vector2<double> RawQuery2(oracle::Ctx& io, int32_t qmode,
        std::vector<Vector2<double>> const& pts, Delaunay2<double> const& del,
        std::vector<size_t> const& order)
    {
        auto const& indices = del.GetIndices();
        int32_t r = io.rawInteger(0, static_cast<int32_t>(order.size()) - 1);
        size_t t = order[static_cast<size_t>(r)];
        Vector2<double> v0 = pts[static_cast<size_t>(indices[3 * t])];
        Vector2<double> v1 = pts[static_cast<size_t>(indices[3 * t + 1])];
        Vector2<double> v2 = pts[static_cast<size_t>(indices[3 * t + 2])];
        Vector2<double> q{ 0.0, 0.0 };
        if (qmode == 0)
        {
            // Strictly inside the triangle; dyadic weights, so the point is
            // exact for a lattice input.
            for (int32_t j = 0; j < 2; ++j)
            {
                q[j] = 0.25 * v0[j] + 0.25 * v1[j] + 0.5 * v2[j];
            }
        }
        else if (qmode == 1)
        {
            // Exactly on the edge <v0,v1>: ToLine returns 0 there.
            for (int32_t j = 0; j < 2; ++j)
            {
                q[j] = 0.5 * (v0[j] + v1[j]);
            }
        }
        else if (qmode == 2)
        {
            // Exactly at a vertex.
            q = v0;
        }
        else if (qmode == 3)
        {
            // Far outside the hull, so the search exits through a hull edge.
            q[0] = 40.0 + static_cast<double>(io.rawInteger(-5, 5));
            q[1] = static_cast<double>(io.rawInteger(-40, 40));
        }
        else
        {
            // Uniform over a box a little larger than the data.
            q[0] = io.raw(-6.0, 6.0);
            q[1] = io.raw(-6.0, 6.0);
        }
        return q;
    }

    void EmitSearch2(oracle::Ctx& io, Delaunay2<double> const& del,
        Vector2<double> const& q, Delaunay2<double>::SearchInfo& info,
        std::vector<size_t> const& rank)
    {
        size_t found = del.GetContainingTriangle(q, info);
        io.outInt(found == Delaunay2<double>::negOne
            ? -1 : static_cast<int32_t>(rank[found]));
        io.outInt(info.numPath);
        for (size_t i = 0; i < info.numPath; ++i)
        {
            io.outInt(rank[info.path[i]]);
        }
        io.outInt(rank[info.finalTriangle]);
        for (size_t j = 0; j < 3; ++j)
        {
            io.outInt(info.finalV[j]);
        }
        io.outInt(rank[info.initialTriangle]);
    }
}

// Delaunay2<T>::GetContainingTriangle with the search started at the
// canonically first triangle on both sides (upstream's default start is
// 'the first triangle of the hash table', which is not a portable choice).
// The query points are aimed at the interior of a triangle, exactly at an
// edge midpoint, exactly at a vertex, far outside the hull, and uniformly
// over a box larger than the data.
ORACLE_CASE("Delaunay2.getContainingTriangle")
{
    int32_t mode = io.index() % 3;
    int32_t qmode = (io.index() / 3) % 5;
    int32_t n = io.integer(4, 10);
    auto pts = SoundPoints2(io, mode, static_cast<size_t>(n));
    Delaunay2<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetNumTriangles());

    auto order = CanonicalOrder<3>(del.GetIndices(), del.GetNumTriangles());
    auto rank = RankOf<3>(order);
    Vector2<double> q = RawQuery2(io, qmode, pts, del, order);
    io.givenVec<2>(q);

    Delaunay2<double>::SearchInfo info{};
    info.initialTriangle = order[0];
    EmitSearch2(io, del, q, info, rank);
}

// The same search with the default SearchInfo, whose initialTriangle is
// negOne and is replaced by 0 inside the query. Exactly 3 input points give
// exactly one triangle, so 'triangle 0' is the same triangle on both sides
// and the hash-table order cannot enter the comparison.
ORACLE_CASE("Delaunay2.getContainingTriangle.defaultStart")
{
    int32_t qmode = io.index() % 5;
    auto pts = SoundPoints2(io, io.index() % 2, 3);
    Delaunay2<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetNumTriangles());

    auto order = CanonicalOrder<3>(del.GetIndices(), del.GetNumTriangles());
    auto rank = RankOf<3>(order);
    Vector2<double> q = RawQuery2(io, qmode, pts, del, order);
    io.givenVec<2>(q);

    Delaunay2<double>::SearchInfo info{};
    EmitSearch2(io, del, q, info, rank);
}

// ---- Delaunay3 -----------------------------------------------------------

namespace
{
    // The whole result of a 3D tetrahedralization, canonicalized. Also covers
    // the GetIndices(t, array) / GetAdjacencies(t, array) accessor overloads
    // and their out-of-range 'false' return.
    void EmitDelaunay3(oracle::Ctx& io, Delaunay3<double>& del, size_t n)
    {
        io.outInt(del.GetDimension());
        io.outInt(del.GetNumVertices());
        io.outInt(del.GetNumUniqueVertices());
        io.outInt(del.GetNumTetrahedra());
        auto const& duplicates = del.GetDuplicates();
        io.outInt(duplicates.size());
        for (auto duplicate : duplicates)
        {
            io.outInt(duplicate);
        }
        io.outInt(n);

        size_t numTetrahedra = del.GetNumTetrahedra();
        auto const& indices = del.GetIndices();
        auto const& adjacencies = del.GetAdjacencies();
        auto order = CanonicalOrder<4>(indices, numTetrahedra);
        auto rank = RankOf<4>(order);
        EmitFeatures<4>(io, indices, adjacencies, order, rank);

        std::vector<size_t> hull{};
        bool hullOk = del.GetHull(hull);
        io.outBool(hullOk);
        EmitSortedTuples<3>(io, hull);

        io.outVec(del.GetLine().origin);
        io.outVec(del.GetLine().direction);
        io.outVec(del.GetPlane().normal);
        io.outReal(del.GetPlane().constant);
        io.outInt(del.GetGraph().GetTetrahedra().size());
        io.outInt(del.GetGraph().GetTriangles().size());

        std::array<int32_t, 4> tetra{ 0, 0, 0, 0 }, adj{ 0, 0, 0, 0 };
        bool tetraOk = del.GetIndices(order[0], tetra);
        io.outBool(tetraOk);
        for (size_t j = 0; j < 4; ++j)
        {
            io.outInt(tetra[j]);
        }
        bool adjOk = del.GetAdjacencies(order[0], adj);
        io.outBool(adjOk);
        for (size_t j = 0; j < 4; ++j)
        {
            io.outInt(adj[j] < 0 ? -1
                : static_cast<int32_t>(rank[static_cast<size_t>(adj[j])]));
        }
        std::array<int32_t, 4> outOfRange{ 0, 0, 0, 0 };
        bool badOk = del.GetIndices(numTetrahedra, outOfRange);
        io.outBool(badOk);
        bool badAdjOk = del.GetAdjacencies(numTetrahedra, outOfRange);
        io.outBool(badAdjOk);
    }
}

// Delaunay3<T>::operator() on inputs where upstream's own floating-point
// dimension classification agrees with exact arithmetic and where no input
// vertex is repeated (see Sound3; repeated vertices are the subject of the
// #283 deviation case). Modes 0-3 are uniform, lattice, cospherical lattice
// and dense lattice.
ORACLE_CASE("Delaunay3.compute")
{
    int32_t mode = io.index() % 4;
    int32_t n = io.integer(4, 10);
    auto pts = SoundPoints3(io, mode, static_cast<size_t>(n));
    Delaunay3<double> del{};
    bool result = del(pts);
    io.outBool(result);
    EmitDelaunay3(io, del, pts.size());
}

// Degenerate inputs that upstream classifies correctly even with epsilon = 0:
// a single repeated point (dimension 0), a collinear set along a coordinate
// axis (dimension 1) and a coplanar set in a plane z = c (dimension 2). In
// both degenerate cases the frame the intrinsics build has exact zeros in the
// right places, so every distance is exactly zero and the strict test holds.
// GetNumVertices() is NOT emitted here; it is the subject of the
// Delaunay3.compute.deviation.numVertices case.
ORACLE_CASE("Delaunay3.compute.lowDimension")
{
    int32_t mode = io.index() % 3;
    int32_t n = io.integer(1, 8);
    std::vector<Vector3<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        double bz = static_cast<double>(io.rawInteger(-3, 3));
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double k0 = (mode == 0 ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
            double k1 = (mode == 2 ? static_cast<double>(io.rawInteger(-3, 3)) : 0.0);
            pts[i][0] = bx + k0;
            pts[i][1] = by + k1;
            pts[i][2] = bz;
        }
        IntrinsicsVector3<double> info(n, pts.data(), 0.0);
        size_t exact = ExactDimension3(pts);
        accepted = (exact <= 2 && static_cast<size_t>(info.dimension) == exact);
    }
    if (!accepted)
    {
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i] = { 0.0, 0.0, 0.0 };
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<3>(pts[i]);
    }

    Delaunay3<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetDimension());
    io.outInt(del.GetNumUniqueVertices());
    io.outInt(del.GetNumTetrahedra());
    io.outInt(del.GetDuplicates().size());
    io.outVec(del.GetLine().origin);
    io.outVec(del.GetLine().direction);
    io.outVec(del.GetPlane().normal);
    io.outReal(del.GetPlane().constant);
}

namespace
{
    Vector3<double> RawQuery3(oracle::Ctx& io, int32_t qmode,
        std::vector<Vector3<double>> const& pts, Delaunay3<double> const& del,
        std::vector<size_t> const& order)
    {
        auto const& indices = del.GetIndices();
        int32_t r = io.rawInteger(0, static_cast<int32_t>(order.size()) - 1);
        size_t t = order[static_cast<size_t>(r)];
        Vector3<double> v0 = pts[static_cast<size_t>(indices[4 * t])];
        Vector3<double> v1 = pts[static_cast<size_t>(indices[4 * t + 1])];
        Vector3<double> v2 = pts[static_cast<size_t>(indices[4 * t + 2])];
        Vector3<double> v3 = pts[static_cast<size_t>(indices[4 * t + 3])];
        Vector3<double> q{ 0.0, 0.0, 0.0 };
        if (qmode == 0)
        {
            // Strictly inside; dyadic weights 1/8, 1/8, 1/4, 1/2.
            for (int32_t j = 0; j < 3; ++j)
            {
                q[j] = 0.125 * v0[j] + 0.125 * v1[j] + 0.25 * v2[j] + 0.5 * v3[j];
            }
        }
        else if (qmode == 1)
        {
            // Exactly on the face <v0,v1,v2>: ToPlane returns 0 there.
            for (int32_t j = 0; j < 3; ++j)
            {
                q[j] = 0.25 * v0[j] + 0.25 * v1[j] + 0.5 * v2[j];
            }
        }
        else if (qmode == 2)
        {
            // Exactly on the edge <v0,v1>.
            for (int32_t j = 0; j < 3; ++j)
            {
                q[j] = 0.5 * (v0[j] + v1[j]);
            }
        }
        else if (qmode == 3)
        {
            // Exactly at a vertex.
            q = v0;
        }
        else if (qmode == 4)
        {
            // Far outside the hull.
            q[0] = 40.0 + static_cast<double>(io.rawInteger(-5, 5));
            q[1] = static_cast<double>(io.rawInteger(-40, 40));
            q[2] = static_cast<double>(io.rawInteger(-40, 40));
        }
        else
        {
            for (int32_t j = 0; j < 3; ++j)
            {
                q[j] = io.raw(-6.0, 6.0);
            }
        }
        return q;
    }

    void EmitSearch3(oracle::Ctx& io, Delaunay3<double> const& del,
        Vector3<double> const& q, Delaunay3<double>::SearchInfo& info,
        std::vector<size_t> const& rank)
    {
        size_t found = del.GetContainingTetrahedron(q, info);
        io.outInt(found == Delaunay3<double>::negOne
            ? -1 : static_cast<int32_t>(rank[found]));
        io.outInt(info.numPath);
        for (size_t i = 0; i < info.numPath; ++i)
        {
            io.outInt(rank[info.path[i]]);
        }
        io.outInt(rank[info.finalTetrahedron]);
        for (size_t j = 0; j < 4; ++j)
        {
            io.outInt(info.finalV[j]);
        }
        io.outInt(rank[info.initialTetrahedron]);
    }
}

// Delaunay3<T>::GetContainingTetrahedron with the search started at the
// canonically first tetrahedron on both sides. The query points are aimed at
// the interior, exactly on a face, exactly on an edge, exactly at a vertex,
// far outside the hull, and uniformly over a box larger than the data.
ORACLE_CASE("Delaunay3.getContainingTetrahedron")
{
    int32_t mode = io.index() % 3;
    int32_t qmode = (io.index() / 3) % 6;
    int32_t n = io.integer(5, 10);
    auto pts = SoundPoints3(io, mode, static_cast<size_t>(n));
    Delaunay3<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetNumTetrahedra());

    auto order = CanonicalOrder<4>(del.GetIndices(), del.GetNumTetrahedra());
    auto rank = RankOf<4>(order);
    Vector3<double> q = RawQuery3(io, qmode, pts, del, order);
    io.givenVec<3>(q);

    Delaunay3<double>::SearchInfo info{};
    info.initialTetrahedron = order[0];
    EmitSearch3(io, del, q, info, rank);
}

// The same search with the default SearchInfo, whose initialTetrahedron is 0.
// Exactly 4 input points give exactly one tetrahedron, so 'tetrahedron 0' is
// the same tetrahedron on both sides.
ORACLE_CASE("Delaunay3.getContainingTetrahedron.defaultStart")
{
    int32_t qmode = io.index() % 6;
    auto pts = SoundPoints3(io, io.index() % 2, 4);
    Delaunay3<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetNumTetrahedra());

    auto order = CanonicalOrder<4>(del.GetIndices(), del.GetNumTetrahedra());
    auto rank = RankOf<4>(order);
    Vector3<double> q = RawQuery3(io, qmode, pts, del, order);
    io.givenVec<3>(q);

    Delaunay3<double>::SearchInfo info{};
    EmitSearch3(io, del, q, info, rank);
}

// Throw parity for the 3D analogues: GetHull() calls
// LogError("The dimension must be 3.") and GetContainingTetrahedron() calls
// LogAssert(mDimension == 3, ...). Every record of this case is a throw
// record.
ORACLE_CASE("Delaunay3.degenerateAccessorsThrow")
{
    int32_t mode = (io.index() / 2) % 3;
    int32_t n = io.integer(1, 8);
    std::vector<Vector3<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        double bz = static_cast<double>(io.rawInteger(-3, 3));
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double k0 = (mode == 0 ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
            double k1 = (mode == 2 ? static_cast<double>(io.rawInteger(-3, 3)) : 0.0);
            pts[i][0] = bx + k0;
            pts[i][1] = by + k1;
            pts[i][2] = bz;
        }
        IntrinsicsVector3<double> info(n, pts.data(), 0.0);
        size_t exact = ExactDimension3(pts);
        accepted = (exact <= 2 && static_cast<size_t>(info.dimension) == exact);
    }
    if (!accepted)
    {
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i] = { 0.0, 0.0, 0.0 };
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<3>(pts[i]);
    }

    Delaunay3<double> del{};
    bool result = del(pts);
    if (io.index() % 2 == 0)
    {
        std::vector<size_t> hull{};
        bool hullOk = del.GetHull(hull);
        io.outBool(hullOk);
    }
    else
    {
        Vector3<double> q{ 0.5, 0.25, 0.125 };
        Delaunay3<double>::SearchInfo info{};
        size_t found = del.GetContainingTetrahedron(q, info);
        io.outInt(found == Delaunay3<double>::negOne ? -1 : static_cast<int32_t>(found));
    }
    io.outBool(result);
    io.outInt(del.GetDimension());
}

// DELIBERATE DEVIATION (#283). Delaunay3<T>::ProcessedVertex hashes AND
// compares its 'location' member, and the lookup key is built as
// ProcessedVertex(mVertices[i], i), so a repeated coordinate triple at a
// later index carries a different location and the find always misses.
// mDuplicates[] stays the identity map, GetNumUniqueVertices() equals the
// input count, and Update(i) runs a second time for a point that is already
// a mesh vertex: it removes that point's incident tetrahedra (the point lies
// on their circumspheres) and rebuilds them against a coincident vertex
// index. The port hashes and compares the vertex only, as Delaunay2 does.
//
// The generator draws point sets with at least one repeated vertex whose
// seed tetrahedron is sound (SoundSeed3), so the only defect on these inputs
// is the duplicate detection.
ORACLE_CASE("Delaunay3.compute.deviation.duplicates")
{
    int32_t mode = io.index() % 2;
    int32_t n = io.integer(5, 9);
    std::vector<Vector3<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        RawPoints3(io, mode == 0 ? 1 : 2, static_cast<size_t>(n), pts);
        int32_t source = io.rawInteger(0, n - 2);
        int32_t target = io.rawInteger(source + 1, n - 1);
        pts[static_cast<size_t>(target)] = pts[static_cast<size_t>(source)];
        accepted = HasDuplicate3(pts) && SoundSeed3(pts);
    }
    if (!accepted)
    {
        Fallback3(static_cast<size_t>(n), pts);
        pts[static_cast<size_t>(n) - 1] = pts[0];
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<3>(pts[i]);
    }

    Delaunay3<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetDimension());
    io.outInt(del.GetNumUniqueVertices());
    io.outInt(del.GetDuplicates().size());
    for (auto duplicate : del.GetDuplicates())
    {
        io.outInt(duplicate);
    }
    io.outInt(del.GetNumTetrahedra());
}

// DELIBERATE DEVIATION (#391), the 3D half. An exactly collinear or exactly
// coplanar set whose frame is not axis aligned has distances to the line or
// the plane that are nonzero by an ulp, so the strict 'maxDistance <= 0' test
// fails and upstream reports dimension 3 with a zero-volume seed tetrahedron
// (it then returns true with zero tetrahedra, or throws from the incremental
// update). The generator accepts only the sets on which upstream's own
// IntrinsicsVector3 disagrees with exact arithmetic, which is the exact
// separator.
ORACLE_CASE("Delaunay3.compute.deviation.epsilon")
{
    int32_t n = io.integer(3, 7);
    std::vector<Vector3<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        // A lattice plane through the origin with normal (a,b,c): the points
        // are u * (b,-a,0) + v * (c,0,-a), which is exactly coplanar.
        int32_t a = io.rawInteger(1, 3);
        int32_t b = io.rawInteger(-3, 3);
        int32_t c = io.rawInteger(-3, 3);
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double u = static_cast<double>(io.rawInteger(-3, 3));
            double v = static_cast<double>(io.rawInteger(-3, 3));
            pts[i][0] = u * static_cast<double>(b) + v * static_cast<double>(c);
            pts[i][1] = u * static_cast<double>(-a);
            pts[i][2] = v * static_cast<double>(-a);
        }
        IntrinsicsVector3<double> info(n, pts.data(), 0.0);
        size_t exact = ExactDimension3(pts);
        accepted = (exact <= 2 && static_cast<size_t>(info.dimension) != exact);
    }
    if (!accepted)
    {
        // The reproduction recorded in docs/UPSTREAM-FINDINGS.md (z = x + 2y),
        // truncated or extended to n points.
        std::array<std::array<double, 3>, 7> const seed
        { {
            { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 }, { 0.0, 1.0, 2.0 },
            { 1.0, 1.0, 3.0 }, { 2.0, 1.0, 4.0 }, { 2.0, 2.0, 6.0 },
            { 3.0, 1.0, 5.0 }
        } };
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i] = { seed[i][0], seed[i][1], seed[i][2] };
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<3>(pts[i]);
    }

    Delaunay3<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetDimension());
    io.outInt(del.GetNumTetrahedra());
}

// DELIBERATE DEVIATION (#283). GetNumVertices() returns mIRVertices.size(),
// which is 0 for dimension 0, 1 or 2 input because operator() returns before
// resizing it, while GetVertices() still returns the caller's array. The port
// returns the stored vertex count on every path. The generator is the
// lowDimension generator, restricted to inputs upstream classifies correctly.
ORACLE_CASE("Delaunay3.compute.deviation.numVertices")
{
    int32_t mode = io.index() % 3;
    int32_t n = io.integer(1, 8);
    std::vector<Vector3<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        double bz = static_cast<double>(io.rawInteger(-3, 3));
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double k0 = (mode == 0 ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
            double k1 = (mode == 2 ? static_cast<double>(io.rawInteger(-3, 3)) : 0.0);
            pts[i][0] = bx + k0;
            pts[i][1] = by + k1;
            pts[i][2] = bz;
        }
        IntrinsicsVector3<double> info(n, pts.data(), 0.0);
        size_t exact = ExactDimension3(pts);
        accepted = (exact <= 2 && static_cast<size_t>(info.dimension) == exact);
    }
    if (!accepted)
    {
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i] = { 0.0, 0.0, 0.0 };
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<3>(pts[i]);
    }

    Delaunay3<double> del{};
    bool result = del(pts);
    io.outBool(result);
    io.outInt(del.GetDimension());
    io.outInt(del.GetNumVertices());
}
