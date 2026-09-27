// Verify group 11 (computational geometry): differential cases for
// ConstrainedDelaunay2.h, ConvexHull3.h, Delaunay2Mesh.h, Delaunay3Mesh.h,
// MinimumAreaBox2.h, MinimumWidthPoints2.h and SeparatePoints2.h.
//
// INSTANTIATIONS. The port implements the *replacement* class templates, so
// the C++ side instantiates the same ones:
//   ConstrainedDelaunay2<double>            (not the deprecated
//                                            <InputType, ComputeType>)
//   ConvexHull3<double>, lgNumThreads = 0   (the port has no threads)
//   Delaunay2Mesh<double>, Delaunay3Mesh<double>
//   MinimumAreaBox2<double, BSRational<UIntegerAP32>>
//       The port is rational-only by design: upstream's own header says
//       "This algorithm guarantees a correct output only when ComputeType is
//       an exact arithmetic type that supports division. In GTE, one such
//       type is BSRational<UIntegerAP32>." The port's BSRational is
//       bigint-backed, which is the same arithmetic. Only
//       MinimumAreaBox2.deviation.floatComputeType instantiates
//       MinimumAreaBox2<double, double>, to show what the choice changes.
//   MinimumWidthPoints2<double>             (Rational = BSRational<UIntegerAP32>
//                                            internally, as the port)
//   SeparatePoints2<double, double>         (ComputeType is vestigial)
//
// COMPARISON POLICY. Everything compared here is either combinatorial
// (indices, adjacencies, dimensions, counts, predicate signs) or computed
// with + - * / sqrt only, either on doubles or on exact rationals that are
// rounded to double at the very end. No libm call is on any compared path,
// so every case is declared { exact: true }.
//
// CANONICALIZATION. ETManifoldMesh, VETManifoldMesh and TSManifoldMesh store
// their features in std::unordered_map, and ConstrainedDelaunay2's inserted
// edges live in a std::unordered_set, so upstream enumerates them in
// hash-table order while the port enumerates them in sorted-key order. The
// feature *sets* are order independent; only the numbering is not. Triangle
// and tetrahedron lists are therefore sorted by their stored vertex tuple
// (ETManifoldMesh::Insert keeps the rotation the caller passed, so the tuple
// is a total key) and adjacency indices are remapped through the same
// permutation, exactly as v09-compgeom does. Vertex and edge lists whose
// order upstream leaves unspecified are sorted. Everything else is emitted
// verbatim.
//
// One io draw per C++ statement: MSVC evaluates function arguments right to
// left, so every generated value goes into a named local before it is used.
#define ORACLE_FAMILY "v11-compgeom"
#include "Oracle.h"

#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/ConstrainedDelaunay2.h>
#include <Mathematics/ConvexHull3.h>
#include <Mathematics/Delaunay2Mesh.h>
#include <Mathematics/Delaunay3Mesh.h>
#include <Mathematics/MinimumAreaBox2.h>
#include <Mathematics/MinimumWidthPoints2.h>
#include <Mathematics/SeparatePoints2.h>

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
    using ExactQ = BSRational<UIntegerAP32>;

    // Sign of DotPerp(B - A, C - A): +1 when <A,B,C> is counterclockwise.
    int32_t ExactOrient2(Vector2<double> const& A, Vector2<double> const& B,
        Vector2<double> const& C)
    {
        Exact ax(A[0]), ay(A[1]), bx(B[0]), by(B[1]), cx(C[0]), cy(C[1]);
        Exact x0 = bx - ax, y0 = by - ay, x1 = cx - ax, y1 = cy - ay;
        Exact det = x0 * y1 - x1 * y0;
        return det.GetSign();
    }

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

    // ---- generator probes (upstream's own control flow) ------------------
    // Delaunay2<T>::operator() builds IntrinsicsVector2<T> with a hardcoded
    // epsilon of 0 and trusts info.dimension, info.extreme[] and
    // info.extremeCCW, all decided by floating-point distances to a
    // *normalized* frame; the port replaces that classification with its own
    // exact ToLine predicate (upstream finding #391). This probe is the exact
    // separator used by v09-compgeom: it accepts exactly those inputs on
    // which upstream's floating-point classification agrees with exact
    // arithmetic, so on an accepted input the port's seed triangle is
    // upstream's seed triangle and the two algorithms are the same algorithm.
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

    // Delaunay3<T> additionally has the duplicate-detection defect (#283), so
    // a 3D point set is "sound" only when no coordinate triple is repeated.
    bool Sound3(std::vector<Vector3<double>> const& pts)
    {
        return !HasDuplicate3(pts) && SoundSeed3(pts);
    }

    // ---- point-set generators (raw, unrecorded) --------------------------
    // The 12 lattice points on the circle x^2 + y^2 = 25 and the 30 lattice
    // points on the sphere x^2 + y^2 + z^2 = 9. Subsets of these are exactly
    // cocircular / cospherical, so the incircle / insphere predicates see an
    // interval that straddles zero and the exact rational fallback decides
    // the sign. Adding the center makes the configuration non-degenerate
    // again, so the mode reaches both.
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
    // mode 2: cocircular lattice about a lattice center, center mixed in.
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

    // Fallbacks for an expired rejection loop: points on a parabola / on the
    // moment curve are in general position, so upstream's floating-point
    // classification is never in doubt there.
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

    // Draw a 2D point set on which upstream's own dimension classification is
    // sound. The loop is capped and every coordinate is redrawn on each
    // attempt. The accepted set is recorded here, once.
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
    // feature whose stored vertex tuple has rank r in lexicographic order.
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

    // Emit a flat list of N-tuples whose element order upstream leaves
    // unspecified (it follows hash-table order). Sorted on both sides.
    template <size_t N, typename T>
    void EmitSortedTuples(oracle::Ctx& io, std::vector<T> const& flat)
    {
        size_t count = flat.size() / N;
        std::vector<std::array<T, N>> tuples(count);
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
}

// ---- ConvexHull3 ---------------------------------------------------------

namespace
{
    // ConvexHull3 classifies every predicate exactly (Colocated, Colinear and
    // ToPlane all fall back to rational arithmetic), so there is no epsilon
    // probe here and no soundness restriction: all of dimension 0, 1, 2 and 3
    // are generated directly.
    //   modes 0-3: the shared 3D generators (uniform, lattice, cospherical
    //              lattice, dense lattice)
    //   mode 4:    exactly collinear lattice points, repeats included
    //   mode 5:    exactly coplanar lattice points in a lattice plane whose
    //              normal is not axis aligned (the hull2 projection path)
    //   mode 6:    one lattice point repeated n times (dimension 0)
    void RawHullPoints3(oracle::Ctx& io, int32_t mode, size_t n,
        std::vector<Vector3<double>>& pts)
    {
        if (mode < 4)
        {
            RawPoints3(io, mode, n, pts);
            return;
        }

        pts.resize(n);
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        double bz = static_cast<double>(io.rawInteger(-3, 3));
        if (mode == 6)
        {
            for (size_t i = 0; i < n; ++i)
            {
                pts[i] = { bx, by, bz };
            }
            return;
        }

        if (mode == 4)
        {
            int32_t dx = 0, dy = 0, dz = 0;
            for (int32_t attempt = 0; attempt < 32 && dx == 0 && dy == 0 && dz == 0;
                ++attempt)
            {
                dx = io.rawInteger(-2, 2);
                dy = io.rawInteger(-2, 2);
                dz = io.rawInteger(-2, 2);
            }
            if (dx == 0 && dy == 0 && dz == 0)
            {
                dx = 1;
                dy = 2;
                dz = -1;
            }
            for (size_t i = 0; i < n; ++i)
            {
                double k = static_cast<double>(io.rawInteger(-3, 3));
                pts[i][0] = bx + k * static_cast<double>(dx);
                pts[i][1] = by + k * static_cast<double>(dy);
                pts[i][2] = bz + k * static_cast<double>(dz);
            }
            return;
        }

        // mode 5: the lattice plane through (bx,by,bz) spanned by (b,-a,0)
        // and (c,0,-a), which is exactly planar for integer a, b, c.
        int32_t a = io.rawInteger(1, 3);
        int32_t b = io.rawInteger(-3, 3);
        int32_t c = io.rawInteger(-3, 3);
        for (size_t i = 0; i < n; ++i)
        {
            double u = static_cast<double>(io.rawInteger(-3, 3));
            double v = static_cast<double>(io.rawInteger(-3, 3));
            pts[i][0] = bx + u * static_cast<double>(b) + v * static_cast<double>(c);
            pts[i][1] = by + u * static_cast<double>(-a);
            pts[i][2] = bz + v * static_cast<double>(-a);
        }
    }
}

// The full public surface of ConvexHull3<T> on the single-threaded path:
// operator()(points, 0), GetDimension(), GetVertices(), GetHull() and
// GetHullMesh(). For dimension 3 the vertex list and the triangle list come
// from std::unordered_map iteration, so they are sorted on both sides (the
// vertex list ascending, the triangle list by its stored vertex triple, which
// TriangleKey<true> has already rotated to put the minimum first). For
// dimension 0, 1 and 2 the lists are ordered results and are emitted verbatim.
//
// The two std facilities on this path are std::sort and std::unique over at
// most 10 indices with the lexicographic Vector3 comparator. MSVC's std::sort
// is an insertion sort, and therefore stable, for at most 32 elements, which
// is what makes the port's stable Array.prototype.sort pick the same
// representative of a repeated point; see the group report.
ORACLE_CASE("ConvexHull3.compute")
{
    int32_t mode = io.index() % 7;
    int32_t n = io.integer(1, 10);
    std::vector<Vector3<double>> pts{};
    RawHullPoints3(io, mode, static_cast<size_t>(n), pts);
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<3>(pts[i]);
    }

    ConvexHull3<double> hull{};
    hull(pts, 0);
    size_t dimension = hull.GetDimension();
    io.outInt(dimension);

    std::vector<size_t> vertices = hull.GetVertices();
    std::vector<size_t> indices = hull.GetHull();
    if (dimension == 3)
    {
        std::sort(vertices.begin(), vertices.end());
        io.outInt(vertices.size());
        for (auto v : vertices)
        {
            io.outInt(v);
        }
        EmitSortedTuples<3, size_t>(io, indices);
    }
    else
    {
        io.outInt(vertices.size());
        for (auto v : vertices)
        {
            io.outInt(v);
        }
        io.outInt(indices.size());
        for (auto v : indices)
        {
            io.outInt(v);
        }
    }

    auto const& mesh = hull.GetHullMesh();
    io.outInt(mesh.GetVertices().size());
    io.outInt(mesh.GetEdges().size());
    io.outInt(mesh.GetTriangles().size());
}

// ---- Delaunay2Mesh -------------------------------------------------------

namespace
{
    // A query point aimed at a configuration whose containing triangle is
    // unique, so that the answer does not depend on where the walk started.
    // Delaunay2Mesh::GetContainingTriangle builds a default SearchInfo, whose
    // initialTriangle is negOne and is replaced by triangle 0 inside the
    // query; "triangle 0" is the first triangle of upstream's hash table and
    // of the port's sorted map, which are different triangles. Only the final
    // answer is emitted, and only for query points that are strictly inside
    // one triangle or strictly outside the hull, where the answer is the same
    // for every start. Query points exactly on an edge or at a vertex are
    // covered by v09-compgeom's Delaunay2.getContainingTriangle cases, which
    // pin the start triangle canonically.
    //   qmode 0: strictly inside the triangle of canonical rank r
    //   qmode 1: far outside the hull
    //   qmode 2: uniform over a box a little larger than the data
    Vector2<double> RawMeshQuery2(oracle::Ctx& io, int32_t qmode,
        std::vector<Vector2<double>> const& pts, Delaunay2<double> const& del,
        std::vector<size_t> const& order)
    {
        Vector2<double> q{ 0.0, 0.0 };
        if (qmode == 0)
        {
            auto const& indices = del.GetIndices();
            int32_t r = io.rawInteger(0, static_cast<int32_t>(order.size()) - 1);
            size_t t = order[static_cast<size_t>(r)];
            Vector2<double> v0 = pts[static_cast<size_t>(indices[3 * t])];
            Vector2<double> v1 = pts[static_cast<size_t>(indices[3 * t + 1])];
            Vector2<double> v2 = pts[static_cast<size_t>(indices[3 * t + 2])];
            for (int32_t j = 0; j < 2; ++j)
            {
                q[j] = 0.25 * v0[j] + 0.25 * v1[j] + 0.5 * v2[j];
            }
        }
        else if (qmode == 1)
        {
            q[0] = 40.0 + static_cast<double>(io.rawInteger(-5, 5));
            q[1] = static_cast<double>(io.rawInteger(-40, 40));
        }
        else
        {
            q[0] = io.raw(-6.0, 6.0);
            q[1] = io.raw(-6.0, 6.0);
        }
        return q;
    }
}

// The full public surface of Delaunay2Mesh<T>: GetNumVertices(),
// GetNumTriangles(), GetVertices(), GetIndices(), GetAdjacencies(),
// GetInvalidIndex(), GetContainingTriangle(P), GetVertices(t, out),
// GetIndices(t, out), GetAdjacencies(t, out) and GetBarycentrics(t, P),
// including the out-of-range 'false' return of the per-triangle accessors.
// The triangle numbering is canonicalized as in v09-compgeom; the stored
// vertex triple of a triangle is the same on both sides because
// TriangleKey<true> rotates the minimum vertex to the front, so the
// barycentric coordinates are taken with respect to the same ordered
// triangle.
ORACLE_CASE("Delaunay2Mesh.query")
{
    int32_t mode = io.index() % 4;
    int32_t qmode = (io.index() / 4) % 3;
    int32_t n = io.integer(3, 10);
    auto pts = SoundPoints2(io, mode, static_cast<size_t>(n));

    Delaunay2<double> del{};
    bool built = del(pts);
    io.outBool(built);

    Delaunay2Mesh<double> mesh(del);
    io.outInt(mesh.GetNumVertices());
    io.outInt(mesh.GetNumTriangles());

    size_t numTriangles = mesh.GetNumTriangles();
    auto const* vertices = mesh.GetVertices();
    for (size_t i = 0; i < static_cast<size_t>(n); ++i)
    {
        io.outVec(vertices[i]);
    }

    std::vector<int32_t> indices(3 * numTriangles), adjacencies(3 * numTriangles);
    std::copy(mesh.GetIndices(), mesh.GetIndices() + 3 * numTriangles, indices.begin());
    std::copy(mesh.GetAdjacencies(), mesh.GetAdjacencies() + 3 * numTriangles,
        adjacencies.begin());
    auto order = CanonicalOrder<3>(indices, numTriangles);
    auto rank = RankOf(order);
    EmitFeatures<3>(io, indices, adjacencies, order, rank);

    io.outInt(mesh.GetInvalidIndex() == Delaunay2<double>::negOne ? -1 : 0);

    Vector2<double> q = RawMeshQuery2(io, qmode, pts, del, order);
    io.givenVec<2>(q);
    size_t found = mesh.GetContainingTriangle(q);
    io.outInt(found == mesh.GetInvalidIndex() ? -1 : static_cast<int32_t>(rank[found]));

    int32_t r = io.integer(0, static_cast<int32_t>(numTriangles) - 1);
    size_t t = order[static_cast<size_t>(r)];

    std::array<Vector2<double>, 3> triVertices{};
    bool vOk = mesh.GetVertices(t, triVertices);
    io.outBool(vOk);
    for (size_t j = 0; j < 3; ++j)
    {
        io.outVec(triVertices[j]);
    }

    std::array<int32_t, 3> triIndices{ 0, 0, 0 };
    bool iOk = mesh.GetIndices(t, triIndices);
    io.outBool(iOk);
    for (size_t j = 0; j < 3; ++j)
    {
        io.outInt(triIndices[j]);
    }

    std::array<int32_t, 3> triAdjacencies{ 0, 0, 0 };
    bool aOk = mesh.GetAdjacencies(t, triAdjacencies);
    io.outBool(aOk);
    for (size_t j = 0; j < 3; ++j)
    {
        io.outInt(triAdjacencies[j] < 0 ? -1
            : static_cast<int32_t>(rank[static_cast<size_t>(triAdjacencies[j])]));
    }

    std::array<double, 3> bary{ 0.0, 0.0, 0.0 };
    bool bOk = mesh.GetBarycentrics(t, q, bary);
    io.outBool(bOk);
    for (size_t j = 0; j < 3; ++j)
    {
        io.outReal(bary[j]);
    }

    // The out-of-range branch of every per-triangle accessor.
    std::array<Vector2<double>, 3> badVertices{};
    bool badV = mesh.GetVertices(numTriangles, badVertices);
    io.outBool(badV);
    for (size_t j = 0; j < 3; ++j)
    {
        io.outVec(badVertices[j]);
    }
    std::array<int32_t, 3> badIndices{ 7, 7, 7 };
    bool badI = mesh.GetIndices(numTriangles, badIndices);
    io.outBool(badI);
    std::array<int32_t, 3> badAdjacencies{ 7, 7, 7 };
    bool badA = mesh.GetAdjacencies(numTriangles, badAdjacencies);
    io.outBool(badA);
    std::array<double, 3> badBary{ 7.0, 7.0, 7.0 };
    bool badB = mesh.GetBarycentrics(numTriangles, q, badBary);
    io.outBool(badB);
    for (size_t j = 0; j < 3; ++j)
    {
        io.outReal(badBary[j]);
    }
}

// Throw parity for the Delaunay2Mesh constructor, which calls
// LogAssert(GetDimension() == 2, "Invalid Delaunay dimension."). Every record
// of this case is a throw record. The inputs are the degenerate sets upstream
// classifies correctly with its hardcoded epsilon of 0 (a repeated point and
// an axis-aligned collinear set), so the dimension itself is not in question.
ORACLE_CASE("Delaunay2Mesh.constructorThrows")
{
    int32_t mode = io.index() % 2;
    int32_t n = io.integer(1, 8);
    std::vector<Vector2<double>> pts(static_cast<size_t>(n));
    for (size_t i = 0; i < pts.size(); ++i)
    {
        double k = (mode == 0 ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
        pts[i][0] = k;
        pts[i][1] = 0.0;
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    Delaunay2<double> del{};
    bool built = del(pts);
    Delaunay2Mesh<double> mesh(del);
    io.outBool(built);
    io.outInt(mesh.GetNumTriangles());
}

// ---- Delaunay3Mesh -------------------------------------------------------

namespace
{
    // The 3D analogue of RawMeshQuery2: only configurations whose containing
    // tetrahedron is unique, because Delaunay3Mesh::GetContainingTetrahedron
    // starts the walk at "tetrahedron 0" of the hash table.
    Vector3<double> RawMeshQuery3(oracle::Ctx& io, int32_t qmode,
        std::vector<Vector3<double>> const& pts, Delaunay3<double> const& del,
        std::vector<size_t> const& order)
    {
        Vector3<double> q{ 0.0, 0.0, 0.0 };
        if (qmode == 0)
        {
            auto const& indices = del.GetIndices();
            int32_t r = io.rawInteger(0, static_cast<int32_t>(order.size()) - 1);
            size_t t = order[static_cast<size_t>(r)];
            Vector3<double> v0 = pts[static_cast<size_t>(indices[4 * t])];
            Vector3<double> v1 = pts[static_cast<size_t>(indices[4 * t + 1])];
            Vector3<double> v2 = pts[static_cast<size_t>(indices[4 * t + 2])];
            Vector3<double> v3 = pts[static_cast<size_t>(indices[4 * t + 3])];
            for (int32_t j = 0; j < 3; ++j)
            {
                q[j] = 0.125 * v0[j] + 0.125 * v1[j] + 0.25 * v2[j] + 0.5 * v3[j];
            }
        }
        else if (qmode == 1)
        {
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
}

// The full public surface of Delaunay3Mesh<T>, mirroring Delaunay2Mesh.query.
ORACLE_CASE("Delaunay3Mesh.query")
{
    int32_t mode = io.index() % 4;
    int32_t qmode = (io.index() / 4) % 3;
    int32_t n = io.integer(4, 10);
    auto pts = SoundPoints3(io, mode, static_cast<size_t>(n));

    Delaunay3<double> del{};
    bool built = del(pts);
    io.outBool(built);

    Delaunay3Mesh<double> mesh(del);
    io.outInt(mesh.GetNumVertices());
    io.outInt(mesh.GetNumTetrahedra());

    size_t numTetrahedra = mesh.GetNumTetrahedra();
    auto const* vertices = mesh.GetVertices();
    for (size_t i = 0; i < static_cast<size_t>(n); ++i)
    {
        io.outVec(vertices[i]);
    }

    std::vector<int32_t> indices(4 * numTetrahedra), adjacencies(4 * numTetrahedra);
    std::copy(mesh.GetIndices(), mesh.GetIndices() + 4 * numTetrahedra, indices.begin());
    std::copy(mesh.GetAdjacencies(), mesh.GetAdjacencies() + 4 * numTetrahedra,
        adjacencies.begin());
    auto order = CanonicalOrder<4>(indices, numTetrahedra);
    auto rank = RankOf(order);
    EmitFeatures<4>(io, indices, adjacencies, order, rank);

    io.outInt(mesh.GetInvalidIndex() == Delaunay3<double>::negOne ? -1 : 0);

    Vector3<double> q = RawMeshQuery3(io, qmode, pts, del, order);
    io.givenVec<3>(q);
    size_t found = mesh.GetContainingTetrahedron(q);
    io.outInt(found == mesh.GetInvalidIndex() ? -1 : static_cast<int32_t>(rank[found]));

    int32_t r = io.integer(0, static_cast<int32_t>(numTetrahedra) - 1);
    size_t t = order[static_cast<size_t>(r)];

    std::array<Vector3<double>, 4> tetVertices{};
    bool vOk = mesh.GetVertices(t, tetVertices);
    io.outBool(vOk);
    for (size_t j = 0; j < 4; ++j)
    {
        io.outVec(tetVertices[j]);
    }

    std::array<int32_t, 4> tetIndices{ 0, 0, 0, 0 };
    bool iOk = mesh.GetIndices(t, tetIndices);
    io.outBool(iOk);
    for (size_t j = 0; j < 4; ++j)
    {
        io.outInt(tetIndices[j]);
    }

    std::array<int32_t, 4> tetAdjacencies{ 0, 0, 0, 0 };
    bool aOk = mesh.GetAdjacencies(t, tetAdjacencies);
    io.outBool(aOk);
    for (size_t j = 0; j < 4; ++j)
    {
        io.outInt(tetAdjacencies[j] < 0 ? -1
            : static_cast<int32_t>(rank[static_cast<size_t>(tetAdjacencies[j])]));
    }

    std::array<double, 4> bary{ 0.0, 0.0, 0.0, 0.0 };
    bool bOk = mesh.GetBarycentrics(t, q, bary);
    io.outBool(bOk);
    for (size_t j = 0; j < 4; ++j)
    {
        io.outReal(bary[j]);
    }

    std::array<Vector3<double>, 4> badVertices{};
    bool badV = mesh.GetVertices(numTetrahedra, badVertices);
    io.outBool(badV);
    for (size_t j = 0; j < 4; ++j)
    {
        io.outVec(badVertices[j]);
    }
    std::array<int32_t, 4> badIndices{ 7, 7, 7, 7 };
    bool badI = mesh.GetIndices(numTetrahedra, badIndices);
    io.outBool(badI);
    std::array<int32_t, 4> badAdjacencies{ 7, 7, 7, 7 };
    bool badA = mesh.GetAdjacencies(numTetrahedra, badAdjacencies);
    io.outBool(badA);
    std::array<double, 4> badBary{ 7.0, 7.0, 7.0, 7.0 };
    bool badB = mesh.GetBarycentrics(numTetrahedra, q, badBary);
    io.outBool(badB);
    for (size_t j = 0; j < 4; ++j)
    {
        io.outReal(badBary[j]);
    }
}

// Throw parity for the Delaunay3Mesh constructor, which calls
// LogAssert(GetDimension() == 3, "Invalid Delaunay dimension."). Every record
// is a throw record. The inputs are the degenerate sets upstream classifies
// correctly with its hardcoded epsilon of 0: a repeated point (dimension 0),
// a collinear set along a coordinate axis (dimension 1) and a coplanar set in
// a plane z = c (dimension 2).
ORACLE_CASE("Delaunay3Mesh.constructorThrows")
{
    int32_t mode = io.index() % 3;
    int32_t n = io.integer(1, 8);
    std::vector<Vector3<double>> pts(static_cast<size_t>(n));
    for (size_t i = 0; i < pts.size(); ++i)
    {
        double k0 = (mode == 0 ? 0.0 : static_cast<double>(io.rawInteger(-3, 3)));
        double k1 = (mode == 2 ? static_cast<double>(io.rawInteger(-3, 3)) : 0.0);
        pts[i][0] = k0;
        pts[i][1] = k1;
        pts[i][2] = 0.0;
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<3>(pts[i]);
    }

    Delaunay3<double> del{};
    bool built = del(pts);
    Delaunay3Mesh<double> mesh(del);
    io.outBool(built);
    io.outInt(mesh.GetNumTetrahedra());
}

// ---- MinimumAreaBox2 -----------------------------------------------------

namespace
{
    // The port instantiates the exact-rational path only, because upstream's
    // own header says a correct output is guaranteed only for an exact
    // ComputeType that supports division and names BSRational<UIntegerAP32>.
    using MAB2 = MinimumAreaBox2<double, ExactQ>;

    int32_t HullDimension2(std::vector<Vector2<double>> const& pts)
    {
        ConvexHull2<double> ch2{};
        ch2(pts);
        return ch2.GetDimension();
    }

    void EmitBox2(oracle::Ctx& io, OrientedBox2<double> const& box)
    {
        io.outVec(box.center);
        io.outVec(box.axis[0]);
        io.outVec(box.axis[1]);
        io.outVec(box.extent);
    }

    // mode 0-3: the shared 2D generators (uniform, lattice, cocircular
    //           lattice, dense lattice)
    // mode 4:   one lattice point repeated n times (hull dimension 0)
    // mode 5:   an axis-aligned square or a 45-degree-rotated square with
    //           extra lattice points mixed in. Every edge of a square gives
    //           a bounding box of the same area, so the minimum is attained
    //           four times and the tie decides the answer. Both sides break
    //           the tie the same way: ComputeBoxForEdgeOrderN seeds minBox
    //           with SmallestBox(n-1, 0) and replaces it only on a strict
    //           'box.area < minBox.area', so the first edge of the rotation
    //           wins; ComputeBoxForEdgeOrderNSqr does the same over the
    //           edges in index order.
    void RawBoxPoints2(oracle::Ctx& io, int32_t mode, size_t n,
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
        if (mode == 4)
        {
            for (size_t i = 0; i < n; ++i)
            {
                pts[i] = { bx, by };
            }
            return;
        }

        double s = static_cast<double>(io.rawInteger(1, 3));
        bool rotated = (io.rawInteger(0, 1) != 0);
        std::array<Vector2<double>, 4> corner{};
        if (rotated)
        {
            corner[0] = { bx + s, by };
            corner[1] = { bx, by + s };
            corner[2] = { bx - s, by };
            corner[3] = { bx, by - s };
        }
        else
        {
            corner[0] = { bx, by };
            corner[1] = { bx + s, by };
            corner[2] = { bx + s, by + s };
            corner[3] = { bx, by + s };
        }
        for (size_t i = 0; i < n; ++i)
        {
            if (i < 4)
            {
                pts[i] = corner[i % 4];
            }
            else
            {
                // Extra points on the square's boundary or inside it, which
                // exercise RemoveCollinearPoints and the support updates.
                int32_t k = io.rawInteger(0, 3);
                double u = 0.5 * static_cast<double>(io.rawInteger(0, 2));
                pts[i][0] = corner[k][0] + u * (corner[(k + 1) % 4][0] - corner[k][0]);
                pts[i][1] = corner[k][1] + u * (corner[(k + 1) % 4][1] - corner[k][1]);
            }
        }
    }

    // Draw a point set whose exact convex-hull dimension is 0 or 2. Hull
    // dimension 1 is excluded from the main cases because upstream's
    // dimension-1 branch reports wrong extreme indices (#328), which the port
    // fixes; see MinimumAreaBox2.deviation.dimension1Extremes.
    std::vector<Vector2<double>> BoxPoints2(oracle::Ctx& io, int32_t mode, size_t n)
    {
        std::vector<Vector2<double>> pts{};
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            RawBoxPoints2(io, mode, n, pts);
            if (HullDimension2(pts) != 1)
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
}

// MinimumAreaBox2 overloads 1 and 2 (the arbitrary-point path), with both
// values of useRotatingCalipers so that ComputeBoxForEdgeOrderN and
// ComputeBoxForEdgeOrderNSqr are both covered. The functor is fresh for every
// record, so upstream's stale-state defect (#328) cannot fire here; it has
// its own deviation case. RemoveCollinearPoints is a no-op on this path
// because ConvexHull2 removes duplicates and guarantees no three consecutive
// collinear hull points, so the port's #286-pattern fix is inert and the two
// implementations are the same computation.
ORACLE_CASE("MinimumAreaBox2.compute")
{
    int32_t mode = io.index() % 6;
    bool useRotatingCalipers = io.boolean();
    int32_t n = io.integer(1, 10);
    auto pts = BoxPoints2(io, mode, static_cast<size_t>(n));

    MAB2 mab{};
    OrientedBox2<double> box = mab(pts, useRotatingCalipers);
    EmitBox2(io, box);
    io.outReal(mab.GetArea());
    auto const& support = mab.GetSupportIndices();
    for (size_t i = 0; i < 4; ++i)
    {
        io.outInt(support[i]);
    }
    auto const& hull = mab.GetHull();
    io.outInt(hull.size());
    for (auto h : hull)
    {
        io.outInt(h);
    }
    io.outInt(mab.GetNumPoints());
    io.outBool(mab.GetPoints() != nullptr);
}

// The hull-dimension-1 branch of overloads 1 and 2. Everything but the
// extreme indices agrees: upstream seeds tmin = tmax = 0 with imin = imax = 0
// and the port seeds them from point 0, but the line origin is one of the
// input points and therefore attains t = 0 exactly, so both compute the same
// tmin and tmax. mHull is not emitted here; it is the subject of
// MinimumAreaBox2.deviation.dimension1Extremes.
ORACLE_CASE("MinimumAreaBox2.compute.dimension1")
{
    bool useRotatingCalipers = io.boolean();
    int32_t n = io.integer(2, 8);
    std::vector<Vector2<double>> pts(static_cast<size_t>(n));
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        int32_t dx = 0, dy = 0;
        for (int32_t k = 0; k < 32 && dx == 0 && dy == 0; ++k)
        {
            dx = io.rawInteger(-2, 2);
            dy = io.rawInteger(-2, 2);
        }
        if (dx == 0 && dy == 0)
        {
            dx = 1;
            dy = 2;
        }
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double k = static_cast<double>(io.rawInteger(-3, 3));
            pts[i][0] = bx + k * static_cast<double>(dx);
            pts[i][1] = by + k * static_cast<double>(dy);
        }
        accepted = (HullDimension2(pts) == 1);
    }
    if (!accepted)
    {
        for (size_t i = 0; i < pts.size(); ++i)
        {
            pts[i][0] = static_cast<double>(i);
            pts[i][1] = 2.0 * static_cast<double>(i);
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    MAB2 mab{};
    OrientedBox2<double> box = mab(pts, useRotatingCalipers);
    EmitBox2(io, box);
    io.outReal(mab.GetArea());
    auto const& support = mab.GetSupportIndices();
    for (size_t i = 0; i < 4; ++i)
    {
        io.outInt(support[i]);
    }
    io.outInt(mab.GetHull().size());
    io.outInt(mab.GetNumPoints());
}

namespace
{
    // Record an integer that was not produced by io.integer().
    void GivenInt(oracle::Ctx& io, int32_t v)
    {
        io.given(static_cast<double>(v));
    }

    // Build a counterclockwise, nondegenerate convex polygon as the convex
    // hull of a lattice point set. ConvexHull2 orders the hull
    // counterclockwise and guarantees no three consecutive collinear hull
    // points, which is exactly overload 3's precondition.
    bool RawConvexPolygon(oracle::Ctx& io, size_t n,
        std::vector<Vector2<double>>& base, std::vector<int32_t>& hull)
    {
        RawPoints2(io, 1, n, base);
        ConvexHull2<double> ch2{};
        ch2(base);
        if (ch2.GetDimension() != 2)
        {
            return false;
        }
        hull = ch2.GetHull();
        return hull.size() >= 3;
    }
}

// MinimumAreaBox2 overload 3, the caller-supplied convex polygon. The four
// submodes are: the points themselves are the polygon (indices == nullptr);
// the polygon is an index subset of the points; the polygon carries exactly
// collinear vertices (every hull vertex doubled and every edge midpoint
// inserted, all integer valued) so that RemoveCollinearPoints actually
// removes something; and the two early-return guards (numPoints < 3 and
// numIndices < 3).
//
// The polygon never has duplicate vertices, so the port's #286-pattern fix of
// RemoveCollinearPoints is inert and the two implementations are the same
// computation; the fix is demonstrated by
// MinimumAreaBox2.deviation.removeCollinear. GetNumPoints() and GetPoints()
// are not emitted here; they are the subject of
// MinimumAreaBox2.deviation.polygonPoints (#402).
ORACLE_CASE("MinimumAreaBox2.computeConvexPolygon")
{
    int32_t submode = io.index() % 4;
    bool useRotatingCalipers = io.boolean();

    std::vector<Vector2<double>> points{};
    std::vector<int32_t> indices{};

    if (submode == 3)
    {
        bool byCount = (io.rawInteger(0, 1) != 0);
        int32_t n = (byCount ? io.rawInteger(1, 2) : io.rawInteger(3, 5));
        points.resize(static_cast<size_t>(n));
        for (size_t i = 0; i < points.size(); ++i)
        {
            points[i][0] = static_cast<double>(io.rawInteger(-3, 3));
            points[i][1] = static_cast<double>(io.rawInteger(-3, 3));
        }
        if (!byCount)
        {
            int32_t numIndices = io.rawInteger(1, 2);
            indices.resize(static_cast<size_t>(numIndices));
            for (size_t i = 0; i < indices.size(); ++i)
            {
                indices[i] = io.rawInteger(0, n - 1);
            }
        }
    }
    else
    {
        std::vector<Vector2<double>> base{};
        std::vector<int32_t> hull{};
        bool accepted = false;
        for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
        {
            base.clear();
            hull.clear();
            size_t n = static_cast<size_t>(io.rawInteger(3, 9));
            accepted = RawConvexPolygon(io, n, base, hull);
        }
        if (!accepted)
        {
            base = { { 0.0, 0.0 }, { 3.0, 0.0 }, { 3.0, 2.0 }, { 0.0, 2.0 } };
            hull = { 0, 1, 2, 3 };
        }

        if (submode == 0)
        {
            points.resize(hull.size());
            for (size_t i = 0; i < hull.size(); ++i)
            {
                points[i] = base[static_cast<size_t>(hull[i])];
            }
        }
        else if (submode == 1)
        {
            points = base;
            indices = hull;
        }
        else
        {
            points.resize(2 * hull.size());
            for (size_t i = 0; i < hull.size(); ++i)
            {
                size_t j = (i + 1) % hull.size();
                Vector2<double> const& p = base[static_cast<size_t>(hull[i])];
                Vector2<double> const& q = base[static_cast<size_t>(hull[j])];
                points[2 * i] = { 2.0 * p[0], 2.0 * p[1] };
                points[2 * i + 1] = { p[0] + q[0], p[1] + q[1] };
            }
        }
    }

    GivenInt(io, static_cast<int32_t>(points.size()));
    for (size_t i = 0; i < points.size(); ++i)
    {
        io.givenVec<2>(points[i]);
    }
    GivenInt(io, static_cast<int32_t>(indices.size()));
    for (size_t i = 0; i < indices.size(); ++i)
    {
        GivenInt(io, indices[i]);
    }

    MAB2 mab{};
    OrientedBox2<double> box = (indices.empty()
        ? mab(static_cast<int32_t>(points.size()), points.data(), 0, nullptr,
            useRotatingCalipers)
        : mab(static_cast<int32_t>(points.size()), points.data(),
            static_cast<int32_t>(indices.size()), indices.data(),
            useRotatingCalipers));
    EmitBox2(io, box);
    io.outReal(mab.GetArea());
    auto const& support = mab.GetSupportIndices();
    for (size_t i = 0; i < 4; ++i)
    {
        io.outInt(support[i]);
    }
    auto const& hull = mab.GetHull();
    io.outInt(hull.size());
    for (auto h : hull)
    {
        io.outInt(h);
    }
}

namespace
{
    // The generator shared by the dimension-1 cases: exactly collinear
    // lattice points whose hull dimension really is 1.
    std::vector<Vector2<double>> Collinear2(oracle::Ctx& io, size_t n)
    {
        std::vector<Vector2<double>> pts(n);
        bool accepted = false;
        for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
        {
            double bx = static_cast<double>(io.rawInteger(-3, 3));
            double by = static_cast<double>(io.rawInteger(-3, 3));
            int32_t dx = 0, dy = 0;
            for (int32_t k = 0; k < 32 && dx == 0 && dy == 0; ++k)
            {
                dx = io.rawInteger(-2, 2);
                dy = io.rawInteger(-2, 2);
            }
            if (dx == 0 && dy == 0)
            {
                dx = 1;
                dy = 2;
            }
            for (size_t i = 0; i < pts.size(); ++i)
            {
                double k = static_cast<double>(io.rawInteger(-3, 3));
                pts[i][0] = bx + k * static_cast<double>(dx);
                pts[i][1] = by + k * static_cast<double>(dy);
            }
            accepted = (HullDimension2(pts) == 1);
        }
        if (!accepted)
        {
            for (size_t i = 0; i < pts.size(); ++i)
            {
                pts[i][0] = static_cast<double>(i);
                pts[i][1] = 2.0 * static_cast<double>(i);
            }
        }
        return pts;
    }
}

// DELIBERATE DEVIATION (#328). On the hull-dimension-1 path, upstream starts
// the t-extremes at tmin = tmax = 0 with imin = imax = 0 "because we know
// that 'origin' is an input vertex". The line origin is points[hull[0]], not
// points[0], so the reported extreme indices are wrong unless hull[0] == 0:
// GetHull() then names an interior point as an extreme. The port seeds both
// extremes from point 0, which gives the same tmin and tmax (the origin is
// still an input point, so t = 0 is attained) and the correct indices.
ORACLE_CASE("MinimumAreaBox2.deviation.dimension1Extremes")
{
    bool useRotatingCalipers = io.boolean();
    int32_t n = io.integer(3, 8);
    auto pts = Collinear2(io, static_cast<size_t>(n));
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    MAB2 mab{};
    OrientedBox2<double> box = mab(pts, useRotatingCalipers);
    io.outReal(box.extent[0]);
    auto const& hull = mab.GetHull();
    io.outInt(hull.size());
    for (auto h : hull)
    {
        io.outInt(h);
    }
}

// DELIBERATE DEVIATION (#328). The degenerate branches of overloads 1 and 2
// return before assigning mArea and mSupportIndices, so a reused functor
// reports the previous data set's area and support indices for a data set
// whose box is a point. The port resets both on every query. The first query
// is a nondegenerate lattice set and the second is a single repeated point,
// whose hull dimension is 0 on both sides.
ORACLE_CASE("MinimumAreaBox2.deviation.staleState")
{
    bool useRotatingCalipers = io.boolean();
    int32_t n = io.integer(4, 9);
    std::vector<Vector2<double>> first{};
    for (int32_t attempt = 0; attempt < 64; ++attempt)
    {
        RawPoints2(io, 1, static_cast<size_t>(n), first);
        if (HullDimension2(first) == 2)
        {
            break;
        }
        first.clear();
    }
    if (first.empty())
    {
        Fallback2(static_cast<size_t>(n), first);
    }
    for (size_t i = 0; i < first.size(); ++i)
    {
        io.givenVec<2>(first[i]);
    }
    double px = io.lattice(-3, 3);
    double py = io.lattice(-3, 3);
    std::vector<Vector2<double>> second(2, Vector2<double>{ px, py });

    MAB2 mab{};
    OrientedBox2<double> box0 = mab(first, useRotatingCalipers);
    io.outReal(box0.extent[0]);
    OrientedBox2<double> box1 = mab(second, useRotatingCalipers);
    EmitBox2(io, box1);
    io.outReal(mab.GetArea());
    auto const& support = mab.GetSupportIndices();
    for (size_t i = 0; i < 4; ++i)
    {
        io.outInt(support[i]);
    }
}

// DELIBERATE DEVIATION (#402). Overloads 1 and 2 begin with
// 'mNumPoints = numPoints; mPoints = points;'; overload 3 only clears mHull,
// so after a caller-supplied polygon GetPoints() is null on a fresh functor
// (or the previous data set's pointer on a reused one) while GetHull() and
// GetSupportIndices() refer to the polygon just processed, and the documented
// mPoints[hull[...]] lookup dereferences null. The port assigns both.
ORACLE_CASE("MinimumAreaBox2.deviation.polygonPoints")
{
    bool useRotatingCalipers = io.boolean();
    std::vector<Vector2<double>> base{};
    std::vector<int32_t> hull{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        base.clear();
        hull.clear();
        size_t n = static_cast<size_t>(io.rawInteger(3, 9));
        accepted = RawConvexPolygon(io, n, base, hull);
    }
    if (!accepted)
    {
        base = { { 0.0, 0.0 }, { 3.0, 0.0 }, { 3.0, 2.0 }, { 0.0, 2.0 } };
        hull = { 0, 1, 2, 3 };
    }
    std::vector<Vector2<double>> polygon(hull.size());
    for (size_t i = 0; i < hull.size(); ++i)
    {
        polygon[i] = base[static_cast<size_t>(hull[i])];
    }
    GivenInt(io, static_cast<int32_t>(polygon.size()));
    for (size_t i = 0; i < polygon.size(); ++i)
    {
        io.givenVec<2>(polygon[i]);
    }

    MAB2 mab{};
    OrientedBox2<double> box = mab(static_cast<int32_t>(polygon.size()),
        polygon.data(), 0, nullptr, useRotatingCalipers);
    io.outReal(box.extent[0]);
    io.outInt(mab.GetNumPoints());
    io.outBool(mab.GetPoints() != nullptr);
}

// DELIBERATE DEVIATION (#286 pattern, recorded for MinimumAreaBox2 under
// #328). RemoveCollinearPoints tests collinearity against the immediately
// preceding edge of the input array, which is the zero-length edge of a
// duplicated polygon vertex; DotPerp is then zero and the genuine corner is
// discarded along with its duplicate. The port compares against the most
// recent nonzero edge, which is identical whenever the polygon has no
// duplicates. The generator duplicates one vertex of a convex polygon with at
// least 4 corners, so upstream's reduced polygon still has at least 3
// vertices and the comparison is box against box rather than throw against
// garbage.
ORACLE_CASE("MinimumAreaBox2.deviation.removeCollinear")
{
    bool useRotatingCalipers = io.boolean();
    std::vector<Vector2<double>> base{};
    std::vector<int32_t> hull{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        base.clear();
        hull.clear();
        size_t n = static_cast<size_t>(io.rawInteger(4, 9));
        accepted = RawConvexPolygon(io, n, base, hull) && hull.size() >= 4;
    }
    if (!accepted)
    {
        base = { { 0.0, 0.0 }, { 3.0, 0.0 }, { 3.0, 2.0 }, { 0.0, 2.0 } };
        hull = { 0, 1, 2, 3 };
    }
    int32_t d = io.rawInteger(0, static_cast<int32_t>(hull.size()) - 1);

    std::vector<Vector2<double>> polygon{};
    for (size_t i = 0; i < hull.size(); ++i)
    {
        Vector2<double> const& p = base[static_cast<size_t>(hull[i])];
        if (static_cast<int32_t>(i) == d)
        {
            polygon.push_back(p);
        }
        polygon.push_back(p);
    }
    GivenInt(io, static_cast<int32_t>(polygon.size()));
    for (size_t i = 0; i < polygon.size(); ++i)
    {
        io.givenVec<2>(polygon[i]);
    }

    MAB2 mab{};
    OrientedBox2<double> box = mab(static_cast<int32_t>(polygon.size()),
        polygon.data(), 0, nullptr, useRotatingCalipers);
    EmitBox2(io, box);
    io.outReal(mab.GetArea());
    auto const& support = mab.GetSupportIndices();
    for (size_t i = 0; i < 4; ++i)
    {
        io.outInt(support[i]);
    }
}

// DESIGN DEVIATION (port note in src/MinimumAreaBox2.ts). Upstream's header
// says a correct output is guaranteed only for an exact ComputeType and
// names BSRational<UIntegerAP32>; the port implements only that path. This
// case runs upstream's floating-point instantiation
// MinimumAreaBox2<double, double> (whose useRotatingCalipers default is
// false) against the port's exact computation on the same uniform points, to
// show what the choice changes: the box is rounded once at the end on the
// port and at every step on upstream's double path. It is not a defect fix
// and the main cases compare the matching exact instantiation.
ORACLE_CASE("MinimumAreaBox2.deviation.floatComputeType")
{
    int32_t n = io.integer(3, 10);
    std::vector<Vector2<double>> pts{};
    for (int32_t attempt = 0; attempt < 64; ++attempt)
    {
        RawPoints2(io, 0, static_cast<size_t>(n), pts);
        if (HullDimension2(pts) == 2)
        {
            break;
        }
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    MinimumAreaBox2<double, double> mab{};
    OrientedBox2<double> box = mab(pts);
    EmitBox2(io, box);
    io.outReal(mab.GetArea());
}

// ---- MinimumWidthPoints2 -------------------------------------------------

// MinimumWidthPoints2<T> overloads 1 and 2, with both values of
// useRotatingCalipers so that the rotating-calipers branch (which minimizes
// an exact BSRational<UIntegerAP32> squared width over the antipodes) and the
// brute-force branch (which minimizes a floating-point width over the edges)
// are both covered. All three hull dimensions are reachable; the degenerate
// branch's opposite frame handedness is an upstream quirk the port preserves
// verbatim, so it is compared here rather than deviating.
//
// The brute-force branch's own duplicate-point corner drop is unreachable
// (its input is always a ConvexHull2 hull, which has no duplicates), so the
// port's fix of it is inert and the two implementations are the same
// computation; the deviation is demonstrated on MinimumAreaBox2, which shares
// the defect and does reach it through the caller-supplied-polygon path.
//
// Ties: the antipode search keeps the first minimum ('sqrWidth < minSqrWidth'
// over the antipode array, which RotatingCalipers builds deterministically)
// and the brute-force search keeps the first minimum over the edges in index
// order, so a square or another tie-rich lattice set resolves identically on
// both sides.
ORACLE_CASE("MinimumWidthPoints2.compute")
{
    int32_t mode = io.index() % 6;
    bool useRotatingCalipers = io.boolean();
    int32_t n = io.integer(3, 10);
    std::vector<Vector2<double>> pts{};
    RawBoxPoints2(io, mode, static_cast<size_t>(n), pts);
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    MinimumWidthPoints2<double> mwp{};
    OrientedBox2<double> box = mwp(pts, useRotatingCalipers);
    EmitBox2(io, box);
}

// MinimumWidthPoints2<T> overload 3 with a nonempty index subset, which
// compacts the selected points and forwards to overload 2. The selection is
// an arbitrary index list (upstream's precondition is only numIndices >= 3;
// it recomputes the convex hull of the compacted points anyway), so duplicate
// and collinear selections are included.
ORACLE_CASE("MinimumWidthPoints2.computeIndexed")
{
    int32_t mode = io.index() % 6;
    bool useRotatingCalipers = io.boolean();
    int32_t n = io.integer(3, 10);
    std::vector<Vector2<double>> pts{};
    RawBoxPoints2(io, mode, static_cast<size_t>(n), pts);
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }
    int32_t numIndices = io.integer(3, 10);
    std::vector<int32_t> indices(static_cast<size_t>(numIndices));
    for (size_t i = 0; i < indices.size(); ++i)
    {
        indices[i] = io.integer(0, n - 1);
    }

    MinimumWidthPoints2<double> mwp{};
    OrientedBox2<double> box = mwp(pts, indices, useRotatingCalipers);
    EmitBox2(io, box);
}

// Throw parity for MinimumWidthPoints2's input validation:
// LogAssert(numPoints >= 3 && points != nullptr) in overloads 1 and 2, and
// LogAssert(numPoints >= 3 && ((indices == nullptr && numIndices == 0) ||
// (indices != nullptr && numIndices >= 3))) in overload 3. Every record is a
// throw record; the generator alternates between too few points and too few
// indices.
ORACLE_CASE("MinimumWidthPoints2.invalidInputThrows")
{
    bool byCount = (io.index() % 2 == 0);
    int32_t n = (byCount ? io.integer(1, 2) : io.integer(3, 6));
    std::vector<Vector2<double>> pts(static_cast<size_t>(n));
    for (size_t i = 0; i < pts.size(); ++i)
    {
        pts[i][0] = io.lattice(-3, 3);
        pts[i][1] = io.lattice(-3, 3);
    }
    // Recorded on both branches so that the record layout is fixed.
    int32_t numIndices = (byCount ? io.integer(0, 0) : io.integer(1, 2));
    std::vector<int32_t> indices(static_cast<size_t>(numIndices));
    for (size_t i = 0; i < indices.size(); ++i)
    {
        indices[i] = io.integer(0, n - 1);
    }

    MinimumWidthPoints2<double> mwp{};
    OrientedBox2<double> box = (byCount
        ? mwp(static_cast<int32_t>(pts.size()), pts.data(), true)
        : mwp(static_cast<int32_t>(pts.size()), pts.data(),
            static_cast<int32_t>(indices.size()), indices.data(), true));
    EmitBox2(io, box);
}

// ---- SeparatePoints2 -----------------------------------------------------

namespace
{
    // A verbatim copy of upstream's side classification for one candidate
    // edge: lineNormal = Perp(Normalize(P1 - P0)),
    // lineConstant = Dot(lineNormal, P0), and the comparison of
    // Dot(lineNormal, Q) against lineConstant. Both OnSameSide and WhichSide
    // use exactly this expression.
    int32_t UpstreamSide(Vector2<double> const& p0, Vector2<double> const& p1,
        Vector2<double> const& q)
    {
        Vector2<double> direction = p1 - p0;
        Normalize(direction);
        Vector2<double> lineNormal = Perp(direction);
        double lineConstant = Dot(lineNormal, p0);
        double c0 = Dot(lineNormal, q);
        return (c0 > lineConstant ? +1 : (c0 < lineConstant ? -1 : 0));
    }

    // The exact sign of the same quantity. Perp(P1 - P0) = (dy, -dx), so
    // Dot(Perp(P1-P0), Q - P0) = -DotPerp(P1 - P0, Q - P0), and normalizing
    // the direction only scales by a positive factor.
    int32_t ExactSide(Vector2<double> const& p0, Vector2<double> const& p1,
        Vector2<double> const& q)
    {
        return -ExactOrient2(p0, p1, q);
    }

    std::vector<Vector2<double>> HullPoints(std::vector<Vector2<double>> const& pts,
        std::vector<int32_t> const& edges)
    {
        std::vector<Vector2<double>> hull(edges.size());
        for (size_t i = 0; i < edges.size(); ++i)
        {
            hull[i] = pts[static_cast<size_t>(edges[i])];
        }
        return hull;
    }

    // Accept a record only when every floating-point side classification the
    // query can perform agrees with the exact one. This is an exact superset
    // of the inputs on which upstream's control flow equals the port's: when
    // all per-point signs agree, OnSameSide's return value agrees (it depends
    // only on whether both a positive and a negative sign occur, which is
    // order independent) and WhichSide's agrees (it is always called with the
    // hull that owns the edge, whose vertices are all on the non-positive
    // side in exact arithmetic, so the first nonzero sign is -1 in either
    // visiting order). A point set whose hull is 0- or 1-dimensional is
    // accepted unconditionally: upstream returns false before any side test.
    bool SeparateSound(std::vector<Vector2<double>> const& pts0,
        std::vector<Vector2<double>> const& pts1)
    {
        ConvexHull2<double> ch0{};
        ch0(pts0);
        if (ch0.GetDimension() != 2)
        {
            return true;
        }
        ConvexHull2<double> ch1{};
        ch1(pts1);
        if (ch1.GetDimension() != 2)
        {
            return true;
        }

        std::vector<Vector2<double>> hull0 = HullPoints(pts0, ch0.GetHull());
        std::vector<Vector2<double>> hull1 = HullPoints(pts1, ch1.GetHull());
        std::array<std::vector<Vector2<double>> const*, 2> hulls{ &hull0, &hull1 };
        for (size_t h = 0; h < 2; ++h)
        {
            auto const& owner = *hulls[h];
            size_t numEdges = owner.size();
            for (size_t j1 = 0, j0 = numEdges - 1; j1 < numEdges; j0 = j1++)
            {
                for (size_t k = 0; k < 2; ++k)
                {
                    for (auto const& q : *hulls[k])
                    {
                        if (UpstreamSide(owner[j0], owner[j1], q)
                            != ExactSide(owner[j0], owner[j1], q))
                        {
                            return false;
                        }
                    }
                }
            }
        }
        return true;
    }

    // The port's algorithm, evaluated with the exact predicate. Used only to
    // classify a draw for the deviation case.
    bool ExactSeparate(std::vector<Vector2<double>> const& pts0,
        std::vector<Vector2<double>> const& pts1)
    {
        ConvexHull2<double> ch0{};
        ch0(pts0);
        if (ch0.GetDimension() != 2)
        {
            return false;
        }
        ConvexHull2<double> ch1{};
        ch1(pts1);
        if (ch1.GetDimension() != 2)
        {
            return false;
        }

        std::vector<Vector2<double>> hull0 = HullPoints(pts0, ch0.GetHull());
        std::vector<Vector2<double>> hull1 = HullPoints(pts1, ch1.GetHull());
        std::array<std::vector<Vector2<double>> const*, 2> hulls{ &hull0, &hull1 };
        for (size_t h = 0; h < 2; ++h)
        {
            auto const& owner = *hulls[h];
            auto const& other = *hulls[1 - h];
            size_t numEdges = owner.size();
            for (size_t j1 = 0, j0 = numEdges - 1; j1 < numEdges; j0 = j1++)
            {
                int32_t posSide = 0, negSide = 0;
                for (auto const& q : other)
                {
                    int32_t s = ExactSide(owner[j0], owner[j1], q);
                    if (s > 0)
                    {
                        ++posSide;
                    }
                    else if (s < 0)
                    {
                        ++negSide;
                    }
                }
                int32_t sideOther = (posSide != 0 && negSide != 0)
                    ? 0 : (posSide != 0 ? +1 : -1);
                if (sideOther != 0)
                {
                    int32_t sideOwner = 0;
                    for (auto const& q : owner)
                    {
                        int32_t s = ExactSide(owner[j0], owner[j1], q);
                        if (s != 0)
                        {
                            sideOwner = s;
                            break;
                        }
                    }
                    if (sideOwner * sideOther <= 0)
                    {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    // A regular k-gon. The trigonometry runs only in the generator; the
    // vertices themselves are recorded, so the replay never calls libm.
    void RawPolygon2(oracle::Ctx& io, size_t k, double cx, double cy,
        double radius, double phase, std::vector<Vector2<double>>& pts)
    {
        (void)io;
        pts.resize(k);
        for (size_t i = 0; i < k; ++i)
        {
            double angle = phase + 6.283185307179586 * static_cast<double>(i)
                / static_cast<double>(k);
            pts[i][0] = cx + radius * std::cos(angle);
            pts[i][1] = cy + radius * std::sin(angle);
        }
    }

    // mode 0: two lattice sets
    // mode 1: two uniform sets
    // mode 2: two regular polygons, overlapping or disjoint
    // mode 3: one of the sets is degenerate (collinear or a repeated point),
    //         which is the 'GetDimension() != 2' early return
    void RawSeparatePair(oracle::Ctx& io, int32_t mode, size_t n0, size_t n1,
        std::vector<Vector2<double>>& pts0, std::vector<Vector2<double>>& pts1)
    {
        if (mode == 0 || mode == 1)
        {
            RawPoints2(io, mode == 0 ? 1 : 0, n0, pts0);
            RawPoints2(io, mode == 0 ? 1 : 0, n1, pts1);
            double shiftX = static_cast<double>(io.rawInteger(-6, 6));
            double shiftY = static_cast<double>(io.rawInteger(-6, 6));
            for (size_t i = 0; i < pts1.size(); ++i)
            {
                pts1[i][0] += shiftX;
                pts1[i][1] += shiftY;
            }
            return;
        }

        if (mode == 2)
        {
            double c0x = io.raw(-2.0, 2.0);
            double c0y = io.raw(-2.0, 2.0);
            double r0 = io.raw(0.5, 2.0);
            double a0 = io.raw(0.0, 6.283185307179586);
            RawPolygon2(io, n0, c0x, c0y, r0, a0, pts0);
            double c1x = io.raw(-3.0, 3.0);
            double c1y = io.raw(-3.0, 3.0);
            double r1 = io.raw(0.5, 2.0);
            double a1 = io.raw(0.0, 6.283185307179586);
            RawPolygon2(io, n1, c1x, c1y, r1, a1, pts1);
            return;
        }

        RawPoints2(io, 1, n0, pts0);
        double bx = static_cast<double>(io.rawInteger(-3, 3));
        double by = static_cast<double>(io.rawInteger(-3, 3));
        bool collinear = (io.rawInteger(0, 1) != 0);
        pts1.resize(n1);
        for (size_t i = 0; i < n1; ++i)
        {
            double k = (collinear ? static_cast<double>(io.rawInteger(-3, 3)) : 0.0);
            pts1[i][0] = bx + k;
            pts1[i][1] = by;
        }
    }
}

// SeparatePoints2<T> on inputs where every floating-point side classification
// upstream performs agrees with the exact one (see SeparateSound), so the two
// implementations follow the same control flow and return the same line. The
// separating line is emitted only when the query returns true: upstream
// writes 'separatingLine' for every candidate edge it tests and therefore
// leaves the caller's line holding the last candidate on a false return,
// while the port returns a default-constructed line there; the upstream
// documentation promises the field only on a true return.
ORACLE_CASE("SeparatePoints2.compute")
{
    int32_t mode = io.index() % 4;
    // Sets of one or two points are always degenerate (the early 'false'
    // return); they are drawn only in mode 3, the degenerate mode, so that
    // modes 0-2 spend their records on the side tests.
    int32_t minCount = (mode == 3 ? 1 : 3);
    int32_t n0 = io.integer(minCount, 8);
    int32_t n1 = io.integer(minCount, 8);
    std::vector<Vector2<double>> pts0{}, pts1{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        pts0.clear();
        pts1.clear();
        RawSeparatePair(io, mode, static_cast<size_t>(n0),
            static_cast<size_t>(n1), pts0, pts1);
        accepted = SeparateSound(pts0, pts1);
    }
    if (!accepted)
    {
        // Two axis-aligned lattice squares far apart: every side test is on
        // exactly representable integers scaled by an axis-aligned unit
        // normal, so no classification is in doubt.
        pts0.assign(static_cast<size_t>(n0), Vector2<double>{ 0.0, 0.0 });
        pts1.assign(static_cast<size_t>(n1), Vector2<double>{ 10.0, 0.0 });
        for (size_t i = 0; i < pts0.size(); ++i)
        {
            pts0[i][0] = static_cast<double>(i % 2);
            pts0[i][1] = static_cast<double>((i / 2) % 2);
        }
        for (size_t i = 0; i < pts1.size(); ++i)
        {
            pts1[i][0] = 10.0 + static_cast<double>(i % 2);
            pts1[i][1] = static_cast<double>((i / 2) % 2);
        }
    }
    for (size_t i = 0; i < pts0.size(); ++i)
    {
        io.givenVec<2>(pts0[i]);
    }
    for (size_t i = 0; i < pts1.size(); ++i)
    {
        io.givenVec<2>(pts1[i]);
    }

    Line2<double> separatingLine{};
    SeparatePoints2<double, double> query{};
    bool separated = query(static_cast<int32_t>(pts0.size()), pts0.data(),
        static_cast<int32_t>(pts1.size()), pts1.data(), separatingLine);
    io.outBool(separated);
    if (separated)
    {
        io.outVec(separatingLine.origin);
        io.outVec(separatingLine.direction);
    }
}

// DELIBERATE DEVIATION (#328). OnSameSide and WhichSide compare
// Dot(Perp(Normalize(P1 - P0)), Q) against the rounded constant
// Dot(Perp(Normalize(P1 - P0)), P0), so the far endpoint P1 of the hull's own
// candidate edge can classify as strictly positive by round-off. WhichSide
// visits the hull vertices starting at the last one, so for the edge
// <hull[n-1], hull[0]> it sees P0 (exactly on the line) and then P1, and a
// round-off there makes it return +1 instead of -1; 'side0 * side1 <= 0' then
// succeeds and two heavily overlapping point sets are reported as separated.
// The port evaluates the side with an exact orientation predicate on the
// unnormalized edge normal, which has the same sign in exact arithmetic.
//
// The generator draws pairs of regular polygons, keeps only those that
// overlap in exact arithmetic, and accepts a draw only when upstream reports
// a separation - that is, only genuinely wrong upstream results.
ORACLE_CASE("SeparatePoints2.deviation.roundoff")
{
    int32_t n0 = io.integer(5, 8);
    int32_t n1 = io.integer(5, 8);
    std::vector<Vector2<double>> pts0{}, pts1{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 512 && !accepted; ++attempt)
    {
        double c0x = io.raw(-0.5, 0.5);
        double c0y = io.raw(-0.5, 0.5);
        double r0 = io.raw(0.8, 1.8);
        double a0 = io.raw(0.0, 6.283185307179586);
        RawPolygon2(io, static_cast<size_t>(n0), c0x, c0y, r0, a0, pts0);
        double c1x = io.raw(-1.2, 1.2);
        double c1y = io.raw(-1.2, 1.2);
        double r1 = io.raw(0.8, 1.8);
        double a1 = io.raw(0.0, 6.283185307179586);
        RawPolygon2(io, static_cast<size_t>(n1), c1x, c1y, r1, a1, pts1);
        if (ExactSeparate(pts0, pts1))
        {
            continue;
        }
        Line2<double> probeLine{};
        SeparatePoints2<double, double> probe{};
        accepted = probe(n0, pts0.data(), n1, pts1.data(), probeLine);
    }
    if (!accepted)
    {
        // The reproduction recorded in docs/UPSTREAM-FINDINGS.md.
        RawPolygon2(io, static_cast<size_t>(n0), 0.0, 0.0, 1.0, 0.3, pts0);
        RawPolygon2(io, static_cast<size_t>(n1), 1.0, 0.5, 1.5, 1.1, pts1);
    }
    for (size_t i = 0; i < pts0.size(); ++i)
    {
        io.givenVec<2>(pts0[i]);
    }
    for (size_t i = 0; i < pts1.size(); ++i)
    {
        io.givenVec<2>(pts1[i]);
    }

    Line2<double> separatingLine{};
    SeparatePoints2<double, double> query{};
    bool separated = query(static_cast<int32_t>(pts0.size()), pts0.data(),
        static_cast<int32_t>(pts1.size()), pts1.data(), separatingLine);
    io.outBool(separated);
}

// ---- ConstrainedDelaunay2 ------------------------------------------------

namespace
{
    using CDT2 = ConstrainedDelaunay2<double>;

    // The Delaunay edges of a triangulation as sorted (min, max) pairs. Used
    // only by the generators to aim a constraint at an existing edge.
    std::vector<std::array<int32_t, 2>> SortedGraphEdges(Delaunay2<double> const& del)
    {
        std::vector<std::array<int32_t, 2>> edges{};
        for (auto const& element : del.GetGraph().GetEdges())
        {
            edges.push_back({ element.first.V[0], element.first.V[1] });
        }
        std::sort(edges.begin(), edges.end());
        return edges;
    }

    // Draw one constraint edge, unrecorded. The two endpoints never map to
    // the same duplicate representative, so upstream's 'Invalid edge.'
    // assert cannot fire (it has its own throw-parity case). The endpoints
    // are raw input indices: when one is a duplicate, Insert substitutes its
    // representative, and that path is covered too.
    //   emode 0: a uniform pair
    //   emode 1: an edge of the unconstrained Delaunay triangulation (already
    //            present unless an earlier constraint removed it)
    //   emode 2: the vertex farthest from a random vertex, which crosses the
    //            most triangles
    //   emode 3: a chain, starting at the previous constraint's end
    std::array<int32_t, 2> RawConstraint(oracle::Ctx& io,
        std::vector<Vector2<double>> const& pts, std::vector<size_t> const& duplicates,
        std::vector<std::array<int32_t, 2>> const& graphEdges, int32_t previousEnd)
    {
        int32_t n = static_cast<int32_t>(pts.size());
        int32_t emode = io.rawInteger(0, 3);
        std::array<int32_t, 2> edge{ 0, 0 };
        if (emode == 1 && !graphEdges.empty())
        {
            int32_t k = io.rawInteger(0, static_cast<int32_t>(graphEdges.size()) - 1);
            edge = graphEdges[static_cast<size_t>(k)];
            if (io.rawInteger(0, 1) != 0)
            {
                std::swap(edge[0], edge[1]);
            }
            return edge;
        }

        int32_t start = (emode == 3 && previousEnd >= 0 ? previousEnd : io.rawInteger(0, n - 1));
        if (emode == 2)
        {
            double best = -1.0;
            for (int32_t j = 0; j < n; ++j)
            {
                Vector2<double> diff = pts[static_cast<size_t>(j)] - pts[static_cast<size_t>(start)];
                double sqrLength = Dot(diff, diff);
                if (sqrLength > best)
                {
                    best = sqrLength;
                    edge = { start, j };
                }
            }
            return edge;
        }

        // A capped rejection loop over the second endpoint; at least one
        // vertex has a different representative because the set has
        // dimension 2.
        edge = { start, start };
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            int32_t other = io.rawInteger(0, n - 1);
            if (duplicates[static_cast<size_t>(other)] != duplicates[static_cast<size_t>(start)])
            {
                edge[1] = other;
                break;
            }
        }
        if (edge[1] == start)
        {
            for (int32_t j = 0; j < n; ++j)
            {
                if (duplicates[static_cast<size_t>(j)] != duplicates[static_cast<size_t>(start)])
                {
                    edge[1] = j;
                    break;
                }
            }
        }
        return edge;
    }

    // The whole constrained triangulation after the insertions: the inserted
    // edge set (a std::unordered_set upstream, sorted here), then the
    // compact triangle arrays after UpdateIndicesAdjacencies, canonicalized
    // exactly as Delaunay2, the hull edges and the graph sizes.
    void EmitCDT2(oracle::Ctx& io, CDT2& cdt)
    {
        std::vector<int32_t> inserted{};
        for (auto const& ekey : cdt.GetInsertedEdges())
        {
            inserted.push_back(ekey.V[0]);
            inserted.push_back(ekey.V[1]);
        }
        EmitSortedTuples<2, int32_t>(io, inserted);

        cdt.UpdateIndicesAdjacencies();
        size_t numTriangles = cdt.GetNumTriangles();
        auto const& indices = cdt.GetIndices();
        auto const& adjacencies = cdt.GetAdjacencies();
        auto order = CanonicalOrder<3>(indices, numTriangles);
        auto rank = RankOf(order);
        EmitFeatures<3>(io, indices, adjacencies, order, rank);

        std::vector<size_t> hull{};
        bool hullOk = cdt.GetHull(hull);
        io.outBool(hullOk);
        EmitSortedTuples<2, size_t>(io, hull);

        io.outInt(cdt.GetGraph().GetTriangles().size());
        io.outInt(cdt.GetGraph().GetEdges().size());
    }

    void EmitPartition(oracle::Ctx& io, std::vector<int32_t> const& partitionedEdge)
    {
        io.outInt(partitionedEdge.size());
        for (auto v : partitionedEdge)
        {
            io.outInt(v);
        }
    }
}

// ConstrainedDelaunay2<T>: both operator() overloads (alternating by
// record), Insert, GetInsertedEdges and the inherited compact-array
// accessors after UpdateIndicesAdjacencies. The point sets are the shared 2D
// generators restricted by Sound2, because the Delaunay2 base carries the
// deliberate #391 deviation (v09-compgeom). One to four constraints are
// inserted in sequence; a later constraint may cross an earlier one, which
// upstream documents as legal ("the second insertion will interfere with the
// retriangulation of the first edge").
//
// Every predicate on this path is exact (interval filter plus rational
// fallback in ToLine, rational pseudosquared distances in SelectSplit), and
// the strip retriangulation is a deterministic bisection, so the whole
// combinatorial output is compared bit for bit. The strip fill is NOT
// constrained-Delaunay (docs/UPSTREAM-FINDINGS.md, ConstrainedDelaunay2
// finding 3, preserved); nothing here asserts that it is. Order: GetLinkEdges
// iterates VETManifoldMesh::Vertex::TAdjacent, a pointer-keyed set, but at
// most two link triangles pass the 'sign0 >= 0 && sign1 <= 0' test, and only
// when the constraint runs through their shared vertex, where both call
// ProcessCoincidentEdge with that vertex; the result is order independent.
ORACLE_CASE("ConstrainedDelaunay2.insert")
{
    int32_t mode = io.index() % 4;
    int32_t n = io.integer(3, 10);
    auto pts = SoundPoints2(io, mode, static_cast<size_t>(n));

    std::vector<size_t> duplicates{};
    std::vector<std::array<int32_t, 2>> graphEdges{};
    {
        Delaunay2<double> probe{};
        probe(pts);
        duplicates = probe.GetDuplicates();
        graphEdges = SortedGraphEdges(probe);
    }
    int32_t numEdges = io.integer(1, 4);
    std::vector<std::array<int32_t, 2>> edges(static_cast<size_t>(numEdges));
    int32_t previousEnd = -1;
    for (size_t k = 0; k < edges.size(); ++k)
    {
        edges[k] = RawConstraint(io, pts, duplicates, graphEdges, previousEnd);
        previousEnd = edges[k][1];
        GivenInt(io, edges[k][0]);
        GivenInt(io, edges[k][1]);
    }

    CDT2 cdt{};
    bool built = (io.index() % 2 == 0 ? cdt(pts) : cdt(pts.size(), pts.data()));
    io.outBool(built);
    for (auto const& edge : edges)
    {
        std::vector<int32_t> partitionedEdge{};
        cdt.Insert(edge, partitionedEdge);
        EmitPartition(io, partitionedEdge);
    }
    EmitCDT2(io, cdt);
}

namespace
{
    // A long constraint through a band of points. The endpoints are
    // base -/+ 5*d for a lattice direction d; the other points sit at
    // base + s*d + h*perp(d) with |s| <= 3, so the Delaunay triangulation
    // has many triangles across the segment and the strip is long, which is
    // where the bisection retriangulation (SelectSplit) does real work.
    //   smode 0: uniform s and h (h bounded away from 0)
    //   smode 1: lattice s and h, h = 0 on a quarter of the points, so the
    //            strip walk meets vertices on the constraint
    //            (ProcessTriangleStrip's querySign == 0 subdivision)
    //   smode 2: every band point on the segment's line except one or two
    //            apexes, so most of the constraint is coincident with
    //            existing edges (the ProcessCoincidentEdge chain)
    void RawStripPoints(oracle::Ctx& io, int32_t smode, size_t n,
        std::vector<Vector2<double>>& pts)
    {
        pts.resize(n);
        double bx = static_cast<double>(io.rawInteger(-2, 2));
        double by = static_cast<double>(io.rawInteger(-2, 2));
        double dx = static_cast<double>(io.rawInteger(1, 2));
        double dy = static_cast<double>(io.rawInteger(-1, 1));
        pts[0] = { bx - 5.0 * dx, by - 5.0 * dy };
        pts[1] = { bx + 5.0 * dx, by + 5.0 * dy };
        size_t numApexes = static_cast<size_t>(io.rawInteger(1, 2));
        for (size_t i = 2; i < n; ++i)
        {
            double s = 0.0, h = 0.0;
            if (smode == 0)
            {
                s = io.raw(-3.0, 3.0);
                h = io.raw(0.25, 2.0);
                if (io.rawInteger(0, 1) != 0)
                {
                    h = -h;
                }
            }
            else if (smode == 1)
            {
                s = static_cast<double>(io.rawInteger(-3, 3));
                h = (io.rawInteger(0, 3) == 0 ? 0.0 : static_cast<double>(io.rawInteger(-2, 2)));
            }
            else
            {
                s = static_cast<double>(io.rawInteger(-4, 4));
                h = (i < 2 + numApexes ? static_cast<double>(io.rawInteger(1, 3))
                    * (i % 2 == 0 ? 1.0 : -1.0) : 0.0);
            }
            pts[i][0] = bx + s * dx - h * dy;
            pts[i][1] = by + s * dy + h * dx;
        }
    }
}

// Insert the long constraint <0, 1>, then optionally a second constraint
// between the two band points farthest apart across the first one, which
// crosses it (the documented interference), then the first one again.
ORACLE_CASE("ConstrainedDelaunay2.insert.strip")
{
    int32_t smode = io.index() % 3;
    int32_t n = io.integer(5, 11);
    std::vector<Vector2<double>> pts{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        RawStripPoints(io, smode, static_cast<size_t>(n), pts);
        accepted = Sound2(pts);
    }
    if (!accepted)
    {
        Fallback2(static_cast<size_t>(n), pts);
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    // The crossing constraint: the lowest-indexed point strictly left of the
    // line <0,1> and the lowest-indexed point strictly right of it, by the
    // exact predicate. -1 when one side is empty.
    int32_t left = -1, right = -1;
    for (int32_t i = 2; i < n; ++i)
    {
        int32_t sign = ExactOrient2(pts[0], pts[1], pts[static_cast<size_t>(i)]);
        if (sign > 0 && left < 0)
        {
            left = i;
        }
        else if (sign < 0 && right < 0)
        {
            right = i;
        }
    }
    bool cross = (left >= 0 && right >= 0 && io.rawInteger(0, 1) != 0);
    GivenInt(io, cross ? 1 : 0);
    GivenInt(io, cross ? left : 0);
    GivenInt(io, cross ? right : 0);

    CDT2 cdt{};
    bool built = cdt(pts);
    io.outBool(built);
    std::vector<int32_t> partitionedEdge{};
    cdt.Insert({ 0, 1 }, partitionedEdge);
    EmitPartition(io, partitionedEdge);
    if (cross)
    {
        std::vector<int32_t> crossing{};
        cdt.Insert({ left, right }, crossing);
        EmitPartition(io, crossing);
        std::vector<int32_t> again{};
        cdt.Insert({ 1, 0 }, again);
        EmitPartition(io, again);
    }
    EmitCDT2(io, cdt);
}

// Throw parity for Insert's LogAssert(edge[0] != edge[1] && in range,
// "Invalid edge."), evaluated after the duplicate substitution: the two
// endpoints are the same index, or two different indices of the same point.
// Every record is a throw record. Out-of-range indices are not drawn here:
// upstream reads duplicates[] before the range check, which is the subject
// of ConstrainedDelaunay2.deviation.duplicatesRead.
ORACLE_CASE("ConstrainedDelaunay2.insert.invalidEdgeThrows")
{
    int32_t n = io.integer(4, 9);
    std::vector<Vector2<double>> pts{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        RawPoints2(io, 1, static_cast<size_t>(n - 1), pts);
        // Repeat one point so that a duplicate pair exists.
        Vector2<double> repeated = pts[static_cast<size_t>(io.rawInteger(0, n - 2))];
        pts.push_back(repeated);
        accepted = Sound2(pts);
    }
    if (!accepted)
    {
        Fallback2(static_cast<size_t>(n - 1), pts);
        pts.push_back(pts[0]);
    }
    for (size_t i = 0; i < pts.size(); ++i)
    {
        io.givenVec<2>(pts[i]);
    }

    std::vector<size_t> duplicates{};
    {
        Delaunay2<double> probe{};
        probe(pts);
        duplicates = probe.GetDuplicates();
    }
    std::array<int32_t, 2> edge{ 0, 0 };
    if (io.rawInteger(0, 1) == 0)
    {
        int32_t v = io.rawInteger(0, n - 1);
        edge = { v, v };
    }
    else
    {
        // The last point duplicates an earlier one.
        edge = { static_cast<int32_t>(duplicates[static_cast<size_t>(n - 1)]), n - 1 };
    }
    GivenInt(io, edge[0]);
    GivenInt(io, edge[1]);

    CDT2 cdt{};
    bool built = cdt(pts);
    std::vector<int32_t> partitionedEdge{};
    cdt.Insert(edge, partitionedEdge);
    io.outBool(built);
    EmitPartition(io, partitionedEdge);
}

// DELIBERATE DEVIATION (#325, ConstrainedDelaunay2 finding 1). operator()
// forwards to Delaunay2<T>::operator(), which never touches mInsertedEdges,
// so a functor reused for a second data set reports the first set's
// constraints from GetInsertedEdges(); those keys index unrelated vertices
// of the new triangulation. The port clears the set in compute(). The
// record: triangulate A, insert one to three constraints, triangulate B
// (with the other operator() overload), emit GetInsertedEdges(), insert one
// constraint into B and emit the partition and GetInsertedEdges() again.
// Every record whose first set received a constraint not also inserted into
// B deviates; the partitions agree.
ORACLE_CASE("ConstrainedDelaunay2.deviation.staleInsertedEdges")
{
    int32_t nA = io.integer(3, 9);
    auto ptsA = SoundPoints2(io, 1, static_cast<size_t>(nA));
    int32_t nB = io.integer(3, 9);
    auto ptsB = SoundPoints2(io, 0, static_cast<size_t>(nB));

    std::vector<size_t> dupA{}, dupB{};
    std::vector<std::array<int32_t, 2>> edgesA0{}, edgesB0{};
    {
        Delaunay2<double> probe{};
        probe(ptsA);
        dupA = probe.GetDuplicates();
        edgesA0 = SortedGraphEdges(probe);
        probe(ptsB);
        dupB = probe.GetDuplicates();
        edgesB0 = SortedGraphEdges(probe);
    }
    int32_t numEdgesA = io.integer(1, 3);
    std::vector<std::array<int32_t, 2>> edgesA(static_cast<size_t>(numEdgesA));
    for (size_t k = 0; k < edgesA.size(); ++k)
    {
        edgesA[k] = RawConstraint(io, ptsA, dupA, edgesA0, -1);
        GivenInt(io, edgesA[k][0]);
        GivenInt(io, edgesA[k][1]);
    }
    std::array<int32_t, 2> edgeB = RawConstraint(io, ptsB, dupB, edgesB0, -1);
    GivenInt(io, edgeB[0]);
    GivenInt(io, edgeB[1]);

    CDT2 cdt{};
    cdt(ptsA);
    for (auto const& edge : edgesA)
    {
        std::vector<int32_t> partitionedEdge{};
        cdt.Insert(edge, partitionedEdge);
    }
    bool built = cdt(ptsB.size(), ptsB.data());
    io.outBool(built);
    std::vector<int32_t> stale{};
    for (auto const& ekey : cdt.GetInsertedEdges())
    {
        stale.push_back(ekey.V[0]);
        stale.push_back(ekey.V[1]);
    }
    EmitSortedTuples<2, int32_t>(io, stale);
    std::vector<int32_t> partitionedEdge{};
    cdt.Insert(edgeB, partitionedEdge);
    EmitPartition(io, partitionedEdge);
    EmitCDT2(io, cdt);
}

// DELIBERATE DEVIATION (#325, ConstrainedDelaunay2 finding 2). Insert
// substitutes 'edge[i] = duplicates[edge[i]]' before its range assert, so an
// out-of-range index is read from beyond the end of mDuplicates. The read is
// made deterministic here: the functor first triangulates a set A of nA
// points whose point j (nB <= j < nA) repeats point r < nB, so A leaves
// mDuplicates[j] = r; it then triangulates the prefix B of nB points, and
// Delaunay2<T>::operator() does 'mDuplicates.clear(); resize(nB)', which
// keeps the capacity and the stale element j. Insert({j, s}) on B then reads
// r, passes the assert, and upstream inserts the unrelated constraint
// <r, s> into B and returns normally. The port validates the raw indices
// first and throws "Invalid edge." Every accepted record deviates (C++
// returns, the port throws).
ORACLE_CASE("ConstrainedDelaunay2.deviation.duplicatesRead")
{
    int32_t nB = io.integer(3, 7);
    int32_t m = io.integer(1, 3);
    std::vector<Vector2<double>> ptsA{};
    int32_t j = -1, s = -1;
    for (int32_t attempt = 0; attempt < 64 && j < 0; ++attempt)
    {
        std::vector<Vector2<double>> ptsB{};
        RawPoints2(io, 1, static_cast<size_t>(nB), ptsB);
        ptsA = ptsB;
        for (int32_t k = 0; k < m; ++k)
        {
            Vector2<double> extra = ptsB[static_cast<size_t>(io.rawInteger(0, nB - 1))];
            if (k > 0 && io.rawInteger(0, 1) != 0)
            {
                extra[0] = static_cast<double>(io.rawInteger(-3, 3));
                extra[1] = static_cast<double>(io.rawInteger(-3, 3));
            }
            ptsA.push_back(extra);
        }
        if (!Sound2(ptsA) || !Sound2(ptsB))
        {
            continue;
        }
        Delaunay2<double> probe{};
        probe(ptsA);
        std::vector<size_t> dupA = probe.GetDuplicates();
        probe(ptsB);
        std::vector<size_t> dupB = probe.GetDuplicates();
        for (int32_t k = nB; k < nB + m && j < 0; ++k)
        {
            size_t r = dupA[static_cast<size_t>(k)];
            if (r < static_cast<size_t>(nB) && dupB[r] == r)
            {
                for (int32_t t = 0; t < nB; ++t)
                {
                    if (dupB[static_cast<size_t>(t)] != r)
                    {
                        j = k;
                        s = t;
                        break;
                    }
                }
            }
        }
    }
    if (j < 0)
    {
        // Points on a parabola (general position) with every extra point
        // repeating point 0, the first extreme processed.
        Fallback2(static_cast<size_t>(nB), ptsA);
        Vector2<double> first = ptsA[0];
        ptsA.resize(static_cast<size_t>(nB + m), first);
        j = nB;
        s = 1;
    }
    for (size_t i = 0; i < ptsA.size(); ++i)
    {
        io.givenVec<2>(ptsA[i]);
    }
    bool swapped = (io.rawInteger(0, 1) != 0);
    std::array<int32_t, 2> edge{ swapped ? s : j, swapped ? j : s };
    GivenInt(io, edge[0]);
    GivenInt(io, edge[1]);

    std::vector<Vector2<double>> ptsB(ptsA.begin(), ptsA.begin() + nB);
    CDT2 cdt{};
    cdt(ptsA);
    bool built = cdt(ptsB);
    io.outBool(built);
    std::vector<int32_t> partitionedEdge{};
    cdt.Insert(edge, partitionedEdge);
    EmitPartition(io, partitionedEdge);
    EmitCDT2(io, cdt);
}
