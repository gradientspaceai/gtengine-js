// Verify group 12 (computational geometry): MinimumVolumeBox3FloatingPoint.h
// and MinimumVolumeBox3Rational.h, the two implementation headers behind the
// MinimumVolumeBox3.h wrapper.
//
// Group 8 (oracle/cpp/cases/v08-compgeom.cpp) already runs both headers
// through the wrapper's four operator() overloads. It compares a fixed mesh
// exactly when upstream's answer is stable under reordering the triangles,
// and a point cloud only by dimension, flags and a toleranced volume, because
// the winning candidate depends on a std::unordered_map iteration order
// (finding #530). This file closes that gap.
//
// WHAT DECIDES THE ANSWER. Neither header calls the C math library: the only
// function is std::sqrt, on doubles and on BSRational (which GTE converts to
// double first). Everything else is + - * / on doubles or exact BSNumber /
// BSRational arithmetic. So there is nothing libm-decided; the tolerance of
// v08 is entirely order-decided, by three unspecified orders:
//   (1) ExtractMeshTopology numbers the edges by iterating ETManifoldMesh's
//       unordered_map mEMap (and the triangles by mTMap); mEdgeIndices and
//       the first-minimum tie break follow that numbering, and
//       RemoveCoplanarTriangleAdjacencies visits the edges in it;
//   (2) ExtractVertexAdjacencies copies each vertex's VAdjacent, an
//       std::unordered_set<int32_t>, into the adjacency pool, and the pool
//       order decides GetExtreme's climb path and which of several tied
//       vertices it returns;
//   (3) for a point cloud, ConvexHull3::GetHull lists the hull triangles in
//       mTMap order, and that order is the mesh insertion order, which
//       decides each edge's V[0]/V[1] orientation (ComputeVolume's assumed
//       minimum support vertex) and T[0]/T[1] order (the N/M roles).
// The port enumerates all three in sorted-key order (std::map order).
//
// THE CONTROL. The subclass Canonical<Base> below runs upstream's own
// protected pipeline step by step, exactly as operator() does, with two of
// its steps replaced by copies that differ only in iterating sorted copies of
// the hash containers (ExtractMeshTopologySorted, ExtractVertexAdjacencies-
// Sorted), and for a point cloud it sorts the ConvexHull3 triangles (which
// are TriangleKey tuples) before building the mesh. Everything else - the
// convex hull, the geometry, the coplanar merge, the aligned candidate, all
// 81 level-curve processors, the four minimizers, ComputeVolume, GetExtreme,
// the threaded candidate search and the rational box - is upstream's code.
// The order is then the port's, and the full box compares bit for bit. Every
// record also runs upstream's raw operator() (all four overloads are used)
// and records, as a diagnostic input the replay ignores, whether the raw
// answer equals the canonical one: the fraction of records on which the
// hash order changes the box is thereby measured on the MSVC build.
//
// DELIBERATE PORT FIXES (docs/UPSTREAM-FINDINGS.md) and how they are kept
// out of the exact cases:
//   #352/#355  dimension-2 Newell loop: the exact dimension-2 case accepts
//              only point sets on which upstream's orthonormal basis equals
//              the corrected one bit for bit (probe below);
//   #405/#426  ComputeVolume's support assumption and GetExtreme's plateau
//              stall (floating point): the main floating-point cases reject
//              records whose canonical upstream box does not contain the
//              points; GetExtreme is compared against an exact maximizer;
//   #355       MinimizerVariableT rounds its parameters to double (rational):
//              the rational cases reject records that reach it, except the
//              minimizer case, whose subclass replaces MinimizerVariableT and
//              so compares its (rounded) arguments exactly.
// Design deviation inherited from MinimumAreaBox2 (port note in
// src/MinimumAreaBox2.ts, v11 MinimumAreaBox2.deviation.floatComputeType):
// the dimension-2 path of both headers calls MinimumAreaBox2<T, T>, the
// floating-point compute type, and the port's MinimumAreaBox2 is exact-only.
// The exact dimension-2 case accepts only projections on which upstream's
// floating-point rectangle equals the exact one bit for bit; the rest is a
// deviation case.
#define ORACLE_FAMILY "v12-compgeom"
#include "Oracle.h"

#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/ConvexHull3.h>
#include <Mathematics/MinimumAreaBox2.h>
#include <Mathematics/MinimumVolumeBox3FloatingPoint.h>
#include <Mathematics/MinimumVolumeBox3Rational.h>
#include <Mathematics/UniqueVerticesSimplices.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>

using namespace gte;

namespace
{
    // Both specializations put exact-arithmetic objects of several KiB on
    // the stack (the rational header asks for a 1 GiB stack reserve), so the
    // cases run on one reused worker thread with a large stack reservation,
    // as in v08. The callable's exception, if any, is rethrown on the calling
    // thread so the harness records the throw.
    struct BigStackCall
    {
        std::function<void()> fn;
        std::exception_ptr error;
    };

    BigStackCall gBigStackCall{};
    HANDLE gBigStackThread = nullptr;
    HANDLE gBigStackRequest = nullptr;
    HANDLE gBigStackDone = nullptr;

    DWORD WINAPI BigStackWorker(LPVOID)
    {
        for (;;)
        {
            WaitForSingleObject(gBigStackRequest, INFINITE);
            try
            {
                gBigStackCall.fn();
            }
            catch (...)
            {
                gBigStackCall.error = std::current_exception();
            }
            SetEvent(gBigStackDone);
        }
    }

    void RunWithBigStack(std::function<void()> fn)
    {
        if (gBigStackThread == nullptr)
        {
            gBigStackRequest = CreateEventW(nullptr, FALSE, FALSE, nullptr);
            gBigStackDone = CreateEventW(nullptr, FALSE, FALSE, nullptr);
            gBigStackThread = CreateThread(nullptr, 512ull * 1024ull * 1024ull,
                BigStackWorker, nullptr, STACK_SIZE_PARAM_IS_A_RESERVATION, nullptr);
            if (gBigStackThread == nullptr)
            {
                throw std::runtime_error("CreateThread failed");
            }
        }

        gBigStackCall.fn = std::move(fn);
        gBigStackCall.error = nullptr;
        SetEvent(gBigStackRequest);
        WaitForSingleObject(gBigStackDone, INFINITE);
        if (gBigStackCall.error)
        {
            std::exception_ptr error = gBigStackCall.error;
            gBigStackCall.error = nullptr;
            std::rethrow_exception(error);
        }
    }

    using MVB3FP = MinimumVolumeBox3<double, int32_t, MVB3FloatingPoint>;
    using MVB3R = MinimumVolumeBox3<double, int32_t, MVB3Rational>;
    using Vec3 = Vector3<double>;

    std::size_t constexpr invalidIndex = std::numeric_limits<std::size_t>::max();

    // The port stores invalid indices as -1.
    int32_t IndexOut(std::size_t i)
    {
        return (i == invalidIndex ? -1 : static_cast<int32_t>(i));
    }
}

namespace
{
    // Upstream's pipeline with the three unspecified orders pinned to sorted
    // key order. See the file comment. Base is MVB3FP or MVB3R; both define
    // the same protected members and steps.
    template <typename Base>
    class Canonical : public Base
    {
    public:
        Canonical(std::size_t numThreads) : Base(numThreads) {}

        // Upstream's operator()(numVertices, vertices, numIndices, indices,
        // lgMaxSample, box, volume) with ExtractMeshTopology and
        // ExtractVertexAdjacencies replaced by their sorted-order copies.
        void RunHull(std::size_t numVertices, Vec3 const* vertices,
            std::size_t numIndices, int32_t const* indices, std::size_t lgMaxSample,
            OrientedBox3<double>& box, double& volume)
        {
            LogAssert(
                numVertices >= 4 && vertices != nullptr &&
                numIndices >= 12 && (numIndices % 3) == 0 && indices != nullptr &&
                lgMaxSample >= 2,
                "Invalid argument.");
            this->GenerateSubdivision(lgMaxSample);
            std::size_t const numTriangles = numIndices / 3;
            VETManifoldMesh mesh{};
            this->CreateMeshTopology(numTriangles, indices, mesh);
            ExtractMeshTopologySorted(mesh);
            ExtractVertexAdjacenciesSorted(mesh);
            this->ExtractMeshGeometry(numVertices, vertices);
            this->RemoveCoplanarTriangleAdjacencies();
            this->ComputeAlignedCandidate();
            this->GetMinimumVolumeCandidate();
            this->GetMinimumVolumeBox(box, volume);
        }

        void RunHull(std::vector<Vec3> const& vertices, std::vector<int32_t> const& indices,
            std::size_t lgMaxSample, OrientedBox3<double>& box, double& volume)
        {
            RunHull(vertices.size(), vertices.data(), indices.size(), indices.data(),
                lgMaxSample, box, volume);
        }

        // Upstream's operator()(numPoints, points, lgMaxSample, box, volume)
        // with the ConvexHull3 triangles sorted before the mesh is built.
        // The triangles are TriangleKey<true> tuples (smallest index first)
        // and UniqueVerticesSimplices renumbers the used vertices in
        // increasing order, so a lexicographic sort of the triples is the
        // port's std::map order.
        std::size_t RunCloud(std::vector<Vec3> const& points, std::size_t lgMaxSample,
            OrientedBox3<double>& box, double& volume)
        {
            LogAssert(points.size() >= 4 && lgMaxSample >= 2, "Invalid argument.");
            std::vector<Vec3> hullVertices{};
            std::vector<int32_t> hullIndices{};
            std::size_t dimension = this->ComputeConvexHull(0, points.size(), points.data(),
                hullVertices, hullIndices, box, volume);
            if (dimension == 3)
            {
                SortTriangles(hullIndices);
                RunHull(hullVertices, hullIndices, lgMaxSample, box, volume);
            }
            return dimension;
        }

        static void SortTriangles(std::vector<int32_t>& indices)
        {
            std::vector<std::array<int32_t, 3>> triangles(indices.size() / 3);
            for (std::size_t t = 0; t < triangles.size(); ++t)
            {
                triangles[t] = { indices[3 * t], indices[3 * t + 1], indices[3 * t + 2] };
            }
            std::sort(triangles.begin(), triangles.end());
            for (std::size_t t = 0; t < triangles.size(); ++t)
            {
                for (std::size_t j = 0; j < 3; ++j)
                {
                    indices[3 * t + j] = triangles[t][j];
                }
            }
        }

        // ---- read access for the replay ----
        auto const& Winner() const { return this->mMinimumVolumeObject; }
        auto const& Aligned() const { return this->mAlignedCandidate; }
        auto const& Edges() const { return this->mEdges; }
        std::size_t NumEdgePairs() const { return this->mEdgeIndices.size(); }
        std::size_t ClimbStart() const { return this->mVClimbStart; }
        std::size_t NumVertices() const { return this->mAdjacentPoolLocation.size(); }
        std::size_t MaxSample() const { return this->mMaxSample; }
        auto const& DomainIndex() const { return this->mDomainIndex; }
        // Floating-point header only (instantiated only there).
        auto const& TVertices() const { return this->mTVertices; }
        auto const& NVertices() const { return this->mNVertices; }

        std::vector<std::size_t> Adjacent(std::size_t v) const
        {
            std::size_t const* adjacent = &this->mAdjacentPool[this->mAdjacentPoolLocation[v]];
            return std::vector<std::size_t>(adjacent + 1, adjacent + 1 + adjacent[0]);
        }

    private:
        // ExtractMeshTopology verbatim, except that eMap and tMap are sorted
        // copies of the mesh's hash maps.
        void ExtractMeshTopologySorted(VETManifoldMesh const& mesh)
        {
            std::map<EdgeKey<false>, ETManifoldMesh::Edge*> eMap{};
            for (auto const& element : mesh.GetEdges())
            {
                eMap.emplace(element.first, element.second.get());
            }
            std::map<TriangleKey<true>, ETManifoldMesh::Triangle*> tMap{};
            for (auto const& element : mesh.GetTriangles())
            {
                tMap.emplace(element.first, element.second.get());
            }

            std::map<ETManifoldMesh::Edge*, std::size_t> edgeIndexMap{};
            std::size_t index = 0;
            for (auto const& element : eMap)
            {
                edgeIndexMap.emplace(element.second, index);
                for (std::size_t j = 0; j < 2; ++j)
                {
                    this->mEdges[index].v[j] = static_cast<std::size_t>(element.second->V[j]);
                }
                ++index;
            }

            std::map<ETManifoldMesh::Triangle*, std::size_t> triangleIndexMap{};
            index = 0;
            for (auto const& element : tMap)
            {
                triangleIndexMap.emplace(element.second, index);
                for (std::size_t j = 0; j < 3; ++j)
                {
                    this->mTriangles[index].v[j] = static_cast<std::size_t>(element.second->V[j]);
                }
                ++index;
            }

            index = 0;
            for (auto const& element : eMap)
            {
                for (std::size_t j = 0; j < 2; ++j)
                {
                    auto tIter = triangleIndexMap.find(element.second->T[j]);
                    this->mEdges[index].t[j] = tIter->second;
                }
                ++index;
            }

            index = 0;
            for (auto const& element : tMap)
            {
                for (std::size_t j = 0; j < 3; ++j)
                {
                    auto eIter = edgeIndexMap.find(element.second->E[j]);
                    this->mTriangles[index].e[j] = eIter->second;
                }
                for (std::size_t j = 0; j < 3; ++j)
                {
                    auto tIter = triangleIndexMap.find(element.second->T[j]);
                    this->mTriangles[index].t[j] = tIter->second;
                }
                ++index;
            }

            for (std::size_t e0 = 0; e0 < this->mEdges.size(); ++e0)
            {
                for (std::size_t e1 = e0 + 1; e1 < this->mEdges.size(); ++e1)
                {
                    this->mEdgeIndices.push_back({ e0, e1 });
                }
            }
        }

        // ExtractVertexAdjacencies verbatim, except that each VAdjacent
        // unordered_set is visited in increasing order.
        void ExtractVertexAdjacenciesSorted(VETManifoldMesh const& mesh)
        {
            std::map<std::int32_t, VETManifoldMesh::Vertex*> sortedVMap{};
            for (auto const& element : mesh.GetVertices())
            {
                sortedVMap.emplace(element.first, element.second.get());
            }

            std::size_t numAdjacentPool = 0;
            for (auto const& element : sortedVMap)
            {
                numAdjacentPool += element.second->VAdjacent.size() + 1;
            }
            this->mAdjacentPool.resize(numAdjacentPool);
            this->mAdjacentPoolLocation.resize(sortedVMap.size());
            std::size_t apIndex = 0, vaIndex = 0;
            for (auto const& element : sortedVMap)
            {
                std::set<std::int32_t> adjacent(element.second->VAdjacent.begin(),
                    element.second->VAdjacent.end());
                this->mAdjacentPoolLocation[vaIndex++] = apIndex;
                this->mAdjacentPool[apIndex++] = adjacent.size();
                for (auto v : adjacent)
                {
                    this->mAdjacentPool[apIndex++] = static_cast<std::size_t>(v);
                }
            }
        }
    };

    using CanonicalFP = Canonical<MVB3FP>;
    using CanonicalR = Canonical<MVB3R>;
}

namespace
{
    // ---- generators ----------------------------------------------------
    //
    // Nothing here is recorded; the callers record the final coordinates.
    // Cloud modes:
    //   0  small lattice [-3,3]^3: exact ties, coplanar hull faces and
    //      collinear hull edges, so RemoveCoplanarTriangleAdjacencies takes
    //      both of its branches and many candidates tie exactly;
    //   1  uniform [-5,5]^3: the generic branch;
    //   2  a lattice cloud rotated by 0.03 to 0.3 degrees (sin/cos applied
    //      here, only the rotated coordinates are recorded): the regime of
    //      #426, where coplanarity is approximate;
    //   3  a random subset of a {0,1,2} x {0..b} x {0..c} grid with the eight
    //      corners forced: many hull vertices in face interiors and on edges;
    //   4  a lattice cloud rotated by the exact Pythagorean rotation
    //      (3/5, 4/5) about a coordinate axis and scaled by 5, so that the
    //      coordinates stay integers and the minimum box is tilted: exact
    //      ties between different edge pairs.
    std::vector<Vec3> DrawCloud(oracle::Ctx& io, int32_t mode, int32_t n)
    {
        std::vector<Vec3> points{};
        if (mode == 3)
        {
            double b = static_cast<double>(io.rawInteger(1, 2));
            double c = static_cast<double>(io.rawInteger(1, 2));
            for (int32_t k = 0; k < 8; ++k)
            {
                points.push_back(Vec3{ (k & 1) ? 2.0 : 0.0, (k & 2) ? b : 0.0,
                    (k & 4) ? c : 0.0 });
            }
            for (int32_t i = 8; i < n + 4; ++i)
            {
                double x = static_cast<double>(io.rawInteger(0, 2));
                double y = static_cast<double>(io.rawInteger(0, static_cast<int32_t>(b)));
                double z = static_cast<double>(io.rawInteger(0, static_cast<int32_t>(c)));
                points.push_back(Vec3{ x, y, z });
            }
            // Shuffle so that the corners are not always first.
            for (std::size_t i = points.size(); i > 1; --i)
            {
                std::size_t j = static_cast<std::size_t>(io.rawInteger(0,
                    static_cast<int32_t>(i) - 1));
                std::swap(points[i - 1], points[j]);
            }
            return points;
        }

        for (int32_t i = 0; i < n; ++i)
        {
            Vec3 p{};
            for (int32_t d = 0; d < 3; ++d)
            {
                p[d] = (mode == 1 ? io.raw(-5.0, 5.0) : static_cast<double>(io.rawInteger(-3, 3)));
            }
            points.push_back(p);
        }

        if (mode == 2 || mode == 4)
        {
            // Mode 4 uses the integer matrix 5 * R with R the rotation by
            // (cos, sin) = (3/5, 4/5), so every product is an exact integer.
            double c = 3.0, s = 4.0, scale = 5.0;
            if (mode == 2)
            {
                double angle = io.raw(0.0005, 0.005);
                c = std::cos(angle);
                s = std::sin(angle);
                scale = 1.0;
            }
            int32_t axis = io.rawInteger(0, 2);
            int32_t a0 = (axis + 1) % 3, a1 = (axis + 2) % 3;
            for (auto& p : points)
            {
                double u = p[a0], v = p[a1];
                p[a0] = c * u - s * v;
                p[a1] = s * u + c * v;
                p[axis] = scale * p[axis];
            }
        }
        return points;
    }

    void GiveCloud(oracle::Ctx& io, std::vector<Vec3> const& points)
    {
        int32_t n = static_cast<int32_t>(points.size());
        io.integer(n, n);
        for (auto const& p : points)
        {
            io.givenVec(p);
        }
    }

    void GiveMesh(oracle::Ctx& io, std::vector<Vec3> const& vertices,
        std::vector<int32_t> const& indices)
    {
        GiveCloud(io, vertices);
        int32_t numIndices = static_cast<int32_t>(indices.size());
        io.integer(numIndices, numIndices);
        for (int32_t i : indices)
        {
            io.integer(i, i);
        }
    }

    // The convex hull mesh of a cloud, built with upstream's ConvexHull3 and
    // UniqueVerticesSimplices exactly as MinimumVolumeBox3's point-cloud
    // query builds it. Returns false unless the hull is 3-dimensional. The
    // triangle order is then scrambled and each triple rotated (orientation
    // preserved), because the vertices-and-indices query takes the caller's
    // order as given; the replay passes the same order.
    bool MakeHullMesh(oracle::Ctx& io, std::vector<Vec3> const& points,
        std::vector<Vec3>& vertices, std::vector<int32_t>& indices)
    {
        ConvexHull3<double> ch3{};
        ch3(points.size(), points.data(), 0);
        if (ch3.GetDimension() != 3)
        {
            return false;
        }
        std::vector<int32_t> source{};
        for (auto index : ch3.GetHull())
        {
            source.push_back(static_cast<int32_t>(index));
        }
        UniqueVerticesSimplices<Vec3, int32_t, 3> uvs{};
        uvs.RemoveDuplicateAndUnusedVertices(points, source, vertices, indices);

        std::size_t const numTriangles = indices.size() / 3;
        std::vector<std::array<int32_t, 3>> triangles(numTriangles);
        for (std::size_t t = 0; t < numTriangles; ++t)
        {
            int32_t r = io.rawInteger(0, 2);
            for (int32_t j = 0; j < 3; ++j)
            {
                triangles[t][j] = indices[3 * t + static_cast<std::size_t>((j + r) % 3)];
            }
        }
        for (std::size_t i = numTriangles; i > 1; --i)
        {
            std::size_t j = static_cast<std::size_t>(io.rawInteger(0,
                static_cast<int32_t>(i) - 1));
            std::swap(triangles[i - 1], triangles[j]);
        }
        for (std::size_t t = 0; t < numTriangles; ++t)
        {
            for (std::size_t j = 0; j < 3; ++j)
            {
                indices[3 * t + j] = triangles[t][j];
            }
        }
        return true;
    }

    // Maximum amount by which a point sticks out of the box, measured along
    // the box axes. Zero when the box contains every point.
    double ContainmentViolation(OrientedBox3<double> const& box,
        std::vector<Vec3> const& points)
    {
        double worst = 0.0;
        for (auto const& p : points)
        {
            Vec3 d = p - box.center;
            for (int32_t i = 0; i < 3; ++i)
            {
                double t = std::fabs(Dot(d, box.axis[i])) - box.extent[i];
                worst = std::max(worst, t);
            }
        }
        return worst;
    }

    double PointScale(std::vector<Vec3> const& points)
    {
        double scale = 1.0;
        for (auto const& p : points)
        {
            for (int32_t i = 0; i < 3; ++i)
            {
                scale = std::max(scale, std::fabs(p[i]));
            }
        }
        return scale;
    }

    bool Contains(OrientedBox3<double> const& box, std::vector<Vec3> const& points)
    {
        return ContainmentViolation(box, points) <= 1e-9 * PointScale(points);
    }

    bool SameBits(double a, double b)
    {
        return std::memcmp(&a, &b, sizeof(double)) == 0;
    }

    bool SameBox(OrientedBox3<double> const& b0, double v0,
        OrientedBox3<double> const& b1, double v1)
    {
        bool same = SameBits(v0, v1);
        for (int32_t i = 0; i < 3; ++i)
        {
            same = same && SameBits(b0.center[i], b1.center[i])
                && SameBits(b0.extent[i], b1.extent[i]);
            for (int32_t j = 0; j < 3; ++j)
            {
                same = same && SameBits(b0.axis[i][j], b1.axis[i][j]);
            }
        }
        return same;
    }

    // v08's invariants of a box: the centre, the sorted extents and the six
    // entries of sum_i extent[i]^2 * axis[i] * axis[i]^T; together they
    // determine the box up to the signs and the order of its axes.
    std::array<double, 12> BoxInvariants(OrientedBox3<double> const& box)
    {
        std::array<double, 12> out{};
        std::array<double, 3> extent{ box.extent[0], box.extent[1], box.extent[2] };
        std::sort(extent.begin(), extent.end());
        for (int32_t i = 0; i < 3; ++i)
        {
            out[i] = box.center[i];
            out[3 + i] = extent[i];
            double w = box.extent[i] * box.extent[i];
            Vec3 const& a = box.axis[i];
            out[6] += w * a[0] * a[0];
            out[7] += w * a[0] * a[1];
            out[8] += w * a[0] * a[2];
            out[9] += w * a[1] * a[1];
            out[10] += w * a[1] * a[2];
            out[11] += w * a[2] * a[2];
        }
        return out;
    }

    // The diagnostic recorded by the canonical cases: how upstream's raw
    // operator() (hash order) compares with the canonical run. 3: the same
    // box bit for bit; 2: the same box up to axis signs and order (the same
    // invariants bit for bit); 1: the same volume, a different box; 0: a
    // different volume.
    int32_t RawAgreement(OrientedBox3<double> const& box, double volume,
        OrientedBox3<double> const& rawBox, double rawVolume)
    {
        if (SameBox(box, volume, rawBox, rawVolume))
        {
            return 3;
        }
        if (!SameBits(volume, rawVolume))
        {
            return 0;
        }
        auto i0 = BoxInvariants(box), i1 = BoxInvariants(rawBox);
        for (std::size_t i = 0; i < i0.size(); ++i)
        {
            if (!SameBits(i0[i], i1[i]))
            {
                return 1;
            }
        }
        return 2;
    }

    void OutBox(oracle::Ctx& io, OrientedBox3<double> const& box, double volume)
    {
        io.outVec(box.center);
        for (int32_t i = 0; i < 3; ++i)
        {
            io.outVec(box.axis[i]);
        }
        for (int32_t i = 0; i < 3; ++i)
        {
            io.outReal(box.extent[i]);
        }
        io.outReal(volume);
    }
}

namespace
{
    // The fields of a Candidate that upstream sets on the path that produced
    // it. For the aligned candidate that is the support indices and the
    // volume (ComputeAlignedCandidate); for the winner, everything.
    void OutAlignedFP(oracle::Ctx& io, CanonicalFP const& query)
    {
        auto const& c = query.Aligned();
        for (int32_t i = 0; i < 3; ++i)
        {
            io.outInt(IndexOut(c.minSupportIndex[i]));
            io.outInt(IndexOut(c.maxSupportIndex[i]));
        }
        io.outReal(c.volume);
    }

    void OutWinnerFP(oracle::Ctx& io, CanonicalFP const& query)
    {
        auto const& c = query.Winner();
        io.outInt(IndexOut(c.edgeIndex[0]));
        io.outInt(IndexOut(c.edgeIndex[1]));
        io.outInt(IndexOut(c.levelCurveProcessorIndex));
        for (int32_t i = 0; i < 2; ++i)
        {
            io.outVec(c.N[i]);
            io.outVec(c.M[i]);
        }
        io.outReal(c.f00);
        io.outReal(c.f10);
        io.outReal(c.f01);
        io.outReal(c.f11);
        for (int32_t i = 0; i < 3; ++i)
        {
            io.outVec(c.axis[i]);
            io.outInt(IndexOut(c.minSupportIndex[i]));
            io.outInt(IndexOut(c.maxSupportIndex[i]));
        }
        io.outReal(c.volume);
    }

    // The mesh topology after RemoveCoplanarTriangleAdjacencies: the edges
    // (vertex and triangle indices), the number of edge pairs, the climb
    // start and every adjacency list, plus the subdivision of the sample
    // domain. Shared by both headers.
    template <typename Query>
    void OutTopology(oracle::Ctx& io, Query const& query)
    {
        auto const& edges = query.Edges();
        io.outInt(static_cast<int32_t>(edges.size()));
        for (auto const& e : edges)
        {
            io.outInt(IndexOut(e.v[0]));
            io.outInt(IndexOut(e.v[1]));
            io.outInt(IndexOut(e.t[0]));
            io.outInt(IndexOut(e.t[1]));
        }
        io.outInt(static_cast<int32_t>(query.NumEdgePairs()));
        io.outInt(IndexOut(query.ClimbStart()));
        for (std::size_t v = 0; v < query.NumVertices(); ++v)
        {
            std::vector<std::size_t> adjacent = query.Adjacent(v);
            io.outInt(static_cast<int32_t>(adjacent.size()));
            for (auto a : adjacent)
            {
                io.outInt(IndexOut(a));
            }
        }
        io.outInt(static_cast<int32_t>(query.MaxSample()));
        for (auto const& item : query.DomainIndex())
        {
            io.outInt(static_cast<int32_t>(item[0]));
            io.outInt(static_cast<int32_t>(item[1]));
            io.outInt(static_cast<int32_t>(item[2]));
        }
    }

    // THE EXACT SEPARATOR for the floating-point header's two deliberate
    // fixes. The port replaces upstream's support vertex along a candidate
    // axis (the hull-edge vertex assumed by ComputeVolume, #405, or the end of
    // GetExtreme's hill climb, #426) only when another vertex projects beyond
    // it by more than a rigorous bound on the rounding error of the two dot
    // products, i.e. only when upstream's vertex is provably not the extreme
    // one in exact arithmetic over the translated vertices mTVertices that
    // the climb uses. This returns true when all six support vertices of
    // the winning candidate are exact extremes along its three axes, so the
    // fixes cannot fire on the winner; then the rational box upstream builds
    // from them contains the hull exactly. A non-winning candidate on which a
    // fix fires only gets a larger volume in the port, so it cannot win there
    // either (barring an exact volume tie, which the deep run would show).
    using ExactNumber = BSNumber<UIntegerAP32>;

    ExactNumber ExactDot(Vec3 const& u, Vec3 const& v)
    {
        ExactNumber sum = ExactNumber(u[0]) * ExactNumber(v[0]);
        sum = sum + ExactNumber(u[1]) * ExactNumber(v[1]);
        return sum + ExactNumber(u[2]) * ExactNumber(v[2]);
    }

    bool SupportsAreExact(CanonicalFP const& query)
    {
        auto const& c = query.Winner();
        auto const& T = query.TVertices();
        for (int32_t i = 0; i < 3; ++i)
        {
            ExactNumber dmin = ExactDot(c.axis[i], T[c.minSupportIndex[i]]);
            ExactNumber dmax = ExactDot(c.axis[i], T[c.maxSupportIndex[i]]);
            for (auto const& v : T)
            {
                ExactNumber d = ExactDot(c.axis[i], v);
                if (d < dmin || d > dmax)
                {
                    return false;
                }
            }
        }
        return true;
    }

    std::vector<Vec3> FallbackCloud()
    {
        return { Vec3{ 0.0, 0.0, 0.0 }, Vec3{ 4.0, 0.0, 0.0 }, Vec3{ 0.0, 3.0, 0.0 },
            Vec3{ 0.0, 0.0, 2.0 }, Vec3{ 1.0, 1.0, 1.0 } };
    }
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.compute.canonical")
{
    RunWithBigStack([&io]() {
    // The point-cloud query (overloads 1 and 2) with the unspecified orders
    // pinned (see the file comment), compared bit for bit: this replaces
    // v08's volume-only comparison with tolerance 2e-2. Restricted to clouds
    // whose hull is 3-dimensional and on which the six support vertices of
    // upstream's winning candidate are exact extremes (SupportsAreExact, the
    // exact separator of #405 and #426; the loop is capped at 24 attempts,
    // redraws the whole cloud and falls back to the cloud whose box misses
    // the points by the least). numThreads > 0 runs upstream's std::thread
    // candidate search; the port is single-threaded and must agree.
    int32_t mode = io.integer(0, 4);
    int32_t lgMaxSample = io.integer(2, 4);
    int32_t n = io.integer(5, 10);
    int32_t numThreads = io.integer(0, 4);

    std::vector<Vec3> best{};
    double bestViolation = std::numeric_limits<double>::max();
    for (int32_t attempt = 0; attempt < 24; ++attempt)
    {
        std::vector<Vec3> points = DrawCloud(io, mode, n);
        OrientedBox3<double> box{};
        double volume = 0.0;
        CanonicalFP probe(0);
        std::size_t dimension = 0;
        try
        {
            dimension = probe.RunCloud(points, static_cast<std::size_t>(lgMaxSample),
                box, volume);
        }
        catch (std::exception const&)
        {
            continue;
        }
        if (dimension != 3)
        {
            continue;
        }
        double violation = ContainmentViolation(box, points) / PointScale(points);
        bool sound = SupportsAreExact(probe);
        if (sound || violation < bestViolation)
        {
            bestViolation = violation;
            best = points;
        }
        if (sound)
        {
            break;
        }
    }
    if (best.empty())
    {
        best = FallbackCloud();
    }
    GiveCloud(io, best);

    // Upstream's raw query, overload 1 or 2. Diagnostic only: whether the
    // hash order changed the answer (RawAgreement) is recorded as an input the replay
    // ignores.
    OrientedBox3<double> rawBox{};
    double rawVolume = 0.0;
    MVB3FP raw(static_cast<std::size_t>(numThreads));
    if (io.index() % 2 == 0)
    {
        raw(best.size(), best.data(), static_cast<std::size_t>(lgMaxSample), rawBox, rawVolume);
    }
    else
    {
        raw(best, static_cast<std::size_t>(lgMaxSample), rawBox, rawVolume);
    }

    OrientedBox3<double> box{};
    double volume = 0.0;
    CanonicalFP query(static_cast<std::size_t>(numThreads));
    std::size_t dimension = query.RunCloud(best, static_cast<std::size_t>(lgMaxSample),
        box, volume);
    int32_t rawAgreement = RawAgreement(box, volume, rawBox, rawVolume);
    io.integer(rawAgreement, rawAgreement);

    io.outInt(static_cast<int32_t>(dimension));
    OutBox(io, box, volume);
    io.outBool(Contains(box, best));
    });
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.computeHull.canonical")
{
    RunWithBigStack([&io]() {
    // The vertices-and-indices query (overloads 3 and 4) on convex hull
    // meshes presented in a scrambled triangle order, with the orders of
    // the hash containers pinned. Besides the box, the whole observable
    // state is compared: the edges, the adjacency lists after the coplanar
    // merge, the climb start, the sample subdivision, the aligned candidate
    // and every field of the winning candidate (edge pair, level-curve
    // processor, normals, bilinear corner values, axes, support indices,
    // volume). Same exact-support restriction as compute.canonical.
    int32_t mode = io.integer(0, 4);
    int32_t lgMaxSample = io.integer(2, 4);
    int32_t n = io.integer(5, 10);
    int32_t numThreads = io.integer(0, 4);

    std::vector<Vec3> bestVertices{};
    std::vector<int32_t> bestIndices{};
    double bestViolation = std::numeric_limits<double>::max();
    for (int32_t attempt = 0; attempt < 24; ++attempt)
    {
        std::vector<Vec3> points = DrawCloud(io, mode, n);
        std::vector<Vec3> vertices{};
        std::vector<int32_t> indices{};
        if (!MakeHullMesh(io, points, vertices, indices))
        {
            continue;
        }
        OrientedBox3<double> box{};
        double volume = 0.0;
        CanonicalFP probe(0);
        try
        {
            probe.RunHull(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);
        }
        catch (std::exception const&)
        {
            continue;
        }
        double violation = ContainmentViolation(box, vertices) / PointScale(vertices);
        bool sound = SupportsAreExact(probe);
        if (sound || violation < bestViolation)
        {
            bestViolation = violation;
            bestVertices = vertices;
            bestIndices = indices;
        }
        if (sound)
        {
            break;
        }
    }
    if (bestVertices.empty())
    {
        MakeHullMesh(io, FallbackCloud(), bestVertices, bestIndices);
    }
    GiveMesh(io, bestVertices, bestIndices);

    OrientedBox3<double> rawBox{};
    double rawVolume = 0.0;
    MVB3FP raw(static_cast<std::size_t>(numThreads));
    if (io.index() % 2 == 0)
    {
        raw(bestVertices.size(), bestVertices.data(), bestIndices.size(), bestIndices.data(),
            static_cast<std::size_t>(lgMaxSample), rawBox, rawVolume);
    }
    else
    {
        raw(bestVertices, bestIndices, static_cast<std::size_t>(lgMaxSample), rawBox, rawVolume);
    }

    OrientedBox3<double> box{};
    double volume = 0.0;
    CanonicalFP query(static_cast<std::size_t>(numThreads));
    query.RunHull(bestVertices, bestIndices, static_cast<std::size_t>(lgMaxSample), box, volume);
    int32_t rawAgreement = RawAgreement(box, volume, rawBox, rawVolume);
    io.integer(rawAgreement, rawAgreement);

    OutBox(io, box, volume);
    io.outBool(Contains(box, bestVertices));
    OutTopology(io, query);
    OutAlignedFP(io, query);
    OutWinnerFP(io, query);
    });
}

// @@APPEND@@
