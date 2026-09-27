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

        // RunHull up to and including RemoveCoplanarTriangleAdjacencies, the
        // state in which GetExtreme is called.
        void Prepare(std::vector<Vec3> const& vertices, std::vector<int32_t> const& indices,
            std::size_t lgMaxSample)
        {
            this->GenerateSubdivision(lgMaxSample);
            VETManifoldMesh mesh{};
            this->CreateMeshTopology(indices.size() / 3, indices.data(), mesh);
            ExtractMeshTopologySorted(mesh);
            ExtractVertexAdjacenciesSorted(mesh);
            this->ExtractMeshGeometry(vertices.size(), vertices.data());
            this->RemoveCoplanarTriangleAdjacencies();
        }

        template <typename Direction, typename Value>
        std::size_t Extreme(Direction const& direction, Value& dMax)
        {
            return this->GetExtreme(direction, dMax);
        }
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

    // Independent reference checks of a returned 3D box: the axes are
    // orthonormal and the volume is the product of the side lengths (both
    // to 1e-9, the box being rounded from rationals).
    bool BoxIsConsistent(OrientedBox3<double> const& box, double volume)
    {
        for (int32_t i = 0; i < 3; ++i)
        {
            for (int32_t j = 0; j < 3; ++j)
            {
                double expected = (i == j ? 1.0 : 0.0);
                if (std::fabs(Dot(box.axis[i], box.axis[j]) - expected) > 1e-9)
                {
                    return false;
                }
            }
        }
        double product = 8.0 * box.extent[0] * box.extent[1] * box.extent[2];
        return std::fabs(product - volume) <= 1e-9 * std::max(1.0, std::fabs(volume));
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

    // The largest amount, over the six supports of the winning candidate,
    // by which upstream's support vertex falls short of the exact extreme
    // minus the port's rounding bound. Positive exactly when one of the
    // port's fixes replaces a support of upstream's winner.
    double WinnerSupportExcess(CanonicalFP const& query)
    {
        auto const& c = query.Winner();
        auto const& T = query.TVertices();
        double excess = -std::numeric_limits<double>::max();
        for (int32_t i = 0; i < 3; ++i)
        {
            Vec3 const& a = c.axis[i];
            double tolerance = 8.0 * std::numeric_limits<double>::epsilon()
                * std::max(std::fabs(a[0]), std::max(std::fabs(a[1]), std::fabs(a[2])));
            double maxL1 = 0.0;
            for (auto const& v : T)
            {
                maxL1 = std::max(maxL1, std::fabs(v[0]) + std::fabs(v[1]) + std::fabs(v[2]));
            }
            tolerance *= maxL1;
            ExactNumber dmin = ExactDot(a, T[c.minSupportIndex[i]]);
            ExactNumber dmax = ExactDot(a, T[c.maxSupportIndex[i]]);
            ExactNumber lo = dmin, hi = dmax;
            for (auto const& v : T)
            {
                ExactNumber d = ExactDot(a, v);
                if (d < lo) { lo = d; }
                if (d > hi) { hi = d; }
            }
            excess = std::max(excess, static_cast<double>(dmin - lo) - tolerance);
            excess = std::max(excess, static_cast<double>(hi - dmax) - tolerance);
        }
        return excess;
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
    // Diagnostic input, ignored by the replay: 1 when the recorded cloud
    // passed the exact-support separator (0 only for a capped fallback).
    int32_t sound = (SupportsAreExact(query) ? 1 : 0);
    io.integer(sound, sound);

    io.outInt(static_cast<int32_t>(dimension));
    OutBox(io, box, volume);
    io.outBool(Contains(box, best));
    io.outBool(BoxIsConsistent(box, volume));
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
    int32_t sound = (SupportsAreExact(query) ? 1 : 0);
    io.integer(sound, sound);

    OutBox(io, box, volume);
    io.outBool(Contains(box, bestVertices));
    io.outBool(BoxIsConsistent(box, volume));
    OutTopology(io, query);
    OutAlignedFP(io, query);
    OutWinnerFP(io, query);
    });
}

namespace
{
    // A search direction for GetExtreme. Kind 0: a signed coordinate axis;
    // kind 1: a lattice direction in [-2,2]^3 (edges and faces of lattice
    // hulls are often perpendicular to it, so ties are exact); kind 2: a
    // uniform unit vector; kind 3: the normalized normal of a hull triangle
    // (upstream's candidate axes are such normals), which makes the whole
    // face a plateau up to rounding. Nothing is recorded here.
    Vec3 DrawDirection(oracle::Ctx& io, int32_t kind, std::vector<Vec3> const& vertices,
        std::vector<int32_t> const& indices)
    {
        Vec3 d{ 0.0, 0.0, 0.0 };
        if (kind == 0)
        {
            d[io.rawInteger(0, 2)] = (io.rawInteger(0, 1) == 0 ? -1.0 : 1.0);
        }
        else if (kind == 1)
        {
            do
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    d[j] = static_cast<double>(io.rawInteger(-2, 2));
                }
            } while (d[0] == 0.0 && d[1] == 0.0 && d[2] == 0.0);
        }
        else if (kind == 2)
        {
            double length = 0.0;
            do
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    d[j] = io.raw(-1.0, 1.0);
                }
                length = Length(d);
            } while (length < 0.1 || length > 1.0);
            Normalize(d);
        }
        else
        {
            std::size_t t = static_cast<std::size_t>(io.rawInteger(0,
                static_cast<int32_t>(indices.size() / 3) - 1));
            Vec3 const& p0 = vertices[indices[3 * t]];
            Vec3 const& p1 = vertices[indices[3 * t + 1]];
            Vec3 const& p2 = vertices[indices[3 * t + 2]];
            d = Cross(p2 - p0, p1 - p0);
            if (io.rawInteger(0, 1) == 0)
            {
                d = -d;
            }
            Normalize(d);
        }
        return d;
    }

    // The port's rounding bound (MinimumVolumeBox3FloatingPoint.ts,
    // climbTolerance): 8 * eps * max_i |d[i]| * max_v L1(mTVertices[v]).
    double ClimbTolerance(std::vector<Vec3> const& T, Vec3 const& d)
    {
        double maxL1 = 0.0;
        for (auto const& v : T)
        {
            maxL1 = std::max(maxL1, std::fabs(v[0]) + std::fabs(v[1]) + std::fabs(v[2]));
        }
        double maxD = std::max(std::fabs(d[0]), std::max(std::fabs(d[1]), std::fabs(d[2])));
        return 8.0 * std::numeric_limits<double>::epsilon() * maxL1 * maxD;
    }

    // How far upstream's GetExtreme vertex falls short of the exact maximum
    // of Dot(d, mTVertices[v]): 0 when it is an exact maximizer.
    double ExtremeShortfall(std::vector<Vec3> const& T, Vec3 const& d, std::size_t vUpstream)
    {
        ExactNumber best = ExactDot(d, T[vUpstream]);
        ExactNumber reached = best;
        for (auto const& v : T)
        {
            ExactNumber value = ExactDot(d, v);
            if (value > best)
            {
                best = value;
            }
        }
        return static_cast<double>(best - reached);
    }
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.getExtreme")
{
    RunWithBigStack([&io]() {
    // GetExtreme, the hill climb every candidate uses, called directly (via
    // the Canonical subclass) on a hull mesh after the coplanar merge, in
    // four directions per record. Restricted to meshes and directions where
    // upstream's vertex is an exact maximizer (the separator of #426; the
    // direction is redrawn up to 16 times, the mesh up to 8 times). The
    // rotated-lattice mode is left to the plateau deviation case.
    int32_t const modes[4] = { 0, 1, 3, 4 };
    int32_t mode = modes[io.integer(0, 3)];
    int32_t n = io.integer(5, 10);

    std::vector<Vec3> vertices{};
    std::vector<int32_t> indices{};
    for (int32_t attempt = 0; attempt < 8 && vertices.empty(); ++attempt)
    {
        MakeHullMesh(io, DrawCloud(io, mode, n), vertices, indices);
    }
    if (vertices.empty())
    {
        MakeHullMesh(io, FallbackCloud(), vertices, indices);
    }
    GiveMesh(io, vertices, indices);

    CanonicalFP query(0);
    query.Prepare(vertices, indices, 2);
    for (int32_t k = 0; k < 4; ++k)
    {
        int32_t kind = io.integer(0, 3);
        Vec3 direction{};
        for (int32_t attempt = 0; attempt < 16; ++attempt)
        {
            direction = DrawDirection(io, kind, vertices, indices);
            double dMax = 0.0;
            std::size_t v = query.Extreme(direction, dMax);
            if (ExtremeShortfall(query.TVertices(), direction, v) == 0.0)
            {
                break;
            }
        }
        io.givenVec(direction);
        double dMax = 0.0;
        std::size_t vMax = query.Extreme(direction, dMax);
        io.outInt(IndexOut(vMax));
        io.outReal(dMax);
    }
    });
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.getExtreme.plateau")
{
    RunWithBigStack([&io]() {
    // Deviation (#426): upstream's strict-improvement climb stops at a
    // vertex that falls short of the exact maximum by more than the rounding
    // bound, so the port's plateau traversal replaces it. The construction
    // follows the finding: a lattice box whose FIRST point lies in the
    // interior of a face (so it becomes hull vertex 0, the climb start),
    // rotated by 0.03 to 0.3 degrees so that the face is coplanar only up to
    // rounding and the exact coplanar merge leaves the point in the graph;
    // the directions are the signed normals of the hull triangles at vertex
    // 0 and the signed axes of upstream's own winning candidate. Capped at 48
    // clouds; each record keeps the direction with the largest shortfall.
    int32_t n = io.integer(6, 10);
    std::vector<Vec3> bestVertices{};
    std::vector<int32_t> bestIndices{};
    Vec3 bestDirection{ 1.0, 0.0, 0.0 };
    double bestExcess = -1.0;
    for (int32_t attempt = 0; attempt < 48 && bestExcess <= 0.0; ++attempt)
    {
        double sx = static_cast<double>(io.rawInteger(2, 4));
        double sy = static_cast<double>(io.rawInteger(2, 4));
        double sz = static_cast<double>(io.rawInteger(1, 4));
        std::vector<Vec3> points{};
        points.push_back(Vec3{ static_cast<double>(io.rawInteger(1, static_cast<int32_t>(sx) - 1)),
            static_cast<double>(io.rawInteger(1, static_cast<int32_t>(sy) - 1)), sz });
        for (int32_t k = 0; k < 8; ++k)
        {
            points.push_back(Vec3{ (k & 1) ? sx : 0.0, (k & 2) ? sy : 0.0, (k & 4) ? sz : 0.0 });
        }
        for (int32_t i = 9; i < n; ++i)
        {
            points.push_back(Vec3{ static_cast<double>(io.rawInteger(0, static_cast<int32_t>(sx))),
                static_cast<double>(io.rawInteger(0, static_cast<int32_t>(sy))),
                static_cast<double>(io.rawInteger(0, static_cast<int32_t>(sz))) });
        }
        double angle = io.raw(0.0005, 0.005);
        double c = std::cos(angle), s = std::sin(angle);
        int32_t axis = io.rawInteger(0, 1);
        for (auto& p : points)
        {
            double u = p[axis], v = p[2];
            p[axis] = c * u - s * v;
            p[2] = s * u + c * v;
        }

        std::vector<Vec3> vertices{};
        std::vector<int32_t> indices{};
        if (!MakeHullMesh(io, points, vertices, indices))
        {
            continue;
        }
        std::vector<Vec3> directions{};
        for (std::size_t t = 0; t < indices.size() / 3; ++t)
        {
            if (indices[3 * t] == 0 || indices[3 * t + 1] == 0 || indices[3 * t + 2] == 0)
            {
                Vec3 const& p0 = vertices[indices[3 * t]];
                Vec3 const& p1 = vertices[indices[3 * t + 1]];
                Vec3 const& p2 = vertices[indices[3 * t + 2]];
                Vec3 normal = Cross(p2 - p0, p1 - p0);
                Normalize(normal);
                directions.push_back(normal);
                directions.push_back(-normal);
            }
        }
        {
            CanonicalFP run(0);
            OrientedBox3<double> box{};
            double volume = 0.0;
            try
            {
                run.RunHull(vertices, indices, 2, box, volume);
                for (int32_t i = 0; i < 3; ++i)
                {
                    Vec3 axis = run.Winner().axis[i];
                    directions.push_back(axis);
                    directions.push_back(-axis);
                }
            }
            catch (std::exception const&)
            {
            }
        }
        CanonicalFP probe(0);
        probe.Prepare(vertices, indices, 2);
        for (auto const& direction : directions)
        {
            double dMax = 0.0;
            std::size_t v = probe.Extreme(direction, dMax);
            double excess = ExtremeShortfall(probe.TVertices(), direction, v)
                - ClimbTolerance(probe.TVertices(), direction);
            if (excess > bestExcess)
            {
                bestExcess = excess;
                bestVertices = vertices;
                bestIndices = indices;
                bestDirection = direction;
            }
        }
    }
    if (bestVertices.empty())
    {
        MakeHullMesh(io, FallbackCloud(), bestVertices, bestIndices);
    }
    GiveMesh(io, bestVertices, bestIndices);
    io.givenVec(bestDirection);
    // Diagnostic input, ignored by the replay: 1 when the shortfall of
    // upstream's vertex exceeds the rounding bound (the port must differ).
    int32_t found = (bestExcess > 0.0 ? 1 : 0);
    io.integer(found, found);

    CanonicalFP query(0);
    query.Prepare(bestVertices, bestIndices, 2);
    double dMax = 0.0;
    std::size_t vMax = query.Extreme(bestDirection, dMax);
    io.outInt(IndexOut(vMax));
    io.outReal(dMax);
    });
}

namespace
{
    // The four virtual minimizers are upstream's documented customization
    // point. This subclass records every call (which minimizer, at which
    // level-curve processor, with which s/t arguments) and can replace any
    // of them by a no-op ('disabled' bit mask: 1 MinimizerConstantS, 2
    // MinimizerConstantT, 4 MinimizerVariableS, 8 MinimizerVariableT), which
    // checks that the port dispatches through the same overridable methods
    // from the same processors. Single-threaded (the counters are not
    // synchronized).
    class MinimizerProbeFP : public CanonicalFP
    {
    public:
        MinimizerProbeFP(int32_t inDisabled) : CanonicalFP(0), disabled(inDisabled) {}

        int32_t disabled;
        std::array<std::size_t, 4> calls{};
        std::array<double, 6> argSum{};
        std::map<std::size_t, std::size_t> processors{};

    protected:
        void MinimizerConstantS(Candidate& c, Candidate& mvc) override
        {
            ++calls[0];
            ++processors[c.levelCurveProcessorIndex];
            if ((disabled & 1) == 0)
            {
                CanonicalFP::MinimizerConstantS(c, mvc);
            }
        }

        void MinimizerConstantT(Candidate& c, Candidate& mvc) override
        {
            ++calls[1];
            ++processors[c.levelCurveProcessorIndex];
            if ((disabled & 2) == 0)
            {
                CanonicalFP::MinimizerConstantT(c, mvc);
            }
        }

        void MinimizerVariableS(double const& sminNumer, double const& smaxNumer,
            double const& sDenom, Candidate& c, Candidate& mvc) override
        {
            ++calls[2];
            ++processors[c.levelCurveProcessorIndex];
            argSum[0] += sminNumer;
            argSum[1] += smaxNumer;
            argSum[2] += sDenom;
            if ((disabled & 4) == 0)
            {
                CanonicalFP::MinimizerVariableS(sminNumer, smaxNumer, sDenom, c, mvc);
            }
        }

        void MinimizerVariableT(double const& tminNumer, double const& tmaxNumer,
            double const& tDenom, Candidate& c, Candidate& mvc) override
        {
            ++calls[3];
            ++processors[c.levelCurveProcessorIndex];
            argSum[3] += tminNumer;
            argSum[4] += tmaxNumer;
            argSum[5] += tDenom;
            if ((disabled & 8) == 0)
            {
                CanonicalFP::MinimizerVariableT(tminNumer, tmaxNumer, tDenom, c, mvc);
            }
        }
    };

    template <typename Probe>
    void OutMinimizerCalls(oracle::Ctx& io, Probe const& probe)
    {
        for (int32_t i = 0; i < 4; ++i)
        {
            io.outInt(static_cast<int32_t>(probe.calls[i]));
        }
        for (int32_t i = 0; i < 6; ++i)
        {
            io.outReal(static_cast<double>(probe.argSum[i]));
        }
        io.outInt(static_cast<int32_t>(probe.processors.size()));
        for (auto const& element : probe.processors)
        {
            io.outInt(static_cast<int32_t>(element.first));
            io.outInt(static_cast<int32_t>(element.second));
        }
    }
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.minimizers")
{
    RunWithBigStack([&io]() {
    // The minimizer overrides (see MinimizerProbeFP) on canonical hull
    // meshes: the call counts per minimizer, the level-curve processors that
    // reached them, the sums of the s/t numerators and denominators passed
    // to the variable minimizers, and the resulting box and winner. Same
    // exact-support restriction as the main cases, evaluated on the probe
    // with the same minimizers disabled.
    int32_t mode = io.integer(0, 4);
    int32_t lgMaxSample = io.integer(2, 4);
    int32_t n = io.integer(5, 10);
    int32_t disabled = io.integer(0, 15);

    std::vector<Vec3> bestVertices{};
    std::vector<int32_t> bestIndices{};
    double bestViolation = std::numeric_limits<double>::max();
    for (int32_t attempt = 0; attempt < 24; ++attempt)
    {
        std::vector<Vec3> vertices{};
        std::vector<int32_t> indices{};
        if (!MakeHullMesh(io, DrawCloud(io, mode, n), vertices, indices))
        {
            continue;
        }
        OrientedBox3<double> box{};
        double volume = 0.0;
        MinimizerProbeFP probe(disabled);
        probe.RunHull(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);
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

    OrientedBox3<double> box{};
    double volume = 0.0;
    MinimizerProbeFP query(disabled);
    query.RunHull(bestVertices, bestIndices, static_cast<std::size_t>(lgMaxSample), box, volume);
    OutMinimizerCalls(io, query);
    OutBox(io, box, volume);
    OutWinnerFP(io, query);
    });
}

namespace
{
    // ---- the dimension-2 path (both headers, identical code) -------------
    //
    // A coplanar cloud. Mode 0: a lattice plane perpendicular to a
    // coordinate axis; mode 1: a tilted lattice plane origin + a*d0 + b*d1;
    // mode 2: an axis-perpendicular plane with uniform in-plane coordinates;
    // mode 3: an axis-perpendicular plane whose lattice points are rotated
    // in-plane by the exact Pythagorean rotation (3/5, 4/5) and scaled by 5,
    // so that the minimum rectangle is tilted and its arithmetic exact.
    std::vector<Vec3> DrawPlanarCloud(oracle::Ctx& io, int32_t mode, int32_t n)
    {
        std::vector<Vec3> points{};
        if (mode == 1)
        {
            Vec3 origin{}, d0{}, d1{};
            for (int32_t j = 0; j < 3; ++j)
            {
                origin[j] = static_cast<double>(io.rawInteger(-3, 3));
            }
            do
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    d0[j] = static_cast<double>(io.rawInteger(-2, 2));
                    d1[j] = static_cast<double>(io.rawInteger(-2, 2));
                }
            } while (Length(Cross(d0, d1)) == 0.0);
            for (int32_t i = 0; i < n; ++i)
            {
                double a = static_cast<double>(io.rawInteger(-3, 3));
                double b = static_cast<double>(io.rawInteger(-3, 3));
                points.push_back(origin + a * d0 + b * d1);
            }
            return points;
        }

        int32_t axis = io.rawInteger(0, 2);
        int32_t a0 = (axis + 1) % 3, a1 = (axis + 2) % 3;
        double level = static_cast<double>(io.rawInteger(-3, 3));
        for (int32_t i = 0; i < n; ++i)
        {
            Vec3 p{};
            p[axis] = level;
            if (mode == 2)
            {
                p[a0] = io.raw(-4.0, 4.0);
                p[a1] = io.raw(-4.0, 4.0);
            }
            else
            {
                double u = static_cast<double>(io.rawInteger(-3, 3));
                double v = static_cast<double>(io.rawInteger(-3, 3));
                p[a0] = (mode == 3 ? 3.0 * u - 4.0 * v : u);
                p[a1] = (mode == 3 ? 4.0 * u + 3.0 * v : v);
            }
            points.push_back(p);
        }
        return points;
    }

    // The two probes of the dimension-2 path. 'basisSame': upstream's
    // orthonormal basis (Newell loop from i1 = 1, #352/#355) equals the
    // port's (loop from i1 = 0) bit for bit. 'rectangleSame': upstream's
    // MinimumAreaBox2<double, double> rectangle of the projected points
    // equals MinimumAreaBox2<double, BSRational<UIntegerAP32>>, the exact
    // instantiation the port implements (bit-identical to it per v11), bit
    // for bit. Both use upstream's own code; the dimension is 2 on return
    // true.
    bool Dimension2Probe(std::vector<Vec3> const& points, bool& basisSame, bool& rectangleSame,
        bool& liftedSame)
    {
        ConvexHull3<double> ch3{};
        ch3(points.size(), points.data(), 0);
        if (ch3.GetDimension() != 2)
        {
            return false;
        }
        auto const& hull = ch3.GetHull();
        std::array<std::array<Vec3, 3>, 2> basis{};
        for (std::size_t k = 0; k < 2; ++k)
        {
            Vec3 normal = Vec3::Zero();
            std::size_t numHull = hull.size();
            for (std::size_t i0 = numHull - 1, i1 = (k == 0 ? 1 : 0); i1 < numHull; i0 = i1++)
            {
                normal += Cross(points[hull[i0]], points[hull[i1]]);
            }
            basis[k][0] = normal;
            ComputeOrthogonalComplement(1, basis[k].data());
        }
        basisSame = true;
        for (std::size_t i = 0; i < 3; ++i)
        {
            for (int32_t j = 0; j < 3; ++j)
            {
                basisSame = basisSame && SameBits(basis[0][i][j], basis[1][i][j]);
            }
        }

        Vec3 origin = points[hull[0]];
        std::vector<Vector2<double>> projection(points.size());
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            Vec3 diff = points[i] - origin;
            projection[i][0] = Dot(basis[0][1], diff);
            projection[i][1] = Dot(basis[0][2], diff);
        }
        MinimumAreaBox2<double, double> floatBox{};
        OrientedBox2<double> r0 = floatBox(static_cast<int32_t>(projection.size()),
            projection.data());
        MinimumAreaBox2<double, BSRational<UIntegerAP32>> exactBox{};
        OrientedBox2<double> r1 = exactBox(static_cast<int32_t>(projection.size()),
            projection.data());
        rectangleSame = true;
        for (int32_t i = 0; i < 2; ++i)
        {
            rectangleSame = rectangleSame && SameBits(r0.center[i], r1.center[i])
                && SameBits(r0.extent[i], r1.extent[i]);
            for (int32_t j = 0; j < 2; ++j)
            {
                rectangleSame = rectangleSame && SameBits(r0.axis[i][j], r1.axis[i][j]);
            }
        }

        // 'liftedSame': the two rectangles lifted into 3D with upstream's
        // lifting expressions give the same box bit for bit. A difference
        // can vanish there (the sign of a zero, or an ulp absorbed by the
        // additions), so the deviation case selects on this observable.
        auto lift = [&origin, &basis](OrientedBox2<double> const& r)
        {
            std::array<Vec3, 3> b = basis[0];
            OrientedBox3<double> box{};
            box.center = origin + r.center[0] * b[1] + r.center[1] * b[2];
            box.axis[0] = r.axis[0][0] * b[1] + r.axis[0][1] * b[2];
            box.axis[1] = r.axis[1][0] * b[1] + r.axis[1][1] * b[2];
            box.axis[2] = b[0];
            box.extent = { r.extent[0], r.extent[1], 0.0 };
            return box;
        };
        liftedSame = SameBox(lift(r0), 0.0, lift(r1), 0.0);
        return true;
    }

    // The dimension-2 record of either header, through upstream's raw
    // operator() (overload 1 or 2; nothing on this path depends on a hash
    // order). 'wantFloatDeviation' selects the deviation generator: basis
    // equal but rectangle different.
    template <typename Query>
    void Dimension2Record(oracle::Ctx& io, bool wantFloatDeviation)
    {
        int32_t mode = io.integer(0, 3);
        int32_t lgMaxSample = io.integer(2, 3);
        int32_t n = io.integer(4, 9);

        // Capped at 32 draws, each redrawing the whole cloud. The main case
        // falls back to a fixed tilted lattice rectangle (below), the
        // deviation case to the last draw whose basis agrees. The diagnostic
        // input the replay ignores is 1 for an accepted draw, 2 for the
        // fallback when it passes both probes, 0 otherwise.
        std::vector<Vec3> best{};
        int32_t accepted = 0;
        for (int32_t attempt = 0; attempt < 32 && accepted == 0; ++attempt)
        {
            std::vector<Vec3> points = DrawPlanarCloud(io, mode, n);
            bool basisSame = false, rectangleSame = false, liftedSame = false;
            if (!Dimension2Probe(points, basisSame, rectangleSame, liftedSame)
                || !basisSame)
            {
                continue;
            }
            // The main case needs the rectangles bit for bit equal, the
            // deviation the lifted 3D boxes different (the observable).
            if (wantFloatDeviation ? !liftedSame : rectangleSame)
            {
                best = points;
                accepted = 1;
            }
            else if (wantFloatDeviation)
            {
                best = points;
            }
        }
        if (best.empty())
        {
            // The lattice rectangle (0,0),(2,0),(2,1),(0,1) plus (1,1) mapped
            // by (u,v) -> (3u-4v, 4u+3v) onto the plane z = 1: both probes
            // pass (the diagnostic records 2 when they do). An axis-aligned
            // rectangle does not do: the floating-point
            // MinimumAreaBox2 returns a -0 axis component where the exact one
            // returns +0.
            best = { Vec3{ 0.0, 0.0, 1.0 }, Vec3{ 6.0, 8.0, 1.0 }, Vec3{ 2.0, 11.0, 1.0 },
                Vec3{ -4.0, 3.0, 1.0 }, Vec3{ -1.0, 7.0, 1.0 } };
            bool basisSame = false, rectangleSame = false, liftedSame = false;
            Dimension2Probe(best, basisSame, rectangleSame, liftedSame);
            accepted = (basisSame && rectangleSame ? 2 : 0);
        }
        GiveCloud(io, best);
        io.integer(accepted, accepted);

        OrientedBox3<double> box{};
        double volume = -1.0;
        Query query(0);
        std::size_t dimension = (io.index() % 2 == 0
            ? query(best.size(), best.data(), static_cast<std::size_t>(lgMaxSample), box, volume)
            : query(best, static_cast<std::size_t>(lgMaxSample), box, volume));
        io.outInt(static_cast<int32_t>(dimension));
        OutBox(io, box, volume);
    }
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.compute.dimension2")
{
    RunWithBigStack([&io]() {
    // Coplanar clouds on which both deliberate differences of the
    // dimension-2 path are inert (Dimension2Probe): upstream's Newell basis
    // equals the corrected one bit for bit (the unclosed loop sums the
    // Newell normal of the hull polygon with hull[0] removed, which for four
    // or more hull vertices has the right direction; its normalization then
    // often agrees, always for an axis-perpendicular plane), and the
    // floating-point MinimumAreaBox2 rectangle equals the exact one.
    Dimension2Record<MVB3FP>(io, false);
    });
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.compute.dimension2.floatComputeType")
{
    RunWithBigStack([&io]() {
    // Deviation (design, inherited from MinimumAreaBox2): the basis agrees,
    // but upstream's MinimumAreaBox2<T, T> floating-point rectangle differs
    // from the exact one the port computes.
    Dimension2Record<MVB3FP>(io, true);
    });
}

namespace
{
    // ---- the rational header ---------------------------------------------
    //
    // Everything is exact until GetMinimumVolumeBox rounds the box, so
    // neither #405 nor #426 exists here. The one deliberate port fix is
    // #355: MinimizerVariableT declares its parameters 'T const&', so
    // upstream samples every t-variable level curve at parameters rounded to
    // double while the port samples them exactly. MinimizerProbeR counts the
    // calls; the main cases reject records that reach MinimizerVariableT, and
    // the minimizer case always replaces it (its rounded arguments are still
    // compared).
    class MinimizerProbeR : public CanonicalR
    {
    public:
        MinimizerProbeR(int32_t inDisabled) : CanonicalR(0), disabled(inDisabled) {}

        int32_t disabled;
        std::array<std::size_t, 4> calls{};
        std::array<double, 6> argSum{};
        std::map<std::size_t, std::size_t> processors{};

    protected:
        void MinimizerConstantS(Candidate& c, Candidate& mvc) override
        {
            ++calls[0];
            ++processors[c.levelCurveProcessorIndex];
            if ((disabled & 1) == 0)
            {
                CanonicalR::MinimizerConstantS(c, mvc);
            }
        }

        void MinimizerConstantT(Candidate& c, Candidate& mvc) override
        {
            ++calls[1];
            ++processors[c.levelCurveProcessorIndex];
            if ((disabled & 2) == 0)
            {
                CanonicalR::MinimizerConstantT(c, mvc);
            }
        }

        void MinimizerVariableS(Number const& sminNumer, Number const& smaxNumer,
            Number const& sDenom, Candidate& c, Candidate& mvc) override
        {
            ++calls[2];
            ++processors[c.levelCurveProcessorIndex];
            argSum[0] += static_cast<double>(sminNumer);
            argSum[1] += static_cast<double>(smaxNumer);
            argSum[2] += static_cast<double>(sDenom);
            if ((disabled & 4) == 0)
            {
                CanonicalR::MinimizerVariableS(sminNumer, smaxNumer, sDenom, c, mvc);
            }
        }

        // Upstream's signature (#355): the arguments arrive rounded to double.
        void MinimizerVariableT(double const& tminNumer, double const& tmaxNumer,
            double const& tDenom, Candidate& c, Candidate& mvc) override
        {
            ++calls[3];
            ++processors[c.levelCurveProcessorIndex];
            argSum[3] += tminNumer;
            argSum[4] += tmaxNumer;
            argSum[5] += tDenom;
            if ((disabled & 8) == 0)
            {
                CanonicalR::MinimizerVariableT(tminNumer, tmaxNumer, tDenom, c, mvc);
            }
        }
    };

    template <typename Vector3Number>
    void OutNVector(oracle::Ctx& io, Vector3Number const& v)
    {
        for (int32_t j = 0; j < 3; ++j)
        {
            io.outReal(static_cast<double>(v[j]));
        }
    }

    template <typename Query>
    void OutAlignedR(oracle::Ctx& io, Query const& query)
    {
        auto const& c = query.Aligned();
        for (int32_t i = 0; i < 3; ++i)
        {
            io.outInt(IndexOut(c.minSupportIndex[i]));
            io.outInt(IndexOut(c.maxSupportIndex[i]));
        }
        io.outReal(static_cast<double>(c.volume));
    }

    // The winner's exact fields, each rounded to double once.
    template <typename Query>
    void OutWinnerR(oracle::Ctx& io, Query const& query)
    {
        auto const& c = query.Winner();
        io.outInt(IndexOut(c.edgeIndex[0]));
        io.outInt(IndexOut(c.edgeIndex[1]));
        io.outInt(IndexOut(c.levelCurveProcessorIndex));
        for (int32_t i = 0; i < 2; ++i)
        {
            OutNVector(io, c.N[i]);
            OutNVector(io, c.M[i]);
        }
        io.outReal(static_cast<double>(c.f00));
        io.outReal(static_cast<double>(c.f10));
        io.outReal(static_cast<double>(c.f01));
        io.outReal(static_cast<double>(c.f11));
        for (int32_t i = 0; i < 3; ++i)
        {
            OutNVector(io, c.axis[i]);
            io.outInt(IndexOut(c.minSupportIndex[i]));
            io.outInt(IndexOut(c.maxSupportIndex[i]));
        }
        io.outReal(static_cast<double>(c.volume));
    }

    // Draw a hull mesh (or a cloud, when 'vertices' is the cloud itself) on
    // which the canonical rational run never reaches MinimizerVariableT.
    // Capped at 8 attempts; falls back to the attempt with the fewest calls.
    template <typename Make>
    void DrawWithoutVariableT(Make const& make, int32_t lgMaxSample,
        std::vector<Vec3>& bestVertices, std::vector<int32_t>& bestIndices)
    {
        std::size_t bestCalls = std::numeric_limits<std::size_t>::max();
        for (int32_t attempt = 0; attempt < 8 && bestCalls > 0; ++attempt)
        {
            std::vector<Vec3> vertices{};
            std::vector<int32_t> indices{};
            if (!make(vertices, indices))
            {
                continue;
            }
            OrientedBox3<double> box{};
            double volume = 0.0;
            MinimizerProbeR probe(0);
            probe.RunHull(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);
            if (probe.calls[3] < bestCalls)
            {
                bestCalls = probe.calls[3];
                bestVertices = vertices;
                bestIndices = indices;
            }
        }
    }
}

ORACLE_CASE("MinimumVolumeBox3Rational.computeHull.canonical")
{
    RunWithBigStack([&io]() {
    // The rational vertices-and-indices query (overloads 3 and 4) on
    // scrambled hull meshes of all five cloud modes, lattice ones included
    // (v08 had to exclude lattice clouds: their many exactly tied candidates
    // are decided by the hash order, which the control pins). Full box and
    // observable state compared, as in the floating-point sibling. Records
    // that reach MinimizerVariableT (#355) are rejected.
    int32_t mode = io.integer(0, 4);
    int32_t lgMaxSample = io.integer(2, 3);
    int32_t n = io.integer(5, 7);
    int32_t numThreads = io.integer(0, 3);

    std::vector<Vec3> vertices{};
    std::vector<int32_t> indices{};
    DrawWithoutVariableT([&io, mode, n](std::vector<Vec3>& v, std::vector<int32_t>& t)
        { return MakeHullMesh(io, DrawCloud(io, mode, n), v, t); },
        lgMaxSample, vertices, indices);
    if (vertices.empty())
    {
        MakeHullMesh(io, FallbackCloud(), vertices, indices);
    }
    GiveMesh(io, vertices, indices);

    OrientedBox3<double> rawBox{};
    double rawVolume = 0.0;
    MVB3R raw(static_cast<std::size_t>(numThreads));
    if (io.index() % 2 == 0)
    {
        raw(vertices.size(), vertices.data(), indices.size(), indices.data(),
            static_cast<std::size_t>(lgMaxSample), rawBox, rawVolume);
    }
    else
    {
        raw(vertices, indices, static_cast<std::size_t>(lgMaxSample), rawBox, rawVolume);
    }

    OrientedBox3<double> box{};
    double volume = 0.0;
    CanonicalR query(static_cast<std::size_t>(numThreads));
    query.RunHull(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);
    int32_t rawAgreement = RawAgreement(box, volume, rawBox, rawVolume);
    io.integer(rawAgreement, rawAgreement);

    OutBox(io, box, volume);
    io.outBool(Contains(box, vertices));
    io.outBool(BoxIsConsistent(box, volume));
    OutTopology(io, query);
    OutAlignedR(io, query);
    OutWinnerR(io, query);
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.compute.canonical")
{
    RunWithBigStack([&io]() {
    // The rational point-cloud query (overloads 1 and 2), canonical order,
    // full box bit for bit; replaces v08's volume comparison with tolerance
    // 1.5e-1. All five cloud modes; records that reach MinimizerVariableT
    // (#355) are rejected (capped at 8 attempts, falling back to the cloud
    // with the fewest calls).
    int32_t mode = io.integer(0, 4);
    int32_t lgMaxSample = io.integer(2, 3);
    int32_t n = io.integer(5, 7);
    int32_t numThreads = io.integer(0, 3);

    std::vector<Vec3> best{};
    std::size_t bestCalls = std::numeric_limits<std::size_t>::max();
    for (int32_t attempt = 0; attempt < 8 && bestCalls > 0; ++attempt)
    {
        std::vector<Vec3> points = DrawCloud(io, mode, n);
        OrientedBox3<double> box{};
        double volume = 0.0;
        MinimizerProbeR probe(0);
        if (probe.RunCloud(points, static_cast<std::size_t>(lgMaxSample), box, volume) != 3)
        {
            continue;
        }
        if (probe.calls[3] < bestCalls)
        {
            bestCalls = probe.calls[3];
            best = points;
        }
    }
    if (best.empty())
    {
        best = FallbackCloud();
    }
    GiveCloud(io, best);

    OrientedBox3<double> rawBox{};
    double rawVolume = 0.0;
    MVB3R raw(static_cast<std::size_t>(numThreads));
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
    CanonicalR query(static_cast<std::size_t>(numThreads));
    std::size_t dimension = query.RunCloud(best, static_cast<std::size_t>(lgMaxSample),
        box, volume);
    int32_t rawAgreement = RawAgreement(box, volume, rawBox, rawVolume);
    io.integer(rawAgreement, rawAgreement);

    io.outInt(static_cast<int32_t>(dimension));
    OutBox(io, box, volume);
    io.outBool(Contains(box, best));
    io.outBool(BoxIsConsistent(box, volume));
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.getExtreme")
{
    RunWithBigStack([&io]() {
    // The exact GetExtreme on hull meshes of all five modes (the rotated
    // mode included: the dot products are exact, so there is no plateau
    // defect), four directions of the four kinds per record. No restriction.
    int32_t mode = io.integer(0, 4);
    int32_t n = io.integer(5, 10);

    std::vector<Vec3> vertices{};
    std::vector<int32_t> indices{};
    for (int32_t attempt = 0; attempt < 8 && vertices.empty(); ++attempt)
    {
        MakeHullMesh(io, DrawCloud(io, mode, n), vertices, indices);
    }
    if (vertices.empty())
    {
        MakeHullMesh(io, FallbackCloud(), vertices, indices);
    }
    GiveMesh(io, vertices, indices);

    CanonicalR query(0);
    query.Prepare(vertices, indices, 2);
    for (int32_t k = 0; k < 4; ++k)
    {
        int32_t kind = io.integer(0, 3);
        Vec3 direction = DrawDirection(io, kind, vertices, indices);
        io.givenVec(direction);
        MVB3R::NVector3 nDirection{};
        for (int32_t j = 0; j < 3; ++j)
        {
            nDirection[j] = MVB3R::Number(direction[j]);
        }
        MVB3R::Number dMax{};
        std::size_t vMax = query.Extreme(nDirection, dMax);
        io.outInt(IndexOut(vMax));
        io.outReal(static_cast<double>(dMax));
    }
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.minimizers")
{
    RunWithBigStack([&io]() {
    // The minimizer overrides of the rational header. MinimizerVariableT is
    // always replaced (bit 8 set): its body is #355, but its callers and its
    // double-rounded arguments are compared exactly. The other three are
    // disabled at random.
    int32_t mode = io.integer(0, 4);
    int32_t lgMaxSample = io.integer(2, 3);
    int32_t n = io.integer(5, 7);
    int32_t disabled = 8 | io.integer(0, 7);

    std::vector<Vec3> vertices{};
    std::vector<int32_t> indices{};
    for (int32_t attempt = 0; attempt < 8 && vertices.empty(); ++attempt)
    {
        MakeHullMesh(io, DrawCloud(io, mode, n), vertices, indices);
    }
    if (vertices.empty())
    {
        MakeHullMesh(io, FallbackCloud(), vertices, indices);
    }
    GiveMesh(io, vertices, indices);

    OrientedBox3<double> box{};
    double volume = 0.0;
    MinimizerProbeR query(disabled);
    query.RunHull(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);
    OutMinimizerCalls(io, query);
    OutBox(io, box, volume);
    OutWinnerR(io, query);
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.compute.dimension2")
{
    RunWithBigStack([&io]() {
    // The rational header's dimension-2 path is the floating-point one
    // verbatim (T = double, MinimumAreaBox2<T, T>); same probes.
    Dimension2Record<MVB3R>(io, false);
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.compute.dimension2.floatComputeType")
{
    RunWithBigStack([&io]() {
    // Deviation (design), as in the floating-point sibling.
    Dimension2Record<MVB3R>(io, true);
    });
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.computeHull.provenViolation")
{
    RunWithBigStack([&io]() {
    // Deviation (#405, #426): hull meshes on which a support vertex of
    // upstream's winning candidate falls short of the exact extreme by more
    // than the port's rounding bound (WinnerSupportExcess > 0), so the
    // port's confined fix replaces it. Selected by the exact separator
    // rather than by containment (v08 nonContaining), so it includes
    // shortfalls too small to fail the 1e-9 containment test (the output
    // boolean). Uniform and rotated-lattice clouds, capped at 64 attempts;
    // the diagnostic input the replay ignores is 1 when the search succeeded.
    int32_t lgMaxSample = io.integer(2, 4);
    int32_t n = io.integer(6, 10);

    std::vector<Vec3> bestVertices{};
    std::vector<int32_t> bestIndices{};
    double bestExcess = -std::numeric_limits<double>::max();
    for (int32_t attempt = 0; attempt < 64 && bestExcess <= 0.0; ++attempt)
    {
        std::vector<Vec3> vertices{};
        std::vector<int32_t> indices{};
        if (!MakeHullMesh(io, DrawCloud(io, (attempt % 2 == 0 ? 1 : 2), n), vertices, indices))
        {
            continue;
        }
        OrientedBox3<double> box{};
        double volume = 0.0;
        CanonicalFP probe(0);
        probe.RunHull(vertices, indices, static_cast<std::size_t>(lgMaxSample), box, volume);
        double excess = WinnerSupportExcess(probe);
        if (excess > bestExcess)
        {
            bestExcess = excess;
            bestVertices = vertices;
            bestIndices = indices;
        }
    }
    GiveMesh(io, bestVertices, bestIndices);
    int32_t found = (bestExcess > 0.0 ? 1 : 0);
    io.integer(found, found);

    OrientedBox3<double> box{};
    double volume = 0.0;
    CanonicalFP query(0);
    query.RunHull(bestVertices, bestIndices, static_cast<std::size_t>(lgMaxSample), box, volume);
    OutBox(io, box, volume);
    io.outBool(Contains(box, bestVertices));
    });
}

namespace
{
    // Reusing one functor for two meshes. CreateMeshTopology resizes mEdges
    // but only reserves mEdgeIndices, and ExtractMeshTopology appends, so the
    // second call processes the first mesh's edge pairs, read as indices into
    // the second mesh's edges, before its own. When the second mesh has fewer
    // edges those indices are out of range (undefined behaviour; not
    // exercised). Otherwise the stale pairs are duplicates processed first:
    // every minimizer they reach is called again (observable through the
    // minimizer overrides, the documented customization point), and when
    // candidates tie exactly the stale order can change the winner. The port
    // resets mEdgeIndices, which is a fresh functor's behaviour.
    //
    // The record runs a MinimizerProbe functor on mesh A, clears its
    // counters, runs it on mesh B (at least as many triangles, hence edges,
    // as A) and emits the counters and the box. Lattice clouds; for the
    // floating-point header mesh B must satisfy the exact-support separator
    // and for the rational one MinimizerVariableT is replaced (#355), so that
    // only the stale pairs separate upstream from the port. The diagnostic
    // input the replay ignores has bit 1 set when the stale pairs changed
    // the box and bit 2 when they changed the minimizer calls, compared
    // with a fresh functor.
    template <typename Probe, typename Sound>
    void ReuseRecord(oracle::Ctx& io, int32_t lgMaxSample, int32_t disabled,
        Sound const& sound)
    {
        std::vector<Vec3> vA{}, vB{};
        std::vector<int32_t> iA{}, iB{};
        for (int32_t attempt = 0; attempt < 32; ++attempt)
        {
            std::vector<Vec3> a{}, b{};
            std::vector<int32_t> ia{}, ib{};
            if (!MakeHullMesh(io, DrawCloud(io, 0, io.rawInteger(4, 6)), a, ia)
                || !MakeHullMesh(io, DrawCloud(io, 0, io.rawInteger(5, 8)), b, ib)
                || ib.size() < ia.size() || !sound(b, ib))
            {
                continue;
            }
            vA = a; iA = ia; vB = b; iB = ib;
            break;
        }
        if (vA.empty())
        {
            MakeHullMesh(io, FallbackCloud(), vA, iA);
            vB = vA;
            iB = iA;
        }
        GiveMesh(io, vA, iA);
        GiveMesh(io, vB, iB);

        std::size_t const lg = static_cast<std::size_t>(lgMaxSample);
        OrientedBox3<double> box{}, fresh{};
        double volume = 0.0, freshVolume = 0.0;
        Probe once(disabled);
        once.RunHull(vB, iB, lg, fresh, freshVolume);
        Probe query(disabled);
        query.RunHull(vA, iA, lg, box, volume);
        query.calls = {};
        query.argSum = {};
        query.processors.clear();
        query.RunHull(vB, iB, lg, box, volume);
        int32_t changed = (SameBox(box, volume, fresh, freshVolume) ? 0 : 1)
            + (query.calls == once.calls && query.argSum == once.argSum
                && query.processors == once.processors ? 0 : 2);
        io.integer(changed, changed);

        OutMinimizerCalls(io, query);
        OutBox(io, box, volume);
    }
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.computeHull.reuse")
{
    RunWithBigStack([&io]() {
    // Deviation (new upstream suspect, see the group report).
    int32_t lgMaxSample = io.integer(2, 3);
    ReuseRecord<MinimizerProbeFP>(io, lgMaxSample, 0,
        [lgMaxSample](std::vector<Vec3> const& v, std::vector<int32_t> const& i)
        {
            OrientedBox3<double> box{};
            double volume = 0.0;
            CanonicalFP probe(0);
            probe.RunHull(v, i, static_cast<std::size_t>(lgMaxSample), box, volume);
            return SupportsAreExact(probe);
        });
    });
}

ORACLE_CASE("MinimumVolumeBox3Rational.computeHull.reuse")
{
    RunWithBigStack([&io]() {
    // The same in the rational header, MinimizerVariableT replaced.
    int32_t lgMaxSample = io.integer(2, 2);
    ReuseRecord<MinimizerProbeR>(io, lgMaxSample, 8,
        [](std::vector<Vec3> const&, std::vector<int32_t> const&) { return true; });
    });
}

namespace
{
    // Throw parity for the LogAssert preconditions of the four operator()
    // overloads, through upstream's raw functor (the asserts fire before any
    // order-dependent step). Kind 0: fewer than 4 points; 1: lgMaxSample < 2
    // on a valid cloud; 2: fewer than 12 indices or a count that is not a
    // multiple of 3; 3: fewer than 4 vertices; 4: lgMaxSample < 2 on a valid
    // mesh. Every record throws.
    template <typename Query>
    void InvalidArgumentRecord(oracle::Ctx& io)
    {
        int32_t kind = io.integer(0, 4);
        int32_t overload = io.integer(0, 1);
        int32_t lgMaxSample = (kind == 1 || kind == 4 ? io.integer(0, 1) : io.integer(2, 3));
        std::vector<Vec3> points{};
        std::vector<int32_t> indices{};
        if (kind <= 1)
        {
            int32_t n = (kind == 0 ? io.integer(1, 3) : io.integer(4, 6));
            for (int32_t i = 0; i < n; ++i)
            {
                points.push_back(io.latticeVec<3>(-3, 3));
            }
        }
        else
        {
            std::vector<Vec3> vertices{};
            MakeHullMesh(io, FallbackCloud(), vertices, indices);
            if (kind == 3)
            {
                vertices.resize(3);
            }
            else if (kind == 2)
            {
                indices.resize(indices.size() - static_cast<std::size_t>(io.rawInteger(1, 2)));
                if (io.rawInteger(0, 1) == 0)
                {
                    indices.resize(9);
                }
            }
            GiveMesh(io, vertices, indices);
            points = vertices;
        }

        OrientedBox3<double> box{};
        double volume = 0.0;
        Query query(0);
        std::size_t lg = static_cast<std::size_t>(lgMaxSample);
        if (kind <= 1)
        {
            std::size_t dimension = (overload == 0
                ? query(points.size(), points.data(), lg, box, volume)
                : query(points, lg, box, volume));
            io.outInt(static_cast<int32_t>(dimension));
        }
        else if (overload == 0)
        {
            query(points.size(), points.data(), indices.size(), indices.data(), lg, box, volume);
        }
        else
        {
            query(points, indices, lg, box, volume);
        }
        OutBox(io, box, volume);
    }
}

ORACLE_CASE("MinimumVolumeBox3FloatingPoint.invalidArgument")
{
    RunWithBigStack([&io]() { InvalidArgumentRecord<MVB3FP>(io); });
}

ORACLE_CASE("MinimumVolumeBox3Rational.invalidArgument")
{
    RunWithBigStack([&io]() { InvalidArgumentRecord<MVB3R>(io); });
}

// @@APPEND@@
