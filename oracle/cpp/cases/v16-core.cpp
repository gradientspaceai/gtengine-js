// Group 16 (core): CurveExtractorSquares, CurveExtractorTriangles,
// IEEEBinary16, MeshSmoother, PolygonTree, PolygonWindingOrder,
// ETNonmanifoldMesh, ImplicitSurface3, MeshCurvature, RevolutionMesh,
// TubeMesh, VETNonmanifoldMesh, VertexCollapseMesh.
//
// Everything is integer or bit manipulation, container bookkeeping or IEEE
// arithmetic with sqrt, and is compared bit for bit, except where a case
// comment names a C math library call (the cos/sin tables of RevolutionMesh
// and TubeMesh, the std:: math wrappers of IEEEBinary16).
//
// Functions.h is compiled with optimization off (the v15 lesson: MSVC /O2
// compiles clamp's 'x >= xmax ? xmax : x' to minsd, which returns the wrong
// zero on a signed-zero tie); IEEEBinary16's clamp wrapper reaches it.
//
// Matrix2x2.h and Matrix3x3.h are included explicitly, so the dependent
// calls Inverse(M) of Mesh::UpdateFrame (2x2) and MeshCurvature (3x3)
// resolve to the closed-form overloads, which is what the port implements
// (ORACLE.md, v34/v43).
#define ORACLE_FAMILY "v16-core"
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <utility>
#include <vector>
#pragma optimize("", off)
#include <Mathematics/Functions.h>
#pragma optimize("", on)
#include "Oracle.h"

#include <Mathematics/BitHacks.h>
#include <Mathematics/IEEEBinary.h>
// Convert32To16/Convert16To32 are private upstream (the port exposes them
// as static methods). Every header IEEEBinary16.h includes is already
// included, so the macro reaches IEEEBinary16.h alone.
#define private public
#include <Mathematics/IEEEBinary16.h>
#undef private

#include <Mathematics/Matrix2x2.h>
#include <Mathematics/Matrix3x3.h>
#include <Mathematics/Vector2.h>
#include <Mathematics/Vector3.h>
#include <Mathematics/CurveExtractorSquares.h>
#include <Mathematics/CurveExtractorTriangles.h>
#include <Mathematics/PolygonTree.h>
#include <Mathematics/PolygonWindingOrder.h>
#include <Mathematics/MeshSmoother.h>
#include <Mathematics/MeshCurvature.h>
#include <Mathematics/ImplicitSurface3.h>
#include <Mathematics/ETNonmanifoldMesh.h>
#include <Mathematics/VETNonmanifoldMesh.h>
#include <Mathematics/BezierCurve.h>
#include <Mathematics/RevolutionMesh.h>
#include <Mathematics/TubeMesh.h>

// VertexCollapseMesh keeps its heap, mesh and helpers private; the
// canonical-order driver below replays DoCollapse through them. Every header
// it includes comes first, so the macro reaches VertexCollapseMesh.h alone.
#include <Mathematics/MinHeap.h>
#include <Mathematics/Polygon2.h>
#include <Mathematics/TriangulateEC.h>
#include <Mathematics/VETManifoldMesh.h>
#define private public
#include <Mathematics/VertexCollapseMesh.h>
#undef private

using namespace gte;

// ================================================================ CurveExtractor

namespace
{
    // An image of size xBound x yBound in one of five modes:
    //   0: values and level in [-2, 2] (many zero corners: every one of the
    //      81 sign patterns of a square, the '0000' and '+-+-' saddle cases)
    //   1: values in [0, 255], level in [0, 255] (the uint8_t image type)
    //   2: values over the int16_t range
    //   3: a '+-+-' checkerboard aimed at the saddle determinant: each value
    //      is +-(1..6) so that f00*f11 == f10*f01 (det == 0) happens often
    //   4: values up to 2^20 in magnitude (int32_t; the port's number
    //      arithmetic is exact while |f| * bound stays below 2^53).
    struct ImageDraw
    {
        int32_t xBound, yBound;
        int mode;
        std::vector<int64_t> pixels;
        int64_t level;
    };

    ImageDraw DrawImage(oracle::Ctx& io)
    {
        ImageDraw d{};
        d.mode = io.integer(0, 4);
        d.xBound = io.integer(2, 6);
        d.yBound = io.integer(2, 6);
        size_t n = static_cast<size_t>(d.xBound) * static_cast<size_t>(d.yBound);
        d.pixels.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            int64_t v = 0;
            switch (d.mode)
            {
            case 0: v = io.rawInteger(-2, 2); break;
            case 1: v = io.rawInteger(0, 255); break;
            case 2: v = io.rawInteger(-32768, 32767); break;
            case 3:
            {
                int32_t x = static_cast<int32_t>(i) % d.xBound;
                int32_t y = static_cast<int32_t>(i) / d.xBound;
                int64_t s = ((x + y) % 2 == 0 ? 1 : -1);
                int64_t m = io.rawInteger(1, 6);
                v = s * m;
                break;
            }
            default: v = io.rawInteger(-(1 << 20), 1 << 20); break;
            }
            d.pixels[i] = static_cast<int64_t>(io.given(static_cast<double>(v)));
        }
        int64_t level = 0;
        switch (d.mode)
        {
        case 0: level = io.rawInteger(-2, 2); break;
        case 1: level = io.rawInteger(0, 255); break;
        case 2: level = io.rawInteger(-32768, 32767); break;
        case 3: level = io.rawInteger(-1, 1); break;
        default: level = io.rawInteger(-(1 << 20), 1 << 20); break;
        }
        d.level = static_cast<int64_t>(io.given(static_cast<double>(level)));
        return d;
    }

    template <typename Vertex, typename Edge>
    void EmitRational(oracle::Ctx& io, std::vector<Vertex> const& vertices,
        std::vector<Edge> const& edges)
    {
        io.outInt(vertices.size());
        for (auto const& v : vertices)
        {
            io.outInt(v.xNumer); io.outInt(v.xDenom);
            io.outInt(v.yNumer); io.outInt(v.yDenom);
        }
        io.outInt(edges.size());
        for (auto const& e : edges) { io.outInt(e.v[0]); io.outInt(e.v[1]); }
    }

    template <typename Edge>
    void EmitReal(oracle::Ctx& io, std::vector<std::array<double, 2>> const& vertices,
        std::vector<Edge> const& edges)
    {
        io.outInt(vertices.size());
        for (auto const& v : vertices) { io.outReal(v[0]); io.outReal(v[1]); }
        io.outInt(edges.size());
        for (auto const& e : edges) { io.outInt(e.v[0]); io.outInt(e.v[1]); }
    }

    // Rational extraction, MakeUnique of it, and the real-valued Extract
    // with and without duplicate removal.
    template <template <typename, typename> class Extractor, typename T>
    void RunExtractor(oracle::Ctx& io, ImageDraw const& d)
    {
        std::vector<T> image(d.pixels.size());
        for (size_t i = 0; i < image.size(); ++i) { image[i] = static_cast<T>(d.pixels[i]); }
        Extractor<T, double> extractor(d.xBound, d.yBound, image.data());
        using Ex = Extractor<T, double>;
        std::vector<typename Ex::Vertex> vertices;
        std::vector<typename Ex::Edge> edges;
        extractor.Extract(static_cast<T>(d.level), vertices, edges);
        EmitRational(io, vertices, edges);
        extractor.MakeUnique(vertices, edges);
        EmitRational(io, vertices, edges);
        for (int remove = 0; remove < 2; ++remove)
        {
            std::vector<std::array<double, 2>> rv;
            std::vector<typename Ex::Edge> re;
            // The derived Extract hides the base overload; call it through
            // the base class, as a caller holding a CurveExtractor does.
            CurveExtractor<T, double>& base = extractor;
            base.Extract(static_cast<T>(d.level), remove != 0, rv, re);
            EmitReal(io, rv, re);
        }
    }

    template <template <typename, typename> class Extractor>
    void RunExtractorByMode(oracle::Ctx& io)
    {
        ImageDraw d = DrawImage(io);
        switch (d.mode)
        {
        case 1: RunExtractor<Extractor, uint8_t>(io, d); break;
        case 2: RunExtractor<Extractor, int16_t>(io, d); break;
        default: RunExtractor<Extractor, int32_t>(io, d); break;
        }
    }
}

ORACLE_CASE("CurveExtractorSquares.extract")
{
    RunExtractorByMode<CurveExtractorSquares>(io);
}

ORACLE_CASE("CurveExtractorTriangles.extract")
{
    RunExtractorByMode<CurveExtractorTriangles>(io);
}
