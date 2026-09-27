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

    ImageDraw DrawImageRaw(oracle::Ctx& io, int forcedMode)
    {
        ImageDraw d{};
        d.mode = (forcedMode >= 0 ? forcedMode : io.rawInteger(0, 4));
        d.xBound = io.rawInteger(2, 6);
        d.yBound = io.rawInteger(2, 6);
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
            d.pixels[i] = v;
        }
        switch (d.mode)
        {
        case 0: d.level = io.rawInteger(-2, 2); break;
        case 1: d.level = io.rawInteger(0, 255); break;
        case 2: d.level = io.rawInteger(-32768, 32767); break;
        case 3: d.level = io.rawInteger(-1, 1); break;
        default: d.level = io.rawInteger(-(1 << 20), 1 << 20); break;
        }
        return d;
    }

    // The number of squares whose corner values F = pixel - level have the
    // sign pattern '+000' or '00+0' (up to a global sign) in upstream's
    // (f00, f10, f11, f01) order: exactly one nonzero corner, at (x,y) or at
    // (x+1,y+1). CurveExtractorSquares emits the wrong pair of square edges
    // for those two patterns (the port fixes it, see the deviation case).
    int CountThreeZeroSquares(ImageDraw const& d)
    {
        int count = 0;
        for (int32_t y = 0; y + 1 < d.yBound; ++y)
        {
            for (int32_t x = 0; x + 1 < d.xBound; ++x)
            {
                int32_t i00 = x + d.xBound * y;
                int64_t f00 = d.pixels[i00] - d.level;
                int64_t f10 = d.pixels[i00 + 1] - d.level;
                int64_t f01 = d.pixels[i00 + d.xBound] - d.level;
                int64_t f11 = d.pixels[i00 + d.xBound + 1] - d.level;
                if (f10 == 0 && f01 == 0 && ((f00 != 0) != (f11 != 0)))
                {
                    ++count;
                }
            }
        }
        return count;
    }

    ImageDraw GivenImage(oracle::Ctx& io, ImageDraw const& d)
    {
        io.given(static_cast<double>(d.mode));
        io.given(static_cast<double>(d.xBound));
        io.given(static_cast<double>(d.yBound));
        for (int64_t p : d.pixels) { io.given(static_cast<double>(p)); }
        io.given(static_cast<double>(d.level));
        return d;
    }

    // 'avoidThreeZero': redraw (at most 64 times, then fall back to an
    // int16_t image, where a zero corner has probability 2^-16) until no
    // square has the '+000'/'00+0' pattern, the inputs where upstream's
    // squares extractor is sound.
    ImageDraw DrawImage(oracle::Ctx& io, bool avoidThreeZero)
    {
        ImageDraw d = DrawImageRaw(io, -1);
        for (int attempt = 0; avoidThreeZero && CountThreeZeroSquares(d) > 0; ++attempt)
        {
            d = DrawImageRaw(io, attempt < 64 ? -1 : 2);
        }
        return GivenImage(io, d);
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
    void RunExtractorOnImage(oracle::Ctx& io, ImageDraw const& d)
    {
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
    // Restricted to images without a '+000' or '00+0' square (see
    // CurveExtractorSquares.extract.threeZeroCorners).
    ImageDraw d = DrawImage(io, true);
    RunExtractorOnImage<CurveExtractorSquares>(io, d);
}

ORACLE_CASE("CurveExtractorSquares.extract.threeZeroCorners")
{
    // Deliberate port fix (new finding of this group): for a square whose
    // only nonzero corner is (x,y) ('+000') or (x+1,y+1) ('00+0'),
    // upstream emits the two square edges incident to that corner, whose
    // points other than the far corners are not on the level set; the
    // port emits the two zero edges. A mode-0 image (values in [-2, 2]) with
    // one such square planted at a random position; every record contains
    // at least one, so every record deviates.
    ImageDraw d = DrawImageRaw(io, 0);
    int32_t sx = io.rawInteger(0, d.xBound - 2);
    int32_t sy = io.rawInteger(0, d.yBound - 2);
    int32_t i00 = sx + d.xBound * sy;
    int64_t nonzero = d.level + (io.rawInteger(0, 1) != 0 ? 1 : -1) * io.rawInteger(1, 2);
    bool atOrigin = (io.rawInteger(0, 1) != 0);
    d.pixels[i00] = (atOrigin ? nonzero : d.level);
    d.pixels[i00 + 1] = d.level;
    d.pixels[i00 + d.xBound] = d.level;
    d.pixels[i00 + d.xBound + 1] = (atOrigin ? d.level : nonzero);
    GivenImage(io, d);
    RunExtractorOnImage<CurveExtractorSquares>(io, d);
}

ORACLE_CASE("CurveExtractorTriangles.extract")
{
    ImageDraw d = DrawImage(io, false);
    RunExtractorOnImage<CurveExtractorTriangles>(io, d);
}

// ================================================================ IEEEBinary16
//
// Encodings travel as integers (a 32-bit pattern is below 2^53), never as
// doubles, so that no NaN payload is lost. Real outputs are the binary32
// values of half numbers widened to double (exact).

namespace
{
    uint32_t RawBits32(oracle::Ctx& io)
    {
        uint32_t hi = static_cast<uint32_t>(io.rawInteger(0, 0xFFFF));
        uint32_t lo = static_cast<uint32_t>(io.rawInteger(0, 0xFFFF));
        return (hi << 16) | lo;
    }

    // A binary32 pattern aimed at every branch of Convert32To16.
    uint32_t DrawConvertBits(oracle::Ctx& io)
    {
        uint32_t const sign = (io.rawInteger(0, 1) != 0 ? 0x80000000u : 0u);
        int mode = io.rawInteger(0, 9);
        switch (mode)
        {
        case 0:  // zero, 32-subnormal
            return sign | (io.rawInteger(0, 1) != 0 ? (RawBits32(io) & 0x007FFFFFu) : 0u);
        case 1:  // the five class thresholds and their neighbours
        {
            static uint32_t const thresholds[5] =
                { 0x33000000u, 0x33800000u, 0x38800000u, 0x477FE000u, 0x477FF000u };
            uint32_t t = thresholds[io.rawInteger(0, 4)];
            return sign | static_cast<uint32_t>(static_cast<int64_t>(t) + io.rawInteger(-2, 2));
        }
        case 2:  // 16-normal range, the 13 dropped bits at, below or above half
        {
            uint32_t biased = static_cast<uint32_t>(io.rawInteger(113, 142));
            uint32_t t16 = static_cast<uint32_t>(io.rawInteger(0, 0x3FF));
            static uint32_t const lows[6] = { 0x1000u, 0x0FFFu, 0x1001u, 0u, 0x1FFFu, 0x0001u };
            int which = io.rawInteger(0, 6);
            uint32_t low = (which < 6 ? lows[which] : (RawBits32(io) & 0x1FFFu));
            return sign | (biased << 23) | (t16 << 13) | low;
        }
        case 3:  // 16-subnormal range, the dropped bits at, below or above half
        {
            uint32_t biased = static_cast<uint32_t>(io.rawInteger(103, 112));
            int32_t p = static_cast<int32_t>(biased) - 127 + 15;
            uint32_t rshift = static_cast<uint32_t>(14 - p);
            uint32_t high = RawBits32(io) & 0x007FFFFFu & ~((1u << rshift) - 1u);
            int64_t delta = io.rawInteger(-1, 1);
            uint32_t low = static_cast<uint32_t>(static_cast<int64_t>(1u << (rshift - 1)) + delta);
            return sign | (biased << 23) | high | low;
        }
        case 4:  // signaling NaN, payload below 2^13 (finding #110) or above
        {
            uint32_t payload = (io.rawInteger(0, 1) != 0
                ? static_cast<uint32_t>(io.rawInteger(1, 0x1FFF))
                : (RawBits32(io) & 0x003FFFFFu) | 0x2000u);
            return sign | 0x7F800000u | payload;
        }
        case 5:  // quiet NaN, default or with payload
            return sign | 0x7FC00000u | (io.rawInteger(0, 1) != 0 ? (RawBits32(io) & 0x003FFFFFu) : 0u);
        case 6:  // infinity, float max
            return sign | (io.rawInteger(0, 1) != 0 ? 0x7F800000u : 0x7F7FFFFFu);
        case 7:  // above 16-max-normal
            return sign | (static_cast<uint32_t>(io.rawInteger(143, 254)) << 23) | (RawBits32(io) & 0x007FFFFFu);
        case 8:  // below 2^-25
            return sign | (static_cast<uint32_t>(io.rawInteger(1, 101)) << 23) | (RawBits32(io) & 0x007FFFFFu);
        default:
            return RawBits32(io);
        }
    }

    uint32_t GivenBits(oracle::Ctx& io, uint32_t u)
    {
        io.given(static_cast<double>(u));
        return u;
    }

    float FloatOfBits(uint32_t u)
    {
        float f;
        std::memcpy(&f, &u, sizeof(f));
        return f;
    }

    uint32_t BitsOfFloat(float f)
    {
        uint32_t u;
        std::memcpy(&u, &f, sizeof(u));
        return u;
    }

    // A half encoding: every class, the class boundaries, random patterns.
    uint16_t DrawHalfBits(oracle::Ctx& io)
    {
        uint16_t const sign = (io.rawInteger(0, 1) != 0 ? 0x8000u : 0u);
        int mode = io.rawInteger(0, 7);
        switch (mode)
        {
        case 0: return static_cast<uint16_t>(sign | static_cast<uint16_t>(io.rawInteger(0, 2)));
        case 1: return static_cast<uint16_t>(sign | static_cast<uint16_t>(0x03FF - io.rawInteger(0, 1)));
        case 2: return static_cast<uint16_t>(sign | static_cast<uint16_t>(0x0400 + io.rawInteger(0, 1)));
        case 3: return static_cast<uint16_t>(sign | static_cast<uint16_t>(0x7BFF - io.rawInteger(0, 1)));
        case 4: return static_cast<uint16_t>(sign | 0x7C00u);
        case 5: return static_cast<uint16_t>(sign | 0x7C00u | static_cast<uint16_t>(io.rawInteger(1, 0x3FF)));
        default: return static_cast<uint16_t>(io.rawInteger(0, 0xFFFF));
        }
    }

    uint16_t GivenHalf(oracle::Ctx& io, uint16_t h)
    {
        io.given(static_cast<double>(h));
        return h;
    }
}

ORACLE_CASE("IEEEBinary16.convert32To16")
{
    // Four patterns per record through the private Convert32To16 and
    // through the float constructor (which the port has no counterpart for
    // on signaling NaNs, so its encoding is emitted for non-NaN inputs).
    for (int k = 0; k < 4; ++k)
    {
        uint32_t bits = GivenBits(io, DrawConvertBits(io));
        uint16_t h = IEEEBinary16::Convert32To16(bits);
        io.outInt(h);
        if ((bits & 0x7FFFFFFFu) <= 0x7F800000u)
        {
            IEEEBinary16 x(FloatOfBits(bits));
            io.outInt(static_cast<uint16_t>(x));
        }
    }
}

ORACLE_CASE("IEEEBinary16.convert16To32")
{
    // Eight encodings per record: the binary32 pattern of Convert16To32,
    // the operator float() and operator double() values (bits of the float
    // for NaNs, since the harness matches any NaN with any NaN).
    for (int k = 0; k < 8; ++k)
    {
        uint16_t h = GivenHalf(io, DrawHalfBits(io));
        uint32_t bits = IEEEBinary16::Convert16To32(h);
        io.outInt(bits);
        IEEEBinary16 x(h);
        float f = static_cast<float>(x);
        double d = static_cast<double>(x);
        io.outInt(BitsOfFloat(f));
        if (d == d) { io.outReal(d); }
    }
}

ORACLE_CASE("IEEEBinary16.fields")
{
    // Every accessor, classifier and neighbour function of the base
    // IEEEBinary<int16_t, uint16_t, 16, 11>, and SetEncoding with in-range
    // and out-of-range fields (uint16_t parameters; the shifted value is
    // truncated to 16 bits on assignment).
    uint16_t h = GivenHalf(io, DrawHalfBits(io));
    IEEEBinary16 x(h);
    io.outInt(static_cast<uint16_t>(x));
    io.outInt(x.GetSign());
    io.outInt(x.GetBiased());
    io.outInt(x.GetTrailing());
    uint16_t sign = 0, biased = 0, trailing = 0;
    x.GetEncoding(sign, biased, trailing);
    io.outInt(sign);
    io.outInt(biased);
    io.outInt(trailing);
    io.outInt(static_cast<int32_t>(x.GetClassification()));
    io.outBool(x.IsZero());
    io.outBool(x.IsSignMinus());
    io.outBool(x.IsSubnormal());
    io.outBool(x.IsNormal());
    io.outBool(x.IsFinite());
    io.outBool(x.IsInfinite());
    io.outBool(x.IsNaN());
    io.outBool(x.IsQuietNaN());
    io.outBool(x.IsSignalingNaN());
    io.outInt(x.GetNextUp());
    io.outInt(x.GetNextDown());
    int mode = io.rawInteger(0, 3);
    uint16_t rs = static_cast<uint16_t>(io.rawInteger(0, 1));
    uint16_t rb = static_cast<uint16_t>(io.rawInteger(0, 31));
    uint16_t rt = static_cast<uint16_t>(io.rawInteger(0, 0x3FF));
    if (mode == 1) { rs = static_cast<uint16_t>(io.rawInteger(2, 7)); }
    if (mode == 2) { rb = static_cast<uint16_t>(io.rawInteger(32, 0xFFFF)); }
    if (mode == 3) { rt = static_cast<uint16_t>(io.rawInteger(0x400, 0xFFFF)); }
    uint16_t s = GivenHalf(io, rs);
    uint16_t b = GivenHalf(io, rb);
    uint16_t t = GivenHalf(io, rt);
    IEEEBinary16 y(static_cast<uint16_t>(0));
    y.SetEncoding(s, b, t);
    io.outInt(static_cast<uint16_t>(y));
}

ORACLE_CASE("IEEEBinary16.fromNumber")
{
    // IEEEBinary16(double): static_cast<float> first (round to nearest
    // even, NaN quieted keeping the high payload bits), then Convert32To16.
    // Doubles on and beside binary32 ties, half ties, both class
    // thresholds, NaNs with payloads, infinities, signed zeros. Read in the
    // replay with io.real() (a NaN keeps its bits).
    int mode = io.rawInteger(0, 4);
    double d = 0.0;
    if (mode == 0)
    {
        // A binary32 value from the Convert32To16 generator, or 1 binary64
        // ulp beside it, or halfway to its binary32 neighbour (a double
        // rounding tie).
        uint32_t f = DrawConvertBits(io);
        double base = static_cast<double>(FloatOfBits(f));
        uint64_t b;
        std::memcpy(&b, &base, sizeof(b));
        int shift = io.rawInteger(0, 2);
        if ((f & 0x7F800000u) != 0x7F800000u)
        {
            if (shift == 1) { b += static_cast<uint64_t>(io.rawInteger(-1, 1)); }
            if (shift == 2) { b += (1ull << 28); }
        }
        std::memcpy(&d, &b, sizeof(d));
    }
    else if (mode == 1)
    {
        uint64_t b = (static_cast<uint64_t>(RawBits32(io)) << 32) | RawBits32(io);
        b = (b & 0x800FFFFFFFFFFFFFull) | 0x7FF0000000000000ull;
        if ((b & 0x000FFFFFFFFFFFFFull) == 0) { b |= 1; }
        std::memcpy(&d, &b, sizeof(d));
    }
    else if (mode == 2)
    {
        d = std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(-30, 17));
        if (io.rawInteger(0, 1) != 0) { d = -d; }
    }
    else if (mode == 3)
    {
        d = (io.rawInteger(0, 1) != 0 ? -0.0 : 0.0);
    }
    else
    {
        d = io.raw(-70000.0, 70000.0);
    }
    d = io.given(d);
    IEEEBinary16 x(d);
    io.outInt(static_cast<uint16_t>(x));
}

ORACLE_CASE("IEEEBinary16.compare")
{
    // The six comparison operators, which compare the binary32 values:
    // -0 == +0, NaN unordered, subnormals and neighbours.
    uint16_t a = GivenHalf(io, DrawHalfBits(io));
    int mode = io.rawInteger(0, 3);
    uint16_t rb = DrawHalfBits(io);
    if (mode == 0) { rb = a; }
    if (mode == 1) { rb = static_cast<uint16_t>(a ^ 0x8000u); }
    if (mode == 2) { rb = IEEEBinary16(a).GetNextUp(); }
    uint16_t b = GivenHalf(io, rb);
    IEEEBinary16 x(a), y(b);
    io.outBool(x == y);
    io.outBool(x != y);
    io.outBool(x < y);
    io.outBool(x <= y);
    io.outBool(x > y);
    io.outBool(x >= y);
}

ORACLE_CASE("IEEEBinary16.arithmetic")
{
    // Unary minus (a sign flip in binary16), the binary operators (binary32
    // results) with half and float operands, and the compound updates
    // (binary32 result converted back to binary16). A float operand is a
    // binary32 value recorded as a double.
    uint16_t a = GivenHalf(io, DrawHalfBits(io));
    uint16_t b = GivenHalf(io, DrawHalfBits(io));
    float f = FloatOfBits(DrawConvertBits(io) & 0xFF7FFFFFu);
    f = static_cast<float>(io.given(static_cast<double>(f)));
    IEEEBinary16 x(a), y(b);
    io.outInt(static_cast<uint16_t>(-x));
    // NaN results: the sign and payload of a NaN produced by x64 SSE
    // arithmetic (default NaN 0xFFC00000, operand NaN quieted) is not
    // reproducible from JavaScript, whose NaN arithmetic results carry no
    // defined bits; a NaN result is emitted as a real (any NaN matches any
    // NaN) and a NaN half as its NaN flag.
    float r[12] = { x + y, x - y, x * y, x / y, x + f, x - f, x * f, x / f,
        f + x, f - x, f * x, f / x };
    for (float v : r) { io.outReal(static_cast<double>(v)); }
    IEEEBinary16 u[8] = { x, x, x, x, x, x, x, x };
    u[0] += y; u[1] -= y; u[2] *= y; u[3] /= y;
    u[4] += f; u[5] -= f; u[6] *= f; u[7] /= f;
    for (auto const& v : u)
    {
        io.outBool(v.IsNaN());
        if (!v.IsNaN()) { io.outInt(static_cast<uint16_t>(v)); }
    }
}

namespace
{
    // A half for the math wrappers: any encoding, a value in [-4, 4], or a
    // multiple of 1/8 in [-2, 2] (exact domain boundaries of acos, asin,
    // atanh, the saturate and clamp ties, integers for ceil/floor/fmod).
    uint16_t DrawMathHalf(oracle::Ctx& io)
    {
        int mode = io.rawInteger(0, 2);
        uint16_t h = 0;
        if (mode == 0) { h = DrawHalfBits(io); }
        else if (mode == 1)
        {
            h = static_cast<uint16_t>(IEEEBinary16(static_cast<float>(io.raw(-4.0, 4.0))));
        }
        else
        {
            h = static_cast<uint16_t>(IEEEBinary16(static_cast<float>(io.rawInteger(-16, 16) / 8.0)));
        }
        return GivenHalf(io, h);
    }

    void OutHalf(oracle::Ctx& io, IEEEBinary16 const& h)
    {
        io.outReal(static_cast<double>(h));
    }
}

ORACLE_CASE("IEEEBinary16.mathExact")
{
    // The wrappers whose binary32 function is correctly rounded or exact:
    // ceil, floor, fabs, sqrt, fmod, frexp, ldexp, and the Functions.h
    // wrappers clamp, isign, saturate, sign, sqr, invsqrt (1/sqrt in
    // binary32: two roundings). Results are emitted as values (the sign of
    // a NaN is x64-specific, see IEEEBinary16.arithmetic).
    IEEEBinary16 x(DrawMathHalf(io));
    IEEEBinary16 y(DrawMathHalf(io));
    IEEEBinary16 z(DrawMathHalf(io));
    int mode = io.rawInteger(0, 3);
    int32_t e = io.rawInteger(-40, 40);
    if (mode == 1) { e = io.rawInteger(1070, 1200) * (io.rawInteger(0, 1) != 0 ? -1 : 1); }
    int32_t exponent = static_cast<int32_t>(io.given(static_cast<double>(e)));
    OutHalf(io, std::ceil(x));
    OutHalf(io, std::floor(x));
    OutHalf(io, std::fabs(x));
    OutHalf(io, std::sqrt(x));
    OutHalf(io, std::fmod(x, y));
    int32_t fexp = 0;
    IEEEBinary16 fr = std::frexp(x, &fexp);
    OutHalf(io, fr);
    if (x.IsFinite()) { io.outInt(fexp); }
    OutHalf(io, std::ldexp(x, exponent));
    IEEEBinary16 lo = (y < z ? y : z), hi = (y < z ? z : y);
    OutHalf(io, clamp(x, lo, hi));
    io.outInt(isign(x));
    OutHalf(io, saturate(x));
    OutHalf(io, sign(x));
    OutHalf(io, sqr(x));
    OutHalf(io, invsqrt(x));
}

ORACLE_CASE("IEEEBinary16.mathLibm")
{
    // The wrappers around binary32 libm functions (MSVC acosf, ..., tanhf;
    // Functions.h atandivpi, atan2divpi, cospi, sinpi, exp10 in float). The
    // port evaluates the binary64 V8 function and rounds to binary32, so a
    // result can differ by one binary32 ulp, which moves the binary16 result
    // by at most one binary16 ulp (2^-10 relative): tolerance 1e-3.
    IEEEBinary16 x(DrawMathHalf(io));
    IEEEBinary16 y(DrawMathHalf(io));
    IEEEBinary16 r[24] = {
        std::acos(x), std::acosh(x), std::asin(x), std::asinh(x), std::atan(x),
        std::atanh(x), std::atan2(x, y), std::cos(x), std::cosh(x), std::exp(x),
        std::exp2(x), std::log(x), std::log2(x), std::log10(x), std::pow(x, y),
        std::sin(x), std::sinh(x), std::tan(x), std::tanh(x),
        atandivpi(x), atan2divpi(x, y), cospi(x), sinpi(x), exp10(x) };
    for (auto const& v : r) { OutHalf(io, v); }
}

// ================================================================ PolygonWindingOrder

ORACLE_CASE("PolygonWindingOrder.operator")
{
    // Star-shaped polygons (vertices sorted by angle around the origin,
    // then reversed half the time) on a small lattice or with real
    // coordinates, collinear runs through the lexicographically smallest
    // vertex, and fully degenerate (collinear) polygons.
    int mode = io.integer(0, 3);
    int n = io.integer(3, 8);
    std::vector<Vector2<double>> polygon(static_cast<size_t>(n));
    if (mode == 3)
    {
        // Collinear: every vertex on one lattice line.
        auto base = io.latticeVec<2>(-3, 3);
        auto dir = io.latticeDir<2>(-2, 2);
        for (auto& p : polygon)
        {
            double k = io.lattice(-3, 3);
            p = base + k * dir;
        }
    }
    else
    {
        struct AngleVertex { double angle; Vector2<double> p; };
        std::vector<AngleVertex> raw(static_cast<size_t>(n));
        for (auto& v : raw)
        {
            if (mode == 0) { v.p = { (double)io.rawInteger(-4, 4), (double)io.rawInteger(-4, 4) }; }
            else if (mode == 1) { v.p = { io.raw(-5.0, 5.0), io.raw(-5.0, 5.0) }; }
            else
            {
                // A collinear run through the lower-left corner (-4,-4).
                int k = io.rawInteger(0, 4);
                v.p = (io.rawInteger(0, 1) != 0 ? Vector2<double>{ -4.0 + k, -4.0 }
                    : Vector2<double>{ (double)io.rawInteger(-4, 4), (double)io.rawInteger(-3, 4) });
            }
            v.angle = std::atan2(v.p[1], v.p[0]);
        }
        std::sort(raw.begin(), raw.end(), [](AngleVertex const& a, AngleVertex const& b)
            { return a.angle < b.angle || (a.angle == b.angle && a.p < b.p); });
        bool reverse = (io.rawInteger(0, 1) != 0);
        for (int i = 0; i < n; ++i)
        {
            polygon[static_cast<size_t>(i)] = raw[static_cast<size_t>(reverse ? n - 1 - i : i)].p;
        }
        for (auto& p : polygon) { io.givenVec(p); }
    }
    PolygonWindingOrder<double> query;
    io.outBool(query(polygon));
}

// ================================================================ PolygonTree

namespace
{
    // A PolygonTreeEx with 1..6 nodes in breadth-first order, each node with
    // a chirality in {-1, 0, +1} and 0..3 triangles over a pool of lattice
    // points (degenerate triangles included), plus the triangle lists the
    // triangulators fill.
    struct TreeDraw
    {
        std::vector<Vector2<double>> points;
        PolygonTreeEx tree;
    };

    TreeDraw DrawTree(oracle::Ctx& io)
    {
        TreeDraw d{};
        int numPoints = io.integer(3, 8);
        d.points.resize(static_cast<size_t>(numPoints));
        for (auto& p : d.points) { p = io.latticeVec<2>(-4, 4); }
        int numNodes = io.integer(1, 6);
        d.tree.nodes.resize(static_cast<size_t>(numNodes));
        // Breadth-first: node k's parent is drawn among the earlier nodes in
        // nondecreasing order, so children of a node are contiguous.
        std::vector<size_t> parent(static_cast<size_t>(numNodes), 0);
        for (int k = 1; k < numNodes; ++k)
        {
            int lo = static_cast<int>(parent[static_cast<size_t>(k - 1)]);
            parent[static_cast<size_t>(k)] = static_cast<size_t>(io.integer(k == 1 ? 0 : lo, k - 1));
        }
        for (int k = 0; k < numNodes; ++k)
        {
            auto& node = d.tree.nodes[static_cast<size_t>(k)];
            node.self = static_cast<size_t>(k);
            node.parent = (k == 0 ? std::numeric_limits<size_t>::max() : parent[static_cast<size_t>(k)]);
            node.minChild = static_cast<size_t>(numNodes);
            node.supChild = static_cast<size_t>(numNodes);
            node.chirality = io.integer(-1, 1);
            int numTriangles = io.integer(0, 3);
            node.triangulation.resize(static_cast<size_t>(numTriangles));
            for (auto& t : node.triangulation)
            {
                for (auto& v : t) { v = io.integer(0, numPoints - 1); }
            }
        }
        for (int k = numNodes - 1; k >= 1; --k)
        {
            auto& p = d.tree.nodes[parent[static_cast<size_t>(k)]];
            p.minChild = static_cast<size_t>(k);
            if (p.supChild == static_cast<size_t>(numNodes) && k + 1 <= numNodes)
            {
                p.supChild = static_cast<size_t>(k + 1);
            }
        }
        for (auto& node : d.tree.nodes)
        {
            if (node.minChild == static_cast<size_t>(numNodes)) { node.supChild = node.minChild; }
        }
        for (size_t k = 0; k < d.tree.nodes.size(); ++k)
        {
            for (auto const& t : d.tree.nodes[k].triangulation)
            {
                d.tree.insideTriangles.push_back(t);
                d.tree.insideNodeIndices.push_back(k);
            }
        }
        return d;
    }

    Vector2<double> DrawTest(oracle::Ctx& io, TreeDraw const& d)
    {
        // A pool point, the midpoint of two pool points (on an edge), or a
        // lattice point.
        int mode = io.rawInteger(0, 2);
        Vector2<double> test{};
        if (mode == 0) { test = d.points[static_cast<size_t>(io.rawInteger(0, static_cast<int>(d.points.size()) - 1))]; }
        else if (mode == 1)
        {
            auto const& a = d.points[static_cast<size_t>(io.rawInteger(0, static_cast<int>(d.points.size()) - 1))];
            auto const& b = d.points[static_cast<size_t>(io.rawInteger(0, static_cast<int>(d.points.size()) - 1))];
            test = 0.5 * (a + b);
        }
        else { test = { (double)io.rawInteger(-4, 4), (double)io.rawInteger(-4, 4) }; }
        return io.givenVec(test);
    }

    void OutIndexOrInvalid(oracle::Ctx& io, size_t i)
    {
        if (i == std::numeric_limits<size_t>::max()) { io.outInt(-1); }
        else { io.outInt(i); }
    }
}

ORACLE_CASE("PolygonTree.getContainingTriangle")
{
    // The three point-containment queries of PolygonTreeEx on lattice
    // triangles with test points at vertices, on edges and elsewhere, so the
    // 'sdot > 0' tests are evaluated at exact zero. The tree search visits
    // the nodes in stack order; the list searches return the first hit.
    TreeDraw d = DrawTree(io);
    Vector2<double> test = DrawTest(io, d);
    int64_t chirality = io.integer(-1, 1);
    auto r0 = d.tree.GetContainingTriangle(test, d.points.data());
    auto r1 = d.tree.GetContainingTriangle(test, d.tree.insideTriangles,
        d.tree.insideNodeIndices, d.points.data());
    size_t r2 = d.tree.GetContainingTriangle(test, d.tree.insideTriangles,
        chirality, d.points.data());
    OutIndexOrInvalid(io, r0.first);
    OutIndexOrInvalid(io, r0.second);
    OutIndexOrInvalid(io, r1.first);
    OutIndexOrInvalid(io, r1.second);
    OutIndexOrInvalid(io, r2);
}

ORACLE_CASE("PolygonTree.getContainingTriangle.invalidArgument")
{
    // LogAssert(triangles.size() == nodeIndices.size()): throw parity.
    TreeDraw d = DrawTree(io);
    Vector2<double> test = DrawTest(io, d);
    int drop = io.integer(0, 1);
    std::vector<size_t> nodeIndices = d.tree.insideNodeIndices;
    if (drop == 1) { nodeIndices.push_back(0); }
    auto r = d.tree.GetContainingTriangle(test, d.tree.insideTriangles,
        nodeIndices, d.points.data());
    OutIndexOrInvalid(io, r.first);
    OutIndexOrInvalid(io, r.second);
}

// ================================================================ surface meshes

namespace
{
    struct SurfaceMesh
    {
        int kind;
        std::vector<Vector3<double>> vertices;
        std::vector<int32_t> indices;
    };

    // Small triangle meshes: 0 tetrahedron, 1 octahedron, 2 cube (closed),
    // 3 a heightfield grid with random diagonals (open), 4 an octahedron
    // subdivided once and projected to a sphere of radius R (closed; the
    // replay checks curvature estimates near 1/R). The closed solids are
    // scaled by an integer and either kept on the lattice (exact normal
    // sums, umbilic vertices of the regular octahedron) or perturbed.
    SurfaceMesh DrawSurfaceMeshRaw(oracle::Ctx& io, int forcedKind)
    {
        SurfaceMesh m{};
        m.kind = (forcedKind >= 0 ? forcedKind : io.rawInteger(0, 4));
        bool perturb = (io.rawInteger(0, 1) != 0);
        double scale = static_cast<double>(io.rawInteger(1, 3));
        auto jitter = [&io, perturb](Vector3<double> const& p)
        {
            if (!perturb) { return p; }
            Vector3<double> q = p;
            for (int j = 0; j < 3; ++j) { q[j] += io.raw(-0.2, 0.2); }
            return q;
        };
        if (m.kind == 0)
        {
            Vector3<double> v[4] = { {1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1} };
            for (auto const& p : v) { m.vertices.push_back(jitter(scale * p)); }
            m.indices = { 0,1,2, 0,3,1, 0,2,3, 1,3,2 };
        }
        else if (m.kind == 1 || m.kind == 4)
        {
            Vector3<double> v[6] = { {1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1} };
            for (auto const& p : v) { m.vertices.push_back(p); }
            std::vector<int32_t> tri = { 0,2,4, 2,1,4, 1,3,4, 3,0,4, 2,0,5, 1,2,5, 3,1,5, 0,3,5 };
            if (m.kind == 1) { m.indices = tri; }
            else
            {
                std::map<std::pair<int32_t, int32_t>, int32_t> mid;
                auto midpoint = [&m, &mid](int32_t a, int32_t b)
                {
                    auto key = std::make_pair(std::min(a, b), std::max(a, b));
                    auto it = mid.find(key);
                    if (it != mid.end()) { return it->second; }
                    int32_t k = static_cast<int32_t>(m.vertices.size());
                    m.vertices.push_back(0.5 * (m.vertices[a] + m.vertices[b]));
                    mid[key] = k;
                    return k;
                };
                for (size_t t = 0; t < tri.size(); t += 3)
                {
                    int32_t a = tri[t], b = tri[t + 1], c = tri[t + 2];
                    int32_t ab = midpoint(a, b), bc = midpoint(b, c), ca = midpoint(c, a);
                    int32_t sub[12] = { a,ab,ca, ab,b,bc, ca,bc,c, ab,bc,ca };
                    m.indices.insert(m.indices.end(), sub, sub + 12);
                }
            }
            double radius = scale * (m.kind == 4 ? io.raw(0.5, 2.0) : 1.0);
            for (auto& p : m.vertices)
            {
                if (m.kind == 4) { Normalize(p); }
                p = jitter(radius * p);
            }
            if (m.kind == 4 && perturb)
            {
                // Keep the sphere samples on the sphere (radial jitter only).
                for (auto& p : m.vertices) { Normalize(p); p = radius * p; }
            }
        }
        else if (m.kind == 2)
        {
            for (int i = 0; i < 8; ++i)
            {
                Vector3<double> p{ (double)((i & 1) * 2 - 1), (double)(((i >> 1) & 1) * 2 - 1),
                    (double)(((i >> 2) & 1) * 2 - 1) };
                m.vertices.push_back(jitter(scale * p));
            }
            m.indices = { 0,2,1, 1,2,3, 4,5,6, 5,7,6, 0,1,4, 1,5,4,
                2,6,3, 3,6,7, 0,4,2, 2,4,6, 1,3,5, 3,7,5 };
        }
        else
        {
            int rows = io.rawInteger(2, 4), cols = io.rawInteger(2, 4);
            bool latticeZ = (io.rawInteger(0, 1) != 0);
            for (int r = 0; r < rows; ++r)
            {
                for (int c = 0; c < cols; ++c)
                {
                    double z = (latticeZ ? (double)io.rawInteger(-1, 1) : io.raw(-1.0, 1.0));
                    m.vertices.push_back({ scale * c, scale * r, z });
                }
            }
            for (int r = 0; r + 1 < rows; ++r)
            {
                for (int c = 0; c + 1 < cols; ++c)
                {
                    int32_t v00 = r * cols + c, v10 = v00 + 1, v01 = v00 + cols, v11 = v01 + 1;
                    if (io.rawInteger(0, 1) != 0)
                    {
                        int32_t t[6] = { v00, v10, v11, v00, v11, v01 };
                        m.indices.insert(m.indices.end(), t, t + 6);
                    }
                    else
                    {
                        int32_t t[6] = { v00, v10, v01, v10, v11, v01 };
                        m.indices.insert(m.indices.end(), t, t + 6);
                    }
                }
            }
        }
        return m;
    }

    SurfaceMesh GivenSurfaceMesh(oracle::Ctx& io, SurfaceMesh const& m)
    {
        io.given(static_cast<double>(m.kind));
        io.given(static_cast<double>(m.vertices.size()));
        for (auto const& p : m.vertices) { io.givenVec(p); }
        io.given(static_cast<double>(m.indices.size()));
        for (int32_t i : m.indices) { io.given(static_cast<double>(i)); }
        return m;
    }

    // Arithmetic-only weights, defined identically in the replay, so that
    // the virtual hooks of MeshSmoother are compared too.
    class WeightedSmoother : public MeshSmoother<double>
    {
    protected:
        virtual bool VertexInfluenced(size_t i, double t) override
        {
            return (i % 3 != 1) || t > 0.75;
        }
        virtual double GetTangentWeight(size_t, double t) override
        {
            return 0.25 + 0.125 * t;
        }
        virtual double GetNormalWeight(size_t i, double t) override
        {
            return (i % 2 == 0 ? -0.0625 : 0.03125) * t;
        }
    };
}

ORACLE_CASE("MeshSmoother.update")
{
    // Three Update(t) steps (t = 0, 0.5, 1) with the default weights or the
    // WeightedSmoother hooks. Some records append a vertex that no triangle
    // references (neighbour count 0: Vector's operator/= sets the mean to
    // zero). After each step: normals, means and the vertex positions.
    SurfaceMesh m = DrawSurfaceMeshRaw(io, -1);
    if (io.rawInteger(0, 3) == 0) { m.vertices.push_back({ io.raw(-1.0, 1.0), 1.5, -0.5 }); }
    GivenSurfaceMesh(io, m);
    bool weighted = io.boolean();
    WeightedSmoother weightedSmoother;
    MeshSmoother<double> plainSmoother;
    MeshSmoother<double>& smoother = (weighted ? static_cast<MeshSmoother<double>&>(weightedSmoother) : plainSmoother);
    std::vector<Vector3<double>> vertices = m.vertices;
    smoother(vertices, m.indices);
    io.outInt(smoother.GetNumVertices());
    io.outInt(smoother.GetNumTriangles());
    for (size_t c : smoother.GetNeighborCounts()) { io.outInt(c); }
    for (int step = 0; step < 3; ++step)
    {
        smoother.Update(0.5 * step);
        for (auto const& v : smoother.GetNormals()) { io.outVec(v); }
        for (auto const& v : smoother.GetMeans()) { io.outVec(v); }
        for (auto const& v : vertices) { io.outVec(v); }
    }
}

ORACLE_CASE("MeshSmoother.invalidInput")
{
    // LogAssert(numVertices >= 3 && numTriangles >= 1): throw parity.
    int numVertices = io.integer(1, 4);
    int numIndices = io.integer(0, 4);
    std::vector<Vector3<double>> vertices(static_cast<size_t>(numVertices), Vector3<double>{ 0.0, 0.0, 0.0 });
    std::vector<int32_t> indices(static_cast<size_t>(numIndices), 0);
    MeshSmoother<double> smoother;
    smoother(vertices, indices);
    io.outInt(smoother.GetNumTriangles());
}

ORACLE_CASE("MeshCurvature.compute")
{
    // Normals, principal curvatures and directions (closed-form 2x2
    // eigensystem, sqrt only). The singularity threshold is 0 (the
    // documented default, which never fires: finding #240), 1e-3, or 1e10
    // (every vertex takes the planar branch). The lattice octahedron has
    // exact umbilics (finding #412: zero directions).
    SurfaceMesh m = GivenSurfaceMesh(io, DrawSurfaceMeshRaw(io, -1));
    int thresholdMode = io.rawInteger(0, 2);
    double threshold = io.given(thresholdMode == 0 ? 0.0 : (thresholdMode == 1 ? 1e-3 : 1e10));
    std::vector<uint32_t> indices(m.indices.begin(), m.indices.end());
    MeshCurvature<double> curvature;
    curvature(m.vertices, indices, threshold);
    for (auto const& v : curvature.GetNormals()) { io.outVec(v); }
    for (double k : curvature.GetMinCurvatures()) { io.outReal(k); }
    for (double k : curvature.GetMaxCurvatures()) { io.outReal(k); }
    for (auto const& v : curvature.GetMinDirections()) { io.outVec(v); }
    for (auto const& v : curvature.GetMaxDirections()) { io.outVec(v); }
}

// ================================================================ ImplicitSurface3

namespace
{
    // F = c0 x^2 + c1 y^2 + c2 z^2 + c3 xy + c4 yz + c5 xz + c6 x + c7 y
    //   + c8 z + c9 + c10 xyz, with every sum grouped left to right as
    // written; the replay defines the same class with the same grouping
    // (ORACLE.md, the v18 precedent), so the base-class algorithms are
    // compared bit for bit.
    class CubicSurface : public ImplicitSurface3<double>
    {
    public:
        CubicSurface(std::array<double, 11> const& c) : mC(c) {}

        virtual double F(Vector3<double> const& p) const override
        {
            double x = p[0], y = p[1], z = p[2];
            return ((((((((((mC[0] * (x * x) + mC[1] * (y * y)) + mC[2] * (z * z))
                + mC[3] * (x * y)) + mC[4] * (y * z)) + mC[5] * (x * z)) + mC[6] * x)
                + mC[7] * y) + mC[8] * z) + mC[9]) + mC[10] * ((x * y) * z));
        }
        virtual double FX(Vector3<double> const& p) const override
        {
            return (((((2.0 * mC[0]) * p[0] + mC[3] * p[1]) + mC[5] * p[2]) + mC[6]) + mC[10] * (p[1] * p[2]));
        }
        virtual double FY(Vector3<double> const& p) const override
        {
            return (((((2.0 * mC[1]) * p[1] + mC[3] * p[0]) + mC[4] * p[2]) + mC[7]) + mC[10] * (p[0] * p[2]));
        }
        virtual double FZ(Vector3<double> const& p) const override
        {
            return (((((2.0 * mC[2]) * p[2] + mC[4] * p[1]) + mC[5] * p[0]) + mC[8]) + mC[10] * (p[0] * p[1]));
        }
        virtual double FXX(Vector3<double> const&) const override { return 2.0 * mC[0]; }
        virtual double FXY(Vector3<double> const& p) const override { return mC[3] + mC[10] * p[2]; }
        virtual double FXZ(Vector3<double> const& p) const override { return mC[5] + mC[10] * p[1]; }
        virtual double FYY(Vector3<double> const&) const override { return 2.0 * mC[1]; }
        virtual double FYZ(Vector3<double> const& p) const override { return mC[4] + mC[10] * p[0]; }
        virtual double FZZ(Vector3<double> const&) const override { return 2.0 * mC[2]; }

    private:
        std::array<double, 11> mC;
    };
}

ORACLE_CASE("ImplicitSurface3.queries")
{
    // Modes: 0 random coefficients and position; 1 the sphere
    // x^2 + y^2 + z^2 - R^2 at a point on it (lattice Pythagorean points or
    // a scaled unit vector; the replay checks both principal curvatures are
    // 1/R); 2 lattice coefficients and lattice position (exact gradients,
    // repeated eigenvalues, axis-aligned normals); 3 a zero gradient (the
    // sphere at its center, or a lattice quadric at its critical point).
    int mode = io.integer(0, 3);
    std::array<double, 11> c{};
    Vector3<double> p{};
    if (mode == 0)
    {
        for (auto& x : c) { x = io.real(-2.0, 2.0); }
        p = io.vec<3>(-2.0, 2.0);
    }
    else if (mode == 1)
    {
        static int const pyth[4][4] = { {1,2,2,3}, {2,3,6,7}, {1,4,8,9}, {2,6,9,11} };
        int k = io.rawInteger(0, 4);
        Vector3<double> q{};
        double radius = 0.0;
        if (k < 4)
        {
            radius = static_cast<double>(pyth[k][3]);
            for (int j = 0; j < 3; ++j)
            {
                q[j] = (io.rawInteger(0, 1) != 0 ? -1.0 : 1.0) * pyth[k][j];
            }
        }
        else
        {
            radius = io.raw(0.5, 3.0);
            Vector3<double> u{ io.raw(-1.0, 1.0), io.raw(-1.0, 1.0), io.raw(-1.0, 1.0) };
            Normalize(u);
            q = radius * u;
        }
        c = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -radius * radius, 0.0 };
        for (auto& x : c) { io.given(x); }
        p = io.givenVec(q);
    }
    else if (mode == 2)
    {
        for (auto& x : c) { x = io.lattice(-2, 2); }
        p = io.latticeVec<3>(-2, 2);
    }
    else
    {
        c = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0 };
        if (io.rawInteger(0, 1) != 0) { c[3] = 1.0; c[0] = 2.0; }
        for (auto& x : c) { io.given(x); }
        p = io.givenVec(Vector3<double>{ 0.0, 0.0, 0.0 });
    }
    double epsilon = io.real(0.0, 1.0);
    CubicSurface surface(c);
    io.outBool(surface.IsOnSurface(p, epsilon));
    io.outReal(surface.F(p));
    io.outVec(surface.GetGradient(p));
    io.outMat(surface.GetHessian(p));
    Vector3<double> t0{}, t1{}, n{};
    surface.GetFrame(p, t0, t1, n);
    io.outVec(t0);
    io.outVec(t1);
    io.outVec(n);
    double k0 = 0.0, k1 = 0.0;
    Vector3<double> d0{}, d1{};
    bool valid = surface.GetPrincipalInformation(p, k0, k1, d0, d1);
    io.outBool(valid);
    io.outReal(k0);
    io.outReal(k1);
    io.outVec(d0);
    io.outVec(d1);
}

// ================================================================ ET/VETNonmanifoldMesh
//
// Every container these classes iterate is ordered by value (std::map keyed
// by EdgeKey/TriangleKey; std::set with WeakPtrLT/SharedPtrLT, which
// compare the pointees' keys; std::set<int32_t>), so every output is
// emitted in upstream's own order and nothing is re-sorted.

namespace
{
    struct MeshOp { int type; int32_t v[3]; };

    void EmitET(oracle::Ctx& io, ETNonmanifoldMesh const& mesh)
    {
        io.outInt(mesh.GetEdges().size());
        for (auto const& e : mesh.GetEdges())
        {
            io.outInt(e.first.V[0]); io.outInt(e.first.V[1]);
            io.outInt(e.second->V[0]); io.outInt(e.second->V[1]);
            io.outInt(e.second->T.size());
            for (auto const& tw : e.second->T)
            {
                auto t = tw.lock();
                for (int j = 0; j < 3; ++j) { io.outInt(t->V[j]); }
            }
        }
        io.outInt(mesh.GetTriangles().size());
        for (auto const& t : mesh.GetTriangles())
        {
            for (int j = 0; j < 3; ++j) { io.outInt(t.first.V[j]); }
            for (int j = 0; j < 3; ++j) { io.outInt(t.second->V[j]); }
            for (int j = 0; j < 3; ++j)
            {
                auto e = t.second->E[j].lock();
                io.outInt(e->V[0]); io.outInt(e->V[1]);
            }
        }
        io.outBool(mesh.IsManifold());
        io.outBool(mesh.IsClosed());
        std::vector<std::vector<std::shared_ptr<ETNonmanifoldMesh::Triangle>>> components;
        mesh.GetComponents(components);
        io.outInt(components.size());
        for (auto const& c : components)
        {
            io.outInt(c.size());
            for (auto const& t : c) { for (int j = 0; j < 3; ++j) { io.outInt(t->V[j]); } }
        }
        std::vector<std::vector<TriangleKey<true>>> keyComponents;
        mesh.GetComponents(keyComponents);
        io.outInt(keyComponents.size());
        for (auto const& c : keyComponents)
        {
            io.outInt(c.size());
            for (auto const& k : c) { for (int j = 0; j < 3; ++j) { io.outInt(k.V[j]); } }
        }
    }

    void EmitVET(oracle::Ctx& io, VETNonmanifoldMesh const& mesh)
    {
        EmitET(io, mesh);
        io.outInt(mesh.GetVertices().size());
        for (auto const& v : mesh.GetVertices())
        {
            io.outInt(v.first);
            io.outInt(v.second->V);
            io.outInt(v.second->VAdjacent.size());
            for (int32_t a : v.second->VAdjacent) { io.outInt(a); }
            io.outInt(v.second->EAdjacent.size());
            for (auto const& e : v.second->EAdjacent) { io.outInt(e->V[0]); io.outInt(e->V[1]); }
            io.outInt(v.second->TAdjacent.size());
            for (auto const& t : v.second->TAdjacent)
            {
                for (int j = 0; j < 3; ++j) { io.outInt(t->V[j]); }
            }
        }
    }

    // A random triangle over the vertex pool 0..poolSize-1 (distinct
    // vertices unless 'degenerate'), or, for removals, usually an existing
    // triangle in one of its three rotations.
    MeshOp DrawMeshOp(oracle::Ctx& io, ETNonmanifoldMesh const& mesh, int poolSize,
        bool allowDegenerate)
    {
        MeshOp op{};
        op.type = (io.rawInteger(0, 9) < 7 ? 0 : 1);
        auto const& tmap = mesh.GetTriangles();
        if (op.type == 1 && !tmap.empty() && io.rawInteger(0, 3) != 0)
        {
            auto it = tmap.begin();
            std::advance(it, io.rawInteger(0, static_cast<int>(tmap.size()) - 1));
            int r = io.rawInteger(0, 2);
            for (int j = 0; j < 3; ++j) { op.v[j] = it->second->V[(j + r) % 3]; }
            return op;
        }
        bool degenerate = allowDegenerate && io.rawInteger(0, 5) == 0;
        op.v[0] = io.rawInteger(0, poolSize - 1);
        do { op.v[1] = io.rawInteger(0, poolSize - 1); } while (!degenerate && op.v[1] == op.v[0]);
        if (degenerate) { op.v[2] = op.v[io.rawInteger(0, 1)]; }
        else
        {
            do { op.v[2] = io.rawInteger(0, poolSize - 1); } while (op.v[2] == op.v[0] || op.v[2] == op.v[1]);
        }
        return op;
    }

    bool IsDegenerate(int32_t const* v)
    {
        return v[0] == v[1] || v[1] == v[2] || v[2] == v[0];
    }

    // True when removing <v0,v1,v2> leaves some vertex of the triangle with
    // no triangle: VETNonmanifoldMesh::Remove's inverted LogAssert fires
    // exactly then (finding, #240).
    bool IsolatesVertex(VETNonmanifoldMesh const& mesh, int32_t const* v)
    {
        auto const& tmap = mesh.GetTriangles();
        if (tmap.find(TriangleKey<true>(v[0], v[1], v[2])) == tmap.end()) { return false; }
        for (int j = 0; j < 3; ++j)
        {
            auto it = mesh.GetVertices().find(v[j]);
            if (it != mesh.GetVertices().end() && it->second->TAdjacent.size() == 1) { return true; }
        }
        return false;
    }

    void GivenOp(oracle::Ctx& io, MeshOp const& op)
    {
        io.given(static_cast<double>(op.type));
        for (int j = 0; j < 3; ++j) { io.given(static_cast<double>(op.v[j])); }
    }
}

ORACLE_CASE("ETNonmanifoldMesh.sequence")
{
    // Up to 16 inserts and removes over a pool of 4..6 vertices: shared,
    // nonmanifold (three or more triangles on an edge), duplicate (the same
    // key in another rotation: Insert returns null), reversed (a different
    // key) and degenerate triangles (repeated vertex, #179: accepted by
    // Insert, one Edge object for two keys). Removing a degenerate triangle
    // throws after corrupting the mesh (#179, preserved), so such removals
    // are turned into inserts here (the throw has its own case). After
    // every operation its result; at the end the whole mesh, a copy
    // (operator=, which reinserts by key: rotated vertices) and Clear.
    int poolSize = io.integer(4, 6);
    int numOps = io.integer(1, 16);
    ETNonmanifoldMesh mesh;
    for (int k = 0; k < numOps; ++k)
    {
        MeshOp op = DrawMeshOp(io, mesh, poolSize, true);
        if (op.type == 1 && IsDegenerate(op.v)) { op.type = 0; }
        GivenOp(io, op);
        if (op.type == 0)
        {
            auto tri = mesh.Insert(op.v[0], op.v[1], op.v[2]);
            io.outBool(tri != nullptr);
            if (tri) { for (int j = 0; j < 3; ++j) { io.outInt(tri->V[j]); } }
        }
        else
        {
            io.outBool(mesh.Remove(op.v[0], op.v[1], op.v[2]));
        }
    }
    EmitET(io, mesh);
    ETNonmanifoldMesh copy(mesh);
    EmitET(io, copy);
    copy.Clear();
    io.outInt(copy.GetEdges().size());
    io.outInt(copy.GetTriangles().size());
}

ORACLE_CASE("ETNonmanifoldMesh.remove.degenerate")
{
    // #179 (preserved): Remove of a triangle with a repeated vertex throws
    // ("Unexpected condition.") on both sides.
    int poolSize = io.integer(4, 6);
    int numOps = io.integer(0, 6);
    ETNonmanifoldMesh mesh;
    for (int k = 0; k < numOps; ++k)
    {
        MeshOp op = DrawMeshOp(io, mesh, poolSize, false);
        op.type = 0;
        GivenOp(io, op);
        mesh.Insert(op.v[0], op.v[1], op.v[2]);
    }
    int32_t a = io.integer(0, poolSize - 1);
    int32_t b = io.integer(0, poolSize - 1);
    int which = io.integer(0, 2);
    int32_t v[3] = { a, a, b };
    if (which == 1) { v[0] = b; v[1] = a; v[2] = a; }
    if (which == 2) { v[0] = a; v[1] = b; v[2] = a; }
    mesh.Insert(v[0], v[1], v[2]);
    bool removed = mesh.Remove(v[0], v[1], v[2]);
    io.outBool(removed);
}

ORACLE_CASE("VETNonmanifoldMesh.sequence")
{
    // As ETNonmanifoldMesh.sequence, plus the vertex adjacency sets. A
    // removal that would leave a vertex of the triangle with no triangle
    // is turned into an insert: upstream's inverted assertion throws there
    // (#240, fixed in the port; see the deviation case). Degenerate
    // triangles are inserted but never removed.
    int poolSize = io.integer(4, 6);
    int numOps = io.integer(1, 20);
    VETNonmanifoldMesh mesh;
    for (int k = 0; k < numOps; ++k)
    {
        MeshOp op = DrawMeshOp(io, mesh, poolSize, true);
        if (op.type == 1 && (IsDegenerate(op.v) || IsolatesVertex(mesh, op.v))) { op.type = 0; }
        GivenOp(io, op);
        if (op.type == 0)
        {
            auto tri = mesh.Insert(op.v[0], op.v[1], op.v[2]);
            io.outBool(tri != nullptr);
        }
        else
        {
            io.outBool(mesh.Remove(op.v[0], op.v[1], op.v[2]));
        }
    }
    EmitVET(io, mesh);
    VETNonmanifoldMesh copy(mesh);
    EmitVET(io, copy);
    copy.Clear();
    io.outInt(copy.GetVertices().size());
}

ORACLE_CASE("VETNonmanifoldMesh.remove.isolatesVertex")
{
    // Deliberate port fix (#240): upstream's Remove asserts
    // 'VAdjacent.size() != 0 || EAdjacent.size() != 0' when a vertex loses
    // its last triangle, which is inverted (both sets have just been
    // emptied), so every such valid removal throws upstream; the port
    // removes the triangle. Some inserts, then the removal of an existing
    // nondegenerate triangle that isolates a vertex.
    int poolSize = io.integer(4, 6);
    int numOps = io.integer(0, 8);
    VETNonmanifoldMesh mesh;
    for (int k = 0; k < numOps; ++k)
    {
        MeshOp op = DrawMeshOp(io, mesh, poolSize, false);
        op.type = 0;
        GivenOp(io, op);
        mesh.Insert(op.v[0], op.v[1], op.v[2]);
    }
    std::vector<std::array<int32_t, 3>> candidates;
    for (auto const& t : mesh.GetTriangles())
    {
        if (IsolatesVertex(mesh, t.second->V.data())) { candidates.push_back(t.second->V); }
    }
    std::array<int32_t, 3> target{ 0, 1, poolSize };
    if (candidates.empty())
    {
        mesh.Insert(target[0], target[1], target[2]);
    }
    else
    {
        target = candidates[static_cast<size_t>(io.rawInteger(0, static_cast<int>(candidates.size()) - 1))];
    }
    io.given(static_cast<double>(candidates.empty() ? 1 : 0));
    for (int32_t x : target) { io.given(static_cast<double>(x)); }
    bool removed = mesh.Remove(target[0], target[1], target[2]);
    io.outBool(removed);
    EmitVET(io, mesh);
}
