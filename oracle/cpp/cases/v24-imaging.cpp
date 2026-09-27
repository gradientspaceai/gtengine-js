// Group 24 (imaging): FastGaussianBlur1, FastGaussianBlur2,
// FastGaussianBlur3, Image, MarchingCubes, PdeFilter,
// AdaptiveSkeletonClimbing2, FastMarch, Histogram, Image2.
//
// Everything here is integer bookkeeping or IEEE-exact arithmetic (+ - * /,
// floor, ceil, conversions), so every case is compared bit for bit.
//
// Upstream functions whose outputs can carry a signed zero are called
// through __declspec(noinline) wrappers (ORACLE.md, v23). Large outputs (the
// 256-entry MarchingCubes table, AdaptiveSkeletonClimbing2 on 33x33 and
// 65x65 images) are emitted partly as 32-bit digests computed identically on
// both sides; each case comment says which.
#define ORACLE_FAMILY "v24-imaging"
#include "Oracle.h"

#include <Mathematics/AdaptiveSkeletonClimbing2.h>
#include <Mathematics/BSRational.h>
#include <Mathematics/FastGaussianBlur1.h>
#include <Mathematics/FastGaussianBlur2.h>
#include <Mathematics/FastGaussianBlur3.h>
#include <Mathematics/FastMarch.h>
#include <Mathematics/Histogram.h>
#include <Mathematics/Image.h>
#include <Mathematics/Image2.h>
#include <Mathematics/MarchingCubes.h>
#include <Mathematics/PdeFilter.h>
#include <Mathematics/UIntegerAP32.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

using namespace gte;

namespace
{
    // ---- digests (mirrored in the replay) ----

    // FNV-1a over 32-bit words.
    struct Digest
    {
        uint32_t h = 2166136261u;
        void word(uint32_t w)
        {
            for (int i = 0; i < 4; ++i)
            {
                h ^= (w >> (8 * i)) & 0xFFu;
                h *= 16777619u;
            }
        }
        void integer(int64_t v) { word(static_cast<uint32_t>(v)); }
        void real(double x)
        {
            uint64_t bits;
            std::memcpy(&bits, &x, sizeof(bits));
            word(static_cast<uint32_t>(bits));
            word(static_cast<uint32_t>(bits >> 32));
        }
    };

    // A string as its length followed by 6 characters per output.
    void OutString(oracle::Ctx& io, std::string const& s)
    {
        io.outInt(s.size());
        for (size_t i = 0; i < s.size(); i += 6)
        {
            int64_t packed = 0;
            for (size_t j = i; j < std::min(s.size(), i + 6); ++j)
            {
                packed = packed * 256 + static_cast<unsigned char>(s[j]);
            }
            io.outInt(packed);
        }
    }

    // ---- FastGaussianBlur ----

    // One recorded double per pixel whatever the mode: uniform, a small
    // lattice, a constant image (the blur must reproduce it exactly) or a
    // single spike on zeros. Integer types stay in [-1000, 1000] so that
    // every output fits int16_t (static_cast of an out-of-range double is
    // undefined behaviour).
    template <typename T>
    std::vector<T> DrawPixels(oracle::Ctx& io, size_t n)
    {
        bool const isInt = std::is_integral<T>::value;
        int mode = io.rawInteger(0, 3);
        double constant = isInt ? io.rawInteger(-1000, 1000) : io.raw(-10.0, 10.0);
        size_t spike = static_cast<size_t>(io.rawInteger(0, static_cast<int>(n) - 1));
        std::vector<T> pixels(n);
        for (size_t i = 0; i < n; ++i)
        {
            double v = 0.0;
            switch (mode)
            {
            case 0: v = isInt ? io.rawInteger(-1000, 1000) : io.raw(-10.0, 10.0); break;
            case 1: v = io.rawInteger(-3, 3); break;
            case 2: v = constant; break;
            default: v = (i == spike ? constant : 0.0); break;
            }
            pixels[i] = static_cast<T>(v);
            io.given(static_cast<double>(pixels[i]));
        }
        return pixels;
    }

    // One recorded double. Integer and half-integer scales put rxps and rxms
    // on the lattice (delta = 0 exactly), scales past the image bound select
    // the boundary branches everywhere, 0 is the degenerate step.
    double DrawScale(oracle::Ctx& io)
    {
        switch (io.rawInteger(0, 5))
        {
        case 0: return io.given(io.raw(0.0, 3.0));
        case 1: return io.given(io.rawInteger(0, 5));
        case 2: return io.given(0.5 * io.rawInteger(0, 9));
        case 3: return io.given(io.raw(3.0, 30.0));
        case 4: return io.given(io.raw(0.0, 1e-3));
        default: return io.given(0.0);
        }
    }

    double DrawLogBase(oracle::Ctx& io)
    {
        switch (io.rawInteger(0, 3))
        {
        case 0: return io.given(io.raw(0.0, 0.5));
        case 1: return io.given(io.rawInteger(-8, 16) / 16.0);
        case 2: return io.given(io.raw(-1.0, 2.0));
        default: return io.given(0.0);
        }
    }

    template <typename T>
    __declspec(noinline) void CallBlur1(int32_t xb, T const* in, T* out, double s, double lb)
    {
        FastGaussianBlur1<T> blur;
        blur.Execute(xb, in, out, s, lb);
    }

    template <typename T>
    __declspec(noinline) void CallBlur2(int32_t xb, int32_t yb, T const* in, T* out,
        double s, double lb)
    {
        FastGaussianBlur2<T> blur;
        blur.Execute(xb, yb, in, out, s, lb);
    }

    template <typename T>
    __declspec(noinline) void CallBlur3(int32_t xb, int32_t yb, int32_t zb, T const* in,
        T* out, double s, double lb)
    {
        FastGaussianBlur3<T> blur;
        blur.Execute(xb, yb, zb, in, out, s, lb);
    }

    template <typename T>
    void Blur1Case(oracle::Ctx& io)
    {
        int32_t xBound = io.integer(1, 24);
        std::vector<T> input = DrawPixels<T>(io, xBound);
        double scale = DrawScale(io);
        double logBase = DrawLogBase(io);
        std::vector<T> output(xBound);
        CallBlur1<T>(xBound, input.data(), output.data(), scale, logBase);
        for (auto v : output) { io.outReal(static_cast<double>(v)); }
    }

    template <typename T>
    void Blur2Case(oracle::Ctx& io)
    {
        int32_t xBound = io.integer(1, 8);
        int32_t yBound = io.integer(1, 8);
        std::vector<T> input = DrawPixels<T>(io, static_cast<size_t>(xBound) * yBound);
        double scale = DrawScale(io);
        double logBase = DrawLogBase(io);
        std::vector<T> output(input.size());
        CallBlur2<T>(xBound, yBound, input.data(), output.data(), scale, logBase);
        for (auto v : output) { io.outReal(static_cast<double>(v)); }
    }

    template <typename T>
    void Blur3Case(oracle::Ctx& io)
    {
        int32_t xBound = io.integer(1, 5);
        int32_t yBound = io.integer(1, 5);
        int32_t zBound = io.integer(1, 5);
        std::vector<T> input = DrawPixels<T>(io,
            static_cast<size_t>(xBound) * yBound * zBound);
        double scale = DrawScale(io);
        double logBase = DrawLogBase(io);
        std::vector<T> output(input.size());
        CallBlur3<T>(xBound, yBound, zBound, input.data(), output.data(), scale, logBase);
        for (auto v : output) { io.outReal(static_cast<double>(v)); }
    }
}

ORACLE_CASE("FastGaussianBlur1.execute.double") { Blur1Case<double>(io); }
ORACLE_CASE("FastGaussianBlur1.execute.float") { Blur1Case<float>(io); }
ORACLE_CASE("FastGaussianBlur1.execute.int32") { Blur1Case<int32_t>(io); }
ORACLE_CASE("FastGaussianBlur1.execute.int16") { Blur1Case<int16_t>(io); }
ORACLE_CASE("FastGaussianBlur2.execute.double") { Blur2Case<double>(io); }
ORACLE_CASE("FastGaussianBlur2.execute.float") { Blur2Case<float>(io); }
ORACLE_CASE("FastGaussianBlur2.execute.int32") { Blur2Case<int32_t>(io); }
ORACLE_CASE("FastGaussianBlur2.execute.int16") { Blur2Case<int16_t>(io); }
ORACLE_CASE("FastGaussianBlur3.execute.double") { Blur3Case<double>(io); }
ORACLE_CASE("FastGaussianBlur3.execute.float") { Blur3Case<float>(io); }
ORACLE_CASE("FastGaussianBlur3.execute.int32") { Blur3Case<int32_t>(io); }
ORACLE_CASE("FastGaussianBlur3.execute.int16") { Blur3Case<int16_t>(io); }

// ---- Image ----

namespace
{
    // A dimension that is positive five times in six; nonpositive dimensions
    // leave the image empty.
    int32_t DrawDimension(oracle::Ctx& io, int32_t maxDim)
    {
        int32_t d = (io.rawInteger(0, 5) == 0 ? io.rawInteger(-2, 0) : io.rawInteger(1, maxDim));
        io.given(d);
        return d;
    }

    // size_t results as signed 64-bit integers: an index that wrapped below
    // zero (a negative coordinate) comes out negative, as in the port.
    int64_t Signed(size_t i) { return static_cast<int64_t>(i); }

    void OutImageShape(oracle::Ctx& io, Image<double> const& image)
    {
        io.outInt(image.GetNumDimensions());
        io.outInt(image.GetDimensions().size());
        for (int32_t d = 0; d < image.GetNumDimensions(); ++d)
        {
            io.outInt(image.GetDimension(d));
            io.outInt(image.GetOffset(d));
            io.outInt(image.GetOffsets()[d]);
        }
        io.outInt(image.GetNumPixels());
        io.outInt(image.GetPixels().size());
    }
}

// Construction (directly or through Reconstruct, over a previously
// reconstructed image), GetIndex on in-range, out-of-range and negative
// coordinates, GetCoordinates on indices up to twice the pixel count,
// operator[] writes and the clamping Get(i) on out-of-range indices
// (including size_t(-1)). The pixels are value-initialized upstream and
// zero-filled in the port; a few are read back before any write.
ORACLE_CASE("Image.access")
{
    int32_t numDims = io.integer(0, 4);
    std::vector<int32_t> dims(numDims);
    for (int32_t d = 0; d < numDims; ++d) { dims[d] = DrawDimension(io, 5); }
    int32_t how = io.integer(0, 2);
    Image<double> image;
    if (how == 0)
    {
        image = Image<double>(dims);
    }
    else
    {
        if (how == 2)
        {
            int32_t firstDims = io.integer(1, 3);
            std::vector<int32_t> first(firstDims);
            for (int32_t d = 0; d < firstDims; ++d) { first[d] = io.integer(1, 4); }
            image.Reconstruct(first);
            image.GetPixels()[0] = 7.0;
        }
        image.Reconstruct(dims);
    }
    OutImageShape(io, image);

    size_t const numPixels = image.GetNumPixels();
    int32_t const nd = image.GetNumDimensions();
    int32_t numQueries = io.integer(1, 6);
    for (int32_t q = 0; q < numQueries; ++q)
    {
        // GetIndex reads coord[0] even for an empty image.
        std::vector<int32_t> coord(std::max(nd, 1));
        for (int32_t d = 0; d < static_cast<int32_t>(coord.size()); ++d)
        {
            int32_t hi = (d < nd ? image.GetDimension(d) : 3);
            coord[d] = io.integer(-1, hi);
        }
        io.outInt(Signed(image.GetIndex(coord.data())));

        if (numPixels > 0)
        {
            int32_t index = io.integer(0, static_cast<int32_t>(2 * numPixels));
            std::vector<int32_t> c(nd);
            image.GetCoordinates(static_cast<size_t>(index), c.data());
            for (auto v : c) { io.outInt(v); }

            int32_t w = io.integer(0, static_cast<int32_t>(numPixels) - 1);
            double value = io.real(-5.0, 5.0);
            io.outReal(image[w]);
            image[w] = value;
            int32_t r = io.integer(-2, static_cast<int32_t>(numPixels) + 1);
            io.outReal(image.Get(static_cast<size_t>(static_cast<int64_t>(r))));
            io.outReal(image[static_cast<size_t>(w)]);
        }
    }
}

// ---- Image2 ----

// Construction (constructor or Reconstruct), both GetIndex and both
// GetCoordinates overloads, operator() in both forms for writes and in-range
// reads, the clamping Get in both forms on coordinates up to two outside the
// image, and the base Image::Get(i) (hidden by Image2's overloads, reached
// through a base reference; the port's one-argument getClamped).
ORACLE_CASE("Image2.access")
{
    int32_t d0 = DrawDimension(io, 7);
    int32_t d1 = DrawDimension(io, 7);
    bool viaReconstruct = io.boolean();
    Image2<double> image;
    if (viaReconstruct)
    {
        image.Reconstruct(3, 2);
        image.Reconstruct(d0, d1);
    }
    else
    {
        image = Image2<double>(d0, d1);
    }
    OutImageShape(io, image);
    size_t const numPixels = image.GetNumPixels();
    if (numPixels == 0)
    {
        return;
    }
    int32_t const dim0 = image.GetDimension(0), dim1 = image.GetDimension(1);
    for (int32_t y = 0; y < dim1; ++y)
    {
        for (int32_t x = 0; x < dim0; ++x)
        {
            if ((x + y) % 2 == 0) { image(x, y) = x + 10.0 * y + 0.5; }
            else { image(std::array<int32_t, 2>{ x, y }) = -x - 10.0 * y - 0.25; }
        }
    }
    int32_t numQueries = io.integer(1, 8);
    for (int32_t q = 0; q < numQueries; ++q)
    {
        int32_t x = io.integer(-2, dim0 + 1);
        int32_t y = io.integer(-2, dim1 + 1);
        std::array<int32_t, 2> coord{ x, y };
        io.outInt(Signed(image.GetIndex(x, y)));
        io.outInt(Signed(image.GetIndex(coord)));
        io.outReal(image.Get(x, y));
        io.outReal(image.Get(coord));
        if (0 <= x && x < dim0 && 0 <= y && y < dim1)
        {
            io.outReal(image(x, y));
            io.outReal(image(coord));
        }
        int32_t index = io.integer(0, static_cast<int32_t>(2 * numPixels));
        int32_t cx = -7, cy = -7;
        image.GetCoordinates(static_cast<size_t>(index), cx, cy);
        io.outInt(cx);
        io.outInt(cy);
        auto c = image.GetCoordinates(static_cast<size_t>(index));
        io.outInt(c[0]);
        io.outInt(c[1]);
        int32_t r = io.integer(-2, static_cast<int32_t>(numPixels) + 1);
        Image<double> const& base = image;
        io.outReal(base.Get(static_cast<size_t>(static_cast<int64_t>(r))));
    }
}

namespace
{
    template <size_t K>
    void OutRelative(oracle::Ctx& io, std::array<int32_t, K> const& nbr)
    {
        for (auto v : nbr) { io.outInt(v); }
    }

    template <size_t K>
    void OutRelative(oracle::Ctx& io, std::array<std::array<int32_t, 2>, K> const& nbr)
    {
        for (auto const& v : nbr) { io.outInt(v[0]); io.outInt(v[1]); }
    }

    template <size_t K>
    void OutAbsolute(oracle::Ctx& io, std::array<size_t, K> const& nbr, bool asUnsigned)
    {
        for (auto v : nbr)
        {
            if (asUnsigned) { io.outReal(static_cast<double>(v)); }
            else { io.outInt(Signed(v)); }
        }
    }

    template <size_t K>
    void OutAbsolute(oracle::Ctx& io, std::array<std::array<size_t, 2>, K> const& nbr,
        bool asUnsigned)
    {
        for (auto const& v : nbr)
        {
            if (asUnsigned)
            {
                io.outReal(static_cast<double>(v[0]));
                io.outReal(static_cast<double>(v[1]));
            }
            else
            {
                io.outInt(Signed(v[0]));
                io.outInt(Signed(v[1]));
            }
        }
    }

    // The relative tables, then the absolute neighbourhoods at 'count'
    // pixels. 'boundaryOnly' puts every pixel on the xmin or ymin edge, where
    // upstream's size_t results wrap (finding #64).
    void NeighborhoodCase(oracle::Ctx& io, bool boundaryOnly, bool asUnsigned)
    {
        int32_t d0 = io.integer(1, 7);
        int32_t d1 = io.integer(1, 7);
        Image2<double> image(d0, d1);
        std::array<int32_t, 4> n4{};
        std::array<int32_t, 8> n8{};
        std::array<int32_t, 4> c4{};
        std::array<int32_t, 9> f9{};
        image.GetNeighborhood(n4);
        image.GetNeighborhood(n8);
        image.GetCorners(c4);
        image.GetFull(f9);
        std::array<std::array<int32_t, 2>, 4> t4{};
        std::array<std::array<int32_t, 2>, 8> t8{};
        std::array<std::array<int32_t, 2>, 4> tc{};
        std::array<std::array<int32_t, 2>, 9> tf{};
        image.GetNeighborhood(t4);
        image.GetNeighborhood(t8);
        image.GetCorners(tc);
        image.GetFull(tf);
        if (!boundaryOnly)  // the wrap case compares only its absolute tables
        {
            OutRelative(io, n4);
            OutRelative(io, n8);
            OutRelative(io, c4);
            OutRelative(io, f9);
            OutRelative(io, t4);
            OutRelative(io, t8);
            OutRelative(io, tc);
            OutRelative(io, tf);
        }

        int32_t count = io.integer(1, boundaryOnly ? 1 : 2);
        for (int32_t k = 0; k < count; ++k)
        {
            int32_t x, y;
            if (boundaryOnly)
            {
                bool onXMin = (io.rawInteger(0, 1) == 0);
                x = static_cast<int32_t>(io.given(onXMin ? 0 : io.rawInteger(0, d0 - 1)));
                y = static_cast<int32_t>(io.given(onXMin ? io.rawInteger(0, d1 - 1) : 0));
            }
            else
            {
                x = io.integer(-1, d0);
                y = io.integer(-1, d1);
            }
            std::array<size_t, 4> a4{};
            std::array<size_t, 8> a8{};
            std::array<size_t, 4> ac{};
            std::array<size_t, 9> af{};
            image.GetNeighborhood(x, y, a4);
            image.GetNeighborhood(x, y, a8);
            image.GetCorners(x, y, ac);
            image.GetFull(x, y, af);
            OutAbsolute(io, a4, asUnsigned);
            OutAbsolute(io, a8, asUnsigned);
            OutAbsolute(io, ac, asUnsigned);
            OutAbsolute(io, af, asUnsigned);
            std::array<std::array<size_t, 2>, 4> b4{};
            std::array<std::array<size_t, 2>, 8> b8{};
            std::array<std::array<size_t, 2>, 4> bc{};
            std::array<std::array<size_t, 2>, 9> bf{};
            image.GetNeighborhood(x, y, b4);
            image.GetNeighborhood(x, y, b8);
            image.GetCorners(x, y, bc);
            image.GetFull(x, y, bf);
            OutAbsolute(io, b4, asUnsigned);
            OutAbsolute(io, b8, asUnsigned);
            OutAbsolute(io, bc, asUnsigned);
            OutAbsolute(io, bf, asUnsigned);
        }
    }
}

// Every neighbourhood table, relative and absolute, 1-D and 2-tuple, at
// pixels inside, on and one outside the border. size_t results are emitted
// as signed 64-bit integers, which is how the port represents them.
ORACLE_CASE("Image2.neighborhoods") { NeighborhoodCase(io, false, false); }

// Finding #64: at a pixel on the xmin or ymin edge upstream's size_t
// neighbour indices and coordinates wrap to SIZE_MAX (and SIZE_MAX - dim0 ...),
// the port's are negative. Emitted as the unsigned values, every record
// deviates.
ORACLE_CASE("Image2.neighborhoods.wrap") { NeighborhoodCase(io, true, true); }

// ---- Histogram ----

namespace
{
    // Tail amounts: generic, exact fractions k/N of the sample count (where
    // fl(k/N)*N can round below k and the truncation drops one), 0, 1 and
    // the out-of-range values that run the loops off either end.
    double DrawTail(oracle::Ctx& io, int32_t total)
    {
        switch (io.rawInteger(0, 4))
        {
        case 0: return io.given(io.raw(0.0, 0.5));
        case 1: return io.given(static_cast<double>(io.rawInteger(0, total)) / std::max(total, 1));
        case 2: return io.given(io.rawInteger(0, 1));
        case 3: return io.given(io.raw(1.0, 1.5));
        default: return io.given(io.raw(-0.5, 0.0));
        }
    }

    void OutHistogram(oracle::Ctx& io, Histogram& h, int32_t total)
    {
        auto const& buckets = h.GetBuckets();
        io.outInt(buckets.size());
        for (auto b : buckets) { io.outInt(b); }
        io.outInt(h.GetExcessLess());
        io.outInt(h.GetExcessGreater());
        for (int32_t k = 0; k < 3; ++k)
        {
            double tail = DrawTail(io, total);
            io.outInt(h.GetLowerTail(tail));
            io.outInt(h.GetUpperTail(tail));
            int32_t lower = -7, upper = -7;
            h.GetTails(tail, lower, upper);
            io.outInt(lower);
            io.outInt(upper);
        }
    }

    // Integer samples: small ranges with many ties, full-range int32 values,
    // a constant image, two values only, and values around [0, numBuckets).
    std::vector<int32_t> DrawIntSamples(oracle::Ctx& io, int32_t numBuckets, int32_t n)
    {
        int mode = io.rawInteger(0, 4);
        int32_t c0 = io.rawInteger(-100, 100), c1 = io.rawInteger(-100, 100);
        std::vector<int32_t> s(n);
        for (int32_t i = 0; i < n; ++i)
        {
            int32_t v;
            switch (mode)
            {
            case 0: v = io.rawInteger(-5, 5); break;
            case 1: v = io.rawInteger(std::numeric_limits<int32_t>::min(),
                std::numeric_limits<int32_t>::max()); break;
            case 2: v = c0; break;
            case 3: v = (io.rawInteger(0, 1) == 0 ? c0 : c1); break;
            default: v = io.rawInteger(-3, numBuckets + 3); break;
            }
            s[i] = static_cast<int32_t>(io.given(v));
        }
        return s;
    }
}

ORACLE_CASE("Histogram.int.noRescaling")
{
    int32_t numBuckets = io.integer(1, 12);
    int32_t n = io.integer(1, 40);
    std::vector<int32_t> samples = DrawIntSamples(io, numBuckets, n);
    Histogram h(numBuckets, n, samples.data(), true);
    OutHistogram(io, h, n);
}

ORACLE_CASE("Histogram.int.rescaled")
{
    int32_t numBuckets = io.integer(1, 12);
    int32_t n = io.integer(1, 40);
    std::vector<int32_t> samples = DrawIntSamples(io, numBuckets, n);
    Histogram h(numBuckets, n, samples.data(), false);
    OutHistogram(io, h, n);
}

namespace
{
    // Real samples: uniform, a lattice with ties, a constant image (both
    // zeros included), dyadic values, values spread over many binades, and
    // a NaN first sample (every comparison with the running extremes is
    // false, so upstream takes the constant branch). NaNs elsewhere and
    // infinities make static_cast<int32_t>(NaN) undefined behaviour and a
    // subnormal range overflows 'mult' (finding #436), so neither is drawn.
    std::vector<double> DrawRealSamples(oracle::Ctx& io, int32_t n, bool asFloat)
    {
        int mode = io.rawInteger(0, 5);
        double c = io.raw(-3.0, 3.0);
        std::vector<double> s(n);
        for (int32_t i = 0; i < n; ++i)
        {
            double v;
            switch (mode)
            {
            case 0: v = io.raw(-10.0, 10.0); break;
            case 1: v = io.rawInteger(-3, 3); break;
            case 2: v = (io.rawInteger(0, 1) == 0 ? 0.0 : -0.0); break;
            case 3: v = io.rawInteger(-64, 64) / 32.0; break;
            case 4: v = std::ldexp(io.raw(-1.0, 1.0), io.rawInteger(-40, 40)); break;
            default: v = (i == 0 ? std::numeric_limits<double>::quiet_NaN() : c); break;
            }
            if (asFloat) { v = static_cast<double>(static_cast<float>(v)); }
            s[i] = io.given(v);
        }
        return s;
    }
}

ORACLE_CASE("Histogram.double")
{
    int32_t numBuckets = io.integer(1, 12);
    int32_t n = io.integer(1, 40);
    std::vector<double> samples = DrawRealSamples(io, n, false);
    Histogram h(numBuckets, n, samples.data());
    OutHistogram(io, h, n);
}

// The float constructor; the samples are recorded as the doubles of the
// floats, which is what upstream computes with.
ORACLE_CASE("Histogram.float")
{
    int32_t numBuckets = io.integer(1, 12);
    int32_t n = io.integer(1, 40);
    std::vector<double> samples = DrawRealSamples(io, n, true);
    std::vector<float> fsamples(n);
    for (int32_t i = 0; i < n; ++i) { fsamples[i] = static_cast<float>(samples[i]); }
    Histogram h(numBuckets, n, fsamples.data());
    OutHistogram(io, h, n);
}

// Histogram(numBuckets) followed by InsertCheck (any value) and Insert
// (in-range values only: Insert does no bounds check).
ORACLE_CASE("Histogram.incremental")
{
    int32_t numBuckets = io.integer(1, 12);
    Histogram h(numBuckets);
    int32_t numOps = io.integer(0, 30);
    for (int32_t k = 0; k < numOps; ++k)
    {
        bool checked = io.boolean();
        if (checked)
        {
            int32_t v = io.integer(-3, numBuckets + 3);
            h.InsertCheck(v);
        }
        else
        {
            int32_t v = io.integer(0, numBuckets - 1);
            h.Insert(v);
        }
    }
    OutHistogram(io, h, numOps);
}

// The LogAssert of every constructor: numBuckets <= 0 or numSamples <= 0
// (a negative bucket count makes the std::vector member throw first; both
// are exceptions).
ORACLE_CASE("Histogram.invalid")
{
    int32_t which = io.integer(0, 4);
    int32_t numBuckets = io.integer(-2, 3);
    int32_t n = io.integer(0, 3);
    std::vector<int32_t> is(std::max(n, 1), 1);
    std::vector<float> fs(std::max(n, 1), 1.0f);
    std::vector<double> ds(std::max(n, 1), 1.0);
    std::unique_ptr<Histogram> h;
    switch (which)
    {
    case 0: h = std::make_unique<Histogram>(numBuckets, n, is.data(), true); break;
    case 1: h = std::make_unique<Histogram>(numBuckets, n, is.data(), false); break;
    case 2: h = std::make_unique<Histogram>(numBuckets, n, fs.data()); break;
    case 3: h = std::make_unique<Histogram>(numBuckets, n, ds.data()); break;
    default: h = std::make_unique<Histogram>(numBuckets); break;
    }
    io.outInt(h->GetBuckets().size());
}

// ---- PdeFilter ----

namespace
{
    // PdeFilter is abstract and its concrete filters belong to group 25. The
    // base has code of its own (the constructor's range scan and scale
    // selection, the accessors, Update's call order); this subclass exposes
    // it, and the replay defines the same subclass over the port.
    class PdeProbe : public PdeFilter<double>
    {
    public:
        PdeProbe(int32_t quantity, double const* data, double borderValue, ScaleType scaleType)
            :
            PdeFilter<double>(quantity, data, borderValue, scaleType)
        {
        }

        double Min() const { return mMin; }
        double Offset() const { return mOffset; }
        double Scale() const { return mScale; }
        std::vector<int32_t> calls;

    protected:
        virtual void OnPreUpdate() override { calls.push_back(1); }
        virtual void OnUpdate() override { calls.push_back(2); }
        virtual void OnPostUpdate() override { calls.push_back(3); }
    };

    __declspec(noinline) PdeProbe* MakePdeProbe(int32_t quantity, double const* data,
        double borderValue, int32_t scaleType)
    {
        return new PdeProbe(quantity, data, borderValue,
            static_cast<PdeFilter<double>::ScaleType>(scaleType));
    }
}

// Data: uniform, all positive, all negative, a lattice with ties, a
// constant image, signed zeros, dyadic values with the PRESERVE_ZERO tie
// max == -min, infinities and NaNs. The scale type 4 is outside the enum's
// named values (the switch falls through, offset and scale stay 0).
ORACLE_CASE("PdeFilter.construct")
{
    int32_t quantity = io.integer(1, 12);
    int mode = io.rawInteger(0, 8);
    double c = io.raw(-5.0, 5.0);
    std::vector<double> data(quantity);
    for (int32_t i = 0; i < quantity; ++i)
    {
        double v;
        switch (mode)
        {
        case 0: v = io.raw(-10.0, 10.0); break;
        case 1: v = io.raw(0.5, 10.0); break;
        case 2: v = io.raw(-10.0, -0.5); break;
        case 3: v = io.rawInteger(-3, 3); break;
        case 4: v = c; break;
        case 5: v = (io.rawInteger(0, 1) == 0 ? 0.0 : -0.0); break;
        case 6: v = (i % 2 == 0 ? 1.0 : -1.0) * io.rawInteger(0, 4) * 0.25; break;
        case 7:
        {
            int k = io.rawInteger(0, 4);
            v = (k == 0 ? std::numeric_limits<double>::infinity()
                : (k == 1 ? -std::numeric_limits<double>::infinity() : io.raw(-3.0, 3.0)));
            break;
        }
        default:
            v = (io.rawInteger(0, 3) == 0 ? std::numeric_limits<double>::quiet_NaN()
                : io.raw(-3.0, 3.0));
            break;
        }
        data[i] = io.given(v);
    }
    double borderValue = io.given(io.rawInteger(0, 3) == 0
        ? std::numeric_limits<double>::max() : io.raw(-2.0, 2.0));
    int32_t scaleType = io.integer(0, 4);
    std::unique_ptr<PdeProbe> filter(MakePdeProbe(quantity, data.data(), borderValue, scaleType));
    io.outInt(filter->GetQuantity());
    io.outReal(filter->GetBorderValue());
    io.outInt(static_cast<int32_t>(filter->GetScaleType()));
    io.outReal(filter->Min());
    io.outReal(filter->Offset());
    io.outReal(filter->Scale());
    io.outReal(filter->GetTimeStep());
    double timeStep = io.real(-1.0, 1.0);
    filter->SetTimeStep(timeStep);
    io.outReal(filter->GetTimeStep());
    int32_t numUpdates = io.integer(0, 3);
    for (int32_t k = 0; k < numUpdates; ++k) { filter->Update(); }
    io.outInt(filter->calls.size());
    for (auto v : filter->calls) { io.outInt(v); }
}

// ---- FastMarch ----

namespace
{
    // FastMarch is abstract; FastMarch2/3 belong to group 25. This 1-D
    // subclass, defined identically in the replay, drives the base's own
    // code: both constructors (seed and zero-speed marking, inverted
    // speeds), the protected MinHeap<size_t, Real> and trial records, and
    // every public classifier. Its arrival time is the smaller valid
    // neighbour time plus the inverse speed (arithmetic only).
    class FastMarch1 : public FastMarch<double>
    {
    public:
        FastMarch1(size_t quantity, std::vector<size_t> const& seeds,
            std::vector<double> const& speeds)
            :
            FastMarch<double>(quantity, seeds, speeds)
        {
            Initialize();
        }

        FastMarch1(size_t quantity, std::vector<size_t> const& seeds, double speed)
            :
            FastMarch<double>(quantity, seeds, speed)
        {
            Initialize();
        }

        virtual void GetBoundary(std::vector<size_t>& boundary) const override
        {
            for (size_t i = 0; i < mQuantity; ++i)
            {
                if (IsBoundary(i))
                {
                    boundary.push_back(i);
                }
            }
        }

        virtual bool IsBoundary(size_t i) const override
        {
            return IsValid(i) && !IsTrial(i)
                && ((i > 0 && IsTrial(i - 1)) || (i + 1 < mQuantity && IsTrial(i + 1)));
        }

        // Upstream FastMarch2/3 call Remove unconditionally; the subclass
        // stops on an empty heap instead (both sides).
        virtual void Iterate() override
        {
            removed = false;
            if (mHeap.GetNumElements() == 0)
            {
                return;
            }
            size_t i = 0;
            double value = 0.0;
            mHeap.Remove(i, value);
            removed = true;
            removedKey = i;
            removedValue = value;
            mTrials[i] = nullptr;
            if (i > 0) { Visit(i - 1); }
            if (i + 1 < mQuantity) { Visit(i + 1); }
        }

        int32_t NumTrials() const { return mHeap.GetNumElements(); }
        double InvSpeed(size_t i) const { return mInvSpeeds[i]; }

        bool removed = false;
        size_t removedKey = 0;
        double removedValue = 0.0;

    private:
        void Initialize()
        {
            for (size_t i = 0; i < mQuantity; ++i)
            {
                if (IsFar(i)
                    && ((i > 0 && IsValid(i - 1) && !IsTrial(i - 1))
                        || (i + 1 < mQuantity && IsValid(i + 1) && !IsTrial(i + 1))))
                {
                    ComputeTime(i);
                    mTrials[i] = mHeap.Insert(i, mTimes[i]);
                }
            }
        }

        void Visit(size_t j)
        {
            if (IsTrial(j))
            {
                ComputeTime(j);
                mHeap.Update(mTrials[j], mTimes[j]);
            }
            else if (IsFar(j))
            {
                ComputeTime(j);
                mTrials[j] = mHeap.Insert(j, mTimes[j]);
            }
        }

        void ComputeTime(size_t i)
        {
            double t = std::numeric_limits<double>::max();
            if (i > 0 && IsValid(i - 1) && mTimes[i - 1] < t) { t = mTimes[i - 1]; }
            if (i + 1 < mQuantity && IsValid(i + 1) && mTimes[i + 1] < t) { t = mTimes[i + 1]; }
            mTimes[i] = t + mInvSpeeds[i];
        }
    };

    void OutFastMarchState(oracle::Ctx& io, FastMarch1 const& march)
    {
        size_t const q = march.GetQuantity();
        io.outInt(march.NumTrials());
        for (size_t i = 0; i < q; ++i)
        {
            io.outReal(march.GetTime(i));
            // IsValid, IsTrial, IsFar, IsZeroSpeed, IsInterior, IsBoundary
            // as bits 0..5.
            io.outInt((march.IsValid(i) ? 1 : 0) | (march.IsTrial(i) ? 2 : 0)
                | (march.IsFar(i) ? 4 : 0) | (march.IsZeroSpeed(i) ? 8 : 0)
                | (march.IsInterior(i) ? 16 : 0) | (march.IsBoundary(i) ? 32 : 0));
        }
        std::vector<size_t> interior, boundary;
        march.GetInterior(interior);
        march.GetBoundary(boundary);
        io.outInt(interior.size());
        for (auto v : interior) { io.outInt(v); }
        io.outInt(boundary.size());
        for (auto v : boundary) { io.outInt(v); }
        double minValue = 0.0, maxValue = 0.0;
        march.GetTimeExtremes(minValue, maxValue);
        io.outReal(minValue);
        io.outReal(maxValue);
    }

    // Seeds (duplicates allowed, possibly on zero-speed pixels) and either
    // per-pixel speeds (positive, zero, negative, NaN: all but the positive
    // ones mark the pixel zero-speed) or one constant speed (lattice speeds
    // give tied arrival times, so the heap's tie order decides the removal
    // sequence; 0 gives infinite times, a negative speed negative ones).
    std::unique_ptr<FastMarch1> DrawFastMarch(oracle::Ctx& io, bool allowNaN)
    {
        int32_t quantity = io.integer(1, 16);
        int32_t numSeeds = io.integer(0, 3);
        std::vector<size_t> seeds(numSeeds);
        for (int32_t k = 0; k < numSeeds; ++k) { seeds[k] = io.integer(0, quantity - 1); }
        bool perPixel = io.boolean();
        if (perPixel)
        {
            std::vector<double> speeds(quantity);
            for (int32_t i = 0; i < quantity; ++i)
            {
                int k = io.rawInteger(0, 9);
                double v = (k == 0 ? 0.0 : (k == 1 ? -1.0 : (k == 2 ? std::nan("")
                    : (k < 6 ? io.rawInteger(1, 4) * 0.5 : io.raw(0.1, 4.0)))));
                speeds[i] = io.given(v);
            }
            return std::make_unique<FastMarch1>(quantity, seeds, speeds);
        }
        int k = io.rawInteger(0, 9);
        double speed = (k == 0 ? 0.0 : (k == 1 ? -2.0 : (k < 6 ? io.rawInteger(1, 4) * 0.5
            : (allowNaN && k == 9 ? std::nan("") : io.raw(0.1, 4.0)))));
        io.given(speed);
        return std::make_unique<FastMarch1>(quantity, seeds, speed);
    }

    void FastMarchCase(oracle::Ctx& io, bool allowNaN)
    {
        std::unique_ptr<FastMarch1> march = DrawFastMarch(io, allowNaN);
        size_t const q = march->GetQuantity();
        for (size_t i = 0; i < q; ++i) { io.outReal(march->InvSpeed(i)); }
        OutFastMarchState(io, *march);
        int32_t numIterations = io.integer(0, 20);
        for (int32_t k = 0; k < numIterations; ++k)
        {
            march->Iterate();
            io.outBool(march->removed);
            if (march->removed)
            {
                io.outInt(march->removedKey);
                io.outReal(march->removedValue);
            }
        }
        OutFastMarchState(io, *march);
    }
}

ORACLE_CASE("FastMarch.march1") { FastMarchCase(io, false); }

// A NaN constant speed gives NaN inverse speeds and NaN arrival times in the
// heap, which is where MinHeap's NaN comparisons matter.
ORACLE_CASE("FastMarch.march1.nanSpeed") { FastMarchCase(io, true); }

// SetTime with arbitrary values (both zeros, max and -max, infinities, NaN,
// negatives) followed by every classifier and GetTimeExtremes, including
// the all-invalid case that returns (max, -max).
ORACLE_CASE("FastMarch.accessors")
{
    std::unique_ptr<FastMarch1> march = DrawFastMarch(io, false);
    size_t const q = march->GetQuantity();
    int32_t numSets = io.integer(0, 8);
    for (int32_t k = 0; k < numSets; ++k)
    {
        size_t i = io.integer(0, static_cast<int32_t>(q) - 1);
        double const big = std::numeric_limits<double>::max();
        double const inf = std::numeric_limits<double>::infinity();
        double const values[9] = { 0.0, -0.0, big, -big, inf, -inf, std::nan(""), -1.5, 2.5 };
        int32_t which = io.rawInteger(0, 11);
        double t = io.given(which < 9 ? values[which] : io.raw(0.0, 10.0));
        march->SetTime(i, t);
    }
    OutFastMarchState(io, *march);
}

// ---- MarchingCubes ----

// The whole 256-entry table: record r covers entries 13*(r mod 20) through
// 13*(r mod 20) + 12, so the 20 committed records cover every entry and
// entries 256..259, where only GetConfigurationType is defined (it returns
// ""). Per entry: numVertices, numTriangles, the 12 vertex pairs packed as
// a*64 + b and the 5 triples packed as (i*64 + j)*64 + k (padding
// included), a digest of the entry's 41 values in the flat GetTable()
// storage, a digest of the GetPrebuiltTable() row, and the configuration
// type name. One more GetConfigurationType call takes an entry far out of
// range.
ORACLE_CASE("MarchingCubes.table")
{
    int32_t first = static_cast<int32_t>(io.given(13 * (io.index() % 20)));
    MarchingCubes<int32_t> mc;
    int32_t const* flat = mc.GetTable();
    auto const& prebuilt = MarchingCubes<int32_t>::GetPrebuiltTable();
    for (int32_t entry = first; entry < first + 13; ++entry)
    {
        if (entry < 256)
        {
            auto const& topology = mc.GetTable(static_cast<size_t>(entry));
            io.outInt(topology.numVertices);
            io.outInt(topology.numTriangles);
            for (auto const& p : topology.vpair) { io.outInt(p[0] * 64 + p[1]); }
            for (auto const& t : topology.itriple) { io.outInt((t[0] * 64 + t[1]) * 64 + t[2]); }
            Digest flatDigest, prebuiltDigest;
            for (int32_t k = 0; k < 41; ++k)
            {
                flatDigest.integer(flat[41 * entry + k]);
                prebuiltDigest.integer(prebuilt[entry][k]);
            }
            io.outInt(flatDigest.h);
            io.outInt(prebuiltDigest.h);
        }
        OutString(io, MarchingCubes<int32_t>::GetConfigurationType(static_cast<size_t>(entry)));
    }
    int32_t far = io.integer(256, 1 << 30);
    OutString(io, MarchingCubes<int32_t>::GetConfigurationType(static_cast<size_t>(far)));
}

// ---- AdaptiveSkeletonClimbing2 ----

namespace
{
    using ASCNumber = BSNumber<UIntegerAP32>;

    // A unit cell whose four edges are all crossed becomes a 1x1 case-15
    // rectangle (a saddle cell never merges: its ymin and ymax crossings run
    // in opposite directions, so neither DoXMerge nor DoYMerge finds it
    // strongly mono). Upstream pairs the four crossings by the sign of
    // det = i00*i11 - i01*i10: det > 0 cuts off corners 00 and 11, det < 0
    // corners 10 and 01, det = 0 adds the plus-sign branch point. The
    // bilinear interpolant's topology is decided by
    // dg = (i00-L)(i11-L) - (i01-L)(i10-L) = det - L*(i00+i11-i01-i10):
    // dg > 0 connects 00 and 11 through the saddle (cut off 10 and 01),
    // dg < 0 cuts off 00 and 11, dg = 0 is the plus sign. Upstream is right
    // iff sign(det) = -sign(dg). Returns -1 for no saddle cell or a corner on
    // the level (GetInterp throws on both sides first), 0 when upstream's
    // pairing is right, 1 when it is not. dg is evaluated exactly.
    int SaddleDefect(int64_t i00, int64_t i10, int64_t i01, int64_t i11, double level)
    {
        bool a00 = static_cast<double>(i00) > level, a10 = static_cast<double>(i10) > level;
        bool a01 = static_cast<double>(i01) > level, a11 = static_cast<double>(i11) > level;
        if (!(a00 == a11 && a10 == a01 && a00 != a10))
        {
            return -1;
        }
        if (static_cast<double>(i00) == level || static_cast<double>(i10) == level
            || static_cast<double>(i01) == level || static_cast<double>(i11) == level)
        {
            return -1;
        }
        int64_t det = i00 * i11 - i01 * i10;
        ASCNumber dg = ASCNumber(det) - ASCNumber(level) * ASCNumber(i00 + i11 - i01 - i10);
        int sdg = dg.GetSign();
        int sdet = (det > 0 ? 1 : (det < 0 ? -1 : 0));
        return (sdet == -sdg ? 0 : 1);
    }

    // [0] saddle cells where upstream is right, [1] where it is not.
    std::array<int, 2> CountSaddles(std::vector<int64_t> const& px, int size, double level)
    {
        std::array<int, 2> count{ 0, 0 };
        for (int y = 0; y + 1 < size; ++y)
        {
            for (int x = 0; x + 1 < size; ++x)
            {
                int64_t i00 = px[x + size * y], i10 = px[x + 1 + size * y];
                int64_t i01 = px[x + size * (y + 1)], i11 = px[x + 1 + size * (y + 1)];
                int d = SaddleDefect(i00, i10, i01, i11, level);
                if (d >= 0) { ++count[d]; }
            }
        }
        return count;
    }

    // An unrecorded image in [lo, hi]: small random integers (ties and
    // saddles), a ramp, a paraboloid, a constant, a signed checkerboard
    // (saddle-rich) and full-range values.
    std::vector<int64_t> RawImage(oracle::Ctx& io, int size, int64_t lo, int64_t hi)
    {
        int mode = io.rawInteger(0, 5);
        int64_t mid = (lo < 0 ? 0 : (lo + hi) / 2);
        int64_t span = hi - mid;
        int R = (io.rawInteger(0, 2) == 0 ? 2 : (io.rawInteger(0, 1) == 0 ? 5 : 100));
        int a = io.rawInteger(-3, 3), b = io.rawInteger(-3, 3), c = io.rawInteger(-5, 5);
        int cx = io.rawInteger(0, size - 1), cy = io.rawInteger(0, size - 1);
        int sgn = (io.rawInteger(0, 1) == 0 ? 1 : -1);
        std::vector<int64_t> px(static_cast<size_t>(size) * size);
        for (int y = 0; y < size; ++y)
        {
            for (int x = 0; x < size; ++x)
            {
                int64_t v;
                switch (mode)
                {
                case 0: v = io.rawInteger(-R, R); break;
                case 1: v = a * x + b * y + c; break;
                case 2: v = sgn * ((x - cx) * (x - cx) + (y - cy) * (y - cy)) + c; break;
                case 3: v = c; break;
                case 4: v = ((x + y) % 2 == 0 ? 1 : -1) * io.rawInteger(1, R) + io.rawInteger(-1, 1); break;
                default: v = static_cast<int64_t>(io.raw(-1.0, 1.0) * static_cast<double>(span)); break;
                }
                px[x + size * y] = std::min(std::max(mid + v, lo), hi);
            }
        }
        return px;
    }

    // An unrecorded level: half-integers (the documented use), reals, the
    // integers of the pixel range (a pixel on the level throws in GetInterp
    // when it ends a crossed edge), 0, values outside the range (nothing
    // extracted), quarter-integers and NaN (every comparison false).
    double RawLevel(oracle::Ctx& io, std::vector<int64_t> const& px)
    {
        auto mm = std::minmax_element(px.begin(), px.end());
        double lo = static_cast<double>(*mm.first), hi = static_cast<double>(*mm.second);
        switch (io.rawInteger(0, 9))
        {
        case 0: case 1: case 2:
            return std::floor(io.raw(lo - 1.0, hi + 1.0)) + 0.5;
        case 3: case 4:
            return io.raw(lo - 0.5, hi + 0.5);
        case 5:
            return std::floor(io.raw(lo, hi + 1.0));
        case 6:
            return 0.0;
        case 7:
            return (io.rawInteger(0, 1) == 0 ? hi + io.raw(0.25, 3.0) : lo - io.raw(0.25, 3.0));
        case 8:
            return std::floor(4.0 * io.raw(lo - 1.0, hi + 1.0)) / 4.0;
        default:
            return (io.rawInteger(0, 3) == 0 ? std::nan("") : 0.5);
        }
    }

    struct ASCResult
    {
        std::vector<std::array<double, 2>> v, vu;
        std::vector<std::array<int32_t, 2>> e, eu;
    };

    template <typename T>
    __declspec(noinline) std::vector<ASCResult> RunASC(int32_t N, std::vector<int64_t> const& px,
        std::vector<double> const& levels, std::vector<int32_t> const& depths)
    {
        std::vector<T> pixels(px.size());
        for (size_t i = 0; i < px.size(); ++i) { pixels[i] = static_cast<T>(px[i]); }
        AdaptiveSkeletonClimbing2<T, double> asc(N, pixels.data());
        std::vector<ASCResult> results(levels.size());
        for (size_t k = 0; k < levels.size(); ++k)
        {
            ASCResult& r = results[k];
            asc.Extract(levels[k], depths[k], r.v, r.e);
            r.vu = r.v;
            r.eu = r.e;
            asc.MakeUnique(r.vu, r.eu);
        }
        return results;
    }

    void OutCurves(oracle::Ctx& io, std::vector<std::array<double, 2>> const& v,
        std::vector<std::array<int32_t, 2>> const& e)
    {
        io.outInt(v.size());
        for (auto const& p : v) { io.outReal(p[0]); io.outReal(p[1]); }
        io.outInt(e.size());
        for (auto const& s : e) { io.outInt(s[0]); io.outInt(s[1]); }
    }

    void OutASCResults(oracle::Ctx& io, std::vector<ASCResult> const& results)
    {
        for (auto const& r : results)
        {
            OutCurves(io, r.v, r.e);
            OutCurves(io, r.vu, r.eu);
        }
    }

    // N = 1, 2, 3 or 4 (the last one in six records: 289 pixels).
    int32_t DrawN(oracle::Ctx& io, int32_t maxN)
    {
        int32_t k = io.rawInteger(0, 5);
        int32_t N = (k < 2 ? 1 : (k < 4 ? 2 : (k == 4 ? 3 : 4)));
        return static_cast<int32_t>(io.given(std::min(N, maxN)));
    }

    // Pixels in [lo, hi] and 1 or 2 levels on which upstream's saddle
    // pairing is sound (see SaddleDefect; the port pairs by the interpolant,
    // the deviation case covers the rest). Everything the acceptance test
    // depends on is redrawn per attempt; after 64 attempts the image falls
    // back to a ramp, which has no saddle cell.
    void DrawSoundImage(oracle::Ctx& io, int32_t N, int64_t lo, int64_t hi,
        std::vector<int64_t>& px, std::vector<double>& levels)
    {
        int const size = (1 << N) + 1;
        int32_t numExtracts = io.integer(1, 2);
        levels.resize(numExtracts);
        bool sound = false;
        for (int attempt = 0; attempt < 64 && !sound; ++attempt)
        {
            px = RawImage(io, size, lo, hi);
            sound = true;
            for (auto& level : levels)
            {
                level = RawLevel(io, px);
                sound = sound && CountSaddles(px, size, level)[1] == 0;
            }
        }
        if (!sound)
        {
            for (int y = 0; y < size; ++y)
            {
                for (int x = 0; x < size; ++x) { px[x + size * y] = std::min<int64_t>(lo + x + y, hi); }
            }
            for (size_t k = 0; k < levels.size(); ++k) { levels[k] = static_cast<double>(lo) + 0.5 + k; }
        }
        for (auto v : px) { io.given(static_cast<double>(v)); }
    }

    template <typename T>
    void ASCCase(oracle::Ctx& io, int32_t N, int64_t lo, int64_t hi)
    {
        std::vector<int64_t> px;
        std::vector<double> levels;
        DrawSoundImage(io, N, lo, hi, px, levels);
        std::vector<int32_t> depths(levels.size());
        for (size_t k = 0; k < levels.size(); ++k)
        {
            io.given(levels[k]);
            depths[k] = io.integer(-2, N + 1);
        }
        std::vector<ASCResult> results = RunASC<T>(N, px, levels, depths);
        OutASCResults(io, results);
    }
}

// int32_t pixels, N = 1..4, 1 or 2 extractions on the same object (the merge
// trees are rebuilt per call), depth -2..N+1 (merging blocked above the
// depth), levels on and between the sample values. Every output list is
// emitted verbatim: upstream's order is the deterministic recursion order of
// GetRectangles, and MakeUnique numbers by first occurrence.
ORACLE_CASE("AdaptiveSkeletonClimbing2.extract")
{
    int32_t N = DrawN(io, 4);
    ASCCase<int32_t>(io, N, std::numeric_limits<int32_t>::min(), std::numeric_limits<int32_t>::max());
}

// The other pixel types upstream allows. uint32_t values stay below 2^31:
// above that, det = i00*i11 - i01*i10 overflows int64_t (undefined
// behaviour; see the group report).
ORACLE_CASE("AdaptiveSkeletonClimbing2.extract.types")
{
    int32_t type = io.integer(0, 4);
    int32_t N = DrawN(io, 3);
    switch (type)
    {
    case 0: ASCCase<int8_t>(io, N, -128, 127); break;
    case 1: ASCCase<uint8_t>(io, N, 0, 255); break;
    case 2: ASCCase<int16_t>(io, N, -32768, 32767); break;
    case 3: ASCCase<uint16_t>(io, N, 0, 65535); break;
    default: ASCCase<uint32_t>(io, N, 0, int64_t(1) << 31); break;
    }
}

namespace
{
    // A saddle-rich image (signed checkerboard with random magnitudes, or
    // random small values) with one level, optionally with a constructed
    // plus-sign cell (-p, p; q, -q around 0: det = 0) or a cell whose level
    // is det/S (the interpolant's saddle exactly on the level).
    void RawSaddleImage(oracle::Ctx& io, int size, std::vector<int64_t>& px, double& level)
    {
        int mode = io.rawInteger(0, 3);
        int R = io.rawInteger(2, 9);
        px.assign(static_cast<size_t>(size) * size, 0);
        for (int y = 0; y < size; ++y)
        {
            for (int x = 0; x < size; ++x)
            {
                px[x + size * y] = (mode == 0 ? io.rawInteger(-R, R)
                    : ((x + y) % 2 == 0 ? 1 : -1) * io.rawInteger(1, R) + io.rawInteger(-2, 2));
            }
        }
        level = std::floor(io.raw(-2.0, 2.0)) + 0.5;
        int cx = io.rawInteger(0, size - 2), cy = io.rawInteger(0, size - 2);
        int64_t& f00 = px[cx + size * cy];
        int64_t& f10 = px[cx + 1 + size * cy];
        int64_t& f01 = px[cx + size * (cy + 1)];
        int64_t& f11 = px[cx + 1 + size * (cy + 1)];
        if (mode == 2)
        {
            int64_t p = io.rawInteger(1, R), q = io.rawInteger(1, R);
            f00 = -p; f10 = p; f01 = q; f11 = -q;
            level = (io.rawInteger(0, 1) == 0 ? 0.0 : io.raw(-0.9, 0.9));
        }
        else if (mode == 3)
        {
            int64_t det = f00 * f11 - f01 * f10;
            int64_t S = f00 + f11 - f01 - f10;
            if (S != 0) { level = static_cast<double>(det) / static_cast<double>(S); }
        }
    }

    // wantDefect = false: at least one saddle cell and upstream sound on all
    // of them. wantDefect = true: at least one saddle cell where upstream's
    // pairing contradicts the interpolant, and no pixel on the level (so
    // neither side throws). Capped; the fallbacks are fixed configurations
    // of the required kind (cell 00 = 10, -1; -1, 10 at level 5.5 is sound,
    // at level 0.5 defective).
    void DrawSaddleImage(oracle::Ctx& io, bool wantDefect, int32_t& N,
        std::vector<int64_t>& px, double& level)
    {
        N = static_cast<int32_t>(io.given(io.index() % 4 == 3 ? 2 : 1));
        int const size = (1 << N) + 1;
        io.integer(1, 1);  // numExtracts, the layout of the other cases
        bool ok = false;
        for (int attempt = 0; attempt < 256 && !ok; ++attempt)
        {
            RawSaddleImage(io, size, px, level);
            auto count = CountSaddles(px, size, level);
            bool onLevel = std::any_of(px.begin(), px.end(),
                [level](int64_t v) { return static_cast<double>(v) == level; });
            ok = !onLevel && (wantDefect ? count[1] > 0 : count[0] > 0 && count[1] == 0);
        }
        if (!ok)
        {
            px.assign(static_cast<size_t>(size) * size, 10);
            px[1] = -1;
            px[size] = -1;
            level = (wantDefect ? 0.5 : 5.5);
        }
        for (auto v : px) { io.given(static_cast<double>(v)); }
    }

    void SaddleCase(oracle::Ctx& io, bool wantDefect)
    {
        int32_t N = 0;
        std::vector<int64_t> px;
        double level = 0.0;
        DrawSaddleImage(io, wantDefect, N, px, level);
        io.given(level);
        std::vector<int32_t> depths{ io.integer(-2, N + 1) };
        std::vector<ASCResult> results = RunASC<int32_t>(N, px, { level }, depths);
        OutASCResults(io, results);
    }
}

// Saddle cells (the case-15 rectangle) on which upstream's pairing agrees
// with the bilinear interpolant: disjoint pairs of both kinds and the
// plus-sign branch point (det = 0 at level 0).
ORACLE_CASE("AdaptiveSkeletonClimbing2.extract.saddle") { SaddleCase(io, false); }

// Every record has a saddle cell on which upstream's det = i00*i11 - i01*i10
// pairing contradicts the bilinear interpolant (it ignores the level and its
// two disjoint pairings are swapped; see the group report). The port pairs
// by dg = det - L*(i00+i11-i01-i10), exactly.
ORACLE_CASE("AdaptiveSkeletonClimbing2.extract.saddlePairing") { SaddleCase(io, true); }

namespace
{
    // The pixel formula of the large case (mirrored in the replay). kind 0:
    // a quadratic a*u^2 + b*v^2 + c*u*v + d*x + e*y with u = x - cx,
    // v = y - cy; kind 1: a*|u| + b*|v| + c*min(|u|,|v|) + d*x + e*y.
    int64_t LargePixel(int64_t kind, std::array<int64_t, 7> const& p, int64_t x, int64_t y)
    {
        int64_t u = x - p[5], v = y - p[6];
        if (kind == 0)
        {
            return p[0] * u * u + p[1] * v * v + p[2] * u * v + p[3] * x + p[4] * y;
        }
        int64_t au = (u < 0 ? -u : u), av = (v < 0 ? -v : v);
        return p[0] * au + p[1] * av + p[2] * std::min(au, av) + p[3] * x + p[4] * y;
    }

    void DigestCurves(oracle::Ctx& io, std::vector<std::array<double, 2>> const& v,
        std::vector<std::array<int32_t, 2>> const& e)
    {
        Digest dv, de;
        for (auto const& p : v) { dv.real(p[0]); dv.real(p[1]); }
        for (auto const& s : e) { de.integer(s[0]); de.integer(s[1]); }
        io.outInt(v.size());
        io.outInt(e.size());
        io.outInt(dv.h);
        io.outInt(de.h);
    }
}

// 33x33 and 65x65 images from a recorded closed-form integer formula, so the
// merge trees are deep and large monotone regions merge into big
// rectangles. The outputs are the counts and FNV-1a digests of the vertex
// bit patterns and edge indices, before and after MakeUnique (hundreds of
// vertices per record). Upstream-sound saddles only (capped rejection,
// fallback: a paraboloid bowl).
ORACLE_CASE("AdaptiveSkeletonClimbing2.extract.large")
{
    int32_t N = static_cast<int32_t>(io.given(io.rawInteger(5, 6)));
    int const size = (1 << N) + 1;
    int64_t kind = 0;
    std::array<int64_t, 7> p{};
    std::vector<int64_t> px(static_cast<size_t>(size) * size);
    double level = 0.0;
    bool sound = false;
    for (int attempt = 0; attempt < 64 && !sound; ++attempt)
    {
        kind = io.rawInteger(0, 1);
        for (int k = 0; k < 3; ++k) { p[k] = io.rawInteger(-3, 3); }
        for (int k = 3; k < 5; ++k) { p[k] = io.rawInteger(-10, 10); }
        for (int k = 5; k < 7; ++k) { p[k] = io.rawInteger(0, size - 1); }
        for (int y = 0; y < size; ++y)
        {
            for (int x = 0; x < size; ++x) { px[x + size * y] = LargePixel(kind, p, x, y); }
        }
        level = RawLevel(io, px);
        // A pixel on the level throws (covered by the main case); here the
        // record would carry nothing.
        sound = CountSaddles(px, size, level)[1] == 0
            && std::none_of(px.begin(), px.end(),
                [level](int64_t v) { return static_cast<double>(v) == level; });
    }
    if (!sound)
    {
        kind = 0;
        p = { 1, 1, 0, 0, 0, size / 2, size / 2 };
        level = 100.5;
    }
    io.given(static_cast<double>(kind));
    for (auto v : p) { io.given(static_cast<double>(v)); }
    io.given(level);
    std::vector<int32_t> depths{ io.integer(-2, N + 1) };
    for (int y = 0; y < size; ++y)
    {
        for (int x = 0; x < size; ++x) { px[x + size * y] = LargePixel(kind, p, x, y); }
    }
    std::vector<ASCResult> results = RunASC<int32_t>(N, px, { level }, depths);
    DigestCurves(io, results[0].v, results[0].e);
    DigestCurves(io, results[0].vu, results[0].eu);
}

// The constructor's LogError for N <= 0 (N = 0 here; a negative N makes
// upstream's 1 << N undefined behaviour). N = 1 records extract normally.
ORACLE_CASE("AdaptiveSkeletonClimbing2.invalid")
{
    int32_t N = io.integer(0, 1);
    int const size = (1 << N) + 1;
    std::vector<int64_t> px(static_cast<size_t>(size) * size);
    for (auto& v : px) { v = static_cast<int64_t>(io.integer(-3, 3)); }
    std::vector<ASCResult> results = RunASC<int32_t>(N, px, { 0.5 }, { -1 });
    OutASCResults(io, results);
}
