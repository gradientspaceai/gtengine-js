// Group 26 (imaging): ImageUtility2, ImageUtility3, SurfaceExtractorCubes,
// SurfaceExtractorMC, SurfaceExtractorTetrahedra.
//
// Everything here is integer bookkeeping or IEEE-exact arithmetic (+ - * /,
// sqrt, floor, conversions); no case calls a C math library function other
// than sqrt/sqrtf, so every case is compared bit for bit.
//
// Binary and small-valued images are recorded packed (several pixels per
// recorded double, see GivenPacked) to keep the golden file small; the
// replay unpacks them identically. Callback sequences of the Draw* functions
// are emitted one packed point per output. Bulky outputs (label images,
// rational vertex arrays before duplicate removal, large extractions) are
// emitted as 32-bit FNV-1a digests computed identically on both sides; each
// case comment says which outputs are digests.
#define ORACLE_FAMILY "v26-imaging"
#include "Oracle.h"

#include <Mathematics/ImageUtility2.h>
#include <Mathematics/ImageUtility3.h>
#include <Mathematics/SurfaceExtractorCubes.h>
#include <Mathematics/SurfaceExtractorMC.h>
#include <Mathematics/SurfaceExtractorTetrahedra.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace gte;

namespace
{
    // ---- digests and packing (mirrored in the replay) ----

    // FNV-1a over 32-bit words (the v24 digest).
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
        void integer(int64_t v)
        {
            word(static_cast<uint32_t>(v));
            word(static_cast<uint32_t>(static_cast<uint64_t>(v) >> 32));
        }
        void real(double x)
        {
            uint64_t bits;
            std::memcpy(&bits, &x, sizeof(bits));
            word(static_cast<uint32_t>(bits));
            word(static_cast<uint32_t>(bits >> 32));
        }
    };

    // Values in [0, 2^bits) packed floor(48 / bits) per double, the first
    // value in the lowest bits.
    std::vector<double> Pack(std::vector<int64_t> const& v, int bits)
    {
        int const per = 48 / bits;
        std::vector<double> out;
        for (size_t i = 0; i < v.size(); i += per)
        {
            int64_t chunk = 0;
            for (size_t j = std::min(v.size(), i + per); j > i; --j)
            {
                chunk = (chunk << bits) | v[j - 1];
            }
            out.push_back(static_cast<double>(chunk));
        }
        return out;
    }

    void GivenPacked(oracle::Ctx& io, std::vector<int64_t> const& v, int bits)
    {
        for (double c : Pack(v, bits)) { io.given(c); }
    }

    void OutPacked(oracle::Ctx& io, std::vector<int64_t> const& v, int bits)
    {
        for (double c : Pack(v, bits)) { io.outReal(c); }
    }

    // A pixel or voxel position as one output: 16 bits per coordinate,
    // offset by 2^15.
    void OutPoint(oracle::Ctx& io, int32_t x, int32_t y, int32_t z = 0)
    {
        int64_t p = static_cast<int64_t>(x + 32768)
            + 65536ll * static_cast<int64_t>(y + 32768)
            + 4294967296ll * static_cast<int64_t>(z + 32768);
        io.outInt(p);
    }

    // Records the callback sequence of a Draw* function: fully when short,
    // as count plus digest when long.
    struct PointRecorder
    {
        std::vector<std::array<int32_t, 3>> points;
        void operator()(int32_t x, int32_t y) { points.push_back({ x, y, 0 }); }
        void operator()(int32_t x, int32_t y, int32_t z) { points.push_back({ x, y, z }); }
        void emit(oracle::Ctx& io, bool full) const
        {
            io.outInt(points.size());
            if (full)
            {
                for (auto const& p : points) { OutPoint(io, p[0], p[1], p[2]); }
            }
            else
            {
                Digest d;
                for (auto const& p : points)
                {
                    d.integer(p[0]);
                    d.integer(p[1]);
                    d.integer(p[2]);
                }
                io.outInt(d.h);
            }
        }
    };

    // ---- binary image generators ----

    // Pixel values in {0, 1} (optionally some 2s) on a d0 x d1 grid. Modes:
    // uniform density, rectangles, discs, one solid square (the #443
    // skeleton case), a ring with a hole, a checkerboard or stripes, empty,
    // full. Nothing is recorded here.
    std::vector<int64_t> RawBinary2(oracle::Ctx& io, int d0, int d1, int mode)
    {
        std::vector<int64_t> pix(static_cast<size_t>(d0) * d1, 0);
        auto set = [&](int x, int y) { if (0 <= x && x < d0 && 0 <= y && y < d1) pix[x + d0 * y] = 1; };
        switch (mode)
        {
        case 0:
        {
            double p = io.raw(0.2, 0.8);
            for (auto& v : pix) { v = (io.raw(0.0, 1.0) < p ? 1 : 0); }
            break;
        }
        case 1:
        {
            int n = io.rawInteger(1, 3);
            for (int k = 0; k < n; ++k)
            {
                int x0 = io.rawInteger(0, d0 - 1), y0 = io.rawInteger(0, d1 - 1);
                int w = io.rawInteger(1, d0), h = io.rawInteger(1, d1);
                for (int y = y0; y < y0 + h; ++y) for (int x = x0; x < x0 + w; ++x) set(x, y);
            }
            break;
        }
        case 2:
        {
            int n = io.rawInteger(1, 2);
            for (int k = 0; k < n; ++k)
            {
                int cx = io.rawInteger(0, d0 - 1), cy = io.rawInteger(0, d1 - 1);
                int r2 = io.rawInteger(1, 12);
                for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
                    if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= r2) set(x, y);
            }
            break;
        }
        case 3:
        {
            int s = io.rawInteger(1, std::max(1, std::min(d0, d1) - 2));
            int x0 = io.rawInteger(1, std::max(1, d0 - 1 - s)), y0 = io.rawInteger(1, std::max(1, d1 - 1 - s));
            for (int y = y0; y < y0 + s; ++y) for (int x = x0; x < x0 + s; ++x) set(x, y);
            break;
        }
        case 4:
        {
            for (int y = 1; y < d1 - 1; ++y) for (int x = 1; x < d0 - 1; ++x) set(x, y);
            int t = io.rawInteger(1, 2);
            for (int y = 1 + t; y < d1 - 1 - t; ++y) for (int x = 1 + t; x < d0 - 1 - t; ++x) pix[x + d0 * y] = 0;
            break;
        }
        case 5:
        {
            int kind = io.rawInteger(0, 2);
            for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
            {
                bool on = (kind == 0 ? ((x + y) & 1) == 0 : kind == 1 ? (x % 3) != 2 : ((x + 2 * y) % 4) < 2);
                if (on) set(x, y);
            }
            break;
        }
        case 6:
            break;
        default:
            std::fill(pix.begin(), pix.end(), 1);
            break;
        }
        return pix;
    }

    void ZeroBorder2(std::vector<int64_t>& pix, int d0, int d1)
    {
        for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
            if (x == 0 || y == 0 || x == d0 - 1 || y == d1 - 1) pix[x + d0 * y] = 0;
    }

    // Records d0, d1 and the packed pixels (2 bits each).
    Image2<int32_t> GivenImage2(oracle::Ctx& io, int d0, int d1, std::vector<int64_t> const& pix)
    {
        GivenPacked(io, pix, 2);
        Image2<int32_t> image(d0, d1);
        for (size_t i = 0; i < pix.size(); ++i) { image[i] = static_cast<int32_t>(pix[i]); }
        return image;
    }

    // A binary image; zeroBorder imposes the ImageUtility2 precondition,
    // withTwos turns some 1s into 2s (the '== 1' / '== 0' tests).
    Image2<int32_t> DrawImage2(oracle::Ctx& io, int dmin, int dmax, bool zeroBorder, bool withTwos)
    {
        int d0 = io.integer(dmin, dmax);
        int d1 = io.integer(dmin, dmax);
        std::vector<int64_t> pix = RawBinary2(io, d0, d1, io.rawInteger(0, 7));
        if (withTwos && io.rawInteger(0, 3) == 0)
        {
            for (auto& v : pix) { if (v == 1 && io.rawInteger(0, 3) == 0) v = 2; }
        }
        if (zeroBorder) { ZeroBorder2(pix, d0, d1); }
        return GivenImage2(io, d0, d1, pix);
    }

    void OutImage2(oracle::Ctx& io, Image2<int32_t> const& image, int bits)
    {
        std::vector<int64_t> v(image.GetNumPixels());
        for (size_t i = 0; i < v.size(); ++i) { v[i] = image[i]; }
        OutPacked(io, v, bits);
    }
}

// ---- ImageUtility2: components and morphology ----

// Outputs: the component count (components.size(), 0 when there are none),
// each component's size and indices (8 bits packed; upstream fills them in
// increasing index order), and the relabelled image (8 bits packed). kind 2
// passes 0..8 caller-specified 1-D offsets from the 3x3 block (duplicates
// and the 0 offset allowed, an empty list included).
ORACLE_CASE("ImageUtility2.getComponents")
{
    int kind = io.integer(0, 2);
    Image2<int32_t> image = DrawImage2(io, 3, 10, true, false);
    std::vector<std::vector<size_t>> components;
    if (kind == 0)
    {
        ImageUtility2::GetComponents<4>(image, components);
    }
    else if (kind == 1)
    {
        ImageUtility2::GetComponents<8>(image, components);
    }
    else
    {
        int n = io.integer(0, 8);
        std::vector<int32_t> nbrs;
        for (int k = 0; k < n; ++k)
        {
            int dx = io.integer(-1, 1);
            int dy = io.integer(-1, 1);
            nbrs.push_back(dx + image.GetDimension(0) * dy);
        }
        ImageUtility2::GetComponents(image, nbrs.size(), nbrs.data(), components);
    }
    io.outInt(components.size());
    for (size_t k = 1; k < components.size(); ++k)
    {
        io.outInt(components[k].size());
        std::vector<int64_t> idx(components[k].begin(), components[k].end());
        OutPacked(io, idx, 8);
    }
    OutImage2(io, image, 8);
}

namespace
{
    // kind 0: N = 4, 1: N = 8, 2: 0..6 caller-specified offsets in
    // [-2,2]^2 (0 offsets is the LogAssert path, a throw record).
    std::vector<std::array<int32_t, 2>> DrawOffsets2(oracle::Ctx& io, int kind)
    {
        std::vector<std::array<int32_t, 2>> nbrs;
        if (kind == 2)
        {
            int n = io.integer(0, 6);
            for (int k = 0; k < n; ++k)
            {
                int dx = io.integer(-2, 2);
                int dy = io.integer(-2, 2);
                nbrs.push_back({ dx, dy });
            }
        }
        return nbrs;
    }

    void Morph2(oracle::Ctx& io, int op)
    {
        int kind = io.integer(0, 2);
        Image2<int32_t> image = DrawImage2(io, 1, 9, false, true);
        bool zeroExterior = (op != 0 ? io.boolean() : false);
        auto nbrs = DrawOffsets2(io, kind);
        Image2<int32_t> out;
        switch (op)
        {
        case 0:
            if (kind == 0) ImageUtility2::Dilate<4>(image, out);
            else if (kind == 1) ImageUtility2::Dilate<8>(image, out);
            else ImageUtility2::Dilate(image, nbrs.size(), nbrs.data(), out);
            break;
        case 1:
            if (kind == 0) ImageUtility2::Erode<4>(image, zeroExterior, out);
            else if (kind == 1) ImageUtility2::Erode<8>(image, zeroExterior, out);
            else ImageUtility2::Erode(image, zeroExterior, nbrs.size(), nbrs.data(), out);
            break;
        case 2:
            if (kind == 0) ImageUtility2::Open<4>(image, zeroExterior, out);
            else if (kind == 1) ImageUtility2::Open<8>(image, zeroExterior, out);
            else ImageUtility2::Open(image, zeroExterior, nbrs.size(), nbrs.data(), out);
            break;
        default:
            if (kind == 0) ImageUtility2::Close<4>(image, zeroExterior, out);
            else if (kind == 1) ImageUtility2::Close<8>(image, zeroExterior, out);
            else ImageUtility2::Close(image, zeroExterior, nbrs.size(), nbrs.data(), out);
            break;
        }
        io.outInt(out.GetDimension(0));
        io.outInt(out.GetDimension(1));
        OutImage2(io, out, 2);
    }

    // Replica of ExtractBoundary's control flow with a step cap: the walk
    // stops only when it returns to its first pixel, which a start inside a
    // component need not do. Used to reject such starts.
    bool BoundaryTerminates(Image2<int32_t> const& image, int32_t x, int32_t y, size_t cap)
    {
        int32_t const d0 = image.GetDimension(0);
        std::vector<int32_t> pix(image.GetNumPixels());
        for (size_t k = 0; k < pix.size(); ++k) { pix[k] = image[k]; }
        size_t i = static_cast<size_t>(x) + static_cast<size_t>(d0) * y;
        for (; i < pix.size(); ++i) { if (pix[i]) break; }
        if (i == pix.size()) return true;
        std::array<int32_t, 8> const dx = { -1, 0, +1, +1, +1, 0, -1, -1 };
        std::array<int32_t, 8> const dy = { -1, -1, -1, 0, +1, +1, +1, 0 };
        int32_t x0 = static_cast<int32_t>(i % d0), y0 = static_cast<int32_t>(i / d0);
        int32_t cx = x0, cy = y0, nx = x0 - 1, ny = y0, dir = 7;
        pix[i] = 2;
        for (size_t step = 0; step < cap; ++step)
        {
            int32_t j = 0, nbr = dir;
            for (j = 0, nbr = dir; j < 8; ++j, nbr = (nbr + 1) % 8)
            {
                nx = cx + dx[nbr];
                ny = cy + dy[nbr];
                if (pix[nx + d0 * ny]) break;
            }
            if (j == 8) return true;
            if (nx == x0 && ny == y0) return true;
            pix[nx + d0 * ny] = 2;
            cx = nx;
            cy = ny;
            dir = (j + 5 + dir) % 8;
        }
        return false;
    }
}

ORACLE_CASE("ImageUtility2.dilate") { Morph2(io, 0); }
ORACLE_CASE("ImageUtility2.erode") { Morph2(io, 1); }
ORACLE_CASE("ImageUtility2.open") { Morph2(io, 2); }
ORACLE_CASE("ImageUtility2.close") { Morph2(io, 3); }

// One to three walks on the same image (the second sees the pixels the
// first marked 2). Each start is redrawn (at most 32 times) until the
// replica terminates; the fallback is the last pixel, which lies on the
// zero border, so the walk finds nothing. Outputs per walk: the return
// value, the boundary length and indices (8 bits packed); then the image.
ORACLE_CASE("ImageUtility2.extractBoundary")
{
    Image2<int32_t> image = DrawImage2(io, 3, 10, true, false);
    int32_t const d0 = image.GetDimension(0), d1 = image.GetDimension(1);
    int calls = io.integer(1, 3);
    for (int c = 0; c < calls; ++c)
    {
        int32_t x = d0 - 1, y = d1 - 1;
        for (int attempt = 0; attempt < 32; ++attempt)
        {
            int32_t tx = io.rawInteger(0, d0 - 1), ty = io.rawInteger(0, d1 - 1);
            if (BoundaryTerminates(image, tx, ty, 4 * image.GetNumPixels() + 16))
            {
                x = tx;
                y = ty;
                break;
            }
        }
        io.given(x);
        io.given(y);
        std::vector<size_t> boundary;
        bool ok = ImageUtility2::ExtractBoundary(x, y, image, boundary);
        io.outBool(ok);
        io.outInt(boundary.size());
        std::vector<int64_t> idx(boundary.begin(), boundary.end());
        OutPacked(io, idx, 8);
    }
    OutImage2(io, image, 2);
}

// ---- ImageUtility2: fills and distance transforms ----

// Colors 0..3; the background color is the seed's color on three records
// of four (a real fill) and the foreground differs from it (equal colors
// never terminate). Seeds one outside the image are the early return.
// Output: the image (2 bits packed).
ORACLE_CASE("ImageUtility2.floodFill4.int")
{
    int32_t d0 = io.integer(1, 9);
    int32_t d1 = io.integer(1, 9);
    std::vector<int64_t> pix = RawBinary2(io, d0, d1, io.rawInteger(0, 7));
    for (auto& v : pix) { if (io.rawInteger(0, 4) == 0) v = io.rawInteger(2, 3); }
    Image2<int32_t> image = GivenImage2(io, d0, d1, pix);
    int32_t x = io.integer(-1, d0);
    int32_t y = io.integer(-1, d1);
    bool inside = (0 <= x && x < d0 && 0 <= y && y < d1);
    int32_t back = (inside && io.rawInteger(0, 3) != 0 ? image(x, y) : io.rawInteger(0, 3));
    io.given(back);
    int32_t fore = (back + io.integer(1, 3)) % 4;
    ImageUtility2::FloodFill4<int32_t>(image, x, y, fore, back);
    OutImage2(io, image, 2);
}

namespace
{
    double const gPalette[5] = { 0.0, 1.0, -0.0, 2.5, std::numeric_limits<double>::quiet_NaN() };
}

// Double pixels from {0, 1, -0, 2.5, NaN}: -0 == +0 fills across signed
// zeros, a NaN background fills only the seed. The foreground is drawn
// until it does not compare equal to the background. Output: every pixel.
ORACLE_CASE("ImageUtility2.floodFill4.double")
{
    int32_t d0 = io.integer(1, 8);
    int32_t d1 = io.integer(1, 8);
    std::vector<int64_t> pix = RawBinary2(io, d0, d1, io.rawInteger(0, 7));
    Image2<double> image(d0, d1);
    for (size_t i = 0; i < pix.size(); ++i)
    {
        int p = (io.rawInteger(0, 3) == 0 ? io.rawInteger(0, 4) : static_cast<int>(pix[i]));
        image[i] = io.given(gPalette[p]);
    }
    int32_t x = io.integer(-1, d0);
    int32_t y = io.integer(-1, d1);
    bool inside = (0 <= x && x < d0 && 0 <= y && y < d1);
    double back = (inside && io.rawInteger(0, 3) != 0 ? image(x, y) : gPalette[io.rawInteger(0, 4)]);
    io.given(back);
    double fore = gPalette[io.rawInteger(0, 4)];
    for (int k = 0; k < 64 && fore == back; ++k) { fore = gPalette[io.rawInteger(0, 4)]; }
    if (fore == back) { fore = (back == 1.0 ? 2.5 : 1.0); }
    io.given(fore);
    ImageUtility2::FloodFill4<double>(image, x, y, fore, back);
    for (size_t i = 0; i < image.GetNumPixels(); ++i) { io.outReal(image[i]); }
}

// The border is zero on half of the records (the precondition); otherwise
// border pixels keep their value 1 and act as distance 1. Outputs: maximum
// distance, its pixel, the image (4 bits packed).
ORACLE_CASE("ImageUtility2.getL1Distance")
{
    int32_t d0 = io.integer(1, 12);
    int32_t d1 = io.integer(1, 12);
    std::vector<int64_t> pix = RawBinary2(io, d0, d1, io.rawInteger(0, 7));
    if (io.boolean()) { ZeroBorder2(pix, d0, d1); }
    Image2<int32_t> image = GivenImage2(io, d0, d1, pix);
    int32_t maxDistance = -7, xMax = -7, yMax = -7;
    ImageUtility2::GetL1Distance(image, maxDistance, xMax, yMax);
    io.outInt(maxDistance);
    io.outInt(xMax);
    io.outInt(yMax);
    OutImage2(io, image, 4);
}

// Any binary image (all-foreground images leave every distance at
// INT32_MAX, sqrtf(2^31)). Outputs: maximum (float), its pixel, the float
// transform. sqrtf is IEEE-exact.
ORACLE_CASE("ImageUtility2.getL2Distance")
{
    int32_t d0 = io.integer(1, 12);
    int32_t d1 = io.integer(1, 12);
    std::vector<int64_t> pix = RawBinary2(io, d0, d1, io.rawInteger(0, 7));
    Image2<int32_t> image = GivenImage2(io, d0, d1, pix);
    float maxDistance = -7.0f;
    int32_t xMax = -7, yMax = -7;
    Image2<float> transform(d0, d1);
    ImageUtility2::GetL2Distance(image, maxDistance, xMax, yMax, transform);
    io.outReal(maxDistance);
    io.outInt(xMax);
    io.outInt(yMax);
    for (size_t i = 0; i < transform.GetNumPixels(); ++i) { io.outReal(transform[i]); }
}

// Large images with a few background pixels, so that distances pass the
// K2..K5 thresholds (13, 31, 49, 72) of the extended neighbourhood checks
// while staying below 100 (the documented exactness range; d <= 70 keeps
// every distance below 99). Recorded: dimensions, whether the border is
// zero, the background pixel indices. Outputs: maximum, its pixel, a
// digest of the transform.
ORACLE_CASE("ImageUtility2.getL2Distance.large")
{
    int32_t d0 = io.integer(20, 70);
    int32_t d1 = io.integer(20, 70);
    Image2<int32_t> image(d0, d1);
    for (size_t i = 0; i < image.GetNumPixels(); ++i) { image[i] = 1; }
    if (io.boolean())
    {
        for (int32_t y = 0; y < d1; ++y) for (int32_t x = 0; x < d0; ++x)
            if (x == 0 || y == 0 || x == d0 - 1 || y == d1 - 1) image(x, y) = 0;
    }
    int k = io.integer(1, (io.rawInteger(0, 1) == 0 ? 4 : 40));
    for (int j = 0; j < k; ++j) { image[io.integer(0, d0 * d1 - 1)] = 0; }
    float maxDistance = -7.0f;
    int32_t xMax = -7, yMax = -7;
    Image2<float> transform(d0, d1);
    ImageUtility2::GetL2Distance(image, maxDistance, xMax, yMax, transform);
    io.outReal(maxDistance);
    io.outInt(xMax);
    io.outInt(yMax);
    Digest d;
    for (size_t i = 0; i < transform.GetNumPixels(); ++i) { d.real(transform[i]); }
    io.outInt(d.h);
}

// Zero-border binary images, including solid squares of side 1..6 (#443:
// even sides vanish, preserved) and rings. Output: the skeleton (1 bit).
ORACLE_CASE("ImageUtility2.getSkeleton")
{
    Image2<int32_t> image = DrawImage2(io, 3, 12, true, false);
    ImageUtility2::GetSkeleton(image);
    OutImage2(io, image, 1);
}

// ---- ImageUtility2: Draw* through a recording callback ----
//
// Each case emits the callback sequence in call order: the count, then one
// packed point per call when the sequence is short (full), or a digest of
// it (the large modes). The mode is recorded.

ORACLE_CASE("ImageUtility2.drawThickPixel")
{
    int32_t x = io.integer(-50, 50);
    int32_t y = io.integer(-50, 50);
    int32_t thick = io.integer(-2, 4);
    PointRecorder rec;
    ImageUtility2::DrawThickPixel(x, y, thick, std::ref(rec));
    rec.emit(io, true);
}

// Modes: random short lines, axis-aligned, exact diagonals, a single point,
// shallow/steep lines with |dx| = |dy| +- 1, long lines (digest).
ORACLE_CASE("ImageUtility2.drawLine")
{
    int mode = io.integer(0, 5);
    int32_t r = (mode == 5 ? 5000 : 20);
    int32_t x0 = io.integer(-r, r);
    int32_t y0 = io.integer(-r, r);
    int32_t x1 = io.integer(-r, r);
    int32_t y1 = io.integer(-r, r);
    int32_t len = io.integer(-12, 12);
    int32_t s = io.integer(-1, 1);
    switch (mode)
    {
    case 1: if (s >= 0) { y1 = y0; } else { x1 = x0; } break;
    case 2: x1 = x0 + len; y1 = y0 + (s >= 0 ? len : -len); break;
    case 3: x1 = x0; y1 = y0; break;
    case 4: x1 = x0 + len; y1 = y0 + (len >= 0 ? len + s : len - s); break;
    default: break;
    }
    PointRecorder rec;
    ImageUtility2::DrawLine(x0, y0, x1, y1, std::ref(rec));
    rec.emit(io, mode != 5);
}

// Radii -2..12; mode 3 draws radii up to 1000 (outline) or 150 (solid),
// emitted as digests.
ORACLE_CASE("ImageUtility2.drawCircle")
{
    int mode = io.integer(0, 3);
    int32_t xc = io.integer(-20, 20);
    int32_t yc = io.integer(-20, 20);
    bool solid = io.boolean();
    int32_t radius = (mode == 3 ? io.integer(0, solid ? 150 : 1000) : io.integer(-2, 12));
    PointRecorder rec;
    ImageUtility2::DrawCircle(xc, yc, radius, solid, std::ref(rec));
    rec.emit(io, mode != 3);
}

ORACLE_CASE("ImageUtility2.drawRectangle")
{
    int32_t xMin = io.integer(-10, 10);
    int32_t yMin = io.integer(-10, 10);
    int32_t xMax = io.integer(-10, 10);
    int32_t yMax = io.integer(-10, 10);
    bool solid = io.boolean();
    PointRecorder rec;
    ImageUtility2::DrawRectangle(xMin, yMin, xMax, yMax, solid, std::ref(rec));
    rec.emit(io, true);
}

// Extents in [-2, 14] (negative extents return early or draw one arc),
// never both zero (#443, the deviation case below); mode 3 draws extents up
// to 100 (digest; larger extents overflow int32 in 4*a^2*(1-y)).
ORACLE_CASE("ImageUtility2.drawEllipse")
{
    int mode = io.integer(0, 3);
    int32_t xc = io.integer(-20, 20);
    int32_t yc = io.integer(-20, 20);
    int32_t hi = (mode == 3 ? 100 : 14);
    int32_t a = io.integer(-2, hi);
    int32_t b = io.integer((a == 0 ? 1 : -2), hi);
    if (mode == 1) { b = a + (b & 1); if (b == 0) b = 1; }
    PointRecorder rec;
    ImageUtility2::DrawEllipse(xc, yc, a, b, std::ref(rec));
    rec.emit(io, mode != 3);
}

namespace
{
    struct ThrowingRecorder
    {
        std::vector<std::array<int32_t, 3>>* points;
        void operator()(int32_t x, int32_t y)
        {
            if (points->size() >= 1000) { throw std::runtime_error("DrawEllipse does not terminate"); }
            points->push_back({ x, y, 0 });
        }
    };
}

// #443 (fixed in the port): with both extents zero upstream never
// terminates (0 <= 0 forever while y decreases). The C++ callback throws
// after 1000 calls, so every record is a throw record; the port visits the
// center once and returns, so every record deviates.
ORACLE_CASE("ImageUtility2.drawEllipse.zeroExtents")
{
    int32_t xc = io.integer(-20, 20);
    int32_t yc = io.integer(-20, 20);
    std::vector<std::array<int32_t, 3>> points;
    ImageUtility2::DrawEllipse(xc, yc, 0, 0, ThrowingRecorder{ &points });
    io.outInt(points.size());
}

// DrawFloodFill4 with callbacks on a C++ int image. Outputs: the sequence
// of setCallback calls (count, points), the number of getCallback calls,
// the final image (2 bits packed).
ORACLE_CASE("ImageUtility2.drawFloodFill4")
{
    int32_t d0 = io.integer(1, 9);
    int32_t d1 = io.integer(1, 9);
    std::vector<int64_t> pix = RawBinary2(io, d0, d1, io.rawInteger(0, 7));
    for (auto& v : pix) { if (io.rawInteger(0, 4) == 0) v = io.rawInteger(2, 3); }
    Image2<int32_t> image = GivenImage2(io, d0, d1, pix);
    int32_t x = io.integer(-1, d0);
    int32_t y = io.integer(-1, d1);
    bool inside = (0 <= x && x < d0 && 0 <= y && y < d1);
    int32_t back = (inside && io.rawInteger(0, 3) != 0 ? image(x, y) : io.rawInteger(0, 3));
    io.given(back);
    int32_t fore = (back + io.integer(1, 3)) % 4;
    PointRecorder sets;
    int64_t numGets = 0;
    ImageUtility2::DrawFloodFill4<int32_t>(x, y, d0, d1, fore, back,
        [&](int32_t px, int32_t py, int32_t value) { sets(px, py); image(px, py) = value; },
        [&](int32_t px, int32_t py) { ++numGets; return image(px, py); });
    sets.emit(io, true);
    io.outInt(numGets);
    OutImage2(io, image, 2);
}

// ---- ImageUtility3 ----

namespace
{
    // Voxel values in {0, 1} on a d0 x d1 x d2 grid: uniform density, boxes,
    // balls, a solid block, a checkerboard, empty, full. Nothing recorded.
    std::vector<int64_t> RawBinary3(oracle::Ctx& io, int d0, int d1, int d2, int mode)
    {
        std::vector<int64_t> pix(static_cast<size_t>(d0) * d1 * d2, 0);
        auto at = [&](int x, int y, int z) -> int64_t& { return pix[x + d0 * (y + d1 * z)]; };
        auto inside = [&](int x, int y, int z) { return 0 <= x && x < d0 && 0 <= y && y < d1 && 0 <= z && z < d2; };
        switch (mode)
        {
        case 0:
        {
            double p = io.raw(0.2, 0.8);
            for (auto& v : pix) { v = (io.raw(0.0, 1.0) < p ? 1 : 0); }
            break;
        }
        case 1:
        {
            int n = io.rawInteger(1, 3);
            for (int k = 0; k < n; ++k)
            {
                int x0 = io.rawInteger(0, d0 - 1), y0 = io.rawInteger(0, d1 - 1), z0 = io.rawInteger(0, d2 - 1);
                int w = io.rawInteger(1, d0), h = io.rawInteger(1, d1), t = io.rawInteger(1, d2);
                for (int z = z0; z < z0 + t; ++z) for (int y = y0; y < y0 + h; ++y) for (int x = x0; x < x0 + w; ++x)
                    if (inside(x, y, z)) at(x, y, z) = 1;
            }
            break;
        }
        case 2:
        {
            int cx = io.rawInteger(0, d0 - 1), cy = io.rawInteger(0, d1 - 1), cz = io.rawInteger(0, d2 - 1);
            int r2 = io.rawInteger(1, 8);
            for (int z = 0; z < d2; ++z) for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
                if ((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz) <= r2) at(x, y, z) = 1;
            break;
        }
        case 3:
        {
            for (int z = 0; z < d2; ++z) for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
                if (((x + y + z) & 1) == 0) at(x, y, z) = 1;
            break;
        }
        case 4:
            break;
        case 5:
            std::fill(pix.begin(), pix.end(), 1);
            break;
        default:
        {
            // Sparse: a few isolated voxels and short rods.
            int n = io.rawInteger(1, 6);
            for (int k = 0; k < n; ++k)
            {
                int x = io.rawInteger(0, d0 - 1), y = io.rawInteger(0, d1 - 1), z = io.rawInteger(0, d2 - 1);
                int axis = io.rawInteger(0, 2), len = io.rawInteger(1, 3);
                for (int j = 0; j < len; ++j)
                {
                    int px = x + (axis == 0 ? j : 0), py = y + (axis == 1 ? j : 0), pz = z + (axis == 2 ? j : 0);
                    if (inside(px, py, pz)) at(px, py, pz) = 1;
                }
            }
            break;
        }
        }
        return pix;
    }

    // borderMode 0: none, 1: zero on all six faces (the header's
    // precondition), 2: zero on the x = 0 face only (where upstream's Dilate
    // is sound, #129).
    Image3<int32_t> DrawImage3(oracle::Ctx& io, int dmin, int dmax, int borderMode, bool withTwos)
    {
        int d0 = io.integer(dmin, dmax);
        int d1 = io.integer(dmin, dmax);
        int d2 = io.integer(dmin, dmax);
        std::vector<int64_t> pix = RawBinary3(io, d0, d1, d2, io.rawInteger(0, 6));
        if (withTwos && io.rawInteger(0, 3) == 0)
        {
            for (auto& v : pix) { if (v == 1 && io.rawInteger(0, 3) == 0) v = 2; }
        }
        for (int z = 0; z < d2; ++z) for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
        {
            bool onBorder = (x == 0 || y == 0 || z == 0 || x == d0 - 1 || y == d1 - 1 || z == d2 - 1);
            if ((borderMode == 1 && onBorder) || (borderMode == 2 && x == 0)) pix[x + d0 * (y + d1 * z)] = 0;
        }
        GivenPacked(io, pix, 2);
        Image3<int32_t> image(d0, d1, d2);
        for (size_t i = 0; i < pix.size(); ++i) { image[i] = static_cast<int32_t>(pix[i]); }
        return image;
    }

    void OutImage3(oracle::Ctx& io, Image3<int32_t> const& image, int bits)
    {
        std::vector<int64_t> v(image.GetNumPixels());
        for (size_t i = 0; i < v.size(); ++i) { v[i] = image[i]; }
        OutPacked(io, v, bits);
    }
}

// kind 0/1/2: N = 6/18/26; kind 3: 0..8 caller-specified 1-D offsets from
// the 3x3x3 block (0 offsets is the LogAssert path, a throw record).
// Outputs as for ImageUtility2.getComponents.
ORACLE_CASE("ImageUtility3.getComponents")
{
    int kind = io.integer(0, 3);
    Image3<int32_t> image = DrawImage3(io, 3, 6, 1, false);
    std::vector<std::vector<size_t>> components;
    if (kind == 0) { ImageUtility3::GetComponents<6>(image, components); }
    else if (kind == 1) { ImageUtility3::GetComponents<18>(image, components); }
    else if (kind == 2) { ImageUtility3::GetComponents<26>(image, components); }
    else
    {
        int n = io.integer(0, 8);
        std::vector<int32_t> nbrs;
        for (int k = 0; k < n; ++k)
        {
            int dx = io.integer(-1, 1);
            int dy = io.integer(-1, 1);
            int dz = io.integer(-1, 1);
            nbrs.push_back(dx + image.GetDimension(0) * (dy + image.GetDimension(1) * dz));
        }
        ImageUtility3::GetComponents(image, nbrs.size(), nbrs.data(), components);
    }
    io.outInt(components.size());
    for (size_t k = 1; k < components.size(); ++k)
    {
        io.outInt(components[k].size());
        std::vector<int64_t> idx(components[k].begin(), components[k].end());
        OutPacked(io, idx, 8);
    }
    OutImage3(io, image, 8);
}

namespace
{
    template <size_t N>
    std::vector<std::array<int32_t, 3>> Neighborhood3(Image3<int32_t> const& image)
    {
        std::array<std::array<int32_t, 3>, N> a{};
        image.GetNeighborhood(a);
        return std::vector<std::array<int32_t, 3>>(a.begin(), a.end());
    }

    // op 0..3: Dilate, Erode, Open, Close. kind 0/1/2: N = 6/18/26, kind 3:
    // 0..6 caller-specified offsets in [-2,2]^3. Close<N> does not compile
    // upstream (#129), so kinds 0..2 of Close call the offset-list Close with
    // the N-neighbourhood, which is what the port's close6/18/26 compute.
    void Morph3(oracle::Ctx& io, int op, int borderMode)
    {
        int kind = io.integer(0, 3);
        Image3<int32_t> image = DrawImage3(io, 1, 6, borderMode, true);
        bool zeroExterior = (op != 0 ? io.boolean() : false);
        std::vector<std::array<int32_t, 3>> nbrs;
        if (kind == 3)
        {
            int n = io.integer(0, 6);
            for (int k = 0; k < n; ++k)
            {
                int dx = io.integer(-2, 2);
                int dy = io.integer(-2, 2);
                int dz = io.integer(-2, 2);
                nbrs.push_back({ dx, dy, dz });
            }
        }
        else if (op == 3)
        {
            nbrs = (kind == 0 ? Neighborhood3<6>(image) : kind == 1 ? Neighborhood3<18>(image) : Neighborhood3<26>(image));
        }
        Image3<int32_t> out;
        switch (op)
        {
        case 0:
            if (kind == 0) ImageUtility3::Dilate<6>(image, out);
            else if (kind == 1) ImageUtility3::Dilate<18>(image, out);
            else if (kind == 2) ImageUtility3::Dilate<26>(image, out);
            else ImageUtility3::Dilate(image, nbrs.size(), nbrs.data(), out);
            break;
        case 1:
            if (kind == 0) ImageUtility3::Erode<6>(image, zeroExterior, out);
            else if (kind == 1) ImageUtility3::Erode<18>(image, zeroExterior, out);
            else if (kind == 2) ImageUtility3::Erode<26>(image, zeroExterior, out);
            else ImageUtility3::Erode(image, zeroExterior, nbrs.size(), nbrs.data(), out);
            break;
        case 2:
            if (kind == 0) ImageUtility3::Open<6>(image, zeroExterior, out);
            else if (kind == 1) ImageUtility3::Open<18>(image, zeroExterior, out);
            else if (kind == 2) ImageUtility3::Open<26>(image, zeroExterior, out);
            else ImageUtility3::Open(image, zeroExterior, nbrs.size(), nbrs.data(), out);
            break;
        default:
            ImageUtility3::Close(image, zeroExterior, nbrs.size(), nbrs.data(), out);
            break;
        }
        io.outInt(out.GetDimension(0));
        io.outInt(out.GetDimension(1));
        io.outInt(out.GetDimension(2));
        OutImage3(io, out, 2);
    }
}

// Main cases: the x = 0 face is zero (upstream's Dilate skips x = 0 sources,
// #129, fixed in the port); the rest of the image is unrestricted, with 2s
// on some records.
ORACLE_CASE("ImageUtility3.dilate") { Morph3(io, 0, 2); }
ORACLE_CASE("ImageUtility3.erode") { Morph3(io, 1, io.index() % 2 == 0 ? 0 : 2); }
ORACLE_CASE("ImageUtility3.open") { Morph3(io, 2, 2); }
ORACLE_CASE("ImageUtility3.close") { Morph3(io, 3, 2); }

// #129 deviation: images with 1s on the x = 0 face, through Dilate, Open
// and Close.
ORACLE_CASE("ImageUtility3.morphology.xmin")
{
    int op = io.integer(0, 2);
    Morph3(io, op == 0 ? 0 : op == 1 ? 2 : 3, 0);
}

// Any binary image. Output: the result (1 bit packed).
ORACLE_CASE("ImageUtility3.computeCDConvex")
{
    Image3<int32_t> image = DrawImage3(io, 1, 6, 0, false);
    ImageUtility3::ComputeCDConvex(image);
    OutImage3(io, image, 1);
}

ORACLE_CASE("ImageUtility3.floodFill6.int")
{
    int32_t d0 = io.integer(1, 5);
    int32_t d1 = io.integer(1, 5);
    int32_t d2 = io.integer(1, 5);
    std::vector<int64_t> pix = RawBinary3(io, d0, d1, d2, io.rawInteger(0, 6));
    for (auto& v : pix) { if (io.rawInteger(0, 4) == 0) v = io.rawInteger(2, 3); }
    GivenPacked(io, pix, 2);
    Image3<int32_t> image(d0, d1, d2);
    for (size_t i = 0; i < pix.size(); ++i) { image[i] = static_cast<int32_t>(pix[i]); }
    int32_t x = io.integer(-1, d0);
    int32_t y = io.integer(-1, d1);
    int32_t z = io.integer(-1, d2);
    bool inside = (0 <= x && x < d0 && 0 <= y && y < d1 && 0 <= z && z < d2);
    int32_t back = (inside && io.rawInteger(0, 3) != 0 ? image(x, y, z) : io.rawInteger(0, 3));
    io.given(back);
    int32_t fore = (back + io.integer(1, 3)) % 4;
    ImageUtility3::FloodFill6<int32_t>(image, x, y, z, fore, back);
    OutImage3(io, image, 2);
}

ORACLE_CASE("ImageUtility3.floodFill6.double")
{
    int32_t d0 = io.integer(1, 4);
    int32_t d1 = io.integer(1, 4);
    int32_t d2 = io.integer(1, 4);
    std::vector<int64_t> pix = RawBinary3(io, d0, d1, d2, io.rawInteger(0, 6));
    Image3<double> image(d0, d1, d2);
    for (size_t i = 0; i < pix.size(); ++i)
    {
        int p = (io.rawInteger(0, 3) == 0 ? io.rawInteger(0, 4) : static_cast<int>(pix[i]));
        image[i] = io.given(gPalette[p]);
    }
    int32_t x = io.integer(-1, d0);
    int32_t y = io.integer(-1, d1);
    int32_t z = io.integer(-1, d2);
    bool inside = (0 <= x && x < d0 && 0 <= y && y < d1 && 0 <= z && z < d2);
    double back = (inside && io.rawInteger(0, 3) != 0 ? image(x, y, z) : gPalette[io.rawInteger(0, 4)]);
    io.given(back);
    double fore = gPalette[io.rawInteger(0, 4)];
    for (int k = 0; k < 64 && fore == back; ++k) { fore = gPalette[io.rawInteger(0, 4)]; }
    if (fore == back) { fore = (back == 1.0 ? 2.5 : 1.0); }
    io.given(fore);
    ImageUtility3::FloodFill6<double>(image, x, y, z, fore, back);
    for (size_t i = 0; i < image.GetNumPixels(); ++i) { io.outReal(image[i]); }
}

// Modes: random, axis-aligned, exact diagonals, a single point, ties between
// the two largest |d| components (the strict '>' selects the axis), long
// lines (digest).
ORACLE_CASE("ImageUtility3.drawLine")
{
    int mode = io.integer(0, 5);
    int32_t r = (mode == 5 ? 3000 : 12);
    int32_t x0 = io.integer(-r, r);
    int32_t y0 = io.integer(-r, r);
    int32_t z0 = io.integer(-r, r);
    int32_t x1 = io.integer(-r, r);
    int32_t y1 = io.integer(-r, r);
    int32_t z1 = io.integer(-r, r);
    int32_t len = io.integer(-10, 10);
    int32_t s = io.integer(0, 5);
    switch (mode)
    {
    case 1: x1 = (s % 3 == 0 ? x1 : x0); y1 = (s % 3 == 1 ? y1 : y0); z1 = (s % 3 == 2 ? z1 : z0); break;
    case 2: x1 = x0 + len; y1 = y0 + ((s & 1) ? len : -len); z1 = z0 + ((s & 2) ? len : -len); break;
    case 3: x1 = x0; y1 = y0; z1 = z0; break;
    case 4:
    {
        // Two components of equal magnitude |len|, the third smaller.
        int32_t small = (len == 0 ? 0 : (len > 0 ? len - 1 : len + 1));
        std::array<int32_t, 3> d = { len, (s & 1) ? len : -len, small };
        std::rotate(d.begin(), d.begin() + (s % 3), d.end());
        x1 = x0 + d[0]; y1 = y0 + d[1]; z1 = z0 + d[2];
        break;
    }
    default: break;
    }
    PointRecorder rec;
    ImageUtility3::DrawLine(x0, y0, z0, x1, y1, z1, std::ref(rec));
    rec.emit(io, mode != 5);
}

// ---- surface extractors: shared generators and emission ----

namespace
{
    // Raw integer voxel values. Modes: uniform in [lo, hi], small lattice,
    // linear (no saddle faces), a quadric, binary, a hyperbolic product
    // (saddle faces).
    std::vector<int64_t> RawVoxels(oracle::Ctx& io, int d0, int d1, int d2, int mode, int64_t lo, int64_t hi)
    {
        std::vector<int64_t> v(static_cast<size_t>(d0) * d1 * d2);
        int64_t a = io.rawInteger(-3, 3), b = io.rawInteger(-3, 3), c = io.rawInteger(-3, 3);
        int64_t e = io.rawInteger(-5, 5);
        int64_t cx = io.rawInteger(0, 2 * d0 - 2), cy = io.rawInteger(0, 2 * d1 - 2), cz = io.rawInteger(0, 2 * d2 - 2);
        int64_t r = io.rawInteger(0, 12);
        for (int z = 0; z < d2; ++z) for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
        {
            int64_t value = 0;
            switch (mode)
            {
            case 0: value = lo + static_cast<int64_t>(io.raw(0.0, 1.0) * static_cast<double>(hi - lo + 1)); break;
            case 1: value = io.rawInteger(-3, 3); break;
            case 2: value = a * x + b * y + c * z + e; break;
            case 3: value = (2 * x - cx) * (2 * x - cx) + (2 * y - cy) * (2 * y - cy) + (2 * z - cz) * (2 * z - cz) - r; break;
            case 4: value = io.rawInteger(0, 1); break;
            default: value = (2 * x - cx) * (2 * y - cy) + a * (2 * z - cz) + e; break;
            }
            v[x + d0 * (y + d1 * z)] = std::min(hi, std::max(lo, value));
        }
        return v;
    }

    // Faces of the voxel grid whose corners alternate in sign about the
    // Cubes level (shifted values 2v - 2L - 1, never zero), split by the
    // sign of SurfaceExtractorCubes' face determinant f00*f11 - f01*f10.
    struct SaddleCount { int nonzero = 0, zero = 0; };

    SaddleCount CubesSaddles(std::vector<int64_t> const& v, int d0, int d1, int d2, int64_t level)
    {
        SaddleCount count;
        auto s = [&](int x, int y, int z) { return 2 * v[x + d0 * (y + d1 * z)] - (2 * level + 1); };
        auto face = [&](int64_t f00, int64_t f10, int64_t f11, int64_t f01)
        {
            if (f00 * f10 < 0 && f10 * f11 < 0 && f11 * f01 < 0)
            {
                if (f00 * f11 - f01 * f10 != 0) ++count.nonzero; else ++count.zero;
            }
        };
        for (int z = 0; z < d2; ++z) for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
        {
            if (y + 1 < d1 && z + 1 < d2) face(s(x, y, z), s(x, y + 1, z), s(x, y + 1, z + 1), s(x, y, z + 1));
            if (x + 1 < d0 && z + 1 < d2) face(s(x, y, z), s(x + 1, y, z), s(x + 1, y, z + 1), s(x, y, z + 1));
            if (x + 1 < d0 && y + 1 < d1) face(s(x, y, z), s(x + 1, y, z), s(x + 1, y + 1, z), s(x, y + 1, z));
        }
        return count;
    }

    std::pair<int64_t, int64_t> MinMax(std::vector<int64_t> const& v)
    {
        auto mm = std::minmax_element(v.begin(), v.end());
        return { *mm.first, *mm.second };
    }

    int BitsFor(int64_t range)
    {
        int b = 1;
        while ((int64_t(1) << b) <= range) { ++b; }
        return b;
    }

    // Records voxel values in [lo, hi] packed as v - lo.
    void GivenVoxels(oracle::Ctx& io, std::vector<int64_t> const& v, int64_t lo, int64_t hi)
    {
        std::vector<int64_t> u;
        for (int64_t x : v) { u.push_back(x - lo); }
        GivenPacked(io, u, BitsFor(hi - lo));
    }

    // Emits, for one extractor and level: the rational extraction (counts,
    // digest), its MakeUnique (counts, digest), the real extraction without
    // duplicate removal (count, digest) and with it (counts; vertices and
    // triangles in full when there are at most 16 triangles and full is
    // requested, else a digest), OrientTriangles (one packed bit per
    // triangle: swapped) and ComputeNormals on the deduplicated mesh (full
    // or digest like the vertices).
    template <typename T, typename Extractor>
    void EmitSurface(oracle::Ctx& io, Extractor& ex, T level, bool sameDir, bool allowFull)
    {
        using Vertex = typename Extractor::Vertex;
        using Triangle = typename Extractor::Triangle;
        SurfaceExtractor<T, double>& base = ex;
        std::vector<Vertex> rv;
        std::vector<Triangle> rt;
        base.Extract(level, rv, rt);
        io.outInt(rv.size());
        io.outInt(rt.size());
        Digest d0;
        for (auto const& v : rv)
        {
            d0.integer(v.xNumer); d0.integer(v.xDenom); d0.integer(v.yNumer);
            d0.integer(v.yDenom); d0.integer(v.zNumer); d0.integer(v.zDenom);
        }
        for (auto const& t : rt) { d0.integer(t.v[0]); d0.integer(t.v[1]); d0.integer(t.v[2]); }
        io.outInt(d0.h);

        base.MakeUnique(rv, rt);
        io.outInt(rv.size());
        io.outInt(rt.size());
        Digest d1;
        for (auto const& v : rv)
        {
            d1.integer(v.xNumer); d1.integer(v.xDenom); d1.integer(v.yNumer);
            d1.integer(v.yDenom); d1.integer(v.zNumer); d1.integer(v.zDenom);
        }
        for (auto const& t : rt) { d1.integer(t.v[0]); d1.integer(t.v[1]); d1.integer(t.v[2]); }
        io.outInt(d1.h);

        std::vector<std::array<double, 3>> xv;
        std::vector<Triangle> xt;
        base.Extract(level, false, xv, xt);
        io.outInt(xv.size());
        Digest d2;
        for (auto const& v : xv) { d2.real(v[0]); d2.real(v[1]); d2.real(v[2]); }
        for (auto const& t : xt) { d2.integer(t.v[0]); d2.integer(t.v[1]); d2.integer(t.v[2]); }
        io.outInt(d2.h);

        base.Extract(level, true, xv, xt);
        io.outInt(xv.size());
        io.outInt(xt.size());
        bool const full = allowFull && xt.size() <= 16;
        Digest d3;
        for (auto const& v : xv)
        {
            if (full) { io.outReal(v[0]); io.outReal(v[1]); io.outReal(v[2]); }
            else { d3.real(v[0]); d3.real(v[1]); d3.real(v[2]); }
        }
        for (auto const& t : xt)
        {
            if (full) { io.outInt(t.v[0]); io.outInt(t.v[1]); io.outInt(t.v[2]); }
            else { d3.integer(t.v[0]); d3.integer(t.v[1]); d3.integer(t.v[2]); }
        }
        if (!full) { io.outInt(d3.h); }

        std::vector<Triangle> ot = xt;
        base.OrientTriangles(xv, ot, sameDir);
        std::vector<int64_t> swapped(ot.size());
        for (size_t t = 0; t < ot.size(); ++t) { swapped[t] = (ot[t].v[1] != xt[t].v[1] ? 1 : 0); }
        OutPacked(io, swapped, 1);

        std::vector<std::array<double, 3>> normals;
        base.ComputeNormals(xv, ot, normals);
        Digest d4;
        for (auto const& n : normals)
        {
            if (full) { io.outReal(n[0]); io.outReal(n[1]); io.outReal(n[2]); }
            else { d4.real(n[0]); d4.real(n[1]); d4.real(n[2]); }
        }
        if (!full) { io.outInt(d4.h); }
    }
}

// ---- SurfaceExtractorCubes ----
//
// Upstream pairs the four crossings of a saddle face the wrong way round on
// every face whose determinant is nonzero (this report, upstream suspect 1;
// fixed in the port). The main cases reject images that have such a face
// (exact integer predicate CubesSaddles, at most 64 redraws, then a linear
// image, which has none); plus-sign faces (determinant 0) are sound and
// have their own case; the deviation case requires a defective face.

namespace
{
    // Draws dimensions, values and level until 'accept' holds (at most 64
    // attempts, then a linear image and a level inside its range), records
    // them, and returns the values and level.
    template <typename Accept>
    std::vector<int64_t> DrawCubesImage(oracle::Ctx& io, int d0, int d1, int d2,
        int64_t lo, int64_t hi, int maxMode, int64_t& level, Accept accept)
    {
        std::vector<int64_t> v;
        bool found = false;
        for (int attempt = 0; attempt < 64 && !found; ++attempt)
        {
            v = RawVoxels(io, d0, d1, d2, io.rawInteger(0, maxMode), lo, hi);
            auto mm = MinMax(v);
            level = std::max(lo, mm.first - 1) + static_cast<int64_t>(
                io.raw(0.0, 1.0) * static_cast<double>(mm.second - std::max(lo, mm.first - 1) + 1));
            found = accept(v, level);
        }
        if (!found)
        {
            v = RawVoxels(io, d0, d1, d2, 2, lo, hi);
            level = MinMax(v).first;
        }
        GivenVoxels(io, v, lo, hi);
        io.given(static_cast<double>(level));
        return v;
    }

    bool NoDefectiveSaddle(std::vector<int64_t> const& v, int d0, int d1, int d2, int64_t level)
    {
        return CubesSaddles(v, d0, d1, d2, level).nonzero == 0;
    }

    template <typename T>
    void CubesTyped(oracle::Ctx& io, int d0, int d1, int d2, std::vector<int64_t> const& v,
        int64_t level, bool sameDir, bool allowFull)
    {
        std::vector<T> voxels(v.size());
        for (size_t i = 0; i < v.size(); ++i) { voxels[i] = static_cast<T>(v[i]); }
        SurfaceExtractorCubes<T, double> ex(d0, d1, d2, voxels.data());
        EmitSurface<T>(io, ex, static_cast<T>(level), sameDir, allowFull);
    }
}

// int32_t images of 2..3 voxels per axis; small lattice, uniform [-50, 50],
// linear, quadric, binary and hyperbolic values; levels from one below the
// minimum to the maximum.
ORACLE_CASE("SurfaceExtractorCubes.extract")
{
    int d0 = io.integer(2, 3);
    int d1 = io.integer(2, 3);
    int d2 = io.integer(2, 3);
    int64_t level = 0;
    auto v = DrawCubesImage(io, d0, d1, d2, -50, 50, 5, level,
        [&](auto const& w, int64_t L) { return NoDefectiveSaddle(w, d0, d1, d2, L); });
    bool sameDir = io.boolean();
    CubesTyped<int32_t>(io, d0, d1, d2, v, level, sameDir, true);
}

// Every pixel type upstream documents. 32-bit values stay within 2^20 so
// that the port's double arithmetic (exact below 2^53) and upstream's int64
// products agree; digests only.
ORACLE_CASE("SurfaceExtractorCubes.extract.types")
{
    int type = io.integer(0, 5);
    int d0 = io.integer(2, 4);
    int d1 = io.integer(2, 4);
    int d2 = io.integer(2, 4);
    int64_t const los[6] = { -128, -32768, -(1 << 20), 0, 0, 0 };
    int64_t const his[6] = { 127, 32767, 1 << 20, 255, 65535, 1 << 20 };
    int64_t level = 0;
    auto v = DrawCubesImage(io, d0, d1, d2, los[type], his[type], 5, level,
        [&](auto const& w, int64_t L) { return NoDefectiveSaddle(w, d0, d1, d2, L); });
    bool sameDir = io.boolean();
    switch (type)
    {
    case 0: CubesTyped<int8_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    case 1: CubesTyped<int16_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    case 2: CubesTyped<int32_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    case 3: CubesTyped<uint8_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    case 4: CubesTyped<uint16_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    default: CubesTyped<uint32_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    }
}

// 4..7 voxels per axis (many cubes, long ear-clipping sequences); digests.
ORACLE_CASE("SurfaceExtractorCubes.extract.large")
{
    int d0 = io.integer(4, 6);
    int d1 = io.integer(4, 6);
    int d2 = io.integer(4, 6);
    int64_t level = 0;
    auto v = DrawCubesImage(io, d0, d1, d2, -200, 200, 5, level,
        [&](auto const& w, int64_t L) { return NoDefectiveSaddle(w, d0, d1, d2, L); });
    bool sameDir = io.boolean();
    CubesTyped<int32_t>(io, d0, d1, d2, v, level, sameDir, false);
}

// Plus-sign faces (four alternating corners with f00*f11 == f01*f10, the
// branch-point path): one face is planted with shifted values u*p, -u*q,
// w*q, -w*p (u, w, p, q odd), the rest drawn and rejected if a defective
// saddle appears; the fallback extrudes the planted face along its normal.
ORACLE_CASE("SurfaceExtractorCubes.extract.plusSign")
{
    int d0 = io.integer(2, 3);
    int d1 = io.integer(2, 3);
    int d2 = io.integer(2, 3);
    int64_t level = 0;
    std::vector<int64_t> v;
    bool found = false;
    for (int attempt = 0; attempt < 64 && !found; ++attempt)
    {
        level = io.rawInteger(-3, 3);
        int64_t const odd[2] = { 1, 3 };
        int64_t u = odd[io.rawInteger(0, 1)], w = odd[io.rawInteger(0, 1)];
        int64_t p = odd[io.rawInteger(0, 1)], q = odd[io.rawInteger(0, 1)];
        int64_t sgn = (io.rawInteger(0, 1) == 0 ? 1 : -1);
        std::array<int64_t, 4> s = { sgn * u * p, -sgn * u * q, sgn * w * q, -sgn * w * p };
        int axis = io.rawInteger(0, 2);
        int dims[3] = { d0, d1, d2 };
        int at = io.rawInteger(0, dims[axis] - 1);
        int b0 = io.rawInteger(0, dims[(axis + 1) % 3] - 2), b1 = io.rawInteger(0, dims[(axis + 2) % 3] - 2);
        // The last attempt gives every voxel the value of the planted corner
        // nearest to it in the face's plane, so that the planted squares are
        // the only alternating faces.
        bool extrude = (attempt == 63);
        v = RawVoxels(io, d0, d1, d2, 1, -30, 30);
        for (int z = 0; z < d2; ++z) for (int y = 0; y < d1; ++y) for (int x = 0; x < d0; ++x)
        {
            int c[3] = { x, y, z };
            int e0 = c[(axis + 1) % 3] - b0, e1 = c[(axis + 2) % 3] - b1;
            bool planted = (c[axis] == at && (e0 == 0 || e0 == 1) && (e1 == 0 || e1 == 1));
            if (!planted && !extrude) continue;
            e0 = std::min(1, std::max(0, e0));
            e1 = std::min(1, std::max(0, e1));
            int j = (e1 == 0 ? (e0 == 0 ? 0 : 1) : (e0 == 1 ? 2 : 3));
            v[x + d0 * (y + d1 * z)] = (s[j] + 2 * level + 1) / 2;
        }
        found = extrude || NoDefectiveSaddle(v, d0, d1, d2, level);
    }
    GivenVoxels(io, v, -30, 30);
    io.given(static_cast<double>(level));
    bool sameDir = io.boolean();
    CubesTyped<int32_t>(io, d0, d1, d2, v, level, sameDir, true);
}

// Upstream suspect 1 of this report: every record has at least one saddle
// face with a nonzero determinant (at most 64 redraws, then the 2x2x2
// image with 3 at (0,0,0) and (0,1,1), 0 elsewhere, level 0).
ORACLE_CASE("SurfaceExtractorCubes.extract.saddle")
{
    int d0 = io.integer(2, 3);
    int d1 = io.integer(2, 3);
    int d2 = io.integer(2, 3);
    int64_t level = 0;
    std::vector<int64_t> v;
    bool found = false;
    for (int attempt = 0; attempt < 64 && !found; ++attempt)
    {
        v = RawVoxels(io, d0, d1, d2, io.rawInteger(0, 5), -20, 20);
        auto mm = MinMax(v);
        level = mm.first - 1 + io.rawInteger(0, static_cast<int>(mm.second - mm.first + 1));
        found = CubesSaddles(v, d0, d1, d2, level).nonzero > 0;
    }
    if (!found)
    {
        v.assign(static_cast<size_t>(d0) * d1 * d2, 0);
        v[0] = 3;
        v[d0 * (1 + d1)] = 3;
        level = 0;
    }
    GivenVoxels(io, v, -20, 20);
    io.given(static_cast<double>(level));
    bool sameDir = io.boolean();
    CubesTyped<int32_t>(io, d0, d1, d2, v, level, sameDir, true);
}

namespace
{
    // One recorded vertex: uniform in the image box widened by 1/4 (points
    // just outside give the zero gradient; slightly negative coordinates
    // truncate to cube 0), an integer point, or a multiple of 1/4. 'allowed'
    // rejects points (at most 256 redraws, then the image center).
    template <typename Allowed>
    std::array<double, 3> DrawPoint(oracle::Ctx& io, int mode, int const* dims, Allowed allowed)
    {
        std::array<double, 3> p{};
        bool ok = false;
        for (int attempt = 0; attempt < 256 && !ok; ++attempt)
        {
            for (int k = 0; k < 3; ++k)
            {
                double hi = static_cast<double>(dims[k] - 1);
                p[k] = (mode == 0 ? io.raw(-0.25, hi + 0.25)
                    : mode == 1 ? static_cast<double>(io.rawInteger(0, dims[k] - 1))
                    : io.rawInteger(-1, 4 * (dims[k] - 1) + 1) / 4.0);
            }
            ok = allowed(p);
        }
        if (!ok)
        {
            for (int k = 0; k < 3; ++k) { p[k] = 0.5 * (dims[k] - 1); }
        }
        for (int k = 0; k < 3; ++k) { io.given(p[k]); }
        return p;
    }

    // Arbitrary vertices and triangles (repeated indices allowed) through
    // OrientTriangles and ComputeNormals. Outputs: one bit per triangle
    // (swapped) and the normals in full.
    template <typename T, typename Allowed>
    void OrientPoints(oracle::Ctx& io, SurfaceExtractor<T, double>& ex, int const* dims, Allowed allowed)
    {
        using Triangle = typename SurfaceExtractor<T, double>::Triangle;
        int nv = io.integer(3, 7);
        int mode = io.integer(0, 2);
        std::vector<std::array<double, 3>> xv;
        for (int i = 0; i < nv; ++i) { xv.push_back(DrawPoint(io, mode, dims, allowed)); }
        int nt = io.integer(1, 6);
        std::vector<Triangle> xt;
        for (int t = 0; t < nt; ++t)
        {
            int i0 = io.integer(0, nv - 1);
            int i1 = io.integer(0, nv - 1);
            int i2 = io.integer(0, nv - 1);
            xt.push_back(Triangle(i0, i1, i2));
        }
        bool sameDir = io.boolean();
        std::vector<Triangle> ot = xt;
        ex.OrientTriangles(xv, ot, sameDir);
        std::vector<int64_t> swapped(ot.size());
        for (size_t t = 0; t < ot.size(); ++t) { swapped[t] = (ot[t].v[1] != xt[t].v[1] ? 1 : 0); }
        OutPacked(io, swapped, 1);
        std::vector<std::array<double, 3>> normals;
        ex.ComputeNormals(xv, ot, normals);
        for (auto const& n : normals) { io.outReal(n[0]); io.outReal(n[1]); io.outReal(n[2]); }
    }
}

// GetGradient (trilinear, from the shifted voxels 2v - 2L - 1 that Extract
// leaves behind) through OrientTriangles on arbitrary points, including
// integer points on lattice images, where the average gradient can be
// exactly orthogonal to the normal and the rounding of the division by 3
// decides the sign.
ORACLE_CASE("SurfaceExtractorCubes.orientTriangles.points")
{
    int dims[3];
    for (int k = 0; k < 3; ++k) { dims[k] = io.integer(2, 4); }
    std::vector<int64_t> v = RawVoxels(io, dims[0], dims[1], dims[2], io.rawInteger(0, 5), -20, 20);
    GivenVoxels(io, v, -20, 20);
    int32_t level = io.integer(-3, 3);
    std::vector<int32_t> voxels(v.begin(), v.end());
    SurfaceExtractorCubes<int32_t, double> ex(dims[0], dims[1], dims[2], voxels.data());
    std::vector<SurfaceExtractorCubes<int32_t, double>::Vertex> rv;
    std::vector<SurfaceExtractorCubes<int32_t, double>::Triangle> rt;
    ex.Extract(level, rv, rt);
    OrientPoints<int32_t>(io, ex, dims, [](auto const&) { return true; });
}

// Constructor precondition: every bound at least 2 (LogAssert). Draws which
// extractor and bounds in 0..3; the throw records carry no outputs.
ORACLE_CASE("SurfaceExtractor.invalidBounds")
{
    int which = io.integer(0, 1);
    int d0 = io.integer(0, 3);
    int d1 = io.integer(0, 3);
    int d2 = io.integer(0, 3);
    std::vector<int32_t> voxels(64, 1);
    if (which == 0)
    {
        SurfaceExtractorCubes<int32_t, double> ex(d0, d1, d2, voxels.data());
    }
    else
    {
        SurfaceExtractorTetrahedra<int32_t, double> ex(d0, d1, d2, voxels.data());
    }
    io.outInt(d0 * d1 * d2);
}

// ---- SurfaceExtractorMC (T = double, IndexType = int32_t) ----
//
// #443 (fixed in the port): the edge interpolation omits 'level'. The main
// cases use level +0 or -0 (where upstream is sound) and compare every
// output; the '.levelTopology' cases use nonzero levels and compare what
// does not depend on the vertex positions; the '.level' cases are the
// deviations.

namespace
{
    using MCExtractor = SurfaceExtractorMC<double, int32_t>;

    // One recorded value: uniform, a small lattice (corners on the level),
    // a multiple of 1/8, or a signed zero.
    double DrawMCValue(oracle::Ctx& io, int mode)
    {
        switch (mode)
        {
        case 0: return io.real(-1.0, 1.0);
        case 1: return io.given(static_cast<double>(io.rawInteger(-2, 2)));
        case 2: return io.given(io.rawInteger(-16, 16) / 8.0);
        default: return io.given(io.rawInteger(0, 1) == 0 ? 0.0 : -0.0);
        }
    }

    // perturb: 0, -0, +-1/1024, +-1e-300 (a perturbation that survives the
    // addition), or uniform.
    double DrawPerturb(oracle::Ctx& io)
    {
        switch (io.rawInteger(0, 4))
        {
        case 0: return io.given(0.0);
        case 1: return io.given(-0.0);
        case 2: return io.given(io.rawInteger(0, 1) == 0 ? 1.0 / 1024.0 : -1.0 / 1024.0);
        case 3: return io.given(io.rawInteger(0, 1) == 0 ? 1e-300 : -1e-300);
        default: return io.real(-0.01, 0.01);
        }
    }

    double DrawLevel(oracle::Ctx& io, bool zero)
    {
        if (zero) { return io.given(io.rawInteger(0, 1) == 0 ? 0.0 : -0.0); }
        switch (io.rawInteger(0, 2))
        {
        case 0: return io.real(-1.0, 1.0);
        case 1: return io.given(io.rawInteger(1, 16) * (io.rawInteger(0, 1) == 0 ? 0.125 : -0.125));
        default: return io.given(static_cast<double>(io.rawInteger(0, 1) == 0 ? 1 : -1));
        }
    }

    __declspec(noinline) bool CallExtractVoxel(MCExtractor const& mc, double level, double perturb,
        std::array<double, 8> const& F, MCExtractor::Mesh& mesh)
    {
        return mc.Extract(level, perturb, F, mesh);
    }

    // Per-voxel extraction. Outputs: valid; when valid the counts, the
    // vertex pairs and triangle triples (4 bits packed) and, if requested,
    // the local vertices.
    void MCVoxel(oracle::Ctx& io, bool zeroLevel, bool withVertices)
    {
        Image3<double> image(2, 2, 2);
        MCExtractor mc(image);
        int mode = io.integer(0, 3);
        std::array<double, 8> F{};
        for (int i = 0; i < 8; ++i) { F[i] = DrawMCValue(io, (mode == 3 && io.rawInteger(0, 1) == 0) ? 1 : mode); }
        double level = DrawLevel(io, zeroLevel);
        double perturb = DrawPerturb(io);
        MCExtractor::Mesh mesh;
        bool valid = CallExtractVoxel(mc, level, perturb, F, mesh);
        io.outBool(valid);
        if (!valid) { return; }
        io.outInt(mesh.topology.numVertices);
        io.outInt(mesh.topology.numTriangles);
        std::vector<int64_t> packed;
        for (int i = 0; i < mesh.topology.numVertices; ++i)
        {
            packed.push_back(mesh.topology.vpair[i][0]);
            packed.push_back(mesh.topology.vpair[i][1]);
        }
        for (int i = 0; i < mesh.topology.numTriangles; ++i)
        {
            for (int j = 0; j < 3; ++j) { packed.push_back(mesh.topology.itriple[i][j]); }
        }
        OutPacked(io, packed, 4);
        if (withVertices)
        {
            for (int i = 0; i < mesh.topology.numVertices; ++i) { io.outVec(mesh.vertices[i]); }
        }
    }
}

ORACLE_CASE("SurfaceExtractorMC.extractVoxel") { MCVoxel(io, true, true); }
ORACLE_CASE("SurfaceExtractorMC.extractVoxel.levelTopology") { MCVoxel(io, false, false); }
ORACLE_CASE("SurfaceExtractorMC.extractVoxel.level") { MCVoxel(io, false, true); }

namespace
{
    // An Image3<double> of 2..3 voxels per axis: per-voxel values as in
    // DrawMCValue, a quadric, or a hyperbolic product (ambiguous faces).
    Image3<double> DrawMCImage(oracle::Ctx& io, int dmin, int dmax)
    {
        int d0 = io.integer(dmin, dmax);
        int d1 = io.integer(dmin, dmax);
        int d2 = io.integer(dmin, dmax);
        int mode = io.integer(0, 5);
        Image3<double> image(d0, d1, d2);
        double cx = io.raw(0.0, d0 - 1.0), cy = io.raw(0.0, d1 - 1.0), cz = io.raw(0.0, d2 - 1.0);
        double r = io.raw(0.25, 2.0);
        for (size_t i = 0; i < image.GetNumPixels(); ++i)
        {
            auto c = image.GetCoordinates(i);
            double x = c[0], y = c[1], z = c[2];
            double value = 0.0;
            switch (mode)
            {
            case 4: value = io.given((x - cx) * (x - cx) + (y - cy) * (y - cy) + (z - cz) * (z - cz) - r * r); break;
            case 5: value = io.given((x - cx) * (y - cy) + 0.25 * (z - cz)); break;
            default: value = DrawMCValue(io, mode); break;
            }
            image[i] = value;
        }
        return image;
    }

    // Whole-image extraction. Outputs: vertex and index counts, then the
    // vertices and indices (in full when there are at most 16 triangles and
    // full is requested, else a digest); when there are vertices, the
    // MakeUnique result (same rule), OrientTriangles (one bit per triangle)
    // and ComputeNormals.
    void MCImage(oracle::Ctx& io, bool zeroLevel, bool topologyOnly, bool allowFull, int dmax)
    {
        Image3<double> image = DrawMCImage(io, 2, dmax);
        double level = DrawLevel(io, zeroLevel);
        double perturb = DrawPerturb(io);
        bool sameDir = io.boolean();
        MCExtractor mc(image);
        std::vector<Vector3<double>> vertices;
        std::vector<int32_t> indices;
        mc.Extract(level, perturb, vertices, indices);
        io.outInt(vertices.size());
        io.outInt(indices.size());
        bool full = allowFull && indices.size() <= 48;
        auto emitMesh = [&](std::vector<Vector3<double>> const& vs, std::vector<int32_t> const& is, bool withVertices)
        {
            Digest d;
            if (withVertices)
            {
                for (auto const& v : vs) { for (int k = 0; k < 3; ++k) { if (full) io.outReal(v[k]); else d.real(v[k]); } }
            }
            for (int32_t i : is) { if (full) io.outInt(i); else d.integer(i); }
            if (!full) { io.outInt(d.h); }
        };
        emitMesh(vertices, indices, !topologyOnly);
        if (topologyOnly || vertices.empty()) { return; }
        mc.MakeUnique(vertices, indices);
        io.outInt(vertices.size());
        emitMesh(vertices, indices, true);
        std::vector<int32_t> oriented = indices;
        mc.OrientTriangles(vertices, oriented, sameDir);
        std::vector<int64_t> swapped(oriented.size() / 3);
        for (size_t t = 0; t < swapped.size(); ++t) { swapped[t] = (oriented[3 * t + 1] != indices[3 * t + 1] ? 1 : 0); }
        OutPacked(io, swapped, 1);
        // ComputeNormals never advances its triangle pointer (this report,
        // upstream suspect 2; fixed in the port): it is compared here only
        // for single-triangle meshes, where upstream is sound, and in the
        // deviation case SurfaceExtractorMC.computeNormals.multiple.
        if (oriented.size() == 3)
        {
            std::vector<Vector3<double>> normals;
            mc.ComputeNormals(vertices, oriented, normals);
            for (auto const& n : normals) { io.outVec(n); }
        }
    }
}

ORACLE_CASE("SurfaceExtractorMC.extract") { MCImage(io, true, false, true, 3); }
ORACLE_CASE("SurfaceExtractorMC.extract.large") { MCImage(io, true, false, false, 6); }
ORACLE_CASE("SurfaceExtractorMC.extract.levelTopology") { MCImage(io, false, true, true, 3); }
ORACLE_CASE("SurfaceExtractorMC.extract.level") { MCImage(io, false, false, false, 3); }

// MakeUnique (UniqueVerticesSimplices::RemoveDuplicateVertices) on
// arbitrary inputs: duplicated vertices, -0 against +0 (equal under
// Vector3's operator<, the first one is kept), and the preconditions (no
// vertices, no indices or a count not divisible by 3, an index out of
// range: LogAssert, throw records). Outputs: the vertices and indices.
ORACLE_CASE("SurfaceExtractorMC.makeUnique")
{
    double const palette[5] = { 0.0, -0.0, 1.0, 0.5, 2.0 };
    int mode = io.integer(0, 4);
    int nv = (mode == 2 ? 0 : io.integer(1, 8));
    std::vector<Vector3<double>> vertices(nv);
    for (auto& v : vertices)
    {
        for (int k = 0; k < 3; ++k) { v[k] = io.given(palette[io.rawInteger(0, 4)]); }
    }
    int ni = (mode == 3 ? io.integer(0, 2) : 3 * io.integer(1, 4));
    std::vector<int32_t> indices(ni);
    for (auto& i : indices) { i = io.integer(0, std::max(0, nv - 1)); }
    if (mode == 4)
    {
        int bad = io.integer(0, ni - 1);
        indices[bad] = (io.boolean() ? nv : -1);
    }
    Image3<double> image(2, 2, 2);
    MCExtractor mc(image);
    mc.MakeUnique(vertices, indices);
    io.outInt(vertices.size());
    for (auto const& v : vertices) { io.outVec(v); }
    for (int32_t i : indices) { io.outInt(i); }
}

// GetGradient (trilinear on the raw image, floor-based cell selection,
// zero outside) through OrientTriangles on arbitrary points and index
// triples; the integer points on lattice images make the average gradient
// exactly orthogonal to the normal, where the rounding of the division by
// 3 decides the sign (Vector3's operator/ multiplies by 1/3).
ORACLE_CASE("SurfaceExtractorMC.orientTriangles.points")
{
    Image3<double> image = DrawMCImage(io, 2, 4);
    int dims[3] = { image.GetDimension(0), image.GetDimension(1), image.GetDimension(2) };
    int nv = io.integer(3, 7);
    int mode = io.integer(0, 2);
    std::vector<Vector3<double>> vertices;
    for (int i = 0; i < nv; ++i)
    {
        auto p = DrawPoint(io, mode, dims, [](auto const&) { return true; });
        vertices.push_back(Vector3<double>{ p[0], p[1], p[2] });
    }
    int nt = io.integer(1, 6);
    std::vector<int32_t> indices;
    for (int t = 0; t < 3 * nt; ++t) { indices.push_back(io.integer(0, nv - 1)); }
    bool sameDir = io.boolean();
    MCExtractor mc(image);
    std::vector<int32_t> oriented = indices;
    mc.OrientTriangles(vertices, oriented, sameDir);
    std::vector<int64_t> swapped(nt);
    for (int t = 0; t < nt; ++t) { swapped[t] = (oriented[3 * t + 1] != indices[3 * t + 1] ? 1 : 0); }
    OutPacked(io, swapped, 1);
    // ComputeNormals only where upstream is sound (one triangle), see
    // MCImage.
    if (nt == 1)
    {
        std::vector<Vector3<double>> normals;
        mc.ComputeNormals(vertices, oriented, normals);
        for (auto const& n : normals) { io.outVec(n); }
    }
}

// Upstream suspect 2 of this report (fixed in the port): ComputeNormals
// reads indices.data() for every triangle without advancing, so only the
// first triangle's normal is accumulated, numTriangles times, at its three
// vertices; every other vertex gets the zero normal. Arbitrary points and
// 2..6 triangles; outputs the normals.
ORACLE_CASE("SurfaceExtractorMC.computeNormals.multiple")
{
    Image3<double> image(2, 2, 2);
    int dims[3] = { 3, 3, 3 };
    int nv = io.integer(3, 7);
    int mode = io.integer(0, 2);
    std::vector<Vector3<double>> vertices;
    for (int i = 0; i < nv; ++i)
    {
        auto p = DrawPoint(io, mode, dims, [](auto const&) { return true; });
        vertices.push_back(Vector3<double>{ p[0], p[1], p[2] });
    }
    int nt = io.integer(2, 6);
    std::vector<int32_t> indices;
    for (int t = 0; t < 3 * nt; ++t) { indices.push_back(io.integer(0, nv - 1)); }
    MCExtractor mc(image);
    std::vector<Vector3<double>> normals;
    mc.ComputeNormals(vertices, indices, normals);
    for (auto const& n : normals) { io.outVec(n); }
}

// ---- SurfaceExtractorTetrahedra ----
//
// Levels on sample values are allowed (the zero-corner cases of
// ProcessTetrahedron). #132 (fixed in the port): GetGradient's test for the
// corner tetrahedron at (1,1,1) of an odd-parity cube is dx + dy + dz >= 0
// (always true) instead of >= 2, so the central tetrahedron 0752 is dead
// code upstream. Extracted vertices always lie on a cube face and never
// reach it; the '.points' cases reject points in the region, the
// '.centralTetra' case is the deviation.

namespace
{
    template <typename T>
    void TetraTyped(oracle::Ctx& io, int d0, int d1, int d2, std::vector<int64_t> const& v,
        int64_t level, bool sameDir, bool allowFull)
    {
        std::vector<T> voxels(v.size());
        for (size_t i = 0; i < v.size(); ++i) { voxels[i] = static_cast<T>(v[i]); }
        SurfaceExtractorTetrahedra<T, double> ex(d0, d1, d2, voxels.data());
        EmitSurface<T>(io, ex, static_cast<T>(level), sameDir, allowFull);
    }

    // Draws values (modes of RawVoxels) and a level from one below the
    // minimum to the maximum; lattice modes put the level on samples.
    std::vector<int64_t> DrawTetraImage(oracle::Ctx& io, int d0, int d1, int d2,
        int64_t lo, int64_t hi, int64_t& level)
    {
        std::vector<int64_t> v = RawVoxels(io, d0, d1, d2, io.rawInteger(0, 5), lo, hi);
        auto mm = MinMax(v);
        int64_t l0 = std::max(lo, mm.first - 1);
        level = l0 + static_cast<int64_t>(io.raw(0.0, 1.0) * static_cast<double>(mm.second - l0 + 1));
        GivenVoxels(io, v, lo, hi);
        io.given(static_cast<double>(level));
        return v;
    }
}

ORACLE_CASE("SurfaceExtractorTetrahedra.extract")
{
    int d0 = io.integer(2, 3);
    int d1 = io.integer(2, 3);
    int d2 = io.integer(2, 3);
    int64_t level = 0;
    auto v = DrawTetraImage(io, d0, d1, d2, -50, 50, level);
    bool sameDir = io.boolean();
    TetraTyped<int32_t>(io, d0, d1, d2, v, level, sameDir, true);
}

ORACLE_CASE("SurfaceExtractorTetrahedra.extract.types")
{
    int type = io.integer(0, 5);
    int d0 = io.integer(2, 4);
    int d1 = io.integer(2, 4);
    int d2 = io.integer(2, 4);
    int64_t const los[6] = { -128, -32768, -(1 << 20), 0, 0, 0 };
    int64_t const his[6] = { 127, 32767, 1 << 20, 255, 65535, 1 << 20 };
    int64_t level = 0;
    auto v = DrawTetraImage(io, d0, d1, d2, los[type], his[type], level);
    bool sameDir = io.boolean();
    switch (type)
    {
    case 0: TetraTyped<int8_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    case 1: TetraTyped<int16_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    case 2: TetraTyped<int32_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    case 3: TetraTyped<uint8_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    case 4: TetraTyped<uint16_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    default: TetraTyped<uint32_t>(io, d0, d1, d2, v, level, sameDir, false); break;
    }
}

ORACLE_CASE("SurfaceExtractorTetrahedra.extract.large")
{
    int d0 = io.integer(4, 6);
    int d1 = io.integer(4, 6);
    int d2 = io.integer(4, 6);
    int64_t level = 0;
    auto v = DrawTetraImage(io, d0, d1, d2, -200, 200, level);
    bool sameDir = io.boolean();
    TetraTyped<int32_t>(io, d0, d1, d2, v, level, sameDir, false);
}

namespace
{
    // True where upstream's GetGradient and the port's differ (#132): an
    // odd-parity cube, outside the corner tetrahedra at (1,0,0), (0,1,0)
    // and (0,0,1), and dx + dy + dz < 2. Same truncation and differences as
    // GetGradient.
    bool InCentralOddTetra(std::array<double, 3> const& p, int const* dims)
    {
        int32_t c[3];
        for (int k = 0; k < 3; ++k)
        {
            c[k] = static_cast<int32_t>(p[k]);
            if (c[k] < 0 || c[k] + 1 >= dims[k]) return false;
        }
        if (((c[0] & 1) ^ (c[1] & 1) ^ (c[2] & 1)) == 0) return false;
        double dx = p[0] - static_cast<double>(c[0]);
        double dy = p[1] - static_cast<double>(c[1]);
        double dz = p[2] - static_cast<double>(c[2]);
        if (dx - dy - dz >= 0.0 || dx - dy + dz <= 0.0 || dx + dy - dz <= 0.0) return false;
        return dx + dy + dz < 2.0;
    }

    class TetraProbe : public SurfaceExtractorTetrahedra<int32_t, double>
    {
    public:
        TetraProbe(int32_t d0, int32_t d1, int32_t d2, int32_t const* voxels)
            : SurfaceExtractorTetrahedra<int32_t, double>(d0, d1, d2, voxels) {}

        std::array<double, 3> Upstream(std::array<double, 3> const& p) { return GetGradient(p); }

        // The gradient of the central tetrahedron 0752 (upstream's own
        // expressions in its dead branch), which the port returns there.
        std::array<double, 3> Central(std::array<double, 3> const& p)
        {
            int32_t x = static_cast<int32_t>(p[0]), y = static_cast<int32_t>(p[1]), z = static_cast<int32_t>(p[2]);
            int32_t i000 = x + mXBound * (y + mYBound * z);
            double f000 = static_cast<double>(mVoxels[i000]);
            double f110 = static_cast<double>(mVoxels[i000 + mXBound + 1]);
            double f101 = static_cast<double>(mVoxels[i000 + mXYBound + 1]);
            double f011 = static_cast<double>(mVoxels[i000 + mXYBound + mXBound]);
            return { 0.5 * (-f000 - f011 + f101 + f110), 0.5 * (-f000 + f011 - f101 + f110),
                0.5 * (-f000 + f011 + f101 - f110) };
        }
    };

    // Whether OrientTriangles decides differently for some triangle when
    // the gradients in the central tetrahedron are the port's.
    bool OrientationDiffers(TetraProbe& probe, int const* dims, std::vector<std::array<double, 3>> const& xv,
        std::vector<std::array<int32_t, 3>> const& tris, bool sameDir)
    {
        for (auto const& indices : tris)
        {
            // The order OrientTriangles sees (the constructor's rotation).
            SurfaceExtractorTetrahedra<int32_t, double>::Triangle t(indices[0], indices[1], indices[2]);
            std::array<double, 3> v[3] = { xv[t.v[0]], xv[t.v[1]], xv[t.v[2]] };
            std::array<double, 3> e1, e2, n;
            for (int k = 0; k < 3; ++k) { e1[k] = v[1][k] - v[0][k]; e2[k] = v[2][k] - v[0][k]; }
            n[0] = e1[1] * e2[2] - e1[2] * e2[1];
            n[1] = e1[2] * e2[0] - e1[0] * e2[2];
            n[2] = e1[0] * e2[1] - e1[1] * e2[0];
            bool swap[2];
            for (int which = 0; which < 2; ++which)
            {
                std::array<double, 3> g[3];
                for (int j = 0; j < 3; ++j)
                {
                    g[j] = (which == 1 && InCentralOddTetra(v[j], dims) ? probe.Central(v[j]) : probe.Upstream(v[j]));
                }
                double avr[3];
                for (int k = 0; k < 3; ++k) { avr[k] = (g[0][k] + g[1][k] + g[2][k]) / 3.0; }
                double dot = avr[0] * n[0] + avr[1] * n[1] + avr[2] * n[2];
                swap[which] = (sameDir ? dot < 0.0 : dot > 0.0);
            }
            if (swap[0] != swap[1]) return true;
        }
        return false;
    }
}

ORACLE_CASE("SurfaceExtractorTetrahedra.orientTriangles.points")
{
    int dims[3];
    for (int k = 0; k < 3; ++k) { dims[k] = io.integer(2, 4); }
    std::vector<int64_t> v = RawVoxels(io, dims[0], dims[1], dims[2], io.rawInteger(0, 5), -20, 20);
    GivenVoxels(io, v, -20, 20);
    int32_t level = io.integer(-3, 3);
    std::vector<int32_t> voxels(v.begin(), v.end());
    SurfaceExtractorTetrahedra<int32_t, double> ex(dims[0], dims[1], dims[2], voxels.data());
    std::vector<SurfaceExtractorTetrahedra<int32_t, double>::Vertex> rv;
    std::vector<SurfaceExtractorTetrahedra<int32_t, double>::Triangle> rt;
    ex.Extract(level, rv, rt);
    OrientPoints<int32_t>(io, ex, dims, [&](auto const& p) { return !InCentralOddTetra(p, dims); });
}

// #132 deviation. Every vertex lies in the central tetrahedron of an
// odd-parity cube, and the record is kept only when some triangle's
// orientation depends on which gradient is used there (at most 64 redraws
// of image, level, points, triangles and direction). Recorded in the
// layout of the '.points' case (point mode 0).
ORACLE_CASE("SurfaceExtractorTetrahedra.orientTriangles.centralTetra")
{
    int dims[3];
    dims[0] = io.integer(3, 4);
    dims[1] = io.integer(2, 4);
    dims[2] = io.integer(2, 4);
    std::vector<int64_t> v;
    int32_t level = 0;
    std::vector<std::array<double, 3>> xv;
    std::vector<std::array<int32_t, 3>> tris;
    bool sameDir = true;
    for (int attempt = 0; attempt < 64; ++attempt)
    {
        v = RawVoxels(io, dims[0], dims[1], dims[2], io.rawInteger(0, 5), -20, 20);
        level = io.rawInteger(-3, 3);
        xv.assign(io.rawInteger(3, 5), {});
        for (auto& p : xv)
        {
            p = { 1.5, 0.5, 0.5 };
            for (int k = 0; k < 256; ++k)
            {
                std::array<double, 3> q = { io.raw(0.0, dims[0] - 1.0), io.raw(0.0, dims[1] - 1.0), io.raw(0.0, dims[2] - 1.0) };
                if (InCentralOddTetra(q, dims)) { p = q; break; }
            }
        }
        tris.assign(io.rawInteger(1, 4), {});
        for (auto& t : tris)
        {
            for (int j = 0; j < 3; ++j) { t[j] = io.rawInteger(0, static_cast<int>(xv.size()) - 1); }
        }
        sameDir = (io.rawInteger(0, 1) == 1);
        std::vector<int32_t> voxels(v.begin(), v.end());
        TetraProbe probe(dims[0], dims[1], dims[2], voxels.data());
        std::vector<TetraProbe::Vertex> rv;
        std::vector<TetraProbe::Triangle> rt;
        probe.Extract(level, rv, rt);
        if (OrientationDiffers(probe, dims, xv, tris, sameDir)) break;
    }
    GivenVoxels(io, v, -20, 20);
    io.given(level);
    std::vector<int32_t> voxels(v.begin(), v.end());
    SurfaceExtractorTetrahedra<int32_t, double> ex(dims[0], dims[1], dims[2], voxels.data());
    std::vector<SurfaceExtractorTetrahedra<int32_t, double>::Vertex> rv;
    std::vector<SurfaceExtractorTetrahedra<int32_t, double>::Triangle> rt;
    ex.Extract(level, rv, rt);
    io.given(static_cast<double>(xv.size()));
    io.given(0.0);
    for (auto const& p : xv) { io.given(p[0]); io.given(p[1]); io.given(p[2]); }
    io.given(static_cast<double>(tris.size()));
    std::vector<SurfaceExtractorTetrahedra<int32_t, double>::Triangle> xt;
    for (auto const& t : tris)
    {
        io.given(t[0]);
        io.given(t[1]);
        io.given(t[2]);
        xt.push_back(SurfaceExtractorTetrahedra<int32_t, double>::Triangle(t[0], t[1], t[2]));
    }
    io.given(sameDir ? 1.0 : 0.0);
    auto ot = xt;
    ex.OrientTriangles(xv, ot, sameDir);
    std::vector<int64_t> swapped(ot.size());
    for (size_t t = 0; t < ot.size(); ++t) { swapped[t] = (ot[t].v[1] != xt[t].v[1] ? 1 : 0); }
    OutPacked(io, swapped, 1);
    std::vector<std::array<double, 3>> normals;
    ex.ComputeNormals(xv, ot, normals);
    for (auto const& n : normals) { io.outReal(n[0]); io.outReal(n[1]); io.outReal(n[2]); }
}
