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

// @@END@@
