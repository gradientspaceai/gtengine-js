// Group 25 (imaging): Image3, PdeFilter1, PdeFilter2, PdeFilter3,
// SurfaceExtractor, TetrahedraRasterizer, CurvatureFlow2, CurvatureFlow3,
// FastMarch2, FastMarch3, GaussianBlur2, GaussianBlur3,
// GradientAnisotropic2, GradientAnisotropic3.
//
// Everything is integer bookkeeping or IEEE-exact arithmetic (+ - * /, sqrt,
// ceil, floor) and is compared bit for bit, except the GradientAnisotropic
// filters, whose conductances call std::exp (tolerance 1e-12 on the values
// that depend on an exp; everything computed before the first exp and every
// cell no update writes is still held to bit identity).
//
// Large arrays (padded PDE buffers, neighbourhood tables, FastMarch times,
// voxel grids) are emitted as 32-bit FNV-1a digests computed identically on
// both sides (NaN canonicalized); each case comment says which. Upstream
// methods whose result can carry a signed zero are called through
// __declspec(noinline) wrappers (ORACLE.md, v23).
#define ORACLE_FAMILY "v25-imaging"
#include "Oracle.h"

#include <Mathematics/CurvatureFlow2.h>
#include <Mathematics/CurvatureFlow3.h>
#include <Mathematics/FastMarch2.h>
#include <Mathematics/FastMarch3.h>
#include <Mathematics/GaussianBlur2.h>
#include <Mathematics/GaussianBlur3.h>
#include <Mathematics/GradientAnisotropic2.h>
#include <Mathematics/GradientAnisotropic3.h>
#include <Mathematics/Image3.h>
#include <Mathematics/PdeFilter1.h>
#include <Mathematics/PdeFilter2.h>
#include <Mathematics/PdeFilter3.h>
#include <Mathematics/SurfaceExtractor.h>
#include <Mathematics/TetrahedraRasterizer.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// Finding #122: upstream GradientAnisotropic2/3::ComputeParameter iterates
// the padded coordinates 1..bound into the unpadded GetUx/GetUy/GetUz, so
// its last row reads Array2::mIndirect1[yBound + 2] (and Array3's
// equivalent) past the end of the row-pointer vector and dereferences it.
// That is undefined behaviour on every construction and every Update, so
// no golden record of the unmodified member can exist. These explicit
// specializations replace exactly that one member function for Real =
// double with the port's loop over 0..bound-1; the constructors, OnPreUpdate,
// OnUpdateSingle and the PdeFilter2/3 bases are upstream's code unchanged.
namespace gte
{
    template <>
    void GradientAnisotropic2<double>::ComputeParameter()
    {
        double gradMagSqr = 0.0;
        for (int32_t y = 0; y < this->mYBound; ++y)
        {
            for (int32_t x = 0; x < this->mXBound; ++x)
            {
                double ux = this->GetUx(x, y);
                double uy = this->GetUy(x, y);
                gradMagSqr += ux * ux + uy * uy;
            }
        }
        gradMagSqr /= static_cast<double>(this->mQuantity);
        mParameter = 1.0 / (mK * mK * gradMagSqr);
        mMHalfParameter = -0.5 * mParameter;
    }

    template <>
    void GradientAnisotropic3<double>::ComputeParameter()
    {
        double gradMagSqr = 0.0;
        for (int32_t z = 0; z < this->mZBound; ++z)
        {
            for (int32_t y = 0; y < this->mYBound; ++y)
            {
                for (int32_t x = 0; x < this->mXBound; ++x)
                {
                    double ux = this->GetUx(x, y, z);
                    double uy = this->GetUy(x, y, z);
                    double uz = this->GetUz(x, y, z);
                    gradMagSqr += ux * ux + uy * uy + uz * uz;
                }
            }
        }
        gradMagSqr /= static_cast<double>(this->mQuantity);
        mParameter = 1.0 / (mK * mK * gradMagSqr);
        mMHalfParameter = -0.5 * mParameter;
    }
}

using namespace gte;

namespace
{
    // ---- digests (mirrored in the replay) ----

    // FNV-1a over 32-bit words. Every NaN hashes as one canonical pattern:
    // MSVC's default NaN is negative, V8's positive, and the harness treats
    // any NaN as equal to any NaN.
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
            if (std::isnan(x)) { bits = 0x7FF8000000000000ull; }
            else { std::memcpy(&bits, &x, sizeof(bits)); }
            word(static_cast<uint32_t>(bits));
            word(static_cast<uint32_t>(bits >> 32));
        }
        void out(oracle::Ctx& io) const { io.outInt(h); }
    };

    int64_t Signed(size_t i) { return static_cast<int64_t>(i); }
}

// ---- Image3 ----

namespace
{
    // Positive five times in six; a nonpositive dimension leaves the image
    // empty.
    int32_t DrawDimension(oracle::Ctx& io, int32_t maxDim)
    {
        int32_t d = (io.rawInteger(0, 5) == 0 ? io.rawInteger(-2, 0) : io.rawInteger(1, maxDim));
        io.given(d);
        return d;
    }

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

// Construction (constructor, or Reconstruct over a 2x3x2 image), both
// GetIndex overloads on coordinates up to two outside the image (negative
// ones wrap in size_t and are emitted signed, which is the port's value),
// both GetCoordinates overloads on indices up to twice the voxel count,
// operator() in both forms for writes and in-range reads, the clamping Get
// in both forms, and the base Image::Get(i) (hidden by Image3's overloads,
// reached through a base reference; the port's one-argument getClamped).
ORACLE_CASE("Image3.access")
{
    int32_t d0 = DrawDimension(io, 5);
    int32_t d1 = DrawDimension(io, 5);
    int32_t d2 = DrawDimension(io, 5);
    bool viaReconstruct = io.boolean();
    Image3<double> image;
    if (viaReconstruct)
    {
        image.Reconstruct(2, 3, 2);
        image.Reconstruct(d0, d1, d2);
    }
    else
    {
        image = Image3<double>(d0, d1, d2);
    }
    OutImageShape(io, image);
    size_t const numPixels = image.GetNumPixels();
    if (numPixels == 0)
    {
        return;
    }
    int32_t const dim0 = image.GetDimension(0), dim1 = image.GetDimension(1);
    int32_t const dim2 = image.GetDimension(2);
    for (int32_t z = 0; z < dim2; ++z)
    {
        for (int32_t y = 0; y < dim1; ++y)
        {
            for (int32_t x = 0; x < dim0; ++x)
            {
                double v = x + 10.0 * y + 100.0 * z;
                if ((x + y + z) % 2 == 0) { image(x, y, z) = v + 0.5; }
                else { image(std::array<int32_t, 3>{ x, y, z }) = -v - 0.25; }
            }
        }
    }
    int32_t numQueries = io.integer(1, 8);
    for (int32_t q = 0; q < numQueries; ++q)
    {
        int32_t x = io.integer(-2, dim0 + 1);
        int32_t y = io.integer(-2, dim1 + 1);
        int32_t z = io.integer(-2, dim2 + 1);
        std::array<int32_t, 3> coord{ x, y, z };
        io.outInt(Signed(image.GetIndex(x, y, z)));
        io.outInt(Signed(image.GetIndex(coord)));
        io.outReal(image.Get(x, y, z));
        io.outReal(image.Get(coord));
        if (0 <= x && x < dim0 && 0 <= y && y < dim1 && 0 <= z && z < dim2)
        {
            io.outReal(image(x, y, z));
            io.outReal(image(coord));
        }
        int32_t index = io.integer(0, static_cast<int32_t>(2 * numPixels));
        int32_t cx = -7, cy = -7, cz = -7;
        image.GetCoordinates(static_cast<size_t>(index), cx, cy, cz);
        io.outInt(cx);
        io.outInt(cy);
        io.outInt(cz);
        auto c = image.GetCoordinates(static_cast<size_t>(index));
        io.outInt(c[0]);
        io.outInt(c[1]);
        io.outInt(c[2]);
        int32_t r = io.integer(-2, static_cast<int32_t>(numPixels) + 1);
        Image<double> const& base = image;
        io.outReal(base.Get(static_cast<size_t>(static_cast<int64_t>(r))));
    }
}

namespace
{
    template <size_t K>
    void DigestTable(oracle::Ctx& io, std::array<int32_t, K> const& nbr)
    {
        Digest d;
        for (auto v : nbr) { d.integer(v); }
        d.out(io);
    }

    template <size_t K>
    void DigestTable(oracle::Ctx& io, std::array<std::array<int32_t, 3>, K> const& nbr)
    {
        Digest d;
        for (auto const& v : nbr) { d.integer(v[0]); d.integer(v[1]); d.integer(v[2]); }
        d.out(io);
    }

    // size_t tables: the digest takes the low 32 bits, which a wrapped
    // SIZE_MAX - k and the port's -1 - k share, so the digest compares the
    // signed values (the port's representation).
    template <size_t K>
    void DigestTable(oracle::Ctx& io, std::array<size_t, K> const& nbr)
    {
        Digest d;
        for (auto v : nbr) { d.integer(Signed(v)); }
        d.out(io);
    }

    template <size_t K>
    void DigestTable(oracle::Ctx& io, std::array<std::array<size_t, 3>, K> const& nbr)
    {
        Digest d;
        for (auto const& v : nbr)
        {
            d.integer(Signed(v[0])); d.integer(Signed(v[1])); d.integer(Signed(v[2]));
        }
        d.out(io);
    }

    // Every relative table (1-D offsets and 3-tuples; 6, 18, 26, corners,
    // full) as a digest, the 6-neighbourhood in full, then the absolute
    // tables at 'count' voxels inside, on and one outside the border: the
    // 1-D 6-neighbourhood in full (signed), every other table as a digest.
    void NeighborhoodCase(oracle::Ctx& io)
    {
        int32_t d0 = io.integer(1, 6);
        int32_t d1 = io.integer(1, 6);
        int32_t d2 = io.integer(1, 6);
        Image3<double> image(d0, d1, d2);
        std::array<int32_t, 6> n6{};
        std::array<int32_t, 18> n18{};
        std::array<int32_t, 26> n26{};
        std::array<int32_t, 8> c8{};
        std::array<int32_t, 27> f27{};
        image.GetNeighborhood(n6);
        image.GetNeighborhood(n18);
        image.GetNeighborhood(n26);
        image.GetCorners(c8);
        image.GetFull(f27);
        std::array<std::array<int32_t, 3>, 6> t6{};
        std::array<std::array<int32_t, 3>, 18> t18{};
        std::array<std::array<int32_t, 3>, 26> t26{};
        std::array<std::array<int32_t, 3>, 8> tc8{};
        std::array<std::array<int32_t, 3>, 27> tf27{};
        image.GetNeighborhood(t6);
        image.GetNeighborhood(t18);
        image.GetNeighborhood(t26);
        image.GetCorners(tc8);
        image.GetFull(tf27);
        for (auto v : n6) { io.outInt(v); }
        DigestTable(io, n18);
        DigestTable(io, n26);
        DigestTable(io, c8);
        DigestTable(io, f27);
        DigestTable(io, t6);
        DigestTable(io, t18);
        DigestTable(io, t26);
        DigestTable(io, tc8);
        DigestTable(io, tf27);

        int32_t count = io.integer(1, 3);
        for (int32_t k = 0; k < count; ++k)
        {
            int32_t x = io.integer(-1, d0);
            int32_t y = io.integer(-1, d1);
            int32_t z = io.integer(-1, d2);
            std::array<size_t, 6> a6{};
            std::array<size_t, 18> a18{};
            std::array<size_t, 26> a26{};
            std::array<size_t, 8> ac8{};
            std::array<size_t, 27> af27{};
            image.GetNeighborhood(x, y, z, a6);
            image.GetNeighborhood(x, y, z, a18);
            image.GetNeighborhood(x, y, z, a26);
            image.GetCorners(x, y, z, ac8);
            image.GetFull(x, y, z, af27);
            for (auto v : a6) { io.outInt(Signed(v)); }
            DigestTable(io, a18);
            DigestTable(io, a26);
            DigestTable(io, ac8);
            DigestTable(io, af27);
            std::array<std::array<size_t, 3>, 6> b6{};
            std::array<std::array<size_t, 3>, 18> b18{};
            std::array<std::array<size_t, 3>, 26> b26{};
            std::array<std::array<size_t, 3>, 8> bc8{};
            std::array<std::array<size_t, 3>, 27> bf27{};
            image.GetNeighborhood(x, y, z, b6);
            image.GetNeighborhood(x, y, z, b18);
            image.GetNeighborhood(x, y, z, b26);
            image.GetCorners(x, y, z, bc8);
            image.GetFull(x, y, z, bf27);
            DigestTable(io, b6);
            DigestTable(io, b18);
            DigestTable(io, b26);
            DigestTable(io, bc8);
            DigestTable(io, bf27);
        }
    }
}

ORACLE_CASE("Image3.neighborhoods") { NeighborhoodCase(io); }

// Finding #64: at a voxel on the xmin, ymin or zmin face upstream's size_t
// neighbour coordinates wrap to SIZE_MAX, the port's are -1. The 3-tuple
// 6-neighbourhood and the 1-D 6-neighbourhood are emitted as the unsigned
// values converted to double; every record deviates (its voxel has a zero
// coordinate). Image3.neighborhoods compares the same tables as signed
// values on the full generator.
ORACLE_CASE("Image3.neighborhoods.wrap")
{
    int32_t d0 = io.integer(1, 6);
    int32_t d1 = io.integer(1, 6);
    int32_t d2 = io.integer(1, 6);
    Image3<double> image(d0, d1, d2);
    int32_t face = io.rawInteger(0, 2);
    int32_t x = static_cast<int32_t>(io.given(face == 0 ? 0 : io.rawInteger(0, d0 - 1)));
    int32_t y = static_cast<int32_t>(io.given(face == 1 ? 0 : io.rawInteger(0, d1 - 1)));
    int32_t z = static_cast<int32_t>(io.given(face == 2 ? 0 : io.rawInteger(0, d2 - 1)));
    std::array<size_t, 6> a6{};
    image.GetNeighborhood(x, y, z, a6);
    for (auto v : a6) { io.outReal(static_cast<double>(v)); }
    std::array<std::array<size_t, 3>, 6> b6{};
    image.GetNeighborhood(x, y, z, b6);
    for (auto const& v : b6)
    {
        io.outReal(static_cast<double>(v[0]));
        io.outReal(static_cast<double>(v[1]));
        io.outReal(static_cast<double>(v[2]));
    }
}

// ---- PDE filters: shared input generators ----

namespace
{
    double const kMax = std::numeric_limits<double>::max();

    // One recorded double per pixel, x fastest. Modes: uniform, all
    // positive, a small lattice with ties, a constant image, signed zeros
    // and dyadic values, an image that depends on x only (a lattice
    // quadratic: its level sets are planes, so the curvature flows leave the
    // interior unchanged, checked in the replay), a single spike on zeros,
    // a signed checkerboard.
    std::vector<double> DrawData(oracle::Ctx& io, int32_t xB, int32_t yB, int32_t zB)
    {
        int mode = io.rawInteger(0, 7);
        double c = io.raw(-5.0, 5.0);
        double a = io.rawInteger(-3, 3), b = io.rawInteger(-2, 2), q = io.rawInteger(-1, 1);
        int32_t n = xB * yB * zB;
        int32_t spike = io.rawInteger(0, n - 1);
        std::vector<double> data(n);
        for (int32_t z = 0, i = 0; z < zB; ++z)
        {
            for (int32_t y = 0; y < yB; ++y)
            {
                for (int32_t x = 0; x < xB; ++x, ++i)
                {
                    double v;
                    switch (mode)
                    {
                    case 0: v = io.raw(-10.0, 10.0); break;
                    case 1: v = io.raw(0.5, 10.0); break;
                    case 2: v = io.rawInteger(-3, 3); break;
                    case 3: v = c; break;
                    case 4: v = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0) * io.rawInteger(0, 3) * 0.25; break;
                    case 5: v = a + b * x + q * x * x; break;
                    case 6: v = (i == spike ? c : 0.0); break;
                    default: v = ((x + y + z) % 2 == 0 ? 1.0 : -1.0) * io.rawInteger(1, 3); break;
                    }
                    data[i] = io.given(v);
                }
            }
        }
        return data;
    }

    // Mask: none (the recorded mode is 0), or one recorded integer per
    // pixel: random 0/1 mostly 1, all ones, a box of ones inside zeros, or
    // values in {-1, 0, 1, 2} (upstream tests the mask as a boolean).
    std::vector<int32_t> DrawMask(oracle::Ctx& io, int32_t xB, int32_t yB, int32_t zB, bool& hasMask)
    {
        int mode = io.integer(0, 4);
        hasMask = (mode != 0);
        std::vector<int32_t> mask;
        if (!hasMask)
        {
            return mask;
        }
        int32_t bx0 = io.rawInteger(0, xB - 1), bx1 = io.rawInteger(bx0, xB - 1);
        int32_t by0 = io.rawInteger(0, yB - 1), by1 = io.rawInteger(by0, yB - 1);
        int32_t bz0 = io.rawInteger(0, zB - 1), bz1 = io.rawInteger(bz0, zB - 1);
        for (int32_t z = 0; z < zB; ++z)
        {
            for (int32_t y = 0; y < yB; ++y)
            {
                for (int32_t x = 0; x < xB; ++x)
                {
                    int32_t m;
                    switch (mode)
                    {
                    case 1: m = (io.rawInteger(0, 3) == 0 ? 0 : 1); break;
                    case 2: m = 1; break;
                    case 3: m = (bx0 <= x && x <= bx1 && by0 <= y && y <= by1
                        && bz0 <= z && z <= bz1 ? 1 : 0); break;
                    default: m = io.rawInteger(-1, 2); break;
                    }
                    mask.push_back(static_cast<int32_t>(io.given(m)));
                }
            }
        }
        return mask;
    }

    // One recorded double: the Neumann marker (double max) half the time,
    // otherwise a Dirichlet value (0, a lattice value or uniform).
    double DrawBorder(oracle::Ctx& io)
    {
        switch (io.rawInteger(0, 5))
        {
        case 0: case 1: case 2: return io.given(kMax);
        case 3: return io.given(0.0);
        case 4: return io.given(io.rawInteger(-2, 2));
        default: return io.given(io.raw(-2.0, 2.0));
        }
    }

    // One recorded integer: the four named scale types, and 4 (outside the
    // enum; upstream's switch falls through and offset = scale = 0) rarely.
    int32_t DrawScaleType(oracle::Ctx& io)
    {
        return static_cast<int32_t>(io.given(io.rawInteger(0, 15) == 0 ? 4 : io.rawInteger(0, 3)));
    }

    double DrawSpacing(oracle::Ctx& io)
    {
        switch (io.rawInteger(0, 4))
        {
        case 0: case 1: return io.given(1.0);
        case 2: return io.given(0.5);
        case 3: return io.given(2.0);
        default: return io.given(io.raw(0.3, 2.0));
        }
    }

    // One recorded double: a stable step, a dyadic one, uniform, zero,
    // negative and a large (unstable) one.
    double DrawTimeStep(oracle::Ctx& io, double stable)
    {
        switch (io.rawInteger(0, 6))
        {
        case 0: case 1: return io.given(stable);
        case 2: return io.given(0.125);
        case 3: return io.given(io.raw(0.0, 0.3));
        case 4: return io.given(0.0);
        case 5: return io.given(-0.1);
        default: return io.given(io.raw(0.5, 2.0));
        }
    }

    // Bounds 1..maxB, with the maximum (5x5, 5x5x5) often.
    int32_t DrawBound(oracle::Ctx& io, int32_t maxB)
    {
        return static_cast<int32_t>(io.given(io.rawInteger(0, 2) == 0 ? maxB : io.rawInteger(1, maxB)));
    }

    using ST = PdeFilter<double>::ScaleType;

    __declspec(noinline) void UpdateFilter(PdeFilter<double>& filter) { filter.Update(); }
}

// ---- PdeFilter1 ----

namespace
{
    // PdeFilter1 has no concrete filter upstream. This subclass (defined
    // identically in the replay) is the explicit heat equation through the
    // base's LookUp3, and exposes the protected state.
    class Heat1 : public PdeFilter1<double>
    {
    public:
        Heat1(int32_t xBound, double xSpacing, double const* data, int32_t const* mask,
            double borderValue, ST scaleType)
            :
            PdeFilter1<double>(xBound, xSpacing, data, mask, borderValue, scaleType)
        {
        }

        double Min() const { return this->mMin; }
        double Offset() const { return this->mOffset; }
        double Scale() const { return this->mScale; }
        int32_t Src() const { return this->mSrc; }
        double Buf(int32_t k, int32_t x) const { return this->mBuffer[k][x]; }
        int32_t MaskAt(int32_t x) const { return this->mMask[x]; }

    protected:
        virtual void OnUpdate(int32_t x) override
        {
            this->LookUp3(x);
            this->mBuffer[this->mDst][x] = this->mUz + this->mTimeStep * this->mInvDxDx
                * (this->mUp - 2.0 * this->mUz + this->mUm);
        }
    };

    __declspec(noinline) Heat1* MakeHeat1(int32_t xBound, double xSpacing, double const* data,
        int32_t const* mask, double borderValue, int32_t scaleType)
    {
        return new Heat1(xBound, xSpacing, data, mask, borderValue, static_cast<ST>(scaleType));
    }

    // The source index, both padded buffers (digest), GetU at every pixel
    // in full, and GetUx, GetUxx, GetMask at every pixel plus the padded
    // mask (one digest).
    void OutState1(oracle::Ctx& io, Heat1 const& f)
    {
        int32_t const xB = f.GetXBound();
        io.outInt(f.Src());
        for (int32_t k = 0; k < 2; ++k)
        {
            Digest d;
            for (int32_t x = 0; x < xB + 2; ++x) { d.real(f.Buf(k, x)); }
            d.out(io);
        }
        Digest d;
        for (int32_t x = 0; x < xB; ++x)
        {
            io.outReal(f.GetU(x));
            d.real(f.GetUx(x));
            d.real(f.GetUxx(x));
            d.integer(f.GetMask(x));
        }
        for (int32_t x = 0; x < xB + 2; ++x) { d.integer(f.MaskAt(x)); }
        d.out(io);
    }
}

// Bounds 1..8, every data mode, masks, Dirichlet and Neumann borders, the
// scale types, spacings and time steps; the state after construction and
// after each of up to four Update calls.
ORACLE_CASE("PdeFilter1.heat")
{
    int32_t xB = DrawBound(io, 8);
    double dx = DrawSpacing(io);
    std::vector<double> data = DrawData(io, xB, 1, 1);
    bool hasMask = false;
    std::vector<int32_t> mask = DrawMask(io, xB, 1, 1, hasMask);
    double borderValue = DrawBorder(io);
    int32_t scaleType = DrawScaleType(io);
    std::unique_ptr<Heat1> f(MakeHeat1(xB, dx, data.data(), hasMask ? mask.data() : nullptr,
        borderValue, scaleType));
    double dt = DrawTimeStep(io, 0.5 * dx * dx);
    f->SetTimeStep(dt);
    io.outInt(f->GetQuantity());
    io.outInt(f->GetXBound());
    io.outReal(f->GetXSpacing());
    io.outReal(f->GetBorderValue());
    io.outInt(static_cast<int32_t>(f->GetScaleType()));
    io.outReal(f->Min());
    io.outReal(f->Offset());
    io.outReal(f->Scale());
    io.outReal(f->GetTimeStep());
    OutState1(io, *f);
    int32_t numSteps = io.integer(0, 4);
    for (int32_t k = 0; k < numSteps; ++k)
    {
        UpdateFilter(*f);
        OutState1(io, *f);
    }
}

// ---- PdeFilter2: GaussianBlur2, CurvatureFlow2, GradientAnisotropic2 ----

namespace
{
    // Exposes the protected PdeFilter2 state of a concrete filter; the
    // replay reads the same members of the port's classes.
    template <class Base>
    class Probe2 : public Base
    {
    public:
        template <typename... A>
        Probe2(A... args) : Base(args...) {}

        double Min() const { return this->mMin; }
        double Offset() const { return this->mOffset; }
        double Scale() const { return this->mScale; }
        int32_t Src() const { return this->mSrc; }
        bool HasMask() const { return this->mHasMask; }
        double Buf(int32_t k, int32_t x, int32_t y) const { return this->mBuffer[k][y][x]; }
        int32_t MaskAt(int32_t x, int32_t y) const { return this->mMask[y][x]; }
        // GradientAnisotropic2 only (instantiated on use).
        double Param() const { return this->mParameter; }
        double MHalfParam() const { return this->mMHalfParameter; }
    };

    using Blur2 = Probe2<GaussianBlur2<double>>;
    using Flow2 = Probe2<CurvatureFlow2<double>>;
    using Aniso2 = Probe2<GradientAnisotropic2<double>>;

    template <class F>
    void OutConstants2(oracle::Ctx& io, F const& f)
    {
        io.outInt(f.GetQuantity());
        io.outInt(f.GetXBound());
        io.outInt(f.GetYBound());
        io.outReal(f.GetXSpacing());
        io.outReal(f.GetYSpacing());
        io.outReal(f.GetBorderValue());
        io.outInt(static_cast<int32_t>(f.GetScaleType()));
        io.outReal(f.Min());
        io.outReal(f.Offset());
        io.outReal(f.Scale());
        io.outReal(f.GetTimeStep());
    }

    // Every public derived quantity at every pixel and the padded mask.
    template <class F>
    void DigestDerived2(oracle::Ctx& io, F const& f)
    {
        int32_t const xB = f.GetXBound(), yB = f.GetYBound();
        Digest d;
        for (int32_t y = 0; y < yB; ++y)
        {
            for (int32_t x = 0; x < xB; ++x)
            {
                d.real(f.GetU(x, y));
                d.real(f.GetUx(x, y));
                d.real(f.GetUy(x, y));
                d.real(f.GetUxx(x, y));
                d.real(f.GetUxy(x, y));
                d.real(f.GetUyy(x, y));
                d.integer(f.GetMask(x, y));
            }
        }
        for (int32_t y = 0; y < yB + 2; ++y)
        {
            for (int32_t x = 0; x < xB + 2; ++x) { d.integer(f.MaskAt(x, y)); }
        }
        d.out(io);
    }

    // Exact filters: the source index, both padded buffers (digest), every
    // derived quantity (digest), and the source image in full when it has at
    // most 9 pixels.
    template <class F>
    void OutState2(oracle::Ctx& io, F const& f)
    {
        int32_t const xB = f.GetXBound(), yB = f.GetYBound();
        io.outInt(f.Src());
        for (int32_t k = 0; k < 2; ++k)
        {
            Digest d;
            for (int32_t y = 0; y < yB + 2; ++y)
            {
                for (int32_t x = 0; x < xB + 2; ++x) { d.real(f.Buf(k, x, y)); }
            }
            d.out(io);
        }
        DigestDerived2(io, f);
        if (xB * yB <= 9)
        {
            for (int32_t y = 0; y < yB; ++y)
            {
                for (int32_t x = 0; x < xB; ++x) { io.outReal(f.GetU(x, y)); }
            }
        }
    }

    // GradientAnisotropic2: the cells an Update writes (interior pixels with
    // no mask or a nonzero mask) of the source buffer in full, compared with
    // a tolerance (exp); every other cell of both padded buffers is never
    // written after construction and goes into an exact digest.
    void OutStateAniso2(oracle::Ctx& io, Aniso2 const& f)
    {
        int32_t const xB = f.GetXBound(), yB = f.GetYBound();
        io.outInt(f.Src());
        Digest d;
        for (int32_t k = 0; k < 2; ++k)
        {
            for (int32_t y = 0; y < yB + 2; ++y)
            {
                for (int32_t x = 0; x < xB + 2; ++x)
                {
                    bool written = (1 <= x && x <= xB && 1 <= y && y <= yB
                        && (!f.HasMask() || f.MaskAt(x, y) != 0));
                    if (!written) { d.real(f.Buf(k, x, y)); }
                    else if (k == f.Src()) { io.outReal(f.Buf(k, x, y)); }
                }
            }
        }
        d.out(io);
        io.outReal(f.Param());
        io.outReal(f.MHalfParam());
    }

    double DrawK(oracle::Ctx& io)
    {
        switch (io.rawInteger(0, 5))
        {
        case 0: case 1: return io.given(1.0);
        case 2: return io.given(0.5);
        case 3: return io.given(2.0);
        case 4: return io.given(0.0);
        default: return io.given(io.raw(0.1, 3.0));
        }
    }

    // kind 0: GaussianBlur2, 1: CurvatureFlow2, 2: GradientAnisotropic2.
    template <class F>
    void Filter2Case(oracle::Ctx& io, int kind, int32_t maxB)
    {
        int32_t xB = DrawBound(io, maxB);
        int32_t yB = DrawBound(io, maxB);
        double dx = DrawSpacing(io);
        double dy = DrawSpacing(io);
        std::vector<double> data = DrawData(io, xB, yB, 1);
        bool hasMask = false;
        std::vector<int32_t> mask = DrawMask(io, xB, yB, 1, hasMask);
        double borderValue = DrawBorder(io);
        int32_t scaleType = DrawScaleType(io);
        int32_t const* maskPtr = (hasMask ? mask.data() : nullptr);
        std::unique_ptr<F> f;
        if constexpr (std::is_same<F, Aniso2>::value)
        {
            double K = DrawK(io);
            f = std::make_unique<F>(xB, yB, dx, dy, data.data(), maskPtr, borderValue,
                static_cast<ST>(scaleType), K);
        }
        else
        {
            f = std::make_unique<F>(xB, yB, dx, dy, data.data(), maskPtr, borderValue,
                static_cast<ST>(scaleType));
        }
        double stable = 0.125 * std::min(dx * dx, dy * dy);
        if constexpr (std::is_same<F, Blur2>::value)
        {
            io.outReal(f->GetMaximumTimeStep());
            stable = f->GetMaximumTimeStep();
        }
        double dt = DrawTimeStep(io, stable);
        f->SetTimeStep(dt);
        OutConstants2(io, *f);
        DigestDerived2(io, *f);
        int32_t numSteps = io.integer(0, 4);
        for (int32_t k = 0; k <= numSteps; ++k)
        {
            if (k > 0) { UpdateFilter(*f); }
            if constexpr (std::is_same<F, Aniso2>::value) { OutStateAniso2(io, *f); }
            else { OutState2(io, *f); }
        }
        (void)kind;
    }
}

// Bounds 1..6 (6x6 often), every data mode, masks, Dirichlet and Neumann
// borders, scale types, spacings and time steps (including
// GetMaximumTimeStep); the state after construction and after each of up to
// four Update calls.
ORACLE_CASE("GaussianBlur2.update") { Filter2Case<Blur2>(io, 0, 6); }

// As GaussianBlur2.update. The x-only data mode puts the uy == 0 and
// denom == 0 branches in every record of that mode.
ORACLE_CASE("CurvatureFlow2.update") { Filter2Case<Flow2>(io, 1, 6); }

// Bounds 1..5; the updated cells are compared with the default tolerance
// 1e-12 (std::exp in the four conductances), everything else exactly; K
// includes 0 (mParameter infinite) and the constant image gives the
// preserved #122 infinity/NaN.
ORACLE_CASE("GradientAnisotropic2.update") { Filter2Case<Aniso2>(io, 2, 5); }

// ---- PdeFilter3: GaussianBlur3, CurvatureFlow3, GradientAnisotropic3 ----

namespace
{
    template <class Base>
    class Probe3 : public Base
    {
    public:
        template <typename... A>
        Probe3(A... args) : Base(args...) {}

        double Min() const { return this->mMin; }
        double Offset() const { return this->mOffset; }
        double Scale() const { return this->mScale; }
        int32_t Src() const { return this->mSrc; }
        bool HasMask() const { return this->mHasMask; }
        double Buf(int32_t k, int32_t x, int32_t y, int32_t z) const
        {
            return this->mBuffer[k][z][y][x];
        }
        int32_t MaskAt(int32_t x, int32_t y, int32_t z) const { return this->mMask[z][y][x]; }
        double Param() const { return this->mParameter; }
        double MHalfParam() const { return this->mMHalfParameter; }
    };

    using Blur3 = Probe3<GaussianBlur3<double>>;
    using Flow3 = Probe3<CurvatureFlow3<double>>;
    using Aniso3 = Probe3<GradientAnisotropic3<double>>;

    template <class F>
    void OutConstants3(oracle::Ctx& io, F const& f)
    {
        io.outInt(f.GetQuantity());
        io.outInt(f.GetXBound());
        io.outInt(f.GetYBound());
        io.outInt(f.GetZBound());
        io.outReal(f.GetXSpacing());
        io.outReal(f.GetYSpacing());
        io.outReal(f.GetZSpacing());
        io.outReal(f.GetBorderValue());
        io.outInt(static_cast<int32_t>(f.GetScaleType()));
        io.outReal(f.Min());
        io.outReal(f.Offset());
        io.outReal(f.Scale());
        io.outReal(f.GetTimeStep());
    }

    template <class F>
    void DigestDerived3(oracle::Ctx& io, F const& f)
    {
        int32_t const xB = f.GetXBound(), yB = f.GetYBound(), zB = f.GetZBound();
        Digest d;
        for (int32_t z = 0; z < zB; ++z)
        {
            for (int32_t y = 0; y < yB; ++y)
            {
                for (int32_t x = 0; x < xB; ++x)
                {
                    d.real(f.GetU(x, y, z));
                    d.real(f.GetUx(x, y, z));
                    d.real(f.GetUy(x, y, z));
                    d.real(f.GetUz(x, y, z));
                    d.real(f.GetUxx(x, y, z));
                    d.real(f.GetUxy(x, y, z));
                    d.real(f.GetUxz(x, y, z));
                    d.real(f.GetUyy(x, y, z));
                    d.real(f.GetUyz(x, y, z));
                    d.real(f.GetUzz(x, y, z));
                    d.integer(f.GetMask(x, y, z));
                }
            }
        }
        for (int32_t z = 0; z < zB + 2; ++z)
        {
            for (int32_t y = 0; y < yB + 2; ++y)
            {
                for (int32_t x = 0; x < xB + 2; ++x) { d.integer(f.MaskAt(x, y, z)); }
            }
        }
        d.out(io);
    }

    template <class F>
    void OutState3(oracle::Ctx& io, F const& f)
    {
        int32_t const xB = f.GetXBound(), yB = f.GetYBound(), zB = f.GetZBound();
        io.outInt(f.Src());
        for (int32_t k = 0; k < 2; ++k)
        {
            Digest d;
            for (int32_t z = 0; z < zB + 2; ++z)
            {
                for (int32_t y = 0; y < yB + 2; ++y)
                {
                    for (int32_t x = 0; x < xB + 2; ++x) { d.real(f.Buf(k, x, y, z)); }
                }
            }
            d.out(io);
        }
        DigestDerived3(io, f);
        if (xB * yB * zB <= 8)
        {
            for (int32_t z = 0; z < zB; ++z)
            {
                for (int32_t y = 0; y < yB; ++y)
                {
                    for (int32_t x = 0; x < xB; ++x) { io.outReal(f.GetU(x, y, z)); }
                }
            }
        }
    }

    void OutStateAniso3(oracle::Ctx& io, Aniso3 const& f)
    {
        int32_t const xB = f.GetXBound(), yB = f.GetYBound(), zB = f.GetZBound();
        io.outInt(f.Src());
        Digest d;
        for (int32_t k = 0; k < 2; ++k)
        {
            for (int32_t z = 0; z < zB + 2; ++z)
            {
                for (int32_t y = 0; y < yB + 2; ++y)
                {
                    for (int32_t x = 0; x < xB + 2; ++x)
                    {
                        bool written = (1 <= x && x <= xB && 1 <= y && y <= yB
                            && 1 <= z && z <= zB && (!f.HasMask() || f.MaskAt(x, y, z) != 0));
                        if (!written) { d.real(f.Buf(k, x, y, z)); }
                        else if (k == f.Src()) { io.outReal(f.Buf(k, x, y, z)); }
                    }
                }
            }
        }
        d.out(io);
        io.outReal(f.Param());
        io.outReal(f.MHalfParam());
    }

    template <class F>
    void Filter3Case(oracle::Ctx& io, int32_t maxB)
    {
        int32_t xB = DrawBound(io, maxB);
        int32_t yB = DrawBound(io, maxB);
        int32_t zB = DrawBound(io, maxB);
        double dx = DrawSpacing(io);
        double dy = DrawSpacing(io);
        double dz = DrawSpacing(io);
        std::vector<double> data = DrawData(io, xB, yB, zB);
        bool hasMask = false;
        std::vector<int32_t> mask = DrawMask(io, xB, yB, zB, hasMask);
        double borderValue = DrawBorder(io);
        int32_t scaleType = DrawScaleType(io);
        int32_t const* maskPtr = (hasMask ? mask.data() : nullptr);
        std::unique_ptr<F> f;
        if constexpr (std::is_same<F, Aniso3>::value)
        {
            double K = DrawK(io);
            f = std::make_unique<F>(xB, yB, zB, dx, dy, dz, data.data(), maskPtr, borderValue,
                static_cast<ST>(scaleType), K);
        }
        else
        {
            f = std::make_unique<F>(xB, yB, zB, dx, dy, dz, data.data(), maskPtr, borderValue,
                static_cast<ST>(scaleType));
        }
        double stable = 0.125 * std::min(std::min(dx * dx, dy * dy), dz * dz);
        if constexpr (std::is_same<F, Blur3>::value)
        {
            io.outReal(f->GetMaximumTimeStep());
            stable = f->GetMaximumTimeStep();
        }
        double dt = DrawTimeStep(io, stable);
        f->SetTimeStep(dt);
        OutConstants3(io, *f);
        DigestDerived3(io, *f);
        int32_t numSteps = io.integer(0, 3);
        for (int32_t k = 0; k <= numSteps; ++k)
        {
            if (k > 0) { UpdateFilter(*f); }
            if constexpr (std::is_same<F, Aniso3>::value) { OutStateAniso3(io, *f); }
            else { OutState3(io, *f); }
        }
    }
}

// Bounds 1..5 (5x5x5 often); as GaussianBlur2.update, up to three steps.
ORACLE_CASE("GaussianBlur3.update") { Filter3Case<Blur3>(io, 5); }

ORACLE_CASE("CurvatureFlow3.update") { Filter3Case<Flow3>(io, 5); }

// Bounds 1..3; tolerance 1e-12 on the updated cells (std::exp in the six
// conductances), exact elsewhere.
ORACLE_CASE("GradientAnisotropic3.update") { Filter3Case<Aniso3>(io, 3); }

// ---- FastMarch2, FastMarch3 ----

namespace
{
    // Exposes the protected heap and speeds; the replay reads the same
    // members of the port's classes. Upstream Iterate removes from the heap
    // unconditionally (on an empty heap it then reads mTrials[SIZE_MAX]), so
    // the cases call it only while NumTrials() > 0.
    template <class Base>
    class MarchProbe : public Base
    {
    public:
        template <typename... A>
        MarchProbe(A... args) : Base(args...) {}

        int32_t NumTrials() const { return this->mHeap.GetNumElements(); }
        bool PeekMin(size_t& key, double& value) const { return this->mHeap.GetMinimum(key, value); }
        double InvSpeed(size_t i) const { return this->mInvSpeeds[i]; }
        typename MinHeap<size_t, double>::Record const* Trial(size_t i) const
        {
            return this->mTrials[i];
        }
    };

    using March2 = MarchProbe<FastMarch2<double>>;
    using March3 = MarchProbe<FastMarch3<double>>;

    // IsValid, IsTrial, IsFar, IsZeroSpeed, IsInterior, IsBoundary as bits.
    template <class M>
    int32_t Flags(M const& m, size_t i)
    {
        return (m.IsValid(i) ? 1 : 0) | (m.IsTrial(i) ? 2 : 0) | (m.IsFar(i) ? 4 : 0)
            | (m.IsZeroSpeed(i) ? 8 : 0) | (m.IsInterior(i) ? 16 : 0) | (m.IsBoundary(i) ? 32 : 0);
    }

    // The number of trials, digests of the times, the flags and the inverse
    // speeds of the included pixels, the trial records (key, value, heap
    // position) as count + digest in key order, GetInterior and GetBoundary
    // as count + digest, and GetTimeExtremes.
    template <class M>
    void OutMarchState(oracle::Ctx& io, M const& m, std::vector<char> const& include)
    {
        size_t const q = m.GetQuantity();
        io.outInt(m.NumTrials());
        Digest times, flags, speeds, trials;
        int32_t numTrials = 0;
        for (size_t i = 0; i < q; ++i)
        {
            if (!include[i]) { continue; }
            times.real(m.GetTime(i));
            flags.integer(Flags(m, i));
            speeds.real(m.InvSpeed(i));
            auto const* record = m.Trial(i);
            if (record)
            {
                ++numTrials;
                trials.integer(static_cast<int64_t>(record->key));
                trials.real(record->value);
                trials.integer(record->index);
            }
        }
        times.out(io);
        flags.out(io);
        speeds.out(io);
        io.outInt(numTrials);
        trials.out(io);
        std::vector<size_t> interior, boundary;
        m.GetInterior(interior);
        m.GetBoundary(boundary);
        Digest di, db;
        for (auto v : interior) { di.integer(static_cast<int64_t>(v)); }
        for (auto v : boundary) { db.integer(static_cast<int64_t>(v)); }
        io.outInt(interior.size());
        di.out(io);
        io.outInt(boundary.size());
        db.out(io);
        double minValue = 0.0, maxValue = 0.0;
        m.GetTimeExtremes(minValue, maxValue);
        io.outReal(minValue);
        io.outReal(maxValue);
    }

    // Removes up to maxIter minima while the heap is not empty (and, when
    // 'stop' says so, before the march would reach an upstream face voxel),
    // emitting each removed key and value.
    template <class M, class Stop>
    void March(oracle::Ctx& io, M& m, int32_t maxIter, Stop stop)
    {
        for (int32_t k = 0; k < maxIter; ++k)
        {
            size_t key = 0;
            double value = 0.0;
            if (m.NumTrials() == 0 || !m.PeekMin(key, value) || stop(key))
            {
                break;
            }
            io.outInt(key);
            io.outReal(value);
            m.Iterate();
        }
        io.outInt(-1);
    }

    // One recorded double per pixel: lattice and uniform positive speeds,
    // 0, -1 and NaN (all three mark the pixel zero speed), or with
    // 'contrast' a mix of 0.25 and 4 that drives the negative-discriminant
    // fallback (#439).
    double DrawPixelSpeed(oracle::Ctx& io, bool contrast)
    {
        if (contrast)
        {
            return io.given(io.rawInteger(0, 1) == 0 ? 0.25 : 4.0);
        }
        int k = io.rawInteger(0, 11);
        double v = (k == 0 ? 0.0 : (k == 1 ? -1.0 : (k == 2 ? std::nan("")
            : (k < 7 ? io.rawInteger(1, 4) * 0.5 : io.raw(0.1, 4.0)))));
        return io.given(v);
    }

    // One recorded double: lattice and uniform positive, 0 (infinite inverse
    // speed), negative and, rarely, NaN.
    double DrawConstantSpeed(oracle::Ctx& io)
    {
        int k = io.rawInteger(0, 15);
        double v = (k == 0 ? 0.0 : (k == 1 ? -2.0 : (k == 2 ? std::nan("")
            : (k < 9 ? io.rawInteger(1, 4) * 0.5 : io.raw(0.1, 4.0)))));
        return io.given(v);
    }
}

namespace
{
    // Grid 1..7 squared, spacings (stored, unused by the method: #439),
    // 0..4 seeds anywhere (duplicates allowed; border seeds are overwritten
    // by the zero-speed border), and a constant speed, per-pixel speeds or
    // the high-contrast per-pixel mix. Every 10th record (repro) is the
    // #439 reproduction: 6x4, speed 1/2, seeds (4,2) and (1,2).
    std::unique_ptr<March2> DrawMarch2(oracle::Ctx& io, bool allowRepro)
    {
        bool repro = allowRepro && (io.index() % 10 == 9);
        int32_t xB = static_cast<int32_t>(io.given(repro ? 6 : io.rawInteger(1, 7)));
        int32_t yB = static_cast<int32_t>(io.given(repro ? 4 : io.rawInteger(1, 7)));
        int32_t q = xB * yB;
        double dx = DrawSpacing(io);
        double dy = DrawSpacing(io);
        int32_t numSeeds = static_cast<int32_t>(io.given(repro ? 2 : io.rawInteger(0, 4)));
        std::vector<size_t> seeds;
        for (int32_t k = 0; k < numSeeds; ++k)
        {
            seeds.push_back(static_cast<size_t>(io.given(repro ? (k == 0 ? 16 : 13)
                : io.rawInteger(0, q - 1))));
        }
        int32_t speedMode = static_cast<int32_t>(io.given(repro ? 0 : io.rawInteger(0, 2)));
        if (speedMode == 0)
        {
            double speed = (repro ? io.given(0.5) : DrawConstantSpeed(io));
            return std::make_unique<March2>(static_cast<size_t>(xB), static_cast<size_t>(yB),
                dx, dy, seeds, speed);
        }
        std::vector<double> speeds(q);
        for (int32_t i = 0; i < q; ++i) { speeds[i] = DrawPixelSpeed(io, speedMode == 2); }
        return std::make_unique<March2>(static_cast<size_t>(xB), static_cast<size_t>(yB),
            dx, dy, seeds, speeds);
    }

    void OutMarch2Constants(oracle::Ctx& io, March2 const& m)
    {
        io.outInt(m.GetQuantity());
        io.outInt(m.GetXBound());
        io.outInt(m.GetYBound());
        io.outReal(m.GetXSpacing());
        io.outReal(m.GetYSpacing());
        io.outInt(m.Index(m.GetXBound() - 1, m.GetYBound() - 1));
    }
}

// The constructor state (every classifier, trial records with their heap
// positions, interior and boundary sets, time extremes), then up to 2q
// removals (key and value of each, so ties in the heap are compared), then
// the final state. The negative-discriminant fallback (#439, preserved) is
// reached by the repro records and the high-contrast speeds; the replay
// counts the ComputeTime branches.
ORACLE_CASE("FastMarch2.march")
{
    std::unique_ptr<March2> m = DrawMarch2(io, true);
    std::vector<char> all(m->GetQuantity(), 1);
    OutMarch2Constants(io, *m);
    OutMarchState(io, *m, all);
    int32_t maxIter = io.integer(0, 2 * static_cast<int32_t>(m->GetQuantity()));
    March(io, *m, maxIter, [](size_t) { return false; });
    OutMarchState(io, *m, all);
}

// SetTime on interior pixels with arbitrary values (both zeros, +-max,
// infinities, NaN, negatives, uniform), the state, a few removals, the
// state again. Border pixels are not set: a valid or far border pixel makes
// upstream's IsBoundary/Iterate index outside the grid.
ORACLE_CASE("FastMarch2.accessors")
{
    std::unique_ptr<March2> m = DrawMarch2(io, false);
    int32_t const xB = static_cast<int32_t>(m->GetXBound());
    int32_t const yB = static_cast<int32_t>(m->GetYBound());
    bool const hasInterior = (xB >= 3 && yB >= 3);
    int32_t numSets = static_cast<int32_t>(io.given(hasInterior ? io.rawInteger(0, 6) : 0));
    for (int32_t k = 0; k < numSets; ++k)
    {
        int32_t x = io.integer(1, xB - 2);
        int32_t y = io.integer(1, yB - 2);
        double const inf = std::numeric_limits<double>::infinity();
        double const values[9] = { 0.0, -0.0, kMax, -kMax, inf, -inf, std::nan(""), -1.5, 2.5 };
        int32_t which = io.rawInteger(0, 11);
        double t = io.given(which < 9 ? values[which] : io.raw(0.0, 10.0));
        m->SetTime(m->Index(x, y), t);
    }
    std::vector<char> all(m->GetQuantity(), 1);
    OutMarchState(io, *m, all);
    int32_t maxIter = io.integer(0, 8);
    March(io, *m, maxIter, [](size_t) { return false; });
    OutMarchState(io, *m, all);
}

namespace
{
    // The number of coordinates of voxel i on the grid boundary. Upstream
    // FastMarch3 marks the voxels with two or three (edges, vertices) zero
    // speed and leaves the faces (exactly one) alone; the port marks all
    // (#121).
    int32_t NumExtremes(size_t i, int32_t xB, int32_t yB, int32_t zB)
    {
        int32_t x = static_cast<int32_t>(i % xB);
        int32_t y = static_cast<int32_t>((i / xB) % yB);
        int32_t z = static_cast<int32_t>(i / (static_cast<size_t>(xB) * yB));
        return (x == 0 || x == xB - 1 ? 1 : 0) + (y == 0 || y == yB - 1 ? 1 : 0)
            + (z == 0 || z == zB - 1 ? 1 : 0);
    }

    bool IsZFace(size_t i, int32_t xB, int32_t yB, int32_t zB)
    {
        int32_t z = static_cast<int32_t>(i / (static_cast<size_t>(xB) * yB));
        return NumExtremes(i, xB, yB, zB) == 1 && (z == 0 || z == zB - 1);
    }

    void OutMarch3Constants(oracle::Ctx& io, March3 const& m)
    {
        io.outInt(m.GetQuantity());
        io.outInt(m.GetXBound());
        io.outInt(m.GetYBound());
        io.outInt(m.GetZBound());
        io.outReal(m.GetXSpacing());
        io.outReal(m.GetYSpacing());
        io.outReal(m.GetZSpacing());
        io.outInt(m.Index(m.GetXBound() - 1, m.GetYBound() - 1, m.GetZBound() - 1));
    }

    // Bounds lo..5; seeds anywhere or (interiorSeeds) strictly inside.
    std::vector<size_t> DrawSeeds3(oracle::Ctx& io, int32_t xB, int32_t yB, int32_t zB,
        int32_t minSeeds, bool interiorSeeds)
    {
        int32_t numSeeds = io.integer(minSeeds, 4);
        std::vector<size_t> seeds;
        for (int32_t k = 0; k < numSeeds; ++k)
        {
            if (interiorSeeds)
            {
                int32_t x = io.integer(1, xB - 2);
                int32_t y = io.integer(1, yB - 2);
                int32_t z = io.integer(1, zB - 2);
                seeds.push_back(static_cast<size_t>(x + xB * (y + yB * z)));
            }
            else
            {
                seeds.push_back(static_cast<size_t>(io.integer(0, xB * yB * zB - 1)));
            }
        }
        return seeds;
    }
}

// Per-voxel speeds with every face voxel (exactly one coordinate on the
// boundary) given a nonpositive speed, so upstream's missing face marking
// (#121) is invisible and the whole state is comparable; bounds 1..5,
// seeds anywhere, the three per-voxel speed mixes; up to 2q removals.
ORACLE_CASE("FastMarch3.march")
{
    int32_t xB = io.integer(1, 5);
    int32_t yB = io.integer(1, 5);
    int32_t zB = io.integer(1, 5);
    int32_t q = xB * yB * zB;
    double dx = DrawSpacing(io);
    double dy = DrawSpacing(io);
    double dz = DrawSpacing(io);
    std::vector<size_t> seeds = DrawSeeds3(io, xB, yB, zB, 0, false);
    bool contrast = io.boolean();
    std::vector<double> speeds(q);
    for (int32_t i = 0; i < q; ++i)
    {
        if (NumExtremes(i, xB, yB, zB) == 1)
        {
            int k = io.rawInteger(0, 2);
            speeds[i] = io.given(k == 0 ? 0.0 : (k == 1 ? -1.0 : std::nan("")));
        }
        else
        {
            speeds[i] = DrawPixelSpeed(io, contrast);
        }
    }
    auto m = std::make_unique<March3>(static_cast<size_t>(xB), static_cast<size_t>(yB),
        static_cast<size_t>(zB), dx, dy, dz, seeds, speeds);
    std::vector<char> all(q, 1);
    OutMarch3Constants(io, *m);
    OutMarchState(io, *m, all);
    int32_t maxIter = io.integer(0, 2 * q);
    March(io, *m, maxIter, [](size_t) { return false; });
    OutMarchState(io, *m, all);
}

// The constant-speed constructor, where upstream leaves the face voxels far
// (#121). Seeds strictly inside; the march stops before removing a voxel
// with a face neighbour (upstream would then compute a face voxel's time
// and index outside the grid for a z face), and every output leaves the
// face voxels out (their time, flags and inverse speed are the #121
// difference). Bounds 3..5.
ORACLE_CASE("FastMarch3.march.constantSpeed")
{
    int32_t xB = io.integer(3, 5);
    int32_t yB = io.integer(3, 5);
    int32_t zB = io.integer(3, 5);
    int32_t q = xB * yB * zB;
    double dx = DrawSpacing(io);
    double dy = DrawSpacing(io);
    double dz = DrawSpacing(io);
    std::vector<size_t> seeds = DrawSeeds3(io, xB, yB, zB, 0, true);
    double speed = DrawConstantSpeed(io);
    auto m = std::make_unique<March3>(static_cast<size_t>(xB), static_cast<size_t>(yB),
        static_cast<size_t>(zB), dx, dy, dz, seeds, speed);
    std::vector<char> include(q, 1);
    for (int32_t i = 0; i < q; ++i) { include[i] = (NumExtremes(i, xB, yB, zB) != 1); }
    OutMarch3Constants(io, *m);
    OutMarchState(io, *m, include);
    auto hasFaceNeighbor = [&](size_t key)
    {
        size_t const n[6] = { key - 1, key + 1, key - xB, key + xB,
            key - static_cast<size_t>(xB) * yB, key + static_cast<size_t>(xB) * yB };
        for (size_t j : n) { if (NumExtremes(j, xB, yB, zB) == 1) { return true; } }
        return false;
    };
    int32_t maxIter = io.integer(0, 2 * q);
    March(io, *m, maxIter, hasFaceNeighbor);
    OutMarchState(io, *m, include);
}

// Finding #121 (fixed in the port): per-voxel speeds with the z faces at
// zero speed and the x and y faces positive. Upstream leaves the x/y face
// voxels far, marches into them (their wrapped neighbours i - 1, i + xB ...
// stay inside the grid because the z faces are never visited, so every
// record is deterministic) and assigns them times; the port has them zero
// speed from the start. Every record deviates. Bounds 3..5, 1..4 interior
// seeds, positive speeds elsewhere, the march runs to completion.
ORACLE_CASE("FastMarch3.faces")
{
    int32_t xB = io.integer(3, 5);
    int32_t yB = io.integer(3, 5);
    int32_t zB = io.integer(3, 5);
    int32_t q = xB * yB * zB;
    std::vector<size_t> seeds = DrawSeeds3(io, xB, yB, zB, 1, true);
    std::vector<double> speeds(q);
    for (int32_t i = 0; i < q; ++i)
    {
        double v = (IsZFace(i, xB, yB, zB) ? 0.0
            : (io.rawInteger(0, 1) == 0 ? io.rawInteger(1, 4) * 0.5 : io.raw(0.1, 4.0)));
        speeds[i] = io.given(v);
    }
    auto m = std::make_unique<March3>(static_cast<size_t>(xB), static_cast<size_t>(yB),
        static_cast<size_t>(zB), 1.0, 1.0, 1.0, seeds, speeds);
    std::vector<char> all(q, 1);
    OutMarchState(io, *m, all);
    March(io, *m, 2 * q, [](size_t) { return false; });
    OutMarchState(io, *m, all);
}

// ---- SurfaceExtractor ----

namespace
{
    using SE = SurfaceExtractor<int32_t, double>;

    // SurfaceExtractor is abstract (its extractors are group 26). This
    // subclass, defined identically in the replay, returns recorded rational
    // vertices and triangles from Extract and an affine gradient
    // (g0 + g1 x, g2 + g3 y, g4 + g5 z) from GetGradient, which exposes the
    // base's own code: Vertex and Triangle, MakeUnique, Convert, the
    // public Extract overload, OrientTriangles and ComputeNormals.
    class Recorded : public SE
    {
    public:
        Recorded(int32_t xB, int32_t yB, int32_t zB, int32_t const* voxels)
            :
            SE(xB, yB, zB, voxels)
        {
        }

        using SE::Extract;

        virtual void Extract(int32_t, std::vector<Vertex>& vertices,
            std::vector<Triangle>& triangles) override
        {
            vertices = rv;
            triangles = rt;
        }

        std::vector<Vertex> rv;
        std::vector<Triangle> rt;
        std::array<double, 6> g{};

    protected:
        virtual std::array<double, 3> GetGradient(std::array<double, 3> const& p) override
        {
            return { g[0] + g[1] * p[0], g[2] + g[3] * p[1], g[4] + g[5] * p[2] };
        }
    };

    int32_t const kVoxels[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

    // One rational coordinate as (numer, denom), both recorded: a base value
    // a/b (a in -2..4, b in 1..2) scaled by m in 1..3 and a sign, so equal
    // rationals come in several representations, negative denominators are
    // normalized, and 0 over a negative denominator probes the sign of the
    // normalized zero numerator.
    void DrawRational(oracle::Ctx& io, int64_t& numer, int64_t& denom)
    {
        int64_t a = io.rawInteger(-2, 4), b = io.rawInteger(1, 2);
        int64_t m = io.rawInteger(1, 3), s = (io.rawInteger(0, 2) == 0 ? -1 : 1);
        numer = static_cast<int64_t>(io.given(static_cast<double>(s * m * a)));
        denom = static_cast<int64_t>(io.given(static_cast<double>(s * m * b)));
    }

    // 1..12 vertices (a third of them re-represent an earlier one) and
    // 0..10 triangles over them (repeated indices allowed).
    void DrawSurface(oracle::Ctx& io, Recorded& r)
    {
        int32_t numVertices = io.integer(1, 12);
        for (int32_t v = 0; v < numVertices; ++v)
        {
            int64_t n[3], d[3];
            bool repeat = (v > 0 && io.rawInteger(0, 2) == 0);
            SE::Vertex const* prev = (repeat ? &r.rv[io.rawInteger(0, v - 1)] : nullptr);
            for (int32_t c = 0; c < 3; ++c)
            {
                if (prev)
                {
                    int64_t pn = (c == 0 ? prev->xNumer : (c == 1 ? prev->yNumer : prev->zNumer));
                    int64_t pd = (c == 0 ? prev->xDenom : (c == 1 ? prev->yDenom : prev->zDenom));
                    int64_t m = io.rawInteger(1, 2) * (io.rawInteger(0, 1) == 0 ? -1 : 1);
                    n[c] = static_cast<int64_t>(io.given(static_cast<double>(m * pn)));
                    d[c] = static_cast<int64_t>(io.given(static_cast<double>(m * pd)));
                }
                else
                {
                    DrawRational(io, n[c], d[c]);
                }
            }
            r.rv.push_back(SE::Vertex(n[0], d[0], n[1], d[1], n[2], d[2]));
        }
        int32_t numTriangles = io.integer(0, 10);
        for (int32_t t = 0; t < numTriangles; ++t)
        {
            int32_t v0 = io.integer(0, numVertices - 1);
            int32_t v1 = io.integer(0, numVertices - 1);
            int32_t v2 = io.integer(0, numVertices - 1);
            r.rt.push_back(SE::Triangle(v0, v1, v2));
        }
    }

    void OutRationals(oracle::Ctx& io, std::vector<SE::Vertex> const& vertices)
    {
        io.outInt(vertices.size());
        for (auto const& v : vertices)
        {
            io.outInt(v.xNumer); io.outInt(v.xDenom);
            io.outInt(v.yNumer); io.outInt(v.yDenom);
            io.outInt(v.zNumer); io.outInt(v.zDenom);
        }
    }

    void OutTriangles(oracle::Ctx& io, std::vector<SE::Triangle> const& triangles)
    {
        io.outInt(triangles.size());
        for (auto const& t : triangles) { io.outInt(t.v[0]); io.outInt(t.v[1]); io.outInt(t.v[2]); }
    }

    void OutPoints(oracle::Ctx& io, std::vector<std::array<double, 3>> const& points)
    {
        io.outInt(points.size());
        for (auto const& p : points) { io.outReal(p[0]); io.outReal(p[1]); io.outReal(p[2]); }
    }

    __declspec(noinline) void CallMakeUnique(Recorded& r, std::vector<SE::Vertex>& v,
        std::vector<SE::Triangle>& t)
    {
        r.MakeUnique(v, t);
    }

    __declspec(noinline) void CallConvert(Recorded& r, std::vector<SE::Vertex> const& v,
        std::vector<std::array<double, 3>>& out)
    {
        r.Convert(v, out);
    }

    __declspec(noinline) void CallExtract(Recorded& r, bool unique,
        std::vector<std::array<double, 3>>& v, std::vector<SE::Triangle>& t)
    {
        r.Extract(0, unique, v, t);
    }

    __declspec(noinline) void CallOrient(Recorded& r, std::vector<std::array<double, 3>>& v,
        std::vector<SE::Triangle>& t, bool sameDir)
    {
        r.OrientTriangles(v, t, sameDir);
    }

    __declspec(noinline) void CallNormals(Recorded& r, std::vector<std::array<double, 3>> const& v,
        std::vector<SE::Triangle> const& t, std::vector<std::array<double, 3>>& n)
    {
        r.ComputeNormals(v, t, n);
    }
}

// The Vertex constructor's sign normalization (fields emitted), operator==
// and operator< on each consecutive pair in both orders, and the Triangle
// constructor's rotation, operator== and operator< likewise.
ORACLE_CASE("SurfaceExtractor.vertexTriangle")
{
    Recorded r(2, 2, 2, kVoxels);
    DrawSurface(io, r);
    OutRationals(io, r.rv);
    for (size_t k = 0; k + 1 < r.rv.size(); ++k)
    {
        io.outBool(r.rv[k] == r.rv[k + 1]);
        io.outBool(r.rv[k] < r.rv[k + 1]);
        io.outBool(r.rv[k + 1] < r.rv[k]);
    }
    OutTriangles(io, r.rt);
    for (size_t k = 0; k + 1 < r.rt.size(); ++k)
    {
        io.outBool(r.rt[k] == r.rt[k + 1]);
        io.outBool(r.rt[k] < r.rt[k + 1]);
        io.outBool(r.rt[k + 1] < r.rt[k]);
    }
}

// MakeUnique on the recorded lists (the rational vertices and index triples
// it keeps, in upstream's first-occurrence numbering), then Convert.
ORACLE_CASE("SurfaceExtractor.makeUnique")
{
    Recorded r(2, 2, 2, kVoxels);
    DrawSurface(io, r);
    std::vector<SE::Vertex> v = r.rv;
    std::vector<SE::Triangle> t = r.rt;
    CallMakeUnique(r, v, t);
    OutRationals(io, v);
    OutTriangles(io, t);
    std::vector<std::array<double, 3>> points;
    CallConvert(r, v, points);
    OutPoints(io, points);
}

// The public Extract(level, removeDuplicateVertices, ...) overload (through
// the subclass's Extract), OrientTriangles with both sameDir values on an
// affine gradient (lattice, dyadic or uniform coefficients), and
// ComputeNormals.
ORACLE_CASE("SurfaceExtractor.extract")
{
    Recorded r(2, 2, 2, kVoxels);
    DrawSurface(io, r);
    bool unique = io.boolean();
    int gmode = io.rawInteger(0, 2);
    for (int32_t i = 0; i < 6; ++i)
    {
        r.g[i] = io.given(gmode == 0 ? io.rawInteger(-2, 2)
            : (gmode == 1 ? 0.25 * io.rawInteger(-8, 8) : io.raw(-2.0, 2.0)));
    }
    bool sameDir = io.boolean();
    std::vector<std::array<double, 3>> vertices;
    std::vector<SE::Triangle> triangles;
    CallExtract(r, unique, vertices, triangles);
    OutPoints(io, vertices);
    OutTriangles(io, triangles);
    CallOrient(r, vertices, triangles, sameDir);
    OutTriangles(io, triangles);
    std::vector<std::array<double, 3>> normals;
    CallNormals(r, vertices, triangles, normals);
    OutPoints(io, normals);
}

// The constructor's LogAssert (every bound must exceed 1); bounds 0..3, so
// about half the records throw on both sides.
ORACLE_CASE("SurfaceExtractor.invalid")
{
    int32_t xB = io.integer(0, 3);
    int32_t yB = io.integer(1, 3);
    int32_t zB = io.integer(1, 3);
    std::vector<int32_t> voxels(27, 0);
    Recorded r(xB, yB, zB, voxels.data());
    io.outInt(1);
}

// ---- TetrahedraRasterizer ----

namespace
{
    using V3 = std::array<double, 3>;
    using T4 = std::array<size_t, 4>;

    // The rasterized grid: its size, the number of voxels set, a digest of
    // the (index, tetrahedron) pairs in index order, and the pairs in full
    // when there are at most 40.
    void OutGrid(oracle::Ctx& io, std::vector<int32_t> const& grid)
    {
        io.outInt(grid.size());
        Digest d;
        int32_t count = 0;
        for (size_t j = 0; j < grid.size(); ++j)
        {
            if (grid[j] != -1) { ++count; d.integer(static_cast<int64_t>(j)); d.integer(grid[j]); }
        }
        io.outInt(count);
        d.out(io);
        if (count <= 40)
        {
            for (size_t j = 0; j < grid.size(); ++j)
            {
                if (grid[j] != -1) { io.outInt(j); io.outInt(grid[j]); }
            }
        }
    }

    // Lattice tetrahedra in [lo, hi]^3 (raw draws; RecordMesh records the
    // result): random (either orientation), random with a flat fourth
    // vertex, or a lattice box split into the six Kuhn tetrahedra (shared
    // faces and the box diagonal put voxel centres exactly on faces).
    void DrawLatticeMesh(oracle::Ctx& io, int32_t lo, int32_t hi,
        std::vector<V3>& vertices, std::vector<T4>& tetra)
    {
        int32_t mode = io.rawInteger(0, 2);
        if (mode == 2)
        {
            int32_t a[3], s[3];
            for (int32_t i = 0; i < 3; ++i)
            {
                a[i] = io.rawInteger(lo, hi - 1);
                s[i] = io.rawInteger(1, hi - a[i]);
            }
            for (int32_t k = 0; k < 8; ++k)
            {
                vertices.push_back({ static_cast<double>(a[0] + (k & 1 ? s[0] : 0)),
                    static_cast<double>(a[1] + (k & 2 ? s[1] : 0)),
                    static_cast<double>(a[2] + (k & 4 ? s[2] : 0)) });
            }
            // Kuhn: the paths 0 -> 7 adding one axis bit at a time.
            size_t const perms[6][2] = { {1,2}, {1,4}, {2,1}, {2,4}, {4,1}, {4,2} };
            for (auto const& p : perms)
            {
                tetra.push_back({ 0, p[0], p[0] | p[1], 7 });
            }
            return;
        }
        int32_t numTetra = io.rawInteger(1, 5);
        for (int32_t t = 0; t < numTetra; ++t)
        {
            for (int32_t k = 0; k < 4; ++k)
            {
                V3 v{ static_cast<double>(io.rawInteger(lo, hi)),
                    static_cast<double>(io.rawInteger(lo, hi)),
                    static_cast<double>(io.rawInteger(lo, hi)) };
                if (mode == 1 && k == 3) { v[2] = vertices[vertices.size() - 3][2]; }
                vertices.push_back(v);
            }
            size_t b = 4 * static_cast<size_t>(t);
            tetra.push_back({ b, b + 1, b + 2, b + 3 });
        }
    }

    // Records the vertex count, the vertices, the tetrahedron count and the
    // index quadruples.
    void RecordMesh(oracle::Ctx& io, std::vector<V3> const& vertices, std::vector<T4> const& tetra)
    {
        io.given(static_cast<double>(vertices.size()));
        for (auto const& v : vertices) { io.given(v[0]); io.given(v[1]); io.given(v[2]); }
        io.given(static_cast<double>(tetra.size()));
        for (auto const& t : tetra)
        {
            for (int32_t j = 0; j < 4; ++j) { io.given(static_cast<double>(t[j])); }
        }
    }

    // Recorded region and bound.
    void RecordRegion(oracle::Ctx& io, V3 const& rmin, V3 const& rmax,
        std::array<size_t, 3> const& bound)
    {
        for (int32_t i = 0; i < 3; ++i) { io.given(rmin[i]); }
        for (int32_t i = 0; i < 3; ++i) { io.given(rmax[i]); }
        for (int32_t i = 0; i < 3; ++i) { io.given(static_cast<double>(bound[i])); }
    }

    __declspec(noinline) void Rasterize(TetrahedraRasterizer<double>& r, size_t numThreads,
        V3 const& rmin, V3 const& rmax, std::array<size_t, 3> const& bound,
        std::vector<int32_t>& grid)
    {
        r(numThreads, rmin, rmax, bound, grid);
    }
}

// Lattice vertices and a region [r0, r0 + b - 1] per axis with bound b, so
// the grid multiplier is exactly 1 and every predicate is evaluated exactly:
// voxel centres on faces, edges and vertices are exact ties. Vertices reach
// one past the region (clipping). numThreads is 0 or 1 (single-threaded;
// one worker thread). The replay checks every record against an exact
// BigInt point-in-tetrahedron reference.
ORACLE_CASE("TetrahedraRasterizer.lattice")
{
    std::array<size_t, 3> bound;
    V3 rmin, rmax;
    for (int32_t i = 0; i < 3; ++i)
    {
        bound[i] = io.rawInteger(2, 9);
        rmin[i] = io.rawInteger(-2, 2);
        rmax[i] = rmin[i] + static_cast<double>(bound[i] - 1);
    }
    int32_t lo = 100, hi = -100;
    for (int32_t i = 0; i < 3; ++i)
    {
        lo = std::min(lo, static_cast<int32_t>(rmin[i]) - 1);
        hi = std::max(hi, static_cast<int32_t>(rmin[i] + bound[i]));
    }
    std::vector<V3> vertices;
    std::vector<T4> tetra;
    DrawLatticeMesh(io, lo, hi, vertices, tetra);
    RecordMesh(io, vertices, tetra);
    RecordRegion(io, rmin, rmax, bound);
    size_t numThreads = io.integer(0, 1);
    TetrahedraRasterizer<double> r(vertices.size(), vertices.data(), tetra.size(), tetra.data());
    std::vector<int32_t> grid;
    Rasterize(r, numThreads, rmin, rmax, bound, grid);
    OutGrid(io, grid);
}

// Uniform vertices in [-1, 2]^3 (random tetrahedra of both orientations),
// a random region inside or straddling them (sometimes with dyadic width,
// so grid vertices fall on voxel centres), bounds 2..10.
ORACLE_CASE("TetrahedraRasterizer.general")
{
    std::vector<V3> vertices;
    std::vector<T4> tetra;
    int32_t numTetra = io.rawInteger(1, 6);
    for (int32_t t = 0; t < numTetra; ++t)
    {
        for (int32_t k = 0; k < 4; ++k)
        {
            vertices.push_back({ io.raw(-1.0, 2.0), io.raw(-1.0, 2.0), io.raw(-1.0, 2.0) });
        }
        size_t b = 4 * static_cast<size_t>(t);
        tetra.push_back({ b, b + 1, b + 2, b + 3 });
    }
    std::array<size_t, 3> bound;
    V3 rmin, rmax;
    bool dyadic = (io.rawInteger(0, 2) == 0);
    for (int32_t i = 0; i < 3; ++i)
    {
        bound[i] = io.rawInteger(2, 10);
        rmin[i] = (dyadic ? -0.25 * io.rawInteger(0, 4) : io.raw(-1.0, 0.5));
        rmax[i] = rmin[i] + (dyadic ? 0.25 * static_cast<double>(bound[i] - 1) : io.raw(0.5, 2.5));
    }
    RecordMesh(io, vertices, tetra);
    RecordRegion(io, rmin, rmax, bound);
    size_t numThreads = io.integer(0, 1);
    TetrahedraRasterizer<double> r(vertices.size(), vertices.data(), tetra.size(), tetra.data());
    std::vector<int32_t> grid;
    Rasterize(r, numThreads, rmin, rmax, bound, grid);
    OutGrid(io, grid);
}

// Two operator() calls on one rasterizer: ClipCullAABBs clips the stored
// boxes in place (#439, preserved), so the second call (a larger or a
// shifted region, recorded) scans the boxes the first call left. Lattice
// meshes; the multipliers are dyadic.
ORACLE_CASE("TetrahedraRasterizer.twice")
{
    std::array<size_t, 3> bound0, bound1;
    V3 rmin0, rmax0, rmin1, rmax1;
    for (int32_t i = 0; i < 3; ++i)
    {
        bound0[i] = io.rawInteger(2, 7);
        rmin0[i] = io.rawInteger(-1, 2);
        rmax0[i] = rmin0[i] + static_cast<double>(bound0[i] - 1);
        bound1[i] = io.rawInteger(2, 9);
        rmin1[i] = rmin0[i] - io.rawInteger(0, 2);
        rmax1[i] = rmin1[i] + 0.5 * static_cast<double>(bound1[i] - 1) * io.rawInteger(1, 2);
    }
    std::vector<V3> vertices;
    std::vector<T4> tetra;
    DrawLatticeMesh(io, -2, 6, vertices, tetra);
    RecordMesh(io, vertices, tetra);
    RecordRegion(io, rmin0, rmax0, bound0);
    RecordRegion(io, rmin1, rmax1, bound1);
    TetrahedraRasterizer<double> r(vertices.size(), vertices.data(), tetra.size(), tetra.data());
    std::vector<int32_t> grid;
    Rasterize(r, 0, rmin0, rmax0, bound0, grid);
    OutGrid(io, grid);
    Rasterize(r, 0, rmin1, rmax1, bound1, grid);
    OutGrid(io, grid);
}

// The multithreaded rasterizer (2..4 std::threads; the port has one code
// path). Upstream's threads write the grid concurrently, so overlapping
// tetrahedra would make the result depend on the schedule; here every
// tetrahedron lies in its own cell of the 2x2x2 split of the region, the
// cells' voxel ranges are disjoint, and the result is deterministic.
ORACLE_CASE("TetrahedraRasterizer.threads")
{
    std::array<size_t, 3> bound;
    V3 rmin, rmax;
    int32_t half[3];
    for (int32_t i = 0; i < 3; ++i)
    {
        bound[i] = io.rawInteger(4, 9);
        rmin[i] = 0.0;
        rmax[i] = static_cast<double>(bound[i] - 1);
        half[i] = static_cast<int32_t>(bound[i]) / 2;
    }
    int32_t cells[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
    for (int32_t k = 7; k > 0; --k) { std::swap(cells[k], cells[io.rawInteger(0, k)]); }
    int32_t numTetra = io.rawInteger(2, 8);
    std::vector<V3> vertices;
    std::vector<T4> tetra;
    for (int32_t t = 0; t < numTetra; ++t)
    {
        int32_t lo[3], hi[3];
        for (int32_t i = 0; i < 3; ++i)
        {
            bool upper = ((cells[t] >> i) & 1) != 0;
            lo[i] = (upper ? half[i] : 0);
            hi[i] = (upper ? static_cast<int32_t>(bound[i]) - 1 : half[i] - 1);
        }
        for (int32_t k = 0; k < 4; ++k)
        {
            vertices.push_back({ static_cast<double>(io.rawInteger(lo[0], hi[0])),
                static_cast<double>(io.rawInteger(lo[1], hi[1])),
                static_cast<double>(io.rawInteger(lo[2], hi[2])) });
        }
        size_t b = 4 * static_cast<size_t>(t);
        tetra.push_back({ b, b + 1, b + 2, b + 3 });
    }
    RecordMesh(io, vertices, tetra);
    RecordRegion(io, rmin, rmax, bound);
    size_t numThreads = io.integer(2, 4);
    TetrahedraRasterizer<double> r(vertices.size(), vertices.data(), tetra.size(), tetra.data());
    std::vector<int32_t> grid;
    Rasterize(r, numThreads, rmin, rmax, bound, grid);
    OutGrid(io, grid);
}

// Throw parity: no vertices or no tetrahedra (constructor LogError), a
// bound component below 2 (operator() LogError), or a valid call.
ORACLE_CASE("TetrahedraRasterizer.invalid")
{
    int32_t which = io.integer(0, 3);
    std::vector<V3> vertices{ { 0.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0 }, { 0.0, 2.0, 0.0 },
        { 0.0, 0.0, 2.0 } };
    std::vector<T4> tetra{ { 0, 1, 2, 3 } };
    if (which == 0) { vertices.clear(); }
    if (which == 1) { tetra.clear(); }
    std::array<size_t, 3> bound{ 3, 3, 3 };
    if (which == 2)
    {
        int32_t axis = io.integer(0, 2);
        int32_t value = io.integer(0, 1);
        bound[axis] = value;
    }
    TetrahedraRasterizer<double> r(vertices.size(), vertices.data(), tetra.size(), tetra.data());
    std::vector<int32_t> grid;
    Rasterize(r, 0, V3{ 0.0, 0.0, 0.0 }, V3{ 2.0, 2.0, 2.0 }, bound, grid);
    OutGrid(io, grid);
}
