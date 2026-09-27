// Group 28 (interpolation): Hermite cubic/quintic polynomials on the unit
// cell and their tensor products, the Akima interpolators (the IntpAkima1
// base through IntpAkimaUniform1, IntpAkimaUniform2, IntpAkimaUniform3) and
// the uniform bilinear and bicubic interpolators.
//
// Everything here is + - * / and a truncating float-to-int conversion for
// the cell index, so every case is declared exact on the TypeScript side.
//
// Generators. Sample values are drawn uniform, on a small integer lattice,
// dyadic, or sparse (zeros of both signs); lattice data produce the equal
// slopes that select the non-generic branches of Akima's ComputeDerivative.
// Query points are drawn inside cells, exactly on nodes, on the upper domain
// boundary, outside the domain on both sides (the clamping and the
// documented-but-absent clamping, #69), and as -0 against a zero minimum
// (std::max(-0, +0) keeps -0, Math.max would not).
//
// IntpAkimaUniform2/3 (#58): upstream's GetFXY/GetFXZ/GetFYZ/GetFXYZ reuse
// the min-boundary one-sided stencils at the max boundaries without the sign
// flip, and the port fixes it. The main cases therefore query only cells
// whose corners avoid every max face (probed with upstream's own clamp and
// lookup expressions), or use additively separable data, on which every
// mixed-derivative stencil is exactly zero; the '.maxBoundary' cases query the
// cells on the max faces and are declared deviations.
#define ORACLE_FAMILY "v28-interpolation"
#include "Oracle.h"

#include <Mathematics/HermiteCubic.h>
#include <Mathematics/HermiteQuintic.h>
#include <Mathematics/HermiteBicubic.h>
#include <Mathematics/HermiteBiquintic.h>
#include <Mathematics/HermiteTricubic.h>
#include <Mathematics/HermiteTriquintic.h>
#include <Mathematics/IntpAkimaUniform1.h>
#include <Mathematics/IntpAkimaUniform2.h>
#include <Mathematics/IntpAkimaUniform3.h>
#include <Mathematics/IntpBicubic2.h>
#include <Mathematics/IntpBilinear2.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

using namespace gte;

namespace
{
    // A sample value, drawn raw (the caller records it).
    //   mode 0: uniform in [-10, 10]
    //   mode 1: integer lattice [-3, 3]
    //   mode 2: dyadic k/8 in [-4, 4]
    //   mode 3: sparse: +0, -0 or a small integer
    double RawValue(oracle::Ctx& io, int mode)
    {
        switch (mode)
        {
        case 0:
            return io.raw(-10.0, 10.0);
        case 1:
            return static_cast<double>(io.rawInteger(-3, 3));
        case 2:
            return io.rawInteger(-32, 32) / 8.0;
        default:
        {
            int k = io.rawInteger(0, 4);
            if (k == 0) { return 0.0; }
            if (k == 1) { return -0.0; }
            if (k == 2) { return static_cast<double>(io.rawInteger(-2, 2)); }
            return 0.0;
        }
        }
    }

    double Value(oracle::Ctx& io, int mode)
    {
        return io.given(RawValue(io, mode));
    }

    // A local coordinate for the unit-cell Hermite polynomials.
    //   0: uniform in [0, 1]   1: dyadic k/16 in [0, 1]
    //   2: an endpoint +0, -0 or 1   3: outside the cell, in [-0.5, 1.5]
    double UnitParam(oracle::Ctx& io)
    {
        int m = io.rawInteger(0, 3);
        double t;
        if (m == 0) { t = io.raw(0.0, 1.0); }
        else if (m == 1) { t = io.rawInteger(0, 16) / 16.0; }
        else if (m == 2)
        {
            int k = io.rawInteger(0, 2);
            t = (k == 0 ? 0.0 : (k == 1 ? -0.0 : 1.0));
        }
        else { t = io.raw(-0.5, 1.5); }
        return io.given(t);
    }

    // ---------------------------------------------------------- Hermite 1D

    std::array<HermiteCubic<double>::Sample, 2> CubicBlocks(oracle::Ctx& io)
    {
        int mode = io.index() % 4;
        std::array<HermiteCubic<double>::Sample, 2> blocks;
        for (size_t b = 0; b < 2; ++b)
        {
            blocks[b].F = Value(io, mode);
            blocks[b].Fx = Value(io, mode);
        }
        return blocks;
    }

    std::array<HermiteQuintic<double>::Sample, 2> QuinticBlocks(oracle::Ctx& io)
    {
        int mode = io.index() % 4;
        std::array<HermiteQuintic<double>::Sample, 2> blocks;
        for (size_t b = 0; b < 2; ++b)
        {
            blocks[b].F = Value(io, mode);
            blocks[b].Fx = Value(io, mode);
            blocks[b].Fxx = Value(io, mode);
        }
        return blocks;
    }
}

ORACLE_CASE("HermiteCubic.evaluate")
{
    auto blocks = CubicBlocks(io);
    HermiteCubic<double> hermite(blocks);
    for (int q = 0; q < 3; ++q)
    {
        double x = UnitParam(io);
        for (size_t order = 0; order <= 4; ++order)
        {
            io.outReal(hermite(order, x));
        }
    }
}

ORACLE_CASE("HermiteCubic.generate")
{
    // Generate on a default-constructed (zero) polynomial, then on a second
    // sample set, overwriting every coefficient.
    HermiteCubic<double> hermite;
    double x = UnitParam(io);
    io.outReal(hermite(0, x));
    auto blocks0 = CubicBlocks(io);
    hermite.Generate(blocks0);
    for (size_t order = 0; order <= 3; ++order) { io.outReal(hermite(order, x)); }
    auto blocks1 = CubicBlocks(io);
    hermite.Generate(blocks1);
    for (size_t order = 0; order <= 3; ++order) { io.outReal(hermite(order, x)); }
}

ORACLE_CASE("HermiteCubic.P")
{
    double t = UnitParam(io);
    for (size_t select = 0; select < 4; ++select)
    {
        for (size_t order = 0; order <= 4; ++order)
        {
            io.outReal(HermiteCubic<double>::P(select, order, t));
        }
    }
}

ORACLE_CASE("HermiteQuintic.evaluate")
{
    auto blocks = QuinticBlocks(io);
    HermiteQuintic<double> hermite(blocks);
    for (int q = 0; q < 3; ++q)
    {
        double x = UnitParam(io);
        for (size_t order = 0; order <= 6; ++order)
        {
            io.outReal(hermite(order, x));
        }
    }
}

ORACLE_CASE("HermiteQuintic.generate")
{
    HermiteQuintic<double> hermite;
    double x = UnitParam(io);
    io.outReal(hermite(0, x));
    auto blocks0 = QuinticBlocks(io);
    hermite.Generate(blocks0);
    for (size_t order = 0; order <= 5; ++order) { io.outReal(hermite(order, x)); }
    auto blocks1 = QuinticBlocks(io);
    hermite.Generate(blocks1);
    for (size_t order = 0; order <= 5; ++order) { io.outReal(hermite(order, x)); }
}

ORACLE_CASE("HermiteQuintic.P")
{
    double t = UnitParam(io);
    for (size_t select = 0; select < 6; ++select)
    {
        for (size_t order = 0; order <= 6; ++order)
        {
            io.outReal(HermiteQuintic<double>::P(select, order, t));
        }
    }
}

// ---------------------------------------------------------- Hermite 2D/3D
//
// Each tensor-product case records the samples (blocks[b0][b1][b2] in index
// order, fields in declaration order), emits every generated coefficient
// c[i][j][k] in index order, then evaluates. The '.manual' cases set the
// public coefficient array directly, as the header comment invites; the
// '.default' cases evaluate the identically zero polynomial.

namespace
{
    std::array<std::array<HermiteBicubic<double>::Sample, 2>, 2> BicubicBlocks(oracle::Ctx& io)
    {
        int mode = io.index() % 4;
        std::array<std::array<HermiteBicubic<double>::Sample, 2>, 2> blocks;
        for (size_t b0 = 0; b0 < 2; ++b0)
        {
            for (size_t b1 = 0; b1 < 2; ++b1)
            {
                auto& s = blocks[b0][b1];
                s.F = Value(io, mode);
                s.Fx = Value(io, mode);
                s.Fy = Value(io, mode);
                s.Fxy = Value(io, mode);
            }
        }
        return blocks;
    }

    std::array<std::array<HermiteBiquintic<double>::Sample, 2>, 2> BiquinticBlocks(oracle::Ctx& io)
    {
        int mode = io.index() % 4;
        std::array<std::array<HermiteBiquintic<double>::Sample, 2>, 2> blocks;
        for (size_t b0 = 0; b0 < 2; ++b0)
        {
            for (size_t b1 = 0; b1 < 2; ++b1)
            {
                auto& s = blocks[b0][b1];
                s.F = Value(io, mode);
                s.Fx = Value(io, mode);
                s.Fy = Value(io, mode);
                s.Fxx = Value(io, mode);
                s.Fxy = Value(io, mode);
                s.Fyy = Value(io, mode);
                s.Fxxy = Value(io, mode);
                s.Fxyy = Value(io, mode);
                s.Fxxyy = Value(io, mode);
            }
        }
        return blocks;
    }

    using TricubicBlocksT = std::array<std::array<std::array<HermiteTricubic<double>::Sample, 2>, 2>, 2>;
    using TriquinticBlocksT = std::array<std::array<std::array<HermiteTriquintic<double>::Sample, 2>, 2>, 2>;

    TricubicBlocksT TricubicBlocks(oracle::Ctx& io)
    {
        int mode = io.index() % 4;
        TricubicBlocksT blocks;
        for (size_t b0 = 0; b0 < 2; ++b0)
        {
            for (size_t b1 = 0; b1 < 2; ++b1)
            {
                for (size_t b2 = 0; b2 < 2; ++b2)
                {
                    auto& s = blocks[b0][b1][b2];
                    s.F = Value(io, mode);
                    s.Fx = Value(io, mode);
                    s.Fy = Value(io, mode);
                    s.Fz = Value(io, mode);
                    s.Fxy = Value(io, mode);
                    s.Fxz = Value(io, mode);
                    s.Fyz = Value(io, mode);
                    s.Fxyz = Value(io, mode);
                }
            }
        }
        return blocks;
    }

    // The 27 fields of HermiteTriquintic::Sample in declaration order.
    std::array<double HermiteTriquintic<double>::Sample::*, 27> const kTriquinticFields =
    { {
        &HermiteTriquintic<double>::Sample::F,
        &HermiteTriquintic<double>::Sample::Fx,
        &HermiteTriquintic<double>::Sample::Fy,
        &HermiteTriquintic<double>::Sample::Fz,
        &HermiteTriquintic<double>::Sample::Fxx,
        &HermiteTriquintic<double>::Sample::Fxy,
        &HermiteTriquintic<double>::Sample::Fxz,
        &HermiteTriquintic<double>::Sample::Fyy,
        &HermiteTriquintic<double>::Sample::Fyz,
        &HermiteTriquintic<double>::Sample::Fzz,
        &HermiteTriquintic<double>::Sample::Fxxy,
        &HermiteTriquintic<double>::Sample::Fxxz,
        &HermiteTriquintic<double>::Sample::Fxyy,
        &HermiteTriquintic<double>::Sample::Fxyz,
        &HermiteTriquintic<double>::Sample::Fxzz,
        &HermiteTriquintic<double>::Sample::Fyyz,
        &HermiteTriquintic<double>::Sample::Fyzz,
        &HermiteTriquintic<double>::Sample::Fxxyy,
        &HermiteTriquintic<double>::Sample::Fxxyz,
        &HermiteTriquintic<double>::Sample::Fxxzz,
        &HermiteTriquintic<double>::Sample::Fxyyz,
        &HermiteTriquintic<double>::Sample::Fxyzz,
        &HermiteTriquintic<double>::Sample::Fyyzz,
        &HermiteTriquintic<double>::Sample::Fxxyyz,
        &HermiteTriquintic<double>::Sample::Fxxyzz,
        &HermiteTriquintic<double>::Sample::Fxyyzz,
        &HermiteTriquintic<double>::Sample::Fxxyyzz
    } };

    TriquinticBlocksT TriquinticBlocks(oracle::Ctx& io)
    {
        int mode = io.index() % 4;
        TriquinticBlocksT blocks;
        for (size_t b0 = 0; b0 < 2; ++b0)
        {
            for (size_t b1 = 0; b1 < 2; ++b1)
            {
                for (size_t b2 = 0; b2 < 2; ++b2)
                {
                    for (auto field : kTriquinticFields)
                    {
                        blocks[b0][b1][b2].*field = Value(io, mode);
                    }
                }
            }
        }
        return blocks;
    }

    template <size_t N>
    void OutCoeff2(oracle::Ctx& io, std::array<std::array<double, N>, N> const& c)
    {
        for (size_t i = 0; i < N; ++i) { for (size_t j = 0; j < N; ++j) { io.outReal(c[i][j]); } }
    }

    template <size_t N>
    void OutCoeff3(oracle::Ctx& io, std::array<std::array<std::array<double, N>, N>, N> const& c)
    {
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                for (size_t k = 0; k < N; ++k) { io.outReal(c[i][j][k]); }
            }
        }
    }

    // All (xOrder, yOrder) in [0, maxOrder]^2, row-major in xOrder.
    template <typename H>
    void OutOrders2(oracle::Ctx& io, H const& hermite, size_t maxOrder, double x, double y)
    {
        for (size_t xo = 0; xo <= maxOrder; ++xo)
        {
            for (size_t yo = 0; yo <= maxOrder; ++yo) { io.outReal(hermite(xo, yo, x, y)); }
        }
    }

    // A drawn order triple, each order in [0, maxOrder].
    template <typename H>
    void OutDrawnOrders3(oracle::Ctx& io, H const& hermite, int maxOrder, int count,
        double x, double y, double z)
    {
        for (int i = 0; i < count; ++i)
        {
            size_t xo = static_cast<size_t>(io.integer(0, maxOrder));
            size_t yo = static_cast<size_t>(io.integer(0, maxOrder));
            size_t zo = static_cast<size_t>(io.integer(0, maxOrder));
            io.outReal(hermite(xo, yo, zo, x, y, z));
        }
    }
}

ORACLE_CASE("HermiteBicubic.evaluate")
{
    auto blocks = BicubicBlocks(io);
    HermiteBicubic<double> hermite(blocks);
    OutCoeff2(io, hermite.c);
    for (int q = 0; q < 2; ++q)
    {
        double x = UnitParam(io);
        double y = UnitParam(io);
        OutOrders2(io, hermite, 4, x, y);
    }
}

ORACLE_CASE("HermiteBicubic.manual")
{
    HermiteBicubic<double> hermite;
    double x0 = UnitParam(io);
    double y0 = UnitParam(io);
    io.outReal(hermite(0, 0, x0, y0));
    int mode = io.index() % 4;
    for (size_t i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < 4; ++j) { hermite.c[i][j] = Value(io, mode); }
    }
    double x = UnitParam(io);
    double y = UnitParam(io);
    OutOrders2(io, hermite, 4, x, y);
    // Generate overwrites every manual coefficient.
    auto blocks = BicubicBlocks(io);
    hermite.Generate(blocks);
    OutCoeff2(io, hermite.c);
}

ORACLE_CASE("HermiteBiquintic.evaluate")
{
    auto blocks = BiquinticBlocks(io);
    HermiteBiquintic<double> hermite(blocks);
    OutCoeff2(io, hermite.c);
    for (int q = 0; q < 2; ++q)
    {
        double x = UnitParam(io);
        double y = UnitParam(io);
        OutOrders2(io, hermite, 6, x, y);
    }
}

ORACLE_CASE("HermiteBiquintic.manual")
{
    HermiteBiquintic<double> hermite;
    double x0 = UnitParam(io);
    double y0 = UnitParam(io);
    io.outReal(hermite(0, 0, x0, y0));
    int mode = io.index() % 4;
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 6; ++j) { hermite.c[i][j] = Value(io, mode); }
    }
    double x = UnitParam(io);
    double y = UnitParam(io);
    OutOrders2(io, hermite, 6, x, y);
    auto blocks = BiquinticBlocks(io);
    hermite.Generate(blocks);
    OutCoeff2(io, hermite.c);
}

ORACLE_CASE("HermiteTricubic.evaluate")
{
    auto blocks = TricubicBlocks(io);
    HermiteTricubic<double> hermite(blocks);
    OutCoeff3(io, hermite.c);
    double x = UnitParam(io);
    double y = UnitParam(io);
    double z = UnitParam(io);
    for (size_t xo = 0; xo <= 3; ++xo)
    {
        for (size_t yo = 0; yo <= 3; ++yo)
        {
            for (size_t zo = 0; zo <= 3; ++zo) { io.outReal(hermite(xo, yo, zo, x, y, z)); }
        }
    }
    // Orders up to 4 (4 in any axis is the zero branch).
    OutDrawnOrders3(io, hermite, 4, 4, x, y, z);
}

ORACLE_CASE("HermiteTricubic.manual")
{
    HermiteTricubic<double> hermite;
    double x0 = UnitParam(io);
    double y0 = UnitParam(io);
    double z0 = UnitParam(io);
    io.outReal(hermite(0, 0, 0, x0, y0, z0));
    int mode = io.index() % 4;
    for (size_t i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            for (size_t k = 0; k < 4; ++k) { hermite.c[i][j][k] = Value(io, mode); }
        }
    }
    double x = UnitParam(io);
    double y = UnitParam(io);
    double z = UnitParam(io);
    OutDrawnOrders3(io, hermite, 4, 16, x, y, z);
}

ORACLE_CASE("HermiteTriquintic.evaluate")
{
    auto blocks = TriquinticBlocks(io);
    HermiteTriquintic<double> hermite(blocks);
    OutCoeff3(io, hermite.c);
    double x = UnitParam(io);
    double y = UnitParam(io);
    double z = UnitParam(io);
    io.outReal(hermite(0, 0, 0, x, y, z));
    OutDrawnOrders3(io, hermite, 6, 20, x, y, z);
}

ORACLE_CASE("HermiteTriquintic.manual")
{
    HermiteTriquintic<double> hermite;
    double x0 = UnitParam(io);
    double y0 = UnitParam(io);
    double z0 = UnitParam(io);
    io.outReal(hermite(0, 0, 0, x0, y0, z0));
    // Integer coefficients; HermiteTriquintic.evaluate covers real ones.
    for (size_t i = 0; i < 6; ++i)
    {
        for (size_t j = 0; j < 6; ++j)
        {
            for (size_t k = 0; k < 6; ++k) { hermite.c[i][j][k] = io.lattice(-3, 3); }
        }
    }
    double x = UnitParam(io);
    double y = UnitParam(io);
    double z = UnitParam(io);
    io.outReal(hermite(0, 0, 0, x, y, z));
    OutDrawnOrders3(io, hermite, 6, 12, x, y, z);
}

// ---------------------------------------------------------- grids

namespace
{
    double RawMin(oracle::Ctx& io)
    {
        int m = io.rawInteger(0, 2);
        if (m == 0) { return io.rawInteger(-16, 16) / 4.0; }
        if (m == 1) { return io.raw(-5.0, 5.0); }
        return 0.0;
    }

    double RawSpacing(oracle::Ctx& io)
    {
        static double const dyadic[7] = { 0.25, 0.5, 1.0, 2.0, 0.75, 1.5, 3.0 };
        int m = io.rawInteger(0, 3);
        if (m == 0) { return dyadic[io.rawInteger(0, 6)]; }
        if (m == 1) { return io.raw(0.1, 3.0); }
        if (m == 2) { return 1.0; }
        return io.rawInteger(1, 8) / 4.0;
    }

    // Validation draws: on even records every parameter is valid; on odd
    // records each parameter is invalid with probability 1/2.
    bool DrawInvalid(oracle::Ctx& io)
    {
        return io.index() % 2 == 1 && io.rawInteger(0, 1) == 0;
    }

    int ValidationBound(oracle::Ctx& io, int lowest, int minValid, int maxValid)
    {
        bool invalid = DrawInvalid(io);
        return (invalid ? io.integer(lowest, minValid - 1) : io.integer(minValid, maxValid));
    }

    double ValidationSpacing(oracle::Ctx& io)
    {
        bool invalid = DrawInvalid(io);
        int k = io.rawInteger(0, 1);
        return io.given(invalid ? (k == 0 ? -1.0 : 0.0) : (k == 0 ? 0.5 : 1.25));
    }

    // Upstream's clamp and cell lookup (IntpAkima1/IntpAkimaUniform2/3):
    // the cell whose lower node is the last node below x, the last cell for
    // x at or above the last interior node.
    int AkimaCell(double x, double min, double spacing, int bound)
    {
        double max = min + spacing * (static_cast<double>(bound) - 1.0);
        x = std::min(std::max(x, min), max);
        for (int i = 0, ip1 = 1; ip1 < bound; ++i, ++ip1)
        {
            if (x < min + spacing * static_cast<double>(ip1)) { return i; }
        }
        return bound - 2;
    }

    // A raw query coordinate on a uniform grid with nodes min + k*spacing,
    // k = 0..bound-1. Inside modes pick a cell in [cellLo, cellHi].
    //   0 uniform inside the cell   1 the cell's lower node
    //   2 dyadic inside the cell    3 below min, up to two spacings
    //   4 above max, up to two spacings   5 exactly max
    //   6 -0 when min is zero, else min
    // Modes 3-6 are replaced by mode 0 when 'outside' is false.
    double RawCoord(oracle::Ctx& io, double min, double spacing, int bound,
        int cellLo, int cellHi, bool outside)
    {
        int m = io.rawInteger(0, outside ? 6 : 2);
        double max = min + spacing * (static_cast<double>(bound) - 1.0);
        int k = io.rawInteger(cellLo, cellHi);
        switch (m)
        {
        case 1: return min + spacing * static_cast<double>(k);
        case 2: return min + spacing * (k + io.rawInteger(0, 7) / 8.0);
        case 3: return min - spacing * io.raw(0.0, 2.0);
        case 4: return max + spacing * io.raw(0.0, 2.0);
        case 5: return max;
        case 6: return (min == 0.0 ? -0.0 : min);
        default: return min + spacing * (k + io.raw(0.0, 1.0));
        }
    }

    // A raw query coordinate whose Akima cell is at most bound-3, so that no
    // corner of its cell lies on the max face (the #58-sound region). Draws
    // below min are allowed: they clamp to min, in cell 0. The acceptance
    // test is upstream's own clamp and lookup; capped, falling back to min.
    double RawSoundCoord(oracle::Ctx& io, double min, double spacing, int bound)
    {
        for (int attempt = 0; attempt < 16; ++attempt)
        {
            double x = RawCoord(io, min, spacing, bound, 0, bound - 3, false);
            if (io.rawInteger(0, 5) == 0)
            {
                x = (min == 0.0 && io.rawInteger(0, 1) == 0 ? -0.0 : min - spacing * io.raw(0.0, 2.0));
            }
            if (AkimaCell(x, min, spacing, bound) <= bound - 3) { return x; }
        }
        return min;
    }

    // A raw query coordinate in the last cell of the Akima lookup (a corner
    // on the max face): inside it, its lower node, exactly max, or above max.
    double RawLastCellCoord(oracle::Ctx& io, double min, double spacing, int bound)
    {
        int m = io.rawInteger(0, 3);
        double max = min + spacing * (static_cast<double>(bound) - 1.0);
        double node = min + spacing * static_cast<double>(bound - 2);
        switch (m)
        {
        case 1: return node;
        case 2: return max;
        case 3: return max + spacing * io.raw(0.0, 2.0);
        default: return min + spacing * (bound - 2 + io.raw(0.0, 1.0));
        }
    }

    // 1D Akima data.
    //   0 uniform   1 lattice [-2, 2] (ties)   2 nondecreasing dyadic steps
    //   3 runs of equal slopes   4 cubic in the node index
    std::vector<double> RawAkimaData1(oracle::Ctx& io, int mode, int n)
    {
        std::vector<double> F(static_cast<size_t>(n));
        double run = 0.0, value = 0.0;
        int a = io.rawInteger(-2, 2), b = io.rawInteger(-2, 2);
        int c = io.rawInteger(-2, 2), d = io.rawInteger(-1, 1);
        for (int i = 0; i < n; ++i)
        {
            double v;
            switch (mode)
            {
            case 0: v = io.raw(-10.0, 10.0); break;
            case 1: v = static_cast<double>(io.rawInteger(-2, 2)); break;
            case 2: value += io.rawInteger(0, 4) / 4.0; v = value; break;
            case 3:
                if (i == 0 || io.rawInteger(0, 2) == 0) { run = io.rawInteger(-2, 2); }
                value += run; v = value; break;
            default: v = a + i * (b + i * (c + i * d)); break;
            }
            F[static_cast<size_t>(i)] = v;
        }
        return F;
    }
}

ORACLE_CASE("IntpAkima1.evaluate.uniform1")
{
    // The abstract IntpAkima1 (clamp, lookup dispatch, Polynomial,
    // ComputeDerivative) through its concrete subclass IntpAkimaUniform1,
    // whose own constructor and Lookup belong to group 29.
    int mode = io.index() % 5;
    int quantity = io.integer(3, 8);
    double xMin = io.given(RawMin(io));
    double xSpacing = io.given(RawSpacing(io));
    std::vector<double> F = RawAkimaData1(io, mode, quantity);
    for (auto& f : F) { f = io.given(f); }
    IntpAkimaUniform1<double> interp(quantity, xMin, xSpacing, F.data());
    io.outInt(interp.GetQuantity());
    io.outReal(interp.GetXMin());
    io.outReal(interp.GetXMax());
    for (int q = 0; q < 5; ++q)
    {
        double x = io.given(RawCoord(io, xMin, xSpacing, quantity, 0, quantity - 2, true));
        io.outReal(interp(x));
        for (int order = -1; order <= 4; ++order) { io.outReal(interp(order, x)); }
    }
}

ORACLE_CASE("IntpAkima1.construct.validation")
{
    // IntpAkima1 asserts quantity >= 3; IntpAkimaUniform1 asserts a positive
    // spacing after the base constructor.
    int quantity = ValidationBound(io, 1, 3, 4);
    double xSpacing = ValidationSpacing(io);
    std::vector<double> F(4);
    for (auto& f : F) { f = io.lattice(-3, 3); }
    IntpAkimaUniform1<double> interp(quantity, 0.0, xSpacing, F.data());
    io.outReal(interp(0.25));
    io.outReal(interp(1, 0.75));
}

// ---------------------------------------------------------- IntpAkimaUniform2

namespace
{
    // Data on an nx-by-ny lattice, row-major F[i + nx*j].
    //   0 uniform   1 lattice [-3, 3]   2 additively separable g(i) + h(j)
    //   3 bilinear a + b*i + c*j + e*i*j   4 plateau {0, 1}   5 cubic in (i, j)
    std::vector<double> RawGrid2(oracle::Ctx& io, int mode, int nx, int ny)
    {
        std::vector<double> g(static_cast<size_t>(nx)), h(static_cast<size_t>(ny));
        for (auto& v : g) { v = io.rawInteger(-3, 3); }
        for (auto& v : h) { v = io.rawInteger(-3, 3); }
        double a = io.rawInteger(-3, 3), b = io.rawInteger(-3, 3);
        double c = io.rawInteger(-3, 3), e = io.rawInteger(-2, 2);
        double p = io.rawInteger(-1, 1), q = io.rawInteger(-1, 1);
        std::vector<double> F(static_cast<size_t>(nx) * static_cast<size_t>(ny));
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                double v;
                switch (mode)
                {
                case 0: v = io.raw(-10.0, 10.0); break;
                case 1: v = io.rawInteger(-3, 3); break;
                case 2: v = g[i] + h[j]; break;
                case 3: v = a + b * i + c * j + e * i * j; break;
                case 4: v = io.rawInteger(0, 1); break;
                default: v = a + b * i + c * j + e * i * j + p * i * i * j + q * i * j * j * j; break;
                }
                F[static_cast<size_t>(i) + static_cast<size_t>(nx) * j] = v;
            }
        }
        return F;
    }

    struct Grid2
    {
        int nx = 0, ny = 0;
        double xMin = 0.0, xSpacing = 0.0, yMin = 0.0, ySpacing = 0.0;
        std::vector<double> F{};
    };

    // Records xBound, yBound, xMin, xSpacing, yMin, ySpacing, F.
    Grid2 RecordGrid2(oracle::Ctx& io, int mode, int minBound, int maxBound)
    {
        Grid2 grid{};
        grid.nx = io.integer(minBound, maxBound);
        grid.ny = io.integer(minBound, maxBound);
        grid.xMin = io.given(RawMin(io));
        grid.xSpacing = io.given(RawSpacing(io));
        grid.yMin = io.given(RawMin(io));
        grid.ySpacing = io.given(RawSpacing(io));
        grid.F = RawGrid2(io, mode, grid.nx, grid.ny);
        for (auto& f : grid.F) { f = io.given(f); }
        return grid;
    }

    void OutAkima2(oracle::Ctx& io, IntpAkimaUniform2<double> const& interp, double x, double y)
    {
        io.outReal(interp(x, y));
        for (int xo = 0; xo <= 3; ++xo)
        {
            for (int yo = 0; yo <= 3; ++yo) { io.outReal(interp(xo, yo, x, y)); }
        }
        io.outReal(interp(4, 0, x, y));
        io.outReal(interp(0, 4, x, y));
        io.outReal(interp(-1, 1, x, y));
    }

    void OutAkima2Header(oracle::Ctx& io, IntpAkimaUniform2<double> const& interp)
    {
        io.outInt(interp.GetQuantity());
        io.outReal(interp.GetXMax());
        io.outReal(interp.GetYMax());
    }
}

ORACLE_CASE("IntpAkimaUniform2.evaluate")
{
    // Queries restricted to cells none of whose corners lies on a max face,
    // where upstream's mixed-derivative estimates are sound (#58).
    static int const modes[5] = { 0, 1, 3, 4, 5 };
    int mode = modes[io.index() % 5];
    Grid2 grid = RecordGrid2(io, mode, 3, 6);
    IntpAkimaUniform2<double> interp(grid.nx, grid.ny, grid.xMin, grid.xSpacing,
        grid.yMin, grid.ySpacing, grid.F.data());
    OutAkima2Header(io, interp);
    for (int q = 0; q < 3; ++q)
    {
        double x = io.given(RawSoundCoord(io, grid.xMin, grid.xSpacing, grid.nx));
        double y = io.given(RawSoundCoord(io, grid.yMin, grid.ySpacing, grid.ny));
        OutAkima2(io, interp, x, y);
    }
}

ORACLE_CASE("IntpAkimaUniform2.evaluate.separable")
{
    // g(x) + h(y) on an integer lattice: every mixed-derivative stencil,
    // defective or not, is exactly zero, so all cells are comparable.
    Grid2 grid = RecordGrid2(io, 2, 3, 6);
    IntpAkimaUniform2<double> interp(grid.nx, grid.ny, grid.xMin, grid.xSpacing,
        grid.yMin, grid.ySpacing, grid.F.data());
    OutAkima2Header(io, interp);
    for (int q = 0; q < 3; ++q)
    {
        double x = io.given(RawCoord(io, grid.xMin, grid.xSpacing, grid.nx, 0, grid.nx - 2, true));
        double y = io.given(RawCoord(io, grid.yMin, grid.ySpacing, grid.ny, 0, grid.ny - 2, true));
        OutAkima2(io, interp, x, y);
    }
}

ORACLE_CASE("IntpAkimaUniform2.evaluate.maxBoundary")
{
    // Deviation (#58): cells with a corner on a max face, generic data.
    static int const modes[4] = { 0, 1, 3, 5 };
    int mode = modes[io.index() % 4];
    Grid2 grid = RecordGrid2(io, mode, 3, 6);
    IntpAkimaUniform2<double> interp(grid.nx, grid.ny, grid.xMin, grid.xSpacing,
        grid.yMin, grid.ySpacing, grid.F.data());
    OutAkima2Header(io, interp);
    for (int q = 0; q < 2; ++q)
    {
        bool lastX = (io.rawInteger(0, 1) == 0);
        double rx = (lastX ? RawLastCellCoord(io, grid.xMin, grid.xSpacing, grid.nx)
            : RawCoord(io, grid.xMin, grid.xSpacing, grid.nx, 0, grid.nx - 2, true));
        double ry = (lastX ? RawCoord(io, grid.yMin, grid.ySpacing, grid.ny, 0, grid.ny - 2, true)
            : RawLastCellCoord(io, grid.yMin, grid.ySpacing, grid.ny));
        double x = io.given(rx);
        double y = io.given(ry);
        OutAkima2(io, interp, x, y);
    }
}

ORACLE_CASE("IntpAkimaUniform2.construct.validation")
{
    int xBound = ValidationBound(io, 2, 3, 4);
    int yBound = ValidationBound(io, 2, 3, 4);
    double xSpacing = ValidationSpacing(io);
    double ySpacing = ValidationSpacing(io);
    std::vector<double> F(static_cast<size_t>(xBound) * static_cast<size_t>(yBound));
    for (auto& f : F) { f = io.lattice(-3, 3); }
    IntpAkimaUniform2<double> interp(xBound, yBound, 0.0, xSpacing, 0.0, ySpacing, F.data());
    // Cell (0, 0) for either valid spacing, away from the #58 faces.
    io.outReal(interp(0.25, 0.125));
}

// ---------------------------------------------------------- IntpAkimaUniform3

namespace
{
    // Data on an nx-by-ny-by-nz lattice, F[i + nx*(j + ny*k)].
    //   0 uniform   1 lattice [-3, 3]   2 additively separable
    //   3 trilinear with integer coefficients   4 plateau {0, 1}
    //   5 cubic in (i, j, k)
    std::vector<double> RawGrid3(oracle::Ctx& io, int mode, int nx, int ny, int nz)
    {
        std::vector<double> g(static_cast<size_t>(nx)), h(static_cast<size_t>(ny));
        std::vector<double> w(static_cast<size_t>(nz));
        for (auto& v : g) { v = io.rawInteger(-3, 3); }
        for (auto& v : h) { v = io.rawInteger(-3, 3); }
        for (auto& v : w) { v = io.rawInteger(-3, 3); }
        std::array<double, 8> a{};
        for (auto& v : a) { v = io.rawInteger(-2, 2); }
        double p = io.rawInteger(-1, 1), q = io.rawInteger(-1, 1);
        std::vector<double> F(static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz));
        for (int k = 0; k < nz; ++k)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int i = 0; i < nx; ++i)
                {
                    double tri = a[0] + a[1] * i + a[2] * j + a[3] * k + a[4] * i * j
                        + a[5] * i * k + a[6] * j * k + a[7] * i * j * k;
                    double v;
                    switch (mode)
                    {
                    case 0: v = io.raw(-10.0, 10.0); break;
                    case 1: v = io.rawInteger(-3, 3); break;
                    case 2: v = g[i] + h[j] + w[k]; break;
                    case 3: v = tri; break;
                    case 4: v = io.rawInteger(0, 1); break;
                    default: v = tri + p * i * i * k + q * j * k * k * i; break;
                    }
                    F[static_cast<size_t>(i) + static_cast<size_t>(nx) * (j + static_cast<size_t>(ny) * k)] = v;
                }
            }
        }
        return F;
    }

    struct Grid3
    {
        int nx = 0, ny = 0, nz = 0;
        double xMin = 0.0, xSpacing = 0.0, yMin = 0.0, ySpacing = 0.0, zMin = 0.0, zSpacing = 0.0;
        std::vector<double> F{};
    };

    // Records xBound, yBound, zBound, xMin, xSpacing, yMin, ySpacing, zMin,
    // zSpacing, F.
    // With 'pow2Spacing' the spacings are powers of two, so that FXYZ's
    // masks invDXDYDZ * ODer[i] * ODer[j] * ODer[k] and their products with
    // integer samples are exact (the 2D-style FXY, FXZ, FYZ stencils scale an
    // exact sum and need no such restriction).
    Grid3 RecordGrid3(oracle::Ctx& io, int mode, bool pow2Spacing = false)
    {
        static double const pow2[5] = { 0.25, 0.5, 1.0, 2.0, 4.0 };
        Grid3 grid{};
        grid.nx = io.integer(3, 4);
        grid.ny = io.integer(3, 4);
        grid.nz = io.integer(3, 4);
        grid.xMin = io.given(RawMin(io));
        grid.xSpacing = io.given(pow2Spacing ? pow2[io.rawInteger(0, 4)] : RawSpacing(io));
        grid.yMin = io.given(RawMin(io));
        grid.ySpacing = io.given(pow2Spacing ? pow2[io.rawInteger(0, 4)] : RawSpacing(io));
        grid.zMin = io.given(RawMin(io));
        grid.zSpacing = io.given(pow2Spacing ? pow2[io.rawInteger(0, 4)] : RawSpacing(io));
        grid.F = RawGrid3(io, mode, grid.nx, grid.ny, grid.nz);
        for (auto& f : grid.F) { f = io.given(f); }
        return grid;
    }

    IntpAkimaUniform3<double> MakeAkima3(Grid3 const& grid)
    {
        return IntpAkimaUniform3<double>(grid.nx, grid.ny, grid.nz, grid.xMin, grid.xSpacing,
            grid.yMin, grid.ySpacing, grid.zMin, grid.zSpacing, grid.F.data());
    }

    void OutAkima3(oracle::Ctx& io, IntpAkimaUniform3<double> const& interp,
        double x, double y, double z)
    {
        io.outReal(interp(x, y, z));
        for (int xo = 0; xo <= 3; ++xo)
        {
            for (int yo = 0; yo <= 3; ++yo)
            {
                for (int zo = 0; zo <= 3; ++zo) { io.outReal(interp(xo, yo, zo, x, y, z)); }
            }
        }
        io.outReal(interp(4, 0, 0, x, y, z));
        io.outReal(interp(0, 0, 4, x, y, z));
        io.outReal(interp(1, -1, 0, x, y, z));
    }

    void OutAkima3Header(oracle::Ctx& io, IntpAkimaUniform3<double> const& interp)
    {
        io.outInt(interp.GetQuantity());
        io.outReal(interp.GetXMax());
        io.outReal(interp.GetYMax());
        io.outReal(interp.GetZMax());
    }
}

ORACLE_CASE("IntpAkimaUniform3.evaluate")
{
    // Queries restricted to cells none of whose corners lies on a max face,
    // where upstream's mixed-derivative estimates are sound (#58).
    static int const modes[5] = { 0, 1, 3, 4, 5 };
    int mode = modes[io.index() % 5];
    Grid3 grid = RecordGrid3(io, mode);
    auto interp = MakeAkima3(grid);
    OutAkima3Header(io, interp);
    for (int q = 0; q < 2; ++q)
    {
        double x = io.given(RawSoundCoord(io, grid.xMin, grid.xSpacing, grid.nx));
        double y = io.given(RawSoundCoord(io, grid.yMin, grid.ySpacing, grid.ny));
        double z = io.given(RawSoundCoord(io, grid.zMin, grid.zSpacing, grid.nz));
        OutAkima3(io, interp, x, y, z);
    }
}

ORACLE_CASE("IntpAkimaUniform3.evaluate.separable")
{
    // g(x) + h(y) + w(z) on an integer lattice with power-of-two spacings:
    // every mixed-derivative stencil is exactly zero, so all cells are
    // comparable. (With other spacings FXYZ keeps a rounding residue of the
    // exact zero, and the port negates that residue at the max faces.)
    Grid3 grid = RecordGrid3(io, 2, true);
    auto interp = MakeAkima3(grid);
    OutAkima3Header(io, interp);
    for (int q = 0; q < 2; ++q)
    {
        double x = io.given(RawCoord(io, grid.xMin, grid.xSpacing, grid.nx, 0, grid.nx - 2, true));
        double y = io.given(RawCoord(io, grid.yMin, grid.ySpacing, grid.ny, 0, grid.ny - 2, true));
        double z = io.given(RawCoord(io, grid.zMin, grid.zSpacing, grid.nz, 0, grid.nz - 2, true));
        OutAkima3(io, interp, x, y, z);
    }
}

ORACLE_CASE("IntpAkimaUniform3.evaluate.maxBoundary")
{
    // Deviation (#58): cells with a corner on at least one max face.
    static int const modes[4] = { 0, 1, 3, 5 };
    int mode = modes[io.index() % 4];
    Grid3 grid = RecordGrid3(io, mode);
    auto interp = MakeAkima3(grid);
    OutAkima3Header(io, interp);
    for (int q = 0; q < 2; ++q)
    {
        // A nonempty subset of the axes is placed in its last cell.
        int subset = io.rawInteger(1, 7);
        double rx = ((subset & 1) ? RawLastCellCoord(io, grid.xMin, grid.xSpacing, grid.nx)
            : RawCoord(io, grid.xMin, grid.xSpacing, grid.nx, 0, grid.nx - 2, true));
        double ry = ((subset & 2) ? RawLastCellCoord(io, grid.yMin, grid.ySpacing, grid.ny)
            : RawCoord(io, grid.yMin, grid.ySpacing, grid.ny, 0, grid.ny - 2, true));
        double rz = ((subset & 4) ? RawLastCellCoord(io, grid.zMin, grid.zSpacing, grid.nz)
            : RawCoord(io, grid.zMin, grid.zSpacing, grid.nz, 0, grid.nz - 2, true));
        double x = io.given(rx);
        double y = io.given(ry);
        double z = io.given(rz);
        OutAkima3(io, interp, x, y, z);
    }
}

ORACLE_CASE("IntpAkimaUniform3.construct.validation")
{
    int xBound = ValidationBound(io, 2, 3, 4);
    int yBound = ValidationBound(io, 2, 3, 4);
    int zBound = ValidationBound(io, 2, 3, 4);
    double xSpacing = ValidationSpacing(io);
    double ySpacing = ValidationSpacing(io);
    double zSpacing = ValidationSpacing(io);
    std::vector<double> F(static_cast<size_t>(xBound * yBound * zBound));
    for (auto& f : F) { f = io.lattice(-3, 3); }
    IntpAkimaUniform3<double> interp(xBound, yBound, zBound, 0.0, xSpacing, 0.0, ySpacing,
        0.0, zSpacing, F.data());
    // Cell (0, 0, 0) for either valid spacing, away from the #58 faces.
    io.outReal(interp(0.25, 0.125, 0.0625));
}

// ---------------------------------------------------------- IntpBicubic2, IntpBilinear2
//
// Both clamp only the cell index; the fractional coordinate extrapolates
// (#69, preserved), so the outside-the-domain query modes compare the
// extrapolation bit for bit. |x - xMin| / xSpacing stays far below 2^31: a
// larger index makes static_cast<int32_t> undefined (MSVC's cvttsd2si gives
// INT_MIN, clamped to cell 0, where Math.trunc gives a huge index clamped to
// the last cell), so such inputs are out of scope.

ORACLE_CASE("IntpBicubic2.evaluate")
{
    static int const modes[6] = { 0, 1, 5, 3, 2, 4 };
    int mode = modes[io.index() % 6];
    int xBound = io.integer(3, 6);
    int yBound = io.integer(3, 6);
    double xMin = io.given(RawMin(io));
    double xSpacing = io.given(RawSpacing(io));
    double yMin = io.given(RawMin(io));
    double ySpacing = io.given(RawSpacing(io));
    bool catmullRom = io.boolean();
    std::vector<double> F = RawGrid2(io, mode, xBound, yBound);
    for (auto& f : F) { f = io.given(f); }
    IntpBicubic2<double> interp(xBound, yBound, xMin, xSpacing, yMin, ySpacing, F.data(), catmullRom);
    io.outInt(interp.GetQuantity());
    io.outReal(interp.GetXMax());
    io.outReal(interp.GetYMax());
    for (int q = 0; q < 3; ++q)
    {
        double x = io.given(RawCoord(io, xMin, xSpacing, xBound, 0, xBound - 2, true));
        double y = io.given(RawCoord(io, yMin, ySpacing, yBound, 0, yBound - 2, true));
        io.outReal(interp(x, y));
        for (int xo = 0; xo <= 4; ++xo)
        {
            for (int yo = 0; yo <= 4; ++yo) { io.outReal(interp(xo, yo, x, y)); }
        }
        io.outReal(interp(-1, 0, x, y));
    }
}

ORACLE_CASE("IntpBicubic2.construct.validation")
{
    int xBound = ValidationBound(io, 1, 3, 4);
    int yBound = ValidationBound(io, 1, 3, 4);
    double xSpacing = ValidationSpacing(io);
    double ySpacing = ValidationSpacing(io);
    bool catmullRom = io.boolean();
    std::vector<double> F(static_cast<size_t>(xBound) * static_cast<size_t>(yBound));
    for (auto& f : F) { f = io.lattice(-3, 3); }
    IntpBicubic2<double> interp(xBound, yBound, 0.0, xSpacing, 0.0, ySpacing, F.data(), catmullRom);
    io.outReal(interp(0.25, 0.5));
}

ORACLE_CASE("IntpBilinear2.evaluate")
{
    static int const modes[6] = { 0, 1, 3, 5, 2, 4 };
    int mode = modes[io.index() % 6];
    int xBound = io.integer(2, 5);
    int yBound = io.integer(2, 5);
    double xMin = io.given(RawMin(io));
    double xSpacing = io.given(RawSpacing(io));
    double yMin = io.given(RawMin(io));
    double ySpacing = io.given(RawSpacing(io));
    std::vector<double> F = RawGrid2(io, mode, xBound, yBound);
    for (auto& f : F) { f = io.given(f); }
    IntpBilinear2<double> interp(xBound, yBound, xMin, xSpacing, yMin, ySpacing, F.data());
    io.outInt(interp.GetQuantity());
    io.outReal(interp.GetXMax());
    io.outReal(interp.GetYMax());
    for (int q = 0; q < 4; ++q)
    {
        double x = io.given(RawCoord(io, xMin, xSpacing, xBound, 0, xBound - 2, true));
        double y = io.given(RawCoord(io, yMin, ySpacing, yBound, 0, yBound - 2, true));
        io.outReal(interp(x, y));
        for (int xo = 0; xo <= 2; ++xo)
        {
            for (int yo = 0; yo <= 2; ++yo) { io.outReal(interp(xo, yo, x, y)); }
        }
        io.outReal(interp(0, -1, x, y));
    }
}

ORACLE_CASE("IntpBilinear2.construct.validation")
{
    int xBound = ValidationBound(io, 1, 2, 3);
    int yBound = ValidationBound(io, 1, 2, 3);
    double xSpacing = ValidationSpacing(io);
    double ySpacing = ValidationSpacing(io);
    std::vector<double> F(static_cast<size_t>(xBound) * static_cast<size_t>(yBound));
    for (auto& f : F) { f = io.lattice(-3, 3); }
    IntpBilinear2<double> interp(xBound, yBound, 0.0, xSpacing, 0.0, ySpacing, F.data());
    io.outReal(interp(0.25, 0.5));
}
