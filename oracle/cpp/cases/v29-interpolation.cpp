// Verify group 29 (interpolation): differential cases for IntpTrilinear3.h,
// IntpTricubic3.h, IntpAkimaUniform1.h, IntpAkimaNonuniform1.h (and through
// them the public surface of IntpAkima1.h), IntpBSplineUniform.h,
// IntpLinearNonuniform2.h, IntpLinearNonuniform3.h,
// IntpQuadraticNonuniform2.h, IntpThinPlateSpline2.h,
// IntpThinPlateSpline3.h, IntpSphere2.h and IntpVectorField2.h.
//
// Almost everything here is + - * / sqrt, so the default is { exact: true }.
// The C math library enters only in IntpThinPlateSpline2 (std::log in the
// kernel) and IntpSphere2::GetSphericalCoordinates (std::atan2, std::acos).
//
// MESHES. The nonuniform interpolators take a duck-typed mesh. Upstream's
// Delaunay2Mesh<T> / Delaunay3Mesh<T> number the simplices in the iteration
// order of ETManifoldMesh / TSManifoldMesh's std::unordered_map, while the
// port numbers them in sorted TriangleKey<true> / TetrahedronKey<true>
// order. Two things depend on that numbering: the default start simplex of
// GetContainingTriangle (a point on a shared edge is reported in whichever
// triangle the walk reaches first) and IntpQuadraticNonuniform2's
// EstimateDerivatives, which sums the triangle normals per vertex in triangle
// order. Cases therefore come in two flavours:
//  - the real Delaunay2Mesh<double> / Delaunay3Mesh<double>, restricted to
//    what is order independent (queries strictly inside a simplex or outside
//    the hull, the derivative constructor of IntpQuadraticNonuniform2);
//  - SortedMesh2 / SortedMesh3 below, a mesh adapter over upstream's own
//    Delaunay2<double> / Delaunay3<double> that presents the simplices in the
//    port's sorted order and starts every search at simplex 0 of that order.
//    Everything else (the triangulation, the exact barycentrics copied from
//    Delaunay2Mesh<T>::GetBarycentrics, the search itself) is upstream's.
// IntpSphere2<T> and IntpVectorField2<T> build a Delaunay2Mesh<T> member
// before their triangulation exists, and that mesh's constructor asserts
// dimension 2, so neither class can be constructed (a 'deviation' case each;
// see the group report). Each gets a "sortedMesh" case that replays its
// constructor and operator() over SortedMesh2 and is compared exactly with
// the port's class.
//
// One io draw per C++ statement: MSVC evaluates function arguments right to
// left, so every generated value goes into a named local before it is used.
#define ORACLE_FAMILY "v29-interpolation"
#include "Oracle.h"

#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/Delaunay2.h>
#include <Mathematics/Delaunay2Mesh.h>
#include <Mathematics/Delaunay3.h>
#include <Mathematics/Delaunay3Mesh.h>
#include <Mathematics/IntpAkimaNonuniform1.h>
#include <Mathematics/IntpAkimaUniform1.h>
#include <Mathematics/IntpBSplineUniform.h>
#include <Mathematics/IntpLinearNonuniform2.h>
#include <Mathematics/IntpLinearNonuniform3.h>
#include <Mathematics/IntpQuadraticNonuniform2.h>
#include <Mathematics/IntpSphere2.h>
#include <Mathematics/IntpThinPlateSpline2.h>
#include <Mathematics/IntpThinPlateSpline3.h>
#include <Mathematics/IntpTricubic3.h>
#include <Mathematics/IntpTrilinear3.h>
#include <Mathematics/IntpVectorField2.h>
#include <Mathematics/TetrahedronKey.h>
#include <Mathematics/TriangleKey.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <vector>

using namespace gte;

namespace
{
    // ---- uniform 3D grids (IntpTrilinear3, IntpTricubic3) ----------------
    // Layout: 3 bounds, then (min, spacing) per axis, then the samples in
    // lexicographic order. Every mode records the same number of doubles.
    struct Grid3
    {
        std::array<int32_t, 3> bound{};
        std::array<double, 3> min{}, spacing{};
        std::vector<double> F;
    };

    // Small integer cubic used by the polynomial-data mode; exactly
    // representable on dyadic grids, so the sample is the exact value.
    double Cubic3(double x, double y, double z)
    {
        return 2.0 + x - 3.0 * y + 0.5 * z + x * y - z * z + 0.25 * x * x * x
            - y * y * z + x * y * z;
    }

    Grid3 MakeGrid3(oracle::Ctx& io, int32_t lo, int32_t hi, bool cubicData)
    {
        Grid3 g;
        for (int32_t d = 0; d < 3; ++d)
        {
            g.bound[d] = io.integer(lo, hi);
        }
        int32_t mode = io.index() % 5;
        for (int32_t d = 0; d < 3; ++d)
        {
            double mn = 0.0, sp = 1.0;
            if (mode == 0 || mode == 3)
            {
                mn = static_cast<double>(io.rawInteger(-2, 2));
                sp = (mode == 3 && io.rawInteger(0, 1) == 1 ? 0.5 : 1.0);
            }
            else if (mode == 1)
            {
                mn = 0.25 * static_cast<double>(io.rawInteger(-8, 8));
                int32_t k = io.rawInteger(0, 2);
                sp = (k == 0 ? 0.25 : (k == 1 ? 0.5 : 2.0));
            }
            else if (mode == 2)
            {
                mn = io.raw(-5.0, 5.0);
                sp = io.raw(0.1, 3.0);
            }
            g.min[d] = io.given(mn);
            g.spacing[d] = io.given(sp);
        }
        size_t q = static_cast<size_t>(g.bound[0]) * g.bound[1] * g.bound[2];
        g.F.resize(q);
        for (int32_t s = 0, i = 0; s < g.bound[2]; ++s)
        {
            for (int32_t r = 0; r < g.bound[1]; ++r)
            {
                for (int32_t c = 0; c < g.bound[0]; ++c, ++i)
                {
                    double f;
                    if (mode == 0)
                    {
                        f = static_cast<double>(io.rawInteger(-3, 3));
                    }
                    else if (mode == 3)
                    {
                        double x = g.min[0] + g.spacing[0] * c;
                        double y = g.min[1] + g.spacing[1] * r;
                        double z = g.min[2] + g.spacing[2] * s;
                        f = (cubicData ? Cubic3(x, y, z) : 1.5 - 2.0 * x + 0.5 * y + 3.0 * z);
                    }
                    else if (mode == 4)
                    {
                        // Signed zeros and small integers.
                        int32_t k = io.rawInteger(-2, 3);
                        f = (k == 3 ? -0.0 : static_cast<double>(k));
                    }
                    else
                    {
                        f = io.raw(-10.0, 10.0);
                    }
                    g.F[i] = io.given(f);
                }
            }
        }
        return g;
    }

    // One query coordinate for an axis with 'bound' samples. Modes: uniform
    // over the domain widened by 1.5 cells, a grid node, a dyadic point in a
    // cell, exactly the maximum, far outside (#69: the cell index is clamped
    // but the fractional coordinate extrapolates), exactly the minimum or -0.
    double GridQuery(oracle::Ctx& io, int32_t bound, double mn, double sp)
    {
        double mx = mn + sp * static_cast<double>(bound - 1);
        int32_t mode = io.rawInteger(0, 6);
        double x;
        if (mode == 0)
        {
            x = io.raw(mn - 1.5 * sp, mx + 1.5 * sp);
        }
        else if (mode == 1)
        {
            x = mn + sp * static_cast<double>(io.rawInteger(0, bound - 1));
        }
        else if (mode == 2)
        {
            double cell = static_cast<double>(io.rawInteger(0, bound - 1));
            double frac = 0.125 * static_cast<double>(io.rawInteger(1, 7));
            x = mn + sp * (cell + frac);
        }
        else if (mode == 3)
        {
            x = mx;
        }
        else if (mode == 4)
        {
            double k = static_cast<double>(io.rawInteger(1, 3)) + 0.25;
            x = (io.rawInteger(0, 1) == 0 ? mn - k * sp : mx + k * sp);
        }
        else if (mode == 5)
        {
            x = mn;
        }
        else
        {
            x = (mn == 0.0 ? -0.0 : mn + sp * 0.5);
        }
        return io.given(x);
    }
}

ORACLE_CASE("IntpTrilinear3.evaluate")
{
    Grid3 g = MakeGrid3(io, 2, 4, false);
    IntpTrilinear3<double> intp(g.bound[0], g.bound[1], g.bound[2], g.min[0],
        g.spacing[0], g.min[1], g.spacing[1], g.min[2], g.spacing[2], g.F.data());
    io.outInt(intp.GetXBound());
    io.outInt(intp.GetYBound());
    io.outInt(intp.GetZBound());
    io.outInt(intp.GetQuantity());
    io.outReal(intp.GetXMin());
    io.outReal(intp.GetXMax());
    io.outReal(intp.GetXSpacing());
    io.outReal(intp.GetYMin());
    io.outReal(intp.GetYMax());
    io.outReal(intp.GetYSpacing());
    io.outReal(intp.GetZMin());
    io.outReal(intp.GetZMax());
    io.outReal(intp.GetZSpacing());
    for (int32_t q = 0; q < 3; ++q)
    {
        double x = GridQuery(io, g.bound[0], g.min[0], g.spacing[0]);
        double y = GridQuery(io, g.bound[1], g.min[1], g.spacing[1]);
        double z = GridQuery(io, g.bound[2], g.min[2], g.spacing[2]);
        io.outReal(intp(x, y, z));
        // Orders 0 and 1 are the two evaluated branches; 2 and -1 take the
        // 'default: return 0' branch of each axis.
        for (int32_t zo = -1; zo <= 2; ++zo)
        {
            for (int32_t yo = -1; yo <= 2; ++yo)
            {
                for (int32_t xo = -1; xo <= 2; ++xo)
                {
                    io.outReal(intp(xo, yo, zo, x, y, z));
                }
            }
        }
    }
}

ORACLE_CASE("IntpTricubic3.evaluate")
{
    Grid3 g = MakeGrid3(io, 4, 5, true);
    bool catmullRom = io.boolean();
    IntpTricubic3<double> intp(g.bound[0], g.bound[1], g.bound[2], g.min[0],
        g.spacing[0], g.min[1], g.spacing[1], g.min[2], g.spacing[2], g.F.data(),
        catmullRom);
    io.outInt(intp.GetXBound());
    io.outInt(intp.GetYBound());
    io.outInt(intp.GetZBound());
    io.outInt(intp.GetQuantity());
    io.outReal(intp.GetXMin());
    io.outReal(intp.GetXMax());
    io.outReal(intp.GetXSpacing());
    io.outReal(intp.GetYMin());
    io.outReal(intp.GetYMax());
    io.outReal(intp.GetYSpacing());
    io.outReal(intp.GetZMin());
    io.outReal(intp.GetZMax());
    io.outReal(intp.GetZSpacing());
    for (int32_t q = 0; q < 2; ++q)
    {
        double x = GridQuery(io, g.bound[0], g.min[0], g.spacing[0]);
        double y = GridQuery(io, g.bound[1], g.min[1], g.spacing[1]);
        double z = GridQuery(io, g.bound[2], g.min[2], g.spacing[2]);
        io.outReal(intp(x, y, z));
        // Orders 0..3 are the evaluated branches; 4 and -1 take 'default'.
        for (int32_t zo = 0; zo <= 3; ++zo)
        {
            for (int32_t yo = 0; yo <= 3; ++yo)
            {
                for (int32_t xo = 0; xo <= 3; ++xo)
                {
                    io.outReal(intp(xo, yo, zo, x, y, z));
                }
            }
        }
        io.outReal(intp(4, 0, 0, x, y, z));
        io.outReal(intp(0, -1, 0, x, y, z));
        io.outReal(intp(1, 2, 4, x, y, z));
    }
}

namespace
{
    // Throw-parity inputs for the grid constructors: bounds around the
    // minimum and spacings that are positive, zero, -0, negative or NaN.
    // The samples are drawn only when there are some, so an empty grid
    // records no samples; the C++ buffer is never empty (a null F would be
    // a second reason to throw).
    Grid3 MakeInvalidGrid3(oracle::Ctx& io, int32_t lo, int32_t hi)
    {
        Grid3 g;
        for (int32_t d = 0; d < 3; ++d)
        {
            g.bound[d] = io.integer(lo, hi);
        }
        for (int32_t d = 0; d < 3; ++d)
        {
            g.min[d] = io.lattice(-2, 2);
            int32_t k = io.rawInteger(0, 9);
            double sp = (k <= 5 ? 0.5 * static_cast<double>(k + 1) : (k == 6 ? 0.0
                : (k == 7 ? -0.0 : (k == 8 ? -1.0 : std::numeric_limits<double>::quiet_NaN()))));
            g.spacing[d] = io.given(sp);
        }
        size_t q = static_cast<size_t>(g.bound[0]) * g.bound[1] * g.bound[2];
        g.F.resize(std::max<size_t>(q, 1), 0.0);
        for (size_t i = 0; i < q; ++i)
        {
            g.F[i] = io.real(-3.0, 3.0);
        }
        return g;
    }
}

ORACLE_CASE("IntpTrilinear3.construct.invalid")
{
    Grid3 g = MakeInvalidGrid3(io, 0, 3);
    IntpTrilinear3<double> intp(g.bound[0], g.bound[1], g.bound[2], g.min[0],
        g.spacing[0], g.min[1], g.spacing[1], g.min[2], g.spacing[2], g.F.data());
    io.outInt(intp.GetQuantity());
    io.outReal(intp.GetXMax());
    io.outReal(intp.GetYMax());
    io.outReal(intp.GetZMax());
    io.outReal(intp(g.min[0] + 0.25, g.min[1] + 0.5, g.min[2] + 0.75));
}

ORACLE_CASE("IntpTricubic3.construct.invalid")
{
    Grid3 g = MakeInvalidGrid3(io, 2, 5);
    bool catmullRom = io.boolean();
    IntpTricubic3<double> intp(g.bound[0], g.bound[1], g.bound[2], g.min[0],
        g.spacing[0], g.min[1], g.spacing[1], g.min[2], g.spacing[2], g.F.data(),
        catmullRom);
    io.outInt(intp.GetQuantity());
    io.outReal(intp.GetXMax());
    io.outReal(intp.GetYMax());
    io.outReal(intp.GetZMax());
    io.outReal(intp(g.min[0] + 0.25, g.min[1] + 0.5, g.min[2] + 0.75));
}

namespace
{
    // ---- Akima 1D --------------------------------------------------------
    // The polynomials are protected; the subclass exposes their coefficients
    // (the port's subclass in the replay does the same), which covers the
    // coefficient construction directly rather than only through Evaluate.
    template <typename Base>
    class AkimaExposed : public Base
    {
    public:
        using Base::Base;
        double Coeff(int32_t i, int32_t k) { return this->mPoly[i][k]; }
    };

    // Samples F for 'q' nodes at spacing-free positions 0..q-1 (the modes
    // only shape the values). Modes: uniform, small lattice (equal
    // consecutive slopes reach every branch of ComputeDerivative), strictly
    // increasing, plateaus with jumps, exactly linear.
    std::vector<double> AkimaSamples(oracle::Ctx& io, int32_t q, int32_t mode)
    {
        std::vector<double> F(q);
        double acc = static_cast<double>(io.rawInteger(-3, 3));
        double slope = static_cast<double>(io.rawInteger(-2, 2)) * 0.5;
        for (int32_t i = 0; i < q; ++i)
        {
            double f;
            if (mode == 0)
            {
                f = io.raw(-5.0, 5.0);
            }
            else if (mode == 1)
            {
                // Lattice values with signed zeros: a -0 sample makes -0
                // slopes and coefficients, and std::max(-0, 0) (returns its
                // first argument) keeps a -0 query, which Math.max would not.
                int32_t k = io.rawInteger(-2, 3);
                f = (k == 3 ? -0.0 : static_cast<double>(k));
            }
            else if (mode == 2)
            {
                acc += io.raw(0.01, 3.0);
                f = acc;
            }
            else if (mode == 3)
            {
                if (io.rawInteger(0, 2) == 0)
                {
                    acc += static_cast<double>(io.rawInteger(-3, 3));
                }
                f = acc;
            }
            else
            {
                f = acc + slope * static_cast<double>(i);
            }
            F[i] = io.given(f);
        }
        return F;
    }

    // A query for an Akima domain [x0, x1] with nodes X (possibly uniform).
    double AkimaQuery(oracle::Ctx& io, std::vector<double> const& X)
    {
        double x0 = X.front(), x1 = X.back();
        size_t n = X.size();
        int32_t mode = io.rawInteger(0, 5);
        double x;
        if (mode == 0)
        {
            double w = x1 - x0;
            x = io.raw(x0 - 0.25 * w, x1 + 0.25 * w);
        }
        else if (mode == 1)
        {
            x = X[static_cast<size_t>(io.rawInteger(0, static_cast<int32_t>(n) - 1))];
        }
        else if (mode == 2)
        {
            size_t i = static_cast<size_t>(io.rawInteger(0, static_cast<int32_t>(n) - 2));
            double frac = 0.125 * static_cast<double>(io.rawInteger(1, 7));
            x = X[i] + frac * (X[i + 1] - X[i]);
        }
        else if (mode == 3)
        {
            x = x1;
        }
        else if (mode == 4)
        {
            x = (x0 == 0.0 ? -0.0 : x0);
        }
        else
        {
            x = (io.rawInteger(0, 1) == 0 ? x0 - 2.0 : x1 + 2.0);
        }
        return io.given(x);
    }

    template <typename Intp>
    void EmitAkima(oracle::Ctx& io, AkimaExposed<Intp>& intp, std::vector<double> const& X)
    {
        int32_t q = intp.GetQuantity();
        io.outInt(q);
        io.outReal(intp.GetXMin());
        io.outReal(intp.GetXMax());
        for (int32_t i = 0; i + 1 < q; ++i)
        {
            for (int32_t k = 0; k < 4; ++k)
            {
                io.outReal(intp.Coeff(i, k));
            }
        }
        for (int32_t j = 0; j < 5; ++j)
        {
            double x = AkimaQuery(io, X);
            io.outReal(intp(x));
            for (int32_t order = -1; order <= 4; ++order)
            {
                io.outReal(intp(order, x));
            }
        }
    }
}

ORACLE_CASE("IntpAkimaUniform1.evaluate")
{
    int32_t q = io.integer(3, 8);
    int32_t mode = io.index() % 5;
    double xMin = 0.0, spacing = 1.0;
    int32_t xmode = io.index() % 3;
    if (xmode == 1)
    {
        xMin = 0.25 * static_cast<double>(io.rawInteger(-8, 8));
        spacing = 0.25 * static_cast<double>(io.rawInteger(1, 8));
    }
    else if (xmode == 2)
    {
        xMin = io.raw(-5.0, 5.0);
        spacing = io.raw(0.05, 2.0);
    }
    xMin = io.given(xMin);
    spacing = io.given(spacing);
    std::vector<double> F = AkimaSamples(io, q, mode);
    AkimaExposed<IntpAkimaUniform1<double>> intp(q, xMin, spacing, F.data());
    std::vector<double> X(q);
    for (int32_t i = 0; i < q; ++i)
    {
        X[i] = xMin + spacing * static_cast<double>(i);
    }
    X[q - 1] = intp.GetXMax();
    io.outReal(intp.GetXSpacing());
    EmitAkima(io, intp, X);
}

ORACLE_CASE("IntpAkimaNonuniform1.evaluate")
{
    int32_t q = io.integer(3, 8);
    int32_t mode = io.index() % 5;
    int32_t xmode = io.index() % 3;
    std::vector<double> X(q);
    double x = (xmode == 2 ? io.raw(-5.0, 5.0) : static_cast<double>(io.rawInteger(-3, 3)));
    for (int32_t i = 0; i < q; ++i)
    {
        X[i] = io.given(x);
        if (xmode == 0)
        {
            x += static_cast<double>(io.rawInteger(1, 3));
        }
        else if (xmode == 1)
        {
            x += 0.125 * static_cast<double>(io.rawInteger(1, 16));
        }
        else
        {
            x += io.raw(0.01, 2.0);
        }
    }
    std::vector<double> F = AkimaSamples(io, q, mode);
    AkimaExposed<IntpAkimaNonuniform1<double>> intp(q, X.data(), F.data());
    EmitAkima(io, intp, X);
}

ORACLE_CASE("IntpAkimaUniform1.construct.invalid")
{
    int32_t q = io.integer(0, 4);
    double xMin = io.lattice(-2, 2);
    int32_t k = io.rawInteger(0, 4);
    double spacing = io.given(k == 0 ? 0.0 : (k == 1 ? -0.0 : (k == 2 ? -0.5 : 0.5 * k)));
    std::vector<double> F(std::max(q, 1), 0.0);
    for (int32_t i = 0; i < q; ++i)
    {
        F[i] = io.lattice(-3, 3);
    }
    IntpAkimaUniform1<double> intp(q, xMin, spacing, F.data());
    io.outReal(intp.GetXMax());
    io.outReal(intp(xMin + 0.75));
}

ORACLE_CASE("IntpAkimaNonuniform1.construct.invalid")
{
    // Modes: too few samples (the base-class assert), strictly increasing
    // (valid), one repeated node, one decreasing node.
    int32_t mode = io.index() % 4;
    int32_t q = io.integer(mode == 0 ? 1 : 3, mode == 0 ? 2 : 6);
    std::vector<double> X(q), F(q);
    int32_t bad = (q > 1 ? io.rawInteger(1, q - 1) : 0);
    double x = 0.0;
    for (int32_t i = 0; i < q; ++i)
    {
        double step = 0.5 * static_cast<double>(io.rawInteger(1, 4));
        x += (i == bad && mode == 2 ? 0.0 : (i == bad && mode == 3 ? -step : step));
        X[i] = io.given(x);
    }
    for (int32_t i = 0; i < q; ++i)
    {
        F[i] = io.lattice(-3, 3);
    }
    IntpAkimaNonuniform1<double> intp(q, X.data(), F.data());
    io.outReal(intp.GetXMax());
    io.outReal(intp(X[0] + 0.5));
}

namespace
{
    // ---- IntpBSplineUniform ----------------------------------------------
    // A Controls adapter over a lexicographic array (dimension 0 fastest)
    // with both accessor forms: the tuple form for the general-dimension
    // code and the component form for the 1/2/3-dimensional specializations.
    struct GridControls
    {
        using Type = double;
        std::vector<int32_t> size;
        std::vector<double> data;

        int32_t GetSize(int32_t dim) const { return size[dim]; }

        double operator()(int32_t const* tuple) const
        {
            size_t n = size.size();
            int32_t index = tuple[n - 1];
            for (size_t d = n - 1; d-- > 0; )
            {
                index = index * size[d] + tuple[d];
            }
            return data[index];
        }

        double operator()(int32_t i0) const { return data[i0]; }
        double operator()(int32_t i0, int32_t i1) const { return data[i0 + size[0] * i1]; }
        double operator()(int32_t i0, int32_t i1, int32_t i2) const
        {
            return data[i0 + size[0] * (i1 + size[1] * i2)];
        }
    };

    // Layout: sizes (ints), then the control values. Modes: uniform,
    // lattice with signed zeros, constant.
    GridControls MakeControls(oracle::Ctx& io, std::vector<int32_t> const& sizes)
    {
        GridControls c;
        c.size = sizes;
        size_t q = 1;
        for (int32_t s : sizes)
        {
            q *= static_cast<size_t>(s);
        }
        c.data.resize(q);
        int32_t mode = io.index() % 3;
        double constant = static_cast<double>(io.rawInteger(-3, 3));
        for (size_t i = 0; i < q; ++i)
        {
            double v;
            if (mode == 0)
            {
                v = io.raw(-5.0, 5.0);
            }
            else if (mode == 1)
            {
                int32_t k = io.rawInteger(-2, 3);
                v = (k == 3 ? -0.0 : static_cast<double>(k));
            }
            else
            {
                v = constant;
            }
            c.data[i] = io.given(v);
        }
        return c;
    }

    // The zero control point passed to the constructor; its sign is the seed
    // of every accumulation.
    double MakeCTZero(oracle::Ctx& io)
    {
        return io.given(io.rawInteger(0, 1) == 0 ? 0.0 : -0.0);
    }

    // A parameter for an axis with 'numControls' controls and 'degree'.
    // Modes: uniform over [tmin - 1, tmax + 1], tmin, tmax, a knot t with
    // dsdt*(t - tmin) near an integer, a dyadic interior value, far outside.
    double BSplineQuery(oracle::Ctx& io, int32_t numControls, int32_t degree)
    {
        double tmin = -0.5, tmax = static_cast<double>(numControls) - 0.5;
        int32_t mode = io.rawInteger(0, 5);
        double t;
        if (mode == 0)
        {
            t = io.raw(tmin - 1.0, tmax + 1.0);
        }
        else if (mode == 1)
        {
            t = tmin;
        }
        else if (mode == 2)
        {
            t = tmax;
        }
        else if (mode == 3)
        {
            double dsdt = (static_cast<double>(numControls) - static_cast<double>(degree))
                / (tmax - tmin);
            double k = static_cast<double>(io.rawInteger(1, numControls - degree - 1));
            t = tmin + k / dsdt;
        }
        else if (mode == 4)
        {
            t = 0.125 * static_cast<double>(io.rawInteger(-3, 8 * numControls - 5));
        }
        else
        {
            t = (io.rawInteger(0, 1) == 0 ? tmin - 3.0 : tmax + 3.0);
        }
        return io.given(t);
    }
}

ORACLE_CASE("IntpBSplineUniform1.evaluate")
{
    int32_t degree = io.integer(1, 5);
    int32_t n = io.integer(degree + 2, degree + 6);
    int32_t cacheMode = io.integer(0, 2);
    GridControls controls = MakeControls(io, { n });
    double ctZero = MakeCTZero(io);
    IntpBSplineUniform<double, GridControls, 1> intp(degree, controls, ctZero, cacheMode);
    io.outInt(intp.GetDegree(0));
    io.outInt(intp.GetNumControls(0));
    io.outReal(intp.GetTMin(0));
    io.outReal(intp.GetTMax(0));
    io.outInt(intp.GetCacheMode());
    // Several evaluations on one object: ON_DEMAND_CACHING fills its tensor
    // lazily, so later calls read what earlier calls cached.
    for (int32_t j = 0; j < 6; ++j)
    {
        double t = BSplineQuery(io, n, degree);
        for (int32_t order = -1; order <= degree + 1; ++order)
        {
            io.outReal(intp.Evaluate({ order }, { t }));
        }
    }
}

ORACLE_CASE("IntpBSplineUniform2.evaluate")
{
    std::array<int32_t, 2> degrees{}, sizes{};
    for (int32_t d = 0; d < 2; ++d)
    {
        degrees[d] = io.integer(1, 3);
    }
    for (int32_t d = 0; d < 2; ++d)
    {
        sizes[d] = io.integer(degrees[d] + 2, degrees[d] + 4);
    }
    int32_t cacheMode = io.integer(0, 2);
    GridControls controls = MakeControls(io, { sizes[0], sizes[1] });
    double ctZero = MakeCTZero(io);
    IntpBSplineUniform<double, GridControls, 2> intp(degrees, controls, ctZero, cacheMode);
    for (int32_t d = 0; d < 2; ++d)
    {
        io.outInt(intp.GetDegree(d));
        io.outInt(intp.GetNumControls(d));
        io.outReal(intp.GetTMin(d));
        io.outReal(intp.GetTMax(d));
    }
    io.outInt(intp.GetCacheMode());
    for (int32_t j = 0; j < 4; ++j)
    {
        double t0 = BSplineQuery(io, sizes[0], degrees[0]);
        double t1 = BSplineQuery(io, sizes[1], degrees[1]);
        for (int32_t o1 = -1; o1 <= degrees[1] + 1; ++o1)
        {
            for (int32_t o0 = -1; o0 <= degrees[0] + 1; ++o0)
            {
                io.outReal(intp.Evaluate({ o0, o1 }, { t0, t1 }));
            }
        }
    }
}

ORACLE_CASE("IntpBSplineUniform3.evaluate")
{
    std::array<int32_t, 3> degrees{}, sizes{};
    for (int32_t d = 0; d < 3; ++d)
    {
        degrees[d] = io.integer(1, 2);
    }
    for (int32_t d = 0; d < 3; ++d)
    {
        sizes[d] = io.integer(degrees[d] + 2, degrees[d] + 3);
    }
    int32_t cacheMode = io.integer(0, 2);
    GridControls controls = MakeControls(io, { sizes[0], sizes[1], sizes[2] });
    double ctZero = MakeCTZero(io);
    IntpBSplineUniform<double, GridControls, 3> intp(degrees, controls, ctZero, cacheMode);
    for (int32_t d = 0; d < 3; ++d)
    {
        io.outInt(intp.GetDegree(d));
        io.outInt(intp.GetNumControls(d));
        io.outReal(intp.GetTMin(d));
        io.outReal(intp.GetTMax(d));
    }
    io.outInt(intp.GetCacheMode());
    for (int32_t j = 0; j < 3; ++j)
    {
        std::array<double, 3> t{};
        for (int32_t d = 0; d < 3; ++d)
        {
            t[d] = BSplineQuery(io, sizes[d], degrees[d]);
        }
        for (int32_t o2 = 0; o2 <= degrees[2] + 1; ++o2)
        {
            for (int32_t o1 = 0; o1 <= degrees[1] + 1; ++o1)
            {
                for (int32_t o0 = 0; o0 <= degrees[0] + 1; ++o0)
                {
                    io.outReal(intp.Evaluate({ o0, o1, o2 }, t));
                }
            }
        }
        io.outReal(intp.Evaluate({ -1, 0, 0 }, t));
    }
}

namespace
{
    // The general-dimension code (IntpBSplineUniformShared::Evaluate*Caching)
    // through the run-time-dimension class IntpBSplineUniform<Real,
    // Controls, 0>, which is the class the port implements.
    void RunGeneralBSpline(oracle::Ctx& io, int32_t N, int32_t minDegree, int32_t maxDegree,
        int32_t maxExtra)
    {
        std::vector<int32_t> degrees(N), sizes(N);
        for (int32_t d = 0; d < N; ++d)
        {
            degrees[d] = io.integer(minDegree, maxDegree);
        }
        for (int32_t d = 0; d < N; ++d)
        {
            sizes[d] = io.integer(degrees[d] + 2, degrees[d] + 2 + maxExtra);
        }
        int32_t cacheMode = io.integer(0, 2);
        GridControls controls = MakeControls(io, sizes);
        double ctZero = MakeCTZero(io);
        IntpBSplineUniform<double, GridControls> intp(degrees, controls, ctZero, cacheMode);
        for (int32_t d = 0; d < N; ++d)
        {
            io.outInt(intp.GetDegree(d));
            io.outInt(intp.GetNumControls(d));
            io.outReal(intp.GetTMin(d));
            io.outReal(intp.GetTMax(d));
        }
        io.outInt(intp.GetCacheMode());
        for (int32_t j = 0; j < 3; ++j)
        {
            std::vector<double> t(N);
            for (int32_t d = 0; d < N; ++d)
            {
                t[d] = BSplineQuery(io, sizes[d], degrees[d]);
            }
            // Every order tuple with entries in [0, degree + 1] (the last
            // value takes the out-of-range early return), odometer order.
            std::vector<int32_t> order(N, 0);
            for (;;)
            {
                io.outReal(intp.Evaluate(order, t));
                int32_t d = 0;
                for (; d < N; ++d)
                {
                    if (++order[d] <= degrees[d] + 1)
                    {
                        break;
                    }
                    order[d] = 0;
                }
                if (d == N)
                {
                    break;
                }
            }
            // The size guard of the run-time-dimension Evaluate: too short
            // an order or t vector returns the zero control point.
            std::vector<int32_t> shortOrder(N - 1, 0);
            std::vector<double> shortT(N - 1, 0.0);
            io.outReal(intp.Evaluate(shortOrder, t));
            io.outReal(intp.Evaluate(std::vector<int32_t>(N, 0), shortT));
            std::vector<int32_t> negative(N, 0);
            negative[N - 1] = -1;
            io.outReal(intp.Evaluate(negative, t));
        }
    }
}

ORACLE_CASE("IntpBSplineUniform.evaluate.general")
{
    // N = 1..3 with degrees 1..3 (1..2 for N = 3), and N = 4 with degree 1.
    int32_t N = 1 + io.index() % 4;
    if (N == 1)
    {
        RunGeneralBSpline(io, 1, 1, 4, 3);
    }
    else if (N == 2)
    {
        RunGeneralBSpline(io, 2, 1, 3, 2);
    }
    else if (N == 3)
    {
        RunGeneralBSpline(io, 3, 1, 2, 1);
    }
    else
    {
        RunGeneralBSpline(io, 4, 1, 1, 0);
    }
}

ORACLE_CASE("IntpBSplineUniform.evaluate.fixedN")
{
    // The compile-time-dimension template IntpBSplineUniform<Real, Controls,
    // 4> runs the same shared code without the size guard.
    std::array<int32_t, 4> degrees{}, sizes{};
    for (int32_t d = 0; d < 4; ++d)
    {
        degrees[d] = io.integer(1, 2);
    }
    for (int32_t d = 0; d < 4; ++d)
    {
        sizes[d] = io.integer(degrees[d] + 2, degrees[d] + 2);
    }
    int32_t cacheMode = io.integer(0, 2);
    GridControls controls = MakeControls(io, { sizes[0], sizes[1], sizes[2], sizes[3] });
    double ctZero = MakeCTZero(io);
    IntpBSplineUniform<double, GridControls, 4> intp(degrees, controls, ctZero, cacheMode);
    for (int32_t j = 0; j < 3; ++j)
    {
        std::array<double, 4> t{};
        for (int32_t d = 0; d < 4; ++d)
        {
            t[d] = BSplineQuery(io, sizes[d], degrees[d]);
        }
        std::array<int32_t, 4> order{};
        for (int32_t d = 0; d < 4; ++d)
        {
            order[d] = io.integer(0, degrees[d]);
        }
        io.outReal(intp.Evaluate(order, t));
        io.outReal(intp.Evaluate({ 0, 0, 0, 0 }, t));
    }
}

// Degree 0 is not covered. Upstream's ComputePowers resizes powerDSDT to
// degree + 1 = 1 element and then writes powerDSDT[1] = dsdt (#135, fixed in
// the port). On this MSVC build the write is not benign: constructing a
// single IntpBSplineUniform<double, Controls, 1> of degree 0 terminates the
// process with STATUS_HEAP_CORRUPTION (0xC0000374) when the heap is next
// checked, so no golden record can exist for any degree-0 path (every
// specialization and the general class call ComputePowers in their
// constructors). The port's guard is pinned by test/IntpBSplineUniform.test.ts.

ORACLE_CASE("IntpBSplineUniform.construct.invalid")
{
    // numControls <= degree + 1 throws in each of the four classes.
    int32_t which = io.integer(0, 3);
    int32_t degree = io.integer(1, 3);
    int32_t n = io.integer(degree, degree + 2);
    int32_t n2 = io.integer(degree + 2, degree + 3);
    int32_t cacheMode = io.integer(0, 2);
    std::vector<int32_t> sizes = (which == 0 || which == 3)
        ? std::vector<int32_t>{ n } : (which == 1 ? std::vector<int32_t>{ n2, n }
        : std::vector<int32_t>{ n2, n2, n });
    GridControls controls = MakeControls(io, sizes);
    double value = 0.0;
    if (which == 0)
    {
        IntpBSplineUniform<double, GridControls, 1> intp(degree, controls, 0.0, cacheMode);
        value = intp.Evaluate({ 0 }, { 0.25 });
    }
    else if (which == 1)
    {
        IntpBSplineUniform<double, GridControls, 2> intp({ degree, degree }, controls, 0.0,
            cacheMode);
        value = intp.Evaluate({ 0, 0 }, { 0.25, 0.25 });
    }
    else if (which == 2)
    {
        IntpBSplineUniform<double, GridControls, 3> intp({ degree, degree, degree }, controls,
            0.0, cacheMode);
        value = intp.Evaluate({ 0, 0, 0 }, { 0.25, 0.25, 0.25 });
    }
    else
    {
        IntpBSplineUniform<double, GridControls> intp({ degree }, controls, 0.0, cacheMode);
        value = intp.Evaluate({ 0 }, { 0.25 });
    }
    io.outReal(value);
}

namespace
{
    // ---- exact predicates (generator probes only) -------------------------
    using Exact = BSNumber<UIntegerAP32>;

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

    // Upstream Delaunay2<T>/Delaunay3<T> classify the seed simplex with
    // IntrinsicsVector2/3 and a hardcoded epsilon of 0 (finding #391, fixed
    // in the port with exact predicates), and Delaunay3<T> never detects a
    // repeated vertex (#283). These probes (from v09-compgeom) accept exactly
    // the point sets on which upstream's seed is the port's seed and, in 3D,
    // no vertex repeats, so both sides build the same triangulation.
    bool Sound2(std::vector<Vector2<double>> const& pts)
    {
        IntrinsicsVector2<double> info(static_cast<int32_t>(pts.size()), pts.data(), 0.0);
        if (info.dimension != 2)
        {
            return false;
        }
        int32_t sign = ExactOrient2(pts[info.extreme[0]], pts[info.extreme[1]],
            pts[info.extreme[2]]);
        return sign != 0 && info.extremeCCW == (sign > 0);
    }

    bool Sound3(std::vector<Vector3<double>> const& pts)
    {
        for (size_t i = 1; i < pts.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                if (pts[i] == pts[j])
                {
                    return false;
                }
            }
        }
        IntrinsicsVector3<double> info(static_cast<int32_t>(pts.size()), pts.data(), 0.0);
        if (info.dimension != 3)
        {
            return false;
        }
        int32_t sign = ExactOrient3(pts[info.extreme[0]], pts[info.extreme[1]],
            pts[info.extreme[2]], pts[info.extreme[3]]);
        return sign != 0 && info.extremeCCW == (sign > 0);
    }

    std::array<std::array<int32_t, 2>, 12> const gCircle25
    { {
        { 5, 0 }, { -5, 0 }, { 0, 5 }, { 0, -5 }, { 3, 4 }, { 3, -4 },
        { -3, 4 }, { -3, -4 }, { 4, 3 }, { 4, -3 }, { -4, 3 }, { -4, -3 }
    } };

    // Planar point sets, recorded as n then the coordinates. Modes: lattice
    // [-3,3] (collinear triples, cocircular quadruples, duplicates),
    // uniform, cocircular lattice points about a lattice center, a
    // perturbed grid. Capped rejection on Sound2 with a fixed fallback.
    std::vector<Vector2<double>> MakePoints2(oracle::Ctx& io, int32_t nmin, int32_t nmax)
    {
        int32_t n = io.integer(nmin, nmax);
        int32_t mode = io.index() % 4;
        std::vector<Vector2<double>> pts(n);
        bool sound = false;
        for (int32_t attempt = 0; attempt < 64 && !sound; ++attempt)
        {
            double cx = static_cast<double>(io.rawInteger(-2, 2));
            double cy = static_cast<double>(io.rawInteger(-2, 2));
            for (int32_t i = 0; i < n; ++i)
            {
                if (mode == 0)
                {
                    pts[i][0] = static_cast<double>(io.rawInteger(-3, 3));
                    pts[i][1] = static_cast<double>(io.rawInteger(-3, 3));
                }
                else if (mode == 1)
                {
                    pts[i][0] = io.raw(-4.0, 4.0);
                    pts[i][1] = io.raw(-4.0, 4.0);
                }
                else if (mode == 2)
                {
                    int32_t k = io.rawInteger(0, 12);
                    pts[i][0] = cx + (k < 12 ? static_cast<double>(gCircle25[k][0]) : 0.0);
                    pts[i][1] = cy + (k < 12 ? static_cast<double>(gCircle25[k][1]) : 0.0);
                }
                else
                {
                    pts[i][0] = static_cast<double>(i % 3) + 0.25 * static_cast<double>(io.rawInteger(-1, 1));
                    pts[i][1] = static_cast<double>(i / 3) + 0.25 * static_cast<double>(io.rawInteger(-1, 1));
                }
            }
            sound = Sound2(pts);
        }
        if (!sound)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                double s = static_cast<double>(i);
                pts[i] = { s, s * s * 0.25 };
            }
        }
        for (int32_t i = 0; i < n; ++i)
        {
            io.givenVec(pts[i]);
        }
        return pts;
    }

    std::vector<Vector3<double>> MakePoints3(oracle::Ctx& io, int32_t nmin, int32_t nmax)
    {
        int32_t n = io.integer(nmin, nmax);
        int32_t mode = io.index() % 3;
        std::vector<Vector3<double>> pts(n);
        bool sound = false;
        for (int32_t attempt = 0; attempt < 64 && !sound; ++attempt)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    pts[i][j] = (mode == 1 ? io.raw(-3.0, 3.0)
                        : static_cast<double>(io.rawInteger(mode == 0 ? -2 : -1, mode == 0 ? 2 : 1)));
                }
            }
            sound = Sound3(pts);
        }
        if (!sound)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                double s = static_cast<double>(i);
                pts[i] = { s, s * s, s * s * s };
            }
        }
        for (int32_t i = 0; i < n; ++i)
        {
            io.givenVec(pts[i]);
        }
        return pts;
    }

    // Sample values at the vertices. Modes: uniform, lattice with signed
    // zeros, affine in (x, y[, z]) with small integer coefficients.
    template <int32_t N>
    std::vector<double> MakeValues(oracle::Ctx& io, std::vector<Vector<N, double>> const& pts)
    {
        int32_t mode = io.index() % 3;
        double c0 = static_cast<double>(io.rawInteger(-3, 3));
        std::array<double, 3> c{};
        for (int32_t j = 0; j < 3; ++j)
        {
            c[j] = 0.5 * static_cast<double>(io.rawInteger(-4, 4));
        }
        std::vector<double> F(pts.size());
        for (size_t i = 0; i < pts.size(); ++i)
        {
            double f;
            if (mode == 0)
            {
                f = io.raw(-5.0, 5.0);
            }
            else if (mode == 1)
            {
                int32_t k = io.rawInteger(-2, 3);
                f = (k == 3 ? -0.0 : static_cast<double>(k));
            }
            else
            {
                f = c0;
                for (int32_t j = 0; j < N; ++j)
                {
                    f += c[j] * pts[i][j];
                }
            }
            F[i] = io.given(f);
        }
        return F;
    }
}

namespace
{
    // ---- mesh adapters -----------------------------------------------------
    // SortedMesh2: upstream's Delaunay2<double>, with its triangles presented
    // in the port's order (sorted by TriangleKey<true> of the stored vertex
    // tuple, i.e. the ETManifoldMesh map key) and every containment search
    // started at triangle 0 of that order. Stored tuples are unchanged.
    // GetBarycentrics is Delaunay2Mesh<T>::GetBarycentrics verbatim.
    class SortedMesh2
    {
    public:
        SortedMesh2(Delaunay2<double> const& del)
            :
            mDel(&del)
        {
            auto const& ind = del.GetIndices();
            auto const& adj = del.GetAdjacencies();
            size_t nt = del.GetNumTriangles();
            mOrder.resize(nt);
            std::iota(mOrder.begin(), mOrder.end(), size_t(0));
            std::sort(mOrder.begin(), mOrder.end(), [&ind](size_t a, size_t b)
            {
                return TriangleKey<true>(ind[3 * a], ind[3 * a + 1], ind[3 * a + 2])
                    < TriangleKey<true>(ind[3 * b], ind[3 * b + 1], ind[3 * b + 2]);
            });
            mRank.resize(nt);
            for (size_t r = 0; r < nt; ++r)
            {
                mRank[mOrder[r]] = r;
            }
            mIndices.resize(3 * nt);
            mAdjacencies.resize(3 * nt);
            for (size_t r = 0; r < nt; ++r)
            {
                for (size_t j = 0; j < 3; ++j)
                {
                    mIndices[3 * r + j] = ind[3 * mOrder[r] + j];
                    int32_t a = adj[3 * mOrder[r] + j];
                    mAdjacencies[3 * r + j] = (a < 0 ? -1 : static_cast<int32_t>(mRank[a]));
                }
            }
        }

        size_t GetNumVertices() const { return mDel->GetNumVertices(); }
        size_t GetNumTriangles() const { return mOrder.size(); }
        Vector2<double> const* GetVertices() const { return mDel->GetVertices(); }
        int32_t const* GetIndices() const { return mIndices.data(); }
        size_t GetInvalidIndex() const { return Delaunay2<double>::negOne; }

        size_t GetContainingTriangle(Vector2<double> const& P) const
        {
            Delaunay2<double>::SearchInfo info{};
            info.initialTriangle = mOrder[0];
            size_t t = mDel->GetContainingTriangle(P, info);
            return t == Delaunay2<double>::negOne ? t : mRank[t];
        }

        bool GetIndices(size_t t, std::array<int32_t, 3>& indices) const
        {
            if (t < mOrder.size())
            {
                for (size_t j = 0; j < 3; ++j)
                {
                    indices[j] = mIndices[3 * t + j];
                }
                return true;
            }
            return false;
        }

        bool GetAdjacencies(size_t t, std::array<int32_t, 3>& adjacencies) const
        {
            if (t < mOrder.size())
            {
                for (size_t j = 0; j < 3; ++j)
                {
                    adjacencies[j] = mAdjacencies[3 * t + j];
                }
                return true;
            }
            return false;
        }

        bool GetVertices(size_t t, std::array<Vector2<double>, 3>& vertices) const
        {
            std::array<int32_t, 3> indices = { 0, 0, 0 };
            if (GetIndices(t, indices))
            {
                for (size_t i = 0; i < 3; ++i)
                {
                    vertices[i] = mDel->GetVertices()[indices[i]];
                }
                return true;
            }
            for (auto& vertex : vertices)
            {
                vertex.MakeZero();
            }
            return false;
        }

        bool GetBarycentrics(size_t t, Vector2<double> const& P, std::array<double, 3>& bary) const
        {
            using Rational = BSRational<UIntegerAP32>;
            std::array<int32_t, 3> indices = { 0, 0, 0 };
            if (GetIndices(t, indices))
            {
                std::array<Vector2<Rational>, 3> rtV{};
                for (size_t i = 0; i < 3; ++i)
                {
                    auto const& V = mDel->GetVertices()[indices[i]];
                    for (int32_t j = 0; j < 2; ++j)
                    {
                        rtV[i][j] = static_cast<Rational>(V[j]);
                    }
                }
                Vector2<Rational> rtP{ P[0], P[1] };
                std::array<Rational, 3> rtBary{};
                if (ComputeBarycentrics(rtP, rtV[0], rtV[1], rtV[2], rtBary))
                {
                    for (size_t i = 0; i < 3; ++i)
                    {
                        bary[i] = static_cast<double>(rtBary[i]);
                    }
                    return true;
                }
            }
            for (auto& b : bary)
            {
                b = 0.0;
            }
            return false;
        }

    private:
        Delaunay2<double> const* mDel;
        std::vector<size_t> mOrder, mRank;
        std::vector<int32_t> mIndices, mAdjacencies;
    };

    // SortedMesh3: the same for Delaunay3<double>, TetrahedronKey<true>.
    class SortedMesh3
    {
    public:
        SortedMesh3(Delaunay3<double> const& del)
            :
            mDel(&del)
        {
            auto const& ind = del.GetIndices();
            size_t nt = del.GetNumTetrahedra();
            mOrder.resize(nt);
            std::iota(mOrder.begin(), mOrder.end(), size_t(0));
            std::sort(mOrder.begin(), mOrder.end(), [&ind](size_t a, size_t b)
            {
                return TetrahedronKey<true>(ind[4 * a], ind[4 * a + 1], ind[4 * a + 2], ind[4 * a + 3])
                    < TetrahedronKey<true>(ind[4 * b], ind[4 * b + 1], ind[4 * b + 2], ind[4 * b + 3]);
            });
            mRank.resize(nt);
            for (size_t r = 0; r < nt; ++r)
            {
                mRank[mOrder[r]] = r;
            }
        }

        size_t GetInvalidIndex() const { return Delaunay3<double>::negOne; }
        size_t Upstream(size_t r) const { return mOrder[r]; }

        size_t GetContainingTetrahedron(Vector3<double> const& P) const
        {
            Delaunay3<double>::SearchInfo info{};
            info.initialTetrahedron = mOrder[0];
            size_t t = mDel->GetContainingTetrahedron(P, info);
            return t == Delaunay3<double>::negOne ? t : mRank[t];
        }

        bool GetIndices(size_t t, std::array<int32_t, 4>& indices) const
        {
            return t < mOrder.size() && mDel->GetIndices(mOrder[t], indices);
        }

        bool GetBarycentrics(size_t t, Vector3<double> const& P, std::array<double, 4>& bary) const
        {
            Delaunay3Mesh<double> mesh(*mDel);
            if (t < mOrder.size())
            {
                return mesh.GetBarycentrics(mOrder[t], P, bary);
            }
            bary.fill(0.0);
            return false;
        }

    private:
        Delaunay3<double> const* mDel;
        std::vector<size_t> mOrder, mRank;
    };
}

namespace
{
    // ---- query points on a triangulation -----------------------------------
    // Dyadic barycentric combinations of the vertices of a random simplex:
    // strictly inside (all weights k/8 > 0), on an edge (or face), at a
    // vertex; plus far outside and uniform over a box covering the hull.
    // With 'strictOnly' the query is redrawn (at most 32 times, then a far
    // outside point) until it is strictly inside a triangle or outside the
    // hull, where the containing triangle does not depend on the search's
    // start (see the file comment).
    Vector2<double> RawQuery2(oracle::Ctx& io, SortedMesh2 const& mesh, int32_t mode)
    {
        int32_t const* ind = mesh.GetIndices();
        auto const* V = mesh.GetVertices();
        size_t t = static_cast<size_t>(io.rawInteger(0, static_cast<int32_t>(mesh.GetNumTriangles()) - 1));
        Vector2<double> V0 = V[ind[3 * t]], V1 = V[ind[3 * t + 1]], V2 = V[ind[3 * t + 2]];
        if (mode == 0)
        {
            double a = static_cast<double>(io.rawInteger(1, 6));
            double b = static_cast<double>(io.rawInteger(1, 7 - static_cast<int32_t>(a)));
            double w0 = a / 8.0, w1 = b / 8.0, w2 = 1.0 - w0 - w1;
            return w0 * V0 + w1 * V1 + w2 * V2;
        }
        if (mode == 1)
        {
            double w0 = static_cast<double>(io.rawInteger(1, 7)) / 8.0;
            return w0 * V0 + (1.0 - w0) * V1;
        }
        if (mode == 2)
        {
            return V0;
        }
        if (mode == 3)
        {
            double s = (io.rawInteger(0, 1) == 0 ? -1.0 : 1.0);
            return Vector2<double>{ s * io.raw(8.0, 10.0), io.raw(-10.0, 10.0) };
        }
        return Vector2<double>{ io.raw(-7.0, 7.0), io.raw(-7.0, 7.0) };
    }

    Vector2<double> MakeQuery2(oracle::Ctx& io, Delaunay2<double> const& del, bool strictOnly)
    {
        SortedMesh2 mesh(del);
        for (int32_t attempt = 0; attempt < 32; ++attempt)
        {
            int32_t mode = io.rawInteger(0, 4);
            Vector2<double> P = RawQuery2(io, mesh, mode);
            if (!strictOnly)
            {
                return io.givenVec(P);
            }
            size_t t = mesh.GetContainingTriangle(P);
            std::array<double, 3> bary{};
            if (t == mesh.GetInvalidIndex()
                || (mesh.GetBarycentrics(t, P, bary) && bary[0] > 0.0 && bary[1] > 0.0 && bary[2] > 0.0))
            {
                return io.givenVec(P);
            }
        }
        return io.givenVec(Vector2<double>{ 20.0, 20.0 });
    }

    Vector3<double> RawQuery3(oracle::Ctx& io, Delaunay3<double> const& del, SortedMesh3 const& mesh, int32_t mode)
    {
        auto const& ind = del.GetIndices();
        auto const* V = del.GetVertices();
        size_t t = mesh.Upstream(static_cast<size_t>(io.rawInteger(0, static_cast<int32_t>(del.GetNumTetrahedra()) - 1)));
        std::array<Vector3<double>, 4> W{};
        for (size_t j = 0; j < 4; ++j)
        {
            W[j] = V[ind[4 * t + j]];
        }
        if (mode <= 2)
        {
            // mode 0: inside (4 positive weights), 1: on a face, 2: on an edge.
            int32_t numPositive = 4 - mode;
            std::array<double, 4> w{ 0.0, 0.0, 0.0, 0.0 };
            int32_t remaining = 8;
            for (int32_t j = 0; j + 1 < numPositive; ++j)
            {
                int32_t k = io.rawInteger(1, remaining - (numPositive - 1 - j));
                w[j] = static_cast<double>(k);
                remaining -= k;
            }
            w[numPositive - 1] = static_cast<double>(remaining);
            Vector3<double> P{ 0.0, 0.0, 0.0 };
            for (size_t j = 0; j < 4; ++j)
            {
                P += (w[j] / 8.0) * W[j];
            }
            return P;
        }
        if (mode == 3)
        {
            return W[0];
        }
        if (mode == 4)
        {
            return Vector3<double>{ io.raw(5.0, 8.0), io.raw(-8.0, 8.0), io.raw(-8.0, 8.0) };
        }
        return Vector3<double>{ io.raw(-3.5, 3.5), io.raw(-3.5, 3.5), io.raw(-3.5, 3.5) };
    }

    Vector3<double> MakeQuery3(oracle::Ctx& io, Delaunay3<double> const& del, bool strictOnly)
    {
        SortedMesh3 mesh(del);
        for (int32_t attempt = 0; attempt < 32; ++attempt)
        {
            int32_t mode = io.rawInteger(0, 5);
            Vector3<double> P = RawQuery3(io, del, mesh, mode);
            if (!strictOnly)
            {
                return io.givenVec(P);
            }
            size_t t = mesh.GetContainingTetrahedron(P);
            std::array<double, 4> bary{};
            if (t == mesh.GetInvalidIndex()
                || (mesh.GetBarycentrics(t, P, bary) && bary[0] > 0.0 && bary[1] > 0.0
                    && bary[2] > 0.0 && bary[3] > 0.0))
            {
                return io.givenVec(P);
            }
        }
        return io.givenVec(Vector3<double>{ 20.0, 20.0, 20.0 });
    }

    // A mesh whose per-triangle GetIndices always fails while the search
    // and the barycentrics work: the input on which IntpLinearNonuniform2
    // discards the failure flag (#135).
    class NoIndicesMesh2 : public SortedMesh2
    {
    public:
        using SortedMesh2::SortedMesh2;
        bool GetIndices(size_t, std::array<int32_t, 3>&) const { return false; }
    };

    class NoIndicesMesh3 : public SortedMesh3
    {
    public:
        using SortedMesh3::SortedMesh3;
        bool GetIndices(size_t, std::array<int32_t, 4>&) const { return false; }
    };
}

ORACLE_CASE("IntpLinearNonuniform2.evaluate")
{
    // The real Delaunay2Mesh<double>; queries strictly inside or outside.
    auto pts = MakePoints2(io, 3, 9);
    auto F = MakeValues<2>(io, pts);
    Delaunay2<double> del{};
    del(pts);
    Delaunay2Mesh<double> mesh(del);
    IntpLinearNonuniform2<double, Delaunay2Mesh<double>> intp(mesh, F.data());
    for (int32_t q = 0; q < 6; ++q)
    {
        Vector2<double> P = MakeQuery2(io, del, true);
        double f = 0.0;
        bool valid = intp(P, f);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
        }
    }
}

ORACLE_CASE("IntpLinearNonuniform2.evaluate.sortedMesh")
{
    // Every query mode, including edges and vertices shared by several
    // triangles, over SortedMesh2.
    auto pts = MakePoints2(io, 3, 9);
    auto F = MakeValues<2>(io, pts);
    Delaunay2<double> del{};
    del(pts);
    SortedMesh2 mesh(del);
    IntpLinearNonuniform2<double, SortedMesh2> intp(mesh, F.data());
    for (int32_t q = 0; q < 6; ++q)
    {
        Vector2<double> P = MakeQuery2(io, del, false);
        double f = 0.0;
        bool valid = intp(P, f);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
        }
    }
}

ORACLE_CASE("IntpLinearNonuniform2.deviation.getIndices")
{
    // #135: upstream ignores GetIndices' false and blends F[0] three times
    // with the barycentrics, returning true; the port returns invalid.
    auto pts = MakePoints2(io, 3, 9);
    auto F = MakeValues<2>(io, pts);
    Delaunay2<double> del{};
    del(pts);
    NoIndicesMesh2 mesh(del);
    IntpLinearNonuniform2<double, NoIndicesMesh2> intp(mesh, F.data());
    for (int32_t q = 0; q < 4; ++q)
    {
        Vector2<double> P = MakeQuery2(io, del, false);
        double f = 0.0;
        bool valid = intp(P, f);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
        }
    }
}

ORACLE_CASE("IntpLinearNonuniform3.evaluate")
{
    // The real Delaunay3Mesh<double>; queries strictly inside or outside.
    auto pts = MakePoints3(io, 4, 8);
    auto F = MakeValues<3>(io, pts);
    Delaunay3<double> del{};
    del(pts);
    Delaunay3Mesh<double> mesh(del);
    IntpLinearNonuniform3<double, Delaunay3Mesh<double>> intp(mesh, F.data());
    for (int32_t q = 0; q < 6; ++q)
    {
        Vector3<double> P = MakeQuery3(io, del, true);
        double f = 0.0;
        bool valid = intp(P, f);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
        }
    }
}

ORACLE_CASE("IntpLinearNonuniform3.evaluate.sortedMesh")
{
    auto pts = MakePoints3(io, 4, 8);
    auto F = MakeValues<3>(io, pts);
    Delaunay3<double> del{};
    del(pts);
    SortedMesh3 mesh(del);
    IntpLinearNonuniform3<double, SortedMesh3> intp(mesh, F.data());
    for (int32_t q = 0; q < 6; ++q)
    {
        Vector3<double> P = MakeQuery3(io, del, false);
        double f = 0.0;
        bool valid = intp(P, f);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
        }
    }
}

ORACLE_CASE("IntpLinearNonuniform3.deviation.getIndices")
{
    // #135, 3D half: upstream blends F[0] four times and returns true.
    auto pts = MakePoints3(io, 4, 8);
    auto F = MakeValues<3>(io, pts);
    Delaunay3<double> del{};
    del(pts);
    NoIndicesMesh3 mesh(del);
    IntpLinearNonuniform3<double, NoIndicesMesh3> intp(mesh, F.data());
    for (int32_t q = 0; q < 4; ++q)
    {
        Vector3<double> P = MakeQuery3(io, del, false);
        double f = 0.0;
        bool valid = intp(P, f);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
        }
    }
}

namespace
{
    // ---- IntpQuadraticNonuniform2 ------------------------------------------
    // Derivative samples: uniform, lattice with signed zeros, or the exact
    // gradient of the affine data mode of MakeValues (constant c[0], c[1]).
    void MakeGradients(oracle::Ctx& io, size_t n, std::vector<double>& FX,
        std::vector<double>& FY)
    {
        int32_t mode = io.index() % 3;
        double gx = 0.5 * static_cast<double>(io.rawInteger(-4, 4));
        double gy = 0.5 * static_cast<double>(io.rawInteger(-4, 4));
        FX.resize(n);
        FY.resize(n);
        for (size_t i = 0; i < n; ++i)
        {
            double fx, fy;
            if (mode == 0)
            {
                fx = io.raw(-3.0, 3.0);
                fy = io.raw(-3.0, 3.0);
            }
            else if (mode == 1)
            {
                int32_t k = io.rawInteger(-2, 3);
                fx = (k == 3 ? -0.0 : static_cast<double>(k));
                fy = static_cast<double>(io.rawInteger(-2, 2));
            }
            else
            {
                fx = gx;
                fy = gy;
            }
            FX[i] = io.given(fx);
            FY[i] = io.given(fy);
        }
    }

    template <typename Intp>
    void EmitQuadratic(oracle::Ctx& io, Intp const& intp, Vector2<double> const& P)
    {
        double f = 0.0, fx = 0.0, fy = 0.0;
        bool valid = intp(P, f, fx, fy);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
            io.outReal(fx);
            io.outReal(fy);
        }
    }

    double MakeSpatialDelta(oracle::Ctx& io)
    {
        int32_t k = io.rawInteger(0, 4);
        return io.given(k == 0 ? 1.0 : (k == 1 ? 0.5 : (k == 2 ? 0.0 : io.raw(0.1, 4.0))));
    }
}

ORACLE_CASE("IntpQuadraticNonuniform2.fromDerivatives")
{
    // The real Delaunay2Mesh<double>. With the derivatives given, the
    // preprocessing is per triangle and order independent; the queries are
    // strictly inside a triangle or outside the hull.
    auto pts = MakePoints2(io, 3, 9);
    auto F = MakeValues<2>(io, pts);
    std::vector<double> FX, FY;
    MakeGradients(io, pts.size(), FX, FY);
    Delaunay2<double> del{};
    del(pts);
    Delaunay2Mesh<double> mesh(del);
    IntpQuadraticNonuniform2<double, Delaunay2Mesh<double>> intp(mesh, F.data(), FX.data(), FY.data());
    for (int32_t q = 0; q < 5; ++q)
    {
        EmitQuadratic(io, intp, MakeQuery2(io, del, true));
    }
}

ORACLE_CASE("IntpQuadraticNonuniform2.fromDerivatives.sortedMesh")
{
    auto pts = MakePoints2(io, 3, 9);
    auto F = MakeValues<2>(io, pts);
    std::vector<double> FX, FY;
    MakeGradients(io, pts.size(), FX, FY);
    Delaunay2<double> del{};
    del(pts);
    SortedMesh2 mesh(del);
    IntpQuadraticNonuniform2<double, SortedMesh2> intp(mesh, F.data(), FX.data(), FY.data());
    for (int32_t q = 0; q < 5; ++q)
    {
        EmitQuadratic(io, intp, MakeQuery2(io, del, false));
    }
}

ORACLE_CASE("IntpQuadraticNonuniform2.fromSpatialDelta.sortedMesh")
{
    auto pts = MakePoints2(io, 3, 9);
    auto F = MakeValues<2>(io, pts);
    double spatialDelta = MakeSpatialDelta(io);
    Delaunay2<double> del{};
    del(pts);
    SortedMesh2 mesh(del);
    IntpQuadraticNonuniform2<double, SortedMesh2> intp(mesh, F.data(), spatialDelta);
    for (int32_t q = 0; q < 5; ++q)
    {
        EmitQuadratic(io, intp, MakeQuery2(io, del, false));
    }
}

// The spatialDelta constructor over the real Delaunay2Mesh<double> is not a
// case: EstimateDerivatives sums the triangle normals per vertex in the
// mesh's triangle order, and upstream's Delaunay2<T> numbering is not even
// reproducible from run to run (3 of 2000 records changed bits between runs
// of the same executable; the order comes from pointer-keyed containers in
// the incremental insertion). Measured against the port it differed by at
// most 9.1e-15 (scaled) on 450 of 11571 outputs. The sortedMesh case above
// is the exact comparison of the same code.

namespace
{
    // SortedMesh2 with one per-triangle accessor failing for every triangle
    // (upstream Delaunay2Mesh-style: false and a zero-filled output). Modes:
    // 0 GetVertices, 1 GetIndices, 2 GetAdjacencies, 3 GetBarycentrics.
    class FailingMesh2 : public SortedMesh2
    {
    public:
        FailingMesh2(Delaunay2<double> const& del, int32_t mode)
            :
            SortedMesh2(del),
            mMode(mode)
        {
        }

        using SortedMesh2::GetVertices;
        using SortedMesh2::GetIndices;

        bool GetVertices(size_t t, std::array<Vector2<double>, 3>& vertices) const
        {
            if (mMode == 0)
            {
                for (auto& v : vertices)
                {
                    v.MakeZero();
                }
                return false;
            }
            return SortedMesh2::GetVertices(t, vertices);
        }

        bool GetIndices(size_t t, std::array<int32_t, 3>& indices) const
        {
            if (mMode == 1)
            {
                return false;
            }
            return SortedMesh2::GetIndices(t, indices);
        }

        bool GetAdjacencies(size_t t, std::array<int32_t, 3>& adjacencies) const
        {
            if (mMode == 2)
            {
                return false;
            }
            return SortedMesh2::GetAdjacencies(t, adjacencies);
        }

        bool GetBarycentrics(size_t t, Vector2<double> const& P, std::array<double, 3>& bary) const
        {
            if (mMode == 3)
            {
                bary.fill(0.0);
                return false;
            }
            return SortedMesh2::GetBarycentrics(t, P, bary);
        }

    private:
        int32_t mMode;
    };
}

ORACLE_CASE("IntpQuadraticNonuniform2.deviation.meshFlags")
{
    // #337: upstream discards every mesh-accessor flag and computes with the
    // zero-filled outputs; the port skips the failed triangle and reports
    // the query invalid (or leaves the triangle's data zero).
    int32_t mode = io.integer(0, 3);
    auto pts = MakePoints2(io, 3, 9);
    auto F = MakeValues<2>(io, pts);
    std::vector<double> FX, FY;
    MakeGradients(io, pts.size(), FX, FY);
    Delaunay2<double> del{};
    del(pts);
    FailingMesh2 mesh(del, mode);
    IntpQuadraticNonuniform2<double, FailingMesh2> intp(mesh, F.data(), FX.data(), FY.data());
    for (int32_t q = 0; q < 4; ++q)
    {
        EmitQuadratic(io, intp, MakeQuery2(io, del, false));
    }
}

namespace
{
    // ListMesh2: an explicit triangle list with explicit adjacencies, the
    // minimal duck-typed mesh. Containment is a scan in triangle order for
    // the first triangle whose exact barycentrics are all nonnegative; a
    // degenerate triangle is never reported. Used to put a collinear
    // triangle next to two proper ones, which Delaunay never produces:
    //   A=V0, B=V1, D=V2 on segment AB, C=V3 above:
    //   t0 = <A,D,C>, t1 = <D,B,C>, t2 = <A,B,D> (collinear).
    // Edge <D,A> of t2 is shared with t0 and <B,D> with t1, so t2's
    // Inscribe center enters both proper triangles' cross-edge points.
    class ListMesh2
    {
    public:
        ListMesh2(std::vector<Vector2<double>> const& vertices, std::vector<int32_t> const& indices,
            std::vector<int32_t> const& adjacencies)
            :
            mVertices(vertices),
            mIndices(indices),
            mAdjacencies(adjacencies)
        {
        }

        size_t GetNumVertices() const { return mVertices.size(); }
        size_t GetNumTriangles() const { return mIndices.size() / 3; }
        Vector2<double> const* GetVertices() const { return mVertices.data(); }
        int32_t const* GetIndices() const { return mIndices.data(); }
        size_t GetInvalidIndex() const { return std::numeric_limits<size_t>::max(); }

        size_t GetContainingTriangle(Vector2<double> const& P) const
        {
            for (size_t t = 0; t < GetNumTriangles(); ++t)
            {
                std::array<double, 3> bary{};
                if (GetBarycentrics(t, P, bary) && bary[0] >= 0.0 && bary[1] >= 0.0 && bary[2] >= 0.0)
                {
                    return t;
                }
            }
            return GetInvalidIndex();
        }

        bool GetIndices(size_t t, std::array<int32_t, 3>& indices) const
        {
            for (size_t j = 0; j < 3; ++j)
            {
                indices[j] = mIndices[3 * t + j];
            }
            return true;
        }

        bool GetAdjacencies(size_t t, std::array<int32_t, 3>& adjacencies) const
        {
            for (size_t j = 0; j < 3; ++j)
            {
                adjacencies[j] = mAdjacencies[3 * t + j];
            }
            return true;
        }

        bool GetVertices(size_t t, std::array<Vector2<double>, 3>& vertices) const
        {
            for (size_t j = 0; j < 3; ++j)
            {
                vertices[j] = mVertices[mIndices[3 * t + j]];
            }
            return true;
        }

        bool GetBarycentrics(size_t t, Vector2<double> const& P, std::array<double, 3>& bary) const
        {
            using Rational = BSRational<UIntegerAP32>;
            std::array<Vector2<Rational>, 3> rtV{};
            for (size_t i = 0; i < 3; ++i)
            {
                auto const& V = mVertices[mIndices[3 * t + i]];
                rtV[i] = { V[0], V[1] };
            }
            Vector2<Rational> rtP{ P[0], P[1] };
            std::array<Rational, 3> rtBary{};
            if (ComputeBarycentrics(rtP, rtV[0], rtV[1], rtV[2], rtBary))
            {
                for (size_t i = 0; i < 3; ++i)
                {
                    bary[i] = static_cast<double>(rtBary[i]);
                }
                return true;
            }
            bary.fill(0.0);
            return false;
        }

    private:
        std::vector<Vector2<double>> mVertices;
        std::vector<int32_t> mIndices, mAdjacencies;
    };
}

ORACLE_CASE("IntpQuadraticNonuniform2.degenerateTriangle")
{
    // Inscribe's failure return is ignored (#337, preserved). For a
    // collinear triangle with a positive perimeter upstream's Inscribe has
    // already stored the perimeter-weighted center when it returns false.
    std::vector<Vector2<double>> V(4);
    V[0] = io.latticeVec<2>(-3, 0);
    V[1] = V[0] + io.latticeDir<2>(-3, 3);
    double k = static_cast<double>(io.integer(1, 3)) * 0.25;
    V[2] = V[0] + k * (V[1] - V[0]);
    Vector2<double> dir = V[1] - V[0];
    double h = io.real(0.5, 3.0);
    V[3] = io.givenVec(V[0] + 0.5 * dir + h * Vector2<double>{ -dir[1], dir[0] });
    std::vector<double> F(4);
    for (size_t i = 0; i < 4; ++i)
    {
        F[i] = io.real(-3.0, 3.0);
    }
    bool useSpatialDelta = io.boolean();
    std::vector<double> FX(4), FY(4);
    for (size_t i = 0; i < 4; ++i)
    {
        FX[i] = io.real(-2.0, 2.0);
        FY[i] = io.real(-2.0, 2.0);
    }
    std::vector<int32_t> indices{ 0, 2, 3, 2, 1, 3, 0, 1, 2 };
    std::vector<int32_t> adjacencies{ 2, 1, -1, 2, -1, 0, -1, 1, 0 };
    ListMesh2 mesh(V, indices, adjacencies);
    std::unique_ptr<IntpQuadraticNonuniform2<double, ListMesh2>> intp;
    if (useSpatialDelta)
    {
        intp = std::make_unique<IntpQuadraticNonuniform2<double, ListMesh2>>(mesh, F.data(), 1.0);
    }
    else
    {
        intp = std::make_unique<IntpQuadraticNonuniform2<double, ListMesh2>>(mesh, F.data(),
            FX.data(), FY.data());
    }
    for (int32_t q = 0; q < 4; ++q)
    {
        double a = static_cast<double>(io.rawInteger(1, 6)) / 8.0;
        double b = static_cast<double>(io.rawInteger(1, 7 - static_cast<int32_t>(a * 8.0))) / 8.0;
        size_t t = static_cast<size_t>(io.rawInteger(0, 1));
        Vector2<double> P = a * V[indices[3 * t]] + b * V[indices[3 * t + 1]]
            + (1.0 - a - b) * V[indices[3 * t + 2]];
        EmitQuadratic(io, *intp, io.givenVec(P));
    }
}

namespace
{
    // ---- thin-plate splines ------------------------------------------------
    // Sample coordinates for one axis, unrecorded. Modes (per record):
    // uniform, small lattice (repeated points make A singular when
    // smooth = 0), and one axis constant (#191: a zero coordinate range gives
    // NaN coordinates when transformToUnitSquare is set).
    std::vector<double> RawTPSCoordinates(oracle::Ctx& io, int32_t n, int32_t mode, bool flat)
    {
        std::vector<double> c(n);
        double constant = static_cast<double>(io.rawInteger(-2, 2));
        for (int32_t i = 0; i < n; ++i)
        {
            if (flat)
            {
                c[i] = constant;
            }
            else if (mode == 0)
            {
                c[i] = io.raw(-3.0, 3.0);
            }
            else
            {
                c[i] = static_cast<double>(io.rawInteger(-3, 3));
            }
        }
        return c;
    }

    std::vector<double> TPSCoordinates(oracle::Ctx& io, int32_t n, int32_t mode, bool flat)
    {
        std::vector<double> c = RawTPSCoordinates(io, n, mode, flat);
        for (double& v : c)
        {
            v = io.given(v);
        }
        return c;
    }

    double RawTPSSmooth(oracle::Ctx& io)
    {
        int32_t k = io.rawInteger(0, 3);
        return (k <= 1 ? 0.0 : (k == 2 ? 0.125 : io.raw(0.001, 2.0)));
    }

    double TPSSmooth(oracle::Ctx& io)
    {
        return io.given(RawTPSSmooth(io));
    }

    double Norm1(GMatrix<double> const& M)
    {
        double best = 0.0;
        for (int32_t c = 0; c < M.GetNumCols(); ++c)
        {
            double sum = 0.0;
            for (int32_t r = 0; r < M.GetNumRows(); ++r)
            {
                sum += std::fabs(M(r, c));
            }
            best = std::max(best, sum);
        }
        return best;
    }

    // Conditioning probe for IntpThinPlateSpline2: cond1(A) * cond1(Q) of
    // the two matrices the constructor inverts, built exactly as upstream
    // builds them (std::minmax_element, t^2 log t^2, GMatrix Inverse,
    // MultiplyATB). Infinity when either inverse fails or is not finite.
    double TPS2Conditioning(std::vector<double> const& X, std::vector<double> const& Y,
        double smooth, bool transform)
    {
        int32_t n = static_cast<int32_t>(X.size());
        std::vector<double> x(X), y(Y);
        if (transform)
        {
            auto ex = std::minmax_element(X.begin(), X.end());
            auto ey = std::minmax_element(Y.begin(), Y.end());
            double xinv = 1.0 / (*ex.second - *ex.first);
            double yinv = 1.0 / (*ey.second - *ey.first);
            for (int32_t i = 0; i < n; ++i)
            {
                x[i] = (X[i] - *ex.first) * xinv;
                y[i] = (Y[i] - *ey.first) * yinv;
            }
        }
        GMatrix<double> A(n, n), B(n, 3);
        for (int32_t r = 0; r < n; ++r)
        {
            for (int32_t c = 0; c < n; ++c)
            {
                double dx = x[r] - x[c], dy = y[r] - y[c];
                double t2 = dx * dx + dy * dy;
                A(r, c) = (r == c ? smooth : (t2 > 0.0 ? t2 * std::log(t2) : 0.0));
            }
            B(r, 0) = 1.0;
            B(r, 1) = x[r];
            B(r, 2) = y[r];
        }
        double const inf = std::numeric_limits<double>::infinity();
        bool invertible = false;
        GMatrix<double> invA = Inverse(A, &invertible);
        if (!invertible)
        {
            return inf;
        }
        GMatrix<double> Q = MultiplyATB(B, invA) * B;
        GMatrix<double> invQ = Inverse(Q, &invertible);
        if (!invertible)
        {
            return inf;
        }
        double cond = Norm1(A) * Norm1(invA) * Norm1(Q) * Norm1(invQ);
        return std::isfinite(cond) ? cond : inf;
    }
}

// std::log in the kernel, compared with a tolerance. MSVC's log and V8's
// Math.log differ in the last bit on about 3.6% of the kernel arguments
// (5594 of 153924 in a 2000-record probe), and the two inverted systems
// amplify that by their conditioning. With MSVC's log values substituted
// for Math.log the port reproduces all 2000 probe records bit for bit,
// invertibility flags included (see the report). The generator therefore
// accepts only records with cond1(A) * cond1(Q) <= 1e4 (capped rejection,
// fallback the best-conditioned candidate), where the measured scaled error
// stays below 1e-12; the replay compares at 1e-11. The flat-axis mode is
// always transformed and takes the #191 NaN path, which is exact.
ORACLE_CASE("IntpThinPlateSpline2.evaluate")
{
    int32_t n = io.integer(3, 8);
    int32_t mode = io.index() % 2;
    bool flatAxis = (io.index() % 7 == 6);
    std::vector<double> X, Y;
    double smooth = 0.0, best = std::numeric_limits<double>::infinity();
    bool transform = false;
    for (int32_t attempt = 0; attempt < 64; ++attempt)
    {
        auto x = RawTPSCoordinates(io, n, mode, false);
        auto y = RawTPSCoordinates(io, n, mode, flatAxis);
        double s = RawTPSSmooth(io);
        bool tr = (flatAxis || io.rawInteger(0, 1) == 1);
        double cond = (flatAxis ? 0.0 : TPS2Conditioning(x, y, s, tr));
        if (attempt == 0 || cond < best)
        {
            best = cond;
            X = x;
            Y = y;
            smooth = s;
            transform = tr;
        }
        if (best <= 1e4)
        {
            break;
        }
    }
    for (double& v : X)
    {
        v = io.given(v);
    }
    for (double& v : Y)
    {
        v = io.given(v);
    }
    std::vector<double> F(n);
    for (int32_t i = 0; i < n; ++i)
    {
        F[i] = io.real(-5.0, 5.0);
    }
    smooth = io.given(smooth);
    io.given(transform ? 1.0 : 0.0);
    IntpThinPlateSpline2<double> tps(n, X.data(), Y.data(), F.data(), smooth, transform);
    io.outBool(tps.IsInitialized());
    for (int32_t i = 0; i < n; ++i)
    {
        io.outReal(tps(X[i], Y[i]));
    }
    for (int32_t q = 0; q < 3; ++q)
    {
        double x = io.real(-4.0, 4.0);
        double y = io.real(-4.0, 4.0);
        io.outReal(tps(x, y));
    }
    io.outReal(tps.ComputeFunctional());
}

// Kernel(t) = -|t|: + - * / sqrt only, exact.
ORACLE_CASE("IntpThinPlateSpline3.evaluate")
{
    int32_t n = io.integer(4, 8);
    int32_t mode = io.index() % 2;
    bool flatAxis = (io.index() % 7 == 6);
    auto X = TPSCoordinates(io, n, mode, false);
    auto Y = TPSCoordinates(io, n, mode, false);
    auto Z = TPSCoordinates(io, n, mode, flatAxis);
    std::vector<double> F(n);
    for (int32_t i = 0; i < n; ++i)
    {
        F[i] = io.real(-5.0, 5.0);
    }
    double smooth = TPSSmooth(io);
    bool transform = io.boolean();
    IntpThinPlateSpline3<double> tps(n, X.data(), Y.data(), Z.data(), F.data(), smooth, transform);
    io.outBool(tps.IsInitialized());
    for (int32_t i = 0; i < n; ++i)
    {
        io.outReal(tps(X[i], Y[i], Z[i]));
    }
    for (int32_t q = 0; q < 3; ++q)
    {
        double x = io.real(-4.0, 4.0);
        double y = io.real(-4.0, 4.0);
        double z = io.real(-4.0, 4.0);
        io.outReal(tps(x, y, z));
    }
    io.outReal(tps.ComputeFunctional());
}

ORACLE_CASE("IntpThinPlateSpline2.construct.invalid")
{
    int32_t n = io.integer(1, 4);
    std::vector<double> X(n), Y(n), F(n);
    for (int32_t i = 0; i < n; ++i)
    {
        X[i] = io.lattice(-3, 3);
        Y[i] = io.lattice(-3, 3);
        F[i] = io.lattice(-3, 3);
    }
    double smooth = io.given(io.rawInteger(0, 1) == 0 ? -0.5 : 0.0);
    IntpThinPlateSpline2<double> tps(n, X.data(), Y.data(), F.data(), smooth, false);
    io.outBool(tps.IsInitialized());
}

ORACLE_CASE("IntpThinPlateSpline3.construct.invalid")
{
    int32_t n = io.integer(1, 5);
    std::vector<double> X(n), Y(n), Z(n), F(n);
    for (int32_t i = 0; i < n; ++i)
    {
        X[i] = io.lattice(-3, 3);
        Y[i] = io.lattice(-3, 3);
        Z[i] = io.lattice(-3, 3);
        F[i] = io.lattice(-3, 3);
    }
    double smooth = io.given(io.rawInteger(0, 1) == 0 ? -0.5 : 0.0);
    IntpThinPlateSpline3<double> tps(n, X.data(), Y.data(), Z.data(), F.data(), smooth, false);
    io.outBool(tps.IsInitialized());
}

namespace
{
    // ---- IntpSphere2 --------------------------------------------------------
    // The wrap-around point set exactly as IntpSphere2<T>'s constructor
    // builds it: the samples, then theta + 2*pi, then theta - 2*pi.
    std::vector<Vector2<double>> WrapAngles(std::vector<double> const& theta,
        std::vector<double> const& phi)
    {
        size_t n = theta.size();
        std::vector<Vector2<double>> wrap(3 * n);
        for (size_t i = 0; i < n; ++i)
        {
            wrap[i] = { theta[i], phi[i] };
        }
        for (size_t i0 = 0, i1 = n, i2 = 2 * n; i0 < n; ++i0, ++i1, ++i2)
        {
            wrap[i1][0] = wrap[i0][0] + static_cast<double>(GTE_C_TWO_PI);
            wrap[i2][0] = wrap[i0][0] - static_cast<double>(GTE_C_TWO_PI);
            wrap[i1][1] = wrap[i0][1];
            wrap[i2][1] = wrap[i0][1];
        }
        return wrap;
    }

    // Samples (theta, phi, F). Modes: uniform angles; angles on a pi/4
    // lattice (repeated phi values: horizontal collinear runs, and theta =
    // -pi and pi, whose wrapped copies coincide); uniform plus the two
    // poles (-pi, 0) and (-pi, pi) the header recommends. Rejection (capped,
    // then a fixed ring) on Sound2 of the wrapped set, so upstream's
    // Delaunay2<T> seed is the port's (#391).
    void MakeSphereSamples(oracle::Ctx& io, std::vector<double>& theta, std::vector<double>& phi,
        std::vector<double>& F)
    {
        int32_t n = io.integer(4, 7);
        int32_t mode = io.index() % 3;
        double const pi = static_cast<double>(GTE_C_PI);
        theta.resize(n);
        phi.resize(n);
        bool sound = false;
        for (int32_t attempt = 0; attempt < 64 && !sound; ++attempt)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                if (mode == 1)
                {
                    theta[i] = 0.25 * pi * static_cast<double>(io.rawInteger(-4, 4));
                    phi[i] = 0.25 * pi * static_cast<double>(io.rawInteger(0, 4));
                }
                else
                {
                    theta[i] = io.raw(-pi, pi);
                    phi[i] = io.raw(0.0, pi);
                }
            }
            if (mode == 2)
            {
                theta[0] = -pi;
                phi[0] = 0.0;
                theta[1] = -pi;
                phi[1] = pi;
            }
            sound = Sound2(WrapAngles(theta, phi));
        }
        if (!sound)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                theta[i] = -pi + 2.0 * pi * static_cast<double>(i) / static_cast<double>(n);
                phi[i] = 0.25 * pi + 0.5 * pi * static_cast<double>(i % 2);
            }
        }
        for (int32_t i = 0; i < n; ++i)
        {
            io.given(theta[i]);
        }
        for (int32_t i = 0; i < n; ++i)
        {
            io.given(phi[i]);
        }
        F.resize(n);
        for (int32_t i = 0; i < n; ++i)
        {
            F[i] = io.real(-3.0, 3.0);
        }
    }

    // Query angles: uniform over the domain, a sample point, a dyadic mix
    // of two samples, outside the phi range of the samples.
    Vector2<double> SphereQuery(oracle::Ctx& io, std::vector<double> const& theta,
        std::vector<double> const& phi)
    {
        double const pi = static_cast<double>(GTE_C_PI);
        int32_t n = static_cast<int32_t>(theta.size());
        int32_t mode = io.rawInteger(0, 3);
        Vector2<double> q;
        if (mode == 0)
        {
            q = { io.raw(-pi, pi), io.raw(0.0, pi) };
        }
        else if (mode == 1)
        {
            int32_t i = io.rawInteger(0, n - 1);
            q = { theta[i], phi[i] };
        }
        else if (mode == 2)
        {
            int32_t i = io.rawInteger(0, n - 1);
            int32_t j = io.rawInteger(0, n - 1);
            q = { 0.5 * (theta[i] + theta[j]), 0.5 * (phi[i] + phi[j]) };
        }
        else
        {
            q = { io.raw(-pi, pi), (io.rawInteger(0, 1) == 0 ? -0.5 : pi + 0.5) };
        }
        return io.givenVec(q);
    }
}

// std::atan2 and std::acos: the default 1e-12 tolerance. The two pole
// branches return constants and are hit exactly by the z = +-1 modes.
ORACLE_CASE("IntpSphere2.getSphericalCoordinates")
{
    int32_t mode = io.index() % 5;
    Vector3<double> v;
    if (mode == 0 || mode == 1)
    {
        v = io.unit<3>();
    }
    else if (mode == 2)
    {
        double s = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        v = io.givenVec(Vector3<double>{ 0.0, 0.0, s });
    }
    else if (mode == 3)
    {
        // Just inside the poles.
        double z = 1.0 - std::ldexp(1.0, -io.rawInteger(20, 52));
        z = (io.rawInteger(0, 1) == 0 ? z : -z);
        double r = std::sqrt(1.0 - z * z);
        v = io.givenVec(Vector3<double>{ r, 0.0, z });
    }
    else
    {
        // Not unit length: |z| >= 1 takes the pole branches.
        v = io.givenVec(Vector3<double>{ io.raw(-1.0, 1.0), io.raw(-1.0, 1.0),
            (io.rawInteger(0, 1) == 0 ? 1.5 : -2.0) });
    }
    double theta = 0.0, phi = 0.0;
    IntpSphere2<double>::GetSphericalCoordinates(v[0], v[1], v[2], theta, phi);
    io.outReal(theta);
    io.outReal(phi);
}

// The real IntpSphere2<double> cannot be constructed: its member
// initializer list constructs mMesh(mDelaunay), and Delaunay2Mesh<T>'s
// constructor asserts mDelaunay->GetDimension() == 2, before the constructor
// body has run mDelaunay(mWrapAngles). A default-constructed Delaunay2<T>
// has dimension 0, so every construction throws "Invalid Delaunay
// dimension." (new finding, see the group report). The port builds the mesh
// after the triangulation. All inputs are drawn before the constructor so
// that the record is complete; every record is a C++ throw.
ORACLE_CASE("IntpSphere2.deviation.constructorThrows")
{
    std::vector<double> theta, phi, F;
    MakeSphereSamples(io, theta, phi, F);
    std::vector<Vector2<double>> queries(5);
    for (auto& a : queries)
    {
        a = SphereQuery(io, theta, phi);
    }
    IntpSphere2<double> intp(static_cast<int32_t>(theta.size()), theta.data(), phi.data(), F.data());
    for (auto const& a : queries)
    {
        double f = 0.0;
        bool valid = intp(a[0], a[1], f);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
        }
    }
}

// IntpSphere2<T>'s constructor and operator() replayed over SortedMesh2:
// the port's IntpSphere2 is compared bit for bit.
ORACLE_CASE("IntpSphere2.evaluate.sortedMesh")
{
    std::vector<double> theta, phi, F;
    MakeSphereSamples(io, theta, phi, F);
    auto wrap = WrapAngles(theta, phi);
    size_t n = theta.size();
    std::vector<double> wrapF(3 * n);
    for (size_t i = 0; i < n; ++i)
    {
        wrapF[i] = F[i];
        wrapF[i + n] = F[i];
        wrapF[i + 2 * n] = F[i];
    }
    Delaunay2<double> del{};
    del(wrap);
    SortedMesh2 mesh(del);
    IntpQuadraticNonuniform2<double, SortedMesh2> intp(mesh, wrapF.data(), 1.0);
    for (int32_t q = 0; q < 5; ++q)
    {
        Vector2<double> a = SphereQuery(io, theta, phi);
        double f = 0.0, thetaDeriv = 0.0, phiDeriv = 0.0;
        bool valid = intp(a, f, thetaDeriv, phiDeriv);
        io.outBool(valid);
        if (valid)
        {
            io.outReal(f);
        }
    }
}

namespace
{
    // ---- IntpVectorField2 ---------------------------------------------------
    std::vector<Vector2<double>> MakeRange(oracle::Ctx& io, size_t n)
    {
        int32_t mode = io.index() % 3;
        std::vector<Vector2<double>> range(n);
        for (size_t i = 0; i < n; ++i)
        {
            Vector2<double> r;
            if (mode == 0)
            {
                r = { io.raw(-5.0, 5.0), io.raw(-5.0, 5.0) };
            }
            else if (mode == 1)
            {
                r = { static_cast<double>(io.rawInteger(-3, 3)), static_cast<double>(io.rawInteger(-3, 3)) };
            }
            else
            {
                r = { 0.0, (io.rawInteger(0, 1) == 0 ? -0.0 : 1.0) };
            }
            range[i] = io.givenVec(r);
        }
        return range;
    }
}

// The real IntpVectorField2<double> cannot be constructed either: the same
// mMesh(mDelaunay) initializer runs Delaunay2Mesh<T>'s dimension assert
// before mDelaunay(numPoints, domain). Every record is a C++ throw.
ORACLE_CASE("IntpVectorField2.deviation.constructorThrows")
{
    auto domain = MakePoints2(io, 3, 9);
    auto range = MakeRange(io, domain.size());
    Delaunay2<double> del{};
    del(domain);
    std::vector<Vector2<double>> queries(5);
    for (auto& P : queries)
    {
        P = MakeQuery2(io, del, false);
    }
    IntpVectorField2<double> intp(static_cast<int32_t>(domain.size()), domain.data(), range.data());
    for (auto const& P : queries)
    {
        Vector2<double> output{ 0.0, 0.0 };
        bool valid = intp(P, output);
        io.outBool(valid);
        if (valid)
        {
            io.outVec(output);
        }
    }
}

// IntpVectorField2<T>'s constructor and operator() replayed over
// SortedMesh2; compared exactly with the port's class, every query mode.
ORACLE_CASE("IntpVectorField2.evaluate.sortedMesh")
{
    auto domain = MakePoints2(io, 3, 9);
    auto range = MakeRange(io, domain.size());
    std::vector<double> xRange(domain.size()), yRange(domain.size());
    for (size_t i = 0; i < domain.size(); ++i)
    {
        xRange[i] = range[i][0];
        yRange[i] = range[i][1];
    }
    Delaunay2<double> del{};
    del(domain);
    SortedMesh2 mesh(del);
    IntpQuadraticNonuniform2<double, SortedMesh2> xInterp(mesh, xRange.data(), 1.0);
    IntpQuadraticNonuniform2<double, SortedMesh2> yInterp(mesh, yRange.data(), 1.0);
    for (int32_t q = 0; q < 5; ++q)
    {
        Vector2<double> P = MakeQuery2(io, del, false);
        Vector2<double> output{ 0.0, 0.0 };
        double xDeriv = 0.0, yDeriv = 0.0;
        bool valid = xInterp(P, output[0], xDeriv, yDeriv)
            && yInterp(P, output[1], xDeriv, yDeriv);
        io.outBool(valid);
        if (valid)
        {
            io.outVec(output);
        }
    }
}

// #337: on a failed query upstream's operator() leaves the caller's output
// untouched (the x-interpolation fails before writing F, and the && skips
// the y-interpolation), where the port returns (0, 0). Since the real class
// cannot be constructed, the demonstration runs upstream's operator()
// expression over the sortedMesh replay; the output is emitted even when the
// query fails, from a caller-initialized sentinel.
ORACLE_CASE("IntpVectorField2.deviation.staleOutput")
{
    auto domain = MakePoints2(io, 3, 9);
    auto range = MakeRange(io, domain.size());
    std::vector<double> xRange(domain.size()), yRange(domain.size());
    for (size_t i = 0; i < domain.size(); ++i)
    {
        xRange[i] = range[i][0];
        yRange[i] = range[i][1];
    }
    Delaunay2<double> del{};
    del(domain);
    SortedMesh2 mesh(del);
    IntpQuadraticNonuniform2<double, SortedMesh2> xInterp(mesh, xRange.data(), 1.0);
    IntpQuadraticNonuniform2<double, SortedMesh2> yInterp(mesh, yRange.data(), 1.0);
    for (int32_t q = 0; q < 3; ++q)
    {
        double s = (io.rawInteger(0, 1) == 0 ? -1.0 : 1.0);
        Vector2<double> P = io.givenVec(Vector2<double>{ s * io.raw(8.0, 10.0), io.raw(-10.0, 10.0) });
        Vector2<double> output{ 123.0, -456.0 };
        double xDeriv = 0.0, yDeriv = 0.0;
        bool valid = xInterp(P, output[0], xDeriv, yDeriv)
            && yInterp(P, output[1], xDeriv, yDeriv);
        io.outBool(valid);
        io.outVec(output);
    }
}
