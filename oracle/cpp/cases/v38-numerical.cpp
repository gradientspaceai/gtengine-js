// Group 38 (numerical): eigensolvers, banded and dense linear solvers,
// numerical integration, the LCP solver and the rounding-mode interval type.
//
// Nearly every path here is + - * / sqrt on doubles, which MSVC (/fp:precise,
// x64 SSE2, no FMA) and V8 both evaluate as IEEE 754 requires, so the
// expectation is bit-identity including the eigensolver iteration counts and
// the LCP pivoting sequence. The two exceptions are stated at their cases:
//   * NISymmetricEigensolver3x3 calls std::acos and std::cos;
//   * one Integration case deliberately uses a libm integrand.
//
// FENV_ACCESS. FPInterval.h switches the hardware rounding mode with
// std::fesetround around every arithmetic operation. Measured on this
// toolchain (MSVC 19.44, x64, /O2): with plain /fp:precise the compiler
// common-subexpression-eliminates the two evaluations of "u / v" that
// straddle the fesetround calls, so both endpoints come back as the
// round-to-nearest value and upstream's intervals are degenerate. With
// "#pragma fenv_access (on)" -- which /fp:precise accepts -- the two
// evaluations are kept and the rounding mode takes effect (1.0/3.0 yields
// 3fd5555555555555 and 3fd5555555555556). The pragma is therefore enabled
// for this translation unit; it was verified that the goldens of every other
// case in this file are byte-identical with and without it.
#pragma fenv_access (on)

#define ORACLE_FAMILY "v38-numerical"
#include "Oracle.h"

#include <Mathematics/BandedMatrix.h>
#include <Mathematics/FPInterval.h>
#include <Mathematics/GaussianElimination.h>
#include <Mathematics/Integration.h>
#include <Mathematics/LCPSolver.h>
#include <Mathematics/SymmetricEigensolver2x2.h>
#include <Mathematics/SymmetricEigensolver3x3.h>
#include <Mathematics/UnsymmetricEigenvalues.h>

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <vector>

using namespace gte;

namespace
{
    // ---------------------------------------------------------------------
    // Shared generators
    // ---------------------------------------------------------------------

    // Exact powers of two, so that a lattice matrix scaled by one of them
    // stays exactly representable and only the exponent range changes.
    double PowerOfTwo(int32_t exponent)
    {
        return std::ldexp(1.0, exponent);
    }

    // A symmetric 3x3 matrix, recorded as the six unique entries a00, a01,
    // a02, a11, a12, a22 (one io draw per statement). The mode is recorded
    // first and every mode records exactly six doubles afterwards.
    //
    // Modes: 0 uniform, 1 small lattice, 2 diagonal, 3 already tridiagonal
    // (a02 = 0), 4 a multiple of the identity (a triple eigenvalue),
    // 5 rank one v*v^T on a lattice (a double eigenvalue 0), 6 the zero
    // matrix, 7 a lattice matrix scaled by 2^k with |k| <= 500 (the in-band
    // extreme scales), 8 signed zeros mixed with lattice entries.
    struct Sym3
    {
        int32_t mode = 0;
        double a00 = 0.0, a01 = 0.0, a02 = 0.0, a11 = 0.0, a12 = 0.0, a22 = 0.0;
    };

    // The unrecorded draw, so that a rejection loop can redraw every
    // quantity its acceptance test depends on. RecordSym3 writes the
    // accepted draw: the mode, then the six entries.
    Sym3 DrawSym3Raw(oracle::Ctx& io, int32_t numModes = 9, int32_t forcedMode = -1)
    {
        int32_t mode = io.rawInteger(0, numModes - 1);
        if (forcedMode >= 0) { mode = forcedMode; }
        Sym3 a{};
        a.mode = mode;
        if (mode == 0)
        {
            a.a00 = io.raw(-5.0, 5.0);
            a.a01 = io.raw(-5.0, 5.0);
            a.a02 = io.raw(-5.0, 5.0);
            a.a11 = io.raw(-5.0, 5.0);
            a.a12 = io.raw(-5.0, 5.0);
            a.a22 = io.raw(-5.0, 5.0);
        }
        else if (mode == 1)
        {
            a.a00 = io.rawInteger(-4, 4);
            a.a01 = io.rawInteger(-4, 4);
            a.a02 = io.rawInteger(-4, 4);
            a.a11 = io.rawInteger(-4, 4);
            a.a12 = io.rawInteger(-4, 4);
            a.a22 = io.rawInteger(-4, 4);
        }
        else if (mode == 2)
        {
            a.a00 = io.rawInteger(-4, 4);
            a.a11 = io.rawInteger(-4, 4);
            a.a22 = io.rawInteger(-4, 4);
        }
        else if (mode == 3)
        {
            a.a00 = io.rawInteger(-4, 4);
            a.a01 = io.rawInteger(-4, 4);
            a.a11 = io.rawInteger(-4, 4);
            a.a12 = io.rawInteger(-4, 4);
            a.a22 = io.rawInteger(-4, 4);
        }
        else if (mode == 4)
        {
            double d = io.rawInteger(-4, 4);
            a.a00 = d;
            a.a11 = d;
            a.a22 = d;
        }
        else if (mode == 5)
        {
            double v0 = io.rawInteger(-3, 3);
            double v1 = io.rawInteger(-3, 3);
            double v2 = io.rawInteger(-3, 3);
            a.a00 = v0 * v0;
            a.a01 = v0 * v1;
            a.a02 = v0 * v2;
            a.a11 = v1 * v1;
            a.a12 = v1 * v2;
            a.a22 = v2 * v2;
        }
        else if (mode == 6)
        {
            // The zero matrix; every entry stays at its zero initializer.
        }
        else if (mode == 7)
        {
            double scale = PowerOfTwo(io.rawInteger(-500, 500));
            a.a00 = scale * io.rawInteger(-4, 4);
            a.a01 = scale * io.rawInteger(-4, 4);
            a.a02 = scale * io.rawInteger(-4, 4);
            a.a11 = scale * io.rawInteger(-4, 4);
            a.a12 = scale * io.rawInteger(-4, 4);
            a.a22 = scale * io.rawInteger(-4, 4);
        }
        else
        {
            // Signed zeros next to lattice entries: the sign of a zero
            // survives multiplication and decides std::max/std::min ties.
            double signedZero[2] = { 0.0, -0.0 };
            a.a00 = signedZero[io.rawInteger(0, 1)];
            a.a01 = io.rawInteger(-2, 2);
            a.a02 = signedZero[io.rawInteger(0, 1)];
            a.a11 = signedZero[io.rawInteger(0, 1)];
            a.a12 = io.rawInteger(-2, 2);
            a.a22 = signedZero[io.rawInteger(0, 1)];
        }

        return a;
    }

    void RecordSym3(oracle::Ctx& io, Sym3 const& a)
    {
        io.given(static_cast<double>(a.mode));
        io.given(a.a00);
        io.given(a.a01);
        io.given(a.a02);
        io.given(a.a11);
        io.given(a.a12);
        io.given(a.a22);
    }

    Sym3 DrawSym3(oracle::Ctx& io, int32_t numModes = 9, int32_t forcedMode = -1)
    {
        Sym3 a = DrawSym3Raw(io, numModes, forcedMode);
        RecordSym3(io, a);
        return a;
    }

    // -----------------------------------------------------------------
    // Probe for gtengine-js issue #42 (UnsymmetricEigenvalues drops the
    // trailing eigenvalue). Upstream's eigenvalue-packing loop runs over
    // i < N-1 and therefore never reports A(N-1,N-1) when the last
    // subdiagonal entry has decoupled, that is when the final
    // mSubdiagonalFlag[N-2] is 0. The port adds the missing 1x1 block, so
    // the two sides report a different NUMBER of eigenvalues on exactly
    // those inputs.
    //
    // mSubdiagonalFlag is private, so this is a verbatim copy of the parts
    // of UnsymmetricEigenvalues<double> that produce it: House, RowHouse,
    // ColHouse, ReduceToUpperHessenberg, FrancisQRStep, GetBlock and the
    // iteration loop of Solve, with the packing and sorting removed. Run on
    // the same input it reproduces upstream's control flow exactly, and
    // TrailingFlag() is the separator: 0 means upstream drops an
    // eigenvalue, 1 means it does not.
    class ProbeUnsymmetric
    {
    public:
        ProbeUnsymmetric(int32_t size, uint32_t maxIterations)
            :
            mSize(size), mSizeM1(size - 1), mMaxIterations(maxIterations),
            mMatrix(static_cast<size_t>(size)* static_cast<size_t>(size), 0.0),
            mX(size, 0.0), mV(size, 0.0), mScaledV(size, 0.0), mW(size, 0.0),
            mFlagStorage(static_cast<size_t>(size) + 1, 0),
            mSubdiagonalFlag(&mFlagStorage[1])
        {
        }

        void Solve(double const* input)
        {
            std::copy(input, input + static_cast<size_t>(mSize) * mSize, mMatrix.begin());
            ReduceToUpperHessenberg();
            std::array<int32_t, 2> block{};
            bool found = GetBlock(block);
            for (uint32_t numIterations = 0; numIterations < mMaxIterations; ++numIterations)
            {
                if (found)
                {
                    FrancisQRStep(block[0], block[1] + 1);
                    found = GetBlock(block);
                }
                else
                {
                    break;
                }
            }
        }

        int32_t TrailingFlag() const
        {
            return mSubdiagonalFlag[mSizeM1 - 1];
        }

    private:
        double const& A(int32_t r, int32_t c) const
        {
            return mMatrix[c + r * static_cast<size_t>(mSize)];
        }

        double& A(int32_t r, int32_t c)
        {
            return mMatrix[c + r * static_cast<size_t>(mSize)];
        }

        void House(int32_t rmin, int32_t rmax)
        {
            double length = 0.0;
            for (int32_t r = rmin; r <= rmax; ++r) { length += mX[r] * mX[r]; }
            length = std::sqrt(length);
            if (length != 0.0)
            {
                double sign = (mX[rmin] >= 0.0 ? 1.0 : -1.0);
                double invDenom = 1.0 / (mX[rmin] + sign * length);
                for (int32_t r = rmin + 1; r <= rmax; ++r) { mV[r] = mX[r] * invDenom; }
            }
            mV[rmin] = 1.0;

            double dot = 1.0;
            for (int32_t r = rmin + 1; r <= rmax; ++r) { dot += mV[r] * mV[r]; }
            double scale = -2.0 / dot;
            for (int32_t r = rmin; r <= rmax; ++r) { mScaledV[r] = scale * mV[r]; }
        }

        void RowHouse(int32_t rmin, int32_t rmax, int32_t cmin, int32_t cmax)
        {
            for (int32_t c = cmin; c <= cmax; ++c)
            {
                mW[c] = 0.0;
                for (int32_t r = rmin; r <= rmax; ++r) { mW[c] += mScaledV[r] * A(r, c); }
            }
            for (int32_t r = rmin; r <= rmax; ++r)
            {
                for (int32_t c = cmin; c <= cmax; ++c) { A(r, c) += mV[r] * mW[c]; }
            }
        }

        void ColHouse(int32_t rmin, int32_t rmax, int32_t cmin, int32_t cmax)
        {
            for (int32_t r = rmin; r <= rmax; ++r)
            {
                mW[r] = 0.0;
                for (int32_t c = cmin; c <= cmax; ++c) { mW[r] += mScaledV[c] * A(r, c); }
            }
            for (int32_t r = rmin; r <= rmax; ++r)
            {
                for (int32_t c = cmin; c <= cmax; ++c) { A(r, c) += mW[r] * mV[c]; }
            }
        }

        void ReduceToUpperHessenberg()
        {
            for (int32_t c = 0, cp1 = 1; c <= mSize - 3; ++c, ++cp1)
            {
                for (int32_t r = cp1; r <= mSizeM1; ++r) { mX[r] = A(r, c); }
                House(cp1, mSizeM1);
                RowHouse(cp1, mSizeM1, c, mSizeM1);
                ColHouse(0, mSizeM1, cp1, mSizeM1);
            }
        }

        void FrancisQRStep(int32_t rmin, int32_t rmax)
        {
            int32_t const i0 = rmax - 1, i1 = rmax;
            double a00 = A(i0, i0);
            double a01 = A(i0, i1);
            double a10 = A(i1, i0);
            double a11 = A(i1, i1);
            double tr = a00 + a11;
            double det = a00 * a11 - a01 * a10;

            int32_t const j0 = rmin, j1 = j0 + 1, j2 = j1 + 1;
            double b00 = A(j0, j0);
            double b01 = A(j0, j1);
            double b10 = A(j1, j0);
            double b11 = A(j1, j1);
            double b21 = A(j2, j1);
            mX[rmin] = b00 * (b00 - tr) + b01 * b10 + det;
            mX[static_cast<size_t>(rmin) + 1] = b10 * (b00 + b11 - tr);
            mX[static_cast<size_t>(rmin) + 2] = b10 * b21;

            House(rmin, rmin + 2);
            RowHouse(rmin, rmin + 2, rmin, rmax);
            ColHouse(rmin, std::min(rmax, rmin + 3), rmin, rmin + 2);

            for (int32_t c = 0, cp1 = 1; c <= mSize - 3; ++c, ++cp1)
            {
                int32_t kmax = std::min(cp1 + 2, mSizeM1);
                for (int32_t r = cp1; r <= kmax; ++r) { mX[r] = A(r, c); }
                House(cp1, kmax);
                RowHouse(cp1, kmax, c, mSizeM1);
                ColHouse(0, mSizeM1, cp1, kmax);
            }
        }

        bool GetBlock(std::array<int32_t, 2>& block)
        {
            for (int32_t i = 0; i < mSizeM1; ++i)
            {
                double a00 = A(i, i);
                double a11 = A(i + 1, i + 1);
                double a21 = A(i + 1, i);
                double sum0 = a00 + a11;
                double sum1 = sum0 + a21;
                mSubdiagonalFlag[i] = (sum1 != sum0 ? 1 : 0);
            }

            for (int32_t i = 0; i < mSizeM1; ++i)
            {
                if (mSubdiagonalFlag[i] == 1)
                {
                    block = { i, -1 };
                    while (i < mSizeM1 && mSubdiagonalFlag[i] == 1) { block[1] = i++; }
                    if (block[1] != block[0]) { return true; }
                }
            }
            return false;
        }

        int32_t mSize, mSizeM1;
        uint32_t mMaxIterations;
        std::vector<double> mMatrix, mX, mV, mScaledV, mW;
        std::vector<int32_t> mFlagStorage;
        int32_t* mSubdiagonalFlag;
    };

    void OutSym3Result(oracle::Ctx& io, std::array<double, 3> const& eval,
        std::array<std::array<double, 3>, 3> const& evec)
    {
        for (size_t i = 0; i < 3; ++i) { io.outReal(eval[i]); }
        for (size_t i = 0; i < 3; ++i)
        {
            for (size_t j = 0; j < 3; ++j) { io.outReal(evec[i][j]); }
        }
    }

    // A general NxN matrix in row-major order, drawn without recording.
    // Modes: 0 uniform, 1 small lattice, 2 the companion matrix of a monic
    // lattice polynomial with prescribed real roots, 3 the companion matrix
    // of a monic lattice polynomial with a complex pair, 4 upper triangular
    // on a lattice (already quasi-triangular, zero subdiagonal), 5 symmetric
    // on a lattice (an all-real spectrum), 6 a repeated real root, 7 the
    // zero matrix.
    std::vector<double> DrawMatrixRaw(oracle::Ctx& io, int32_t n, int32_t mode)
    {
        std::vector<double> m(static_cast<size_t>(n) * n, 0.0);
        auto at = [&m, n](int32_t r, int32_t c) -> double& { return m[c + n * r]; };

        if (mode == 0)
        {
            for (auto& e : m) { e = io.raw(-5.0, 5.0); }
        }
        else if (mode == 1)
        {
            for (auto& e : m) { e = io.rawInteger(-4, 4); }
        }
        else if (mode == 2 || mode == 3 || mode == 6)
        {
            // Monic polynomial x^n + c[n-1] x^{n-1} + ... + c[0] built from
            // prescribed roots, then its companion matrix. Mode 2 uses
            // distinct real lattice roots, mode 6 repeats the first root and
            // mode 3 multiplies in a quadratic with no real roots.
            std::vector<double> poly(1, 1.0);
            int32_t remaining = n;
            if (mode == 3)
            {
                double b = io.rawInteger(-3, 3);
                double c = io.rawInteger(1, 6);
                // x^2 - 2*b*x + (b*b + c) has the complex roots b +- i*sqrt(c).
                std::vector<double> quad{ b * b + c, -2.0 * b, 1.0 };
                std::vector<double> next(poly.size() + 2, 0.0);
                for (size_t i = 0; i < poly.size(); ++i)
                {
                    for (size_t j = 0; j < 3; ++j) { next[i + j] += poly[i] * quad[j]; }
                }
                poly = next;
                remaining -= 2;
            }
            double firstRoot = io.rawInteger(-4, 4);
            for (int32_t k = 0; k < remaining; ++k)
            {
                double root = firstRoot;
                if (k > 0) { root = (mode == 6 ? firstRoot : io.rawInteger(-4, 4)); }
                std::vector<double> next(poly.size() + 1, 0.0);
                for (size_t i = 0; i < poly.size(); ++i)
                {
                    next[i + 1] += poly[i];
                    next[i] -= poly[i] * root;
                }
                poly = next;
            }
            for (int32_t r = 0; r < n; ++r) { at(r, n - 1) = -poly[r]; }
            for (int32_t r = 1; r < n; ++r) { at(r, r - 1) = 1.0; }
        }
        else if (mode == 4)
        {
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = r; c < n; ++c) { at(r, c) = io.rawInteger(-4, 4); }
            }
        }
        else if (mode == 5)
        {
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = r; c < n; ++c)
                {
                    double value = io.rawInteger(-4, 4);
                    at(r, c) = value;
                    at(c, r) = value;
                }
            }
        }
        // mode 7 leaves the zero matrix.
        return m;
    }

    void RecordMatrix(oracle::Ctx& io, std::vector<double> const& m)
    {
        for (double e : m) { io.given(e); }
    }
}

// =======================================================================
// SymmetricEigensolver2x2.h
// =======================================================================

// The only entry point. Arithmetic is + - * / sqrt plus std::max and
// std::fabs, so the result is bit-identical.
ORACLE_CASE("SymmetricEigensolver2x2.solve")
{
    // Modes: 0 uniform, 1 lattice, 2 a00 == a11 and a01 == 0 (the
    // maxAbsComp == 0 branch, which is also the repeated-eigenvalue case),
    // 3 a01 == 0 with distinct diagonal, 4 signed zeros, 5 a lattice matrix
    // scaled by a power of two.
    int32_t mode = io.integer(0, 5);
    double a00 = 0.0, a01 = 0.0, a11 = 0.0;
    if (mode == 0)
    {
        a00 = io.raw(-5.0, 5.0);
        a01 = io.raw(-5.0, 5.0);
        a11 = io.raw(-5.0, 5.0);
    }
    else if (mode == 1)
    {
        a00 = io.rawInteger(-4, 4);
        a01 = io.rawInteger(-4, 4);
        a11 = io.rawInteger(-4, 4);
    }
    else if (mode == 2)
    {
        a00 = io.rawInteger(-4, 4);
        a11 = a00;
    }
    else if (mode == 3)
    {
        a00 = io.rawInteger(-4, 4);
        a11 = io.rawInteger(-4, 4);
    }
    else if (mode == 4)
    {
        double signedZero[2] = { 0.0, -0.0 };
        a00 = signedZero[io.rawInteger(0, 1)];
        a01 = signedZero[io.rawInteger(0, 1)];
        a11 = signedZero[io.rawInteger(0, 1)];
    }
    else
    {
        double scale = PowerOfTwo(io.rawInteger(-500, 500));
        a00 = scale * io.rawInteger(-4, 4);
        a01 = scale * io.rawInteger(-4, 4);
        a11 = scale * io.rawInteger(-4, 4);
    }
    io.given(a00);
    io.given(a01);
    io.given(a11);
    int32_t sortType = io.integer(-1, 1);

    std::array<double, 2> eval{};
    std::array<std::array<double, 2>, 2> evec{};
    SymmetricEigensolver2x2<double> solver;
    solver(a00, a01, a11, sortType, eval, evec);
    io.outReal(eval[0]);
    io.outReal(eval[1]);
    io.outReal(evec[0][0]);
    io.outReal(evec[0][1]);
    io.outReal(evec[1][0]);
    io.outReal(evec[1][1]);
}

// =======================================================================
// SymmetricEigensolver3x3.h
// =======================================================================

// SortEigenstuff is a public class of the header used by both solvers; it
// only permutes and negates, so it is exercised directly with lattice
// eigenvalues that produce plenty of ties across the six orderings.
ORACLE_CASE("SortEigenstuff.sort")
{
    std::array<double, 3> eval{};
    eval[0] = io.lattice(-3, 3);
    eval[1] = io.lattice(-3, 3);
    eval[2] = io.lattice(-3, 3);
    std::array<std::array<double, 3>, 3> evec{};
    for (size_t i = 0; i < 3; ++i)
    {
        for (size_t j = 0; j < 3; ++j) { evec[i][j] = io.lattice(-2, 2); }
    }
    bool isRotation = io.boolean();
    int32_t sortType = io.integer(-1, 1);
    SortEigenstuff<double>()(sortType, isRotation, eval, evec);
    OutSym3Result(io, eval, evec);
}

// The iterative solver. Both 'aggressive' settings and all three sort types
// are drawn. The generator stays inside the magnitude band on which the
// port's GetCosSin rescaling (issue #379) is inactive, so the port evaluates
// the upstream expression unchanged; the out-of-band scales are the subject
// of the deviation case below.
ORACLE_CASE("SymmetricEigensolver3x3.solve")
{
    Sym3 a = DrawSym3(io);
    bool aggressive = io.boolean();
    int32_t sortType = io.integer(-1, 1);
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    SymmetricEigensolver3x3<double> solver;
    int32_t iterations = solver(a.a00, a.a01, a.a02, a.a11, a.a12, a.a22,
        aggressive, sortType, eval, evec);
    io.outInt(iterations);
    OutSym3Result(io, eval, evec);
}

// Deviation: gtengine-js issue #379. Upstream's GetCosSin evaluates
// sqrt(u*u + v*v) with no rescaling, so for max(|u|,|v|) outside
// [2^-511, 2^511] the squares overflow to infinity or underflow into the
// subnormals and the returned 2-tuple is (0,0) or has lost most of its
// mantissa. The port rescales by max(|u|,|v|) exactly outside that band.
// The generator scales a lattice matrix by 2^k with 520 <= |k| <= 1000 and
// forces the off-diagonal entries that reach the first GetCosSin call to be
// nonzero, so the defect fires on every record.
ORACLE_CASE("SymmetricEigensolver3x3.solve.outOfBandScale")
{
    int32_t magnitude = io.rawInteger(520, 1000);
    int32_t sign = io.rawInteger(0, 1);
    double scale = PowerOfTwo(sign == 0 ? magnitude : -magnitude);
    int32_t nonzero0 = io.rawInteger(1, 4);
    int32_t nonzero1 = io.rawInteger(1, 4);
    int32_t flip0 = io.rawInteger(0, 1);
    int32_t flip1 = io.rawInteger(0, 1);
    double a00 = scale * io.rawInteger(-4, 4);
    double a01 = scale * io.rawInteger(-4, 4);
    double a02 = scale * (flip0 == 0 ? nonzero0 : -nonzero0);
    double a11 = scale * io.rawInteger(-4, 4);
    double a12 = scale * (flip1 == 0 ? nonzero1 : -nonzero1);
    double a22 = scale * io.rawInteger(-4, 4);
    io.given(a00);
    io.given(a01);
    io.given(a02);
    io.given(a11);
    io.given(a12);
    io.given(a22);
    bool aggressive = io.boolean();
    int32_t sortType = io.integer(-1, 1);
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    SymmetricEigensolver3x3<double> solver;
    int32_t iterations = solver(a00, a01, a02, a11, a12, a22,
        aggressive, sortType, eval, evec);
    io.outInt(iterations);
    OutSym3Result(io, eval, evec);
}

// The non-iterative solver. This is the only eigensolver in the group that
// calls the C math library: the eigenvalues come from std::acos and
// std::cos, whose MSVC and V8 implementations may differ by one ulp, and the
// eigenvectors are built from those eigenvalues. The comparison therefore
// carries a measured tolerance.
//
// The branch on 'halfDet' that selects which eigenvector is computed first
// is decided before the acos call and from + - * / alone, so it is the same
// on both sides. The remaining libm-decided control flow is the choice of
// the largest of d0, d1, d2 in ComputeEigenvector0 and of the largest row of
// M in ComputeEigenvector1; both are stable when the eigenvalues are well
// separated, so the generator rejects matrices whose eigenvalue gaps are
// below 1e-3 of the spectral radius, using upstream's own eigenvalues as the
// probe. Every quantity the test depends on is redrawn inside the loop, the
// loop is capped and the fallback is the best candidate seen.
ORACLE_CASE("NISymmetricEigensolver3x3.solve")
{
    NISymmetricEigensolver3x3<double> solver;
    Sym3 best{};
    double bestGap = -1.0;
    for (int32_t attempt = 0; attempt < 64; ++attempt)
    {
        // Modes 6 (zero) and 2 (diagonal) take the norm == 0 path, which is
        // covered by its own exact case, so they are excluded here.
        Sym3 candidate = DrawSym3Raw(io, 9);
        if (candidate.mode == 2 || candidate.mode == 6) { candidate.mode = 0; }
        if (candidate.mode == 0)
        {
            candidate.a00 = io.raw(-5.0, 5.0);
            candidate.a01 = io.raw(-5.0, 5.0);
            candidate.a02 = io.raw(-5.0, 5.0);
            candidate.a11 = io.raw(-5.0, 5.0);
            candidate.a12 = io.raw(-5.0, 5.0);
            candidate.a22 = io.raw(-5.0, 5.0);
        }

        std::array<double, 3> probeEval{};
        std::array<std::array<double, 3>, 3> probeEvec{};
        solver(candidate.a00, candidate.a01, candidate.a02, candidate.a11,
            candidate.a12, candidate.a22, 1, probeEval, probeEvec);
        double radius = std::max(std::fabs(probeEval[0]), std::fabs(probeEval[2]));
        double denom = std::max(1.0, radius);
        double gap0 = (probeEval[1] - probeEval[0]) / denom;
        double gap1 = (probeEval[2] - probeEval[1]) / denom;
        double gap = std::min(gap0, gap1);
        if (gap > bestGap)
        {
            bestGap = gap;
            best = candidate;
        }
        if (gap > 1e-3) { break; }
    }
    RecordSym3(io, best);
    int32_t sortType = io.integer(-1, 1);
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    solver(best.a00, best.a01, best.a02, best.a11, best.a12, best.a22,
        sortType, eval, evec);
    OutSym3Result(io, eval, evec);
}

// The two paths of the non-iterative solver that never reach the C math
// library: the zero matrix (maxAbsElement == 0) and a diagonal matrix
// (norm == 0). Both are bit-identical.
ORACLE_CASE("NISymmetricEigensolver3x3.solve.degenerate")
{
    // Mode 2 is diagonal, mode 6 is the zero matrix.
    int32_t pick = io.rawInteger(0, 3);
    Sym3 a = DrawSym3Raw(io, 9, pick == 0 ? 6 : 2);
    RecordSym3(io, a);
    int32_t sortType = io.integer(-1, 1);
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    NISymmetricEigensolver3x3<double> solver;
    solver(a.a00, a.a01, a.a02, a.a11, a.a12, a.a22, sortType, eval, evec);
    OutSym3Result(io, eval, evec);
}

// =======================================================================
// UnsymmetricEigenvalues.h
// =======================================================================

namespace
{
    void OutUnsymmetric(oracle::Ctx& io, int32_t n, uint32_t numIterations,
        UnsymmetricEigenvalues<double>& solver)
    {
        io.outInt(numIterations);
        uint32_t numEigenvalues = 0;
        std::vector<double> eigenvalues(static_cast<size_t>(n), 0.0);
        solver.GetEigenvalues(numEigenvalues, eigenvalues.data());
        io.outInt(numEigenvalues);
        for (uint32_t i = 0; i < numEigenvalues; ++i) { io.outReal(eigenvalues[i]); }
    }
}

// The QR solver. Everything is + - * / sqrt, and the eigenvalue sort is
// std::sort over at most six values, which MSVC implements as a (stable)
// insertion sort for at most 32 elements, matching the port's stable
// Array.prototype.sort. The generator rejects the inputs on which upstream
// drops its trailing eigenvalue (issue #42), using the verbatim
// ProbeUnsymmetric copy of upstream's own control flow as the separator;
// those inputs are the subject of the deviation case below.
ORACLE_CASE("UnsymmetricEigenvalues.solve")
{
    int32_t n = 0;
    uint32_t maxIterations = 0;
    std::vector<double> best;
    bool found = false;
    for (int32_t attempt = 0; attempt < 64; ++attempt)
    {
        n = io.rawInteger(3, 6);
        maxIterations = static_cast<uint32_t>(io.rawInteger(1, 64));
        int32_t mode = io.rawInteger(0, 7);
        std::vector<double> candidate = DrawMatrixRaw(io, n, mode);
        ProbeUnsymmetric probe(n, maxIterations);
        probe.Solve(candidate.data());
        if (!found)
        {
            best = candidate;
        }
        if (probe.TrailingFlag() == 1)
        {
            best = candidate;
            found = true;
            break;
        }
    }
    io.given(static_cast<double>(n));
    io.given(static_cast<double>(maxIterations));
    RecordMatrix(io, best);
    int32_t sortType = io.integer(-1, 1);

    UnsymmetricEigenvalues<double> solver(n, maxIterations);
    uint32_t numIterations = solver.Solve(best.data(), sortType);
    OutUnsymmetric(io, n, numIterations, solver);
}

// Deviation: gtengine-js issue #42. Upstream's packing loop stops at
// i < N-1, so when the final subdiagonal entry has decoupled the 1x1 block
// at (N-1,N-1) is never reported and one real eigenvalue is silently lost.
// The generator keeps only the inputs where upstream's own flags say the
// trailing block decoupled, so the port reports one eigenvalue more than
// upstream on every record.
ORACLE_CASE("UnsymmetricEigenvalues.solve.trailingBlock")
{
    int32_t n = 0;
    uint32_t maxIterations = 0;
    std::vector<double> best;
    for (int32_t attempt = 0; attempt < 64; ++attempt)
    {
        n = io.rawInteger(3, 6);
        maxIterations = static_cast<uint32_t>(io.rawInteger(8, 64));
        // Modes 2, 4, 5 and 6 have an all-real spectrum and converge to a
        // fully quasi-triangular form with a trailing 1x1 block.
        int32_t const realModes[4] = { 2, 4, 5, 6 };
        int32_t mode = realModes[io.rawInteger(0, 3)];
        std::vector<double> candidate = DrawMatrixRaw(io, n, mode);
        ProbeUnsymmetric probe(n, maxIterations);
        probe.Solve(candidate.data());
        best = candidate;
        if (probe.TrailingFlag() == 0) { break; }
    }
    io.given(static_cast<double>(n));
    io.given(static_cast<double>(maxIterations));
    RecordMatrix(io, best);
    int32_t sortType = io.integer(-1, 1);

    UnsymmetricEigenvalues<double> solver(n, maxIterations);
    uint32_t numIterations = solver.Solve(best.data(), sortType);
    OutUnsymmetric(io, n, numIterations, solver);
}

// The constructor rejects size < 3 and maxIterations == 0 by leaving the
// solver at size 0; Solve then returns 0 and GetEigenvalues reports none.
ORACLE_CASE("UnsymmetricEigenvalues.solve.invalidSize")
{
    int32_t n = io.rawInteger(0, 3);
    int32_t rawMax = io.rawInteger(0, 2);
    // n == 3 is a valid size, so pair it only with the invalid budget 0.
    uint32_t maxIterations = static_cast<uint32_t>(n >= 3 ? 0 : rawMax);
    io.given(static_cast<double>(n));
    io.given(static_cast<double>(maxIterations));
    int32_t storage = std::max(n, 1);
    std::vector<double> m(static_cast<size_t>(storage) * storage, 0.0);
    for (auto& e : m) { e = io.real(-4.0, 4.0); }
    int32_t sortType = io.integer(-1, 1);

    UnsymmetricEigenvalues<double> solver(n, maxIterations);
    uint32_t numIterations = solver.Solve(m.data(), sortType);
    io.outInt(numIterations);
    uint32_t numEigenvalues = 0;
    std::vector<double> eigenvalues(static_cast<size_t>(storage), 0.0);
    solver.GetEigenvalues(numEigenvalues, eigenvalues.data());
    io.outInt(numEigenvalues);
    for (uint32_t i = 0; i < numEigenvalues; ++i) { io.outReal(eigenvalues[i]); }
}

// The documented non-convergence example of UPSTREAM-FINDINGS (issue #476):
// A = [[0,8,0],[8,0,8],[0,8,0]] has the well separated real eigenvalues 0
// and +-8*sqrt(2), but FrancisQRStep has no exceptional shift and the
// iteration cycles, burning the whole budget and reporting no eigenvalues.
// The behaviour is preserved by the port, so the two sides agree. The
// off-diagonal magnitude and the iteration budget are drawn.
ORACLE_CASE("UnsymmetricEigenvalues.solve.cycling")
{
    double t = 1.0;
    uint32_t maxIterations = 1;
    for (int32_t attempt = 0; attempt < 32; ++attempt)
    {
        t = static_cast<double>(io.rawInteger(1, 20));
        maxIterations = static_cast<uint32_t>(io.rawInteger(1, 64));
        std::vector<double> candidate{ 0.0, t, 0.0, t, 0.0, t, 0.0, t, 0.0 };
        ProbeUnsymmetric probe(3, maxIterations);
        probe.Solve(candidate.data());
        // Only a still-coupled trailing block keeps issue #42 out of this
        // case; the cycling example satisfies it for every t seen.
        if (probe.TrailingFlag() == 1) { break; }
    }
    io.given(t);
    io.given(static_cast<double>(maxIterations));
    int32_t sortType = io.integer(-1, 1);
    std::vector<double> m{ 0.0, t, 0.0, t, 0.0, t, 0.0, t, 0.0 };
    UnsymmetricEigenvalues<double> solver(3, maxIterations);
    uint32_t numIterations = solver.Solve(m.data(), sortType);
    OutUnsymmetric(io, 3, numIterations, solver);
}
