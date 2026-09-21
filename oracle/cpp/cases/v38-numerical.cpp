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

// Two more exact reproductions of the preserved non-convergence defect
// (issue #476): FrancisQRStep has no exceptional shift, so the iteration can
// cycle on a well separated real spectrum. Variant 0 is a nearly symmetric
// 3x3 found by a CI run of the port; variant 1 is the same matrix with a00
// shifted by 1e-13, which converges in 33 iterations and is the contrast;
// variant 2 is the zero-diagonal matrix of the UPSTREAM-FINDINGS entry. The
// probe keeps issue #42 (the dropped trailing eigenvalue) out of the case by
// falling back to the cycling variant whenever the trailing block decouples,
// which is why variant 1 survives only with a budget below its 33
// iterations: once it converges, upstream drops A(2,2) and the port does
// not, which is the subject of the deviation case above.
ORACLE_CASE("UnsymmetricEigenvalues.solve.nonConvergence")
{
    double const cycling[9] =
    {
        -8.0, -7.99999999999999, 0.0,
        -7.99999999999999, 0.250000000000008, 7.99999999999999,
        0.0, 7.99999999999999, -8.0
    };
    double const zeroDiagonal[9] = { 0.0, 8.0, 0.0, 8.0, 0.0, 8.0, 0.0, 8.0, 0.0 };
    uint32_t maxIterations = static_cast<uint32_t>(io.rawInteger(1, 64));
    int32_t variant = io.rawInteger(0, 2);
    std::vector<double> m(9, 0.0);
    for (int32_t attempt = 0; attempt < 3; ++attempt)
    {
        if (variant == 2)
        {
            for (int32_t i = 0; i < 9; ++i) { m[i] = zeroDiagonal[i]; }
        }
        else
        {
            for (int32_t i = 0; i < 9; ++i) { m[i] = cycling[i]; }
            if (variant == 1) { m[0] = cycling[0] + 1e-13; }
        }
        ProbeUnsymmetric probe(3, maxIterations);
        probe.Solve(m.data());
        if (probe.TrailingFlag() == 1) { break; }
        variant = 0;
    }
    io.given(static_cast<double>(variant));
    io.given(static_cast<double>(maxIterations));
    for (double e : m) { io.given(e); }
    int32_t sortType = io.integer(-1, 1);
    UnsymmetricEigenvalues<double> solver(3, maxIterations);
    uint32_t numIterations = solver.Solve(m.data(), sortType);
    OutUnsymmetric(io, 3, numIterations, solver);
}

// =======================================================================
// BandedMatrix.h
// =======================================================================

namespace
{
    // A banded matrix of the given shape. Modes: 0 uniform entries,
    // 1 lattice entries, 2 symmetric positive definite built as B*B^T from a
    // lower-banded lattice B (so CholeskyFactor succeeds), 3 symmetric on a
    // lattice but not definite, 4 a symmetric lattice matrix with a zero
    // diagonal entry (the CholeskyFactor and ComputeInverse failure
    // returns). Only the stored band entries are drawn; the layout is
    // diagonal first, then each lower band, then each upper band, which is
    // what RecordBanded writes.
    struct Banded
    {
        int32_t size = 0, numLBands = 0, numUBands = 0, mode = 0;
        std::vector<double> d;
        std::vector<std::vector<double>> l, u;
    };

    void ResizeBanded(Banded& b)
    {
        b.d.assign(static_cast<size_t>(b.size), 0.0);
        b.l.clear();
        b.u.clear();
        for (int32_t k = 0; k < b.numLBands; ++k)
        {
            b.l.emplace_back(static_cast<size_t>(b.size - 1 - k), 0.0);
        }
        for (int32_t k = 0; k < b.numUBands; ++k)
        {
            b.u.emplace_back(static_cast<size_t>(b.size - 1 - k), 0.0);
        }
    }

    Banded DrawBandedRaw(oracle::Ctx& io, int32_t mode, int32_t minSize = 2,
        int32_t maxSize = 6)
    {
        Banded b;
        b.mode = mode;
        b.size = io.rawInteger(minSize, maxSize);
        int32_t bands = io.rawInteger(0, b.size - 1);
        b.numLBands = bands;
        b.numUBands = bands;
        if (mode == 0 || mode == 1)
        {
            // The general (possibly non-symmetric) shape: draw the band
            // counts independently.
            b.numLBands = io.rawInteger(0, b.size - 1);
            b.numUBands = io.rawInteger(0, b.size - 1);
        }
        ResizeBanded(b);

        if (mode == 0 || mode == 1)
        {
            bool lattice = (mode == 1);
            for (auto& e : b.d) { e = lattice ? io.rawInteger(-4, 4) : io.raw(-5.0, 5.0); }
            for (auto& band : b.l)
            {
                for (auto& e : band) { e = lattice ? io.rawInteger(-4, 4) : io.raw(-5.0, 5.0); }
            }
            for (auto& band : b.u)
            {
                for (auto& e : band) { e = lattice ? io.rawInteger(-4, 4) : io.raw(-5.0, 5.0); }
            }
        }
        else if (mode == 2)
        {
            // A = B*B^T with B lower triangular and banded with 'bands'
            // subdiagonals, so A is symmetric positive definite and banded
            // with the same band count. The diagonal of B is nonzero on a
            // lattice, so A is exactly representable.
            int32_t const n = b.size;
            std::vector<std::vector<double>> B(static_cast<size_t>(n),
                std::vector<double>(static_cast<size_t>(n), 0.0));
            for (int32_t r = 0; r < n; ++r)
            {
                B[r][r] = io.rawInteger(1, 4);
                for (int32_t c = std::max(0, r - bands); c < r; ++c)
                {
                    B[r][c] = io.rawInteger(-3, 3);
                }
            }
            auto entry = [&B, n](int32_t r, int32_t c)
            {
                double sum = 0.0;
                for (int32_t k = 0; k < n; ++k) { sum += B[r][k] * B[c][k]; }
                return sum;
            };
            for (int32_t r = 0; r < n; ++r) { b.d[r] = entry(r, r); }
            for (int32_t k = 0; k < bands; ++k)
            {
                for (int32_t r = 0; r + k + 1 < n; ++r)
                {
                    double value = entry(r + k + 1, r);
                    b.l[k][r] = value;
                    b.u[k][r] = value;
                }
            }
        }
        else
        {
            // Symmetric on a lattice, indefinite. Mode 4 zeroes one diagonal
            // entry so that the pivot tests fail.
            for (auto& e : b.d) { e = io.rawInteger(-4, 4); }
            for (int32_t k = 0; k < bands; ++k)
            {
                for (size_t r = 0; r < b.l[k].size(); ++r)
                {
                    double value = io.rawInteger(-3, 3);
                    b.l[k][r] = value;
                    b.u[k][r] = value;
                }
            }
            if (mode == 4)
            {
                b.d[io.rawInteger(0, b.size - 1)] = 0.0;
            }
        }
        return b;
    }

    // Records the shape (size, numLBands, numUBands, mode) and then every
    // stored entry: the diagonal, the lower bands, the upper bands.
    void RecordBanded(oracle::Ctx& io, Banded const& b)
    {
        io.given(static_cast<double>(b.size));
        io.given(static_cast<double>(b.numLBands));
        io.given(static_cast<double>(b.numUBands));
        io.given(static_cast<double>(b.mode));
        for (double e : b.d) { io.given(e); }
        for (auto const& band : b.l)
        {
            for (double e : band) { io.given(e); }
        }
        for (auto const& band : b.u)
        {
            for (double e : band) { io.given(e); }
        }
    }

    BandedMatrix<double> MakeBanded(Banded const& b)
    {
        BandedMatrix<double> matrix(b.size, b.numLBands, b.numUBands);
        matrix.GetDBand() = b.d;
        for (int32_t k = 0; k < b.numLBands; ++k) { matrix.GetLBands()[k] = b.l[k]; }
        for (int32_t k = 0; k < b.numUBands; ++k) { matrix.GetUBands()[k] = b.u[k]; }
        return matrix;
    }

    // Every stored entry of the matrix, in the order RecordBanded uses.
    void OutBanded(oracle::Ctx& io, BandedMatrix<double> const& matrix)
    {
        for (double e : matrix.GetDBand()) { io.outReal(e); }
        for (auto const& band : matrix.GetLBands())
        {
            for (double e : band) { io.outReal(e); }
        }
        for (auto const& band : matrix.GetUBands())
        {
            for (double e : band) { io.outReal(e); }
        }
    }
}

// The element accessors, including the in-range indices outside the stored
// bands and the out-of-range indices, both of which upstream services from
// its mZero scratch member. Writes through that reference are discarded.
ORACLE_CASE("BandedMatrix.accessor")
{
    int32_t mode = io.rawInteger(0, 1);
    Banded b = DrawBandedRaw(io, mode);
    RecordBanded(io, b);
    BandedMatrix<double> matrix = MakeBanded(b);
    io.outInt(matrix.GetSize());

    // Read every entry of the full square, plus the four out-of-range
    // probes.
    for (int32_t r = -1; r <= b.size; ++r)
    {
        for (int32_t c = -1; c <= b.size; ++c)
        {
            io.outReal(matrix(r, c));
        }
    }

    // Write a distinctive value at every index and read the whole matrix
    // back: entries outside the bands must be unchanged.
    for (int32_t r = -1; r <= b.size; ++r)
    {
        for (int32_t c = -1; c <= b.size; ++c)
        {
            matrix(r, c) = static_cast<double>(1 + r * (b.size + 2) + c);
        }
    }
    OutBanded(io, matrix);
}

// The Cholesky factorization. Mode 2 is positive definite and succeeds;
// modes 3 and 4 exercise the 'diagonal <= 0' failure return, and the
// non-symmetric modes 0 and 1 reach it too. The factored matrix is emitted
// whatever the return value, because upstream leaves the partial
// factorization in place.
ORACLE_CASE("BandedMatrix.choleskyFactor")
{
    int32_t mode = io.rawInteger(0, 4);
    Banded b = DrawBandedRaw(io, mode);
    RecordBanded(io, b);
    BandedMatrix<double> matrix = MakeBanded(b);
    bool success = matrix.CholeskyFactor();
    io.outBool(success);
    OutBanded(io, matrix);
}

// SolveSystem with a vector right-hand side. Note that upstream's
// vector-valued SolveLower/SolveUpper divide by the pivot while the
// matrix-valued ones multiply by its reciprocal, so the two overloads are
// different computations and both are covered.
ORACLE_CASE("BandedMatrix.solveSystem.vector")
{
    int32_t mode = io.rawInteger(0, 4);
    Banded b = DrawBandedRaw(io, mode);
    RecordBanded(io, b);
    std::vector<double> rhs(static_cast<size_t>(b.size), 0.0);
    for (auto& e : rhs) { e = io.real(-6.0, 6.0); }
    BandedMatrix<double> matrix = MakeBanded(b);
    bool success = matrix.SolveSystem(rhs.data());
    io.outBool(success);
    for (double e : rhs) { io.outReal(e); }
    OutBanded(io, matrix);
}

// SolveSystem with a matrix right-hand side, in both storage orders.
ORACLE_CASE("BandedMatrix.solveSystem.matrix")
{
    int32_t mode = io.rawInteger(0, 4);
    Banded b = DrawBandedRaw(io, mode);
    RecordBanded(io, b);
    int32_t numBColumns = io.integer(1, 3);
    bool rowMajor = io.boolean();
    std::vector<double> rhs(static_cast<size_t>(b.size) * numBColumns, 0.0);
    for (auto& e : rhs) { e = io.real(-6.0, 6.0); }
    BandedMatrix<double> matrix = MakeBanded(b);
    bool success = rowMajor
        ? matrix.SolveSystem<true>(rhs.data(), numBColumns)
        : matrix.SolveSystem<false>(rhs.data(), numBColumns);
    io.outBool(success);
    for (double e : rhs) { io.outReal(e); }
    OutBanded(io, matrix);
}

// ComputeInverse, in both storage orders, including the 'diag == 0' failure
// return. Upstream leaves 'inverse' in a partially reduced state when it
// fails, and the caller is told not to use it, so only the success flag is
// emitted on that path.
ORACLE_CASE("BandedMatrix.computeInverse")
{
    int32_t mode = io.rawInteger(0, 4);
    Banded b = DrawBandedRaw(io, mode);
    RecordBanded(io, b);
    bool rowMajor = io.boolean();
    std::vector<double> inverse(static_cast<size_t>(b.size) * b.size, 0.0);
    BandedMatrix<double> matrix = MakeBanded(b);
    bool success = rowMajor
        ? matrix.ComputeInverse<true>(inverse.data())
        : matrix.ComputeInverse<false>(inverse.data());
    io.outBool(success);
    if (success)
    {
        for (double e : inverse) { io.outReal(e); }
    }
    // ComputeInverse works on a copy, so the matrix itself is unchanged.
    OutBanded(io, matrix);
}

// UPSTREAM DEFECT found by this oracle (both sides agree bit for bit, and
// both are wrong). ComputeInverse runs Gaussian elimination with NO pivoting
// and rejects a matrix only when the pivot is exactly zero. When a leading
// principal minor of the banded matrix vanishes, the pivot is zero in exact
// arithmetic but a round-off residue in floating point, so the test passes,
// the multiplier is about 1e16 and the reported inverse is wrong in the
// second digit while the return value is 'true'. The 6x6 matrix below is a
// deep-run witness (2 of 1605 integer records of the case above): its fourth
// leading principal minor is 0 and its determinant is -656, so it is
// invertible, yet the reported inverse differs from the exact rational
// inverse by 4.4e-2 relative. The port preserves the behaviour, so the case
// is compared bit for bit; see oracle/reports/v38-numerical.md.
ORACLE_CASE("BandedMatrix.computeInverse.zeroLeadingMinor")
{
    Banded b;
    b.size = 6;
    b.numLBands = 2;
    b.numUBands = 4;
    b.mode = 9;
    ResizeBanded(b);
    // A power-of-two scale keeps every entry exactly representable, so the
    // elimination reproduces the same cancellation at every scale.
    double s = PowerOfTwo(io.rawInteger(-8, 8));
    double const dBand[6] = { 1, 2, 0, 1, 2, -2 };
    double const lBand0[5] = { 4, -4, 0, -3, 4 };
    double const lBand1[4] = { 0, 2, 4, -2 };
    double const uBand0[5] = { 0, -2, -2, 4, -3 };
    double const uBand1[4] = { 3, -1, -2, -2 };
    double const uBand2[3] = { -3, -1, 3 };
    double const uBand3[2] = { 3, 3 };
    for (int32_t i = 0; i < 6; ++i) { b.d[i] = s * dBand[i]; }
    for (int32_t i = 0; i < 5; ++i) { b.l[0][i] = s * lBand0[i]; }
    for (int32_t i = 0; i < 4; ++i) { b.l[1][i] = s * lBand1[i]; }
    for (int32_t i = 0; i < 5; ++i) { b.u[0][i] = s * uBand0[i]; }
    for (int32_t i = 0; i < 4; ++i) { b.u[1][i] = s * uBand1[i]; }
    for (int32_t i = 0; i < 3; ++i) { b.u[2][i] = s * uBand2[i]; }
    for (int32_t i = 0; i < 2; ++i) { b.u[3][i] = s * uBand3[i]; }
    RecordBanded(io, b);
    bool rowMajor = io.boolean();
    std::vector<double> inverse(36, 0.0);
    BandedMatrix<double> matrix = MakeBanded(b);
    bool success = rowMajor
        ? matrix.ComputeInverse<true>(inverse.data())
        : matrix.ComputeInverse<false>(inverse.data());
    io.outBool(success);
    for (double e : inverse) { io.outReal(e); }
}

// The constructor rejects size <= 0 and band counts outside [0, size) by
// collapsing the matrix to size 0; every accessor then reports zero.
ORACLE_CASE("BandedMatrix.invalidShape")
{
    int32_t size = io.integer(-1, 3);
    int32_t numLBands = io.integer(-1, 4);
    int32_t numUBands = io.integer(-1, 4);
    BandedMatrix<double> matrix(size, numLBands, numUBands);
    io.outInt(matrix.GetSize());
    io.outInt(static_cast<int32_t>(matrix.GetDBand().size()));
    io.outInt(static_cast<int32_t>(matrix.GetLBands().size()));
    io.outInt(static_cast<int32_t>(matrix.GetUBands().size()));
    io.outBool(matrix.CholeskyFactor());
    for (int32_t r = 0; r < 3; ++r)
    {
        for (int32_t c = 0; c < 3; ++c) { io.outReal(matrix(r, c)); }
    }
}

// =======================================================================
// GaussianElimination.h
// =======================================================================

namespace
{
    // An NxN matrix in the row-major order that GTE_USE_ROW_MAJOR selects
    // (the port's default and the only order covered here). Modes:
    // 0 uniform, 1 small lattice, 2 a permutation matrix, 3 exactly singular
    // (a repeated row), 4 nearly singular (a row perturbed by 2^-40),
    // 5 diagonal with mixed magnitudes, 6 the zero matrix, 7 a row of zeros.
    std::vector<double> DrawSquareRaw(oracle::Ctx& io, int32_t n, int32_t mode)
    {
        std::vector<double> m(static_cast<size_t>(n) * n, 0.0);
        auto at = [&m, n](int32_t r, int32_t c) -> double& { return m[c + n * r]; };
        if (mode == 0)
        {
            for (auto& e : m) { e = io.raw(-5.0, 5.0); }
        }
        else if (mode == 1 || mode == 3 || mode == 4)
        {
            for (auto& e : m) { e = io.rawInteger(-4, 4); }
            if (mode != 1 && n >= 2)
            {
                int32_t source = io.rawInteger(0, n - 1);
                int32_t target = (source + 1 + io.rawInteger(0, n - 2)) % n;
                double perturbation = (mode == 4 ? PowerOfTwo(-40) : 0.0);
                for (int32_t c = 0; c < n; ++c)
                {
                    at(target, c) = at(source, c) + (c == 0 ? perturbation : 0.0);
                }
            }
        }
        else if (mode == 2)
        {
            std::vector<int32_t> permutation(static_cast<size_t>(n));
            for (int32_t i = 0; i < n; ++i) { permutation[i] = i; }
            for (int32_t i = n - 1; i > 0; --i)
            {
                std::swap(permutation[i], permutation[io.rawInteger(0, i)]);
            }
            for (int32_t r = 0; r < n; ++r) { at(r, permutation[r]) = 1.0; }
        }
        else if (mode == 5)
        {
            for (int32_t r = 0; r < n; ++r)
            {
                at(r, r) = PowerOfTwo(io.rawInteger(-60, 60)) * io.rawInteger(1, 4);
            }
        }
        else if (mode == 7)
        {
            for (auto& e : m) { e = io.rawInteger(-4, 4); }
            int32_t zeroRow = io.rawInteger(0, n - 1);
            for (int32_t c = 0; c < n; ++c) { at(zeroRow, c) = 0.0; }
        }
        // mode 6 leaves the zero matrix.
        return m;
    }
}

// The full operator(): the inverse, the determinant, the vector solve and
// the matrix solve in one call, which is the widest path through the
// routine. Singular inputs take the 'maxValue == zero' early return, which
// zero-fills every output.
ORACLE_CASE("GaussianElimination.compute")
{
    int32_t n = io.rawInteger(1, 5);
    int32_t mode = io.rawInteger(0, 7);
    std::vector<double> m = DrawSquareRaw(io, n, mode);
    io.given(static_cast<double>(n));
    io.given(static_cast<double>(mode));
    for (double e : m) { io.given(e); }
    int32_t numCols = io.integer(1, 3);
    std::vector<double> b(static_cast<size_t>(n), 0.0);
    for (auto& e : b) { e = io.real(-6.0, 6.0); }
    std::vector<double> c(static_cast<size_t>(n) * numCols, 0.0);
    for (auto& e : c) { e = io.real(-6.0, 6.0); }

    std::vector<double> inverse(static_cast<size_t>(n) * n, 0.0);
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    std::vector<double> y(static_cast<size_t>(n) * numCols, 0.0);
    double determinant = 0.0;
    GaussianElimination<double> solver;
    bool invertible = solver(n, m.data(), inverse.data(), determinant,
        b.data(), x.data(), c.data(), numCols, y.data());
    io.outBool(invertible);
    io.outReal(determinant);
    for (double e : inverse) { io.outReal(e); }
    for (double e : x) { io.outReal(e); }
    for (double e : y) { io.outReal(e); }
}

// The narrower calls: no inverse wanted (upstream then eliminates in a local
// scratch buffer and skips the final un-permutation of the rows), and only
// one of the two right-hand sides supplied.
ORACLE_CASE("GaussianElimination.compute.partial")
{
    int32_t n = io.rawInteger(1, 5);
    int32_t mode = io.rawInteger(0, 7);
    std::vector<double> m = DrawSquareRaw(io, n, mode);
    io.given(static_cast<double>(n));
    io.given(static_cast<double>(mode));
    for (double e : m) { io.given(e); }
    // selector: 0 determinant only, 1 vector solve, 2 matrix solve,
    // 3 inverse and determinant only.
    int32_t selector = io.integer(0, 3);
    int32_t numCols = io.integer(1, 3);
    std::vector<double> b(static_cast<size_t>(n), 0.0);
    for (auto& e : b) { e = io.real(-6.0, 6.0); }
    std::vector<double> c(static_cast<size_t>(n) * numCols, 0.0);
    for (auto& e : c) { e = io.real(-6.0, 6.0); }

    std::vector<double> inverse(static_cast<size_t>(n) * n, 0.0);
    std::vector<double> x(static_cast<size_t>(n), 0.0);
    std::vector<double> y(static_cast<size_t>(n) * numCols, 0.0);
    double determinant = 0.0;
    GaussianElimination<double> solver;
    bool invertible = false;
    if (selector == 0)
    {
        invertible = solver(n, m.data(), nullptr, determinant,
            nullptr, nullptr, nullptr, 0, nullptr);
    }
    else if (selector == 1)
    {
        invertible = solver(n, m.data(), nullptr, determinant,
            b.data(), x.data(), nullptr, 0, nullptr);
    }
    else if (selector == 2)
    {
        invertible = solver(n, m.data(), nullptr, determinant,
            nullptr, nullptr, c.data(), numCols, y.data());
    }
    else
    {
        invertible = solver(n, m.data(), inverse.data(), determinant,
            nullptr, nullptr, nullptr, 0, nullptr);
    }
    io.outBool(invertible);
    io.outReal(determinant);
    if (selector == 3)
    {
        for (double e : inverse) { io.outReal(e); }
    }
    if (selector == 1)
    {
        for (double e : x) { io.outReal(e); }
    }
    if (selector == 2)
    {
        for (double e : y) { io.outReal(e); }
    }
}

// Preserved upstream defect (gtengine-js issue #375). A matrix whose entries
// are all of subnormal magnitude passes the pivot test (the largest entry is
// nonzero), but 1 / pivot overflows to infinity and the row operations
// produce infinity * 0 = NaN, so the call reports invertible = true with NaN
// entries. The port preserves this deliberately, so the two sides are
// compared bit for bit; the discrete 'invertible' output keeps the record
// from being vacuous when every real output is NaN.
ORACLE_CASE("GaussianElimination.compute.subnormal")
{
    int32_t n = io.rawInteger(2, 4);
    std::vector<double> m(static_cast<size_t>(n) * n, 0.0);
    for (auto& e : m)
    {
        e = PowerOfTwo(-1070) * io.rawInteger(-8, 8);
    }
    io.given(static_cast<double>(n));
    for (double e : m) { io.given(e); }
    std::vector<double> inverse(static_cast<size_t>(n) * n, 0.0);
    double determinant = 0.0;
    GaussianElimination<double> solver;
    bool invertible = solver(n, m.data(), inverse.data(), determinant,
        nullptr, nullptr, nullptr, 0, nullptr);
    io.outBool(invertible);
    io.outReal(determinant);
    for (double e : inverse) { io.outReal(e); }
}

// The LogError precondition. Upstream rejects numRows <= 0, a B without an
// X, a C without a Y and a C with numCols < 1; the port replaces the pointer
// pairs with optional inputs, so only the two conditions it can express are
// covered here: numRows <= 0 and a C with numCols < 1. Upstream throws and
// the port must throw on the same records.
ORACLE_CASE("GaussianElimination.compute.invalidInput")
{
    int32_t which = io.integer(0, 1);
    int32_t n = io.integer(-1, 3);
    std::vector<double> m(static_cast<size_t>(std::max(n, 1)) * std::max(n, 1), 0.0);
    for (auto& e : m) { e = io.real(-4.0, 4.0); }
    std::vector<double> c(static_cast<size_t>(std::max(n, 1)), 0.0);
    for (auto& e : c) { e = io.real(-4.0, 4.0); }
    std::vector<double> y(c.size(), 0.0);
    double determinant = 0.0;
    GaussianElimination<double> solver;
    bool invertible = false;
    if (which == 0)
    {
        // numRows <= 0 when n <= 0; otherwise a legitimate call.
        invertible = solver(n, m.data(), nullptr, determinant,
            nullptr, nullptr, nullptr, 0, nullptr);
    }
    else
    {
        // C with numCols < 1.
        invertible = solver(n, m.data(), nullptr, determinant,
            nullptr, nullptr, c.data(), 0, y.data());
    }
    io.outBool(invertible);
    io.outReal(determinant);
}

// =======================================================================
// Integration.h
// =======================================================================

namespace
{
    // The integrands are written identically on both sides so that the
    // comparison is of the quadrature rule, not of the integrand. Kind 0 is
    // a polynomial in Horner form, kind 1 is a rational function whose
    // denominator is bounded away from zero, and kind 2 is the only libm
    // integrand in the group (std::exp), used by its own case.
    struct Integrand
    {
        int32_t kind = 0;
        int32_t degree = 0;
        std::vector<double> coefficient;
    };

    Integrand DrawIntegrand(oracle::Ctx& io, int32_t kind)
    {
        Integrand f;
        f.kind = kind;
        f.degree = io.integer(0, 5);
        f.coefficient.resize(static_cast<size_t>(f.degree) + 1);
        for (auto& e : f.coefficient) { e = io.real(-3.0, 3.0); }
        return f;
    }

    double Evaluate(Integrand const& f, double t)
    {
        double result = f.coefficient[f.degree];
        for (int32_t i = f.degree - 1; i >= 0; --i)
        {
            result = result * t + f.coefficient[i];
        }
        if (f.kind == 1)
        {
            result = result / (t * t + 1.0);
        }
        else if (f.kind == 2)
        {
            result = std::exp(-t * t) * result;
        }
        return result;
    }
}

ORACLE_CASE("Integration.trapezoidRule")
{
    int32_t numSamples = io.integer(2, 24);
    double a = io.real(-3.0, 3.0);
    double b = io.real(-3.0, 3.0);
    int32_t kind = io.integer(0, 1);
    Integrand f = DrawIntegrand(io, kind);
    auto integrand = [&f](double t) { return Evaluate(f, t); };
    io.outReal(Integration<double>::TrapezoidRule(numSamples, a, b, integrand));
}

ORACLE_CASE("Integration.romberg")
{
    int32_t order = io.integer(1, 10);
    double a = io.real(-3.0, 3.0);
    double b = io.real(-3.0, 3.0);
    int32_t kind = io.integer(0, 1);
    Integrand f = DrawIntegrand(io, kind);
    auto integrand = [&f](double t) { return Evaluate(f, t); };
    io.outReal(Integration<double>::Romberg(order, a, b, integrand));
}

// ComputeQuadratureInfo builds the Legendre polynomial of the requested
// degree and finds its roots with RootsPolynomial::Find, which is bisection
// on + - * / alone, then solves for the coefficients with the subset-product
// recursion. All of it is exact arithmetic.
ORACLE_CASE("Integration.computeQuadratureInfo")
{
    int32_t degree = io.integer(2, 9);
    std::vector<double> roots, coefficients;
    Integration<double>::ComputeQuadratureInfo(degree, roots, coefficients);
    io.outInt(static_cast<int32_t>(roots.size()));
    for (double e : roots) { io.outReal(e); }
    io.outInt(static_cast<int32_t>(coefficients.size()));
    for (double e : coefficients) { io.outReal(e); }
}

ORACLE_CASE("Integration.gaussianQuadrature")
{
    int32_t degree = io.integer(2, 9);
    double a = io.real(-3.0, 3.0);
    double b = io.real(-3.0, 3.0);
    int32_t kind = io.integer(0, 1);
    Integrand f = DrawIntegrand(io, kind);
    std::vector<double> roots, coefficients;
    Integration<double>::ComputeQuadratureInfo(degree, roots, coefficients);
    auto integrand = [&f](double t) { return Evaluate(f, t); };
    io.outReal(Integration<double>::GaussianQuadrature(roots, coefficients,
        a, b, integrand));
    io.outInt(static_cast<int32_t>(roots.size()));
    for (double e : roots) { io.outReal(e); }
    for (double e : coefficients) { io.outReal(e); }
}

// The one case with a libm integrand: std::exp is accurate to within an ulp
// in both runtimes but is not required to be correctly rounded, so the
// quadrature sums differ in the last bits. Compared with a tolerance.
ORACLE_CASE("Integration.libmIntegrand")
{
    int32_t order = io.integer(1, 8);
    int32_t numSamples = io.integer(2, 24);
    double a = io.real(-2.0, 2.0);
    double b = io.real(-2.0, 2.0);
    Integrand f = DrawIntegrand(io, 2);
    auto integrand = [&f](double t) { return Evaluate(f, t); };
    io.outReal(Integration<double>::TrapezoidRule(numSamples, a, b, integrand));
    io.outReal(Integration<double>::Romberg(order, a, b, integrand));
}

// =======================================================================
// LCPSolver.h
// =======================================================================

namespace
{
    // The LCP w = q + M*z. Modes: 0 uniform M and q, 1 lattice M and q,
    // 2 M positive definite (M = A*A^T + I on a lattice, so a solution
    // always exists), 3 M = 0, 4 M with only nonpositive entries and q < 0
    // (provably infeasible: w = q + M*z < 0 for every z >= 0), 5 q >= 0 (the
    // trivial solution), 6 a lattice M with a zero row.
    struct Lcp
    {
        int32_t n = 0, mode = 0;
        std::vector<double> q, m;
    };

    Lcp DrawLcpRaw(oracle::Ctx& io, int32_t mode, int32_t minN = 1, int32_t maxN = 6)
    {
        Lcp lcp;
        lcp.mode = mode;
        lcp.n = io.rawInteger(minN, maxN);
        int32_t const n = lcp.n;
        lcp.q.assign(static_cast<size_t>(n), 0.0);
        lcp.m.assign(static_cast<size_t>(n) * n, 0.0);
        auto at = [&lcp, n](int32_t r, int32_t c) -> double& { return lcp.m[c + n * r]; };

        if (mode == 0)
        {
            for (auto& e : lcp.m) { e = io.raw(-4.0, 4.0); }
            for (auto& e : lcp.q) { e = io.raw(-4.0, 4.0); }
        }
        else if (mode == 1 || mode == 6)
        {
            for (auto& e : lcp.m) { e = io.rawInteger(-3, 3); }
            for (auto& e : lcp.q) { e = io.rawInteger(-4, 4); }
            if (mode == 6)
            {
                int32_t zeroRow = io.rawInteger(0, n - 1);
                for (int32_t c = 0; c < n; ++c) { at(zeroRow, c) = 0.0; }
            }
        }
        else if (mode == 2)
        {
            std::vector<double> a(static_cast<size_t>(n) * n, 0.0);
            for (auto& e : a) { e = io.rawInteger(-2, 2); }
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = 0; c < n; ++c)
                {
                    double sum = 0.0;
                    for (int32_t k = 0; k < n; ++k)
                    {
                        sum += a[k + n * r] * a[k + n * c];
                    }
                    at(r, c) = sum + (r == c ? 1.0 : 0.0);
                }
            }
            for (auto& e : lcp.q) { e = io.rawInteger(-5, 5); }
        }
        else if (mode == 3)
        {
            for (auto& e : lcp.q) { e = io.rawInteger(-4, 4); }
        }
        else if (mode == 4)
        {
            for (auto& e : lcp.m) { e = -static_cast<double>(io.rawInteger(0, 3)); }
            for (auto& e : lcp.q) { e = -static_cast<double>(io.rawInteger(1, 5)); }
        }
        else
        {
            for (auto& e : lcp.m) { e = io.rawInteger(-3, 3); }
            for (auto& e : lcp.q) { e = io.rawInteger(0, 5); }
        }
        return lcp;
    }

    void RecordLcp(oracle::Ctx& io, Lcp const& lcp)
    {
        io.given(static_cast<double>(lcp.n));
        io.given(static_cast<double>(lcp.mode));
        for (double e : lcp.q) { io.given(e); }
        for (double e : lcp.m) { io.given(e); }
    }

    void RunLcp(oracle::Ctx& io, Lcp const& lcp, int32_t maxIterations)
    {
        LCPSolver<double> solver(lcp.n);
        if (maxIterations != 0) { solver.SetMaxIterations(maxIterations); }
        io.outInt(solver.GetMaxIterations());
        std::vector<double> w(static_cast<size_t>(lcp.n), 0.0);
        std::vector<double> z(static_cast<size_t>(lcp.n), 0.0);
        LCPSolverShared<double>::Result result{};
        bool success = solver.Solve(lcp.q, lcp.m, w, z, &result);
        io.outBool(success);
        io.outInt(static_cast<int32_t>(result));
        io.outInt(solver.GetNumIterations());
        for (double e : w) { io.outReal(e); }
        for (double e : z) { io.outReal(e); }
    }
}

// Lemke's method over the whole generator: solvable, infeasible, trivially
// solvable and degenerate inputs. Every step is + - * / and the pivot choice
// is a lexicographic comparison of the perturbation polynomials, so the
// pivoting sequence, the iteration count and the solution are bit-identical.
ORACLE_CASE("LCPSolver.solve")
{
    int32_t mode = io.rawInteger(0, 6);
    Lcp lcp = DrawLcpRaw(io, mode);
    RecordLcp(io, lcp);
    // 0 keeps the default n*n budget; the others are explicit budgets, and
    // a negative one restores the default through SetMaxIterations.
    int32_t maxIterations = io.integer(-2, 40);
    RunLcp(io, lcp, maxIterations);
}

// M positive definite: the LCP always has a solution, so the solver must
// report HAS_TRIVIAL_SOLUTION or HAS_NONTRIVIAL_SOLUTION.
ORACLE_CASE("LCPSolver.solve.positiveDefinite")
{
    Lcp lcp = DrawLcpRaw(io, 2);
    RecordLcp(io, lcp);
    RunLcp(io, lcp, 0);
}

// q >= 0 takes the trivial-solution early return before any pivoting.
ORACLE_CASE("LCPSolver.solve.trivial")
{
    Lcp lcp = DrawLcpRaw(io, 5);
    RecordLcp(io, lcp);
    RunLcp(io, lcp, 0);
}

// M <= 0 entrywise with q < 0 is infeasible: w = q + M*z is negative for
// every z >= 0, so the driving variable can never leave the dictionary and
// the solver returns NO_SOLUTION. The magnitudes are ordinary, which keeps
// the case clear of the subnormal-pivot defect of issue #476.
ORACLE_CASE("LCPSolver.solve.noSolution")
{
    Lcp lcp = DrawLcpRaw(io, 4);
    RecordLcp(io, lcp);
    RunLcp(io, lcp, 0);
}

// A budget of one or two iterations forces FAILED_TO_CONVERGE on inputs that
// would otherwise need more pivots.
ORACLE_CASE("LCPSolver.solve.maxIterations")
{
    int32_t mode = io.rawInteger(0, 2);
    Lcp lcp = DrawLcpRaw(io, mode, 3, 6);
    RecordLcp(io, lcp);
    int32_t maxIterations = io.integer(1, 2);
    RunLcp(io, lcp, maxIterations);
}

// UPSTREAM DEFECT found by this oracle (both sides agree bit for bit, and
// both are wrong). The header's own comment anticipates the mechanism -- "it
// is possible that theoretically mAugmented[r][driving] is zero but rounding
// errors cause it to be slightly negative" -- and hopes the outcome is a
// FAILED_TO_CONVERGE. On the small integer LCP below (q < 0, M <= 0
// entrywise, so w = q + M*z < 0 for every z >= 0 and the problem is provably
// infeasible) the outcome is worse: at iteration 6 the ratio test accepts a
// pivot of -1.1102230246251565e-16, the reciprocal is about 1e16, the
// dictionary loses all its digits, and after 19 iterations the artificial
// variable happens to leave the basis, so the solver returns
// HAS_NONTRIVIAL_SOLUTION with z = (16, 4, 2/9, 16, 4, 0) and a residual
// |w - q - M*z| of about 80. This is the round-off failure recorded as
// gtengine-js issue #476 for a SUBNORMAL pivot, reached here with ordinary
// integer data. Found once in 2000 records of LCPSolver.solve.noSolution.
ORACLE_CASE("LCPSolver.solve.roundoffPivot")
{
    // A power-of-two scale keeps every entry exactly representable and
    // reproduces the same pivot sequence.
    double s = PowerOfTwo(io.rawInteger(-8, 8));
    double const q6[6] = { -5, -3, -2, -4, -1, -5 };
    double const m6[36] =
    {
        -1, -2, -2,  0, -3, -3,
        -2,  0, -3, -1, -2, -3,
        -3, -3, -2,  0, -2, -2,
        -2, -3,  0, -1, -2, -2,
        -3,  0, -2, -2,  0, -3,
         0, -3, -2, -2, -3, -2
    };
    Lcp lcp;
    lcp.n = 6;
    lcp.mode = 9;
    lcp.q.resize(6);
    lcp.m.resize(36);
    for (int32_t i = 0; i < 6; ++i) { lcp.q[i] = s * q6[i]; }
    for (int32_t i = 0; i < 36; ++i) { lcp.m[i] = s * m6[i]; }
    RecordLcp(io, lcp);
    int32_t maxIterations = io.integer(0, 40);
    RunLcp(io, lcp, maxIterations);
}

// The INVALID_INPUT return of the dynamic solver: a q or an M with fewer
// elements than the dimension demands. (The other INVALID_INPUT the port
// reports, construction with n <= 0, cannot be compared: upstream leaves its
// member pointers null and Solve dereferences them, which is an access
// violation rather than a catchable exception.)
ORACLE_CASE("LCPSolver.solve.invalidInput")
{
    int32_t n = io.integer(1, 4);
    int32_t shortfall = io.integer(1, 2);
    bool shortenQ = io.boolean();
    size_t qSize = static_cast<size_t>(shortenQ ? std::max(0, n - shortfall) : n);
    size_t mSize = static_cast<size_t>(shortenQ ? n * n
        : std::max(0, n * n - shortfall));
    std::vector<double> q(qSize, 0.0), m(mSize, 0.0);
    for (auto& e : q) { e = io.real(-4.0, 4.0); }
    for (auto& e : m) { e = io.real(-4.0, 4.0); }
    LCPSolver<double> solver(n);
    std::vector<double> w(static_cast<size_t>(n), 0.0);
    std::vector<double> z(static_cast<size_t>(n), 0.0);
    LCPSolverShared<double>::Result result{};
    bool success = solver.Solve(q, m, w, z, &result);
    io.outBool(success);
    io.outInt(static_cast<int32_t>(result));
    for (double e : w) { io.outReal(e); }
    for (double e : z) { io.outReal(e); }
}

// =======================================================================
// FPInterval.h
// =======================================================================
//
// Upstream computes each endpoint under std::fesetround(FE_DOWNWARD) or
// FE_UPWARD. JavaScript has no rounding-mode control, so the port emulates
// directed rounding: it computes the round-to-nearest value, proves
// exactness with TwoSum/TwoProduct where it can, and otherwise steps one
// representable value outward. That is always an enclosure but is up to one
// ulp wider than upstream whenever round-to-nearest happened to round the
// right way.
//
// The exact cases below therefore use dyadic operands on which every
// operation is exactly representable: both sides then produce the same
// endpoints bit for bit, which is what tests the branch structure and the
// port's exactness proofs. The wider behaviour is the subject of the
// declared deviation case at the end.

namespace
{
    using Interval = FPInterval<double>;

    // m * 2^-k with |m| <= 64 and k <= 4: every sum, difference and product
    // of two such values is exactly representable as a double.
    double DrawDyadic(oracle::Ctx& io)
    {
        int32_t m = io.rawInteger(-64, 64);
        int32_t k = io.rawInteger(0, 4);
        return static_cast<double>(m) * PowerOfTwo(-k);
    }

    // +-2^k with |k| <= 8, or exactly zero: division by such a value is
    // exact for a dyadic numerator.
    double DrawPowerOfTwoOrZero(oracle::Ctx& io, int32_t signClass)
    {
        if (signClass == 0) { return 0.0; }
        double magnitude = PowerOfTwo(io.rawInteger(-8, 8));
        return signClass > 0 ? magnitude : -magnitude;
    }

    // An interval whose sign class is 0 (0 <= e0 <= e1), 1 (e0 <= e1 <= 0)
    // or 2 (e0 < 0 < e1). Ordered endpoints are the documented invariant.
    std::array<double, 2> DrawIntervalRaw(oracle::Ctx& io, int32_t signClass)
    {
        double x = std::fabs(DrawDyadic(io));
        double y = std::fabs(DrawDyadic(io));
        double lo = std::min(x, y);
        double hi = std::max(x, y);
        if (signClass == 0) { return { lo, hi }; }
        if (signClass == 1) { return { -hi, -lo }; }
        return { -hi - 1.0, lo + 1.0 };
    }

    // The same three sign classes with power-of-two endpoints, plus the two
    // one-sided classes with a zero endpoint that the division operators
    // single out: 3 is [0, +2^k] and 4 is [-2^k, 0].
    std::array<double, 2> DrawDivisorRaw(oracle::Ctx& io, int32_t signClass)
    {
        double x = PowerOfTwo(io.rawInteger(-8, 8));
        double y = PowerOfTwo(io.rawInteger(-8, 8));
        double lo = std::min(x, y);
        double hi = std::max(x, y);
        if (signClass == 0) { return { lo, hi }; }
        if (signClass == 1) { return { -hi, -lo }; }
        if (signClass == 2) { return { -hi, hi }; }
        if (signClass == 3) { return { 0.0, hi }; }
        return { -hi, 0.0 };
    }

    void OutInterval(oracle::Ctx& io, Interval const& w)
    {
        io.outReal(w[0]);
        io.outReal(w[1]);
    }
}

// The leaf-node operations Add/Sub/Mul/Div(u, v) on two raw variables,
// including Div's division-by-zero return of the whole real line.
ORACLE_CASE("FPInterval.leaf")
{
    double u = DrawDyadic(io);
    io.given(u);
    int32_t divisorClass = io.integer(-1, 1);
    double v = DrawPowerOfTwoOrZero(io, divisorClass);
    io.given(v);
    OutInterval(io, Interval::Add(u, v));
    OutInterval(io, Interval::Sub(u, v));
    OutInterval(io, Interval::Mul(u, v));
    OutInterval(io, Interval::Div(u, v));
}

// The interior-node helpers: the four-argument Add/Sub/Mul, Mul2,
// Reciprocal, ReciprocalDown, ReciprocalUp and Reals. ReciprocalDown and
// ReciprocalUp deliberately produce an infinite endpoint, which is exactly
// where the port must not reuse the nextafter-based widening of SWInterval.
ORACLE_CASE("FPInterval.internal")
{
    int32_t uClass = io.integer(0, 2);
    std::array<double, 2> u = DrawIntervalRaw(io, uClass);
    io.given(u[0]);
    io.given(u[1]);
    // Classes 3 and 4 have a zero endpoint, so Div, Reciprocal,
    // ReciprocalDown and ReciprocalUp evaluate 1/+-0 and produce infinite
    // endpoints. Those are exact results, not rounded ones, in every
    // rounding mode.
    int32_t vClass = io.integer(0, 4);
    std::array<double, 2> v = DrawDivisorRaw(io, vClass);
    io.given(v[0]);
    io.given(v[1]);

    OutInterval(io, Interval::Add(u[0], u[1], v[0], v[1]));
    OutInterval(io, Interval::Sub(u[0], u[1], v[0], v[1]));
    OutInterval(io, Interval::Mul(u[0], u[1], v[0], v[1]));
    OutInterval(io, Interval::Mul2(u[0], u[1], v[0], v[1]));
    OutInterval(io, Interval::Div(u[0], u[1], v[0], v[1]));
    OutInterval(io, Interval::Reciprocal(v[0], v[1]));
    OutInterval(io, Interval::ReciprocalDown(v[1]));
    OutInterval(io, Interval::ReciprocalUp(v[0]));
    OutInterval(io, Interval::Reals());
    Interval fromEndpoints(v);
    OutInterval(io, fromEndpoints);
    Interval degenerate(u[0]);
    OutInterval(io, degenerate);
    Interval defaultConstructed{};
    OutInterval(io, defaultConstructed);
}

// The free operators. The nine sign-class combinations of operator*(u, v)
// and the five of operator/(u, v) are all reached, and the scalar overloads
// are covered on both sides of the operand.
ORACLE_CASE("FPInterval.operators")
{
    int32_t uClass = io.integer(0, 2);
    std::array<double, 2> uRaw = DrawIntervalRaw(io, uClass);
    io.given(uRaw[0]);
    io.given(uRaw[1]);
    int32_t vClass = io.integer(0, 2);
    std::array<double, 2> vRaw = DrawIntervalRaw(io, vClass);
    io.given(vRaw[0]);
    io.given(vRaw[1]);
    int32_t divisorClass = io.integer(0, 4);
    std::array<double, 2> dRaw = DrawDivisorRaw(io, divisorClass);
    io.given(dRaw[0]);
    io.given(dRaw[1]);
    int32_t scalarClass = io.integer(-1, 1);
    double scalar = DrawPowerOfTwoOrZero(io, scalarClass);
    io.given(scalar);

    Interval u(uRaw);
    Interval v(vRaw);
    Interval d(dRaw);

    OutInterval(io, +u);
    OutInterval(io, -u);
    OutInterval(io, u + v);
    OutInterval(io, u + scalar);
    OutInterval(io, scalar + u);
    OutInterval(io, u - v);
    OutInterval(io, u - scalar);
    OutInterval(io, scalar - u);
    OutInterval(io, u * v);
    OutInterval(io, u * scalar);
    OutInterval(io, scalar * u);
    OutInterval(io, u / d);
    OutInterval(io, u / scalar);
    OutInterval(io, scalar / d);
}

// ProductLowerBound and ProductUpperBound. The header requires the caller to
// have selected the rounding mode, which is what the fesetround calls below
// do. Both functions return the wrong bound when u and v both straddle zero
// (gtengine-js issue #75); the port preserves that deliberately, so the two
// sides are compared bit for bit and that branch is exercised on purpose.
ORACLE_CASE("FPInterval.productBounds")
{
    int32_t uClass = io.integer(0, 2);
    std::array<double, 2> u = DrawIntervalRaw(io, uClass);
    io.given(u[0]);
    io.given(u[1]);
    int32_t vClass = io.integer(0, 2);
    std::array<double, 2> v = DrawIntervalRaw(io, vClass);
    io.given(v[0]);
    io.given(v[1]);

    auto saveMode = std::fegetround();
    std::fesetround(FE_DOWNWARD);
    double lower = Interval::ProductLowerBound(u, v);
    std::fesetround(FE_UPWARD);
    double upper = Interval::ProductUpperBound(u, v);
    std::fesetround(saveMode);
    io.outReal(lower);
    io.outReal(upper);
}

// Deviation: the port emulates directed rounding, upstream has it in
// hardware. On operands where the operation is not exactly representable the
// port widens by one ulp in both directions, while upstream keeps the
// round-to-nearest value on whichever side it rounded toward. The port's
// interval therefore contains upstream's, and the endpoints differ by at
// most one ulp; see the PORT DEVIATION note at the top of src/FPInterval.ts.
// The operands are ordinary reals, so almost every operation is inexact.
ORACLE_CASE("FPInterval.directedRounding")
{
    double u0 = io.raw(-8.0, 8.0);
    double u1 = u0 + io.raw(0.0, 8.0);
    io.given(u0);
    io.given(u1);
    double v0 = io.raw(0.25, 8.0);
    double v1 = v0 + io.raw(0.0, 8.0);
    io.given(v0);
    io.given(v1);
    Interval u(u0, u1);
    Interval v(v0, v1);
    OutInterval(io, Interval::Add(u0, v0));
    OutInterval(io, Interval::Sub(u0, v0));
    OutInterval(io, Interval::Mul(u0, v0));
    OutInterval(io, Interval::Div(u0, v0));
    OutInterval(io, u + v);
    OutInterval(io, u - v);
    OutInterval(io, u * v);
    OutInterval(io, u / v);
    OutInterval(io, Interval::Reciprocal(v0, v1));
    // 1/+0 is +infinity in every rounding mode; the port's emulated
    // round-down steps that to MAX_VALUE, which is still a valid lower
    // bound but not upstream's endpoint.
    OutInterval(io, Interval::Reciprocal(-v1, 0.0));
    OutInterval(io, Interval::ReciprocalDown(0.0));
    OutInterval(io, Interval::ReciprocalUp(0.0));
}
