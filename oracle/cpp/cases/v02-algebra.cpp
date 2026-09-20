// Verify group 2 (algebra): differential cases for ConvertCoordinates.h,
// Matrix2x2.h, Matrix3x3.h, Matrix4x4.h, Quaternion.h, Projection.h,
// Rotation.h, RotationEstimate.h and Transform.h.
//
// Neither GTE_USE_COL_MAJOR nor GTE_USE_VEC_MAT is defined, so the row-major /
// GTE_USE_MAT_VEC default that the port implements is the configuration
// measured.
//
// Overload resolution. This translation unit includes Matrix2x2.h, Matrix3x3.h
// and Matrix4x4.h (directly and through Transform.h and Projection.h), so
// every unqualified Inverse/Adjoint/Determinant/Trace of a 2x2, 3x3 or 4x4
// matrix resolves to the closed forms of those headers and NOT to the
// Gaussian-elimination templates of Matrix.h, which group 1 measures. That
// includes the calls inside Transform.h (GetHInverse, Inverse) and inside
// Hyperellipsoid::FromCoefficients, which Projection.h reaches.
//
// libm. Many entry points here call sin, cos, acos, asin or atan2. Cases whose
// compared path never touches libm are declared { exact: true }; the others
// carry a tolerance and say which function is responsible. Where a conversion
// mixes the two (an axis computed with sqrt alone, an angle computed with
// acos), the arithmetic-only outputs are emitted with the harness's
// outRealExact so they are still required to be bit-identical.
//
// Generator modes. Every mode of a given helper records exactly the same
// number of doubles, so the TypeScript replay reads a fixed layout and never
// needs to know which mode produced a record:
//   0  small integer lattice: exact arithmetic, exact zeros and ties,
//      exactly singular matrices
//   1  uniform
//   2  a constructed exactly singular matrix (last row an integer combination
//      of the others)
//   3  wild magnitudes: signed zeros, subnormals on both sides of the
//      1/x-overflow threshold, values near the overflow limit
// Rotation-valued inputs have their own modes (coordinate axes, signed
// permutation matrices, angles at 0 and +-pi/2, unit quaternions with w of
// either sign), listed at the helpers.
#define ORACLE_FAMILY "v02-algebra"
#include "Oracle.h"

#include <Mathematics/ConvertCoordinates.h>
#include <Mathematics/Matrix2x2.h>
#include <Mathematics/Matrix3x3.h>
#include <Mathematics/Matrix4x4.h>
#include <Mathematics/Projection.h>
#include <Mathematics/Quaternion.h>
#include <Mathematics/Rotation.h>
#include <Mathematics/RotationEstimate.h>
#include <Mathematics/Transform.h>

#include <array>
#include <cmath>
#include <cstdint>

using namespace gte;

namespace
{
    // ---- scalar populations -------------------------------------------

    // A value drawn from the "wild" population. Only unrecorded raw draws are
    // used here; the caller records the single resulting double.
    double WildRaw(oracle::Ctx& io)
    {
        double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        int32_t k = io.rawInteger(0, 9);
        double mantissa = io.raw(1.0, 2.0);
        switch (k)
        {
        case 0:  return sign * 0.0;                          // +-0
        case 1:  return sign * std::ldexp(1.0, -1074);       // smallest subnormal
        case 2:  return sign * std::ldexp(mantissa, -1060);  // 1/x overflows
        case 3:  return sign * std::ldexp(mantissa, -1023);  // 1/x is finite
        case 4:  return sign * std::ldexp(mantissa, -1022);  // smallest normal
        case 5:  return sign * std::ldexp(mantissa, 500);
        case 6:  return sign * std::ldexp(mantissa, -500);
        case 7:  return sign * static_cast<double>(io.rawInteger(0, 4));
        case 8:  return sign * std::ldexp(mantissa, 52);     // integer-valued
        default: return sign * io.raw(0.0, 4.0);
        }
    }

    // One recorded scalar in the requested mode.
    double Scalar(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 0)
        {
            return io.lattice(-3, 3);
        }
        if (mode == 1)
        {
            return io.real(-10.0, 10.0);
        }
        if (mode == 2)
        {
            return io.lattice(-4, 4);
        }
        return io.given(WildRaw(io));
    }

    // ---- matrix helpers -----------------------------------------------

    // Record an R-by-C matrix as R*C inputs, row-major.
    template <int32_t R, int32_t C>
    Matrix<R, C, double> GivenMat(oracle::Ctx& io, Matrix<R, C, double> const& m)
    {
        for (int32_t r = 0; r < R; ++r)
        {
            for (int32_t c = 0; c < C; ++c)
            {
                io.given(m(r, c));
            }
        }
        return m;
    }

    // A square matrix in the requested mode; N*N recorded doubles either way.
    // Mode 2 is an exactly singular matrix: the last row is an integer
    // combination of the preceding rows of a small integer matrix.
    template <int32_t N>
    Matrix<N, N, double> SquareMat(oracle::Ctx& io, int32_t mode)
    {
        Matrix<N, N, double> m{};
        if (mode == 2)
        {
            for (int32_t r = 0; r + 1 < N; ++r)
            {
                for (int32_t c = 0; c < N; ++c)
                {
                    m(r, c) = static_cast<double>(io.rawInteger(-3, 3));
                }
            }
            std::array<double, N> w{};
            for (int32_t r = 0; r + 1 < N; ++r)
            {
                w[r] = static_cast<double>(io.rawInteger(-2, 2));
            }
            for (int32_t c = 0; c < N; ++c)
            {
                double sum = 0.0;
                for (int32_t r = 0; r + 1 < N; ++r)
                {
                    sum += w[r] * m(r, c);
                }
                m(N - 1, c) = sum;
            }
            return GivenMat(io, m);
        }

        for (int32_t r = 0; r < N; ++r)
        {
            for (int32_t c = 0; c < N; ++c)
            {
                m(r, c) = Scalar(io, mode);
            }
        }
        return m;
    }

    // A vector in the requested mode; N recorded doubles either way.
    template <int32_t N>
    Vector<N, double> Vec(oracle::Ctx& io, int32_t mode)
    {
        Vector<N, double> v{};
        for (int32_t i = 0; i < N; ++i)
        {
            v[i] = Scalar(io, mode);
        }
        return v;
    }
}

// ===================================================================
// Matrix2x2.h
// ===================================================================

// Inverse (both the matrix and the invertibility flag), Adjoint, Determinant
// and Trace of a 2x2 matrix. Arithmetic only.
ORACLE_CASE("Matrix2x2.inverse")
{
    Matrix2x2<double> M = SquareMat<2>(io, io.index() % 4);
    bool invertible = false;
    Matrix2x2<double> inv = Inverse(M, &invertible);
    io.outBool(invertible);
    io.outMat(inv);
    io.outMat(Adjoint(M));
    io.outReal(Determinant(M));
    io.outReal(Trace(M));
    // The no-flag overload must agree with the flag overload.
    io.outMat(Inverse(M));
}

// MakeRotation and GetRotationAngle. MakeRotation calls cos and sin,
// GetRotationAngle calls atan2, so this case carries the default tolerance.
ORACLE_CASE("Matrix2x2.rotation")
{
    // Mode 0 pins the angle at values where cos and sin are exactly
    // representable on both runtimes (0 and +-pi/2 from the double nearest
    // pi/2 are not exact, so the mode exists to reach the axis-aligned
    // rotations, not to avoid libm).
    double angle = 0.0;
    if (io.index() % 4 == 0)
    {
        int32_t k = io.rawInteger(-2, 2);
        angle = io.given(static_cast<double>(k) * GTE_C_HALF_PI);
    }
    else
    {
        angle = io.real(-7.0, 7.0);
    }
    Matrix2x2<double> rotation{};
    MakeRotation(angle, rotation);
    io.outMat(rotation);
    io.outReal(GetRotationAngle(rotation));

    // GetRotationAngle of an arbitrary recorded matrix (the caller is
    // responsible for it being a rotation; atan2 is defined regardless).
    Matrix2x2<double> other = SquareMat<2>(io, io.index() % 4);
    io.outReal(GetRotationAngle(other));
}

// DoTransform(M,V), DoTransform(A,B), SetBasis and GetBasis.
ORACLE_CASE("Matrix2x2.doTransform")
{
    int32_t mode = io.index() % 4;
    Matrix2x2<double> A = SquareMat<2>(io, mode);
    Matrix2x2<double> B = SquareMat<2>(io, mode);
    Vector2<double> V = Vec<2>(io, mode);
    Vector2<double> W = Vec<2>(io, mode);
    int32_t i = io.integer(0, 1);
    io.outVec(DoTransform(A, V));
    io.outMat(DoTransform(A, B));
    io.outVec(GetBasis(A, i));
    Matrix2x2<double> C = A;
    SetBasis(C, i, W);
    io.outMat(C);
}

// ===================================================================
// Matrix3x3.h
// ===================================================================

ORACLE_CASE("Matrix3x3.inverse")
{
    Matrix3x3<double> M = SquareMat<3>(io, io.index() % 4);
    bool invertible = false;
    Matrix3x3<double> inv = Inverse(M, &invertible);
    io.outBool(invertible);
    io.outMat(inv);
    io.outMat(Adjoint(M));
    io.outReal(Determinant(M));
    io.outReal(Trace(M));
    io.outMat(Inverse(M));
}

ORACLE_CASE("Matrix3x3.doTransform")
{
    int32_t mode = io.index() % 4;
    Matrix3x3<double> A = SquareMat<3>(io, mode);
    Matrix3x3<double> B = SquareMat<3>(io, mode);
    Vector3<double> V = Vec<3>(io, mode);
    Vector3<double> W = Vec<3>(io, mode);
    int32_t i = io.integer(0, 2);
    io.outVec(DoTransform(A, V));
    io.outMat(DoTransform(A, B));
    io.outVec(GetBasis(A, i));
    Matrix3x3<double> C = A;
    SetBasis(C, i, W);
    io.outMat(C);
}

// ===================================================================
// Matrix4x4.h
// ===================================================================

ORACLE_CASE("Matrix4x4.inverse")
{
    Matrix4x4<double> M = SquareMat<4>(io, io.index() % 4);
    bool invertible = false;
    Matrix4x4<double> inv = Inverse(M, &invertible);
    io.outBool(invertible);
    io.outMat(inv);
    io.outMat(Adjoint(M));
    io.outReal(Determinant(M));
    io.outReal(Trace(M));
    io.outMat(Inverse(M));
}

ORACLE_CASE("Matrix4x4.doTransform")
{
    int32_t mode = io.index() % 4;
    Matrix4x4<double> A = SquareMat<4>(io, mode);
    Matrix4x4<double> B = SquareMat<4>(io, mode);
    Vector4<double> V = Vec<4>(io, mode);
    Vector4<double> W = Vec<4>(io, mode);
    int32_t i = io.integer(0, 3);
    io.outVec(DoTransform(A, V));
    io.outMat(DoTransform(A, B));
    io.outVec(GetBasis(A, i));
    Matrix4x4<double> C = A;
    SetBasis(C, i, W);
    io.outMat(C);
}

namespace
{
    // ---- quaternion and rotation helpers ------------------------------
    //
    // These build their values from unrecorded raw draws so that a case can
    // reject a candidate before recording it. The caller records the accepted
    // value with GivenQuat / GivenMat, which is also the only thing the
    // TypeScript replay reads.

    Quaternion<double> GivenQuat(oracle::Ctx& io, Quaternion<double> const& q)
    {
        for (int32_t i = 0; i < 4; ++i)
        {
            io.given(q[i]);
        }
        return q;
    }

    // A unit-length quaternion, unrecorded. The modes are
    //   0  exact: the eight axis-aligned quaternions (+-1 in one slot) and the
    //      sixteen (+-1/2,+-1/2,+-1/2,+-1/2); all exactly unit length, and
    //      they make the matrix->quaternion trace split land on exact ties
    //   1  a normalized uniform 4-tuple (w of either sign)
    //   2  the same with w forced negative
    //   3  near-identity (angle ~ 2^-e) or near-pi (angle ~ pi - 2^-e) for
    //      e in [1,30]
    Quaternion<double> UnitQuatRaw(oracle::Ctx& io, int32_t mode)
    {
        Quaternion<double> q{};
        if (mode == 0)
        {
            int32_t k = io.rawInteger(0, 23);
            if (k < 8)
            {
                q[k % 4] = (k < 4 ? 1.0 : -1.0);
            }
            else
            {
                int32_t bits = k - 8;
                for (int32_t i = 0; i < 4; ++i)
                {
                    q[i] = (((bits >> i) & 1) != 0 ? -0.5 : 0.5);
                }
            }
            return q;
        }

        if (mode == 3)
        {
            Vector3<double> u{};
            double len = 0.0;
            do
            {
                for (int32_t i = 0; i < 3; ++i) { u[i] = io.raw(-1.0, 1.0); }
                len = Length(u);
            }
            while (len < 0.1 || len > 1.0);
            Normalize(u);
            double s = std::ldexp(1.0, -io.rawInteger(1, 30));
            bool nearPi = (io.rawInteger(0, 1) != 0);
            Quaternion<double> p{};
            if (nearPi)
            {
                p = Quaternion<double>(u[0], u[1], u[2], s);
            }
            else
            {
                p = Quaternion<double>(s * u[0], s * u[1], s * u[2], 1.0);
            }
            Normalize(p);
            return p;
        }

        double len = 0.0;
        do
        {
            for (int32_t i = 0; i < 4; ++i) { q[i] = io.raw(-1.0, 1.0); }
            len = std::sqrt(Dot(q, q));
        }
        while (len < 0.1 || len > 1.0);
        Normalize(q);
        if (mode == 2 && q[3] > 0.0)
        {
            for (int32_t i = 0; i < 4; ++i) { q[i] = -q[i]; }
        }
        return q;
    }

    // A quaternion in the requested mode, unrecorded. Mode 0 is a small
    // integer lattice (including the zero quaternion, which Inverse and
    // Normalize special-case); the other modes are the unit-length
    // populations.
    Quaternion<double> QuatRaw(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 0)
        {
            Quaternion<double> q{};
            for (int32_t i = 0; i < 4; ++i)
            {
                q[i] = static_cast<double>(io.rawInteger(-2, 2));
            }
            return q;
        }
        return UnitQuatRaw(io, mode);
    }

    // A 3x3 rotation matrix written into the upper-left block of an N-by-N
    // identity, unrecorded. Mode 0 is a signed permutation matrix of
    // determinant +1 (entries exactly 0 and +-1, which is where the Euler
    // gimbal-lock tests r(i,j) == +-1 fire); the other modes are the
    // quaternion populations converted to a matrix.
    template <int32_t N>
    Matrix<N, N, double> RotMatRaw(oracle::Ctx& io, int32_t mode)
    {
        Matrix<N, N, double> r{};
        r.MakeIdentity();
        if (mode == 0)
        {
            // The 24 rotations of the octahedral group: a permutation of the
            // three axes times a sign pattern whose product matches the
            // permutation's parity, so the determinant is +1.
            static int32_t const perm[6][3] =
            { {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0} };
            static int32_t const permParity[6] = { 0, 1, 1, 0, 0, 1 };
            int32_t p = io.rawInteger(0, 5);
            int32_t s = io.rawInteger(0, 3);
            int32_t sign0 = ((s & 1) != 0 ? -1 : 1);
            int32_t sign1 = ((s & 2) != 0 ? -1 : 1);
            int32_t sign2 = (sign0 * sign1) * (permParity[p] != 0 ? -1 : 1);
            int32_t sgn[3] = { sign0, sign1, sign2 };
            for (int32_t c = 0; c < 3; ++c)
            {
                for (int32_t rr = 0; rr < 3; ++rr) { r(rr, c) = 0.0; }
                r(perm[p][c], c) = static_cast<double>(sgn[c]);
            }
            return r;
        }

        Quaternion<double> q = UnitQuatRaw(io, mode);
        Rotation<N, double> rotation(q);
        r = rotation;
        return r;
    }

    // The axis triple of an Euler-angle factorization. 'which' in [0,11]
    // selects one of the twelve valid orders, six with three distinct axes
    // and six with the first and last repeated. 'which' == 12 is an invalid
    // order (axis[1] == axis[0]), which the conversions report as INVALID.
    void EulerAxes(int32_t which, int32_t& i0, int32_t& i1, int32_t& i2)
    {
        static int32_t const axes[13][3] =
        {
            {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0},
            {0,1,0}, {0,2,0}, {1,0,1}, {1,2,1}, {2,0,2}, {2,1,2},
            {1,1,2}
        };
        i0 = axes[which][0];
        i1 = axes[which][1];
        i2 = axes[which][2];
    }
}

// ===================================================================
// Quaternion.h
// ===================================================================

// Unary operators, the linear-algebraic operators, Dot, Length, Normalize,
// the six comparisons and the five special quaternions. Arithmetic only
// (Length and Normalize use sqrt).
ORACLE_CASE("Quaternion.arithmetic")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q0 = GivenQuat(io, QuatRaw(io, mode));
    Quaternion<double> q1 = GivenQuat(io, QuatRaw(io, mode));
    double s = Scalar(io, mode);

    Quaternion<double> negated = -q0;
    io.outVec(Vector<4, double>{ negated[0], negated[1], negated[2], negated[3] });
    Quaternion<double> plus = +q0;
    io.outVec(Vector<4, double>{ plus[0], plus[1], plus[2], plus[3] });
    Quaternion<double> sum = q0 + q1;
    io.outVec(Vector<4, double>{ sum[0], sum[1], sum[2], sum[3] });
    Quaternion<double> dif = q0 - q1;
    io.outVec(Vector<4, double>{ dif[0], dif[1], dif[2], dif[3] });
    Quaternion<double> qs = q0 * s;
    io.outVec(Vector<4, double>{ qs[0], qs[1], qs[2], qs[3] });
    Quaternion<double> sq = s * q0;
    io.outVec(Vector<4, double>{ sq[0], sq[1], sq[2], sq[3] });
    Quaternion<double> qd = q0 / s;
    io.outVec(Vector<4, double>{ qd[0], qd[1], qd[2], qd[3] });

    io.outReal(Dot(q0, q1));
    io.outReal(Length(q0));
    Quaternion<double> unit = q0;
    double length = Normalize(unit);
    io.outReal(length);
    io.outVec(Vector<4, double>{ unit[0], unit[1], unit[2], unit[3] });

    io.outBool(q0 == q1);
    io.outBool(q0 != q1);
    io.outBool(q0 < q1);
    io.outBool(q0 <= q1);
    io.outBool(q0 > q1);
    io.outBool(q0 >= q1);

    Quaternion<double> zero = Quaternion<double>::Zero();
    io.outVec(Vector<4, double>{ zero[0], zero[1], zero[2], zero[3] });
    Quaternion<double> qi = Quaternion<double>::I();
    io.outVec(Vector<4, double>{ qi[0], qi[1], qi[2], qi[3] });
    Quaternion<double> qj = Quaternion<double>::J();
    io.outVec(Vector<4, double>{ qj[0], qj[1], qj[2], qj[3] });
    Quaternion<double> qk = Quaternion<double>::K();
    io.outVec(Vector<4, double>{ qk[0], qk[1], qk[2], qk[3] });
    Quaternion<double> qid = Quaternion<double>::Identity();
    io.outVec(Vector<4, double>{ qid[0], qid[1], qid[2], qid[3] });
}

// The Hamilton product (noncommutative), Conjugate and Inverse, including the
// zero quaternion that Inverse special-cases.
ORACLE_CASE("Quaternion.multiply")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q0 = GivenQuat(io, QuatRaw(io, mode));
    Quaternion<double> q1 = GivenQuat(io, QuatRaw(io, mode));

    Quaternion<double> p01 = q0 * q1;
    io.outVec(Vector<4, double>{ p01[0], p01[1], p01[2], p01[3] });
    Quaternion<double> p10 = q1 * q0;
    io.outVec(Vector<4, double>{ p10[0], p10[1], p10[2], p10[3] });
    Quaternion<double> conj = Conjugate(q0);
    io.outVec(Vector<4, double>{ conj[0], conj[1], conj[2], conj[3] });
    Quaternion<double> inv = Inverse(q0);
    io.outVec(Vector<4, double>{ inv[0], inv[1], inv[2], inv[3] });
    Quaternion<double> invZero = Inverse(Quaternion<double>::Zero());
    io.outVec(Vector<4, double>{ invZero[0], invZero[1], invZero[2], invZero[3] });
}

// Rotate(q, Vector<3>). Arithmetic only.
ORACLE_CASE("Quaternion.rotate.3d")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q = GivenQuat(io, QuatRaw(io, mode));
    Vector3<double> u = Vec<3>(io, mode);
    io.outVec(Rotate(q, u));
}

// Rotate(q, Vector<4>), the homogeneous form. The 4-tuple has a zero last
// component on half the records (the documented use) and an arbitrary one on
// the rest.
ORACLE_CASE("Quaternion.rotate.4d")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q = GivenQuat(io, QuatRaw(io, mode));
    Vector4<double> u{};
    for (int32_t i = 0; i < 3; ++i)
    {
        u[i] = Scalar(io, mode);
    }
    bool affine = (io.index() % 2 == 0);
    u[3] = (affine ? io.given(0.0) : Scalar(io, mode));
    io.outVec(Rotate(q, u));
}

// Quaternion::Slerp. ChebyshevRatiosUsingCosAngle calls acos and sin unless
// cosA >= 1 (the exact arm, which mode 0 reaches by drawing q1 == q0), so the
// case carries the default tolerance.
ORACLE_CASE("Quaternion.slerp")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q0raw = UnitQuatRaw(io, mode);
    Quaternion<double> q1raw = (io.rawInteger(0, 3) == 0 ? q0raw : UnitQuatRaw(io, mode));
    Quaternion<double> q0 = GivenQuat(io, q0raw);
    Quaternion<double> q1 = GivenQuat(io, q1raw);
    double t = io.real(0.0, 1.0);
    Quaternion<double> q = Slerp(t, q0, q1);
    io.outVec(Vector<4, double>{ q[0], q[1], q[2], q[3] });
}

// Quaternion::SlerpR. The documented preprocessing makes Dot(q0,q1) >= 0, so
// q1 is negated here when the dot product is negative; the recorded q1 is the
// preprocessed one. acos and sin unless the dot product is >= 1.
ORACLE_CASE("Quaternion.slerpR")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q0raw = UnitQuatRaw(io, mode);
    Quaternion<double> q1raw = (io.rawInteger(0, 3) == 0 ? q0raw : UnitQuatRaw(io, mode));
    if (Dot(q0raw, q1raw) < 0.0)
    {
        for (int32_t i = 0; i < 4; ++i) { q1raw[i] = -q1raw[i]; }
    }
    Quaternion<double> q0 = GivenQuat(io, q0raw);
    Quaternion<double> q1 = GivenQuat(io, q1raw);
    double t = io.real(0.0, 1.0);
    Quaternion<double> q = SlerpR(t, q0, q1);
    io.outVec(Vector<4, double>{ q[0], q[1], q[2], q[3] });
}

// Quaternion::SlerpRP, the preprocessed form that takes cosA as an argument.
// cosA is the preprocessed Dot(q0,q1) and is recorded as an input.
ORACLE_CASE("Quaternion.slerpRP")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q0raw = UnitQuatRaw(io, mode);
    Quaternion<double> q1raw = (io.rawInteger(0, 3) == 0 ? q0raw : UnitQuatRaw(io, mode));
    if (Dot(q0raw, q1raw) < 0.0)
    {
        for (int32_t i = 0; i < 4; ++i) { q1raw[i] = -q1raw[i]; }
    }
    Quaternion<double> q0 = GivenQuat(io, q0raw);
    Quaternion<double> q1 = GivenQuat(io, q1raw);
    double cosA = io.given(Dot(q0raw, q1raw));
    double t = io.real(0.0, 1.0);
    Quaternion<double> q = SlerpRP(t, q0, q1, cosA);
    io.outVec(Vector<4, double>{ q[0], q[1], q[2], q[3] });
}

// Quaternion::SlerpRPH, the half-angle form. The midpoint qh and cosAH come
// from the documented preprocessing and are recorded as inputs, so the only
// libm on the compared path is inside ChebyshevRatiosUsingCosAngle. The
// parameter t is drawn across 1/2 so that both arms of 'twoT <= 1' run.
ORACLE_CASE("Quaternion.slerpRPH")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q0raw = UnitQuatRaw(io, mode);
    Quaternion<double> q1raw = (io.rawInteger(0, 3) == 0 ? q0raw : UnitQuatRaw(io, mode));
    if (Dot(q0raw, q1raw) < 0.0)
    {
        for (int32_t i = 0; i < 4; ++i) { q1raw[i] = -q1raw[i]; }
    }
    double cosAraw = Dot(q0raw, q1raw);
    double cosAHraw = std::sqrt((1.0 + cosAraw) / 2.0);
    Quaternion<double> qhraw = (q0raw + q1raw) / (2.0 * cosAHraw);
    Quaternion<double> q0 = GivenQuat(io, q0raw);
    Quaternion<double> q1 = GivenQuat(io, q1raw);
    Quaternion<double> qh = GivenQuat(io, qhraw);
    double cosAH = io.given(cosAHraw);
    double t = io.real(0.0, 1.0);
    Quaternion<double> q = SlerpRPH(t, q0, q1, qh, cosAH);
    io.outVec(Vector<4, double>{ q[0], q[1], q[2], q[3] });
}

// Throw parity: ChebyshevRatiosUsingCosAngle calls LogError when the angle is
// pi, which SlerpRP reaches for antipodal quaternions (cosA == -1). Half the
// records use cosA == -1 exactly and half a valid cosA, so the case measures
// both arms of the guard.
ORACLE_CASE("Quaternion.slerpRP.invalidAngle")
{
    Quaternion<double> q0 = GivenQuat(io, UnitQuatRaw(io, 1));
    Quaternion<double> q1raw = q0;
    for (int32_t i = 0; i < 4; ++i) { q1raw[i] = -q1raw[i]; }
    Quaternion<double> q1 = GivenQuat(io, q1raw);
    double cosA = io.given(io.index() % 2 == 0 ? -1.0 : 0.25);
    double t = io.real(0.0, 1.0);
    Quaternion<double> q = SlerpRP(t, q0, q1, cosA);
    io.outVec(Vector<4, double>{ q[0], q[1], q[2], q[3] });
}

namespace
{
    // Every Rotation case records its input in 3D and runs the conversion for
    // both N = 3 and N = 4, so both instantiations are measured on every
    // record with a single fixed input layout. A Rotation<4> matrix is the
    // 3x3 rotation in the upper-left block of the identity, which is exactly
    // what the conversions produce and consume.
    Matrix4x4<double> Lift3To4(Matrix3x3<double> const& r3)
    {
        Matrix4x4<double> r4{};
        r4.MakeIdentity();
        for (int32_t r = 0; r < 3; ++r)
        {
            for (int32_t c = 0; c < 3; ++c)
            {
                r4(r, c) = r3(r, c);
            }
        }
        return r4;
    }

    // Emit an axis-angle pair: the axis is arithmetic only (a difference of
    // matrix entries and a Normalize, or a quaternion component divided by a
    // sqrt), so it is required to be bit-identical; the angle comes from acos.
    template <int32_t N>
    void OutAxisAngle(oracle::Ctx& io, AxisAngle<N, double> const& a)
    {
        for (int32_t i = 0; i < N; ++i)
        {
            io.outReal(a.axis[i]);
        }
        io.outReal(a.angle);
    }

    void OutEulerAngles(oracle::Ctx& io, EulerAngles<double> const& e)
    {
        for (int32_t i = 0; i < 3; ++i) { io.outInt(e.axis[i]); }
        for (int32_t i = 0; i < 3; ++i) { io.outReal(e.angle[i]); }
        io.outInt(static_cast<int32_t>(e.result));
    }

    // The port validates the axis that upstream's 'angle in (0,pi)' arm
    // extracts from R - Transpose(R) and falls back to the symmetric formula
    // when the extraction underflowed (issue #374). This is that predicate,
    // written with upstream's own expressions; cases that must stay on inputs
    // where upstream is sound reject the candidates it rejects.
    bool AxisExtractionIsSound(Matrix3x3<double> const& r)
    {
        double trace = r(0, 0) + r(1, 1) + r(2, 2);
        double cs = 0.5 * (trace - 1.0);
        cs = std::max(std::min(cs, 1.0), -1.0);
        double angle = std::acos(cs);
        if (angle <= 0.0 || angle >= GTE_C_PI)
        {
            // The identity and the symmetric arms do not extract an axis from
            // the antisymmetric part, so they are always sound.
            return true;
        }
        Vector3<double> axis{};
        axis[0] = r(2, 1) - r(1, 2);
        axis[1] = r(0, 2) - r(2, 0);
        axis[2] = r(1, 0) - r(0, 1);
        Normalize(axis);
        return std::fabs(Dot(axis, axis) - 1.0) <= 1e-12;
    }

    // A rotation matrix on which upstream's axis extraction is sound. The
    // rejection loop is capped; the fallback is a coordinate-axis rotation by
    // pi/2, which takes the well-conditioned arm.
    Matrix3x3<double> SoundRotMatRaw(oracle::Ctx& io, int32_t mode)
    {
        Matrix3x3<double> r{};
        for (int32_t attempt = 0; attempt < 32; ++attempt)
        {
            r = RotMatRaw<3>(io, mode);
            if (AxisExtractionIsSound(r))
            {
                return r;
            }
        }
        r.MakeIdentity();
        r(1, 1) = 0.0;
        r(1, 2) = -1.0;
        r(2, 1) = 1.0;
        r(2, 2) = 0.0;
        return r;
    }
}

// ===================================================================
// Rotation.h
// ===================================================================

// Quaternion -> matrix, for N = 3 and N = 4. Arithmetic only.
ORACLE_CASE("Rotation.quaternionToMatrix")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q = GivenQuat(io, UnitQuatRaw(io, mode));
    Rotation<3, double> r3(q);
    Matrix3x3<double> m3 = r3;
    io.outMat(m3);
    Rotation<4, double> r4(q);
    Matrix4x4<double> m4 = r4;
    io.outMat(m4);
}

// Matrix -> quaternion, for N = 3 and N = 4. Arithmetic only (sqrt). The
// four arms of the trace-based case split are r22 <= 0 crossed with
// dif10 <= 0 / sum10 <= 0; mode 0 (signed permutation matrices) evaluates
// those tests at exact ties.
ORACLE_CASE("Rotation.matrixToQuaternion")
{
    int32_t mode = io.index() % 4;
    Matrix3x3<double> m3 = GivenMat(io, RotMatRaw<3>(io, mode));
    Rotation<3, double> r3(m3);
    Quaternion<double> q3 = r3;
    io.outVec(Vector<4, double>{ q3[0], q3[1], q3[2], q3[3] });
    Matrix4x4<double> m4 = Lift3To4(m3);
    Rotation<4, double> r4(m4);
    Quaternion<double> q4 = r4;
    io.outVec(Vector<4, double>{ q4[0], q4[1], q4[2], q4[3] });
}

// Axis-angle -> matrix, for N = 3 and N = 4. Calls cos and sin.
ORACLE_CASE("Rotation.axisAngleToMatrix")
{
    Vector3<double> axis = io.unit<3>();
    double angle = io.real(-7.0, 7.0);
    AxisAngle<3, double> a3(axis, angle);
    Rotation<3, double> r3(a3);
    Matrix3x3<double> m3 = r3;
    io.outMat(m3);
    Vector4<double> axis4{ axis[0], axis[1], axis[2], 0.0 };
    AxisAngle<4, double> a4(axis4, angle);
    Rotation<4, double> r4(a4);
    Matrix4x4<double> m4 = r4;
    io.outMat(m4);
}

// Matrix -> axis-angle, for N = 3 and N = 4. The axis is arithmetic only, the
// angle comes from acos. The generator rejects the inputs on which upstream's
// extraction underflows (issue #374), which the deviation case below covers.
ORACLE_CASE("Rotation.matrixToAxisAngle")
{
    int32_t mode = io.index() % 4;
    Matrix3x3<double> m3 = GivenMat(io, SoundRotMatRaw(io, mode));
    Rotation<3, double> r3(m3);
    AxisAngle<3, double> a3 = r3;
    OutAxisAngle(io, a3);
    Matrix4x4<double> m4 = Lift3To4(m3);
    Rotation<4, double> r4(m4);
    AxisAngle<4, double> a4 = r4;
    OutAxisAngle(io, a4);
}

// Axis-angle -> quaternion, for N = 3 and N = 4. Calls sin and cos.
ORACLE_CASE("Rotation.axisAngleToQuaternion")
{
    Vector3<double> axis = io.unit<3>();
    double angle = io.real(-7.0, 7.0);
    AxisAngle<3, double> a3(axis, angle);
    Rotation<3, double> r3(a3);
    Quaternion<double> q3 = r3;
    io.outVec(Vector<4, double>{ q3[0], q3[1], q3[2], q3[3] });
    Vector4<double> axis4{ axis[0], axis[1], axis[2], 0.0 };
    AxisAngle<4, double> a4(axis4, angle);
    Rotation<4, double> r4(a4);
    Quaternion<double> q4 = r4;
    io.outVec(Vector<4, double>{ q4[0], q4[1], q4[2], q4[3] });
}

// Quaternion -> axis-angle, for N = 3 and N = 4. The axis is arithmetic only
// (a division by sqrt), the angle comes from acos. Mode 0 includes the zero
// quaternion arm, where the axis is Unit(0) and the angle 0.
ORACLE_CASE("Rotation.quaternionToAxisAngle")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q = GivenQuat(io, QuatRaw(io, mode));
    Rotation<3, double> r3(q);
    AxisAngle<3, double> a3 = r3;
    OutAxisAngle(io, a3);
    Rotation<4, double> r4(q);
    AxisAngle<4, double> a4 = r4;
    OutAxisAngle(io, a4);
}

namespace
{
    // One Euler angle, unrecorded. Mode 0 is an exact multiple of pi/2 (the
    // gimbal-lock configurations), mode 2 is a multiple of pi/2 perturbed by
    // 2^-e for e in [10,40], and the rest are uniform on (-pi,pi).
    double EulerAngleRaw(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 0)
        {
            return static_cast<double>(io.rawInteger(-2, 2)) * GTE_C_HALF_PI;
        }
        if (mode == 2)
        {
            double base = static_cast<double>(io.rawInteger(-2, 2)) * GTE_C_HALF_PI;
            return base + std::ldexp(io.raw(-1.0, 1.0), -io.rawInteger(10, 40));
        }
        return io.raw(-GTE_C_PI, GTE_C_PI);
    }

    // Convert(Matrix, Quaternion) branches on r22 <= 0 and then on
    // r11 - r00 <= 0 or r11 + r00 <= 0. When the matrix is itself the output
    // of sin/cos, a 1-ulp difference between the MSVC runtime and V8 can move
    // one of those quantities across zero and select a different branch,
    // which no tolerance repairs. Composite cases that build the matrix from
    // angles therefore reject candidates whose branch quantities are within
    // 'margin' of zero; the branches are covered exactly by
    // Rotation.matrixToQuaternion, whose input matrix is recorded and
    // therefore identical on both sides.
    bool QuatBranchesClear(Matrix3x3<double> const& r, double margin)
    {
        return std::fabs(r(2, 2)) > margin
            && std::fabs(r(1, 1) - r(0, 0)) > margin
            && std::fabs(r(1, 1) + r(0, 0)) > margin;
    }

    // The same for Convert(Matrix, EulerAngles), whose gimbal-lock tests are
    // r(i2,i0) < 1 / > -1 (distinct axes) or r(i2,i2) < 1 / > -1 (repeated
    // first and last axis). Staying away from +-1 also keeps atan2 and asin
    // well conditioned, since the remaining entries of the row or column then
    // have squared sum at least 2*margin. The gimbal-lock arms themselves are
    // covered exactly by Rotation.matrixToEulerAngles.
    bool EulerBranchesClear(Matrix3x3<double> const& r, int32_t i0, int32_t i1,
        int32_t i2, double margin)
    {
        if (i1 == i0 || i1 == i2)
        {
            // The INVALID arm reads no matrix entry.
            return true;
        }
        double v = (i0 != i2 ? r(i2, i0) : r(i2, i2));
        return std::fabs(v) < 1.0 - margin;
    }
}

// Matrix -> Euler angles, for N = 3 and N = 4, for all twelve valid axis
// orders and one invalid one. The input matrix is recorded, so the
// gimbal-lock tests r(i,j) < 1 / > -1 are evaluated on identical doubles on
// both sides; mode 0 (signed permutation matrices) makes them exact ties, so
// the UNIQUE, NOT_UNIQUE_SUM and NOT_UNIQUE_DIF arms are all reached
// deterministically. Calls atan2, asin and acos.
ORACLE_CASE("Rotation.matrixToEulerAngles")
{
    int32_t mode = io.index() % 4;
    int32_t which = io.integer(0, 12);
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which, i0, i1, i2);
    Matrix3x3<double> m3 = GivenMat(io, RotMatRaw<3>(io, mode));
    Rotation<3, double> r3(m3);
    OutEulerAngles(io, r3(i0, i1, i2));
    Matrix4x4<double> m4 = Lift3To4(m3);
    Rotation<4, double> r4(m4);
    OutEulerAngles(io, r4(i0, i1, i2));
}

// Euler angles -> matrix, for N = 3 and N = 4. Three axis-angle conversions
// (cos and sin) and two matrix products; no branch depends on a libm result.
ORACLE_CASE("Rotation.eulerAnglesToMatrix")
{
    int32_t mode = io.index() % 4;
    int32_t which = io.integer(0, 12);
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which, i0, i1, i2);
    double a0 = io.given(EulerAngleRaw(io, mode));
    double a1 = io.given(EulerAngleRaw(io, mode));
    double a2 = io.given(EulerAngleRaw(io, mode));
    EulerAngles<double> e(i0, i1, i2, a0, a1, a2);
    Rotation<3, double> r3(e);
    Matrix3x3<double> m3 = r3;
    io.outMat(m3);
    Rotation<4, double> r4(e);
    Matrix4x4<double> m4 = r4;
    io.outMat(m4);
}

// Quaternion -> Euler angles, for N = 3 and N = 4. The quaternion is
// recorded, the intermediate matrix is pure arithmetic, so the gimbal-lock
// tests see identical doubles on both sides and no margin is needed.
ORACLE_CASE("Rotation.quaternionToEulerAngles")
{
    int32_t mode = io.index() % 4;
    int32_t which = io.integer(0, 12);
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which, i0, i1, i2);
    Quaternion<double> q = GivenQuat(io, UnitQuatRaw(io, mode));
    Rotation<3, double> r3(q);
    OutEulerAngles(io, r3(i0, i1, i2));
    Rotation<4, double> r4(q);
    OutEulerAngles(io, r4(i0, i1, i2));
}

// Euler angles -> quaternion, for N = 3 and N = 4. The intermediate matrix
// comes from sin and cos, so the candidate angles are rejected when the
// matrix -> quaternion branch quantities are within 1e-6 of zero.
ORACLE_CASE("Rotation.eulerAnglesToQuaternion")
{
    int32_t mode = io.index() % 4;
    int32_t which = io.integer(0, 12);
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which, i0, i1, i2);
    double a0 = 0.0, a1 = 0.0, a2 = 0.0;
    for (int32_t attempt = 0; attempt < 32; ++attempt)
    {
        a0 = EulerAngleRaw(io, mode);
        a1 = EulerAngleRaw(io, mode);
        a2 = EulerAngleRaw(io, mode);
        EulerAngles<double> candidate(i0, i1, i2, a0, a1, a2);
        Rotation<3, double> probe(candidate);
        Matrix3x3<double> pm = probe;
        if (QuatBranchesClear(pm, 1e-6))
        {
            break;
        }
    }
    io.given(a0);
    io.given(a1);
    io.given(a2);
    EulerAngles<double> e(i0, i1, i2, a0, a1, a2);
    Rotation<3, double> r3(e);
    Quaternion<double> q3 = r3;
    io.outVec(Vector<4, double>{ q3[0], q3[1], q3[2], q3[3] });
    Rotation<4, double> r4(e);
    Quaternion<double> q4 = r4;
    io.outVec(Vector<4, double>{ q4[0], q4[1], q4[2], q4[3] });
}

// Axis-angle -> Euler angles, for N = 3 and N = 4. The pair goes through a
// quaternion (sin, cos) and a matrix, so the candidate is rejected when the
// matrix -> Euler gimbal-lock tests are within 1e-6 of +-1.
ORACLE_CASE("Rotation.axisAngleToEulerAngles")
{
    int32_t which = io.integer(0, 12);
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which, i0, i1, i2);
    Vector3<double> axis{ 1.0, 0.0, 0.0 };
    double angle = 0.0;
    for (int32_t attempt = 0; attempt < 32; ++attempt)
    {
        double len = 0.0;
        do
        {
            for (int32_t i = 0; i < 3; ++i) { axis[i] = io.raw(-1.0, 1.0); }
            len = Length(axis);
        }
        while (len < 0.1 || len > 1.0);
        Normalize(axis);
        angle = io.raw(-7.0, 7.0);
        AxisAngle<3, double> candidate(axis, angle);
        Rotation<3, double> probe(candidate);
        Matrix3x3<double> pm = probe;
        if (EulerBranchesClear(pm, i0, i1, i2, 1e-6))
        {
            break;
        }
    }
    io.givenVec(axis);
    io.given(angle);
    AxisAngle<3, double> a3(axis, angle);
    Rotation<3, double> r3(a3);
    OutEulerAngles(io, r3(i0, i1, i2));
    Vector4<double> axis4{ axis[0], axis[1], axis[2], 0.0 };
    AxisAngle<4, double> a4(axis4, angle);
    Rotation<4, double> r4(a4);
    OutEulerAngles(io, r4(i0, i1, i2));
}

// Euler angles -> axis-angle, for N = 3 and N = 4. The chain is
// angles -> matrix -> quaternion -> axis-angle, so the candidate is rejected
// both when the matrix -> quaternion branch quantities are within 1e-6 of
// zero and when |q[3]| is within 1e-3 of 1. The second condition is a
// conditioning bound, not a branch: the angle is 2*acos(q[3]), and a 1-ulp
// difference in the intermediate q[3] near +-1 moves the angle by
// 1.1e-16/sqrt(2*(1-|q3|)), which is 5.5e-14 relative at 1-|q3| = 1e-3 and
// unbounded as |q3| approaches 1.
ORACLE_CASE("Rotation.eulerAnglesToAxisAngle")
{
    int32_t mode = io.index() % 4;
    int32_t which = io.integer(0, 12);
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which, i0, i1, i2);
    double a0 = 0.0, a1 = 0.0, a2 = 0.0;
    for (int32_t attempt = 0; attempt < 32; ++attempt)
    {
        a0 = EulerAngleRaw(io, mode);
        a1 = EulerAngleRaw(io, mode);
        a2 = EulerAngleRaw(io, mode);
        EulerAngles<double> candidate(i0, i1, i2, a0, a1, a2);
        Rotation<3, double> probe(candidate);
        Matrix3x3<double> pm = probe;
        Quaternion<double> pq = probe;
        if (QuatBranchesClear(pm, 1e-6) && std::fabs(pq[3]) < 1.0 - 1e-3)
        {
            break;
        }
    }
    io.given(a0);
    io.given(a1);
    io.given(a2);
    EulerAngles<double> e(i0, i1, i2, a0, a1, a2);
    Rotation<3, double> r3(e);
    AxisAngle<3, double> aa3 = r3;
    OutAxisAngle(io, aa3);
    Rotation<4, double> r4(e);
    AxisAngle<4, double> aa4 = r4;
    OutAxisAngle(io, aa4);
}

// The four identity conversions: a Rotation returns its own representation
// unchanged. Arithmetic-free, but it pins the constructors and the cached
// members.
ORACLE_CASE("Rotation.passthrough")
{
    int32_t mode = io.index() % 4;
    int32_t which = io.integer(0, 12);
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which, i0, i1, i2);
    Matrix3x3<double> m3 = GivenMat(io, RotMatRaw<3>(io, mode));
    Quaternion<double> q = GivenQuat(io, UnitQuatRaw(io, mode));
    Vector3<double> axis = io.unit<3>();
    double angle = io.real(-7.0, 7.0);
    double a0 = io.given(EulerAngleRaw(io, mode));
    double a1 = io.given(EulerAngleRaw(io, mode));
    double a2 = io.given(EulerAngleRaw(io, mode));

    Rotation<3, double> fromMatrix(m3);
    Matrix3x3<double> outMatrix = fromMatrix;
    io.outMat(outMatrix);

    Rotation<3, double> fromQuaternion(q);
    Quaternion<double> outQuaternion = fromQuaternion;
    io.outVec(Vector<4, double>{ outQuaternion[0], outQuaternion[1],
        outQuaternion[2], outQuaternion[3] });

    AxisAngle<3, double> a3(axis, angle);
    Rotation<3, double> fromAxisAngle(a3);
    AxisAngle<3, double> outAxisAngle = fromAxisAngle;
    OutAxisAngle(io, outAxisAngle);

    EulerAngles<double> e(i0, i1, i2, a0, a1, a2);
    Rotation<3, double> fromEuler(e);
    OutEulerAngles(io, fromEuler(i0, i1, i2));
}

// Deliberate port deviation, issue #374: upstream's 'angle in (0,pi)' arm
// extracts the axis from R - Transpose(R), whose entries are 2*sin(angle)*a_i.
// The inputs here are a coordinate-axis rotation by pi (an exact matrix with
// entries 0 and +-1) whose trace is lifted by 2^-50, so that acos reports an
// angle just below pi and the antisymmetric arm runs, plus an antisymmetric
// perturbation of about 2^-k with k in [513,600]. Upstream then squares axis
// components of about 2^(1-k): for k in [513,538] those squares are subnormal
// and lose most of their significand, so Normalize returns an axis that is
// not unit length (the 0.707 case of the finding), and for k >= 539 they
// underflow to zero, so Normalize returns the zero vector. Either way the
// reported axis is not a rotation axis, and converting it back does not
// reproduce the input rotation. The port validates the normalized axis and
// falls back to the symmetric formula, which is exact here. Random rotation
// matrices never reach this regime, which is why the main case's probe -- the
// same predicate -- accepts essentially every candidate.
ORACLE_CASE("Rotation.matrixToAxisAngle.nearPi")
{
    int32_t d = io.rawInteger(0, 2);
    int32_t k = io.rawInteger(513, 600);
    double e0 = std::ldexp(io.raw(1.0, 2.0), -k);
    double e1 = std::ldexp(io.raw(1.0, 2.0), -k);
    double e2 = std::ldexp(io.raw(1.0, 2.0), -k);
    double delta = std::ldexp(1.0, -50);
    Matrix3x3<double> m3{};
    m3.MakeZero();
    for (int32_t i = 0; i < 3; ++i)
    {
        m3(i, i) = (i == d ? 1.0 : -1.0);
    }
    m3((d + 1) % 3, (d + 1) % 3) += delta;
    m3(2, 1) = e0;
    m3(1, 2) = -e0;
    m3(0, 2) = e1;
    m3(2, 0) = -e1;
    m3(1, 0) = e2;
    m3(0, 1) = -e2;
    GivenMat(io, m3);
    Rotation<3, double> r3(m3);
    AxisAngle<3, double> a3 = r3;
    OutAxisAngle(io, a3);
    Matrix4x4<double> m4 = Lift3To4(m3);
    Rotation<4, double> r4(m4);
    AxisAngle<4, double> a4 = r4;
    OutAxisAngle(io, a4);
}

// Deliberate port deviation, issue #225: upstream's operator()(i0,i1,i2)
// writes the requested axis indices into the cached EulerAngles and then does
// nothing at all when the Rotation was built from Euler angles, so asking an
// Euler-sourced Rotation for a *different* factorization returns the original
// angles relabelled with the new axes, which is a different rotation. The
// port recomputes through the rotation matrix, as every other source type
// does. The two axis orders here are always different.
ORACLE_CASE("Rotation.eulerAngles.refactorization")
{
    int32_t mode = io.index() % 4;
    int32_t which0 = io.integer(0, 11);
    int32_t step = io.integer(1, 11);
    int32_t which1 = (which0 + step) % 12;
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which0, i0, i1, i2);
    int32_t j0 = 0, j1 = 0, j2 = 0;
    EulerAxes(which1, j0, j1, j2);
    double a0 = io.given(EulerAngleRaw(io, mode));
    double a1 = io.given(EulerAngleRaw(io, mode));
    double a2 = io.given(EulerAngleRaw(io, mode));
    EulerAngles<double> e(i0, i1, i2, a0, a1, a2);
    Rotation<3, double> r3(e);
    OutEulerAngles(io, r3(j0, j1, j2));
}

// ===================================================================
// RotationEstimate.h
// ===================================================================

namespace
{
#define ORACLE_ROTC_DISPATCH(NAME)                                  \
    double NAME##Dispatch(double t, int32_t degree)                 \
    {                                                               \
        switch (degree)                                             \
        {                                                           \
        case 4:  return NAME##Estimate<double, 4>(t);               \
        case 6:  return NAME##Estimate<double, 6>(t);               \
        case 8:  return NAME##Estimate<double, 8>(t);               \
        case 10: return NAME##Estimate<double, 10>(t);              \
        case 12: return NAME##Estimate<double, 12>(t);              \
        case 14: return NAME##Estimate<double, 14>(t);              \
        default: return NAME##Estimate<double, 16>(t);              \
        }                                                           \
    }                                                               \
    double NAME##MaxErrorDispatch(int32_t degree)                   \
    {                                                               \
        switch (degree)                                             \
        {                                                           \
        case 4:  return Get##NAME##EstimateMaxError<double, 4>();   \
        case 6:  return Get##NAME##EstimateMaxError<double, 6>();   \
        case 8:  return Get##NAME##EstimateMaxError<double, 8>();   \
        case 10: return Get##NAME##EstimateMaxError<double, 10>();  \
        case 12: return Get##NAME##EstimateMaxError<double, 12>();  \
        case 14: return Get##NAME##EstimateMaxError<double, 14>();  \
        default: return Get##NAME##EstimateMaxError<double, 16>();  \
        }                                                           \
    }

    ORACLE_ROTC_DISPATCH(RotC0)
    ORACLE_ROTC_DISPATCH(RotC1)
    ORACLE_ROTC_DISPATCH(RotC2)
    ORACLE_ROTC_DISPATCH(RotC3)
    ORACLE_ROTC_DISPATCH(RotC4)
#undef ORACLE_ROTC_DISPATCH

    void RotationEstimateDispatch(Vector3<double> const& p, int32_t degree,
        Matrix3x3<double>& R)
    {
        switch (degree)
        {
        case 4:  RotationEstimate<double, 4>(p, R); break;
        case 6:  RotationEstimate<double, 6>(p, R); break;
        case 8:  RotationEstimate<double, 8>(p, R); break;
        case 10: RotationEstimate<double, 10>(p, R); break;
        case 12: RotationEstimate<double, 12>(p, R); break;
        case 14: RotationEstimate<double, 14>(p, R); break;
        default: RotationEstimate<double, 16>(p, R); break;
        }
    }

    void RotationDerivativeEstimateDispatch(Vector3<double> const& p,
        int32_t degree, std::array<Matrix3x3<double>, 3>& Rder)
    {
        switch (degree)
        {
        case 4:  RotationDerivativeEstimate<double, 4>(p, Rder); break;
        case 6:  RotationDerivativeEstimate<double, 6>(p, Rder); break;
        case 8:  RotationDerivativeEstimate<double, 8>(p, Rder); break;
        case 10: RotationDerivativeEstimate<double, 10>(p, Rder); break;
        case 12: RotationDerivativeEstimate<double, 12>(p, Rder); break;
        case 14: RotationDerivativeEstimate<double, 14>(p, Rder); break;
        default: RotationDerivativeEstimate<double, 16>(p, Rder); break;
        }
    }

    void RotationAndDerivativeEstimateDispatch(Vector3<double> const& p,
        int32_t degree, Matrix3x3<double>& R,
        std::array<Matrix3x3<double>, 3>& Rder)
    {
        switch (degree)
        {
        case 4:  RotationAndDerivativeEstimate<double, 4>(p, R, Rder); break;
        case 6:  RotationAndDerivativeEstimate<double, 6>(p, R, Rder); break;
        case 8:  RotationAndDerivativeEstimate<double, 8>(p, R, Rder); break;
        case 10: RotationAndDerivativeEstimate<double, 10>(p, R, Rder); break;
        case 12: RotationAndDerivativeEstimate<double, 12>(p, R, Rder); break;
        case 14: RotationAndDerivativeEstimate<double, 14>(p, R, Rder); break;
        default: RotationAndDerivativeEstimate<double, 16>(p, R, Rder); break;
        }
    }

    // The argument of the rotc_k estimates, documented as t in [0,pi]. Mode 0
    // pins the ends (0, pi and the small lattice), mode 3 is wild, and the
    // rest are uniform on [0,pi] or, for a quarter of the records, outside
    // the documented range, where the polynomials are still defined.
    double RotCArgument(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 0)
        {
            int32_t k = io.rawInteger(0, 4);
            double table[5] = { 0.0, 1.0, 2.0, GTE_C_HALF_PI, GTE_C_PI };
            return io.given(table[k]);
        }
        if (mode == 3)
        {
            return io.given(WildRaw(io));
        }
        if (io.index() % 8 == 5)
        {
            return io.real(-8.0, 8.0);
        }
        return io.real(0.0, GTE_C_PI);
    }
}

// The five polynomial estimates and the five max-error tables, for all seven
// degrees. Polynomial evaluation only; the tables are constants.
ORACLE_CASE("RotationEstimate.rotCEstimate")
{
    int32_t mode = io.index() % 4;
    int32_t degree = 4 + 2 * io.integer(0, 6);
    double t = RotCArgument(io, mode);
    io.outReal(RotC0Dispatch(t, degree));
    io.outReal(RotC1Dispatch(t, degree));
    io.outReal(RotC2Dispatch(t, degree));
    io.outReal(RotC3Dispatch(t, degree));
    io.outReal(RotC4Dispatch(t, degree));
    io.outReal(RotC0MaxErrorDispatch(degree));
    io.outReal(RotC1MaxErrorDispatch(degree));
    io.outReal(RotC2MaxErrorDispatch(degree));
    io.outReal(RotC3MaxErrorDispatch(degree));
    io.outReal(RotC4MaxErrorDispatch(degree));
}

// RotationEstimate, RotationDerivativeEstimate and the combined function.
// Arithmetic only: Length uses sqrt and the estimates are polynomials.
ORACLE_CASE("RotationEstimate.rotationEstimate")
{
    int32_t mode = io.index() % 4;
    int32_t degree = 4 + 2 * io.integer(0, 6);
    Vector3<double> p = Vec<3>(io, mode);

    Matrix3x3<double> R{};
    RotationEstimateDispatch(p, degree, R);
    io.outMat(R);

    std::array<Matrix3x3<double>, 3> Rder{};
    RotationDerivativeEstimateDispatch(p, degree, Rder);
    for (int32_t i = 0; i < 3; ++i) { io.outMat(Rder[i]); }

    Matrix3x3<double> Rboth{};
    std::array<Matrix3x3<double>, 3> RderBoth{};
    RotationAndDerivativeEstimateDispatch(p, degree, Rboth, RderBoth);
    io.outMat(Rboth);
    for (int32_t i = 0; i < 3; ++i) { io.outMat(RderBoth[i]); }
}

// ===================================================================
// Transform.h
// ===================================================================

namespace
{
    // A nonzero scalar; SetScale and SetUniformScale assert on a zero scale.
    // Records one double.
    double NonzeroScalar(oracle::Ctx& io, int32_t mode)
    {
        double s = 0.0;
        do
        {
            if (mode == 0) { s = static_cast<double>(io.rawInteger(-3, 3)); }
            else if (mode == 1) { s = io.raw(-10.0, 10.0); }
            else if (mode == 2) { s = static_cast<double>(io.rawInteger(-4, 4)); }
            else { s = WildRaw(io); }
        }
        while (s == 0.0);
        return io.given(s);
    }

    // A 4x4 matrix for SetMatrix. Mode 2 makes the upper-left 3x3 block
    // exactly singular (its third row an integer combination of the first
    // two), which is what decides whether the homogeneous matrix is
    // invertible: UpdateHMatrix copies only that block and the translation,
    // so the last row of H is always (0,0,0,1). Sixteen recorded doubles,
    // row-major, in either mode.
    Matrix4x4<double> GeneralMat4(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 2)
        {
            Matrix4x4<double> m{};
            for (int32_t r = 0; r < 2; ++r)
            {
                for (int32_t c = 0; c < 3; ++c)
                {
                    m(r, c) = static_cast<double>(io.rawInteger(-3, 3));
                }
            }
            double w0 = static_cast<double>(io.rawInteger(-2, 2));
            double w1 = static_cast<double>(io.rawInteger(-2, 2));
            for (int32_t c = 0; c < 3; ++c)
            {
                m(2, c) = w0 * m(0, c) + w1 * m(1, c);
            }
            for (int32_t r = 0; r < 4; ++r)
            {
                m(r, 3) = static_cast<double>(io.rawInteger(-3, 3));
            }
            for (int32_t c = 0; c < 3; ++c)
            {
                m(3, c) = static_cast<double>(io.rawInteger(-3, 3));
            }
            return GivenMat(io, m);
        }
        return SquareMat<4>(io, mode);
    }

    // The channels that every Transform defines.
    void OutTransform(oracle::Ctx& io, Transform<double> const& t)
    {
        io.outBool(t.IsIdentity());
        io.outBool(t.IsRSMatrix());
        io.outBool(t.IsUniformScale());
        io.outMat(t.GetHMatrix());
        io.outMat(t.GetHInverse());
        io.outReal(t.GetNorm());
    }

    // The full accessor set. The rotation-only accessors are emitted only
    // when the is-rsmatrix hint holds, since otherwise they assert; the
    // emitted count therefore follows from the flags emitted first. The
    // matrix -> quaternion conversion is arithmetic only, so this dump stays
    // bit-exact as long as the transform's matrix channel is a recorded
    // input rather than a libm result.
    void OutTransformFull(oracle::Ctx& io, Transform<double> const& t)
    {
        OutTransform(io, t);
        io.outMat(t.GetMatrix());
        io.outVec(t.GetTranslation());
        io.outVec(t.GetTranslationW0());
        io.outVec(t.GetTranslationW1());
        if (t.IsRSMatrix())
        {
            io.outMat(t.GetRotation());
            Matrix3x3<double> rotate3{};
            t.GetRotation(rotate3);
            io.outMat(rotate3);
            io.outVec(t.GetScale());
            io.outVec(t.GetScaleW1());
            Quaternion<double> q{};
            t.GetRotation(q);
            io.outVec(Vector<4, double>{ q[0], q[1], q[2], q[3] });
        }
        if (t.IsUniformScale())
        {
            io.outReal(t.GetUniformScale());
        }
    }

    // Apply the translation channel through one of the three SetTranslation
    // overloads, selected by the record index.
    void ApplyTranslation(oracle::Ctx& io, Transform<double>& t,
        Vector3<double> const& tr)
    {
        int32_t which = io.index() % 3;
        if (which == 0)
        {
            t.SetTranslation(tr[0], tr[1], tr[2]);
        }
        else if (which == 1)
        {
            t.SetTranslation(tr);
        }
        else
        {
            Vector4<double> tr4{ tr[0], tr[1], tr[2], 1.0 };
            t.SetTranslation(tr4);
        }
    }

    // Apply the scale channel, selected by the record index:
    //   0  leave the unit scale alone
    //   1  SetScale(s0,s1,s2), the non-uniform path
    //   2  SetUniformScale(us)
    //   3  SetScale(Vector3) followed by MakeUnitScale
    void ApplyScale(oracle::Ctx& io, Transform<double>& t, double s0, double s1,
        double s2, double us)
    {
        int32_t which = io.index() % 4;
        if (which == 1)
        {
            t.SetScale(s0, s1, s2);
        }
        else if (which == 2)
        {
            t.SetUniformScale(us);
        }
        else if (which == 3)
        {
            Vector3<double> s{ s0, s1, s2 };
            t.SetScale(s);
            t.MakeUnitScale();
        }
    }
}

// The rotation channel set from a recorded rotation matrix (the 4x4 overload
// on even records, the 3x3 overload on odd ones), the translation channel
// through each of its three overloads and the scale channel through each of
// its forms. Arithmetic only.
ORACLE_CASE("Transform.rotationMatrix")
{
    int32_t mode = io.index() % 4;
    Matrix3x3<double> rot3 = GivenMat(io, RotMatRaw<3>(io, mode));
    Vector3<double> tr = Vec<3>(io, mode);
    double s0 = NonzeroScalar(io, mode);
    double s1 = NonzeroScalar(io, mode);
    double s2 = NonzeroScalar(io, mode);
    double us = NonzeroScalar(io, mode);

    Transform<double> t{};
    OutTransform(io, t);
    if (io.index() % 2 == 0)
    {
        Matrix4x4<double> rot4 = Lift3To4(rot3);
        t.SetRotation(rot4);
    }
    else
    {
        t.SetRotation(rot3);
    }
    ApplyTranslation(io, t, tr);
    ApplyScale(io, t, s0, s1, s2, us);
    OutTransformFull(io, t);
    t.MakeIdentity();
    OutTransform(io, t);
}

// SetRotation(Quaternion) and GetRotation(Quaternion). Both conversions are
// arithmetic only.
ORACLE_CASE("Transform.rotationQuaternion")
{
    int32_t mode = io.index() % 4;
    Quaternion<double> q = GivenQuat(io, UnitQuatRaw(io, mode));
    Vector3<double> tr = Vec<3>(io, mode);
    double us = NonzeroScalar(io, mode);

    Transform<double> t{};
    t.SetRotation(q);
    ApplyTranslation(io, t, tr);
    t.SetUniformScale(us);
    OutTransformFull(io, t);
}

// The general-matrix path: SetMatrix clears the is-rsmatrix and
// is-uniform-scale hints, GetNorm switches to the max-row-sum norm and
// GetHInverse takes the Matrix4x4 closed-form Inverse of the whole
// homogeneous matrix. Mode 2 supplies exactly singular matrices, for which
// that Inverse is the zero matrix. Arithmetic only.
ORACLE_CASE("Transform.setMatrix")
{
    int32_t mode = io.index() % 4;
    Matrix4x4<double> m = GeneralMat4(io, mode);
    Vector3<double> tr = Vec<3>(io, mode);

    Transform<double> t{};
    t.SetMatrix(m);
    ApplyTranslation(io, t, tr);
    OutTransformFull(io, t);
}

// SetRotation(AxisAngle<3>) and SetRotation(AxisAngle<4>), with both
// GetRotation(AxisAngle<N>) forms. Calls sin, cos and acos. The angle is kept
// in [0.3, 2.8] in absolute value so that the matrix -> axis-angle
// conversion stays away from the 'angle == 0' and 'angle == pi' arms, whose
// selection would otherwise depend on a 1-ulp difference in cos.
ORACLE_CASE("Transform.rotationAxisAngle")
{
    Vector3<double> axis = io.unit<3>();
    double magnitude = io.real(0.3, 2.8);
    double sign = io.given(io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    double angle = magnitude * sign;
    Vector3<double> tr = Vec<3>(io, 1);

    Transform<double> t{};
    if (io.index() % 2 == 0)
    {
        AxisAngle<3, double> a3(axis, angle);
        t.SetRotation(a3);
    }
    else
    {
        Vector4<double> axis4{ axis[0], axis[1], axis[2], 0.0 };
        AxisAngle<4, double> a4(axis4, angle);
        t.SetRotation(a4);
    }
    ApplyTranslation(io, t, tr);
    OutTransform(io, t);

    AxisAngle<3, double> out3{};
    t.GetRotation(out3);
    OutAxisAngle(io, out3);
    AxisAngle<4, double> out4{};
    t.GetRotation(out4);
    OutAxisAngle(io, out4);
}

// SetRotation(EulerAngles) and GetRotation(EulerAngles). Calls sin, cos,
// atan2, asin and acos. The angles are rejected when the gimbal-lock tests of
// the requested factorization would be within 1e-6 of +-1 on the matrix that
// sin and cos produce.
ORACLE_CASE("Transform.rotationEulerAngles")
{
    int32_t mode = io.index() % 4;
    int32_t which = io.integer(0, 11);
    int32_t i0 = 0, i1 = 0, i2 = 0;
    EulerAxes(which, i0, i1, i2);
    double a0 = 0.0, a1 = 0.0, a2 = 0.0;
    for (int32_t attempt = 0; attempt < 32; ++attempt)
    {
        a0 = EulerAngleRaw(io, mode);
        a1 = EulerAngleRaw(io, mode);
        a2 = EulerAngleRaw(io, mode);
        EulerAngles<double> candidate(i0, i1, i2, a0, a1, a2);
        Rotation<3, double> probe(candidate);
        Matrix3x3<double> pm = probe;
        if (EulerBranchesClear(pm, i0, i1, i2, 1e-6))
        {
            break;
        }
    }
    io.given(a0);
    io.given(a1);
    io.given(a2);
    Vector3<double> tr = Vec<3>(io, 1);

    EulerAngles<double> e(i0, i1, i2, a0, a1, a2);
    Transform<double> t{};
    t.SetRotation(e);
    ApplyTranslation(io, t, tr);
    OutTransform(io, t);

    EulerAngles<double> out{};
    out.axis = { i0, i1, i2 };
    t.GetRotation(out);
    OutEulerAngles(io, out);
}

namespace
{
    // Build a transform from recorded channels. 'kind' selects the structure:
    //   0  identity
    //   1  rotation + translation, unit scale (RS, uniform)
    //   2  rotation + translation + uniform scale (RS, uniform)
    //   3  rotation + translation + non-uniform scale (RS, not uniform)
    //   4  general matrix + translation (not RS)
    // The caller draws all the channels first, so every kind records the same
    // number of doubles.
    Transform<double> MakeTransform(int32_t kind, Matrix3x3<double> const& rot3,
        Matrix4x4<double> const& general, Vector3<double> const& tr,
        double s0, double s1, double s2, double us)
    {
        Transform<double> t{};
        if (kind == 0)
        {
            return t;
        }
        if (kind == 4)
        {
            t.SetMatrix(general);
        }
        else
        {
            t.SetRotation(rot3);
        }
        t.SetTranslation(tr[0], tr[1], tr[2]);
        if (kind == 2)
        {
            t.SetUniformScale(us);
        }
        else if (kind == 3)
        {
            t.SetScale(s0, s1, s2);
        }
        return t;
    }
}

// Transform::Inverse. The uniform-scale arm reassigns the channels (rotation,
// uniform scale, translation); every other state falls back to the
// Matrix4x4 closed-form Inverse of the homogeneous matrix. The M channel of
// the general arm is a deliberate port deviation and is compared by the
// separate case below, so it is not emitted here. Arithmetic only.
ORACLE_CASE("Transform.inverse")
{
    int32_t mode = io.index() % 4;
    int32_t kind = io.integer(0, 4);
    Matrix3x3<double> rot3 = GivenMat(io, RotMatRaw<3>(io, mode));
    Matrix4x4<double> general = GeneralMat4(io, mode);
    Vector3<double> tr = Vec<3>(io, mode);
    double s0 = NonzeroScalar(io, mode);
    double s1 = NonzeroScalar(io, mode);
    double s2 = NonzeroScalar(io, mode);
    double us = NonzeroScalar(io, mode);

    Transform<double> t = MakeTransform(kind, rot3, general, tr, s0, s1, s2, us);
    io.outMat(t.GetHInverse());
    Transform<double> inv = t.Inverse();
    OutTransform(io, inv);
    io.outVec(inv.GetTranslation());
    if (inv.IsRSMatrix())
    {
        io.outMat(inv.GetRotation());
        io.outVec(inv.GetScale());
    }
    if (inv.IsUniformScale())
    {
        io.outReal(inv.GetUniformScale());
    }
}

// Deliberate port deviation, issue #265: the general arm of Transform::Inverse
// passes the full affine 4x4 to SetMatrix, so the inverse's M channel carries
// the translation in its last column and a last row that is whatever the
// closed-form Inverse produced, violating GetMatrix's documented
// {{M,0},{0,1}} structure. The port zeroes the last column and row before
// storing. The generator uses the general path with a nonzero translation, so
// the M channel differs on essentially every record.
ORACLE_CASE("Transform.inverse.matrixChannel")
{
    int32_t mode = io.index() % 2;
    Matrix4x4<double> general = GeneralMat4(io, mode);
    Vector3<double> tr = Vec<3>(io, mode);
    Transform<double> t{};
    t.SetMatrix(general);
    t.SetTranslation(tr[0], tr[1], tr[2]);
    Transform<double> inv = t.Inverse();
    io.outMat(inv.GetMatrix());
}

// Deliberate port deviation, issue #265: the RS arms of GetHInverse write only
// the upper-left 3x3 block and the last column of mInvHMatrix, assuming the
// last row still holds (0,0,0,1) from the constructor. After the general arm
// has assigned Inverse(mHMatrix) for a singular mHMatrix -- the zero matrix --
// that assumption is false, so a later RS state reports an inverse whose last
// row is (0,0,0,0). The port writes the last row explicitly. The generator
// forces exactly that sequence: a singular general matrix, a GetHInverse that
// caches the zero matrix, then a rotation and a second GetHInverse.
ORACLE_CASE("Transform.getHInverse.lastRow")
{
    Matrix4x4<double> singular = GeneralMat4(io, 2);
    Matrix3x3<double> rot3 = GivenMat(io, RotMatRaw<3>(io, io.index() % 4));
    Vector3<double> tr = Vec<3>(io, 0);
    Transform<double> t{};
    t.SetMatrix(singular);
    io.outMat(t.GetHInverse());
    t.SetRotation(rot3);
    t.SetTranslation(tr[0], tr[1], tr[2]);
    io.outMat(t.GetHInverse());
}

// The operator* overloads: M*V, V^T*M, A*B for transforms (all five
// structures against all five, chosen per record), Matrix4x4*Transform and
// Transform*Matrix4x4. Arithmetic only.
ORACLE_CASE("Transform.multiply")
{
    int32_t mode = io.index() % 4;
    int32_t kindA = io.integer(0, 4);
    int32_t kindB = io.integer(0, 4);
    Matrix3x3<double> rotA = GivenMat(io, RotMatRaw<3>(io, mode));
    Matrix3x3<double> rotB = GivenMat(io, RotMatRaw<3>(io, mode));
    Matrix4x4<double> generalA = GeneralMat4(io, mode);
    Matrix4x4<double> generalB = GeneralMat4(io, mode);
    Vector3<double> trA = Vec<3>(io, mode);
    Vector3<double> trB = Vec<3>(io, mode);
    double s0 = NonzeroScalar(io, mode);
    double s1 = NonzeroScalar(io, mode);
    double s2 = NonzeroScalar(io, mode);
    double us = NonzeroScalar(io, mode);
    Vector4<double> V = Vec<4>(io, mode);
    Matrix4x4<double> K = SquareMat<4>(io, mode);

    Transform<double> A = MakeTransform(kindA, rotA, generalA, trA, s0, s1, s2, us);
    Transform<double> B = MakeTransform(kindB, rotB, generalB, trB, s2, s0, s1, us);
    // The transform product is formed before anything is emitted: when both
    // factors are rotation-scale and A's scale is uniform, it calls SetScale
    // or SetUniformScale on the product of the two scales, which asserts if
    // that product underflows to zero -- 21 of the 2000 deep-run records, all
    // from the wild scale population. A record that throws records no
    // outputs, so no output may precede the throwing call.
    Transform<double> product = A * B;
    io.outVec(A * V);
    io.outVec(V * A);
    OutTransform(io, product);
    io.outVec(product.GetTranslation());
    io.outMat(K * A);
    io.outMat(A * K);
}

// Throw parity for the LogAssert preconditions. The scenario is selected by
// the record index: a zero component in SetScale, SetScale on a transform
// that is not rotation-scale, GetUniformScale on a non-uniform transform and
// a valid sequence that returns normally.
ORACLE_CASE("Transform.assertions")
{
    int32_t which = io.index() % 4;
    Matrix3x3<double> rot3 = GivenMat(io, RotMatRaw<3>(io, io.index() % 4));
    Matrix4x4<double> general = SquareMat<4>(io, 1);
    double s = io.real(1.0, 3.0);

    Transform<double> t{};
    t.SetRotation(rot3);
    if (which == 0)
    {
        t.SetScale(0.0, s, s);
    }
    else if (which == 1)
    {
        t.SetMatrix(general);
        t.SetScale(s, s, s);
    }
    else if (which == 2)
    {
        t.SetScale(s, s + 1.0, s + 2.0);
        io.outReal(t.GetUniformScale());
    }
    else
    {
        t.SetUniformScale(s);
    }
    OutTransform(io, t);
}

// ===================================================================
// ConvertCoordinates.h
// ===================================================================

namespace
{
    template <int32_t N>
    void ConvertCoordinatesBody(oracle::Ctx& io, int32_t mode, bool vU, bool vV)
    {
        Matrix<N, N, double> U = SquareMat<N>(io, mode);
        Matrix<N, N, double> V = SquareMat<N>(io, mode);
        Vector<N, double> X = Vec<N>(io, mode);
        Vector<N, double> Y = Vec<N>(io, mode);
        Matrix<N, N, double> A = SquareMat<N>(io, mode);
        Matrix<N, N, double> B = SquareMat<N>(io, mode);

        ConvertCoordinates<N, double> convert;
        io.outBool(convert(U, vU, V, vV));
        io.outMat(convert.GetC());
        io.outMat(convert.GetInverseC());
        io.outBool(convert.IsVectorOnRightU());
        io.outBool(convert.IsVectorOnRightV());
        io.outBool(convert.IsRightHandedU());
        io.outBool(convert.IsRightHandedV());
        io.outVec(convert.UToV(X));
        io.outVec(convert.VToU(Y));
        io.outMat(convert.UToV(A));
        io.outMat(convert.VToU(B));
    }
}

// The change-of-basis operator(), the six accessors, the vector conversions
// and the matrix conversions for all four combinations of the
// vector-on-the-right conventions, in dimensions 2, 3 and 4. Mode 2 supplies
// exactly singular U and V, for which operator() returns false and leaves the
// identity in place. GaussianElimination is arithmetic only.
ORACLE_CASE("ConvertCoordinates.convert")
{
    int32_t mode = io.index() % 4;
    int32_t n = io.integer(2, 4);
    bool vU = io.boolean();
    bool vV = io.boolean();
    if (n == 2)
    {
        ConvertCoordinatesBody<2>(io, mode, vU, vV);
    }
    else if (n == 3)
    {
        ConvertCoordinatesBody<3>(io, mode, vU, vV);
    }
    else
    {
        ConvertCoordinatesBody<4>(io, mode, vU, vV);
    }
}

// ===================================================================
// Projection.h
// ===================================================================

namespace
{
    // An orthonormal frame, unrecorded: the columns of a rotation matrix.
    void Frame(oracle::Ctx& io, int32_t mode, std::array<Vector3<double>, 3>& f)
    {
        Matrix3x3<double> r = RotMatRaw<3>(io, mode);
        for (int32_t c = 0; c < 3; ++c)
        {
            f[c] = r.GetCol(c);
        }
    }
}

// Orthogonal projection of an ellipse onto a line. Arithmetic only (sqrt).
ORACLE_CASE("Projection.projectEllipse2")
{
    Ellipse2<double> ellipse{};
    ellipse.center = io.vec<2>(-10.0, 10.0);
    Vector2<double> axis0 = io.unit<2>();
    // GTE's Perp(x,y) is (y,-x), so {axis0, Perp(axis0)} is left-handed;
    // negate it to make the ellipse frame right-handed.
    Vector2<double> axis1 = io.givenVec(-Perp(axis0));
    ellipse.axis[0] = axis0;
    ellipse.axis[1] = axis1;
    ellipse.extent = io.vec<2>(0.25, 4.0);
    Vector2<double> origin = io.vec<2>(-10.0, 10.0);
    Vector2<double> direction = io.unit<2>();
    Line2<double> line(origin, direction);
    double smin = 0.0, smax = 0.0;
    Project(ellipse, line, smin, smax);
    io.outReal(smin);
    io.outReal(smax);
}

// Orthogonal projection of an ellipsoid onto a line. Arithmetic only (sqrt).
ORACLE_CASE("Projection.projectEllipsoid3")
{
    int32_t mode = io.index() % 4;
    std::array<Vector3<double>, 3> f{};
    Frame(io, mode, f);
    Ellipsoid3<double> ellipsoid{};
    ellipsoid.center = io.vec<3>(-10.0, 10.0);
    ellipsoid.axis[0] = io.givenVec(f[0]);
    ellipsoid.axis[1] = io.givenVec(f[1]);
    ellipsoid.axis[2] = io.givenVec(f[2]);
    ellipsoid.extent = io.vec<3>(0.25, 4.0);
    Vector3<double> origin = io.vec<3>(-10.0, 10.0);
    Vector3<double> direction = io.unit<3>();
    Line3<double> line(origin, direction);
    double smin = 0.0, smax = 0.0;
    Project(ellipsoid, line, smin, smax);
    io.outReal(smin);
    io.outReal(smax);
}

namespace
{
    void OutEllipse2(oracle::Ctx& io, Ellipse2<double> const& e)
    {
        io.outVec(e.center);
        io.outVec(e.axis[0]);
        io.outVec(e.axis[1]);
        io.outVec(e.extent);
    }

    // The configuration required by PerspectiveProject: an ellipsoid strictly
    // between the eyepoint E and the view plane Dot(N,X) = constant. The
    // ellipsoid centre is placed at E + (maxExtent + d)*N with d in [3,8] and
    // the plane at near distance maxExtent + d + gap with gap in [1,6], so the
    // documented precondition (the projection of the ellipsoid onto the line
    // E + s*N has smin > 0, and the whole ellipsoid is in front of the plane)
    // holds by construction. The centre, the near distance and the plane
    // constant are recorded, so the replay never recomputes them.
    struct ProjectionSetup
    {
        Ellipsoid3<double> ellipsoid;
        Vector3<double> E, N, U, V;
        double constant, n;
    };

    ProjectionSetup MakeProjectionSetup(oracle::Ctx& io, int32_t mode)
    {
        ProjectionSetup s{};
        std::array<Vector3<double>, 3> frame{};
        Frame(io, mode, frame);
        s.N = io.givenVec(frame[0]);
        s.U = io.givenVec(frame[1]);
        s.V = io.givenVec(frame[2]);
        std::array<Vector3<double>, 3> eframe{};
        Frame(io, (mode + 1) % 4, eframe);
        s.ellipsoid.axis[0] = io.givenVec(eframe[0]);
        s.ellipsoid.axis[1] = io.givenVec(eframe[1]);
        s.ellipsoid.axis[2] = io.givenVec(eframe[2]);
        s.ellipsoid.extent = io.vec<3>(0.25, 2.0);
        s.E = io.vec<3>(-5.0, 5.0);
        double maxExtent = std::max(std::max(s.ellipsoid.extent[0],
            s.ellipsoid.extent[1]), s.ellipsoid.extent[2]);
        double d = io.raw(3.0, 8.0);
        double gap = io.raw(1.0, 6.0);
        s.ellipsoid.center = io.givenVec(s.E + (maxExtent + d) * s.N);
        s.n = io.given(maxExtent + d + gap);
        s.constant = io.given(Dot(s.N, s.E) + s.n);
        return s;
    }
}

// Perspective projection of an ellipsoid onto a plane, the precomputed-basis
// overload. Arithmetic only: the conic coefficients are + - * /, and
// Ellipse2's FromCoefficients runs a SymmetricEigensolver whose Householder
// and QL steps use nothing but sqrt and fabs. Inside FromCoefficients,
// Inverse(A) for the 2x2 A resolves to the closed form of Matrix2x2.h, which
// this translation unit makes visible -- not to the Gaussian-elimination
// template of Matrix.h. That is measurable here: before the port's matching
// fix landed (PR #508) the ellipse centre, which is the only output that uses
// the inverse, differed by 1 to 2 ulps on 12 of the 20 committed records
// while every axis and extent already agreed bit for bit. Upstream discards
// the bool that FromCoefficients returns, so it is not an output here; the
// TypeScript replay throws when the port's flag is false, which would show up
// as a disagreement.
ORACLE_CASE("Projection.perspectiveProject.basis")
{
    int32_t mode = io.index() % 4;
    ProjectionSetup s = MakeProjectionSetup(io, mode);
    Ellipse2<double> ellipse{};
    PerspectiveProject(s.ellipsoid, s.E, s.N, s.U, s.V, s.n, ellipse);
    OutEllipse2(io, ellipse);
}

// The single-plane overload, which derives U, V and the near distance from
// the plane with ComputeOrthogonalComplement.
ORACLE_CASE("Projection.perspectiveProject.plane")
{
    int32_t mode = io.index() % 4;
    ProjectionSetup s = MakeProjectionSetup(io, mode);
    Plane3<double> plane(s.N, s.constant);
    Ellipse2<double> ellipse{};
    PerspectiveProject(s.ellipsoid, s.E, plane, ellipse);
    OutEllipse2(io, ellipse);
}

// ===================================================================
// Matrix4x4.h, the special matrices
// ===================================================================

namespace
{
    // A 4-tuple for the special-matrix builders; four recorded doubles in
    // every mode. Mode 1 is the documented shape (a point with w = 1, or a
    // unit-length 3-direction lifted with w = 0); the other modes are the
    // general populations, which also produce the degenerate Dot(N,D) == 0
    // and signed zeros.
    Vector4<double> Special4(oracle::Ctx& io, int32_t mode, bool isPoint)
    {
        if (mode != 1)
        {
            return Vec<4>(io, mode);
        }
        Vector4<double> v{};
        if (isPoint)
        {
            Vector3<double> p = io.vec<3>(-10.0, 10.0);
            v = Vector4<double>{ p[0], p[1], p[2], 0.0 };
        }
        else
        {
            Vector3<double> u = io.unit<3>();
            v = Vector4<double>{ u[0], u[1], u[2], 0.0 };
        }
        v[3] = io.given(isPoint ? 1.0 : 0.0);
        return v;
    }
}

// The oblique projection onto a plane. Arithmetic only.
ORACLE_CASE("Matrix4x4.makeObliqueProjection")
{
    int32_t mode = io.index() % 4;
    Vector4<double> origin = Special4(io, mode, true);
    Vector4<double> normal = Special4(io, mode, false);
    Vector4<double> direction = Special4(io, mode, false);
    io.outMat(MakeObliqueProjection(origin, normal, direction));
}

// The perspective projection onto a plane. Arithmetic only. The last column
// of the matrix reads entries of the same row that were assigned earlier in
// the function, so the assignment order is part of the result.
ORACLE_CASE("Matrix4x4.makePerspectiveProjection")
{
    int32_t mode = io.index() % 4;
    Vector4<double> origin = Special4(io, mode, true);
    Vector4<double> normal = Special4(io, mode, false);
    Vector4<double> eye = Special4(io, mode, true);
    io.outMat(MakePerspectiveProjection(origin, normal, eye));
}

// The reflection through a plane. Arithmetic only. Six entries are copies of
// entries assigned earlier, so the assignment order is part of the result.
ORACLE_CASE("Matrix4x4.makeReflection")
{
    int32_t mode = io.index() % 4;
    Vector4<double> origin = Special4(io, mode, true);
    Vector4<double> normal = Special4(io, mode, false);
    io.outMat(MakeReflection(origin, normal));
}

// Transform::Identity, the static instance, next to a default-constructed
// transform. Arithmetic-free but it pins the documented initial state.
ORACLE_CASE("Transform.identity")
{
    Transform<double> constructed{};
    OutTransformFull(io, constructed);
    Transform<double> identity = Transform<double>::Identity();
    OutTransformFull(io, identity);
    // The implicit conversion 'operator Matrix4x4 const&' is GetHMatrix.
    Matrix4x4<double> const& h = identity;
    io.outMat(h);
}
