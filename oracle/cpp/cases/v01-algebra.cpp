// Verify group 1 (algebra): differential cases for Vector.h, Vector2.h,
// Vector3.h, Vector4.h, GVector.h, Matrix.h and GMatrix.h.
//
// This group is arithmetic only: every compared path is + - * / sqrt fabs and
// comparisons, so every case is declared exact (bit-identical) on the
// TypeScript side. sqrt is correctly rounded in both the MSVC runtime and V8,
// and Length/Normalize/UnitPerp/UnitCross are the only callers.
//
// Only Matrix.h, GMatrix.h, Vector.h, Vector2.h, Vector3.h, Vector4.h,
// GVector.h are included (through Oracle.h), so Inverse/Determinant resolve to
// the Gaussian-elimination templates of Matrix.h/GMatrix.h and not to the
// closed forms of Matrix2x2.h / Matrix3x3.h / Matrix4x4.h, which belong to
// group 2. Neither GTE_USE_COL_MAJOR nor GTE_USE_VEC_MAT is defined, so the
// row-major / MAT_VEC default that the port implements is the one measured.
//
// Generator modes. Every mode records exactly one double per component, so the
// TypeScript replay reads a fixed layout and never needs to know the mode:
//   0  small integer lattice: exact arithmetic, ties, exact zeros, parallel
//      and linearly dependent configurations, exactly singular matrices
//   1  uniform
//   2  wild magnitudes: signed zeros, subnormals on both sides of the
//      1/x-overflow threshold, values near the overflow and underflow limits
// Several cases add a fourth, constructed mode (linearly dependent vectors,
// collinear point sets, exactly singular and nearly singular matrices).
#define ORACLE_FAMILY "v01-algebra"
#include "Oracle.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

using namespace gte;

namespace
{
    // A value drawn from the "wild" population. Only unrecorded raw draws are
    // used here; the caller records the single resulting double.
    double WildRaw(oracle::Ctx& io)
    {
        double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        int32_t k = io.rawInteger(0, 11);
        double mantissa = io.raw(1.0, 2.0);
        switch (k)
        {
        case 0:  return sign * 0.0;                       // +-0
        case 1:  return sign * std::ldexp(1.0, -1074);    // smallest subnormal
        case 2:  return sign * std::ldexp(mantissa, -1060);  // 1/x overflows
        case 3:  return sign * std::ldexp(mantissa, -1023);  // 1/x is finite
        case 4:  return sign * std::ldexp(mantissa, -1022);  // smallest normal
        case 5:  return sign * std::ldexp(mantissa, 1020);   // near overflow
        case 6:  return sign * std::ldexp(mantissa, 500);
        case 7:  return sign * std::ldexp(mantissa, -500);
        case 8:  return sign * static_cast<double>(io.rawInteger(0, 4));
        case 9:  return sign * mantissa;
        case 10: return sign * std::ldexp(mantissa, 52);  // integer-valued
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
        return io.given(WildRaw(io));
    }

    template <int32_t N>
    Vector<N, double> Vec(oracle::Ctx& io, int32_t mode)
    {
        Vector<N, double> v;
        for (int32_t i = 0; i < N; ++i)
        {
            v[i] = Scalar(io, mode);
        }
        return v;
    }

    // A 4-tuple whose last component is exactly 0 when 'affine' is set (the
    // affine representation of a 3D vector). Records four doubles either way.
    Vector4<double> AffineVec4(oracle::Ctx& io, int32_t mode, bool affine)
    {
        Vector4<double> v;
        for (int32_t i = 0; i < 3; ++i)
        {
            v[i] = Scalar(io, mode);
        }
        v[3] = (affine ? io.given(0.0) : Scalar(io, mode));
        return v;
    }

    GVector<double> GVec(oracle::Ctx& io, int32_t size, int32_t mode)
    {
        GVector<double> v(size);
        for (int32_t i = 0; i < size; ++i)
        {
            v[i] = Scalar(io, mode);
        }
        return v;
    }

    // A nonzero scalar. GVector's and GMatrix's operator/= call LogError for a
    // zero divisor (unlike Vector.h's and Matrix.h's, which yield the zero
    // tuple), and the port shares Vector.h's division, so the GVector and
    // GMatrix division cases stay away from zero. Records one double.
    double NonzeroScalar(oracle::Ctx& io, int32_t mode)
    {
        double s = 0.0;
        do
        {
            if (mode == 0)
            {
                s = static_cast<double>(io.rawInteger(-3, 3));
            }
            else if (mode == 1)
            {
                s = io.raw(-10.0, 10.0);
            }
            else
            {
                s = WildRaw(io);
            }
        }
        while (s == 0.0);
        return io.given(s);
    }

    template <int32_t R, int32_t C>
    Matrix<R, C, double> Mat(oracle::Ctx& io, int32_t mode)
    {
        Matrix<R, C, double> m;
        for (int32_t r = 0; r < R; ++r)
        {
            for (int32_t c = 0; c < C; ++c)
            {
                m(r, c) = Scalar(io, mode);
            }
        }
        return m;
    }

    GMatrix<double> GMat(oracle::Ctx& io, int32_t numRows, int32_t numCols, int32_t mode)
    {
        GMatrix<double> m(numRows, numCols);
        for (int32_t r = 0; r < numRows; ++r)
        {
            for (int32_t c = 0; c < numCols; ++c)
            {
                m(r, c) = Scalar(io, mode);
            }
        }
        return m;
    }

    // A unit-length 3-tuple. Mode 0 is an exact coordinate axis (so
    // FastComputeOrthogonalComplement sees v2[2] exactly 0, +1 or -1), mode 1
    // is a normalized lattice direction and mode 2 is a uniform unit vector.
    // All three record exactly three doubles.
    Vector3<double> UnitVec3(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 0)
        {
            int32_t d = io.rawInteger(0, 2);
            double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
            Vector3<double> v{ 0.0, 0.0, 0.0 };
            v[d] = sign;
            return io.givenVec(v);
        }
        if (mode == 1)
        {
            Vector3<double> v{ 0.0, 0.0, 0.0 };
            double len = 0.0;
            do
            {
                for (int32_t i = 0; i < 3; ++i)
                {
                    v[i] = static_cast<double>(io.rawInteger(-4, 4));
                }
                len = Length(v);
            }
            while (len == 0.0);
            Normalize(v);
            return io.givenVec(v);
        }
        return io.unit<3>();
    }

    // A square matrix in the requested mode. Modes 0 and 1 are the shared
    // lattice and uniform populations, mode 2 builds an exactly singular
    // matrix (the last row is an integer combination of the others), mode 3
    // perturbs that combination by a single tiny power of two so that full
    // pivoting ends on a minute (sometimes subnormal) pivot, and mode 4 is the
    // wild population, which reaches the denormal-pivot NaN of
    // GaussianElimination (upstream issue #375, preserved by the port).
    template <int32_t N>
    Matrix<N, N, double> SquareMat(oracle::Ctx& io, int32_t mode)
    {
        Matrix<N, N, double> m;
        if (mode != 2 && mode != 3)
        {
            for (int32_t r = 0; r < N; ++r)
            {
                for (int32_t c = 0; c < N; ++c)
                {
                    m(r, c) = Scalar(io, mode == 4 ? 2 : mode);
                }
            }
            return m;
        }

        for (int32_t r = 0; r + 1 < N; ++r)
        {
            for (int32_t c = 0; c < N; ++c)
            {
                m(r, c) = static_cast<double>(io.rawInteger(-4, 4));
            }
        }
        for (int32_t c = 0; c < N; ++c)
        {
            m(N - 1, c) = 0.0;
        }
        for (int32_t r = 0; r + 1 < N; ++r)
        {
            double k = static_cast<double>(io.rawInteger(-2, 2));
            for (int32_t c = 0; c < N; ++c)
            {
                m(N - 1, c) += k * m(r, c);
            }
        }
        if (mode == 3)
        {
            int32_t c = io.rawInteger(0, N - 1);
            int32_t e = io.rawInteger(-1070, -40);
            m(N - 1, c) += std::ldexp(1.0, e);
        }
        for (int32_t r = 0; r < N; ++r)
        {
            for (int32_t c = 0; c < N; ++c)
            {
                io.given(m(r, c));
            }
        }
        return m;
    }
}

// ---------------------------------------------------------------------------
// Vector.h
// ---------------------------------------------------------------------------

ORACLE_CASE("Vector.arithmetic")
{
    // Unary operator- and operator+, operator+/-, operator* and operator/ by a
    // scalar (including the zero divisor, which upstream special-cases to the
    // zero vector), and the componentwise operator* / operator/.
    int32_t mode = io.index() % 3;
    auto v0 = Vec<4>(io, mode);
    auto v1 = Vec<4>(io, mode);
    double s = Scalar(io, mode);
    io.outVec(-v0);
    io.outVec(+v0);
    io.outVec(v0 + v1);
    io.outVec(v0 - v1);
    io.outVec(v0 * s);
    io.outVec(s * v0);
    io.outVec(v0 / s);
    io.outVec(v0 * v1);
    io.outVec(v0 / v1);
}

ORACLE_CASE("Vector.dot")
{
    int32_t mode = io.index() % 3;
    auto v0 = Vec<4>(io, mode);
    auto v1 = Vec<4>(io, mode);
    io.outReal(Dot(v0, v1));
}

ORACLE_CASE("Vector.dot.signedZero")
{
    // Every product v0[i]*v1[i] is a signed zero, so the seed of the
    // accumulation decides the sign of the result: Vector.h's Dot starts from
    // v0[0]*v1[0] (all-negative products give -0) while GVector.h's starts
    // from the literal 0, which turns a leading -0 into +0.
    Vector<4, double> v0, v1;
    for (int32_t i = 0; i < 4; ++i)
    {
        double zero = (io.rawInteger(0, 1) == 0 ? 0.0 : -0.0);
        v0[i] = io.given(zero);
    }
    for (int32_t i = 0; i < 4; ++i)
    {
        double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        double magnitude = static_cast<double>(io.rawInteger(1, 3));
        v1[i] = io.given(sign * magnitude);
    }
    io.outReal(Dot(v0, v1));
    io.outReal(Dot(v1, v0));
}

ORACLE_CASE("Vector.length")
{
    int32_t mode = io.index() % 3;
    auto v = Vec<4>(io, mode);
    io.outReal(Length(v, false));
    io.outReal(Length(v, true));
}

ORACLE_CASE("Vector.length.subnormal")
{
    // Upstream issue #370, preserved by the port: the robust paths of Length
    // and Normalize rescale with 'v /= maxAbsComp', which Vector.h implements
    // as a multiplication by 1/maxAbsComp. When the largest component is
    // subnormal and below about 2^-1024 that reciprocal overflows to infinity
    // and the zero components become 0 * inf = NaN, so the robust length is
    // NaN while the plain one is correct. Between 2^-1074 and 2^-1024 the
    // reciprocal is still finite and both paths agree; both regimes occur
    // here, as does the all-zero vector.
    Vector<4, double> v;
    for (int32_t i = 0; i < 4; ++i)
    {
        int32_t kind = io.rawInteger(0, 3);
        double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        double mantissa = io.raw(1.0, 2.0);
        int32_t exponent = io.rawInteger(-1074, -1000);
        double x = (kind == 0 ? sign * 0.0 : sign * std::ldexp(mantissa, exponent));
        v[i] = io.given(x);
    }
    Vector<4, double> plain = v;
    Vector<4, double> robust = v;
    io.outReal(Length(v, false));
    io.outReal(Length(v, true));
    io.outReal(Normalize(plain, false));
    io.outVec(plain);
    io.outReal(Normalize(robust, true));
    io.outVec(robust);
}

ORACLE_CASE("Vector.normalize")
{
    int32_t mode = io.index() % 3;
    auto v = Vec<4>(io, mode);
    Vector<4, double> plain = v;
    Vector<4, double> robust = v;
    double lengthPlain = Normalize(plain, false);
    double lengthRobust = Normalize(robust, true);
    io.outReal(lengthPlain);
    io.outVec(plain);
    io.outReal(lengthRobust);
    io.outVec(robust);
}

ORACLE_CASE("Vector.orthonormalize")
{
    // Mode 3 makes the inputs exactly linearly dependent, so Gram-Schmidt
    // produces a zero-length vector and the zero-length branch of Normalize
    // fires with minLength exactly 0.
    int32_t numInputs = io.integer(1, 4);
    bool robust = io.boolean();
    int32_t mode = io.index() % 4;
    std::array<Vector<4, double>, 4> v{};
    if (mode == 3)
    {
        for (int32_t i = 0; i < 4; ++i)
        {
            Vector<4, double> w;
            w.MakeZero();
            if (i == 0)
            {
                for (int32_t k = 0; k < 4; ++k)
                {
                    w[k] = static_cast<double>(io.rawInteger(-3, 3));
                }
            }
            else
            {
                for (int32_t j = 0; j < i; ++j)
                {
                    double k = static_cast<double>(io.rawInteger(-2, 2));
                    w = w + v[j] * k;
                }
            }
            v[i] = io.givenVec(w);
        }
    }
    else
    {
        for (int32_t i = 0; i < 4; ++i)
        {
            v[i] = Vec<4>(io, mode);
        }
    }
    double minLength = Orthonormalize<4, double>(numInputs, v.data(), robust);
    io.outReal(minLength);
    for (int32_t i = 0; i < 4; ++i)
    {
        io.outVec(v[i]);
    }
}

ORACLE_CASE("Vector.getOrthogonal")
{
    int32_t n = io.integer(2, 4);
    bool unitLength = io.boolean();
    int32_t mode = io.index() % 3;
    if (n == 2)
    {
        auto v = Vec<2>(io, mode);
        io.outVec(GetOrthogonal(v, unitLength));
    }
    else if (n == 3)
    {
        auto v = Vec<3>(io, mode);
        io.outVec(GetOrthogonal(v, unitLength));
    }
    else
    {
        auto v = Vec<4>(io, mode);
        io.outVec(GetOrthogonal(v, unitLength));
    }
}

ORACLE_CASE("Vector.computeExtremes")
{
    // numVectors == 0 is the "invalid input" arm: the return value is false
    // and vmin/vmax keep the zero values they were given.
    int32_t numVectors = io.integer(0, 5);
    int32_t mode = io.index() % 3;
    std::vector<Vector<3, double>> v(static_cast<size_t>(numVectors));
    for (int32_t i = 0; i < numVectors; ++i)
    {
        v[static_cast<size_t>(i)] = Vec<3>(io, mode);
    }
    Vector<3, double> vmin, vmax;
    vmin.MakeZero();
    vmax.MakeZero();
    bool valid = ComputeExtremes<3, double>(numVectors, v.data(), vmin, vmax);
    io.outBool(valid);
    io.outVec(vmin);
    io.outVec(vmax);
}

ORACLE_CASE("Vector.liftProject")
{
    int32_t inject = io.integer(0, 4);
    int32_t reject = io.integer(0, 3);
    int32_t mode = io.index() % 3;
    double last = Scalar(io, mode);
    auto v = Vec<4>(io, mode);
    io.outVec(HLift(v, last));
    io.outVec(HProject(v));
    io.outVec(Lift(v, inject, last));
    io.outVec(Project(v, reject));
}

ORACLE_CASE("Vector.special")
{
    // MakeZero/MakeOnes/MakeUnit and the Zero/Ones/Unit factories, plus the
    // six comparison operators (std::array lexicographic order).
    int32_t d = io.integer(-1, 4);
    int32_t mode = io.index() % 3;
    auto v0 = Vec<4>(io, mode);
    auto v1 = Vec<4>(io, mode);
    Vector<4, double> zero, ones, unit;
    zero.MakeZero();
    ones.MakeOnes();
    unit.MakeUnit(d);
    io.outVec(zero);
    io.outVec(ones);
    io.outVec(unit);
    io.outVec(Vector<4, double>::Zero());
    io.outVec(Vector<4, double>::Ones());
    io.outVec(Vector<4, double>::Unit(d));
    io.outInt(v0.GetSize());
    io.outBool(v0 == v1);
    io.outBool(v0 != v1);
    io.outBool(v0 < v1);
    io.outBool(v0 <= v1);
    io.outBool(v0 > v1);
    io.outBool(v0 >= v1);
    io.outBool(v0 == v0);
    io.outBool(v0 <= v0);
}

// ---------------------------------------------------------------------------
// Vector2.h
// ---------------------------------------------------------------------------

ORACLE_CASE("Vector2.perp")
{
    int32_t mode = io.index() % 3;
    bool robust = io.boolean();
    auto v0 = Vec<2>(io, mode);
    auto v1 = Vec<2>(io, mode);
    io.outVec(Perp(v0));
    io.outVec(UnitPerp(v0, robust));
    io.outReal(DotPerp(v0, v1));
    io.outReal(DotPerp(v1, v0));
}

ORACLE_CASE("Vector2.computeOrthogonalComplement")
{
    // numInputs == 0 exercises the "invalid input" arm, which returns 0 and
    // leaves both vectors untouched.
    int32_t numInputs = io.integer(0, 1);
    bool robust = io.boolean();
    int32_t mode = io.index() % 3;
    std::array<Vector2<double>, 2> v{};
    v[0] = Vec<2>(io, mode);
    v[1] = Vec<2>(io, mode);
    double minLength = ComputeOrthogonalComplement<double>(numInputs, v.data(), robust);
    io.outReal(minLength);
    io.outVec(v[0]);
    io.outVec(v[1]);
}

ORACLE_CASE("Vector2.computeBarycentrics")
{
    // Mode 3 puts v2 and p exactly on the line through v0 and v1, so the
    // determinant DotPerp(diff[0], diff[1]) is exactly zero and the
    // "linearly dependent" arm fires even with epsilon = 0.
    int32_t epsilonMode = io.index() % 3;
    double epsilon = (epsilonMode == 0 ? io.given(0.0)
        : (epsilonMode == 1 ? io.real(0.0, 0.5) : io.given(-1.0)));
    int32_t mode = io.index() % 4;
    Vector2<double> v0, v1, v2, p;
    if (mode == 3)
    {
        Vector2<double> a, b;
        for (int32_t k = 0; k < 2; ++k)
        {
            a[k] = static_cast<double>(io.rawInteger(-4, 4));
        }
        for (int32_t k = 0; k < 2; ++k)
        {
            b[k] = static_cast<double>(io.rawInteger(-4, 4));
        }
        double t2 = static_cast<double>(io.rawInteger(-3, 3));
        double tp = static_cast<double>(io.rawInteger(-3, 3));
        v0 = io.givenVec(a);
        v1 = io.givenVec(b);
        v2 = io.givenVec(a + (b - a) * t2);
        p = io.givenVec(a + (b - a) * tp);
    }
    else
    {
        v0 = Vec<2>(io, mode);
        v1 = Vec<2>(io, mode);
        v2 = Vec<2>(io, mode);
        p = Vec<2>(io, mode);
    }
    std::array<double, 3> bary{};
    bool valid = ComputeBarycentrics(p, v0, v1, v2, bary, epsilon);
    io.outBool(valid);
    for (int32_t i = 0; i < 3; ++i)
    {
        io.outReal(bary[i]);
    }
}

ORACLE_CASE("Vector2.intrinsicsVector2")
{
    // Mode 2 places every point on one lattice line, so the dimension-1 arm
    // fires with maxDistance exactly 0. A negative epsilon exercises the
    // "invalid input" arm, in which every member keeps its zero default.
    int32_t numVectors = io.integer(1, 6);
    int32_t epsilonMode = io.index() % 3;
    double epsilon = (epsilonMode == 0 ? io.given(0.0)
        : (epsilonMode == 1 ? io.real(0.0, 0.25) : io.given(-1.0)));
    int32_t mode = io.index() % 3;
    std::vector<Vector2<double>> v(static_cast<size_t>(numVectors));
    if (mode == 2)
    {
        Vector2<double> a, d;
        for (int32_t k = 0; k < 2; ++k)
        {
            a[k] = static_cast<double>(io.rawInteger(-4, 4));
        }
        for (int32_t k = 0; k < 2; ++k)
        {
            d[k] = static_cast<double>(io.rawInteger(-3, 3));
        }
        for (int32_t i = 0; i < numVectors; ++i)
        {
            double t = static_cast<double>(io.rawInteger(-4, 4));
            v[static_cast<size_t>(i)] = io.givenVec(a + d * t);
        }
    }
    else
    {
        for (int32_t i = 0; i < numVectors; ++i)
        {
            v[static_cast<size_t>(i)] = Vec<2>(io, mode);
        }
    }
    IntrinsicsVector2<double> intr(numVectors, v.data(), epsilon);
    io.outReal(intr.epsilon);
    io.outInt(intr.dimension);
    io.outReal(intr.min[0]);
    io.outReal(intr.min[1]);
    io.outReal(intr.max[0]);
    io.outReal(intr.max[1]);
    io.outReal(intr.maxRange);
    io.outVec(intr.origin);
    io.outVec(intr.direction[0]);
    io.outVec(intr.direction[1]);
    io.outInt(intr.extreme[0]);
    io.outInt(intr.extreme[1]);
    io.outInt(intr.extreme[2]);
    io.outBool(intr.extremeCCW);
}

// ---------------------------------------------------------------------------
// Vector3.h
// ---------------------------------------------------------------------------

ORACLE_CASE("Vector3.cross")
{
    int32_t mode = io.index() % 3;
    bool robust = io.boolean();
    auto v0 = Vec<3>(io, mode);
    auto v1 = Vec<3>(io, mode);
    auto v2 = Vec<3>(io, mode);
    io.outVec(Cross(v0, v1));
    io.outVec(UnitCross(v0, v1, robust));
    io.outReal(DotCross(v0, v1, v2));
    io.outReal(DotCross(v2, v0, v1));
}

ORACLE_CASE("Vector3.cross.4d")
{
    // Cross/UnitCross/DotCross also accept 4-tuples (affine vectors with
    // w = 0); the fourth component of the result is always 0.
    int32_t mode = io.index() % 3;
    bool robust = io.boolean();
    bool affine = io.boolean();
    auto v0 = AffineVec4(io, mode, affine);
    auto v1 = AffineVec4(io, mode, affine);
    auto v2 = AffineVec4(io, mode, affine);
    io.outVec(Cross(v0, v1));
    io.outVec(UnitCross(v0, v1, robust));
    io.outReal(DotCross(v0, v1, v2));
}

ORACLE_CASE("Vector3.computeOrthogonalComplement")
{
    // numInputs == 0 is the "invalid input" arm. Mode 3 makes v[1] an exact
    // integer multiple of v[0], so the Gram-Schmidt step produces a
    // zero-length vector.
    int32_t numInputs = io.integer(0, 2);
    bool robust = io.boolean();
    int32_t mode = io.index() % 4;
    std::array<Vector3<double>, 3> v{};
    if (mode == 3)
    {
        Vector3<double> a;
        for (int32_t k = 0; k < 3; ++k)
        {
            a[k] = static_cast<double>(io.rawInteger(-3, 3));
        }
        double scale = static_cast<double>(io.rawInteger(-3, 3));
        v[0] = io.givenVec(a);
        v[1] = io.givenVec(a * scale);
        v[2] = io.givenVec(a * scale);
    }
    else
    {
        v[0] = Vec<3>(io, mode);
        v[1] = Vec<3>(io, mode);
        v[2] = Vec<3>(io, mode);
    }
    double minLength = ComputeOrthogonalComplement<double>(numInputs, v.data(), robust);
    io.outReal(minLength);
    for (int32_t i = 0; i < 3; ++i)
    {
        io.outVec(v[i]);
    }
}

ORACLE_CASE("Vector3.fastComputeOrthogonalComplement")
{
    // Upstream documents v2 as unit length. Mode 0 supplies exact coordinate
    // axes, so v2[2] is exactly 0, +1 or -1 and both arms of the v2[2] >= 0
    // test are taken at the boundary.
    int32_t mode = io.index() % 3;
    auto v2 = UnitVec3(io, mode);
    Vector3<double> v0, v1;
    v0.MakeZero();
    v1.MakeZero();
    FastComputeOrthogonalComplement(v2, v0, v1);
    io.outVec(v0);
    io.outVec(v1);
}

ORACLE_CASE("Vector3.computeBarycentrics")
{
    // Mode 3 puts v3 and p in the plane of v0, v1, v2 on a lattice, so the
    // determinant DotCross(diff[0], diff[1], diff[2]) is exactly zero.
    int32_t epsilonMode = io.index() % 3;
    double epsilon = (epsilonMode == 0 ? io.given(0.0)
        : (epsilonMode == 1 ? io.real(0.0, 0.5) : io.given(-1.0)));
    int32_t mode = io.index() % 4;
    Vector3<double> v0, v1, v2, v3, p;
    if (mode == 3)
    {
        Vector3<double> a, e0, e1;
        for (int32_t k = 0; k < 3; ++k)
        {
            a[k] = static_cast<double>(io.rawInteger(-3, 3));
        }
        for (int32_t k = 0; k < 3; ++k)
        {
            e0[k] = static_cast<double>(io.rawInteger(-3, 3));
        }
        for (int32_t k = 0; k < 3; ++k)
        {
            e1[k] = static_cast<double>(io.rawInteger(-3, 3));
        }
        double s1 = static_cast<double>(io.rawInteger(-2, 2));
        double t1 = static_cast<double>(io.rawInteger(-2, 2));
        double s2 = static_cast<double>(io.rawInteger(-2, 2));
        double t2 = static_cast<double>(io.rawInteger(-2, 2));
        double s3 = static_cast<double>(io.rawInteger(-2, 2));
        double t3 = static_cast<double>(io.rawInteger(-2, 2));
        double sp = static_cast<double>(io.rawInteger(-2, 2));
        double tp = static_cast<double>(io.rawInteger(-2, 2));
        v0 = io.givenVec(a);
        v1 = io.givenVec(a + e0 * s1 + e1 * t1);
        v2 = io.givenVec(a + e0 * s2 + e1 * t2);
        v3 = io.givenVec(a + e0 * s3 + e1 * t3);
        p = io.givenVec(a + e0 * sp + e1 * tp);
    }
    else
    {
        v0 = Vec<3>(io, mode);
        v1 = Vec<3>(io, mode);
        v2 = Vec<3>(io, mode);
        v3 = Vec<3>(io, mode);
        p = Vec<3>(io, mode);
    }
    std::array<double, 4> bary{};
    bool valid = ComputeBarycentrics(p, v0, v1, v2, v3, bary, epsilon);
    io.outBool(valid);
    for (int32_t i = 0; i < 4; ++i)
    {
        io.outReal(bary[i]);
    }
}

ORACLE_CASE("Vector3.intrinsicsVector3")
{
    // Mode 2 places the points on one lattice line (dimension 1) and mode 3
    // in one lattice plane (dimension 2), both with maxDistance exactly 0.
    int32_t numVectors = io.integer(1, 6);
    int32_t epsilonMode = io.index() % 3;
    double epsilon = (epsilonMode == 0 ? io.given(0.0)
        : (epsilonMode == 1 ? io.real(0.0, 0.25) : io.given(-1.0)));
    int32_t mode = io.index() % 4;
    std::vector<Vector3<double>> v(static_cast<size_t>(numVectors));
    if (mode == 2 || mode == 3)
    {
        Vector3<double> a, d0, d1;
        for (int32_t k = 0; k < 3; ++k)
        {
            a[k] = static_cast<double>(io.rawInteger(-4, 4));
        }
        for (int32_t k = 0; k < 3; ++k)
        {
            d0[k] = static_cast<double>(io.rawInteger(-3, 3));
        }
        for (int32_t k = 0; k < 3; ++k)
        {
            d1[k] = (mode == 3 ? static_cast<double>(io.rawInteger(-3, 3)) : 0.0);
        }
        for (int32_t i = 0; i < numVectors; ++i)
        {
            double s = static_cast<double>(io.rawInteger(-4, 4));
            double t = static_cast<double>(io.rawInteger(-4, 4));
            v[static_cast<size_t>(i)] = io.givenVec(a + d0 * s + d1 * t);
        }
    }
    else
    {
        for (int32_t i = 0; i < numVectors; ++i)
        {
            v[static_cast<size_t>(i)] = Vec<3>(io, mode);
        }
    }
    IntrinsicsVector3<double> intr(numVectors, v.data(), epsilon);
    io.outReal(intr.epsilon);
    io.outInt(intr.dimension);
    for (int32_t k = 0; k < 3; ++k)
    {
        io.outReal(intr.min[k]);
    }
    for (int32_t k = 0; k < 3; ++k)
    {
        io.outReal(intr.max[k]);
    }
    io.outReal(intr.maxRange);
    io.outVec(intr.origin);
    io.outVec(intr.direction[0]);
    io.outVec(intr.direction[1]);
    io.outVec(intr.direction[2]);
    for (int32_t k = 0; k < 4; ++k)
    {
        io.outInt(intr.extreme[k]);
    }
    io.outBool(intr.extremeCCW);
}

// ---------------------------------------------------------------------------
// Vector4.h
// ---------------------------------------------------------------------------

ORACLE_CASE("Vector4.hyperCross")
{
    int32_t mode = io.index() % 3;
    bool robust = io.boolean();
    auto v0 = Vec<4>(io, mode);
    auto v1 = Vec<4>(io, mode);
    auto v2 = Vec<4>(io, mode);
    auto v3 = Vec<4>(io, mode);
    io.outVec(HyperCross(v0, v1, v2));
    io.outVec(UnitHyperCross(v0, v1, v2, robust));
    io.outReal(DotHyperCross(v0, v1, v2, v3));
    io.outReal(DotHyperCross(v3, v0, v1, v2));
}

ORACLE_CASE("Vector4.computeOrthogonalComplement")
{
    // numInputs == 0 is the "invalid input" arm. Mode 3 supplies affine
    // vectors (w = 0), the documented use case of the maxIndex == 3 branch.
    int32_t numInputs = io.integer(0, 3);
    bool robust = io.boolean();
    int32_t mode = io.index() % 4;
    std::array<Vector4<double>, 4> v{};
    for (int32_t i = 0; i < 4; ++i)
    {
        v[static_cast<size_t>(i)] = AffineVec4(io, mode == 3 ? 0 : mode, mode == 3);
    }
    double minLength = ComputeOrthogonalComplement<double>(numInputs, v.data(), robust);
    io.outReal(minLength);
    for (int32_t i = 0; i < 4; ++i)
    {
        io.outVec(v[static_cast<size_t>(i)]);
    }
}

ORACLE_CASE("Vector4.computeOrthogonalComplement.zeroMiddle")
{
    // Upstream defect #87 (preserved by the port): for v[0] = (a, 0, 0, b)
    // with |b| > |a| the maxIndex == 3 branch sets v[1] = (0, +v[0][2],
    // -v[0][1], 0) = 0, so Orthonormalize normalizes the zero vector and the
    // returned minimum length is 0. The port preserves the defect, so this
    // case compares bit for bit rather than being a deviation.
    bool robust = io.boolean();
    double a = static_cast<double>(io.rawInteger(-3, 3));
    double b = static_cast<double>(io.rawInteger(4, 9));
    double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
    std::array<Vector4<double>, 4> v{};
    v[0] = io.givenVec(Vector4<double>{ a, 0.0, 0.0, sign * b });
    v[1] = io.givenVec(Vector4<double>{ 0.0, 0.0, 0.0, 0.0 });
    v[2] = io.givenVec(Vector4<double>{ 0.0, 0.0, 0.0, 0.0 });
    v[3] = io.givenVec(Vector4<double>{ 0.0, 0.0, 0.0, 0.0 });
    double minLength = ComputeOrthogonalComplement<double>(1, v.data(), robust);
    io.outReal(minLength);
    for (int32_t i = 0; i < 4; ++i)
    {
        io.outVec(v[static_cast<size_t>(i)]);
    }
}

// ---------------------------------------------------------------------------
// GVector.h
//
// The port merges GVector into Vector (GVector extends Vector and reuses its
// free functions), so these cases measure the shared implementation against
// upstream's separate GVector.h one. Four documented merge deviations are kept
// out of the generators and listed in the group report: division by zero,
// HProject/Project of a size-1 tuple, Dot/Length of an empty tuple, and the
// LogError arms of Orthonormalize/ComputeExtremes for invalid input.
// ---------------------------------------------------------------------------

ORACLE_CASE("GVector.arithmetic")
{
    int32_t size = io.integer(1, 5);
    int32_t mode = io.index() % 3;
    auto v0 = GVec(io, size, mode);
    auto v1 = GVec(io, size, mode);
    double s = NonzeroScalar(io, mode);
    io.outGVec(-v0);
    io.outGVec(+v0);
    io.outGVec(v0 + v1);
    io.outGVec(v0 - v1);
    io.outGVec(v0 * s);
    io.outGVec(s * v0);
    io.outGVec(v0 / s);
}

ORACLE_CASE("GVector.dot")
{
    // GVector.h's Dot seeds the accumulation with the literal 0 where
    // Vector.h's seeds it with v0[0]*v1[0]; the two differ only when every
    // product is -0, which the signed zeros of the wild mode reach.
    int32_t size = io.integer(1, 5);
    int32_t mode = io.index() % 3;
    auto v0 = GVec(io, size, mode);
    auto v1 = GVec(io, size, mode);
    io.outReal(Dot(v0, v1));
}

ORACLE_CASE("GVector.dot.signedZero")
{
    // Every product is a signed zero, so the accumulation seed decides the
    // sign of the result. Aimed squarely at the difference above.
    int32_t size = io.integer(1, 4);
    GVector<double> v0(size), v1(size);
    for (int32_t i = 0; i < size; ++i)
    {
        double zero = (io.rawInteger(0, 1) == 0 ? 0.0 : -0.0);
        v0[i] = io.given(zero);
    }
    for (int32_t i = 0; i < size; ++i)
    {
        double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
        double magnitude = static_cast<double>(io.rawInteger(1, 3));
        v1[i] = io.given(sign * magnitude);
    }
    io.outReal(Dot(v0, v1));
    io.outReal(Dot(v1, v0));
}

ORACLE_CASE("GVector.lengthNormalize")
{
    int32_t size = io.integer(1, 5);
    int32_t mode = io.index() % 3;
    auto v = GVec(io, size, mode);
    GVector<double> plain = v;
    GVector<double> robust = v;
    io.outReal(Length(v, false));
    io.outReal(Length(v, true));
    double lengthPlain = Normalize(plain, false);
    double lengthRobust = Normalize(robust, true);
    io.outReal(lengthPlain);
    io.outGVec(plain);
    io.outReal(lengthRobust);
    io.outGVec(robust);
}

ORACLE_CASE("GVector.orthonormalize")
{
    int32_t size = io.integer(2, 5);
    int32_t numInputs = io.integer(1, size);
    bool robust = io.boolean();
    int32_t mode = io.index() % 4;
    std::vector<GVector<double>> v(static_cast<size_t>(size));
    if (mode == 3)
    {
        // Exactly linearly dependent inputs.
        for (int32_t i = 0; i < size; ++i)
        {
            GVector<double> w(size);
            w.MakeZero();
            if (i == 0)
            {
                for (int32_t k = 0; k < size; ++k)
                {
                    w[k] = static_cast<double>(io.rawInteger(-3, 3));
                }
            }
            else
            {
                for (int32_t j = 0; j < i; ++j)
                {
                    double k = static_cast<double>(io.rawInteger(-2, 2));
                    w = w + v[static_cast<size_t>(j)] * k;
                }
            }
            for (int32_t k = 0; k < size; ++k)
            {
                io.given(w[k]);
            }
            v[static_cast<size_t>(i)] = w;
        }
    }
    else
    {
        for (int32_t i = 0; i < size; ++i)
        {
            v[static_cast<size_t>(i)] = GVec(io, size, mode);
        }
    }
    double minLength = Orthonormalize<double>(numInputs, v.data(), robust);
    io.outReal(minLength);
    for (int32_t i = 0; i < size; ++i)
    {
        io.outGVec(v[static_cast<size_t>(i)]);
    }
}

ORACLE_CASE("GVector.computeExtremes")
{
    int32_t size = io.integer(1, 4);
    int32_t numVectors = io.integer(1, 5);
    int32_t mode = io.index() % 3;
    std::vector<GVector<double>> v(static_cast<size_t>(numVectors));
    for (int32_t i = 0; i < numVectors; ++i)
    {
        v[static_cast<size_t>(i)] = GVec(io, size, mode);
    }
    GVector<double> vmin(size), vmax(size);
    bool valid = ComputeExtremes<double>(numVectors, v.data(), vmin, vmax);
    io.outBool(valid);
    io.outGVec(vmin);
    io.outGVec(vmax);
}

ORACLE_CASE("GVector.liftProject")
{
    int32_t size = io.integer(2, 5);
    int32_t inject = io.integer(0, size);
    int32_t reject = io.integer(0, size - 1);
    int32_t mode = io.index() % 3;
    double last = Scalar(io, mode);
    auto v = GVec(io, size, mode);
    io.outGVec(HLift(v, last));
    io.outGVec(HProject(v));
    io.outGVec(Lift(v, inject, last));
    io.outGVec(Project(v, reject));
}

ORACLE_CASE("GVector.special")
{
    // MakeZero/MakeUnit, the Zero/Unit factories, SetSize (std::vector::resize
    // semantics: the leading elements survive, new ones are value-initialized)
    // and the six std::vector comparison operators, which unlike Vector's
    // accept tuples of different sizes.
    int32_t size0 = io.integer(1, 5);
    int32_t size1 = io.integer(1, 5);
    int32_t d = io.integer(-1, 5);
    int32_t newSize = io.integer(0, 7);
    int32_t mode = io.index() % 3;
    auto v0 = GVec(io, size0, mode);
    auto v1 = GVec(io, size1, mode);
    GVector<double> zero(size0), unit(size0);
    zero.MakeZero();
    unit.MakeUnit(d);
    io.outGVec(zero);
    io.outGVec(unit);
    io.outGVec(GVector<double>::Zero(size0));
    io.outGVec(GVector<double>::Unit(size0, d));
    io.outInt(v0.GetSize());
    io.outBool(v0 == v1);
    io.outBool(v0 != v1);
    io.outBool(v0 < v1);
    io.outBool(v0 <= v1);
    io.outBool(v0 > v1);
    io.outBool(v0 >= v1);
    GVector<double> resized = v0;
    resized.SetSize(newSize);
    io.outGVec(resized);
}

// ---------------------------------------------------------------------------
// Matrix.h
// ---------------------------------------------------------------------------

namespace
{
    template <int32_t R, int32_t C>
    void MatArithmetic(oracle::Ctx& io, int32_t mode)
    {
        auto m0 = Mat<R, C>(io, mode);
        auto m1 = Mat<R, C>(io, mode);
        double s = Scalar(io, mode);
        io.outMat(-m0);
        io.outMat(+m0);
        io.outMat(m0 + m1);
        io.outMat(m0 - m1);
        io.outMat(m0 * s);
        io.outMat(s * m0);
        io.outMat(m0 / s);
        io.outReal(L1Norm(m0));
        io.outReal(L2Norm(m0));
        io.outReal(LInfinityNorm(m0));
    }

    template <int32_t R, int32_t K, int32_t C>
    void MatMultiply(oracle::Ctx& io, int32_t mode)
    {
        auto A = Mat<R, K>(io, mode);
        auto B = Mat<K, C>(io, mode);
        auto Bt = Mat<C, K>(io, mode);
        auto At = Mat<K, R>(io, mode);
        auto v = Vec<K>(io, mode);
        auto w = Vec<R>(io, mode);
        io.outMat(A * B);
        io.outMat(MultiplyAB(A, B));
        io.outMat(MultiplyABT(A, Bt));
        io.outMat(MultiplyATB(At, B));
        io.outMat(MultiplyATBT(At, Bt));
        io.outVec(A * v);
        io.outVec(w * A);
        io.outMat(Transpose(A));
    }

    template <int32_t R, int32_t C>
    void MatDiagonal(oracle::Ctx& io, int32_t mode)
    {
        auto M = Mat<R, C>(io, mode);
        auto dc = Vec<C>(io, mode);
        auto dr = Vec<R>(io, mode);
        io.outMat(MultiplyMD(M, dc));
        io.outMat(MultiplyDM(dr, M));
        io.outMat(OuterProduct(dr, dc));
        Matrix<C, C, double> diag;
        MakeDiagonal(dc, diag);
        io.outMat(diag);
    }

    template <int32_t N>
    void MatInverse(oracle::Ctx& io, int32_t mode)
    {
        auto M = SquareMat<N>(io, mode);
        bool invertible = false;
        Matrix<N, N, double> invM = Inverse(M, &invertible);
        io.outBool(invertible);
        io.outMat(invM);
        io.outReal(Determinant(M));
        io.outMat(Inverse(M));
    }

    // Upstream issue #375, preserved by the port: for a matrix whose entries
    // are all of denormal magnitude the pivot search finds a nonzero pivot,
    // but 1/pivot overflows to infinity and the row operations produce
    // inf * 0 = NaN, so Inverse returns NaN entries while reporting
    // invertible = true and Determinant can be NaN as well.
    template <int32_t N>
    void DenormalInverse(oracle::Ctx& io)
    {
        Matrix<N, N, double> M;
        for (int32_t r = 0; r < N; ++r)
        {
            for (int32_t c = 0; c < N; ++c)
            {
                double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
                double mantissa = io.raw(1.0, 2.0);
                int32_t exponent = io.rawInteger(-1074, -1023);
                M(r, c) = io.given(sign * std::ldexp(mantissa, exponent));
            }
        }
        bool invertible = false;
        Matrix<N, N, double> invM = Inverse(M, &invertible);
        io.outBool(invertible);
        io.outMat(invM);
        io.outReal(Determinant(M));
    }

    template <int32_t R, int32_t C>
    void MatAccess(oracle::Ctx& io, int32_t mode)
    {
        int32_t r = io.integer(0, R - 1);
        int32_t c = io.integer(0, C - 1);
        int32_t ur = io.integer(-1, R);
        int32_t uc = io.integer(-1, C);
        int32_t flat = io.integer(0, R * C - 1);
        auto M = Mat<R, C>(io, mode);
        auto M2 = Mat<R, C>(io, mode);
        auto row = Vec<C>(io, mode);
        auto col = Vec<R>(io, mode);
        Matrix<R, C, double> A = M;
        A.SetRow(r, row);
        A.SetCol(c, col);
        io.outMat(A);
        io.outVec(M.GetRow(r));
        io.outVec(M.GetCol(c));
        io.outReal(M(r, c));
        io.outReal(M[flat]);
        Matrix<R, C, double> zero, unit, identity;
        zero.MakeZero();
        unit.MakeUnit(ur, uc);
        identity.MakeIdentity();
        io.outMat(zero);
        io.outMat(unit);
        io.outMat(identity);
        io.outMat(Matrix<R, C, double>::Zero());
        io.outMat(Matrix<R, C, double>::Unit(ur, uc));
        io.outMat(Matrix<R, C, double>::Identity());
        io.outBool(M == M2);
        io.outBool(M != M2);
        io.outBool(M < M2);
        io.outBool(M <= M2);
        io.outBool(M > M2);
        io.outBool(M >= M2);
        io.outBool(M == M);
    }
}

ORACLE_CASE("Matrix.arithmetic")
{
    // Unary operator- and operator+, operator+/-, operator* and operator/ by a
    // scalar (the zero divisor yields the zero matrix) and the three norms.
    int32_t shape = io.integer(0, 2);
    int32_t mode = io.index() % 3;
    if (shape == 0)
    {
        MatArithmetic<2, 3>(io, mode);
    }
    else if (shape == 1)
    {
        MatArithmetic<3, 3>(io, mode);
    }
    else
    {
        MatArithmetic<4, 2>(io, mode);
    }
}

ORACLE_CASE("Matrix.multiply")
{
    int32_t shape = io.integer(0, 2);
    int32_t mode = io.index() % 3;
    if (shape == 0)
    {
        MatMultiply<2, 3, 4>(io, mode);
    }
    else if (shape == 1)
    {
        MatMultiply<3, 3, 3>(io, mode);
    }
    else
    {
        MatMultiply<4, 2, 3>(io, mode);
    }
}

ORACLE_CASE("Matrix.multiplyDiagonal")
{
    int32_t shape = io.integer(0, 2);
    int32_t mode = io.index() % 3;
    if (shape == 0)
    {
        MatDiagonal<2, 3>(io, mode);
    }
    else if (shape == 1)
    {
        MatDiagonal<3, 3>(io, mode);
    }
    else
    {
        MatDiagonal<4, 2>(io, mode);
    }
}

ORACLE_CASE("Matrix.inverse")
{
    // Only Matrix.h is included, so Inverse and Determinant resolve to the
    // Gaussian-elimination templates for every N (the closed forms of
    // Matrix2x2.h / Matrix3x3.h / Matrix4x4.h belong to group 2). Mode 2 is
    // exactly singular, mode 3 nearly singular and mode 4 wild, which reaches
    // the denormal-pivot NaN of upstream issue #375 that the port preserves.
    int32_t n = io.integer(2, 5);
    int32_t mode = io.index() % 5;
    if (n == 2)
    {
        MatInverse<2>(io, mode);
    }
    else if (n == 3)
    {
        MatInverse<3>(io, mode);
    }
    else if (n == 4)
    {
        MatInverse<4>(io, mode);
    }
    else
    {
        MatInverse<5>(io, mode);
    }
}

ORACLE_CASE("Matrix.inverse.denormal")
{
    int32_t n = io.integer(2, 4);
    if (n == 2)
    {
        DenormalInverse<2>(io);
    }
    else if (n == 3)
    {
        DenormalInverse<3>(io);
    }
    else
    {
        DenormalInverse<4>(io);
    }
}

ORACLE_CASE("Matrix.access")
{
    int32_t shape = io.integer(0, 2);
    int32_t mode = io.index() % 3;
    if (shape == 0)
    {
        MatAccess<2, 3>(io, mode);
    }
    else if (shape == 1)
    {
        MatAccess<3, 3>(io, mode);
    }
    else
    {
        MatAccess<4, 2>(io, mode);
    }
}

ORACLE_CASE("Matrix.hliftProject")
{
    int32_t n = io.integer(2, 4);
    int32_t mode = io.index() % 3;
    if (n == 2)
    {
        auto M = Mat<2, 2>(io, mode);
        io.outMat(HLift(M));
        io.outMat(HProject(M));
    }
    else if (n == 3)
    {
        auto M = Mat<3, 3>(io, mode);
        io.outMat(HLift(M));
        io.outMat(HProject(M));
    }
    else
    {
        auto M = Mat<4, 4>(io, mode);
        io.outMat(HLift(M));
        io.outMat(HProject(M));
    }
}

// ---------------------------------------------------------------------------
// GMatrix.h
//
// As for GVector, the port merges GMatrix into Matrix. The documented merge
// deviations kept out of these generators are division by zero (upstream
// throws) and LInfinityNorm's seed, which differs from Matrix.h's only when
// element 0 is NaN. Both are listed in the group report.
// ---------------------------------------------------------------------------

namespace
{
    // Runtime-sized counterpart of SquareMat, with the same five modes.
    GMatrix<double> SquareGMat(oracle::Ctx& io, int32_t n, int32_t mode)
    {
        GMatrix<double> m(n, n);
        if (mode != 2 && mode != 3)
        {
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = 0; c < n; ++c)
                {
                    m(r, c) = Scalar(io, mode == 4 ? 2 : mode);
                }
            }
            return m;
        }

        for (int32_t r = 0; r + 1 < n; ++r)
        {
            for (int32_t c = 0; c < n; ++c)
            {
                m(r, c) = static_cast<double>(io.rawInteger(-4, 4));
            }
        }
        for (int32_t c = 0; c < n; ++c)
        {
            m(n - 1, c) = 0.0;
        }
        for (int32_t r = 0; r + 1 < n; ++r)
        {
            double k = static_cast<double>(io.rawInteger(-2, 2));
            for (int32_t c = 0; c < n; ++c)
            {
                m(n - 1, c) += k * m(r, c);
            }
        }
        if (mode == 3)
        {
            int32_t c = io.rawInteger(0, n - 1);
            int32_t e = io.rawInteger(-1070, -40);
            m(n - 1, c) += std::ldexp(1.0, e);
        }
        for (int32_t r = 0; r < n; ++r)
        {
            for (int32_t c = 0; c < n; ++c)
            {
                io.given(m(r, c));
            }
        }
        return m;
    }
}

ORACLE_CASE("GMatrix.arithmetic")
{
    int32_t numRows = io.integer(1, 4);
    int32_t numCols = io.integer(1, 4);
    int32_t mode = io.index() % 3;
    auto m0 = GMat(io, numRows, numCols, mode);
    auto m1 = GMat(io, numRows, numCols, mode);
    double s = NonzeroScalar(io, mode);
    io.outGMat(-m0);
    io.outGMat(+m0);
    io.outGMat(m0 + m1);
    io.outGMat(m0 - m1);
    io.outGMat(m0 * s);
    io.outGMat(s * m0);
    io.outGMat(m0 / s);
    io.outReal(L1Norm(m0));
    io.outReal(L2Norm(m0));
    io.outReal(LInfinityNorm(m0));
}

ORACLE_CASE("GMatrix.multiply")
{
    int32_t r = io.integer(1, 4);
    int32_t k = io.integer(1, 4);
    int32_t c = io.integer(1, 4);
    int32_t mode = io.index() % 3;
    auto A = GMat(io, r, k, mode);
    auto B = GMat(io, k, c, mode);
    auto Bt = GMat(io, c, k, mode);
    auto At = GMat(io, k, r, mode);
    auto v = GVec(io, k, mode);
    auto w = GVec(io, r, mode);
    io.outGMat(A * B);
    io.outGMat(MultiplyAB(A, B));
    io.outGMat(MultiplyABT(A, Bt));
    io.outGMat(MultiplyATB(At, B));
    io.outGMat(MultiplyATBT(At, Bt));
    io.outGVec(A * v);
    io.outGVec(w * A);
    io.outGMat(Transpose(A));
}

ORACLE_CASE("GMatrix.multiplyDiagonal")
{
    int32_t numRows = io.integer(1, 4);
    int32_t numCols = io.integer(1, 4);
    int32_t mode = io.index() % 3;
    auto M = GMat(io, numRows, numCols, mode);
    auto dc = GVec(io, numCols, mode);
    auto dr = GVec(io, numRows, mode);
    io.outGMat(MultiplyMD(M, dc));
    io.outGMat(MultiplyDM(dr, M));
    io.outGMat(OuterProduct(dr, dc));
    GMatrix<double> diag(numRows, numCols);
    MakeDiagonal(dc, diag);
    io.outGMat(diag);
}

ORACLE_CASE("GMatrix.inverse")
{
    // GMatrix.h's MakeDiagonal is the nonsquare-tolerant version the port
    // shares; Inverse and Determinant are the same Gaussian elimination as
    // Matrix.h's.
    int32_t n = io.integer(1, 5);
    int32_t mode = io.index() % 5;
    auto M = SquareGMat(io, n, mode);
    bool invertible = false;
    GMatrix<double> invM = Inverse(M, &invertible);
    io.outBool(invertible);
    io.outGMat(invM);
    io.outReal(Determinant(M));
    io.outGMat(Inverse(M));
}

ORACLE_CASE("GMatrix.access")
{
    int32_t numRows = io.integer(1, 4);
    int32_t numCols = io.integer(1, 4);
    int32_t r = io.integer(0, numRows - 1);
    int32_t c = io.integer(0, numCols - 1);
    int32_t flat = io.integer(0, numRows * numCols - 1);
    int32_t newRows = io.integer(0, 5);
    int32_t newCols = io.integer(0, 5);
    int32_t mode = io.index() % 3;
    auto M = GMat(io, numRows, numCols, mode);
    auto M2 = GMat(io, numRows, numCols, mode);
    auto row = GVec(io, numCols, mode);
    auto col = GVec(io, numRows, mode);
    GMatrix<double> A = M;
    A.SetRow(r, row);
    A.SetCol(c, col);
    io.outGMat(A);
    io.outGVec(M.GetRow(r));
    io.outGVec(M.GetCol(c));
    io.outReal(M(r, c));
    io.outReal(M[flat]);
    io.outInt(M.GetNumRows());
    io.outInt(M.GetNumCols());
    io.outInt(M.GetNumElements());
    GMatrix<double> zero(numRows, numCols), unit(numRows, numCols),
        identity(numRows, numCols);
    zero.MakeZero();
    unit.MakeUnit(r, c);
    identity.MakeIdentity();
    io.outGMat(zero);
    io.outGMat(unit);
    io.outGMat(identity);
    io.outGMat(GMatrix<double>::Zero(numRows, numCols));
    io.outGMat(GMatrix<double>::Unit(numRows, numCols, r, c));
    io.outGMat(GMatrix<double>::Identity(numRows, numCols));
    io.outBool(M == M2);
    io.outBool(M != M2);
    io.outBool(M < M2);
    io.outBool(M <= M2);
    io.outBool(M > M2);
    io.outBool(M >= M2);
    io.outBool(M == M);
    GMatrix<double> resized = M;
    resized.SetSize(newRows, newCols);
    io.outGMat(resized);
}

ORACLE_CASE("GMatrix.access.mismatchedSizes")
{
    // Upstream GMatrix's comparisons are false for every relation but != when
    // the dimensions differ (issue #88, preserved by the port), and its
    // element access throws for an out-of-range index. Both sides must agree,
    // including which records throw.
    int32_t rows0 = io.integer(1, 3);
    int32_t cols0 = io.integer(1, 3);
    int32_t rows1 = io.integer(1, 3);
    int32_t cols1 = io.integer(1, 3);
    int32_t mode = io.index() % 3;
    auto M0 = GMat(io, rows0, cols0, mode);
    auto M1 = GMat(io, rows1, cols1, mode);
    io.outBool(M0 == M1);
    io.outBool(M0 != M1);
    io.outBool(M0 < M1);
    io.outBool(M0 <= M1);
    io.outBool(M0 > M1);
    io.outBool(M0 >= M1);
}

ORACLE_CASE("GMatrix.access.invalidIndex")
{
    // Throw parity: GMatrix's range-checked operator() calls LogError for an
    // out-of-range index, where the fixed-size Matrix accessor is unchecked.
    // About half of these records are throw records.
    int32_t numRows = io.integer(1, 3);
    int32_t numCols = io.integer(1, 3);
    int32_t r = io.integer(-1, 3);
    int32_t c = io.integer(-1, 3);
    int32_t mode = io.index() % 3;
    auto M = GMat(io, numRows, numCols, mode);
    io.outReal(M(r, c));
}
