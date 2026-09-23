// Verify group 42 (primitives): differential cases for the geometric value
// types AlignedBox.h, CanonicalBox.h, Halfspace.h, Hypersphere.h, Line.h,
// OrientedBox.h, Ray.h, Rectangle.h, Segment.h, SegmentMesh.h, Triangle.h,
// Arc2.h, Capsule.h, Circle3.h, Cylinder3.h, Ellipse3.h, Frustum3.h,
// Lozenge3.h, Parallelepiped3.h, Parallelogram2.h, Polyhedron3.h, Sector2.h
// and Torus3.h.
//
// Almost every compared path here is + - * / sqrt fabs and comparisons, so
// almost every case is declared exact (bit-identical) on the TypeScript side.
// The exceptions call the C math library and are named in their case comment:
//   Sector2.setAngle        std::cos, std::sin
//   Torus3.evaluate         std::cos, std::sin
//   Torus3.getParameters    std::atan2
// Sector2::Contains reads the stored cosAngle rather than recomputing it, so
// the containment case sets the member directly from a recorded input and is
// exact.
//
// Generator modes (selected by io.index() so that every record layout is
// fixed and the replay never has to know the mode):
//   0  small integer lattice: exact arithmetic, ties, exact zeros, coincident
//      and degenerate primitives
//   1  uniform
//   2  wild magnitudes: signed zeros, subnormals, 2^+-500, near overflow
//   3  wild magnitudes with NaN components (comparison cases only, because a
//      record whose real outputs are all NaN tests nothing: the harness
//      treats any NaN as equal to any NaN)
//
// One io draw per C++ statement: MSVC evaluates function arguments right to
// left, so every generated value goes into a named local before it is used.
#define ORACLE_FAMILY "v42-primitives"
#include "Oracle.h"

#include <Mathematics/AlignedBox.h>
#include <Mathematics/Arc2.h>
#include <Mathematics/CanonicalBox.h>
#include <Mathematics/Capsule.h>
#include <Mathematics/Circle3.h>
#include <Mathematics/Cylinder3.h>
#include <Mathematics/Ellipse3.h>
#include <Mathematics/Frustum3.h>
#include <Mathematics/Halfspace.h>
#include <Mathematics/Hypersphere.h>
#include <Mathematics/Line.h>
#include <Mathematics/Lozenge3.h>
#include <Mathematics/OrientedBox.h>
#include <Mathematics/Parallelepiped3.h>
#include <Mathematics/Parallelogram2.h>
#include <Mathematics/Polyhedron3.h>
#include <Mathematics/Ray.h>
#include <Mathematics/Rectangle.h>
#include <Mathematics/Sector2.h>
#include <Mathematics/Segment.h>
#include <Mathematics/SegmentMesh.h>
#include <Mathematics/Torus3.h>
#include <Mathematics/Triangle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

using namespace gte;

namespace
{
    // ---- unrecorded raw draws -------------------------------------------

    // A value drawn from the "wild" population. Only unrecorded raw draws are
    // used here; the caller records the resulting double.
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
        case 3:  return sign * std::ldexp(mantissa, -1022);  // smallest normal
        case 4:  return sign * std::ldexp(mantissa, 1020);   // near overflow
        case 5:  return sign * std::ldexp(mantissa, 500);
        case 6:  return sign * std::ldexp(mantissa, -500);
        case 7:  return sign * static_cast<double>(io.rawInteger(0, 4));
        case 8:  return sign * std::ldexp(mantissa, 52);     // integer-valued
        default: return sign * mantissa;
        }
    }

    // One unrecorded scalar in the requested mode. Mode 3 adds NaN.
    double RawScalar(oracle::Ctx& io, int32_t mode)
    {
        if (mode == 0)
        {
            return static_cast<double>(io.rawInteger(-3, 3));
        }
        if (mode == 1)
        {
            return io.raw(-8.0, 8.0);
        }
        if (mode == 2)
        {
            return WildRaw(io);
        }
        if (io.rawInteger(0, 3) == 0)
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return WildRaw(io);
    }

    // ---- recorded generators --------------------------------------------

    // One recorded scalar. Records one double.
    double Scalar(oracle::Ctx& io, int32_t mode)
    {
        return io.given(RawScalar(io, mode));
    }

    // A recorded vector. Records N doubles.
    template <int32_t N>
    Vector<N, double> Vec(oracle::Ctx& io, int32_t mode)
    {
        Vector<N, double> v{};
        for (int32_t i = 0; i < N; ++i)
        {
            v[i] = RawScalar(io, mode);
        }
        return io.givenVec(v);
    }

    // A recorded nonnegative extent vector; mode 0 makes zero extents common
    // (a box degenerate in that direction). Records N doubles.
    template <int32_t N>
    Vector<N, double> Extent(oracle::Ctx& io, int32_t mode)
    {
        Vector<N, double> e{};
        for (int32_t i = 0; i < N; ++i)
        {
            e[i] = std::fabs(RawScalar(io, mode));
        }
        return io.givenVec(e);
    }

    // A right-handed orthonormal frame of R^3 from unrecorded draws; the
    // caller records the axes it actually uses. GTE's convention is
    // axis[2] = Cross(axis[0], axis[1]), so this frame has determinant +1.
    std::array<Vector3<double>, 3> RawFrame3(oracle::Ctx& io)
    {
        std::array<Vector3<double>, 3> f{};
        f[0] = { 1.0, 0.0, 0.0 };
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            Vector3<double> candidate{};
            for (int32_t i = 0; i < 3; ++i) { candidate[i] = io.raw(-1.0, 1.0); }
            double len = Length(candidate);
            if (0.25 <= len && len <= 1.0) { f[0] = candidate; break; }
        }
        Normalize(f[0]);
        ComputeOrthogonalComplement<double>(1, f.data());
        return f;
    }

    // A right-handed orthonormal frame of R^2 from unrecorded draws. GTE's
    // Perp(x,y) = (y,-x) is the clockwise rotation, so the right-handed
    // partner of U = (x,y) is (-y,x) = -Perp(U), for which
    // DotPerp(axis[0], axis[1]) = 1 > 0.
    std::array<Vector2<double>, 2> RawFrame2(oracle::Ctx& io)
    {
        std::array<Vector2<double>, 2> f{};
        f[0] = { 1.0, 0.0 };
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            Vector2<double> candidate{};
            for (int32_t i = 0; i < 2; ++i) { candidate[i] = io.raw(-1.0, 1.0); }
            double len = Length(candidate);
            if (0.25 <= len && len <= 1.0) { f[0] = candidate; break; }
        }
        Normalize(f[0]);
        f[1] = { -f[0][1], f[0][0] };
        return f;
    }

    // An axis-aligned box. Some records are flat in every direction
    // (min == max) and some are inverted (min > max): upstream documents
    // min[i] <= max[i] but neither the constructor nor GetCenteredForm nor
    // GetVertices enforces it. Records 2*N doubles, min then max.
    template <int32_t N>
    AlignedBox<N, double> Box(oracle::Ctx& io, int32_t mode)
    {
        AlignedBox<N, double> box{};
        int32_t shape = io.rawInteger(0, 5);
        for (int32_t i = 0; i < N; ++i)
        {
            double x = RawScalar(io, mode);
            double y = RawScalar(io, mode);
            if (shape == 0) { y = x; }
            else if (shape == 1) { if (x < y) { std::swap(x, y); } }
            else if (x > y) { std::swap(x, y); }
            box.min[i] = x;
            box.max[i] = y;
        }
        io.givenVec(box.min);
        io.givenVec(box.max);
        return box;
    }

    // Two endpoints that coincide on some records (a zero-length segment,
    // whose GetCenteredForm normalizes the zero vector). Records 2*N doubles.
    template <int32_t N>
    Segment<N, double> Seg(oracle::Ctx& io, int32_t mode)
    {
        Segment<N, double> s{};
        bool coincident = (io.rawInteger(0, 4) == 0);
        for (int32_t i = 0; i < N; ++i) { s.p[0][i] = RawScalar(io, mode); }
        for (int32_t i = 0; i < N; ++i)
        {
            s.p[1][i] = (coincident ? s.p[0][i] : RawScalar(io, mode));
        }
        io.givenVec(s.p[0]);
        io.givenVec(s.p[1]);
        return s;
    }

    // ---- comparison operators -------------------------------------------

    // Draws the K member scalars of two objects of the same type, recording
    // the K of the first object and then the K of the second. A random prefix
    // of the second object's members is an exact copy of the first's (or, for
    // a zero, its negation), which is what makes the lexicographic member-by-
    // member ordering of upstream's operator< observable: without ties the
    // first member alone decides every comparison. Records 2*K doubles in
    // every mode.
    void ComparePair(oracle::Ctx& io, int32_t K,
        std::vector<double>& a, std::vector<double>& b)
    {
        int32_t mode = io.index() % 4;
        int32_t tie = io.rawInteger(0, K);
        a.resize(static_cast<size_t>(K));
        b.resize(static_cast<size_t>(K));
        for (int32_t i = 0; i < K; ++i)
        {
            a[static_cast<size_t>(i)] = RawScalar(io, mode);
        }
        for (int32_t i = 0; i < K; ++i)
        {
            double x = a[static_cast<size_t>(i)];
            if (i < tie)
            {
                // An exact tie, or the other signed zero: +0 == -0 is true
                // while neither is less than the other, so a -0/+0 tie must
                // not change any of the six relations.
                bool flip = (x == 0.0 && io.rawInteger(0, 2) == 0);
                b[static_cast<size_t>(i)] = (flip ? -x : x);
            }
            else
            {
                b[static_cast<size_t>(i)] = RawScalar(io, mode);
            }
        }
        for (int32_t i = 0; i < K; ++i) { io.given(a[static_cast<size_t>(i)]); }
        for (int32_t i = 0; i < K; ++i) { io.given(b[static_cast<size_t>(i)]); }
    }

    // The six comparison operators of A against B, the two reversed relations
    // (upstream defines operator> as the reversed operator< and operator<= as
    // its negation, which differ for NaN) and the four self-relations.
    template <typename T>
    void EmitComparisons(oracle::Ctx& io, T const& A, T const& B)
    {
        io.outBool(A == B);
        io.outBool(A != B);
        io.outBool(A < B);
        io.outBool(A <= B);
        io.outBool(A > B);
        io.outBool(A >= B);
        io.outBool(B < A);
        io.outBool(B <= A);
        io.outBool(A == A);
        io.outBool(A < A);
        io.outBool(A <= A);
        io.outBool(A >= A);
    }

    // Fills v[0..N-1] from src[offset..offset+N-1].
    template <int32_t N>
    Vector<N, double> Unpack(std::vector<double> const& src, int32_t offset)
    {
        Vector<N, double> v{};
        for (int32_t i = 0; i < N; ++i)
        {
            v[i] = src[static_cast<size_t>(offset + i)];
        }
        return v;
    }
}

// ---------------------------------------------------------------------------
// Comparison operators. Every type of this group that defines them gets a case
// with NaN components, -0/+0 ties and lexicographic prefix ties. The objects
// are built by assigning the public members rather than through the value
// constructors, because Parallelogram2 and Parallelepiped3 assert a
// right-handed basis, Sector2's constructor calls std::cos/std::sin and
// Frustum3's calls Update(); none of those touch the compared members.
// ---------------------------------------------------------------------------

ORACLE_CASE("AlignedBox.comparisons.3d")
{
    std::vector<double> a, b;
    ComparePair(io, 6, a, b);
    AlignedBox3<double> A{}, B{};
    A.min = Unpack<3>(a, 0);
    A.max = Unpack<3>(a, 3);
    B.min = Unpack<3>(b, 0);
    B.max = Unpack<3>(b, 3);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("CanonicalBox.comparisons.2d")
{
    std::vector<double> a, b;
    ComparePair(io, 2, a, b);
    CanonicalBox2<double> A{}, B{};
    A.extent = Unpack<2>(a, 0);
    B.extent = Unpack<2>(b, 0);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Halfspace.comparisons.3d")
{
    std::vector<double> a, b;
    ComparePair(io, 4, a, b);
    Halfspace3<double> A{}, B{};
    A.normal = Unpack<3>(a, 0);
    A.constant = a[3];
    B.normal = Unpack<3>(b, 0);
    B.constant = b[3];
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Hypersphere.comparisons.3d")
{
    std::vector<double> a, b;
    ComparePair(io, 4, a, b);
    Sphere3<double> A{}, B{};
    A.center = Unpack<3>(a, 0);
    A.radius = a[3];
    B.center = Unpack<3>(b, 0);
    B.radius = b[3];
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Line.comparisons.3d")
{
    std::vector<double> a, b;
    ComparePair(io, 6, a, b);
    Line3<double> A{}, B{};
    A.origin = Unpack<3>(a, 0);
    A.direction = Unpack<3>(a, 3);
    B.origin = Unpack<3>(b, 0);
    B.direction = Unpack<3>(b, 3);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Ray.comparisons.2d")
{
    std::vector<double> a, b;
    ComparePair(io, 4, a, b);
    Ray2<double> A{}, B{};
    A.origin = Unpack<2>(a, 0);
    A.direction = Unpack<2>(a, 2);
    B.origin = Unpack<2>(b, 0);
    B.direction = Unpack<2>(b, 2);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("OrientedBox.comparisons.3d")
{
    std::vector<double> a, b;
    ComparePair(io, 15, a, b);
    OrientedBox3<double> A{}, B{};
    A.center = Unpack<3>(a, 0);
    for (int32_t d = 0; d < 3; ++d) { A.axis[d] = Unpack<3>(a, 3 + 3 * d); }
    A.extent = Unpack<3>(a, 12);
    B.center = Unpack<3>(b, 0);
    for (int32_t d = 0; d < 3; ++d) { B.axis[d] = Unpack<3>(b, 3 + 3 * d); }
    B.extent = Unpack<3>(b, 12);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Rectangle.comparisons.3d")
{
    std::vector<double> a, b;
    ComparePair(io, 11, a, b);
    Rectangle3<double> A{}, B{};
    A.center = Unpack<3>(a, 0);
    for (int32_t d = 0; d < 2; ++d) { A.axis[d] = Unpack<3>(a, 3 + 3 * d); }
    A.extent = Unpack<2>(a, 9);
    B.center = Unpack<3>(b, 0);
    for (int32_t d = 0; d < 2; ++d) { B.axis[d] = Unpack<3>(b, 3 + 3 * d); }
    B.extent = Unpack<2>(b, 9);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Segment.comparisons.3d")
{
    std::vector<double> a, b;
    ComparePair(io, 6, a, b);
    Segment3<double> A{}, B{};
    A.p[0] = Unpack<3>(a, 0);
    A.p[1] = Unpack<3>(a, 3);
    B.p[0] = Unpack<3>(b, 0);
    B.p[1] = Unpack<3>(b, 3);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Triangle.comparisons.2d")
{
    // B goes through the std::array constructor and A through the
    // three-vertex one, so both are covered.
    std::vector<double> a, b;
    ComparePair(io, 6, a, b);
    Triangle2<double> A(Unpack<2>(a, 0), Unpack<2>(a, 2), Unpack<2>(a, 4));
    std::array<Vector2<double>, 3> bv{};
    for (int32_t k = 0; k < 3; ++k) { bv[k] = Unpack<2>(b, 2 * k); }
    Triangle2<double> B(bv);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Arc2.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 7, a, b);
    Arc2<double> A{}, B{};
    A.center = Unpack<2>(a, 0);
    A.radius = a[2];
    A.end[0] = Unpack<2>(a, 3);
    A.end[1] = Unpack<2>(a, 5);
    B.center = Unpack<2>(b, 0);
    B.radius = b[2];
    B.end[0] = Unpack<2>(b, 3);
    B.end[1] = Unpack<2>(b, 5);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Capsule.comparisons.3d")
{
    std::vector<double> a, b;
    ComparePair(io, 7, a, b);
    Capsule3<double> A{}, B{};
    A.segment.p[0] = Unpack<3>(a, 0);
    A.segment.p[1] = Unpack<3>(a, 3);
    A.radius = a[6];
    B.segment.p[0] = Unpack<3>(b, 0);
    B.segment.p[1] = Unpack<3>(b, 3);
    B.radius = b[6];
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Circle3.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 7, a, b);
    Circle3<double> A{}, B{};
    A.center = Unpack<3>(a, 0);
    A.normal = Unpack<3>(a, 3);
    A.radius = a[6];
    B.center = Unpack<3>(b, 0);
    B.normal = Unpack<3>(b, 3);
    B.radius = b[6];
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Cylinder3.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 8, a, b);
    Cylinder3<double> A{}, B{};
    A.axis.origin = Unpack<3>(a, 0);
    A.axis.direction = Unpack<3>(a, 3);
    A.radius = a[6];
    A.height = a[7];
    B.axis.origin = Unpack<3>(b, 0);
    B.axis.direction = Unpack<3>(b, 3);
    B.radius = b[6];
    B.height = b[7];
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Ellipse3.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 14, a, b);
    Ellipse3<double> A{}, B{};
    A.center = Unpack<3>(a, 0);
    A.normal = Unpack<3>(a, 3);
    A.axis[0] = Unpack<3>(a, 6);
    A.axis[1] = Unpack<3>(a, 9);
    A.extent = Unpack<2>(a, 12);
    B.center = Unpack<3>(b, 0);
    B.normal = Unpack<3>(b, 3);
    B.axis[0] = Unpack<3>(b, 6);
    B.axis[1] = Unpack<3>(b, 9);
    B.extent = Unpack<2>(b, 12);
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Frustum3.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 16, a, b);
    Frustum3<double> A{}, B{};
    A.origin = Unpack<3>(a, 0);
    A.dVector = Unpack<3>(a, 3);
    A.uVector = Unpack<3>(a, 6);
    A.rVector = Unpack<3>(a, 9);
    A.dMin = a[12];
    A.dMax = a[13];
    A.uBound = a[14];
    A.rBound = a[15];
    B.origin = Unpack<3>(b, 0);
    B.dVector = Unpack<3>(b, 3);
    B.uVector = Unpack<3>(b, 6);
    B.rVector = Unpack<3>(b, 9);
    B.dMin = b[12];
    B.dMax = b[13];
    B.uBound = b[14];
    B.rBound = b[15];
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Lozenge3.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 12, a, b);
    Lozenge3<double> A{}, B{};
    A.rectangle.center = Unpack<3>(a, 0);
    A.rectangle.axis[0] = Unpack<3>(a, 3);
    A.rectangle.axis[1] = Unpack<3>(a, 6);
    A.rectangle.extent = Unpack<2>(a, 9);
    A.radius = a[11];
    B.rectangle.center = Unpack<3>(b, 0);
    B.rectangle.axis[0] = Unpack<3>(b, 3);
    B.rectangle.axis[1] = Unpack<3>(b, 6);
    B.rectangle.extent = Unpack<2>(b, 9);
    B.radius = b[11];
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Parallelepiped3.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 12, a, b);
    Parallelepiped3<double> A{}, B{};
    A.center = Unpack<3>(a, 0);
    for (int32_t d = 0; d < 3; ++d) { A.axis[d] = Unpack<3>(a, 3 + 3 * d); }
    B.center = Unpack<3>(b, 0);
    for (int32_t d = 0; d < 3; ++d) { B.axis[d] = Unpack<3>(b, 3 + 3 * d); }
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Parallelogram2.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 6, a, b);
    Parallelogram2<double> A{}, B{};
    A.center = Unpack<2>(a, 0);
    for (int32_t d = 0; d < 2; ++d) { A.axis[d] = Unpack<2>(a, 2 + 2 * d); }
    B.center = Unpack<2>(b, 0);
    for (int32_t d = 0; d < 2; ++d) { B.axis[d] = Unpack<2>(b, 2 + 2 * d); }
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Sector2.comparisons")
{
    // Upstream compares vertex, radius, direction and angle but not the
    // derived cosAngle/sinAngle, so those members are left at their defaults.
    std::vector<double> a, b;
    ComparePair(io, 6, a, b);
    Sector2<double> A{}, B{};
    A.vertex = Unpack<2>(a, 0);
    A.radius = a[2];
    A.direction = Unpack<2>(a, 3);
    A.angle = a[5];
    B.vertex = Unpack<2>(b, 0);
    B.radius = b[2];
    B.direction = Unpack<2>(b, 3);
    B.angle = b[5];
    EmitComparisons(io, A, B);
}

ORACLE_CASE("Torus3.comparisons")
{
    std::vector<double> a, b;
    ComparePair(io, 14, a, b);
    Torus3<double> A{}, B{};
    A.center = Unpack<3>(a, 0);
    A.direction0 = Unpack<3>(a, 3);
    A.direction1 = Unpack<3>(a, 6);
    A.normal = Unpack<3>(a, 9);
    A.radius0 = a[12];
    A.radius1 = a[13];
    B.center = Unpack<3>(b, 0);
    B.direction0 = Unpack<3>(b, 3);
    B.direction1 = Unpack<3>(b, 6);
    B.normal = Unpack<3>(b, 9);
    B.radius0 = b[12];
    B.radius1 = b[13];
    EmitComparisons(io, A, B);
}

// ---------------------------------------------------------------------------
// Boxes, rectangles, segments: the centered forms and the corner
// reconstructions. Every one of these is a sum of scaled axes whose
// association order is observable, so all are compared bit for bit.
// ---------------------------------------------------------------------------

ORACLE_CASE("AlignedBox.centeredFormAndVertices.3d")
{
    int32_t mode = io.index() % 3;
    auto box = Box<3>(io, mode);
    Vector3<double> center{}, extent{};
    box.GetCenteredForm(center, extent);
    io.outVec(center);
    io.outVec(extent);
    std::array<Vector3<double>, 8> vertex{};
    box.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("AlignedBox.centeredFormAndVertices.2d")
{
    int32_t mode = io.index() % 3;
    auto box = Box<2>(io, mode);
    Vector2<double> center{}, extent{};
    box.GetCenteredForm(center, extent);
    io.outVec(center);
    io.outVec(extent);
    std::array<Vector2<double>, 4> vertex{};
    box.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("CanonicalBox.getVertices.3d")
{
    int32_t mode = io.index() % 3;
    auto extent = Extent<3>(io, mode);
    CanonicalBox3<double> box(extent);
    std::array<Vector3<double>, 8> vertex{};
    box.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("CanonicalBox.getVertices.2d")
{
    int32_t mode = io.index() % 3;
    auto extent = Extent<2>(io, mode);
    CanonicalBox2<double> box(extent);
    std::array<Vector2<double>, 4> vertex{};
    box.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("OrientedBox.getVertices.3d")
{
    // GetVertices forms product[d] = extent[d]*axis[d] once and then
    // accumulates center +- product[0] +- product[1] +- product[2] left to
    // right, which is what the exact comparison pins down.
    int32_t mode = io.index() % 3;
    auto center = Vec<3>(io, mode);
    auto frame = RawFrame3(io);
    std::array<Vector3<double>, 3> axis{};
    axis[0] = io.givenVec(frame[0]);
    axis[1] = io.givenVec(frame[1]);
    axis[2] = io.givenVec(frame[2]);
    auto extent = Extent<3>(io, mode);
    OrientedBox3<double> box(center, axis, extent);
    std::array<Vector3<double>, 8> vertex{};
    box.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("OrientedBox.getVertices.2d")
{
    int32_t mode = io.index() % 3;
    auto center = Vec<2>(io, mode);
    auto frame = RawFrame2(io);
    std::array<Vector2<double>, 2> axis{};
    axis[0] = io.givenVec(frame[0]);
    axis[1] = io.givenVec(frame[1]);
    auto extent = Extent<2>(io, mode);
    OrientedBox2<double> box(center, axis, extent);
    std::array<Vector2<double>, 4> vertex{};
    box.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("Rectangle.getVertices.3d")
{
    int32_t mode = io.index() % 3;
    auto center = Vec<3>(io, mode);
    auto frame = RawFrame3(io);
    std::array<Vector3<double>, 2> axis{};
    axis[0] = io.givenVec(frame[0]);
    axis[1] = io.givenVec(frame[1]);
    auto extent = Extent<2>(io, mode);
    Rectangle3<double> rectangle(center, axis, extent);
    std::array<Vector3<double>, 4> vertex{};
    rectangle.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("Rectangle.getVertices.2d")
{
    int32_t mode = io.index() % 3;
    auto center = Vec<2>(io, mode);
    auto frame = RawFrame2(io);
    std::array<Vector2<double>, 2> axis{};
    axis[0] = io.givenVec(frame[0]);
    axis[1] = io.givenVec(frame[1]);
    auto extent = Extent<2>(io, mode);
    Rectangle2<double> rectangle(center, axis, extent);
    std::array<Vector2<double>, 4> vertex{};
    rectangle.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("Segment.getCenteredForm.3d")
{
    int32_t mode = io.index() % 3;
    auto segment = Seg<3>(io, mode);
    Vector3<double> center{}, direction{};
    double extent = 0.0;
    segment.GetCenteredForm(center, direction, extent);
    io.outVec(center);
    io.outVec(direction);
    io.outReal(extent);
    // Round trip: upstream's own note says the endpoints need not come back.
    Segment3<double> round(center, direction, extent);
    io.outVec(round.p[0]);
    io.outVec(round.p[1]);
}

ORACLE_CASE("Segment.getCenteredForm.2d")
{
    int32_t mode = io.index() % 3;
    auto segment = Seg<2>(io, mode);
    Vector2<double> center{}, direction{};
    double extent = 0.0;
    segment.GetCenteredForm(center, direction, extent);
    io.outVec(center);
    io.outVec(direction);
    io.outReal(extent);
    Segment2<double> round(center, direction, extent);
    io.outVec(round.p[0]);
    io.outVec(round.p[1]);
}

ORACLE_CASE("Segment.setCenteredForm.3d")
{
    // SetCenteredForm from independent inputs, including a non-unit and a
    // zero direction, and the array constructor.
    int32_t mode = io.index() % 3;
    auto center = Vec<3>(io, mode);
    auto direction = Vec<3>(io, mode);
    double extent = Scalar(io, mode);
    Segment3<double> segment{};
    segment.SetCenteredForm(center, direction, extent);
    io.outVec(segment.p[0]);
    io.outVec(segment.p[1]);
    std::array<Vector3<double>, 2> p{ segment.p[0], segment.p[1] };
    Segment3<double> copy(p);
    io.outVec(copy.p[0]);
    io.outVec(copy.p[1]);
}

// ---------------------------------------------------------------------------
// Parallelepiped3 and Parallelogram2: right-handed, not necessarily orthogonal
// bases; the corner reconstruction accumulates left to right.
// ---------------------------------------------------------------------------

namespace
{
    // Three axes with DotCross(a0,a1,a2) > 0, which the Parallelepiped3
    // constructor asserts. Mode 0 draws integer axes and rejects the
    // left-handed and degenerate triples; the other modes shear a right-handed
    // orthonormal frame, which keeps the determinant equal to s0*s1*s2 in
    // exact arithmetic but is rejected when round-off makes it nonpositive.
    // Every rejection loop is capped and redraws every quantity it tests.
    // Records 9 doubles.
    std::array<Vector3<double>, 3> PpdAxes(oracle::Ctx& io, int32_t mode)
    {
        auto f = RawFrame3(io);
        std::array<Vector3<double>, 3> a{ f[0], f[1], f[2] };
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            std::array<Vector3<double>, 3> c{};
            if (mode == 0)
            {
                for (int32_t d = 0; d < 3; ++d)
                {
                    for (int32_t i = 0; i < 3; ++i)
                    {
                        c[d][i] = static_cast<double>(io.rawInteger(-3, 3));
                    }
                }
            }
            else
            {
                double s0 = std::fabs(RawScalar(io, mode));
                double s1 = std::fabs(RawScalar(io, mode));
                double s2 = std::fabs(RawScalar(io, mode));
                double t0 = RawScalar(io, mode);
                double t1 = RawScalar(io, mode);
                double t2 = RawScalar(io, mode);
                c[0] = s0 * f[0];
                c[1] = t0 * f[0] + s1 * f[1];
                c[2] = t1 * f[0] + t2 * f[1] + s2 * f[2];
            }
            if (DotCross(c[0], c[1], c[2]) > 0.0) { a = c; break; }
        }
        io.givenVec(a[0]);
        io.givenVec(a[1]);
        io.givenVec(a[2]);
        return a;
    }

    // Two axes with DotPerp(a0,a1) > 0, which the Parallelogram2 constructor
    // asserts. Records 4 doubles.
    std::array<Vector2<double>, 2> PgmAxes(oracle::Ctx& io, int32_t mode)
    {
        auto f = RawFrame2(io);
        std::array<Vector2<double>, 2> a{ f[0], f[1] };
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            std::array<Vector2<double>, 2> c{};
            if (mode == 0)
            {
                for (int32_t d = 0; d < 2; ++d)
                {
                    for (int32_t i = 0; i < 2; ++i)
                    {
                        c[d][i] = static_cast<double>(io.rawInteger(-3, 3));
                    }
                }
            }
            else
            {
                double s0 = std::fabs(RawScalar(io, mode));
                double s1 = std::fabs(RawScalar(io, mode));
                double t0 = RawScalar(io, mode);
                c[0] = s0 * f[0];
                c[1] = t0 * f[0] + s1 * f[1];
            }
            if (DotPerp(c[0], c[1]) > 0.0) { a = c; break; }
        }
        io.givenVec(a[0]);
        io.givenVec(a[1]);
        return a;
    }
}

ORACLE_CASE("Parallelepiped3.getVertices")
{
    int32_t mode = io.index() % 3;
    auto center = Vec<3>(io, mode);
    auto axis = PpdAxes(io, mode);
    Parallelepiped3<double> ppd(center, axis);
    std::array<Vector3<double>, 8> vertex{};
    ppd.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("Parallelepiped3.rightHandedAssert")
{
    // Throw parity for 'The axes must form a right-handed basis.'. The axes
    // are unconstrained, so many records are left-handed or degenerate and the
    // constructor throws.
    int32_t mode = io.index() % 3;
    auto center = Vec<3>(io, mode);
    auto a0 = Vec<3>(io, mode);
    auto a1 = Vec<3>(io, mode);
    auto a2 = Vec<3>(io, mode);
    std::array<Vector3<double>, 3> axis{ a0, a1, a2 };
    Parallelepiped3<double> ppd(center, axis);
    std::array<Vector3<double>, 8> vertex{};
    ppd.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("Parallelogram2.getVertices")
{
    int32_t mode = io.index() % 3;
    auto center = Vec<2>(io, mode);
    auto axis = PgmAxes(io, mode);
    Parallelogram2<double> pgm(center, axis);
    std::array<Vector2<double>, 4> vertex{};
    pgm.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

ORACLE_CASE("Parallelogram2.rightHandedAssert")
{
    int32_t mode = io.index() % 3;
    auto center = Vec<2>(io, mode);
    auto a0 = Vec<2>(io, mode);
    auto a1 = Vec<2>(io, mode);
    std::array<Vector2<double>, 2> axis{ a0, a1 };
    Parallelogram2<double> pgm(center, axis);
    std::array<Vector2<double>, 4> vertex{};
    pgm.GetVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
}

// ---------------------------------------------------------------------------
// Frustum3: the derived ratios of Update() and the eight corners.
// ---------------------------------------------------------------------------

namespace
{
    // An unrecorded strictly positive scalar (capped rejection).
    double RawPositive(oracle::Ctx& io, int32_t mode)
    {
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            double c = std::fabs(RawScalar(io, mode));
            if (c > 0.0) { return c; }
        }
        return 1.0;
    }
}

ORACLE_CASE("Frustum3.updateAndComputeVertices")
{
    // The default frustum has R = Unit(0), U = Unit(1) and
    // D = Unit(2) = Cross(R, U), so {rVector, uVector, dVector} is the
    // right-handed frame; the generator keeps that orientation.
    int32_t mode = io.index() % 3;
    auto origin = Vec<3>(io, mode);
    auto frame = RawFrame3(io);
    auto rVector = io.givenVec(frame[0]);
    auto uVector = io.givenVec(frame[1]);
    auto dVector = io.givenVec(frame[2]);
    double n0 = RawPositive(io, mode);
    double n1 = RawPositive(io, mode);
    double dMin = io.given(std::min(n0, n1));
    double dMax = io.given(std::max(n0, n1));
    double uBound = io.given(RawPositive(io, mode));
    double rBound = io.given(RawPositive(io, mode));
    Frustum3<double> frustum(origin, dVector, uVector, rVector,
        dMin, dMax, uBound, rBound);
    io.outReal(frustum.GetDRatio());
    io.outReal(frustum.GetMTwoUF());
    io.outReal(frustum.GetMTwoRF());
    std::array<Vector3<double>, 8> vertex{};
    frustum.ComputeVertices(vertex);
    for (auto const& v : vertex) { io.outVec(v); }
    // Update() again after changing the four scalars, which is the documented
    // way to use the class.
    double uBound2 = io.given(RawPositive(io, mode));
    frustum.uBound = uBound2;
    frustum.Update();
    io.outReal(frustum.GetDRatio());
    io.outReal(frustum.GetMTwoUF());
    io.outReal(frustum.GetMTwoRF());
}

// ---------------------------------------------------------------------------
// Arc2 and Sector2. Both are arithmetic only on the compared path: Arc2 uses
// Length, std::fabs and DotPerp, and Sector2::Contains reads the stored
// cosAngle, which this case records as an input rather than computing it with
// std::cos.
// ---------------------------------------------------------------------------

namespace
{
    // The twelve lattice points at distance exactly 5 from the origin, in
    // counterclockwise order. sqrt(25) is exact, so a point C + kCircle5[i]
    // with an integer center C satisfies |P - C| - 5 == 0 exactly and the
    // epsilon test of Arc2::Contains is evaluated at equality.
    Vector2<double> const kCircle5[12] =
    {
        { 5.0, 0.0 }, { 4.0, 3.0 }, { 3.0, 4.0 }, { 0.0, 5.0 },
        { -3.0, 4.0 }, { -4.0, 3.0 }, { -5.0, 0.0 }, { -4.0, -3.0 },
        { -3.0, -4.0 }, { 0.0, -5.0 }, { 3.0, -4.0 }, { 4.0, -3.0 }
    };
}

ORACLE_CASE("Arc2.contains")
{
    // Mode 0 places the center, both endpoints and the test point on the
    // radius-5 lattice circle, so |P-C| == radius exactly and
    // DotPerp(P-E0, E1-E0) is an exact integer that is zero whenever P is E0,
    // E1 or the third lattice point on the chord. epsilon is negative on a
    // sixth of the records: upstream's comment claims a negative epsilon
    // behaves as if zero was passed, but ||P-C| - r| >= 0 makes the test fail
    // for every P (UPSTREAM-FINDINGS, issue #155, preserved by the port).
    int32_t mode = io.index() % 3;
    Vector2<double> center{}, e0{}, e1{}, p{};
    double radius = 0.0;
    if (mode == 0)
    {
        for (int32_t i = 0; i < 2; ++i)
        {
            center[i] = static_cast<double>(io.rawInteger(-3, 3));
        }
        e0 = center + kCircle5[io.rawInteger(0, 11)];
        e1 = center + kCircle5[io.rawInteger(0, 11)];
        p = center + kCircle5[io.rawInteger(0, 11)];
        radius = 5.0;
    }
    else
    {
        for (int32_t i = 0; i < 2; ++i) { center[i] = RawScalar(io, mode); }
        for (int32_t i = 0; i < 2; ++i) { e0[i] = RawScalar(io, mode); }
        for (int32_t i = 0; i < 2; ++i) { e1[i] = RawScalar(io, mode); }
        for (int32_t i = 0; i < 2; ++i) { p[i] = RawScalar(io, mode); }
        // Half the records put P exactly on the circle by taking the radius
        // from the point, so the epsilon branch is reached in both directions.
        radius = (io.rawInteger(0, 1) == 0 ? Length(p - center)
            : std::fabs(RawScalar(io, mode)));
    }
    io.givenVec(center);
    io.given(radius);
    io.givenVec(e0);
    io.givenVec(e1);
    io.givenVec(p);
    int32_t epsKind = io.integer(0, 3);
    double epsilon = (epsKind == 0 ? -1.0 : (epsKind == 1 ? 0.0
        : (epsKind == 2 ? 1e-9 : 0.5)));
    Arc2<double> arc(center, radius, e0, e1);
    io.outBool(arc.Contains(p));
    io.outBool(arc.Contains(p, epsilon));
}

ORACLE_CASE("Sector2.contains")
{
    // Mode 0 builds exact boundary cases: the direction is an axis direction,
    // the test point is at distance exactly 5 from the vertex and cosAngle is
    // one of the dyadic-times-small-integer values for which 5*cosAngle is
    // exactly representable, so both 'length <= radius' and
    // 'Dot(direction,diff) >= length*cosAngle' are evaluated at equality.
    int32_t mode = io.index() % 3;
    Vector2<double> vertex{}, direction{};
    std::array<Vector2<double>, 3> p{};
    double radius = 0.0, cosAngle = 0.0, sinAngle = 0.0, angle = 0.0;
    if (mode == 0)
    {
        static double const kCos[9] =
            { -1.0, -0.8, -0.6, -0.5, 0.0, 0.5, 0.6, 0.8, 1.0 };
        for (int32_t i = 0; i < 2; ++i)
        {
            vertex[i] = static_cast<double>(io.rawInteger(-3, 3));
        }
        int32_t d = io.rawInteger(0, 3);
        direction = { (d == 0 ? 1.0 : (d == 2 ? -1.0 : 0.0)),
            (d == 1 ? 1.0 : (d == 3 ? -1.0 : 0.0)) };
        for (int32_t k = 0; k < 3; ++k)
        {
            p[k] = vertex + kCircle5[io.rawInteger(0, 11)];
        }
        radius = static_cast<double>(io.rawInteger(4, 6));
        cosAngle = kCos[io.rawInteger(0, 8)];
        sinAngle = 0.0;
        angle = 0.0;
    }
    else
    {
        auto frame = RawFrame2(io);
        for (int32_t i = 0; i < 2; ++i) { vertex[i] = RawScalar(io, mode); }
        direction = frame[0];
        for (int32_t k = 0; k < 3; ++k)
        {
            for (int32_t i = 0; i < 2; ++i) { p[k][i] = RawScalar(io, mode); }
        }
        radius = std::fabs(RawScalar(io, mode));
        cosAngle = io.raw(-1.0, 1.0);
        sinAngle = io.raw(-1.0, 1.0);
        angle = RawScalar(io, mode);
    }
    Sector2<double> sector{};
    sector.vertex = io.givenVec(vertex);
    sector.radius = io.given(radius);
    sector.direction = io.givenVec(direction);
    sector.angle = io.given(angle);
    sector.cosAngle = io.given(cosAngle);
    sector.sinAngle = io.given(sinAngle);
    for (int32_t k = 0; k < 3; ++k) { io.givenVec(p[k]); }
    for (int32_t k = 0; k < 3; ++k) { io.outBool(sector.Contains(p[k])); }
}

ORACLE_CASE("Sector2.setAngle")
{
    // std::cos and std::sin: the MSVC runtime and V8 may differ by 1 ulp, so
    // this is the only Sector2 case compared with a tolerance.
    double angle = io.real(-7.0, 7.0);
    Sector2<double> sector{};
    sector.SetAngle(angle);
    io.outReal(sector.angle);
    io.outReal(sector.cosAngle);
    io.outReal(sector.sinAngle);
}

// ---------------------------------------------------------------------------
// Torus3. Evaluate calls std::cos and std::sin and GetParameters calls
// std::atan2, so both are compared with the default scaled tolerance.
// ---------------------------------------------------------------------------

ORACLE_CASE("Torus3.evaluate")
{
    // maxOrder 3 exercises the preserved upstream quirk that the second-order
    // block is guarded by 'maxOrder == 2' rather than '>= 2': for maxOrder 3
    // only jet[0..2] are written (UPSTREAM-FINDINGS, issue #484).
    int32_t mode = io.index() % 2;
    auto center = Vec<3>(io, mode);
    auto frame = RawFrame3(io);
    auto direction0 = io.givenVec(frame[0]);
    auto direction1 = io.givenVec(frame[1]);
    auto normal = io.givenVec(frame[2]);
    double radius0 = io.given(RawPositive(io, mode));
    double radius1 = io.given(RawPositive(io, mode));
    double u = io.real(-7.0, 7.0);
    double v = io.real(-7.0, 7.0);
    int32_t maxOrder = io.integer(0, 3);
    Torus3<double> torus(center, direction0, direction1, normal,
        radius0, radius1);
    std::array<Vector3<double>, 6> jet{};
    torus.Evaluate(u, v, static_cast<uint32_t>(maxOrder), jet.data());
    int32_t count = (maxOrder == 0 ? 1 : (maxOrder == 2 ? 6 : 3));
    io.outInt(count);
    for (int32_t i = 0; i < count; ++i) { io.outVec(jet[i]); }
}

ORACLE_CASE("Torus3.getParameters")
{
    // Half the records take X off the torus and half take a point produced by
    // the parametric form; the point itself is recorded, so the replay never
    // repeats the cos/sin that made it.
    int32_t mode = io.index() % 2;
    auto center = Vec<3>(io, mode);
    auto frame = RawFrame3(io);
    auto direction0 = io.givenVec(frame[0]);
    auto direction1 = io.givenVec(frame[1]);
    auto normal = io.givenVec(frame[2]);
    double radius0 = io.given(RawPositive(io, mode));
    double radius1 = io.given(RawPositive(io, mode));
    Torus3<double> torus(center, direction0, direction1, normal,
        radius0, radius1);
    Vector3<double> X{};
    if (io.rawInteger(0, 1) == 0)
    {
        double rawU = io.raw(-7.0, 7.0);
        double rawV = io.raw(-7.0, 7.0);
        std::array<Vector3<double>, 6> jet{};
        torus.Evaluate(rawU, rawV, 0, jet.data());
        X = jet[0];
    }
    else
    {
        for (int32_t i = 0; i < 3; ++i) { X[i] = RawScalar(io, mode); }
    }
    io.givenVec(X);
    double u = 0.0, v = 0.0;
    torus.GetParameters(X, u, v);
    io.outReal(u);
    io.outReal(v);
}

// ---------------------------------------------------------------------------
// Cylinder3: the finite/infinite height sentinel helpers.
// ---------------------------------------------------------------------------

ORACLE_CASE("Cylinder3.finiteInfinite")
{
    int32_t mode = io.index() % 3;
    auto origin = Vec<3>(io, mode);
    auto frame = RawFrame3(io);
    auto direction = io.givenVec(frame[2]);
    double radius = Scalar(io, mode);
    double height = Scalar(io, mode);
    Line3<double> axis(origin, direction);
    Cylinder3<double> cylinder(axis, radius, height);
    io.outReal(cylinder.height);
    io.outBool(cylinder.IsFinite());
    io.outBool(cylinder.IsInfinite());
    double newHeight = Scalar(io, mode);
    cylinder.MakeFiniteCylinder(newHeight);
    io.outReal(cylinder.height);
    io.outBool(cylinder.IsFinite());
    cylinder.MakeInfiniteCylinder();
    io.outReal(cylinder.height);
    io.outBool(cylinder.IsFinite());
    io.outBool(cylinder.IsInfinite());
}

// ---------------------------------------------------------------------------
// Polyhedron3: the three geometric queries and the construction validation.
// ComputeVertexAverage ends with 'average /= n', which Vector.h implements as
// a multiplication by 1/n, so the accumulation order (ascending unique index,
// std::set order) and the reciprocal are both observable. ComputeVolume ends
// with a scalar 'volume /= 6', a true division.
// ---------------------------------------------------------------------------

namespace
{
    // The four faces of the tetrahedron <0,1,2,3>, a closed surface for any
    // four points, so the divergence-theorem volume is |det|/6 exactly.
    int32_t const kTetraFaces[12] =
        { 0, 1, 2,  0, 2, 3,  0, 3, 1,  1, 3, 2 };

    // The twelve faces of the box whose corner i has bit pattern b2b1b0.
    int32_t const kBoxFaces[36] =
    {
        0, 2, 3,  0, 3, 1,   // z = min
        4, 5, 7,  4, 7, 6,   // z = max
        0, 1, 5,  0, 5, 4,   // y = min
        2, 6, 7,  2, 7, 3,   // y = max
        0, 4, 6,  0, 6, 2,   // x = min
        1, 3, 7,  1, 7, 5    // x = max
    };
}

ORACLE_CASE("Polyhedron3.queries")
{
    // shape 0: a tetrahedron on four pool vertices (plus unused pool entries,
    // which the class supports); shape 1: the eight corners of an aligned box
    // with its twelve triangles; shape 2: an arbitrary index list over the
    // pool, which the queries also accept.
    int32_t mode = io.index() % 3;
    int32_t shape = io.integer(0, 2);
    int32_t extra = io.integer(0, 2);
    int32_t poolSize = (shape == 0 ? 4 : (shape == 1 ? 8 : 4)) + extra;
    auto pool = std::make_shared<std::vector<Vector3<double>>>(poolSize);
    if (shape == 1)
    {
        Vector3<double> lo{}, hi{};
        for (int32_t i = 0; i < 3; ++i)
        {
            double x = RawScalar(io, mode);
            double y = RawScalar(io, mode);
            lo[i] = std::min(x, y);
            hi[i] = std::max(x, y);
        }
        for (int32_t k = 0; k < 8; ++k)
        {
            for (int32_t d = 0; d < 3; ++d)
            {
                (*pool)[k][d] = ((k & (1 << d)) != 0 ? hi[d] : lo[d]);
            }
        }
        for (int32_t k = 8; k < poolSize; ++k)
        {
            for (int32_t d = 0; d < 3; ++d) { (*pool)[k][d] = RawScalar(io, mode); }
        }
    }
    else
    {
        for (int32_t k = 0; k < poolSize; ++k)
        {
            for (int32_t d = 0; d < 3; ++d) { (*pool)[k][d] = RawScalar(io, mode); }
        }
    }
    for (int32_t k = 0; k < poolSize; ++k) { io.givenVec((*pool)[k]); }

    std::vector<int32_t> indices;
    if (shape == 0)
    {
        indices.assign(kTetraFaces, kTetraFaces + 12);
    }
    else if (shape == 1)
    {
        indices.assign(kBoxFaces, kBoxFaces + 36);
    }
    else
    {
        int32_t numTriangles = io.rawInteger(4, 7);
        indices.resize(static_cast<size_t>(3 * numTriangles));
        for (auto& index : indices)
        {
            index = io.rawInteger(0, poolSize - 1);
        }
    }
    io.given(static_cast<double>(indices.size()));
    int32_t numIndices = static_cast<int32_t>(indices.size());
    for (auto index : indices) { io.given(static_cast<double>(index)); }
    bool counterClockwise = io.boolean();

    Polyhedron3<double> polyhedron(pool, numIndices, indices.data(),
        counterClockwise);
    io.outBool(static_cast<bool>(polyhedron));
    io.outBool(polyhedron.CounterClockwise());
    io.outInt(polyhedron.GetVertices().size());
    io.outInt(polyhedron.GetUniqueIndices().size());
    for (int32_t index : polyhedron.GetUniqueIndices()) { io.outInt(index); }
    io.outInt(polyhedron.GetIndices().size());
    io.outVec(polyhedron.ComputeVertexAverage());
    io.outReal(polyhedron.ComputeSurfaceArea());
    io.outReal(polyhedron.ComputeVolume());
}

ORACLE_CASE("Polyhedron3.invalidInput")
{
    // numIndices < 12 and numIndices % 3 != 0 both fail construction; the
    // object then reports invalid, has no indices and returns zero from the
    // three queries.
    int32_t mode = io.index() % 3;
    int32_t poolSize = io.integer(4, 6);
    std::vector<Vector3<double>> pool(static_cast<size_t>(poolSize));
    for (auto& vertex : pool)
    {
        for (int32_t d = 0; d < 3; ++d) { vertex[d] = RawScalar(io, mode); }
    }
    for (auto const& vertex : pool) { io.givenVec(vertex); }
    int32_t numIndices = io.integer(3, 15);
    std::vector<int32_t> indices(static_cast<size_t>(numIndices));
    for (auto& index : indices) { index = io.integer(0, poolSize - 1); }
    bool counterClockwise = io.boolean();
    auto shared = std::make_shared<std::vector<Vector3<double>>>(pool);
    Polyhedron3<double> polyhedron(shared, numIndices, indices.data(),
        counterClockwise);
    io.outBool(static_cast<bool>(polyhedron));
    io.outBool(polyhedron.CounterClockwise());
    io.outInt(polyhedron.GetUniqueIndices().size());
    io.outInt(polyhedron.GetIndices().size());
    io.outVec(polyhedron.ComputeVertexAverage());
    io.outReal(polyhedron.ComputeSurfaceArea());
    io.outReal(polyhedron.ComputeVolume());
}

// ---------------------------------------------------------------------------
// SegmentMesh: the four constructors, the index lists they derive and the
// three LogAssert conditions (throw parity).
// ---------------------------------------------------------------------------

namespace
{
    // Runs one SegmentMesh construction for dimension N. All inputs are drawn
    // before the constructor runs so that a throwing record still records the
    // full input layout.
    template <int32_t N>
    void SegmentMeshCase(oracle::Ctx& io)
    {
        int32_t mode = io.index() % 3;
        int32_t kind = io.integer(0, 4);
        int32_t numVertices = io.integer(1, 6);
        std::vector<Vector<N, double>> vertices(static_cast<size_t>(numVertices));
        for (auto& vertex : vertices)
        {
            for (int32_t d = 0; d < N; ++d) { vertex[d] = RawScalar(io, mode); }
        }
        for (auto const& vertex : vertices) { io.givenVec(vertex); }
        bool validate = io.boolean();
        int32_t numPairs = io.integer(0, 4);
        std::vector<std::array<size_t, 2>> indices(static_cast<size_t>(numPairs));
        for (auto& pair : indices)
        {
            // Out-of-range indices are drawn only when the caller asks for
            // validation. Without validation upstream stores the raw size_t
            // and never looks at it, so an out-of-range value there is
            // undefined use of the class rather than a disagreement to
            // measure.
            int32_t hi = (validate ? numVertices : numVertices - 1);
            int32_t i0 = io.integer(0, std::max(hi, 0));
            int32_t i1 = io.integer(0, std::max(hi, 0));
            pair[0] = static_cast<size_t>(i0);
            pair[1] = static_cast<size_t>(i1);
        }

        SegmentMesh<N, double> mesh{};
        if (kind == 1)
        {
            mesh = SegmentMesh<N, double>(vertices);
        }
        else if (kind == 2)
        {
            mesh = SegmentMesh<N, double>(vertices, true);
        }
        else if (kind == 3)
        {
            mesh = SegmentMesh<N, double>(vertices, false);
        }
        else if (kind == 4)
        {
            mesh = SegmentMesh<N, double>(vertices, indices, validate);
        }

        io.outInt(static_cast<int32_t>(mesh.GetTopology()));
        io.outInt(mesh.GetVertices().size());
        for (auto const& vertex : mesh.GetVertices()) { io.outVec(vertex); }
        io.outInt(mesh.GetIndices().size());
        for (auto const& pair : mesh.GetIndices())
        {
            io.outInt(pair[0]);
            io.outInt(pair[1]);
        }
    }
}

ORACLE_CASE("SegmentMesh.constructors.2d")
{
    SegmentMeshCase<2>(io);
}

ORACLE_CASE("SegmentMesh.constructors.3d")
{
    SegmentMeshCase<3>(io);
}

// ---------------------------------------------------------------------------
// The default constructors of every type of the group. These are the only
// paths that build the canonical axes and extents, and several of them differ
// from what the upstream class comment claims (Cylinder3's default axis is
// Line3's (1,0,0), not the documented (0,0,1)).
// ---------------------------------------------------------------------------

namespace
{
    template <int32_t N>
    void EmitTemplatedDefaults(oracle::Ctx& io)
    {
        AlignedBox<N, double> box{};
        io.outVec(box.min);
        io.outVec(box.max);
        CanonicalBox<N, double> canonicalBox{};
        io.outVec(canonicalBox.extent);
        Halfspace<N, double> halfspace{};
        io.outVec(halfspace.normal);
        io.outReal(halfspace.constant);
        Hypersphere<N, double> hypersphere{};
        io.outVec(hypersphere.center);
        io.outReal(hypersphere.radius);
        Line<N, double> line{};
        io.outVec(line.origin);
        io.outVec(line.direction);
        Ray<N, double> ray{};
        io.outVec(ray.origin);
        io.outVec(ray.direction);
        OrientedBox<N, double> orientedBox{};
        io.outVec(orientedBox.center);
        for (int32_t d = 0; d < N; ++d) { io.outVec(orientedBox.axis[d]); }
        io.outVec(orientedBox.extent);
        Rectangle<N, double> rectangle{};
        io.outVec(rectangle.center);
        io.outVec(rectangle.axis[0]);
        io.outVec(rectangle.axis[1]);
        io.outVec(rectangle.extent);
        Segment<N, double> segment{};
        io.outVec(segment.p[0]);
        io.outVec(segment.p[1]);
        Triangle<N, double> triangle{};
        for (int32_t k = 0; k < 3; ++k) { io.outVec(triangle.v[k]); }
        Capsule<N, double> capsule{};
        io.outVec(capsule.segment.p[0]);
        io.outVec(capsule.segment.p[1]);
        io.outReal(capsule.radius);
    }
}

ORACLE_CASE("Primitives.defaultConstructors")
{
    int32_t n = io.integer(2, 4);
    if (n == 2) { EmitTemplatedDefaults<2>(io); }
    else if (n == 3) { EmitTemplatedDefaults<3>(io); }
    else { EmitTemplatedDefaults<4>(io); }

    Arc2<double> arc{};
    io.outVec(arc.center);
    io.outReal(arc.radius);
    io.outVec(arc.end[0]);
    io.outVec(arc.end[1]);
    Circle3<double> circle{};
    io.outVec(circle.center);
    io.outVec(circle.normal);
    io.outReal(circle.radius);
    Cylinder3<double> cylinder{};
    io.outVec(cylinder.axis.origin);
    io.outVec(cylinder.axis.direction);
    io.outReal(cylinder.radius);
    io.outReal(cylinder.height);
    io.outBool(cylinder.IsFinite());
    Ellipse3<double> ellipse{};
    io.outVec(ellipse.center);
    io.outVec(ellipse.normal);
    io.outVec(ellipse.axis[0]);
    io.outVec(ellipse.axis[1]);
    io.outVec(ellipse.extent);
    Frustum3<double> frustum{};
    io.outVec(frustum.origin);
    io.outVec(frustum.dVector);
    io.outVec(frustum.uVector);
    io.outVec(frustum.rVector);
    io.outReal(frustum.dMin);
    io.outReal(frustum.dMax);
    io.outReal(frustum.uBound);
    io.outReal(frustum.rBound);
    io.outReal(frustum.GetDRatio());
    io.outReal(frustum.GetMTwoUF());
    io.outReal(frustum.GetMTwoRF());
    std::array<Vector3<double>, 8> frustumVertex{};
    frustum.ComputeVertices(frustumVertex);
    for (auto const& v : frustumVertex) { io.outVec(v); }
    Lozenge3<double> lozenge{};
    io.outVec(lozenge.rectangle.center);
    io.outVec(lozenge.rectangle.axis[0]);
    io.outVec(lozenge.rectangle.axis[1]);
    io.outVec(lozenge.rectangle.extent);
    io.outReal(lozenge.radius);
    Parallelepiped3<double> parallelepiped{};
    io.outVec(parallelepiped.center);
    for (int32_t d = 0; d < 3; ++d) { io.outVec(parallelepiped.axis[d]); }
    std::array<Vector3<double>, 8> ppdVertex{};
    parallelepiped.GetVertices(ppdVertex);
    for (auto const& v : ppdVertex) { io.outVec(v); }
    Parallelogram2<double> parallelogram{};
    io.outVec(parallelogram.center);
    for (int32_t d = 0; d < 2; ++d) { io.outVec(parallelogram.axis[d]); }
    std::array<Vector2<double>, 4> pgmVertex{};
    parallelogram.GetVertices(pgmVertex);
    for (auto const& v : pgmVertex) { io.outVec(v); }
    Sector2<double> sector{};
    io.outVec(sector.vertex);
    io.outReal(sector.radius);
    io.outVec(sector.direction);
    io.outReal(sector.angle);
    io.outReal(sector.cosAngle);
    io.outReal(sector.sinAngle);
    Torus3<double> torus{};
    io.outVec(torus.center);
    io.outVec(torus.direction0);
    io.outVec(torus.direction1);
    io.outVec(torus.normal);
    io.outReal(torus.radius0);
    io.outReal(torus.radius1);
    SegmentMesh<3, double> mesh{};
    io.outInt(static_cast<int32_t>(mesh.GetTopology()));
    io.outInt(mesh.GetVertices().size());
    io.outInt(mesh.GetIndices().size());
}
