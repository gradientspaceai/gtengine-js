// Group 43 (primitives): Hyperellipsoid, Hyperplane, RectangleManager,
// AlignedBoxBV, Cone, Polygon2, RectangleMesh, RectanglePatchMesh,
// Tetrahedron3, the four AlignedBox/OrientedBox trees of points, segments and
// triangles and the two bounding-volume classes.
//
// NOTE on overload resolution (v34): Hyperellipsoid::FromCoefficients calls
// Inverse(A, &invertible). Matrix2x2.h and Matrix3x3.h declare more
// specialized closed-form Inverse overloads than the Gaussian-elimination
// template of Matrix.h, and the call is dependent, so the overload chosen
// depends on what is visible in the translation unit at the point of
// instantiation. This file includes Matrix2x2.h and Matrix3x3.h explicitly
// (Cone.h and ContOrientedBox3.h pull Matrix3x3.h in anyway), so the closed
// form is used for N = 2 and N = 3, which is what the port assumes.
#define ORACLE_FAMILY "v43-primitives"
#include "Oracle.h"

#include <Mathematics/Matrix2x2.h>
#include <Mathematics/Matrix3x3.h>
#include <Mathematics/Hyperellipsoid.h>
#include <Mathematics/Hyperplane.h>
#include <Mathematics/Cone.h>
#include <Mathematics/Polygon2.h>
#include <Mathematics/Tetrahedron3.h>
#include <Mathematics/RectangleManager.h>
#include <Mathematics/RectangleMesh.h>
#include <Mathematics/RectanglePatchMesh.h>
#include <Mathematics/AlignedBoxBV.h>
#include <Mathematics/OrientedBoxBV.h>
#include <Mathematics/AlignedBoxTreeOfPoints.h>
#include <Mathematics/AlignedBoxTreeOfSegments.h>
#include <Mathematics/AlignedBoxTreeOfTriangles.h>
#include <Mathematics/OrientedBoxTreeOfPoints.h>
#include <Mathematics/OrientedBoxTreeOfSegments.h>
#include <Mathematics/OrientedBoxTreeOfTriangles.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <set>
#include <vector>

using namespace gte;

namespace
{
    // A right-handed orthonormal frame. Every mode records exactly 9 doubles.
    // Mode 0 is an exactly representable signed permutation of the standard
    // basis, which produces the exact-arithmetic configurations (axis-aligned
    // ellipsoids, boxes whose axes are the coordinate axes).
    std::array<Vector3<double>, 3> MakeFrame3(oracle::Ctx& io)
    {
        std::array<Vector3<double>, 3> b{};
        if (io.rawInteger(0, 3) == 0)
        {
            int32_t p = io.rawInteger(0, 2);
            double s = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
            b[0].MakeUnit(p);
            b[0] = s * b[0];
            b[1].MakeUnit((p + 1) % 3);
            b[2] = Cross(b[0], b[1]);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int32_t i = 0; i < 3; ++i) { b[0][i] = io.raw(-1.0, 1.0); }
                len = Length(b[0]);
            }
            while (len < 0.25 || len > 1.0);
            Normalize(b[0]);
            ComputeOrthogonalComplement(1, b.data());
        }
        io.givenVec(b[0]);
        io.givenVec(b[1]);
        io.givenVec(b[2]);
        return b;
    }

    // A right-handed orthonormal frame in 2D: {u, -Perp(u)}. Records 4
    // doubles.
    std::array<Vector2<double>, 2> MakeFrame2(oracle::Ctx& io)
    {
        std::array<Vector2<double>, 2> b{};
        if (io.rawInteger(0, 3) == 0)
        {
            int32_t p = io.rawInteger(0, 1);
            double s = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
            b[0].MakeUnit(p);
            b[0] = s * b[0];
            b[1] = -Perp(b[0]);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int32_t i = 0; i < 2; ++i) { b[0][i] = io.raw(-1.0, 1.0); }
                len = Length(b[0]);
            }
            while (len < 0.25 || len > 1.0);
            Normalize(b[0]);
            ComputeOrthogonalComplement(1, b.data());
        }
        io.givenVec(b[0]);
        io.givenVec(b[1]);
        return b;
    }

    // Positive extents; every 4th record uses small integers so the algebra
    // is exact.
    template <int32_t N>
    Vector<N, double> MakeExtent(oracle::Ctx& io)
    {
        Vector<N, double> e{};
        bool lattice = (io.rawInteger(0, 3) == 0);
        for (int32_t i = 0; i < N; ++i)
        {
            e[i] = (lattice ? io.lattice(1, 4) : io.real(0.25, 4.0));
        }
        return e;
    }
}

// ---------------------------------------------------------------- Hyperellipsoid

ORACLE_CASE("Hyperellipsoid.getM.2d")
{
    auto center = io.vec<2>(-4.0, 4.0);
    auto axis = MakeFrame2(io);
    auto extent = MakeExtent<2>(io);
    Ellipse2<double> E(center, axis, extent);
    Matrix<2, 2, double> M{}, MInverse{};
    E.GetM(M);
    E.GetMInverse(MInverse);
    io.outMat(M);
    io.outMat(MInverse);
}

ORACLE_CASE("Hyperellipsoid.getM.3d")
{
    auto center = io.vec<3>(-4.0, 4.0);
    auto axis = MakeFrame3(io);
    auto extent = MakeExtent<3>(io);
    Ellipsoid3<double> E(center, axis, extent);
    Matrix<3, 3, double> M{}, MInverse{};
    E.GetM(M);
    E.GetMInverse(MInverse);
    io.outMat(M);
    io.outMat(MInverse);
}

ORACLE_CASE("Hyperellipsoid.toCoefficients.2d")
{
    auto center = io.vec<2>(-4.0, 4.0);
    auto axis = MakeFrame2(io);
    auto extent = MakeExtent<2>(io);
    Ellipse2<double> E(center, axis, extent);
    Matrix<2, 2, double> A{};
    Vector<2, double> B{};
    double C{};
    E.ToCoefficients(A, B, C);
    std::array<double, 6> coeff{};
    E.ToCoefficients(coeff);
    io.outMat(A);
    io.outVec(B);
    io.outReal(C);
    for (double k : coeff) { io.outReal(k); }
}

ORACLE_CASE("Hyperellipsoid.toCoefficients.3d")
{
    auto center = io.vec<3>(-4.0, 4.0);
    auto axis = MakeFrame3(io);
    auto extent = MakeExtent<3>(io);
    Ellipsoid3<double> E(center, axis, extent);
    Matrix<3, 3, double> A{};
    Vector<3, double> B{};
    double C{};
    E.ToCoefficients(A, B, C);
    std::array<double, 10> coeff{};
    E.ToCoefficients(coeff);
    io.outMat(A);
    io.outVec(B);
    io.outReal(C);
    for (double k : coeff) { io.outReal(k); }
}

ORACLE_CASE("Hyperellipsoid.fromCoefficients.2d")
{
    // Round trip: E -> coefficients -> F. FromCoefficients runs the
    // closed-form Inverse of Matrix2x2.h and the N-dimensional
    // SymmetricEigensolver, both of which use only + - * / and sqrt.
    auto center = io.vec<2>(-4.0, 4.0);
    auto axis = MakeFrame2(io);
    auto extent = MakeExtent<2>(io);
    Ellipse2<double> E(center, axis, extent);
    std::array<double, 6> coeff{};
    E.ToCoefficients(coeff);
    Ellipse2<double> F{};
    bool valid = F.FromCoefficients(coeff);
    io.outBool(valid);
    if (valid)
    {
        io.outVec(F.center);
        io.outVec(F.axis[0]);
        io.outVec(F.axis[1]);
        io.outVec(F.extent);
    }
}

ORACLE_CASE("Hyperellipsoid.fromCoefficients.3d")
{
    auto center = io.vec<3>(-4.0, 4.0);
    auto axis = MakeFrame3(io);
    auto extent = MakeExtent<3>(io);
    Ellipsoid3<double> E(center, axis, extent);
    std::array<double, 10> coeff{};
    E.ToCoefficients(coeff);
    Ellipsoid3<double> F{};
    bool valid = F.FromCoefficients(coeff);
    io.outBool(valid);
    if (valid)
    {
        io.outVec(F.center);
        io.outVec(F.axis[0]);
        io.outVec(F.axis[1]);
        io.outVec(F.axis[2]);
        io.outVec(F.extent);
    }
}

ORACLE_CASE("Hyperellipsoid.fromCoefficientsABC.2d")
{
    // Small-integer (A, B, C). 'definite' makes A diagonally dominant and
    // positive definite so the accepting path is reached; the other half
    // reaches the rejection branches (a singular A, a zero right side, a
    // nonpositive eigenvalue).
    bool definite = io.boolean();
    double a00 = (definite ? io.lattice(2, 4) : io.lattice(-3, 3));
    double a01 = (definite ? io.lattice(-1, 1) : io.lattice(-3, 3));
    double a11 = (definite ? io.lattice(2, 4) : io.lattice(-3, 3));
    auto B = io.latticeVec<2>(-3, 3);
    double C = (definite ? io.lattice(-4, -1) : io.lattice(-3, 3));
    Matrix<2, 2, double> A{};
    A(0, 0) = a00;
    A(0, 1) = a01;
    A(1, 0) = a01;
    A(1, 1) = a11;
    Ellipse2<double> F{};
    bool valid = F.FromCoefficients(A, B, C);
    io.outBool(valid);
    if (valid)
    {
        io.outVec(F.center);
        io.outVec(F.axis[0]);
        io.outVec(F.axis[1]);
        io.outVec(F.extent);
    }
}

ORACLE_CASE("Hyperellipsoid.fromCoefficientsABC.3d")
{
    // See the 2d case: 'definite' reaches the accepting path.
    bool definite = io.boolean();
    double a00 = (definite ? io.lattice(3, 5) : io.lattice(-3, 3));
    double a01 = (definite ? io.lattice(-1, 1) : io.lattice(-3, 3));
    double a02 = (definite ? io.lattice(-1, 1) : io.lattice(-3, 3));
    double a11 = (definite ? io.lattice(3, 5) : io.lattice(-3, 3));
    double a12 = (definite ? io.lattice(-1, 1) : io.lattice(-3, 3));
    double a22 = (definite ? io.lattice(3, 5) : io.lattice(-3, 3));
    auto B = io.latticeVec<3>(-3, 3);
    double C = (definite ? io.lattice(-4, -1) : io.lattice(-3, 3));
    Matrix<3, 3, double> A{};
    A(0, 0) = a00; A(0, 1) = a01; A(0, 2) = a02;
    A(1, 0) = a01; A(1, 1) = a11; A(1, 2) = a12;
    A(2, 0) = a02; A(2, 1) = a12; A(2, 2) = a22;
    Ellipsoid3<double> F{};
    bool valid = F.FromCoefficients(A, B, C);
    io.outBool(valid);
    if (valid)
    {
        io.outVec(F.center);
        io.outVec(F.axis[0]);
        io.outVec(F.axis[1]);
        io.outVec(F.axis[2]);
        io.outVec(F.extent);
    }
}

ORACLE_CASE("Hyperellipsoid.compare.2d")
{
    // Draws from a tiny set so that ties in the center, the axes and the
    // extents are common and every branch of operator< is exercised.
    Ellipse2<double> E0{}, E1{};
    E0.center = io.latticeVec<2>(-1, 1);
    E0.axis[0] = io.latticeVec<2>(-1, 1);
    E0.axis[1] = io.latticeVec<2>(-1, 1);
    E0.extent = io.latticeVec<2>(-1, 1);
    E1.center = io.latticeVec<2>(-1, 1);
    E1.axis[0] = io.latticeVec<2>(-1, 1);
    E1.axis[1] = io.latticeVec<2>(-1, 1);
    E1.extent = io.latticeVec<2>(-1, 1);
    io.outBool(E0 == E1);
    io.outBool(E0 != E1);
    io.outBool(E0 < E1);
    io.outBool(E0 <= E1);
    io.outBool(E0 > E1);
    io.outBool(E0 >= E1);
}

// -------------------------------------------------------------------- Hyperplane

ORACLE_CASE("Hyperplane.construct.3d")
{
    auto normal = io.unit<3>();
    double constant = io.real(-5.0, 5.0);
    auto origin = io.vec<3>(-5.0, 5.0);
    Plane3<double> defaultPlane{};
    Plane3<double> fromConstant(normal, constant);
    Plane3<double> fromOrigin(normal, origin);
    io.outVec(defaultPlane.normal);
    io.outVec(defaultPlane.origin);
    io.outReal(defaultPlane.constant);
    io.outVec(fromConstant.normal);
    io.outVec(fromConstant.origin);
    io.outReal(fromConstant.constant);
    io.outVec(fromOrigin.normal);
    io.outVec(fromOrigin.origin);
    io.outReal(fromOrigin.constant);
}

ORACLE_CASE("Hyperplane.fromPoints.3d")
{
    // The N = 3 specialization is UnitCross, which is arithmetic plus sqrt.
    // Every 3rd record uses a lattice so that collinear (degenerate) and
    // axis-aligned triples occur.
    std::array<Vector3<double>, 3> p{};
    if (io.index() % 3 == 0)
    {
        p[0] = io.latticeVec<3>(-3, 3);
        p[1] = io.latticeVec<3>(-3, 3);
        p[2] = io.latticeVec<3>(-3, 3);
    }
    else
    {
        p[0] = io.vec<3>(-5.0, 5.0);
        p[1] = io.vec<3>(-5.0, 5.0);
        p[2] = io.vec<3>(-5.0, 5.0);
    }
    Plane3<double> plane(p);
    io.outVec(plane.normal);
    io.outVec(plane.origin);
    io.outReal(plane.constant);
}

ORACLE_CASE("Hyperplane.compare.3d")
{
    Hyperplane<3, double> h0{}, h1{};
    h0.normal = io.latticeVec<3>(-1, 1);
    h0.origin = io.latticeVec<3>(-1, 1);
    h0.constant = io.lattice(-1, 1);
    h1.normal = io.latticeVec<3>(-1, 1);
    h1.origin = io.latticeVec<3>(-1, 1);
    h1.constant = io.lattice(-1, 1);
    io.outBool(h0 == h1);
    io.outBool(h0 != h1);
    io.outBool(h0 < h1);
    io.outBool(h0 <= h1);
    io.outBool(h0 > h1);
    io.outBool(h0 >= h1);
}

// Deliberate port deviation, issue #217. Upstream routes N != 3 through
// ComputeFromPoints, which builds SingularValueDecomposition(N, N-1, 32) and
// calls Solve(&edge[0], -1). For N = 2 the SVD constructor's
// LogAssert(mNumCols >= 2) fires; for N >= 4 Solve's LogAssert(multiplier > 0)
// fires. So Hyperplane<N,T> cannot be built from points at all upstream and
// every record throws. The port computes the orthogonal complement.
ORACLE_CASE("Hyperplane.fromPoints.deviation.2d")
{
    std::array<Vector2<double>, 2> p{};
    p[0] = io.latticeVec<2>(-3, 3);
    p[1] = io.latticeVec<2>(-3, 3);
    Hyperplane<2, double> plane(p);
    io.outVec(plane.normal);
    io.outVec(plane.origin);
    io.outReal(plane.constant);
}

ORACLE_CASE("Hyperplane.fromPoints.deviation.4d")
{
    std::array<Vector<4, double>, 4> p{};
    p[0] = io.latticeVec<4>(-3, 3);
    p[1] = io.latticeVec<4>(-3, 3);
    p[2] = io.latticeVec<4>(-3, 3);
    p[3] = io.latticeVec<4>(-3, 3);
    Hyperplane<4, double> plane(p);
    io.outVec(plane.normal);
    io.outVec(plane.origin);
    io.outReal(plane.constant);
}

// -------------------------------------------------------------------------- Cone

ORACLE_CASE("Cone.setAngle")
{
    // cos, sin and tan come from the C runtime here and from V8 in the port,
    // so this case measures libm agreement. The angle itself is exact.
    auto origin = io.vec<3>(-4.0, 4.0);
    auto direction = io.unit<3>();
    double angle = io.real(0.05, 1.5);
    Ray3<double> ray(origin, direction);
    Cone3<double> cone(ray, angle);
    io.outReal(cone.angle);  // replayed with outRealExact
    io.outReal(cone.cosAngle);
    io.outReal(cone.sinAngle);
    io.outReal(cone.tanAngle);
    io.outReal(cone.cosAngleSqr);
    io.outReal(cone.sinAngleSqr);
    io.outReal(cone.invSinAngle);
}

ORACLE_CASE("Cone.heights")
{
    // The four cone types and the height predicates. The angle is recorded as
    // an input and no libm-derived quantity reaches an output, so the case is
    // arithmetic-only.
    auto origin = io.vec<3>(-4.0, 4.0);
    auto direction = io.unit<3>();
    double angle = io.real(0.05, 1.5);
    int32_t mode = io.integer(0, 4);
    double h0 = io.lattice(0, 4);
    double h1 = io.lattice(0, 4);
    auto probe = io.latticeVec<4>(-1, 5);
    Ray3<double> ray(origin, direction);
    Cone3<double> cone(ray, angle);
    if (mode == 1)
    {
        cone.MakeInfiniteTruncatedCone(h0);
    }
    else if (mode == 2)
    {
        cone.MakeFiniteCone(h1 > 0.0 ? h1 : 1.0);
    }
    else if (mode == 3)
    {
        cone.MakeConeFrustum(h0, h0 + (h1 > 0.0 ? h1 : 1.0));
    }
    else if (mode == 4)
    {
        cone.MakeInfiniteCone();
    }
    io.outReal(cone.GetMinHeight());
    io.outReal(cone.GetMaxHeight());
    io.outBool(cone.IsFinite());
    io.outBool(cone.IsInfinite());
    for (int32_t i = 0; i < 4; ++i)
    {
        io.outBool(cone.HeightInRange(probe[i]));
        io.outBool(cone.HeightLessThanMin(probe[i]));
        io.outBool(cone.HeightGreaterThanMax(probe[i]));
    }
}

ORACLE_CASE("Cone.heights.throwParity")
{
    // MakeInfiniteTruncatedCone, MakeFiniteCone and MakeConeFrustum assert
    // their height preconditions, and SetAngle asserts 0 < angle < pi/2. The
    // draws straddle each threshold so both the throwing and the
    // non-throwing sides occur.
    auto origin = io.vec<3>(-2.0, 2.0);
    auto direction = io.unit<3>();
    double angle = io.real(-0.5, 2.0);
    int32_t mode = io.integer(0, 2);
    double h0 = io.lattice(-2, 2);
    double h1 = io.lattice(-2, 2);
    Ray3<double> ray(origin, direction);
    Cone3<double> cone(ray, angle);
    if (mode == 0)
    {
        cone.MakeInfiniteTruncatedCone(h0);
    }
    else if (mode == 1)
    {
        cone.MakeFiniteCone(h1);
    }
    else
    {
        cone.MakeConeFrustum(h0, h1);
    }
    io.outReal(cone.GetMinHeight());
    io.outReal(cone.GetMaxHeight());
}

ORACLE_CASE("Cone.compare")
{
    // Ray, angle and the two heights are drawn from a tiny set so that ties
    // reach every branch of operator<.
    auto origin0 = io.latticeVec<3>(-1, 1);
    auto direction0 = io.latticeVec<3>(-1, 1);
    double angle0 = io.real(0.5, 0.5000000000000001);
    double h0min = io.lattice(0, 1);
    double h0max = io.lattice(2, 3);
    auto origin1 = io.latticeVec<3>(-1, 1);
    auto direction1 = io.latticeVec<3>(-1, 1);
    double angle1 = io.real(0.5, 0.5000000000000001);
    double h1min = io.lattice(0, 1);
    double h1max = io.lattice(2, 3);
    Cone3<double> cone0(Ray3<double>(origin0, direction0), angle0, h0min, h0max);
    Cone3<double> cone1(Ray3<double>(origin1, direction1), angle1, h1min, h1max);
    io.outBool(cone0 == cone1);
    io.outBool(cone0 != cone1);
    io.outBool(cone0 < cone1);
    io.outBool(cone0 <= cone1);
    io.outBool(cone0 > cone1);
    io.outBool(cone0 >= cone1);
}

ORACLE_CASE("Cone.createMesh")
{
    // std::cos/std::sin/std::tan decide the vertex positions, so the reals
    // carry the default tolerance. The mesh structure (numExtra, the unique
    // vertex count and the index array) depends only on exact equality of
    // identically computed vertices, so it is emitted as integers and must
    // match.
    auto origin = io.vec<3>(-2.0, 2.0);
    auto direction = io.unit<3>();
    double angle = io.real(0.2, 0.8);
    double hMin = io.real(0.5, 2.0);
    double hMax = io.real(2.5, 4.0);
    bool inscribed = io.boolean();
    Ray3<double> ray(origin, direction);
    Cone3<double> cone(ray, angle, hMin, hMax);
    std::vector<Vector3<double>> vertices{};
    std::vector<int32_t> indices{};
    cone.CreateMesh(3, inscribed, vertices, indices);
    io.outInt(vertices.size());
    io.outInt(indices.size());
    for (int32_t i : indices) { io.outInt(i); }
    for (auto const& v : vertices) { io.outVec(v); }
}

// ---------------------------------------------------------------------- Polygon2

namespace
{
    // The 16 lattice points on the boundary of [-2,2]^2 in counterclockwise
    // order, used for a convex polygon that has collinear vertex triples.
    std::array<Vector2<double>, 16> const& SquarePerimeter()
    {
        static std::array<Vector2<double>, 16> const s =
        { {
            { -2.0, -2.0 }, { -1.0, -2.0 }, { 0.0, -2.0 }, { 1.0, -2.0 },
            { 2.0, -2.0 }, { 2.0, -1.0 }, { 2.0, 0.0 }, { 2.0, 1.0 },
            { 2.0, 2.0 }, { 1.0, 2.0 }, { 0.0, 2.0 }, { -1.0, 2.0 },
            { -2.0, 2.0 }, { -2.0, 1.0 }, { -2.0, 0.0 }, { -2.0, -1.0 }
        } };
        return s;
    }
}

ORACLE_CASE("Polygon2.queries")
{
    // Five generator modes; each records the same number of doubles.
    //   0: lattice points ordered by angle about their average (simple)
    //   1: uniform points ordered by angle about their average (simple)
    //   2: lattice points in draw order (often self-intersecting)
    //   3: convex, on the boundary of a square, with collinear triples
    //   4: a duplicated index, which fails the constructor
    int32_t n = io.integer(3, 7);
    int32_t mode = io.integer(0, 4);
    std::vector<Vector2<double>> pool(static_cast<size_t>(n) + 1);
    std::vector<int32_t> order(static_cast<size_t>(n) + 1);
    if (mode == 3)
    {
        auto const& perimeter = SquarePerimeter();
        std::vector<int32_t> pick{};
        for (int32_t i = 0; i < 16; ++i) { pick.push_back(i); }
        while (static_cast<int32_t>(pick.size()) > n + 1)
        {
            pick.erase(pick.begin() + io.rawInteger(0,
                static_cast<int32_t>(pick.size()) - 1));
        }
        for (size_t i = 0; i < pool.size(); ++i)
        {
            pool[i] = perimeter[static_cast<size_t>(pick[i])];
        }
    }
    else
    {
        for (size_t i = 0; i < pool.size(); ++i)
        {
            if (mode == 1)
            {
                pool[i][0] = io.raw(-4.0, 4.0);
                pool[i][1] = io.raw(-4.0, 4.0);
            }
            else
            {
                pool[i][0] = static_cast<double>(io.rawInteger(-4, 4));
                pool[i][1] = static_cast<double>(io.rawInteger(-4, 4));
            }
        }
    }
    for (auto const& p : pool) { io.givenVec(p); }

    for (size_t i = 0; i < order.size(); ++i)
    {
        order[i] = static_cast<int32_t>(i);
    }
    if (mode == 0 || mode == 1)
    {
        // Order the first n pool points by angle about their average, which
        // yields a simple star-shaped polygon.
        Vector2<double> average{ 0.0, 0.0 };
        for (int32_t i = 0; i < n; ++i) { average += pool[static_cast<size_t>(i)]; }
        average /= static_cast<double>(n);
        std::stable_sort(order.begin(), order.begin() + n,
            [&pool, &average](int32_t a, int32_t b)
            {
                Vector2<double> da = pool[static_cast<size_t>(a)] - average;
                Vector2<double> db = pool[static_cast<size_t>(b)] - average;
                return std::atan2(da[1], da[0]) < std::atan2(db[1], db[0]);
            });
    }
    std::vector<int32_t> indices(order.begin(), order.begin() + n);
    if (mode == 4)
    {
        indices[static_cast<size_t>(n) - 1] = indices[0];
    }
    for (int32_t i : indices) { io.given(static_cast<double>(i)); }

    Polygon2<double> polygon(pool.data(), n, indices.data(), true);
    io.outBool(static_cast<bool>(polygon));
    io.outInt(polygon.GetVertices().size());
    for (int32_t v : polygon.GetVertices()) { io.outInt(v); }
    io.outInt(polygon.GetIndices().size());
    for (int32_t v : polygon.GetIndices()) { io.outInt(v); }
    io.outBool(polygon.CounterClockwise());
    io.outVec(polygon.ComputeVertexAverage());
    io.outReal(polygon.ComputePerimeterLength());
    io.outReal(polygon.ComputeArea());
    io.outBool(polygon.IsSimple());
    io.outBool(polygon.IsConvex());
}

ORACLE_CASE("Polygon2.queries.clockwise")
{
    // The same queries with counterClockwise = false, which flips the sign
    // used by IsConvexInternal. The vertices are ordered clockwise by angle.
    int32_t n = io.integer(4, 7);
    std::vector<Vector2<double>> pool(static_cast<size_t>(n));
    for (size_t i = 0; i < pool.size(); ++i)
    {
        pool[i][0] = static_cast<double>(io.rawInteger(-4, 4));
        pool[i][1] = static_cast<double>(io.rawInteger(-4, 4));
    }
    for (auto const& p : pool) { io.givenVec(p); }

    std::vector<int32_t> indices(static_cast<size_t>(n));
    for (size_t i = 0; i < indices.size(); ++i)
    {
        indices[i] = static_cast<int32_t>(i);
    }
    Vector2<double> average{ 0.0, 0.0 };
    for (auto const& p : pool) { average += p; }
    average /= static_cast<double>(n);
    std::stable_sort(indices.begin(), indices.end(),
        [&pool, &average](int32_t a, int32_t b)
        {
            Vector2<double> da = pool[static_cast<size_t>(a)] - average;
            Vector2<double> db = pool[static_cast<size_t>(b)] - average;
            return std::atan2(da[1], da[0]) > std::atan2(db[1], db[0]);
        });
    for (int32_t i : indices) { io.given(static_cast<double>(i)); }

    Polygon2<double> polygon(pool.data(), n, indices.data(), false);
    io.outBool(static_cast<bool>(polygon));
    io.outVec(polygon.ComputeVertexAverage());
    io.outReal(polygon.ComputePerimeterLength());
    io.outReal(polygon.ComputeArea());
    io.outBool(polygon.IsSimple());
    io.outBool(polygon.IsConvex());
}

// ------------------------------------------------------------- Tetrahedron3

namespace
{
    // Every 3rd record is a lattice tetrahedron, which produces exact
    // arithmetic, zero-volume (coplanar) configurations and repeated
    // vertices.
    Tetrahedron3<double> MakeTetrahedron(oracle::Ctx& io)
    {
        std::array<Vector3<double>, 4> v{};
        if (io.index() % 3 == 0)
        {
            v[0] = io.latticeVec<3>(-3, 3);
            v[1] = io.latticeVec<3>(-3, 3);
            v[2] = io.latticeVec<3>(-3, 3);
            v[3] = io.latticeVec<3>(-3, 3);
        }
        else
        {
            v[0] = io.vec<3>(-5.0, 5.0);
            v[1] = io.vec<3>(-5.0, 5.0);
            v[2] = io.vec<3>(-5.0, 5.0);
            v[3] = io.vec<3>(-5.0, 5.0);
        }
        return Tetrahedron3<double>(v);
    }
}

ORACLE_CASE("Tetrahedron3.normals")
{
    auto tetra = MakeTetrahedron(io);
    for (size_t i = 0; i < 4; ++i) { io.outVec(tetra.ComputeFaceNormal(i)); }
    for (size_t i = 0; i < 6; ++i) { io.outVec(tetra.ComputeEdgeNormal(i)); }
    for (size_t i = 0; i < 4; ++i) { io.outVec(tetra.ComputeVertexNormal(i)); }
}

ORACLE_CASE("Tetrahedron3.computeCentroid")
{
    auto tetra = MakeTetrahedron(io);
    io.outVec(tetra.ComputeCentroid());
}

ORACLE_CASE("Tetrahedron3.getPlanes")
{
    // Only 'normal' and 'constant' are emitted here: GetPlanes never writes
    // 'origin', which the deviation case below covers (issue #268).
    auto tetra = MakeTetrahedron(io);
    std::array<Plane3<double>, 4> plane{};
    tetra.GetPlanes(plane);
    for (size_t i = 0; i < 4; ++i)
    {
        io.outVec(plane[i].normal);
        io.outReal(plane[i].constant);
    }
}

// Deliberate port deviation, issue #268. GetPlanes assigns 'normal' and
// 'constant' but never 'origin', so the caller's array keeps whatever was in
// it. The array here is filled with a recognizable plane first, exactly as a
// caller reusing a buffer would; upstream leaves those stale origins in
// place, the port returns planes whose three members are consistent.
ORACLE_CASE("Tetrahedron3.getPlanes.deviation")
{
    auto tetra = MakeTetrahedron(io);
    auto staleNormal = io.unit<3>();
    double staleConstant = io.real(-3.0, 3.0);
    std::array<Plane3<double>, 4> plane{};
    for (size_t i = 0; i < 4; ++i)
    {
        plane[i] = Plane3<double>(staleNormal, staleConstant);
    }
    tetra.GetPlanes(plane);
    for (size_t i = 0; i < 4; ++i)
    {
        io.outVec(plane[i].origin);
    }
}

ORACLE_CASE("Tetrahedron3.tables")
{
    size_t face = static_cast<size_t>(io.integer(0, 3));
    size_t edge = static_cast<size_t>(io.integer(0, 5));
    size_t vertex = static_cast<size_t>(io.integer(0, 3));
    for (size_t i : Tetrahedron3<double>::GetFaceIndices(face)) { io.outInt(i); }
    for (size_t i : Tetrahedron3<double>::GetAllFaceIndices()) { io.outInt(i); }
    for (size_t i : Tetrahedron3<double>::GetEdgeIndices(edge)) { io.outInt(i); }
    for (size_t i : Tetrahedron3<double>::GetAllEdgeIndices()) { io.outInt(i); }
    for (size_t i : Tetrahedron3<double>::GetEdgeAugmented(edge)) { io.outInt(i); }
    for (size_t i : Tetrahedron3<double>::GetVertexAugmented(vertex)) { io.outInt(i); }
}

ORACLE_CASE("Tetrahedron3.compare")
{
    std::array<Vector3<double>, 4> a{}, b{};
    a[0] = io.latticeVec<3>(-1, 1);
    a[1] = io.latticeVec<3>(-1, 1);
    a[2] = io.latticeVec<3>(-1, 1);
    a[3] = io.latticeVec<3>(-1, 1);
    b[0] = io.latticeVec<3>(-1, 1);
    b[1] = io.latticeVec<3>(-1, 1);
    b[2] = io.latticeVec<3>(-1, 1);
    b[3] = io.latticeVec<3>(-1, 1);
    Tetrahedron3<double> t0(a), t1(b);
    io.outBool(t0 == t1);
    io.outBool(t0 != t1);
    io.outBool(t0 < t1);
    io.outBool(t0 <= t1);
    io.outBool(t0 > t1);
    io.outBool(t0 >= t1);
}

ORACLE_CASE("Tetrahedron3.defaultConstruct")
{
    // The default tetrahedron and the two value constructors.
    auto v0 = io.latticeVec<3>(-2, 2);
    auto v1 = io.latticeVec<3>(-2, 2);
    auto v2 = io.latticeVec<3>(-2, 2);
    auto v3 = io.latticeVec<3>(-2, 2);
    Tetrahedron3<double> canonical{};
    Tetrahedron3<double> fromFour(v0, v1, v2, v3);
    for (size_t i = 0; i < 4; ++i) { io.outVec(canonical.v[i]); }
    for (size_t i = 0; i < 4; ++i) { io.outVec(fromFour.v[i]); }
    io.outVec(canonical.ComputeCentroid());
}

// --------------------------------------------------------- RectangleManager

namespace
{
    // Small integer rectangles so that shared endpoint values (ties in
    // Endpoint::operator<) and touching rectangles are common.
    void DrawRectangles(oracle::Ctx& io, std::vector<AlignedBox2<double>>& r)
    {
        for (size_t i = 0; i < r.size(); ++i)
        {
            double x0 = io.lattice(-4, 4);
            double y0 = io.lattice(-4, 4);
            double w = io.lattice(0, 4);
            double h = io.lattice(0, 4);
            r[i].min = { x0, y0 };
            r[i].max = { x0 + w, y0 + h };
        }
    }

    void OutOverlap(oracle::Ctx& io, std::set<EdgeKey<false>> const& overlap)
    {
        io.outInt(overlap.size());
        for (auto const& e : overlap)
        {
            io.outInt(e.V[0]);
            io.outInt(e.V[1]);
        }
    }
}

ORACLE_CASE("RectangleManager.initialize")
{
    int32_t n = io.integer(2, 8);
    std::vector<AlignedBox2<double>> rectangles(static_cast<size_t>(n));
    DrawRectangles(io, rectangles);
    RectangleManager<double> manager(rectangles);
    OutOverlap(io, manager.GetOverlap());
    for (int32_t i = 0; i < n; ++i)
    {
        AlignedBox2<double> box{};
        manager.GetRectangle(i, box);
        io.outVec(box.min);
        io.outVec(box.max);
    }
}

ORACLE_CASE("RectangleManager.update")
{
    // Move every rectangle by a small integer offset and re-run the
    // incremental update. The endpoints stay nearly sorted, which is the
    // regime the insertion sort is written for, and the transpositions drive
    // the incremental insert/erase of the overlap set.
    int32_t n = io.integer(2, 6);
    std::vector<AlignedBox2<double>> rectangles(static_cast<size_t>(n));
    DrawRectangles(io, rectangles);
    RectangleManager<double> manager(rectangles);
    OutOverlap(io, manager.GetOverlap());

    for (int32_t pass = 0; pass < 2; ++pass)
    {
        for (int32_t i = 0; i < n; ++i)
        {
            double dx = io.lattice(-2, 2);
            double dy = io.lattice(-2, 2);
            AlignedBox2<double> box{};
            manager.GetRectangle(i, box);
            box.min[0] += dx;
            box.max[0] += dx;
            box.min[1] += dy;
            box.max[1] += dy;
            manager.SetRectangle(i, box);
        }
        manager.Update();
        OutOverlap(io, manager.GetOverlap());
    }
}

// ------------------------------------------------------------ RectangleMesh

namespace
{
    Rectangle3<double> MakeRectangle3(oracle::Ctx& io)
    {
        auto center = io.vec<3>(-3.0, 3.0);
        auto frame = MakeFrame3(io);
        auto extent = MakeExtent<2>(io);
        std::array<Vector3<double>, 2> axis{ frame[0], frame[1] };
        return Rectangle3<double>(center, axis, extent);
    }
}

ORACLE_CASE("RectangleMesh.construct")
{
    // No texture-coordinate channel is supplied, so RectangleMesh allocates
    // the default coordinates; wantDynamicTangentSpaceUpdate is false, so
    // InitializeNormals runs.
    int32_t numRows = io.integer(2, 5);
    int32_t numCols = io.integer(2, 5);
    auto rectangle = MakeRectangle3(io);
    MeshDescription description(MeshTopology::RECTANGLE,
        static_cast<uint32_t>(numRows), static_cast<uint32_t>(numCols));
    std::vector<Vector3<double>> positions(description.numVertices);
    std::vector<Vector3<double>> normals(description.numVertices);
    std::vector<uint32_t> indices(3 * static_cast<size_t>(description.numTriangles));
    description.vertexAttributes.push_back(VertexAttribute("position",
        positions.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("normal",
        normals.data(), sizeof(Vector3<double>)));
    description.indexAttribute = IndexAttribute(indices.data(), sizeof(uint32_t));
    RectangleMesh<double> mesh(description, rectangle);
    io.outInt(description.numVertices);
    io.outInt(description.numTriangles);
    for (auto const& p : positions) { io.outVec(p); }
    for (auto const& nrm : normals) { io.outVec(nrm); }
    for (uint32_t i : indices) { io.outInt(i); }
    io.outVec(mesh.GetRectangle().center);
}

ORACLE_CASE("RectangleMesh.tcoords")
{
    // A client-supplied texture-coordinate channel, which is what
    // InitializeTCoords writes and InitializePositions reads back.
    int32_t numRows = io.integer(2, 5);
    int32_t numCols = io.integer(2, 5);
    auto rectangle = MakeRectangle3(io);
    bool wantCCW = io.boolean();
    MeshDescription description(MeshTopology::RECTANGLE,
        static_cast<uint32_t>(numRows), static_cast<uint32_t>(numCols));
    description.wantCCW = wantCCW;
    std::vector<Vector3<double>> positions(description.numVertices);
    std::vector<Vector2<double>> tcoords(description.numVertices);
    std::vector<uint32_t> indices(3 * static_cast<size_t>(description.numTriangles));
    description.vertexAttributes.push_back(VertexAttribute("position",
        positions.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("tcoord",
        tcoords.data(), sizeof(Vector2<double>)));
    description.indexAttribute = IndexAttribute(indices.data(), sizeof(uint32_t));
    RectangleMesh<double> mesh(description, rectangle);
    for (auto const& t : tcoords) { io.outVec(t); }
    for (auto const& p : positions) { io.outVec(p); }
    for (uint32_t i : indices) { io.outInt(i); }
}

// Deliberate port deviation, issue #268. RectangleMesh::InitializeFrame
// hardcodes tangent = (1,0,0) and bitangent = (0,1,0) for every vertex
// regardless of the rectangle's axes, so the tangent, bitangent, dpdu and
// dpdv channels are not in the rectangle's plane. The port uses the
// rectangle's own orthonormal axes; the normal channel is unaffected and is
// compared in the main case above.
ORACLE_CASE("RectangleMesh.frame.deviation")
{
    int32_t numRows = io.integer(2, 4);
    int32_t numCols = io.integer(2, 4);
    auto rectangle = MakeRectangle3(io);
    MeshDescription description(MeshTopology::RECTANGLE,
        static_cast<uint32_t>(numRows), static_cast<uint32_t>(numCols));
    description.wantDynamicTangentSpaceUpdate = true;
    std::vector<Vector3<double>> positions(description.numVertices);
    std::vector<Vector3<double>> normals(description.numVertices);
    std::vector<Vector3<double>> tangents(description.numVertices);
    std::vector<Vector3<double>> bitangents(description.numVertices);
    std::vector<Vector3<double>> dpdus(description.numVertices);
    std::vector<Vector3<double>> dpdvs(description.numVertices);
    std::vector<uint32_t> indices(3 * static_cast<size_t>(description.numTriangles));
    description.vertexAttributes.push_back(VertexAttribute("position",
        positions.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("normal",
        normals.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("tangent",
        tangents.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("bitangent",
        bitangents.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("dpdu",
        dpdus.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("dpdv",
        dpdvs.data(), sizeof(Vector3<double>)));
    description.indexAttribute = IndexAttribute(indices.data(), sizeof(uint32_t));
    RectangleMesh<double> mesh(description, rectangle);
    io.outBool(description.allowUpdateFrame);
    for (auto const& v : normals) { io.outVec(v); }
    for (auto const& v : tangents) { io.outVec(v); }
    for (auto const& v : bitangents) { io.outVec(v); }
    for (auto const& v : dpdus) { io.outVec(v); }
    for (auto const& v : dpdvs) { io.outVec(v); }
}

// ------------------------------------------------------- RectanglePatchMesh

namespace
{
    // RectanglePatchMesh's only collaborator is an abstract
    // ParametricSurface. A quadratic graph surface keeps the whole mesh
    // construction arithmetic-only (ORACLE.md, the v18 precedent). The
    // expression grouping below is reproduced verbatim in the replay.
    class QuadraticGraphSurface : public ParametricSurface<3, double>
    {
    public:
        QuadraticGraphSurface(double umin, double umax, double vmin,
            double vmax, std::array<double, 6> const& a)
            :
            ParametricSurface<3, double>(umin, umax, vmin, vmax, true),
            mA(a)
        {
            mConstructed = true;
        }

        virtual void Evaluate(double u, double v, uint32_t order,
            Vector3<double>* jet) const override
        {
            double z = ((((mA[0] + mA[1] * u) + mA[2] * v) + mA[3] * (u * u))
                + mA[4] * (u * v)) + mA[5] * (v * v);
            jet[0] = { u, v, z };
            if (order >= 1)
            {
                double zu = (mA[1] + (mA[3] + mA[3]) * u) + mA[4] * v;
                double zv = (mA[2] + mA[4] * u) + (mA[5] + mA[5]) * v;
                jet[1] = { 1.0, 0.0, zu };
                jet[2] = { 0.0, 1.0, zv };
            }
            if (order >= 2)
            {
                jet[3] = { 0.0, 0.0, mA[3] + mA[3] };
                jet[4] = { 0.0, 0.0, mA[4] };
                jet[5] = { 0.0, 0.0, mA[5] + mA[5] };
            }
        }

    private:
        std::array<double, 6> mA;
    };
}

ORACLE_CASE("RectanglePatchMesh.construct")
{
    int32_t numRows = io.integer(2, 5);
    int32_t numCols = io.integer(2, 5);
    double umin = io.real(-2.0, 0.0);
    double umax = io.real(0.5, 2.0);
    double vmin = io.real(-2.0, 0.0);
    double vmax = io.real(0.5, 2.0);
    std::array<double, 6> a{};
    for (size_t i = 0; i < 6; ++i) { a[i] = io.real(-1.5, 1.5); }
    auto surface = std::make_shared<QuadraticGraphSurface>(umin, umax, vmin,
        vmax, a);
    MeshDescription description(MeshTopology::RECTANGLE,
        static_cast<uint32_t>(numRows), static_cast<uint32_t>(numCols));
    std::vector<Vector3<double>> positions(description.numVertices);
    std::vector<Vector3<double>> normals(description.numVertices);
    std::vector<Vector2<double>> tcoords(description.numVertices);
    std::vector<uint32_t> indices(3 * static_cast<size_t>(description.numTriangles));
    description.vertexAttributes.push_back(VertexAttribute("position",
        positions.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("normal",
        normals.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("tcoord",
        tcoords.data(), sizeof(Vector2<double>)));
    description.indexAttribute = IndexAttribute(indices.data(), sizeof(uint32_t));
    RectanglePatchMesh<double> mesh(description, surface);
    for (auto const& t : tcoords) { io.outVec(t); }
    for (auto const& p : positions) { io.outVec(p); }
    for (auto const& nrm : normals) { io.outVec(nrm); }
    for (uint32_t i : indices) { io.outInt(i); }
}

ORACLE_CASE("RectanglePatchMesh.frame")
{
    // wantDynamicTangentSpaceUpdate with a full tangent-space request runs
    // InitializeFrame, which derives the frame from the surface derivatives
    // (RectanglePatchMesh does not share RectangleMesh's hardcoded frame).
    int32_t numRows = io.integer(2, 4);
    int32_t numCols = io.integer(2, 4);
    double umin = io.real(-2.0, 0.0);
    double umax = io.real(0.5, 2.0);
    double vmin = io.real(-2.0, 0.0);
    double vmax = io.real(0.5, 2.0);
    std::array<double, 6> a{};
    for (size_t i = 0; i < 6; ++i) { a[i] = io.real(-1.5, 1.5); }
    auto surface = std::make_shared<QuadraticGraphSurface>(umin, umax, vmin,
        vmax, a);
    MeshDescription description(MeshTopology::RECTANGLE,
        static_cast<uint32_t>(numRows), static_cast<uint32_t>(numCols));
    description.wantDynamicTangentSpaceUpdate = true;
    std::vector<Vector3<double>> positions(description.numVertices);
    std::vector<Vector3<double>> normals(description.numVertices);
    std::vector<Vector3<double>> tangents(description.numVertices);
    std::vector<Vector3<double>> bitangents(description.numVertices);
    std::vector<Vector3<double>> dpdus(description.numVertices);
    std::vector<Vector3<double>> dpdvs(description.numVertices);
    std::vector<uint32_t> indices(3 * static_cast<size_t>(description.numTriangles));
    description.vertexAttributes.push_back(VertexAttribute("position",
        positions.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("normal",
        normals.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("tangent",
        tangents.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("bitangent",
        bitangents.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("dpdu",
        dpdus.data(), sizeof(Vector3<double>)));
    description.vertexAttributes.push_back(VertexAttribute("dpdv",
        dpdvs.data(), sizeof(Vector3<double>)));
    description.indexAttribute = IndexAttribute(indices.data(), sizeof(uint32_t));
    RectanglePatchMesh<double> mesh(description, surface);
    io.outBool(description.allowUpdateFrame);
    for (auto const& v : normals) { io.outVec(v); }
    for (auto const& v : tangents) { io.outVec(v); }
    for (auto const& v : bitangents) { io.outVec(v); }
    for (auto const& v : dpdus) { io.outVec(v); }
    for (auto const& v : dpdvs) { io.outVec(v); }
}

// ------------------------------------------- AlignedBoxBV and OrientedBoxBV

ORACLE_CASE("AlignedBoxBV.queries")
{
    // useDefault keeps the default-constructed box, which is [-1,1]^3.
    bool useDefault = io.boolean();
    auto boxMin = io.latticeVec<3>(-3, 3);
    auto boxExtent = io.latticeVec<3>(0, 3);
    auto P = io.vec<3>(-5.0, 5.0);
    auto D = io.unit<3>();
    auto R = io.vec<3>(-5.0, 5.0);
    AlignedBoxBV<double> bv{};
    if (!useDefault)
    {
        bv.box.min = boxMin;
        bv.box.max = boxMin + boxExtent;
    }
    Vector3<double> origin{}, direction{};
    bv.GetSplittingAxis(origin, direction);
    io.outVec(bv.box.min);
    io.outVec(bv.box.max);
    io.outVec(origin);
    io.outVec(direction);
    io.outBool(AlignedBoxBV<double>::IntersectLine(P, D, bv));
    io.outBool(AlignedBoxBV<double>::IntersectRay(P, D, bv));
    io.outBool(AlignedBoxBV<double>::IntersectSegment(P, R, bv));
}

ORACLE_CASE("OrientedBoxBV.queries")
{
    bool useDefault = io.boolean();
    auto center = io.vec<3>(-3.0, 3.0);
    auto frame = MakeFrame3(io);
    auto extent = MakeExtent<3>(io);
    auto P = io.vec<3>(-5.0, 5.0);
    auto D = io.unit<3>();
    auto R = io.vec<3>(-5.0, 5.0);
    OrientedBoxBV<double> bv{};
    if (!useDefault)
    {
        bv.box.center = center;
        bv.box.axis = frame;
        bv.box.extent = extent;
    }
    Vector3<double> origin{}, direction{};
    bv.GetSplittingAxis(origin, direction);
    io.outVec(bv.box.center);
    io.outVec(bv.box.axis[0]);
    io.outVec(bv.box.axis[1]);
    io.outVec(bv.box.axis[2]);
    io.outVec(bv.box.extent);
    io.outVec(origin);
    io.outVec(direction);
    io.outBool(OrientedBoxBV<double>::IntersectLine(P, D, bv));
    io.outBool(OrientedBoxBV<double>::IntersectRay(P, D, bv));
    io.outBool(OrientedBoxBV<double>::IntersectSegment(P, R, bv));
}

// -------------------------------------------------- Bounding-volume trees

namespace
{
    std::int64_t NodeIndexOut(std::size_t i)
    {
        return i == std::numeric_limits<std::size_t>::max()
            ? -1 : static_cast<std::int64_t>(i);
    }

    // Vertices for the tree primitives. Every mode records 3 * count
    // doubles.
    //   0: small lattice     3: collinear on a lattice line
    //   1: uniform           4: coplanar (z = 0) lattice
    //   2: all equal         5: two distinct lattice values (duplicates)
    std::vector<Vector3<double>> MakeVertices(oracle::Ctx& io, int32_t count,
        int32_t mode)
    {
        std::vector<Vector3<double>> v(static_cast<size_t>(count));
        Vector3<double> base{}, step{};
        for (int32_t i = 0; i < 3; ++i)
        {
            base[i] = static_cast<double>(io.rawInteger(-3, 3));
            step[i] = static_cast<double>(io.rawInteger(-2, 2));
        }
        for (size_t i = 0; i < v.size(); ++i)
        {
            switch (mode)
            {
            case 0:
                for (int32_t k = 0; k < 3; ++k)
                {
                    v[i][k] = static_cast<double>(io.rawInteger(-4, 4));
                }
                break;
            case 1:
                for (int32_t k = 0; k < 3; ++k) { v[i][k] = io.raw(-5.0, 5.0); }
                break;
            case 2:
                v[i] = base;
                break;
            case 3:
                v[i] = base + static_cast<double>(io.rawInteger(-3, 3)) * step;
                break;
            case 4:
                v[i][0] = static_cast<double>(io.rawInteger(-4, 4));
                v[i][1] = static_cast<double>(io.rawInteger(-4, 4));
                v[i][2] = 0.0;
                break;
            default:
                for (int32_t k = 0; k < 3; ++k)
                {
                    v[i][k] = (io.rawInteger(0, 1) == 0 ? base[k] : base[k] + step[k]);
                }
                break;
            }
        }
        for (auto const& p : v) { io.givenVec(p); }
        return v;
    }

    struct Linear
    {
        Vector3<double> P, D, R;
    };

    // The linear component for the tree queries. Half the records aim
    // through a primitive vertex so the hit rate is high. Records 9 doubles.
    Linear MakeLinear(oracle::Ctx& io, std::vector<Vector3<double>> const& targets)
    {
        Linear L{};
        Vector3<double> d{};
        double len = 0.0;
        do
        {
            for (int32_t i = 0; i < 3; ++i) { d[i] = io.raw(-1.0, 1.0); }
            len = Length(d);
        }
        while (len < 0.25 || len > 1.0);
        Normalize(d);
        if (io.rawInteger(0, 1) == 0 && !targets.empty())
        {
            auto const& c = targets[static_cast<size_t>(io.rawInteger(0,
                static_cast<int32_t>(targets.size()) - 1))];
            L.P = c - 3.0 * d;
            L.R = c + 3.0 * d;
        }
        else
        {
            for (int32_t i = 0; i < 3; ++i) { L.P[i] = io.raw(-6.0, 6.0); }
            for (int32_t i = 0; i < 3; ++i) { L.R[i] = io.raw(-6.0, 6.0); }
            d = L.R - L.P;
            Normalize(d);
        }
        L.D = d;
        io.givenVec(L.P);
        io.givenVec(L.D);
        io.givenVec(L.R);
        return L;
    }

    template <typename Tree>
    void OutTreeStructure(oracle::Ctx& io, Tree const& tree)
    {
        io.outInt(tree.GetHeight());
        auto const& partition = tree.GetPartition();
        io.outInt(partition.size());
        for (auto p : partition) { io.outInt(p); }
        auto const& nodes = tree.GetNodes();
        io.outInt(nodes.size());
        for (auto const& n : nodes)
        {
            io.outInt(NodeIndexOut(n.minIndex));
            io.outInt(NodeIndexOut(n.maxIndex));
            io.outInt(NodeIndexOut(n.leftChild));
            io.outInt(NodeIndexOut(n.rightChild));
        }
    }

    template <typename Tree>
    void OutAlignedBoxes(oracle::Ctx& io, Tree const& tree)
    {
        for (auto const& n : tree.GetNodes())
        {
            io.outVec(n.boundingVolume.box.min);
            io.outVec(n.boundingVolume.box.max);
        }
    }

    // 'skipLeafExtent' omits the extents of leaf nodes, which
    // OrientedBoxTreeOfTriangles collapses along the wrong axis upstream
    // (issue #343); its own deviation case compares them.
    template <typename Tree>
    void OutOrientedBoxes(oracle::Ctx& io, Tree const& tree, bool skipLeafExtent)
    {
        for (auto const& n : tree.GetNodes())
        {
            auto const& box = n.boundingVolume.box;
            io.outVec(box.center);
            io.outVec(box.axis[0]);
            io.outVec(box.axis[1]);
            io.outVec(box.axis[2]);
            bool isLeaf = (n.leftChild == Tree::Node::invalid);
            if (!(skipLeafExtent && isLeaf)) { io.outVec(box.extent); }
        }
    }

    template <typename Tree>
    void OutTreeQueries(oracle::Ctx& io, Tree& tree, Linear const& L)
    {
        std::vector<std::size_t> nodeIndices{};
        tree.Execute(Tree::LINE_QUERY, L.P, L.D, nodeIndices);
        io.outInt(nodeIndices.size());
        for (auto i : nodeIndices) { io.outInt(i); }
        tree.Execute(Tree::RAY_QUERY, L.P, L.D, nodeIndices);
        io.outInt(nodeIndices.size());
        for (auto i : nodeIndices) { io.outInt(i); }
        tree.Execute(Tree::SEGMENT_QUERY, L.P, L.R, nodeIndices);
        io.outInt(nodeIndices.size());
        for (auto i : nodeIndices) { io.outInt(i); }
    }
}

namespace
{
    std::vector<std::array<std::size_t, 2>> MakeSegments(oracle::Ctx& io,
        int32_t numSegments, int32_t numVertices)
    {
        std::vector<std::array<std::size_t, 2>> segments(
            static_cast<size_t>(numSegments));
        for (auto& s : segments)
        {
            s[0] = static_cast<size_t>(io.rawInteger(0, numVertices - 1));
            s[1] = static_cast<size_t>(io.rawInteger(0, numVertices - 1));
        }
        for (auto const& s : segments)
        {
            io.given(static_cast<double>(s[0]));
            io.given(static_cast<double>(s[1]));
        }
        return segments;
    }

    // Distinct index triples: two triangles with the same three indices are
    // hit at bit-identical parameters, which upstream's Intersection set
    // collapses (issue #167). The deviation case below is aimed at exactly
    // that; the main cases keep the triples distinct.
    std::vector<std::array<std::size_t, 3>> MakeTriangles(oracle::Ctx& io,
        int32_t numTriangles, int32_t numVertices)
    {
        std::vector<std::array<std::size_t, 3>> triangles{};
        for (int32_t attempt = 0;
            attempt < 200 && static_cast<int32_t>(triangles.size()) < numTriangles;
            ++attempt)
        {
            std::array<std::size_t, 3> t{};
            t[0] = static_cast<size_t>(io.rawInteger(0, numVertices - 1));
            t[1] = static_cast<size_t>(io.rawInteger(0, numVertices - 1));
            t[2] = static_cast<size_t>(io.rawInteger(0, numVertices - 1));
            std::array<std::size_t, 3> sorted = t;
            std::sort(sorted.begin(), sorted.end());
            bool duplicate = false;
            for (auto const& u : triangles)
            {
                std::array<std::size_t, 3> other = u;
                std::sort(other.begin(), other.end());
                if (other == sorted) { duplicate = true; break; }
            }
            if (!duplicate) { triangles.push_back(t); }
        }
        while (static_cast<int32_t>(triangles.size()) < numTriangles)
        {
            // Fallback when the vertex count is too small to supply enough
            // distinct triples; the records are still comparable.
            triangles.push_back(triangles.back());
        }
        for (auto const& t : triangles)
        {
            io.given(static_cast<double>(t[0]));
            io.given(static_cast<double>(t[1]));
            io.given(static_cast<double>(t[2]));
        }
        return triangles;
    }

    // True when two triangles are hit at bit-identical parameters, the
    // configuration upstream's std::set<Intersection> collapses. The probe is
    // an exact superset of the reported hits (it tests every triangle, not
    // only the ones in the reported leaves), which is allowed: it can only
    // narrow the generator.
    bool HasCoincidentParameters(std::uint32_t queryType,
        Vector3<double> const& A, Vector3<double> const& B,
        std::vector<Vector3<double>> const& V,
        std::vector<std::array<std::size_t, 3>> const& T)
    {
        std::vector<double> params{};
        for (auto const& t : T)
        {
            Triangle3<double> tri(V[t[0]], V[t[1]], V[t[2]]);
            if (queryType == 0)
            {
                FIQuery<double, Line3<double>, Triangle3<double>> q{};
                auto r = q(Line3<double>(A, B), tri);
                if (r.intersect) { params.push_back(r.parameter); }
            }
            else if (queryType == 1)
            {
                FIQuery<double, Ray3<double>, Triangle3<double>> q{};
                auto r = q(Ray3<double>(A, B), tri);
                if (r.intersect) { params.push_back(r.parameter); }
            }
            else
            {
                FIQuery<double, Segment3<double>, Triangle3<double>> q{};
                auto r = q(Segment3<double>(A, B), tri);
                if (r.intersect) { params.push_back(r.parameter); }
            }
        }
        for (size_t i = 0; i < params.size(); ++i)
        {
            for (size_t j = i + 1; j < params.size(); ++j)
            {
                if (params[i] == params[j]) { return true; }
            }
        }
        return false;
    }

    template <typename Tree>
    void OutTriangleTreeQueries(oracle::Ctx& io, Tree& tree, Linear const& L,
        std::vector<Vector3<double>> const& V,
        std::vector<std::array<std::size_t, 3>> const& T)
    {
        for (std::uint32_t queryType = 0; queryType < 3; ++queryType)
        {
            Vector3<double> A = L.P;
            Vector3<double> B = (queryType == 2 ? L.R : L.D);
            std::vector<std::size_t> nodeIndices{};
            std::set<typename Tree::Intersection> intersections{};
            tree.Execute(queryType, A, B, nodeIndices, intersections);
            io.outInt(nodeIndices.size());
            for (auto i : nodeIndices) { io.outInt(i); }
            bool coincident = HasCoincidentParameters(queryType, A, B, V, T);
            io.outBool(coincident);
            if (!coincident)
            {
                io.outInt(intersections.size());
                for (auto const& it : intersections)
                {
                    io.outInt(it.triangleIndex);
                    io.outReal(it.parameter);
                    io.outVec(it.point);
                }
            }
        }
    }
}

ORACLE_CASE("AlignedBoxTreeOfPoints.create")
{
    int32_t count = io.integer(1, 8);
    int32_t mode = io.integer(0, 5);
    int32_t height = io.integer(-1, 4);
    auto vertices = MakeVertices(io, count, mode);
    auto L = MakeLinear(io, vertices);
    AlignedBoxTreeOfPoints<double> tree{};
    if (height < 0)
    {
        tree.Create(vertices);
    }
    else
    {
        tree.Create(vertices, static_cast<std::size_t>(height));
    }
    OutTreeStructure(io, tree);
    OutAlignedBoxes(io, tree);
    OutTreeQueries(io, tree, L);
}

ORACLE_CASE("OrientedBoxTreeOfPoints.create")
{
    int32_t count = io.integer(1, 8);
    int32_t mode = io.integer(0, 5);
    int32_t height = io.integer(-1, 4);
    auto vertices = MakeVertices(io, count, mode);
    auto L = MakeLinear(io, vertices);
    OrientedBoxTreeOfPoints<double> tree{};
    if (height < 0)
    {
        tree.Create(vertices);
    }
    else
    {
        tree.Create(vertices, static_cast<std::size_t>(height));
    }
    OutTreeStructure(io, tree);
    OutOrientedBoxes(io, tree, false);
    OutTreeQueries(io, tree, L);
}

ORACLE_CASE("AlignedBoxTreeOfSegments.create")
{
    int32_t numVertices = io.integer(2, 8);
    int32_t numSegments = io.integer(1, 6);
    int32_t mode = io.integer(0, 5);
    int32_t height = io.integer(-1, 3);
    auto vertices = MakeVertices(io, numVertices, mode);
    auto segments = MakeSegments(io, numSegments, numVertices);
    auto L = MakeLinear(io, vertices);
    AlignedBoxTreeOfSegments<double> tree{};
    if (height < 0)
    {
        tree.Create(vertices, segments);
    }
    else
    {
        tree.Create(vertices, segments, static_cast<std::size_t>(height));
    }
    OutTreeStructure(io, tree);
    OutAlignedBoxes(io, tree);
    OutTreeQueries(io, tree, L);
    for (auto const& c : tree.GetCentroids()) { io.outVec(c); }
}

ORACLE_CASE("OrientedBoxTreeOfSegments.create")
{
    int32_t numVertices = io.integer(2, 8);
    int32_t numSegments = io.integer(1, 6);
    int32_t mode = io.integer(0, 5);
    int32_t height = io.integer(-1, 3);
    auto vertices = MakeVertices(io, numVertices, mode);
    auto segments = MakeSegments(io, numSegments, numVertices);
    auto L = MakeLinear(io, vertices);
    OrientedBoxTreeOfSegments<double> tree{};
    if (height < 0)
    {
        tree.Create(vertices, segments);
    }
    else
    {
        tree.Create(vertices, segments, static_cast<std::size_t>(height));
    }
    OutTreeStructure(io, tree);
    OutOrientedBoxes(io, tree, false);
    OutTreeQueries(io, tree, L);
    for (auto const& c : tree.GetCentroids()) { io.outVec(c); }
}

ORACLE_CASE("AlignedBoxTreeOfTriangles.create")
{
    int32_t numVertices = io.integer(3, 8);
    int32_t numTriangles = io.integer(1, 6);
    int32_t mode = io.integer(0, 5);
    int32_t height = io.integer(-1, 3);
    auto vertices = MakeVertices(io, numVertices, mode);
    auto triangles = MakeTriangles(io, numTriangles, numVertices);
    auto L = MakeLinear(io, vertices);
    AlignedBoxTreeOfTriangles<double> tree{};
    if (height < 0)
    {
        tree.Create(vertices, triangles);
    }
    else
    {
        tree.Create(vertices, triangles, static_cast<std::size_t>(height));
    }
    OutTreeStructure(io, tree);
    OutAlignedBoxes(io, tree);
    OutTriangleTreeQueries(io, tree, L, vertices, triangles);
    for (auto const& c : tree.GetCentroids()) { io.outVec(c); }
}

ORACLE_CASE("OrientedBoxTreeOfTriangles.create")
{
    // The extents of the leaf boxes are omitted here: upstream collapses a
    // leaf box along its largest axis (issue #343). The deviation case below
    // compares them.
    int32_t numVertices = io.integer(3, 8);
    int32_t numTriangles = io.integer(1, 6);
    int32_t mode = io.integer(0, 5);
    int32_t height = io.integer(-1, 3);
    auto vertices = MakeVertices(io, numVertices, mode);
    auto triangles = MakeTriangles(io, numTriangles, numVertices);
    auto L = MakeLinear(io, vertices);
    OrientedBoxTreeOfTriangles<double> tree{};
    if (height < 0)
    {
        tree.Create(vertices, triangles);
    }
    else
    {
        tree.Create(vertices, triangles, static_cast<std::size_t>(height));
    }
    OutTreeStructure(io, tree);
    OutOrientedBoxes(io, tree, true);
    OutTriangleTreeQueries(io, tree, L, vertices, triangles);
    for (auto const& c : tree.GetCentroids()) { io.outVec(c); }
}

// Deliberate port deviation, issue #343.
// OrientedBoxTreeOfTriangles::ComputeLeafBoundingVolume scans for the
// smallest |extent| but its last comparison is 'absExtent > minAbsExtent',
// so the *largest* extent is selected whenever extent[2] is not the smallest
// and the leaf box is collapsed along its longest axis. The port uses '<'.
ORACLE_CASE("OrientedBoxTreeOfTriangles.leafExtent.deviation")
{
    int32_t numVertices = io.integer(3, 8);
    int32_t numTriangles = io.integer(1, 6);
    int32_t mode = io.integer(0, 5);
    auto vertices = MakeVertices(io, numVertices, mode);
    auto triangles = MakeTriangles(io, numTriangles, numVertices);
    OrientedBoxTreeOfTriangles<double> tree{};
    tree.Create(vertices, triangles);
    for (auto const& n : tree.GetNodes())
    {
        if (n.leftChild == OrientedBoxTreeOfTriangles<double>::Node::invalid)
        {
            io.outVec(n.boundingVolume.box.extent);
        }
    }
}

// Deliberate port deviation, issue #167. BVTreeOfTriangles::Execute collects
// hits in a std::set<Intersection> whose operator< compares only the
// parameter, so two triangles hit at the same parameter are set-equivalent
// and all but one are dropped. Two triangles with the same three vertex
// indices are hit at bit-identical parameters. The port orders by
// (parameter, triangleIndex) and keeps both.
ORACLE_CASE("BVTreeOfTriangles.coincident.deviation")
{
    int32_t numVertices = io.integer(3, 6);
    auto vertices = MakeVertices(io, numVertices, 0);
    int32_t i0 = io.integer(0, numVertices - 1);
    int32_t i1 = io.integer(0, numVertices - 1);
    int32_t i2 = io.integer(0, numVertices - 1);
    std::vector<std::array<std::size_t, 3>> triangles(2);
    triangles[0] = { static_cast<size_t>(i0), static_cast<size_t>(i1),
        static_cast<size_t>(i2) };
    triangles[1] = triangles[0];
    // Aim the line through the triangle's centroid.
    Vector3<double> centroid = (vertices[triangles[0][0]]
        + vertices[triangles[0][1]] + vertices[triangles[0][2]]) / 3.0;
    Vector3<double> d{};
    double len = 0.0;
    do
    {
        for (int32_t i = 0; i < 3; ++i) { d[i] = io.raw(-1.0, 1.0); }
        len = Length(d);
    }
    while (len < 0.25 || len > 1.0);
    Normalize(d);
    Vector3<double> P = centroid - 3.0 * d;
    io.givenVec(P);
    io.givenVec(d);
    AlignedBoxTreeOfTriangles<double> tree{};
    tree.Create(vertices, triangles);
    std::vector<std::size_t> nodeIndices{};
    std::set<AlignedBoxTreeOfTriangles<double>::Intersection> intersections{};
    tree.Execute(AlignedBoxTreeOfTriangles<double>::LINE_QUERY, P, d,
        nodeIndices, intersections);
    io.outInt(intersections.size());
    for (auto const& it : intersections)
    {
        io.outInt(it.triangleIndex);
        io.outReal(it.parameter);
        io.outVec(it.point);
    }
}

// __END__
