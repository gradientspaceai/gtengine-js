// Smoke family: a few cases across categories that exercise every part of
// the harness (random and lattice generators, variable-length outputs,
// libm-dependent code, iteration counts).
#define ORACLE_FAMILY "smoke"
#include "Oracle.h"

#include <Mathematics/DistSegmentSegment.h>
#include <Mathematics/DistPointTriangle.h>
#include <Mathematics/IntrLine3Sphere3.h>
#include <Mathematics/IntrRay3Sphere3.h>
#include <Mathematics/SymmetricEigensolver3x3.h>
#include <Mathematics/Rotation.h>

using namespace gte;

namespace
{
    Segment3<double> MakeSegment3(oracle::Ctx& io)
    {
        // Every 3rd record is on a small lattice so that parallel,
        // intersecting and degenerate configurations occur.
        Segment3<double> s;
        if (io.index() % 3 == 0)
        {
            s.p[0] = io.latticeVec<3>(-2, 2);
            s.p[1] = io.latticeVec<3>(-2, 2);
        }
        else
        {
            s.p[0] = io.vec<3>(-10.0, 10.0);
            s.p[1] = io.vec<3>(-10.0, 10.0);
        }
        return s;
    }
}

ORACLE_CASE("DistSegmentSegment.compute.3d")
{
    auto s0 = MakeSegment3(io);
    auto s1 = MakeSegment3(io);
    DCPQuery<double, Segment3<double>, Segment3<double>> query;
    auto r = query(s0, s1);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

ORACLE_CASE("DistSegmentSegment.computeRobust.3d")
{
    auto s0 = MakeSegment3(io);
    auto s1 = MakeSegment3(io);
    DCPQuery<double, Segment3<double>, Segment3<double>> query;
    auto r = query.ComputeRobust(s0, s1);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    io.outReal(r.parameter[0]);
    io.outReal(r.parameter[1]);
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

ORACLE_CASE("DistPointTriangle.compute.3d")
{
    Vector3<double> point;
    Triangle3<double> tri;
    if (io.index() % 3 == 0)
    {
        point = io.latticeVec<3>(-3, 3);
        for (int i = 0; i < 3; ++i) { tri.v[i] = io.latticeVec<3>(-3, 3); }
    }
    else
    {
        point = io.vec<3>(-10.0, 10.0);
        for (int i = 0; i < 3; ++i) { tri.v[i] = io.vec<3>(-10.0, 10.0); }
    }
    DCPQuery<double, Vector3<double>, Triangle3<double>> query;
    auto r = query(point, tri);
    io.outReal(r.distance);
    io.outReal(r.sqrDistance);
    for (int i = 0; i < 3; ++i) { io.outReal(r.barycentric[i]); }
    io.outVec(r.closest[0]);
    io.outVec(r.closest[1]);
}

ORACLE_CASE("IntrLine3Sphere3.find")
{
    // One io draw per statement: C++ leaves the evaluation order of function
    // arguments unspecified (MSVC goes right to left).
    auto origin = io.vec<3>(-3.0, 3.0);
    auto direction = io.unit<3>();
    auto center = io.vec<3>(-3.0, 3.0);
    double radius = io.real(0.5, 4.0);
    Line3<double> line(origin, direction);
    Sphere3<double> sphere(center, radius);
    FIQuery<double, Line3<double>, Sphere3<double>> query;
    auto r = query(line, sphere);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

ORACLE_CASE("IntrRay3Sphere3.find")
{
    auto origin = io.vec<3>(-3.0, 3.0);
    auto direction = io.unit<3>();
    auto center = io.vec<3>(-3.0, 3.0);
    double radius = io.real(0.5, 4.0);
    Ray3<double> ray(origin, direction);
    Sphere3<double> sphere(center, radius);
    FIQuery<double, Ray3<double>, Sphere3<double>> query;
    auto r = query(ray, sphere);
    io.outBool(r.intersect);
    io.outInt(r.numIntersections);
    for (int i = 0; i < r.numIntersections; ++i)
    {
        io.outReal(r.parameter[i]);
        io.outVec(r.point[i]);
    }
}

ORACLE_CASE("SymmetricEigensolver3x3.solve")
{
    double a00 = io.real(-5.0, 5.0), a01 = io.real(-5.0, 5.0), a02 = io.real(-5.0, 5.0);
    double a11 = io.real(-5.0, 5.0), a12 = io.real(-5.0, 5.0), a22 = io.real(-5.0, 5.0);
    bool aggressive = io.boolean();
    int sortType = io.integer(-1, 1);
    std::array<double, 3> eval{};
    std::array<std::array<double, 3>, 3> evec{};
    SymmetricEigensolver3x3<double> solver;
    int32_t iterations = solver(a00, a01, a02, a11, a12, a22, aggressive, sortType, eval, evec);
    io.outInt(iterations);
    for (int i = 0; i < 3; ++i) { io.outReal(eval[i]); }
    for (int i = 0; i < 3; ++i) { for (int j = 0; j < 3; ++j) { io.outReal(evec[i][j]); } }
}

ORACLE_CASE("Rotation.axisAngleToMatrixAndQuaternion")
{
    // sin/cos come from the C runtime here and from V8 in the port, so this
    // case measures libm agreement rather than bit-exact arithmetic.
    auto axis = io.unit<3>();
    double angle = io.real(-3.0, 3.0);
    AxisAngle<3, double> aa(axis, angle);
    Rotation<3, double> rotation(aa);
    Matrix<3, 3, double> m = rotation;
    Quaternion<double> q = rotation;
    io.outMat(m);
    for (int i = 0; i < 4; ++i) { io.outReal(q[i]); }
}
