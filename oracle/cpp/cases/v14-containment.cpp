// Verify group 14 (containment): differential cases for the Cont* headers
// listed in plan/verify-groups.json group 14.
//
// Everything in this group is + - * / sqrt fabs and comparisons, including
// the covariance fits (ApprGaussian2/3, ApprOrthogonalLine3) and the
// SymmetricEigensolver2x2 / SymmetricEigensolver3x3 they run, so every
// ordinary case is declared exact on the TypeScript side. The trigonometry
// that Cone::SetAngle derives from the angle is applied to an *unrecorded*
// draw and the six derived values are recorded as inputs, so the C math
// library never enters a compared computation; the same holds for the sin/cos
// used to build orthonormal frames, unit directions and arc endpoints.
//
// Generators alternate an exactly representable small-lattice mode, a uniform
// mode and constructed degenerate modes (all points equal, collinear,
// coplanar, points exactly on a container boundary, polygon queries through a
// vertex and along an edge, nested / disjoint / identical inputs for the merge
// functions), so the rank-deficient and failure branches are reached.
//
// Frames are built RIGHT-HANDED (2D: axis1 = -Perp(axis0), because GTE's
// Perp(x,y) = (y,-x) is the clockwise one; 3D: axis2 = Cross(axis0, axis1),
// det = +1). Separate ".leftHanded" cases feed det = -1 frames to the
// oriented-box, ellipse and ellipsoid queries on purpose, to see whether any
// of them silently depends on handedness.
#define ORACLE_FAMILY "v14-containment"
#include "Oracle.h"

#include <Mathematics/ContAlignedBox.h>
#include <Mathematics/ContAlignedBox2Arc2.h>
#include <Mathematics/ContCapsule3.h>
#include <Mathematics/ContCircle2.h>
#include <Mathematics/ContCone.h>
#include <Mathematics/ContCylinder3.h>
#include <Mathematics/ContEllipse2.h>
#include <Mathematics/ContEllipse2MinCR.h>
#include <Mathematics/ContEllipsoid3.h>
#include <Mathematics/ContEllipsoid3MinCR.h>
#include <Mathematics/ContLozenge3.h>
#include <Mathematics/ContOrientedBox2.h>
#include <Mathematics/ContOrientedBox3.h>
#include <Mathematics/ContPointInPolygon2.h>
#include <Mathematics/ContPointInPolyhedron3.h>
#include <Mathematics/ContScribeCircle2.h>
#include <Mathematics/ContScribeCircle3Sphere3.h>
#include <Mathematics/ContSphere3.h>
#include <Mathematics/ContTetrahedron3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <limits>
#include <random>
#include <vector>

using namespace gte;

namespace
{
    // ---- shared generators -------------------------------------------------
    // Every helper records a fixed number of doubles in every mode, so the
    // TypeScript replay reads a fixed layout.

    // mode 0 is an exactly representable lattice point, mode 1 is uniform.
    template <int N>
    Vector<N, double> Pt(oracle::Ctx& io, int mode, int lat, double range)
    {
        return (mode == 0 ? io.latticeVec<N>(-lat, lat) : io.vec<N>(-range, range));
    }

    // A point cloud of 'count' points; records count*N doubles in every kind.
    //   0 lattice cloud
    //   1 uniform cloud
    //   2 every point identical (rank 0)
    //   3 collinear on a lattice line (rank 1)
    //   4 coplanar on a lattice plane (rank 2); for N = 2 this is kind 3
    template <int N>
    std::vector<Vector<N, double>> Cloud(oracle::Ctx& io, int kind, int count,
        int lat, double range)
    {
        std::vector<Vector<N, double>> points(static_cast<size_t>(count));
        if (kind == 0 || kind == 1)
        {
            for (int i = 0; i < count; ++i)
            {
                points[static_cast<size_t>(i)] = Pt<N>(io, kind, lat, range);
            }
            return points;
        }

        Vector<N, double> base{}, u{}, v{};
        for (int k = 0; k < N; ++k)
        {
            base[k] = static_cast<double>(io.rawInteger(-lat, lat));
            u[k] = static_cast<double>(io.rawInteger(-lat, lat));
            v[k] = static_cast<double>(io.rawInteger(-lat, lat));
        }
        for (int i = 0; i < count; ++i)
        {
            Vector<N, double> p = base;
            if (kind == 3 || (kind == 4 && N < 3))
            {
                double s = static_cast<double>(io.rawInteger(-lat, lat));
                p = base + s * u;
            }
            else if (kind == 4)
            {
                double s = static_cast<double>(io.rawInteger(-lat, lat));
                double t = static_cast<double>(io.rawInteger(-lat, lat));
                p = base + s * u + t * v;
            }
            points[static_cast<size_t>(i)] = io.givenVec(p);
        }
        return points;
    }

    // Right-handed orthonormal 2D frame (det = +1): axis1 = -Perp(axis0),
    // because GTE's Perp(x,y) = (y,-x) is the clockwise one. Unrecorded.
    std::array<Vector2<double>, 2> DrawFrame2(oracle::Ctx& io, int mode)
    {
        Vector2<double> a0{};
        if (mode == 0)
        {
            int k = io.rawInteger(0, 3);
            a0.MakeZero();
            a0[k % 2] = (k < 2 ? 1.0 : -1.0);
        }
        else
        {
            double len = 0.0;
            do
            {
                a0[0] = io.raw(-1.0, 1.0);
                a0[1] = io.raw(-1.0, 1.0);
                len = Length(a0);
            } while (len < 0.1 || len > 1.0);
            Normalize(a0);
        }
        std::array<Vector2<double>, 2> axis{};
        axis[0] = a0;
        axis[1] = -Perp(a0);
        return axis;
    }

    // Right-handed orthonormal 3D frame (det = +1). Unrecorded.
    std::array<Vector3<double>, 3> DrawFrame3(oracle::Ctx& io, int mode)
    {
        std::array<Vector3<double>, 3> v{};
        if (mode == 0)
        {
            int const perm[6][3] =
            {
                { 0, 1, 2 }, { 0, 2, 1 }, { 1, 0, 2 },
                { 1, 2, 0 }, { 2, 0, 1 }, { 2, 1, 0 }
            };
            int p = io.rawInteger(0, 5);
            for (int i = 0; i < 2; ++i)
            {
                double sign = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0);
                v[i].MakeZero();
                v[i][perm[p][i]] = sign;
            }
            v[2] = Cross(v[0], v[1]);
        }
        else
        {
            double len = 0.0;
            do
            {
                for (int i = 0; i < 3; ++i) { v[0][i] = io.raw(-1.0, 1.0); }
                len = Length(v[0]);
            } while (len < 0.1 || len > 1.0);
            Normalize(v[0]);
            ComputeOrthogonalComplement(1, v.data());
        }
        return v;
    }

    // Records 4 doubles.
    std::array<Vector2<double>, 2> Frame2(oracle::Ctx& io, int mode, bool leftHanded)
    {
        std::array<Vector2<double>, 2> v = DrawFrame2(io, mode);
        if (leftHanded) { v[1] = -v[1]; }
        std::array<Vector2<double>, 2> axis{};
        axis[0] = io.givenVec(v[0]);
        axis[1] = io.givenVec(v[1]);
        return axis;
    }

    // Records 9 doubles.
    std::array<Vector3<double>, 3> Frame3(oracle::Ctx& io, int mode, bool leftHanded)
    {
        std::array<Vector3<double>, 3> v = DrawFrame3(io, mode);
        if (leftHanded) { v[2] = -v[2]; }
        std::array<Vector3<double>, 3> axis{};
        axis[0] = io.givenVec(v[0]);
        axis[1] = io.givenVec(v[1]);
        axis[2] = io.givenVec(v[2]);
        return axis;
    }

    // An unrecorded point: lattice in mode 0, uniform otherwise.
    template <int N>
    Vector<N, double> DrawPt(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> p{};
        for (int k = 0; k < N; ++k)
        {
            p[k] = (mode == 0 ? static_cast<double>(io.rawInteger(-lat, lat))
                : io.raw(-range, range));
        }
        return p;
    }

    // A signed coordinate axis (an exact unit vector). Records N doubles.
    template <int N>
    Vector<N, double> AxisUnit(oracle::Ctx& io)
    {
        int k = io.rawInteger(0, 2 * N - 1);
        Vector<N, double> v{};
        v.MakeZero();
        v[k % N] = (k < N ? 1.0 : -1.0);
        return io.givenVec(v);
    }

    // A unit-length direction: an exact coordinate axis in lattice mode, a
    // uniformly distributed unit vector otherwise. Records N doubles.
    template <int N>
    Vector<N, double> Dir(oracle::Ctx& io, int mode)
    {
        return (mode == 0 ? AxisUnit<N>(io) : io.unit<N>());
    }

    // Positive extents; records N doubles.
    template <int N>
    Vector<N, double> Extent(oracle::Ctx& io, int mode, int lat, double range)
    {
        Vector<N, double> e{};
        for (int i = 0; i < N; ++i)
        {
            e[i] = (mode == 0 ? static_cast<double>(io.integer(1, lat))
                : io.real(0.25, range));
        }
        return e;
    }

    // Records the box as min then max (2*N doubles) in every mode, so that a
    // constructed box and a drawn box have the same layout.
    template <int N>
    AlignedBox<N, double> ABox(oracle::Ctx& io, int mode, int lat, double range)
    {
        AlignedBox<N, double> box{};
        for (int i = 0; i < N; ++i)
        {
            double center = (mode == 0 ? static_cast<double>(io.rawInteger(-lat, lat))
                : io.raw(-range, range));
            double extent = (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
                : io.raw(0.25, range));
            box.min[i] = center - extent;
            box.max[i] = center + extent;
        }
        io.givenVec(box.min);
        io.givenVec(box.max);
        return box;
    }

    template <int N>
    OrientedBox<N, double> OBox(oracle::Ctx& io, int mode, int lat, double range,
        bool leftHanded)
    {
        OrientedBox<N, double> box{};
        box.center = Pt<N>(io, mode, lat, range);
        if constexpr (N == 2) { box.axis = Frame2(io, mode, leftHanded); }
        else { box.axis = Frame3(io, mode, leftHanded); }
        box.extent = Extent<N>(io, mode, lat, range);
        return box;
    }

    template <int N>
    Hyperellipsoid<N, double> Ellip(oracle::Ctx& io, int mode, int lat, double range,
        bool leftHanded)
    {
        Hyperellipsoid<N, double> e{};
        e.center = Pt<N>(io, mode, lat, range);
        if constexpr (N == 2) { e.axis = Frame2(io, mode, leftHanded); }
        else { e.axis = Frame3(io, mode, leftHanded); }
        e.extent = Extent<N>(io, mode, lat, range);
        return e;
    }

    template <int N>
    Hypersphere<N, double> Sph(oracle::Ctx& io, int mode, int lat, double range)
    {
        Hypersphere<N, double> s{};
        s.center = Pt<N>(io, mode, lat, range);
        s.radius = (mode == 0 ? io.lattice(0, lat) : io.real(0.25, range));
        return s;
    }

    Capsule3<double> Cap3(oracle::Ctx& io, int mode, int lat, double range)
    {
        Capsule3<double> c{};
        Vector3<double> p0 = Pt<3>(io, mode, lat, range);
        Vector3<double> p1 = Pt<3>(io, mode, lat, range);
        c.segment = Segment3<double>(p0, p1);
        c.radius = (mode == 0 ? io.lattice(1, lat) : io.real(0.25, range));
        return c;
    }
}

namespace
{
    // A point placed relative to an aligned box so that InContainer's '<' and
    // '>' are evaluated at exact equality: 'choice' selects, per coordinate,
    // box.min[i], box.max[i], the midpoint or a point outside. Records N
    // doubles.
    template <int N>
    Vector<N, double> BoxProbe(oracle::Ctx& io, AlignedBox<N, double> const& box)
    {
        Vector<N, double> p{};
        for (int i = 0; i < N; ++i)
        {
            switch (io.rawInteger(0, 4))
            {
            case 0: p[i] = box.min[i]; break;
            case 1: p[i] = box.max[i]; break;
            case 2: p[i] = 0.5 * (box.min[i] + box.max[i]); break;
            case 3: p[i] = box.min[i] - 1.0; break;
            default: p[i] = box.max[i] + 1.0; break;
            }
        }
        return io.givenVec(p);
    }

    template <int N>
    void OutABox(oracle::Ctx& io, AlignedBox<N, double> const& box)
    {
        io.outVec(box.min);
        io.outVec(box.max);
    }
}

// ============================== ContAlignedBox ==========================

ORACLE_CASE("ContAlignedBox.getContainer.2d")
{
    int kind = io.index() % 5;
    int count = io.integer(1, 8);
    auto points = Cloud<2>(io, kind, count, 5, 8.0);
    AlignedBox2<double> box{};
    bool ok = GetContainer(count, points.data(), box);
    io.outBool(ok);
    OutABox<2>(io, box);
}

ORACLE_CASE("ContAlignedBox.getContainer.3d")
{
    int kind = io.index() % 5;
    int count = io.integer(1, 8);
    auto points = Cloud<3>(io, kind, count, 5, 8.0);
    AlignedBox3<double> box{};
    bool ok = GetContainer(count, points.data(), box);
    io.outBool(ok);
    OutABox<3>(io, box);
}

ORACLE_CASE("ContAlignedBox.inContainer.2d")
{
    int mode = io.index() % 2;
    auto box = ABox<2>(io, mode, 4, 5.0);
    auto p = BoxProbe<2>(io, box);
    io.outBool(InContainer(p, box));
}

ORACLE_CASE("ContAlignedBox.inContainer.3d")
{
    int mode = io.index() % 2;
    auto box = ABox<3>(io, mode, 4, 5.0);
    auto p = BoxProbe<3>(io, box);
    io.outBool(InContainer(p, box));
}

ORACLE_CASE("ContAlignedBox.mergeContainers.3d")
{
    // Every fourth record merges a box with itself (identical inputs) and
    // every fourth shrinks the second box into the first (nested).
    int mode = io.index() % 2;
    int rel = (io.index() / 2) % 4;
    auto box0 = ABox<3>(io, mode, 4, 5.0);
    AlignedBox3<double> box1{};
    if (rel == 0)
    {
        box1 = box0;
        io.givenVec(box1.min);
        io.givenVec(box1.max);
    }
    else if (rel == 1)
    {
        for (int i = 0; i < 3; ++i)
        {
            double c = 0.5 * (box0.min[i] + box0.max[i]);
            box1.min[i] = 0.5 * (box0.min[i] + c);
            box1.max[i] = 0.5 * (box0.max[i] + c);
        }
        io.givenVec(box1.min);
        io.givenVec(box1.max);
    }
    else
    {
        box1 = ABox<3>(io, mode, 4, 5.0);
    }
    AlignedBox3<double> merge{};
    bool ok = MergeContainers(box0, box1, merge);
    io.outBool(ok);
    OutABox<3>(io, merge);
}

ORACLE_CASE("ContAlignedBox.mergeContainers.signedZero")
{
    // std::min(a,b) is (b < a ? b : a) and std::max(a,b) is (a < b ? b : a),
    // so both return their FIRST argument when the arguments are +0 and -0.
    // Math.min / Math.max order -0 below +0, so this case pins the port's
    // stdMin / stdMax.
    AlignedBox3<double> box0{}, box1{};
    for (int i = 0; i < 3; ++i)
    {
        int s0 = io.integer(0, 1);
        int s1 = io.integer(0, 1);
        box0.min[i] = (s0 == 0 ? 0.0 : -0.0);
        box0.max[i] = (s1 == 0 ? 0.0 : -0.0);
    }
    for (int i = 0; i < 3; ++i)
    {
        int s0 = io.integer(0, 1);
        int s1 = io.integer(0, 1);
        box1.min[i] = (s0 == 0 ? 0.0 : -0.0);
        box1.max[i] = (s1 == 0 ? 0.0 : -0.0);
    }
    AlignedBox3<double> merge{};
    bool ok = MergeContainers(box0, box1, merge);
    io.outBool(ok);
    OutABox<3>(io, merge);
}

// ======================= ContCircle2 / ContSphere3 ======================

namespace
{
    // A point placed relative to a hypersphere so that 'Length(diff) <= r' is
    // evaluated at exact equality: the probe is center + r*axis for a signed
    // coordinate axis (exact), the center itself, or a point outside.
    // Records N doubles.
    template <int N>
    Vector<N, double> SphereProbe(oracle::Ctx& io, Hypersphere<N, double> const& s)
    {
        Vector<N, double> p = s.center;
        int choice = io.rawInteger(0, 3);
        if (choice != 0)
        {
            int k = io.rawInteger(0, 2 * N - 1);
            double sign = (k < N ? 1.0 : -1.0);
            double t = (choice == 1 ? s.radius
                : (choice == 2 ? 0.5 * s.radius : s.radius + 1.0));
            p[k % N] += sign * t;
        }
        return io.givenVec(p);
    }

    template <int N>
    void OutSphere(oracle::Ctx& io, Hypersphere<N, double> const& s)
    {
        io.outVec(s.center);
        io.outReal(s.radius);
    }

    // Two spheres in a chosen relationship: 0 identical, 1 the second nested
    // inside the first, 2 the first nested inside the second, 3 independent.
    // Records 2*(N+1) doubles.
    template <int N>
    void SpherePair(oracle::Ctx& io, int mode, int rel, int lat, double range,
        Hypersphere<N, double>& s0, Hypersphere<N, double>& s1)
    {
        s0 = Sph<N>(io, mode, lat, range);
        if (rel == 3)
        {
            s1 = Sph<N>(io, mode, lat, range);
            return;
        }
        if (rel == 0)
        {
            // Identical: lenSqr = 0 and rDiff = 0, so rDiffSqr >= lenSqr.
            s1 = s0;
        }
        else if (rel == 1)
        {
            // s1 strictly inside s0, off-centre.
            s1.center = s0.center;
            s1.center[io.rawInteger(0, N - 1)] += 0.25 * s0.radius;
            s1.radius = 0.5 * s0.radius;
        }
        else
        {
            // s0 strictly inside s1.
            s1.center = s0.center;
            s1.radius = 2.0 * s0.radius + 1.0;
        }
        io.givenVec(s1.center);
        io.given(s1.radius);
    }
}

ORACLE_CASE("ContCircle2.getContainer")
{
    int kind = io.index() % 5;
    int count = io.integer(1, 8);
    auto points = Cloud<2>(io, kind, count, 5, 8.0);
    Circle2<double> circle{};
    bool ok = GetContainer(count, points.data(), circle);
    io.outBool(ok);
    OutSphere<2>(io, circle);
}

ORACLE_CASE("ContCircle2.inContainer")
{
    int mode = io.index() % 2;
    auto circle = Sph<2>(io, mode, 4, 5.0);
    auto p = SphereProbe<2>(io, circle);
    io.outBool(InContainer(p, circle));
}

ORACLE_CASE("ContCircle2.mergeContainers")
{
    int mode = io.index() % 2;
    int rel = (io.index() / 2) % 4;
    Circle2<double> c0{}, c1{}, merge{};
    SpherePair<2>(io, mode, rel, 4, 5.0, c0, c1);
    bool ok = MergeContainers(c0, c1, merge);
    io.outBool(ok);
    OutSphere<2>(io, merge);
}

ORACLE_CASE("ContSphere3.getContainer")
{
    int kind = io.index() % 5;
    int count = io.integer(1, 8);
    auto points = Cloud<3>(io, kind, count, 5, 8.0);
    Sphere3<double> sphere{};
    bool ok = GetContainer(count, points.data(), sphere);
    io.outBool(ok);
    OutSphere<3>(io, sphere);
}

ORACLE_CASE("ContSphere3.inContainer")
{
    int mode = io.index() % 2;
    auto sphere = Sph<3>(io, mode, 4, 5.0);
    auto p = SphereProbe<3>(io, sphere);
    io.outBool(InContainer(p, sphere));
}

ORACLE_CASE("ContSphere3.mergeContainers")
{
    int mode = io.index() % 2;
    int rel = (io.index() / 2) % 4;
    Sphere3<double> s0{}, s1{}, merge{};
    SpherePair<3>(io, mode, rel, 4, 5.0, s0, s1);
    bool ok = MergeContainers(s0, s1, merge);
    io.outBool(ok);
    OutSphere<3>(io, merge);
}

// =========================== ContPointInPolygon2 ========================

namespace
{
    // Edge directions with strictly increasing polar angle in [0, pi).
    double const kEdgeDir[12][2] =
    {
        {  1.0,  0.0 }, {  3.0,  1.0 }, {  2.0,  1.0 }, {  1.0,  1.0 },
        {  1.0,  2.0 }, {  1.0,  3.0 }, {  0.0,  1.0 }, { -1.0,  3.0 },
        { -1.0,  2.0 }, { -1.0,  1.0 }, { -2.0,  1.0 }, { -3.0,  1.0 }
    };

    // A strictly convex, counterclockwise polygon with 2*m vertices: m edge
    // vectors of strictly increasing angle followed by their negations, which
    // closes the polygon exactly. Mode 0 scales each direction by a positive
    // integer (an exact lattice polygon); mode 1 scales by a positive real.
    // Records 4*m doubles.
    std::vector<Vector2<double>> ConvexPoly(oracle::Ctx& io, int m, int mode, int lat)
    {
        int maxStep = 11 / (m - 1);
        int step = io.rawInteger(1, maxStep);
        int start = io.rawInteger(0, 11 - (m - 1) * step);

        std::vector<Vector2<double>> edge(static_cast<size_t>(m));
        for (int k = 0; k < m; ++k)
        {
            int idx = start + k * step;
            double scale = (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
                : io.raw(0.5, 3.0));
            edge[static_cast<size_t>(k)][0] = scale * kEdgeDir[idx][0];
            edge[static_cast<size_t>(k)][1] = scale * kEdgeDir[idx][1];
        }

        Vector2<double> v{};
        v[0] = (mode == 0 ? static_cast<double>(io.rawInteger(-lat, lat))
            : io.raw(-4.0, 4.0));
        v[1] = (mode == 0 ? static_cast<double>(io.rawInteger(-lat, lat))
            : io.raw(-4.0, 4.0));

        std::vector<Vector2<double>> poly;
        poly.reserve(static_cast<size_t>(2 * m));
        for (int k = 0; k < m; ++k)
        {
            poly.push_back(io.givenVec(v));
            v = v + edge[static_cast<size_t>(k)];
        }
        for (int k = 0; k < m; ++k)
        {
            poly.push_back(io.givenVec(v));
            v = v - edge[static_cast<size_t>(k)];
        }
        return poly;
    }

    // A test point for a polygon query. The constructed choices put the point
    // exactly on a vertex, exactly on an edge midpoint or quarter point, and
    // exactly on the horizontal ray through a vertex, which is where the
    // '<' / '<=' asymmetry of Contains is decided. Records 2 doubles.
    Vector2<double> PolyProbe(oracle::Ctx& io, std::vector<Vector2<double>> const& poly,
        int lat, double range)
    {
        int n = static_cast<int>(poly.size());
        Vector2<double> p{};
        switch (io.rawInteger(0, 5))
        {
        case 0:
            p = poly[static_cast<size_t>(io.rawInteger(0, n - 1))];
            break;
        case 1:
        {
            int i = io.rawInteger(0, n - 1);
            int j = (i + 1) % n;
            p = 0.5 * (poly[static_cast<size_t>(i)] + poly[static_cast<size_t>(j)]);
            break;
        }
        case 2:
        {
            int i = io.rawInteger(0, n - 1);
            int j = (i + 1) % n;
            p = poly[static_cast<size_t>(i)]
                + 0.25 * (poly[static_cast<size_t>(j)] - poly[static_cast<size_t>(i)]);
            break;
        }
        case 3:
        {
            // Same y as a vertex, arbitrary x: the horizontal ray passes
            // exactly through a vertex.
            int i = io.rawInteger(0, n - 1);
            p[0] = static_cast<double>(io.rawInteger(-2 * lat, 2 * lat));
            p[1] = poly[static_cast<size_t>(i)][1];
            break;
        }
        case 4:
            p[0] = static_cast<double>(io.rawInteger(-2 * lat, 2 * lat));
            p[1] = static_cast<double>(io.rawInteger(-2 * lat, 2 * lat));
            break;
        default:
            p[0] = io.raw(-range, range);
            p[1] = io.raw(-range, range);
            break;
        }
        return io.givenVec(p);
    }
}

ORACLE_CASE("ContPointInPolygon2.contains.convex")
{
    int mode = io.index() % 2;
    int m = 2 + (io.index() / 2) % 4;
    io.integer(m, m);
    auto poly = ConvexPoly(io, m, mode, 4);
    auto p = PolyProbe(io, poly, 4, 6.0);
    PointInPolygon2<double> pip(static_cast<int32_t>(poly.size()), poly.data());
    io.outBool(pip.Contains(p));
}

ORACLE_CASE("ContPointInPolygon2.contains.nonconvex")
{
    // An L-shaped hexagon (counterclockwise), which the ray-counting query
    // must handle and the convex queries must not be given.
    int mode = io.index() % 2;
    double a = (mode == 0 ? static_cast<double>(io.rawInteger(2, 6)) : io.raw(1.0, 6.0));
    double b = (mode == 0 ? static_cast<double>(io.rawInteger(2, 6)) : io.raw(1.0, 6.0));
    double ox = (mode == 0 ? static_cast<double>(io.rawInteger(-4, 4)) : io.raw(-4.0, 4.0));
    double oy = (mode == 0 ? static_cast<double>(io.rawInteger(-4, 4)) : io.raw(-4.0, 4.0));
    std::vector<Vector2<double>> poly(6);
    poly[0] = { ox, oy };
    poly[1] = { ox + 2.0 * a, oy };
    poly[2] = { ox + 2.0 * a, oy + b };
    poly[3] = { ox + a, oy + b };
    poly[4] = { ox + a, oy + 2.0 * b };
    poly[5] = { ox, oy + 2.0 * b };
    for (auto const& q : poly) { io.givenVec(q); }
    auto p = PolyProbe(io, poly, 6, 10.0);
    PointInPolygon2<double> pip(6, poly.data());
    io.outBool(pip.Contains(p));
}

ORACLE_CASE("ContPointInPolygon2.containsConvexOrderN")
{
    int mode = io.index() % 2;
    int m = 2 + (io.index() / 2) % 4;
    io.integer(m, m);
    auto poly = ConvexPoly(io, m, mode, 4);
    auto p = PolyProbe(io, poly, 4, 6.0);
    PointInPolygon2<double> pip(static_cast<int32_t>(poly.size()), poly.data());
    io.outBool(pip.ContainsConvexOrderN(p));
}

ORACLE_CASE("ContPointInPolygon2.containsConvexOrderLogN")
{
    int mode = io.index() % 2;
    int m = 2 + (io.index() / 2) % 4;
    io.integer(m, m);
    auto poly = ConvexPoly(io, m, mode, 4);
    auto p = PolyProbe(io, poly, 4, 6.0);
    PointInPolygon2<double> pip(static_cast<int32_t>(poly.size()), poly.data());
    io.outBool(pip.ContainsConvexOrderLogN(p));
}

ORACLE_CASE("ContPointInPolygon2.containsQuadrilateral")
{
    // Every third record uses a hexagon, for which upstream returns false
    // because the vertex count is not four.
    int mode = io.index() % 2;
    int m = (io.index() % 3 == 2 ? 3 : 2);
    io.integer(m, m);
    auto poly = ConvexPoly(io, m, mode, 4);
    auto p = PolyProbe(io, poly, 4, 6.0);
    PointInPolygon2<double> pip(static_cast<int32_t>(poly.size()), poly.data());
    io.outBool(pip.ContainsQuadrilateral(p));
}

// =========================== ContAlignedBox2Arc2 ========================

ORACLE_CASE("ContAlignedBox2Arc2.getContainer")
{
    // Mode 0 places the arc endpoints on exact lattice points of a radius-5
    // circle, so the four axis points C+(+-r,0), C+(0,+-r) are themselves
    // possible endpoints and Arc2::Contains is evaluated at DotPerp == 0.
    // Records center (2), radius (1), end[0] (2), end[1] (2).
    int mode = io.index() % 2;
    Arc2<double> arc{};
    if (mode == 0)
    {
        double const pt[12][2] =
        {
            {  5.0,  0.0 }, {  4.0,  3.0 }, {  3.0,  4.0 }, {  0.0,  5.0 },
            { -3.0,  4.0 }, { -4.0,  3.0 }, { -5.0,  0.0 }, { -4.0, -3.0 },
            { -3.0, -4.0 }, {  0.0, -5.0 }, {  3.0, -4.0 }, {  4.0, -3.0 }
        };
        Vector2<double> center{ static_cast<double>(io.rawInteger(-4, 4)),
            static_cast<double>(io.rawInteger(-4, 4)) };
        int i0 = io.rawInteger(0, 11);
        int i1 = io.rawInteger(0, 11);
        arc.center = io.givenVec(center);
        arc.radius = io.given(5.0);
        arc.end[0] = io.givenVec(Vector2<double>{ center[0] + pt[i0][0],
            center[1] + pt[i0][1] });
        arc.end[1] = io.givenVec(Vector2<double>{ center[0] + pt[i1][0],
            center[1] + pt[i1][1] });
    }
    else
    {
        Vector2<double> center{ io.raw(-4.0, 4.0), io.raw(-4.0, 4.0) };
        double radius = io.raw(0.5, 4.0);
        double t0 = io.raw(0.0, 6.2831853071795862);
        double t1 = io.raw(0.0, 6.2831853071795862);
        arc.center = io.givenVec(center);
        arc.radius = io.given(radius);
        arc.end[0] = io.givenVec(Vector2<double>{
            center[0] + radius * std::cos(t0), center[1] + radius * std::sin(t0) });
        arc.end[1] = io.givenVec(Vector2<double>{
            center[0] + radius * std::cos(t1), center[1] + radius * std::sin(t1) });
    }
    AlignedBox2<double> box{};
    bool ok = GetContainer(arc, box);
    io.outBool(ok);
    OutABox<2>(io, box);
}

// ============================== ContCapsule3 ============================

namespace
{
    void OutCapsule(oracle::Ctx& io, Capsule3<double> const& c)
    {
        io.outVec(c.segment.p[0]);
        io.outVec(c.segment.p[1]);
        io.outReal(c.radius);
    }
}

ORACLE_CASE("ContCapsule3.getContainer")
{
    int kind = io.index() % 5;
    int count = io.integer(2, 8);
    auto points = Cloud<3>(io, kind, count, 5, 8.0);
    Capsule3<double> capsule{};
    bool ok = GetContainer(count, points.data(), capsule);
    io.outBool(ok);
    OutCapsule(io, capsule);
}

ORACLE_CASE("ContCapsule3.inContainer.point")
{
    // Every second record places the point exactly at a segment endpoint
    // offset by the radius along a coordinate axis, where 'distance <= radius'
    // is evaluated at equality.
    int mode = io.index() % 2;
    auto capsule = Cap3(io, mode, 4, 5.0);
    Vector3<double> p{};
    int choice = io.rawInteger(0, 3);
    if (choice == 0)
    {
        int e = io.rawInteger(0, 1);
        int k = io.rawInteger(0, 5);
        p = capsule.segment.p[e];
        p[k % 3] += (k < 3 ? capsule.radius : -capsule.radius);
    }
    else if (choice == 1)
    {
        p = 0.5 * (capsule.segment.p[0] + capsule.segment.p[1]);
    }
    else
    {
        for (int k = 0; k < 3; ++k)
        {
            p[k] = (mode == 0 ? static_cast<double>(io.rawInteger(-6, 6))
                : io.raw(-7.0, 7.0));
        }
    }
    io.givenVec(p);
    io.outBool(InContainer(p, capsule));
}

ORACLE_CASE("ContCapsule3.inContainer.sphere")
{
    int mode = io.index() % 2;
    auto capsule = Cap3(io, mode, 4, 5.0);
    auto sphere = Sph<3>(io, mode, 4, 5.0);
    io.outBool(InContainer<double>(sphere, capsule));
}

ORACLE_CASE("ContCapsule3.inContainer.capsule")
{
    int mode = io.index() % 2;
    auto testCapsule = Cap3(io, mode, 4, 5.0);
    auto capsule = Cap3(io, mode, 4, 5.0);
    io.outBool(InContainer<double>(testCapsule, capsule));
}

ORACLE_CASE("ContCapsule3.mergeContainers")
{
    // rel 0 merges a capsule with itself (both containment tests succeed),
    // rel 1 nests the second capsule inside the first, rel 2 and 3 are
    // independent capsules, which is the general branch.
    int mode = io.index() % 2;
    int rel = (io.index() / 2) % 4;
    auto c0 = Cap3(io, mode, 4, 5.0);
    Capsule3<double> c1{};
    if (rel == 0)
    {
        c1 = c0;
    }
    else if (rel == 1)
    {
        Vector3<double> center = 0.5 * (c0.segment.p[0] + c0.segment.p[1]);
        c1.segment.p[0] = 0.5 * (c0.segment.p[0] + center);
        c1.segment.p[1] = 0.5 * (c0.segment.p[1] + center);
        c1.radius = 0.25 * c0.radius;
    }
    else
    {
        c1 = Cap3(io, mode, 4, 5.0);
    }
    if (rel < 2)
    {
        io.givenVec(c1.segment.p[0]);
        io.givenVec(c1.segment.p[1]);
        io.given(c1.radius);
    }
    Capsule3<double> merge{};
    bool ok = MergeContainers(c0, c1, merge);
    io.outBool(ok);
    OutCapsule(io, merge);
}

// ============================== ContCylinder3 ===========================

ORACLE_CASE("ContCylinder3.getContainer")
{
    int kind = io.index() % 5;
    int count = io.integer(2, 8);
    auto points = Cloud<3>(io, kind, count, 5, 8.0);
    Cylinder3<double> cylinder{};
    bool ok = GetContainer(count, points.data(), cylinder);
    io.outBool(ok);
    io.outVec(cylinder.axis.origin);
    io.outVec(cylinder.axis.direction);
    io.outReal(cylinder.radius);
    io.outReal(cylinder.height);
}

ORACLE_CASE("ContCylinder3.inContainer")
{
    // Mode 0 gives an exact coordinate-axis cylinder, and the probe places the
    // point exactly on the flat cap (|zProj|*2 == height) or exactly on the
    // lateral surface (Length(xyProj) == radius).
    int mode = io.index() % 2;
    Cylinder3<double> cylinder{};
    cylinder.axis.origin = Pt<3>(io, mode, 4, 5.0);
    cylinder.axis.direction = Dir<3>(io, mode);
    cylinder.radius = (mode == 0 ? io.lattice(1, 4) : io.real(0.25, 4.0));
    cylinder.height = (mode == 0 ? io.lattice(1, 6) : io.real(0.25, 6.0));

    Vector3<double> basis[3];
    basis[0] = cylinder.axis.direction;
    ComputeOrthogonalComplement(1, basis);
    Vector3<double> p{};
    switch (io.rawInteger(0, 3))
    {
    case 0:
        p = cylinder.axis.origin + (0.5 * cylinder.height) * cylinder.axis.direction;
        break;
    case 1:
        p = cylinder.axis.origin + cylinder.radius * basis[1];
        break;
    case 2:
        p = cylinder.axis.origin - (0.5 * cylinder.height) * cylinder.axis.direction
            + cylinder.radius * basis[2];
        break;
    default:
        for (int k = 0; k < 3; ++k)
        {
            p[k] = (mode == 0 ? static_cast<double>(io.rawInteger(-6, 6))
                : io.raw(-7.0, 7.0));
        }
        break;
    }
    io.givenVec(p);
    io.outBool(InContainer(p, cylinder));
}

// =============================== ContLozenge3 ===========================

namespace
{
    // The fields of a lozenge other than the rectangle centre, which is the
    // one the port deliberately computes differently (issue #174).
    void OutLozengeNoCenter(oracle::Ctx& io, Lozenge3<double> const& l)
    {
        io.outReal(l.radius);
        io.outVec(l.rectangle.axis[0]);
        io.outVec(l.rectangle.axis[1]);
        io.outReal(l.rectangle.extent[0]);
        io.outReal(l.rectangle.extent[1]);
    }
}

ORACLE_CASE("ContLozenge3.getContainer")
{
    // The rectangle CENTRE is not emitted here: upstream centres the
    // rectangle on a corner of the fitted parameter interval while the port
    // uses the midpoint (docs/UPSTREAM-FINDINGS.md ContLozenge3.h, issue
    // #174). Everything else agrees bit for bit. The centre is compared in
    // ContLozenge3.getContainer.sphereCenter (where the two formulas coincide
    // mathematically) and its deviation is demonstrated in
    // ContLozenge3.getContainer.cornerDeviation.
    int kind = io.index() % 5;
    int count = io.integer(2, 8);
    auto points = Cloud<3>(io, kind, count, 5, 8.0);
    Lozenge3<double> lozenge{};
    bool ok = GetContainer(count, points.data(), lozenge);
    io.outBool(ok);
    OutLozengeNoCenter(io, lozenge);
}

ORACLE_CASE("ContLozenge3.getContainer.cornerDeviation")
{
    // Same generator, emitting only the rectangle centre. Upstream's
    //   center = box.center + aMin*axis[2] + bMin*axis[1]
    // is off by the rectangle's own half-extents; the port uses the interval
    // midpoints in every branch (issue #174).
    int kind = io.index() % 2;
    int count = io.integer(3, 8);
    auto points = Cloud<3>(io, kind, count, 5, 8.0);
    Lozenge3<double> lozenge{};
    GetContainer(count, points.data(), lozenge);
    io.outVec(lozenge.rectangle.center);
}

ORACLE_CASE("ContLozenge3.getContainer.sphereBranch")
{
    // Coincident points are the only configuration that reaches upstream's
    // "container is a sphere" branch, where its own formula is
    //   box.center + (0.5*(aMin+aMax))*axis[2] + (0.5*(bMin+bMax))*axis[1],
    // which is the expression the port uses in every branch. (The extreme-w
    // points have radical == 0, so aMin >= aMax forces them to share a single
    // u, and because box.center is the mean that shared value is 0; a
    // collinear cloud therefore lands in the capsule branch, not this one.)
    // Both offsets are zero here, so the case pins the branch and the centre
    // but not the accumulation order of the sum; that is covered by the
    // regression test in test/ContLozenge3.test.ts.
    int count = io.integer(3, 8);
    auto points = Cloud<3>(io, 2, count, 5, 8.0);
    Lozenge3<double> lozenge{};
    bool ok = GetContainer(count, points.data(), lozenge);
    io.outBool(ok);
    io.outVec(lozenge.rectangle.center);
    OutLozengeNoCenter(io, lozenge);
}

ORACLE_CASE("ContLozenge3.inContainer")
{
    int mode = io.index() % 2;
    Lozenge3<double> lozenge{};
    lozenge.rectangle.center = Pt<3>(io, mode, 4, 5.0);
    auto frame = Frame3(io, mode, false);
    lozenge.rectangle.axis[0] = frame[0];
    lozenge.rectangle.axis[1] = frame[1];
    lozenge.rectangle.extent = Extent<2>(io, mode, 4, 4.0);
    lozenge.radius = (mode == 0 ? io.lattice(0, 3) : io.real(0.0, 3.0));

    Vector3<double> p{};
    switch (io.rawInteger(0, 3))
    {
    case 0:
        // Exactly on the rectangle corner offset by the radius along the
        // rectangle normal, where 'distance <= radius' is an equality.
        p = lozenge.rectangle.center
            + lozenge.rectangle.extent[0] * lozenge.rectangle.axis[0]
            + lozenge.rectangle.extent[1] * lozenge.rectangle.axis[1]
            + lozenge.radius * frame[2];
        break;
    case 1:
        p = lozenge.rectangle.center + lozenge.radius * frame[2];
        break;
    default:
        for (int k = 0; k < 3; ++k)
        {
            p[k] = (mode == 0 ? static_cast<double>(io.rawInteger(-7, 7))
                : io.raw(-8.0, 8.0));
        }
        break;
    }
    io.givenVec(p);
    io.outBool(InContainer(p, lozenge));
}

// ====================== ContOrientedBox2 / ContOrientedBox3 =============

namespace
{
    template <int N>
    void OutOBox(oracle::Ctx& io, OrientedBox<N, double> const& box)
    {
        io.outVec(box.center);
        for (int i = 0; i < N; ++i) { io.outVec(box.axis[i]); }
        io.outVec(box.extent);
    }

    // A point placed relative to an oriented box or ellipsoid frame so that
    // the '|coeff| > extent' / 'Length(standardized) <= 1' tests are evaluated
    // at exact equality. Records N doubles.
    template <int N>
    Vector<N, double> FrameProbe(oracle::Ctx& io, Vector<N, double> const& center,
        std::array<Vector<N, double>, N> const& axis, Vector<N, double> const& extent,
        int mode)
    {
        Vector<N, double> p = center;
        int choice = io.rawInteger(0, 3);
        if (choice == 0)
        {
            // A face centre: |coeff| == extent exactly along one axis.
            int k = io.rawInteger(0, N - 1);
            int s = io.rawInteger(0, 1);
            p = center + (s == 0 ? extent[k] : -extent[k]) * axis[k];
        }
        else if (choice == 1)
        {
            // A corner.
            for (int k = 0; k < N; ++k)
            {
                p = p + (io.rawInteger(0, 1) == 0 ? extent[k] : -extent[k]) * axis[k];
            }
        }
        else
        {
            for (int k = 0; k < N; ++k)
            {
                p[k] = (mode == 0 ? static_cast<double>(io.rawInteger(-8, 8))
                    : io.raw(-9.0, 9.0));
            }
        }
        return io.givenVec(p);
    }
}

ORACLE_CASE("ContOrientedBox2.getContainer")
{
    int kind = io.index() % 5;
    int count = io.integer(1, 8);
    auto points = Cloud<2>(io, kind, count, 5, 8.0);
    OrientedBox2<double> box{};
    bool ok = GetContainer(count, points.data(), box);
    io.outBool(ok);
    OutOBox<2>(io, box);
}

ORACLE_CASE("ContOrientedBox2.inContainer")
{
    int mode = io.index() % 2;
    auto box = OBox<2>(io, mode, 4, 5.0, false);
    auto p = FrameProbe<2>(io, box.center, box.axis, box.extent, mode);
    io.outBool(InContainer(p, box));
}

ORACLE_CASE("ContOrientedBox2.inContainer.leftHanded")
{
    // det(axis) = -1. InContainer projects onto each axis independently, so
    // it must not depend on the handedness of the frame; this case records
    // what upstream actually does with one.
    int mode = io.index() % 2;
    auto box = OBox<2>(io, mode, 4, 5.0, true);
    auto p = FrameProbe<2>(io, box.center, box.axis, box.extent, mode);
    io.outBool(InContainer(p, box));
}

ORACLE_CASE("ContOrientedBox2.mergeContainers")
{
    int mode = io.index() % 2;
    int rel = (io.index() / 2) % 3;
    auto box0 = OBox<2>(io, mode, 4, 5.0, false);
    OrientedBox2<double> box1{};
    if (rel == 0)
    {
        box1 = box0;
        io.givenVec(box1.center);
        io.givenVec(box1.axis[0]);
        io.givenVec(box1.axis[1]);
        io.givenVec(box1.extent);
    }
    else
    {
        box1 = OBox<2>(io, mode, 4, 5.0, false);
    }
    OrientedBox2<double> merge{};
    bool ok = MergeContainers(box0, box1, merge);
    io.outBool(ok);
    OutOBox<2>(io, merge);
}

ORACLE_CASE("ContOrientedBox3.getContainer")
{
    int kind = io.index() % 5;
    int count = io.integer(1, 8);
    auto points = Cloud<3>(io, kind, count, 5, 8.0);
    OrientedBox3<double> box{};
    bool ok = GetContainer(count, points.data(), box);
    io.outBool(ok);
    OutOBox<3>(io, box);
}

ORACLE_CASE("ContOrientedBox3.inContainer")
{
    int mode = io.index() % 2;
    auto box = OBox<3>(io, mode, 4, 5.0, false);
    auto p = FrameProbe<3>(io, box.center, box.axis, box.extent, mode);
    io.outBool(InContainer(p, box));
}

ORACLE_CASE("ContOrientedBox3.inContainer.leftHanded")
{
    int mode = io.index() % 2;
    auto box = OBox<3>(io, mode, 4, 5.0, true);
    auto p = FrameProbe<3>(io, box.center, box.axis, box.extent, mode);
    io.outBool(InContainer(p, box));
}

ORACLE_CASE("ContOrientedBox3.mergeContainers")
{
    int mode = io.index() % 2;
    int rel = (io.index() / 2) % 3;
    auto box0 = OBox<3>(io, mode, 4, 5.0, false);
    OrientedBox3<double> box1{};
    if (rel == 0)
    {
        box1 = box0;
        io.givenVec(box1.center);
        io.givenVec(box1.axis[0]);
        io.givenVec(box1.axis[1]);
        io.givenVec(box1.axis[2]);
        io.givenVec(box1.extent);
    }
    else
    {
        box1 = OBox<3>(io, mode, 4, 5.0, false);
    }
    OrientedBox3<double> merge{};
    bool ok = MergeContainers(box0, box1, merge);
    io.outBool(ok);
    OutOBox<3>(io, merge);
}

ORACLE_CASE("ContOrientedBox3.mergeContainers.leftHanded")
{
    // Upstream's merge converts each axis triple to a quaternion with
    // Rotation<3,Real>, which assumes det = +1 (docs/UPSTREAM-FINDINGS.md
    // ContOrientedBox3.h, issue #292). The defect is preserved in the port, so
    // the two sides are still compared bit for bit here.
    int mode = io.index() % 2;
    auto box0 = OBox<3>(io, mode, 4, 5.0, true);
    auto box1 = OBox<3>(io, mode, 4, 5.0, true);
    OrientedBox3<double> merge{};
    bool ok = MergeContainers(box0, box1, merge);
    io.outBool(ok);
    OutOBox<3>(io, merge);
}

// ======================= ContEllipse2 / ContEllipsoid3 ==================

namespace
{
    template <int N>
    void OutEllip(oracle::Ctx& io, Hyperellipsoid<N, double> const& e)
    {
        io.outVec(e.center);
        for (int i = 0; i < N; ++i) { io.outVec(e.axis[i]); }
        io.outVec(e.extent);
    }
}

ORACLE_CASE("ContEllipse2.getContainer")
{
    // The degenerate clouds (kinds 2 to 4) give the fitter an exactly zero
    // eigenvalue, which upstream turns into an infinite extent
    // (docs/UPSTREAM-FINDINGS.md ContEllipse2.h, issue #292); the defect is
    // preserved in the port, so infinity and NaN are compared bit for bit.
    int kind = io.index() % 5;
    int count = io.integer(1, 8);
    auto points = Cloud<2>(io, kind, count, 5, 8.0);
    Ellipse2<double> ellipse{};
    bool ok = GetContainer(count, points.data(), ellipse);
    io.outBool(ok);
    OutEllip<2>(io, ellipse);
}

ORACLE_CASE("ContEllipse2.inContainer")
{
    int mode = io.index() % 2;
    auto ellipse = Ellip<2>(io, mode, 4, 5.0, false);
    auto p = FrameProbe<2>(io, ellipse.center, ellipse.axis, ellipse.extent, mode);
    io.outBool(InContainer(p, ellipse));
}

ORACLE_CASE("ContEllipse2.inContainer.leftHanded")
{
    int mode = io.index() % 2;
    auto ellipse = Ellip<2>(io, mode, 4, 5.0, true);
    auto p = FrameProbe<2>(io, ellipse.center, ellipse.axis, ellipse.extent, mode);
    io.outBool(InContainer(p, ellipse));
}

ORACLE_CASE("ContEllipse2.mergeContainers")
{
    int mode = io.index() % 2;
    auto e0 = Ellip<2>(io, mode, 4, 5.0, false);
    auto e1 = Ellip<2>(io, mode, 4, 5.0, false);
    Ellipse2<double> merge{};
    bool ok = MergeContainers(e0, e1, merge);
    io.outBool(ok);
    OutEllip<2>(io, merge);
}

ORACLE_CASE("ContEllipsoid3.getContainer")
{
    int kind = io.index() % 5;
    int count = io.integer(1, 8);
    auto points = Cloud<3>(io, kind, count, 5, 8.0);
    Ellipsoid3<double> ellipsoid{};
    bool ok = GetContainer(count, points.data(), ellipsoid);
    io.outBool(ok);
    OutEllip<3>(io, ellipsoid);
}

ORACLE_CASE("ContEllipsoid3.inContainer")
{
    int mode = io.index() % 2;
    auto e = Ellip<3>(io, mode, 4, 5.0, false);
    auto p = FrameProbe<3>(io, e.center, e.axis, e.extent, mode);
    io.outBool(InContainer(p, e));
}

ORACLE_CASE("ContEllipsoid3.inContainer.leftHanded")
{
    int mode = io.index() % 2;
    auto e = Ellip<3>(io, mode, 4, 5.0, true);
    auto p = FrameProbe<3>(io, e.center, e.axis, e.extent, mode);
    io.outBool(InContainer(p, e));
}

ORACLE_CASE("ContEllipsoid3.mergeContainers")
{
    int mode = io.index() % 2;
    auto e0 = Ellip<3>(io, mode, 4, 5.0, false);
    auto e1 = Ellip<3>(io, mode, 4, 5.0, false);
    Ellipsoid3<double> merge{};
    bool ok = MergeContainers(e0, e1, merge);
    io.outBool(ok);
    OutEllip<3>(io, merge);
}

ORACLE_CASE("ContEllipsoid3.mergeContainers.leftHanded")
{
    // As ContOrientedBox3.mergeContainers.leftHanded: the quaternion average
    // assumes det = +1 and the defect is preserved in the port.
    int mode = io.index() % 2;
    auto e0 = Ellip<3>(io, mode, 4, 5.0, true);
    auto e1 = Ellip<3>(io, mode, 4, 5.0, true);
    Ellipsoid3<double> merge{};
    bool ok = MergeContainers(e0, e1, merge);
    io.outBool(ok);
    OutEllip<3>(io, merge);
}

// =================================== ContCone ===========================

namespace
{
    // A cone whose six angle-derived constants are recorded as inputs, so the
    // C math library never enters the compared computation. 'kind' selects
    // infinite (0), infinite truncated (1), finite (2) or frustum (3).
    // Records origin (N), direction (N), angle and its six derived values,
    // then minHeight and maxHeight (-1 means "infinite").
    template <int N>
    Cone<N, double> Cn(oracle::Ctx& io, int mode, int lat, double range, int kind)
    {
        Cone<N, double> cone{};
        cone.ray.origin = Pt<N>(io, mode, lat, range);
        cone.ray.direction = Dir<N>(io, mode);
        double angle = io.raw(0.15, 1.35);
        cone.SetAngle(angle);
        io.given(cone.angle);
        io.given(cone.cosAngle);
        io.given(cone.sinAngle);
        io.given(cone.tanAngle);
        io.given(cone.cosAngleSqr);
        io.given(cone.sinAngleSqr);
        io.given(cone.invSinAngle);
        double minHeight = 0.0, maxHeight = -1.0;
        if (kind == 1 || kind == 3)
        {
            minHeight = (mode == 0 ? static_cast<double>(io.rawInteger(1, lat))
                : io.raw(0.25, range));
        }
        if (kind == 2 || kind == 3)
        {
            maxHeight = minHeight + (mode == 0
                ? static_cast<double>(io.rawInteger(1, lat)) : io.raw(0.25, range));
        }
        io.given(minHeight);
        io.given(maxHeight);
        if (maxHeight < 0.0) { cone.MakeInfiniteTruncatedCone(minHeight); }
        else { cone.MakeConeFrustum(minHeight, maxHeight); }
        return cone;
    }
}

ORACLE_CASE("ContCone.inContainer.3d")
{
    int mode = io.index() % 2;
    int kind = (io.index() / 2) % 4;
    auto cone = Cn<3>(io, mode, 4, 5.0, kind);
    Vector3<double> p{};
    int choice = io.rawInteger(0, 2);
    if (choice == 0)
    {
        // Exactly on the cone axis at a height equal to an interval endpoint,
        // so HeightInRange is evaluated at equality.
        double h = (kind == 2 || kind == 3 ? cone.GetMaxHeight() : cone.GetMinHeight());
        p = cone.ray.origin + h * cone.ray.direction;
    }
    else
    {
        for (int k = 0; k < 3; ++k)
        {
            p[k] = (mode == 0 ? static_cast<double>(io.rawInteger(-6, 6))
                : io.raw(-7.0, 7.0));
        }
    }
    io.givenVec(p);
    io.outBool(InContainer(p, cone));
}

ORACLE_CASE("ContCone.inContainer.2d")
{
    int mode = io.index() % 2;
    int kind = (io.index() / 2) % 4;
    auto cone = Cn<2>(io, mode, 4, 5.0, kind);
    Vector2<double> p{};
    for (int k = 0; k < 2; ++k)
    {
        p[k] = (mode == 0 ? static_cast<double>(io.rawInteger(-6, 6))
            : io.raw(-7.0, 7.0));
    }
    io.givenVec(p);
    io.outBool(InContainer(p, cone));
}

// ===================== ContScribeCircle2 / ...Circle3Sphere3 ============

ORACLE_CASE("ContScribeCircle2.circumscribe")
{
    // kinds 2 and 3 are coincident and collinear vertices, for which the
    // 2x2 system is singular and upstream returns false with the circle
    // untouched; nothing is emitted for it then.
    int kind = io.index() % 4;
    auto v = Cloud<2>(io, kind, 3, 5, 8.0);
    Circle2<double> circle{};
    bool ok = Circumscribe(v[0], v[1], v[2], circle);
    io.outBool(ok);
    if (ok) { OutSphere<2>(io, circle); }
}

ORACLE_CASE("ContScribeCircle2.inscribe")
{
    int kind = io.index() % 4;
    auto v = Cloud<2>(io, kind, 3, 5, 8.0);
    Circle2<double> circle{};
    bool ok = Inscribe(v[0], v[1], v[2], circle);
    io.outBool(ok);
    if (ok) { OutSphere<2>(io, circle); }
}

ORACLE_CASE("ContScribeCircle3Sphere3.circumscribeCircle3")
{
    int kind = io.index() % 5;
    auto v = Cloud<3>(io, kind, 3, 5, 8.0);
    Circle3<double> circle{};
    bool ok = Circumscribe(v[0], v[1], v[2], circle);
    io.outBool(ok);
    if (ok)
    {
        io.outVec(circle.center);
        io.outVec(circle.normal);
        io.outReal(circle.radius);
    }
}

ORACLE_CASE("ContScribeCircle3Sphere3.circumscribeSphere3")
{
    int kind = io.index() % 5;
    auto v = Cloud<3>(io, kind, 4, 5, 8.0);
    Sphere3<double> sphere{};
    bool ok = Circumscribe(v[0], v[1], v[2], v[3], sphere);
    io.outBool(ok);
    if (ok) { OutSphere<3>(io, sphere); }
}

ORACLE_CASE("ContScribeCircle3Sphere3.inscribeCircle3")
{
    int kind = io.index() % 5;
    auto v = Cloud<3>(io, kind, 3, 5, 8.0);
    Circle3<double> circle{};
    bool ok = Inscribe(v[0], v[1], v[2], circle);
    io.outBool(ok);
    if (ok)
    {
        io.outVec(circle.center);
        io.outVec(circle.normal);
        io.outReal(circle.radius);
    }
}

ORACLE_CASE("ContScribeCircle3Sphere3.inscribeSphere3")
{
    int kind = io.index() % 5;
    auto v = Cloud<3>(io, kind, 4, 5, 8.0);
    Sphere3<double> sphere{};
    bool ok = Inscribe(v[0], v[1], v[2], v[3], sphere);
    io.outBool(ok);
    if (ok) { OutSphere<3>(io, sphere); }
}

// ============================= ContTetrahedron3 =========================

ORACLE_CASE("ContTetrahedron3.inContainer")
{
    int kind = io.index() % 5;
    auto v = Cloud<3>(io, kind, 4, 4, 6.0);
    Tetrahedron3<double> tetra{};
    for (int i = 0; i < 4; ++i) { tetra.v[i] = v[static_cast<size_t>(i)]; }
    Vector3<double> p{};
    switch (io.rawInteger(0, 4))
    {
    case 0:
        // A vertex, where every triple scalar product of the three incident
        // faces is exactly zero.
        p = tetra.v[io.rawInteger(0, 3)];
        break;
    case 1:
    {
        // A face centroid: exactly on one face plane.
        int a = io.rawInteger(0, 3);
        int b = (a + 1) % 4;
        int c = (a + 2) % 4;
        p = (tetra.v[a] + tetra.v[b] + tetra.v[c]) / 3.0;
        break;
    }
    case 2:
    {
        // An edge midpoint.
        int a = io.rawInteger(0, 3);
        int b = (a + 1 + io.rawInteger(0, 2)) % 4;
        p = 0.5 * (tetra.v[a] + tetra.v[b]);
        break;
    }
    default:
        for (int k = 0; k < 3; ++k)
        {
            p[k] = (kind == 1 ? io.raw(-7.0, 7.0)
                : static_cast<double>(io.rawInteger(-5, 5)));
        }
        break;
    }
    io.givenVec(p);
    io.outBool(InContainer(p, tetra));
}

// ============================ ContEllipse2MinCR =========================

namespace
{
    struct MinCR2Draw
    {
        Vector2<double> C{};
        std::array<Vector2<double>, 2> axis{};
        std::vector<Vector2<double>> points{};
    };

    // Unrecorded draw. 'forceOnAxis' puts one point exactly on the first
    // ellipse axis, which makes that point's constraint coefficient A[1]
    // exactly zero: those are the vertical constraint lines on which
    // upstream's MaxProduct divides by zero.
    MinCR2Draw DrawMinCR2(oracle::Ctx& io, int mode, bool forceOnAxis)
    {
        MinCR2Draw d{};
        d.C = DrawPt<2>(io, mode, 4, 5.0);
        d.axis = DrawFrame2(io, mode);
        int count = io.rawInteger(3, 7);
        d.points.resize(static_cast<size_t>(count));
        for (int i = 0; i < count; ++i)
        {
            d.points[static_cast<size_t>(i)] = DrawPt<2>(io, mode, 5, 6.0);
        }
        if (forceOnAxis)
        {
            int i = io.rawInteger(0, count - 1);
            double t = static_cast<double>(io.rawInteger(1, 5));
            d.points[static_cast<size_t>(i)] = d.C + t * d.axis[0];
        }
        return d;
    }

    // Records C (2), the two columns of R (4), the point count and the points.
    void RecordMinCR2(oracle::Ctx& io, MinCR2Draw const& d)
    {
        io.givenVec(d.C);
        io.givenVec(d.axis[0]);
        io.givenVec(d.axis[1]);
        io.given(static_cast<double>(d.points.size()));
        for (auto const& p : d.points) { io.givenVec(p); }
    }

    // true when upstream returns a non-finite component for this draw. The
    // division by a zero A[iYMin][1] yields NaN when the numerator is also
    // zero and +-infinity otherwise; both are the defect.
    bool MinCR2IsNaN(MinCR2Draw const& d, double D[2])
    {
        Matrix2x2<double> R{};
        R.SetCol(0, d.axis[0]);
        R.SetCol(1, d.axis[1]);
        D[0] = 0.0;
        D[1] = 0.0;
        try
        {
            ContEllipse2MinCR<double> query;
            query(static_cast<int32_t>(d.points.size()), d.points.data(), d.C, R, D);
        }
        catch (std::exception const&) { return false; }
        return !std::isfinite(D[0]) || !std::isfinite(D[1]);
    }

    void RunMinCR2(oracle::Ctx& io, int mode, bool forceOnAxis, bool wantNaN,
        int attempts)
    {
        MinCR2Draw best = DrawMinCR2(io, mode, forceOnAxis);
        double D[2] = { 0.0, 0.0 };
        bool isNaN = MinCR2IsNaN(best, D);
        for (int attempt = 1; attempt < attempts && isNaN != wantNaN; ++attempt)
        {
            MinCR2Draw d = DrawMinCR2(io, mode, forceOnAxis);
            double tD[2] = { 0.0, 0.0 };
            bool tNaN = MinCR2IsNaN(d, tD);
            if (tNaN == wantNaN || attempt == attempts - 1)
            {
                best = d;
                D[0] = tD[0];
                D[1] = tD[1];
                isNaN = tNaN;
            }
        }
        RecordMinCR2(io, best);
        Matrix2x2<double> R{};
        R.SetCol(0, best.axis[0]);
        R.SetCol(1, best.axis[1]);
        double out[2] = { 0.0, 0.0 };
        ContEllipse2MinCR<double> query;
        query(static_cast<int32_t>(best.points.size()), best.points.data(),
            best.C, R, out);
        io.outReal(out[0]);
        io.outReal(out[1]);
    }
}

ORACLE_CASE("ContEllipse2MinCR.compute")
{
    // Inputs on which upstream is sound: the accepted draw is one whose result
    // has no NaN component. Upstream divides by A[iYMin][1] when the walk
    // steps onto a vertical constraint line, which is exactly where the port's
    // fix applies (docs/UPSTREAM-FINDINGS.md ContEllipse2MinCR.h, issue #234).
    // The rejection loop is capped at 24 attempts and redraws every quantity.
    int mode = io.index() % 2;
    RunMinCR2(io, mode, false, false, 24);
}

ORACLE_CASE("ContEllipse2MinCR.compute.verticalLineDeviation")
{
    // The complementary set: draws on which upstream returns a NaN component
    // because it evaluates (1 - a0*x0)/b0 with b0 == 0. The port evaluates the
    // shared hull vertex on the previous line of the walk instead (issue
    // #234), so these records deviate by construction.
    RunMinCR2(io, 0, true, true, 64);
}

// =========================== ContEllipsoid3MinCR ========================

namespace
{
    // A replica of upstream's facet/edge walk, used only to CLASSIFY a draw
    // before the real query is run on it. It reproduces upstream's control
    // flow exactly (including the mt19937 jitter, drawn from the same
    // std::uniform_real_distribution in this translation unit) and reports
    //   0 : upstream terminates and never sees a negative numerator
    //   1 : upstream's LogAssert(numer >= 0) fires, so upstream throws
    //   2 : the walk does not terminate; upstream overflows the stack, which
    //       is not a catchable C++ exception, so such draws must be rejected
    //       rather than recorded.
    // Over 50000 random integer clouds the classification agreed with the real
    // query on every draw (0 => no throw, 1 => throw).
    struct MinCR3Walk
    {
        std::vector<Vector3<double>> A{};
        int status = 0;
        int budget = 0;

        void Edge(int plane0, int plane1, double D[3])
        {
            if (--budget < 0) { status = 2; return; }
            double xDir = A[plane0][1] * A[plane1][2] - A[plane1][1] * A[plane0][2];
            double yDir = A[plane0][2] * A[plane1][0] - A[plane1][2] * A[plane0][0];
            double zDir = A[plane0][0] * A[plane1][1] - A[plane1][0] * A[plane0][1];
            double a0 = D[0] * D[1] * zDir + D[0] * D[2] * yDir + D[1] * D[2] * xDir;
            double a1 = 2.0 * (D[2] * xDir * yDir + D[1] * xDir * zDir
                + D[0] * yDir * zDir);
            double a2 = 3.0 * (xDir * yDir * zDir);
            double tFinal;
            if (a2 != 0.0)
            {
                double invA2 = 1.0 / a2;
                double discr = a1 * a1 - 4.0 * a0 * a2;
                discr = std::sqrt(std::max(discr, 0.0));
                tFinal = -0.5 * (a1 + discr) * invA2;
                if (a1 + 2.0 * a2 * tFinal > 0.0)
                {
                    tFinal = 0.5 * (-a1 + discr) * invA2;
                }
            }
            else if (a1 != 0.0) { tFinal = -a0 / a1; }
            else if (a0 != 0.0)
            {
                double fmax = std::numeric_limits<double>::max();
                tFinal = (a0 >= 0.0 ? fmax : -fmax);
            }
            else { return; }
            if (tFinal < 0.0)
            {
                tFinal = -tFinal; xDir = -xDir; yDir = -yDir; zDir = -zDir;
            }
            double tMax = tFinal;
            int plane2 = -1;
            int n = static_cast<int>(A.size());
            for (int i = 0; i < n; ++i)
            {
                if (i == plane0 || i == plane1) { continue; }
                double norDotDir = A[i][0] * xDir + A[i][1] * yDir + A[i][2] * zDir;
                if (norDotDir <= 0.0) { continue; }
                double numer = 1.0 - A[i][0] * D[0] - A[i][1] * D[1] - A[i][2] * D[2];
                if (numer < 0.0) { status = 1; return; }
                double t = numer / norDotDir;
                if (0 <= t && t < tMax) { plane2 = i; tMax = t; }
            }
            D[0] += tMax * xDir; D[1] += tMax * yDir; D[2] += tMax * zDir;
            if (tMax == tFinal) { return; }
            if (tMax > 0.0) { Facet(plane2, D); return; }
        }

        void Facet(int plane0, double D[3])
        {
            if (--budget < 0) { status = 2; return; }
            double tFinal, xDir, yDir, zDir;
            if (A[plane0][0] > 0.0 && A[plane0][1] > 0.0 && A[plane0][2] > 0.0)
            {
                double oneThird = 1.0 / 3.0;
                xDir = oneThird / A[plane0][0] - D[0];
                yDir = oneThird / A[plane0][1] - D[1];
                zDir = oneThird / A[plane0][2] - D[2];
                tFinal = 1.0;
            }
            else
            {
                tFinal = std::numeric_limits<double>::max();
                xDir = (A[plane0][0] > 0.0 ? 0.0 : 1.0);
                yDir = (A[plane0][1] > 0.0 ? 0.0 : 1.0);
                zDir = (A[plane0][2] > 0.0 ? 0.0 : 1.0);
            }
            double tMax = tFinal;
            int plane1 = -1;
            int n = static_cast<int>(A.size());
            for (int i = 0; i < n; ++i)
            {
                if (i == plane0) { continue; }
                double norDotDir = A[i][0] * xDir + A[i][1] * yDir + A[i][2] * zDir;
                if (norDotDir <= 0.0) { continue; }
                double numer = 1.0 - A[i][0] * D[0] - A[i][1] * D[1] - A[i][2] * D[2];
                if (numer < 0.0) { status = 1; return; }
                double t = numer / norDotDir;
                if (0 <= t && t < tMax) { plane1 = i; tMax = t; }
            }
            D[0] += tMax * xDir; D[1] += tMax * yDir; D[2] += tMax * zDir;
            if (tMax == 1.0) { return; }
            if (tMax > 0.0) { Facet(plane1, D); return; }
            Edge(plane0, plane1, D);
        }
    };

    struct MinCR3Draw
    {
        Vector3<double> C{};
        std::array<Vector3<double>, 3> axis{};
        std::vector<Vector3<double>> points{};
    };

    MinCR3Draw DrawMinCR3(oracle::Ctx& io, int mode, bool flatten)
    {
        MinCR3Draw d{};
        d.C = DrawPt<3>(io, mode, 4, 5.0);
        d.axis = DrawFrame3(io, mode);
        int count = io.rawInteger(4, 8);
        d.points.resize(static_cast<size_t>(count));
        for (int i = 0; i < count; ++i)
        {
            d.points[static_cast<size_t>(i)] = DrawPt<3>(io, mode, 5, 6.0);
            if (flatten)
            {
                // Remove the third coordinate in the frame, so every A[i][2]
                // is exactly zero before the jitter and the whole of
                // D[2] = 1/zmax comes from the jitter.
                Vector3<double>& p = d.points[static_cast<size_t>(i)];
                p = p - Dot(p - d.C, d.axis[2]) * d.axis[2];
            }
        }
        return d;
    }

    void RecordMinCR3(oracle::Ctx& io, MinCR3Draw const& d)
    {
        io.givenVec(d.C);
        io.givenVec(d.axis[0]);
        io.givenVec(d.axis[1]);
        io.givenVec(d.axis[2]);
        io.given(static_cast<double>(d.points.size()));
        for (auto const& p : d.points) { io.givenVec(p); }
    }

    int ClassifyMinCR3(MinCR3Draw const& d)
    {
        Matrix3x3<double> R{};
        R.SetCol(0, d.axis[0]);
        R.SetCol(1, d.axis[1]);
        R.SetCol(2, d.axis[2]);
        MinCR3Walk w;
        w.A.resize(d.points.size());
        for (size_t i = 0; i < d.points.size(); ++i)
        {
            Vector3<double> diff = d.points[i] - d.C;
            Vector3<double> prod = diff * R;
            w.A[i] = prod * prod;
        }
        std::mt19937 mte;
        std::uniform_real_distribution<double> rnd(0.0, 1.0);
        double maxJitter = 1e-12;
        int n = static_cast<int>(w.A.size());
        for (int i = 0; i < n; ++i)
        {
            w.A[i][0] += maxJitter * rnd(mte);
            w.A[i][1] += maxJitter * rnd(mte);
            w.A[i][2] += maxJitter * rnd(mte);
        }
        int plane = -1;
        double zmax = 0.0;
        for (int i = 0; i < n; ++i)
        {
            if (w.A[i][2] > zmax) { zmax = w.A[i][2]; plane = i; }
        }
        if (plane == -1) { return 1; }
        double D[3] = { 0.0, 0.0, 1.0 / zmax };
        w.budget = 4 * n + 64;
        w.Facet(plane, D);
        return w.status;
    }

    void RunMinCR3(oracle::Ctx& io, int mode, bool flatten, int wantStatus,
        int attempts)
    {
        MinCR3Draw best{};
        bool have = false;
        for (int attempt = 0; attempt < attempts; ++attempt)
        {
            MinCR3Draw d = DrawMinCR3(io, mode, flatten);
            int status = ClassifyMinCR3(d);
            if (status == wantStatus) { best = d; have = true; break; }
            // Status 2 draws crash upstream and can never be recorded, so the
            // fallback is always a status 0 draw.
            if (!have && status == 0) { best = d; have = true; }
        }
        if (!have)
        {
            // A cloud that is known to be status 0: four points in general
            // position around the centre.
            best.C = Vector3<double>{ 0.0, 0.0, 0.0 };
            best.axis[0] = Vector3<double>{ 1.0, 0.0, 0.0 };
            best.axis[1] = Vector3<double>{ 0.0, 1.0, 0.0 };
            best.axis[2] = Vector3<double>{ 0.0, 0.0, 1.0 };
            best.points = { Vector3<double>{ 1.0, 2.0, 3.0 },
                Vector3<double>{ -2.0, 1.0, 4.0 },
                Vector3<double>{ 3.0, -4.0, 1.0 },
                Vector3<double>{ -1.0, -2.0, -5.0 } };
        }
        RecordMinCR3(io, best);
        Matrix3x3<double> R{};
        R.SetCol(0, best.axis[0]);
        R.SetCol(1, best.axis[1]);
        R.SetCol(2, best.axis[2]);
        double D[3] = { 0.0, 0.0, 0.0 };
        ContEllipsoid3MinCR<double> query;
        query(static_cast<int32_t>(best.points.size()), best.points.data(),
            best.C, R, D);
        io.outReal(D[0]);
        io.outReal(D[1]);
        io.outReal(D[2]);
    }
}

ORACLE_CASE("ContEllipsoid3MinCR.compute")
{
    // Inputs on which upstream is sound (classifier status 0). Every fourth
    // record flattens the cloud into the plane of the first two axes, so that
    // every A[i][2] is exactly zero and D[2] = 1/zmax is decided entirely by
    // the mt19937 jitter; that is the configuration in which the exact
    // std::generate_canonical construction matters.
    int mode = io.index() % 2;
    bool flatten = (io.index() / 2) % 2 == 1;
    RunMinCR3(io, mode, flatten, 0, 40);
}

ORACLE_CASE("ContEllipsoid3MinCR.compute.assertDeviation")
{
    // Classifier status 1: an already-active constraint plane has a slack of
    // about -1e-16 and upstream's LogAssert(numer >= 0) fires, contradicting
    // the adjacent comment ("some numerical error may make this a small
    // negative number. In that case set tmax = 0"). The port clamps as the
    // comment prescribes (docs/UPSTREAM-FINDINGS.md ContEllipsoid3MinCR.h,
    // issue #409), so it returns where upstream throws.
    int mode = io.index() % 2;
    RunMinCR3(io, mode, false, 1, 40);
}

// ========================= ContPointInPolyhedron3 =======================

namespace
{
    struct PolyFace
    {
        std::vector<int32_t> indices{};
        std::vector<int32_t> triangles{};
        Vector3<double> normal{};
        double constant = 0.0;
    };

    struct Polyhedron
    {
        std::vector<Vector3<double>> vertices{};
        std::vector<PolyFace> faces{};
    };

    // Orient one face so that its normal points away from 'inside' and its
    // vertices are counterclockwise when viewed from outside, which is the
    // convention PointInPolyhedron3 documents.
    void OrientFace(Polyhedron& poly, PolyFace& face, Vector3<double> const& inside)
    {
        Vector3<double> const& v0 = poly.vertices[static_cast<size_t>(face.indices[0])];
        Vector3<double> const& v1 = poly.vertices[static_cast<size_t>(face.indices[1])];
        Vector3<double> const& v2 = poly.vertices[static_cast<size_t>(face.indices[2])];
        Vector3<double> n = UnitCross(v1 - v0, v2 - v0);
        if (Dot(n, inside - v0) > 0.0)
        {
            n = -n;
            std::reverse(face.indices.begin(), face.indices.end());
        }
        face.normal = n;
        face.constant = Dot(n, poly.vertices[static_cast<size_t>(face.indices[0])]);
    }

    // Fill face.triangles with a triangle fan of the (already oriented) face.
    void FanFace(PolyFace& face)
    {
        face.triangles.clear();
        for (size_t k = 1; k + 1 < face.indices.size(); ++k)
        {
            face.triangles.push_back(face.indices[0]);
            face.triangles.push_back(face.indices[k]);
            face.triangles.push_back(face.indices[k + 1]);
        }
    }

    // A tetrahedron (shape 0) or an axis-aligned box (shape 1). 'quads' keeps
    // the box faces as quadrilaterals; otherwise they are split into
    // triangles. Unrecorded.
    Polyhedron DrawPolyhedron(oracle::Ctx& io, int shape, bool quads)
    {
        Polyhedron poly{};
        if (shape == 0)
        {
            for (int attempt = 0; attempt < 16; ++attempt)
            {
                poly.vertices.clear();
                for (int i = 0; i < 4; ++i)
                {
                    poly.vertices.push_back(DrawPt<3>(io, 0, 4, 4.0));
                }
                double vol = DotCross(poly.vertices[1] - poly.vertices[0],
                    poly.vertices[2] - poly.vertices[0],
                    poly.vertices[3] - poly.vertices[0]);
                if (std::fabs(vol) > 0.5) { break; }
                if (attempt == 15)
                {
                    poly.vertices = { Vector3<double>{ 0.0, 0.0, 0.0 },
                        Vector3<double>{ 4.0, 0.0, 0.0 },
                        Vector3<double>{ 0.0, 4.0, 0.0 },
                        Vector3<double>{ 0.0, 0.0, 4.0 } };
                }
            }
            int32_t const tri[4][3] = { { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 }, { 1, 2, 3 } };
            for (auto const& t : tri)
            {
                PolyFace f{};
                f.indices = { t[0], t[1], t[2] };
                poly.faces.push_back(f);
            }
        }
        else
        {
            Vector3<double> c = DrawPt<3>(io, 0, 3, 3.0);
            double e[3];
            for (int k = 0; k < 3; ++k)
            {
                e[k] = static_cast<double>(io.rawInteger(1, 4));
            }
            for (int i = 0; i < 8; ++i)
            {
                Vector3<double> v{};
                for (int k = 0; k < 3; ++k)
                {
                    v[k] = c[k] + ((i >> k) & 1 ? e[k] : -e[k]);
                }
                poly.vertices.push_back(v);
            }
            int32_t const quad[6][4] =
            {
                { 0, 2, 3, 1 }, { 4, 5, 7, 6 }, { 0, 1, 5, 4 },
                { 2, 6, 7, 3 }, { 0, 4, 6, 2 }, { 1, 3, 7, 5 }
            };
            for (auto const& q : quad)
            {
                if (quads)
                {
                    PolyFace f{};
                    f.indices = { q[0], q[1], q[2], q[3] };
                    poly.faces.push_back(f);
                }
                else
                {
                    PolyFace f0{}, f1{};
                    f0.indices = { q[0], q[1], q[2] };
                    f1.indices = { q[0], q[2], q[3] };
                    poly.faces.push_back(f0);
                    poly.faces.push_back(f1);
                }
            }
        }

        Vector3<double> inside{ 0.0, 0.0, 0.0 };
        for (auto const& v : poly.vertices) { inside = inside + v; }
        inside = inside / static_cast<double>(poly.vertices.size());
        for (auto& f : poly.faces)
        {
            OrientFace(poly, f, inside);
            FanFace(f);
        }
        return poly;
    }

    void RecordPolyhedron(oracle::Ctx& io, Polyhedron const& poly)
    {
        io.given(static_cast<double>(poly.vertices.size()));
        for (auto const& v : poly.vertices) { io.givenVec(v); }
        io.given(static_cast<double>(poly.faces.size()));
        for (auto const& f : poly.faces)
        {
            io.given(static_cast<double>(f.indices.size()));
            for (int32_t i : f.indices) { io.given(static_cast<double>(i)); }
            io.givenVec(f.normal);
            io.given(f.constant);
            io.given(static_cast<double>(f.triangles.size()));
            for (int32_t i : f.triangles) { io.given(static_cast<double>(i)); }
        }
    }

    // Builds upstream's Face array. Face::indices is std::array<int32_t,3>, so
    // only the first three indices of each face fit: that truncation IS the
    // defect of issue #343, which the port fixes with a dynamic array.
    std::vector<PointInPolyhedron3<double>::Face> UpstreamFaces(Polyhedron const& poly)
    {
        std::vector<PointInPolyhedron3<double>::Face> faces(poly.faces.size());
        for (size_t i = 0; i < poly.faces.size(); ++i)
        {
            for (size_t k = 0; k < 3; ++k)
            {
                faces[i].indices[k] = poly.faces[i].indices[k];
            }
            faces[i].plane.normal = poly.faces[i].normal;
            faces[i].plane.constant = poly.faces[i].constant;
            faces[i].triangles = poly.faces[i].triangles;
        }
        return faces;
    }

    // Records the number of rays (1 or 3) and their unit directions.
    std::vector<Vector3<double>> RayDirections(oracle::Ctx& io, int mode)
    {
        int numRays = 1 + 2 * io.integer(0, 1);
        std::vector<Vector3<double>> dirs;
        for (int j = 0; j < numRays; ++j) { dirs.push_back(Dir<3>(io, mode)); }
        return dirs;
    }

    // A test point: the centroid, a vertex, a face centroid or a lattice point.
    Vector3<double> PolyhedronProbe(oracle::Ctx& io, Polyhedron const& poly)
    {
        Vector3<double> p{};
        switch (io.rawInteger(0, 3))
        {
        case 0:
            for (auto const& v : poly.vertices) { p = p + v; }
            p = p / static_cast<double>(poly.vertices.size());
            break;
        case 1:
            p = poly.vertices[static_cast<size_t>(
                io.rawInteger(0, static_cast<int>(poly.vertices.size()) - 1))];
            break;
        case 2:
        {
            PolyFace const& f = poly.faces[static_cast<size_t>(
                io.rawInteger(0, static_cast<int>(poly.faces.size()) - 1))];
            for (int32_t i : f.indices)
            {
                p = p + poly.vertices[static_cast<size_t>(i)];
            }
            p = p / static_cast<double>(f.indices.size());
            break;
        }
        default:
            p = DrawPt<3>(io, 0, 6, 6.0);
            break;
        }
        return io.givenVec(p);
    }

    void RunPolyhedron(oracle::Ctx& io, int shape, bool quads,
        PointInPolyhedron3<double>::FaceType type, uint32_t method)
    {
        int mode = io.index() % 2;
        Polyhedron poly = DrawPolyhedron(io, shape, quads);
        RecordPolyhedron(io, poly);
        std::vector<Vector3<double>> dirs = RayDirections(io, mode);
        Vector3<double> p = PolyhedronProbe(io, poly);
        std::vector<PointInPolyhedron3<double>::Face> faces = UpstreamFaces(poly);
        PointInPolyhedron3<double> query(type, static_cast<int32_t>(poly.vertices.size()),
            poly.vertices.data(), static_cast<int32_t>(faces.size()), faces.data(),
            static_cast<int32_t>(dirs.size()), dirs.data(), method);
        io.outBool(query.Contains(p));
    }
}

ORACLE_CASE("ContPointInPolyhedron3.contains.triangle")
{
    RunPolyhedron(io, io.index() % 2, false,
        PointInPolyhedron3<double>::TRIANGLE, 0);
}

ORACLE_CASE("ContPointInPolyhedron3.contains.convex0")
{
    RunPolyhedron(io, io.index() % 2, false,
        PointInPolyhedron3<double>::CONVEX, 0);
}

ORACLE_CASE("ContPointInPolyhedron3.contains.convex1")
{
    RunPolyhedron(io, io.index() % 2, false,
        PointInPolyhedron3<double>::CONVEX, 1);
}

ORACLE_CASE("ContPointInPolyhedron3.contains.convex2")
{
    RunPolyhedron(io, io.index() % 2, false,
        PointInPolyhedron3<double>::CONVEX, 2);
}

ORACLE_CASE("ContPointInPolyhedron3.contains.simple0")
{
    // Method 0 for simple faces reads Face::triangles, a std::vector, so the
    // quadrilateral faces of the box are NOT truncated and upstream and the
    // port agree.
    RunPolyhedron(io, 1, true, PointInPolyhedron3<double>::SIMPLE, 0);
}

ORACLE_CASE("ContPointInPolyhedron3.contains.simple1")
{
    RunPolyhedron(io, io.index() % 2, false,
        PointInPolyhedron3<double>::SIMPLE, 1);
}

ORACLE_CASE("ContPointInPolyhedron3.contains.unsupported")
{
    // SIMPLE with method >= 2 is not implemented and upstream returns a silent
    // 'false' (preserved in the port).
    RunPolyhedron(io, io.index() % 2, false,
        PointInPolyhedron3<double>::SIMPLE, 2);
}

ORACLE_CASE("ContPointInPolyhedron3.contains.convex0QuadDeviation")
{
    // Quadrilateral faces: upstream's Face::indices is std::array<int32_t,3>
    // while ContainsC0 reads indices.size() as the face's vertex count, so the
    // trifan covers only the first three vertices of each quad and half of the
    // box surface is missing (docs/UPSTREAM-FINDINGS.md
    // ContPointInPolyhedron3.h, issue #343). The port keeps the full index
    // list.
    RunPolyhedron(io, 1, true, PointInPolyhedron3<double>::CONVEX, 0);
}

ORACLE_CASE("ContPointInPolyhedron3.contains.convex12QuadDeviation")
{
    // The same truncation in SharedContains: only the first three vertices of
    // each quad are projected, so the point-in-polygon test runs against a
    // triangle. Alternates the O(N) and O(log N) polygon tests.
    RunPolyhedron(io, 1, true, PointInPolyhedron3<double>::CONVEX,
        1 + static_cast<uint32_t>(io.index() % 2));
}

ORACLE_CASE("ContPointInPolyhedron3.contains.simple1QuadDeviation")
{
    RunPolyhedron(io, 1, true, PointInPolyhedron3<double>::SIMPLE, 1);
}

