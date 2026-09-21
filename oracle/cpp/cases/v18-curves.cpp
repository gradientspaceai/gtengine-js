// Verify group 18 (curves): differential cases for NaturalSplineCurve.h,
// RiemannianGeodesic.h and EllipsoidGeodesic.h.
//
// Comparison policy per header:
//
//  * NaturalSplineCurve is + - * / only (the closed variant additionally runs
//    Gaussian elimination with full pivoting through LinearSystem::Solve), so
//    every case is compared bit for bit.
//  * RiemannianGeodesic is an abstract class whose only concrete subclass in
//    GTE is EllipsoidGeodesic, whose metric calls sin and cos. To compare the
//    geodesic algorithm itself - the steepest-descent search of Refine, the
//    subdivision bookkeeping of ComputeGeodesic, the Trapezoid Rule
//    integrators, the metric inverse and the Christoffel symbols of the second
//    kind - bit for bit, this file defines its own concrete subclass
//    PolySurfaceGeodesic over the polynomial graph z = u^2*v + u*v^2. Its
//    metric and Christoffel symbols are polynomials, so the whole algorithm,
//    including every comparison that decides a branch, is exact IEEE
//    arithmetic and the cases are declared { exact: true }.
//  * EllipsoidGeodesic's own cases call sin and cos, so they carry a
//    tolerance. Refine's steepest-descent search picks the argument of the
//    minimum over 2*searchSamples+1 sampled lengths and ComputeSegmentLength's
//    values are libm-derived, so that choice is not comparable at any
//    tolerance. The cases that reach it either set searchSamples = 0, which
//    makes the search a single sample at tRay = 0 and leaves the midpoint
//    exactly where Subdivide put it (no libm value reaches the output, so the
//    case is exact), or run a probe that rejects the inputs on which the
//    winning sample is not separated from the others; see
//    EllipsoidSearchSeparation.
//
// Upstream computes the derived parameters mIntegralStep, mSearchStep and
// mDerivativeFactor in the constructor only, so a case that changes
// integralSamples or searchSamples afterwards is comparing upstream's stale
// derived values. That quirk is preserved by the port (which offers
// updateDerivedParameters() as an opt-in escape hatch that these cases never
// call) and is exercised on purpose: integralSamples is varied while
// mIntegralStep stays 1/15.
#define ORACLE_FAMILY "v18-curves"
#include "Oracle.h"

#include <Mathematics/EllipsoidGeodesic.h>
#include <Mathematics/NaturalSplineCurve.h>
#include <Mathematics/RiemannianGeodesic.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <vector>

using namespace gte;

namespace
{
    // ---- NaturalSplineCurve generators ------------------------------------
    // Records 1 + 3*n + n doubles: the point count, the points and the
    // strictly increasing times, in every mode.
    struct SplineInput
    {
        std::vector<Vector3<double>> points;
        std::vector<double> times;
    };

    SplineInput MakeSpline(oracle::Ctx& io, int minPoints, int maxPoints)
    {
        SplineInput s{};
        int n = io.integer(minPoints, maxPoints);
        s.points.resize(static_cast<size_t>(n));
        s.times.resize(static_cast<size_t>(n));

        int pointMode = io.index() % 2;
        for (int i = 0; i < n; ++i)
        {
            s.points[static_cast<size_t>(i)] = (pointMode == 0
                ? io.vec<3>(-4.0, 4.0) : io.latticeVec<3>(-4, 4));
        }

        int timeMode = (io.index() / 2) % 3;
        if (timeMode == 0)
        {
            // Uniformly spaced knots.
            double t0 = io.raw(-2.0, 2.0);
            double h = io.raw(0.5, 2.0);
            for (int i = 0; i < n; ++i)
            {
                s.times[static_cast<size_t>(i)] =
                    io.given(t0 + static_cast<double>(i) * h);
            }
        }
        else if (timeMode == 1)
        {
            // Nonuniformly spaced knots.
            double t = io.raw(-2.0, 2.0);
            for (int i = 0; i < n; ++i)
            {
                s.times[static_cast<size_t>(i)] = io.given(t);
                t = t + io.raw(0.25, 2.0);
            }
        }
        else
        {
            // Integer knots, which make dt and the tridiagonal recurrence
            // exact and reach the degenerate configurations.
            double t = static_cast<double>(io.rawInteger(-3, 3));
            for (int i = 0; i < n; ++i)
            {
                s.times[static_cast<size_t>(i)] = io.given(t);
                t = t + static_cast<double>(io.rawInteger(1, 3));
            }
        }
        return s;
    }

    // Five evaluation parameters, recorded in this order: one below the first
    // knot, one above the last knot, one exactly at a knot, and a pair
    // straddling an interior knot (which the replay uses to check C2
    // continuity independently). The straddling pair requires n >= 3; for
    // n == 2 the interior knot is replaced by the midpoint of the single
    // segment.
    std::vector<double> MakeEvalTimes(oracle::Ctx& io, SplineInput const& s)
    {
        size_t n = s.times.size();
        double tmin = s.times[0];
        double tmax = s.times[n - 1];
        std::vector<double> t(5);
        t[0] = io.given(tmin - io.raw(0.1, 2.0));
        t[1] = io.given(tmax + io.raw(0.1, 2.0));
        size_t k = static_cast<size_t>(io.rawInteger(0, static_cast<int>(n) - 1));
        t[2] = io.given(s.times[k]);
        double interior = (n >= 3
            ? s.times[static_cast<size_t>(io.rawInteger(1, static_cast<int>(n) - 2))]
            : 0.5 * (tmin + tmax));
        double delta = io.raw(1e-7, 1e-6);
        t[3] = io.given(interior - delta);
        t[4] = io.given(interior + delta);
        return t;
    }

    void EmitSplineJets(oracle::Ctx& io, NaturalSplineCurve<3, double> const& curve,
        std::vector<double> const& times, int order)
    {
        for (size_t j = 0; j < times.size(); ++j)
        {
            Vector3<double> jet[5]{};
            curve.Evaluate(times[j], static_cast<uint32_t>(order), jet);
            for (int i = 0; i <= order; ++i) { io.outVec(jet[i]); }
        }
    }
}

// ---- NaturalSplineCurve ----------------------------------------------------

ORACLE_CASE("NaturalSplineCurve.free.evaluate")
{
    SplineInput s = MakeSpline(io, 3, 6);
    std::vector<double> t = MakeEvalTimes(io, s);
    int order = io.integer(0, 4);
    NaturalSplineCurve<3, double> curve(true, static_cast<int32_t>(s.points.size()),
        s.points.data(), s.times.data());
    io.outInt(curve.GetNumPoints());
    EmitSplineJets(io, curve, t, order);
}

// The closed spline solves a numPoints-by-numPoints linear system whose last
// row wraps around. Upstream writes that row with three plain assignments to
// columns numSegments-1, 0 and 1, which are distinct only when
// numSegments >= 3; the port accumulates instead (issue #295). This case is
// restricted to numPoints >= 4 so that upstream is sound.
ORACLE_CASE("NaturalSplineCurve.closed.evaluate")
{
    SplineInput s = MakeSpline(io, 4, 6);
    std::vector<double> t = MakeEvalTimes(io, s);
    int order = io.integer(0, 4);
    NaturalSplineCurve<3, double> curve(false, static_cast<int32_t>(s.points.size()),
        s.points.data(), s.times.data());
    io.outInt(curve.GetNumPoints());
    EmitSplineJets(io, curve, t, order);
}

// The deliberate port fix of issue #295: for numPoints 2 and 3 the wrap-around
// columns coincide and upstream's third assignment discards the first term.
ORACLE_CASE("NaturalSplineCurve.closed.wrapRow.deviation")
{
    SplineInput s = MakeSpline(io, 2, 3);
    std::vector<double> t = MakeEvalTimes(io, s);
    int order = io.integer(0, 4);
    NaturalSplineCurve<3, double> curve(false, static_cast<int32_t>(s.points.size()),
        s.points.data(), s.times.data());
    io.outInt(curve.GetNumPoints());
    EmitSplineJets(io, curve, t, order);
}

ORACLE_CASE("NaturalSplineCurve.clamped.evaluate")
{
    SplineInput s = MakeSpline(io, 3, 6);
    std::vector<double> t = MakeEvalTimes(io, s);
    int order = io.integer(0, 4);
    Vector3<double> derivative0 = io.vec<3>(-3.0, 3.0);
    Vector3<double> derivative1 = io.vec<3>(-3.0, 3.0);
    NaturalSplineCurve<3, double> curve(static_cast<int32_t>(s.points.size()),
        s.points.data(), s.times.data(), derivative0, derivative1);
    io.outInt(curve.GetNumPoints());
    EmitSplineJets(io, curve, t, order);
}

// Two points is the documented minimum; the free and clamped constructors
// reduce to a single segment.
ORACLE_CASE("NaturalSplineCurve.twoPoints.evaluate")
{
    SplineInput s = MakeSpline(io, 2, 2);
    std::vector<double> t = MakeEvalTimes(io, s);
    int order = io.integer(0, 3);
    bool isFree = io.boolean();
    Vector3<double> derivative0 = io.vec<3>(-3.0, 3.0);
    Vector3<double> derivative1 = io.vec<3>(-3.0, 3.0);
    if (isFree)
    {
        NaturalSplineCurve<3, double> curve(true,
            static_cast<int32_t>(s.points.size()), s.points.data(), s.times.data());
        io.outInt(curve.GetNumPoints());
        EmitSplineJets(io, curve, t, order);
    }
    else
    {
        NaturalSplineCurve<3, double> curve(static_cast<int32_t>(s.points.size()),
            s.points.data(), s.times.data(), derivative0, derivative1);
        io.outInt(curve.GetNumPoints());
        EmitSplineJets(io, curve, t, order);
    }
}

// Throw parity for LogAssert(numPoints >= 2). ParametricCurve's constructor
// runs first with numSegments = 0, which allocates a one-element time array
// and copies one time, so nothing is read out of bounds; the assert in
// NaturalSplineCurve then fires. numPoints = 0 is not generated: it would make
// ParametricCurve allocate a vector of SIZE_MAX elements.
ORACLE_CASE("NaturalSplineCurve.construct.tooFewPoints")
{
    Vector3<double> point = io.vec<3>(-4.0, 4.0);
    double time = io.real(-2.0, 2.0);
    NaturalSplineCurve<3, double> curve(true, 1, &point, &time);
    io.outInt(curve.GetNumPoints());
}

// ---- RiemannianGeodesic on a polynomial surface ----------------------------
// The graph z = u^2*v + u*v^2 of R^2 into R^3. With P_u = (1, 0, du) and
// P_v = (0, 1, dv), the metric is the Gram matrix of (P_u, P_v), which is
// positive definite everywhere (its determinant is 1 + du^2 + dv^2), and the
// Christoffel symbols of the first kind are Gamma_{k,ij} = Dot(P_ij, P_k),
// upstream's convention. Everything is polynomial, so the entire geodesic
// algorithm runs in exact IEEE arithmetic on both sides, including the
// comparisons that decide its branches.

namespace
{
    class PolySurfaceGeodesic : public RiemannianGeodesic<double>
    {
    public:
        PolySurfaceGeodesic(int32_t dimension)
            :
            RiemannianGeodesic<double>(dimension)
        {
        }

        virtual ~PolySurfaceGeodesic() = default;

        // Public access to the protected helpers, in the order
        // ComputeIntegrand calls them, plus ComputeMetricDerivative, which
        // the algorithm itself never calls.
        bool ProbeTensors(GVector<double> const& point)
        {
            ComputeMetric(point);
            ComputeChristoffel1(point);
            bool invertible = ComputeMetricInverse();
            ComputeChristoffel2();
            ComputeMetricDerivative();
            return invertible;
        }

        double ProbeIntegrand(GVector<double> const& pos, GVector<double> const& der)
        {
            return ComputeIntegrand(pos, der);
        }

        GMatrix<double> const& Metric() const { return mMetric; }
        GMatrix<double> const& MetricInverse() const { return mMetricInverse; }
        GMatrix<double> const& Christoffel1(int32_t i) const { return mChristoffel1[i]; }
        GMatrix<double> const& Christoffel2(int32_t i) const { return mChristoffel2[i]; }
        GMatrix<double> const& MetricDerivative(int32_t i) const { return mMetricDerivative[i]; }

        virtual void ComputeMetric(GVector<double> const& point) override
        {
            double u = point[0];
            double v = point[1];
            double du = 2.0 * u * v + v * v;
            double dv = u * u + 2.0 * u * v;
            mMetric(0, 0) = 1.0 + du * du;
            mMetric(0, 1) = du * dv;
            mMetric(1, 0) = mMetric(0, 1);
            mMetric(1, 1) = 1.0 + dv * dv;
        }

        virtual void ComputeChristoffel1(GVector<double> const& point) override
        {
            double u = point[0];
            double v = point[1];
            double du = 2.0 * u * v + v * v;
            double dv = u * u + 2.0 * u * v;
            double duu = 2.0 * v;
            double duv = 2.0 * u + 2.0 * v;
            double dvv = 2.0 * u;
            mChristoffel1[0](0, 0) = duu * du;
            mChristoffel1[0](0, 1) = duv * du;
            mChristoffel1[0](1, 0) = mChristoffel1[0](0, 1);
            mChristoffel1[0](1, 1) = dvv * du;
            mChristoffel1[1](0, 0) = duu * dv;
            mChristoffel1[1](0, 1) = duv * dv;
            mChristoffel1[1](1, 0) = mChristoffel1[1](0, 1);
            mChristoffel1[1](1, 1) = dvv * dv;
        }
    };

    // A parameter point of the polynomial surface. Two recorded doubles.
    GVector<double> PolyPoint(oracle::Ctx& io, int mode)
    {
        GVector<double> p(2);
        if (mode == 0)
        {
            p[0] = io.real(-2.0, 2.0);
            p[1] = io.real(-2.0, 2.0);
        }
        else
        {
            p[0] = io.lattice(-2, 2);
            p[1] = io.lattice(-2, 2);
        }
        return p;
    }

    // A polyline of 'n' parameter points whose consecutive points differ by at
    // least 0.5 in one coordinate, so the length and curvature integrands
    // never see a zero difference vector. Records 2*n doubles.
    std::vector<GVector<double>> PolyPath(oracle::Ctx& io, int mode, int n)
    {
        std::vector<GVector<double>> path(static_cast<size_t>(n));
        double u = (mode == 0 ? io.raw(-2.0, 2.0)
            : static_cast<double>(io.rawInteger(-2, 2)));
        double v = (mode == 0 ? io.raw(-2.0, 2.0)
            : static_cast<double>(io.rawInteger(-2, 2)));
        for (int i = 0; i < n; ++i)
        {
            path[static_cast<size_t>(i)].SetSize(2);
            path[static_cast<size_t>(i)][0] = io.given(u);
            path[static_cast<size_t>(i)][1] = io.given(v);
            if (mode == 0)
            {
                u = u + (io.raw(0.0, 1.0) < 0.5 ? io.raw(0.5, 1.5) : io.raw(-1.5, -0.5));
                v = v + (io.raw(0.0, 1.0) < 0.5 ? io.raw(0.5, 1.5) : io.raw(-1.5, -0.5));
            }
            else
            {
                u = u + static_cast<double>(io.rawInteger(1, 2));
                v = v + static_cast<double>(io.rawInteger(-2, -1));
            }
        }
        return path;
    }

    void EmitMatrix2(oracle::Ctx& io, GMatrix<double> const& m)
    {
        for (int32_t r = 0; r < 2; ++r)
        {
            for (int32_t c = 0; c < 2; ++c) { io.outReal(m(r, c)); }
        }
    }

    void EmitGVec(oracle::Ctx& io, GVector<double> const& v)
    {
        for (int32_t i = 0; i < v.GetSize(); ++i) { io.outReal(v[i]); }
    }
}

// ComputeMetricInverse (Gaussian elimination with full pivoting on the 2x2
// GMatrix), ComputeChristoffel2 and ComputeMetricDerivative. The last one is
// never called by the algorithm and preserves an upstream quirk: it computes
// 2*Gamma_{d,i0i1}, not dg_{i0i1}/dx_d (issue #295).
ORACLE_CASE("RiemannianGeodesic.poly.tensors")
{
    int mode = io.index() % 2;
    GVector<double> point = PolyPoint(io, mode);
    PolySurfaceGeodesic manifold(2);
    bool invertible = manifold.ProbeTensors(point);
    io.outBool(invertible);
    EmitMatrix2(io, manifold.Metric());
    EmitMatrix2(io, manifold.MetricInverse());
    EmitMatrix2(io, manifold.Christoffel1(0));
    EmitMatrix2(io, manifold.Christoffel1(1));
    EmitMatrix2(io, manifold.Christoffel2(0));
    EmitMatrix2(io, manifold.Christoffel2(1));
    EmitMatrix2(io, manifold.MetricDerivative(0));
    EmitMatrix2(io, manifold.MetricDerivative(1));
}

ORACLE_CASE("RiemannianGeodesic.poly.computeIntegrand")
{
    int mode = io.index() % 2;
    std::vector<GVector<double>> path = PolyPath(io, mode, 2);
    PolySurfaceGeodesic manifold(2);
    GVector<double> der = path[1] - path[0];
    double curvature = manifold.ProbeIntegrand(path[0], der);
    io.outReal(curvature);
}

ORACLE_CASE("RiemannianGeodesic.poly.computeSegmentLength")
{
    int mode = io.index() % 2;
    std::vector<GVector<double>> path = PolyPath(io, mode, 2);
    int integralSamples = io.integer(2, 16);
    PolySurfaceGeodesic manifold(2);
    manifold.integralSamples = integralSamples;
    double length = manifold.ComputeSegmentLength(path[0], path[1]);
    io.outReal(length);
}

ORACLE_CASE("RiemannianGeodesic.poly.computeTotalLength")
{
    int mode = io.index() % 2;
    int quantity = io.integer(2, 5);
    std::vector<GVector<double>> path = PolyPath(io, mode, quantity);
    PolySurfaceGeodesic manifold(2);
    double length = manifold.ComputeTotalLength(quantity, path);
    io.outReal(length);
}

ORACLE_CASE("RiemannianGeodesic.poly.computeSegmentCurvature")
{
    int mode = io.index() % 2;
    std::vector<GVector<double>> path = PolyPath(io, mode, 2);
    PolySurfaceGeodesic manifold(2);
    double curvature = manifold.ComputeSegmentCurvature(path[0], path[1]);
    io.outReal(curvature);
}

ORACLE_CASE("RiemannianGeodesic.poly.computeTotalCurvature")
{
    int mode = io.index() % 2;
    int quantity = io.integer(2, 5);
    std::vector<GVector<double>> path = PolyPath(io, mode, quantity);
    PolySurfaceGeodesic manifold(2);
    double curvature = manifold.ComputeTotalCurvature(quantity, path);
    io.outReal(curvature);
}

// The steepest-descent search. searchSamples and derivativeStep are set after
// construction, so upstream's mSearchStep is 1/32 and its mDerivativeFactor is
// 5000 whatever the case chooses; that stale-derived-parameter quirk is
// preserved and compared.
ORACLE_CASE("RiemannianGeodesic.poly.refine")
{
    int mode = io.index() % 2;
    std::vector<GVector<double>> path = PolyPath(io, mode, 2);
    GVector<double> mid(2);
    mid[0] = io.real(-2.0, 2.0);
    mid[1] = io.real(-2.0, 2.0);
    int searchSamples = io.integer(1, 4);
    double searchRadius = io.real(0.25, 2.0);
    double derivativeStep = io.real(1e-5, 1e-3);
    PolySurfaceGeodesic manifold(2);
    manifold.searchSamples = searchSamples;
    manifold.searchRadius = searchRadius;
    manifold.derivativeStep = derivativeStep;
    bool changed = manifold.Refine(path[0], mid, path[1]);
    io.outBool(changed);
    EmitGVec(io, mid);
}

ORACLE_CASE("RiemannianGeodesic.poly.subdivide")
{
    int mode = io.index() % 2;
    std::vector<GVector<double>> path = PolyPath(io, mode, 2);
    int searchSamples = io.integer(1, 4);
    double searchRadius = io.real(0.25, 2.0);
    PolySurfaceGeodesic manifold(2);
    manifold.searchSamples = searchSamples;
    manifold.searchRadius = searchRadius;
    GVector<double> mid(2);
    bool changed = manifold.Subdivide(path[0], mid, path[1]);
    io.outBool(changed);
    EmitGVec(io, mid);
}

// The full multiresolution algorithm. The direct segment length is emitted
// beside the total length of the polyline so that the replay can check the
// algorithm's own contract (the refined polyline is no longer than the
// straight parameter-space segment) without a second implementation.
ORACLE_CASE("RiemannianGeodesic.poly.computeGeodesic")
{
    int mode = io.index() % 2;
    std::vector<GVector<double>> path = PolyPath(io, mode, 2);
    int subdivisions = io.integer(1, 2);
    int refinements = io.integer(0, 2);
    int searchSamples = io.integer(1, 2);
    double searchRadius = io.real(0.25, 2.0);
    PolySurfaceGeodesic manifold(2);
    manifold.subdivisions = subdivisions;
    manifold.refinements = refinements;
    manifold.searchSamples = searchSamples;
    manifold.searchRadius = searchRadius;
    int32_t quantity = 0;
    std::vector<GVector<double>> geodesic;
    manifold.ComputeGeodesic(path[0], path[1], quantity, geodesic);
    io.outInt(quantity);
    for (int32_t i = 0; i < quantity; ++i) { EmitGVec(io, geodesic[static_cast<size_t>(i)]); }
    io.outReal(manifold.ComputeTotalLength(quantity, geodesic));
    io.outReal(manifold.ComputeSegmentLength(path[0], path[1]));
    io.outReal(manifold.ComputeTotalCurvature(quantity, geodesic));
}

// refineCallback and the three progress accessors. Subdivide installs a no-op
// callback for the Refine it performs itself and restores the caller's
// afterwards, so the recorded trace has one entry per Refine of the
// refinement loop and none from the subdivision pass.
ORACLE_CASE("RiemannianGeodesic.poly.refineCallback")
{
    int mode = io.index() % 2;
    std::vector<GVector<double>> path = PolyPath(io, mode, 2);
    int subdivisions = io.integer(1, 2);
    int refinements = io.integer(0, 2);
    int searchSamples = io.integer(1, 2);
    PolySurfaceGeodesic manifold(2);
    manifold.subdivisions = subdivisions;
    manifold.refinements = refinements;
    manifold.searchSamples = searchSamples;
    std::vector<int32_t> trace;
    manifold.refineCallback = [&manifold, &trace]()
    {
        trace.push_back(manifold.GetSubdivisionStep());
        trace.push_back(manifold.GetRefinementStep());
        trace.push_back(manifold.GetCurrentQuantity());
    };
    int32_t quantity = 0;
    std::vector<GVector<double>> geodesic;
    manifold.ComputeGeodesic(path[0], path[1], quantity, geodesic);
    io.outInt(static_cast<int32_t>(trace.size()));
    for (size_t i = 0; i < trace.size(); ++i) { io.outInt(trace[i]); }
    io.outInt(manifold.GetSubdivisionStep());
    io.outInt(manifold.GetRefinementStep());
    io.outInt(manifold.GetCurrentQuantity());
}

// Throw parity for LogAssert(dimension >= 2).
ORACLE_CASE("RiemannianGeodesic.construct.invalidDimension")
{
    int dimension = io.integer(0, 1);
    PolySurfaceGeodesic manifold(dimension);
    io.outInt(manifold.GetDimension());
}

// ---- EllipsoidGeodesic -----------------------------------------------------
// sin and cos in ComputeMetric and ComputeChristoffel1, so these cases carry
// a tolerance. The v parameter stays in [0.6, 2.5], away from the poles where
// sin(v) vanishes and the metric degenerates, so that LogAssert(qForm > 0)
// cannot be decided by libm rounding.

namespace
{
    // Every third record is a sphere (equal extents), which the replay checks
    // against the closed-form sphere metric. Three recorded doubles.
    void EllipsoidExtents(oracle::Ctx& io, int mode, double& a, double& b, double& c)
    {
        if (mode == 0)
        {
            double r = io.raw(0.5, 3.0);
            a = io.given(r);
            b = io.given(r);
            c = io.given(r);
        }
        else
        {
            a = io.real(0.5, 3.0);
            b = io.real(0.5, 3.0);
            c = io.real(0.5, 3.0);
        }
    }

    // A polyline of parameter points (u, v) whose consecutive points are at
    // least 0.4 apart in u and in v. Records 2*n doubles.
    std::vector<GVector<double>> EllipsoidPath(oracle::Ctx& io, int n)
    {
        std::vector<GVector<double>> path(static_cast<size_t>(n));
        double u = io.raw(-3.0, 3.0);
        double v = io.raw(0.8, 2.2);
        double du = io.raw(0.4, 1.0);
        double dv = io.raw(0.15, 0.3);
        for (int i = 0; i < n; ++i)
        {
            path[static_cast<size_t>(i)].SetSize(2);
            path[static_cast<size_t>(i)][0] = io.given(u + static_cast<double>(i) * du);
            path[static_cast<size_t>(i)][1] = io.given(v
                + static_cast<double>(i % 2) * dv);
        }
        return path;
    }
}

ORACLE_CASE("EllipsoidGeodesic.computePosition")
{
    int mode = io.index() % 3;
    double a = 0.0, b = 0.0, c = 0.0;
    EllipsoidExtents(io, mode, a, b, c);
    std::vector<GVector<double>> path = EllipsoidPath(io, 1);
    EllipsoidGeodesic<double> eg(a, b, c);
    Vector3<double> position = eg.ComputePosition(path[0]);
    io.outVec(position);
}

ORACLE_CASE("EllipsoidGeodesic.computeSegmentLength")
{
    int mode = io.index() % 3;
    double a = 0.0, b = 0.0, c = 0.0;
    EllipsoidExtents(io, mode, a, b, c);
    std::vector<GVector<double>> path = EllipsoidPath(io, 2);
    int integralSamples = io.integer(2, 16);
    EllipsoidGeodesic<double> eg(a, b, c);
    eg.integralSamples = integralSamples;
    io.outReal(eg.ComputeSegmentLength(path[0], path[1]));
}

ORACLE_CASE("EllipsoidGeodesic.computeTotalLength")
{
    int mode = io.index() % 3;
    double a = 0.0, b = 0.0, c = 0.0;
    EllipsoidExtents(io, mode, a, b, c);
    int quantity = io.integer(2, 5);
    std::vector<GVector<double>> path = EllipsoidPath(io, quantity);
    EllipsoidGeodesic<double> eg(a, b, c);
    io.outReal(eg.ComputeTotalLength(quantity, path));
}

ORACLE_CASE("EllipsoidGeodesic.computeSegmentCurvature")
{
    int mode = io.index() % 3;
    double a = 0.0, b = 0.0, c = 0.0;
    EllipsoidExtents(io, mode, a, b, c);
    std::vector<GVector<double>> path = EllipsoidPath(io, 2);
    EllipsoidGeodesic<double> eg(a, b, c);
    io.outReal(eg.ComputeSegmentCurvature(path[0], path[1]));
}

ORACLE_CASE("EllipsoidGeodesic.computeTotalCurvature")
{
    int mode = io.index() % 3;
    double a = 0.0, b = 0.0, c = 0.0;
    EllipsoidExtents(io, mode, a, b, c);
    int quantity = io.integer(2, 5);
    std::vector<GVector<double>> path = EllipsoidPath(io, quantity);
    EllipsoidGeodesic<double> eg(a, b, c);
    io.outReal(eg.ComputeTotalCurvature(quantity, path));
}

// searchSamples = 0 leaves the steepest-descent search with the single sample
// tRay = 0, whose candidate point is bit-identical to the midpoint Subdivide
// produced. Whichever way the libm-decided comparison goes, the midpoint that
// leaves Refine is the arithmetic midpoint, so no libm value reaches the path
// and the case is exact. It compares ComputeGeodesic's own bookkeeping: the
// vertex count, the interleaving copy and the subdivision order.
ORACLE_CASE("EllipsoidGeodesic.computeGeodesic.noSearch")
{
    int mode = io.index() % 3;
    double a = 0.0, b = 0.0, c = 0.0;
    EllipsoidExtents(io, mode, a, b, c);
    std::vector<GVector<double>> path = EllipsoidPath(io, 2);
    int subdivisions = io.integer(1, 3);
    int refinements = io.integer(0, 2);
    EllipsoidGeodesic<double> eg(a, b, c);
    eg.subdivisions = subdivisions;
    eg.refinements = refinements;
    eg.searchSamples = 0;
    int32_t quantity = 0;
    std::vector<GVector<double>> geodesic;
    eg.ComputeGeodesic(path[0], path[1], quantity, geodesic);
    io.outInt(quantity);
    for (int32_t i = 0; i < quantity; ++i)
    {
        EmitGVec(io, geodesic[static_cast<size_t>(i)]);
    }
}

namespace
{
    // Replicates the sampled lengths of RiemannianGeodesic::Refine and returns
    // the smallest relative gap between the length the search selects and any
    // other sampled length, including the length at the current midpoint (the
    // value that decides Refine's boolean result). A negative return means the
    // configuration throws or is degenerate. mDerivativeFactor and mSearchStep
    // are upstream's constructor-time values, 0.5/1e-4 and 1/32, because the
    // case changes searchSamples only after construction.
    double EllipsoidSearchSeparation(EllipsoidGeodesic<double>& eg,
        GVector<double> const& end0, GVector<double> const& mid,
        GVector<double> const& end1, int searchSamples, double searchRadius)
    {
        double const derivativeStep = 1e-4;
        double const derivativeFactor = 0.5 / derivativeStep;
        double const searchStep = 1.0 / 32.0;

        GVector<double> temp = mid;
        GVector<double> gradient(2);
        for (int32_t i = 0; i < 2; ++i)
        {
            temp[i] = mid[i] + derivativeStep;
            gradient[i] = eg.ComputeSegmentLength(end0, temp);
            gradient[i] += eg.ComputeSegmentLength(temp, end1);
            temp[i] = mid[i] - derivativeStep;
            gradient[i] -= eg.ComputeSegmentLength(end0, temp);
            gradient[i] -= eg.ComputeSegmentLength(temp, end1);
            temp[i] = mid[i];
            gradient[i] *= derivativeFactor;
        }

        double oldLength = eg.ComputeSegmentLength(end0, mid);
        oldLength += eg.ComputeSegmentLength(mid, end1);

        double multiplier = searchStep * searchRadius;
        std::vector<double> lengths;
        for (int32_t i = -searchSamples; i <= searchSamples; ++i)
        {
            double tRay = multiplier * static_cast<double>(i);
            GVector<double> pRay = mid - tRay * gradient;
            double value = eg.ComputeSegmentLength(end0, pRay);
            value += eg.ComputeSegmentLength(end1, pRay);
            lengths.push_back(value);
        }

        // Refine keeps the first strictly smaller sample, so the winner is the
        // first index attaining the running minimum.
        double best = oldLength;
        size_t bestIndex = lengths.size();
        for (size_t i = 0; i < lengths.size(); ++i)
        {
            if (lengths[i] < best) { best = lengths[i]; bestIndex = i; }
        }

        double gap = std::numeric_limits<double>::max();
        for (size_t i = 0; i < lengths.size(); ++i)
        {
            if (i == bestIndex) { continue; }
            double d = std::fabs(lengths[i] - best) / std::max(1.0, std::fabs(best));
            gap = std::min(gap, d);
        }
        if (bestIndex != lengths.size())
        {
            double d = std::fabs(oldLength - best) / std::max(1.0, std::fabs(oldLength));
            gap = std::min(gap, d);
        }
        return gap;
    }
}

// Refine on the real ellipsoid. The search picks the argument of the minimum
// over 2*searchSamples+1 libm-derived lengths, which no tolerance can rescue,
// so the generator rejects (capped at 24 attempts, redrawing every quantity)
// the configurations where the winning length is within a relative 1e-9 of any
// other sampled length or of the length at the current midpoint. On the
// accepted inputs both builds select the same sample and the same boolean, and
// the emitted point differs only by the rounding of the gradient, which the
// case's 1e-12 tolerance covers.
ORACLE_CASE("EllipsoidGeodesic.refine.separated")
{
    bool sphere = (io.index() % 3 == 0);
    double ba = 1.0, bb = 1.0, bc = 1.0, bRadius = 1.0, bestGap = -1.0;
    int bSamples = 1;
    GVector<double> bEnd0(2), bMid(2), bEnd1(2);
    bEnd0[0] = 0.0; bEnd0[1] = 1.0;
    bMid[0] = 0.5; bMid[1] = 1.25;
    bEnd1[0] = 1.0; bEnd1[1] = 1.5;
    for (int attempt = 0; attempt < 24; ++attempt)
    {
        double a = io.raw(0.5, 3.0);
        double b = (sphere ? a : io.raw(0.5, 3.0));
        double c = (sphere ? a : io.raw(0.5, 3.0));
        GVector<double> end0(2), mid(2), end1(2);
        end0[0] = io.raw(-3.0, 3.0);
        end0[1] = io.raw(0.8, 2.3);
        end1[0] = end0[0] + io.raw(0.5, 1.5);
        end1[1] = end0[1] + io.raw(-0.4, 0.4);
        mid[0] = 0.5 * (end0[0] + end1[0]) + io.raw(-0.1, 0.1);
        mid[1] = 0.5 * (end0[1] + end1[1]) + io.raw(-0.1, 0.1);
        int samples = io.rawInteger(1, 3);
        double radius = io.raw(0.25, 2.0);
        if (mid[1] < 0.6 || mid[1] > 2.5 || end1[1] < 0.6 || end1[1] > 2.5)
        {
            continue;
        }
        double gap = -1.0;
        EllipsoidGeodesic<double> probe(a, b, c);
        try
        {
            gap = EllipsoidSearchSeparation(probe, end0, mid, end1, samples, radius);
        }
        catch (std::exception const&)
        {
            gap = -1.0;
        }
        if (gap > bestGap)
        {
            bestGap = gap;
            ba = a; bb = b; bc = c;
            bEnd0 = end0; bMid = mid; bEnd1 = end1;
            bSamples = samples; bRadius = radius;
        }
        if (gap > 1e-9) { break; }
    }

    io.given(ba);
    io.given(bb);
    io.given(bc);
    io.given(bEnd0[0]);
    io.given(bEnd0[1]);
    io.given(bMid[0]);
    io.given(bMid[1]);
    io.given(bEnd1[0]);
    io.given(bEnd1[1]);
    io.given(static_cast<double>(bSamples));
    io.given(bRadius);

    EllipsoidGeodesic<double> eg(ba, bb, bc);
    eg.searchSamples = bSamples;
    eg.searchRadius = bRadius;
    GVector<double> mid = bMid;
    bool changed = eg.Refine(bEnd0, mid, bEnd1);
    io.outBool(changed);
    EmitGVec(io, mid);
}
