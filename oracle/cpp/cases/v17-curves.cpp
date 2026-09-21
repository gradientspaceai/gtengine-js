// Group 17 (curves): B-spline basis functions, parametric curves and the
// frames, samplers and polyline tools built on them.
//
// Almost everything here is a pipeline of + - * / and sqrt, so the cases are
// declared exact on the TypeScript side. The exceptions, each called out at
// the case, are:
//   * FrenetFrame2/3::GetCurvature and ImplicitCurve2::GetCurvature call
//     std::pow(x, 1.5),
//   * SampleCircularArc calls std::acos, but only to decide how many samples
//     to produce; the generator rejects draws whose sample count sits within
//     0.1 of an integer boundary so that the (libm-dependent) truncation is
//     the same on both sides, and the sample values themselves are computed
//     with + - * / and sqrt alone.
//
// Division note: Vector.h's operator/=(Vector, Real) multiplies by
// 1/scalar (and yields the zero vector for scalar == 0), so an upstream
// 'v / s' is NOT a componentwise division. The port's div() has the same
// semantics; the cases below reach every such site (NURBS homogeneous
// divides, TCB derivative scaling, natural-spline derivative scaling,
// SampleCircularArc, NURBSCircle, NURBSSphere).
#define ORACLE_FAMILY "v17-curves"
#include "Oracle.h"

#include <Mathematics/Arc2.h>
#include <Mathematics/BasisFunction.h>
#include <Mathematics/BezierCurve.h>
#include <Mathematics/CLODPolyline.h>
#include <Mathematics/DarbouxFrame.h>
#include <Mathematics/FrenetFrame.h>
#include <Mathematics/ImplicitCurve2.h>
#include <Mathematics/MassSpringCurve.h>
#include <Mathematics/NURBSCircle.h>
#include <Mathematics/NURBSCurve.h>
#include <Mathematics/NURBSSphere.h>
#include <Mathematics/NaturalCubicSpline.h>
#include <Mathematics/NaturalQuinticSpline.h>
#include <Mathematics/ParametricCurve.h>
#include <Mathematics/PolylineOffset.h>
#include <Mathematics/ReparameterizeByArclength.h>
#include <Mathematics/SampleCircularArc.h>
#include <Mathematics/TCBSplineCurve.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

using namespace gte;

namespace
{
    // ---------------------------------------------------------- basis input

    // A drawn (not yet recorded) BasisFunctionInput. The knot values are
    // dyadic so that the knot differences and the derived domain are exact.
    struct BasisSpec
    {
        int32_t numControls = 0;
        int32_t degree = 0;
        bool uniform = false;
        bool periodic = false;
        std::vector<UniqueKnot<double>> knots{};
    };

    // Modes:
    //   0  open uniform, built exactly as BasisFunctionInput(n, d) builds it
    //   1  open nonuniform, dyadic interior knots of multiplicity 1
    //   2  open nonuniform with a repeated interior knot (multiplicity <= d)
    //   3  periodic, uniform dyadic knots
    //   4  periodic, nonuniform dyadic knots
    BasisSpec DrawBasis(oracle::Ctx& io, int32_t minDegree, int32_t maxDegree,
        int32_t mode)
    {
        BasisSpec s{};
        s.degree = io.rawInteger(minDegree, maxDegree);
        s.periodic = (mode >= 3);
        s.numControls = s.degree + 1 + io.rawInteger(0, 5);

        if (!s.periodic)
        {
            // Total multiplicity is numControls + degree + 1; the boundary
            // knots take degree + 1 each, leaving numControls - degree - 1
            // for the interior.
            int32_t interiorTotal = s.numControls - s.degree - 1;
            std::vector<int32_t> mult{};
            if (mode == 2 && interiorTotal >= 2)
            {
                // Put multiplicity 2 (at most degree) on one interior knot.
                int32_t remaining = interiorTotal;
                int32_t repeated = (s.degree >= 2 ? 2 : 1);
                mult.push_back(repeated);
                remaining -= repeated;
                for (int32_t i = 0; i < remaining; ++i) { mult.push_back(1); }
            }
            else
            {
                for (int32_t i = 0; i < interiorTotal; ++i) { mult.push_back(1); }
            }

            s.uniform = (mode == 0);
            if (mode == 0)
            {
                // BasisFunctionInput(numControls, degree): t in [0,1] with
                // equally spaced interior knots i / (numUniqueKnots - 1).
                int32_t numUniqueKnots = s.numControls - s.degree + 1;
                s.knots.resize(static_cast<size_t>(numUniqueKnots));
                s.knots.front().t = 0.0;
                s.knots.front().multiplicity = s.degree + 1;
                for (int32_t i = 1; i <= numUniqueKnots - 2; ++i)
                {
                    s.knots[i].t = i / (numUniqueKnots - 1.0);
                    s.knots[i].multiplicity = 1;
                }
                s.knots.back().t = 1.0;
                s.knots.back().multiplicity = s.degree + 1;
            }
            else
            {
                double u = 0.0;
                s.knots.push_back({ u, s.degree + 1 });
                for (size_t i = 0; i < mult.size(); ++i)
                {
                    u += io.rawInteger(1, 3) * 0.125;
                    s.knots.push_back({ u, mult[i] });
                }
                u += io.rawInteger(1, 3) * 0.125;
                s.knots.push_back({ u, s.degree + 1 });
            }
        }
        else
        {
            // A periodic basis has mNumControls = numControls + degree, so
            // the knot count is numControls + 2*degree + 1. Give every
            // unique knot multiplicity 1.
            s.uniform = (mode == 3);
            int32_t numUniqueKnots = s.numControls + 2 * s.degree + 1;
            s.knots.resize(static_cast<size_t>(numUniqueKnots));
            double u = 0.0;
            for (int32_t i = 0; i < numUniqueKnots; ++i)
            {
                s.knots[i].t = u;
                s.knots[i].multiplicity = 1;
                u += (mode == 3 ? 0.25 : io.rawInteger(1, 3) * 0.125);
            }
        }
        return s;
    }

    void RecordBasis(oracle::Ctx& io, BasisSpec const& s)
    {
        io.given(static_cast<double>(s.numControls));
        io.given(static_cast<double>(s.degree));
        io.given(s.uniform ? 1.0 : 0.0);
        io.given(s.periodic ? 1.0 : 0.0);
        io.given(static_cast<double>(s.knots.size()));
        for (auto const& k : s.knots)
        {
            io.given(k.t);
            io.given(static_cast<double>(k.multiplicity));
        }
    }

    BasisFunctionInput<double> ToInput(BasisSpec const& s)
    {
        BasisFunctionInput<double> input{};
        input.numControls = s.numControls;
        input.degree = s.degree;
        input.uniform = s.uniform;
        input.periodic = s.periodic;
        input.numUniqueKnots = static_cast<int32_t>(s.knots.size());
        input.uniqueKnots = s.knots;
        return input;
    }

    BasisSpec MakeBasis(oracle::Ctx& io, int32_t minDegree, int32_t maxDegree,
        int32_t numModes = 5)
    {
        int32_t mode = io.rawInteger(0, numModes - 1);
        BasisSpec s = DrawBasis(io, minDegree, maxDegree, mode);
        RecordBasis(io, s);
        return s;
    }

    // A parameter for a curve on [tmin, tmax]. The modes visit the domain
    // endpoints, the interior knots (including repeated ones), dyadic
    // midpoints, parameters just outside the domain (clamping) and
    // parameters far outside (periodic wrapping).
    double DrawParameter(oracle::Ctx& io, double tmin, double tmax,
        std::vector<UniqueKnot<double>> const* knots)
    {
        double tlength = tmax - tmin;
        switch (io.rawInteger(0, 7))
        {
        case 0: return tmin;
        case 1: return tmax;
        case 2:
            if (knots != nullptr && knots->size() > 2)
            {
                size_t i = static_cast<size_t>(io.rawInteger(1,
                    static_cast<int32_t>(knots->size()) - 2));
                return (*knots)[i].t;
            }
            return tmin + 0.5 * tlength;
        case 3: return tmin + tlength * (io.rawInteger(1, 7) * 0.125);
        case 4: return tmin - 0.25 * tlength;
        case 5: return tmax + 0.25 * tlength;
        case 6: return tmin - tlength * (io.rawInteger(1, 5) + 0.375);
        default: return tmin + tlength * io.raw(0.0, 1.0);
        }
    }

    // ------------------------------------------------------- control points

    // Modes: 0 uniform cloud, 1 small integer lattice, 2 lattice with a
    // repeated point, 3 all points coincident (zero-speed curve).
    std::vector<Vector<3, double>> DrawControls3(oracle::Ctx& io, int32_t count,
        int32_t mode)
    {
        std::vector<Vector<3, double>> P(static_cast<size_t>(count));
        if (mode == 0)
        {
            for (int32_t i = 0; i < count; ++i)
            {
                for (int32_t j = 0; j < 3; ++j) { P[i][j] = io.raw(-5.0, 5.0); }
            }
        }
        else if (mode == 1)
        {
            for (int32_t i = 0; i < count; ++i)
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    P[i][j] = static_cast<double>(io.rawInteger(-4, 4));
                }
            }
        }
        else if (mode == 2)
        {
            for (int32_t i = 0; i < count; ++i)
            {
                for (int32_t j = 0; j < 3; ++j)
                {
                    P[i][j] = static_cast<double>(io.rawInteger(-3, 3));
                }
            }
            int32_t a = io.rawInteger(0, count - 1);
            int32_t b = io.rawInteger(0, count - 1);
            P[a] = P[b];
        }
        else
        {
            Vector<3, double> q{};
            for (int32_t j = 0; j < 3; ++j)
            {
                q[j] = static_cast<double>(io.rawInteger(-3, 3));
            }
            for (int32_t i = 0; i < count; ++i) { P[i] = q; }
        }
        return P;
    }

    std::vector<Vector<3, double>> RecordControls3(oracle::Ctx& io,
        std::vector<Vector<3, double>> const& P)
    {
        for (auto const& p : P) { io.givenVec(p); }
        return P;
    }

    std::vector<Vector<3, double>> MakeControls3(oracle::Ctx& io, int32_t count,
        int32_t numModes = 4)
    {
        int32_t mode = io.rawInteger(0, numModes - 1);
        return RecordControls3(io, DrawControls3(io, count, mode));
    }

    // Positive weights: exact small integers, dyadic fractions or uniform.
    std::vector<double> MakeWeights(oracle::Ctx& io, int32_t count)
    {
        std::vector<double> w(static_cast<size_t>(count));
        int32_t mode = io.rawInteger(0, 2);
        for (int32_t i = 0; i < count; ++i)
        {
            if (mode == 0) { w[i] = static_cast<double>(io.rawInteger(1, 4)); }
            else if (mode == 1) { w[i] = io.rawInteger(1, 8) * 0.125; }
            else { w[i] = io.raw(0.25, 4.0); }
        }
        for (int32_t i = 0; i < count; ++i) { io.given(w[i]); }
        return w;
    }

    // Strictly increasing times with dyadic differences.
    std::vector<double> MakeTimes(oracle::Ctx& io, int32_t count)
    {
        std::vector<double> t(static_cast<size_t>(count));
        int32_t mode = io.rawInteger(0, 1);
        double u = static_cast<double>(io.rawInteger(-4, 4));
        for (int32_t i = 0; i < count; ++i)
        {
            t[i] = u;
            u += (mode == 0 ? 1.0 : io.rawInteger(1, 6) * 0.25);
        }
        for (int32_t i = 0; i < count; ++i) { io.given(t[i]); }
        return t;
    }
}

// ============================================================ BasisFunction

ORACLE_CASE("BasisFunction.create")
{
    auto spec = MakeBasis(io, 1, 5);
    BasisFunction<double> bf(ToInput(spec));
    io.outInt(bf.GetNumControls());
    io.outInt(bf.GetDegree());
    io.outInt(bf.GetNumUniqueKnots());
    io.outInt(bf.GetNumKnots());
    io.outReal(bf.GetMinDomain());
    io.outReal(bf.GetMaxDomain());
    io.outBool(bf.IsOpen());
    io.outBool(bf.IsUniform());
    io.outBool(bf.IsPeriodic());
    auto const* unique = bf.GetUniqueKnots();
    for (int32_t i = 0; i < bf.GetNumUniqueKnots(); ++i)
    {
        io.outReal(unique[i].t);
        io.outInt(unique[i].multiplicity);
    }
    auto const* knots = bf.GetKnots();
    for (int32_t i = 0; i < bf.GetNumKnots(); ++i) { io.outReal(knots[i]); }
}

ORACLE_CASE("BasisFunction.evaluate")
{
    auto spec = MakeBasis(io, 1, 5);
    BasisFunction<double> bf(ToInput(spec));
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, bf.GetMinDomain(), bf.GetMaxDomain(),
        &spec.knots));

    int32_t imin = 0, imax = 0;
    bf.Evaluate(t, static_cast<uint32_t>(order), imin, imax);
    io.outInt(imin);
    io.outInt(imax);
    int32_t const numValues = bf.GetNumControls() + bf.GetDegree();
    for (int32_t o = 0; o <= order; ++o)
    {
        for (int32_t i = 0; i < numValues; ++i)
        {
            io.outReal(bf.GetValue(static_cast<uint32_t>(o), i));
        }
    }
}

// A second Evaluate on the same object: upstream only overwrites the current
// support, so GetValue outside [minIndex, maxIndex] returns values left by
// the first call (docs/UPSTREAM-FINDINGS.md, issue #415, preserved).
ORACLE_CASE("BasisFunction.evaluate.repeated")
{
    auto spec = MakeBasis(io, 1, 5);
    BasisFunction<double> bf(ToInput(spec));
    double t0 = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    double t1 = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    double tmin = bf.GetMinDomain(), tmax = bf.GetMaxDomain();
    double u0 = tmin + (tmax - tmin) * t0;
    double u1 = tmin + (tmax - tmin) * t1;

    int32_t imin = 0, imax = 0;
    bf.Evaluate(u0, 3, imin, imax);
    io.outInt(imin);
    io.outInt(imax);
    bf.Evaluate(u1, 3, imin, imax);
    io.outInt(imin);
    io.outInt(imax);
    int32_t const numValues = bf.GetNumControls() + bf.GetDegree();
    for (int32_t o = 0; o <= 3; ++o)
    {
        for (int32_t i = 0; i < numValues; ++i)
        {
            io.outReal(bf.GetValue(static_cast<uint32_t>(o), i));
        }
    }
}

// The BasisFunctionInput(numControls, degree) convenience constructor.
ORACLE_CASE("BasisFunctionInput.openUniform")
{
    int32_t degree = io.integer(1, 5);
    int32_t numControls = degree + io.integer(1, 7);
    BasisFunctionInput<double> input(numControls, degree);
    io.outInt(input.numControls);
    io.outInt(input.degree);
    io.outBool(input.uniform);
    io.outBool(input.periodic);
    io.outInt(input.numUniqueKnots);
    for (int32_t i = 0; i < input.numUniqueKnots; ++i)
    {
        io.outReal(input.uniqueKnots[i].t);
        io.outInt(input.uniqueKnots[i].multiplicity);
    }
}

// GetValue rejects an out-of-range index and an order of 4 or more.
ORACLE_CASE("BasisFunction.getValue.invalid")
{
    auto spec = MakeBasis(io, 1, 3);
    BasisFunction<double> bf(ToInput(spec));
    int32_t which = io.integer(0, 1);
    int32_t imin = 0, imax = 0;
    bf.Evaluate(bf.GetMinDomain(), 1, imin, imax);
    if (which == 0)
    {
        io.outReal(bf.GetValue(0, bf.GetNumControls() + bf.GetDegree()));
    }
    else
    {
        io.outReal(bf.GetValue(4, 0));
    }
}

// ============================================================== BezierCurve

ORACLE_CASE("BezierCurve.evaluate")
{
    int32_t degree = io.integer(2, 5);
    auto controls = MakeControls3(io, degree + 1);
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));

    BezierCurve<3, double> curve(degree, controls.data());
    std::array<Vector<3, double>, 4> jet{};
    curve.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    io.outInt(curve.GetDegree());
    io.outInt(curve.GetNumControls());
    io.outBool(static_cast<bool>(curve));
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
    auto const* stored = curve.GetControls();
    for (int32_t i = 0; i < curve.GetNumControls(); ++i) { io.outVec(stored[i]); }
}

// order >= SUP_ORDER returns a zero jet.
ORACLE_CASE("BezierCurve.evaluate.invalidOrder")
{
    int32_t degree = io.integer(2, 4);
    auto controls = MakeControls3(io, degree + 1);
    double t = io.real(0.0, 1.0);

    BezierCurve<3, double> curve(degree, controls.data());
    std::array<Vector<3, double>, 4> jet{};
    curve.Evaluate(t, 4, jet.data());
    for (int32_t i = 0; i < 4; ++i) { io.outVec(jet[i]); }
}

// ========================================================== ParametricCurve
// The shared machinery is exercised through BezierCurve (one segment) and,
// for the multiple-segment paths, through NaturalCubicSpline below.

ORACLE_CASE("ParametricCurve.differentialGeometry")
{
    int32_t degree = io.integer(2, 5);
    auto controls = MakeControls3(io, degree + 1);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));

    BezierCurve<3, double> curve(degree, controls.data());
    io.outReal(curve.GetTMin());
    io.outReal(curve.GetTMax());
    io.outInt(curve.GetNumSegments());
    io.outVec(curve.GetPosition(t));
    io.outVec(curve.GetTangent(t));
    io.outReal(curve.GetSpeed(t));
}

ORACLE_CASE("ParametricCurve.setTimeInterval")
{
    int32_t degree = io.integer(2, 4);
    auto controls = MakeControls3(io, degree + 1);
    double tmin = io.real(-3.0, 0.0);
    double tmax = io.real(0.5, 3.0);
    double t = io.real(-1.0, 2.0);

    BezierCurve<3, double> curve(degree, controls.data());
    curve.SetTimeInterval(tmin, tmax);
    io.outReal(curve.GetTMin());
    io.outReal(curve.GetTMax());
    auto const* times = curve.GetTimes();
    io.outReal(times[0]);
    io.outReal(times[1]);
    // The domain of the Bezier evaluation is unchanged by SetTimeInterval.
    io.outVec(curve.GetPosition(t));
}

// GetLength, GetTotalLength: Romberg integration of the speed. The Romberg
// order is drawn small to keep the deep run fast; Romberg uses only + - * /.
ORACLE_CASE("ParametricCurve.getLength")
{
    int32_t degree = io.integer(2, 5);
    auto controls = MakeControls3(io, degree + 1);
    int32_t rombergOrder = io.integer(2, 5);
    double t0 = io.real(-0.25, 1.25);
    double t1 = io.real(-0.25, 1.25);

    BezierCurve<3, double> curve(degree, controls.data());
    curve.SetRombergOrder(rombergOrder);
    io.outReal(curve.GetLength(t0, t1));
    io.outReal(curve.GetTotalLength());
    io.outReal(curve.GetLength(t1, t0));
    io.outReal(curve.GetTotalLength());
}

ORACLE_CASE("ParametricCurve.getTime")
{
    int32_t degree = io.integer(2, 4);
    auto controls = MakeControls3(io, degree + 1);
    int32_t rombergOrder = io.integer(2, 4);
    int32_t maxBisections = io.integer(4, 16);
    double fraction = io.real(-0.25, 1.25);

    BezierCurve<3, double> curve(degree, controls.data());
    curve.SetRombergOrder(rombergOrder);
    curve.SetMaxBisections(static_cast<uint32_t>(maxBisections));
    double total = curve.GetTotalLength();
    double length = fraction * total;
    io.outReal(total);
    io.outReal(curve.GetTime(length));
}

ORACLE_CASE("ParametricCurve.subdivideByTime")
{
    int32_t degree = io.integer(2, 5);
    auto controls = MakeControls3(io, degree + 1);
    int32_t numPoints = io.integer(2, 6);

    BezierCurve<3, double> curve(degree, controls.data());
    std::vector<Vector<3, double>> points(static_cast<size_t>(numPoints));
    curve.SubdivideByTime(numPoints, points.data());
    for (auto const& p : points) { io.outVec(p); }
}

ORACLE_CASE("ParametricCurve.subdivideByLength")
{
    int32_t degree = io.integer(2, 4);
    auto controls = MakeControls3(io, degree + 1);
    int32_t rombergOrder = io.integer(2, 4);
    int32_t maxBisections = io.integer(4, 12);
    int32_t numPoints = io.integer(2, 5);

    BezierCurve<3, double> curve(degree, controls.data());
    curve.SetRombergOrder(rombergOrder);
    curve.SetMaxBisections(static_cast<uint32_t>(maxBisections));
    std::vector<Vector<3, double>> points(static_cast<size_t>(numPoints));
    curve.SubdivideByLength(numPoints, points.data());
    for (auto const& p : points) { io.outVec(p); }
}

// =============================================================== NURBSCurve

ORACLE_CASE("NURBSCurve.evaluate")
{
    auto spec = MakeBasis(io, 1, 4);
    auto controls = MakeControls3(io, spec.numControls);
    auto weights = MakeWeights(io, spec.numControls);
    int32_t order = io.integer(0, 3);

    NURBSCurve<3, double> curve(ToInput(spec), controls.data(), weights.data());
    double t = io.given(DrawParameter(io, curve.GetTMin(), curve.GetTMax(),
        &spec.knots));
    std::array<Vector<3, double>, 4> jet{};
    curve.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    io.outReal(curve.GetTMin());
    io.outReal(curve.GetTMax());
    io.outInt(curve.GetNumControls());
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

ORACLE_CASE("NURBSCurve.evaluate.invalidOrder")
{
    auto spec = MakeBasis(io, 1, 3);
    auto controls = MakeControls3(io, spec.numControls);
    auto weights = MakeWeights(io, spec.numControls);

    NURBSCurve<3, double> curve(ToInput(spec), controls.data(), weights.data());
    double t = io.given(DrawParameter(io, curve.GetTMin(), curve.GetTMax(),
        &spec.knots));
    std::array<Vector<3, double>, 4> jet{};
    curve.Evaluate(t, 4, jet.data());
    for (int32_t i = 0; i < 4; ++i) { io.outVec(jet[i]); }
}

// The control-point and weight accessors, including the out-of-range
// branches that return element 0.
ORACLE_CASE("NURBSCurve.accessors")
{
    auto spec = MakeBasis(io, 1, 3);
    auto controls = MakeControls3(io, spec.numControls);
    auto weights = MakeWeights(io, spec.numControls);
    int32_t setIndex = io.integer(-2, spec.numControls + 1);
    auto newControl = io.latticeVec<3>(-6, 6);
    double newWeight = io.real(0.5, 3.0);
    int32_t getIndex = io.integer(-2, spec.numControls + 1);

    NURBSCurve<3, double> curve(ToInput(spec), controls.data(), weights.data());
    curve.SetControl(setIndex, newControl);
    curve.SetWeight(setIndex, newWeight);
    io.outVec(curve.GetControl(getIndex));
    io.outReal(curve.GetWeight(getIndex));
    auto const* stored = curve.GetControls();
    auto const* storedW = curve.GetWeights();
    for (int32_t i = 0; i < curve.GetNumControls(); ++i)
    {
        io.outVec(stored[i]);
        io.outReal(storedW[i]);
    }
}

// Deferred controls and weights: the null-pointer constructor path.
ORACLE_CASE("NURBSCurve.deferred")
{
    auto spec = MakeBasis(io, 1, 3);
    NURBSCurve<3, double> curve(ToInput(spec), nullptr, nullptr);
    io.outInt(curve.GetNumControls());
    auto const* stored = curve.GetControls();
    auto const* storedW = curve.GetWeights();
    for (int32_t i = 0; i < curve.GetNumControls(); ++i)
    {
        io.outVec(stored[i]);
        io.outReal(storedW[i]);
    }
}

// ============================================================== NURBSCircle

namespace
{
    // A circular arc whose endpoints lie exactly on the circle. Mode 0 uses
    // the twelve lattice directions of the radius-5 circle (so the endpoints
    // are exact), mode 1 uses unrecorded angle draws whose cosine and sine
    // are recorded with givenVec. Both modes record center, radius and the
    // two endpoints, and both guarantee DotPerp(P0, P2) > 0, that is, a
    // counterclockwise arc subtending less than pi radians.
    constexpr int32_t kCircleLattice[12][2] =
    {
        { 5, 0 }, { 4, 3 }, { 3, 4 }, { 0, 5 }, { -3, 4 }, { -4, 3 },
        { -5, 0 }, { -4, -3 }, { -3, -4 }, { 0, -5 }, { 3, -4 }, { 4, -3 }
    };

    Arc2<double> MakeArc(oracle::Ctx& io, double radiusLo, double radiusHi)
    {
        Arc2<double> arc{};
        int32_t mode = io.rawInteger(0, 1);
        Vector2<double> center{}, d0{}, d1{};
        double radius = 0.0;
        for (int32_t attempt = 0; attempt < 32; ++attempt)
        {
            if (mode == 0)
            {
                center[0] = static_cast<double>(io.rawInteger(-4, 4));
                center[1] = static_cast<double>(io.rawInteger(-4, 4));
                radius = 5.0 * io.rawInteger(1, 3);
                int32_t a = io.rawInteger(0, 11);
                int32_t span = io.rawInteger(1, 5);
                int32_t b = (a + span) % 12;
                d0 = { kCircleLattice[a][0] / 5.0, kCircleLattice[a][1] / 5.0 };
                d1 = { kCircleLattice[b][0] / 5.0, kCircleLattice[b][1] / 5.0 };
            }
            else
            {
                center[0] = io.raw(-4.0, 4.0);
                center[1] = io.raw(-4.0, 4.0);
                radius = io.raw(radiusLo, radiusHi);
                double a0 = io.raw(0.0, 6.2831853071795862);
                double a1 = a0 + io.raw(0.1, 3.0);
                d0 = { std::cos(a0), std::sin(a0) };
                d1 = { std::cos(a1), std::sin(a1) };
            }
            if (DotPerp(d0, d1) > 0.0) { break; }
        }
        arc.center = io.givenVec(center);
        arc.radius = io.given(radius);
        arc.end[0] = io.givenVec(Vector2<double>{
            center[0] + radius * d0[0], center[1] + radius * d0[1] });
        arc.end[1] = io.givenVec(Vector2<double>{
            center[0] + radius * d1[0], center[1] + radius * d1[1] });
        return arc;
    }

    // Evaluate a NURBS curve at a drawn parameter and emit the jet.
    void EmitCurveJet(oracle::Ctx& io, NURBSCurve<2, double> const& curve,
        double t, int32_t order)
    {
        std::array<Vector<2, double>, 4> jet{};
        curve.Evaluate(t, static_cast<uint32_t>(order), jet.data());
        io.outReal(curve.GetTMin());
        io.outReal(curve.GetTMax());
        for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
    }
}

ORACLE_CASE("NURBSCircle.quarterDegree2")
{
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    NURBSQuarterCircleDegree2<double> curve;
    EmitCurveJet(io, curve, t, order);
}

ORACLE_CASE("NURBSCircle.quarterDegree4")
{
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    NURBSQuarterCircleDegree4<double> curve;
    EmitCurveJet(io, curve, t, order);
}

ORACLE_CASE("NURBSCircle.halfDegree3")
{
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    NURBSHalfCircleDegree3<double> curve;
    EmitCurveJet(io, curve, t, order);
}

ORACLE_CASE("NURBSCircle.fullDegree3")
{
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    NURBSFullCircleDegree3<double> curve;
    EmitCurveJet(io, curve, t, order);
}

ORACLE_CASE("NURBSCircle.arcDegree2")
{
    auto arc = MakeArc(io, 0.5, 4.0);
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    NURBSCircularArcDegree2<double> curve(arc);
    io.outReal(curve.GetWeight(0));
    io.outReal(curve.GetWeight(1));
    io.outReal(curve.GetWeight(2));
    for (int32_t i = 0; i < 3; ++i) { io.outVec(curve.GetControl(i)); }
    EmitCurveJet(io, curve, t, order);
}

// ============================================================== NURBSSphere

namespace
{
    // A point of the triangular domain {u >= 0, v >= 0, u + v <= 1} of
    // NURBSEighthSphereDegree4. The modes visit the three corners, the three
    // edges, dyadic interior points and uniform interior points.
    void DrawUV(oracle::Ctx& io, double& u, double& v)
    {
        int32_t mode = io.rawInteger(0, 7);
        switch (mode)
        {
        case 0: u = 0.0; v = 0.0; break;
        case 1: u = 1.0; v = 0.0; break;
        case 2: u = 0.0; v = 1.0; break;
        case 3: u = 0.0; v = io.rawInteger(1, 7) * 0.125; break;
        case 4: u = io.rawInteger(1, 7) * 0.125; v = 0.0; break;
        case 5:
        {
            double a = io.rawInteger(1, 7) * 0.125;
            u = a; v = 1.0 - a;
            break;
        }
        case 6:
        {
            int32_t a = io.rawInteger(1, 6);
            int32_t b = io.rawInteger(1, 7 - a);
            u = a * 0.125; v = b * 0.125;
            break;
        }
        default:
        {
            double a = io.raw(0.0, 1.0);
            double b = io.raw(0.0, 1.0 - a);
            u = a; v = b;
            break;
        }
        }
    }
}

ORACLE_CASE("NURBSSphere.eighthDegree4")
{
    int32_t maxOrder = io.integer(0, 2);
    double u = 0.0, v = 0.0;
    DrawUV(io, u, v);
    io.given(u);
    io.given(v);

    NURBSEighthSphereDegree4<double> patch;
    std::array<Vector<3, double>, 6> values{};
    patch.Evaluate(u, v, static_cast<uint32_t>(maxOrder), values.data());
    int32_t const numValues = (maxOrder == 0 ? 1 : (maxOrder == 1 ? 3 : 6));
    for (int32_t i = 0; i < numValues; ++i) { io.outVec(values[i]); }
}

ORACLE_CASE("NURBSSphere.halfDegree3")
{
    int32_t order = io.integer(0, 2);
    double u = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    double v = io.given(DrawParameter(io, 0.0, 1.0, nullptr));

    NURBSHalfSphereDegree3<double> surface;
    std::array<Vector<3, double>, 6> jet{};
    surface.Evaluate(u, v, static_cast<uint32_t>(order), jet.data());
    int32_t const numValues = (order == 0 ? 1 : (order == 1 ? 3 : 6));
    for (int32_t i = 0; i < numValues; ++i) { io.outVec(jet[i]); }
}

ORACLE_CASE("NURBSSphere.fullDegree3")
{
    int32_t order = io.integer(0, 2);
    double u = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    double v = io.given(DrawParameter(io, 0.0, 1.0, nullptr));

    NURBSFullSphereDegree3<double> surface;
    std::array<Vector<3, double>, 6> jet{};
    surface.Evaluate(u, v, static_cast<uint32_t>(order), jet.data());
    int32_t const numValues = (order == 0 ? 1 : (order == 1 ? 3 : 6));
    for (int32_t i = 0; i < numValues; ++i) { io.outVec(jet[i]); }
}

// ============================================================= DarbouxFrame
// The surface is a NURBS sphere from this group, so the frame and the
// principal information have an independent reference: the frame is
// orthonormal and the principal curvatures of a unit sphere are +-1.

ORACLE_CASE("DarbouxFrame.compute")
{
    int32_t which = io.integer(0, 1);
    double u = io.real(0.05, 0.95);
    double v = io.real(0.05, 0.95);

    std::shared_ptr<ParametricSurface<3, double>> surface;
    if (which == 0)
    {
        surface = std::make_shared<NURBSHalfSphereDegree3<double>>();
    }
    else
    {
        surface = std::make_shared<NURBSFullSphereDegree3<double>>();
    }
    DarbouxFrame3<double> frame(surface);
    Vector3<double> position{}, tangent0{}, tangent1{}, normal{};
    frame(u, v, position, tangent0, tangent1, normal);
    io.outVec(position);
    io.outVec(tangent0);
    io.outVec(tangent1);
    io.outVec(normal);
}

ORACLE_CASE("DarbouxFrame.getPrincipalInformation")
{
    int32_t which = io.integer(0, 1);
    double u = io.real(0.05, 0.95);
    double v = io.real(0.05, 0.95);

    std::shared_ptr<ParametricSurface<3, double>> surface;
    if (which == 0)
    {
        surface = std::make_shared<NURBSHalfSphereDegree3<double>>();
    }
    else
    {
        surface = std::make_shared<NURBSFullSphereDegree3<double>>();
    }
    DarbouxFrame3<double> frame(surface);
    double curvature0 = 0.0, curvature1 = 0.0;
    Vector3<double> direction0{}, direction1{};
    frame.GetPrincipalInformation(u, v, curvature0, curvature1, direction0,
        direction1);
    io.outReal(curvature0);
    io.outReal(curvature1);
    io.outVec(direction0);
    io.outVec(direction1);
}

// =========================================================== TCBSplineCurve

namespace
{
    struct TCBInput
    {
        std::vector<Vector<3, double>> point{};
        std::vector<double> time{}, tension{}, continuity{}, bias{}, lambda{};
        Vector<3, double> firstOut{}, lastIn{};
        bool hasFirstOut = false, hasLastIn = false;
    };

    // Draws and records: count, controls, times, tension, continuity, bias,
    // the lambda flag and (when set) the lambdas, then the two optional
    // boundary tangents preceded by their flags.
    TCBInput MakeTCBInput(oracle::Ctx& io, bool forceUniform)
    {
        TCBInput in{};
        int32_t n = io.integer(2, 6);
        int32_t controlMode = (forceUniform ? 0 : io.rawInteger(0, 2));
        in.point = RecordControls3(io, DrawControls3(io, n, controlMode));
        in.time = MakeTimes(io, n);
        in.tension.resize(static_cast<size_t>(n));
        in.continuity.resize(static_cast<size_t>(n));
        in.bias.resize(static_cast<size_t>(n));
        for (int32_t i = 0; i < n; ++i) { in.tension[i] = io.real(-1.0, 1.0); }
        for (int32_t i = 0; i < n; ++i) { in.continuity[i] = io.real(-1.0, 1.0); }
        for (int32_t i = 0; i < n; ++i) { in.bias[i] = io.real(-1.0, 1.0); }
        int32_t useLambda = io.integer(0, 1);
        if (useLambda == 1)
        {
            in.lambda.resize(static_cast<size_t>(n));
            for (int32_t i = 0; i < n; ++i) { in.lambda[i] = io.real(0.25, 2.0); }
        }
        int32_t useFirst = io.integer(0, 1);
        if (useFirst == 1)
        {
            in.firstOut = io.latticeVec<3>(-3, 3);
            in.hasFirstOut = true;
        }
        int32_t useLast = io.integer(0, 1);
        if (useLast == 1)
        {
            in.lastIn = io.latticeVec<3>(-3, 3);
            in.hasLastIn = true;
        }
        return in;
    }
}

ORACLE_CASE("TCBSplineCurve.evaluate")
{
    // Lambda scaling divides by the sum of the two Kochanek-Bartels tangent
    // lengths, which is zero only for coincident neighbours; the uniform
    // control-point mode is forced when lambda is used so that the documented
    // division by zero (issue #415, preserved) stays out of this case. See
    // TCBSplineCurve.degenerateLambda.
    auto in = MakeTCBInput(io, false);
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, in.time.front(), in.time.back(),
        nullptr));

    TCBSplineCurve<3, double> curve(in.point, in.time, in.tension,
        in.continuity, in.bias, in.lambda,
        in.hasFirstOut ? &in.firstOut : nullptr,
        in.hasLastIn ? &in.lastIn : nullptr);
    io.outInt(static_cast<int32_t>(curve.GetNumKeyFrames()));
    for (auto const& v : curve.GetInTangents()) { io.outVec(v); }
    for (auto const& v : curve.GetOutTangents()) { io.outVec(v); }
    std::array<Vector<3, double>, 4> jet{};
    curve.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

// Constructed degenerate: coincident neighbours plus lambda make both
// Kochanek-Bartels tangents vanish, so 'common' is 0/0 and the tangents and
// two segments' coefficients become NaN (docs/UPSTREAM-FINDINGS.md issue
// #415, preserved). The NaN pattern is emitted as booleans, because the
// harness treats any NaN as equal to any NaN.
ORACLE_CASE("TCBSplineCurve.degenerateLambda")
{
    int32_t n = io.integer(3, 5);
    int32_t k = io.integer(1, n - 2);
    std::vector<Vector<3, double>> point(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i)
    {
        for (int32_t j = 0; j < 3; ++j)
        {
            point[i][j] = static_cast<double>(io.rawInteger(-3, 3));
        }
    }
    // Make P[k-1] == P[k+1] so that both tangents at key k vanish.
    point[k + 1] = point[k - 1];
    for (int32_t i = 0; i < n; ++i) { io.givenVec(point[i]); }
    std::vector<double> time(static_cast<size_t>(n));
    for (int32_t i = 0; i < n; ++i) { time[i] = static_cast<double>(i); }
    std::vector<double> zeros(static_cast<size_t>(n), 0.0);
    std::vector<double> lambda(static_cast<size_t>(n), 1.0);
    double t = io.real(0.0, static_cast<double>(n - 1));

    TCBSplineCurve<3, double> curve(point, time, zeros, zeros, zeros, lambda,
        nullptr, nullptr);
    for (auto const& v : curve.GetInTangents())
    {
        for (int32_t j = 0; j < 3; ++j) { io.outBool(std::isnan(v[j])); }
    }
    for (auto const& v : curve.GetOutTangents())
    {
        for (int32_t j = 0; j < 3; ++j) { io.outBool(std::isnan(v[j])); }
    }
    std::array<Vector<3, double>, 4> jet{};
    curve.Evaluate(t, 3, jet.data());
    for (int32_t i = 0; i < 4; ++i)
    {
        for (int32_t j = 0; j < 3; ++j) { io.outBool(std::isnan(jet[i][j])); }
    }
}

// The port sets mConstructed; upstream never does (issue #182, fixed in the
// port). Every record deviates.
ORACLE_CASE("TCBSplineCurve.isConstructed")
{
    auto in = MakeTCBInput(io, true);
    TCBSplineCurve<3, double> curve(in.point, in.time, in.tension,
        in.continuity, in.bias, in.lambda,
        in.hasFirstOut ? &in.firstOut : nullptr,
        in.hasLastIn ? &in.lastIn : nullptr);
    io.outBool(static_cast<bool>(curve));
}

// ======================================== NaturalCubicSpline / QuinticSpline

namespace
{
    // Interpolation points for a natural spline. 'closed' forces the last
    // point to equal the first, which is what a closed spline means.
    std::vector<Vector<3, double>> MakeSplinePoints(oracle::Ctx& io,
        int32_t count, bool closed)
    {
        int32_t mode = io.rawInteger(0, 3);
        auto P = DrawControls3(io, count, mode);
        if (closed) { P[static_cast<size_t>(count) - 1] = P[0]; }
        return RecordControls3(io, P);
    }
}

ORACLE_CASE("NaturalCubicSpline.free")
{
    int32_t n = io.integer(3, 6);
    auto f0 = MakeSplinePoints(io, n, false);
    auto times = MakeTimes(io, n);
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, times.front(), times.back(), nullptr));

    NaturalCubicSpline<3, double> spline(true, n, f0.data(), times.data());
    io.outInt(spline.GetNumSegments());
    for (auto const& poly : spline.GetPolynomials())
    {
        for (auto const& c : poly) { io.outVec(c); }
    }
    std::array<Vector<3, double>, 4> jet{};
    spline.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

ORACLE_CASE("NaturalCubicSpline.closed")
{
    // n = 3 is the minimum point count (issue #295 covers the sibling
    // NaturalSplineCurve; this class solves the closed system differently).
    int32_t n = io.integer(3, 6);
    auto f0 = MakeSplinePoints(io, n, true);
    auto times = MakeTimes(io, n);
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, times.front(), times.back(), nullptr));

    NaturalCubicSpline<3, double> spline(false, n, f0.data(), times.data());
    io.outInt(spline.GetNumSegments());
    for (auto const& poly : spline.GetPolynomials())
    {
        for (auto const& c : poly) { io.outVec(c); }
    }
    std::array<Vector<3, double>, 4> jet{};
    spline.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

ORACLE_CASE("NaturalCubicSpline.clamped")
{
    int32_t n = io.integer(3, 6);
    auto f0 = MakeSplinePoints(io, n, false);
    auto times = MakeTimes(io, n);
    auto d0 = io.latticeVec<3>(-4, 4);
    auto d1 = io.latticeVec<3>(-4, 4);
    int32_t order = io.integer(0, 3);
    double t = io.given(DrawParameter(io, times.front(), times.back(), nullptr));

    NaturalCubicSpline<3, double> spline(n, f0.data(), times.data(), d0, d1);
    io.outInt(spline.GetNumSegments());
    for (auto const& poly : spline.GetPolynomials())
    {
        for (auto const& c : poly) { io.outVec(c); }
    }
    std::array<Vector<3, double>, 4> jet{};
    spline.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

// The multiple-segment paths of ParametricCurve::GetLength (a partial first
// segment, whole interior segments and a partial last segment) and of
// GetTime, reached through a natural cubic spline.
ORACLE_CASE("ParametricCurve.getLength.multiSegment")
{
    int32_t n = io.integer(3, 5);
    auto f0 = MakeSplinePoints(io, n, false);
    auto times = MakeTimes(io, n);
    int32_t rombergOrder = io.integer(2, 4);
    int32_t maxBisections = io.integer(4, 12);
    double a = io.real(0.0, 1.0);
    double b = io.real(0.0, 1.0);

    NaturalCubicSpline<3, double> spline(true, n, f0.data(), times.data());
    spline.SetRombergOrder(rombergOrder);
    spline.SetMaxBisections(static_cast<uint32_t>(maxBisections));
    double tmin = spline.GetTMin(), tmax = spline.GetTMax();
    double t0 = tmin + (tmax - tmin) * a;
    double t1 = tmin + (tmax - tmin) * b;
    io.outReal(spline.GetLength(t0, t1));
    io.outReal(spline.GetTotalLength());
    io.outReal(spline.GetTime(0.5 * spline.GetTotalLength()));
    // Exact segment boundaries.
    io.outReal(spline.GetLength(times.front(), times[1]));
    io.outReal(spline.GetLength(times[1], times.back()));
}

namespace
{
    std::vector<Vector<3, double>> MakeSplineDerivatives(oracle::Ctx& io,
        int32_t count)
    {
        int32_t mode = io.rawInteger(0, 1);
        std::vector<Vector<3, double>> D(static_cast<size_t>(count));
        for (int32_t i = 0; i < count; ++i)
        {
            for (int32_t j = 0; j < 3; ++j)
            {
                D[i][j] = (mode == 0
                    ? static_cast<double>(io.rawInteger(-3, 3))
                    : io.raw(-3.0, 3.0));
            }
        }
        return RecordControls3(io, D);
    }
}

ORACLE_CASE("NaturalQuinticSpline.free")
{
    int32_t n = io.integer(3, 6);
    auto f0 = MakeSplinePoints(io, n, false);
    auto f1 = MakeSplineDerivatives(io, n);
    auto times = MakeTimes(io, n);
    int32_t order = io.integer(0, 5);
    double t = io.given(DrawParameter(io, times.front(), times.back(), nullptr));

    NaturalQuinticSpline<3, double> spline(true, n, f0.data(), f1.data(),
        times.data());
    io.outInt(spline.GetNumSegments());
    for (auto const& poly : spline.GetPolynomials())
    {
        for (auto const& c : poly) { io.outVec(c); }
    }
    std::array<Vector<3, double>, 6> jet{};
    spline.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

ORACLE_CASE("NaturalQuinticSpline.closed")
{
    int32_t n = io.integer(3, 6);
    auto f0 = MakeSplinePoints(io, n, true);
    auto f1 = MakeSplineDerivatives(io, n);
    auto times = MakeTimes(io, n);
    int32_t order = io.integer(0, 5);
    double t = io.given(DrawParameter(io, times.front(), times.back(), nullptr));

    NaturalQuinticSpline<3, double> spline(false, n, f0.data(), f1.data(),
        times.data());
    io.outInt(spline.GetNumSegments());
    for (auto const& poly : spline.GetPolynomials())
    {
        for (auto const& c : poly) { io.outVec(c); }
    }
    std::array<Vector<3, double>, 6> jet{};
    spline.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

ORACLE_CASE("NaturalQuinticSpline.clamped")
{
    int32_t n = io.integer(3, 6);
    auto f0 = MakeSplinePoints(io, n, false);
    auto f1 = MakeSplineDerivatives(io, n);
    auto times = MakeTimes(io, n);
    auto d0 = io.latticeVec<3>(-4, 4);
    auto d1 = io.latticeVec<3>(-4, 4);
    int32_t order = io.integer(0, 5);
    double t = io.given(DrawParameter(io, times.front(), times.back(), nullptr));

    NaturalQuinticSpline<3, double> spline(n, f0.data(), f1.data(),
        times.data(), d0, d1);
    io.outInt(spline.GetNumSegments());
    for (auto const& poly : spline.GetPolynomials())
    {
        for (auto const& c : poly) { io.outVec(c); }
    }
    std::array<Vector<3, double>, 6> jet{};
    spline.Evaluate(t, static_cast<uint32_t>(order), jet.data());
    for (int32_t i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

// ============================================================= FrenetFrame

namespace
{
    std::vector<Vector<2, double>> MakeControls2(oracle::Ctx& io, int32_t count)
    {
        int32_t mode = io.rawInteger(0, 3);
        std::vector<Vector<2, double>> P(static_cast<size_t>(count));
        if (mode == 3)
        {
            // All control points coincident: the curve has zero speed
            // everywhere, which is the degenerate branch of GetCurvature and
            // of the tangent normalization.
            Vector<2, double> q{ static_cast<double>(io.rawInteger(-3, 3)),
                static_cast<double>(io.rawInteger(-3, 3)) };
            for (int32_t i = 0; i < count; ++i) { P[i] = q; }
        }
        else
        {
            for (int32_t i = 0; i < count; ++i)
            {
                for (int32_t j = 0; j < 2; ++j)
                {
                    P[i][j] = (mode == 0 ? io.raw(-5.0, 5.0)
                        : static_cast<double>(io.rawInteger(-4, 4)));
                }
            }
        }
        for (auto const& p : P) { io.givenVec(p); }
        return P;
    }
}

ORACLE_CASE("FrenetFrame2.compute")
{
    int32_t degree = io.integer(2, 5);
    auto controls = MakeControls2(io, degree + 1);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));

    auto curve = std::make_shared<BezierCurve<2, double>>(degree, controls.data());
    FrenetFrame2<double> frame(curve);
    Vector2<double> position{}, tangent{}, normal{};
    frame(t, position, tangent, normal);
    io.outVec(position);
    io.outVec(tangent);
    io.outVec(normal);
}

// GetCurvature calls std::pow(speedSqr, 1.5); the case carries a tolerance.
ORACLE_CASE("FrenetFrame2.getCurvature")
{
    int32_t degree = io.integer(2, 5);
    auto controls = MakeControls2(io, degree + 1);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));

    auto curve = std::make_shared<BezierCurve<2, double>>(degree, controls.data());
    FrenetFrame2<double> frame(curve);
    io.outReal(frame.GetCurvature(t));
}

ORACLE_CASE("FrenetFrame3.compute")
{
    int32_t degree = io.integer(2, 5);
    auto controls = MakeControls3(io, degree + 1);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));

    auto curve = std::make_shared<BezierCurve<3, double>>(degree, controls.data());
    FrenetFrame3<double> frame(curve);
    Vector3<double> position{}, tangent{}, normal{}, binormal{};
    frame(t, position, tangent, normal, binormal);
    io.outVec(position);
    io.outVec(tangent);
    io.outVec(normal);
    io.outVec(binormal);
}

// GetCurvature calls std::pow(speedSqr, 1.5); GetTorsion is arithmetic only
// and is therefore required to be bit-identical (outRealExact on the
// TypeScript side).
ORACLE_CASE("FrenetFrame3.getCurvatureAndTorsion")
{
    int32_t degree = io.integer(3, 5);
    auto controls = MakeControls3(io, degree + 1);
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));

    auto curve = std::make_shared<BezierCurve<3, double>>(degree, controls.data());
    FrenetFrame3<double> frame(curve);
    io.outReal(frame.GetCurvature(t));
    io.outReal(frame.GetTorsion(t));
}

// =========================================================== ImplicitCurve2

namespace
{
    // A general conic F(x,y) = c0 + c1*x + c2*y + c3*x^2 + c4*x*y + c5*y^2.
    // The TypeScript replay defines the same subclass with the same
    // expression order.
    class ConicCurve2 : public ImplicitCurve2<double>
    {
    public:
        explicit ConicCurve2(std::array<double, 6> const& c) : mC(c) {}

        virtual double F(Vector2<double> const& p) const override
        {
            return mC[0] + mC[1] * p[0] + mC[2] * p[1] + mC[3] * p[0] * p[0]
                + mC[4] * p[0] * p[1] + mC[5] * p[1] * p[1];
        }

        virtual double FX(Vector2<double> const& p) const override
        {
            return mC[1] + 2.0 * mC[3] * p[0] + mC[4] * p[1];
        }

        virtual double FY(Vector2<double> const& p) const override
        {
            return mC[2] + mC[4] * p[0] + 2.0 * mC[5] * p[1];
        }

        virtual double FXX(Vector2<double> const&) const override
        {
            return 2.0 * mC[3];
        }

        virtual double FXY(Vector2<double> const&) const override
        {
            return mC[4];
        }

        virtual double FYY(Vector2<double> const&) const override
        {
            return 2.0 * mC[5];
        }

    private:
        std::array<double, 6> mC;
    };

    std::array<double, 6> MakeConic(oracle::Ctx& io)
    {
        std::array<double, 6> c{};
        int32_t mode = io.rawInteger(0, 2);
        for (int32_t i = 0; i < 6; ++i)
        {
            if (mode == 0) { c[i] = io.raw(-3.0, 3.0); }
            else if (mode == 1) { c[i] = static_cast<double>(io.rawInteger(-3, 3)); }
            else { c[i] = (i >= 3 ? static_cast<double>(io.rawInteger(0, 2)) : 0.0); }
        }
        for (int32_t i = 0; i < 6; ++i) { io.given(c[i]); }
        return c;
    }
}

ORACLE_CASE("ImplicitCurve2.evaluate")
{
    auto c = MakeConic(io);
    auto p = io.vec<2>(-4.0, 4.0);
    double epsilon = io.real(0.0, 1.0);

    ConicCurve2 curve(c);
    io.outReal(curve.F(p));
    io.outReal(curve.FX(p));
    io.outReal(curve.FY(p));
    io.outReal(curve.FXX(p));
    io.outReal(curve.FXY(p));
    io.outReal(curve.FYY(p));
    io.outBool(curve.IsOnCurve(p, epsilon));
    io.outVec(curve.GetGradient(p));
    io.outMat(curve.GetHessian(p));
    Vector2<double> tangent{}, normal{};
    curve.GetFrame(p, tangent, normal);
    io.outVec(tangent);
    io.outVec(normal);
}

// GetCurvature calls std::pow(fx^2 + fy^2, 1.5); the case carries a
// tolerance. The boolean is exact: it is decided by 'denom == 0', and
// std::pow and Math.pow both return exactly 0 only for a zero base.
ORACLE_CASE("ImplicitCurve2.getCurvature")
{
    auto c = MakeConic(io);
    auto p = io.vec<2>(-4.0, 4.0);

    ConicCurve2 curve(c);
    double curvature = 0.0;
    bool valid = curve.GetCurvature(p, curvature);
    io.outBool(valid);
    io.outReal(curvature);
}

// ================================================ ReparameterizeByArclength

ORACLE_CASE("ReparameterizeByArclength.getT")
{
    int32_t degree = io.integer(2, 4);
    auto controls = MakeControls3(io, degree + 1);
    int32_t rombergOrder = io.integer(2, 4);
    int32_t useBisection = io.integer(0, 1);
    double fraction = io.real(-0.2, 1.2);

    auto curve = std::make_shared<BezierCurve<3, double>>(degree, controls.data());
    curve->SetRombergOrder(rombergOrder);
    ReparameterizeByArclength<3, double> reparam(curve);
    io.outReal(reparam.GetTMin());
    io.outReal(reparam.GetTMax());
    io.outReal(reparam.GetTotalArclength());
    double s = fraction * reparam.GetTotalArclength();
    auto output = reparam.GetT(s, useBisection == 1);
    io.outReal(output.t);
    io.outReal(output.f);
    io.outInt(static_cast<int32_t>(output.numIterations));
}

// =========================================================== SampleCircularArc

namespace
{
    // The number of samples is static_cast<size_t>(radius * angle / k) with
    // angle from std::acos, so it is decided by the C math library. MSVC and
    // V8 may differ by one ulp in acos, which would change the truncation
    // only when radius*angle/k sits at an integer boundary; the probe below
    // rejects those draws (fractional part outside [0.15, 0.85]) and bounds
    // the sample count. Everything else in SampleCircularArc is + - * / and
    // sqrt, so the case is exact.
    bool ProbeArc(Vector2<double> const& center, double radius,
        Vector2<double> const& e0, Vector2<double> const& e1, bool wantLarge)
    {
        Vector2<double> P0 = (e0 - center) / radius;
        Vector2<double> P2 = (e1 - center) / radius;
        double dot = std::max(-1.0, std::min(Dot(P0, P2), 1.0));
        double angle = std::acos(dot);
        double dotPerp = DotPerp(P0, P2);
        double value = 0.0;
        if (dotPerp >= 0.0)
        {
            if (wantLarge) { return false; }
            value = (dot >= 0.0 ? radius * angle : radius * angle / 2.0);
        }
        else
        {
            if (!wantLarge) { return false; }
            double a = static_cast<double>(GTE_C_TWO_PI) - angle;
            value = (dot <= 0.0 ? radius * a / 3.0 : radius * a / 4.0);
        }
        double frac = value - std::floor(value);
        return frac >= 0.15 && frac <= 0.85 && value >= 1.0 && value <= 8.0;
    }

    Arc2<double> MakeSampledArc(oracle::Ctx& io, bool wantLarge)
    {
        Vector2<double> center{ 0.0, 0.0 }, e0{}, e1{};
        double radius = 4.0;
        double lo = (wantLarge ? 3.2 : 0.2);
        double hi = (wantLarge ? 6.2 : 3.0);
        bool accepted = false;
        for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
        {
            center[0] = io.raw(-3.0, 3.0);
            center[1] = io.raw(-3.0, 3.0);
            radius = io.raw(1.0, 8.0);
            double a0 = io.raw(0.0, static_cast<double>(GTE_C_TWO_PI));
            double a1 = a0 + io.raw(lo, hi);
            e0 = { center[0] + radius * std::cos(a0),
                center[1] + radius * std::sin(a0) };
            e1 = { center[0] + radius * std::cos(a1),
                center[1] + radius * std::sin(a1) };
            accepted = ProbeArc(center, radius, e0, e1, wantLarge);
        }
        if (!accepted)
        {
            // Deterministic fallback: a quarter circle (or its complement)
            // of radius 4, for which radius*angle/k is 6.283... .
            center = { 0.0, 0.0 };
            radius = 4.0;
            e0 = { 4.0, 0.0 };
            e1 = (wantLarge ? Vector2<double>{ 0.0, -4.0 }
                : Vector2<double>{ 0.0, 4.0 });
        }

        Arc2<double> arc{};
        arc.center = io.givenVec(center);
        arc.radius = io.given(radius);
        arc.end[0] = io.givenVec(e0);
        arc.end[1] = io.givenVec(e1);
        return arc;
    }
}

ORACLE_CASE("SampleCircularArc.compute")
{
    auto arc = MakeSampledArc(io, false);
    SampleCircularArc<double> sampler;
    std::vector<Vector2<double>> points;
    sampler(arc, points);
    io.outInt(static_cast<int32_t>(points.size()));
    for (auto const& p : points) { io.outVec(p); }
}

// Arcs larger than pi reach SampleArc3/SampleArc4, whose split directions
// are exact only for the midpoint construction (docs/UPSTREAM-FINDINGS.md,
// issue #183, preserved): the samples lie on the circle but are not
// angularly ordered. The port reproduces upstream bit for bit.
ORACLE_CASE("SampleCircularArc.compute.largeArc")
{
    auto arc = MakeSampledArc(io, true);
    SampleCircularArc<double> sampler;
    std::vector<Vector2<double>> points;
    sampler(arc, points);
    io.outInt(static_cast<int32_t>(points.size()));
    for (auto const& p : points) { io.outVec(p); }
}

// ============================================================ PolylineOffset

namespace
{
    // Reject polylines with an exact direction reversal at a vertex: the
    // offset there is (d / (1 + Dot(N0, N1))) * (N0 + N1) = (d/0) * 0 = NaN
    // on both sides, which tests nothing.
    bool ProbePolyline(std::vector<Vector2<double>> const& V, bool isOpen)
    {
        size_t const n = V.size();
        if (n < 3) { return true; }
        std::vector<Vector2<double>> N(isOpen ? n - 1 : n);
        for (size_t i0 = 0, i1 = 1; i1 < n; i0 = i1++)
        {
            Vector2<double> D = V[i1] - V[i0];
            Normalize(D);
            N[i0] = Perp(D);
        }
        if (!isOpen)
        {
            Vector2<double> D = V[0] - V[n - 1];
            Normalize(D);
            N[n - 1] = Perp(D);
        }
        for (size_t i0 = 0, i1 = 1; i1 < N.size(); i0 = i1++)
        {
            if (1.0 + Dot(N[i0], N[i1]) == 0.0) { return false; }
        }
        if (!isOpen)
        {
            if (1.0 + Dot(N[N.size() - 1], N[0]) == 0.0) { return false; }
            if (1.0 + Dot(N[N.size() - 2], N[N.size() - 1]) == 0.0) { return false; }
        }
        return true;
    }

    // Modes: 0 uniform, 1 small lattice (parallel and collinear runs occur),
    // 2 convex lattice polygon (all turns the same way), 3 a lattice polyline
    // with an explicitly collinear run.
    std::vector<Vector2<double>> MakePolyline2(oracle::Ctx& io, int32_t n,
        bool isOpen)
    {
        std::vector<Vector2<double>> V(static_cast<size_t>(n));
        int32_t mode = io.rawInteger(0, 3);
        bool accepted = false;
        for (int32_t attempt = 0; attempt < 48 && !accepted; ++attempt)
        {
            if (mode == 0)
            {
                for (int32_t i = 0; i < n; ++i)
                {
                    V[i][0] = io.raw(-6.0, 6.0);
                    V[i][1] = io.raw(-6.0, 6.0);
                }
            }
            else if (mode == 1)
            {
                for (int32_t i = 0; i < n; ++i)
                {
                    V[i][0] = static_cast<double>(io.rawInteger(-4, 4));
                    V[i][1] = static_cast<double>(io.rawInteger(-4, 4));
                }
            }
            else if (mode == 2)
            {
                // A convex polygon: the twelve lattice directions of the
                // radius-5 circle, taken in increasing order.
                int32_t start = io.rawInteger(0, 11);
                for (int32_t i = 0; i < n; ++i)
                {
                    int32_t k = (start + i) % 12;
                    V[i][0] = static_cast<double>(kCircleLattice[k][0]);
                    V[i][1] = static_cast<double>(kCircleLattice[k][1]);
                }
            }
            else
            {
                Vector2<double> base{
                    static_cast<double>(io.rawInteger(-3, 3)),
                    static_cast<double>(io.rawInteger(-3, 3)) };
                Vector2<double> dir{
                    static_cast<double>(io.rawInteger(-2, 2)),
                    static_cast<double>(io.rawInteger(1, 3)) };
                for (int32_t i = 0; i < n; ++i)
                {
                    if (i < n / 2)
                    {
                        V[i][0] = base[0] + i * dir[0];
                        V[i][1] = base[1] + i * dir[1];
                    }
                    else
                    {
                        V[i][0] = static_cast<double>(io.rawInteger(-5, 5));
                        V[i][1] = static_cast<double>(io.rawInteger(-5, 5));
                    }
                }
            }
            accepted = ProbePolyline(V, isOpen);
        }
        if (!accepted)
        {
            for (int32_t i = 0; i < n; ++i)
            {
                V[i][0] = static_cast<double>(kCircleLattice[i % 12][0]);
                V[i][1] = static_cast<double>(kCircleLattice[i % 12][1]);
            }
        }
        for (auto const& v : V) { io.givenVec(v); }
        return V;
    }
}

ORACLE_CASE("PolylineOffset.execute")
{
    int32_t isOpen = io.integer(0, 1);
    int32_t n = io.integer(isOpen == 1 ? 2 : 3, 8);
    auto vertices = MakePolyline2(io, n, isOpen == 1);
    double distance = io.real(0.125, 3.0);
    int32_t which = io.integer(0, 2);

    PolylineOffset<double> offset(vertices, isOpen == 1);
    std::vector<Vector2<double>> right{}, left{};
    bool offsetRight = (which != 1);
    bool offsetLeft = (which != 0);
    offset.Execute(distance, offsetRight, right, offsetLeft, left);
    io.outInt(static_cast<int32_t>(right.size()));
    for (auto const& p : right) { io.outVec(p); }
    io.outInt(static_cast<int32_t>(left.size()));
    for (auto const& p : left) { io.outVec(p); }
}

// Execute rejects a nonpositive distance and a request for neither polyline;
// the constructor rejects too few vertices.
ORACLE_CASE("PolylineOffset.throwParity")
{
    int32_t which = io.integer(0, 2);
    int32_t isOpen = io.integer(0, 1);
    int32_t n = (which == 2 ? io.integer(1, 2) : io.integer(3, 5));
    auto vertices = MakePolyline2(io, n, isOpen == 1);
    double distance = (which == 0 ? io.real(-2.0, 0.0) : io.real(0.5, 2.0));

    PolylineOffset<double> offset(vertices, isOpen == 1);
    std::vector<Vector2<double>> right{}, left{};
    offset.Execute(distance, which != 1, right, false, left);
    io.outInt(static_cast<int32_t>(right.size()));
    for (auto const& p : right) { io.outVec(p); }
}

// ============================================================= CLODPolyline
// Only closed polylines and the two-vertex open polyline are generated.
// For an open polyline with three or more vertices upstream's ComputeEdges
// builds the surviving edge from collapses[0] + 1 == numVertices and
// ReorderVertices then reads permute[numVertices], one past the end of the
// vector: undefined behaviour, so those inputs cannot be compared against the
// real build at all (docs/UPSTREAM-FINDINGS.md, issue #182; the port fixes it
// and test/CLODPolyline.test.ts pins the fixed behaviour).

ORACLE_CASE("CLODPolyline.construct")
{
    int32_t closed = io.integer(0, 1);
    int32_t n = (closed == 1 ? io.integer(3, 8) : 2);
    auto vertices = MakePolyline2(io, n, closed == 0);

    CLODPolyline<2, double> polyline(vertices, closed == 1);
    io.outInt(polyline.GetNumVertices());
    io.outBool(polyline.GetClosed());
    io.outInt(polyline.GetMinLevelOfDetail());
    io.outInt(polyline.GetMaxLevelOfDetail());
    io.outInt(polyline.GetLevelOfDetail());
    io.outInt(polyline.GetNumEdges());
    for (auto const& v : polyline.GetVertices()) { io.outVec(v); }
    auto const& edges = polyline.GetEdges();
    for (int32_t i = 0; i < 2 * polyline.GetNumEdges(); ++i)
    {
        io.outInt(edges[i]);
    }
}

ORACLE_CASE("CLODPolyline.setLevelOfDetail")
{
    int32_t n = io.integer(4, 8);
    auto vertices = MakePolyline2(io, n, false);
    int32_t lod = io.integer(1, n + 1);
    int32_t lod2 = io.integer(1, n + 1);

    CLODPolyline<2, double> polyline(vertices, true);
    polyline.SetLevelOfDetail(lod);
    io.outInt(polyline.GetLevelOfDetail());
    io.outInt(polyline.GetNumEdges());
    auto const& edges = polyline.GetEdges();
    for (int32_t i = 0; i < 2 * polyline.GetNumEdges(); ++i)
    {
        io.outInt(edges[i]);
    }
    polyline.SetLevelOfDetail(lod2);
    io.outInt(polyline.GetLevelOfDetail());
    io.outInt(polyline.GetNumEdges());
    for (int32_t i = 0; i < 2 * polyline.GetNumEdges(); ++i)
    {
        io.outInt(edges[i]);
    }
}

// =========================================================== MassSpringCurve

ORACLE_CASE("MassSpringCurve.update")
{
    int32_t numParticles = io.integer(3, 5);
    double step = io.real(0.01, 0.2);
    int32_t numSteps = io.integer(1, 3);
    int32_t immovable = io.integer(-1, numParticles - 1);

    std::vector<double> mass(static_cast<size_t>(numParticles));
    for (int32_t i = 0; i < numParticles; ++i) { mass[i] = io.real(0.5, 3.0); }
    std::vector<Vector<3, double>> position(static_cast<size_t>(numParticles));
    for (int32_t i = 0; i < numParticles; ++i)
    {
        // Separated along x so that no two particles coincide (a zero-length
        // spring divides by zero in Acceleration).
        double dx = io.real(0.5, 1.5);
        double py = io.real(-1.0, 1.0);
        double pz = io.real(-1.0, 1.0);
        position[i][0] = static_cast<double>(i) + dx;
        position[i][1] = py;
        position[i][2] = pz;
    }
    std::vector<Vector<3, double>> velocity(static_cast<size_t>(numParticles));
    for (int32_t i = 0; i < numParticles; ++i)
    {
        velocity[i] = io.vec<3>(-1.0, 1.0);
    }
    int32_t numSprings = numParticles - 1;
    std::vector<double> constant(static_cast<size_t>(numSprings));
    for (int32_t i = 0; i < numSprings; ++i) { constant[i] = io.real(0.5, 4.0); }
    std::vector<double> restLength(static_cast<size_t>(numSprings));
    for (int32_t i = 0; i < numSprings; ++i) { restLength[i] = io.real(0.25, 2.0); }
    double time = io.real(0.0, 1.0);

    MassSpringCurve<3, double> system(numParticles, step);
    for (int32_t i = 0; i < numParticles; ++i)
    {
        system.SetMass(i, i == immovable ? std::numeric_limits<double>::max()
            : mass[i]);
        system.SetPosition(i, position[i]);
        system.SetVelocity(i, velocity[i]);
    }
    for (int32_t i = 0; i < numSprings; ++i)
    {
        system.SetConstant(i, constant[i]);
        system.SetLength(i, restLength[i]);
    }

    io.outInt(system.GetNumParticles());
    io.outInt(system.GetNumSprings());
    io.outReal(system.GetStep());
    double t = time;
    for (int32_t s = 0; s < numSteps; ++s)
    {
        system.Update(t);
        t += step;
    }
    for (int32_t i = 0; i < numParticles; ++i)
    {
        io.outReal(system.GetMass(i));
        io.outVec(system.GetPosition(i));
        io.outVec(system.GetVelocity(i));
    }
    for (int32_t i = 0; i < numSprings; ++i)
    {
        io.outReal(system.GetConstant(i));
        io.outReal(system.GetLength(i));
    }
}
