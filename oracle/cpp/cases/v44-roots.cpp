// Group 44 (category roots-polynomials): CubicRootsQR.h, QuarticRootsQR.h,
// Polynomial1.h, PolynomialCurve.h, PolynomialRoot.h, RootsBisection.h,
// RootsBisection1.h, RootsBisection2.h, RootsBrentsMethod.h, RootsLinear.h,
// RootsQuadratic.h, RootsPolynomial.h, RootsGeneralPolynomial.h.
//
// Everything here is arithmetic only (+ - * /, comparisons, std::sqrt,
// std::fabs, std::max) except RootsPolynomial's low-degree closed forms, which
// call std::pow, std::atan2, std::cos and std::sin. The arithmetic cases are
// compared bit-for-bit; the closed-form ones carry a libm tolerance.
//
// Instantiations. The port has no templates, so each case instantiates the
// C++ side the way the port does:
//   * RootsLinear / RootsQuadratic expose both the T = double and the
//     T = BSRational<UIntegerAP32> instantiations, and the classification is
//     exact rational in both.
//   * RootsPolynomial is templated on a separate 'Rational' type for the
//     classification; the port instantiates it with 'number', so the C++ side
//     is called as RootsPolynomial<double>::SolveQuadratic<double>(...).
//   * RootsBisection1/RootsBisection2 are instantiated with T = double only
//     (the arbitrary-precision constructors are not ported).
//
// Deliberate port fixes demonstrated by 'deviation' cases:
//   * RootsBisection1 with maxIterations == 1 leaves its outputs unwritten
//     (issue #84); RootsBisection2 then reports a previous call's values.
//   * PolynomialCurve never sets mConstructed, so operator bool is false for
//     every valid curve (issue #319).
//   * RootsGeneralPolynomial runs on a zero-padded rational polynomial when
//     the input has high-order zero coefficients (issue #340, item 4).
#define ORACLE_FAMILY "v44-roots"
#include "Oracle.h"

#include <Mathematics/CubicRootsQR.h>
#include <Mathematics/Polynomial1.h>
#include <Mathematics/PolynomialCurve.h>
#include <Mathematics/PolynomialRoot.h>
#include <Mathematics/QuarticRootsQR.h>
#include <Mathematics/RootsBisection.h>
#include <Mathematics/RootsBisection1.h>
#include <Mathematics/RootsBisection2.h>
#include <Mathematics/RootsBrentsMethod.h>
#include <Mathematics/RootsGeneralPolynomial.h>
#include <Mathematics/RootsLinear.h>
#include <Mathematics/RootsPolynomial.h>
#include <Mathematics/RootsQuadratic.h>

#include <array>
#include <cmath>
#include <map>
#include <vector>

using namespace gte;

namespace
{
    using Rational = BSRational<UIntegerAP32>;

    // ---- generators -------------------------------------------------------

    void ExpandRoots(double lead, std::vector<double> const& r, std::vector<double>& c)
    {
        c.assign(r.size() + 1, 0.0);
        c[0] = 1.0;
        size_t degree = 0;
        for (double root : r)
        {
            ++degree;
            c[degree] = c[degree - 1];
            for (size_t i = degree - 1; i >= 1; --i)
            {
                c[i] = c[i - 1] - root * c[i];
            }
            c[0] = -root * c[0];
        }
        for (double& v : c) { v *= lead; }
    }

    double NonzeroRawInteger(oracle::Ctx& io, int lo, int hi)
    {
        int v = 0;
        do { v = io.rawInteger(lo, hi); } while (v == 0);
        return static_cast<double>(v);
    }

    // Eight generator modes, cycled on io.index() % 8, each producing
    // degree+1 coefficients (index 0 is the constant term). 'monic' forces the
    // leading coefficient to 1.
    void DrawCoefficients(oracle::Ctx& io, int degree, bool monic, int mode,
        std::vector<double>& c)
    {
        size_t const n = static_cast<size_t>(degree) + 1;
        size_t const free = n - (monic ? 1 : 0);
        c.assign(n, 0.0);
        std::vector<double> roots;
        switch (mode)
        {
        case 0:
            // Small integer lattice.
            for (size_t i = 0; i < free; ++i)
            {
                c[i] = static_cast<double>(io.rawInteger(-5, 5));
            }
            break;
        case 1:
            // Uniform reals.
            for (size_t i = 0; i < free; ++i) { c[i] = io.raw(-10.0, 10.0); }
            break;
        case 2:
            // Distinct integer roots.
            for (int i = 0; i < degree; ++i)
            {
                roots.push_back(static_cast<double>(io.rawInteger(-6, 6)));
            }
            ExpandRoots(monic ? 1.0 : NonzeroRawInteger(io, -4, 4), roots, c);
            break;
        case 3:
        {
            // Integer roots with repeated values: double, triple and
            // quadruple roots, and two double roots.
            double a = static_cast<double>(io.rawInteger(-4, 4));
            double b = static_cast<double>(io.rawInteger(-4, 4));
            for (int i = 0; i < degree; ++i)
            {
                roots.push_back(i % 2 == 0 ? a : (io.rawInteger(0, 1) == 0 ? a : b));
            }
            ExpandRoots(monic ? 1.0 : NonzeroRawInteger(io, -4, 4), roots, c);
            break;
        }
        case 4:
        {
            // A common power-of-ten scale.
            double s = std::pow(10.0, static_cast<double>(io.rawInteger(-100, 100)));
            for (size_t i = 0; i < free; ++i) { c[i] = io.raw(-10.0, 10.0) * s; }
            if (!monic) { c[n - 1] = io.raw(-10.0, 10.0) * s; }
            break;
        }
        case 5:
            // Per-coefficient scales.
            for (size_t i = 0; i < free; ++i)
            {
                c[i] = io.raw(-10.0, 10.0)
                    * std::pow(10.0, static_cast<double>(io.rawInteger(-30, 30)));
            }
            break;
        case 6:
            // Leading coefficient zero (degree drop), sometimes two.
            for (size_t i = 0; i < free; ++i)
            {
                c[i] = static_cast<double>(io.rawInteger(-5, 5));
            }
            if (!monic)
            {
                c[n - 1] = 0.0;
                if (io.rawInteger(0, 1) == 0 && n >= 3) { c[n - 2] = 0.0; }
            }
            else if (free >= 2)
            {
                c[free - 1] = 0.0;
            }
            break;
        case 7:
        default:
            // Zero-valued roots, including a signed zero constant term.
            for (size_t i = 0; i < free; ++i) { c[i] = io.raw(-10.0, 10.0); }
            if (!monic) { c[n - 1] = io.raw(-10.0, 10.0); }
            c[0] = (io.rawInteger(0, 1) == 0 ? 0.0 : -0.0);
            if (free >= 2 && io.rawInteger(0, 1) == 0) { c[1] = 0.0; }
            if (free >= 3 && io.rawInteger(0, 2) == 0) { c[2] = 0.0; }
            break;
        }
        if (monic) { c[n - 1] = 1.0; }
    }

    // Draw and record degree+1 (or 'degree' for monic) coefficients.
    void DrawAndRecord(oracle::Ctx& io, int degree, bool monic, double* c)
    {
        std::vector<double> v;
        DrawCoefficients(io, degree, monic, io.index() % 8, v);
        int count = degree + (monic ? 0 : 1);
        for (int i = 0; i < count; ++i) { c[i] = io.given(v[i]); }
    }

    // A Polynomial1 whose degree is itself a recorded input, constructed the
    // way the port does: the degree-argument constructor (which zero-fills and
    // does NOT eliminate leading zeros) followed by per-coefficient
    // assignment, so the degree-drop mode keeps its leading zeros.
    // 'nonzeroLeading' is used for divisors, whose leading coefficient
    // upstream inverts.
    Polynomial1<double> DrawPoly(oracle::Ctx& io, int maxDegree, bool nonzeroLeading)
    {
        int degree = io.integer(0, maxDegree);
        std::vector<double> v;
        DrawCoefficients(io, degree, false, io.index() % 8, v);
        if (nonzeroLeading && v[static_cast<size_t>(degree)] == 0.0)
        {
            v[static_cast<size_t>(degree)] = 1.0;
        }
        Polynomial1<double> p(static_cast<uint32_t>(degree));
        for (int i = 0; i <= degree; ++i) { p[i] = io.given(v[i]); }
        return p;
    }

    void EmitPoly(oracle::Ctx& io, Polynomial1<double> const& p)
    {
        uint32_t degree = p.GetDegree();
        io.outInt(degree);
        for (uint32_t i = 0; i <= degree; ++i) { io.outReal(p[i]); }
    }
}

// ---------------------------------------------------------------------------
// Polynomial1.h
// ---------------------------------------------------------------------------

// Every binary and unary operation, plus the comparison operators, which
// std::vector implements lexicographically.
ORACLE_CASE("Polynomial1.arithmetic")
{
    auto p0 = DrawPoly(io, 4, false);
    auto p1 = DrawPoly(io, 4, false);
    double s = io.real(-5.0, 5.0);
    double nonzero = NonzeroRawInteger(io, -6, 6);
    io.given(nonzero);
    EmitPoly(io, +p0);
    EmitPoly(io, -p0);
    EmitPoly(io, p0 + p1);
    EmitPoly(io, p0 - p1);
    EmitPoly(io, p0 * p1);
    EmitPoly(io, p0 + s);
    EmitPoly(io, s + p0);
    EmitPoly(io, p0 - s);
    EmitPoly(io, s - p0);
    EmitPoly(io, p0 * s);
    EmitPoly(io, s * p0);
    EmitPoly(io, p0 / nonzero);
    io.outBool(p0 == p1);
    io.outBool(p0 != p1);
    io.outBool(p0 < p1);
    io.outBool(p0 <= p1);
    io.outBool(p0 > p1);
    io.outBool(p0 >= p1);
}

// Horner evaluation, the derivative chain, the inversion and the translation.
ORACLE_CASE("Polynomial1.evaluateAndTransforms")
{
    auto p = DrawPoly(io, 5, false);
    double t0 = io.real(-3.0, 3.0);
    double t1 = io.lattice(-4, 4);
    io.outReal(p(t0));
    io.outReal(p(t1));
    auto d1 = p.GetDerivative();
    EmitPoly(io, d1);
    EmitPoly(io, d1.GetDerivative());
    EmitPoly(io, p.GetInversion());
    EmitPoly(io, p.GetTranslation(t0));
    EmitPoly(io, p.GetTranslation(t1));
}

// The in-place mutators. EliminateLeadingZeros matters for the degree-drop
// generator mode; SetDegree grows with zeros and shrinks by truncation.
ORACLE_CASE("Polynomial1.mutators")
{
    auto p = DrawPoly(io, 5, false);
    int newDegree = io.integer(0, 6);
    auto q = p;
    q.EliminateLeadingZeros();
    EmitPoly(io, q);
    auto r = p;
    r.MakeMonic();
    EmitPoly(io, r);
    auto t = p;
    t.SetDegree(static_cast<uint32_t>(newDegree));
    EmitPoly(io, t);
    auto u = p;
    u.SetCoefficients(2.5);
    EmitPoly(io, u);
}

// Euclidean division, including the degree(P) < degree(D) path.
ORACLE_CASE("Polynomial1.divide")
{
    auto p = DrawPoly(io, 5, false);
    auto d = DrawPoly(io, 3, true);
    Polynomial1<double> quotient{}, remainder{};
    p.Divide(d, quotient, remainder);
    EmitPoly(io, quotient);
    EmitPoly(io, remainder);
}

ORACLE_CASE("Polynomial1.greatestCommonDivisor")
{
    auto p0 = DrawPoly(io, 4, false);
    auto p1 = DrawPoly(io, 4, false);
    EmitPoly(io, GreatestCommonDivisor(p0, p1));
    EmitPoly(io, GreatestCommonDivisor(p1, p0));
}

namespace
{
    // A verbatim copy of SquareFreeFactorization with an iteration cap. The
    // upstream do-while loop is unbounded and does not terminate whenever the
    // floating-point GreatestCommonDivisor fails to reduce the degree of b
    // (issue #83): for f = (t-2)^2 (t+2)(t+3)(t-4) it spins forever. The port
    // caps the iterations at degree(f) + 1, the exact-arithmetic bound, and
    // throws. The generator uses this probe with the same cap, so the main
    // case covers only the inputs on which upstream terminates within the
    // port's bound; there is no way to record a non-terminating upstream run.
    bool SquareFreeTerminates(Polynomial1<double> const& f)
    {
        uint32_t const cap = f.GetDegree() + 1;
        Polynomial1<double> fder = f.GetDerivative();
        Polynomial1<double> a, b, c, d, q, r;
        a = GreatestCommonDivisor(f, fder);
        f.Divide(a, b, r);
        fder.Divide(a, c, r);
        d = c - b.GetDerivative();
        uint32_t iteration = 0;
        do
        {
            if (++iteration > cap) { return false; }
            a = GreatestCommonDivisor(b, d);
            b.Divide(a, q, r);
            b = q;
            d.Divide(a, c, r);
            d = c - b.GetDerivative();
        } while (b.GetDegree() > 0);
        return true;
    }
}

ORACLE_CASE("Polynomial1.squareFreeFactorization")
{
    // Rejection: only inputs on which upstream's unbounded loop terminates
    // within the port's iteration cap. The loop is capped and redraws the
    // whole polynomial; the fallback is the monic linear polynomial t, whose
    // factorization terminates in one iteration.
    int degree = io.integer(1, 4);
    std::vector<double> v;
    bool ok = false;
    Polynomial1<double> f(static_cast<uint32_t>(degree));
    for (int attempt = 0; attempt < 64 && !ok; ++attempt)
    {
        DrawCoefficients(io, degree, true, io.index() % 8, v);
        for (int i = 0; i <= degree; ++i) { f[i] = v[i]; }
        ok = SquareFreeTerminates(f);
    }
    if (!ok)
    {
        v.assign(static_cast<size_t>(degree) + 1, 0.0);
        v[1] = 1.0;
        v[static_cast<size_t>(degree)] = 1.0;
    }
    for (int i = 0; i <= degree; ++i) { f[i] = io.given(v[i]); }
    std::vector<Polynomial1<double>> factors;
    SquareFreeFactorization(f, factors);
    io.outInt(factors.size());
    for (auto const& factor : factors) { EmitPoly(io, factor); }
}

// ---------------------------------------------------------------------------
// PolynomialCurve.h
// ---------------------------------------------------------------------------

ORACLE_CASE("PolynomialCurve.evaluate")
{
    double tmin = io.real(-2.0, 0.0);
    double tmax = io.real(0.5, 3.0);
    std::array<Polynomial1<double>, 3> components{};
    components[0] = DrawPoly(io, 4, false);
    components[1] = DrawPoly(io, 4, false);
    components[2] = DrawPoly(io, 4, false);
    double t = io.real(-2.0, 3.0);
    uint32_t order = static_cast<uint32_t>(io.integer(0, 3));
    PolynomialCurve<3, double> curve(tmin, tmax, components);
    for (int i = 0; i < 3; ++i)
    {
        EmitPoly(io, curve.GetPolynomial(i));
        EmitPoly(io, curve.GetDer1Polynomial(i));
        EmitPoly(io, curve.GetDer2Polynomial(i));
        EmitPoly(io, curve.GetDer3Polynomial(i));
    }
    std::array<Vector<3, double>, 4> jet{};
    curve.Evaluate(t, order, jet.data());
    for (uint32_t k = 0; k <= order; ++k) { io.outVec(jet[k]); }
}

// Deviation: upstream never sets mConstructed, so ParametricCurve's
// operator bool reports failure for every valid polynomial curve. The port
// sets it. UPSTREAM-FINDINGS PolynomialCurve.h, issue #319.
ORACLE_CASE("PolynomialCurve.deviation.constructed")
{
    double tmin = io.real(-2.0, 0.0);
    double tmax = io.real(0.5, 3.0);
    std::array<Polynomial1<double>, 3> components{};
    components[0] = DrawPoly(io, 2, false);
    components[1] = DrawPoly(io, 2, false);
    components[2] = DrawPoly(io, 2, false);
    PolynomialCurve<3, double> curve(tmin, tmax, components);
    io.outBool(static_cast<bool>(curve));
}

// ---------------------------------------------------------------------------
// PolynomialRoot.h, RootsBisection.h, RootsBisection1.h, RootsBisection2.h,
// RootsBrentsMethod.h
//
// The target functions of every bisection and Brent case are polynomials
// evaluated by Horner from the recorded coefficients, written identically on
// both sides, so they use only multiplication and addition.
// ---------------------------------------------------------------------------

namespace
{
    // Horner, identical on both sides: c[n] is the leading coefficient.
    double HornerAt(std::vector<double> const& c, double t)
    {
        size_t i = c.size() - 1;
        double result = c[i];
        while (i-- > 0) { result = t * result + c[i]; }
        return result;
    }

    // A quintic with integer roots in [-3,3], recorded as 6 coefficients, so
    // that the bracketing endpoints can be placed exactly on a root, exactly
    // between two roots, or outside every root.
    std::vector<double> DrawBracketPoly(oracle::Ctx& io, int degree)
    {
        std::vector<double> v;
        DrawCoefficients(io, degree, false, io.index() % 8, v);
        if (v[static_cast<size_t>(degree)] == 0.0) { v[static_cast<size_t>(degree)] = 1.0; }
        for (size_t i = 0; i < v.size(); ++i) { v[i] = io.given(v[i]); }
        return v;
    }
}

ORACLE_CASE("PolynomialRoot.bisect")
{
    auto c = DrawBracketPoly(io, 3);
    double xMin = io.lattice(-6, 0);
    double xMax = io.lattice(1, 6);
    int signFMin = io.integer(0, 1) == 0 ? -1 : +1;
    int signFMax = -signFMin;
    auto F = [&c](double x) { return HornerAt(c, x); };
    double a = xMin, b = xMax;
    PolynomialRootBisect<double>(F, signFMin, signFMax, a, b);
    io.outReal(a);
    io.outReal(b);
}

// The comparison operators of PolynomialRoot compare only the root estimate.
ORACLE_CASE("PolynomialRoot.compare")
{
    double x0 = io.lattice(-3, 3);
    double x1 = io.lattice(-3, 3);
    int m0 = io.integer(1, 4);
    int m1 = io.integer(1, 4);
    PolynomialRoot<double> r0(x0, static_cast<size_t>(m0));
    PolynomialRoot<double> r1(x1, static_cast<size_t>(m1));
    io.outBool(r0 == r1);
    io.outBool(r0 < r1);
    io.outBool(r1 < r0);
    PolynomialRoot<double> def{};
    io.outReal(def.x);
    io.outInt(def.m);
}

ORACLE_CASE("RootsBisection.find")
{
    // The endpoint ranges overlap, so t1 <= t0 occurs and the invalid-interval
    // return 0 is covered along with the f0*f1 > 0 rejection.
    auto c = DrawBracketPoly(io, 4);
    double t0 = io.lattice(-6, 1);
    double t1 = io.lattice(-1, 6);
    int maxIterations = io.integer(0, 60);
    auto F = [&c](double t) { return HornerAt(c, t); };
    double root = 0.0;
    uint32_t iterations = RootsBisection<double>::Find(F, t0, t1,
        static_cast<uint32_t>(maxIterations), root);
    io.outInt(iterations);
    io.outReal(root);
}

// The overload that takes f0 = F(t0) and f1 = F(t1). Passing the signs rather
// than the values is the documented use, so the case draws both forms.
ORACLE_CASE("RootsBisection.findWithValues")
{
    auto c = DrawBracketPoly(io, 4);
    double t0 = io.lattice(-6, 1);
    double t1 = io.lattice(-1, 6);
    int maxIterations = io.integer(0, 60);
    bool useSigns = io.boolean();
    auto F = [&c](double t) { return HornerAt(c, t); };
    double f0 = F(t0), f1 = F(t1);
    if (useSigns)
    {
        f0 = (f0 > 0.0 ? 1.0 : (f0 < 0.0 ? -1.0 : 0.0));
        f1 = (f1 > 0.0 ? 1.0 : (f1 < 0.0 ? -1.0 : 0.0));
    }
    double root = 0.0;
    uint32_t iterations = RootsBisection<double>::Find(F, t0, t1, f0, f1,
        static_cast<uint32_t>(maxIterations), root);
    io.outInt(iterations);
    io.outReal(root);
}

ORACLE_CASE("RootsBisection1.find")
{
    auto c = DrawBracketPoly(io, 4);
    double tMin = io.lattice(-6, 0);
    double tMax = io.lattice(1, 6);
    // maxIterations >= 2: with 1 the outputs are left unwritten upstream,
    // which the deviation case below covers.
    int maxIterations = io.integer(2, 60);
    bool passValues = io.boolean();
    auto F = [&c](double const& t) { return HornerAt(c, t); };
    RootsBisection1<double> bisector(static_cast<uint32_t>(maxIterations));
    double tRoot = 0.0, fAtTRoot = 0.0;
    uint32_t iterations;
    if (passValues)
    {
        double fMin = F(tMin), fMax = F(tMax);
        iterations = bisector(F, tMin, tMax, fMin, fMax, tRoot, fAtTRoot);
    }
    else
    {
        iterations = bisector(F, tMin, tMax, tRoot, fAtTRoot);
    }
    io.outInt(iterations);
    io.outReal(tRoot);
    io.outReal(fAtTRoot);
}

// Deviation: with maxIterations == 1 the bisection loop never runs and
// upstream returns tRoot and fAtTRoot unassigned, so the caller reads back
// whatever it passed in. The port returns explicit zeros. The case seeds the
// output variables with a sentinel, which is what a caller reusing a variable
// sees. UPSTREAM-FINDINGS RootsBisection1.h item 1, issue #84.
ORACLE_CASE("RootsBisection1.deviation.maxIterationsOne")
{
    auto c = DrawBracketPoly(io, 3);
    double tMin = io.lattice(-6, -1);
    double tMax = io.lattice(1, 6);
    double sentinel = io.real(1.0, 2.0);
    auto F = [&c](double const& t) { return HornerAt(c, t); };
    RootsBisection1<double> bisector(1);
    double tRoot = sentinel, fAtTRoot = sentinel;
    uint32_t iterations = bisector(F, tMin, tMax, tRoot, fAtTRoot);
    io.outInt(iterations);
    io.outReal(tRoot);
    io.outReal(fAtTRoot);
}

// RootsBisection2 on the rectangle [-4,4]x[-4,4] with
//   F(x,y) = x^3 + a*x + b*y + c,  G(x,y) = y^3 + d*y + e*x + f
// and |a|..|f| <= 5. Then |a*x + b*y + c| <= 45 < 64 at every corner, so
// F(-4,y) < 0 < F(4,y) and G(x,-4) < 0 < G(x,4) hold for every draw and the
// preconditions of the bisector are satisfied without rejection.
ORACLE_CASE("RootsBisection2.find")
{
    double a = io.real(-5.0, 5.0);
    double b = io.real(-5.0, 5.0);
    double c = io.real(-5.0, 5.0);
    double d = io.real(-5.0, 5.0);
    double e = io.real(-5.0, 5.0);
    double f = io.real(-5.0, 5.0);
    int xMaxIterations = io.integer(2, 40);
    int yMaxIterations = io.integer(2, 40);
    auto F = [&a, &b, &c](double const& x, double const& y)
        { return x * x * x + a * x + b * y + c; };
    auto G = [&d, &e, &f](double const& x, double const& y)
        { return y * y * y + d * y + e * x + f; };
    RootsBisection2<double> bisector(static_cast<uint32_t>(xMaxIterations),
        static_cast<uint32_t>(yMaxIterations));
    double xRoot = 0.0, yRoot = 0.0, fAtRoot = 0.0, gAtRoot = 0.0;
    uint32_t iterations = bisector(F, G, -4.0, 4.0, -4.0, 4.0, xRoot, yRoot,
        fAtRoot, gAtRoot);
    io.outInt(iterations);
    io.outReal(xRoot);
    io.outReal(yRoot);
    io.outReal(fAtRoot);
    io.outReal(gAtRoot);
    io.outBool(bisector.NoGuaranteeForRootBound());
}

// Deviation: with xMaxIterations == 1 the x-bisector never writes its outputs
// (see RootsBisection1 item 1), so upstream's members keep the values of the
// PREVIOUS call while the return value claims a successful bisection. The
// first call has F(xMin, y) exactly zero, which makes RootsBisection1 write
// mXRoot = xMin; the second call then reads that stale xMin where the port
// reports 0. UPSTREAM-FINDINGS RootsBisection2.h item 2, issues #84, #152.
ORACLE_CASE("RootsBisection2.deviation.staleOutputs")
{
    double a = io.lattice(-5, 5);
    double c2 = io.lattice(-5, 5);
    double d = io.lattice(-5, 5);
    double e = io.lattice(-5, 5);
    double f = io.lattice(-5, 5);
    // c1 makes F1(-4, y) = 0 exactly (integer arithmetic).
    double c1 = -(-64.0 + a * -4.0);
    io.given(c1);
    auto G = [&d, &e, &f](double const& x, double const& y)
        { return y * y * y + d * y + e * x + f; };
    RootsBisection2<double> bisector(1, 20);
    double xRoot = 0.0, yRoot = 0.0, fAtRoot = 0.0, gAtRoot = 0.0;

    auto F1 = [&a, &c1](double const& x, double const&)
        { return x * x * x + a * x + c1; };
    uint32_t it1 = bisector(F1, G, -4.0, 4.0, -4.0, 4.0, xRoot, yRoot, fAtRoot, gAtRoot);
    io.outInt(it1);
    io.outReal(xRoot);
    io.outReal(fAtRoot);

    auto F2 = [&a, &c2](double const& x, double const&)
        { return x * x * x + a * x + c2; };
    uint32_t it2 = bisector(F2, G, -4.0, 4.0, -4.0, 4.0, xRoot, yRoot, fAtRoot, gAtRoot);
    io.outInt(it2);
    io.outReal(xRoot);
    io.outReal(fAtRoot);
}

// Brent's method. The parameter draws reach every validation rejection
// (t1 <= t0, maxIterations 0, a positive negFTolerance, a negative
// posFTolerance / stepTTolerance / convTTolerance) and every termination
// criterion (an endpoint inside the function tolerance, the consecutive-float
// bisection exit, the function tolerance at an interpolated point, and the
// subinterval tolerance).
ORACLE_CASE("RootsBrentsMethod.find")
{
    auto c = DrawBracketPoly(io, 4);
    double t0 = io.lattice(-6, 1);
    double t1 = io.lattice(-1, 6);
    int maxIterations = io.integer(0, 80);
    int mode = io.integer(0, 5);
    double scale = io.real(0.0, 1.0);
    double negFTolerance = 0.0, posFTolerance = 0.0;
    double stepTTolerance = 0.0, convTTolerance = 0.0;
    switch (mode)
    {
    case 0: break;                                        // all tolerances 0
    case 1: negFTolerance = -scale; posFTolerance = scale; break;
    case 2: stepTTolerance = scale; break;
    case 3: convTTolerance = scale; break;
    case 4: negFTolerance = scale; break;                 // invalid (> 0)
    default: posFTolerance = -scale; break;               // invalid (< 0)
    }
    auto F = [&c](double t) { return HornerAt(c, t); };
    double root = 0.0;
    bool found = RootsBrentsMethod<double>::Find(F, t0, t1,
        static_cast<uint32_t>(maxIterations), negFTolerance, posFTolerance,
        stepTTolerance, convTTolerance, root);
    io.outBool(found);
    io.outReal(root);
}

// Throw parity: RootsBisection1 asserts tMin < tMax.
ORACLE_CASE("RootsBisection1.throwParity")
{
    auto c = DrawBracketPoly(io, 2);
    double tMin = io.lattice(-3, 3);
    double tMax = io.lattice(-3, 3);
    auto F = [&c](double const& t) { return HornerAt(c, t); };
    RootsBisection1<double> bisector(20);
    double tRoot = 0.0, fAtTRoot = 0.0;
    uint32_t iterations = bisector(F, tMin, tMax, tRoot, fAtTRoot);
    io.outInt(iterations);
    io.outReal(tRoot);
    io.outReal(fAtTRoot);
}

// ---------------------------------------------------------------------------
// CubicRootsQR.h, QuarticRootsQR.h
// ---------------------------------------------------------------------------

ORACLE_CASE("CubicRootsQR.solve")
{
    double c[3]{};
    DrawAndRecord(io, 3, true, c);
    int maxIterations = io.integer(0, 64);
    uint32_t numRoots = 0;
    std::array<double, 3> roots{};
    uint32_t iterations = CubicRootsQR<double>()(static_cast<uint32_t>(maxIterations),
        c[0], c[1], c[2], numRoots, roots);
    io.outInt(iterations);
    io.outInt(numRoots);
    for (double r : roots) { io.outReal(r); }
}

ORACLE_CASE("CubicRootsQR.solveMatrix")
{
    CubicRootsQR<double>::Matrix A{};
    for (int r = 0; r < 3; ++r)
    {
        for (int cc = 0; cc < 3; ++cc) { A[r][cc] = io.real(-4.0, 4.0); }
    }
    int maxIterations = io.integer(0, 64);
    uint32_t numRoots = 0;
    std::array<double, 3> roots{};
    uint32_t iterations = CubicRootsQR<double>()(static_cast<uint32_t>(maxIterations),
        A, numRoots, roots);
    io.outInt(iterations);
    io.outInt(numRoots);
    for (double r : roots) { io.outReal(r); }
    for (int r = 0; r < 3; ++r)
    {
        for (int cc = 0; cc < 3; ++cc) { io.outReal(A[r][cc]); }
    }
}

ORACLE_CASE("QuarticRootsQR.solve")
{
    double c[4]{};
    DrawAndRecord(io, 4, true, c);
    int maxIterations = io.integer(0, 64);
    uint32_t numRoots = 0;
    std::array<double, 4> roots{};
    uint32_t iterations = QuarticRootsQR<double>()(static_cast<uint32_t>(maxIterations),
        c[0], c[1], c[2], c[3], numRoots, roots);
    io.outInt(iterations);
    io.outInt(numRoots);
    for (double r : roots) { io.outReal(r); }
}

ORACLE_CASE("QuarticRootsQR.solveMatrix")
{
    QuarticRootsQR<double>::Matrix A{};
    for (int r = 0; r < 4; ++r)
    {
        for (int cc = 0; cc < 4; ++cc) { A[r][cc] = io.real(-4.0, 4.0); }
    }
    int maxIterations = io.integer(0, 64);
    uint32_t numRoots = 0;
    std::array<double, 4> roots{};
    uint32_t iterations = QuarticRootsQR<double>()(static_cast<uint32_t>(maxIterations),
        A, numRoots, roots);
    io.outInt(iterations);
    io.outInt(numRoots);
    for (double r : roots) { io.outReal(r); }
    for (int r = 0; r < 4; ++r)
    {
        for (int cc = 0; cc < 4; ++cc) { io.outReal(A[r][cc]); }
    }
}

// ---------------------------------------------------------------------------
// RootsLinear.h, RootsQuadratic.h
// ---------------------------------------------------------------------------

namespace
{
    void EmitRoots(oracle::Ctx& io, PolynomialRoot<double> const* roots, size_t numRoots)
    {
        io.outInt(numRoots);
        for (size_t i = 0; i < numRoots; ++i)
        {
            io.outReal(roots[i].x);
            io.outInt(roots[i].m);
        }
    }

    void EmitRationalRoots(oracle::Ctx& io, PolynomialRoot<Rational> const* roots,
        size_t numRoots)
    {
        io.outInt(numRoots);
        for (size_t i = 0; i < numRoots; ++i)
        {
            io.outReal(static_cast<double>(roots[i].x));
            io.outInt(roots[i].m);
        }
    }
}

ORACLE_CASE("RootsLinear.solve")
{
    double c[2]{};
    DrawAndRecord(io, 1, false, c);
    std::array<PolynomialRoot<double>, 1> roots{};
    size_t numRoots = RootsLinear<double>::Solve(c[0], c[1], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsLinear.solveMonic")
{
    double c[1]{};
    DrawAndRecord(io, 1, true, c);
    std::array<PolynomialRoot<double>, 1> roots{};
    size_t numRoots = RootsLinear<double>::Solve(c[0], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsLinear.solveRational")
{
    double c[2]{};
    DrawAndRecord(io, 1, false, c);
    std::array<PolynomialRoot<Rational>, 1> gRoots{};
    size_t numRoots = RootsLinear<Rational>::Solve(Rational(c[0]), Rational(c[1]),
        gRoots.data());
    EmitRationalRoots(io, gRoots.data(), numRoots);
    std::array<PolynomialRoot<Rational>, 1> mRoots{};
    size_t numMonic = RootsLinear<Rational>::Solve(Rational(c[0]), mRoots.data());
    EmitRationalRoots(io, mRoots.data(), numMonic);
}

ORACLE_CASE("RootsQuadratic.solve.bisection")
{
    double c[3]{};
    DrawAndRecord(io, 2, false, c);
    std::array<PolynomialRoot<double>, 2> roots{};
    size_t numRoots = RootsQuadratic<double>::Solve(true, c[0], c[1], c[2], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

// The closed form uses std::sqrt only, which is correctly rounded in both
// runtimes, so this case is exact too.
ORACLE_CASE("RootsQuadratic.solve.closedForm")
{
    double c[3]{};
    DrawAndRecord(io, 2, false, c);
    std::array<PolynomialRoot<double>, 2> roots{};
    size_t numRoots = RootsQuadratic<double>::Solve(false, c[0], c[1], c[2], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuadratic.solveMonic")
{
    double c[2]{};
    DrawAndRecord(io, 2, true, c);
    bool useBisection = io.boolean();
    std::array<PolynomialRoot<double>, 2> roots{};
    size_t numRoots = RootsQuadratic<double>::Solve(useBisection, c[0], c[1], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuadratic.solveDepressed")
{
    double d0 = io.real(-10.0, 10.0);
    bool useBisection = io.boolean();
    std::array<PolynomialRoot<double>, 2> roots{};
    size_t numRoots = RootsQuadratic<double>::Solve(useBisection, d0, roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuadratic.computeDepressedRoots")
{
    double d0 = io.lattice(-8, 8);
    bool useBisection = io.boolean();
    std::array<PolynomialRoot<Rational>, 2> roots{};
    size_t numRoots = RootsQuadratic<double>::ComputeDepressedRoots(useBisection,
        Rational(d0), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuadratic.solveRational")
{
    double c[3]{};
    DrawAndRecord(io, 2, false, c);
    bool useBisection = io.boolean();
    std::array<PolynomialRoot<Rational>, 2> gRoots{};
    size_t numRoots = RootsQuadratic<Rational>::Solve(useBisection, Rational(c[0]),
        Rational(c[1]), Rational(c[2]), gRoots.data());
    EmitRationalRoots(io, gRoots.data(), numRoots);
    std::array<PolynomialRoot<Rational>, 2> mRoots{};
    size_t numMonic = RootsQuadratic<Rational>::Solve(useBisection, Rational(c[0]),
        Rational(c[1]), mRoots.data());
    EmitRationalRoots(io, mRoots.data(), numMonic);
    std::array<PolynomialRoot<Rational>, 2> dRoots{};
    size_t numDepressed = RootsQuadratic<Rational>::Solve(useBisection, Rational(c[0]),
        dRoots.data());
    EmitRationalRoots(io, dRoots.data(), numDepressed);
}

// ---------------------------------------------------------------------------
// RootsPolynomial.h
//
// The port instantiates the classification type with 'number', so the C++
// side is RootsPolynomial<double> with the template parameter Rational
// deduced as double. The low-degree Solve* functions call std::pow,
// std::atan2, std::cos and std::sin, so their root VALUES carry the libm
// tolerance; the GetRootInfo* functions and Find use only arithmetic and
// std::sqrt and are exact.
// ---------------------------------------------------------------------------

namespace
{
    void EmitRootMap(oracle::Ctx& io, std::map<double, int32_t> const& rmMap)
    {
        io.outInt(rmMap.size());
        for (auto const& rm : rmMap)
        {
            io.outReal(rm.first);
            io.outInt(rm.second);
        }
    }

    // The leading coefficient must be nonzero for the Solve*/GetRootInfo*
    // preconditions.
    void DrawNonzeroLeading(oracle::Ctx& io, int degree, double* c)
    {
        std::vector<double> v;
        DrawCoefficients(io, degree, false, io.index() % 8, v);
        if (v[static_cast<size_t>(degree)] == 0.0) { v[static_cast<size_t>(degree)] = 1.0; }
        for (int i = 0; i <= degree; ++i) { c[i] = io.given(v[i]); }
    }

    // Two rejections for the Solve* cases.
    //
    // 1. A non-finite root. The classifiers are computed in 'double' here (the
    //    port's instantiation), so a cancelling denominator such as
    //    9*c1^2 - 2*c2*a1 can round to exactly zero and produce an infinite or
    //    NaN root. Inserting a NaN key into std::map is undefined behaviour:
    //    std::less<double> is not a strict weak ordering on NaN, and MSVC's
    //    lower_bound then treats NaN as EQUAL to whatever key it lands on, so
    //    the insert is silently dropped. For
    //    (0, 0, 0, -2.085220648210411, -5.801142781619129) upstream reports
    //    one root, -infinity of multiplicity 2, where the port reports three
    //    (-infinity and two NaNs). That is upstream undefined behaviour, not a
    //    port defect, so those inputs stay out of the generator.
    //
    // 2. Cancellation. The cubic and quartic closed forms build roots as
    //    differences of libm-derived terms, so a root much smaller than the
    //    largest one inherits the 1 ulp libm disagreement amplified by
    //    max|root| / |root|. Records whose smallest NONZERO root is more than
    //    1e6 times smaller than the largest are rejected, which bounds the
    //    relative disagreement by about 1e-15 * 1e6 = 1e-9, the tolerance the
    //    replay uses. Exactly zero roots come from the exact branches and are
    //    ignored by the test.
    bool MapComparable(std::map<double, int32_t> const& rmMap, bool requireSeparation)
    {
        double maxAbs = 0.0, minAbs = 0.0;
        bool haveMin = false;
        for (auto const& rm : rmMap)
        {
            if (!std::isfinite(rm.first)) { return false; }
            double a = std::fabs(rm.first);
            if (a > maxAbs) { maxAbs = a; }
            if (a != 0.0 && (!haveMin || a < minAbs)) { minAbs = a; haveMin = true; }
        }
        if (!requireSeparation || !haveMin || maxAbs == 0.0) { return true; }
        return minAbs * 1.0e6 >= maxAbs;
    }
}

namespace
{
    // Draw coefficients with a nonzero leading term until upstream's own
    // result passes MapComparable. The loop is capped and redraws every
    // coefficient; the fallback is the polynomial (x-1)(x-2)...(x-degree),
    // whose roots are exact and well separated.
    template <typename Solver>
    void DrawComparable(oracle::Ctx& io, int degree, bool requireSeparation,
        Solver solve, double* c)
    {
        std::vector<double> v;
        bool ok = false;
        std::map<double, int32_t> rmMap;
        for (int attempt = 0; attempt < 64 && !ok; ++attempt)
        {
            DrawCoefficients(io, degree, false, io.index() % 8, v);
            if (v[static_cast<size_t>(degree)] == 0.0)
            {
                v[static_cast<size_t>(degree)] = 1.0;
            }
            rmMap.clear();
            solve(v.data(), rmMap);
            ok = MapComparable(rmMap, requireSeparation);
        }
        if (!ok)
        {
            std::vector<double> roots;
            for (int i = 1; i <= degree; ++i) { roots.push_back(static_cast<double>(i)); }
            ExpandRoots(1.0, roots, v);
        }
        for (int i = 0; i <= degree; ++i) { c[i] = io.given(v[i]); }
    }
}

// The depressed quadratic uses std::sqrt only, which is correctly rounded in
// both runtimes, so this case is exact.
ORACLE_CASE("RootsPolynomial.solveQuadratic")
{
    double c[3]{};
    DrawComparable(io, 2, false, [](double const* v, std::map<double, int32_t>& m)
        { RootsPolynomial<double>::SolveQuadratic(v[0], v[1], v[2], m); }, c);
    std::map<double, int32_t> rmMap;
    RootsPolynomial<double>::SolveQuadratic(c[0], c[1], c[2], rmMap);
    EmitRootMap(io, rmMap);
}

ORACLE_CASE("RootsPolynomial.solveCubic")
{
    double c[4]{};
    DrawComparable(io, 3, true, [](double const* v, std::map<double, int32_t>& m)
        { RootsPolynomial<double>::SolveCubic(v[0], v[1], v[2], v[3], m); }, c);
    std::map<double, int32_t> rmMap;
    RootsPolynomial<double>::SolveCubic(c[0], c[1], c[2], c[3], rmMap);
    EmitRootMap(io, rmMap);
}

ORACLE_CASE("RootsPolynomial.solveQuartic")
{
    double c[5]{};
    DrawComparable(io, 4, true, [](double const* v, std::map<double, int32_t>& m)
        { RootsPolynomial<double>::SolveQuartic(v[0], v[1], v[2], v[3], v[4], m); }, c);
    std::map<double, int32_t> rmMap;
    RootsPolynomial<double>::SolveQuartic(c[0], c[1], c[2], c[3], c[4], rmMap);
    EmitRootMap(io, rmMap);
}

ORACLE_CASE("RootsPolynomial.getRootInfo")
{
    double c[5]{};
    DrawNonzeroLeading(io, 4, c);
    std::vector<int32_t> info;
    RootsPolynomial<double>::GetRootInfoQuadratic(c[0], c[1], c[2], info);
    io.outInt(info.size());
    for (int32_t m : info) { io.outInt(m); }
    RootsPolynomial<double>::GetRootInfoCubic(c[0], c[1], c[2], c[3], info);
    io.outInt(info.size());
    for (int32_t m : info) { io.outInt(m); }
    RootsPolynomial<double>::GetRootInfoQuartic(c[0], c[1], c[2], c[3], c[4], info);
    io.outInt(info.size());
    for (int32_t m : info) { io.outInt(m); }
}

ORACLE_CASE("RootsPolynomial.find")
{
    int degree = io.integer(0, 6);
    std::vector<double> v;
    DrawCoefficients(io, degree, false, io.index() % 8, v);
    std::vector<double> c(v.size());
    for (int i = 0; i <= degree; ++i) { c[i] = io.given(v[i]); }
    int maxIterations = io.integer(1, 128);
    std::vector<double> roots(static_cast<size_t>(degree) + 1, 0.0);
    int32_t numRoots = RootsPolynomial<double>::Find(degree, c.data(),
        static_cast<uint32_t>(maxIterations), roots.data());
    io.outInt(numRoots);
    for (int32_t i = 0; i < numRoots; ++i) { io.outReal(roots[i]); }
}

ORACLE_CASE("RootsPolynomial.findBounded")
{
    int degree = io.integer(1, 5);
    std::vector<double> v;
    DrawCoefficients(io, degree, false, io.index() % 8, v);
    std::vector<double> c(v.size());
    for (int i = 0; i <= degree; ++i) { c[i] = io.given(v[i]); }
    double tmin = io.lattice(-6, 1);
    double tmax = io.lattice(-1, 6);
    int maxIterations = io.integer(0, 128);
    double root = 0.0;
    bool found = RootsPolynomial<double>::Find(degree, c.data(), tmin, tmax,
        static_cast<uint32_t>(maxIterations), root);
    io.outBool(found);
    io.outReal(root);
}

// ---------------------------------------------------------------------------
// RootsGeneralPolynomial.h
// ---------------------------------------------------------------------------

ORACLE_CASE("RootsGeneralPolynomial.solve")
{
    // The leading coefficient is forced nonzero: with high-order zeros the
    // port and upstream deliberately differ (see the deviation case).
    int degree = io.integer(1, 5);
    std::vector<double> v;
    DrawCoefficients(io, degree, false, io.index() % 8, v);
    if (v[static_cast<size_t>(degree)] == 0.0) { v[static_cast<size_t>(degree)] = 1.0; }
    std::vector<double> p(v.size());
    for (int i = 0; i <= degree; ++i) { p[i] = io.given(v[i]); }
    bool useThreading = io.boolean();
    std::vector<double> roots;
    RootsGeneralPolynomial<double>::Solve(p, useThreading, roots);
    io.outInt(roots.size());
    for (double r : roots) { io.outReal(r); }
}

ORACLE_CASE("RootsGeneralPolynomial.solveRational")
{
    int degree = io.integer(1, 4);
    std::vector<double> v;
    DrawCoefficients(io, degree, false, io.index() % 8, v);
    if (v[static_cast<size_t>(degree)] == 0.0) { v[static_cast<size_t>(degree)] = 1.0; }
    std::vector<Rational> rP(v.size());
    for (int i = 0; i <= degree; ++i) { rP[i] = Rational(io.given(v[i])); }
    std::vector<Rational> rRoots;
    RootsGeneralPolynomial<double>::Solve(rP, false, rRoots);
    io.outInt(rRoots.size());
    for (auto const& r : rRoots) { io.outReal(static_cast<double>(r)); }
}

// Deviation: with high-order zero coefficients upstream allocates the
// rational coefficient vector at the UNTRIMMED size and assigns only the
// trimmed entries, so the solver runs on a zero-padded polynomial: the
// leading coefficient is zero, the Cauchy bound includes the monic 1 of the
// real leading term, and the derivative recursion descends through phantom
// degrees whose base case is bisected rather than solved exactly. The port
// allocates degree+1 entries. UPSTREAM-FINDINGS RootsGeneralPolynomial.h
// item 4, issue #340.
ORACLE_CASE("RootsGeneralPolynomial.deviation.zeroPadding")
{
    int degree = io.integer(2, 4);
    int numZeros = io.integer(1, 2);
    // The effect is a last-bit one on most draws, so the generator rejects
    // until the padded run really differs from the unpadded one, which is
    // exactly what the port computes. The loop is capped and redraws every
    // coefficient; the fallback is the last draw.
    std::vector<double> v, padded, rootsPadded, rootsPlain;
    for (int attempt = 0; attempt < 200; ++attempt)
    {
        DrawCoefficients(io, degree, false, io.index() % 8, v);
        if (v[static_cast<size_t>(degree)] == 0.0) { v[static_cast<size_t>(degree)] = 1.0; }
        padded.assign(static_cast<size_t>(degree) + 1 + numZeros, 0.0);
        for (int i = 0; i <= degree; ++i) { padded[i] = v[i]; }
        v.resize(static_cast<size_t>(degree) + 1);
        RootsGeneralPolynomial<double>::Solve(padded, false, rootsPadded);
        RootsGeneralPolynomial<double>::Solve(v, false, rootsPlain);
        if (rootsPadded != rootsPlain) { break; }
    }
    std::vector<double> p(static_cast<size_t>(degree) + 1 + numZeros, 0.0);
    for (int i = 0; i <= degree; ++i) { p[i] = io.given(v[i]); }
    std::vector<double> roots;
    RootsGeneralPolynomial<double>::Solve(p, false, roots);
    io.outInt(roots.size());
    for (double r : roots) { io.outReal(r); }
}
