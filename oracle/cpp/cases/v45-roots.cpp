// Group 45 (category roots-polynomials): RootsCubic.h, RootsQuartic.h.
//
// Both headers classify roots and multiplicities with exact rational
// arithmetic (BSRational<UIntegerAP32>) and compute root estimates either by
// bisection (useBisection = true, arithmetic only: +, *, /, comparisons and
// gte::FMA) or by a closed form (useBisection = false, which calls the C math
// library through std::pow, std::atan2, std::cos and std::sin on BSRational
// operands; those overloads convert to double, call libm and convert back).
//
//   * every "bisection" case is compared bit-for-bit (exact: true)
//   * every "closedForm" case is compared with a tolerance, because MSVC and
//     V8 disagree by an ulp on pow/atan2/cos/sin
//
// The port deliberately fixes three upstream defects recorded in
// docs/UPSTREAM-FINDINGS.md under RootsCubic.h/RootsQuartic.h (issue #340):
//   1. the cubic's invalid bisection bound max(1,|d0|,|d1|),
//   2. the quartic's incomplete "two complex-conjugate pairs" test,
//   3. the quartic's stale reads of a reused square-root array.
// The main cases reject the inputs on which upstream is defective with the
// probes below, which are verbatim copies of upstream's own control flow, and
// each defect gets its own 'deviation' case.
#define ORACLE_FAMILY "v45-roots"
#include "Oracle.h"

#include <Mathematics/RootsCubic.h>
#include <Mathematics/RootsQuartic.h>
#include <array>
#include <cmath>
#include <vector>

using namespace gte;

namespace
{
    using Rational = BSRational<UIntegerAP32>;

    // ---- generators -------------------------------------------------------
    //
    // Every mode fills 'c' with exactly degree+1 values (c[0] is the constant
    // term) using unrecorded raw draws; the caller records them with
    // io.given(), one statement per value.

    // Expand lead * prod_i (x - r[i]) into c[0..r.size()].
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

    // Multiply the monic quadratics (x^2 + p0*x + q0) and (x^2 + p1*x + q1)
    // by 'lead'. Integer p,q with p^2 - 4q < 0 give a complex-conjugate pair,
    // so this mode reaches the "no real roots" and "two real roots"
    // classifications with exactly zero rounding error.
    void ExpandQuadraticFactors(double lead, double p0, double q0, double p1, double q1,
        std::vector<double>& c)
    {
        c.assign(5, 0.0);
        c[0] = lead * (q0 * q1);
        c[1] = lead * (p0 * q1 + p1 * q0);
        c[2] = lead * (q0 + q1 + p0 * p1);
        c[3] = lead * (p0 + p1);
        c[4] = lead;
    }

    double NonzeroRawInteger(oracle::Ctx& io, int lo, int hi)
    {
        int v = 0;
        do { v = io.rawInteger(lo, hi); } while (v == 0);
        return static_cast<double>(v);
    }

    // The three entry points of every Roots* class. GENERAL records all
    // degree+1 coefficients, MONIC records the low 'degree' ones (the leading
    // coefficient is 1 by construction) and DEPRESSED records the low
    // degree-1 ones (the leading coefficient is 1 and the next one is 0).
    enum class Form { GENERAL, MONIC, DEPRESSED };

    // Integer roots that sum to zero, so that the expansion is a depressed
    // polynomial with integer coefficients.
    void DepressedRoots(std::vector<double>& roots)
    {
        double sum = 0.0;
        for (size_t i = 0; i + 1 < roots.size(); ++i) { sum += roots[i]; }
        roots[roots.size() - 1] = -sum;
    }

    // Draw degree+1 coefficients (index 0 is the constant term). 'mode'
    // selects the generator; every mode produces the same number of values.
    void DrawCoefficients(oracle::Ctx& io, int degree, Form form, int mode,
        std::vector<double>& c)
    {
        size_t const n = static_cast<size_t>(degree) + 1;
        bool const monic = (form != Form::GENERAL);
        bool const depressed = (form == Form::DEPRESSED);
        size_t const free = n - (monic ? 1 : 0) - (depressed ? 1 : 0);
        c.assign(n, 0.0);
        if (monic) { c[n - 1] = 1.0; }
        std::vector<double> roots;
        switch (mode)
        {
        case 0:
            // Small integer lattice: exact arithmetic, frequent degenerate
            // classifications (zero coefficients, zero discriminant).
            for (size_t i = 0; i < free; ++i)
            {
                c[i] = static_cast<double>(io.rawInteger(-5, 5));
            }
            break;
        case 1:
            // Uniform reals: the generic branch.
            for (size_t i = 0; i < free; ++i) { c[i] = io.raw(-10.0, 10.0); }
            break;
        case 2:
            // Distinct integer roots: all real, all multiplicity 1.
            for (int i = 0; i < degree; ++i)
            {
                roots.push_back(static_cast<double>(io.rawInteger(-6, 6)));
            }
            if (depressed) { DepressedRoots(roots); }
            ExpandRoots(monic ? 1.0 : NonzeroRawInteger(io, -4, 4), roots, c);
            break;
        case 3:
        {
            // Integer roots with a multiplicity pattern: double, triple,
            // quadruple roots and two double roots. The coefficients are
            // integers, so the rational discriminant is exactly zero and the
            // exact classifier branches are reached at exact equality.
            double a = static_cast<double>(io.rawInteger(-4, 4));
            double b = static_cast<double>(io.rawInteger(-4, 4));
            int pattern = io.rawInteger(0, degree == 3 ? 1 : 2);
            if (degree == 3)
            {
                // (x-a)^2 (x-b) and (x-a)^3; depressed: b = -2a and a = 0.
                if (pattern == 0) { roots = { a, a, depressed ? -2.0 * a : b }; }
                else { roots = { depressed ? 0.0 : a, depressed ? 0.0 : a,
                    depressed ? 0.0 : a }; }
            }
            else
            {
                // (x-a)^2 (x-b)^2, (x-a)^3 (x-b) and (x-a)^4; depressed:
                // b = -a, b = -3a and a = 0.
                if (pattern == 0) { roots = { a, a, depressed ? -a : b,
                    depressed ? -a : b }; }
                else if (pattern == 1) { roots = { a, a, a,
                    depressed ? -3.0 * a : b }; }
                else { roots = { depressed ? 0.0 : a, depressed ? 0.0 : a,
                    depressed ? 0.0 : a, depressed ? 0.0 : a }; }
            }
            ExpandRoots(monic ? 1.0 : NonzeroRawInteger(io, -4, 4), roots, c);
            break;
        }
        case 4:
        {
            // A common power-of-ten scale. For the general form, scaling every
            // coefficient leaves the roots unchanged but drives the rational
            // classifiers to extreme exponents; for the monic forms only the
            // low coefficients are scaled, which moves the roots instead.
            double s = std::pow(10.0, static_cast<double>(io.rawInteger(-100, 100)));
            for (size_t i = 0; i < free; ++i) { c[i] = io.raw(-10.0, 10.0) * s; }
            if (!monic) { c[n - 1] = io.raw(-10.0, 10.0) * s; }
            break;
        }
        case 5:
            // Per-coefficient scales: badly scaled polynomials.
            for (size_t i = 0; i < free; ++i)
            {
                c[i] = io.raw(-10.0, 10.0)
                    * std::pow(10.0, static_cast<double>(io.rawInteger(-30, 30)));
            }
            break;
        case 6:
            // Degree drop for the general form (the leading coefficient is
            // zero, sometimes two of them, which hands the polynomial to the
            // lower-degree solver); for the monic forms, zero interior
            // coefficients.
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
            // Zero-valued roots, including the signed zero -0.0 and a run of
            // zero low-order coefficients.
            for (size_t i = 0; i < free; ++i) { c[i] = io.raw(-10.0, 10.0); }
            if (!monic) { c[n - 1] = io.raw(-10.0, 10.0); }
            c[0] = (io.rawInteger(0, 1) == 0 ? 0.0 : -0.0);
            if (free >= 2 && io.rawInteger(0, 1) == 0) { c[1] = 0.0; }
            if (free >= 3 && io.rawInteger(0, 2) == 0) { c[2] = 0.0; }
            break;
        }
        if (monic) { c[n - 1] = 1.0; }
        if (depressed) { c[n - 2] = 0.0; }
    }

    // ---- probes -----------------------------------------------------------
    //
    // Verbatim copies of the upstream control flow that decides whether the
    // port's deliberate fixes change the answer. ORACLE.md allows a verbatim
    // copy of the upstream routine as a probe.

    // RootsCubic<T>::ComputeClassifiers, monic form.
    void CubicClassifiers(Rational const& rM0, Rational const& rM1, Rational const& rM2,
        Rational& rD0, Rational& rD1, Rational& rM2Div3)
    {
        rM2Div3 = Rational(1, 3) * rM2;
        rD0 = rM0 - rM2Div3 * (rM1 - Rational(2) * rM2Div3 * rM2Div3);
        rD1 = rM1 - rM2 * rM2Div3;
    }

    // RootsQuartic<T>::ComputeClassifiers, monic form.
    void QuarticClassifiers(Rational const& rM0, Rational const& rM1, Rational const& rM2,
        Rational const& rM3, Rational& rD0, Rational& rD1, Rational& rD2, Rational& rM3Div4)
    {
        rM3Div4 = Rational(1, 4) * rM3;
        Rational rM3Div4Sqr = rM3Div4 * rM3Div4;
        rD0 = rM0 - rM3Div4 * (rM1 - rM3Div4 * (rM2 - Rational(3) * rM3Div4Sqr));
        rD1 = rM1 - Rational(2) * rM3Div4 * (rM2 - Rational(4) * rM3Div4Sqr);
        rD2 = rM2 - Rational(6) * rM3Div4Sqr;
    }

    Rational CubicDelta(Rational const& rD0, Rational const& rD1)
    {
        return Rational(-27) * rD0 * rD0 + Rational(-4) * rD1 * rD1 * rD1;
    }

    Rational QuarticDelta(Rational const& rD0, Rational const& rD1, Rational const& rD2)
    {
        Rational rD0sqr = rD0 * rD0, rD1sqr = rD1 * rD1, rD2sqr = rD2 * rD2;
        return rD1sqr * (Rational(-27) * rD1sqr +
            Rational(4) * rD2 * (Rational(36) * rD0 - rD2sqr)) +
            Rational(16) * rD0 * (rD2sqr * (rD2sqr - Rational(8) * rD0) +
                Rational(16) * rD0sqr);
    }

    // True when upstream's bisection bound b = max(1,|d0|,|d1|) really does
    // bracket the single real root of x^3 + d1*x + d0, that is, when
    // PolynomialRootBisect's own sign preconditions hold. This is exactly the
    // condition under which the port (which falls back to Cauchy's bound
    // 1 + max(|d0|,|d1|) only when the upstream bound fails) is bit-identical
    // to upstream. UPSTREAM-FINDINGS RootsCubic.h item 1, issue #340.
    bool CubicBisectionSound(Rational const& rD0, Rational const& rD1)
    {
        if (rD0.GetSign() == 0 || rD1.GetSign() == 0) { return true; }
        if (CubicDelta(rD0, rD1).GetSign() >= 0) { return true; }
        double d0 = rD0;
        double d1 = rD1;
        double b = std::max(1.0, std::max(std::fabs(d0), std::fabs(d1)));
        auto F = [&d0, &d1](double x) { return gte::FMA(x, gte::FMA(x, x, d1), d0); };
        return F(-b) < 0.0 && F(b) > 0.0;
    }

    // ---- conditioning of the closed forms ---------------------------------
    //
    // ComputeDepressedRootsClosedForm builds its roots out of std::pow,
    // std::atan2, std::cos and std::sin results, each of which MSVC and V8 may
    // round differently in the last bit. A root that is the difference of two
    // nearly equal such terms amplifies that 1e-16 relative difference by the
    // cancellation factor max|term| / |root|: for x^3 + d1*x + d0 with roots
    // near +-9.4e16 and -11.3, the middle root of the closed form differs
    // between the two libms in its leading digit while every digit of the
    // outer roots agrees.
    //
    // The closed-form cases therefore accept a record only when that factor
    // stays below CANCELLATION_LIMIT, which bounds the relative disagreement
    // by about 1e-15 * 1e6 = 1e-9, the tolerance the replay uses. The
    // ill-conditioned regime is not dropped: it gets its own case that emits
    // the root count and the multiplicities, which exact rational arithmetic
    // decides and which are therefore comparable there too.
    double const CANCELLATION_LIMIT = 1.0e6;

    bool WellConditioned(std::vector<double> const& terms,
        PolynomialRoot<double> const* roots, size_t numRoots)
    {
        double scale = 0.0;
        for (double t : terms) { scale = std::max(scale, std::fabs(t)); }
        for (size_t i = 0; i < numRoots; ++i)
        {
            scale = std::max(scale, std::fabs(roots[i].x));
        }
        if (scale == 0.0) { return true; }
        for (size_t i = 0; i < numRoots; ++i)
        {
            if (!(std::fabs(roots[i].x) * CANCELLATION_LIMIT >= scale)) { return false; }
        }
        return true;
    }

    // The intermediate magnitudes of ComputeDepressedRootsClosedForm whose
    // cancellation the roots inherit. Exact branches (rational arithmetic or
    // std::sqrt, which is correctly rounded on both sides) contribute nothing.
    void CubicClosedFormTerms(Rational const& rD0, Rational const& rD1,
        std::vector<double>& terms)
    {
        int32_t signD0 = rD0.GetSign(), signD1 = rD1.GetSign();
        if (signD0 == 0 || signD1 == 0) { return; }
        Rational rDelta = CubicDelta(rD0, rD1);
        int32_t signDelta = rDelta.GetSign();
        if (signDelta == 0) { return; }
        Rational const r1Div3(1, 3);
        if (signDelta > 0)
        {
            // The roots are 2*t0, -t0-t1 and -t0+t1, so the roots themselves
            // already carry the term magnitudes (|t0| = |2*t0| / 2).
            return;
        }
        Rational rSqrtNegDeltaDiv27 = std::sqrt(-rDelta / Rational(27));
        Rational rD1Div3 = rD1 * r1Div3;
        Rational rCbrt = (signD0 < 0
            ? std::pow(Rational(1, 2) * (-rD0 + rSqrtNegDeltaDiv27), r1Div3)
            : -std::pow(Rational(1, 2) * (rD0 + rSqrtNegDeltaDiv27), r1Div3));
        if (rCbrt.GetSign() == 0) { return; }
        terms.push_back(static_cast<double>(rCbrt));
        terms.push_back(static_cast<double>(rD1Div3 / rCbrt));
    }

    // The value upstream reads out of rQRoots[1] after
    //   RootsQuadratic<Rational>::ComputeDepressedRoots(useBisection, -v, rQRoots)
    // when v is strictly positive (two real roots are written).
    Rational SqrtViaQuadratic(bool useBisection, Rational const& v)
    {
        std::array<PolynomialRoot<Rational>, 2> rQRoots{};
        RootsQuadratic<Rational>::ComputeDepressedRoots(useBisection, -v, rQRoots.data());
        return rQRoots[1].x;
    }

    // True when the port's RootsQuartic reproduces upstream bit for bit on
    // the depressed quartic d0 + d1*x + d2*x^2 + x^4. It replays upstream's
    // branch structure, including the single reused square-root array whose
    // index 1 upstream reads even when the call wrote nothing (item 3), and
    // rejects the inputs of items 1 to 3 of the #340 findings.
    bool QuarticDepressedSound(bool useBisection,
        Rational const& rD0, Rational const& rD1, Rational const& rD2)
    {
        int32_t signD0 = rD0.GetSign(), signD1 = rD1.GetSign();
        if (signD0 == 0)
        {
            if (signD1 == 0) { return true; }
            // The depressed cubic (rD1, rD2) supplies the other roots.
            return !useBisection || CubicBisectionSound(rD1, rD2);
        }

        if (signD1 == 0)
        {
            // SolveBiquadratic. The only unsound configuration is
            // rS + sqrt(rT) == 0: upstream then divides by zero and reads a
            // stale rQRoots[1].
            Rational rS = Rational(-0.5) * rD2;
            Rational rT = rS * rS - rD0;
            if (rT.GetSign() <= 0) { return true; }
            Rational rSPsqrtT = rS + SqrtViaQuadratic(useBisection, rT);
            return rSPsqrtT.GetSign() != 0;
        }

        Rational rDelta = QuarticDelta(rD0, rD1, rD2);
        int32_t signDelta = rDelta.GetSign();
        if (signDelta == 0) { return true; }
        if (signDelta > 0 && rD2.GetSign() > 0) { return true; }
        if (signDelta > 0)
        {
            // Item 2: upstream falls through and reports four roots for a
            // quartic with two complex-conjugate pairs.
            if ((Rational(4) * rD0 - rD2 * rD2).GetSign() >= 0) { return false; }
        }

        Rational rD1sqr = rD1 * rD1;
        Rational rM0 = Rational(0.125) * (Rational(4) * rD0 * rD2 - rD1sqr);
        Rational rM1 = -rD0;
        Rational rM2 = Rational(-0.5) * rD2;
        Rational rM2Div3 = Rational(1, 3) * rM2;
        Rational rC0 = rM0 - rM2Div3 * (rM1 - Rational(2) * rM2Div3 * rM2Div3);
        Rational rC1 = rM1 - rM2 * rM2Div3;
        if (useBisection && !CubicBisectionSound(rC0, rC1)) { return false; }

        std::array<PolynomialRoot<Rational>, 3> rCRoots{};
        RootsCubic<Rational>::ComputeDepressedRoots(useBisection, rC0, rC1, rCRoots.data());
        Rational rT = rCRoots[0].x - rM2Div3;

        // 'state' is upstream's rQRoots[1], which starts at zero and is
        // overwritten only by a call whose argument is strictly negative.
        Rational state(0);
        Rational rAlphaSqr = Rational(2) * rT - rD2;
        if (rAlphaSqr.GetSign() > 0) { state = SqrtViaQuadratic(useBisection, rAlphaSqr); }

        Rational rArg = rT * rT - rD0;
        if (rArg.GetSign() > 0) { state = SqrtViaQuadratic(useBisection, rArg); }
        else if (state.GetSign() != 0) { return false; }
        Rational rBeta = (signD1 > 0 ? Rational(1) : Rational(-1)) * state;
        if (rArg.GetSign() <= 0) { rBeta = Rational(0); }

        Rational rDiscr0 = rAlphaSqr - Rational(4) * (rT + rBeta);
        bool ok0 = (rDiscr0.GetSign() > 0 || state.GetSign() == 0);
        if (rDiscr0.GetSign() > 0) { state = SqrtViaQuadratic(useBisection, rDiscr0); }

        Rational rDiscr1 = rAlphaSqr - Rational(4) * (rT - rBeta);
        bool ok1 = (rDiscr1.GetSign() > 0 || state.GetSign() == 0);

        // For a negative discriminant only one of the two square roots is
        // used, so only that one has to agree.
        if (signDelta > 0) { return ok0 && ok1; }
        return signD1 > 0 ? ok1 : ok0;
    }

    // Soundness of RootsCubic<T>::Solve's dispatch. The quadratic and linear
    // sub-solvers are always bit-identical, so only the depressed cubic is
    // checked.
    bool CubicSolveSound(bool useBisection, double g0, double g1, double g2, double g3)
    {
        if (!useBisection) { return true; }
        if (g3 == 0.0 || g0 == 0.0) { return true; }
        Rational rD0{}, rD1{}, rM2Div3{};
        Rational rG3(g3);
        CubicClassifiers(Rational(g0) / rG3, Rational(g1) / rG3, Rational(g2) / rG3,
            rD0, rD1, rM2Div3);
        return CubicBisectionSound(rD0, rD1);
    }

    bool MonicCubicSolveSound(bool useBisection, double m0, double m1, double m2)
    {
        if (!useBisection) { return true; }
        if (m0 == 0.0) { return true; }
        Rational rD0{}, rD1{}, rM2Div3{};
        CubicClassifiers(Rational(m0), Rational(m1), Rational(m2), rD0, rD1, rM2Div3);
        return CubicBisectionSound(rD0, rD1);
    }

    // Closed-form comparability: the classification is exact on both sides,
    // but a root built by cancelling two libm-derived terms is not comparable
    // at 1e-9. These run upstream's own solver and measure the cancellation.
    // The depressed cubic's closed form calls std::pow if and only if d0 is
    // nonzero and the discriminant is nonzero; every other branch is exact
    // rational arithmetic or std::sqrt, which is correctly rounded on both
    // sides.
    bool CubicDepressedInexact(Rational const& rD0, Rational const& rD1)
    {
        return rD0.GetSign() != 0 && CubicDelta(rD0, rD1).GetSign() != 0;
    }

    bool DepressedCubicComparable(double d0, double d1)
    {
        Rational rD0(d0), rD1(d1);
        if (!CubicDepressedInexact(rD0, rD1)) { return true; }
        std::vector<double> terms;
        CubicClosedFormTerms(rD0, rD1, terms);
        std::array<PolynomialRoot<double>, 3> r{};
        size_t n = RootsCubic<double>::Solve(false, d0, d1, r.data());
        return WellConditioned(terms, r.data(), n);
    }

    bool MonicCubicComparable(double m0, double m1, double m2)
    {
        if (m0 == 0.0) { return true; }
        Rational rD0{}, rD1{}, rM2Div3{};
        CubicClassifiers(Rational(m0), Rational(m1), Rational(m2), rD0, rD1, rM2Div3);
        if (!CubicDepressedInexact(rD0, rD1)) { return true; }
        std::vector<double> terms{ static_cast<double>(rM2Div3) };
        CubicClosedFormTerms(rD0, rD1, terms);
        std::array<PolynomialRoot<double>, 3> r{};
        size_t n = RootsCubic<double>::Solve(false, m0, m1, m2, r.data());
        return WellConditioned(terms, r.data(), n);
    }

    bool CubicComparable(double g0, double g1, double g2, double g3)
    {
        if (g3 == 0.0 || g0 == 0.0) { return true; }
        Rational rD0{}, rD1{}, rM2Div3{}, rG3(g3);
        CubicClassifiers(Rational(g0) / rG3, Rational(g1) / rG3, Rational(g2) / rG3,
            rD0, rD1, rM2Div3);
        if (!CubicDepressedInexact(rD0, rD1)) { return true; }
        std::vector<double> terms{ static_cast<double>(rM2Div3) };
        CubicClosedFormTerms(rD0, rD1, terms);
        std::array<PolynomialRoot<double>, 3> r{};
        size_t n = RootsCubic<double>::Solve(false, g0, g1, g2, g3, r.data());
        return WellConditioned(terms, r.data(), n);
    }

    // The depressed quartic's closed form is inexact exactly when it solves a
    // resolvent (or sub-) cubic whose own closed form calls std::pow; every
    // other step is exact rational arithmetic or std::sqrt. 'terms' collects
    // the intermediate magnitudes whose cancellation the roots inherit.
    bool QuarticDepressedInexact(Rational const& rD0, Rational const& rD1,
        Rational const& rD2, std::vector<double>& terms)
    {
        int32_t signD0 = rD0.GetSign(), signD1 = rD1.GetSign();
        if (signD0 == 0)
        {
            if (signD1 == 0) { return false; }
            return CubicDepressedInexact(rD1, rD2)
                ? (CubicClosedFormTerms(rD1, rD2, terms), true) : false;
        }
        if (signD1 == 0) { return false; }
        Rational rDelta = QuarticDelta(rD0, rD1, rD2);
        int32_t signDelta = rDelta.GetSign();
        if (signDelta == 0) { return false; }
        if (signDelta > 0 && rD2.GetSign() > 0) { return false; }

        Rational rD1sqr = rD1 * rD1;
        Rational rM0 = Rational(0.125) * (Rational(4) * rD0 * rD2 - rD1sqr);
        Rational rM1 = -rD0;
        Rational rM2 = Rational(-0.5) * rD2;
        Rational rM2Div3 = Rational(1, 3) * rM2;
        Rational rC0 = rM0 - rM2Div3 * (rM1 - Rational(2) * rM2Div3 * rM2Div3);
        Rational rC1 = rM1 - rM2 * rM2Div3;
        if (!CubicDepressedInexact(rC0, rC1)) { return false; }
        CubicClosedFormTerms(rC0, rC1, terms);

        std::array<PolynomialRoot<Rational>, 3> rCRoots{};
        RootsCubic<Rational>::ComputeDepressedRoots(false, rC0, rC1, rCRoots.data());
        Rational rT = rCRoots[0].x - rM2Div3;
        terms.push_back(static_cast<double>(rT));
        terms.push_back(static_cast<double>(rM2Div3));
        Rational rAlphaSqr = Rational(2) * rT - rD2;
        if (rAlphaSqr.GetSign() > 0)
        {
            terms.push_back(static_cast<double>(SqrtViaQuadratic(false, rAlphaSqr)));
        }
        Rational rArg = rT * rT - rD0;
        if (rArg.GetSign() > 0)
        {
            terms.push_back(static_cast<double>(SqrtViaQuadratic(false, rArg)));
        }
        return true;
    }

    bool QuarticDepressedComparable(double d0, double d1, double d2)
    {
        std::vector<double> terms;
        if (!QuarticDepressedInexact(Rational(d0), Rational(d1), Rational(d2), terms))
        {
            return true;
        }
        std::array<PolynomialRoot<double>, 4> r{};
        size_t n = RootsQuartic<double>::Solve(false, d0, d1, d2, r.data());
        return WellConditioned(terms, r.data(), n);
    }

    bool MonicQuarticComparable(double m0, double m1, double m2, double m3)
    {
        if (m0 == 0.0) { return MonicCubicComparable(m1, m2, m3); }
        Rational rD0{}, rD1{}, rD2{}, rM3Div4{};
        QuarticClassifiers(Rational(m0), Rational(m1), Rational(m2), Rational(m3),
            rD0, rD1, rD2, rM3Div4);
        std::vector<double> terms{ static_cast<double>(rM3Div4) };
        if (!QuarticDepressedInexact(rD0, rD1, rD2, terms)) { return true; }
        std::array<PolynomialRoot<double>, 4> r{};
        size_t n = RootsQuartic<double>::Solve(false, m0, m1, m2, m3, r.data());
        return WellConditioned(terms, r.data(), n);
    }

    bool QuarticComparable(double g0, double g1, double g2, double g3, double g4)
    {
        if (g4 == 0.0) { return CubicComparable(g0, g1, g2, g3); }
        if (g0 == 0.0) { return CubicComparable(g1, g2, g3, g4); }
        Rational rD0{}, rD1{}, rD2{}, rM3Div4{}, rG4(g4);
        QuarticClassifiers(Rational(g0) / rG4, Rational(g1) / rG4, Rational(g2) / rG4,
            Rational(g3) / rG4, rD0, rD1, rD2, rM3Div4);
        std::vector<double> terms{ static_cast<double>(rM3Div4) };
        if (!QuarticDepressedInexact(rD0, rD1, rD2, terms)) { return true; }
        std::array<PolynomialRoot<double>, 4> r{};
        size_t n = RootsQuartic<double>::Solve(false, g0, g1, g2, g3, g4, r.data());
        return WellConditioned(terms, r.data(), n);
    }

    // The exact condition of UPSTREAM-FINDINGS RootsQuartic.h item 2: a
    // positive discriminant with P = 8*d2 >= 0 or D = 64*d0 - 16*d2^2 >= 0
    // means two complex-conjugate pairs, but upstream returns early only for
    // d2 > 0 and otherwise reports four spurious roots.
    bool QuarticItem2Defect(Rational const& rD0, Rational const& rD1, Rational const& rD2)
    {
        if (rD0.GetSign() == 0 || rD1.GetSign() == 0) { return false; }
        if (QuarticDelta(rD0, rD1, rD2).GetSign() <= 0) { return false; }
        if (rD2.GetSign() > 0) { return false; }
        return (Rational(4) * rD0 - rD2 * rD2).GetSign() >= 0;
    }

    bool QuarticSolveSound(bool useBisection, double g0, double g1, double g2, double g3,
        double g4)
    {
        if (g4 == 0.0) { return CubicSolveSound(useBisection, g0, g1, g2, g3); }
        if (g0 == 0.0)
        {
            if (g1 == 0.0) { return true; }
            return CubicSolveSound(useBisection, g1, g2, g3, g4);
        }
        Rational rD0{}, rD1{}, rD2{}, rM3Div4{};
        Rational rG4(g4);
        QuarticClassifiers(Rational(g0) / rG4, Rational(g1) / rG4, Rational(g2) / rG4,
            Rational(g3) / rG4, rD0, rD1, rD2, rM3Div4);
        return QuarticDepressedSound(useBisection, rD0, rD1, rD2);
    }

    bool MonicQuarticSolveSound(bool useBisection, double m0, double m1, double m2, double m3)
    {
        if (m0 == 0.0)
        {
            if (m1 == 0.0) { return true; }
            return MonicCubicSolveSound(useBisection, m1, m2, m3);
        }
        Rational rD0{}, rD1{}, rD2{}, rM3Div4{};
        QuarticClassifiers(Rational(m0), Rational(m1), Rational(m2), Rational(m3),
            rD0, rD1, rD2, rM3Div4);
        return QuarticDepressedSound(useBisection, rD0, rD1, rD2);
    }

    // ---- recording --------------------------------------------------------

    // Every rejection loop is capped and redraws all four (five) coefficients
    // on each attempt. Mode 3 builds the polynomial from repeated integer
    // roots, whose rational discriminant is exactly zero, so it is sound for
    // every solver; it is the fallback when the cap is reached.
    int NumRecorded(int degree, Form form)
    {
        return degree + 1 - (form == Form::GENERAL ? 0 : (form == Form::MONIC ? 1 : 2));
    }

    template <typename Predicate>
    void DrawAndRecord(oracle::Ctx& io, int degree, Form form, Predicate sound, double* c)
    {
        std::vector<double> v;
        int mode = io.index() % 8;
        bool ok = false;
        for (int attempt = 0; attempt < 64 && !ok; ++attempt)
        {
            DrawCoefficients(io, degree, form, mode, v);
            ok = sound(v);
        }
        if (!ok) { DrawCoefficients(io, degree, form, 3, v); }
        for (int i = 0; i < NumRecorded(degree, form); ++i) { c[i] = io.given(v[i]); }
    }

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

// ---------------------------------------------------------------------------
// RootsCubic.h
// ---------------------------------------------------------------------------

// The general cubic, bisection. Arithmetic only (BSRational exact arithmetic,
// gte::FMA, comparisons), so the comparison is bit-for-bit.
ORACLE_CASE("RootsCubic.solve.bisection")
{
    double c[4]{};
    DrawAndRecord(io, 3, Form::GENERAL, [](std::vector<double> const& v)
        { return CubicSolveSound(true, v[0], v[1], v[2], v[3]); }, c);
    std::array<PolynomialRoot<double>, 3> roots{};
    size_t numRoots = RootsCubic<double>::Solve(true, c[0], c[1], c[2], c[3], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

// The general cubic, closed form. std::pow, std::atan2, std::cos and std::sin
// on the depressed cubic make the root values libm-dependent; the root count
// and the multiplicities are decided by exact rational arithmetic.
ORACLE_CASE("RootsCubic.solve.closedForm")
{
    double c[4]{};
    DrawAndRecord(io, 3, Form::GENERAL, [](std::vector<double> const& v)
        { return CubicComparable(v[0], v[1], v[2], v[3]); }, c);
    std::array<PolynomialRoot<double>, 3> roots{};
    size_t numRoots = RootsCubic<double>::Solve(false, c[0], c[1], c[2], c[3], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsCubic.solveMonic.bisection")
{
    double c[3]{};
    DrawAndRecord(io, 3, Form::MONIC, [](std::vector<double> const& v)
        { return MonicCubicSolveSound(true, v[0], v[1], v[2]); }, c);
    std::array<PolynomialRoot<double>, 3> roots{};
    size_t numRoots = RootsCubic<double>::Solve(true, c[0], c[1], c[2], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsCubic.solveMonic.closedForm")
{
    double c[3]{};
    DrawAndRecord(io, 3, Form::MONIC, [](std::vector<double> const& v)
        { return MonicCubicComparable(v[0], v[1], v[2]); }, c);
    std::array<PolynomialRoot<double>, 3> roots{};
    size_t numRoots = RootsCubic<double>::Solve(false, c[0], c[1], c[2], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsCubic.solveDepressed.bisection")
{
    double c[2]{};
    DrawAndRecord(io, 3, Form::DEPRESSED, [](std::vector<double> const& v)
        { return CubicBisectionSound(Rational(v[0]), Rational(v[1])); }, c);
    std::array<PolynomialRoot<double>, 3> roots{};
    size_t numRoots = RootsCubic<double>::Solve(true, c[0], c[1], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsCubic.solveDepressed.closedForm")
{
    double c[2]{};
    DrawAndRecord(io, 3, Form::DEPRESSED, [](std::vector<double> const& v)
        { return DepressedCubicComparable(v[0], v[1]); }, c);
    std::array<PolynomialRoot<double>, 3> roots{};
    size_t numRoots = RootsCubic<double>::Solve(false, c[0], c[1], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

// The public ComputeDepressedRoots, entered with rational classifiers. The
// roots stay rational, so the only rounding is the final conversion for the
// golden file.
ORACLE_CASE("RootsCubic.computeDepressedRoots.bisection")
{
    double c[2]{};
    DrawAndRecord(io, 3, Form::DEPRESSED, [](std::vector<double> const& v)
        { return CubicBisectionSound(Rational(v[0]), Rational(v[1])); }, c);
    std::array<PolynomialRoot<Rational>, 3> roots{};
    size_t numRoots = RootsCubic<double>::ComputeDepressedRoots(true, Rational(c[0]),
        Rational(c[1]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsCubic.computeDepressedRoots.closedForm")
{
    double c[2]{};
    DrawAndRecord(io, 3, Form::DEPRESSED, [](std::vector<double> const& v)
        { return DepressedCubicComparable(v[0], v[1]); }, c);
    std::array<PolynomialRoot<Rational>, 3> roots{};
    size_t numRoots = RootsCubic<double>::ComputeDepressedRoots(false, Rational(c[0]),
        Rational(c[1]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

// The T = Rational instantiation: the inverse transform x - m2/3 is exact and
// the roots are converted to double only for the golden file.
ORACLE_CASE("RootsCubic.solveRational.bisection")
{
    double c[4]{};
    DrawAndRecord(io, 3, Form::GENERAL, [](std::vector<double> const& v)
        { return CubicSolveSound(true, v[0], v[1], v[2], v[3]); }, c);
    std::array<PolynomialRoot<Rational>, 3> roots{};
    size_t numRoots = RootsCubic<Rational>::Solve(true, Rational(c[0]), Rational(c[1]),
        Rational(c[2]), Rational(c[3]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsCubic.solveMonicRational.closedForm")
{
    double c[3]{};
    DrawAndRecord(io, 3, Form::MONIC, [](std::vector<double> const& v)
        { return MonicCubicComparable(v[0], v[1], v[2]); }, c);
    std::array<PolynomialRoot<Rational>, 3> roots{};
    size_t numRoots = RootsCubic<Rational>::Solve(false, Rational(c[0]), Rational(c[1]),
        Rational(c[2]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsCubic.solveDepressedRational.bisection")
{
    double c[2]{};
    DrawAndRecord(io, 3, Form::DEPRESSED, [](std::vector<double> const& v)
        { return CubicBisectionSound(Rational(v[0]), Rational(v[1])); }, c);
    std::array<PolynomialRoot<Rational>, 3> roots{};
    size_t numRoots = RootsCubic<Rational>::Solve(true, Rational(c[0]), Rational(c[1]),
        roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

// The ill-conditioned regime that the closed-form cases reject: badly scaled
// cubics whose middle root is the difference of two nearly equal libm-derived
// terms. Only the root count and the multiplicities are emitted; exact
// rational arithmetic decides both, so they are comparable here too.
ORACLE_CASE("RootsCubic.solveDepressed.closedForm.illConditioned")
{
    double c[2]{};
    DrawAndRecord(io, 3, Form::DEPRESSED, [](std::vector<double> const& v)
        { return !DepressedCubicComparable(v[0], v[1]); }, c);
    std::array<PolynomialRoot<double>, 3> roots{};
    size_t numRoots = RootsCubic<double>::Solve(false, c[0], c[1], roots.data());
    io.outInt(numRoots);
    for (size_t i = 0; i < numRoots; ++i) { io.outInt(roots[i].m); }
}

// Deviation: UPSTREAM-FINDINGS RootsCubic.h item 1, issue #340. The
// 'signDelta < 0' bisection branch bounds the single real root of
// x^3 + d1*x + d0 by max(1,|d0|,|d1|), which is not a root bound; the port
// uses Cauchy's bound 1 + max(|d0|,|d1|). |x|^3 <= |d1||x| + |d0| forces
// |x| <= sqrt(2b) for b = max(1,|d0|,|d1|), so the bound can only fail for
// b < 2 and the generator draws in [-2,2].
ORACLE_CASE("RootsCubic.solveDepressed.deviation.bisectionBound")
{
    double d0 = 0.0, d1 = 0.0;
    for (int attempt = 0; attempt < 400; ++attempt)
    {
        d0 = io.raw(-2.0, 2.0);
        d1 = io.raw(-2.0, 2.0);
        if (!CubicBisectionSound(Rational(d0), Rational(d1))) { break; }
    }
    io.given(d0);
    io.given(d1);
    std::array<PolynomialRoot<double>, 3> roots{};
    size_t numRoots = RootsCubic<double>::Solve(true, d0, d1, roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

// ---------------------------------------------------------------------------
// RootsQuartic.h
// ---------------------------------------------------------------------------

ORACLE_CASE("RootsQuartic.solve.bisection")
{
    double c[5]{};
    DrawAndRecord(io, 4, Form::GENERAL, [](std::vector<double> const& v)
        { return QuarticSolveSound(true, v[0], v[1], v[2], v[3], v[4]); }, c);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(true, c[0], c[1], c[2], c[3], c[4],
        roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solveMonic.bisection")
{
    double c[4]{};
    DrawAndRecord(io, 4, Form::MONIC, [](std::vector<double> const& v)
        { return MonicQuarticSolveSound(true, v[0], v[1], v[2], v[3]); }, c);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(true, c[0], c[1], c[2], c[3], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solveDepressed.bisection")
{
    double c[3]{};
    DrawAndRecord(io, 4, Form::DEPRESSED, [](std::vector<double> const& v)
        { return QuarticDepressedSound(true, Rational(v[0]), Rational(v[1]),
            Rational(v[2])); }, c);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(true, c[0], c[1], c[2], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.computeDepressedRoots.bisection")
{
    double c[3]{};
    DrawAndRecord(io, 4, Form::DEPRESSED, [](std::vector<double> const& v)
        { return QuarticDepressedSound(true, Rational(v[0]), Rational(v[1]),
            Rational(v[2])); }, c);
    std::array<PolynomialRoot<Rational>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::ComputeDepressedRoots(true, Rational(c[0]),
        Rational(c[1]), Rational(c[2]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solveRational.bisection")
{
    double c[5]{};
    DrawAndRecord(io, 4, Form::GENERAL, [](std::vector<double> const& v)
        { return QuarticSolveSound(true, v[0], v[1], v[2], v[3], v[4]); }, c);
    std::array<PolynomialRoot<Rational>, 4> roots{};
    size_t numRoots = RootsQuartic<Rational>::Solve(true, Rational(c[0]), Rational(c[1]),
        Rational(c[2]), Rational(c[3]), Rational(c[4]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solveDepressedRational.bisection")
{
    double c[3]{};
    DrawAndRecord(io, 4, Form::DEPRESSED, [](std::vector<double> const& v)
        { return QuarticDepressedSound(true, Rational(v[0]), Rational(v[1]),
            Rational(v[2])); }, c);
    std::array<PolynomialRoot<Rational>, 4> roots{};
    size_t numRoots = RootsQuartic<Rational>::Solve(true, Rational(c[0]), Rational(c[1]),
        Rational(c[2]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solveMonicRational.closedForm")
{
    double c[4]{};
    DrawAndRecord(io, 4, Form::MONIC, [](std::vector<double> const& v)
        { return MonicQuarticSolveSound(false, v[0], v[1], v[2], v[3])
            && MonicQuarticComparable(v[0], v[1], v[2], v[3]); }, c);
    std::array<PolynomialRoot<Rational>, 4> roots{};
    size_t numRoots = RootsQuartic<Rational>::Solve(false, Rational(c[0]), Rational(c[1]),
        Rational(c[2]), Rational(c[3]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solve.closedForm")
{
    double c[5]{};
    DrawAndRecord(io, 4, Form::GENERAL, [](std::vector<double> const& v)
        { return QuarticSolveSound(false, v[0], v[1], v[2], v[3], v[4])
            && QuarticComparable(v[0], v[1], v[2], v[3], v[4]); }, c);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(false, c[0], c[1], c[2], c[3], c[4],
        roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solveMonic.closedForm")
{
    double c[4]{};
    DrawAndRecord(io, 4, Form::MONIC, [](std::vector<double> const& v)
        { return MonicQuarticSolveSound(false, v[0], v[1], v[2], v[3])
            && MonicQuarticComparable(v[0], v[1], v[2], v[3]); }, c);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(false, c[0], c[1], c[2], c[3],
        roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solveDepressed.closedForm")
{
    double c[3]{};
    DrawAndRecord(io, 4, Form::DEPRESSED, [](std::vector<double> const& v)
        { return QuarticDepressedSound(false, Rational(v[0]), Rational(v[1]),
            Rational(v[2])) && QuarticDepressedComparable(v[0], v[1], v[2]); }, c);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(false, c[0], c[1], c[2], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.computeDepressedRoots.closedForm")
{
    double c[3]{};
    DrawAndRecord(io, 4, Form::DEPRESSED, [](std::vector<double> const& v)
        { return QuarticDepressedSound(false, Rational(v[0]), Rational(v[1]),
            Rational(v[2])) && QuarticDepressedComparable(v[0], v[1], v[2]); }, c);
    std::array<PolynomialRoot<Rational>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::ComputeDepressedRoots(false, Rational(c[0]),
        Rational(c[1]), Rational(c[2]), roots.data());
    EmitRationalRoots(io, roots.data(), numRoots);
}

// SolveBiquadratic: d1 is exactly zero, so d0 + d2*x^2 + x^4 factors through
// two nested square roots. The generator reaches all four classifications
// (four real roots, two real roots and a complex pair, two complex pairs, two
// real roots of multiplicity 2) by building (x^2 - a)(x^2 - b) from integers.
namespace
{
    void DrawBiquadratic(oracle::Ctx& io, bool useBisection, double* c)
    {
        double d0 = 0.0, d2 = 0.0;
        int mode = io.index() % 4;
        for (int attempt = 0; attempt < 64; ++attempt)
        {
            if (mode == 0)
            {
                double a = static_cast<double>(io.rawInteger(-4, 4));
                double b = static_cast<double>(io.rawInteger(-4, 4));
                d0 = a * b;
                d2 = -(a + b);
            }
            else if (mode == 1)
            {
                double a = static_cast<double>(io.rawInteger(-4, 4));
                d0 = a * a;
                d2 = -2.0 * a;
            }
            else if (mode == 2)
            {
                d0 = static_cast<double>(io.rawInteger(-5, 5));
                d2 = static_cast<double>(io.rawInteger(-5, 5));
            }
            else
            {
                d0 = io.raw(-10.0, 10.0);
                d2 = io.raw(-10.0, 10.0);
            }
            if (QuarticDepressedSound(useBisection, Rational(d0), Rational(0.0),
                Rational(d2)))
            {
                break;
            }
        }
        c[0] = io.given(d0);
        c[1] = io.given(0.0);
        c[2] = io.given(d2);
    }
}

ORACLE_CASE("RootsQuartic.solveBiquadratic.bisection")
{
    double c[3]{};
    DrawBiquadratic(io, true, c);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(true, c[0], c[1], c[2], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

ORACLE_CASE("RootsQuartic.solveBiquadratic.closedForm")
{
    double c[3]{};
    DrawBiquadratic(io, false, c);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(false, c[0], c[1], c[2], roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

// Deviation: UPSTREAM-FINDINGS RootsQuartic.h item 2, issue #340.
// (x^2 + p*x + q)(x^2 - p*x + r) with p^2 < 4q, p^2 < 4r and q + r <= p^2 has
// two complex-conjugate pairs, a positive discriminant and d2 = q + r - p^2
// <= 0, so upstream's 'd2 > 0' early-out misses it and it reports four roots
// with large residuals. The port returns none.
ORACLE_CASE("RootsQuartic.solveDepressed.deviation.complexPairs")
{
    double d0 = 0.0, d1 = 0.0, d2 = 0.0;
    for (int attempt = 0; attempt < 400; ++attempt)
    {
        double p = io.raw(1.0, 4.0);
        double pSqr = p * p;
        double q = io.raw(0.26 * pSqr, 0.5 * pSqr);
        double r = io.raw(0.26 * pSqr, pSqr - q);
        if (!(r > 0.25 * pSqr) || q == r) { continue; }
        d0 = q * r;
        d1 = p * (r - q);
        d2 = q + r - pSqr;
        if (QuarticItem2Defect(Rational(d0), Rational(d1), Rational(d2))) { break; }
    }
    io.given(d0);
    io.given(d1);
    io.given(d2);
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(true, d0, d1, d2, roots.data());
    EmitRoots(io, roots.data(), numRoots);
}

// Deviation: UPSTREAM-FINDINGS RootsQuartic.h item 3, issue #340. Upstream
// reads rQRoots[1] after every square-root extraction, but
// ComputeDepressedRoots writes index 1 only for a strictly negative argument.
// beta = sqrt(T^2 - d0) is mathematically zero exactly when d1 is zero, so a
// d1 of about 1e-200 makes T^2 - d0 round to a value whose sign is pure
// round-off: on the non-positive side upstream reuses alpha = sqrt(2T - d2)
// as beta, which is order 1, and the roots it reports are wrong. The port
// reads zero there.
ORACLE_CASE("RootsQuartic.solveDepressed.deviation.staleSqrt")
{
    double d0 = 0.0, d1 = 0.0, d2 = 0.0;
    for (int attempt = 0; attempt < 400; ++attempt)
    {
        d0 = io.raw(-4.0, 4.0);
        d2 = io.raw(-4.0, 4.0);
        d1 = (io.rawInteger(0, 1) == 0 ? 1.0 : -1.0) * io.raw(1.0e-200, 1.0e-199);
        Rational rD0(d0), rD1(d1), rD2(d2);
        if (QuarticItem2Defect(rD0, rD1, rD2)) { continue; }
        if (!QuarticDepressedSound(false, rD0, rD1, rD2)) { break; }
    }
    io.given(d0);
    io.given(d1);
    io.given(d2);
    // The closed form is used so that the invalid bisection bound of item 1
    // cannot be the cause of the disagreement.
    std::array<PolynomialRoot<double>, 4> roots{};
    size_t numRoots = RootsQuartic<double>::Solve(false, d0, d1, d2, roots.data());
    EmitRoots(io, roots.data(), numRoots);
}
