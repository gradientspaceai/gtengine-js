// Group 39 (numerical): Minimize1, MinimizeN, the four Ode solvers,
// SingularValueDecomposition, SymmetricEigensolver (NxN), RemezAlgorithm,
// CholeskyDecomposition and LDLTDecomposition.
//
// Everything in this family is + - * / sqrt except RemezAlgorithm's initial
// Chebyshev nodes, which call std::cos. Those angles form a fixed finite set
// (j * (GTE_C_HALF_PI / degree) for odd j < 2*degree), and the MSVC runtime
// and V8 were measured to agree bit for bit on every one of them for
// degree = 1..8 (the measurement is reported in oracle/reports/
// v39-numerical.md), so the Remez cases are compared exactly as well.
//
// IMPORTANT: this translation unit must NOT include Matrix2x2.h, Matrix3x3.h
// or Matrix4x4.h. OdeImplicitEuler::Update and BlockLDLTDecomposition::Factor
// call the unqualified Inverse(...) on a Matrix<N,N,double>; with those
// headers visible, ADL at the point of instantiation picks their closed-form
// overloads instead of Matrix.h's Gaussian elimination (see v34). The port
// implements the Gaussian-elimination form, and the OracleInverseIsGaussian
// case below verifies on every record that this TU really resolves to it.
#define ORACLE_FAMILY "v39-numerical"
#include "Oracle.h"

#include <Mathematics/CholeskyDecomposition.h>
#include <Mathematics/GaussianElimination.h>
#include <Mathematics/LDLTDecomposition.h>
#include <Mathematics/Minimize1.h>
#include <Mathematics/MinimizeN.h>
#include <Mathematics/OdeEuler.h>
#include <Mathematics/OdeImplicitEuler.h>
#include <Mathematics/OdeMidpoint.h>
#include <Mathematics/OdeRungeKutta4.h>
#include <Mathematics/RemezAlgorithm.h>
#include <Mathematics/SingularValueDecomposition.h>
#include <Mathematics/SymmetricEigensolver.h>

#include <cmath>
#include <cstdint>
#include <vector>

using namespace gte;

namespace
{
    // ------------------------------------------------------------ utilities

    // A power of two 2^k with |k| <= 60, exactly representable, used to drive
    // the generators far from unit scale without leaving the normal range.
    double DrawScale(oracle::Ctx& io, int32_t maxExponent)
    {
        int32_t k = io.integer(-maxExponent, maxExponent);
        return std::ldexp(1.0, k);
    }

    // ------------------------------------------- objective functions (1-D)
    //
    // Arithmetic only, so that the minimizer's control flow is comparable:
    // a libm call inside F would let a 1 ulp difference send the search to a
    // different local minimum, which no tolerance can repair.
    //
    // Every mode records exactly one integer and five doubles.
    struct Objective1
    {
        int32_t kind;
        double c[5];
    };

    Objective1 DrawObjective1(oracle::Ctx& io)
    {
        Objective1 f{};
        f.kind = io.integer(0, 3);
        switch (f.kind)
        {
        case 0:
            // Quartic with uniform coefficients, evaluated by Horner.
            for (int32_t i = 0; i < 5; ++i) { f.c[i] = io.real(-3.0, 3.0); }
            break;
        case 1:
            // Quartic on a small lattice: many exact ties and flat pieces.
            for (int32_t i = 0; i < 5; ++i) { f.c[i] = io.lattice(-3, 3); }
            break;
        case 2:
        {
            // Even quartic c0 + c2*t^2 + c4*t^4. On an interval symmetric
            // about 0 the two endpoints have exactly equal values, which is
            // the configuration of issue #298.
            double c0 = io.real(-2.0, 2.0);
            double c2 = io.real(-4.0, 4.0);
            double c4 = io.real(0.25, 3.0);
            f.c[0] = c0;
            f.c[1] = io.given(0.0);
            f.c[2] = c2;
            f.c[3] = io.given(0.0);
            f.c[4] = c4;
            break;
        }
        default:
        {
            // Rational: c0 + c1 / ((t - c2)^2 + c3), c3 > 0. Flat far from
            // c2 and sharply peaked near it.
            double c0 = io.real(-2.0, 2.0);
            double c1 = io.real(-6.0, 6.0);
            double c2 = io.real(-2.0, 2.0);
            double c3 = io.real(0.05, 2.0);
            f.c[0] = c0;
            f.c[1] = c1;
            f.c[2] = c2;
            f.c[3] = c3;
            f.c[4] = io.given(0.0);
            break;
        }
        }
        return f;
    }

    double Evaluate1(Objective1 const& f, double t)
    {
        if (f.kind == 3)
        {
            double d = t - f.c[2];
            return f.c[0] + f.c[1] / (d * d + f.c[3]);
        }
        double result = f.c[4];
        result = result * t + f.c[3];
        result = result * t + f.c[2];
        result = result * t + f.c[1];
        result = result * t + f.c[0];
        return result;
    }

    // ------------------------------------------------- Minimize1 (#298) probe
    //
    // A verbatim copy of Minimize1<double> that can run either upstream's
    // GetBracketedMinimum (mFixed == false) or the port's (mFixed == true).
    // The port differs in exactly two places:
    //   (a) when |denom| <= mEpsilon upstream returns, abandoning the search,
    //       where the port refines around the middle sample and continues;
    //   (b) upstream compares the parabola vertex with the middle sample
    //       exactly, where the port treats them as coincident when they are
    //       within 8 ulps of the bracket.
    // Running both and comparing (tMin, fMin) bit for bit is the exact
    // observable separator for issue #298: the main cases keep the records
    // where the two agree, the deviation case keeps the rest.
    class ProbeMinimize1
    {
    public:
        ProbeMinimize1(std::function<double(double)> const& F,
            int32_t maxSubdivisions, int32_t maxBisections, bool fixed,
            double epsilon = 1e-08, double tolerance = 1e-04)
            :
            mFunction(F),
            mMaxSubdivisions(maxSubdivisions),
            mMaxBisections(maxBisections),
            mTMin(0.0),
            mFMin(0.0),
            mEpsilon(std::max(epsilon, 0.0)),
            mTolerance(std::max(tolerance, 0.0)),
            mFixed(fixed)
        {
        }

        void GetMinimum(double t0, double t1, double tInitial, double& tMin, double& fMin)
        {
            mTMin = std::numeric_limits<double>::max();
            mFMin = std::numeric_limits<double>::max();

            double f0 = mFunction(t0);
            if (f0 < mFMin) { mTMin = t0; mFMin = f0; }

            double fInitial = mFunction(tInitial);
            if (fInitial < mFMin) { mTMin = tInitial; mFMin = fInitial; }

            double f1 = mFunction(t1);
            if (f1 < mFMin) { mTMin = t1; mFMin = f1; }

            if (((fInitial < f0) && (f1 >= fInitial)) ||
                ((f1 > fInitial) && (f0 >= fInitial)))
            {
                GetBracketedMinimum(t0, f0, tInitial, fInitial, t1, f1);
            }
            else
            {
                Subdivide(t0, f0, tInitial, fInitial, mMaxSubdivisions);
                Subdivide(tInitial, fInitial, t1, f1, mMaxSubdivisions);
            }

            tMin = mTMin;
            fMin = mFMin;
        }

    private:
        void Subdivide(double t0, double f0, double t1, double f1, int32_t subdivisionsRemaining)
        {
            if (subdivisionsRemaining-- == 0) { return; }

            double tm = 0.5 * (t0 + t1);
            double fm = mFunction(tm);
            if (fm < mFMin) { mTMin = tm; mFMin = fm; }

            if (((fm < f0) && (f1 >= fm)) || ((f1 > fm) && (f0 >= fm)))
            {
                GetBracketedMinimum(t0, f0, tm, fm, t1, f1);
            }
            else
            {
                Subdivide(t0, f0, tm, fm, subdivisionsRemaining);
                Subdivide(tm, fm, t1, f1, subdivisionsRemaining);
            }
        }

        // Upstream's "the vertex is at the middle sample" case, which the
        // port also uses for the degenerate parabola.
        void RefineAroundMidpoint(double& t0, double& f0, double& tm, double& fm,
            double& t1, double& f1)
        {
            double const half = 0.5;
            double tm0 = half * (t0 + tm);
            double fm0 = mFunction(tm0);
            double tm1 = half * (tm + t1);
            double fm1 = mFunction(tm1);

            if (fm0 < fm)
            {
                if (fm1 < fm)
                {
                    if (fm0 < fm1) { t1 = tm; f1 = fm; tm = tm0; fm = fm0; }
                    else { t0 = tm; f0 = fm; tm = tm1; fm = fm1; }
                }
                else { t1 = tm; f1 = fm; tm = tm0; fm = fm0; }
            }
            else if (fm0 > fm)
            {
                if (fm1 < fm) { t0 = tm; f0 = fm; tm = tm1; fm = fm1; }
                else { t0 = tm0; f0 = fm0; t1 = tm1; f1 = fm1; }
            }
            else
            {
                if (fm1 < fm) { t0 = tm; f0 = fm; tm = tm1; fm = fm1; }
                else { t0 = tm0; f0 = fm0; t1 = tm1; f1 = fm1; }
            }
        }

        void GetBracketedMinimum(double t0, double f0, double tm, double fm, double t1, double f1)
        {
            double const two = 2.0;
            double const half = 0.5;

            for (int32_t i = 0; i < mMaxBisections; ++i)
            {
                if (fm < mFMin) { mTMin = tm; mFMin = fm; }

                double dt10 = t1 - t0;
                double dtBound = two * mTolerance * std::fabs(tm) + mEpsilon;
                if (dt10 <= dtBound) { break; }

                double dt0m = t0 - tm;
                double dt1m = t1 - tm;
                double df0m = f0 - fm;
                double df1m = f1 - fm;
                double tmp0 = dt0m * df1m;
                double tmp1 = dt1m * df0m;
                double denom = tmp1 - tmp0;
                if (std::fabs(denom) <= mEpsilon)
                {
                    if (!mFixed) { return; }
                    RefineAroundMidpoint(t0, f0, tm, fm, t1, f1);
                    continue;
                }

                double tv = tm + half * (dt1m * tmp1 - dt0m * tmp0) / denom;
                tv = std::max(t0, std::min(tv, t1));
                double fv = mFunction(tv);
                if (fv < mFMin) { mTMin = tv; mFMin = fv; }

                double tvBound = mFixed
                    ? 8.0 * std::numeric_limits<double>::epsilon()
                        * std::max(std::fabs(t0), std::fabs(t1))
                    : 0.0;

                if (tv < tm - tvBound)
                {
                    if (fv < fm) { t1 = tm; f1 = fm; tm = tv; fm = fv; }
                    else { t0 = tv; f0 = fv; }
                }
                else if (tv > tm + tvBound)
                {
                    if (fv < fm) { t0 = tm; f0 = fm; tm = tv; fm = fv; }
                    else { t1 = tv; f1 = fv; }
                }
                else
                {
                    RefineAroundMidpoint(t0, f0, tm, fm, t1, f1);
                }
            }
        }

        std::function<double(double)> mFunction;
        int32_t mMaxSubdivisions;
        int32_t mMaxBisections;
        double mTMin, mFMin;
        double mEpsilon, mTolerance;
        bool mFixed;
    };

    // True when upstream's Minimize1 and the port's produce different
    // results on this configuration.
    bool Min1Differs(std::function<double(double)> const& F, double t0, double t1,
        double tInitial, int32_t maxSubdivisions, int32_t maxBisections,
        double epsilon, double tolerance)
    {
        ProbeMinimize1 up(F, maxSubdivisions, maxBisections, false, epsilon, tolerance);
        double tUp = 0.0, fUp = 0.0;
        up.GetMinimum(t0, t1, tInitial, tUp, fUp);
        ProbeMinimize1 port(F, maxSubdivisions, maxBisections, true, epsilon, tolerance);
        double tPort = 0.0, fPort = 0.0;
        port.GetMinimum(t0, t1, tInitial, tPort, fPort);
        return std::memcmp(&tUp, &tPort, sizeof(double)) != 0
            || std::memcmp(&fUp, &fPort, sizeof(double)) != 0;
    }

    // ------------------------------------------------------ Minimize1 draws

    struct Min1Draw
    {
        Objective1 f;
        double t0, t1, tInitial;
        int32_t maxSubdivisions, maxBisections;
        double epsilon, tolerance;
    };

    // Unrecorded draw. 'mode' cycles the coefficient family and the interval
    // construction; every mode produces the same 12 recorded values.
    Min1Draw RawMin1(oracle::Ctx& io, int32_t mode)
    {
        Min1Draw d{};
        d.f.kind = mode % 4;
        switch (d.f.kind)
        {
        case 0:
            for (int32_t i = 0; i < 5; ++i) { d.f.c[i] = io.raw(-3.0, 3.0); }
            break;
        case 1:
            for (int32_t i = 0; i < 5; ++i)
            {
                d.f.c[i] = static_cast<double>(io.rawInteger(-3, 3));
            }
            break;
        case 2:
            d.f.c[0] = io.raw(-2.0, 2.0);
            d.f.c[1] = 0.0;
            d.f.c[2] = io.raw(-4.0, 4.0);
            d.f.c[3] = 0.0;
            d.f.c[4] = io.raw(0.25, 3.0);
            break;
        default:
            d.f.c[0] = io.raw(-2.0, 2.0);
            d.f.c[1] = io.raw(-6.0, 6.0);
            d.f.c[2] = io.raw(-2.0, 2.0);
            d.f.c[3] = io.raw(0.05, 2.0);
            d.f.c[4] = 0.0;
            break;
        }

        double halfWidth = io.raw(0.5, 4.0);
        double center = io.raw(-2.0, 2.0);
        if (d.f.kind == 2)
        {
            // Symmetric about 0: F(t0) and F(t1) are exactly equal, which is
            // the bracket-collapse configuration of issue #298.
            d.t0 = -halfWidth;
            d.t1 = halfWidth;
        }
        else if (mode >= 4)
        {
            // A lattice interval, where the midpoints of the bisection are
            // exactly representable.
            d.t0 = static_cast<double>(io.rawInteger(-4, 0));
            d.t1 = static_cast<double>(io.rawInteger(1, 5));
        }
        else
        {
            d.t0 = center - halfWidth;
            d.t1 = center + halfWidth;
        }

        int32_t guess = io.rawInteger(0, 3);
        if (guess == 0) { d.tInitial = 0.5 * (d.t0 + d.t1); }
        else if (guess == 1) { d.tInitial = d.t0; }
        else if (guess == 2) { d.tInitial = d.t1; }
        else { d.tInitial = d.t0 + (d.t1 - d.t0) * io.raw(0.0, 1.0); }

        d.maxSubdivisions = io.rawInteger(1, 8);
        d.maxBisections = io.rawInteger(1, 48);
        d.epsilon = io.raw(1e-12, 1e-6);
        d.tolerance = io.raw(1e-8, 1e-3);
        return d;
    }

    void RecordMin1(oracle::Ctx& io, Min1Draw const& d)
    {
        io.given(static_cast<double>(d.f.kind));
        for (int32_t i = 0; i < 5; ++i) { io.given(d.f.c[i]); }
        io.given(d.t0);
        io.given(d.t1);
        io.given(d.tInitial);
        io.given(static_cast<double>(d.maxSubdivisions));
        io.given(static_cast<double>(d.maxBisections));
        io.given(d.epsilon);
        io.given(d.tolerance);
    }

    // True when upstream's Minimize1 and the port's report a different
    // minimum on this input; see ProbeMinimize1.
    bool Min1IsDefective(Min1Draw const& d)
    {
        Objective1 f = d.f;
        return Min1Differs([&f](double t) { return Evaluate1(f, t); },
            d.t0, d.t1, d.tInitial, d.maxSubdivisions, d.maxBisections,
            d.epsilon, d.tolerance);
    }

    // A configuration on which upstream is sound, used as the capped
    // rejection loop's fallback: F(t) = t^2 on [-1,2] with a well-resolved
    // parabola vertex.
    Min1Draw SoundMin1()
    {
        Min1Draw d{};
        d.f.kind = 0;
        d.f.c[0] = 0.0; d.f.c[1] = 0.0; d.f.c[2] = 1.0; d.f.c[3] = 0.0; d.f.c[4] = 0.0;
        d.t0 = -1.0;
        d.t1 = 2.0;
        d.tInitial = 0.25;
        d.maxSubdivisions = 8;
        d.maxBisections = 32;
        d.epsilon = 1e-8;
        d.tolerance = 1e-4;
        return d;
    }

    // ------------------------------------------ objective functions (N-D)
    //
    // Arithmetic only, for the same reason as Objective1.
    struct ObjectiveN
    {
        int32_t kind;
        int32_t dimensions;
        std::vector<double> w, a;
        double g;
    };

    double EvaluateN(ObjectiveN const& f, double const* x)
    {
        double result = 0.0;
        if (f.kind == 0)
        {
            // Anisotropic quadratic with a cross term: the cross term is
            // what Powell's conjugate directions are supposed to handle.
            for (int32_t i = 0; i < f.dimensions; ++i)
            {
                double t = x[i] - f.a[i];
                result += f.w[i] * t * t;
            }
            for (int32_t i = 0; i + 1 < f.dimensions; ++i)
            {
                result += f.g * x[i] * x[i + 1];
            }
        }
        else
        {
            // A Rosenbrock-style curved valley.
            for (int32_t i = 0; i + 1 < f.dimensions; ++i)
            {
                double u = x[i + 1] - x[i] * x[i];
                double v = f.a[i] - x[i];
                result += f.w[i] * u * u + v * v;
            }
        }
        return result;
    }

    // ------------------------------------------------- MinimizeN (#146) probe
    //
    // A verbatim copy of MinimizeN<double> with two switches: 'stray' selects
    // upstream's 'mDConjIndex = 0' (true) or the port's omission of it
    // (false), and the embedded ProbeMinimize1 reports whether the port's fix
    // of issue #298 would change any of the line searches.
    class ProbeMinimizeN
    {
    public:
        ProbeMinimizeN(int32_t dimensions, std::function<double(double const*)> const& F,
            int32_t maxLevel, int32_t maxBracket, int32_t maxIterations, bool stray,
            bool fixedMin1, double epsilon = 1e-06)
            :
            mDimensions(dimensions),
            mFunction(F),
            mMaxIterations(maxIterations),
            mEpsilon(0.0),
            mDirections(static_cast<size_t>(dimensions) + 1),
            mDConjIndex(dimensions),
            mDCurrIndex(0),
            mTCurr(dimensions),
            mTSave(dimensions),
            mFCurr(0.0),
            mStray(stray),
            mMinimizer([this](double t) { return mFunction(&(mTCurr + t * mDirections[mDCurrIndex])[0]); },
                maxLevel, maxBracket, fixedMin1)
        {
            mEpsilon = (epsilon > 0.0 ? epsilon : 0.0);
            for (auto& direction : mDirections) { direction.SetSize(dimensions); }
        }

        void GetMinimum(double const* t0, double const* t1, double const* tInitial,
            double* tMin, double& fMin)
        {
            size_t numBytes = mDimensions * sizeof(double);
            mFCurr = mFunction(tInitial);
            std::memcpy(&mTSave[0], tInitial, numBytes);
            std::memcpy(&mTCurr[0], tInitial, numBytes);

            for (int32_t i = 0; i < mDimensions; ++i) { mDirections[i].MakeUnit(i); }

            double ell0 = 0.0, ell1 = 0.0, ellMin = 0.0;
            for (int32_t iter = 0; iter < mMaxIterations; ++iter)
            {
                for (int32_t i = 0; i < mDimensions; ++i)
                {
                    mDCurrIndex = i;
                    ComputeDomain(t0, t1, ell0, ell1);
                    mMinimizer.GetMinimum(ell0, ell1, 0.0, ellMin, mFCurr);
                    mTCurr += ellMin * mDirections[i];
                }

                mDirections[mDConjIndex] = mTCurr - mTSave;
                double length = Length(mDirections[mDConjIndex]);
                if (length <= mEpsilon) { break; }

                mDirections[mDConjIndex] /= length;

                mDCurrIndex = mDConjIndex;
                ComputeDomain(t0, t1, ell0, ell1);
                mMinimizer.GetMinimum(ell0, ell1, 0.0, ellMin, mFCurr);
                mTCurr += ellMin * mDirections[mDCurrIndex];

                // UPSTREAM: mDConjIndex = 0. The port omits it (issue #146).
                if (mStray) { mDConjIndex = 0; }
                for (int32_t i = 0, ip1 = 1; i < mDimensions; ++i, ++ip1)
                {
                    mDirections[i] = mDirections[ip1];
                }

                mTSave = mTCurr;
            }

            std::memcpy(tMin, &mTCurr[0], numBytes);
            fMin = mFCurr;
        }

    private:
        void ComputeDomain(double const* t0, double const* t1, double& ell0, double& ell1)
        {
            ell0 = -std::numeric_limits<double>::max();
            ell1 = +std::numeric_limits<double>::max();

            for (int32_t i = 0; i < mDimensions; ++i)
            {
                double value = mDirections[mDCurrIndex][i];
                if (value != 0.0)
                {
                    double b0 = t0[i] - mTCurr[i];
                    double b1 = t1[i] - mTCurr[i];
                    double inv = 1.0 / value;
                    if (value > 0.0)
                    {
                        b0 *= inv;
                        if (b0 > ell0) { ell0 = b0; }
                        b1 *= inv;
                        if (b1 < ell1) { ell1 = b1; }
                    }
                    else
                    {
                        b0 *= inv;
                        if (b0 < ell1) { ell1 = b0; }
                        b1 *= inv;
                        if (b1 > ell0) { ell0 = b1; }
                    }
                }
            }

            if (ell0 > 0.0) { ell0 = 0.0; }
            if (ell1 < 0.0) { ell1 = 0.0; }
        }

        int32_t mDimensions;
        std::function<double(double const*)> mFunction;
        int32_t mMaxIterations;
        double mEpsilon;
        std::vector<GVector<double>> mDirections;
        int32_t mDConjIndex;
        int32_t mDCurrIndex;
        GVector<double> mTCurr;
        GVector<double> mTSave;
        double mFCurr;
        bool mStray;
        ProbeMinimize1 mMinimizer;
    };

    // ----------------------------------------------------- MinimizeN draws

    struct MinNDraw
    {
        ObjectiveN f;
        std::vector<double> t0, t1, tInitial;
        int32_t maxLevel, maxBracket, maxIterations;
        double epsilon;
    };

    MinNDraw RawMinN(oracle::Ctx& io, int32_t mode, int32_t minIterations, int32_t maxIterations)
    {
        MinNDraw d{};
        d.f.dimensions = io.rawInteger(2, 4);
        d.f.kind = mode % 2;
        int32_t n = d.f.dimensions;
        d.f.w.resize(n);
        d.f.a.resize(n);
        bool lattice = (mode >= 2);
        for (int32_t i = 0; i < n; ++i)
        {
            d.f.w[i] = lattice ? static_cast<double>(io.rawInteger(1, 4)) : io.raw(0.25, 4.0);
        }
        for (int32_t i = 0; i < n; ++i)
        {
            d.f.a[i] = lattice ? static_cast<double>(io.rawInteger(-2, 2)) : io.raw(-1.5, 1.5);
        }
        d.f.g = lattice ? static_cast<double>(io.rawInteger(-1, 1)) : io.raw(-1.0, 1.0);

        d.t0.resize(n);
        d.t1.resize(n);
        d.tInitial.resize(n);
        for (int32_t i = 0; i < n; ++i) { d.t0[i] = io.raw(-3.0, -1.0); }
        for (int32_t i = 0; i < n; ++i) { d.t1[i] = io.raw(1.0, 3.0); }
        for (int32_t i = 0; i < n; ++i)
        {
            d.tInitial[i] = d.t0[i] + (d.t1[i] - d.t0[i]) * io.raw(0.05, 0.95);
        }

        d.maxLevel = io.rawInteger(2, 6);
        d.maxBracket = io.rawInteger(8, 32);
        d.maxIterations = io.rawInteger(minIterations, maxIterations);
        d.epsilon = io.raw(1e-8, 1e-4);
        return d;
    }

    void RecordMinN(oracle::Ctx& io, MinNDraw const& d)
    {
        int32_t n = d.f.dimensions;
        io.given(static_cast<double>(n));
        io.given(static_cast<double>(d.f.kind));
        for (int32_t i = 0; i < n; ++i) { io.given(d.f.w[i]); }
        for (int32_t i = 0; i < n; ++i) { io.given(d.f.a[i]); }
        io.given(d.f.g);
        for (int32_t i = 0; i < n; ++i) { io.given(d.t0[i]); }
        for (int32_t i = 0; i < n; ++i) { io.given(d.t1[i]); }
        for (int32_t i = 0; i < n; ++i) { io.given(d.tInitial[i]); }
        io.given(static_cast<double>(d.maxLevel));
        io.given(static_cast<double>(d.maxBracket));
        io.given(static_cast<double>(d.maxIterations));
        io.given(d.epsilon);
    }

    // The capped rejection loops' fallback: a Rosenbrock valley in two
    // variables with a single iteration, on which both probes report that
    // upstream and the port agree. (An exactly quadratic objective is a poor
    // fallback: the parabola fit is exact, the next bracket is degenerate and
    // the issue #298 branch fires.)
    MinNDraw SoundMinN()
    {
        MinNDraw d{};
        d.f.kind = 1;
        d.f.dimensions = 2;
        d.f.w = { 3.0, 1.0 };
        d.f.a = { 0.5, -0.5 };
        d.f.g = 0.0;
        d.t0 = { -2.0, -2.0 };
        d.t1 = { 2.0, 2.0 };
        d.tInitial = { -1.25, 0.75 };
        d.maxLevel = 4;
        d.maxBracket = 16;
        d.maxIterations = 1;
        d.epsilon = 1e-6;
        return d;
    }

    // Run the probe as upstream (stray 'mDConjIndex = 0', unfixed Minimize1)
    // and as the port (neither), and compare the reported minima bit for
    // bit. That is the exact observable separator for the two deliberate
    // fixes that reach MinimizeN, issues #146 and #298.
    bool MinNDiffers(MinNDraw const& d)
    {
        ObjectiveN f = d.f;
        auto F = [&f](double const* x) { return EvaluateN(f, x); };
        int32_t n = d.f.dimensions;
        std::vector<double> tUp(n, 0.0), tPort(n, 0.0);
        double fUp = 0.0, fPort = 0.0;

        ProbeMinimizeN up(n, F, d.maxLevel, d.maxBracket, d.maxIterations,
            true, false, d.epsilon);
        up.GetMinimum(d.t0.data(), d.t1.data(), d.tInitial.data(), tUp.data(), fUp);
        ProbeMinimizeN port(n, F, d.maxLevel, d.maxBracket, d.maxIterations,
            false, true, d.epsilon);
        port.GetMinimum(d.t0.data(), d.t1.data(), d.tInitial.data(), tPort.data(), fPort);

        if (std::memcmp(&fUp, &fPort, sizeof(double)) != 0) { return true; }
        for (int32_t i = 0; i < n; ++i)
        {
            if (std::memcmp(&tUp[i], &tPort[i], sizeof(double)) != 0) { return true; }
        }
        return false;
    }

    // ----------------------------------------------------- Ode systems
    //
    // dx/dt = F(t,x) with F arithmetic only, written identically on both
    // sides. 'kind' 0 is linear, F[i] = b[i]*t + sum_j A(i,j)*x[j]; 'kind' 1
    // is quadratic, F[i] = b[i]*t + sum_j A(i,j)*x[j]*x[j]. The Jacobian
    // DF(r,c) is A(r,c) and 2*A(r,c)*x[c] respectively.
    template <int32_t N>
    struct OdeSystem
    {
        int32_t kind;
        Matrix<N, N, double> A;
        Vector<N, double> b;

        Vector<N, double> F(double t, Vector<N, double> const& x) const
        {
            Vector<N, double> f{};
            for (int32_t i = 0; i < N; ++i)
            {
                double s = b[i] * t;
                for (int32_t j = 0; j < N; ++j)
                {
                    s += (kind == 0 ? A(i, j) * x[j] : A(i, j) * x[j] * x[j]);
                }
                f[i] = s;
            }
            return f;
        }

        Matrix<N, N, double> DF(double t, Vector<N, double> const& x) const
        {
            (void)t;
            Matrix<N, N, double> df{};
            for (int32_t r = 0; r < N; ++r)
            {
                for (int32_t c = 0; c < N; ++c)
                {
                    df(r, c) = (kind == 0 ? A(r, c) : 2.0 * A(r, c) * x[c]);
                }
            }
            return df;
        }
    };

    // Draw and record a system. 'mode' selects the entry family:
    // 0 uniform, 1 small lattice (exact arithmetic), 2 diagonally dominant
    // and stiff (large negative diagonal), 3 nilpotent/singular (a zero row,
    // which makes I - tDelta*DF singular for the implicit solver when
    // combined with the right tDelta).
    template <int32_t N>
    OdeSystem<N> DrawOdeSystem(oracle::Ctx& io, int32_t mode, int32_t kind)
    {
        OdeSystem<N> sys{};
        sys.kind = kind;
        for (int32_t r = 0; r < N; ++r)
        {
            for (int32_t c = 0; c < N; ++c)
            {
                double value = 0.0;
                if (mode == 0) { value = io.real(-2.0, 2.0); }
                else if (mode == 1) { value = io.lattice(-3, 3); }
                else if (mode == 2)
                {
                    value = (r == c) ? io.real(-40.0, -10.0) : io.real(-0.5, 0.5);
                }
                else
                {
                    value = (r == 0) ? io.given(0.0) : io.lattice(-2, 2);
                }
                sys.A(r, c) = value;
            }
        }
        for (int32_t i = 0; i < N; ++i) { sys.b[i] = io.real(-1.0, 1.0); }
        return sys;
    }

    // Run 'numSteps' Update calls of one solver, changing tDelta between the
    // first and second step so that SetTDelta/GetTDelta are covered too.
    // solverId: 0 Euler, 1 Midpoint, 2 RungeKutta4, 3 ImplicitEuler.
    template <int32_t N>
    void RunOde(oracle::Ctx& io, int32_t solverId, int32_t mode, int32_t kind)
    {
        OdeSystem<N> sys = DrawOdeSystem<N>(io, mode, kind);
        double tDelta = (mode == 2 ? io.real(0.001, 0.05) : io.real(0.01, 0.5));
        double newTDelta = io.real(-0.25, 0.25);
        double t = io.real(-1.0, 1.0);
        Vector<N, double> x{};
        for (int32_t i = 0; i < N; ++i)
        {
            x[i] = (mode == 1 || mode == 3) ? io.lattice(-3, 3) : io.real(-2.0, 2.0);
        }
        int32_t numSteps = io.integer(1, 4);

        OdeSystem<N> const& s = sys;
        auto F = [&s](double tIn, Vector<N, double> const& xIn) { return s.F(tIn, xIn); };
        auto DF = [&s](double tIn, Vector<N, double> const& xIn) { return s.DF(tIn, xIn); };

        OdeEuler<double, Vector<N, double>> euler(tDelta, F);
        OdeMidpoint<double, Vector<N, double>> midpoint(tDelta, F);
        OdeRungeKutta4<double, Vector<N, double>> rk4(tDelta, F);
        OdeImplicitEuler<double, Vector<N, double>, Matrix<N, N, double>> implicit(tDelta, F, DF);
        OdeSolver<double, Vector<N, double>>* solver = nullptr;
        if (solverId == 0) { solver = &euler; }
        else if (solverId == 1) { solver = &midpoint; }
        else if (solverId == 2) { solver = &rk4; }
        else { solver = &implicit; }

        io.outReal(solver->GetTDelta());
        for (int32_t step = 0; step < numSteps; ++step)
        {
            double tOut = 0.0;
            Vector<N, double> xOut{};
            solver->Update(t, x, tOut, xOut);
            io.outReal(tOut);
            io.outVec(xOut);
            t = tOut;
            x = xOut;
            if (step == 0)
            {
                solver->SetTDelta(newTDelta);
                io.outReal(solver->GetTDelta());
            }
        }
    }

    // ---------------------------------------------- dense matrix generators

    // An MxN matrix drawn by 'mode' and recorded entry by entry (row major).
    // Modes: 0 uniform, 1 small lattice, 2 rank deficient (the last column is
    // an integer combination of the others), 3 diagonal with repeated
    // entries, 4 the zero matrix, 5 a lattice matrix scaled by 2^k.
    std::vector<double> DrawMatrix(oracle::Ctx& io, int32_t numRows, int32_t numCols,
        int32_t mode)
    {
        std::vector<double> A(static_cast<size_t>(numRows) * static_cast<size_t>(numCols), 0.0);
        auto at = [&A, numCols](int32_t r, int32_t c) -> double& {
            return A[static_cast<size_t>(c) + static_cast<size_t>(numCols) * static_cast<size_t>(r)];
        };

        if (mode == 0)
        {
            for (auto& value : A) { value = io.raw(-4.0, 4.0); }
        }
        else if (mode == 1)
        {
            for (auto& value : A) { value = static_cast<double>(io.rawInteger(-3, 3)); }
        }
        else if (mode == 2)
        {
            for (int32_t r = 0; r < numRows; ++r)
            {
                for (int32_t c = 0; c + 1 < numCols; ++c)
                {
                    at(r, c) = static_cast<double>(io.rawInteger(-3, 3));
                }
            }
            std::vector<double> w(static_cast<size_t>(numCols) - 1, 0.0);
            for (auto& value : w) { value = static_cast<double>(io.rawInteger(-2, 2)); }
            for (int32_t r = 0; r < numRows; ++r)
            {
                double s = 0.0;
                for (int32_t c = 0; c + 1 < numCols; ++c) { s += w[c] * at(r, c); }
                at(r, numCols - 1) = s;
            }
        }
        else if (mode == 3)
        {
            // Diagonal with a repeated value, which ties the singular values
            // and exercises the sort's treatment of equal keys.
            double repeated = static_cast<double>(io.rawInteger(1, 4));
            for (int32_t c = 0; c < numCols; ++c)
            {
                at(c, c) = (c < numCols / 2) ? repeated : repeated + 1.0;
            }
        }
        else if (mode == 4)
        {
            // Already zero.
            (void)io.rawInteger(0, 1);
        }
        else
        {
            double scale = std::ldexp(1.0, io.rawInteger(-200, 200));
            for (auto& value : A)
            {
                value = static_cast<double>(io.rawInteger(-3, 3)) * scale;
            }
        }

        for (auto& value : A) { io.given(value); }
        return A;
    }

    // An NxN symmetric matrix. Modes: 0 uniform, 1 lattice, 2 diagonal,
    // 3 tridiagonal, 4 block diagonal (two decoupled blocks), 5 rank one,
    // 6 zero, 7 diagonal with repeated entries, 8 a lattice matrix scaled by
    // a power of two. Unrecorded: the caller records the accepted draw.
    std::vector<double> DrawSymmetricRaw(oracle::Ctx& io, int32_t n, int32_t mode)
    {
        std::vector<double> A(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
        auto at = [&A, n](int32_t r, int32_t c) -> double& {
            return A[static_cast<size_t>(c) + static_cast<size_t>(n) * static_cast<size_t>(r)];
        };
        double scale = (mode == 8) ? std::ldexp(1.0, io.rawInteger(-200, 200)) : 1.0;

        if (mode == 5)
        {
            std::vector<double> v(static_cast<size_t>(n), 0.0);
            for (auto& value : v) { value = static_cast<double>(io.rawInteger(-3, 3)); }
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = 0; c < n; ++c) { at(r, c) = v[r] * v[c]; }
            }
        }
        else if (mode == 6)
        {
            (void)io.rawInteger(0, 1);
        }
        else if (mode == 7)
        {
            double repeated = static_cast<double>(io.rawInteger(1, 4));
            for (int32_t r = 0; r < n; ++r)
            {
                at(r, r) = (r < n / 2) ? repeated : repeated + 1.0;
            }
        }
        else
        {
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = r; c < n; ++c)
                {
                    bool keep = true;
                    if (mode == 2) { keep = (r == c); }
                    else if (mode == 3) { keep = (c - r <= 1); }
                    else if (mode == 4) { keep = ((r < n / 2) == (c < n / 2)); }
                    double value = 0.0;
                    if (keep)
                    {
                        value = (mode == 0) ? io.raw(-4.0, 4.0)
                            : static_cast<double>(io.rawInteger(-3, 3)) * scale;
                    }
                    at(r, c) = value;
                    at(c, r) = value;
                }
            }
        }

        return A;
    }

    // A verbatim copy of SymmetricEigensolver<double>::Tridiagonalize that
    // returns true as soon as a Householder step is degenerate
    // (length == 0). That is an exact characterization of the inputs on
    // which the port's fix of issue #80 changes the eigenvector matrix: the
    // fix only replaces the stored reflection parameter of such a step.
    // 'matrix' is a copy, since the routine overwrites it.
    bool EigenDecouples(std::vector<double> matrix, int32_t size)
    {
        std::vector<double> v(size, 0.0), p(size, 0.0), w(size, 0.0);
        int32_t r, c;
        for (int32_t i = 0, ip1 = 1; i < size - 2; ++i, ++ip1)
        {
            double length = 0.0;
            for (r = 0; r < ip1; ++r) { v[r] = 0.0; }
            for (r = ip1; r < size; ++r)
            {
                double vr = matrix[r + static_cast<size_t>(size) * i];
                v[r] = vr;
                length += vr * vr;
            }
            double vdv = 1.0;
            length = std::sqrt(length);
            if (length > 0.0)
            {
                double v1 = v[ip1];
                double sgn = (v1 >= 0.0 ? 1.0 : -1.0);
                double invDenom = 1.0 / (v1 + sgn * length);
                v[ip1] = 1.0;
                for (r = ip1 + 1; r < size; ++r)
                {
                    v[r] *= invDenom;
                    vdv += v[r] * v[r];
                }
            }
            else
            {
                return true;
            }

            double invvdv = 1.0 / vdv;
            double twoinvvdv = invvdv * 2.0;
            double pdvtvdv = 0.0;
            for (r = i; r < size; ++r)
            {
                p[r] = 0.0;
                for (c = i; c < r; ++c)
                {
                    p[r] += matrix[r + static_cast<size_t>(size) * c] * v[c];
                }
                for (/**/; c < size; ++c)
                {
                    p[r] += matrix[c + static_cast<size_t>(size) * r] * v[c];
                }
                p[r] *= twoinvvdv;
                pdvtvdv += p[r] * v[r];
            }

            pdvtvdv *= invvdv;
            for (r = i; r < size; ++r) { w[r] = p[r] - pdvtvdv * v[r]; }

            for (r = i; r < size; ++r)
            {
                double vr = v[r], wr = w[r];
                double offset = vr * wr * 2.0;
                matrix[r + static_cast<size_t>(size) * r] -= offset;
                for (c = r + 1; c < size; ++c)
                {
                    offset = vr * w[c] + wr * v[c];
                    matrix[c + static_cast<size_t>(size) * r] -= offset;
                }
            }

            matrix[i + static_cast<size_t>(size) * ip1] = twoinvvdv;
            for (r = ip1 + 1; r < size; ++r)
            {
                matrix[i + static_cast<size_t>(size) * r] = v[r];
            }
        }
        return false;
    }

    // ------------------------------------------------- SPD-ish generators
    //
    // Modes: 0 SPD from an integer lower-triangular L0 (A = L0*L0^T is exact
    // in integers), 1 SPD from a uniform L0, 2 positive semidefinite (one
    // diagonal entry of L0 is zero), 3 indefinite (a uniform symmetric
    // matrix), 4 the zero matrix, 5 an integer symmetric matrix whose (0,0)
    // entry is zero, which fails on the first pivot.
    std::vector<double> DrawSpd(oracle::Ctx& io, int32_t n, int32_t mode)
    {
        std::vector<double> A(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
        auto at = [&A, n](int32_t r, int32_t c) -> double& {
            return A[static_cast<size_t>(c) + static_cast<size_t>(n) * static_cast<size_t>(r)];
        };

        if (mode <= 2)
        {
            std::vector<double> L(static_cast<size_t>(n) * static_cast<size_t>(n), 0.0);
            int32_t zeroRow = (mode == 2) ? io.rawInteger(0, n - 1) : -1;
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = 0; c <= r; ++c)
                {
                    double value = 0.0;
                    if (c == r)
                    {
                        if (r == zeroRow) { value = 0.0; }
                        else if (mode == 1) { value = io.raw(0.5, 3.0); }
                        else { value = static_cast<double>(io.rawInteger(1, 3)); }
                    }
                    else
                    {
                        value = (mode == 1) ? io.raw(-2.0, 2.0)
                            : static_cast<double>(io.rawInteger(-2, 2));
                    }
                    L[static_cast<size_t>(c) + static_cast<size_t>(n) * r] = value;
                }
            }
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = 0; c < n; ++c)
                {
                    double s = 0.0;
                    for (int32_t k = 0; k < n; ++k)
                    {
                        s += L[static_cast<size_t>(k) + static_cast<size_t>(n) * r]
                            * L[static_cast<size_t>(k) + static_cast<size_t>(n) * c];
                    }
                    at(r, c) = s;
                }
            }
        }
        else if (mode == 3)
        {
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = r; c < n; ++c)
                {
                    double value = io.raw(-3.0, 3.0);
                    at(r, c) = value;
                    at(c, r) = value;
                }
            }
        }
        else if (mode == 4)
        {
            (void)io.rawInteger(0, 1);
        }
        else
        {
            for (int32_t r = 0; r < n; ++r)
            {
                for (int32_t c = r; c < n; ++c)
                {
                    double value = (r == 0 && c == 0) ? 0.0
                        : static_cast<double>(io.rawInteger(-3, 3));
                    at(r, c) = value;
                    at(c, r) = value;
                }
            }
        }

        for (double value : A) { io.given(value); }
        return A;
    }

    // Upstream's Solve returns std::numeric_limits<size_t>::max() (or
    // 0xFFFFFFFF) when the iteration does not converge, which exceeds 2^53.
    // Both sides emit -1 for it instead.
    void OutIterations(oracle::Ctx& io, size_t iterations, size_t invalid)
    {
        io.outInt(iterations == invalid ? -1 : static_cast<int64_t>(iterations));
    }

    void RunOdeDispatch(oracle::Ctx& io, int32_t solverId)
    {
        int32_t n = io.integer(2, 4);
        int32_t mode = io.integer(0, 3);
        int32_t kind = io.integer(0, 1);
        if (n == 2) { RunOde<2>(io, solverId, mode, kind); }
        else if (n == 3) { RunOde<3>(io, solverId, mode, kind); }
        else { RunOde<4>(io, solverId, mode, kind); }
    }

    // The same comparison for a second GetMinimum call on the same object,
    // where upstream's never-reset mDConjIndex matters (issue #146, part 2).
    bool MinNReuseDiffers(MinNDraw const& d, std::vector<double> const& second)
    {
        ObjectiveN f = d.f;
        auto F = [&f](double const* x) { return EvaluateN(f, x); };
        int32_t n = d.f.dimensions;
        std::vector<double> tUp(n, 0.0), tPort(n, 0.0), scratch(n, 0.0);
        double fUp = 0.0, fPort = 0.0, fScratch = 0.0;

        ProbeMinimizeN up(n, F, d.maxLevel, d.maxBracket, d.maxIterations,
            true, false, d.epsilon);
        up.GetMinimum(d.t0.data(), d.t1.data(), d.tInitial.data(), scratch.data(), fScratch);
        up.GetMinimum(d.t0.data(), d.t1.data(), second.data(), tUp.data(), fUp);
        ProbeMinimizeN port(n, F, d.maxLevel, d.maxBracket, d.maxIterations,
            false, true, d.epsilon);
        port.GetMinimum(d.t0.data(), d.t1.data(), d.tInitial.data(), scratch.data(), fScratch);
        port.GetMinimum(d.t0.data(), d.t1.data(), second.data(), tPort.data(), fPort);

        if (std::memcmp(&fUp, &fPort, sizeof(double)) != 0) { return true; }
        for (int32_t i = 0; i < n; ++i)
        {
            if (std::memcmp(&tUp[i], &tPort[i], sizeof(double)) != 0) { return true; }
        }
        return false;
    }
}

// =========================================================== Minimize1
//
// F is arithmetic only (Horner on a quartic, or one division for the
// rational family), so the whole search - the V-shape tests, the bracket
// updates and the parabola vertex - is IEEE-exact on both sides.

ORACLE_CASE("Minimize1.getMinimum")
{
    // Rejection sampling: keep only inputs on which the port's fix of
    // issue #298 does not change the iteration (see ProbeMinimize1).
    Min1Draw d{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        d = RawMin1(io, (io.index() + attempt) % 8);
        accepted = !Min1IsDefective(d);
    }
    if (!accepted) { d = SoundMin1(); }
    RecordMin1(io, d);

    Objective1 f = d.f;
    Minimize1<double> minimizer([&f](double t) { return Evaluate1(f, t); },
        d.maxSubdivisions, d.maxBisections, d.epsilon, d.tolerance);
    double tMin = 0.0, fMin = 0.0;
    minimizer.GetMinimum(d.t0, d.t1, d.tInitial, tMin, fMin);
    io.outReal(tMin);
    io.outReal(fMin);
    io.outReal(minimizer.GetEpsilon());
    io.outReal(minimizer.GetTolerance());
}

// The two-argument overload, whose initial guess is (t0+t1)/2.
ORACLE_CASE("Minimize1.getMinimum.defaultGuess")
{
    Min1Draw d{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        d = RawMin1(io, (io.index() + attempt) % 8);
        d.tInitial = 0.5 * (d.t0 + d.t1);
        accepted = !Min1IsDefective(d);
    }
    if (!accepted) { d = SoundMin1(); d.tInitial = 0.5 * (d.t0 + d.t1); }
    RecordMin1(io, d);

    Objective1 f = d.f;
    Minimize1<double> minimizer([&f](double t) { return Evaluate1(f, t); },
        d.maxSubdivisions, d.maxBisections, d.epsilon, d.tolerance);
    double tMin = 0.0, fMin = 0.0;
    minimizer.GetMinimum(d.t0, d.t1, tMin, fMin);
    io.outReal(tMin);
    io.outReal(fMin);
}

// The default epsilon and tolerance of the constructor, plus SetEpsilon and
// SetTolerance with possibly negative arguments (std::max(x, 0) clamps them
// to zero).
ORACLE_CASE("Minimize1.setEpsilonTolerance")
{
    Min1Draw d{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        d = RawMin1(io, (io.index() + attempt) % 8);
        d.epsilon = 1e-08;
        d.tolerance = 1e-04;
        accepted = !Min1IsDefective(d);
    }
    if (!accepted) { d = SoundMin1(); }
    RecordMin1(io, d);
    double newEpsilon = io.real(-1.0, 1.0);
    double newTolerance = io.real(-1.0, 1.0);

    Objective1 f = d.f;
    Minimize1<double> minimizer([&f](double t) { return Evaluate1(f, t); },
        d.maxSubdivisions, d.maxBisections);
    io.outReal(minimizer.GetEpsilon());
    io.outReal(minimizer.GetTolerance());
    double tMin = 0.0, fMin = 0.0;
    minimizer.GetMinimum(d.t0, d.t1, d.tInitial, tMin, fMin);
    io.outReal(tMin);
    io.outReal(fMin);
    minimizer.SetEpsilon(newEpsilon);
    minimizer.SetTolerance(newTolerance);
    io.outReal(minimizer.GetEpsilon());
    io.outReal(minimizer.GetTolerance());
}

// Issue #298: at equal-value bracket endpoints the vertex of the
// interpolating parabola lands within round-off of the middle sample, where
// upstream's exact 'tv == tm' test fails and the bracket collapses, or the
// parabola goes degenerate and upstream abandons the search. The port
// resolves the vertex to within 8 ulps of the bracket and continues.
ORACLE_CASE("Minimize1.getMinimum.deviation")
{
    Min1Draw d{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        // Mode 2 is the even quartic on an interval symmetric about 0.
        d = RawMin1(io, 2);
        d.maxBisections = io.rawInteger(8, 48);
        accepted = Min1IsDefective(d);
    }
    RecordMin1(io, d);

    Objective1 f = d.f;
    Minimize1<double> minimizer([&f](double t) { return Evaluate1(f, t); },
        d.maxSubdivisions, d.maxBisections, d.epsilon, d.tolerance);
    double tMin = 0.0, fMin = 0.0;
    minimizer.GetMinimum(d.t0, d.t1, d.tInitial, tMin, fMin);
    io.outReal(tMin);
    io.outReal(fMin);
}

// Issue #298, second part: with the initial guess at an endpoint the first
// bracket is {(t0,f0),(t1,f1),(t1,f1)}, the interpolating parabola is exactly
// degenerate (denom == 0), and upstream returns after three function
// evaluations, reporting the better endpoint. The port refines around the
// middle sample instead and finds the interior minimum. F(t) = c*t^2 here,
// whose minimum is 0 at t = 0, an interior point of every drawn interval, so
// each deviating record is a demonstrably wrong upstream answer.
ORACLE_CASE("Minimize1.getMinimum.endpointGuessDeviation")
{
    double t0 = 0.0, t1 = 0.0, tInitial = 0.0, c2 = 0.0;
    int32_t maxSubdivisions = 0, maxBisections = 0;
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        t0 = static_cast<double>(io.rawInteger(-4, -1));
        t1 = static_cast<double>(io.rawInteger(1, 4));
        tInitial = (io.rawInteger(0, 1) == 0) ? t0 : t1;
        c2 = io.raw(0.25, 2.0);
        maxSubdivisions = io.rawInteger(1, 6);
        maxBisections = io.rawInteger(4, 40);
        accepted = Min1Differs([c2](double t) { return c2 * t * t; },
            t0, t1, tInitial, maxSubdivisions, maxBisections, 1e-08, 1e-04);
    }
    io.given(t0);
    io.given(t1);
    io.given(tInitial);
    io.given(c2);
    io.given(static_cast<double>(maxSubdivisions));
    io.given(static_cast<double>(maxBisections));

    Minimize1<double> minimizer([c2](double t) { return c2 * t * t; },
        maxSubdivisions, maxBisections);
    double tMin = 0.0, fMin = 0.0;
    minimizer.GetMinimum(t0, t1, tInitial, tMin, fMin);
    io.outReal(tMin);
    io.outReal(fMin);
}

// Throw parity: the constructor asserts maxSubdivisions > 0 and
// maxBisections > 0, and GetMinimum asserts t0 <= tInitial <= t1.
ORACLE_CASE("Minimize1.invalidInput")
{
    // The point of the case is throw parity, so the configurations that do
    // NOT throw are restricted to those on which upstream's Minimize1 and
    // the port's agree. Without that restriction the generator lands on
    // issue #298 often: with tInitial at an endpoint the initial bracket is
    // degenerate, |denom| is exactly zero, upstream abandons the search and
    // returns the endpoint while the port finds the interior minimum.
    int32_t maxSubdivisions = 0, maxBisections = 0;
    double t0 = 0.0, t1 = 0.0, tInitial = 0.0, c2 = 0.0;
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        maxSubdivisions = io.rawInteger((io.index() + attempt) % 2 == 0 ? 1 : -1, 3);
        maxBisections = io.rawInteger((io.index() + attempt) % 3 == 0 ? 1 : -1, 3);
        t0 = static_cast<double>(io.rawInteger(-3, 0));
        t1 = static_cast<double>(io.rawInteger(1, 4));
        tInitial = ((io.index() + attempt) % 5 == 0) ? io.raw(-5.0, 6.0)
            : static_cast<double>(io.rawInteger(-3, 4));
        c2 = io.raw(0.25, 2.0);
        bool willThrow = (maxSubdivisions <= 0 || maxBisections <= 0
            || !(t0 <= tInitial && tInitial <= t1));
        accepted = willThrow || !Min1Differs([c2](double t) { return c2 * t * t; },
            t0, t1, tInitial, maxSubdivisions, maxBisections, 1e-08, 1e-04);
    }
    io.given(static_cast<double>(maxSubdivisions));
    io.given(static_cast<double>(maxBisections));
    io.given(t0);
    io.given(t1);
    io.given(tInitial);
    io.given(c2);

    Minimize1<double> minimizer([c2](double t) { return c2 * t * t; },
        maxSubdivisions, maxBisections);
    double tMin = 0.0, fMin = 0.0;
    minimizer.GetMinimum(t0, t1, tInitial, tMin, fMin);
    io.outReal(tMin);
    io.outReal(fMin);
}

// ================================================================ Ode
//
// Four explicit Update steps of small arithmetic-only systems. Every solver
// gets the same generator so that a formula difference shows in one of them
// and not the others.

ORACLE_CASE("OdeEuler.update")
{
    RunOdeDispatch(io, 0);
}

ORACLE_CASE("OdeMidpoint.update")
{
    RunOdeDispatch(io, 1);
}

ORACLE_CASE("OdeRungeKutta4.update")
{
    RunOdeDispatch(io, 2);
}

ORACLE_CASE("OdeImplicitEuler.update")
{
    RunOdeDispatch(io, 3);
}

// I - tDelta*DF exactly singular: Inverse returns the zero matrix with no
// warning, so the step is xOut = xIn and tOut = tIn + tDelta. The system has
// a zero first row, DF is then singular for kind 0, and tDelta is chosen so
// that I - tDelta*DF is singular as well: with a11 the (1,1) entry of DF and
// the rest of row/column 0 zero, tDelta = 1/a11 makes row 1 dependent.
ORACLE_CASE("OdeImplicitEuler.update.singularJacobian")
{
    // A 2x2 system with F[0] = b0*t and F[1] = b1*t + a*x[1], so that
    // DF = [[0,0],[0,a]] and I - tDelta*DF is singular exactly at
    // tDelta = 1/a.
    double a = io.lattice(1, 4);
    double b0 = io.real(-1.0, 1.0);
    double b1 = io.real(-1.0, 1.0);
    double t = io.real(-1.0, 1.0);
    double x0 = io.lattice(-3, 3);
    double x1 = io.lattice(-3, 3);
    double tDelta = io.given(1.0 / a);
    int32_t numSteps = io.integer(1, 3);

    auto F = [a, b0, b1](double tIn, Vector<2, double> const& xIn)
    {
        Vector<2, double> f{};
        f[0] = b0 * tIn;
        f[1] = b1 * tIn + a * xIn[1];
        return f;
    };
    auto DF = [a](double tIn, Vector<2, double> const& xIn)
    {
        (void)tIn; (void)xIn;
        Matrix<2, 2, double> df{};
        df(0, 0) = 0.0; df(0, 1) = 0.0;
        df(1, 0) = 0.0; df(1, 1) = a;
        return df;
    };

    OdeImplicitEuler<double, Vector<2, double>, Matrix<2, 2, double>> solver(tDelta, F, DF);
    Vector<2, double> x{};
    x[0] = x0;
    x[1] = x1;
    for (int32_t step = 0; step < numSteps; ++step)
    {
        double tOut = 0.0;
        Vector<2, double> xOut{};
        solver.Update(t, x, tOut, xOut);
        io.outReal(tOut);
        io.outVec(xOut);
        t = tOut;
        x = xOut;
    }
}

// This translation unit must resolve the unqualified Inverse(Matrix<N,N,T>)
// of OdeImplicitEuler and BlockLDLTDecomposition to Matrix.h's
// Gaussian-elimination template, not to a closed-form overload of
// Matrix2x2.h / Matrix3x3.h / Matrix4x4.h (see v34). The case compares the
// two on every record; the port implements the Gaussian-elimination form, so
// the booleans must all be true and the matrices must agree bit for bit.
ORACLE_CASE("Matrix.inverse.isGaussianElimination")
{
    Matrix<3, 3, double> M{};
    for (int32_t r = 0; r < 3; ++r)
    {
        for (int32_t c = 0; c < 3; ++c) { M(r, c) = io.real(-4.0, 4.0); }
    }
    Matrix<2, 2, double> M2{};
    for (int32_t r = 0; r < 2; ++r)
    {
        for (int32_t c = 0; c < 2; ++c) { M2(r, c) = io.real(-4.0, 4.0); }
    }

    bool invertible3 = false, invertible2 = false;
    Matrix<3, 3, double> inv3 = Inverse(M, &invertible3);
    Matrix<2, 2, double> inv2 = Inverse(M2, &invertible2);

    Matrix<3, 3, double> ge3{};
    Matrix<2, 2, double> ge2{};
    double det = 0.0;
    bool geInvertible3 = GaussianElimination<double>()(3, &M[0], &ge3[0], det,
        nullptr, nullptr, nullptr, 0, nullptr);
    bool geInvertible2 = GaussianElimination<double>()(2, &M2[0], &ge2[0], det,
        nullptr, nullptr, nullptr, 0, nullptr);

    bool same = (invertible3 == geInvertible3) && (invertible2 == geInvertible2);
    for (int32_t i = 0; i < 9 && same; ++i)
    {
        same = (std::memcmp(&inv3[i], &ge3[i], sizeof(double)) == 0);
    }
    for (int32_t i = 0; i < 4 && same; ++i)
    {
        same = (std::memcmp(&inv2[i], &ge2[i], sizeof(double)) == 0);
    }
    io.outBool(same);
    io.outBool(invertible3);
    io.outMat(inv3);
    io.outBool(invertible2);
    io.outMat(inv2);
}

// ========================================= SingularValueDecomposition
//
// Bidiagonalization plus Golub-Kahan steps: + - * / sqrt only. The sort of
// the singular values uses std::sort with a comparison on the value alone,
// so tied values leave the U/V column permutation unspecified; MSVC's
// std::sort is an insertion sort (hence stable) for at most 32 elements and
// numCols <= 4 here, and the port's Array.prototype.sort is stable, so the
// two agree. The .tiedSingularValues case is the one that ties them.

namespace
{
    void RunSvd(oracle::Ctx& io, int32_t mode)
    {
        int32_t numCols = io.integer(2, 4);
        int32_t numRows = io.integer(numCols, 5);
        int32_t maxIterations = io.integer(1, 48);
        double multiplier = io.real(0.5, 24.0);
        std::vector<double> A = DrawMatrix(io, numRows, numCols, mode);
        int32_t uIndex = io.integer(0, numRows - 1);
        int32_t vIndex = io.integer(0, numCols - 1);

        SingularValueDecomposition<double> svd(
            static_cast<size_t>(numRows), static_cast<size_t>(numCols),
            static_cast<size_t>(maxIterations));
        size_t iterations = svd.Solve(A.data(), multiplier);
        OutIterations(io, iterations, std::numeric_limits<size_t>::max());

        std::vector<double> singular(static_cast<size_t>(numCols), 0.0);
        svd.GetSingularValues(singular.data());
        for (double value : singular) { io.outReal(value); }

        std::vector<double> U(static_cast<size_t>(numRows) * numRows, 0.0);
        svd.GetU(U.data());
        for (double value : U) { io.outReal(value); }

        std::vector<double> V(static_cast<size_t>(numCols) * numCols, 0.0);
        svd.GetV(V.data());
        for (double value : V) { io.outReal(value); }

        std::vector<double> S(static_cast<size_t>(numRows) * numCols, 0.0);
        svd.GetS(S.data());
        for (double value : S) { io.outReal(value); }

        std::vector<double> uColumn(static_cast<size_t>(numRows), 0.0);
        svd.GetUColumn(static_cast<size_t>(uIndex), uColumn.data());
        for (double value : uColumn) { io.outReal(value); }

        std::vector<double> vColumn(static_cast<size_t>(numCols), 0.0);
        svd.GetVColumn(static_cast<size_t>(vIndex), vColumn.data());
        for (double value : vColumn) { io.outReal(value); }

        io.outReal(svd.GetSingularValue(static_cast<size_t>(vIndex)));
    }
}

ORACLE_CASE("SingularValueDecomposition.solve")
{
    RunSvd(io, io.integer(0, 1));
}

ORACLE_CASE("SingularValueDecomposition.solve.rankDeficient")
{
    RunSvd(io, 2);
}

ORACLE_CASE("SingularValueDecomposition.solve.tiedSingularValues")
{
    RunSvd(io, 3);
}

ORACLE_CASE("SingularValueDecomposition.solve.zeroMatrix")
{
    RunSvd(io, 4);
}

ORACLE_CASE("SingularValueDecomposition.solve.extremeScale")
{
    RunSvd(io, 5);
}

// Throw parity: the constructor asserts numCols >= 2, numRows >= numCols and
// maxIterations > 0, and Solve asserts multiplier > 0. Nothing is emitted
// before the calls that can throw.
ORACLE_CASE("SingularValueDecomposition.invalidInput")
{
    int32_t numRows = io.integer(1, 5);
    int32_t numCols = io.integer(1, 5);
    int32_t maxIterations = io.integer(0, 8);
    double multiplier = io.real(-2.0, 8.0);
    std::vector<double> A(static_cast<size_t>(numRows) * static_cast<size_t>(numCols), 0.0);
    for (auto& value : A) { value = io.real(-4.0, 4.0); }

    SingularValueDecomposition<double> svd(
        static_cast<size_t>(numRows), static_cast<size_t>(numCols),
        static_cast<size_t>(maxIterations));
    size_t iterations = svd.Solve(A.data(), multiplier);
    OutIterations(io, iterations, std::numeric_limits<size_t>::max());
}

// Throw parity for the accessor range checks (issue #478). Exactly one
// accessor is called per record, so a record either throws or produces its
// whole output.
ORACLE_CASE("SingularValueDecomposition.invalidIndex")
{
    int32_t numCols = io.integer(2, 4);
    int32_t numRows = io.integer(numCols, 5);
    int32_t maxIterations = io.integer(8, 32);
    int32_t accessor = io.integer(0, 2);
    int32_t index = io.integer(-1, 6);
    std::vector<double> A(static_cast<size_t>(numRows) * static_cast<size_t>(numCols), 0.0);
    for (auto& value : A) { value = io.real(-4.0, 4.0); }

    SingularValueDecomposition<double> svd(
        static_cast<size_t>(numRows), static_cast<size_t>(numCols),
        static_cast<size_t>(maxIterations));
    svd.Solve(A.data(), 8.0);
    if (accessor == 0)
    {
        std::vector<double> uColumn(static_cast<size_t>(numRows), 0.0);
        svd.GetUColumn(static_cast<size_t>(index), uColumn.data());
        for (double value : uColumn) { io.outReal(value); }
    }
    else if (accessor == 1)
    {
        std::vector<double> vColumn(static_cast<size_t>(numCols), 0.0);
        svd.GetVColumn(static_cast<size_t>(index), vColumn.data());
        for (double value : vColumn) { io.outReal(value); }
    }
    else
    {
        io.outReal(svd.GetSingularValue(static_cast<size_t>(index)));
    }
}

// ============================================== SymmetricEigensolver (NxN)

namespace
{
    // UPSTREAM HANG (see oracle/reports/v39-numerical.md). When Solve does
    // not converge it never calls ComputePermutation, so mPermutation stays
    // all zeros and GetEigenvectors' cycle walk
    //   while ((next = mPermutation[current]) != start)
    // never leaves index 0. Both the MSVC build and the port loop forever.
    // Every case that calls GetEigenvectors therefore requires a converged
    // Solve; SymmetricEigensolver.solve.nonConvergence covers the other
    // branch without touching GetEigenvectors.
    bool EigenConverges(std::vector<double> const& A, int32_t size, int32_t maxIterations,
        int32_t sortType)
    {
        SymmetricEigensolver<double> probe(size, static_cast<uint32_t>(maxIterations));
        return probe.Solve(A.data(), sortType) != 0xFFFFFFFFu;
    }

    void RunSymmetricEigensolver(oracle::Ctx& io, bool wantDecoupled,
        int32_t modeLo, int32_t modeHi)
    {
        int32_t size = 0, mode = 0, maxIterations = 0, sortType = 0;
        std::vector<double> A;
        bool accepted = false;
        for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
        {
            size = io.rawInteger(wantDecoupled ? 3 : 2, 6);
            mode = io.rawInteger(modeLo, modeHi);
            A = DrawSymmetricRaw(io, size, mode);
            maxIterations = io.rawInteger(8, 64);
            sortType = io.rawInteger(-1, 1);
            accepted = (EigenDecouples(A, size) == wantDecoupled)
                && EigenConverges(A, size, maxIterations, sortType);
        }
        if (!accepted)
        {
            // Fallbacks: a 3x3 lattice matrix whose first subcolumn is
            // nonzero (no degenerate step) or a diagonal one (every step
            // degenerate). Both converge.
            size = 3;
            maxIterations = 32;
            sortType = 1;
            A.assign(9, 0.0);
            if (wantDecoupled)
            {
                A[0] = 1.0; A[4] = 2.0; A[8] = 3.0;
            }
            else
            {
                A[0] = 1.0; A[1] = 2.0; A[2] = 3.0;
                A[3] = 2.0; A[4] = 4.0; A[5] = 1.0;
                A[6] = 3.0; A[7] = 1.0; A[8] = 5.0;
            }
        }

        io.given(static_cast<double>(size));
        for (double value : A) { io.given(value); }
        io.given(static_cast<double>(maxIterations));
        io.given(static_cast<double>(sortType));
        int32_t index = io.integer(0, size - 1);

        SymmetricEigensolver<double> solver(size, static_cast<uint32_t>(maxIterations));
        uint32_t iterations = solver.Solve(A.data(), sortType);
        OutIterations(io, iterations, 0xFFFFFFFFu);
        io.outInt(solver.GetEigenvectorMatrixType());

        std::vector<double> eigenvalues(static_cast<size_t>(size), 0.0);
        solver.GetEigenvalues(eigenvalues.data());
        for (double value : eigenvalues) { io.outReal(value); }

        std::vector<double> eigenvectors(static_cast<size_t>(size) * size, 0.0);
        solver.GetEigenvectors(eigenvectors.data());
        for (double value : eigenvectors) { io.outReal(value); }
        io.outInt(solver.GetEigenvectorMatrixType());

        std::vector<double> eigenvector(static_cast<size_t>(size), 0.0);
        solver.GetEigenvector(index, eigenvector.data());
        for (double value : eigenvector) { io.outReal(value); }
        io.outReal(solver.GetEigenvalue(index));
    }
}

ORACLE_CASE("SymmetricEigensolver.solve")
{
    RunSymmetricEigensolver(io, false, 0, 8);
}

// Issue #80: when the subcolumn below the subdiagonal is already zero the
// Householder reflection actually applied is the identity, but upstream
// stores 2/Dot(v,v) = 2 as its parameter, and GetEigenvectors /
// GetEigenvector rebuild I - 2*e*e^T from it. The port stores 0 there, so
// the rebuilt reflection is the identity. The eigenvalues are unaffected;
// the eigenvector matrix is not.
ORACLE_CASE("SymmetricEigensolver.solve.decoupledDeviation")
{
    RunSymmetricEigensolver(io, true, 2, 8);
}

// The non-convergence branch: Solve exhausts its budget and returns
// 0xFFFFFFFF. GetEigenvectors is deliberately NOT called here - after a
// failed Solve it loops forever on both sides; see the case-file comment on
// EigenConverges and the group report.
ORACLE_CASE("SymmetricEigensolver.solve.nonConvergence")
{
    int32_t size = 0, mode = 0, maxIterations = 0, sortType = 0;
    std::vector<double> A;
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 64 && !accepted; ++attempt)
    {
        size = io.rawInteger(3, 6);
        mode = io.rawInteger(0, 8);
        A = DrawSymmetricRaw(io, size, mode);
        maxIterations = io.rawInteger(1, 3);
        sortType = io.rawInteger(-1, 1);
        accepted = !EigenConverges(A, size, maxIterations, sortType);
    }
    if (!accepted)
    {
        // A 3x3 unreduced tridiagonal matrix with a single iteration.
        size = 3;
        maxIterations = 1;
        sortType = 1;
        A = { 1.0, 2.0, 0.0, 2.0, 4.0, 3.0, 0.0, 3.0, 5.0 };
    }

    io.given(static_cast<double>(size));
    for (double value : A) { io.given(value); }
    io.given(static_cast<double>(maxIterations));
    io.given(static_cast<double>(sortType));
    int32_t index = io.integer(0, size - 1);

    SymmetricEigensolver<double> solver(size, static_cast<uint32_t>(maxIterations));
    uint32_t iterations = solver.Solve(A.data(), sortType);
    OutIterations(io, iterations, 0xFFFFFFFFu);
    io.outInt(solver.GetEigenvectorMatrixType());
    std::vector<double> eigenvalues(static_cast<size_t>(size), 0.0);
    solver.GetEigenvalues(eigenvalues.data());
    for (double value : eigenvalues) { io.outReal(value); }
    std::vector<double> eigenvector(static_cast<size_t>(size), 0.0);
    solver.GetEigenvector(index, eigenvector.data());
    for (double value : eigenvector) { io.outReal(value); }
    io.outReal(solver.GetEigenvalue(index));
}

// A size of 1 or 0, or a zero iteration budget, leaves the solver inactive:
// Solve returns 0, the accessors do nothing and GetEigenvalue returns
// DBL_MAX.
ORACLE_CASE("SymmetricEigensolver.solve.invalidSize")
{
    int32_t size = io.integer(0, 2);
    // A valid size gets a zero budget, so the solver is inactive on every
    // record; a valid size with a small budget can fail to converge, which
    // is the deviation case above (issue #517), not this one.
    int32_t maxIterations = (size >= 2 ? io.integer(0, 0) : io.integer(0, 4));
    int32_t index = io.integer(0, 1);
    std::vector<double> A(static_cast<size_t>(size) * static_cast<size_t>(size), 0.0);
    for (auto& value : A) { value = io.lattice(-3, 3); }

    SymmetricEigensolver<double> solver(size, static_cast<uint32_t>(maxIterations));
    uint32_t iterations = solver.Solve(A.data(), 0);
    OutIterations(io, iterations, 0xFFFFFFFFu);
    io.outInt(solver.GetEigenvectorMatrixType());
    io.outReal(solver.GetEigenvalue(index));
}

// ======================================================= RemezAlgorithm
//
// F and F' are arithmetic only. The only libm call on upstream's path is the
// std::cos of ComputeInitialXNodes, whose arguments form a fixed finite set
// (j * (GTE_C_HALF_PI / degree) for odd j < 2*degree); MSVC and V8 were
// measured to agree bit for bit on every one of them for degree = 1..8, so
// these cases are compared exactly.
//
// Upstream's ComputePCoefficients reads poly[i] past the coefficient array
// when the Newton product loses its leading coefficient to exact
// cancellation (issue #147, fixed in the port by returning 0 there). That
// needs u[degree] - E*v[degree] to be exactly zero, which random draws never
// produce; upstream's out-of-range read stays out of these generators.

namespace
{
    struct RemezFunction
    {
        int32_t kind;
        double c[6];

        double F(double x) const
        {
            if (kind == 0)
            {
                double n = c[0] + c[1] * x;
                double d = c[2] + x * x;
                return n / d;
            }
            double result = c[5];
            result = result * x + c[4];
            result = result * x + c[3];
            result = result * x + c[2];
            result = result * x + c[1];
            result = result * x + c[0];
            return result;
        }

        double FDer(double x) const
        {
            if (kind == 0)
            {
                double n = c[0] + c[1] * x;
                double d = c[2] + x * x;
                return (c[1] * d - n * (2.0 * x)) / (d * d);
            }
            double result = 5.0 * c[5];
            result = result * x + 4.0 * c[4];
            result = result * x + 3.0 * c[3];
            result = result * x + 2.0 * c[2];
            result = result * x + c[1];
            return result;
        }
    };

    // 'fixedKind' is 0 or 1, or -1 to draw the kind. Either way exactly one
    // integer is recorded for it.
    void RunRemez(oracle::Ctx& io, int32_t fixedKind, int32_t degreeLo, int32_t degreeHi,
        int32_t maxRemezLo, int32_t maxRemezHi)
    {
        RemezFunction f{};
        int32_t kind = (fixedKind < 0) ? io.integer(0, 1) : io.integer(fixedKind, fixedKind);
        f.kind = kind;
        if (kind == 0)
        {
            f.c[0] = io.real(-2.0, 2.0);
            f.c[1] = io.real(-2.0, 2.0);
            f.c[2] = io.real(0.5, 4.0);
            f.c[3] = io.given(0.0);
            f.c[4] = io.given(0.0);
            f.c[5] = io.given(0.0);
        }
        else
        {
            for (int32_t i = 0; i < 6; ++i) { f.c[i] = io.real(-1.5, 1.5); }
        }
        double xMin = io.real(-2.0, -0.25);
        double xMax = io.real(0.25, 2.0);
        int32_t degree = io.integer(degreeLo, degreeHi);
        int32_t maxRemezIterations = io.integer(maxRemezLo, maxRemezHi);
        int32_t maxBisectionIterations = io.integer(4, 64);
        int32_t maxBracketIterations = io.integer(1, 8);

        RemezAlgorithm<double> remez;
        size_t iterations = remez.Execute(
            [&f](double const& x) { return f.F(x); },
            [&f](double const& x) { return f.FDer(x); },
            xMin, xMax, static_cast<size_t>(degree),
            static_cast<size_t>(maxRemezIterations),
            static_cast<size_t>(maxBisectionIterations),
            static_cast<size_t>(maxBracketIterations));
        OutIterations(io, iterations, std::numeric_limits<size_t>::max());
        io.outReal(remez.GetEstimatedMaxError());
        for (double value : remez.GetCoefficients()) { io.outReal(value); }
        for (double value : remez.GetXNodes()) { io.outReal(value); }
        for (double value : remez.GetErrors()) { io.outReal(value); }
    }
}

// Degree 1 has no libm dependence at all: its single Chebyshev cosine is at
// pi/2 and ComputeInitialXNodes overwrites it with an exact zero.
ORACLE_CASE("RemezAlgorithm.execute.degreeOne")
{
    RunRemez(io, -1, 1, 1, 1, 8);
}

ORACLE_CASE("RemezAlgorithm.execute.rational")
{
    RunRemez(io, 0, 2, 5, 1, 6);
}

ORACLE_CASE("RemezAlgorithm.execute.polynomial")
{
    RunRemez(io, 1, 2, 5, 1, 6);
}

// A single Remez iteration, which stops before the exchange and reports the
// interpolation data of the initial Chebyshev nodes.
ORACLE_CASE("RemezAlgorithm.execute.singleIteration")
{
    RunRemez(io, -1, 2, 5, 1, 1);
}

// Throw parity: Execute asserts xMin < xMax, degree > 0 and the three
// iteration counts > 0; GetXExtreme asserts that the partition subinterval
// brackets a root of E'.
ORACLE_CASE("RemezAlgorithm.execute.invalidInput")
{
    double xMin = io.lattice(-2, 1);
    double xMax = io.lattice(-1, 2);
    int32_t degree = io.integer(0, 3);
    int32_t maxRemezIterations = io.integer(0, 3);
    int32_t maxBisectionIterations = io.integer(0, 8);
    int32_t maxBracketIterations = io.integer(0, 3);
    double c0 = io.real(-2.0, 2.0);
    double c1 = io.real(-2.0, 2.0);
    double c2 = io.real(0.5, 4.0);

    RemezFunction f{};
    f.kind = 0;
    f.c[0] = c0;
    f.c[1] = c1;
    f.c[2] = c2;
    RemezAlgorithm<double> remez;
    size_t iterations = remez.Execute(
        [&f](double const& x) { return f.F(x); },
        [&f](double const& x) { return f.FDer(x); },
        xMin, xMax, static_cast<size_t>(degree),
        static_cast<size_t>(maxRemezIterations),
        static_cast<size_t>(maxBisectionIterations),
        static_cast<size_t>(maxBracketIterations));
    OutIterations(io, iterations, std::numeric_limits<size_t>::max());
    io.outReal(remez.GetEstimatedMaxError());
    for (double value : remez.GetCoefficients()) { io.outReal(value); }
}

// ================================================= CholeskyDecomposition

namespace
{
    // Fill a GMatrix from a row-major array.
    GMatrix<double> ToGMatrix(std::vector<double> const& A, int32_t n)
    {
        GMatrix<double> M(n, n);
        for (int32_t r = 0; r < n; ++r)
        {
            for (int32_t c = 0; c < n; ++c)
            {
                M(r, c) = A[static_cast<size_t>(c) + static_cast<size_t>(n) * r];
            }
        }
        return M;
    }

    void OutGMatrixEntries(oracle::Ctx& io, GMatrix<double> const& M)
    {
        for (int32_t r = 0; r < M.GetNumRows(); ++r)
        {
            for (int32_t c = 0; c < M.GetNumCols(); ++c) { io.outReal(M(r, c)); }
        }
    }
}

// Factor, then SolveLower and SolveUpper on the factored matrix. Factor
// modifies only the lower triangle, so the whole matrix is emitted.
ORACLE_CASE("CholeskyDecomposition.factorAndSolve")
{
    int32_t n = io.integer(1, 5);
    int32_t mode = io.integer(0, 5);
    std::vector<double> A = DrawSpd(io, n, mode);
    std::vector<double> B(static_cast<size_t>(n), 0.0);
    for (auto& value : B) { value = io.lattice(-4, 4); }

    GMatrix<double> M = ToGMatrix(A, n);
    CholeskyDecomposition<double> decomposer(n);
    bool success = decomposer.Factor(M);
    io.outBool(success);
    OutGMatrixEntries(io, M);
    if (success)
    {
        GVector<double> Y(n);
        for (int32_t i = 0; i < n; ++i) { Y[i] = B[i]; }
        decomposer.SolveLower(M, Y);
        for (int32_t i = 0; i < n; ++i) { io.outReal(Y[i]); }
        decomposer.SolveUpper(M, Y);
        for (int32_t i = 0; i < n; ++i) { io.outReal(Y[i]); }
    }
}

// The compile-time specialization CholeskyDecomposition<double, N>, which
// the port collapsed into the same run-time class. N is 3 here.
ORACLE_CASE("CholeskyDecomposition.factorAndSolve.fixedSize")
{
    int32_t mode = io.integer(0, 5);
    std::vector<double> A = DrawSpd(io, 3, mode);
    Vector<3, double> B{};
    for (int32_t i = 0; i < 3; ++i) { B[i] = io.lattice(-4, 4); }

    Matrix<3, 3, double> M{};
    for (int32_t r = 0; r < 3; ++r)
    {
        for (int32_t c = 0; c < 3; ++c) { M(r, c) = A[static_cast<size_t>(c) + 3 * r]; }
    }
    CholeskyDecomposition<double, 3> decomposer;
    bool success = decomposer.Factor(M);
    io.outBool(success);
    io.outMat(M);
    if (success)
    {
        Vector<3, double> Y = B;
        decomposer.SolveLower(M, Y);
        io.outVec(Y);
        decomposer.SolveUpper(M, Y);
        io.outVec(Y);
    }
}

// Throw parity: the run-time class raises LogError when the matrix is not
// N-by-N or the vector does not have N components.
ORACLE_CASE("CholeskyDecomposition.invalidSize")
{
    // Half the records use the valid shape, so both branches are covered.
    int32_t n = io.integer(1, 4);
    bool valid = io.boolean();
    int32_t rows = valid ? io.integer(n, n) : io.integer(1, 4);
    int32_t cols = valid ? io.integer(n, n) : io.integer(1, 4);
    int32_t vectorSize = valid ? io.integer(n, n) : io.integer(1, 4);
    std::vector<double> A(static_cast<size_t>(rows) * static_cast<size_t>(cols), 0.0);
    for (auto& value : A) { value = io.lattice(1, 4); }
    std::vector<double> B(static_cast<size_t>(vectorSize), 0.0);
    for (auto& value : B) { value = io.lattice(-3, 3); }

    GMatrix<double> M(rows, cols);
    for (int32_t r = 0; r < rows; ++r)
    {
        for (int32_t c = 0; c < cols; ++c)
        {
            M(r, c) = A[static_cast<size_t>(c) + static_cast<size_t>(cols) * r];
        }
    }
    CholeskyDecomposition<double> decomposer(n);
    bool success = decomposer.Factor(M);
    GVector<double> Y(vectorSize);
    for (int32_t i = 0; i < vectorSize; ++i) { Y[i] = B[i]; }
    decomposer.SolveLower(M, Y);
    std::vector<double> afterLower(static_cast<size_t>(vectorSize), 0.0);
    for (int32_t i = 0; i < vectorSize; ++i) { afterLower[i] = Y[i]; }
    decomposer.SolveUpper(M, Y);

    // Emitted only once every call that can throw has returned, so a record
    // either throws or produces its whole output.
    io.outBool(success);
    for (double value : afterLower) { io.outReal(value); }
    for (int32_t i = 0; i < vectorSize; ++i) { io.outReal(Y[i]); }
}

// ============================================ BlockCholeskyDecomposition

namespace
{
    // The compile-time specialization, which indexes inside a block with
    // block(row, col) and is the behavior the port implements.
    template <int32_t B, int32_t Nb>
    void RunBlockCholesky(oracle::Ctx& io, std::vector<double> const& A,
        std::vector<double> const& rhs)
    {
        int32_t const dim = B * Nb;
        BlockCholeskyDecomposition<double, B, Nb> decomposer;
        typename BlockCholeskyDecomposition<double, B, Nb>::BlockMatrix M{};
        for (int32_t r = 0; r < dim; ++r)
        {
            for (int32_t c = 0; c < dim; ++c)
            {
                decomposer.Set(M, r, c, A[static_cast<size_t>(c) + static_cast<size_t>(dim) * r]);
            }
        }
        // Get() over the whole matrix, before any modification.
        for (int32_t r = 0; r < dim; ++r)
        {
            for (int32_t c = 0; c < dim; ++c) { io.outReal(decomposer.Get(M, r, c)); }
        }

        bool success = decomposer.Factor(M);
        io.outBool(success);
        for (int32_t r = 0; r < dim; ++r)
        {
            for (int32_t c = 0; c < dim; ++c) { io.outReal(decomposer.Get(M, r, c)); }
        }

        if (success)
        {
            typename BlockCholeskyDecomposition<double, B, Nb>::BlockVector Y{};
            for (int32_t r = 0, k = 0; r < Nb; ++r)
            {
                for (int32_t i = 0; i < B; ++i, ++k) { Y[r][i] = rhs[k]; }
            }
            decomposer.SolveLower(M, Y);
            for (int32_t r = 0; r < Nb; ++r)
            {
                for (int32_t i = 0; i < B; ++i) { io.outReal(Y[r][i]); }
            }
            decomposer.SolveUpper(M, Y);
            for (int32_t r = 0; r < Nb; ++r)
            {
                for (int32_t i = 0; i < B; ++i) { io.outReal(Y[r][i]); }
            }
        }
    }

    // The run-time specialization BlockCholeskyDecomposition<Real,0,0>,
    // whose in-block scalar offset uses the block-level stride NumBlocks
    // (issue #209). It agrees with the compile-time one exactly when
    // BlockSize == NumBlocks; NumBlocks < BlockSize keeps every index in
    // range while making the two differ.
    void RunBlockCholeskyRuntime(oracle::Ctx& io, int32_t B, int32_t Nb,
        std::vector<double> const& A, std::vector<double> const& rhs)
    {
        int32_t const dim = B * Nb;
        BlockCholeskyDecomposition<double, 0, 0> decomposer(B, Nb);
        std::vector<GMatrix<double>> M(static_cast<size_t>(Nb) * Nb);
        for (auto& block : M) { block.SetSize(B, B); }
        for (int32_t r = 0; r < dim; ++r)
        {
            for (int32_t c = 0; c < dim; ++c)
            {
                decomposer.Set(M, r, c, A[static_cast<size_t>(c) + static_cast<size_t>(dim) * r]);
            }
        }

        bool success = decomposer.Factor(M);
        io.outBool(success);
        for (size_t b = 0; b < M.size(); ++b)
        {
            for (int32_t r = 0; r < B; ++r)
            {
                for (int32_t c = 0; c < B; ++c) { io.outReal(M[b](r, c)); }
            }
        }

        if (success)
        {
            std::vector<GVector<double>> Y(static_cast<size_t>(Nb));
            for (int32_t r = 0, k = 0; r < Nb; ++r)
            {
                Y[r].SetSize(B);
                for (int32_t i = 0; i < B; ++i, ++k) { Y[r][i] = rhs[k]; }
            }
            decomposer.SolveLower(M, Y);
            for (int32_t r = 0; r < Nb; ++r)
            {
                for (int32_t i = 0; i < B; ++i) { io.outReal(Y[r][i]); }
            }
            decomposer.SolveUpper(M, Y);
            for (int32_t r = 0; r < Nb; ++r)
            {
                for (int32_t i = 0; i < B; ++i) { io.outReal(Y[r][i]); }
            }
        }
    }
}

ORACLE_CASE("BlockCholeskyDecomposition.factorAndSolve")
{
    // (blockSize, numBlocks) pairs whose product is at most 6.
    static int32_t const shapes[][2] = {
        { 1, 2 }, { 1, 3 }, { 1, 4 }, { 2, 1 }, { 2, 2 }, { 2, 3 },
        { 3, 1 }, { 3, 2 }, { 1, 1 }
    };
    int32_t shape = io.integer(0, 8);
    int32_t B = shapes[shape][0], Nb = shapes[shape][1];
    int32_t dim = B * Nb;
    int32_t mode = io.integer(0, 5);
    std::vector<double> A = DrawSpd(io, dim, mode);
    std::vector<double> rhs(static_cast<size_t>(dim), 0.0);
    for (auto& value : rhs) { value = io.lattice(-4, 4); }

    if (B == 1 && Nb == 1) { RunBlockCholesky<1, 1>(io, A, rhs); }
    else if (B == 1 && Nb == 2) { RunBlockCholesky<1, 2>(io, A, rhs); }
    else if (B == 1 && Nb == 3) { RunBlockCholesky<1, 3>(io, A, rhs); }
    else if (B == 1) { RunBlockCholesky<1, 4>(io, A, rhs); }
    else if (B == 2 && Nb == 1) { RunBlockCholesky<2, 1>(io, A, rhs); }
    else if (B == 2 && Nb == 2) { RunBlockCholesky<2, 2>(io, A, rhs); }
    else if (B == 2) { RunBlockCholesky<2, 3>(io, A, rhs); }
    else if (Nb == 1) { RunBlockCholesky<3, 1>(io, A, rhs); }
    else { RunBlockCholesky<3, 2>(io, A, rhs); }
}

// Issue #209: the run-time BlockCholeskyDecomposition<Real,0,0> addresses
// scalars inside a BlockSize-by-BlockSize block with its block-level helper
// GetIndex(row,col) = col + row*NumBlocks, so SolveLower, SolveUpper,
// LowerTriangularSolver and SubtractiveUpdate read and write the wrong
// entries whenever BlockSize != NumBlocks. The port follows the
// compile-time specialization, which uses block(row,col). NumBlocks <
// BlockSize is chosen so that every one of upstream's indices is still
// inside the block's storage (no out-of-range read).
ORACLE_CASE("BlockCholeskyDecomposition.runtimeStrideDeviation")
{
    // NumBlocks >= 2 (otherwise the inner loops that carry the defect are
    // empty) and BlockSize > NumBlocks (otherwise the two indexings agree,
    // and BlockSize < NumBlocks would read out of the block's storage).
    static int32_t const shapes[][2] = { { 3, 2 }, { 4, 2 } };
    int32_t shape = io.integer(0, 1);
    int32_t B = shapes[shape][0], Nb = shapes[shape][1];
    int32_t dim = B * Nb;
    int32_t mode = io.integer(0, 1);
    std::vector<double> A = DrawSpd(io, dim, mode);
    std::vector<double> rhs(static_cast<size_t>(dim), 0.0);
    for (auto& value : rhs) { value = io.lattice(-4, 4); }

    RunBlockCholeskyRuntime(io, B, Nb, A, rhs);
}

// Throw parity for the run-time constructor: blockSize > 0 and numBlocks > 0.
ORACLE_CASE("BlockCholeskyDecomposition.invalidSize")
{
    int32_t B = io.integer(-1, 2);
    int32_t Nb = io.integer(-1, 2);
    double value = io.real(1.0, 4.0);

    BlockCholeskyDecomposition<double, 0, 0> decomposer(B, Nb);
    io.outInt(decomposer.NumDimensions);
    io.outReal(value);
}

// ==================================================== LDLTDecomposition

ORACLE_CASE("LDLTDecomposition.factorAndSolve")
{
    int32_t n = io.integer(1, 5);
    int32_t mode = io.integer(0, 5);
    std::vector<double> A = DrawSpd(io, n, mode);
    std::vector<double> B(static_cast<size_t>(n), 0.0);
    for (auto& value : B) { value = io.lattice(-4, 4); }

    GMatrix<double> M = ToGMatrix(A, n);
    GVector<double> rhs(n);
    for (int32_t i = 0; i < n; ++i) { rhs[i] = B[i]; }

    LDLTDecomposition<double> decomposer(n);
    GMatrix<double> L{}, D{};
    bool success = decomposer.Factor(M, L, D);
    io.outBool(success);
    OutGMatrixEntries(io, L);
    OutGMatrixEntries(io, D);
    if (success)
    {
        GVector<double> X(n);
        decomposer.Solve(L, D, rhs, X);
        for (int32_t i = 0; i < n; ++i) { io.outReal(X[i]); }
    }

    // The factoring overload of Solve.
    GVector<double> X2(n);
    bool success2 = decomposer.Solve(M, rhs, X2);
    io.outBool(success2);
    if (success2)
    {
        for (int32_t i = 0; i < n; ++i) { io.outReal(X2[i]); }
    }
}

// The compile-time specialization LDLTDecomposition<double, N>, N = 4.
ORACLE_CASE("LDLTDecomposition.factorAndSolve.fixedSize")
{
    int32_t mode = io.integer(0, 5);
    std::vector<double> A = DrawSpd(io, 4, mode);
    Vector<4, double> B{};
    for (int32_t i = 0; i < 4; ++i) { B[i] = io.lattice(-4, 4); }

    Matrix<4, 4, double> M{};
    for (int32_t r = 0; r < 4; ++r)
    {
        for (int32_t c = 0; c < 4; ++c) { M(r, c) = A[static_cast<size_t>(c) + 4 * r]; }
    }
    LDLTDecomposition<double, 4> decomposer;
    Matrix<4, 4, double> L{}, D{};
    bool success = decomposer.Factor(M, L, D);
    io.outBool(success);
    io.outMat(L);
    io.outMat(D);
    Vector<4, double> X{};
    bool success2 = decomposer.Solve(M, B, X);
    io.outBool(success2);
    if (success2) { io.outVec(X); }
}

// Throw parity: the run-time class asserts N > 0 and the matrix and vector
// sizes.
ORACLE_CASE("LDLTDecomposition.invalidSize")
{
    int32_t n = io.integer(0, 4);
    bool valid = io.boolean();
    int32_t rows = (valid && n > 0) ? io.integer(n, n) : io.integer(1, 4);
    int32_t vectorSize = (valid && n > 0) ? io.integer(n, n) : io.integer(1, 4);
    std::vector<double> A(static_cast<size_t>(rows) * static_cast<size_t>(rows), 0.0);
    for (auto& value : A) { value = io.lattice(1, 4); }
    std::vector<double> B(static_cast<size_t>(vectorSize), 0.0);
    for (auto& value : B) { value = io.lattice(-3, 3); }

    LDLTDecomposition<double> decomposer(n);
    GMatrix<double> M = ToGMatrix(A, rows);
    GVector<double> rhs(vectorSize);
    for (int32_t i = 0; i < vectorSize; ++i) { rhs[i] = B[i]; }
    GVector<double> X{};
    bool success = decomposer.Solve(M, rhs, X);
    io.outBool(success);
    // Upstream leaves X untouched when the factoring fails, so nothing but
    // the boolean is defined on that path.
    if (success)
    {
        for (int32_t i = 0; i < X.GetSize(); ++i) { io.outReal(X[i]); }
    }
}

// =============================================== BlockLDLTDecomposition

namespace
{
    void RunBlockLdlt(oracle::Ctx& io, int32_t B, int32_t Nb,
        std::vector<double> const& A, std::vector<double> const& rhs)
    {
        int32_t const dim = B * Nb;
        BlockLDLTDecomposition<double> decomposer(B, Nb);
        GMatrix<double> full = ToGMatrix(A, dim);
        GVector<double> fullRhs(dim);
        for (int32_t i = 0; i < dim; ++i) { fullRhs[i] = rhs[i]; }

        std::vector<GMatrix<double>> MBlock{};
        decomposer.Convert(full, MBlock);
        std::vector<GVector<double>> BBlock{};
        decomposer.Convert(fullRhs, BBlock);

        // Convert back and Get() over the whole matrix: a round trip that
        // must reproduce the input.
        GMatrix<double> roundTrip{};
        decomposer.Convert(MBlock, roundTrip);
        OutGMatrixEntries(io, roundTrip);
        for (int32_t r = 0; r < dim; ++r)
        {
            for (int32_t c = 0; c < dim; ++c)
            {
                double value = 0.0;
                decomposer.Get(MBlock, r, c, value);
                io.outReal(value);
            }
        }

        std::vector<GMatrix<double>> L{}, D{};
        bool success = decomposer.Factor(MBlock, L, D);
        io.outBool(success);
        for (auto const& block : L) { OutGMatrixEntries(io, block); }
        for (auto const& block : D) { OutGMatrixEntries(io, block); }

        std::vector<GVector<double>> X{};
        bool success2 = decomposer.Solve(MBlock, BBlock, X);
        io.outBool(success2);
        if (success2)
        {
            for (auto const& block : X)
            {
                for (int32_t i = 0; i < block.GetSize(); ++i) { io.outReal(block[i]); }
            }
        }
    }
}

ORACLE_CASE("BlockLDLTDecomposition.factorAndSolve")
{
    static int32_t const shapes[][2] = {
        { 1, 2 }, { 1, 3 }, { 2, 2 }, { 2, 3 }, { 3, 2 }, { 1, 1 }, { 2, 1 }, { 3, 3 }
    };
    int32_t shape = io.integer(0, 7);
    int32_t B = shapes[shape][0], Nb = shapes[shape][1];
    int32_t dim = B * Nb;
    int32_t mode = io.integer(0, 5);
    std::vector<double> A = DrawSpd(io, dim, mode);
    std::vector<double> rhs(static_cast<size_t>(dim), 0.0);
    for (auto& value : rhs) { value = io.lattice(-4, 4); }

    RunBlockLdlt(io, B, Nb, A, rhs);
}

// Set() over the whole matrix, then Get() to read it back, plus the vector
// conversions with BlockSize == NumBlocks (where upstream's Convert
// assertion happens to accept the input).
ORACLE_CASE("BlockLDLTDecomposition.getSetConvert")
{
    int32_t B = io.integer(1, 3);
    // BlockSize == NumBlocks, the one shape on which upstream's Convert
    // assertion accepts a block vector.
    int32_t Nb = static_cast<int32_t>(io.given(static_cast<double>(B)));
    int32_t dim = B * Nb;
    std::vector<double> A(static_cast<size_t>(dim) * static_cast<size_t>(dim), 0.0);
    for (auto& value : A) { value = io.lattice(-5, 5); }
    std::vector<double> v(static_cast<size_t>(dim), 0.0);
    for (auto& value : v) { value = io.lattice(-5, 5); }

    BlockLDLTDecomposition<double> decomposer(B, Nb);
    std::vector<GMatrix<double>> MBlock(static_cast<size_t>(Nb) * Nb);
    for (auto& block : MBlock) { block.SetSize(B, B); block.MakeZero(); }
    for (int32_t r = 0; r < dim; ++r)
    {
        for (int32_t c = 0; c < dim; ++c)
        {
            decomposer.Set(MBlock, r, c, A[static_cast<size_t>(c) + static_cast<size_t>(dim) * r]);
        }
    }
    for (int32_t r = 0; r < dim; ++r)
    {
        for (int32_t c = 0; c < dim; ++c)
        {
            double value = 0.0;
            decomposer.Get(MBlock, r, c, value);
            io.outReal(value);
        }
    }

    GVector<double> full(dim);
    for (int32_t i = 0; i < dim; ++i) { full[i] = v[i]; }
    std::vector<GVector<double>> VBlock{};
    decomposer.Convert(full, VBlock);
    GVector<double> back{};
    decomposer.Convert(VBlock, back);
    for (int32_t i = 0; i < back.GetSize(); ++i) { io.outReal(back[i]); }
}

// Issue #209: BlockLDLTDecomposition<T>::Convert(BlockVector, GVector)
// verifies each block vector with 'current.GetSize() == NumBlocks', but a
// block vector has BlockSize components, so upstream raises LogAssert for
// every valid input with BlockSize != NumBlocks. The port checks blockSize
// and converts.
ORACLE_CASE("BlockLDLTDecomposition.convertBlockToVector.deviation")
{
    int32_t B = io.integer(1, 4);
    // NumBlocks != BlockSize, which is where upstream's assertion rejects a
    // valid block vector.
    int32_t raw = io.rawInteger(1, 3);
    int32_t Nb = static_cast<int32_t>(io.given(static_cast<double>(raw >= B ? raw + 1 : raw)));
    int32_t dim = B * Nb;
    std::vector<double> v(static_cast<size_t>(dim), 0.0);
    for (auto& value : v) { value = io.lattice(-5, 5); }

    BlockLDLTDecomposition<double> decomposer(B, Nb);
    GVector<double> full(dim);
    for (int32_t i = 0; i < dim; ++i) { full[i] = v[i]; }
    std::vector<GVector<double>> VBlock{};
    decomposer.Convert(full, VBlock);
    GVector<double> back{};
    decomposer.Convert(VBlock, back);
    io.outInt(back.GetSize());
    for (int32_t i = 0; i < back.GetSize(); ++i) { io.outReal(back[i]); }
}

// The compile-time specialization BlockLDLTDecomposition<double, B, Nb>,
// B = 2, Nb = 2, including its Convert overloads.
ORACLE_CASE("BlockLDLTDecomposition.factorAndSolve.fixedSize")
{
    int32_t mode = io.integer(0, 5);
    std::vector<double> A = DrawSpd(io, 4, mode);
    std::vector<double> rhs(4, 0.0);
    for (auto& value : rhs) { value = io.lattice(-4, 4); }

    Matrix<4, 4, double> full{};
    for (int32_t r = 0; r < 4; ++r)
    {
        for (int32_t c = 0; c < 4; ++c) { full(r, c) = A[static_cast<size_t>(c) + 4 * r]; }
    }
    Vector<4, double> fullRhs{};
    for (int32_t i = 0; i < 4; ++i) { fullRhs[i] = rhs[i]; }

    BlockLDLTDecomposition<double, 2, 2> decomposer;
    BlockLDLTDecomposition<double, 2, 2>::BlockMatrix MBlock{};
    decomposer.Convert(full, MBlock);
    BlockLDLTDecomposition<double, 2, 2>::BlockVector BBlock{};
    decomposer.Convert(fullRhs, BBlock);

    Matrix<4, 4, double> roundTrip{};
    decomposer.Convert(MBlock, roundTrip);
    io.outMat(roundTrip);

    BlockLDLTDecomposition<double, 2, 2>::BlockMatrix L{}, D{};
    bool success = decomposer.Factor(MBlock, L, D);
    io.outBool(success);
    for (int32_t r = 0; r < 2; ++r)
    {
        for (int32_t c = 0; c < 2; ++c) { io.outMat(L[r][c]); }
    }
    for (int32_t r = 0; r < 2; ++r)
    {
        for (int32_t c = 0; c < 2; ++c) { io.outMat(D[r][c]); }
    }

    BlockLDLTDecomposition<double, 2, 2>::BlockVector X{};
    bool success2 = decomposer.Solve(MBlock, BBlock, X);
    io.outBool(success2);
    if (success2)
    {
        Vector<4, double> flat{};
        decomposer.Convert(X, flat);
        io.outVec(flat);
    }
}

// Throw parity for the run-time constructor: blockSize > 0, numBlocks > 0.
ORACLE_CASE("BlockLDLTDecomposition.invalidSize")
{
    int32_t B = io.integer(-1, 2);
    int32_t Nb = io.integer(-1, 2);
    double value = io.real(1.0, 4.0);

    BlockLDLTDecomposition<double> decomposer(B, Nb);
    io.outInt(decomposer.NumDimensions);
    io.outReal(value);
}

// =========================================================== MinimizeN

ORACLE_CASE("MinimizeN.getMinimum")
{
    // Accept only inputs on which neither of the port's two deliberate
    // fixes changes the trajectory: the Minimize1 bracket fix (#298) must
    // not fire in any line search, and dropping upstream's stray
    // 'mDConjIndex = 0' (#146) must not change the result.
    MinNDraw d{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        d = RawMinN(io, (io.index() + attempt) % 4, 1, 6);
        accepted = !MinNDiffers(d);
    }
    if (!accepted) { d = SoundMinN(); }
    RecordMinN(io, d);

    ObjectiveN f = d.f;
    MinimizeN<double> minimizer(d.f.dimensions,
        [&f](double const* x) { return EvaluateN(f, x); },
        d.maxLevel, d.maxBracket, d.maxIterations, d.epsilon);
    std::vector<double> tMin(d.f.dimensions, 0.0);
    double fMin = 0.0;
    minimizer.GetMinimum(d.t0.data(), d.t1.data(), d.tInitial.data(), tMin.data(), fMin);
    for (double value : tMin) { io.outReal(value); }
    io.outReal(fMin);
    io.outReal(minimizer.GetEpsilon());
}

// Issue #146: upstream's stray 'mDConjIndex = 0' writes the conjugate
// direction into slot 0, where the cycling shift immediately overwrites it,
// and slot 'dimensions' keeps the first iteration's conjugate direction
// forever. The direction set degenerates and Powell's method stalls. The
// records kept here are exactly those on which the assignment changes the
// reported minimum.
ORACLE_CASE("MinimizeN.getMinimum.deviation")
{
    MinNDraw d{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        d = RawMinN(io, (io.index() + attempt) % 4, 3, 10);
        accepted = MinNDiffers(d);
    }
    RecordMinN(io, d);

    ObjectiveN f = d.f;
    MinimizeN<double> minimizer(d.f.dimensions,
        [&f](double const* x) { return EvaluateN(f, x); },
        d.maxLevel, d.maxBracket, d.maxIterations, d.epsilon);
    std::vector<double> tMin(d.f.dimensions, 0.0);
    double fMin = 0.0;
    minimizer.GetMinimum(d.t0.data(), d.t1.data(), d.tInitial.data(), tMin.data(), fMin);
    for (double value : tMin) { io.outReal(value); }
    io.outReal(fMin);
}

// A second GetMinimum call on the same object. Upstream never resets
// mDConjIndex, so the second call starts with whatever the first left
// behind (issue #146, second part); the port's mDConjIndex is invariant, so
// the second call behaves like the first. The records kept here are those
// where the two disagree.
ORACLE_CASE("MinimizeN.getMinimum.reuseDeviation")
{
    MinNDraw d{};
    std::vector<double> second;
    double shift = 0.0;
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        d = RawMinN(io, (io.index() + attempt) % 4, 2, 8);
        shift = io.raw(0.05, 0.4);
        second.assign(d.f.dimensions, 0.0);
        for (int32_t i = 0; i < d.f.dimensions; ++i)
        {
            second[i] = d.tInitial[i] + shift * (d.t1[i] - d.tInitial[i]);
        }
        accepted = MinNReuseDiffers(d, second);
    }
    RecordMinN(io, d);
    io.given(shift);

    ObjectiveN f = d.f;
    MinimizeN<double> minimizer(d.f.dimensions,
        [&f](double const* x) { return EvaluateN(f, x); },
        d.maxLevel, d.maxBracket, d.maxIterations, d.epsilon);
    std::vector<double> tMin(d.f.dimensions, 0.0);
    double fMin = 0.0;
    minimizer.GetMinimum(d.t0.data(), d.t1.data(), d.tInitial.data(), tMin.data(), fMin);

    std::vector<double> tMin2(d.f.dimensions, 0.0);
    double fMin2 = 0.0;
    minimizer.GetMinimum(d.t0.data(), d.t1.data(), second.data(), tMin2.data(), fMin2);
    for (double value : tMin2) { io.outReal(value); }
    io.outReal(fMin2);
}

// SetEpsilon clamps a nonpositive argument to zero, which turns off the
// "the position did not change" break and lets every iteration run.
ORACLE_CASE("MinimizeN.setEpsilon")
{
    MinNDraw d{};
    bool accepted = false;
    for (int32_t attempt = 0; attempt < 32 && !accepted; ++attempt)
    {
        d = RawMinN(io, (io.index() + attempt) % 4, 1, 2);
        // The recorded 'epsilon' is the argument of SetEpsilon here, so the
        // probe (which clamps it exactly as MinimizeN does) sees the value
        // the search will really use.
        d.epsilon = io.raw(-1.0, 1.0);
        accepted = !MinNDiffers(d);
    }
    if (!accepted) { d = SoundMinN(); }
    RecordMinN(io, d);

    ObjectiveN f = d.f;
    MinimizeN<double> minimizer(d.f.dimensions,
        [&f](double const* x) { return EvaluateN(f, x); },
        d.maxLevel, d.maxBracket, d.maxIterations, -1.0);
    io.outReal(minimizer.GetEpsilon());
    minimizer.SetEpsilon(d.epsilon);
    io.outReal(minimizer.GetEpsilon());
    std::vector<double> tMin(d.f.dimensions, 0.0);
    double fMin = 0.0;
    minimizer.GetMinimum(d.t0.data(), d.t1.data(), d.tInitial.data(), tMin.data(), fMin);
    for (double value : tMin) { io.outReal(value); }
    io.outReal(fMin);
}
