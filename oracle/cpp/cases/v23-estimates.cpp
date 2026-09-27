// Verify group 23 (estimates): differential cases for Exp2Estimate.h,
// ExpEstimate.h, Log2Estimate.h, LogEstimate.h, SinEstimate.h,
// CosEstimate.h, TanEstimate.h, ATanEstimate.h, ACosEstimate.h,
// ASinEstimate.h, SqrtEstimate.h, InvSqrtEstimate.h, ChebyshevRatio.h,
// ChebyshevRatioEstimate.h, Slerp.h and SlerpEstimate.h.
//
// Instantiation: T = double, the only instantiation the port implements (it
// takes the degree as a runtime argument). Every record evaluates EVERY
// degree the header's static_assert admits, through a fold expression over an
// index sequence, so every template instantiation is compared on every record.
//
// libm. The polynomial evaluations are + - * and std::sqrt; the range
// reductions call std::floor, std::frexp, std::ldexp, std::remainder and
// std::fabs, all exactly specified by IEEE 754 / C99. Every estimate case is
// therefore { exact: true }. Only ChebyshevRatio.h and Slerp.h call libm
// (std::sin, std::acos); those cases carry a tolerance, and their exactly
// evaluated branches (angle 0, cosAngle >= 1) are emitted where the replay
// holds them to bit identity with outRealExact.
//
// Scalar generators. Each record draws one value from EACH of eight modes
// (uniform, dyadic lattice, domain endpoints and their neighbours, values
// approaching an endpoint, outside the domain, a "wild" population of signed
// zeros, subnormals, DBL_MIN, DBL_MAX, infinities and NaN), so every
// committed record exercises every population. Each mode records exactly one
// double per sample.
#define ORACLE_FAMILY "v23-estimates"
#include "Oracle.h"

#include <Mathematics/ACosEstimate.h>
#include <Mathematics/ASinEstimate.h>
#include <Mathematics/ATanEstimate.h>
#include <Mathematics/ChebyshevRatio.h>
#include <Mathematics/ChebyshevRatioEstimate.h>
#include <Mathematics/CosEstimate.h>
#include <Mathematics/Exp2Estimate.h>
#include <Mathematics/ExpEstimate.h>
#include <Mathematics/InvSqrtEstimate.h>
#include <Mathematics/Log2Estimate.h>
#include <Mathematics/LogEstimate.h>
#include <Mathematics/SinEstimate.h>
#include <Mathematics/Slerp.h>
#include <Mathematics/SlerpEstimate.h>
#include <Mathematics/SqrtEstimate.h>
#include <Mathematics/TanEstimate.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>

using namespace gte;

namespace
{
    int32_t constexpr NUM_MODES = 8;
    double const INF = std::numeric_limits<double>::infinity();
    double const QNAN = std::numeric_limits<double>::quiet_NaN();
    double const MIN_SUB = std::numeric_limits<double>::denorm_min();
    double const MIN_NORM = std::numeric_limits<double>::min();
    double const MAX_DBL = std::numeric_limits<double>::max();

    // x moved by k ulps (k < 0 moves down).
    double Ulps(double x, int32_t k)
    {
        for (; k > 0; --k) { x = std::nextafter(x, INF); }
        for (; k < 0; ++k) { x = std::nextafter(x, -INF); }
        return x;
    }

    double RawSign(oracle::Ctx& io)
    {
        return io.rawInteger(0, 1) == 0 ? 1.0 : -1.0;
    }

    // Signed zeros, subnormals, the normal limits, infinities, NaN, small
    // integers and ordinary values. Unrecorded; the caller records it.
    double WildRaw(oracle::Ctx& io)
    {
        double s = RawSign(io);
        switch (io.rawInteger(0, 11))
        {
        case 0:  return s * 0.0;
        case 1:  return s * MIN_SUB;
        case 2:  return s * std::ldexp(io.raw(1.0, 2.0), io.rawInteger(-1074, -1023));
        case 3:  return s * MIN_NORM;
        case 4:  return s * MAX_DBL;
        case 5:  return s * INF;
        case 6:  return QNAN;
        case 7:  return s * std::ldexp(io.raw(1.0, 2.0), io.rawInteger(-1022, 1023));
        case 8:  return s * static_cast<double>(io.rawInteger(0, 8));
        case 9:  return s * std::ldexp(1.0, io.rawInteger(-1074, 1023));
        case 10: return s * Ulps(1.0, io.rawInteger(-2, 2));
        default: return s * io.raw(0.0, 4.0);
        }
    }

    // A value in the closed interval [a, b] or around it, by mode:
    //   0, 7 uniform in [a, b]
    //   1    dyadic a + (b - a) * k/64, k in [0, 64]
    //   2    the endpoints and their neighbours within 2 ulps, and -0 when
    //        a is 0
    //   3    approaching a: a + (b - a) * 2^-k * u, k in [1, 1074]
    //   4    approaching b: b - (b - a) * 2^-k * u, k in [1, 60]
    //   5    uniform in [a - (b - a)/2, b + (b - a)/2] (outside the domain)
    //   6    wild
    double IntervalRaw(oracle::Ctx& io, int32_t mode, double a, double b)
    {
        double w = b - a;
        switch (mode)
        {
        case 1:
            return a + w * (static_cast<double>(io.rawInteger(0, 64)) / 64.0);
        case 2:
        {
            int32_t which = io.rawInteger(0, 2);
            int32_t k = io.rawInteger(-2, 2);
            if (which == 2) { return a == 0.0 ? -0.0 : Ulps(a, k); }
            return Ulps(which == 0 ? a : b, k);
        }
        case 3:
            return a + w * std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(1, 1074));
        case 4:
            return b - w * std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(1, 60));
        case 5:
            return io.raw(a - 0.5 * w, b + 0.5 * w);
        case 6:
            return WildRaw(io);
        default:
            return io.raw(a, b);
        }
    }

    // A value in [-L, L] or around it, by mode:
    //   0, 7 uniform in [-L, L]
    //   1    lattice k*L/16, k in [-16, 16]
    //   2    +-L and their neighbours within 2 ulps, +-0
    //   3    tiny: +-2^-k * u, k in [1, 1074]
    //   4    uniform in [-2L, 2L] (outside the domain)
    //   5    wild
    //   6    near the endpoints: +-(L - L * 2^-k * u), k in [1, 60]
    double SymmetricRaw(oracle::Ctx& io, int32_t mode, double L)
    {
        double s = RawSign(io);
        switch (mode)
        {
        case 1:
            return L * (static_cast<double>(io.rawInteger(-16, 16)) / 16.0);
        case 2:
        {
            int32_t which = io.rawInteger(0, 1);
            int32_t k = io.rawInteger(-2, 2);
            return which == 0 ? s * 0.0 : s * Ulps(L, k);
        }
        case 3:
            return s * std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(1, 1074));
        case 4:
            return io.raw(-2.0 * L, 2.0 * L);
        case 5:
            return WildRaw(io);
        case 6:
            return s * (L - L * std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(1, 60)));
        default:
            return io.raw(-L, L);
        }
    }
}

namespace
{
    // An argument for the periodic range reductions (period P, branch
    // thresholds at multiples of Q), by mode:
    //   0    uniform in [-4P, 4P]
    //   1    k*Q, k in [-64, 64], moved by up to 2 ulps: straddles the
    //        r > halfPi / r < -halfPi / |y| <= quarterPi comparisons
    //   2    exact remainder ties odd*(P/2), odd in {+-1, +-3, +-5, +-7}
    //        (exact products: the significand of pi ends in three zero
    //        bits), moved by at most 1 ulp
    //   3    huge: +-u*2^k, k in [0, 1023]
    //   4    +-integers below 2^31 and +-2^k exactly, k in [0, 1023]
    //   5    wild
    //   6    uniform in [-1e6, 1e6]
    //   7    uniform in [-P/2 - 1/2, P/2 + 1/2]
    double PeriodicRaw(oracle::Ctx& io, int32_t mode, double P, double Q)
    {
        double s = RawSign(io);
        switch (mode)
        {
        case 1:
        {
            double x = static_cast<double>(io.rawInteger(-64, 64)) * Q;
            return Ulps(x, io.rawInteger(-2, 2));
        }
        case 2:
        {
            double odd = static_cast<double>(2 * io.rawInteger(0, 3) + 1);
            double x = s * odd * (0.5 * P);
            return Ulps(x, io.rawInteger(-1, 1));
        }
        case 3:
            return s * std::ldexp(io.raw(1.0, 2.0), io.rawInteger(0, 1023));
        case 4:
            if (io.rawInteger(0, 1) == 0)
            {
                return s * static_cast<double>(io.rawInteger(0, 2147483647));
            }
            return s * std::ldexp(1.0, io.rawInteger(0, 1023));
        case 5:
            return WildRaw(io);
        case 6:
            return io.raw(-1e6, 1e6);
        case 7:
            return io.raw(-0.5 * P - 0.5, 0.5 * P + 0.5);
        default:
            return io.raw(-4.0 * P, 4.0 * P);
        }
    }

    // An argument for Exp2EstimateRR (scale 1) or ExpEstimateRR (scale ln 2),
    // x = scale * g, by mode of g:
    //   0    uniform in [-10, 10]
    //   1    integers in [-1080, 1030] moved by up to 2 ulps (y = x - floor(x)
    //        is 0 or 1 - ulp; 2^k boundaries, subnormal and overflowing
    //        results)
    //   2    uniform in [-1080, -1015] (subnormal results, underflow to 0)
    //   3    uniform in [1015, 1030] (overflow)
    //   4    +-u*2^k, k in [-1074, 29]
    //   5    specials: +-0, +-subnormals, +-DBL_MIN, +-inf, NaN,
    //        +-(2^30 - 1), +-(2^30 - 1/2)
    //   6    dyadic k/64 in [-4, 4]
    //   7    uniform in [-100, 100]
    // |g| stays below 2^30: static_cast<int32_t> of floor(x) outside the
    // int32_t range is undefined behaviour (MSVC produces INT_MIN, so
    // Exp2EstimateRR(3e9) returns 0 instead of +inf); the deviation case
    // Exp2Estimate.estimateRR.hugeArgument covers that range. For
    // ExpEstimateRR the bound keeps x*(1/ln 2) far inside the range as well.
    double ExpRRRaw(oracle::Ctx& io, int32_t mode, double scale)
    {
        double s = RawSign(io);
        double g = 0.0;
        switch (mode)
        {
        case 1:
            g = Ulps(static_cast<double>(io.rawInteger(-1080, 1030)), io.rawInteger(-2, 2));
            break;
        case 2:
            g = io.raw(-1080.0, -1015.0);
            break;
        case 3:
            g = io.raw(1015.0, 1030.0);
            break;
        case 4:
            g = s * std::ldexp(io.raw(1.0, 2.0), io.rawInteger(-1074, 29));
            break;
        case 5:
            switch (io.rawInteger(0, 7))
            {
            case 0:  g = s * 0.0; break;
            case 1:  g = s * MIN_SUB; break;
            case 2:  g = s * std::ldexp(io.raw(1.0, 2.0), io.rawInteger(-1074, -1023)); break;
            case 3:  g = s * MIN_NORM; break;
            case 4:  g = s * INF; break;
            case 5:  g = QNAN; break;
            case 6:  g = s * (1073741824.0 - 1.0); break;
            default: g = s * (1073741824.0 - 0.5); break;
            }
            break;
        case 6:
            g = static_cast<double>(io.rawInteger(-256, 256)) / 64.0;
            break;
        case 7:
            g = io.raw(-100.0, 100.0);
            break;
        default:
            g = io.raw(-10.0, 10.0);
            break;
        }
        return scale == 1.0 ? g : g * scale;
    }

    // A positive argument (Log2/Log/Sqrt/InvSqrt EstimateRR), by mode:
    //   0    u*2^k, k in [-1074, 1023] (log-uniform over all doubles)
    //   1    2^k, k in [-1074, 1023], moved by up to 2 ulps (even and odd
    //        exponents, the frexp boundaries)
    //   2    uniform in [0, 4]
    //   3    subnormal: u*2^-1022, u in [0, 1)
    //   4    1 +- u*2^-k, k in [1, 60]
    //   5    specials outside the domain: +-0 (SqrtEstimateRR(0), issue #57),
    //        MIN_SUB, DBL_MAX, +-inf, NaN, negative values
    //   6    uniform in [1, 2]
    //   7    u*2^k, k in [-30, 30]
    double PositiveRaw(oracle::Ctx& io, int32_t mode)
    {
        double s = RawSign(io);
        switch (mode)
        {
        case 1:
        {
            double x = std::ldexp(1.0, io.rawInteger(-1074, 1023));
            return Ulps(x, io.rawInteger(-2, 2));
        }
        case 2:
            return io.raw(0.0, 4.0);
        case 3:
            return std::ldexp(io.raw(0.0, 1.0), -1022);
        case 4:
            return 1.0 + s * std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(1, 60));
        case 5:
            switch (io.rawInteger(0, 6))
            {
            case 0:  return s * 0.0;
            case 1:  return MIN_SUB;
            case 2:  return MAX_DBL;
            case 3:  return s * INF;
            case 4:  return QNAN;
            default: return -std::ldexp(io.raw(1.0, 2.0), io.rawInteger(-1074, 1023));
            }
        case 6:
            return io.raw(1.0, 2.0);
        case 7:
            return std::ldexp(io.raw(1.0, 2.0), io.rawInteger(-30, 30));
        default:
            return std::ldexp(io.raw(1.0, 2.0), io.rawInteger(-1074, 1023));
        }
    }

    // An argument for ATanEstimateRR, by mode:
    //   0    uniform in [-4, 4]
    //   1    +-1 moved by up to 2 ulps (the |x| <= 1 comparison)
    //   2    +-u*2^k, k in [-1074, 1023] (1/x subnormal for the huge ones)
    //   3    wild
    //   4    k/8, k in [-32, 32]
    //   5    uniform in [-1, 1]
    //   6    uniform in [-100, 100]
    //   7    64/k, k in [-64, 64] \ {0}
    double ATanRRRaw(oracle::Ctx& io, int32_t mode)
    {
        double s = RawSign(io);
        switch (mode)
        {
        case 1:
            return s * Ulps(1.0, io.rawInteger(-2, 2));
        case 2:
            return s * std::ldexp(io.raw(1.0, 2.0), io.rawInteger(-1074, 1023));
        case 3:
            return WildRaw(io);
        case 4:
            return static_cast<double>(io.rawInteger(-32, 32)) / 8.0;
        case 5:
            return io.raw(-1.0, 1.0);
        case 6:
            return io.raw(-100.0, 100.0);
        case 7:
            return s * 64.0 / static_cast<double>(io.rawInteger(1, 64));
        default:
            return io.raw(-4.0, 4.0);
        }
    }
}

// ---------------------------------------------------------------------------
// Every degree of every scalar estimate. ALL_DEGREES(NAME, CALL) defines
// NAME(io, x, std::index_sequence<I...>) that emits CALL for every I in
// order (a comma fold evaluates left to right); CALL is parenthesized so the
// comma of the template argument list survives the macro.
// ---------------------------------------------------------------------------
#define ALL_DEGREES(NAME, CALL) \
    template <size_t... I> \
    void NAME(oracle::Ctx& io, double x, std::index_sequence<I...>) \
    { (io.outReal(CALL), ...); }

#define ALL_MAX_ERRORS(NAME, CALL) \
    template <size_t... I> \
    void NAME(oracle::Ctx& io, std::index_sequence<I...>) \
    { (io.outReal(CALL), ...); }

namespace
{
    // Degrees 1..7.
    ALL_DEGREES(AllExp2, (Exp2Estimate<double, I + 1>(x)))
    ALL_DEGREES(AllExp2RR, (Exp2EstimateRR<double, I + 1>(x)))
    ALL_MAX_ERRORS(AllExp2Err, (GetExp2EstimateMaxError<double, I + 1>()))
    ALL_DEGREES(AllExp, (ExpEstimate<double, I + 1>(x)))
    ALL_DEGREES(AllExpRR, (ExpEstimateRR<double, I + 1>(x)))
    ALL_MAX_ERRORS(AllExpErr, (GetExpEstimateMaxError<double, I + 1>()))
    // Degrees 1..8.
    ALL_DEGREES(AllLog2, (Log2Estimate<double, I + 1>(x)))
    ALL_DEGREES(AllLog2RR, (Log2EstimateRR<double, I + 1>(x)))
    ALL_MAX_ERRORS(AllLog2Err, (GetLog2EstimateMaxError<double, I + 1>()))
    ALL_DEGREES(AllLog, (LogEstimate<double, I + 1>(x)))
    ALL_DEGREES(AllLogRR, (LogEstimateRR<double, I + 1>(x)))
    ALL_MAX_ERRORS(AllLogErr, (GetLogEstimateMaxError<double, I + 1>()))
    // Degrees 3, 5, ..., 11.
    ALL_DEGREES(AllSin, (SinEstimate<double, 2 * I + 3>(x)))
    ALL_DEGREES(AllSinRR, (SinEstimateRR<double, 2 * I + 3>(x)))
    ALL_MAX_ERRORS(AllSinErr, (GetSinEstimateMaxError<double, 2 * I + 3>()))
    // Degrees 2, 4, ..., 10.
    ALL_DEGREES(AllCos, (CosEstimate<double, 2 * I + 2>(x)))
    ALL_DEGREES(AllCosRR, (CosEstimateRR<double, 2 * I + 2>(x)))
    ALL_MAX_ERRORS(AllCosErr, (GetCosEstimateMaxError<double, 2 * I + 2>()))
    // Degrees 3, 5, ..., 13.
    ALL_DEGREES(AllTan, (TanEstimate<double, 2 * I + 3>(x)))
    ALL_DEGREES(AllTanRR, (TanEstimateRR<double, 2 * I + 3>(x)))
    ALL_MAX_ERRORS(AllTanErr, (GetTanEstimateMaxError<double, 2 * I + 3>()))
    ALL_DEGREES(AllATan, (ATanEstimate<double, 2 * I + 3>(x)))
    ALL_DEGREES(AllATanRR, (ATanEstimateRR<double, 2 * I + 3>(x)))
    ALL_MAX_ERRORS(AllATanErr, (GetATanEstimateMaxError<double, 2 * I + 3>()))
    // Degrees 1..8.
    ALL_DEGREES(AllACos, (ACosEstimate<double, I + 1>(x)))
    ALL_MAX_ERRORS(AllACosErr, (GetACosEstimateMaxError<double, I + 1>()))
    ALL_DEGREES(AllASin, (ASinEstimate<double, I + 1>(x)))
    ALL_MAX_ERRORS(AllASinErr, (GetASinEstimateMaxError<double, I + 1>()))
    ALL_DEGREES(AllSqrt, (SqrtEstimate<double, I + 1>(x)))
    ALL_DEGREES(AllSqrtRR, (SqrtEstimateRR<double, I + 1>(x)))
    ALL_MAX_ERRORS(AllSqrtErr, (GetSqrtEstimateMaxError<double, I + 1>()))
    ALL_DEGREES(AllInvSqrt, (InvSqrtEstimate<double, I + 1>(x)))
    ALL_DEGREES(AllInvSqrtRR, (InvSqrtEstimateRR<double, I + 1>(x)))
    ALL_MAX_ERRORS(AllInvSqrtErr, (GetInvSqrtEstimateMaxError<double, I + 1>()))
    // Degrees 1..16 and 1..12.
    ALL_MAX_ERRORS(AllChbErr, (GetChebyshevRatioEstimateMaxError<double, I + 1>()))
    ALL_MAX_ERRORS(AllChbErrR, (GetChebyshevRatioEstimateRMaxError<double, I + 1>()))

    using Seq5 = std::make_index_sequence<5>;
    using Seq6 = std::make_index_sequence<6>;
    using Seq7 = std::make_index_sequence<7>;
    using Seq8 = std::make_index_sequence<8>;
}

// One case per entry point. Each record draws one argument per generator
// mode (NUM_MODES recorded doubles) and emits the estimate at every degree
// for each of them.
#define ESTIMATE_CASE(CASENAME, GENERATE, ALL, SEQ) \
    ORACLE_CASE(CASENAME) \
    { \
        for (int32_t mode = 0; mode < NUM_MODES; ++mode) \
        { \
            double x = io.given(GENERATE); \
            ALL(io, x, SEQ{}); \
        } \
    }

#define MAX_ERROR_CASE(CASENAME, ALL, SEQ) \
    ORACLE_CASE(CASENAME) \
    { \
        ALL(io, SEQ{}); \
    }

// Exp2Estimate.h: x in [0, 1]; Exp2EstimateRR: any x (below 2^30 in
// magnitude, see ExpRRRaw).
ESTIMATE_CASE("Exp2Estimate.estimate", IntervalRaw(io, mode, 0.0, 1.0), AllExp2, Seq7)
ESTIMATE_CASE("Exp2Estimate.estimateRR", ExpRRRaw(io, mode, 1.0), AllExp2RR, Seq7)
MAX_ERROR_CASE("Exp2Estimate.getMaxError", AllExp2Err, Seq7)

// ExpEstimate.h: x in [0, ln 2]; ExpEstimateRR: any x.
ESTIMATE_CASE("ExpEstimate.estimate", IntervalRaw(io, mode, 0.0, GTE_C_LN_2), AllExp, Seq7)
ESTIMATE_CASE("ExpEstimate.estimateRR", ExpRRRaw(io, mode, GTE_C_LN_2), AllExpRR, Seq7)
MAX_ERROR_CASE("ExpEstimate.getMaxError", AllExpErr, Seq7)

// Log2Estimate.h: x in [1, 2]; Log2EstimateRR: x > 0 (frexp).
ESTIMATE_CASE("Log2Estimate.estimate", IntervalRaw(io, mode, 1.0, 2.0), AllLog2, Seq8)
ESTIMATE_CASE("Log2Estimate.estimateRR", PositiveRaw(io, mode), AllLog2RR, Seq8)
MAX_ERROR_CASE("Log2Estimate.getMaxError", AllLog2Err, Seq8)

// LogEstimate.h. GetLogEstimateMaxError returns the log2 bound (issue #57,
// preserved).
ESTIMATE_CASE("LogEstimate.estimate", IntervalRaw(io, mode, 1.0, 2.0), AllLog, Seq8)
ESTIMATE_CASE("LogEstimate.estimateRR", PositiveRaw(io, mode), AllLogRR, Seq8)
MAX_ERROR_CASE("LogEstimate.getMaxError", AllLogErr, Seq8)

// SinEstimate.h: x in [-pi/2, pi/2]; SinEstimateRR: std::remainder(x, 2pi).
ESTIMATE_CASE("SinEstimate.estimate", SymmetricRaw(io, mode, GTE_C_HALF_PI), AllSin, Seq5)
ESTIMATE_CASE("SinEstimate.estimateRR", PeriodicRaw(io, mode, GTE_C_TWO_PI, GTE_C_HALF_PI), AllSinRR, Seq5)
MAX_ERROR_CASE("SinEstimate.getMaxError", AllSinErr, Seq5)

// CosEstimate.h.
ESTIMATE_CASE("CosEstimate.estimate", SymmetricRaw(io, mode, GTE_C_HALF_PI), AllCos, Seq5)
ESTIMATE_CASE("CosEstimate.estimateRR", PeriodicRaw(io, mode, GTE_C_TWO_PI, GTE_C_HALF_PI), AllCosRR, Seq5)
MAX_ERROR_CASE("CosEstimate.getMaxError", AllCosErr, Seq5)

// TanEstimate.h: x in [-pi/4, pi/4]; TanEstimateRR: std::remainder(x, pi)
// and the |y| <= pi/4 split.
ESTIMATE_CASE("TanEstimate.estimate", SymmetricRaw(io, mode, GTE_C_QUARTER_PI), AllTan, Seq6)
ESTIMATE_CASE("TanEstimate.estimateRR", PeriodicRaw(io, mode, GTE_C_PI, GTE_C_QUARTER_PI), AllTanRR, Seq6)
MAX_ERROR_CASE("TanEstimate.getMaxError", AllTanErr, Seq6)

// ATanEstimate.h: x in [-1, 1]; ATanEstimateRR: any x.
ESTIMATE_CASE("ATanEstimate.estimate", SymmetricRaw(io, mode, 1.0), AllATan, Seq6)
ESTIMATE_CASE("ATanEstimate.estimateRR", ATanRRRaw(io, mode), AllATanRR, Seq6)
MAX_ERROR_CASE("ATanEstimate.getMaxError", AllATanErr, Seq6)

// ACosEstimate.h and ASinEstimate.h: x in [0, 1].
ESTIMATE_CASE("ACosEstimate.estimate", IntervalRaw(io, mode, 0.0, 1.0), AllACos, Seq8)
MAX_ERROR_CASE("ACosEstimate.getMaxError", AllACosErr, Seq8)
ESTIMATE_CASE("ASinEstimate.estimate", IntervalRaw(io, mode, 0.0, 1.0), AllASin, Seq8)
MAX_ERROR_CASE("ASinEstimate.getMaxError", AllASinErr, Seq8)

// SqrtEstimate.h and InvSqrtEstimate.h: x in [1, 2]; the RR variants take
// x >= 0 (x > 0 for InvSqrt) through frexp/ldexp.
ESTIMATE_CASE("SqrtEstimate.estimate", IntervalRaw(io, mode, 1.0, 2.0), AllSqrt, Seq8)
ESTIMATE_CASE("SqrtEstimate.estimateRR", PositiveRaw(io, mode), AllSqrtRR, Seq8)
MAX_ERROR_CASE("SqrtEstimate.getMaxError", AllSqrtErr, Seq8)
ESTIMATE_CASE("InvSqrtEstimate.estimate", IntervalRaw(io, mode, 1.0, 2.0), AllInvSqrt, Seq8)
ESTIMATE_CASE("InvSqrtEstimate.estimateRR", PositiveRaw(io, mode), AllInvSqrtRR, Seq8)
MAX_ERROR_CASE("InvSqrtEstimate.getMaxError", AllInvSqrtErr, Seq8)

// Exp2EstimateRR and ExpEstimateRR for floor(x*scale) >= 2^31. Upstream
// converts floor(x) with static_cast<int32_t>, which is undefined behaviour
// out of range; MSVC's cvttsd2si yields INT_MIN, so std::ldexp(poly, INT_MIN)
// returns +0 where 2^x overflows to +inf. The port keeps the double exponent
// and returns +inf. x*scale stays finite (k <= 1021), so y = x - floor(x) is
// a number and the result is decided by the exponent alone.
ORACLE_CASE("Exp2Estimate.estimateRR.hugeArgument")
{
    double x = io.given(std::ldexp(io.raw(1.0, 2.0), io.rawInteger(31, 1021)));
    AllExp2RR(io, x, Seq7{});
    AllExpRR(io, x, Seq7{});
}

// ---------------------------------------------------------------------------
// ChebyshevRatio.h (std::sin, std::acos) and ChebyshevRatioEstimate.h
// (arithmetic only).
// ---------------------------------------------------------------------------
namespace
{
    // An interpolation parameter, by mode (mode taken modulo 4):
    //   0    uniform in [0, 1]
    //   1    dyadic k/16, k in [0, 16]
    //   2    0, -0, 1, 1/2 and 1/2 +- 1 ulp (the SlerpRPH-style 2t <= 1
    //        split of the midpoint overloads)
    //   3    uniform in [-1/4, 5/4] (extrapolation)
    double TRaw(oracle::Ctx& io, int32_t mode)
    {
        switch (mode % 4)
        {
        case 1:
            return static_cast<double>(io.rawInteger(0, 16)) / 16.0;
        case 2:
            switch (io.rawInteger(0, 4))
            {
            case 0:  return 0.0;
            case 1:  return -0.0;
            case 2:  return 1.0;
            case 3:  return 0.5;
            default: return Ulps(0.5, io.rawInteger(0, 1) == 0 ? -1 : 1);
            }
        case 3:
            return io.raw(-0.25, 1.25);
        default:
            return io.raw(0.0, 1.0);
        }
    }

    // An angle in [0, pi), by mode:
    //   0, 5, 7  uniform in [0, pi)
    //   1        +-0 (the angle == 0 branch; -0 == 0)
    //   2        tiny: u*2^-k, k in [1, 1074]
    //   3        approaching pi from below: pi - pi*u*2^-k, k in [1, 52],
    //            or pi - 1..3 ulps
    //   4        k*pi/16, k in [0, 15]
    //   6        pi/2 moved by up to 2 ulps
    double AngleRaw(oracle::Ctx& io, int32_t mode)
    {
        switch (mode)
        {
        case 1:
            return io.rawInteger(0, 1) == 0 ? 0.0 : -0.0;
        case 2:
            return std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(1, 1074));
        case 3:
            if (io.rawInteger(0, 1) == 0)
            {
                return GTE_C_PI - GTE_C_PI * std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(1, 52));
            }
            return Ulps(GTE_C_PI, -io.rawInteger(1, 3));
        case 4:
            return GTE_C_PI * (static_cast<double>(io.rawInteger(0, 15)) / 16.0);
        case 6:
            return Ulps(GTE_C_HALF_PI, io.rawInteger(-2, 2));
        default:
            return io.raw(0.0, GTE_C_PI);
        }
    }

    // A cosine for the UsingCosAngle variants, by mode:
    //   0, 7 uniform in [-1 + 2^-20, 1): the angle stays at least 1.4e-3
    //        below pi, where a 1-ulp difference in acos moves the ratio by
    //        at most 3.2e-13 relative (the near-pi regime has its own case)
    //   1    exactly 1 (the angle-0 branch)
    //   2    above 1: 1 + u, 1 + 1..2 ulps, +inf (angle-0 branch)
    //   3    approaching 1: 1 - u*2^-k, k in [1, 53]
    //   4    NaN (neither comparison holds: the angle-0 branch returns t)
    //   5    dyadic k/16, k in [-15, 15]
    //   6    +-0
    double CosRaw(oracle::Ctx& io, int32_t mode)
    {
        double const lo = -1.0 + std::ldexp(1.0, -20);
        switch (mode)
        {
        case 1:
            return 1.0;
        case 2:
            switch (io.rawInteger(0, 2))
            {
            case 0:  return 1.0 + io.raw(0.0, 1.0);
            case 1:  return Ulps(1.0, io.rawInteger(1, 2));
            default: return INF;
            }
        case 3:
            return 1.0 - std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(1, 53));
        case 4:
            return QNAN;
        case 5:
            return static_cast<double>(io.rawInteger(-15, 15)) / 16.0;
        case 6:
            return io.rawInteger(0, 1) == 0 ? 0.0 : -0.0;
        default:
            return io.raw(lo, 1.0);
        }
    }
}

// std::sin: tolerance 1e-12 in the replay, except the angle == 0 branch
// (returns t), which the replay holds to bit identity.
ORACLE_CASE("ChebyshevRatio.ratio")
{
    for (int32_t mode = 0; mode < NUM_MODES; ++mode)
    {
        double t = io.given(TRaw(io, mode + io.index()));
        double angle = io.given(AngleRaw(io, mode));
        io.outReal(ChebyshevRatio(t, angle));
    }
}

ORACLE_CASE("ChebyshevRatio.ratios")
{
    for (int32_t mode = 0; mode < NUM_MODES; ++mode)
    {
        double t = io.given(TRaw(io, mode + io.index()));
        double angle = io.given(AngleRaw(io, mode));
        std::array<double, 2> f = ChebyshevRatios(t, angle);
        io.outReal(f[0]);
        io.outReal(f[1]);
    }
}

// std::acos, std::sin: tolerance 1e-12, except the branch !(cosAngle < 1).
ORACLE_CASE("ChebyshevRatio.ratioUsingCosAngle")
{
    for (int32_t mode = 0; mode < NUM_MODES; ++mode)
    {
        double t = io.given(TRaw(io, mode + io.index()));
        double cosAngle = io.given(CosRaw(io, mode));
        io.outReal(ChebyshevRatioUsingCosAngle(t, cosAngle));
    }
}

ORACLE_CASE("ChebyshevRatio.ratiosUsingCosAngle")
{
    for (int32_t mode = 0; mode < NUM_MODES; ++mode)
    {
        double t = io.given(TRaw(io, mode + io.index()));
        double cosAngle = io.given(CosRaw(io, mode));
        std::array<double, 2> f = ChebyshevRatiosUsingCosAngle(t, cosAngle);
        io.outReal(f[0]);
        io.outReal(f[1]);
    }
}

// cosAngle in (-1, -1 + 2^-20]: the angle is within 1.4e-3 of pi and the
// ratio is ill-conditioned in acos: a 1-ulp change dA of the angle moves
// the ratio by dA/(pi - A) relative. The replay compares with a per-output
// tolerance of 4 ulps(pi)/(pi - A) derived from the recorded cosine.
ORACLE_CASE("ChebyshevRatio.usingCosAngle.nearPi")
{
    for (int32_t mode = 0; mode < NUM_MODES; ++mode)
    {
        double t = io.given(TRaw(io, mode + io.index()));
        double cosAngle = io.given(-1.0 + std::ldexp(io.raw(0.5, 1.0), -io.rawInteger(20, 53)));
        io.outReal(ChebyshevRatioUsingCosAngle(t, cosAngle));
        std::array<double, 2> f = ChebyshevRatiosUsingCosAngle(t, cosAngle);
        io.outReal(f[0]);
        io.outReal(f[1]);
    }
}

// LogError("Invalid angle."): an angle outside [0, pi) (negative, pi itself,
// above pi, +inf, NaN) or a cosine at or below -1. One call per record,
// selected by a recorded integer; about one record in six is valid so that
// the non-throwing side of each test is present too.
ORACLE_CASE("ChebyshevRatio.throwParity")
{
    int32_t which = io.integer(0, 3);
    double t = io.real(0.0, 1.0);
    double value = 0.0;
    int32_t kind = io.rawInteger(0, 5);
    if (which < 2)
    {
        switch (kind)
        {
        case 0:  value = -io.raw(0.0, 4.0); break;
        case 1:  value = -MIN_SUB; break;
        case 2:  value = GTE_C_PI; break;
        case 3:  value = io.rawInteger(0, 1) == 0 ? io.raw(GTE_C_PI, 10.0) : INF; break;
        case 4:  value = QNAN; break;
        default: value = Ulps(GTE_C_PI, -1); break;
        }
    }
    else
    {
        switch (kind)
        {
        case 0:  value = -1.0; break;
        case 1:  value = -1.0 - io.raw(0.0, 4.0); break;
        case 2:  value = Ulps(-1.0, -1); break;
        case 3:  value = -INF; break;
        case 4:  value = -MAX_DBL; break;
        default: value = Ulps(-1.0, 1); break;
        }
    }
    double arg = io.given(value);
    std::array<double, 2> f{};
    switch (which)
    {
    case 0:  f[0] = ChebyshevRatio(t, arg); break;
    case 1:  f = ChebyshevRatios(t, arg); break;
    case 2:  f[0] = ChebyshevRatioUsingCosAngle(t, arg); break;
    default: f = ChebyshevRatiosUsingCosAngle(t, arg); break;
    }
    io.outReal(f[0]);
    io.outReal(f[1]);
}

namespace
{
    template <size_t... I>
    void AllChbEst(oracle::Ctx& io, double t, double x, std::index_sequence<I...>)
    {
        ((io.outReal(ChebyshevRatioEstimate<double, I + 1>(t, x)[0]),
          io.outReal(ChebyshevRatioEstimate<double, I + 1>(t, x)[1])), ...);
    }

    template <size_t... I>
    void AllChbEstR(oracle::Ctx& io, double t, double x, std::index_sequence<I...>)
    {
        ((io.outReal(ChebyshevRatioEstimateR<double, I + 1>(t, x)[0]),
          io.outReal(ChebyshevRatioEstimateR<double, I + 1>(t, x)[1])), ...);
    }
}

// Arithmetic only: exact. x = cos(A) in [0, 1] (A in [0, pi/2]); the modes
// of IntervalRaw include x = 0, -0, 1, values approaching 1 (A -> 0), values
// outside [0, 1] and the wild population. Every degree 1..16 on every sample.
ORACLE_CASE("ChebyshevRatioEstimate.estimate")
{
    for (int32_t mode = 0; mode < NUM_MODES; ++mode)
    {
        double t = io.given(TRaw(io, mode + io.index()));
        double x = io.given(IntervalRaw(io, mode, 0.0, 1.0));
        AllChbEst(io, t, x, std::make_index_sequence<16>{});
    }
}

// Exact. x = cos(A) in [cos(pi/4), 1] (A in [0, pi/4]); every degree 1..12.
ORACLE_CASE("ChebyshevRatioEstimate.estimateR")
{
    for (int32_t mode = 0; mode < NUM_MODES; ++mode)
    {
        double t = io.given(TRaw(io, mode + io.index()));
        double x = io.given(IntervalRaw(io, mode, GTE_C_INV_SQRT_2, 1.0));
        AllChbEstR(io, t, x, std::make_index_sequence<12>{});
    }
}

MAX_ERROR_CASE("ChebyshevRatioEstimate.getMaxError", AllChbErr, std::make_index_sequence<16>)
MAX_ERROR_CASE("ChebyshevRatioEstimate.getMaxErrorR", AllChbErrR, std::make_index_sequence<12>)

// ---------------------------------------------------------------------------
// Slerp.h and SlerpEstimate.h, for N = 2, 3, 4 on every record.
// ---------------------------------------------------------------------------
namespace
{
    template <int32_t N>
    using Arr = std::array<double, N>;

    template <int32_t N>
    double RawDot(Arr<N> const& a, Arr<N> const& b)
    {
        // Upstream's accumulation: seed 0, i ascending.
        double d = 0.0;
        for (int32_t i = 0; i < N; ++i) { d += a[i] * b[i]; }
        return d;
    }

    template <int32_t N>
    Arr<N> RawNormalized(Arr<N> v)
    {
        double len = std::sqrt(RawDot<N>(v, v));
        for (int32_t i = 0; i < N; ++i) { v[i] /= len; }
        return v;
    }

    template <int32_t N>
    Arr<N> RawUnit(oracle::Ctx& io)
    {
        Arr<N> v{};
        for (int32_t attempt = 0; attempt < 64; ++attempt)
        {
            for (int32_t i = 0; i < N; ++i) { v[i] = io.raw(-1.0, 1.0); }
            double len = std::sqrt(RawDot<N>(v, v));
            if (len >= 0.1 && len <= 1.0) { return RawNormalized<N>(v); }
        }
        v.fill(0.0);
        v[0] = 1.0;
        return v;
    }

    // A signed coordinate axis whose other components are +-0 at random.
    template <int32_t N>
    Arr<N> RawAxis(oracle::Ctx& io, int32_t axis)
    {
        Arr<N> v{};
        for (int32_t i = 0; i < N; ++i) { v[i] = RawSign(io) * 0.0; }
        v[axis] = RawSign(io);
        return v;
    }

    // A unit vector with small rational components: (+-1/2)^4 for N = 4,
    // a signed permutation of (3/5, 4/5) or (1/3, 2/3, 2/3) otherwise.
    template <int32_t N>
    Arr<N> RawFraction(oracle::Ctx& io)
    {
        Arr<N> v{};
        if (N == 4)
        {
            for (int32_t i = 0; i < N; ++i) { v[i] = 0.5 * RawSign(io); }
            return v;
        }
        std::array<double, 3> c = (N == 2)
            ? std::array<double, 3>{ 3.0 / 5.0, 4.0 / 5.0, 0.0 }
            : std::array<double, 3>{ 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0 };
        int32_t shift = io.rawInteger(0, N - 1);
        for (int32_t i = 0; i < N; ++i) { v[i] = RawSign(io) * c[(i + shift) % N]; }
        return v;
    }

    template <int32_t N>
    Arr<N> RawPerturbed(oracle::Ctx& io, Arr<N> const& q, double sign, int32_t kmin, int32_t kmax)
    {
        Arr<N> r = RawUnit<N>(io);
        double scale = std::ldexp(1.0, -io.rawInteger(kmin, kmax));
        Arr<N> v{};
        for (int32_t i = 0; i < N; ++i) { v[i] = sign * q[i] + scale * r[i]; }
        return RawNormalized<N>(v);
    }

    // A pair of unit vectors, by mode:
    //   0  independent uniform unit vectors
    //   1  signed coordinate axes with +-0 elsewhere (dot 0 or 1)
    //   2  q1 = q0 (dot within an ulp of 1: either side of the cosA < 1 test)
    //   3  obtuse: angle uniform in (pi/2, pi - 0.1]
    //   4  near zero angle: q1 = normalize(q0 + 2^-k r), k in [10, 40]
    //   5  small rational components (dot -1/2, 0, 1/2, ... for N = 4)
    //   6  near antipodal: q1 = normalize(-q0 + 2^-k r), k in [1, 30]; only
    //      when nearPi is set (the midpoint overloads)
    // Without nearPi a pair whose dot is below -0.995 (angle above
    // pi - 0.1) has q1 negated: Slerp's ratios are ill-conditioned in acos
    // there (a 1-ulp angle change moves the result by ulp/(pi - A)^2), and
    // dot = -1 is the invalid antipodal configuration of throwParity.
    template <int32_t N>
    void RawPair(oracle::Ctx& io, int32_t mode, bool nearPi, Arr<N>& q0, Arr<N>& q1)
    {
        switch (nearPi ? mode % 7 : mode % 6)
        {
        case 1:
            q0 = RawAxis<N>(io, io.rawInteger(0, N - 1));
            q1 = RawAxis<N>(io, io.rawInteger(0, N - 1));
            break;
        case 2:
            q0 = RawUnit<N>(io);
            q1 = q0;
            break;
        case 3:
        {
            q0 = RawUnit<N>(io);
            Arr<N> w = RawUnit<N>(io);
            double d = RawDot<N>(w, q0);
            for (int32_t i = 0; i < N; ++i) { w[i] -= d * q0[i]; }
            w = RawNormalized<N>(w);
            double A = io.raw(GTE_C_HALF_PI, GTE_C_PI - 0.1);
            double c = std::cos(A), s = std::sin(A);
            for (int32_t i = 0; i < N; ++i) { q1[i] = c * q0[i] + s * w[i]; }
            q1 = RawNormalized<N>(q1);
            break;
        }
        case 4:
            q0 = RawUnit<N>(io);
            q1 = RawPerturbed<N>(io, q0, 1.0, 10, 40);
            break;
        case 5:
            q0 = RawFraction<N>(io);
            q1 = RawFraction<N>(io);
            break;
        case 6:
            q0 = RawUnit<N>(io);
            q1 = RawPerturbed<N>(io, q0, -1.0, 1, 30);
            break;
        default:
            q0 = RawUnit<N>(io);
            q1 = RawUnit<N>(io);
            break;
        }
        double const limit = nearPi ? -1.0 : -0.995;
        if (RawDot<N>(q0, q1) <= limit)
        {
            for (int32_t i = 0; i < N; ++i) { q1[i] = -q1[i]; }
        }
    }

    template <int32_t N>
    Arr<N> GivenArr(oracle::Ctx& io, Arr<N> const& v)
    {
        for (int32_t i = 0; i < N; ++i) { io.given(v[i]); }
        return v;
    }

    template <int32_t N>
    void OutArr(oracle::Ctx& io, Arr<N> const& v)
    {
        for (int32_t i = 0; i < N; ++i) { io.outReal(v[i]); }
    }
}

namespace
{
    // The three overloads share their input construction. Every record runs
    // N = 2, 3 and 4; the pair mode cycles with the record index and N.
    // The documented preprocessing (negate q1 when the dot is negative, so
    // the angle is at most pi/2) is applied on odd records of the cosA
    // overloads, whose cosA is then recorded as the upstream dot of the
    // recorded pair.
    // Upstream is called through non-inlined wrappers. Slerp and
    // SlerpEstimate accumulate result.fill(0); result[i] += f0*a[i] + f1*b[i],
    // and IEEE 754 makes 0 + (-0) = +0. Where MSVC 19.44 (/O2 /fp:precise)
    // inlined Slerp<double, 2> into the case body it knew result[i] == 0,
    // dropped the addition (mulpd/addpd with no add of the seed) and so
    // returned -0 in slots where both products are -0, contrary to its own
    // documented /fp:precise treatment of -0.0; the out-of-line
    // instantiation keeps the addsd. Which slots the optimizer rewrote
    // depended on the inlining context (N = 2 dot overload, N = 4 cosA
    // overload in the first build), so without the wrappers the sign of those
    // zeros is not a property of the upstream source. The replay holds every
    // such slot to bit identity, which verifies that the wrappers keep the
    // seed on every record.
    template <int32_t N>
    __declspec(noinline) Arr<N> CallSlerp(double t, Arr<N> const& a, Arr<N> const& b)
    {
        return Slerp<double, N>(t, a, b);
    }

    template <int32_t N>
    __declspec(noinline) Arr<N> CallSlerp(double t, Arr<N> const& a, Arr<N> const& b, double c)
    {
        return Slerp<double, N>(t, a, b, c);
    }

    template <int32_t N>
    __declspec(noinline) Arr<N> CallSlerp(double t, Arr<N> const& a, Arr<N> const& b,
        Arr<N> const& h, double c)
    {
        return Slerp<double, N>(t, a, b, h, c);
    }

    template <int32_t N, size_t D>
    __declspec(noinline) Arr<N> CallSlerpEst(double t, Arr<N> const& a, Arr<N> const& b)
    {
        return SlerpEstimate<double, N, D>(t, a, b);
    }

    template <int32_t N, size_t D>
    __declspec(noinline) Arr<N> CallSlerpEst(double t, Arr<N> const& a, Arr<N> const& b, double c)
    {
        return SlerpEstimate<double, N, D>(t, a, b, c);
    }

    template <int32_t N, size_t D>
    __declspec(noinline) Arr<N> CallSlerpEst(double t, Arr<N> const& a, Arr<N> const& b,
        Arr<N> const& h, double c)
    {
        return SlerpEstimate<double, N, D>(t, a, b, h, c);
    }

    template <int32_t N>
    struct PairInputs
    {
        Arr<N> q0, q1;
        double t;
    };

    template <int32_t N>
    PairInputs<N> DrawPair(oracle::Ctx& io, bool nearPi, bool preprocess)
    {
        Arr<N> a{}, b{};
        RawPair<N>(io, io.index() + N, nearPi, a, b);
        if (preprocess && RawDot<N>(a, b) < 0.0)
        {
            for (int32_t i = 0; i < N; ++i) { b[i] = -b[i]; }
        }
        PairInputs<N> in{};
        in.q0 = GivenArr<N>(io, a);
        in.q1 = GivenArr<N>(io, b);
        in.t = io.given(TRaw(io, io.rawInteger(0, 3)));
        return in;
    }

    // Midpoint preprocessing as documented in Slerp.h: cosAH =
    // sqrt((1 + cosA)/2), qh = (q0 + q1)/(2*cosAH). Records qh, then cosAH.
    template <int32_t N>
    void DrawMidpoint(oracle::Ctx& io, PairInputs<N> const& in, Arr<N>& qh, double& cosAH)
    {
        double cosA = RawDot<N>(in.q0, in.q1);
        double c = std::sqrt((1.0 + cosA) / 2.0);
        Arr<N> h{};
        for (int32_t i = 0; i < N; ++i) { h[i] = (in.q0[i] + in.q1[i]) / (2.0 * c); }
        qh = GivenArr<N>(io, h);
        cosAH = io.given(c);
    }

    template <int32_t N>
    void SlerpBody(oracle::Ctx& io)
    {
        PairInputs<N> in = DrawPair<N>(io, false, false);
        OutArr<N>(io, CallSlerp<N>(in.t, in.q0, in.q1));
    }

    template <int32_t N>
    void SlerpCosBody(oracle::Ctx& io)
    {
        PairInputs<N> in = DrawPair<N>(io, false, io.index() % 2 == 1);
        double cosA = io.given(RawDot<N>(in.q0, in.q1));
        OutArr<N>(io, CallSlerp<N>(in.t, in.q0, in.q1, cosA));
    }

    template <int32_t N>
    void SlerpMidBody(oracle::Ctx& io)
    {
        PairInputs<N> in = DrawPair<N>(io, true, false);
        Arr<N> qh{};
        double cosAH = 0.0;
        DrawMidpoint<N>(io, in, qh, cosAH);
        OutArr<N>(io, CallSlerp<N>(in.t, in.q0, in.q1, qh, cosAH));
    }

    template <int32_t N, size_t... I>
    void AllSlerpEst(oracle::Ctx& io, PairInputs<N> const& in, std::index_sequence<I...>)
    {
        (OutArr<N>(io, CallSlerpEst<N, I + 1>(in.t, in.q0, in.q1)), ...);
    }

    template <int32_t N, size_t... I>
    void AllSlerpEstCos(oracle::Ctx& io, PairInputs<N> const& in, double cosA,
        std::index_sequence<I...>)
    {
        (OutArr<N>(io, CallSlerpEst<N, I + 1>(in.t, in.q0, in.q1, cosA)), ...);
    }

    template <int32_t N, size_t... I>
    void AllSlerpEstMid(oracle::Ctx& io, PairInputs<N> const& in, Arr<N> const& qh,
        double cosAH, std::index_sequence<I...>)
    {
        (OutArr<N>(io, CallSlerpEst<N, I + 1>(in.t, in.q0, in.q1, qh, cosAH)), ...);
    }

    template <int32_t N>
    void SlerpEstBody(oracle::Ctx& io)
    {
        PairInputs<N> in = DrawPair<N>(io, false, false);
        AllSlerpEst<N>(io, in, std::make_index_sequence<16>{});
    }

    template <int32_t N>
    void SlerpEstCosBody(oracle::Ctx& io)
    {
        PairInputs<N> in = DrawPair<N>(io, false, io.index() % 2 == 1);
        double cosA = io.given(RawDot<N>(in.q0, in.q1));
        AllSlerpEstCos<N>(io, in, cosA, std::make_index_sequence<16>{});
    }

    template <int32_t N>
    void SlerpEstMidBody(oracle::Ctx& io)
    {
        PairInputs<N> in = DrawPair<N>(io, true, false);
        Arr<N> qh{};
        double cosAH = 0.0;
        DrawMidpoint<N>(io, in, qh, cosAH);
        AllSlerpEstMid<N>(io, in, qh, cosAH, std::make_index_sequence<16>{});
    }
}

// std::acos, std::sin through ChebyshevRatiosUsingCosAngle: tolerance 1e-12,
// with the angle-0 branch (dot >= 1) and the +-0 components (both inputs
// zero in that slot) held to bit identity by the replay.
ORACLE_CASE("Slerp.slerp")
{
    SlerpBody<2>(io);
    SlerpBody<3>(io);
    SlerpBody<4>(io);
}

ORACLE_CASE("Slerp.slerpCosAngle")
{
    SlerpCosBody<2>(io);
    SlerpCosBody<3>(io);
    SlerpCosBody<4>(io);
}

// Angles up to (not including) pi: the midpoint's half angle stays below
// pi/2, which is the point of the overload.
ORACLE_CASE("Slerp.slerpMidpoint")
{
    SlerpMidBody<2>(io);
    SlerpMidBody<3>(io);
    SlerpMidBody<4>(io);
}

// LogError("Invalid angle.") from ChebyshevRatiosUsingCosAngle: antipodal
// pairs (dot exactly -1), cosA <= -1 and cosAH <= -1. N = 4. One record in
// four per overload is valid.
ORACLE_CASE("Slerp.throwParity")
{
    int32_t which = io.integer(0, 2);
    int32_t valid = io.integer(0, 3);
    Arr<4> a = RawAxis<4>(io, io.rawInteger(0, 3));
    Arr<4> b = a;
    if (which == 0 && valid != 0)
    {
        for (int32_t i = 0; i < 4; ++i) { b[i] = -b[i]; }
    }
    Arr<4> q0 = GivenArr<4>(io, a);
    Arr<4> q1 = GivenArr<4>(io, b);
    double t = io.real(0.0, 1.0);
    double c = 0.5;
    if (valid != 0)
    {
        switch (io.rawInteger(0, 2))
        {
        case 0:  c = -1.0; break;
        case 1:  c = -1.0 - io.raw(0.0, 4.0); break;
        default: c = -INF; break;
        }
    }
    double cosArg = io.given(c);
    Arr<4> qh = GivenArr<4>(io, a);
    Arr<4> r{};
    switch (which)
    {
    case 0:  r = CallSlerp<4>(t, q0, q1); break;
    case 1:  r = CallSlerp<4>(t, q0, q1, cosArg); break;
    default: r = CallSlerp<4>(t, q0, q1, qh, cosArg); break;
    }
    OutArr<4>(io, r);
}

// Arithmetic only (ChebyshevRatioEstimate): exact, every degree 1..16 for
// every N. The dot overload also receives obtuse pairs (mode 3), outside its
// documented [0, pi/2] range, which it still evaluates deterministically.
ORACLE_CASE("SlerpEstimate.slerpEstimate")
{
    SlerpEstBody<2>(io);
    SlerpEstBody<3>(io);
    SlerpEstBody<4>(io);
}

ORACLE_CASE("SlerpEstimate.slerpEstimateCosAngle")
{
    SlerpEstCosBody<2>(io);
    SlerpEstCosBody<3>(io);
    SlerpEstCosBody<4>(io);
}

ORACLE_CASE("SlerpEstimate.slerpEstimateMidpoint")
{
    SlerpEstMidBody<2>(io);
    SlerpEstMidBody<3>(io);
    SlerpEstMidBody<4>(io);
}
