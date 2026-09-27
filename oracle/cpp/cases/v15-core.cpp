// Group 15 (core): Array2, Array3, Array4, Constants, DCPQuery, FIQuery,
// HashCombine, IEEEBinary, LexicoArray2, Logger, MinHeap, TIQuery,
// TypeTraits, BitHacks, CurveExtractor, Functions, MeshStaticManifold2,
// MeshStaticManifold3, UniqueVerticesTriangles.
//
// Everything here is integer or bit manipulation, comparisons, IEEE-exact
// arithmetic (including std::fma, which is correctly rounded), or container
// bookkeeping, and is compared bit for bit. The one exception is
// Functions.libm (std::atan, std::atan2, std::cos, std::sin, std::exp).
//
// Bit patterns (IEEEBinary encodings, 64-bit BitHacks inputs, hashes) are
// recorded and emitted as two 32-bit halves, high half first, so that every
// integer stays below 2^53 and no NaN payload travels through a double.
//
// Where the sign of zero can reach an output the upstream function is called
// through a __declspec(noinline) wrapper (see ORACLE.md, v23).
//
// Functions.h is compiled with optimization off. MSVC 19.44 /O2 compiles
// clamp's inner conditional (x >= xmax ? xmax : x) to 'minsd xmax, x', which
// returns x instead of xmax when x and xmax are zeros of opposite sign, so
// clamp(+0, -2, -0) is +0 where the C++ source says -0 (and
// clamp(-0, NaN, +0) is -0 where it says +0). The noinline wrapper does not
// help, because the out-of-line copy is compiled the same way; /fp:strict
// does not either; /Od and this pragma give the source's value. Measured on
// the deep run: without the pragma 5 of 2000 Functions.clamp records differ
// from the source semantics, all of them signed-zero ties, and no other
// case of the family changes. The standard headers Functions.h includes
// come first so that the pragma reaches the GTE code alone; no header
// included by Oracle.h includes Functions.h.
#define ORACLE_FAMILY "v15-core"
#include <cmath>
#include <cstdint>
#pragma optimize("", off)
#include <Mathematics/Functions.h>
#pragma optimize("", on)
#include "Oracle.h"

// MinHeap keeps its pointer array private; the cases compare the whole
// internal state, so the header is included with 'private' opened up. Every
// standard header it uses is included first, so the macro reaches
// MinHeap.h alone.
#include <cstdint>
#include <vector>
#define private public
#include <Mathematics/MinHeap.h>
#undef private

#include <Mathematics/Array2.h>
#include <Mathematics/Array3.h>
#include <Mathematics/Array4.h>
#include <Mathematics/BitHacks.h>
#include <Mathematics/BSNumber.h>
#include <Mathematics/BSRational.h>
#include <Mathematics/Constants.h>
#include <Mathematics/CurveExtractor.h>
#include <Mathematics/DCPQuery.h>
#include <Mathematics/FIQuery.h>
#include <Mathematics/Functions.h>
#include <Mathematics/HashCombine.h>
#include <Mathematics/IEEEBinary.h>
#include <Mathematics/LexicoArray2.h>
#include <Mathematics/Logger.h>
#include <Mathematics/MeshStaticManifold2.h>
#include <Mathematics/MeshStaticManifold3.h>
#include <Mathematics/QFNumber.h>
#include <Mathematics/TIQuery.h>
#include <Mathematics/TypeTraits.h>
#include <Mathematics/UIntegerAP32.h>
// The precondition checks of UniqueVerticesTriangles are compiled in (the
// port's 'validate' flag); every main case feeds valid inputs, so they only
// matter for the UniqueVerticesTriangles.validate throw-parity case.
#define GTL_VALIDATE_UNIQUE_VERTICES_TRIANGLES
#include <Mathematics/UniqueVerticesTriangles.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <set>
#include <utility>

using namespace gte;

namespace
{
    // ------------------------------------------------------------ bits

    double FromBits(uint64_t b)
    {
        double d;
        std::memcpy(&d, &b, sizeof(d));
        return d;
    }

    uint64_t ToBits(double d)
    {
        uint64_t b;
        std::memcpy(&b, &d, sizeof(b));
        return b;
    }

    uint32_t RawU32(oracle::Ctx& io)
    {
        uint32_t hi = static_cast<uint32_t>(io.rawInteger(0, 0xFFFF));
        uint32_t lo = static_cast<uint32_t>(io.rawInteger(0, 0xFFFF));
        return (hi << 16) | lo;
    }

    uint64_t RawU64(oracle::Ctx& io)
    {
        uint64_t hi = RawU32(io);
        uint64_t lo = RawU32(io);
        return (hi << 32) | lo;
    }

    uint32_t GivenU32(oracle::Ctx& io, uint32_t u)
    {
        io.given(static_cast<double>(u));
        return u;
    }

    uint64_t GivenU64(oracle::Ctx& io, uint64_t b)
    {
        io.given(static_cast<double>(static_cast<uint32_t>(b >> 32)));
        io.given(static_cast<double>(static_cast<uint32_t>(b & 0xFFFFFFFFull)));
        return b;
    }

    void OutU64(oracle::Ctx& io, uint64_t b)
    {
        io.outInt(static_cast<uint32_t>(b >> 32));
        io.outInt(static_cast<uint32_t>(b & 0xFFFFFFFFull));
    }

    // A 32-bit pattern aimed at the bit-hack branches: 0, 1, powers of two
    // and their neighbours, all-ones, high-bit-set values, sparse and dense
    // random patterns.
    uint32_t DrawU32(oracle::Ctx& io)
    {
        int mode = io.rawInteger(0, 9);
        int k = io.rawInteger(0, 31);
        uint32_t p = 1u << k;
        switch (mode)
        {
        case 0: return static_cast<uint32_t>(io.rawInteger(0, 2));
        case 1: return p;
        case 2: return p - 1u;
        case 3: return p + 1u;
        case 4: return 0xFFFFFFFFu - static_cast<uint32_t>(io.rawInteger(0, 2));
        case 5: return 0x80000000u | RawU32(io);
        case 6: return RawU32(io) & RawU32(io) & RawU32(io);
        case 7: return static_cast<uint32_t>(io.rawInteger(0, 1000));
        case 8: return p | (1u << io.rawInteger(0, 31));
        default: return RawU32(io);
        }
    }

    uint64_t DrawU64(oracle::Ctx& io)
    {
        int mode = io.rawInteger(0, 7);
        int k = io.rawInteger(0, 63);
        uint64_t p = 1ull << k;
        switch (mode)
        {
        case 0: return static_cast<uint64_t>(io.rawInteger(0, 2));
        case 1: return p;
        case 2: return p - 1ull;
        case 3: return p + 1ull;
        case 4: return ~0ull - static_cast<uint64_t>(io.rawInteger(0, 2));
        case 5: return static_cast<uint64_t>(DrawU32(io)) << 32;
        case 6: return static_cast<uint64_t>(DrawU32(io));
        default: return RawU64(io);
        }
    }
}

namespace
{
    // ------------------------------------------------------------ doubles

    double Pow2(int k)
    {
        return std::ldexp(1.0, k);
    }

    // A binary64 pattern from every IEEE class: signed zeros, the smallest
    // and largest subnormals and normals, random subnormals and normals,
    // powers of two, infinities, quiet NaNs (default and with payload),
    // signaling NaNs (payload 1, maximal, random), the neighbours of every
    // class boundary and arbitrary patterns.
    uint64_t DrawBits64(oracle::Ctx& io)
    {
        uint64_t const sign = (io.rawInteger(0, 1) != 0 ? 0x8000000000000000ull : 0ull);
        uint64_t const trailingMask = 0x000FFFFFFFFFFFFFull;
        uint64_t const expMask = 0x7FF0000000000000ull;
        uint64_t const quiet = 0x0008000000000000ull;
        int mode = io.rawInteger(0, 13);
        switch (mode)
        {
        case 0: return sign;                                         // +-0
        case 1: return sign | static_cast<uint64_t>(io.rawInteger(1, 3));  // tiny subnormal
        case 2: return sign | (trailingMask - static_cast<uint64_t>(io.rawInteger(0, 2)));
        case 3: return sign | (RawU64(io) & trailingMask);             // subnormal or 0
        case 4: return sign | (0x0010000000000000ull + static_cast<uint64_t>(io.rawInteger(0, 2)));
        case 5: return sign | (0x7FEFFFFFFFFFFFFFull - static_cast<uint64_t>(io.rawInteger(0, 2)));
        case 6:
        {
            uint64_t biased = static_cast<uint64_t>(io.rawInteger(1, 2046));
            return sign | (biased << 52) | (RawU64(io) & trailingMask);
        }
        case 7:
        {
            uint64_t biased = static_cast<uint64_t>(io.rawInteger(1, 2046));
            return sign | (biased << 52);                               // power of two
        }
        case 8: return sign | expMask;                                  // +-inf
        case 9: return sign | expMask | quiet;                          // default qNaN
        case 10: return sign | expMask | quiet | (RawU64(io) & (quiet - 1ull));
        case 11:
        {
            int which = io.rawInteger(0, 2);
            uint64_t payload = (which == 0 ? 1ull : (which == 1 ? quiet - 1ull
                : (RawU64(io) & (quiet - 1ull))));
            if (payload == 0) { payload = 1; }
            return sign | expMask | payload;                            // sNaN
        }
        case 12: return sign | (expMask - static_cast<uint64_t>(io.rawInteger(0, 1)));
        default: return RawU64(io);
        }
    }

    // The binary32 counterpart of DrawBits64.
    uint32_t DrawBits32(oracle::Ctx& io)
    {
        uint32_t const sign = (io.rawInteger(0, 1) != 0 ? 0x80000000u : 0u);
        uint32_t const trailingMask = 0x007FFFFFu;
        uint32_t const expMask = 0x7F800000u;
        uint32_t const quiet = 0x00400000u;
        int mode = io.rawInteger(0, 13);
        switch (mode)
        {
        case 0: return sign;
        case 1: return sign | static_cast<uint32_t>(io.rawInteger(1, 3));
        case 2: return sign | (trailingMask - static_cast<uint32_t>(io.rawInteger(0, 2)));
        case 3: return sign | (RawU32(io) & trailingMask);
        case 4: return sign | (0x00800000u + static_cast<uint32_t>(io.rawInteger(0, 2)));
        case 5: return sign | (0x7F7FFFFFu - static_cast<uint32_t>(io.rawInteger(0, 2)));
        case 6:
        {
            uint32_t biased = static_cast<uint32_t>(io.rawInteger(1, 254));
            return sign | (biased << 23) | (RawU32(io) & trailingMask);
        }
        case 7:
        {
            uint32_t biased = static_cast<uint32_t>(io.rawInteger(1, 254));
            return sign | (biased << 23);
        }
        case 8: return sign | expMask;
        case 9: return sign | expMask | quiet;
        case 10: return sign | expMask | quiet | (RawU32(io) & (quiet - 1u));
        case 11:
        {
            int which = io.rawInteger(0, 2);
            uint32_t payload = (which == 0 ? 1u : (which == 1 ? quiet - 1u
                : (RawU32(io) & (quiet - 1u))));
            if (payload == 0) { payload = 1; }
            return sign | expMask | payload;
        }
        case 12: return sign | (expMask - static_cast<uint32_t>(io.rawInteger(0, 1)));
        default: return RawU32(io);
        }
    }

    // A double for the scalar functions: uniform, small lattice, dyadic
    // (exact ties with clamp bounds and saturate's 0 and 1), signed zeros,
    // infinities, NaN, huge and tiny magnitudes, arbitrary patterns.
    double DrawScalar(oracle::Ctx& io)
    {
        double s = (io.rawInteger(0, 1) != 0 ? -1.0 : 1.0);
        int mode = io.rawInteger(0, 9);
        switch (mode)
        {
        case 0: return io.raw(-2.0, 2.0);
        case 1: return static_cast<double>(io.rawInteger(-3, 3));
        case 2: return io.rawInteger(-8, 16) / 8.0;
        case 3: return s * 0.0;
        case 4: return s * std::numeric_limits<double>::infinity();
        case 5: return std::numeric_limits<double>::quiet_NaN();
        case 6: return s * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(500, 1023));
        case 7: return s * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(-1074, -500));
        case 8: return s * (1.0 - std::ldexp(1.0, -io.rawInteger(1, 53)));
        default: return FromBits(DrawBits64(io));
        }
    }

    // A finite double over the whole binary64 range for fma: normal values
    // over a wide exponent range, dyadic mantissas, subnormals, signed zeros.
    double DrawFinite(oracle::Ctx& io, int lo, int hi)
    {
        double s = (io.rawInteger(0, 1) != 0 ? -1.0 : 1.0);
        int mode = io.rawInteger(0, 5);
        switch (mode)
        {
        case 0: return s * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(lo, hi));
        case 1: return s * std::ldexp(1.0 + io.rawInteger(0, 15) / 16.0, io.rawInteger(lo, hi));
        case 2: return s * std::ldexp(Pow2(io.rawInteger(1, 53)) - 1.0, io.rawInteger(lo, hi) - 52);
        case 3: return io.raw(-4.0, 4.0);
        case 4: return s * FromBits(RawU64(io) & 0x000FFFFFFFFFFFFFull);
        default: return static_cast<double>(io.rawInteger(-9, 9));
        }
    }
}

// ================================================================ Constants

ORACLE_CASE("Constants.values")
{
    io.outReal(GTE_C_PI);
    io.outReal(GTE_C_HALF_PI);
    io.outReal(GTE_C_QUARTER_PI);
    io.outReal(GTE_C_TWO_PI);
    io.outReal(GTE_C_INV_PI);
    io.outReal(GTE_C_INV_TWO_PI);
    io.outReal(GTE_C_INV_HALF_PI);
    io.outReal(GTE_C_DEG_TO_RAD);
    io.outReal(GTE_C_RAD_TO_DEG);
    io.outReal(GTE_C_SQRT_2);
    io.outReal(GTE_C_INV_SQRT_2);
    io.outReal(GTE_C_LN_2);
    io.outReal(GTE_C_INV_LN_2);
    io.outReal(GTE_C_LN_10);
    io.outReal(GTE_C_INV_LN_10);
}

// ================================================================ Functions
//
// Only the double overloads are ported (the float overloads are listed under
// "Not covered"). stdMax / stdMin are port-only additions documented as the
// ports of std::max / std::min; they are compared with MSVC's std::max and
// std::min, which are (a < b ? b : a) and (b < a ? b : a).

namespace
{
    __declspec(noinline) double CallSign(double x) { return gte::sign(x); }
    __declspec(noinline) int32_t CallISign(double x) { return gte::isign(x); }
    __declspec(noinline) double CallSaturate(double x) { return gte::saturate(x); }
    __declspec(noinline) double CallSqr(double x) { return gte::sqr(x); }
    __declspec(noinline) double CallInvSqrt(double x) { return gte::invsqrt(x); }
    __declspec(noinline) double CallClamp(double x, double a, double b) { return gte::clamp(x, a, b); }
    __declspec(noinline) double CallStdMax(double a, double b) { return std::max(a, b); }
    __declspec(noinline) double CallStdMin(double a, double b) { return std::min(a, b); }
    __declspec(noinline) double CallFMA(double u, double v, double w) { return gte::FMA(u, v, w); }
    __declspec(noinline) double CallSOP(double u, double v, double w, double z) { return gte::RobustSOP(u, v, w, z); }
    __declspec(noinline) double CallDOP(double u, double v, double w, double z) { return gte::RobustDOP(u, v, w, z); }
}

ORACLE_CASE("Functions.arithmetic")
{
    // sign, isign, saturate, sqr, invsqrt (1/std::sqrt is IEEE exact).
    double x = io.given(DrawScalar(io));
    io.outReal(CallSign(x));
    io.outInt(CallISign(x));
    io.outReal(CallSaturate(x));
    io.outReal(CallSqr(x));
    io.outReal(CallInvSqrt(x));
}

ORACLE_CASE("Functions.clamp")
{
    // Lattice and dyadic draws make x equal to a bound; NaN, signed zeros
    // and xmin > xmax reach every arm of the nested conditional.
    double x = io.given(DrawScalar(io));
    double xmin = io.given(DrawScalar(io));
    double xmax = io.given(DrawScalar(io));
    io.outReal(CallClamp(x, xmin, xmax));
}

ORACLE_CASE("Functions.clamp.signedZero")
{
    // Targeted at the ties the deep run found. Half of the records are the
    // miscompiled configuration itself: x a zero, xmax the other zero and
    // xmin below or NaN, so the source returns xmax and 'minsd' returned x
    // (see the pragma at the top of the file). The other half draw x, xmin
    // and xmax from {+0, -0, NaN, -1, 1}.
    static double const special[] = { 0.0, -0.0, std::numeric_limits<double>::quiet_NaN(),
        -1.0, 1.0 };
    double rx = special[io.rawInteger(0, 4)];
    double rmin = special[io.rawInteger(0, 4)];
    double rmax = special[io.rawInteger(0, 4)];
    if (io.rawInteger(0, 1) == 0)
    {
        rx = (io.rawInteger(0, 1) == 0 ? 0.0 : -0.0);
        rmax = -rx;
        rmin = (io.rawInteger(0, 1) == 0 ? -1.0 : std::numeric_limits<double>::quiet_NaN());
    }
    double x = io.given(rx);
    double xmin = io.given(rmin);
    double xmax = io.given(rmax);
    io.outReal(CallClamp(x, xmin, xmax));
}

ORACLE_CASE("Functions.stdMaxMin")
{
    static double const special[] = { 0.0, -0.0, 1.0, -1.0, 2.0,
        std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN() };
    int ia = io.rawInteger(0, 9);
    int ib = io.rawInteger(0, 9);
    double a = io.given(ia < 8 ? special[ia] : DrawScalar(io));
    double b = io.given(ib < 8 ? special[ib] : DrawScalar(io));
    io.outReal(CallStdMax(a, b));
    io.outReal(CallStdMin(a, b));
}

ORACLE_CASE("Functions.libm")
{
    // std::atan, std::atan2, std::cos, std::sin, std::exp: MSVC's CRT and V8
    // may each be 1 ulp off, so this is the one tolerance case of the file.
    // Every result is O(1) or is compared relative to its own magnitude
    // (exp10), so the default 1e-12 is not loosened. The special arguments
    // (signed zeros, infinities, NaN, multiples of 1/2 for cospi/sinpi) are
    // in the draw because libm returns exact values there.
    int mode = io.rawInteger(0, 3);
    double x = 0.0, y = 0.0;
    if (mode == 0)
    {
        x = io.given(io.raw(-4.0, 4.0));
        y = io.given(io.raw(-4.0, 4.0));
    }
    else if (mode == 1)
    {
        x = io.given(io.rawInteger(-16, 16) / 4.0);
        y = io.given(io.rawInteger(-16, 16) / 4.0);
    }
    else if (mode == 2)
    {
        x = io.given(DrawScalar(io));
        y = io.given(DrawScalar(io));
    }
    else
    {
        x = io.given(io.raw(-300.0, 300.0));
        y = io.given(io.raw(-1.0e6, 1.0e6));
    }
    io.outReal(gte::atandivpi(x));
    io.outReal(gte::atan2divpi(y, x));
    io.outReal(gte::cospi(x));
    io.outReal(gte::sinpi(x));
    io.outReal(gte::exp10(x));
}

ORACLE_CASE("Functions.fma")
{
    // std::fma is correctly rounded (one rounding of the exact u*v + w). The
    // port is a software fma over BigInt. Modes: independent draws over a
    // wide exponent range (products that overflow, underflow into the
    // subnormals, or are absorbed by w); TwoProduct cancellation w = -u*v
    // (the result is the exact rounding error of the product, often
    // subnormal); near-cancellation w = -u*v * (1 + small); a lattice.
    int mode = io.rawInteger(0, 4);
    double u = 0.0, v = 0.0, w = 0.0;
    if (mode == 0)
    {
        u = io.given(DrawFinite(io, -600, 600));
        v = io.given(DrawFinite(io, -600, 600));
        w = io.given(DrawFinite(io, -1074, 1023));
    }
    else if (mode == 1 || mode == 2)
    {
        double ru = DrawFinite(io, mode == 1 ? -60 : -560, mode == 1 ? 60 : -500);
        double rv = DrawFinite(io, mode == 1 ? -60 : -560, mode == 1 ? 60 : -500);
        double rw = -(ru * rv);
        if (io.rawInteger(0, 1) != 0)
        {
            rw = rw * (1.0 + std::ldexp(static_cast<double>(io.rawInteger(-8, 8)), -50));
        }
        u = io.given(ru);
        v = io.given(rv);
        w = io.given(rw);
    }
    else if (mode == 3)
    {
        u = io.given(DrawFinite(io, 500, 1023));
        v = io.given(DrawFinite(io, 0, 700));
        w = io.given(DrawFinite(io, 900, 1023));
    }
    else
    {
        u = io.given(static_cast<double>(io.rawInteger(-7, 7)));
        v = io.given(static_cast<double>(io.rawInteger(-7, 7)));
        w = io.given(static_cast<double>(io.rawInteger(-49, 49)));
    }
    io.outReal(CallFMA(u, v, w));
}

ORACLE_CASE("Functions.fma.ties")
{
    // Exact halfway cases of the single rounding, and their neighbours:
    // w has ulp 2^e and u*v = (2t+1) * 2^(e-1) exactly, so u*v + w lies
    // exactly halfway between two doubles (ties to even). The modes cover
    // normal w, subnormal w, the overflow threshold (MAX_NORMAL plus half an
    // ulp rounds to infinity because MAX_NORMAL's significand is odd), and a
    // product 2^-20 ulp above or below the tie.
    int mode = io.rawInteger(0, 3);
    double s0 = (io.rawInteger(0, 1) != 0 ? -1.0 : 1.0);
    double s1 = (io.rawInteger(0, 1) != 0 ? -1.0 : 1.0);
    double t = static_cast<double>(2 * io.rawInteger(0, 20) + 1);
    double ru = 0.0, rv = 0.0, rw = 0.0;
    if (mode == 0)
    {
        int e = io.rawInteger(-900, 900);
        double j = static_cast<double>(io.rawInteger(0, 1 << 30));
        rw = s0 * std::ldexp(Pow2(52) + j, e);
        ru = s1 * t;
        rv = Pow2(e - 1);
    }
    else if (mode == 1)
    {
        double j = static_cast<double>(2 * io.rawInteger(0, 1 << 20) + io.rawInteger(0, 1));
        rw = s0 * std::ldexp(j, -1074);
        ru = s1 * t * Pow2(-537);
        rv = Pow2(-538);
    }
    else if (mode == 2)
    {
        rw = s0 * std::numeric_limits<double>::max();
        ru = s1 * (io.rawInteger(0, 1) != 0 ? 1.0 : (1.0 - Pow2(-52)));
        rv = Pow2(970);
    }
    else
    {
        int e = io.rawInteger(-900, 900);
        double j = static_cast<double>(io.rawInteger(0, 1 << 30));
        rw = s0 * std::ldexp(Pow2(52) + j, e);
        ru = s1 * (t * Pow2(20) + (io.rawInteger(0, 1) != 0 ? 1.0 : -1.0));
        rv = Pow2(e - 21);
    }
    double u = io.given(ru);
    double v = io.given(rv);
    double w = io.given(rw);
    io.outReal(CallFMA(u, v, w));
}

ORACLE_CASE("Functions.fma.special")
{
    // Signed zeros, infinities, NaN, MAX_NORMAL and the smallest subnormal
    // in every operand position: 0*inf, inf - inf, a finite product that
    // overflows against an infinite w, exact zero sums and their sign.
    static double const special[] = { 0.0, -0.0, 1.0, -1.0,
        std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::max(),
        -std::numeric_limits<double>::max(), std::numeric_limits<double>::denorm_min(),
        -std::numeric_limits<double>::denorm_min(), 0.5 };
    int iu = io.rawInteger(0, 15);
    int iv = io.rawInteger(0, 15);
    int iw = io.rawInteger(0, 15);
    double u = io.given(iu < 12 ? special[iu] : DrawFinite(io, -600, 600));
    double v = io.given(iv < 12 ? special[iv] : DrawFinite(io, -600, 600));
    double w = io.given(iw < 12 ? special[iw] : DrawFinite(io, -1074, 1023));
    io.outReal(CallFMA(u, v, w));
}

namespace
{
    // Four operands for RobustSOP / RobustDOP: independent draws, a lattice,
    // and u*v = +-w*z up to rounding so that the sum cancels and the
    // compensation term decides the result.
    std::array<double, 4> DrawProducts(oracle::Ctx& io, bool isDOP)
    {
        int mode = io.rawInteger(0, 3);
        std::array<double, 4> a{};
        if (mode == 0)
        {
            for (auto& x : a) { x = DrawFinite(io, -200, 200); }
        }
        else if (mode == 1)
        {
            for (auto& x : a) { x = static_cast<double>(io.rawInteger(-9, 9)); }
        }
        else
        {
            a[2] = DrawFinite(io, -60, 60);
            a[3] = DrawFinite(io, -60, 60);
            a[0] = a[2] * (1.0 + std::ldexp(static_cast<double>(io.rawInteger(-4, 4)), -52));
            a[1] = a[3] * (isDOP ? 1.0 : -1.0);
            if (mode == 3)
            {
                a[0] = io.raw(-4.0, 4.0);
                a[1] = (isDOP ? 1.0 : -1.0) * (a[2] * a[3]) / a[0];
            }
        }
        return a;
    }
}

ORACLE_CASE("Functions.robustSOP")
{
    auto r = DrawProducts(io, false);
    double u = io.given(r[0]);
    double v = io.given(r[1]);
    double w = io.given(r[2]);
    double z = io.given(r[3]);
    io.outReal(CallSOP(u, v, w, z));
}

ORACLE_CASE("Functions.robustDOP")
{
    auto r = DrawProducts(io, true);
    double u = io.given(r[0]);
    double v = io.given(r[1]);
    double w = io.given(r[2]);
    double z = io.given(r[3]);
    io.outReal(CallDOP(u, v, w, z));
}

// ================================================================ BitHacks
//
// GTE_THROW_ON_BITHACKS_ERROR is not defined (upstream's default), so the
// int32_t / int64_t overloads reinterpret negative inputs and zero inputs
// return the documented "invalid" 0. Both are deterministic and compared.
// The port has one 32-bit function per operation (ToUint32 of its argument)
// and one bigint function for the 64-bit overloads; the replay calls it with
// the unsigned value and with the signed reinterpretation.

namespace
{
    void EmitBits32(oracle::Ctx& io, uint32_t u)
    {
        int32_t i = static_cast<int32_t>(u);
        io.outBool(BitHacks::IsPowerOfTwo(u));
        io.outBool(BitHacks::IsPowerOfTwo(i));
        io.outInt(BitHacks::Log2OfPowerOfTwo(u));
        io.outInt(BitHacks::Log2OfPowerOfTwo(i));
        io.outInt(BitHacks::GetLeadingBit(u));
        io.outInt(BitHacks::GetLeadingBit(i));
        io.outInt(BitHacks::GetTrailingBit(u));
        io.outInt(BitHacks::GetTrailingBit(i));
        io.outInt(BitHacks::RoundUpToPowerOfTwo(u));
        io.outInt(BitHacks::RoundDownToPowerOfTwo(u));
    }

    void EmitBits64(oracle::Ctx& io, uint64_t u)
    {
        int64_t i = static_cast<int64_t>(u);
        io.outInt(BitHacks::GetLeadingBit(u));
        io.outInt(BitHacks::GetLeadingBit(i));
        io.outInt(BitHacks::GetTrailingBit(u));
        io.outInt(BitHacks::GetTrailingBit(i));
    }
}

ORACLE_CASE("BitHacks.bits32")
{
    uint32_t u = GivenU32(io, DrawU32(io));
    EmitBits32(io, u);
}

ORACLE_CASE("BitHacks.bits64")
{
    uint64_t u = GivenU64(io, DrawU64(io));
    EmitBits64(io, u);
}

ORACLE_CASE("BitHacks.powerSweep")
{
    // Every power of two 2^k and its neighbours 2^k - 1, 2^k + 1 for four
    // consecutive k; the 32-bit functions see the low 32 bits.
    int part = io.integer(0, 15);
    for (int k = 4 * part; k < 4 * part + 4; ++k)
    {
        uint64_t p = 1ull << k;
        uint64_t values[3] = { p - 1ull, p, p + 1ull };
        for (uint64_t v : values)
        {
            EmitBits64(io, v);
            EmitBits32(io, static_cast<uint32_t>(v & 0xFFFFFFFFull));
        }
    }
}

// ================================================================ IEEEBinary

namespace
{
    template <typename Binary>
    void EmitEncoding(oracle::Ctx& io, typename Binary::UIntType e)
    {
        if constexpr (sizeof(e) == 8) { OutU64(io, e); }
        else { io.outInt(e); }
    }

    // Every accessor, classifier and neighbour function of one encoding.
    template <typename Binary>
    void EmitFields(oracle::Ctx& io, Binary const& x)
    {
        using UInt = typename Binary::UIntType;
        EmitEncoding<Binary>(io, x.encoding);
        io.outInt(static_cast<uint32_t>(x.GetSign()));
        io.outInt(static_cast<uint32_t>(x.GetBiased()));
        EmitEncoding<Binary>(io, x.GetTrailing());
        UInt sign = 0, biased = 0, trailing = 0;
        x.GetEncoding(sign, biased, trailing);
        io.outInt(static_cast<uint32_t>(sign));
        io.outInt(static_cast<uint32_t>(biased));
        EmitEncoding<Binary>(io, trailing);
        io.outInt(static_cast<int32_t>(x.GetClassification()));
        io.outBool(x.IsZero());
        io.outBool(x.IsSignMinus());
        io.outBool(x.IsSubnormal());
        io.outBool(x.IsNormal());
        io.outBool(x.IsFinite());
        io.outBool(x.IsInfinite());
        io.outBool(x.IsNaN());
        io.outBool(x.IsQuietNaN());
        io.outBool(x.IsSignalingNaN());
        EmitEncoding<Binary>(io, x.GetNextUp());
        EmitEncoding<Binary>(io, x.GetNextDown());
        Binary y(sign, biased, trailing);
        EmitEncoding<Binary>(io, y.encoding);
    }
}

ORACLE_CASE("IEEEBinary64.fields")
{
    uint64_t bits = GivenU64(io, DrawBits64(io));
    IEEEBinary64 x(bits);
    EmitFields(io, x);
    // The Float member of the union. A NaN is emitted as its encoding above
    // (the harness matches any NaN with any NaN), so it is skipped here.
    if (!x.IsNaN()) { io.outReal(x.number); }
}

ORACLE_CASE("IEEEBinary32.fields")
{
    uint32_t bits = GivenU32(io, DrawBits32(io));
    IEEEBinary32 x(bits);
    EmitFields(io, x);
    if (!x.IsNaN()) { io.outReal(static_cast<double>(x.number)); }
}

ORACLE_CASE("IEEEBinary64.fromNumber")
{
    // IEEEBinary64(double): the union's Float member. The double travels as
    // a recorded double, NaN payloads included (V8 keeps the bits of a
    // double it does not compute with).
    double d = io.given(FromBits(DrawBits64(io)));
    IEEEBinary64 x(d);
    OutU64(io, x.encoding);
    io.outInt(static_cast<int32_t>(x.GetClassification()));
}

ORACLE_CASE("IEEEBinary32.fromNumber")
{
    // IEEEBinary32(float) after the double-to-float conversion a caller
    // holding a double performs (round to nearest even, overflow to
    // infinity, NaN quieted with its high payload bits). The port's
    // fromNumber rounds through a Float32Array.
    int mode = io.rawInteger(0, 3);
    double d = 0.0;
    if (mode == 0)
    {
        d = io.given(FromBits(DrawBits64(io)));
    }
    else if (mode == 1)
    {
        // Halfway between two binary32 values (ties to even) or 1 ulp of
        // binary64 away from the tie.
        uint32_t f = DrawBits32(io) & 0x7F7FFFFFu;
        float ff;
        std::memcpy(&ff, &f, sizeof(ff));
        double lo = static_cast<double>(ff);
        uint64_t b = ToBits(lo) + (1ull << 28);
        b = b + static_cast<uint64_t>(io.rawInteger(-1, 1));
        d = io.given(FromBits(b));
    }
    else if (mode == 2)
    {
        d = io.given(std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(-160, 130)));
    }
    else
    {
        d = io.given(DrawScalar(io));
    }
    IEEEBinary32 x(static_cast<float>(d));
    io.outInt(x.encoding);
    io.outInt(static_cast<int32_t>(x.GetClassification()));
}

ORACLE_CASE("IEEEBinary64.setEncoding")
{
    // IEEEBinary(sign, biased, trailing) with in-range fields and with
    // out-of-range fields: SetEncoding shifts and ORs in UInt arithmetic,
    // so an oversized field wraps modulo 2^64 or spills into the next field.
    int mode = io.rawInteger(0, 3);
    uint64_t rs = static_cast<uint64_t>(io.rawInteger(0, 1));
    uint64_t rb = static_cast<uint64_t>(io.rawInteger(0, 2047));
    uint64_t rt = RawU64(io) & 0x000FFFFFFFFFFFFFull;
    if (mode == 1) { rs = static_cast<uint64_t>(io.rawInteger(2, 7)); }
    if (mode == 2) { rb = static_cast<uint64_t>(RawU32(io)); }
    if (mode == 3) { rt = RawU64(io); }
    uint64_t sign = static_cast<uint64_t>(GivenU32(io, static_cast<uint32_t>(rs)));
    uint64_t biased = static_cast<uint64_t>(GivenU32(io, static_cast<uint32_t>(rb)));
    uint64_t trailing = GivenU64(io, rt);
    IEEEBinary64 x(sign, biased, trailing);
    OutU64(io, x.encoding);
}

ORACLE_CASE("IEEEBinary32.setEncoding")
{
    int mode = io.rawInteger(0, 3);
    uint32_t rs = static_cast<uint32_t>(io.rawInteger(0, 1));
    uint32_t rb = static_cast<uint32_t>(io.rawInteger(0, 255));
    uint32_t rt = RawU32(io) & 0x007FFFFFu;
    if (mode == 1) { rs = static_cast<uint32_t>(io.rawInteger(2, 7)); }
    if (mode == 2) { rb = RawU32(io); }
    if (mode == 3) { rt = RawU32(io); }
    uint32_t sign = GivenU32(io, rs);
    uint32_t biased = GivenU32(io, rb);
    uint32_t trailing = GivenU32(io, rt);
    IEEEBinary32 x(sign, biased, trailing);
    io.outInt(x.encoding);
}

namespace
{
    template <typename Binary>
    void EmitConstants(oracle::Ctx& io)
    {
        io.outInt(Binary::NUM_ENCODING_BITS);
        io.outInt(Binary::NUM_EXPONENT_BITS);
        io.outInt(Binary::NUM_SIGNIFICAND_BITS);
        io.outInt(Binary::NUM_TRAILING_BITS);
        io.outInt(Binary::EXPONENT_BIAS);
        io.outInt(Binary::MAX_BIASED_EXPONENT);
        io.outInt(Binary::MIN_SUB_EXPONENT);
        io.outInt(Binary::MIN_EXPONENT);
        io.outInt(Binary::SIGN_SHIFT);
        EmitEncoding<Binary>(io, Binary::SIGN_MASK);
        EmitEncoding<Binary>(io, Binary::NOT_SIGN_MASK);
        EmitEncoding<Binary>(io, Binary::TRAILING_MASK);
        EmitEncoding<Binary>(io, Binary::EXPONENT_MASK);
        EmitEncoding<Binary>(io, Binary::NAN_QUIET_MASK);
        EmitEncoding<Binary>(io, Binary::NAN_PAYLOAD_MASK);
        EmitEncoding<Binary>(io, Binary::MAX_TRAILING);
        EmitEncoding<Binary>(io, Binary::SUP_TRAILING);
        EmitEncoding<Binary>(io, Binary::POS_ZERO);
        EmitEncoding<Binary>(io, Binary::NEG_ZERO);
        EmitEncoding<Binary>(io, Binary::MIN_SUBNORMAL);
        EmitEncoding<Binary>(io, Binary::MAX_SUBNORMAL);
        EmitEncoding<Binary>(io, Binary::MIN_NORMAL);
        EmitEncoding<Binary>(io, Binary::MAX_NORMAL);
        EmitEncoding<Binary>(io, Binary::POS_INFINITY);
        EmitEncoding<Binary>(io, Binary::NEG_INFINITY);
    }
}

ORACLE_CASE("IEEEBinary.constants")
{
    EmitConstants<IEEEBinary32>(io);
    EmitConstants<IEEEBinary64>(io);
}

// ================================================================ HashCombine
//
// The port's documented contract (porting-status.json, the header comment of
// src/HashCombine.ts) is "deterministic, not equal to any C++ library": it
// hashes a number by folding its binary64 bits (-0 folded as +0) into 32
// bits and combines in 32-bit arithmetic, where MSVC hashes a double with
// 64-bit FNV-1a (after mapping -0 to +0) and combines in 64-bit size_t.
//   * HashCombine.combine compares what is comparable: upstream's combine
//     step, instantiated with a key type whose std::hash is the port's fold,
//     from a 32-bit seed. One step from a seed below 2^32 agrees with the
//     port modulo 2^32 (the carries out of bit 31 never flow back down);
//     later steps do not, because (seed >> 2) brings the carried bits in.
//   * HashCombine.equalities compares the equality structure of HashValue
//     with MSVC's real std::hash<double>: which argument lists hash alike.
//   * HashCombine.hashValue.msvc is a deviation case: the raw values.

namespace
{
    struct PortFolded
    {
        size_t h;
    };
}

namespace std
{
    template <>
    struct hash<PortFolded>
    {
        size_t operator()(PortFolded const& x) const noexcept { return x.h; }
    };
}

namespace
{
    uint32_t PortFold(double x)
    {
        uint64_t b = (x == 0.0 ? 0ull : ToBits(x));
        return static_cast<uint32_t>(b & 0xFFFFFFFFull) ^ static_cast<uint32_t>(b >> 32);
    }
}

ORACLE_CASE("HashCombine.combine")
{
    uint32_t seed = GivenU32(io, DrawU32(io));
    double x = io.given(DrawScalar(io));
    size_t s = static_cast<size_t>(seed);
    HashCombine(s, PortFolded{ static_cast<size_t>(PortFold(x)) });
    io.outInt(static_cast<uint32_t>(s & 0xFFFFFFFFull));
}

namespace
{
    // A variant of x that hashes alike (the same bits, or the other zero)
    // or differently (a neighbour, another NaN payload, a sign flip).
    double Variant(oracle::Ctx& io, double x)
    {
        int mode = io.rawInteger(0, 4);
        uint64_t b = ToBits(x);
        switch (mode)
        {
        case 0: return x;
        case 1: return (x == 0.0 ? FromBits(b ^ 0x8000000000000000ull) : x);
        case 2: return FromBits(b ^ 1ull);
        case 3: return FromBits(b ^ 0x8000000000000000ull);
        default: return (std::isnan(x) ? FromBits(b ^ 0x0000000000000100ull) : x);
        }
    }
}

ORACLE_CASE("HashCombine.equalities")
{
    // HashValue(a, b) == HashValue(c, d), with (c, d) a variant of (a, b)
    // or of (b, a). Signed zeros hash alike on both sides (MSVC's
    // std::hash<double> maps -0 to +0); distinct bit patterns, NaN payloads
    // included, hash differently.
    int zeroBias = io.rawInteger(0, 2);
    double ra = (zeroBias == 0 ? 0.0 : FromBits(DrawBits64(io)));
    double rb = (zeroBias == 1 ? -0.0 : DrawScalar(io));
    bool swap = (io.rawInteger(0, 3) == 0);
    double rc = Variant(io, swap ? rb : ra);
    double rd = Variant(io, swap ? ra : rb);
    double a = io.given(ra);
    double b = io.given(rb);
    double c = io.given(rc);
    double d = io.given(rd);
    size_t h0 = HashValue(a, b);
    size_t h1 = HashValue(c, d);
    size_t h2 = HashValue(a);
    size_t h3 = HashValue(c);
    io.outBool(h0 == h1);
    io.outBool(h2 == h3);
}

ORACLE_CASE("HashCombine.equalities.signalingNaN")
{
    // Targeted at the deep-run finding: a signaling NaN against its quiet
    // twin, against the signed zeros, and against itself. MSVC hashes the
    // raw bits; the port folds them, so the comparison holds only while the
    // signaling NaN reaches the port unquieted (the replay reads it with
    // io.real(), not through a JS array).
    static uint64_t const snans[] = { 0xFFF7FFFFFFFFFFFFull, 0x7FF7FFFFFFFFFFFFull,
        0x7FF0000000000001ull, 0xFFF0000000000001ull, 0x7FF4000000000000ull };
    uint64_t sBits = snans[io.rawInteger(0, 4)];
    int mode = io.rawInteger(0, 3);
    uint64_t cBits = (mode == 0 ? sBits | 0x0008000000000000ull
        : (mode == 1 ? 0ull : (mode == 2 ? 0x8000000000000000ull : sBits)));
    double a = io.given(FromBits(sBits));
    double c = io.given(FromBits(cBits));
    io.outBool(HashValue(a) == HashValue(c));
    io.outBool(HashValue(a, c) == HashValue(c, a));
}

ORACLE_CASE("HashCombine.hashValue.msvc")
{
    // Deviation: the raw HashValue(a, b, c) of MSVC's 64-bit std::hash<double>
    // against the port's 32-bit stand-in (not meant to be equal).
    double a = io.given(DrawScalar(io));
    double b = io.given(DrawScalar(io));
    double c = io.given(DrawScalar(io));
    size_t h = HashValue(a, b, c);
    OutU64(io, static_cast<uint64_t>(h));
}

// ================================================================ TypeTraits

ORACLE_CASE("TypeTraits.traits")
{
    // Compile-time traits, compared with the port's runtime predicates on a
    // value of the corresponding port type. Order: double, BSNumber,
    // BSRational, QFNumber<double, 1>, QFNumber<BSRational, 1>.
    using BSN = BSNumber<UIntegerAP32>;
    using BSR = BSRational<UIntegerAP32>;
    io.outBool(is_arbitrary_precision<double>::value);
    io.outBool(has_division_operator<double>::value);
    io.outBool(is_arbitrary_precision<BSN>::value);
    io.outBool(has_division_operator<BSN>::value);
    io.outBool(is_arbitrary_precision<BSR>::value);
    io.outBool(has_division_operator<BSR>::value);
    io.outBool(is_arbitrary_precision<QFNumber<double, 1>>::value);
    io.outBool(has_division_operator<QFNumber<double, 1>>::value);
    io.outBool(is_arbitrary_precision<QFNumber<BSR, 1>>::value);
    io.outBool(has_division_operator<QFNumber<BSR, 1>>::value);
}

// ================================================================ Logger

ORACLE_CASE("Logger.logAssert")
{
    // LogAssert throws std::runtime_error when the condition is false; the
    // port's logAssert throws Error. Emitted only when it does not throw.
    bool condition = io.boolean();
    LogAssert(condition, "oracle assertion");
    io.outBool(condition);
}

ORACLE_CASE("Logger.logError")
{
    int unused = io.integer(0, 1);
    if (unused >= 0)
    {
        LogError("oracle error");
    }
    io.outInt(unused);
}

// ================================================================ MinHeap
//
// MinHeap<int32_t, double>, keys = insertion ids. A record is a random
// sequence of Insert (including inserts into a full heap, which return
// nullptr), Remove (including from an empty heap), Update (increase,
// decrease, equal value; of live records, removed records and the nullptr
// of a failed insert) and GetMinimum. Values come from a small set so that
// ties, including -0 against +0, decide the sift directions. After every
// operation the case emits its result, GetNumElements and IsValid; at the
// end it emits the whole internal state: for each slot of mPointers the
// identity of the record (its position in mRecords) and its index member,
// and for the live slots the key and value.

namespace
{
    using Heap = MinHeap<int32_t, double>;

    double DrawHeapValue(oracle::Ctx& io, bool allowNaN)
    {
        int mode = io.rawInteger(0, allowNaN ? 9 : 8);
        switch (mode)
        {
        case 0: case 1: case 2: case 3:
            return static_cast<double>(io.rawInteger(0, 3));
        case 4: return 0.0;
        case 5: return -0.0;
        case 6: return io.raw(-1.0, 4.0);
        case 7: return static_cast<double>(io.rawInteger(-2, 6)) * 0.5;
        case 8: return static_cast<double>(io.rawInteger(0, 1));
        default: return std::numeric_limits<double>::quiet_NaN();
        }
    }

    void RunHeapSequence(oracle::Ctx& io, bool allowNaN)
    {
        int maxElements = io.integer(0, 12);
        Heap heap(maxElements);
        std::vector<Heap::Record*> handles;
        int numOps = io.integer(1, 40);
        for (int op = 0; op < numOps; ++op)
        {
            int opcode = io.integer(0, 6);
            if (opcode <= 2)
            {
                double value = io.given(DrawHeapValue(io, allowNaN));
                int32_t key = static_cast<int32_t>(handles.size());
                Heap::Record* record = heap.Insert(key, value);
                handles.push_back(record);
                io.outInt(record != nullptr ? record->key : -1);
            }
            else if (opcode == 3)
            {
                int32_t key = -1;
                double value = 0.0;
                bool ok = heap.Remove(key, value);
                io.outBool(ok);
                if (ok) { io.outInt(key); io.outReal(value); }
            }
            else if (opcode <= 5)
            {
                if (!handles.empty())
                {
                    int which = io.integer(0, static_cast<int>(handles.size()) - 1);
                    double value = io.given(DrawHeapValue(io, allowNaN));
                    heap.Update(handles[which], value);
                }
            }
            else
            {
                int32_t key = -1;
                double value = 0.0;
                bool ok = heap.GetMinimum(key, value);
                io.outBool(ok);
                if (ok) { io.outInt(key); io.outReal(value); }
            }
            io.outInt(heap.GetNumElements());
            io.outBool(heap.IsValid());
        }
        for (int i = 0; i < maxElements; ++i)
        {
            Heap::Record const* record = heap.mPointers[i];
            io.outInt(static_cast<int32_t>(record - heap.mRecords.data()));
            io.outInt(record->index);
            if (i < heap.GetNumElements())
            {
                io.outInt(record->key);
                io.outReal(record->value);
            }
        }
    }
}

ORACLE_CASE("MinHeap.sequence")
{
    RunHeapSequence(io, false);
}

ORACLE_CASE("MinHeap.sequence.nan")
{
    // The same sequences with NaN values: upstream compares with '<' and
    // '<=' directly, and every comparison with a NaN is false.
    RunHeapSequence(io, true);
}

// ================================================================ Array2/3/4
//
// Caller-owned storage (the constructor that takes T* objects) is aliased:
// reads and writes through a[i1][i0] (a[i2][i1][i0], a[i3][i2][i1][i0])
// land at i0 + b0*(i1 + b1*(i2 + b2*i3)). The owned constructor
// value-initializes its std::vector<double> (zeros); the port's generic
// storage is not zero-filled (a documented port note), which the
// Array2.owned.zeroInit deviation case shows.

namespace
{
    std::vector<double> DrawObjects(oracle::Ctx& io, size_t n)
    {
        std::vector<double> objects(n);
        for (size_t i = 0; i < n; ++i)
        {
            objects[i] = io.real(-10.0, 10.0);
        }
        return objects;
    }
}

ORACLE_CASE("Array2.access")
{
    size_t b0 = static_cast<size_t>(io.integer(1, 5));
    size_t b1 = static_cast<size_t>(io.integer(1, 5));
    std::vector<double> objects = DrawObjects(io, b0 * b1);
    Array2<double> a(b0, b1, objects.data());
    io.outInt(a.GetBound0());
    io.outInt(a.GetBound1());
    int numOps = io.integer(1, 12);
    for (int op = 0; op < numOps; ++op)
    {
        int i0 = io.integer(0, static_cast<int>(b0) - 1);
        int i1 = io.integer(0, static_cast<int>(b1) - 1);
        bool write = io.boolean();
        if (write) { a[i1][i0] = io.real(-10.0, 10.0); }
        else { io.outReal(a[i1][i0]); }
    }
    for (double x : objects) { io.outReal(x); }
}

ORACLE_CASE("Array2.owned")
{
    // Owned storage, every element written before it is read.
    size_t b0 = static_cast<size_t>(io.integer(1, 5));
    size_t b1 = static_cast<size_t>(io.integer(1, 5));
    Array2<double> a(b0, b1);
    for (size_t i1 = 0; i1 < b1; ++i1)
    {
        for (size_t i0 = 0; i0 < b0; ++i0)
        {
            a[static_cast<int32_t>(i1)][i0] = io.real(-10.0, 10.0);
        }
    }
    int numReads = io.integer(1, 8);
    for (int r = 0; r < numReads; ++r)
    {
        int i0 = io.integer(0, static_cast<int>(b0) - 1);
        int i1 = io.integer(0, static_cast<int>(b1) - 1);
        io.outReal(a[i1][i0]);
    }
}

ORACLE_CASE("Array2.owned.zeroInit")
{
    // Deviation: reads of never-written owned elements are 0 upstream.
    size_t b0 = static_cast<size_t>(io.integer(1, 5));
    size_t b1 = static_cast<size_t>(io.integer(1, 5));
    Array2<double> a(b0, b1);
    int i0 = io.integer(0, static_cast<int>(b0) - 1);
    int i1 = io.integer(0, static_cast<int>(b1) - 1);
    io.outReal(a[i1][i0]);
}

ORACLE_CASE("Array3.access")
{
    size_t b0 = static_cast<size_t>(io.integer(1, 4));
    size_t b1 = static_cast<size_t>(io.integer(1, 4));
    size_t b2 = static_cast<size_t>(io.integer(1, 4));
    std::vector<double> objects = DrawObjects(io, b0 * b1 * b2);
    Array3<double> a(b0, b1, b2, objects.data());
    io.outInt(a.GetBound0());
    io.outInt(a.GetBound1());
    io.outInt(a.GetBound2());
    int numOps = io.integer(1, 12);
    for (int op = 0; op < numOps; ++op)
    {
        int i0 = io.integer(0, static_cast<int>(b0) - 1);
        int i1 = io.integer(0, static_cast<int>(b1) - 1);
        int i2 = io.integer(0, static_cast<int>(b2) - 1);
        bool write = io.boolean();
        if (write) { a[i2][i1][i0] = io.real(-10.0, 10.0); }
        else { io.outReal(a[i2][i1][i0]); }
    }
    for (double x : objects) { io.outReal(x); }
}

ORACLE_CASE("Array4.access")
{
    size_t b0 = static_cast<size_t>(io.integer(1, 3));
    size_t b1 = static_cast<size_t>(io.integer(1, 3));
    size_t b2 = static_cast<size_t>(io.integer(1, 3));
    size_t b3 = static_cast<size_t>(io.integer(1, 3));
    std::vector<double> objects = DrawObjects(io, b0 * b1 * b2 * b3);
    Array4<double> a(b0, b1, b2, b3, objects.data());
    io.outInt(a.GetBound0());
    io.outInt(a.GetBound1());
    io.outInt(a.GetBound2());
    io.outInt(a.GetBound3());
    int numOps = io.integer(1, 12);
    for (int op = 0; op < numOps; ++op)
    {
        int i0 = io.integer(0, static_cast<int>(b0) - 1);
        int i1 = io.integer(0, static_cast<int>(b1) - 1);
        int i2 = io.integer(0, static_cast<int>(b2) - 1);
        int i3 = io.integer(0, static_cast<int>(b3) - 1);
        bool write = io.boolean();
        if (write) { a[i3][i2][i1][i0] = io.real(-10.0, 10.0); }
        else { io.outReal(a[i3][i2][i1][i0]); }
    }
    for (double x : objects) { io.outReal(x); }
}

ORACLE_CASE("Array3.owned")
{
    size_t b0 = static_cast<size_t>(io.integer(1, 4));
    size_t b1 = static_cast<size_t>(io.integer(1, 4));
    size_t b2 = static_cast<size_t>(io.integer(1, 4));
    Array3<double> a(b0, b1, b2);
    for (size_t i2 = 0; i2 < b2; ++i2)
    {
        for (size_t i1 = 0; i1 < b1; ++i1)
        {
            for (size_t i0 = 0; i0 < b0; ++i0)
            {
                a[static_cast<int32_t>(i2)][i1][i0] = io.real(-10.0, 10.0);
            }
        }
    }
    int i0 = io.integer(0, static_cast<int>(b0) - 1);
    int i1 = io.integer(0, static_cast<int>(b1) - 1);
    int i2 = io.integer(0, static_cast<int>(b2) - 1);
    io.outReal(a[i2][i1][i0]);
}

ORACLE_CASE("Array4.owned")
{
    size_t b0 = static_cast<size_t>(io.integer(1, 3));
    size_t b1 = static_cast<size_t>(io.integer(1, 3));
    size_t b2 = static_cast<size_t>(io.integer(1, 3));
    size_t b3 = static_cast<size_t>(io.integer(1, 3));
    Array4<double> a(b0, b1, b2, b3);
    for (size_t i3 = 0; i3 < b3; ++i3)
    {
        for (size_t i2 = 0; i2 < b2; ++i2)
        {
            for (size_t i1 = 0; i1 < b1; ++i1)
            {
                for (size_t i0 = 0; i0 < b0; ++i0)
                {
                    a[static_cast<int32_t>(i3)][i2][i1][i0] = io.real(-10.0, 10.0);
                }
            }
        }
    }
    int i0 = io.integer(0, static_cast<int>(b0) - 1);
    int i1 = io.integer(0, static_cast<int>(b1) - 1);
    int i2 = io.integer(0, static_cast<int>(b2) - 1);
    int i3 = io.integer(0, static_cast<int>(b3) - 1);
    io.outReal(a[i3][i2][i1][i0]);
}

// ================================================================ LexicoArray2
//
// Modes 0/1: the run-time-dimension specializations (row-major, column-
// major); modes 2-5: compile-time dimensions 3x4 and 2x5 in both orders.
// The port has one class with a run-time flag and run-time dimensions.

namespace
{
    template <typename Lexico>
    void RunLexico(oracle::Ctx& io, Lexico& a, int numRows, int numCols,
        std::vector<double>& storage)
    {
        io.outInt(a.GetNumRows());
        io.outInt(a.GetNumCols());
        int numOps = io.integer(1, 12);
        for (int op = 0; op < numOps; ++op)
        {
            int r = io.integer(0, numRows - 1);
            int c = io.integer(0, numCols - 1);
            bool write = io.boolean();
            if (write) { a(r, c) = io.real(-10.0, 10.0); }
            else { io.outReal(a(r, c)); }
        }
        for (double x : storage) { io.outReal(x); }
    }
}

ORACLE_CASE("LexicoArray2.access")
{
    int mode = io.integer(0, 5);
    int numRows = 0, numCols = 0;
    if (mode <= 1)
    {
        numRows = io.integer(1, 5);
        numCols = io.integer(1, 5);
    }
    else
    {
        numRows = (mode <= 3 ? 3 : 2);
        numCols = (mode <= 3 ? 4 : 5);
    }
    std::vector<double> storage = DrawObjects(io, static_cast<size_t>(numRows * numCols));
    if (mode == 0) { LexicoArray2<true, double> a(numRows, numCols, storage.data()); RunLexico(io, a, numRows, numCols, storage); }
    else if (mode == 1) { LexicoArray2<false, double> a(numRows, numCols, storage.data()); RunLexico(io, a, numRows, numCols, storage); }
    else if (mode == 2) { LexicoArray2<true, double, 3, 4> a(storage.data()); RunLexico(io, a, numRows, numCols, storage); }
    else if (mode == 3) { LexicoArray2<false, double, 3, 4> a(storage.data()); RunLexico(io, a, numRows, numCols, storage); }
    else if (mode == 4) { LexicoArray2<true, double, 2, 5> a(storage.data()); RunLexico(io, a, numRows, numCols, storage); }
    else { LexicoArray2<false, double, 2, 5> a(storage.data()); RunLexico(io, a, numRows, numCols, storage); }
}

// ================================================================ MeshStaticManifold2
//
// Meshes are drawn raw and recorded as numVertices, numThreads, the triangle
// count and the index triples. Generators (the 'kind' passed in):
//   grid:     an nx-by-ny grid (1..4 cells each way), each cell split along a
//             random diagonal, counterclockwise, each triangle kept with
//             probability 0.55..1 (boundaries, holes, bow-tie vertices),
//             then relabelled by a random permutation, each triangle rotated
//             cyclically, the triangle order shuffled and 0..2 unused
//             vertices appended;
//   closed:   a tetrahedron or octahedron surface (no boundary edges);
//   defective: a grid plus one defect: a duplicated triangle, a reversed
//             triangle, a third triangle on an edge, or a triangle with a
//             repeated vertex.
// The upstream multithreaded adjacency update (numThreads >= 2) is drawn as
// well; the port always runs the single-threaded loop. Invalid indices
// (SIZE_MAX upstream, Number.MAX_SAFE_INTEGER in the port) are emitted as -1.
// Every output of the class is an ordered std::vector (insertion order,
// depth-first order, std::map key order), so nothing is re-sorted.

namespace
{
    enum class MeshKind { grid, closed, defective, any };

    struct Mesh2
    {
        size_t numVertices = 0;
        std::vector<std::array<size_t, 3>> triangles;
    };

    template <typename T>
    void Shuffle(oracle::Ctx& io, std::vector<T>& v)
    {
        for (size_t i = v.size(); i > 1; --i)
        {
            size_t j = static_cast<size_t>(io.rawInteger(0, static_cast<int>(i) - 1));
            std::swap(v[i - 1], v[j]);
        }
    }

    // Relabel, rotate, shuffle and append unused vertices.
    void Scramble2(oracle::Ctx& io, Mesh2& mesh)
    {
        std::vector<size_t> perm(mesh.numVertices);
        for (size_t i = 0; i < perm.size(); ++i) { perm[i] = i; }
        Shuffle(io, perm);
        for (auto& tri : mesh.triangles)
        {
            int r = io.rawInteger(0, 2);
            std::array<size_t, 3> t = tri;
            for (int k = 0; k < 3; ++k) { tri[k] = perm[t[(k + r) % 3]]; }
        }
        Shuffle(io, mesh.triangles);
        mesh.numVertices += static_cast<size_t>(io.rawInteger(0, 2));
    }

    Mesh2 GridMesh2(oracle::Ctx& io)
    {
        Mesh2 mesh;
        int nx = io.rawInteger(1, 4);
        int ny = io.rawInteger(1, 4);
        double keep = io.raw(0.55, 1.0);
        mesh.numVertices = static_cast<size_t>((nx + 1) * (ny + 1));
        auto id = [nx](int x, int y) { return static_cast<size_t>(x + (nx + 1) * y); };
        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                size_t a = id(x, y), b = id(x + 1, y), c = id(x + 1, y + 1), d = id(x, y + 1);
                std::array<std::array<size_t, 3>, 2> pair{};
                if (io.rawInteger(0, 1) == 0) { pair = { { { a, b, c }, { a, c, d } } }; }
                else { pair = { { { a, b, d }, { b, c, d } } }; }
                for (auto const& tri : pair)
                {
                    if (io.raw(0.0, 1.0) < keep) { mesh.triangles.push_back(tri); }
                }
            }
        }
        if (mesh.triangles.empty()) { mesh.triangles.push_back({ id(0, 0), id(1, 0), id(1, 1) }); }
        return mesh;
    }

    Mesh2 ClosedMesh2(oracle::Ctx& io)
    {
        Mesh2 mesh;
        if (io.rawInteger(0, 1) == 0)
        {
            mesh.numVertices = 4;
            mesh.triangles = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
        }
        else
        {
            // Octahedron: 0 = +x, 1 = -x, 2 = +y, 3 = -y, 4 = +z, 5 = -z.
            mesh.numVertices = 6;
            mesh.triangles = { { 0, 2, 4 }, { 2, 1, 4 }, { 1, 3, 4 }, { 3, 0, 4 },
                { 2, 0, 5 }, { 1, 2, 5 }, { 3, 1, 5 }, { 0, 3, 5 } };
        }
        return mesh;
    }

    void AddDefect2(oracle::Ctx& io, Mesh2& mesh)
    {
        auto const tri = mesh.triangles[static_cast<size_t>(
            io.rawInteger(0, static_cast<int>(mesh.triangles.size()) - 1))];
        int defect = io.rawInteger(0, 3);
        if (defect == 0) { mesh.triangles.push_back(tri); }
        else if (defect == 1) { mesh.triangles.push_back({ tri[0], tri[2], tri[1] }); }
        else if (defect == 2)
        {
            size_t extra = mesh.numVertices++;
            mesh.triangles.push_back({ tri[1], tri[0], extra });
            mesh.triangles.push_back({ tri[0], tri[1], extra });
        }
        else { mesh.triangles.push_back({ tri[0], tri[0], tri[1] }); }
    }

    Mesh2 DrawMesh2(oracle::Ctx& io, MeshKind kind)
    {
        if (kind == MeshKind::any)
        {
            int k = io.rawInteger(0, 5);
            kind = (k <= 2 ? MeshKind::grid : (k == 3 ? MeshKind::closed : MeshKind::defective));
        }
        Mesh2 mesh = (kind == MeshKind::closed ? ClosedMesh2(io) : GridMesh2(io));
        if (kind == MeshKind::defective) { AddDefect2(io, mesh); }
        Scramble2(io, mesh);
        return mesh;
    }

    // Record the mesh: numVertices, triangle count, triples.
    void GivenMesh2(oracle::Ctx& io, Mesh2 const& mesh)
    {
        io.given(static_cast<double>(mesh.numVertices));
        io.given(static_cast<double>(mesh.triangles.size()));
        for (auto const& tri : mesh.triangles)
        {
            for (size_t k = 0; k < 3; ++k) { io.given(static_cast<double>(tri[k])); }
        }
    }

    void OutIndex(oracle::Ctx& io, size_t i)
    {
        if (i == std::numeric_limits<size_t>::max()) { io.outInt(-1); }
        else { io.outInt(i); }
    }
}

ORACLE_CASE("MeshStaticManifold2.construct")
{
    Mesh2 mesh = DrawMesh2(io, MeshKind::any);
    GivenMesh2(io, mesh);
    size_t numThreads = static_cast<size_t>(io.integer(0, 4));
    MeshStaticManifold2 m(mesh.numVertices, mesh.triangles, numThreads);
    io.outInt(m.GetMinNumTrianglesAtVertex());
    io.outInt(m.GetMaxNumTrianglesAtVertex());
    for (auto const& tri : m.GetTriangles())
    {
        for (size_t k = 0; k < 3; ++k) { io.outInt(tri[k]); }
    }
    for (auto const& adj : m.GetAdjacents())
    {
        for (size_t k = 0; k < 3; ++k) { OutIndex(io, adj[k]); }
    }
    for (auto const& vertex : m.GetVertices())
    {
        io.outInt(vertex.GetNumAdjacents());
        auto const* tuple = vertex.GetAdjacents();
        for (size_t j = 0; j < vertex.GetNumAdjacents(); ++j, ++tuple)
        {
            for (size_t k = 0; k < 4; ++k) { OutIndex(io, (*tuple)[k]); }
        }
    }
}

namespace
{
    // The first triangle (in input order) that contains the directed edge
    // <v0,v1>, or SIZE_MAX. Independent of the class: a scan of the input.
    size_t DirectedEdgeTriangle(Mesh2 const& mesh, size_t v0, size_t v1)
    {
        for (size_t t = 0; t < mesh.triangles.size(); ++t)
        {
            auto const& tri = mesh.triangles[t];
            for (size_t k = 0; k < 3; ++k)
            {
                if (tri[k] == v0 && tri[(k + 1) % 3] == v1) { return t; }
            }
        }
        return std::numeric_limits<size_t>::max();
    }

    // A query vertex pair: an existing edge in either direction, two random
    // indices in [-1, numVertices] (-1 becomes SIZE_MAX, numVertices is out
    // of range), or v0 == v1.
    std::array<int, 2> DrawEdgeQuery(oracle::Ctx& io, Mesh2 const& mesh)
    {
        int mode = io.rawInteger(0, 5);
        int n = static_cast<int>(mesh.numVertices);
        if (mode <= 2)
        {
            auto const& tri = mesh.triangles[static_cast<size_t>(
                io.rawInteger(0, static_cast<int>(mesh.triangles.size()) - 1))];
            int k = io.rawInteger(0, 2);
            int a = static_cast<int>(tri[k]);
            int b = static_cast<int>(tri[(k + 1) % 3]);
            return (mode == 0 ? std::array<int, 2>{ b, a } : std::array<int, 2>{ a, b });
        }
        if (mode <= 4)
        {
            return { io.rawInteger(-1, n), io.rawInteger(-1, n) };
        }
        int v = io.rawInteger(0, n - 1);
        return { v, v };
    }

    size_t AsIndex(int i)
    {
        return (i < 0 ? std::numeric_limits<size_t>::max() : static_cast<size_t>(i));
    }
}

ORACLE_CASE("MeshStaticManifold2.queries")
{
    // EdgeExists, and GetAdjacentTriangles where upstream is sound: its
    // return value always, adj0 when the directed edge <v0,v1> exists (then
    // upstream's adj0 is that edge's triangle, the documented contract).
    // Upstream's adj1, and its adj0 when only <v1,v0> exists, are the
    // defect #66; see the .adjacentTriangles deviation case.
    Mesh2 mesh = DrawMesh2(io, MeshKind::any);
    GivenMesh2(io, mesh);
    MeshStaticManifold2 m(mesh.numVertices, mesh.triangles, 0);
    int numQueries = io.integer(1, 8);
    for (int q = 0; q < numQueries; ++q)
    {
        auto e = DrawEdgeQuery(io, mesh);
        int i0 = static_cast<int>(io.given(static_cast<double>(e[0])));
        int i1 = static_cast<int>(io.given(static_cast<double>(e[1])));
        size_t v0 = AsIndex(i0), v1 = AsIndex(i1);
        io.outBool(m.EdgeExists(v0, v1));
        size_t adj0 = 0, adj1 = 0;
        bool found = m.GetAdjacentTriangles(v0, v1, adj0, adj1);
        io.outBool(found);
        bool directed = (i0 >= 0 && i1 >= 0 && i0 != i1
            && v0 < mesh.numVertices && v1 < mesh.numVertices
            && DirectedEdgeTriangle(mesh, v0, v1) != std::numeric_limits<size_t>::max());
        io.outBool(directed);
        if (directed) { OutIndex(io, adj0); }
    }
}

ORACLE_CASE("MeshStaticManifold2.adjacentTriangles")
{
    // Deviation #66: upstream returns the adjacency of the edge opposite v0
    // as adj1 (and a valid adj0 when only <v1,v0> exists). Manifold grid
    // meshes; the query edge is redrawn (capped) until upstream's pair
    // differs from the documented contract computed by a scan of the input:
    // adj0 = triangle of <v0,v1>, adj1 = triangle of <v1,v0>.
    Mesh2 mesh = DrawMesh2(io, MeshKind::grid);
    MeshStaticManifold2 probe(mesh.numVertices, mesh.triangles, 0);
    std::array<int, 2> best{ 0, 1 };
    for (int attempt = 0; attempt < 64; ++attempt)
    {
        auto const& tri = mesh.triangles[static_cast<size_t>(
            io.rawInteger(0, static_cast<int>(mesh.triangles.size()) - 1))];
        int k = io.rawInteger(0, 2);
        bool reverse = (io.rawInteger(0, 1) != 0);
        size_t a = tri[k], b = tri[(k + 1) % 3];
        size_t v0 = reverse ? b : a, v1 = reverse ? a : b;
        size_t adj0 = 0, adj1 = 0;
        probe.GetAdjacentTriangles(v0, v1, adj0, adj1);
        best = { static_cast<int>(v0), static_cast<int>(v1) };
        if (adj0 != DirectedEdgeTriangle(mesh, v0, v1) || adj1 != DirectedEdgeTriangle(mesh, v1, v0))
        {
            break;
        }
    }
    GivenMesh2(io, mesh);
    MeshStaticManifold2 m(mesh.numVertices, mesh.triangles, 0);
    size_t v0 = static_cast<size_t>(io.given(static_cast<double>(best[0])));
    size_t v1 = static_cast<size_t>(io.given(static_cast<double>(best[1])));
    size_t adj0 = 0, adj1 = 0;
    io.outBool(m.GetAdjacentTriangles(v0, v1, adj0, adj1));
    OutIndex(io, adj0);
    OutIndex(io, adj1);
}

ORACLE_CASE("MeshStaticManifold2.components")
{
    Mesh2 mesh = DrawMesh2(io, MeshKind::any);
    GivenMesh2(io, mesh);
    MeshStaticManifold2 m(mesh.numVertices, mesh.triangles, 0);
    std::vector<std::vector<size_t>> components;
    m.GetComponents(components);
    io.outInt(components.size());
    for (auto const& component : components)
    {
        io.outInt(component.size());
        for (size_t t : component) { io.outInt(t); }
    }
}

ORACLE_CASE("MeshStaticManifold2.boundaryPolygons")
{
    // Manifold inputs only (grid meshes with holes and bow-tie vertices,
    // closed surfaces): on a nonmanifold mesh the boundary walk can cycle
    // forever or trip its LogAssert. The polygons come in std::map key
    // order of their first boundary edge.
    int k = io.rawInteger(0, 3);
    Mesh2 mesh = DrawMesh2(io, k == 0 ? MeshKind::closed : MeshKind::grid);
    GivenMesh2(io, mesh);
    bool duplicateEndpoints = io.boolean();
    MeshStaticManifold2 m(mesh.numVertices, mesh.triangles, 0);
    std::vector<std::vector<size_t>> polygons;
    m.GetBoundaryPolygons(polygons, duplicateEndpoints);
    io.outInt(polygons.size());
    for (auto const& polygon : polygons)
    {
        io.outInt(polygon.size());
        for (size_t v : polygon) { io.outInt(v); }
    }
}

ORACLE_CASE("MeshStaticManifold2.invalidInput")
{
    // LogAssert(numVertices >= 3 && triangles.size() > 0): throw parity.
    size_t numVertices = static_cast<size_t>(io.integer(0, 4));
    bool empty = io.boolean();
    std::vector<std::array<size_t, 3>> triangles;
    if (!empty) { triangles.push_back({ 0, 1, 2 }); }
    MeshStaticManifold2 m(numVertices, triangles, 0);
    io.outInt(m.GetMaxNumTrianglesAtVertex());
}

// ================================================================ MeshStaticManifold3
//
// Tetrahedral meshes: the Freudenthal (Kuhn) triangulation of an
// nx-by-ny-by-nz block of cubes (1..2 each way, six tetrahedra per cube,
// conforming across cubes), each tetrahedron oriented positively (GTE's
// canonical tetrahedron (0,0,0),(1,0,0),(0,1,0),(0,0,1) has positive
// determinant), kept with probability 0.5..1, relabelled, reordered by a
// random even permutation of its vertices (orientation preserving) and
// shuffled. Defective meshes add a duplicated, a reversed (odd permutation)
// or a degenerate (repeated vertex) tetrahedron.

namespace
{
    struct Mesh3
    {
        size_t numVertices = 0;
        std::vector<std::array<size_t, 4>> tetrahedra;
    };

    Mesh3 DrawMesh3(oracle::Ctx& io, bool allowDefects)
    {
        Mesh3 mesh;
        int n[3] = { io.rawInteger(1, 2), io.rawInteger(1, 2), io.rawInteger(1, 2) };
        double keep = io.raw(0.5, 1.0);
        mesh.numVertices = static_cast<size_t>((n[0] + 1) * (n[1] + 1) * (n[2] + 1));
        auto id = [&n](std::array<int, 3> const& p)
        {
            return static_cast<size_t>(p[0] + (n[0] + 1) * (p[1] + (n[1] + 1) * p[2]));
        };
        static int const perms[6][3] = { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 },
            { 0, 2, 1 }, { 2, 1, 0 }, { 1, 0, 2 } };
        for (int z = 0; z < n[2]; ++z)
        {
            for (int y = 0; y < n[1]; ++y)
            {
                for (int x = 0; x < n[0]; ++x)
                {
                    for (int p = 0; p < 6; ++p)
                    {
                        std::array<int, 3> c{ x, y, z };
                        std::array<size_t, 4> tet{};
                        tet[0] = id(c);
                        for (int s = 0; s < 3; ++s)
                        {
                            ++c[perms[p][s]];
                            tet[s + 1] = id(c);
                        }
                        // Permutations 3..5 are odd: swap to orient positively.
                        if (p >= 3) { std::swap(tet[2], tet[3]); }
                        if (io.raw(0.0, 1.0) < keep) { mesh.tetrahedra.push_back(tet); }
                    }
                }
            }
        }
        if (mesh.tetrahedra.empty()) { mesh.tetrahedra.push_back({ id({ 0, 0, 0 }), id({ 1, 0, 0 }), id({ 1, 1, 0 }), id({ 1, 1, 1 }) }); }
        if (allowDefects && io.rawInteger(0, 2) == 0)
        {
            auto const tet = mesh.tetrahedra[static_cast<size_t>(
                io.rawInteger(0, static_cast<int>(mesh.tetrahedra.size()) - 1))];
            int defect = io.rawInteger(0, 2);
            if (defect == 0) { mesh.tetrahedra.push_back(tet); }
            else if (defect == 1) { mesh.tetrahedra.push_back({ tet[0], tet[1], tet[3], tet[2] }); }
            else { mesh.tetrahedra.push_back({ tet[0], tet[0], tet[1], tet[2] }); }
        }
        // Relabel, reorder by an even permutation, shuffle, append unused.
        std::vector<size_t> perm(mesh.numVertices);
        for (size_t i = 0; i < perm.size(); ++i) { perm[i] = i; }
        Shuffle(io, perm);
        static int const even[12][4] = {
            { 0, 1, 2, 3 }, { 0, 2, 3, 1 }, { 0, 3, 1, 2 }, { 1, 0, 3, 2 },
            { 1, 2, 0, 3 }, { 1, 3, 2, 0 }, { 2, 0, 1, 3 }, { 2, 1, 3, 0 },
            { 2, 3, 0, 1 }, { 3, 0, 2, 1 }, { 3, 1, 0, 2 }, { 3, 2, 1, 0 } };
        for (auto& tet : mesh.tetrahedra)
        {
            int r = io.rawInteger(0, 11);
            std::array<size_t, 4> t = tet;
            for (int k = 0; k < 4; ++k) { tet[k] = perm[t[even[r][k]]]; }
        }
        Shuffle(io, mesh.tetrahedra);
        mesh.numVertices += static_cast<size_t>(io.rawInteger(0, 2));
        return mesh;
    }

    void GivenMesh3(oracle::Ctx& io, Mesh3 const& mesh)
    {
        io.given(static_cast<double>(mesh.numVertices));
        io.given(static_cast<double>(mesh.tetrahedra.size()));
        for (auto const& tet : mesh.tetrahedra)
        {
            for (size_t k = 0; k < 4; ++k) { io.given(static_cast<double>(tet[k])); }
        }
    }

    // A face query: a face of a tetrahedron in one of its six vertex orders,
    // three random indices in [-1, numVertices], or a repeated index.
    std::array<int, 3> DrawFaceQuery(oracle::Ctx& io, Mesh3 const& mesh)
    {
        int mode = io.rawInteger(0, 4);
        int n = static_cast<int>(mesh.numVertices);
        if (mode <= 2)
        {
            auto const& tet = mesh.tetrahedra[static_cast<size_t>(
                io.rawInteger(0, static_cast<int>(mesh.tetrahedra.size()) - 1))];
            auto const& f = MeshStaticManifold3::face[static_cast<size_t>(io.rawInteger(0, 3))];
            std::array<int, 3> q{ static_cast<int>(tet[f[0]]), static_cast<int>(tet[f[1]]),
                static_cast<int>(tet[f[2]]) };
            int r = io.rawInteger(0, 5);
            if (r >= 3) { std::swap(q[1], q[2]); }
            std::rotate(q.begin(), q.begin() + (r % 3), q.end());
            return q;
        }
        if (mode == 3)
        {
            return { io.rawInteger(-1, n), io.rawInteger(-1, n), io.rawInteger(-1, n) };
        }
        int v = io.rawInteger(0, n - 1);
        return { v, v, io.rawInteger(0, n - 1) };
    }
}

ORACLE_CASE("MeshStaticManifold3.construct")
{
    Mesh3 mesh = DrawMesh3(io, true);
    GivenMesh3(io, mesh);
    size_t numThreads = static_cast<size_t>(io.integer(0, 4));
    MeshStaticManifold3 m(mesh.numVertices, mesh.tetrahedra, numThreads);
    io.outInt(m.GetMinNumTetrahedraAtVertex());
    io.outInt(m.GetMaxNumTetrahedraAtVertex());
    for (auto const& tet : m.GetTetrahedra())
    {
        for (size_t k = 0; k < 4; ++k) { io.outInt(tet[k]); }
    }
    for (auto const& adj : m.GetAdjacents())
    {
        for (size_t k = 0; k < 4; ++k) { OutIndex(io, adj[k]); }
    }
    for (auto const& vertex : m.GetVertices())
    {
        io.outInt(vertex.GetNumAdjacents());
        auto const* tuple = vertex.GetAdjacents();
        for (size_t j = 0; j < vertex.GetNumAdjacents(); ++j, ++tuple)
        {
            for (size_t k = 0; k < 5; ++k) { OutIndex(io, (*tuple)[k]); }
        }
    }
}

ORACLE_CASE("MeshStaticManifold3.queries")
{
    // FaceExists and the return value of GetAdjacentTetrahedra. Upstream's
    // adj0/adj1 are the defect #66 (a vertex index returned as a
    // tetrahedron index); see the .adjacentTetrahedra deviation case.
    Mesh3 mesh = DrawMesh3(io, true);
    GivenMesh3(io, mesh);
    MeshStaticManifold3 m(mesh.numVertices, mesh.tetrahedra, 0);
    int numQueries = io.integer(1, 8);
    for (int q = 0; q < numQueries; ++q)
    {
        auto f = DrawFaceQuery(io, mesh);
        int i0 = static_cast<int>(io.given(static_cast<double>(f[0])));
        int i1 = static_cast<int>(io.given(static_cast<double>(f[1])));
        int i2 = static_cast<int>(io.given(static_cast<double>(f[2])));
        size_t v0 = AsIndex(i0), v1 = AsIndex(i1), v2 = AsIndex(i2);
        io.outBool(m.FaceExists(v0, v1, v2));
        size_t adj0 = 0, adj1 = 0;
        io.outBool(m.GetAdjacentTetrahedra(v0, v1, v2, adj0, adj1));
    }
}

ORACLE_CASE("MeshStaticManifold3.adjacentTetrahedra")
{
    // Deviation #66: on existing faces upstream returns element [2] of the
    // 5-tuple (a vertex index) as adj0 and element [3] (the tetrahedron) as
    // adj1. Every record queries a face of the mesh.
    Mesh3 mesh = DrawMesh3(io, false);
    GivenMesh3(io, mesh);
    MeshStaticManifold3 m(mesh.numVertices, mesh.tetrahedra, 0);
    auto const& tet = mesh.tetrahedra[static_cast<size_t>(
        io.integer(0, static_cast<int>(mesh.tetrahedra.size()) - 1))];
    auto const& f = MeshStaticManifold3::face[static_cast<size_t>(io.integer(0, 3))];
    size_t adj0 = 0, adj1 = 0;
    io.outBool(m.GetAdjacentTetrahedra(tet[f[0]], tet[f[1]], tet[f[2]], adj0, adj1));
    OutIndex(io, adj0);
    OutIndex(io, adj1);
}

ORACLE_CASE("MeshStaticManifold3.invalidInput")
{
    // LogAssert(numVertices >= 4 && tetrahedra.size() > 0): throw parity.
    size_t numVertices = static_cast<size_t>(io.integer(0, 5));
    bool empty = io.boolean();
    std::vector<std::array<size_t, 4>> tetrahedra;
    if (!empty) { tetrahedra.push_back({ 0, 1, 2, 3 }); }
    MeshStaticManifold3 m(numVertices, tetrahedra, 0);
    io.outInt(m.GetMaxNumTetrahedraAtVertex());
}

// ================================================================ UniqueVerticesTriangles
//
// VertexType = Vector3<double> (the port's default key joins the components,
// so its equivalence is component-wise ==, like Vector's lexicographic
// operator< as a std::map key), plus a scalar VertexType = double case. The
// vertex lists repeat pool points exactly, as a -0/+0 variant (equivalent
// upstream: the first-inserted representative is kept, sign included) and
// as a 1-ulp neighbour (distinct). Both overloads of each operation (flat
// indices and index triples) are called; outputs are ordered vectors.

namespace
{
    std::vector<Vector3<double>> DrawVertexList(oracle::Ctx& io, size_t count)
    {
        int poolSize = io.rawInteger(1, 6);
        std::vector<Vector3<double>> pool(static_cast<size_t>(poolSize));
        for (auto& p : pool)
        {
            bool lattice = (io.rawInteger(0, 2) != 0);
            for (int k = 0; k < 3; ++k)
            {
                p[k] = (lattice ? static_cast<double>(io.rawInteger(-1, 1)) : io.raw(-1.0, 1.0));
            }
        }
        std::vector<Vector3<double>> list(count);
        for (auto& v : list)
        {
            v = pool[static_cast<size_t>(io.rawInteger(0, poolSize - 1))];
            int variant = io.rawInteger(0, 5);
            int k = io.rawInteger(0, 2);
            if (variant == 4 && v[k] == 0.0) { v[k] = -v[k]; }
            if (variant == 5) { v[k] = FromBits(ToBits(v[k]) + 1ull); }
        }
        return list;
    }

    std::vector<Vector3<double>> GivenVertexList(oracle::Ctx& io, std::vector<Vector3<double>> const& list)
    {
        io.given(static_cast<double>(list.size()));
        for (auto const& v : list) { io.givenVec(v); }
        return list;
    }

    std::vector<int32_t> DrawIndices(oracle::Ctx& io, size_t numVertices)
    {
        int numTriangles = io.integer(1, 5);
        std::vector<int32_t> indices(3 * static_cast<size_t>(numTriangles));
        for (auto& i : indices) { i = io.integer(0, static_cast<int>(numVertices) - 1); }
        return indices;
    }

    std::vector<std::array<int32_t, 3>> ToTriples(std::vector<int32_t> const& indices)
    {
        std::vector<std::array<int32_t, 3>> triples(indices.size() / 3);
        for (size_t t = 0; t < triples.size(); ++t)
        {
            for (size_t j = 0; j < 3; ++j) { triples[t][j] = indices[3 * t + j]; }
        }
        return triples;
    }

    template <typename V>
    void EmitUnique(oracle::Ctx& io, std::vector<V> const& vertices, std::vector<int32_t> const& indices)
    {
        io.outInt(vertices.size());
        for (auto const& v : vertices)
        {
            if constexpr (std::is_same<V, double>::value) { io.outReal(v); }
            else { io.outVec(v); }
        }
        io.outInt(indices.size());
        for (int32_t i : indices) { io.outInt(i); }
    }

    template <typename V>
    void EmitUnique(oracle::Ctx& io, std::vector<V> const& vertices,
        std::vector<std::array<int32_t, 3>> const& triangles)
    {
        std::vector<int32_t> flat;
        for (auto const& t : triangles) { flat.insert(flat.end(), t.begin(), t.end()); }
        EmitUnique(io, vertices, flat);
    }
}

ORACLE_CASE("UniqueVerticesTriangles.generateIndexedTriangles")
{
    size_t numTriangles = static_cast<size_t>(io.rawInteger(1, 5));
    auto inVertices = GivenVertexList(io, DrawVertexList(io, 3 * numTriangles));
    UniqueVerticesTriangles<Vector3<double>> uvt;
    std::vector<Vector3<double>> outVertices;
    std::vector<int32_t> outIndices;
    uvt.GenerateIndexedTriangles(inVertices, outVertices, outIndices);
    EmitUnique(io, outVertices, outIndices);
    std::vector<Vector3<double>> outVertices3;
    std::vector<std::array<int32_t, 3>> outTriangles;
    uvt.GenerateIndexedTriangles(inVertices, outVertices3, outTriangles);
    EmitUnique(io, outVertices3, outTriangles);
}

namespace
{
    enum class UVTOp { duplicate, unused, both };

    void RunUVT(oracle::Ctx& io, UVTOp op)
    {
        size_t numVertices = static_cast<size_t>(io.rawInteger(1, 9));
        auto inVertices = GivenVertexList(io, DrawVertexList(io, numVertices));
        auto inIndices = DrawIndices(io, numVertices);
        auto inTriangles = ToTriples(inIndices);
        UniqueVerticesTriangles<Vector3<double>> uvt;
        std::vector<Vector3<double>> outVertices, outVertices3;
        std::vector<int32_t> outIndices;
        std::vector<std::array<int32_t, 3>> outTriangles;
        if (op == UVTOp::duplicate)
        {
            uvt.RemoveDuplicateVertices(inVertices, inIndices, outVertices, outIndices);
            uvt.RemoveDuplicateVertices(inVertices, inTriangles, outVertices3, outTriangles);
        }
        else if (op == UVTOp::unused)
        {
            uvt.RemoveUnusedVertices(inVertices, inIndices, outVertices, outIndices);
            uvt.RemoveUnusedVertices(inVertices, inTriangles, outVertices3, outTriangles);
        }
        else
        {
            uvt.RemoveDuplicateAndUnusedVertices(inVertices, inIndices, outVertices, outIndices);
            uvt.RemoveDuplicateAndUnusedVertices(inVertices, inTriangles, outVertices3, outTriangles);
        }
        EmitUnique(io, outVertices, outIndices);
        EmitUnique(io, outVertices3, outTriangles);
    }
}

ORACLE_CASE("UniqueVerticesTriangles.removeDuplicateVertices")
{
    RunUVT(io, UVTOp::duplicate);
}

ORACLE_CASE("UniqueVerticesTriangles.removeUnusedVertices")
{
    RunUVT(io, UVTOp::unused);
}

ORACLE_CASE("UniqueVerticesTriangles.removeDuplicateAndUnusedVertices")
{
    RunUVT(io, UVTOp::both);
}

ORACLE_CASE("UniqueVerticesTriangles.validate")
{
    // Throw parity of the GTL_VALIDATE_UNIQUE_VERTICES_TRIANGLES checks:
    // empty inputs, vertex counts that are not a multiple of 3, index lists
    // whose length is not a multiple of 3, indices -1 and numVertices. The
    // operation is drawn; on success only the output sizes are emitted.
    // Half of the records are drawn valid, the other half unconstrained.
    int op = io.integer(0, 5);
    bool valid = io.boolean();
    int nv = 0, ni = 0;
    if (valid)
    {
        nv = (op == 0 || op == 5 ? 3 * io.rawInteger(1, 2) : io.rawInteger(1, 7));
        ni = 3 * io.rawInteger(1, 2);
    }
    else
    {
        nv = io.rawInteger(0, 7);
        ni = io.rawInteger(0, 7);
        if (op >= 3) { ni -= ni % 3; }
    }
    size_t numVertices = static_cast<size_t>(io.given(static_cast<double>(nv)));
    std::vector<Vector3<double>> inVertices(numVertices);
    for (size_t i = 0; i < numVertices; ++i)
    {
        inVertices[i] = { static_cast<double>(i), 0.0, 0.0 };
    }
    size_t numIndices = static_cast<size_t>(io.given(static_cast<double>(ni)));
    std::vector<int32_t> inIndices(numIndices);
    for (auto& i : inIndices)
    {
        int r = (valid ? io.rawInteger(0, nv - 1) : io.rawInteger(-1, nv));
        i = static_cast<int32_t>(io.given(static_cast<double>(r)));
    }
    UniqueVerticesTriangles<Vector3<double>> uvt;
    std::vector<Vector3<double>> outVertices;
    std::vector<int32_t> outIndices;
    std::vector<std::array<int32_t, 3>> outTriangles;
    if (op == 0) { uvt.GenerateIndexedTriangles(inVertices, outVertices, outIndices); }
    else if (op == 1) { uvt.RemoveDuplicateVertices(inVertices, inIndices, outVertices, outIndices); }
    else if (op == 2) { uvt.RemoveUnusedVertices(inVertices, inIndices, outVertices, outIndices); }
    else if (op == 3) { uvt.RemoveDuplicateVertices(inVertices, ToTriples(inIndices), outVertices, outTriangles); }
    else if (op == 4) { uvt.RemoveUnusedVertices(inVertices, ToTriples(inIndices), outVertices, outTriangles); }
    else { uvt.GenerateIndexedTriangles(inVertices, outVertices, outTriangles); }
    io.outInt(outVertices.size());
    io.outInt(outIndices.size() + 3 * outTriangles.size());
}

ORACLE_CASE("UniqueVerticesTriangles.scalar")
{
    // VertexType = double: std::map<double, int32_t>, -0 and +0 equivalent.
    size_t numVertices = static_cast<size_t>(io.integer(1, 9));
    std::vector<double> inVertices(numVertices);
    for (auto& v : inVertices)
    {
        int mode = io.rawInteger(0, 3);
        double x = static_cast<double>(io.rawInteger(-2, 2));
        if (mode == 0 && x == 0.0) { x = -0.0; }
        if (mode == 1) { x = FromBits(ToBits(x) + 1ull); }
        v = io.given(x);
    }
    auto inIndices = DrawIndices(io, numVertices);
    UniqueVerticesTriangles<double> uvt;
    std::vector<double> outVertices;
    std::vector<int32_t> outIndices;
    uvt.RemoveDuplicateAndUnusedVertices(inVertices, inIndices, outVertices, outIndices);
    EmitUnique(io, outVertices, outIndices);
    std::vector<double> genVertices;
    std::vector<int32_t> genIndices;
    std::vector<double> triple(inVertices.begin(), inVertices.begin() + 3 * (numVertices / 3));
    if (!triple.empty())
    {
        uvt.GenerateIndexedTriangles(triple, genVertices, genIndices);
    }
    EmitUnique(io, genVertices, genIndices);
}

// ================================================================ CurveExtractor
//
// The abstract base has public functions of its own: the Vertex and Edge
// structs (constructors, operator==, operator<), MakeUnique, Convert and the
// non-virtual Extract(level, removeDuplicateVertices, vertices, edges). They
// are driven through a case-file subclass whose pure-virtual Extract returns
// a recorded rational vertex/edge list built with the protected AddEdge, as
// the concrete extractors (group 16) do. Mirrored in the replay.

namespace
{
    class ReplayExtractor : public CurveExtractor<int32_t, double>
    {
    public:
        using CurveExtractor<int32_t, double>::Extract;

        ReplayExtractor(int32_t xBound = 2, int32_t yBound = 2)
            :
            CurveExtractor<int32_t, double>(xBound, yBound, msPixels)
        {
        }

        void AddRecordedEdge(std::array<int64_t, 8> const& e)
        {
            AddEdge(mVertices, mEdges, e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]);
        }

        void AddRecordedVertex(std::array<int64_t, 4> const& v)
        {
            AddVertex(mVertices, v[0], v[1], v[2], v[3]);
        }

        virtual void Extract(int32_t, std::vector<Vertex>& vertices,
            std::vector<Edge>& edges) override
        {
            vertices = mVertices;
            edges = mEdges;
        }

    private:
        static int32_t const msPixels[16];
        std::vector<Vertex> mVertices;
        std::vector<Edge> mEdges;
    };

    int32_t const ReplayExtractor::msPixels[16] = {};

    // A rational coordinate numer/denom with numerator and denominator of
    // the same sign (the extractors' invariant): 0..6 over 1..4, in reduced
    // or unreduced form, sometimes both negated.
    std::array<int64_t, 2> DrawRational(oracle::Ctx& io)
    {
        int64_t d = io.rawInteger(1, 4);
        int64_t n = io.rawInteger(0, 6);
        int64_t scale = io.rawInteger(1, 3);
        int64_t sign = (io.rawInteger(0, 3) == 0 ? -1 : 1);
        return { sign * scale * n, sign * scale * d };
    }
}

ORACLE_CASE("CurveExtractor.extract")
{
    // Edges between points of a small pool, each point written in several
    // rational forms (x/y, kx/ky, -x/-y), so vertices repeat under ==, edges
    // repeat in both directions (finding #362: a reversed duplicate survives
    // MakeUnique) and degenerate edges (both endpoints equal) occur.
    ReplayExtractor extractor;
    int poolSize = io.rawInteger(1, 5);
    std::vector<std::array<std::array<int64_t, 2>, 2>> pool(static_cast<size_t>(poolSize));
    for (auto& p : pool) { p = { DrawRational(io), DrawRational(io) }; }
    int numEdges = io.integer(0, 8);
    for (int e = 0; e < numEdges; ++e)
    {
        std::array<int64_t, 8> r{};
        for (int end = 0; end < 2; ++end)
        {
            auto const& p = pool[static_cast<size_t>(io.rawInteger(0, poolSize - 1))];
            for (int c = 0; c < 2; ++c)
            {
                int64_t k = io.rawInteger(1, 3) * (io.rawInteger(0, 3) == 0 ? -1 : 1);
                r[4 * end + 2 * c] = k * p[c][0];
                r[4 * end + 2 * c + 1] = k * p[c][1];
            }
        }
        for (auto& x : r) { x = static_cast<int64_t>(io.given(static_cast<double>(x))); }
        extractor.AddRecordedEdge(r);
    }
    int numSingles = io.integer(0, 2);
    for (int s = 0; s < numSingles; ++s)
    {
        auto const& p = pool[static_cast<size_t>(io.rawInteger(0, poolSize - 1))];
        std::array<int64_t, 4> v{ p[0][0], p[0][1], p[1][0], p[1][1] };
        for (auto& x : v) { x = static_cast<int64_t>(io.given(static_cast<double>(x))); }
        extractor.AddRecordedVertex(v);
    }
    bool removeDuplicates = io.boolean();
    std::vector<std::array<double, 2>> vertices;
    std::vector<ReplayExtractor::Edge> edges;
    extractor.Extract(0, removeDuplicates, vertices, edges);
    io.outInt(vertices.size());
    for (auto const& v : vertices) { io.outReal(v[0]); io.outReal(v[1]); }
    io.outInt(edges.size());
    for (auto const& e : edges) { io.outInt(e.v[0]); io.outInt(e.v[1]); }
}

ORACLE_CASE("CurveExtractor.vertexEdge")
{
    // Vertex construction (sign normalization), == and < in both
    // directions, with components up to 2^31 so that the cross products
    // reach 2^62 (the port switches to bigint above 2^53); Edge construction
    // (index sort), == and <.
    using V = ReplayExtractor::Vertex;
    using E = ReplayExtractor::Edge;
    std::array<int64_t, 8> r{};
    bool big = (io.rawInteger(0, 1) != 0);
    bool related = (io.rawInteger(0, 1) != 0);
    for (int i = 0; i < 4; ++i)
    {
        int64_t d = (big ? static_cast<int64_t>(io.rawInteger(1, 0x7FFFFFFF)) : io.rawInteger(1, 4));
        int64_t n = (big ? static_cast<int64_t>(io.rawInteger(0, 0x7FFFFFFF)) : io.rawInteger(0, 6));
        int64_t s = (io.rawInteger(0, 3) == 0 ? -1 : 1);
        r[2 * i] = s * n;
        r[2 * i + 1] = s * d;
    }
    if (related)
    {
        // The second vertex equals the first as rationals, or differs in y.
        r[4] = r[0]; r[5] = r[1];
        if (io.rawInteger(0, 1) != 0) { r[6] = r[2]; r[7] = r[3]; }
    }
    for (auto& x : r) { x = static_cast<int64_t>(io.given(static_cast<double>(x))); }
    V a(r[0], r[1], r[2], r[3]);
    V b(r[4], r[5], r[6], r[7]);
    io.outInt(a.xNumer); io.outInt(a.xDenom); io.outInt(a.yNumer); io.outInt(a.yDenom);
    io.outInt(b.xNumer); io.outInt(b.xDenom); io.outInt(b.yNumer); io.outInt(b.yDenom);
    io.outBool(a == b);
    io.outBool(a < b);
    io.outBool(b < a);
    int i0 = io.integer(0, 5);
    int i1 = io.integer(0, 5);
    int j0 = io.integer(0, 5);
    int j1 = io.integer(0, 5);
    E e0(i0, i1);
    E e1(j0, j1);
    io.outInt(e0.v[0]); io.outInt(e0.v[1]);
    io.outBool(e0 == e1);
    io.outBool(e0 < e1);
    io.outBool(e1 < e0);
}

ORACLE_CASE("CurveExtractor.invalidBounds")
{
    // The protected constructor's LogAssert(xBound > 1 && yBound > 1 &&
    // inputPixels != nullptr): throw parity over bounds 0..3 (the pixel
    // array always holds enough samples, the port's extra length check).
    int32_t xBound = io.integer(0, 3);
    int32_t yBound = io.integer(0, 3);
    ReplayExtractor extractor(xBound, yBound);
    io.outInt(xBound * yBound);
}
