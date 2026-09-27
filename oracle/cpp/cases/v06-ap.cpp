// Group 6 (arbitrary precision; the group is alphabetical): BSRational,
// APConversion, APInterval, BSPPolygon2, BSplineCurve, BSplineGeodesic,
// BSplineReduction (ArbitraryPrecision.h is an umbrella header with no code).
//
// BSRational is instantiated with UIntegerAP32. The port replaces the
// UInteger layer with bigint, so a BSRational is compared by VALUE and
// canonical form, never by storage layout: OutBSR emits the numerator and the
// denominator the way v05's OutBSN emits a BSNumber (sign, biased exponent,
// exponent, bit count, the 32-bit words of the odd integer or, above 8 words,
// the top three, the bottom two and a 48-bit FNV-1a digest of all words, then
// operator double and operator float of that BSNumber), followed by the
// BSRational's own operator double and operator float. OutBSN is copied from
// v05-ap.cpp.
//
// Everything is + - * / sqrt, comparisons and integer/bit manipulation except
// the std:: math overloads of BSRational that round-trip through a libm call;
// those live in the one tolerance case BSRational.std.libm. BSplineGeodesic
// is arithmetic only (its metric is the B-spline surface, the geodesic code
// calls only std::sqrt and a Gaussian-elimination inverse).
#define ORACLE_FAMILY "v06-ap"
#include "Oracle.h"

// BSPPolygon2 keeps its vertex/edge arrays and its BSP tree private. The
// cases compare the whole state (the edge array can be longer than the edge
// map, finding #169), so the header is included with 'private' opened up.
// Every standard and GTE header it uses is included first, so the macro
// reaches BSPPolygon2.h alone.
#include <Mathematics/Logger.h>
#include <Mathematics/Vector2.h>
#include <Mathematics/EdgeKey.h>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <utility>
#include <vector>
#define private public
#include <Mathematics/BSPPolygon2.h>
#undef private

#include <Mathematics/BSNumber.h>
#include <Mathematics/BSRational.h>
#include <Mathematics/UIntegerAP32.h>
#include <Mathematics/ArbitraryPrecision.h>
#include <Mathematics/APConversion.h>
#include <Mathematics/APInterval.h>
#include <Mathematics/BSplineCurve.h>
#include <Mathematics/BSplineSurface.h>
#include <Mathematics/BSplineGeodesic.h>
#include <Mathematics/BSplineReduction.h>

#include <array>
#include <cfenv>
#include <cstring>
#include <limits>
#include <string>

using namespace gte;
using BSN = BSNumber<UIntegerAP32>;
using BSR = BSRational<UIntegerAP32>;

namespace
{
    // ------------------------------------------------------------ doubles
    // FromBits, RawBits64, RawSign, Pow2, DrawDouble and DrawNormal are
    // copied from v05-ap.cpp.

    double FromBits(uint64_t b)
    {
        double d;
        std::memcpy(&d, &b, sizeof(d));
        return d;
    }

    uint64_t RawBits64(oracle::Ctx& io)
    {
        uint64_t b = 0;
        for (int i = 0; i < 4; ++i)
        {
            uint64_t piece = static_cast<uint64_t>(io.rawInteger(0, 0xFFFF));
            b = (b << 16) | piece;
        }
        return b;
    }

    double RawSign(oracle::Ctx& io)
    {
        return io.rawInteger(0, 1) != 0 ? -1.0 : 1.0;
    }

    double Pow2(int k)
    {
        return std::ldexp(1.0, k);
    }

    // A finite double over the whole binary64 range, drawn (not recorded):
    // small integers, powers of two from 2^-1074 to 2^1023, odd and all-ones
    // mantissas, arbitrary bit patterns, subnormals, signed zeros, the top
    // and bottom normal binades and binary32 values.
    double DrawDouble(oracle::Ctx& io)
    {
        double s = RawSign(io);
        int mode = io.rawInteger(0, 10);
        switch (mode)
        {
        case 0:
            return io.raw(-10.0, 10.0);
        case 1:
            return static_cast<double>(io.rawInteger(-20, 20));
        case 2:
            return s * Pow2(io.rawInteger(-1074, 1023));
        case 3:
        {
            double m = 2.0 * io.rawInteger(0, 1 << 19) + 1.0;
            int k = io.rawInteger(-1094, 1003);
            return s * std::ldexp(m, k);
        }
        case 4:
        {
            int m = io.rawInteger(1, 53);
            int k = io.rawInteger(-1000, 900);
            return s * std::ldexp(Pow2(m) - 1.0, k);
        }
        case 5:
        {
            double d = 1.0;
            for (int guard = 0; guard < 16; ++guard)
            {
                d = FromBits(RawBits64(io));
                if (std::isfinite(d)) { break; }
                d = 1.0;
            }
            return d;
        }
        case 6:
        {
            uint64_t t = RawBits64(io) & ((1ull << 52) - 1ull);
            if (t == 0) { t = 1; }
            return s * FromBits(t);
        }
        case 7:
            return s * 0.0;
        case 8:
        {
            double m = 1.0 + io.raw(0.0, 1.0);
            bool top = (io.rawInteger(0, 1) != 0);
            int k = (top ? io.rawInteger(1000, 1023) : io.rawInteger(-1022, -1000));
            return s * std::ldexp(m, k);
        }
        case 9:
            return static_cast<double>(static_cast<float>(io.raw(-1.0e6, 1.0e6)));
        default:
        {
            double m = 2.0 * io.rawInteger(0, (1 << 23) - 1) + 1.0;
            int k = io.rawInteger(-173, 104);
            return s * std::ldexp(m, k);
        }
        }
    }

    // A nonzero finite double (DrawDouble with the zeros redrawn; capped).
    double DrawNonzero(oracle::Ctx& io)
    {
        for (int guard = 0; guard < 32; ++guard)
        {
            double d = DrawDouble(io);
            if (d != 0.0) { return d; }
        }
        return 1.0;
    }

    // A moderate double: small lattice values, dyadics and uniform values in
    // [-8, 8], with the zeros kept. Rationals built from these stay small, so
    // chains of arithmetic do not explode the bit counts.
    double DrawModerate(oracle::Ctx& io)
    {
        switch (io.rawInteger(0, 3))
        {
        case 0: return static_cast<double>(io.rawInteger(-9, 9));
        case 1: return io.rawInteger(-64, 64) * 0.125;
        case 2: return io.raw(-8.0, 8.0);
        default: return RawSign(io) * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(-40, 40));
        }
    }

    // A normal double a = +-(1 + f) * 2^e with e in [lo, hi].
    double DrawNormal(oracle::Ctx& io, int lo, int hi)
    {
        double s = RawSign(io);
        bool dyadic = (io.rawInteger(0, 1) != 0);
        double m = (dyadic ? 1.0 + io.rawInteger(0, 15) / 16.0 : 1.0 + io.raw(0.0, 1.0));
        int k = io.rawInteger(lo, hi);
        return s * std::ldexp(m, k);
    }
}

namespace
{
    // ----------------------------------------------------------- emission

    // Emit a BSNumber by value and canonical form (copied from v05-ap.cpp).
    void OutBSN(oracle::Ctx& io, BSN const& x)
    {
        io.outInt(x.GetSign());
        io.outInt(x.GetBiasedExponent());
        io.outInt(x.GetExponent());
        auto const& u = x.GetUInteger();
        int32_t numBits = u.GetNumBits();
        io.outInt(numBits);
        int32_t numWords = (numBits + 31) / 32;
        auto const& bits = u.GetBits();
        if (numWords <= 8)
        {
            for (int32_t i = 0; i < numWords; ++i) { io.outInt(bits[i]); }
        }
        else
        {
            io.outInt(bits[numWords - 1]);
            io.outInt(bits[numWords - 2]);
            io.outInt(bits[numWords - 3]);
            io.outInt(bits[1]);
            io.outInt(bits[0]);
            uint64_t h = 0xCBF29CE484222325ull;
            for (int32_t i = 0; i < numWords; ++i)
            {
                for (int32_t k = 0; k < 4; ++k)
                {
                    h ^= static_cast<uint64_t>((bits[i] >> (8 * k)) & 0xFFu);
                    h *= 0x100000001B3ull;
                }
            }
            io.outInt(h & ((1ull << 48) - 1ull));
        }
        io.outReal(static_cast<double>(x));
        io.outReal(static_cast<double>(static_cast<float>(x)));
    }

    // Emit a valid BSRational: numerator, denominator, then the rational's
    // own operator double and operator float (Convert with FE_TONEAREST to
    // 53 and 24 bits, then the BSNumber conversion).
    void OutBSR(oracle::Ctx& io, BSR const& r)
    {
        OutBSN(io, r.GetNumerator());
        OutBSN(io, r.GetDenominator());
        io.outReal(static_cast<double>(r));
        io.outReal(static_cast<double>(static_cast<float>(r)));
    }

    // ------------------------------------------------------- BSRational draws

    // The half ulp of a at the binary64 grid (normal or subnormal), the
    // binary32 grid, used to put a rational on or next to a rounding tie.
    int HalfUlpExponent(double a, int grid)
    {
        if (grid == 0)
        {
            int ea = (a == 0.0 ? -1022 : std::ilogb(a));
            return (ea >= -1022 ? ea : -1022) - 53;
        }
        int ef = (a == 0.0 ? -126 : std::ilogb(a));
        return (ef >= -126 ? ef : -126) - 24;
    }

    // A BSRational built from recorded doubles. The mode is recorded first;
    // the replay reads the same doubles and applies the same operations.
    //   0  BSR(a)                         1 double, any finite
    //   1  BSR(BSN(a), BSN(b))            2 nonzero doubles
    //   2  BSR(a) / BSR(b)                2 doubles, b nonzero
    //   3  BSR(a) + BSR(b)                2 doubles (cancellation included)
    //   4  BSR(a) * BSR(b) - BSR(c)       3 moderate doubles
    //   5  BSR(int64 p, int64 q)          2 integer doubles, |p|,|q| < 2^45,
    //                                     often with a common odd factor
    //   6  BSR(a) + ldexp(BSR(hs), e) [+ ldexp(BSR(cs), e - j) / BSR(3)]
    //                                     a on the binary64 (normal or
    //                                     subnormal) or binary32 grid, then
    //                                     e, hs, j, cs: +-half ulp (a tie)
    //                                     and, for j > 0, a non-dyadic tie
    //                                     breaker
    BSR DrawBSR(oracle::Ctx& io, int lo, int hi)
    {
        int mode = io.integer(lo, hi);
        switch (mode)
        {
        case 0:
        {
            double a = io.given(DrawDouble(io));
            return BSR(a);
        }
        case 1:
        {
            double a = io.given(DrawNonzero(io));
            double b = io.given(DrawNonzero(io));
            return BSR(BSN(a), BSN(b));
        }
        case 2:
        {
            double a = io.given(DrawDouble(io));
            double b = io.given(DrawNonzero(io));
            return BSR(a) / BSR(b);
        }
        case 3:
        {
            double a = io.given(DrawDouble(io));
            double b = io.given(io.rawInteger(0, 3) == 0 ? -a : DrawDouble(io));
            return BSR(a) + BSR(b);
        }
        case 4:
        {
            double a = io.given(DrawModerate(io));
            double b = io.given(DrawModerate(io));
            double c = io.given(DrawModerate(io));
            return BSR(a) * BSR(b) - BSR(c);
        }
        case 5:
        {
            double k = static_cast<double>(2 * io.rawInteger(0, 20) + 1);
            double p = static_cast<double>(io.rawInteger(-1000000, 1000000));
            double q = static_cast<double>(io.rawInteger(1, 1000000));
            if (io.rawInteger(0, 2) == 0) { p = std::ldexp(p, io.rawInteger(0, 20)); }
            if (io.rawInteger(0, 2) == 0) { q = std::ldexp(q, io.rawInteger(0, 20)); }
            double pk = io.given(p * k);
            double qk = io.given(RawSign(io) * q * k);
            return BSR(static_cast<int64_t>(pk), static_cast<int64_t>(qk));
        }
        default:
        {
            int grid = io.rawInteger(0, 2);
            double a;
            if (grid == 0) { a = DrawNormal(io, -1021, 1022); }
            else if (grid == 1) { a = static_cast<double>(static_cast<float>(DrawNormal(io, -125, 126))); }
            else { a = RawSign(io) * std::ldexp(static_cast<double>(io.rawInteger(1, 1 << 20)), -1074); }
            a = io.given(a);
            int e = static_cast<int>(io.given(static_cast<double>(HalfUlpExponent(a, grid == 1 ? 1 : 0))));
            double hs = io.given(RawSign(io));
            int j = io.integer(0, 60);
            double cs = io.given(RawSign(io));
            BSR r = BSR(a) + std::ldexp(BSR(hs), e);
            if (j > 0)
            {
                r = r + std::ldexp(BSR(cs), e - j) / BSR(3);
            }
            return r;
        }
        }
    }

    // The second operand of a comparison or an arithmetic operation, drawn
    // relative to the first so that ties and near ties are reached.
    //   0  independent            1  y = x               2  y = -x
    //   3  y = x + BSR(e)         4  y = x * BSR(2^k)
    //   5  x and y the same value with different representations:
    //      x = BSR(p*k, q*k) and y = BSR(p, q) (k odd), x replaced
    //   6  y = x + BSR(1) / BSR(3 * 2^k) (a non-dyadic neighbour)
    void DrawPairBSR(oracle::Ctx& io, BSR& x, BSR& y)
    {
        x = DrawBSR(io, 0, 6);
        int mode = io.integer(0, 6);
        switch (mode)
        {
        case 0: y = DrawBSR(io, 0, 6); break;
        case 1: y = x; break;
        case 2: y = -x; break;
        case 3:
        {
            double e = io.given(RawSign(io) * Pow2(io.rawInteger(-80, 10)));
            y = x + BSR(e);
            break;
        }
        case 4:
        {
            double s = io.given(Pow2(io.rawInteger(-40, 40)));
            y = x * BSR(s);
            break;
        }
        case 5:
        {
            double k = static_cast<double>(2 * io.rawInteger(1, 20) + 1);
            double p = io.given(static_cast<double>(io.rawInteger(-5000, 5000)));
            double q = io.given(static_cast<double>(io.rawInteger(1, 5000)));
            double kk = io.given(k);
            x = BSR(static_cast<int64_t>(p * kk), static_cast<int64_t>(q * kk));
            y = BSR(static_cast<int64_t>(p), static_cast<int64_t>(q));
            break;
        }
        default:
        {
            double t = io.given(3.0 * Pow2(io.rawInteger(-60, 4)));
            y = x + BSR(1) / BSR(t);
            break;
        }
        }
    }

    // A 64-bit integer pattern (unrecorded), as v05's integer case draws
    // them: small values, powers of two, 2^k +- 1, extremes, random bits.
    uint64_t DrawIntegerBits(oracle::Ctx& io)
    {
        switch (io.rawInteger(0, 5))
        {
        case 0: return static_cast<uint64_t>(static_cast<int64_t>(io.rawInteger(-20, 20)));
        case 1: return 1ull << io.rawInteger(0, 63);
        case 2:
        {
            uint64_t p = 1ull << io.rawInteger(1, 63);
            return (io.rawInteger(0, 1) != 0 ? p + 1ull : p - 1ull);
        }
        case 3:
        {
            static uint64_t const extremes[8] = {
                0x7FFFFFFFull, 0x80000001ull, 0xFFFFFFFFull, 0xFFFFFFFF80000001ull,
                0x7FFFFFFFFFFFFFFFull, 0x8000000000000001ull, 0xFFFFFFFFFFFFFFFFull,
                (1ull << 53) + 1ull };
            return extremes[io.rawInteger(0, 7)];
        }
        default: return RawBits64(io);
        }
    }

    // Reduce a pattern to the range of the integer kind (0 int32_t,
    // 1 uint32_t, 2 int64_t, 3 uint64_t). INT32_MIN and INT64_MIN are
    // excluded: BSNumber evaluates '-number' for them (signed overflow).
    uint64_t ToKind(uint64_t bits, int kind)
    {
        if (kind == 0)
        {
            int32_t v = static_cast<int32_t>(static_cast<uint32_t>(bits));
            if (v == std::numeric_limits<int32_t>::min()) { v = v + 1; }
            return static_cast<uint64_t>(static_cast<int64_t>(v));
        }
        if (kind == 1) { return bits & 0xFFFFFFFFull; }
        if (kind == 2 && bits == 0x8000000000000000ull) { return 0x8000000000000001ull; }
        return bits;
    }

    void RecordBits(oracle::Ctx& io, uint64_t bits)
    {
        io.given(static_cast<double>(bits >> 32));
        io.given(static_cast<double>(bits & 0xFFFFFFFFull));
    }

    void RecordString(oracle::Ctx& io, std::string const& s)
    {
        io.given(static_cast<double>(s.size()));
        for (char c : s) { io.given(static_cast<double>(static_cast<unsigned char>(c))); }
    }
}

// ============================================================= BSRational

// BSRational(double) and BSRational(double, double): every finite double,
// +-inf and NaN (BSNumber's graceful exit: +-2^1024 and 0), zero and
// negative denominators (the sign moves to the numerator; this pair
// constructor does NOT normalize the denominator's exponent), and a zero
// denominator (LogAssert; throw parity).
ORACLE_CASE("BSRational.construct.double")
{
    int kind = io.integer(0, 3);
    auto special = [&io](double v)
    {
        int k = io.rawInteger(0, 15);
        if (k == 0) { return RawSign(io) * std::numeric_limits<double>::infinity(); }
        if (k == 1) { return std::numeric_limits<double>::quiet_NaN(); }
        return v;
    };
    double n = io.given(special(DrawDouble(io)));
    if (kind == 0)
    {
        BSR r(n);
        OutBSR(io, r);
        return;
    }
    double dv = (kind == 3 ? RawSign(io) * 0.0 : special(DrawNonzero(io)));
    if (std::isnan(dv)) { dv = 3.0; }
    double d = io.given(dv);
    BSR r(n, d);
    OutBSR(io, r);
}

// BSRational(float) and BSRational(float, float): binary32 values incl.
// subnormals, FLT_MAX, inf, NaN (recorded exactly as doubles).
ORACLE_CASE("BSRational.construct.float")
{
    int kind = io.integer(0, 3);
    auto drawFloat = [&io](bool nonzero)
    {
        float s = static_cast<float>(RawSign(io));
        switch (io.rawInteger(0, 7))
        {
        case 0: return s * std::numeric_limits<float>::infinity();
        case 1: return nonzero ? 3.0f : std::numeric_limits<float>::quiet_NaN();
        case 2: return s * std::ldexp(static_cast<float>(io.rawInteger(1, (1 << 23) - 1)), -149);
        case 3: return s * std::numeric_limits<float>::max();
        case 4: return s * std::numeric_limits<float>::min();
        case 5: return nonzero ? s * 0.5f : s * 0.0f;
        default:
        {
            float f = static_cast<float>(DrawDouble(io));
            return (nonzero && f == 0.0f ? 1.5f : f);
        }
        }
    };
    float n = static_cast<float>(io.given(static_cast<double>(drawFloat(false))));
    if (kind == 0)
    {
        BSR r(n);
        OutBSR(io, r);
        return;
    }
    float dv = (kind == 3 ? static_cast<float>(RawSign(io)) * 0.0f : drawFloat(true));
    float d = static_cast<float>(io.given(static_cast<double>(dv)));
    BSR r(n, d);
    OutBSR(io, r);
}

// BSRational(int32_t), (uint32_t), (int64_t), (uint64_t) and the four pair
// constructors. Each value is recorded as the high and low 32-bit words of
// its 64-bit pattern. A zero denominator asserts (throw parity).
ORACLE_CASE("BSRational.construct.integer")
{
    int kind = io.integer(0, 3);
    bool pair = io.boolean();
    uint64_t nb = ToKind(DrawIntegerBits(io), kind);
    RecordBits(io, nb);
    BSR r;
    if (!pair)
    {
        if (kind == 0) { r = BSR(static_cast<int32_t>(static_cast<uint32_t>(nb))); }
        else if (kind == 1) { r = BSR(static_cast<uint32_t>(nb)); }
        else if (kind == 2) { r = BSR(static_cast<int64_t>(nb)); }
        else { r = BSR(static_cast<uint64_t>(nb)); }
        OutBSR(io, r);
        return;
    }
    uint64_t db = ToKind(io.rawInteger(0, 9) == 0 ? 0ull : DrawIntegerBits(io), kind);
    RecordBits(io, db);
    if (kind == 0)
    {
        r = BSR(static_cast<int32_t>(static_cast<uint32_t>(nb)), static_cast<int32_t>(static_cast<uint32_t>(db)));
    }
    else if (kind == 1) { r = BSR(static_cast<uint32_t>(nb), static_cast<uint32_t>(db)); }
    else if (kind == 2) { r = BSR(static_cast<int64_t>(nb), static_cast<int64_t>(db)); }
    else { r = BSR(static_cast<uint64_t>(nb), static_cast<uint64_t>(db)); }
    OutBSR(io, r);
}

namespace
{
    // A BSNumber from recorded doubles: BSN(a), BSN(a) + BSN(b) or
    // BSN(a) * BSN(b) (mode recorded first).
    BSN DrawBSN(oracle::Ctx& io, bool nonzero)
    {
        int mode = io.integer(0, 2);
        double a = io.given(nonzero ? DrawNonzero(io) : DrawDouble(io));
        double b = io.given(nonzero ? DrawNonzero(io) : DrawDouble(io));
        if (mode == 0) { return BSN(a); }
        if (mode == 1)
        {
            BSN s = BSN(a) + BSN(b);
            return (nonzero && s.GetSign() == 0 ? BSN(a) : s);
        }
        return BSN(a) * BSN(b);
    }
}

// BSRational(BSNumber) and BSRational(BSNumber, BSNumber) with a nonzero
// numerator (the pair constructor moves the denominator's exponent into the
// numerator); a zero denominator asserts. A zero numerator with a pair is
// finding #168 and has its own deviation case.
ORACLE_CASE("BSRational.construct.bsnumber")
{
    int kind = io.integer(0, 2);
    BSN n = DrawBSN(io, kind != 0);
    if (kind == 0)
    {
        BSR r(n);
        OutBSR(io, r);
        return;
    }
    BSN d = (io.integer(0, 9) == 0 ? BSN(0) : DrawBSN(io, true));
    BSR r(n, d);
    OutBSR(io, r);
}

// Finding #168 (port fixed): BSRational(BSNumber, BSNumber) adjusts the
// numerator's biased exponent even when the numerator is zero, producing a
// zero with a nonzero biased exponent (visible through GetNumerator; the
// value, comparisons and arithmetic are unaffected). The denominator is
// drawn with a nonzero exponent, so every record deviates in the numerator's
// biased exponent.
ORACLE_CASE("BSRational.construct.bsnumber.zeroNumerator")
{
    double z = io.given(RawSign(io) * 0.0);
    double dv = DrawNonzero(io);
    for (int guard = 0; guard < 32 && BSN(dv).GetExponent() == 0; ++guard) { dv = DrawNonzero(io); }
    if (BSN(dv).GetExponent() == 0) { dv = 6.0; }
    double d = io.given(dv);
    BSR r = BSR(BSN(z), BSN(d));
    // The value is unaffected: equality with zero, order and arithmetic
    // agree with the canonical zero (emitted first, so they are compared
    // before the deviating biased exponent).
    io.outBool(r == BSR(0));
    io.outBool(BSR(0) == r);
    io.outBool(r < BSR(1));
    io.outBool(BSR(-1) < r);
    OutBSR(io, r + BSR(d));
    OutBSR(io, r);
}

namespace
{
    std::string Digits(oracle::Ctx& io, int n, bool nonzeroLead)
    {
        std::string s;
        for (int i = 0; i < n; ++i)
        {
            int lo = (i == 0 && nonzeroLead ? 1 : 0);
            s.push_back(static_cast<char>('0' + io.rawInteger(lo, 9)));
        }
        return s;
    }

    bool AllZeroDigits(std::string const& s)
    {
        for (char c : s)
        {
            if (c >= '1' && c <= '9') { return false; }
        }
        return true;
    }

    // A string for BSRational(std::string). Valid forms "x", "x.y", ".y",
    // "x." with an optional sign; x is "0" or has no leading zero, y may
    // have leading and trailing zeros; 1..60 digits on either side. Strings
    // whose value is zero are finding #168 (their own deviation case) and a
    // one-character non-digit integer part is finding #95 (its own
    // deviation case); both are redrawn here. Invalid forms (throw parity):
    // "", "+", "-", leading zeros in x, non-digits in x or y, a second '.',
    // exponent notation.
    std::string DrawRationalString(oracle::Ctx& io)
    {
        static char const* const junk = "abexyzE/:;= _,+-";
        for (int guard = 0; guard < 64; ++guard)
        {
            int signKind = io.rawInteger(0, 2);
            std::string sign = (signKind == 0 ? "" : (signKind == 1 ? "+" : "-"));
            int mode = io.rawInteger(0, 9);
            std::string x, y, s;
            int nx = (io.rawInteger(0, 3) == 0 ? io.rawInteger(20, 60) : io.rawInteger(1, 12));
            int ny = (io.rawInteger(0, 3) == 0 ? io.rawInteger(20, 60) : io.rawInteger(1, 12));
            x = (io.rawInteger(0, 4) == 0 ? std::string("0") : Digits(io, nx, true));
            y = Digits(io, ny, false);
            switch (mode)
            {
            case 0: case 1: s = x; break;
            case 2: case 3: case 4: s = x + "." + y; break;
            case 5: s = "." + y; break;
            case 6: s = x + "."; break;
            case 7:
            {
                // Invalid forms, including a second sign in front of a
                // multi-character integer part, which ConvertToInteger
                // rejects ("--5" leaves "-5").
                static char const* const forms[12] = { "", "+", "-", "1e5", "1.5e3", "2.5.1",
                    "--5", "+-2.5", "-+12", "++0.5", "-+0", "+-37.25" };
                return forms[io.rawInteger(0, 11)];
            }
            case 8:
            {
                // A leading zero in a multi-digit integer part.
                s = "0" + Digits(io, io.rawInteger(1, 6), false) + (io.rawInteger(0, 1) != 0 ? "." + y : "");
                break;
            }
            default:
            {
                // A non-digit inside the integer part (length >= 2) or the
                // fractional part.
                bool inFraction = (io.rawInteger(0, 1) != 0);
                std::string& t = (inFraction ? y : x);
                if (!inFraction && t.size() < 2) { t = "1" + t; }
                int pos = io.rawInteger(inFraction ? 0 : 1, static_cast<int>(t.size()) - 1);
                t[static_cast<size_t>(pos)] = junk[io.rawInteger(0, 15)];
                s = x + "." + y;
                break;
            }
            }
            s = sign + s;
            // Redraw the #168 zero strings and the #95 single non-digit
            // integer parts (neither can occur in modes 7..9 except by the
            // junk character, which leaves a multi-character part).
            if (mode <= 6)
            {
                bool zero = AllZeroDigits(s);
                if (zero) { continue; }
            }
            return s;
        }
        return "1.5";
    }
}

ORACLE_CASE("BSRational.construct.string")
{
    std::string s = DrawRationalString(io);
    RecordString(io, s);
    BSR r(s);
    OutBSR(io, r);
}

// Finding #168 (port fixed): the string constructor calls
// mNumerator.SetSign(sign) unconditionally, so every string whose value is
// zero ("0", "-0.0", "+.000", "0.", ...) produces a numerator with sign +-1
// and a zero integer. Only the signs are emitted: converting or comparing
// such a number reads integer words that do not exist.
ORACLE_CASE("BSRational.construct.string.zero")
{
    int signKind = io.rawInteger(0, 2);
    std::string sign = (signKind == 0 ? "" : (signKind == 1 ? "+" : "-"));
    std::string zeros(static_cast<size_t>(io.rawInteger(1, 5)), '0');
    std::string s;
    switch (io.rawInteger(0, 3))
    {
    case 0: s = "0"; break;
    case 1: s = "0." + zeros; break;
    case 2: s = "." + zeros; break;
    default: s = "0."; break;
    }
    s = sign + s;
    RecordString(io, s);
    BSR r(s);
    io.outInt(r.GetSign());
    io.outInt(r.GetNumerator().GetSign());
    io.outInt(r.GetNumerator().GetUInteger().GetNumBits());
}

// Finding #95 (port fixed), inherited from BSNumber::ConvertToInteger: an
// integer part of one non-digit character is not validated, so upstream
// reads 'c' - '0' as its value ("a.5" is 49.5); the port asserts. Every
// record deviates structurally.
ORACLE_CASE("BSRational.construct.string.singleChar")
{
    static char const* const junk = "abcxyzXYZ/:;= _eE,!#$%&()*<>?@[]^`{|}~";
    int signKind = io.rawInteger(0, 2);
    std::string s = (signKind == 0 ? "" : (signKind == 1 ? "+" : "-"));
    s.push_back(junk[io.rawInteger(0, 37)]);
    if (io.rawInteger(0, 1) != 0)
    {
        s += "." + Digits(io, io.rawInteger(1, 4), true);
    }
    RecordString(io, s);
    BSR r(s);
    OutBSR(io, r);
}

// == != < <= > >= in both orders, on independent pairs, equal values with
// equal and with different representations, negations, near neighbours.
ORACLE_CASE("BSRational.compare")
{
    BSR x, y;
    DrawPairBSR(io, x, y);
    io.outBool(x == y);
    io.outBool(x != y);
    io.outBool(x < y);
    io.outBool(x <= y);
    io.outBool(x > y);
    io.outBool(x >= y);
    io.outBool(y == x);
    io.outBool(y < x);
    io.outBool(y <= x);
    io.outInt(x.GetSign());
    io.outInt(y.GetSign());
}

// Unary + and -, binary + - * / in both orders, the compound assignments,
// Negate, SetSign (on nonzero values) and GetSign.
ORACLE_CASE("BSRational.arithmetic")
{
    BSR x, y;
    DrawPairBSR(io, x, y);
    OutBSR(io, +x);
    OutBSR(io, -x);
    OutBSR(io, x + y);
    OutBSR(io, y + x);
    OutBSR(io, x - y);
    OutBSR(io, y - x);
    OutBSR(io, x * y);
    OutBSR(io, y * x);
    if (y.GetSign() != 0) { OutBSR(io, x / y); }
    if (x.GetSign() != 0) { OutBSR(io, y / x); }
    BSR z = x;
    z += y;
    OutBSR(io, z);
    z -= x;
    OutBSR(io, z);
    z *= x;
    OutBSR(io, z);
    if (y.GetSign() != 0)
    {
        z /= y;
        OutBSR(io, z);
    }
    BSR w = y;
    w.Negate();
    OutBSR(io, w);
    if (w.GetSign() != 0)
    {
        int s = io.integer(0, 1) != 0 ? 1 : -1;
        w.SetSign(s);
        OutBSR(io, w);
    }
}

// Division by a zero BSRational asserts (throw parity), for x / 0 and x /= 0
// with a zero built in several ways.
ORACLE_CASE("BSRational.divide.zero")
{
    BSR x = DrawBSR(io, 0, 6);
    int how = io.integer(0, 3);
    BSR zero;
    if (how == 1) { zero = BSR(-0.0); }
    else if (how == 2) { zero = x - x; }
    else if (how == 3) { zero = BSR(0) * x; }
    if (io.boolean())
    {
        OutBSR(io, x / zero);
    }
    else
    {
        x /= zero;
        OutBSR(io, x);
    }
}

namespace
{
    int const gRoundingModes[5] = { FE_TONEAREST, FE_DOWNWARD, FE_TOWARDZERO, FE_UPWARD, 12345 };
}

// Convert(BSRational, precision, mode, BSNumber&), the BSRational& overload
// and Convert(BSRational, mode, double&) / (float&), for all four <cfenv>
// modes and an unsupported mode (LogError for a nonzero input), precision
// 53, 24, 1..8, 9..120, 121..300 and <= 0 (LogError). Inputs include exact
// binary64/binary32 ties, tie +- a non-dyadic breaker and subnormal grids.
// Everything is computed before anything is emitted (throw parity).
ORACLE_CASE("BSRational.convert")
{
    BSR x = DrawBSR(io, 0, 6);
    int pmode = io.rawInteger(0, 5);
    int p = 0;
    switch (pmode)
    {
    case 0: p = 53; break;
    case 1: p = 24; break;
    case 2: p = io.rawInteger(1, 8); break;
    case 3: p = io.rawInteger(9, 120); break;
    case 4: p = io.rawInteger(121, 300); break;
    default: p = (io.rawInteger(0, 3) == 0 ? io.rawInteger(-3, 0) : 64); break;
    }
    int precision = io.integer(p, p);
    int m = io.integer(0, 4);
    int mode = gRoundingModes[m];
    BSN outN;
    BSR outR;
    double outD = 0.0;
    float outF = 0.0f;
    Convert(x, precision, mode, outN);
    Convert(x, precision, mode, outR);
    Convert(x, mode, outD);
    Convert(x, mode, outF);
    OutBSN(io, outN);
    OutBSR(io, outR);
    io.outReal(outD);
    io.outReal(static_cast<double>(outF));
}

// Convert(BSRational, mode, double&) / (float&) and operator double / float
// where the result is subnormal or beyond the largest finite value.
// Upstream rounds twice: Convert to 53 (24) bits with the requested mode,
// then BSNumber's conversion rounds to nearest onto the subnormal grid, or
// to infinity. Preserved (a new upstream suspect, see the report): the
// double operator is up to 1 ulp off in the subnormal range, and the
// directed modes return results on the wrong side of x there (FE_UPWARD
// gives 0 for a positive x below 2^-1075) and +-inf instead of the largest
// finite value for FE_DOWNWARD / FE_TOWARDZERO beyond it. Compared bit for
// bit. x = sign * (p / q) * 2^e with small p and q (q odd: not dyadic) or
// q = 1 and a dyadic tail far below the leading bit.
ORACLE_CASE("BSRational.convert.subnormalAndOverflow")
{
    int range = io.integer(0, 3);
    double p = io.given(static_cast<double>(io.rawInteger(1, 1 << 20)));
    double q = io.given(io.rawInteger(0, 2) == 0 ? 1.0 : static_cast<double>(2 * io.rawInteger(1, 50) + 1));
    int e = 0;
    switch (range)
    {
    case 0: e = io.integer(-1100, -1040); break;   // binary64 subnormals
    case 1: e = io.integer(-175, -140); break;     // binary32 subnormals
    case 2: e = io.integer(1004, 1010); break;     // binary64 overflow
    default: e = io.integer(108, 114); break;      // binary32 overflow
    }
    int tail = io.integer(0, 80);
    double s = io.given(RawSign(io));
    int m = io.integer(0, 3);
    BSR x = std::ldexp(BSR(p) / BSR(q), e);
    if (tail > 0) { x = x + std::ldexp(BSR(1), e - tail); }
    x = BSR(s) * x;
    int mode = gRoundingModes[m];
    double d = 0.0;
    float f = 0.0f;
    Convert(x, mode, d);
    Convert(x, mode, f);
    io.outReal(d);
    io.outReal(static_cast<double>(f));
    io.outReal(static_cast<double>(x));
    io.outReal(static_cast<double>(static_cast<float>(x)));
}

// The std:: and gte:: overloads that are exact or correctly rounded: fabs,
// frexp and ldexp (exact on the rational), floor, ceil, sqrt, fmod,
// remainder (exact IEEE operations on the converted doubles), clamp,
// invsqrt, isign, saturate, sign, sqr, FMA, RobustSOP and RobustDOP.
ORACLE_CASE("BSRational.std.exact")
{
    BSR x = DrawBSR(io, 0, 6);
    BSR y = DrawBSR(io, 0, 6);
    BSR z = DrawBSR(io, 0, 6);
    BSR w = DrawBSR(io, 0, 6);
    int k = io.integer(-1100, 1100);
    OutBSR(io, std::fabs(x));
    int32_t exponent = 0;
    BSR fr = std::frexp(x, &exponent);
    OutBSR(io, fr);
    io.outInt(exponent);
    OutBSR(io, std::ldexp(x, k));
    OutBSR(io, std::floor(x));
    OutBSR(io, std::ceil(x));
    OutBSR(io, std::sqrt(x));
    OutBSR(io, std::fmod(x, y));
    OutBSR(io, std::remainder(x, y));
    OutBSR(io, gte::clamp(x, y, z));
    OutBSR(io, gte::invsqrt(x));
    io.outInt(gte::isign(x));
    OutBSR(io, gte::saturate(x));
    OutBSR(io, gte::sign(x));
    OutBSR(io, gte::sqr(x));
    OutBSR(io, FMA(x, y, z));
    OutBSR(io, RobustSOP(x, y, z, w));
    OutBSR(io, RobustDOP(x, y, z, w));
}

// std::remainder and std::fmod with quotients far beyond 2^53, where only an
// exact quotient gives the IEEE result (the port reuses BSNumber's bigint
// remainder), plus the infinite and zero divisors (NaN -> 0, x).
ORACLE_CASE("BSRational.std.remainder")
{
    double a = io.given(DrawNormal(io, -40, 1000));
    double b = io.given(DrawNormal(io, -60, 60));
    double c = io.given(io.rawInteger(0, 2) == 0 ? RawSign(io) * 0.0 : DrawNormal(io, 1023, 1023));
    BSR x = BSR(a) / BSR(3);
    BSR y(b);
    BSR u(c);
    OutBSR(io, std::remainder(x, y));
    OutBSR(io, std::remainder(y, x));
    OutBSR(io, std::fmod(x, y));
    OutBSR(io, std::remainder(y, u));
    OutBSR(io, std::fmod(y, u));
    OutBSR(io, std::remainder(y, u * BSR(4)));
}

// The overloads that round-trip through a C math library call: acos acosh
// asin asinh atan atanh atan2 cos cosh exp exp2 log log2 log10 pow sin sinh
// tan tanh (std::) and atandivpi atan2divpi cospi exp10 sinpi (gte::, built
// on std::atan, atan2, cos, exp, sin; the port's exp2 is Math.pow(2, x)).
// The MSVC runtime and V8 may differ by an ulp: compared with the default
// scaled tolerance 1e-12 on the converted double. The arguments are ratios
// p/q of small integers (not dyadic, so the conversion to double rounds) in
// [-12, 12], plus signed zeros.
ORACLE_CASE("BSRational.std.libm")
{
    BSR v[2];
    for (int i = 0; i < 2; ++i)
    {
        int mode = io.rawInteger(0, 3);
        double p = 0.0, q = 1.0;
        switch (mode)
        {
        case 0: p = static_cast<double>(io.rawInteger(-99, 99)); q = 100.0 + io.rawInteger(0, 2); break;
        case 1: p = static_cast<double>(io.rawInteger(-3600, 3600)); q = 300.0 + io.rawInteger(0, 3); break;
        case 2: p = static_cast<double>(io.rawInteger(10, 400)); q = 7.0 + io.rawInteger(0, 3); break;
        default: p = RawSign(io) * 0.0; q = 3.0; break;
        }
        double pp = io.given(p);
        double qq = io.given(q);
        v[i] = BSR(pp) / BSR(qq);
    }
    BSR const& x = v[0];
    BSR const& y = v[1];
    io.outReal(static_cast<double>(std::acos(x)));
    io.outReal(static_cast<double>(std::acosh(x)));
    io.outReal(static_cast<double>(std::asin(x)));
    io.outReal(static_cast<double>(std::asinh(x)));
    io.outReal(static_cast<double>(std::atan(x)));
    io.outReal(static_cast<double>(std::atanh(x)));
    io.outReal(static_cast<double>(std::atan2(y, x)));
    io.outReal(static_cast<double>(std::cos(x)));
    io.outReal(static_cast<double>(std::cosh(x)));
    io.outReal(static_cast<double>(std::exp(x)));
    io.outReal(static_cast<double>(std::exp2(x)));
    io.outReal(static_cast<double>(std::log(x)));
    io.outReal(static_cast<double>(std::log2(x)));
    io.outReal(static_cast<double>(std::log10(x)));
    io.outReal(static_cast<double>(std::pow(x, y)));
    io.outReal(static_cast<double>(std::sin(x)));
    io.outReal(static_cast<double>(std::sinh(x)));
    io.outReal(static_cast<double>(std::tan(x)));
    io.outReal(static_cast<double>(std::tanh(x)));
    io.outReal(static_cast<double>(gte::atandivpi(x)));
    io.outReal(static_cast<double>(gte::atan2divpi(y, x)));
    io.outReal(static_cast<double>(gte::cospi(x)));
    io.outReal(static_cast<double>(gte::exp10(x)));
    io.outReal(static_cast<double>(gte::sinpi(x)));
}

// =========================================================== APConversion

namespace
{
    using APC = APConversion<BSR>;

    // A nonnegative rational square, recorded: the mode, then two or three
    // doubles.
    //   0  (p/q)^2 for small integers (an exact square: a is rational)
    //   1  p/q, small integers (usually not a square)
    //   2  ldexp(p/q, k) with |k| <= 300 (huge and tiny exponents)
    //   3  a square of a dyadic double, 2^(2k) or an odd integer squared
    //   4  zero
    BSR DrawSquare(oracle::Ctx& io)
    {
        int mode = io.integer(0, 4);
        double p = io.given(static_cast<double>(io.rawInteger(1, 2000)));
        double q = io.given(static_cast<double>(io.rawInteger(1, 60)));
        double k = io.given(static_cast<double>(io.rawInteger(-300, 300)));
        BSR r = BSR(p) / BSR(q);
        switch (mode)
        {
        case 0: return r * r;
        case 1: return r;
        case 2: return std::ldexp(r, static_cast<int32_t>(k));
        case 3:
        {
            BSR s = std::ldexp(BSR(p), static_cast<int32_t>(k) / 2);
            return s * s;
        }
        default: return BSR(0);
        }
    }

    // Precision and iteration limits: small, moderate and large.
    int DrawPrecision(oracle::Ctx& io)
    {
        static double const table[10] = { 1, 2, 3, 5, 8, 16, 24, 53, 64, 100 };
        return static_cast<int>(io.given(table[io.rawInteger(0, 9)]));
    }

    uint32_t DrawMaxIterations(oracle::Ctx& io)
    {
        static double const table[8] = { 1, 1, 2, 3, 4, 8, 16, 32 };
        return static_cast<uint32_t>(io.given(table[io.rawInteger(0, 7)]));
    }

    // ---- Probes for finding #280 item 1. AmBBisectionExhausted replicates
    // the bisection loops of APConversion::EstimateAmB (upstream's own
    // control flow) and reports that one ran out of iterations; the helpers
    // are copies of APConversion's private PreprocessSqr, GetMinOfSqrt and
    // GetMaxOfSqrt. CmpAmB is an exact comparison of t with sqrt(a2) -
    // sqrt(b2). AmBDefective is the observable symptom the port fixes: an
    // exhausted bisection, aSqr >= bSqr, and upstream's bracket misses a - b.

    void ProbePreprocess(BSR const& aSqr, BSR& rSqr, int32_t& exponentA)
    {
        int32_t exponentASqr;
        rSqr = std::frexp(aSqr, &exponentASqr);
        if (exponentASqr & 1)
        {
            exponentA = (exponentASqr - 1) / 2;
            rSqr = std::ldexp(rSqr, 1);
        }
        else
        {
            exponentA = exponentASqr / 2;
        }
    }

    BSR ProbeMinOfSqrt(BSR const& rSqr, int32_t exponent)
    {
        double lowerRSqr = 0.0;
        Convert(rSqr, FE_DOWNWARD, lowerRSqr);
        BSR aMin = std::nextafter(std::sqrt(lowerRSqr), -std::numeric_limits<double>::max());
        return std::ldexp(aMin, exponent);
    }

    BSR ProbeMaxOfSqrt(BSR const& rSqr, int32_t exponent)
    {
        double upperRSqr = 0.0;
        Convert(rSqr, FE_UPWARD, upperRSqr);
        BSR aMax = std::nextafter(std::sqrt(upperRSqr), +std::numeric_limits<double>::max());
        return std::ldexp(aMax, exponent);
    }

    // Returns true when a bisection loop of EstimateAmB runs out of iterations.
    // Only the bisection part is replicated; it decides staleness.
    bool AmBBisectionExhausted(BSR const& aSqr, BSR const& bSqr, int32_t precision, uint32_t maxIterations)
    {
        BSR const zero(0), three(3), five(5);
        BSR threshold = std::ldexp(BSR(1), -precision);
        BSR a2tb2 = aSqr * bSqr, a2pb2 = aSqr + bSqr, a2mb2 = aSqr - bSqr;
        BSR a2mb2Sqr = a2mb2 * a2mb2;
        BSR twoa2pb2 = std::ldexp(a2pb2, 1);
        BSR uSqr, vSqr;
        int32_t exponentA, exponentB;
        ProbePreprocess(aSqr, uSqr, exponentA);
        ProbePreprocess(bSqr, vSqr, exponentB);
        BSR signSecDer = a2mb2Sqr - five * a2tb2;
        if (signSecDer > zero)
        {
            BSR tMin = ProbeMinOfSqrt(uSqr, exponentA) - ProbeMaxOfSqrt(vSqr, exponentB);
            if (tMin < zero) { tMin = zero; }
            if (three * tMin * tMin - a2pb2 >= zero) { return false; }
            BSR tMax = ProbeMaxOfSqrt(uSqr, exponentA) - ProbeMinOfSqrt(vSqr, exponentB);

            for (uint32_t iterate = 1; iterate <= maxIterations; ++iterate)
            {
                if (tMax - tMin < threshold) { return false; }
                BSR tMid = std::ldexp(tMin + tMax, -1);
                BSR tMidSqr = tMid * tMid;
                if (three * tMidSqr - a2pb2 >= zero)
                {
                    BSR f = tMidSqr * (tMidSqr - twoa2pb2) + a2mb2Sqr;
                    if (f >= zero) { return false; }
                    tMax = tMid;
                    Convert(tMax, 2 * precision, FE_UPWARD, tMax);
                }
                else
                {
                    tMin = tMid;
                    Convert(tMin, 2 * precision, FE_DOWNWARD, tMin);

                }
            }
            return true;
        }
        if (signSecDer < zero)
        {
            BSR tMax = ProbeMaxOfSqrt(uSqr, exponentA) - ProbeMinOfSqrt(vSqr, exponentB);
            if (three * tMax * tMax - a2pb2 <= zero) { return false; }
            BSR tMin = ProbeMinOfSqrt(uSqr, exponentA) - ProbeMaxOfSqrt(vSqr, exponentB);

            for (uint32_t iterate = 1; iterate <= maxIterations; ++iterate)
            {
                if (tMax - tMin < threshold) { return false; }
                BSR tMid = std::ldexp(tMin + tMax, -1);
                BSR tMidSqr = tMid * tMid;
                if (three * tMidSqr - a2pb2 <= zero)
                {
                    BSR f = tMidSqr * (tMidSqr - twoa2pb2) + a2mb2Sqr;
                    if (f <= zero) { return false; }
                    tMin = tMid;
                    Convert(tMin, 2 * precision, FE_DOWNWARD, tMin);
                }
                else
                {
                    tMax = tMid;
                    Convert(tMax, 2 * precision, FE_UPWARD, tMax);

                }
            }
            return true;
        }
        return false;
    }

    // The exact sign of t - (sqrt(a2) - sqrt(b2)) for a2, b2 >= 0: with
    // u = t + b, u < 0 exactly when t < 0 and t^2 > b2; otherwise
    // sign(u - a) = sign(2*t*b - r) with r = a2 - b2 - t^2.
    int CmpAmB(BSR const& t, BSR const& a2, BSR const& b2)
    {
        BSR const zero(0);
        BSR tSqr = t * t;
        if (t < zero && tSqr > b2) { return -1; }
        BSR r = a2 - b2 - tSqr;
        int sx = (b2.GetSign() == 0 ? 0 : t.GetSign());
        int sr = r.GetSign();
        if (sx == 0) { return -sr; }
        if (sr == 0) { return sx; }
        if (sx != sr) { return sx; }
        BSR x2 = BSR(4) * tSqr * b2;
        BSR r2 = r * r;
        int c = (x2 < r2 ? -1 : (r2 < x2 ? 1 : 0));
        return sx > 0 ? c : -c;
    }

    bool AmBDefective(BSR const& aSqr, BSR const& bSqr, int32_t precision, uint32_t maxIterations)
    {
        if (aSqr < bSqr || !AmBBisectionExhausted(aSqr, bSqr, precision, maxIterations))
        {
            return false;
        }
        try
        {
            APConversion<BSR> apc(precision, maxIterations);
            BSR tMin, tMax;
            apc.EstimateAmB(aSqr, bSqr, tMin, tMax);
            return CmpAmB(tMin, aSqr, bSqr) > 0 || CmpAmB(tMax, aSqr, bSqr) < 0;
        }
        catch (std::exception const&)
        {
            return false;
        }
    }

    // aSqr/bSqr near the ratio (7 + 3*sqrt(5))/2 at which f''(a - b) = 0,
    // where both bisection loops of EstimateAmB run: a continued-fraction
    // convergent P/Q of that ratio (terms 6; 1, 5, 1, 5, ...) scaled by an
    // integer, recorded as P*s and Q*s.
    void DrawNearRatio(oracle::Ctx& io, BSR& aSqr, BSR& bSqr)
    {
        int n = io.rawInteger(1, 26);
        double p0 = 1.0, q0 = 0.0, p1 = 6.0, q1 = 1.0;
        for (int i = 1; i <= n; ++i)
        {
            double term = (i % 2 == 1 ? 1.0 : 5.0);
            double p2 = term * p1 + p0, q2 = term * q1 + q0;
            p0 = p1; q0 = q1; p1 = p2; q1 = q2;
        }
        double s = static_cast<double>(io.rawInteger(1, 7));
        double P = io.given(p1 * s);
        double Q = io.given(q1 * s);
        aSqr = BSR(P);
        bSqr = BSR(Q);
    }
}

// EstimateSqrt(aSqr, aMin, aMax) and EstimateSqrt(aSqr, a) over exact
// squares, non-squares, huge and tiny exponents, zero (finding #280 item 3:
// [0, 2^-1074]) and a negative input (the NaN from std::sqrt becomes a zero
// rational and aSqr / aMax asserts; throw parity).
ORACLE_CASE("APConversion.estimateSqrt")
{
    BSR aSqr = DrawSquare(io);
    if (io.integer(0, 15) == 0) { aSqr = -aSqr - BSR(1); }
    int32_t precision = DrawPrecision(io);
    uint32_t maxIterations = DrawMaxIterations(io);
    APC apc(precision, maxIterations);
    BSR aMin, aMax, a;
    uint32_t n0 = apc.EstimateSqrt(aSqr, aMin, aMax);
    uint32_t n1 = apc.EstimateSqrt(aSqr, a);
    io.outInt(n0);
    OutBSR(io, aMin);
    OutBSR(io, aMax);
    io.outInt(n1);
    OutBSR(io, a);
}

// EstimateApB(aSqr, bSqr, tMin, tMax): bounds on sqrt(aSqr) + sqrt(bSqr).
ORACLE_CASE("APConversion.estimateApB")
{
    BSR aSqr = DrawSquare(io);
    BSR bSqr = DrawSquare(io);
    int32_t precision = DrawPrecision(io);
    uint32_t maxIterations = DrawMaxIterations(io);
    APC apc(precision, maxIterations);
    BSR tMin, tMax;
    uint32_t n = apc.EstimateApB(aSqr, bSqr, tMin, tMax);
    io.outInt(n);
    OutBSR(io, tMin);
    OutBSR(io, tMax);
}

namespace
{
    // Inputs for EstimateAmB with aSqr >= bSqr (the undocumented
    // precondition of finding #280 item 2; its violation has its own case):
    //   0  two drawn squares, swapped into order
    //   1  near the ratio where f''(a - b) = 0 (both bisection loops)
    //   2  aSqr = bSqr, and bSqr = 0
    void DrawAmB(oracle::Ctx& io, BSR& aSqr, BSR& bSqr)
    {
        int mode = io.integer(0, 2);
        if (mode == 0)
        {
            aSqr = DrawSquare(io);
            bSqr = DrawSquare(io);
            if (aSqr < bSqr) { std::swap(aSqr, bSqr); }
        }
        else if (mode == 1)
        {
            DrawNearRatio(io, aSqr, bSqr);
        }
        else
        {
            aSqr = DrawSquare(io);
            bSqr = (io.boolean() ? aSqr : BSR(0));
        }
    }
}

// EstimateAmB(aSqr, bSqr, tMin, tMax) on inputs where upstream is sound:
// records on which a bisection loop runs out of iterations and upstream's
// bracket then misses a - b (finding #280 item 1, fixed in the port; the
// probe AmBDefective is exactly the port's replacement condition) are
// redrawn; the rest (every branch, early returns, exhausted bisections with
// a valid bracket, Newton exhaustion) is compared. aSqr = bSqr = 0 reaches
// the LogError on a zero second derivative (throw parity).
ORACLE_CASE("APConversion.estimateAmB")
{
    BSR aSqr, bSqr;
    int32_t precision = 1;
    uint32_t maxIterations = 1;
    for (int guard = 0; guard < 64; ++guard)
    {
        oracle::Ctx probe = io;
        DrawAmB(probe, aSqr, bSqr);
        precision = DrawPrecision(probe);
        maxIterations = DrawMaxIterations(probe);
        if (!AmBDefective(aSqr, bSqr, precision, maxIterations) || guard == 63)
        {
            io = probe;
            break;
        }
        io.rawInteger(0, 1);
    }
    APC apc(precision, maxIterations);
    BSR tMin, tMax;
    uint32_t n = apc.EstimateAmB(aSqr, bSqr, tMin, tMax);
    io.outInt(n);
    OutBSR(io, tMin);
    OutBSR(io, tMax);
}

// EstimateAmB with aSqr < bSqr, violating the undocumented precondition
// (finding #280 item 2, preserved): the f'' > 0 block clamps tMin to zero
// and returns the inverted bracket [0, aMax - bMin] at once; the f'' < 0
// block runs as usual. Compared bit for bit (both sides preserve it).
ORACLE_CASE("APConversion.estimateAmB.aLessThanB")
{
    BSR aSqr, bSqr;
    int32_t precision = 1;
    uint32_t maxIterations = 1;
    for (int guard = 0; guard < 64; ++guard)
    {
        oracle::Ctx probe = io;
        aSqr = DrawSquare(probe);
        bSqr = DrawSquare(probe);
        precision = DrawPrecision(probe);
        maxIterations = DrawMaxIterations(probe);
        if (aSqr < bSqr || guard == 63)
        {
            io = probe;
            break;
        }
        io.rawInteger(0, 1);
    }
    APC apc(precision, maxIterations);
    BSR tMin, tMax;
    uint32_t n = apc.EstimateAmB(aSqr, bSqr, tMin, tMax);
    io.outInt(n);
    OutBSR(io, tMin);
    OutBSR(io, tMax);
}

// Finding #280 item 1 (port fixed): when a bisection loop of EstimateAmB
// runs out of iterations, upstream runs Newton's method from outside its
// basin (after an update in the rounding branch also with a stale square)
// and can return a bracket that misses a - b; the port returns the
// bisection bracket instead, exactly when upstream's bracket is wrong.
// Inputs near the ratio (7 + 3*sqrt(5))/2 with few iterations; only records
// on which upstream's bracket is wrong (AmBDefective: the observable
// symptom, checked exactly) are kept (capped; the fallback is the last
// candidate).
ORACLE_CASE("APConversion.estimateAmB.bisectionExhausted")
{
    BSR aSqr, bSqr;
    int32_t precision = 1;
    uint32_t maxIterations = 1;
    for (int guard = 0; guard < 256; ++guard)
    {
        oracle::Ctx probe = io;
        DrawNearRatio(probe, aSqr, bSqr);
        precision = static_cast<int32_t>(probe.given(static_cast<double>(probe.rawInteger(30, 110))));
        maxIterations = static_cast<uint32_t>(probe.given(static_cast<double>(probe.rawInteger(1, 4))));
        if (AmBDefective(aSqr, bSqr, precision, maxIterations) || guard == 255)
        {
            io = probe;
            break;
        }
        io.rawInteger(0, 1);
    }
    APC apc(precision, maxIterations);
    BSR tMin, tMax;
    uint32_t n = apc.EstimateAmB(aSqr, bSqr, tMin, tMax);
    io.outInt(n);
    OutBSR(io, tMin);
    OutBSR(io, tMax);
}

// Estimate(QFN1 q, qMin, qMax) and Estimate(q, qEstimate) for
// q = x + y*sqrt(d): d = 0 or y = 0 (no iteration, [x, x]), y > 0, y < 0 (the
// interval is reflected), and d < 0 (EstimateSqrt of a negative number;
// throw parity).
ORACLE_CASE("APConversion.estimate")
{
    BSR x = DrawBSR(io, 4, 5);
    int ymode = io.integer(0, 3);
    BSR y = (ymode == 0 ? BSR(0) : DrawBSR(io, 4, 5));
    int dmode = io.integer(0, 5);
    BSR d = (dmode == 0 ? BSR(0) : DrawSquare(io));
    if (dmode == 5) { d = -d - BSR(1); }
    int32_t precision = DrawPrecision(io);
    uint32_t maxIterations = DrawMaxIterations(io);
    APC apc(precision, maxIterations);
    QFNumber<BSR, 1> q(x, y, d);
    BSR qMin, qMax, qEstimate;
    uint32_t n0 = apc.Estimate(q, qMin, qMax);
    uint32_t n1 = apc.Estimate(q, qEstimate);
    io.outInt(n0);
    OutBSR(io, qMin);
    OutBSR(io, qMax);
    io.outInt(n1);
    OutBSR(io, qEstimate);
}

// The constructor, SetPrecision, SetMaxIterations and the getters, with
// invalid arguments (precision <= 0, maxIterations = 0: LogAssert, throw
// parity); the new threshold is observed through EstimateSqrt(2).
ORACLE_CASE("APConversion.accessors")
{
    int32_t p0 = io.integer(-1, 40);
    uint32_t m0 = static_cast<uint32_t>(io.integer(0, 6));
    int32_t p1 = io.integer(-1, 60);
    uint32_t m1 = static_cast<uint32_t>(io.integer(0, 6));
    APC apc(p0, m0);
    int32_t g0 = apc.GetPrecision();
    uint32_t h0 = apc.GetMaxIterations();
    BSR aMin0, aMax0;
    uint32_t n0 = apc.EstimateSqrt(BSR(2), aMin0, aMax0);
    apc.SetPrecision(p1);
    apc.SetMaxIterations(m1);
    BSR aMin1, aMax1;
    uint32_t n1 = apc.EstimateSqrt(BSR(2), aMin1, aMax1);
    io.outInt(g0);
    io.outInt(h0);
    io.outInt(n0);
    OutBSR(io, aMin0);
    OutBSR(io, aMax0);
    io.outInt(apc.GetPrecision());
    io.outInt(apc.GetMaxIterations());
    io.outInt(n1);
    OutBSR(io, aMin1);
    OutBSR(io, aMax1);
}

// ============================================================= APInterval

namespace
{
    using API = APInterval<BSR>;

    // Emit an interval endpoint. A finite endpoint is a valid BSRational
    // (OutBSR). An infinity sentinel of finding #280 (SetSign(+-2) on a zero
    // rational, possibly multiplied by a finite number) is emitted as its
    // sign, its numerator's sign and biased exponent and its denominator;
    // not its numerator's bit count and exponent (upstream's UInteger::Mul
    // of a zero-bit operand leaves a count of zero words, the bigint port
    // has no words), and not its conversions (they read words that do not
    // exist).
    void OutAPE(oracle::Ctx& io, BSR const& e)
    {
        int32_t s = e.GetSign();
        io.outInt(s);
        if (s >= 2 || s <= -2)
        {
            io.outInt(e.GetNumerator().GetSign());
            io.outInt(e.GetNumerator().GetBiasedExponent());
            OutBSN(io, e.GetDenominator());
            return;
        }
        OutBSR(io, e);
    }

    void OutAPI(oracle::Ctx& io, API const& w)
    {
        OutAPE(io, w[0]);
        OutAPE(io, w[1]);
    }

    // A scalar endpoint: a moderate rational (DrawBSR modes 4..5) or zero.
    BSR DrawAPScalar(oracle::Ctx& io)
    {
        return (io.integer(0, 5) == 0 ? BSR(0) : DrawBSR(io, 4, 5));
    }

    // An interval [e0, e1] with e0 <= e1, recorded. The shapes reach every
    // sign branch of operator* and operator/:
    //   0 [a, b] general  1 [0, b]  2 [a, 0]  3 [a, a]  4 [a, b] with
    //   a < 0 < b  5 [0, 0]  6 both positive  7 both negative
    API DrawAPI(oracle::Ctx& io, int hi)
    {
        int shape = io.integer(0, hi);
        BSR a = DrawBSR(io, 4, 5);
        BSR b = DrawBSR(io, 4, 5);
        BSR lo = (a < b ? a : b), up = (a < b ? b : a);
        BSR fa = std::fabs(a), fb = std::fabs(b) + BSR(1);
        switch (shape)
        {
        case 0: return API(lo, up);
        case 1: return API(BSR(0), fb);
        case 2: return API(-fb, BSR(0));
        case 3: return API(a, a);
        case 4: return API(-fa - BSR(1), fb);
        case 5: return API(BSR(0), BSR(0));
        case 6: return API(fa + BSR(1), fa + fb + BSR(1));
        default: return API(-fa - fb - BSR(1), -fa - BSR(1));
        }
    }
}

// The leaf operations Add, Sub, Mul, Div(u, v) (Div by zero returns Reals).
ORACLE_CASE("APInterval.leaf")
{
    BSR u = DrawAPScalar(io);
    BSR v = DrawAPScalar(io);
    OutAPI(io, API::Add(u, v));
    OutAPI(io, API::Sub(u, v));
    OutAPI(io, API::Mul(u, v));
    OutAPI(io, API::Div(u, v));
    OutAPI(io, API::Div(v, u));
}

// The internal four-argument Add, Sub, Mul, Mul2, Div, and Reciprocal,
// ReciprocalDown, ReciprocalUp (nonzero arguments), Reals.
ORACLE_CASE("APInterval.internal")
{
    BSR u0 = DrawAPScalar(io);
    BSR u1 = DrawAPScalar(io);
    BSR v0 = DrawBSR(io, 4, 5);
    BSR v1 = DrawBSR(io, 4, 5);
    if (v0.GetSign() == 0) { v0 = BSR(3); }
    if (v1.GetSign() == 0) { v1 = BSR(-5); }
    OutAPI(io, API::Add(u0, u1, v0, v1));
    OutAPI(io, API::Sub(u0, u1, v0, v1));
    OutAPI(io, API::Mul(u0, u1, v0, v1));
    OutAPI(io, API::Mul2(u0, u1, v0, v1));
    OutAPI(io, API::Div(u0, u1, v0, v1));
    OutAPI(io, API::Reciprocal(v0, v1));
    OutAPI(io, API::ReciprocalDown(v0));
    OutAPI(io, API::ReciprocalUp(v1));
    OutAPI(io, API::Reals());
}

// Every non-class operator on interval/interval, interval/scalar and
// scalar/interval, the compound assignments, operator[], GetEndpoints and
// the constructors. Division by an interval with a zero endpoint returns
// the infinity sentinels of finding #280 (a sentinel only ever appears as
// the right operand of a product here: a sentinel on the left makes
// UIntegerALU32::Mul read word 0 of an empty array); their negation is
// emitted too. A [0, 0] divisor asserts and has its own case.
ORACLE_CASE("APInterval.operators")
{
    API u = DrawAPI(io, 7);
    API v = DrawAPI(io, 7);
    BSR s = DrawAPScalar(io);
    bool vZero = (v[0].GetSign() == 0 && v[1].GetSign() == 0);
    bool uZero = (u[0].GetSign() == 0 && u[1].GetSign() == 0);

    API d;
    OutAPI(io, d);
    OutAPI(io, API(s));
    std::array<BSR, 2> ends = u.GetEndpoints();
    OutAPI(io, API(ends));
    API c(u);
    OutAPI(io, c);
    OutAPE(io, u[0]);
    OutAPE(io, u[1]);

    OutAPI(io, +u);
    OutAPI(io, -u);
    OutAPI(io, s + u);
    OutAPI(io, u + s);
    OutAPI(io, u + v);
    OutAPI(io, s - u);
    OutAPI(io, u - s);
    OutAPI(io, u - v);
    OutAPI(io, s * u);
    OutAPI(io, u * s);
    OutAPI(io, u * v);
    OutAPI(io, v * u);
    OutAPI(io, u / s);
    if (!vZero)
    {
        API q = u / v;
        OutAPI(io, q);
        OutAPI(io, -q);
        OutAPI(io, s / v);
        OutAPI(io, -(s / v));
    }
    if (!uZero)
    {
        OutAPI(io, v / u);
    }

    API w = u;
    w += s;
    OutAPI(io, w);
    w += v;
    OutAPI(io, w);
    w -= s;
    OutAPI(io, w);
    w -= v;
    OutAPI(io, w);
    w *= s;
    OutAPI(io, w);
    w *= v;
    OutAPI(io, w);
    w /= s;
    OutAPI(io, w);
    bool wFinite = (std::abs(w[0].GetSign()) < 2 && std::abs(w[1].GetSign()) < 2);
    if (!vZero && wFinite)
    {
        w /= v;
        OutAPI(io, w);
    }
}

// Division by the interval [0, 0] (finding #280, preserved): the documented
// result is Reals(), but the v[0] == 0 branch calls ReciprocalDown(v[1]) and
// 1/0 asserts. Interval / interval, scalar / interval and /= all throw.
ORACLE_CASE("APInterval.divide.zeroInterval")
{
    API u = DrawAPI(io, 7);
    BSR s = DrawAPScalar(io);
    int how = io.integer(0, 2);
    API zero(BSR(0), BSR(0));
    if (how == 0) { OutAPI(io, u / zero); }
    else if (how == 1) { OutAPI(io, s / zero); }
    else
    {
        u /= zero;
        OutAPI(io, u);
    }
}

// ============================================================ BSPPolygon2

namespace
{
    using Poly = BSPPolygon2<double>;

    // An edge soup: vertex coordinates and directed edges (index pairs).
    struct Soup
    {
        std::vector<Vector2<double>> v;
        std::vector<std::array<int, 2>> e;
    };

    void AddLoop(Soup& s, std::vector<Vector2<double>> const& loop)
    {
        int base = static_cast<int>(s.v.size());
        int n = static_cast<int>(loop.size());
        for (auto const& p : loop) { s.v.push_back(p); }
        for (int i = 0; i < n; ++i) { s.e.push_back({ base + i, base + (i + 1) % n }); }
    }

    // A lattice shape (unrecorded), counterclockwise so that the interior is
    // on the negative side of the directed edges:
    //   0 rectangle  1 triangle  2 L  3 U  4 square with a square hole
    //   5 star-shaped polygon from 5..7 lattice points sorted by angle
    //   6 diamond
    Soup DrawShape(oracle::Ctx& io, int shape, bool uniform)
    {
        Soup s;
        auto P = [](double x, double y) { return Vector2<double>{ x, y }; };
        switch (shape)
        {
        case 0:
        {
            double w = io.rawInteger(1, 4), h = io.rawInteger(1, 4);
            AddLoop(s, { P(0, 0), P(w, 0), P(w, h), P(0, h) });
            break;
        }
        case 1:
        {
            for (int guard = 0; guard < 32; ++guard)
            {
                Vector2<double> a = P(io.rawInteger(-3, 3), io.rawInteger(-3, 3));
                Vector2<double> b = P(io.rawInteger(-3, 3), io.rawInteger(-3, 3));
                Vector2<double> c = P(io.rawInteger(-3, 3), io.rawInteger(-3, 3));
                double det = DotPerp(b - a, c - a);
                if (det != 0.0 || guard == 31)
                {
                    if (det == 0.0) { a = P(0, 0); b = P(2, 0); c = P(0, 2); }
                    if (det < 0.0) { std::swap(b, c); }
                    AddLoop(s, { a, b, c });
                    break;
                }
            }
            break;
        }
        case 2:
            AddLoop(s, { P(0, 0), P(3, 0), P(3, 1), P(1, 1), P(1, 3), P(0, 3) });
            break;
        case 3:
            AddLoop(s, { P(0, 0), P(3, 0), P(3, 3), P(2, 3), P(2, 1), P(1, 1), P(1, 3), P(0, 3) });
            break;
        case 4:
            AddLoop(s, { P(0, 0), P(4, 0), P(4, 4), P(0, 4) });
            AddLoop(s, { P(1, 1), P(1, 3), P(3, 3), P(3, 1) });
            break;
        case 5:
        {
            int n = io.rawInteger(5, 7);
            std::vector<Vector2<double>> pts;
            for (int i = 0; i < n; ++i)
            {
                Vector2<double> p = P(io.rawInteger(-3, 3), io.rawInteger(-3, 3));
                bool dup = false;
                for (auto const& q : pts) { dup = dup || (q == p); }
                if (!dup) { pts.push_back(p); }
            }
            // Sorting by angle about a point inside the hull (the centroid,
            // nudged off the lattice) gives a simple star-shaped polygon; a
            // fixed point, possibly outside the hull, on a quarter of the
            // records gives self-intersecting soups (garbage in: both sides
            // must still agree).
            Vector2<double> c = P(0.1, 0.05);
            if (io.rawInteger(0, 3) != 0 && !pts.empty())
            {
                c = P(0.01, 0.003);
                for (auto const& p : pts) { c += p / static_cast<double>(pts.size()); }
            }
            std::sort(pts.begin(), pts.end(), [&c](Vector2<double> const& a, Vector2<double> const& b)
                { return std::atan2(a[1] - c[1], a[0] - c[0]) < std::atan2(b[1] - c[1], b[0] - c[0]); });
            if (pts.size() < 3) { pts = { P(0, 0), P(2, 0), P(0, 2) }; }
            AddLoop(s, pts);
            break;
        }
        default:
        {
            double r = io.rawInteger(1, 3);
            AddLoop(s, { P(0, -r), P(r, 0), P(0, r), P(-r, 0) });
            break;
        }
        }
        // Scale, translate; the uniform mode perturbs every coordinate.
        double scale = io.rawInteger(1, 2);
        double tx = io.rawInteger(-3, 3), ty = io.rawInteger(-3, 3);
        for (auto& p : s.v)
        {
            p[0] = scale * p[0] + tx;
            p[1] = scale * p[1] + ty;
            if (uniform)
            {
                p[0] += io.raw(-0.3, 0.3);
                p[1] += io.raw(-0.3, 0.3);
            }
        }
        return s;
    }

    void RecordSoup(oracle::Ctx& io, Soup const& s)
    {
        io.given(static_cast<double>(s.v.size()));
        for (auto const& p : s.v) { io.given(p[0]); io.given(p[1]); }
        io.given(static_cast<double>(s.e.size()));
        for (auto const& e : s.e) { io.given(e[0]); io.given(e[1]); }
    }

    // Insert the soup edge by edge (vertices first, deduplicated by the
    // polygon's vertex map), then Finalize.
    void BuildPoly(Poly& poly, Soup const& s, bool finalize)
    {
        for (auto const& e : s.e)
        {
            int v0 = poly.InsertVertex(s.v[e[0]]);
            int v1 = poly.InsertVertex(s.v[e[1]]);
            poly.InsertEdge(Poly::Edge(v0, v1));
        }
        if (finalize) { poly.Finalize(); }
    }

    void OutTree(oracle::Ctx& io, Poly::BSPTree2 const* node)
    {
        io.outInt(node->mCoincident.size());
        for (auto const& e : node->mCoincident) { io.outInt(e.V[0]); io.outInt(e.V[1]); }
        io.outBool(node->mPosChild != nullptr);
        if (node->mPosChild) { OutTree(io, node->mPosChild.get()); }
        io.outBool(node->mNegChild != nullptr);
        if (node->mNegChild) { OutTree(io, node->mNegChild.get()); }
    }

    // The whole state: vertices, the edge-map size and the edge array (which
    // can be longer, finding #169), the BSP tree.
    void OutPoly(oracle::Ctx& io, Poly const& poly)
    {
        io.outInt(poly.GetNumVertices());
        for (int i = 0; i < poly.GetNumVertices(); ++i)
        {
            io.outReal(poly.GetVertex(i)[0]);
            io.outReal(poly.GetVertex(i)[1]);
        }
        io.outInt(poly.GetNumEdges());
        io.outInt(poly.mEArray.size());
        for (auto const& e : poly.mEArray) { io.outInt(e.V[0]); io.outInt(e.V[1]); }
        io.outBool(poly.mTree != nullptr);
        if (poly.mTree) { OutTree(io, poly.mTree.get()); }
    }
}

namespace
{
    // The epsilon of a polygon: 0 (most records), 1e-10 or 0.125.
    double DrawEpsilon(oracle::Ctx& io)
    {
        int k = io.rawInteger(0, 5);
        return io.given(k <= 3 ? 0.0 : (k == 4 ? 1e-10 : 0.125));
    }

    // The first polygon: a lattice shape, sometimes with perturbed
    // coordinates.
    Soup DrawSoupP(oracle::Ctx& io)
    {
        int shape = io.rawInteger(0, 6);
        bool uniform = (io.rawInteger(0, 4) == 0);
        return DrawShape(io, shape, uniform);
    }

    // The second polygon, related to the first:
    //   0 independent shape   1 identical   2 translated by a lattice
    //   offset (shared and overlapping edges, shared vertices)
    //   3 reversed orientation   4 scaled about a vertex (nested)
    Soup DrawSoupQ(oracle::Ctx& io, Soup const& p)
    {
        int rel = io.rawInteger(0, 4);
        if (rel == 0) { return DrawSoupP(io); }
        Soup q = p;
        if (rel == 2)
        {
            double dx = io.rawInteger(-4, 4), dy = io.rawInteger(-4, 4);
            for (auto& v : q.v) { v[0] += dx; v[1] += dy; }
        }
        else if (rel == 3)
        {
            for (auto& e : q.e) { std::swap(e[0], e[1]); }
        }
        else if (rel == 4)
        {
            Vector2<double> c = p.v[0];
            double f = (io.rawInteger(0, 1) != 0 ? 0.5 : 2.0);
            for (auto& v : q.v) { v = c + f * (v - c); }
        }
        return q;
    }

    // Query points (recorded): the polygon's vertices, edge midpoints,
    // lattice points and uniform points.
    std::vector<Vector2<double>> DrawQueries(oracle::Ctx& io, Soup const& p, Soup const& q)
    {
        std::vector<Vector2<double>> pts;
        for (int i = 0; i < 8; ++i)
        {
            Vector2<double> x;
            switch (io.rawInteger(0, 4))
            {
            case 0: x = p.v[io.rawInteger(0, static_cast<int>(p.v.size()) - 1)]; break;
            case 1: x = q.v[io.rawInteger(0, static_cast<int>(q.v.size()) - 1)]; break;
            case 2:
            {
                auto const& e = p.e[io.rawInteger(0, static_cast<int>(p.e.size()) - 1)];
                x = 0.5 * (p.v[e[0]] + p.v[e[1]]);
                break;
            }
            case 3: x = { io.rawInteger(-8, 8) * 0.5, io.rawInteger(-8, 8) * 0.5 }; break;
            default: x = { io.raw(-8.0, 8.0), io.raw(-8.0, 8.0) }; break;
            }
            pts.push_back(x);
            io.given(x[0]);
            io.given(x[1]);
        }
        return pts;
    }

    // One Boolean operation (0 ~P, 1 P & Q, 2 P | Q, 3 P - Q, 4 P ^ Q), the
    // whole state of the result and its point locations. An empty result
    // throws (finding #169, preserved): Finalize asserts a nonempty edge
    // list. Nothing is emitted before the operation.
    void PolyBoolean(oracle::Ctx& io, int op)
    {
        double eps = DrawEpsilon(io);
        Soup sp = DrawSoupP(io);
        RecordSoup(io, sp);
        Soup sq = DrawSoupQ(io, sp);
        RecordSoup(io, sq);
        std::vector<Vector2<double>> pts = DrawQueries(io, sp, sq);
        Poly P(eps), Q(eps);
        BuildPoly(P, sp, true);
        BuildPoly(Q, sq, true);
        Poly R(eps);
        switch (op)
        {
        case 0: R = ~P; break;
        case 1: R = P & Q; break;
        case 2: R = P | Q; break;
        case 3: R = P - Q; break;
        default: R = P ^ Q; break;
        }
        OutPoly(io, R);
        for (auto const& x : pts) { io.outInt(R.PointLocation(x)); }
    }
}

ORACLE_CASE("BSPPolygon2.negation") { PolyBoolean(io, 0); }
ORACLE_CASE("BSPPolygon2.intersection") { PolyBoolean(io, 1); }
ORACLE_CASE("BSPPolygon2.union") { PolyBoolean(io, 2); }
ORACLE_CASE("BSPPolygon2.difference") { PolyBoolean(io, 3); }
ORACLE_CASE("BSPPolygon2.exclusiveOr") { PolyBoolean(io, 4); }

// Construction (InsertVertex returns existing indices for repeated
// vertices, InsertEdge for repeated edges), Finalize (the BSP tree splits
// edges, which inserts vertices and edges), the accessors, PointLocation,
// copy construction and assignment.
ORACLE_CASE("BSPPolygon2.construct")
{
    double eps = DrawEpsilon(io);
    Soup sp = DrawSoupP(io);
    if (io.rawInteger(0, 3) == 0)
    {
        Soup extra = DrawShape(io, io.rawInteger(0, 6), false);
        int base = static_cast<int>(sp.v.size());
        for (auto const& v : extra.v) { sp.v.push_back(v); }
        for (auto const& e : extra.e) { sp.e.push_back({ e[0] + base, e[1] + base }); }
    }
    RecordSoup(io, sp);
    std::vector<Vector2<double>> pts = DrawQueries(io, sp, sp);
    Poly P(eps);
    std::vector<int> indices;
    for (auto const& e : sp.e)
    {
        int v0 = P.InsertVertex(sp.v[e[0]]);
        int v1 = P.InsertVertex(sp.v[e[1]]);
        indices.push_back(v0);
        indices.push_back(v1);
        indices.push_back(P.InsertEdge(Poly::Edge(v0, v1)));
    }
    // Finalize can throw (SplitEdge's "Edge does not exist" after an edge
    // was left unmapped, finding #169), so everything is computed before
    // anything is emitted; the state before Finalize is kept in a copy.
    Poly pre(P);
    P.Finalize();
    std::vector<int> locations;
    for (auto const& x : pts) { locations.push_back(P.PointLocation(x)); }
    Poly C(P);
    Poly D(0.5);
    D = C;
    for (int i : indices) { io.outInt(i); }
    OutPoly(io, pre);
    OutPoly(io, P);
    for (int loc : locations) { io.outInt(loc); }
    OutPoly(io, D);
    for (int i = 0; i < D.GetNumEdges(); ++i)
    {
        io.outInt(D.GetEdge(i).V[0]);
        io.outInt(D.GetEdge(i).V[1]);
    }
}

// Finding #169 item 2 (preserved): SplitEdge(v0, v1, vmid) erases <v0,v1>,
// rewrites the edge-array entry to <v0,vmid> and std::map::insert's it,
// which is a silent no-op when <v0,vmid> is already an edge. The soup has a
// splitter edge whose line crosses the edge A->C exactly at the existing
// vertex B, and an edge A->B: after Finalize the edge array holds an entry
// the map does not index (GetNumEdges < the array size). The lattice
// placement is random (translation, scale, one of four orientations), plus
// a closing edge set and an optional second loop; the state, the point
// locations and a Boolean operation with a square are compared.
ORACLE_CASE("BSPPolygon2.splitEdge.unmapped")
{
    double eps = DrawEpsilon(io);
    Soup s;
    auto P = [](double x, double y) { return Vector2<double>{ x, y }; };
    // Splitter <S0,S1> (vertical through x = 2), then A->C, A->B, and edges
    // closing the figure.
    s.v = { P(2, -2), P(2, 2), P(0, 0), P(4, 0), P(2, 0), P(4, 3), P(0, 3) };
    s.e = { { 0, 1 }, { 2, 3 }, { 2, 4 }, { 3, 5 }, { 5, 6 }, { 6, 2 } };
    if (io.rawInteger(0, 1) != 0) { s.e.push_back({ 1, 0 }); }
    int rot = io.rawInteger(0, 3);
    double scale = io.rawInteger(1, 3);
    double tx = io.rawInteger(-3, 3), ty = io.rawInteger(-3, 3);
    for (auto& p : s.v)
    {
        Vector2<double> r = p;
        for (int k = 0; k < rot; ++k) { r = P(-r[1], r[0]); }
        p = P(scale * r[0] + tx, scale * r[1] + ty);
    }
    RecordSoup(io, s);
    Soup sq;
    double w = io.rawInteger(1, 5);
    AddLoop(sq, { P(tx, ty), P(tx + w, ty), P(tx + w, ty + w), P(tx, ty + w) });
    RecordSoup(io, sq);
    std::vector<Vector2<double>> pts = DrawQueries(io, s, sq);
    int op = io.integer(0, 3);
    Poly A(eps), B(eps);
    BuildPoly(A, s, true);
    BuildPoly(B, sq, true);
    OutPoly(io, A);
    for (auto const& x : pts) { io.outInt(A.PointLocation(x)); }
    // The overlapping collinear edges A->B and A->C often make the Boolean
    // operation assert (GetCoPartition's interval ordering, or an empty
    // result); the exception is caught and reported as a flag so that the
    // state above is compared on every record.
    Poly R(eps);
    bool threw = false;
    try
    {
        switch (op)
        {
        case 0: R = A & B; break;
        case 1: R = A | B; break;
        case 2: R = A - B; break;
        default: R = B - A; break;
        }
    }
    catch (std::exception const&)
    {
        threw = true;
    }
    io.outBool(threw);
    if (!threw)
    {
        OutPoly(io, R);
        for (auto const& x : pts) { io.outInt(R.PointLocation(x)); }
    }
}

// Throw parity: a degenerate edge (InsertEdge asserts), PointLocation and
// the Boolean operators before Finalize (LogAssert "Tree must exist"),
// Finalize of an empty polygon (BSPTree2 asserts a nonempty edge list), and
// the intersection of two disjoint polygons (empty result, finding #169).
ORACLE_CASE("BSPPolygon2.invalid")
{
    int kind = io.integer(0, 4);
    Soup sp = DrawShape(io, io.rawInteger(0, 6), false);
    RecordSoup(io, sp);
    Poly P(0.0);
    if (kind == 0)
    {
        BuildPoly(P, sp, false);
        int v = P.InsertVertex(sp.v[0]);
        P.InsertEdge(Poly::Edge(v, v));
        io.outInt(v);
    }
    else if (kind == 1)
    {
        BuildPoly(P, sp, false);
        io.outInt(P.PointLocation(sp.v[0]));
    }
    else if (kind == 2)
    {
        BuildPoly(P, sp, false);
        Poly N = ~P;
        io.outInt(N.GetNumEdges());
    }
    else if (kind == 3)
    {
        P.Finalize();
        io.outInt(P.GetNumEdges());
    }
    else
    {
        BuildPoly(P, sp, true);
        Soup far = sp;
        for (auto& v : far.v) { v[0] += 100.0; }
        Poly F(0.0);
        BuildPoly(F, far, true);
        Poly I = P & F;
        io.outInt(I.GetNumEdges());
    }
}

// ================================================================ BSplines

namespace
{
    // BasisSpec, DrawBasis, RecordBasis, ToInput, DrawParameter and
    // DrawPoints are copied from v05-ap.cpp. Knot values are dyadic so knot
    // differences and the domain are exact.
    //   0 open uniform (as BasisFunctionInput(n, d) builds it)
    //   1 open nonuniform, 2 open nonuniform with a repeated interior knot,
    //   3 periodic uniform, 4 periodic nonuniform
    struct BasisSpec
    {
        int32_t numControls = 0;
        int32_t degree = 0;
        bool uniform = false;
        bool periodic = false;
        std::vector<UniqueKnot<double>> knots{};
    };

    BasisSpec DrawBasis(oracle::Ctx& io, int32_t minDegree, int32_t maxDegree,
        int32_t maxExtra, int32_t maxMode = 4)
    {
        int32_t mode = io.rawInteger(0, maxMode);
        BasisSpec s{};
        s.degree = io.rawInteger(minDegree, maxDegree);
        s.periodic = (mode >= 3);
        s.numControls = s.degree + 1 + io.rawInteger(0, maxExtra);
        if (!s.periodic)
        {
            int32_t interiorTotal = s.numControls - s.degree - 1;
            std::vector<int32_t> mult{};
            if (mode == 2 && interiorTotal >= 2)
            {
                int32_t repeated = (s.degree >= 2 ? 2 : 1);
                mult.push_back(repeated);
                for (int32_t i = 0; i < interiorTotal - repeated; ++i) { mult.push_back(1); }
            }
            else
            {
                for (int32_t i = 0; i < interiorTotal; ++i) { mult.push_back(1); }
            }
            s.uniform = (mode == 0);
            if (mode == 0)
            {
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

    // A parameter on [tmin, tmax] (unrecorded): the endpoints, an interior
    // knot, dyadic interior points, just outside (clamping), far outside
    // (periodic wrapping) and uniform.
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
                size_t i = static_cast<size_t>(io.rawInteger(1, static_cast<int32_t>(knots->size()) - 2));
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

    // 'count' points of dimension 'dim' (flat, unrecorded). Modes: 0 uniform
    // [-5,5], 1 integer lattice [-3,3], 2 all points equal, 3 collinear
    // (p + t*dir with dyadic t), 4 all zero with signed zeros.
    std::vector<double> DrawPoints(oracle::Ctx& io, int32_t count, int32_t dim, int32_t mode)
    {
        std::vector<double> P(static_cast<size_t>(count) * dim);
        std::vector<double> base(dim), dir(dim);
        for (int32_t j = 0; j < dim; ++j)
        {
            base[j] = static_cast<double>(io.rawInteger(-3, 3));
            dir[j] = static_cast<double>(io.rawInteger(-2, 2));
        }
        for (int32_t i = 0; i < count; ++i)
        {
            double t = io.rawInteger(-8, 8) * 0.25;
            for (int32_t j = 0; j < dim; ++j)
            {
                double& v = P[static_cast<size_t>(i) * dim + j];
                switch (mode)
                {
                case 0: v = io.raw(-5.0, 5.0); break;
                case 1: v = static_cast<double>(io.rawInteger(-3, 3)); break;
                case 2: v = base[j]; break;
                case 3: v = base[j] + t * dir[j]; break;
                default: v = RawSign(io) * 0.0; break;
                }
            }
        }
        return P;
    }

    template <int N>
    std::vector<Vector<N, double>> DrawControlsN(oracle::Ctx& io, int32_t count, int32_t mode)
    {
        std::vector<double> P = DrawPoints(io, count, N, mode);
        std::vector<Vector<N, double>> C(static_cast<size_t>(count));
        for (int32_t i = 0; i < count; ++i)
        {
            for (int j = 0; j < N; ++j) { C[i][j] = io.given(P[static_cast<size_t>(i) * N + j]); }
        }
        return C;
    }
}

// =========================================================== BSplineCurve

namespace
{
    template <int N>
    Vector<N, double> DrawLatticeVector(oracle::Ctx& io)
    {
        Vector<N, double> v;
        for (int j = 0; j < N; ++j) { v[j] = io.lattice(-4, 4); }
        return v;
    }

    template <int N>
    void CurveBody(oracle::Ctx& io, BasisSpec const& s, int cmode)
    {
        int32_t n = s.numControls;
        std::vector<Vector<N, double>> C;
        if (cmode <= 4) { C = DrawControlsN<N>(io, n, cmode); }
        BasisFunctionInput<double> input = ToInput(s);
        BSplineCurve<N, double> curve(input, cmode <= 4 ? C.data() : nullptr);
        for (int k = 0; k < 3; ++k)
        {
            int i = io.integer(-1, n);
            Vector<N, double> value = DrawLatticeVector<N>(io);
            curve.SetControl(i, value);
        }
        BasisFunction<double> const& basis = curve.GetBasisFunction();
        io.outBool(static_cast<bool>(curve));
        io.outReal(curve.GetTMin());
        io.outReal(curve.GetTMax());
        io.outInt(curve.GetNumSegments());
        io.outInt(curve.GetNumControls());
        io.outInt(basis.GetNumKnots());
        io.outInt(basis.GetDegree());
        io.outBool(basis.IsOpen());
        io.outBool(basis.IsUniform());
        io.outBool(basis.IsPeriodic());
        io.outReal(basis.GetMinDomain());
        io.outReal(basis.GetMaxDomain());
        for (int k = 0; k < 3; ++k)
        {
            int i = io.integer(-1, n);
            io.outVec(curve.GetControl(i));
        }
        io.outVec(curve.GetControls()[n - 1]);
        static double const orders[6] = { 0, 1, 2, 3, 4, 5 };
        for (int k = 0; k < 5; ++k)
        {
            double t = io.given(DrawParameter(io, curve.GetTMin(), curve.GetTMax(), &s.knots));
            int order = static_cast<int>(io.given(orders[io.rawInteger(0, 5)]));
            std::array<Vector<N, double>, 4> jet{};
            curve.Evaluate(t, static_cast<uint32_t>(order), jet.data());
            for (int i = 0; i < 4; ++i) { io.outVec(jet[i]); }
        }
        double t = io.given(DrawParameter(io, curve.GetTMin(), curve.GetTMax(), &s.knots));
        io.outVec(curve.GetPosition(t));
        io.outVec(curve.GetTangent(t));
        io.outReal(curve.GetSpeed(t));
        double tl = io.given(DrawParameter(io, curve.GetTMin(), curve.GetTMax(), &s.knots));
        io.outReal(curve.GetLength(curve.GetTMin(), tl));
        double total = curve.GetTotalLength();
        io.outReal(total);
        double fraction = io.given(io.rawInteger(0, 8) * 0.125);
        io.outReal(curve.GetTime(fraction * total));
    }
}

// BSplineCurve<N, double> for N = 1..4 over open uniform, open nonuniform
// (with a repeated interior knot) and periodic (uniform and nonuniform)
// bases of degree 1..5; controls present (uniform, lattice, equal,
// collinear, zero) or deferred (null pointer: zero controls); SetControl and
// GetControl in and out of range; the BasisFunction and ParametricCurve
// accessors; Evaluate for orders 0..3 and >= SUP_ORDER (a zero jet) at the
// domain ends, knots, dyadic, outside (clamped or wrapped) and uniform
// parameters; GetPosition, GetTangent, GetSpeed, GetLength, GetTotalLength
// and GetTime (Romberg integration and the Newton/bisection inversion).
ORACLE_CASE("BSplineCurve.evaluate")
{
    int N = io.integer(1, 4);
    BasisSpec s = DrawBasis(io, 1, 5, 5);
    RecordBasis(io, s);
    int cmode = io.integer(0, 5);
    if (N == 1) { CurveBody<1>(io, s, cmode); }
    else if (N == 2) { CurveBody<2>(io, s, cmode); }
    else if (N == 3) { CurveBody<3>(io, s, cmode); }
    else { CurveBody<4>(io, s, cmode); }
}

// ======================================================= BSplineReduction

namespace
{
    template <int N>
    void ReductionBody(oracle::Ctx& io, int numIn, int degree, double fraction, int cmode)
    {
        std::vector<Vector<N, double>> in = DrawControlsN<N>(io, numIn, cmode);
        BSplineReduction<N, double> reduction;
        std::vector<Vector<N, double>> out;
        reduction(in, degree, fraction, out);
        io.outInt(out.size());
        for (auto const& c : out) { io.outVec(c); }
        // Reuse the object (its knot arrays are resized) with the reversed
        // controls and fraction 1/2.
        std::reverse(in.begin(), in.end());
        std::vector<Vector<N, double>> out2;
        reduction(in, degree, 0.5, out2);
        io.outInt(out2.size());
        for (auto const& c : out2) { io.outVec(c); }
    }
}

// operator()(inControls, degree, fraction, outControls) for N = 1..3,
// 2..12 input controls, degree 1..min(numIn - 1, 4), and fractions: uniform
// in [0, 1], exact multiples k/numIn (the truncation boundary), 1 and above
// (a copy), 0 and negative (clamped to degree + 1 controls). The Gram
// matrix entries come from Integration::Romberg(8) and the inverse from
// BandedMatrix::ComputeInverse (finding #169: the Romberg error across
// knots is preserved; compared bit for bit).
ORACLE_CASE("BSplineReduction.reduce")
{
    int N = io.integer(1, 3);
    int numIn = io.integer(2, 12);
    int maxDegree = (numIn - 1 < 4 ? numIn - 1 : 4);
    int degree = io.integer(1, maxDegree);
    double f;
    switch (io.rawInteger(0, 5))
    {
    case 0: case 1: f = io.raw(0.0, 1.0); break;
    case 2: f = static_cast<double>(io.rawInteger(0, numIn)) / numIn; break;
    case 3: f = (io.rawInteger(0, 1) != 0 ? 1.0 : 1.5); break;
    case 4: f = (io.rawInteger(0, 1) != 0 ? 0.0 : -0.5); break;
    default: f = std::nextafter(static_cast<double>(io.rawInteger(1, numIn)) / numIn, 0.0); break;
    }
    double fraction = io.given(f);
    int cmode = io.integer(0, 4);
    if (N == 1) { ReductionBody<1>(io, numIn, degree, fraction, cmode); }
    else if (N == 2) { ReductionBody<2>(io, numIn, degree, fraction, cmode); }
    else { ReductionBody<3>(io, numIn, degree, fraction, cmode); }
}

// The LogAssert on invalid input: fewer than 2 controls, degree < 1,
// degree >= number of controls.
ORACLE_CASE("BSplineReduction.invalid")
{
    int kind = io.integer(0, 2);
    int numIn = (kind == 0 ? io.integer(0, 1) : io.integer(2, 6));
    int degree = (kind == 1 ? io.integer(-1, 0) : io.integer(numIn, numIn + 2));
    std::vector<Vector<2, double>> in = DrawControlsN<2>(io, numIn, 1);
    BSplineReduction<2, double> reduction;
    std::vector<Vector<2, double>> out;
    reduction(in, degree, 0.5, out);
    io.outInt(out.size());
}

// ======================================================== BSplineGeodesic

// BSplineGeodesic's metric and Christoffel symbols come from the B-spline
// surface (arithmetic only); RiemannianGeodesic adds std::sqrt and the
// Gaussian-elimination Inverse of GMatrix. Unlike EllipsoidGeodesic (v18) no
// C math library call decides the control flow, so the whole algorithm,
// including the argmin of Refine's line search, is compared bit for bit.
// The tuning parameters are drawn small (the derived steps stay those of
// the defaults, as upstream computes them only in the constructor).

namespace
{
    // A graph-like 3D B-spline surface over two bases (recorded): control
    // (i0, i1) is (i0 * sx + p, i1 * sy + q, z) with |p|, |q| < sx/4, so the
    // partial derivatives are independent and the metric is positive
    // definite on the domain; z is uniform or on a lattice.
    struct GeoSurface
    {
        BasisSpec s0, s1;
        std::vector<Vector3<double>> C;
    };

    GeoSurface DrawGeoSurface(oracle::Ctx& io, int maxMode)
    {
        GeoSurface g;
        g.s0 = DrawBasis(io, 1, 3, 2, maxMode);
        RecordBasis(io, g.s0);
        g.s1 = DrawBasis(io, 1, 3, 2, maxMode);
        RecordBasis(io, g.s1);
        int n0 = g.s0.numControls, n1 = g.s1.numControls;
        double sx = io.rawInteger(1, 3) * 0.5, sy = io.rawInteger(1, 3) * 0.5;
        bool lattice = (io.rawInteger(0, 1) != 0);
        g.C.resize(static_cast<size_t>(n0) * n1);
        for (int i1 = 0; i1 < n1; ++i1)
        {
            for (int i0 = 0; i0 < n0; ++i0)
            {
                Vector3<double>& c = g.C[i0 + static_cast<size_t>(n0) * i1];
                c[0] = io.given(i0 * sx + (lattice ? 0.0 : io.raw(-0.2, 0.2) * sx));
                c[1] = io.given(i1 * sy + (lattice ? 0.0 : io.raw(-0.2, 0.2) * sy));
                c[2] = io.given(lattice ? static_cast<double>(io.rawInteger(-2, 2)) : io.raw(-1.0, 1.0));
            }
        }
        return g;
    }

    GVector<double> DrawDomainPoint(oracle::Ctx& io, BSplineSurface<3, double> const& surface,
        GeoSurface const& g)
    {
        GVector<double> p(2);
        p[0] = io.given(DrawParameter(io, surface.GetUMin(), surface.GetUMax(), &g.s0.knots));
        p[1] = io.given(DrawParameter(io, surface.GetVMin(), surface.GetVMax(), &g.s1.knots));
        return p;
    }

    void DrawTuning(oracle::Ctx& io, BSplineGeodesic<double>& geo)
    {
        static double const steps[3] = { 1e-4, 1e-3, 1e-2 };
        static double const radii[3] = { 0.5, 1.0, 2.0 };
        geo.integralSamples = io.integer(2, 8);
        geo.searchSamples = io.integer(1, 6);
        geo.derivativeStep = io.given(steps[io.rawInteger(0, 2)]);
        geo.subdivisions = io.integer(0, 3);
        geo.refinements = io.integer(0, 2);
        geo.searchRadius = io.given(radii[io.rawInteger(0, 2)]);
    }

    void OutGVec2(oracle::Ctx& io, GVector<double> const& v)
    {
        io.outInt(v.GetSize());
        for (int i = 0; i < v.GetSize(); ++i) { io.outReal(v[i]); }
    }
}

// ComputeGeodesic(end0, end1, quantity, path) with the refine callback
// (its calls and the progress accessors it can read), then ComputeTotalLength
// and ComputeTotalCurvature of the path. Open bases only (the metric stays
// positive definite).
ORACLE_CASE("BSplineGeodesic.computeGeodesic")
{
    GeoSurface g = DrawGeoSurface(io, 2);
    BasisFunctionInput<double> inputs[2] = { ToInput(g.s0), ToInput(g.s1) };
    BSplineSurface<3, double> surface(inputs, g.C.data());
    BSplineGeodesic<double> geo(surface);
    DrawTuning(io, geo);
    GVector<double> end0 = DrawDomainPoint(io, surface, g);
    GVector<double> end1 = DrawDomainPoint(io, surface, g);
    int calls = 0;
    int64_t progress = 0;
    geo.refineCallback = [&]()
    {
        ++calls;
        progress = progress * 7 + geo.GetSubdivisionStep() * 100 + geo.GetRefinementStep() * 10
            + geo.GetCurrentQuantity();
        progress %= 1000000007;
    };
    int32_t quantity = 0;
    std::vector<GVector<double>> path;
    geo.ComputeGeodesic(end0, end1, quantity, path);
    double length = geo.ComputeTotalLength(quantity, path);
    double curvature = geo.ComputeTotalCurvature(quantity, path);
    io.outInt(geo.GetDimension());
    io.outInt(quantity);
    io.outInt(path.size());
    for (auto const& p : path) { OutGVec2(io, p); }
    io.outInt(calls);
    io.outInt(progress);
    io.outInt(geo.GetSubdivisionStep());
    io.outInt(geo.GetRefinementStep());
    io.outInt(geo.GetCurrentQuantity());
    io.outReal(length);
    io.outReal(curvature);
}

// ComputeSegmentLength, ComputeSegmentCurvature, Subdivide and Refine on
// one segment (the callback is suppressed in Subdivide and called once in
// Refine). Periodic bases included.
ORACLE_CASE("BSplineGeodesic.segment")
{
    GeoSurface g = DrawGeoSurface(io, 4);
    BasisFunctionInput<double> inputs[2] = { ToInput(g.s0), ToInput(g.s1) };
    BSplineSurface<3, double> surface(inputs, g.C.data());
    BSplineGeodesic<double> geo(surface);
    DrawTuning(io, geo);
    GVector<double> p0 = DrawDomainPoint(io, surface, g);
    GVector<double> p1 = DrawDomainPoint(io, surface, g);
    GVector<double> m = DrawDomainPoint(io, surface, g);
    int calls = 0;
    geo.refineCallback = [&calls]() { ++calls; };
    double length = geo.ComputeSegmentLength(p0, p1);
    double curvature = geo.ComputeSegmentCurvature(p0, p1);
    GVector<double> mid(2);
    bool changed0 = geo.Subdivide(p0, mid, p1);
    bool changed1 = geo.Refine(p0, m, p1);
    io.outReal(length);
    io.outReal(curvature);
    io.outBool(changed0);
    OutGVec2(io, mid);
    io.outBool(changed1);
    OutGVec2(io, m);
    io.outInt(calls);
}

// Throw parity: a segment of zero length (the quadratic form is 0 and
// LogAssert(qForm > 0) fires), and a surface with deferred (zero) controls.
ORACLE_CASE("BSplineGeodesic.degenerate")
{
    GeoSurface g = DrawGeoSurface(io, 2);
    int kind = io.integer(0, 1);
    BasisFunctionInput<double> inputs[2] = { ToInput(g.s0), ToInput(g.s1) };
    BSplineSurface<3, double> surface(inputs, kind == 0 ? g.C.data() : nullptr);
    BSplineGeodesic<double> geo(surface);
    GVector<double> p0 = DrawDomainPoint(io, surface, g);
    GVector<double> p1 = (kind == 0 ? p0 : DrawDomainPoint(io, surface, g));
    if (kind == 1 && p1 == p0) { p1[0] += 0.125; }
    io.outReal(geo.ComputeSegmentLength(p0, p1));
}
