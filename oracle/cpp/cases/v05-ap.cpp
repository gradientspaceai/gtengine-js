// Group 5 (arbitrary precision; the group is alphabetical): BSNumber,
// BSPrecision, QFNumber, SWInterval, BSplineCurveFit, BSplineSurface,
// BSplineSurfaceFit, BSplineVolume.
//
// BSNumber is instantiated with UIntegerAP32, the arbitrary-precision
// storage. The port replaces the UInteger layer with bigint, so a BSNumber is
// compared by VALUE and canonical form, never by storage layout: OutBSN emits
// the sign, the biased exponent, the exponent, the bit count of the odd
// integer, its 32-bit words (all of them up to 8 words; beyond that the top
// three, the bottom two and a 48-bit FNV-1a digest of all words) and both
// floating-point conversions, operator double and operator float.
//
// Everything is + - * / sqrt, comparisons and integer/bit manipulation except
// the std:: math overloads of BSNumber that round-trip through a libm call;
// those live in the one tolerance case BSNumber.std.libm.
#define ORACLE_FAMILY "v05-ap"
#include "Oracle.h"

#include <Mathematics/BSNumber.h>
#include <Mathematics/UIntegerAP32.h>
#include <Mathematics/BSPrecision.h>
#include <Mathematics/QFNumber.h>
#include <Mathematics/SWInterval.h>
#include <Mathematics/BSplineCurveFit.h>
#include <Mathematics/BSplineSurface.h>
#include <Mathematics/BSplineSurfaceFit.h>
#include <Mathematics/BSplineVolume.h>

#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

using namespace gte;
using BSN = BSNumber<UIntegerAP32>;

namespace
{
    // ------------------------------------------------------------ doubles

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

    // A finite double over the whole binary64 range, drawn (not recorded).
    // The modes visit the values the conversions branch on: small integers,
    // exact powers of two from 2^-1074 to 2^1023, odd and all-ones
    // mantissas, arbitrary bit patterns, subnormals, signed zeros, the top
    // and bottom normal binades and binary32 values including binary32
    // subnormals and values just below and above FLT_MAX.
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
            // (2^m - 1) * 2^k: rounding the all-ones mantissa carries into
            // the next power of two.
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
    // ----------------------------------------------------------- BSNumber

    // A BSNumber built from recorded doubles. The mode is recorded first; the
    // replay reads the same number of doubles and applies the same
    // operations in the same order.
    //   0  BSN(a)                    1 double
    //   1  BSN(a) + BSN(b)           2 doubles, any exponent gap (<= 2100 bits)
    //   2  BSN(a) * BSN(b)           2 doubles (<= 106 bits)
    //   3  BSN(a) * BSN(b) + BSN(c)  3 doubles
    //   4  BSN(a) + BSN(b)           moderate gap, 2..13 words
    //   5  a + h + c                 h = +-half ulp of a (binary64 tie), c = 0
    //                                or a tiny tie breaker
    //   6  f + h + c                 the same for a binary32 value f
    //   7  BSN(a) * BSN(b)           a subnormal (binary64 or binary32 grid),
    //                                b in {1/8, 1/4, 3/8, 1/2, 3/4, 5/4, 3/2}:
    //                                values between and halfway between the
    //                                smallest subnormals
    //   8  M + h + c                 M = DBL_MAX or FLT_MAX, h in {+-half ulp,
    //                                +ulp}: the overflow boundary
    BSN DrawBSN(oracle::Ctx& io, int lo, int hi)
    {
        int mode = io.integer(lo, hi);
        switch (mode)
        {
        case 0:
        {
            double a = io.given(DrawDouble(io));
            return BSN(a);
        }
        case 1:
        case 2:
        {
            double a = io.given(DrawDouble(io));
            double b = io.given(DrawDouble(io));
            return (mode == 1 ? BSN(a) + BSN(b) : BSN(a) * BSN(b));
        }
        case 3:
        {
            double a = io.given(DrawDouble(io));
            double b = io.given(DrawDouble(io));
            double c = io.given(DrawDouble(io));
            return BSN(a) * BSN(b) + BSN(c);
        }
        case 4:
        {
            double a = io.given(DrawNormal(io, -60, 60));
            double b = io.given(DrawNormal(io, -360, 60));
            return BSN(a) + BSN(b);
        }
        case 5:
        {
            double a = io.given(DrawNormal(io, -1021, 1022));
            int e = std::ilogb(a) - 53;
            double hs = RawSign(io);
            double h = io.given(hs * Pow2(e));
            int j = io.rawInteger(0, 60);
            double cs = RawSign(io);
            double c = io.given(j == 0 || e - j < -1074 ? 0.0 : cs * Pow2(e - j));
            return BSN(a) + BSN(h) + BSN(c);
        }
        case 6:
        {
            double fv;
            if (io.rawInteger(0, 3) == 0)
            {
                double m = static_cast<double>(io.rawInteger(1, (1 << 23) - 1));
                fv = RawSign(io) * std::ldexp(m, -149);
            }
            else
            {
                fv = static_cast<double>(static_cast<float>(DrawNormal(io, -126, 127)));
                if (!std::isfinite(fv)) { fv = static_cast<double>(std::numeric_limits<float>::max()); }
            }
            double f = io.given(fv);
            int ef = std::ilogb(f);
            int e = (ef >= -126 ? ef : -126) - 24;
            double hs = RawSign(io);
            double h = io.given(hs * Pow2(e));
            int j = io.rawInteger(0, 60);
            double cs = RawSign(io);
            double c = io.given(j == 0 ? 0.0 : cs * Pow2(e - j));
            return BSN(f) + BSN(h) + BSN(c);
        }
        case 7:
        {
            static double const factors[7] = { 0.125, 0.25, 0.375, 0.5, 0.75, 1.25, 1.5 };
            int grid = (io.rawInteger(0, 1) != 0 ? -1074 : -149);
            double m = static_cast<double>(io.rawInteger(1, 8));
            double a = io.given(RawSign(io) * std::ldexp(m, grid));
            double b = io.given(factors[io.rawInteger(0, 6)]);
            return BSN(a) * BSN(b);
        }
        default:
        {
            bool dbl = (io.rawInteger(0, 1) != 0);
            double top = (dbl ? std::numeric_limits<double>::max()
                : static_cast<double>(std::numeric_limits<float>::max()));
            int e = (dbl ? 970 : 103);
            double M = io.given(RawSign(io) * top);
            int hk = io.rawInteger(0, 2);
            double hm = (hk == 0 ? 1.0 : (hk == 1 ? -1.0 : 2.0));
            double h = io.given((M < 0.0 ? -hm : hm) * Pow2(e));
            int j = io.rawInteger(0, 60);
            double cs = RawSign(io);
            double c = io.given(j == 0 ? 0.0 : cs * Pow2(e - j));
            return BSN(M) + BSN(h) + BSN(c);
        }
        }
    }

    // Emit a BSNumber by value and canonical form (see the file comment).
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
}

// =============================================================== BSNumber

// BSNumber(double): every finite double plus the infinities and NaNs, which
// upstream maps to +-2^1024 and 0 (the graceful exit when
// GTE_THROW_ON_CONVERT_FROM_INFINITY_OR_NAN is not defined).
ORACLE_CASE("BSNumber.construct.double")
{
    int kind = io.rawInteger(0, 11);
    double v = DrawDouble(io);
    if (kind == 0) { v = RawSign(io) * std::numeric_limits<double>::infinity(); }
    if (kind == 1) { v = std::numeric_limits<double>::quiet_NaN(); }
    double d = io.given(v);
    BSN x(d);
    OutBSN(io, x);
}

// BSNumber(float): binary32 values including subnormals, FLT_MIN, FLT_MAX,
// the infinities and NaN. The float is recorded exactly as a double.
ORACLE_CASE("BSNumber.construct.float")
{
    float f = 0.0f;
    int kind = io.rawInteger(0, 7);
    float s = static_cast<float>(RawSign(io));
    switch (kind)
    {
    case 0: f = s * std::numeric_limits<float>::infinity(); break;
    case 1: f = std::numeric_limits<float>::quiet_NaN(); break;
    case 2: f = s * std::ldexp(static_cast<float>(io.rawInteger(1, (1 << 23) - 1)), -149); break;
    case 3: f = s * std::numeric_limits<float>::max(); break;
    case 4: f = s * std::numeric_limits<float>::min(); break;
    case 5: f = s * 0.0f; break;
    default: f = static_cast<float>(DrawDouble(io)); break;
    }
    double fd = io.given(static_cast<double>(f));
    BSN x(static_cast<float>(fd));
    OutBSN(io, x);
}

// BSNumber(int32_t), (uint32_t), (int64_t), (uint64_t). The value is recorded
// as the high and low 32-bit words of its 64-bit two's complement form.
// INT32_MIN and INT64_MIN are excluded: the constructors evaluate '-number',
// which is signed overflow (undefined behaviour) for them.
ORACLE_CASE("BSNumber.construct.integer")
{
    int kind = io.integer(0, 3);
    uint64_t bits = 0;
    int mode = io.rawInteger(0, 5);
    switch (mode)
    {
    case 0: bits = static_cast<uint64_t>(static_cast<int64_t>(io.rawInteger(-20, 20))); break;
    case 1: bits = 1ull << io.rawInteger(0, 63); break;
    case 2:
    {
        uint64_t p = 1ull << io.rawInteger(1, 63);
        bits = (io.rawInteger(0, 1) != 0 ? p + 1ull : p - 1ull);
        break;
    }
    case 3:
    {
        static uint64_t const extremes[8] = {
            0x7FFFFFFFull, 0x80000001ull, 0xFFFFFFFFull, 0xFFFFFFFF80000001ull,
            0x7FFFFFFFFFFFFFFFull, 0x8000000000000001ull, 0xFFFFFFFFFFFFFFFFull,
            (1ull << 53) + 1ull };
        bits = extremes[io.rawInteger(0, 7)];
        break;
    }
    default: bits = RawBits64(io); break;
    }
    if (kind == 0)
    {
        int32_t v = static_cast<int32_t>(static_cast<uint32_t>(bits));
        if (v == std::numeric_limits<int32_t>::min()) { v = v + 1; }
        bits = static_cast<uint64_t>(static_cast<int64_t>(v));
    }
    else if (kind == 1)
    {
        bits &= 0xFFFFFFFFull;
    }
    else if (kind == 2 && bits == 0x8000000000000000ull)
    {
        bits = 0x8000000000000001ull;
    }
    io.given(static_cast<double>(bits >> 32));
    io.given(static_cast<double>(bits & 0xFFFFFFFFull));
    BSN x;
    if (kind == 0) { x = BSN(static_cast<int32_t>(static_cast<uint32_t>(bits))); }
    else if (kind == 1) { x = BSN(static_cast<uint32_t>(bits)); }
    else if (kind == 2) { x = BSN(static_cast<int64_t>(bits)); }
    else { x = BSN(static_cast<uint64_t>(bits)); }
    OutBSN(io, x);
}

namespace
{
    // A string for the BSNumber(std::string) constructor, recorded as its
    // length followed by the character codes.
    //   0  valid: optional sign, 1..45 digits with a nonzero leading digit
    //   1  "0", "+0", "-0"
    //   2  leading zero in a multi-digit number (asserts)
    //   3  a non-digit inside a multi-character number (asserts)
    //   4  "", "+", "-" (asserts)
    //   5  decimal point or exponent notation (asserts; BSRational's syntax)
    //   6  long valid numbers, 46..80 digits, near 2^64, 2^128 and 10^k
    // A single non-digit character (finding #95) has its own deviation case.
    std::string DrawNumberString(oracle::Ctx& io)
    {
        static char const* const junk = "abcxyzXYZ./:;= _eE,";
        int mode = io.rawInteger(0, 6);
        int signKind = io.rawInteger(0, 2);
        std::string sign = (signKind == 0 ? "" : (signKind == 1 ? "+" : "-"));
        std::string digits;
        switch (mode)
        {
        case 0:
        case 6:
        {
            int n = (mode == 0 ? io.rawInteger(1, 45) : io.rawInteger(46, 80));
            digits.push_back(static_cast<char>('1' + io.rawInteger(0, 8)));
            for (int i = 1; i < n; ++i) { digits.push_back(static_cast<char>('0' + io.rawInteger(0, 9))); }
            if (mode == 6 && io.rawInteger(0, 2) == 0)
            {
                static char const* const specials[4] = {
                    "18446744073709551616", "18446744073709551615",
                    "340282366920938463463374607431768211456",
                    "100000000000000000000000000000000000000000000000000000000000" };
                digits = specials[io.rawInteger(0, 3)];
            }
            break;
        }
        case 1:
            digits = "0";
            break;
        case 2:
        {
            int n = io.rawInteger(2, 12);
            digits = "0";
            for (int i = 1; i < n; ++i) { digits.push_back(static_cast<char>('0' + io.rawInteger(0, 9))); }
            break;
        }
        case 3:
        {
            int n = io.rawInteger(2, 12);
            digits.push_back(static_cast<char>('1' + io.rawInteger(0, 8)));
            for (int i = 1; i < n; ++i) { digits.push_back(static_cast<char>('0' + io.rawInteger(0, 9))); }
            int pos = io.rawInteger(1, n - 1);
            digits[pos] = junk[io.rawInteger(0, 18)];
            break;
        }
        case 4:
            digits = "";
            break;
        default:
        {
            static char const* const forms[6] = { "12.5", "1e5", "3.0", "0.25", "7E2", "1.5e-3" };
            digits = forms[io.rawInteger(0, 5)];
            break;
        }
        }
        return sign + digits;
    }

    void RecordString(oracle::Ctx& io, std::string const& s)
    {
        io.given(static_cast<double>(s.size()));
        for (char c : s) { io.given(static_cast<double>(static_cast<unsigned char>(c))); }
    }
}

ORACLE_CASE("BSNumber.construct.string")
{
    std::string s = DrawNumberString(io);
    RecordString(io, s);
    BSN x(s);
    OutBSN(io, x);
}

// Finding #95 (port fixed): ConvertToInteger validates the characters only
// when the string has more than one character, so upstream accepts a
// one-character non-digit ('x' becomes 72) while the port asserts. Every
// record deviates structurally (C++ returns, the port throws).
ORACLE_CASE("BSNumber.construct.string.singleChar")
{
    static char const* const junk = "abcxyzXYZ./:;= _eE,!#$%&()*<>?@[]^`{|}~";
    int signKind = io.rawInteger(0, 2);
    std::string s = (signKind == 0 ? "" : (signKind == 1 ? "+" : "-"));
    s.push_back(junk[io.rawInteger(0, 38)]);
    RecordString(io, s);
    BSN x(s);
    OutBSN(io, x);
}

namespace
{
    // The second operand of a comparison or an arithmetic operation, drawn
    // relative to the first so that ties and the left-aligned word
    // comparison of UIntegerALU32::operator< are reached.
    //   0  independent
    //   1  y = x
    //   2  y = -x
    //   3  y = x + e, e = +-2^k with k at most 120 bits below x's exponent
    //   4  y = x * 2^k, the same odd integer at another exponent
    //   5  x = a*c and y = b*c where a and b share their binary exponent but
    //      have odd integers of different lengths (x is replaced)
    void DrawPair(oracle::Ctx& io, BSN& x, BSN& y)
    {
        x = DrawBSN(io, 0, 8);
        int mode = io.integer(0, 5);
        switch (mode)
        {
        case 0:
            y = DrawBSN(io, 0, 8);
            break;
        case 1:
            y = x;
            break;
        case 2:
            y = -x;
            break;
        case 3:
        {
            int k = (x.GetSign() != 0 ? x.GetExponent() - io.rawInteger(0, 120) : io.rawInteger(-1074, 1023));
            k = (k < -1074 ? -1074 : (k > 1023 ? 1023 : k));
            double e = io.given(RawSign(io) * Pow2(k));
            y = x + BSN(e);
            break;
        }
        case 4:
        {
            double p = io.given(Pow2(io.rawInteger(-3, 3)));
            y = x * BSN(p);
            break;
        }
        default:
        {
            int E = io.rawInteger(-40, 40);
            int j0 = io.rawInteger(1, 30);
            int j1 = io.rawInteger(1, 30);
            double m0 = static_cast<double>(io.rawInteger(0, (1 << j0) - 1));
            double m1 = static_cast<double>(io.rawInteger(0, (1 << j1) - 1));
            double s = RawSign(io);
            double a = io.given(s * std::ldexp(1.0 + std::ldexp(m0, -j0), E));
            double b = io.given(s * std::ldexp(1.0 + std::ldexp(m1, -j1), E));
            double c = io.given(DrawNormal(io, -30, 30));
            x = BSN(a) * BSN(c);
            y = BSN(b) * BSN(c);
            break;
        }
        }
    }
}

// operator==, !=, <, <=, >, >= in both orders.
ORACLE_CASE("BSNumber.compare")
{
    BSN x, y;
    DrawPair(io, x, y);
    io.outBool(x == y);
    io.outBool(x != y);
    io.outBool(x < y);
    io.outBool(x <= y);
    io.outBool(x > y);
    io.outBool(x >= y);
    io.outBool(y < x);
    io.outBool(y <= x);
}

// Unary + and -, binary + - *, the compound assignments and Negate().
ORACLE_CASE("BSNumber.arithmetic")
{
    BSN x, y;
    DrawPair(io, x, y);
    OutBSN(io, +x);
    OutBSN(io, -x);
    OutBSN(io, x + y);
    OutBSN(io, x - y);
    OutBSN(io, y - x);
    OutBSN(io, x * y);
    BSN z = x;
    z += y;
    OutBSN(io, z);
    z -= x;
    OutBSN(io, z);
    z *= y;
    OutBSN(io, z);
    BSN w = x;
    w.Negate();
    OutBSN(io, w);
}

// operator double and operator float at their rounding boundaries: binary64
// and binary32 ties (both parities, with and without a tie breaker), values
// between the smallest subnormals, and the overflow boundary.
ORACLE_CASE("BSNumber.conversions.boundary")
{
    BSN x = DrawBSN(io, 5, 8);
    OutBSN(io, x);
    OutBSN(io, -x);
}

// SetSign, SetBiasedExponent, SetExponent and the matching getters. The sign
// is only flipped on a nonzero number: SetSign(+-1) on zero creates an
// invalid number whose conversion reads UInteger words that do not exist.
ORACLE_CASE("BSNumber.accessors")
{
    BSN x = DrawBSN(io, 0, 8);
    int e = io.integer(-1200, 1200);
    int be = io.integer(-1200, 1200);
    BSN y = x;
    y.SetExponent(e);
    OutBSN(io, y);
    BSN z = x;
    z.SetBiasedExponent(be);
    OutBSN(io, z);
    BSN w = x;
    w.SetSign(-x.GetSign());
    OutBSN(io, w);
}

namespace
{
    int const gRoundingModes[5] = { FE_TONEAREST, FE_DOWNWARD, FE_TOWARDZERO, FE_UPWARD, 12345 };
}

// Convert(input, precision, roundingMode, output). The precision is drawn
// relative to the input's bit count: anywhere in [1, numBits + 2], exactly
// numBits - 1 (a single remainder bit, the round-to-nearest tie), numBits,
// 1, small, and invalid (<= 0, LogError). Mode index 4 is not a <cfenv> mode
// (LogError once rounding is actually needed).
ORACLE_CASE("BSNumber.convert")
{
    BSN x = DrawBSN(io, 0, 8);
    int numBits = x.GetUInteger().GetNumBits();
    int p = 1;
    switch (io.rawInteger(0, 5))
    {
    case 0: p = io.rawInteger(1, numBits + 2); break;
    case 1: p = (numBits > 1 ? numBits - 1 : 1); break;
    case 2: p = (numBits > 0 ? numBits : 1); break;
    case 3: p = 1; break;
    case 4: p = io.rawInteger(1, 24); break;
    default: p = io.rawInteger(-3, 0); break;
    }
    int precision = static_cast<int>(io.given(static_cast<double>(p)));
    int modeIndex = io.integer(0, 4);
    BSN output;
    Convert(x, precision, gRoundingModes[modeIndex], output);
    OutBSN(io, output);
}

// The std:: and gte:: overloads that are exact or correctly rounded: fabs,
// frexp and ldexp (exact on the BSNumber), floor, ceil, sqrt, fmod,
// remainder (exact IEEE operations on the converted doubles), clamp,
// invsqrt, isign, saturate, sign, sqr, FMA, RobustSOP and RobustDOP.
ORACLE_CASE("BSNumber.std.exact")
{
    BSN x = DrawBSN(io, 0, 8);
    BSN y = DrawBSN(io, 0, 8);
    BSN z = DrawBSN(io, 0, 8);
    BSN w = DrawBSN(io, 0, 8);
    int k = io.integer(-1100, 1100);
    OutBSN(io, std::fabs(x));
    int32_t exponent = 0;
    BSN fr = std::frexp(x, &exponent);
    OutBSN(io, fr);
    io.outInt(exponent);
    OutBSN(io, std::ldexp(x, k));
    OutBSN(io, std::floor(x));
    OutBSN(io, std::ceil(x));
    OutBSN(io, std::sqrt(x));
    OutBSN(io, std::fmod(x, y));
    OutBSN(io, std::remainder(x, y));
    OutBSN(io, gte::clamp(x, y, z));
    OutBSN(io, gte::invsqrt(x));
    io.outInt(gte::isign(x));
    OutBSN(io, gte::saturate(x));
    OutBSN(io, gte::sign(x));
    OutBSN(io, gte::sqr(x));
    OutBSN(io, FMA(x, y, z));
    OutBSN(io, RobustSOP(x, y, z, w));
    OutBSN(io, RobustDOP(x, y, z, w));
}

// std::remainder with quotients far beyond 2^53, where only an exact
// quotient gives the IEEE result (the port computes it with bigint).
ORACLE_CASE("BSNumber.std.remainder")
{
    double a = io.given(DrawNormal(io, -40, 1000));
    double b = io.given(DrawNormal(io, -60, 60));
    BSN x(a), y(b);
    OutBSN(io, std::remainder(x, y));
    OutBSN(io, std::remainder(y, x));
    OutBSN(io, std::fmod(x, y));
}

// The std:: and gte:: overloads that round-trip through a C math library
// call: acos acosh asin asinh atan atanh atan2 cos cosh exp exp2 log log2
// log10 pow sin sinh tan tanh (std::) and atandivpi atan2divpi cospi exp10
// sinpi (gte::, built on std::atan, std::atan2, std::cos, std::exp and
// std::sin). The MSVC runtime and V8 may differ by an ulp; compared with the
// default scaled tolerance 1e-12 on the converted double. Arguments are kept
// in [-12, 12] (plus [-1,1], [1,40] and signed zeros) so that no result is
// ill-conditioned in its argument.
ORACLE_CASE("BSNumber.std.libm")
{
    double v[2];
    for (int i = 0; i < 2; ++i)
    {
        int mode = io.rawInteger(0, 4);
        double r = 0.0;
        switch (mode)
        {
        case 0: r = io.raw(-1.0, 1.0); break;
        case 1: r = io.raw(-12.0, 12.0); break;
        case 2: r = io.raw(1.0, 40.0); break;
        case 3: r = static_cast<double>(io.rawInteger(-4, 4)) * 0.25; break;
        default: r = RawSign(io) * 0.0; break;
        }
        v[i] = io.given(r);
    }
    BSN x(v[0]), y(v[1]);
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

// ============================================================ BSPrecision

namespace
{
    void OutParameters(oracle::Ctx& io, BSPrecision::Parameters const& p)
    {
        io.outInt(p.minExponent);
        io.outInt(p.maxExponent);
        io.outInt(p.maxBits);
        io.outInt(p.maxWords);
    }

    void OutBSP(oracle::Ctx& io, BSPrecision const& p)
    {
        OutParameters(io, p.bsn);
        OutParameters(io, p.bsr);
    }

    // A BSPrecision leaf, recorded as a kind followed by its arguments.
    //   0  BSPrecision(Type) for each of the six types
    //   1  BSPrecision(minExponent, maxExponent, maxBits), wide lattice
    //      (maxBits may be <= 0, which GetMaxWords truncates toward zero)
    //   2  the same on a small lattice, where operator+'s two branches and
    //      its carry tests tie often
    //   3  the finding #366 operands: IS_DOUBLE, IS_FLOAT, (-17, 20, 24),
    //      (-5, 20, 12)
    BSPrecision DrawBSP(oracle::Ctx& io)
    {
        int kind = io.integer(0, 3);
        if (kind == 0)
        {
            int t = io.integer(0, 5);
            return BSPrecision(static_cast<BSPrecision::Type>(t));
        }
        if (kind == 1)
        {
            int mn = io.integer(-1100, 50);
            int mx = io.integer(-50, 1100);
            int bits = io.integer(-40, 2200);
            return BSPrecision(mn, mx, bits);
        }
        if (kind == 2)
        {
            int mn = io.integer(-20, 0);
            int mx = io.integer(0, 24);
            int bits = io.integer(1, 30);
            return BSPrecision(mn, mx, bits);
        }
        int which = io.integer(0, 3);
        switch (which)
        {
        case 0: return BSPrecision(BSPrecision::Type::IS_DOUBLE);
        case 1: return BSPrecision(BSPrecision::Type::IS_FLOAT);
        case 2: return BSPrecision(-17, 20, 24);
        default: return BSPrecision(-5, 20, 12);
        }
    }

    // A random expression tree of depth <= 'depth' over DrawBSP leaves. The
    // operator code is recorded before the operands.
    BSPrecision DrawBSPExpression(oracle::Ctx& io, int depth)
    {
        int op = (depth == 0 ? 0 : io.integer(0, 6));
        if (op == 0)
        {
            return DrawBSP(io);
        }
        BSPrecision a = DrawBSPExpression(io, depth - 1);
        BSPrecision b = DrawBSPExpression(io, depth - 1);
        switch (op)
        {
        case 1: return a + b;
        case 2: return a - b;
        case 3: return a * b;
        case 4: return a / b;
        case 5: return a == b;
        default: return a < b;
        }
    }
}

// Every operator, both orders of the non-symmetric operator+, and the
// constructors. Finding #366 (operator+ pairs one operand's maxExponent with
// the other's minExponent) is preserved and compared bit for bit.
ORACLE_CASE("BSPrecision.operators")
{
    BSPrecision a = DrawBSP(io);
    BSPrecision b = DrawBSP(io);
    OutBSP(io, a);
    OutBSP(io, b);
    OutBSP(io, a + b);
    OutBSP(io, b + a);
    OutBSP(io, a - b);
    OutBSP(io, a * b);
    OutBSP(io, a / b);
    OutBSP(io, a == b);
    OutBSP(io, a != b);
    OutBSP(io, a < b);
    OutBSP(io, a <= b);
    OutBSP(io, a > b);
    OutBSP(io, a >= b);
    OutBSP(io, BSPrecision());
    io.outInt(BSPrecision::Parameters().GetMaxWords());
}

// Expression trees of depth <= 3, the way PrimalQuery2/3 chain the operators.
ORACLE_CASE("BSPrecision.expression")
{
    BSPrecision r = DrawBSPExpression(io, 3);
    OutBSP(io, r);
}

// =============================================================== QFNumber

namespace
{
    // Coefficient modes of the quadratic-field generators:
    //   0  uniform [-4, 4]
    //   1  integers [-5, 5]
    //   2  {-2, -1, -0, +0, 1, 2}
    double DrawQFCoefficient(oracle::Ctx& io, int mode)
    {
        switch (mode)
        {
        case 0: return io.raw(-4.0, 4.0);
        case 1: return static_cast<double>(io.rawInteger(-5, 5));
        default:
        {
            static double const values[6] = { -2.0, -1.0, -0.0, 0.0, 1.0, 2.0 };
            return values[io.rawInteger(0, 5)];
        }
        }
    }

    // A d-value: uniform [0, 10], an integer in [0, 12] (perfect squares
    // included), or a signed zero.
    double DrawQFD(oracle::Ctx& io, int mode)
    {
        switch (mode)
        {
        case 0: return io.raw(0.0, 10.0);
        case 1: return static_cast<double>(io.rawInteger(0, 12));
        default: return RawSign(io) * 0.0;
        }
    }

    // An unrecorded QFNumber<double, N>; d[0] is the innermost d-value.
    template <size_t N>
    QFNumber<double, N> RawQF(oracle::Ctx& io, double const* d, int mode)
    {
        if constexpr (N == 1)
        {
            double x0 = DrawQFCoefficient(io, mode);
            double x1 = DrawQFCoefficient(io, mode);
            return QFNumber<double, 1>(x0, x1, d[0]);
        }
        else
        {
            QFNumber<double, N - 1> x0 = RawQF<N - 1>(io, d, mode);
            QFNumber<double, N - 1> x1 = RawQF<N - 1>(io, d, mode);
            return QFNumber<double, N>(x0, x1, d[N - 1]);
        }
    }

    // Record x[0], x[1], d recursively (the replay reads the same order).
    template <size_t N>
    void RecordQF(oracle::Ctx& io, QFNumber<double, N> const& q)
    {
        if constexpr (N == 1)
        {
            io.given(q.x[0]);
            io.given(q.x[1]);
        }
        else
        {
            RecordQF<N - 1>(io, q.x[0]);
            RecordQF<N - 1>(io, q.x[1]);
        }
        io.given(q.d);
    }

    template <size_t N>
    void OutQF(oracle::Ctx& io, QFNumber<double, N> const& q)
    {
        if constexpr (N == 1)
        {
            io.outReal(q.x[0]);
            io.outReal(q.x[1]);
        }
        else
        {
            OutQF<N - 1>(io, q.x[0]);
            OutQF<N - 1>(io, q.x[1]);
        }
        io.outReal(q.d);
    }

    // A pair of QFNumber<double, N> that share their d-values, recorded.
    //   0..2  coefficient modes 0..2 with independent q0, q1
    //   3     equal values in different representations: the outer d is a
    //         perfect square k^2 and q1 = (q0.x[0] + k*T, q0.x[1] - T)
    //   4     the same x[1] (the d == 0 || x[1] == x[1] branch)
    //   5     q1 = q0 +- one unit in x[0] after mode 3 (near equality)
    //   6     mismatched outer d-values (arithmetic uses q0.d)
    template <size_t N>
    void DrawQFPair(oracle::Ctx& io, int maxMode, QFNumber<double, N>& q0,
        QFNumber<double, N>& q1, double& s)
    {
        int mode = io.rawInteger(0, maxMode);
        int cmode = (mode <= 2 ? mode : 1);
        double d[N];
        for (size_t i = 0; i < N; ++i)
        {
            d[i] = DrawQFD(io, mode == 2 ? io.rawInteger(1, 2) : (mode == 0 ? 0 : 1));
        }
        int k = io.rawInteger(1, 4);
        if (mode == 3 || mode == 5) { d[N - 1] = static_cast<double>(k * k); }
        q0 = RawQF<N>(io, d, cmode);
        q1 = RawQF<N>(io, d, cmode);
        if (mode == 3 || mode == 5)
        {
            if constexpr (N == 1)
            {
                double t = static_cast<double>(io.rawInteger(-3, 3));
                q1 = QFNumber<double, 1>(q0.x[0] + static_cast<double>(k) * t, q0.x[1] - t, d[0]);
                if (mode == 5) { q1.x[0] = q1.x[0] + RawSign(io); }
            }
            else
            {
                QFNumber<double, N - 1> T = RawQF<N - 1>(io, d, 1);
                q1 = QFNumber<double, N>(q0.x[0] + T * static_cast<double>(k), q0.x[1] - T, d[N - 1]);
                if (mode == 5) { q1.x[0] = q1.x[0] + RawSign(io); }
            }
        }
        else if (mode == 4)
        {
            q1.x[1] = q0.x[1];
        }
        else if (mode == 6)
        {
            q1.d = q0.d + 1.0 + static_cast<double>(io.rawInteger(0, 3));
        }
        int smode = io.rawInteger(0, 2);
        s = (smode == 0 ? io.raw(-4.0, 4.0) : (smode == 1 ? static_cast<double>(io.rawInteger(-3, 3))
            : RawSign(io) * 0.0));
        RecordQF<N>(io, q0);
        RecordQF<N>(io, q1);
        io.given(s);
    }

    // Every arithmetic operator: unary + and -, q+q, q+s, s+q, q-q, q-s,
    // s-q, q*q, q*s, s*q, q/q, q/s, s/q, then the compound assignments in
    // sequence.
    template <size_t N>
    void QFArithmetic(oracle::Ctx& io)
    {
        QFNumber<double, N> q0, q1;
        double s = 0.0;
        DrawQFPair<N>(io, 6, q0, q1, s);
        OutQF<N>(io, +q0);
        OutQF<N>(io, -q0);
        OutQF<N>(io, q0 + q1);
        OutQF<N>(io, q0 + s);
        OutQF<N>(io, s + q0);
        OutQF<N>(io, q0 - q1);
        OutQF<N>(io, q0 - s);
        OutQF<N>(io, s - q0);
        OutQF<N>(io, q0 * q1);
        OutQF<N>(io, q0 * s);
        OutQF<N>(io, s * q0);
        OutQF<N>(io, q0 / q1);
        OutQF<N>(io, q0 / s);
        OutQF<N>(io, s / q0);
        QFNumber<double, N> c = q0;
        c += q1;
        OutQF<N>(io, c);
        c += s;
        OutQF<N>(io, c);
        c -= q1;
        OutQF<N>(io, c);
        c -= s;
        OutQF<N>(io, c);
        c *= q1;
        OutQF<N>(io, c);
        c *= s;
        OutQF<N>(io, c);
        c /= q1;
        OutQF<N>(io, c);
        c /= s;
        OutQF<N>(io, c);
    }

    template <size_t N>
    void QFCompare(oracle::Ctx& io)
    {
        QFNumber<double, N> q0, q1;
        double s = 0.0;
        DrawQFPair<N>(io, 5, q0, q1, s);
        io.outBool(q0 == q1);
        io.outBool(q0 != q1);
        io.outBool(q0 < q1);
        io.outBool(q0 <= q1);
        io.outBool(q0 > q1);
        io.outBool(q0 >= q1);
        io.outBool(q1 < q0);
        io.outBool(q1 == q0);
    }
}

ORACLE_CASE("QFNumber.arithmetic.n1") { QFArithmetic<1>(io); }
ORACLE_CASE("QFNumber.arithmetic.n2") { QFArithmetic<2>(io); }
ORACLE_CASE("QFNumber.arithmetic.n3") { QFArithmetic<3>(io); }
ORACLE_CASE("QFNumber.compare.n1") { QFCompare<1>(io); }
ORACLE_CASE("QFNumber.compare.n2") { QFCompare<2>(io); }
ORACLE_CASE("QFNumber.compare.n3") { QFCompare<3>(io); }

// The constructors: QFNumber(), QFNumber(d), QFNumber(x0, x1, d) and
// QFNumber(std::array, d) for N = 1, and the coefficient and array
// constructors for N = 2. (The N >= 2 default and QFNumber(d) constructors
// create depth-(N-1) zero coefficients; the port's coefficients carry their
// depth at run time, so it builds those with the coefficient constructor.)
ORACLE_CASE("QFNumber.construct")
{
    double d0 = io.real(0.0, 10.0);
    double d1 = io.real(0.0, 10.0);
    double a = io.real(-4.0, 4.0);
    double b = io.real(-4.0, 4.0);
    OutQF<1>(io, QFNumber<double, 1>());
    OutQF<1>(io, QFNumber<double, 1>(d0));
    OutQF<1>(io, QFNumber<double, 1>(a, b, d0));
    OutQF<1>(io, QFNumber<double, 1>(std::array<double, 2>{ b, a }, d0));
    QFNumber<double, 1> inner(a, b, d0);
    OutQF<2>(io, QFNumber<double, 2>(inner, -inner, d1));
    OutQF<2>(io, QFNumber<double, 2>(std::array<QFNumber<double, 1>, 2>{ -inner, inner }, d1));
}

// ============================================================= SWInterval

namespace
{
    using SWI = SWInterval<double>;

    void OutSWI(oracle::Ctx& io, SWI const& w)
    {
        io.outReal(w[0]);
        io.outReal(w[1]);
    }

    // A scalar for the interval generators: uniform, integers, signed zeros,
    // subnormals, values near the overflow threshold, infinities.
    double DrawSWScalar(oracle::Ctx& io)
    {
        double s = RawSign(io);
        switch (io.rawInteger(0, 7))
        {
        case 0: return io.raw(-10.0, 10.0);
        case 1: return static_cast<double>(io.rawInteger(-4, 4));
        case 2: return s * 0.0;
        case 3: return s * static_cast<double>(io.rawInteger(1, 16)) * Pow2(-1074);
        case 4: return s * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(1000, 1023));
        case 5: return s * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(-1022, -1015));
        case 6: return (io.rawInteger(0, 7) == 0 ? s * std::numeric_limits<double>::infinity() : io.raw(-1.0e3, 1.0e3));
        default: return s * Pow2(io.rawInteger(-8, 8));
        }
    }

    // An interval [e0, e1], e0 <= e1, recorded as e0 then e1.
    //   0 positive, 1 negative, 2 straddling zero, 3 a zero endpoint (+0 or
    //   -0 on either side, [0,0], [-0,+0]), 4 degenerate [a,a], 5 integer
    //   lattice [-3,3], 6 near the overflow threshold, 7 subnormal, 8 an
    //   infinite endpoint.
    SWI DrawSWI(oracle::Ctx& io)
    {
        double e0 = 0.0, e1 = 0.0;
        int mode = io.rawInteger(0, 8);
        switch (mode)
        {
        case 0: e0 = io.raw(0.1, 10.0); e1 = e0 + io.raw(0.0, 5.0); break;
        case 1: e1 = -io.raw(0.1, 10.0); e0 = e1 - io.raw(0.0, 5.0); break;
        case 2: e0 = -io.raw(0.1, 10.0); e1 = io.raw(0.1, 10.0); break;
        case 3:
        {
            double z0 = RawSign(io) * 0.0;
            double z1 = RawSign(io) * 0.0;
            switch (io.rawInteger(0, 2))
            {
            case 0: e0 = z0; e1 = io.raw(0.1, 10.0); break;
            case 1: e0 = -io.raw(0.1, 10.0); e1 = z1; break;
            default: e0 = z0; e1 = z1; break;
            }
            break;
        }
        case 4: e0 = io.raw(-10.0, 10.0); e1 = e0; break;
        case 5:
        {
            int a = io.rawInteger(-3, 3);
            int b = io.rawInteger(-3, 3);
            e0 = static_cast<double>(a < b ? a : b);
            e1 = static_cast<double>(a < b ? b : a);
            break;
        }
        case 6:
        {
            double s = RawSign(io);
            double a = s * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(1015, 1023));
            double b = s * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(1015, 1023));
            e0 = (a < b ? a : b);
            e1 = (a < b ? b : a);
            break;
        }
        case 7:
        {
            double a = static_cast<double>(io.rawInteger(-16, 16)) * Pow2(-1074);
            double b = static_cast<double>(io.rawInteger(-16, 16)) * Pow2(-1074);
            e0 = (a < b ? a : b);
            e1 = (a < b ? b : a);
            break;
        }
        default:
        {
            double inf = std::numeric_limits<double>::infinity();
            if (io.rawInteger(0, 1) != 0) { e0 = -inf; e1 = io.raw(-10.0, 10.0); }
            else { e0 = io.raw(-10.0, 10.0); e1 = inf; }
            break;
        }
        }
        double r0 = io.given(e0);
        double r1 = io.given(e1);
        return SWI(r0, r1);
    }
}

// The leaf-node operations Add, Sub, Mul, Div(u, v), including Div by +-0
// (the whole real line).
ORACLE_CASE("SWInterval.leaf")
{
    double u = io.given(DrawSWScalar(io));
    double v = io.given(DrawSWScalar(io));
    OutSWI(io, SWI::Add(u, v));
    OutSWI(io, SWI::Sub(u, v));
    OutSWI(io, SWI::Mul(u, v));
    OutSWI(io, SWI::Div(u, v));
    OutSWI(io, SWI::Div(v, u));
}

// The internal four-argument operations, Mul2, Reciprocal, ReciprocalDown,
// ReciprocalUp and Reals, on arbitrary (not necessarily ordered) operands.
ORACLE_CASE("SWInterval.internal")
{
    double u0 = io.given(DrawSWScalar(io));
    double u1 = io.given(DrawSWScalar(io));
    double v0 = io.given(DrawSWScalar(io));
    double v1 = io.given(DrawSWScalar(io));
    OutSWI(io, SWI::Add(u0, u1, v0, v1));
    OutSWI(io, SWI::Sub(u0, u1, v0, v1));
    OutSWI(io, SWI::Mul(u0, u1, v0, v1));
    OutSWI(io, SWI::Mul2(u0, u1, v0, v1));
    OutSWI(io, SWI::Div(u0, u1, v0, v1));
    OutSWI(io, SWI::Reciprocal(v0, v1));
    OutSWI(io, SWI::ReciprocalDown(v0));
    OutSWI(io, SWI::ReciprocalUp(v1));
    OutSWI(io, SWI::Reals());
}

// The non-class operators on intervals and scalars in every sign
// configuration, the compound assignments in sequence, and the accessors.
ORACLE_CASE("SWInterval.operators")
{
    SWI u = DrawSWI(io);
    SWI v = DrawSWI(io);
    double s = io.given(DrawSWScalar(io));
    OutSWI(io, +u);
    OutSWI(io, -u);
    OutSWI(io, u + v);
    OutSWI(io, u + s);
    OutSWI(io, s + u);
    OutSWI(io, u - v);
    OutSWI(io, u - s);
    OutSWI(io, s - u);
    OutSWI(io, u * v);
    OutSWI(io, v * u);
    OutSWI(io, u * s);
    OutSWI(io, s * u);
    OutSWI(io, u / v);
    OutSWI(io, u / s);
    OutSWI(io, s / u);
    SWI c = u;
    c += v;
    OutSWI(io, c);
    c += s;
    OutSWI(io, c);
    c -= v;
    OutSWI(io, c);
    c -= s;
    OutSWI(io, c);
    c *= v;
    OutSWI(io, c);
    c *= s;
    OutSWI(io, c);
    c /= v;
    OutSWI(io, c);
    c /= s;
    OutSWI(io, c);
    std::array<double, 2> ends = u.GetEndpoints();
    io.outReal(ends[0]);
    io.outReal(ends[1]);
    OutSWI(io, SWI());
    OutSWI(io, SWI(s));
    OutSWI(io, SWI(std::array<double, 2>{ v[0], v[1] }));
}

// std::nextafter(x, -max) and std::nextafter(x, +max) through Mul(x, 1),
// whose product is x exactly: signed zeros, the smallest and largest
// subnormals, the smallest normal, +-1, +-DBL_MAX, the infinities (finding
// #50: nextafter pulls them back to +-DBL_MAX), NaN and arbitrary bit
// patterns.
ORACLE_CASE("SWInterval.nextafter")
{
    static uint64_t const specials[10] = {
        0x0000000000000000ull, 0x0000000000000001ull, 0x0000000000000002ull,
        0x000FFFFFFFFFFFFFull, 0x0010000000000000ull, 0x3FF0000000000000ull,
        0x7FEFFFFFFFFFFFFFull, 0x7FF0000000000000ull, 0x7FF8000000000000ull,
        0x3FEFFFFFFFFFFFFFull };
    uint64_t bits = 0;
    if (io.rawInteger(0, 2) != 0)
    {
        bits = specials[io.rawInteger(0, 9)];
    }
    else
    {
        bits = RawBits64(io);
    }
    if (io.rawInteger(0, 1) != 0) { bits |= 0x8000000000000000ull; }
    double x = io.given(FromBits(bits));
    OutSWI(io, SWI::Mul(x, 1.0));
}

// Finding #50 (preserved): a bound that overflows to an infinity is pulled
// back to +-DBL_MAX by the nextafter widening, and dividing by an interval
// with a zero endpoint yields a finite MAX bound. Both implementations do
// it, so the case is compared bit for bit.
ORACLE_CASE("SWInterval.overflow")
{
    double a = io.given(RawSign(io) * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(1015, 1023)));
    double b = io.given(RawSign(io) * std::ldexp(1.0 + io.raw(0.0, 1.0), io.rawInteger(1015, 1023)));
    double t = io.given(RawSign(io) * static_cast<double>(io.rawInteger(1, 64)) * Pow2(-1074));
    OutSWI(io, SWI::Add(a, b));
    OutSWI(io, SWI::Mul(a, b));
    OutSWI(io, SWI::Div(a, t));
    OutSWI(io, SWI::Div(8.881784197001252e-16, -std::numeric_limits<double>::denorm_min()));
    SWI p(1.0, 2.0), q(0.0, 4.0);
    OutSWI(io, p / q);
    OutSWI(io, SWI(a, a) * SWI(b, b));
}

// ================================================================ BSplines

namespace
{
    // A drawn BasisFunctionInput, as in v17-curves.cpp. Knot values are
    // dyadic so knot differences and the domain are exact.
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
        int32_t maxExtra)
    {
        int32_t mode = io.rawInteger(0, 4);
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

    void RecordPoints(oracle::Ctx& io, std::vector<double> const& P)
    {
        for (double v : P) { io.given(v); }
    }
}

// ======================================================== BSplineCurveFit

// The least-squares fit (A^T*A banded, the Cholesky solve of
// BandedMatrix::SolveSystem<true>, the X0*P accumulation and the end-point
// override), every accessor, and Evaluate/GetPosition for orders 0..3 at
// parameters inside, at the knots and outside [0,1] (clamped). One third of
// the records use the minimum number of samples numControls + degree + 1.
// Sample modes: uniform, lattice, all samples equal, collinear, signed
// zeros, and a polynomial curve of degree <= degree (which the spline space
// contains, so the fit must reproduce it; checked against the polynomial
// outside the harness).
ORACLE_CASE("BSplineCurveFit.fit")
{
    int dimension = io.integer(1, 4);
    int degree = io.integer(1, 5);
    int numControls = io.integer(degree + 1, degree + 6);
    int minSamples = numControls + degree + 1;
    int ns = (io.rawInteger(0, 2) == 0 ? minSamples : io.rawInteger(minSamples, minSamples + 10));
    int numSamples = static_cast<int>(io.given(static_cast<double>(ns)));
    int pmode = io.rawInteger(0, 5);
    std::vector<double> P;
    if (pmode <= 4)
    {
        P = DrawPoints(io, numSamples, dimension, pmode);
    }
    else
    {
        P.resize(static_cast<size_t>(numSamples) * dimension);
        for (int j = 0; j < dimension; ++j)
        {
            std::vector<double> c(static_cast<size_t>(degree) + 1);
            for (auto& ck : c) { ck = static_cast<double>(io.rawInteger(-3, 3)); }
            for (int i = 0; i < numSamples; ++i)
            {
                double t = static_cast<double>(i) / static_cast<double>(numSamples - 1);
                double value = 0.0;
                for (int k = degree; k >= 0; --k) { value = value * t + c[k]; }
                P[static_cast<size_t>(i) * dimension + j] = value;
            }
        }
    }
    RecordPoints(io, P);
    BSplineCurveFit<double> fit(dimension, numSamples, P.data(), degree, numControls);
    io.outInt(fit.GetDimension());
    io.outInt(fit.GetNumSamples());
    io.outInt(fit.GetDegree());
    io.outInt(fit.GetNumControls());
    io.outReal(fit.GetSampleData()[0]);
    double const* controls = fit.GetControlData();
    for (int i = 0; i < dimension * numControls; ++i) { io.outReal(controls[i]); }
    io.outReal(fit.GetBasis().GetMinDomain());
    io.outReal(fit.GetBasis().GetMaxDomain());
    io.outInt(fit.GetBasis().GetNumKnots());
    std::vector<double> value(dimension);
    for (int k = 0; k < 6; ++k)
    {
        double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
        int order = io.integer(0, 3);
        fit.Evaluate(t, static_cast<uint32_t>(order), value.data());
        for (int j = 0; j < dimension; ++j) { io.outReal(value[j]); }
    }
    double t = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
    fit.GetPosition(t, value.data());
    for (int j = 0; j < dimension; ++j) { io.outReal(value[j]); }
}

// Throw parity for the constructor preconditions and for Evaluate with
// order 4 (BasisFunction::Evaluate asserts order <= 3).
//   0 dimension 0, 1 degree 0, 2 degree >= numControls, 3 fewer samples
//   than numControls + degree + 1, 4 a valid fit evaluated at order 4
ORACLE_CASE("BSplineCurveFit.invalid")
{
    int mode = io.integer(0, 4);
    int d = io.rawInteger(1, 3);
    int n = d + 1 + io.rawInteger(0, 3);
    int s = n + d + 1 + io.rawInteger(0, 3);
    int dim = io.rawInteger(1, 3);
    switch (mode)
    {
    case 0: dim = 0; break;
    case 1: d = 0; break;
    case 2: n = d - io.rawInteger(0, 1); if (n < 1) { n = 1; } break;
    case 3: s = n + d - io.rawInteger(0, 2); break;
    default: break;
    }
    int dimension = static_cast<int>(io.given(static_cast<double>(dim)));
    int degree = static_cast<int>(io.given(static_cast<double>(d)));
    int numControls = static_cast<int>(io.given(static_cast<double>(n)));
    int numSamples = static_cast<int>(io.given(static_cast<double>(s)));
    std::vector<double> P = DrawPoints(io, numSamples, dimension, 1);
    RecordPoints(io, P);
    P.push_back(0.0);
    BSplineCurveFit<double> fit(dimension, numSamples, P.data(), degree, numControls);
    std::vector<double> value(dimension);
    fit.Evaluate(0.5, 4, value.data());
    io.outReal(value[0]);
}

// ====================================================== BSplineSurfaceFit

namespace
{
    // Sample modes 0..4 as DrawPoints; mode 5 samples a bilinear polynomial
    // in (u, v) = (i0/(numSamples0-1), i1/(numSamples1-1)), which every
    // degree >= 1 spline space contains.
    std::vector<double> DrawSurfaceSamples(oracle::Ctx& io, int ns0, int ns1)
    {
        int pmode = io.rawInteger(0, 5);
        if (pmode <= 4)
        {
            return DrawPoints(io, ns0 * ns1, 3, pmode);
        }
        std::vector<double> P(static_cast<size_t>(ns0) * ns1 * 3);
        for (int j = 0; j < 3; ++j)
        {
            double c[4];
            for (int k = 0; k < 4; ++k) { c[k] = static_cast<double>(io.rawInteger(-3, 3)); }
            for (int i1 = 0; i1 < ns1; ++i1)
            {
                double v = static_cast<double>(i1) / static_cast<double>(ns1 - 1);
                for (int i0 = 0; i0 < ns0; ++i0)
                {
                    double u = static_cast<double>(i0) / static_cast<double>(ns0 - 1);
                    P[(static_cast<size_t>(i0) + static_cast<size_t>(ns0) * i1) * 3 + j] =
                        c[0] + c[1] * u + c[2] * v + c[3] * u * v;
                }
            }
        }
        return P;
    }

    std::vector<Vector3<double>> ToVector3(std::vector<double> const& P)
    {
        std::vector<Vector3<double>> V(P.size() / 3);
        for (size_t i = 0; i < V.size(); ++i)
        {
            V[i] = { P[3 * i], P[3 * i + 1], P[3 * i + 2] };
        }
        return V;
    }
}

// The tensor-product least-squares fit Q = X0*P*X1^T, the accessors and
// GetPosition inside, at and outside [0,1]^2. numSamples ranges down to
// numControls (the square, interpolating case).
ORACLE_CASE("BSplineSurfaceFit.fit")
{
    int degree0 = io.integer(1, 3);
    int numControls0 = io.integer(degree0 + 2, degree0 + 4);
    int numSamples0 = io.integer(numControls0, numControls0 + 3);
    int degree1 = io.integer(1, 3);
    int numControls1 = io.integer(degree1 + 2, degree1 + 4);
    int numSamples1 = io.integer(numControls1, numControls1 + 3);
    std::vector<double> P = DrawSurfaceSamples(io, numSamples0, numSamples1);
    RecordPoints(io, P);
    std::vector<Vector3<double>> samples = ToVector3(P);
    BSplineSurfaceFit<double> fit(degree0, numControls0, numSamples0,
        degree1, numControls1, numSamples1, samples.data());
    for (int d = 0; d < 2; ++d)
    {
        io.outInt(fit.GetNumSamples(d));
        io.outInt(fit.GetDegree(d));
        io.outInt(fit.GetNumControls(d));
        io.outReal(fit.GetBasis(d).GetMinDomain());
        io.outReal(fit.GetBasis(d).GetMaxDomain());
    }
    io.outVec(fit.GetSampleData()[0]);
    Vector3<double> const* controls = fit.GetControlData();
    for (int i = 0; i < numControls0 * numControls1; ++i) { io.outVec(controls[i]); }
    for (int k = 0; k < 5; ++k)
    {
        double u = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
        double v = io.given(DrawParameter(io, 0.0, 1.0, nullptr));
        io.outVec(fit.GetPosition(u, v));
    }
}

// Throw parity for the constructor preconditions: 0 degree0 = 0,
// 1 degree0 + 1 >= numControls0, 2 numControls0 > numSamples0, 3..5 the
// same for the second dimension.
ORACLE_CASE("BSplineSurfaceFit.invalid")
{
    int mode = io.integer(0, 5);
    int d[2], n[2], s[2];
    for (int k = 0; k < 2; ++k)
    {
        d[k] = io.rawInteger(1, 2);
        n[k] = d[k] + 2 + io.rawInteger(0, 2);
        s[k] = n[k] + io.rawInteger(0, 2);
    }
    int k = mode / 3;
    switch (mode % 3)
    {
    case 0: d[k] = 0; break;
    case 1: n[k] = d[k] + 1 - io.rawInteger(0, 1); break;
    default: s[k] = n[k] - io.rawInteger(1, 2); break;
    }
    int degree0 = static_cast<int>(io.given(static_cast<double>(d[0])));
    int numControls0 = static_cast<int>(io.given(static_cast<double>(n[0])));
    int numSamples0 = static_cast<int>(io.given(static_cast<double>(s[0])));
    int degree1 = static_cast<int>(io.given(static_cast<double>(d[1])));
    int numControls1 = static_cast<int>(io.given(static_cast<double>(n[1])));
    int numSamples1 = static_cast<int>(io.given(static_cast<double>(s[1])));
    std::vector<double> P = DrawPoints(io, numSamples0 * numSamples1, 3, 1);
    RecordPoints(io, P);
    std::vector<Vector3<double>> samples = ToVector3(P);
    samples.push_back(Vector3<double>{ 0.0, 0.0, 0.0 });
    BSplineSurfaceFit<double> fit(degree0, numControls0, numSamples0,
        degree1, numControls1, numSamples1, samples.data());
    io.outInt(fit.GetNumControls(0));
}

// ========================================================= BSplineSurface

namespace
{
    // Controls of dimension N (count 'count'), mode 0..4 as DrawPoints, or
    // mode 5 = none (the null pointer: zero controls), recorded.
    template <int N>
    std::vector<Vector<N, double>> DrawControls(oracle::Ctx& io, int32_t count, int32_t mode)
    {
        std::vector<Vector<N, double>> C;
        if (mode <= 4)
        {
            std::vector<double> P = DrawPoints(io, count, N, mode);
            RecordPoints(io, P);
            C.resize(static_cast<size_t>(count));
            for (int32_t i = 0; i < count; ++i)
            {
                for (int j = 0; j < N; ++j) { C[i][j] = P[static_cast<size_t>(i) * N + j]; }
            }
        }
        return C;
    }

    template <int N>
    Vector<N, double> DrawLatticeVector(oracle::Ctx& io)
    {
        Vector<N, double> v;
        for (int j = 0; j < N; ++j) { v[j] = io.lattice(-4, 4); }
        return v;
    }

    // Evaluation orders: 0..3 compute a jet (order 3 evaluates the basis at
    // order 3 and still fills entries 0..5), 6 and 7 are >= SUP_ORDER and
    // return a zero jet. Orders 4 and 5 assert in BasisFunction and have
    // their own throw-parity case.
    int DrawSurfaceOrder(oracle::Ctx& io)
    {
        static int const orders[6] = { 0, 1, 2, 3, 6, 7 };
        return static_cast<int>(io.given(static_cast<double>(orders[io.rawInteger(0, 5)])));
    }

    template <int N>
    void SurfaceBody(oracle::Ctx& io, BasisSpec const& s0, BasisSpec const& s1, int cmode)
    {
        int32_t n0 = s0.numControls, n1 = s1.numControls;
        std::vector<Vector<N, double>> C = DrawControls<N>(io, n0 * n1, cmode);
        BasisFunctionInput<double> inputs[2] = { ToInput(s0), ToInput(s1) };
        BSplineSurface<N, double> surface(inputs, cmode <= 4 ? C.data() : nullptr);
        for (int k = 0; k < 3; ++k)
        {
            int i0 = io.integer(-1, n0);
            int i1 = io.integer(-1, n1);
            Vector<N, double> value = DrawLatticeVector<N>(io);
            surface.SetControl(i0, i1, value);
        }
        io.outBool(static_cast<bool>(surface));
        io.outReal(surface.GetUMin());
        io.outReal(surface.GetUMax());
        io.outReal(surface.GetVMin());
        io.outReal(surface.GetVMax());
        io.outBool(surface.IsRectangular());
        io.outInt(surface.GetNumControls(0));
        io.outInt(surface.GetNumControls(1));
        io.outReal(surface.GetBasisFunction(0).GetMaxDomain());
        io.outReal(surface.GetBasisFunction(1).GetMinDomain());
        for (int k = 0; k < 3; ++k)
        {
            int i0 = io.integer(-1, n0);
            int i1 = io.integer(-1, n1);
            io.outVec(surface.GetControl(i0, i1));
        }
        io.outVec(surface.GetControls()[n0 * n1 - 1]);
        for (int k = 0; k < 5; ++k)
        {
            double u = io.given(DrawParameter(io, surface.GetUMin(), surface.GetUMax(), &s0.knots));
            double v = io.given(DrawParameter(io, surface.GetVMin(), surface.GetVMax(), &s1.knots));
            int order = DrawSurfaceOrder(io);
            std::array<Vector<N, double>, 6> jet{};
            surface.Evaluate(u, v, static_cast<uint32_t>(order), jet.data());
            for (int i = 0; i < 6; ++i) { io.outVec(jet[i]); }
        }
        double u = io.given(DrawParameter(io, surface.GetUMin(), surface.GetUMax(), &s0.knots));
        double v = io.given(DrawParameter(io, surface.GetVMin(), surface.GetVMax(), &s1.knots));
        io.outVec(surface.GetPosition(u, v));
        io.outVec(surface.GetUTangent(u, v));
        io.outVec(surface.GetVTangent(u, v));
    }
}

// BSplineSurface<N, double> for N = 1, 2, 3 over open uniform, open
// nonuniform (with a repeated interior knot) and periodic (uniform and
// nonuniform) bases: construction with and without controls, SetControl
// and GetControl in and out of range, the ParametricSurface accessors,
// Evaluate for every order, GetPosition, GetUTangent and GetVTangent.
ORACLE_CASE("BSplineSurface.evaluate")
{
    int N = io.integer(1, 3);
    BasisSpec s0 = DrawBasis(io, 1, 3, 3);
    RecordBasis(io, s0);
    BasisSpec s1 = DrawBasis(io, 1, 3, 3);
    RecordBasis(io, s1);
    int cmode = io.integer(0, 5);
    if (N == 1) { SurfaceBody<1>(io, s0, s1, cmode); }
    else if (N == 2) { SurfaceBody<2>(io, s0, s1, cmode); }
    else { SurfaceBody<3>(io, s0, s1, cmode); }
}

// Evaluate with order 4 or 5 (below SUP_ORDER = 6, above BasisFunction's
// limit of 3): BasisFunction::Evaluate asserts.
ORACLE_CASE("BSplineSurface.evaluate.invalidOrder")
{
    BasisSpec s0 = DrawBasis(io, 1, 2, 2);
    RecordBasis(io, s0);
    BasisSpec s1 = DrawBasis(io, 1, 2, 2);
    RecordBasis(io, s1);
    std::vector<Vector<2, double>> C = DrawControls<2>(io, s0.numControls * s1.numControls, 1);
    BasisFunctionInput<double> inputs[2] = { ToInput(s0), ToInput(s1) };
    BSplineSurface<2, double> surface(inputs, C.data());
    int order = io.integer(4, 5);
    std::array<Vector<2, double>, 6> jet{};
    surface.Evaluate(0.5 * (surface.GetUMin() + surface.GetUMax()), surface.GetVMin(),
        static_cast<uint32_t>(order), jet.data());
    io.outVec(jet[0]);
}

// ========================================================== BSplineVolume

namespace
{
    // Orders 0..3 compute a jet (order 3 fills entries 0..9), 10 and 11 are
    // >= SUP_ORDER and return a zero jet; 4..9 assert in BasisFunction and
    // have their own throw-parity case.
    int DrawVolumeOrder(oracle::Ctx& io)
    {
        static int const orders[6] = { 0, 1, 2, 3, 10, 11 };
        return static_cast<int>(io.given(static_cast<double>(orders[io.rawInteger(0, 5)])));
    }

    template <int N>
    void VolumeBody(oracle::Ctx& io, BasisSpec const* s, int cmode)
    {
        int32_t n0 = s[0].numControls, n1 = s[1].numControls, n2 = s[2].numControls;
        std::vector<Vector<N, double>> C = DrawControls<N>(io, n0 * n1 * n2, cmode);
        BasisFunctionInput<double> inputs[3] = { ToInput(s[0]), ToInput(s[1]), ToInput(s[2]) };
        BSplineVolume<N, double> volume(inputs, cmode <= 4 ? C.data() : nullptr);
        for (int k = 0; k < 3; ++k)
        {
            int i0 = io.integer(-1, n0);
            int i1 = io.integer(-1, n1);
            int i2 = io.integer(-1, n2);
            Vector<N, double> value = DrawLatticeVector<N>(io);
            volume.SetControl(i0, i1, i2, value);
        }
        io.outBool(static_cast<bool>(volume));
        for (int d = 0; d < 3; ++d)
        {
            io.outReal(volume.GetMinDomain(d));
            io.outReal(volume.GetMaxDomain(d));
            io.outInt(volume.GetNumControls(d));
            io.outInt(volume.GetBasisFunction(d).GetNumKnots());
        }
        for (int k = 0; k < 3; ++k)
        {
            int i0 = io.integer(-1, n0);
            int i1 = io.integer(-1, n1);
            int i2 = io.integer(-1, n2);
            io.outVec(volume.GetControl(i0, i1, i2));
        }
        io.outVec(volume.GetControls()[n0 * n1 * n2 - 1]);
        for (int k = 0; k < 4; ++k)
        {
            double t[3];
            for (int d = 0; d < 3; ++d)
            {
                t[d] = io.given(DrawParameter(io, volume.GetMinDomain(d), volume.GetMaxDomain(d), &s[d].knots));
            }
            int order = DrawVolumeOrder(io);
            std::array<Vector<N, double>, 10> jet{};
            volume.Evaluate(t[0], t[1], t[2], static_cast<uint32_t>(order), jet.data());
            for (int i = 0; i < 10; ++i) { io.outVec(jet[i]); }
        }
    }
}

// BSplineVolume<N, double> for N = 1, 2, 3 with three independent bases of
// the five kinds, controls present or deferred, SetControl/GetControl in
// and out of range, the domain and basis accessors, and Evaluate for every
// order.
ORACLE_CASE("BSplineVolume.evaluate")
{
    int N = io.integer(1, 3);
    BasisSpec s[3];
    for (int d = 0; d < 3; ++d)
    {
        s[d] = DrawBasis(io, 1, 2, 2);
        RecordBasis(io, s[d]);
    }
    int cmode = io.integer(0, 5);
    if (N == 1) { VolumeBody<1>(io, s, cmode); }
    else if (N == 2) { VolumeBody<2>(io, s, cmode); }
    else { VolumeBody<3>(io, s, cmode); }
}

// Evaluate with 4 <= order <= 9 (below SUP_ORDER = 10, above BasisFunction's
// limit of 3): BasisFunction::Evaluate asserts.
ORACLE_CASE("BSplineVolume.evaluate.invalidOrder")
{
    BasisSpec s[3];
    for (int d = 0; d < 3; ++d)
    {
        s[d] = DrawBasis(io, 1, 1, 1);
        RecordBasis(io, s[d]);
    }
    std::vector<Vector<1, double>> C = DrawControls<1>(io,
        s[0].numControls * s[1].numControls * s[2].numControls, 1);
    BasisFunctionInput<double> inputs[3] = { ToInput(s[0]), ToInput(s[1]), ToInput(s[2]) };
    BSplineVolume<1, double> volume(inputs, C.data());
    int order = io.integer(4, 9);
    std::array<Vector<1, double>, 10> jet{};
    volume.Evaluate(volume.GetMinDomain(0), volume.GetMaxDomain(1), volume.GetMinDomain(2),
        static_cast<uint32_t>(order), jet.data());
    io.outVec(jet[0]);
}
