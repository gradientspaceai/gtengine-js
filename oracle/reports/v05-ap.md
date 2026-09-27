# Verify group 5 (arbitrary precision) against the MSVC build of upstream GTE

Family `v05-ap`, 35 cases, 20 golden records each. The group is
alphabetical rather than thematic: `BSNumber.h`, `BSPrecision.h`,
`QFNumber.h`, `SWInterval.h`, `BSplineCurveFit.h`, `BSplineSurface.h`,
`BSplineSurfaceFit.h`, `BSplineVolume.h`.

Deep run `npm run oracle:deep -- 2000 v05-ap`: 70 000 records (9703 of them
throw records, all with throw parity), **all 36 tests pass** (35 cases plus
the claim check). 33 cases are `{ exact: true }`, one is compared with the
default tolerance (`BSNumber.std.libm`, C math library) and one is a declared
deviation (`#95`). Outside the deviation case the deep run compares
2 158 091 floating-point outputs, of which 2 155 699 are bit-identical; the
2392 others all belong to `BSNumber.std.libm` (max scaled error 2.22e-16).
Every discrete output (signs, exponents, bit counts, 32-bit words, digests,
booleans, counts) matches exactly.

**One port defect was found and fixed** (`QFNumber.scalarDiv`, a sign of zero
for `N >= 2`; see below).

## How BSNumber is compared

The port replaces upstream's `UInteger` layer (`UIntegerAP32`,
`UIntegerFP32`, `UIntegerALU32`, all intentionally not ported) with `bigint`,
so a `BSNumber` is compared by value and canonical form, never by storage:
the C++ side instantiates `BSNumber<UIntegerAP32>` and `OutBSN` emits

* `GetSign()`, `GetBiasedExponent()`, `GetExponent()`,
* `GetUInteger().GetNumBits()` and the 32-bit words of the odd integer
  (all of them up to 8 words; beyond that the top three, the bottom two and
  a 48-bit FNV-1a digest of every word, so the goldens stay small while a
  one-bit difference anywhere is still detected),
* `operator double` and `operator float` (widened exactly to double).

The replay extracts the same words from the port's `bigint`. Every `BSNumber`
result in the family, from every case, therefore also tests both conversions.
No helper was added to `Oracle.h` or `harness.ts`; the encoding lives in the
case file and the replay.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `BSNumber.h` | `BSNumber.construct.double` (every finite double incl. subnormals and signed zeros, plus +-inf and NaN, the graceful-exit path), `.construct.float` (binary32 subnormals, `FLT_MIN`, `FLT_MAX`, inf, NaN), `.construct.integer` (`int32_t`, `uint32_t`, `int64_t`, `uint64_t`: small, powers of two, 2^k +- 1, 2^53 + 1, the extremes, full-width random), `.construct.string` (valid numbers of 1..80 digits incl. 2^64, 2^64 - 1, 2^128, 10^59, `"0"`, `"+0"`, `"-0"`; throw parity for leading zeros, embedded non-digits, `""`, `"+"`, `"-"`, decimal-point and exponent notation) | exact | pass; string: 1122 of 2000 throw on both sides |
| | `BSNumber.construct.string.singleChar` | deviation (#95) | 2000 of 2000 records deviate |
| | `BSNumber.compare` (`==`, `!=`, `<`, `<=`, `>`, `>=`, both orders) | exact | pass; 748 / 450 / 802 records with x < y / x = y / x > y; 346 records compare two numbers of equal sign and exponent whose odd integers have different lengths (the left-aligned word comparison of `UIntegerALU32::operator<`) |
| | `BSNumber.arithmetic` (unary `+` and `-`, binary `+ - *` in both orders, `+=`, `-=`, `*=`, `Negate()`) | exact | pass; 426 of 2000 sums exceed 256 bits (multi-word, up to ~2100 bits) |
| | `BSNumber.conversions.boundary` (`operator double` / `operator float` at binary64 and binary32 ties, just above and below them, between and halfway between the smallest subnormals, and at the overflow boundary) | exact | pass; binary64: 106 exact ties, 572 above, 618 below, 704 exact, 140 overflow to inf, 64 round to `DBL_TRUE_MIN`, 30 round to zero; binary32: 84 ties, 607 overflow |
| | `BSNumber.accessors` (`SetExponent`, `SetBiasedExponent`, `SetSign` and the getters) | exact | pass |
| | `BSNumber.convert` (`Convert` with `FE_TONEAREST`, `FE_DOWNWARD`, `FE_TOWARDZERO`, `FE_UPWARD` and a non-`<cfenv>` mode; precision anywhere in [1, numBits + 2], numBits - 1, numBits, 1, small, and <= 0) | exact | pass; 581 throw records (precision <= 0, or the invalid mode when rounding is needed), 479 no-ops, 60 exact nearest-mode ties, 141 carries into the next power of two (`RoundUp` shift) |
| | `BSNumber.std.exact` (`std::fabs`, `frexp`, `ldexp`, `floor`, `ceil`, `sqrt`, `fmod`, `remainder`; `gte::clamp`, `invsqrt`, `isign`, `saturate`, `sign`, `sqr`, `FMA`, `RobustSOP`, `RobustDOP`) | exact | pass |
| | `BSNumber.std.remainder` (`std::remainder` and `std::fmod` with quotients up to 2^1060) | exact | pass |
| | `BSNumber.std.libm` (`std::acos acosh asin asinh atan atanh atan2 cos cosh exp exp2 log log2 log10 pow sin sinh tan tanh`, `gte::atandivpi atan2divpi cospi exp10 sinpi`) | tolerance 1e-12 (C math library; arguments in [-12, 12]) | pass; 45 608 of 48 000 bit-identical, max scaled error 2.22e-16 |
| `BSPrecision.h` | `BSPrecision.operators` (the three constructors, `Parameters()`, `GetMaxWords` incl. `maxBits <= 0`, `+` in both orders, `- * / == != < <= > >=`; leaves from the six `Type`s, wide and small lattices, and the #366 operands) | exact | pass |
| | `BSPrecision.expression` (random operator trees of depth <= 3, the way `PrimalQuery2/3` chain them) | exact | pass |
| `QFNumber.h` | `QFNumber.arithmetic.n1`, `.n2`, `.n3` (`QFNumber<double, 1..3>`: unary `+ -`, `q+q`, `q+s`, `s+q`, `q-q`, `q-s`, `s-q`, `q*q`, `q*s`, `s*q`, `q/q`, `q/s`, `s/q`, `+= -= *= /=` with both operand kinds) | exact | pass (after the fix below) |
| | `QFNumber.compare.n1`, `.n2`, `.n3` (`== != < <= > >=`) | exact | pass; N = 1 branch split: 674 `d == 0 || x1 equal`, 215 / 486 `x1 >` early-out / squared test, 206 / 419 `x1 <` early-out / squared test |
| | `QFNumber.construct` | exact | pass |
| `SWInterval.h` | `SWInterval.leaf` (`Add Sub Mul Div(u, v)` incl. division by +-0) | exact | pass |
| | `SWInterval.internal` (four-argument `Add Sub Mul Div`, `Mul2`, `Reciprocal`, `ReciprocalDown`, `ReciprocalUp`, `Reals`) | exact | pass |
| | `SWInterval.operators` (every non-class operator on interval/interval, interval/scalar, scalar/interval, the compound assignments, `operator[]`, `GetEndpoints`, the constructors) | exact | pass; all nine sign branches of `operator*` reached (82 to 273 records each); 161 / 94 / 358 divisors with a zero lower endpoint / zero upper endpoint / straddling zero |
| | `SWInterval.nextafter` (`std::nextafter(x, +-max)` through `Mul(x, 1)`: +-0, +-`DBL_TRUE_MIN`, the largest subnormal, `DBL_MIN`, +-1, +-`DBL_MAX`, +-inf, NaN, random bit patterns) | exact | pass |
| | `SWInterval.overflow` (#50, preserved) | exact | pass |
| `BSplineCurveFit.h` | `BSplineCurveFit.fit` (dimension 1..4, degree 1..5, the minimum sample count on a third of the records; uniform, lattice, all-equal, collinear, signed-zero and polynomial samples; every accessor, `Evaluate` for orders 0..3 inside, at and outside [0,1], `GetPosition`) | exact | pass |
| | `BSplineCurveFit.invalid` (throw parity: dimension 0, degree 0, degree >= numControls, too few samples, `Evaluate` at order 4) | exact | 2000 of 2000 throw on both sides |
| `BSplineSurfaceFit.h` | `BSplineSurfaceFit.fit` (degrees 1..3, numSamples down to numControls; the same sample modes plus a bilinear polynomial; accessors, control data, `GetPosition` inside and outside [0,1]^2) | exact | pass |
| | `BSplineSurfaceFit.invalid` (throw parity for the four preconditions) | exact | 2000 of 2000 throw on both sides |
| `BSplineSurface.h` | `BSplineSurface.evaluate` (`N` = 1, 2, 3; open uniform, open nonuniform, repeated interior knot, periodic uniform, periodic nonuniform bases; controls present or deferred; `SetControl`/`GetControl` in and out of range; `GetControls`, `GetBasisFunction`, `GetNumControls`; the `ParametricSurface` accessors, `Evaluate` for orders 0, 1, 2, 3 and >= `SUP_ORDER`, `GetPosition`, `GetUTangent`, `GetVTangent`) | exact | pass |
| | `BSplineSurface.evaluate.invalidOrder` (orders 4, 5: `BasisFunction` asserts) | exact | 2000 of 2000 throw on both sides |
| `BSplineVolume.h` | `BSplineVolume.evaluate` (the same over three bases; `operator bool`, domains, `Evaluate` for orders 0..3 and >= `SUP_ORDER`) | exact | pass |
| | `BSplineVolume.evaluate.invalidOrder` (orders 4..9) | exact | 2000 of 2000 throw on both sides |

## Port defects fixed

**`src/QFNumber.ts`, `scalarDiv` (the port of `operator/(T s, QFNumber q)`)
negated after dividing.** Upstream writes

```cpp
auto x1 = -(s * q.x[1]) / denom;
```

and unary minus binds tighter than `/`, so the numerator is negated *before*
the division. The port computed `-((q.x[1] * s) / denom)`. For `N = 1` the two
are the same computation (IEEE division is sign symmetric). For `N >= 2` the
division is itself a quadratic-field division whose numerator coefficients are
differences `a*b - c*d*d'`; when the two products are equal, `(-a)*b -
(-c)*d*d'` is `+0` while `-(a*b - c*d*d')` is `-0`, so the port flipped the
sign of an exact zero coefficient of `x[1]`. Error size: the sign bit of a
zero coefficient; the value is otherwise identical (nonzero coefficients are
exact negations either way). Reach before the fix: 626 of 2000 deep-run
records of `QFNumber.arithmetic.n2` and 821 of 2000 of `.n3` (7 of 20 golden
records each). Fixed by negating the product before the division, which is
upstream's expression; regression test in `test/QFNumber.test.ts`
("negates the numerator before dividing, as upstream does") on
`q = ((1, 0, 1), (0, 1, 1), 0)`, `s = 1`, which also shows that the two
groupings differ on that input.

Nothing else in the eight headers differs from upstream on any of the 70 000
deep-run records: the bigint replacement of `UIntegerAP32` (`ShiftRightToOdd`,
`RoundUp`, `GetPrefix`, the left-aligned `operator<`), the round-to-nearest
conversions to `double` and `float`, `Convert`, the `std::nextafter`
emulation, the banded Cholesky solves and accumulation orders of both fits and
the tensor-product evaluation of surfaces and volumes are bit-identical.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `BSNumber.construct.string.singleChar` | [#95](https://github.com/gradientspaceai/gtengine-js/issues/95) | 2000 of 2000 records. `ConvertToInteger` validates characters only when the string has more than one, so upstream accepts a one-character non-digit (optionally signed) and returns `digit - '0'` as its value (`"x"` is 72); the port asserts. The main string case excludes exactly these strings and agrees on every record, including all its throw records. |

Preserved upstream defects compared bit for bit:
[#366](https://github.com/gradientspaceai/gtengine-js/issues/366)
(`BSPrecision::operator+` pairs one operand's `maxExponent` with the other's
`minExponent`; `BSPrecision.operators` draws the finding's operands
`IS_DOUBLE + IS_FLOAT` and `(-17,20,24) + (-5,20,12)` in both orders) and
[#50](https://github.com/gradientspaceai/gtengine-js/issues/50)
(`SWInterval.overflow`, and the huge and subnormal modes of the other
interval cases).

## Independent-reference checks (agreement is not correctness)

Run once over the 2000-record deep-run outputs of the C++ build (which the port
reproduces bit for bit), outside the harness:

| check | result |
| --- | --- |
| Every `BSNumber` output of `construct.double/float/integer/string`, `conversions.boundary`, `arithmetic`, `accessors`, `convert` and the exact members of `std.exact` (`fabs`, `frexp` incl. its exponent, `ldexp`, `FMA`, `RobustSOP`, `RobustDOP`) against exact dyadic arithmetic over `BigInt`: canonical sign, biased exponent, bit count and words, and both conversions against an independent round-to-nearest-even to binary64 and binary32 | 0 disagreements. `Convert` is checked against an independent model of all four rounding modes and its throw conditions |
| `BSNumber.compare` against the sign of the exact difference | 0 disagreements |
| String validity: throw iff the string does not match `[+-]?(0|[1-9][0-9]*)` | 0 disagreements |
| `SWInterval` enclosure of the exact result (BigInt rationals; divisions by cross multiplication) for the leaf operations and `u+v`, `u-v`, `u*v`, `u/v` over all four endpoint combinations | every violation is finding #50: 53 leaf products and 230 leaf quotients that overflow, 140 interval products, 1 sum, 1 difference, 81 interval quotients. Six of the quotients have no `DBL_MAX` in the output at all: dividing `[-6.4e-323, 1.5e-323]` by `[2e-323, 3e-323]` returns `[-1.15e-14, 2.66e-15]` instead of an enclosure of `[-3.2, 0.75]`, because `Reciprocal` overflows to `[+inf, +inf]`, `nextafter` pulls both ends back to `DBL_MAX`, and the product with the subnormal numerator is then tiny. This is #50's mechanism, but the error is 14 orders of magnitude, not an ulp; proposed as an addition to the #50 entry below |
| `SWInterval.nextafter` against the adjacent binary64 values computed from the bit pattern | 0 disagreements over 1851 non-NaN records (NaN stays NaN on 149) |
| `QFNumber.compare.n1` against the exact sign of `(a0-b0) + (a1-b1)*sqrt(d)` | 0 disagreements, including the 327 non-lattice records where upstream's `diff` is computed in rounded doubles |
| `BSplineCurveFit`: interior control points against a dense least-squares solve (independent Cox-de Boor basis, normal equations by Gauss-Jordan), end control points against the first and last samples | worst scaled difference 1.7e-13; ends bit-identical to the samples on every record |
| `BSplineSurfaceFit`: every control point against `X0 P X1^T` with dense `X_d = (A^T A)^{-1} A^T` | worst scaled difference 3.6e-14 |
| `BSplineSurface` position and first derivatives, `BSplineVolume` position, against an independent Cox-de Boor evaluation with upstream's periodic wrapping and open clamping | worst 8.9e-17 (13 101 surface components) and 4.4e-17 (10 567 volume components) |

## Not covered

| header / entry point | reason |
| --- | --- |
| `BSNumber` over `UIntegerFP32<N>` | Fixed-precision storage is not ported (`porting-status.json`: `UIntegerFP32.h` omitted); the value semantics are the same as `UIntegerAP32`'s except for the capacity checks, which the port drops. |
| `BSNumber::Write`, `Read` | Serialize the C++ word layout; not ported (no word layout in the bigint port). |
| `BSNumber::IsValid` | Compiled only under `GTE_VALIDATE_BSNUMBER`, which the reference build does not define; the port's `isValid()` is covered by `test/BSNumber.test.ts`. |
| `BSNumber(char const*)` | Delegates to `BSNumber(std::string)`, which is covered. |
| `BSNumber(int32_t)` with `INT32_MIN`, `BSNumber(int64_t)` with `INT64_MIN` | Upstream evaluates `-number`, which is signed overflow (undefined behaviour) for these values; excluded from the generator (see the suspects). MSVC happens to produce the right value. |
| `BSNumber::SetSign(+-1)` on zero | Creates an invalid number whose conversion reads `UInteger` words that do not exist (out-of-range read); cannot produce a golden record. |
| `GTE_THROW_ON_CONVERT_FROM_INFINITY_OR_NAN`, `GTE_ASSERT_ON_QFNUMBER_MISMATCHED_D`, `GTE_THROW_ON_INVALID_SWINTERVAL` | Configuration macros the reference build does not define; the default (graceful) paths are covered, including mismatched `d` values in the QFNumber arithmetic case. |
| `QFNumber<BSRational, N>` and other non-`double` `T` | The port implements the number-based path only (`src/QFNumber.ts` header). |
| `QFNumber<T, N>()` and `QFNumber<T, N>(d)` for `N >= 2` | The port's coefficients carry their depth at run time; these constructors have no port equivalent beyond the coefficient constructor, which is covered. |
| `SWInterval<float>` | The port maps C++ floating point to binary64 only. |
| `BSRational` | Belongs to verify group 6 (it is used here only indirectly, nowhere in a compared path). |

## Upstream bug suspects

New ones only.

**1. `Convert(BSNumber const& input, ..., BSNumber& output)` is not alias
safe (result-corrupting when aliased; minor/doc otherwise).** The function
calls `output.GetUInteger().SetNumBits(precision)` and
`SetAllBitsToZero()` before it reads `input.GetUInteger()`. When `output` and
`input` are the same object, the input's bits are zeroed first. Measured on
the MSVC build: with `x = BSN(12345.678) * BSN(3.25)`,
`Convert(x, 10, FE_TONEAREST, y)` gives `40128`, but
`Convert(x, 10, FE_TONEAREST, x)` gives `3.637978807091713e-12` (bit count 1).
No upstream caller aliases this overload: `APConversion.h` writes
`Convert(aMax, 2 * mPrecision, FE_UPWARD, aMax)` eight times, but on
`BSRational`, whose overload converts into a local `BSNumber` first and is
alias safe. The hazard is undocumented and the sibling overload invites the
pattern. *Port:* not applicable (`convertBSNumber` returns a new object).

**2. `BSNumber(int32_t)` and `BSNumber(int64_t)` evaluate `-number` for the
most negative value**, which is signed overflow (undefined behaviour). MSVC
x64 happens to produce the correct magnitude (`BSN(INT32_MIN)` is
`-2^31`, `BSN(INT64_MIN)` is `-2^63`), but an optimizer is entitled not to.
Suggested fix: `magnitude = 0u - static_cast<uint32_t>(number)` (and the
64-bit analogue). *Port:* not affected (`fromNumber` / `fromBigInt`).
