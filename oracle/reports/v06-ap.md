# Verify group 6 (arbitrary precision) against the MSVC build of upstream GTE

Family `v06-ap`, 41 cases, 20 golden records each. The group is
alphabetical rather than thematic: `BSRational.h`, `APConversion.h`,
`APInterval.h`, `BSPPolygon2.h`, `BSplineCurve.h`, `BSplineGeodesic.h`,
`BSplineReduction.h` and `ArbitraryPrecision.h` (an umbrella header with no
code).

Deep run `npm run oracle:deep -- 2000 v06-ap`: 82 000 records (16 321 of
them throw records, all with throw parity), **all 42 tests pass** (41 cases
plus the claim check). 36 cases are `{ exact: true }`, one is compared with
the default tolerance (`BSRational.std.libm`, C math library) and four are
declared deviations (#95, #168 twice, #280). Outside the deviation cases the
deep run compares 2 213 849 floating-point outputs, of which 2 211 357 are
bit-identical; the 2492 others all belong to `BSRational.std.libm` (max
scaled error 2.22e-16). Every discrete output (signs, exponents, bit counts,
32-bit words, digests, iteration counts, polygon indices, BSP trees, point
locations, booleans) matches exactly.

**Four port defects were found and fixed** (the `BSRational` scalar pair
constructors, `BSRational.fromString` with two signs, the dot-product seed of
`RiemannianGeodesic`, and the incomplete #280 fix of
`APConversion.estimateAmB`; see below).

## How BSRational is compared

The port replaces upstream's `UInteger` layer with `bigint` (v05). The C++
side instantiates `BSRational<UIntegerAP32>`; `OutBSR` emits the numerator and
the denominator with v05's `OutBSN` (sign, biased exponent, exponent, bit
count, the 32-bit words of the odd integer, or above 8 words the top three,
the bottom two and a 48-bit FNV-1a digest, then `operator double` and
`operator float` of that `BSNumber`), followed by the rational's own
`operator double` and `operator float`. Every `BSRational` output therefore
also tests both conversions (Convert to 53/24 bits, then the `BSNumber`
conversion). An `APInterval` endpoint is emitted the same way, except an
infinity sentinel of finding #280 (`SetSign(+-2)` on a zero rational, possibly
multiplied by a finite number): its sign, its numerator's sign and biased
exponent and its denominator; not its numerator's bit count (upstream's
`UInteger::Mul` of a zero-bit operand leaves `numBits - 1` zero words, the
port a zero `bigint`) and not its conversions (they would read words that do
not exist). `BSPPolygon2.h` is included with `private` opened up, so the edge
array (which can be longer than the edge map, #169) and the whole BSP tree are
compared; the replay reads the port's private fields the same way.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `BSRational.h` | `BSRational.construct.double` (every finite double, +-inf and NaN through `BSNumber`'s graceful exit, pairs with negative, infinite and zero denominators), `.construct.float` (binary32 subnormals, `FLT_MAX`, inf, NaN), `.construct.integer` (`int32_t`, `uint32_t`, `int64_t`, `uint64_t`, single and pairs: small, powers of two, 2^k +- 1, 2^53 + 1, the extremes, full-width random, zero denominators), `.construct.bsnumber` (`BSRational(BSNumber)`, `(BSNumber, BSNumber)` with a nonzero numerator), `.construct.string` (signs, `x`, `x.y`, `.y`, `x.`, 1..60 digits on each side, leading and trailing zeros in the fraction; throw parity for `""`, `"+"`, `"-"`, exponents, a second `.`, leading zeros in `x`, non-digits, and a second sign) | exact | pass; 501 / 499 / 134 / 130 / 649 of 2000 throw on both sides |
| | `BSRational.construct.bsnumber.zeroNumerator`, `.construct.string.zero`, `.construct.string.singleChar` | deviation (#168, #168, #95) | 2000 of 2000 records deviate in each |
| | `BSRational.compare` (`== != < <= > >=`, both orders, on independent pairs, equal values with equal and with different representations `pk/qk` vs `p/q`, negations, dyadic and non-dyadic neighbours) | exact | pass; 856 / 590 / 554 records with x < y / x = y / x > y |
| | `BSRational.arithmetic` (unary `+ -`, `+ - * /` both orders, `+= -= *= /=`, `Negate`, `SetSign`, `GetSign`), `.divide.zero` (throw parity for `/` and `/=` by zeros built four ways) | exact | pass; `divide.zero` 2000 of 2000 throw on both sides |
| | `BSRational.convert` (`Convert` to `BSNumber` and to `BSRational` at precision 53, 24, 1..300 and <= 0, `Convert` to `double` and `float`, all four `<cfenv>` modes and an unsupported mode; inputs on and next to binary64/binary32 ties, subnormal grids), `.convert.subnormalAndOverflow` (the same conversions for results in the binary64/binary32 subnormal ranges and beyond the largest finite values) | exact | pass; 427 throw records (precision <= 0, unsupported mode) |
| | `BSRational.std.exact` (`std::fabs frexp ldexp floor ceil sqrt fmod remainder`, `gte::clamp invsqrt isign saturate sign sqr`, `FMA`, `RobustSOP`, `RobustDOP`), `.std.remainder` (quotients up to 2^1060, zero and infinite divisors) | exact | pass |
| | `BSRational.std.libm` (`std::acos acosh asin asinh atan atanh atan2 cos cosh exp exp2 log log2 log10 pow sin sinh tan tanh`, `gte::atandivpi atan2divpi cospi exp10 sinpi` on non-dyadic arguments) | tolerance 1e-12 (C math library; the port's `exp2` is `Math.pow(2, x)`) | pass; 45 508 of 48 000 bit-identical, max 2.22e-16 |
| `APConversion.h` | `APConversion.estimateSqrt` (both overloads; exact squares, non-squares, exponents to +-300, zero, a negative input: throw parity), `.estimateApB`, `.estimateAmB` (both f'' branches, early returns, exhausted bisections whose bracket stays valid, Newton exhaustion; the LogError on `aSqr = bSqr = 0`), `.estimateAmB.aLessThanB` (the violated precondition of #280 item 2), `.estimate` (both overloads over `QFNumber<BSRational, 1>`: d = 0 or y = 0, y > 0, y < 0, d < 0 throws), `.accessors` (constructor, `SetPrecision`, `SetMaxIterations`, getters, invalid arguments) | exact | pass; estimateAmB branch split: f'' > 0 795 converged / 369 bisection exhausted, f'' < 0 539 / 123, 174 throw; estimate: 931 / 410 / 402 / 257 |
| | `APConversion.estimateAmB.bisectionExhausted` | deviation (#280) | 2000 of 2000 records deviate |
| `APInterval.h` (over `BSRational`) | `APInterval.leaf` (`Add Sub Mul Div(u, v)`, division by zero returns `Reals`), `.internal` (four-argument `Add Sub Mul Mul2 Div`, `Reciprocal`, `ReciprocalDown`, `ReciprocalUp`, `Reals`), `.operators` (every non-class operator on interval/interval, interval/scalar, scalar/interval, compound assignments, `operator[]`, `GetEndpoints`, all constructors; divisors with a zero endpoint produce sentinels, whose negation is emitted too) | exact | pass; all nine sign branches of `operator*` reached (71 to 471 records each); 873 divisors without zero, 873 with a zero endpoint or straddling zero |
| | `APInterval.divide.zeroInterval` (#280: `[0, 0]` asserts instead of returning `Reals`) | exact | 2000 of 2000 throw on both sides |
| `BSPPolygon2.h` | `BSPPolygon2.negation`, `.intersection`, `.union`, `.difference`, `.exclusiveOr` (lattice rectangles, triangles, L, U, square with hole, star polygons, diamonds, perturbed copies; the second polygon independent, identical, translated by a lattice offset (shared edges and vertices), reversed or scaled about a vertex (nested); epsilon 0, 1e-10, 0.125; the whole result state, BSP tree and eight point locations), `.construct` (`InsertVertex`/`InsertEdge` return values, the state before and after `Finalize`, `PointLocation`, copy construction and assignment, `GetEdge`), `.splitEdge.unmapped` (#169: every record leaves an edge unmapped), `.invalid` (throw parity) | exact | pass; empty results throw on 713 / 352 / 551 / 1102 of 2000 records (intersection / union / difference / xor, #169); `Finalize` splits edges on 980 construct records; SplitEdge's "Edge does not exist" fires on 1 |
| `BSplineCurve.h` | `BSplineCurve.evaluate` (`N` = 1..4, degree 1..5, open uniform, open nonuniform with a repeated interior knot, periodic uniform and nonuniform; controls present or deferred; `SetControl`/`GetControl` in and out of range; the `BasisFunction` and `ParametricCurve` accessors; `Evaluate` for orders 0..3 and >= `SUP_ORDER`; `GetPosition`, `GetTangent`, `GetSpeed`, `GetLength`, `GetTotalLength`, `GetTime`) | exact | pass |
| `BSplineReduction.h` | `BSplineReduction.reduce` (`N` = 1..3, 2..12 controls, degree 1..4, fractions uniform, k/n, 1, 1.5, 0, -0.5 and just below k/n; the object reused), `.invalid` (throw parity) | exact | pass; 406 reduced, 884 clamped to degree + 1, 710 copies |
| `BSplineGeodesic.h` | `BSplineGeodesic.computeGeodesic` (with the refine callback and the progress accessors, then `ComputeTotalLength` and `ComputeTotalCurvature`), `.segment` (`ComputeSegmentLength`, `ComputeSegmentCurvature`, `Subdivide`, `Refine`; periodic bases included), `.degenerate` (throw parity: a zero-length segment, zero controls) | exact | pass; paths of 2, 3, 5, 9 points; `Refine` moved the midpoint on 1136 of 1956 records |

`BSplineGeodesic` is arithmetic only: the metric and the Christoffel symbols
come from `BSplineSurface` (whose periodic wrap calls `std::fmod`, exact), and
`RiemannianGeodesic` adds `std::sqrt` and the Gaussian-elimination `Inverse`
of `GMatrix`. No C math library call decides its control flow, so unlike
`EllipsoidGeodesic` (v18) the whole algorithm, including the argmin of
`Refine`'s line search, is compared bit for bit on the real class. The tuning
parameters are drawn small (`integralSamples` 2..8, `searchSamples` 1..6,
`subdivisions` 0..3, `refinements` 0..2); as upstream, the derived steps stay
those of the defaults (the port's `updateDerivedParameters` is not called).

## Port defects fixed

**1. `src/BSRational.ts`: the scalar pair constructors normalized the
denominator.** Upstream's `(float, float)`, `(double, double)`,
`(int32_t, int32_t)`, `(uint32_t, uint32_t)`, `(int64_t, int64_t)` and
`(uint64_t, uint64_t)` constructors only move the sign of a negative
denominator to the numerator; only the `(BSNumber, BSNumber)` constructor
moves the denominator's exponent into the numerator. The port routed
`fromNumber(n, d)`, `fromFloat32(n, d)` and `fromBigInt(n, d)` through
`fromBSNumber`, so `fromNumber(3, 4)` stored 3*2^-2 / 1 where upstream stores
3 / 4 (denominator exponent 2). Visible through `getNumerator` /
`getDenominator` only: the value, the comparisons, the conversions and every
arithmetic result are identical (each operator ends in the normalizing
constructor, which cancels the difference). Reach before the fix: 6, 11 and 3
of the 20 golden records of `construct.double`, `.float`, `.integer`, and
every case whose generator used an integer pair (`arithmetic`, `std.exact`,
`APConversion.estimate`, `APInterval.operators`). Fixed with a
`fromScalarPair` body that does exactly what upstream does; regression test in
`test/BSRational.test.ts` ("normalizes the denominator exponent to zero in
fromBSNumber only"), and the port's own test helper no longer asserts the
normalization for these constructors.

**2. `src/BSRational.ts`: `fromString` accepted a second sign.** Upstream
strips one sign and hands the integer part to `BSNumber::ConvertToInteger`,
which has no sign handling: `"--5"` leaves `"-5"`, a multi-character string
that does not start with a nonzero digit, and asserts. The port called
`BSNumber.fromString`, which accepts a sign, and returned -5 for `"--5"`
(likewise `"+-2.5"`, `"-+12"`, `"++0.5"`, `"-+0"`, `"+-37.25"`): a value where
upstream throws. Found on 8 of 200 records of a first deep run (added to the
generator's invalid forms). Fixed by rejecting a leading sign of the integer
part before `BSNumber.fromString` (a lone sign character such as `"+-"` is a
one-character integer part, finding #95, on which the port already asserts);
regression test "rejects a second sign as upstream does".

**3. `src/RiemannianGeodesic.ts`: the dot-product seed.** Upstream's vectors
are `GVector`s, whose `Dot` accumulates from the literal 0; the port built
`diff` (in `computeSegmentLength` and `computeSegmentCurvature`) and `acc` (in
`computeIntegrand`) with the base `Vector` functions `sub` and `mul`, so
`dot` used `Vector.h`'s seed, the first product. On a surface whose
Christoffel symbols vanish along the segment and a direction with negative
components, `qForm1` became -0 instead of +0, the acceleration (-0, -0), and
the curvature `sqrt(-0) = -0` where upstream returns +0: the sign of a zero
curvature, deep-run record 129 of `computeGeodesic` (1 of 2000). Fixed by
wrapping both vectors as `GVector`s; regression test in
`test/BSplineGeodesic.test.ts` ("returns +0 curvature on a plane for any
segment direction"), which fails on the old source. The file belongs to v18;
`EllipsoidGeodesic` shares it and its tests still pass.

**4. `src/APConversion.ts`: the #280 fix of `estimateAmB` was incomplete.**
Finding #280 item 1 says that when a bisection loop of `EstimateAmB` runs out
of iterations after updating the bisected endpoint in its rounding branch,
upstream feeds Newton's method a stale square, and that "all 145 [failing
brackets] become valid brackets once the squares are recomputed". The port
recomputed the square. The oracle shows that this is not the defect's cause:

* Substituting upstream's expression (removing the recompute) makes the port
  bit-identical to C++ on all 2000 deep-run records of the deviation case,
  which proves the recompute was the only difference.
* On those 2000 records (selected because upstream's bracket misses a - b,
  checked exactly), the recompute fix still returned a bracket missing a - b
  on 755.
* On 3 of 1826 records of the main case upstream returns a wrong bracket
  although no square is stale (the bisection only moved the other endpoint):
  aSqr = 701408733, bSqr = 102334155, precision 100, maxIterations 1 returns
  [t, t] with t = 1.99805758446352...*2^13 > a - b.

The cause is that after an exhausted bisection Newton's method starts outside
the basin in which it converges monotonically; the bisection bracket itself
is always valid (every step keeps tMin <= a - b <= tMax and rounds outward).
The port now evaluates upstream's expressions unchanged (the recompute is
removed) and replaces the result by the bisection bracket only when a
bisection ran out of iterations, aSqr >= bSqr (item 2, the precondition, stays
preserved), and the returned bracket provably misses a - b (exact rational
comparisons, no square roots). Result on the deep run: all 2000 deviation
records and all 1826 main records return valid brackets; the port differs
from upstream on exactly the records where upstream is wrong (3 main records
are redrawn by the generator's exact probe), and agrees on the 6 main records
whose bisection ran out but whose Newton bracket happens to be valid. Size of
the upstream error: the bracket misses a - b by up to its own width (both
endpoints on the same side). Regression test in `test/APConversion.test.ts`
("returns a valid bracket whenever the bisection runs out of iterations"):
on all six inputs both upstream's expression and the previous port fix return
a bracket that misses a - b (the first three have no stale square, the last
three were among the 755).

Nothing else in the seven headers differs from upstream on any of the 82 000
deep-run records: the bigint replacement of the `UInteger` arithmetic under
`BSRational` (constructors, cross-multiplied comparisons, the four operators,
`Convert` in every mode, `frexp`, `ldexp`, the string parser), the interval
endpoint selection including the sentinel products, the BSP construction with
its edge splits and `std::map` orders, the B-spline evaluation, the Romberg
Gram matrices and banded inverse of the reduction, and the geodesic
subdivision and refinement are bit-identical.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `BSRational.construct.bsnumber.zeroNumerator` | [#168](https://github.com/gradientspaceai/gtengine-js/issues/168) item 1 | 2000 of 2000 records. `BSRational(BSN(+-0), BSN(d))` with d of nonzero exponent: upstream's numerator keeps sign 0 but gets the biased exponent `-exponent(d)`; the port keeps the canonical zero. The case first emits `r == 0`, `0 == r`, `r < 1`, `-1 < r` and `r + d`, which agree on every record: the invalid encoding is only visible through `GetNumerator`. |
| `BSRational.construct.string.zero` | #168 item 2 | 2000 of 2000 records. Every zero-valued string (`"0"`, `"+0.000"`, `".0"`, `"0."`, `"-0.0"`, ...) gets a numerator of sign +1 or -1 and no bits upstream; the port keeps sign 0. Only the signs and the bit count are emitted (converting or comparing upstream's number would read words that do not exist). |
| `BSRational.construct.string.singleChar` | [#95](https://github.com/gradientspaceai/gtengine-js/issues/95), inherited from `BSNumber::ConvertToInteger` | 2000 of 2000 records: `"a.5"` is 49.5 upstream, the port asserts. |
| `APConversion.estimateAmB.bisectionExhausted` | [#280](https://github.com/gradientspaceai/gtengine-js/issues/280) item 1 (port fix corrected, above) | 2000 of 2000 records: inputs near the ratio (7 + 3*sqrt(5))/2 with 1..4 iterations, kept only when upstream's bracket misses a - b (the observable symptom, checked exactly on the C++ side). Upstream's bracket is wrong on all of them by construction; the port's is right on all of them (checked exactly). |

Preserved upstream behaviour compared bit for bit: #280 item 2 (`EstimateAmB`
with aSqr < bSqr: 1904 of 2000 brackets inverted, `[0, aMax - bMin]`), #280
item 3 (`EstimateSqrt(0)` = `[0, 2^-1074]`, 382 records), the #280 `APInterval`
items (sentinel arithmetic, `[0, 0]` divisor asserts), #169 (empty Boolean
results throw; `SplitEdge` leaves an edge unmapped on all 2000
`splitEdge.unmapped` records and on 3 `construct` records; the Romberg Gram
error of `BSplineReduction`), and the two-step rounding of `BSRational`
conversions (new suspect 1 below).

## Independent-reference checks (agreement is not correctness)

Run once over the 2000-record deep-run outputs (the port reproduces them bit
for bit, so the checks were run on the port's objects), outside the harness:

| check | result |
| --- | --- |
| Every `BSRational` built by `construct.double`, `construct.string` and `arithmetic` (x, -x, x+y, x-y, x*y, x/y), and `fabs`, `frexp` (value and range [1/2,1)), `ldexp`, `FMA`, `RobustSOP`, `RobustDOP`, `floor`, `ceil` of `std.exact`, against exact BigInt rationals built independently from the recorded doubles | 0 value mismatches |
| `BSRational.compare` against the sign of the exact difference | 0 mismatches |
| String validity: throw iff the string does not match `[+-]?((0\|[1-9][0-9]*)(\.[0-9]*)?\|\.[0-9]+)` | 0 disagreements |
| `Convert(x, p, mode, BSNumber&)` against an independent rounding of the exact rational to p bits in each mode | 0 mismatches (1548 records) |
| `operator double` / `operator float` and `Convert(x, mode, double&/float&)` against correct IEEE rounding (subnormals, overflow) | every difference is explained by upstream's two-step rounding (new suspect 1): `operator double` is 1 ulp off in the subnormal range on 18 of 2000 `arithmetic` operands, 2..12 of 2000 in the other checks; in `convert`, 16 / 21 / 20 double and 78 / 63 / 86 float results of the down / toward-zero / up modes are on the wrong side of x |
| `APConversion`: brackets against exact comparisons of the endpoints with sqrt(aSqr), sqrt(aSqr) + sqrt(bSqr), sqrt(aSqr) - sqrt(bSqr) (squared forms, no square roots), and the width < 2^-precision when the iterations sufficed | `estimateSqrt`: 0 violations except `EstimateSqrt(0)` = `[0, 2^-1074]` (382 records, #280 item 3); `estimateApB`: 0 violations; `estimateAmB`: 0 violations after the fix (3 before, upstream); `aLessThanB`: 1904 inverted brackets (#280 item 2) |
| `APInterval.operators` against exact endpoint arithmetic: `u+v`, `u-v`, `u*v` and `s*u` as [min, max] of the exact endpoint combinations, `u/v` for 0 not in v; for divisors with a zero endpoint, every sentinel on the left must be -inf and on the right +inf, and no finite endpoint may cut into the true range u * [1/b, +inf) or u * (-inf, 1/a] | 0 violations |
| `BSplineCurve` positions and derivatives of orders 1..3 against an independent de Boor evaluation (knots expanded from the unique knots, periodic wrap, clamping, control-point difference formula) | worst scaled difference 4.4e-16 over 16 861 jet entries |
| `BSplineReduction` against the exact L2 projection (Gram matrices integrated exactly span by span with 6-point Gauss-Legendre, dense solve) | worst difference 0.137 of the control scale; >= 1e-2 on 486 of 1290 reduced records (degree 2..4 mostly): the Romberg error of #169, much larger in the controls than the "~1e-3 entry error" the finding quotes |
| `BSPPolygon2` Boolean results: point locations of off-boundary query points (distance > 1e-6 + 2*epsilon from every input edge, both inputs simple) against the winding-number predicate combined by the operation | negation 0 of 5897 wrong; intersection 68 of 3729, union 33 of 4888, difference 84 of 4236, exclusive or 83 of 2583 wrong (new suspect 2) |
| `BSplineGeodesic` | not checked independently beyond the zero-curvature regression; `RiemannianGeodesic` was checked in v18 |

## Not covered

| header / entry point | reason |
| --- | --- |
| `ArbitraryPrecision.h` | An umbrella header of includes; no computational entry point (the port's module re-exports). |
| `BSRational` over `UIntegerFP32<N>` | Fixed-precision storage is not ported (`UIntegerFP32.h` omitted); the "maximum precision exceeded" LogError of `Convert` cannot fire for `UIntegerAP32`. |
| `BSRational::Write`, `Read` | Serialize the C++ word layout; not ported. |
| `GTE_BINARY_SCIENTIFIC_SHOW_DOUBLE` (`mValue`) | Debug macro, not defined in the reference build. |
| `BSRational(int32_t)` / `(int64_t)` with `INT32_MIN` / `INT64_MIN` | `BSNumber` evaluates `-number` (signed overflow, v05 suspect 2); excluded. |
| `BSRational(char const*)` | Delegates to `BSRational(std::string)`, covered. |
| `BSRational::SetSign(+-1)` on a zero rational | Creates an invalid number whose conversion reads integer words that do not exist; `SetSign` is covered on nonzero values. |
| `APInterval<BSNumber>`, `APConversion<Rational>` for other types | The port implements the `BSRational` instantiation only (header comments); `QFN2` is a dead alias. |
| `APInterval` arithmetic with a sentinel as the left operand of a product, or in a sum | Upstream's `UIntegerALU32::Mul` reads word 0 of the sentinel's empty word array (undefined behaviour; an access violation, not an exception), so no golden record can be produced; `Reals() * s`, `Reals() + v` and `(u / [0, b]) / v` are therefore not generated (suspect 3). |
| `GTE_THROW_ON_INVALID_APINTERVAL` | Not defined in the reference build; the port's `throwOnInvalid` flag stays false. |
| `BSPPolygon2::Print` | Debug macro. |
| `BSplineCurve` of degree 0 or other `BasisFunction` asserts | `BasisFunction` preconditions belong to its own group. |
| `BSplineGeodesic` with the default tuning (7 subdivisions, 8 refinements, 32 search samples) | About 10^6 surface evaluations per record; the same code paths run with the small parameters above. |
| `float` instantiations | The port maps C++ floating point to binary64 only. |

## Upstream bug suspects

New ones, and corrections to existing entries.

**1. `BSRational` conversions to `double` and `float` round twice
(result-corrupting for the directed modes; minor for round-to-nearest).**
`Convert(BSRational, mode, FPType&)` (and `operator double` / `operator
float`, which call it with `FE_TONEAREST`) first converts to a `BSNumber` of
53 (24) bits with the requested mode and then applies `BSNumber`'s own
conversion, which always rounds to nearest onto the subnormal grid and to
infinity beyond the largest finite value. So (a) `operator double` is 1 ulp
off in the binary64 subnormal range whenever the 53-bit rounding lands on a
tie of the subnormal grid (x ~ 1.0773*2^-1055 gives `0x...89e36` instead of
`0x...89e37`; 18 of 2000 `arithmetic` operands); (b) the directed modes
return results on the wrong side of x there: `FE_UPWARD` of a positive
x < 2^-1075 returns 0, `FE_DOWNWARD` of a negative one returns -0,
`FE_TOWARDZERO` of 1.99999809*2^-145 to `float` returns 2^-144 * 1.0 >
x; (c) beyond the largest finite value `FE_DOWNWARD` and `FE_TOWARDZERO`
return +inf for positive x (and `FE_UPWARD`, `FE_TOWARDZERO` return -inf for
negative x) instead of the largest finite value. Measured on the targeted
case `BSRational.convert.subnormalAndOverflow` (2000 records): 14 to 115
wrong-side results per mode and range. The only upstream caller of the
directed double conversions, `APConversion::GetMinOfSqrt/GetMaxOfSqrt`,
converts values in [1/2, 2), where the conversions are correct. *Port:*
preserved (bit-identical); a fix would round once, at the subnormal quantum,
with the requested mode, and saturate per mode at overflow.

**2. `BSPPolygon2` Boolean results misclassify points far from the
polygons (result-corrupting).** With `Real = double` the intersection point
of an edge with a splitting line is computed separately along each split
edge (`v0 + t*(v1 - v0)`), and the same geometric point comes out as
different doubles along different edges; `InsertVertex` matches vertices
exactly, so the result polygon gets near-duplicate vertices joined by
sub-ulp edges. When `Finalize` builds the result's BSP tree, such an edge's
line has an essentially arbitrary direction and splits the plane wrongly.
Example (epsilon 0): P = diamond (-3,-3), (0,0), (-3,3), (-6,0); Q = hexagon
(-9,-6), (-3,-2), (-1,2), (-7,6), (-9,6), (-7,4). `P & Q` has the vertices
(-3.6,-2.4), (-3.5999999999999996,-2.4) and (-3.5999999999999996,
-2.4000000000000004) and the 4e-16-long edges between them, and
`PointLocation` reports (100, 0) and (7.98, -1.83) inside the intersection.
Over the deep run 0.7 % to 3.2 % of off-boundary query points get the wrong
side in the results of `&`, `|`, `-` and `^` (68 / 33 / 84 / 83 points);
negation, which does not split edges, is always right. Some of the failures
involve coincident edges (P with itself, P with its reversal) and were not
root-caused individually. *Port:* preserved (bit-identical). An exact
`Real` (BSRational) would remove the mechanism shown.

**3. `APInterval` crashes on its own infinite intervals (undefined
behaviour; extends the #280 `APInterval` entry).** The `SetSign(+-2)`
sentinels have an empty word array, and `UIntegerALU32::Mul` reads word 0 of
its left operand unconditionally, so any product with a sentinel as the left
operand reads out of bounds: `APInterval<BSRational>::Reals() * BSRational(2)`
is an access violation (segmentation fault) on the MSVC build, as is
dividing an interval with an infinite endpoint by an interval that does not
straddle zero. A sentinel as
the right operand (`u / [0, b]`, `s / [a, 0]`) is harmless and produces a
product whose numerator keeps `numBits(u) - 1` zero words. *Port:* not
affected (bigint arithmetic); `isInfinite` detects the sentinels.

**4. Correction to #280 item 1 (`EstimateAmB`).** The stale square is not the
cause: after an exhausted bisection Newton's method starts outside its basin,
and the bracket can miss a - b with or without a stale square (3 main-case
records with none; 755 of 2000 records still wrong after recomputing the
square). Returning the bisection bracket when it ran out of iterations is
always valid. See "Port defects fixed" 4.

**5. Corrections to the #168 entry (`BSRational.h`).** Item 1: the invalid
zero is not reachable from `BSRational(0.0, 1024.0)` (the scalar pair
constructors do not touch the exponents) nor from arithmetic (every operator
returns `BSRational(0)` for a zero numerator); it arises only from
`BSRational(BSNumber, BSNumber)` with a zero numerator. It compares equal to
zero (`operator==` tests the signs first) and does not disturb arithmetic,
only `GetNumerator()` shows it. Item 2 concerns every zero-valued string,
not only `"-0.0"`: `"0"`, `"0.0"`, `".0"`, `"+0"` get sign +1 (a "positive
zero" that compares unequal to zero and whose conversion reads integer
words that do not exist).

**6. Addition to the #169 `BSplineReduction` entry (magnitude).** Against the
exact L2 projection the reduced control points are off by up to 0.137 of the
control-point scale, >= 1e-2 on 486 of 1290 reduced records (degree >= 2
mostly); the ~1e-3 Gram-entry error is amplified by the inverse.
