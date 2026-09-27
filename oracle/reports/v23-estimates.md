# Group 23 (`v23-estimates`) — C++ oracle report

52 cases, 20 records each in `oracle/golden/v23-estimates.txt` (1040 records,
29 of them C++ exceptions from the two throw-parity cases). 41 cases are
declared `{ exact: true }`, 10 carry the default `1e-12` tolerance because
their upstream path calls `std::sin` or `std::acos` (with the arithmetic-only
outputs inside them still held to bit identity through `outRealExact`), and 1
is a `deviation` case.

Deep run `npm run oracle:deep -- 2000 v23-estimates`: 52 cases, 104000
records, 3194 C++ exceptions (all matched by the port), **4 644 558
floating-point outputs compared**. The 4 416 000 outputs of the exact cases
are **100.0000 % bit-identical**; the 200 558 outputs of the tolerance cases
are 96.96 % bit-identical with a **worst scaled error of 8.9e-16** (two ulps,
`Slerp.slerpCosAngle`); the deviation case deviates on 2000 of 2000 records.
All 55 tests pass.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 19.44.35207
x64 `/O2 /fp:precise`. `T = double` throughout: the port implements only the
double instantiation and takes the degree as a runtime argument, so each
record evaluates **every degree the header's `static_assert` admits** (a comma
fold over an index sequence), and every template instantiation is compared on
every record.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `Exp2Estimate.h` | `estimate`, `estimateRR`, `getMaxError` (degrees 1..7) | exact | pass |
| `Exp2Estimate.h` / `ExpEstimate.h` | `estimateRR.hugeArgument` | deviation | 2000/2000 deviate |
| `ExpEstimate.h` | `estimate`, `estimateRR`, `getMaxError` (1..7) | exact | pass |
| `Log2Estimate.h` | `estimate`, `estimateRR`, `getMaxError` (1..8) | exact | pass |
| `LogEstimate.h` | `estimate`, `estimateRR`, `getMaxError` (1..8) | exact | pass |
| `SinEstimate.h` | `estimate`, `estimateRR`, `getMaxError` (3, 5, ..., 11) | exact | pass |
| `CosEstimate.h` | `estimate`, `estimateRR`, `getMaxError` (2, 4, ..., 10) | exact | pass |
| `TanEstimate.h` | `estimate`, `estimateRR`, `getMaxError` (3, 5, ..., 13) | exact | pass |
| `ATanEstimate.h` | `estimate`, `estimateRR`, `getMaxError` (3, 5, ..., 13) | exact | pass |
| `ACosEstimate.h` | `estimate`, `getMaxError` (1..8) | exact | pass |
| `ASinEstimate.h` | `estimate`, `getMaxError` (1..8) | exact | pass |
| `SqrtEstimate.h` | `estimate`, `estimateRR`, `getMaxError` (1..8) | exact | pass |
| `InvSqrtEstimate.h` | `estimate`, `estimateRR`, `getMaxError` (1..8) | exact | pass |
| `ChebyshevRatio.h` | `ratio`, `ratios`, `ratioUsingCosAngle`, `ratiosUsingCosAngle`, `usingCosAngle.nearPi`, `throwParity` | tolerance (`sin`, `acos`); exact on the angle-0 branches | pass |
| `ChebyshevRatioEstimate.h` | `estimate` (1..16), `estimateR` (1..12), `getMaxError`, `getMaxErrorR` | exact | pass |
| `Slerp.h` | `slerp`, `slerpCosAngle`, `slerpMidpoint`, `throwParity` (N = 2, 3, 4) | tolerance (`acos`, `sin`); exact on the angle-0 branch and the +-0 slots | pass |
| `SlerpEstimate.h` | `slerpEstimate`, `slerpEstimateCosAngle`, `slerpEstimateMidpoint` (N = 2, 3, 4, degrees 1..16) | exact | pass |

Every computational public entry point of the 16 headers is reached: each
`XEstimate`, each `XEstimateRR`, each `GetXEstimateMaxError`, the four
`ChebyshevRatio*` functions, `ChebyshevRatioEstimate`,
`ChebyshevRatioEstimateR` and their two max-error getters, and the three
overloads each of `Slerp` and `SlerpEstimate`.

### Generators and the branches they reach

Each record draws one argument from each of eight modes, so every committed
record covers every population: uniform in the domain; a dyadic lattice
(`k/64`, `k*L/16`); the domain endpoints and their neighbours within 2 ulps
(and `-0` at a zero endpoint); values approaching an endpoint geometrically
(down to subnormal offsets); uniform values outside the domain; and a wild
population (signed zeros, `denorm_min`, random subnormals, `DBL_MIN`,
`DBL_MAX`, `+-inf`, NaN, powers of two across the whole exponent range,
`1 +- 1..2 ulps`). The range-reduced variants have their own populations:

- **`Sin/Cos/TanEstimateRR`**: `k*Q` for the branch thresholds (`Q = pi/2`,
  resp. `pi/4`) moved by up to 2 ulps; the exact `std::remainder` ties
  `odd*(P/2)` for `odd` in `{+-1, +-3, +-5, +-7}` (exact products: the
  significand of `pi` ends in three zero bits); huge arguments up to
  `DBL_MAX`; integers below `2^31`; `+-2^k`.
- **`Exp2/ExpEstimateRR`**: integers in `[-1080, 1030]` moved by up to 2 ulps
  (`y = 0` and `y = 1 - ulp`), subnormal results, underflow to zero, overflow,
  `|x| < 2^30` otherwise (see the deviation below).
- **`Log2/Log/Sqrt/InvSqrtEstimateRR`**: log-uniform over all positive doubles
  including subnormals, exact powers of two of both exponent parities moved by
  up to 2 ulps, `+-0`, negative values, `+-inf`, NaN.
- **`ATanEstimateRR`**: `+-1` moved by up to 2 ulps, `64/k`, magnitudes up to
  `DBL_MAX` (subnormal `1/x`), NaN.
- **`ChebyshevRatio`**: angle `+-0` (the `angle == 0` branch, `-0 == 0`),
  subnormal angles, `pi - 1..3 ulps`, `pi/2 +- 2 ulps`; cosines `1`, above 1,
  `+inf`, NaN (both comparisons fail, so NaN returns `t`), `1 - 2^-k`.
- **`Slerp`/`SlerpEstimate` pairs**: independent unit vectors; signed
  coordinate axes whose other components are random `+-0`; `q1 = q0` (dot
  within an ulp of 1, on either side of `cosA < 1`); obtuse pairs with the
  angle in `(pi/2, pi - 0.1]`; `q1 = normalize(q0 + 2^-k r)`; small rational
  components (`(+-1/2)^4`, signed permutations of `(3/5, 4/5)` and
  `(1/3, 2/3, 2/3)`); near-antipodal pairs `normalize(-q0 + 2^-k r)` for the
  midpoint overloads. `t` includes `0`, `-0`, `1`, `1/2`, `1/2 +- 1 ulp` (the
  `2t <= 1` split) and extrapolation in `[-1/4, 5/4]`. The cosA overloads
  apply the documented "negate q1 when the dot is negative" preprocessing on
  odd records.

Branch histogram of the 2000-record deep run (recomputed independently in the
replay; the `std::remainder` branches with an exact BigInt remainder):

| entry point | branches |
| --- | --- |
| `SinEstimateRR` | `r > pi/2` 3771 (372 ties), `r < -pi/2` 3929 (357 ties), `|r| <= pi/2` 7264, non-finite 307 |
| `CosEstimateRR` | `r > pi/2` 3825 (326 ties), `r < -pi/2` 3963 (356 ties), `|r| <= pi/2` 7239, non-finite 291 |
| `TanEstimateRR` | `y > pi/4` 4087 (360 ties), `y < -pi/4` 3936 (367 ties), `|y| <= pi/4` 6893, non-finite 357; `r > pi/2` / `r < -pi/2` never (dead, #57) |
| `Exp2EstimateRR` | normal 11658 (+703 with `y = 0`), subnormal 1678 (+11), underflow to 0 283 (+116), overflow 937 (+122), non-finite 492 |
| `Log2EstimateRR` | `p` even 6810, odd 5370, subnormal even 992 / odd 1394, `+0` 134, `-0` 143, negative 551, non-finite 606 |
| `SqrtEstimateRR` | `p` even 6669, odd 5510, subnormal 2405, `+0` 144, `-0` 121, negative 581, non-finite 570 |
| `ATanEstimateRR` | `|x| <= 1` 6447, `x > 1` 4668, `x < -1` 4730, NaN (else arm) 155 |
| `ChebyshevRatioUsingCosAngle` | `(-1, 1)` 10000, `>= 1` 4000, NaN 2000 |
| `Slerp` (dot) | acos branch 4394, angle-0 branch 1606 |
| `Slerp` (midpoint) | acos branch `2t <= 1` 2410 / `2t > 1` 1855, angle-0 branch 949 / 786, `dot < -0.99` 727 |
| `SlerpEstimate` (dot) | angle `<= pi/2` 4241, `> pi/2` 1759 |

The replay's last test fails if one of the required branches (the RR arms,
ties, non-finite and `+-0` arguments, the angle-0 and NaN branches of the
ratios, both sides of `2t <= 1`) is missing from the replayed golden set. The
committed 20-record goldens reach every branch in the table, including 2 to
4 remainder ties per reduction; only the `subnormal, y = 0` combination of
`Exp2EstimateRR`, a `+0` argument of `LogEstimateRR` and a `-0` argument of
`InvSqrtEstimateRR` are left to the deep run.

### Where a tolerance is used, and why

Every estimate is `+ - * /`, `std::sqrt`, `std::fabs`, comparisons and the
IEEE-exact library functions `std::floor`, `std::frexp`, `std::ldexp` and
`std::remainder`, so all 34 scalar cases, the four `ChebyshevRatioEstimate`
cases and the three `SlerpEstimate` cases are `exact: true` and bit-identical
on every deep-run record. That includes the port's own `ldexp` (two-step
scaling), `frexp` (exponent-field read) and exact long-division `remainder`,
at subnormal, overflowing and `DBL_MAX`-sized arguments and at the remainder
ties, and every published max-error constant.

The ten tolerance cases call libm on the compared path:

- `ChebyshevRatio.ratio`, `ratios` — `std::sin`. Worst deep-run scaled error
  4.1e-16.
- `ChebyshevRatio.ratioUsingCosAngle`, `ratiosUsingCosAngle` — `std::acos`,
  `std::sin`. The generator keeps `cosAngle >= -1 + 2^-20`, where the angle is
  at least 1.4e-3 below pi and a 1-ulp difference in `acos` moves the ratio by
  at most 3.2e-13 relative (`dA/(pi - A)`). Worst 4.4e-16.
- `ChebyshevRatio.usingCosAngle.nearPi` — the same functions for
  `cosAngle` in `(-1, -1 + 2^-20]`, compared per output with the conditioning
  tolerance `max(1e-12, 4 * 2^-51 / sqrt(2(1 + cosAngle)))` (up to 1.7e-7 at
  the smallest gap). Measured worst: **3.9e-16**; the two runtimes' `acos`
  never differed by enough to show the conditioning.
- `ChebyshevRatio.throwParity` — valid records call `sin`/`acos`. Worst 1.8e-16.
- `Slerp.slerp`, `slerpCosAngle`, `slerpMidpoint`, `throwParity` —
  `ChebyshevRatiosUsingCosAngle`. The dot overloads keep the angle at most
  `pi - 0.1` (a pair with dot below -0.995 has `q1` negated), where a 1-ulp
  angle difference moves the result by about `ulp/(pi - A)^2 = 4.4e-14`; the
  midpoint overload reaches the near-antipodal pairs through its half angle.
  Worst 8.9e-16.

Inside those cases the exactly computed outputs are held to bit identity with
`outRealExact`: the `angle == 0` / `!(cosAngle < 1)` branches (which return
`t` and `{1 - t, t}`), and every Slerp slot in which both combined inputs are
`+-0` (the result is a signed zero whatever the libm ratios are). No branch
in this group is decided by a libm value: the ratio branches compare the
recorded angle or cosine, and the Slerp branches compare the dot product or
the recorded cosine, which are exact arithmetic on both sides.

## Port defects fixed

**`Slerp.h` and `SlerpEstimate.h` (all six functions): the zero seed of the
result accumulation.** Upstream writes

```cpp
std::array<T, N> result{};
result.fill(zero);
for (size_t i = 0; i < N; ++i) { result[i] += f[0] * q0[i] + f[1] * q1[i]; }
```

and the port assigned `result[i] = f[0] * q0[i] + f[1] * q1[i]`. The two
differ exactly when both products are `-0` (both input components `-0` with
nonnegative ratios, or a `-0` input times a zero ratio): IEEE 754 gives
`0 + (-0) = +0`, so upstream returns `+0` where the port returned `-0`. It is
the "accumulation seed" class of ORACLE.md (v01, `GVector` dot). Size: the
sign of zero only, on 305 (`Slerp.slerp`), 295 (`SlerpEstimate.slerpEstimate`)
and 287 (`SlerpEstimate.slerpEstimateCosAngle`) of 2000 deep-run records;
replaying the committed goldens through the unfixed port fails all six
Slerp/SlerpEstimate main cases, on 17 of their 120 records. The port now
fills with `0` and accumulates with `+=` exactly as upstream does
(`src/Slerp.ts`, `src/SlerpEstimate.ts`, with a comment);
`test/Slerp.test.ts` and `test/SlerpEstimate.test.ts` pin `+0` in both
branches, every overload and every degree, and require the direct assignment
to give `-0` on the same inputs (both tests fail on the old code).

**Reference-build caveat found on the way (MSVC, not upstream).** In the
first build the case bodies called `Slerp<double, N>` directly, and MSVC
19.44 inlined it into `SlerpBody<2>`: knowing `result[i] == 0.0` after the
`fill`, it emitted `mulpd`/`addpd` for `f0*q0 + f1*q1` and stored that,
**dropping the `0.0 +`**, which is not value-safe for `-0` and contradicts the
documented `/fp:precise` rule that `-0.0` is processed per IEEE 754. The
out-of-line `Slerp<double, 2>` in the same object keeps
`addsd xmm1, QWORD PTR [rbx]`. The golden of that build therefore had `-0` in
2 of the zero slots (N = 2 dot overload, N = 4 cosA overload) and `+0` in all
others, so the sign in those slots was a property of the inlining decision,
not of the upstream source. The case file now calls upstream through
`__declspec(noinline)` wrappers; regenerating changed exactly 11 outputs, all
`-0 -> +0` in Slerp zero slots (9 of them in `throwParity`, where the default
tolerance had hidden them), and the replay holds every zero slot to bit
identity on all 2000 x 6 Slerp/SlerpEstimate deep records, which verifies that
the wrappers keep the seed. Any other group whose upstream code zero-fills an
accumulator and whose case body lets MSVC inline it can have the same
problem.

## Deliberate deviations demonstrated

**`Exp2EstimateRR` / `ExpEstimateRR` for `floor(x) >= 2^31`**
(`Exp2Estimate.estimateRR.hugeArgument`, new, see the suspects below).
Upstream computes `static_cast<int32_t>(static_cast<double>(p))` with
`p = floor(x)`; outside the `int32_t` range the conversion is undefined
behaviour and MSVC's `cvttsd2si` produces `INT_MIN`, so
`std::ldexp(poly, INT_MIN)` returns `+0` where `2^x` overflows to `+inf`.
The port keeps the exponent as a double and returns `+inf`. The case draws
`x = u*2^k`, `k` in `[31, 1021]` (so `x * (1/ln 2)` stays finite and `y` is a
number) and emits all seven degrees of both functions: 2000 of 2000 records
deviate, every output `+0` (C++) against `+inf` (port). The main
`estimateRR` generators keep `|x| < 2^30`. The port was already correct; the
oracle made the difference visible, and this group adds the explanation to
`src/Exp2Estimate.ts` and a regression test (`test/Exp2Estimate.test.ts`,
`+inf` for `x >= 2^31`, `+0` for `x < -2^31`).

The six preserved #57 findings (ASin/ACos max-error digits, LogEstimate's
log2 bound, `SqrtEstimateRR(0)`, the TanEstimateRR comment, the
ChebyshevRatioEstimate R-variant comment, the Slerp header formula) are
either values the port reproduces or comments only; they agree bit for bit
(`SqrtEstimateRR(+-0)` is in the committed goldens, the max-error tables are
compared constant by constant) and need no deviation case. The port's
corrected comments change no computation.

## Independent reference checks

Run in the replay on every record, over the committed goldens and the
2000-record deep run (ORACLE.md, "Agreement is not correctness"); the truth is
`Math.*` in double precision.

- **Every estimate inside its documented domain is within its published
  `GetMaxError` bound** of the true function, for every degree:
  `Exp2`/`Exp` on `[0, 1]`/`[0, ln 2]`, `Log2`/`Log`/`Sqrt`/`InvSqrt` on
  `[1, 2]`, `Sin`/`Cos` on `[-pi/2, pi/2]`, `Tan` on `[-pi/4, pi/4]`, `ATan` on
  `[-1, 1]`, `ACos`/`ASin` on `[0, 1]`; and the range-reduced variants as a
  relative bound (`Exp2RR`, `ExpRR`, `SqrtRR`: the reduced absolute bound;
  `InvSqrtRR`: `sqrt(2)` times it), absolutely (`Log2RR`, `LogRR`, `ATanRR` on
  all reals), or for `|x| <= 1000` (`SinRR`, `CosRR`, beyond which the double
  `2*pi` of the reduction is no longer the true period to the bound's
  precision) and for `|y| <= pi/4` (`TanRR`). About 2.0 million checks. The
  worst `|error|/bound` is 0.99996 to 1.00000 for every non-RR estimate and
  degree (the published minimax bounds are attained, not loose), 0.71 to 0.99 for the RR
  variants, and **0.693 for `LogEstimate`**: its published bound is the log2
  bound (#57), and the extra `logTight` check confirms the true error reaches
  `ln 2` times it (ratio 1.00000). The only values above a bound are 14
  `ChebyshevRatioEstimateR` checks at degrees 9 to 12 (worst ratio 1.0034 at
  degree 12, whose bound is 3.3e-14, i.e. an excess of 1.1e-16) and one
  `Log2EstimateRR` check at degree 5 (ratio 1.000000): evaluation rounding at
  the bound's attained maximum, inside the stated 1e-15 slack.
- **`ChebyshevRatioEstimate` / `ChebyshevRatioEstimateR`** are within their
  bounds of `sin(tA)/sin(A)`, `A = acos(x)`, for `t` in `[0, 1]` and `x` in
  `[0, 1]` resp. `[cos(pi/4), 1]` (23 828 resp. 22 490 samples per degree).
- **`ChebyshevRatio` identity** `f(t, A) sin(A) = sin(t A)` holds to
  `1.1e-16 * max(1, |sin(tA)|)` (51 629 checks) for all four functions.
- **`Slerp`**: for `t` in `[0, 1]` the result is a unit vector (worst
  `||r| - 1|` 1.2e-14) whose angle from `q0` is `t*A` and from `q1` is
  `(1 - t)*A` (worst 8.8e-15), measured with `2*atan2(|a - b|, |a + b|)`. For
  the midpoint overload the bound is scaled by the conditioning of the
  documented preprocessing (`1e-12 + 4e-16/cosAH^2`: `1 + cosA` cancels near
  the antipode, and `cosAH`, `qh` carry a relative error of about
  `2^-52/(4 cosAH^2)`); with that scaling the worst ratio is 0.071. With a
  flat 1e-12 bound, 5 near-antipodal records (`cosAH` about 2e-3) reached
  4.2e-12 — the error of the recorded `qh`, not of `Slerp`.
- **`SlerpEstimate`** is within `GetChebyshevRatioEstimateMaxError(D)` times
  `|a_i| + |b_i|` of the true slerp for every degree and component, on the
  pairs inside the documented angle range (dot overloads: `dot >= 0`;
  midpoint: `cosAH >= 1e-3`, compared on the half arc it evaluates).

## Sensitivity

The cases were checked to discriminate the plausible alternative forms of the
suspect computations, recomputed from the 2000-record deep run and compared
with the C++ outputs (records of 2000 that would have differed):

| computation | alternative | records differing |
| --- | --- | --- |
| `ChebyshevRatioEstimate` `b` coefficient `u*(i+1)/(2i+3)` | `u*((i+1)/(2i+3))` | 1593 |
| `ChebyshevRatioEstimate` `term *= (b - a*sqr)*y` | `term = term*(b - a*sqr)*y` | 1265 |
| `ChebyshevRatioEstimate` `a` coefficient `u/((i+1)(2i+3))` | `u/(i+1)/(2i+3)` | 6 |
| `SinEstimateRR` exact `std::remainder` | `x - round(x/2pi)*2pi` | 2000 |
| `TanEstimateRR` exact `std::remainder` (ties to even) | `fmod`, then shift into `[-pi/2, pi/2]` | 341 |
| `Exp2EstimateRR` `std::ldexp` | `poly * 2**p` | 28 |
| `Exp2Estimate` Horner | power form `sum c_i x^i` | 2000 |
| `Slerp`, `SlerpEstimate` zero seed | direct assignment | 287 to 305 per case |

`SinEstimateRR` and `CosEstimateRR` cannot tell `std::remainder` from an exact
`fmod` shifted into `[-pi, pi]`: the two differ only at the ties `r = +-pi`,
where `pi - r` and `-pi - r` are both `+0` (0 records); `TanEstimateRR` ties
do discriminate (341 records).

## Not covered

- **Invalid degrees.** Upstream rejects them with `static_assert` at compile
  time; the port's runtime `logAssert` has no C++ counterpart to compare.
- **Mismatched or too-small dimensions** of the `Slerp`/`SlerpEstimate`
  arrays: compile-time in upstream (`std::array<T, N>`, `static_assert(N >=
  2)`), runtime asserts in the port.
- **`T = float`** (and any other `T`): the port implements the double
  instantiation only.
- **`Exp2EstimateRR`/`ExpEstimateRR` with `floor(x*scale) < -2^31`**: the
  same undefined conversion; MSVC's `INT_MIN` happens to give the right `+0`,
  and the port returns `+0` as well (pinned by `test/Exp2Estimate.test.ts`),
  but a golden of undefined behaviour is not recorded.
- **The `frexp` exponent of `+-inf` and NaN** is unspecified (MSVC returns
  `-1`, the port computes 1025). The RR cases do feed `+-inf` and NaN and
  agree bit for bit, because `poly` is then infinite or NaN and the exponent
  cannot reach the output.
- **`TanEstimateRR`'s `r > pi/2` / `r < -pi/2` arms** are unreachable (#57);
  the replay asserts that no record reaches them.

## Upstream bug suspects

**`Exp2EstimateRR` converts `floor(x)` to `int32_t` without a range check
(undefined behaviour; wrong result `+0` on MSVC for `x >= 2^31`).** New, not
in `docs/UPSTREAM-FINDINGS.md`. `Exp2Estimate.h`:

```cpp
T p = std::floor(x);
T y = x - p;
T poly = Exp2Estimate<T, Degree>(y);
int32_t power = static_cast<int32_t>(static_cast<double>(p));
T result = std::ldexp(poly, power);
```

The header documents `x` as "any real number". For `floor(x)` outside
`[-2^31, 2^31 - 1]` the conversion is undefined ([conv.fpint]); MSVC x64
produces `INT_MIN`, so `Exp2EstimateRR<double, 3>(3e9)` returns `0` (probe
and `Exp2Estimate.estimateRR.hugeArgument`, 2000 of 2000 records) where
`2^x` overflows to `+inf`; `ExpEstimateRR` inherits it for
`x >= 2^31 ln 2`. Below `-2^31` the same conversion happens to produce the
correct `+0`. A fix would clamp `p` to, e.g., `[-1100, 1100]` before the
conversion (every `ldexp` result is already saturated there). The port keeps
the exponent as a double (`+inf`/`+0`), documented in `src/Exp2Estimate.ts`.
