# Group 45 (roots-polynomials): `RootsCubic.h`, `RootsQuartic.h`

28 cases, 2000 records each in the deep run (`npm run oracle:deep -- 2000
v45-roots`, 58 s wall, all 29 tests green). 15 cases are bit-for-bit
(`exact: true`), 10 carry the libm tolerance 1e-9, 3 are deliberate
`deviation` cases.

Both headers classify roots and multiplicities with exact rational arithmetic
(`BSRational<UIntegerAP32>`) and then estimate the root values either by
bisection (`useBisection = true`) or by a closed form (`useBisection =
false`). The bisection path uses only `+ - * /`, comparisons and `gte::FMA`,
so it is held to bit-identity; the closed form calls `std::pow`,
`std::atan2`, `std::cos` and `std::sin` through the `BSRational` overloads
that convert to `double`, call libm and convert back, so it carries a
tolerance. The `SolveBiquadratic` path and the `d0 = 0` branches use only
`std::sqrt`, which is correctly rounded in both runtimes, and are compared
bit-for-bit even in the closed-form variant.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `RootsCubic.h` | `solve`, `solveMonic`, `solveDepressed`, `computeDepressedRoots` (bisection + closed form), `solveRational`, `solveDepressedRational` (bisection), `solveMonicRational` (closed form), `closedForm.illConditioned` | 7 exact, 5 tol 1e-9 | pass |
| `RootsCubic.h` deviation | `solveDepressed.deviation.bisectionBound` | deviation, exact | 2000/2000 deviate |
| `RootsQuartic.h` | `solve`, `solveMonic`, `solveDepressed`, `computeDepressedRoots` (bisection + closed form), `solveRational`, `solveDepressedRational` (bisection), `solveMonicRational` (closed form), `solveBiquadratic` (bisection + closed form) | 8 exact, 5 tol 1e-9 | pass |
| `RootsQuartic.h` deviations | `solveDepressed.deviation.complexPairs`, `solveDepressed.deviation.staleSqrt` | deviation | 2000/2000 and 313/2000 deviate |

Every computational public entry point of both headers is covered for both
template instantiations the port exposes: `T = double` (`solve`,
`solveMonic`, `solveDepressed`), `T = Rational` (`solveRational`,
`solveMonicRational`, `solveDepressedRational`) and the `T`-independent
`ComputeDepressedRoots`. The private `HasZeroValuedRoots` and
`ComputeClassifiers` helpers are reached through the public entry points by
the `g0 = 0` / degree-drop generator modes.

### Generators

Eight modes cycle on `io.index() % 8`, each recording the same number of
doubles (degree+1 for the general form, degree for monic, degree-1 for
depressed):

0. small integer lattice `[-5, 5]`;
1. uniform reals `[-10, 10]`;
2. expanded from distinct integer roots (all real, all multiplicity 1);
3. expanded from integer roots with a multiplicity pattern — `(x-a)^2(x-b)`,
   `(x-a)^3`, `(x-a)^2(x-b)^2`, `(x-a)^3(x-b)`, `(x-a)^4`, with the depressed
   forms constructed so that the roots sum to zero. The coefficients are
   integers, so the rational discriminant is *exactly* zero and the exact
   classifier branches run at exact equality;
4. a common power-of-ten scale `10^k`, `k` in `[-100, 100]`;
5. per-coefficient scales `10^k`, `k` in `[-30, 30]`;
6. leading coefficient(s) zero (degree drop into the cubic/quadratic/linear
   solvers);
7. zero-valued roots, including a `-0.0` constant term and runs of zero
   low-order coefficients.

The biquadratic cases build `(x^2 - a)(x^2 - b)` from integers so that all
four classifications of `SolveBiquadratic` (four real roots, two real roots
and a complex pair, two complex pairs, two real roots of multiplicity 2) are
reached exactly.

## Port defects fixed

None. Every bisection case is bit-identical to the MSVC build on all 2000
records (0 inexact real outputs out of ~4000 per case), including the
`1e±100` scale modes, the lattice modes and the `T = Rational`
instantiations.

## Deliberate deviations demonstrated

All three are the `RootsCubic.h` / `RootsQuartic.h` entries of
`docs/UPSTREAM-FINDINGS.md`, issue
[#340](https://github.com/gradientspaceai/gtengine-js/issues/340). The main
cases reject the inputs on which upstream is defective with probes that are
verbatim copies of upstream's own control flow (`CubicBisectionSound`,
`QuarticDepressedSound`, `QuarticItem2Defect` in the case file).

1. **Invalid cubic bisection bound** (item 1). The `signDelta < 0` branch
   bisects `x^3 + d1*x + d0` over `[-b, b]` with `b = max(1,|d0|,|d1|)`,
   which is not a root bound; the port uses Cauchy's `1 + max(|d0|,|d1|)`.
   `|x|^3 <= |d1||x| + |d0|` gives `|x| <= sqrt(2b)`, so the bound can only
   fail for `b < 2`, and the deviation generator draws `d0, d1` in `[-2, 2]`.
   All 2000 records deviate, by up to 0.24 relative.

   **The fix is confined.** A separate probe compiled against upstream
   (40000 draws over four generator modes, 38048 of them with a sound
   upstream bound) bisected each depressed cubic over both intervals: the
   two runs returned *identical* `[xMin, xMax]` pairs on all 38048 sound
   records and differed on all 1952 unsound ones. `gte::FMA` rounds
   `x*(x^2+d1) + d0` once, so the evaluation is exactly zero only at a true
   double root of the floating-point function and both bisections converge
   on the same pair of adjacent doubles. This is why
   `RootsCubic.solveDepressed.bisection` and the quartic cases that go
   through the resolvent cubic are bit-identical.

2. **Incomplete complex-pair classification** (item 2). `if (signDelta > 0 &&
   rD2.GetSign() > 0) return 0;` is only half the criterion. The deviation
   generator builds `(x^2 + p*x + q)(x^2 - p*x + r)` with `p^2 < 4q`,
   `p^2 < 4r` and `q + r <= p^2`: two complex-conjugate pairs, positive
   discriminant, `d2 = q + r - p^2 <= 0`. Upstream reports four roots, the
   port none; all 2000 records deviate structurally.

3. **Stale reads of the reused square-root array** (item 3). `beta =
   sqrt(T^2 - d0)` is mathematically zero exactly when `d1` is zero, so the
   generator draws `|d1|` near `1e-200`: `T^2 - d0` is then pure round-off
   and its sign decides whether upstream writes index 1 of `rQRoots` or
   re-reads `alpha = sqrt(2T - d2)` from the previous extraction. 313 of
   2000 records deviate by order 1. The case is compared at 1e-6 rather than
   bit-for-bit so that the libm disagreement of the resolvent cubic
   (measured maximum 6.7e-9) is not counted as a deviation; it runs the
   closed form so that item 1 cannot be the cause.

## Tolerances

`1e-9` on the ten closed-form cases; measured maximum over the deep run
**6.53e-11** (`RootsCubic.solveDepressed.closedForm`). Cause: `std::pow`,
`std::atan2`, `std::cos` and `std::sin` differ between the MSVC runtime and
V8 in the last bit, and the depressed cubic's three-root formula builds its
roots as `-t0 - t1` and `-t0 + t1`, which amplifies that by the cancellation
factor `max|t| / |root|`. The generators of the closed-form cases reject
records whose cancellation factor exceeds `1e6` (`CANCELLATION_LIMIT`), which
bounds the relative disagreement by about `1e-15 * 1e6 = 1e-9`; the
acceptance threshold and the tolerance were chosen together. Without it the
deep run showed disagreements of 1.9e-2 and worse: for
`x^3 - 8.86e33 x + 1.62e33` the outer roots `+-9.4e16` agree to every digit
while the middle root is `-11.257` in MSVC and `0` in V8.

The rejected regime is not dropped:
`RootsCubic.solveDepressed.closedForm.illConditioned` inverts the predicate
and emits the root count and the multiplicities, which exact rational
arithmetic decides; those agree on all 2000 ill-conditioned records.

`1e-6` on `RootsQuartic.solveDepressed.deviation.staleSqrt`, justified above.

## Independent reference checks

Run once over the 2000-record deep run with exact BigInt rational arithmetic
(the scripts are throwaway; the numbers are reproduced here).

* **Root count versus a Sturm chain.** For ten cases the polynomial was
  rebuilt exactly from the recorded coefficients, its square-free part taken
  over the rationals and the number of distinct real roots counted by sign
  variations of the Sturm chain at `-inf` and `+inf`. Over 19971 records the
  reported root count equalled the Sturm count **every time** (0 mismatches),
  and the multiplicities never summed to more than the degree.
* **Residual `|p(r)|` relative to `max_i |a_i| |r|^i`.** Over the six
  generator modes that are not extreme-scale: worst 3.9e-10
  (`RootsQuartic.solve.bisection`), typically below 1e-12.
* **Lattice constructions.** For the modes built from integer roots, `p(r)`
  is *exactly* zero over the rationals for every reported root in about 90%
  of records (448/500, 453/500, 461/500 ...), and the multiplicities sum to
  the degree in 500/500. The remainder are the records where the depressed
  root is irrational as a double even though the original roots are integers.

## Not covered

Nothing. Both headers are fully ported and every public entry point has a
case.

## Upstream bug suspects

**The inverse depression transform loses all significance when a root is much
smaller than `m2/3` (or `m3/4`).** New, not in
`docs/UPSTREAM-FINDINGS.md`. `RootsCubic<T>::Solve` shifts the roots by
`rM2Div3` in exact rational arithmetic, but the depressed root it shifts is
only a double-precision bisection (or closed-form) estimate, so the shift
subtracts two nearly equal quantities whose difference is below the
resolution of the estimate.

Explicit input (deep-run record 29 of `RootsCubic.solve.bisection`):

```
g = (2.996630553652356e-25, -496.23907294959224,
     6.901865787578214e-15, -9.341590182409976e-6)
```

`m2/3 = -2.4627733440124166e-10` and the depressed root is
`-2.462773344012417e-10`; upstream returns `x = -1.946996467443576e-26`,
whereas the true root is `+6.0387e-28` (the linearisation `-g0/g1`, since the
cubic and quadratic terms are 10 orders of magnitude smaller there). The
residual `|p(x)|` is 1.03 times the coefficient scale — that is, `x` is not
a root at all, and even its sign is wrong. The same shape appears in
`RootsQuartic` (`m3/4`) and in the closed-form variants.

This is a conditioning limitation with no exact defective condition to guard
on, so per PORTING.md the port preserves it; it is pinned by the
bit-for-bit `RootsCubic.solve.bisection` and `RootsQuartic.solve.bisection`
cases, whose generator modes 4 and 5 reach it on about 3% of records. It is
invisible to the classification: the root *count* and the multiplicities stay
correct (the Sturm check above found no mismatch), only the root values are
wrong, so callers that use the count are unaffected.
