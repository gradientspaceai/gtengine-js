# Group 44 (roots-polynomials): the 13 remaining root and polynomial headers

40 cases, 2000 records each in the deep run (`npm run oracle:deep -- 2000
v44-roots`, 87 s wall, all 41 tests green). 34 cases are bit-for-bit
(`exact: true`), 2 carry the libm tolerance 1e-9, 4 are deliberate
`deviation` cases.

Almost everything in this group is `+ - * /`, comparisons, `std::sqrt`,
`std::fabs` and `std::max`, so the expectation is bit-identity. The only libm
users are `RootsPolynomial`'s low-degree closed forms (`std::pow`,
`std::atan2`, `std::cos`, `std::sin` in `SolveDepressedCubic` and, through
it, `SolveDepressedQuartic`).

Instantiations follow the port: `RootsLinear`/`RootsQuadratic` are exercised
in both the `T = double` and the `T = BSRational<UIntegerAP32>` forms;
`RootsPolynomial` is instantiated with its `Rational` template parameter equal
to `double`, which is what the port implements (docs/API.md, "Root finders
report different things under similar names"); `RootsBisection1` and
`RootsBisection2` only in the floating-point form, the only one ported.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `Polynomial1.h` | `arithmetic` (all 12 operators plus the 6 comparisons), `evaluateAndTransforms` (Horner, two derivatives, inversion, two translations), `mutators` (`EliminateLeadingZeros`, `MakeMonic`, `SetDegree`, `SetCoefficients`), `divide`, `greatestCommonDivisor` (both argument orders), `squareFreeFactorization` | 6 exact | pass |
| `PolynomialCurve.h` | `evaluate` (the four stored polynomials per component and the jet at orders 0..3) | exact | pass |
| `PolynomialCurve.h` deviation | `deviation.constructed` | deviation | 2000/2000 |
| `PolynomialRoot.h` | `bisect` (all four `PolynomialRootBisect` exits), `compare` (`operator==`, `operator<`, the default constructor) | 2 exact | pass |
| `RootsBisection.h` | `find`, `findWithValues` (both overloads, iteration counts included) | 2 exact | pass |
| `RootsBisection1.h` | `find` (both overloads), `throwParity` (the `tMin < tMax` assert, 1129 of 2000 records throw on both sides) | 2 exact | pass |
| `RootsBisection1.h` deviation | `deviation.maxIterationsOne` | deviation | 1302/2000 |
| `RootsBisection2.h` | `find` | exact | pass |
| `RootsBisection2.h` deviation | `deviation.staleOutputs` | deviation | 2000/2000 |
| `RootsBrentsMethod.h` | `find` (every validation rejection and every termination criterion) | exact | pass |
| `CubicRootsQR.h` | `solve`, `solveMatrix` (iteration counts, root counts, all root slots and the modified matrix) | 2 exact | pass |
| `QuarticRootsQR.h` | `solve`, `solveMatrix` | 2 exact | pass |
| `RootsLinear.h` | `solve`, `solveMonic`, `solveRational` (both rational overloads) | 3 exact | pass |
| `RootsQuadratic.h` | `solve.bisection`, `solve.closedForm`, `solveMonic`, `solveDepressed`, `computeDepressedRoots`, `solveRational` (all three rational overloads) | 6 exact | pass |
| `RootsPolynomial.h` | `solveQuadratic` (exact: `std::sqrt` only), `solveCubic`, `solveQuartic`, `getRootInfo` (all three), `find`, `findBounded` | 4 exact, 2 tol 1e-9 | pass |
| `RootsGeneralPolynomial.h` | `solve`, `solveRational` | 2 exact | pass |
| `RootsGeneralPolynomial.h` deviation | `deviation.zeroPadding` | deviation | 851/2000 |

### Generators

Eight modes cycle on `io.index() % 8`: small integer lattice; uniform reals;
expanded from distinct integer roots; expanded from integer roots with
repeats (double, triple and quadruple roots); a common `10^k` scale with `k`
in `[-100, 100]`; per-coefficient `10^k` scales with `k` in `[-30, 30]`;
leading coefficient(s) zero; and zero-valued roots with a `-0.0` constant
term. Polynomial degrees are themselves recorded inputs, so the degree-drop
and identically-zero paths are covered. The bisection and Brent cases draw
their interval endpoints from overlapping lattice ranges, so `t1 <= t0`, an
exact endpoint root, no sign change, and `maxIterations` exhausted all occur.
The target functions are Horner evaluations of the recorded coefficients,
written character for character the same on both sides.

## Port defects fixed

None. Every arithmetic case is bit-identical to the MSVC build on all 2000
records: 0 inexact real outputs out of 78772 (`Polynomial1.arithmetic`),
63189 (`PolynomialCurve.evaluate`), 40000 (`QuarticRootsQR.solveMatrix`) and
so on, including the `1e±100` scale modes, signed-zero coefficients and the
`T = Rational` instantiations.

## Deliberate deviations demonstrated

1. **`RootsBisection1` with `maxIterations == 1`** (issue
   [#84](https://github.com/gradientspaceai/gtengine-js/issues/84)). The loop
   `for (iteration = 2; iteration <= mMaxIterations; ++iteration)` never runs,
   so upstream returns `tRoot` and `fAtTRoot` unassigned and the caller reads
   back whatever it passed in. The case seeds both with a sentinel in
   `[1, 2]`, which is what a caller reusing a variable sees; the port returns
   explicit zeros. 1302 of 2000 records deviate (the rest are the records
   where the endpoint signs agree or an endpoint is an exact root, which
   upstream writes explicitly).
2. **`RootsBisection2` stale outputs** (issues #84 and
   [#152](https://github.com/gradientspaceai/gtengine-js/issues/152)). Built
   on the same defect: the case constructs one bisector with
   `xMaxIterations = 1`, calls it twice, and makes the first call's
   `F(xMin, y)` exactly zero so that `RootsBisection1` writes `mXRoot = xMin`.
   The second call leaves the member untouched while returning a successful
   iteration count, so upstream reports the first call's root. All 2000
   records deviate.
3. **`RootsGeneralPolynomial` zero padding** (issue
   [#340](https://github.com/gradientspaceai/gtengine-js/issues/340), item 4).
   `std::vector<Rational> rP(p.size())` is allocated at the untrimmed size
   while only `0..degree` are assigned, so the solver runs on a zero-padded
   polynomial. Two consequences: the Cauchy bound takes the maximum over
   `rP[0..size-2]`, which now includes the monic `1` of the real leading term,
   and the derivative recursion descends through phantom degrees whose base
   case `degree == 1` is never reached with the real linear derivative, so
   that derivative root is bisected in floating point instead of being
   returned as the exact rational `-rP[0]/rP[1]`. The effect is a last-bit one
   on most draws, so the case's generator rejects until the padded run really
   differs from the unpadded one (which is what the port computes); 851 of
   2000 records then deviate, by up to 1.9e-15 relative. Without that
   rejection only about 1% of records differ, which is why the upstream
   finding records that "values survived in every configuration exercised".
4. **`PolynomialCurve` never sets `mConstructed`** (issue
   [#319](https://github.com/gradientspaceai/gtengine-js/issues/319)). Both
   upstream constructors leave `ParametricCurve::mConstructed` false, so
   `operator bool` reports failure for every valid curve. All 2000 records
   deviate (upstream `false`, port `true`); `Evaluate` is unaffected, which
   the separate `PolynomialCurve.evaluate` case confirms bit-for-bit.

## Tolerances

`1e-9` on `RootsPolynomial.solveCubic` and `RootsPolynomial.solveQuartic`.
Cause: `SolveDepressedCubic` computes its three roots as
`rhoPowThird*cos(theta/3)` plus or minus `sqrt(3)*rhoPowThird*sin(theta/3)`
with `rhoPowThird = pow(rhoSqr, 1/6)` and `theta = atan2(...)`, and MSVC and
V8 round `pow`, `atan2`, `cos` and `sin` differently in the last bit.
Measured maximum over the deep run: **1.55e-15** (cubic) and **1.33e-15**
(quartic), so the tolerance has six orders of headroom; it is set at 1e-9
because the generator's acceptance threshold (below) bounds the error at
about `1e-15 * 1e6 = 1e-9` rather than at the measured value.

Generator acceptance (chosen together with that tolerance): a record is
redrawn when upstream's own root map has a root whose magnitude is more than
`1e6` times smaller than the largest root, since such a root is the
difference of two nearly equal libm-derived terms. Without it the deep run
showed differences of 1.95e-3.

`RootsPolynomial.solveQuadratic` is `exact: true`: its depressed solver calls
only `std::sqrt`, which is correctly rounded in both runtimes.

## Independent reference checks

Run once over the 2000-record deep run with exact BigInt rational arithmetic.

* **`Polynomial1.divide` satisfies `P = Q*D + R` exactly**: the rational
  residual of `Q*D + R - P`, relative to the largest coefficient magnitude, is
  at most **2.8e-16** over all 2000 records.
* **QR eigenvalue residuals**: over the six non-extreme-scale generator modes
  and excluding polynomials with a zero constant term (where the relative
  residual of a root at 0 is meaningless), `|p(r)| / max_i |a_i||r|^i` is at
  most **3.95e-13** for `CubicRootsQR` (1757 roots) and **1.23e-12** for
  `QuarticRootsQR` (1962 roots).
* **`RootsGeneralPolynomial.solve`**: residual at most **1.20e-15** over 2154
  roots, and over 1500 Sturm comparisons the number of reported roots never
  exceeded the exact number of distinct real roots (it was smaller in 28,
  which is the documented behaviour: a root of even multiplicity has no sign
  change and is not reported).
* **`RootsPolynomial.find`**: the reported root values on the lattice
  constructions are exactly the construction's roots (`(x+3)^3` gives
  `-3, -3, -3`; `(x-3)^4(x+1)·(-3)` gives `-1, 3, 3, 3, 3, 3`), but a root of
  odd multiplicity `m` is reported `m` times, once per bounding subinterval
  whose endpoint value is exactly zero. In 45 of 264 well-converged records
  the count therefore exceeds the Sturm count of distinct real roots. This is
  upstream's documented deprecation reason and the port's class comment says
  the same; `RootsGeneralPolynomial` never does it.

## Not covered

* **`RootsBisection1`/`RootsBisection2` arbitrary-precision instantiation.**
  The port only has the floating-point one (`porting-status.json`: "floating-point
  instantiation only; arbitrary-precision path deferred until BSNumber exists"),
  so the `precision` constructors, `RoundInitial` and the `std::ldexp`-based
  `RoundAverage` have no counterpart to compare against.
* **`RootsGeneralPolynomial` with `useThreading = true` in the top-level
  call.** It is drawn and recorded, and the results agree, because upstream's
  threaded branch collects the subinterval roots in the same order as the
  sequential one after joining; the port ignores the flag. The `std::thread`
  scheduling itself is not comparable and nothing depends on it.
* **`Polynomial1::SquareFreeFactorization` on inputs where upstream does not
  terminate.** The upstream `do { ... } while (b.GetDegree() > 0)` loop is
  unbounded and spins forever whenever the floating-point
  `GreatestCommonDivisor` fails to reduce the degree of `b` (issue
  [#83](https://github.com/gradientspaceai/gtengine-js/issues/83), reproduced
  with `f = (t-2)^2 (t+2)(t+3)(t-4)`). A non-terminating run cannot be
  recorded, so the generator rejects those inputs with a verbatim copy of the
  routine carrying the port's own iteration cap, `degree(f) + 1`.
* **`Polynomial1`'s compound assignment operators** (`+=`, `-=`, `*=`, `/=`).
  Upstream defines each as the corresponding binary operation followed by
  assignment and the port omits them for that reason, so the binary cases
  cover them.

## Upstream bug suspects

**`RootsPolynomial`'s `std::map<Real, int32_t>` receives NaN keys, which is
undefined behaviour.** New, not in `docs/UPSTREAM-FINDINGS.md` (the file
records the *conditioning* of `RootsPolynomial` with `Rational = double`, not
this).

With the template parameter `Rational = double` — upstream's own documented
instantiation and the one the port implements — the `delta == 0` branch of
`SolveDepressedQuartic` divides by `9*c1^2 - 2*c2*a1`. That denominator is a
difference of two nearly equal products and can round to exactly zero, making
`root0` infinite and the derived `root1`, `root2` NaN. Those are then inserted
into `std::map<Real, int32_t>`, whose comparator `std::less<double>` is not a
strict weak ordering on NaN.

Explicit input (deep-run record 631 of `RootsPolynomial.solveQuartic`), the
quartic `x^3 * (-2.085220648210411 - 5.801142781619129*x)`:

```
p = (0, 0, 0, -2.085220648210411, -5.801142781619129)
```

MSVC's `lower_bound` finds neither `NaN < key` nor `key < NaN`, treats the NaN
as a duplicate of the key it landed on and drops the insert, so upstream
reports **one** root, `-infinity` of multiplicity 2. The port's explicit
ordered-array replacement for `std::map` appends both NaNs and reports
**three**. Neither answer is meaningful; the standard gives no behaviour at
all for the C++ one. The generator rejects records whose root map contains a
non-finite key, and the case comment records the input. A fix would need a
finiteness guard upstream does not have, so the port preserves its own
behaviour.

**`CubicRootsQR`/`QuarticRootsQR` lose all significance on badly scaled
companion matrices.** Characterised, not filed: the Francis iteration forms
`A[0][0]*A[0][0] + A[0][1]*A[1][0] - tr*A[0][0] + det`, so companion entries
of order `1e92` square to `1e184` and the uncoupling test
`tr + A[i][j] == tr` fires on the first iteration. For
`x^4 + 8.9e92 x^3 + 4.8e92 x^2 + 8.9e92 x + 5.8e92` upstream returns
`0.6666666666666667`, `-8.8e75`, `-8.9e92`, `0` after 0 iterations, of which
only `-8.9e92` is a root. The port reproduces it bit for bit. This is a
conditioning limitation with no exact separator, so it is preserved; the
reference check above reports the residuals with those two generator modes
excluded.
