# Group 38 (numerical) — C++ oracle report

Family `v38-numerical`, 39 cases: 34 compared bit for bit (`exact: true`),
2 with a measured tolerance, 3 declared deviations. The deep run
(`npm run oracle:deep -- 2000 v38-numerical`, 78 000 records, 2.8 s wall)
passes; outside the deviation cases it compares 815 669 floating-point
outputs, 99.87 % of them bit-identical to the MSVC build, and the only
inexact ones come from the two libm cases.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `SymmetricEigensolver2x2.h` | `SymmetricEigensolver2x2.solve` | exact | 12 000 reals, 100 % |
| `SymmetricEigensolver3x3.h` | `SortEigenstuff.sort` | exact | 24 000 reals, 100 % |
| | `SymmetricEigensolver3x3.solve` | exact | 24 000 reals, 100 % |
| | `SymmetricEigensolver3x3.solve.outOfBandScale` | deviation (#379) | 2000 of 2000 records deviate |
| | `NISymmetricEigensolver3x3.solve` | tolerance 1e-12 (`std::acos`, `std::cos`) | 24 000 reals, 96.74 % exact, max scaled error 1.44e-15 |
| | `NISymmetricEigensolver3x3.solve.degenerate` | exact | 24 000 reals, 100 % |
| `UnsymmetricEigenvalues.h` | `UnsymmetricEigenvalues.solve` | exact | 4306 reals, 100 % |
| | `UnsymmetricEigenvalues.solve.trailingBlock` | deviation (#42) | 2000 of 2000 records deviate |
| | `UnsymmetricEigenvalues.solve.invalidSize` | exact | discrete only |
| | `UnsymmetricEigenvalues.solve.cycling` | exact | discrete only |
| | `UnsymmetricEigenvalues.solve.nonConvergence` | exact | discrete only |
| `BandedMatrix.h` | `BandedMatrix.accessor` | exact | 101 918 reals, 100 % |
| | `BandedMatrix.choleskyFactor` | exact | 24 302 reals, 100 % |
| | `BandedMatrix.solveSystem.vector` | exact | 32 054 reals, 100 % |
| | `BandedMatrix.solveSystem.matrix` | exact | 40 395 reals, 100 % |
| | `BandedMatrix.computeInverse` | exact | 53 567 reals, 100 % |
| | `BandedMatrix.computeInverse.zeroLeadingMinor` | exact | 72 000 reals, 100 % |
| | `BandedMatrix.invalidShape` | exact | 18 000 reals, 100 % |
| `GaussianElimination.h` | `GaussianElimination.compute` | exact | 41 731 reals, 100 % |
| | `GaussianElimination.compute.partial` | exact | 11 925 reals, 100 % |
| | `GaussianElimination.compute.subnormal` | exact | 21 059 reals, 100 % |
| | `GaussianElimination.compute.invalidInput` | exact, throw parity | 1398 of 2000 records throw on both sides |
| `Integration.h` | `Integration.trapezoidRule` | exact | 2000 reals, 100 % |
| | `Integration.romberg` | exact | 2000 reals, 100 % |
| | `Integration.computeQuadratureInfo` | exact | 22 038 reals, 100 % |
| | `Integration.gaussianQuadrature` | exact | 23 950 reals, 100 % |
| | `Integration.libmIntegrand` | tolerance 1e-12 (`std::exp`) | 4000 reals, 92.3 % exact, max scaled error 6.42e-16 |
| `LCPSolver.h` | `LCPSolver.solve` | exact | 13 970 reals, 100 % |
| | `LCPSolver.solve.positiveDefinite` | exact | 13 934 reals, 100 % |
| | `LCPSolver.solve.trivial` | exact | 13 682 reals, 100 % |
| | `LCPSolver.solve.noSolution` | exact | 13 986 reals, 100 % |
| | `LCPSolver.solve.maxIterations` | exact | 18 172 reals, 100 % |
| | `LCPSolver.solve.roundoffPivot` | exact | 24 000 reals, 100 % |
| | `LCPSolver.solve.invalidInput` | exact | 10 078 reals, 100 % |
| `FPInterval.h` | `FPInterval.leaf` | exact | 16 000 reals, 100 % |
| | `FPInterval.internal` | exact | 48 000 reals, 100 % |
| | `FPInterval.operators` | exact | 56 000 reals, 100 % |
| | `FPInterval.productBounds` | exact | 4000 reals, 100 % |
| | `FPInterval.directedRounding` | deviation (emulated rounding) | 2000 of 2000 records deviate, max 1 ulp |
| `OdeSolver.h` | — | — | see "Not covered" |

Generator coverage (deep-run histograms). The symmetric 3x3 generator spreads
its nine modes evenly (214–245 records each: uniform, lattice, diagonal,
already tridiagonal, a multiple of the identity, rank one, the zero matrix, a
lattice matrix scaled by 2^k with |k| <= 500, signed zeros); the iterative
solver's iteration count ranges over 0, 4–9 and 11–16. `UnsymmetricEigenvalues`
covers n = 3…6, reports 0–6 eigenvalues and exhausts its budget on 579 of 2000
records. `BandedMatrix.choleskyFactor` succeeds on 454 of 2000 and
`computeInverse` on 1584 of 2000; `GaussianElimination` reports invertible on
1279 of 2000. The LCP cases reach every `Result` code: 557
`HAS_TRIVIAL_SOLUTION`, 365 `HAS_NONTRIVIAL_SOLUTION`, 1042 `NO_SOLUTION` and
36 `FAILED_TO_CONVERGE` in the general case, plus 2000 `INVALID_INPUT` in its
own case.

## FENV_ACCESS: how the C++ side really exercises the rounding mode

`FPInterval.h` switches the hardware rounding mode with `std::fesetround`
around every arithmetic operation. Measured on this toolchain (MSVC 19.44,
x64, `/O2`):

* with plain `/fp:precise` the compiler common-subexpression-eliminates the
  two evaluations of `u + v` that straddle the `fesetround` calls, so both
  endpoints come back as the round-to-nearest value and **upstream's intervals
  are degenerate**: `Add`, `Sub`, `Mul`, `Div` and `Reciprocal` all return
  `[x, x]`;
* `#pragma fenv_access (on)`, which `/fp:precise` accepts, restores the
  intended behaviour for all of them except one: `1.0/3.0` becomes
  `[3fd5555555555555, 3fd5555555555556]`.

The case file therefore carries `#pragma fenv_access (on)`. It was verified
by regenerating the goldens with and without the pragma that **only the
`FPInterval` cases change**; every other case in the file is byte-identical,
so the pragma does not perturb the rest of the group.

The one exception: even with the pragma, `/O2` still CSEs the two divisions
of the two-argument leaf `Div(u, v)` (whose bodies are the textually identical
`u / v`), so that overload alone returns a degenerate interval where
`/fp:strict` and `/Od` return the correct directed pair. The `FPInterval.leaf`
case divides by exact powers of two, where the quotient is exact and the CSE
is invisible; the four-argument `Div(u0, u1, v0, v1)` is unaffected and is
covered by `FPInterval.internal` and `FPInterval.directedRounding`.

## Port defects fixed

**`src/FPInterval.ts`, the emulated round toward negative infinity signed an
exact zero wrongly.** JavaScript cannot set the rounding mode, so the port
computes the round-to-nearest value, proves exactness with TwoSum/TwoProduct
where it can, and otherwise steps one representable value outward. The
exactness proof was complete for the magnitude but not for the sign: under
`FE_DOWNWARD` an exactly zero sum or difference is `-0`, not `+0` (IEEE
754-2019 6.3), and the only exception is the sum of two `+0` operands. The
MSVC build confirms it directly: with `FE_DOWNWARD`, `1 - 1`, `(+0) + (-0)`,
`(+0) - (+0)` and `(-0) - (-0)` all give `-0` while `(+0) + (+0)` gives `+0`,
and multiplication and division keep the exclusive-or sign rule in every
mode. The port returned `+0`, so the lower endpoint of `[1,1] - [1,1]`
differed from the C++ build in its sign bit. Error size: a sign bit on an
exact zero, which `Object.is` and the oracle harness both distinguish. Fixed
by `signedZeroDown` in `addDown`/`subDown`; two existing expectations in
`test/FPInterval.test.ts` that asserted `+0` were updated, and a regression
test covers all six sign combinations plus the multiplication and division
cases that must NOT change.

**`src/FPInterval.ts`, the emulated rounding widened infinities that are
exact results.** `1/0`, `x * infinity` and `infinity + x` are produced by
IEEE without rounding, so every rounding mode returns them unchanged; only an
*overflow* of finite operands is a rounded result. The port's exactness tests
answered "not proven exact" for every non-finite result, so
`FPInterval.reciprocal(-4, 0)` returned `MAX_VALUE` where upstream's
`FE_DOWNWARD` returns `+infinity`. Error size: one endpoint,
`1.7976931348623157e308` versus `Infinity`. `FPInterval` produces infinite
endpoints routinely through `reciprocalDown`, `reciprocalUp` and `reals`, and
those flow through `mul`, so this is reached on ordinary inputs. Fixed by
answering "exact" in `sumIsExact`, `productIsExact` and `quotientIsExact`
whenever an operand is non-finite (or the divisor is zero), which leaves the
overflow widening — the thing that reproduces upstream's `MAX_VALUE` lower
bound for `MAX_VALUE + MAX_VALUE` — untouched. Both behaviours are pinned by
regression tests.

No other port defect was found: every other case is bit-identical on every
one of its 78 000 deep-run records, including the eigensolver iteration
counts, the LCP pivot sequences and iteration counts, and the Cholesky and
Gaussian-elimination factorizations.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `SymmetricEigensolver3x3.solve.outOfBandScale` | [#379](https://github.com/gradientspaceai/gtengine-js/issues/379) | 2000 of 2000 records. Upstream's `GetCosSin` evaluates `sqrt(u*u + v*v)` with no rescaling, so for a lattice matrix scaled by 2^k with 520 <= \|k\| <= 1000 the squares overflow or underflow and upstream returns `(0,0)` instead of a unit 2-tuple; only 4.9 % of its real outputs still agree, and the worst scaled difference is 1.83 (upstream's eigenvectors have length 0). The main case stays inside \|k\| <= 500, where the port's rescaling is inactive and evaluates the upstream expression unchanged — 24 000 outputs, all bit-identical. |
| `UnsymmetricEigenvalues.solve.trailingBlock` | [#42](https://github.com/gradientspaceai/gtengine-js/issues/42) | 2000 of 2000 records. Upstream's eigenvalue-packing loop runs `i < N-1` and never reports `A(N-1,N-1)` when the final subdiagonal has decoupled, so the port reports one eigenvalue more and the counts differ structurally. The separator is upstream's own `mSubdiagonalFlag[N-2]`; because it is private, the case file carries `ProbeUnsymmetric`, a verbatim copy of `House`, `RowHouse`, `ColHouse`, `ReduceToUpperHessenberg`, `FrancisQRStep`, `GetBlock` and the iteration loop of `Solve`, and both the main case (flag 1) and the deviation case (flag 0) select on it. |
| `FPInterval.directedRounding` | the PORT DEVIATION note at the top of `src/FPInterval.ts` | 2000 of 2000 records; 33 % of the 48 000 endpoints differ, every one of them by exactly one ulp (max scaled error 4.43e-16). Upstream takes directed rounding from `std::fesetround`; the port widens whenever it cannot prove the operation exact, so its interval contains upstream's. The four exact `FPInterval` cases use dyadic operands on which every operation is exactly representable and are bit-identical there, which is what tests the branch structure and the exactness proofs. |

The preserved upstream defects of these headers are simply compared bit for
bit and are reached on purpose:
[#375](https://github.com/gradientspaceai/gtengine-js/issues/375)
(`GaussianElimination.compute.subnormal`: a matrix of subnormal magnitude is
reported invertible with a NaN inverse — the `invertible` boolean keeps the
record from being vacuous when every real output is NaN),
[#75](https://github.com/gradientspaceai/gtengine-js/issues/75)
(`FPInterval.productBounds` deliberately draws the both-straddle-zero branch
where `ProductLowerBound` returns a positive lower bound), and
[#476](https://github.com/gradientspaceai/gtengine-js/issues/476)
(`UnsymmetricEigenvalues.solve.cycling` and `.nonConvergence`: `FrancisQRStep`
has no exceptional shift, so all three inputs burn their whole budget and
report no eigenvalues on both sides).

## Independent-reference checks (agreement is not correctness)

Every main case's 2000-record deep-run output was checked once against a
reference computed outside both implementations. All of these run on the C++
outputs, which the port reproduces bit for bit, so they check both sides.

| check | result |
| --- | --- |
| `SymmetricEigensolver3x3`: `A v = lambda v`, `V V^T = I`, `det V = 1` | worst scaled residual 1.37e-15, worst \|V V^T - I\| 1.11e-15, worst \|det V - 1\| 1.33e-15 (279 records skipped: the zero matrix and the scaled modes whose entries overflow the residual product) |
| `NISymmetricEigensolver3x3`: same | worst 3.25e-15 / 1.33e-15 / 1.33e-15, nothing skipped |
| `SymmetricEigensolver2x2`: same in 2D | worst scaled residual 4.72e-16, worst \|det V - 1\| 2.22e-16 |
| `UnsymmetricEigenvalues`: backward error `sigma_min(A - lambda I) / \|\|A\|\|_F`, by a one-sided Jacobi SVD (forming `A^T A` floors the measurement at sqrt(eps) and gives a misleading 1e-8) | worst 5.59e-15 over 4306 eigenvalues of the main case and 1.69e-15 over 5798 of the deviation case: every reported eigenvalue is an exact eigenvalue of a matrix within 6e-15 of the input |
| `GaussianElimination`: `M M^-1 = I`, and the determinant against the **exact** integer determinant (Bareiss over `BigInt`) on the 623 lattice records | worst \|M M^-1 - I\| 3.43e-16; worst determinant error 3.7e-13 relative to `max(1, \|det\|)`, consistent with a product of five pivots of intermediate magnitude ~4^5; no record reported singular that was not exactly singular |
| `LCPSolver`: `w = q + M z`, `w >= 0`, `z >= 0`, `w . z = 0` | main case 0 violations of 922 solved records (worst residual 1.18e-15); positive-definite case 0 of 1837; trivial case 0 of 2000; infeasible case **1 violation of 1 solved record** — see the suspects below |
| `LCPSolver`, positive definite `M`: the LCP is always solvable, yet 163 of 2000 records return `FAILED_TO_CONVERGE` | not a defect: all 163 solve, with valid solutions, when the iteration budget is raised from the default `n*n` to 20 000. The header calls the `n*n` default "chosen arbitrarily". |
| `Integration.gaussianQuadrature` against the closed-form polynomial integral, restricted to the 964 records where the rule is exact (`2*degree - 1 >= polynomial degree`) | worst relative error 1.31e-12 |
| `Integration.romberg` against the closed form, restricted to the 637 records with `order >= degree + 2` | worst relative error 6.77e-15 |
| `Integration.trapezoidRule` on linear integrands, where the rule is exact | worst relative error 3.44e-15 |
| `BandedMatrix.solveSystem.vector`: residual `A x - b` against the recorded matrix | mode 2 (symmetric positive definite, the only inputs the Cholesky routine is defined for) worst 3.46e-15 over 400 records; mode 3 (symmetric indefinite) 6.41e-16. The six mode-0 records that pass `CholeskyFactor` with a non-symmetric matrix have residuals up to 5.4, which is a precondition violation of the routine, not a defect. |
| `BandedMatrix.computeInverse`: `A A^-1 = I`, and for the 1605 integer records the exact rational inverse (`BigInt` fractions) | 1603 correct; **2 wrong** — see the suspects below |

## Not covered

| header / entry point | reason |
| --- | --- |
| `OdeSolver.h` | Abstract base class. Its only members are the trivial accessors `SetTDelta`/`GetTDelta` and the pure virtual `Update`; it has no computational entry point of its own, and the concrete solvers that implement `Update` (`OdeEuler`, `OdeImplicitEuler`, `OdeMidpoint`, `OdeRungeKutta4`) belong to other verify groups. |
| `FPInterval::operator+=`, `-=`, `*=`, `/=` | Not ported: intervals are immutable in the port, and upstream defines each as `u = u OP v`, so they add no arithmetic over the binary operators, which are covered. |
| `FPInterval<float>` | The port maps C++ floating point to IEEE binary64 only (PORTING.md). |
| `FPInterval::operator=` and the copy constructor | Value semantics in C++; the port's `clone()` is a plain field copy with no arithmetic. The other four constructors are covered by `FPInterval.internal`. |
| `GaussianElimination` with `GTE_USE_COL_MAJOR` | The port defaults to row major, which is the GTE default, and the group brief restricts this header to it. The column-major path differs only in `LexicoArray2` indexing, which group 26 covers. |
| `GaussianElimination`'s "`B` without `X`" and "`C` without `Y`" preconditions | Upstream signals "compute this output" with a pointer pair and raises `LogError` when only one of the pair is non-null. The port takes an options object, where the pair cannot be half-supplied. The two preconditions the port can express, `numRows <= 0` and a `C` with `numCols < 1`, are covered with throw parity. |
| `LCPSolver<T, n>` (compile-time dimension) | No TypeScript analogue; the port is upstream's dynamic `LCPSolver<T>`. Both derived classes only allocate storage and delegate to the identical `LCPSolverShared<T>::Solve`, which is what the cases exercise. |
| `LCPSolverShared(int32_t n, T const& zero, T const& one)` | The constructor exists to select the right zero and one for `QFNumber`/rational instantiations; only the floating-point instantiation is ported. |
| `LCPSolver` constructed with `n <= 0` ([#76](https://github.com/gradientspaceai/gtengine-js/issues/76), fixed in the port) | Cannot be compared against the C++ build: upstream leaves its member pointers null and `Solve` dereferences them, which is an access violation (`0xC0000409`-class process death), not a catchable `std::exception`. The oracle driver would not survive the record. The port's `INVALID_INPUT` return is covered by `test/LCPSolver.test.ts`. |
| `BandedMatrix`'s `const` accessor overloads | Identical bodies to the non-`const` ones; the port has a single `get`/`set` pair, which the `BandedMatrix.accessor` case exercises over the whole square plus the out-of-range probes. |

## Upstream bug suspects

Two new ones. In both, the port and the MSVC build agree bit for bit and both
are wrong; the port preserves the behaviour, and each has a pinned regression
test and a targeted oracle case so the goldens carry the reproduction.

**1. `BandedMatrix::ComputeInverse` returns a wrong inverse, and `true`, when
a leading principal minor vanishes (result-corrupting).** The routine
eliminates with **no pivoting** and rejects a matrix only on `diag == 0`.
When a leading principal minor of the banded matrix is exactly zero, the
pivot is zero in exact arithmetic but a round-off residue in floating point
(the elimination has already divided by earlier pivots, so the entries are no
longer integers and the cancellation is inexact), the test passes, the
multiplier is about 1e16 and the result is garbage with a `true` return.

*Reproduction* (deep-run record 792 of `BandedMatrix.computeInverse`; lower
bandwidth 2, upper bandwidth 4):

```
 1   0   3  -3   3   0
 4   2  -2  -1  -1   3
 0  -4   0  -2  -2   3
 0   2   0   1   4  -2
 0   0   4  -3   2  -3
 0   0   0  -2   4  -2
```

Its determinant is -656 and its fourth leading principal minor is 0. The
reported first row of the inverse is
`(-0.7434619, 0.4358654, 0.7960934, 1.1295357, 0.7889220, -0.4839663)`; the
exact inverse's first row is `(-31, 18, 33, 48, 129/4, -159/8)/41`, that is
`(-0.7560975, 0.4390243, 0.8048780, 1.1707317, 0.7865853, -0.4847560)`: a
relative error of 4.4e-2, wrong in the second digit. Dense
`GaussianElimination` with full pivoting gets the same matrix right to 1e-12,
so the loss is the missing pivoting, not the conditioning. A second record
(665, a 5x5 with bands 3/3) is **exactly singular** and still gets
`ComputeInverse` returning `true`. Frequency: 2 of the 1605 integer records
of the deep run; 418 of them have some vanishing leading minor, and the other
416 are caught because their pivot is still exactly zero in floating point.

*Port: preserved*, by the same judgment already recorded for
`GaussianElimination`'s denormal pivots
([#375](https://github.com/gradientspaceai/gtengine-js/issues/375)). There is
no exact condition separating "this pivot is a round-off residue of an exact
zero" from "this pivot is legitimately small", so replacing `diag == 0` with
a tolerance would turn correct inverses into "not invertible"; and adding
pivoting would change the result of every input and destroy the band
structure the routine relies on. Pinned by
`test/BandedMatrix.test.ts`, "upstream (preserved): computeInverse returns a
wrong inverse with a true return when a leading principal minor vanishes",
and by the oracle case `BandedMatrix.computeInverse.zeroLeadingMinor`.

**2. `LCPSolverShared::Solve` reports a solution for a provably infeasible
LCP when a round-off-sized pivot passes the ratio test (result-corrupting).**
The header's own comment anticipates the mechanism — "it is possible that
theoretically `mAugmented[r][driving]` is zero but rounding errors cause it
to be slightly negative" — and says the hoped-for outcome is a failure to
converge, hence `NO_SOLUTION`. It is not always.

*Reproduction* (deep-run record 501 of `LCPSolver.solve.noSolution`), all
integers, `q < 0` and `M <= 0` entrywise, so `w = q + M z <= q < 0` for every
`z >= 0` and the problem has no solution:

```
q = (-5, -3, -2, -4, -1, -5)
M = [[-1,-2,-2, 0,-3,-3], [-2, 0,-3,-1,-2,-3], [-3,-3,-2, 0,-2,-2],
     [-2,-3, 0,-1,-2,-2], [-3, 0,-2,-2, 0,-3], [ 0,-3,-2,-2,-3,-2]]
```

An instrumented replay of the pivoting shows the failure exactly: iterations
0 to 5 use pivots of magnitude 1 to 3 and keep the dictionary consistent;
iteration 6 accepts `Augmented(r, driving) = -1.1102230246251565e-16`, the
round-off residue of an exactly zero entry, whose reciprocal is
-9.0e15, and the dictionary residual jumps to 5.4e16. After a few more pivots
the artificial variable happens to leave the basis, so at iteration 19 the
solver returns `HAS_NONTRIVIAL_SOLUTION` with
`w = (0,0,0,0,0,2)`, `z = (16, 4, 2/9, 16, 4, 0)` and
`|w - q - M z|` between 41 and 81. Raising the budget to 4096 does not help:
the damage is done by the pivot, not the budget. Frequency: 1 of 2000 records
of that case; the same generator at `n <= 5` (the existing property test in
`test/LCPSolver.test.ts`) never reproduces it.

This is the ordinary-magnitude version of the subnormal-pivot failure
recorded as [#476](https://github.com/gradientspaceai/gtengine-js/issues/476)
("a subnormal pivot turns a provably infeasible problem into
`HAS_NONTRIVIAL_SOLUTION` with infinite `z`"), and it is a stronger statement:
no subnormals and no extreme magnitudes are needed, only small integers.

*Port: preserved.* A guard on the pivot magnitude would be a tolerance with
no principled threshold, and it would change the pivot sequence of every
problem. Pinned by `test/LCPSolver.test.ts`, "upstream (preserved): a
round-off-sized pivot turns a provably infeasible LCP into a reported
solution", and by the oracle case `LCPSolver.solve.roundoffPivot`, which
reproduces it at every power-of-two scale (11 of its 20 golden records return
the bogus `HAS_NONTRIVIAL_SOLUTION`, the other 9 run out of the drawn budget
first).

**Not a defect, for the record.** The `FE_DOWNWARD` common-subexpression
elimination described above is an MSVC codegen limitation, not an upstream
bug, and `FPInterval.h` has no way to defend against it; callers who need
rigorous bounds must compile with `/fp:strict` or
`#pragma fenv_access (on)`. And the 163 positive-definite LCPs that return
`FAILED_TO_CONVERGE` are the arbitrary `n*n` default budget, not a failure of
the method: every one of them solves correctly with a larger budget.
