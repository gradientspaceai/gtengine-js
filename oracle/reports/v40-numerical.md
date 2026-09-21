# Verify group 40 (numerical) against the MSVC build of upstream GTE

Family `v40-numerical`, 22 cases, 20 golden records each, deep run
`npm run oracle:deep -- 2000 v40-numerical` (44000 records, 4000 of them
throw records): **every case passes**, and the deviation case deviates as
required. Wall time of the deep run, generation and replay together: about
2 seconds.

Everything in this group is `+ - * / sqrt`: the closed-form 2x2/3x3/4x4
inverses, Gaussian elimination with full pivoting, the tridiagonal
recurrences, conjugate gradient (whose only `sqrt` is in the residual test)
and the Cholesky-based normal-equation solves of the two minimizers. The
residual functions `F` and the Jacobians `J` that drive the minimizers are
written in `oracle/cpp/cases/v40-numerical.cpp` and in
`test/oracle/v40-numerical.oracle.test.ts` as the same polynomial
expressions in the same association. **No case needs a tolerance**; all 20
non-throw cases are declared `{ exact: true }` and every real output is
bit-identical to the MSVC value.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `LinearSystem.h` | `solve.2x2`, `solve.3x3`, `solve.4x4` | exact | pass, 2000 records each |
| `LinearSystem.h` | `solve.nxn` (Gaussian elimination, Nx1 right-hand side), `solve.nxm` (NxM) | exact | pass |
| `LinearSystem.h` | `solveTridiagonal`, `solveConstantTridiagonal` | exact | pass |
| `LinearSystem.h` | `solveSymmetricCG.dense`, `solveSymmetricCG.sparse`, `solveSymmetricCG.zeroRHS` | exact | pass |
| `GaussNewtonMinimizer.h` | `minimize.rosenbrock`, `minimize.circleFit`, `minimize.jPlus`, `minimize.rankDeficient`, `minimize.zeroResidual`, `construct.invalidDimensions` | exact | pass |
| `LevenbergMarquardtMinimizer.h` | `minimize.rosenbrock`, `minimize.circleFit`, `minimize.jPlus`, `minimize.rankDeficient`, `construct.invalidDimensions` | exact | pass |
| `LevenbergMarquardtMinimizer.h` | `minimize.staleResidual.deviation` | deviation (#261) | 1857 of 2000 records deviate |

Both minimizer constructors are covered: the `JFunction` form (upstream
builds `J^T*J` with `MultiplyATB` and `-J^T*F` as `-(mF * mJ)`) and the
`JPlusFunction` form (the caller supplies both). Both `Result` structs are
compared field by field and bit for bit: `numIterations`, `numAdjustments`,
`converged`, `minError`, `minErrorDifference` (which is `DBL_MAX` on the
paths where no step ever reduced the error), `minUpdateLength` and every
component of `minLocation`.

### Branches the generators reach

Every matrix generator mixes uniform draws, a small integer lattice and an
aimed singular construction whose last row is a linear combination of the
earlier rows, and all three modes record the same number of doubles in
row-major order, so the replay does not need to know the mode.

* `solve.2x2` / `3x3` / `4x4` / `nxn` / `nxm`: the "not invertible" return is
  taken on 5-25% of records (`solve.2x2`: 561 of 2000).
* `solveTridiagonal` / `solveConstantTridiagonal`: mode 2 forces
  `diagonal[0] == 0` (the early return) and mode 3 constructs an exactly zero
  pivot in the recurrence; about half the deep-run records return `false`.
* `solveSymmetricCG`: the tolerance is drawn from four regimes (exactly 0,
  1e-8, 1e-3..1e-1 and 1e12) and `maxIterations` from 0..6, so the returned
  iteration count covers the immediate break at iteration 1, convergence
  inside the loop, and the exhausted loop that returns `maxIterations + 1`.
* Minimizers: `rankDeficient` makes `J^T*J` singular in two of three modes,
  so `CholeskyDecomposition::Factor` fails and the minimizer returns on its
  first iteration on 631 of 2000 records; `zeroResidual` starts at an exactly
  zero residual, which no later error can beat; `minimize.jPlus` for
  Levenberg-Marquardt reaches `numAdjustments` 1, 2 and 3 on 117 of 2000
  records, so the inner lambda-adjustment loop and the Gauss-Newton fallback
  after it are exercised with full agreement on the path where upstream is
  sound.

## Port defects fixed

None. Every arithmetic output of all 20 ordinary cases was bit-identical to
the MSVC build on the first run, over 44000 deep-run records. In particular
the accumulation orders that have produced defects in other groups are
already correct here: `LinearSystem`'s private `Dot` seeds with the literal
`0` (`src/LinearSystem.ts`, `dotArray`), `Mul` accumulates `P[row] +=
matA(row, col) * X[col]` in row-then-column order, `MultiplyATB` sweeps
`(r, c)` outside and the common index inside, `Cholesky` divides element by
element rather than multiplying by a reciprocal, and the sparse
conjugate-gradient path sorts its entries by `(row, col)` to reproduce
`std::map` iteration order.

The C++ overload question of ORACLE.md's triage list was checked
explicitly: `LinearSystem.h` includes `Matrix2x2.h`, `Matrix3x3.h` and
`Matrix4x4.h`, so `Inverse(A, &invertible)` inside `Solve` resolves to the
closed-form overloads, not to the Gaussian-elimination template of
`Matrix.h`. The port calls `inverse2x2` / `inverse3x3` / `inverse4x4`, which
matches, and the 2000-record agreement of `solve.2x2`, `solve.3x3` and
`solve.4x4` confirms it against the real build.

## Deliberate deviations demonstrated

`LevenbergMarquardtMinimizer::DoIteration` builds `-J^T*F` from the member
`mF`, which holds `F` at the previously *rejected* candidate whenever the
inner lambda-adjustment loop, or the Gauss-Newton fallback after it, runs
`DoIteration` more than once for the same `pCurrent`. The port re-evaluates
`F(pCurrent)` when the stored residual is stale (issue
[#261](https://github.com/gradientspaceai/gtengine-js/issues/261), recorded
as "fixed" in `docs/UPSTREAM-FINDINGS.md`).

A `DoIteration` is repeated exactly when the inner loop increments
`numAdjustments`, so the port and upstream agree exactly on the runs whose
`result.numAdjustments` is 0 for every prefix of the iteration. The C++ side
runs upstream itself with `maxIterations = 1, 2, ...` and records the longest
sound prefix, so `maxIterations` is an aimed input rather than a rejection
loop (`SoundLMIterations`, modelled on the `ApprCone3` cases of group 4).
The main cases `minimize.rosenbrock`, `minimize.circleFit` and
`minimize.rankDeficient` use that prefix and agree bit for bit on 2000
records each, which is the evidence that the fix is confined to the defective
inputs. `minimize.staleResidual.deviation` records the first prefix length at
which upstream repeats a `DoIteration`; 1857 of its 2000 records disagree
with the port, most of them already in the discrete `numIterations` /
`numAdjustments` outputs.

Two further entries of `docs/UPSTREAM-FINDINGS.md` for this group are
*preserved* defects and are therefore compared bit for bit rather than given
a deviation case:

* `SolveSymmetricCG` with `B = 0` computes `alpha = 0/0` and overwrites the
  exact solution `X = 0` with NaN. `solveSymmetricCG.zeroRHS` pins this:
  because the harness treats any NaN as equal to any NaN, the case emits the
  iteration count as an integer and the NaN pattern as booleans, so the
  record is not vacuous. All 2000 deep-run records produce an all-NaN `X` on
  both sides.
* `Result::numIterations` is the raw loop counter, so an exhausted loop
  reports `maxIterations + 1`; this is visible throughout the histograms
  above and is compared exactly.
* The Gauss-Newton fallback of `LevenbergMarquardtMinimizer` uses the
  maximally inflated lambda and never divides it back down. It is reached by
  `minimize.jPlus`, where it agrees bit for bit.

## Independent-reference checks

Agreement is not correctness (issue #507), so the deep-run outputs of the
main cases were checked once against references written from scratch (a
partial-pivot Gaussian elimination unrelated to GTE's full-pivot code, and
direct residual evaluation):

| check | records | median | maximum |
| --- | --- | --- | --- |
| `solve.2x2` relative `\|A X - B\|` | 1439 invertible | 3.7e-17 | 2.0e-15 |
| `solve.3x3` | 1765 | 7.6e-17 | 2.0e-14 |
| `solve.4x4` | 1842 | 9.7e-17 | 1.9e-14 |
| `solve.nxn` | 1767 | 6.3e-17 | 4.1e-16 |
| `solveTridiagonal` | 970 solved | 9.1e-17 | 1.9e-2 (see below) |
| `solveConstantTridiagonal` | 972 solved | 9.7e-17 | 6.1e-16 |
| `solveSymmetricCG.dense` vs the independent direct solve | 256 converged | 1.5e-16 | 4.4e-16 |

* Gauss-Newton on Rosenbrock reaches the known minimum `(b, b^2)`: over the
  861 records with `maxIterations >= 4` the worst deviation is 4.4e-16, and
  the median is exactly 0.
* No record of any minimizer case ends with `minError` larger than the error
  at the starting point; the circle-fit cases strictly reduce it on 1630 of
  2000 (Gauss-Newton) and 1878 of 2000 (Levenberg-Marquardt) records.

Two of the 970 solved `solveTridiagonal` records have a relative residual of
3.7e-3 and 1.9e-2. Both are lattice systems that are *exactly singular in
rational arithmetic* (for record 189, `diagonal = [3,2,3,3,4,1]`,
`subdiagonal = [2,-2,1,-2,0]`, `superdiagonal = [2,-1,-1,-1,1]`, whose third
pivot `3 - (-2)*(-1)/(2/3)` is zero in exact arithmetic) but whose pivot
rounds to a tiny nonzero value, so the pivot-free recurrence returns `true`
with a garbage solution. This is a property of the algorithm upstream
chooses, is reproduced bit for bit by the port, and is not a new finding: it
is the ordinary reason a tridiagonal solver without pivoting is documented as
requiring a diagonally dominant matrix. The lattice generator finding these
configurations is exactly what it is there for.

## Not covered

* **`Inverse`** has no entry point in any of this group's three headers; the
  closed-form `Inverse` of `Matrix2x2.h` / `Matrix3x3.h` / `Matrix4x4.h` is
  covered indirectly through `LinearSystem::Solve` (see above), and the
  Gaussian-elimination `Inverse` belongs to `GaussianElimination.h`, which is
  another group's header.
* **The column-major storage order.** Upstream selects it at compile time
  with `GTE_USE_COL_MAJOR`, and the oracle build does not define it, so only
  the row-major path (GTE's default) is compared. The port takes the order as
  a `rowMajor` argument; the column-major branch of `LexicoArray2` is covered
  by the port's own tests, not by this oracle.
* **`SolveTridiagonal` / `SolveConstantTridiagonal` with `N = 0`.** Upstream
  never validates `N`; `static_cast<size_t>(N) - 1` wraps to `SIZE_MAX` and
  `B[0]` is read out of bounds. That is undefined behaviour, so it stays out
  of the generators. The port asserts `N >= 1` instead (issue #261, "fixed
  (assert)"), a deviation that cannot be demonstrated against a build whose
  behaviour is undefined.
* **`GaussNewtonMinimizer` / `LevenbergMarquardtMinimizer` with a large
  number of residuals.** The cases keep `numF` at 3-8 so that a record stays
  a few tens of doubles, per ORACLE.md's size rule. Nothing in the two
  headers depends on `numF` beyond the loop bounds of `MultiplyATB` and the
  vector-matrix product.

## Upstream bug suspects

No new ones. Everything observed in this group is already recorded in
`docs/UPSTREAM-FINDINGS.md` under issue #261: the stale residual (fixed in
the port and demonstrated here), the `alpha = 0/0` NaN for a zero
right-hand side, the raw-loop-counter `numIterations`, the unconditional
`pCurrent = pNext` in Gauss-Newton that accepts error-increasing steps, the
never-reduced lambda of the Gauss-Newton fallback, and the unvalidated `N`
of the tridiagonal solvers.
