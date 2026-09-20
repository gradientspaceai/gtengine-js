# Group 3 (`v03-approximation`) — C++ oracle report

75 cases, 20 records each in `oracle/golden/v03-approximation.txt` (1500
records). 70 of the cases are ordinary comparisons and 5 are deliberate
`deviation` cases; 45 of the golden records are recorded C++ assertion
failures that the port reproduces (the `ApprPolynomialSpecial*` constructor
asserts and `ApproximateCurveByArcs`'s `numArcs >= 1`). Every ordinary case is
declared `{ exact: true }`; they carry 10644 floating-point outputs in the
goldens, **all 10644 bit-identical** to the MSVC build of upstream GTE. No
case needed a tolerance.

Deep run `npm run oracle:deep -- 2000 v03-approximation`: 75 cases, 150000
records, 4735 recorded C++ exceptions, **1 064 007 floating-point outputs
compared over the 70 ordinary cases, 100.000 % bit-identical** (worst scaled
error 0), and every deviation case deviating. Wall time 40 s (generation plus
replay). All discrete outputs (the success flags, `GetMinimumRequired`, the
parameter-array sizes, the iteration counts of `FitUsingLengths`, the arc and
time counts of `ApproximateCurveByArcs`, the RANSAC consensus indices) are
compared exactly and agreed on every record of every ordinary case.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise`. No `src/` file was changed by this group.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `ApprQuery.h` | `fit.all`, `fit.range`, `fit.indexed`, `fit.numIndices`, `fitIndexed.direct`, `ransac.tooFew`, `ransac.minimum` | exact | pass |
| `ApprCircle2.h` | `fitUsingSquaredLengths`, `fitUsingSquaredLengths.cocircular`, `fitUsingLengths`, `fitUsingLengths.epsilon` | exact | pass |
| `ApprGaussian2.h` | `fit`, `fitIndexed`, `copyParameters` | exact | pass |
| `ApprGaussian3.h` | `fit`, `fitIndexed`, `copyParameters` | exact | pass |
| `ApprHeightLine2.h` | `fit`, `fit.verticalData`, `copyParameters` | exact | pass |
| `ApprHeightPlane3.h` | `fit`, `fit.collinearXY`, `copyParameters` | exact | pass |
| `ApprOrthogonalLine2.h` | `fit`, `fitIndexed`, `copyParameters` | exact | pass |
| `ApprOrthogonalLine3.h` | `fit`, `fitIndexed`, `copyParameters` | exact | pass |
| `ApprOrthogonalPlane3.h` | `fit`, `fitIndexed`, `copyParameters` | exact | pass |
| `ApprParallelLines2.h` | `fit`, `fit.deviation` | exact (+1 deviation) | pass |
| `ApprSphere3.h` | `fitUsingSquaredLengths`, `fitUsingSquaredLengths.cospherical`, `fitUsingLengths`, `fitUsingLengths.epsilon` | exact | pass |
| `ApprCurveByArcs.h` | `compute`, `compute.collinear`, `compute.epsilon`, `compute.assert` | exact | pass |
| `ApprPolynomial2.h` | `fit`, `fitIndexed`, `copyParameters` | exact | pass |
| `ApprPolynomial3.h` | `fit`, `fitIndexed`, `copyParameters` | exact | pass |
| `ApprPolynomial4.h` | `fit`, `fitIndexed`, `copyParameters` | exact | pass |
| `ApprPolynomialSpecial2.h` | `fit`, `fitIndexed`, `copyParameters`, `constructor.assert` | exact | pass |
| `ApprPolynomialSpecial3.h` | `fit`, `fitIndexed`, `copyParameters`, `constructor.affineAssert` | exact | pass |
| `ApprPolynomialSpecial4.h` | `fit`, `fitIndexed`, `copyParameters`, `constructor.sizeAssert` | exact | pass |
| `ApprQuadratic2.h` | `ApprQuadratic2.fit`, `.fit.cocircular`, `.fit.decoupledDeviation`, `ApprQuadraticCircle2.fit`, `.fit.cocircular`, `.fit.decoupledDeviation` | exact (+2 deviations) | pass |
| `ApprQuadratic3.h` | `ApprQuadratic3.fit`, `.fit.cospherical`, `.fit.decoupledDeviation`, `ApprQuadraticSphere3.fit`, `.fit.cospherical`, `.fit.decoupledDeviation` | exact (+2 deviations) | pass |

Every computational public entry point the port implements is covered: all
four `ApprQuery::Fit` overloads and `FitIndexed` (the `(count, pointer)` and
`std::vector` forms of `Fit(observations)` build the same identity index list
and collapse to one function in the port), `GetMinimumRequired`, `Error`,
`CopyParameters`, the two deterministic `RANSAC` paths, both `FitUsing*`
variants of the circle and the sphere with and without the
`initialCenterIsAverage` option and with and without a positive `epsilon`,
`GetParameters` and every domain accessor of the polynomial fitters,
`Evaluate` for each of them, `ApproximateCurveByArcs`, and both classes of
each `ApprQuadratic` header.

### Why every case is exact

No fitting path in this group reaches the C math library except through
`std::sqrt`. `SymmetricEigensolver2x2`, the **iterative**
`SymmetricEigensolver3x3` (which `ApprGaussian3`, `ApprOrthogonalLine3` and
`ApprOrthogonalPlane3` call — the non-iterative `NISymmetricEigensolver3x3`,
which uses `std::acos` and `std::cos`, is a different class that no fitter in
this group calls), the NxN `SymmetricEigensolver` (Householder plus
implicit-shift QL, behind `ApprQuadratic*`), `GaussianElimination` (behind
`Inverse(GMatrix)` in the polynomial fitters) and
`RootsPolynomial::Find` (Cauchy bound plus recursive bisection, behind
`ApprParallelLines2`) use only `+ - * /`, `sqrt` and `fabs`. The
`ApprCurveByArcs` path adds Romberg integration of `Length(X'(t))` and a
bisection for `GetTime`, again arithmetic and `sqrt` only.

### Generators and the branches they reach

Every point-set helper alternates a uniform cloud with exactly representable
lattice sets and with constructions aimed at a specific branch, and records
the count, the generator mode and then the coordinates. The distributions
below are from the 2000-record deep run.

- **2D point sets** (`MakePoints2`): uniform; a small integer lattice; all
  points coincident; collinear on an integer lattice line; and exactly
  cocircular (the twelve lattice points of the radius-5 circle about an
  integer centre). **3D point sets** (`MakePoints3`) add a coplanar lattice
  mode and use the 30 lattice points of the radius-5 sphere.
- **Rank-deficient branches.** `ApprCircle2.fitUsingSquaredLengths` reports
  `det == 0` on 969 of 2000 records and `ApprSphere3` on 1026;
  `ApprHeightLine2.fit` reaches its `covar00 > 0` failure branch on 629
  records (the dedicated `fit.verticalData` case makes `covar00` exactly zero
  on every record by giving every sample the same abscissa) and
  `ApprHeightPlane3.fit` its `det != 0` failure branch on 942 (the
  `fit.collinearXY` case makes the 2x2 determinant exactly zero by placing the
  (x,y) projections on a lattice line).
- **Eigenvalue ties.** The uniqueness flags of the orthogonal fitters are
  `eval[0] < eval[1]` (2D) and `eval[1] < eval[2]` / `eval[0] < eval[1]` (3D),
  so the coincident and lattice modes, which produce exactly equal
  eigenvalues, reach the `false` branch: 577 of 2000 for
  `ApprOrthogonalLine2`, 516 for `ApprOrthogonalLine3` and 574 for
  `ApprOrthogonalPlane3`.
- **Iteration counts.** `ApprCircle2.fitUsingLengths` draws `maxIterations` in
  `[0,6]` and reaches every count from 0 to 6; the `.epsilon` variant draws
  `maxIterations` in `[1,12]` and a positive `epsilon`, so the convergence
  `break` fires and the counts spread over 0 to 12 (665 records stop at 0,
  that is, on the first iteration). `ApprSphere3` behaves the same way.
- **Index plumbing.** `MakeIndices` produces identity, reversed, repeated
  (`i % 2`) and random index lists, so every `FitIndexed` case sees
  out-of-order and duplicated indices. `ApprQuery.fit.range` draws `imin` and
  `imax` independently, so the `imin <= imax` rejection fires on about 44 % of
  the records (583 of 2000 fits succeed); `ApprQuery.fit.numIndices` draws
  `numIndices` up to `indices.size() + 2`, which exercises the
  `std::min(numIndices, indices.size())` clamp.
- **Polynomial fitters.** The observation helpers add a degenerate-x mode
  (every sample shares its abscissa), a degenerate-w mode, and an exact
  polynomial-sample mode. `ApprPolynomial2/3/4` return `hasNonzero == false`
  on 684, 499 and 639 of 2000 records (a singular Vandermonde matrix, which
  `Inverse` reports by returning the zero matrix). The special-polynomial
  fitters divide by the sample range in `Transform` with no guard, so the
  degenerate-domain modes produce `Infinity` scales and NaN coefficients on
  928, 927 and 863 of 2000 records — the documented upstream defect, which
  the port preserves and which agrees bit for bit here. `MakeDegrees` builds
  strictly increasing degree lists with gaps of 1 to 3, so the term sets are
  genuinely sparse.
- **Curve by arcs.** The curve is a `BezierCurve<2,double>` of degree 2 to 4
  with lattice control points. With `epsilon = 0` the test
  `fabs(det) >= epsilon` is always true, so `compute` takes the division
  branch on all 3950 arcs — that is upstream's issue #163 defect, preserved.
  `compute.collinear` puts the control points on a lattice line and draws
  `epsilon` positive on half its records: 1970 of 4010 arcs take the
  `numeric_limits::max()` line-segment sentinel. `compute.epsilon` draws
  `epsilon` in `[0,20]` on ordinary curves and reaches the sentinel on 4341 of
  4913 arcs. `compute.assert` draws `numArcs` in `{0,1}`, so half its records
  are the `LogAssert` throw-parity check.
- **Quadratic fits.** The cocircular and cospherical constructions put the
  samples exactly on a conic or quadric, which drives the minimum eigenvalue
  to the round-off floor and makes the `std::max(eigenvalue, 0)` clamp
  reachable.

## Port defects fixed

None. The port's transcription of all 20 headers is bit-exact against the
MSVC build on every one of the 1 064 007 floating-point outputs of the deep
run. Every disagreement the oracle found was a deliberate, documented port fix
of an upstream defect; the two of them are demonstrated by the five
`deviation` cases below.

## Agreement is not correctness: independent-reference spot checks

Bit-identity only says the port does what upstream does. On the 2000-record
deep run the following outputs were additionally checked against a reference
that does not use the fitted result at all. The reference values are computed
from the recorded inputs, so they check upstream and the port at once.

- **`ApprCircle2.fitUsingSquaredLengths.cocircular`.** The samples lie exactly
  on the radius-5 circle about an integer centre, so the least-squares circle
  must be that circle. The integer centre is recovered from the samples by
  exhaustive search (the only lattice point at squared distance exactly 25
  from every sample). On the 1920 records whose determinant is nonzero the
  worst centre error is **2.1e-14** and the worst `|radius - 5|` is
  **1.6e-14**. The 80 skipped records are the ones where the determinant is
  zero or the samples collapse to fewer than three distinct points.
- **`ApprSphere3.fitUsingSquaredLengths.cospherical`.** The same check in 3D,
  restricted to the 1764 records whose sample covariance has determinant above
  20 (otherwise the samples are coplanar and the sphere centre is not
  determined by them at all, even though upstream's determinant is only
  numerically, not exactly, zero). Worst centre error **1.9e-14**, worst
  `|radius - 5|` **9.8e-15**.
- **`ApprPolynomial2.fit`.** On the 48 records whose generator mode samples an
  exact cubic on a lattice, with `degree == 3` and at least four distinct
  abscissae, the fitted polynomial must reproduce every sample. Worst relative
  residual **2.7e-13**.
- **`ApprOrthogonalLine3.fit`.** On the 248 records whose samples are collinear
  on a lattice line, the fitted direction must be parallel to that line.
  Worst `|d x D| / |d|` is **2.4e-16`.

These checks were run once from the deep-run output with a throwaway script;
they are not part of the committed test suite (the invariant-style property
tests of the verification wave already cover them continuously).

## Deliberate deviations demonstrated

| case | records disagreeing (golden / deep) | record of the decision |
| --- | --- | --- |
| `ApprParallelLines2.fit.deviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `ApprParallelLines2.h` items 1-4; issue [#91](https://github.com/gradientspaceai/gtengine-js/issues/91) |
| `ApprQuadratic2.fit.decoupledDeviation` | 5/20, 579/2000 | `docs/UPSTREAM-FINDINGS.md` `SymmetricEigensolver.h` `Tridiagonalize`; issue [#80](https://github.com/gradientspaceai/gtengine-js/issues/80) |
| `ApprQuadratic3.fit.decoupledDeviation` | 3/20, 282/2000 | the same |
| `ApprQuadraticCircle2.fit.decoupledDeviation` | 3/20, 79/2000 | the same |
| `ApprQuadraticSphere3.fit.decoupledDeviation` | 9/20, 349/2000 | the same |

### `SymmetricEigensolver::Tridiagonalize`, reached through `ApprQuadratic*`

All four `ApprQuadratic` entry points minimise `C^T M C` over unit `C` by
taking the eigenvector of the smallest eigenvalue of `M` from the NxN
`SymmetricEigensolver`. When a Householder step is degenerate — the subcolumn
below the subdiagonal is already zero, so `length == 0`, `mVVector` stays the
zero vector and the reflection actually applied to the matrix is the
identity — upstream still stores the reflection parameter `2/Dot(v,v) == 2`,
and `GetEigenvector` rebuilds `H = I - 2*e*e^T` instead of the identity. The
eigenvalues are unaffected (the stored value lands in the lower triangle,
which the diagonal copy never reads); `Q` is not. The port stores 0 instead.

This group is where the oracle first exercised that fix against the real
build: `ApprQuadratic3.fit` disagreed on one of its first 20 records, on the
single-point input `(0,1,3)`, whose `M = V*V^T` is rank 1 and decouples at the
first Householder step. Both sides return a legitimate eigenvector of the
9-fold eigenvalue 0 there — they are negatives of each other on the nonzero
components — so the symptom is a different (not a wrong) eigenvector; the
already-recorded reproduction in `docs/UPSTREAM-FINDINGS.md`, a block-diagonal
`A` with `max|Q*D*Q^T - A| = 4.0`, is the wrong-answer one.

**Probe.** `oracle/cpp/cases/v03-approximation.cpp` carries a verbatim copy of
upstream's `Tridiagonalize` (`EigenDecouples`) and verbatim copies of the four
matrix assemblies (`BuildQuadratic2M`, `BuildQuadraticCircle2M`,
`BuildQuadratic3M`, `BuildQuadraticSphere3M`, including `ApprQuadratic3`'s
documented omission of the trailing `M(0,0) = 1`). `EigenDecouples` reports
whether any Householder step is degenerate, which is an **exact**
characterisation of the inputs on which the two implementations can differ:
the fix changes only the value stored for such a step. The six main cases
reject those inputs by rejection sampling (capped at 24 attempts, redrawing
the whole point set — count, mode and every coordinate — on each attempt,
which is everything the predicate depends on, and falling back to the last
candidate); the four `.decoupledDeviation` cases keep only them, with point
counts of 1 to 5 (2D) or 1 to 6 (3D) so that `M` is rank-deficient.
`ApprQuadraticCircle2` deviates least often (79 of 2000) because its matrix is
only 4x4, so there are just two Householder steps to be degenerate.

### `ApprParallelLines2`

The port fixes four defects of upstream's `Fit` (issue #91): `ComputeF` writes
`a30[1] = -3` where the equivalent evaluation in `UpdateParameters` uses
`-3*Z12`; `Fit` accepts roots with `sigma^2 > 1`, which cannot come from a
unit direction; the `f1 == 0` branch sets `gamma = sqrt(sigma)` instead of
`sqrt(1 - sigma^2)`; and `Fit` reads `Polynomial1::operator[]` past the end of
its coefficient vector whenever the arithmetic operators eliminate leading
zeros.

The fourth of these is undefined behaviour in C++ (MSVC reads heap memory in
release), so **both** cases reject the draws that reach it, by checking the
degrees of upstream's own `f0`, `f1` and `h`. The first one is a wrong
coefficient, not a reassociation, so upstream's `h(sigma^2)` is a different
polynomial from the corrected one on essentially every input and the two fits
agree only when neither polynomial's roots improve the seed
`(sigma, gamma) = (0, 1)` — about 4 % of the draws. The main case therefore
needs a 400-attempt cap to find an agreeing draw (it redraws the point set and
`maxIterations` on every attempt); it is a narrow but genuine window in which
the shared arithmetic — the moments, `RootsPolynomial::Find`, the
`UpdateParameters` error function and the final projection of `C` onto `U` —
is compared bit for bit. The deviation case takes the complement and deviates
on every record.

The case file carries a copy of upstream's `Fit` with the first three fixes
switchable (`parallel::FitFixed`), so the split is decided by comparing
upstream's own `(C, V, radius)` with the corrected one bit for bit rather than
by a proxy. The corrected copy was validated against the port on 20 records
down to the coefficients of `f0` and `f1`, the reduced polynomials, the root
list and the final result: all bit-identical.

Upstream defects recorded for this group's headers that the port *preserves*
need no deviation case; the oracle compares them bit for bit and they agree:

- `ApprCircle2.h` / `ApprSphere3.h`'s inverted `initialCenterIsAverage`
  documentation (#92) — behaviour only, both cases drive both settings.
- `ApprQuery.h`'s `RANSAC` re-seeding and its `true`-without-fitting path
  (#378) — see "Not covered" for the shuffling path.
- `ApprCurveByArcs.h`'s always-true `fabs(det) >= epsilon` test with the
  default `epsilon = 0` (#163): `compute` divides by a near-zero determinant
  on every arc and both sides produce the same large or infinite centres.
- `ApprPolynomialSpecial2/3/4.h`'s per-axis strictly-increasing assert, which
  rejects the affine model `{1, x, y}` the header documents as admissible, and
  its unguarded division by the sample range (#163, #380). The
  `constructor.affineAssert` case pins the first as a throw-parity case; the
  NaN models produced by the second are compared bit for bit.
- `ApprQuadratic3.h`'s missing trailing `M(0,0) = 1` and
  `ApprQuadraticCircle2` / `ApprQuadraticSphere3`'s unguarded division by the
  last eigenvector component (#163, #380): the cocircular and cospherical
  cases reach centres of `(-Infinity, NaN)` with a reported measure of 0, and
  the recorded bits agree.

## Not covered

All 20 headers of group 3 are implemented by the port and every computational
entry point has a case. Three deliberate restrictions are worth recording.

- **`ApprQuery::RANSAC`'s shuffling path.** The loop calls
  `std::shuffle(candidates.begin(), candidates.end(), std::default_random_engine())`.
  Both the engine (MSVC makes `default_random_engine` a `mt19937`) and
  `std::shuffle`'s consumption of it are implementation-defined, so the
  permutation is not a property MSVC and V8 can be expected to agree on, and
  the port's own generator (a `minstd`-style Lehmer engine with Fisher-Yates)
  deliberately differs — `src/ApprQuery.ts` says so. The two paths that do not
  shuffle are covered: `numObservations < GetMinimumRequired()`
  (`ransac.tooFew`, an immediate `false` with `bestConsensus` untouched) and
  `numObservations == GetMinimumRequired()` (`ransac.minimum`, which returns
  the identity consensus and `bestModel.Fit(observations)`). The upstream
  defects of the shuffling path (#378: the engine is re-seeded every
  iteration, so `numIterations` has no observable effect; and `RANSAC` can
  return `true` with `bestModel` never fitted) are preserved in the port and
  documented, but they cannot be demonstrated against the C++ build.
- **`ApprQuery::ValidIndices` with `GTE_APPR_QUERY_VALIDATE_INDICES`
  defined.** Upstream's index validation is a compile-time toggle that is off
  in the default build, and defining it would change every case in the
  translation unit. The port turns it into the runtime flag
  `ApprQuery.validateIndices`, whose default (`false`) is what the oracle
  compares. The validating branch has its own unit tests in
  `test/ApprQuery.test.ts`.
- **Empty index lists.** `FitIndexed` with zero indices divides by
  `numIndices == 0`. `ApprHeightLine2` and `ApprHeightPlane3` reach GTE's
  `Vector::operator/=`, which zeroes the vector for a zero divisor, while the
  other fitters compute `1/0 = Infinity` and multiply, producing NaN
  covariances that are then fed to an eigensolver. The generators draw at
  least one index, so this configuration is not sampled: the NaN path makes
  every real output NaN, which the harness cannot distinguish (any NaN matches
  any NaN), so the records would be vacuous. The port reproduces upstream's
  divisor semantics explicitly and `test/ApprHeightLine2.test.ts` pins them.

## Upstream bug suspects

None new for this group's headers. One observation worth recording, about a
header outside the group that these fitters reach:

- **`SymmetricEigensolver::ComputePermutation` sorts with `std::sort`,** which
  is not stable, so the order of tied eigenvalues — and therefore which
  eigenvector `GetEigenvector(0)` returns for a degenerate smallest
  eigenvalue — is unspecified. The port uses `Array.prototype.sort`, which is
  stable. For MSVC the two happen to agree on every input this group reaches,
  because `std::sort` falls back to insertion sort for ranges of at most 32
  elements and the matrices here are 4x4 to 10x10; the agreement is therefore
  an accident of the implementation, not a guarantee. `SingularValueDecomposition`
  has the same pattern and it is already recorded under item 5 of the
  `SymmetricEigensolver.h` entry in `docs/UPSTREAM-FINDINGS.md`; the
  `SymmetricEigensolver` instance of it is not, and should be added there.
  It is not a wrong answer — every tie order gives a valid eigenbasis — but it
  makes the returned eigenvector platform-dependent, which matters for a
  caller such as `ApprQuadratic2` that hands the eigenvector straight back as
  conic coefficients.
