# Verify group 4 (approximation) against the MSVC-built upstream GTE

Family `v04-approximation`, 36 cases in `oracle/cpp/cases/v04-approximation.cpp`
and `test/oracle/v04-approximation.oracle.test.ts`. Deep run:
`npm run oracle:deep -- 2000 v04-approximation`, 72000 records (8033 of them
upstream `LogAssert` throws, all matched by the port), 65 s wall time, all 37
tests green.

The group is ten headers of least-squares fitters and iterative minimizers.
Three of them (`ApprParabola2`, `ApprParaboloid3`, `ApprGreatCircle3`) and two
of `ApprCylinder3`'s three point-fitting modes use only `+ - * / sqrt`, so
they are declared `exact` and agree bit for bit on every record. The rest
reach the C math library, either directly (`atan2` in `ApprGreatArc3` and
`ApprCone3::ComputeInitialCone`, `pow` in `ApprEllipseByArcs`, `sin`/`cos` in
`ApprTorus3`, `ApprCylinder3`'s hemisphere search and
`ApprCone3EllipseAndPoints`) or through `RootsPolynomial::SolveCubic`
(`ApprEllipse2`, `ApprEllipsoid3`, `ApprTorus3`), and are compared with a
tolerance.

Three levers keep libm out of code where it would decide control flow rather
than only perturb a value:

* `ApprCone3`'s iterated cases feed the minimizer a caller-supplied cone whose
  angle is one of the 24 dyadic values `j/16`. Every `k/64` for `k = 1..100`
  was compared between the MSVC `std::cos` and V8's `Math.cos`; only `k = 25`
  and `k = 54` differ, and neither is a multiple of 4, so `cos(j/16)` is
  bit-identical on both platforms. With the cosine fixed, the whole
  Gauss-Newton and Levenberg-Marquardt iteration (Jacobian, `MultiplyATB`,
  Cholesky factor and the two triangular solves) is arithmetic, and every
  output but the final `acos`-derived angle is held to bit-identity with
  `io.outRealExact`. Both cases agree to the last bit on all 2000 records.
* `ApprGreatArc3` picks the largest gap between `atan2` angles. The generator
  places the samples either in one narrow cluster (the wrap-around gap wins)
  or in two clusters separated by a 2-radian hole (an interior gap wins), so
  the winning gap leads by more than a radian and no rounding can change it.
* `ApprCylinder3`'s hemisphere and mesh searches take an argmin over
  `cos`/`sin`-generated directions. An even number of theta samples puts two
  candidates on perpendicular coordinate planes whose projected measures are
  exactly equal for a symmetric lattice point set; the first version of the
  case drew 2 to 4 samples and one record in 2000 picked a different winner on
  each side. Restricting the count to 3 or 5 removes the tie and both cases
  agree within 6e-15.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `ApprParabola2.h` | `fit`, `fitRobust`, `fit.throw`, `fitRobust.throw` | exact | bit-identical, 8000 records |
| `ApprParaboloid3.h` | `fit`, `fitRobust`, `fit.throw`, `fitRobust.throw` | exact | bit-identical, 8000 records |
| `ApprGreatCircle3.h` | `ApprGreatCircle3.compute` | exact | bit-identical, 2000 records |
| `ApprGreatCircle3.h` (`ApprGreatArc3`) | `ApprGreatArc3.compute` | tolerance (`atan2`) | bit-identical, 2000 records |
| `ApprEllipseByArcs.h` | `approximate`, `approximate.deviation` | tolerance (`pow`) | bit-identical, 4000 records |
| `ApprEllipse2.h` | `compute.box`, `compute.ellipse` | tolerance (`SolveCubic`) | max 1.8e-15 |
| `ApprEllipsoid3.h` | `compute.box`, `compute.ellipsoid`, `compute.degenerateBox` | tolerance (`SolveCubic`) | max 1.4e-15 |
| `ApprCylinder3.h` | `compute.eigenIndex`, `compute.specifiedAxis`, `compute.throw` | exact | bit-identical, 6000 records |
| `ApprCylinder3.h` | `compute.hemisphere`, `computeMesh` | tolerance (`sin`, `cos`) | max 5.6e-15 |
| `ApprCone3.h` | `gaussNewton.initialGuess`, `levenbergMarquardt.initialGuess` | tolerance, every field but the angle `outRealExact` | max 2.2e-16 (the `acos` angle only) |
| `ApprCone3.h` | `gaussNewton.computeInitialCone`, `levenbergMarquardt.computeInitialCone` | tolerance (`atan2`, `cos`, `acos`) | max 5.9e-14 |
| `ApprCone3.h` | `levenbergMarquardt.staleResidual.deviation` | deviation | 1028 of 2000 records deviate |
| `ApprTorus3.h` | `compute` | tolerance (`SolveCubic`) | max 1.4e-14 |
| `ApprTorus3.h` | `gaussNewton.computeInitialTorus`, `levenbergMarquardt.computeInitialTorus` | tolerance (`sin`, `cos`, `atan2`, `acos`) | max 2.4e-14 |
| `ApprTorus3.h` | `gaussNewton.initialGuess` | tolerance `1e-7` | max 1.25e-8 |
| `ApprTorus3.h` | `levenbergMarquardt.initialGuess` | tolerance `1e-10` | max 9.1e-12 |
| `ApprCone3EllipseAndPoints.h` | `fit`, `fit.throw` | tolerance (`sin`, `cos` through `Minimize1`) | max 7.0e-14 |
| `ApprCone3EllipseAndPoints.h` | `fit.deviation` | deviation | 1367 of 2000 records deviate |
| `ApprCone3EllipseAndPoints.h` (`ApprCone3ExtractEllipses`) | `extract` | tolerance, structural outputs only | bit-identical, 2000 records |

Generators alternate a small integer lattice, a uniform cloud and samples
placed on the target surface (a parabola or paraboloid with dyadic
coefficients at integer abscissas, points on a circle, cylinder, cone, torus
or ellipse), plus the degenerate sets that reach the failure branches: too few
points for the `LogAssert` (four throw-parity cases, 1100 to 1400 throwing
records each), all abscissas equal so the moment matrix is singular, one to
three points so the initial oriented box has a zero extent, a zero cylinder
axis, `numArcs < 2` and `a == b` for `ApprEllipseByArcs`, an invalid `Control`
and an empty point set for `ApprCone3EllipseAndPoints`.

## Port defects fixed

**1. `src/ApprParabola2.ts`, `src/ApprParaboloid3.ts`, `src/ApprCone3.ts`: the
sample average was computed with a division instead of a multiplication by the
reciprocal.** Upstream writes `average /= tNumPoints` (and `center /=
tNumPoints` in `ApprCone3::ComputeInitialCone`), and `Vector.h`'s
`operator/=(Vector&, Real)` multiplies each component by `1/scalar`. The port
divided each component. `a/b` and `a*(1/b)` differ in the last bit for most
`b`: for a coordinate sum of `-3` over 10 samples, `-3/10` is
`bfd3333333333333` and `-3*(1/10)` is `bfd3333333333334`. The error is 1 ulp
in the average, which `ApprParabola2::FitRobust` and
`ApprParaboloid3::FitRobust` carry into every centred coordinate (4 of 20 and
6 of 20 golden records disagreed, with scaled errors up to 1.2e-14 in the
coefficients), and which `ApprCone3::ComputeInitialCone` amplifies to 6.1e-10
in the reported `minError` of the initial cone. Regression tests:
`test/ApprParabola2.test.ts`, `test/ApprParaboloid3.test.ts` and
`test/ApprCone3.test.ts` each run the port on a point set whose coordinate sum
makes the two groupings differ and require the upstream grouping; the
`ApprCone3` test pins both the value the multiply grouping produces and the
one the divide grouping produced, so it fails if the fix is reverted.

## Deliberate deviations demonstrated

**`ApprEllipseByArcs.approximate.deviation` (issue #322).** The
intermediate-arc loop ignores `Circumscribe`'s return value, so a degenerate
point triple leaves `circle` holding the previous arc's centre and radius and
the function still reports success; the port propagates the failure. The
defect is reached when two consecutive interior points coincide, which happens
when `a` and `b` are adjacent doubles: the curvature range is then about
`3 * 2^-52` relative and `pow(ab/curv, 2/3)` rounds to the same double for
consecutive `i`. 628 of 2000 records deviate. Ordinary semi-axes never reach
it: 50000 uniform and 50000 extreme-aspect-ratio draws through the same
expression chain produced no degenerate triple, so the main case keeps the
full generator.

**`ApprCone3.levenbergMarquardt.staleResidual.deviation` (issue #261).**
`LevenbergMarquardtMinimizer::DoIteration` builds `-J^T*F` from the member
`mF`, which holds `F` at the previously *rejected* candidate whenever the
inner lambda-adjustment loop runs `DoIteration` more than once for the same
`pCurrent`; the port re-evaluates `F(pCurrent)`. `DoIteration` is repeated
exactly when the inner loop increments `numAdjustments`, so both the main case
and the deviation case choose `maxIterations` by running upstream itself with
`maxIterations = 1, 2, ...` and reading `result.numAdjustments`: the main case
keeps the longest prefix on which it stays 0 (an aimed construction, not a
rejection loop), the deviation case the first prefix on which it is not. 1028
of 2000 records deviate; the main case agrees to the last bit on all 2000.
`ApprTorus3.levenbergMarquardt.initialGuess` uses the same construction. The
first version of the Levenberg-Marquardt case, without it, disagreed on 9 of
20 golden records, including the iteration counts.

**`ApprCone3EllipseAndPoints.fit.deviation` (issue #349).** `ComputeCone`
computes `b / a` without validating the fitted ellipse extents and `Fit`
divides the accumulated error by `points.size()` with no empty-set guard; the
port asserts on both. Only the extent guard is observable: with an empty point
set the error function is NaN for every theta, upstream's own "Failed to find
fitted cone" assertion fires and the two sides agree on throwing (the
`fit.throw` case covers that, 1281 throwing records of 2000). A zero or
negative extent instead lets upstream return a degenerate cone, which the port
refuses: 1367 of 2000 records deviate.

## Not covered

* **`ApprEllipse2.h` and `ApprEllipsoid3.h`, the discarded `GetContainer`
  failure flag (issues #224 and #322).** The port raises an error where
  upstream reads the oriented box that `GetContainer` left untouched. No
  `deviation` case exists, because in the default upstream configuration the
  defect is unreachable: `GTE_APPR_QUERY_VALIDATE_INDICES` is not defined, so
  `ApprQuery::ValidIndices` returns `true` unconditionally,
  `ApprGaussian2/3::FitIndexed` never reports failure and `GetContainer` never
  returns `false`. A standalone probe compiled against the upstream headers
  confirms it: `ApprGaussian2<double>().Fit(0, nullptr)` returns `true`. The
  one input that would exercise the discarded flag, the empty point set,
  instead makes `GetContainer` dereference `points[0]` and terminates the
  process with an access violation, so it cannot be recorded. See the upstream
  suspects below.
* **`ApprCylinder3.h`'s multithreaded paths** (`ComputeMultiThreaded`,
  `FitToMeshMultiThreaded`): excluded by the task; the port runs them
  single-threaded anyway, and the single-threaded partition differs from the
  multithreaded one (upstream re-evaluates the north pole in every thread),
  which is already recorded in `docs/UPSTREAM-FINDINGS.md` under issue #224.
* **The centre, axes and extents of the ellipses that
  `ApprCone3ExtractEllipses::Extract` returns.** The case compares the plane
  count, every plane normal and constant, the index lists and the ellipse
  normals, all of which come from the oriented-box tree, `ApprGaussian3` and
  the iterative `SymmetricEigensolver3x3` and agree bit for bit on all 2000
  records. The remaining fields come from `ApprEllipse2` run with upstream's
  hard-coded 1024 iterations of a two-step gradient descent, each step solving
  a cubic with `pow`, `atan2`, `sin` and `cos`. On six to eight samples that
  iterate does not converge (it returns extents an order of magnitude larger
  than the generating ones), so the two sides wander apart: with those fields
  emitted, 247 of 2000 records differed, by up to 0.4 relative. No tolerance
  makes those numbers a test. `ApprEllipse2` itself is compared directly by
  its own cases, at 0 and 1.8e-15 relative error for 0 to 2 iterations.

## Independent-reference check

`ORACLE.md` asks for one check that the agreeing outputs are also right. Three
were made on the 2000-record deep run, by recomputing an independent quantity
from the recorded inputs and the recorded outputs. All three cases are
declared `exact` and agree bit for bit on every record, so the numbers below
are the port's as well as upstream's.

* `ApprParabola2.fit`, the 455 records whose samples lie exactly on a parabola
  `y = u0 x^2 + u1 x + u2` (dyadic coefficients, integer abscissas) and have
  at least three distinct abscissas: the fitted polynomial reproduces every
  sample to within 4.5e-11 absolute.
* `ApprParaboloid3.fit`, the 465 analogous records with at least six distinct
  `(x,y)`: 1.5e-12 absolute.
* `ApprCylinder3.compute.specifiedAxis`, the 666 records whose samples lie on
  a cylinder and whose specified axis is the generating axis: the distance of
  every sample from the *fitted* axis matches the fitted radius to within 3.9
  per cent, consistent with the 3 per cent radial perturbation the generator
  applies.

The parabola check first reported a residual of 143 on record 1062. The input
there is `{(2, 10.25), (3, 20.75), (3, 20.75)}`: only two distinct abscissas,
so the 3x3 moment matrix is rank deficient, its determinant merely rounds to
something nonzero, `Inverse(A, &invertible)` reports success and both
implementations return the same garbage coefficients. That is the behaviour of
a determinant-only rank test, not a defect of either side; the check now skips
those records. The generators were also checked against the library's
conventions: every frame is built right-handed (`RawFrame3` uses
`ComputeOrthogonalComplement`, whose `{N, U, V}` satisfies `N = U x V`, and
the 2D frames use `axis[1] = (-s, c) = -Perp(axis[0])`, since GTE's `Perp` is
the clockwise rotation), all directions are unit length and all extents are
positive except where a case is about violating that.

## Upstream bug suspects

**1. `ContOrientedBox2.h` and `ContOrientedBox3.h`: `GetContainer`
dereferences `points[0]` for an empty point set.** `ApprGaussian2/3::Fit`
returns `true` unconditionally in the default configuration
(`ApprQuery::ValidIndices` is compiled out unless
`GTE_APPR_QUERY_VALIDATE_INDICES` is defined), so `GetContainer` enters the
branch that reads `points[0]` to start the min/max scan and the process dies
with an access violation. The documented failure mode of
`ApprEllipse2`/`ApprEllipsoid3` (issue #224, #322, "reads the
default-constructed box after ignoring `GetContainer`'s failure flag") is
therefore not what happens in a default build: the flag is never false, and
the empty input crashes inside `GetContainer` instead. Reproduction: a
translation unit that calls `GetContainer(0, std::vector<Vector2<double>>{}
.data(), box)`.

**2. `ApprCone3EllipseAndPoints.h`: two paths of
`ApprCone3ExtractEllipses::Extract` terminate the process.** Both were hit
while building this group's generator, which is why the case selects its
`boxExtentEpsilon` with a replica of upstream's own plane location and
association (`PlanePointCounts` in the case file) instead of probing `Extract`
directly.
  * When `LocatePlanes` finds no flat box, `mPlanes` is empty,
    `AssociatePointsWithPlanes` leaves `minJ` at `size_t(-1)` and
    `mIndices[minJ]` is an out-of-bounds write. This is already recorded under
    issue #349; what is new is that it is an immediate access violation, not a
    silent corruption.
  * When a plane ends up with no supporting points, `ComputeEllipse` runs
    `ApprEllipse2` over an empty vector, which reaches suspect 1 above and
    also kills the process. Issue #349 describes this case as producing
    garbage ellipses; on a default build of the current upstream it crashes.

**3. `ApprEllipse2.h`: the 1024-iteration fit that
`ApprCone3ExtractEllipses::ComputeEllipse` hard-codes does not converge on
small point sets.** On six to eight samples taken exactly from an ellipse with
semi-axes in [2,4] and [0.5,1], `ApprEllipse2` returns extents of 16 to 30.
The two-step gradient descent is not a reliable fitter at that sample count,
and `ComputeEllipse` gives the caller no way to tune it. This is a quality
observation rather than a coding defect, and it is the reason the extract case
cannot compare those fields.

