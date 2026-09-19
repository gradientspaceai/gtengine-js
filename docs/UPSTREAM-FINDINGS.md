# Defects found in Geometric Tools Engine / Mathematics while porting it to TypeScript

This report collects every defect that a complete, test-driven port of the
GTE `Mathematics` headers found in the upstream C++ sources. It is written for
the maintainer of <https://github.com/davideberly/GeometricTools> and assumes no
knowledge of the port.

## What produced these findings

`gtengine-js` is a line-by-line TypeScript port of the header-only
`GTE/Mathematics` library. The baseline is upstream commit `d29e7758ae26`
(2026-08-13); every file, line number and reproduction below refers to that
commit unless the text says otherwise.

The findings came out of two passes:

1. **Porting.** Each header was translated function by function, with the
   upstream algorithm and its numerical behaviour preserved deliberately, and
   covered by unit tests. Anything that made the tests disagree with an
   independently derived expectation was investigated rather than papered over.
2. **Independent verification.** Every ported file was then re-read against the
   upstream header by a second reviewer who had not written the port, using a
   fixed translation-hazard checklist and property-based tests (random inputs
   checked against invariants, brute-force oracles, exact-arithmetic oracles,
   and mutation testing that reverts a fix to confirm a test fails without it).
   This pass confirmed, sharpened, corrected or refuted the porting-pass claims
   and found a further tranche of defects that unit tests had missed.

Because the port is in TypeScript, several classes of C++ undefined behaviour
show up as concrete, observable failures (an out-of-bounds read returns
`undefined` and poisons arithmetic instead of silently reading adjacent memory),
which is why a number of latent C++ issues below have precise reproductions.

## How to read this report

Every finding is tracked as an issue in the port's repository, labelled
`upstream-bug`. The issue holds the original analysis; its comments hold the
verification-pass confirmations, sharper measurements, corrections and
"now fixed in the port" status updates. Issue links are given for every row.

Severity is one of:

| Code | Meaning |
| --- | --- |
| **RC** | Result-corrupting: silently wrong numeric or geometric output, silent state corruption, a hang, or C++ undefined behaviour on reachable input. |
| **WR** | Wrong but recoverable: the failure is loud (an assert fires, an exception is thrown, a NaN propagates to an obvious place) or the wrong value provably never escapes the function. |
| **minor** | Dead code, unreachable branches, contract or API inconsistencies, performance, nondeterminism with mathematically identical results. |
| **doc** | Comment or documentation only; the code is correct. |

Port status is one of:

- **fixed**: the port deviates from upstream and implements the corrected
  behaviour, with a regression test that fails against the upstream formulation.
- **preserved**: the port reproduces the upstream behaviour deliberately, with a
  test that pins it, because the fix would change a public contract, require
  inventing a tolerance or a new algorithm, or because the upstream behaviour is
  arguably intentional.
- **corrected** / **dropped**: a documentation or dead-code item; the port
  carries the corrected comment, or omits the dead code, and the behaviour is
  unchanged.
- **fixed (assert)** / **fixed (throws)**: the port adds the guard that upstream
  is missing, turning undefined behaviour into a diagnosed error; behaviour on
  valid input is unchanged.
- **deviates**: the port deliberately departs from upstream on a point where
  upstream's behaviour is itself undefined or untranslatable, and says so in the
  source.
- **n/a**: the code path is not ported (for example the `GTE_USE_VEC_MAT`
  branches, or arbitrary-precision instantiations that were not needed).

Nothing here has been reported upstream before; this document is the report.

## Counts

- **493 distinct findings** across **157** tracked issues (one issue
  frequently holds several findings in related files).
- By severity: **242 result-corrupting**, **14 wrong but
  recoverable**, **164 minor**, **73 documentation**.
- By port status: **255 fixed or corrected in the port** (of which 154 are code
  fixes with regression tests, 22 are added guards or asserts where upstream has
  undefined behaviour, 64 are comment corrections, 11 are dead-code removals and
  4 are documented deliberate deviations), **229 preserved deliberately**, and
  **9 not ported** (the `GTE_USE_VEC_MAT` branches, dead code that cannot
  compile, and two arbitrary-precision paths).
- **288 distinct upstream headers** are implicated.

Nine claims made during the porting pass were later corrected, sharpened or
withdrawn by the verification pass; they are listed in
[Claims withdrawn or corrected](#claims-withdrawn-or-corrected) and the
corrected form is what appears above.

## Summary table

Issue links point at <https://github.com/gradientspaceai/gtengine-js/issues>. "Port" values: *fixed* = the port deviates and implements the corrected behaviour with a regression test; *preserved* = the port reproduces upstream deliberately, pinned by a test; *n/a* = the code path is not ported.

| Upstream file | Function / location | Symptom | Sev | Port | Issue |
| --- | --- | --- | --- | --- | --- |
| `AdaptiveSkeletonClimbing2.h` | `LinearMergeTree::GetEdge`, `GetRectangle` | `GetEdge` can never return -1, so the `!= -1` guards are dead code | minor | preserved | [#52](https://github.com/gradientspaceai/gtengine-js/issues/52) |
| `AdaptiveSkeletonClimbing2.h` | constructor comment | comment says `N >= 0` is accepted; the code rejects `N <= 0` | doc | preserved | [#52](https://github.com/gradientspaceai/gtengine-js/issues/52) |
| `AdaptiveSkeletonClimbing3.h` | root `Merge` call | the root call discards `Merge`'s return, so a whole-image monobox is dropped and the mesh comes back empty | RC | preserved | [#194](https://github.com/gradientspaceai/gtengine-js/issues/194) |
| `AdaptiveSkeletonClimbing3.h` | `Get{X,Y,Z}{Min,Max}EdgesM` | `GetZeroBase` may return -1 and index the image at -1 (latent) | minor | preserved | [#194](https://github.com/gradientspaceai/gtengine-js/issues/194) |
| `AdaptiveSkeletonClimbing3.h` | `GetVertices` | narrows `Real` box corners through `static_cast<float>` | minor | preserved | [#194](https://github.com/gradientspaceai/gtengine-js/issues/194) |
| `AlignedBox.h` | `GetVertices` comment | typos "whern", `vertex[i][d = max[d]` | doc | corrected | [#78](https://github.com/gradientspaceai/gtengine-js/issues/78) |
| `AlignedBoxBV.h` | `GetSplittingAxis` | dead store `maxExtent = extents[2]` | minor | dropped | [#268](https://github.com/gradientspaceai/gtengine-js/issues/268) |
| `APConversion.h` | `EstimateAmB` | `tMinSqr`/`tMaxSqr` are not recomputed in the bisection rounding branch, so the Newton bound on exhaustion uses a stale square | RC | fixed | [#280](https://github.com/gradientspaceai/gtengine-js/issues/280) |
| `APConversion.h` | `EstimateAmB` | undocumented `aSqr >= bSqr` precondition; otherwise the bracket is inverted | doc | preserved | [#280](https://github.com/gradientspaceai/gtengine-js/issues/280) |
| `APConversion.h` | `EstimateSqrt` | `EstimateSqrt(0)` returns `[0, 5e-324]`, violating the documented strict lower bound; no negative-input guard | minor | preserved | [#280](https://github.com/gradientspaceai/gtengine-js/issues/280) |
| `APConversion.h` | header comment | `r^2 - 7*r^2 + 1 = 0` should read `r^2 - 7*r + 1 = 0` | doc | corrected | [#280](https://github.com/gradientspaceai/gtengine-js/issues/280) |
| `APInterval.h` | `SetSign(+-2)` sentinels | the infinity sentinels violate the `BSNumber` invariant and compare smaller than finite values of the same sign | RC | preserved | [#280](https://github.com/gradientspaceai/gtengine-js/issues/280) |
| `APInterval.h` | `operator/` | divisor `[0,0]` throws instead of returning the documented `Reals()` | WR | preserved | [#280](https://github.com/gradientspaceai/gtengine-js/issues/280) |
| `APInterval.h` | `QFN2` | dead alias | minor | n/a | [#280](https://github.com/gradientspaceai/gtengine-js/issues/280) |
| `ApprCircle2.h`, `ApprSphere3.h` | `FitUsingLengths` comment | the `initialCenterIsAverage` description is inverted relative to the code | doc | corrected | [#92](https://github.com/gradientspaceai/gtengine-js/issues/92) |
| `ApprCone3.h` | `ComputeInitialCone` | unguarded `rRange/hRange` and `rMax/tanAngle` divisions | RC | preserved | [#271](https://github.com/gradientspaceai/gtengine-js/issues/271) |
| `ApprCone3.h` | `ComputeInitialCone` | the `ApprHeightLine2::Fit` failure return is discarded, giving a zeroed initial cone | RC | preserved | [#271](https://github.com/gradientspaceai/gtengine-js/issues/271) |
| `ApprCone3EllipseAndPoints.h` | `Fit` | divides by `points.size()` with no empty-set guard | RC | fixed (assert) | [#349](https://github.com/gradientspaceai/gtengine-js/issues/349) |
| `ApprCone3EllipseAndPoints.h` | `ComputeCone` | computes `b / a` without validating the fitted ellipse extents | RC | fixed (assert) | [#349](https://github.com/gradientspaceai/gtengine-js/issues/349) |
| `ApprCone3EllipseAndPoints.h` | `Control` | duplicated default member initializers for three of six fields | minor | preserved | [#349](https://github.com/gradientspaceai/gtengine-js/issues/349) |
| `ApprCone3EllipseAndPoints.h` | `AssociatePointsWithPlanes` | indexes `mIndices[size_t(-1)]` when `mPlanes` is empty | RC | fixed | [#349](https://github.com/gradientspaceai/gtengine-js/issues/349) |
| `ApprCone3EllipseAndPoints.h` | `ComputeEllipse` | runs on spurious planes with an empty index set (a three-point tree node always has a flat box) | RC | fixed | [#349](https://github.com/gradientspaceai/gtengine-js/issues/349) |
| `ApprCone3EllipseAndPoints.h` | `LocatePlanes` | planes with only 1-2 supporting points give an ellipse with infinite extents | RC | preserved | [#349](https://github.com/gradientspaceai/gtengine-js/issues/349) |
| `ApprCone3EllipseAndPoints.h` | `ProcessPlane` | compares a plane constant (a distance) against `mCosAngleEpsilon` (a cosine) | minor | preserved | [#349](https://github.com/gradientspaceai/gtengine-js/issues/349) |
| `ApprCone3EllipseAndPoints.h` | error reporting | `sqrt(error)/n` instead of `sqrt(error/n)` | minor | preserved | [#349](https://github.com/gradientspaceai/gtengine-js/issues/349) |
| `ApprCurveByArcs.h` | `fabs(det) >= epsilon` | `epsilon` defaults to 0 so the test is always true; a colinear triple divides by zero instead of taking the documented sentinel path | RC | preserved | [#163](https://github.com/gradientspaceai/gtengine-js/issues/163) |
| `ApprCurveByArcs.h` | `numArcs == 1` | the bisection chord collapses on a closed curve, giving a NaN or meaningless arc | RC | preserved | [#163](https://github.com/gradientspaceai/gtengine-js/issues/163) |
| `ApprCurveByArcs.h` | comment | `{t0,t1} = {times[2*i], times[2*i+1]}` should be `times[2*i+2]` | doc | corrected | [#163](https://github.com/gradientspaceai/gtengine-js/issues/163) |
| `ApprCylinder3.h` | `FinishCylinder` | `hmax` initialised to 0 instead of `-max()` | minor | preserved | [#224](https://github.com/gradientspaceai/gtengine-js/issues/224) |
| `ApprCylinder3.h` | multithreaded paths | the partition starts at `j = 0` (north pole re-evaluated) and leaves threads idle when `numThreads > numPhiSamples` | minor | preserved | [#224](https://github.com/gradientspaceai/gtengine-js/issues/224) |
| `ApprCylinder3.h` | class and mesh docs | stale "cone axis/height/center" text; RMS claimed where the mean squared error is returned | doc | corrected | [#224](https://github.com/gradientspaceai/gtengine-js/issues/224) |
| `ApprEllipse2.h` | `operator()` | `GetContainer`'s bool is ignored; on failure `box` is read uninitialised | RC | fixed (throws) | [#224](https://github.com/gradientspaceai/gtengine-js/issues/224) |
| `ApprEllipse2.h` | header comment | the `H(t)` paragraph says "smallest G-value" and `M - T * dF/dC`; `UpdateMatrix` names a renamed variable | doc | corrected | [#224](https://github.com/gradientspaceai/gtengine-js/issues/224) |
| `ApprEllipse2.h` | `operator()` | dead store: the initial `ErrorFunction` result is overwritten whenever `numIterations > 0` | minor | preserved | [#224](https://github.com/gradientspaceai/gtengine-js/issues/224) |
| `ApprEllipseByArcs.h` | intermediate-arc loop | `Circumscribe`'s failure is ignored, so the previous arc's circle is stored and success is still reported | RC | fixed | [#322](https://github.com/gradientspaceai/gtengine-js/issues/322) |
| `ApprEllipseByArcs.h` | header comment | calls the outputs "returned input arrays"; writes `numArc` for `numArcs` | doc | corrected | [#322](https://github.com/gradientspaceai/gtengine-js/issues/322) |
| `ApprEllipsoid3.h` | initial oriented box | `GetContainer`'s return is ignored and the default-constructed box is read | RC | fixed (assert) | [#322](https://github.com/gradientspaceai/gtengine-js/issues/322) |
| `ApprEllipsoid3.h` | header comment | "smallest G-value", `dF/dC`, typos "definitess" and "ellisoid" | doc | corrected | [#322](https://github.com/gradientspaceai/gtengine-js/issues/322) |
| `ApprGaussian3.h` | `GetMinimumRequired` | returns 2, so RANSAC seeds a 3D Gaussian from a rank-1 covariance | minor | preserved | [#378](https://github.com/gradientspaceai/gtengine-js/issues/378) |
| `ApprOrthogonalLine3.h` | header comment | says "the minimum eigenvalue is unique"; the code correctly tests the maximum | doc | corrected | [#380](https://github.com/gradientspaceai/gtengine-js/issues/380) |
| `ApprParabola2.h`, `ApprParaboloid3.h` | all `Fit` and `FitRobust` | the reported "mean square error" is `sqrt(sum)/n`, neither MSE nor RMS | minor | preserved | [#271](https://github.com/gradientspaceai/gtengine-js/issues/271) |
| `ApprParaboloid3.h` | comment | "the i-index into A(2,2) is 35" should reference `A(5,5)` | doc | corrected | [#271](https://github.com/gradientspaceai/gtengine-js/issues/271) |
| `ApprParallelLines2.h` | `ComputeF` | `a30[1] = -mR3;` is missing the factor `Z12`; the fit differs on every random input | RC | fixed | [#91](https://github.com/gradientspaceai/gtengine-js/issues/91) |
| `ApprParallelLines2.h` | `Fit` | roots with `sigma^2 > 1` are not rejected; non-unit direction and NaN radius | RC | fixed | [#91](https://github.com/gradientspaceai/gtengine-js/issues/91) |
| `ApprParallelLines2.h` | `Fit` | `gamma = sqrt(sigma)` in the `f1 == 0` branch violates `gamma^2 + sigma^2 = 1` | RC | fixed | [#91](https://github.com/gradientspaceai/gtengine-js/issues/91) |
| `ApprParallelLines2.h` | `Fit` | out-of-range `Polynomial1::operator[]` reads after leading-zero elimination | RC | fixed | [#91](https://github.com/gradientspaceai/gtengine-js/issues/91) |
| `ApprParallelLines2.h` | `Fit` | the returned direction is never renormalised (`Dot(V,V)` observed at 1.0019 and 1.8e169) | RC | preserved | [#380](https://github.com/gradientspaceai/gtengine-js/issues/380) |
| `ApprPolynomialSpecial3.h`, `ApprPolynomialSpecial4.h` | constructor | asserts each degree list is separately increasing, rejecting the documented distinct-pair term sets | minor | preserved | [#163](https://github.com/gradientspaceai/gtengine-js/issues/163) |
| `ApprPolynomialSpecial2/3/4.h` | `Transform` | divides by the sample range unguarded; constant data yields a NaN model reported as success | RC | preserved | [#380](https://github.com/gradientspaceai/gtengine-js/issues/380) |
| `ApprQuadratic2.h`, `ApprQuadratic3.h` | comment | says `M = (sum V)(sum V)^T` (rank 1); the code accumulates `sum V V^T` | doc | corrected | [#163](https://github.com/gradientspaceai/gtengine-js/issues/163) |
| `ApprQuadratic3.h` | `Fit` | missing the trailing `M(0,0) = 1` its three siblings have | minor | preserved | [#163](https://github.com/gradientspaceai/gtengine-js/issues/163) |
| `ApprQuadraticCircle2.h`, `ApprQuadraticSphere3.h` | `Fit` | divide by the last eigenvector component unguarded; NaN center and radius reported with measure 0 | RC | preserved | [#380](https://github.com/gradientspaceai/gtengine-js/issues/380) |
| `ApprQuery.h` | `RANSAC` | the random engine is re-seeded inside the iteration loop, so `numIterations` has no observable effect | RC | preserved | [#378](https://github.com/gradientspaceai/gtengine-js/issues/378) |
| `ApprQuery.h` | `RANSAC` | can return `true` without ever fitting `bestModel` | RC | preserved | [#378](https://github.com/gradientspaceai/gtengine-js/issues/378) |
| `ApprTorus3.h` | `operator()` | unguarded `1/b0` and unguarded `SolveCubic` leading coefficient | RC | preserved | [#271](https://github.com/gradientspaceai/gtengine-js/issues/271) |
| `ApprTorus3.h` | header comment | defines `a2` twice (the third occurrence should be `a0`) | doc | corrected | [#271](https://github.com/gradientspaceai/gtengine-js/issues/271) |
| `ApprTorus3.h` | `Rational f0..f3` | dead code: the intended exact root classification never happens | minor | preserved | [#271](https://github.com/gradientspaceai/gtengine-js/issues/271) |
| `Arc2.h` | `Contains(P, epsilon)` | a negative epsilon cannot behave "as if zero"; the function returns false for every point | minor | preserved | [#155](https://github.com/gradientspaceai/gtengine-js/issues/155) |
| `Array4.h` | `SetPointers` comment | "else 'other' is an empty Array3." copy-pasted from `Array3.h` | doc | corrected | [#363](https://github.com/gradientspaceai/gtengine-js/issues/363) |
| `ASinEstimate.h` | `C_ASIN_EST_MAX_ERROR` | disagrees with `C_ACOS_EST_MAX_ERROR` in trailing digits although provably the same quantity | minor | preserved | [#57](https://github.com/gradientspaceai/gtengine-js/issues/57) |
| `BasisFunction.h` | `GetIndex` | the `t <= tmin` and `t >= tmax` shortcuts ignore knot multiplicity; all-zero basis at a degenerate endpoint | RC | preserved | [#415](https://github.com/gradientspaceai/gtengine-js/issues/415) |
| `BasisFunction.h` | `GetValue` | doc promises zeros outside `[minIndex, maxIndex]`; `Evaluate` leaves stale nonzero values | doc | preserved | [#415](https://github.com/gradientspaceai/gtengine-js/issues/415) |
| `BasisFunction.h` | `Create` | validates interior multiplicities against `d+1` where the class doc requires `d`; no total-multiplicity check | minor | preserved | [#415](https://github.com/gradientspaceai/gtengine-js/issues/415) |
| `BitHacks.h` | header comment | `GetTrailingBit(10) = 2` is wrong (it is 1); a hex prefix on binary digits | doc | corrected | [#67](https://github.com/gradientspaceai/gtengine-js/issues/67), [#363](https://github.com/gradientspaceai/gtengine-js/issues/363) |
| `BSNumber.h` | `ConvertToInteger` | the format asserts sit inside `if (number.size() > 1)`, so one-character garbage is accepted | RC | fixed | [#95](https://github.com/gradientspaceai/gtengine-js/issues/95) |
| `BSPPolygon2.h` | `operator&`, `operator-` | `Finalize()` asserts a nonempty edge list, so an empty intersection throws | WR | preserved | [#169](https://github.com/gradientspaceai/gtengine-js/issues/169) |
| `BSPPolygon2.h` | `SplitEdge` | `std::map::insert` is a silent no-op on an existing key, leaving an edge index unmapped | RC | preserved | [#169](https://github.com/gradientspaceai/gtengine-js/issues/169) |
| `BSPPolygon2.h`, `BSPTree2.h` | class comments | the two point-classification sign conventions contradict; `BSPTree2`'s is inverted | doc | corrected | [#388](https://github.com/gradientspaceai/gtengine-js/issues/388) |
| `BSPrecision.h` | `operator+` | pairs one operand's `maxExponent` with the other's `minExponent`; under-estimates `maxBits` for mixed sets and is asymmetric | RC | preserved | [#366](https://github.com/gradientspaceai/gtengine-js/issues/366) |
| `BSplineReduction.h` | Gram matrix | `Integration::Romberg` is applied across knots, giving ~1e-3 entry error | RC | preserved | [#169](https://github.com/gradientspaceai/gtengine-js/issues/169) |
| `BSplineSurface.h` | header comment (L10-18) | stale copy-pasted comment describing `BSplineReduction` | doc | dropped | [#96](https://github.com/gradientspaceai/gtengine-js/issues/96) |
| `BSRational.h` | two-`BSNumber` constructor | the unconditional exponent adjustment produces an invalid zero encoding | RC | fixed | [#168](https://github.com/gradientspaceai/gtengine-js/issues/168) |
| `BSRational.h` | string constructor | unconditional `SetSign(sign)` yields an invalid negative zero for `-0.0` | RC | fixed | [#168](https://github.com/gradientspaceai/gtengine-js/issues/168) |
| `BSRational.h` | string constructor | the `x.` branch is unreachable | minor | dropped | [#168](https://github.com/gradientspaceai/gtengine-js/issues/168) |
| `BVTree.h` | `GetLeafIndices` | a leaf is never tested against its own bounding volume, only its parent's | minor | preserved | [#103](https://github.com/gradientspaceai/gtengine-js/issues/103) |
| `BVTreeOfTriangles.h` | `Execute`, `Intersection::operator<` | the `std::set` orders by `parameter` only, so coincident hits are silently dropped | RC | fixed | [#167](https://github.com/gradientspaceai/gtengine-js/issues/167) |
| `BVTreeOfTriangles.h` | `IntersectSegmentTriangle` | reports the centered-form parameter where `BVTree.h` documents `t` in `[0,1]` | minor | preserved | [#387](https://github.com/gradientspaceai/gtengine-js/issues/387) |
| `CholeskyDecomposition.h` | run-time `BlockCholeskyDecomposition` | the block-level `GetIndex` stride `NumBlocks` is used for the scalar offset inside blocks | RC | fixed | [#209](https://github.com/gradientspaceai/gtengine-js/issues/209) |
| `CholeskyDecomposition.h` | run-time class | comments promise an `N > 0` check that does not exist | minor | fixed (assert) | [#478](https://github.com/gradientspaceai/gtengine-js/issues/478) |
| `ChebyshevRatioEstimate.h` | R-variant comment | misstates the x-domain | doc | preserved | [#57](https://github.com/gradientspaceai/gtengine-js/issues/57) |
| `CircleThroughPointSpecifiedTangentAndRadius.h` | header comment | `Perp(N) = (-n1,n0)` is the wrong convention; wrong base point; bisector-origin copy-paste | doc | corrected | [#101](https://github.com/gradientspaceai/gtengine-js/issues/101) |
| `CLODPolyline.h` | `ComputeEdges` | the open-polyline branch reads `permute[numVertices]`, one past the end | RC | fixed | [#182](https://github.com/gradientspaceai/gtengine-js/issues/182) |
| `Cone.h` | `CreateMesh` | `tNumExtra` is an undocumented, non-scale-invariant replacement for the original density rule | minor | preserved | [#268](https://github.com/gradientspaceai/gtengine-js/issues/268) |
| `Cone.h` | `CreateMesh` | no lower bound on `numMinVertices` | minor | fixed (assert) | [#268](https://github.com/gradientspaceai/gtengine-js/issues/268) |
| `ConformalMapGenus0.h` | Laplacian build (~L89) | `element.second->T[1]` is dereferenced unconditionally; null for a boundary edge | RC | fixed (assert) | [#318](https://github.com/gradientspaceai/gtengine-js/issues/318) |
| `ConformalMapGenus0.h` | bisection (~L375) | `LogAssert(fmin > 0 && fmax < 0)` fires on reachable valid input | WR | preserved | [#318](https://github.com/gradientspaceai/gtengine-js/issues/318) |
| `ConformalMapGenus0.h` | CG solve (~L155/175) | convergence flag off by one | minor | preserved | [#318](https://github.com/gradientspaceai/gtengine-js/issues/318) |
| `ConformalMapGenus0.h` | bisection (~L386) | dead store to `fmax` | minor | dropped | [#318](https://github.com/gradientspaceai/gtengine-js/issues/318) |
| `ConstrainedDelaunay2.h` | `operator()` | `mInsertedEdges` is never cleared, so a reused object inherits stale edge keys | RC | fixed | [#325](https://github.com/gradientspaceai/gtengine-js/issues/325) |
| `ConstrainedDelaunay2.h` | `Insert` | reads `duplicates[edge[i]]` before range-checking `edge[i]` | RC | fixed | [#325](https://github.com/gradientspaceai/gtengine-js/issues/325) |
| `ConstrainedDelaunay2.h` | `Retriangulate` | the strip fill is not constrained-Delaunay; nothing says so | doc | preserved | [#325](https://github.com/gradientspaceai/gtengine-js/issues/325) |
| `ContAlignedBox.h` | `MergeContainers` | the `bool` return is vestigial (only `true` is possible) | minor | preserved | [#106](https://github.com/gradientspaceai/gtengine-js/issues/106) |
| `ContAlignedBox2Arc2.h` | `GetContainer` | `0 < numPoints && numPoints <= 6` is always true; trailing `return false` unreachable | minor | dropped | [#174](https://github.com/gradientspaceai/gtengine-js/issues/174) |
| `ContCircle2.h`, `ContSphere3.h` | `GetContainer` | reads `points[0]` and divides by `numPoints` before any count check | RC | fixed (throws) | [#106](https://github.com/gradientspaceai/gtengine-js/issues/106) |
| `ContCone.h` | includes | relies on transitive includes for `Vector`/`Dot` | minor | preserved | [#292](https://github.com/gradientspaceai/gtengine-js/issues/292) |
| `ContEllipse2.h`, `ContEllipsoid3.h` | `GetContainer` | only strictly negative eigenvalues are adjusted; an exactly zero eigenvalue gives an infinite extent | RC | preserved | [#292](https://github.com/gradientspaceai/gtengine-js/issues/292) |
| `ContEllipse2.h`, `ContEllipsoid3.h` | `MergeContainers` | the merged ellipse does not contain its inputs | RC | preserved | [#292](https://github.com/gradientspaceai/gtengine-js/issues/292) |
| `ContEllipse2.h`, `ContEllipsoid3.h` | `InContainer` | rejects the extremal input point that defines the container by one ulp | minor | preserved | [#409](https://github.com/gradientspaceai/gtengine-js/issues/409) |
| `ContEllipse2MinCR.h` | `MaxProduct` | divides by zero whenever the walk steps onto a vertical constraint line | RC | fixed | [#234](https://github.com/gradientspaceai/gtengine-js/issues/234) |
| `ContEllipsoid3MinCR.h` | `FindFacetMax`, `FindEdgeMax` | `LogAssert(numer >= 0)` contradicts the adjacent comment; without it the walk leaves the feasible region | RC | fixed | [#234](https://github.com/gradientspaceai/gtengine-js/issues/234) |
| `ContEllipsoid3MinCR.h` | `FindFacetMax`, `FindEdgeMax` | the mutual recursion does not terminate for ~0.3% of random integer clouds | RC | fixed | [#409](https://github.com/gradientspaceai/gtengine-js/issues/409) |
| `ContEllipsoid3MinCR.h` | `FindFacetMax` | degenerate branch can recurse with `plane0 = plane1 = -1`, indexing `A[-1]` | minor | preserved | [#409](https://github.com/gradientspaceai/gtengine-js/issues/409) |
| `ContEllipsoid3MinCR.h` | facet/edge walk | stalls at a polytope vertex on ~4% of clouds, giving a containing but non-minimal ellipsoid | minor | preserved | [#234](https://github.com/gradientspaceai/gtengine-js/issues/234) |
| `ContLozenge3.h` | `GetContainer` | the rectangle is centered on a corner of the parameter interval, so the lozenge does not contain the points | RC | fixed | [#174](https://github.com/gradientspaceai/gtengine-js/issues/174) |
| `ContOrientedBox3.h` | `GetContainer` | reads `points[0]` with no count guard | RC | fixed (assert) | [#292](https://github.com/gradientspaceai/gtengine-js/issues/292) |
| `ContOrientedBox3.h`, `ContEllipsoid3.h` | `MergeContainers` | `Rotation<3,Real>(rot)` assumes `det = +1`; eigenvector frames are not guaranteed right-handed | RC | preserved | [#292](https://github.com/gradientspaceai/gtengine-js/issues/292) |
| `ContOrientedBox3.h` | `MergeContainers` | seeds `pmin`/`pmax` at zero rather than the first projected vertex | minor | preserved | [#292](https://github.com/gradientspaceai/gtengine-js/issues/292) |
| `ContPointInPolygon2.h` | constructor comment | claims the class stores a copy of the points; the member is a bare pointer | doc | corrected | [#106](https://github.com/gradientspaceai/gtengine-js/issues/106) |
| `ContPointInPolyhedron3.h` | `Face::indices` | declared `std::array<int32_t,3>` while `indices.size()` is read as the true vertex count; faces are truncated | RC | fixed | [#343](https://github.com/gradientspaceai/gtengine-js/issues/343) |
| `ContPointInPolyhedron3.h` | `mNumPoints` | dead member | minor | dropped | [#343](https://github.com/gradientspaceai/gtengine-js/issues/343) |
| `ContPointInPolyhedron3.h` | `Contains` | returns a silent `false` for unsupported (type, method) pairs | minor | preserved | [#343](https://github.com/gradientspaceai/gtengine-js/issues/343) |
| `ContScribeCircle2.h` | `Inscribe` | writes a degenerate circle into the output before returning `false` | minor | fixed | [#292](https://github.com/gradientspaceai/gtengine-js/issues/292) |
| `ConvertCoordinates.h` | affine example comment | the tuples documented as `B.GetCol(0..3)` are the rows of `B` | doc | corrected | [#160](https://github.com/gradientspaceai/gtengine-js/issues/160) |
| `ConvertCoordinates.h` | 3D example comment | a 4-tuple used as a 3D basis vector; the snippet declares `cs`/`sn` then uses `c`/`s` | doc | corrected | [#160](https://github.com/gradientspaceai/gtengine-js/issues/160) |
| `ConvexHull2.h` | `GetTangent` | silently returns whatever indices it last held if the bounding loop expires | RC | preserved | [#277](https://github.com/gradientspaceai/gtengine-js/issues/277) |
| `ConvexHull3.h` | `GetHull()` comment | claims `E = T/2` satisfies Euler's formula; for a closed triangle mesh `E = 3T/2` | doc | corrected | [#325](https://github.com/gradientspaceai/gtengine-js/issues/325) |
| `ConvexHull3.h` | `SelectSplit` | binds scratch-pool slots by const reference and then mutates them | minor | n/a | [#325](https://github.com/gradientspaceai/gtengine-js/issues/325) |
| `ConvexPolyhedron3.h` | constructor | validates `indices.size() >= 12` but never `% 3 == 0`, unlike `Polyhedron3.h` | minor | preserved | [#175](https://github.com/gradientspaceai/gtengine-js/issues/175) |
| `CubicRootsQR.h`, `QuarticRootsQR.h` | `GetQuadraticRoots` | real roots of even multiplicity are silently dropped (rounded-negative discriminant) | RC | preserved | [#488](https://github.com/gradientspaceai/gtengine-js/issues/488) |
| `CurvatureFlow2.h` | line 48 | mixed-derivative coefficient is 0.5 where the curvature-flow numerator has 2 | RC | preserved | [#123](https://github.com/gradientspaceai/gtengine-js/issues/123) |
| `CurveExtractor.h` | `MakeUnique` | remapped edge indices are not re-sorted, so reversed duplicates survive dedup | minor | preserved | [#67](https://github.com/gradientspaceai/gtengine-js/issues/67), [#362](https://github.com/gradientspaceai/gtengine-js/issues/362) |
| `Cylinder3.h` | class comment | claims the default constructor sets the axis to (0,0,1); `Line3<T>()` gives (1,0,0) | doc | corrected | [#155](https://github.com/gradientspaceai/gtengine-js/issues/155) |
| `Delaunay2.h` | `GetNumVertices` | returns 0 after degenerate input while `GetVertices` still returns the caller's pointer | RC | fixed | [#277](https://github.com/gradientspaceai/gtengine-js/issues/277) |
| `Delaunay2.h` | `GetContainingTriangle` comment | "inside all four edges" is a 3D copy-paste | doc | corrected | [#277](https://github.com/gradientspaceai/gtengine-js/issues/277) |
| `Delaunay2.h`, `Delaunay3.h` | non-void functions | fall off the end after a bare `LogError` | minor | preserved | [#277](https://github.com/gradientspaceai/gtengine-js/issues/277), [#283](https://github.com/gradientspaceai/gtengine-js/issues/283) |
| `Delaunay2.h` | deprecated `Update` | comments out its own `LogError` diagnostic | minor | n/a | [#277](https://github.com/gradientspaceai/gtengine-js/issues/277) |
| `Delaunay2.h` (L907), `Delaunay3.h` (L978) | `IntrinsicsVector2/3` construction | hardcoded `epsilon = 0` misclassifies exactly collinear or coplanar input | RC | fixed | [#391](https://github.com/gradientspaceai/gtengine-js/issues/391) |
| `Delaunay3.h` | `ProcessedVertex` (~L1486) | hashes and compares `location`, so duplicate detection can never match | RC | fixed | [#283](https://github.com/gradientspaceai/gtengine-js/issues/283) |
| `Delaunay3.h` | `GetNumVertices` | returns 0 for dimension 0/1/2 input | RC | fixed | [#283](https://github.com/gradientspaceai/gtengine-js/issues/283) |
| `Delaunay3.h` | comment (~L1021) | "IntrinsicsVector2{T}" should read `IntrinsicsVector3<T>` | doc | corrected | [#283](https://github.com/gradientspaceai/gtengine-js/issues/283) |
| `DistAlignedBoxAlignedBox.h` | result comment | "any choice of P0 and P1" is false on an overlap axis | doc | corrected | [#418](https://github.com/gradientspaceai/gtengine-js/issues/418) |
| `DistCircle2Circle2.h` | `operator()` | swaps arguments by radius while `DoQuery` writes its own first argument into `closest[j][0]` | minor | preserved | [#118](https://github.com/gradientspaceai/gtengine-js/issues/118) |
| `DistCircle3Circle3.h` | root finder | `p7` roundoff splits the double roots of `phi`, so mirror-symmetric configurations miss the true minimum | RC | fixed | [#331](https://github.com/gradientspaceai/gtengine-js/issues/331) |
| `DistCircle3Circle3.h` | candidate selection | `rRoots.size() > 1` used where the surviving-candidate count is meant; fabricates a pair at the origin | RC | fixed | [#331](https://github.com/gradientspaceai/gtengine-js/issues/331) |
| `DistCircle3Circle3.h` | candidate selection | unguarded `candidates[0]` read when nothing survives filtering | RC | fixed (assert) | [#331](https://github.com/gradientspaceai/gtengine-js/issues/331) |
| `DistCircle3Circle3.h` | candidate array | fixed 16 elements with no bound check | RC | fixed | [#331](https://github.com/gradientspaceai/gtengine-js/issues/331) |
| `DistCircle3Circle3.h` | `SCPolynomial` | dead private class | minor | n/a | [#331](https://github.com/gradientspaceai/gtengine-js/issues/331) |
| `DistCircle3Circle3.h` | `PrepareCircles` | the `normal[2] < 0` alignment is a no-op for normals orthogonal to z; anti-parallel normals take the polynomial path | RC | fixed | [#431](https://github.com/gradientspaceai/gtengine-js/issues/431) |
| `DistCircle3Circle3.h` | `sn = -p6/p7` | lands off the unit circle at the double roots that matter | RC | fixed | [#431](https://github.com/gradientspaceai/gtengine-js/issues/431) |
| `DistCircle3Circle3.h` | concentric/coaxial input | `phi` and `p6` are perfect squares, the sign-change bisection finds nothing and a default result is returned | RC | fixed | [#431](https://github.com/gradientspaceai/gtengine-js/issues/431) |
| `DistCircle3Circle3.h` | `p6 == p7 == 0` | throws "Unexpected degree for p6" for coaxial circles in parallel planes | WR | fixed | [#431](https://github.com/gradientspaceai/gtengine-js/issues/431) |
| `DistCircle3Circle3.h` | `DoQueryParallelPlanes` | projects with the rounded rotated normal, so coaxial circles report both closest points on the axis | RC | fixed | [#442](https://github.com/gradientspaceai/gtengine-js/issues/442) |
| `DistCircle3Circle3.h` | `PrepareCircles` | divides by an underflowed length for subnormal centre offsets | minor | preserved | [#431](https://github.com/gradientspaceai/gtengine-js/issues/431) |
| `DistCircle3Circle3.h` | result contract | at most two closest pairs; `numClosestPairs` is not symmetric under argument swap | minor | preserved | [#431](https://github.com/gradientspaceai/gtengine-js/issues/431) |
| `DistLine2Triangle2.h` | line 340 | the edge parameter denominator can round to exactly 0 while the normal components are nonzero; NaN result | RC | fixed | [#441](https://github.com/gradientspaceai/gtengine-js/issues/441) |
| `DistLine2Triangle2.h` | `NoCommonPoints` | the `abs(D)`-scaled distance is a dead assignment (recomputed before return) | minor | preserved | [#118](https://github.com/gradientspaceai/gtengine-js/issues/118) |
| `DistLine3CanonicalBox3.h` | `sqrDistance` accumulation | cancellation makes it negative for a flat box, so `sqrt` gives NaN | RC | fixed | [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) |
| `DistLine3Circle3.h` | `PDFSection422` | `tauHat` is missing the division by `a2` | RC | fixed | [#247](https://github.com/gradientspaceai/gtengine-js/issues/247) |
| `DistLine3Circle3.h` | header comment | `G'(t)` numerator given as `a1*a2`; it is `a1*a3` | doc | corrected | [#247](https://github.com/gradientspaceai/gtengine-js/issues/247) |
| `DistLine3Circle3.h` | `a3 == 0` path | NaN intercept and an inverted bisection bracket; throws "Invalid ordering of t-interval endpoints" | WR | fixed | [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) |
| `DistLine3Circle3.h` | `Finalize` | normalises a possibly-zero projection and reports the circle centre at distance 0 | RC | fixed | [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) |
| `DistLine3Circle3.h` | `Bisect`, `PDFSection422` back-substitution | for near-perpendicular lines the bracket collapses and `t = tau + s` then cancels every significant digit | RC | fixed | [#421](https://github.com/gradientspaceai/gtengine-js/issues/421), [#495](https://github.com/gradientspaceai/gtengine-js/issues/495) |
| `DistLine3Circle3.h` | `Execute`, `PDFSection421`, `PDFSection422` | `Dot(NxM, NxM)` underflows to exactly 0 while `NxM != 0`, so every output is NaN | RC | fixed | [#495](https://github.com/gradientspaceai/gtengine-js/issues/495) |
| `DistLine3OrientedBox3.h` | `operator()` | `closest[0]` is written in world space and then transformed again as if it were box-frame | RC | port already correct | [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) |
| `DistLine3OrientedBox3.h` | `operator()` | dead assignment to `result.closest[0]` | minor | preserved | [#187](https://github.com/gradientspaceai/gtengine-js/issues/187) |
| `DistLine3Rectangle3.h`, `DistLine3Triangle3.h` | `sqrDistance == invalid` | a squared distance can never equal the `-1` sentinel, so the test is dead | minor | preserved | [#187](https://github.com/gradientspaceai/gtengine-js/issues/187) |
| `DistLineSegment.h` | file comment | swaps the documented roles of `closest[0]` and `closest[1]` ("stoed" typo included) | doc | corrected | [#118](https://github.com/gradientspaceai/gtengine-js/issues/118) |
| `DistLineLine.h`, `DistLineRay.h`, `DistLineSegment.h`, `DistRayRay.h`, `DistRaySegment.h`, `DistSegmentSegment.h` | `det > 0` | the parallelism test is not reliable; exactly parallel pairs take the nonparallel branch | RC | preserved | [#418](https://github.com/gradientspaceai/gtengine-js/issues/418) |
| `DistOrientedBox3Cone3.h` | `operator()` | the documented finite-`hmax` precondition is never validated; an infinite cone gives an infeasible LCP | RC | fixed (assert) | [#298](https://github.com/gradientspaceai/gtengine-js/issues/298) |
| `DistOrientedBox3Cone3.h`, `LCPSolver.h` | 10-D LCP | the Hessian is a Gram matrix of five vectors in R^3, so Lemke returns a non-solution for ~1 in 400 configurations | RC | preserved | [#431](https://github.com/gradientspaceai/gtengine-js/issues/431) |
| `DistOrientedBox3Cone3.h` | angle sweep | exact `==` in the minimum assertion; silent `MAX_VALUE` on LCP failure; by-value `Z`; redundant containment test | minor | preserved | [#298](https://github.com/gradientspaceai/gtengine-js/issues/298) |
| `DistPoint3Cylinder3.h` | `operator()` | tests `height == max()` while `Cylinder3` uses `height = -1`, so infinite-cylinder queries throw | WR | fixed | [#187](https://github.com/gradientspaceai/gtengine-js/issues/187) |
| `DistPoint3Frustum3.h` | LF/UF edge cases | 2 of 10 assignments do not clamp the free far-edge coordinate; the closest point can lie outside the frustum | RC | fixed | [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) |
| `DistPoint3Tetrahedron3.h` | comments | stale "not unit length"; wrong trailing comment in `GetClosestRpmp` | doc | corrected | [#298](https://github.com/gradientspaceai/gtengine-js/issues/298) |
| `DistPointHyperellipsoid.h` | `Bisector` | the lower bracket rounds to exactly -1 near the smallest-extent axis and the result loses all precision | RC | preserved | [#424](https://github.com/gradientspaceai/gtengine-js/issues/424) |
| `DistRaySegment.h` | parallel branch, regions 1 and 5 | the ray parameter is never clamped to `s0 >= 0` | RC | fixed | [#126](https://github.com/gradientspaceai/gtengine-js/issues/126) |
| `DistSegment2Circle2.h` | contained segment | zeroes the whole result: `numClosestPairs = 0` and `distance = 0` for a segment strictly inside | RC | preserved | [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) |
| `DistSegment2Circle2.h` | boundary tests | half-open comparisons leave an off-segment parameter when an endpoint is exactly on the circle | minor | preserved | [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) |
| `DistSegment2Circle2.h` | interval branch chain | measure-zero gaps fall through to the final `else` | minor | preserved | [#187](https://github.com/gradientspaceai/gtengine-js/issues/187) |
| `DistSegmentSegment.h` | `ComputeIntersection` | the out-of-range ratio is replaced by `1/2` on a false premise; `ComputeRobust` then searches the wrong edge | RC | fixed | [#418](https://github.com/gradientspaceai/gtengine-js/issues/418) |
| `DistSegmentSegment.h` | `ComputeMinimumParameters` comment | gives `H(z)` with the wrong endpoint difference | doc | corrected | [#126](https://github.com/gradientspaceai/gtengine-js/issues/126) |
| `DistTriangle3CanonicalBox3.h` | degenerate triangle | divides by zero; the NaN fallback happens to reach the right answer | minor | preserved | [#298](https://github.com/gradientspaceai/gtengine-js/issues/298) |
| `ETManifoldMesh.h`, `TSManifoldMesh.h` | `Insert` | graceful nonmanifold rejection leaves phantom edges/faces referencing a triangle that is never added | RC | preserved | [#212](https://github.com/gradientspaceai/gtengine-js/issues/212) |
| `ETManifoldMesh.h`, `TSManifoldMesh.h` | `Insert` (throwing path) | the same entries are left pointing at a `unique_ptr` that dies during unwinding | RC | preserved | [#472](https://github.com/gradientspaceai/gtengine-js/issues/472) |
| `ETManifoldMesh.h` | `Insert` | the "nonmanifold" guard actually enforces consistent orientation and reports a misleading message | minor | preserved | [#212](https://github.com/gradientspaceai/gtengine-js/issues/212) |
| `ETManifoldMesh.h` | `GetBoundaryPolygon` | `std::map::operator[]` inserts a null entry into the map `GetBoundaryPolygons` is iterating | RC | fixed (assert) | [#212](https://github.com/gradientspaceai/gtengine-js/issues/212) |
| `ETNonmanifoldMesh.h` | `Insert` / `Remove` | a degenerate triangle aliases two edges; `Remove` throws after already erasing the shared edge | RC | preserved | [#179](https://github.com/gradientspaceai/gtengine-js/issues/179) |
| `ETNonmanifoldMesh.h` | `Insert` comment | "the (bad) triangle will not be part of the mesh" covers only `mTMap` | doc | corrected | [#179](https://github.com/gradientspaceai/gtengine-js/issues/179) |
| `ExtremalQuery3BSP.h` | `InsertArc` | arcs are compared only against node great circles, not the accumulated region; ~1% of icosahedron queries are wrong | RC | preserved | [#290](https://github.com/gradientspaceai/gtengine-js/issues/290) |
| `ExtremalQuery3BSP.h` | `GetTreeDepth` | reports the DFS stack high-water mark, not the tree depth | minor | preserved | [#290](https://github.com/gradientspaceai/gtengine-js/issues/290) |
| `ExtremalQuery3BSP.h` | construction | no check of `VETManifoldMesh::Insert` failure or a null `Edge::T[1]` | RC | fixed (assert) | [#290](https://github.com/gradientspaceai/gtengine-js/issues/290) |
| `FastGaussianBlur2.h`, `FastGaussianBlur3.h` | cleanup | 2D clears `mInput`/`mOutput` on exit, 3D leaves them dangling | minor | preserved | [#436](https://github.com/gradientspaceai/gtengine-js/issues/436) |
| `FastMarch.h` | `GetTimeExtremes` | redundantly re-tests one element | minor | preserved | [#52](https://github.com/gradientspaceai/gtengine-js/issues/52) |
| `FastMarch2.h`, `FastMarch3.h` | `ComputeTime` | the negative-discriminant fallback takes the larger neighbour time instead of the Godunov minimum | RC | preserved | [#439](https://github.com/gradientspaceai/gtengine-js/issues/439) |
| `FastMarch2.h`, `FastMarch3.h` | grid spacing | spacings are stored but never used by the numerical method | minor | preserved | [#439](https://github.com/gradientspaceai/gtengine-js/issues/439) |
| `FastMarch3.h` | `Initialize` (L232-341) | the six boundary faces are omitted from zero-speed marking, so face voxels index off-grid | RC | fixed | [#121](https://github.com/gradientspaceai/gtengine-js/issues/121) |
| `FPInterval.h` | `ProductLowerBound`, `ProductUpperBound` | wrong bounds in the both-straddle-zero branch | RC | preserved | [#75](https://github.com/gradientspaceai/gtengine-js/issues/75) |
| `FPInterval.h` | arithmetic | NaN endpoints lose the enclosure property | RC | preserved | [#476](https://github.com/gradientspaceai/gtengine-js/issues/476) |
| `GaussianElimination.h` | pivoting | denormal-magnitude matrices give NaN entries with `invertible = true` | RC | preserved | [#375](https://github.com/gradientspaceai/gtengine-js/issues/375) |
| `GaussNewtonMinimizer.h` | iteration | the unconditional `pCurrent = pNext` accepts error-increasing steps | minor | preserved | [#261](https://github.com/gradientspaceai/gtengine-js/issues/261) |
| `GenerateMeshUV.h` | `AssignBoundaryTextureCoordinatesSquare` | corner selection has no tolerance; a boundary vertex on a quarter mark is displaced a full edge | RC | preserved | [#262](https://github.com/gradientspaceai/gtengine-js/issues/262) |
| `GenerateMeshUV.h` | `mInteriorEdges` | a `std::set<Edge*>` iterated in pointer order makes the floating-point sums run-dependent | minor | fixed | [#262](https://github.com/gradientspaceai/gtengine-js/issues/262) |
| `GenerateMeshUV.h` | boundary detection | a closed surface leaves `mBoundaryStart == INT32_MAX` and indexes out of bounds | RC | fixed (assert) | [#262](https://github.com/gradientspaceai/gtengine-js/issues/262) |
| `GenerateMeshUV.h` | multithreaded path | unguarded `numV / mNumThreads` integer division | minor | n/a | [#262](https://github.com/gradientspaceai/gtengine-js/issues/262) |
| `GenerateMeshUV.h` | misc | `!= -1.0f` float literal against a `(Real)-1` sentinel; missing `break`; stale disk-boundary comment | minor | corrected | [#262](https://github.com/gradientspaceai/gtengine-js/issues/262) |
| `GMatrix.h` | comparison operators | mixed dimensions are not a strict weak ordering, so `std::map`/`std::set` keys are UB | minor | preserved | [#88](https://github.com/gradientspaceai/gtengine-js/issues/88) |
| `GradientAnisotropic2.h`, `GradientAnisotropic3.h` | `ComputeParameter` | padded coordinates fed to unpadded accessors: shifted window plus a read past the buffer | RC | fixed | [#122](https://github.com/gradientspaceai/gtengine-js/issues/122) |
| `GradientAnisotropic2.h`, `GradientAnisotropic3.h` | `ComputeParameter` | a constant image gives `mParameter = Infinity` and then NaN | RC | preserved | [#122](https://github.com/gradientspaceai/gtengine-js/issues/122) |
| `GradientAnisotropic2.h`, `GradientAnisotropic3.h` | `OnPreUpdate` | never calls the base implementation, so a masked Neumann filter never refreshes its mask border | minor | preserved | [#439](https://github.com/gradientspaceai/gtengine-js/issues/439) |
| `GradientAnisotropic2.h`, `GradientAnisotropic3.h` | divergence term | divided by the spacing once, not twice | minor | preserved | [#439](https://github.com/gradientspaceai/gtengine-js/issues/439) |
| `GVector.h`, `Vector.h` | `Length`, `Normalize` (robust) | rescaling by `1/maxAbsComp` overflows for a subnormal maximum, giving NaN | RC | preserved | [#370](https://github.com/gradientspaceai/gtengine-js/issues/370) |
| `HelmertTransformation7.h` | `A = sum(u v^T)` | computed and never read (dead Procrustes leftover) | minor | n/a | [#262](https://github.com/gradientspaceai/gtengine-js/issues/262) |
| `HelmertTransformation7.h` | `LogAssert(mNumPoints >= 7)` | confuses "7 parameters" with "7 points"; 3 correspondences suffice | minor | preserved | [#262](https://github.com/gradientspaceai/gtengine-js/issues/262) |
| `HelmertTransformation7.h` | `UpdateF` | `0.0` double literal in template code; stray `;;` | minor | corrected | [#262](https://github.com/gradientspaceai/gtengine-js/issues/262) |
| `Histogram.h` | bucket index (L80-87, 126-133, 171-178) | the maximum sample lands in bucket `B-2` about 7% of the time | minor | preserved | [#436](https://github.com/gradientspaceai/gtengine-js/issues/436) |
| `Histogram.h` | rescaled path | a subnormal sample range overflows `mult` and writes outside `mBuckets` | RC | preserved | [#436](https://github.com/gradientspaceai/gtengine-js/issues/436) |
| `Hyperellipsoid.h` | `ToCoefficients` | divides by `maxValue` with no zero check; `maxIndex` computed but never read | minor | preserved | [#217](https://github.com/gradientspaceai/gtengine-js/issues/217) |
| `Hyperplane.h` | `ComputeFromPoints` | passes `-1` as the SVD `multiplier`, which asserts for every `N != 3` | WR | fixed | [#217](https://github.com/gradientspaceai/gtengine-js/issues/217) |
| `Hyperplane.h` | `ComputeFromPoints` | `N = 2` constructs `SingularValueDecomposition(2,1,32)`, which itself asserts | WR | fixed | [#217](https://github.com/gradientspaceai/gtengine-js/issues/217) |
| `IEEEBinary16.h` | `Convert32To16` | NaNs with payload below 2^13 convert to infinity | RC | preserved | [#110](https://github.com/gradientspaceai/gtengine-js/issues/110) |
| `Image2.h` (L482-515), `Image3.h` (~L690-760) | `GetNeighborhood`, `GetCorners`, `GetFull` | `size_t` coordinates wrap to `SIZE_MAX` at a boundary, contradicting the adjacent documentation | RC | deviates | [#64](https://github.com/gradientspaceai/gtengine-js/issues/64) |
| `ImageUtility2.h` | `DrawEllipse` | infinite loop for `xExtent == yExtent == 0` | RC | fixed | [#443](https://github.com/gradientspaceai/gtengine-js/issues/443) |
| `ImageUtility2.h` | `GetSkeleton` | solid even-sided squares skeletonise to nothing | RC | preserved | [#443](https://github.com/gradientspaceai/gtengine-js/issues/443) |
| `ImageUtility2.h` | `DrawLine` | never updates `maxValue` in its `dy > maxValue` branch | minor | preserved | [#129](https://github.com/gradientspaceai/gtengine-js/issues/129) |
| `ImageUtility2.h`, `ImageUtility3.h` | neighbourhood helpers | read neighbours with no range test, relying wholly on the zero-boundary precondition | minor | preserved | [#129](https://github.com/gradientspaceai/gtengine-js/issues/129) |
| `ImageUtility3.h` | `Dilate` | the innermost loop starts at `i0 = 1`, so `x = 0` voxels never act as dilation sources | RC | fixed | [#129](https://github.com/gradientspaceai/gtengine-js/issues/129) |
| `ImageUtility3.h` | `Close<N>` | uncompilable dead code (2D static_assert, nonexistent two-argument `Image3` constructor) | minor | fixed | [#129](https://github.com/gradientspaceai/gtengine-js/issues/129) |
| `IncrementalDelaunay2.h` | `GetHull` | dereferences `edges.begin()` on an empty edge map | RC | fixed | [#290](https://github.com/gradientspaceai/gtengine-js/issues/290) |
| `IncrementalDelaunay2.h` | `GetHull` | the unbounded boundary walk never returns to the start for collinear input | RC | fixed | [#290](https://github.com/gradientspaceai/gtengine-js/issues/290) |
| `IncrementalDelaunay2.h` | `RetriangulateBoundaryRemovalPolygon` | the `numPolygon == 2` branch tests a condition that can never hold | minor | fixed | [#290](https://github.com/gradientspaceai/gtengine-js/issues/290) |
| `IncrementalDelaunay2.h` | comments | claim `GetNumVertices`/`GetVertices` exclude the supervertices; they do not | doc | corrected | [#290](https://github.com/gradientspaceai/gtengine-js/issues/290) |
| `IncrementalDelaunay2.h` | `DoEarClipping` | skips `RPPolygon::Remove` on its early exit, leaving `GetNumActive()` one too large | minor | preserved | [#290](https://github.com/gradientspaceai/gtengine-js/issues/290) |
| `InscribedFixedAspectRectInQuad.h` | `Execute` | three edge normals in one quadrant share a case index and the `alpha` assertion fires on a solvable quad | WR | preserved | [#395](https://github.com/gradientspaceai/gtengine-js/issues/395) |
| `InscribedFixedAspectRectInQuad.h` | `Execute` | a degenerate feasible interval gives "Unexpected interval intersection type" | WR | preserved | [#395](https://github.com/gradientspaceai/gtengine-js/issues/395) |
| `InscribedFixedAspectRectInQuad.h` | comment | `quad[(i + 1) % 3]` for a quadrilateral | doc | corrected | [#101](https://github.com/gradientspaceai/gtengine-js/issues/101) |
| `Integration.h` | `Romberg` | latent `int32_t` overflow at `order >= 17` and `1u << (n-1)` allocation growth | minor | preserved | [#476](https://github.com/gradientspaceai/gtengine-js/issues/476) |
| `IntpAkimaUniform2.h`, `IntpAkimaUniform3.h` | `GetFXY` (and 3D `GetFXZ`, `GetFYZ`, `GetFXYZ`) | max-boundary stencils reuse the min-boundary coefficients on a reversed index sequence, so the sign is inverted | RC | fixed | [#58](https://github.com/gradientspaceai/gtengine-js/issues/58) |
| `IntpBicubic2.h`, `IntpBilinear2.h`, `IntpTricubic3.h`, `IntpTrilinear3.h` | `operator()` | documented input clamping clamps only the cell index; the fractional coordinate extrapolates | RC | preserved | [#69](https://github.com/gradientspaceai/gtengine-js/issues/69) |
| `IntpBicubic2.h`, `IntpTricubic3.h` | `mBlend[0][3]` initialiser | stray double semicolon | minor | corrected | [#69](https://github.com/gradientspaceai/gtengine-js/issues/69) |
| `IntpBilinear2.h` | constructor comment | claims a 3x3 block minimum where the assertion and algorithm need 2x2 | doc | corrected | [#69](https://github.com/gradientspaceai/gtengine-js/issues/69) |
| `IntpBSplineUniform.h` | `ComputePowers` | writes `powerDSDT[1]` unconditionally after `resize(degree+1)`; overruns for degree 0 | RC | fixed | [#135](https://github.com/gradientspaceai/gtengine-js/issues/135) |
| `IntpBSplineUniform.h` | `A` extraction | reads `Q[k][col]` past the degree after leading-zero elimination | minor | fixed | [#135](https://github.com/gradientspaceai/gtengine-js/issues/135) |
| `IntpBSplineUniform.h` | `ComputeBlendingMatrix` | dead local `sm1` with a comment claiming it is used | minor | dropped | [#135](https://github.com/gradientspaceai/gtengine-js/issues/135) |
| `IntpLinearNonuniform2.h`, `IntpLinearNonuniform3.h` | `operator()` | `GetIndices`' failure flag is discarded and the zero-initialised index tuple is used | RC | deviates | [#135](https://github.com/gradientspaceai/gtengine-js/issues/135) |
| `IntpQuadraticNonuniform2.h` | all mesh accessors | every failure flag is discarded; zero index triples blend unrelated samples | RC | fixed | [#337](https://github.com/gradientspaceai/gtengine-js/issues/337) |
| `IntpQuadraticNonuniform2.h` | `ProcessTriangles` | ignores `Inscribe`'s failure return; a degenerate triangle gets centre (0,0) | minor | preserved | [#337](https://github.com/gradientspaceai/gtengine-js/issues/337) |
| `IntpQuadraticNonuniform2.h` | `operator()` fallback | consults barycentrics even when `ComputeBarycentrics` failed; a degenerate subtriangle always wins | RC | preserved | [#337](https://github.com/gradientspaceai/gtengine-js/issues/337) |
| `IntpQuadraticNonuniform2.h` | `ComputeCoefficients` | unguarded divisions; stale "circumscribing circle" comment | minor | preserved | [#337](https://github.com/gradientspaceai/gtengine-js/issues/337) |
| `IntpThinPlateSpline2.h`, `IntpThinPlateSpline3.h` | WARNING comment (L14-15) | the invariance claim is wrong in both directions | doc | corrected | [#191](https://github.com/gradientspaceai/gtengine-js/issues/191) |
| `IntpThinPlateSpline2.h`, `IntpThinPlateSpline3.h` | domain transform | unguarded division by a zero coordinate range gives NaN coordinates | minor | preserved | [#191](https://github.com/gradientspaceai/gtengine-js/issues/191) |
| `IntpThinPlateSpline2.h`, `IntpThinPlateSpline3.h` | `ComputeFunctional` | discontinuous at `lambda = 0`; values are not comparable across smoothing settings | minor | preserved | [#191](https://github.com/gradientspaceai/gtengine-js/issues/191) |
| `IntpVectorField2.h` | `operator()` | the `&&` short-circuit leaves `output[1]` stale when the x-interpolation fails | RC | fixed | [#337](https://github.com/gradientspaceai/gtengine-js/issues/337) |
| `IntrAlignedBox3Cone3.h` | `BoxFullyInConeSlab` | stale adjacency bits survive into a later clipping query on the same object; false negatives | RC | fixed | [#301](https://github.com/gradientspaceai/gtengine-js/issues/301) |
| `IntrAlignedBox3Cone3.h`, `IntrOrientedBox3Cone3.h` | header TODO | still written for the retired `maxHeight = max()` infinite-cone representation | doc | preserved | [#334](https://github.com/gradientspaceai/gtengine-js/issues/334) |
| `IntrAlignedBox3OrientedBox3.h` vs `IntrOrientedBox3OrientedBox3.h` | parallel-pair cutoff | one file uses `>= cutoff`, the other `> cutoff`, for the same separating-axis algorithm | minor | preserved | [#450](https://github.com/gradientspaceai/gtengine-js/issues/450) |
| `IntrAlignedBox3Sphere3.h` | `FIQuery::operator()` | `contactPoint += boxCenter` is unconditional, so a no-contact result reports the box centre | minor | preserved | [#250](https://github.com/gradientspaceai/gtengine-js/issues/250) |
| `IntrAlignedBox3Sphere3.h` | `DoQueryRayRoundedFace`, `DoQuery` | first-wins probe selection over pieces of the Minkowski sum reports contacts late or misses them entirely | RC | fixed | [#458](https://github.com/gradientspaceai/gtengine-js/issues/458), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrAreaEllipse2Ellipse2.h` | `mZero`, `mOne`, `mTwo`, `mPi`, `mTwoPi` | declared with no constructor and never assigned; every area reads indeterminate values | RC | fixed | [#301](https://github.com/gradientspaceai/gtengine-js/issues/301) |
| `IntrAreaEllipse2Ellipse2.h` | class comment | "axes are not required to be normalized" is false for the polar angles, the `Area4` ordering and the full-ellipse area | doc | fixed | [#301](https://github.com/gradientspaceai/gtengine-js/issues/301), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrCanonicalBox3Cylinder3.h` | `DoQueryNoZeros`, `(U1, -D)` block | sign typo puts a segment endpoint on the wrong box edge; false negatives | RC | fixed | [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrCapsule3Capsule3.h` | `operator()` | uses the non-robust segment-segment distance, so a strictly contained capsule is reported as no intersection | RC | preserved | [#455](https://github.com/gradientspaceai/gtengine-js/issues/455) |
| `IntrCircle2Circle2.h` | TI vs FI | TI tests solid disks, FI intersects the curves; nested circles give TI `true` and FI `false` | doc | preserved | [#450](https://github.com/gradientspaceai/gtengine-js/issues/450) |
| `IntrConvexMesh3Plane3.h` | `GetIntersectionPolygon` | the coplanar-face predecessor map is read in index order instead of traversed as a cycle | RC | fixed | [#301](https://github.com/gradientspaceai/gtengine-js/issues/301) |
| `IntrConvexMesh3Plane3.h` | map lookups | `eiVMap[EdgeKey(...)]` default-inserts 0; a `find()` result is dereferenced without an `end()` check | RC | fixed (assert) | [#301](https://github.com/gradientspaceai/gtengine-js/issues/301) |
| `IntrConvexPolygonHyperplane.h` | TI `NEGATIVE_SIDE_VERTEX` comment | says "vertex or an edge"; the branch handles only the single-vertex case | doc | corrected | [#250](https://github.com/gradientspaceai/gtengine-js/issues/250) |
| `IntrCylinder3Cylinder3.h`, `IntrHalfspace3Cylinder3.h`, `IntrRay3Cylinder3.h`, `IntrSegment3Cylinder3.h`, `IntrTriangle3Cylinder3.h` | `operator()` | no `IsFinite()` guard, so the `height = -1` infinite sentinel silently yields a negative half-height | RC | fixed (assert) | [#197](https://github.com/gradientspaceai/gtengine-js/issues/197), [#206](https://github.com/gradientspaceai/gtengine-js/issues/206), [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrCylinder3Cylinder3.h` | constructor comment | says `phi[j] = pi*j/numPhi`; the code uses `GTE_C_HALF_PI` | doc | corrected | [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrEllipse2Ellipse2.h` | `CaseE4NotZero`, `divisor == 0` | both symmetric points are written to the same slot, so one real intersection is lost | RC | fixed | [#250](https://github.com/gradientspaceai/gtengine-js/issues/250) |
| `IntrEllipse2Ellipse2.h` | TIQuery | drops the `c_i = 0` terms of `f(s)`, losing the two pole critical points | RC | fixed | [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) |
| `IntrEllipse2Ellipse2.h` | quartic root finding | concentric congruent ellipses rotated exactly 90 degrees report no intersection | RC | preserved | [#250](https://github.com/gradientspaceai/gtengine-js/issues/250) |
| `IntrEllipse2Ellipse2.h` | TI bracket, quartic | the bracket collapses onto its pole for nearly concentric ellipses; near-double roots are lost or invented | RC | preserved | [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) |
| `IntrEllipse2Ellipse2.h` | `CaseE4ZeroE2NotZeroE3Zero`, `GetRoots` | `<=` where the siblings use `<`; dead `fval = F(s)` stores | minor | preserved | [#250](https://github.com/gradientspaceai/gtengine-js/issues/250) |
| `IntrEllipsoid3Ellipsoid3.h` | valid-pair analysis, `d0 > d1 = d2` | `param[1].second += param[0].second` folds the coefficient of the distinct eigenvalue into the repeated one (counting it twice) and drops `param[2].second`; `f(s)` is built with the wrong coefficients and the classification is wrong | RC | fixed | [#503](https://github.com/gradientspaceai/gtengine-js/issues/503) |
| `IntrEllipsoid3Ellipsoid3.h` | `GetRoots` | the ad-hoc `epsilon = 0.001` bracketing asserts fire for ordinary close-centre input | WR | preserved | [#255](https://github.com/gradientspaceai/gtengine-js/issues/255), [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) |
| `IntrEllipsoid3Ellipsoid3.h` | `Matrix3x3 D0` | dead store | minor | dropped | [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrHalfspace2Polygon2.h` | FIQuery | returns `intersect = true` with an empty polygon when the input lies entirely inside | RC | preserved | [#139](https://github.com/gradientspaceai/gtengine-js/issues/139) |
| `IntrHalfspace3Cylinder3.h` | `root` computation | `max((T)1, (T)1 - x*x)` is always 1, so the documented `sqrt(1 - Dot(N,W)^2)` is never computed; false positives | RC | fixed | [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrHalfspace3Segment3.h` | FIQuery, case (1,1,0) | sets `numPoints = 1` and returns only the crossing point where its own table documents the clipped sub-segment | RC | preserved | [#139](https://github.com/gradientspaceai/gtengine-js/issues/139) |
| `IntrIntervals.h` | dynamic FIQuery (~L485) | the left-approach branch uses the moved left endpoint instead of the contacting right endpoint | RC | fixed | [#62](https://github.com/gradientspaceai/gtengine-js/issues/62) |
| `IntrIntervals.h` | ~L446 | provable dead store | minor | preserved | [#62](https://github.com/gradientspaceai/gtengine-js/issues/62) |
| `IntrLine2Ray2.h`, `IntrLine2Segment2.h` | zero-length input | a zero-direction supporting line makes `IntrLine2Line2` report "same line" even for an off-line point | RC | preserved | [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrLine2Segment2.h` | parameter test | `line1Parameter[0] >= 0 && line1Parameter[1] <= 1` is redundant; both entries hold the same value | minor | preserved | [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrLine2SegmentMesh2.h` | record construction | copies the `+-max()` sentinel into `lineParameter` while storing the real endpoint as `point` | RC | preserved | [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) |
| `IntrLine3Capsule3.h` | hemisphere roots | `intersect` is never set when only one hemisphere root is accepted | RC | fixed | [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) |
| `IntrLine3Capsule3.h` | cap-junction plane | overlapping acceptance regions accept a junction root twice and discard the true far endpoint | RC | fixed | [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) |
| `IntrLine3Capsule3.h` | zero-length capsule segment | FI returns `[-radius, +radius]` regardless of the line | minor | preserved | [#200](https://github.com/gradientspaceai/gtengine-js/issues/200) |
| `IntrLine3Cone3.h` | `CaseC2NotZeroDiscrZero` | the vertex test uses only the U-component of `(P - V) + t*U = 0` | RC | fixed | [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrLine3Cone3.h` | `CaseC2NotZeroDiscrZero` | the vertex branch is gated on an exact floating-point equality that essentially never holds | RC | fixed | [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrLine3Cone3.h` | `CaseC2NotZeroDiscrPos` Block 3 | assumes `c2 > 0`; a through-vertex line has `discr` as a cancelling difference and is misclassified | RC | fixed | [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrLine3Cone3.h` | `Result::Convert` | the member templates cannot compile (`QFNumber` has no conversion to `Real`) | minor | n/a | [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrLine3Cone3.h` | `isRayNegative` doc | contradicts itself | doc | corrected | [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrLine3Cylinder3.h` | `operator()` | no infinite-cylinder branch at all; assumes a unit direction | minor | preserved | [#197](https://github.com/gradientspaceai/gtengine-js/issues/197), [#200](https://github.com/gradientspaceai/gtengine-js/issues/200) |
| `IntrLine3Rectangle3.h`, `IntrSegment3Rectangle3.h` | `Result::rectCoord` | a 3-tuple for a 2-coordinate rectangle; component 2 is always zero | minor | preserved | [#141](https://github.com/gradientspaceai/gtengine-js/issues/141) |
| `IntrLine3Rectangle3.h` | header comment | `X = C + sum_{i=0}^2 s[i]*W[i]` should be `i = 0..1`; branches say "Line and triangle are parallel" | doc | corrected | [#141](https://github.com/gradientspaceai/gtengine-js/issues/141) |
| `IntrOrientedBox2OrientedBox2.h` | TI vs FI | TI uses closed-box separation, FI requires a strictly positive vertex distance; touching boxes disagree | minor | preserved | [#141](https://github.com/gradientspaceai/gtengine-js/issues/141) |
| `IntrOrientedBox2Sector2.h` | boundary clipping | `polygon = std::move(hpResult.polygon)` wipes the working polygon when a clip was a no-op | RC | fixed | [#200](https://github.com/gradientspaceai/gtengine-js/issues/200) |
| `IntrOrientedBox2Sector2.h` | wedge model | modelled as two halfplanes, which is wrong for half-angle > pi/2 (including `Sector2`'s own default) | RC | preserved | [#200](https://github.com/gradientspaceai/gtengine-js/issues/200) |
| `IntrPlane3Cylinder3.h` | CIRCLE/ELLIPSE label | chosen by an exact extent comparison | minor | preserved | [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrPlane3Cylinder3.h` | `Result` docs | refer to a nonexistent `type = NONE` and to "all zero members" for value-initialised members | doc | corrected | [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrPlane3Plane3.h` | parallel test | `abs(Dot(N0,N1)) >= 1` is exact; identical planes are frequently classified transverse with `invDet ~ 1e16` | RC | preserved | [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrPlane3Triangle3.h` | FIQuery | uses a raw uninitialised `Real s[3]` where TIQuery uses an initialised array | minor | corrected | [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrRay2Arc2.h`, `IntrSegment2Arc2.h` | `operator()` | solid-disk clipping feeds an off-circle point into `Arc2::Contains`, reporting hits that are not on the arc | RC | fixed | [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrRay2Ray2.h` | TI vs FI | opposite collinear rays with a shared origin give `numIntersections` 1 (TI) vs 2 (FI) | minor | preserved | [#200](https://github.com/gradientspaceai/gtengine-js/issues/200), [#455](https://github.com/gradientspaceai/gtengine-js/issues/455) |
| `IntrSegment2AlignedBox2.h` | `cdeParameter` | self-assignment no-op | minor | dropped | [#203](https://github.com/gradientspaceai/gtengine-js/issues/203) |
| `IntrSegment2Arc2.h` | `parameter[0] = 0` | duplicated assignment | minor | preserved | [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrSegment2OrientedBox2.h` | FIQuery | mixes box-frame and world-frame coordinates when building `result.point[i]` | RC | fixed | [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrSegment2OrientedBox2.h` | `cdeParameter` | self-assignment no-op | minor | dropped | [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrSegment2OrientedBox2.h` | degenerate segment | a point segment inside the box reports `numIntersections = 2`; the aligned sibling reports 1 | minor | preserved | [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrSegment2Segment2.h` | FIQuery | `segment1Parameter = overlap - t` has the wrong sign for antiparallel collinear segments | RC | fixed | [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) |
| `IntrSegment2Segment2.h` | `Result` documentation | the point identity and the two ascending-order claims are mutually incompatible | doc | preserved | [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) |
| `IntrSegment2Segment2.h` | `Exact` | `line0Parameter` index mix; a zero-length segment is classified collinear | minor | preserved | [#203](https://github.com/gradientspaceai/gtengine-js/issues/203) |
| `IntrSegment3Cylinder3.h` | `GetCenteredForm` | no zero-length-segment guard; a degenerate segment inside the cylinder returns `intersect = true` with NaN | RC | preserved | [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrSegment3Ellipsoid3.h` | TIQuery | a segment fully inside the ellipsoid reports no intersection, contradicting the FI query in the same header | RC | fixed | [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrSegment3Ellipsoid3.h` | comments | stale `a0`, `a1`, `-a1/2`, `a3*e` descriptions | doc | corrected | [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrSegment3Sphere3.h` | TIQuery | a segment strictly inside the sphere reports no intersection | RC | fixed | [#203](https://github.com/gradientspaceai/gtengine-js/issues/203) |
| `IntrSphere3Cone3.h` | FI query | `result.point` omits the cone vertex `V` | RC | fixed | [#307](https://github.com/gradientspaceai/gtengine-js/issues/307) |
| `IntrSphere3Cone3.h` | `DoQueryFiniteCone` | calls `GetMaxHeight()` twice redundantly | minor | preserved | [#307](https://github.com/gradientspaceai/gtengine-js/issues/307) |
| `IntrSphere3Sphere3.h` | FIQuery type-4 branch | the internal-tangency contact point is the antipode of the true one | RC | fixed | [#203](https://github.com/gradientspaceai/gtengine-js/issues/203) |
| `IntrSphere3Triangle3.h` | `intersectionType = +1` | returns the sphere centre rather than the surface contact point | minor | preserved | [#203](https://github.com/gradientspaceai/gtengine-js/issues/203) |
| `IntrSphere3Triangle3.h` | arbitrary-precision path | the `closest[j]` loop runs out of bounds | RC | n/a | [#203](https://github.com/gradientspaceai/gtengine-js/issues/203) |
| `IntrTetrahedron3Tetrahedron3.h` | edge-edge phase | `abs(Dot(E0,E1)) < cutoff` compares unnormalised dots against a cosine cutoff, skipping the whole phase for long edges | RC | fixed | [#307](https://github.com/gradientspaceai/gtengine-js/issues/307) |
| `IntrTetrahedron3Tetrahedron3.h` | edge-edge phase | the separation test is a plane-side test, not projection-interval disjointness; argument-order dependent | RC | fixed | [#307](https://github.com/gradientspaceai/gtengine-js/issues/307) |
| `IntrTetrahedron3Tetrahedron3.h` | face-normal phase | undocumented positive-orientation precondition; a negatively oriented tetrahedron gives inward normals | doc | preserved | [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrTriangle2Triangle2.h` | `WhichSide` | relies on C++ int truthiness via `--negative` | minor | preserved | [#307](https://github.com/gradientspaceai/gtengine-js/issues/307) |
| `IntrTriangle3CanonicalBox3.h`, `IntrTriangle3*Box3.h` | comment | "at most 7 vertices" for the clipped polygon; the bound is 9 and octagons occur | doc | corrected | [#307](https://github.com/gradientspaceai/gtengine-js/issues/307), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrTriangle3Cylinder3.h` | `DiskOverlapsPolygon` (L451-487) | a degenerate projected polygon gives `positive == negative == 0`, reported as containment at any distance | RC | fixed | [#206](https://github.com/gradientspaceai/gtengine-js/issues/206), [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) |
| `IntrTriangle3Cylinder3.h` | comments | "on the tob of the slab"; case 4c/4d parameterises edge `<V2,V0>` the other way round | doc | corrected | [#206](https://github.com/gradientspaceai/gtengine-js/issues/206) |
| `IntrTriangle3Triangle3.h` | moving-triangle overloads | `fabs(Dot(N0,N1)) < 1` compares an area-scaled quantity against 1; the wrong branch is taken both ways | RC | fixed | [#334](https://github.com/gradientspaceai/gtengine-js/issues/334) |
| `IntrTriangle3Triangle3.h` | contact configuration | `ContactSide` and the contact-time configurations are maintained through ~150 lines and never read | minor | preserved | [#334](https://github.com/gradientspaceai/gtengine-js/issues/334) |
| `IntrTriangle3Triangle3.h` | `GetCoplanarIntersection`, `IntersectsSegment`, `ContainsPoint` | divide by `normal[lookup[2]]` with no guard; a zero-area triangle emits NaN/Infinity | RC | preserved | [#334](https://github.com/gradientspaceai/gtengine-js/issues/334) |
| `IntrTriangle3Triangle3.h` | file comment | documents a moving-triangle FI query the header does not contain | doc | corrected | [#307](https://github.com/gradientspaceai/gtengine-js/issues/307) |
| `IsPlanarGraph.h` | `operator()` | after flagging invalid vertex indices the code proceeds into position lookups with them | RC | fixed | [#44](https://github.com/gradientspaceai/gtengine-js/issues/44) |
| `IsPlanarGraph.h` | `InvalidSegmentIntersection` | the collinear branch parameterises by the first segment, so a zero-length first segment collapses to "no overlap" | RC | preserved | [#472](https://github.com/gradientspaceai/gtengine-js/issues/472) |
| `LCPSolver.h` | `Solve` | constructing the dynamic solver with `n <= 0` leaves null members and `Solve` dereferences them | RC | fixed | [#76](https://github.com/gradientspaceai/gtengine-js/issues/76) |
| `LCPSolver.h` | Lemke pivot | a subnormal pivot turns a provably infeasible problem into `HAS_NONTRIVIAL_SOLUTION` with infinite `z` | RC | preserved | [#476](https://github.com/gradientspaceai/gtengine-js/issues/476) |
| `LDLTDecomposition.h` | `BlockLDLTDecomposition::Convert` | asserts `GetSize() == NumBlocks` where block vectors have `BlockSize` components | WR | fixed | [#209](https://github.com/gradientspaceai/gtengine-js/issues/209) |
| `LDLTDecomposition.h` | documentation | claims positive-definite input and positive `D`; `Factor` fails only on an exactly zero pivot | doc | preserved | [#209](https://github.com/gradientspaceai/gtengine-js/issues/209) |
| `LevenbergMarquardtMinimizer.h` | `DoIteration` | forms `mNegJTF` from `F` at the previous rejected candidate while `J` is at `pCurrent` | RC | fixed | [#261](https://github.com/gradientspaceai/gtengine-js/issues/261) |
| `LevenbergMarquardtMinimizer.h` | Gauss-Newton fallback | uses the maximally inflated lambda and never divides it back down; comment and code disagree | minor | preserved | [#261](https://github.com/gradientspaceai/gtengine-js/issues/261) |
| `LevenbergMarquardtMinimizer.h`, `GaussNewtonMinimizer.h` | `Result::numIterations` | is the raw loop counter, so an exhausted loop reports `maxIterations + 1` | minor | preserved | [#261](https://github.com/gradientspaceai/gtengine-js/issues/261) |
| `LieGroupsAlgebras.h` | `LieSO3::Log`, angle = pi | scales by `pi/sqrt(2)` where the correct factor is `pi`; `Exp(Log(Y)) != Y` for every pi rotation | RC | fixed | [#313](https://github.com/gradientspaceai/gtengine-js/issues/313) |
| `LieGroupsAlgebras.h` | `LieSO3::Log`, near pi | the generic branch amplifies round-off by ~1e8 when the trace lands just above -1 | RC | preserved | [#313](https://github.com/gradientspaceai/gtengine-js/issues/313) |
| `LinearSystem.h` | `SolveTridiagonal`, `SolveConstantTridiagonal` | never validate `N`; `N = 0` wraps the temporary size to `SIZE_MAX` | RC | fixed (assert) | [#261](https://github.com/gradientspaceai/gtengine-js/issues/261) |
| `LinearSystem.h` | `SolveSymmetricCG` | recomputes the loop-invariant `Dot(N, B, B)` every iteration | minor | fixed | [#261](https://github.com/gradientspaceai/gtengine-js/issues/261) |
| `LinearSystem.h` | `SolveSymmetricCG` | with `B = 0` computes `alpha = 0/0` and overwrites the exact solution with NaN | RC | preserved | [#261](https://github.com/gradientspaceai/gtengine-js/issues/261) |
| `LogEstimate.h` | `GetLogEstimateMaxError` | returns the log2 bound, loose by a factor of ~1.443, unlike `ExpEstimate.h` | minor | preserved | [#57](https://github.com/gradientspaceai/gtengine-js/issues/57) |
| `Matrix4x4.h` | `MakePerspectiveProjection`, `GTE_USE_VEC_MAT` | `M(3,0)` and `M(3,1)` read entries before they are assigned; the last row is wrong | RC | n/a | [#160](https://github.com/gradientspaceai/gtengine-js/issues/160) |
| `Mesh.h` | `ComputeIndices`, SPHERE | the south-pole fan uses the `numCols` stride where rows are `numCols + 1` wide | RC | fixed | [#220](https://github.com/gradientspaceai/gtengine-js/issues/220) |
| `Mesh.h` | `ComputeIndices`, SPHERE | the second pole fan is wound opposite to the body (inward normal) | RC | fixed | [#240](https://github.com/gradientspaceai/gtengine-js/issues/240) |
| `Mesh.h` family | `allowUpdateFrame` | re-checked only inside `if (!mTCoords)`, so a client-supplied tcoord channel skips the stricter requirement | minor | preserved | [#268](https://github.com/gradientspaceai/gtengine-js/issues/268) |
| `MeshCurvature.h` | principal directions | both candidate vectors are zero at an exact umbilic, so both directions come back `(0,0,0)` | RC | preserved | [#412](https://github.com/gradientspaceai/gtengine-js/issues/412) |
| `MeshCurvature.h` | `maxAbs < singularityThreshold` | never fires at the documented threshold of zero | minor | preserved | [#240](https://github.com/gradientspaceai/gtengine-js/issues/240) |
| `MeshStaticManifold2.h` | `GetAdjacentTriangles` | returns the neighbour across the edge opposite `v0`, not the queried edge; documented case 3 is unreachable | RC | fixed | [#66](https://github.com/gradientspaceai/gtengine-js/issues/66), [#361](https://github.com/gradientspaceai/gtengine-js/issues/361) |
| `MeshStaticManifold3.h` | `GetAdjacentTetrahedra` | reads the 2D tuple indices on a 5-wide record, returning a vertex index as a tetrahedron index | RC | fixed | [#66](https://github.com/gradientspaceai/gtengine-js/issues/66), [#361](https://github.com/gradientspaceai/gtengine-js/issues/361) |
| `MeshStaticManifold2.h` | `GetComponents` comment | describes a `range[]` output the function does not have | doc | corrected | [#361](https://github.com/gradientspaceai/gtengine-js/issues/361) |
| `MinimalCycleBasis.h` | detach-time `DepthFirstSearch` | visited flags persist between searches, flattening nested cycle forests at depth >= 2 | RC | fixed | [#310](https://github.com/gradientspaceai/gtengine-js/issues/310) |
| `MinimalCycleBasis.h` | `GetClockwiseMost` | the non-convex branch uses `< 0` where five siblings and the published pseudocode use `<= 0` | minor | preserved | [#310](https://github.com/gradientspaceai/gtengine-js/issues/310) |
| `MinimalCycleBasis.h` | `Extract`, `ExtractConnectedComponents` | isolated vertices go unreported and a lone-edge component is reported as nothing | minor | preserved | [#310](https://github.com/gradientspaceai/gtengine-js/issues/310) |
| `Minimize1.h` | `GetBracketedMinimum` | at equal-value endpoints the exact `tv == tm` test is missed and the bracket collapses | RC | fixed | [#298](https://github.com/gradientspaceai/gtengine-js/issues/298) |
| `MinimizeN.h` | `GetMinimum` | the stray `mDConjIndex = 0` destroys Powell's direction-set update; the method stalls | RC | fixed | [#146](https://github.com/gradientspaceai/gtengine-js/issues/146) |
| `MinimizeN.h` | `GetMinimum` | `mDConjIndex` is never reset, so a second call on the same object is broken from its first iteration | RC | fixed | [#146](https://github.com/gradientspaceai/gtengine-js/issues/146) |
| `MinimizeN.h` | algorithm | basic Powell can stall in a proper subspace while reporting convergence | minor | preserved | [#478](https://github.com/gradientspaceai/gtengine-js/issues/478) |
| `MinimumAreaBox2.h` | `RemoveCollinearPoints` | a duplicate point's zero-length edge also removes the next genuine corner | RC | fixed | [#328](https://github.com/gradientspaceai/gtengine-js/issues/328) |
| `MinimumAreaBox2.h` | `dimension == 1` branch | `imin`/`imax` are left at 0 although the degenerate line's origin is `points[hull[0]]` | RC | fixed | [#328](https://github.com/gradientspaceai/gtengine-js/issues/328) |
| `MinimumAreaBox2.h` | degenerate branches | return stale `mArea`/`mSupportIndices` from a previous query | RC | fixed | [#328](https://github.com/gradientspaceai/gtengine-js/issues/328) |
| `MinimumAreaBox2.h` | convex-polygon overload | never assigns `mNumPoints`/`mPoints`, so `GetPoints()` is null or stale | RC | fixed | [#402](https://github.com/gradientspaceai/gtengine-js/issues/402) |
| `MinimumAreaBox2.h` | overload 4 | contradicts overload 3 on empty polygon indices; stale accessor comments | minor | corrected | [#328](https://github.com/gradientspaceai/gtengine-js/issues/328) |
| `MinimumAreaCircle2.h`, `MinimumVolumeSphere3.h` | exception-trap fallback | `GetContainer` is called with the unique-point count but the full point array | RC | fixed | [#286](https://github.com/gradientspaceai/gtengine-js/issues/286) |
| `MinimumAreaCircle2.h`, `MinimumVolumeSphere3.h` | `UpdateSupport{2,3,4}` | rewrite the support set unconditionally, so a rejected update hides a point and the result does not contain it | RC | preserved | [#399](https://github.com/gradientspaceai/gtengine-js/issues/399) |
| `MinimumSpanningTree.h` | `ExtractMinimumSpanningTree` | an empty edge list writes `records[0]` on a zero-length vector | RC | fixed | [#74](https://github.com/gradientspaceai/gtengine-js/issues/74) |
| `MinimumSpanningTree.h` | header comment | describes a spanning tree with one sentinel; a disconnected graph gives a forest | doc | corrected | [#74](https://github.com/gradientspaceai/gtengine-js/issues/74) |
| `MinimumSpanningTree.h` | back edges | filtered by remapped index before `ConvertToOriginalIndices`, so the reported pairs are unordered in caller labels | doc | preserved | [#472](https://github.com/gradientspaceai/gtengine-js/issues/472) |
| `MinimumVolumeBox3FloatingPoint.h` | `ComputeConvexHull`, dimension 2 | the Newell-normal loop drops the wrap-around term; a triangular hull gives an exactly zero normal | RC | fixed | [#352](https://github.com/gradientspaceai/gtengine-js/issues/352) |
| `MinimumVolumeBox3FloatingPoint.h` | `ComputeVolume` | assumes a hull-edge vertex realises the `axis[0]`/`axis[1]` minima; false in floating point | RC | fixed | [#405](https://github.com/gradientspaceai/gtengine-js/issues/405) |
| `MinimumVolumeBox3FloatingPoint.h` | `GetExtreme` | the strict-improvement hill climb stalls on a floating-point plateau; degenerate non-containing box | RC | fixed | [#426](https://github.com/gradientspaceai/gtengine-js/issues/426) |
| `MinimumVolumeBox3FloatingPoint.h` | misc | duplicated `axis[2]` computation; compute-type literal in `MinimizerConstantS`; stale level-curve comments | minor | corrected | [#352](https://github.com/gradientspaceai/gtengine-js/issues/352) |
| `MinimumVolumeBox3Rational.h` | `MinimizerVariableT` (~L1177) | declares `T const&` parameters where the sibling uses `Number const&`; every exact sample is rounded to double | RC | fixed | [#355](https://github.com/gradientspaceai/gtengine-js/issues/355) |
| `MinimumVolumeBox3Rational.h` | dimension 2 (~L336) | the same Newell wrap-around omission as the floating-point sibling | RC | fixed | [#355](https://github.com/gradientspaceai/gtengine-js/issues/355) |
| `MinimumVolumeBox3Rational.h` | `ComputeVolume` | recomputes `axis[2]` from unchanged operands; `&&indices != nullptr` spacing | minor | corrected | [#355](https://github.com/gradientspaceai/gtengine-js/issues/355) |
| `MinimumWidthPoints2.h` | `ComputeMinWidth` brute-force branch | carries the same duplicate-point corner drop (currently unreachable) | minor | fixed | [#328](https://github.com/gradientspaceai/gtengine-js/issues/328) |
| `MinimumWidthPoints2.h` | frames | inconsistent handedness between the degenerate and general branches | minor | preserved | [#328](https://github.com/gradientspaceai/gtengine-js/issues/328) |
| `NaturalSplineCurve.h` | `CreateClosed` | the wrap-around row uses three plain assignments; for a 3-point closed spline the third overwrites the first | RC | fixed | [#295](https://github.com/gradientspaceai/gtengine-js/issues/295) |
| `NaturalSplineCurve.h` | constructors | `LogAssert` runs after the base constructor has already resized on a negative count and dereferenced `times` | RC | fixed | [#295](https://github.com/gradientspaceai/gtengine-js/issues/295) |
| `NaturalSplineCurve.h` | `CreateFree` | `storageSize` over-allocates by `N-1` reals | minor | preserved | [#295](https://github.com/gradientspaceai/gtengine-js/issues/295) |
| `NearestNeighborQuery.h` | traversal stack | `maxLevel == 32` is allowed but the stack has 32 entries, permitting depth 33 | RC | preserved (harmless in the port: JS arrays grow) | [#48](https://github.com/gradientspaceai/gtengine-js/issues/48) |
| `NearestNeighborQuery.h` | accessors | do not reject negative indices | RC | fixed | [#48](https://github.com/gradientspaceai/gtengine-js/issues/48) |
| `NURBSCircle.h` | `NURBSHalfCircleDegree3` comment | says `x >= 0`; the control points trace `y >= 0` | doc | corrected | [#415](https://github.com/gradientspaceai/gtengine-js/issues/415) |
| `OBBTree.h` | height clamp | `height > 31` is clamped to 31 and then preallocates `2^32 - 1` nodes | RC | preserved | [#274](https://github.com/gradientspaceai/gtengine-js/issues/274) |
| `OBBTreeOfPoints.h`, `OBBTreeOfSegments.h`, `OBBTreeOfTriangles.h` | `ComputeInteriorBox` | seeds `pmin = pmax = Zero()` instead of the first projection, forcing the box to contain the frame origin | minor | preserved | [#103](https://github.com/gradientspaceai/gtengine-js/issues/103), [#274](https://github.com/gradientspaceai/gtengine-js/issues/274) |
| `OBBTreeOfSegments.h` | frame construction | unguarded `Normalize` on a zero-length segment collapses to a zero frame | minor | preserved | [#103](https://github.com/gradientspaceai/gtengine-js/issues/103) |
| `OBBTreeOfTriangles.h` | `Execute` | the same coincident-hit drop as `BVTreeOfTriangles.h` | RC | fixed | [#167](https://github.com/gradientspaceai/gtengine-js/issues/167) |
| `OBBTreeOfTriangles.h` | `Execute`, `IntersectSegmentTriangle` | leaf boxes are never tested against the linear component; unguarded `Length(Q - P)` division | minor | preserved | [#274](https://github.com/gradientspaceai/gtengine-js/issues/274) |
| `OrientedBox.h` | class comment | writes box coordinates with only two terms for an N-dimensional class | doc | corrected | [#78](https://github.com/gradientspaceai/gtengine-js/issues/78) |
| `OrientedBoxTreeOfTriangles.h` | `ComputeLeafBoundingVolume` | the smallest-extent scan's last comparison is `>`, so the largest extent is zeroed and the leaf box collapses | RC | fixed | [#343](https://github.com/gradientspaceai/gtengine-js/issues/343) |
| `Parallelepiped3.h`, `Parallelogram2.h` | `GetVertices` comment | documented counterclockwise; the code emits bit-pattern order | doc | corrected | [#155](https://github.com/gradientspaceai/gtengine-js/issues/155) |
| `Parallelepiped3.h` | handedness assert | `DotCross(...) > 0` accepts a numerically degenerate (repeated-axis) basis | minor | preserved | [#484](https://github.com/gradientspaceai/gtengine-js/issues/484) |
| `ParametricCurve.h` | `GetTime` (~L213) | integrates across knots as a single Romberg call while `GetLength` splits at knots | RC | preserved | [#113](https://github.com/gradientspaceai/gtengine-js/issues/113) |
| `ParametricCurve.h` | lazy init (~L152/L199) | length-zero sentinels re-run the full segment-length initialisation on every call | minor | preserved | [#113](https://github.com/gradientspaceai/gtengine-js/issues/113) |
| `PdeFilter.h`, `PdeFilter1/2/3.h` | `ScaleType::NONE` | still subtracts the data minimum, so "as is" shifts the image and a constant image stores as zeros | RC | preserved | [#60](https://github.com/gradientspaceai/gtengine-js/issues/60) |
| `PdeFilter1/2/3.h` | `OnPreUpdate` | only the Neumann mask border is refreshed, so the image border keeps construction-time values after the first swap | RC | preserved | [#60](https://github.com/gradientspaceai/gtengine-js/issues/60) |
| `PdeFilter.h` | border assignment | `mBorderValue` is written verbatim into offset/scale-transformed buffers | minor | preserved | [#60](https://github.com/gradientspaceai/gtengine-js/issues/60) |
| `PlanarMesh.h` | first constructor | half-constructs on duplicate triangles (zero vertex count, null pointer, unset query) | RC | preserved | [#256](https://github.com/gradientspaceai/gtengine-js/issues/256) |
| `PlanarMesh.h` | `find` results, `Contains` | map lookups dereferenced without an `end()` check; `Contains` lacks the range check its siblings have | RC | fixed (assert) | [#256](https://github.com/gradientspaceai/gtengine-js/issues/256) |
| `Polygon2.h` | constructor | cannot reject an empty vertex pool | minor | fixed (assert) | [#268](https://github.com/gradientspaceai/gtengine-js/issues/268) |
| `Polyhedron3.h` | class comment | says `numIndices` must be 6 or larger; the constructor requires >= 12 | doc | corrected | [#155](https://github.com/gradientspaceai/gtengine-js/issues/155) |
| `Polyhedron3.h` | geometric queries | `auto vertexPool = GetVertices();` deep-copies the whole pool on every call | minor | preserved | [#155](https://github.com/gradientspaceai/gtengine-js/issues/155) |
| `Polynomial1.h` | `SquareFreeFactorization` | the trailing `do/while` is unbounded and the exact-zero GCD test never lets the degree drop: infinite loop | RC | fixed | [#83](https://github.com/gradientspaceai/gtengine-js/issues/83) |
| `Polynomial1.h` | `GreatestCommonDivisor` | the exact-zero remainder test is unusable at double precision; genuine common factors are missed | RC | preserved | [#83](https://github.com/gradientspaceai/gtengine-js/issues/83) |
| `Polynomial1.h` | `GreatestCommonDivisor` | asymmetric for degree-0 inputs: `gcd(c,0) = c` but `gcd(0,c) = 0` | minor | preserved | [#488](https://github.com/gradientspaceai/gtengine-js/issues/488) |
| `PolynomialCurve.h` | constructors | never set `mConstructed`, so `operator bool` reports failure for every valid curve | RC | fixed | [#182](https://github.com/gradientspaceai/gtengine-js/issues/182) |
| `PolynomialCurve.h` | `Evaluate` | `if (order == 3)` where lower orders use `>=`, so `order > 3` leaves `jet[3]` untouched | minor | preserved | [#182](https://github.com/gradientspaceai/gtengine-js/issues/182) |
| `PolynomialRoot.h` | `PolynomialRootBisect` | the sign-mismatch branches collapse the interval, so an exact endpoint root is never reported as bracketed | minor | preserved | [#319](https://github.com/gradientspaceai/gtengine-js/issues/319) |
| `PolynomialRoot.h` | `operator==`, comment | ignores the multiplicity `m`; typo `// x is the root estimate)and m` | minor | preserved | [#319](https://github.com/gradientspaceai/gtengine-js/issues/319) |
| `PolylineOffset.h` | member initializer list (~L61-73) | `vertices.size() - 1` wraps for an empty open polyline, so allocation throws before the intended `LogAssert` | WR | fixed | [#112](https://github.com/gradientspaceai/gtengine-js/issues/112) |
| `PrimalQuery2.h` | four-argument `ToLine` | squares `P - V0` instead of `V1 - V0`, swapping collinear orders 0 and +2 | RC | preserved | [#100](https://github.com/gradientspaceai/gtengine-js/issues/100) |
| `PrimalQuery2.h`, `PrimalQuery3.h` | precision tables | the BSRational N-value tables were computed with a removed BSPrecision API and no longer match | doc | preserved | [#43](https://github.com/gradientspaceai/gtengine-js/issues/43) |
| `PrimalQuery2.h` | `order` table, `ToCircumcircle` comment | mislabels +2, omits +3; "involves three calls of ToLine" is false | doc | corrected | [#101](https://github.com/gradientspaceai/gtengine-js/issues/101) |
| `Projection.h` | `PerspectiveProject` | discards the `bool` from `Ellipse2::FromCoefficients` | minor | fixed | [#225](https://github.com/gradientspaceai/gtengine-js/issues/225) |
| `QuadricSurface.h` | `GetClassification` (~L200) | symmetrizes `A` from the upper triangle while `F`/`FX`/`FY`/`FZ` use the stored matrix verbatim | RC | preserved | [#318](https://github.com/gradientspaceai/gtengine-js/issues/318) |
| `QuadricSurface.h` | ~L210, ~L195 | dead `rS12`; the "Sturm sequences" comment describes Descartes' rule of signs | minor | corrected | [#318](https://github.com/gradientspaceai/gtengine-js/issues/318) |
| `Quaternion.h` | `Slerp` comment | the two `sin` coefficients are swapped, so the formula yields `q1` at `t = 0` | doc | corrected | [#160](https://github.com/gradientspaceai/gtengine-js/issues/160) |
| `RectangleMesh.h` | `InitializeFrame` | hardcodes `tangent = (1,0,0)`, `bitangent = (0,1,0)` regardless of the rectangle's axes | RC | fixed | [#268](https://github.com/gradientspaceai/gtengine-js/issues/268) |
| `RemezAlgorithm.h` | `ComputeXExtremes` | pins the outer nodes to `xMin`/`xMax`, so a non-minimax fixed point is reported as success | RC | preserved | [#478](https://github.com/gradientspaceai/gtengine-js/issues/478) |
| `RemezAlgorithm.h` | `ComputePCoefficients` | reads `poly[i]` past `mCoefficient` after leading-zero elimination | RC | fixed | [#147](https://github.com/gradientspaceai/gtengine-js/issues/147) |
| `RemezAlgorithm.h` | `maxBracketIterations` | validated and stored but never read; the QFLS comments describe a search that is now plain bisection | minor | preserved | [#147](https://github.com/gradientspaceai/gtengine-js/issues/147) |
| `RemezAlgorithm.h` | `Execute` | leaves `GetXNodes()` one exchange ahead of `GetErrors()` and `GetCoefficients()` | minor | preserved | [#478](https://github.com/gradientspaceai/gtengine-js/issues/478) |
| `ReparameterizeByArclength.h` | `DoNewtonsMethod` | falls through to `tMid - fMid/dfdt` after the `dfdt == 0` bisection step, producing an infinite quotient | minor | preserved | [#183](https://github.com/gradientspaceai/gtengine-js/issues/183) |
| `RevolutionMesh.h` | DISK texture coordinates | the declared `origin{0.5, 0.5}` is never added to the ring, putting the hub at a corner | RC | preserved | [#412](https://github.com/gradientspaceai/gtengine-js/issues/412) |
| `RevolutionMesh.h` | one-row SPHERE | divides by `numRows - 1` although `MeshDescription` clamps SPHERE to `numRows >= 1`; NaN tcoords | RC | preserved | [#412](https://github.com/gradientspaceai/gtengine-js/issues/412) |
| `RiemannianGeodesic.h` | `ComputeMetricDerivative` | computes `2*Gamma_{d,i0i1}`, not `dg_{i0i1}/dx_d` | minor | preserved | [#295](https://github.com/gradientspaceai/gtengine-js/issues/295) |
| `RiemannianGeodesic.h` | derived parameters | computed only in the constructor although the inputs are public and documented as tweakable | minor | fixed (accessor) | [#295](https://github.com/gradientspaceai/gtengine-js/issues/295) |
| `RiemannianGeodesic.h` | `ComputeGeodesic` | `LogAssert(subdivisions < 32)` admits 31, where `1 << subdivisions` is signed overflow | minor | preserved | [#295](https://github.com/gradientspaceai/gtengine-js/issues/295) |
| `RigidBody.h` | `RigidBodyState::SetQOrientation` | builds `mROrientation` from the raw parameter, not the normalized member; wrong on every RK4 stage | RC | fixed | [#313](https://github.com/gradientspaceai/gtengine-js/issues/313) |
| `RigidBody.h` | `SetBodyInertia`, `GetAngularVelocity` | inconsistent `IsMovable()` guard; stray `;;` | minor | preserved | [#313](https://github.com/gradientspaceai/gtengine-js/issues/313) |
| `RootsBisection1.h` | `operator()` | with `maxIterations == 1` the loop never runs and `tRoot`/`fAtTRoot` are returned unassigned | RC | fixed | [#84](https://github.com/gradientspaceai/gtengine-js/issues/84) |
| `RootsBisection2.h` | inner bisector | with `maxIterations == 1` returns 2 without writing its outputs, so members keep a previous call's values | RC | deviates | [#84](https://github.com/gradientspaceai/gtengine-js/issues/84) |
| `RootsBisection2.h` | `XFunction` | `mNoGuaranteeForRootBound = ...` assigns instead of accumulating, erasing earlier y-failures | minor | preserved | [#152](https://github.com/gradientspaceai/gtengine-js/issues/152) |
| `RootsBisection2.h` | `operator()` | the y-outputs hold the last `XFunction` call's values and need not match the returned x-root | RC | preserved | [#152](https://github.com/gradientspaceai/gtengine-js/issues/152) |
| `RootsCubic.h` | `signDelta < 0` branch | bisects over `[-b, b]` with `b = max(1, abs(d0), abs(d1))`, which is not a valid root bound | RC | fixed | [#340](https://github.com/gradientspaceai/gtengine-js/issues/340) |
| `RootsCubic.h` | closed form, comment | the `sin(theta/3) <= 0` reorder is unreachable; the second Samuelson subinterval is mislabelled | minor | corrected | [#340](https://github.com/gradientspaceai/gtengine-js/issues/340) |
| `RootsGeneralPolynomial.h` | `operator()` | the rational polynomial is allocated at the untrimmed size, so the solver runs on a zero-padded, non-monic polynomial | minor | fixed | [#319](https://github.com/gradientspaceai/gtengine-js/issues/319) |
| `RootsGeneralPolynomial.h` | `Bisect` | dead stores `tPMax`/`tPMin = tPAtRoot` | minor | dropped | [#319](https://github.com/gradientspaceai/gtengine-js/issues/319) |
| `RootsQuartic.h` | four-real-root classification | `signDelta > 0 && rD2 > 0` is only half the criterion; two complex pairs are reported as four real roots | RC | fixed | [#340](https://github.com/gradientspaceai/gtengine-js/issues/340) |
| `RootsQuartic.h` | square-root extraction | a shared two-element array is reused and index 1 is read unconditionally, returning a stale root | RC | fixed | [#340](https://github.com/gradientspaceai/gtengine-js/issues/340) |
| `RotatingCalipers.h` | `CreatePolygon` | collinearity is tested against the immediately preceding edge, so a duplicate point discards the next real corner | RC | fixed | [#286](https://github.com/gradientspaceai/gtengine-js/issues/286) |
| `RotatingCalipers.h` | `CreatePolygon` | reads `vertices.back()` and `vertices[1]` before the size assert | RC | fixed (assert) | [#286](https://github.com/gradientspaceai/gtengine-js/issues/286) |
| `RotatingCalipers.h` | `ComputeAntipodes` | duplicates one edge and omits another when caliper angles tie | minor | preserved | [#286](https://github.com/gradientspaceai/gtengine-js/issues/286) |
| `Rotation.h` | `operator()(i0,i1,i2)` | the `IS_EULER_ANGLES` arm is a bare `break` after the axes have been overwritten; returns relabelled angles | RC | fixed | [#225](https://github.com/gradientspaceai/gtengine-js/issues/225) |
| `Rotation.h` | `Convert(Matrix, AxisAngle)` | near angle pi the antisymmetric part underflows, giving a zero or non-unit axis | RC | fixed | [#374](https://github.com/gradientspaceai/gtengine-js/issues/374) |
| `Rotation.h` | `Convert(Matrix, AxisAngle)` | slightly further from pi, absorbed round-off gives an axis wrong by up to 5e-3 | RC | preserved | [#374](https://github.com/gradientspaceai/gtengine-js/issues/374) |
| `Rotation.h` | `Convert(Matrix, EulerAngles)` | the gimbal-lock test uses exact `r(i,j) == +-1` | RC | preserved | [#374](https://github.com/gradientspaceai/gtengine-js/issues/374) |
| `RotationEstimate.h` | `C_ROTC*_EST_MAX_ERROR` | the published maximum-error tables understate the true error at high degrees and for all of `rotc4` | doc | preserved | [#225](https://github.com/gradientspaceai/gtengine-js/issues/225) |
| `RotationEstimate.h` | `C_ROTC4_EST_COEFF` degree 14 | the last coefficient's exponent looks like a typo (degree 14 is worse than degree 12) | minor | preserved | [#225](https://github.com/gradientspaceai/gtengine-js/issues/225) |
| `SampleCircularArc.h` | `SampleArc3`, `SampleArc4` | the trisector/quadsector directions are exact only for the 1:1 midpoint construction | RC | preserved | [#183](https://github.com/gradientspaceai/gtengine-js/issues/183) |
| `SegmentMesh.h` | `CONTIGUOUS_CLOSED` | the constructor stores `S[i] = {i-1, i}`, a rotation of the documented `{i, (i+1) % L}` | doc | preserved | [#78](https://github.com/gradientspaceai/gtengine-js/issues/78) |
| `SegmentMesh.h` | `DISJOINT` | the doc requires an even vertex count; the constructor silently drops the last vertex for odd counts | minor | preserved | [#78](https://github.com/gradientspaceai/gtengine-js/issues/78) |
| `SeparatePoints2.h` | `OnSameSide`, `WhichSide` | the side test uses a rounded normal and constant, so a candidate edge's own endpoint classifies positive | RC | fixed | [#328](https://github.com/gradientspaceai/gtengine-js/issues/328) |
| `SeparatePoints2.h` | hull scan | each hull vertex is classified twice | minor | preserved | [#328](https://github.com/gradientspaceai/gtengine-js/issues/328) |
| `SeparatePoints3.h` | `OnSameSide`, `WhichSide` | the same rounding defect; 101 of 400 overlapping cloud pairs were reported separated | RC | fixed | [#348](https://github.com/gradientspaceai/gtengine-js/issues/348) |
| `SeparatePoints3.h` | cross-product-axis loop | sets `normal`/`constant` while `origin` keeps the last face plane, the inconsistent state `Hyperplane.h` warns about | RC | fixed | [#348](https://github.com/gradientspaceai/gtengine-js/issues/348) |
| `SingularValueDecomposition.h` | singular-value sort | `std::sort` is unstable, so the U/V permutation for tied singular values is unspecified | minor | fixed | [#478](https://github.com/gradientspaceai/gtengine-js/issues/478) |
| `Slerp.h` | header formula | second numerator term is `sin(theta)*q1`; it must be `sin(t*theta)*q1` | doc | corrected | [#57](https://github.com/gradientspaceai/gtengine-js/issues/57), [#427](https://github.com/gradientspaceai/gtengine-js/issues/427) |
| `SortPointsOnCircle.h` | `LessThanByGeometry` | not a strict weak ordering when a point coincides with the sort centre | minor | preserved | [#394](https://github.com/gradientspaceai/gtengine-js/issues/394) |
| `SplitMeshByPlane.h` | `SplitTriangle*` | read the edge map with `std::map::operator[]`, default-constructing index 0 for a missing key | RC | fixed (assert) | [#310](https://github.com/gradientspaceai/gtengine-js/issues/310) |
| `SplitMeshByPlane.h` | `ClassifyTriangles`, `mEMap` | narrowing cast; the `Vector3` half of the map value is never read | minor | preserved | [#310](https://github.com/gradientspaceai/gtengine-js/issues/310) |
| `SqrtEstimate.h` | `SqrtEstimateRR(0)` | returns a meaningless nonzero value despite the documented `x >= 0` domain | minor | preserved | [#57](https://github.com/gradientspaceai/gtengine-js/issues/57) |
| `StaticVETManifoldMesh2.h` | `GetAdjacentTriangles` | the reversed-direction fallback returns R then L without swapping; documented case 3 unreachable | RC | fixed | [#66](https://github.com/gradientspaceai/gtengine-js/issues/66), [#472](https://github.com/gradientspaceai/gtengine-js/issues/472) |
| `StaticVTSManifoldMesh3.h` | `GetAdjacentTetrahedra` | the same missing swap after `SortFace` | RC | fixed | [#66](https://github.com/gradientspaceai/gtengine-js/issues/66), [#472](https://github.com/gradientspaceai/gtengine-js/issues/472) |
| `SurfaceExtractor.h` | `MakeUnique` | keeps rotated duplicate triangles and emits non-canonical rotations | minor | preserved | [#439](https://github.com/gradientspaceai/gtengine-js/issues/439) |
| `SurfaceExtractorMC.h` | `Extract` | the edge interpolation omits `level`; vertices are misplaced for any nonzero level | RC | fixed | [#443](https://github.com/gradientspaceai/gtengine-js/issues/443) |
| `SurfaceExtractorMC.h` | 15-case table | not face-consistent: a shared face with alternating signs is resolved differently by each voxel | RC | preserved | [#443](https://github.com/gradientspaceai/gtengine-js/issues/443) |
| `SurfaceExtractorTetrahedra.h` | `GetGradient` (~L296) | the odd-parity branch condition `dx + dy + dz >= 0` is unconditionally true; the plane is at 2 | RC | fixed | [#132](https://github.com/gradientspaceai/gtengine-js/issues/132) |
| `SWInterval.h` | outward rounding | `std::nextafter(value, +-max)` pulls an infinite bound back to `+-MAX_VALUE`, destroying enclosure | RC | preserved | [#50](https://github.com/gradientspaceai/gtengine-js/issues/50), [#367](https://github.com/gradientspaceai/gtengine-js/issues/367) |
| `SymmetricEigensolver.h` | `Tridiagonalize` | stores the reflection parameter for a degenerate Householder step, corrupting `Q` | RC | fixed | [#80](https://github.com/gradientspaceai/gtengine-js/issues/80) |
| `SymmetricEigensolver3x3.h` | `GetCosSin` | lacks the `maxAbsComp` rescaling its 2x2 sibling documents; underflows for extreme matrix scales | RC | fixed | [#379](https://github.com/gradientspaceai/gtengine-js/issues/379) |
| `TanEstimate.h` | `TanEstimateRR` comment | claims `r` in `[-pi, pi]`; `remainder(x, pi)` gives `[-pi/2, pi/2]`, leaving two branches dead | doc | preserved | [#57](https://github.com/gradientspaceai/gtengine-js/issues/57) |
| `TCBSplineCurve.h` | constructor | never sets `mConstructed`, so `operator bool` reports failure for every valid curve | RC | fixed | [#182](https://github.com/gradientspaceai/gtengine-js/issues/182) |
| `TCBSplineCurve.h` | `ComputeInteriorTangents` | the lambda pass divides by zero when both Kochanek-Bartels tangents vanish; NaN coefficients | RC | preserved | [#415](https://github.com/gradientspaceai/gtengine-js/issues/415) |
| `Tetrahedron3.h` | `GetPlanes` | sets only `normal` and `constant`, leaving `Plane3::origin` at the world origin | RC | fixed | [#268](https://github.com/gradientspaceai/gtengine-js/issues/268) |
| `TetrahedraRasterizer.h` | `ClipCullAABBs` | clips the stored boxes in place, so a second call with a larger region scans the clipped box | RC | preserved | [#439](https://github.com/gradientspaceai/gtengine-js/issues/439) |
| `Torus3.h` | `GetParameters` | documents `u, v` in `[0, 2*pi)` but returns `atan2` values in `[-pi, pi]` | doc | preserved | [#455](https://github.com/gradientspaceai/gtengine-js/issues/455) |
| `Torus3.h` | `Evaluate` | guards the second-order derivatives with `maxOrder == 2` rather than `>= 2` | minor | preserved | [#484](https://github.com/gradientspaceai/gtengine-js/issues/484) |
| `Transform.h` | `GetHInverse` | the RS branches assume the last row is still `(0,0,0,1)`; after a singular general state it is `(0,0,0,0)` | RC | fixed | [#265](https://github.com/gradientspaceai/gtengine-js/issues/265) |
| `Transform.h` | `Inverse()` | stores the full affine matrix in the M channel, violating `GetMatrix`'s documented block structure | RC | fixed | [#265](https://github.com/gradientspaceai/gtengine-js/issues/265) |
| `Transform.h` | misc | `SetRotation(AxisAngle<3>)` lifts with `w = 1`; dead `Invert3x3`; `mIsUniformScale` never reset; inverted doc | minor | preserved | [#265](https://github.com/gradientspaceai/gtengine-js/issues/265) |
| `Triangle.h` | comment | "sets the/ vertices" | doc | corrected | [#78](https://github.com/gradientspaceai/gtengine-js/issues/78) |
| `TriangulateCDT.h` | `RemapPolygonTree` | a later duplicate overwrites the first occurrence's remapping, so a never-passed index is reported | minor | preserved | [#348](https://github.com/gradientspaceai/gtengine-js/issues/348), [#405](https://github.com/gradientspaceai/gtengine-js/issues/405) |
| `TriangulateCDT.h` | `ConstrainedTriangulate` | re-inserts triangles the graph copy already contains | minor | preserved | [#348](https://github.com/gradientspaceai/gtengine-js/issues/348) |
| `TriangulateEC.h` | `VertexList::RemoveR` | removing the only reflex vertex leaves `mRLast` set, breaking the class's own emptiness invariant | minor | preserved | [#175](https://github.com/gradientspaceai/gtengine-js/issues/175) |
| `TriangulateEC.h` | `DoEarClipping` | closes the ear ring without checking `mEFirst != -1` | minor | preserved | [#175](https://github.com/gradientspaceai/gtengine-js/issues/175) |
| `TubeMesh.h` | `UpdatePositions` | the closed-tube fixup uses the `numCols` stride and stops one vertex short, cracking the seam | RC | fixed | [#240](https://github.com/gradientspaceai/gtengine-js/issues/240) |
| `TubeMesh.h` | ring sampling | rings have only `numCols - 1` distinct angles; a closed tube drops the last ring | minor | preserved | [#240](https://github.com/gradientspaceai/gtengine-js/issues/240) |
| `UIntegerALU32.h` | `RoundUp` comment | the return value is the trailing-zero shift count, not the 0/1 carry the comments imply | doc | corrected | [#96](https://github.com/gradientspaceai/gtengine-js/issues/96) |
| `UnsymmetricEigenvalues.h` | eigenvalue packing | the loop bound drops `A(N-1,N-1)` when the last row decouples as a 1x1 block | RC | fixed | [#42](https://github.com/gradientspaceai/gtengine-js/issues/42) |
| `UnsymmetricEigenvalues.h` | `FrancisQRStep` | no exceptional shift, so the iteration can cycle on a well-separated real spectrum | RC | preserved | [#476](https://github.com/gradientspaceai/gtengine-js/issues/476) |
| `UnsymmetricEigenvalues.h` | `Solve` doc | claims a `0xFFFFFFFF` non-convergence return that does not exist | doc | corrected | [#476](https://github.com/gradientspaceai/gtengine-js/issues/476) |
| `Vector4.h` | `ComputeOrthogonalComplement` | the `maxIndex == 3` branch yields the zero vector whenever components 1 and 2 vanish | WR | preserved | [#87](https://github.com/gradientspaceai/gtengine-js/issues/87) |
| `VEManifoldMesh.h` | `Insert` | writes the new edge into `mEMap` before the nonmanifold check, leaving a phantom edge on failure | RC | preserved | [#73](https://github.com/gradientspaceai/gtengine-js/issues/73) |
| `VertexCollapseMesh.h` | `Collapsed` comment, `VCM_NO_MORE_ALLOWED` | the comment describes a restore that never happens; the code never returns that status | doc | corrected | [#295](https://github.com/gradientspaceai/gtengine-js/issues/295), [#412](https://github.com/gradientspaceai/gtengine-js/issues/412) |
| `VertexCollapseMesh.h` | `TriangulateLink`, `Collapsed` | floating-point ear clipping of a link with collinear vertices returns a duplicate triangle; the old fan is removed before the failure is noticed, so `DoCollapse` returns false with a corrupted mesh | RC | fixed | [#498](https://github.com/gradientspaceai/gtengine-js/issues/498) |
| `VertexCollapseMesh.h` | `DoCollapse` | `record.vertex = 0x80000000` relies on implementation-defined conversion pre-C++20 | minor | preserved | [#295](https://github.com/gradientspaceai/gtengine-js/issues/295) |
| `VETNonmanifoldMesh.h` | `Remove` | the assertion is inverted and fires for every well-formed mesh | RC | fixed | [#240](https://github.com/gradientspaceai/gtengine-js/issues/240) |
| `VTSManifoldMesh.h` | `Remove` | over-erases `VAdjacent`, dropping adjacencies still contributed by surviving faces | RC | fixed | [#256](https://github.com/gradientspaceai/gtengine-js/issues/256) |

## Findings by upstream header

Headers are in alphabetical order. Where one defect spans a family of files the
section is filed under the alphabetically first member and names the others.
Purely cosmetic comment items are listed compactly; anything result-corrupting
carries its cause and a concrete reproduction.

### `AdaptiveSkeletonClimbing2.h`, `FastMarch.h`

- `LinearMergeTree::GetEdge` can never return -1, so the `!= -1` guards in
  `GetRectangle` are dead code.
- The constructor comment says `N >= 0` is accepted; the code rejects `N <= 0`.
- `FastMarch::GetTimeExtremes` redundantly re-tests one element.

All minor, preserved. Issue [#52](https://github.com/gradientspaceai/gtengine-js/issues/52).

### `AdaptiveSkeletonClimbing3.h`

**Symptom.** A whole-image monobox is silently dropped and the extractor returns
an empty mesh.

**Cause.** `Merge(v, LX, LY, LZ, x0, y0, z0, stride, depth)` adds a fully merged
box on behalf of its *parent* and signals "I merged completely" by returning
`true`. The root invocation `Merge(0, 0, 0, 0, 0, 0, 0, mTwoPowerN, depth)`
discards that return value, so when the entire image merges into one box nobody
calls `AddBox` for it.

**Reproduction.** The ramp `F(x,y,z) = x` on a 9^3 image (`N = 8`), extracted at
level 3.5 with `depth <= 0`, yields 0 boxes and 0 triangles although the planar
surface `x = 3.5` is present. `AdaptiveSkeletonClimbing2` avoids this because
`GetRectangles` is a separate pass starting from the root.

Also in this header (both minor): the guard
`merge->IsZeroEdge(L) || merge->HasZeroSubedge(L)` admits a `CFG_MULT` node whose
`mZeroBases[L]` is -1, which would index the image at -1 along that axis
(extensive probing never reached it); and `GetVertices` narrows `Real` box
corners through `static_cast<float>`.

Issue [#194](https://github.com/gradientspaceai/gtengine-js/issues/194). Port: preserved.

### `APConversion.h`

**1. `EstimateAmB` uses stale bisection squares (result-corrupting).**

**Cause.** In the bisection loop's rounding branch `tMin`/`tMax` are updated but
`tMinSqr`/`tMaxSqr` are not recomputed. When the loop exhausts its iteration
budget the subsequent Newton bound is computed from a square that does not
correspond to its `t`.

**Reproduction.** Using continued-fraction convergents of `(7 + 3*sqrt(5))/2`,
where `a^2 - 3ab + b^2` is within an ulp of zero, 145 of 304 (convergent,
maxIterations) pairs return a `[tMin, tMax]` that does not bracket `a - b`; for
example `a^2 = 10983760033`, `b^2 = 1602508992`. All 145 become valid brackets
once the squares are recomputed on update, which is provably a no-op on every
internally consistent path.

**Suggested fix.** Recompute `tMinSqr`/`tMaxSqr` wherever `tMin`/`tMax` are
assigned in that branch. Port: fixed.

**2-4 (minor/doc, preserved).** `EstimateAmB` has an undocumented
`aSqr >= bSqr` precondition (otherwise the zero-clamp produces the inverted
bracket `[max, min]`); `EstimateSqrt(0)` returns `[0, 5e-324]`, violating the
documented strict lower bound `aMin < a`, and there is no negative-input guard;
the comment `r^2 - 7*r^2 + 1 = 0` should read `r^2 - 7*r + 1 = 0`.

Issue [#280](https://github.com/gradientspaceai/gtengine-js/issues/280).

### `APInterval.h`

**Symptom.** The `SetSign(+-2)` infinity sentinels order incorrectly against
finite values.

**Cause.** They encode zero magnitude with a nonzero sign, which violates the
`BSNumber`/`BSRational` invariant: `IsValid()` rejects them, `posinf < 1` is
true, and converting a sentinel produces `0.333...`.

Also: `operator/` with divisor `[0, 0]` reaches `1/0` (an assert inside the
rational reciprocal) instead of returning the documented whole-line `Reals()`
that the zero-containing-divisor branch produces for `[0, x]`. `QFN2` is a dead
alias.

Issue [#280](https://github.com/gradientspaceai/gtengine-js/issues/280). Port: preserved (fixing would change `BSNumber` semantics); an `isInfinite()` predicate was added so callers can detect indeterminate results.

### `ApprCircle2.h`, `ApprSphere3.h`

The `FitUsingLengths` doc comment describing `initialCenterIsAverage` is
inverted relative to the code (`if (initialCenterIsAverage) { center = average; }`).
Documentation only. Issue [#92](https://github.com/gradientspaceai/gtengine-js/issues/92).

### `ApprCone3.h`

`ComputeInitialCone` divides by `rRange/hRange` and `rMax/tanAngle` with no
guard (all points at one height, or a right-angle estimate), and discards
`ApprHeightLine2::Fit`'s return value, so a failed height-line fit silently
produces a zeroed initial cone that is handed to the minimizers. Both
preserved. Issue [#271](https://github.com/gradientspaceai/gtengine-js/issues/271).

### `ApprCone3EllipseAndPoints.h`

**1. `ApprCone3ExtractEllipses::AssociatePointsWithPlanes` indexes with `size_t(-1)`.**
When `mPlanes` is empty the minimising index stays at its sentinel and
`mIndices[minJ]` is an out-of-bounds access on the first point. Reachable with
fewer than three input points, because no tree node then satisfies
`maxIndex >= minIndex + 2` and no plane is ever located. Port: returns early.

**2. `ComputeEllipse` runs on planes with an empty index set (result-corrupting).**
Any tree node holding exactly three points has a flat box (three points are
always coplanar), so `LocatePlanes` readily emits planes that later receive no
points; upstream then runs `ApprGaussian3::FitIndexed` with zero indices and
`ApprEllipse2` over an empty vector. Reproduction: three exact circular cone
sections produced 3 correct planes plus 2 spurious empty ones. Port: discards
point-less planes, keeping the plane, index and ellipse arrays parallel.

**3. Planes with only 1-2 supporting points survive (result-corrupting, preserved).**
Reproduction: two exact circular cone sections of 19 points each at `z = 20` and
`z = 40` with `boxExtentEpsilon = 1e-6` give 3 planes: the two real ones (19 and
18 points) plus a spurious tilted plane `n = (0.511268, -0.782555, 0.355264)`
holding 1 point, whose ellipse has `extent = (Infinity, Infinity)`. Randomized
runs also show the opposite failure, a real section missed. Damage is contained:
such an ellipse makes the error function NaN and both implementations then throw
"Failed to find fitted cone."

**4-7 (fitter side).** `Fit` divides by `points.size()` with no empty-set guard
(NaN error function, garbage cone); `ComputeCone` computes `b / a` without
validating the fitted ellipse extents (Inf/NaN on a degenerate ellipse); both
guarded in the port. `Control` duplicates default member initializers for three
of six fields. `ProcessPlane` compares plane constants (a distance) against
`mCosAngleEpsilon` (a cosine), so one epsilon serves two incompatible units.
The error metric is `sqrt(error)/n` rather than `sqrt(error/n)` (argmin
unaffected).

Issue [#349](https://github.com/gradientspaceai/gtengine-js/issues/349).

### `ApprCurveByArcs.h`

`if (fabs(det) >= epsilon)` with `epsilon` defaulting to 0 is always true, so an
exactly colinear triple divides by zero instead of taking the documented
`numeric_limits::max()` sentinel path. Separately, with `numArcs == 1` on a
closed curve the chord the bisection needs collapses; combined with the
always-true determinant test the arc comes back NaN or meaningless (radius 1.414
for the unit circle, centred at (3,-1)). A stale comment writes
`{t0,t1} = {times[2*i], times[2*i+1]}` where the code uses `times[2*i+2]`.

Issue [#163](https://github.com/gradientspaceai/gtengine-js/issues/163). Port: preserved.

### `ApprCylinder3.h`

Minor and documentation only, all preserved: `FinishCylinder` initialises `hmax`
to 0 instead of `-max()` while `hmin` gets `+max()` (harmless because the points
are mean-centred, so the projections straddle 0); the multithreaded partition
starts at `j = 0`, re-evaluating the north pole `numThetaSamples` times, and
`numThreads > numPhiSamples` leaves all but the last worker idle; the mesh
`operator()` docs say "cone axis/height/center" (stale copy-paste from
`ApprCone3.h`), and the class doc calls the return an RMS error while `G(...)`
returns the mean squared error.

Issue [#224](https://github.com/gradientspaceai/gtengine-js/issues/224).

### `ApprEllipse2.h`

The `bool` returned by `GetContainer` is ignored, so on failure `box` is read
uninitialised (undefined behaviour in C++; the port throws). Documentation:
the `H(t)` paragraph says "smallest G-value" where it means H-value and
`M - T * dF/dC` where it means `dF/dM`; `UpdateMatrix` refers to a renamed
`negGradM`. Dead store: the initial `ErrorFunction` result is overwritten
whenever `numIterations > 0`.

Issue [#224](https://github.com/gradientspaceai/gtengine-js/issues/224).

### `ApprEllipseByArcs.h`

**Symptom.** A returned arc silently repeats the previous arc's centre and
radius, and the function still reports success.

**Cause.** The intermediate-arc loop ignores `Circumscribe`'s return value. On a
collinear point triple `Circumscribe` returns `false` and leaves `circle`
untouched, so the loop stores whatever was there from the previous iteration.
The two end-arc computations do check the flag, which is the evidence this is an
oversight.

**Suggested fix.** Propagate the failure as the end-arc code already does.
Port: fixed. Documentation in the same header calls the outputs "returned input
... arrays" and writes `numArc` for `numArcs`.

Issue [#322](https://github.com/gradientspaceai/gtengine-js/issues/322).

### `ApprEllipsoid3.h`

`GetContainer`'s return value is ignored for the initial oriented box, so the
default-constructed box is read when `ApprGaussian3::Fit` fails (same pattern as
`ApprEllipse2.h`). Header comment errors: "smallest G-value" should be H-value,
`M' = M - T * dF/dC` should be `dF/dM`, typos "definitess" and "ellisoid".

Issue [#322](https://github.com/gradientspaceai/gtengine-js/issues/322).

### `ApprParabola2.h`, `ApprParaboloid3.h`

All four `Fit`/`FitRobust` variants compute
`*meanSquareError = std::sqrt(totalSqrError) / numPoints`, which is `sqrt(sum)/n`:
neither the mean squared error `sum/n` nor the RMS error `sqrt(sum/n)`. The
value is monotone in the residual sum, so the fitted coefficients are
unaffected, but it is not the documented quantity. `ApprParaboloid3.h` also
carries a stale comment, "the i-index into A(2,2) is 35", which should reference
`A(5,5)`.

Issue [#271](https://github.com/gradientspaceai/gtengine-js/issues/271). Port: preserved.

### `ApprParallelLines2.h`

**1. `ComputeF`: `a30[1] = -mR3;` is missing the factor `Z12` (result-corrupting).**

**Cause.** It should be `-mR3 * data.Z12`, per the equivalent evaluation in
`UpdateParameters` (`A30 = -sigma*(3*Z12 + (Z30 - 3*Z12)*sigma^2)`) and the
expansion of `S30 = average((-sigma*x + gamma*y)^3)`.

**Reproduction.** Over 200 random 12-point sets the fit differs from the
corrected version in all 200 trials; the worst direction difference is 1.88, and
one trial returns the un-improved seed.

**2. `Fit` does not reject roots with `sigma^2 > 1`.** `(gamma, sigma)` is unit
length, so those are spurious roots of the squared polynomial `h`, and their
"error" can go negative. On 14 samples on `y = +-3` upstream returns direction
`(0, 1.7320509)`, which is not unit, and radius NaN. Port filters
`sigmaSqr > 0 && sigmaSqr <= 1`.

**3. `Fit`: `gamma = std::sqrt(sigma)` in the `f1 == 0` branch** violates
`gamma^2 + sigma^2 = 1`; it should be `sqrt(1 - sigmaSqr)`. Hard to reach,
because the moment conditions that kill `f1` also kill `f0`.

**4. `Fit` reads `Polynomial1::operator[]` out of range.** `f0[2*i]`,
`f1[2*i+1]` and `h[2*i]` assume degrees 8, 7 and 16, but the arithmetic
operators eliminate leading zeros, so cancellation makes these index past the
end of the `std::vector`. Port reads missing coefficients as 0.

**5. The returned direction is never renormalised (result-corrupting, preserved).**
`gamma` comes from the polynomial quotient, so `Dot(V, V)` was observed at 1.0019
on ordinary data and 1.8e169 on near-degenerate data. Distinct cause from item 2.

Issues [#91](https://github.com/gradientspaceai/gtengine-js/issues/91) (items 1-4, all fixed in the port) and [#380](https://github.com/gradientspaceai/gtengine-js/issues/380) (item 5).

### `ApprPolynomialSpecial2/3/4.h`

`ApprPolynomialSpecial3/4`'s constructor rejects term sets the header promises
to accept. The header documents "distinct pairs `<p[i], q[i]>`", but the
constructor asserts that each degree list is *separately* strictly increasing.
That rejects the affine model `{1, x, y}` (`xDegrees = {0,1,0}`,
`yDegrees = {0,0,1}`); only "diagonal" term sets such as `{1, x*y, x^2*y^2}`
are admissible. It looks like the correct 1-D check from `ApprPolynomialSpecial2`
was copied per axis instead of being replaced by a tuple-distinctness test.

Separately, `Transform` divides by the sample range with no guard, so constant
height data (`max == min`) yields a NaN model with no failure indication.

Issues [#163](https://github.com/gradientspaceai/gtengine-js/issues/163), [#380](https://github.com/gradientspaceai/gtengine-js/issues/380). Port: preserved.

### `ApprQuadratic2.h`, `ApprQuadratic3.h`, `ApprQuadraticCircle2.h`, `ApprQuadraticSphere3.h`

- Comment says `M = (sum_i V[i])(sum_i V[i])^T` (rank 1); the code accumulates
  `sum_i V[i] V[i]^T`. Documentation only.
- `ApprQuadratic3.h` is missing the trailing `M(0,0) = 1` its three siblings
  have (harmless, the division already leaves exactly 1).
- `ApprQuadraticCircle2` and `ApprQuadraticSphere3` divide by the last
  eigenvector component with no guard. When it is zero the centre becomes
  `(-Infinity, NaN)` and the radius NaN, while the returned measure is `0`, that
  is, reported as a perfect fit. Result-corrupting, preserved.

Issues [#163](https://github.com/gradientspaceai/gtengine-js/issues/163), [#380](https://github.com/gradientspaceai/gtengine-js/issues/380).

### `ApprQuery.h`

**1. `RANSAC` re-seeds the random engine inside the iteration loop.**
`std::default_random_engine` is constructed per iteration with the same seed, so
`std::shuffle` produces the same permutation every time. When an iteration swaps
nothing, every subsequent iteration is bit-identical: `numIterations` has no
observable effect (results identical for 1, 2, 10 and 50 iterations).

**2. `RANSAC` can return `true` without ever fitting `bestModel`.** When
`numRequiredForGoodFit <= GetMinimumRequired()` the consensus branch that
assigns `bestModel` is skipped, yet the function still reports success, leaving
the model default-initialised.

**3. Minor.** `ApprGaussian3::GetMinimumRequired()` returns 2, so RANSAC seeds a
3D Gaussian from a rank-1 covariance.

Issue [#378](https://github.com/gradientspaceai/gtengine-js/issues/378). Port: preserved.

### `ApprTorus3.h`

`operator()` divides by `b0` and passes `f3` as a `SolveCubic` leading
coefficient, both unguarded; on degenerate input the resulting NaNs fail the
`v > 0` / `u > v` validity tests, so the caller sees failure rather than a
garbage torus. The header comment defines `a2` twice (the third occurrence
should be `a0`). The `Rational f0..f3` locals are dead code: the right-hand
sides are computed in `Real`, `Real` to `BSRational` is exact, and
`SolveCubic(Real const&)` converts straight back by round-to-nearest, so the
intended exact root classification never happens.

Issue [#271](https://github.com/gradientspaceai/gtengine-js/issues/271). Port: preserved.

### `Arc2.h`

`Contains(P, epsilon)`'s comment promises that a negative epsilon "behaves as if
a value of zero was passed", but `std::fabs(length - radius) <= epsilon` can
never hold for `epsilon < 0`, so the function returns false for every point.

Issue [#155](https://github.com/gradientspaceai/gtengine-js/issues/155). Port: preserved.

### `ASinEstimate.h`, `ChebyshevRatioEstimate.h`, `LogEstimate.h`, `SqrtEstimate.h`, `TanEstimate.h`, `Slerp.h`

A group of documentation and estimate-table findings, all preserved (issues
[#57](https://github.com/gradientspaceai/gtengine-js/issues/57) and [#427](https://github.com/gradientspaceai/gtengine-js/issues/427)):

- `TanEstimateRR`'s comment claims `r` in `[-pi, pi]`, but `remainder(x, pi)`
  yields `[-pi/2, pi/2]`, leaving two branches dead (verified including the tie
  cases).
- `ChebyshevRatioEstimate`'s R-variant comment misstates the x-domain.
- `SqrtEstimateRR(0)` produces a meaningless nonzero value despite the
  documented `x >= 0` constraint.
- `GetLogEstimateMaxError` forwards straight to `GetLog2EstimateMaxError`,
  reporting the log2 bound. Since
  `log(x) - LogEstimate(x) = log(2)*(log2(x) - Log2Estimate(x))`, the tight bound
  is about 1.443 times smaller. It is still a valid upper bound; `ExpEstimate.h`
  by contrast carries its own scaled error table.
- `C_ASIN_EST_MAX_ERROR` disagrees with `C_ACOS_EST_MAX_ERROR` in the trailing
  digits although the two are provably the same quantity (degree 8:
  `3.5952707963527e-9` vs `3.5952707477805e-9`, relative difference 1.4e-8).
  This looks like two separate minimax solver runs.
- `Slerp.h`'s header formula is written
  `[sin((1-t)*theta)*q0 + sin(theta)*q1] / sin(theta)`; the second numerator term
  must be `sin(t*theta)*q1`, otherwise the formula does not reproduce `q1` at
  `t = 1`. The code is correct.

Every other published max-error bound in the estimate group holds when
re-measured on a dense grid.

### `BasisFunction.h`

**Symptom.** `GetIndex` returns an index for which every `invD` is zero, so the
basis evaluates to all zeros and partition of unity fails exactly at a domain
endpoint.

**Cause.** The `t <= tmin` and `t >= tmax` shortcuts ignore knot multiplicities.
For a floating knot vector with `knots[degree] == knots[degree+1]` the returned
index does not satisfy `knots[i] <= t < knots[i+1]`.

**Reproduction.** `numControls = 3`, `degree = 1`, unique knots `(0,1)`, `(1,2)`,
`(2,2)`.

Related: `Create` validates interior multiplicities against `d+1` while the class
documentation requires `m[j] <= d`; the looser check is what makes the above
reachable. Nothing validates that the multiplicities sum to
`numControls + degree + 1`, so a short input leaves `mKnots` partly
default-initialised. `GetValue`'s comment promises zero outside
`[minIndex, maxIndex]`, but `Evaluate` only overwrites the current support, so
stale nonzero values are returned.

Issue [#415](https://github.com/gradientspaceai/gtengine-js/issues/415). Port: preserved.

### `BSNumber.h`

**Symptom.** `BSNumber("x")` silently yields 72 and `BSNumber("+/")` yields -1.

**Cause.** In `ConvertToInteger` the format assertions sit inside
`if (number.size() > 1)`, so a one-character non-digit string is never
validated and the character's ASCII value minus `'0'` flows straight into the
digit accumulator.

**Suggested fix.** Move the digit-range assertion outside the size test.
Port: fixed. Issue [#95](https://github.com/gradientspaceai/gtengine-js/issues/95).

### `BSPPolygon2.h`, `BSPTree2.h`

**1. `operator&` / `operator-` throw on an empty result.** `Finalize()` is called
unconditionally and asserts `edges.size() > 0`, so an empty intersection throws
instead of returning an empty polygon; `A - A` throws too.

**2. `SplitEdge` uses `std::map::insert`,** which is a silent no-op if the
shortened edge is already a key. The edge index `eIndex` is then left unmapped
and `GetNumEdges()` falls below `mEArray.size()`.

**3. Documentation contradiction.** The `BSPPolygon2` class comment states point
classification as `(-1 inside, 0 on polygon, +1 outside)`; `BSPTree2::PointLocation`
states `(-1 outside, 0 on, +1 inside)`. With the negative-side interior
convention the class comment is correct and the `BSPTree2` comment is inverted.

Issues [#169](https://github.com/gradientspaceai/gtengine-js/issues/169), [#388](https://github.com/gradientspaceai/gtengine-js/issues/388). Port: preserved (items 1-2), documentation corrected (item 3).

### `BSPrecision.h`

**Symptom.** `operator+`'s `maxBits` is under-estimated for heterogeneous
operand sets, and the operation is not symmetric.

**Cause.** In the `bsp0.bsn.maxExponent >= bsp1.bsn.maxExponent` branch the bit
count is `bsp0.bsn.maxExponent - bsp1.bsn.minExponent + 1`, pairing *one*
operand's `maxExponent` with the *other* operand's `minExponent`. That is
correct only when the operand with the larger `maxExponent` also has the larger
`minExponent`, in particular when both operands come from the same set (the
documented usage).

**Reproduction.**
- Under-estimate: `IS_DOUBLE + IS_FLOAT` reports `minExponent = -1074`,
  `maxExponent = 1023`, `maxBits = 1173` (`1023 - (-149) + 1`). But `2^-1074`
  (double set) plus `2^127` (float set) is an exact sum whose odd mantissa has
  `127 + 1074 + 1 = 1202` bits. The reported `maxBits` is also inconsistent with
  the reported exponent span, `1023 - (-1074) + 1 = 2098`.
- Asymmetry: with `a = (-17, 20, 24)` and `b = (-5, 20, 12)`, `a + b` gives
  `maxBits = 27` while `b + a` gives `maxBits = 38`.

A fuzz over about 400k random members of the sets confirms the `mul` bound holds
in general and the `add` bound holds whenever both operands come from the same
set; every failure involves two different sets. Note that the `N` word counts
published in the `PrimalQuery2.h`/`PrimalQuery3.h` comments are computed through
chains that mix operand sets, so they may be under-estimates for the worst case.

Issue [#366](https://github.com/gradientspaceai/gtengine-js/issues/366). Port: preserved (the port's arbitrary precision is `bigint`-backed, so nothing allocates from these numbers).

### `BSplineReduction.h`

**Symptom.** The reduction is accurate only to about 1e-2 on inputs that lie
exactly in the output space.

**Cause.** The Gram-matrix entries `integral N_i N_j` are computed with
`Integration<Real>::Romberg(8, ...)` over the whole support of each product,
straddling knots where the integrand's derivative is discontinuous. Romberg's
Richardson extrapolation assumes smoothness, so the entries carry about 1e-3
error.

**Reproduction.** A straight degree-1 spline, which lies exactly in the output
space, is reproduced only to ~1e-2, and the rows of `A^-1 B` sum to 1 only to
that accuracy. A per-knot-span Gauss or Romberg integration, where the integrand
is polynomial on each span, recovers the exact controls to 1e-5.

Issue [#169](https://github.com/gradientspaceai/gtengine-js/issues/169). Port: preserved, with the accuracy documented in the header.

### `BSRational.h`

**1. The two-`BSNumber` constructor produces an invalid zero.**
`mNumerator.mBiasedExponent -= mDenominator.GetExponent()` runs unconditionally,
so a zero numerator becomes a `BSNumber` with sign 0, uinteger 0, but a nonzero
biased exponent: an invalid encoding that compares unequal to canonical zero.
Reachable from `BSRational(0.0, 1024.0)` and from any arithmetic whose numerator
cancels.

**2. The string constructor's `SetSign(sign)` is unconditional,** so `-0.0`
yields an invalid negative zero (same class as the `BSNumber` finding above).

**3. The string constructor's `x.` branch is dead.** `decimal < fpNumber.size()`
is always true when `find` succeeded, so the branch is never taken; the `x.y`
branch produces the same result.

Issue [#168](https://github.com/gradientspaceai/gtengine-js/issues/168). Port: items 1 and 2 fixed (guarded on a nonzero numerator, matching the `BSNumber` guard), item 3 dropped.

### `BVTree.h`, `BVTreeOfTriangles.h`, `OBBTree.h`, `OBBTreeOfPoints.h`, `OBBTreeOfSegments.h`, `OBBTreeOfTriangles.h`

**1. `BVTreeOfTriangles::Execute` and `OBBTreeOfTriangles::Execute` drop coincident hits (result-corrupting).**
Hits are collected in a `std::set<Intersection>` whose `operator<` compares only
`parameter`:

```cpp
bool operator<(Intersection const& other) const
{
    return parameter < other.parameter;
}
```

Two triangles hit at the same parameter, the routine case of a line, ray or
segment passing through a shared edge or vertex of a closed mesh, are
set-equivalent, so all but one are silently dropped, and which `triangleIndex`
survives depends on the tree partition and traversal order.

**Suggested fix.** Order lexicographically by `(parameter, triangleIndex)`, or
use `std::multiset`. Port: fixed in both files.

**2. `BVTree::GetLeafIndices` never tests a leaf's own bounding volume.**
Only interior nodes get the `linearBoundaryVolumeQuery` test; the leaf branch is
an unconditional `nodeIndices.push_back(nodeIndex)`. A leaf is therefore
reported whenever its *parent's* volume is hit, and a height-0 tree reports its
root for every query, contradicting the documented "leaf nodes whose bounding
volumes are intersected". Only work is wasted, because the derived classes
re-test each reported leaf's primitive exactly. Preserved.

**3. `BVTreeOfTriangles::IntersectSegmentTriangle` reports the wrong parameter scale.**
`BVTree.h` documents `SEGMENT_QUERY` results as the parameter `t` of
`(1-t)*P + t*Q`, `0 <= t <= 1`, but the function forwards the
`FIQuery<Segment3, Triangle3>` output unchanged, and that query reports `s` of
the centered form `C + s*D` with `|s| <= e`. The older
`OBBTreeOfTriangles.h` does convert:

```cpp
triResult.parameter = result.parameter / Length(Q - P) + static_cast<T>(0.5);
```

so the newer file lost the conversion and its three query types now report
parameters on different scales. Ordering along the segment is unaffected.
Preserved.

**4. `ComputeInteriorBox` seeds `pmin = pmax = Vector3<T>::Zero()`** in
`OBBTreeOfPoints`, `OBBTreeOfSegments` and `OBBTreeOfTriangles` instead of
seeding from the first projection. This is still a correct min/max:
`pmin <= 0 <= pmax` is a loop invariant and the projections always straddle zero
because the projection origin is the mean of the node's centroids. The only
consequence is that the box is forced to contain that mean. Preserved.

**5. `OBBTree.h`: `height > 31` is clamped to 31 and then preallocates `2^32 - 1` nodes,**
an out-of-memory allocation in practice (hit while writing tests; unusable in
C++ as well). A height cap derived from the input size, or lazy node allocation,
would fix it. Preserved.

**6. Minor.** `OBBTreeOfTriangles::Execute` never tests leaf boxes against the
linear component (correct results, wasted work);
`OBBTreeOfTriangles::IntersectSegmentTriangle` divides by `Length(Q - P)`
unguarded; `OBBTreeOfSegments` runs an unguarded `Normalize` on a zero-length
segment, which collapses to a zero frame.

Issues [#103](https://github.com/gradientspaceai/gtengine-js/issues/103), [#167](https://github.com/gradientspaceai/gtengine-js/issues/167), [#274](https://github.com/gradientspaceai/gtengine-js/issues/274), [#387](https://github.com/gradientspaceai/gtengine-js/issues/387).

### `CholeskyDecomposition.h`, `LDLTDecomposition.h`

**1. Run-time `BlockCholeskyDecomposition<Real, 0, 0>` uses the wrong scalar offset inside blocks (result-corrupting).**
The block-level helper `GetIndex(row, col) = col + row * NumBlocks` is also used
for the *scalar* offset inside a `BlockSize x BlockSize` block
(`GMatrix::operator[]`) in `SolveLower`, `SolveUpper`, `LowerTriangularSolver`
and `SubtractiveUpdate`. The correct in-block offset is `col + row * BlockSize`.
Whenever `BlockSize != NumBlocks` the run-time class reads and writes the wrong
scalars, and for `BlockSize < NumBlocks` it indexes past the end of the block's
`std::vector` (for example `BlockSize = 2, NumBlocks = 3` reaches index 4 of a
4-element block). The compile-time specialization is the same code written
correctly as `Lrc(i,j)` / `Acc(j,i)` / `Ark(j,i)`. Port: fixed; tests at
`(2,3)` and `(3,2)` require the block factor to equal the unblocked one.

**2. `BlockLDLTDecomposition<T>::Convert(BlockVector, GVector&)` asserts the wrong dimension.**
The per-block check is `current.GetSize() == NumBlocks`; block vectors have
`BlockSize` components, as the class docs, the sibling `Convert` and `Solve`'s
own checks all state, so a valid `(BlockSize = 3, NumBlocks = 2)` block vector
can never be converted back. Port: fixed.

**3. Documentation.** `LDLTDecomposition.h`'s comments claim positive-definite
input and positive `D`, but `Factor` fails only on an exactly zero pivot, so
symmetric indefinite matrices factor fine with negative `D` entries. The
run-time `CholeskyDecomposition` comments say "Ensure that N > 0 at run time"
but there is no check (the LDLT class has one).

Issues [#209](https://github.com/gradientspaceai/gtengine-js/issues/209), [#478](https://github.com/gradientspaceai/gtengine-js/issues/478).

### `CLODPolyline.h`

**Symptom.** After all collapses, the surviving edge of an open polyline is
built from an out-of-range index and every reduced level of detail is unusable.

**Cause.**

```cpp
edges[0] = collapses[0];
edges[1] = collapses[0] + 1;
```

Both endpoints of an open polyline receive `std::numeric_limits<Real>::max()`
weight, and `MinHeap` breaks that tie such that `collapses[0] == numVertices - 1`.
Then `edges[1] == numVertices`, and `ReorderVertices` reads
`permute[numVertices]`, one past the end.

**Suggested fix (applied in the port).**
`vIndex = (collapses[0] != numVertices - 1 ? collapses[0] : collapses[1]);`
then `edges = { vIndex, vIndex + 1 }`.

Issue [#182](https://github.com/gradientspaceai/gtengine-js/issues/182). Port: fixed.

### `Cone.h`

`CreateMesh`'s `tNumExtra = 0.5*(1 + rMax)/(1 + rMin) - 1` is an undocumented,
non-scale-invariant replacement for the original `0.5*rMax/rMin - 1` (which is
infinite at `hMin == 0`); it affects only sample density. There is also no lower
bound on `numMinVertices`.

Issue [#268](https://github.com/gradientspaceai/gtengine-js/issues/268). Port: density rule preserved, `numMinVertices >= 3` asserted.

### `ConformalMapGenus0.h`

**1. Null dereference on open meshes (~L89).** `element.second->T[j]` for
`j = 0, 1` is dereferenced unconditionally while building the Laplacian, and
`T[1]` is null for a boundary edge. The closed-mesh precondition is documented
but never checked. Port: asserts closed-manifold input.

**2. `LogAssert(fmin > 0 && fmax < 0)` (~L375) can fire on valid input.**
`fmin = poly2(0)` is exactly zero whenever a puncture-triangle vertex lands on
the plane origin, and the origin is the subtracted centroid computed a few lines
earlier, so the configuration is reachable rather than pathological. No
obviously correct alternative bracket presents itself. Preserved.

**3. Convergence flag off by one (~L155/175).** `iterations >= maxIterations`
reports non-convergence when the conjugate-gradient solve converged on exactly
the last permitted iteration. Output coordinates are unaffected.

**4. Dead store** to `fmax` inside the bisection loop (~L386).

Issue [#318](https://github.com/gradientspaceai/gtengine-js/issues/318).

### `ConstrainedDelaunay2.h`, `ConvexHull3.h`

**1. `ConstrainedDelaunay2<T>::operator()` never clears `mInsertedEdges`.**
A second data set processed by the same object inherits the first set's
inserted-edge keys, which then index unrelated vertices in the new
triangulation. Port: cleared at the start of compute.

**2. `ConstrainedDelaunay2<T>::Insert` reads `duplicates[edge[i]]` before range-checking `edge[i]`.**
Only the substituted values are asserted, so an out-of-range input index is
undefined behaviour the assert cannot catch. The deprecated specialization has
the mirror-image flaw: it checks the raw indices but never substitutes
duplicates, so a duplicate vertex index trips "Failed to find vertex in graph."
Port: raw indices validated first.

**3. `ConstrainedDelaunay2<T>::Retriangulate` is not Delaunay inside the strip.**
The strip is filled by a minimum-pseudosquared-distance bisection rather than
the empty-circumcircle rule, so the output is a valid constrained triangulation
that is not constrained-Delaunay near inserted edges. Neither the class name nor
the header comments say so. Preserved.

**4. `ConvexHull3::GetHull()` documentation** claims `E = T/2` satisfies
`V - E + T = 2`; for a closed triangle mesh the identity is `E = 3T/2` (cube
hull: `8 - 18 + 12 = 2`, not `8 - 6 + 12 = 14`).

**5. `ConvexHull3::SelectSplit`** binds `crV2x`/`crV2y` as `const&` into the
scratch pool and then mutates those pool slots inside the candidate loop, relying
on the references to observe the mutation. Correct, but invisible at the call
site.

Issue [#325](https://github.com/gradientspaceai/gtengine-js/issues/325).

### `ContAlignedBox.h`, `ContAlignedBox2Arc2.h`, `ContCircle2.h`, `ContCone.h`, `ContPointInPolygon2.h`, `ContSphere3.h`

- `ContCircle2::GetContainer` and `ContSphere3::GetContainer` execute
  `circle.center = points[0]` before any count check and then `/= (Real)numPoints`:
  an out-of-bounds read plus a divide by zero for an empty `std::vector` via the
  convenience overload. The same pattern recurs in `ContCapsule3`,
  `ContCylinder3`, `ContLozenge3` and `ContOrientedBox2`. Port: throws.
- `ContPointInPolygon2`'s constructor comment, "The class object stores a copy of
  'points'", contradicts the member `Vector2<Real> const* mPoints`, a bare
  pointer with no copy.
- `ContAlignedBox::MergeContainers` and the circle/sphere `GetContainer` can only
  return `true`; the `bool` return is vestigial.
- `ContAlignedBox2Arc2`'s `0 < numPoints && numPoints <= 6` is always true and
  the trailing `return false` is unreachable.
- `ContCone.h` relies on transitive includes for `Vector` and `Dot`.

Issues [#106](https://github.com/gradientspaceai/gtengine-js/issues/106), [#174](https://github.com/gradientspaceai/gtengine-js/issues/174), [#292](https://github.com/gradientspaceai/gtengine-js/issues/292).

### `ContEllipse2.h`, `ContEllipsoid3.h`, `ContOrientedBox3.h`, `ContScribeCircle2.h`

**1. `GetContainer` turns a zero eigenvalue into an infinite extent.** The
comment says nonpositive eigenvalues are adjusted, but only strictly negative
ones are; an exactly zero eigenvalue (collinear or coplanar input) divides by
zero, the extent becomes infinite, and `InContainer` then accepts every point.

**2. `MergeContainers` does not contain its inputs.** The result is inscribed in
the bounding box of the projected intervals, so it generally excludes parts of
the input ellipses. Counterexample: merging unit circles centred at `(-1,0)` and
`(1,0)` yields extents `(2,1)`, which excludes `(-1,1)`, a point on the first
circle. A correct merge is a different algorithm, so this is preserved, but
callers relying on the containment contract should know it does not hold.

**3. `ContOrientedBox3::MergeContainers` and `ContEllipsoid3::MergeContainers`
assume `det(rot) = +1`.** Eigenvector frames from `GetContainer` are not
guaranteed right-handed, so round-tripping a fitted container through merge can
produce a garbage frame (quaternion conversion of a reflection).

**4. `ContOrientedBox3::GetContainer` reads `points[0]` with no count guard,**
reachable because `ApprGaussian3::Fit(0, nullptr)` can report success.

**5. `ContOrientedBox3::MergeContainers` seeds `pmin`/`pmax` at zero** rather
than the first projected vertex; sound only because the provisional centre lies
between the boxes.

**6. `ContScribeCircle2::Inscribe` writes a degenerate circle into the output
before returning `false`.**

**7. `InContainer` rejects the defining extremal point by one ulp.**
`extent[j] = sqrt(maxValue / D[j])` puts the extremal point exactly on the
boundary and `Length(standardized) <= 1` has no tolerance (observed
`Q = 1.0000000000000004`).

Issues [#292](https://github.com/gradientspaceai/gtengine-js/issues/292), [#409](https://github.com/gradientspaceai/gtengine-js/issues/409).

### `ContEllipse2MinCR.h`

**Symptom.** `MaxProduct` returns a NaN component.

**Cause.** In the `0.5 < A[iYMin][0] * x0` branch,
`D[1] = (1 - a0 * D[0]) / A[iYMin][1]` divides by zero whenever the walk steps
onto a vertical constraint line, which is exactly when that branch is entered
(`a0 * x0 = 1 > 1/2`). Any input point lying on the second ellipse axis produces
such a line.

**Reproduction.** `C = (0,0)`, `R = I`, points `(3,0), (0,2), (2,1.4), (-2.5,-1)`
gives `D = (1/9, NaN)`; the correct answer is `(1/9, 1/4)`.

**Suggested fix (applied in the port).** Evaluate the shared vertex's `y` on the
previous line of the walk. Issue [#234](https://github.com/gradientspaceai/gtengine-js/issues/234). Port: fixed.

### `ContEllipsoid3MinCR.h`

**1. `LogAssert(numer >= 0)` in `FindFacetMax` / `FindEdgeMax` contradicts its own comment.**
The adjacent comment says "some numerical error may make this a small negative
number. In that case set tmax = 0". An already-active plane routinely has slack
about -1e-16, so upstream either aborts on ordinary random clouds or, with
assertions disabled, gets a slightly negative `t`, fails the `0 <= t` test,
ignores the blocking plane and walks straight out of the feasible region; slack
was measured reaching -8.1 with non-terminating recursion. Port: clamps `numer`
to 0 as the comment prescribes.

**2. The facet/edge walk does not terminate (result-corrupting).**
With three or more planes active and every candidate step zero or numerically
negligible, the mutual recursion cycles with the point unchanged until the stack
is exhausted; about 0.3% of random integer clouds. Traced example:
`points = {(1,-5,2), (2,-2,2), (3,-2,-1), (-1,1,5), (0,3,-2)}`, `C = 0`, `R = I`
repeats `facet(0) -> facet(2) -> facet(5) -> edge(5,2) -> facet(0)` forever with
`D` bit-identical. The upstream assert in item 1 usually masks this; a variant
with neither assert nor clamp still overflows the stack, so the cycle belongs to
the algorithm. Port: bounds the walk to `4*N + 64` steps and stops at the
current (feasible) point; over 6000 random clouds, 5947 results are bit-identical,
26 previously crashing clouds now return, and none is infeasible.

**3. Latent out-of-bounds.** In `FindFacetMax`'s degenerate `else` branch, if no
plane blocks, `tMax` stays at `numeric_limits::max()` and the code recurses with
`plane0 = plane1 = -1`, indexing `A[-1]`. Made essentially unreachable by the
1e-12 jitter.

**4. The walk is a heuristic (preserved).** Even with the clamp,
`FindEdgeMax`'s "tmax == 0, so return" case stalls at a polytope vertex on about
4% of random clouds (worst product 0.096 of the true maximum), giving a
containing but non-minimal ellipsoid. Containment always holds, so this is a
fit-quality limit of the algorithm rather than a correctness bug.

Issues [#234](https://github.com/gradientspaceai/gtengine-js/issues/234), [#409](https://github.com/gradientspaceai/gtengine-js/issues/409).

### `ContLozenge3.h`

**Symptom.** The returned lozenge does not contain the input points.

**Cause.** `GetContainer` centres the fitted rectangle on a *corner* of the
fitted parameter interval instead of its midpoint:

```cpp
lozenge.rectangle.center = box.center + aMin * box.axis[2] + bMin * box.axis[1];
lozenge.rectangle.extent[0] = 0.5 * (aMax - aMin);
lozenge.rectangle.extent[1] = 0.5 * (bMax - bMin);
```

`Rectangle<3,Real>` is a *centred* rectangle (`R(s0,s1) = C + s0*A0 + s1*A1`,
`|si| <= ei`), so the result is shifted by its own half-extents. The
degenerate-direction branches in the same function correctly use
`0.5*(bMin + bMax)`, so the code is internally inconsistent; it looks like a
leftover from the Wild Magic corner-origin rectangle (stale origin wording also
survives in the `Rectangle.h` and `Lozenge3.h` comments).

**Reproduction.** Planar grid `{-3..3} x {-2..2} x {0}` yields `radius = 0`,
`extent = (3, 2)`, `center = (-3, -2, 0)`, leaving the input point `(3, 2, 0)`
at distance 3 from the container.

**Suggested fix.** Use interval midpoints in every branch. Port: fixed; worst
containment excess over random clouds is about 1e-15.

Issue [#174](https://github.com/gradientspaceai/gtengine-js/issues/174).

### `ContPointInPolyhedron3.h`

**Symptom.** A unit cube described by its six quadrilateral faces reports its own
centre `(0.5, 0.5, 0.5)` as OUTSIDE for CONVEX methods 0, 1 and 2.

**Cause.** `Face::indices` is declared `std::array<int32_t, 3>`, but `ContainsC0`
and `SharedContains` read `face->indices.size()` as the face's true vertex count.
This looks like a regression from an older `std::vector<int32_t>`: any face with
more than 3 vertices is silently truncated to its first three.

**Suggested fix.** Use a dynamic index array. Port: fixed; six (type, method)
combinations now agree with analytic tests for the cube and a non-convex L-prism.

Minor in the same header: `mNumPoints` is a dead member, and `Contains` returns a
silent `false` for unsupported (type, method) pairs (for example SIMPLE with
method >= 2), indistinguishable from "outside".

Issue [#343](https://github.com/gradientspaceai/gtengine-js/issues/343).

### `ConvexHull2.h`

`GetTangent` can silently return non-tangent indices. If the `size0 + size1`
bounding loop expires without the tangency tests passing, the function returns
whatever indices it last held, corrupting the hull merge with no diagnostic.
Preserved (adding a throw would change behaviour in cases where upstream happens
to succeed); randomized monotone-chain cross-checks never triggered it.

Issue [#277](https://github.com/gradientspaceai/gtengine-js/issues/277).

### `ConvexPolyhedron3.h`, `TriangulateEC.h`

- `ConvexPolyhedron3`'s constructor validates `indices.size() >= 12` but never
  `indices.size() % 3 == 0`, unlike `Polyhedron3.h` on the same data; a trailing
  partial triangle is silently dropped.
- `TriangulateEC::VertexList::RemoveR`: removing the *only* reflex vertex sets
  `mRFirst = -1` but leaves `mRLast` pointing at the removed vertex, breaking the
  "empty iff both are -1" invariant the class's own `LogAssert` relies on.
  Appears unreachable in practice.
- `TriangulateEC::DoEarClipping` closes the ear ring with
  `V(mEFirst).ePrev = mELast` without checking `mEFirst != -1`, so malformed
  input surfaces as an "Index out of range" assertion rather than a meaningful
  error.
- Stale comments: "modifty"; "triangles that are simple" where polygons are
  meant.

Issue [#175](https://github.com/gradientspaceai/gtengine-js/issues/175). Port: preserved.

### `CubicRootsQR.h`, `QuarticRootsQR.h`

**Symptom.** Real roots of even multiplicity are silently dropped.

**Cause.** `GetQuadraticRoots` emits roots from the deflated 2x2 block only when
`discriminant >= 0`. An even-multiplicity root gives a discriminant that is
exactly zero in exact arithmetic, so rounding pushes it negative and the block
contributes nothing, indistinguishable from a genuine complex pair.

**Reproduction.** `(x-6)^2 (x+1)` reports only `-1`; `(x-1)^2 (x-2)^2` reports
only the pair near 2.

A fix would mean inventing a discriminant tolerance, so this is preserved and
documented. Issue [#488](https://github.com/gradientspaceai/gtengine-js/issues/488).

### `CurvatureFlow2.h`

Line 48 computes the update numerator with `- 0.5*uxy*ux*uy`, but the
mean-curvature-flow numerator `uxx*uy^2 - 2*uxy*ux*uy + uyy*ux^2` has coefficient
2, and `CurvatureFlow3`'s own expansion uses 2 for its mixed terms. The 2D
filter therefore under-weights the mixed derivative by a factor of 4.

Issue [#123](https://github.com/gradientspaceai/gtengine-js/issues/123). Port: preserved, pinned by a hand-computed test.

### `CurveExtractor.h`, `SurfaceExtractor.h`, `BitHacks.h`, `Array4.h`

- `CurveExtractor::MakeUnique` does not re-sort remapped edges. `Edge` sorts its
  two indices at construction and `operator<` compares positionally, but
  `MakeUnique` overwrites `edge.v[i]` in place without restoring order, so an
  edge remapping to `(1,0)` is a distinct `std::map` key from `(0,1)` and
  duplicates survive dedup.
- `SurfaceExtractor::MakeUnique` has the same shape: it keeps rotated duplicate
  triangles and emits non-canonical rotations.
- `BitHacks.h`'s header comment states `GetTrailingBit(10) = 2`. 10 is `0b1010`,
  whose lowest set bit is at index 1, and the code correctly returns 1. The same
  sentence writes the value as `0x0000000000001010`, mixing a hex prefix with
  binary digits; `GetLeadingBit(10) = 3` is correct.
- `Array4.h`'s private `SetPointers(Array4 const&)` ends with
  `// else 'other' is an empty Array3.`, copy-pasted from `Array3.h`.

Issues [#67](https://github.com/gradientspaceai/gtengine-js/issues/67), [#362](https://github.com/gradientspaceai/gtengine-js/issues/362), [#363](https://github.com/gradientspaceai/gtengine-js/issues/363), [#439](https://github.com/gradientspaceai/gtengine-js/issues/439). Port: preserved.
### `Cylinder3.h`, `Polyhedron3.h`, `Parallelepiped3.h`, `Parallelogram2.h`, `Torus3.h`, `SegmentMesh.h`, `AlignedBox.h`, `OrientedBox.h`, `Triangle.h`

Primitive-level contract and documentation findings, all preserved:

- `Cylinder3.h`: the comment claims the default constructor sets the axis to
  `(0,0,1)`; it uses `Line3<T>()`, whose direction is `(1,0,0)`.
- `Polyhedron3.h`: the class comment says `numIndices` "must be 6 or larger"; the
  constructor requires `>= 12`. Separately, `auto vertexPool = GetVertices();` in
  all three geometric queries deduces by value from a `const&` return,
  deep-copying the entire vertex pool on every call; `auto const&` was surely
  intended.
- `Parallelepiped3.h` / `Parallelogram2.h`: `GetVertices` is documented as
  "counterclockwise order" but emits bit-pattern order (0, 1, 3, 2 around the
  boundary). (`Rectangle::GetVertices` emits the same bit-pattern order and its
  comment is correct; see "Claims withdrawn or corrected".)
- `Parallelepiped3.h`: the right-handedness assert `DotCross(...) > 0` accepts a
  numerically degenerate basis, because a repeated axis rounds to a tiny positive
  triple product. `Parallelogram2`'s `DotPerp(A0, A0)` cancels exactly, so only
  3D is affected.
- `Torus3.h`: `GetParameters` documents `u, v` in `[0, 2*pi)` but returns
  `atan2` values in `[-pi, pi]`. Harmless for `Evaluate` (periodic) but visible
  through `IntrLine3Torus3FIResult.torusParameter`. `Evaluate` guards the
  second-order derivatives with `maxOrder == 2` rather than `>= 2`, so
  `maxOrder = 3` silently returns only the first-order jet.
- `SegmentMesh.h`: the `CONTIGUOUS_CLOSED` enum doc says `S[i] = {i, (i+1) % L}`,
  but the constructor loop stores `S[i] = {i-1, i}`, the same segment set rotated
  so the wrap-around `<V[L-1], V[0]>` lands at index 0 instead of `L-1`.
  `DISJOINT` documents "L must be an even number" but the constructor asserts
  only `size() >= 2` and builds `floor(L/2)` segments, silently dropping the last
  vertex for odd `L`.
- Comment typos: `AlignedBox.h` `GetVertices` ("whern", `vertex[i][d = max[d]`);
  `Triangle.h` ("sets the/ vertices"); `OrientedBox.h`'s class comment writes box
  coordinates with only two terms for an N-dimensional class.

Issues [#78](https://github.com/gradientspaceai/gtengine-js/issues/78), [#155](https://github.com/gradientspaceai/gtengine-js/issues/155), [#455](https://github.com/gradientspaceai/gtengine-js/issues/455), [#484](https://github.com/gradientspaceai/gtengine-js/issues/484).

### `Delaunay2.h`, `Delaunay3.h`

**1. Hardcoded `epsilon = 0` misclassifies degenerate input (result-corrupting).**
`Delaunay2.h` (L907) and `Delaunay3.h` (L978) construct `IntrinsicsVector2/3`
with `static_cast<T>(0)`. The intrinsics measure point distances against a
*normalized* frame, so an exactly degenerate point set has distances that are
nonzero by an ulp and the strict `maxDistance <= 0` test fails: collinear input
is classified as dimension 2 and coplanar or collinear input as dimension 3 or
2. The downstream `LogAssert("The tetrahedron should not be degenerate.")` only
checks that the mesh returned a feature, not that it has area or volume.

**Reproduction (small integers).** `Delaunay2` on `(-3,-9), (0,0)` returns
success with indices `[0, 1, 1]`, a degenerate seed triangle; more collinear
points can throw. `Delaunay3` on the five coplanar points
`(0,0,0), (1,0,1), (0,1,2), (1,1,3), (2,1,4)` returns success with zero
tetrahedra; other inputs throw "Attempt to create nonmanifold mesh." or
"Unexpected insertion failure."

**Suggested fix.** Classify the dimension with the classes' own exact `ToLine` /
`ToPlane` predicates. Port: fixed; inputs upstream already handles produce
bit-identical results.

**2. `Delaunay3<T>` duplicate detection can never match (result-corrupting).**
`ProcessedVertex` (~L1486) hashes *and* compares the `location` field:

```cpp
return HashValue(v.vertex[0], v.vertex[1], v.vertex[2], v.location);
...
return v0.vertex == v1.vertex && v0.location == v1.location;
```

The lookup key is constructed as `ProcessedVertex(mVertices[i], i)`, so a
repeated coordinate triple at a later index carries a different `location` and
the find always misses. `mDuplicates[]` then remains the identity map,
`GetNumUniqueVertices()` equals `mNumVertices`, and `Update(i)` runs for a point
already in the mesh, removing that point's incident tetrahedra (the point lies
on their circumspheres) and rebuilding them against a coincident vertex index:
zero-volume tetrahedra or an insertion assert. `Delaunay2.h` (~L1331) has the
correct form, hashing and comparing the vertex only. Port: fixed.

**3. `GetNumVertices()` disagrees with `GetVertices()` after degenerate input.**
For 0- or 1-dimensional input (2D) or dimension 0/1/2 (3D), `operator()` clears
`mIRVertices` and returns early, so `GetNumVertices()` reports 0 while
`GetVertices()` still returns the caller's non-empty vertex pointer. A client
iterating `GetVertices()` bounded by `GetNumVertices()` silently sees nothing.
Port: fixed (returns the stored vertex count).

**4-6 (minor and doc).** Several non-void functions fall off the end after a bare
`LogError` (formally undefined behaviour if `LogError` ever returns). The
deprecated `Delaunay2<InputType, ComputeType>::Update` comments out its own
`LogError` diagnostic, silently swallowing an internal-state error. Stale
comments: "inside all four edges" in the 2D `GetContainingTriangle`;
"constructed by the call to IntrinsicsVector2{T}" in `Delaunay3.h` (~L1021),
which should read `IntrinsicsVector3<T>`.

Issues [#277](https://github.com/gradientspaceai/gtengine-js/issues/277), [#283](https://github.com/gradientspaceai/gtengine-js/issues/283), [#391](https://github.com/gradientspaceai/gtengine-js/issues/391).

### `DistCircle3Circle3.h`

Six result-corrupting defects, all fixed in the port, plus three preserved
quirks. This is the single most affected file in the library.

**1. Mirror-symmetric configurations miss the true minimum.** The solver squares
`H(cs, sn) = p6(cs) + sn * p7(cs)` into `phi = p6^2 - (1 - cs^2) * p7^2`.
Whenever the configuration is mirror-symmetric, `p7` is identically zero in
exact arithmetic and `phi`'s minimizing roots are double roots. Roundoff in
`PrepareCircles` leaves `p7` at about 1e-16 instead of 0, splitting each double
root into two simple roots about 1e-32 apart, which the bisection-based root
finder cannot resolve, so the global minimizer is never evaluated as a
candidate. Reproduction: the unit circle in the `z = 0` plane versus the unit
circle in the `x = 3` plane returns `sqrt(5) = 2.2360679...` where the true
distance is `2.1622776601683795`. Fix applied in the port: additionally pair the
roots of `p6` with `sn = +-sqrt(1 - cs^2)` (upstream's own `p7 == 0` rule), which
is safe because every candidate is then a genuine point of circle1 whose
distance is evaluated exactly. Caution for any fix: at a root of `p6`, upstream's
quotient `sn = -p6/p7` must NOT be used; it yields `sn = 0` and a non-unit
`(cs, sn)`, that is, a point off the circle, which does under-report.

**2. `PrepareCircles`'s normal alignment is a no-op for normals orthogonal to z.**
It aligns with `if (normal[2] < 0) negate`, so parallel-plane circles with
anti-parallel normals in the xy-orthogonal case are answered by the polynomial
path: 3.7385 instead of 5.1851 on the pinned case.

**3. `sn = -p6(cs)/p7(cs)` is on the unit circle only for an exact root of `phi`.**
The roots that matter are double roots, so the reported point can be 0.43 inside
its circle, with a distance 0.22 below the true minimum.

**4. Concentric and coaxial circles make `phi` and `p6` perfect squares,** so the
sign-change bisection finds nothing and upstream reads a default `ClosestInfo`,
reporting distance 0 at the shared axis point (0 instead of 0.3 for radii 0.2 and
0.5 in perpendicular planes).

**5. `p6 == p7 == 0`** for coaxial circles in parallel planes that slip past the
exact `normal[2] < 1` test, and for a zero-radius circle on the axis: upstream
reports "Unexpected degree for p6".

**6. `DoQueryParallelPlanes` projects with the rounded rotated normal.**
`compProj = D - Dot(N,D)*N` uses the rotated, only approximately unit normal, so
for coaxial circles the residue `(1 - |N|^2) Dot(N,D) N` points along the axis;
`d > 0` then selects the inside-circle branch with `U` along the axis and both
closest points are reported on the common axis (the distance is still correct).
Reproduction: circle0 centre `(0,0,0)`, normal `(0,0,-1)`, radius 0.2; circle1
centre `(0,0,-0.001)`, normal `(0,0,0.9999999999999999)`, radius 0.2 gives
`distance = 0.0009999999999999996` (correct) with closest points `(0,0,-0.2)` and
`(0,0,-0.201)`, which are 0.2 away from their circles. A sweep of 20000 random
coaxial pairs found 2540 bad results, 2521 of them through upstream's own
`circle0.normal[2] >= 1` branch, about 13% of coaxial configurations. Port fix:
use the exact common normal `(0,0,1)` of the prepared frame at both call sites.

**7-10 (structural).** `rRoots.size() > 1` is used where the surviving-candidate
count is meant, fabricating a second "closest pair" at the origin when only one
candidate exists; `candidates[0]` is read with no guard when no candidate
survives filtering; the `candidates` array is a fixed 16 elements with no bound
check; `SCPolynomial` is a dead private class.

**11-13 (preserved).** `PrepareCircles` divides by an underflowed length for
subnormal centre offsets (NaN); at most two closest pairs are ever reported; and
`numClosestPairs` is not symmetric under argument swap, although the distance is.

Issues [#331](https://github.com/gradientspaceai/gtengine-js/issues/331), [#431](https://github.com/gradientspaceai/gtengine-js/issues/431), [#442](https://github.com/gradientspaceai/gtengine-js/issues/442).

### `DistLine3Circle3.h`, `DistRay3Circle3.h`, `DistSegment3Circle3.h`

**1. `PDFSection422` computes `tauHat` without dividing by `a2` (result-corrupting).**

```cpp
T tauHat = std::sqrt(std::fabs(std::pow(a1 * a3, twoThirds) - a3));
```

**Cause.** With `Phi(t) = (t + a0) - a1*t/(a2*t^2 + a3)^{1/2}` and
`G(t) = a1*t/(a2*t^2 + a3)^{1/2}`, the derivative is
`G'(t) = a1*a3/(a2*t^2 + a3)^{3/2}`. (The header comment states
`G'(t) = a1*a2/(a2*t^2 + a3)^{3/2}`, which is also incorrect: the numerator is
`a1*a3`.) Solving `G'(tauHat) = 1` gives `a2*tauHat^2 + a3 = (a1*a3)^{2/3}`,
hence

```
tauHat = sqrt(|(a1 * a3)^{2/3} - a3| / a2)
```

Since `a2 = |Cross(N, M)|^2` it equals `sin^2(angle(N, M))` for unit `N` and `M`,
so `a2 != 1` in essentially every non-perpendicular configuration. `tauHat`
serves both as a bisection bracket endpoint and directly as a reported critical
parameter (`critical.parameter[...] = +-tauHat + s`), so the `G'(0) > 1` branch
can misplace critical points and report a single spurious critical point,
missing the global minimum.

**Reproduction.** A randomized search over 29,119 three-critical-point
configurations found upstream errors up to 2.9: circle radius 3.0447, normal
`(0.3617, -0.8460, 0.3917)`, line `(3.2987, 0.9600, -0.4038) + t*(-5.1028, 1.0799, 0.4136)`
gives 2.9236 where the true distance is 0.0131. With the `/a2` the maximum error
over that family is exactly 0. The failure needs `a2 = |N x M|^2 > 1`, that is a
long direction vector, which is how the segment-circle query calls it.
`DistRay3Circle3.h` and `DistSegment3Circle3.h` inherit the error through their
critical-point clamping.

**2. `a3 == 0` gives a NaN intercept and an inverted bisection bracket.**
A line meeting the circle axis away from the centre and not parallel to the
plane throws "Invalid ordering of t-interval endpoints"; for example the unit
circle in `z = 0` with the line `(1,0,6) + t*(1,0,1)`. The port delegates to
`PDFSection421`, whose closed form is provably the `a3 = 0` solution.

**3. `Finalize` normalizes a possibly-zero projection,** returning the circle
centre as the closest circle point with distance 0. Reachable by round-off for a
line along the circle axis.

**4. Near-perpendicular lines lose every significant digit, and `Bisect`'s bracket collapses.**
The first symptom is that `Bisect`'s bracket endpoints round together and
`LogAssert(tMin < tMax)` throws. Removing that assertion is not enough: the
returned root still carries no significant digits in the same regime. Worked
example: circle `C = (0,0,0)`, `N = (-1,0,0)`, `r = 0.2`, segment
`(-3.443461197292675, 7.999999999999849, 0) -> (0, 7.999999999999964, 0)`
reports 7.800046 where the truth is 7.8.

**Cause.** The loss is not in the bisection but in `PDFSection422`'s
back-substitution `t = tau + s`. Here `s` grows like `1/|NxM|^2`, and every
bracket is an interval of width `r*|NxM|/Dot(M,M)` with an endpoint at `-a0`, so
`tau` approaches `-s` to far below `ulp(s)` and the sum cancels every digit. On
the worked example `s = -66719994479561.64`, the bracket width is 2e-15 and
`ulp(a0) = 0.0148`, so `tau + s` evaluates to `0.9921875` where the true value is
1.

**Suggested fix.** Use the algebraically identical
`t = G(tau) - Dot(M,D)/Dot(M,M)`, which never forms the cancelling sum. Applied
in the port: 987 of 4000 near-perpendicular draws were wrong by more than 1e-9
relative (worst 1.64) before, 0 of 4000 after (worst 1e-15), and
well-conditioned inputs change only at rounding level, at most 1.8e-15 absolute.
This supersedes the earlier reading of the defect as a bisection-bracket
problem, and the earlier "only half fixed" status.

**5. `Dot(NxM, NxM)` underflows to exactly 0 while `NxM != 0`, so every output is NaN.**
`Execute` branches on `NxM != 0`, but `PDFSection421` and `PDFSection422` then
divide by `Dot(NxM, NxM)`, which underflows to exactly 0 for
`|NxM| < ~1.5e-162`. Reproduction: circle of radius 0.2 in the plane `x = 0`
versus the line `(0, 4, 0) + t*(-0.05, 0, -5e-324)`. Port: fixed, by routing
such directions to the perpendicular branch, which is the exact limit of the
other two.

Issues [#247](https://github.com/gradientspaceai/gtengine-js/issues/247), [#421](https://github.com/gradientspaceai/gtengine-js/issues/421), [#495](https://github.com/gradientspaceai/gtengine-js/issues/495).

### `DistLine2Triangle2.h`

**Symptom.** The query returns `distance = NaN`, `parameter = NaN` and
barycentric coordinates `(+Infinity, -Infinity, 0)`, although the vertex sign
classification correctly reports opposite signs.

**Cause.** Line 340 computes the edge parameter as
`s = DotPerp(D, P - V[i0]) / DotPerp(D, V[i1] - V[i0])`. Algebraically the
numerator equals `ncomp[i0]` and the denominator equals `ncomp[i0] - ncomp[i1]`
(the same normal components used for the sign test), but rounding the difference
`V[i1] - V[i0]` first can make the denominator exactly 0 while `ncomp[i0]` and
`ncomp[i1]` are tiny nonzero values of opposite sign.

**Reproduction.** Line origin `(-1.7086881202421092e-14, 1.9960937490641202)`,
direction `(-2.5613545596609416e-14, 2.9921808228827933)`; triangle
`(-6.848127733652727e-14, 7.999999999999985)`, `(0, 0)`,
`(0.49999999999995026, 5.820550528229657)`. Then
`ncomp = (-2.5e-29, 6.3e-30, 1.496)` with signs `(-,+,+)` and the denominator
rounds to 0.

**Suggested fix (applied in the port).** `s = ncomp[i0] / (ncomp[i0] - ncomp[i1])`,
which the sign condition guarantees is finite and in `[0,1]`. The port uses it
only when upstream's denominator `DotPerp(D, V[i1] - V[i0])` is exactly 0 and
evaluates upstream's quotient otherwise: the two are algebraically equal but
round differently, and the C++ oracle (family `v19-distance`) showed that the
unconditional replacement moved 12.6% of ordinary results by a few ulps.

Minor in the same file: `NoCommonPoints` sets
`result.distance = min |Dot(Perp(D), V[i] - P)|`, which is scaled by `|D|` and so
wrong for non-unit directions, but `operator()` recomputes the distance from
`closest[0] - closest[1]` before returning, so the wrong value never escapes.

Issues [#441](https://github.com/gradientspaceai/gtengine-js/issues/441), [#118](https://github.com/gradientspaceai/gtengine-js/issues/118).

### `DistLine3CanonicalBox3.h`

**Symptom.** `distance = NaN` for a flat (zero-extent) box; the defect propagates
to every line, ray and segment versus aligned or oriented box query.

**Cause.** `sqrDistance` is accumulated in the cancellation-prone form
`pme^2 + tmp^2 + PpE^2 + delta*parameter`.

**Reproduction.** Extent `(0, 3.999999999999937, 3.999999999970607)`, line origin
`(0, -7.499019495845319, 0)`, direction `(-7.56e-9, -0.976, 0.216)` gives
`sqrDistance = -1.78e-15` and `sqrt` returns NaN. The same incremental form loses
about half the mantissa for grazing lines, returning exactly 0 where the true
distance is 1.2e-7.

**Suggested fix.** Clamp negative round-off to 0 (applied in the port), or
reformulate the accumulation. Issue [#421](https://github.com/gradientspaceai/gtengine-js/issues/421).

### `DistLine3OrientedBox3.h`, `DistRay3OrientedBox3.h`, `DistSegment3OrientedBox3.h`

Upstream writes the world-space line point into `closest[0]` before the loop that
reads `closest[i]` as box-frame coordinates, so it is transformed a second time.
Box centre `(10,0,0)`, identity axes, extent `(1,1,1)`, line
`(10,5,0) + t*(0,0,1)` gives `(20,5,1)` instead of `(10,5,1)`. The defect
propagates from the line query to the ray and segment queries. A dead assignment
to `result.closest[0]` in the same function (overwritten by the rotated canonical
point, equal values) is a separate cosmetic item.

Issues [#421](https://github.com/gradientspaceai/gtengine-js/issues/421), [#187](https://github.com/gradientspaceai/gtengine-js/issues/187).

### `DistLineLine.h`, `DistLineRay.h`, `DistLineSegment.h`, `DistRayRay.h`, `DistRaySegment.h`, `DistSegmentSegment.h`

**1. `DistRaySegment` never clamps the ray parameter to `s0 >= 0` (result-corrupting).**
This affects the parallel branch (both cases) and regions 1 and 5 of the
nonparallel branch: reaching those regions only guarantees that the
unconstrained line/line `s0` is nonnegative, and since
`s0(t) = -(a01*t + b0)/a00` is linear in `t`, clamping `t` to 0 or 1 can drive
`s0` negative. Upstream then returns a "closest ray point" behind the ray origin
and an under-reported distance.

Witnesses: parallel, ray `(0,0,0) + t(1,0,0)` versus segment `(-5,1,0)` to
`(-3,1,0)`, upstream distance 1, truth `sqrt(10)`. Nonparallel, well conditioned
(determinant about 2.79, region 1): ray origin
`(-1.9852731227874756, -2.484701693058014, -1.566794514656067)`, direction
`(-2.7053112387657166, -0.4477686882019043, -1.6644186973571777)`, segment
`(1.403347671031952, -2.980962038040161, 2.2564467592164874)` to
`(-1.6737277507781982, -2.9632678627967834, 0.3852156177163124)`, upstream
`s0 = -0.377` and `sqrDistance` 2.675 versus the true 4.136. `DistRayRay.h`
clamps its analogous branches, which supports this being an oversight. Port:
fixed by convexity (if `s0 < 0` the minimum lies on the face `s0 = 0`, where
`s1 = clamp(-b1/a11, 0, 1)`), verified over 6000 randomized configurations
against a brute-force reference.

**2. `DistSegmentSegment::ComputeIntersection`'s `1/2` fallback picks the wrong edge (result-corrupting).**
An out-of-range `f/b` ratio is replaced by `1/2` on the premise that numerator
and denominator are then both nearly zero. The premise is false whenever the line
`dR/ds = 0` passes through a domain corner: the ratio is `1 + 1 ulp` with both
values about 48.8, and `ComputeRobust` then searches the wrong edge. Witness:
`segment0 = <(-0.8119320124387741,-3,0), (7,0,0)>`,
`segment1 = <(0,1.9484716467559338,0), (7,0,6)>`, sine of the angle between
directions 0.79, upstream returns 4.6265 instead of 3.5103. Port: clamp to the
nearest endpoint of `[0,1]`; a 300000-configuration scan went from 12 wrong
minima (worst excess 1.12) to none (worst 2.7e-14). Note that this changes
results for every consumer of `ComputeRobust`.

**3. `det = max(a00*a11 - a01*a01, 0) > 0` does not reliably detect parallelism (result-corrupting, preserved).**
Present in `DistLineLine`, `DistLineRay`, `DistLineSegment`, `DistRayRay`,
`DistRaySegment` and `DistSegmentSegment::operator()`. Witness: ray
`(0,-8,1) + t*(0,-0.2563283918425441,0)` versus segment
`<(0,0,2), (0,-6.4209048338234425,2)>`, both along y, returns 8.0623 instead of
1.8691. In a 170000-configuration degenerate-prone scan each query failed in 136
to 254 cases, worst excess about 16 units; every failing case had exactly
parallel directions. A repair needs a relative tolerance, that is an algorithm
change; upstream's own answer is `ComputeRobust`. `IntrCapsule3Capsule3.h`
inherits the defect: for a capsule strictly containing another the distance comes
back as 0.5375 where the truth is 0, so the query reports no intersection.

**4. Documentation.** `DistLineSegment.h`'s file comment states the segment point
is `closest[0]` and the line point `closest[1]` ("stoed" typo included); the code
does the opposite. `DistSegmentSegment.h`'s `ComputeMinimumParameters` comment
states `H(z) = (end[1][1] - end[1][0]) * dR/dt(...)` while the code correctly
uses `end[1][1] - end[0][1]`. `DistAlignedBoxAlignedBox.h`'s "Any choice of P0 in
closest[0] and any P1 in closest[1]" is false: on an overlap axis both closest
boxes carry the same nondegenerate interval and only coordinate-matched pairs
realize the distance.

Issues [#126](https://github.com/gradientspaceai/gtengine-js/issues/126), [#418](https://github.com/gradientspaceai/gtengine-js/issues/418), [#455](https://github.com/gradientspaceai/gtengine-js/issues/455), [#118](https://github.com/gradientspaceai/gtengine-js/issues/118).

### `DistOrientedBox3Cone3.h`, `Minimize1.h`, `LCPSolver.h`

**1. `Minimize1::GetBracketedMinimum` collapses the bracket at equal-value endpoints (result-corrupting).**
The box-cone query sweeps the quad angle over `[-pi/2, +pi/2]`, whose two
endpoint slices are the same point set, so `F(-pi/2) == F(+pi/2)` to round-off
for every input. When `F` is V-shaped there, `Minimize1` fits a parabola whose
vertex lands a few ulps off the bracket midpoint; the exact `tv == tm` equality
branch is missed, the asymmetric branch collapses the bracket to a degenerate
interval, and the `|denom| <= epsilon` test returns on the next iteration.

**Reproduction.** Frustum `V = (0,0,0)`, `D = (0,0,1)`, angle `pi/6`, `hmin = 1`,
`hmax = 4`, box near `(-1.386, 2.672, 3.307)`: reported distance
`0.362774519113053` versus the true `0.0775455784704954`, after only four
evaluations of `F`. Raising `maxSubdivisions` or shrinking `epsilon` does not
help.

**Suggested fix (applied in the port).** Make the `tv` versus `tm` test
resolution-aware rather than exact, and let the degenerate-parabola case fall
into upstream's own neighbourhood examination. Either change alone leaves inputs
broken. Well-conditioned searches are bit-for-bit unchanged; the consumers
`MinimizeN`, `RemezAlgorithm` and `ApprCone3EllipseAndPoints` are unaffected.

**2. The infinite-cone precondition is never validated.** `Cone` encodes an
infinite cone as `maxHeight = -1`, so an infinite cone silently produces the
infeasible LCP constraint `z3 + z4 <= -1`. Port: asserts finiteness.

**3. The 10-D LCP is degenerate (result-corrupting, preserved).** The
box-quadrilateral Hessian is the Gram matrix of five vectors in R^3, so the LCP
is rank deficient. For about 1 in 400 configurations Lemke's method returns
`HAS_NONTRIVIAL_SOLUTION` with `w != q + M*z` (`w[6] = 0` reported where
`q + Mz` gives -0.71); the query then returns a "box point" outside the box and a
"frustum point" below `hmin`, with a distance below the true minimum.

**4. `LCPSolver<T>::Solve` null-dereferences when constructed with `n <= 0`.**
The dynamic solver leaves `mPoly`, `mQMin`, `mAugmented`, `mVarBasic` and
`mVarNonbasic` null; `Solve`'s size check passes (`0 > 0` is false) and execution
reaches `Copy(mPoly[0], mQMin)`. Port: returns `INVALID_INPUT`, the
classification upstream already uses for unusable input.

**5. A subnormal pivot turns an infeasible problem into a "solution".**
With `M <= 0`, `q < 0` (provably infeasible) a subnormal pivot overflows and
yields `HAS_NONTRIVIAL_SOLUTION` with infinite `z`. Preserved.

**6-10 (minor).** Exact `==` in the minimum-distance assertion of the angle
sweep; a silent `MAX_VALUE` result on LCP failure; a stale "not unit length"
comment in `DistPoint3Tetrahedron3.h`; a wrong trailing comment in
`GetClosestRpmp`; by-value `Z` parameters; a redundant containment test;
non-`const` `j0`/`j1` index arrays; `DistTriangle3CanonicalBox3.h` divides by
zero for a degenerate triangle, where the NaN fallback happens to reach the right
answer.

Issues [#298](https://github.com/gradientspaceai/gtengine-js/issues/298), [#431](https://github.com/gradientspaceai/gtengine-js/issues/431), [#76](https://github.com/gradientspaceai/gtengine-js/issues/76), [#476](https://github.com/gradientspaceai/gtengine-js/issues/476).

### `DistPoint3Cylinder3.h`

**Symptom.** The infinite-cylinder query is unreachable; a cylinder built through
the documented API throws.

**Cause.** The query selects the infinite-cylinder path with

```cpp
if (cylinder.height == std::numeric_limits<T>::max())
```

but `Cylinder3.h` documents and implements the infinite sentinel as
`height = -1` (`MakeInfiniteCylinder` / `IsInfinite`). The finite branch is taken
instead and immediately hits
`LogAssert(cylinder.height > 0, "The cylinder must have a positive height.")`.
Version skew supports this reading: `Cylinder3.h` is `8.0.2025.05.10`,
`DistPoint3Cylinder3.h` is `8.0.2026.08.08`.

**Suggested fix.** Dispatch on `cylinder.IsInfinite()`. Port: fixed.
Issue [#187](https://github.com/gradientspaceai/gtengine-js/issues/187).

### `DistPoint3Frustum3.h`

The LF and UF edge cases do not clamp the free far-edge coordinate; 2 of the 10
assignments are unbounded. With `dMin = 1`, `dMax = 2`, `uBound = rBound = 1`, a
point at `(r,u,d) = (100,3,0)` gives `closest = (2,3,2)`, which is outside the
frustum, at distance 98.0204 instead of the LUF vertex `(2,2,2)` at 98.0255.

Issue [#421](https://github.com/gradientspaceai/gtengine-js/issues/421). Port: fixed.

### `DistPointHyperellipsoid.h`

**Symptom.** For a query point just off the smallest-extent axis, `Bisector`
returns Infinity or a badly wrong distance.

**Cause.** The root satisfies `s > -pSqr[last] = -1` and the lower bracket is
`smin = z[last] - 1`. When `z[last]` (the query point's coordinate along the
smallest-extent axis divided by that extent) is tiny but nonzero, `smin` rounds
to exactly -1, the root lies below the double spacing there, and
`x[i] = pSqr[i]*y[i]/(s + pSqr[i])` cancels catastrophically.

**Reproduction.** Ellipse extents `(3, 2)`, query point `(0, y)`: `y = 0` gives 2
(correct); `y = 1e-16` gives Infinity; `1e-15` gives 1.5012; `1e-14` gives
1.9581; `1e-12` gives 1.99982. The correct distance is 2 throughout.

A double-precision fix means reformulating the bisection variable, so this is
preserved and documented. Issue [#424](https://github.com/gradientspaceai/gtengine-js/issues/424).

### `DistSegment2Circle2.h`, `DistCircle2Circle2.h`, `DistLine3Rectangle3.h`, `DistLine3Triangle3.h`

- `DistSegment2Circle2` zeroes the whole result for a segment strictly inside the
  circle: `numClosestPairs = 0` (the header says 1 or 2) and `distance = 0`
  although the segment does not touch the circle. `DistSegment2Arc2.h` depends on
  the 0 to fall through.
- `DistSegment2Circle2`'s half-open boundary tests (`t1 < one`, `t0 > zero`)
  leave an off-segment parameter when an endpoint lies exactly on the circle
  (observed `[0,4]` and `[-3,1]`).
- `DistSegment2Circle2`'s interval branch chain has measure-zero gaps:
  `t0 == 0 && t1 > 1` and `t0 < 0 && t1 == 1` fall through to the final `else`
  that assumes `0 <= t0 < t1 <= 1`, retaining a pair whose parameter lies outside
  the segment.
- `DistCircle2Circle2::operator()` swaps its arguments when
  `circle1.radius > circle0.radius`, but `DoQuery` writes the point on *its own*
  first argument into `closest[j][0]`. Passing the smaller circle first therefore
  puts the circle1 point in `closest[j][0]`, contradicting the file comment. The
  distance is unaffected.
- `DistLine3Triangle3.h` and `DistLine3Rectangle3.h` test
  `result.sqrDistance == invalid` with `invalid = -1`; a squared distance can
  never equal the sentinel, so the test is dead.

Issues [#421](https://github.com/gradientspaceai/gtengine-js/issues/421), [#187](https://github.com/gradientspaceai/gtengine-js/issues/187), [#118](https://github.com/gradientspaceai/gtengine-js/issues/118). Port: preserved.

### `ETManifoldMesh.h`, `TSManifoldMesh.h`, `ETNonmanifoldMesh.h`, `VEManifoldMesh.h`, `VETNonmanifoldMesh.h`, `VTSManifoldMesh.h`

**1. `VEManifoldMesh::Insert` leaves a rejected edge in the map.** `Insert`
writes the new edge into `mEMap` *before* performing the nonmanifold check.
Under `ThrowOnNonmanifoldInsertion(false)`, the "continue gracefully" path, a
failed insert returns null but leaves the rejected edge in the map, with any
already-processed vertex pointing at it. The mesh is silently corrupted and
subsequent queries see a phantom edge.

**2. `ETManifoldMesh::Insert` / `TSManifoldMesh::Insert` leak the partially built
feature (the same shape one level up).** When `ThrowOnNonmanifoldInsertion(false)`
is active and `Insert` returns `nullptr` mid-loop, edges and faces created
earlier in the same call stay in `mEMap`/`mFMap` referencing a phantom triangle
or tetrahedron that is never added to the feature map, and an existing neighbour
can be left with `T[j]` / `S[j]` pointing at that phantom. On the *throwing*
path the same entries point at `newTri`/`newTetra`, whose `unique_ptr` dies as
the exception unwinds, so they dangle. The header's claim that "the (bad)
triangle will not be part of the mesh" holds for `mTMap`, not for `mEMap`.

**3. `ETManifoldMesh::Insert`'s "nonmanifold" guard actually enforces consistent
orientation.** It checks only `edge->T[0]` (never `T[1]`), compares vertex order
rather than counting adjacent triangles, and never `break`s. A
manifold-but-inconsistently-oriented insertion is therefore rejected with a
misleading "nonmanifold" message, which makes `IsOriented()` and
`MakeConsistentChirality()` reachable only after
`ThrowOnNonmanifoldInsertion(false)`.

**4. `ETManifoldMesh::GetBoundaryPolygon` uses `boundaryEdges[vEdge]`**
(`std::map::operator[]`), which inserts a null-triangle entry for a non-boundary
key while `GetBoundaryPolygons` is iterating that same map: a null dereference on
a malformed mesh.

**5. `ETNonmanifoldMesh::Remove` throws and corrupts the mesh after a
degenerate-triangle `Insert`.** `Insert` accepts a triangle with repeated
vertices. For `(v0, v0, v1)` the edge keys `<v1,v0>` and `<v0,v1>` collide, so
`tri->E[1]` and `tri->E[2]` are the *same* `Edge` object and the second
`edge->T.insert(tri)` is a no-op. `Remove` then iterates `i = 0,1,2` with

```cpp
size_t numRemoved = edge->T.erase(tri);
LogAssert(numRemoved > 0, "Unexpected condition.");
```

The third iteration sees `numRemoved == 0` and throws, after the shared edge has
already been erased from `mEMap` and with the triangle still present in `mTMap`.
Either a distinct-vertex guard in `Insert` or a tolerant `Remove` would fix it.

**6. `VETNonmanifoldMesh::Remove` has an inverted assertion.**
`LogAssert(VAdjacent.size() != 0 || EAdjacent.size() != 0, "Malformed mesh.")`
fires for every well-formed mesh: by the time it runs both sets have just been
emptied, so even `Insert(0,1,2); Remove(0,1,2)` throws. Port: asserts both are
empty.

**7. `VTSManifoldMesh::Remove` over-erases `VAdjacent` (result-corrupting).**
On removing a tetrahedron, upstream erases the other two vertices of every
destroyed face from each vertex's `VAdjacent`. Unlike the edge-based VET mesh a
vertex pair can lie on several faces, so adjacencies still contributed by
surviving faces get dropped. With tetrahedra `<0,1,2,3>` and `<0,1,2,4>` sharing
face `<0,1,2>`, removing `<0,1,2,3>` leaves vertex 0 with `VAdjacent = {4}`
instead of `{1,2,4}`. The end-of-loop assertion does not catch this. Port:
fixed, by restoring the insert-time invariant after the destroyed faces are
removed.

Issues [#73](https://github.com/gradientspaceai/gtengine-js/issues/73), [#179](https://github.com/gradientspaceai/gtengine-js/issues/179), [#212](https://github.com/gradientspaceai/gtengine-js/issues/212), [#240](https://github.com/gradientspaceai/gtengine-js/issues/240), [#256](https://github.com/gradientspaceai/gtengine-js/issues/256), [#472](https://github.com/gradientspaceai/gtengine-js/issues/472).

### `ExtremalQuery3BSP.h`

**Symptom.** The BSP tree answers some extremal-vertex queries with the wrong
vertex.

**Cause.** During construction, `InsertArc` descends the tree comparing only the
arc's two endpoints against each visited node's great circle, never against the
accumulated spherical region that node governs. A node can therefore end up
governing a region its arc does not touch.

**Reproduction.** Measured against brute-force max-dot over 5000 random
directions: tetrahedron 0 wrong, octahedron 0 wrong (also 0 of 6000 over 200
random rotations each), regular icosahedron 49 of 5000 (0.98%), subdivided
icosahedron 188 of 5000 (3.76%). Sharpest statement: the 42-vertex subdivided
icosahedron answers 2 of its own 42 vertex directions with the wrong vertex. One
failure was verified geometrically: the leaf's spherical region lies entirely
inside vertex 9's Gauss cell while the tree answers vertex 4.

A prototype repair that clips arcs against the accumulated region during descent
fixed 42-vertex polytopes but not the icosahedron, so a correct implementation
appears to need a redesigned construction. Preserved in the port with a pointer
to `ExtremalQuery3PRJ` for exact results.

Also: `GetTreeDepth` reports the DFS stack high-water mark, not the tree depth;
and the class never checks `VETManifoldMesh::Insert` failure or a boundary edge's
null `Edge::T[1]`, so `std::map::operator[]` silently maps null to a zero face
normal and the dot-sign classification is garbage for nonmanifold or open input.

Issue [#290](https://github.com/gradientspaceai/gtengine-js/issues/290).

### `FastMarch.h`, `FastMarch2.h`, `FastMarch3.h`

**1. `FastMarch3::Initialize` omits the six boundary faces (result-corrupting).**
Lines 232-341 mark only the vertices and edges of the grid boundary as zero
speed, while the 2D sibling marks its entire border and the comment says the
whole boundary is meant to be excluded. Face voxels stay "far", get pulled into
the front, and then index off-grid (`i - mXYBound` for a `z = 0` voxel wraps in
`size_t`). Port: fixed, all six faces marked.

**2. `ComputeTime`'s negative-discriminant fallback takes the LARGER neighbour time.**
`mTimes[i] = (diff >= 0 ? xConst : yConst)` is the maximum; the Godunov fallback
is `min + 1/speed`. Combined with `IsValid` counting trial pixels as upwind data,
a recomputation can raise an already-correct time. Reproduction: a 6x4 grid,
speed 1/2, seeds `(4,2)` and `(1,2)`; pixel `(2,2)`, a 4-neighbour of a seed,
starts at the correct 2 and is raised to `2*(1 + sqrt(2)/2) = 3.414` when `(3,2)`
is accepted. Consequence: the computed times are not the discrete-eikonal
solution, are not monotone in the seed set, and break the L1 bound. Preserved.

**3. Minor.** `FastMarch2/3` store grid spacings that the numerical method never
uses. `FastMarch::GetTimeExtremes` redundantly re-tests one element.

Issues [#121](https://github.com/gradientspaceai/gtengine-js/issues/121), [#439](https://github.com/gradientspaceai/gtengine-js/issues/439), [#52](https://github.com/gradientspaceai/gtengine-js/issues/52).

### `FPInterval.h`

**1. `ProductLowerBound` / `ProductUpperBound` are wrong in the both-straddle-zero branch.**
The branch returns `u[0]*v[0]` as the lower bound and `u[1]*v[1]` as the upper
bound. `u[0]*v[0]` is positive (negative times negative), so it cannot bound a
product set containing negatives. The correct values are
`min(u[0]*v[1], u[1]*v[0])` and `max(u[0]*v[0], u[1]*v[1])`, which is what
`Mul2`/`operator*` computes. Evidence: `u = v = [-1,1]` has true range `[-1,1]`
but `ProductLowerBound` returns `+1`; `u = v = [-4,1]` has true range `[-4,16]`
but `ProductUpperBound` returns `1`. Every other branch agrees with `operator*`.
Nothing in GTE calls either function.

**2. NaN endpoints lose the enclosure.** `FPInterval(0,0) / FPInterval(MIN_VALUE, 1)`
and `(MAX + MAX) * 0` both give NaN bounds.

Note: `FPInterval.h` does *not* share `SWInterval.h`'s overflow defect, because it
uses `fesetround` and directed rounding leaves an infinity alone.

Issues [#75](https://github.com/gradientspaceai/gtengine-js/issues/75), [#476](https://github.com/gradientspaceai/gtengine-js/issues/476). Port: preserved.

### `GaussianElimination.h`

For a matrix whose entries are all of denormal magnitude, `Inverse()` returns NaN
entries while reporting `invertible = true`, and `Determinant()` returns NaN. The
pivot magnitude comparisons succeed (the pivot is nonzero) but the reciprocal
`1 / pivot` overflows to infinity and subsequent row operations produce
`inf * 0 = NaN`. A robust implementation would either report non-invertible when
`1/pivot` is not finite, or scale rows before elimination.

The port preserves this deliberately: full pivoting takes pivots in decreasing
magnitude, so a subnormal pivot does not imply an underflowing determinant.
`diag(1e300, 1e-320)` has determinant `9.99988867182683e-21`, which is returned
correctly today although its inverse overflows; reporting non-invertible would
replace a correct, representable determinant with 0 in every dependent query.

Issue [#375](https://github.com/gradientspaceai/gtengine-js/issues/375).

### `GaussNewtonMinimizer.h`, `LevenbergMarquardtMinimizer.h`, `LinearSystem.h`

**1. `LevenbergMarquardtMinimizer::DoIteration` builds the step from a stale residual (result-corrupting).**
`ComputeLinearSystemInputs` recomputes `mJ` at `pCurrent` but forms
`mNegJTF = -(mF * mJ)` with `mF` still holding `F` evaluated at the *previous
rejected candidate* `pNext`. Any repeat of the lambda-adjustment loop for the
same `pCurrent`, precisely the mechanism that distinguishes LM from
Gauss-Newton, therefore solves a system with mismatched `J` and `F`.

**Reproduction.** `E(x) = atan(x)^2` from `x0 = 2`: the J-path gives first
iterate `1.3693536108807434` while the equivalent J-plus path (immune, since it
derives everything from `pCurrent`) gives `-0.7678717944852251`, and the
trajectories diverge from there. `GaussNewtonMinimizer.h` is not affected.
Port: fixed by re-evaluating `F(pCurrent)` when stale; the two paths then agree
to 12 or more digits.

**2. `LinearSystem::SolveTridiagonal` / `SolveConstantTridiagonal` never validate `N`.**
They read `diagonal[0]`/`B[0]` unconditionally and allocate
`std::vector<Real> tmp(static_cast<size_t>(N) - 1)`, so `N = 0` wraps the size to
`SIZE_MAX`. Port: asserts `N >= 1`.

**3. `SolveSymmetricCG` with `B = 0`** computes `alpha = 0/0 = NaN` on the first
iteration and overwrites the already-exact solution `X = 0` with NaN. It also
recomputes the loop-invariant `Dot(N, B, B)` every iteration.

**4-6 (minor, preserved).** `Result::numIterations` is the raw loop counter, so
an exhausted loop reports `maxIterations + 1` and `maxIterations = 0` reports 1.
The LM "fall back to a Gauss-Newton iterate" branch actually uses the maximally
inflated lambda (not `lambda = 0`) and never divides it back down, so the damping
persists into later outer iterations; the comment and code disagree.
`GaussNewtonMinimizer`'s unconditional `pCurrent = pNext` accepts
error-increasing steps, with the side effect that `converged` is usually `false`
for exactly linear models.

Issue [#261](https://github.com/gradientspaceai/gtengine-js/issues/261).

### `GenerateMeshUV.h`, `HelmertTransformation7.h`

**1. `AssignBoundaryTextureCoordinatesSquare` selects square corners without tolerance (result-corrupting).**
It locates the four corners with `std::lower_bound(..., 0.25/0.50/0.75)` on
normalized accumulated arc length and then continues with `for (++i; ...)`,
assuming the corner consumed exactly that index. When a boundary vertex sits
mathematically on a quarter mark but its accumulated arc length rounds just below
it, the corner lands on the following vertex and one boundary vertex is displaced
by a full edge.

**Reproduction.** A uniform 6x6 grid on `[0,1]^2`, where the total arc length
accumulates to `4.000000000000001`: the boundary vertex at `(1, 1/5)` receives uv
`(1, 0)`. The same symptom was pinned on a 4x4 planar grid, where the corner arc
length is one ulp below 0.25; translating the grid by 4 in y makes it go away.
Preserved, because the map remains a monotone parametrization of the square
boundary and an epsilon would change the output for every mesh.

**2. Nondeterminism.** `mInteriorEdges` is a `std::set<Edge*>` iterated in
pointer order, so the vertex-weight assembly order varies run to run. Results are
order-independent in exact arithmetic, but the floating-point sums differ. The
port iterates in mesh edge-key order.

**3. No boundary validation.** A closed surface leaves
`mBoundaryStart == INT32_MAX` and indexes out of bounds.

**4-8 (minor).** `!= -1.0f` compares a float literal against a `(Real)-1`
sentinel in template code; a stale comment claims the first disk boundary vertex
is `(1/2, 0)` while the code assigns `(1, 1/2)`; a missing `break` in the weight
write-back loop causes a harmless extra scan; the multithreaded path divides by
`mNumThreads` unguarded.

**9-11. `HelmertTransformation7.h`.** A full 3x3 matrix `A = sum(u_i v_i^T)` is
computed and never read, a dead SVD/Procrustes leftover; `A(r, c) = 0.0` uses a
non-generic `double` literal in template code and `UpdateF` has a stray `;;`;
`LogAssert(mNumPoints >= 7)` appears to confuse "7 parameters" with "7 points",
since a 7-parameter similarity fit needs only 3 non-collinear correspondences.

`PolyhedralMassProperties.h` was reviewed on the same pass with no suspects; it
reproduces the analytic cube, box and tetrahedron values exactly.

Issue [#262](https://github.com/gradientspaceai/gtengine-js/issues/262).

### `GMatrix.h`, `SortPointsOnCircle.h`

Both provide comparators that are not strict weak orderings, which is undefined
behaviour for the standard containers and algorithms that consume them:

- `GMatrix`'s comparison operators return `false` for all of `<`, `<=`, `>`, `>=`
  when dimensions differ, while `!=` is `true`. Mixed-dimension `GMatrix` keys in
  `std::map`/`std::set` are therefore UB. Harmless in GTE's own usage.
- `SortPointsOnCircle::LessThanByGeometry` is not a strict weak ordering when a
  point coincides with the sort centre. With `W = (0,0)` the comparator says `W`
  is "not less" than every point with `y >= 0` and "greater" than every point with
  `y < 0`, so `W` is equivalent to points that are themselves strictly ordered.

Issues [#88](https://github.com/gradientspaceai/gtengine-js/issues/88), [#394](https://github.com/gradientspaceai/gtengine-js/issues/394). Port: preserved.

### `GradientAnisotropic2.h`, `GradientAnisotropic3.h`

`ComputeParameter` feeds padded coordinates to the unpadded `GetUx`/`GetUy`/`GetUz`:
the average is taken over a window shifted by one pixel, and the final sample
reads one element past the padded buffer. Verified against a hand-computed
average (1/0.625 for the `u = x` ramp). Port: iterates `0 <= x < mXBound`.

On a constant image the gradient average is zero, so `mParameter = Infinity` and
`exp(-Infinity*0) = NaN`. No guard exists upstream; none was added.

Minor: `OnPreUpdate` never calls the base-class `OnPreUpdate`, so a masked filter
with Neumann conditions never refreshes its mask border; the divergence term is
divided by the spacing once rather than twice, which is dimensionally consistent
only at unit spacing.

Issues [#122](https://github.com/gradientspaceai/gtengine-js/issues/122), [#439](https://github.com/gradientspaceai/gtengine-js/issues/439).

### `Histogram.h`

**1. The maximum sample misses the last bucket about 7% of the time.**
`fl((B-1)/d) * d` can round below `B-1`, so `static_cast<int32_t>` truncates the
maximum sample into bucket `B-2`. Measured 140053 of 2000000 (7.0%) over random
`(B, min, max)` triples. This is an IEEE rounding property, not a translation
error.

**2. A subnormal sample range writes outside `mBuckets`.** `mult` overflows to
infinity, `static_cast<int32_t>(inf)` is undefined behaviour, and
`++mBuckets[index]` is unbounded on the rescaled path.

Also noted: `FastGaussianBlur2.h` clears `mInput`/`mOutput` on exit while
`FastGaussianBlur3.h` leaves them dangling.

Issue [#436](https://github.com/gradientspaceai/gtengine-js/issues/436). Port: preserved.

### `Hyperplane.h`, `Hyperellipsoid.h`

**1. `ComputeFromPoints` passes `-1` as the SVD multiplier.**
`SingularValueDecomposition::Solve`'s second parameter is `Real multiplier`
(default 8), and `Solve` begins with
`LogAssert(input != nullptr && multiplier > zero, ...)`, so `svd.Solve(&edge[0], -1)`
throws for every input. The `-1` looks like a stale `sortType` argument from an
older signature. Port: uses the default multiplier.

**2. `N = 2` constructs `SingularValueDecomposition(2, 1, 32)`,** which trips that
constructor's `LogAssert(mNumCols >= 2 && ...)`. So `Hyperplane<2,T>` cannot be
constructed from points at all upstream, even with item 1 fixed. Port: routes
`N = 2` through GTE's own `ComputeOrthogonalComplement` for `Vector2`
(`v[1] = -Perp(v[0])`, orthonormalized), which is the one-dimensional orthogonal
complement the SVD path is meant to compute; pinned across `N = 2..5`.

**3. `Hyperellipsoid::ToCoefficients`** divides `coeff[i] /= maxValue` with no
zero check and computes a `maxIndex` that is never read. `maxValue == 0` requires
every extent to be zero and yields NaN or Inf rather than a crash.

Issue [#217](https://github.com/gradientspaceai/gtengine-js/issues/217).

### `IEEEBinary16.h`

**Symptom.** `Convert32To16` turns some NaNs into infinity, violating IEEE-754's
requirement that a NaN convert to a NaN.

**Cause.**

```cpp
uint16_t maskPayload = static_cast<uint16_t>(trailing32 >> CONVERSION_TRAILING_SHIFT);
return sign16 | F16::EXPONENT_MASK | maskPayload;
```

`CONVERSION_TRAILING_SHIFT` is 13, so only the high 9 of the 23 trailing bits
survive. Any NaN with payload below `2^13`, for example `0x7F800001`, yields
`maskPayload == 0`, and `sign16 | 0x7C00 | 0` is the encoding of INFINITY. Only
signaling NaNs with tiny payloads are affected; quiet NaNs always carry
`0x00400000`, which maps to the 16-bit quiet bit.

Issue [#110](https://github.com/gradientspaceai/gtengine-js/issues/110). Port: preserved, pinned by an exact property over all 65536 encodings.

### `Image2.h`, `Image3.h`

In `Image3.h` (~L690-760) and `Image2.h` (L482-515) the absolute-coordinate
neighborhood getters compute

```cpp
nbr[i][0] = static_cast<size_t>(x) + inbr[i][0];
```

where `inbr[i][0]` is a negative `int32_t` for the minus-direction offsets. For a
boundary pixel or voxel (say `x == 0`) the addition wraps to `SIZE_MAX`,
contradicting the adjacent documentation, which says no clamping is used when the
point is on the boundary. That wording implies the caller receives out-of-range
but meaningful coordinates, not wrapped ones. The port deviates deliberately:
coordinates are plain signed numbers, so out-of-bounds neighbours come back as
negative values as the documentation intends.

Issue [#64](https://github.com/gradientspaceai/gtengine-js/issues/64).

### `ImageUtility2.h`, `ImageUtility3.h`

**1. `ImageUtility3::Dilate` skips `x = 0` sources.** The innermost loop is
`for (std::int32_t i0 = 1; i0 < dim0; ++i0)`, so a foreground voxel at `x = 0`
never acts as a dilation source. The 2D counterpart uses `i0 = 0`. Invisible
under the header's zero-boundary precondition but wrong for general binary input.
Port: fixed.

**2. `ImageUtility3::Close<N>` is uncompilable dead code.** It carries the 2D
`static_assert(N == 4 || N == 8, ...)` and constructs
`Image3<std::int32_t> temp(GetDimension(0), GetDimension(1))`, a two-argument
`Image3` constructor that does not exist. Any instantiation fails to compile.
Port: implemented as `close6`/`close18`/`close26` with correct neighbourhoods.

**3. `ImageUtility2::DrawEllipse` hangs for `xExtent == yExtent == 0`.**
`for (; yExtSqr * x <= xExtSqr * y; ++x)` evaluates `0 <= 0` forever. Port:
visits the degenerate point once.

**4. `ImageUtility2::GetSkeleton` empties solid even-sided squares.** 2x2 and 4x4
blocks skeletonise to nothing, contradicting the header's "connectivity ... is
preserved"; 1x1, 3x3 and 5x5 behave. Preserved.

**5. Minor.** `ImageUtility2::DrawLine` never updates `maxValue` in its
`if (dy > maxValue)` branch (the 3D version does assign it).
`GetComponents`/`ExtractBoundary`/`GetSkeleton` read neighbours with no range
test, relying wholly on the documented zero-boundary precondition.

Issues [#129](https://github.com/gradientspaceai/gtengine-js/issues/129), [#443](https://github.com/gradientspaceai/gtengine-js/issues/443).

### `IncrementalDelaunay2.h`

- `GetHull` dereferences `edges.begin()` on an empty edge map.
- `GetHull`'s unbounded `while (vNext != vStart)` walk loops forever, or writes
  past the output, for collinear input: the supervertex triangles contribute a
  boundary path that ends in a 2-cycle which never returns to the start.
  Reproduced; the port uses a bounded walk.
- `RetriangulateBoundaryRemovalPolygon`'s `numPolygon == 2` branch tests
  `polygon[0] == vRemovalIndex`, which can never hold, because `polygon` contains
  the *opposite* vertices. Instrumented: unreachable in practice, since the
  `mVertexIndexMap.size() == 4` path catches the case first.
- Comments claim `GetNumVertices`/`GetVertices` exclude the three supervertices;
  they do not.
- `DoEarClipping` skips `RPPolygon::Remove` on its early exit, leaving
  `GetNumActive()` one too large. Harmless on every reachable path.

Issue [#290](https://github.com/gradientspaceai/gtengine-js/issues/290).

### `InscribedFixedAspectRectInQuad.h`

`Execute` fails on solvable convex quadrilaterals in two ways. The LP optimum was
verified against an independent 16-constraint vertex enumeration.

**(a) Shared case index.** Constraint normals that share a case index `j` lie in
a common plane of `(u, v, w)` space, so when three edge normals fall in the same
quadrant the `alpha` assertion fires. Smallest lattice example: the quad
`<(0,3), (0,2), (1,1), (3,0)>`.

**(b) Degenerate feasible interval.** When the optimum touches all four edges the
feasible interval collapses to a point, and the untoleranced endpoint comparison
can order it empty in both pairings, giving "Unexpected interval intersection
type". Reproduced with a quad whose optimum width is about 0.6565.

The header comment also writes `quad[(i + 1) % 3]` for a quadrilateral (should be
`% 4`).

Issues [#395](https://github.com/gradientspaceai/gtengine-js/issues/395), [#101](https://github.com/gradientspaceai/gtengine-js/issues/101). Port: preserved.

### `IntpAkimaUniform2.h`, `IntpAkimaUniform3.h`

**Symptom.** Upstream does not reproduce a bilinear function in the boundary
cells: for `f = x*y` on a 3x3 grid, `evaluate(1.5, 0.75)` returns `1.1015625`
where the exact value is `1.125`, a first-order error.

**Cause.** `GetFXY` reuses the min-boundary one-sided difference coefficients at
the max boundaries without negating the sign. The high-index boundary stencil
runs on the reversed index sequence `(n-1, n-2, n-3)` with the forward
coefficients `(-3, 4, -1)/(2h)`, which computes `-dF/dx`. Recovered `FXY` is `+1`
everywhere except the last row and column, where it is `-1`.

The same slip exists in `IntpAkimaUniform3.h`: `GetFXY`, `GetFXZ`, `GetFYZ` and
`GetFXYZ` build the min-boundary masks from outer products of
`ODer = {-1.5, 2, -0.5}` and reuse the identical coefficients at max boundaries
with reflected sample indices. Reflection reverses the differentiation direction,
so the mask needs a sign flip per reflected axis and never gets one: the mixed
second derivatives are negated on single max faces (double-max corners cancel
back to correct), and `FXYZ` is negated wherever an odd number of max boundaries
is involved.

There is no 1D analogue: `IntpAkimaUniform1` and `IntpAkimaNonuniform1` use the
textbook Akima slope extrapolation with correct indices and signs at both ends.

Issue [#58](https://github.com/gradientspaceai/gtengine-js/issues/58). Port: fixed in both files; tests pin exact bilinear and trilinear reproduction.

### `IntpBicubic2.h`, `IntpBilinear2.h`, `IntpTricubic3.h`, `IntpTrilinear3.h`

**Symptom.** The documented input clamping does not happen; the interpolants
extrapolate outside the domain.

**Cause.** All four headers say "The functions clamp the inputs to
`xmin <= x <= xmax` ...", but only the cell index `ix` is clamped; the fractional
coordinate `xIndex - ix` keeps its out-of-range value.

**Consequences.** Trilinear and bilinear extrapolate linearly below the minimum;
above the maximum they accidentally hold the boundary value, because the
collapsed stencil's blend weights sum to 1. Tricubic and bicubic extrapolate an
unbounded cubic on both sides, and far outside the domain even a constant field
is destroyed by round-off (weights reach about 1e7 and cancellation fails).
Reproduction: for `f = x + 2y` on `[0,2]^2`, `IntpBilinear2(-3, 1)` returns `-1`
rather than the clamped `2`.

`IntpAkima1.h` and `IntpAkimaUniform2.h` genuinely clamp, so this reads as
copy-pasted stale comment text.

Cosmetic in the same family: `IntpTricubic3.h` and `IntpBicubic2.h` have a stray
double semicolon in the `mBlend[0][3]` initializer; `IntpBilinear2.h`'s
constructor comment claims a 3x3 block minimum while the assertion and algorithm
need only 2x2.

Issue [#69](https://github.com/gradientspaceai/gtengine-js/issues/69). Port: preserved.

### `IntpBSplineUniform.h`, `IntpLinearNonuniform2.h`, `IntpLinearNonuniform3.h`

- `ComputePowers` overruns its buffer for degree 0: `resize(degree+1)` is
  followed by an unconditional `powerDSDT[1] = dsdt;`, and `GetKey` later reads
  `[1]`. Degree 0 is clearly meant to work, since `ComputeBlendingMatrix`
  special-cases it.
- Latent out-of-range read extracting `A`: `Q[k][col]` runs `col = 0..degree`,
  but `Q[k]` comes from `Polynomial1` arithmetic that calls
  `EliminateLeadingZeros`.
- Dead local `Polynomial1<Real> sm1 = {-1, 1};` in `ComputeBlendingMatrix`, with
  a comment claiming it is used; the code uses `GetTranslation(1)`.
- `IntpLinearNonuniform2/3` discard `GetIndices`' `bool` return, so on failure the
  value-initialized `{0,0,0}` tuple is used and `true` is returned to the caller.

Issue [#135](https://github.com/gradientspaceai/gtengine-js/issues/135).

### `IntpQuadraticNonuniform2.h`, `IntpVectorField2.h`

**1. Every mesh-accessor failure flag is discarded (result-corrupting).**
Each `GetVertices` / `GetIndices` / `GetAdjacencies` / `GetBarycentrics` call
ignores the returned `bool` and computes with value-initialized zero arrays on
failure. Zero index triples then index sample 0, silently blending unrelated
samples. Port: fixed (failed triangles are skipped and the result is marked
invalid).

**2. The "closest subtriangle" fallback consults `bary` even when
`ComputeBarycentrics` returned false.** The value-initialized `(0,0,0)` is a
corner of the `[0,1]^3` cube, so a degenerate subtriangle always wins the
fallback and yields `F = 0` with a `1/0` gradient. Preserved, since skipping
invalid candidates would leave the selector outside its intended 1..6 range.

**3. `IntpVectorField2::operator()`'s `&&` short-circuit leaves a half-written output.**
A failed x-interpolation skips the y-query, so `output[1]` retains the caller's
stale value while `output[0]` may already be overwritten. Port: fixed.

**4-5 (minor).** `ProcessTriangles` ignores `Inscribe`'s failure return, so a
degenerate triangle silently receives centre `(0, 0)` (Delaunay input contains no
degenerate triangles). `ComputeCoefficients` and `ComputeCrossEdgeIntersections`
divide by `alpha`/`beta`/`gamma` and `invDet` unguarded; `ComputeCoefficients`'
comment says "circumscribing circle" where the code uses the inscribed-circle
centre.

Issue [#337](https://github.com/gradientspaceai/gtengine-js/issues/337).

### `IntpThinPlateSpline2.h`, `IntpThinPlateSpline3.h`

**1. The WARNING comment's invariance claim is wrong in both directions**
(lines 14-15: the interpolation "is invariant to translations and rotations of
(x,y) but not to scaling").

- The *classical* thin-plate spline **is** invariant to uniform scaling. In 3D
  `Kernel(t) = -|t|` is homogeneous of degree 1. In 2D
  `Kernel(s*r) = s^2*Kernel(r) + 2*s^2*log(s)*r^2`, and the extra `r^2` term
  contributes only a constant under the side conditions
  `sum(a_i) = sum(a_i x_i) = sum(a_i y_i) = 0`, absorbed by the affine part.
  Verified numerically to about 1e-15 for `s` in `{0.05, 4, 1000}`.
- The *transformed* spline, the one the WARNING actually describes, with the
  per-axis unit-square or unit-cube rescaling, is **not** rotation-invariant,
  because independent per-axis scaling is not a similarity. Verified with a
  non-square bounding box.

**2. Unguarded division by a zero coordinate range in the domain transform.**
A flat sample set (all `y` equal, say) gives `mYInvRange = Infinity`, so
`0 * Infinity = NaN` for every transformed coordinate. The NaN propagates into
`A`, the inverse reports non-invertible, and `IsInitialized()` is `false`: the
right outcome, but by NaN accident rather than a check.

**3. `ComputeFunctional` is discontinuous in the smoothing parameter.** It
returns `a^T A a` when `lambda == 0` but `lambda * (a^T A a)` when `lambda > 0`,
so the value is not comparable across smoothing settings and jumps as
`lambda -> 0+`.

Recorded but not filed as a defect: singular sample configurations are only
detected when the floating-point inverse fails, so a coplanar 3D sample set passes
`Inverse(Q)` on the rounded matrix and the spline reports initialized with
ill-determined affine coefficients. There is no rank test.

Issue [#191](https://github.com/gradientspaceai/gtengine-js/issues/191). Port: preserved.
### `IntrAlignedBox3Cone3.h`, `IntrOrientedBox3Cone3.h`

**Symptom.** A reused query object reports false negatives.

**Cause.** `BoxFullyInConeSlab` copies 12 edges into `mCandidateEdges` and sets
`mNumCandidateEdges = 12` without clearing `mAdjacencyMatrix`, orphaning
adjacency bits at indices `>= 12`. A later clipping query on the same object then
has `InsertEdge` silently drop edges.

**Reproduction.** A 4000-case randomized reused-versus-fresh comparison shows 89
disagreements with the upstream logic and 0 once the adjacency matrix is cleared.
Port: fixed.

Both headers also still carry upstream's own TODO noting they were written for
the old `maxHeight = max()` infinite-cone representation and never rewritten for
the redesigned `Cone` class.

Issues [#301](https://github.com/gradientspaceai/gtengine-js/issues/301), [#334](https://github.com/gradientspaceai/gtengine-js/issues/334).

### `IntrAlignedBox3Sphere3.h`, `IntrAlignedBox2Circle2.h`

**1. `DoQueryRayRoundedFace` accepts the first probe (result-corrupting).**
It accepts the first rounded-edge probe, including its vertex fallback, and never
tries the other edge. Any piece of the Minkowski sum gives only an upper bound on
the contact time, so first-wins selection is unsound. There is also no case for
the rounded edge parallel to the entry-face normal, and `VertexSeparated` probes
only the rounded vertex.

**Reproduction.** Against a Minkowski-sum oracle over 3,000,000 random
configurations, the upstream logic missed 674 contacts and reported 14,689 late,
worst lateness 1.90. Two individual cases were pinned first: a contact reported
0.11 late with 0.014 penetration, and a contact missed entirely with penetration
5.7e-3 at closest approach (observed through `IntrOrientedBox3Sphere3`). Port:
fixed; zero misses, zero late, zero spurious, and every one of the 2,784,307
cases where the old time was already bit-identical to the oracle is unchanged.

**2. `FIQuery::operator()` reports `contactPoint = boxCenter` for a no-contact result.**
`result.contactPoint += boxCenter` is applied unconditionally, so a no-contact
result returns the box centre instead of the documented `(0,0,0)`. The 2D sibling
`IntrAlignedBox2Circle2.h` guards the same translation with
`if (result.intersectionType != 0)`. Preserved, since the field is documented
invalid when there is no contact.

Issues [#458](https://github.com/gradientspaceai/gtengine-js/issues/458), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465), [#250](https://github.com/gradientspaceai/gtengine-js/issues/250).

### `IntrAlignedBox3OrientedBox3.h`, `IntrOrientedBox3OrientedBox3.h`

The two files implement the same separating-axis algorithm but disagree on the
parallel-pair cutoff: `absDot01[i][j] >= cutoff` in the aligned file and
`> cutoff` in the oriented one. At the default `epsilon = 0` (cutoff exactly 1) a
pair whose rounded `|Dot|` lands on `1.0` is treated as parallel by one and not
by the other; the oriented version then evaluates three nearly-zero cross-product
axes, exactly the case the epsilon mechanism exists to avoid. No result-changing
case was produced in testing.

Issue [#450](https://github.com/gradientspaceai/gtengine-js/issues/450). Port: both spellings preserved.

### `IntrAreaEllipse2Ellipse2.h`

**1. Five member variables are never initialized (result-corrupting).**
`T mZero, mOne, mTwo, mPi, mTwoPi;` are declared with no constructor and never
assigned anywhere in the header, so every computed area, every `pi` scaling and
the `dtheta <= mPi` branch read indeterminate values. Port: fixed with literals
and pi constants.

**2. The documented "axes are not required to be normalized" is false.**
The matrix `M` is length-invariant, but the polar angles
`atan2(Dot(axis[1], X), Dot(axis[0], X))`, the point ordering in `Area4`, and the
`pi * extent[0] * extent[1]` full-ellipse area are not. (The FI standard form
divides by `|U|^2`, so the claim holds only for a common axis scale.) Port:
normalizes private axis copies, with no behaviour change for unit axes.

Issues [#301](https://github.com/gradientspaceai/gtengine-js/issues/301), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465).

### `IntrCanonicalBox3Cylinder3.h`, `IntrHalfspace3Cylinder3.h`, `IntrCylinder3Cylinder3.h`, `IntrLine3Cylinder3.h`, `IntrRay3Cylinder3.h`, `IntrSegment3Cylinder3.h`, `IntrTriangle3Cylinder3.h`

**1. `IntrCanonicalBox3Cylinder3::DoQueryNoZeros`, the `(U1, -D)` block, has a sign typo (false negatives).**

```cpp
P1 = { +E[0] - s3p / D[0], -E[1], -E[2] };
```

`s3p` is derived from vertex `V3 = (+E[0], +E[1], -E[2])`, so the `y` component
should be `+E[1]`. The wrong sign puts the segment endpoint on the wrong box
edge, inflating the projected distance, so intersecting configurations are
reported as disjoint. Verified empirically: restoring the upstream sign fails
three of the port's tests, including a brute-force case with a sampled cylinder
point provably inside the box. Port: fixed.

**2. `IntrHalfspace3Cylinder3`'s `max((T)1, ...)` clamp makes `root` unconditionally 1 (false positives).**

```cpp
T root = std::sqrt(std::max((T)1, (T)1 - absNdW * absNdW));
```

`max(1, 1 - x^2)` is always 1, so the term the comment documents as
`sqrt(1 - Dot(N,W)^2)` is never computed; the clamp should be `max((T)0, ...)`.
Cylinders tilted toward the plane normal are reported as intersecting when they
are not. Port: fixed.

**3. No `IsFinite()` check for the `height = -1` infinite sentinel.**
`IntrHalfspace3Cylinder3`, `IntrCylinder3Cylinder3`, `IntrRay3Cylinder3`,
`IntrSegment3Cylinder3` and `IntrTriangle3Cylinder3` all treat an infinite
cylinder as finite; `IntrTriangle3Cylinder3` in particular gets an empty slab
(`-h/2 = 0.5 > -0.5 = h/2`) and silently returns garbage. `IntrLine3Cylinder3` has
no infinite-cylinder branch at all, so an infinite cylinder silently yields no
intersection. Port: asserts finiteness using upstream's own message from
`IntrCanonicalBox3Cylinder3.h`.

**4. `IntrSegment3Cylinder3` has no zero-length-segment guard.** `GetCenteredForm`
normalizes `p[1] - p[0]`, `Normalize` leaves the zero vector zero, and the
`FIQuery<Line3, Cylinder3>::DoQuery` it delegates to then divides by a zero
quadratic coefficient; a zero-length segment inside the cylinder returns
`intersect = true` with NaN parameters and points (outside is rejected cleanly).
The same unit-direction assumption sits in `IntrLine3Cylinder3.h`. Preserved.

**5. `IntrTriangle3Cylinder3::DiskOverlapsPolygon` reports containment for every degenerate projected polygon (result-corrupting).**
Lines 451-487 count the signs of `DotPerp(Q[i0], Q[i0] - Q[i1])` and conclude the
origin is inside when `positive == 0 || negative == 0`. When the projected
polygon has collapsed to a point or to a segment through the origin, every
`DotPerp` is zero, both counts are zero, and the query reports an intersection
regardless of distance. This is not confined to degenerate triangles: **any**
nondegenerate triangle whose plane contains the cylinder axis projects to a
segment through the origin and is reported as intersecting however far away it
is. Point-degenerate reproduction: triangle `(4,0,-4), (4,0,4), (4,0,-4)` versus
a cylinder with centre `(0,0,0)`, axis `(0,0,1)`, `r = 1`, `h = 2` returns
`intersect = true`. Port: fixed.

**6. Comments.** `IntrCylinder3Cylinder3`'s constructor comment says
`phi[j] = pi*j/numPhi` while the code uses `GTE_C_HALF_PI` (the code is correct
for hemisphere sampling). `IntrTriangle3Cylinder3` writes "on the tob of the
slab" at lines 178, 186, 213 and 246, and its case 4c/4d comment parameterizes
edge `<V2,V0>` as `V0 + t*(V2-V0)` while the formulas and code use
`V2 + t*(V0-V2)`.

Issues [#197](https://github.com/gradientspaceai/gtengine-js/issues/197), [#206](https://github.com/gradientspaceai/gtengine-js/issues/206), [#255](https://github.com/gradientspaceai/gtengine-js/issues/255), [#458](https://github.com/gradientspaceai/gtengine-js/issues/458), [#200](https://github.com/gradientspaceai/gtengine-js/issues/200).

### `IntrCircle2Circle2.h`

TI and FI answer different questions: TI is `Length(C0 - C1) <= R0 + R1` (solid
disks) while FI intersects the curves. Nested circles give TI `true` and FI
`false` with `numIntersections = 0`. This is undocumented upstream. It looks
deliberate rather than a coding error, so the port preserves and documents it.
Note that `IntrSphere3Sphere3.h` does **not** share the split; see "Claims
withdrawn or corrected".

Issue [#450](https://github.com/gradientspaceai/gtengine-js/issues/450).

### `IntrConvexMesh3Plane3.h`, `IntrConvexPolygonHyperplane.h`

**Symptom.** For a coplanar face, `GetIntersectionPolygon` returns an arbitrary
permutation of the boundary.

**Cause.** `polygonIndices` is a predecessor map indexed by vertex index, but it
is read back in index order instead of being traversed as a cycle. It is also
sized by edge count while indexed by vertex index, which is an out-of-bounds
write when the face has a non-boundary vertex, and `-1` entries would index out
of bounds.

**Reproduction.** A relabeled cube top face yields side lengths
`1.414, 1, 1.414, 1` with upstream's readback versus the correct `1, 1, 1, 1`.

**Suggested fix.** Walk the cycle. Port: fixed.

Also in that file: `eiVMap[EdgeKey(...)]` default-inserts 0 on a missing key, and
the split-polygon traversal dereferences `find()` without an `end()` check.
Unreachable for well-formed convex input; guarded with asserts in the port.

`IntrConvexPolygonHyperplane.h`'s TIQuery `NEGATIVE_SIDE_VERTEX` branch comment
says "vertex or an edge"; the branch handles only the single-vertex case.

Issues [#301](https://github.com/gradientspaceai/gtengine-js/issues/301), [#250](https://github.com/gradientspaceai/gtengine-js/issues/250).

### `IntrEllipse2Ellipse2.h`, `IntrEllipsoid3Ellipsoid3.h`

**1. `CaseE4NotZero`'s `divisor == 0` branch loses one of two real intersection points.**
The two symmetric points `(y0, +sqrt(lambda/d1))` and `(y0, -sqrt(lambda/d1))`
are both written to `result.points[result.numPoints]`, but `result.numPoints` is
incremented only once, so the second store overwrites the first. The sibling
handlers `CaseE4ZeroE2ZeroE3NotZero` and `CaseE4ZeroE2ZeroE3Zero` increment
correctly for each stored point.

**Reproduction (plain doubles).** `M0 = I` centred at the origin and
`M1 = [[3, 1/2], [1/2, 1/2]]` centred at `(0, 1/2)`: all intermediates are
dyadic, the quartic root `y0 = 1/2` is exact, and `divisor` is exactly zero.
Upstream reports 3 intersection points where the correct answer is 4; the lost
and kept points are `(1/2, +sqrt(3)/2)` and `(1/2, -sqrt(3)/2)`. Port: fixed.

**2. The TIQuery drops the `c_i = 0` terms of `f(s)` (result-corrupting).**
The two critical points at the pole `s = 1/d_j` are lost. The unit circle versus
the ellipse with extents `(1/2, 2)` centred at `(0, 1/2)` is reported
`ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0` although they overlap. This affects any two
axis-aligned ellipses whose centres differ along one axis. Port: fixed.

**3. Floating-point cliffs (preserved).** The FI query returns `numPoints = 0`
and `intersect = false` for two concentric congruent ellipses rotated exactly 90
degrees (for example `x^2/9 + y^2 = 1` versus `x^2 + y^2/9 = 1`, which intersect
in 4 points); finding the standard form with exact matrices recovers all 4, and
perturbing the angle by 1e-6 recovers them too, so the root cause is quartic-root
conditioning rather than a logic error. Downstream, `IntrAreaEllipse2Ellipse2`
then reports `pi*a*b` (the containment case) instead of the true lens area. On
the TI side the root bracket collapses onto its pole for nearly concentric
ellipses (throwing "Unexpected condition."), the quartic loses nearly double
roots, and with the early exit disabled it invents two spurious roots.

**4. `IntrEllipsoid3Ellipsoid3::GetRoots` bracketing asserts are reachable.**
The bracketing uses an ad-hoc `epsilon = 0.001` (upstream's own comment:
`// TODO: What role does epsilon play?`) backed by `LogAssert`s. The query throws
"Unexpected condition." for ordinary input when the centres are close, erratically
in the offset: 0.002 and 0.0001 classify, 0.001 and 1e-6 throw. Preserved, with a
deterministic reproduction pinned.

**5. `IntrEllipsoid3Ellipsoid3`'s `d0 > d1 = d2` branch folds the wrong coefficient (result-corrupting).**
The valid-pair analysis reduces `f(s) = sum_i d_i*c_i/(d_i*s - 1)^2 - 1` by
folding the coefficients that share an eigenvalue, because two equal `d` values
give the single term `d*(c_i + c_j)/(d*s - 1)^2`. The branches for
`d0 = d1 > d2` and `d0 = d1 = d2` do that correctly. The branch for
`d0 > d1 = d2` writes

```cpp
if (param[0].second > (T)0) { valid.push_back(param[0]); }
param[1].second += param[0].second;    // should be param[2].second
if (param[1].second > (T)0) { valid.push_back(param[1]); }
```

folding in the `c` of the *distinct* eigenvalue `d0` (counted twice, since
`param[0]` is pushed as well) and never using `param[2].second`. `GetRoots`
then solves the wrong `f(s)` and the classification is wrong.

**Reproduction (plain doubles, all dyadic).** Ellipsoid0 is the unit ball at the
origin; ellipsoid1 is axis-aligned, centred at `(1/2, 1/4, 0)`, with extents
`(1/4, 1, 1)`. Both frames are the identity and the extents are powers of two,
so `M2 = diag(16, 1, 1)` is exact and its two trailing eigenvalues are exactly
equal. Upstream reports `ELLIPSOID0_CONTAINS_ELLIPSOID1`, although
`(1/2, 5/4, 0)` lies on ellipsoid1 and is `sqrt(29)/4 = 1.3462...` from the
origin; the correct answer is `ELLIPSOIDS_INTERSECTING`. Over 1620
configurations reaching this branch, dense sampling of ellipsoid1's surface
against ellipsoid0's quadratic form disagrees with upstream on 83 and with the
corrected expression on none. The C++ oracle confirms it against the real MSVC
build: `oracle/cpp/cases/v33-intersection.cpp`'s
`IntrEllipsoid3Ellipsoid3.test.equalEigenvaluesDeviation` deviates on 2000 of
2000 records, upstream reporting containment on 1260 and separation on 637 of
them. Port: fixed (`param[1][1] += param[2][1]`). Found by the C++ oracle wave.

**6. Minor.** `CaseE4ZeroE2NotZeroE3Zero` uses `test0 <= test1` where the other
three handlers use `test0 < test1`; `GetRoots` contains dead `fval = F(s)` stores;
`IntrEllipsoid3Ellipsoid3` has a dead `Matrix3x3 D0`.

Issues [#250](https://github.com/gradientspaceai/gtengine-js/issues/250), [#458](https://github.com/gradientspaceai/gtengine-js/issues/458), [#255](https://github.com/gradientspaceai/gtengine-js/issues/255), [#461](https://github.com/gradientspaceai/gtengine-js/issues/461), [#503](https://github.com/gradientspaceai/gtengine-js/issues/503).

### `IntrHalfspace2Polygon2.h`, `IntrHalfspace3Segment3.h`, `IntrOrientedBox2Sector2.h`

**1. `IntrHalfspace2Polygon2`'s FIQuery leaves `Result::polygon` empty in the fully-inside case.**
`negative == 0` returns `intersect = true` with no vertices, contradicting the
field's own doc comment. Likely deliberate copy-avoidance, but it is a trap for
callers.

**2. `IntrOrientedBox2Sector2` discards the clipped polygon when no clipping is needed (result-corrupting).**
After each boundary clip the header does

```cpp
polygon = std::move(hpResult.polygon);
```

and `FIQuery<Halfspace<2,T>, std::vector<Vector2<T>>>` returns `intersect = true`
with an **empty** polygon when the input lies entirely inside the halfspace (item
1). The assignment therefore wipes the working polygon and the query reports no
intersection. This is not a boundary case: the second clip's guard is computed
from the *original* box, so the first clip's output frequently lies wholly inside
the second halfplane. A 300-pair randomized sweep hit it 9 times with the box
demonstrably overlapping the sector. Port: keeps the polygon when the halfspace
query returns an empty result with `intersect = true`.

**3. `IntrOrientedBox2Sector2` is only valid for half-angles <= pi/2.**
The wedge is modelled as the intersection of two halfplanes, which is wrong for a
non-convex sector; `Sector2`'s own default half-angle of pi (a full disk) is
mishandled. Preserved.

**4. `IntrHalfspace3Segment3`'s FIQuery contradicts its own table.**
The doc table says `1 1 0  segment (clipped)` with `numPoints` 0/1/2, but the
code sets `result.numPoints = 1` before the `numPositive` split and returns only
the crossing point; the clipped intersection should be the sub-segment from the
positive endpoint to the crossing point, that is 2 points. The sibling
`IntrHalfspace3Triangle3.h` sets `numPoints` per sub-case and does not have this
defect. Nothing else upstream instantiates this specialization.

Issues [#139](https://github.com/gradientspaceai/gtengine-js/issues/139), [#200](https://github.com/gradientspaceai/gtengine-js/issues/200).

### `IntrIntervals.h`

**Symptom.** The dynamic find-intersection query reports a contact point outside
the second interval.

**Cause.** At about line 485 (interval0 initially left of interval1) the contact
point is computed as `interval0[0] + firstTime*speed0`, the moved LEFT endpoint,
instead of `interval0[1] + firstTime*speed0`, the endpoint that actually makes
contact. The mirrored right-approach branch at about line 503 is correct, so the
query is asymmetric.

**Reproduction.** `[0,1]` moving at speed 2 versus a static `[5,7]` reports
contact at 4 instead of 5.

Port: fixed; restoring the upstream expression fails three tests. A provable dead
store at about line 446 (`intersect = (a0 >= a1)`) is a separate minor item.

Issue [#62](https://github.com/gradientspaceai/gtengine-js/issues/62).

### `IntrLine2Ray2.h`, `IntrLine2Segment2.h`, `IntrLine2SegmentMesh2.h`, `IntrSegment2Segment2.h`, `IntrSegment2AlignedBox2.h`

- A zero-length segment gives a zero-direction supporting line, so
  `IntrLine2Line2` classifies the pair as "same line" and
  `IntrLine2Segment2`/`IntrLine2Ray2` report the collinear case even when the
  point is off the line. The same shape appears in `IntrSegment2Segment2`'s
  `Exact` path, producing a spurious collinear classification.
- `IntrLine2SegmentMesh2` copies `IntrLine2Segment2`'s `+-max()` sentinel into
  `lineParameter` while storing the real endpoint as `point`, so records are
  internally inconsistent. `IntrRay2SegmentMesh2`, which filters on
  `lineParameter >= 0`, then drops the reachable endpoint and keeps the one behind
  the ray.
- `IntrSegment2Segment2`'s FIQuery computes `segment1Parameter = overlap - t`,
  which has the wrong sign for antiparallel collinear segments. Port: fixed.
- `IntrSegment2Segment2`'s `Result` documentation is self-contradictory: the point
  identity and both ascending orderings cannot all hold.
- Minor: `IntrLine2Segment2`'s `line1Parameter[0] >= 0 && line1Parameter[1] <= 1`
  is redundant because both entries hold the same value; `IntrSegment2Segment2`'s
  `Exact` mixes `line0Parameter[0]` and `[1]` (harmless, same value);
  `IntrSegment2AlignedBox2` contains `cdeParameter = cdeParameter;`.

Issues [#197](https://github.com/gradientspaceai/gtengine-js/issues/197), [#203](https://github.com/gradientspaceai/gtengine-js/issues/203), [#458](https://github.com/gradientspaceai/gtengine-js/issues/458), [#461](https://github.com/gradientspaceai/gtengine-js/issues/461).

### `IntrLine3Capsule3.h`

**1. `intersect` is never set when only one hemisphere root is accepted.**
A line tangent to a cap yields `intersect = false` with `numIntersections = 1`,
and the ray and segment clips, which are guarded on `intersect`, drop the
contact. `IntrLine3Cylinder3` sets the flag for one root, so the capsule query is
the outlier. Port: fixed.

**2. The query collapses to a degenerate interval on a cap-junction plane.**
The acceptance regions (wall `|z| <= e`, bottom `z <= -e`, top `z >= e`) overlap
on the junction circles `z = +-e`, so a root there is accepted twice and the
early return at two roots discards the true far endpoint. A segment starting 0.25
inside the capsule is reported as a miss while the TI query says true. Port:
fixed; over 20000 random configurations away from junction circles nothing
changes, and with junction-aimed lines the 654 substantive changes are exactly
the cases an oracle says the old code got wrong.

**3. Minor (preserved).** The FI query degenerates on a zero-length capsule
segment: the centered-form direction is zero, the orthogonal-complement basis
degenerates, and `parameter` is `[-radius, +radius]` regardless of the line. The
distance-based TI query is unaffected.

Issues [#461](https://github.com/gradientspaceai/gtengine-js/issues/461), [#200](https://github.com/gradientspaceai/gtengine-js/issues/200).

### `IntrLine3Cone3.h`

Three result-corrupting defects around lines that pass through the cone vertex.

**1. `CaseC2NotZeroDiscrZero`'s vertex test uses only the U-component.**
The test `t*UdU + UdPmV == 0` is only the axis component of the full condition
`(P - V) + t*U = 0`, so a tangency at a point `X != V` takes the vertex branch
and clamps with height 0 instead of `X`'s height. A truncated cone can therefore
reject a tangent point that lies inside its height range. Port: fixed with the
full componentwise test.

**2. The same branch is gated on an exact floating-point equality that essentially never holds.**
A line through the cone vertex therefore reports a point or an empty set instead
of the intersection ray, because control falls through to the tangency branch.
Reproduction through the segment query: cone apex
`(-2.966273275177028, 0, 0)`, axis `(1,0,0)`, angle 0.15, infinite; segment
`(0,0,0) -> (-0.6862336043960249, 0, 0)`, exactly on the axis and fully inside,
returns type "empty" with `intersect = false`.

**3. `CaseC2NotZeroDiscrPos` Block 3 assumes `c2 > 0`.**
For a line through the cone vertex `Q(t) = c2*(t - tv)^2` has a double root, so
`discr = c1*c1 - c0*c2` is exactly 0 in exact arithmetic and a cancelling
difference with no significant digits in floating point. Rounded above zero with
`c2 < 0` the code calls `SetRayClamp` and reports a segment reaching the cone's
max height where the true set is the vertex alone (recorded case: far endpoint
3.19 from the vertex, `G = -2.94`); rounded below zero it reports nothing for a
line along the axis.

**Suggested fix (applied in the port).** Classify the through-vertex line by the
sign of `c2` rather than by the rounded discriminant. Bit-for-bit comparison over
400k configurations: no changes for plain random lines, and the through-vertex
family went from 29796 containment violations to 0. Lines merely *near* the
vertex, and directions numerically on the cone boundary, remain preserved because
there are no significant digits there.

**4. Minor.** `Result::Convert`'s member templates cannot compile, because
`QFNumber` has no conversion to `Real`; dead code as shipped. The
`isRayNegative` documentation contradicts itself, and `SetRayNegative` carries a
`-1 // +infinity` comment.

Issues [#304](https://github.com/gradientspaceai/gtengine-js/issues/304), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465).

### `IntrLine3Rectangle3.h`, `IntrSegment3Rectangle3.h`

`FIQuery::Result::rectCoord` is `std::array<T,3>` for a 2-coordinate rectangle;
component 2 is always zero (copy-paste from `triangleBary`). Preserved to keep
the result shape. Comment defects: the header says
`X = C + sum_{i=0}^2 s[i]*W[i]` (should be `i = 0..1`) and the zero-`Dot(D,N)`
branches say "Line and triangle are parallel".

Issue [#141](https://github.com/gradientspaceai/gtengine-js/issues/141).

### `IntrOrientedBox2OrientedBox2.h`

TI and FI disagree on edge-touching boxes. TI uses closed-box separation (`>`),
so touching boxes report `intersect = true`; FI's `Outside` requires a strictly
positive vertex distance, so the same pair reports `intersect = false` with an
empty polygon. Comment typos in the same file: "the kast clip vertex", "closed
halfspacedefined by the line".

Issue [#141](https://github.com/gradientspaceai/gtengine-js/issues/141). Port: preserved.

### `IntrPlane3Plane3.h`, `IntrPlane3Cylinder3.h`, `IntrPlane3Triangle3.h`

- `IntrPlane3Plane3`'s parallel test `|Dot(N0,N1)| >= 1` is exact. Only about
  half of normalized doubles satisfy `Dot(N,N) === 1`, so identical planes are
  frequently classified transverse with `invDet` around 1e16. Visible through
  `IntrPlane3Circle3`.
- `IntrPlane3Cylinder3` chooses the CIRCLE versus ELLIPSE label by an exact
  extent comparison; its `Result` docs refer to a nonexistent `type = NONE` and
  to "all zero members" for value-initialized `Line3`/`Ellipse3`; and
  `GetEllipseOfIntersection` discards `FromCoefficients`' status.
- `IntrPlane3Triangle3`'s FIQuery uses a raw uninitialized `Real s[3]` where the
  TIQuery uses an initialized `std::array`. Harmless as used.

Issues [#465](https://github.com/gradientspaceai/gtengine-js/issues/465), [#304](https://github.com/gradientspaceai/gtengine-js/issues/304), [#255](https://github.com/gradientspaceai/gtengine-js/issues/255).

### `IntrRay2Arc2.h`, `IntrSegment2Arc2.h`

**Symptom.** The query reports an arc hit at a point that is not on the arc, and
not even on the circle.

**Cause.** These queries reuse the ray/segment-versus-circle machinery, which
clips against the *solid disk*. With the ray origin (or a segment endpoint)
inside the disk, the clipped "intersection point" is that origin or endpoint
itself. `Arc2::Contains`'s one-argument form assumes its input lies on the
circle, so it happily accepts it.

**Suggested fix (applied in the port).** Intersect with the circular curve
(line-circle roots filtered by `t >= 0` or `|t| <= extent`), the same technique
upstream itself uses in its `SegmentMesh2` queries. `IntrLine2Arc2` is not
affected, because it intersects the curve rather than the disk.

Minor: `IntrSegment2Arc2` has a duplicated `parameter[0] = 0` assignment.

Issue [#304](https://github.com/gradientspaceai/gtengine-js/issues/304).

### `IntrRay2Ray2.h`

TI and FI disagree for opposite collinear rays sharing an origin: TI reports
`numIntersections = 1` via an explicit `t == 0` branch, FI reports 2 with the
degenerate segment `[0,0]`. Both agree on `intersect`.

Issues [#200](https://github.com/gradientspaceai/gtengine-js/issues/200), [#455](https://github.com/gradientspaceai/gtengine-js/issues/455). Port: preserved.

### `IntrSegment2OrientedBox2.h`

**Symptom.** Every reported intersection point is wrong for a rotated box.

**Cause.** The FIQuery computes
`result.point[i] = box.center + (segOrigin + parameter[i] * segDirection)` where
`segOrigin` and `segDirection` are in the box frame (built from
`Dot(diff, box.axis[i])`) while `box.center` is world-space. This is correct only
when the box axes are the identity, which is why the identical line in
`IntrSegment2AlignedBox2.h` is fine.

**Reproduction.** A 45-degree-rotated box produces points like `(+-1, +-1)` where
the correct answer is `(+-sqrt(2), 0)`.

Port: fixed (world-space evaluation). Minor in the same file:
`result.cdeParameter = result.cdeParameter;` is a no-op; a degenerate (point)
segment inside the box reports `numIntersections = 2` while the aligned-box
sibling reports 1 for the same input.

Issue [#255](https://github.com/gradientspaceai/gtengine-js/issues/255).

### `IntrSegment3Ellipsoid3.h`, `IntrSegment3Sphere3.h`, `IntrSphere3Sphere3.h`, `IntrSphere3Triangle3.h`

**1. `IntrSphere3Sphere3`'s FIQuery type-4 branch (internal tangency) reports the antipode.**

```cpp
result.point = sphere1.center + r1 * C1mC0;   // should be sphere1.center - r1 * C1mC0
```

With `sphere0 = ((2,0,0), 1)` and `sphere1 = ((0,0,0), 3)` the query returns
`(-3,0,0)`; the actual tangency point is `(3,0,0)`. The sibling type-6 branch
uses the correct sign. Port: fixed.

**2. `IntrSegment3Sphere3`'s TIQuery misses contained segments.** It only looks
for the segment crossing the sphere surface, so a segment strictly inside the
solid sphere is reported as not intersecting, contradicting the FIQuery in the
same file and `IntrRay3Sphere3.h`'s TIQuery, which handles "origin inside"
explicitly. Port: fixed.

**3. `IntrSegment3Ellipsoid3`'s TIQuery has the same defect** (a segment fully
inside the ellipsoid reports no intersection). Port: fixed the same way
(`qm <= 0 || qp <= 0` tested first). Stale coefficient comments (`a0`, `a1`,
`-a1/2`, `a3*e`) in the same header.

**4. `IntrSphere3Triangle3` returns the sphere *centre*** rather than the surface
contact point for `intersectionType = +1`, inconsistent with the `-1` path. Its
arbitrary-precision path also runs the `closest[j]` loop out of bounds (not
ported; floating-point instantiation only).

Issues [#203](https://github.com/gradientspaceai/gtengine-js/issues/203), [#304](https://github.com/gradientspaceai/gtengine-js/issues/304).

### `IntrSphere3Cone3.h`

The FI query drops the cone vertex from the reported point:
`result.point = t * (cosAngle*D + tmp*B)` omits `V`, although the file's own
reduction parameterizes the ray as `X(t) = V + t*D`. The reported contact point
is wrong whenever the cone vertex is not the origin. Port: fixed (`V + t*D`),
confirmed against an exact centre-to-solid-cone distance for all four cone kinds.
Minor: `DoQueryFiniteCone` calls `GetMaxHeight()` twice redundantly.

Issue [#307](https://github.com/gradientspaceai/gtengine-js/issues/307).

### `IntrTetrahedron3Tetrahedron3.h`

**1. The edge-edge parallelism test compares unnormalised dot products against a cosine cutoff.**
`|Dot(E0, E1)| < cutoff` uses raw edge vectors against a cutoff in `[0, 1]`, so
for edges longer than about one unit the entire edge-edge phase is effectively
skipped. A faithful transcription disagreed with the corrected version on 32 of
4000 random tetrahedron pairs, always as a false "intersecting" report. Port:
normalized comparison.

**2. The edge-edge separation test is not a separating-axis test.**
Upstream tests which side of a plane through one edge endpoint the other
tetrahedron lies on, instead of testing projection-interval disjointness on the
cross-product axis. It misses genuine separations (observed gaps up to 1.12
units, not round-off) and makes the query argument-order dependent: 26 of 2500
random pairs answered differently for `(A, B)` versus `(B, A)`. Port: interval
disjointness plus an explicit zero-cross-product guard; the result agrees with an
independent SAT implementation on all 2500 pairs in both argument orders, and
with a 44-axis brute-force SAT over 166k tie-free integer configurations with
zero disagreements.

**3. Undocumented positive-orientation precondition.** A negatively oriented
tetrahedron gives inward face normals, and a contained tetrahedron is then
reported as separated.

Issues [#307](https://github.com/gradientspaceai/gtengine-js/issues/307), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465).

### `IntrTriangle2Triangle2.h`, `IntrTriangle3CanonicalBox3.h`, `IntrTriangle3Triangle3.h`

**1. Both moving-triangle overloads of `IntrTriangle3Triangle3` use a meaningless parallel test.**
The test is `std::fabs(Dot(N0, N1)) < 1` with `N0` and `N1` the *unnormalized*
edge cross products, so the magnitude scales with the product of the triangle
areas and the comparison against 1 is meaningless. The consequences run both
ways: large non-parallel triangles take the coplanar branch, skipping the `N1`
axis and all nine `E0[i0] x E1[i1]` axes, and can report contact for triangles
that never touch; small parallel triangles take the non-parallel branch, where
every `UnitCross(E0[i0], E1[i1])` is zero. The stationary query in the same
header uses the correct criterion, `|Cross(N0, N1)|^2 > 0`. Port: adopts that
criterion for the moving queries, checked by asserting that the zero-velocity
dynamic query agrees with the stationary query over 1500 random pairs.

**2. The contact-configuration machinery is dead code.** `ContactSide side` and
the contact-time `Configuration tcfg0/tcfg1` are threaded through every
`FindOverlap` call and carefully maintained (including an explicit vertex sort,
about 150 lines) but never read: the contact set is recomputed by the stationary
query at `tFirst`, which makes `FindOverlap` behaviorally identical to
`TestOverlap`.

**3. Unguarded divisions for degenerate triangles.** `GetCoplanarIntersection`,
`IntersectsSegment` and `ContainsPoint` divide by `normal[lookup[2]]` with no
guard; a zero-area triangle0 takes the `numZero == 3` branch and emits NaN or
Infinity coordinates. Upstream's own sign-count table lists `(0,0,3)` as
including the degenerate case but never guards it.

**4. Documentation and convention.** `IntrTriangle3CanonicalBox3.h` (and the
other `IntrTriangle3*Box3.h` files) say the clipped polygon has "at most 7
vertices"; the bound is 9 (a hexagonal plane section plus 3 triangle edges), and
an octagon case exists. The file comment of `IntrTriangle3Triangle3.h` documents
a moving-triangle FI query the header does not contain.
`IntrTriangle2Triangle2::WhichSide` relies on C++ int truthiness via
`--negative`. Finally, the triangle-triangle queries and the tetrahedron
face-normal phase treat measure-zero contact, a shared edge or vertex, as
separation (`WhichSide` classifies `[0, b]`, `b > 0`, as strictly positive); this
is upstream's convention but is not stated.

Issues [#334](https://github.com/gradientspaceai/gtengine-js/issues/334), [#307](https://github.com/gradientspaceai/gtengine-js/issues/307), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465).

### `IsPlanarGraph.h`

**1. Out-of-bounds reads when edges reference invalid vertices.** When input
edges reference out-of-range vertex indices the code sets
`IPG_EDGES_WITH_INVALID_VERTICES` but then proceeds into position lookups using
those indices. Port: returns the accumulated flags early.

**2. `InvalidSegmentIntersection` is asymmetric for degenerate edges.** The
collinear branch parameterizes by the first segment, so a zero-length *first*
segment collapses to "no overlap". With positions `0 = (0,0)`, `1 = (2,0)`,
`2 = (1,0)`: `[[0,1],[2,2]]` gives `DEGENERATE | INVALID_INTERSECTIONS` while
`[[2,2],[0,1]]` gives `DEGENERATE` only.

Issues [#44](https://github.com/gradientspaceai/gtengine-js/issues/44), [#472](https://github.com/gradientspaceai/gtengine-js/issues/472).
### `LieGroupsAlgebras.h`, `RigidBody.h`

**1. `LieSO3::Log`'s angle = pi branch is off by `1/sqrt(2)`.**
Upstream scales the normalized row of `Y + I` by `GTE_C_PI * GTE_C_INV_SQRT_2`.
In this branch `Y + I = (2/pi^2) x x^T`, so each row is parallel to `x` and after
normalization the correct factor is simply `pi`. Equivalently, for
`R = 2 n n^T - I` each row of `R + I` is parallel to `n`. The extra `1/sqrt(2)`
yields a log of length `2.2214...` instead of `pi`, so `Exp(Log(Y)) != Y` for
every rotation by pi: `Log(diag(-1,-1,1))` returns `(0, 0, 2.2214)` instead of
`(0, 0, pi)`. The error propagates into `LieSE3::Log`. Port: fixed.

**2. Knife edge at angle = pi (preserved).** The pi branch requires
`(Trace(Y) - 1)/2 <= -1`; when trace round-off lands just above -1, the generic
branch multiplies round-off-level `Y - Y^T` entries by a factor of about 1e8 and
returns garbage (length about 4e-8 instead of about pi). This is inherent
ill-conditioning of the formulation.

**3. `RigidBodyState::SetQOrientation` builds the rotation matrix from the un-normalized input.**
Upstream executes
`mQOrientation = q; if (normalize) { Normalize(mQOrientation); } mROrientation = Rotation<3, T>(qOrientation);`,
constructing the matrix from the raw parameter rather than the normalized member.
`RigidBody::Update` calls this with `q + h * dq/dt` and `normalize = true` at
every RK4 stage, so the orientation matrix and the world inertia tensors
`R J R^T` are derived from a non-unit quaternion throughout integration, scaling
them by `|q|^2` and `|q|^4`. Port: fixed (matrix from the normalized member); a
regression test fails against the upstream formulation.

**4. Minor.** `SetBodyInertia` has an inconsistent `IsMovable()` guard around
`UpdateWorldInertialQuantities()` compared with `SetQOrientation` and
`SetROrientation` (harmless: immovable bodies carry zero body inertia).
`RigidBody::GetAngularVelocity` has a stray `;;`.

Issue [#313](https://github.com/gradientspaceai/gtengine-js/issues/313).

### `Matrix4x4.h`, `ConvertCoordinates.h`, `Quaternion.h`

**1. `MakePerspectiveProjection`'s `GTE_USE_VEC_MAT` branch reads unassigned entries.**
`M(3,0) = -(M(0,0)*eye[0] + M(0,1)*eye[1] + M(0,2)*eye[2])` reads `M(0,1)` and
`M(0,2)` *before* they are assigned (still 0 from `MakeZero`), so it evaluates to
`-M(0,0)*eye[0]`. By transposition it should read `M(0,0)`, `M(1,0)`, `M(2,0)`.
`M(3,1)` has the same defect (`M(1,2)` unassigned), and `M(3,2)` uses
assigned-but-untransposed indices. The last row is wrong under `VEC_MAT`. Not
reachable in the port, which implements only the `MAT_VEC` branch; that branch
was verified property-wise (projected points land on the plane and on the eye
ray).

**2. `ConvertCoordinates.h`'s affine example.** The four tuples documented as
`B.GetCol(0..3)` are actually the ROWS of `B`, because the example sets
`vectorOnRightV = false` and `UToV(Matrix)` transposes. The port reproduces
`{1,0,0,0}, {0,c,s,0}, {0,-s,c,0}, {2.0,-0.9,-2.6,1}` exactly as `getRow(i)`. The
3D example comments use a 4-tuple `(1,0,0,0)` as a basis vector in a 3D example
and declare `double cs = 0.6, sn = 0.8;` then use `c`/`s`.

**3. `Quaternion.h`'s `Slerp` comment.**
`slerp(t,q0,q1) = [sin(t*theta)*q0 + sin((1-t)*theta)*q1]/sin(theta)` yields `q1`
at `t = 0`; the coefficients are swapped relative to the implementation (and
relative to `Slerp.h`, whose own comment has a different error).

Issue [#160](https://github.com/gradientspaceai/gtengine-js/issues/160).

### `Mesh.h`, `RevolutionMesh.h`, `RectangleMesh.h`, `TubeMesh.h`, `MeshCurvature.h`

**1. `Mesh::ComputeIndices`, SPHERE topology: the south-pole fan uses the wrong row stride.**

```cpp
v0 = (mDescription.numRows - 1) * mDescription.numCols;   // wrong
```

Each row of a SPHERE mesh holds `numCols + 1` vertices (the seam duplicate),
which is why the same function uses `rIncrement = numCols + 1`, and
`RevolutionMesh.h`'s SPHERE initializers confirm the layout. The first vertex of
the last row is therefore `rMax * rIncrement`; the two expressions agree only
when `numRows == 1`. With `numRows = 4, numCols = 5` upstream gets `v0 = 15`
instead of 18, so the fan straddles two rows, reuses an interior seam vertex,
never closes, and leaves the last row's vertices unreferenced by any triangle.
Port: fixed.

**2. The second pole fan is wound opposite to the rest of the mesh (distinct from item 1).**
Both pole fans emit `(v0, pole, v1)`. The first pole sits before row 0, where
that order is consistent with the grid triangles (it plays the role of
`(v1, v3, v2)` with `r = -1`). The second pole sits *after* the last row, where
the consistent order is `(v0, v1, pole)` (the role of `(v0, v1, v2)` with
`r = rMax`, derived by collapsing the virtual row `rMax + 1`). Verified on a
`RevolutionMesh` sphere: exactly `numCols` directed edges, all on the last ring,
are traversed twice in the same direction, the second pole's computed normal
points inward (`dot(n, p)/r = -1`), the last-ring normals degrade to about 0.79
of unit length, and the trailing seam vertex ends up with a zero normal. Port:
fixed.

**3. `TubeMesh::UpdatePositions`'s closed-tube fixup uses the wrong row stride.**
`i1 = col + numCols * (numRows - 1)`, but a CYLINDER row is
`rIncrement = numCols + 1` vertices wide, so the closing copy overwrites an
interior ring and leaves the last ring unclosed (`numRows = 6, numCols = 7`
writes vertices 35..41 instead of 40..47). The loop also stops at
`col < numCols`, leaving the trailing seam vertex stale and cracking the seam.
Port: uses `rIncrement * (numRows - 1)` and copies `col <= numCols`.

**4. `RectangleMesh::InitializeFrame` hardcodes the tangent frame.**
`tangent = (1,0,0)`, `bitangent = (0,1,0)` regardless of the rectangle's axes, so
for any rectangle not in the `z = 0` plane the tangent, bitangent, dpdu and dpdv
channels are not orthogonal to the normal (the normal itself is computed
correctly). Upstream carries its own `// TODO: Are tangent and bitangent
correct?` and a stale `// bitangent = Cross(normal, tangent)` comment the code
does not implement. Port: `tangent = rectangle.axis[0]`,
`bitangent = rectangle.axis[1]`, `normal = UnitCross(tangent, bitangent)`
(identical normal values). `RectanglePatchMesh.h` was checked and does not share
the defect, since it derives the frame from surface derivatives.

**5. `RevolutionMesh`'s DISK texture coordinates never add the declared origin.**
`Vector2<Real> origin{0.5, 0.5}` is used only for the hub; the ring is written as
`radius * (cos, sin)`, covering `[-0.5, 0.5]^2` with the hub at a corner.
`std::min(radius, 0.5)` is vestigial.

**6. `RevolutionMesh`'s one-row SPHERE gets NaN texture coordinates.**
`MeshDescription` clamps SPHERE to `numRows >= 1` but the SPHERE branch divides by
`numRows - 1`.

**7. `MeshCurvature` returns zero-length principal directions at an umbilic vertex.**
With the shape operator `S = k*I` exactly, both candidate vectors
`W0 = {S(0,1), k - S(0,0)}` and `W1 = {k - S(1,1), S(1,0)}` are zero,
`Dot(W0,W0) >= Dot(W1,W1)` is `0 >= 0`, `Normalize` leaves the zero vector, and
both principal directions come back `(0,0,0)`. Reproducible on the unit
icosahedron. The nearby `DWTrnZero` branch returns `U`, `V` for the analogous
planar case. Near-umbilic vertices get `dmin == dmax`, which are then not
orthogonal.

**8. Minor (preserved).** `TubeMesh` rings have only `numCols - 1` distinct
angles, giving one degenerate quad column with a zero normal, and a closed
`TubeMesh` discards the ring at `t = (numRows-1)/numRows`, making the last band
double-width. `MeshCurvature`'s `maxAbs < singularityThreshold` never fires at
the documented threshold of zero. The `Mesh`-family `allowUpdateFrame` re-check
sits inside `if (!mTCoords)`, so a client-supplied tcoord channel skips the
stricter `mNormals` requirement.

Issues [#220](https://github.com/gradientspaceai/gtengine-js/issues/220), [#240](https://github.com/gradientspaceai/gtengine-js/issues/240), [#268](https://github.com/gradientspaceai/gtengine-js/issues/268), [#412](https://github.com/gradientspaceai/gtengine-js/issues/412).

### `MeshStaticManifold2.h`, `MeshStaticManifold3.h`, `StaticVETManifoldMesh2.h`, `StaticVTSManifoldMesh3.h`

**1. `MeshStaticManifold2::GetAdjacentTriangles` returns the adjacency of the wrong edge.**
The 4-tuple stored at `vertex[v0]` for triangle `<v0,v1,v2>` is `{v1, v2, t, a}`,
where `a` is set by `UpdateAdjacencyForEdge` to the triangle adjacent across the
edge *opposite* `v0`, namely `(v1,v2)`, not across the queried edge `(v0,v1)`.
`GetAdjacentTriangles` returns elements `[2]` and `[3]` of a single
`GetDirectedEdge` lookup, so `adj1` is a different edge's neighbour.

Independent evidence: the function's documented case 3, "(invalid, valid)", is
unreachable in the upstream code, because when only `<v1,v0>` exists the code
takes the second branch and returns a valid triangle index as `adj0`.

**Reproduction.** Two triangles `<0,1,2>` and `<0,2,3>` tiling a square: for the
shared diagonal `{0,2}` upstream returns `(1, invalid)`, reporting an interior
edge as a boundary edge.

**2. `MeshStaticManifold3::GetAdjacentTetrahedra` has the same bug plus an index off-by-one.**
The 3D record is a *5*-tuple `{v1, v2, v3, tet, adj}`, but the function reads
`[2]`/`[3]`, the indices the 2D code uses on its 4-tuple, so `adj0` is actually
`v3`, a **vertex** index returned as a tetrahedron index.
`UpdateAdjacencyForFace` in the same header correctly uses `[3]`/`[4]`. Even
shifted to `[3]`/`[4]` it is still wrong for the reason in item 1, and case 3 is
likewise unreachable.

**3. `StaticVETManifoldMesh2::GetAdjacentTriangles` fails to swap on the fallback path.**
The stored triple `<AV, LT, RT>` does describe the queried edge, but after the
`<v0,v1>` lookup fails the code searches `<v1,v0>` and returns its `[1]`/`[2]`
components verbatim, which relative to the queried direction are R then L. The
documented case 3 is again unreachable; it collapses into case 2's signature.

**4. `StaticVTSManifoldMesh3::GetAdjacentTetrahedra` has the same defect one dimension up.**
Records are `<AV0, AV1, LS, RS>`; after `SortFace`, the fallback lookup of the
opposite orientation `<u0,u2,u1>` returns `([2],[3])` unswapped, against the same
documented four-case contract.

**Suggested fix (applied in the port).** Perform both directed-edge or
ordered-face lookups and return the documented contract (`adj0` = L-simplex of
`<v0,v1[,v2]>`, `adj1` = L-simplex of the reversed edge or face). All four
documented cases then become reachable; exhaustive property tests check the
result against a brute-force map.

**5. Documentation.** The comment above `GetComponents` in
`MeshStaticManifold2.h` describes a `range[]` vector that the function does not
have; the signature returns `std::vector<std::vector<size_t>>` directly.

Issues [#66](https://github.com/gradientspaceai/gtengine-js/issues/66), [#361](https://github.com/gradientspaceai/gtengine-js/issues/361), [#472](https://github.com/gradientspaceai/gtengine-js/issues/472).

### `MinimalCycleBasis.h`, `SplitMeshByPlane.h`

**1. Persistent `Vertex::visited` flags flatten nested cycle forests.**
The detach-time `DepthFirstSearch(clone, component)` calls run after
`ExtractConnectedComponents` has reset the visited flags, but the detach-time
searches never reset them among themselves, so a subgraph whose vertices were
already marked in an earlier detach search collapses to just the clone vertex.
Emulating upstream faithfully confirms it: for a doubly nested square
configuration, upstream reports the innermost square as a *sibling* of the square
that contains it, plus two spurious empty-cycle nodes. Affects nesting depth
`>= 2` only. Port: fixed (flags reset per search).

**2-4 (minor, preserved).** `GetClockwiseMost`'s non-convex branch uses
`GetSignDet(...) < 0` where the other five convexity updates in the header and
the published pseudocode use `<= 0` (differs only for exactly parallel
directions). `Extract` returns early when the graph has no edges, so isolated
vertices go unreported. `ExtractConnectedComponents` requires
`adjacents.size() >= 2`, so a lone-edge component is reported as nothing at all,
neither a filament nor isolated vertices.

**5. `SplitMeshByPlane::SplitTriangle*` read the edge map with `std::map::operator[]`,**
which default-constructs index 0 for an absent key and silently corrupts the
output mesh on any internal inconsistency. Provably unreachable for valid input;
the port asserts instead. Cosmetic: a `static_cast<int32_t>(indices.size() / 3)`
narrowing in `ClassifyTriangles`, and the `Vector3` half of the `mEMap` pair value
is never read.

Issue [#310](https://github.com/gradientspaceai/gtengine-js/issues/310).

### `MinimizeN.h`, `RemezAlgorithm.h`

**1. `MinimizeN`'s stray `mDConjIndex = 0` destroys Powell's direction-set update.**
The constructor sets `mDConjIndex(dimensions)`, the slot Powell's method stores
the new conjugate direction in. But each iteration executes `mDConjIndex = 0;`
immediately before the direction-cycling loop:

```cpp
mDConjIndex = 0;
for (int32_t i = 0, ip1 = 1; i < mDimensions; ++i, ++ip1)
{
    mDirections[i] = mDirections[ip1];
}
```

From iteration 2 on, the conjugate direction is written to `mDirections[0]` and
immediately overwritten by the shift, while `mDirections[mDimensions]` is never
updated again: the shift copies the FIRST iteration's conjugate direction into
slot `mDimensions-1` forever, and the last two direction slots become identical.
The direction set degenerates and the method stalls on any non-separable
objective.

**Reproduction.** Rosenbrock from `(-1.2, 1)` on `[-2,2]^2`: upstream freezes at
`f = 3.0548` at `(-0.7436, 0.5650)` after 4 iterations and never improves through
256; without the stray assignment the same run reaches `f = 8.0e-12` at
`(1.0000, 1.0000)` by 16 iterations. Additionally, `mDConjIndex` is never reset
in `GetMinimum`, so a second call on the same object is broken from its first
iteration even before the loop bug bites. Port: the stray assignment is omitted.

**2. `MinimizeN` is basic Powell,** whose direction set can become linearly
dependent and stall in a proper subspace while reporting convergence; confirmed
on a 3-D bowl where every line search is exact yet the search stops at `f = 3.02`.
This is an algorithmic limitation, distinct from item 1.

**3. `RemezAlgorithm::ComputeXExtremes` pins the outer nodes to `xMin`/`xMax`,**
so a function with an inflection point inside the interval converges to a
non-minimax fixed point that is reported as success. `sin(2x + 0.3)` on
`[-1, 1]`, degree 1: reported error 0.33346, true error of the returned
polynomial 0.49382, actual minimax 0.36192. Preserved, since a fix is a different
node-exchange step and there is no in-library consumer.

**4. `RemezAlgorithm::ComputePCoefficients` reads out of range.** It reads
`poly[i]` up to `mPCoefficients.size()-1`, but `Polynomial1`'s operators call
`EliminateLeadingZeros`, so an exactly-zero leading coefficient drops the degree
and the read indexes past `mCoefficient`. The port substitutes 0 for an absent
coefficient, which is mathematically correct and identical to upstream in the
non-degenerate case.

**5. Dead parameter and stale comments.** `maxBracketIterations` is validated and
stored in `mMaxBracketIterations` but never read. The comment above
`ComputeXExtremes` claims a "quadratic-fit line-search (QFLS)", but `GetXExtreme`
is plain bisection bounded by `mMaxBisectionIterations`, and `ComputePartition`'s
comment repeats the claim. It looks as if a bracketing search was replaced by
bisection with the parameter and comments left behind. `Execute` also leaves
`GetXNodes()` one exchange ahead of `GetErrors()` and `GetCoefficients()`.

Issues [#146](https://github.com/gradientspaceai/gtengine-js/issues/146), [#147](https://github.com/gradientspaceai/gtengine-js/issues/147), [#478](https://github.com/gradientspaceai/gtengine-js/issues/478).

### `MinimumAreaBox2.h`, `MinimumWidthPoints2.h`, `RotatingCalipers.h`

**1. `RotatingCalipers::CreatePolygon` discards a genuine corner after a duplicate point.**
Collinearity is tested against the immediately preceding edge, so a duplicate
point's zero-length edge zeroes the `DotPerp` and removes the *next real corner*
as well. Example: the square `(0,0), (0,0), (1,0), (1,1), (1,1), (0,1)` retains
only 2 of its 4 corners and then trips `LogAssert(indices.size() >= 3)`.
`MinimumAreaBox2::RemoveCollinearPoints` carries the same defect, reachable via
the caller-supplied-polygon path, and `MinimumWidthPoints2::ComputeMinWidth`'s
brute-force branch carries it in a currently unreachable form. Port: compares
against the most recent nonzero edge (identical behaviour when no duplicates are
present).

**2. `CreatePolygon` reads `vertices.back()` and `vertices[1]` before any size check;**
the size assert only runs afterwards.

**3. `ComputeAntipodes` emits a duplicated edge and omits another when caliper angles tie (preserved).**
On the unit square the output is `(v2, e3->0), (v2, e0->1), (v1, e2->3), (v1, e3->0)`:
edge `(1, 2)` never appears and `e3->0` appears twice. Reproduced on a second
asymmetric example, the trapezoid `(-1,4), (0,-4), (1,0), (0,8)`. Every emitted
pair is still genuinely antipodal, and min/max-derived quantities (width,
diameter) are unaffected. Flipping the `AngleLessThan` tie from `>` to `>=` does
not repair it, and a real fix would change both the traversal and the output size
contract.

**4. `MinimumAreaBox2`'s `dimension == 1` branch reports wrong extreme indices.**
`imin`/`imax` are left at 0 although the degenerate line's origin is
`points[hull[0]]`, so `GetHull()` reports an interior point as an extreme unless
`hull[0] == 0`. Reproduction: t-values `[4, 0, 7]` make upstream report
`imin = 0` for a point with `t = 4`.

**5. The degenerate branches return stale `mArea`/`mSupportIndices` from a previous query.**

**6. The convex-polygon overload never assigns `mNumPoints`/`mPoints`.**
Overloads 1 and 2 begin with `mNumPoints = numPoints; mPoints = points;`;
overload 3 only does `mHull.clear()`. A caller supplying its own convex polygon
therefore gets `GetPoints() == nullptr` on a fresh functor, or the previous data
set's pointer on a reused one, while `GetHull()` and `GetSupportIndices()` refer
to the polygon just processed. The documented `mPoints[hull[...]]` lookup then
dereferences null or indexes the wrong array. Port: fixed.

**7. Minor.** Overload 4 contradicts overload 3 on empty polygon indices;
`MinimumWidthPoints2` uses inconsistent frame handedness between its degenerate
and general branches (harmless, zero extent); stale `GetSupportIndices` and
`ComputeAngles` comments.

Issues [#286](https://github.com/gradientspaceai/gtengine-js/issues/286), [#328](https://github.com/gradientspaceai/gtengine-js/issues/328), [#402](https://github.com/gradientspaceai/gtengine-js/issues/402).

### `MinimumAreaCircle2.h`, `MinimumVolumeSphere3.h`

**1. The trapped-failure fallback bounds only a prefix of the input.**
In the exception-trap branch, `GetContainer(numPoints, points, minimal)` is
called with `numPoints` *after* it has been overwritten by the unique-point
count, while `points` still refers to the full input array. The fallback
"bounding" circle or sphere is then computed over only a prefix and can fail to
bound the point set; with a duplicated point it can return the radius-0 circle
about that point. Port: the complete array is passed.

**2. An untrapped floating-point failure returns success with a non-containing circle or sphere.**
`UpdateSupport{2,3,4}` rewrite `mSupport`/`mNumSupport` unconditionally, while
`operator()` keeps the returned circle only when it is larger. A rejected update
therefore leaves a point hidden inside the support set (or, in 3D, simply never
rechecked). Two triggers:

1. a cocircular quadruple whose squared radius differs by one ulp between two
   triples;
2. a nearly collinear triple whose 2x2 determinant is a cancellation residue
   instead of exact zero, so `LinearSystem::Solve`'s `det != 0` test returns
   garbage (a radius-0 circle) instead of the `radius = max()` sentinel.

**Reproductions.** A 2D nine-point set misses `(2,-6)` by 0.23; a 2D collinear
four-point set; a 3D seven-point set misses a support point by 0.084; a 3D
eleven-point set misses a non-support point by 0.40. Rates: about 1e-4 on 2D
lattices, about 0.05% on 3D `[-2,2]^3` lattices, none observed for random real
coordinates.

Preserved in the port: a final containment check would false-trap about 45% of
correct results (the returned radius is a square root, so re-squaring loses the
tie), and restoring the support set on rejection does not repair either
reproduction. Instantiating upstream with `ComputeType = BSRational` would
eliminate the defect.

Issues [#286](https://github.com/gradientspaceai/gtengine-js/issues/286), [#399](https://github.com/gradientspaceai/gtengine-js/issues/399).

### `MinimumSpanningTree.h`

- An empty `edges[]` yields `numVertices == 0`, and
  `ExtractMinimumSpanningTree` executes `records[0] = heap.Insert(...)` on a
  zero-length vector: an out-of-bounds write. `Execute` has no nonempty
  precondition even with `validateInputs` enabled. Port: guards and returns an
  empty tree.
- Documentation: the header describes the output as a spanning *tree* with "the
  first element ... a sentinel", but a disconnected graph yields a spanning
  *forest* with one sentinel per component.
- Documentation: back edges are ordered by remapped indices, not labels, because
  the `v0 < v1` filter runs before `ConvertToOriginalIndices`, so reported back
  edges are neither sorted nor smaller-first in the caller's labels.

Issues [#74](https://github.com/gradientspaceai/gtengine-js/issues/74), [#472](https://github.com/gradientspaceai/gtengine-js/issues/472).

### `MinimumVolumeBox3FloatingPoint.h`, `MinimumVolumeBox3Rational.h`

**1. The `dimension == 2` Newell-normal loop is not closed (both files).**
`for (i0 = numHull - 1, i1 = 1; i1 < numHull; i0 = i1++)` drops the wrap-around
term involving `hull[0]` from the cross-product sum. For a coplanar point set
whose hull is a triangle the sum degenerates to `Cross(P2, P1) + Cross(P1, P2)`,
which is exactly zero, so the plane normal is exactly zero and the returned box
is garbage for *every* such input. Port: fixed by starting at `i1 = 0`; the
regression test fails with the upstream loop.

**2. `MinimizerVariableT` declares floating-point parameters inside the exact pipeline (Rational file).**
Lines about 1177-1179 declare `T const& tminNumer, T const& tmaxNumer, T const& tDenom`
where the sibling `MinimizerVariableS` (line about 1138) was correctly changed to
`Number const&`. Every caller passes `Number` (`BSNumber`) expressions, and
`BSNumber::operator double()` is implicit, so the code compiles and silently
rounds every exact numerator and denominator to `double` before converting back
into `std::vector<Number>`. Every `t`-variable level curve in the rational
pipeline is therefore sampled at rounded parameters, defeating the purpose of the
exact instantiation. Port: fixed.

**3. `ComputeVolume` assumes a hull-edge vertex realizes the `axis[0]`/`axis[1]` minima (floating-point file).**

```cpp
candidate.minSupportIndex[0] = mEdges[candidate.edgeIndex[0]].v[0];
pmin[0] = Dot(candidate.axis[0], mTVertices[candidate.minSupportIndex[0]]);
```

That is true in exact arithmetic but false in floating point. When
`d = f00*f11 - f10*f01` is at the rounding level, the endpoint sample of
`MinimizerVariableS/T` computes `q0` and `q1` as cancellation noise (both exactly
zero in exact arithmetic; the residue is `fl(fl(f00-f10) - f00) != -f10`), so
`q1*M0 - q0*M1` is a random direction, the volume is underestimated, and that
candidate wins. In 5 of 300 random 9-point integer clouds the box missed a point
by up to 3.26 on a coordinate range of 10. Port: the two minima use the existing
`GetExtreme` hill climb, identical when the assumption holds and correct when it
does not (about 19% slower). The Rational pipeline is unaffected.

**4. `GetExtreme`'s hill climb stalls on a floating-point plateau (floating-point file).**
It uses a strict comparison and stalls at a hull vertex lying in the relative
interior of a face, because `RemoveCoplanarTriangleAdjacencies` removes such a
vertex only under the exact `BSNumber` test, which rotated coordinates defeat.
Reproduction: the 8-point cloud
`(3,1,1) (0,0,0) (0,-1,-3) (2,-1,3) (-2,0,1) (0,0,0) (0,-4,0) (4,3,-1)`
gives a correct box (volume 116.667). Rotating it by about 0.1 degree about y
with the frame
`(0.9999983994110458, 0, -0.0017891828711411968), (0, 1, 0), (0.0017891828711411968, 0, 0.9999983994110458)`
gives, at `lgMaxSample` 2, 3 and 4, a box with extents
`(2.4807, 8.5e-17, 3.3741)`, volume 5.7e-15, and a maximum containment violation
of 3.90 units. Simulating upstream's `ComputeVolume` verbatim reproduces it
(volume 7.13e-16, violation 3.8966). `MinimumVolumeBox3Rational` is unaffected
(exact dots). Port: after the upstream climb, flood through vertices within
`8*eps*max|d|*L1` of the best value, keeping strict improvements; 4500 queries
over 1500 rotated clouds then show zero containment failures, against 4 without
the fix.

**5. Minor and documentation.** `MinimizerConstantS` uses
`T const half = static_cast<Number>(0.5)` where the siblings use
`static_cast<T>(0.5)`; `ComputeVolume` computes
`axis[2] = Cross(axis[0], axis[1])` twice with no intervening change (both
files); `RemoveCoplanarTriangleAdjacencies` has a redundant `numAdjacent = 0`;
stale comments in five level-curve processors (0x09's `smax = 0` should read 1;
0x69 has the wrong denominator in the `1 - ta` derivation; typos in 0x50, 0x40,
0x96); `&&indices != nullptr` is missing a space at about line 168 of the
rational file. All 79 level-curve processors and both 81-entry dispatch tables
were compared mechanically between the two files with zero semantic differences.

Issues [#352](https://github.com/gradientspaceai/gtengine-js/issues/352), [#355](https://github.com/gradientspaceai/gtengine-js/issues/355), [#405](https://github.com/gradientspaceai/gtengine-js/issues/405), [#426](https://github.com/gradientspaceai/gtengine-js/issues/426).

### `NaturalSplineCurve.h`, `ParametricCurve.h`, `PolynomialCurve.h`, `TCBSplineCurve.h`, `SampleCircularArc.h`, `ReparameterizeByArclength.h`, `NURBSCircle.h`

**1. `NaturalSplineCurve::CreateClosed`'s wrap row overwrite breaks C1 for 3-point closed splines.**
The wrap-around matrix row is written with three plain assignments to columns
`numSm1`, `0`, `1`. Those columns are distinct only for `numSegments >= 3`; for a
3-point closed spline (`numSm1 == 1`) the third assignment overwrites the first,
silently dropping the `dt[numSm1] * C[numSm1]` term, so the resulting "closed"
spline is not C1 at the wrap. Port: fixed by zero-filling the row and
accumulating with `+=`, which is bit-identical for `numSegments >= 3`.

**2. `NaturalSplineCurve`'s constructors validate after damage.** `LogAssert`
runs after the `ParametricCurve(numPoints - 1, times)` base constructor has
already resized on a negative count and dereferenced a possibly-null `times`.
Port: validates before delegating. `CreateFree`'s `storageSize` accounting does
not match the subsequent pointer arithmetic; it over-allocates by `N-1` reals by
accident.

**3. `TCBSplineCurve` and `PolynomialCurve` never set `mConstructed = true`,**
so `ParametricCurve::operator bool()` reports failure for every valid curve.
`Evaluate` is unaffected. Port: fixed in both. `PolynomialCurve::Evaluate` also
uses `if (order == 3)` for the third derivative while the lower orders use `>=`,
so `order > 3` leaves `jet[3]` untouched (the documented precondition is
`order <= 3`).

**4. `TCBSplineCurve::ComputeInteriorTangents` divides by zero in the lambda pass.**
When both Kochanek-Bartels tangents at a key frame vanish (for example
`P2 == P0`, uniform times, zero TCB), `0 * infinity` makes both tangents NaN and
poisons two segments' coefficients. Without lambda the same input gives zero
tangents.

**5. `ParametricCurve::GetTime` ignores the segment decomposition `GetLength` honors.**
`GetLength` splits at knots and uses cached segment lengths, but `GetTime`
(around line 213) builds `F(t) = Romberg(mTime.front(), t, speed) - length` as a
single integration across the whole interval. On any multi-segment curve with a
speed or derivative discontinuity at a knot, exactly the case the multi-segment
constructor exists for, Richardson extrapolation assumes smoothness and converges
poorly, making `GetTime` inconsistent with `GetLength`: measured arc length
12.0048 versus the exact 12 on a piecewise-linear curve. It also re-integrates
from `tmin` on every bisection step. This affects arc-length reparameterization
of Bezier, B-spline and NURBS curves at knot boundaries.

**6. `ParametricCurve`'s lazy-init sentinels use length zero** (around lines 152
and 199). `GetLength` re-runs the full segment-length initialisation whenever
`mSegmentLength[0] == 0`, and `GetTotalLength` whenever
`mAccumulatedLength.back() == 0`. A curve with a legal zero-length first segment
(a repeated B-spline knot) therefore re-integrates everything on every call.
Performance only; values are identical.

**7. `SampleCircularArc::SampleArc3` / `SampleArc4` split directions are wrong for arcs larger than the midpoint case.**
The "trisector" and "quadsector" directions are formed as `-(2*P0 + P2)`,
`-(3*P0 + P2)` and so on, which are exact only for the 1:1 midpoint
construction. For a 270-degree arc the split directions land at 153.43 and 116.57
degrees instead of 90 and 180 degrees, so the sub-arcs run backward and exceed
the pi/2 precondition of the base sampler; the resulting samples lie exactly on
the circle but are not angularly ordered. A full-circle arc additionally yields
NaN. Arcs of at most pi radians are correct and monotone. Preserved, since a fix
would mean inventing new math.

**8. `ReparameterizeByArclength::DoNewtonsMethod`** falls through to
`tMid - fMid/dfdt` after the `dfdt == 0` bisection step. `fMid` is nonzero there,
so the quotient is `+-Infinity` and the out-of-interval test recomputes the same
midpoint. The final result is correct, by a confusing route.

**9. `NURBSCircle::NURBSHalfCircleDegree3`'s comment names the wrong half.**
It says `x >= 0`; the control points trace `y >= 0`.

Issues [#113](https://github.com/gradientspaceai/gtengine-js/issues/113), [#182](https://github.com/gradientspaceai/gtengine-js/issues/182), [#183](https://github.com/gradientspaceai/gtengine-js/issues/183), [#295](https://github.com/gradientspaceai/gtengine-js/issues/295), [#415](https://github.com/gradientspaceai/gtengine-js/issues/415).

### `NearestNeighborQuery.h`

`NearestNeighborQuery` allows `maxLevel == 32` but uses a fixed 32-entry
traversal stack, permitting depth 33; the suppressed MSVC C28020 warning in the
header points at exactly this. `GetRectangle`-style accessors also do not reject
negative indices (an out-of-bounds read in C++; the port returns null). The
header's own TODO acknowledges that the query is approximate.

Port status for the stack bound: preserved. The port keeps upstream's structure
(a 32-entry array allocated at construction, with `maxLevel <= 32` asserted), but
a JavaScript array grows on write, so an overflow cannot corrupt memory. The
condition was judged unreachable in practice in any case: it needs on the order
of 2^32 sites.

Issue [#48](https://github.com/gradientspaceai/gtengine-js/issues/48).

### `OrientedBoxTreeOfTriangles.h`

**Symptom.** The tree culls leaves it should visit, so `Execute` silently drops
intersections.

**Cause.** In `ComputeLeafBoundingVolume`, the "find the smallest extent" scan's
last comparison is `>` instead of `<`:

```cpp
absExtent = std::fabs(box.extent[2]);
if (absExtent > minAbsExtent) { minAbsExtent = absExtent; minIndex = 2; }
box.extent[minIndex] = static_cast<T>(0);
```

This is a copy-paste slip from the correct "find the maximum" scan in
`OrientedBoxTreeOfSegments.h`. Whenever `extent[2]` is not the smallest, the
*largest* extent is zeroed, the leaf box collapses along its longest axis and no
longer contains its own triangle.

**How often.** `GetContainer` returns axes ordered by increasing variance, so
`extent[2]` is the largest for essentially every triangle: upstream zeroes the
largest extent every time, not occasionally. Restoring the upstream comparison
fails 7 of the port's tests.

Issue [#343](https://github.com/gradientspaceai/gtengine-js/issues/343). Port: fixed.

### `PdeFilter.h`, `PdeFilter1.h`, `PdeFilter2.h`, `PdeFilter3.h`

**1. `ScaleType::NONE` still subtracts the data minimum,** so "as is" actually
shifts the image, and a constant image stores as all zeros because the
`min == max` branch leaves `mMin` at the data value. Reproduces in
`PdeFilter1`, `PdeFilter2` (see `PdeFilter2.h:59`) and `PdeFilter3`.

**2. The Neumann image border is stale after the first buffer swap.**
`AssignNeumannImageBorder` runs only in the constructor and `OnPreUpdate`
refreshes only the *mask* border, so after the first swap the boundary derivative
estimates use construction-time duplicated values instead of tracking the
evolving interior.

**3. Minor.** The Dirichlet and Neumann border assignments write `mBorderValue`
verbatim into buffers that hold `mOffset`/`mScale`-transformed data, so with any
non-identity `ScaleType` the border value is in different units than the interior
image.

Issue [#60](https://github.com/gradientspaceai/gtengine-js/issues/60). Port: preserved, pinned by tests.

### `PlanarMesh.h`

- The first constructor half-constructs on duplicate triangles: it silently
  leaves `mNumVertices = 0`, a null vertex pointer and an unset `mQuery`. The
  adjacent comment, plus its own `// TODO: Fix this comment`, is stale.
- `tmap.find(key)` and `mTriIndexMap.find(key)->second` are dereferenced without
  an `end()` check.
- `Contains` has no range check on `triangle`, although every sibling accessor
  does, giving an out-of-bounds read in the documented exhaustive-search usage.

Issue [#256](https://github.com/gradientspaceai/gtengine-js/issues/256). Port: asserts; no behaviour change for valid input.

### `Polynomial1.h`

**1. `SquareFreeFactorization` never terminates for floating-point `Real`.**
The trailing `do { ... } while (b.GetDegree() > 0);` is unbounded, and
`GreatestCommonDivisor` tests remainders against exactly zero, so a few-ULP
rounding error makes it report the constant 1 for a pair with a genuine common
factor; then `b = b / 1 = b` never loses degree.

**Reproduction (exact integer coefficients).**
`f = (t-2)^2 (t+2)(t+3)(t-4) = [-96, 40, 36, -14, -3, 1]`. `gcd(f, f')` comes out
as `t - 1.9999999999999998`; the next remainder is
`[4.26e-14, 3.55e-15, -2.44e-15]` instead of 0; `b` is then pinned at
`[48, 4, -16, -1, 1]` while `d` grows linearly forever. Confirmed by tracing 40
iterations and by an unbounded run that had to be killed after 2 minutes. Port:
caps iterations at `degree(f) + 1`, the exact-arithmetic bound of one iteration
per multiplicity, and throws; the algorithm is otherwise untouched.

**2. `GreatestCommonDivisor`'s exact-zero remainder test is unusable at double precision.**
Same root cause, non-fatal: for `f = (t-1)(t-2)^2 = [-4, 8, -5, 1]`, `gcd(f, f')`
should be `t - 2` but computes as 1, and square-free factorization returns `f` as
a single factor. A fix would need a tolerance parameter upstream does not have,
so the port does not invent one.

**3. `GreatestCommonDivisor` is asymmetric for degree-0 inputs.**
`gcd(c, 0) = c` but `gcd(0, c) = 0`, because the tie in
`p0.GetDegree() >= p1.GetDegree()` goes to the first argument. Only the
degree-0/degree-0 case is affected.

Issues [#83](https://github.com/gradientspaceai/gtengine-js/issues/83), [#488](https://github.com/gradientspaceai/gtengine-js/issues/488).

### `PolylineOffset.h`

The member initializer list sizes the arrays before `LogAssert` can fire
(around lines 61-73):

```cpp
mDirections(isOpen ? vertices.size() - 1 : vertices.size()),
```

`vertices.size()` is `size_t`, so an empty input with `isOpen == true` wraps to
`SIZE_MAX` and the vector constructions throw `length_error`/`bad_alloc` before
the intended diagnostic. The `resize` calls in the body are then dead; they
always resize to the size the members already have.

Issue [#112](https://github.com/gradientspaceai/gtengine-js/issues/112). Port: fixed (validates first, so an empty input throws the documented "Invalid number of polyline vertices.").

### `PrimalQuery2.h`, `PrimalQuery3.h`

**1. The four-argument `ToLine(test, v0, v1, order)` swaps collinear orders 0 and +2.**
The collinear branch computes:

```cpp
Real x0x0 = x0 * x0;          // x0,y0 = P - V0
Real y0y0 = y0 * y0;
Real sqrLength = x0x0 + y0y0; // |P-V0|^2, but |V1-V0|^2 is intended
```

It squares `P - V0` instead of `V1 - V0`. With `P - V0 = c*(V1 - V0)`, the test
`dot > sqrLength` reduces to `c > c^2`, which is true for `0 < c < 1`, so `order`
becomes +2 ("beyond the segment") for points *interior* to `[V0, V1]`, and 0
("interior") for points beyond `V1`. `P = V1` (+1) and the `dot <= 0` cases stay
correct by coincidence. The sibling `ToLineExtended` performs the same test
correctly against `|Q1-Q0|^2`, which is the evidence this is a slip. Nothing in
`GTE/Mathematics` calls this overload.

**2. The BSRational precision tables in the comments are stale.**
The N-value tables in the `PrimalQuery2.h` and `PrimalQuery3.h` comments were
computed with an older `BSPrecision` API (they cite a since-removed four-argument
bool constructor) and do not match the current 8.0 formulas; for example `ToLine`
float is documented as 35 where the current formulas give 70. Note also the
`BSPrecision::operator+` finding above: the published N values are computed
through chains that mix operand sets, so they may be under-estimates.

**3. Comment defects.** The `order` table mislabels +2 as relative to `V0` and
omits +3; `ToCircumcircle`'s "involves three calls of ToLine" is false, since it
is an inline 3x3 determinant, and its own N-values contradict the claim.

Issues [#100](https://github.com/gradientspaceai/gtengine-js/issues/100), [#43](https://github.com/gradientspaceai/gtengine-js/issues/43), [#101](https://github.com/gradientspaceai/gtengine-js/issues/101), [#366](https://github.com/gradientspaceai/gtengine-js/issues/366).

### `QuadricSurface.h`

`GetClassification` (around line 200) symmetrizes `A` from the upper triangle
while `F`, `FX`, `FY` and `FZ` use the stored matrix verbatim, so an asymmetric
input classifies one quadric and evaluates another. Minor in the same file:
`rS12` (around line 210) is computed and never used, and the comment around line
195 claims Sturm sequences where the code implements Descartes' rule of signs
(still exact for a real symmetric matrix).

Issue [#318](https://github.com/gradientspaceai/gtengine-js/issues/318).

### `RiemannianGeodesic.h`

- `ComputeMetricDerivative` computes
  `mChristoffel1[d](i0,i1) + mChristoffel1[d](i1,i0)`, which equals
  `2*Gamma_{d,i0i1}`, not `dg_{i0i1}/dx_d`, under the
  `Gamma_{k,ij} = Dot(P_ij, P_k)` convention that `EllipsoidGeodesic.h` uses and
  `ComputeChristoffel2` requires. `mMetricDerivative` is never read anywhere in
  GTE, so the quirk is preserved rather than reinterpreted.
- `mIntegralStep`, `mSearchStep` and `mDerivativeFactor` are computed only in the
  constructor although their inputs are public and documented as tweakable;
  changing `integralSamples` afterwards makes `ComputeSegmentLength` return
  roughly double the true length. The port adds an
  `updateDerivedParameters()` escape hatch; defaults are unchanged.
- `ComputeGeodesic`'s `LogAssert(subdivisions < 32)` admits 31, for which
  `1 << subdivisions` is signed overflow; negative values are not guarded.

Issue [#295](https://github.com/gradientspaceai/gtengine-js/issues/295).

### `RootsBisection1.h`, `RootsBisection2.h`

**1. `RootsBisection1` returns uninitialised outputs when `maxIterations == 1`.**
With `mMaxIterations == 1` the loop
`for (iteration = 2; iteration <= mMaxIterations; ++iteration)` never runs and
both output references `tRoot` and `fAtTRoot` are returned unassigned, so the
caller reads back whatever it passed in. Port: returns explicit zeros, matching
the zeros upstream returns on its `sign0 == sign1` path.

**2. `RootsBisection2` has the same analogue with a stale-value twist.**
With `xMaxIterations == 1` (or `yMaxIterations == 1`) the inner bisector returns
2 without ever writing its root and value outputs, so upstream's members keep a
stale value from a *previous* call while the return value claims a successful
bisection.

**3. `mNoGuaranteeForRootBound = (numYIterations == 0);` inside `XFunction` assigns instead of accumulating.**
The x-bisector calls `XFunction` repeatedly, so an early failed y-bisection is
erased by any later successful one. Only the x-status is OR-ed in (`|=`), which
strongly suggests accumulation was intended for y too, plus a `false` reset at
the top of `operator()`. It does not affect the computed roots.

**4. The y-outputs need not correspond to the returned x-root.**
`mYRoot` and `mGAtRoot` hold whatever the LAST `XFunction` call produced. Since
`RootsBisection1` evaluates `F(tMin)` then `F(tMax)` before testing signs, an
exact root at `xMin` returns `xRoot = xMin` while `yRoot` came from `x = xMax`,
so `gAtRoot != G(xRoot, yRoot)`. The same mismatch occurs when the x-precondition
fails.

Issues [#84](https://github.com/gradientspaceai/gtengine-js/issues/84), [#152](https://github.com/gradientspaceai/gtengine-js/issues/152).

### `RootsCubic.h`, `RootsQuartic.h`, `RootsGeneralPolynomial.h`, `PolynomialRoot.h`

**1. `RootsCubic` uses an invalid root bound in the bisection branch.**
The `signDelta < 0` branch bisects over `[-b, b]` with
`b = max(1, |d0|, |d1|)`, which is not a valid root bound for monic
`x^3 + d1*x + d0`. Example: `x^3 - 0.9x - 0.9` has its only real root at
`1.26863...`, but `b = 1`; `PolynomialRootBisect` sees the wrong sign at `xMax`
and returns the endpoint, so upstream reports `1.0`, an absolute error of 0.27.
`RootsQuartic` inherits the defect through the resolvent cubic. Port: Cauchy's
bound `1 + max(|d0|, |d1|)`.

**2. `RootsQuartic`'s complex-pair classification is incomplete.**
`if (signDelta > 0 && rD2.GetSign() > 0) { return 0; }` is only half the
criterion: a positive discriminant means four real roots iff `P = 8*d2 < 0` AND
`D = 64*d0 - 16*d2^2 < 0`; otherwise there are two complex pairs. Found by
fuzzing against `RootsGeneralPolynomial`/`RootsPolynomial`: the quartic with
coefficients
`(2.3588208672590554, 1.996267369017005, -0.8864723690785468, -2.481981527991593, 1.4308462475892156)`
has no real roots, but upstream reports four, with residuals 1.6 and 2.2 at the
"roots". Port: `signDelta > 0 && (d2 >= 0 || 4*d0 - d2^2 >= 0)` classifies the
complex-pair case.

**3. `RootsQuartic` reads stale values from a reused root array.**
A single `std::array<PolynomialRoot<Rational>, 2>` is reused for every
square-root extraction and index 1 is read unconditionally, but
`ComputeDepressedRoots` writes index 1 only for a strictly negative argument, so
a zero argument leaves index 1 holding the previous extraction's value
(`sqrt(0)` should be 0). Latent alone; combined with item 2 it produces the
duplicated bogus roots above.

**4. `RootsGeneralPolynomial` runs on a zero-padded rational polynomial.**
After trimming high-order zero coefficients,
`std::vector<Rational> rP(p.size())` (and `rPMonic(rP.size())`, plus the
`rPMonic = rP` copy) allocate to the *untrimmed* size while only entries
`0..degree` are assigned, so the solver runs on a zero-padded polynomial whose
leading coefficient is zero, violating the assumption behind the Cauchy bound and
the degree recursion. Values survived in every configuration exercised, since
padding does not change `p(x)`, but the recursion burns extra
arbitrary-precision levels on phantom degrees.

**5. `PolynomialRootBisect` collapses the interval for endpoint roots.**
The sign-mismatch branches set `xMax = xMin` or `xMin = xMax`, so an exact
endpoint root (true sign 0 at an endpoint) is never reported as a bracketed root,
only as the collapsed point. Documented fallback behaviour.

**6. Minor and documentation.** The cubic closed form's `sin(theta/3) <= 0`
reorder is unreachable, since `atan2(positive, x)` lies in `(0, pi)`; the cubic
bisection labels its second Samuelson subinterval `[-2s, s]` instead of
`[-s, s]`; `RootsQuadratic<T>` versus `RootsQuadratic<Rational>` qualification of
the same static is inconsistent; `RootsGeneralPolynomial::Bisect` has dead stores
`tPMax/tPMin = tPAtRoot`; `PolynomialRoot.h` has the typo comment
`// x is the root estimate)and m ...` and an `operator==` that ignores the
multiplicity `m` (consistent with `operator<` for `std::set` usage).

Also characterised, not filed separately: `RootsPolynomial` with
`Rational = double`, upstream's own template parameter, misclassifies near-double
roots both ways (`-x(x+1)^2` loses its double root; `(1, -2/3, 1/9)` gains a
spurious double root at 3). This is the conditioning `IntrEllipse2Ellipse2`
inherits.

Issues [#340](https://github.com/gradientspaceai/gtengine-js/issues/340), [#319](https://github.com/gradientspaceai/gtengine-js/issues/319), [#488](https://github.com/gradientspaceai/gtengine-js/issues/488).

### `Rotation.h`, `RotationEstimate.h`, `Projection.h`

**1. `operator()(i0, i1, i2)` returns relabelled angles for a re-factorization.**
The operator writes the requested axis indices into `mEulerAngles.axis` *before*
the `switch`, and the `IS_EULER_ANGLES` arm is a bare `break;`. A `Rotation`
constructed from Euler angles that is asked for a *different* axis order
therefore returns the original angles relabelled with the new axes, which is a
different rotation. Port: returns the cache unchanged when the requested axes
match the stored ones (bit-identical to upstream) and otherwise recomputes
through the matrix, as every other source type does.

**2. `Convert(Matrix<3,3,T>, AxisAngle<3,T>&)` underflows near angle = pi.**
For angles near pi the axis is taken from the antisymmetric part
`R - Transpose(R)`, whose entries are `2 sin(theta) a_i`. Near pi these are tiny
and, once the products fall into the denormal range, the extracted axis is
non-unit. For the rotation matrix of quaternion `(0, 0, -1, 0)` the antisymmetric
part is exactly zero and the axis returned is the zero vector; for nearby inputs
the "unit" axis has length about 0.707. Converting the result back to a matrix
yields `-I` instead of `R`. Port: validates the normalized axis and falls back to
the symmetric formula (the largest diagonal entry of `(R + Transpose(R))/2`)
when the antisymmetric extraction is degenerate; the check is on the result, so
no input upstream handles correctly changes.

**3. Residual (preserved).** Slightly further from pi the antisymmetric entries
are small but not denormal; round-off absorbed into them gives a unit axis that
is wrong by up to about 5e-3. Fixing this needs an input threshold, a numerical
policy change.

**4. The Euler factorization gimbal-lock test uses exact `+-1`.**
`Convert(Matrix, EulerAngles)` tests `r(i,j) == +-1` to detect the NOT_UNIQUE
case. One ulp inside the boundary the UNIQUE branch runs on pure round-off and
returns angles whose recomposition differs from `R` by up to 2.0 in Frobenius
norm.

**5. `RotationEstimate`'s `C_ROTC*_EST_MAX_ERROR` tables understate the real error.**
Recomputed `max |poly(t) - rotc_k(t)|` on `[0, pi]` in 60-digit decimal
arithmetic (Maclaurin-series references, no cancellation) using the header's own
coefficients. Degrees 4 to 12 of `rotc0`/`rotc1`/`rotc2` and 4 to 10 of `rotc3`
match the table to every printed digit. `rotc0`/`rotc1`/`rotc2` at degrees 14 and
16 and `rotc3` at degrees 12, 14 and 16 are off by factors of 3 to 1004 (for
example `rotc0` degree 16: tabulated `6.80012e-16`, true `6.82690e-13`), and
**every** `rotc4` row is wrong, by factors of 607 to 1.5e12. Several `rotc4`
entries (`1.16e-16`, `1.50e-19`, `1.59e-22`) are below what a `double` can
represent as an absolute error of a quantity of magnitude 1/6. The constants are
returned verbatim by the `GetRotC*EstimateMaxError` functions and used in no
computation. All 245 `RotationEstimate` constants were independently verified to
match the header byte for byte, so this is a table-generation problem, not a
transcription one.

**6. `C_ROTC4_EST_COEFF` degree 14, last coefficient `-2.753603733509153125e-14`, looks like an exponent typo.**
As printed, degree 14 is *worse* than degree 12 (`2.26e-07` versus `1.24e-08`),
the only non-monotone step across all five families. Changing `e-14` to `e-15`
drops the error to `3.28e-10` and restores monotonicity. `rotc4` is unused by the
header's rotation-estimate functions.

**7. Minor.** `Projection.h`'s single-plane `PerspectiveProject` discards the
`bool` from `Ellipse2::FromCoefficients`.

Issues [#225](https://github.com/gradientspaceai/gtengine-js/issues/225), [#374](https://github.com/gradientspaceai/gtengine-js/issues/374).

### `SeparatePoints2.h`, `SeparatePoints3.h`, `TriangulateCDT.h`

**1. Round-off in the side tests yields false-positive separations (both files).**
The side test compares `Dot(Perp(Normalize(P1 - P0)), P)` (2D) or the
`Plane3({P0,P1,P2})` signed distance (3D) against a rounded plane constant, so
the far endpoint of the hull's *own* candidate face can classify as strictly
positive by round-off; `WhichSide` then returns `+1` instead of `-1` and heavily
overlapping point sets are reported as separated.

**Reproductions.** 2D: hexagons centred `(0,0)` with radius 1 and phase 0.3, and
`(1,0.5)` with radius 1.5 and phase 1.1; the offending signed distance is
`+4.44e-16`. 3D: running the verbatim upstream floating-point algorithm against
an exact-arithmetic implementation over 400 randomized pairs of *overlapping*
point clouds, upstream reported a separation in **101 cases**, and in all 101 the
returned plane does not separate (for example `points0` signed distances span
`[-1.989, +eps]` while `points1` spans `[-1.777, -0.640]`). An independent
upstream-faithful implementation written later gave 70 false separations in 396
overlapping trials, so the rate is in the tens of percent either way.

**Suggested fix (applied in the port).** Use an exact orientation predicate on
the unnormalized normal, `Dot(Cross(P1-P0, P2-P0), Q-P0)` in 3D.

**2. `SeparatePoints3`'s cross-product-axis loop leaves `Plane3` inconsistent.**
It sets `normal` and `constant` directly while `origin` retains whichever face
plane was built last, the inconsistent state `Hyperplane.h`'s own comments warn
against. Port: constructs via normal plus constant.

**3. `TriangulateCDT::RemapPolygonTree` overwrites a remapping.**
`remapping[iter->second] = node.polygon[i]` overwrites the first occurrence's
remapping with a later duplicate's input index, contradicting the adjacent
comment, so the restored polygon reports an index the caller never passed.
Coordinates are identical, so output geometry is unaffected.
`ConstrainedTriangulate` also re-inserts triangles the graph copy already
contains.

**4. Minor.** `SeparatePoints2` classifies each hull vertex twice redundantly.

Issues [#328](https://github.com/gradientspaceai/gtengine-js/issues/328), [#348](https://github.com/gradientspaceai/gtengine-js/issues/348), [#405](https://github.com/gradientspaceai/gtengine-js/issues/405).

### `SingularValueDecomposition.h`, `SymmetricEigensolver.h`, `SymmetricEigensolver3x3.h`, `UnsymmetricEigenvalues.h`

**1. `SymmetricEigensolver::Tridiagonalize` stores a reflection parameter for a degenerate Householder step.**
When the subcolumn below the subdiagonal is already zero, `length == 0`, so
`mVVector` stays the zero vector (`v1 = 1` is inside the `if (length > 0)` guard)
and the reflection actually applied to the matrix is the identity. But `vdv`
remains 1, so the code still stores the reflection parameter:

```cpp
mMatrix[i + static_cast<size_t>(mSize) * ip1] = twoinvvdv;   // == 2
```

`GetEigenvectors()`/`GetEigenvector()` rebuild the reflection from that parameter
with the *implied* `v[i+1] = 1`, reconstructing `H = I - 2*e_{i+1}*e_{i+1}^T`
instead of the identity: a spurious sign flip of row `i+1` of `Q`. Eigenvalues
stay correct, because the stored value lands in the lower triangle, which the
diagonal and superdiagonal copy never reads; only `Q` is corrupted.

**Reproduction.** `A = diag({{2,1},{1,2}}, {{5,2},{2,5}})` gives correct
eigenvalues 1, 3, 3, 7 but returns `(0,0,-1,1)/sqrt(2)` for eigenvalue 7 instead
of `(0,0,1,1)/sqrt(2)`, with `max|Q*D*Q^T - A| = 4.0`. Triggers on
block-diagonal input, already-tridiagonal input, or an isolated eigenvalue.
Port: stores 0 so the rebuilt reflection is the identity; eigenvalues and the
tridiagonal matrix remain bit-for-bit identical to upstream.

**2. `SymmetricEigensolver3x3::GetCosSin` lacks the `maxAbsComp` rescaling its 2x2 sibling documents as necessary.**
For a covariance matrix whose off-diagonal entries are around 1e-160, as arise
from orthonormal-frame-style inputs with subnormal components, the squares
underflow. Measured: with `A = 1e160*B` the eigenvectors have length 0 and the
eigenvalues have no correct digits; with `A = 1e-170*B` the eigenvalues are wrong
by 36%. (The first observation was milder: eigenvalues wrong in the 7th
significant digit and eigenvectors of length 0.9999998833.) Port: rescales
`(u, v)` by `max(|u|, |v|)` only when that maximum leaves `[2^-511, 2^511]`, so
in-band inputs evaluate the upstream expression bit-identically; the error is
then at most 6e-16 everywhere. The 2x2 sibling is already scale-equivariant over
1e-300 to 1e300.

**3. `UnsymmetricEigenvalues` drops a trailing eigenvalue.**
The eigenvalue-packing loop iterates `i < mSizeM1` and never reports
`A(N-1,N-1)` when the final row decouples as a 1x1 block, so one real eigenvalue
is silently dropped for matrices that converge fully triangular. Port: adds a
guarded post-loop check replicating the 1x1 case for index `N-1`.

**4. `UnsymmetricEigenvalues::FrancisQRStep` has no exceptional shift, so the iteration can cycle.**
`A = [[0,8,0],[8,0,8],[0,8,0]]` has eigenvalues 0 and `+-8*sqrt(2)`, all real and
well separated; `Solve` burns its whole budget (16, 256 and 4096 iterations all
return the full limit) and reports none. The `Solve` documentation also claims it
returns `0xFFFFFFFF` on non-convergence; it returns `numIterations` on every path.

**5. `SingularValueDecomposition` sorts with unstable `std::sort`,** so the U/V
column permutation for tied singular values is unspecified. `SingularValueDecomposition`
has no analogue of item 1: it stores full Householder vectors and replays them
consistently. `SymmetricEigensolver::Solve`'s block scan does not have item 3's
trailing-block defect.

Issues [#42](https://github.com/gradientspaceai/gtengine-js/issues/42), [#80](https://github.com/gradientspaceai/gtengine-js/issues/80), [#379](https://github.com/gradientspaceai/gtengine-js/issues/379), [#476](https://github.com/gradientspaceai/gtengine-js/issues/476), [#478](https://github.com/gradientspaceai/gtengine-js/issues/478).

### `SurfaceExtractorMC.h`, `SurfaceExtractorTetrahedra.h`, `TetrahedraRasterizer.h`

**1. `SurfaceExtractorMC::Extract` omits `level` from the edge interpolation.**
Upstream writes `vertex[index] = F[j0] / (F[j0] - F[j1]);` where it must be
`(F[j0] - level) / (F[j0] - F[j1])`. This is correct only for `level == 0`;
otherwise vertices are misplaced and can leave the voxel. Port: shifts only the
numerators and keeps the unperturbed `F`, so `perturb` still affects only the
sign classification.

**2. The marching-cubes table is not face-consistent (preserved).** A shared face
with alternating corner signs is resolved differently by each voxel, leaving a
hole; the minimal 3x2x2 case produces 6 triangles with 8 boundary edges. This is
inherent to the 15-case table; the sibling extractors are unaffected. (The
256-entry articulation table was script-diffed and matches.)

**3. `SurfaceExtractorTetrahedra::GetGradient`'s odd-parity branch condition is unconditionally true.**
At about line 296:

```cpp
else if (dx + dy + dz >= (Real)0)   { /* tetra 6572 */ }
else                                { /* tetra 0752 */ }
```

`dx`, `dy`, `dz` are fractional offsets within the voxel, all in `[0,1]`, so
`dx + dy + dz >= 0` is unconditionally true: the 6572 corner tetrahedron absorbs
the central one and the 0752 branch is dead code. The plane through `(1,1,0)`,
`(1,0,1)`, `(0,1,1)` cutting off corner `(1,1,1)` is `dx + dy + dz = 2`; the
even-parity branch immediately below uses the analogous planes at 1 and is
correct. The wrong gradient corrupts `OrientTriangles` for any point of an
odd-parity voxel in the central tetrahedron. Port: condition changed to `>= 2`.

**4. `TetrahedraRasterizer::ClipCullAABBs` clips the stored boxes in place,** so a
second `operator()` with a larger region scans the box the first call left.

Issues [#443](https://github.com/gradientspaceai/gtengine-js/issues/443), [#132](https://github.com/gradientspaceai/gtengine-js/issues/132), [#439](https://github.com/gradientspaceai/gtengine-js/issues/439).

### `SWInterval.h`

Every bound is widened with `std::nextafter(value, -max)` or
`std::nextafter(value, +max)`. When the round-to-nearest computation of a bound
overflows to an infinity, `nextafter` steps it *back* toward `max`, so the bound
becomes `+-MAX_VALUE` and the interval no longer contains the exact result.
`SWInterval<double>::Div(8.881784197001252e-16, -DBL_TRUE_MIN)` returns
`[-MAX, -MAX]` although the exact quotient is about `-1.78e308`. The same
mechanism makes the reciprocal of an interval with a zero endpoint finite:
`[1,2] / [0,4]` returns an upper bound of `MAX_VALUE` rather than `+infinity`,
because `Mul` widens `2 * infinity` with `nextafter(infinity, +max)`.

`FPInterval.h` does not share this defect: it uses `fesetround`, and directed
rounding leaves an infinity alone.

Issues [#50](https://github.com/gradientspaceai/gtengine-js/issues/50), [#367](https://github.com/gradientspaceai/gtengine-js/issues/367). Port: preserved (the port's own interval type uses `getNextUp`/`getNextDown`, which preserve enclosure, rather than reusing these).

### `Tetrahedron3.h`

`GetPlanes` sets only `normal` and `constant`, leaving `Plane3::origin` at
`(0,0,0)`, so `Dot(normal, origin) != constant`: the returned planes violate the
`Hyperplane` invariant that the file's own comments warn about. Port: constructs
via normal plus constant (normal and constant values bit-identical to upstream).

Issue [#268](https://github.com/gradientspaceai/gtengine-js/issues/268).

### `Transform.h`

**1. `GetHInverse` can return a matrix with a corrupt last row.**
The RS branches write only the upper-left 3x3 block and last column of
`mInvHMatrix`, assuming the last row still holds `(0,0,0,1)`. But the non-RS
branch assigns `mInvHMatrix = Inverse(mHMatrix)` wholesale, and for a singular
`mHMatrix`, `Inverse` returns the zero matrix. A transform that goes through a
singular `SetHMatrix` or general state and later switches back to an RS transform
then reports an "inverse" whose last row is `(0,0,0,0)`. Port: fixed by writing
the last row explicitly.

**2. `Inverse()` stores the translation in the M channel.**
`inverse.SetMatrix(invMatrix)` passes the full affine 4x4 whose last column is
`-M^{-1} * T`, violating `GetMatrix`'s documented `{{M, 0}, {0, 1}}` block
structure. `GetHMatrix` and `operator*` happen to be immune, because they rebuild
from channels, but the `GetMatrix` accessor hands back a wrong matrix. Port:
fixed by zeroing the last column and row before storing.

**3-6 (minor and doc, preserved).** `SetRotation(AxisAngle<3>)` uses
`HLift(axis, 1)` where `AxisAngle<4>` is documented to need `(x, y, z, 0)`
(observationally harmless, since the matrix conversion ignores component 3);
`Invert3x3` is dead code and assigns `0.0f`/`1.0f` float literals into `Real`
entries, which would be lossy for exact-arithmetic instantiations; the doc
comment `"M = S*R (!GTE_USE_VEC_MAT)"` should read `GTE_USE_VEC_MAT`, since as
written both clauses name the same configuration; `SetRotation` never resets
`mIsUniformScale`, so after a `SetMatrix` a subsequent `SetRotation` reports
`IsUniformScale() == false` even with scale `(1,1,1)` (a hint only, both paths
compute correct results).

Issue [#265](https://github.com/gradientspaceai/gtengine-js/issues/265).

### `Vector.h`, `GVector.h`, `Vector4.h`

**1. Robust `Length` and `Normalize` return NaN for subnormal inputs.**
`Length(v, robust = true)` and `Normalize(v, robust = true)` return NaN when the
largest-magnitude component is subnormal, while the non-robust path returns the
correct (tiny or zero) value. For `v = (0, 5e-324, 0)` the robust length is NaN
and the plain length is 0. Cause: the robust path rescales with
`v /= maxAbsComp`, implemented as multiplication by `1 / maxAbsComp`; when
`maxAbsComp` is subnormal that reciprocal overflows to infinity, and the zero
components then produce `0 * inf = NaN`. Dividing each component directly by
`maxAbsComp` would avoid the overflow.

**2. `Vector4::ComputeOrthogonalComplement`'s `maxIndex == 3` branch can return the zero vector.**
The branch sets `v[1] = (0, +v[0][2], -v[0][1], 0)`, which is the zero vector
whenever components 1 and 2 vanish; for example `v[0] = (1,0,0,5)` yields no
basis at all and the function returns 0. The generic fallback described by the
branch's own comment would give `(0,0,-5,0)`. The branch's stated motivation, 3D
affine vectors with `w = 0`, can never trigger it, since `v[0][3] = 0` can only
be the max component for the zero vector.

Issues [#370](https://github.com/gradientspaceai/gtengine-js/issues/370), [#87](https://github.com/gradientspaceai/gtengine-js/issues/87). Port: preserved.

### `VertexCollapseMesh.h`

`Collapsed`'s comment claims the mesh is "restored" to its state before the
`Remove(...)` call; the check runs before any removal. `VCM_NO_MORE_ALLOWED` is
declared and documented but never returned, and the `VCM_DEFERRED` documentation
omits the nonmanifold-diagonal case. `DoCollapse`'s
`record.vertex = 0x80000000` relies on implementation-defined unsigned-to-int32
conversion pre-C++20.

**Invalid link triangulation corrupts the mesh (fixed in the port).**
`TriangulateLink` uses `TriangulateEC<Real, Real>`. For a projected link with
collinear vertices the floating-point ear clipping can return an invalid
triangulation: on a 6x6 grid that is flat except `z[1] = 0.5999999330626311`,
`z[4] = -0.10950932320884953`, the 16th collapse has the link
`0 1 2 3 4 11 17 23 29` and receives `<11,17,29>` twice with the edges `17-23`,
`23-29` uncovered. `Collapsed` removes the 9 old triangles, inserts 6, finds a
missing link edge and returns `VCM_UNEXPECTED_ERROR` with the mesh already
modified. The port checks, before removing anything, that the inserted triangles
are topologically a triangulation of the link (`n-2` triangles, link edges used
once, all other edges twice) and defers the vertex otherwise.

Issues [#295](https://github.com/gradientspaceai/gtengine-js/issues/295), [#412](https://github.com/gradientspaceai/gtengine-js/issues/412), [#498](https://github.com/gradientspaceai/gtengine-js/issues/498).
### Remaining single-item findings

- `AlignedBoxBV.h`: `GetSplittingAxis` has a dead store `maxExtent = extents[2]`.
- `ApprOrthogonalLine3.h`: the header comment says "the minimum eigenvalue is
  unique"; the code correctly tests the maximum.
- `BSplineSurface.h`: lines 10-18 carry a stale copy-pasted header comment
  describing `BSplineReduction` (curve control-point reduction, L2 integral
  norm), which has nothing to do with this class.
- `Polygon2.h`: the constructor cannot reject an empty vertex pool.
- `UIntegerALU32.h`: `RoundUp`, as used by `BSNumber::Convert`, returns the
  trailing-zero shift count, not the 0/1 carry the surrounding comments imply.
  The behaviour is correct; the comment is misleading.

Issues [#96](https://github.com/gradientspaceai/gtengine-js/issues/96), [#268](https://github.com/gradientspaceai/gtengine-js/issues/268), [#380](https://github.com/gradientspaceai/gtengine-js/issues/380).

## Systematic patterns

The findings are not independent. Twelve recurring classes account for most of
them, and each is worth a targeted sweep of the library rather than a per-file
fix.

### 1. Exact floating-point equality or zero used as a degeneracy gate

By far the largest class. A predicate that is exactly true in exact arithmetic is
tested with `==`, `>= 1`, `<= 0` or an epsilon that defaults to 0, and rounding
then routes control into a branch written for the non-degenerate case.

Examples: `Delaunay2`/`Delaunay3`'s hardcoded `epsilon = 0` in the intrinsics
classification ([#391](https://github.com/gradientspaceai/gtengine-js/issues/391));
`IntrLine3Cone3`'s through-vertex test and its `discr` sign
([#465](https://github.com/gradientspaceai/gtengine-js/issues/465));
`IntrPlane3Plane3`'s `|Dot(N0,N1)| >= 1`, which only about half of normalized
doubles satisfy ([#465](https://github.com/gradientspaceai/gtengine-js/issues/465));
`Rotation`'s Euler gimbal-lock test `r(i,j) == +-1`
([#374](https://github.com/gradientspaceai/gtengine-js/issues/374));
`Minimize1`'s `tv == tm` at equal-value bracket endpoints
([#298](https://github.com/gradientspaceai/gtengine-js/issues/298));
`Polynomial1::GreatestCommonDivisor`'s exact-zero remainder test
([#83](https://github.com/gradientspaceai/gtengine-js/issues/83)) and
`CubicRootsQR`/`QuarticRootsQR`'s `discriminant >= 0`
([#488](https://github.com/gradientspaceai/gtengine-js/issues/488));
`ApprCurveByArcs`'s `fabs(det) >= epsilon` with `epsilon` defaulting to 0
([#163](https://github.com/gradientspaceai/gtengine-js/issues/163));
`DistCircle3Circle3`'s `normal[2] < 1` and its perfect-square `phi`
([#431](https://github.com/gradientspaceai/gtengine-js/issues/431),
[#442](https://github.com/gradientspaceai/gtengine-js/issues/442));
`MinimumAreaCircle2`'s reliance on `LinearSystem::Solve`'s `det != 0`
([#399](https://github.com/gradientspaceai/gtengine-js/issues/399));
`GenerateMeshUV`'s untoleranced `lower_bound` on accumulated arc length
([#262](https://github.com/gradientspaceai/gtengine-js/issues/262));
`ContEllipse2`/`ContEllipsoid3`'s `Length <= 1` containment test with no
tolerance ([#409](https://github.com/gradientspaceai/gtengine-js/issues/409));
the `>= cutoff` versus `> cutoff` divergence between the aligned and oriented
box-box queries ([#450](https://github.com/gradientspaceai/gtengine-js/issues/450));
`SeparatePoints2`/`SeparatePoints3`'s side tests against a rounded plane constant
([#328](https://github.com/gradientspaceai/gtengine-js/issues/328),
[#348](https://github.com/gradientspaceai/gtengine-js/issues/348)).

Where the library already has an exact predicate (`PrimalQuery2::ToLine`,
`ToPlane`, `BSNumber` arithmetic), routing the classification through it fixes
the problem without changing any non-degenerate result. That is what the port did
for `Delaunay2/3` and `SeparatePoints2/3`.

### 2. First-wins candidate selection where all candidates must be compared

A search accepts the first candidate that succeeds, although each candidate only
bounds the answer from one side.

`IntrAlignedBox3Sphere3::DoQueryRayRoundedFace` accepts the first rounded-edge
probe, but any piece of the Minkowski sum gives only an upper bound on the
contact time, so contacts are reported late or missed
([#458](https://github.com/gradientspaceai/gtengine-js/issues/458),
[#465](https://github.com/gradientspaceai/gtengine-js/issues/465)).
`MinimumVolumeBox3FloatingPoint::GetExtreme`'s strict-improvement hill climb
stops at the first plateau ([#426](https://github.com/gradientspaceai/gtengine-js/issues/426)),
and `ComputeVolume` assumes a hull-edge vertex realizes the axis minima
([#405](https://github.com/gradientspaceai/gtengine-js/issues/405)).
`RotatingCalipers::CreatePolygon` tests collinearity against the immediately
preceding edge rather than the most recent nonzero one
([#286](https://github.com/gradientspaceai/gtengine-js/issues/286)).
`BVTreeOfTriangles::Execute` keeps one hit per parameter value because the
`std::set` key is a single field
([#167](https://github.com/gradientspaceai/gtengine-js/issues/167)).
`ConvexHull2::GetTangent` returns whatever indices it last held when its bounding
loop expires ([#277](https://github.com/gradientspaceai/gtengine-js/issues/277)).

### 3. `Vector::operator/` zero-divisor semantics: a hazard for translators

GTE's `Vector::operator/` and `operator/=` do not divide componentwise; they
multiply by the reciprocal, and `Vector.h` special-cases a zero divisor to yield
the zero vector. Two consequences:

- **For upstream.** The reciprocal formulation overflows for a subnormal divisor.
  `Length(v, robust = true)` and `Normalize(v, robust = true)` rescale with
  `v /= maxAbsComp`, so a subnormal `maxAbsComp` gives `1/maxAbsComp = infinity`
  and the zero components become `0 * inf = NaN`:
  `Length((0, 5e-324, 0), true)` is NaN while the plain length is 0
  ([#370](https://github.com/gradientspaceai/gtengine-js/issues/370)). Dividing
  componentwise would avoid it. The same reciprocal idiom underlies
  `GaussianElimination`'s denormal-pivot NaN
  ([#375](https://github.com/gradientspaceai/gtengine-js/issues/375)).
- **For anyone re-implementing GTE.** A naive componentwise translation of
  `v / d` changes behaviour at `d == 0` from "zero vector" to "NaN vector", and
  the difference only shows up on degenerate geometry. This was the single most
  frequent translation defect found in the port's own verification pass (eight
  separate sites, in `MeshSmoother`, `TCBSplineCurve`, `SampleCircularArc`,
  `NURBSSphere`, `IntrSphere3Triangle3`, `IntrAlignedBox3Cone3`,
  `IntrConvexMesh3Plane3` and `GenerateMeshUV`). It is worth a sentence in
  `Vector.h`'s own documentation.

### 4. Published error and precision tables that do not match the code

`RotationEstimate.h`'s `C_ROTC*_EST_MAX_ERROR` tables understate the true maximum
error for `rotc0`/`rotc1`/`rotc2` at degrees 14 and 16, for `rotc3` at degrees
12, 14 and 16, and for the entire `rotc4` table, by factors from 3 to 1.5e12;
several `rotc4` entries are below double representability for the quantity they
bound ([#225](https://github.com/gradientspaceai/gtengine-js/issues/225)).
`LogEstimate::GetLogEstimateMaxError` returns the log2 bound, loose by a factor
of about 1.443 ([#57](https://github.com/gradientspaceai/gtengine-js/issues/57)).
`ASinEstimate` and `ACosEstimate` publish different values for a provably
identical quantity ([#57](https://github.com/gradientspaceai/gtengine-js/issues/57)).
The `PrimalQuery2.h`/`PrimalQuery3.h` N-value tables were computed with a removed
`BSPrecision` API ([#43](https://github.com/gradientspaceai/gtengine-js/issues/43)),
and `BSPrecision::operator+` itself under-estimates `maxBits` for mixed operand
sets, which is exactly how those chains are built
([#366](https://github.com/gradientspaceai/gtengine-js/issues/366)).

Every other published estimate bound in the library holds when re-measured on a
dense grid, so the problem is confined to these tables.

### 5. Members that are never initialized or never assigned

`IntrAreaEllipse2Ellipse2` declares `T mZero, mOne, mTwo, mPi, mTwoPi;` with no
constructor and never assigns them, so every computed area reads indeterminate
values ([#301](https://github.com/gradientspaceai/gtengine-js/issues/301)).
`TCBSplineCurve` and `PolynomialCurve` never set `mConstructed`, so
`operator bool` reports failure for every valid curve
([#182](https://github.com/gradientspaceai/gtengine-js/issues/182)).
`MinimumAreaBox2`'s convex-polygon overload never assigns `mNumPoints`/`mPoints`
([#402](https://github.com/gradientspaceai/gtengine-js/issues/402)).
`RootsBisection1` returns its output references unassigned when
`maxIterations == 1` ([#84](https://github.com/gradientspaceai/gtengine-js/issues/84)).
`ApprEllipse2` and `ApprEllipsoid3` read a default-constructed box after ignoring
`GetContainer`'s failure flag
([#224](https://github.com/gradientspaceai/gtengine-js/issues/224),
[#322](https://github.com/gradientspaceai/gtengine-js/issues/322)).
`IntrPlane3Triangle3`'s FIQuery uses a raw uninitialized `Real s[3]`
([#255](https://github.com/gradientspaceai/gtengine-js/issues/255)).

### 6. `size_t` wrap-around and unsigned underflow

`PolylineOffset`'s member initializer computes `vertices.size() - 1` for an empty
open polyline, so allocation throws before the intended diagnostic
([#112](https://github.com/gradientspaceai/gtengine-js/issues/112)).
`Image2`/`Image3`'s neighbourhood getters add a negative `int32_t` offset to a
`size_t` coordinate and wrap to `SIZE_MAX` at the boundary
([#64](https://github.com/gradientspaceai/gtengine-js/issues/64)).
`LinearSystem::SolveTridiagonal` allocates `N - 1` entries without validating
`N` ([#261](https://github.com/gradientspaceai/gtengine-js/issues/261)).
`MinimumSpanningTree` writes `records[0]` on a zero-length vector
([#74](https://github.com/gradientspaceai/gtengine-js/issues/74)).
`FastMarch3` computes `i - mXYBound` for a `z = 0` voxel
([#121](https://github.com/gradientspaceai/gtengine-js/issues/121)).
`ApprCone3ExtractEllipses` indexes with `size_t(-1)` when the plane set is empty
([#349](https://github.com/gradientspaceai/gtengine-js/issues/349)).

### 7. Status returns that are discarded

A `bool` or enum saying "I failed" is ignored and the value-initialized output is
used as if valid. `IntpQuadraticNonuniform2` discards *every* mesh-accessor flag,
so zero index triples silently blend sample 0
([#337](https://github.com/gradientspaceai/gtengine-js/issues/337));
`IntpLinearNonuniform2/3` discard `GetIndices`
([#135](https://github.com/gradientspaceai/gtengine-js/issues/135));
`ApprEllipseByArcs` discards `Circumscribe`, storing the previous arc
([#322](https://github.com/gradientspaceai/gtengine-js/issues/322));
`ApprCone3` discards `ApprHeightLine2::Fit`
([#271](https://github.com/gradientspaceai/gtengine-js/issues/271));
`AdaptiveSkeletonClimbing3`'s root call discards `Merge`'s "I merged completely"
signal, producing an empty mesh
([#194](https://github.com/gradientspaceai/gtengine-js/issues/194));
`ExtremalQuery3BSP` discards `VETManifoldMesh::Insert` failure
([#290](https://github.com/gradientspaceai/gtengine-js/issues/290));
`Projection` discards `Ellipse2::FromCoefficients`
([#225](https://github.com/gradientspaceai/gtengine-js/issues/225)).

### 8. State that survives into the next query on a reused object

`ConstrainedDelaunay2` never clears `mInsertedEdges`
([#325](https://github.com/gradientspaceai/gtengine-js/issues/325));
`IntrAlignedBox3Cone3` never clears `mAdjacencyMatrix`
([#301](https://github.com/gradientspaceai/gtengine-js/issues/301));
`MinimumAreaBox2`'s degenerate branches return the previous query's `mArea` and
`mSupportIndices` ([#328](https://github.com/gradientspaceai/gtengine-js/issues/328));
`TetrahedraRasterizer::ClipCullAABBs` clips the stored boxes in place
([#439](https://github.com/gradientspaceai/gtengine-js/issues/439));
`MinimizeN` never resets `mDConjIndex` in `GetMinimum`
([#146](https://github.com/gradientspaceai/gtengine-js/issues/146));
`RootsBisection2` reports stale roots from a previous call
([#84](https://github.com/gradientspaceai/gtengine-js/issues/84),
[#152](https://github.com/gradientspaceai/gtengine-js/issues/152)).

### 9. Copy-paste between dimensional or structural siblings

A 2D routine is lifted to 3D, or one tree/mesh class is lifted from another, and
one index, stride, sign or bound is not adjusted.

`IntpAkimaUniform2`/`3` reuse min-boundary stencil coefficients at the max
boundaries ([#58](https://github.com/gradientspaceai/gtengine-js/issues/58));
`MeshStaticManifold3::GetAdjacentTetrahedra` uses the 2D tuple indices on a
5-wide record ([#66](https://github.com/gradientspaceai/gtengine-js/issues/66));
`ImageUtility3::Dilate` starts at `i0 = 1` where the 2D version starts at 0, and
`Close<N>` still carries the 2D `static_assert`
([#129](https://github.com/gradientspaceai/gtengine-js/issues/129));
`Mesh::ComputeIndices` and `TubeMesh::UpdatePositions` use the `numCols` stride
where rows are `numCols + 1` wide
([#220](https://github.com/gradientspaceai/gtengine-js/issues/220),
[#240](https://github.com/gradientspaceai/gtengine-js/issues/240));
`OrientedBoxTreeOfTriangles::ComputeLeafBoundingVolume` keeps the maximum scan's
comparison while looking for the minimum
([#343](https://github.com/gradientspaceai/gtengine-js/issues/343));
both `MinimumVolumeBox3` variants drop the same Newell wrap-around term
([#352](https://github.com/gradientspaceai/gtengine-js/issues/352),
[#355](https://github.com/gradientspaceai/gtengine-js/issues/355));
the run-time `BlockCholeskyDecomposition` uses the block stride for the in-block
offset ([#209](https://github.com/gradientspaceai/gtengine-js/issues/209));
`CurvatureFlow2` uses 0.5 where `CurvatureFlow3` uses 2
([#123](https://github.com/gradientspaceai/gtengine-js/issues/123));
`ApprPolynomialSpecial3/4` apply the 1-D degree check per axis
([#163](https://github.com/gradientspaceai/gtengine-js/issues/163));
`IntrSegment2OrientedBox2` keeps `IntrSegment2AlignedBox2`'s world-space point
formula ([#255](https://github.com/gradientspaceai/gtengine-js/issues/255)).

In almost every case the *other* sibling is correct, which is what makes these
easy to confirm and cheap to fix.

### 10. Sentinel and contract skew between a primitive and its queries

`Cylinder3` documents and implements the infinite sentinel as `height = -1`, but
`DistPoint3Cylinder3` tests `height == numeric_limits<T>::max()`, so the infinite
path is unreachable and the finite path throws
([#187](https://github.com/gradientspaceai/gtengine-js/issues/187)).
`IntrHalfspace3Cylinder3`, `IntrCylinder3Cylinder3`, `IntrRay3Cylinder3`,
`IntrSegment3Cylinder3` and `IntrTriangle3Cylinder3` have no finiteness check at
all, so `-1` becomes a negative half-height
([#197](https://github.com/gradientspaceai/gtengine-js/issues/197),
[#206](https://github.com/gradientspaceai/gtengine-js/issues/206),
[#255](https://github.com/gradientspaceai/gtengine-js/issues/255)).
`DistOrientedBox3Cone3` never validates `Cone`'s finite-`hmax` precondition
([#298](https://github.com/gradientspaceai/gtengine-js/issues/298)), and
`IntrAlignedBox3Cone3`/`IntrOrientedBox3Cone3` still carry a TODO for the retired
infinite-cone representation
([#334](https://github.com/gradientspaceai/gtengine-js/issues/334)). The
file-version headers make the skew visible: `Cylinder3.h` is `8.0.2025.05.10`,
`DistPoint3Cylinder3.h` is `8.0.2026.08.08`.

Related documentation drift: four interpolators claim to clamp their inputs and
clamp only the cell index
([#69](https://github.com/gradientspaceai/gtengine-js/issues/69));
`ParametricCurve::GetTime` ignores the segment decomposition `GetLength` honors
([#113](https://github.com/gradientspaceai/gtengine-js/issues/113));
`BasisFunction::GetValue` promises zeros it does not write
([#415](https://github.com/gradientspaceai/gtengine-js/issues/415)).

### 11. Cancellation-prone accumulations used as gates

`DistLine3CanonicalBox3` accumulates `sqrDistance` incrementally, which can go
negative for a flat box and make `sqrt` return NaN, and loses about half the
mantissa for grazing lines
([#421](https://github.com/gradientspaceai/gtengine-js/issues/421)).
`IntrLine3Cone3`'s `discr = c1*c1 - c0*c2` is exactly zero for a through-vertex
line, so in floating point it is a difference with no significant digits and its
*sign* selects the case analysis
([#465](https://github.com/gradientspaceai/gtengine-js/issues/465)).
`DistCircle3Circle3`'s `phi = p6^2 - (1 - cs^2) p7^2` splits double roots
([#331](https://github.com/gradientspaceai/gtengine-js/issues/331)).
`DistPointHyperellipsoid::Bisector`'s bracket rounds to exactly -1 and the root
then cancels catastrophically
([#424](https://github.com/gradientspaceai/gtengine-js/issues/424)).
`MinimumVolumeBox3FloatingPoint`'s endpoint samples compute `q0` and `q1` as
cancellation noise, giving a random axis direction
([#405](https://github.com/gradientspaceai/gtengine-js/issues/405)).
`SymmetricEigensolver3x3::GetCosSin` squares un-rescaled entries and underflows
([#379](https://github.com/gradientspaceai/gtengine-js/issues/379)).

### 12. Container misuse: `operator[]`, `insert`, and weak orderings

`std::map::operator[]` is used for lookup, silently inserting a default entry:
`ETManifoldMesh::GetBoundaryPolygon` (while `GetBoundaryPolygons` iterates the
same map, [#212](https://github.com/gradientspaceai/gtengine-js/issues/212)),
`SplitMeshByPlane::SplitTriangle*`
([#310](https://github.com/gradientspaceai/gtengine-js/issues/310)),
`IntrConvexMesh3Plane3`
([#301](https://github.com/gradientspaceai/gtengine-js/issues/301)),
`ExtremalQuery3BSP` ([#290](https://github.com/gradientspaceai/gtengine-js/issues/290)).
`find()` results are dereferenced without an `end()` check in `PlanarMesh`
([#256](https://github.com/gradientspaceai/gtengine-js/issues/256)) and
`IntrConvexMesh3Plane3`. `std::map::insert` is used where assignment was meant,
so `BSPPolygon2::SplitEdge` silently no-ops and leaves an edge unmapped
([#169](https://github.com/gradientspaceai/gtengine-js/issues/169)). And three
comparators are not strict weak orderings, which is undefined behaviour for the
containers and algorithms that consume them: `GMatrix`'s mixed-dimension
comparisons ([#88](https://github.com/gradientspaceai/gtengine-js/issues/88)),
`SortPointsOnCircle::LessThanByGeometry` when a point equals the sort centre
([#394](https://github.com/gradientspaceai/gtengine-js/issues/394)), and
`BVTreeOfTriangles::Intersection::operator<`, which makes distinct hits
equivalent ([#167](https://github.com/gradientspaceai/gtengine-js/issues/167)).

A related family is partial mutation on a failure path: `VEManifoldMesh`,
`ETManifoldMesh`, `TSManifoldMesh` and `ETNonmanifoldMesh` all leave edges or
faces in their maps referencing a feature that was never added, or that has been
destroyed by exception unwinding
([#73](https://github.com/gradientspaceai/gtengine-js/issues/73),
[#179](https://github.com/gradientspaceai/gtengine-js/issues/179),
[#212](https://github.com/gradientspaceai/gtengine-js/issues/212),
[#472](https://github.com/gradientspaceai/gtengine-js/issues/472)).

## Claims withdrawn or corrected

The verification pass re-derived every porting-pass claim. These were changed;
the corrected form is what appears above.

1. **`Rectangle::GetVertices` is not "genuinely counterclockwise".** Issue
   [#155](https://github.com/gradientspaceai/gtengine-js/issues/155) item 3
   originally said that `Parallelepiped3`/`Parallelogram2` emit bit-pattern order
   "unlike `Rectangle::GetVertices`, which is genuinely CCW". That parenthetical
   is wrong: `Rectangle` emits the same bit-pattern order
   (`vertex[0] = C - sum`, `[1] = C + dif`, `[2] = C - dif`, `[3] = C + sum`).
   Only the `Parallel*` comments are wrong; `Rectangle`'s comment is correct.

2. **`IntrSphere3Sphere3` does not share `IntrCircle2Circle2`'s TI/FI split.**
   Issue [#450](https://github.com/gradientspaceai/gtengine-js/issues/450) item 2
   originally named `IntrSphere3Sphere3.h` as a sibling instance. It is not:
   `IntrSphere3Sphere3`'s FI query reports containment as `intersect = true`
   (types 3 to 6), so its TI and FI agree everywhere. The circle case stands, and
   `IntrRay2Ray2` (collinear opposite rays, `numIntersections` 1 versus 2) is the
   other genuine TI/FI divergence.

3. **`IntrTriangle3Cylinder3::DiskOverlapsPolygon` is not "harmless for valid input".**
   Issue [#206](https://github.com/gradientspaceai/gtengine-js/issues/206)
   item 1 described the point-degenerate polygon case as harmless because a
   nondegenerate triangle cannot project to a point. That is true but incomplete:
   any nondegenerate triangle whose plane contains the cylinder axis projects to
   a *segment through the origin*, where both sign counts are also zero, and is
   reported as intersecting at any distance. The corrected characterisation is
   result-corrupting; see
   [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) item 2.

4. **`IntrEllipsoid3Ellipsoid3::GetRoots`'s asserts are not merely latent.**
   Issue [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) item 6
   recorded the ad-hoc `epsilon = 0.001` bracketing as a "robustness suspect"
   that randomized testing never tripped. It is reachable for ordinary
   close-centre input, erratically in the offset (0.002 and 0.0001 classify;
   0.001 and 1e-6 throw "Unexpected condition."); see
   [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) item 4.

5. **The `OBBTree` family's `pmin = pmax = 0` seed is not merely "conservative".**
   Issue [#103](https://github.com/gradientspaceai/gtengine-js/issues/103)
   item 2 called the resulting boxes conservatively large. Sharpened: it is a
   *correct* min/max. `pmin <= 0 <= pmax` is a loop invariant and the projections
   always straddle zero, because the projection origin is the mean of the node's
   centroids; the only consequence is that the box is forced to contain that
   mean.

6. **`OrientedBoxTreeOfTriangles` zeroes the largest extent every time, not sometimes.**
   Issue [#343](https://github.com/gradientspaceai/gtengine-js/issues/343)
   item 1 said "whenever `extent[2]` is not the smallest". Sharpened:
   `GetContainer` returns axes ordered by increasing variance, so `extent[2]` is
   the largest for essentially every triangle.

7. **`IntpAkimaUniform2`'s boundary defect is a first-order error, not a
   cosmetic one.** Issue [#58](https://github.com/gradientspaceai/gtengine-js/issues/58)
   originally recorded it as affecting only cross-derivative estimates, with
   sample pass-through and C1 continuity unaffected, and the port preserved it.
   Corrected: upstream does not reproduce a bilinear function in the boundary
   cells at all (`f(1.5, 0.75) = 1.1015625` instead of `1.125`). The port now
   fixes it, in both the 2D and 3D headers.

8. **`BVTreeOfTriangles`'s coincident-hit drop is fixed, not preserved, in the port.**
   Issue [#167](https://github.com/gradientspaceai/gtengine-js/issues/167)'s
   "Port status: preserved" line is stale: `OBBTreeOfTriangles` was fixed first
   and `BVTreeOfTriangles` was aligned on it, so both ports now use lexicographic
   `(parameter, triangleIndex)` ordering. The upstream defect is unchanged.

9. **Measurements sharpened, not withdrawn.** `ExtremalQuery3BSP`'s icosahedron
   error rate was first reported as "~1.2% of directions"; a larger run gives
   49/5000 = 0.98% for the regular icosahedron and 188/5000 = 3.76% for its
   subdivision ([#290](https://github.com/gradientspaceai/gtengine-js/issues/290)).
   `SeparatePoints3`'s false-separation count was 101 of 400 in the first
   measurement and 70 of 396 in an independent re-implementation
   ([#348](https://github.com/gradientspaceai/gtengine-js/issues/348)); both
   measure the same upstream algorithm on different random samples.

Two items were investigated and found **not** to be defects, and are recorded
here so they are not re-reported: `ApprEllipseByArcs::UpdateMatrix` uses `a[i]`
on the diagonal versus `2*a[i]` off-diagonal and drops the global factor 2 of
`dF/dM`, which is a deliberate symmetric pairing absorbed by the subsequent
normalization ([#322](https://github.com/gradientspaceai/gtengine-js/issues/322));
and `IntpThinPlateSpline`'s lack of a rank test for singular sample
configurations is a design limit rather than a coding error
([#191](https://github.com/gradientspaceai/gtengine-js/issues/191)).
