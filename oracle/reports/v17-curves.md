# Group 17 (curves) — `v17-curves`

51 cases, 20 records each in `oracle/golden/v17-curves.txt`; 47 declared
`exact` (bit-identity required), 3 compared with the default scaled tolerance
because upstream calls `std::pow`, 1 declared a deliberate `deviation`.

Deep run: `npm run oracle:deep -- 2000 v17-curves` — 102000 records
(3817 of them recorded C++ throws), **all 52 tests pass**, wall time 4.7 s
(generation plus replay). Over the deep run **every real output of every
`exact` case was bit-identical**, and so were the three tolerance cases: the
statistics report 0 inexact real outputs in the whole family.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `BasisFunction.h` | `BasisFunction.create`, `.evaluate`, `.evaluate.repeated`, `.getValue.invalid`, `BasisFunctionInput.openUniform` | exact | pass |
| `ParametricCurve.h` | `ParametricCurve.differentialGeometry`, `.setTimeInterval`, `.getLength`, `.getLength.multiSegment`, `.getTime`, `.subdivideByTime`, `.subdivideByLength` | exact | pass |
| `BezierCurve.h` | `BezierCurve.evaluate`, `.evaluate.invalidOrder` | exact | pass |
| `NURBSCurve.h` | `NURBSCurve.evaluate`, `.evaluate.invalidOrder`, `.accessors`, `.deferred` | exact | pass |
| `NURBSCircle.h` | `NURBSCircle.quarterDegree2`, `.quarterDegree4`, `.halfDegree3`, `.fullDegree3`, `.arcDegree2` | exact | pass |
| `NURBSSphere.h` | `NURBSSphere.eighthDegree4`, `.halfDegree3`, `.fullDegree3` | exact | pass |
| `DarbouxFrame.h` | `DarbouxFrame.compute`, `.getPrincipalInformation` | exact | pass |
| `FrenetFrame.h` | `FrenetFrame2.compute`, `FrenetFrame3.compute` | exact | pass |
| `FrenetFrame.h` | `FrenetFrame2.getCurvature`, `FrenetFrame3.getCurvatureAndTorsion` | tolerance (`std::pow(speedSqr, 1.5)`); the torsion output is held to bit-identity with `outRealExact` | pass |
| `ImplicitCurve2.h` | `ImplicitCurve2.evaluate` | exact | pass |
| `ImplicitCurve2.h` | `ImplicitCurve2.getCurvature` | tolerance (`std::pow(fx^2+fy^2, 1.5)`); the `valid` flag compares exactly | pass |
| `NaturalCubicSpline.h` | `NaturalCubicSpline.free`, `.closed`, `.clamped` | exact | pass |
| `NaturalQuinticSpline.h` | `NaturalQuinticSpline.free`, `.closed`, `.clamped` | exact | pass (after the fix below) |
| `TCBSplineCurve.h` | `TCBSplineCurve.evaluate`, `.degenerateLambda` | exact | pass |
| `TCBSplineCurve.h` | `TCBSplineCurve.isConstructed` | deviation | 2000/2000 records deviate |
| `ReparameterizeByArclength.h` | `ReparameterizeByArclength.getT` | exact, including the iteration count | pass |
| `SampleCircularArc.h` | `SampleCircularArc.compute`, `.compute.largeArc` | exact (see the `std::acos` note) | pass |
| `PolylineOffset.h` | `PolylineOffset.execute`, `.throwParity` | exact | pass |
| `CLODPolyline.h` | `CLODPolyline.construct`, `.setLevelOfDetail` | exact | pass |
| `MassSpringCurve.h` | `MassSpringCurve.update` (1–3 RK4 steps) | exact | pass |

Generators. Bases are drawn in five modes — open uniform (exactly as
`BasisFunctionInput(n, d)` builds it), open nonuniform with dyadic knots, open
with a repeated interior knot, periodic uniform and periodic nonuniform — with
degrees 1..5; the deep run splits 371/797/412/420 over
open-uniform / open-nonuniform / periodic-uniform / periodic-nonuniform.
Evaluation parameters visit `tmin`, `tmax`, an exact interior knot (including a
repeated one), dyadic interior points, parameters just outside the domain
(clamping) and parameters several periods away (wrapping). Control points
alternate a uniform cloud, a small integer lattice, a lattice with a deliberately
repeated point and an all-coincident set (zero-speed curves, which are the
degenerate branch of `GetTangent`, `FrenetFrame*` and `GetCurvature`).
Weights are drawn as small integers, dyadic fractions or uniform reals; times
have dyadic differences. Polylines are drawn uniform, on a lattice, as a convex
lattice polygon and with an explicit collinear run.

`std::acos` note. `SampleCircularArc` truncates `radius * angle / k` to a
sample count, and `angle` comes from `std::acos`. MSVC and V8 may differ by one
ulp there, which would change the count only at an integer boundary, so the
generator rejects draws whose value has a fractional part outside `[0.15, 0.85]`
(and bounds the value to `[1, 8]`, keeping the records small). Everything else in
the sampler is `+ - * /` and `sqrt`, so the case is exact; the deep run produced
sample counts 1..14 for arcs of at most pi and 3..28 for larger arcs, with no
count disagreement in 4000 records.

## Port defects fixed

**`src/NaturalQuinticSpline.ts`, `backSubstitute`: accumulation order of the
four-term back-substitution sum.** Upstream writes

```cpp
poly[2] = invR(0,0)*B[j0] + invR(0,1)*B[j1] + invR(0,2)*B[j2] + invR(0,3)*B[j3];
```

which C++ accumulates strictly left to right, `((a+b)+c)+d`. The port grouped
the terms pairwise, `(a+b)+(c+d)`. The two differ in the last bits whenever
`B[j2]` and `B[j3]` are both nonzero, which is the case for every closed and
every clamped spline; a *free* spline has `B[j2] = 0` exactly (its `boundary0`
is the zero vector and the `ell20` row reduction does not touch it), which is
why the defect was invisible there and why the port's own property tests never
caught it. The error is amplified by the backward recursion over the segments:
`NaturalQuinticSpline.closed` and `.clamped` disagreed on **20 of 20** golden
records with scaled errors up to **2.3e-13** on the polynomial coefficients and
the jets. After the fix all 2000 deep records of all three quintic cases are
bit-identical. Regression test:
`test/NaturalQuinticSplineAssociation.test.ts`, which pins the MSVC values for a
closed and a clamped spline bit for bit and records the pre-fix (pairwise)
values to show that the two groupings really are different computations on that
input.

**Faithfulness alignment (no observable change on valid inputs).**
`NaturalCubicSpline.evaluate` and `NaturalQuinticSpline.evaluate` scaled the
derivative jets with `mul(v, 1 / denom)`. Upstream writes `v / denom`, i.e.
`operator/(Vector, Real)`, which multiplies by the reciprocal *and yields the
zero vector when the scalar is zero*. For `denom != 0` the two are bit-identical
(which is why the oracle cannot distinguish them — the spline contract requires
strictly increasing times), but a zero-length segment would give infinities
instead of upstream's zero vector. Both now use `div()`, matching
`TCBSplineCurve`, `SampleCircularArc`, `NURBSCircle` and `NURBSSphere`, which
already had the right semantics.

Every other vector-by-scalar division on these paths was checked against the
exact upstream operator and was already correct: the NURBS homogeneous divides
(`NURBSCurve`, `NURBSSurface`, `NURBSSphere`) compute `invW = 1/w` and multiply,
exactly as upstream does; `NURBSCircle`, `SampleCircularArc` and
`NURBSEighthSphereDegree4` use the `operator/` semantics; `TCBSplineCurve`
already used `div()`.

## Deliberate deviations demonstrated

`TCBSplineCurve.isConstructed` (issue
[#182](https://github.com/gradientspaceai/gtengine-js/issues/182)). Upstream's
`TCBSplineCurve` constructor never sets `mConstructed`, so `operator bool()`
reports failure for every valid curve; the port sets it. All 2000 deep records
deviate. The curve's values are unaffected, which the neighbouring
`TCBSplineCurve.evaluate` case (exact, 2000/2000 bit-identical) shows.

Preserved upstream defects, compared bit for bit rather than deviated from:

* `BasisFunction::GetIndex`'s endpoint shortcuts and the stale values left in
  the jet by `Evaluate` (issue #415) — `BasisFunction.evaluate.repeated`
  evaluates twice on one object and emits the whole jet, including the stale
  entries outside the current support.
* `TCBSplineCurve::ComputeInteriorTangents`'s division by zero in the lambda
  pass (issue #415) — `TCBSplineCurve.degenerateLambda` constructs
  `P[k-1] == P[k+1]` with uniform times and zero TCB parameters and emits the
  NaN *pattern* as booleans, because the harness treats any NaN as equal to any
  NaN.
* `SampleCircularArc::SampleArc3`/`SampleArc4`'s split directions for arcs
  larger than pi (issue #183) — `SampleCircularArc.compute.largeArc`.
* `ReparameterizeByArclength::DoNewtonsMethod`'s fall-through to
  `tMid - fMid/dfdt` after a `dfdt == 0` bisection step (issue #183) — reached
  by `ReparameterizeByArclength.getT` with `useBisection = false`.
* `ParametricCurve::GetTime`'s single Romberg integration across knots (issue
  #113) — `ParametricCurve.getLength.multiSegment` compares `GetTime` on a
  multi-segment natural spline bit for bit.
* `PolylineOffset`'s constructor validation order (issue #112): the port
  validates before sizing, upstream sizes first. Both throw on the same inputs,
  so `PolylineOffset.throwParity` is a normal throw-parity case (1817 of 2000
  deep records are throws) and not a deviation.

## Independent-reference checks

Run once against the port (which the oracle shows is bit-identical to upstream
on these paths), so the numbers below are upstream's as well:

| check | worst residual |
| --- | --- |
| basis functions nonnegative and summing to 1 (400 open bases) | `5.6e-16`, no negative value |
| jet[1] versus a central difference (Bezier / NURBS / cubic / quintic, 200 each) | `9.4e-9` / `2.2e-8` / `1.6e-9` / `6.7e-9` (the difference quotient's own error) |
| `NURBSCircle` quarter/half/full: `abs(X(t)) - 1` | `2.2e-16` |
| `NURBSCircularArcDegree2`: `abs(X - C) - r` (relative, 100 arcs) | `1.2e-15` |
| `NURBSSphere` eighth/half/full: `abs(X(u,v)) - 1` | `3.3e-16` |
| natural cubic spline interpolates its samples; derivative and second-derivative jumps at the joints (exact one-sided values from the coefficients) | `4.4e-15` interpolation, `1.6e-14` C1, `4.3e-14` C2 |
| natural quintic spline interpolates `f0` and `f1` | `2.4e-14` and `1.1e-13` |
| `GetTotalLength` versus a 20000-segment chord sum | `1.4e-6` relative — the chord sum is the less accurate of the two |
| `GetTime` inverts `GetLength` | `2.9e-16` relative |
| Frenet frames orthonormal and right-handed (200 curves) | `4.4e-15`; `min det(T,N,B) = 0.9999999999999994` |
| Darboux principal curvatures of the unit sphere are `+-1` | `1.3e-8` near the degenerate pole rows of the NURBS sphere, `2.2e-16` away from them — conditioning of the parameterization, not an error |
| offset vertices at the offset distance from the incident edges | `8.9e-16` |
| `SampleCircularArc` samples on the circle, for arcs of at most pi also angularly ordered | `5.4e-16`; 100/100 monotone |

No convention violation of the kind that produced issue #507 was found: the
generated arcs are counterclockwise with `DotPerp(P0, P2) > 0`, the polylines
use GTE's `Perp(x,y) = (y,-x)` normals, and the drawn weights and times satisfy
the documented preconditions.

## Not covered

* **`CLODPolyline` with an open polyline of three or more vertices.** For every
  such input upstream's `ComputeEdges` writes `edges[1] = collapses[0] + 1`,
  which the `MinHeap` tie-break makes equal to `numVertices`, and
  `ReorderVertices` then evaluates `permute[numVertices]` — a read one past the
  end of a `std::vector<int32_t>`. That is undefined behaviour, so there is no
  upstream answer to compare against and no honest `deviation` case can be
  built: the recorded value would be whatever heap bytes follow the vector in
  that particular build. Issue
  [#182](https://github.com/gradientspaceai/gtengine-js/issues/182); the port's
  fix is pinned by `test/CLODPolyline.test.ts`. The generator therefore uses
  closed polylines (3 to 8 vertices) plus the open two-vertex polyline, which
  upstream returns early for.
* `PolylineOffset::GetDirections`/`GetNormals` are port-only accessors
  (upstream keeps the arrays private); they are exercised by the
  independent-reference check above, not by a case.
* `ParametricCurve::Evaluate(Real, uint32_t, Real*)`, the `reinterpret_cast`
  overload, has no TypeScript analogue and is not ported.
* `MassSpringCurve::ExternalAcceleration` is the zero default; a nonzero
  override would be test code on both sides rather than upstream code, so only
  the default (which `Acceleration` calls on every Runge-Kutta stage) is
  covered.

## Upstream bug suspects

None new. Everything observed in this group is already recorded in
`docs/UPSTREAM-FINDINGS.md` (issues #112, #113, #182, #183, #415).

One observation worth recording for future groups rather than as an upstream
defect: `std::pow(x, 1.5)` in the MSVC runtime and `Math.pow(x, 1.5)` in V8
agreed **bit for bit on all 6000 deep-run arguments** of
`FrenetFrame2.getCurvature`, `FrenetFrame3.getCurvatureAndTorsion` and
`ImplicitCurve2.getCurvature`. The cases keep the default tolerance anyway,
since neither library guarantees it.
