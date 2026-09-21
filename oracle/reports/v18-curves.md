# Verify group 18 (curves) against the MSVC build of upstream GTE

Family `v18-curves`, 24 cases, 20 golden records each, deep run
`npm run oracle:deep -- 2000 v18-curves` (48000 records, 4000 of them throw
records): **every case passes**, and the deviation case deviates on all 2000
records. Wall time of the deep run, generation and replay together: about 12
seconds.

Twenty of the twenty-four cases are declared `{ exact: true }`, including the
whole geodesic algorithm; the four that call `sin`/`cos` carry a tolerance,
and their measured maxima are reported below.

## How the geodesic algorithm is compared bit for bit

`RiemannianGeodesic` is abstract, and GTE's only concrete subclass is
`EllipsoidGeodesic`, whose metric calls `sin` and `cos`. That matters because
`Refine`'s steepest-descent search selects the *argument of the minimum* over
`2*searchSamples+1` sampled lengths, and a search that compares libm-derived
quantities is not comparable at any tolerance: one ulp elsewhere sends the
midpoint to a different sample.

The case file therefore defines its own concrete subclass,
`PolySurfaceGeodesic`, over the polynomial graph `z = u^2*v + u*v^2`. Its
metric is the Gram matrix of the tangent vectors, its Christoffel symbols of
the first kind are `Dot(P_ij, P_k)` (upstream's convention), and both are
polynomials. Everything downstream - `ComputeSegmentLength`,
`ComputeTotalLength`, `ComputeIntegrand`, `ComputeSegmentCurvature`,
`ComputeTotalCurvature`, `ComputeMetricInverse` (Gaussian elimination with
full pivoting), `ComputeChristoffel2`, `ComputeMetricDerivative`, `Refine`,
`Subdivide`, `ComputeGeodesic` and `refineCallback` - is then `+ - * / sqrt`,
so every comparison that decides a branch is exact on both sides and the
cases are bit-for-bit. The same subclass is written in the TypeScript replay.

`EllipsoidGeodesic` then only has to be compared where it differs from the
generic algorithm, which is the metric and the Christoffel symbols.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `NaturalSplineCurve.h` | `free.evaluate`, `clamped.evaluate`, `twoPoints.evaluate` | exact | pass, 2000 records each |
| `NaturalSplineCurve.h` | `closed.evaluate` (numPoints >= 4) | exact | pass |
| `NaturalSplineCurve.h` | `closed.wrapRow.deviation` (numPoints 2, 3) | deviation (#295) | 2000 of 2000 records deviate |
| `NaturalSplineCurve.h` | `construct.tooFewPoints` | throw parity | pass, 2000 throw records |
| `RiemannianGeodesic.h` | `poly.tensors` (metric inverse, Christoffel 2, metric derivative) | exact | pass |
| `RiemannianGeodesic.h` | `poly.computeIntegrand`, `poly.computeSegmentLength`, `poly.computeTotalLength`, `poly.computeSegmentCurvature`, `poly.computeTotalCurvature` | exact | pass |
| `RiemannianGeodesic.h` | `poly.refine`, `poly.subdivide`, `poly.computeGeodesic`, `poly.refineCallback` | exact | pass |
| `RiemannianGeodesic.h` | `construct.invalidDimension` | throw parity | pass, 2000 throw records |
| `EllipsoidGeodesic.h` | `computePosition`, `computeSegmentLength`, `computeTotalLength`, `computeSegmentCurvature`, `computeTotalCurvature` | tolerance 1e-12 (sin, cos) | pass |
| `EllipsoidGeodesic.h` | `computeGeodesic.noSearch` | exact | pass |
| `EllipsoidGeodesic.h` | `refine.separated` | tolerance 1e-10 (sin, cos, amplified) | pass |

`Evaluate` is driven at five parameters per record - one below the first
knot, one above the last knot, one exactly at a knot, and a pair straddling an
interior knot - with the jet order drawn over 0..4, so `GetKeyInfo`'s three
branches and every early exit of `Evaluate`'s order ladder are reached. The
knots are drawn uniformly spaced, nonuniformly spaced or on an integer
lattice; the points are drawn uniformly or on a lattice.

`RiemannianGeodesic.poly.refineCallback` records `GetSubdivisionStep`,
`GetRefinementStep` and `GetCurrentQuantity` from inside the callback during
`ComputeGeodesic`; the trace lengths in the deep run range over 0, 3, 6, 12
and 24 entries, and the absence of entries from the subdivision pass confirms
that `Subdivide` installs and restores the no-op callback.

## Port defects fixed

None. Every arithmetic output of all 20 ordinary cases was bit-identical to
the MSVC build on the first run, over 48000 deep-run records.

Two fidelity questions from ORACLE.md's triage list were checked explicitly
against the real build and came out clean:

* The port's shared `dot` switches on `GVector`'s `dotAccumulatesFromZero`
  flag, so `Vector.h`'s first-product seed and `GVector.h`'s literal-zero seed
  are both reproduced. `RiemannianGeodesic.ts` builds its difference vectors
  with `sub()` from `Vector.ts`, which returns the base `Vector` type, so those
  `dot` calls use the first-product seed where upstream's `GVector` argument
  would use the zero seed. The two differ only when the first product is `-0`
  *and* the remaining sum is `±0`, which makes the quadratic form zero and
  fires upstream's own `LogAssert(qForm > 0)` on both sides; 48000 records
  produced no disagreement.
* `Inverse(mMetric, &mMetricInverseExists)` in `ComputeMetricInverse` takes a
  `GMatrix`, so it resolves to `GMatrix.h`'s Gaussian-elimination `Inverse`,
  not to the closed-form `Matrix2x2.h` overload that `LinearSystem.h` brings
  into scope. The port calls `inverse()` from `Matrix.ts`, which is Gaussian
  elimination, and `poly.tensors` confirms it bit for bit on 2000 records.

## Deliberate deviations demonstrated

`NaturalSplineCurve::CreateClosed` builds the wrap-around row of its linear
system with three plain assignments, to columns `numSegments-1`, `0` and `1`.
Those columns are distinct only when `numSegments >= 3`. For
`numSegments == 2` the columns `numSegments-1` and `1` coincide, so the third
assignment discards the `dt[numSegments-1]` term; for `numSegments == 1` the
columns `numSegments-1` and `0` coincide and the second assignment discards
the `dt[numSegments-1]` term. The wrap equation is the ordinary C2 joint
condition `dt[i-1]*c[i-1] + 2*(dt[i-1]+dt[i])*c[i] + dt[i]*c[i+1] = rhs` at
`i = 0` with `i-1` taken modulo the number of segments, so when two of those
indices name the same unknown their coefficients must *add*. The port
accumulates with `+=` into a zero-filled matrix, which is identical to
upstream whenever the columns are distinct (issue
[#295](https://github.com/gradientspaceai/gtengine-js/issues/295)).

`closed.evaluate` therefore draws `numPoints >= 4` and agrees bit for bit on
2000 records, which is the evidence that the fix is confined to the defective
inputs; `closed.wrapRow.deviation` draws `numPoints` 2 or 3 and disagrees on
all 2000 of its records, with relative differences up to 2.0.

Preserved upstream defects of this group, compared bit for bit rather than
deviated from:

* `RiemannianGeodesic::ComputeMetricDerivative` computes `2*Gamma_{d,i0i1}`,
  not `dg_{i0i1}/dx_d`. It is emitted by `poly.tensors` and matches exactly.
  (It is never called by the algorithm itself; the case reaches it through a
  public probe added to the test subclass on both sides.)
* The derived parameters `mIntegralStep`, `mSearchStep` and
  `mDerivativeFactor` are computed in the constructor only, although
  `integralSamples`, `searchSamples`, `derivativeStep` and `searchRadius` are
  public and documented as tweakable. The cases change those fields after
  construction on purpose - `integralSamples` is drawn over 2..16 while
  `mIntegralStep` stays `1/15`, and `searchSamples` is changed while
  `mSearchStep` stays `1/32` - and the port reproduces the stale values
  exactly. The port's `updateDerivedParameters()` escape hatch is an addition
  with no upstream counterpart, so no case calls it; calling it would be a
  deviation, and it is covered by the port's own unit tests instead.
* `LogAssert(subdivisions < 32)` admits 31, where `1 << subdivisions` is
  signed overflow. The generators stay far below that (`subdivisions <= 3`)
  and the assert's boundary is not exercised, because upstream's behaviour
  there is undefined.
* `CreateFree`'s `storageSize` over-allocates by `N-1` reals. It has no
  observable effect and the port does not reproduce the packing at all (it
  uses four separate arrays), so there is nothing to compare.

## Tolerances and their proven causes

| case | tolerance | measured maximum over 2000 records | cause |
| --- | --- | --- | --- |
| `EllipsoidGeodesic.computePosition` | 1e-12 (default) | 3.75e-16 | `sin`, `cos` |
| `EllipsoidGeodesic.computeSegmentLength` | 1e-12 | 3.69e-16 | `sin`, `cos` in `ComputeMetric` |
| `EllipsoidGeodesic.computeTotalLength` | 1e-12 | 3.11e-16 | as above |
| `EllipsoidGeodesic.computeSegmentCurvature` | 1e-12 | 3.33e-16 | as above |
| `EllipsoidGeodesic.computeTotalCurvature` | 1e-12 | 3.33e-16 | as above |
| `EllipsoidGeodesic.refine.separated` | 1e-10 | 3.57e-13 | the centered-difference gradient multiplies a libm-level difference of two nearly equal lengths by `mDerivativeFactor = 0.5/1e-4 = 5000` |

`EllipsoidGeodesic.refine.separated` is the only case whose *control flow*
depends on libm values: `Refine` keeps the first sampled length that strictly
beats the running minimum, and the boolean it returns compares that minimum
with the length at the incoming midpoint. The C++ generator replicates the
sampled lengths (`EllipsoidSearchSeparation`) and rejects, with a cap of 24
attempts and a full redraw of every quantity per attempt, any configuration
whose winning length is within a relative `1e-9` of another sample or of the
length at the midpoint. On the accepted inputs both builds select the same
sample and the same boolean; what remains is the rounding of the gradient,
three orders of magnitude below the case tolerance and ten orders below the
distance between two adjacent search samples, so a flipped argument of the
minimum would still fail the case. The acceptance threshold (1e-9) and the
comparison tolerance (1e-10) are chosen together: the band between them is
empty because a same-index result differs by at most 3.6e-13.

`EllipsoidGeodesic.computeGeodesic.noSearch` sets `searchSamples = 0`, which
leaves the search with the single sample `tRay = 0`. Its candidate point is
`mid - 0*gradient`, bit-identical to the midpoint `Subdivide` computed, so
whichever way the libm-decided comparison goes the point that leaves `Refine`
is the arithmetic midpoint and no libm value reaches the path. The case is
therefore exact and compares `ComputeGeodesic`'s own bookkeeping - the vertex
count, the interleaving copy, the subdivision order and the refinement loop -
on the real ellipsoid.

## Independent-reference checks

The deep-run outputs of the main cases were checked once against references
written from scratch (no port, no GTE):

| check | records | median | maximum |
| --- | --- | --- | --- |
| `free/closed/clamped/twoPoints`: `X(t > tmax)` equals the last interpolation point | 2000 each | 2.2e-16 | 7.1e-15 |
| C2 joints: the jump in position, velocity and acceleration across an interior knot, divided by `2.5*delta*(bound on the next derivative)` | 2382-2460 samples each | 0.80 | 0.80 (never above 1) |
| a sphere's `computeSegmentLength` against the closed-form metric `diag(r^2 sin^2 v, r^2)`, derived from the geometry rather than from GTE's `Dot(dP/du, dP/du)` | 667 | 0 | 4.4e-16 |
| `poly.tensors`: the emitted metric equals the Gram matrix of the tangent vectors of `z = u^2 v + u v^2` | 2000 | 0 | 0 |
| `poly.tensors`: `g * g^-1 - I` | 2000 | 0 | 1.4e-14 |
| `poly.computeGeodesic`: the length of the refined polyline divided by the length of the straight parameter-space segment | 2000 | 0.9979 | 1.00068 |

The geodesic ratio is below 1 on 1863 of 2000 records, down to 0.647, which
is the algorithm doing its job. The 137 records slightly above 1 exceed it by
at most 0.07%: splitting a segment and applying the 16-sample Trapezoid Rule
to each half is a different quadrature from one rule over the whole segment,
and that discretization difference is `O(h^2)` with either sign. It is not a
failure of the refinement.

## Not covered

* **`ParametricCurve`'s inherited interface** (`GetTMin`, `GetTotalLength`,
  `GetSpeed`, `GetTime(length)`, `SubdivideByTime`, `SubdivideByLength` and
  the Romberg integration behind them). Those belong to `ParametricCurve.h`,
  which is another verify group's header; only the members
  `NaturalSplineCurve.h` itself defines (`Evaluate`, `GetNumPoints`,
  `GetPoints`, the three construction paths) are covered here.
* **`NaturalSplineCurve` with `numPoints = 0`.** `ParametricCurve`'s
  constructor runs before the `LogAssert` and would build a vector of
  `SIZE_MAX` elements, so the case is not generated. `numPoints = 1` is
  generated and gives throw parity (the base constructor allocates a
  one-element time array, copies one time and reads nothing out of bounds).
* **`RiemannianGeodesic` with `subdivisions >= 31`.** `1 << subdivisions` is
  signed overflow, which is undefined behaviour.
* **`updateDerivedParameters()`**, which has no upstream counterpart (see
  above).
* **`EllipsoidGeodesic` near the poles.** The generators keep the `v`
  parameter inside `[0.6, 2.6]`. At `sin(v) = 0` the metric is singular,
  `Inverse` reports non-invertibility and `LogAssert(qForm > 0)` fires, all
  decided by the rounding of `sin`, so throw parity cannot be promised there.
  The degenerate branch of `ComputeMetricInverse` is instead covered exactly
  by the polynomial manifold's own cases, whose metric is never singular, and
  the invertibility flag is emitted by `poly.tensors`.
* **`EllipsoidGeodesic::ComputeGeodesic` with a live search.** The search is
  run and compared through `EllipsoidGeodesic.refine.separated`, where the
  probe can bound a single `Refine`. Bounding the separation of every one of
  the tens of `Refine` calls inside a full `ComputeGeodesic` would mean
  re-implementing the whole algorithm in the case file; `ComputeGeodesic`'s
  own bookkeeping is covered exactly, on the ellipsoid by the `noSearch` case
  and with a live search by `RiemannianGeodesic.poly.computeGeodesic`.

## Upstream bug suspects

No new ones. Everything observed is already recorded in
`docs/UPSTREAM-FINDINGS.md` under issue #295: the wrap-around row of
`CreateClosed` (fixed in the port and demonstrated here), the `LogAssert`
that runs after the base constructor, the over-allocating `storageSize` of
`CreateFree`, `ComputeMetricDerivative` computing `2*Gamma`, the derived
parameters that are never recomputed, and the `subdivisions < 32` assert.
