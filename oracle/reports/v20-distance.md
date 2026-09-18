# Oracle group 20 — `v20-distance`

Differential testing of the 19 headers of verify group 20 against the upstream
GTE build (MSVC 194435215, x64, `/O2 /fp:precise`, upstream commit
`d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`). See [ORACLE.md](../../ORACLE.md).

- `oracle/cpp/cases/v20-distance.cpp` — 28 cases
- `oracle/golden/v20-distance.txt` — 560 committed records (20 per case)
- `test/oracle/v20-distance.oracle.test.ts` — the TypeScript replay

**Result.** 22 ordinary cases, 440 records, 6 further cases demonstrate
deliberate fixes of upstream defects. 21 of the 22 ordinary cases are declared
`{ exact: true }` and are **100% bit-identical** to the C++ build, at 20
records per case and at 2000 records per case in the deep run. The one
exception is `DistLine3Circle3.compute`, which is compared with the default
scaled tolerance `1e-12`; its maximum scaled error over 2000 deep-run records
is `2.3e-14` and 69% of its outputs are still bit-identical. The reason is
given under "Comparison policy" below. No case is skipped.

Deep run: `npm run oracle:deep -- 2000 v20` passes (56 000 records).

Three port defects were found and fixed, all of them a reassociated
accumulation of a change of basis; before the fix they made 11 of 20, 20 of 20
and 8 of 20 records of three different cases disagree with the C++ build.

Every generator mixes five modes selected from `io.index()`: uniform random
(modes 0, 1), a small integer lattice (modes 2, 3 — exact arithmetic, and
therefore exactly parallel directions, exactly zero direction components in the
box frame, zero box extents, arc endpoints exactly on the circle and integer
radii) and a construction aimed at a specific branch (mode 4 — the zero
direction, a point exactly at a circle center or on a circle axis or on a
cylinder axis, a line aimed at a point of a rectangle or triangle, a segment
strictly inside a circle, a ray aimed away from the primitive). Every mode
records the same number of doubles, so the replay reads the inputs without
knowing which mode produced them.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `DistLine2Arc2.h` | `compute` | exact | pass |
| `DistLine2OrientedBox2.h` | `compute` | exact | pass |
| `DistLine3AlignedBox3.h` | `compute` | exact | pass |
| `DistLine3Circle3.h` | `compute.axis`, `compute`, `compute.axis.deviation`, `compute.tauHat`, `compute.nearPerpendicular` | exact / tol 1e-12 / deviation | pass |
| `DistLine3OrientedBox3.h` | `compute`, `compute.deviation` | exact / deviation | pass |
| `DistLine3Rectangle3.h` | `compute` | exact | pass |
| `DistLine3Triangle3.h` | `compute` | exact | pass |
| `DistPoint2Arc2.h` | `compute` | exact | pass |
| `DistPoint3Circle3.h` | `compute` | exact | pass |
| `DistPoint3Cylinder3.h` | `compute`, `compute.deviation` | exact / deviation | pass |
| `DistPoint3Frustum3.h` | `compute`, `compute.deviation` | exact / deviation | pass |
| `DistPointAlignedBox.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistPointOrientedBox.h` | `compute.2d`, `compute.3d` | exact | pass |
| `DistRay2Circle2.h` | `compute` | exact | pass |
| `DistRay2Triangle2.h` | `compute` | exact | pass |
| `DistRay3CanonicalBox3.h` | `compute` | exact | pass |
| `DistSegment2Circle2.h` | `compute` | exact | pass |
| `DistSegment2Triangle2.h` | `compute` | exact | pass |
| `DistSegment3CanonicalBox3.h` | `compute` | exact | pass |

Entry points covered: every `operator()` of every header of the group. None of
these headers exposes a second public method, and the port adds no exported
free function for any of them, so `operator()` is the complete public surface.
The private helpers `DistPoint3Cylinder3::DoQueryInfiniteCylinder` /
`DoQueryFiniteCylinder`, `DistPoint3Circle3::ComputeOrthogonalBasis` and
`DistLine3Circle3::Execute` / `PDFSection41x` / `PDFSection42x` / `Finalize` /
`Bisect` are reached through `operator()`; the cases below name which case
reaches which. (`DistLine3Circle3::Execute` is exported by the port as
`distLine3Circle3Execute` for the ray-circle and segment-circle queries, which
belong to other groups; it differs from `operator()` only in also returning the
critical points, and `operator()` is a one-line call to it.)

Branch histogram of the deep run (2000 records per case), the counts that show
the non-generic branches are reached:

- `distance == 0` (the linear component meets the solid): 557 of 2000 for
  `DistLine2OrientedBox2`, 126 for `DistLine3AlignedBox3`, 207 for
  `DistLine3OrientedBox3`, 498 for `DistLine3Rectangle3`, 549 for
  `DistLine3Triangle3`, 190 for `DistRay3CanonicalBox3`, 487 for
  `DistSegment3CanonicalBox3`, 496 for `DistRay2Triangle2`, 676 for
  `DistSegment2Triangle2`, 53 for `DistPoint3Frustum3` (the point is inside the
  frustum), 225 for `DistPoint3Cylinder3` (inside the solid cylinder).
- The endpoint branches of the ray and segment queries: `parameter == 0` in
  1218 of 2000 records of `DistRay2Triangle2`, 1214 of `DistRay3CanonicalBox3`,
  856 of `DistSegment2Triangle2` and 552 of `DistSegment3CanonicalBox3`;
  `parameter == 1` in 67, 69, 368 and 249 records respectively. The clamps are
  exercised on both sides.
- `numClosestPairs`: `DistLine2Arc2` 1537/463 for one/two pairs,
  `DistRay2Circle2` 1819/181, `DistSegment2Circle2` 1394/144 plus **462 records
  with zero pairs** (the preserved upstream quirk that resets the whole result
  for a segment strictly inside the circle), `DistLine3Circle3.compute.axis`
  1801/199.
- The `equidistant` branches: 422 of 2000 for `DistPoint2Arc2` (the point at
  the arc center), 416 for `DistPoint3Circle3` (the point on the circle axis,
  where the scaled projection is exactly the zero vector), 440 for
  `DistLine3Circle3.compute.axis` (PDFSection411).
- `DistLine2OrientedBox2` and `DistLine3AlignedBox3` dispatch on the number of
  nonzero direction components in the box frame. Mode 4 passes the exactly zero
  direction (`DoQuery0D`) and the lattice modes use the coordinate axes as box
  axes, so lattice directions keep their exactly zero components (`DoQuery1D`,
  `DoQuery2D`).
- `DistPoint3Frustum3`'s test point is drawn in frustum coordinates over a
  range that covers all ten Voronoi regions, in the octant-folded form the
  algorithm uses.

## Comparison policy

Every case but one is `{ exact: true }`: the paths in these headers are
`+ - * /`, `sqrt`, `fabs`, `min`, `max` and comparisons only.

`DistLine3Circle3.compute` (PDFSection422) is the exception and uses the
default scaled tolerance `1e-12`. Two properties of that path are intrinsic to
it:

1. it is iterative — `RootsBisection1` halves the bracket up to 4096 times; and
2. the port recovers the critical line parameter as
   `t = G(tau) - Dot(M,D)/Dot(M,M)` instead of upstream's `t = tau + s`.

The second is the documented fix for issue
[#421](https://github.com/gradientspaceai/gtengine-js/issues/421) item 3. It is
a **conditioning** fix, not a branch fix: `s` grows like `1/|NxM|^2` while every
bisection bracket has width `r*|NxM|/Dot(M,M)`, so `tau` approaches `-s` to far
below `ulp(s)` and upstream's sum loses its significant digits continuously as
the line approaches the axis of the circle. There is no exact predicate that
separates the inputs upstream gets right from the ones it gets wrong, so the
substitution cannot be guarded the way the `DistLine2Triangle2` fix of group 19
was narrowed — any guard would be an accuracy threshold, and a guard loose
enough to keep the documented worked example (`7.8000463921561485` where the
truth is `7.8`, a relative error of 6e-6) would fire on essentially every
record anyway.

The residual was root-caused rather than assumed: substituting upstream's
`t = tau + s` back into `src/DistLine3Circle3.ts` makes **every record of this
case bit-identical** to the C++ build, which proves the reassociation is the
only source of the difference and that nothing else on the path deviates. The
generator keeps only configurations on which upstream's own reported distance
agrees with an independent reference (the minimum over the circle of the
point-to-line distance, resolved by a dense scan plus golden-section
refinement) to `1e-9` relative, so the case measures the two expressions where
upstream is sound; the regime in which upstream is actually wrong has its own
`compute.nearPerpendicular` deviation case.

## Port defects fixed

All three are the same defect: upstream accumulates a change of basis left to
right and the port had grouped the basis terms first. The sums are
algebraically equal and differ only in the last bits, but the difference is
real and the oracle measures it. Each fix is a pure reassociation back to
upstream's order, and each is pinned by a new regression test that runs the
port beside both groupings over 500 constructed configurations, asserts the
port equals the left-to-right grouping bit for bit, and asserts that the two
groupings really do differ on the sample so the test cannot pass vacuously.

### `src/DistLine2OrientedBox2.ts` — closest points mapped back out of order

Upstream writes

```cpp
result.closest[i] = box.center + temp[i][0] * box.axis[0]
    + temp[i][1] * box.axis[1];
```

which C++ evaluates as `((center + t0*a0) + t1*a1)`. The port grouped the two
axis terms first. **11 of 20** committed records of
`DistLine2OrientedBox2.compute` disagreed, in `distance`, `sqrDistance` and
both closest points, by up to `5.7e-16` scaled. Regression test:
`test/DistLine2OrientedBox2.test.ts`, `DistLine2OrientedBox2 accumulation
order`.

### `src/DistPoint3Cylinder3.ts` — closest point mapped back out of order

Upstream writes `cylinder.axis.origin + c[0]*basis[1] + c[1]*basis[2] +
c[2]*basis[0]`, that is `(((origin + c0*b1) + c1*b2) + c2*b0)`; the port
grouped the three basis terms first. **20 of 20** committed records of
`DistPoint3Cylinder3.compute` disagreed in `closest[1]`. Regression test:
`test/DistPoint3Cylinder3.test.ts`, `DistPoint3Cylinder3 accumulation order`.

### `src/DistPoint3Frustum3.ts` — closest point mapped back out of order

Upstream writes `frustum.origin + closest[0]*rVector + closest[1]*uVector +
closest[2]*dVector`, that is `(((origin + c0*r) + c1*u) + c2*d)`; the port
grouped the three axis terms first. **8 of 20** committed records of
`DistPoint3Frustum3.compute` disagreed in `closest[1]`. Regression test:
`test/DistPoint3Frustum3.test.ts`, `DistPoint3Frustum3 accumulation order`.

After these three fixes, all 440 committed records and all 44 000 deep-run
records of the 21 `exact` cases agree bit for bit, and
`DistLine3Circle3.compute` agrees to `2.3e-14` scaled.

## Deliberate deviations demonstrated

Each of these cases generates only inputs on which upstream is defective; the
harness requires the port to disagree with the C++ build on at least one
record, so the fix is demonstrated against the real build rather than against a
TypeScript re-derivation. All six deviate on 20 of 20 committed records.

| case | records disagreeing | upstream defect |
| --- | ---: | --- |
| `DistLine3OrientedBox3.compute.deviation` | 20 of 20 | `operator()` assigns `result.closest[0] = line.origin + parameter * line.direction` (world space) *before* the loop that reads `result.closest[i]` as box-frame coordinates and maps them to the world, so upstream returns `box.center + sum_j (world line point)[j] * box.axis[j]` as the closest line point. The port drops the premature assignment. `distance`, `sqrDistance`, `parameter` and `closest[1]` come from the canonical-box query and are unaffected, so the main `compute` case emits those four and is bit-identical; this case emits `closest[0]` alone. Max scaled error 1.8. |
| `DistPoint3Cylinder3.compute.deviation` | 20 of 20 (C++ throws on all 20) | `operator()` detects an infinite cylinder with `cylinder.height == std::numeric_limits<T>::max()`, but `Cylinder3::MakeInfiniteCylinder` sets and documents `height = -1`. The sentinel therefore falls into the finite branch and trips `LogAssert(height > 0)`: upstream throws on every infinite cylinder built the documented way. The port uses `Cylinder3::isInfinite()` and answers the infinite-cylinder query. The generator passes `height = -1`. |
| `DistPoint3Frustum3.compute.deviation` | 20 of 20 | Two of the ten far-edge assignments do not clamp the free coordinate of the edge: in the `test[2] <= dmin` group and in the `dmin < test[2] < dmax` group, the `rEdgeDot >= 0` then `rdDot >= maxRDDot` branch reaches the LF-edge with `test[1] > umax`, and symmetrically the `uEdgeDot` branch reaches the UF-edge with `test[0] > rmax`. The reported closest point is off the end of the edge, outside the frustum, and the distance is too small. The generator replicates the branch analysis in C++ and keeps only the records that reach one of those two assignments with the coordinate actually past the end. The other eight assignments are bounded by their own branch conditions, where the port's clamp is a no-op — which is why the main case is bit-identical. Max scaled error 1.2. |
| `DistLine3Circle3.compute.axis.deviation` | 20 of 20 | `Finalize` normalizes the in-plane component of a critical line point without checking that it is nonzero; `Vector::Normalize` leaves a zero vector at zero, so upstream reports the circle **center** as the closest circle point and the distance to the center instead of to the circle (issue [#421](https://github.com/gradientspaceai/gtengine-js/issues/421)). PDFSection421 reaches it exactly on the integer lattice whenever a critical parameter puts the line point on the axis of the circle. Since every circle point is then equidistant from that line point, the port's corrected distance `sqrt(h^2 + r^2)` always exceeds upstream's `|h|`, so the two agree unless upstream *selected* the defective critical point; the generator uses that exact symptom (upstream reports `circularClosest[j] == circle.center`) as its selector, and the main `compute.axis` case rejects it. Max scaled error 2.0. |
| `DistLine3Circle3.compute.tauHat` | 20 of 20 | `PDFSection422` solves `G'(tauHat) = 1` as `tauHat = sqrt(|(a1*a3)^{2/3} - a3|)`, omitting the division by `a2`, so `tauHat` is a factor `sqrt(a2)` too large (issue [#247](https://github.com/gradientspaceai/gtengine-js/issues/247)). The generator keeps only records for which `tauHat` can influence the answer — `a1 > sqrt(a3)` and `abs(a0)` no larger than either intercept, since above both intercepts the two implementations take the same single-critical-point branch with the same bracket — *and* for which upstream's reported distance is wrong by more than `1e-9` relative to the independent reference. The direction vector is long (`a2 = |NxM|^2 > 1`), which is the regime the defect needs and the way the segment-circle query calls this file. Max scaled error 2.0. |
| `DistLine3Circle3.compute.nearPerpendicular` | 20 of 20 | The back-substitution `t = tau + s` cancels away the significant digits of the critical parameter for a line nearly perpendicular to the plane of the circle (issue [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) item 3, see "Comparison policy"). The generator perturbs the direction off the circle normal by `eps` drawn log-uniformly in `[10^-8.5, 10^-6]`, rejects the configurations for which upstream's bisection bracket collapses entirely (upstream throws there, a different symptom of the same finding), and keeps only those for which upstream's reported distance is wrong by more than `1e-9` relative to the independent reference. Max scaled error 2.4e-8. |

## Not covered

Nothing. All 19 headers of the group are implemented by the port and all have
cases. No header in the group is pure data.

Two narrower notes, not gaps in entry-point coverage:

- `DistPoint3Cylinder3`'s infinite-cylinder branch is reachable in upstream
  only by setting `height` to exactly `std::numeric_limits<double>::max()`,
  which no `Cylinder3` method produces. It is covered through the deviation
  case, which passes the documented `height = -1` sentinel; the port answers
  the infinite-cylinder query there and upstream throws. For
  `height = DBL_MAX` the two agree: upstream takes the infinite branch and the
  port takes the finite branch with `halfHeight = DBL_MAX/2`, which never
  clamps, so the results are identical.
- `DistPointAlignedBox` and `DistPointOrientedBox` are `N`-dimensional
  templates and are covered in 2D and 3D. `DistLine3Circle3`'s two-closest-pair
  and `equidistant` results are produced by the closed-form branches and are
  covered by `compute.axis` (199 and 440 records of 2000); the generic
  `compute` case reaches them with probability zero because they require two
  bisected critical distances to be exactly equal.

## Upstream bug suspects

None new. The five defects exercised above are already recorded in
`docs/UPSTREAM-FINDINGS.md`, as is the preserved `DistSegment2Circle2` quirk
that reports `numClosestPairs = 0` with `distance = 0` for a segment strictly
inside the circle; the deep run confirms that quirk against the real build on
462 of 2000 records, where the port reproduces upstream bit for bit.

One observation worth recording, though it is a consequence of an already
documented defect rather than a new one: `DistPoint3Cylinder3.h` and
`Cylinder3.h` disagree about the infinite-cylinder sentinel in a way that makes
the upstream infinite-cylinder path unreachable from the upstream API. Every
call built with `MakeInfiniteCylinder` throws.
