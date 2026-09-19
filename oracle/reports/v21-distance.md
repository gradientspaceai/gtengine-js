# Oracle group 21 — `v21-distance`

Differential testing of the 21 headers of verify group 21 against the upstream
GTE build (MSVC 194435215, x64, `/O2 /fp:precise`, upstream commit
`d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`). See [ORACLE.md](../../ORACLE.md).

- `oracle/cpp/cases/v21-distance.cpp` — 38 cases
- `oracle/golden/v21-distance.txt` — 760 committed records (20 per case)
- `test/oracle/v21-distance.oracle.test.ts` — the TypeScript replay

**Result.** 34 ordinary cases, 680 records, 5295 floating-point outputs. 32 of
the 34 are declared `{ exact: true }` and are **100.00% bit-identical** to the
C++ build (4975 outputs at 20 records per case, 63 000 records in the deep
run). The two exceptions are `DistRay3Circle3.compute` and
`DistSegment3Circle3.compute`, compared with the default scaled tolerance
`1e-12`; their maximum scaled error over 2000 deep-run records each is
`4.4e-15` and `2.9e-15`, and the reason is root-caused below. 4 further cases
demonstrate deliberate fixes of upstream defects; each deviates on every one of
its records. No case is skipped.

Deep run: `npm run oracle:deep -- 2000 v21` passes (76 000 records).

Two port defects were found, both the same defect: a reassociated left-to-right
accumulation of a change of basis. Before the fix they made 3 of 20 records of
three different cases disagree with the C++ build. One of the two is in
`DistLine2OrientedBox2`, which group 20 found independently and fixed; this
branch is rebased on that fix and carries only the other one.

Every generator mixes five modes selected from `io.index()`: uniform random
(modes 0, 1), a small integer lattice (modes 2, 3 — exact arithmetic, and
therefore exactly zero direction and normal components, zero box extents, equal
hyperellipsoid extents, arc endpoints exactly on the circle with integer
coordinates, integer radii and coplanar tetrahedra) and a construction aimed at
a specific branch (mode 4 — a ray or segment aimed at a point of the primitive
with the primitive behind the origin about half the time, a query point exactly
on a hyperplane, on a parallelogram edge, on a hyperellipsoid or inside a
tetrahedron, the exactly zero plane normal). Every mode records the same number
of doubles, so the replay reads the inputs without knowing which mode produced
them.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `DistPlane3CanonicalBox3.h` | `compute` | exact | pass |
| `DistPoint2Parallelogram2.h` | `compute`, `compute.degenerate`, `getMinimizer` | exact | pass |
| `DistPoint3ConvexPolyhedron3.h` | `compute.tetrahedron`, `compute.box`, `compute.noPlanes`, `setMaxLCPIterations` | exact | pass |
| `DistPointHyperellipsoid.h` | `compute.2d`, `compute.3d`, `computeAxisAligned.2d`, `computeAxisAligned.3d` | exact | pass |
| `DistPointHyperplane.h` | `compute.2d`, `compute.3d`, `compute.4d` | exact | pass |
| `DistRay2AlignedBox2.h` | `compute` | exact | pass |
| `DistRay2Arc2.h` | `compute` | exact | pass |
| `DistRay2OrientedBox2.h` | `compute` | exact | pass |
| `DistRay3AlignedBox3.h` | `compute` | exact | pass |
| `DistRay3Circle3.h` | `compute.axis`, `compute`, `compute.axis.deviation` | exact / tol 1e-12 / deviation | pass |
| `DistRay3OrientedBox3.h` | `compute`, `compute.deviation` | exact / deviation | pass |
| `DistRay3Rectangle3.h` | `compute` | exact | pass |
| `DistRay3Triangle3.h` | `compute` | exact | pass |
| `DistSegment2AlignedBox2.h` | `compute` | exact | pass |
| `DistSegment2Arc2.h` | `compute` | exact | pass |
| `DistSegment2OrientedBox2.h` | `compute` | exact | pass |
| `DistSegment3AlignedBox3.h` | `compute` | exact | pass |
| `DistSegment3Circle3.h` | `compute.axis`, `compute.twoPairs`, `compute`, `compute.axis.deviation` | exact / tol 1e-12 / deviation | pass |
| `DistSegment3OrientedBox3.h` | `compute`, `compute.deviation` | exact / deviation | pass |
| `DistSegment3Rectangle3.h` | `compute` | exact | pass |
| `DistSegment3Triangle3.h` | `compute` | exact | pass |

Entry points covered: every `operator()` of every header of the group, plus

- `DistPoint2Parallelogram2::GetMinimizer`, which is a public member function
  upstream and a public method of the port (`getMinimizer`), covered by its own
  case as well as through `operator()`;
- `DistPoint3ConvexPolyhedron3`'s two other public entry points, the
  `numTriangles` constructor argument (`compute.box` builds the query with
  `numTriangles = 12`, the cached-solver path; the other three cases use the
  `numTriangles = 0` per-query path) and `SetMaxLCPIterations`;
- both public `operator()` overloads of `DistPointHyperellipsoid`, the general
  one and the axis-aligned one that takes the extents directly
  (`computeAxisAligned` in the port).

The port exports no additional free function for any header of this group, and
no header of the group is pure data. The private helpers reached through
`operator()` are `DistPlane3CanonicalBox3::DoQuery3D` / `DoQuery2D` /
`DoQuery1D` / `DoQuery0D`, `DistPointHyperellipsoid::SqrDistance` /
`SqrDistanceSpecial` / `Bisector`, and the `Execute` /
`HasOneCriticalPoint` / `HasTwoCriticalPoints` / `RayOriginClosest` /
`SegmentEndpointClosest` / `SelectClosestPoint` helpers of the ray-circle and
segment-circle queries; the notes below name which case reaches which.

Branch histogram of the deep run (2000 records per case), the counts that show
the non-generic branches are reached:

- The endpoint clamps of the wrappers. `parameter == 0` (the ray origin or the
  first segment endpoint) in 920 records of `DistRay2AlignedBox2`, 940 of
  `DistRay3AlignedBox3`, 898 of `DistRay2OrientedBox2`, 953 of
  `DistRay3OrientedBox3`, 1008 of `DistRay3Rectangle3`, 1011 of
  `DistRay3Triangle3`, 549/490/514/443/561/471 of the six segment cases;
  `parameter == 1` (the second segment endpoint) in 605, 446, 663, 517, 421 and
  322 records of the segment cases and in 84, 58, 85, 62, 72 and 29 records of
  the ray cases (where it is the ordinary interior solution, not a clamp).
- `distance == 0` (the linear component meets the solid or the surface): 486,
  297, 476, 305 for the four aligned-box cases, 587, 681, 399, 451 for the four
  oriented-box cases, 330/292 for the rectangles, 265/289 for the triangles,
  832 for `DistPlane3CanonicalBox3` (the plane cuts the box), 428 for
  `DistPoint2Parallelogram2` (the point is inside), 755/453 for the two arc
  cases.
- `DistPlane3CanonicalBox3` dispatches on the number of strictly positive
  components of the octant-folded normal: 883 records reach `DoQuery3D`, 288
  `DoQuery2D`, 429 `DoQuery1D` and 400 the zero-normal `DoQuery0D`, which
  upstream calls a "low-probability event".
- `DistPoint2Parallelogram2.getMinimizer` lands on the domain boundary of
  `[-1,1]^2` in both coordinates: `K[0]` is exactly `-1` in 622 records and
  exactly `+1` in 617, `K[1]` in 622 and 654; the interior (the "query point is
  inside" arm) accounts for the rest. `compute.degenerate` throws on 2000 of
  2000 records, and the port throws on the same records.
- `DistPointHyperplane` puts the query point exactly on the hyperplane
  (`signedDistance == 0`) in 253, 240 and 237 records of the 2D, 3D and 4D
  cases.
- `DistPointHyperellipsoid` reaches `distance == 0` (the point is exactly on
  the hyperellipsoid, the `sumZSqr == 1` early return of `Bisector`) in 122,
  53, 214 and 148 records of the four cases; mode 4 zeroes one coordinate
  exactly in every record, which is the `y[i] == 0` branch of
  `SqrDistanceSpecial` and, when it is the smallest-extent coordinate, the
  subhyperellipsoid branch.
- `DistPoint3ConvexPolyhedron3`: the LCP solver succeeds on 2000 of 2000
  tetrahedron records (443 with the point inside, distance below `1e-9`) and on
  1962 of 2000 box records (198 inside, 38 LCP failures, which exercise the
  `queryIsSuccessful == false` path with a nonzero iteration count);
  `compute.noPlanes` takes the `numTriangles == 0` early return on all 2000
  records and `setMaxLCPIterations` makes the solver give up on 1993 of 2000.
  `numLCPIterations` is emitted as an integer and compares exactly everywhere,
  so the port's Lemke pivoting follows upstream's step for step.
- `DistRay2Arc2` / `DistSegment2Arc2`: `numClosestPairs` is 2 in 147 and 93
  records (both circle closest points on the arc), 1 in the rest; the
  endpoint/arc sort arm (no circle closest point on the arc) is reached
  whenever the ray or segment misses the arc.
- `DistRay3Circle3.compute.axis` / `DistSegment3Circle3.compute.axis` reach the
  three closed-form branches: `equidistant` is true in 440 and 447 records
  (PDFSection411, the linear component along the axis of the circle), and the
  lattice mode puts the origin on the axis (PDFSection421).
  `DistSegment3Circle3.compute.twoPairs` returns `numClosestPairs == 2` on
  2000 of 2000 records.

## Comparison policy

32 of the 34 ordinary cases are `{ exact: true }`: their paths are `+ - * /`,
`sqrt`, `fabs`, `min`, `max` and comparisons only. The arc generators call
`cos`, `sin` and `atan2` to place the arc endpoints and the aim point, but
those values are recorded as inputs with `io.givenVec`, so no libm result
enters the computation under test.

`DistRay3Circle3.compute` and `DistSegment3Circle3.compute` use the default
scaled tolerance `1e-12`. Both reach `DistLine3Circle3`'s PDFSection422, which
is iterative (`RootsBisection1`) and in which the port recovers the critical
line parameter as `t = G(tau) - Dot(M,D)/Dot(M,M)` instead of upstream's
`t = tau + s`. That is the documented conditioning fix of issue
[#421](https://github.com/gradientspaceai/gtengine-js/issues/421) item 3; group
20 measured and reported the same residual for `DistLine3Circle3.compute`
itself.

The residual was root-caused rather than assumed. Substituting upstream's
`t = tau + s` back into `src/DistLine3Circle3.ts` makes **every output of both
cases bit-identical** to the C++ build over the full 2000-record deep run
(`bitident 100.00`, `maxRelErr 0`), which proves the reassociation is the only
source of the difference and that nothing else on the ray and segment paths
deviates. The generators keep only configurations on which upstream's own
reported distance agrees, to `1e-9` relative, with an independent C++-side
reference (the minimum over the circle of the point-to-clamped-line distance,
resolved by a 2048-sample scan plus a golden-section refinement) and on which
the `tauHat` defect of issue
[#247](https://github.com/gradientspaceai/gtengine-js/issues/247) cannot
influence the answer, so the cases measure the two expressions where upstream
is sound. The regimes in which upstream is actually wrong belong to group 20's
`DistLine3Circle3.compute.tauHat` and `compute.nearPerpendicular` cases, which
exercise the same code through the line query.

The closed-form branches of the two headers are measured bit for bit by
`compute.axis` and `compute.twoPairs`, which reach them through the ray and
segment wrappers.

## Port defects fixed

Both are the same defect: upstream accumulates a change of basis left to right
and the port had grouped the basis terms first. The sums are algebraically
equal and differ only in the last bits, but the difference is real and the
oracle measures it.

### `src/DistPoint2Parallelogram2.ts` — the closest point mapped back out of order

Upstream writes

```cpp
result.closest[1] = pgm.center + K[0] * pgm.axis[0] + K[1] * pgm.axis[1];
```

which C++ evaluates as `((C + K0*A0) + K1*A1)`. The port grouped the two axis
terms first. **3 of 20** committed records of
`DistPoint2Parallelogram2.compute` disagreed, in `distance`, `sqrDistance` and
`closest[1]`, by up to `8.6e-16` scaled.

Regression test: `test/DistPoint2Parallelogram2.test.ts`, "accumulates the
closest point left to right, as upstream does". It runs the port beside both
groupings over 500 constructed configurations, asserts the port equals the
left-to-right grouping bit for bit, and asserts that the two groupings really
do differ on the sample, so the test cannot pass vacuously. Reverting the
`src/` change makes it fail.

### `src/DistLine2OrientedBox2.ts` — the same defect, fixed by group 20

`DistRay2OrientedBox2` and `DistSegment2OrientedBox2` call
`DistLine2OrientedBox2`, whose closest-point reconstruction

```cpp
result.closest[i] = box.center + temp[i][0] * box.axis[0]
    + temp[i][1] * box.axis[1];
```

had the same reassociation. **3 of 20** committed records of each of
`DistRay2OrientedBox2.compute` and `DistSegment2OrientedBox2.compute`
disagreed, in `distance`, `sqrDistance` and both closest points.

Group 20 found this independently (`DistLine2OrientedBox2.compute`, 11 of 20
records) and fixed it with the identical one-hunk `src/` change plus a
regression test in `test/DistLine2OrientedBox2.test.ts`. Group 21 was developed
from the commit before that landed and reproduced the same fix byte for byte;
after rebasing onto `main` with group 20 merged, this branch no longer carries
the hunk and group 20's regression test is the one that pins the order. The
finding is recorded here because two different headers of this group reach the
defect and because it is the second independent confirmation of it.

After these two fixes, all 640 committed records and all 63 000 deep-run
records of the 32 `exact` cases agree bit for bit.

## Deliberate deviations demonstrated

Each of these cases generates only inputs on which upstream is defective; the
harness requires the port to disagree with the C++ build on at least one
record, so the fix is demonstrated against the real build rather than against a
TypeScript re-derivation. All four deviate on 20 of 20 committed records and on
2000 of 2000 deep-run records.

| case | records disagreeing | upstream defect |
| --- | ---: | --- |
| `DistRay3OrientedBox3.compute.deviation` | 20 of 20 | `DistLine3OrientedBox3::operator()` assigns `result.closest[0] = line.origin + parameter * line.direction` (world space) *before* the loop that reads `result.closest[i]` as box-frame coordinates and maps them to the world, so upstream returns `box.center + sum_j (world line point)[j] * box.axis[j]` as the closest linear point (issue [#421](https://github.com/gradientspaceai/gtengine-js/issues/421)). The ray query inherits it whenever it returns the line-box output. The port drops the premature assignment. Max scaled error 2.0. |
| `DistSegment3OrientedBox3.compute.deviation` | 20 of 20 | The same defect, inherited by the segment query on the branch that returns the line-box output. Max scaled error 2.0. |
| `DistRay3Circle3.compute.axis.deviation` | 20 of 20 | `DistLine3Circle3::Finalize` normalizes the in-plane component of a critical line point without checking that it is nonzero; `Vector::Normalize` leaves a zero vector at zero, so upstream reports the circle **center** as the closest circle point and the distance to the center instead of to the circle (issue [#421](https://github.com/gradientspaceai/gtengine-js/issues/421)). PDFSection421 reaches it on the integer lattice whenever a critical parameter puts the line point on the axis of the circle, and the ray query keeps that critical point. Max scaled error 2.0. |
| `DistSegment3Circle3.compute.axis.deviation` | 20 of 20 | The same defect, reached through the segment query. Max scaled error 2.0. |

The selectors are the observable symptoms, not replicas of the branch analysis.
For the two oriented-box cases the generator keeps only the records for which
upstream's own `parameter` is strictly inside the domain of the linear
component, which is exactly when the query returns the line-box output and the
doubly-transformed `closest[0]` appears; when the parameter is clamped the
query writes the ray origin or a segment endpoint into `closest[0]` itself and
the port agrees, so keeping those records would make the deviation look
narrower than it is. For the two circle cases the selector is
`circularClosest[j] == circle.center`, which is upstream's symptom, and the
ordinary `compute.axis` cases reject exactly those records.

The two ordinary oriented-box cases in 3D emit `distance`, `sqrDistance`,
`parameter` and `closest[1]` — every field of the result that upstream computes
correctly on that path — and are bit-identical on all 2000 deep-run records
each; `closest[0]` is the one corrupted field and it has the deviation cases to
itself. This is the arrangement group 20 used for `DistLine3OrientedBox3`.

## Not covered

Nothing. All 21 headers of the group are implemented by the port and all have
cases. No header in the group is pure data.

Four narrower notes, not gaps in entry-point coverage:

- `numClosestPairs == 2` for `DistRay3Circle3` is not covered by an exact case.
  It can only arise from the `t0 >= 0` pass-through of a two-pair line result,
  and a line result has two pairs only when two bisected critical distances are
  exactly equal — values produced by the bisection whose back-substitution the
  port deliberately computes differently (see "Comparison policy"). A first
  attempt at a `DistRay3Circle3.compute.twoPairs` case disagreed on 6 of 20
  records, 2 of them structurally (`numClosestPairs` 1 versus 2), so it was
  removed rather than hidden behind a tolerance that cannot apply to an integer
  output. The corresponding segment arm is covered:
  `DistSegment3Circle3.compute.twoPairs` reaches the tie arm of
  `SelectClosestPoint` from the branch `t0 < 0` and `t1 >= 1`, where the two
  candidates are the segment endpoints, so the output comes from two
  point-circle queries alone and is exact on 2000 of 2000 records. The
  construction is symmetric about the circle center under `X -> 2C - X`, which
  maps the circle to itself, and the segment is short next to the circle so
  that both critical parameters are far outside `[0,1]` and the branch choice
  does not depend on their last bits.
- `DistPointHyperellipsoid`'s documented ill-conditioning near the
  smallest-extent axis (issue
  [#424](https://github.com/gradientspaceai/gtengine-js/issues/424)) is
  *preserved* in the port, not fixed, so there is no deviation case: the port
  runs upstream's bisection unchanged and the four cases agree bit for bit,
  including on the inputs where both are equally wrong.
- `DistPoint2Parallelogram2.compute.degenerate` exists because
  `Parallelogram2`'s constructor asserts `DotPerp(axis[0], axis[1]) > 0`. It
  feeds left-handed and exactly collinear axis pairs; the C++ build throws on
  every record and the port throws on the same records, which is the parity
  ORACLE.md rule 9 asks for. The ordinary `compute` case swaps a left-handed
  pair and redraws a degenerate one, so it never hits the assert.
- `DistPointHyperplane` is an `N`-dimensional template and is covered in 2D, 3D
  and 4D. `DistPointHyperellipsoid` is covered in 2D and 3D; its `N - 1`-sized
  intermediate vectors make `N = 2` the interesting edge case and it is
  included.

## Upstream bug suspects

None new. The three defects exercised above (the `DistLine3OrientedBox3`
`closest[0]` double transform, `DistLine3Circle3::Finalize`'s unchecked
normalization and the `t = tau + s` back-substitution) are already recorded in
`docs/UPSTREAM-FINDINGS.md`, as is the preserved `DistPointHyperellipsoid`
conditioning defect.

Two observations worth recording, neither a new defect:

- `Parallelogram2`'s right-handedness assert makes every singular axis pair
  unusable, so `DistPoint2Parallelogram2`'s `Inverse(B) == 0` path (which
  upstream guards against by returning the zero matrix, giving `Z = 0` and the
  parallelogram center as the closest point) is unreachable from the upstream
  API. It is not a defect, but it means the guard is dead code.
- `DistRay2Arc2` and `DistSegment2Arc2` sort their three and four candidate
  items with `std::sort` on a comparator that only orders `sqrDistance`, so the
  relative order of equal-distance items is unspecified; the port uses
  `Array.prototype.sort`, which is stable. For ranges this small the MSVC
  implementation is an insertion sort, which is stable too, and the deep run
  found no disagreement in 4000 records. The port note in
  `src/DistRay2Arc2.ts` already records the reliance.
