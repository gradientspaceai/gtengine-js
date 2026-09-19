# Group 33 (`v33-intersection`) — C++ oracle report

68 cases, 20 records each in `oracle/golden/v33-intersection.txt` (1360 records).
59 of the cases are ordinary comparisons, 9 are deliberate `deviation` cases and
two of the ordinary cases (`IntrEllipsoid3Ellipsoid3.test.bracketAssert`,
`IntrOrientedBox3Cylinder3.test.infinite`) record a C++ assertion failure that
the port must reproduce. The 59 ordinary cases carry 3703 floating-point outputs
in the goldens, **all 3703 bit-identical** to the MSVC build of upstream GTE.

Deep run `npm run oracle:deep -- 2000 v33`: 68 cases, 136000 records, 4000
recorded C++ exceptions (the two throw-parity cases), 360527 floating-point
outputs compared over the 59 ordinary cases, **100.000 % bit-identical**
(worst scaled error 0), and every deviation case deviating. All discrete
outputs (`numIntersections`, `isLine`, `isInterior`, the ellipsoid
`Classification`, the segment-mesh record count and index pairs) are compared
exactly and agreed on every record of every ordinary case.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise`.

## Coverage

Every computational entry point of every header in the group is covered: each
`TIQuery::operator()` and `FIQuery::operator()`, and each protected `DoQuery`
helper that the header declares and the port exports as a free function (called
directly on both sides through a small `Expose` subclass in C++). The headers
that only *inherit* a `DoQuery` from an aligned-box or line base
(`IntrRay2OrientedBox2`, `IntrRay3OrientedBox3`, `IntrSegment2OrientedBox2`,
`IntrSegment3OrientedBox3`) declare none of their own; the port likewise reuses
the aligned-box module functions, which groups 31 and 32 cover directly.

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `IntrEllipsoid3Ellipsoid3.h` | `test`, `test.concentric`, `test.equalEigenvalues`, `test.equalEigenvaluesDeviation`, `test.bracketAssert` | exact (discrete outputs; libm inside) (+1 deviation) | pass |
| `IntrHalfspace3Ellipsoid3.h` | `test`, `test.tangent` | exact | pass |
| `IntrLine2Arc2.h` | `test`, `find`, `find.throughLatticePoint` | exact | pass |
| `IntrLine2SegmentMesh2.h` | `find`, `find.collinear` | exact | pass |
| `IntrLine3Ellipsoid3.h` | `test`, `find`, `doQuery.fi`, `find.latticeTangent` | exact | pass |
| `IntrOrientedBox3Cylinder3.h` | `test`, `test.infinite`, `test.edgeTypoDeviation` | exact (+1 deviation) | pass |
| `IntrPlane3Plane3.h` | `test`, `find`, `find.parallel` | exact | pass |
| `IntrPlane3Triangle3.h` | `test`, `find`, `find.onPlane` | exact | pass |
| `IntrRay2Circle2.h` | `test`, `find`, `doQuery.fi`, `find.latticeTangent` | exact | pass |
| `IntrRay2OrientedBox2.h` | `test`, `find`, `find.throughPoint` | exact | pass |
| `IntrRay3Capsule3.h` | `test`, `find`, `doQuery.fi`, `test.clampDeviation`, `doQuery.capTangentDeviation` | exact (+2 deviations) | pass |
| `IntrRay3Cylinder3.h` | `find`, `doQuery.fi`, `find.latticeAimed`, `find.infiniteDeviation` | exact (+1 deviation) | pass |
| `IntrRay3OrientedBox3.h` | `test`, `find`, `find.throughPoint` | exact | pass |
| `IntrSegment2Circle2.h` | `test`, `find`, `doQuery.fi`, `find.latticeTangent` | exact | pass |
| `IntrSegment2OrientedBox2.h` | `test`, `find`, `find.throughPoint`, `find.originBox`, `find.degenerate`, `find.frameDeviation` | exact (+1 deviation) | pass |
| `IntrSegment3Capsule3.h` | `test`, `find`, `doQuery.fi`, `doQuery.junctionDeviation` | exact (+1 deviation) | pass |
| `IntrSegment3Cylinder3.h` | `find`, `doQuery.fi`, `find.latticeAimed`, `find.degenerateSegment`, `find.infiniteDeviation` | exact (+1 deviation) | pass |
| `IntrSegment3OrientedBox3.h` | `test`, `find`, `find.throughPoint` | exact | pass |
| `IntrSphere3Frustum3.h` | `test`, `test.clampDeviation` | exact (+1 deviation) | pass |

Only `IntrEllipsoid3Ellipsoid3` reaches the C math library on the compared
path (`SymmetricEigensolver3x3` uses `atan2`/`cos`/`sqrt`, `RootsBisection`
isolates the roots of `f(s)`). Its two outputs are discrete, so the case is
still declared exact and agreed on every one of the 6000 non-throwing deep-run
records. Everywhere else the `sin`/`cos` used to build frames, arcs and unit
directions, and the `sqrt` inside `Normalize`, are applied to *unrecorded*
draws; only the resulting frame, endpoint or unit vector is recorded as an
input, so libm never enters the compared computation.

### Generators and the branches they reach

Each generator alternates an exactly representable small-lattice mode with a
uniform mode; most add a constructed mode aimed at a specific branch. The
distributions below are from the 2000-record deep run.

- **Circles and arcs.** The circle of radius 5 about an integer centre carries
  twelve lattice points. With an axis-parallel unit direction and a lattice
  offset `a` in `{0, +-3, +-4, +-5}` the discriminant is exactly `25 - a^2`, so
  the roots are integers, the intersection points are exactly those lattice
  points and `Arc2::Contains` evaluates `DotPerp` of integer vectors — exactly
  zero when the hit is an arc endpoint. `a = +-5` is the tangent branch.
  `IntrLine2Arc2.find.throughLatticePoint` reports 1 or 2 arc hits on 1470 of
  2000 records against 957 for the uniform variant, and its
  `numIntersections = 1` count (837) is where exactly one of the two circle
  points lies on the arc. `IntrRay2Circle2.find.latticeTangent` and
  `IntrSegment2Circle2.find.latticeTangent` place a root exactly at the ray
  origin (respectively at the segment centre), so the `[0,+infinity)` and
  `[-e,e]` interval clips are evaluated at exact equality; they hit on 1574 and
  1822 of 2000 records against 635 and 774 for the uniform variants.
- **Oriented boxes.** The `throughPoint` cases use a lattice box (integer
  centre and extents, axes a signed permutation of the coordinate axes) and aim
  at a boundary point whose box coordinates are `+-E[f]` on one axis and
  `(k/2)*E[i]` on the others, so every world/box change of basis is exact and
  the Liang-Barsky `<=` tests are evaluated at exact equality. The hit rate
  rises from 780 to 1412 (2D ray), 341 to 1325 (3D ray), 592 to 2000 (3D
  segment) and 976 to 2000 (2D segment) of 2000 records. Single-point contacts
  occur on 316 (2D ray), 377 (3D ray) and 135 (3D segment) records. The 2D
  *segment* query cannot report `numIntersections = 1`: upstream's
  `IntrSegment2AlignedBox2::DoQuery` explicitly rewrites a collapsed clip to 2
  ("ensure the caller computes 2 points of intersection for a degenerate line
  segment"), so the collapsed clip shows instead as
  `cdeParameter[0] == cdeParameter[1]`, which happens on 127 records of
  `find.throughPoint` and 95 of the uniform `find`.
- **Planes.** The lattice mode puts the plane normal on a signed coordinate
  axis with an integer constant, which makes `|Dot(N0,N1)|` exactly 0 or 1 and
  the coplanarity test exact; `IntrPlane3Plane3.find.parallel` uses
  `N1 = +-N0` with a uniform unit normal, and reaches the coplanar branch on 421
  of 2000 records. The remaining 705 records with `N1 = +-N0` are classified
  *transverse*, which is the documented upstream quirk that only about half of
  normalized doubles satisfy `Dot(N,N) == 1` (issue #465); the port reproduces
  it bit for bit, `invDet` around `1e16` included.
- **Planes versus triangles.** With the plane normal on a coordinate axis, an
  integer constant and each triangle vertex offset by an integer in `[-1,1]`,
  every signed distance is an exact integer.
  `IntrPlane3Triangle3.find.onPlane` reaches `numIntersections` 1 (451), 2
  (1307, both the interior-segment and the edge sub-branch) and 3 (91, the
  triangle in the plane) against 218/891/7 for the uniform variant.
- **Ellipsoids versus lines and halfspaces.** With a signed-permutation frame,
  an integer centre and extents that are powers of two, `M = R*diag(1/e^2)*R^T`
  and `MInverse` are exact. `IntrLine3Ellipsoid3.find.latticeTangent` makes the
  discriminant exactly `(1 - k^2/4)/e_i^2`, so `k = +-2` is the tangent branch
  (790 of 2000 records report `numIntersections = 1`, against 53 for the
  uniform variant) and `k = 0` gives integer parameters.
  `IntrHalfspace3Ellipsoid3.test.tangent` makes `tmax` exactly zero on a third
  of its records.
- **Ellipsoid versus ellipsoid.** The uniform and lattice modes reach all four
  classifications (separated 1150, intersecting 833, ellipsoid0 contains
  ellipsoid1 6, ellipsoid1 contains ellipsoid0 11). `test.concentric` forces
  `K2 == 0`, so the closed-form common-centre path is taken (it reaches the
  other three classifications, 1727/136/137). `test.equalEigenvalues` uses two
  identity frames and power-of-two extents, which makes `M2` exactly diagonal
  and therefore two or three of the `d`-values exactly equal: that is what
  reaches the `d0 > d1 = d2`, `d0 = d1 > d2` and `d0 = d1 = d2` branches of the
  valid-pair analysis and with them the one- and two-argument `GetRoots`
  overloads, which the three-argument one would otherwise hide.
- **Cylinders.** `find.latticeAimed` uses a coordinate-axis cylinder of radius
  5 with an even integer height and aims at a lattice point of the wall circle
  at height `(k/2)*h`, so the wall roots and the `|z| <= h/2` clip are exact and
  `k = +-2` puts the hit exactly on an end disk: 869 (ray) and 1063 (segment) of
  2000 records hit, against 347 and 518 for the uniform variants.
  `IntrSegment3Cylinder3.find.degenerateSegment` drives the documented
  zero-length-segment path, where both sides report `intersect` with NaN
  parameters on 954 of 2000 records and the recorded bits agree.
- **Capsules.** The main `find` and `doQuery` cases reject, by rejection
  sampling against an independent reference interval, the inputs on which the
  delegate `IntrLine3Capsule3` is defective (below); they still reach
  `numIntersections` 1 (5-9 records) and 2 (about 610). The TI cases draw the
  capsule radius as a multiple of the computed distance, so both answers occur
  about equally often (1033/967 for the ray, 974/1026 for the segment).
- **Segment meshes.** All four `SegmentMesh2` constructors are used
  (`DISJOINT`, `CONTIGUOUS_OPEN`, `CONTIGUOUS_CLOSED`, `INDEXED`); the indexed
  mode also produces degenerate segments with equal index pairs, which is the
  documented zero-length-segment path of `IntrLine2Segment2`. The uniform case
  reports 0 to 7 intersection records per line. `find.collinear` puts four of
  the six vertices on the query line itself, which reaches the collinear branch
  exactly and records upstream's `+-max()` line-parameter sentinel; it reports
  2 to 12 records per line, and records with equal line parameters (two
  segments meeting at a shared vertex, and the two endpoints of a coincident
  segment) occur throughout, which is exactly where upstream's unstable
  `std::sort` leaves the order unspecified.
- **Sphere versus frustum.** The test point is built in frustum coordinates so
  that every Voronoi region of `DistPoint3Frustum3` is reachable, and the
  sphere radius is drawn as a multiple of the computed distance: 1167 of 2000
  records report an intersection.

## Port defects fixed

None in the port's own transcription: every value disagreement found by the
oracle was a deliberate, documented port fix of an upstream defect (nine of
them, listed below), and no case needed a tolerance — all 360527 floating-point
outputs of the 59 ordinary cases are bit-identical over the deep run.

One *new upstream defect* was found while building the
`test.equalEigenvalues` case and is now fixed in `src/`, per PORTING.md's rule
that result-corrupting upstream defects are corrected rather than preserved.

### `IntrEllipsoid3Ellipsoid3`: the `d0 > d1 = d2` fold (result-corrupting)

The valid-pair analysis folds the coefficients of
`f(s) = sum_i d_i*c_i/(d_i*s - 1)^2 - 1` that share an eigenvalue, because two
equal `d` values give the single term `d*(c_i + c_j)/(d*s - 1)^2`. The branches
for `d0 = d1 > d2` and `d0 = d1 = d2` do that correctly. The branch for
`d0 > d1 = d2` writes

```cpp
if (param[0].second > (T)0) { valid.push_back(param[0]); }
param[1].second += param[0].second;    // upstream; should be param[2].second
if (param[1].second > (T)0) { valid.push_back(param[1]); }
```

which folds in the `c` of the *distinct* eigenvalue `d0` — counting it twice,
since `param[0]` is also pushed — and never uses `param[2].second`. `f(s)` is
then built with the wrong coefficients, `GetRoots` returns the wrong roots, the
minimum and maximum squared distances are wrong and the classification is
wrong.

**Reproduction (dyadic, human-checkable).** Ellipsoid0 is the unit ball at the
origin; ellipsoid1 is axis-aligned, centred at `(1/2, 1/4, 0)`, with extents
`(1/4, 1, 1)`. Both frames are the identity and every extent is a power of two,
so `M2 = diag(16, 1, 1)` is exact and its two trailing eigenvalues are exactly
equal: the `d0 > d1 = d2` branch is taken, and the centre offset has a
component along the distinct axis (`1/2`) and one in the repeated eigenplane
(`1/4`), which is what makes the two folds differ.

- upstream (and the port before the fix): `ELLIPSOID0_CONTAINS_ELLIPSOID1`
- correct: `ELLIPSOIDS_INTERSECTING`

`(1/2, 5/4, 0)` is exactly on ellipsoid1 (`((1/2-1/2)/(1/4))^2 +
((5/4-1/4)/1)^2 = 1`) and its distance from the origin is `sqrt(29)/4 =
1.3462...`, so ellipsoid1 is plainly not contained in the unit ball.

**How the correct answer was computed.** Independently of the query: sample
ellipsoid1's surface densely (400 x 800 points) and evaluate ellipsoid0's
quadratic form `sum_i (Dot(X - C0, U0_i)/e0_i)^2 - 1` at each sample. All
values negative means ellipsoid0 contains ellipsoid1; all positive means
separated, or ellipsoid1 contains ellipsoid0 when `C0` is inside ellipsoid1;
mixed signs means intersecting. Over a grid of 1620 configurations that reach
this branch (ellipsoid0 a ball of radius 1 or 2, ellipsoid1 with extents
`(bd, br, br)` drawn from powers of two with `bd < br`, centre offset
`(dx, dy, 0)` over nine dyadic values each), **upstream's expression disagrees
with the sampled classification on 83 and the corrected expression on none**,
counting only configurations whose sampled extremes are at least `1e-2` away
from zero.

**Fix.** `src/IntrEllipsoid3Ellipsoid3.ts` now writes
`param[1][1] += param[2][1];`, with a `KNOWN UPSTREAM DEFECT (fixed here)` note
in the file header. `test/IntrEllipsoid3Ellipsoid3.test.ts` gained two
regression tests — the pinned dyadic reproduction above, and the 1620-point
grid compared against the surface sampling — both of which fail on upstream's
expression and pass on the corrected one.

## Deliberate deviations demonstrated

| case | records disagreeing (golden / deep) | record of the decision |
| --- | --- | --- |
| `IntrEllipsoid3Ellipsoid3.test.equalEigenvaluesDeviation` | 20/20, 2000/2000 | the `d0 > d1 = d2` fold above; `src/IntrEllipsoid3Ellipsoid3.ts` header note |
| `IntrSegment2OrientedBox2.find.frameDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrSegment2OrientedBox2.h` FIQuery; issue [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrOrientedBox3Cylinder3.test.edgeTypoDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrCanonicalBox3Cylinder3.h` `DoQueryNoZeros` `(U1,-D)` sign typo; issue [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrRay3Cylinder3.find.infiniteDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` cylinder headers, missing `IsFinite` guard; issues [#197](https://github.com/gradientspaceai/gtengine-js/issues/197), [#206](https://github.com/gradientspaceai/gtengine-js/issues/206), [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrSegment3Cylinder3.find.infiniteDeviation` | 20/20, 2000/2000 | the same |
| `IntrRay3Capsule3.test.clampDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `DistRaySegment.h` parallel branch and regions 1 and 5; issue [#126](https://github.com/gradientspaceai/gtengine-js/issues/126) |
| `IntrRay3Capsule3.doQuery.capTangentDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrLine3Capsule3.h` hemisphere roots; issue [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) |
| `IntrSegment3Capsule3.doQuery.junctionDeviation` | 15/20, 1531/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrLine3Capsule3.h` cap-junction plane; issue [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) |
| `IntrSphere3Frustum3.test.clampDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `DistPoint3Frustum3.h` LF/UF edge cases; issue [#421](https://github.com/gradientspaceai/gtengine-js/issues/421) |

Each deviation is confined to inputs on which upstream is actually defective,
and the corresponding main case is generated by rejection sampling against an
explicit predicate (or, where the predicate is a closed form, by construction)
rather than by narrowing the input ranges:

- **`IntrEllipsoid3Ellipsoid3`.** The probe replicates upstream's own control
  flow down to the valid-pair analysis — `K2`, `M2`, the eigensolver, the
  `std::greater` sort — and reports whether the record reaches the
  `d0 > d1 = d2` branch with `param[0].second != param[2].second`. That is an
  *exact* characterisation: both expressions start from the same
  `param[1].second`, so the two sums agree bit for bit exactly when the added
  values do. `test` and `test.equalEigenvalues` reject it (and still reach all
  four classifications, and all three equal-eigenvalue sub-branches). The
  deviation case additionally keeps only the records on which upstream's
  classification disagrees with the dense surface sampling by a margin of at
  least `1e-2`, so every one of its records is a genuine misclassification
  rather than a different intermediate: on the 2000 deep-run records upstream
  reports `ELLIPSOID0_CONTAINS_ELLIPSOID1` 1260 times, `ELLIPSOIDS_SEPARATED`
  637 times and `ELLIPSOIDS_INTERSECTING` 103 times, where the port (and the
  sampling) report `ELLIPSOIDS_INTERSECTING` 1897 times and
  `ELLIPSOID1_CONTAINS_ELLIPSOID0` 103 times.
- **`IntrSegment2OrientedBox2`.** The FI query builds
  `result.point[i] = box.center + (segOrigin + parameter[i] * segDirection)`,
  adding a box-frame vector to the world-space box centre; the port evaluates
  the world-space centred form instead. Upstream is *mathematically* right only
  when the box axes are the coordinate axes, and even then its expression adds
  and subtracts `box.center`, which rounds. The set on which the two are
  bit-identical is therefore `box.center == 0` with the identity axes, which
  is what `find.originBox` uses: it compares every field including `point`.
  `find` covers a general rotated box on every other field (`intersect`,
  `numIntersections`, both `parameter` and both `cdeParameter`),
  `find.throughPoint` does the same for the exact lattice construction, and
  `find.degenerate` compares the zero-extent branch in full, because the points
  it stores are the segment endpoints themselves and involve no change of
  basis. `find.frameDeviation` sends a segment through the centre of a rotated
  box, so every record has a reported point and every record deviates.
- **`IntrOrientedBox3Cylinder3`.** The query is a thin wrapper around
  `IntrCanonicalBox3Cylinder3`, so it inherits the `(U1,-D)` sign typo of
  `DoQueryNoZeros` and the port's correction of it. The case file carries the
  same corrected `DoQueryNoZeros` that `oracle/cpp/cases/v31-intersection.cpp`
  and `v32-intersection.cpp` use, applied to the transformed cylinder and the
  canonical box the wrapper builds; the main case keeps the inputs where the
  corrected answer equals upstream's, the deviation case those where it does
  not.
- **`IntrRay3Cylinder3`, `IntrSegment3Cylinder3`.** `Cylinder3` marks an
  infinite cylinder with the sentinel `height = -1`, which upstream feeds into
  the finite-cylinder formulas (an empty slab, `-h/2 = 0.5 > -0.5 = h/2`); the
  port asserts. The deviation cases pass that sentinel, so the port throws where
  upstream returns. Upstream's answer is not even constant: it reports an
  intersection on 200 (ray) and 320 (segment) of 2000 records.
- **`IntrRay3Capsule3` TI.** The query is
  `DistRaySegment(ray, capsule.segment).distance <= capsule.radius`, and
  upstream never clamps the ray parameter to `s0 >= 0` in its parallel branch
  and in regions 1 and 5, so it measures to a point behind the ray origin. The
  port clamps, which restricts the feasible set and can therefore only increase
  the distance. The main case keeps only the inputs for which upstream reports
  `parameter[0] >= 0` — exactly where the clamp is a no-op, the same predicate
  group 19 uses for the distance query itself — and draws the radius as a
  multiple of the distance. The deviation case biases the segment behind the
  ray origin, keeps the records with `parameter[0] < 0`, and sets the radius to
  *exactly* upstream's distance: upstream then reports an intersection
  (`d <= d`) and the port, whose distance is strictly larger, reports none.
- **`IntrRay3Capsule3` / `IntrSegment3Capsule3` FI.** Both delegate to
  `IntrLine3Capsule3`, whose `intersect` flag is never set when only one
  hemisphere root is accepted and whose overlapping acceptance regions accept a
  cap-junction root twice; the port fixes both. The main `find` and
  `doQuery.fi` cases reject, with group 31's reference interval (the min and
  max over the boundary candidates of the union of the finite cylinder and the
  two end balls), the inputs on which upstream's delegate disagrees with that
  reference — the segment cases probe with the centred form the query itself
  will compute from the recorded endpoints, not with the draws they were built
  from. The two deviation cases are group 31's constructions carried through
  the ray and segment clips: a ray tangent to a cap pole (upstream leaves
  `intersect` false, so the ray clip, guarded on that flag, drops the contact)
  and a segment meeting the wall exactly on the junction circle `z = -e`
  (upstream collapses the interval to a point). 469 of the 2000 junction
  records agree anyway, and the construction says exactly which: the second
  wall root lies at capsule height `-e + 2*sqrt(3)*r`, which is inside the wall
  band `|z| <= e` precisely when `r <= e/sqrt(3)`. For those `(r,e)` pairs — 2
  of the 9 the generator draws, so about 22 % of the records — the count
  reaches two from the wall alone, upstream's early return is correct and both
  sides report the same interval.
- **`IntrSphere3Frustum3`.** The query is
  `DistPoint3Frustum3(sphere.center, frustum).distance <= sphere.radius`, and
  two of upstream's ten far-edge assignments copy the test coordinate straight
  through instead of clamping it to the end of the edge, so upstream's closest
  point can lie outside the frustum and its distance is smaller than the true
  one. The main case rejects group 20's exact predicate for reaching one of
  those two assignments with an out-of-range coordinate; the deviation case
  keeps only it and again sets the radius to exactly upstream's distance.

Upstream defects recorded for this group's headers that the port *preserves*
need no deviation case; the oracle compares them bit for bit and they agree:
`IntrPlane3Plane3`'s exact `|Dot(N0,N1)| >= 1` parallel test (#465, see the
`find.parallel` histogram above), `IntrPlane3Triangle3`'s uninitialised
`Real s[3]` (#255, a correction with no behavioural effect as used),
`IntrLine2SegmentMesh2`'s `+-max()` sentinel copied into `lineParameter` while
the real endpoint is stored as `point` (#461),
`IntrSegment2OrientedBox2`'s `cdeParameter = cdeParameter` no-op and its
`numIntersections = 2` for a contained point segment where the aligned sibling
reports 1 (#255), `IntrSegment3Cylinder3`'s missing zero-length-segment guard
(#197, the `find.degenerateSegment` case), `IntrEllipsoid3Ellipsoid3`'s dead
`Matrix3x3 D0` (#255) and its reachable `GetRoots` bracketing asserts (#255,
#461, the `test.bracketAssert` throw-parity case; its `d0 > d1 = d2` fold is
fixed, not preserved — see above), and `Arc2::Contains`'s
negative-epsilon comment (#155, not on the one-argument path these queries use).

## Not covered

Nothing. All 19 headers of group 33 are implemented by the port and every entry
point has a case. Two deliberate restrictions are worth recording:

- `IntrEllipsoid3Ellipsoid3`'s main `test` and `test.equalEigenvalues` cases
  reject, by running the query in a `try`/`catch` on the C++ side, the
  configurations on which upstream's ad-hoc `epsilon = 0.001` bracketing
  asserts fire. Those inputs are not "upstream is defective and the port is
  fixed" — the port preserves the throw — but the sign of `F` at a bracket
  endpoint there is decided by round-off, so whether the assert fires is not a
  property MSVC and V8 can be expected to agree on. The dedicated
  `test.bracketAssert` case pins the one reproduction the port's own test file
  pins, where both sides throw; it agreed on all 2000 deep-run records.
- The `find.throughPoint` and `find.originBox` cases use an oriented box whose
  axes are a signed permutation of the coordinate axes (respectively the
  identity). That is the only way to make the world/box change of basis exact,
  which is what puts the Liang-Barsky comparisons on exact equality; the
  uniform `find` and `test` cases use fully rotated frames.

## Upstream bug suspects

None new for this group's headers. Two observations worth recording:

- `IntrSegment2AlignedBox2::DoQuery` (group 32's header, reached here through
  `IntrSegment2OrientedBox2`) rewrites `numIntersections == 1` to `2` with the
  comment "ensure the caller computes 2 points of intersection for a degenerate
  line segment representing a single point". The rewrite is unconditional, so a
  segment that touches the box in a single point reports two coincident points,
  whereas the `Ray2`/`Line2` siblings report one. The port reproduces it; it is
  consistent with the already-recorded quirk that a contained point segment
  reports `numIntersections = 2` (#255), and it is a documentation trap rather
  than a wrong answer.

The `IntrEllipsoid3Ellipsoid3` `d0 > d1 = d2` fold, found by this group, is
result-corrupting and is therefore *fixed* in the port rather than listed as a
suspect; see "Port defects fixed" above for the reproduction, the independent
reference and the regression tests. It needs a new
`docs/UPSTREAM-FINDINGS.md` entry (RC, fixed); it is *not* the already-recorded
`GetRoots` bracketing problem.
