# Group 31 (`v31-intersection`) — C++ oracle report

61 cases, 20 records each in `oracle/golden/v31-intersection.txt` (1220 records).
1100 records are ordinary comparisons and carry 2644 floating-point outputs;
2627 of them are bit-identical to the MSVC build of upstream GTE. The 17 that
are not all belong to `IntrLine3Torus3.find`, the one case whose compared path
calls the C math library (worst scaled error 3.7e-14 in the goldens). The
remaining 120 records belong to the six `deviation` cases, which are expected
to disagree.

Deep run `npm run oracle:deep -- 2000 v31-intersection`: 61 cases, 122000
records, 0 C++ exceptions, 258255 floating-point outputs compared over the 55
ordinary cases, **99.11 % bit-identical**, and every case within its declared
comparison. Only two cases are not bit-identical on every output:

| case | bit-identical | worst scaled error | why |
| --- | --- | --- | --- |
| `IntrLine3Torus3.find` | 5247 / 7548 | 1.72e-12 | `RootsPolynomial::SolveQuartic` (cbrt, trigonometric resolvent) and `Torus3::GetParameters` (atan2) |
| `IntrCylinder3Cylinder3.test` | 5998 / 6000 | 2.22e-16 (1 ulp) | `cos`/`sin` of the sampled hemisphere directions |

The other 53 ordinary cases are bit-identical on all 244707 of their
floating-point outputs over 110000 deep-run records.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise`.

## Coverage

Every computational entry point of every header in the group is covered: each
`TIQuery::operator()` and `FIQuery::operator()` overload, and each protected
`DoQuery` helper (the port exports these as free functions, so they are called
directly on both sides through a small `Expose` subclass in C++).

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `IntrArc2Arc2.h` | `find` | exact | pass |
| `IntrCanonicalBox3Cylinder3.h` | `test`, `test.oneZero`, `test.twoZeros`, `test.edgeTypoDeviation` | exact (+1 deviation) | pass |
| `IntrCapsule3Capsule3.h` | `test` | exact | pass |
| `IntrCircle2Arc2.h` | `find` | exact | pass |
| `IntrCylinder3Cylinder3.h` | `test`, `test.parallel`, `test.infiniteDeviation` | tol 1e-12 / exact (+1 deviation) | pass |
| `IntrDisk2Sector2.h` | `test` | exact | pass |
| `IntrHalfspace3Capsule3.h` | `test` | exact | pass |
| `IntrHalfspace3Cylinder3.h` | `test.perpendicular`, `test.rootDeviation` | exact (+1 deviation) | pass |
| `IntrLine2Circle2.h` | `test`, `find`, `doQuery.fi`, `find.tangent` | exact | pass |
| `IntrLine2OrientedBox2.h` | `test`, `find` | exact | pass |
| `IntrLine2Ray2.h` | `test`, `find` | exact | pass |
| `IntrLine2Segment2.h` | `test`, `find` | exact | pass |
| `IntrLine3Capsule3.h` | `test`, `find`, `doQuery.fi`, `doQuery.capTangentDeviation`, `doQuery.junctionDeviation` | exact (+2 deviations) | pass |
| `IntrLine3Cylinder3.h` | `find`, `doQuery.fi`, `find.axisAligned` | exact | pass |
| `IntrLine3OrientedBox3.h` | `test`, `find` | exact | pass |
| `IntrLine3Torus3.h` | `find` | tol 1e-9 (libm) | pass |
| `IntrOrientedBox2Sector2.h` | `test`, `test.clipDeviation` | exact (+1 deviation) | pass |
| `IntrOrientedBox3Frustum3.h` | `test` | exact | pass |
| `IntrRay2AlignedBox2.h` | `test`, `find`, `doQuery.ti`, `doQuery.fi` | exact | pass |
| `IntrRay2Ray2.h` | `test`, `find`, `test.collinear`, `find.collinear` | exact | pass |
| `IntrRay2Segment2.h` | `test`, `find`, `test.collinear`, `find.collinear` | exact | pass |
| `IntrRay2Triangle2.h` | `test`, `find`, `doQuery.fi`, `find.throughPoint` | exact | pass |
| `IntrRay3AlignedBox3.h` | `test`, `find`, `doQuery.ti`, `doQuery.fi` | exact | pass |
| `IntrRay3Rectangle3.h` | `test`, `find`, `find.throughPoint` | exact | pass |

Only `IntrLine3Torus3.find` and `IntrCylinder3Cylinder3.test` reach the C math
library on the compared path. Everywhere else the `sin`/`cos` used to build
frames and sector angles, and the `sqrt` inside `Normalize`, are applied to
*unrecorded* draws; only the resulting frame, unit vector or cosine/sine pair
is recorded as an input, so libm never enters the compared computation.

### Generators and the branches they reach

Each generator alternates an exactly representable small-lattice mode with a
uniform mode; several add a constructed mode aimed at a specific branch. The
distributions below are from the 2000-record deep run.

- **Lines, rays and segments in 2D.** Lattice origins with integer directions
  make the `DotPerp(D0,D1) == 0` parallel test of `IntrLine2Line2` fire
  exactly, so the collinear "same line" branches of `IntrLine2Ray2` and
  `IntrLine2Segment2` (`numIntersections == max int32`) are reached (14 and 36
  records). The four `*.collinear` cases put both primitives on one lattice
  line: for `IntrRay2Ray2` the same-direction (999/1004 records),
  opposite-direction-overlapping (398/552) and disjoint (453/444) branches all
  occur, and the TI variant reports the `t == 0` shared-origin case as one
  intersection on 150 records where the FI variant reports two (finding #200,
  preserved and now compared against the C++ build on both sides). For
  `IntrRay2Segment2` the collinear overlap (1109/984), the single-point touch
  (87/199) and the disjoint case (804/817) all occur.
- **Rays versus boxes and triangles.** `IntrRay2AlignedBox2`,
  `IntrRay3AlignedBox3` and `IntrRay2Triangle2` reach `numIntersections`
  0, 1 and 2 for both the query and its `DoQuery` helper. `find.throughPoint`
  aims the ray at `v0 + (k1*E1 + k2*E2)/4` on an integer triangle with an
  integer direction and origin `target - t*direction`, `t` a multiple of 1/2 in
  `[-2,3]`, so the clip against `[0,+infinity)` is evaluated at exactly zero
  (ray origin on the triangle) and the barycentric comparisons at exact
  equality.
- **Rays versus rectangles.** The `throughPoint` variant aims at
  `C + (a/2)*e0*W0 + (b/2)*e1*W1` with `a, b` in `[-3,3]` on an axis-aligned
  lattice rectangle, so the corner and edge cases of `|rectCoord| <= extent`
  are exact; 529 of 2000 records hit, against 145 for the uniform variant.
- **Lines versus circles.** `find.tangent` puts the line origin at
  `C + (R+d)*e_i + s*e_j` with the axis direction `e_j`, which makes the
  discriminant exactly `R^2 - (R+d)^2`: the tangent branch
  (`numIntersections == 1`) occurs on 622 of 2000 records, against 212 for the
  uniform variant.
- **Lines versus cylinders.** `find.axisAligned` uses a coordinate-axis
  cylinder axis and a coordinate-axis line direction, which reaches the
  otherwise unreachable parallel branch (`|Dot(W,D)| == 1`) and perpendicular
  branch (`D[2] == 0`) exactly, and places the line at radial distance
  `r + d`, `d` in `{-1,0,1}`, so the wall discriminant is exactly zero on a
  third of the records.
- **Arcs.** Mode 0 puts both arcs on the lattice circle of radius 5 about a
  common centre, whose 12 integer points make the arcs cocircular and their
  endpoints compare exactly. All eight `IntrArc2Arc2::Configuration` values
  occur (`NO_INTERSECTION` 805, `COCIRCULAR_ONE_ARC` 566, `TWO_ARCS` 161,
  `NONCOCIRCULAR_ONE_POINT` 183, `TWO_POINTS` 120, `COCIRCULAR_ONE_POINT` 65,
  `COCIRCULAR_TWO_POINTS` 6, `ONE_POINT_ONE_ARC` 94). `IntrCircle2Arc2`
  reaches all four `numIntersections` values, including the `max int32`
  "arc is on the circle" sentinel (1000 records).
- **Sectors.** The cosine and sine of the half-angle are recorded directly
  (mode 0 uses the exact pairs `(0,1)`, `(0.6,0.8)` and `(0.8,0.6)`), so the
  queries are exact and the half-angle stays inside the documented `(0,pi/2]`
  convexity precondition. `IntrDisk2Sector2` reports both answers on roughly
  40/60 of records; `IntrOrientedBox2Sector2` on 46/54.
- **Cylinders.** `IntrCylinder3Cylinder3.test.parallel` uses coordinate-axis
  directions, which makes `Length(Cross(W0,W1))` exactly zero and reaches the
  parallel branch, the coincident-origin early return and both parallel
  separating axes. The general case draws `numTheta` in `[2,6]` and `numPhi`
  in `[2,5]` and always passes `numThreads = 1` (see below); it reports
  separation on 961 of 2000 records, both from the four analytic axes and from
  the sampled hemisphere directions.
- **Canonical box versus cylinder.** The three main cases fix the number of
  zero components of the cylinder axis direction at 0, 1 and 2, which selects
  `DoQueryNoZeros`, `DoQueryOneZero` and `DoQueryTwoZeros` respectively; each
  reports both answers on roughly half its records.
- **Box versus frustum.** The box centre is placed relative to the frustum
  frame at `E + (dMin + t*(dMax-dMin))*D + a*u*U + b*r*R` with `t` in
  `[-2,3]` and `a, b` in `[-6,6]`, which gives 369 intersecting records (on
  which every one of the 26 separating axes is evaluated, since the query
  returns `true` only after all of them pass) and 1631 separated records
  spread over the early rejections.

## Port defects fixed

**`src/IntrIntervals.ts`: `Math.max`/`Math.min` are not `std::max`/`std::min`
for signed zeros.** Upstream clips an interval against a semi-infinite one with
`std::max(finite[0], a)` and `std::min(finite[1], a)`, which are specified as
`(a < b ? b : a)` and `(b < a ? b : a)` and therefore return the *first*
argument when the two compare equal. `Math.max` and `Math.min` instead order
`-0` below `+0`, so `Math.max(-0, +0)` is `+0` where `std::max(-0, +0)` is
`-0`, and `Math.min(+0, -0)` is `-0` where `std::min(+0, -0)` is `+0`.

The oracle caught it through `IntrRay2Triangle2`, whose `DoQuery` clips the
line-triangle parameter interval against `[0,+infinity)`: for a ray whose first
contact parameter is a negative zero the port returned `+0` where the C++ build
returned `-0`. 15 of the 6000 deep-run records of the three
`IntrRay2Triangle2` cases disagreed. The fix replaces the four call sites with
module-private `stdMax`/`stdMin` helpers that use upstream's comparison form;
`test/IntrIntervals.test.ts` gains a regression test that pins the sign of zero
on all four of them. The size of the error is zero in magnitude but the sign of
zero is observable: it is returned to the caller as a ray or segment parameter
and it propagates through `origin + t * direction`.

No other disagreement in the group was a port defect.

## Deliberate deviations demonstrated

| case | records disagreeing (golden / deep) | record of the decision |
| --- | --- | --- |
| `IntrOrientedBox2Sector2.test.clipDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrOrientedBox2Sector2.h` boundary clipping; issue [#200](https://github.com/gradientspaceai/gtengine-js/issues/200) |
| `IntrHalfspace3Cylinder3.test.rootDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrHalfspace3Cylinder3.h` `root` computation; issue [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrLine3Capsule3.doQuery.capTangentDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrLine3Capsule3.h` hemisphere roots; issue [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) item 2 |
| `IntrLine3Capsule3.doQuery.junctionDeviation` | 17/20, 1555/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrLine3Capsule3.h` cap-junction plane; issue [#461](https://github.com/gradientspaceai/gtengine-js/issues/461) item 3 |
| `IntrCylinder3Cylinder3.test.infiniteDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` cylinder headers, missing `IsFinite` guard; issue [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrCanonicalBox3Cylinder3.test.edgeTypoDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrCanonicalBox3Cylinder3.h` `DoQueryNoZeros` `(U1,-D)` sign typo; issue [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |

Each deviation is confined to inputs on which upstream is actually defective,
and the corresponding main case is generated by rejection sampling against a
predicate that makes the restriction explicit rather than by narrowing the
input ranges:

- **`IntrOrientedBox2Sector2`.** `IntrHalfspace2Polygon2` reports
  `intersect = true` with an *empty* polygon when no clipping is necessary, and
  upstream's `polygon = std::move(hpResult.polygon)` then throws the working
  polygon away. The probe `ProbeOb2S2` replicates upstream's control flow and
  reports whether either clip is such a no-op. The main case keeps only the
  inputs where neither is (upstream sound); the deviation case keeps only the
  inputs where one is *and* keeping the polygon changes the answer. The second
  clip is the reachable one: its `prjcen1 >= -radius1` guard is computed from
  the original box while the polygon being clipped is the output of the first
  clip, so the no-op happens on an open set of configurations, not only on a
  boundary.
- **`IntrHalfspace3Cylinder3`.** Upstream's
  `root = sqrt(max(1, 1 - Dot(N,W)^2))` always selects 1, so the documented
  `r*sqrt(1-Dot(N,W)^2)` term is never computed. Upstream is sound exactly when
  `Dot(N,W)` is zero, so the main case makes the halfspace normal orthogonal to
  the cylinder axis — in lattice mode by construction from coordinate axes
  (`Dot` exactly 0), in uniform mode from an orthonormal frame (`Dot` within
  one ulp of 0, for which `1 - Dot^2` rounds to exactly 1). The deviation case
  places the halfspace constant so that upstream's overestimated interval
  reaches the plane and the corrected one does not.
- **`IntrLine3Capsule3`.** The case file carries an independent reference for
  the intersection interval: the solid capsule is the union of the solid finite
  cylinder and the two solid end balls, a line meets a convex set in an
  interval, and every endpoint of that interval is a root of one of the three
  quadratics, so the interval is `[min,max]` over the candidate roots. Upstream
  accumulates roots one at a time, returns as soon as two are accepted and sets
  `intersect` only in the cylinder-wall branches, so it disagrees with the
  reference on two families of inputs. The `find` and `doQuery.fi` generators
  reject those (and the zero-length capsule segment, whose frame is undefined —
  a documented upstream limitation the port preserves, and on which both sides
  produce NaN). The two deviation cases construct the families exactly: a line
  tangent to the bottom end sphere at its pole (the only accepted root comes
  from a hemisphere, so upstream leaves `intersect` false while copying the
  parameter), and a line meeting the wall exactly on the cap-junction circle
  `z = -e` (the wall test and the bottom-ball test both accept that one root,
  the count reaches two and the early return collapses the interval to a
  point). The junction case deviates on 17 of 20 golden records: for the three
  `(radius, extent)` combinations with `3.46*r - e <= e` the second wall root
  is also accepted, the wall branch returns with both true endpoints and
  upstream is right.
- **`IntrCanonicalBox3Cylinder3`.** The case file carries upstream's
  `DoQueryNoZeros` with the `(U1,-D)` sign typo corrected (the only line that
  differs is the `E[1]` sign of the second `P1`; the rest is copied verbatim).
  The three main cases keep only the inputs where the corrected answer equals
  upstream's; the deviation case keeps only the inputs where it does not.
  Upstream reports "no intersection" on 1441 of those 2000 records.
- **`IntrCylinder3Cylinder3`.** `Cylinder3` marks an infinite cylinder with the
  sentinel `height = -1`, which upstream feeds into the finite-cylinder
  formulas; the port asserts `IsFinite` instead, so it throws where the C++
  build returns a number.

All other upstream defects recorded for this group's headers are *preserved* by
the port and therefore need no deviation case; the oracle compares them bit for
bit and they agree: `IntrCapsule3Capsule3` non-robust segment-segment distance
(#455), `IntrLine2Ray2`/`IntrLine2Segment2` zero-direction supporting line and
the redundant `line1Parameter` test (#197), `IntrLine3Capsule3` zero-length
capsule segment (#200), `IntrLine3Cylinder3` missing infinite-cylinder branch
and unit-direction assumption (#197, #200), `IntrOrientedBox2Sector2` halfplane
wedge model for half-angles above pi/2 (#200), `IntrRay2Ray2` TI/FI
disagreement for collinear opposite rays sharing an origin (#200, #455),
`IntrCylinder3Cylinder3`'s stale `phi[j]` comment (#197).

## Not covered

Nothing. All 24 headers of group 31 are implemented by the port and every entry
point has a case. Two deliberate restrictions of the *inputs* are worth
recording:

- `IntrCylinder3Cylinder3` is always constructed with `numThreads = 1`. The
  upstream multithreaded path reports whichever thread found a separating
  direction first rather than the first direction in scan order, so its output
  is not deterministic and cannot be compared; the port documents that it
  always runs the single-threaded loop, which visits the same directions in the
  same order.
- The sector cases keep the half-angle inside the documented `(0, pi/2]`
  convexity precondition of `IntrDisk2Sector2` and `IntrOrientedBox2Sector2`.
  `Sector2`'s own default half-angle of pi violates it (finding #200,
  preserved), so it is not exercised here.

`IntrLine3Torus3` emits `numIntersections` as an exact integer: the quartic
root count agreed on all 2000 deep-run records, so the root classification of
`RootsPolynomial::SolveQuartic` did not diverge on any input generated here.
Its real outputs use a 1e-9 tolerance, justified above.

## Upstream bug suspects

None new. Everything observed in this group is already in
`docs/UPSTREAM-FINDINGS.md`.
