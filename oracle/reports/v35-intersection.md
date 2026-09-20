# Group 35 (`v35-intersection`) — C++ oracle report

56 cases, 20 records each in `oracle/golden/v35-intersection.txt` (1120 records).
48 of the cases are ordinary comparisons and 8 are deliberate `deviation`
cases. No case records a C++ exception: every generator respects the
preconditions of the queries it drives. The 48 ordinary cases carry 6847
floating-point outputs in the goldens, **all 6847 bit-identical** to the MSVC
build of upstream GTE.

Deep run `npm run oracle:deep -- 2000 v35-intersection`: 56 cases, 112000
records, 0 recorded C++ exceptions, 685081 floating-point outputs compared over
the 48 ordinary cases, **100.000 % bit-identical** (worst scaled error 0), and
all 8 deviation cases deviating on every one of their 2000 records. Every
discrete output (`numIntersections`, the line/cone result `type`, the
tetrahedron `separating` pair, the polygon vertex counts and the booleans) is
compared exactly and agreed on every record of every ordinary case.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise`.

## Coverage

Every computational entry point of every header in the group is covered: each
`TIQuery::operator()` (including the `epsilon` and moving-triangle overloads),
each `FIQuery::operator()` (including the moving-triangle overload), each
protected `DoQuery` the port exports as a free function (called directly on
both sides through a small `Expose` subclass in C++), and the one port-only
free function that upstream writes inline three times
(`intrTriangle3BoxFacePlanes`).

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `IntrRay3Plane3.h` | `test`, `find`, `doQuery.fi`, `test.throughLatticePoint`, `find.throughLatticePoint` | exact | pass |
| `IntrSegment3Plane3.h` | `test`, `find`, `doQuery.fi`, `test.throughLatticePoint`, `find.throughLatticePoint` | exact | pass |
| `IntrTriangle2Triangle2.h` | `test`, `find`, `test.sharedVertex`, `find.sharedVertex` | exact | pass |
| `IntrTetrahedron3Tetrahedron3.h` | `test`, `test.smallEdges`, `test.edgeCutoffDeviation`, `test.edgeSeparationDeviation` | exact (+2 deviations) | pass |
| `IntrTriangle3CanonicalBox3.h` | `test`, `find`, `test.onFaces`, `find.onFaces` | exact | pass |
| `IntrTriangle3AlignedBox3.h` | `test`, `find`, `test.onFaces`, `find.onFaces` | exact | pass |
| `IntrTriangle3OrientedBox3.h` | `test`, `find`, `test.onFaces`, `find.onFaces`, `facePlanes` | exact | pass |
| `IntrTriangle3Triangle3.h` | `test`, `find`, `test.onPlane`, `find.onPlane`, `test.coplanar`, `test.moving`, `find.moving`, `test.moving.coplanar`, `find.moving.coplanar`, `test.moving.parallelDeviation`, `find.moving.parallelDeviation` | exact (+2 deviations) | pass |
| `IntrSphere3Cone3.h` | `test.infiniteCone`, `test.infiniteTruncatedCone`, `test.finiteCone`, `test.coneFrustum`, `find`, `find.originVertex`, `find.vertexDeviation` | exact (+1 deviation) | pass |
| `IntrRay3Cone3.h` | `find`, `find.throughVertexDeviation` | exact (+1 deviation) | pass |
| `IntrSegment3Cone3.h` | `find`, `find.throughInterior`, `find.throughVertexDeviation` | exact (+1 deviation) | pass |
| `IntrOrientedBox3Cone3.h` | `test`, `test.staleAdjacencyDeviation` | exact (+1 deviation) | pass |

No case reaches the C math library on the compared path. The `sin`/`cos` used
to build orthonormal frames and unit directions, and the `cos`/`sin`/`tan` of
`Cone3::SetAngle`, are applied on the C++ side to *unrecorded* draws; only the
resulting frame, unit vector or trigonometric value is recorded as an input.
The remaining `SetAngle` members are exact functions of the recorded three
(`cosAngleSqr = cosAngle*cosAngle`, `sinAngleSqr = sinAngle*sinAngle`,
`invSinAngle = 1/sinAngle`), and both sides derive them with those same
operations. `Cone3::angle` is read by no query of this group and is set to zero
on both sides. Everything else in the group is `+ - * / sqrt fabs min max` and
comparisons, so every ordinary case is declared `{ exact: true }` and no case
needed a tolerance.

Two output representations are canonicalised identically on both sides, and
both are documented port conventions rather than computed values:

- `IntrTetrahedron3Tetrahedron3`'s `separating` entries use
  `std::numeric_limits<size_t>::max()` upstream and `Number.MAX_SAFE_INTEGER`
  in the port. The C++ value is not representable as a double at all, so both
  sides emit `-1` for the sentinel and the index otherwise.
- `IntrLine3Plane3`'s "the linear component lies in the plane" sentinel,
  `std::numeric_limits<int32_t>::max()`, is representable and is compared as
  is; it occurs on 398 of the 2000 `IntrRay3Plane3.find.throughLatticePoint`
  records and 1295 of the `IntrSegment3Plane3.find.throughLatticePoint` ones.

### Generators and the branches they reach

Each generator alternates an exactly representable small-lattice mode with a
uniform mode; most add a constructed mode aimed at a specific branch. The
distributions below are from the 2000-record deep run.

- **Planes.** The lattice mode puts the plane normal on a signed coordinate
  axis with an integer constant, so every signed distance of an integer point
  is an exact integer. `find.throughLatticePoint` places the ray origin at
  `Q - k*U` for an integer point `Q` of the plane, an integer direction `U` and
  an integer `k`, so the intersection parameter is exactly `k` and the ray
  query's `parameter < 0` clamp is evaluated at exact equality (939 hits, 663
  misses and 398 line-in-plane records, against 967/1013/20 for the uniform
  variant). The segment construction uses a signed coordinate axis and a
  power-of-two length `L`, which makes `GetCenteredForm` exact (the reciprocal
  `1/L` is exact and the unit direction is a coordinate axis) and the extent
  exactly `L/2`, so `|t| <= extent` is evaluated at exact equality for the
  endpoint-on-plane placements.
- **Triangles in 2D.** Both queries require counterclockwise vertices, so the
  generator computes `DotPerp` and swaps two vertices when the orientation is
  negative (a zero-area triangle has no orientation and is redrawn — a
  precondition, not a narrowing). `find` reports polygons of 3 (412), 4 (613),
  5 (207) and 6 (20) vertices and the empty set 748 times.
  `find.sharedVertex` builds triangle1 on a vertex of triangle0, so `WhichSide`
  evaluates a projection of exactly zero and the clipping query meets its
  vertex-on-line configurations.
- **Tetrahedra.** The generator enforces the (undocumented, #465) positive
  orientation precondition by swapping two vertices when
  `Dot(Cross(v1-v0,v2-v0), v3-v0)` is negative. `test` reaches all eight
  face-normal separations (1236 records) and intersection (564), and
  `test.smallEdges` additionally reaches upstream's *edge-edge* separations —
  the only regime in which upstream's unnormalised cutoff lets that phase run
  at all. The `epsilon` argument is drawn as 0 half the time and in
  `(0, 0.25)` otherwise.
- **Triangles against boxes.** The lattice mode uses integer extents and
  triangle vertices drawn as `centre + sum_j k_j * extent_j * axis_j` with
  `k_j` in `{-2..2}`, so vertices land exactly on faces, edges and corners and
  every separating-axis and clipping comparison is evaluated at exact equality;
  for the oriented box the lattice axes are a signed permutation of the
  coordinate axes, which makes the change of basis exact. The `onFaces` hit
  rates are 1691/1713/1707 of 2000 (canonical/aligned/oriented) against
  1393/934/966 for the uniform variants, and the inside polygon reaches 3 to 7
  vertices routinely, 8 in `IntrTriangle3CanonicalBox3.find` and **9** in
  `IntrTriangle3OrientedBox3.find` — which confirms on the real build the
  documentation correction recorded for these headers (upstream says "at most
  7"; the bound is 9).
- **Triangle against triangle, stationary.** `find.onPlane` puts triangle0 in
  an axis plane and offsets each vertex of triangle1 by an integer in
  `[-1,1]`, so the signed distances are exact integers and every row of the
  header's `(numNegative, numPositive, numZero)` table is reached: the
  intersection is empty on 1371 records, a point on 137, a segment on 462 and a
  polygon of 3 to 5 vertices on 30. Two of the 2000 records have a zero-area
  triangle0, which is upstream's unguarded division by `normal[lookup[2]]`
  (#334, preserved by the port): both sides emit the same NaN coordinates.
  `test.coplanar` puts both triangles in the same axis plane, so `Cross(N0,N1)`
  is exactly zero and the coplanar branch of the separating-axis test is taken.
- **Triangle against triangle, moving.** The relative velocity is aimed from
  the centroid of triangle1 at the centroid of triangle0 with noise, which
  brings the two together on 1295 of 2000 records (a purely random relative
  velocity almost never does). `test.moving.coplanar`/`find.moving.coplanar`
  use the exactly coplanar lattice pair. Both main moving cases reject the
  inputs on which upstream's parallel test and the port's disagree (below).
- **Cones.** The cone kind (infinite, infinite truncated, finite, frustum) is
  selected explicitly, and each of `IntrSphere3Cone3`'s four private helpers
  has its own case with a balanced true/false split (810/1190, 707/1293,
  560/1440 and 622/1378). In lattice mode the cone axis is a signed coordinate
  axis with an integer vertex and integer heights and the ray origin and
  direction are integers, so every height is an exact integer and the
  height-range clamps are evaluated at exact equality. `IntrRay3Cone3.find`
  reports all reachable result kinds (empty 1566, point 13, segment 270,
  positive ray 151; `isRayNegative` cannot survive the ray clip).
  `IntrSegment3Cone3.find.throughInterior` sends the segment through a point of
  the solid cone, so it reports a segment on all 2000 records and drives
  clipping blocks 21 to 31, where the uniform `find` mostly reports the empty
  set.
- **Box against cone.** The main case uses a fresh query object per record —
  the defect below needs a reused one.

## Port defects fixed

None. Every value disagreement the oracle found was a deliberate, documented
port fix of an upstream defect (the eight listed below), and no case needed a
tolerance: all 685081 floating-point outputs of the 48 ordinary cases are
bit-identical over the deep run, including the change-of-basis chains
(`IntrTriangle3OrientedBox3`'s `Dot(box.axis[i], v - box.center)`,
`IntrOrientedBox3Cone3`'s transform of the cone into the box frame,
`IntrSegment3Plane3`'s `GetCenteredForm`, and the `origin + t * direction`
reconstructions of every FI query).

No new result-corrupting upstream defect was found in this group.

## Deliberate deviations demonstrated

| case | records disagreeing (golden / deep) | record of the decision |
| --- | --- | --- |
| `IntrTetrahedron3Tetrahedron3.test.edgeCutoffDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrTetrahedron3Tetrahedron3.h` edge-edge cutoff; issue [#307](https://github.com/gradientspaceai/gtengine-js/issues/307) |
| `IntrTetrahedron3Tetrahedron3.test.edgeSeparationDeviation` | 20/20, 2000/2000 | the same entry, the separation test |
| `IntrTriangle3Triangle3.test.moving.parallelDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrTriangle3Triangle3.h` moving-triangle overloads; issue [#334](https://github.com/gradientspaceai/gtengine-js/issues/334) |
| `IntrTriangle3Triangle3.find.moving.parallelDeviation` | 20/20, 2000/2000 | the same |
| `IntrSphere3Cone3.find.vertexDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrSphere3Cone3.h` FI query; issue [#307](https://github.com/gradientspaceai/gtengine-js/issues/307) |
| `IntrRay3Cone3.find.throughVertexDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrLine3Cone3.h` `CaseC2NotZeroDiscrZero` and `CaseC2NotZeroDiscrPos` Block 3; issues [#304](https://github.com/gradientspaceai/gtengine-js/issues/304), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrSegment3Cone3.find.throughVertexDeviation` | 20/20, 2000/2000 | the same |
| `IntrOrientedBox3Cone3.test.staleAdjacencyDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrAlignedBox3Cone3.h` `BoxFullyInConeSlab`; issue [#301](https://github.com/gradientspaceai/gtengine-js/issues/301) |

Each deviation is confined to inputs on which upstream is actually defective,
and each corresponding main case rejects exactly those inputs by an explicit
predicate evaluated on upstream's own control flow rather than by narrowing the
input ranges.

- **`IntrTetrahedron3Tetrahedron3`.** The port corrects two defects of the
  edge-edge phase: upstream compares `|Dot(E0,E1)|` of *unnormalised* edge
  vectors against the cosine cutoff `1 - epsilon`, and its separation test asks
  which side of a plane through one edge endpoint the other tetrahedron lies on
  instead of testing projection-interval disjointness on the cross-product
  axis. The probe `TetraFixed` is the port's algorithm transcribed into the
  case file; the main cases keep the records on which it agrees with upstream
  in both `intersect` and `separating`, and the deviation cases keep the ones
  on which it does not. The two defects are separated: the cutoff case uses
  large tetrahedra, where almost every `|Dot(E0,E1)|` exceeds the cutoff and
  upstream skips the phase entirely (it reports `intersecting` on 1932 of 2000
  records where the port finds an edge-edge separating axis); the separation
  case uses tetrahedra small enough that `TetraCutoffAgrees` confirms the two
  cutoff tests decide all 36 edge pairs identically, so the only remaining
  difference is the separation test (upstream reports `intersecting` on 682 of
  2000 records where the port separates, and a different separating pair on the
  rest).
- **`IntrTriangle3Triangle3`, moving overloads.** Upstream decides "the
  triangles are parallel" with `fabs(Dot(N0,N1)) < 1` on the unnormalised edge
  cross products, so the comparison against 1 scales with the product of the
  triangle areas; the port uses the criterion upstream's own *stationary* query
  uses, `|Cross(N0,N1)|^2 > 0`. The exact predicate
  `(|Dot(N0,N1)| < 1) == (|Cross(N0,N1)|^2 > 0)`, evaluated on the translated
  triangles as the query itself does, says whether the two take the same
  branch; when they do, the whole computation is identical. The main moving
  cases require it (small uniform triangles, where `|Dot(N0,N1)| < 1` and the
  triangles are not parallel, and the exactly coplanar lattice pair, where
  `Cross(N0,N1)` is exactly zero and `|Dot(N0,N1)| >= 1`); the deviation cases
  use triangles large enough that `|Dot(N0,N1)| >= 1` although `Cross(N0,N1)`
  is nonzero, and keep only the records on which upstream's answer differs from
  the port's. Upstream then skips the `N1` axis and all nine
  `E0[i0] x E1[i1]` axes and reports contact for triangles that never touch:
  all 2000 deviation records have upstream `intersect = true`.
- **`IntrSphere3Cone3` FI.** Upstream assigns
  `point = t * (cosAngle*D + tmp*B)`, which omits the cone vertex `V` although
  the reduction it comes from parameterises the ray as `X(t) = V + t*D`; the
  port assigns `V + t*D`. Adding an exactly zero vector is bit-neutral here
  (`t >= 0` always, so no component is a negative zero), so the two expressions
  are bit-identical exactly when `V` is the origin: `find.originVertex` pins
  `V = (0,0,0)` and compares the point, the general `find` compares
  `intersect` on the full generator (the only other field), and
  `find.vertexDeviation` draws a cone whose vertex is not the origin and keeps
  the records that reach the third branch with a reported point — all 2000 of
  them deviate.
- **`IntrRay3Cone3`, `IntrSegment3Cone3`.** Both derive from the line/cone
  FIQuery and clip its t-interval, so they inherit the port's corrections of
  `IntrLine3Cone3` for a line through the cone vertex. `LineConePortDeviates`
  in the case file evaluates, on upstream's own quantities and after the
  direction flip `DoQuery` performs, exactly the condition under which the port
  takes a different path: the hoisted through-vertex test
  `(P - V) + tv*U == 0` succeeds, or `discr < 0` with `c2 > 0`, or `discr > 0`
  with `c2 < 0` in Block 3, or `discr == 0` where the full vector vertex test
  and upstream's U-component test disagree, or `discr == 0` with `c2 > 0` away
  from the vertex. The main cases reject it. Every one of those conditions
  needs the line to contain the cone vertex exactly, or the rounded
  discriminant to take a sign the exact one cannot, so the rejection fires
  only on the (rare) lattice records where the ray origin happens to be
  `V + k*U`; the main generators are not narrowed. The deviation cases
  construct the defect
  exactly: an integer cone vertex `V`, an integer direction `U` and an origin
  at `V + k*U` make `(P - V) + tv*U` exactly the zero vector, so the quadratic
  has a double root. They additionally require `c2 > 0` (the direction is
  interior to the cone's angular region, so the true intersection is the ray of
  nonnegative heights from the vertex), the vertex parameter inside the clipped
  range, a cone minimum height of zero, and — the condition that actually
  separates defective from sound inputs — that upstream's *rounded*
  discriminant is not exactly zero. When the rounding happens to give exactly
  zero, upstream reaches its own vertex branch and computes what the port
  computes; when it does not, its case analysis reports a point, the empty set
  or a segment reaching the cone's maximum height where the answer is a ray.
  With that condition all 2000 records of each case deviate (without it, only
  about 25 % and 7 % did, which is why the condition is in the generator rather
  than in the report as a residue).
- **`IntrOrientedBox3Cone3`.** The query is a thin wrapper around
  `IntrAlignedBox3Cone3`, from which it derives, so it inherits that query's
  stale-adjacency defect: `BoxFullyInConeSlab` copies the twelve box edges over
  `mCandidateEdges` and sets `mNumCandidateEdges = 12` without touching
  `mAdjacencyMatrix`, so the bits an earlier clipping query set for the clipped
  vertices (indices 8..31) are never cleared; a later `ClearCandidates` visits
  only the twelve box edges, `InsertEdge` then treats those candidate edges as
  already present and drops them, and the query reports a false negative. The
  defect needs a query object that is *reused*, so the main case uses a fresh
  object per record and the deviation case runs three calls on one object: a
  clipping configuration (a finite cone whose maximum-height plane passes
  through the box centre, with the axis offset so the quick acceptance does not
  fire), a fully-inside-the-slab configuration (an infinite cone whose vertex
  is below the box), and the first configuration again. The generator keeps the
  records on which the third answer differs from the answer a fresh object
  gives for the same inputs: all 2000 deep-run records report
  `true, true, false` from upstream where the port reports `true, true, true`.

Upstream defects recorded for this group's headers that the port *preserves*
need no deviation case; the oracle compares them bit for bit and they agree:
`IntrTriangle2Triangle2::WhichSide`'s reliance on C++ int truthiness via
`--negative` (#307, a transcription note with no behavioural difference),
`IntrTriangle3Triangle3`'s dead `ContactSide`/`Configuration` machinery (#334,
computed on both sides and read by neither) and its unguarded division by
`normal[lookup[2]]` for a zero-area triangle0 (#334, reached by 2 of the 2000
`find.onPlane` records, where both sides emit the same NaN coordinates), the
`IntrAlignedBox3Cone3`/`IntrOrientedBox3Cone3` header TODO about the retired
infinite-cone representation (#334, documentation only), the undocumented
positive-orientation precondition of the tetrahedron face-normal phase (#465,
which the generator honours), and `IntrLine3Cone3::Result::Convert`, whose
member templates cannot compile upstream (#304) and which the port implements
with the obvious intent; it is not instantiated by any query and has no C++
counterpart to compare against.

## Not covered

Nothing. All 12 headers of group 35 are implemented by the port and every
computational entry point has a case. Three deliberate restrictions are worth
recording:

- The `onFaces`, `throughLatticePoint` and lattice modes use signed
  coordinate axes or signed permutations of them. That is the only way to make
  the world/box and world/plane changes of basis exact, which is what puts the
  clipping and height comparisons on exact equality; the uniform modes of the
  same cases use fully rotated frames.
- `IntrTriangle3Triangle3`'s main moving cases run on small triangles
  (`|Dot(N0,N1)| < 1`) or on exactly coplanar lattice triangles. That is not a
  narrowing chosen for convenience: it is precisely the set on which upstream's
  parallel test agrees with the port's, expressed as an exact predicate and
  enforced by rejection. Large non-parallel triangles are the deviation case.
- `IntrLine3Cone3::Result::Convert` (reached through the ray and segment
  results) has no upstream instantiation, so there is nothing to compare; the
  port's implementation is covered by its own unit tests.

## Upstream bug suspects

None new for this group's headers. One observation worth recording:

- `IntrTriangle3Triangle3`'s moving `FIQuery` recomputes the contact set by
  calling the *stationary* query on the triangles moved to `tFirst`. On a
  first contact the two triangles touch in measure zero, and the stationary
  query's `WhichSide`-style tests treat a touching configuration as a
  separation, so the moving query frequently returns `intersect = true` with an
  empty `intersection`: 1522 of the 2000 `find.moving` records and all 2000
  `find.moving.parallelDeviation` records. Both sides reproduce it bit for
  bit; it is a consequence of upstream's documented convention (contact of
  measure zero is separation) rather than a new defect, but it makes the FI
  moving query's contact set of limited use.
