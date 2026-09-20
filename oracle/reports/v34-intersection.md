# Group 34 (`v34-intersection`) — C++ oracle report

65 cases, 20 records each in `oracle/golden/v34-intersection.txt` (1300 records).
57 of the cases are ordinary comparisons and 8 are deliberate `deviation` cases.
56 of the 57 ordinary cases are declared `exact` (bit-identical to the MSVC
build of upstream GTE on every floating-point output); the remaining one,
`IntrAreaEllipse2Ellipse2.compute`, reaches `atan2`, `atan`, `sin`, `cos` and
`RootsPolynomial::SolveQuartic` and is compared with a tolerance.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise`.

## Coverage

Every computational entry point of every header in the group is covered: each
`TIQuery::operator()` and `FIQuery::operator()`, each protected `DoQuery`
helper and each protected static helper that the port exports as a free
function (called directly on both sides through a small `Expose` subclass in
C++), and the `AreaEllipse2Ellipse2::operator()` of the non-query class.

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `IntrAlignedBox3Cone3.h` | `test`, `test.straddle`, `computeBoxHeightInterval`, `coneAxisIntersectsBox`, `hasPointInsideCone`, `test.staleAdjacencyDeviation` | exact (+1 deviation) | pass |
| `IntrAreaEllipse2Ellipse2.h` | `compute`, `compute.uninitializedDeviation` | tol 1e-8 (libm) (+1 deviation) | pass |
| `IntrConvexMesh3Plane3.h` | `find.tetrahedron`, `find.configurationOnly`, `find.tangent`, `find.box`, `find.coplanarFaceDeviation` | exact (+1 deviation) | pass |
| `IntrLine3Cone3.h` | `find`, `doQuery.fi`, `find.vertexDeviation` | exact (+1 deviation) | pass |
| `IntrLine3Plane3.h` | `test`, `find`, `doQuery.fi`, `find.parallel` | exact | pass |
| `IntrOrientedBox2Circle2.h` | `test`, `find`, `find.aimed` | exact | pass |
| `IntrOrientedBox2Cone2.h` | `test`, `test.aimed` | exact | pass |
| `IntrOrientedBox3Sphere3.h` | `test`, `find`, `find.aimed`, `find.probeDeviation` | exact (+1 deviation) | pass |
| `IntrPlane3Capsule3.h` | `test`, `test.touching` | exact | pass |
| `IntrPlane3Circle3.h` | `test`, `find`, `find.parallel`, `find.tangent` | exact | pass |
| `IntrPlane3Cylinder3.h` | `test`, `test.infinite`, `find`, `find.parallel` | exact | pass |
| `IntrPlane3Ellipsoid3.h` | `test`, `test.tangent` | exact | pass |
| `IntrPlane3OrientedBox3.h` | `test`, `test.touching` | exact | pass |
| `IntrPlane3Sphere3.h` | `test`, `find`, `find.tangent` | exact | pass |
| `IntrRay2Arc2.h` | `test`, `find`, `find.insideDeviation` | exact (+1 deviation) | pass |
| `IntrRay2SegmentMesh2.h` | `find`, `find.throughVertex` | exact | pass |
| `IntrRay3Ellipsoid3.h` | `test`, `find`, `doQuery.fi`, `find.tangent` | exact | pass |
| `IntrSegment2Arc2.h` | `test`, `find`, `find.insideDeviation` | exact (+1 deviation) | pass |
| `IntrSegment2SegmentMesh2.h` | `find`, `find.throughVertex` | exact | pass |
| `IntrSegment3Ellipsoid3.h` | `test`, `test.containedDeviation`, `find`, `doQuery.fi`, `find.tangent` | exact (+1 deviation) | pass |

Only `IntrAreaEllipse2Ellipse2` reaches the C math library on the compared
path. Everywhere else the `sin`/`cos` used to build orthonormal frames, unit
directions and arc endpoints, and the trigonometry of `Cone::SetAngle`, are
applied to *unrecorded* draws: the resulting frame, unit vector, endpoint or
`cosAngle`/`sinAngle`/`tanAngle`/`cosAngleSqr`/`sinAngleSqr`/`invSinAngle` is
recorded as an input and the replay assigns it, so libm never enters the
compared computation.

### Two non-obvious instantiation choices

- **`IntrConvexMesh3Plane3`** `static_assert`s that `Real` is an
  arbitrary-precision type with division, so `FIQuery<double, ...>` does not
  compile. The port instantiates the query for `number`, and its header
  comment says the floating-point behaviour is what a C++ instantiation with
  `double` would do if the assert were removed. The case file therefore
  defines `v34ap::APD`, a transparent wrapper around `double` whose every
  operation is the corresponding IEEE `double` operation, and specializes
  `gte::_is_arbitrary_precision_internal` and
  `gte::_has_division_operator_internal` for it. That instantiates upstream's
  own algorithm on plain doubles without touching the traits for `double`
  itself and without editing the upstream header, so the comparison is
  bit-for-bit against upstream's code.
- **`IntrAreaEllipse2Ellipse2`** cannot be used as its own reference: upstream
  declares `T mZero, mOne, mTwo, mPi, mTwoPi;` with no constructor and never
  assigns them, so every area it computes reads indeterminate values. The main
  case compares the port against `FixedAreaEllipse2Ellipse2` in the case file,
  which is upstream's code verbatim with exactly the port's two documented
  corrections applied (the five constants initialized, and the private axis
  copies normalized). The unmodified upstream class is exercised by the
  deviation case.

### Generators and the branches they reach

Each generator alternates an exactly representable small-lattice mode with a
uniform mode; several add a constructed mode aimed at a specific branch. The
distributions below are from the committed 20-record goldens and the
2000-record deep run.

- **Planes against spheres, ellipsoids, boxes and capsules.** The `touching`,
  `tangent` variants put the plane normal on a coordinate axis and the object
  on the lattice, so the signed distance and the projection radius are exact
  integers and `distance <= radius`, `distance == radius` and
  `distance <= sqrt(N^T M^{-1} N)` are evaluated at exact equality. Both
  answers occur in every boolean case (worst split 3/17).
- **Plane against circle and cylinder.** `IntrPlane3Plane3`'s parallel test is
  the exact `|Dot(N0,N1)| >= 1`, which uniform normals essentially never
  satisfy; the `find.parallel` variants use the same signed coordinate axis for
  both normals, which reaches the coplanar branch (the whole circle,
  `numIntersections` = SIZE_MAX) and the parallel-but-distinct branch.
  `IntrPlane3Cylinder3.find.parallel` reaches all three branches of
  `GetLinesOfIntersection` (two lines, the tangent single line, no
  intersection) and `find` reaches `ELLIPSE` and `CIRCLE` with their trim
  lines.
- **Rays and segments against arcs.** The lattice mode puts the arc endpoints
  on the twelve lattice points of the circle of radius 5, so the side-of-line
  test of `Arc2::Contains` is exact; the ray or segment is aimed at a point of
  that circle. Both `numIntersections` 1 and 2 occur.
- **Segment meshes.** The mesh is an open or closed contiguous polysegment of
  3 to 6 lattice or uniform vertices. The `throughVertex` variants aim an
  integer-direction line at a mesh vertex, which makes the two incident
  segments report exactly equal line parameters — the ties that `std::sort`
  reorders (see the canonicalization note below).
- **Rays and segments against ellipsoids.** The aimed generator sends the
  component at a point at relative radius 0.5 to 1.5 in the ellipsoid's own
  frame, which raises the hit rate from about 5 % to about 45 %. The `tangent`
  variants use an axis-aligned lattice ellipsoid whose extents are powers of
  two, so `M` is exact and the discriminant is exactly zero for the tangency
  offset; the single-root branch occurs on 3 and 2 of the 20 golden records.
- **Line against cone.** The four cone kinds (infinite, infinite truncated,
  finite, frustum) are visited in turn. The lattice mode uses a
  coordinate-axis cone axis and an integer line direction, which reaches
  `c2 == 0` (the line parallel to a cone ruling) and the `SetSegmentClamp` /
  `SetRayClamp` interval clips at exact equality. All five result types
  (`isEmpty`, `isPoint`, `isSegment`, `isRayPositive`, `isRayNegative`) occur.
- **Box against cone.** `test.straddle` centres the box on one of the two
  height planes with a coordinate-axis cone axis and integer heights, so the
  projections are exact integers and the `0` (Z) face configurations of the
  81-entry table occur; it splits 10/10. The `test` generator places the box
  centre on the cone axis with an offset comparable to the box size, which
  splits 5/15 instead of the 18/2 a uniform box gives.
- **Convex mesh against plane.** The mesh is a lattice or uniform tetrahedron
  (4 vertices, 4 faces) or box (8 vertices, 12 faces) with outward
  counterclockwise faces. `find.tangent` uses all-integer plane coefficients
  through one vertex or one edge, which reaches `CFG_POS_SIDE_VERTEX`,
  `CFG_NEG_SIDE_VERTEX`, `CFG_POS_SIDE_EDGE` and `CFG_NEG_SIDE_EDGE`;
  `find.configurationOnly` additionally reaches `CFG_POS_SIDE_POLYGON` and
  `CFG_NEG_SIDE_POLYGON` (it requests no intersection data, so the defective
  `GetIntersectionPolygon` is not run); `find.tetrahedron` and `find.box`
  reach `CFG_SPLIT` on 17 and 18 of 20 records.
- **Ellipse areas.** The two ellipses are drawn with centres within a short
  offset of each other, which reaches `ELLIPSES_ARE_SEPARATED`,
  `E0_CONTAINS_E1`, `E1_CONTAINS_E0`, `ONE_CHORD_REGION` and
  `FOUR_CHORD_REGION`.

### Outputs canonicalized on both sides

- `IntrRay2SegmentMesh2` and `IntrSegment2SegmentMesh2` delegate to
  `FIQuery<Line2,SegmentMesh2>`, which orders its output with `std::sort` on
  the line parameter alone. `std::sort` is not stable and
  `Array.prototype.sort` is, so records with equal line parameters — which the
  `throughVertex` generators produce on purpose — may come out in either
  order. Both sides re-sort the reported intersections by the total key
  (parameter, indexPair[0], indexPair[1], meshSegmentParameter, point) before
  emitting them.
- `FIQuery<Plane3,Circle3>::Result::numIntersections` uses `SIZE_MAX` for "the
  whole circle" and `FIQuery<Ellipse2,Ellipse2>::Result::numPoints` uses it for
  "the ellipses are the same". Neither is representable as a double (the port
  uses `Number.MAX_SAFE_INTEGER`), so both sides emit a small code instead of
  the raw value.

## Port defects fixed

Four, all found by this oracle. Each has a regression test that runs the port
beside the two candidate computations, requires bit-identity with upstream's
and requires the two candidates to differ on at least one input.

| file | what differed | root cause | size |
| --- | --- | --- | --- |
| `src/Hyperellipsoid.ts` | `fromCoefficientsABC`'s centre, axes and extents | `Inverse(A, &invertible)` resolves to the *closed-form* `Inverse` of `Matrix2x2.h`/`Matrix3x3.h`, which is more specialized than the Gaussian-elimination template of `Matrix.h`; the port called the generic one | up to 1.2e-14 scaled, on 4 of 20 `IntrPlane3Cylinder3.find` records |
| `src/IntrOrientedBox2Circle2.ts` | `contactPoint` | upstream accumulates `box.center + s0*P0*axis0 + s1*P1*axis1` left to right; the port added the two basis terms first | 1 ulp, on 2 of 20 records |
| `src/IntrOrientedBox3Sphere3.ts` | `contactPoint` | the same, with three basis terms | 1 ulp, on 3 of 20 records |
| `src/IntrPlane3Cylinder3.ts` | `ellipse.center` | the same, with two basis terms | 1 ulp, on 4 of 20 records |

The `Hyperellipsoid` defect is in a file of another verify group but it is
reachable only through this group's `IntrPlane3Cylinder3`, so it is fixed here.
It is worth a note for the owner of `Hyperellipsoid.h`: upstream's
`Hyperellipsoid.h` includes only `Matrix.h`, so which `Inverse` overload a C++
translation unit actually picks depends on whether `Matrix2x2.h` or
`Matrix3x3.h` is also included at the point of instantiation. The port has no
such ambiguity — every module is always available — and the specialized
overloads exist precisely to be used, so the port dispatches on the dimension.

## Deliberate deviations demonstrated

| case | records disagreeing (golden / deep) | record of the decision |
| --- | --- | --- |
| `IntrOrientedBox3Sphere3.find.probeDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrAlignedBox3Sphere3.h` `DoQueryRayRoundedFace`/`DoQuery`; issues [#458](https://github.com/gradientspaceai/gtengine-js/issues/458), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrRay2Arc2.find.insideDeviation` | see below | `docs/UPSTREAM-FINDINGS.md` `IntrRay2Arc2.h` `operator()`; issue [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrSegment2Arc2.find.insideDeviation` | see below | `docs/UPSTREAM-FINDINGS.md` `IntrSegment2Arc2.h` `operator()`; issue [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrSegment3Ellipsoid3.test.containedDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrSegment3Ellipsoid3.h` TIQuery; issue [#304](https://github.com/gradientspaceai/gtengine-js/issues/304) |
| `IntrLine3Cone3.find.vertexDeviation` | see below | `docs/UPSTREAM-FINDINGS.md` `IntrLine3Cone3.h`; issues [#304](https://github.com/gradientspaceai/gtengine-js/issues/304), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrAlignedBox3Cone3.test.staleAdjacencyDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrAlignedBox3Cone3.h` `BoxFullyInConeSlab`; issue [#301](https://github.com/gradientspaceai/gtengine-js/issues/301) |
| `IntrAreaEllipse2Ellipse2.compute.uninitializedDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrAreaEllipse2Ellipse2.h`; issue [#301](https://github.com/gradientspaceai/gtengine-js/issues/301) |
| `IntrConvexMesh3Plane3.find.coplanarFaceDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrConvexMesh3Plane3.h` `GetIntersectionPolygon`; issue [#301](https://github.com/gradientspaceai/gtengine-js/issues/301) |

Each deviation is confined to inputs on which upstream is actually defective,
and the corresponding main case rejects those inputs by an explicit predicate
evaluated on upstream's own control flow rather than by narrowing the input
ranges:

- **`IntrOrientedBox3Sphere3`.** The find-intersection query is a thin wrapper
  around the protected `DoQuery` of `FIQuery<AlignedBox3,Sphere3>`, so it
  inherits that query's first-wins probe selection and the port's correction of
  it. The case file carries the corrected `DoQuery` that
  `oracle/cpp/cases/v32-intersection.cpp` uses, applied to the transformed
  inputs the oriented wrapper builds; the two main cases keep only the inputs
  where the corrected answer equals upstream's bit for bit.
- **`IntrRay2Arc2` / `IntrSegment2Arc2`.** Upstream reuses the ray/circle and
  segment/circle queries, which clip against the *solid disk*: when a root of
  the line/circle quadratic falls outside the ray or segment parameter range,
  the clip replaces it by the range endpoint, which is not on the circle, and
  `Arc2::Contains`'s one-argument form accepts it anyway. The port intersects
  the circular curve and filters by the parameter range instead. The probe is
  the port's own computation, written out in the case file, compared field by
  field with `memcmp` against upstream's answer; the main cases reject the rays
  and segments where the two differ, and the deviation cases place the ray
  origin (respectively one segment endpoint) strictly inside the disk.
- **`IntrSegment3Ellipsoid3`.** Upstream's `qm * qp <= 0` test misses the case
  where both endpoint values of `Q` are negative, that is where the whole
  segment lies strictly inside the solid ellipsoid; it then falls through to
  `qm > 0 && |a1| < a2*e` and reports no intersection, contradicting its own FI
  query. The predicate `Q(-e) < 0 && Q(e) < 0`, evaluated on upstream's own
  expressions after upstream's `discr < 0` early exit and on the centred form
  the query itself uses, is exactly the set on which the port deviates.
- **`IntrLine3Cone3`.** A line through the cone vertex makes
  `Q(t) = c2*(t-tv)^2` a double root, so `discr = c1*c1 - c0*c2` is zero in
  exact arithmetic and a cancelling difference with no significant digits in
  floating point; upstream's case analysis is then decided by the rounded sign
  and reports a point, the empty set, or a segment reaching the cone's maximum
  height where the answer is a ray or the vertex alone. The probe follows
  upstream's control flow on the direction `DoQuery` actually uses and reports
  the four places where the port's correction can change the answer: the line
  containing the vertex exactly, `discr < 0` with `c2 > 0`, `discr > 0` with
  `h[0] < 0 < h[1]` and `c2 < 0`, and `discr == 0` with `c2 > 0` away from the
  vertex. The first of those is a *superset* of the deviating inputs (upstream
  sometimes agrees there), which narrows the main generator but cannot hide a
  disagreement. The deviation case places the line origin at `V - t*U` with an
  integer direction `U`, so `tv` is recovered exactly and the vertex is on the
  line in exact arithmetic.
- **`IntrAlignedBox3Cone3`.** `BoxFullyInConeSlab` copies the twelve box edges
  into `mCandidateEdges` and sets `mNumCandidateEdges = 12` without clearing
  `mAdjacencyMatrix`, so the bits an earlier clipping query set for the clipped
  vertices (indices >= 8) survive; the next `ClearCandidates` visits only the
  twelve box edges, and the next clipping query then silently drops those
  candidate edges through `InsertEdge` and reports a false negative. The
  deviation case runs three queries on ONE query object — a clipping
  configuration, a fully-in-slab configuration and the first configuration
  again — and emits all three answers; the first two agree and the third does
  not. The generator rejects the triples on which the reused object and a fresh
  object agree. Every other case constructs a fresh query object per record, as
  the port does, so none of them is affected.
- **`IntrAreaEllipse2Ellipse2`.** The five never-assigned members are read on
  every path. Value-initializing the query object (`AreaEllipse2Ellipse2<double>
  query{}`; the class has no user-provided constructor, so `{}` zero-initializes
  it) makes them deterministically zero, which is the only way to get a
  reproducible golden out of code that otherwise reads indeterminate stack
  values; the golden file is byte-identical across repeated generations. With
  `mTwo = 0` the half-area `AB/mTwo` is infinite, with `mPi = 0` every
  containment area is zero and with `mTwoPi = 0` the `atan2` wrap never
  happens, so 16 of the 20 golden records report `Infinity` or `NaN`; the four
  finite ones are the containment and separation branches, where the area is
  `mPi * AB = 0`, and their configuration is decided by `qform > mOne`, that is
  `qform > 0`, so it too is wrong.
- **`IntrConvexMesh3Plane3`.** `GetIntersectionPolygon` builds
  `polygonIndices` as a predecessor map indexed by vertex index and then reads
  it back in index order instead of traversing it as a cycle. For the
  four-vertex box face the boundary cycle is `0 -> 2 -> 3 -> 1` after
  `RemoveDuplicateAndUnusedVertices` repacks the vertices, and upstream's
  readback is `(1,3,0,2)`, which is not a rotation of it. The effect is not
  confined to faces with four or more vertices: for a three-vertex face the
  readback is the *reversed* cycle, so the main `find.tangent` case
  deliberately omits the face-containing plane and the `find.tetrahedron` and
  `find.box` generators reject the planes that contain a face (upstream's own
  `numZero >= 3` bit of the configuration, obtained from a
  `REQ_CONFIGURATION_ONLY` probe).

All other upstream defects recorded for this group's headers are *preserved* by
the port and therefore need no deviation case; the oracle compares them bit for
bit and they agree: `IntrPlane3Cylinder3`'s CIRCLE-versus-ELLIPSE label chosen
by the exact extent comparison `ellipse2.extent[0] != ellipse2.extent[1]` and
its discarded `FromCoefficients` status (#465, #304), `IntrPlane3Plane3`'s exact
`|Dot(N0,N1)| >= 1` parallel test seen through `IntrPlane3Circle3` and
`IntrPlane3Cylinder3` (#465), `IntrSegment2Arc2`'s duplicated
`parameter[0] = 0` assignment (#304, a no-op), `IntrLine2SegmentMesh2`'s
`+-max()` sentinel copied into `lineParameter` and the resulting
`IntrRay2SegmentMesh2` filtering (#461, reached only for a mesh segment
collinear with the query line), and the `IntrAlignedBox3Cone3` /
`IntrOrientedBox3Cone3` TODO about the retired infinite-cone representation
(#334) — the port substitutes `Number.MAX_VALUE` for
`std::numeric_limits<double>::max()` exactly as upstream does, and the
`test`, `test.straddle` and `coneAxisIntersectsBox` cases all include infinite
cones.

## Not covered

All 20 headers of group 34 are implemented by the port and every computational
entry point has a case. Two deliberate restrictions are worth recording:

- `FIQuery<Real,Line3,Cone3>::Result::Convert`, whose port is
  `intrLine3Cone3Convert` / `intrLine3Cone3ConvertPoint`, has no case: the
  upstream member templates cannot be instantiated at all (`QFNumber` has no
  conversion operator to `Real`), so there is nothing to compare against
  (`docs/UPSTREAM-FINDINGS.md` `IntrLine3Cone3.h`, issue
  [#304](https://github.com/gradientspaceai/gtengine-js/issues/304)). The port
  implements the obvious intent, `x[0] + x[1]*sqrt(d)`.
- `SetRayNegative` of the same header is dead code upstream (`DoQuery` converts
  `isRayPositive` to `isRayNegative` by assigning the type), and the port omits
  it. The `isRayNegative` result type itself is reached by the `find` cases.

## Upstream bug suspects

None new for this group's headers.

One observation about a *shared* file, for the group that owns it, is recorded
above under "Port defects fixed": `Hyperellipsoid::FromCoefficients` calls
`Inverse(A, &invertible)` while `Hyperellipsoid.h` includes only `Matrix.h`, so
which overload a C++ translation unit resolves to — the closed-form
`Matrix2x2`/`Matrix3x3` one or the Gaussian-elimination template — depends on
whether `Matrix2x2.h` or `Matrix3x3.h` happens to be included at the point of
instantiation. Two translation units of the same program can therefore compute
different ellipse parameters for the same input. It is not result-corrupting
(the two inverses agree to rounding), but it is an ODR-adjacent hazard worth an
`#include` in `Hyperellipsoid.h`.
