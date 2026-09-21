# Group 14 (`v14-containment`) — C++ oracle report

70 cases, 20 records each in `oracle/golden/v14-containment.txt` (1400 records).
64 of the cases are ordinary comparisons and 6 are deliberate `deviation`
cases. **All 64 ordinary cases are declared `exact`**: every floating-point
output is bit-identical to the MSVC build of upstream GTE on every record.
There is no tolerance case in this group — nothing on a compared path reaches
the C math library.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise`.

Deep run: `npm run oracle:deep -- 2000 v14-containment` passes (70 cases,
140000 records, 4.2 s wall including the replay).

## Coverage

Every computational entry point of every header in the group is covered: each
`GetContainer`, each `InContainer` overload (distinguished by the contained
type), each `MergeContainers`, the four `Circumscribe` / `Inscribe` functions,
all four `PointInPolygon2` queries, the two `Cont*MinCR` functors, and all six
supported `(FaceType, method)` pairs of `PointInPolyhedron3::Contains` plus its
silent-`false` unsupported pair.

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `ContAlignedBox.h` | `getContainer.2d`, `getContainer.3d`, `inContainer.2d`, `inContainer.3d`, `mergeContainers.3d`, `mergeContainers.signedZero` | exact | pass |
| `ContAlignedBox2Arc2.h` | `getContainer` | exact | pass |
| `ContCapsule3.h` | `getContainer`, `inContainer.point`, `inContainer.sphere`, `inContainer.capsule`, `mergeContainers` | exact | pass |
| `ContCircle2.h` | `getContainer`, `inContainer`, `mergeContainers` | exact | pass |
| `ContCone.h` | `inContainer.2d`, `inContainer.3d` | exact | pass |
| `ContCylinder3.h` | `getContainer`, `inContainer` | exact | pass |
| `ContEllipse2.h` | `getContainer`, `inContainer`, `inContainer.leftHanded`, `mergeContainers` | exact | pass |
| `ContEllipse2MinCR.h` | `compute`, `compute.verticalLineDeviation` | exact (+1 deviation) | pass |
| `ContEllipsoid3.h` | `getContainer`, `inContainer`, `inContainer.leftHanded`, `mergeContainers`, `mergeContainers.leftHanded` | exact | pass |
| `ContEllipsoid3MinCR.h` | `compute`, `compute.assertDeviation` | exact (+1 deviation) | pass |
| `ContLozenge3.h` | `getContainer`, `getContainer.sphereBranch`, `inContainer`, `getContainer.cornerDeviation` | exact (+1 deviation) | pass |
| `ContOrientedBox2.h` | `getContainer`, `inContainer`, `inContainer.leftHanded`, `mergeContainers` | exact | pass |
| `ContOrientedBox3.h` | `getContainer`, `inContainer`, `inContainer.leftHanded`, `mergeContainers`, `mergeContainers.leftHanded` | exact | pass |
| `ContPointInPolygon2.h` | `contains.convex`, `contains.nonconvex`, `containsConvexOrderN`, `containsConvexOrderLogN`, `containsQuadrilateral` | exact | pass |
| `ContPointInPolyhedron3.h` | `contains.triangle`, `contains.convex0`, `contains.convex1`, `contains.convex2`, `contains.simple0`, `contains.simple1`, `contains.unsupported`, three `*QuadDeviation` cases | exact (+3 deviations) | pass |
| `ContScribeCircle2.h` | `circumscribe`, `inscribe` | exact | pass |
| `ContScribeCircle3Sphere3.h` | `circumscribeCircle3`, `circumscribeSphere3`, `inscribeCircle3`, `inscribeSphere3` | exact | pass |
| `ContSphere3.h` | `getContainer`, `inContainer`, `mergeContainers` | exact | pass |
| `ContTetrahedron3.h` | `inContainer` | exact | pass |

### Why every case is exact

The whole group is `+ - * / sqrt fabs` and comparisons, including the
`ApprGaussian2` / `ApprGaussian3` / `ApprOrthogonalLine3` covariance fits and
the `SymmetricEigensolver2x2` / `SymmetricEigensolver3x3` they run (the
iterative solver, not the `acos`/`cos`-based `NISymmetricEigensolver3x3`), the
`Matrix2x2` / `Matrix3x3` closed-form inverses that `LinearSystem::Solve` uses
for the scribe functions, `Projection.h`'s ellipse and ellipsoid projections,
and the quaternion round trip in the oriented-box and ellipsoid merges. The
only libm in the group is `Cone::SetAngle`; the six values it derives from the
angle are computed on an *unrecorded* draw and recorded as inputs, which the
replay assigns, so `cos`/`sin`/`tan` never enter a compared computation. The
same applies to the `cos`/`sin` used to build unit directions, orthonormal
frames and non-lattice arc endpoints.

### Generators

Every generator alternates a small exact lattice mode with a uniform mode, and
adds constructed modes aimed at the non-generic branches:

- point clouds in five kinds — lattice, uniform, all points identical (rank 0),
  collinear on a lattice line (rank 1) and coplanar on a lattice plane
  (rank 2), which is what drives `ApprGaussian` to exactly zero eigenvalues and
  reaches upstream's infinite-extent behaviour in `ContEllipse2` /
  `ContEllipsoid3`;
- boundary points built by exact construction so that `<=` / `>` is evaluated
  at equality: a box face or corner, `center + r * e_k` on a sphere, the
  capsule end cap, the flat cap and lateral surface of a cylinder, the lozenge
  corner offset by the radius along the rectangle normal, the cone height
  interval endpoints, and the box/ellipsoid face centres and corners;
- `ContPointInPolygon2` uses strictly convex counterclockwise lattice polygons
  with 4 to 10 vertices, built from edge vectors of strictly increasing angle
  followed by their negations (so the polygon closes exactly), plus a
  non-convex L-hexagon for `Contains`; the query point is a vertex, an edge
  midpoint, an edge quarter point, a point on the horizontal ray through a
  vertex (where the `<` / `<=` asymmetry of `Contains` is decided), a lattice
  point or a uniform point;
- `ContAlignedBox2Arc2` places the arc endpoints on the twelve exact lattice
  points of a radius-5 circle, so the four axis points `C±(r,0)`, `C±(0,r)`
  are themselves possible endpoints and `Arc2::Contains` is evaluated at
  `DotPerp == 0`;
- the merge cases mix identical, nested and disjoint inputs;
  `ContAlignedBox.mergeContainers.signedZero` feeds `±0` to `std::min` /
  `std::max` so that the port's `stdMin` / `stdMax` is pinned (`Math.min` /
  `Math.max` order `-0` below `+0`);
- all frames are right-handed (2D `axis1 = -Perp(axis0)`, because GTE's
  `Perp(x,y) = (y,-x)` is the clockwise one; 3D `axis2 = Cross(axis0, axis1)`,
  `det = +1`). Five extra `.leftHanded` cases feed `det = -1` frames to the
  oriented-box, ellipse and ellipsoid queries on purpose.

Branch histogram over the deep run: every boolean case splits, the worst being
`ContPointInPolyhedron3.contains.simple1QuadDeviation` at 24 % / 76 % and
`ContPointInPolygon2.contains.nonconvex` at 20 % / 80 %; most are between
35 % and 65 %. `ContPointInPolyhedron3.contains.unsupported` is 100 % `false`
by construction (that is the point of the case). The `getContainer` and
`mergeContainers` cases return `true` on every record because upstream's
`bool` is vestigial there.

## Port defects fixed

### 1. `src/ContEllipsoid3MinCR.ts` — the plane jitter used the wrong `std::generate_canonical`

`MaxProduct` jitters every constraint coefficient by
`1e-12 * std::uniform_real_distribution<double>(0,1)(std::mt19937())`. The port
reproduced the combination the C++ standard describes and libstdc++/libc++
compute, `(g0 + g1 * 2^32) / 2^64` evaluated in floating point. The MSVC STL,
the reference build, takes the power-of-two shortcut in
`std::generate_canonical<double, 53>` instead: it forms the 53-bit integer
`S = (g0 >> 11) + (g1 << 21)` and returns `S * 2^-53`, **truncating** the low
11 bits of the first draw where the portable formula **rounds** them. The first
value of the stream is `0x3fc1574f7b6848dc` under MSVC and
`0x3fc1574f7b6848de` under the portable formula, 2 ulps apart.

The difference is invisible while the jitter is added to a coefficient of
order 1, but it *is* the whole coefficient when an input point lies on a
coordinate plane of the ellipsoid frame (`A[i][k] = 0` before the jitter), and
`D[2] = 1/zmax` then differs in its last bits. Observed: 3 of 20 golden records
disagreed, by 1 to 2 ulps (largest scaled error 3.4e-16, on
`D[2] ≈ 3.57e11`). The port now uses the MSVC construction and all 2000 deep
records agree bit for bit. The `ContEllipsoid3MinCR.compute` generator flattens
every other cloud into the plane of the first two axes precisely so that this
path is exercised.

Regression test: `test/ContEllipsoid3MinCR.test.ts`, "the plane jitter matches
the reference std::mt19937 stream" — it pins the first canonical value to the
hex MSVC prints, shows that the two constructions differ (so the test is not
vacuous), and checks the full `D` for a one-point cloud at the centre, where
the result is a pure function of the first three jitter values.

### 2. `src/ContLozenge3.ts` — change-of-basis accumulation order

Upstream writes the rectangle centre as
`box.center + a * box.axis[2] + b * box.axis[1]`, which C++ accumulates left to
right, `((C + a*A2) + b*A1)`. The port grouped the two basis terms first,
`C + (a*A2 + b*A1)` — the eighth instance of this defect class in the wave. The
two differ in the last bits on essentially every input.

This one is *not* observable through the oracle, and the report says so
plainly: the port deliberately replaces upstream's corner offsets by the
interval midpoints (issue #174, below), so the only branch in which upstream
evaluates the same coefficients as the port is its "container is a sphere"
branch — and that branch is only reachable with `a = b = 0`. (The extreme-`w`
points have `radical == 0`, so `aMin >= aMax` forces them to share a single `u`,
and because `box.center` is the mean of the cloud that shared value is 0; a
collinear cloud therefore lands in the *capsule* branch, not the sphere one.)
`ContLozenge3.getContainer.sphereBranch` pins the branch and the centre with
zero offsets; the accumulation order is pinned by the regression test
`test/ContLozenge3.test.ts`, "the rectangle centre accumulates the basis terms
left to right", which runs the port beside both groupings over 400 clouds,
requires bit-identity with upstream's grouping and requires the two groupings
to differ somewhere.

## Deliberate deviations demonstrated

| case | issue | deviating records (of 2000, deep run) |
| --- | --- | --- |
| `ContLozenge3.getContainer.cornerDeviation` | [#174](https://github.com/gradientspaceai/gtengine-js/issues/174) — upstream centres the fitted rectangle on a *corner* of the parameter interval while `Rectangle<3>` is a centred primitive, so the lozenge does not contain its own points | 2000 |
| `ContEllipse2MinCR.compute.verticalLineDeviation` | [#234](https://github.com/gradientspaceai/gtengine-js/issues/234) — `D[1] = (1 - a0*x0)/A[iYMin][1]` divides by zero whenever the walk steps onto a vertical constraint line, returning NaN or ±infinity | 1999 |
| `ContEllipsoid3MinCR.compute.assertDeviation` | [#409](https://github.com/gradientspaceai/gtengine-js/issues/409) — `LogAssert(numer >= 0)` contradicts its own adjacent comment and fires on an already-active plane with a slack of about -1e-16; upstream throws where the port clamps and returns | 2000 (upstream throws on every record) |
| `ContPointInPolyhedron3.contains.convex0QuadDeviation` | [#343](https://github.com/gradientspaceai/gtengine-js/issues/343) — `Face::indices` is `std::array<int32_t,3>` while `ContainsC0` reads `indices.size()` as the vertex count, so a quad face is truncated to its first triangle | 956 |
| `ContPointInPolyhedron3.contains.convex12QuadDeviation` | #343, the same truncation in `SharedContains` (methods 1 and 2 alternate) | 661 |
| `ContPointInPolyhedron3.contains.simple1QuadDeviation` | #343, the same truncation for SIMPLE method 1 | 952 |

Each deviation case is confined to the inputs on which upstream is defective,
and the corresponding main case is broad:

- `ContLozenge3.getContainer` runs the full five-kind cloud generator and
  compares every field of the lozenge except the rectangle centre (the split by
  output field that ORACLE.md prescribes for a re-derivation);
  `ContLozenge3.getContainer.sphereBranch` compares the centre too on the one
  branch where upstream's formula and the port's coincide.
- `ContEllipse2MinCR.compute` rejects, with a capped 24-attempt loop that
  redraws every quantity, any draw whose upstream result has a non-finite
  component; the deviation case is the complement, driven by forcing one input
  point exactly onto the first ellipse axis so that its `A[1]` is exactly zero.
- `ContEllipsoid3MinCR.compute` classifies each draw with a replica of
  upstream's facet/edge walk (in the case file, next to the case) before
  running the real query, and accepts only draws on which upstream terminates
  and never sees a negative numerator; the deviation case accepts only draws on
  which upstream's assert fires. Over 50000 random integer clouds the
  classifier agreed with the real query on every draw.
- The polyhedron main cases use faces with exactly three vertices, which
  upstream's fixed-size `indices` array represents faithfully; the deviation
  cases use the quadrilateral faces of a box. `contains.simple0` uses the same
  quad box faces as a *main* case, because method 0 for simple faces reads
  `Face::triangles`, a `std::vector`, and is therefore not truncated.

## Independent-reference checks

Run once over the 2000-record deep run (the checking script is not committed;
it is a straightforward replay of the recorded inputs through the port plus a
geometric predicate). Because the main cases are bit-identical, checking the
port's outputs checks upstream's.

| check | result |
| --- | --- |
| every input point is inside its `getContainer` result | holds for aligned box (exactly), circle (0), sphere (0), oriented box 2D/3D (≤ 1.4e-14), capsule (≤ 1.9e-14), cylinder (≤ 2.8e-14), **lozenge (≤ 2.7e-14)** and ellipse/ellipsoid (≤ 2.3e-16 on `Q`) |
| a merged container contains both inputs | holds for the aligned box (exactly), circle and sphere (≤ 2.7e-15), oriented box 2D/3D over all vertices (≤ 3.6e-15) and capsule over both end spheres (≤ 5.4e-15) |
| circumscribed circles/spheres pass through their vertices | relative residual ≤ 3.1e-14 (circle2 1.9e-15, circle3 3.1e-14, sphere3 6.4e-15) |
| inscribed circles/spheres are tangent to the sides | relative residual ≤ 7.0e-14 (circle2 1.9e-15, sphere3 7.0e-14) |
| point-in-polygon matches a crossing-number reference off the boundary | 3927 checks, 0 disagreements, across `Contains` (convex and non-convex), `ContainsConvexOrderN` and `ContainsConvexOrderLogN` |
| point-in-polyhedron matches an exact half-space test on the box | 390 checks with general-position rays, 0 disagreements |
| `Cont*MinCR` results are feasible (`sum A_i D_i <= 1` for every point) | worst excess 6.3e-12 (2D), 0 (3D) |

Two notes from those checks, neither a defect:

- **`ContEllipse2.mergeContainers` really does not contain its inputs.** 3812 of
  4000 input ellipses in the deep run have a boundary point outside the merged
  ellipse. This is exactly the behaviour documented under issue #292
  ("the merged ellipse does not contain its inputs"), preserved in the port,
  and it is the reason the merge cases here do not assert containment.
- **`PointInPolyhedron3` with an axis-aligned ray on an axis-aligned box.** Two
  of 865 records disagreed with the half-space test until records whose ray
  direction is a coordinate axis were excluded; with general-position rays the
  agreement is 390/390. Upstream documents this ("If the ray intersects an edge
  or a vertex, then the counting must be handled differently … You should
  choose these rays randomly"), so this is the algorithm's documented
  limitation surfacing through a deliberately degenerate generator, not a bug.
  The generator keeps those records: both sides agree on them bit for bit, and
  they are the interesting inputs for the parity counting.

## Not covered

- **`ContEllipsoid3MinCR.h`, the non-terminating facet/edge walk**
  (issue #409, item 2). The oracle cannot record it: on such a cloud upstream
  recurses until the stack is exhausted, which is `STATUS_STACK_OVERFLOW` and
  not a catchable C++ exception, so the driver would die rather than write a
  throw record. Measured frequency on random integer clouds of 4 to 8 points:
  0.2 % (108 of 50000). Those draws are rejected by the classifier described
  above, and the port's bounded walk is covered by the existing property tests
  in `test/ContEllipsoid3MinCR.test.ts`.
- **The `std::vector<Vector>` convenience overloads** of `GetContainer` in
  `ContCircle2.h`, `ContSphere3.h`, `ContOrientedBox2.h` and
  `ContOrientedBox3.h`. They forward to the `(count, pointer)` overload with no
  arithmetic of their own, and the port has a single array-taking function, so
  there is nothing separate to compare.
- **`PointInPolyhedron3::Contains` for `(SIMPLE, method >= 2)`** is covered as
  the silent-`false` case `contains.unsupported`; there is no computation
  behind it.

## Upstream bug suspects

No new ones. Everything this group ran into is already in
`docs/UPSTREAM-FINDINGS.md`: the `Cont*` empty-input family (#106), the
`ContLozenge3` corner centring (#174), the `ContEllipse2MinCR` division by
zero (#234), the `ContEllipse2`/`ContEllipsoid3`/`ContOrientedBox3` merge and
zero-eigenvalue behaviour (#292), the `ContPointInPolyhedron3` face truncation
(#343) and the `ContEllipsoid3MinCR` assert and non-termination (#409). The
left-handed-frame cases confirmed what #292 predicts: `InContainer` for
oriented boxes, ellipses and ellipsoids is handedness-independent (it projects
onto each axis separately) and both sides agree bit for bit on `det = -1`
frames, while `MergeContainers` for `OrientedBox3` and `Ellipsoid3` runs a
reflection through `Rotation<3,Real>`'s matrix-to-quaternion conversion and
produces a meaningless frame — identically on both sides, so it is compared
bit for bit and left alone.

The MSVC-versus-portable `std::generate_canonical` difference (port defect 1
above) is a standard-library portability hazard rather than a GTE bug: any
`ContEllipsoid3MinCR` result depends on which standard library built it, in
the last bits generally and materially when an input point lies on a
coordinate plane of the ellipsoid frame.
