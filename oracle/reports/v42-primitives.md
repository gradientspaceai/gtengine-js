# Verify group 42 (primitives) against the MSVC build of upstream GTE

Family `v42-primitives`, 48 cases, 20 golden records each, deep run
`npm run oracle:deep -- 2000 v42-primitives` (96000 records, 2986 of them
throw records): **every case passes**. Wall time of the deep run, generation
and replay together: about 9 seconds on the first run, 2.5 seconds warm.

This group is the library's geometric value types. Apart from three entry
points that call the C math library (`Sector2::SetAngle`, `Torus3::Evaluate`
and `Torus3::GetParameters`), every compared path is `+ - * / sqrt fabs` and
comparisons, so **45 of the 48 cases are declared `{ exact: true }` and every
real output of those cases is bit-identical to the MSVC value**. The three
libm cases carry the default scaled tolerance of `1e-12`; their measured
worst scaled error over the deep run is `8.9e-16`.

**No port defect was found.** Every association order, accumulation seed and
division form in `src/` already matches upstream bit for bit, so this branch
contains no `src/` change.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `AlignedBox.h` | `comparisons.3d`, `centeredFormAndVertices.3d`, `centeredFormAndVertices.2d` | exact | pass, 2000 records each |
| `CanonicalBox.h` | `comparisons.2d`, `getVertices.3d`, `getVertices.2d` | exact | pass |
| `Halfspace.h` | `comparisons.3d` (+ defaults) | exact | pass |
| `Hypersphere.h` | `comparisons.3d` (+ defaults) | exact | pass |
| `Line.h` | `comparisons.3d` (+ defaults) | exact | pass |
| `OrientedBox.h` | `comparisons.3d`, `getVertices.3d`, `getVertices.2d` | exact | pass |
| `Ray.h` | `comparisons.2d` (+ defaults) | exact | pass |
| `Rectangle.h` | `comparisons.3d`, `getVertices.3d`, `getVertices.2d` | exact | pass |
| `Segment.h` | `comparisons.3d`, `getCenteredForm.3d`, `getCenteredForm.2d`, `setCenteredForm.3d` | exact | pass |
| `SegmentMesh.h` | `constructors.2d`, `constructors.3d` | exact | pass, 447 / 449 throw records |
| `Triangle.h` | `comparisons.2d` (both constructors) | exact | pass |
| `Arc2.h` | `comparisons`, `contains` (both overloads) | exact | pass |
| `Capsule.h` | `comparisons.3d` (+ defaults) | exact | pass |
| `Circle3.h` | `comparisons` (+ defaults) | exact | pass |
| `Cylinder3.h` | `comparisons`, `finiteInfinite` | exact | pass |
| `Ellipse3.h` | `comparisons` (+ defaults) | exact | pass |
| `Frustum3.h` | `comparisons`, `updateAndComputeVertices` | exact | pass |
| `Lozenge3.h` | `comparisons` (+ defaults) | exact | pass |
| `Parallelepiped3.h` | `comparisons`, `getVertices`, `rightHandedAssert` | exact | pass, 1034 throw records |
| `Parallelogram2.h` | `comparisons`, `getVertices`, `rightHandedAssert` | exact | pass, 1056 throw records |
| `Polyhedron3.h` | `queries`, `invalidInput` | exact | pass |
| `Sector2.h` | `comparisons`, `contains` | exact | pass |
| `Sector2.h` | `setAngle` | tolerance 1e-12 (`std::cos`, `std::sin`) | pass, worst 1.1e-16 |
| `Torus3.h` | `comparisons` | exact | pass |
| `Torus3.h` | `evaluate` | tolerance 1e-12 (`std::cos`, `std::sin`) | pass, worst 8.9e-16 |
| `Torus3.h` | `getParameters` | tolerance 1e-12 (`std::atan2`) | pass, worst 2.2e-16 |
| all 23 headers | `Primitives.defaultConstructors` | exact | pass |

Every computational public entry point of the group is covered. The default
constructor of every type is in `Primitives.defaultConstructors`, which runs
the dimension-templated types at N = 2, 3 and 4 and the fixed-dimension types
once, and also evaluates `Frustum3::Update` / `ComputeVertices`,
`Parallelepiped3::GetVertices` and `Parallelogram2::GetVertices` on the
default state. The value constructors are exercised by the cases that use
them; `Segment`'s four constructors and `Triangle`'s two are each covered
(`Segment(std::array)` in `setCenteredForm.3d`, `Triangle(std::array)` in
`comparisons.2d`).

### Comparison operators

All 21 types that define comparison operators get a case that emits twelve
booleans: `A == B`, `A != B`, `A < B`, `A <= B`, `A > B`, `A >= B`, the two
reversed relations `B < A` and `B <= A`, and the four self-relations
`A == A`, `A < A`, `A <= A`, `A >= A`. Discrete outputs always compare
exactly, so none of these records is vacuous even when every member is NaN.

The generator draws the K members of A, then the K members of B with a random
prefix copied exactly from A (or replaced by the other signed zero). Without
those prefix ties the first member decides every comparison and the member
*order* of `operator<` is never observed. Over the deep run of, for example,
`OrientedBox.comparisons.3d`: 2000 records, 500 with a NaN member, 869 with a
`-0` member, `A == B` true on 92, `A < B` on 945, `A <= B` on 1094,
`A == A` on 1504 (NaN makes a box unequal to itself), `A < A` on 0 and
`A <= A` / `A >= A` on 2000 (upstream writes `<=` as `!(b < a)`, which is
true for NaN). The other twenty comparison cases have the same shape.

The member order matters and is checked: `OrientedBox` compares centre, then
the axis array, then the extent; `Rectangle` the same with two axes and a
2-vector extent; `Lozenge3` the whole rectangle then the radius; `Cylinder3`
the axis `Line3` then radius then height; `Sector2` compares `angle` but not
the derived `cosAngle`/`sinAngle`; `Parallelepiped3` and `Parallelogram2`
compare centre then axes and nothing else. `Segment` and `Triangle` use
`std::array`'s relational operators directly (so `<=` is `!(b < a)` at the
array level), which the port reproduces through its `comparePoints` /
`compareVertices` helpers.

The objects of the comparison cases are built by assigning public members
rather than through the value constructors, because `Parallelogram2` and
`Parallelepiped3` assert a right-handed basis, `Sector2`'s constructor calls
`std::cos`, and `Frustum3`'s calls `Update()`; none of those touch the
compared members.

### Branches the generators reach

Three populations per case, selected by `io.index()`: a small integer lattice
(mode 0), uniform (mode 1) and "wild" magnitudes (mode 2: signed zeros,
subnormals, `2^-1074`, `2^+-500`, `2^1020`). The comparison cases add mode 3,
wild magnitudes with about 25% NaN components. Every mode records the same
number of doubles, so the replay never needs to know which one ran.

* `AlignedBox`: a sixth of the records are flat (`min == max` on every axis)
  and a sixth are inverted (`min > max`), which upstream documents against but
  does not enforce.
* `CanonicalBox`, `OrientedBox`, `Rectangle`: zero extents are common on the
  lattice mode, so the degenerate box is covered.
* `Segment`: a fifth of the records have `p[0] == p[1]`, so
  `GetCenteredForm` normalizes the zero vector and takes the
  `length > 0` false arm (extent 0, direction zeroed).
* `Parallelepiped3` / `Parallelogram2`: the `getVertices` cases use a capped
  rejection loop that redraws every axis it tests until
  `DotCross(a0,a1,a2) > 0` (mode 0 draws integer axes, the other modes shear a
  right-handed orthonormal frame); the `rightHandedAssert` cases draw the axes
  unconstrained, and about 52% of the deep-run records throw.
* `Arc2.contains`: on mode 0 the centre is an integer point, the radius is 5
  and both endpoints and the test point are drawn from the twelve lattice
  points at distance exactly 5 from the origin, so `|P - C| - r` is exactly
  zero (1060 of 2000 deep-run records) and `DotPerp(P-E0, E1-E0)` is an exact
  integer that is exactly zero on 216 records. `epsilon` is drawn from
  `{-1, 0, 1e-9, 0.5}`; the negative value exercises the preserved upstream
  defect (issue #155). `Contains(P)` is true on 1093 of 2000 records and
  `Contains(P, epsilon)` on 527.
* `Sector2.contains`: mode 0 puts the vertex on the integer lattice, the
  direction on a coordinate axis, the points at distance exactly 5 and
  `cosAngle` in `{-1, -0.8, -0.6, -0.5, 0, 0.5, 0.6, 0.8, 1}`, for which
  `5*cosAngle` is exactly representable. Both `length <= radius` and
  `Dot(direction,diff) >= length*cosAngle` are therefore evaluated at exact
  equality: 790 of the 6000 deep-run point tests hit one of the two
  boundaries exactly.
* `Cylinder3.finiteInfinite`: `IsFinite` is true on 1054 of 2000 records
  before the sentinel is applied, and `MakeFiniteCylinder` is called with a
  negative height (which it ignores) on about half.
* `Polyhedron3.queries`: three shapes — a tetrahedron over four pool
  vertices with up to two unused pool entries, the eight corners of an
  axis-aligned box with its twelve triangles, and an arbitrary index list.
  414 of 2000 records have a unique-index count that is not a power of two,
  which is what makes `average /= n` (a multiplication by `1/n`)
  distinguishable from a per-component division.
* `Polyhedron3.invalidInput`: `numIndices` is drawn from 3..15, so 1697 of
  2000 records fail construction (`< 12` or not a multiple of 3) and report
  invalid with zero queries.
* `SegmentMesh.constructors`: `kind` selects the default constructor,
  DISJOINT, CONTIGUOUS_OPEN, CONTIGUOUS_CLOSED and INDEXED; the vertex count
  is drawn from 1..6 so the `>= 2` assert fires, the index-pair count from
  0..4 so the INDEXED `>= 1` assert fires, and when `validateIndices` is true
  the indices can be out of range so the third assert fires. About 22% of the
  deep-run records throw, with all three asserts represented.

### Association orders the cases actually pin down

A case that agrees bit for bit only proves something if a plausible
alternative spelling would have disagreed. Each of the suspect accumulations
was recomputed from the recorded inputs in the alternative grouping and
compared with the C++ outputs over the 2000-record deep run:

| computation | alternative | records that differ |
| --- | --- | --- |
| `OrientedBox::GetVertices` 3d, `((C -+ p0) -+ p1) -+ p2` | axis terms summed first | 1295 of 2000 |
| `OrientedBox::GetVertices` 2d | axis terms summed first | 938 of 2000 |
| `Parallelepiped3::GetVertices`, `((C -+ a0) -+ a1) -+ a2` | axis terms summed first | 711 of 2000 |
| `Parallelogram2::GetVertices` | axis terms summed first | 518 of 2000 |
| `Rectangle::GetVertices`, via `sum = p0+p1` and `dif = p0-p1` | axes added one at a time | 946 of 2000 |
| `Frustum3::ComputeVertices`, `(dS -+ uS) -+ rS` | `uS +- rS` grouped first | 1112 of 2000 |
| `Polyhedron3::ComputeVertexAverage`, `average *= 1/n` | per-component `/ n` | 270 of 2000 |

`AlignedBox::GetCenteredForm` is the one place where the alternative is not
distinguishable, and correctly so: upstream writes `(max + min) * half` and
`Vector.h`'s `operator/` would evaluate `(max + min) / 2` as a multiplication
by the exactly representable `0.5`, so the two are the same computation. The
case records that, rather than claiming a discrimination it does not have.

### Independent reference checks

Agreement between the port and the MSVC build is not correctness, so the
**C++ deep-run outputs** were checked against independent references
(2000 records per case; records whose own terms span more than about 12
decades are skipped, because a double-precision reference cannot second-guess
`centre +- extent` when the two differ by more than `2^53`):

| check | checked | violations | worst |
| --- | ---: | ---: | --- |
| `AlignedBox` 2d/3d: centre and extent are `(max +- min)/2`, each corner is `min` or `max` per axis, the corner mean is the centre | 2000 / 2000 | 0 | 7.9e-16 |
| `CanonicalBox` 2d/3d: every corner component is `+-extent[d]` with the bit-pattern sign | 2000 / 2000 | 0 | exact |
| `OrientedBox` 2d/3d, `Rectangle` 2d/3d: corner mean is the centre, corner 0 is `C - sum(e_a A_a)` and the opposite corner is `C + sum(e_a A_a)` | 1334-1340 | 0 | 1.3e-15 |
| `Parallelepiped3`, `Parallelogram2`: every corner is `C +- a0 +- a1 (+- a2)` and the corner mean is the centre | 1334 / 1339 | 0 | 1.5e-16 |
| `Frustum3`: all eight corners satisfy `Dot(P-E,D) = n` or `f`, `\|Dot(P-E,U)\| = (u/n)Dot(P-E,D)` and `\|Dot(P-E,R)\| = (r/n)Dot(P-E,D)`; the three derived ratios match `f/n`, `-2uf`, `-2rf` | 1334 | 0 | 3.2e-14 |
| `Torus3::Evaluate`: the returned position satisfies `[\|P-C\|^2+r0^2-r1^2]^2 - 4 r0^2 [\|P-C\|^2 - Dot(N,P-C)^2] = 0` | 2000 | 0 | 2.0e-15 |
| `Torus3::Evaluate`: `dX/du` and `dX/dv` against central differences of the parametric form | 1498 | 0 | 6.1e-11 |
| `Torus3::GetParameters`: `cos(u)` and `sin(v)` reproduce the recorded projections | 1998 | 0 | 5.0e-16 |
| `Segment` 2d/3d: the direction is unit length, the centre is the midpoint, `2*extent` is `\|p1-p0\|` and the round trip returns the endpoints | 1081 / 1080 | 0 | 8.8e-16 |
| `Polyhedron3::ComputeVolume` against the exact BigInt `\|sum DotCross\|/6` of the integer-lattice records | 667 | 0 | exact |
| `Polyhedron3::ComputeSurfaceArea` against `2(ab+bc+ca)` of the axis-aligned box records | 459 | 0 | 5.3e-16 |
| `Arc2::Contains` against the counterclockwise angle sweep from E0 to E1, on records where all three points are exactly on the circle and the test point is strictly off the chord | 484 | 0 | exact |
| `Sector2::Contains` against `acos` of the normalized dot product, on evaluations strictly off both boundaries | 3212 | 0 | exact |

All generated frames are right-handed: `RawFrame3` uses
`ComputeOrthogonalComplement`, whose third vector is `Cross(v0, v1)`, and
`RawFrame2` pairs `U = (x,y)` with `(-y,x) = -Perp(U)`, for which
`DotPerp(axis[0], axis[1]) = 1 > 0`. The `Frustum3` generator assigns
`rVector`, `uVector`, `dVector` in that order, matching the default frustum's
`R = Unit(0)`, `U = Unit(1)`, `D = Unit(2) = Cross(R, U)`.

## Port defects fixed

None. Every case is bit-identical to the MSVC build on every record of the
deep run apart from the three libm cases, whose worst scaled deviation is
`8.9e-16`. No `src/` file is changed by this branch.

## Deliberate deviations demonstrated

None; the group has no deliberate port fix. The three behavioural quirks of
these headers recorded in `docs/UPSTREAM-FINDINGS.md` are all *preserved* by
the port, and each is compared bit for bit rather than declared a deviation:

* `Arc2::Contains(P, epsilon)` with a negative epsilon returns `false` for
  every point, contrary to the comment's claim that it "behaves as if a value
  of zero was passed" (issue #155). `Arc2.contains` draws `epsilon = -1` on a
  quarter of its records and the two sides agree on all of them.
* `Torus3::Evaluate` guards the second-order block with `maxOrder == 2`
  rather than `>= 2`, so `maxOrder = 3` returns only the position and the two
  first derivatives (issue #484). `Torus3.evaluate` draws `maxOrder` from
  0..3 and emits the jet length first, so the quirk is pinned as a discrete
  output.
* `SegmentMesh`'s CONTIGUOUS_CLOSED constructor stores `S[i] = {i-1, i}`, a
  rotation of the documented `{i, (i+1) % L}`, and its DISJOINT constructor
  silently drops the last vertex for an odd vertex count (issue #78). Both
  index lists are emitted element by element and agree.

## Not covered

* `Polyhedron3::GetVertexPool` returns a `std::shared_ptr`, which has no
  observable value on the TypeScript side beyond nullness; that is covered by
  `operator bool` / `isValid()`. `GetVertices()` is covered through its
  element count in `Polyhedron3.queries`.
* `SegmentMesh`'s INDEXED constructor with `validateIndices = false` is only
  given in-range indices. Upstream stores the raw `size_t`, so a negative
  index would wrap to a value the port cannot represent; feeding one is
  undefined use of the class rather than a disagreement to measure.
* The comparison generators keep infinities out of the populations. They are
  not needed to separate the six relations, and an infinite member would make
  the arithmetic cases' outputs NaN, which the harness treats as universally
  equal.

## Upstream bug suspects

Two observations, neither of them result-corrupting on inputs the classes
document, so both are preserved by the port.

1. **`Parallelepiped3::GetVertices` and `Parallelogram2::GetVertices` are not
   `const`.** Every other `GetVertices` in this group
   (`AlignedBox`, `CanonicalBox`, `OrientedBox`, `Rectangle`) is declared
   `const`. The two newer classes declare
   `void GetVertices(std::array<Vector3<T>, 8>& vertices)` with no
   `const`, so the corners of a `Parallelepiped3 const&` cannot be computed
   without a copy. This is an API defect only: the function writes nothing to
   the object. The port has no `const` methods, so it does not reproduce it,
   and the oracle cases construct non-const objects. Not previously recorded
   in `docs/UPSTREAM-FINDINGS.md`.

2. **`Arc2::Contains` and `Sector2::Contains` use the non-robust `Length`,
   which overflows.** `Length(v)` is `sqrt(Dot(v,v))`, so a component larger
   than about `1.3e154` squares to `+inf` and the length is reported as
   infinite. `Sector2::Contains` then evaluates `inf <= radius` and returns
   `false` for a point that is well inside the sector. Explicit input from the
   deep run (record 344 of `Sector2.contains`):
   `vertex = (1.210067641901647e307, -0)`, `radius = 1.3116945087447793e307`,
   `direction = (0.899340699556372, -0.4372485633154847)`,
   `cosAngle = -0.07291478776539106`,
   `P = (1.9993459159770843e307, -5550230175632790)`. The true distance is
   `7.892782740754372e306 < radius` and `Dot(direction, P - vertex)` is
   `7.098e306`, far above `length*cosAngle = -5.755e305`, so the point is
   inside; upstream returns `false`. The same holds for `Arc2::Contains`.
   Both functions take no `robust` flag, so there is no exact separator to
   guard on and no way to fix this without changing the result on ordinary
   inputs: the port preserves the upstream expression and the oracle compares
   it bit for bit (the independent reference check skips these magnitudes and
   reports how many records it skipped). `Vector.h`'s *robust* `Length` has
   its own recorded defect for subnormal inputs (issue #370).
