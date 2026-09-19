# Oracle group 22 — `v22-distance`

Differential testing of the 18 headers of verify group 22 against the upstream
GTE build (MSVC 194435215, x64, `/O2 /fp:precise`, upstream commit
`d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`). See [ORACLE.md](../../ORACLE.md).

- `oracle/cpp/cases/v22-distance.cpp` — 28 cases
- `oracle/golden/v22-distance.txt` — 560 committed records (20 per case)
- `test/oracle/v22-distance.oracle.test.ts` — the TypeScript replay

**Result.** 23 ordinary cases and 5 cases that demonstrate deliberate fixes of
upstream defects. 21 of the 23 ordinary cases are declared `{ exact: true }`
and are **100% bit-identical** to the C++ build, at 20 records per case and at
2000 records per case in the deep run. The two exceptions are the two
`DistCircle3Circle3` cases on the general polynomial path, which are compared
with a measured scaled tolerance of `1e-10`; their root cause is established
by substitution below. **No port defect was found in this group** — every header
that does not inherit a documented deliberate fix agreed bit for bit on the
first run.

Deep run: `npm run oracle:deep -- 2000 v22` passes (56 000 records). The 21
exact cases contribute 412 818 real outputs, every one bit-identical.

Every generator mixes five modes selected from `io.index()`: uniform random
(modes 0, 1), a small integer lattice (modes 2, 3 — exact arithmetic, hence
coplanar and touching configurations, zero extents, degenerate triangles and
tetrahedra, and box axes that are the coordinate axes so that lattice inputs
keep exactly zero components in the box frame) and a construction aimed at
contact (mode 4 — the plane through the box center, the rectangle centered in
the box, a triangle vertex inside the box, two boxes sharing a center, the
query point at the tetrahedron centroid or at the parallelepiped center,
concentric circles). Every mode records the same number of doubles, so the
replay reads the inputs without knowing which mode produced them.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `DistAlignedBox3OrientedBox3.h` | `compute` | exact | pass |
| `DistCircle3Circle3.h` | `compute`, `compute.identityRotation`, `compute.parallelPlanes`, `compute.doubleRoot.deviation`, `compute.antiParallel.deviation`, `compute.coaxial.deviation` | tol 1e-10 / exact / deviation | pass |
| `DistOrientedBox3Cone3.h` | `compute.control`, `compute.deviation`, `compute.infinite` | exact / deviation | pass |
| `DistOrientedBox3OrientedBox3.h` | `compute` | exact | pass |
| `DistPlane3AlignedBox3.h` | `compute` | exact | pass |
| `DistPlane3OrientedBox3.h` | `compute` | exact | pass |
| `DistPoint3Parallelepiped3.h` | `compute`, `getMinimizer`, `ctor.leftHanded` | exact | pass |
| `DistPoint3Tetrahedron3.h` | `compute` | exact | pass |
| `DistRectangle3AlignedBox3.h` | `compute` | exact | pass |
| `DistRectangle3CanonicalBox3.h` | `compute` | exact | pass |
| `DistRectangle3OrientedBox3.h` | `compute` | exact | pass |
| `DistRectangle3Rectangle3.h` | `compute` | exact | pass |
| `DistTetrahedron3Tetrahedron3.h` | `compute`, `compute.nested` | exact | pass |
| `DistTriangle3AlignedBox3.h` | `compute` | exact | pass |
| `DistTriangle3CanonicalBox3.h` | `compute` | exact | pass |
| `DistTriangle3OrientedBox3.h` | `compute` | exact | pass |
| `DistTriangle3Rectangle3.h` | `compute` | exact | pass |
| `DistTriangle3Triangle3.h` | `compute` | exact | pass |

Entry points covered: every `operator()` of every header of the group, both
overloads of `DistOrientedBox3Cone3::operator()` (the default controls in
`compute.deviation`, the `Control` overload in `compute.control`) and the one
public member function the group exposes besides `operator()`,
`DistPoint3Parallelepiped3::GetMinimizer`, which has its own case because its
27-region dispatch tree is not reachable from `operator()` alone with a
uniform generator. `Parallelepiped3`'s right-handedness assertion gets the
throw-parity case `ctor.leftHanded` (C++ throws on 20 of 20 records and the
port must throw on the same ones). The port exports no free function for any
header of this group, and none of these headers is pure data.

## Comparison policy

21 of the 23 ordinary cases are `{ exact: true }`: their paths are
`+ - * /`, `sqrt`, `fabs`, `min`, `max` and comparisons only. That includes
`DistOrientedBox3Cone3.compute.control` and
`DistCircle3Circle3.compute.parallelPlanes`, both of which are arranged so
that no libm call other than `sqrt` is evaluated (see below).

`DistCircle3Circle3.compute` and `DistCircle3Circle3.compute.identityRotation`
use a scaled tolerance of `1e-10`. The reason is the port's documented
fix of issue [#331](https://github.com/gradientspaceai/gtengine-js/issues/331):
upstream recovers the sine of a critical angle as `sn = -p6(cs)/p7(cs)`, which
equals `+/-sqrt(1 - cs^2)` only at an *exact* root of `phi`; the port uses both
signs of `sqrt(1 - cs^2)` and adds further candidates. This is a
**conditioning** fix — there is no exact predicate separating the roots a
floating-point bisection locates well from the ones it does not — so the main
case keeps the configurations on which upstream's own answer is right and
measures the residual, and the regime where upstream is actually wrong has its
own deviation case.

The generator's acceptance threshold and the comparison tolerance are tied
together deliberately. A record is kept only when upstream's own reported
distance agrees with an independent reference (the minimum over circle1 of the
point-to-circle0 distance, resolved by a 4096-sample scan plus 200
golden-section steps, accurate to about `1e-13` relative) to `1e-12`
relative; when nothing reaches that threshold the record keeps the best of
the 200 candidates rather than the last. On such a record upstream is right to
`1e-12` and the port is at least as right, so a disagreement much larger than
that band is a finding rather than the conditioning of the substituted
expression. The tolerance is set above the acceptance threshold because the
closest-point coordinates are more sensitive than the distance the generator
screens on.

Both numbers were arrived at by measurement, not by taste. An earlier draft
accepted at `1e-9` and compared at `1e-12`; the deep run correctly rejected
it, with 10 of 2000 records of `compute.identityRotation` differing by up to
`2.7e-10` — exactly the band the acceptance threshold had left open.
Tightening the acceptance to `1e-12`, and keeping the best rather than the
last candidate when nothing reaches it, brought the worst residual over 2000
records to `7.6e-12` for `compute` and `2.3e-12` for
`compute.identityRotation`. The `1e-10` tolerance is an order of magnitude
above those, which is headroom for the conditioning rather than for an
unexplained difference.

The root cause was established rather than assumed. Temporarily substituting
upstream's `sn = -p6(cs)/p7(cs)` and upstream's candidate set back into
`src/DistCircle3Circle3.ts` makes **every real output of both cases
bit-identical** to the C++ build (160 of 160 in each over the 20 records per
case the experiment ran on, maximum scaled error 0), and makes
`compute.doubleRoot.deviation` stop deviating entirely. Nothing else on
the path differs. `compute.identityRotation` exists to separate the two
possible sources: it forces the larger circle's normal to be exactly `(0,0,1)`,
so `PrepareCircles`'s `std::acos(1)` is `0` and `Rotation` of the zero
axis-angle is exactly the identity matrix, which leaves `sqrt` as the only
libm call on the path. Its residual is therefore the substitution alone.

Measured residuals over the 20 committed records per case:
`compute.identityRotation` `3.44e-15` scaled with 132 of 160 real outputs
still bit-identical, `compute` `1.25e-15` with 110 of 160 bit-identical. Over
2000 deep-run records per case: `compute` `7.58e-12` with 11331 of 16000 real
outputs bit-identical, `compute.identityRotation` `2.27e-12` with 12443 of
16000 bit-identical.

### Why `DistOrientedBox3Cone3` needs a constructed `Control`

This query minimizes `F(angle)`, the distance between the box and the planar
quadrilateral cut from the frustum at that angle, with `Minimize1`. Two
separate things make the default search incomparable, and they are independent
of `DistOrientedBox3Cone3.h` itself:

1. The port fixes two defects of `Minimize1::GetBracketedMinimum`
   (UPSTREAM-FINDINGS `Minimize1.h` item 1, issue
   [#298](https://github.com/gradientspaceai/gtengine-js/issues/298)). For this
   query the defective regime is the *normal* one, because the two endpoint
   slices of the sweep are the same point set, so `F(-pi/2)` and `F(+pi/2)`
   agree to a few ulps and the parabola vertex lands within round-off of the
   middle sample. 16 of 20 committed records of `compute.deviation` differ,
   by up to `0.89` scaled.
2. `F` calls `std::cos` and `std::sin`. **MSVC's and V8's differ by one ulp on
   a couple of percent of their arguments**, and `sin(+/-pi/4)` — one of the
   angles the subdivision visits — is one of them: MSVC returns
   `0x3fe6a09e667f3bcd`, V8 returns `0x3fe6a09e667f3bcc`. A perturbation of
   `F` at one sampled angle flips a strict comparison in the search and sends
   it to a different local minimum, so the difference is not a rounding
   difference a tolerance could absorb.

`compute.control` removes both. It exercises the `Control` overload with
`maxSubdivisions = 1`, `maxBisections = 1`, `epsilon = 10` and
`tolerance = 0`, and its generator keeps only configurations whose initial
polyline `{(-pi/2,f0),(0,fm),(pi/2,f1)}` is V-shaped — evaluated in C++ by a
probe that is upstream's own `DoBoxQuadQuery` copied verbatim. `GetMinimum`
then evaluates `F` at exactly `-pi/2`, `0` and `+pi/2` and calls
`GetBracketedMinimum`, whose first act is the convergence test
`pi <= 2*0*|tm| + 10`, so it returns without a further evaluation. The answer
is the smallest of three box-quadrilateral distances at three fixed angles at
which the two libm implementations agree bit for bit (`cos(0) = 1`,
`sin(0) = 0`, `sin(+/-pi/2) = +/-1`,
`cos(+/-pi/2) = 6.123233995736766e-17`). What the case measures is then
`DistOrientedBox3Cone3`'s own arithmetic — the 5x5 quadratic form, the 10-D
LCP and the mapping of the LCP solution back to the two closest points — and
it is bit-identical on all 20 committed and all 2000 deep-run records (14 000
real outputs).

The default-controls path is therefore covered only as a deviation. That is
stated here rather than papered over: this group does not establish that the
port's box-cone answer with the default controls matches upstream's, only
that the per-angle subproblem does and that the difference comes from the
`Minimize1` fix.

## Port defects fixed

None. No output of any header of this group disagreed with the C++ build for
a reason other than a deliberate, already documented fix; in particular no
reassociated change of basis was found. The three headers of this group that
map closest points out of the box frame — `DistPlane3OrientedBox3`,
`DistRectangle3OrientedBox3` and `DistTriangle3OrientedBox3` — already
accumulate `closest[i] += result.closest[i][j] * box.axis[j]` left to right
as upstream does, and all three are bit-identical.

## Deliberate deviations demonstrated

Each of these cases generates only inputs on which upstream is defective; the
harness requires the port to disagree with the C++ build on at least one
record, so the fix is demonstrated against the real build rather than against
a TypeScript re-derivation. Two of the five select by construction rather than
by a probe: an anti-parallel pair of normals orthogonal to the z-axis, and a
coaxial pair, always reach the misclassified branch, and there the reported
*distance* is often right while the reported closest points are not, so the
distance is not a usable selector.

| case | records disagreeing | upstream defect |
| --- | ---: | --- |
| `DistOrientedBox3Cone3.compute.deviation` | 16 of 20 | `Minimize1::GetBracketedMinimum` collapses its bracket whenever the parabola vertex lands within round-off of the middle sample, which is the normal situation for this query, and then returns on the degenerate parabola of the next iteration (issue [#298](https://github.com/gradientspaceai/gtengine-js/issues/298)). Max scaled error `0.89` over the committed records, `2.0` over 2000 deep-run records, where 1269 of 2000 deviate. |
| `DistOrientedBox3Cone3.compute.infinite` | 20 of 20 (the port throws on all 20) | The header documents `hmax < infinity` but `operator()` never checks it. `Cone` encodes an infinite cone as `maxHeight = -1`, which upstream feeds into the LCP as the infeasible constraint `z3 + z4 <= -1` and answers anyway. The port asserts the precondition (issue [#298](https://github.com/gradientspaceai/gtengine-js/issues/298)). The generator builds the cone with `Cone(ray, angle, minHeight)`, the documented infinite-truncated-cone constructor. |
| `DistCircle3Circle3.compute.doubleRoot.deviation` | 14 of 20 | A mirror-symmetric configuration, constructed directly in the frame `PrepareCircles` produces (circle1 the unit circle in `z = 0`, circle0's center on `x = 0` and its normal with an exactly zero x-component). The reflection `x -> -x` maps each circle to itself, so `p6` and `p7` acquire a factor `cs^2` and `phi` a root of multiplicity four at `cs = 0` — the critical angle that matters. A sign-change bisection cannot find an even-multiplicity root, so upstream never considers it and evaluates `sn = -p6/p7` at some other root, off the unit circle (issue [#331](https://github.com/gradientspaceai/gtengine-js/issues/331)). The generator aims the construction at the regime where the missed critical angle is the global minimum (circle0 offset in +y with its plane nearly perpendicular to circle1's), and keeps, of 24 candidates, the one whose reported distance is furthest from the independent reference. Max scaled error `9.2e-9` over the committed records and `1.4e-2` over 2000 deep-run records, where 1585 of 2000 records deviate. Random sampling never reaches this regime: the symmetry has to be constructed. |
| `DistCircle3Circle3.compute.antiParallel.deviation` | 7 of 20 | `PrepareCircles` aligns normals by negating any whose z-component is negative, which does nothing when both normals are orthogonal to the z-axis. Two circles in parallel planes with anti-parallel normals therefore arrive with `circle0.normal = (0,0,-1)`, the test `circle0.normal[2] < 1` misclassifies them as skew and the general polynomial path answers the wrong question (issue [#431](https://github.com/gradientspaceai/gtengine-js/issues/431)). Max scaled error `1.8` over the committed records, `2.0` over 2000 deep-run records, where 769 of 2000 deviate. |
| `DistCircle3Circle3.compute.coaxial.deviation` | 19 of 20 (C++ throws on one of them, where the port answers) | Coaxial circles in parallel planes. Whenever the rotation leaves the transformed `circle0.normal[2]` a few ulps below one, upstream takes the polynomial path, where `a1 = 0` makes `phi` and `p6` perfect squares: the sign-change bisection finds no root, upstream reads a default-constructed candidate and reports distance 0 at the origin — or `p6` and `p7` vanish identically and it throws "Unexpected degree for p6" (issue [#431](https://github.com/gradientspaceai/gtengine-js/issues/431)). The port answers with the parallel-planes code and the exact common normal (issue [#442](https://github.com/gradientspaceai/gtengine-js/issues/442)). Max scaled error `1.9` over the committed records, `2.0` over 2000 deep-run records, where 1770 of 2000 deviate and upstream throws on 11. |

## Branch histogram

Counts from the deep run (2000 records per case), the ones that show the
non-generic branches are reached.

- **The solids meet, so the distance is exactly zero**, in 2000 records:
  854 for `DistPlane3AlignedBox3`, 1047 for `DistPlane3OrientedBox3`, 711 for
  `DistRectangle3CanonicalBox3`, 510 for `DistRectangle3AlignedBox3`, 581 for
  `DistRectangle3OrientedBox3`, 1311 for `DistTriangle3CanonicalBox3`, 583 for
  `DistTriangle3AlignedBox3`, 725 for `DistTriangle3OrientedBox3`, 492 for
  `DistRectangle3Rectangle3`, 545 for `DistTriangle3Rectangle3`, 860 for
  `DistTriangle3Triangle3`, 681 for `DistOrientedBox3OrientedBox3`, 570 for
  `DistAlignedBox3OrientedBox3`, 427 for `DistPoint3Tetrahedron3` (the point
  inside the solid, where no face is visible and the `invalid` sentinel
  survives the loop), 1414 for `DistTetrahedron3Tetrahedron3` (the
  `foundZeroDistance` early exit), 950 for its `.nested` case (the centroid
  containment branch, the only one that reports zero without a zero face-pair
  distance) and 469 for `DistPoint3Parallelepiped3` (the `Rzzz` region, where
  the minimizer is `Z` itself).
- **Degenerate inputs**: `DistPoint3Tetrahedron3` reports a NaN distance on 4
  of 2000 records and `DistTetrahedron3Tetrahedron3` and
  `DistTriangle3Triangle3` on one each — coplanar or collinear lattice inputs
  where upstream normalizes a zero cross product or divides by a zero squared
  normal length. The port reproduces every one of them bit for bit, including
  the NaN.
- `DistPoint3Parallelepiped3.getMinimizer` visits all 27 regions of its
  dispatch tree (`io.index() % 27`), and within the interior arm the
  coordinate is exactly `-1` or exactly `+1` about half the time, so the
  `Z[i] <= PosOne()` comparisons are evaluated at equality.
- `DistPoint3Parallelepiped3.ctor.leftHanded`: C++ throws on 2000 of 2000
  records and the port throws on the same ones.
- `DistCircle3Circle3.compute.parallelPlanes` reaches all five branches of
  `DoQueryParallelPlanes`: **803 of 2000 records report two closest pairs**
  (the overlapping-circles branch), 1197 report one, **410 report
  `equidistant`** (the concentric branches, where `U` is rebuilt from
  `GetOrthogonal`) and 75 have distance exactly zero (coplanar overlapping
  circles).
- `DistCircle3Circle3.compute` and `.identityRotation` report a single closest
  pair on every record, which is what a non-symmetric configuration gives;
  the two-pair result is covered by `parallelPlanes` above.
- `DistCircle3Circle3.compute.coaxial.deviation`: upstream **throws** on 11 of
  2000 records ("Unexpected degree for p6") and returns a defaulted
  zero-distance result on others; the port answers every one of them.
- `DistOrientedBox3Cone3.compute.deviation` finds a zero distance (the box
  meets the frustum) on 54 of 2000 records and `compute.control` on 6; the
  `.control` generator's rejection of non-V-shaped polylines biases it away
  from deep overlap, which is stated here because it is a real narrowing of
  that case.

## Upstream bug suspects

None new. Every defect exercised above is already recorded in
`docs/UPSTREAM-FINDINGS.md`.

One measurement is worth recording because it is a property of the two
runtimes rather than of either library, and it limits what any oracle case
over a libm-driven search can assert: **MSVC's `std::sin`/`std::cos` and V8's
`Math.sin`/`Math.cos` differ by one ulp on a couple of percent of their
arguments.** In a probe of 180 arguments spread over `[-pi/2, pi/2]` (360
cosine and sine values) five results differed, among them `sin(pi/4)` (`0x3fe6a09e667f3bcd` in MSVC,
`0x3fe6a09e667f3bcc` in V8) — an argument any dyadic subdivision of
`[-pi/2, pi/2]` visits. Neither value is wrong; both are within half an ulp of
the true sine. Any case whose control flow compares two libm-derived
quantities therefore has to be arranged so that the comparison is not
decided by that last bit, as `DistOrientedBox3Cone3.compute.control` is.

## Not covered

Nothing. All 18 headers of the group are implemented by the port and all have
cases.

Two narrower notes, not gaps in entry-point coverage:

- `DistOrientedBox3Cone3` with the **default** controls is compared only as a
  deviation, for the reasons given under "Comparison policy". The `Control`
  overload carries the bit-exact comparison of the query's own arithmetic.
- `DistCircle3Circle3`'s second closest pair is emitted by
  `compute.parallelPlanes` (overlapping circles always report two pairs),
  `compute` and `compute.identityRotation`. The three deviation cases emit
  the first pair only, because there the port's candidate set is deliberately
  a superset of upstream's and a tie it detects need not be a tie upstream
  detects, which would turn a value difference into a structural one and hide
  the measured deviation.
