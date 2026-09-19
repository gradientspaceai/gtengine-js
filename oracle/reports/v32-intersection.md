# Group 32 (`v32-intersection`) — C++ oracle report

64 cases, 20 records each in `oracle/golden/v32-intersection.txt` (1280 records).
55 of the cases are ordinary comparisons, 9 are deliberate `deviation` cases and
one case (`IntrAlignedBox3Cylinder3.test.infinite`) records a C++ assertion
failure that the port must reproduce. The 55 ordinary cases carry 3683
floating-point outputs in the goldens, of which 3666 are bit-identical to the
MSVC build of upstream GTE; the 17 that are not all belong to the two
`IntrEllipse2Ellipse2` find-intersection cases (worst scaled error 9.1e-14 in
the goldens).

Deep run `npm run oracle:deep -- 2000 v32`: 64 cases, 128000 records, 2000
recorded C++ exceptions (the `test.infinite` case), 367448 floating-point
outputs compared over the 55 ordinary cases, **99.371 % bit-identical**, and
every case within its declared comparison. Only two cases are not bit-identical
on every output:

| case | bit-identical | worst scaled error | why |
| --- | --- | --- | --- |
| `IntrEllipse2Ellipse2.find` | 14739 / 16000 | 2.04e-10 | `RootsPolynomial::SolveQuartic` (`pow`, `cos`, `atan2` in the resolvent cubic) |
| `IntrEllipse2Ellipse2.findStandardForm` | 14951 / 16000 | 4.26e-10 | the same |

The other 53 ordinary cases are bit-identical on all 335448 of their
floating-point outputs over 106000 deep-run records. All discrete outputs
(`numIntersections`, `numPoints`, `isTransverse`, `configuration`,
`intersectionType`, the sphere/sphere `type`, the ellipse classification) are
compared exactly and agreed on every record of every case.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise`.

## Coverage

Every computational entry point of every header in the group is covered: each
`TIQuery::operator()` and `FIQuery::operator()` overload (including the named
`Exact` variants and the standard-form overload), each protected `DoQuery`
helper that the port exports as a free function (called directly on both sides
through a small `Expose` subclass in C++), and the public helpers
`GetStandardForm` and both `ComputeAlignedBox` overloads.

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `IntrAlignedBox2Circle2.h` | `test`, `find`, `find.aimed` | exact | pass |
| `IntrAlignedBox3Cylinder3.h` | `test`, `test.infinite`, `test.edgeTypoDeviation` | exact (+1 deviation) | pass |
| `IntrAlignedBox3Sphere3.h` | `test`, `find`, `find.aimed`, `find.probeDeviation` | exact (+1 deviation) | pass |
| `IntrConvexPolygonHyperplane.h` | `test.2d`, `find.2d`, `test.3d`, `find.3d` | exact | pass |
| `IntrEllipse2Ellipse2.h` | `test`, `find`, `findStandardForm`, `getStandardForm`, `test.poleDeviation`, `find.divisorDeviation` | exact / tol 1e-8 (libm) (+2 deviations) | pass |
| `IntrRay3Sphere3.h` | `test`, `find`, `doQuery.fi`, `find.tangent` | exact | pass |
| `IntrSegment2AlignedBox2.h` | `test`, `find`, `find.degenerate`, `doQuery.ti`, `doQuery.fi` | exact | pass |
| `IntrSegment2Segment2.h` | `test`, `testExact`, `find`, `findExact`, `test.collinear`, `testExact.collinear`, `find.collinear`, `findExact.collinear`, `find.antiparallelDeviation` | exact (+1 deviation) | pass |
| `IntrSegment2Triangle2.h` | `test`, `find`, `doQuery.fi`, `find.throughPoint` | exact | pass |
| `IntrSegment3AlignedBox3.h` | `test`, `find`, `doQuery.ti`, `doQuery.fi`, `find.throughPoint` | exact | pass |
| `IntrSegment3Rectangle3.h` | `test`, `find`, `find.throughPoint` | exact | pass |
| `IntrSegment3Sphere3.h` | `test`, `find`, `doQuery.fi`, `doQuery.fi.tangent`, `test.containedDeviation` | exact (+1 deviation) | pass |
| `IntrSphere3Sphere3.h` | `test`, `find`, `find.tangent`, `find.internalTangentDeviation` | exact (+1 deviation) | pass |
| `IntrSphere3Triangle3.h` | `find`, `find.aimed` | exact | pass |
| `IntrTriangle3Cylinder3.h` | `test`, `test.degeneratePolygonDeviation`, `test.infiniteDeviation` | exact (+2 deviations) | pass |

Only the two `IntrEllipse2Ellipse2` find-intersection cases reach the C math
library on the compared path. Everywhere else the `sin`/`cos` used to build
frames, convex polygons and unit directions, and the `sqrt` inside `Normalize`,
are applied to *unrecorded* draws; only the resulting frame, vertex or unit
vector is recorded as an input, so libm never enters the compared computation.

### Generators and the branches they reach

Each generator alternates an exactly representable small-lattice mode with a
uniform mode; several add a constructed mode aimed at a specific branch. The
distributions below are from the 2000-record deep run.

- **Rays, segments and spheres.** `IntrRay3Sphere3.find.tangent` and
  `IntrSegment3Sphere3.doQuery.fi.tangent` place the origin at
  `C + (R+d)*e_i + s*e_j` with the unit direction `e_j`, which makes the
  discriminant exactly `R^2 - (R+d)^2`: the tangent branch
  (`numIntersections == 1`) occurs on 408 and 571 of 2000 records against 41
  and 47 for the uniform variants, and `s = 0` puts a root exactly at the ray
  origin (respectively at the segment centre), so the `[0,+infinity)` and
  `[-e,e]` interval clips are evaluated at exact equality.
- **Segments versus boxes.** `IntrSegment3AlignedBox3.find.throughPoint` aims
  the segment at `C + (k/2)*e_i` on a lattice box with an integer direction and
  dyadic endpoint parameters, which raises the hit rate from 444 to 1960 of
  2000 records and makes 916 of them single-point contacts (the Liang-Barsky
  clip collapsing to `t0 == t1`). `IntrSegment2AlignedBox2.find.degenerate`
  drives the zero-length-segment branch, which upstream answers with a
  containment test: 448 of 2000 records report the contained point.
- **Segments versus triangles and rectangles.**
  `IntrSegment2Triangle2.find.throughPoint` aims at `v0 + (k1*E1 + k2*E2)/4` on
  an integer triangle: 236 of 2000 records report a single-point intersection
  (against 144 for the uniform variant), which is the case where the segment
  touches a vertex or an edge exactly.
  `IntrSegment3Rectangle3.find.throughPoint` aims at
  `C + (a/2)*e0*W0 + (b/2)*e1*W1` on an axis-aligned lattice rectangle, so the
  corner and edge cases of `|rectCoord| <= extent` are exact; 611 of 2000
  records hit, against 242 for the uniform variant.
- **Segments versus segments.** The four `*.collinear` cases put both segments
  on one lattice line, which reaches the "same line" branch of `IntrLine2Line2`
  exactly: the overlap-in-a-segment (934/2000), single-point touch (306) and
  disjoint (760) branches all occur for `find.collinear`, and the `Exact`
  variants likewise. The endpoint order of the second segment is flipped on
  half the records of the `Exact` and `test` variants, which is what makes the
  two centred directions antiparallel.
- **Spheres versus spheres.** `find.tangent` uses `C1 = C0 + d*e_k` on the
  lattice, so `|C1-C0|^2` is exactly `d^2`: external tangency (type 1, 654
  records), internal tangency with sphere1 inside sphere0 (type 6, 668) and
  identical spheres (type 4 with a zero direction, 678) all occur exactly. The
  uniform `find` case reaches types 0, 1, 2, 3, 5 and 6; type 4 is the
  deviation case.
- **Dynamic queries (box/circle, box/sphere, sphere/triangle).** The `aimed`
  variants send the moving object towards a dyadic point of the target's
  boundary, which reaches the rounded face, edge and vertex branches of the
  Minkowski sum: `IntrAlignedBox3Sphere3.find.aimed` reports a future contact
  on 1566 of 2000 records against 256 for the uniform variant, and
  `IntrSphere3Triangle3.find.aimed` on 950 against 249. All three
  `intersectionType` values occur in every one of these cases.
- **Convex polygon versus hyperplane.** The polygon is a cyclically increasing
  subset of the twelve lattice points on the circle of radius 5 (exactly
  convex, exact arithmetic) or a set of increasing angles on a circle; in 3D
  the lattice mode puts it in a coordinate plane so the plane normal is exact
  too. The hyperplane is drawn uniformly, or built through one polygon vertex,
  or from the perpendicular of one polygon edge, or (3D) as the plane of the
  polygon. All nine `Configuration` values are reached: in 3D `SPLIT` 455,
  `CONTAINED` 311, `NEGATIVE_SIDE_EDGE` 196, `POSITIVE_SIDE_VERTEX` 187,
  `NEGATIVE_SIDE_VERTEX` 134, `NEGATIVE_SIDE_STRICT` 130, `POSITIVE_SIDE_STRICT`
  125, `INVALID_POLYGON` 401 and `POSITIVE_SIDE_EDGE` the rest; in 2D the same
  except `CONTAINED`, which a nondegenerate 2D polygon cannot reach.
- **Triangle versus cylinder.** The lattice mode uses a coordinate-axis
  cylinder direction and an even height, so the projected z-coordinates are
  exact integers and `h/2` is an integer: the clip cases 1a, 1b, 2a, 2b (the
  triangle touching a slab plane in a point or a segment) are reached along
  with every polygon clip case. Both answers occur (1193/807).
- **Ellipses.** The lattice mode uses coordinate axes, lattice centres and
  extents that are powers of two, so the standard-form matrices are exact and
  the common-centre branch (`K == 0`), the `d0 == d1` branch and the
  `e4 == 0` case analysis are reached; the uniform mode rotates both ellipses,
  which is what makes `e4 != 0` and drives `CaseE4NotZero`. All eight
  `Classification` values occur in `test` (`ELLIPSES_SEPARATED` 1257,
  `ELLIPSES_OVERLAP` 682, then the six containment and tangency values), and
  `find` reaches `numPoints` 0, 1, 2, 3, 4 and the "identical ellipses"
  sentinel.

## Port defects fixed

None. Every disagreement found in this group was either a deliberate,
documented port fix of an upstream defect (nine of them, listed below) or C
math library rounding in the quartic root finder. No case required a change to
`src/`.

## Deliberate deviations demonstrated

| case | records disagreeing (golden / deep) | record of the decision |
| --- | --- | --- |
| `IntrSegment3Sphere3.test.containedDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrSegment3Sphere3.h` TIQuery; issue [#203](https://github.com/gradientspaceai/gtengine-js/issues/203) |
| `IntrSphere3Sphere3.find.internalTangentDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrSphere3Sphere3.h` FIQuery type-4 branch; issue [#203](https://github.com/gradientspaceai/gtengine-js/issues/203) |
| `IntrSegment2Segment2.find.antiparallelDeviation` | 19/20, 1926/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrSegment2Segment2.h` FIQuery; issue [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) |
| `IntrAlignedBox3Cylinder3.test.edgeTypoDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrCanonicalBox3Cylinder3.h` `DoQueryNoZeros` `(U1,-D)` sign typo; issue [#197](https://github.com/gradientspaceai/gtengine-js/issues/197) |
| `IntrAlignedBox3Sphere3.find.probeDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrAlignedBox3Sphere3.h` `DoQueryRayRoundedFace`/`DoQuery`; issues [#458](https://github.com/gradientspaceai/gtengine-js/issues/458), [#465](https://github.com/gradientspaceai/gtengine-js/issues/465) |
| `IntrTriangle3Cylinder3.test.degeneratePolygonDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrTriangle3Cylinder3.h` `DiskOverlapsPolygon`; issues [#206](https://github.com/gradientspaceai/gtengine-js/issues/206), [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) |
| `IntrTriangle3Cylinder3.test.infiniteDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` cylinder headers, missing `IsFinite` guard; issues [#197](https://github.com/gradientspaceai/gtengine-js/issues/197), [#206](https://github.com/gradientspaceai/gtengine-js/issues/206), [#255](https://github.com/gradientspaceai/gtengine-js/issues/255) |
| `IntrEllipse2Ellipse2.test.poleDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrEllipse2Ellipse2.h` TIQuery; issue [#458](https://github.com/gradientspaceai/gtengine-js/issues/458) |
| `IntrEllipse2Ellipse2.find.divisorDeviation` | 20/20, 2000/2000 | `docs/UPSTREAM-FINDINGS.md` `IntrEllipse2Ellipse2.h` `CaseE4NotZero`, `divisor == 0`; issue [#250](https://github.com/gradientspaceai/gtengine-js/issues/250) |

Each deviation is confined to inputs on which upstream is actually defective,
and the corresponding main case is generated by rejection sampling against an
explicit predicate (or, where the predicate is a closed form, by construction)
rather than by narrowing the input ranges:

- **`IntrSegment3Sphere3`.** Upstream's `qm * qp <= 0` test misses the case
  where *both* endpoint values of `Q` are negative, that is, where the whole
  segment lies strictly inside the solid sphere; it then falls through to
  `qm > 0 && |a1| < e` and reports no intersection, contradicting its own FI
  query. The predicate `Q(-e) < 0 && Q(e) < 0` (evaluated on upstream's own
  expressions, after upstream's `discr < 0` early exit) is exactly the set on
  which the port deviates; the main `test` case rejects it and the deviation
  case draws both endpoints inside the ball of radius `radius/2`.
- **`IntrSphere3Sphere3`.** The type-4 branch (internal tangency with
  `r0 <= r1`) reports `C1 + r1*(C1-C0)/|C1-C0|`, the antipode of the contact
  point. The main `find` case rejects `|C1-C0|^2 == (r0-r1)^2 && r0 <= r1`
  exactly; the deviation case constructs it on the lattice with
  `C1 = C0 + d*e_k` and `r1 = r0 + d`. Identical spheres also take the type-4
  branch, but there `C1-C0` is the zero vector and the two expressions agree,
  so they stay in the main `find.tangent` case.
- **`IntrSegment2Segment2`.** The collinear branch of the centred-form FI query
  reports `segment1Parameter[i] = overlap[i] - t`, a parameter measured along
  *segment0*'s direction, so `C1 + segment1Parameter[i] * D1` is not the
  reported point when the two centred directions are antiparallel. The probe
  replicates upstream's control flow (line-line query, then the interval clip)
  and reports whether the same-line branch is taken with
  `Dot(D0,D1) < 0` and a nonempty overlap; the main `find` case rejects it and
  the deviation case builds two overlapping collinear segments with the second
  one's endpoints in reverse order. 74 of the 2000 deviation records agree
  anyway: there `DotPerp` of the (normalized) direction with the difference of
  the centres rounds to a nonzero value, so upstream classifies the two lines
  as parallel-but-distinct and both sides report no intersection.
- **`IntrAlignedBox3Cylinder3`.** The query is a thin wrapper around
  `IntrCanonicalBox3Cylinder3`, so it inherits the `(U1,-D)` sign typo of
  `DoQueryNoZeros` and the port's correction of it (first demonstrated by
  group 31). The case file carries the same corrected `DoQueryNoZeros` that
  `oracle/cpp/cases/v31-intersection.cpp` uses, applied to the translated
  cylinder and the canonical box the wrapper builds; the main case keeps the
  inputs where the corrected answer equals upstream's, the deviation case those
  where it does not. Without the rejection the main case disagreed on 3 of 2000
  deep-run records.
- **`IntrAlignedBox3Sphere3`.** Upstream probes at most one rounded piece of the
  Minkowski sum per case and accepts the first probe that reports a contact, so
  it reports a contact later than the true first contact or misses one
  entirely; the port probes every candidate piece and keeps the earliest. The
  case file carries the port's corrected `DoQuery` (upstream's code with the
  three marked changes), and the two main cases keep only the inputs where the
  corrected answer equals upstream's bit for bit. On the deviation case
  upstream misses the contact altogether on 14 of 2000 records and reports a
  later contact time on the rest.
- **`IntrTriangle3Cylinder3`.** `DiskOverlapsPolygon` reports containment of the
  origin when *every* `DotPerp(Q[i0], Q[i0]-Q[i1])` is zero, which is not
  containment but degeneracy: all polygon vertices lie on one line through the
  origin. Since `DotPerp(A, A-B) = -DotPerp(A,B)`, "every term is zero" means
  consecutive vertices are parallel, hence zero polygon area, and every polygon
  the query builds is a convex combination of the three projected triangle
  vertices; the main case therefore rejects a zero-area projected triangle,
  which is an exact predicate and a superset of the deviating inputs. The
  deviation case puts the triangle in a plane containing the cylinder axis at
  radial distance `r + 1` or more, so upstream reports an intersection for a
  triangle that is nowhere near the cylinder. The second deviation is the
  missing `IsFinite` guard: `Cylinder3` marks an infinite cylinder with the
  sentinel `height = -1`, which upstream feeds into the finite-cylinder
  formulas (an empty slab), where the port asserts.
- **`IntrEllipse2Ellipse2`.** The TI query drops the `c_i = 0` term of `f(s)`
  when exactly one of `c0`, `c1` is zero, which loses the two critical points at
  the pole `s = 1/d_j`. The probe replicates upstream's control flow down to the
  `valid` list (the eigensolver included) and reports whether a term is dropped;
  the main `test` case rejects it. The deviation case is the documented example
  family: the unit circle against an axis-aligned ellipse whose centre is offset
  along one coordinate axis, narrow across the offset and tall along it, for
  which upstream reports containment for ellipses that in fact overlap. The FI
  deviation is the `divisor == 0` branch of `CaseE4NotZero`, where upstream
  writes both symmetric points to `result.points[numPoints]` without
  incrementing between the two writes. Reaching it needs `e2 + e4*y0` to be
  exactly zero at a quartic root, which is not something a random generator ever
  hits, so the case constructs it: with ellipse0 a circle of radius `2^m`
  centred at `K` and ellipse1 concentric with it, with the exact non-unit axes
  `(p,q)` and `(-q,p)` (`p^2+q^2 = 5`) and extents `2^(m+1)`, `2^(m-1)`, the
  standard-form matrices satisfy `M1(1,1) = M0(1,1)` exactly, which makes `e0`,
  `e1` and `e2` all zero; the quartic is then `y^2 * (c4*y^2 + c2)` and `y0 = 0`
  is a double root with `divisor = 0`. `GetStandardForm` divides by
  `|axis[0]|^2`, which is what makes the non-unit axes legitimate and keeps the
  matrices exact.

All other upstream defects recorded for this group's headers are *preserved* by
the port and therefore need no deviation case; the oracle compares them bit for
bit and they agree: `IntrAlignedBox3Sphere3`'s unconditional
`contactPoint += boxCenter` (#250, which is why the no-contact records of the FI
cases report the box centre on both sides), `IntrSegment2AlignedBox2`'s
`cdeParameter = cdeParameter` self-assignment (#203, dropped as a no-op),
`IntrSegment2Segment2`'s `Exact` `line0Parameter` index mix and its
self-contradictory `Result` documentation (#203, #458),
`IntrSegment3Rectangle3`'s three-component `rectCoord` (#141),
`IntrSphere3Triangle3`'s sphere-centre contact point for
`intersectionType = +1` (#203), and
`IntrConvexPolygonHyperplane`'s `NEGATIVE_SIDE_VERTEX` comment (#250, a
documentation-only correction).

## Not covered

Nothing. All 15 headers of group 32 are implemented by the port and every entry
point has a case. Three deliberate restrictions of the *inputs* are worth
recording:

- `IntrSphere3Triangle3`'s arbitrary-precision `operator()` overload (the
  `QFNumber` one) is not ported — PORTING.md ports only the floating-point
  instantiation unless a dependent file needs the exact path, and no file does —
  so it has no case. Upstream's own comment records that its `closest[j]` loop
  runs out of bounds (#203), which is a second reason not to port it.
- The `IntrEllipse2Ellipse2.find.divisorDeviation` case uses non-unit ellipse
  axes. `Hyperellipsoid` documents unit-length axes, but `GetStandardForm`
  divides by `|axis[0]|^2` precisely so that non-unit axes work, and that is
  what makes the exact construction of the `divisor == 0` branch possible. Every
  other ellipse case uses unit-length axes.
- The `IntrConvexPolygonHyperplane` cases whose hyperplane is built from a
  polygon edge use the integer perpendicular of that edge as the normal rather
  than a unit vector, so that the two endpoint heights are exactly zero and the
  `*_EDGE` configurations are reached. Both queries compare only the signs of
  `Dot(N,X) - c`, which a positive scaling of `(N,c)` leaves unchanged.

## Upstream bug suspects

None new for this group's headers. One observation about a *shared* file, for
the group that owns it:

- `RootsPolynomial.h`'s quartic solver is called by `IntrEllipse2Ellipse2`. The
  port (`src/RootsPolynomial.ts`) writes `Math.sqrt(Math.max(alphaSqr, 0))` and
  three sibling expressions where upstream writes `std::sqrt(std::max(alphaSqr,
  (T)0))`. `Math.max(-0, 0)` is `+0` where `std::max(-0, 0)` is `-0`, so the two
  can disagree in the sign of a zero, exactly the class of defect group 31 found
  in `IntrIntervals`. No record of this group's deep run reaches it (the 2310
  non-bit-identical outputs are all of size 1e-12 or larger and come from the
  resolvent cubic), so nothing was changed here; the owner of
  `RootsPolynomial.h` should check it.
