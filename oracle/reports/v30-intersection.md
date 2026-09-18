# Group 30 (`v30-intersection`) — C++ oracle report

61 cases, 20 records each in `oracle/golden/v30-intersection.txt` (1220 records).
1200 records are ordinary comparisons and carry 1684 floating-point outputs;
**all 1684 are bit-identical to the MSVC build of upstream GTE**. The remaining
20 records belong to the single `deviation` case, which is expected to disagree.

Deep run `npm run oracle:deep -- 2000 v30-intersection`: 61 cases, 122000
records, 0 C++ exceptions, 169507 floating-point outputs compared, **100 %
bit-identical**, no disagreement of any kind. The deviation case disagreed on
all 2000 of its records, as intended.

Every case is declared `{ exact: true }`. None of these queries calls the C math
library on the compared path: the arithmetic is `+ - * / sqrt fabs min max` and
comparisons only. The `sin`/`cos` used to build oriented-box and rectangle
frames, and the `sqrt` inside `Normalize`, are applied to *unrecorded* draws;
only the resulting frame or unit vector is recorded as an input, so libm never
enters the compared computation.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise`.

## Coverage

Every computational entry point of every header in the group is covered: each
`TIQuery::operator()` and `FIQuery::operator()` overload, and each protected
`DoQuery` helper (the port exports these as free functions, so they are called
directly on both sides through a small `Expose` subclass in C++).

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `IntrIntervals.h` | `test.finiteFinite`, `test.finiteSemiInfinite`, `test.semiInfiniteSemiInfinite`, `test.dynamic`, `find.finiteFinite`, `find.finiteSemiInfinite`, `find.semiInfiniteSemiInfinite`, `findDynamic.soundBranches`, `findDynamic.leftApproachDeviation` | exact (+1 deviation) | pass |
| `IntrAlignedBox2AlignedBox2.h` | `test`, `find` | exact | pass |
| `IntrAlignedBox2OrientedBox2.h` | `test` | exact | pass |
| `IntrAlignedBox3AlignedBox3.h` | `test`, `find` | exact | pass |
| `IntrAlignedBox3OrientedBox3.h` | `test` | exact | pass |
| `IntrCircle2Circle2.h` | `test`, `find` | exact | pass |
| `IntrHalfspace2Polygon2.h` | `find` | exact | pass |
| `IntrHalfspace3OrientedBox3.h` | `test` | exact | pass |
| `IntrHalfspace3Segment3.h` | `test`, `find` | exact | pass |
| `IntrHalfspace3Sphere3.h` | `test` | exact | pass |
| `IntrHalfspace3Triangle3.h` | `test`, `find` | exact | pass |
| `IntrLine2AlignedBox2.h` | `test`, `find`, `doQuery.ti`, `doQuery.fi` | exact | pass |
| `IntrLine2Line2.h` | `test`, `find` | exact | pass |
| `IntrLine2Triangle2.h` | `test`, `find`, `doQuery.fi` | exact | pass |
| `IntrLine3AlignedBox3.h` | `test`, `find`, `doQuery.ti`, `doQuery.fi` | exact | pass |
| `IntrLine3Rectangle3.h` | `test`, `find`, `test.throughPoint`, `find.throughPoint` | exact | pass |
| `IntrLine3Sphere3.h` | `test`, `find`, `doQuery.fi`, `test.tangent`, `find.tangent` | exact | pass |
| `IntrLine3Triangle3.h` | `test`, `find`, `test.throughPoint`, `find.throughPoint` | exact | pass |
| `IntrOrientedBox2OrientedBox2.h` | `test`, `find` | exact | pass |
| `IntrOrientedBox3OrientedBox3.h` | `test` | exact | pass |
| `IntrRay3Triangle3.h` | `test`, `find`, `test.throughPoint`, `find.throughPoint` | exact | pass |
| `IntrSegment3Triangle3.h` | `test`, `find`, `test.throughPoint`, `find.throughPoint` | exact | pass |

### Generators and the branches they reach

Each generator alternates an exactly representable small-lattice mode with a
uniform mode; several add a constructed mode aimed at a specific branch. The
distributions below are from the 2000-record deep run.

- **Boxes.** Lattice corners give touching and coincident boxes; oriented-box
  frames in lattice mode are signed permutations of the coordinate axes, so the
  parallel-axis short circuit (`existsParallelPair`) and the axis-aligned
  separating axes are reached exactly. `IntrOrientedBox2OrientedBox2.test`
  reports all four `separating` values; the 3D queries report separating
  indices from both the face axes (`separating[1] == -1`) and the nine cross
  products. The `epsilon` argument is drawn from `[-0.25, 0.6]`, which covers
  the clamp of a negative epsilon and a cutoff loose enough to declare nearly
  parallel random frames parallel.
- **Intervals.** Lattice endpoints give touching (`numIntersections == 1`,
  `isPoint`) and coincident intervals; all five `Result::type` values and both
  semi-infinite orientations occur.
- **Circles.** A third mode builds externally tangent, internally tangent and
  identical circles on a lattice, so `find` returns 0, 1, 2 and the
  `numeric_limits<int32_t>::max()` "same circle" sentinel (647/447/279 records
  of the four outcomes).
- **Halfspaces.** In lattice mode the halfspace constant is `Dot(normal, P)`
  for a vertex/centre `P` of the other primitive, so the `distance == 0`
  branches of the clippers are reached exactly. `IntrHalfspace2Polygon2` sees
  clipped polygons of 3 to 8 vertices, the fully-outside case and the
  fully-inside case (upstream's documented empty-polygon-with-`intersect`
  result, finding #139, is reproduced by the port).
- **Lines vs boxes.** Lattice origins and integer directions produce lines
  through box corners and along box faces: `numIntersections` 0, 1 and 2 all
  occur for both the 2D and 3D queries and for both `DoQuery` helpers.
- **Lines, rays and segments vs triangles and rectangles.** Uniform inputs hit
  only a few percent of the time and never land on a boundary, so the
  `*.throughPoint` cases aim the line at
  `v0 + (k1*E1 + k2*E2)/4` with `k1, k2` in `[0, 5]` on an integer triangle.
  All the arithmetic is exact, so the barycentric outputs land on
  0, 1/4, 1/2, 3/4 and 1, and upstream's boundary comparisons `b1 >= 0`,
  `b2 >= 0` and `b1 + b2 <= 1` are evaluated at exact equality. The ray variant
  draws the target parameter from `[-2, 3]` so the `QdN >= 0` test is evaluated
  at exactly zero (ray origin on the triangle); the segment variant places an
  endpoint exactly on the triangle (the `|t| == extent` boundary) and sometimes
  degenerates to a point, which upstream treats as parallel. The rectangle
  variant aims at `C + (a/2)*e0*W0 + (b/2)*e1*W1` with `a, b` in `[-3, 3]`, so
  the `|W1dDxQ| == extent[0]*|DdN|` edge and corner cases are exact.
- **Lines vs spheres.** `IntrLine3Sphere3.*.tangent` puts the line origin at
  `C + (r+d)*e_i + s*e_j` with an axis direction `e_j`, which makes the
  discriminant exactly `r^2 - (r+d)^2`: negative, exactly zero or positive with
  no rounding. The tangent branch (`numIntersections == 1`) is otherwise
  unreachable and is now covered in the committed goldens.

## Port defects fixed

None. No disagreement was found in 122000 deep-run records, so no `src/` file
was modified by this group.

## Deliberate deviations demonstrated

| case | records disagreeing | record of the decision |
| --- | --- | --- |
| `IntrIntervals.findDynamic.leftApproachDeviation` | 20 of 20 golden, 2000 of 2000 deep | `docs/UPSTREAM-FINDINGS.md`, `IntrIntervals.h` dynamic FIQuery; issue [#62](https://github.com/gradientspaceai/gtengine-js/issues/62) |

When `interval0` lies strictly to the left of `interval1` and the intervals
approach each other, upstream reports `interval0[0] + firstTime * speed0` as the
contact point, which is the moved *left* endpoint of `interval0`; the contact
point is the moved *right* endpoint, `interval0[1] + firstTime * speed0`, which
is also what upstream itself computes in the mirrored "interval0 on the right"
branch. The port keeps the fix. The case constructs a non-degenerate
`interval0` strictly left of `interval1` with `speed0 > speed1`, so upstream and
the port differ on `overlap[0]` and `overlap[1]` on every record — the C++ build
confirms the defect is real and still present upstream. The main dynamic case
(`findDynamic.soundBranches`) swaps the intervals whenever that branch would be
taken, so it compares exactly.

All other upstream defects recorded for this group's headers are *preserved* by
the port (`IntrHalfspace2Polygon2` empty polygon #139, `IntrHalfspace3Segment3`
case (1,1,0) #139, `IntrOrientedBox2OrientedBox2` TI/FI touching disagreement
#141, `IntrLine3Rectangle3::Result::rectCoord[2]` #141,
`IntrAlignedBox3OrientedBox3` vs `IntrOrientedBox3OrientedBox3` cutoff
comparison #450, `IntrCircle2Circle2` TI/FI disk-vs-curve #450,
`IntrIntervals` dead store #62), so they need no deviation case: the oracle
compares them bit for bit and they agree.

## Not covered

Nothing. All 22 headers of group 30 are implemented by the port and every
entry point has a case.

`IntrLine3Rectangle3` `Result::rectCoord[2]` and the `FIQuery` result fields of
the aligned-box/aligned-box queries on the non-intersecting path are emitted
because upstream's `Result` constructor defines them on every path (0 and the
default `AlignedBox` `[-1, 1]` respectively), so they are part of the compared
contract.

## Upstream bug suspects

None new. Everything observed in this group is already in
`docs/UPSTREAM-FINDINGS.md`.
