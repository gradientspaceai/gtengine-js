# v08-compgeom (verify group 8, computational geometry)

41 cases: 35 exact, 2 with a measured tolerance, 4 deliberate deviations.
`oracle/cpp/cases/v08-compgeom.cpp` and
`test/oracle/v08-compgeom.oracle.test.ts`.

Almost everything in this group is `+ - * / sqrt` and comparisons on doubles,
so the expectation is bit-identity. The three places where the C math library
enters are pinned rather than tolerated:

* `SortPointsOnCircle::ByAngle` compares `std::atan2` values, and the result is
  a permutation, which no tolerance can rescue. The generator rejects point
  pairs whose directions from the sort centre are distinct but nearly
  parallel: with `|DotPerp(W0,W1)| >= 1e-6 * |W0| * |W1|` the angular gap is
  ten orders of magnitude above an ulp of `atan2`. Exactly coincident `W`
  vectors are allowed, because `atan2` of identical arguments is identical
  within one library and the comparator then falls through to its arithmetic
  tie-break.
* `InscribedFixedAspectRectInQuad::Execute` uses `std::atan2` only to pick the
  quadrant index `j = floor((2/pi) * angle)` of each inner edge normal; the
  angle never reaches the result. The lattice quads of the generator have no
  axis-parallel edge, so every normal is at least `atan(1/12)` from a quadrant
  boundary and `j` is decided by the geometry.
* `MinimumVolumeBox3` clouds are built by rotating lattice clouds with
  `sin`/`cos`, but only the final rotated coordinates are recorded as inputs.

## Coverage

| header | cases | comparison | deep run (2000 records/case) |
| --- | --- | --- | --- |
| `DisjointIntervals.h` | `insertRemove`, `operators` | exact | pass |
| `DisjointRectangles.h` | `insertRemove`, `operators` | exact | pass |
| `SortPointsOnCircle.h` | `byAngleAndByGeometry`, `byGeometry.ties` | exact | pass |
| `CircleThroughTwoPointsSpecifiedRadius.h` | `compute` | exact | pass |
| `CircleThroughPointSpecifiedTangentAndRadius.h` | `compute` | exact | pass |
| `ConvexHullSimplePolygon.h` | `compute` | exact | pass |
| `ConvexPolyhedron3.h` | `construct` | exact | pass |
| `PrimalQuery2.h` | `toLine`, `toLineWithOrder`, `toTriangle`, `toCircumcircle`, `toLineExtended` | exact | pass |
| `PrimalQuery3.h` | `toPlane`, `toTetrahedron`, `toCircumsphere` | exact | pass |
| `ExtremalQuery3.h`, `ExtremalQuery3PRJ.h` | `getExtremeVertices` | exact | pass |
| `NearestNeighborQuery.h` | `build`, `build.tiedCoordinates`, `build.maxLevelAssert`, `findNeighbors.all`, `findNeighbors.limited` | exact | pass |
| `BoxManager.h` | `initializeAndUpdate` | exact | pass |
| `InscribedFixedAspectRectInQuad.h` | `execute` | exact | pass |
| `TriangulateEC.h` | `triangulate`, `triangulatePolygon`, `triangulateWithHole`, `triangulateWithHoles`, `triangulateTree` | exact | pass |
| `MinimumVolumeBox3.h` (FloatingPoint) | `compute`, `compute.lowDimension`, `computeHull`, + 2 deviations | `computeHull` and `lowDimension` exact, `compute` volume 2e-2 | pass |
| `MinimumVolumeBox3.h` (Rational) | `compute`, `compute.lowDimension`, `computeHull`, + 2 deviations | `computeHull` and `lowDimension` exact, `compute` volume 1.5e-1 | pass |

The two `compute` cases (the arbitrary-point-cloud query) emit the hull
dimension, the minimum volume and the reference booleans, not the oriented box:
see "Not covered". The two `computeHull` cases (the vertices-and-indices query,
where the mesh is fixed) emit the full box and are bit-identical over 2000
records each.

The deep run (`npm run oracle:deep -- 2000 v08-compgeom`) passes: 41 cases,
2000 records each, generation 203 s and replay 60 s on this machine.

Every computational public entry point the port implements is covered. Both
overloads of every `PrimalQuery2` / `PrimalQuery3` predicate (the vertex-index
one and the `Vector` one) are exercised in the same record; the four Boolean
operators of `DisjointIntervals` and `DisjointRectangles` are compared against
a brute-force membership oracle as well as against C++; `BoxManager` is driven
through `Initialize`, `SetBox`, `Update`, `GetBox`, `GetOverlap` and a second
`Initialize`; `TriangulateEC` is exercised through all five queries.

### Generators

Each case mixes a small-lattice mode, a uniform mode and constructions aimed at
particular branches:

* `DisjointIntervals` / `DisjointRectangles` use a `[-3,3]` lattice, so
  touching, nested, coincident and empty (`xmin >= xmax`) requests are
  constant; the interval and rectangle sets are also combined with the empty
  set.
* `PrimalQuery2` / `PrimalQuery3` use the 12 lattice points of the circle
  `x^2+y^2 = 25` and the 30 lattice points of the sphere `x^2+y^2+z^2 = 9`, so
  cocircular and cospherical quadruples occur constantly, plus collinear and
  coplanar runs `V0 + a*dir0 + b*dir1`. With coordinates bounded by 10, every
  product and sum inside the determinants is an integer below `2^53`, so the
  predicates are evaluated at exact zero wherever the configuration is
  degenerate.
* `ConvexHullSimplePolygon` and `TriangulateEC` walk the boundary lattice
  points of an axis-aligned rectangle, which puts long collinear runs on every
  edge so that `WhichSide` returns 0 and the `>= 0` / `<= 0` branches are both
  taken.
* `ExtremalQuery3PRJ` uses lattice directions that are frequently perpendicular
  to a face or parallel to an edge of the tetrahedron, cube and octahedron, so
  the strict comparisons see exact ties.
* `NearestNeighborQuery.findNeighbors.all` allows distance ties on purpose (the
  neighbour limit is at least the number of sites, so the retained set cannot
  depend on the visiting order) and the heap pops them by decreasing site
  index.

### Independent reference checks

* **PrimalQuery, exact BigInt predicates.** `oracle/reports` does not carry the
  script, but the deep-run outputs of `PrimalQuery2.toLine`,
  `PrimalQuery2.toCircumcircle`, `PrimalQuery3.toPlane` and
  `PrimalQuery3.toCircumsphere` were recomputed with exact `BigInt`
  determinants on every record whose coordinates are all integers. On the
  committed goldens: 488 predicates checked, 0 disagreements; on the deep run:
  48280 predicates checked, 0 disagreements. Both the C++ build and the port
  therefore get the exact sign on the lattice inputs, including the exactly
  cocircular and cospherical configurations.
* **Interval and rectangle set operations vs brute force.** Every record of
  `DisjointIntervals.operators` and `DisjointRectangles.operators` emits a
  boolean that compares membership in the computed result with the Boolean
  combination of the memberships in the operands, sampled on a half-integer
  lattice covering the whole support, plus a "the reported intervals are
  nonempty, ascending and separated" boolean. All true on both sides.
* **Box overlaps vs brute force.** `BoxManager` emits a boolean comparing the
  reported overlap set with the pairwise brute-force test after every
  `Initialize` and every `Update`. All true on both sides, so the incremental
  update never drifts from the full sweep.
* **Hull validity.** `ConvexHullSimplePolygon` emits a boolean that the hull is
  convex, counterclockwise, of positive area and contains every polygon vertex.
* **Triangulation validity.** Every `TriangulateEC` case emits a boolean that
  every produced triangle is counterclockwise or degenerate and that the total
  signed area equals the signed area of the region (outer polygon minus holes),
  computed with a rational shoelace sum that is exact on the lattice inputs.
  The simple-polygon cases additionally check the triangle count `n - 2`.
* **Nearest neighbours vs brute force.** Both `findNeighbors` cases emit a
  boolean comparing the returned index set with a brute-force scan.
* **Extreme vertices vs brute force.** `ExtremalQuery3PRJ` emits a boolean that
  the returned vertices realize the maximum and the minimum of
  `Dot(direction, V)` over the polytope vertices, computed without the centroid
  shift the query uses.
* **Minimum-volume box containment.** Every `MinimumVolumeBox3` case emits a
  boolean that the returned box contains the input points, and the
  floating-point case also checks that the reported volume is the product of
  the box side lengths.
* **Circle constructions.** Both circle cases emit a boolean that each returned
  circle passes through the given points and, for the tangent construction,
  that the line is tangent to it.
* **Inscribed rectangle.** The case emits a boolean that all four rectangle
  corners lie on the inner side of all four quad edges.

All of these booleans are `true` on both sides for every record, at both the
committed and the deep size, with one deliberate exception: a handful of
`MinimumVolumeBox3FloatingPoint` records where upstream's own box does not
contain the points, which is what the `nonContaining` deviation case is for.

## Port defects fixed

**`src/MinimumVolumeBox3FloatingPoint.ts`, `computeVolume`: the deliberate fix
of upstream issue #405 was not confined.** Upstream's `ComputeVolume` assumes
that a vertex of the hull edge realizes the minimum along `axis[0]` and
`axis[1]`:

```cpp
candidate.minSupportIndex[0] = mEdges[candidate.edgeIndex[0]].v[0];
pmin[0] = Dot(candidate.axis[0], mTVertices[candidate.minSupportIndex[0]]);
```

The port replaced both lines by a `getExtreme` hill climb unconditionally. The
hill climb returns a vertex with the same or a smaller *double* projection, but
not necessarily the same vertex: among vertices whose double projections are
equal it can pick another one. `minSupportIndex` is then fed to the exact
rational `getMinimumVolumeBox`, where two vertices with equal double
projections can have different *exact* projections, so the box centre and the
extents came out an ulp away from upstream's on ordinary inputs. Error size: up
to 6e-16 relative on the centre and the extents, and occasionally a different
winning candidate and a 0.8 percent different volume.

The fix now evaluates upstream's expression first and only replaces the index
and the projection when the hill climb finds a **strictly smaller** projection,
which is exactly the condition under which upstream's assumption is false.
Regression tests in `test/MinimumVolumeBox3FloatingPoint.test.ts`:

* a probe subclass verifies, over more than a thousand candidates, that the
  minimum support index is upstream's edge vertex on every candidate where that
  vertex really realizes the minimum (0 failures), while the defect itself
  occurs among those candidates;
* a second test runs the port beside a subclass that restores upstream's
  verbatim `ComputeVolume` and requires that the two agree bit for bit whenever
  no candidate violated the assumption, that they differ somewhere (so the test
  fails if the fix is removed), and that upstream's box fails containment on
  some of the 300 sampled clouds.

No other port defect was found in this group: every other case was bit-identical
to the MSVC build on the first run.

## Deliberate deviations demonstrated

| case | finding | what deviates |
| --- | --- | --- |
| `MinimumVolumeBox3FloatingPoint.compute.coplanar` | #352 | the dimension-2 Newell normal loop `for (i0 = numHull - 1, i1 = 1; ...)` drops the wrap-around cross-product term, so the plane normal is wrong for every coplanar point set and exactly zero when the hull is a triangle; the port starts at `i1 = 0` |
| `MinimumVolumeBox3Rational.compute.coplanar` | #355 | the same omission in the rational file |
| `MinimumVolumeBox3FloatingPoint.compute.nonContaining` | #405, #426 | clouds on which upstream's box does not contain the input points, the shared observable of the `ComputeVolume` support assumption and of the `GetExtreme` plateau stall; the generator keeps the worst offender it finds in 48 capped attempts |
| `MinimumVolumeBox3Rational.compute.variableT` | #355 | clouds on which upstream reaches `MinimizerVariableT`, whose `tminNumer`, `tmaxNumer` and `tDenom` are declared `T const&` instead of `Number const&`, so the exact `BSNumber` expressions its callers pass are rounded to `double` and the whole `t`-variable level curve is sampled at rounded parameters; the port samples them exactly. A subclass in the case file counts the calls, which is the exact separator |

Deviating records in the deep run: 1167 / 2000 (`FloatingPoint.compute.coplanar`),
1433 / 2000 (`FloatingPoint.compute.nonContaining`), 1185 / 2000
(`Rational.compute.coplanar`), 321 / 2000 (`Rational.compute.variableT`).

`InscribedFixedAspectRectInQuad`'s two failures (#395: the `alpha` assertion on
convex quads whose normals share a quadrant, and the untoleranced degenerate
feasible interval) are *preserved* in the port, so they are compared as throw
parity inside the ordinary case rather than as deviations; the deep run
contains such records and the port throws on exactly the same ones.
`PrimalQuery2`'s four-argument `ToLine` collinear-order defect (#100) and
`SortPointsOnCircle`'s non-strict-weak-order comparator (#394) are likewise
preserved and compared bit for bit.

## Not covered

* **`MinimumVolumeBox3`, the oriented box of the point-cloud query.** See the
  upstream suspect below. Where the mesh is fixed (`computeHull`) the box is
  emitted in a form invariant under the axis signs and the axis order, the
  acceptance test "upstream's own answer does not change when the same mesh is
  presented in any other cyclic order or reversed" is recorded as an input so
  that the replay takes the same branch, and the comparison is then exact over
  2000 records. For an arbitrary point cloud (`compute`) the hull itself comes
  out of `ConvexHull3`, whose vertex ordering is a second container-dependent
  input to the same search, and a probe over all cyclic reorderings of the
  input points still left a few records in 2000 on which upstream and the port
  found different (both valid, both containing) boxes. Those two cases
  therefore compare the hull dimension and the reference booleans exactly and
  the minimum volume with a measured tolerance: the floating-point case is
  bit-identical on 1997 of 2000 records with a largest scaled difference of
  8.3e-3 (tolerance 2e-2), the rational case on 1994 of 2000 with a largest
  scaled difference of 9.3e-2 (tolerance 1.5e-1).
* **`MinimumVolumeBox3Rational`, `MinimizerVariableT`.** Upstream declares its
  `tminNumer`, `tmaxNumer` and `tDenom` parameters as `T const&` where the
  sibling `MinimizerVariableS` uses `Number const&`, so every `t`-variable level
  curve is sampled at parameters rounded to `double`; the port samples them
  exactly (#355). The rational `compute` case therefore rejects clouds on which
  `MinimizerVariableT` is reached at all (a subclass in the case file counts the
  calls), and the `compute.variableT` deviation case aims at them.
* **`MinimumVolumeBox3Rational` on lattice clouds.** Its arithmetic is exact, so
  many candidates of a lattice cloud tie exactly and the winner is decided by
  the enumeration order. The case draws uniform clouds only.
* **`NearestNeighborQuery` with tied split coordinates.** `std::nth_element`
  only guarantees the partition postconditions; which of the elements equal to
  the median end up on each side is unspecified, and the port uses its own
  quickselect. Measured: emitting the leaf contents makes 11 of 20 records
  disagree, and so does the split value of any internal node below the root,
  because that node's subrange then holds a different multiset. The
  `build.tiedCoordinates` case compares only what the site counts determine
  (depth, largest leaf size, node count, the node ranges, and the root split,
  which is an order statistic over the whole array). The main cases give every
  axis pairwise distinct coordinates, which makes the whole tree determined.
* **`NearestNeighborQuery.findNeighbors` with a neighbour limit below the
  candidate count and equidistant candidates.** Upstream replaces the heap top
  only on a strict improvement, so the retained set depends on the visiting
  order, hence on the `nth_element` permutation. The `findNeighbors.limited`
  generator rejects query points with equidistant sites (capped, falling back to
  the point with the fewest ties seen); `findNeighbors.all` covers distance ties
  with a limit that retains every candidate.
* **`SortPointsOnCircle` with a point coincident with the sort centre.**
  `LessThanByGeometry` is then not a strict weak ordering, which is undefined
  behaviour for `std::sort` (#394). Such inputs are kept out of every
  generator.
* **`ExtremalQuery3BSP.h`** is not in this group (it belongs to group 82's
  file list only through `ExtremalQuery3PRJ`); it is not exercised here.
* **`MinimumVolumeBox3GPU`** is not ported and has no case.
* **`ConvexPolyhedron3`'s `alignedBox` when `GenerateAlignedBox` is not
  called** is the `AlignedBox3` default (`min = -1`, `max = +1`) on both sides
  and is compared as such.

## Upstream bug suspects

**1. `MinimumVolumeBox3FloatingPoint.h` / `MinimumVolumeBox3Rational.h`: the
result depends on a `std::unordered_map` iteration order.** `ExtractMeshTopology`
numbers the mesh edges and triangles by iterating `ETManifoldMesh`'s `mEMap`
and `mTMap`, which are `std::unordered_map`. `mEdgeIndices` is then every
ordered pair `(e0, e1)` with `e0 < e1` **in that numbering**, and
`ProcessEdgePair` is not symmetric in its two edges: the first supplies the
`N` normals and the second the `M` normals, and the level-curve branch and the
sampled axes differ between the two roles. Both the set of candidates examined
and the tie-break ("the first candidate with a strictly smaller volume wins")
therefore depend on a hash-container iteration order, which no specification
fixes and which differs between standard libraries.

Measured with the probe in the case file, which runs upstream's own query twice
on the same mesh with the triangles presented in a different order:

* **Rational pipeline, fixed mesh.** Presenting the same tetrahedron, cube or
  octahedron with its triangles rotated changes the reported box on
  **653 of 2000** generated meshes (33 percent). The rational pipeline compares
  candidate volumes exactly, so this is not rounding: a different one of
  several exactly-tied candidates is reported.
* **Floating-point pipeline, fixed mesh.** The same probe (widened to all
  cyclic rotations and the reversal) plus the containment check rejects
  **1168 of 2000** meshes. That figure mixes the order dependence with the
  #405/#426 containment failures and is not separated here.
* **Point clouds.** On small random clouds the reported *minimum volume*
  changes by up to 0.8 percent under a reordering of the input points, and the
  rational pipeline on lattice clouds was affected on roughly 9 percent of the
  clouds drawn before the generator was restricted to uniform ones.

This is not a rounding effect: the rational pipeline compares candidate volumes
exactly, so what changes is which of several exactly-tied candidates is
reported. The port cannot reproduce MSVC's bucket order and should not try.
Consequence for a user: the oriented box returned by `MinimumVolumeBox3` is not
reproducible across standard library implementations, and neither is the
reported minimum volume when the sampled candidates tie. A fix would be to
number the edges and triangles by a deterministic key order (the mesh already
has `EdgeKey`/`TriangleKey` comparators) instead of by hash order, and to make
`ProcessEdgePair` symmetric or to enumerate both orders of every pair.

**2. `InscribedFixedAspectRectInQuad.h` does not compile on its own.** It uses
`GTE_C_TWO_PI` and `GTE_C_INV_HALF_PI` but includes only `IntrIntervals.h`,
`Vector2.h`, `<array>`, `<cmath>`, `<cstddef>` and `<utility>`; `Constants.h`
is missing. Including the header first in a translation unit is a hard compile
error (`C2065: 'GTE_C_TWO_PI': undeclared identifier`). The case file works
around it by including `Mathematics/Constants.h` first. Minor, but it is a
defect of the header's self-containment.

No other new upstream defect was found: everything else in the group agreed
bit for bit with the MSVC build and with the independent references.
