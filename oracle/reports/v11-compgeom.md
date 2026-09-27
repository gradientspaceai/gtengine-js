# Verify group 11 (computational geometry) against the MSVC build of upstream GTE

Family `v11-compgeom`, 23 cases, 20 golden records each. Deep run
`npm run oracle:deep -- 2000 v11-compgeom` (46000 records, 8000 of them throw
records): **every case passes** (24 tests including the coverage check,
replay 5.3 s, generation and replay together under 10 s).

The group is `ConstrainedDelaunay2.h`, `ConvexHull3.h`, `Delaunay2Mesh.h`,
`Delaunay3Mesh.h`, `MinimumAreaBox2.h`, `MinimumWidthPoints2.h` and
`SeparatePoints2.h`. The C++ side instantiates exactly what the port
implements:

| upstream instantiation | why |
| --- | --- |
| `ConstrainedDelaunay2<double>` | the replacement class; the deprecated `<InputType, ComputeType>` specialization is not ported |
| `ConvexHull3<double>`, `lgNumThreads = 0` | the port has only the single-threaded path |
| `Delaunay2Mesh<double>`, `Delaunay3Mesh<double>` | the replacement classes |
| `MinimumAreaBox2<double, BSRational<UIntegerAP32>>` | the port is exact-only by design; upstream's own header says a correct result is guaranteed only for an exact ComputeType and names this one. The floating-point instantiation is compared separately, as a design deviation |
| `MinimumWidthPoints2<double>` | `Rational = BSRational<UIntegerAP32>` inside, as the port |
| `SeparatePoints2<double, double>` | `ComputeType` is vestigial |

No compared path calls the C math library: outputs are combinatorial
(dimensions, index lists, adjacencies, partitions, predicate results) or are
computed with `+ - * / sqrt` on doubles, or exactly on rationals that are
rounded to double once. **All 23 cases are declared `{ exact: true }`; every
real output of every non-deviation case is bit-identical to the MSVC value**
(259664 real outputs in the deep run, 259664 exact; no tolerance anywhere in
the family).

**No port defect was found.** This branch contains no `src/` change.

## Canonicalization

`ETManifoldMesh`, `VETManifoldMesh` and `TSManifoldMesh` store features in
`std::unordered_map`, and `ConstrainedDelaunay2::mInsertedEdges` is a
`std::unordered_set`; the port enumerates them in sorted-key order. As in
v09, triangle and tetrahedron lists are sorted by their stored vertex tuple
(`Insert` keeps the caller's rotation, so the tuple is a total key) with the
adjacencies remapped through the same permutation; `ConvexHull3`'s dimension-3
vertex list and face triples, `GetHull()` edge pairs and the inserted-edge set
are sorted. Ordered results (`ConvexHull3`'s dimension 0-2 lists, CDT
partitions, `MinimumAreaBox2::GetHull`, support indices) are emitted verbatim.
The containment queries of `Delaunay2Mesh` / `Delaunay3Mesh` start the walk at
"feature 0", which is a different feature on the two sides, so their query
points are aimed at configurations whose answer does not depend on the start
(strictly inside the canonically chosen simplex, far outside, or uniform);
v09's `getContainingTriangle` cases pin the start and cover edges and
vertices.

`ConstrainedDelaunay2::GetLinkEdges` iterates `Vertex::TAdjacent`, a
pointer-keyed set. At most two link triangles pass the
`sign0 >= 0 && sign1 <= 0` test, and only when the constraint runs through
their shared vertex, where both call `ProcessCoincidentEdge` with that vertex,
so the result does not depend on that order. The strip is removed from and
re-inserted into the graph in an order that affects only hash-table layout.
`ConvexHull3`'s `std::set<Triangle*> visited` is used for membership only, and
the visible region is a connected set independent of the seed.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `ConstrainedDelaunay2.h` | `ConstrainedDelaunay2.insert` (both `operator()` overloads, alternating; one to four `Insert` calls per record with uniform, already-present, farthest-vertex and chained constraints, endpoints given by duplicate indices included; `GetInsertedEdges`; then `UpdateIndicesAdjacencies`, `GetNumTriangles`, `GetIndices`, `GetAdjacencies`, `GetHull`, `GetGraph`) | exact | pass, 2000 records |
| | `ConstrainedDelaunay2.insert.strip` (a long constraint through a band of points, with vertices on the constraint and a collinear chain; optionally a second constraint crossing the first, then the first again) | exact | pass, 2000 records |
| | `ConstrainedDelaunay2.insert.invalidEdgeThrows` (throw parity for "Invalid edge." after the duplicate substitution) | exact | pass, 2000 throw records |
| | `ConstrainedDelaunay2.deviation.staleInsertedEdges` (#325), `ConstrainedDelaunay2.deviation.duplicatesRead` (#325) | deviation | pass, 2000/2000 deviate each |
| `ConvexHull3.h` | `ConvexHull3.compute` (`operator()(points, 0)`, `GetDimension`, `GetVertices`, `GetHull`, `GetHullMesh`; dimension 0 / 1 / 2 / 3 reached on 464 / 412 / 347 / 777 records) | exact | pass, 2000 records |
| `Delaunay2Mesh.h` | `Delaunay2Mesh.query` (`GetNumVertices`, `GetNumTriangles`, `GetVertices`, `GetIndices`, `GetAdjacencies`, `GetInvalidIndex`, `GetContainingTriangle`, the per-triangle `GetVertices` / `GetIndices` / `GetAdjacencies` / `GetBarycentrics` in range and out of range) | exact | pass, 2000 records |
| | `Delaunay2Mesh.constructorThrows` (the dimension assert) | exact | pass, 2000 throw records |
| `Delaunay3Mesh.h` | `Delaunay3Mesh.query`, `Delaunay3Mesh.constructorThrows` (the same surface in 3D) | exact | pass, 2x2000 records |
| `MinimumAreaBox2.h` | `MinimumAreaBox2.compute` (overloads 1/2 with both `useRotatingCalipers` values, so `ComputeBoxForEdgeOrderN` and `...OrderNSqr`; box, `GetArea`, `GetSupportIndices`, `GetHull`, `GetNumPoints`, `GetPoints`; hull dimension 0 and 2, square ties) | exact | pass, 2000 records |
| | `MinimumAreaBox2.compute.dimension1` (hull dimension 1, everything but the extreme indices) | exact | pass, 2000 records |
| | `MinimumAreaBox2.computeConvexPolygon` (overload 3: the points as polygon, an index subset, exactly collinear polygon vertices that `RemoveCollinearPoints` removes, the two early returns) | exact | pass, 2000 records |
| | `...deviation.dimension1Extremes`, `...deviation.staleState`, `...deviation.polygonPoints`, `...deviation.removeCollinear` (#328, #402, #286 pattern) and `...deviation.floatComputeType` (design) | deviation | pass; 1427, 2000, 2000, 1966, 1939 of 2000 deviate |
| `MinimumWidthPoints2.h` | `MinimumWidthPoints2.compute` (overloads 1/2, both `useRotatingCalipers` values: the exact antipode search and the floating-point edge search; hull dimension 0, 1, 2) | exact | pass, 2000 records |
| | `MinimumWidthPoints2.computeIndexed` (overload 4 and through it 3: arbitrary index lists with duplicates) | exact | pass, 2000 records |
| | `MinimumWidthPoints2.invalidInputThrows` (both `LogAssert`s) | exact | pass, 2000 throw records |
| `SeparatePoints2.h` | `SeparatePoints2.compute` (lattice, uniform, regular polygons, degenerate sets; 1154 separated / 265 overlapping / 581 degenerate) | exact | pass, 2000 records |
| | `SeparatePoints2.deviation.roundoff` (#328) | deviation | pass, 2000/2000 deviate |

### Generators and the branches they reach

Point sets have at most 11 points. The shared 2D and 3D modes are v09's:
uniform doubles, lattice `[-3,3]`, cocircular / cospherical lattice subsets
with their center, and dense lattice `[-2,2]` (duplicates and collinear or
coplanar subsets). Every triangulating case draws from v09's `Sound2` /
`Sound3` probe, the exact separator for the Delaunay2/3 findings #391 and
#283, inside a 64-attempt rejection loop that redraws every coordinate and
falls back to a parabola / moment-curve set. Every loop in the file is capped.

`ConstrainedDelaunay2`'s branches, counted on the deep-run inputs of the two
main cases by instrumenting the port (which agrees with the build on every
record):

| branch | calls |
| --- | --- |
| `ProcessTriangleStrip`, edge consumed | 3133 |
| `ProcessTriangleStrip`, vertex on the constraint (`querySign == 0`, subdivided) | 508 |
| `ProcessCoincidentEdge` | 4455 |
| constraint already present when inserted | 3634 + 64 |
| constraint endpoint given by a non-representative duplicate index | 390 |
| `SelectSplit` with more than one candidate (the exact pseudo-distance minimum) | 4056 |
| `SelectSplit` with a single candidate | 7605 |
| records with two properly crossing constraints | 224 + 733 |

`ConvexHull3` adds exactly collinear, exactly coplanar (a lattice plane with a
non-axis-aligned normal, the projection path) and all-equal modes.
`MinimumAreaBox2` / `MinimumWidthPoints2` add an axis-aligned or 45-degree
lattice square with boundary and interior points, where the minimum is
attained four times and the first-minimum tie break decides the answer.
`SeparatePoints2` draws two point sets per record; its main case accepts a
draw only when every floating-point side classification upstream performs
agrees with the exact one (the exact superset of the inputs where upstream's
control flow equals the port's), and sets of one or two points are drawn only
in the degenerate mode.

## Port defects fixed

None. Every association order, division form, seed and branch in
`src/ConstrainedDelaunay2.ts`, `src/ConvexHull3.ts`, `src/Delaunay2Mesh.ts`,
`src/Delaunay3Mesh.ts`, `src/MinimumAreaBox2.ts`, `src/MinimumWidthPoints2.ts`
and `src/SeparatePoints2.ts` already matches upstream bit for bit on all 46000
deep-run records, so this branch carries no `src/` change.

## Deliberate deviations demonstrated

1. **`ConstrainedDelaunay2.deviation.staleInsertedEdges` (#325, finding 1),
   2000/2000.** `operator()` forwards to `Delaunay2<T>::operator()`, which
   never touches `mInsertedEdges`. The record triangulates A, inserts one to
   three constraints, triangulates B with the same functor (the pointer
   overload), and emits `GetInsertedEdges()`: upstream still reports A's
   constraints, whose keys index unrelated vertices of B; the port reports
   none. A constraint then inserted into B gets the same partition and the
   same triangulation on both sides; only the edge set differs.
2. **`ConstrainedDelaunay2.deviation.duplicatesRead` (#325, finding 2),
   2000/2000.** `Insert` evaluates `duplicates[edge[i]]` before its range
   assert. The out-of-range read is made deterministic: the functor first
   triangulates a set A of `nB + m` points whose point `j >= nB` repeats a
   point `r < nB` (so A leaves `mDuplicates[j] = r`), then triangulates the
   prefix B of `nB` points; `mDuplicates.clear(); resize(nB)` keeps the
   capacity and the stale element. `Insert({j, s})` on B reads `r`, passes the
   assert and inserts the unrelated constraint `<r, s>`; upstream returns a
   partition and a modified triangulation, the port throws "Invalid edge."
   The generator accepts only draws where `r` is a representative vertex of B
   and `s` a different one, so every record is a silent wrong result.
3. **`MinimumAreaBox2.deviation.dimension1Extremes` (#328), 1427/2000.** The
   residue is exactly explained: over the 2000 deep-run records the two sides
   agree precisely on the 573 whose point 0 coincides with the line origin
   `points[hull[0]]` (checked record by record; no exception), which is the
   only situation in which upstream's `imin = imax = 0` seed names the right
   point.
4. **`MinimumAreaBox2.deviation.staleState` (#328), 2000/2000.**
5. **`MinimumAreaBox2.deviation.polygonPoints` (#402), 2000/2000.**
6. **`MinimumAreaBox2.deviation.removeCollinear` (#286 pattern, #328),
   1966/2000.** The root cause was confirmed record by record: on all 2000
   deep-run records upstream's output (box, area, support indices) is
   bit-identical to the port run on the polygon with the duplicated corner
   removed *entirely* (both copies: `P[d]` is dropped because the next edge
   has length zero, `P[d+1]` because the previous one has). The 34 agreeing
   records are those where removing that corner changes no output.
7. **`MinimumAreaBox2.deviation.floatComputeType` (design, port note),
   1939/2000.** Upstream's `MinimumAreaBox2<double, double>` (brute-force
   default) against the port's exact computation on uniform points. Most
   records differ in the last bits of the box. One committed record differs
   in the box itself with the area equal to 1e-16: the hull is a triangle
   whose angles are all acute, where every edge gives a bounding rectangle
   of exactly twice the triangle's area, so the minimum is an exact
   three-way tie; the exact port takes the first edge, upstream's rounded
   areas pick another. That is the input class for which upstream's header
   recommends an exact `ComputeType`.
8. **`SeparatePoints2.deviation.roundoff` (#328), 2000/2000.** Every record
   is an exactly overlapping pair of regular polygons that upstream reports
   as separated.

The confinement of every fix is checked by the broad main cases, not only
asserted: the main `ConstrainedDelaunay2`, `MinimumAreaBox2` and
`SeparatePoints2` cases agree bit for bit on all their deep-run records.

## Independent reference checks (agreement is not correctness)

The deep-run inputs of the main cases were re-checked against references that
share no code with either implementation (exact predicates on the port's
`BSNumber`, brute force in double). Because the port agrees with the build on
every one of these records, each check covers both sides.

* **ConstrainedDelaunay2** (4000 records): every output triangle is exactly
  counterclockwise; the exact sum of the triangle areas equals that of the
  unconstrained Delaunay triangulation of the same points (same region, no
  gap, no overlap); every sub-edge of every partition is a graph edge when no
  two constraints of the record properly cross, and every sub-edge of the last
  constraint is one when they do. No violation.
* **ConstrainedDelaunay2 finding 3 on the real build.** An exact in-circle test
  across every unconstrained interior edge finds a locally non-Delaunay edge
  in 113 of 2000 `insert` records and 314 of 2000 `insert.strip` records, on
  both sides identically: the strip fill is a valid constrained triangulation
  that is not constrained-Delaunay, as recorded (preserved).
* **ConvexHull3** (777 dimension-3 records): no input point is strictly
  outside any hull face by the exact orientation predicate, the hull mesh
  satisfies `V - E + T = 2`, and no point that is a strict coordinate extreme
  is missing from `GetVertices()`.
* **MinimumAreaBox2** (1339 dimension-2 records): the box contains every
  point (to 1e-12 relative, the box being rounded from rationals), and its
  area is within 8.9e-16 relative of a brute-force minimum over all hull
  edges in double, never above it.
* **MinimumWidthPoints2** (1640 dimension-2 records): the box contains every
  point and `2 * min(extent)` is within 5.6e-16 relative of a brute-force
  minimum over hull edges of the maximum distance to the edge's line.
* **SeparatePoints2** (2000 records): on all 1154 `true` records the returned
  line has all of one set on one closed side and all of the other on the
  other side, up to a residual below 1.6e-11 (the line is rounded); every
  one of the 265 nondegenerate
  `false` records has exactly overlapping hulls (a vertex strictly inside the
  other hull, two properly crossing edges, or a centroid strictly inside);
  the other 581 `false` records have a hull of dimension below 2.
* **Delaunay2Mesh / Delaunay3Mesh** (2 x 2000 records): every returned
  containing simplex has barycentrics >= 0 for the query; the barycentrics of
  the drawn simplex sum to 1 within 1.2e-13 and reconstruct the query within
  3.1e-13 relative (queries go out to 40 units from data of size 4).

## Sensitivity of the cases

Measured by mutating `src/` and replaying the deep run:

* `MinimumAreaBox2`'s axis `Normalize` written as a componentwise division
  instead of upstream's multiplication by `1/length`: 165 of 2000
  `compute` and 392 of 2000 `computeConvexPolygon` records disagree.
* `MinimumWidthPoints2`'s center `origin + e0*a0 + s1*a1` regrouped as
  `origin + (e0*a0 + s1*a1)`: 686 of 2000 `compute` and 522 of 2000
  `computeIndexed` records disagree.
* `SeparatePoints2`'s line direction normalized by division: 126 of 2000
  records disagree.
* `ConstrainedDelaunay2::SelectSplit`'s strict `<` relaxed to `<=` (the tie
  break between equal pseudo-distances): 52 `insert` and 109 `insert.strip`
  records disagree, so the combinatorial output pins the tie rule.

Provably indistinguishable, and therefore not claimed: the order in which the
two strip polygons are retriangulated (it changes only insertion order into a
set; 0 records), `ComputePSD`'s `dot1020 <= 0` against `< 0` (at zero both
branches compute the same product; 0 records), the association of the
rational center in `ConvertTo` (exact), and `0.5 * (tmin + tmax) * d` against
`0.5 * ((tmin + tmax) * d)` in the dimension-1 branches (scaling by 0.5 is
exact).

## Not covered

* **`ConvexHull3`'s multithreaded path (`lgNumThreads > 0`).** Not ported
  (the port has no threads; porting note in `src/ConvexHull3.ts`).
* **The deprecated `ConstrainedDelaunay2<InputType, ComputeType>`,
  `Delaunay2Mesh<InputType, ComputeType, RationalType>` and
  `Delaunay3Mesh<...>` specializations.** Not ported. The deprecated CDT's
  mirror-image flaw (it range-checks raw indices but never substitutes
  duplicates) is recorded under #325 and has no port counterpart.
* **`MinimumAreaBox2` overload 4 (`std::vector` points plus `std::vector`
  indices).** Deliberately not ported (port note): on an empty index vector
  it forwards to the convex-hull path, contradicting overload 3 (#328,
  minor). Overload 3 is covered.
* **`ConstrainedDelaunay2::Insert` with an index outside the array on a
  freshly constructed functor, or after a dimension 0/1 `operator()`.**
  Upstream then indexes an empty or short `mDuplicates` and dereferences
  unallocated memory, which is an access violation or a nondeterministic heap
  read rather than an exception, so no golden record can be written; the port
  throws. The deterministic form of the same defect is
  `ConstrainedDelaunay2.deviation.duplicatesRead`.
* **`MinimumWidthPoints2`'s brute-force duplicate-corner fix (#328).** Inert:
  the branch's input is always a `ConvexHull2` hull, which has no duplicates,
  so the port and upstream are the same computation there (covered by the
  main case); the defect is demonstrated through `MinimumAreaBox2`.
* **Input sets larger than 32 points.** See the suspect below. The generators
  are capped at 11 points, which also keeps records small and the exact
  fallbacks cheap.
* **#352 (the unclosed Newell loop).** It belongs to
  `MinimumVolumeBox3FloatingPoint.h`, not to `ConvexHull3.h`; `ConvexHull3`
  has no Newell normal. Nothing to demonstrate in this group.

## Upstream bug suspects

**`ConvexHull3<Real>::operator()` reports an unspecified index for a
duplicated hull vertex (new; the 3D sibling of v09's `ConvexHull2` finding).**
The points are ordered by `std::sort(sorted, lessThanPoints)` and deduplicated
by `std::unique(…, equalPoints)`, where the comparator compares the *points*.
Repeated points are equivalent under it, `std::sort` leaves equivalent
elements in an unspecified order, and which index survives - and is reported
by `GetVertices()` and `GetHull()` - is unspecified. It is invisible for at
most 32 points (MSVC's insertion-sort threshold, stable by accident), which is
the only reason the port's stable sort agrees on this family's cases. A
scratch case with 33 to 60 lattice points in `[-2,2]^3` disagreed on the
vertex indices of 13 of 20 records, and agreed on all 20 when the vertices and
faces were compared by coordinates, so the hull geometry is unaffected. Example
(57 points):

```
(1,-1,1) (1,2,1) (2,-1,-2) (-1,1,-2) (-2,2,1) (-2,1,1) (0,0,2) (2,1,-2)
(0,-1,-2) (-1,2,0) (-2,-2,-2) (0,-1,-1) (-2,1,2) (-1,-2,2) (-2,0,-1)
(1,-1,-2) (1,0,-1) (-2,-2,0) (1,0,-2) (0,-1,2) (2,-2,2) (-1,-1,0) (-2,2,-1)
(-1,-1,-2) (2,-2,-1) (2,0,0) (1,-2,2) (-1,0,-2) (-1,-2,2) (0,-2,1) (-1,0,-2)
(-1,-1,2) (1,1,0) (1,1,2) (-1,-1,-1) (-2,1,-1) (2,0,-2) (0,-2,-2) (1,-2,2)
(-1,-2,1) (0,2,2) (1,-1,-1) (1,-1,0) (0,0,2) (-1,2,2) (-2,-1,-1) (-1,2,2)
(2,0,1) (0,-2,1) (1,0,0) (-2,2,-1) (1,2,1) (2,1,-1) (0,1,-2) (1,1,-2)
(0,1,-2) (0,1,-2)
```

The MSVC build reports the hull vertex `(1,2,1)` as index 51, the port as
index 1 (`points[1] == points[51]`). `MinimumAreaBox2::GetHull()` inherits
v09's `ConvexHull2` instance of the same issue (the same scratch run with
33 to 60 lattice points in `[-3,3]^2` disagreed on the `GetHull()` indices of
11 of 20 records; the box is computed from hull coordinates, which are the
same). There is no upstream value to match (it is MSVC's introsort
partitioning), so the port is left alone.

No other new suspect. The independent checks above hold on every deep-run
record, so on the inputs where upstream is sound (every recorded finding
excluded) both implementations produce correct constrained triangulations,
hulls, boxes, widths and separations.
