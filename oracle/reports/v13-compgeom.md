# Verify group 13 (computational geometry) against the MSVC build of upstream GTE

Family `v13-compgeom`, 6 cases, 20 golden records each. Deep run
`npm run oracle:deep -- 2000 v13-compgeom` (12000 records, 4000 of them throw
records): **every case passes** (8 tests including the coverage check and
the reference-check test; replay 10-13 s, generation and replay together
about 1 min 45 s). Two deep runs of the final generator are byte-identical.

The group is `SeparatePoints3.h` and `TriangulateCDT.h`. The C++ side
instantiates exactly what the port implements:

| upstream instantiation | why |
| --- | --- |
| `SeparatePoints3<double, double>` | `ComputeType` is vestigial (the current `ConvexHull3` is templated on the input type only); the port drops it |
| `TriangulateCDT<double>` | the replacement class; the deprecated `TriangulateCDT<InputType, ComputeType>` specialization is not ported |

No compared path calls the C math library: outputs are combinatorial
(booleans, indices, counts, node structure, triangle lists) or computed with
`+ - * / sqrt` on doubles (the separating plane is `UnitCross` and `Dot`).
**All 6 cases are declared `{ exact: true }`; every real output of every
non-deviation case is bit-identical to the MSVC value** (2888 real outputs in
the deep run of `SeparatePoints3.compute`, 2888 exact; no tolerance anywhere
in the family).

**No port defect was found.** This branch contains no `src/` change.

## Canonicalization and order dependence

`TriangulateCDT` fills `PolygonTreeEx::allTriangles` from
`ConstrainedDelaunay2::GetIndices()` and `outsideTriangles` from
`ETManifoldMesh::GetTriangles()`, both `std::unordered_map` orders; the port
enumerates them in sorted-key order. Both lists are sorted by their restored
vertex tuple on both sides. Everything else is emitted verbatim:
`Node::triangulation` comes from a `std::set<TriangleKey<true>>` (defined
order, which the port reproduces with a sorted key set), and the interior,
exterior and inside lists are built from the node triangulations in node
order. The seed triangle chosen in `ClassifyDFS` (`edge->T[0]` or `T[1]`)
does not depend on the graph copy's insertion order: exactly one of the two
lies on the node's side.

`SeparatePoints3` visits the faces of `ConvexHull3::GetHull()`. v11 showed
the hull's triangle *set* equals the port's for at most 32 points (the sets
here have at most 8). Its *order* is `std::unordered_map` order, and for
more than four points it is not even reproducible within one process:
`ConvexHull3::Hull3` starts each visible-region search from
`Vertex::TAdjacent`, a `std::unordered_set<Triangle*>`, so the insertion
history of the hull mesh depends on heap addresses. Measured with two calls
on the same input in one process, 2000 uniform sets each: 8-point sets give
a different `GetHull()` order on 850 (identical triangle set on all 2000),
4-point sets on none; `SeparatePoints3` returns a different plane on 95 of
2000 pairs. Upstream returns the *first* separating face in that order, and
`WhichSide` returns the *first* nonzero sign in that order, so the returned
plane and, on the inputs of #348 finding 1, the boolean itself depend on the
run. See "Upstream bug suspects". The cases are built so that no golden
depends on that order:

* The main case accepts an input only when upstream's float answer provably
  equals the port's exact answer in *every* face order (the predicate
  `TraceSeparate3`, below), and emits the plane only when it is the same in
  every order.
* Both deviation cases use a set of exactly four points where the order
  matters (a tetrahedron's hull is built by four fixed insertions and never
  searched). A first version of `deviation.edgeAxisOrigin` with up to 8
  points in the second set produced 7 different records in two identical
  deep runs; with the tetrahedron the deep runs are byte-identical.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `SeparatePoints3.h` | `SeparatePoints3.compute` (`operator()`: the degenerate early return, both face loops, the edge-edge loop, no separation; separated, touching, overlapping, coplanar and degenerate clouds) | exact | pass, 2000 records |
| | `SeparatePoints3.deviation.roundoff` (#348 finding 1) | deviation | pass, 2000 of 2000 deviate |
| | `SeparatePoints3.deviation.edgeAxisOrigin` (#348 finding 2) | deviation | pass, 1855 of 2000 deviate |
| `TriangulateCDT.h` | `TriangulateCDT.compute` (both `operator()` overloads, alternating; every `PolygonTreeEx` member: per node `self`, `chirality`, `parent`, `minChild`, `supChild`, the output polygon with edge splits, the triangulation; `interiorTriangles` / `interiorNodeIndices`, `exteriorTriangles` / `exteriorNodeIndices`, `insideTriangles` / `insideNodeIndices`, `outsideTriangles`, `allTriangles`) | exact | pass, 2000 records |
| | `TriangulateCDT.compute.invalidThrows` ("Invalid argument." and "Invalid polygon tree.") | exact | pass, 2000 throw records |
| | `TriangulateCDT.compute.hullSharedEdgeThrows` (new suspect below) | exact | pass, 2000 throw records |

### SeparatePoints3: generator, soundness predicate and branches

Each record draws two sets of 4 to 8 points (1 to 8 in the degenerate mode)
in one of six modes by record index: lattice sets in `[-3,3]^3` shifted by a
lattice vector in `[-6,6]^3`; uniform sets in `[-4,4]^3` with an integer
shift; crossed tetrahedra with skew edges at height `d` in `{-1,0,1,2}`
(overlap, touching at the crossing point, edge-edge separation), lattice or
jittered by 0.2, under a random signed axis permutation, with edge midpoints
as extra points; points on two spheres with a gap in `[-0.5,0.5]`; lattice
boxes sharing a face, an edge or a vertex, a unit apart or overlapping by a
unit, on a third of the records under a common integer shear (non-axis
faces); and degenerate sets (all equal, collinear, a lattice plane with a
non-axis normal, at most three points, both sets in one plane). The set
sizes are recorded after the draw (the crossed-tetrahedra mode may swap the
sets; recording the drawn counts first was a harness error of the first
version, found on the first replay).

`TraceSeparate3` classifies every candidate face both ways: the port's exact
test (`OnSameSide` of the other hull is +1; `WhichSide` of the owner is -1
exactly for an outward normal), and upstream's float test in every order.
`OnSameSide`'s value is order independent; the float `WhichSide` can be any
nonzero sign that occurs among the owner's vertices, so a face's float test
can pass in some order when the other hull's float side is nonzero and the
owner's float signs contain its negative (or are all zero), and passes in
every order when they do not contain it. An input is accepted when

* some face separates exactly, some face passes the float test in every
  order, and every face that can pass in a loop upstream can reach separates
  exactly (the second condition makes the boolean agree; without the third,
  upstream can return `true` with a non-separating plane, and the first deep
  run without it showed exactly that on 240 of 2000 records through the
  reference check); or
* no face separates exactly, no face can pass the float test in any order,
  and along the edge-edge loop (a `std::set`, defined order) the float and
  exact `side0 * side1 < 0` agree until the loop returns.

The plane is emitted (normal and constant, plus the origin on a face stage)
only when it is *determined*: the returning face loop has exactly one face
that can pass, it separates exactly and passes in every order (and for the
second loop no face of the first can pass), or the stage is edge-edge. The
stage and the flag are recorded as inputs after the points and only steer
which outputs are emitted. The C++ side checks upstream against the trace on
every record (same boolean; same bits for a determined plane) and would
record a throw on a mismatch; none occurred.

Acceptance per draw (capped at 64 draws, everything redrawn, fallback: two
far-apart unit cubes): lattice 35%, uniform 13%, crossed tetrahedra 52%,
spheres 11%, boxes 94%, degenerate 95%; 2 fallbacks in 2000 records. The
rejected draws are the inputs on which upstream's answer depends on its hull
order (a face whose own vertex rounds to the wrong side).

Deep-run histogram of `SeparatePoints3.compute` (2000 records):

| stage | records |
| --- | --- |
| a hull of dimension below 3 (hull dimension pairs 0/0 to 2/3 all reached) | 352 |
| separated by a face of hull 0 | 1241 |
| separated by a face of hull 1 | 83 |
| separated by an edge-edge plane | 113 |
| not separated, both hulls 3D (overlapping) | 211 |

1437 records are separated, 229 of them touching (the closed hulls meet,
checked exactly); 563 are not. The plane is compared bit for bit on 461
separated records (determined); on the other 976 separated records upstream
can return more than one face depending on its hull order (all of them
genuine separators, by the predicate), and only the boolean and the
reference check are compared.

### TriangulateCDT: generator and branches

The pool and the polygon tree are recorded (pool size and points, then the
tree in preorder: polygon size, indices, child count). Six modes by record
index: one simple lattice polygon (rectangles with every boundary lattice
point as a vertex, and 12-direction stars whose angular gaps are closed to at
most 150 degrees); a 12-direction outer polygon with one or two star holes;
concentric nesting two to four levels deep plus an extra hole and an extra
island; coincident configurations (a hole sharing a vertex with the outer
polygon, a hole vertex in the interior of an outer edge, a hole sharing the
bottom edge of a U-shaped polygon's notch entirely or in part, two holes
sharing an edge or a vertex, an island touching its hole at a vertex, a
pinched polygon visiting one vertex twice, collinear vertices on both outer
polygon and hole); the same with about half of the references replaced by
duplicated pool points (so a shared vertex is referenced through two indices
with equal coordinates, including -0 for +0); and non-lattice star polygons
three levels deep. Every record adds 0 to 3 unused pool points and shuffles
the pool, so the order in which `RemapPolygonTree` meets the points differs
from the pool order. The unique points passed to `ConstrainedDelaunay2` are
drawn under v09/v11's `Sound2` predicate (the #391 restriction inherited from
`Delaunay2`), computed on exactly the point sequence `RemapPolygonTree`
builds.

Deep-run histogram of `TriangulateCDT.compute` (2000 records):

| branch | records |
| --- | --- |
| tree with 1 / 2 / 3 / 4 / 5 / 6 nodes | 536 / 715 / 492 / 137 / 87 / 33 |
| exterior triangles (holes) | 1464 |
| outside triangles (non-convex root) | 726 |
| an output polygon longer than its input (`Insert` split an edge at a vertex) | 140 |
| duplicate coordinates referenced (`RemapPolygonTree`'s duplicate branch) | 128 |
| a restored polygon index differs from the input (#348 finding 3, preserved) | 109 |
| a polygon repeating an index (pinched polygon) | 39 |
| -0 among the referenced coordinates | 147 |
| unused pool points | 1587 |
| `allTriangles` count below 10 / 10-19 / 20 or more | 531 / 507 / 962 |

One record of an early deep run was a self-intersecting star (an angular gap
wider than 180 degrees), which upstream documents as unsupported: both sides
threw "Unexpected condition.". The generator now closes such gaps.

## Port defects fixed

None. Every association order, comparison and container order in
`src/SeparatePoints3.ts` and `src/TriangulateCDT.ts` matches upstream bit for
bit on all deep-run records, so this branch carries no `src/` change.

Sensitivity, measured by mutating `src/` and replaying the deep run:

| mutation | records that disagree |
| --- | --- |
| edge-edge plane constant through the other endpoint of edge 0 | 23 of 2000 `SeparatePoints3.compute` |
| face plane built from the rotated triple `(P1, P2, P0)` | 225 of 2000 |
| edge set iterated in descending order | 53 of 2000 |
| `RemapPolygonTree` quirk "fixed" (first occurrence kept) | 128 of 2000 `TriangulateCDT.compute` |
| exterior triangles stored as `(V2, V1, V0)` instead of `(V0, V2, V1)` | 1464 of 2000 |
| region triangles extracted in reverse key order | 1987 of 2000 |

## Deliberate deviations demonstrated

1. **`SeparatePoints3.deviation.roundoff` (#348 finding 1), 2000 of 2000.**
   Set 0 is a tetrahedron (reproducible hull order), set 1 a cloud of 4 to 8
   points, both in `[-1,1]^3` with a shift of at most 0.5. A draw is accepted
   only when (a) the port's exact algorithm reports no separation, (b) the
   closed hulls meet, exhibited exactly by the independent reference,
   (c) a verbatim replica of upstream's first face loop, on the same
   deterministic order, passes the float test and upstream reports a
   separation, and (d) upstream's plane leaves points of one set strictly on
   both sides (residual above 1e-12). Every record is a genuinely wrong
   upstream answer; the port reports `false` on all of them. The fallback is
   a recorded accepted draw.
2. **`SeparatePoints3.deviation.edgeAxisOrigin` (#348 finding 2), 1855 of
   2000.** Crossed tetrahedra separated (or touching) only by the edge-edge
   plane, with set 1 a tetrahedron because the stale origin is that of the
   last face of hull 1 in hull order. `separated`, the normal and the
   constant agree bit for bit; the origin differs. Residue explained record by
   record: the 145 agreeing records all have plane constant 0 (the plane
   passes through the world origin, so the port's origin is `0 * normal`) and
   upstream's stale face plane also passes through it with the same zero
   signs.

The fixes are confined: the main case agrees bit for bit on every one of its
records, including 461 determined planes and 113 edge-edge planes.

## Independent reference checks (agreement is not correctness)

Each main case ends with a reference check that each side computes on its own
output with exact arithmetic (BSNumber on the C++ side, bigints on the common
dyadic scale of `test/helpers/exact.ts` on the TypeScript side), and the
replay's last test requires that the port's check never failed. On the final
deep run all checks hold on both sides:

* **SeparatePoints3** (2000 records): on all 1437 `true` records the plane is
  unit length and puts one set on each closed side, largest residual
  3.4e-16 of the largest coordinate; on every one of the 211 nondegenerate
  `false` records the closed hulls meet (a point of one set in the other
  hull, or a hull edge meeting a hull triangle, exactly, coplanar contacts
  included); the other 352 `false` records have a hull of dimension below 3.
* **TriangulateCDT** (2000 records): every node triangle has the node's
  chirality strictly (exact orientation); per node, the exact signed areas of
  its triangles sum to the node polygon's shoelace area plus its children's
  (rational shoelace, so the triangulation covers exactly the region between
  the polygon and its children); every triangle's centroid has winding
  number equal to the chirality in its node polygon and 0 in every child
  polygon (no triangle in a hole, none in an island), and every outside
  triangle's centroid has winding number 0 in the root polygon; every edge of
  every output polygon (after splits) is an edge of `allTriangles`;
  interior plus exterior equals inside in order, and inside plus outside
  equals `allTriangles` as a multiset of oriented triangles.

The throw-parity cases carry no reference; the port's messages were checked
on sample records ("Invalid argument.", "Invalid polygon tree.",
"Unexpected condition.").

## Not covered

* **The deprecated `TriangulateCDT<InputType, ComputeType>`.** Not ported
  (port note in `src/TriangulateCDT.ts`); its body is the replacement's with
  `int32_t` counts and an explicit `ConstrainedDelaunay2<InputType,
  ComputeType>`.
* **`SeparatePoints3` with a `ComputeType` other than `double`.** Vestigial
  upstream; the instantiation computes nothing different.
* **`PolygonTreeEx::GetContainingTriangle` and the rest of `PolygonTree.h`.**
  Group 16, not this group.
* **`SeparatePoints3` inputs on which upstream's result depends on its hull
  order.** Rejected by the main generator (see the acceptance rates) because
  no reproducible golden exists for them; the defect they carry is
  demonstrated in `deviation.roundoff` on the reproducible (tetrahedron) form.
  Likewise, on the 976 separated records with several separating faces the
  plane is compared only through the reference check.
* **`TriangulateCDT` on a tree whose referenced points are collinear.**
  `Delaunay2::operator()` returns at dimension 1 before sizing
  `mDuplicates`, and `ConstrainedDelaunay2::Insert` then reads
  `duplicates[edge[0]]` from an empty vector, an access violation or heap
  read rather than an exception, so no golden can be written (v11's #325
  finding 2 in a new caller). The port throws.
* **The pointer overload with `numInputPoints` below the largest referenced
  index.** Upstream reads the caller's array past `numInputPoints` and
  writes `remapping` out of range; not generable. The port has only the
  array form.
* **Self-intersecting polygon trees** (documented unsupported by upstream).
  Seen once, above; both sides threw.
* **More than 32 points per set.** `ConvexHull3`'s duplicate-index finding
  (#325) and the heap-address order both apply; the sets here have at most 8
  points. `TriangulateCDT`'s pools reach 79 points, but nothing on its path
  (`Delaunay2`, `ConstrainedDelaunay2`, the manifold meshes) sorts.

## Upstream bug suspects

1. **`ConvexHull3<Real>::operator()` returns its dimension-3 hull in a
   heap-address-dependent order, and `SeparatePoints3` inherits it (new;
   minor for `ConvexHull3`, see below for `SeparatePoints3`).**
   `ConvexHull3::Hull3` seeds each visible-region search with the first
   visible triangle found by iterating `Vertex::TAdjacent`, which is a
   `std::unordered_set<Triangle*>` keyed by pointer, and removes and inserts
   hull triangles in the resulting breadth-first order; `GetHull()` and
   `GetVertices()` then enumerate the mesh's `std::unordered_map`s, whose
   layout depends on that history. Two calls on the same 8 uniform points in
   one process return the same triangle set in a different order on 850 of
   2000 sets (4 points: never, the initial tetrahedron is inserted in a fixed
   sequence). `SeparatePoints3` returns the first separating face in that
   order, so two calls return different (valid) planes on 95 of 2000 pairs,
   and the #348 false separation fires or not depending on whether
   `WhichSide` meets the rounded face vertex before a clearly negative one:
   the defect is not reproducible run to run. The port enumerates the hull in
   sorted-key order (deterministic); nothing to change. Proposed text for
   `docs/UPSTREAM-FINDINGS.md` (under `ConvexHull3.h` / `SeparatePoints3.h`,
   issue #325 or #348):

   > **`ConvexHull3::operator()` returns the dimension-3 hull in a
   > heap-address-dependent order (minor).** `Hull3` seeds each visible-region
   > search from `Vertex::TAdjacent`, a `std::unordered_set<Triangle*>`, so
   > the hull mesh's insertion history, and with it the `std::unordered_map`
   > order that `GetHull()` and `GetVertices()` enumerate, varies between
   > calls on the same input in one process (850 of 2000 8-point sets; never
   > for 4 points). The triangle set is the same. `SeparatePoints3` returns
   > the first separating face in this order, so its plane differs between
   > two calls on 95 of 2000 separated pairs, and the #348 false separation
   > occurs only in some orders. Found by the C++ oracle of group 13. Port:
   > sorted-key order, deterministic.

2. **`TriangulateCDT<T>::operator()` throws on a documented-supported
   coincident edge-edge configuration (new; WR).** When a hole shares an edge,
   or part of an edge, with its parent polygon and that edge lies on the
   convex hull of the referenced points, `ClassifyDFS` processes the hole
   first and removes its triangles; `ETManifoldMesh::Remove` deletes the
   shared edge because no triangle remains on its other side, and the parent's
   pass fails `LogAssert(eiter != emap.end(), "Unexpected condition.")`.
   Example: outer square `(0,0),(6,0),(6,6),(0,6)`, hole `(0,2),(0,4),(2,3)`
   (or `(0,0),(0,6),(2,3)`, the full side). The header states "The algorithm
   supports coincident vertex-edge and coincident edge-edge
   configurations"; the same configuration with the shared edge inside the
   hull (a hole under the notch of a U-shaped polygon) succeeds and is
   covered by `TriangulateCDT.compute`. Loud, so WR. A fix would skip
   polygon edges that are no longer in the graph (they have no triangle on
   the node's side): `if (eiter == emap.end()) { continue; }`. The port
   preserves the assert (throw parity, `TriangulateCDT.compute.hullSharedEdgeThrows`,
   2000 of 2000 records throw on both sides); whether to fix is left to the
   orchestrator. Proposed text for `docs/UPSTREAM-FINDINGS.md` (under
   `TriangulateCDT.h`, issue #348):

   > **4. `TriangulateCDT::ClassifyDFS` rejects a shared edge on the convex
   > hull (WR).** A hole sharing an edge, or part of one, with its parent
   > polygon where that edge is on the convex hull of the referenced points:
   > the hole's triangles are extracted first, `ETManifoldMesh::Remove`
   > deletes the edge (no triangle left on the other side), and the parent's
   > `LogAssert(eiter != emap.end())` throws "Unexpected condition." although
   > the header documents coincident edge-edge support. Example: square
   > `(0,0),(6,0),(6,6),(0,6)` with hole `(0,2),(0,4),(2,3)`. Interior shared
   > edges work. Suggested fix: skip polygon edges absent from the graph.
   > Found by the C++ oracle of group 13. Port: preserved (throw parity).

No other new suspect. The two preserved minors of #348 finding 3 (the
`RemapPolygonTree` overwrite, visible on 109 deep-run records, and the
redundant re-insertion in `ConstrainedTriangulate`) agree bit for bit.
