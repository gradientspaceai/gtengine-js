# Verify group 10 (computational geometry) against the MSVC build of upstream GTE

Family `v10-compgeom`, 18 cases, 20 golden records each. Deep run
`npm run oracle:deep -- 2000 v10-compgeom` (36000 records, 7719 of them throw
records): **every case passes**. Wall time of generation and replay
together: 7 s.

The group is `ExtremalQuery3BSP.h`, `IncrementalDelaunay2.h`,
`MinimumAreaCircle2.h`, `MinimumVolumeSphere3.h` and `RotatingCalipers.h`.
Nothing on any compared path calls the C math library except `sqrt` (the
circle and sphere radii, `UnitCross` of the BSP face normals), which IEEE 754
rounds exactly. **All 14 non-deviation cases are declared `{ exact: true }`;
every one of their 342318 deep-run real outputs is bit-identical to the MSVC
value. No tolerance appears anywhere in the family.**

Two port defects were found and fixed (a 1 ulp association error in
`MinimumVolumeSphere3` and holes in the array `IncrementalDelaunay2.getHull`
returned for a collinear triangulation), and two new upstream findings are
reported below (`ExtremalQuery3BSP` is nondeterministic from run to run;
`IncrementalDelaunay2::GetHull` silently returns garbage for collinear input
where issue #290 only described the infinite loop).

## What is comparable, and how

**The `std::shuffle` of `MinimumAreaCircle2` / `MinimumVolumeSphere3`.**
`operator()` permutes the unique points with `std::shuffle(permuted, mDRE)`
and `mDRE` is a default-constructed `std::default_random_engine`, which is
`std::mt19937` on MSVC. The port uses a Lehmer generator (`16807 * x mod
2^31-1`, seed 1) with its own Fisher-Yates loop. ORACLE.md lists such paths
as not comparable. They are comparable after all: the C++ side makes
upstream run *the port's* permutation.

* `PortPermutation(m)` reproduces the permutation the port's first
  `compute()` applies to `m` unique points (the replay constructs a fresh
  query per record).
* A private `std::default_random_engine` runs in lockstep with the query's
  `mDRE`: both start in the default state, `mDRE` serves exactly one
  `std::shuffle` of `m` elements per `operator()` call, and the engine
  consumption of `std::shuffle` depends only on the length, so shuffling a
  local index array once per call yields exactly the permutation that call
  applied.
* The query is called until the tracked permutation equals the port's, and
  *that* call's result is recorded. Before the shuffle both sides hold the
  same array (sorted unique indices; MSVC's `std::sort` is an insertion
  sort, hence stable, below 32 elements, and the port sorts with an explicit
  index tie-break), so from the shuffle on both run the same computation on
  the same point order, and every output - center, radius, support count and
  the support indices *in upstream's order* - is compared bit for bit.

The expected number of calls is `m!`, so the unique count is capped at 6 (720
calls on average; the loop is capped at 20000, and the probability of not
meeting a given permutation of 6 is below `e^-27`). The earlier draft of this
group accepted only point sets whose result did not change over 48 upstream
shuffles; that filter excluded exactly the interesting inputs. Measured on
the deep run by re-running each recorded input through 24 further port
permutations, the result of **529 of 2000** `MinimumAreaCircle2.compute`
records, **1835 of 2000** `.circumcircle` records and **865 of 2000**
`MinimumVolumeSphere3.compute` records depends on the permutation (in 147,
3 and 134 of them the support *set* changes, in the rest the last bits of
the center or radius, because `ExactCircle3` / `ExactSphere4` are not
symmetric in their arguments). The port defect below is only visible on such
records.

**`ExtremalQuery3BSP`.** The BSP tree is built from `VETManifoldMesh`, whose
`GetEdges()` / `GetVertices()` are `std::unordered_map` and whose
`Vertex::TAdjacent` is an `std::unordered_set<Triangle*>`; the port reads
them through sorted accessors. The trees differ (and upstream's differs from
one construction to the next, see the first upstream finding), so only the
query answer is compared, and only on directions whose extreme vertex is
unique for `D` and `-D` (a tie is decided by the leaf region, which is tree
shape) and on which upstream's own answer is the exact argmax (a guard for
issue #290; it never fired on the polytopes drawn, see below).
`GetFaceNormals()` is compared bit for bit.

**`IncrementalDelaunay2`.** `mGraph` is a `VETManifoldMesh`, so the triangle
numbering of `GetTriangles()` / `GetAdjacencies()` / `GetTriangle(t)` is
hash-table order. The cases canonicalize exactly as v09 does for `Delaunay2`:
triangles sorted by their stored vertex tuple, adjacency indices remapped
through the same permutation. The triangle *set* is order independent: the
insertion polygon is the connected set of triangles whose circumcircle
contains the point, which does not depend on the start triangle even for a
point on an edge, and `DeleteRemovalPolygon` orders the removal polygon
through a `std::map`. `GetContainingTriangle` starts at a recorded canonical
rank that each side maps to its own index, and the whole `SearchInfo`
(`initialTriangle`, `finalTriangle`, `finalV`, `numPath`, `path`) is emitted
translated back to canonical ranks, so query points on a shared edge or at a
vertex are compared too (3980 of the 8000 deep-run queries).

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `RotatingCalipers.h` | `RotatingCalipers.computeAntipodes` (the only public entry point; strictly convex polygons from cocircular lattice points, lattice hulls and real hulls) | exact | pass, 2000 |
| | `….computeAntipodes.collinear` (exact edge midpoints inserted, the collinear-removal path of `CreatePolygon`) | exact | pass, 2000 |
| | `….computeAntipodes.collinearThrows` (2 to 6 collinear points: throw parity for the `indices.size() >= 3` assert) | exact | pass, 2000 throw records |
| | `….computeAntipodes.deviation.duplicate` (#286) | deviation | 2000 of 2000 deviate |
| `MinimumAreaCircle2.h` | `MinimumAreaCircle2.compute` (`operator()`, `GetNumSupport`, `GetSupport`; uniform, lattice, cocircular lattice and dense lattice modes, 1 to 8 points, at most 6 unique) | exact | pass, 2000 |
| | `….compute.circumcircle` (acute lattice triangle plus points strictly inside its circumcircle: three-point support on every record) | exact | pass, 2000 |
| | `….compute.empty` (throw parity for no points) | exact | pass, 2000 throw records |
| | `….compute.deviation.trappedFallback` (#286) | deviation | 1998 of 2000 deviate |
| `MinimumVolumeSphere3.h` | `MinimumVolumeSphere3.compute` (as above plus a cospherical lattice mode, the 30 lattice points of `x^2+y^2+z^2 = 9`) | exact | pass, 2000 |
| | `….compute.empty` | exact | pass, 2000 throw records |
| | `….compute.deviation.trappedFallback` (#286) | deviation | 2000 of 2000 deviate |
| `ExtremalQuery3BSP.h` | `ExtremalQuery3BSP.getExtremeVertices` (constructor, `GetExtremeVertices`, `GetFaceNormals` of the base; tetrahedron and bipyramid under integer linear maps, the regular octahedron; lattice and real directions) | exact | pass, 2000 |
| `IncrementalDelaunay2.h` | `IncrementalDelaunay2.insert` (`Insert` return values, `GetNumVertices`, `GetNumTriangles`, `GetTriangles`, `GetAdjacencies`, `GetHull`) | exact | pass, 2000 |
| | `….remove` (`Remove` of vertices, non-vertices and removed positions, then re-insertion, which gets a fresh index) | exact | pass, 2000 |
| | `….getContainingTriangle` (`GetTriangulation`, `GetContainingTriangle` with the whole `SearchInfo`, `GetTriangle`, `GetAdjacent`, both in and out of range) | exact | pass, 2000 |
| | `….finalizeTriangulation` (`FinalizeTriangulation`, the state after it, a second call, `Insert` / `Remove` after it) | exact | pass, 2000 |
| | `….domainAsserts` (constructor, `Insert`, `Remove` domain preconditions; `Remove` after finalization) | exact | pass, 2000 (1500 throw records) |
| | `….getHull.deviation.collinear` (#290) | deviation | 1151 of 2000 deviate |

Half of the `IncrementalDelaunay2` records use real rather than lattice
domain bounds, so the supervertex arithmetic of the constructor
(`xMin - dx`, `xMin + 5 * dx`, ...) is inexact and `GetTriangulation`'s
vertex list pins its evaluation order.

### Degeneracy actually reached (deep run, 2000 records per case)

| measurement | value |
| --- | --- |
| `MinimumAreaCircle2.compute` support size 1 / 2 / 3 | 319 / 1176 / 505 |
| `MinimumVolumeSphere3.compute` support size 1 / 2 / 3 / 4 | 313 / 854 / 667 / 166 |
| `.circumcircle` support size 3 | 2000 |
| `RotatingCalipers` records with a repeated antipode edge (the preserved #286 tie quirk) | 51 (main), 54 (`.collinear`) |
| `IncrementalDelaunay2.insert` triangles per record | 4 to 18, 20936 in all |
| `IncrementalDelaunay2.remove` successful removals / fresh re-insertion indices | 3022 / 1781 |
| `getContainingTriangle` queries found / on an edge or vertex / outside the hull | 7045 / 3980 / 955 |
| `finalizeTriangulation` Delaunay triangles per record | 1 to 10, 7750 in all |

## Port defects fixed

1. **`MinimumVolumeSphere3.exactSphere4` summed the center pairwise.**
   Upstream writes `X[0]*P0 + X[1]*P1 + X[2]*P2 + x3*P3`, which accumulates
   left to right, `((a + b) + c) + d`; the port computed
   `(a + b) + (c + d)`. The center of a four-point support was 1 ulp off on
   the first record that reached it (center x `0.018733754116927015` where
   the MSVC build returns `0.01873375411692702`, scaled error 3.5e-18); with
   the pairwise grouping restored, 50 of the 2000 deep-run records of
   `MinimumVolumeSphere3.compute` disagree. It went unnoticed before because
   the earlier permutation-invariance filter rejected the records whose last
   bits depend on the argument order, which is every record where it shows.
   Fixed in `src/MinimumVolumeSphere3.ts`; the regression test in
   `test/MinimumVolumeSphere3.test.ts` pins the MSVC center, radius and
   support order for that input.

2. **`IncrementalDelaunay2.getHull` returned an array with holes.** For a
   finalized triangulation of three or more collinear points the edge map
   collected from the supervertex triangles is not one cycle, and the walk
   can return to its start before visiting every edge. Upstream then returns
   `hull[]` with `numEdges` entries whose unwritten tail keeps the zeros of
   the `resize` (supervertex 0; for the points `(4,4) (2,2) (0,0)` in the
   domain `(-6,-8)-(6,6)` the MSVC build returns `[7, 8, 0]`). The port's
   deliberate #290 fix caught only the other outcome (the walk never
   returning) and here returned `[7, 8, <hole>]`, an array with an
   `undefined` entry, which is neither upstream's value nor a hull. The port
   now reports the degenerate triangulation after the walk as well
   (`i + 1 === numEdges`, which always holds for a nondegenerate
   triangulation, so the check changes nothing where upstream is sound: the
   `insert`, `remove` and `finalizeTriangulation` cases call `GetHull` on
   every record and agree bit for bit). 697 of the 2000 deep-run records of
   the deviation case take this path. Regression test in
   `test/IncrementalDelaunay2.test.ts`.

## Deliberate deviations demonstrated

1. **`RotatingCalipers.computeAntipodes.deviation.duplicate` (#286), 2000 of
   2000 deviate.** `CreatePolygon` tests collinearity against the
   immediately preceding edge of the input array, so a duplicated vertex
   gives a zero edge and the next genuine corner is discarded; the port
   compares against the most recent nonzero edge. On a convex polygon with
   one vertex repeated upstream answers with one corner fewer, and in 219
   records that leaves fewer than three corners and upstream's `LogAssert`
   fires, while the port returns the antipodes of the whole polygon.

2. **`MinimumAreaCircle2.compute.deviation.trappedFallback` (#286), 1998 of
   2000.** In the trapped-failure branch upstream calls
   `GetContainer(numPoints, points)` after `numPoints` was overwritten with
   the unique count, so the fallback circle bounds only a prefix of the input
   array; the port bounds all of it. The inputs are catalogued lattice sets
   (4 to 6 points) on which the port traps, presented in a random order with
   1 to 3 copies of one point prefixed; the set of distinct points and hence
   the control flow is unchanged, and the lockstep permutation makes upstream
   trap at the same point (the C++ `success` is `false` on all 2000 records).
   The 2 agreeing records are coincidences of the data: the prefix
   `(0,4) (0,4) (3,-2) (2,0) (0,4)` has the same centroid `(1,2)` as the
   whole array and contains its farthest point `(3,-2)`.

3. **`MinimumVolumeSphere3.compute.deviation.trappedFallback` (#286), 2000
   of 2000**, the same in 3D (catalogued lattice sets of 5 and 6 points in
   `[-2,2]^3`; the C++ side traps on every record).

4. **`IncrementalDelaunay2.getHull.deviation.collinear` (#290), 1151 of
   2000.** Collinear lattice points, finalized. The C++ side runs a verbatim
   copy of `GetHull`'s edge collection and walk as a probe and calls the
   *real* upstream `GetHull` whenever that is safe (the walk returns to its
   start within `numEdges` steps): 697 records, upstream's padded output as
   described under port defect 2. Otherwise upstream would write past the end
   of `hull[]` and never terminate, and the record emits what upstream writes
   up to the end of `hull[]`: 454 records. The port throws on both. The 849
   agreeing records are the residue with one distinct point (111: no hull
   edge at all; upstream would dereference `edges.begin()` of an empty map,
   which is not executed, and the replica and the port both report an empty
   hull) or two distinct points (738: a 2-cycle through both edges, on which
   upstream is sound).

The confinement of the fixes is shown by the main cases, which run broad
generators and agree bit for bit: `MinimumAreaCircle2.compute` and
`MinimumVolumeSphere3.compute` accept trapped failures without duplicates
(upstream's prefix is then the whole array), and the `CreatePolygon` fix is
inactive in the `.collinear` case, whose inserted midpoints are exactly
collinear but not duplicates.

## Independent reference checks (agreement is not correctness)

The deep-run records were re-checked outside both implementations, with
exact BigInt predicates over the dyadic inputs where the question is
combinatorial:

* **Minimum circle / sphere** (6000 records): every successful result
  contains every input point (to 1e-9; the radius is a `sqrt`), every
  support point lies on the boundary, and no ball through 2 to `d+1` of the
  points is smaller than the result (brute force). One record fails:
  `.circumcircle` record 945, points `(4,4) (-6,3) (-3,-3) (-2,3) (-4,4)
  (-3,5)`, where both implementations return `success = true` and the circle
  `(0,1)`, radius 5 through four cocircular lattice points, which misses
  `(-6,3)` by 1.32 and whose support list `(5,2,1)` names that point. This is
  the preserved defect #399 (cocircular input; a rejected update rewrites the
  support set), reproduced bit for bit on both sides; not new.
* **`IncrementalDelaunay2`** (`insert`, `finalizeTriangulation`, `remove`,
  6000 records): every triangle is exactly counterclockwise, no present
  vertex is strictly inside any circumcircle (exact in-circle determinant),
  adjacency is symmetric, the hull turns left, the `remove` return values
  follow an independently replayed vertex numbering, and before
  finalization the exact triangle areas sum to the exact area of the domain
  rectangle.
* **`getContainingTriangle`** (8000 queries): a returned triangle contains
  the query point by the exact orientation tests; an `invalid` answer is
  outside the rectangle.
* **`RotatingCalipers`** (20357 antipodes of the two main cases): for every
  antipode the whole polygon lies on the inner side of the antipodal edge and
  the vertex is at maximum distance from its line (exact).
* **`ExtremalQuery3BSP`**: the main case only accepts directions whose
  answer is the exact brute-force argmax, and the port agrees on all of them.
  Separately, over 300 polytopes of the drawn families plus sheared
  octahedra and 400 random directions each (120000 queries), the port's
  answer was the exact argmax every time; see the #290 measurement below.

## Sensitivity of the cases

The combinatorial outputs (antipodes, triangulations, search paths, support
sets) are decided by exact rational predicates on both sides and cannot move
with a change of floating-point association. The floating-point outputs were
tested by mutating `src/` and replaying the deep run:

| mutation | deep records that disagree |
| --- | --- |
| `exactSphere4` center pairwise (the defect) | 50 (`MinimumVolumeSphere3.compute`) |
| `exactSphere4` `tmp` as `a + (b + c)` | 19 |
| `exactSphere3` (MVS) center as `a + (b + c)` | 268 |
| `exactCircle3` center as `a + (b + c)` | 205 (`compute`), 656 (`.circumcircle`) |
| `exactCircle3` `x2 = 1 - (X0 + X1)` | 195, 604 |
| `exactCircle2` center `0.5*P0 + 0.5*P1` | 0 |

The last row cannot be discriminated by any input short of underflow:
`0.5 * (a + b)` and `0.5*a + 0.5*b` are the same double because halving is
exact.

## Not covered

* **`ExtremalQuery3BSP::GetNumNodes`, `GetTreeDepth`.** The trees differ
  (container order), and upstream's tree is not even reproducible from one
  construction to the next (first upstream finding below), so there is no
  upstream value to match.
* **`ExtremalQuery3BSP` on sheared octahedra, and the directions on which
  #290 fires.** A draft `deviation` case selected directions on which the
  MSVC build answers a non-extreme vertex (the port, whose tree differs, was
  right on all of them), but which directions those are depends on heap
  addresses, so the case regenerated different goldens on every run and was
  dropped. The main case draws the tetrahedron and the triangular bipyramid
  under integer linear maps and the regular octahedron, on which the #290
  guard never fired (0 wrong answers in 80400 + 10800 measured directions).
  Directions whose extreme vertex is not unique are excluded as well: the
  tie is resolved by the leaf region, which is tree shape.
* **The #290 construction assert (`VETManifoldMesh::Insert` failure, a null
  `Edge::T[1]`; the port asserts).** Upstream cannot be executed on such
  input: `SortAdjacentTriangles` follows `tri->T[prev]` into a null pointer
  for an open mesh (an access violation, not an exception), so no
  `deviation` case can record it.
* **`MinimumAreaCircle2` / `MinimumVolumeSphere3` with more than 6 unique
  points.** The lockstep reproduction of the port's permutation costs `m!`
  upstream calls on average; 7 unique points would cost 5040 per record.
  Nothing in the algorithm changes above 6 points (the support is at most 3
  or 4 points), and the four-point sphere support is reached on 166 deep-run
  records. Above 32 points MSVC's `std::sort` is no longer stable and the
  representative index of a repeated point would be unspecified (the v09
  `ConvexHull2` finding).
* **`MinimumAreaCircle2` / `MinimumVolumeSphere3` with a null pointer and a
  positive count.** The port's API takes an array; the empty call is covered
  (`.empty`).
* **`RotatingCalipers::ComputeAntipodes` with 0 or 1 vertex.** Upstream
  reads `vertices.back()` / `vertices[1]` out of bounds before its size
  assert (#286, port asserts up front); 2 to 6 collinear points are the
  throw-parity case.
* **`IncrementalDelaunay2::GetHull` on an empty edge map.** Upstream
  dereferences `edges.begin()` of an empty map (#290, the port returns an
  empty hull); only the replica is run there (111 residue records of the
  deviation case).
* **`IncrementalDelaunay2::GetContainingTriangle` with the default start**
  (`initialTriangle` left invalid). The walk then starts at triangle 0 of
  hash-table numbering; every query in the case starts at an explicit
  canonical triangle instead, which covers the same code.
* **`IncrementalDelaunay2::GetGraph`, `GetVertices`.** The graph object is
  not emitted; its triangles are compared through `GetTriangulation` (which
  copies `mVertices` and lists every graph triangle, supervertex triangles
  included), and `GetVertices` is the same vertex array.
* **`IncrementalDelaunay2` internals with preserved quirks.**
  `DoEarClipping` leaving `GetNumActive()` one too large and the unreachable
  `numPolygon == 2` test of `RetriangulateBoundaryRemovalPolygon` (#290) are
  private and do not reach any public output. The boundary branch of
  `Remove` is reached only through `FinalizeTriangulation`: a point strictly
  inside the rectangle is never adjacent to a supervertex.

## Upstream bug suspects

**1. `ExtremalQuery3BSP` is nondeterministic: its tree and its wrong
answers depend on heap addresses (new).** `CreateSphericalBisectors` sorts
the triangles around each vertex with `SortAdjacentTriangles`, which starts
at `*tAdj.begin()` of `Vertex::TAdjacent`, an
`std::unordered_set<Triangle*>` whose iteration order MSVC takes from the
pointer hash. The bisector arcs, and hence the tree, therefore depend on
where the allocator placed the triangles. Measured with the MSVC build on
200 integer shears of the octahedron, each constructed once and then rebuilt
8 times in the same process between unrelated allocations: every one of the
200 polytopes got a different node count on some rebuild (counts from 11 to
28 nodes), and 55 of them answered some directions differently from the
first construction (3871 of 320000 answers; two different answers for a
direction with a unique extreme vertex mean at least one of them is wrong).
Two runs of the oracle generator produced different goldens for a case that
filtered directions by upstream's own answer. Consequences: `GetNumNodes`,
`GetTreeDepth` and the wrong answers of #290 are not reproducible; a
program that builds the same query twice can get two different answers.
The port reads the same containers in sorted order and is deterministic.

**2. #290 is not confined to icosahedra (sharpened measurement).** The
finding records 0 wrong answers for the octahedron, also over 200 random
rotations. Under integer linear maps with positive determinant (shears,
which are not rotations) the MSVC build answers 142 of 28800 random
directions (0.49%) of the sheared octahedra with a vertex that is not the
exact argmax, while the unsheared octahedron (10800 directions), the
tetrahedron and the triangular bipyramid, sheared or not (80400), gave none.
The port, whose tree comes from sorted containers, answered all 120000
directions of this measurement correctly; that is luck of the tree shape,
not a fix, since the construction defect is preserved.

**3. `IncrementalDelaunay2::GetHull` returns garbage for collinear input
without looping (extends #290).** The finding says the walk "loops forever,
or writes past the output". It can also return to its start early, and then
upstream returns normally with `hull[]` padded by the zeros of the `resize`,
i.e. a "hull" that contains supervertex 0: domain `(-6,-8)-(6,6)`, points
`(4,4) (2,2) (0,0)`, `FinalizeTriangulation()`, `GetHull` returns
`[7, 8, 0]` (the MSVC build). In the deviation case's deep run this is the
outcome on 697 records and the unbounded walk on 454. The port now throws on
both (port defect 2 above).

No other new suspect was found. The independent checks above hold on every
deep-run record except the one #399 record, so on the compared inputs both
implementations return genuine minimum circles and spheres, genuine
antipodes, and genuine Delaunay triangulations.
