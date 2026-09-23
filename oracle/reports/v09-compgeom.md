# Verify group 9 (computational geometry) against the MSVC build of upstream GTE

Family `v09-compgeom`, 16 cases, 20 golden records each. Deep run
`npm run oracle:deep -- 2000 v09-compgeom` (32000 records, 6860 of them throw
records): **every case passes**. Wall time of generation and replay together:
3.7 s.

The group is `ConvexHull2.h`, `Delaunay2.h` and `Delaunay3.h`. The port
implements the *replacement* class templates `ConvexHull2<T>`,
`Delaunay2<T>` and `Delaunay3<T>`, which embed their own
SWInterval-filtered / BSNumber-exact predicates; the deprecated
`Delaunay2<InputType, ComputeType>` and `Delaunay3<InputType, ComputeType>`
specializations that route through `PrimalQuery2/3` are not ported
(porting-status.json, and the header comments of `src/Delaunay2.ts` and
`src/Delaunay3.ts`). The C++ side instantiates `ConvexHull2<double>`,
`Delaunay2<double>` and `Delaunay3<double>` with their own default rational
types, so both sides run the same filter and the same exact fallback.

Nothing on any compared path calls the C math library: the outputs are
combinatorial (dimensions, index lists, adjacency lists, duplicate maps,
hull lists, search paths) plus a handful of values computed with
`+ - * / sqrt` only (the dimension-1 line direction through `Normalize`, the
dimension-2 plane normal through `UnitCross`, and the plane constant through
`Dot`). **All 16 cases are declared `{ exact: true }`; every real output of
every non-deviation case is bit-identical to the MSVC value** (56000 real
outputs in the deep run, 56000 exact, no tolerance anywhere in the family).

**No port defect was found.** This branch contains no `src/` change.

## Canonicalization

`ETManifoldMesh` and `TSManifoldMesh` hold their features in
`std::unordered_map`, so upstream enumerates triangles and tetrahedra in
hash-table order while the port enumerates them in sorted-key order. The
*set* of features and each feature's stored vertex tuple are order
independent (`Insert` stores `V` in the rotation the caller passed and each
feature is inserted exactly once, so the tuple is a total key), but the
numbering is not. Every case therefore

* sorts the triangles / tetrahedra by their stored vertex tuple on both
  sides and emits them in that order, with the adjacency indices remapped
  through the same permutation;
* sorts the `GetHull()` edge pairs (2D) and face triples (3D)
  lexicographically, which upstream's own comment says are unordered;
* starts `GetContainingTriangle` / `GetContainingTetrahedron` at the
  canonically first feature on both sides, because upstream's default start
  is "whatever `mTMap.begin()` yields", and then emits the whole search path
  translated through the same permutation.

Everything else - the returned `bool`, the dimension, all counts, the
duplicate map, `finalV[]` (which holds vertex indices, not feature indices),
the line, the plane, the graph feature counts - is emitted verbatim. The
`ConvexHull2::GetHull()` index list is **not** canonicalized: its order is
the divide-and-conquer merge order and is part of the result.

## Coverage

| header | cases | comparison | deep-run result |
| --- | --- | --- | --- |
| `ConvexHull2.h` | `ConvexHull2.compute` (`operator()`, `GetDimension`, `GetLine`, `GetNumPoints`, `GetNumUniquePoints`, `GetHull`; 6 generator modes reaching hull dimension 0, 1 and 2) | exact | pass, 2000 records |
| `Delaunay2.h` | `Delaunay2.compute` (`operator()`, `GetDimension`, `GetNumVertices`, `GetNumUniqueVertices`, `GetNumTriangles`, `GetDuplicates`, `GetIndices`, `GetAdjacencies`, `GetIndices(t,·)`, `GetAdjacencies(t,·)` in range and out of range, `GetHull`, `GetLine`, `GetGraph`) | exact | pass, 2000 records |
| | `Delaunay2.compute.lowDimension` (dimension 0 and 1 results and their line) | exact | pass |
| | `Delaunay2.getContainingTriangle`, `…defaultStart` (`GetContainingTriangle` + the whole `SearchInfo`: `numPath`, `path`, `finalTriangle`, `finalV`, `initialTriangle`; query points inside a triangle, on an edge, at a vertex, far outside, uniform) | exact | pass, 2x2000 records |
| | `Delaunay2.degenerateAccessorsThrow` (throw parity for `GetHull` and `GetContainingTriangle` on a dimension 0/1 result) | exact | pass, 2000 throw records |
| | `Delaunay2.compute.deviation.numVertices` (#277), `Delaunay2.compute.deviation.epsilon` (#391) | deviation | pass, 2000/2000 records deviate each |
| `Delaunay3.h` | `Delaunay3.compute` (`operator()`, `GetDimension`, `GetNumVertices`, `GetNumUniqueVertices`, `GetNumTetrahedra`, `GetDuplicates`, `GetIndices`, `GetAdjacencies`, `GetIndices(t,·)`, `GetAdjacencies(t,·)` in range and out of range, `GetHull`, `GetLine`, `GetPlane`, `GetGraph`) | exact | pass, 2000 records |
| | `Delaunay3.compute.lowDimension` (dimension 0, 1 and 2 results with their line and plane) | exact | pass |
| | `Delaunay3.getContainingTetrahedron`, `…defaultStart` (`GetContainingTetrahedron` + the whole `SearchInfo`; query points inside, on a face, on an edge, at a vertex, far outside, uniform) | exact | pass, 2x2000 records |
| | `Delaunay3.degenerateAccessorsThrow` (throw parity for `GetHull` and `GetContainingTetrahedron` on a dimension 0/1/2 result) | exact | pass, 2000 throw records |
| | `Delaunay3.compute.deviation.duplicates` (#283), `…deviation.epsilon` (#391), `…deviation.numVertices` (#283) | deviation | pass, 2000/2000 records deviate each |

### Generators

Point sets have at most 10 points (the exact fallback is affordable and the
records stay small). The modes are selected by `io.index()`:

* **uniform** doubles in `[-4,4]`, the generic interval branch;
* **lattice** `[-3,3]`, exact arithmetic with frequent collinear triples;
* **cocircular / cospherical lattice**: subsets of the 12 lattice points on
  `x^2+y^2=25` and the 30 lattice points on `x^2+y^2+z^2=9`, translated by a
  lattice vector, mixed with the center and with nearby lattice points;
* **dense lattice** `[-2,2]`, where duplicates and degeneracies are common.

`ConvexHull2` adds an exactly collinear mode and an all-points-equal mode.
The degenerate cases build axis-aligned collinear sets and `z = c` coplanar
sets, which the intrinsics classify exactly.

The main `Delaunay2` / `Delaunay3` cases use a **capped rejection loop**
(64 attempts, every coordinate redrawn per attempt, fallback to a parabola /
moment-curve point set) around a probe built from upstream's own control
flow: construct `IntrinsicsVector2/3<double>` with upstream's hardcoded
`epsilon = 0`, require the reported dimension to be 2 (3), require the exact
orientation of the extreme triangle (tetrahedron) to be nonzero, and require
`info.extremeCCW` to agree with that exact sign. That is the exact separator
for finding #391: on an accepted input the port's seed simplex *is*
upstream's seed simplex and the two implementations run the same algorithm.
`Delaunay3` additionally rejects repeated vertices, which is the exact
separator for finding #283.

### Degeneracy actually reached (deep run, 2000 records per case)

| measurement | value |
| --- | --- |
| `ConvexHull2` hull dimension 0 / 1 / 2 | 499 / 437 / 1064 |
| `ConvexHull2` records with a duplicated input point | 948 |
| `Delaunay2` triangles per record | mean 4.92, from 1 to 14 |
| `Delaunay2` records with a duplicated input vertex | 707 |
| `Delaunay2` records with an exactly collinear triple | 760 |
| `Delaunay2` records with an exactly cocircular quadruple | 873 |
| `Delaunay2` predicate calls falling back to exact rational arithmetic | 3023 of 68379 (4.4%) |
| `Delaunay3` tetrahedra per record | mean 7.63, from 1 to 24 |
| `Delaunay3` records with an exactly coplanar quadruple | 961 |
| `Delaunay3` records with an exactly cospherical quintuple | 506 |
| `getContainingTriangle` inside / outside the hull | 1282 / 718 (2498 / 1502 over both cases) |
| `getContainingTetrahedron` inside / outside the hull | 1336 / 664 (2673 / 1327 over both cases) |

The fallback rate was measured by subclassing the port's `Delaunay2` and
re-running `toLine`'s and `toCircumcircle`'s interval prologue in the
override; an exactly zero determinant always leaves the interval straddling
zero, so every collinear / cocircular configuration in the table above does
reach the exact rational branch on both sides.

## Port defects fixed

None. Every association order, accumulation seed, division form and branch
in `src/ConvexHull2.ts`, `src/Delaunay2.ts` and `src/Delaunay3.ts` already
matches upstream bit for bit on all 32000 deep-run records, so this branch
carries no `src/` change.

## Deliberate deviations demonstrated

All five deviation cases deviate on **2000 of 2000** deep-run records.

1. **`Delaunay2.compute.deviation.numVertices` (#277).**
   `Delaunay2<T>::GetNumVertices()` returns `mIRVertices.size()`, and
   `operator()` clears `mIRVertices` and returns before resizing it when the
   intrinsic dimension is 0 or 1, so upstream reports 0 input vertices while
   `GetVertices()` still hands back the caller's non-empty array. The port
   returns the stored vertex count on every path. The generator is confined
   to the degenerate inputs upstream classifies *correctly* (all points
   equal, or collinear along a coordinate axis, verified against an exact
   dimension computation), so the accessor is the only difference:
   `operator()`'s return value and `GetDimension()` agree on every record.

2. **`Delaunay3.compute.deviation.numVertices` (#283).** The same defect in
   3D, for dimension 0, 1 and 2 input.

3. **`Delaunay2.compute.deviation.epsilon` (#391).** `Delaunay2<T>` builds
   `IntrinsicsVector2<T>` with a hardcoded `epsilon = 0`, and the intrinsics
   measure perpendicular distances against a *normalized* direction. An
   exactly collinear set whose direction is not axis aligned therefore has
   distances that are nonzero by an ulp, the strict `maxDistance <= 0` test
   fails, and the set is reported as 2-dimensional with an exactly
   degenerate seed triangle: upstream returns `true` with a degenerate
   triangle list, or throws from the incremental update (1411 of the 2000
   deep-run records are throw records). The port classifies the dimension
   with its own exact `toLine` predicate and returns dimension 1. The
   generator draws exactly collinear lattice sets and accepts only those on
   which upstream's own `IntrinsicsVector2` reports dimension 2, which is
   the exact separator: when it reports 1 the port agrees with it bit for
   bit, and that regime is the `Delaunay2.compute.lowDimension` case.

4. **`Delaunay3.compute.deviation.epsilon` (#391).** The 3D half, on exactly
   coplanar lattice sets `u*(b,-a,0) + v*(c,0,-a)`; 1449 of 2000 deep-run
   records are throw records.

5. **`Delaunay3.compute.deviation.duplicates` (#283).**
   `Delaunay3<T>::ProcessedVertex` hashes *and* compares its `location`
   member and the lookup key is built as `ProcessedVertex(mVertices[i], i)`,
   so a repeated coordinate triple at a later index carries a different
   location and the find always misses. `mDuplicates[]` stays the identity
   map, `GetNumUniqueVertices()` equals the input count, and `Update(i)` runs
   a second time for a point that is already a mesh vertex. The port hashes
   and compares the vertex only, as the 2D class correctly does. The
   generator forces one repeated vertex into a point set whose seed
   tetrahedron is sound, so duplicate detection is the only defect on these
   inputs.

The confinement of the port's #391 fix is checked by the main cases, not
only asserted: `Delaunay2.compute`, `Delaunay3.compute`, both containment
searches and both `lowDimension` cases run the broad generators and agree
bit for bit on all 16000 records, so the fix changes nothing where upstream
is sound.

## Independent reference checks (agreement is not correctness)

The deep-run inputs of the main cases were re-checked against exact
predicates evaluated with the port's own `BSNumber` (arbitrary precision,
`+ - *` only), independently of both implementations. Every check passed on
every record:

* **ConvexHull2** (1064 dimension-2 records): every consecutive triple of
  the returned hull turns strictly counterclockwise (so the polygon is
  convex and has no collinear or repeated vertex), and no input point is
  strictly right of any hull edge.
* **Delaunay2** (2000 records): every output triangle is exactly
  counterclockwise; no input vertex is strictly inside any output
  triangle's circumcircle (the Delaunay property, exact in-circle
  determinant); the adjacency across face `j` shares the reversed edge and
  is symmetric; the exact sum of the triangle areas equals the exact area of
  the convex hull of the same points, so the triangulation covers the hull
  with no gap and no overlap; and every edge of `GetHull()` has every input
  point on its left or on it.
* **Delaunay3** (2000 records): every output tetrahedron is exactly
  positively oriented; no input vertex is strictly inside any circumsphere
  (exact in-sphere determinant); each adjacency shares exactly one face; and
  every triangle of `GetHull()` has every input point on its inner side, so
  the returned faces are outward oriented and bound the point set.
* **containment searches** (4 x 2000 records): whenever a simplex index is
  returned, the exact orientation predicates put the query point inside or
  on that simplex; the rest return "outside".

No violation was found, so the oracle's agreement on these cases is
agreement about a result that is independently correct.

## Sensitivity of the cases

Almost every output of this family is combinatorial and is decided by an
*exact* predicate, so it cannot be moved by a change of association order in
floating point. What the cases can and cannot discriminate was measured
rather than assumed, by mutating `src/` and re-running:

* **seed-simplex orientation.** Flipping the sign test that orients the seed
  triangle in `src/Delaunay2.ts`
  (`toLineSign < 0 ? [e0, e1, e2] : [e0, e2, e1]`) makes
  `Delaunay2.compute` disagree on 20 of 20 golden records and both
  containment cases on 19-20 of 20. The combinatorial outputs are therefore
  fully load bearing.
* **interval-filter association (cannot be discriminated).** Regrouping the
  interval accumulation of `toCircumcircle`
  (`x0c0.add(x1c1).add(x2c2)` to `x0c0.add(x1c1.add(x2c2))`) changes no
  output on any of the 32000 deep-run records, as it must: the SWInterval
  filter is only a conservative accelerator and the exact rational branch
  decides every sign the filter leaves open. No case in this family - and no
  case that could be written against the public API - can discriminate the
  association order inside the interval filter, nor the association order
  inside the rational expression trees, which is exact.
* **the one non-combinatorial output.** `ConvexHull2`'s dimension-1
  direction is `Normalize(p1 - p0)`, which `Vector.h` writes as
  `v /= length`, i.e. a multiplication by `1/length`. A componentwise
  `v[i] / length` would have produced a different double on **67 of the 437
  dimension-1 deep-run records**, so that case pins the division form. The
  same `Normalize` and the `UnitCross` / `Dot` of the dimension-2 plane are
  pinned by `Delaunay2.compute.lowDimension` and
  `Delaunay3.compute.lowDimension`.

## Not covered

* **`ConvexHull2<T>::operator()(int32_t, Vector2<T> const*)`,
  `Delaunay2/3<T>::operator()(size_t, Vector<T> const*)`,
  `GetPoints()` / `GetVertices()`.** The pointer-plus-count overloads are
  subsumed by the array overload in the port (a documented porting choice),
  and `GetPoints` / `GetVertices` return the caller's own array unchanged.
  Nothing is computed on these paths.
* **Input sets larger than 32 points.** Above MSVC's insertion-sort
  threshold `std::sort` is no longer stable, and `ConvexHull2::operator()`
  sorts indices by an equivalence-inducing comparator, so upstream's hull
  index list stops being reproducible; see the suspect below. The
  generators are capped at 10 points (which also keeps the exact fallback
  affordable and the golden records small).
* **`ConvexHull2::GetTangent`'s expired-bounding-loop path (#277,
  preserved).** The loop is bounded by `size0 + size1` and, when it expires,
  the function silently returns whatever indices it last held. The port
  preserves that. It was never observed to fire: over 2000 deep-run records
  including all-equal, collinear, duplicate-heavy and cocircular inputs, the
  hull was verified convex and enclosing on every dimension-2 record, which
  it would not be if the merge had been fed non-tangent indices. No
  generator that reaches it is known.
* **`Delaunay2/3::GetHull`'s "there must be at least one triangle
  (tetrahedron)" `LogError`.** Unreachable: `mDimension == 2` (3) implies a
  nonempty mesh, whose boundary always has an unshared face.
* **`GetBarycentrics`.** It lives in `Delaunay2Mesh.h` / `Delaunay3Mesh.h`,
  which belong to verify group 11, not to this group's three headers.
  `Delaunay2.h` / `Delaunay3.h` themselves expose no barycentric entry
  point.
* **`GetHullMesh`.** No such member exists in `ConvexHull2.h`.
* **The deprecated `Delaunay2<InputType, ComputeType>` and
  `Delaunay3<InputType, ComputeType>` specializations.** Not ported
  (porting-status.json); upstream states they will be removed.

## Upstream bug suspects

**`ConvexHull2<T>::operator()` reports an unspecified index for a duplicated
hull vertex (new).** The hull is built from
`std::sort(mHull, lessThanPoints)` followed by
`std::unique(…, equalPoints)`, where `lessThanPoints` compares the *points*,
not the indices. Repeated points are therefore equivalent under the
comparator, `std::sort` leaves equivalent elements in an unspecified order,
and which index of a repeated point survives `std::unique` - and hence which
index `GetHull()` reports for that hull vertex - is unspecified. The hull
*geometry* is unaffected.

It is invisible for at most 32 points, because MSVC's `std::sort` falls back
to an insertion sort there and is stable by accident; that is the only
reason the port's stable `Array.prototype.sort` agrees with the C++ build on
this family's committed cases. Above the threshold the two disagree. A
scratch case with 40 lattice points in `[-3,3]^2` disagreed on 11 of 20
records, for instance on

```
(-1,-1) (0,2) (2,0) (-1,-2) (3,-3) (3,0) (-2,3) (3,-1) (-1,-2) (1,-1)
(0,-2) (3,1) (-1,1) (2,-1) (2,0) (-1,-1) (1,1) (-2,1) (0,1) (3,-1)
(3,3) (-2,3) (-1,3) (-2,-3) (-2,-1) (3,-2) (-2,-1) (3,-3) (0,0) (3,1)
(0,-3) (-3,1) (-1,0) (-1,-1) (-1,-1) (-2,2) (-2,-2) (-2,-3) (-2,1) (0,1)
```

where the MSVC build reports hull index 27 and the port reports 4, and
`points[4] == points[27] == (3,-3)`. This is unspecified behaviour rather
than a wrong answer, and there is no "upstream value" to match (it is
MSVC's introsort partitioning), so the port is left alone and the
generators are capped below the threshold. Same class as the v03 finding on
`SymmetricEigensolver3x3::ComputePermutation`. A caller who needs a
deterministic index should be told to deduplicate the input first.

No other new suspect was found: the independent exact checks above hold on
every deep-run record of every main case, so on the inputs where upstream's
dimension classification and duplicate detection are sound, both
implementations produce a genuinely correct Delaunay triangulation,
tetrahedralization and convex hull.
