# Group 43 (primitives) — differential test against the MSVC build

Family `v43-primitives`, 48 cases: 39 ordinary `exact`, 7 `deviation` (also
`exact`), 2 with the default scaled tolerance. Deep run
(`npm run oracle:deep -- 2000 v43-primitives`): 96 000 records, 5 550 of them
C++ throws, **all 49 tests pass**, 5.1 s wall time.

The case file makes `Matrix2x2.h` and `Matrix3x3.h` visible, so
`Hyperellipsoid::FromCoefficients`'s `Inverse(A, &invertible)` resolves to the
closed-form overloads for `N = 2` and `N = 3`, which is the configuration the
port assumes (v34, issue #217 item 4). `Cone.h` and `ContOrientedBox3.h` pull
`Matrix3x3.h` in regardless.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `Hyperellipsoid.h` | `getM.2d/3d` (GetM + GetMInverse), `toCoefficients.2d/3d` (both overloads), `fromCoefficients.2d/3d`, `fromCoefficientsABC.2d/3d`, `defaultConstruct`, `compare.2d`, `fromCoefficients.decoupledDeviation.3d` | exact | pass |
| `Hyperplane.h` | `construct.3d` (default, normal+constant, normal+origin), `fromPoints.3d`, `compare.3d`, `fromPoints.deviation.2d/4d` | exact | pass |
| `Cone.h` | `construct` (all four constructors), `setAngle`, `heights` (all four cone types, `HeightInRange` / `HeightLessThanMin` / `HeightGreaterThanMax` / `IsFinite` / `IsInfinite`), `heights.throwParity`, `compare`, `createMesh` | `construct`, `heights`, `heights.throwParity`, `compare` exact; `setAngle` and `createMesh` tolerance (`cos`, `sin`, `tan`) | pass |
| `Polygon2.h` | `queries` (5 generator modes), `queries.clockwise` | exact | pass |
| `Tetrahedron3.h` | `normals` (face, edge, vertex), `computeCentroid`, `getPlanes`, `tables`, `compare`, `defaultConstruct`, `getPlanes.deviation` | exact | pass |
| `RectangleManager.h` | `initialize`, `update` (two move passes) | exact | pass |
| `RectangleMesh.h` | `construct`, `tcoords`, `frame.deviation` | exact | pass |
| `RectanglePatchMesh.h` | `construct`, `frame` | exact | pass |
| `AlignedBoxBV.h` | `queries` (default box `[-1,1]^3`, `GetSplittingAxis`, the three static intersection queries) | exact | pass |
| `OrientedBoxBV.h` | `queries` | exact | pass |
| `AlignedBoxTreeOfPoints.h` | `create` (structure, partition, node boxes, the three `Execute` query types) | exact | pass |
| `AlignedBoxTreeOfSegments.h` | `create` | exact | pass |
| `AlignedBoxTreeOfTriangles.h` | `create` (+ the `Execute` intersection lists) | exact | pass |
| `OrientedBoxTreeOfPoints.h` | `create` | exact | pass |
| `OrientedBoxTreeOfSegments.h` | `create` | exact | pass |
| `OrientedBoxTreeOfTriangles.h` | `create`, `leafExtent.deviation` | exact | pass |

Generators mix uniform draws, small lattices, exactly representable signed
permutation frames, collinear and coplanar point sets, all-equal points,
duplicate primitives, zero-volume tetrahedra, self-intersecting and
collinear-vertex polygons, polygons whose constructor fails on a duplicate
index, touching rectangles that share endpoint values, and meshes of 2x2 up to
5x5 samples. Tree heights include the "build the whole tree" sentinel and
explicit heights 0..4, so the early-stop branch of `BuildTree` is covered.

## Tolerances

Only the two cases whose upstream path calls the C math library carry a
tolerance; both use the harness default of `1e-12` and both are far inside it
over the 2 000-record deep run.

| case | libm call | measured max scaled error | bit-identical outputs |
| --- | --- | --- | --- |
| `Cone.setAngle` | `std::cos`, `std::sin`, `std::tan` in `SetAngle` | 2.932e-16 | 13 629 / 14 000 |
| `Cone.createMesh` | the same three, through `GenerateInscribed` / `GenerateCircumscribed` and `tanAngle` | 8.882e-16 | 70 738 / 71 940 |

`Cone.setAngle` emits `angle` itself with `outRealExact`, and `createMesh`
emits `numExtra` (through the vertex count), the index array and the unique
vertex count as integers, so a `tan` rounding difference that moved
`ceil(tNumExtra)`, or a `cos`/`sin` difference that changed which vertices
`UniqueVerticesSimplices` considers equal, would fail loudly rather than hide
under the tolerance. Neither happened in 2 000 records.

## Port defects fixed

None. Every arithmetic-only case agreed bit for bit with the MSVC build on the
first replay, including the oriented-box fits (covariance plus
`SymmetricEigensolver3x3`), the tree partitions and the mesh channels.

## Deliberate deviations demonstrated

| case | issue | what upstream does |
| --- | --- | --- |
| `Hyperplane.fromPoints.deviation.2d` | #217 | `Hyperplane<2,T>(std::array<Vector2,2>)` builds `SingularValueDecomposition(2, 1, 32)`, whose `LogAssert(mNumCols >= 2)` fires. All 2 000 deep records throw; the port computes the complement with `ComputeOrthogonalComplement` for `Vector2`. |
| `Hyperplane.fromPoints.deviation.4d` | #217 | `svd.Solve(&edge[0], -1)` trips `LogAssert(multiplier > 0)`. All 2 000 deep records throw; the port uses the documented default multiplier. |
| `Hyperellipsoid.fromCoefficients.decoupledDeviation.3d` | #80 | `SymmetricEigensolver::Tridiagonalize` stores the reflection parameter `2/Dot(v,v) == 2` for a *degenerate* Householder step (the subcolumn below the subdiagonal is already zero, so the reflection actually applied was the identity). `GetEigenvectors` rebuilds `H = I - 2*e*e^T`, which flips the sign of one row of `Q`, so the recovered axes are not the axes of the quadric. An axis-aligned ellipsoid gives a diagonal `M`, which is exactly that case. |
| `Tetrahedron3.getPlanes.deviation` | #268 | `GetPlanes` assigns `normal` and `constant` but never `origin`, so the caller's array keeps whatever was there (the case pre-fills it, as a caller reusing a buffer would). The port returns planes whose three members are consistent; `normal` and `constant` are compared in the main case and agree bit for bit. |
| `RectangleMesh.frame.deviation` | #268 | `InitializeFrame` hardcodes `tangent = (1,0,0)`, `bitangent = (0,1,0)` for every vertex regardless of the rectangle's axes. The port uses the rectangle's own orthonormal axes. The `normal` channel is unaffected and is compared in `RectangleMesh.construct`. |
| `OrientedBoxTreeOfTriangles.leafExtent.deviation` | #343 | The smallest-extent scan of `ComputeLeafBoundingVolume` ends with `absExtent > minAbsExtent`, so the *largest* extent is zeroed and the leaf box collapses along its longest axis. `ApprGaussian3` returns extents in increasing order, so upstream selects index 2 on essentially every triangle. |
| `BVTreeOfTriangles.coincident.deviation` | #167 | `Execute` collects hits in a `std::set<Intersection>` ordered by `parameter` alone, so two triangles hit at the same parameter are set-equivalent and all but one are dropped. Two triangles with the same three vertex indices are hit at bit-identical parameters. The port orders by `(parameter, triangleIndex)`. |

Deviating records over the 2 000-record deep run: `#217` 2 000/2 000 (both
cases), `#80` 2 000/2 000, `#268` (planes) 2 000/2 000, `#268` (mesh frame)
1 914/2 000, `#343` 1 649/2 000, `#167` 828/2 000. The residues are explained:
the mesh-frame case agrees when the rectangle's axes happen to be exactly
`(1,0,0)` and `(0,1,0)` (one of the eight signed permutation frames the
generator draws); `#343` agrees when `extent[2]` really is the smallest, so
both sides zero the same index; `#167` agrees when the line misses the
duplicated triangle, in which case both sides report no hits.

Both defects that leak into an ordinary case are excluded by an exact predicate
evaluated on upstream's own control flow:

* `Hyperellipsoid.fromCoefficients.3d` and `.fromCoefficientsABC.3d` reject
  inputs whose `M` decouples, using `EigenDecouples`, a verbatim copy of
  upstream's `Tridiagonalize` (the v03 precedent). `N = 2` needs no rejection:
  `Tridiagonalize`'s loop runs `mSize - 2` times, so it is empty.
* `AlignedBoxTreeOfTriangles.create` and `OrientedBoxTreeOfTriangles.create`
  emit an upstream-computed `coincident` flag (a verbatim probe that runs the
  same `FIQuery` over *every* triangle, an exact superset of the reported
  hits) and emit the intersection list only when it is false.
* `OrientedBoxTreeOfTriangles.create` omits the extents of leaf nodes; the
  centre, the axes and every interior node's extents are still compared.

## Independent reference checks

Run once over the 2 000-record deep run (the port and the MSVC build agree bit
for bit on all these cases, so checking the recorded outputs checks both). The
checker is plain JavaScript and uses no gtengine-js code.

* **Tree containment.** For every node of every tree, every primitive vertex
  in `[minIndex, maxIndex]` lies inside that node's bounding volume — aligned
  boxes exactly, oriented boxes to `1e-9 * scale`. Node ranges nest, so this
  also proves containment in every ancestor. 2 000 records each for the point
  and triangle trees, aligned and oriented.
* **Tree queries vs brute force.** A leaf every one of whose ancestors' boxes
  is hit with a margin is reported; a leaf with a definitely-missed ancestor is
  not. (The comparison needs the margin because GTE's line/box test is a
  separating-axis test and the reference is the slab method; they disagree only
  at the rounding boundary.) For the triangle trees, every triangle the line
  crosses well inside its interior (Möller–Trumbore, `u, w > 0.02`,
  `u + w < 0.98`, clearly non-parallel) is among the reported intersections,
  and every reported point lies on the line and in its triangle's plane.
  452 and 441 of the 2 000 records produce line hits.
* **Polygon2 vs exact rational arithmetic.** On the 1 188 lattice records,
  `ComputeArea` equals the exact BigInt shoelace area, and `IsConvex` agrees
  with the exact BigInt sign test of the turns in both directions (for
  `numIndices > 3`; a 3-index polygon short-circuits to "convex" upstream
  without looking at the turns).
* **Tetrahedron planes.** On the 1 983 non-degenerate records, each
  `GetPlanes` plane contains all three of its face's vertices to `1e-9 * scale`
  and the opposite vertex is on the negative side.
* **Hyperellipsoid round trip.** On all 2 000 valid records,
  `FromCoefficients(ToCoefficients(E))` reproduces `M = sum U_d U_d^T / e_d^2`
  and the centre to `1e-9`.
* **Rectangle-mesh frames.** On all 2 000 records, every normal is unit
  length, orthogonal to both rectangle axes and equal to
  `axis[0] x axis[1]`; every position lies in the rectangle's plane and within
  its extents.

## Not covered

* `BVTree::SplitPoints` with **more than 32 primitives in a node range.** The
  port replaces `std::nth_element` with a full stable sort. A direct MSVC probe
  (2 000 trials per size, distinct and tied projections) shows that MSVC's
  `std::nth_element` leaves the range in exactly stable-sorted order for
  `n <= 32` and differs from it on 97–100 % of inputs for `n >= 33`; MSVC's
  `std::sort` behaves the same way. The cases here keep at most 8 primitives,
  so the agreement observed in the deep run is *the accident of the
  insertion-sort path*, not a property of the port. A tree built over more than
  32 primitives in one node would have a different `mPartition` order, which
  feeds the covariance sums of the oriented-box fits, and on tied projections
  even a different left/right *membership*. Reproducing MSVC's introselect
  exactly would hard-code one standard library's unspecified order into the
  port; it is recorded here instead. The same reasoning applies to
  `RectangleManager::Initialize`'s `std::sort` over `2n` endpoints (the cases
  use at most 8 rectangles, so 16 endpoints).
* `Cone::CreateMesh` is exercised with `numMinVertices = 3` only, to keep the
  committed records small.
* `Mesh::Update` / `RectanglePatchMesh::UpdatePositions`, `UpdateNormals` and
  `UpdateFrame` re-run the same `Initialize*` routines the construction cases
  already compare, so they are covered indirectly only.
* `Polygon2`'s empty-vertex-pool guard (the port's `#268` assert) cannot be
  reached from C++ at all: upstream has only a pointer and cannot tell an empty
  pool from a populated one, so there is nothing to compare.

## Upstream bug suspects

New findings, not already in `docs/UPSTREAM-FINDINGS.md`:

1. **`BVTree::SplitPoints` depends on the unspecified order of
   `std::nth_element`.** The left child receives `info[0 .. medianIndex]` and
   the right child receives the rest *in reverse*, both in whatever order
   `nth_element` happened to leave them. The order within each partition is
   unspecified by the standard, and on tied projections — which the class's own
   use cases produce routinely (duplicate primitives, lattice data, a
   splitting axis orthogonal to the spread) — even the left/right *membership*
   is unspecified. Two conforming standard libraries therefore build different
   trees, with different node bounding volumes, from the same input. This is
   the same class of finding as `SymmetricEigensolver::ComputePermutation`'s
   unstable `std::sort` (issue #478), one level up. `std::stable_sort` on
   `(projection, centroidIndex)` would make it deterministic. Measured above
   with an MSVC probe: identical to a stable sort for `n <= 32`, different for
   `n >= 33`.
2. **`Polygon2::IsSimple` and `IsConvex` return `true` for a degenerate
   triangle.** Both short-circuit on `numIndices == 3` without looking at the
   vertices, so three collinear or coincident points are reported as a simple
   convex polygon. Minor: the documented precondition is that the caller
   supplies a simple polygon. Observed on the lattice records of
   `Polygon2.queries`.
3. **`RectangleManager::Initialize`'s y-overlap test is inconsistent with the
   `TIQuery` used by `Update`.** The sweep tests
   `r0.max[1] >= r1.min[1] && r0.min[1] <= r1.max[1]` on the y-interval only
   and relies on the x-sweep for the x-interval, while the incremental update
   calls `TIQuery<AlignedBox2, AlignedBox2>`. The two agree for the closed
   intervals used here, so no disagreement was observed; it is noted only
   because the initial and incremental paths do not share a predicate.
