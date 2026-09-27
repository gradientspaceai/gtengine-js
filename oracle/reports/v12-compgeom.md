# Verify group 12 (computational geometry) against the MSVC build of upstream GTE

Family `v12-compgeom`, 18 cases, 20 golden records each.
`oracle/cpp/cases/v12-compgeom.cpp`, `test/oracle/v12-compgeom.oracle.test.ts`.
Deep run `npm run oracle:deep -- 2000 v12-compgeom` (36000 records, 4000 of
them throw records): **every case passes** (19 tests including the coverage
check; generation about 5 minutes, replay 89 s). Two deep generations are
byte-identical.

The group is `MinimumVolumeBox3FloatingPoint.h` and `MinimumVolumeBox3Rational.h`,
the two implementation headers behind the `MinimumVolumeBox3.h` wrapper that
v08 already drives through all four `operator()` overloads. This group audits
and closes v08's gaps: v08 compared a point cloud only by dimension, flags and
a volume with tolerance 2e-2 (floating point) and 1.5e-1 (rational), and a
fixed mesh only when upstream's answer was stable under reordering.

**All 12 non-deviation cases are declared `{ exact: true }` and every real
output is bit-identical to the MSVC value** (508000 real outputs in the deep
run, 508000 exact). **No tolerance anywhere in the family**, and v08's two
tolerance cases are superseded by exact ones (`compute.canonical` of each
header, full box bit for bit).

One port defect was found and fixed: the v08 confinement of the #405/#426
fixes was still too broad (see "Port defects fixed").

## What decides the answer: order, not libm

Neither header calls the C math library. The only function is `std::sqrt`, on
doubles and on `BSRational` (GTE's `std::sqrt(BSRational)` converts to double
first, as the port does). Everything else is `+ - * /` on doubles or exact
`BSNumber` / `BSRational` arithmetic. So v08's tolerance was entirely
order-decided, by three unspecified orders:

1. `ExtractMeshTopology` numbers the edges and triangles by iterating
   `ETManifoldMesh`'s `std::unordered_map`s. `mEdgeIndices` (every pair
   `e0 < e1`), the first-minimum tie break and the edge visiting order of
   `RemoveCoplanarTriangleAdjacencies` follow that numbering, and
   `ProcessEdgePair` is not symmetric in its two edges (#530).
2. `ExtractVertexAdjacencies` copies each vertex's `VAdjacent`, an
   `std::unordered_set<int32_t>`, into the adjacency pool; the pool order
   decides `GetExtreme`'s climb path and which of several tied vertices it
   returns.
3. For a point cloud, `ConvexHull3::GetHull()` lists the hull triangles in
   hash order (and, per v13, in an order that depends on heap addresses);
   that order is the mesh insertion order, which fixes each edge's
   `V[0]/V[1]` (the vertex `ComputeVolume` assumes to be the minimum) and
   `T[0]/T[1]` (the `N`/`M` roles).

The port enumerates all three in sorted-key order.

### The control

`Canonical<Base>` in the case file runs upstream's own protected pipeline
step by step, exactly as `operator()` does, with two steps replaced by copies
that differ only in iterating sorted copies of the hash containers
(`ExtractMeshTopologySorted`, `ExtractVertexAdjacenciesSorted`), and for a
point cloud sorts the `ConvexHull3` triangles (which are `TriangleKey<true>`
tuples; `UniqueVerticesSimplices` renumbers the used vertices monotonically,
so a lexicographic sort is the port's `std::map` order). Everything else - the
convex hull, the exact geometry, the coplanar merge, the aligned candidate,
the level-curve processors, the four minimizers, `ComputeVolume`,
`GetExtreme`, the `std::thread` candidate search and the rational box - is
upstream's code. The replay calls the port's public `compute` /
`computeHull` unchanged. With the order pinned, the full box, and on the hull
cases every observable internal (edges, adjacency lists after the coplanar
merge, climb start, sample subdivision, aligned candidate, every field of the
winning candidate), compare bit for bit on all 2000 records.

How much the order matters on the MSVC build (the hull cases record
`RawAgreement` between upstream's raw `operator()` and the canonical run as a
diagnostic input; deep run, 2000 records each):

| header | same box bit for bit | same box up to axis signs/order | same volume, different box | different volume |
| --- | ---: | ---: | ---: | ---: |
| floating point, `computeHull` | 1107 | 345 | 93 | 455 |
| rational, `computeHull` | 939 | 353 | 367 | 341 |

So on 17 to 23 percent of small hull meshes the reported minimum volume itself
depends on the hash order, in the rational pipeline too (different edge roles
sample different candidates; the minimum over the samples then differs). A
single run of the point-cloud path, where ConvexHull3's order adds to it,
gave a different volume on 752 of 2000 floating-point and 317 of 2000
rational clouds; that path is not recorded as a diagnostic because the
heap-address order would make the golden irreproducible. This quantifies #530.

## Coverage

| header | cases | comparison | deep run (2000 records/case) |
| --- | --- | --- | --- |
| `MinimumVolumeBox3FloatingPoint.h` | `compute.canonical` (overloads 1/2 path, all 5 cloud modes, `numThreads` 0-4) | exact | pass |
| | `computeHull.canonical` (overloads 3/4 path, scrambled triangle order and rotation, box + all internals, raw overload 3 or 4 run for `RawAgreement`) | exact | pass |
| | `getExtreme` (`GetExtreme` via the subclass, 4 directions of 4 kinds per record) | exact | pass |
| | `minimizers` (the four virtual minimizers overridden: call counts, processors that reach them, summed s/t arguments, any subset replaced by a no-op) | exact | pass |
| | `compute.dimension2` (raw overloads 1/2) | exact | pass |
| | `invalidArgument` (throw parity, all `LogAssert`s of the four overloads) | exact | pass, 2000 throw records |
| | `getExtreme.plateau` (#426), `computeHull.provenViolation` (#405/#426), `compute.dimension2.floatComputeType` (design), `computeHull.reuse` (new suspect) | deviation | pass; 1969, 1698, 1995, 1994 of 2000 deviate |
| `MinimumVolumeBox3Rational.h` | `compute.canonical`, `computeHull.canonical`, `getExtreme` (exact; rotated clouds included), `minimizers` (`MinimizerVariableT` always replaced), `compute.dimension2`, `invalidArgument` | exact | pass (2000 throw records for `invalidArgument`) |
| | `compute.dimension2.floatComputeType` (design), `computeHull.reuse` (new suspect) | deviation | pass; 1996, 1989 of 2000 deviate |
| `MinimumVolumeBox3.h` (wrapper) | covered by v08 (`compute`, `compute.lowDimension`, `computeHull` and the #352/#355/#405/#426 deviations for both headers) | v08 | v08 |

Dimension 0 and 1 of both headers: v08's `compute.lowDimension` cases already
compare the full box exactly on both headers; not duplicated. Dimension 2:
new here (v08 had it only as the #352/#355 deviation).

Every public entry point the port implements is covered: the four
`operator()` overloads (C++ side: overloads 1/2 in `compute.dimension2` and
`invalidArgument`, 3/4 in `computeHull.canonical`, `invalidArgument` and
v08; the canonical cloud and hull runs replicate their bodies), the
constructor's thread count (`numThreads` 0-4 on the floating-point cases,
0-3 on the rational ones: upstream's `std::thread` candidate search, the
port's single thread, bit-identical), and the protected customization surface
a subclass sees: `GetExtreme`, the four minimizers and the `Candidate` /
topology state (`mEdges`, `mEdgeIndices`, adjacency pool, `mVClimbStart`,
`mMaxSample`, `mDomainIndex`, `mAlignedCandidate`, `mMinimumVolumeObject`).
`ComputeVolume`, `ProcessEdgePair`, `Pair` and the level-curve processors are
not virtual; they are observable through the winning candidate's fields
(`levelCurveProcessorIndex`, `edgeIndex`, `N`, `M`, `f00..f11`, axes, support
indices, volume) and through the minimizer overrides, which see
`c.levelCurveProcessorIndex`: 66 distinct processor indices reach a
minimizer in each header over the deep run, all four minimizers are reached
(floating point: `VariableS` on 2000 records, `VariableT` 1560, `ConstantS`
743, `ConstantT` 752).

### Generators

Cloud modes (`DrawCloud`): 0 lattice `[-3,3]^3` (exact ties, coplanar faces,
collinear hull edges: both branches of `RemoveCoplanarTriangleAdjacencies`);
1 uniform `[-5,5]^3`; 2 a lattice cloud rotated by 0.03 to 0.3 degrees (the
#426 regime; the sin/cos is applied in the generator, only coordinates are
recorded); 3 a random subset of a `{0,1,2} x {0..b} x {0..c}` grid with the
corners forced (face-interior and edge-interior hull points); 4 a lattice
cloud mapped by the integer matrix `5R`, `R` the rotation by `(3/5, 4/5)`
(tilted minimum box, exact arithmetic, exact ties between edge pairs). v08
excluded lattice clouds from the rational case because their exact ties are
decided by the hash order; with the order pinned they are all included.
Hull meshes for the vertices-and-indices overloads come from upstream's
`ConvexHull3` + `UniqueVerticesSimplices` (as the cloud query builds them),
with the triangles sorted (for reproducibility) and then shuffled and each
triple rotated. Directions for `getExtreme`: signed axes, lattice directions,
uniform unit vectors and normalized face normals (a whole face is then a
plateau up to rounding). Dimension 2 (`DrawPlanarCloud`): lattice planes
perpendicular to an axis, tilted lattice planes `o + a*d0 + b*d1`, uniform
in-plane coordinates, and in-plane Pythagorean rotations.

### Restrictions of the exact cases and their separators

* **Floating point, #405 and #426.** Accepted only when all six support
  vertices of upstream's winning candidate are exact extremes along its three
  axes over the translated vertices the climb uses (`SupportsAreExact`,
  exact `BSNumber` dot products). That is exactly the condition under which
  the port's (now confined) fixes cannot touch the winner; it also proves
  that upstream's rational box contains the hull. Capped at 24 attempts;
  1993 and 1995 of 2000 records pass, the rest are capped fallbacks, which
  agree too. `getExtreme` requires upstream's vertex to be an exact maximizer
  per direction (16 redraws).
* **Rational, #355.** Records whose canonical run reaches
  `MinimizerVariableT` are rejected (a probe subclass counts the calls;
  capped at 8, fallback the fewest calls). The `minimizers` case instead
  always replaces `MinimizerVariableT` by a recorder: its callers and its
  (double-rounded) arguments are then compared exactly.
* **Dimension 2 (#352/#355 and the MinimumAreaBox2 compute type).**
  `Dimension2Probe` runs upstream's hull, both Newell loops and both
  `MinimumAreaBox2` instantiations: accepted when upstream's basis equals the
  corrected one bit for bit and `MinimumAreaBox2<double, double>`'s rectangle
  equals `MinimumAreaBox2<double, BSRational<UIntegerAP32>>`'s (the exact
  instantiation v11 showed the port reproduces). 1660 / 1634 of 2000 draws
  pass; the rest use a fixed tilted-lattice fallback that passes both probes
  (recorded as a diagnostic). Note on #352: for four or more hull vertices the
  unclosed loop sums the Newell normal of the hull polygon with `hull[0]`
  removed, which has the right direction; only its magnitude differs, and
  after normalization the basis is often bit-identical (always for a plane
  perpendicular to an axis). The loop is wrong only in exact terms for those;
  for a triangular hull the normal is exactly zero.

### Independent reference checks (agreement is not correctness)

* Exact containment: on every accepted floating-point record the winner's
  support vertices are exact extremes (above), so upstream's rational box -
  and the port's, which is bit-identical - contains the hull exactly. The
  rational pipeline is exact by construction.
* Every `compute.canonical` / `computeHull.canonical` record of both headers
  emits two booleans computed on each side from its own output: the
  (double) box contains every input point to 1e-9 relative, and the axes are
  orthonormal and the volume equals `8 e0 e1 e2`, both to 1e-9. All 8000
  records: true on both sides.
* The `getExtreme.plateau` and `provenViolation` diagnostics are the exact
  shortfall of upstream's vertex against the exact extreme; the records on
  which the port differs are exactly those with a positive excess over the
  rounding bound (1969 = 1969, 1698 = 1698), so the port's two fixes fire
  precisely when upstream is provably not extreme and never otherwise.

## Port defects fixed

**`src/MinimumVolumeBox3FloatingPoint.ts`: the v08 confinement of the #405
and #426 fixes was still too broad.** Checked here as asked, on the first
golden it failed: 8 of 20 `compute.canonical` and 11 of 20
`computeHull.canonical` records disagreed with upstream.

* What differed: the box centre, extents and volume by 1 to 4 ulp (scaled
  error up to 3.3e-16), and on those records the winner's `minSupportIndex`
  (for example `[2, 1]` against upstream's `[1, 5]`).
* Root cause: v08 made `computeVolume` replace upstream's assumed minimum
  vertex (the hull-edge vertex, #405) whenever the hill climb found a
  *strictly smaller double* projection, and `getExtreme`'s plateau traversal
  (#426) took any vertex with a *strictly larger double* dot product. When a
  candidate axis is the rounded normal of a hull face, the other vertices of
  that face project within an ulp of the edge vertex, so on ordinary uniform
  clouds both fired on rounding noise. The support index feeds the exact
  rational `GetMinimumVolumeBox`, where the two vertices have different exact
  projections, so the box moved. v08's lattice-only probe could not see it
  (lattice face normals give exact ties).
* Fix: both replacements now require the new vertex to be beyond
  `climbTolerance(direction) = 8 eps max|d_i| max_v L1(v)`, a rigorous bound on
  the difference of the rounding errors of two three-term dot products, which
  proves that upstream's vertex is not extreme in exact arithmetic;
  otherwise upstream's vertex and value are returned. `getExtreme` keeps the
  traversal but returns upstream's climb result unless the traversal's best
  exceeds it by more than the bound.
* After the first narrowing the deep run showed 1 + 4 of 2000 remaining
  disagreements; all were records on which a support of upstream's winner is
  short of the exact extreme by more than the bound (a proven, if tiny,
  #405 violation, 5e-15 relative in volume), which is the deliberate fix
  acting. The main cases now select on that exact separator and
  `computeHull.provenViolation` demonstrates it.
* Regression tests (`test/MinimumVolumeBox3FloatingPoint.test.ts`): an
  upstream-verbatim subclass (`GetExtreme` and `ComputeVolume`) is compared
  with the port on 150 uniform clouds: bit-identical whenever no support
  comparison is violated beyond the bound, with noise-level violations
  present and some proven-violation clouds changed; and the oracle record
  itself. Setting the bound to zero (the v08 guard) fails both tests. v08's
  confinement test still passes; its "different box" count fell from 118 to
  4 of 300 lattice clouds (the visibly broken ones).

**`src/MinimumVolumeBox3FloatingPoint.ts`, `src/MinimumVolumeBox3Rational.ts`:**
comment only, documenting that `createMeshTopology` resets `mEdgeIndices`
(the port's behaviour on functor reuse, see the suspect below).

## Deliberate deviations demonstrated

| case | finding | deviating records (deep) | residue |
| --- | --- | ---: | --- |
| `MinimumVolumeBox3FloatingPoint.getExtreme.plateau` | #426 | 1969 / 2000 | the 31 records whose capped search found no shortfall beyond the bound (diagnostic 0); one-for-one |
| `MinimumVolumeBox3FloatingPoint.computeHull.provenViolation` | #405, #426 | 1698 / 2000 | the 302 records whose search found none; one-for-one. 45 of the deviating boxes still contain the points to 1e-9: proven but tiny shortfalls |
| `MinimumVolumeBox3FloatingPoint.compute.dimension2.floatComputeType` | design (port's `MinimumAreaBox2` is exact-only) | 1995 / 2000 | the 5 records whose search found no differing lifted box |
| `MinimumVolumeBox3Rational.compute.dimension2.floatComputeType` | design | 1996 / 2000 | the 4 such records |
| `MinimumVolumeBox3FloatingPoint.computeHull.reuse` | new suspect | 1994 / 2000 | the 6 records where no stale pair reaches a minimizer and the box is unchanged |
| `MinimumVolumeBox3Rational.computeHull.reuse` | new suspect | 1989 / 2000 | the 11 such records |

The plateau construction follows #426: a lattice box whose *first* point lies
inside a face (it becomes hull vertex 0, the climb start), rotated slightly
so the exact coplanar merge leaves it in the graph; directions are the signed
normals of the triangles at vertex 0 and the axes of upstream's own winner.
The `floatComputeType` deviations select on the observable (the rectangle
lifted into 3D with upstream's expressions), because a difference in the sign
of a zero rectangle component alone vanishes in the lift. v08's
`compute.coplanar` (#352/#355), `compute.nonContaining` (#405/#426) and
`compute.variableT` (#355) deviations are not duplicated.

## Not covered

* **`MinimumVolumeBox3GPU`**: not ported.
* **Upstream's raw point-cloud box** (overloads 1/2 with the hash and
  heap-address orders): not comparable by construction; the canonical run
  replaces it and the table above measures the difference on the MSVC build.
* **Reuse with a second mesh that has fewer edges**: upstream then indexes
  `mEdges` out of range (undefined behaviour); not exercised.
* **Minimizer overrides with `numThreads > 0`**: upstream calls the virtual
  minimizers concurrently from its `std::thread`s, so a counting override is
  a data race; the `minimizers` cases run single-threaded.
* **Rational records that reach `MinimizerVariableT` with upstream's body**:
  #355, v08's `compute.variableT` deviation.

## Upstream bug suspects

**1. `MinimumVolumeBox3FloatingPoint.h`, `MinimumVolumeBox3Rational.h`: a
reused functor processes the previous mesh's edge pairs (new, minor; UB when
the new mesh is smaller).** `CreateMeshTopology` resizes `mEdges` but only
`reserve`s `mEdgeIndices`, and `ExtractMeshTopology` appends every pair
`(e0, e1)`. Calling `operator()` a second time on the same object therefore
processes the first mesh's pairs, read as indices into the second mesh's
edges, before the second mesh's own. When the second mesh has at least as
many edges these are duplicates: every minimizer they reach is called again
(1994 of 2000 floating-point and 1989 of 2000 rational records in the
`reuse` cases, observed through the minimizer overrides) and, when
candidates tie exactly, the stale order can change which one wins (5 and 2
records changed the box). When the second mesh has fewer edges the stale
indices are out of range: undefined behaviour. Fix: `mEdgeIndices.clear()`
in `CreateMeshTopology`. Port: resets the list (a fresh functor's
behaviour); now documented in both files.

**2. `MinimumVolumeBox3FloatingPoint.h`, `MinimumVolumeBox3Rational.h`,
dimension 2: `MinimumAreaBox2<T, T>` (minor, design).** The dimension-2 path
instantiates `MinimumAreaBox2` with the floating-point compute type, for
which `MinimumAreaBox2.h` itself says a correct result is not guaranteed. The
lifted box differs from the exact rectangle's on most coplanar inputs (1995
of 2000 deviation records drawn), usually in the last bits. Not a new defect
class (v11 `MinimumAreaBox2.deviation.floatComputeType`); recorded because the
MVB3 port inherits the exact computation.

**3. #530, quantified.** On the MSVC build the reported minimum *volume*
depends on the `std::unordered_map`/`unordered_set` iteration order on 455 of
2000 floating-point and 341 of 2000 rational hull meshes (table above); v08
measured it only under triangle rotations.
