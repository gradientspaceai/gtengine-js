# Verify group 27 (imaging) against the MSVC build of upstream GTE

Family `v27-imaging`, one header, `AdaptiveSkeletonClimbing3.h` (2905
lines), 10 cases, 20 golden records each (golden file 551 KB). Deep run
`npm run oracle:deep -- 2000 v27-imaging` (20 000 records, 1492 of them
throw records, all in the throw-parity case): **every case passes** (12
tests: the 10 cases, the independent-reference test and the coverage check;
generation and replay together about 50 s). Outside the deviation case the
deep run compares 657 330 floating-point outputs, **all bit-identical** to
the MSVC build; every other output is an integer, a boolean or a 32-bit
FNV-1a digest (vertex bit patterns, triangle indices, boxes) and compares
exactly. The header is arithmetic-only (`+ - * /`, `std::sqrt` in
`ComputeNormals`, `std::floor` in the face subdivision, truncating
conversions in `GetGradient`; no libm call), so every case is declared
`{ exact: true }` and there is no tolerance anywhere in the family. The
goldens and a 2000-record deep file were each generated twice,
byte-identical.

**One upstream defect was found and fixed in the port**: the 3-D instance of
#544. The four-crossing face cases of the six `Get{X,Y,Z}{Min,Max}EdgesS`
pair the face's level-set points by the level-free `det = f00*f11 - f01*f10`
with the two disjoint pairings swapped, exactly as `AdaptiveSkeletonClimbing2`
does. It is demonstrated by a `deviation` case (2000 of 2000 records) and a
C++-side independent check that finds upstream's own mesh wrong on 62 200 of
86 605 four-crossing faces of that case. No translation defect was found.
Six further upstream suspects (cracks between merged boxes, `depth > N`,
`fixBoundary`, in-face ear-clipping diagonals, holes at plus-sign branch
points, `OrientTriangles` consistency) are preserved and described below;
they were found by the independent checks, on which both sides agree bit
for bit.

## Coverage

| header | case | comparison | deep run |
| --- | --- | --- | --- |
| `AdaptiveSkeletonClimbing3.h` | `.extract` (int32, N = 1, 2, one or two extractions) | exact | 76 854 reals + digests, 100 % |
| | `.extract.types` (int8, uint8, int16, uint16, uint32 < 2^31) | exact | 64 755 reals + digests, 100 % |
| | `.extract.saddle` (sound four-crossing faces, plus signs at level 0) | exact | 163 317 reals, 100 % |
| | `.extract.saddlePairing` | deviation (#544 in 3-D) | 2000 of 2000 records deviate |
| | `.extract.large` (9^3, 17^3, 33^3 from a recorded formula) | exact | 49 455 reals + digests, 100 % |
| | `.extract.closedLarge` (level sets inside 9^3, 17^3) | exact | 94 758 reals + digests, 100 % |
| | `.extract.closed` (level sets inside 3^3, 5^3) | exact | 118 722 reals + digests, 100 % |
| | `.extract.rootMonobox` (#194, preserved) | exact | discrete only (empty meshes), 100 % |
| | `.meshOps` (`MakeUnique`, `OrientTriangles`, `ComputeNormals` on arbitrary meshes) | exact | 88 416 reals, 100 % |
| | `.invalid` (N = 0, null voxels) | exact, throw parity | 1492 records throw on both sides; 1053 reals |

Every extraction emits the monobox list (upstream reads it back through the
public `PrintBoxes`, the replay through the port's private `mBoxes`), the raw
vertices and triangles of `Extract`, the vertices and triangles after
`MakeUnique`, the triangles after `OrientTriangles` (random `sameDir`), the
`ComputeNormals` result and the three counts of the C++-side face check
(below). Lists of up to 24 elements are emitted in full, longer ones as a
digest after their count. Upstream's order is deterministic (the octree
recursion of `Merge`, `Tessellate`'s box order, the ear-clipping order,
`MakeUnique`'s first-occurrence numbering through `std::map`), so nothing is
re-sorted: the verbatim order is a stronger comparison than a canonical sort.
Upstream is called through `__declspec(noinline)` wrappers.

### How the cases reach the branches

* **Images.** Random small integers (ties, saddles, isolated voxels), linear
  ramps (strongly monotone: the whole image merges), paraboloid bowls and
  domes (spheres), constants, signed checkerboards (saddle-rich),
  full-range values, a constant with one or two isolated voxels, quadrics
  with cross terms (hyperboloids: saddles, tunnels), two spheres; for the
  large cases also octahedra/boxes (`|u|` forms), tori (tunnels) and cubics.
  `.closed` surrounds random or checkerboard interiors by a constant border
  so that the level set never reaches the border.
* **Levels.** Half-integers (the documented use), reals, integers of the
  voxel range (a voxel on the level: upstream does not assert in 3-D,
  `SetLevel` counts the voxel as below), 0 and -0, outside the range,
  quarter-integers, NaN and both infinities, and for `.saddle` levels equal
  to `det / S` (the face saddle on the level).
* **Depth -2..N+1 and `fixBoundary`** (one record in four), one or two
  extractions on the same object.
* **Face pairing.** The main cases are restricted by the exact predicate
  `sign(det) = -sign(dg)` (below; `BSNumber` on the C++ side, capped
  rejection with a saddle-free fallback), so they compare upstream's own
  answer wherever it is sound; `.saddle` puts at least one sound
  four-crossing face in every record, `.saddlePairing` at least one
  defective one. Record 0 of both is a pinned reproduction (the image of
  suspect 1).
* **Branch histogram** of the deep-run replay (temporary counters in the
  port, removed): x/y/z merges 2.7M/5.3M/10.7M; unit-box face types 0,
  3, 5, 6, 9, 10, 12, 15 = 780k, 262k, 262k, 264k, 263k, 261k, 291k, 159k;
  merged-box face types 0, 3, 5, 6, 9, 10, 12 = 164k, 48k, 47k, 46k, 46k,
  46k, 44k (no merged face ever has four crossings, so its `LogError` and the
  unit-face `LogError` for odd crossing counts are unreachable: a corner is
  either above or not, so a unit face always has 0, 2 or 4 crossings); face
  subdivision 13 870 times, endpoint swap taken 4114 times; fan
  triangulation `RemoveTrianglesSE` 11 180 boxes, ear clipping of merged
  boxes 62 066; face pairing (upstream sign, sign of dg): (+,-) 26 588,
  (-,+) 25 947, (0,0) 2453 sound, and (+,+) 47 877, (-,-) 47 971, (0,+)
  2517, (0,-) 2965, (-,0) 1143, (+,0) 1083 defective.
* **`meshOps`** feeds `MakeUnique` duplicates with flipped zero signs
  (`std::map` equivalence merges +0 and -0, the first survives), rotations
  and duplicates of triangles (kept apart), degenerate triangles and unused
  vertices (zero normals), and vertices on a quarter lattice in
  [-1.5, 2^N + 1], so `GetGradient`'s truncating conversion takes the
  -1 < x < 0 branch and the zero gradient outside [0, 2^N) is reached.

### Independent references (replay; deep-run counts)

Checked on every extraction whose level is finite and no voxel value (the
documented precondition): 14 966 extractions.

* **Level set, exact.** Every vertex of the unique mesh is (a) on a unit grid
  edge whose end values straddle the level, within `(a + 4) 2^-53` of the
  exact root `a + (L - f0)/(f1 - f0)` (BigInt; 587 443 vertices), (b) a
  plus-sign branch point on a four-crossing face whose bilinear saddle value
  is exactly the level, at (u of the v = 0 crossing, v of the u = 0
  crossing) (1666), or (c) a fan centroid strictly inside a merged box
  (9084). No exception.
* **No degenerate triangle** (exact cross product): 0 of all.
* **Face saddles (the #544 check).** For every four-crossing unit face the
  mesh must carry the two cuts of the interpolant's pairing (`dg > 0`: the
  cuts around corners 10 and 01; `dg < 0`: around 00 and 11) and no branch
  point inside the face, or, where `dg = 0`, a branch point adjacent to all
  four crossings. A face with all four cuts is ambiguous (the ear clipping of
  a box put diagonals into the face that complete the other pairing,
  suspect 5). The port: 0 wrong of 98 208 faces in all cases (2750
  ambiguous). The same check runs in C++ on upstream's own mesh and is
  emitted: 0 wrong in every sound case, and in `.saddlePairing` upstream is
  wrong on 62 200 of 86 605 faces (2266 ambiguous; 1984 of 2000 records
  have a wrong face, the other 16 have their defective faces only with a
  corner on the level, which the check skips, or ambiguous).
* **Closed meshes.** When no border voxel is on the other side of the level
  than the rest of the border, the raw triangles (both copies of coincident
  triangles counted) must use every edge an even number of times. 3344
  such meshes: 3198 closed; the 146 open ones are all explained: 65 cracks
  where every border edge lies in a grid plane whose boxes on the two sides
  are different rectangles (suspect 2), 11 at plus-sign branch points
  (suspect 6), 70 at depth > N (suspect 3). Meshes from unit boxes only
  (depth = N or `fixBoundary`) without branch points are always closed.
  8 meshes have an edge on four triangles and 4 closed ones are opened by
  `MakeUnique` (suspect 5).
* **Depth > N.** In all 1606 such extractions without `fixBoundary`, every
  box is a unit box with a four-crossing face (suspect 3).
* **#194.** All 2000 `.rootMonobox` extractions have a nonempty level set and
  an empty mesh with no box, on both sides.
* **`GetZeroBase` -1 (#194, latent).** The replay counted the
  `Get{X,Y,Z}Interp` calls at a negative coordinate: 0 in the deep run.
* **`OrientTriangles`.** Of 3190 closed 2-manifold unique meshes, 1052 are
  not consistently oriented afterwards (6838 directed edges used twice,
  282 012 triangles); orienting with the gradient of the cell that contains
  each triangle instead of the cell at `trunc(vertex)` leaves 169. Reported
  as an observation (suspect 7), not asserted.

## Port defects fixed

None of translation. The one change to `src/` fixes an upstream defect
(suspect 1): `facePairing` in `src/AdaptiveSkeletonClimbing3.ts` computes
`det` as before (BigInt, matching `int64_t`), then the exact sign of
`dg = det - L*(f00 + f11 - f01 - f10)` (level written as `m 2^e`, the helper
of `AdaptiveSkeletonClimbing2`), and keeps upstream's choice exactly when
`sign(det) = -sign(dg)`, otherwise uses `-sign(dg)`. It changes nothing on
sound faces (the main cases agree bit for bit on 657 330 reals) and the
pairing of every defective face. Regression tests in
`test/AdaptiveSkeletonClimbing3.test.ts`: the pinned image at levels 0.5 and
2.5, a `dg = 0` face with `det = -1` (branch point), exactness at voxel
values near 1e9 (`det` is not representable as a difference of rounded
double products), and the existing plus-sign test moved to level 0 (where
`det = dg = 0`). The first oracle generation of `.invalid` drew images
`a(x + y + z) + b yz` unfiltered and 9 of 2000 deep records deviated: the
face 0, 1; 1, 0 of `a = 1, b = -2` has its saddle exactly at level 0.5
(`det = -1`, `dg = 0`), where the port inserts the plus sign. That was a
harness error (the generator was not restricted to sound faces); it now is.

**Sensitivity.** Recomputed from the deep-run outputs: `Get*Interp` as
`x + (L - f0) * (1 / (f1 - f0))` would change 51 139 of 587 443 edge
vertices (6451 of 14 966 extractions); `ComputeNormals` dividing by a
multiplication with `1 / length` would change the normals of 7824
extractions.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `AdaptiveSkeletonClimbing3.extract.saddlePairing` | #544 (3-D face cases), this report, suspect 1 | 2000 of 2000 records: each has a four-crossing face where upstream's pairing contradicts the interpolant (exact predicate, rejection capped at 256, fixed defective fallback). The C++ face check finds upstream wrong on 62 200 faces; the port's check finds 0. |

Preserved findings reached and agreeing bit for bit: #194 (root monobox:
2000 of 2000 `.rootMonobox` records empty on both sides; the latent
`GetZeroBase` -1 read never occurs; the `static_cast<float>` of the integer
box corners in `GetVertices` is exact below 2^24 and so unobservable for any
image that fits in memory).

## Not covered

| header / entry point | reason |
| --- | --- |
| `AdaptiveSkeletonClimbing3<T, float>` | the port has one number type (`Real` = double); only `Real = double` is compared. |
| `uint32_t` voxels at or above 2^31 | `det = f00*f11 - f01*f10` overflows `int64_t` (undefined behaviour, as #544's minor item for the 2-D class); the port computes with BigInt. |
| `N < 0`, huge `N` | `1 << N` is undefined behaviour; memory. `N = 0` is the throw-parity case. |
| levels equal to a voxel value | compared bit for bit (upstream does not assert in 3-D), but excluded from the independent geometric checks: the documented precondition is a non-integer level, and crossings then sit on voxels. |
| `PrintBoxes` output format | not ported (debugging `ostream` output); the box list it prints is compared through a parse. |
| `GetGradient` (through `OrientTriangles`) at coordinates outside `int32_t` or NaN | the conversion to `int32_t` is undefined behaviour. |
| `MakeUnique` / `OrientTriangles` / `ComputeNormals` with out-of-range indices or NaN coordinates | out-of-range reads; NaN keys break `std::map`'s strict weak ordering. |

## Upstream bug suspects

**1. The four-crossing face is paired by the level-free `det`, with the
pairings swapped (#544 in 3-D; result-corrupting; fixed in the port).** The
six `Get{X,Y,Z}{Min,Max}EdgesS`, case 15, read the face corners as f00, f10
(first face coordinate + 1), f01 (second + 1), f11 (lines 1503-1545,
1596-1638, 1689-1731, 1782-1824, 1875-1917, 1968-2012) and join the two
crossings next to corner 00 and the two next to corner 11 when
`det = f00*f11 - f01*f10 > 0`, which cuts those corners off; `det < 0` cuts
off 10 and 01; `det = 0` inserts the plus-sign branch point. The trilinear
interpolant restricted to the face is bilinear, and its topology is decided
by `dg = (f00-L)(f11-L) - (f01-L)(f10-L) = det - L*(f00+f11-f01-f10)`:
`dg > 0` connects 00 and 11 through the face, so 10 and 01 must be cut off;
`dg = 0` is the plus sign. Upstream is right exactly when
`sign(det) = -sign(dg)`; at level 0 it is wrong on every non-degenerate
saddle face, and its plus sign is right only at level 0. Reproduction: the
3^3 image that is -1 everywhere except 4 at (0,0,0), (1,1,0), (2,2,0), level
0.5, depth 1: the face z = 0 of the voxel at the origin has corners 4, -1;
-1, 4, `det = 15`, `S = 10`, saddle value 1.5 > 0.5, so the corners of value
4 are connected through the face; upstream's mesh (record 0 of
`.saddlePairing`) carries the cuts (0.7,0,0)-(0,0.7,0) and (1,0.3,0)-(0.3,1,0)
around the corners of value 4 and neither cut around a corner of value -1.
At level 2.5 (saddle below the level) upstream is right. Measured: upstream
wrong on 62 200 of 86 605 faces of the deviation case by its own mesh. Port:
`facePairing` decides by the exact sign of `dg` and keeps upstream's choice
wherever it is right.

**2. Merged boxes crack where their common faces are cut into different
rectangles (result-corrupting; preserved).** A merged box subdivides its
face polylines (`Get*EdgesM`) only at the zero sub-edges of its own
merge-tree nodes (`IsZeroEdge || HasZeroSubedge` along the face's interior
grid lines). When the boxes on the two sides of a grid plane cut it into
different rectangles, an edge of one side runs across the face of the other
side and is not a sub-edge of it, so the two sides triangulate the plane
with different polylines and the mesh has holes. Reproduction: the 5^3 image
with border -9 and interior (z = 1..3, y = 1..3, x = 1..3)
`8 -6 -3 / 6 5 -8 / 3 5 8`, `8 -4 -4 / 5 -4 -3 / 1 6 0`,
`3 -8 4 / -5 -6 1 / 1 7 1`, level -7.5, depth 1: below z = 2 the boxes are
x 2..4 by y 0..1 and 1..2, above they are x 2..3 and 3..4 by y 0..2; the
lower side draws (2,0.3,2)-(3.7,1,2)-(3.75,2,2), the upper side
(2,0.3,2)-(3,0.3,2)-(3.75,2,2), and the four edges are on one triangle each.
Deep run: 65 of 3344 interior level sets (whose meshes should be closed)
crack this way; every open edge of every one is on such a plane (checked).

**3. `depth > N` drops every mergeable leaf (result-corrupting; preserved;
the #194 mechanism one level down).** `Merge` adds a merged child box on
behalf of its parent, and only at nodes with `depth <= 1`. With a user depth
D the parents of the leaves have depth `D - N + 1`; for `D > N` that is at
least 2, so a leaf that returns "mergeable" is neither merged nor added and
only leaves with a four-crossing face (which add themselves) survive.
Reproduction: `round(9 - |p - (4,4,4)|^2)` on 9^3, level 0.5: depth 3 gives
512 boxes and the sphere, depth 4 gives 0 boxes and an empty mesh. Deep run:
in all 1606 extractions with `depth > N` and no `fixBoundary` the surviving
boxes are exactly such leaves. `depth` is undocumented; `N` is the full
resolution.

**4. `fixBoundary` fixes every voxel, not only the boundary voxels (doc,
minor; preserved).** The constructor comment says image boundary voxels are
not allowed to merge; the leaf branch of `Merge` adds a unit box for every
leaf when `mFixBoundary` is set, so nothing merges anywhere and the depth
has no effect (same image: 512 boxes at depths -1, 0 and 4).

**5. The ear clipping puts diagonals and triangles into voxel faces (minor
to result-corrupting; preserved).** `RemoveTrianglesEC` clips the
lowest-index degree-2 vertex regardless of geometry, and when the two
neighbours of the clipped vertex lie on one face the diagonal lies in that
face. Seen three ways: four-crossing faces carrying all four cuts (26 of
7334 sound faces of `.saddle`, 245 of 2537 of `.closed`; in the pinned image
above the port's box gets a flat quadrilateral on the face), meshes with an
edge on four triangles where the boxes on both sides chose the same
in-face diagonal (8), and flat triangles that both boxes put into their
common face, which `MakeUnique` then keeps once, opening a closed mesh (4).

**6. Plus-sign branch points leave holes (result-corrupting; preserved).**
A branch point is a vertex of degree 4 in the box wireframe; the ear
clipping only clips degree-2 vertices and the triangulation around the
branch point is incomplete. Reproduction: 5^3, border -3, interior -1 except
the z = 2 layer `(x + y)` even 1 / odd -1 (every face of that layer has
`det = dg = 0` at level 0), level 0, depth 2: 24 edges are on an odd number
of triangles (the branch edges on one, face cuts on three). Deep run: 11
interior level sets. The port's fix inserts the plus sign also where
`det != 0` and `dg = 0` (the interpolant's level set is the plus sign there),
so it inherits this behaviour there.

**7. `OrientTriangles` does not orient a closed mesh consistently
(observation, minor; preserved).** Each triangle is flipped by the sign of
its normal against the average `GetGradient` of its vertices, and
`GetGradient` evaluates the trilinear gradient of the cell at
`trunc(position)`, which for a vertex on a cell face or edge is a
neighbouring cell's interpolant, whose normal derivatives differ. 1052 of
3190 closed 2-manifold meshes of the deep run come out inconsistently
oriented; with the gradient of the cell containing each triangle 169 do.
