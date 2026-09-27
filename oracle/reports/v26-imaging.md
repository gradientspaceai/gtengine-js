# Verify group 26 (imaging) against the MSVC build of upstream GTE

Family `v26-imaging`, 55 cases, 20 golden records each (golden file 636 KB).
Deep run `npm run oracle:deep -- 2000 v26-imaging` (110 000 records, 7439 of
them throw records): **every case passes** (59 tests: the 55 cases, the
coverage check and three independent-reference tests; generation and replay
together about 95 s; 2000 of the throw records are the C++ side of the
`DrawEllipse` deviation). Outside the seven deviation cases the deep run
compares 819 821 floating-point outputs (packed pixel words included), **all
bit-identical** to the MSVC build; every
other output is an integer, a boolean, a packed pixel word or a 32-bit
digest and compares exactly. The only C math function on any path is
`sqrt`/`sqrtf` (IEEE-exact), so every case is declared `{ exact: true }`
and there is no tolerance anywhere in the family. The goldens are
deterministic (regenerated, byte-identical).

Headers: `ImageUtility2.h`, `ImageUtility3.h`, `SurfaceExtractorCubes.h`,
`SurfaceExtractorMC.h`, `SurfaceExtractorTetrahedra.h` (the inherited
`SurfaceExtractor.h` members `MakeUnique`, `Convert`, `OrientTriangles` and
`ComputeNormals` are exercised through the Cubes and Tetrahedra extractors;
the base class itself is group 25).

**Two upstream defects were found by this group**: `SurfaceExtractorCubes`
pairs the crossings of every saddle face the wrong way round (fixed in the
port), and `SurfaceExtractorMC::ComputeNormals` never advances its triangle
pointer (the port already iterated correctly; the difference is now
documented, pinned and demonstrated). **Three translation defects were
fixed** (below), one earlier finding is refuted (#443 item 2: the marching-cubes table *is*
face-consistent), and one upstream limitation is documented
(`SurfaceExtractorCubes`' geometry-free ear clipping).

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `ImageUtility2.h` | `.getComponents` (N = 4, 8, caller offsets), `.dilate`, `.erode`, `.open`, `.close` (each N = 4, 8, caller offsets), `.extractBoundary`, `.floodFill4.int`, `.floodFill4.double`, `.getL1Distance`, `.getL2Distance`, `.getL2Distance.large`, `.getSkeleton`, `.drawThickPixel`, `.drawLine`, `.drawCircle`, `.drawRectangle`, `.drawEllipse`, `.drawFloodFill4` | exact | 190 478 reals (packed words included) + integers, 100 %; 373 throw records (empty offset lists) |
| | `.drawEllipse.zeroExtents` | deviation (#443) | 2000 of 2000 records deviate |
| `ImageUtility3.h` | `.getComponents` (N = 6, 18, 26, caller offsets), `.dilate`, `.erode`, `.open`, `.close` (each N = 6, 18, 26, caller offsets), `.computeCDConvex`, `.floodFill6.int`, `.floodFill6.double`, `.drawLine` | exact | 90 014 reals + integers, 100 %; 276 throw records |
| | `.morphology.xmin` | deviation (#129) | 317 of 2000 records deviate (residue below) |
| `SurfaceExtractorCubes.h` | `.extract`, `.extract.types`, `.extract.large`, `.extract.plusSign`, `.orientTriangles.points`, `.orientTriangles.thirds`, `.constructor.invalid` | exact | 246 480 reals + digests, 100 %; 1756 throw records |
| | `.extract.saddle` | deviation (this report, suspect 1) | 2000 of 2000 records deviate |
| `SurfaceExtractorMC.h` | `.extractVoxel`, `.extractVoxel.levelTopology`, `.extractVoxel.signedZero`, `.extract`, `.extract.large`, `.extract.levelTopology`, `.makeUnique`, `.orientTriangles.points`, `.orientTriangles.thirds` | exact | 184 213 reals + integers, 100 %; 1200 throw records |
| | `.extractVoxel.level`, `.extract.level` | deviation (#443 item 1) | 1352 and 1172 of 2000 deviate (residue below) |
| | `.computeNormals.multiple` | deviation (this report, suspect 2) | 1432 of 2000 deviate (residue below) |
| `SurfaceExtractorTetrahedra.h` | `.extract`, `.extract.types`, `.extract.large`, `.orientTriangles.points`, `.constructor.invalid` | exact | 108 636 reals + digests, 100 %; 1779 throw records |
| | `.orientTriangles.centralTetra` | deviation (#132) | 2000 of 2000 records deviate |

### How the cases reach the branches

* **Images** are drawn in eight modes (uniform density, rectangles, discs,
  one solid square of side 1..10, a ring, checkerboards and stripes, empty,
  full; in 3D boxes, balls, a checkerboard, rods), with 2s sprinkled into the
  morphology inputs (the `== 1` / `== 0` tests), and are recorded packed
  (several pixels per recorded double). The zero-border precondition of
  `GetComponents`, `ExtractBoundary`, `GetSkeleton` is imposed where the
  header requires it; the morphology, fill and distance cases also run on
  images without it (`GetL1Distance` then treats border pixels as distance
  1). Caller-specified neighbour lists include duplicates, the zero offset,
  offsets of reach 2 and the empty list (the `LogAssert` of 2D
  `Dilate`/`Erode`, 3D `Erode`/`GetComponents`: throw parity).
* **Label and index outputs** are emitted in full, packed: the component
  lists (upstream fills them in increasing index order; the labels are
  assigned in raster order of the first pixel, so the DFS stack order cannot
  reach an output), the relabelled image, the boundary walk, every fill.
  `ExtractBoundary` runs one to three walks on the same image (later walks
  see the pixels earlier ones marked 2), from starts inside components, in
  holes and on the border; a replica with a step cap would reject starts on
  which the walk does not return to its first pixel, and in a census of
  130 118 starts on 3000 images it never did.
* **Fills** use `int32_t` and `double` pixels; the double palette
  `{0, 1, -0, 2.5, NaN}` fills across `-0 == +0`, a NaN background fills only
  the seed, and the foreground is redrawn until it does not compare equal to
  the background (equal colours never terminate). Seeds one outside the image
  take the early return.
* **Distance transforms:** `.getL2Distance.large` uses 20..70 x 20..70
  images with 1..40 background pixels, so the extended neighbourhoods behind
  the thresholds 13^2, 31^2, 49^2, 72^2 run while every distance stays below
  the documented 100; all-foreground images give `sqrtf(float(INT32_MAX))`.
* **Draw functions** are fed through a recording callback on both sides and
  compared call by call (the sequence, not the set): short sequences in
  full, long ones (lines of 10 000 pixels, circles of radius 1000, solid
  circles of radius 150, ellipses with extents up to 100) as count plus
  digest. Modes aim at axis-aligned, diagonal, single-point and
  `|dx| = |dy| +- 1` lines, ties between the two largest 3D components (the
  strict `>` selecting the step axis), negative radii and extents, inverted
  rectangles. Extents stay at most 100, far inside the `int32_t` range of
  the decision terms (`4 * a^2 * (1 - y)` is of order `4 a^2 b`).
* **SurfaceExtractorCubes/Tetrahedra:** images of 2..3 voxels per axis
  (full emission when the deduplicated mesh has at most 16 triangles),
  4..6 per axis (digests), and every documented voxel type (`int8_t` to
  `uint32_t`, 32-bit values within 2^20 so that upstream's `int64_t`
  products and the port's doubles are both exact); values small lattice,
  uniform, linear, a quadric, binary and a hyperbolic product; levels from
  one below the minimum to the maximum (Tetrahedra: 3204 of the 6000 deep
  extraction records have samples on the level, i.e. the zero-corner cases of
  `ProcessTetrahedron`). Each record emits the rational extraction, its
  `MakeUnique`, the real extraction with and without duplicate removal,
  `OrientTriangles` (one bit per triangle) and `ComputeNormals`. The Cubes
  main cases reject images with a saddle face of nonzero determinant (exact
  integer predicate, at most 64 redraws, then a linear image, which has
  none; see suspect 1); `.extract.plusSign` plants faces with
  `f00*f11 == f01*f10` (8003 plus-sign faces in the deep run).
  `.orientTriangles.points` drives `GetGradient` on arbitrary points (inside,
  just outside, slightly negative, integer and quarter-integer points).
* **SurfaceExtractorMC** (`T = double`, `IndexType = int32_t`): per-voxel
  `Extract` with values uniform, lattice (corners on the level), eighths and
  signed zeros, perturbations `0`, `-0`, `+-1/1024`, `+-1e-300` and uniform;
  whole-image extraction on the same values plus a quadric and a hyperbolic
  product (ambiguous faces, 1586 deep records); `MakeUnique` on arbitrary
  vertex lists with duplicates and `-0` against `+0` (the first kept) and on
  the `UniqueVerticesSimplices` preconditions (no vertices, index counts 0, 1,
  2, an index out of range: 1200 throw records). The per-voxel `Extract`,
  whose outputs carry signed zeros, is called through a
  `__declspec(noinline)` wrapper.

### Independent references (run on the port's outputs of every replayed record)

* **Image utilities:** components partition the 1-pixels, the image holds
  the labels, each component is connected and no two are adjacent (4/8/6/18/
  26); every morphology output equals a direct definition written from the
  output pixel's side, and for binary inputs with symmetric neighbourhoods
  `Erode ⊆ A ⊆ Dilate`, `Open ⊆ A`, `A ⊆ Close` (outside `zeroExterior`);
  boundary walks start at the first nonzero pixel after the start and take
  8-neighbour steps over foreground pixels; fills equal the 4/6-connected
  region of the background colour grown from the seed; `GetL1Distance`
  equals the city-block distance to the nearest background pixel or the
  exterior, and its maximum and location follow the grass-fire loop;
  `GetL2Distance` equals the exact Euclidean distance transform rounded to
  float on every deep record (in `.large` the maximum distance exceeds 13,
  31, 49 and 72 on 1272, 371, 116 and 6 of 2000 records, up to 87.7: the
  "exact below 100" claim holds there); skeletons are subsets of the input.
  Lines: endpoints, one pixel per major-axis step, every minor coordinate
  within 1/2 of the ideal line (exact integers). Circles: every pixel within
  one unit of the radius, the 8-fold symmetry, solid circles contain every
  pixel at distance at most `r - 1`. Ellipses: the 3x3 block of every pixel
  straddles `b^2 X^2 + a^2 Y^2 = a^2 b^2`, 4-fold symmetry. Rectangles and
  thick pixels: exact pixel sets. `ComputeCDConvex`: a voxel is 0 exactly when
  it and everything beyond it in some axis direction is 0.
* **Skeleton (#443, preserved):** of the deep records whose foreground is
  one solid square, every even side vanished (side 2: 76 of 76, 4: 31, 6: 11,
  8: 6, 10: 3) and no odd side did (1: 139, 3: 29, 5: 13, 7: 8, 9: 5),
  identically on both sides.
* **Extracted surfaces (exact, BigInt):** every Cubes vertex is a zero of the
  trilinear interpolant of `2v - 2L - 1` (edge crossings and plus-sign branch
  points alike); every Tetrahedra triangle lies in one tetrahedron of the
  parity-alternating decomposition, on the zero set of its linear
  interpolant (where swapped zero-corner case bodies would show, cf. v16
  #546; none did); every saddle face of a Cubes mesh (4285 in the deep run)
  contains the two segments that cut off the corners on the other side of
  the bilinear saddle value (this check fails on upstream's output, suspect
  1); meshes are closed (every edge of two triangles except on the image's
  boundary planes) for Tetrahedra without samples on the level and for Cubes
  without four-crossing faces. MC: every vertex lies on its cube edge at the
  crossing of the linear interpolant (residual at rounding level), and every
  deduplicated mesh without a voxel on the level (3226 deep records, 1586 of
  them with ambiguous faces) is closed and consistently oriented (no
  directed edge twice).
* **MarchingCubes table:** on each of the 192 ambiguous face instances of
  the 256 entries, the triangles contain exactly two segments in the face
  and both cut off the negative corners (`F < level`); the resolution depends
  on the face alone, so the table is face-consistent (suspect 4).

## Port defects fixed

| file | what differed | root cause | size |
| --- | --- | --- | --- |
| `src/ImageUtility3.ts` `dilate` | the port threw on an empty offset list; upstream returns a copy of the input | the port added the `numNeighbors > 0` assertion of the 2D functions and the 3D `Erode`, which upstream's 3D `Dilate` does not have | throw parity (1 of 20 golden records) |
| `src/SurfaceExtractorMC.ts` `orientTriangles` | a triangle oriented the other way on 3 of 2000 deep `.extract.large` records | upstream divides the `Vector3` gradient sum by 3 with `Vector.h`'s `operator/`, which multiplies by `1/3`; the port divided each component by 3. They differ in the last bit, and where the average gradient is orthogonal to the normal in exact arithmetic (integer gradients) the rounding decides the sign of the dot product | a discrete orientation flip; `.orientTriangles.thirds` is aimed at it (1613 of 2000 deep records flip with the old formula) |
| `src/SurfaceExtractorMC.ts` `extractVoxel` | `-0` against `+0` in vertex coordinates on 44 of 2000 deep `.extractVoxel` records (level -0) | the #443 fix evaluates `(F[j0] - level) / (F[j0] - F[j1])`, and `-0 - (-0)` is `+0`; at level 0 upstream is correct, so the port now evaluates upstream's `F[j0] / (F[j0] - F[j1])` there and the shifted numerator only for a nonzero level | sign of zero; `.extractVoxel.signedZero` is aimed at it |

Two further changes fix upstream defects found by this group (suspects 1
and 2). Regression tests: `test/ImageUtility3.test.ts`,
`test/SurfaceExtractorCubes.test.ts` (the two existing face-ambiguity tests
encoded upstream's inverted pairing and were corrected; a new test pins the
header example), `test/SurfaceExtractorMC.test.ts` (the reciprocal, the zero
sign, `ComputeNormals`; each fails on the previous port). The existing test
"an ambiguous shared face leaves a hole" was corrected (suspect 4).

**Sensitivity.** Replaying the deep run with the base class's average
gradient computed as `sum * (1/3)` instead of upstream's scalar `/ 3`
changes 2000 of 2000 `SurfaceExtractorCubes.orientTriangles.thirds` records
and 1, 1 and 11 records of Cubes `.orientTriangles.points`, Tetrahedra
`.extract` and `.extract.large`; the old MC formula changes 1613 of 2000
`.orientTriangles.thirds` records (the others exhausted their 8192 redraws
without a qualifying draw).

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `ImageUtility2.drawEllipse.zeroExtents` | [#443](https://github.com/gradientspaceai/gtengine-js/issues/443) | 2000 of 2000: upstream never terminates (the C++ callback throws after 1000 calls, a throw record); the port visits the centre once. |
| `ImageUtility3.morphology.xmin` | [#129](https://github.com/gradientspaceai/gtengine-js/issues/129) | 317 of 2000 through `Dilate`, `Open` and `Close` on images with 1s on the x = 0 face. In the residue (1628 records and 55 throw records) no 1-voxel on the x = 0 face reaches a voxel that another source does not also reach (the generator's random modes do not aim at the face, and `Open` erodes first); the direct-definition check passes on the port for every record, and upstream's result differs from the definition on exactly the deviating ones. |
| `SurfaceExtractorCubes.extract.saddle` | this report, suspect 1 | 2000 of 2000; each record has a saddle face with a nonzero determinant (at most 64 redraws, then the header example). The saddle check passes on the port's output for every record. |
| `SurfaceExtractorMC.extractVoxel.level`, `.extract.level` | [#443](https://github.com/gradientspaceai/gtengine-js/issues/443) item 1 | 1352 and 1172 of 2000; the rest have no vertex (all corners on one side, or a corner on the level with `perturb` 0). `.levelTopology` compares the tables, counts and indices at nonzero levels and agrees everywhere. |
| `SurfaceExtractorMC.computeNormals.multiple` | this report, suspect 2 | 1432 of 2000. In 567 of the 568 agreeing records the triangles after the first, apart from copies of the first, add zero at every vertex (degenerate triangles, frequent with repeated indices among 3 to 7 points, or a triangle cancelled by its reversal), so both sides normalise multiples of the first normal; one more agrees by coincidence of the normalised sums. The main cases compare `ComputeNormals` on single-triangle meshes. |
| `SurfaceExtractorTetrahedra.orientTriangles.centralTetra` | [#132](https://github.com/gradientspaceai/gtengine-js/issues/132) | 2000 of 2000; every vertex lies in the central tetrahedron of an odd-parity cube and each record is kept only when some triangle's orientation depends on that gradient. Extracted vertices always lie on a cube face and never reach the region, so `.extract` agrees; `.orientTriangles.points` rejects the region. |

Preserved findings reached and agreeing bit for bit: #443 (`GetSkeleton`
on even squares, statistics above), #129 (`DrawLine`'s unused `maxValue`
assignment; the 3D `Close<N>` does not compile, so the port's
`close6/18/26` are compared against upstream's offset-list `Close` with the
same neighbourhoods). #439 (`SurfaceExtractor::MakeUnique` keeps rotated
duplicate triangles) belongs to group 25's base class; `MakeUnique` is
compared bit for bit here through both extractors.

## Not covered

| header / entry point | reason |
| --- | --- |
| `ImageUtility2/3` `GetComponents`, `ExtractBoundary`, `GetSkeleton` on images violating the zero-border precondition | the neighbour reads have no range test (#129, preserved): upstream reads outside the pixel array, undefined behaviour. |
| `GetComponents`, `ExtractBoundary` output accumulation | upstream appends to (or leaves unchanged) the caller's vectors; the port returns fresh arrays. The cases pass fresh vectors. |
| `ImageUtility3::Close<N>` | does not compile upstream (#129); compared through the offset-list `Close`. |
| `FloodFill4/6`, `DrawFloodFill4` with equal foreground and background colours | never terminate (the stack overruns its allocation); excluded by the generators. |
| `DrawEllipse` extents above about 800, `DrawCircle`/`DrawLine` coordinates near the `int32_t` range | signed overflow in the decision terms (`4 a^2 b` exceeds 2^31), undefined behaviour. |
| `GetL2Distance` with distances of 100 or more | outside the documented exactness range. |
| `SurfaceExtractorCubes/Tetrahedra` with 32-bit voxel values beyond 2^20 | upstream's `int64_t` products can overflow and the port's doubles lose exactness above 2^53 (documented port limit). |
| `SurfaceExtractorMC` with `T = float` or another `IndexType` | the port has one number type. |
| `SurfaceExtractorMC::MakeUnique` with NaN vertices | `std::map` needs a strict weak ordering; NaN breaks it. |
| `SurfaceExtractor.h` itself (`Vertex`, `Triangle`, the base's own cases) | group 25; only reached here through the Cubes and Tetrahedra extractors. |

## Upstream bug suspects

**1. `SurfaceExtractorCubes.h` `Get{X,Y,Z}{Min,Max}Edges`, case 15: the
saddle-face pairing is inverted** (result-corrupting; fixed in the port).
When all four edges of a voxel face are crossed and the face determinant
`det = f00*f11 - f01*f10` of the shifted values `2v - 2L - 1` is nonzero,
the `det > 0` branch joins the two crossings next to corner `f00` and the
two next to `f11`, and `det < 0` the ones next to `f10` and `f01`. For the
bilinear interpolant of the face the saddle value has the sign of
`det * f00` (the denominator `f00 + f11 - f01 - f10` has the sign of
`f00`), so for `det > 0` the corners `f00` and `f11` are connected through
the saddle and it is `f10` and `f01` that must be cut off: upstream is wrong
on every such face, in all six functions. Unlike v24's
`AdaptiveSkeletonClimbing2` (#544) the determinant does include the level;
only the two branch bodies are swapped. Reproduction: the 2x2x2 image with
3 at (0,0,0) and (0,1,1) and 0 elsewhere, level 0: face x = 0 has shifted
values 5, -1, 5, -1 and centre value 2 > 0, but upstream returns two caps
around the positive corners, whose chord through (0, 0.4167, 0.4167) passes
where the (shifted) interpolant is 2.08 instead of 0. Adjacent
voxels make the same choice on a shared face, so the mesh stays closed and
only its topology is wrong. The port swaps the branches; the main cases
reject images with such a face (exact integer predicate on the C++ side) and
agree bit for bit.

**2. `SurfaceExtractorMC.h` `ComputeNormals`: the triangle pointer is never
advanced** (result-corrupting). `IndexType const* triangle = indices.data();`
is read in the loop over the triangles without `triangle += 3` (which
`OrientTriangles` has), so the first triangle's normal is added
`numTriangles` times at its three vertices and every other vertex gets the
zero normal. Reproduction: two disjoint triangles (0,0,0), (1,0,0), (0,1,0)
and (0,0,2), (0,1,2), (0,0,3): upstream returns (0,0,1) for the first three
vertices and (0,0,0) for the last three instead of (1,0,0). The port already
iterated over every triangle (an undocumented deviation until now); it is
now documented in `src/SurfaceExtractorMC.ts`, pinned by a regression test
and demonstrated by `.computeNormals.multiple`.

**3. `SurfaceExtractorCubes.h` `VETable::RemoveTriangles`: geometry-free ear
clipping makes non-manifold meshes on voxels with a four-crossing face**
(quality, preserved). `Remove` always takes the lowest-numbered vertex of
degree 2 as the next ear, so on a voxel whose wireframe has four vertices on
one face (a saddle or plus-sign face) a fan diagonal can lie inside that
face, and the two voxels sharing the face choose their diagonals
independently: edges of three or four triangles and single-triangle edges
inside the image. In the deep run 443 of the 4399 Cubes records with a
four-crossing face are affected; the upstream pairing has the same defect
(in the saddle case 43 records with odd edge counts and 71 with four-fold
edges, against 21 and 36 after the fix of suspect 1). Meshes without a
four-crossing face are always closed.

**4. Correction to [#443](https://github.com/gradientspaceai/gtengine-js/issues/443)
item 2 ("the 15-case table is not face-consistent ... leaving a hole").**
The `MarchingCubes` table resolves each of its 192 ambiguous face instances
by cutting off the negative corners, so two voxels sharing an ambiguous face
always agree. The pinned 3x2x2 test's 8 single-triangle edges all lie on the
image's boundary planes y = 0, y = 1, z = 0, z = 1; the face x = 1 is closed
(the left voxel caps the negative corners, the right one joins the positive
corners by a tube). 3226 deep meshes without a voxel on the level, 1586 of
them with ambiguous faces, are closed and consistently oriented. The fixed
rule is a disambiguation (it can disagree with the trilinear interpolant's
saddle), not a hole.

**5. Minor.** `ImageUtility3`'s offset-list `Dilate` lacks the
`numNeighbors > 0` assertion that the 2D functions and the 3D `Erode` and
`GetComponents` have. `ImageUtility2::GetL1Distance` reports a maximum
distance of 1 at (0,0) for an image without foreground (the loop's first
pass decides `maxDistance = 1`).
