# Group 16 (core) — `v16-core`

34 cases, 20 records each in `oracle/golden/v16-core.txt`: 31 declared
`exact` (5 of them `deviation`s) and 3 with a tolerance (libm); integer,
boolean and bit-pattern outputs always compare exactly. Headers:
`CurveExtractorSquares.h`, `CurveExtractorTriangles.h`, `IEEEBinary16.h`,
`MeshSmoother.h`, `PolygonTree.h`, `PolygonWindingOrder.h`,
`ETNonmanifoldMesh.h`, `ImplicitSurface3.h`, `MeshCurvature.h`,
`RevolutionMesh.h`, `TubeMesh.h`, `VETNonmanifoldMesh.h`,
`VertexCollapseMesh.h`.

Deep run: `npm run oracle:deep -- 2000 v16-core`: 68 000 records (7641 of
them recorded C++ throws), **all 35 tests pass** (34 cases and the claim
check), about 85 s including generation (most of it the #498 search of the
`VertexCollapseMesh` deviation case). Outside the five deviation cases the
deep run compares 6 508 777 real outputs; **6 448 084 are bit-identical**,
and all 60 693 others are in the two libm cases `RevolutionMesh.libm` (max
scaled error 2.1e-13) and `TubeMesh.libm` (2.2e-14). `IEEEBinary16.mathLibm`
is declared with a tolerance but was bit-identical on all 48 000 outputs.
The goldens are deterministic: the committed goldens and a 2000-record deep
file were each generated twice, byte-identical (including the
`VertexCollapseMesh.doCollapse` case, which runs the unmodified upstream
class over its pointer-keyed unordered sets).

**Four port defects were fixed** (three in `IEEEBinary16`, one new upstream
defect in `CurveExtractorSquares` fixed per PORTING.md) and five deliberate
fixes are demonstrated as deviations.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `CurveExtractorSquares.h` | `CurveExtractorSquares.extract` (rational `Extract`, `MakeUnique`, real `Extract` with and without duplicate removal) | exact | pass, 109 910 reals bit-identical |
| | `CurveExtractorSquares.extract.threeZeroCorners` | deviation (new finding, below) | 2000 of 2000 records deviate |
| `CurveExtractorTriangles.h` | `CurveExtractorTriangles.extract` (same four outputs) | exact | pass, 135 110 reals |
| `IEEEBinary16.h` | `.convert32To16`, `.convert16To32`, `.fields`, `.fromNumber`, `.compare`, `.arithmetic`, `.mathExact` | exact | pass (after the fixes below) |
| | `.mathLibm` (24 std:: and Functions.h wrappers) | tolerance 1e-3 per output: binary32 libm (MSVC `acosf`, ...) against V8's binary64 libm rounded to binary32; one binary32 ulp moves a binary16 result by at most one binary16 ulp (2^-10 relative) | pass, 48 000 of 48 000 bit-identical |
| `PolygonWindingOrder.h` | `PolygonWindingOrder.operator` | exact | pass |
| `PolygonTree.h` | `PolygonTree.getContainingTriangle` (all three `PolygonTreeEx` queries), `.invalidArgument` (throw parity, 1009 throws) | exact | pass |
| `MeshSmoother.h` | `MeshSmoother.update` (three `Update(t)` steps, default and overridden hooks), `.invalidInput` (throw parity, 1618 throws) | exact | pass, 513 810 reals |
| `MeshCurvature.h` | `MeshCurvature.compute` | exact | pass, 203 599 reals |
| `ImplicitSurface3.h` | `ImplicitSurface3.queries` (`IsOnSurface`, `F`, `GetGradient`, `GetHessian`, `GetFrame`, `GetPrincipalInformation`) | exact | pass, 60 000 reals |
| `ETNonmanifoldMesh.h` | `ETNonmanifoldMesh.sequence`, `.remove.degenerate` (throw parity, #179 preserved, 2000 throws) | exact | pass |
| `VETNonmanifoldMesh.h` | `VETNonmanifoldMesh.sequence` | exact | pass |
| | `VETNonmanifoldMesh.remove.isolatesVertex` | deviation (#240) | 2000 of 2000 |
| `RevolutionMesh.h` (+ `Mesh.h` base) | `RevolutionMesh.construct` (CYLINDER, TORUS, DISK), `.sphere` (positions, tcoords), `.invalid` (unsupported topology; missing index/position channel throws, 1014 throws) | exact | pass, 823 579 reals |
| | `RevolutionMesh.libm` | tolerance 1e-12 (std::cos/std::sin of the column angles) | pass, 97.92 % bit-identical, max 2.1e-13 |
| | `RevolutionMesh.sphere.indices` | deviation (#220, #240: `Mesh.h` pole fans) | 2000 of 2000 |
| `TubeMesh.h` | `TubeMesh.construct` (open tubes) | exact | pass, 715 466 reals |
| | `TubeMesh.libm` | tolerance 1e-12 (std::cos/std::sin) | pass, 98.80 % bit-identical, max 2.2e-14 |
| | `TubeMesh.closed` | deviation (#240: closed-tube fixup stride) | 2000 of 2000 |
| `VertexCollapseMesh.h` | `.doCollapse` (unmodified class, discrete outputs), `.doCollapse.canonical` (every heap weight bit for bit), `.invalidInput` | exact | pass, 225 537 reals |
| | `.doCollapse.invalidLink` | deviation (#498) | 2000 of 2000 |

### How the cases reach the branches

* **Curve extractors.** Images of 2..6 x 2..6 pixels in five modes: values
  and level in [-2, 2] (zero corners: every sign pattern), `uint8_t`
  images, the full `int16_t` range, a `+-+-` checkerboard of +-(1..6)
  aimed at the saddle determinant, and `int32_t` values up to 2^20. Deep
  run: 162 / 1743 / 1804 saddle squares with `det == 0` / `> 0` / `< 0`
  (squares), 173 / 1702 / 1651 (triangles' images). The preserved
  `MakeUnique` finding (#67, #362: reversed duplicate edges survive) is
  reached on 71 (squares) and 117 (triangles) of 2000 deep records, and
  agrees. The replay checks, with bigint rationals, that every extracted
  vertex lies on the level set of the interpolant the extractor contours
  (bilinear per square; linear per triangle with the parity-alternating
  diagonal): that check found the new upstream defect below.
* **IEEEBinary16.** Encodings travel as integers. `convert32To16` draws
  from ten modes: zero and binary32 subnormals, the five class thresholds
  (2^-25, 2^-24, 2^-14, max normal, the max-normal/infinity midpoint) +-2
  encodings, 16-normal and 16-subnormal values whose dropped bits are
  exactly half, half +-1 or arbitrary (ties to even both ways), signaling
  NaNs with payload below 2^13 (finding #110: infinity, preserved) and
  above, quiet NaNs, infinities, overflow and underflow ranges. The replay
  checks every result against an independent round-to-nearest-even model
  and every half against its 16 -> 32 -> 16 round trip. `fromNumber` covers
  the double -> float -> half double rounding (doubles on and 1 ulp beside
  binary32 ties), NaN payloads read with `io.real()`. `fields` covers
  every accessor, classifier and `GetNextUp`/`GetNextDown` of the base
  template and `SetEncoding` with fields out of range (`uint16_t`
  truncation). `arithmetic` covers unary minus, the twelve binary
  operators with half and float operands and the eight compound updates.
* **Polygons.** `PolygonWindingOrder`: star-shaped lattice and real
  polygons in both orientations, collinear runs through the lower-left
  vertex, fully collinear polygons (`dotPerp == 0` returns false). The
  replay checks, on simple lattice polygons whose lower-left turn is
  strict, that the answer is the sign of the exact shoelace area.
  `PolygonTreeEx`: breadth-first trees of 1..6 nodes with chirality -1, 0,
  +1 and degenerate lattice triangles, test points at vertices, on edge
  midpoints and elsewhere, so every `sdot > 0` is evaluated at exact zero.
* **Surface meshes** (`MeshSmoother`, `MeshCurvature`): tetrahedron,
  octahedron, cube (integer-scaled lattice or perturbed), heightfield grids
  with lattice or real heights and random diagonals, and a once-subdivided
  octahedron projected to a sphere of radius R. `MeshSmoother` runs three
  `Update(t)` steps with the default weights or a subclass overriding all
  three hooks (defined identically in the replay), and one record in four
  has a vertex no triangle references (Vector's `operator/=` zeroes its
  mean). Sensitivity: upstream's `mMeans[i] /= n` is a multiplication by
  `1/n`; a per-component `x / n` would differ from it on 964 of the 2000
  deep records in the first step alone. `MeshCurvature` uses the thresholds
  0 (the documented default, #240: never fires), 1e-3 and 1e10 (every vertex
  planar); the unperturbed regular octahedron is the umbilic configuration
  of #412 (both directions zero, preserved). The replay checks that the
  unperturbed sphere samples estimate both curvatures within [0.5, 1.5]/R.
* **ImplicitSurface3** through a cubic polynomial subclass with every sum
  grouped explicitly, mirrored in the replay (the v18 precedent): random
  coefficients, the sphere of radius R at Pythagorean lattice points or a
  scaled unit vector (the replay checks both principal curvatures are 1/R
  to 1e-12), lattice coefficients and positions (repeated eigenvalues of
  `SymmetricEigensolver2x2`), and zero gradients (the `false` branch).
* **ET/VETNonmanifoldMesh.** Up to 16 (20) inserts and removes over 4..6
  vertices: shared and nonmanifold edges, the same key in another rotation
  (`Insert` returns null), reversed triangles and degenerate triangles
  (#179). Every container these classes iterate is ordered by value
  (`std::map` of `EdgeKey`/`TriangleKey`, `std::set` with
  `WeakPtrLT`/`SharedPtrLT`, which compare the pointees' keys), so every
  output is emitted in upstream's own order: the edge map with each edge's
  stored vertex order and triangle set, the triangle map with each
  triangle's edges, `IsManifold`, `IsClosed`, both `GetComponents`
  overloads (depth-first order), the vertex map with `VAdjacent`,
  `EAdjacent`, `TAdjacent`, then a copy (operator=, which reinserts by key,
  rotating vertices) and `Clear`.
* **RevolutionMesh / TubeMesh** use `BezierCurve` profiles and medial
  curves (arithmetic only, bit-identical in v17), a quadratic radial
  function mirrored in the replay, zero (Frenet) and lattice up vectors,
  uniform and arc-length sampling, both windings and every channel
  combination (normals, client texture coordinates, the four tangent-space
  channels with and without `wantDynamicTangentSpaceUpdate`, i.e.
  `UpdateFrame` versus `UpdateNormals`). The only libm calls are the
  column-angle tables. A stand-alone probe compared MSVC's and V8's `cos` and
  `sin` of `c * (1/n) * 2pi` and of `2pi * c / n` (DISK texture coordinates)
  for n = 1..80: the exact cases draw the column count from the n on which
  every table entry agrees bit for bit (revolution n in {3,4,5,6,7,9,10,14};
  tube numCols - 1 in {2,3,4,5,6,7,9,10,14,18}), so the whole mesh, frames
  included, is compared exactly; the `.libm` cases take any column count.
  The case file includes `Matrix2x2.h`, so `Mesh::UpdateFrame`'s dependent
  `Inverse(mUTU[i])` resolves to the closed-form 2x2 inverse the port
  implements (the v34/v43 overload lesson).
* **VertexCollapseMesh.** `VETManifoldMesh`'s vertex map and adjacency sets
  are unordered (`unordered_map<int32_t>`, `unordered_set<Triangle*>`,
  `unordered_set<int32_t>`); their order decides the heap insertion order
  (weight ties), the accumulation order of `ComputeWeight` (the last bits
  of every weight and vertex normal) and the order of `Record::removed`.
  The port iterates sorted (vertices by index, triangles by key, adjacent
  vertices ascending). `.doCollapse` runs the unmodified class on meshes
  with real coordinates and compares the discrete outputs (removed keys
  sorted): 2000 of 2000 agree. `.doCollapse.canonical` opens the class's
  private members and replays its constructor heap and `DoCollapse` in the
  port's order around upstream's own `TriangulateLink` and `Collapsed` (the
  v12 `Canonical<Base>` precedent), on real and lattice heightfields
  (exact weight ties, collinear links, deferred vertices) and closed
  solids, and emits the weight of every vertex in the heap after every
  step: 225 537 weights bit-identical. It is restricted to meshes on which
  no step meets the #498 defect (probed with the same driver).

## Port defects fixed

| file | what differed | root cause | size |
| --- | --- | --- | --- |
| `src/IEEEBinary16.ts` `sinpi`, `cospi`, `exp10`, `atandivpi`, `atan2divpi`, `invsqrt` | `sinpi(65504)` was -0 in the port, -0.00544 in C++ | the `IEEEBinary16` wrappers of `Functions.h` call the **float** overloads upstream (`std::sin(x * static_cast<float>(GTE_C_PI))` with a binary32 product; `1.0f / std::sqrt(x)`), the port called the double overloads | up to 5.7e-3 absolute for `sinpi`/`cospi` near \|x\| = 65504; 102 of 2000 deep `mathLibm` records disagreed before the fixes (together with `pow` below) |
| `src/IEEEBinary16.ts` `ldexp` | `ldexp(inf, -1100)`, `ldexp(0, 1100)` were NaN (C++: inf, 0) | `x * Math.pow(2, e)` forms `inf * 0` once `2^e` leaves binary64; `std::ldexp` returns zero, infinity and NaN unchanged | NaN instead of the value; 13 of 2000 deep `mathExact` records |
| `src/IEEEBinary16.ts` `pow` | `pow(1, NaN)` and `pow(-1, +-inf)` were NaN (C++: 1) | C's `pow` (C99 F.10.4.4) and `Math.pow` differ on exactly these cases | NaN instead of 1; 5 of 2000 deep `mathLibm` records |
| `src/CurveExtractorSquares.ts` `'+000'`, `'00+0'` | see "Upstream bug suspects" 1: an upstream defect found by this group's level-set check, fixed per PORTING.md | the two case bodies are swapped upstream | two wrong edges per affected square |

The `ldexp` fix clamps the exponent to [-300, 300] (every binary16 times
2^+-300 already overflows or underflows binary32, so the product stays
exact and `Math.fround` rounds it correctly); the `pow` fix returns 1 for
`x == 1` and for `x == -1, y = +-inf`. Regression tests are in
`test/IEEEBinary16.test.ts` and `test/CurveExtractorSquares.test.ts`.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `CurveExtractorSquares.extract.threeZeroCorners` | new finding of this group (suspect 1) | 2000 of 2000 records: a mode-0 image with one `'+000'` or `'00+0'` square planted; upstream emits the two edges incident to the nonzero corner, the port the two zero edges. The main case rejects images containing such a square (at most 64 redraws, then an `int16_t` image) and agrees everywhere; the port passes the level-set check on every record of both cases. |
| `VETNonmanifoldMesh.remove.isolatesVertex` | [#240](https://github.com/gradientspaceai/gtengine-js/issues/240) (UPSTREAM-FINDINGS, mesh item 6) | 2000 of 2000: upstream's inverted `LogAssert` throws on every valid removal that leaves a vertex of the triangle without triangles; the port removes it. The main case turns such removals into inserts and agrees everywhere. |
| `RevolutionMesh.sphere.indices` | [#220](https://github.com/gradientspaceai/gtengine-js/issues/220), [#240](https://github.com/gradientspaceai/gtengine-js/issues/240) (`Mesh.h` items 1 and 2) | 2000 of 2000: the second pole fan's start vertex (row stride `numCols` instead of `numCols + 1`) and its winding; the normals and frames computed from the indices follow. `RevolutionMesh.sphere` compares the positions and texture coordinates of every SPHERE mesh exactly (including the NaN texture coordinates of one-row spheres, #412, preserved). |
| `TubeMesh.closed` | [#240](https://github.com/gradientspaceai/gtengine-js/issues/240) (`TubeMesh` item 3) | 2000 of 2000: the closed-tube fixup's stride and its missing seam vertex; counts, flags, texture coordinates and indices agree. |
| `VertexCollapseMesh.doCollapse.invalidLink` | [#498](https://github.com/gradientspaceai/gtengine-js/issues/498) | 2000 of 2000: every record has a step on which upstream's `TriangulateEC` returns an invalid triangulation of a collinear link **and** `Collapsed` goes on to modify the mesh (an invalid triangulation whose repeated edges already carry two triangles is deferred by upstream's own diagonal test, exactly like the port, and agrees; selecting on the invalid triangulation alone left only 124 of 2000 deviating). Every fourth record is the #498 mesh; the others are flat grids with two raised vertices or lattice grids found by the probe (151 distinct meshes in the deep run; draws that find none within 64 attempts fall back to the #498 mesh). |

Preserved findings reached and agreeing: #110 (`Convert32To16`, NaN payload
below 2^13 becomes infinity), #179 (degenerate-triangle aliasing: `Insert`
accepts it, `Remove` throws on both sides), #67/#362 (reversed duplicate
edges after `MakeUnique`), #412 (DISK texture coordinates without the
origin; one-row SPHERE NaN texture coordinates; umbilic principal
directions), #240 (`TubeMesh` rings with `numCols - 1` distinct angles;
`MeshCurvature`'s zero threshold), #295 (`record.vertex = 0x80000000`,
emitted as INT32_MIN on both sides).

## Not covered

| header / entry point | reason |
| --- | --- |
| `IEEEBinary16.h`: sign and payload of a NaN **produced** by arithmetic or a libm wrapper | x64 SSE returns the default NaN `0xFFC00000` or a quieted operand; a NaN computed in JavaScript has no reproducible bits. Such results are compared as NaN (the harness matches any NaN with any NaN), and a NaN half result is emitted as its NaN flag. NaN **inputs** are covered bit for bit through `convert32To16`, `convert16To32`, `fields` (integers) and `fromNumber` (`io.real()`). |
| `IEEEBinary16.h`: the float constructor on binary32 signaling NaNs | the port has no binary32 argument type (a JS number holding a binary32 sNaN cannot exist: widening quiets it); the conversion itself is compared through the exposed `convert32To16`. |
| `IEEEBinary16.h`: `std::frexp` exponent of infinity and NaN | unspecified by the C standard; the exponent is compared for finite arguments only. |
| `CurveExtractor*.h`: pixel magnitudes where the products or the saddle determinant exceed 2^53 | the port's documented limitation (number instead of `int64_t`; `int32_t` images near their range also overflow `int64_t` upstream in `det`). Images over the full `int16_t` range and `int32_t` values up to 2^20 are covered. |
| `ETNonmanifoldMesh.h`, `VETNonmanifoldMesh.h`: custom `ECreator`/`TCreator`/`VCreator` | only the default creators are exercised; the callbacks are plumbing. |
| `Mesh.h` `UpdateFrame` without texture coordinates | unreachable from `RevolutionMesh` and `TubeMesh`, which always allocate default texture coordinates; `Mesh.h` belongs to group 46. |
| `RevolutionMesh.libm` on TORUS meshes with 2 rows or a quadratic profile | ill-conditioned inputs, not a port difference: with two rows (t = tmin and tmax coincide) or a closed quadratic profile (it retraces itself, which the documented precondition excludes) the TORUS surface folds onto itself and `UpdateFrame`'s per-vertex sums cancel, so a 1-ulp `cos`/`sin` difference flips the frame. Measured on the deep run: 27 of 2000 records disagreed (all in the frame and normal channels of such meshes, positions within 1e-12); 4 remained with the row restriction alone, 24 with the degree restriction alone, none with both. The exact case covers two-row tori bit for bit; all TORUS profiles are now cubic or quartic. |
| `VertexCollapseMesh.h` with lattice meshes through the **unmodified** class | exact weight ties are then decided by `unordered_map`/`unordered_set` iteration order, which the port does not reproduce (suspect 2); the canonical-order case compares them. |
| `PolygonTree` (the class) | pure data (a polygon and its children); `TriangulateEC` and `TriangulateCDT` (groups 8 and 13) exercise it. |

## Upstream bug suspects

**1. `CurveExtractorSquares.h`: the `'+000'` and `'00+0'` cases emit the
wrong pair of square edges** (result-corrupting; found by this group's
level-set reference). In upstream's `(f00, f10, f11, f01)` sign notation,
`'+000'` is a square whose only nonzero corner is `(x, y)`. The zero set of
the bilinear interpolant `f00 (1-u)(1-v)` is the right edge `u = 1` and the
top edge `v = 1`, but the code emits `(x,yp)-(x,y)` and `(x,y)-(xp,y)`, the
left and bottom edges, which meet at the nonzero corner. `'00+0'` (only
`(xp, yp)` nonzero, zero set: left and bottom edges) emits the right and top
edges. The two bodies are swapped: each is the other's correct answer. The
other two three-zero cases, `'0+00'` and `'000+'`, are correct. Example: the
2x2 image `{f00, f10, f01, f11} = {1, 0, 0, 0}` at level 0 yields the edges
`(0,1)-(0,0)` and `(0,0)-(1,0)`, both through `(0,0)`, where the image is 1.
Port: fixed (bodies swapped), regression tests in
`test/CurveExtractorSquares.test.ts`, deviation case
`CurveExtractorSquares.extract.threeZeroCorners`.

**2. `VertexCollapseMesh.h`: the decimation depends on unordered-container
iteration order** (determinism, minor). `VCVertex::ComputeWeight` sums over
`TAdjacent` (`unordered_set<Triangle*>`, pointer-hashed) and `VAdjacent`
(`unordered_set<int32_t>`), so the weights and vertex normals depend in the
last bits on heap addresses and bucket layout; the constructor inserts the
vertices into the min-heap in `unordered_map` order, which decides exact
weight ties (common on symmetric or lattice meshes); and `Record::removed`
lists the triangles in `TAdjacent` order. Two runs of the same build agreed
here, but another standard library, allocator or address layout can
collapse a different vertex. The port iterates in sorted order.
