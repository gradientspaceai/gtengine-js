# Verify group 29 (interpolation) against the MSVC build of upstream GTE

Family `v29-interpolation`, 36 cases, 20 golden records each
(`oracle/golden/v29-interpolation.txt`, 822 KB). 29 cases are declared
`{ exact: true }`, 2 carry a tolerance because upstream calls the C math
library, 5 are `deviation` cases.

Deep run: `npm run oracle:deep -- 2000 v29-interpolation` - 72000 records
(15482 of them C++ throws), **all 37 tests pass** (36 cases plus the
"every golden case has a replay" check), replay 53 s. Over the deep run every
real output of every `exact` case was bit-identical (1910276 real outputs);
the two tolerance cases had 4655 and 152 inexact outputs, with the maxima
listed below. Goldens are reproducible: three generations of 2000 records per
case were byte-identical.

The group's headers are `IntpTricubic3.h`, `IntpTrilinear3.h`,
`IntpAkimaNonuniform1.h`, `IntpAkimaUniform1.h` (and through them the public
surface of `IntpAkima1.h`), `IntpBSplineUniform.h`,
`IntpLinearNonuniform2.h`, `IntpLinearNonuniform3.h`,
`IntpThinPlateSpline2.h`, `IntpThinPlateSpline3.h`,
`IntpQuadraticNonuniform2.h`, `IntpSphere2.h` and `IntpVectorField2.h`.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `IntpTrilinear3.h` | `IntpTrilinear3.evaluate` (all 13 accessors; `operator()(x,y,z)`; `operator()(xo,yo,zo,x,y,z)` for every order triple in {-1,0,1,2}^3, i.e. both evaluated branches and the `default: return 0` of each axis; 3 queries per record) | exact | pass, 408000 reals |
| | `IntpTrilinear3.construct.invalid` (bounds 0..3, spacing positive / 0 / -0 / negative / NaN) | exact, throw parity | pass, 1938 throws |
| `IntpTricubic3.h` | `IntpTricubic3.evaluate` (Catmull-Rom and B-spline blending; accessors; value; all 64 order triples in {0..3}^3 plus orders 4 and -1) | exact | pass, 290000 reals |
| | `IntpTricubic3.construct.invalid` | exact, throw parity | pass, 1947 throws |
| `IntpTrilinear3.h`, `IntpTricubic3.h` | `IntpTrilinear3.evaluate.signedZeroSamples` (every sample +0 or -0, so the zero seeds of `P`, `Q`, `R` and `result` decide the sign; each derivative emitted from an inlined and from a `__declspec(noinline)` call, per the v23 seed-folding lesson; both agree, all outputs +0) | exact | pass, 36000 reals |
| `IntpAkimaUniform1.h`, `IntpAkima1.h` | `IntpAkimaUniform1.evaluate` (`GetQuantity`, `GetXMin`, `GetXMax`, `GetXSpacing`; every polynomial coefficient through a subclass that exposes the protected `mPoly`; `operator()(x)` and `operator()(order,x)` for orders -1..4 at 5 queries: nodes, dyadic in-cell points, both ends, -0, outside (clamped)) | exact | pass, 112392 reals |
| | `IntpAkimaUniform1.construct.invalid` (quantity 0..4, spacing 0 / -0 / negative) | exact, throw parity | pass, 1678 throws |
| `IntpAkimaNonuniform1.h`, `IntpAkima1.h` | `IntpAkimaNonuniform1.evaluate` (the same outputs; nodes on integer, dyadic and uniform increments, 3..8 samples) | exact | pass, 109748 reals |
| | `IntpAkimaNonuniform1.construct.invalid` (too few samples, repeated node, decreasing node) | exact, throw parity | pass, 1500 throws |
| `IntpBSplineUniform.h` | `IntpBSplineUniform1.evaluate`, `IntpBSplineUniform2.evaluate`, `IntpBSplineUniform3.evaluate` (the three specializations; degrees 1..5 / 1..3 / 1..2; all three cache modes; `GetDegree`, `GetNumControls`, `GetTMin`, `GetTMax`, `GetCacheMode`; every order tuple from -1 (or 0) to degree+1; several evaluations per object so ON_DEMAND_CACHING reads what earlier calls cached; `ctZero` = +0 or -0) | exact | pass, 558579 reals |
| | `IntpBSplineUniform.evaluate.general` (the run-time-dimension class `IntpBSplineUniform<Real, Controls, 0>`, i.e. `IntpBSplineUniformShared::EvaluateNoCaching` / `EvaluateCaching`, `ComputeTensor`, `InitializeTensors`, `GetRowIndex`, `GetIndex`, for N = 1..4; every order tuple; the size guard for short `order` / `t`; a negative order) | exact | pass, 245281 reals |
| | `IntpBSplineUniform.evaluate.fixedN` (the compile-time `IntpBSplineUniform<Real, Controls, 4>`) | exact | pass |
| | `IntpBSplineUniform.construct.invalid` (numControls <= degree + 1 in all four classes) | exact, throw parity | pass, 1343 throws |
| `IntpLinearNonuniform2.h` | `IntpLinearNonuniform2.evaluate` (the real `Delaunay2Mesh<double>`, queries strictly inside a triangle or outside the hull) | exact | pass |
| | `IntpLinearNonuniform2.evaluate.sortedMesh` (every query class: interior, edge, vertex, outside) | exact | pass |
| | `IntpLinearNonuniform2.deviation.getIndices` | deviation (#135) | 1955 of 2000 deviate |
| `IntpLinearNonuniform3.h` | `IntpLinearNonuniform3.evaluate`, `.evaluate.sortedMesh` (interior, face, edge, vertex, outside) | exact | pass |
| | `IntpLinearNonuniform3.deviation.getIndices` | deviation (#135) | 1966 of 2000 deviate |
| `IntpQuadraticNonuniform2.h` | `IntpQuadraticNonuniform2.fromDerivatives` (real `Delaunay2Mesh<double>`, strict queries), `.fromDerivatives.sortedMesh`, `.fromSpatialDelta.sortedMesh` (spatialDelta 1, 0.5, 0, uniform); F, FX and FY | exact | pass |
| | `IntpQuadraticNonuniform2.degenerateTriangle` (a collinear triangle next to two proper ones; both constructors) | exact | pass after the fix below |
| | `IntpQuadraticNonuniform2.deviation.meshFlags` | deviation (#337) | 1940 of 2000 deviate |
| `IntpThinPlateSpline2.h` | `IntpThinPlateSpline2.evaluate` (`IsInitialized`, `operator()` at every sample and at 3 other points, `ComputeFunctional`; smoothing 0 and > 0; with and without the unit-square transform; a flat axis, #191) | tolerance 1e-11 (`std::log`, see below) | pass, max scaled error 7.96e-13 |
| | `IntpThinPlateSpline2.construct.invalid` (fewer than 3 points, negative smoothing) | exact, throw parity | pass, 1488 throws |
| `IntpThinPlateSpline3.h` | `IntpThinPlateSpline3.evaluate` (the same outputs; kernel `-|t|`; every 5th record has all samples -0, so the zero-seeded sums decide the signs) | exact | pass, 20048 reals |
| | `IntpThinPlateSpline3.construct.invalid` | exact, throw parity | pass, 1588 throws |
| `IntpSphere2.h` | `IntpSphere2.getSphericalCoordinates` (unit vectors, exactly the poles, just inside the poles, non-unit vectors with abs(z) >= 1) | tolerance 1e-12 (`std::atan2`, `std::acos`) | pass, 3848 of 4000 bit-identical, max 2.2e-16 |
| | `IntpSphere2.evaluate.sortedMesh` (the constructor and `operator()` replayed over the sorted mesh, compared with the port's class) | exact | pass |
| | `IntpSphere2.deviation.constructorThrows` | deviation (new finding) | 2000 of 2000 deviate |
| `IntpVectorField2.h` | `IntpVectorField2.evaluate.sortedMesh` (replay, both components) | exact | pass |
| | `IntpVectorField2.deviation.constructorThrows` | deviation (new finding) | 2000 of 2000 deviate |
| | `IntpVectorField2.deviation.staleOutput` | deviation (#337) | 2000 of 2000 deviate |

### Meshes

The nonuniform interpolators take a duck-typed mesh. Upstream's
`Delaunay2Mesh<T>` / `Delaunay3Mesh<T>` number the simplices in the iteration
order of `ETManifoldMesh` / `TSManifoldMesh`'s `std::unordered_map`, and the
port in sorted `TriangleKey<true>` / `TetrahedronKey<true>` order. Two things
depend on the numbering: the start simplex of the containment walk (a point on
a shared edge is reported in whichever triangle the walk reaches first) and
`IntpQuadraticNonuniform2::EstimateDerivatives`, which sums the triangle
normals per vertex in triangle order. Worse, upstream's numbering is not
reproducible between runs of the same executable: `Delaunay2<T>` keeps a
`std::unordered_set<Triangle*>` (`TrianglePtrSet`, Delaunay2.h:1351), so the
insertion history, and with it the map's iteration order, follows heap
addresses. A first version of the case file that chose query triangles by
upstream index produced different golden inputs on two consecutive runs.

The cases therefore use two meshes:

* the real `Delaunay2Mesh<double>` / `Delaunay3Mesh<double>`, restricted to
  what is independent of the numbering: queries strictly inside a simplex
  (all exact barycentrics positive) or outside the hull, and the derivative
  constructor of `IntpQuadraticNonuniform2`, whose preprocessing is per
  triangle;
* `SortedMesh2` / `SortedMesh3` in the case file: upstream's own
  `Delaunay2<double>` / `Delaunay3<double>` with the simplices presented in
  the port's sorted order and every search started at simplex 0 of that order;
  the barycentrics are `Delaunay2Mesh<T>::GetBarycentrics` verbatim (exact
  `BSRational`), the search is upstream's `GetContainingTriangle` /
  `GetContainingTetrahedron`. The replay uses the port's own
  `Delaunay2Mesh` / `Delaunay3Mesh`, which already number in sorted order, so
  the exact agreement of the sortedMesh cases also checks that claim.

Point sets pass v09's capped probe (`Sound2` / `Sound3`: upstream's
`IntrinsicsVector2/3` seed agrees with the exact orientation, no repeated 3D
vertex), which keeps findings #391 and #283 of `Delaunay2/3` out of this
group. Query simplices are chosen in sorted order.

### Generators and the branches they reach (deep run, 2000 records per case)

| measurement | value |
| --- | --- |
| `IntpTrilinear3` query coordinates below min / above max / on a grid node | 2092 / 2930 / 8322 of 18000 |
| grid modes | integer lattice, dyadic min and spacing, uniform, affine (trilinear) or cubic (tricubic) data, signed-zero samples; spacings 0.25..3 |
| `ComputeDerivative` branches (uniform / nonuniform): `s1==s2`, generic, `s2==s3`, `s0==s1`, both | 3298/1819, 7171/8605, 302/246, 301/247, 26/20 |
| `IntpBSplineUniform1` `GetKey`: t <= tmin / interior / t >= tmax | 3263 / 5570 / 3167; cache modes 698 / 655 / 647; `ctZero = -0` in 1022 records |
| `IntpLinearNonuniform2.sortedMesh` queries interior / edge / vertex / outside | 2844 / 2106 / 2347 / 4703 |
| `IntpLinearNonuniform3.sortedMesh` queries interior / face / edge / vertex / outside | 2690 / 1415 / 1576 / 1984 / 4335 |
| `IntpQuadraticNonuniform2` valid queries (fromDerivatives / fromSpatialDelta, sortedMesh) | 6103 / 6132 of 10000 |
| `IntpQuadraticNonuniform2` closest-subtriangle fallback (#337 item 2, preserved) | 321 and 316 of those valid queries, compared exactly; no degenerate subtriangle occurred |
| `IntpThinPlateSpline2` / `3`: initialized, smoothing 0, transformed | 1715 / 1712, 940 / 974, 1114 / 987 |
| `IntpSphere2.sortedMesh` / `IntpVectorField2.sortedMesh` valid queries | 7018 / 6145 of 10000 |

## Tolerances

**`IntpThinPlateSpline2.evaluate`: 1e-11, with a conditioning restriction.**
The kernel is `t^2 log(t^2)`. MSVC's `std::log` and V8's `Math.log` differ in
the last bit on 5594 of the 153924 kernel arguments of a 2000-record probe
(3.6%), and the constructor inverts `A = M + lambda I` and
`Q = B^T A^-1 B`, which amplifies that by their conditioning. Before the
restriction the deep run had 745 of 2000 records not bit-identical, errors up
to 1.6e-8, and six records whose `IsInitialized` differed: exactly singular
configurations (collinear samples, so `Q` is singular in exact arithmetic;
repeated samples with smoothing 0) in which Gaussian elimination's
"invertible" is decided by whether a pivot rounds to exactly zero.

Root cause proven by substitution: a temporary probe case (not committed)
recorded upstream's outputs together with every `std::log` argument and value
its computation uses; replacing `Math.log` by a lookup of the MSVC values made
the port reproduce **all 2000 probe records bit for bit**, including the
invertibility flips (1302 of 2000 were bit-identical with V8's log). The
port's thin-plate arithmetic (`GMatrix` inverse, `MultiplyATB`, the
accumulations) is therefore upstream's; `IntpThinPlateSpline3`, which runs the
same solver with the arithmetic-only kernel `-|t|`, is exact on all 2000
records, ill-conditioned ones included.

The generator now accepts only records with `cond1(A) * cond1(Q) <= 1e4`,
computed on the C++ side from upstream's own matrices (capped rejection, 64
attempts, fallback the best-conditioned candidate). On the deep run the
measured maximum scaled error is 7.96e-13 (4655 of 18984 real outputs
inexact); the replay compares at 1e-11. Records with conditioning between 1e4
and 1e5 reached 1.04e-11 in an earlier run, which is why the threshold is 1e4.
The flat-axis mode is always transformed and takes the #191 NaN path, which is
exact (both sides report not initialized).

**`IntpSphere2.getSphericalCoordinates`: default 1e-12** (`std::atan2`,
`std::acos`). 3848 of 4000 outputs bit-identical, maximum 2.2e-16. The pole
branches (`z >= 1`, `z <= -1`) return constants and are exact.

No other case needs a tolerance. The one other disagreement class met,
`IntpQuadraticNonuniform2` over the real `Delaunay2Mesh<double>` with the
spatialDelta constructor (at most 9.1e-15 on 450 of 11571 outputs, caused by
the normal-summation order), is not a case because its goldens are not
reproducible between runs; see "Not covered".

## Port defects fixed

**`src/IntpQuadraticNonuniform2.ts`, `processTriangles`: the Inscribe center of
a collinear triangle.** Upstream passes a value-initialized `Circle2` to
`Inscribe` and ignores the returned `bool`. `Inscribe` (ContScribeCircle2.h)
stores `circle.center = len21 * v0 + len20 * v1 + len10 * v2` (edge lengths
normalized by the perimeter) *before* it tests the radius, so a collinear
triangle with a positive perimeter keeps that weighted point and only a
triangle with three coincident vertices keeps `(0,0)`. The port called
`inscribeCircle2`, which returns `null` for every degenerate triangle, and left
`(0,0)` in both situations (its comment, and the summary of finding #337, said
upstream leaves `(0,0)` for any degenerate triangle). A triangle's center
enters its neighbours' cross-edge intersections and Bezier coefficients, so on
a mesh where a collinear triangle shares edges with proper ones the
neighbours' values and gradients were wrong: `IntpQuadraticNonuniform2.
degenerateTriangle` disagreed on 4 of 20 golden records with scaled errors up
to 0.67. A Delaunay triangulation never contains such a triangle, but the mesh
is duck-typed and a hand-built mesh can. The fix computes the center exactly
as `Inscribe` leaves it (`inscribedCenter`: `inscribeCircle2`, and when that
fails the same weighted combination, or `(0,0)` for a zero perimeter); all
2000 deep records are now bit-identical. Regression test:
`test/IntpQuadraticNonuniform2DegenerateTriangle.test.ts` pins record 7 bit for
bit, asserts the pre-fix values differ, and checks the collinear center and
the coincident `(0,0)`.

No other port defect was found; every accumulation order, seed, division
form and `std::max`/`std::min` substitute on the covered paths already
matched (see "Sensitivity").

## Deliberate deviations demonstrated

* **`IntpLinearNonuniform2.deviation.getIndices`,
  `IntpLinearNonuniform3.deviation.getIndices` (#135).** Upstream discards
  `GetIndices`' `bool` and blends `F[0]` with the barycentrics, returning
  `true`; the port returns invalid. The failure is reachable only through a
  mesh whose `GetIndices` fails while the search and `GetBarycentrics` succeed
  (`NoIndicesMesh2/3` in the case file, the same wrapper in the replay); with
  a consistent mesh `GetBarycentrics` fails first. 1955 and 1966 of 2000
  records deviate; the rest have all queries outside the hull.
* **`IntpQuadraticNonuniform2.deviation.meshFlags` (#337).** One per-triangle
  accessor (`GetVertices`, `GetIndices`, `GetAdjacencies` or
  `GetBarycentrics`, by record) fails for every triangle; upstream computes
  with the zero-filled outputs, the port skips the triangle or reports the
  query invalid. 1940 of 2000 records deviate.
* **`IntpVectorField2.deviation.staleOutput` (#337).** On a failed query
  upstream leaves the caller's output untouched and the port returns `(0,0)`;
  2000 of 2000 records deviate. On the finding's wording: the x-interpolation
  fails before it writes `F`, and the y-interpolation shares the mesh and
  fails on the same points, so upstream never writes half of the output;
  *both* components stay stale. Because the real class cannot be constructed
  (next item), the case runs upstream's `operator()` expression
  (`(*mXInterp)(...) && (*mYInterp)(...)`) over the sortedMesh replay.
* **`IntpSphere2.deviation.constructorThrows`,
  `IntpVectorField2.deviation.constructorThrows` (new finding, below).**
  Every construction of the real upstream classes throws; the port's classes
  work. 2000 of 2000 records deviate in each. The sortedMesh cases show that
  the port's classes compute exactly what upstream's constructor and
  `operator()` compute once the mesh is built after the triangulation.

Preserved upstream behaviour compared bit for bit: the #69 clamping
(extrapolation of the boundary cell outside the domain; trilinear holds the
boundary value above the maximum), the #191 NaN coordinates of a zero
coordinate range, the #191 `ComputeFunctional` discontinuity at `lambda = 0`
(smoothing 0 and positive smoothing both covered), the closest-subtriangle
fallback of `IntpQuadraticNonuniform2`, `Inscribe`'s ignored failure (with the
fix above), and the missing rank test of the thin-plate splines (one deep
`IntpThinPlateSpline3` record with a repeated sample and smoothing 0 reports
initialized and misses its samples by 36, identically on both sides).

## Independent reference checks (agreement is not correctness)

Run once on the deep-run outputs, which are the C++ values:

| check | result |
| --- | --- |
| trilinear on affine data, queries in [min, max) | value 1.8e-15, gradient 2.2e-15 |
| tricubic Catmull-Rom at grid nodes (211 queries) | reproduces the sample to 9.4e-16 |
| tricubic Catmull-Rom on quadratic data (fresh 6^3 grid, port) | 1.2e-14; d/dx against a central difference (h = 1e-3) 5.6e-12 |
| tricubic Catmull-Rom on cubic data | **not reproduced**: 9.6e-2 on a cubic with unit coefficients. Catmull-Rom estimates the node tangent as `(f(i+1) - f(i-1))/2`, which is exact for quadratics only; this is the scheme's precision, not a defect |
| Akima (nonuniform) at the nodes | 5.3e-15; d/dx against a central difference with h = 1e-3 of the local interval 2.8e-5 relative (the quotient's own truncation) |
| Akima on strictly increasing data | **not monotone** in 245 of 400 records (overshoot near sharp slope changes, e.g. X = -1, -0.75, 1.125, 2.375, F = 0.556, 2.247, 3.912, 3.982 rises to 4.045 before the last node). Akima's weights do not preserve monotonicity; the formula matches the published one |
| B-spline with constant controls (partition of unity), 3996 values | 9.1e-14; derivatives against central differences (degree >= 2) 1.1e-9 |
| linear nonuniform on affine data, valid queries (2D / 3D) | 9.8e-15 / 1.8e-15 |
| thin-plate spline, smoothing 0, at the samples | 2D 2.7e-13 over 792 splines; 3D the rank-deficient record above (36), every other record below 1e-9 |
| trilinear derivative at exactly `xMax` | 0 (the clamped cell collapses to one sample), see the suspects |

## Sensitivity

Recomputed from the recorded inputs with one grouping changed, against the
C++ outputs of the deep run (the replica with upstream's grouping reproduces
every record first):

* `IntpLinearNonuniform2`: `b0*F0 + (b1*F1 + b2*F2)` instead of
  `(b0*F0 + b1*F1) + b2*F2` would differ on 245 of 2000 records.
* `IntpTrilinear3`: `P*(Q*(R*F))` instead of `((P*Q)*R)*F` would differ on 490
  of 2000 records.
* `IntpAkima1::Polynomial`: the power form instead of Horner's would differ on
  1070 of 2000 `IntpAkimaNonuniform1` records, and `Math.max`/`Math.min` in
  place of `std::max`/`std::min` for the clamp on 9 (the -0 query of the
  signed-zero lattice mode).

The B-spline `ctZero` sign is drawn per record (+0 / -0, 1022 records with
-0) and the controls include signed zeros, so the accumulation seed is
exercised in both signs; how many records a literal-`0` seed would change was
not measured.

## Not covered

* **`IntpBSplineUniform` with a degree-0 axis.** Upstream's `ComputePowers`
  resizes `powerDSDT` to one element and writes `powerDSDT[1]` (#135, fixed in
  the port). On this MSVC x64 build that write is not benign: constructing a
  single `IntpBSplineUniform<double, Controls, 1>` of degree 0 terminates the
  process with `STATUS_HEAP_CORRUPTION` (0xC0000374), verified with a
  one-case probe. Every class calls `ComputePowers` in its constructor, so no
  degree-0 record can exist and no `deviation` case is possible. The port's
  guard is covered by `test/IntpBSplineUniform.test.ts`. The `A` extraction
  read past the degree (#135) is unreachable on every covered degree (the
  leading coefficient of each `Q` polynomial is `(-1)^k C(d,k)/d!`, never
  zero), so it cannot be demonstrated either.
* **`IntpQuadraticNonuniform2` spatialDelta constructor over the real
  `Delaunay2Mesh<double>`.** Its outputs depend on upstream's triangle
  numbering, which is not reproducible between runs (see "Meshes"; 3 of 2000
  records changed bits between two runs), so it cannot have a golden. It was
  within 9.1e-15 of the port; the same code is compared exactly over
  `SortedMesh2`.
* **The real `IntpSphere2<T>` / `IntpVectorField2<T>` `operator()`.** They
  cannot be constructed (below); covered by the replays.
* **`IntpThinPlateSpline2` with `cond1(A) * cond1(Q) > 1e4`.** The rounding of
  `std::log` decides the result there, even `IsInitialized`; the substitution
  experiment above shows the port is exact given MSVC's log values.
* **The deprecated `IntpSphere2<InputType, ComputeType, RationalType>` and
  `IntpVectorField2<InputType, ComputeType, RationalType>`**, not ported
  (their `Delaunay2` / `Delaunay2Mesh` specializations are not ported).
* **`GetF()` / `GetX()` accessors** return the caller's array; nothing is
  computed. `IntpAkima1` has no public `GetPolynomial`; its polynomials are
  compared through a subclass on both sides.
* **Upstream undefined inputs**: grid queries whose index overflows `int32_t`
  in `static_cast<int32_t>(xIndex)` are kept out of the generators (queries
  stay within four cells of the domain).

## Upstream bug suspects

1. **`IntpSphere2<T>` and `IntpVectorField2<T>` cannot be constructed (new,
   result-corrupting: every call throws).** Both declare
   `Delaunay2<T> mDelaunay; TriangleMesh mMesh;` and initialize
   `mMesh(mDelaunay)` in the member initializer list, but compute the
   triangulation only in the constructor body (`mDelaunay(mWrapAngles)`,
   `mDelaunay(numPoints, domain)`). `Delaunay2Mesh<T>`'s constructor
   `LogAssert`s `mDelaunay->GetDimension() == 2`, and a default-constructed
   `Delaunay2<T>` has dimension 0, so every construction throws "Invalid
   Delaunay dimension." (2000 of 2000 deep records for each class). The
   deprecated three-parameter specializations have the same initializer order
   and their `Delaunay2Mesh` has the same assert, so by reading they fail the
   same way (not run: not ported). The port builds the mesh after the
   triangulation.
2. **`Delaunay2<T>`'s triangle numbering depends on heap addresses (new for
   this header; same class as #325 and #290).** `Delaunay2<T>` uses
   `std::unordered_set<Triangle*>` (`TrianglePtrSet`) during insertion, so the
   iteration order of the triangle map, and therefore `GetIndices()`,
   `Delaunay2Mesh`'s numbering and the default start of
   `GetContainingTriangle`, can change between runs of one executable. For
   `IntpQuadraticNonuniform2`'s spatialDelta constructor this makes the
   estimated derivatives, and the interpolated values, vary in the last bits
   from run to run (observed on 3 of 2000 records).
3. **`IntpTrilinear3`'s derivatives are 0 at exactly the domain maximum
   (minor, extends #69).** At `x = xMax` the cell index clamps to the last
   sample, both stencil entries become that sample, and the first derivative
   is `0` instead of the last cell's slope. Preserved and compared.
4. **Finding #135's degree-0 overrun is fatal on MSVC** (evidence for the
   existing entry): 0xC0000374 on construction, see "Not covered".
5. **Finding #337's wording**: `Inscribe` leaves `(0,0)` only for a zero
   perimeter (a collinear triangle gets the weighted point), and
   `IntpVectorField2` leaves both output components stale rather than half of
   the output.
