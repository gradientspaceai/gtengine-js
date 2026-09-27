# Verify group 25 (imaging) against the MSVC build of upstream GTE

Family `v25-imaging`, 25 cases, 20 golden records each (golden file 716 KB).
Deep run `npm run oracle:deep -- 2000 v25-imaging` (50 000 records, 3062 of
them throw records): **every case passes** (27 tests including the coverage
check and the independent-reference test; generation and replay together
about 10 s). Outside the two deviation cases the deep run compares 506 379
floating-point outputs, 503 468 of them bit-identical; the other 2911 are
all updated cells of the two GradientAnisotropic filters, whose conductances
call `std::exp`, and agree to 6.2e-15 (scaled). Every other output is an
integer, a boolean or a 32-bit digest and compares exactly. The goldens and
the deep run are deterministic (regenerated, byte-identical, including the
multithreaded rasterizer case).

Headers: `Image3.h`, `PdeFilter1.h`, `PdeFilter2.h`, `PdeFilter3.h`,
`SurfaceExtractor.h`, `TetrahedraRasterizer.h`, `CurvatureFlow2.h`,
`CurvatureFlow3.h`, `FastMarch2.h`, `FastMarch3.h`, `GaussianBlur2.h`,
`GaussianBlur3.h`, `GradientAnisotropic2.h`, `GradientAnisotropic3.h`.

**One port defect was found and fixed** (`SurfaceExtractor::Vertex`
produced a negative zero numerator, so `Convert` returned -0 where upstream
returns +0). Two deliberate port fixes (#64, #121) are demonstrated by
`deviation` cases; the third (#122) cannot be, because upstream's defective
code is undefined behaviour on every call (below).

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `Image3.h` | `Image3.access`, `Image3.neighborhoods` | exact | 26 644 reals + indices and digests, 100 % |
| | `Image3.neighborhoods.wrap` | deviation (#64) | 2000 of 2000 records deviate |
| `PdeFilter1.h` | `PdeFilter1.heat` (case-file subclass) | exact | 46 438 reals + digests, 100 % |
| `PdeFilter2.h` | through the three 2-D filters | exact | (below) |
| `PdeFilter3.h` | through the three 3-D filters | exact | (below) |
| `GaussianBlur2.h` | `GaussianBlur2.update` | exact | 25 321 reals + digests, 100 % |
| `CurvatureFlow2.h` | `CurvatureFlow2.update` | exact | 22 712 reals + digests, 100 % |
| `GradientAnisotropic2.h` | `GradientAnisotropic2.update` | 1e-12 on updated cells (`exp`), exact elsewhere | 84 583 reals, 82 950 bit-identical, max 6.2e-15 |
| `GaussianBlur3.h` | `GaussianBlur3.update` | exact | 20 142 reals + digests, 100 % |
| `CurvatureFlow3.h` | `CurvatureFlow3.update` | exact | 18 185 reals + digests, 100 % |
| `GradientAnisotropic3.h` | `GradientAnisotropic3.update` | 1e-12 on updated cells (`exp`), exact elsewhere | 72 435 reals, 71 157 bit-identical, max 6.2e-15 |
| `FastMarch2.h` | `FastMarch2.march`, `FastMarch2.accessors` | exact | 26 187 reals + digests, 100 % |
| `FastMarch3.h` | `FastMarch3.march`, `.march.contrast`, `.march.constantSpeed` | exact | 68 443 reals + digests, 100 % |
| | `FastMarch3.faces` | deviation (#121) | 2000 of 2000 records deviate |
| `SurfaceExtractor.h` | `.vertexTriangle`, `.makeUnique`, `.extract` | exact | 95 289 reals + integers, 100 % (after the fix) |
| | `SurfaceExtractor.invalid` | exact, throw parity | 1557 records throw on both sides |
| `TetrahedraRasterizer.h` | `.lattice`, `.general`, `.twice`, `.threads` | exact | discrete only (grid pairs and digests), 100 % |
| | `TetrahedraRasterizer.invalid` | exact, throw parity | 1505 records throw on both sides |

The only tolerance in the family is the default 1e-12 on the cells a
GradientAnisotropic `Update` writes (they depend on the `std::exp`
conductances, MSVC's and V8's `exp` differ by an ulp on a few arguments).
The same cases hold every never-written cell of both padded buffers, the
construction state and `mParameter` at construction to bit identity
(`outRealExact` and exact digests).

### How the cases reach the branches

* **Image3.** Dimensions -2..5 (a nonpositive one leaves the image empty),
  constructor or `Reconstruct` over a 2x3x2 image, both `GetIndex` overloads
  on coordinates up to two outside (negative ones wrap in `size_t` and are
  emitted signed, the port's value), both `GetCoordinates` overloads up to
  twice the voxel count, `operator()` in both forms, the clamping `Get` in
  both forms for reads and writes, and the base `Image::Get(i)`. Every
  neighbourhood table (6, 18, 26, corners, full; 1-D offsets and 3-tuples;
  relative and absolute) on 1..6 cubed images at voxels inside, on and one
  outside the border; the tables are emitted as digests of the signed
  values, the 6-neighbourhoods in full.
* **PDE filters.** `PdeFilter1` has no concrete filter upstream; a
  case-file subclass (explicit heat equation through `LookUp3`) is mirrored
  in the replay. The 2-D and 3-D bases are reached through the six concrete
  filters, each wrapped in a probe subclass that exposes the protected state
  (`mMin`, `mOffset`, `mScale`, `mSrc`, both padded buffers, the padded
  mask). Bounds 1..6 (2-D, 6x6 often), 1..5 (3-D, 5x5x5 often), 1..8 (1-D),
  1..5 and 1..3 for the anisotropic filters; spacings 1, 0.5, 2 and uniform;
  data uniform, positive, lattice with ties, constant, signed zeros and
  dyadic, x-only lattice quadratics, a spike, a signed checkerboard; masks
  none, random, all ones, a box, and values in {-1, 0, 1, 2}; Neumann
  (`borderValue = DBL_MAX`) and Dirichlet borders (0, lattice, uniform), so
  every image- and mask-border routine runs; the four scale types and the
  unnamed value 4; time steps stable (`GetMaximumTimeStep` for the blurs),
  dyadic, uniform, 0, negative and unstable; K = 0, 0.5, 1, 2 and uniform.
  After construction and after each of up to four `Update` calls: the source
  index, digests of both padded buffers, and a digest of every public derived
  quantity (`GetU`, `GetUx` ... `GetUzz`, `GetMask`) at every pixel.
* **FastMarch2/3.** Probe subclasses expose the heap, the trial records and
  the inverse speeds. Grids 1..7 squared and 1..5 cubed; 0..4 seeds anywhere
  (duplicates; border seeds are overwritten by the zero-speed border);
  constant speeds (lattice, uniform, 0, negative, NaN) and per-pixel speeds
  (lattice, uniform, 0, -1, NaN) and a high-contrast mix of 0.25 and 4. Each
  state emits the time, classifier and inverse-speed digests, the trial
  records (key, value, heap position) in key order, `GetInterior`,
  `GetBoundary` and `GetTimeExtremes`; each removal emits the heap minimum's
  key and value, so heap ties are compared. Upstream `Iterate` is called only
  on a nonempty heap (it indexes `mTrials[SIZE_MAX]` otherwise). ComputeTime
  branch histogram of the deep run (counted in the replay from the state
  before each call): one term 15 733, two terms with a nonnegative
  discriminant 17 612 and with the negative-discriminant fallback (#439)
  4338, three terms 19 005 and 11 012. Every tenth `FastMarch2.march` record
  is #439's reproduction (6x4, speed 1/2, seeds (4,2) and (1,2)).
  `FastMarch3.march` gives every face voxel a nonpositive speed so that
  upstream's missing face marking (#121) is invisible and the whole state
  compares; `.march.constantSpeed` uses the constructor where upstream leaves
  the faces far, seeds strictly inside, stops before the march would reach a
  face voxel and leaves the face voxels out of the outputs.
* **SurfaceExtractor** is abstract (its extractors are group 26). A case-file
  subclass, mirrored in the replay, returns recorded rational vertices and
  triangles from `Extract` and an affine gradient from `GetGradient`, which
  exposes the base's own code: the `Vertex` sign normalization, both
  comparison operators of `Vertex` and `Triangle`, the `Triangle` rotation,
  `MakeUnique`, `Convert`, the public `Extract(level, bool, ...)` overload,
  `OrientTriangles` (both `sameDir`) and `ComputeNormals`. Rational
  coordinates are drawn as a base value scaled by 1..3 with a random sign,
  a third of the vertices re-represent an earlier vertex, and a third of the
  triangles rotate an earlier one through rationally equal vertices (181
  deep records keep two rotations of one triangle, #439, preserved).
* **TetrahedraRasterizer.** `.lattice`: region `[r0, r0 + b - 1]` with bound
  b, so the grid multiplier is exactly 1 and every predicate is exact; random
  lattice tetrahedra of both orientations, flat ones, and a lattice box split
  into the six Kuhn tetrahedra (shared faces and the box diagonal put voxel
  centres exactly on faces), vertices one past the region (clipping).
  `.general`: uniform vertices, uniform or dyadic regions, bounds 2..10.
  `.twice`: a second `operator()` over the boxes the first call clipped in
  place (#439, preserved). `.threads`: 2..4 `std::thread`s with every
  tetrahedron in its own cell of the region's 2x2x2 split (see suspect 1 for
  why overlapping tetrahedra cannot be compared multithreaded). The other
  cases pass `numThreads` 0 or 1. Grids are emitted as the count, a digest
  of the (index, tetrahedron) pairs and the pairs in full up to 40.

### Independent references (replay, deep-run counts)

* Image3: every 1-D absolute table equals `GetIndex` of the matching 3-tuple
  table entry by entry, and the relative 3-tuples are the documented offsets
  (3978 voxels).
* GaussianBlur2/3 with a Neumann border and no mask preserve the interior sum
  on the first step (the ghost cells mirror the boundary and the second
  differences telescope): 282 records. On the second step the sum drifts in
  123 of 206 records, because the ghost cells are never refreshed: #60 item 2
  (preserved) seen independently of the oracle.
* CurvatureFlow2/3 leave an image whose level sets are planes (Neumann, no
  mask, data depending on x only) unchanged at every step: 204 records. With
  a small step the flow never increased the total variation (109 records).
* FastMarch: every removed pixel has a known 4-/6-neighbour when it is
  removed (30 434 removals) and a removed time never changes afterwards. On
  the #439 repro records pixel (2,2) starts at 2 and ends at
  2 (1 + sqrt(2)/2) (200 records), the defect #439 describes. The acceptance
  order is not monotone: 552 of 31 214 finite removals are smaller than the
  previous one (see the proposed note under #439).
* SurfaceExtractor: every `Triangle` is a rotation of its input with the
  minimal index first; after `MakeUnique` no two vertices are equal as
  rationals, the triples are distinct and every triangle names the same
  rational points as an input triangle (1811 records with triangles);
  `OrientTriangles`' decision agrees with the exact sign of
  `(grad(v0) + grad(v1) + grad(v2)) . ((v1 - v0) x (v2 - v0))` in dyadic
  arithmetic on all 7961 triangles.
* TetrahedraRasterizer: on the lattice records the grid equals an exact
  BigInt point-in-tetrahedron rasterization (last containing tetrahedron
  wins) on every record; on the other 6000 records upstream's grid equals the
  exact test on its own floating-point grid vertices and clipped boxes (75 325
  voxels, no difference).

## Port defects fixed

**`SurfaceExtractor.ts`, `SurfaceExtractorVertex` constructor: negative zero
numerator.** Upstream stores `int64_t` numerators and denominators and
normalizes `(n, d)` with `d <= 0` to `(-n, -d)`; an integer has no negative
zero. The port negated a zero numerator to `-0`, so `convert()` (and the
public `extract()`) returned `-0` coordinates where upstream returns `+0`,
for every vertex with a zero numerator over a negative denominator; a
caller's `0 * d` with `d < 0` produced a `-0` numerator directly as well.
Whether the group-26 extractors, which build their vertices as
`(n, d)` pairs from voxel differences, produce such vertices is for
group 26 to measure; the base class now cannot. Fixed by adding `+ 0`
to the stored numerators (and the negated denominators), which is exactly
`int64_t`'s behaviour. Size: the sign of zero only; every one of the 3149
differing outputs of the deep run (929 `makeUnique` and 1003 `extract`
records of 2000) was a `-0` against `+0`, no triangle or orientation
decision changed. Regression test in `test/SurfaceExtractor.test.ts` (fails
without the fix).

No other translation defect was found.

**Sensitivity.** Recomputed from the deep-run inputs with the port patched
to the plausible alternative (a temporary test, not committed), the cases
would have caught: `mUzz + dt * (numer / denom)` in `CurvatureFlow2`
(401 of 2000 records), `numer0 + (numer1 + numer2)` in `CurvatureFlow3`
(617), `mUzz + (dt * uxx + dt * uyy)` in `GaussianBlur2` (597),
`uxx + (uyy + uzz)` in `GaussianBlur3` (708), `average * (1 / count)` in
`PdeFilter2`'s Neumann mask border (185, 194 and 117 records of the three
2-D filters), `normal * (1 / length)` in `ComputeNormals` (644), and
`(sum + sqrt(discr)) * (1 / 3)` in `FastMarch3::ComputeTime` (8, 1415 and 23
records of the three FastMarch3 cases). `(Real)0.5 * (sum + sqrt(discr))`
and `(sum + sqrt(discr)) / 2` are one computation and cannot be
distinguished.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `Image3.neighborhoods.wrap` | [#64](https://github.com/gradientspaceai/gtengine-js/issues/64) | 2000 of 2000 records: at a voxel on the xmin, ymin or zmin face upstream's `size_t` 3-tuple neighbours wrap to `SIZE_MAX` (emitted unsigned), the port's are -1. `Image3.neighborhoods` compares the same tables as signed values on the full generator and agrees everywhere. |
| `FastMarch3.faces` | [#121](https://github.com/gradientspaceai/gtengine-js/issues/121) | 2000 of 2000 records: with the z faces at zero speed and the x/y faces positive, upstream leaves the x/y face voxels far, marches into them (their wrapped neighbours `i - 1`, `i + xBound` stay inside the grid because the z faces are never visited, so the records are deterministic) and gives them times; the port has them at zero speed from the start. `FastMarch3.march` (faces nonpositive) and `.march.constantSpeed` (march stopped before a face) agree everywhere. |

**#122 (GradientAnisotropic2/3 `ComputeParameter`) has no deviation case.**
Upstream iterates the padded coordinates `1..bound` into the unpadded
`GetUx`/`GetUy`/`GetUz`. For the last row, `GetUy(x, yBound)` reads
`F[yBound + 2]`, i.e. `Array2::mIndirect1[yBound + 2]`, one past the end of
the row-pointer vector, and dereferences whatever pointer lies there (in 3-D
`GetUz(x, y, zBound)` does the same with `Array3`'s plane pointers). That is
undefined behaviour on every construction and every `Update` (which calls
`ComputeParameter` again), so no golden record of the unmodified member can
exist. The case file replaces exactly that member for `Real = double` with an
explicit specialization implementing the port's loop over `0..bound-1`; the
constructors, `OnPreUpdate`, `OnUpdateSingle` and the `PdeFilter2/3` bases
are upstream's code, and they agree with the port to 6.2e-15 (exact before
the first `exp`). The preserved #122 behaviour (a constant image gives
`mParameter = +Infinity` and NaN updates) is reached and agrees.

Preserved findings reached and agreeing bit for bit: #60 (all three items;
item 2 also seen by the sum reference above), #123 (`CurvatureFlow2`'s 0.5),
#439 (FastMarch fallback and unused spacings, `GradientAnisotropic`'s single
spacing division and missing base `OnPreUpdate`, `MakeUnique`'s rotated
duplicates, `ClipCullAABBs` in place), #52 (`GetTimeExtremes`).

## Not covered

| header / entry point | reason |
| --- | --- |
| `Image3` copy and move construction and assignment | not ported (C++ value semantics). |
| `Image3` with `GTE_THROW_ON_IMAGE3_ERRORS` | compile-time option the port does not offer; the default (unchecked) build is compared. |
| `Image3::Get` on an empty image | `mPixels[0]` of an empty vector is undefined behaviour. |
| `GradientAnisotropic2/3::ComputeParameter` as upstream wrote it (#122) | undefined behaviour on every call (out-of-range row-pointer read and dereference); replaced by an explicit specialization in the case file, see above. |
| `PdeFilter*`, `FastMarch*`, `TetrahedraRasterizer` with `Real`/`T` = `float` | the port is double only. |
| PDE filters with a zero bound | the base constructor reads `data[0]` of an empty array. |
| `FastMarch2/3::Iterate` on an empty heap | upstream then reads `mTrials[SIZE_MAX]` (suspect 2); the port returns. No golden record can exist. |
| `FastMarch2/3` with a valid or far border pixel (`SetTime` on the border), and `FastMarch3` marches that reach a z face | upstream indexes outside the grid (`i - 1`, `i - mXYBound` wrap in `size_t`). |
| `SurfaceExtractor::Vertex` with a zero denominator | outside the contract (the extractors divide by nonzero voxel differences); upstream's `operator<` is then not a strict weak ordering, and the port's reduced-fraction key differs. |
| `SurfaceExtractor::Vertex` with numerators or denominators beyond 2^26 | upstream's `int64_t` cross products versus the port's doubles, exact below 2^53 (documented port limit); the cases use values up to 12. |
| `SurfaceExtractor<T, Real>` with `T` other than `int32_t` | `T` only types `mInputVoxels`, which the base's methods never read. |
| `TetrahedraRasterizer` multithreaded with overlapping tetrahedra | concurrent writes to one grid element; the result depends on the schedule (suspect 1). |
| `TetrahedraRasterizer` with `regionMin >= regionMax` or tetrahedra whose clipped grid box is negative | `static_cast<size_t>` of a negative double is undefined behaviour. |
| the level-surface extractors (`SurfaceExtractorMC`, `...Cubes`, `...Tetrahedra`) | group 26. |

## Upstream bug suspects

**1. `TetrahedraRasterizer.h` `MultiThreadedRasterizer`: concurrent writes to
shared voxels** (undefined behaviour, schedule-dependent result; minor). The
threads rasterize disjoint ranges of tetrahedra into one `grid` with plain
`int32_t` stores. `PointInTetrahedron` is inclusive (`> 0` rejects), so a
voxel centre on a face shared by two tetrahedra of a mesh, the normal input,
is written by both; when the tetrahedra fall into different threads, those
are unsynchronized writes to one object (a data race, undefined behaviour in
the C++ memory model), and in practice the voxel gets whichever index is
stored last, where the single-threaded path always leaves the larger index.
The documentation recommends `numThreads > 0` without a caveat. The port has
only the single-threaded path. The oracle's multithreaded case keeps every
tetrahedron in its own cell of the region so that no voxel is shared.

**2. `FastMarch2.h`, `FastMarch3.h` `Iterate`: an exhausted heap is not
detected** (undefined behaviour; minor). `mHeap.Remove(i, value)` returns
`false` on an empty heap and leaves `i = 0`; `Iterate` ignores the result,
sets `mTrials[0] = nullptr` and then tests `IsTrial(i - 1)`, i.e.
`mTrials[SIZE_MAX]`. The heap is protected and the classes offer no public
"done" test (a caller would have to scan every pixel with `IsTrial`), so a
fixed-count loop that runs past the last trial pixel reads out of range. The
port returns from `iterate()` on an empty heap (port note in
`FastMarch2.ts`/`FastMarch3.ts`); the oracle calls upstream's `Iterate` only
while the heap is nonempty.

### Proposed notes for existing findings

* **#122** (mechanism, for the `GradientAnisotropic2.h`,
  `GradientAnisotropic3.h` entry): the out-of-range access is not only one
  element past the padded buffer. `GetUy(x, yBound)` reads
  `F[yBound + 2][x + 1]`, i.e. `Array2::mIndirect1[yBound + 2]`, one past the
  end of the row-pointer vector, and dereferences that pointer for every `x`
  of the last row (3-D: `GetUz(x, y, zBound)` with `Array3`'s plane
  pointers). This happens in the constructor and in every `Update`
  (`OnPreUpdate` recomputes the parameter), so upstream's filters are
  undefined behaviour on every use. The C++ oracle (group 25) compares
  upstream's filters with only `ComputeParameter` replaced (explicit
  specialization) and finds the rest bit-identical up to `std::exp` rounding
  (6.2e-15).
* **#439** (`FastMarch2.h`, `FastMarch3.h`): the acceptance order is not
  monotone. Reproduction: `FastMarch2` 5x4, spacing 1, seed (1,1), interior
  speeds 4, 4, 1 (row y = 1) and 0.25, 0.5, 0.5 (row y = 2). The heap
  removes (2,1) at 0.25, (3,1) at 1.25, (3,2) at 2.9557 and then (2,2) at
  2.0149: when (2,1) was accepted, (2,2) took the trial time 4 of (1,2) as
  upwind data, the discriminant was negative and the fallback kept the larger
  term 4 (the correct time is 0.25 + 2 = 2.25); only after (3,2) was accepted
  did the recomputation bring it below the time already accepted. The oracle's
  deep run sees 552 decreases in 31 214 removals, identically on both sides.
