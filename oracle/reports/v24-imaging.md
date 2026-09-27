# Verify group 24 (imaging) against the MSVC build of upstream GTE

Family `v24-imaging`, 33 cases, 20 golden records each (golden file 752 KB).
Deep run `npm run oracle:deep -- 2000 v24-imaging` (66 000 records, 2591 of
them throw records): **every case passes** (37 tests including the coverage
check and three independent-reference tests; generation and replay together
about 10 s). Outside the two deviation cases the deep run compares 898 213
floating-point outputs, **all bit-identical** to the MSVC build; every other
output is an integer, a boolean or a 32-bit digest and compares exactly. No
case calls the C math library, so every case is declared `{ exact: true }`
and there is no tolerance anywhere in the family. The goldens and the deep
run are deterministic (regenerated, byte-identical).

Headers: `FastGaussianBlur1.h`, `FastGaussianBlur2.h`, `FastGaussianBlur3.h`,
`Image.h`, `MarchingCubes.h`, `PdeFilter.h`, `AdaptiveSkeletonClimbing2.h`,
`FastMarch.h`, `Histogram.h`, `Image2.h`.

**One upstream defect was found and fixed in the port**
(`AdaptiveSkeletonClimbing2` pairs the crossings of a saddle cell by a
determinant that ignores the level and swaps the two pairings), demonstrated
by a `deviation` case. No translation defect was found in the port.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `FastGaussianBlur1.h` | `.execute.double`, `.float`, `.int32`, `.int16` | exact | 100 261 reals, 100 % |
| `FastGaussianBlur2.h` | `.execute.double`, `.float`, `.int32`, `.int16` | exact | 162 631 reals, 100 % |
| `FastGaussianBlur3.h` | `.execute.double`, `.float`, `.int32`, `.int16` | exact | 217 615 reals, 100 % |
| `Image.h` | `Image.access` | exact | 10 911 reals + indices, 100 % |
| `Image2.h` | `Image2.access`, `Image2.neighborhoods` | exact | 21 400 reals + indices, 100 % |
| | `Image2.neighborhoods.wrap` | deviation (#64) | 2000 of 2000 records deviate |
| `Histogram.h` | `.int.noRescaling`, `.int.rescaled`, `.double`, `.float`, `.incremental` | exact | discrete only, 100 % |
| | `Histogram.invalid` | exact, throw parity | 1202 records throw on both sides |
| `PdeFilter.h` | `PdeFilter.construct` | exact | 12 000 reals, 100 % |
| `FastMarch.h` | `FastMarch.march1`, `.march1.nanSpeed`, `FastMarch.accessors` | exact | 149 997 reals, 100 % |
| `MarchingCubes.h` | `MarchingCubes.table` (all 256 entries in the 20 committed records) | exact | discrete only, 100 % |
| `AdaptiveSkeletonClimbing2.h` | `.extract`, `.extract.types`, `.extract.saddle`, `.extract.large` | exact | 219 926 reals + digests, 100 %; 374 throw records |
| | `.invalid` | exact, throw parity | 1015 records throw on both sides |
| | `.extract.saddlePairing` | deviation (this report) | 2000 of 2000 records deviate |

### How the cases reach the branches

* **FastGaussianBlur.** The upstream algorithm is the explicit
  scale-space step `center + logBase * (xsum [+ ysum [+ zsum]])` with linear
  interpolation at `x +- scale` (not a recursive Deriche filter). Scales are
  drawn uniform, integer and half-integer (`delta = 0` exactly), past the
  image bound (the boundary branches everywhere), tiny, and 0; `logBase`
  uniform, dyadic, negative and 0; images uniform, lattice, constant and a
  single spike; bounds 1..24, 1..8 squared, 1..5 cubed. Every `T` upstream
  documents is instantiated (`double`, `float`, `int32_t`, `int16_t`); the
  port writes into the matching typed array, so `static_cast<T>` is the
  typed-array store (round to float, truncate to integer). Integer images stay
  in [-1000, 1000] so the result fits `int16_t` (out-of-range conversion is
  undefined behaviour). Independent reference in the replay: a constant image
  is returned unchanged on every record that has one (exact zero second
  differences).
* **Image, Image2.** Dimensions 0..4 (Image) with nonpositive dimensions one
  time in six (the image stays empty), construction directly or by
  `Reconstruct` over a previous image, `GetIndex` on in-range, out-of-range
  and negative coordinates, `GetCoordinates` up to twice the pixel count,
  `operator[]`/`operator()` in both forms, the clamping `Get` in all forms
  (including `Image::Get(size_t(-1))` reached through a base reference, since
  `Image2`'s overloads hide it), reads of value-initialized pixels, and every
  neighbourhood table (relative and absolute, 1-D and 2-tuple) at pixels
  inside, on and one outside the border. `size_t` results are emitted as
  signed 64-bit integers, which is how the port represents them; the wrapped
  values of #64 are compared as such in the deviation case.
* **Histogram.** Integer samples with ties, full-range `int32_t`, constant
  and two-valued images and values around `[0, numBuckets)`; real samples
  uniform, lattice, signed zeros, dyadic, spread over 80 binades and with a
  NaN first sample (every comparison with the running extremes is false, so
  upstream takes the constant branch); the float constructor on float values.
  Tail amounts: generic, exact fractions `k/N` (where `fl(k/N)*N` truncates
  to `k - 1`), 0, 1, above 1 (`GetLowerTail` returns `numBuckets`,
  `GetUpperTail` -1) and negative. The maximum sample lands in bucket `B-2`
  (finding #436, preserved) on 48, 36 and 28 deep records of the double, float
  and rescaled-int cases, identically on both sides. Independent reference:
  the buckets plus the two excess counts equal the sample count on every
  record.
* **PdeFilter** is abstract (its filters are group 25). A case-file subclass,
  mirrored in the replay, exposes the base's own code: the range scan, the
  four scale types plus the unnamed value 4 (the switch falls through,
  offset and scale stay 0), the constant image, signed zeros, the
  `PRESERVE_ZERO` tie `max == -min`, infinities and NaNs, the accessors and
  the order of `Update`'s three hooks.
* **FastMarch** is abstract (`FastMarch2/3` are group 25). A 1-D subclass,
  defined identically on both sides, drives both constructors (seeds on
  zero-speed pixels, duplicate seeds, zero, negative and NaN per-pixel speeds,
  constant speeds 0 and negative), the protected `MinHeap<size_t, Real>` and
  trial records, every classifier, `GetInterior`, `GetBoundary` and
  `GetTimeExtremes` (including the all-invalid `(max, -max)` answer and the
  redundant re-test of finding #52), with lattice speeds so that tied arrival
  times make the heap's tie order decide the removal sequence.
  `FastMarch.accessors` sets arbitrary times (both zeros, `+-max`,
  infinities, NaN, negatives) first.
* **MarchingCubes.** Record `r` covers entries `13*(r mod 20)` onwards, so
  the committed goldens hold the whole table: per entry the counts, the 12
  vertex pairs and 5 triples (padding included, packed), a digest of the
  entry's 41 values in the flat `GetTable()` storage and of the
  `GetPrebuiltTable()` row, and the configuration name; entries 256..259 and
  one far out of range check `GetConfigurationType`'s `""`. Independent
  reference (all 256 entries): the vertex pairs are exactly the cube edges
  whose corners differ in sign, no directed triangle edge is used twice (a
  consistently oriented surface), and the prebuilt table equals the
  constructed one.
* **AdaptiveSkeletonClimbing2.** Images of 3x3 to 17x17 recorded pixels
  (random small values, ramps, paraboloids, constants, signed checkerboards,
  full-range values), levels between samples, on samples (221 + 153 deep
  records throw on both sides, every one of them with the level on a sample
  value, i.e. `GetInterp`'s assertion), 0, NaN, outside the range and on
  quarter-integers, depth -2..N+1, one or two extractions on the same
  object, and every pixel type upstream allows (`int8_t` to `uint32_t`).
  `GetComponents`' default `LogError` (a rectangle type other than 3, 5, 6,
  9, 10, 12, 15) was never reached.
  `.extract.large` builds 33x33 and 65x65 images from a recorded
  closed-form formula (deep merge trees, big merged rectangles) and emits
  counts and FNV-1a digests of the vertex bit patterns and edge indices; all
  other cases emit the vertex and edge lists verbatim, before and after
  `MakeUnique` (upstream's order is the deterministic `GetRectangles`
  recursion and first-occurrence numbering, so nothing is re-sorted).
  `.extract.saddle` puts sound saddle cells (both disjoint pairings and the
  plus sign at level 0) in every record. Independent reference over every
  extraction: each vertex on a grid edge interpolates to the level, after
  `MakeUnique` every border vertex has degree 1 and every interior vertex 2
  (4 for a branch point), and every saddle-cell segment cuts off a corner on
  the other side of the level than the interpolant's saddle value.

## Port defects fixed

None of translation. The one change to `src/` fixes an upstream defect (next
section).

**Sensitivity.** Recomputed from the deep-run inputs: associating
`xsum + (ysum + zsum)` in `FastGaussianBlur3` would change 122 of 2000
records, distributing `logBase` over the sums 448 (2-D) and 559 (3-D) of
2000; `(B-1)*(s-min)/d` in place of `mult*(s-min)` changes a bucket on 55 of
2000 `Histogram.double` records. The pre-#543 `MinHeap` (`a <= b` derived as
`!(b < a)`) disagrees on 9 of the 20 committed `FastMarch.march1.nanSpeed`
records (39 of 2000 in the first deep run of the mixed generator, before the
v15 fix reached `main`), so the NaN semantics of v15 are now pinned through
`FastMarch` too.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `Image2.neighborhoods.wrap` | [#64](https://github.com/gradientspaceai/gtengine-js/issues/64) | 2000 of 2000 records: at a pixel on the xmin or ymin edge upstream's `size_t` absolute neighbours wrap to `SIZE_MAX` (emitted unsigned), the port's are negative. `Image2.neighborhoods` compares the same tables as signed values on the full generator and agrees everywhere. |
| `AdaptiveSkeletonClimbing2.extract.saddlePairing` | this report (upstream suspect 1) | 2000 of 2000 records; each has a saddle cell on which upstream's pairing contradicts the bilinear interpolant (rejection on the exact predicate, capped at 256 attempts with a fixed defective fallback). The replay's independent saddle check passes on the port's output for every record. |

Preserved findings reached and agreeing bit for bit: #436 (Histogram maximum
in bucket `B-2`), #52 (`GetEdge` never returns -1; `GetTimeExtremes`
re-test), #60 (the `PdeFilter` base keeps `mMin` at the data value on a
constant image). #436's FastGaussianBlur2/3 pointer cleanup and #60's
Neumann/border items concern state that no public function exposes or the
group-25 filters.

## Not covered

| header / entry point | reason |
| --- | --- |
| `FastGaussianBlur*` with a negative or huge scale | a negative scale reads `input[-1]`, a scale beyond `INT32_MAX` makes `static_cast<int32_t>(floor(...))` undefined; results outside `T`'s range are undefined for integer `T`. |
| `Image`/`Image2` copy and move construction and assignment | not ported (C++ value semantics). |
| `Image::Get` on an empty image | `mPixels.front()` of an empty vector is undefined behaviour. |
| `Image2` with `GTE_THROW_ON_IMAGE2_ERRORS` | compile-time option the port does not offer; the default (unchecked) build is compared. |
| `Histogram` with a NaN other than first, infinities, a subnormal range, `Insert` out of range | `static_cast<int32_t>(NaN)` and unbounded `++mBuckets[index]` are undefined behaviour (#436 item 2); `Insert` is documented unchecked. |
| `PdeFilter1/2/3`, `FastMarch2/3` | group 25; only the bases' own code is covered here, through case-file subclasses. |
| `MarchingCubes` with `IndexType` other than `int32_t` | the port's index type is `number`; the table values are below 64. |
| `AdaptiveSkeletonClimbing2` with `N < 0` | `1 << N` is undefined behaviour; `N = 0` is the throw-parity case. |
| `AdaptiveSkeletonClimbing2<uint32_t>` pixels above 2^31 | `det = i00*i11 - i01*i10` overflows `int64_t` (suspect 2). |
| `AdaptiveSkeletonClimbing2::PrintRectangles` | private debugging output, not ported. |

## Upstream bug suspects

**1. `AdaptiveSkeletonClimbing2.h` `GetComponents`, case 15: the saddle-cell
pairing ignores the level and is swapped** (result-corrupting; fixed in the
port). With all four edges of a unit cell crossed, upstream pairs the four
level-set points by the sign of `det = i00*i11 - i01*i10`: `det > 0` joins
the xmin point to the ymin point and the xmax point to the ymax point, which
cuts off corners 00 and 11. For the bilinear interpolant the decision is the
sign of `dg = (i00-L)(i11-L) - (i01-L)(i10-L) = det - L*(i00+i11-i01-i10)`:
`dg > 0` puts the saddle on the side of corners 00 and 11, which are then
connected through the cell, so it is corners 10 and 01 that must be cut off.
Upstream is right exactly when `sign(det) = -sign(dg)`; at level 0 (where
`det = dg`, the only level at which its `det = 0` plus-sign branch is
meaningful) it is wrong on every non-degenerate saddle. Reproduction: the
3x3 image `10 -1 -1 / -1 10 10 / -1 10 10`, level 0.5, depth 1: cell (0,0)
has center value 4.5 > 0.5 and `det = 99 > 0`, and upstream returns the
segments (0, 0.8636)-(0.8636, 0) and (1, 0.1364)-(0.1364, 1), which cut
off the two positive corners; the first passes through (0.43, 0.43), where
the interpolant is 4.6. At level 5.5 (saddle below the level) upstream's
answer is right. The port computes `dg` exactly (BigInt, level decomposed
as `m * 2^e`) and keeps upstream's choice wherever it is right; the main
cases are restricted by the same exact predicate (BSNumber on the C++ side)
and agree bit for bit, including 5 + 2 sound plus-sign cells and hundreds of
sound saddles. The same `det` construction appears in
`AdaptiveSkeletonClimbing3.h` (face cases, L1518, 1611, 1704, 1797, ...),
which belongs to group 26 and was not examined here.

**2. `AdaptiveSkeletonClimbing2.h`: `int64_t det` overflows for `uint32_t`
images** (undefined behaviour, minor). The class documents `uint32_t` pixels;
`i00*i11` exceeds `INT64_MAX` once two pixels exceed about 3.04e9. The port
computes with BigInt. The oracle keeps `uint32_t` pixels below 2^31.
