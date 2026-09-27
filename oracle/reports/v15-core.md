# Verify group 15 (core) against the MSVC build of upstream GTE

Family `v15-core`, 57 cases, 20 golden records each. Deep run
`npm run oracle:deep -- 2000 v15-core` (114 000 records, 8564 of them throw
records): **every case passes** (58 tests including the coverage check;
generation and replay together about 4 s). Outside the 4 deviation cases the
deep run compares 436 154 floating-point outputs, 435 691 of them (99.89 %)
bit-identical to the MSVC build; the 463 others all come from the one libm
case (`Functions.libm`, max scaled error 2.2e-16). Every other output is an
integer, a boolean or a bit pattern and compares exactly. The goldens are
deterministic (regenerated twice, byte-identical).

The group is infrastructure: `Array2.h`, `Array3.h`, `Array4.h`,
`Constants.h`, `DCPQuery.h`, `FIQuery.h`, `HashCombine.h`, `IEEEBinary.h`,
`LexicoArray2.h`, `Logger.h`, `MinHeap.h`, `TIQuery.h`, `TypeTraits.h`,
`BitHacks.h`, `CurveExtractor.h`, `Functions.h`, `MeshStaticManifold2.h`,
`MeshStaticManifold3.h`, `UniqueVerticesTriangles.h`.

**Three port defects were found and fixed** (IEEEBinary64 field wrap, MinHeap
NaN comparisons, QFNumber type traits), **one MSVC code-generation defect**
was found in the reference build itself (`gte::clamp` on signed-zero ties,
worked around in the case file), and two new minor upstream suspects are
reported below.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `Constants.h` | `Constants.values` (all 15 `GTE_C_*`) | exact | 30 000 reals, 100 % |
| `Functions.h` | `Functions.arithmetic` (`sign`, `isign`, `saturate`, `sqr`, `invsqrt`) | exact | 8000 reals, 100 % |
| | `Functions.clamp` | exact | 2000 reals, 100 % |
| | `Functions.clamp.signedZero` (targeted, see below) | exact | 2000 reals, 100 % |
| | `Functions.stdMaxMin` (port-only `stdMax`/`stdMin` against MSVC `std::max`/`std::min`) | exact | 4000 reals, 100 % |
| | `Functions.libm` (`atandivpi`, `atan2divpi`, `cospi`, `sinpi`, `exp10`) | tolerance 1e-12 (`std::atan`, `std::atan2`, `std::cos`, `std::sin`, `std::exp`) | 10 000 reals, 95.37 % exact, max scaled error 2.22e-16 |
| | `Functions.fma`, `Functions.fma.ties`, `Functions.fma.special` | exact | 6000 reals, 100 % |
| | `Functions.robustSOP`, `Functions.robustDOP` | exact | 4000 reals, 100 % |
| `BitHacks.h` | `BitHacks.bits32`, `BitHacks.bits64`, `BitHacks.powerSweep` | exact | discrete only, 100 % |
| `IEEEBinary.h` | `IEEEBinary64.fields`, `IEEEBinary32.fields` | exact | 3129 reals + bit patterns, 100 % |
| | `IEEEBinary64.fromNumber`, `IEEEBinary32.fromNumber` | exact | bit patterns, 100 % |
| | `IEEEBinary64.setEncoding`, `IEEEBinary32.setEncoding` | exact | bit patterns, 100 % (after the fix below) |
| | `IEEEBinary.constants` | exact | discrete only, 100 % |
| `HashCombine.h` | `HashCombine.combine`, `HashCombine.equalities`, `HashCombine.equalities.signalingNaN` | exact | discrete only, 100 % |
| | `HashCombine.hashValue.msvc` | deviation (port contract) | 2000 of 2000 records deviate |
| `TypeTraits.h` | `TypeTraits.traits` | exact | discrete only, 100 % (after the fix below) |
| `Logger.h` | `Logger.logAssert`, `Logger.logError` | exact, throw parity | 969 + 2000 records throw on both sides |
| `MinHeap.h` | `MinHeap.sequence`, `MinHeap.sequence.nan` | exact | 32 328 reals + every internal index, 100 % (after the fix below) |
| `Array2.h` | `Array2.access`, `Array2.owned` | exact | 33 680 reals, 100 % |
| | `Array2.owned.zeroInit` | deviation (port note) | 2000 of 2000 records deviate |
| `Array3.h` | `Array3.access`, `Array3.owned` | exact | 40 461 reals, 100 % |
| `Array4.h` | `Array4.access`, `Array4.owned` | exact | 40 210 reals, 100 % |
| `LexicoArray2.h` | `LexicoArray2.access` (all four specializations) | exact | 27 042 reals, 100 % |
| `MeshStaticManifold2.h` | `.construct`, `.queries`, `.components`, `.boundaryPolygons`, `.invalidInput` | exact | discrete only, 100 %; 1563 throw records |
| | `.adjacentTriangles` | deviation (#66) | 2000 of 2000 records deviate |
| `MeshStaticManifold3.h` | `.construct`, `.queries`, `.invalidInput` | exact | discrete only, 100 %; 1677 throw records |
| | `.adjacentTetrahedra` | deviation (#66) | 2000 of 2000 records deviate |
| `UniqueVerticesTriangles.h` | `.generateIndexedTriangles`, `.removeDuplicateVertices`, `.removeUnusedVertices`, `.removeDuplicateAndUnusedVertices`, `.scalar` | exact | 170 164 reals, 100 % |
| | `.validate` | exact, throw parity | 878 of 2000 records throw on both sides |
| `CurveExtractor.h` | `CurveExtractor.extract`, `CurveExtractor.vertexEdge`, `CurveExtractor.invalidBounds` | exact | 23 140 reals, 100 %; 1477 throw records |

56 cases are declared `{ exact: true }` (the 4 deviation cases included);
`Functions.libm` is the only tolerance case, and its tolerance is the default
1e-12 because the five functions are one libm call each, scaled by an exact
constant.

### How the cases reach the branches

* **Bit patterns.** `IEEEBinary*.fields` draw from every IEEE class: signed
  zeros, the smallest and largest subnormals and normals and their
  neighbours, random subnormals and normals, powers of two, infinities,
  default quiet NaNs, quiet NaNs with payload and signaling NaNs with payload
  1, the maximal payload and random payloads. The deep run hits all ten
  `Classification` values for both widths (59 to 407 records each).
  Encodings travel as two 32-bit integers, never as doubles.
* **Neighbours.** `GetNextUp`/`GetNextDown` are also checked in the replay
  against an independent model (the encoding mapped to a monotone signed
  integer, stepped by one and mapped back), and the classification against a
  BigInt decoder.
* **BitHacks.** 32-bit inputs are 0, 1, 2, powers of two and their +-1
  neighbours, all-ones, high-bit-set, sparse and dense random patterns; each
  is fed to both the `uint32_t` and the `int32_t` overload.
  `BitHacks.powerSweep` walks every `2^k - 1, 2^k, 2^k + 1` for k = 0..63
  (four k per record) through the 64-bit functions and, on their low 32 bits,
  every 32-bit function. The replay checks each result against BigInt
  references (highest and lowest set bit, power-of-two test, rounding up and
  down).
* **fma.** Independent draws over a wide exponent range (overflow,
  underflow into the subnormals, absorption), TwoProduct cancellation
  `w = -u*v` and near-cancellation, exact halfway cases of the single
  rounding for normal and subnormal `w`, the overflow threshold
  (`MAX_NORMAL + 2^970` rounds to infinity), products 2^-20 ulp beside a tie,
  and every combination of signed zeros, infinities, NaN, `MAX_NORMAL` and
  `denorm_min`. Deep-run result classes: 3665 normal, 698 subnormal, 744
  infinite, 398 NaN, 495 zero results of both signs.
* **MinHeap.** Random sequences (up to 40 operations, capacity 0..12) of
  inserts (5803 of them into a full heap), removes (1391 from an empty heap),
  updates that raise, lower or keep a value (2207 through the null record of
  a failed insert, many through records already removed), and `GetMinimum`,
  over values with many ties and both zeros. After every operation the case
  emits the result, `GetNumElements` and `IsValid`; at the end the whole
  `mPointers` permutation (as record identities) and every `index` member.
  The replay checks every reported minimum against a live multiset.
* **Meshes.** Grid triangulations with random diagonals, triangles removed
  with probability up to 0.45 (777 of 2000 boundary-polygon records contain a
  bow-tie vertex), random relabelling, cyclic rotation and triangle order;
  closed tetrahedron and octahedron surfaces; defective meshes (a duplicated,
  a reversed, a third triangle on an edge, a repeated vertex). Tetrahedral
  meshes are Freudenthal triangulations of up to 2x2x2 cubes, positively
  oriented, thinned, relabelled by even permutations and shuffled, plus
  defective ones. `numThreads` 0..4 is drawn, so upstream's multithreaded
  adjacency update is compared with the port's single-threaded one. The
  replay checks the component partition and the boundary polygons (every
  boundary directed edge exactly once, in its triangle's direction).

### What is compared, and how

* No output of this group has an order upstream leaves unspecified: the
  mesh classes return insertion-ordered vectors, depth-first orders and
  `std::map` key orders (the port sorts its boundary-edge keys to reproduce
  the map), and `UniqueVerticesTriangles` numbers vertices by first
  occurrence. Everything is emitted verbatim, nothing is re-sorted.
* `SIZE_MAX` (upstream `invalid`) and `Number.MAX_SAFE_INTEGER` (the port's)
  are both emitted as -1. A query index of -1 is `SIZE_MAX` in C++; the
  port's explicit `0 <= v` guards take the same branch.
* `MinHeap` internals: the case file opens `private` for `MinHeap.h` alone
  (as v06 did for `BSPPolygon2.h`) and emits each slot's record as its
  position in `mRecords`; the replay reads the port's private arrays.
* `CurveExtractor` is abstract and its concrete extractors belong to group
  16, but the base has public code of its own (the `Vertex` and `Edge`
  structs, `MakeUnique`, `Convert`, the non-virtual `Extract`). A case-file
  subclass whose pure-virtual `Extract` returns a recorded rational
  vertex/edge list (built with the protected `AddEdge`/`AddVertex`, as the
  extractors do) drives all of it; the replay mirrors it. Vertices repeat in
  several rational forms (`x/y`, `kx/ky`, `-x/-y`) so the `std::map`
  equivalence, the sign normalization and the preserved finding #362
  (reversed duplicate edges survive `MakeUnique`) are all reached: 257 deep
  records keep a reversed duplicate, identically on both sides.
  `CurveExtractor.vertexEdge` draws components up to 2^31 so the products
  reach 2^62, where the port's `compareProducts` switches to BigInt.
* `HashCombine`: the port's documented contract (porting-status note, header
  of `src/HashCombine.ts`) is "deterministic, not equal to any C++
  library": 32-bit seeds, and a double hashed by folding its binary64 bits
  (-0 as +0) where MSVC's `std::hash<double>` is 64-bit FNV-1a over the
  bytes (after mapping -0 to +0) and `size_t` is 64-bit. So the raw values
  are not comparable, and `HashCombine.hashValue.msvc` is a deviation case
  that shows it. What is comparable is compared exactly:
  `HashCombine.combine` runs upstream's `HashCombine` template with a key
  type whose `std::hash` is the port's fold, from a 32-bit seed. One step
  from a seed below 2^32 agrees with the port modulo 2^32 (carries out of
  bit 31 never come back down), so the low 32 bits are emitted; a second
  step would not agree, because `seed >> 2` brings the carried bits in.
  `HashCombine.equalities` compares the equality structure with MSVC's real
  `std::hash<double>`: which argument lists hash alike (the other zero, the
  same bits) and which differ (1-bit neighbours, sign flips, other NaN
  payloads, swapped arguments). The replay also recomputes the port's hash
  with BigInt.
* `TypeTraits` is compile-time upstream; `TypeTraits.traits` emits
  `is_arbitrary_precision` and `has_division_operator` for `double`,
  `BSNumber<UIntegerAP32>`, `BSRational<UIntegerAP32>`,
  `QFNumber<double, 1>` and `QFNumber<BSRational<UIntegerAP32>, 1>` and the
  replay evaluates the port's runtime predicates on values of those types.
* `UniqueVerticesTriangles` is instantiated with `Vector3<double>` (the
  port's default key joins the components, so both equivalences are
  component-wise `==`) and with `double`. Vertices repeat exactly, as a
  -0/+0 variant (equivalent on both sides; the first-inserted representative
  is kept, sign included) and as a 1-ulp neighbour (distinct). The case file
  defines `GTL_VALIDATE_UNIQUE_VERTICES_TRIANGLES`, so the precondition
  checks run on the C++ side and `.validate` compares throw parity with the
  port's `validate = true`; every other case feeds valid inputs.

## The reference build miscompiles `gte::clamp` on signed-zero ties

The first deep run disagreed on 5 of 2000 `Functions.clamp` records, all of
them signed zeros: for `clamp(+0, -2, -0)` the C++ source returns `xmax`,
i.e. -0, and so does the port, but the MSVC build returned +0 (and -0 where
the source says +0 for `clamp(-0, NaN, +0)`). The disassembly of the case's
out-of-line wrapper shows why: MSVC 19.44 `/O2` compiles the inner
conditional `x >= xmax ? xmax : x` to `minsd xmm2, xmm0` (xmax, x), which
returns `x` whenever `xmax < x` is false, including when the operands are
equal, so on a tie between two zeros of opposite sign it returns the wrong
zero. The transformation is exact for every other input (a NaN in either
operand gives the source's answer). A stand-alone probe reproduces it with
`/O2 /fp:precise` and with `/O2 /fp:strict`; `/Od` gives the source's value.
The `__declspec(noinline)` wrapper that fixed v23's seed folding does not
help, because the out-of-line copy is compiled the same way.

The case file therefore includes `Functions.h` between
`#pragma optimize("", off)` and `#pragma optimize("", on)` (after the
standard headers it needs, before anything else that includes it). The deep
run was regenerated with and without the pragma: **exactly those 5
`Functions.clamp` records change and no other case of the family changes**,
the libm case included. `Functions.clamp.signedZero` aims half of its
records at the miscompiled configuration (7 of the 20 committed records
return the zero of the other sign), so the committed goldens hold the port
to the source's semantics there.

Exposure elsewhere: other families compile `Functions.h` at `/O2`, so in
their goldens a `gte::clamp` whose `x` ties `xmax` as the other zero, with
`x > xmin`, returns the wrong zero. Upstream calls `clamp` in
`DistLine2AlignedBox2.h` (bounds `-extent, +extent`: the tie needs
`extent = 0`, where `x <= xmin` already takes the first branch, so it is not
reachable there) and in the `clamp` overloads of `BSNumber.h` and
`BSRational.h` (through `double`).

## Port defects fixed

**`src/IEEEBinary.ts`, `IEEEBinary64.setEncoding` did not wrap.** Upstream
shifts and ORs in `uint64_t`, so a field outside its documented range wraps
modulo 2^64 (`sign = 2` contributes nothing, `biased = 4096` falls off the
top) or spills into its neighbour (a trailing field wider than 52 bits ORs
into the exponent and sign). The port built the same expression in
unbounded `bigint`, so `fromParts(2, 0, 5n)` had the "encoding" 2^64 + 5 and
the object violated its own `uint64` invariant; 11 of the first 20 records
disagreed (the port's high 32-bit word reached 2^51).
`IEEEBinary32` already wrapped through `>>> 0`. Fixed with
`BigInt.asUintN(64, ...)`, which is the identity on in-range fields;
regression tests in `test/IEEEBinary.test.ts`.

**`src/MinHeap.ts`, `a <= b` was derived as `!(b < a)`.** Upstream's
`MinHeap<KeyType, double>` compares with the built-in `<=`, which is false
whenever an operand is NaN, while `!(b < a)` is then true. A NaN value
therefore stopped sifting where upstream's keeps moving: 7 of 20 records of
`MinHeap.sequence.nan` disagreed (a different key at the root, `IsValid`
answering differently). NaN weights are reachable in the heap's users
(`CLODPolyline`, `FastMarch`, `MinimumSpanningTree`, `VertexCollapseMesh`) on
degenerate geometry. Fixed by using the built-in `<` and `<=` when no
comparator is passed; a custom `lessThan` still gets `!lessThan(b, a)`, which
is exactly `IncrementalDelaunay2`'s `RPWeight::operator<=` upstream, and may
now be paired with an explicit `lessEqual`. Regression tests in
`test/MinHeap.test.ts` (traced by hand against upstream's code; the NaN-free
behaviour is unchanged, `MinHeap.sequence` was bit-identical before and
after).

**`src/QFNumber.ts`, the type traits of `QFNumber` were true where
upstream's are false.** The port had `QFNumber` implement the
`ArbitraryPrecisionNumber` marker with `hasDivisionOperator = true`, citing
"upstream marks QFNumber as arbitrary precision in TypeTraits.h". Only a
comment there says so; `QFNumber.h` never specializes
`_is_arbitrary_precision_internal` or `_has_division_operator_internal`, and
the MSVC build reports `false`/`false` for both `QFNumber<double, 1>` and
`QFNumber<BSRational<UIntegerAP32>, 1>` (20 of 20 records disagreed). No port
code dispatches on these predicates for a `QFNumber`, so the fix changes no
numeric result: `QFNumber` no longer implements the marker, and
`isArbitraryPrecision(q)` / `hasDivisionOperator(q)` are false. The comments
in `src/QFNumber.ts` and `src/TypeTraits.ts` and the expectation in
`test/QFNumber.test.ts` were updated (`QFNumber.h` belongs to group 5).

No other port defect was found. Sensitivity of the exact cases, measured on
the deep-run inputs: an unfused `u*v + w` would differ from the MSVC `fma`
on 398, 252 and 20 records of the three fma cases; a plain `u*v +- w*z` from
`RobustSOP`/`RobustDOP` on 715 and 693; `Math.min(Math.max(x, xmin), xmax)`
from `clamp` on 1280 records of `Functions.clamp.signedZero`; and
`Math.max`/`Math.min` from `std::max`/`std::min` on 258 records of
`Functions.stdMaxMin`.

## A harness trap: V8 quiets signaling NaNs stored into a JS array

The first deep run also disagreed on 2 records of `HashCombine.equalities`:
`HashValue(fff7ffffffffffff) == HashValue(-0)` was true in the port and false
in C++. The replay read its inputs with `io.reals(4)`, which pushes them into
a JS array; V8 silences a signaling NaN stored into a double-elements array
(it keeps the hole-NaN pattern unique), so the port saw `ffffffffffffffff`,
whose 32-bit fold `hi ^ lo` is 0, the hash of zero. A double read straight
from the golden's `Float64Array` and passed as an argument keeps its bits
(checked in node). The replay now reads NaN-carrying inputs one by one with
`io.real()`, and `HashCombine.equalities.signalingNaN` pins signaling NaNs
against their quiet twins, the zeros and themselves (bit patterns
`fff7ffffffffffff`, `7ff7ffffffffffff`, `7ff0000000000001`,
`fff0000000000001`, `7ff4000000000000`); it agrees on every record, so the
port's `hashValue(...values)` rest array keeps the bits too.

## Deliberate deviations demonstrated

| case | record of the decision | what deviates |
| --- | --- | --- |
| `MeshStaticManifold2.adjacentTriangles` | [#66](https://github.com/gradientspaceai/gtengine-js/issues/66) (UPSTREAM-FINDINGS, `MeshStaticManifold2.h` item 1) | 2000 of 2000 records. Upstream returns element `[3]` of the `<v0,v1>` tuple, the neighbour across the edge *opposite* `v0`, as `adj1`, and a valid `adj0` when only `<v1,v0>` exists. The query edge is redrawn (at most 64 times) until upstream's pair differs from the documented contract computed by a scan of the input; the replay checks the port against that same scan on every record. The main `.queries` case compares upstream's return value always and `adj0` whenever `<v0,v1>` exists (4374 of 8957 deep queries), where upstream is sound, and agrees everywhere. |
| `MeshStaticManifold3.adjacentTetrahedra` | [#66](https://github.com/gradientspaceai/gtengine-js/issues/66) (item 2) | 2000 of 2000 records: upstream returns element `[2]` of the 5-tuple (a vertex index) as `adj0`. The main `.queries` case compares `FaceExists` and the return value (5464 of 8987 deep queries hit a face) and agrees everywhere. |
| `HashCombine.hashValue.msvc` | porting-status note of `HashCombine.h` (32-bit seeds, the port's own `std::hash<double>`) | 2000 of 2000 records: the raw `HashValue(a, b, c)` of MSVC's 64-bit FNV-1a `std::hash<double>` against the port's 32-bit fold. Everything comparable is compared exactly by the three other `HashCombine` cases. |
| `Array2.owned.zeroInit` | the NOTE in the constructor comment of `src/Array2.ts` (also `Array3`, `Array4`) | 2000 of 2000 records: upstream's owned `std::vector<double>` is value-initialized, so a never-written element reads 0; the port's generic storage reads `undefined`. The `.owned` cases write every element first and agree everywhere. |

The preserved findings of these headers are compared bit for bit and are
reached on purpose: #362 (`CurveExtractor::MakeUnique` keeps reversed
duplicate edges; 4 of the 20 committed records, 257 of 2000 deep ones) and
the `BitHacks.h` comment finding #67 (the code, not the comment, defines
`GetTrailingBit`; both sides agree on every input and the replay checks every
result against a BigInt lowest-set-bit reference). Finding #110 concerns
`IEEEBinary16.h`, which is in group 16, not here.

## Not covered

| header / entry point | reason |
| --- | --- |
| `DCPQuery.h`, `FIQuery.h`, `TIQuery.h` | empty primary templates (result-struct plumbing); nothing to execute. Every Dist/Intr family exercises them. |
| `TypeTraits.h` SFINAE aliases (`TraitSelector`, `IsFPType`, `IsAPType`, `IsDivisionType`, `IsNotDivisionType`) | compile-time only; the port omits them (porting-status note). The two traits are covered. |
| `Logger.h` message text | the port throws `Error(message)` without upstream's `file(function,line): ` prefix (documented); the harness compares throw parity, not messages. `GTE_ASSERT_INDIRECT` / `GTE_ERROR_INDIRECT` and `GTE_NO_EXCEPTIONS` have no port. |
| `Functions.h` float overloads | not ported (PORTING.md: `number` is binary64). `Functions.h` has no `RoundNearest`/`RoundUp`/`RoundDown`, `IsPowerOfTwo` or `Log2OfPowerOfTwo`; the last two are in `BitHacks.h` and covered there. `BitHacks.h` has no popcount. |
| `Functions.h` `GTE_DISCARD_FMA` | a compile-time switch the port does not offer (it always fuses); the default build is compared. |
| `BitHacks.h` `GTE_THROW_ON_BITHACKS_ERROR` | off by default and not ported; the default behaviour on zero and negative inputs is deterministic and is compared (nothing in `BitHacks.h` is undefined behaviour: `1 << 31` is defined in C++17). |
| `MinHeap.h` copy constructor and assignment | not ported. |
| `MinHeap.h` `Reset` of a used heap followed by `Update` through an old record | upstream keeps its `mRecords` storage (old pointers stay valid when the size does not grow), the port allocates new records; using a pre-`Reset` handle is outside the documented contract ("constant for the life of the min-heap"). `Reset` through construction is covered. |
| `Array2/3/4.h` copy/move constructors and assignment | not ported (porting-status note). |
| `MeshStaticManifold2::GetBoundaryPolygons` on nonmanifold meshes | the documented precondition is a manifold mesh; on an edge shared by three triangles or inconsistent orientation the fan walk can cycle forever, which an upstream loop cannot be capped against. Manifold meshes with bow-tie vertices are covered. Its `LogAssert(i1 < 3)` is not reached: every neighbour recorded by the class shares the walked edge, so on manifold input it cannot fire. |
| `UniqueVerticesTriangles.h` with NaN coordinates | a NaN breaks the strict weak ordering of the `std::map` key (undefined behaviour). |
| `CurveExtractor.h` with products beyond 2^63 or zero denominators | signed `int64_t` overflow in `operator<` is undefined behaviour; a zero denominator is outside the extractors' invariant. |
| `IEEEBinary.h` NaN inputs to `fromNumber` through a JS array | see the harness trap above; `fromNumber` itself is covered with every NaN class through `io.real()`. |

## Upstream bug suspects

**1. `Constants.h`: `GTE_C_INV_SQRT_2` and `GTE_C_INV_LN_10` are one ulp below
the correctly rounded values** (precision, minor). Checked against 180-digit
values computed with BigInt (Machin's formula for pi, integer square root,
`atanh` series for the logarithms): 13 of the 15 constants are the correctly
rounded doubles, but `0.7071067811865475` (1/sqrt(2) is
0.70710678118654752440..., nearest double 0.7071067811865476, which is
`Math.SQRT1_2`) and `0.43429448190325176` (1/ln 10 is
0.43429448190325182765..., nearest double 0.4342944819032518) are the next
double down: both literals truncate the 17-digit expansion instead of
rounding it. The port copies the literals (its header comment already notes
`GTE_C_INV_SQRT_2 !== Math.SQRT1_2`), and the oracle confirms all 15 are
bit-identical to the MSVC build.

**2. `TypeTraits.h`: the comment says `is_arbitrary_precision` is true for
`QFNumber`, but no header specializes it** (doc). "The value is 'true' for
BSNumber, BSRational and QFNumber, all implemented in their header files";
`QFNumber.h` has no specialization of `_is_arbitrary_precision_internal` or
`_has_division_operator_internal`, so both traits are `false` for every
`QFNumber<T, N>` (measured). Either the comment or `QFNumber.h` is wrong;
code selecting an arbitrary-precision path with `IsAPType<QFNumber<...>>`
gets the floating-point one. The port follows the code (see the fix above).

**3. `MinHeap.h`: `Update` on a record that was already removed silently
corrupts the heap** (API hazard, minor). `Remove` parks the old root's record
in the slot just past the live range with `index = last`, and the pointer
`Insert` returned stays valid. `Update(thatRecord, smaller)` sifts it up
through live slots, moving a live record out of the live range: with values
1, 2, 3 inserted and the 1 removed, `Update(record1, 0)` makes the removed
record the root again and drops the 2 from the heap for good (the next two
`Remove` calls return 0 and 3). Nothing in the header forbids passing a
removed record. Both sides agree (the deep run passes such updates on many
records); the port preserves it.

**Toolchain, not upstream:** MSVC 19.44 `/O2` (also with `/fp:strict`)
compiles `gte::clamp`'s `x >= xmax ? xmax : x` to a `minsd` that returns the
wrong zero on a signed-zero tie; see the section above.
