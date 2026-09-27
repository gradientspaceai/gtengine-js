# Group 28 (interpolation) — `v28-interpolation`

31 cases, 20 records each in `oracle/golden/v28-interpolation.txt`. 29 are
declared `exact` (bit-identity required) and 2 are declared deliberate
`deviation`s (issue #58). No case needs a tolerance: every path in the group
is `+ - * /` plus the truncating float-to-int conversion of a cell index.

Deep run: `npm run oracle:deep -- 2000 v28-interpolation`: 62000 records
(4535 of them recorded C++ throws), **all 32 tests pass** (31 cases and the
claim check), wall time 5.8 s including generation. Over the deep run
**every one of the 3,270,714 real outputs of the 29 exact cases is
bit-identical** to the MSVC build.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `HermiteCubic.h` | `HermiteCubic.evaluate` (orders 0-4), `.generate` (zero polynomial, then two `Generate` calls on one object), `.P` (all 4x5 `P(select, order, t)` table entries) | exact | pass |
| `HermiteQuintic.h` | `HermiteQuintic.evaluate` (orders 0-6), `.generate`, `.P` (all 6x7 entries) | exact | pass |
| `HermiteBicubic.h` | `HermiteBicubic.evaluate` (every generated `c`, orders 0..4 x 0..4), `.manual` (zero polynomial, public `c` set directly, then `Generate` overwriting it) | exact | pass |
| `HermiteBiquintic.h` | `HermiteBiquintic.evaluate` (orders 0..6 x 0..6), `.manual` | exact | pass |
| `HermiteTricubic.h` | `HermiteTricubic.evaluate` (every `c`, all 64 order triples 0..3, drawn triples up to 4), `.manual` | exact | pass |
| `HermiteTriquintic.h` | `HermiteTriquintic.evaluate` (all 216 `c`, drawn order triples 0..6), `.manual` | exact | pass |
| `IntpAkima1.h` | `IntpAkima1.evaluate.uniform1` (both `operator()` overloads, orders -1..4, clamping on both sides), `.construct.validation`, `.evaluate.signedZero` | exact | pass |
| `IntpAkimaUniform2.h` | `IntpAkimaUniform2.evaluate`, `.evaluate.separable`, `.evaluate.signedZero`, `.construct.validation` | exact | pass |
| `IntpAkimaUniform2.h` | `IntpAkimaUniform2.evaluate.maxBoundary` | deviation (#58) | 1884 of 2000 records deviate |
| `IntpAkimaUniform3.h` | `IntpAkimaUniform3.evaluate`, `.evaluate.separable`, `.evaluate.signedZero`, `.construct.validation` | exact | pass |
| `IntpAkimaUniform3.h` | `IntpAkimaUniform3.evaluate.maxBoundary` | deviation (#58) | 2000 of 2000 records deviate |
| `IntpBicubic2.h` | `IntpBicubic2.evaluate` (value, orders 0..4 x 0..4 and -1, Catmull-Rom and B-spline blending), `.construct.validation` | exact | pass |
| `IntpBilinear2.h` | `IntpBilinear2.evaluate` (value, orders 0..2 x 0..2 and -1), `.construct.validation` | exact | pass |

Every Akima/bicubic/bilinear case also emits `GetQuantity` and the computed
`GetXMax`/`GetYMax`/`GetZMax`.

Generators. Samples are uniform, integer lattice, dyadic or sparse (`+0`,
`-0`, small integers) for the Hermite classes; the grids add additively
separable, bilinear/trilinear with integer coefficients, plateau `{0, 1}` and
cubic lattice data. Grid minima are dyadic, uniform or zero; spacings dyadic
(`0.25 .. 3`), powers of two, `1` or uniform in `[0.1, 3]`. Grids run from the
minimum sizes (Akima 3x3 and 3x3x3, bicubic 3x3, bilinear 2x2) up to 6x6 and
4x4x4. Query coordinates (deep run, `IntpBicubic2.evaluate`, x and y
together): 3264 inside a cell, 2985 exactly on a node, 1689 exactly at the
maximum, 1766 above it, 1704 below the minimum, 592 equal to `-0` against a
`+0` minimum; `IntpBilinear2.evaluate` is split alike. Catmull-Rom and
B-spline blending 1000 records each. The unit-cell Hermite parameters are
uniform, dyadic `k/16`, `+0`, `-0`, `1`, or outside the cell in
`[-0.5, 1.5]`. Akima's `ComputeDerivative` branches over the 1D deep run:
`s1 == s2` 2397, all different 7871, `s2 == s3` 312, `s0 == s1` 334,
`s0 == s1 && s2 == s3` 110 windows. The validation cases draw every
parameter valid on even records and each one invalid with probability 1/2 on
odd records; 751 / 930 / 977 / 949 / 928 of 2000 deep records throw.

## Port defects fixed

None. Every output of every exact case agreed bit for bit on the first
complete run and over the deep run. The one `src/` change is documentation:
the header of `src/IntpAkimaUniform3.ts` still said the #58 sign slip was
"preserved" although the code (like `IntpAkimaUniform2.ts`) fixes it; the note
now says fixed.

Sensitivity (ORACLE.md, "a group that finds no port defect"). Each plausible
alternative below was substituted into the port and the deep run replayed;
the count is the records of 2000 on which the substitute disagrees with MSVC:

| alternative | first case to catch it |
| --- | --- |
| `HermiteTricubic` / `HermiteBicubic` evaluate `c*(xy*z)`, `c*(x*y)` | 1806 / 1743 |
| `HermiteTricubic` `v111` regrouped; `HermiteTriquintic` `v222` tail pairwise | 471 / 496 |
| `HermiteQuintic` `P3D0` as `(1-t)*t*(1-t)*t*t` | 1146 (`HermiteBiquintic.evaluate`) |
| `IntpAkimaUniform2::Construct` `A(2,2)` pairwise | 411 |
| `IntpAkimaUniform2` value operator as a power sum instead of Horner | 915 |
| `IntpAkimaUniform2::GetFX` slope `/ spacing` instead of `* (1/spacing)` | 334 |
| `IntpAkimaUniform3` polynomial `c*(x*y*z)` | 1867 |
| `IntpAkima1` `Polynomial` order 1 `x*(3*c3)` | 110 |
| `IntpBicubic2` / `IntpBilinear2` `P*(Q*F)` | 1480 / 375 |
| `IntpAkima1` clamp with `Math.max`/`Math.min` (either overload) | 658 / 781 (`.signedZero`) |
| `IntpAkimaUniform2` value operator clamp with `Math.max`/`Math.min` | 314 (`.signedZero`) |

The `.signedZero` cases exist for the last two rows: `std::max(-0, +0)` is
`-0`, and the sign survives only through Horner's `c0 + dx*(...)` with a
`-0` sample at the lower corner, which the general generators hit too rarely.
Provably indistinguishable, and therefore not a claim of these cases: the
clamp of `IntpAkimaUniform2`'s derivative operator and of both
`IntpAkimaUniform3` operators (their power sums are seeded with `+0` and the
coordinate enters only through products, so a zero's sign never reaches the
result), and `Math.trunc(-0.5) = -0` in `IntpBicubic2`/`IntpBilinear2` (`P`,
`Q` and the result are `+0`-seeded sums of the same kind).

## Deliberate deviations demonstrated

`IntpAkimaUniform2::GetFXY` and `IntpAkimaUniform3::GetFXY/GetFXZ/GetFYZ/GetFXYZ`
reuse the min-boundary one-sided stencils at the max boundaries on reflected
indices without the sign flip (issue
[#58](https://github.com/gradientspaceai/gtengine-js/issues/58)); the port
negates exactly those stencils. The main cases query only cells none of
whose corners lies on a max face, selected with upstream's own clamp and
lookup expressions (capped at 16 attempts, falling back to the minimum), or
use additively separable data with power-of-two spacings, on which every
mixed stencil is exactly zero (with other spacings `FXYZ`'s masks leave a
rounding residue of that zero, which the port negates like any other value;
the first run showed it). The `.maxBoundary` cases query the cells on the max
faces:

* 2D: 1884 of 2000 records deviate. The 116 that agree all have additively
  separable data (every 2x2 mixed difference zero), where upstream is sound.
* 3D: 2000 of 2000 deviate.
* Magnitude, against the exact bilinear/trilinear function on the lattice
  data that are bilinear/trilinear: upstream's worst relative error is
  `2.5e-2` (2D, 1116 queries) and `2.1e-1` (3D, 1122 queries); the port's is
  `1.8e-15` and `5.1e-15`. In the sound cells both sides reproduce
  bilinear/trilinear data to `1.0e-15` / `1.1e-15`.
* Confinement: reverting the port's 29 sign sites (16 negated stencils, 13
  subtracted `FXYZ` mask terms) makes both `.maxBoundary` cases agree with
  MSVC on all 2000 records while every other case still passes. So the fix
  changes nothing but those values, and the port's remaining max-face
  expressions (the `(xmax, ymax)` corner of `FXY` and the corners with an
  even number of max directions) are upstream's bit for bit.

Preserved and compared bit for bit rather than deviated from (issue
[#69](https://github.com/gradientspaceai/gtengine-js/issues/69)):
`IntpBicubic2`/`IntpBilinear2` clamp only the cell index and extrapolate with
the fractional coordinate (the below/above query modes); the stray `;;` in
`IntpBicubic2`'s B-spline `mBlend` initialiser is harmless and the port's
`1/6, -3/6, ...` constants match (1000 B-spline records per deep run).

## Independent-reference checks

Run once against the port, which the oracle shows is bit-identical to
upstream on these paths (for the Akima max faces both sides are reported
above):

| check | result |
| --- | --- |
| Hermite corner conditions: `H` and its derivatives at the corners equal the samples (2000 random sets each; triquintic 300 sets, all 27 fields at 8 corners) | cubic `1.8e-15`, quintic `8.2e-14`, bicubic `5.5e-14`, biquintic `8.3e-12`, tricubic `3.1e-13` (relative); triquintic at most `1.01 eps * sum|c P P P|`, i.e. pure cancellation in coefficients up to ~1e4 |
| `IntpAkima1` reproduces the samples at the nodes (4473 node queries) | exact at the 3059 queries on a node below the last; `1.1e-13` absolute at the 1414 on the upper end node, which is evaluated as the right end of the last cell |
| `IntpAkima1` first derivative against a central difference, `h = 1e-3 * spacing` | within 5x of `h^2 |f'''|/6 + 4 eps |f|/h` |
| `IntpAkimaUniform2` `Fxy` against a central difference of `Fx` | within 5x of the same bound |
| `IntpAkimaUniform2`/`3` reproduce bilinear/trilinear data | see the deviation section: exact to round-off except upstream's max faces |
| `IntpBicubic2` Catmull-Rom: samples at nodes (dyadic spacing) | exact |
| Catmull-Rom reproduces biquadratic data in cells whose 4x4 stencil is inside the grid | `8.3e-15` relative |
| Catmull-Rom on bicubic data | not reproduced (`0.39` relative): Catmull-Rom is exact through degree 2 only, as designed |
| B-spline blending reproduces bilinear data in the same cells | `2.3e-15`; biquadratic data are not reproduced (an approximating spline) |
| `IntpBicubic2` `(1,0)`/`(0,1)` against central differences | `3.6e-5` relative with `h = 1e-3 * spacing` (truncation); the `mInvXSpacing` scaling is right for every spacing |
| `IntpBilinear2` bilinear data inside the domain; below the minimum; above the maximum | `2.6e-15`; linear extrapolation (`2.7e-15`); holds the boundary value (`3.6e-15`), as #69 describes |
| `IntpAkima1` on nondecreasing data, 64 samples per cell | 351 of 729 records overshoot; every one violates the Fritsch-Carlson monotonicity region (an extrapolated end derivative of the wrong sign, or `alpha^2 + beta^2 > 9`) |

The last row answers "is Akima monotone on monotone data": it is not, and
upstream does not claim it; the Akima derivative is a weighted mean of the
neighbouring slopes and its boundary slopes are linear extrapolations, which
can be negative. This is a property of the method, not a defect.

## Not covered

* **`IntpAkima1`'s concrete subclasses.** The abstract base is reached
  through `IntpAkimaUniform1` only. `IntpAkimaUniform1`'s constructor and
  `Lookup`, and `IntpAkimaNonuniform1`, belong to group 29; they ran here as
  the base's vehicle and agree, but group 29 owns their coverage.
* **Trivial accessors** that return a constructor argument (`GetXBound`,
  `GetXMin`, `GetXSpacing`, `GetF`, ...) are not emitted; `GetQuantity` and the
  computed `GetXMax`/`GetYMax`/`GetZMax` are.
* **`HermiteCubic`/`HermiteQuintic` coefficients.** `c` is private upstream
  (see the suspects below), so these two are covered through `Generate` only;
  the tensor-product classes also cover a directly assigned `c`.
* **`HermiteCubic::P` / `HermiteQuintic::P` with `select` out of range**
  index past the function table (undefined behaviour); not generated.
* **Cell indices of magnitude `2^31` or more, and NaN coordinates, in
  `IntpBicubic2`/`IntpBilinear2`.** `static_cast<int32_t>(xIndex)` is
  undefined there (suspect 1 below), so the generators stay within two spacings
  of the domain.
* The **order-dependent sign of zero** in the four provably indistinguishable
  places named under "Port defects fixed".

## Upstream bug suspects

New findings only; #58 and #69 are already in `docs/UPSTREAM-FINDINGS.md`.

1. **`IntpBilinear2.h`, `IntpBicubic2.h` (and the same line in
   `IntpTrilinear3.h`, `IntpTricubic3.h`): the cell index conversion is
   undefined for large or NaN coordinates.** `int32_t ix =
   static_cast<int32_t>(xIndex)` is undefined behaviour when `xIndex` is NaN or
   outside the `int32_t` range. MSVC's `cvttsd2si` returns `INT_MIN`, which the
   following `ix < 0` clamp turns into cell 0, so a query far to the right is
   evaluated on the *first* cell. Measured with the MSVC build, `f = x + 2y` on
   `[0,2]^2` with unit spacing: `IntpBilinear2(2147483646, 1) = 4` (the held
   boundary value, #69) but `IntpBilinear2(2147483648, 1) = 2147483650`, and
   `IntpBicubic2` jumps from `4.95e27` to `-4.95e27` across the same point. The
   port's `Math.trunc` followed by the same clamp selects the last cell, which
   agrees with upstream for every `|xIndex| < 2^31`. Severity: minor (input
   range); port: defined behaviour, no change.
2. **Hermite headers, documentation and access.** `HermiteCubic.h` and
   `HermiteQuintic.h` declare `c` under `private:` right below the comment "Set
   the coefficients manually as desired", while the four tensor-product
   classes make it public; `HermiteQuintic.h` calls itself "Hermite cubic
   interpolation ... globally C1-continuous" (the quintic lattice interpolant
   is C2, as the biquintic/triquintic comments say) and lists `hermite(5, x)`
   as `Hxxxx`; `HermiteBiquintic.h` describes "6x6x6 coefficients ... voxels at
   (x,y,z)" and `HermiteTriquintic.h` "Hermite cubic interpolation ... 4x4x4
   coefficients"; the evaluation comments name classes that do not exist
   (`IntpHermiteCubic2<T>`, `IntpHermiteCubic3<T>`, `HermiteCubic{2,3}`).
   Severity: doc. The formulas themselves are right: the corner conditions
   above hold to round-off for all six classes.
3. **`IntpBicubic2.h` asserts a 3x3 minimum** with the comment "needed to
   construct the estimates of the boundary derivatives", copied from the Akima
   headers; the bicubic interpolator estimates no derivatives, and its clamped
   4x4 stencil works on a 2x2 grid. This is the bicubic twin of #69's
   `IntpBilinear2` comment item. Severity: doc; port: preserved (asserts 3x3).
