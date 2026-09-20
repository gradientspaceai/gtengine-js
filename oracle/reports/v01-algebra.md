# Group 1 (`v01-algebra`) — C++ oracle report

46 cases, 20 records each in `oracle/golden/v01-algebra.txt` (920 records, 18 of
them C++ exceptions from the one throw-parity case). The 902 ordinary records
carry 15792 floating-point outputs and **every one of them is bit-identical** to
the MSVC build of upstream GTE. Every case is declared `{ exact: true }`: this
group is `+ - * /`, `sqrt`, `fabs` and comparisons only, and `sqrt` is correctly
rounded in both the MSVC runtime and V8. There are no tolerances and no
`deviation` cases.

Deep run `npm run oracle:deep -- 2000 v01-algebra`: 46 cases, 92000 records,
1675 C++ exceptions (all in `GMatrix.access.invalidIndex`, matched by the port),
**1602881 floating-point outputs compared, 100.0000 % bit-identical**, worst
scaled error 0. All 46 cases are bit-identical on every output.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise` (no `/arch:AVX2`, so no FMA contraction). Neither
`GTE_USE_COL_MAJOR` nor `GTE_USE_VEC_MAT` is defined, so the row-major /
`MAT_VEC` default that the port implements is the configuration measured. Only
`Vector.h`, `Vector2.h`, `Vector3.h`, `Vector4.h`, `GVector.h`, `Matrix.h` and
`GMatrix.h` are included, so `Inverse`, `Determinant` and `Adjoint` resolve to
the Gaussian-elimination templates of `Matrix.h`/`GMatrix.h` rather than to the
closed forms of `Matrix2x2.h`/`Matrix3x3.h`/`Matrix4x4.h`, which belong to
group 2. That matters here because the port's `inverse()`/`determinant()` always
take the Gaussian path, so this is the comparison that is meaningful.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `Vector.h` | `arithmetic`, `dot`, `dot.signedZero`, `length`, `length.subnormal`, `normalize`, `orthonormalize`, `getOrthogonal`, `computeExtremes`, `liftProject`, `special` | exact | pass |
| `Vector2.h` | `perp`, `computeOrthogonalComplement`, `computeBarycentrics`, `intrinsicsVector2` | exact | pass |
| `Vector3.h` | `cross`, `cross.4d`, `computeOrthogonalComplement`, `fastComputeOrthogonalComplement`, `computeBarycentrics`, `intrinsicsVector3` | exact | pass |
| `Vector4.h` | `hyperCross`, `computeOrthogonalComplement`, `computeOrthogonalComplement.zeroMiddle` | exact | pass |
| `GVector.h` | `arithmetic`, `dot`, `dot.signedZero`, `lengthNormalize`, `orthonormalize`, `computeExtremes`, `liftProject`, `special` | exact | pass |
| `Matrix.h` | `arithmetic`, `multiply`, `multiplyDiagonal`, `inverse`, `inverse.denormal`, `access`, `hliftProject` | exact | pass |
| `GMatrix.h` | `arithmetic`, `multiply`, `multiplyDiagonal`, `inverse`, `access`, `access.mismatchedSizes`, `access.invalidIndex` | exact | pass |

Entry points reached, by header:

- **`Vector.h`** — unary `operator+`/`operator-`, `operator+`, `operator-`,
  `operator*(v, s)`, `operator*(s, v)`, `operator/(v, s)` (including the zero
  divisor, which upstream special-cases to the zero vector), componentwise
  `operator*(v0, v1)` and `operator/(v0, v1)`, `Dot`, `Length` (plain and
  robust), `Normalize` (plain and robust), `Orthonormalize`, `GetOrthogonal`
  (`N` = 2, 3, 4, both `unitLength` values), `ComputeExtremes` (including the
  invalid-input arm), `HLift`, `HProject`, `Lift`, `Project`, `MakeZero`,
  `MakeOnes`, `MakeUnit`, `Zero`, `Ones`, `Unit`, `GetSize` and the six
  comparison operators.
- **`Vector2.h`** — `Perp`, `UnitPerp` (both `robust` values), `DotPerp`,
  `ComputeOrthogonalComplement` (both the `numInputs == 1` and the
  invalid-input arm), `ComputeBarycentrics` (both arms, three epsilon regimes)
  and the whole `IntrinsicsVector2` constructor.
- **`Vector3.h`** — `Cross`, `UnitCross`, `DotCross` for `N` = 3 and for `N` = 4
  (affine, `w = 0`, and general 4-tuples), `ComputeOrthogonalComplement`
  (`numInputs` 0, 1 and 2, both arms of the `|v0[0]| > |v0[1]|` test),
  `FastComputeOrthogonalComplement` (both arms of `v2[2] >= 0`, with exact
  coordinate axes so the test is evaluated at `0`, `+1` and `-1`),
  `ComputeBarycentrics` and the whole `IntrinsicsVector3` constructor.
- **`Vector4.h`** — `HyperCross`, `UnitHyperCross`, `DotHyperCross`,
  `ComputeOrthogonalComplement` for `numInputs` 0, 1, 2 and 3, reaching all
  three `maxIndex` arms of the `numInputs == 1` step and all three of the
  `numInputs == 2` step.
- **`GVector.h`** — the same set on the resizable tuple, plus `SetSize`,
  `GetSize` and the `std::vector` comparison operators on tuples of equal and
  of different sizes.
- **`Matrix.h`** — unary `operator+`/`operator-`, `operator+`, `operator-`,
  `operator*`/`operator/` by a scalar (including the zero divisor), `L1Norm`,
  `L2Norm`, `LInfinityNorm`, `Inverse` (both overloads), `Determinant`,
  `Transpose`, `M*V`, `V^T*M`, `A*B`, `MultiplyAB`, `MultiplyABT`,
  `MultiplyATB`, `MultiplyATBT`, `MultiplyMD`, `MultiplyDM`, `OuterProduct`,
  `MakeDiagonal`, `HLift`, `HProject`, `SetRow`, `SetCol`, `GetRow`, `GetCol`,
  `operator()`, the 1D `operator[]`, `MakeZero`, `MakeUnit` (valid and invalid
  indices), `MakeIdentity`, `Zero`, `Unit`, `Identity` and the six comparison
  operators.
- **`GMatrix.h`** — the same set on the resizable table, plus `SetSize`,
  `GetSize`, `GetNumRows`, `GetNumCols`, `GetNumElements`, the range-checked
  `operator()` (throw parity) and the dimension-aware comparisons.

### Generators and the branches they reach

Every generator mixes three populations, one recorded double per component in
each, so the TypeScript replay reads a fixed layout and never needs to know
which mode produced a record:

- **small integer lattice** (`[-3, 3]` or `[-4, 4]`) — exact arithmetic, exact
  zeros, ties, parallel and linearly dependent configurations, exactly singular
  matrices;
- **uniform** `[-10, 10]`;
- **wild** — signed zeros, the smallest subnormal `2^-1074`, subnormals on both
  sides of the `1/x` overflow threshold (`2^-1060` overflows, `2^-1023` does
  not), the smallest normal, values near overflow (`2^1020`) and integer-valued
  doubles at `2^52`.

Several cases add a fourth, constructed mode: exactly linearly dependent vector
sets, collinear and coplanar lattice point sets, exactly singular matrices (last
row an integer combination of the others), and nearly singular matrices (the
same with a single added power of two between `2^-1070` and `2^-40`, so full
pivoting ends on a minute pivot).

Distributions from the 2000-record deep run:

| case | branch histogram |
| --- | --- |
| `Matrix.inverse` | `invertible = false` 500, `invertible = true` 1441, `invertible = true` with NaN inverse entries 59 |
| `GMatrix.inverse` | `invertible = false` 511, `true` 1448, `true` with NaN entries 41 |
| `Matrix.inverse.denormal` | all-denormal matrices: `invertible = false` 1064, `true` with a finite inverse 100, `true` with NaN inverse entries 836, of which 232 also have a NaN determinant |
| `Vector2.intrinsicsVector2` | `dimension` 0: 887, 1: 262, 2: 851 |
| `Vector3.intrinsicsVector3` | `dimension` 0: 878, 1: 312, 2: 371, 3: 439 |
| `Vector2.computeBarycentrics` | linearly dependent 447, independent 1553 |
| `Vector3.computeBarycentrics` | linearly dependent 503, independent 1497 |
| `Vector.orthonormalize` | `minLength` exactly 0 in 269 records |
| `Vector2/3/4.computeOrthogonalComplement` | `minLength` exactly 0 in 1059 / 909 / 669 records |
| `Vector.length.subnormal` | robust `Length` NaN in 484 records (issue #370), finite in 1516 |
| `Vector.length`, `GVector.lengthNormalize` | robust `Length` NaN in 3 records each, from the wild mode alone |
| `Vector.arithmetic` | zero divisor in 169 records |
| `Vector.computeExtremes` | invalid-input arm in 332 records |
| `GMatrix.access.invalidIndex` | 1675 throw records, 325 in-range |
| `GMatrix.access.mismatchedSizes` | different dimensions in 1880 records |

In the committed 20-record goldens the same branches occur: 6 singular
`Matrix.inverse` records and 1 NaN-with-`invertible = true`, 8 of the 20
`Matrix.inverse.denormal` records NaN with `invertible = true` (2 of them with a
NaN determinant), 6 of the 20 `Vector.length.subnormal` records with a NaN
robust length, `dimension` 0/1/2 and 0/1/2/3 for the two `Intrinsics` cases, and
18 throw records for `GMatrix.access.invalidIndex`.

## Port defects fixed

**`src/Vector.ts` — `dot` used `Vector.h`'s accumulation seed for `GVector`.**

Upstream keeps two copies of the geometric free functions. `Vector.h`'s `Dot`
seeds the accumulation with the first product,

```cpp
Real dot = v0[0] * v1[0];
for (int32_t i = 1; i < N; ++i) { dot += v0[i] * v1[i]; }
```

while `GVector.h`'s seeds it with the literal zero,

```cpp
Real dot(0);
for (int32_t i = 0; i < v0.GetSize(); ++i) { dot += v0[i] * v1[i]; }
```

The port merged `GVector` into `Vector` (`GVector extends Vector` and reuses its
free functions), so `dot()` was `Vector.h`'s for both. The two seeds are equal in
exact arithmetic but not in floating point: `0 + x` is `x` for every `x` except
`-0`, where it is `+0`. They therefore disagree, in the sign of a zero, exactly
when every product `v0[i] * v1[i]` is `-0` — for example
`Dot((0, 0, 0), (-1, -2, -3))`, which upstream `GVector.h` reports as `+0` and
the port reported as `-0`. Size of the error: none in magnitude, one sign bit;
it propagates, because `x - (-0)` and `x - (+0)` differ when `x` is `-0`
(reached through `Orthonormalize`'s `v[i] -= v[j] * dot`) and because a later
division by the result differs in the sign of its infinity.

Found by `GVector.dot.signedZero` (6 of 20 committed records disagreed) and by
one record of `GVector.dot`. Fixed by giving `Vector` a `dotAccumulatesFromZero`
getter, `false` on `Vector` and overridden to `true` on `GVector`, and letting
`dot()` switch on it; `length`, `normalize` and `orthonormalize` inherit the
correct seed through `dot`. `Vector.h`'s own behaviour is unchanged, which the
`Vector.dot` / `Vector.dot.signedZero` cases confirm against the real build.

Regression test: `test/GVector.test.ts`, "Dot accumulation seed (Vector.h vs
GVector.h)". It runs both candidate groupings side by side, asserts that they
really do differ on at least one input (`-0` versus `+0`), asserts that the port
matches `GVector.h` bit for bit for `GVector` arguments and `Vector.h` bit for
bit for `Vector` arguments, and property-checks that the two agree everywhere
else.

After the fix the full suite is green (561 files, 11140 tests) and no other
family's golden replay changed.

## Deliberate deviations demonstrated

None. Both upstream defects that this group's headers carry are *preserved* by
the port, so they are ordinary bit-for-bit cases rather than `deviation` cases:

- **Issue [#370](https://github.com/gradientspaceai/gtengine-js/issues/370)** —
  the robust `Length`/`Normalize` rescale with `v /= maxAbsComp`, implemented as
  a multiplication by `1/maxAbsComp`; for a subnormal maximum below about
  `2^-1024` that reciprocal overflows and the zero components become
  `0 * inf = NaN`. `Vector.length.subnormal` generates both regimes and the port
  reproduces upstream's NaN exactly (484 of 2000 deep-run records).
- **Issue [#375](https://github.com/gradientspaceai/gtengine-js/issues/375)** —
  `GaussianElimination` on an all-denormal matrix returns NaN inverse entries
  while reporting `invertible = true`. `Matrix.inverse.denormal` reproduces it,
  including the records where the determinant is also NaN, and the `invertible`
  flag agrees on every record.
- **Issue [#87](https://github.com/gradientspaceai/gtengine-js/issues/87)** —
  `Vector4::ComputeOrthogonalComplement`'s `maxIndex == 3` branch yields the
  zero vector whenever components 1 and 2 vanish.
  `Vector4.computeOrthogonalComplement.zeroMiddle` pins it: all 2000 deep-run
  records return `minLength = 0` on both sides.
- **Issue [#88](https://github.com/gradientspaceai/gtengine-js/issues/88)** —
  `GMatrix`'s comparisons are `false` for `<`, `<=`, `>` and `>=` when the
  dimensions differ while `!=` is `true`, which is not a strict weak ordering.
  `GMatrix.access.mismatchedSizes` compares that behaviour directly.

## Not covered

| header / entry point | reason |
| --- | --- |
| `AxisAngle.h` | pure data: two members (`axis`, `angle`) and two constructors, no computation. The `static_assert(N == 3 \|\| N == 4)` is a compile-time check with no runtime counterpart. |
| `EulerAngles.h` | pure data: the `EulerResult` enum plus three members (`axis`, `angle`, `result`) and two constructors, no computation. |
| `Vector.h` / `GVector.h` constructors from `std::array` and `std::initializer_list` | value-copy constructors with no arithmetic; the truncate-and-zero-fill behaviour of the initializer-list form is covered by `test/Vector.test.ts`. |
| `GVector::operator/` and `GMatrix::operator/` with a zero divisor | upstream calls `LogError("Division by zero.")` where `Vector.h`/`Matrix.h` yield the zero tuple. The port shares `Vector.h`'s division for both (documented in `src/GVector.ts` and `src/GMatrix.ts`), so throw parity does not hold and the generators use a nonzero divisor. The fixed-size cases do include zero divisors. |
| `GVector::HProject` / `GVector::Project` of a size-1 tuple | upstream returns an empty `GVector`; the shared `hproject`/`project` throw, as `Vector.h`'s `static_assert(N >= 2)` does. Generators use sizes `>= 2`. |
| `GVector::Dot` / `Length` of an empty tuple | upstream `Dot` returns 0 while `Length` reads `v[0]` out of bounds (undefined behaviour). Generators use sizes `>= 1`. |
| `GVector::Orthonormalize` / `ComputeExtremes` invalid-input arms | upstream calls `LogError`; the shared `orthonormalize` returns 0 and `computeExtremes` returns `null`, as `Vector.h`'s bool-returning versions do. Generators use valid input. The fixed-size `Vector.computeExtremes` case does cover the invalid arm. |
| `LInfinityNorm` with a NaN in element 0 | `Matrix.h` seeds the maximum with `fabs(M[0])` and `GMatrix.h` with 0; the port shares `GMatrix.h`'s (documented in `src/Matrix.ts`). The two differ only when element 0 is NaN, which no generator here produces, since every drawn input is finite. Left as a documented port deviation rather than fixed: unlike the `Dot` seed, it is unreachable from finite inputs. |
| `Matrix`/`GMatrix`/`Vector`/`GVector` comparison operators as `std::map` keys | container behaviour, not computation. |

## Upstream bug suspects

None new. Everything this group found is already recorded in
`docs/UPSTREAM-FINDINGS.md` (issues #87, #88, #370, #375). The one new defect
found is on the port's side and is described under "Port defects fixed".

One observation worth recording for a future upstream note, already implied by
the `Vector.h` / `GVector.h` duplication: the two headers' `Dot` implementations
are not interchangeable in floating point, and neither are `Matrix.h`'s and
`GMatrix.h`'s `L1Norm`, `L2Norm` and `LInfinityNorm`. `L1Norm` and `L2Norm`
happen to agree on every input (their seeds `fabs(M[0])` and `M[0]*M[0]` are
never `-0`, and `0 + x` is `x` for every non-negative `x`), but `Dot` and
`LInfinityNorm` do not. A library that documents `GMatrix` as the run-time-sized
counterpart of `Matrix` would be better off with one implementation.
