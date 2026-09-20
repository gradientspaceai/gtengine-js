# Group 2 (`v02-algebra`) — C++ oracle report

52 cases, 20 records each in `oracle/golden/v02-algebra.txt` (1040 records, 25 of
them C++ exceptions from the two throw-parity cases). 30 cases are declared
`{ exact: true }`, 18 carry the default `1e-12` tolerance because their upstream
path calls `sin`, `cos`, `acos`, `asin` or `atan2`, and 4 are `deviation` cases
demonstrating deliberate port fixes of upstream defects.

Deep run `npm run oracle:deep -- 2000 v02-algebra`: 52 cases, 104000 records,
2521 C++ exceptions (all matched by the port), **2 778 728 floating-point
outputs compared**. The 2 352 728 outputs of the exact cases are **100.0000 %
bit-identical**; the 426 000 outputs of the tolerance cases are 97.2758 %
bit-identical with a **worst scaled error of 3.35e-15**, three orders of
magnitude inside the tolerance. All 53 tests pass.

Upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, MSVC 194435215 x64
`/O2 /fp:precise` (no `/arch:AVX2`, so no FMA contraction). Neither
`GTE_USE_COL_MAJOR` nor `GTE_USE_VEC_MAT` is defined, so the row-major /
`GTE_USE_MAT_VEC` default that the port implements is the configuration
measured.

**Overload resolution.** The translation unit includes `Matrix2x2.h`,
`Matrix3x3.h` and `Matrix4x4.h` (directly and through `Transform.h` and
`Projection.h`), so every unqualified `Inverse`/`Adjoint`/`Determinant`/`Trace`
of a 2x2, 3x3 or 4x4 matrix resolves to those headers' closed forms and not to
the Gaussian-elimination templates of `Matrix.h` that group 1 measures. That
includes the calls inside `Transform.h` (`GetHInverse`, `Inverse`) and inside
`Hyperellipsoid::FromCoefficients`, which `Projection.h` reaches — and it is
measurable: before group 34's matching port fix landed (PR #508), the ellipse
centre of `Projection.perspectiveProject.basis`, the only output that uses that
inverse, differed by 1 to 2 ulps on 12 of the 20 committed records while every
axis and extent already agreed bit for bit. After the rebase both
`perspectiveProject` cases are bit-exact.

## Coverage

| header | cases | comparison | deep run |
| --- | --- | --- | --- |
| `Matrix2x2.h` | `inverse`, `doTransform`, `rotation` | exact, exact, tolerance (`cos`, `sin`, `atan2`) | pass |
| `Matrix3x3.h` | `inverse`, `doTransform` | exact | pass |
| `Matrix4x4.h` | `inverse`, `doTransform`, `makeObliqueProjection`, `makePerspectiveProjection`, `makeReflection` | exact | pass |
| `Quaternion.h` | `arithmetic`, `multiply`, `rotate.3d`, `rotate.4d` | exact | pass |
| `Quaternion.h` (slerp) | `slerp`, `slerpR`, `slerpRP`, `slerpRPH`, `slerpRP.invalidAngle` | tolerance (`acos`, `sin`) | pass |
| `Rotation.h` | `quaternionToMatrix`, `matrixToQuaternion`, `passthrough` | exact | pass |
| `Rotation.h` | `axisAngleToMatrix`, `matrixToAxisAngle`, `axisAngleToQuaternion`, `quaternionToAxisAngle`, `matrixToEulerAngles`, `eulerAnglesToMatrix`, `quaternionToEulerAngles`, `eulerAnglesToQuaternion`, `axisAngleToEulerAngles`, `eulerAnglesToAxisAngle` | tolerance | pass |
| `Rotation.h` | `matrixToAxisAngle.nearPi`, `eulerAngles.refactorization` | deviation | pass |
| `RotationEstimate.h` | `rotCEstimate`, `rotationEstimate` | exact | pass |
| `Transform.h` | `identity`, `rotationMatrix`, `rotationQuaternion`, `setMatrix`, `inverse`, `multiply`, `assertions` | exact | pass |
| `Transform.h` | `rotationAxisAngle`, `rotationEulerAngles` | tolerance | pass |
| `Transform.h` | `inverse.matrixChannel`, `getHInverse.lastRow` | deviation | pass |
| `ConvertCoordinates.h` | `convert` | exact | pass |
| `Projection.h` | `projectEllipse2`, `projectEllipsoid3`, `perspectiveProject.basis`, `perspectiveProject.plane` | exact | pass |

Entry points reached, by header:

- **`Matrix2x2.h`** — `MakeRotation`, `GetRotationAngle`, `Inverse` (both the
  flag and the no-flag overload), `Adjoint`, `Determinant`, `Trace`,
  `DoTransform(M,V)`, `DoTransform(A,B)`, `SetBasis`, `GetBasis`.
- **`Matrix3x3.h`** — `Inverse` (both overloads), `Adjoint`, `Determinant`,
  `Trace`, both `DoTransform` overloads, `SetBasis`, `GetBasis`.
- **`Matrix4x4.h`** — the same set plus `MakeObliqueProjection`,
  `MakePerspectiveProjection` (whose last column reads entries of the same row
  assigned earlier in the function, so the assignment order is part of the
  result) and `MakeReflection` (six entries are copies of earlier ones).
- **`Quaternion.h`** — unary `operator+`/`operator-`, `operator+`, `operator-`,
  `operator*(q,s)`, `operator*(s,q)`, `operator/(q,s)` including the zero
  divisor, the Hamilton product `operator*(q0,q1)` in both orders, `Dot`,
  `Length`, `Normalize` (including the zero quaternion), `Conjugate`,
  `Inverse` (including the zero quaternion), `Rotate` for `Vector<3>` and
  `Vector<4>` (affine and general 4-tuples), `Slerp`, `SlerpR`, `SlerpRP`,
  `SlerpRPH` (both arms of `twoT <= 1`), the six comparisons and `Zero`, `I`,
  `J`, `K`, `Identity`.
- **`Rotation.h`** — all four constructors, all four conversion operators and
  all twelve `Convert` overloads, for `N = 3` and `N = 4` on every record:
  matrix <-> quaternion, matrix <-> axis-angle, quaternion <-> axis-angle,
  matrix <-> Euler angles, quaternion <-> Euler angles, axis-angle <-> Euler
  angles, plus the four identity (pass-through) conversions. All twelve valid
  Euler axis orders and one invalid one, all four `EulerResult` codes.
- **`RotationEstimate.h`** — `RotC0Estimate` … `RotC4Estimate` and
  `GetRotC0EstimateMaxError` … `GetRotC4EstimateMaxError` for all seven
  degrees (4, 6, 8, 10, 12, 14, 16), `RotationEstimate`,
  `RotationDerivativeEstimate` and `RotationAndDerivativeEstimate`.
- **`Transform.h`** — the default constructor, `Identity()`, the implicit
  `operator Matrix4x4 const&`, `MakeIdentity`, `MakeUnitScale`, `IsIdentity`,
  `IsRSMatrix`, `IsUniformScale`, all five `SetRotation` overloads (4x4, 3x3,
  quaternion, `AxisAngle<3>`, `AxisAngle<4>`, Euler angles), `SetMatrix`, all
  three `SetTranslation` overloads, all three `SetScale` overloads,
  `SetUniformScale`, all five `GetRotation` forms, `GetMatrix`,
  `GetTranslation`, `GetTranslationW0`, `GetTranslationW1`, `GetScale`,
  `GetScaleW1`, `GetUniformScale`, `GetNorm` (both the RS and the
  max-row-sum arm), `GetHMatrix`, `GetHInverse` (uniform, non-uniform and
  general arms, including a singular homogeneous matrix), `Inverse` (both
  arms) and all five `operator*` overloads.
- **`ConvertCoordinates.h`** — `operator()` (both the invertible and the
  early-exit arm), `GetC`, `GetInverseC`, `IsVectorOnRightU`,
  `IsVectorOnRightV`, `IsRightHandedU`, `IsRightHandedV`, `UToV`/`VToU` for
  vectors and for matrices in all four `vectorOnRight` combinations, in
  dimensions 2, 3 and 4.
- **`Projection.h`** — both `Project` overloads and both `PerspectiveProject`
  overloads (the precomputed-basis one and the single-plane one, which derives
  `U`, `V` and the near distance with `ComputeOrthogonalComplement`).

### Generators and the branches they reach

Matrix- and vector-valued inputs mix four populations, one recorded double per
component in each, so the TypeScript replay reads a fixed layout: a small
integer lattice, uniform, a constructed exactly singular matrix (the last row
an integer combination of the others), and a "wild" population of signed
zeros, the smallest subnormal, subnormals on both sides of the `1/x` overflow
threshold and values near the overflow limit.

Rotation-valued inputs have their own modes, all recorded as the derived
values the query actually uses: the eight axis-aligned unit quaternions and
the sixteen `(±1/2, ±1/2, ±1/2, ±1/2)`; a normalized uniform 4-tuple; the same
with `w` forced negative; and near-identity or near-pi quaternions with
`2^-e`, `e` in `[1, 30]`. Rotation matrices come either from those quaternions
or, in mode 0, from the 24 signed permutation matrices of determinant `+1`,
whose entries are exactly `0` and `±1` — which is what makes the Euler
gimbal-lock tests `r(i,j) == ±1` fire at exact equality. All frames are
right-handed: rotation frames are the columns of a determinant `+1` matrix,
and the 2D ellipse frame uses `-Perp(axis0)`, since GTE's `Perp(x,y) = (y,-x)`
makes `{axis0, Perp(axis0)}` left-handed.

Distributions from the 2000-record deep run:

| case | branch histogram |
| --- | --- |
| `Matrix2x2.inverse` | invertible 1333, singular 667 |
| `Matrix3x3.inverse` | invertible 1360, singular 640 |
| `Matrix4x4.inverse` | invertible 1405, singular 595 |
| `Rotation.matrixToQuaternion` | all four arms of the trace split: x 652, y 422, z 389, w 537 |
| `Rotation.matrixToAxisAngle` | angle 0 (identity arm) 49, angle in (0,pi) 1731, angle pi (symmetric arm) 220 |
| `Rotation.matrixToEulerAngles` | UNIQUE 1674, NOT_UNIQUE_SUM 93, NOT_UNIQUE_DIF 71, INVALID 162 |
| `Rotation.quaternionToEulerAngles` | UNIQUE 1657, NOT_UNIQUE_SUM 85, NOT_UNIQUE_DIF 87, INVALID 171 |
| `Quaternion.slerp` | `acos`/`sin` arm 1554, the exact `cosA >= 1` arm 446 |
| `Quaternion.slerpRP.invalidAngle` | 1000 throws (angle pi), 1000 normal returns |
| `ConvertCoordinates.convert` | invertible 1387, singular 613; n = 2/3/4: 633/680/687; the four `vectorOnRight` combinations 481/550/454/515 |
| `Transform.inverse` | the five structures (identity, RS unit scale, RS uniform, RS non-uniform, general) 406/397/379/387/431 |
| `Transform.multiply` | `kindA` 394/416/376/418/375, `kindB` 387/416/357/377/442, 21 throws |
| `Transform.rotationMatrix` | uniform scale 1500, non-uniform 500 |
| `Transform.assertions` | 1500 throws (zero scale, non-RS `SetScale`, non-uniform `GetUniformScale`), 500 normal returns |
| `RotationEstimate.rotCEstimate` | all seven degrees, 248 to 320 records each |
| `Rotation.matrixToAxisAngle.nearPi` | upstream returns a zero axis 1363, a non-unit axis 446, a sound axis 191 |

The committed 20-record goldens reach the same branches: both invertibility
arms of every `Inverse`, all four arms of the matrix-to-quaternion split, all
four `EulerResult` codes, both arms of `SlerpRPH`, both `PerspectiveProject`
overloads, the throw arms of both throw-parity cases and 19 of the 20
`nearPi` records deviating.

### Where a tolerance is used, and why

libm is kept out of the compared path wherever possible: unit quaternions,
rotation matrices, the `SlerpRP` cosine and the `SlerpRPH` midpoint and
half-angle cosine are all *recorded as inputs*, so both sides start from the
same doubles. That is why 30 of the 52 cases are exact. The remaining 18 call
libm on the compared path itself:

- `Matrix2x2.rotation` — `cos`, `sin` (`MakeRotation`) and `atan2`
  (`GetRotationAngle`). Worst deep-run scaled error 2.22e-16, one ulp.
- `Quaternion.slerp`, `slerpR`, `slerpRP`, `slerpRPH`,
  `slerpRP.invalidAngle` — `ChebyshevRatiosUsingCosAngle` calls `acos` and
  `sin`. Worst 3.33e-16.
- The ten libm `Rotation` conversions and the two libm `Transform` rotation
  channels — `cos`, `sin`, `acos`, `asin`, `atan2`. Worst 3.35e-15
  (`Rotation.axisAngleToEulerAngles`: a quaternion built with `sin`/`cos`, a
  matrix built from it and then three `atan2`/`asin` calls, so about four
  roundings accumulate on a quantity of magnitude 1).

Every one of these errors is at the size of a one-ulp libm difference times a
small growth factor; none approaches the 1e-12 tolerance. The proof that they
are libm differences and not formula differences is structural: where the same
conversion is driven from a *recorded* matrix or quaternion instead of from
angles, the very same code is bit-identical. `Rotation.matrixToQuaternion`,
`Rotation.quaternionToMatrix` and `Transform.rotationQuaternion` are exact;
`Rotation.matrixToAxisAngle` and `Rotation.quaternionToAxisAngle` emit the
axis — the part computed with `sqrt` alone — through the harness's new
`outRealExact`, so it is required to be bit-identical while only the `acos`
angle carries the tolerance, and both pass at 100 % on that axis over 2000
records.

**libm-decided control flow** is never compared. `Convert(Matrix, Quaternion)`
branches on `r22 <= 0` and then on `r11 ∓ r00 <= 0`, and
`Convert(Matrix, EulerAngles)` on `r(i,j) < 1 / > -1`; a 1-ulp difference in
`cos` could move either across its threshold, which no tolerance repairs. The
cases whose matrix comes from *recorded* input feed both sides identical
doubles and therefore take identical branches — that is where the branch
coverage in the table above comes from. The composite cases that build the
matrix from angles (`eulerAnglesToQuaternion`, `eulerAnglesToAxisAngle`,
`axisAngleToEulerAngles`, `Transform.rotationEulerAngles`) reject candidates
whose branch quantities are within `1e-6` of the threshold, ten orders of
magnitude above a 1-ulp move; every rejection loop is capped at 32 attempts and
redraws every quantity its acceptance test depends on.
`Rotation.eulerAnglesToAxisAngle` additionally rejects `|q[3]| > 1 - 1e-3`,
which is a conditioning bound rather than a branch: the angle is `2*acos(q[3])`
and a 1-ulp difference in the intermediate `q[3]` moves it by
`1.1e-16/sqrt(2(1-|q3|))`, which is 5.5e-14 relative at the threshold and
unbounded as `|q3|` approaches 1. `Convert(Matrix, AxisAngle)`'s own branches
are *not* decided by libm despite calling `acos`: `acos(cs) > 0` holds exactly
when `cs < 1` and `acos(cs) < pi` exactly when `cs > -1` on both runtimes, so
the arms follow from exact comparisons on `cs`.

### Independent cross-check of the results

Agreement between the port and upstream would not reveal an error they share,
so the committed goldens' *outputs* were also checked against invariants, in
plain JavaScript with no gtengine-js import (worst residual over the 20
records of each case): `Rotation.quaternionToMatrix` produces matrices with
`R^T R = I` and `det R = 1` (8.9e-16); `Rotation.matrixToQuaternion` returns a
unit quaternion that reproduces the input matrix (4.4e-16);
`Rotation.matrixToAxisAngle` returns a unit axis whose axis-angle matrix
reproduces the input (4.2e-13); `Quaternion.rotate.3d` equals `R(q)*u`
(1.8e-15); `Matrix3x3.inverse` satisfies `M * M^{-1} = I` whenever the
invertibility flag is set (1.1e-14); and `Transform.rotationMatrix` satisfies
`H * H^{-1} = I` (5.3e-15). The last three are restricted to well-scaled
records, because the deliberate "wild" population destroys such identities by
conditioning rather than by a wrong result.

## Port defects fixed

**`src/Vector.ts` / `src/Quaternion.ts` — `normalize` used `Vector.h`'s scalar
division for quaternions.**

Upstream keeps two copies of the scalar `operator/=`. `Vector.h`'s forms the
reciprocal once and multiplies,

```cpp
Real invScalar = (Real)1 / scalar;
for (int32_t i = 0; i < N; ++i) { v[i] *= invScalar; }
```

while `Quaternion.h`'s divides each component,

```cpp
for (int32_t i = 0; i < 4; ++i) { q[i] /= scalar; }
```

Both headers write `Normalize` as `v /= length`, so the two `Normalize`s are
different computations: `a * (1/b)` and `a / b` differ in the last bit for most
`b`, since the reciprocal is rounded once before the multiplication. The port
merged `Quaternion` into `Vector` (`Quaternion extends Vector` and reuses its
free functions) and documented `normalize(q)` as "exactly upstream's
Normalize", which was true of `Dot` and `Length` but not of this. Size of the
error: one ulp of a unit-length component, on 4 of the 20 committed records of
`Quaternion.arithmetic` (and on every record where the reciprocal's rounding
is not absorbed). It propagates through every user of a normalized quaternion.

Fixed by giving `Vector` a `divideByScalarUsesReciprocal` getter, `true` on
`Vector` and overridden to `false` on `Quaternion`, and letting the non-robust
`normalize()` switch on it — the same shape as the `dotAccumulatesFromZero`
hook that group 1 added for `GVector`'s `Dot` seed. Only the non-robust arm
consults it: `Quaternion.h` has no robust `Normalize`, so the robust arm stays
`Vector.h`'s. `Vector.h`'s own behaviour is unchanged, which the exact
`v01-algebra` cases continue to confirm.

Regression test: `test/Quaternion.test.ts`, "Normalize division (Vector.h vs
Quaternion.h)". It runs both candidate computations side by side, asserts that
they really do differ on at least one input, asserts that the port matches
`Quaternion.h` bit for bit for a `Quaternion` and `Vector.h` bit for bit for a
plain `Vector`, checks that the flag follows the type, and property-checks
bit-identity with `Quaternion.h`'s form over random quaternions.

After the fix the full suite is green (564 files, 11328 tests) and no other
family's golden replay changed.

## Deliberate deviations demonstrated

All four are recorded in `docs/UPSTREAM-FINDINGS.md`; each `deviation` case is
confined to inputs on which upstream is actually defective, and the
corresponding main case rejects those inputs by an explicit predicate or avoids
the state that triggers them.

- **`Rotation.matrixToAxisAngle.nearPi`** (issue
  [#374](https://github.com/gradientspaceai/gtengine-js/issues/374)) — for an
  angle just below pi, upstream extracts the axis from `R - Transpose(R)`,
  whose entries are `2 sin(angle) a_i`. The case builds a coordinate-axis
  rotation by pi (exact entries `0` and `±1`) whose trace is lifted by `2^-50`,
  so `acos` reports an angle just below pi and the antisymmetric arm runs, plus
  an antisymmetric perturbation of about `2^-k` with `k` in `[513, 600]`. For
  `k` in `[513, 538]` the squares of the extracted components are subnormal and
  lose most of their significand, so `Normalize` returns a non-unit axis (446
  of 2000 deep-run records); for `k >= 539` they underflow to zero and it
  returns the zero vector (1363 records). Either way upstream's "axis" is not a
  rotation axis. The port validates the normalized axis and falls back to the
  symmetric formula. 19 of the 20 committed records deviate. The main case
  `Rotation.matrixToAxisAngle` rejects exactly the candidates this predicate
  rejects — which over 2000 random rotation matrices is none, since the regime
  is not reachable by sampling — and is bit-identical on every output, so the
  fix is confined to the defect.
- **`Rotation.eulerAngles.refactorization`** (issue
  [#225](https://github.com/gradientspaceai/gtengine-js/issues/225)) —
  `operator()(i0,i1,i2)` writes the requested axis indices into the cached
  `EulerAngles` and then does nothing at all when the `Rotation` was built from
  Euler angles, so asking for a *different* factorization returns the original
  angles relabelled with the new axes, which is a different rotation. The port
  recomputes through the rotation matrix, as every other source type does. The
  case always requests a different order; 20 of 20 committed records deviate.
  `Rotation.passthrough` covers the arm where the requested order matches the
  stored one, which is bit-identical to upstream.
- **`Transform.inverse.matrixChannel`** (issue
  [#265](https://github.com/gradientspaceai/gtengine-js/issues/265)) — the
  general arm of `Transform::Inverse` passes the full affine 4x4 to
  `SetMatrix`, so the inverse's M channel carries the translation in its last
  column and whatever the closed-form `Inverse` produced in its last row,
  violating `GetMatrix`'s documented `{{M,0},{0,1}}` structure. The port zeroes
  them. 20 of 20 records deviate. The main case `Transform.inverse` emits every
  other channel of the inverse — flags, `H`, `H^{-1}`, translation, norm,
  rotation, scale — on all five transform structures and is bit-identical, so
  the fix touches only the accessor it is about.
- **`Transform.getHInverse.lastRow`** (issue #265) — the RS arms of
  `GetHInverse` write only the upper-left 3x3 block and the last column of
  `mInvHMatrix`, assuming the last row still holds `(0,0,0,1)` from the
  constructor. After the general arm has assigned `Inverse(mHMatrix)` for a
  singular `mHMatrix` — the zero matrix — that assumption is false, so a later
  RS state reports an inverse whose last row is `(0,0,0,0)`. The case forces
  exactly that sequence (a general matrix with a singular 3x3 block, a
  `GetHInverse` that caches the zero matrix, then a rotation and a second
  `GetHInverse`) and 20 of 20 records deviate, in the single entry `(3,3)`.

The preserved defects of these headers are ordinary bit-for-bit cases:
`Convert(Matrix, AxisAngle)`'s absorbed round-off away from the underflow
regime (an axis wrong by up to 5e-3, which the port keeps because fixing it
would need an input threshold) and `Convert(Matrix, EulerAngles)`'s exact
`±1` gimbal-lock test are both reproduced exactly, the latter by the signed
permutation matrices that evaluate it at exact equality.

## Not covered

| header / entry point | reason |
| --- | --- |
| `Quaternion::Log` / `Exp` | not present in this upstream version of `Quaternion.h`; the header has no logarithm or exponential. |
| `Transform::Invert3x3` | private, dead code — nothing in the header calls it. Recorded as a preserved minor finding under issue #265; it also assigns `0.0f`/`1.0f` float literals into `Real` entries, which no `double` instantiation can observe. |
| `Matrix2x2.h`/`3x3`/`4x4` `DoTransform`, `SetBasis`, `GetBasis` under `GTE_USE_VEC_MAT` | the port implements only `GTE_USE_MAT_VEC`, the GTE default, and the oracle is built in that configuration. The `#if` branches are the other configuration, not another code path of this one. |
| `Rotation`'s `static_assert(N == 3 \|\| N == 4)` | a compile-time check; the port's runtime `logAssert` counterpart has no C++ record to compare against. |
| `ConvertCoordinates::UToV`/`VToU` for `N = 1` | upstream's `GaussianElimination` and the port both handle it, but a 1x1 change of basis has no branch the larger dimensions do not; dimensions 2, 3 and 4 are generated. |
| `Projection::PerspectiveProject` for a violated precondition | the documented precondition is that the ellipsoid lies strictly between the eyepoint and the view plane; upstream then returns a conic that is not an ellipse and `FromCoefficients` returns `false`, leaving the ellipse members *undefined*. Upstream discards that flag, so there is nothing well-defined to compare. The generator satisfies the precondition by construction and the TypeScript replay throws if the port's flag is ever `false` (it never is). |
| `AxisAngle.h`, `EulerAngles.h` | pure data, covered by group 1's report. |

## Upstream bug suspects

None new. Everything this group's headers carry is already recorded in
`docs/UPSTREAM-FINDINGS.md` (issues #160, #225, #265, #374). The one new defect
the oracle found is on the port's side and is described under "Port defects
fixed".

One observation worth recording, and the exact analogue of group 1's note about
`Vector.h`'s and `GVector.h`'s `Dot`: `Vector.h`'s and `Quaternion.h`'s scalar
`operator/=` are not interchangeable in floating point, so `Normalize` means two
different computations in the two headers even though the source line is
identical. A library that presents `Quaternion` as a 4-tuple with the usual
vector operations would be better off with one implementation.
