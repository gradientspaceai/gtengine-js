# API conventions

gtengine-js mirrors the structure of upstream `GTE/Mathematics` one file per
header, so David Eberly's [documentation](https://www.geometrictools.com/Documentation/Documentation.html)
and the upstream source comments apply directly. This page records how C++
idioms map onto the TypeScript surface so that the whole library reads the
same way. The rules below are stable; anything not covered here may still
change before 1.0.

## Modules and names

- One module per upstream header, same basename: `Vector3.h` becomes
  `src/Vector3.ts`, importable from the package root
  (`import { cross } from 'gtengine-js'`). Every exported symbol is unique
  across the library, so the root import never collides.
- Classes keep upstream names exactly (`DCPQuery`, `IntrRay3Sphere3`,
  `BSNumber`). Methods and free functions are `camelCase`
  (`ComputeRoots` → `computeRoots`). Enums are `PascalCase` and file-qualified
  (`IntrIntervalsFIResultType`).
- Where upstream overloads one name on a type (`GetContainer` for boxes,
  circles, spheres, ...), the port appends the type:
  `getContainerAlignedBox`, `inContainerCircle2`, `mergeContainersSphere3`.
- Fixed-size matrix helpers carry their size uniformly: `inverse3x3`,
  `determinant4x4`, `makeRotation2x2`.
- Upstream mesh-topology records keep their single-letter member names
  (`V`, `E`, `T`, `S` for vertex/edge/triangle/tetrahedron index arrays, and
  the `VAdjacent`/`EAdjacent`/`TAdjacent` sets) because that is how the
  Geometric Tools documentation refers to them. `RigidBodyContact` likewise
  keeps `A`, `B`, `P`, `N`.

## Numbers

- Upstream `template <typename T>` numeric code is ported for `T = double`
  as plain `number`. `float`-specific paths use `Math.fround` only where
  bit-exactness matters (`IEEEBinary`, `BitHacks`).
- Arbitrary precision (`BSNumber`, `BSRational`, `QFNumber`, `APInterval`)
  is built on `bigint`. The C++ `UInteger` storage backends are not ported.
  Algorithms that upstream offers in both floating-point and exact flavours
  expose the exact flavour only where a dependent file needs it; the manifest
  `porting-status.json` records each decision.
- Integer semantics are preserved where upstream relies on them (`| 0`,
  `Math.trunc`, `>>> 0`).

## Vectors and matrices

- `Vector` and `Matrix` have a runtime dimension (`new Vector(3)`,
  `new Matrix(rows, cols)`, `Matrix.fromArray(...)`), so one class serves every upstream `N`.
  `Vector2`/`Vector3`/`Vector4` modules add the dimension-specific helpers
  (`perp`, `cross`, `dotCross`, barycentrics) as free functions on `Vector`.
- Element access is `v.get(i)` / `v.set(i, x)`; `v.values` (a `number[]`) is
  public for hot loops. `v.size` is the dimension.
- Arithmetic is free functions returning new objects: `add`, `sub`, `negate`,
  `mul` (scalar or matrix), `div`, `compMul`, `compDiv`, `dot`, `cross`,
  `length`, `normalize`. Matrix and quaternion variants that would collide are
  suffixed: `addMatrix`, `mulQuaternion`.
- `normalize` and `orthonormalize` mutate their argument in place and return
  the length, exactly as upstream. Everything else is non-mutating unless the
  upstream signature is an output parameter.
- `Quaternion` and `GVector` extend `Vector`, so the base free functions apply.

## Queries

Upstream expresses geometric queries as function objects; the port keeps the
class per query and names `operator()` by kind:

| upstream | port | result |
| --- | --- | --- |
| `DCPQuery<T, A, B>::operator()` | `new DistAB().compute(a, b)` | `DistABResult` (`distance`, `sqrDistance`, `closest`, parameters) |
| `TIQuery<T, A, B>::operator()` | `new IntrABTI().test(a, b)` | `IntrABTIResult` (`intersect`) |
| `FIQuery<T, A, B>::operator()` | `new IntrABFI().find(a, b)` | `IntrABFIResult` (`intersect`, `numIntersections`, `point`, ...) |

- An `Intr*` header that defines both a test and a find query becomes two
  classes with `TI` and `FI` suffixes. Extra `operator()` overloads get
  descriptive names (`findDynamic`, `testFiniteSemiInfinite`,
  `computeRobust`).
- `Result` types are plain interfaces returned by value; they never alias the
  query object's internal state, so a result may be kept across calls.
- Protected `DoQuery` helpers that upstream shares between line/ray/segment
  queries are exported module functions named `<file>DoQuery`
  (`intrLine3Sphere3DoQuery`) so the ray and segment queries can reuse them
  without inheritance.
- Query classes are stateless unless upstream has tunable state
  (`setMaxLCPIterations`, `useConjugateGradient`), and may be reused freely.

## Construction and results

- Default constructors produce a fully initialised zero object (upstream
  leaves members uninitialised). Ambiguous C++ constructor overloads become
  static factories: `AlignedBox.fromMinMax`, `Line.fromOriginDirection`,
  `Triangle.fromVertices`, `Hypersphere.fromCenterRadius`, `BSNumber.fromNumber`.
- C++ output parameters and `std::pair` returns become object literals with
  named fields (`{ inverse, invertible }`). Vestigial always-`true` boolean
  returns are dropped.
- Objects are mutable value types with `clone()`. Assignment aliases; call
  `clone()` where you would have relied on C++ copy semantics.

## Errors

- Upstream `LogError`/`LogAssert`/exceptions throw a plain `Error` carrying
  the upstream message text (see `src/Logger.ts`). There is no error subclass
  hierarchy and nothing is logged to the console.
- Dimension mismatches between vectors or matrices throw.

## Upstream fidelity

- Algorithms and numerical behaviour follow upstream commit
  `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`. Where the port deliberately
  departs from upstream to fix a defect, the source comment says so and the
  defect is documented as an
  [`upstream-bug` issue](https://github.com/gradientspaceai/gtengine-js/issues?q=label%3Aupstream-bug).
- `porting-status.json` lists every upstream header with its status
  (`ported` / `omitted`) and, once the verification wave has reviewed it, the
  `verified` group id.
