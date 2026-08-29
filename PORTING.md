# Porting conventions: GTE/Mathematics (C++) → TypeScript

Upstream baseline: `davideberly/GeometricTools` commit
`d29e7758ae2615e5e37da3eb573b7bf90ee94e9b`, folder `GTE/Mathematics` only.

## File mapping

- One `src/<Name>.ts` per upstream `<Name>.h`, same basename
  (`Vector3.h` → `src/Vector3.ts`). Flat directory, mirroring upstream.
- Every source file is re-exported from `src/index.ts` (alphabetical order).
- Tests live in `test/<Name>.test.ts` (vitest).
- Every ported file starts with the upstream banner adapted to the port:

```ts
// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) <Name>.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
```

- `porting-status.json` maps every upstream header to
  `pending | ported | omitted`, with a `notes` field for omissions and
  deliberate deviations. Update it in the same commit as the port.

## Type mapping

| C++ | TypeScript |
| --- | --- |
| `template <typename Real/T>` (numeric) | `number` (IEEE f64 ≈ C++ `double`) |
| `float`-specific paths | `Math.fround` only where bit-exactness matters |
| `std::array<T, N>`, `std::vector<T>` | `T[]` (use `Float64Array`/`Int32Array` for large hot buffers) |
| `std::map` / `std::set` | `Map` / `Set` (mind ordering: use sorted arrays where iteration order matters) |
| `int32_t`, `size_t` indices | `number` (keep `| 0` / `Math.trunc` semantics where the C++ relies on integer division) |
| `uint32_t`/`uint64_t` bit manipulation | `>>> 0` idioms or `bigint` (e.g. `IEEEBinary`, `BitHacks`) |
| arbitrary precision (`BSNumber`, `BSRational`, `UInteger*`) | ported on top of `bigint` |
| `std::pair`/output reference params | returned object literals with named fields |

## API mapping

- Class/struct names keep upstream names exactly (`DCPQuery`, `IntrRay3Sphere3`, …).
- Methods and free functions become `camelCase` (`ComputeRoots` → `computeRoots`).
- `operator()` on query objects: `DCPQuery::operator()` → `compute(...)`,
  `TIQuery::operator()` → `test(...)`, `FIQuery::operator()` → `find(...)`.
- Arithmetic operators become named functions/methods: `add`, `sub`, `negate`,
  `mul` (scalar or matrix), `div`, `equals`, plus the specialized upstream free
  functions (`Dot`, `Cross`, `Normalize`, …) as camelCase module functions
  (`dot`, `cross`, `normalize`).
- Template dimension parameters (`Vector<N, Real>`) → runtime dimension
  (`new Vector(n)`), with `Vector2/3/4`-style helpers matching upstream
  specialization headers.
- C++ template specializations selected by type traits: port only the
  floating-point instantiation unless the rational/exact path is required by a
  dependent file; note the decision in `porting-status.json`.
- `GTE_ASSERT` / `LogError` / exceptions → helpers in `src/Logger.ts` that
  `throw new Error(...)` with the same message text.

## Semantics rules

- Preserve upstream algorithms and numerical behavior; do not "improve" math
  while porting. Refactors are allowed only where C++ idioms don't exist in TS,
  and must be behavior-preserving.
- Preserve upstream comments that explain the mathematics (adapted, not
  verbatim-required); drop C++-mechanical comments.
- Watch for: integer division and modulo on negatives, `std::numeric_limits`
  (→ `Number.EPSILON`, `Number.MAX_VALUE`, etc. — check exact constant),
  in-place mutation vs copy (C++ value semantics: assignment copies; TS objects
  alias — clone explicitly where upstream copies).
- Deterministic iteration: where upstream relies on `std::map`/`std::set`
  ordering, replicate with explicit sorting.

## Omissions

C++-runtime-specific utilities are omitted (recorded in
`porting-status.json`): thread/atomic wrappers, smart-pointer comparators, and
similar files with no mathematical content. File I/O utilities are ported using
`ArrayBuffer`/`DataView` where they carry geometric value (e.g. STL files).

## Quality gates (every PR)

1. `npm run typecheck` passes (strict).
2. `npm test` passes; each ported file has a test file with meaningful cases:
   known values, mathematical identities, degenerate inputs, and where
   practical randomized cross-checks against an independent computation.
3. No new runtime dependencies. devDependencies only, and only with good reason.
4. `porting-status.json` and `src/index.ts` updated in the same PR.
