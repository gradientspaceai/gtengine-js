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
  `pending | ported | omitted`. It and the generated `src/index.ts` are
  maintained by the orchestrator on `main` (`npm run update:status`,
  `npm run gen:index`) — **port PRs must not touch either file**, so parallel
  PRs never conflict. Record deviations from upstream in the PR description.

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

## Established precedents (follow these; set by the foundation port)

- **Global export uniqueness**: `src/index.ts` star-exports every file, so every
  exported symbol must be unique across the whole library. `npm run gen:index`
  fails on duplicates. Strategies, in order of preference: (1) reuse via
  subclassing like upstream (`GVector extends Vector`, `Quaternion` over
  `Vector`-of-4) so base free functions (`dot`, `normalize`, …) apply; (2) add
  overloads inside the file that owns the name only if you are porting that
  file; (3) prefix with the type context (`IEEEClassification`, not
  `Classification`).
- **Cont\* naming**: upstream overloads `GetContainer`/`InContainer`/
  `MergeContainers` per bounding volume; the port suffixes each with the
  container type - `getContainerAlignedBox`, `inContainerCircle2`,
  `mergeContainersSphere3`, ... (set by B38). Vestigial always-`true` bool
  returns are dropped in favor of returning the container.
- **Dist\* shape** (set by B41): one class per file, named after the file,
  `implements DCPQuery<Type0, Type1, Result>` with `compute()`; nested Result
  structs become exported `<ClassName>Result` interfaces; private static
  helpers become module-private functions; `friend`-granted protected helpers
  export as `<fileNameCamel>DoQuery`; dimension aliases (`DCPPointLine`, ...)
  are dropped since runtime-dimension `Vector` serves every N.
- **Query bases are interfaces**: `DCPQuery<Type0, Type1, Result>` with
  `compute()`, `TIQuery` with `test()`, `FIQuery` with `find()`. Upstream Intr*
  headers containing both TI and FI specializations become two classes with
  suffixes: `IntrRay3Sphere3TI` / `IntrRay3Sphere3FI`, result types
  `...TIResult` / `...FIResult`. Dist* headers export one class named after the
  file. When an upstream Intr* class has multiple `operator()` overloads, only
  the canonical two-argument query keeps `test`/`find`; the others get
  descriptive method names (e.g. `findDynamic`, `testFiniteSemiInfinite` — see
  `src/IntrIntervals.ts`). Result-kind integer constants become a
  file-qualified exported enum (e.g. `IntrIntervalsFIResultType`).
- **Vector conventions**: runtime dimension over `number[]`; element access
  `get(i)`/`set(i, x)` (public `values` for hot loops); operators are module
  functions `negate/add/sub/mul/div/compMul/compDiv` returning new vectors;
  `normalize`/`orthonormalize` mutate in place and return the length (as
  upstream); comparisons are methods (`equals`, `lessThan`, … lexicographic);
  dimension mismatches throw via `logAssert`; `HLift/HProject/Lift/Project` →
  `hlift/hproject/lift/project`.
- **Errors**: use `logAssert(cond, msg)` / `logError(msg)` from `src/Logger.ts`
  (they throw `Error`; no file/line prefix — JS stacks carry that).
- **Arbitrary precision marker**: `BSNumber`/`BSRational`/`QFNumber` implement
  `ArbitraryPrecisionNumber` from `src/TypeTraits.ts`
  (`isArbitraryPrecision: true`, `hasDivisionOperator`).
- **Multi-dim arrays**: `Array2/3/4` use flat storage, `a.get(i0, i1)` with
  index order matching constructor bound order (C++ `a[i1][i0]` → `a.get(i0, i1)`).
- **Constructors**: default constructors zero-fill (upstream leaves
  uninitialized); ambiguous C++ constructor overloads become static factories
  (`fromEncoding/fromNumber/…`).

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

## Upstream bugs

If while porting you find a suspected bug in the upstream C++ (wrong math,
undefined behavior, unreachable code, stale documentation), do all of:

1. Port the fix if the bug would corrupt results in TS (document the fix in a
   code comment referencing upstream), or preserve the quirk faithfully if it
   is harmless — judgment call, explained in the PR.
2. Describe it in the PR body under a heading exactly titled
   `## Upstream bug suspects`, with file, line context, evidence, and whether
   the port fixes or preserves it.

The orchestrator files each confirmed suspect as a GitHub issue labeled
`upstream-bug` — the canonical list of upstream findings.

## Quality gates (every PR)

1. `npm run typecheck` passes (strict).
2. `npm test` passes; each ported file has a test file with meaningful cases:
   known values, mathematical identities, degenerate inputs, and where
   practical randomized cross-checks against an independent computation.
3. No new runtime dependencies. devDependencies only, and only with good reason.
4. PRs touch only `src/<PortedName>.ts` and `test/<PortedName>.test.ts` files —
   never `src/index.ts`, `porting-status.json`, or config files.
