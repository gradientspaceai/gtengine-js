# gtengine-js

A TypeScript port of the **Mathematics** component of the
[Geometric Tools Engine (GTE)](https://www.geometrictools.com/) by David Eberly
([upstream repository](https://github.com/davideberly/GeometricTools)).

This is a machine-assisted port of the ~560 headers in
`GTE/Mathematics`: vector/matrix algebra, geometric primitives, distance
queries, intersection queries, containment, approximation/fitting,
interpolation, curves and surfaces, root finding, computational geometry,
and arbitrary-precision arithmetic.

Port baseline: upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b` (2026-08-13).

## Status

**The port is complete.** All 561 headers in upstream `GTE/Mathematics` are
accounted for: 550 are ported to `src/` and 11 are intentionally omitted (C++
runtime utilities such as thread/mutex wrappers, `std::atomic` helpers and the
32-bit `UInteger` storage backends, which are superseded by the `bigint`-backed
`BSNumber`/`BSRational` implementation). Every ported file has a Vitest test
file; the suite currently runs 550 test files and ~6700 tests.

See [docs/API.md](docs/API.md) for the public API conventions,
[PORTING.md](PORTING.md) for porting rules, [VERIFYING.md](VERIFYING.md) for the
verification wave in progress, and
[porting-status.json](porting-status.json) for a machine-readable per-file
manifest. The port was done in 135 dependency-ordered batches; the closed
[`port-batch` issues](https://github.com/gradientspaceai/gtengine-js/issues?q=label%3Aport-batch)
record what went into each one.

### Upstream defects found while porting

Porting each header with tests uncovered a number of bugs in the upstream C++.
Each is documented as an
[`upstream-bug` issue](https://github.com/gradientspaceai/gtengine-js/issues?q=label%3Aupstream-bug)
with the upstream file, line numbers, a reproducing input and how the port
handles it. Where the defect was clear-cut (wrong formula, dropped loop term,
stale state, uninitialised member) the port implements the corrected behaviour
and keeps a regression test; where the upstream behaviour is arguably
intentional or the fix is a design change, the port preserves upstream
semantics and the issue notes the caveat. These issues have not yet been
reported to the upstream project.

The API is not yet stable: names, parameter conventions and module layout may
still change before a 1.0 release.

## Usage

```ts
import { Vector3, /* ... */ } from 'gtengine-js';
```

Not yet published to npm; consume via git until the port stabilizes.

## Development

```
npm install
npm run typecheck
npm test
```

## License

Licensed under the [Boost Software License 1.0](LICENSE), the same license as
the upstream Geometric Tools Engine.

Upstream copyright: David Eberly, Geometric Tools, Redmond WA 98052,
distributed under the Boost Software License, Version 1.0.
