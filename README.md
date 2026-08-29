# gtengine-js

A TypeScript port of the **Mathematics** component of the
[Geometric Tools Engine (GTE)](https://www.geometrictools.com/) by David Eberly
([upstream repository](https://github.com/davideberly/GeometricTools)).

This is a work-in-progress, machine-assisted port of the ~560 headers in
`GTE/Mathematics`: vector/matrix algebra, geometric primitives, distance
queries, intersection queries, containment, approximation/fitting,
interpolation, curves and surfaces, root finding, computational geometry,
and arbitrary-precision arithmetic.

Port baseline: upstream commit `d29e7758ae2615e5e37da3eb573b7bf90ee94e9b` (2026-08-13).

## Status

The port is in progress. Progress is tracked in the
[issues](https://github.com/gradientspaceai/gtengine-js/issues); each issue is a
batch of upstream headers ported in dependency order. See
[PORTING.md](PORTING.md) for conventions and [porting-status.json](porting-status.json)
for a machine-readable per-file status manifest.

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
