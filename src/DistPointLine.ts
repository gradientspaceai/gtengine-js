// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointLine.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a line in nD.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The input point is stored in closest[0]. The closest point on the line is
// stored in closest[1].
//
// Port notes (this file is part of the first Dist* batch; see DCPQuery.ts for
// the family conventions):
//
// * Upstream defines the specialization
//   'DCPQuery<T, Vector<N,T>, Line<N,T>>' (aliased DCPPointLine). TypeScript
//   has no template specialization, so the port exports a concrete class
//   named after the file, DistPointLine, implementing DCPQuery.
// * The nested 'Result' struct becomes the exported interface
//   DistPointLineResult; src/index.ts star-exports every file, so result
//   types are file-qualified.
// * 'operator()' becomes 'compute(...)' per PORTING.md.
// * The dimension N is a runtime property of the Vector/Line arguments, so
//   the upstream aliases DCPPoint2Line2/DCPPoint3Line3 have no port; the one
//   class serves every dimension.

import type { DCPQuery } from './DCPQuery.js';
import type { Line } from './Line.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

export interface DistPointLineResult {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // closest[0] is the input point, closest[1] is the closest line point.
    closest: [Vector, Vector];
}

export class DistPointLine
    implements DCPQuery<Vector, Line, DistPointLineResult> {
    compute(point: Vector, line: Line): DistPointLineResult {
        let diff = sub(point, line.origin);
        const parameter = dot(line.direction, diff)
            / dot(line.direction, line.direction);
        const closest0 = point.clone();
        const closest1 = add(line.origin, mul(parameter, line.direction));
        diff = sub(closest0, closest1);
        const sqrDistance = dot(diff, diff);

        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            parameter,
            closest: [closest0, closest1]
        };
    }
}
