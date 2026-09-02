// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContScribeCircle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Circumscribed and inscribed circles for a triangle in 2D.
//
// Port notes: upstream overloads Circumscribe/Inscribe on the container type
// written through the output reference. Following the Cont* naming precedent,
// the port suffixes each function with that container type
// (circumscribeCircle2, inscribeCircle2). Upstream returns 'true' when the
// circle was constructed and writes it through an output parameter; the port
// returns the circle, or null when the input points are linearly dependent.

import { Hypersphere, type Circle2 } from './Hypersphere';
import { LinearSystem } from './LinearSystem';
import { logAssert } from './Logger';
import { Matrix } from './Matrix';
import { Vector, add, dot, length, mul, sub } from './Vector';
import { dotPerp } from './Vector2';

// Circle circumscribing a triangle. Returns null when the input points are
// linearly dependent (a degenerate triangle).
export function circumscribeCircle2(v0: Vector, v1: Vector,
    v2: Vector): Circle2 | null {
    logAssert(v0.size === 2 && v1.size === 2 && v2.size === 2,
        'circumscribeCircle2: points must be 2D.');

    const e10 = sub(v1, v0);
    const e20 = sub(v2, v0);
    const A = Matrix.zero(2, 2);
    A.setRow(0, e10);
    A.setRow(1, e20);
    const B = Vector.fromArray([0.5 * dot(e10, e10), 0.5 * dot(e20, e20)]);
    const { X, invertible } = LinearSystem.solve2x2(A, B);
    if (invertible) {
        const circle = new Hypersphere(2);
        circle.center = add(v0, X);
        circle.radius = length(X);
        return circle;
    }
    return null;
}

// Circle inscribing a triangle. Returns null when the perimeter is zero or
// when the triangle is degenerate (zero inradius).
export function inscribeCircle2(v0: Vector, v1: Vector,
    v2: Vector): Circle2 | null {
    logAssert(v0.size === 2 && v1.size === 2 && v2.size === 2,
        'inscribeCircle2: points must be 2D.');

    const d10 = sub(v1, v0);
    const d20 = sub(v2, v0);
    const d21 = sub(v2, v1);
    let len10 = length(d10);
    let len20 = length(d20);
    let len21 = length(d21);
    const perimeter = len10 + len20 + len21;
    if (perimeter > 0) {
        const inv = 1 / perimeter;
        len10 *= inv;
        len20 *= inv;
        len21 *= inv;

        // The incenter is the perimeter-weighted average of the vertices,
        // each weighted by the length of the opposite edge.
        const circle = new Hypersphere(2);
        circle.center = add(add(mul(len21, v0), mul(len20, v1)),
            mul(len10, v2));

        // The inradius is Area/s, where s is the semiperimeter and the area
        // is |DotPerp(d10,d20)|/2.
        circle.radius = inv * Math.abs(dotPerp(d10, d20));
        return circle.radius > 0 ? circle : null;
    }
    return null;
}
