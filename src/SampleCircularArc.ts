// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SampleCircularArc.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithm for representing an arc as a NURBS curve is described in
//   https://www.geometrictools.com/Documentation/NURBSCircleSphere.pdf
// The SampleCircularArc class generates points on an arc. The arc must be
// counterclockwise ordered. The number of returned points is the approximate
// length of the arc. This is useful for 2D applications where you want to
// draw an arc. The alternative is to derive an algorithm for integer-based
// pixel selection similar to Bresenham's algorithm for a full circle.
//
// Port notes: 'operator()' becomes 'compute'. Upstream resizes a
// caller-supplied 'std::vector<Vector2<T>>&'; the port returns the points.
// The subarc helpers write into a subrange of that vector through a pointer;
// the port passes the destination array and a starting offset. The C++
// 'static_cast<std::size_t>(radius * angle)' truncates toward zero, so the
// port uses Math.trunc; the sample count is therefore floor(radius*angle),
// floor(radius*angle/2) per half, and so on.
//
// Sampling semantics: each subarc is sampled at u = i/numPoints for
// 0 <= i < numPoints, so the first point of the arc is included but the final
// endpoint arc.end[1] is not. An arc short enough that the count truncates to
// zero produces no points.

import { Arc2 } from './Arc2.js';
import { GTE_C_TWO_PI } from './Constants.js';
import { Vector, dot, normalize } from './Vector.js';
import { dotPerp, perp } from './Vector2.js';

export class SampleCircularArc {
    compute(arc: Arc2): Vector[] {
        // Translate and scale the arc to the unit circle centered at the
        // origin. Compute the angle subtended by the arc.
        const twoPi = GTE_C_TWO_PI;
        const P0 = Vector.fromArray([
            (arc.end[0].values[0] - arc.center.values[0]) / arc.radius,
            (arc.end[0].values[1] - arc.center.values[1]) / arc.radius
        ]);
        const P2 = Vector.fromArray([
            (arc.end[1].values[0] - arc.center.values[0]) / arc.radius,
            (arc.end[1].values[1] - arc.center.values[1]) / arc.radius
        ]);
        const d = Math.max(-1, Math.min(dot(P0, P2), 1));
        const angle = Math.acos(d);

        // Decompose the arc into subarcs, each with subtended angle in
        // (0,pi/2].
        const dp = dotPerp(P0, P2);
        if (dp >= 0) {
            // The subtended angle is in [0,pi].
            if (d >= 0) {
                // The subtended angle is in [0,pi/2].
                return this.sampleArc1(arc.center, arc.radius, P0, P2, angle);
            }
            else {
                // The subtended angle is in [pi/2,pi].
                return this.sampleArc2(arc.center, arc.radius, P0, P2, angle);
            }
        }
        else {
            // The subtended angle is in [pi,2*pi].
            if (d <= 0) {
                // The subtended angle is in [pi,3*pi/2].
                return this.sampleArc3(arc.center, arc.radius, P0, P2,
                    twoPi - angle);
            }
            else {
                // The subtended angle is in [3*pi/2,2*pi].
                return this.sampleArc4(arc.center, arc.radius, P0, P2,
                    twoPi - angle);
            }
        }
    }

    private sampleArc1(center: Vector, radius: number, P0: Vector, P2: Vector,
        angle: number): Vector[] {
        const numPoints = Math.trunc(radius * angle);
        const points = new Array<Vector>(numPoints);
        this.sampleAcuteArc(center, radius, P0, P2, numPoints, points, 0);
        return points;
    }

    private sampleArc2(center: Vector, radius: number, P0: Vector, P2: Vector,
        angle: number): Vector[] {
        const numPoints = Math.trunc(radius * angle / 2);
        const points = new Array<Vector>(2 * numPoints);
        const bisector0 = addVec(P0, P2);
        normalize(bisector0);
        this.sampleAcuteArc(center, radius, P0, bisector0, numPoints, points, 0);
        this.sampleAcuteArc(center, radius, bisector0, P2, numPoints, points,
            numPoints);
        return points;
    }

    private sampleArc3(center: Vector, radius: number, P0: Vector, P2: Vector,
        angle: number): Vector[] {
        const numPoints = Math.trunc(radius * angle / 3);
        const points = new Array<Vector>(3 * numPoints);
        const trisector0 = negScaleAdd(2, P0, 1, P2);
        normalize(trisector0);
        const trisector1 = negScaleAdd(1, P0, 2, P2);
        normalize(trisector1);
        this.sampleAcuteArc(center, radius, P0, trisector0, numPoints, points, 0);
        this.sampleAcuteArc(center, radius, trisector0, trisector1, numPoints,
            points, numPoints);
        this.sampleAcuteArc(center, radius, trisector1, P2, numPoints, points,
            2 * numPoints);
        return points;
    }

    private sampleArc4(center: Vector, radius: number, P0: Vector, P2: Vector,
        angle: number): Vector[] {
        const numPoints = Math.trunc(radius * angle / 4);
        const points = new Array<Vector>(4 * numPoints);
        const quadsector0 = negScaleAdd(3, P0, 1, P2);
        normalize(quadsector0);
        const quadsector1 = negScaleAdd(1, P0, 1, P2);
        normalize(quadsector1);
        const quadsector2 = negScaleAdd(1, P0, 3, P2);
        normalize(quadsector2);
        this.sampleAcuteArc(center, radius, P0, quadsector0, numPoints, points, 0);
        this.sampleAcuteArc(center, radius, quadsector0, quadsector1, numPoints,
            points, numPoints);
        this.sampleAcuteArc(center, radius, quadsector1, quadsector2, numPoints,
            points, 2 * numPoints);
        this.sampleAcuteArc(center, radius, quadsector2, P2, numPoints, points,
            3 * numPoints);
        return points;
    }

    // The preconditions are:
    // 1. The arc is on the unit circle centered at the origin.
    // 2. The angle subtended by the arc is in the interval (0,pi/2].
    // 3. 'points' has enough storage from 'offset' for the requested number
    //    of samples.
    private sampleAcuteArc(center: Vector, radius: number, P0: Vector,
        P2: Vector, numPoints: number, points: Vector[], offset: number): void {
        // The ordered points of the arc are {P0,P1,P2}.
        const diff = Vector.fromArray([
            P2.values[0] - P0.values[0], P2.values[1] - P0.values[1]
        ]);
        const denom = dotPerp(P0, P2);
        const pp = perp(diff);
        const P1 = Vector.fromArray([
            pp.values[0] / denom, pp.values[1] / denom
        ]);

        // Compute the NURBS weights for the parameterization. The weights
        // w1 = 1 and w2 = w0.
        const w0 = Math.sqrt(2 * (dot(P1, P1) - 1) / (1 - dot(P0, P2)));

        // Compute the NURBS control points for the parameterization.
        const C0x = center.values[0] + radius * P0.values[0];
        const C0y = center.values[1] + radius * P0.values[1];
        const C1x = center.values[0] + radius * P1.values[0];
        const C1y = center.values[1] + radius * P1.values[1];
        const C2x = center.values[0] + radius * P2.values[0];
        const C2y = center.values[1] + radius * P2.values[1];

        // Compute the samples for u in [0,1).
        for (let i = 0; i < numPoints; ++i) {
            const u = i / numPoints;
            const onemu = 1 - u;
            const k0 = w0 * onemu * onemu;
            const k1 = 2 * u * onemu;
            const k2 = w0 * u * u;
            const k = k0 + k1 + k2;
            points[offset + i] = Vector.fromArray([
                (k0 * C0x + k1 * C1x + k2 * C2x) / k,
                (k0 * C0y + k1 * C1y + k2 * C2y) / k
            ]);
        }
    }
}

function addVec(v0: Vector, v1: Vector): Vector {
    return Vector.fromArray([
        v0.values[0] + v1.values[0], v0.values[1] + v1.values[1]
    ]);
}

// Compute -(a0 * v0 + a1 * v1) for 2D vectors.
function negScaleAdd(a0: number, v0: Vector, a1: number, v1: Vector): Vector {
    return Vector.fromArray([
        -(a0 * v0.values[0] + a1 * v1.values[0]),
        -(a0 * v0.values[1] + a1 * v1.values[1])
    ]);
}
