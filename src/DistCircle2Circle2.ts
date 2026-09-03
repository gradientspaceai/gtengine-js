// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistCircle2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between two circles in 2D. The circles are considered
// to be curves, not solid disks.
//
// The circles are C[i] + r[i] * U(s[i]) for i in {0,1}, where C[i] is the
// center, r[i] > 0 is the radius, and U(s[i]) = (cos(s[i]), sin(s[i])) for
// s[i] in [0,2*pi). The circles are concentric when C[0] = C[1]. The circles
// are cocircular if they are concentric and r[0] = r[1].
//
// The number of pairs of closest points is result.numClosestPairs, which is
// 1 or 2.
//
// If result.numClosestPairs is 1, the circle[i] closest points are
// result.closest[0][i]. The possible geometric configurations are
//   1. The circles are strictly separated with positive distance.
//   2. The circles are separated but tangent at a point.
//   3. One circle is strictly inside the other circle with positive distance.
//   4. One circle is inside the other circle but tangent at a point.
//
// If result.numClosestPairs is 2 and the circles are not concentric, there
// are 2 pairs of circle[i] closest points, which are the intersection points
// of the circles.
//
// If the circles are concentric or cocircular, there are infinitely many
// pairs of closest points. The distance is |r[1]-r[0]|. The reported pairs
// are (C[0] - r[0] * (1,0), C[1] - r[1] * (1,0)) and
// (C[0] + r[0] * (1,0), C[1] + r[1] * (1,0)). The number of pairs reported
// is 2.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Circle2<T>, Circle2<T>>' becomes the
// class DistCircle2Circle2 with the result type DistCircle2Circle2Result. The
// private helper 'DoQuery' becomes the private method doQuery.

import type { DCPQuery } from './DCPQuery.js';
import type { Circle2 } from './Hypersphere.js';
import { Vector, add, length, mul, normalize, sub } from './Vector.js';
import { perp } from './Vector2.js';

export interface DistCircle2Circle2Result {
    distance: number;
    sqrDistance: number;

    // The number of pairs of closest points, 1 or 2.
    numClosestPairs: number;

    // closest[j][i] is the point on circle[i] of the j-th closest pair.
    closest: [[Vector, Vector], [Vector, Vector]];

    concentric: boolean;
    cocircular: boolean;
}

function defaultResult(): DistCircle2Circle2Result {
    return {
        distance: 0,
        sqrDistance: 0,
        numClosestPairs: 0,
        closest: [
            [new Vector(2), new Vector(2)],
            [new Vector(2), new Vector(2)]
        ],
        concentric: false,
        cocircular: false
    };
}

export class DistCircle2Circle2
    implements DCPQuery<Circle2, Circle2, DistCircle2Circle2Result> {
    compute(circle0: Circle2, circle1: Circle2): DistCircle2Circle2Result {
        const result = defaultResult();
        if (circle0.radius >= circle1.radius) {
            this.doQuery(circle0, circle1, result);
        }
        else {
            this.doQuery(circle1, circle0, result);
        }
        return result;
    }

    // The query requires circle0.radius >= circle1.radius.
    private doQuery(circle0: Circle2, circle1: Circle2,
        result: DistCircle2Circle2Result): void {
        if (circle0.center.notEquals(circle1.center)) {
            const delta = sub(circle1.center, circle0.center);
            const lenDelta = length(delta);
            const rSum = circle0.radius + circle1.radius;
            let distance = lenDelta - rSum;
            if (distance >= 0) {
                // Configurations 1 or 2 for the case numClosestPairs = 1.
                // Configuration 2 occurs when lenDelta equals rSum.
                result.distance = distance;
                result.sqrDistance = distance * distance;
                result.numClosestPairs = 1;
                normalize(delta);
                result.closest[0][0] =
                    add(circle0.center, mul(circle0.radius, delta));
                if (distance > 0) {
                    result.closest[0][1] =
                        sub(circle1.center, mul(circle1.radius, delta));
                }
                else {
                    result.closest[0][1] = result.closest[0][0].clone();
                }
            }
            else {
                const rDif = circle0.radius - circle1.radius;
                distance = rDif - lenDelta;
                if (distance >= 0) {
                    result.distance = distance;
                    result.sqrDistance = distance * distance;
                    result.numClosestPairs = 1;
                    normalize(delta);
                    result.closest[0][0] =
                        add(circle0.center, mul(circle0.radius, delta));
                    if (distance > 0) {
                        result.closest[0][1] =
                            add(circle1.center, mul(circle1.radius, delta));
                    }
                    else {
                        result.closest[0][1] = result.closest[0][0].clone();
                    }
                }
                else {
                    // Let D = C1 - C0. The circles intersect at the points
                    // X = C0 + u * D + v * Perp(D) for some u in (0,1). We
                    // know r0^2 = |X-C0|^2 and r1^2 = |X-C1|^2, which leads to
                    //   r0^2 = (u^2 + v^2) * |D|^2
                    //   r1^2 = ((u-1)^2 + v^2) * |D|^2
                    // The solutions are u = (1 + (r0^2 - r1^2)/|D|^2)/2 and
                    // v = +/- (r0^2 / |D|^2 - u^2).
                    const rSumDivLen = rSum / lenDelta;
                    const rDifDivLen = rDif / lenDelta;
                    const r0DivLen = circle0.radius / lenDelta;
                    const u = 0.5 * (1 + rSumDivLen * rDifDivLen);
                    const v = Math.sqrt(
                        Math.max(r0DivLen * r0DivLen - u * u, 0));

                    result.distance = 0;
                    result.sqrDistance = 0;
                    result.numClosestPairs = 2;
                    const temp0 = add(circle0.center, mul(u, delta));
                    const temp1 = mul(v, perp(delta));
                    result.closest[0][0] = add(temp0, temp1);
                    result.closest[0][1] = result.closest[0][0].clone();
                    result.closest[1][0] = sub(temp0, temp1);
                    result.closest[1][1] = result.closest[1][0].clone();
                }
            }
        }
        else {
            result.distance = Math.abs(circle0.radius - circle1.radius);
            result.sqrDistance = result.distance * result.distance;
            result.numClosestPairs = 2;
            const offset0 = Vector.fromArray([circle0.radius, 0]);
            const offset1 = Vector.fromArray([circle1.radius, 0]);
            result.closest[0][0] = sub(circle0.center, offset0);
            result.closest[0][1] = sub(circle1.center, offset1);
            result.closest[1][0] = add(circle0.center, offset0);
            result.closest[1][1] = add(circle1.center, offset1);
            result.concentric = true;
            result.cocircular = (circle0.radius === circle1.radius);
        }
    }
}
