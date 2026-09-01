// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrCircle2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: see IntrIntervals.ts for the Intr* precedent. The circles are
// the ported Hypersphere with dimension 2.

import { Hypersphere } from './Hypersphere';
import { Vector, add, sub, mul, dot, length } from './Vector';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';

// The port of std::numeric_limits<int32_t>::max(), which upstream uses as the
// 'numIntersections' value meaning "the circles are the same".
const INT32_MAX = 2147483647;

// The result of IntrCircle2Circle2TI.test.
export interface IntrCircle2Circle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrCircle2Circle2TIResult {
    return { intersect: false };
}

// The result of IntrCircle2Circle2FI.find.
export interface IntrCircle2Circle2FIResult {
    intersect: boolean;

    // The number of intersections is 0, 1, 2 or 2147483647 (the port of
    // std::numeric_limits<int32_t>::max()). When 1, the circles are tangent
    // and intersect in a single point. When 2, the circles have two
    // transverse intersection points. When 2147483647, the circles are the
    // same.
    numIntersections: number;

    // Valid only when numIntersections is 1 or 2.
    point: [Vector, Vector];

    // Valid only when numIntersections is 2147483647.
    circle: Hypersphere;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrCircle2Circle2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        point: [Vector.zero(2), Vector.zero(2)],
        circle: Hypersphere.fromCenterRadius(Vector.zero(2), 0)
    };
}

export class IntrCircle2Circle2TI implements
    TIQuery<Hypersphere, Hypersphere, IntrCircle2Circle2TIResult> {

    test(circle0: Hypersphere, circle1: Hypersphere): IntrCircle2Circle2TIResult {
        const result = defaultTIResult();
        const diff = sub(circle0.center, circle1.center);
        result.intersect = (length(diff) <= circle0.radius + circle1.radius);
        return result;
    }
}

export class IntrCircle2Circle2FI implements
    FIQuery<Hypersphere, Hypersphere, IntrCircle2Circle2FIResult> {

    find(circle0: Hypersphere, circle1: Hypersphere): IntrCircle2Circle2FIResult {
        // The two circles are |X-C0| = R0 and |X-C1| = R1. Define U = C1 - C0
        // and V = Perp(U) where Perp(x,y) = (y,-x). Note that Dot(U,V) = 0
        // and |V|^2 = |U|^2. The intersection points X can be written in the
        // form X = C0+s*U+t*V and X = C1+(s-1)*U+t*V. Squaring the circle
        // equations and substituting these formulas into them yields
        //   R0^2 = (s^2 + t^2)*|U|^2
        //   R1^2 = ((s-1)^2 + t^2)*|U|^2.
        // Subtracting and solving for s yields
        //   s = ((R0^2-R1^2)/|U|^2 + 1)/2
        // Then replace in the first equation and solve for t^2
        //   t^2 = (R0^2/|U|^2) - s^2.
        // In order for there to be solutions, the right-hand side must be
        // nonnegative. Some algebra leads to the condition for existence of
        // solutions,
        //   (|U|^2 - (R0+R1)^2)*(|U|^2 - (R0-R1)^2) <= 0.
        // This reduces to
        //   |R0-R1| <= |U| <= |R0+R1|.
        // If |U| = |R0-R1|, then the circles are side-by-side and just
        // tangent. If |U| = |R0+R1|, then the circles are nested and just
        // tangent. If |R0-R1| < |U| < |R0+R1|, then the two circles intersect
        // in two points.

        const result = defaultFIResult();

        const U = sub(circle1.center, circle0.center);
        const USqrLen = dot(U, U);
        const R0 = circle0.radius, R1 = circle1.radius;
        const R0mR1 = R0 - R1;
        if (USqrLen === 0 && R0mR1 === 0) {
            // Circles are the same.
            result.intersect = true;
            result.numIntersections = INT32_MAX;
            result.circle = circle0.clone();
            return result;
        }

        const R0mR1Sqr = R0mR1 * R0mR1;
        if (USqrLen < R0mR1Sqr) {
            // The circles do not intersect.
            result.intersect = false;
            result.numIntersections = 0;
            return result;
        }

        const R0pR1 = R0 + R1;
        const R0pR1Sqr = R0pR1 * R0pR1;
        if (USqrLen > R0pR1Sqr) {
            // The circles do not intersect.
            result.intersect = false;
            result.numIntersections = 0;
            return result;
        }

        if (USqrLen < R0pR1Sqr) {
            if (R0mR1Sqr < USqrLen) {
                const invUSqrLen = 1 / USqrLen;
                const s = 0.5 * ((R0 * R0 - R1 * R1) * invUSqrLen + 1);
                const tmp = add(circle0.center, mul(s, U));

                // In theory, discr is nonnegative. However, numerical
                // round-off errors can make it slightly negative. Clamp it to
                // zero.
                let discr = R0 * R0 * invUSqrLen - s * s;
                if (discr < 0) {
                    discr = 0;
                }
                const t = Math.sqrt(discr);
                const V = Vector.fromArray([U.values[1], -U.values[0]]);
                result.point[0] = sub(tmp, mul(t, V));
                result.point[1] = add(tmp, mul(t, V));
                result.numIntersections = (t > 0 ? 2 : 1);
            } else {
                // |U| = |R0-R1|, circles are tangent.
                result.numIntersections = 1;
                result.point[0] = add(circle0.center, mul(R0 / R0mR1, U));
            }
        } else {
            // |U| = |R0+R1|, circles are tangent.
            result.numIntersections = 1;
            result.point[0] = add(circle0.center, mul(R0 / R0pR1, U));
        }

        // The circles intersect in 1 or 2 points.
        result.intersect = true;
        return result;
    }
}
