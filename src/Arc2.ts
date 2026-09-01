// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The circle containing the arc is represented as |X-C| = r where C is the
// center and r is the radius. The arc is defined by two points E0 and E1 on
// the circle so that E1 is obtained from E0 by traversing counterclockwise.
// The application is responsible for ensuring that E0 and E1 are on the
// circle and that they are properly ordered.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their Vector arguments,
// comparison methods). The class is not templated on the dimension upstream,
// so the default constructor takes no arguments and builds 2D vectors. The
// C++ 'std::array<Vector2<T>, 2> end' becomes a Vector[] of length 2. The
// two 'Contains' overloads become 'contains' (with tolerance) and
// 'containsOnCircle' (assumes P is on the circle).

import { logAssert } from './Logger';
import { Vector, sub, length } from './Vector';
import { dotPerp } from './Vector2';

export class Arc2 {
    // Public member access.
    center: Vector;
    radius: number;
    end: Vector[];

    // The port of the default constructor, which sets the center to (0,0),
    // radius to 1, end0 to (1,0), and end1 to (0,1).
    constructor() {
        this.center = new Vector(2);
        this.radius = 1;
        this.end = [Vector.unit(2, 0), Vector.unit(2, 1)];
    }

    // The port of 'Arc2(C, r, E0, E1)'. The vectors are copied, matching C++
    // value semantics.
    static fromCenterRadiusEnds(c: Vector, r: number, e0: Vector,
        e1: Vector): Arc2 {
        logAssert(c.size === 2 && e0.size === 2 && e1.size === 2,
            'Arc2: mismatched sizes.');
        const arc = new Arc2();
        arc.center = c.clone();
        arc.radius = r;
        arc.end = [e0.clone(), e1.clone()];
        return arc;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Arc2 {
        return Arc2.fromCenterRadiusEnds(this.center, this.radius,
            this.end[0], this.end[1]);
    }

    // Test whether P is on the arc.
    //
    // Formulated for real arithmetic, |P-C| - r = 0 is necessary for P to be
    // on the circle of the arc. If P is on the circle, then P is on the arc
    // from E0 to E1 when it is on the side of the line containing E0 with
    // normal Perp(E1-E0) where Perp(u,v) = (v,-u). This test works for any
    // angle between E0-C and E1-C, even if the angle is larger than or equal
    // to pi radians.
    //
    // Formulated for floating-point types, rounding errors cause |P-C| - r
    // rarely to be 0 when P is on (or numerically near) the circle. To allow
    // for this, choose a small and nonnegative tolerance epsilon. The test
    // concludes that P is on the circle when ||P-C| - r| <= epsilon;
    // otherwise, P is not on the circle. If P is on the circle (in the
    // epsilon-tolerance sense), the side-of-line test of the previous
    // paragraph is applied.
    //
    // NOTE: The upstream comment claims that a negative epsilon behaves as if
    // zero was passed. It does not: ||P-C| - r| is nonnegative, so the test
    // ||P-C| - r| <= epsilon fails for every P when epsilon < 0 and the
    // function returns false. The port preserves the upstream code behavior,
    // not the comment.
    contains(p: Vector, epsilon: number): boolean {
        const len = length(sub(p, this.center));
        if (Math.abs(len - this.radius) <= epsilon) {
            return this.containsOnCircle(p);
        }
        return false;
    }

    // The port of the single-argument 'Contains'. This function assumes P is
    // on the circle containing the arc (with possibly a small amount of
    // floating-point rounding error).
    containsOnCircle(p: Vector): boolean {
        const diffPE0 = sub(p, this.end[0]);
        const diffE1E0 = sub(this.end[1], this.end[0]);
        return dotPerp(diffPE0, diffE1E0) >= 0;
    }

    // Comparisons to support sorted containers.
    equals(arc: Arc2): boolean {
        return this.center.equals(arc.center)
            && this.radius === arc.radius
            && this.end[0].equals(arc.end[0])
            && this.end[1].equals(arc.end[1]);
    }

    notEquals(arc: Arc2): boolean {
        return !this.equals(arc);
    }

    lessThan(arc: Arc2): boolean {
        if (this.center.lessThan(arc.center)) {
            return true;
        }

        if (this.center.greaterThan(arc.center)) {
            return false;
        }

        if (this.radius < arc.radius) {
            return true;
        }

        if (this.radius > arc.radius) {
            return false;
        }

        if (this.end[0].lessThan(arc.end[0])) {
            return true;
        }

        if (this.end[0].greaterThan(arc.end[0])) {
            return false;
        }

        return this.end[1].lessThan(arc.end[1]);
    }

    lessThanOrEqual(arc: Arc2): boolean {
        return !arc.lessThan(this);
    }

    greaterThan(arc: Arc2): boolean {
        return arc.lessThan(this);
    }

    greaterThanOrEqual(arc: Arc2): boolean {
        return !this.lessThan(arc);
    }
}
