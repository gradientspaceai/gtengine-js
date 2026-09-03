// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Ellipse3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The plane containing the ellipse is Dot(N,X-C) = 0 where X is any point in
// the plane, C is the ellipse center, and N is a unit-length normal to the
// plane. Vectors A0, A1, and N form an orthonormal right-handed set. The
// ellipse in the plane is parameterized by
// X = C + e0*cos(t)*A0 + e1*sin(t)*A1, where A0 is the major axis, A1 is the
// minor axis, and e0 and e1 are the extents along those axes. The angle t is
// in [-pi,pi) and e0 >= e1 > 0.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their Vector arguments,
// comparison methods). The class is not templated on the dimension upstream,
// so the default constructor takes no arguments and builds 3D vectors. The
// C++ 'std::array<Vector3<Real>, 2> axis' becomes a Vector[] of length 2 and
// its comparisons are lexicographic over the elements, as std::array's are.
// 'extent' remains a 2-component Vector.

import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

export class Ellipse3 {
    // Public member access.
    center: Vector;
    normal: Vector;
    axis: Vector[];
    extent: Vector;

    // The port of the default constructor, which sets the center to (0,0,0),
    // A0 to (1,0,0), A1 to (0,1,0), the normal to (0,0,1), e0 to 1 and e1
    // to 1.
    constructor() {
        this.center = new Vector(3);
        this.normal = Vector.unit(3, 2);
        this.axis = [Vector.unit(3, 0), Vector.unit(3, 1)];
        this.extent = Vector.filled(2, 1);
    }

    // The port of 'Ellipse3(inCenter, inNormal, inAxis, inExtent)'. The
    // vectors are copied, matching C++ value semantics.
    static fromCenterNormalAxisExtent(inCenter: Vector, inNormal: Vector,
        inAxis: readonly Vector[], inExtent: Vector): Ellipse3 {
        logAssert(inCenter.size === 3 && inNormal.size === 3,
            'Ellipse3: mismatched sizes.');
        logAssert(inAxis.length === 2 && inAxis[0].size === 3
            && inAxis[1].size === 3, 'Ellipse3: mismatched sizes.');
        logAssert(inExtent.size === 2, 'Ellipse3: mismatched sizes.');
        const ellipse = new Ellipse3();
        ellipse.center = inCenter.clone();
        ellipse.normal = inNormal.clone();
        ellipse.axis = [inAxis[0].clone(), inAxis[1].clone()];
        ellipse.extent = inExtent.clone();
        return ellipse;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Ellipse3 {
        return Ellipse3.fromCenterNormalAxisExtent(this.center, this.normal,
            this.axis, this.extent);
    }

    // Comparisons to support sorted containers.
    equals(ellipse: Ellipse3): boolean {
        return this.center.equals(ellipse.center)
            && this.normal.equals(ellipse.normal)
            && this.axis[0].equals(ellipse.axis[0])
            && this.axis[1].equals(ellipse.axis[1])
            && this.extent.equals(ellipse.extent);
    }

    notEquals(ellipse: Ellipse3): boolean {
        return !this.equals(ellipse);
    }

    lessThan(ellipse: Ellipse3): boolean {
        if (this.center.lessThan(ellipse.center)) {
            return true;
        }

        if (this.center.greaterThan(ellipse.center)) {
            return false;
        }

        if (this.normal.lessThan(ellipse.normal)) {
            return true;
        }

        if (this.normal.greaterThan(ellipse.normal)) {
            return false;
        }

        if (this.axis[0].lessThan(ellipse.axis[0])) {
            return true;
        }

        if (this.axis[0].greaterThan(ellipse.axis[0])) {
            return false;
        }

        if (this.axis[1].lessThan(ellipse.axis[1])) {
            return true;
        }

        if (this.axis[1].greaterThan(ellipse.axis[1])) {
            return false;
        }

        return this.extent.lessThan(ellipse.extent);
    }

    lessThanOrEqual(ellipse: Ellipse3): boolean {
        return !ellipse.lessThan(this);
    }

    greaterThan(ellipse: Ellipse3): boolean {
        return ellipse.lessThan(this);
    }

    greaterThanOrEqual(ellipse: Ellipse3): boolean {
        return !this.lessThan(ellipse);
    }
}
