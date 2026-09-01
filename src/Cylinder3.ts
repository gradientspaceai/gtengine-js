// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The cylinder axis is a line. The origin of the cylinder is chosen to be the
// line origin. The cylinder wall is at a distance R units from the axis. An
// infinite cylinder has infinite height. A finite cylinder has center C at
// the line origin and has a finite height H. The segment for the finite
// cylinder has endpoints C-(H/2)*D and C+(H/2)*D where D is a unit-length
// direction of the line.
//
// NOTE: Some of the geometric queries involve infinite cylinders. To support
// exact arithmetic, it is necessary to avoid infinity/max sentinels. Instead,
// the queries require you to set the infinite cylinder 'height' to -1.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their arguments, comparison
// methods). The class is not templated on the dimension upstream, so the
// default constructor takes no arguments and builds a 3D axis line.

import { logAssert } from './Logger';
import { Line } from './Line';

export class Cylinder3 {
    // Public member access.
    axis: Line;
    radius: number;
    height: number;

    // The port of the default constructor, which sets the axis to the
    // default 3D line (origin (0,0,0), direction (1,0,0)), the radius to 1
    // and the height to 1.
    //
    // NOTE: The upstream comment claims the default constructor sets the
    // axis to (0,0,1), but Line3<T>() has direction (1,0,0). The port
    // preserves the code behavior, not the comment.
    constructor() {
        this.axis = new Line(3);
        this.radius = 1;
        this.height = 1;
    }

    // The port of 'Cylinder3(inAxis, inRadius, inHeight)'. The line is
    // copied, matching C++ value semantics.
    static fromAxisRadiusHeight(inAxis: Line, inRadius: number,
        inHeight: number): Cylinder3 {
        logAssert(inAxis.dimension === 3, 'Cylinder3: mismatched sizes.');
        const cylinder = new Cylinder3();
        cylinder.axis = inAxis.clone();
        cylinder.radius = inRadius;
        cylinder.height = inHeight;
        return cylinder;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Cylinder3 {
        return Cylinder3.fromAxisRadiusHeight(this.axis, this.radius,
            this.height);
    }

    // Please read the NOTE at the beginning of this file about setting the
    // 'height' member for infinite cylinders.
    makeInfiniteCylinder(): void {
        this.height = -1;
    }

    makeFiniteCylinder(inHeight: number): void {
        if (inHeight >= 0) {
            this.height = inHeight;
        }
    }

    isFinite(): boolean {
        return this.height >= 0;
    }

    isInfinite(): boolean {
        return this.height < 0;
    }

    // Comparisons to support sorted containers.
    equals(cylinder: Cylinder3): boolean {
        return this.axis.equals(cylinder.axis)
            && this.radius === cylinder.radius
            && this.height === cylinder.height;
    }

    notEquals(cylinder: Cylinder3): boolean {
        return !this.equals(cylinder);
    }

    lessThan(cylinder: Cylinder3): boolean {
        if (this.axis.lessThan(cylinder.axis)) {
            return true;
        }

        if (this.axis.greaterThan(cylinder.axis)) {
            return false;
        }

        if (this.radius < cylinder.radius) {
            return true;
        }

        if (this.radius > cylinder.radius) {
            return false;
        }

        return this.height < cylinder.height;
    }

    lessThanOrEqual(cylinder: Cylinder3): boolean {
        return !cylinder.lessThan(this);
    }

    greaterThan(cylinder: Cylinder3): boolean {
        return cylinder.lessThan(this);
    }

    greaterThanOrEqual(cylinder: Cylinder3): boolean {
        return !this.lessThan(cylinder);
    }
}
