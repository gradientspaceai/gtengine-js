// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Sector2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A solid sector is the intersection of a disk and a 2D cone. The disk has
// center C, radius R, and contains points X for which |X-C| <= R. The 2D cone
// has vertex C, unit-length axis direction D, angle A in (0,pi) measured from
// D, and contains points X for which Dot(D,(X-C)/|X-C|) >= cos(A). Sector
// points X satisfy both inequality constraints.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their Vector arguments,
// comparison methods). The class is not templated on the dimension upstream,
// so the default constructor takes no arguments and builds 2D vectors.

import { GTE_C_PI } from './Constants.js';
import { logAssert } from './Logger.js';
import { Vector, sub, dot, length } from './Vector.js';

export class Sector2 {
    // The cosine and sine of the angle are used in queries, so all of angle,
    // cos(angle) and sin(angle) are stored. If you set 'angle' via the public
    // members, you must set all to be consistent. You can also call
    // setAngle(...) to ensure consistency.
    vertex: Vector;
    radius: number;
    direction: Vector;
    angle: number;
    cosAngle: number;
    sinAngle: number;

    // The port of the default constructor, which sets the vertex to (0,0),
    // the radius to 1, the axis direction to (1,0), and the angle to pi, all
    // of which define a disk.
    constructor() {
        this.vertex = new Vector(2);
        this.radius = 1;
        this.direction = Vector.unit(2, 0);
        this.angle = GTE_C_PI;
        this.cosAngle = -1;
        this.sinAngle = 0;
    }

    // The port of 'Sector2(inVertex, inRadius, inDirection, inAngle)'. The
    // vectors are copied, matching C++ value semantics.
    static fromVertexRadiusDirectionAngle(inVertex: Vector, inRadius: number,
        inDirection: Vector, inAngle: number): Sector2 {
        logAssert(inVertex.size === 2 && inDirection.size === 2,
            'Sector2: mismatched sizes.');
        const sector = new Sector2();
        sector.vertex = inVertex.clone();
        sector.radius = inRadius;
        sector.direction = inDirection.clone();
        sector.setAngle(inAngle);
        return sector;
    }

    // A deep copy (the port of C++ copy construction/assignment). The stored
    // cos/sin are copied rather than recomputed so that a hand-set
    // inconsistent state survives the copy, as it does in C++.
    clone(): Sector2 {
        const sector = new Sector2();
        sector.vertex = this.vertex.clone();
        sector.radius = this.radius;
        sector.direction = this.direction.clone();
        sector.angle = this.angle;
        sector.cosAngle = this.cosAngle;
        sector.sinAngle = this.sinAngle;
        return sector;
    }

    // Set the angle, cos(angle) and sin(angle) simultaneously.
    setAngle(inAngle: number): void {
        this.angle = inAngle;
        this.cosAngle = Math.cos(inAngle);
        this.sinAngle = Math.sin(inAngle);
    }

    // Test whether P is in the sector.
    contains(p: Vector): boolean {
        const diff = sub(p, this.vertex);
        const len = length(diff);
        return len <= this.radius
            && dot(this.direction, diff) >= len * this.cosAngle;
    }

    // Comparisons to support sorted containers. Note that upstream compares
    // only vertex, radius, direction and angle (not the derived cos/sin).
    equals(sector: Sector2): boolean {
        return this.vertex.equals(sector.vertex)
            && this.radius === sector.radius
            && this.direction.equals(sector.direction)
            && this.angle === sector.angle;
    }

    notEquals(sector: Sector2): boolean {
        return !this.equals(sector);
    }

    lessThan(sector: Sector2): boolean {
        if (this.vertex.lessThan(sector.vertex)) {
            return true;
        }

        if (this.vertex.greaterThan(sector.vertex)) {
            return false;
        }

        if (this.radius < sector.radius) {
            return true;
        }

        if (this.radius > sector.radius) {
            return false;
        }

        if (this.direction.lessThan(sector.direction)) {
            return true;
        }

        if (this.direction.greaterThan(sector.direction)) {
            return false;
        }

        return this.angle < sector.angle;
    }

    lessThanOrEqual(sector: Sector2): boolean {
        return !sector.lessThan(this);
    }

    greaterThan(sector: Sector2): boolean {
        return sector.lessThan(this);
    }

    greaterThanOrEqual(sector: Sector2): boolean {
        return !this.lessThan(sector);
    }
}
