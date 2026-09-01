// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Torus3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A torus with origin (0,0,0), outer radius r0 and inner radius r1 (with
// r0 >= r1) is defined implicitly as follows. The point P0 = (x,y,z) is on
// the torus. Its projection onto the xy-plane is P1 = (x,y,0). The circular
// cross section of the torus that contains the projection has radius r0 and
// center P2 = r0*(x,y,0)/sqrt(x^2+y^2). The triangle <P0,P1,P2> is a right
// triangle with right angle at P1. The hypotenuse <P0,P2> has length r1, leg
// <P1,P2> has length z and leg <P0,P1> has length |r0 - sqrt(x^2+y^2)|. The
// Pythagorean theorem says z^2 + |r0 - sqrt(x^2+y^2)|^2 = r1^2. This can be
// algebraically manipulated to
//   (x^2 + y^2 + z^2 + r0^2 - r1^2)^2 - 4 * r0^2 * (x^2 + y^2) = 0
//
// A parametric form is
//   x = (r0 + r1 * cos(v)) * cos(u)
//   y = (r0 + r1 * cos(v)) * sin(u)
//   z = r1 * sin(v)
// for u in [0,2*pi) and v in [0,2*pi).
//
// Generally, let the torus center be C with plane of symmetry containing C
// and having directions D0 and D1. The axis of symmetry is the line
// containing C and having direction N (the plane normal). The radius from the
// center of the torus is r0 and the radius of the tube of the torus is r1. A
// point P may be written as P = C + x*D0 + y*D1 + z*N, where matrix [D0 D1 N]
// is orthonormal and has determinant 1. Thus, x = Dot(D0,P-C),
// y = Dot(D1,P-C) and z = Dot(N,P-C). The implicit form is
//      [|P-C|^2 + r0^2 - r1^2]^2 - 4*r0^2*[|P-C|^2 - (Dot(N,P-C))^2] = 0
// Observe that D0 and D1 are not present in the equation, which is to be
// expected by the symmetry. The parametric form is
//      P(u,v) = C + (r0 + r1*cos(v))*(cos(u)*D0 + sin(u)*D1) + r1*sin(v)*N
// for u in [0,2*pi) and v in [0,2*pi).
//
// In the class Torus3, the members are 'center' C, 'direction0' D0,
// 'direction1' D1, 'normal' N, 'radius0' r0 and 'radius1' r1.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their Vector arguments,
// comparison methods). The class is not templated on the dimension upstream,
// so the default constructor takes no arguments and builds 3D vectors.
// 'Evaluate' returns the jet array rather than filling a caller-supplied
// buffer, and 'GetParameters' returns { u, v } instead of writing to output
// references.

import { logAssert } from './Logger';
import { Vector, add, sub, mul, dot, negate } from './Vector';

export class Torus3 {
    // Public member access.
    center: Vector;
    direction0: Vector;
    direction1: Vector;
    normal: Vector;
    radius0: number;
    radius1: number;

    // The port of the default constructor, which sets the center to (0,0,0),
    // direction0 to (1,0,0), direction1 to (0,1,0), the normal to (0,0,1),
    // radius0 to 2 and radius1 to 1.
    constructor() {
        this.center = new Vector(3);
        this.direction0 = Vector.unit(3, 0);
        this.direction1 = Vector.unit(3, 1);
        this.normal = Vector.unit(3, 2);
        this.radius0 = 2;
        this.radius1 = 1;
    }

    // The port of the six-argument constructor. The vectors are copied,
    // matching C++ value semantics.
    static fromCenterFrameRadii(inCenter: Vector, inDirection0: Vector,
        inDirection1: Vector, inNormal: Vector, inRadius0: number,
        inRadius1: number): Torus3 {
        logAssert(inCenter.size === 3 && inDirection0.size === 3
            && inDirection1.size === 3 && inNormal.size === 3,
            'Torus3: mismatched sizes.');
        const torus = new Torus3();
        torus.center = inCenter.clone();
        torus.direction0 = inDirection0.clone();
        torus.direction1 = inDirection1.clone();
        torus.normal = inNormal.clone();
        torus.radius0 = inRadius0;
        torus.radius1 = inRadius1;
        return torus;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Torus3 {
        return Torus3.fromCenterFrameRadii(this.center, this.direction0,
            this.direction1, this.normal, this.radius0, this.radius1);
    }

    // Evaluation of the surface. The function supports derivative calculation
    // through order 2; that is, maxOrder <= 2 is required. If you want only
    // the position, pass in maxOrder of 0. If you want the position and
    // first-order derivatives, pass in maxOrder of 1, and so on. The returned
    // 'jet' values are ordered as: position X; first-order derivatives dX/du,
    // dX/dv; second-order derivatives d2X/du2, d2X/dudv, d2X/dv2.
    evaluate(u: number, v: number, maxOrder: number): Vector[] {
        const jet: Vector[] = [];

        // Compute position.
        const csu = Math.cos(u);
        const snu = Math.sin(u);
        const csv = Math.cos(v);
        const snv = Math.sin(v);
        const r1csv = this.radius1 * csv;
        const r1snv = this.radius1 * snv;
        const r0pr1csv = this.radius0 + r1csv;
        const combo0 = add(mul(csu, this.direction0),
            mul(snu, this.direction1));
        const r0pr1csvcombo0 = mul(r0pr1csv, combo0);
        const r1snvnormal = mul(r1snv, this.normal);
        jet[0] = add(add(this.center, r0pr1csvcombo0), r1snvnormal);

        if (maxOrder >= 1) {
            // Compute first-order derivatives.
            const combo1 = add(mul(-snu, this.direction0),
                mul(csu, this.direction1));
            jet[1] = mul(r0pr1csv, combo1);
            jet[2] = add(mul(-r1snv, combo0), mul(r1csv, this.normal));

            if (maxOrder === 2) {
                // Compute second-order derivatives.
                jet[3] = negate(r0pr1csvcombo0);
                jet[4] = mul(-r1snv, combo1);
                jet[5] = sub(mul(-r1csv, combo0), r1snvnormal);
            }
        }

        return jet;
    }

    // Reverse lookup of parameters from position.
    getParameters(x: Vector): { u: number, v: number } {
        const delta = sub(x, this.center);

        // (r0 + r1*cos(v))*cos(u)
        const dot0 = dot(this.direction0, delta);

        // (r0 + r1*cos(v))*sin(u)
        const dot1 = dot(this.direction1, delta);

        // r1*sin(v)
        const dot2 = dot(this.normal, delta);

        // r1*cos(v)
        const r1csv = Math.sqrt(dot0 * dot0 + dot1 * dot1) - this.radius0;

        return { u: Math.atan2(dot1, dot0), v: Math.atan2(dot2, r1csv) };
    }

    // Comparisons to support sorted containers.
    equals(torus: Torus3): boolean {
        return this.center.equals(torus.center)
            && this.direction0.equals(torus.direction0)
            && this.direction1.equals(torus.direction1)
            && this.normal.equals(torus.normal)
            && this.radius0 === torus.radius0
            && this.radius1 === torus.radius1;
    }

    notEquals(torus: Torus3): boolean {
        return !this.equals(torus);
    }

    lessThan(torus: Torus3): boolean {
        if (this.center.lessThan(torus.center)) {
            return true;
        }

        if (this.center.greaterThan(torus.center)) {
            return false;
        }

        if (this.direction0.lessThan(torus.direction0)) {
            return true;
        }

        if (this.direction0.greaterThan(torus.direction0)) {
            return false;
        }

        if (this.direction1.lessThan(torus.direction1)) {
            return true;
        }

        if (this.direction1.greaterThan(torus.direction1)) {
            return false;
        }

        if (this.normal.lessThan(torus.normal)) {
            return true;
        }

        if (this.normal.greaterThan(torus.normal)) {
            return false;
        }

        if (this.radius0 < torus.radius0) {
            return true;
        }

        if (this.radius0 > torus.radius0) {
            return false;
        }

        return this.radius1 < torus.radius1;
    }

    lessThanOrEqual(torus: Torus3): boolean {
        return !torus.lessThan(this);
    }

    greaterThan(torus: Torus3): boolean {
        return torus.lessThan(this);
    }

    greaterThanOrEqual(torus: Torus3): boolean {
        return !this.lessThan(torus);
    }
}
