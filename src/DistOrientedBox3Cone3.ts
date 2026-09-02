// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistOrientedBox3Cone3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between an oriented box and a cone frustum. The
// frustum is part of a single-sided cone with heights measured along the
// axis direction. The single-sided cone heights h satisfy
// 0 <= h <= infinity. The cone frustum has heights that satisfy
// 0 <= hmin < h <= hmax < infinity. The algorithm is described in
// https://www.geometrictools.com/Documentation/DistanceBox3Cone3.pdf
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, OrientedBox3<T>, Cone3<T>>' becomes
// the class DistOrientedBox3Cone3. The nested structs Control and Result
// become the exported interfaces DistOrientedBox3Cone3Control and
// DistOrientedBox3Cone3Result; the upstream default-constructed Control is
// available as defaultDistOrientedBox3Cone3Control(). The private member
// function DoBoxQuadQuery becomes a module-private function. The upstream
// 'Matrix<5,5,T>' objects A and D and the 'Vector<5,T>' objects b and e are
// plain flat number arrays here because they are only used to fill the
// LCP matrix and vector.
//
// Known upstream limitation, preserved by the port. The angle parameter is
// swept over [-pi/2,+pi/2] and the quadrilateral at -pi/2 is the same set as
// the one at +pi/2, so F(-pi/2) and F(+pi/2) agree to within round-off for
// every input. When F is V-shaped over that bracket, Minimize1 fits a
// parabola whose vertex lands at the midpoint of the bracket up to
// round-off. Minimize1 has an exact 'vertex == midpoint' branch that handles
// the symmetric case, but a vertex that is merely near the midpoint takes
// the asymmetric branch instead, which collapses the bracket to a degenerate
// interval and stops. The reported distance is then the distance for some
// valid pair of points, but not necessarily the global minimum. See the
// 'Upstream bug suspects' section of the port PR.

import type { Cone3 } from './Cone';
import type { DCPQuery } from './DCPQuery';
import { LCPSolver } from './LCPSolver';
import { logAssert } from './Logger';
import { Minimize1 } from './Minimize1';
import type { OrientedBox3 } from './OrientedBox';
import { GTE_C_HALF_PI } from './Constants';
import { Vector, add, dot, length, mul, sub } from './Vector';
import { computeOrthogonalComplement3 } from './Vector3';

// Parameters used internally for controlling the minimizer.
export interface DistOrientedBox3Cone3Control {
    maxSubdivisions: number;
    maxBisections: number;
    epsilon: number;
    tolerance: number;
}

// The port of the upstream default-constructed Control.
export function defaultDistOrientedBox3Cone3Control():
    DistOrientedBox3Cone3Control {
    return {
        maxSubdivisions: 8,
        maxBisections: 128,
        epsilon: 1e-08,
        tolerance: 1e-04
    };
}

// The output of the query, which is the distance between the objects and a
// pair of closest points, one from each object.
export interface DistOrientedBox3Cone3Result {
    distance: number;
    boxClosestPoint: Vector;
    coneClosestPoint: Vector;
}

// The result of the box-quadrilateral subproblem for a fixed angle. Upstream
// returns these through reference parameters.
interface BoxQuadResult {
    distance: number;
    boxClosestPoint: Vector;
    quadClosestPoint: Vector;
}

// Compute the distance between the box and the planar quadrilateral obtained
// by slicing the cone frustum with the half plane at the specified angle.
// The problem is posed as a linear complementarity problem of dimension 10.
function doBoxQuadQuery(lcp: LCPSolver, box: OrientedBox3, cone: Cone3,
    coneW0: Vector, coneW1: Vector, quadAngle: number): BoxQuadResult {
    let K = box.center.clone();
    const ell = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        K = sub(K, mul(box.extent.values[i], box.axis[i]));
        ell.values[i] = 2 * box.extent.values[i];
    }

    const cs = Math.cos(quadAngle);
    const sn = Math.sin(quadAngle);
    const term = mul(cone.tanAngle, add(mul(cs, coneW0), mul(sn, coneW1)));
    const G: [Vector, Vector] = [
        sub(cone.ray.direction, term),
        add(cone.ray.direction, term)
    ];

    // A is the 5-by-5 symmetric matrix of the quadratic form, stored
    // row-major.
    const A = new Array<number>(25).fill(0);
    const setA = (r: number, c: number, value: number): void => {
        A[c + 5 * r] = value;
    };
    const getA = (r: number, c: number): number => A[c + 5 * r];

    setA(0, 0, 1);
    setA(0, 1, 0);
    setA(0, 2, 0);
    setA(0, 3, -dot(box.axis[0], G[0]));
    setA(0, 4, -dot(box.axis[0], G[1]));
    setA(1, 0, getA(0, 1));
    setA(1, 1, 1);
    setA(1, 2, 0);
    setA(1, 3, -dot(box.axis[1], G[0]));
    setA(1, 4, -dot(box.axis[1], G[1]));
    setA(2, 0, getA(0, 2));
    setA(2, 1, getA(1, 2));
    setA(2, 2, 1);
    setA(2, 3, -dot(box.axis[2], G[0]));
    setA(2, 4, -dot(box.axis[2], G[1]));
    setA(3, 0, getA(0, 3));
    setA(3, 1, getA(1, 3));
    setA(3, 2, getA(2, 3));
    setA(3, 3, dot(G[0], G[0]));
    setA(3, 4, dot(G[0], G[1]));
    setA(4, 0, getA(0, 4));
    setA(4, 1, getA(1, 4));
    setA(4, 2, getA(2, 4));
    setA(4, 3, getA(3, 4));
    setA(4, 4, dot(G[1], G[1]));

    const KmV = sub(K, cone.ray.origin);
    const b = [
        dot(box.axis[0], KmV),
        dot(box.axis[1], KmV),
        dot(box.axis[2], KmV),
        -dot(G[0], KmV),
        -dot(G[1], KmV)
    ];

    // D is the 5-by-5 matrix of the linear inequality constraints, stored
    // row-major.
    const D = new Array<number>(25).fill(0);
    const getD = (r: number, c: number): number => D[c + 5 * r];
    D[0 + 5 * 0] = -1;
    D[1 + 5 * 1] = -1;
    D[2 + 5 * 2] = -1;
    D[3 + 5 * 3] = +1;
    D[4 + 5 * 3] = +1;
    D[3 + 5 * 4] = -1;
    D[4 + 5 * 4] = -1;

    const e = [
        -ell.values[0],
        -ell.values[1],
        -ell.values[2],
        cone.getMinHeight(),
        -cone.getMaxHeight()
    ];

    const q = new Array<number>(10).fill(0);
    for (let i = 0, ip5 = 5; i < 5; ++i, ++ip5) {
        q[i] = b[i];
        q[ip5] = -e[i];
    }

    const M = new Array<number>(100).fill(0);
    for (let r = 0, rp5 = 5; r < 5; ++r, ++rp5) {
        for (let c = 0, cp5 = 5; c < 5; ++c, ++cp5) {
            M[c + 10 * r] = getA(r, c);
            M[c + 10 * rp5] = getD(r, c);
            M[cp5 + 10 * r] = -getD(c, r);
            M[cp5 + 10 * rp5] = 0;
        }
    }

    const output = lcp.solve(q, M);
    if (output.success) {
        let boxClosestPoint = K;
        for (let i = 0; i < 3; ++i) {
            boxClosestPoint = add(boxClosestPoint,
                mul(output.z[i], box.axis[i]));
        }

        let quadClosestPoint = cone.ray.origin.clone();
        for (let i = 0, ip3 = 3; i < 2; ++i, ++ip3) {
            quadClosestPoint = add(quadClosestPoint, mul(output.z[ip3], G[i]));
        }

        return {
            distance: length(sub(boxClosestPoint, quadClosestPoint)),
            boxClosestPoint,
            quadClosestPoint
        };
    }

    return {
        distance: Number.MAX_VALUE,
        boxClosestPoint: new Vector(3),
        quadClosestPoint: new Vector(3)
    };
}

export class DistOrientedBox3Cone3
    implements DCPQuery<OrientedBox3, Cone3, DistOrientedBox3Cone3Result> {
    private mLCP: LCPSolver = new LCPSolver(10);

    // The default minimizer controls are reasonable choices generally, in
    // which case you can call compute(box, cone). If your application
    // requires specialized controls, pass them as the third argument.
    compute(box: OrientedBox3, cone: Cone3,
        inControl?: DistOrientedBox3Cone3Control):
        DistOrientedBox3Cone3Result {
        const control = inControl !== undefined
            ? { ...inControl } : defaultDistOrientedBox3Cone3Control();

        // Port addition. The header documents the precondition
        // 0 <= hmin < h <= hmax < infinity, but upstream does not check it.
        // An infinite cone stores maxHeight as the sentinel -1, which would
        // silently produce a nonsensical constraint (e[4] = +1) rather than
        // an error, so the port reports the precondition violation.
        logAssert(cone.getMaxHeight() !== -1,
            'The cone must be a frustum with finite maximum height.');

        // Compute a basis for the cone coordinate system.
        const basis: Vector[] = [cone.ray.direction.clone(), new Vector(3),
            new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const coneW0 = basis[1];
        const coneW1 = basis[2];

        const result: DistOrientedBox3Cone3Result = {
            distance: Number.MAX_VALUE,
            boxClosestPoint: new Vector(3),
            coneClosestPoint: new Vector(3)
        };

        const F = (angle: number): number => {
            const bqResult = doBoxQuadQuery(this.mLCP, box, cone, coneW0,
                coneW1, angle);

            if (bqResult.distance < result.distance) {
                result.distance = bqResult.distance;
                result.boxClosestPoint = bqResult.boxClosestPoint;
                result.coneClosestPoint = bqResult.quadClosestPoint;
            }

            return bqResult.distance;
        };

        const minimizer = new Minimize1(F, control.maxSubdivisions,
            control.maxBisections, control.epsilon, control.tolerance);
        const angle0 = -GTE_C_HALF_PI;
        const angle1 = +GTE_C_HALF_PI;
        const minResult = minimizer.getMinimum(angle0, angle1, 0);
        logAssert(minResult.fMin === result.distance,
            'Unexpected mismatch in minimum distance.');

        return result;
    }
}
