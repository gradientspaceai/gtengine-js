// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrPlane3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The intersection queries between a plane and a cylinder (finite or
// infinite) are described in
// https://www.geometrictools.com/Documentation/IntersectionCylinderPlane.pdf
//
// The plane is Dot(N, X - P) = 0, where P is a point on the plane and N is a
// nonzero vector that is not necessarily unit length.
//
// The cylinder is (X - C)^T * (I - W * W^T) * (X - C) = r^2, where C is the
// center, W is the axis direction and r > 0 is the radius. The cylinder has
// height h. In the intersection queries, an infinite cylinder is specified by
// setting h = -1 (call cylinder.makeInfiniteCylinder()). This avoids the
// problem of setting the height to the maximum or infinite floating-point
// value, which does not work for exact rational types.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrPlane3Cylinder3TI
// and IntrPlane3Cylinder3FI. The nested 'Result::Type' enum class becomes the
// exported enum IntrPlane3Cylinder3FIResultType (a nested enum is exported
// with a file-qualified name because src/index.ts star-exports every file).
// The private static helpers become module-private functions.

import { Cylinder3 } from './Cylinder3.js';
import { Ellipse3 } from './Ellipse3.js';
import type { FIQuery } from './FIQuery.js';
import { Hyperellipsoid } from './Hyperellipsoid.js';
import { Hyperplane } from './Hyperplane.js';
import type { Plane3 } from './Hyperplane.js';
import { IntrPlane3Plane3FI } from './IntrPlane3Plane3.js';
import { Line } from './Line.js';
import type { Line3 } from './Line.js';
import { logAssert } from './Logger.js';
import { Matrix, mulMatrix, outerProduct, subMatrix } from './Matrix.js';
import { Vector, add, dot, length, mul, sub } from './Vector.js';
import { computeOrthogonalComplement3, cross } from './Vector3.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrPlane3Cylinder3TI.test.
export interface IntrPlane3Cylinder3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrPlane3Cylinder3TIResult():
    IntrPlane3Cylinder3TIResult {
    return { intersect: false };
}

// The type of intersection reported by IntrPlane3Cylinder3FI.
export enum IntrPlane3Cylinder3FIResultType {
    // The cylinder and plane are separated.
    noIntersection = 0,

    // The plane is tangent to the cylinder direction.
    singleLine = 1,

    // The cylinder direction is parallel to the plane and the plane cuts
    // through the cylinder in two lines.
    parallelLines = 2,

    // The cylinder direction is perpendicular to the plane.
    circle = 3,

    // The cylinder direction is not parallel to the plane. When the direction
    // is perpendicular to the plane, the intersection is a circle which is an
    // ellipse with equal extents.
    ellipse = 4
}

// The result of IntrPlane3Cylinder3FI.find. The members are set according to
// 'type'.
//
// type = noIntersection
//   intersect = false
//   line[0,1] and ellipse are default constructed
//
// type = singleLine
//   intersect = true
//   line[0] is valid; line[1] and ellipse are default constructed
//
// type = parallelLines
//   intersect = true
//   line[0] and line[1] are valid; ellipse is default constructed
//
// type = circle
//   intersect = true
//   ellipse is valid (with extent[0] = extent[1]); line[0,1] are default
//   constructed
//
// type = ellipse
//   intersect = true
//   ellipse is valid; line[0,1] are default constructed
export interface IntrPlane3Cylinder3FIResult {
    intersect: boolean;
    type: IntrPlane3Cylinder3FIResultType;
    line: [Line3, Line3];
    ellipse: Ellipse3;

    // Trim lines when the cylinder is finite. They are computed when the
    // plane and infinite cylinder intersect. If there is no intersection, the
    // trim lines are default constructed.
    trimLine: [Line3, Line3];
}

// The port of the upstream FIQuery::Result default constructor. Upstream
// value-initializes the Line3 and Ellipse3 members, so the port uses their
// default constructors.
export function defaultIntrPlane3Cylinder3FIResult():
    IntrPlane3Cylinder3FIResult {
    return {
        intersect: false,
        type: IntrPlane3Cylinder3FIResultType.noIntersection,
        line: [new Line(3), new Line(3)],
        ellipse: new Ellipse3(),
        trimLine: [new Line(3), new Line(3)]
    };
}

// Test-intersection query for a plane and a cylinder in 3D.
export class IntrPlane3Cylinder3TI implements
    TIQuery<Plane3, Cylinder3, IntrPlane3Cylinder3TIResult> {

    // For an infinite cylinder, call cylinder.makeInfiniteCylinder(), which
    // sets the height to -1. For a finite cylinder, set cylinder.height > 0.
    test(plane: Plane3, cylinder: Cylinder3): IntrPlane3Cylinder3TIResult {
        logAssert(plane.dimension === 3,
            'IntrPlane3Cylinder3TI: mismatched sizes.');
        const result = defaultIntrPlane3Cylinder3TIResult();

        // Convenient names.
        const P = plane.origin;
        const N = plane.normal;
        const C = cylinder.axis.origin;
        const W = cylinder.axis.direction;
        const r = cylinder.radius;
        const h = cylinder.height;

        if (cylinder.isInfinite()) {
            if (dot(N, W) !== 0) {
                // The cylinder direction and plane are not parallel.
                result.intersect = true;
            }
            else {
                // The cylinder direction and plane are parallel.
                const dotNCmP = dot(N, sub(C, P));
                result.intersect = (Math.abs(dotNCmP) <= r);
            }
        }
        else {  // the cylinder is finite
            const dotNCmP = dot(N, sub(C, P));
            const dotNW = dot(N, W);
            const crossNW = cross(N, W);
            const lhs = Math.abs(dotNCmP);
            const rhs = r * length(crossNW) + 0.5 * h * Math.abs(dotNW);
            result.intersect = (lhs <= rhs);
        }

        return result;
    }
}

// Find-intersection query for a plane and a cylinder in 3D.
export class IntrPlane3Cylinder3FI implements
    FIQuery<Plane3, Cylinder3, IntrPlane3Cylinder3FIResult> {

    find(plane: Plane3, cylinder: Cylinder3): IntrPlane3Cylinder3FIResult {
        logAssert(plane.dimension === 3,
            'IntrPlane3Cylinder3FI: mismatched sizes.');
        const result = defaultIntrPlane3Cylinder3FIResult();

        const tiQuery = new IntrPlane3Cylinder3TI();
        const tiOutput = tiQuery.test(plane, cylinder);
        if (tiOutput.intersect) {
            const dotNW = dot(plane.normal, cylinder.axis.direction);
            if (dotNW !== 0) {
                // The cylinder direction is not parallel to the plane. The
                // intersection is an ellipse or circle.
                getEllipseOfIntersection(plane, cylinder, result);
                getTrimLines(plane, cylinder, result.trimLine);
            }
            else {
                // The cylinder direction is parallel to the plane. There are
                // no trim lines for this geometric configuration.
                getLinesOfIntersection(plane, cylinder, result);
            }
        }

        return result;
    }
}

// The cylinder is infinite and its direction is not parallel to the plane.
function getEllipseOfIntersection(plane: Plane3, cylinder: Cylinder3,
    result: IntrPlane3Cylinder3FIResult): void {
    // Convenient names.
    const P = plane.origin;
    const N = plane.normal;
    const C = cylinder.axis.origin;
    const W = cylinder.axis.direction;
    const r = cylinder.radius;

    // Compute a right-handed orthonormal basis {N,A,B}. The plane is spanned
    // by A and B. Upstream binds A and B as references into the basis array
    // before the call; the port reads the array after the call, because
    // computeOrthogonalComplement3 assigns new vectors into the slots.
    const basis: Vector[] = [N.clone(), Vector.zero(3), Vector.zero(3)];
    computeOrthogonalComplement3(1, basis);
    const A = basis[1];
    const B = basis[2];

    // Compute the projection matrix M = I - W * W^T.
    const M = subMatrix(Matrix.identity(3, 3), outerProduct(W, W));

    // Compute the coefficients of the quadratic equation
    // c00 + c10*x + c01*y + c20*x^2 + c11*x*y + c02*y^2 = 0.
    const PmC = sub(P, C);
    const MtPmC = mulMatrix(M, PmC) as Vector;
    const MtA = mulMatrix(M, A) as Vector;
    const MtB = mulMatrix(M, B) as Vector;
    const coefficients: number[] = [
        dot(PmC, MtPmC) - r * r,
        2 * dot(A, MtPmC),
        2 * dot(B, MtPmC),
        dot(A, MtA),
        2 * dot(A, MtB),
        dot(B, MtB)
    ];

    // Compute the 2D ellipse parameters in plane coordinates.
    const ellipse2 = new Hyperellipsoid(2);
    ellipse2.fromCoefficients(coefficients);

    // Lift the 2D ellipse/circle to the 3D ellipse/circle.
    result.intersect = true;
    result.type = (ellipse2.extent.values[0] !== ellipse2.extent.values[1]
        ? IntrPlane3Cylinder3FIResultType.ellipse
        : IntrPlane3Cylinder3FIResultType.circle);
    result.ellipse.center = add(plane.origin,
        add(mul(ellipse2.center.values[0], A),
            mul(ellipse2.center.values[1], B)));
    result.ellipse.normal = plane.normal.clone();
    result.ellipse.axis[0] = add(mul(ellipse2.axis[0].values[0], A),
        mul(ellipse2.axis[0].values[1], B));
    result.ellipse.axis[1] = add(mul(ellipse2.axis[1].values[0], A),
        mul(ellipse2.axis[1].values[1], B));
    result.ellipse.extent = ellipse2.extent.clone();
}

// The cylinder is infinite and its direction is parallel to the plane.
function getLinesOfIntersection(plane: Plane3, cylinder: Cylinder3,
    result: IntrPlane3Cylinder3FIResult): void {
    // Convenient names.
    const P = plane.origin;
    const N = plane.normal;
    const C = cylinder.axis.origin;
    const W = cylinder.axis.direction;
    const r = cylinder.radius;

    const CmP = sub(C, P);
    const dotNCmP = dot(N, CmP);
    const ellSqr = r * r - dotNCmP * dotNCmP;  // r^2 - d^2
    if (ellSqr > 0) {
        // The plane cuts through the cylinder in two lines.
        result.intersect = true;
        result.type = IntrPlane3Cylinder3FIResultType.parallelLines;
        const projC = sub(C, mul(dotNCmP, N));
        const crsNW = cross(N, W);
        const ell = Math.sqrt(ellSqr);
        result.line[0].origin = sub(projC, mul(ell, crsNW));
        result.line[0].direction = W.clone();
        result.line[1].origin = add(projC, mul(ell, crsNW));
        result.line[1].direction = W.clone();
    }
    else if (ellSqr < 0) {
        // The cylinder does not intersect the plane.
        result.intersect = false;
        result.type = IntrPlane3Cylinder3FIResultType.noIntersection;
    }
    else {  // ellSqr = 0
        // The plane is tangent to the cylinder.
        result.intersect = true;
        result.type = IntrPlane3Cylinder3FIResultType.singleLine;
        result.line[0].origin = sub(C, mul(dotNCmP, N));
        result.line[0].direction = W.clone();
    }
}

// Compute the intersections of the plane with the two cylinder end planes.
function getTrimLines(plane: Plane3, cylinder: Cylinder3,
    trimLine: [Line3, Line3]): void {
    // Compute the cylinder end planes.
    const C = cylinder.axis.origin;
    const D = cylinder.axis.direction;
    const h = cylinder.height;
    const offset = mul(0.5 * h, D);

    const ppQuery = new IntrPlane3Plane3FI();

    const endPlaneNeg = Hyperplane.fromNormalOrigin(D, sub(C, offset));
    let ppResult = ppQuery.find(plane, endPlaneNeg);
    trimLine[0] = ppResult.line;

    const endPlanePos = Hyperplane.fromNormalOrigin(D, add(C, offset));
    ppResult = ppQuery.find(plane, endPlanePos);
    trimLine[1] = ppResult.line;
}
