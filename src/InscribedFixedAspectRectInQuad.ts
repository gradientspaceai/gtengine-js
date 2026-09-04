// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) InscribedFixedAspectRectInQuad.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the maximum-area, fixed-aspect-ratio, and axis-aligned rectangle
// inscribed in a convex quadrilateral. The algorithm is described in
// https://www.geometrictools.com/Documentation/MaximumAreaAspectRectangle.pdf
//
// KNOWN UPSTREAM LIMITATION (preserved): Execute solves the linear program by
// intersecting the constraint planes 0 and 2 and clipping the resulting line
// with the constraints 1 and 3 (falling back to the complementary pairing when
// that clip is empty). It asserts alpha1 != 0 and alpha3 != 0, i.e. that the
// line is not parallel to the other two planes. That assumption fails for
// legitimate convex quadrilaterals. For a fixed case index j, the constraint
// normal is a linear function of the 2D edge normal (n0,n1):
//   j = 0: (n0, n1, 0)                  j = 1: (n0, n1, n0)
//   j = 2: (n0, n1, n0 + n1/r)          j = 3: (n0, n1, n1/r)
// so all constraints with the same j lie in a common plane through the origin
// of (u,v,w)-space. Whenever three of the four edge normals fall in the same
// quadrant (in particular normals 0, 2 and one of 1 or 3), those three
// constraint normals are coplanar, the corresponding alpha is exactly zero and
// the assert fires with "Unexpected condition". The smallest lattice example
// is the counterclockwise convex quad <(0,3),(0,2),(1,1),(3,0)> (normals 0, 1
// and 2 all in quadrant 0), which has a well-defined maximum inscribed
// rectangle for every aspect ratio. The port preserves the upstream behavior
// (a thrown Error) rather than redesigning the solver.
//
// A second, round-off-driven failure of the same kind: when the maximum
// rectangle touches all four edges, the feasible interval along the line of
// planes 0 and 2 degenerates to a single point. The endpoints -beta1/alpha1
// and -beta3/alpha3 are compared without a tolerance, so one ulp of error can
// order them the wrong way; the interval is then reported empty, the
// complementary pairing is empty for the same reason, and Execute reports
// "Unexpected interval intersection type" for a solvable problem. Both
// failures are pinned by tests in test/InscribedFixedAspectRectInQuad.test.ts.
//
// Port notes: the upstream static Execute has the output parameters
// 'rectOrigin', 'rectWidth' and 'rectHeight' and returns a 'bool'; the port
// returns an object with the fields 'isUnique', 'rectOrigin', 'rectWidth' and
// 'rectHeight'. The upstream 'std::pair<Vector3<T>, T>' constraints become
// objects with the fields 'normal' and 'constant'. The upstream nested type
// alias IIQuery (FIQuery for two intervals) becomes IntrIntervalsFI and the
// semiinfinite-semiinfinite operator() overload is the ported method
// findSemiInfiniteSemiInfinite.

import { GTE_C_INV_HALF_PI, GTE_C_TWO_PI } from './Constants.js';
import { IntrIntervalsFI, IntrIntervalsFIResultType } from './IntrIntervals.js';
import { logAssert, logError } from './Logger.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { perp } from './Vector2.js';
import { cross } from './Vector3.js';

export interface InscribedFixedAspectRectInQuadResult {
    // The value is 'true' when there is a unique solution or 'false' when
    // there are infinitely many solutions.
    isUnique: boolean;

    // The rectangle origin (u,v).
    rectOrigin: Vector;

    // The rectangle width w.
    rectWidth: number;

    // The rectangle height h = w / aspectRatio.
    rectHeight: number;
}

// A plane of the form Dot(normal,X) = constant.
interface Constraint {
    normal: Vector;
    constant: number;
}

interface LineFromPlanes {
    isLine: boolean;
    origin: Vector;
    direction: Vector;
}

export class InscribedFixedAspectRectInQuad {
    // The caller is responsible for the 'quad' vertices occurring in
    // counterclockwise order. The output 'rectOrigin' is (u,v), the
    // 'rectWidth' is w, and the 'rectHeight' is h. The rectangle vertices are
    // (u,v), (u + w, v), (u + w, v + h), and (u, v + h) in counterclockwise
    // order.
    static execute(quad: readonly Vector[], aspectRatio: number):
        InscribedFixedAspectRectInQuadResult {
        logAssert(quad.length === 4, 'The quadrilateral must have 4 vertices.');

        let isUnique = false;

        // The i-th edge lies on a line with origin quad[i] and non-unit
        // direction edges[i] = quad[(i + 1) % 4] - quad[i]. The lines
        // containing the edges have these inner-pointing normal vectors.
        const normals: Vector[] = [
            perp(sub(quad[0], quad[1])),
            perp(sub(quad[1], quad[2])),
            perp(sub(quad[2], quad[3])),
            perp(sub(quad[3], quad[0]))
        ];

        // Compute the 4 linear inequality constraints of the form
        // Dot(N[i], R[floor(2*angle[i]/pi)] - V[i]) >= 0, where V[i] is a
        // quad vertex and N[i] is a corresponding normal. The angle[i] is the
        // angle formed by N[i] with the positive x-axis and is in [0,2*pi).
        // Each constraint is written as Dot((c0,c1,c2),(u,v,w)) + c3 >= 0. In
        // the comments for rect[], r is the aspect ratio w/h.
        const zero = 0;
        const twoPi = GTE_C_TWO_PI;
        const invTwoDivPi = GTE_C_INV_HALF_PI;

        const constraints: Constraint[] = [];
        for (let i = 0; i < 4; ++i) {
            let angle = Math.atan2(normals[i].values[1], normals[i].values[0]);
            if (angle < zero) {
                angle += twoPi;
            }

            const j = Math.floor(invTwoDivPi * angle);
            const normal = Vector.fromArray([normals[i].values[0], normals[i].values[1], 0]);
            const constant = dot(normals[i], quad[i]);
            if (j === 0) {
                // rect[0] = (u, v)
                normal.values[2] = zero;
            }
            else if (j === 1) {
                // rect[1] = (u, v) + (w, 0)
                normal.values[2] = normals[i].values[0];
            }
            else if (j === 2) {
                // rect[2] = (u, v) + (w, w / r)
                normal.values[2] = normals[i].values[0] + normals[i].values[1] / aspectRatio;
            }
            else { // j === 3
                // rect[3] = (u, v) + (0, w / r)
                normal.values[2] = normals[i].values[1] / aspectRatio;
            }
            constraints.push({ normal, constant });
        }

        // Solve the linear programming problem.
        const iiQuery = new IntrIntervalsFI();
        let line = InscribedFixedAspectRectInQuad.findIntersection(
            constraints[0].normal, constraints[0].constant,
            constraints[2].normal, constraints[2].constant);
        logAssert(line.isLine, 'Unexpected condition.');

        const alpha1 = dot(constraints[1].normal, line.direction);
        const beta1 = dot(constraints[1].normal, line.origin) - constraints[1].constant;
        const alpha3 = dot(constraints[3].normal, line.direction);
        const beta3 = dot(constraints[3].normal, line.origin) - constraints[3].constant;
        logAssert(alpha1 !== zero && alpha3 !== zero, 'Unexpected condition.');

        const end1 = -beta1 / alpha1;
        const isPositiveInfinite1 = (alpha1 > zero);
        const end3 = -beta3 / alpha3;
        const isPositiveInfinite3 = (alpha3 > zero);
        let iiResult = iiQuery.findSemiInfiniteSemiInfinite(
            end1, isPositiveInfinite1, end3, isPositiveInfinite3);

        let rectOrigin = Vector.zero(2);
        let rectWidth = 0;
        let rectHeight = 0;

        // Select the solution with the larger w. The rectangle height is
        // determined by the width and the aspect ratio.
        const assignFromSolution = (solution: Vector): void => {
            rectOrigin = Vector.fromArray([solution.values[0], solution.values[1]]);
            rectWidth = solution.values[2];
            rectHeight = rectWidth / aspectRatio;
        };

        if (iiResult.type === IntrIntervalsFIResultType.isFinite) {
            const solution0 = add(mul(line.direction, iiResult.overlap[0]), line.origin);
            const solution1 = add(mul(line.direction, iiResult.overlap[1]), line.origin);
            if (solution0.values[2] > solution1.values[2]) {
                assignFromSolution(solution0);
            }
            else {
                assignFromSolution(solution1);
            }

            isUnique = (solution0.values[2] !== solution1.values[2]);
        }
        else if (iiResult.type === IntrIntervalsFIResultType.isPoint) {
            assignFromSolution(add(mul(line.direction, iiResult.overlap[0]), line.origin));
            isUnique = true;
        }
        else if (iiResult.type === IntrIntervalsFIResultType.isEmpty) {
            line = InscribedFixedAspectRectInQuad.findIntersection(
                constraints[1].normal, constraints[1].constant,
                constraints[3].normal, constraints[3].constant);
            logAssert(line.isLine, 'Unexpected condition.');

            const alpha0 = dot(constraints[0].normal, line.direction);
            const beta0 = dot(constraints[0].normal, line.origin) - constraints[0].constant;
            const alpha2 = dot(constraints[2].normal, line.direction);
            const beta2 = dot(constraints[2].normal, line.origin) - constraints[2].constant;
            logAssert(alpha0 !== zero && alpha2 !== zero, 'Unexpected condition.');

            const end0 = -beta0 / alpha0;
            const isPositiveInfinite0 = (alpha0 > zero);
            const end2 = -beta2 / alpha2;
            const isPositiveInfinite2 = (alpha2 > zero);
            iiResult = iiQuery.findSemiInfiniteSemiInfinite(
                end0, isPositiveInfinite0, end2, isPositiveInfinite2);
            if (iiResult.type === IntrIntervalsFIResultType.isFinite) {
                const solution0 = add(mul(line.direction, iiResult.overlap[0]), line.origin);
                const solution1 = add(mul(line.direction, iiResult.overlap[1]), line.origin);
                if (solution0.values[2] > solution1.values[2]) {
                    assignFromSolution(solution0);
                }
                else {
                    assignFromSolution(solution1);
                }

                isUnique = (solution0.values[2] !== solution1.values[2]);
            }
            else if (iiResult.type === IntrIntervalsFIResultType.isPoint) {
                assignFromSolution(add(mul(line.direction, iiResult.overlap[0]), line.origin));
                isUnique = true;
            }
            else {
                logError('Unexpected interval intersection type.');
            }
        }
        else {
            logError('Unexpected interval intersection type.');
        }

        return { isUnique, rectOrigin, rectWidth, rectHeight };
    }

    private static findIntersection(normal0: Vector, constant0: number,
        normal1: Vector, constant1: number): LineFromPlanes {
        // The intersection line is of the form
        // t * Cross(normal0, normal1) + a0 * normal0 + a1 * normal1
        const direction = cross(normal0, normal1);
        if (direction.values[0] !== 0 || direction.values[1] !== 0
            || direction.values[2] !== 0) {
            const dotN0N0 = dot(normal0, normal0);
            const dotN0N1 = dot(normal0, normal1);
            const dotN1N1 = dot(normal1, normal1);
            const det = dot(direction, direction);
            const a0 = (dotN1N1 * constant0 - dotN0N1 * constant1) / det;
            const a1 = (dotN0N0 * constant1 - dotN0N1 * constant0) / det;
            const origin = add(mul(normal0, a0), mul(normal1, a1));
            return { isLine: true, origin, direction };
        }
        else {
            return { isLine: false, origin: Vector.zero(3), direction: Vector.zero(3) };
        }
    }
}
