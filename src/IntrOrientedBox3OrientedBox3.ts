// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box to be a solid.
//
// The test-intersection query uses the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The set of potential separating directions includes the 3 face normals of
// box0, the 3 face normals of box1, and 9 directions, each of which is the
// cross product of an edge of box0 and an edge of box1.
//
// The separating axes involving cross products of edges has numerical
// robustness problems when the two edges are nearly parallel. The cross
// product of the edges is nearly the zero vector, so normalization of the
// cross product may produce unit-length directions that are not close to the
// true direction. Such a pair of edges occurs when an object0 face normal N0
// and an object1 face normal N1 are nearly parallel. In this case, you may
// skip the edge-edge directions, which is equivalent to projecting the
// objects onto the plane with normal N0 and applying a 2D separating axis
// test. The ability to do so involves choosing a small nonnegative epsilon.
// It is used to determine whether two face normals, one from each object, are
// nearly parallel: |Dot(N0,N1)| >= 1 - epsilon. If the epsilon input to the
// test(...) function is negative, it is clamped to zero.
//
// The pair of integers 'separating', say, (i0,i1), identifies the axes that
// reported separation; there may be more than one but only one is reported.
// If the separating axis is a face normal N[i0] of object0, then (i0,-1) is
// returned. If the axis is a face normal N[i1], then (-1,i1) is returned. If
// the axis is a cross product of edges, Cross(N[i0],N[i1]), then (i0,i1) is
// returned. If 'intersect' is true, the separating[] values are invalid
// because there is no separation.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): upstream has
// only the TIQuery specialization for this pair, which becomes the class
// IntrOrientedBox3OrientedBox3TI with the result type
// IntrOrientedBox3OrientedBox3TIResult.

import type { TIQuery } from './TIQuery';
import type { OrientedBox } from './OrientedBox';
import { dot, sub } from './Vector';

// The result of IntrOrientedBox3OrientedBox3TI queries.
export interface IntrOrientedBox3OrientedBox3TIResult {
    intersect: boolean;

    // The pair of separating axis indices, valid only when 'intersect' is
    // false. See the file comments for the interpretation.
    separating: [number, number];
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrOrientedBox3OrientedBox3TIResult {
    return { intersect: false, separating: [0, 0] };
}

// Test-intersection query for two solid oriented boxes in 3D using the
// 15-axis separating-axis test.
export class IntrOrientedBox3OrientedBox3TI implements TIQuery<OrientedBox, OrientedBox, IntrOrientedBox3OrientedBox3TIResult> {
    test(box0: OrientedBox, box1: OrientedBox,
        epsilon: number = 0): IntrOrientedBox3OrientedBox3TIResult {
        const result = defaultTIResult();

        // Convenience variables.
        const C0 = box0.center;
        const A0 = box0.axis;
        const E0 = box0.extent.values;
        const C1 = box1.center;
        const A1 = box1.axis;
        const E1 = box1.extent.values;

        epsilon = Math.max(epsilon, 0);
        const cutoff = 1 - epsilon;
        let existsParallelPair = false;

        // Compute difference of box centers.
        const D = sub(C1, C0);

        // dot01[i][j] = Dot(A0[i],A1[j])
        const dot01 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

        // |dot01[i][j]|
        const absDot01 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

        // Dot(D, A0[i])
        const dotDA0 = [0, 0, 0];

        // interval radii and distance between centers
        let r0: number, r1: number, r: number;

        // r0 + r1
        let r01: number;

        // Test for separation on the axis C0 + t*A0[0].
        for (let i = 0; i < 3; ++i) {
            dot01[0][i] = dot(A0[0], A1[i]);
            absDot01[0][i] = Math.abs(dot01[0][i]);
            if (absDot01[0][i] > cutoff) {
                existsParallelPair = true;
            }
        }
        dotDA0[0] = dot(D, A0[0]);
        r = Math.abs(dotDA0[0]);
        r1 = E1[0] * absDot01[0][0] + E1[1] * absDot01[0][1] + E1[2] * absDot01[0][2];
        r01 = E0[0] + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 0;
            result.separating[1] = -1;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[1].
        for (let i = 0; i < 3; ++i) {
            dot01[1][i] = dot(A0[1], A1[i]);
            absDot01[1][i] = Math.abs(dot01[1][i]);
            if (absDot01[1][i] > cutoff) {
                existsParallelPair = true;
            }
        }
        dotDA0[1] = dot(D, A0[1]);
        r = Math.abs(dotDA0[1]);
        r1 = E1[0] * absDot01[1][0] + E1[1] * absDot01[1][1] + E1[2] * absDot01[1][2];
        r01 = E0[1] + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 1;
            result.separating[1] = -1;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[2].
        for (let i = 0; i < 3; ++i) {
            dot01[2][i] = dot(A0[2], A1[i]);
            absDot01[2][i] = Math.abs(dot01[2][i]);
            if (absDot01[2][i] > cutoff) {
                existsParallelPair = true;
            }
        }
        dotDA0[2] = dot(D, A0[2]);
        r = Math.abs(dotDA0[2]);
        r1 = E1[0] * absDot01[2][0] + E1[1] * absDot01[2][1] + E1[2] * absDot01[2][2];
        r01 = E0[2] + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 2;
            result.separating[1] = -1;
            return result;
        }

        // Test for separation on the axis C0 + t*A1[0].
        r = Math.abs(dot(D, A1[0]));
        r0 = E0[0] * absDot01[0][0] + E0[1] * absDot01[1][0] + E0[2] * absDot01[2][0];
        r01 = r0 + E1[0];
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = -1;
            result.separating[1] = 0;
            return result;
        }

        // Test for separation on the axis C0 + t*A1[1].
        r = Math.abs(dot(D, A1[1]));
        r0 = E0[0] * absDot01[0][1] + E0[1] * absDot01[1][1] + E0[2] * absDot01[2][1];
        r01 = r0 + E1[1];
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = -1;
            result.separating[1] = 1;
            return result;
        }

        // Test for separation on the axis C0 + t*A1[2].
        r = Math.abs(dot(D, A1[2]));
        r0 = E0[0] * absDot01[0][2] + E0[1] * absDot01[1][2] + E0[2] * absDot01[2][2];
        r01 = r0 + E1[2];
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = -1;
            result.separating[1] = 2;
            return result;
        }

        // At least one pair of box axes was parallel, so the separation is
        // effectively in 2D. The edge-edge axes do not need to be tested.
        if (existsParallelPair) {
            // The result.separating[] values are invalid because there is no
            // separation.
            result.intersect = true;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[0]xA1[0].
        r = Math.abs(dotDA0[2] * dot01[1][0] - dotDA0[1] * dot01[2][0]);
        r0 = E0[1] * absDot01[2][0] + E0[2] * absDot01[1][0];
        r1 = E1[1] * absDot01[0][2] + E1[2] * absDot01[0][1];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 0;
            result.separating[1] = 0;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[0]xA1[1].
        r = Math.abs(dotDA0[2] * dot01[1][1] - dotDA0[1] * dot01[2][1]);
        r0 = E0[1] * absDot01[2][1] + E0[2] * absDot01[1][1];
        r1 = E1[0] * absDot01[0][2] + E1[2] * absDot01[0][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 0;
            result.separating[1] = 1;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[0]xA1[2].
        r = Math.abs(dotDA0[2] * dot01[1][2] - dotDA0[1] * dot01[2][2]);
        r0 = E0[1] * absDot01[2][2] + E0[2] * absDot01[1][2];
        r1 = E1[0] * absDot01[0][1] + E1[1] * absDot01[0][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 0;
            result.separating[1] = 2;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[1]xA1[0].
        r = Math.abs(dotDA0[0] * dot01[2][0] - dotDA0[2] * dot01[0][0]);
        r0 = E0[0] * absDot01[2][0] + E0[2] * absDot01[0][0];
        r1 = E1[1] * absDot01[1][2] + E1[2] * absDot01[1][1];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 1;
            result.separating[1] = 0;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[1]xA1[1].
        r = Math.abs(dotDA0[0] * dot01[2][1] - dotDA0[2] * dot01[0][1]);
        r0 = E0[0] * absDot01[2][1] + E0[2] * absDot01[0][1];
        r1 = E1[0] * absDot01[1][2] + E1[2] * absDot01[1][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 1;
            result.separating[1] = 1;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[1]xA1[2].
        r = Math.abs(dotDA0[0] * dot01[2][2] - dotDA0[2] * dot01[0][2]);
        r0 = E0[0] * absDot01[2][2] + E0[2] * absDot01[0][2];
        r1 = E1[0] * absDot01[1][1] + E1[1] * absDot01[1][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 1;
            result.separating[1] = 2;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[2]xA1[0].
        r = Math.abs(dotDA0[1] * dot01[0][0] - dotDA0[0] * dot01[1][0]);
        r0 = E0[0] * absDot01[1][0] + E0[1] * absDot01[0][0];
        r1 = E1[1] * absDot01[2][2] + E1[2] * absDot01[2][1];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 2;
            result.separating[1] = 0;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[2]xA1[1].
        r = Math.abs(dotDA0[1] * dot01[0][1] - dotDA0[0] * dot01[1][1]);
        r0 = E0[0] * absDot01[1][1] + E0[1] * absDot01[0][1];
        r1 = E1[0] * absDot01[2][2] + E1[2] * absDot01[2][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 2;
            result.separating[1] = 1;
            return result;
        }

        // Test for separation on the axis C0 + t*A0[2]xA1[2].
        r = Math.abs(dotDA0[1] * dot01[0][2] - dotDA0[0] * dot01[1][2]);
        r0 = E0[0] * absDot01[1][2] + E0[1] * absDot01[0][2];
        r1 = E1[0] * absDot01[2][1] + E1[1] * absDot01[2][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating[0] = 2;
            result.separating[1] = 2;
            return result;
        }

        // The result.separating[] values are invalid because there is no
        // separation.
        result.intersect = true;
        return result;
    }
}
