// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrAlignedBox3OrientedBox3.h
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
// true direction. Such a pair of edges occurs when a box0 face normal N0 and
// a box1 face normal N1 are nearly parallel. In this case, you may skip the
// edge-edge directions, which is equivalent to projecting the boxes onto the
// plane with normal N0 and applying a 2D separating axis test. The ability
// to do so involves choosing a small nonnegative epsilon. It is used to
// determine whether two face normals, one from each box, are nearly
// parallel: |Dot(N0,N1)| >= 1 - epsilon. If the epsilon input to test(...)
// is negative, it is clamped to zero.
//
// The pair of integers 'separating', say, (i0,i1), identifies the axis that
// reported separation; there may be more than one but only one is reported.
// If the separating axis is a face normal N[i0] of the aligned box0 in
// dimension i0, then (i0,-1) is returned. If the axis is a face normal
// box1.axis[i1], then (-1,i1) is returned. If the axis is a cross product of
// edges, Cross(N[i0],box1.axis[i1]), then (i0,i1) is returned. If 'intersect'
// is true, the separating[] values are invalid because there is no
// separation.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only a TIQuery specialization for this pair of primitives. The optional
// epsilon argument of 'operator()' is an optional third argument of test().

import { AlignedBox } from './AlignedBox';
import { OrientedBox } from './OrientedBox';
import { sub, dot } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrAlignedBox3OrientedBox3TI.test.
export interface IntrAlignedBox3OrientedBox3TIResult {
    intersect: boolean;

    // Valid only when 'intersect' is false. See the file comments for the
    // meaning of the pair.
    separating: [number, number];
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrAlignedBox3OrientedBox3TIResult {
    return { intersect: false, separating: [0, 0] };
}

export class IntrAlignedBox3OrientedBox3TI implements
    TIQuery<AlignedBox, OrientedBox, IntrAlignedBox3OrientedBox3TIResult> {

    test(box0: AlignedBox, box1: OrientedBox, epsilon: number = 0):
        IntrAlignedBox3OrientedBox3TIResult {
        const result = defaultTIResult();

        // Get the centered form of the aligned box. The axes are implicitly
        // A0[0] = (1,0,0), A0[1] = (0,1,0) and A0[2] = (0,0,1).
        const { center: C0, extent: E0 } = box0.getCenteredForm();

        // Convenience variables.
        const C1 = box1.center;
        const A1 = box1.axis;
        const E1 = box1.extent;
        const e0 = E0.values;
        const e1 = E1.values;

        epsilon = Math.max(epsilon, 0);
        const cutoff = 1 - epsilon;
        let existsParallelPair = false;

        // Compute the difference of box centers.
        const Dvec = sub(C1, C0);
        const D = Dvec.values;

        // dot01[i][j] = Dot(A0[i],A1[j]) = A1[j][i]
        const dot01: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

        // |dot01[i][j]|
        const absDot01: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

        // interval radii and distance between centers
        let r0 = 0, r1 = 0, r = 0;

        // r0 + r1
        let r01 = 0;

        // Test for separation on the axis C0 + t*A0[0].
        for (let i = 0; i < 3; ++i) {
            dot01[0][i] = A1[i].values[0];
            absDot01[0][i] = Math.abs(A1[i].values[0]);
            if (absDot01[0][i] >= cutoff) {
                existsParallelPair = true;
            }
        }
        r = Math.abs(D[0]);
        r1 = e1[0] * absDot01[0][0] + e1[1] * absDot01[0][1] + e1[2] * absDot01[0][2];
        r01 = e0[0] + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [0, -1];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[1].
        for (let i = 0; i < 3; ++i) {
            dot01[1][i] = A1[i].values[1];
            absDot01[1][i] = Math.abs(A1[i].values[1]);
            if (absDot01[1][i] >= cutoff) {
                existsParallelPair = true;
            }
        }
        r = Math.abs(D[1]);
        r1 = e1[0] * absDot01[1][0] + e1[1] * absDot01[1][1] + e1[2] * absDot01[1][2];
        r01 = e0[1] + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [1, -1];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[2].
        for (let i = 0; i < 3; ++i) {
            dot01[2][i] = A1[i].values[2];
            absDot01[2][i] = Math.abs(A1[i].values[2]);
            if (absDot01[2][i] >= cutoff) {
                existsParallelPair = true;
            }
        }
        r = Math.abs(D[2]);
        r1 = e1[0] * absDot01[2][0] + e1[1] * absDot01[2][1] + e1[2] * absDot01[2][2];
        r01 = e0[2] + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [2, -1];
            return result;
        }

        // Test for separation on the axis C0 + t*A1[0].
        r = Math.abs(dot(Dvec, A1[0]));
        r0 = e0[0] * absDot01[0][0] + e0[1] * absDot01[1][0] + e0[2] * absDot01[2][0];
        r01 = r0 + e1[0];
        if (r > r01) {
            result.intersect = false;
            result.separating = [-1, 0];
            return result;
        }

        // Test for separation on the axis C0 + t*A1[1].
        r = Math.abs(dot(Dvec, A1[1]));
        r0 = e0[0] * absDot01[0][1] + e0[1] * absDot01[1][1] + e0[2] * absDot01[2][1];
        r01 = r0 + e1[1];
        if (r > r01) {
            result.intersect = false;
            result.separating = [-1, 1];
            return result;
        }

        // Test for separation on the axis C0 + t*A1[2].
        r = Math.abs(dot(Dvec, A1[2]));
        r0 = e0[0] * absDot01[0][2] + e0[1] * absDot01[1][2] + e0[2] * absDot01[2][2];
        r01 = r0 + e1[2];
        if (r > r01) {
            result.intersect = false;
            result.separating = [-1, 2];
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
        r = Math.abs(D[2] * dot01[1][0] - D[1] * dot01[2][0]);
        r0 = e0[1] * absDot01[2][0] + e0[2] * absDot01[1][0];
        r1 = e1[1] * absDot01[0][2] + e1[2] * absDot01[0][1];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [0, 0];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[0]xA1[1].
        r = Math.abs(D[2] * dot01[1][1] - D[1] * dot01[2][1]);
        r0 = e0[1] * absDot01[2][1] + e0[2] * absDot01[1][1];
        r1 = e1[0] * absDot01[0][2] + e1[2] * absDot01[0][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [0, 1];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[0]xA1[2].
        r = Math.abs(D[2] * dot01[1][2] - D[1] * dot01[2][2]);
        r0 = e0[1] * absDot01[2][2] + e0[2] * absDot01[1][2];
        r1 = e1[0] * absDot01[0][1] + e1[1] * absDot01[0][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [0, 2];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[1]xA1[0].
        r = Math.abs(D[0] * dot01[2][0] - D[2] * dot01[0][0]);
        r0 = e0[0] * absDot01[2][0] + e0[2] * absDot01[0][0];
        r1 = e1[1] * absDot01[1][2] + e1[2] * absDot01[1][1];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [1, 0];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[1]xA1[1].
        r = Math.abs(D[0] * dot01[2][1] - D[2] * dot01[0][1]);
        r0 = e0[0] * absDot01[2][1] + e0[2] * absDot01[0][1];
        r1 = e1[0] * absDot01[1][2] + e1[2] * absDot01[1][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [1, 1];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[1]xA1[2].
        r = Math.abs(D[0] * dot01[2][2] - D[2] * dot01[0][2]);
        r0 = e0[0] * absDot01[2][2] + e0[2] * absDot01[0][2];
        r1 = e1[0] * absDot01[1][1] + e1[1] * absDot01[1][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [1, 2];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[2]xA1[0].
        r = Math.abs(D[1] * dot01[0][0] - D[0] * dot01[1][0]);
        r0 = e0[0] * absDot01[1][0] + e0[1] * absDot01[0][0];
        r1 = e1[1] * absDot01[2][2] + e1[2] * absDot01[2][1];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [2, 0];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[2]xA1[1].
        r = Math.abs(D[1] * dot01[0][1] - D[0] * dot01[1][1]);
        r0 = e0[0] * absDot01[1][1] + e0[1] * absDot01[0][1];
        r1 = e1[0] * absDot01[2][2] + e1[2] * absDot01[2][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [2, 1];
            return result;
        }

        // Test for separation on the axis C0 + t*A0[2]xA1[2].
        r = Math.abs(D[1] * dot01[0][2] - D[0] * dot01[1][2]);
        r0 = e0[0] * absDot01[1][2] + e0[1] * absDot01[0][2];
        r1 = e1[0] * absDot01[2][1] + e1[1] * absDot01[2][0];
        r01 = r0 + r1;
        if (r > r01) {
            result.intersect = false;
            result.separating = [2, 2];
            return result;
        }

        // The result.separating[] values are invalid because there is no
        // separation.
        result.intersect = true;
        return result;
    }
}
