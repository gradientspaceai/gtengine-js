// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrCanonicalBox3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query is for finite cylinders. The cylinder and box are considered to
// be solids. The cylinder has center C, unit-length axis direction D, radius
// r and height h. The canonical box has center at the origin and extents E.
//
// The abstract algorithm clips the canonical box by the planes of the
// cylinder end disks to obtain a convex polyhedron Q. This polyhedron is
// projected to a convex polygon P in the plane Dot(D, X - C) = 0. The
// cylinder axis projects to C. The box and cylinder intersect when
// Distance(C,P) <= r. If C is inside or on P, the distance is 0 and there is
// an intersection. If C is outside P, the distance is the minimum of the
// distances from C to the edges of P.
//
// The implementation is described in
//   https://www.geometrictools.com/Documentation/IntersectionBoxCylinder.pdf
// and avoids using a generic convex hull algorithm. Reductions in dimension
// occur based on the number of 0-valued components of the cylinder axis
// direction.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrCanonicalBox3Cylinder3TI.
// The private static helpers become module-private functions. The upstream
// private static CylinderAxisIntersectsBox3D is dead code (never called; the
// same computation is inlined at the top of DoQueryNoZeros) and is omitted.

import type { CanonicalBox } from './CanonicalBox.js';
import type { Cylinder3 } from './Cylinder3.js';
import { logAssert, logError } from './Logger.js';
import { Vector, dot } from './Vector.js';
import { computeOrthogonalComplement3 } from './Vector3.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrCanonicalBox3Cylinder3TI.test.
export interface IntrCanonicalBox3Cylinder3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrCanonicalBox3Cylinder3TIResult {
    return { intersect: false };
}

// Test whether the box is outside the slab contained by the planes of the
// cylinder end disks. This is accomplished by computing the interval of
// projection of the box onto the cylinder axis.
function boxIsOutsideCylinderSlab(box: CanonicalBox,
    cylinder: Cylinder3): boolean {
    const C = cylinder.axis.origin;
    const D = cylinder.axis.direction;
    const absD = Vector.fromArray([
        Math.abs(D.values[0]),
        Math.abs(D.values[1]),
        Math.abs(D.values[2])
    ]);
    const hDiv2 = 0.5 * cylinder.height;
    const E = box.extent;

    // Let g be the interval center, p be the interval radius and h be the
    // cylinder height. The culling test is g - p > h/2 (box above the slab)
    // or g + p < -h/2 (box below the slab). The tests can be rewritten as
    // g > p + h/2 or -g > p + h/2. In turn these are combined to
    // |g| > p + h/2.
    const intervalCenter = -dot(D, C);  // Dot(D, boxCenter-cylCenter)
    const intervalRadius = dot(E, absD);
    return Math.abs(intervalCenter) > intervalRadius + hDiv2;
}

// Compute the squared distance from (0,0) to the projection of the segment
// <P0,P1>. The projection plane has origin C and is spanned by the
// orthonormal vectors W0 and W1.
function computeSqrDistance(P0: readonly number[], P1: readonly number[],
    C: readonly number[], W0: readonly number[],
    W1: readonly number[]): number {
    const p0mc = [P0[0] - C[0], P0[1] - C[1], P0[2] - C[2]];
    const p1mc = [P1[0] - C[0], P1[1] - C[1], P1[2] - C[2]];
    const q0x = W0[0] * p0mc[0] + W0[1] * p0mc[1] + W0[2] * p0mc[2];
    const q0y = W1[0] * p0mc[0] + W1[1] * p0mc[1] + W1[2] * p0mc[2];
    const q1x = W0[0] * p1mc[0] + W0[1] * p1mc[1] + W0[2] * p1mc[2];
    const q1y = W1[0] * p1mc[0] + W1[1] * p1mc[1] + W1[2] * p1mc[2];

    const dirX = q1x - q0x, dirY = q1y - q0y;
    let s = dirX * q1x + dirY * q1y;
    if (s <= 0) {
        return q1x * q1x + q1y * q1y;
    } else {
        s = dirX * q0x + dirY * q0y;
        if (s >= 0) {
            return q0x * q0x + q0y * q0y;
        } else {
            s /= dirX * dirX + dirY * dirY;
            const cx = q0x - s * dirX, cy = q0y - s * dirY;
            return cx * cx + cy * cy;
        }
    }
}

function doQueryTwoZeros(i: readonly number[], C: readonly number[], r: number,
    E: readonly number[]): boolean {
    // The 2-tuple (C[i[1]], C[i[2]]) is the projected cylinder axis. The
    // 2-tuple (E[i[1]], E[i[2]]) is the extent of the projected canonical
    // box, which is an axis-aligned rectangle.
    const absC1 = Math.abs(C[i[1]]), absC2 = Math.abs(C[i[2]]);
    const E1 = E[i[1]], E2 = E[i[2]];

    // Test whether the cylinder axis and canonical box intersect.
    if (absC1 <= E1 && absC2 <= E2) {
        return true;
    }

    // Compute the squared distance from the projected cylinder axis to the
    // projected canonical box.
    let sqrDistance = 0;
    let delta = absC1 - E1;
    if (delta > 0) {
        sqrDistance += delta * delta;
    }

    delta = absC2 - E2;
    if (delta > 0) {
        sqrDistance += delta * delta;
    }

    return sqrDistance <= r * r;
}

function doQueryOneZero(i: readonly number[], C: readonly number[],
    D: readonly number[], r: number, hDiv2: number,
    E: readonly number[]): boolean {
    const c0 = C[i[0]], c1 = C[i[1]], c2 = C[i[2]];
    const d0 = D[i[0]], d1 = D[i[1]];
    const e0 = E[i[0]], e1 = E[i[1]], e2 = E[i[2]];
    const e0pc0 = e0 + c0, e0mc0 = e0 - c0;
    const e1pc1 = e1 + c1, e1mc1 = e1 - c1;

    // Test whether the cylinder axis and canonical box intersect.
    const absC2 = Math.abs(c2);
    if (absC2 <= e2) {
        const negEmCDivD = [-e0pc0 / d0, -e1pc1 / d1];
        const posEmCDivD = [e0mc0 / d0, e1mc1 / d1];
        const lower = Math.max(Math.max(negEmCDivD[0], negEmCDivD[1]), -hDiv2);
        const upper = Math.min(Math.min(posEmCDivD[0], posEmCDivD[1]), hDiv2);
        if (lower <= upper) {
            return true;
        }
    }

    // Compute the squared distance from the projected cylinder axis (a point)
    // to the projected convex polyhedron (a rectangle).
    let sMin = 0;
    const tHat = d1 * e1mc1 - d0 * e0pc0;
    if (-hDiv2 <= tHat) {
        if (tHat <= hDiv2) {
            sMin = -(d0 * e1mc1 + d1 * e0pc0);
        } else {
            // tHat > +h/2
            sMin = -(e0pc0 + d0 * hDiv2) / d1;
        }
    } else {
        // tHat < -h/2
        sMin = -(e1mc1 + d1 * hDiv2) / d0;
    }

    let sMax = 0;
    const tBar = d0 * e0mc0 - d1 * e1pc1;
    if (-hDiv2 <= tBar) {
        if (tBar <= hDiv2) {
            sMax = d0 * e1pc1 + d1 * e0mc0;
        } else {
            // tBar > +h/2
            sMax = (e1pc1 + d1 * hDiv2) / d0;
        }
    } else {
        // tBar < -h/2
        sMax = (e0mc0 + d0 * hDiv2) / d1;
    }

    logAssert(sMin < sMax, 'The s-interval is invalid, which is unexpected.');

    let sqrDistance = 0;
    if (0 < sMin) {
        sqrDistance += sMin * sMin;
    } else if (sMax < 0) {
        sqrDistance += sMax * sMax;
    }

    const delta = absC2 - e2;
    if (delta > 0) {
        sqrDistance += delta * delta;
    }

    return sqrDistance <= r * r;
}

function doQueryNoZeros(C: readonly number[], D: readonly number[], r: number,
    hDiv2: number, E: readonly number[]): boolean {
    // Test whether the cylinder axis and canonical box intersect.
    const negEmCDivD = [
        (-E[0] - C[0]) / D[0],
        (-E[1] - C[1]) / D[1],
        (-E[2] - C[2]) / D[2]
    ];

    const posEmCDivD = [
        (E[0] - C[0]) / D[0],
        (E[1] - C[1]) / D[1],
        (E[2] - C[2]) / D[2]
    ];

    const max01 = Math.max(negEmCDivD[0], negEmCDivD[1]);
    const max23 = Math.max(negEmCDivD[2], -hDiv2);
    let lower = Math.max(max01, max23);
    const min01 = Math.min(posEmCDivD[0], posEmCDivD[1]);
    const min23 = Math.min(posEmCDivD[2], hDiv2);
    let upper = Math.min(min01, min23);
    if (lower <= upper) {
        return true;
    }

    // Compute t[i] = Dot(D, V[i] - C) for box vertices V[i]. These are used
    // in computing the intervals associated with extreme edges.
    const dotDC = D[0] * C[0] + D[1] * C[1] + D[2] * C[2];
    const d0e0 = D[0] * E[0], d1e1 = D[1] * E[1], d2e2 = D[2] * E[2];
    const t1 = +d0e0 - d1e1 - d2e2 - dotDC, s1p = t1 + hDiv2, s1n = t1 - hDiv2;
    const t2 = -d0e0 + d1e1 - d2e2 - dotDC, s2p = t2 + hDiv2, s2n = t2 - hDiv2;
    const t3 = +d0e0 + d1e1 - d2e2 - dotDC, s3p = t3 + hDiv2, s3n = t3 - hDiv2;
    const t4 = -d0e0 - d1e1 + d2e2 - dotDC, s4p = t4 + hDiv2, s4n = t4 - hDiv2;
    const t5 = +d0e0 - d1e1 + d2e2 - dotDC, s5p = t5 + hDiv2, s5n = t5 - hDiv2;
    const t6 = -d0e0 + d1e1 + d2e2 - dotDC, s6p = t6 + hDiv2, s6n = t6 - hDiv2;

    // Compute an orthonormal basis containing D.
    const basis: Vector[] = [Vector.fromArray([D[0], D[1], D[2]]),
        new Vector(3), new Vector(3)];
    computeOrthogonalComplement3(1, basis);
    const W0 = basis[1].values;
    const W1 = basis[2].values;

    const sqrRadius = r * r;
    let sqrDistance = 0;
    let P0: number[] = [0, 0, 0];
    let P1: number[] = [0, 0, 0];

    // (U0, -U1)
    lower = (s1p >= 0 ? -E[2] : -E[2] - s1p / D[2]);
    upper = (s5n <= 0 ? +E[2] : +E[2] - s5n / D[2]);
    if (lower <= upper) {
        P0 = [+E[0], -E[1], lower];
        P1 = [+E[0], -E[1], upper];
        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (U1, -U0)
    lower = (s2p >= 0 ? -E[2] : -E[2] - s2p / D[2]);
    upper = (s6n <= 0 ? +E[2] : +E[2] - s6n / D[2]);
    if (lower <= upper) {
        P0 = [-E[0], +E[1], lower];
        P1 = [-E[0], +E[1], upper];
        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (U0, -U2)
    lower = (s1p >= 0 ? -E[1] : -E[1] - s1p / D[1]);
    upper = (s3n <= 0 ? +E[1] : +E[1] - s3n / D[1]);
    if (lower <= upper) {
        P0 = [+E[0], lower, -E[2]];
        P1 = [+E[0], upper, -E[2]];
        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (U2, -U0)
    lower = (s4p >= 0 ? -E[1] : -E[1] - s4p / D[1]);
    upper = (s6n <= 0 ? +E[1] : +E[1] - s6n / D[1]);
    if (lower <= upper) {
        P0 = [-E[0], lower, +E[2]];
        P1 = [-E[0], upper, +E[2]];
        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (U1, -U2)
    lower = (s2p >= 0 ? -E[0] : -E[0] - s2p / D[0]);
    upper = (s3n <= 0 ? +E[0] : +E[0] - s3n / D[0]);
    if (lower <= upper) {
        P0 = [lower, +E[1], -E[2]];
        P1 = [upper, +E[1], -E[2]];
        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (U2, -U1)
    lower = (s4p >= 0 ? -E[0] : -E[0] - s4p / D[0]);
    upper = (s5n <= 0 ? +E[0] : +E[0] - s5n / D[0]);
    if (lower <= upper) {
        P0 = [lower, -E[1], +E[2]];
        P1 = [upper, -E[1], +E[2]];
        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (U0, -D)
    lower = (s3p >= 0 ? -E[2] : -E[2] - s3p / D[2]);
    upper = (s5p <= 0 ? +E[2] : +E[2] - s5p / D[2]);
    if (lower <= upper) {
        if (s3p >= 0) {
            P0 = [+E[0], +E[1] - s3p / D[1], -E[2]];
        } else {
            P0 = [+E[0], +E[1], -E[2] - s3p / D[2]];
        }

        if (s5p <= 0) {
            P1 = [+E[0], -E[1] - s5p / D[1], +E[2]];
        } else {
            P1 = [+E[0], -E[1], +E[2] - s5p / D[2]];
        }

        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (D, -U0)
    lower = (s2n >= 0 ? -E[2] : -E[2] - s2n / D[2]);
    upper = (s4n <= 0 ? +E[2] : +E[2] - s4n / D[2]);
    if (lower <= upper) {
        if (s2n >= 0) {
            P0 = [-E[0], +E[1] - s2n / D[1], -E[2]];
        } else {
            P0 = [-E[0], +E[1], -E[2] - s2n / D[2]];
        }

        if (s4n <= 0) {
            P1 = [-E[0], -E[1] - s4n / D[1], +E[2]];
        } else {
            P1 = [-E[0], -E[1], +E[2] - s4n / D[2]];
        }

        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (U1, -D)
    lower = (s6p >= 0 ? -E[0] : -E[0] - s6p / D[0]);
    upper = (s3p <= 0 ? +E[0] : +E[0] - s3p / D[0]);
    if (lower <= upper) {
        if (s6p >= 0) {
            P0 = [-E[0], +E[1], +E[2] - s6p / D[2]];
        } else {
            P0 = [-E[0] - s6p / D[0], +E[1], +E[2]];
        }

        if (s3p <= 0) {
            P1 = [+E[0], +E[1], -E[2] - s3p / D[2]];
        } else {
            // Upstream bug (fixed here): the C++ writes
            //   P1 = { +E[0] - s3p / D[0], -E[1], -E[2] }
            // but s3p is derived from t3, the projection of the box vertex
            // V3 = (+E[0], +E[1], -E[2]). Every other branch of every other
            // edge in this function builds the point from the vertex that
            // generated its s-value, replacing exactly one coordinate; the
            // sign of E[1] here is inconsistent with V3 and with the matching
            // 'lower/upper' clip of this same edge (which clips the x
            // coordinate). Compare the (D,-U1) block, which uses -E[1]
            // consistently because it is built from V1 = (+E[0],-E[1],-E[2]).
            P1 = [+E[0] - s3p / D[0], +E[1], -E[2]];
        }

        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (D, -U1)
    lower = (s4n >= 0 ? -E[0] : -E[0] - s4n / D[0]);
    upper = (s1n <= 0 ? +E[0] : +E[0] - s1n / D[0]);
    if (lower <= upper) {
        if (s4n >= 0) {
            P0 = [-E[0], -E[1], +E[2] - s4n / D[2]];
        } else {
            P0 = [-E[0] - s4n / D[0], -E[1], +E[2]];
        }

        if (s1n <= 0) {
            P1 = [+E[0], -E[1], -E[2] - s1n / D[2]];
        } else {
            P1 = [+E[0] - s1n / D[0], -E[1], -E[2]];
        }

        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (U2, -D)
    lower = (s5p >= 0 ? -E[1] : -E[1] - s5p / D[1]);
    upper = (s6p <= 0 ? +E[1] : +E[1] - s6p / D[1]);
    if (lower <= upper) {
        if (s5p >= 0) {
            P0 = [+E[0] - s5p / D[0], -E[1], +E[2]];
        } else {
            P0 = [+E[0], -E[1] - s5p / D[1], +E[2]];
        }

        if (s6p <= 0) {
            P1 = [-E[0] - s6p / D[0], +E[1], +E[2]];
        } else {
            P1 = [-E[0], E[1] - s6p / D[1], +E[2]];
        }

        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    // (D, -U2)
    lower = (s1n >= 0 ? -E[1] : -E[1] - s1n / D[1]);
    upper = (s2n <= 0 ? +E[1] : +E[1] - s2n / D[1]);
    if (lower <= upper) {
        if (s1n >= 0) {
            P0 = [+E[0] - s1n / D[0], -E[1], -E[2]];
        } else {
            P0 = [+E[0], -E[1] - s1n / D[1], -E[2]];
        }

        if (s2n <= 0) {
            P1 = [-E[0] - s2n / D[0], +E[1], -E[2]];
        } else {
            P1 = [-E[0], E[1] - s2n / D[1], -E[2]];
        }

        sqrDistance = computeSqrDistance(P0, P1, C, W0, W1);
        if (sqrDistance <= sqrRadius) {
            return true;
        }
    }

    return false;
}

export class IntrCanonicalBox3Cylinder3TI implements
    TIQuery<CanonicalBox, Cylinder3, IntrCanonicalBox3Cylinder3TIResult> {

    test(box: CanonicalBox, cylinder: Cylinder3):
        IntrCanonicalBox3Cylinder3TIResult {
        logAssert(box.dimension === 3,
            'IntrCanonicalBox3Cylinder3TI: mismatched sizes.');
        logAssert(cylinder.isFinite(),
            'Infinite cylinders are not yet supported.');

        // The result.intersect is initially false.
        const result = defaultTIResult();

        if (boxIsOutsideCylinderSlab(box, cylinder)) {
            // The box does not intersect the slab, so it does not intersect
            // the cylinder.
            result.intersect = false;
            return result;
        }

        // Apply reflections to obtain a cylinder whose axis direction is in
        // the first octant (positive- or zero-valued components). The
        // reflections applied to the canonical box do not require any
        // computational changes.
        const C = cylinder.axis.origin.values.slice();
        const D = cylinder.axis.direction.values.slice();
        const r = cylinder.radius;
        const hDiv2 = 0.5 * cylinder.height;
        const E = box.extent.values;
        for (let i = 0; i < 3; ++i) {
            if (D[i] < 0) {
                C[i] = -C[i];
                D[i] = -D[i];
            }
        }

        // D is now in the first octant. The box vertices are
        //   V[0] = (-E[0],-E[1],-E[2]), V[4] = (-E[0],-E[1],+E[2])
        //   V[1] = (+E[0],-E[1],-E[2]), V[5] = (+E[0],-E[1],+E[2])
        //   V[2] = (-E[0],+E[1],-E[2]), V[6] = (-E[0],+E[1],+E[2])
        //   V[3] = (+E[0],+E[1],-E[2]), V[7] = (+E[0],+E[1],+E[2])
        if (D[0] > 0) {
            if (D[1] > 0) {
                if (D[2] > 0) {
                    // (+,+,+)
                    result.intersect = doQueryNoZeros(C, D, r, hDiv2, E);
                } else {
                    // (+,+,0)
                    result.intersect = doQueryOneZero([0, 1, 2], C, D, r,
                        hDiv2, E);
                }
            } else {
                if (D[2] > 0) {
                    // (+,0,+)
                    result.intersect = doQueryOneZero([2, 0, 1], C, D, r,
                        hDiv2, E);
                } else {
                    // (+,0,0)
                    result.intersect = doQueryTwoZeros([0, 1, 2], C, r, E);
                }
            }
        } else {
            if (D[1] > 0) {
                if (D[2] > 0) {
                    // (0,+,+)
                    result.intersect = doQueryOneZero([1, 2, 0], C, D, r,
                        hDiv2, E);
                } else {
                    // (0,+,0)
                    result.intersect = doQueryTwoZeros([1, 2, 0], C, r, E);
                }
            } else {
                if (D[2] > 0) {
                    // (0,0,+)
                    result.intersect = doQueryTwoZeros([2, 0, 1], C, r, E);
                } else {
                    // (0,0,0)
                    logError('The cylinder direction cannot be (0,0,0).');
                }
            }
        }

        return result;
    }
}
