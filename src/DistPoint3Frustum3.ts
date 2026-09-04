// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint3Frustum3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance from a point to an orthogonal frustum. The algorithm
// is described in
// https://www.geometrictools.com/Documentation/DistancePointToFrustum.pdf
//
// The input point is stored in the member closest[0]. The frustum point
// closest to it is stored in the member closest[1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector3<T>, Frustum3<T>>' becomes the
// class DistPoint3Frustum3 with the result type DistPoint3Frustum3Result.

import type { DCPQuery } from './DCPQuery.js';
import type { Frustum3 } from './Frustum3.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

export interface DistPoint3Frustum3Result {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest frustum point.
    closest: [Vector, Vector];
}

// Upstream bug (fixed here): the LF-edge and UF-edge cases assign the free
// coordinate of the far edge straight from the (sign-folded) test point,
// which is correct only where the surrounding case analysis already bounds
// it. Two of the ten assignments are not so bounded: in the test[2] <= dmin
// group and in the dmin < test[2] < dmax group, the branch guarded by
// rEdgeDot >= 0 followed by rdDot >= maxRDDot reaches the LF-edge with
// test[1] > umax (and symmetrically the uEdgeDot branch reaches the UF-edge
// with test[0] > rmax). The reported closest point is then off the end of the
// edge and outside the frustum, and the distance is too small.
//
// Example: an orthogonal frustum with dMin = 1, dMax = 2, uBound = rBound = 1
// and the test point at frustum coordinates (r,u,d) = (100,3,0). Upstream
// returns closest = (2,3,2) with distance 98.0204; the frustum only reaches
// |u| <= 2 at d = 2, and the true closest point is the LUF vertex (2,2,2) at
// distance 98.0255.
//
// The closest point on a segment is the projection clamped to its endpoints,
// so clamping the free coordinate to the far half-extent repairs those cases
// and is a no-op for the eight that were already bounded. The test point has
// been folded into the octant with nonnegative R and U coordinates, so only
// the upper bound needs to be applied.
function clampToFarEdge(coordinate: number, halfExtent: number): number {
    return coordinate <= halfExtent ? coordinate : halfExtent;
}

export class DistPoint3Frustum3
    implements DCPQuery<Vector, Frustum3, DistPoint3Frustum3Result> {
    compute(point: Vector, frustum: Frustum3): DistPoint3Frustum3Result {
        // Compute coordinates of the point with respect to the frustum
        // coordinate system.
        const diffOrigin = sub(point, frustum.origin);
        const test = [
            dot(diffOrigin, frustum.rVector),
            dot(diffOrigin, frustum.uVector),
            dot(diffOrigin, frustum.dVector)
        ];

        // Perform calculations in the octant with nonnegative R and U
        // coordinates.
        let rSignChange = false;
        if (test[0] < 0) {
            rSignChange = true;
            test[0] = -test[0];
        }

        let uSignChange = false;
        if (test[1] < 0) {
            uSignChange = true;
            test[1] = -test[1];
        }

        // Frustum derived parameters.
        const rmin = frustum.rBound;
        const rmax = frustum.getDRatio() * rmin;
        const umin = frustum.uBound;
        const umax = frustum.getDRatio() * umin;
        const dmin = frustum.dMin;
        const dmax = frustum.dMax;
        const rminSqr = rmin * rmin;
        const uminSqr = umin * umin;
        const dminSqr = dmin * dmin;
        const minRDDot = rminSqr + dminSqr;
        const minUDDot = uminSqr + dminSqr;
        const minRUDDot = rminSqr + minUDDot;
        const maxRDDot = frustum.getDRatio() * minRDDot;
        const maxUDDot = frustum.getDRatio() * minUDDot;
        const maxRUDDot = frustum.getDRatio() * minRUDDot;

        // The algorithm computes the closest point in all cases by
        // determining in which Voronoi region of the vertices, edges and
        // faces of the frustum the test point lives.
        const closest = [0, 0, 0];
        let rDot: number, uDot: number, rdDot: number, udDot: number;
        let rudDot: number, rEdgeDot: number, uEdgeDot: number, t: number;
        if (test[2] >= dmax) {
            if (test[0] <= rmax) {
                if (test[1] <= umax) {
                    // F-face
                    closest[0] = test[0];
                    closest[1] = test[1];
                    closest[2] = dmax;
                }
                else {
                    // UF-edge
                    closest[0] = clampToFarEdge(test[0], rmax);
                    closest[1] = umax;
                    closest[2] = dmax;
                }
            }
            else {
                if (test[1] <= umax) {
                    // LF-edge
                    closest[0] = rmax;
                    closest[1] = clampToFarEdge(test[1], umax);
                    closest[2] = dmax;
                }
                else {
                    // LUF-vertex
                    closest[0] = rmax;
                    closest[1] = umax;
                    closest[2] = dmax;
                }
            }
        }
        else if (test[2] <= dmin) {
            if (test[0] <= rmin) {
                if (test[1] <= umin) {
                    // N-face
                    closest[0] = test[0];
                    closest[1] = test[1];
                    closest[2] = dmin;
                }
                else {
                    udDot = umin * test[1] + dmin * test[2];
                    if (udDot >= maxUDDot) {
                        // UF-edge
                        closest[0] = clampToFarEdge(test[0], rmax);
                        closest[1] = umax;
                        closest[2] = dmax;
                    }
                    else if (udDot >= minUDDot) {
                        // U-face
                        uDot = dmin * test[1] - umin * test[2];
                        t = uDot / minUDDot;
                        closest[0] = test[0];
                        closest[1] = test[1] - t * dmin;
                        closest[2] = test[2] + t * umin;
                    }
                    else {
                        // UN-edge
                        closest[0] = test[0];
                        closest[1] = umin;
                        closest[2] = dmin;
                    }
                }
            }
            else {
                if (test[1] <= umin) {
                    rdDot = rmin * test[0] + dmin * test[2];
                    if (rdDot >= maxRDDot) {
                        // LF-edge
                        closest[0] = rmax;
                        closest[1] = clampToFarEdge(test[1], umax);
                        closest[2] = dmax;
                    }
                    else if (rdDot >= minRDDot) {
                        // L-face
                        rDot = dmin * test[0] - rmin * test[2];
                        t = rDot / minRDDot;
                        closest[0] = test[0] - t * dmin;
                        closest[1] = test[1];
                        closest[2] = test[2] + t * rmin;
                    }
                    else {
                        // LN-edge
                        closest[0] = rmin;
                        closest[1] = test[1];
                        closest[2] = dmin;
                    }
                }
                else {
                    rudDot = rmin * test[0] + umin * test[1] + dmin * test[2];
                    rEdgeDot = umin * rudDot - minRUDDot * test[1];
                    if (rEdgeDot >= 0) {
                        rdDot = rmin * test[0] + dmin * test[2];
                        if (rdDot >= maxRDDot) {
                            // LF-edge
                            closest[0] = rmax;
                            closest[1] = clampToFarEdge(test[1], umax);
                            closest[2] = dmax;
                        }
                        else if (rdDot >= minRDDot) {
                            // L-face
                            rDot = dmin * test[0] - rmin * test[2];
                            t = rDot / minRDDot;
                            closest[0] = test[0] - t * dmin;
                            closest[1] = test[1];
                            closest[2] = test[2] + t * rmin;
                        }
                        else {
                            // LN-edge
                            closest[0] = rmin;
                            closest[1] = test[1];
                            closest[2] = dmin;
                        }
                    }
                    else {
                        uEdgeDot = rmin * rudDot - minRUDDot * test[0];
                        if (uEdgeDot >= 0) {
                            udDot = umin * test[1] + dmin * test[2];
                            if (udDot >= maxUDDot) {
                                // UF-edge
                                closest[0] = clampToFarEdge(test[0], rmax);
                                closest[1] = umax;
                                closest[2] = dmax;
                            }
                            else if (udDot >= minUDDot) {
                                // U-face
                                uDot = dmin * test[1] - umin * test[2];
                                t = uDot / minUDDot;
                                closest[0] = test[0];
                                closest[1] = test[1] - t * dmin;
                                closest[2] = test[2] + t * umin;
                            }
                            else {
                                // UN-edge
                                closest[0] = test[0];
                                closest[1] = umin;
                                closest[2] = dmin;
                            }
                        }
                        else {
                            if (rudDot >= maxRUDDot) {
                                // LUF-vertex
                                closest[0] = rmax;
                                closest[1] = umax;
                                closest[2] = dmax;
                            }
                            else if (rudDot >= minRUDDot) {
                                // LU-edge
                                t = rudDot / minRUDDot;
                                closest[0] = t * rmin;
                                closest[1] = t * umin;
                                closest[2] = t * dmin;
                            }
                            else {
                                // LUN-vertex
                                closest[0] = rmin;
                                closest[1] = umin;
                                closest[2] = dmin;
                            }
                        }
                    }
                }
            }
        }
        else {
            rDot = dmin * test[0] - rmin * test[2];
            uDot = dmin * test[1] - umin * test[2];
            if (rDot <= 0) {
                if (uDot <= 0) {
                    // The point is inside the frustum.
                    closest[0] = test[0];
                    closest[1] = test[1];
                    closest[2] = test[2];
                }
                else {
                    udDot = umin * test[1] + dmin * test[2];
                    if (udDot >= maxUDDot) {
                        // UF-edge
                        closest[0] = clampToFarEdge(test[0], rmax);
                        closest[1] = umax;
                        closest[2] = dmax;
                    }
                    else {
                        // U-face
                        t = uDot / minUDDot;
                        closest[0] = test[0];
                        closest[1] = test[1] - t * dmin;
                        closest[2] = test[2] + t * umin;
                    }
                }
            }
            else {
                if (uDot <= 0) {
                    rdDot = rmin * test[0] + dmin * test[2];
                    if (rdDot >= maxRDDot) {
                        // LF-edge
                        closest[0] = rmax;
                        closest[1] = clampToFarEdge(test[1], umax);
                        closest[2] = dmax;
                    }
                    else {
                        // L-face
                        t = rDot / minRDDot;
                        closest[0] = test[0] - t * dmin;
                        closest[1] = test[1];
                        closest[2] = test[2] + t * rmin;
                    }
                }
                else {
                    rudDot = rmin * test[0] + umin * test[1] + dmin * test[2];
                    rEdgeDot = umin * rudDot - minRUDDot * test[1];
                    if (rEdgeDot >= 0) {
                        rdDot = rmin * test[0] + dmin * test[2];
                        if (rdDot >= maxRDDot) {
                            // LF-edge
                            closest[0] = rmax;
                            closest[1] = clampToFarEdge(test[1], umax);
                            closest[2] = dmax;
                        }
                        else {
                            // L-face; here rdDot >= minRDDot.
                            t = rDot / minRDDot;
                            closest[0] = test[0] - t * dmin;
                            closest[1] = test[1];
                            closest[2] = test[2] + t * rmin;
                        }
                    }
                    else {
                        uEdgeDot = rmin * rudDot - minRUDDot * test[0];
                        if (uEdgeDot >= 0) {
                            udDot = umin * test[1] + dmin * test[2];
                            if (udDot >= maxUDDot) {
                                // UF-edge
                                closest[0] = clampToFarEdge(test[0], rmax);
                                closest[1] = umax;
                                closest[2] = dmax;
                            }
                            else {
                                // U-face; here udDot >= minUDDot.
                                t = uDot / minUDDot;
                                closest[0] = test[0];
                                closest[1] = test[1] - t * dmin;
                                closest[2] = test[2] + t * umin;
                            }
                        }
                        else {
                            if (rudDot >= maxRUDDot) {
                                // LUF-vertex
                                closest[0] = rmax;
                                closest[1] = umax;
                                closest[2] = dmax;
                            }
                            else {
                                // LU-edge; here rudDot >= minRUDDot.
                                t = rudDot / minRUDDot;
                                closest[0] = t * rmin;
                                closest[1] = t * umin;
                                closest[2] = t * dmin;
                            }
                        }
                    }
                }
            }
        }

        // The difference is computed in the folded octant, before the signs
        // are restored.
        const diff = [
            test[0] - closest[0],
            test[1] - closest[1],
            test[2] - closest[2]
        ];

        // Convert back to the original quadrant.
        if (rSignChange) {
            closest[0] = -closest[0];
        }

        if (uSignChange) {
            closest[1] = -closest[1];
        }

        // Convert back to the original coordinates.
        const closest1 = add(frustum.origin,
            add(mul(closest[0], frustum.rVector),
                add(mul(closest[1], frustum.uVector),
                    mul(closest[2], frustum.dVector))));

        const sqrDistance = diff[0] * diff[0] + diff[1] * diff[1]
            + diff[2] * diff[2];
        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            closest: [point.clone(), closest1]
        };
    }
}
