// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPlane3CanonicalBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a plane and a solid canonical box in 3D.
//
// The plane is defined by Dot(N, X - P) = 0, where P is the plane origin and
// N is a unit-length normal for the plane.
//
// The canonical box has center at the origin and is aligned with the
// coordinate axes. The extents are E = (e[0],e[1],e[2]). A box point is
// Y = (y[0],y[1],y[2]) with |y[i]| <= e[i] for all i.
//
// The closest point on the plane is stored in closest[0]. The closest point
// on the box is stored in closest[1]. When there are infinitely many choices
// for the pair of closest points, only one of them is returned.
//
// TODO (upstream): Modify to support non-unit-length N.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Plane3<T>, CanonicalBox3<T>>' becomes
// the class DistPlane3CanonicalBox3 with the result type
// DistPlane3CanonicalBox3Result. The private static helpers DoQuery3D,
// DoQuery2D, DoQuery1D and DoQuery0D become the module-private functions
// doQuery3D, doQuery2D, doQuery1D and doQuery0D. Where upstream assigns one
// closest point to the other (C++ value semantics), the port clones.

import type { CanonicalBox3 } from './CanonicalBox.js';
import type { DCPQuery } from './DCPQuery.js';
import { clamp } from './Functions.js';
import type { Plane3 } from './Hyperplane.js';
import { Vector, add, dot, mul, negate, sub } from './Vector.js';

export interface DistPlane3CanonicalBox3Result {
    distance: number;
    sqrDistance: number;

    // closest[0] is on the plane, closest[1] is on the box.
    closest: [Vector, Vector];
}

function doQuery3D(origin: Vector, normal: Vector, extent: Vector,
    result: DistPlane3CanonicalBox3Result): void {
    const dmin = -dot(normal, add(extent, origin));
    if (dmin >= 0) {
        result.closest[0] = sub(negate(extent), mul(dmin, normal));
        result.closest[1] = negate(extent);
    }
    else {
        const dmax = dot(normal, sub(extent, origin));
        if (dmax <= 0) {
            result.closest[0] = sub(extent, mul(dmax, normal));
            result.closest[1] = extent.clone();
        }
        else {
            // t = dmin / (dmin - dmax) in [0,1], compute s = 2*t-1
            const s = 2 * dmin / (dmin - dmax) - 1;
            result.closest[0] = mul(s, extent);
            result.closest[1] = result.closest[0].clone();
        }
    }
}

function doQuery2D(i0: number, i1: number, i2: number, origin: Vector,
    normal: Vector, extent: Vector,
    result: DistPlane3CanonicalBox3Result): void {
    const e = extent.values;
    const n = normal.values;
    const p = origin.values;
    const c0 = result.closest[0].values;
    const c1 = result.closest[1].values;

    const dmin = -(n[i0] * (e[i0] + p[i0]) + n[i1] * (e[i1] + p[i1]));

    if (dmin >= 0) {
        c0[i0] = -e[i0] - dmin * n[i0];
        c0[i1] = -e[i1] - dmin * n[i1];
        c0[i2] = e[i2];
        c1[i0] = -e[i0];
        c1[i1] = -e[i1];
        c1[i2] = e[i2];
    }
    else {
        const dmax = n[i0] * (e[i0] - p[i0]) + n[i1] * (e[i1] - p[i1]);

        if (dmax <= 0) {
            c0[i0] = e[i0] - dmax * n[i0];
            c0[i1] = e[i1] - dmax * n[i1];
            c0[i2] = e[i2];
            c1[i0] = e[i0];
            c1[i1] = e[i1];
            c1[i2] = e[i2];
        }
        else {
            // t = dmin / (dmin - dmax) in [0,1], compute s = 2*t-1
            const s = 2 * dmin / (dmin - dmax) - 1;
            c0[i0] = s * e[i0];
            c0[i1] = s * e[i1];
            c0[i2] = e[i2];
            result.closest[1] = result.closest[0].clone();
        }
    }
}

function doQuery1D(i0: number, i1: number, i2: number, origin: Vector,
    extent: Vector, result: DistPlane3CanonicalBox3Result): void {
    const e = extent.values;
    const p = origin.values;
    const c0 = result.closest[0].values;
    const c1 = result.closest[1].values;
    c0[i0] = p[i0];
    c0[i1] = e[i1];
    c0[i2] = e[i2];
    c1[i0] = clamp(p[i0], -e[i0], e[i0]);
    c1[i1] = e[i1];
    c1[i2] = e[i2];
}

function doQuery0D(origin: Vector, extent: Vector,
    result: DistPlane3CanonicalBox3Result): void {
    const e = extent.values;
    const p = origin.values;
    result.closest[0] = origin.clone();
    const c1 = result.closest[1].values;
    c1[0] = clamp(p[0], -e[0], e[0]);
    c1[1] = clamp(p[1], -e[1], e[1]);
    c1[2] = clamp(p[2], -e[2], e[2]);
}

export class DistPlane3CanonicalBox3
    implements DCPQuery<Plane3, CanonicalBox3, DistPlane3CanonicalBox3Result> {
    compute(plane: Plane3, box: CanonicalBox3):
        DistPlane3CanonicalBox3Result {
        const result: DistPlane3CanonicalBox3Result = {
            distance: 0,
            sqrDistance: 0,
            closest: [new Vector(3), new Vector(3)]
        };

        // Copies are made so that we can transform the plane normal to the
        // first octant (nonnegative components) using reflections.
        const origin = mul(plane.constant, plane.normal);
        const normal = plane.normal.clone();
        const reflect: boolean[] = [false, false, false];
        for (let i = 0; i < 3; ++i) {
            if (normal.values[i] < 0) {
                origin.values[i] = -origin.values[i];
                normal.values[i] = -normal.values[i];
                reflect[i] = true;
            }
        }

        // Compute the plane-box closest points.
        const n = normal.values;
        if (n[0] > 0) {
            if (n[1] > 0) {
                if (n[2] > 0) {
                    // The normal signs are (+,+,+).
                    doQuery3D(origin, normal, box.extent, result);
                }
                else {
                    // The normal signs are (+,+,0).
                    doQuery2D(0, 1, 2, origin, normal, box.extent, result);
                }
            }
            else {
                if (n[2] > 0) {
                    // The normal signs are (+,0,+).
                    doQuery2D(0, 2, 1, origin, normal, box.extent, result);
                }
                else {
                    // The normal signs are (+,0,0). The closest box point is
                    // (x0,e1,e2) where x0 = clamp(p0,[-e0,e0]). The closest
                    // plane point is (p0,e1,e2).
                    doQuery1D(0, 1, 2, origin, box.extent, result);
                }
            }
        }
        else {
            if (n[1] > 0) {
                if (n[2] > 0) {
                    // The normal signs are (0,+,+).
                    doQuery2D(1, 2, 0, origin, normal, box.extent, result);
                }
                else {
                    // The normal signs are (0,+,0). The closest box point is
                    // (e0,x1,e2) where x1 = clamp(p1,[-e1,e1]). The closest
                    // plane point is (e0,p1,e2).
                    doQuery1D(1, 2, 0, origin, box.extent, result);
                }
            }
            else {
                if (n[2] > 0) {
                    // The normal signs are (0,0,+). The closest box point is
                    // (e0,e1,x2) where x2 = clamp(p2,[-e2,e2]). The closest
                    // plane point is (e0,e1,p2).
                    doQuery1D(2, 0, 1, origin, box.extent, result);
                }
                else {
                    // The normal signs are (0,0,0). Execute the DCP query for
                    // the plane origin and the canonical box. This is a
                    // low-probability event. Upstream passes the stored
                    // plane.origin here rather than the local (reflected)
                    // origin; with a zero normal no reflections were applied,
                    // so the port keeps the upstream expression.
                    doQuery0D(plane.origin, box.extent, result);
                }
            }
        }

        // Undo the reflections. The origin and normal are not consumed, so
        // these do not need to be reflected. However, the closest points are
        // consumed.
        for (let i = 0; i < 3; ++i) {
            if (reflect[i]) {
                for (let j = 0; j < 2; ++j) {
                    result.closest[j].values[i] = -result.closest[j].values[i];
                }
            }
        }

        const diff = sub(result.closest[0], result.closest[1]);
        result.sqrDistance = dot(diff, diff);
        result.distance = Math.sqrt(result.sqrDistance);
        return result;
    }
}
