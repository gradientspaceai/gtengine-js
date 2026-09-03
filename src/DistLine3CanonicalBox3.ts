// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine3CanonicalBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a canonical box in 3D.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The canonical box has center at the origin and is aligned with the
// coordinate axes. The extents are E = (e[0],e[1],e[2]). A box point is
// Y = (y[0],y[1],y[2]) with |y[i]| <= e[i] for all i.
//
// The closest point on the line is stored in closest[0] with parameter t. The
// closest point on the box is stored in closest[1]. When there are infinitely
// many choices for the pair of closest points, only one of them is returned.
//
// The doQueryND functions are described in Section 10.9.4 Linear Component to
// Oriented Bounding Box of
//    Geometric Tools for Computer Graphics,
//    Philip J. Schneider and David H. Eberly,
//    Morgan Kaufmann, San Francisco CA, 2002
//
// Upstream TODO: the code in DistLine2AlignedBox2.h effectively uses the same
// approach, although in 2D, and is cleaner than this 3D code.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line3<T>, CanonicalBox3<T>>' becomes
// the class DistLine3CanonicalBox3 with the result type
// DistLine3CanonicalBox3Result. The private static helpers Face, DoQuery3D,
// DoQuery2D, DoQuery1D and DoQuery0D become module-private functions. As
// upstream does, they mutate the 'origin' vector in place: on return it holds
// the closest box point (in the reflected coordinate frame).

import type { CanonicalBox3 } from './CanonicalBox.js';
import type { DCPQuery } from './DCPQuery.js';
import type { Line3 } from './Line.js';
import { Vector, add, mul } from './Vector.js';

export interface DistLine3CanonicalBox3Result {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // closest[0] is on the line, closest[1] is on the box.
    closest: [Vector, Vector];
}

function face(i0: number, i1: number, i2: number, origin: Vector,
    direction: Vector, PmE: Vector, extent: Vector,
    result: DistLine3CanonicalBox3Result): void {
    const d = direction.values;
    const o = origin.values;
    const e = extent.values;
    const pme = PmE.values;
    const PpE = [o[0] + e[0], o[1] + e[1], o[2] + e[2]];

    if (d[i0] * PpE[i1] >= d[i1] * pme[i0]) {
        if (d[i0] * PpE[i2] >= d[i2] * pme[i0]) {
            // v[i1] >= -e[i1], v[i2] >= -e[i2] (distance = 0)
            o[i0] = e[i0];
            o[i1] -= d[i1] * pme[i0] / d[i0];
            o[i2] -= d[i2] * pme[i0] / d[i0];
            result.parameter = -pme[i0] / d[i0];
        }
        else {
            // v[i1] >= -e[i1], v[i2] < -e[i2]
            let lenSqr = d[i0] * d[i0] + d[i2] * d[i2];
            let tmp = lenSqr * PpE[i1]
                - d[i1] * (d[i0] * pme[i0] + d[i2] * PpE[i2]);
            if (tmp <= 2 * lenSqr * e[i1]) {
                const t = tmp / lenSqr;
                lenSqr += d[i1] * d[i1];
                tmp = PpE[i1] - t;
                const delta = d[i0] * pme[i0] + d[i1] * tmp + d[i2] * PpE[i2];
                result.parameter = -delta / lenSqr;
                result.sqrDistance += pme[i0] * pme[i0] + tmp * tmp
                    + PpE[i2] * PpE[i2] + delta * result.parameter;

                o[i0] = e[i0];
                o[i1] = t - e[i1];
                o[i2] = -e[i2];
            }
            else {
                lenSqr += d[i1] * d[i1];
                const delta = d[i0] * pme[i0] + d[i1] * pme[i1]
                    + d[i2] * PpE[i2];
                result.parameter = -delta / lenSqr;
                result.sqrDistance += pme[i0] * pme[i0] + pme[i1] * pme[i1]
                    + PpE[i2] * PpE[i2] + delta * result.parameter;

                o[i0] = e[i0];
                o[i1] = e[i1];
                o[i2] = -e[i2];
            }
        }
    }
    else {
        if (d[i0] * PpE[i2] >= d[i2] * pme[i0]) {
            // v[i1] < -e[i1], v[i2] >= -e[i2]
            let lenSqr = d[i0] * d[i0] + d[i1] * d[i1];
            let tmp = lenSqr * PpE[i2]
                - d[i2] * (d[i0] * pme[i0] + d[i1] * PpE[i1]);
            if (tmp <= 2 * lenSqr * e[i2]) {
                const t = tmp / lenSqr;
                lenSqr += d[i2] * d[i2];
                tmp = PpE[i2] - t;
                const delta = d[i0] * pme[i0] + d[i1] * PpE[i1] + d[i2] * tmp;
                result.parameter = -delta / lenSqr;
                result.sqrDistance += pme[i0] * pme[i0] + PpE[i1] * PpE[i1]
                    + tmp * tmp + delta * result.parameter;

                o[i0] = e[i0];
                o[i1] = -e[i1];
                o[i2] = t - e[i2];
            }
            else {
                lenSqr += d[i2] * d[i2];
                const delta = d[i0] * pme[i0] + d[i1] * PpE[i1]
                    + d[i2] * pme[i2];
                result.parameter = -delta / lenSqr;
                result.sqrDistance += pme[i0] * pme[i0] + PpE[i1] * PpE[i1]
                    + pme[i2] * pme[i2] + delta * result.parameter;

                o[i0] = e[i0];
                o[i1] = -e[i1];
                o[i2] = e[i2];
            }
        }
        else {
            // v[i1] < -e[i1], v[i2] < -e[i2]
            let lenSqr = d[i0] * d[i0] + d[i2] * d[i2];
            let tmp = lenSqr * PpE[i1]
                - d[i1] * (d[i0] * pme[i0] + d[i2] * PpE[i2]);
            if (tmp >= 0) {
                // v[i1]-edge is closest
                if (tmp <= 2 * lenSqr * e[i1]) {
                    const t = tmp / lenSqr;
                    lenSqr += d[i1] * d[i1];
                    tmp = PpE[i1] - t;
                    const delta = d[i0] * pme[i0] + d[i1] * tmp
                        + d[i2] * PpE[i2];
                    result.parameter = -delta / lenSqr;
                    result.sqrDistance += pme[i0] * pme[i0] + tmp * tmp
                        + PpE[i2] * PpE[i2] + delta * result.parameter;

                    o[i0] = e[i0];
                    o[i1] = t - e[i1];
                    o[i2] = -e[i2];
                }
                else {
                    lenSqr += d[i1] * d[i1];
                    const delta = d[i0] * pme[i0] + d[i1] * pme[i1]
                        + d[i2] * PpE[i2];
                    result.parameter = -delta / lenSqr;
                    result.sqrDistance += pme[i0] * pme[i0] + pme[i1] * pme[i1]
                        + PpE[i2] * PpE[i2] + delta * result.parameter;

                    o[i0] = e[i0];
                    o[i1] = e[i1];
                    o[i2] = -e[i2];
                }
                return;
            }

            lenSqr = d[i0] * d[i0] + d[i1] * d[i1];
            tmp = lenSqr * PpE[i2]
                - d[i2] * (d[i0] * pme[i0] + d[i1] * PpE[i1]);
            if (tmp >= 0) {
                // v[i2]-edge is closest
                if (tmp <= 2 * lenSqr * e[i2]) {
                    const t = tmp / lenSqr;
                    lenSqr += d[i2] * d[i2];
                    tmp = PpE[i2] - t;
                    const delta = d[i0] * pme[i0] + d[i1] * PpE[i1]
                        + d[i2] * tmp;
                    result.parameter = -delta / lenSqr;
                    result.sqrDistance += pme[i0] * pme[i0] + PpE[i1] * PpE[i1]
                        + tmp * tmp + delta * result.parameter;

                    o[i0] = e[i0];
                    o[i1] = -e[i1];
                    o[i2] = t - e[i2];
                }
                else {
                    lenSqr += d[i2] * d[i2];
                    const delta = d[i0] * pme[i0] + d[i1] * PpE[i1]
                        + d[i2] * pme[i2];
                    result.parameter = -delta / lenSqr;
                    result.sqrDistance += pme[i0] * pme[i0] + PpE[i1] * PpE[i1]
                        + pme[i2] * pme[i2] + delta * result.parameter;

                    o[i0] = e[i0];
                    o[i1] = -e[i1];
                    o[i2] = e[i2];
                }
                return;
            }

            // (v[i1],v[i2])-corner is closest
            lenSqr += d[i2] * d[i2];
            const delta = d[i0] * pme[i0] + d[i1] * PpE[i1] + d[i2] * PpE[i2];
            result.parameter = -delta / lenSqr;
            result.sqrDistance += pme[i0] * pme[i0] + PpE[i1] * PpE[i1]
                + PpE[i2] * PpE[i2] + delta * result.parameter;

            o[i0] = e[i0];
            o[i1] = -e[i1];
            o[i2] = -e[i2];
        }
    }
}

function doQuery3D(origin: Vector, direction: Vector, extent: Vector,
    result: DistLine3CanonicalBox3Result): void {
    const d = direction.values;
    const PmE = Vector.fromArray([
        origin.values[0] - extent.values[0],
        origin.values[1] - extent.values[1],
        origin.values[2] - extent.values[2]
    ]);
    const pme = PmE.values;
    const prodDxPy = d[0] * pme[1];
    const prodDyPx = d[1] * pme[0];

    if (prodDyPx >= prodDxPy) {
        const prodDzPx = d[2] * pme[0];
        const prodDxPz = d[0] * pme[2];
        if (prodDzPx >= prodDxPz) {
            // line intersects x = e0
            face(0, 1, 2, origin, direction, PmE, extent, result);
        }
        else {
            // line intersects z = e2
            face(2, 0, 1, origin, direction, PmE, extent, result);
        }
    }
    else {
        const prodDzPy = d[2] * pme[1];
        const prodDyPz = d[1] * pme[2];
        if (prodDzPy >= prodDyPz) {
            // line intersects y = e1
            face(1, 2, 0, origin, direction, PmE, extent, result);
        }
        else {
            // line intersects z = e2
            face(2, 0, 1, origin, direction, PmE, extent, result);
        }
    }
}

function doQuery2D(i0: number, i1: number, i2: number, origin: Vector,
    direction: Vector, extent: Vector,
    result: DistLine3CanonicalBox3Result): void {
    const o = origin.values;
    const d = direction.values;
    const e = extent.values;

    const PmE0 = o[i0] - e[i0];
    const PmE1 = o[i1] - e[i1];
    const prod0 = d[i1] * PmE0;
    const prod1 = d[i0] * PmE1;

    if (prod0 >= prod1) {
        // line intersects P[i0] = e[i0]
        o[i0] = e[i0];

        const PpE1 = o[i1] + e[i1];
        const delta = prod0 - d[i0] * PpE1;
        if (delta >= 0) {
            const lenSqr = d[i0] * d[i0] + d[i1] * d[i1];
            result.sqrDistance += delta * delta / lenSqr;
            o[i1] = -e[i1];
            result.parameter = -(d[i0] * PmE0 + d[i1] * PpE1) / lenSqr;
        }
        else {
            o[i1] -= prod0 / d[i0];
            result.parameter = -PmE0 / d[i0];
        }
    }
    else {
        // line intersects P[i1] = e[i1]
        o[i1] = e[i1];

        const PpE0 = o[i0] + e[i0];
        const delta = prod1 - d[i1] * PpE0;
        if (delta >= 0) {
            const lenSqr = d[i0] * d[i0] + d[i1] * d[i1];
            result.sqrDistance += delta * delta / lenSqr;
            o[i0] = -e[i0];
            result.parameter = -(d[i0] * PpE0 + d[i1] * PmE1) / lenSqr;
        }
        else {
            o[i0] -= prod1 / d[i1];
            result.parameter = -PmE1 / d[i1];
        }
    }

    if (o[i2] < -e[i2]) {
        const delta = o[i2] + e[i2];
        result.sqrDistance += delta * delta;
        o[i2] = -e[i2];
    }
    else if (o[i2] > e[i2]) {
        const delta = o[i2] - e[i2];
        result.sqrDistance += delta * delta;
        o[i2] = e[i2];
    }
}

function doQuery1D(i0: number, i1: number, i2: number, origin: Vector,
    direction: Vector, extent: Vector,
    result: DistLine3CanonicalBox3Result): void {
    const o = origin.values;
    const e = extent.values;

    result.parameter = (e[i0] - o[i0]) / direction.values[i0];

    o[i0] = e[i0];
    for (const i of [i1, i2]) {
        if (o[i] < -e[i]) {
            const delta = o[i] + e[i];
            result.sqrDistance += delta * delta;
            o[i] = -e[i];
        }
        else if (o[i] > e[i]) {
            const delta = o[i] - e[i];
            result.sqrDistance += delta * delta;
            o[i] = e[i];
        }
    }
}

function doQuery0D(origin: Vector, extent: Vector,
    result: DistLine3CanonicalBox3Result): void {
    const o = origin.values;
    const e = extent.values;
    for (let i = 0; i < 3; ++i) {
        if (o[i] < -e[i]) {
            const delta = o[i] + e[i];
            result.sqrDistance += delta * delta;
            o[i] = -e[i];
        }
        else if (o[i] > e[i]) {
            const delta = o[i] - e[i];
            result.sqrDistance += delta * delta;
            o[i] = e[i];
        }
    }
}

export class DistLine3CanonicalBox3
    implements DCPQuery<Line3, CanonicalBox3, DistLine3CanonicalBox3Result> {
    compute(line: Line3, box: CanonicalBox3): DistLine3CanonicalBox3Result {
        const result: DistLine3CanonicalBox3Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            closest: [new Vector(3), new Vector(3)]
        };

        // Copies are made so that the line direction can be transformed to
        // the first octant (nonnegative components) using reflections. The
        // line parameter t is invariant under these reflections, so the
        // closest line point can be computed from the original line.
        const origin = line.origin.clone();
        const direction = line.direction.clone();
        const reflect: boolean[] = [false, false, false];
        for (let i = 0; i < 3; ++i) {
            if (direction.values[i] < 0) {
                origin.values[i] = -origin.values[i];
                direction.values[i] = -direction.values[i];
                reflect[i] = true;
            }
        }

        // Compute the line-box distance and closest points. The doQueryND
        // calls compute result.parameter and update result.sqrDistance. The
        // result.distance is computed after the specialized queries. The
        // result.closest[] points are computed afterwards.
        if (direction.values[0] > 0) {
            if (direction.values[1] > 0) {
                if (direction.values[2] > 0) {  // (+,+,+)
                    doQuery3D(origin, direction, box.extent, result);
                }
                else {  // (+,+,0)
                    doQuery2D(0, 1, 2, origin, direction, box.extent, result);
                }
            }
            else {
                if (direction.values[2] > 0) {  // (+,0,+)
                    doQuery2D(0, 2, 1, origin, direction, box.extent, result);
                }
                else {  // (+,0,0)
                    doQuery1D(0, 1, 2, origin, direction, box.extent, result);
                }
            }
        }
        else {
            if (direction.values[1] > 0) {
                if (direction.values[2] > 0) {  // (0,+,+)
                    doQuery2D(1, 2, 0, origin, direction, box.extent, result);
                }
                else {  // (0,+,0)
                    doQuery1D(1, 0, 2, origin, direction, box.extent, result);
                }
            }
            else {
                if (direction.values[2] > 0) {  // (0,0,+)
                    doQuery1D(2, 0, 1, origin, direction, box.extent, result);
                }
                else {  // (0,0,0)
                    doQuery0D(origin, box.extent, result);
                }
            }
        }

        // Undo the reflections applied previously.
        for (let i = 0; i < 3; ++i) {
            if (reflect[i]) {
                origin.values[i] = -origin.values[i];
            }
        }

        result.distance = Math.sqrt(result.sqrDistance);

        // Compute the closest point on the line.
        result.closest[0] =
            add(line.origin, mul(result.parameter, line.direction));

        // Compute the closest point on the box. The 'origin' copy is modified
        // by the doQueryND functions to hold the closest box point.
        result.closest[1] = origin;
        return result;
    }
}
