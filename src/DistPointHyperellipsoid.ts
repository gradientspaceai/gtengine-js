// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointHyperellipsoid.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance from a point to a hyperellipsoid in nD. The
// hyperellipsoid is considered to be a closed surface, not a solid. In 2D,
// this is a point-ellipse distance query. In 3D, this is a point-ellipsoid
// distance query. The following document describes the algorithm.
//   https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
// The hyperellipsoid can have arbitrary center and orientation; that is, it
// does not have to be axis-aligned with center at the origin.
//
// The input point is stored in closest[0]. The closest point on the
// hyperellipsoid is stored in closest[1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, Hyperellipsoid<N,T>>'
// becomes the class DistPointHyperellipsoid with the result type
// DistPointHyperellipsoidResult. Upstream has two 'operator()' overloads; the
// canonical point-hyperellipsoid query keeps the name 'compute' and the
// overload that takes only the extents (an axis-aligned hyperellipsoid
// centered at the origin) becomes 'computeAxisAligned', following the
// precedent for multiple upstream query overloads. The private helpers
// SqrDistance, SqrDistanceSpecial and Bisector become the module-private
// functions sqrDistance, sqrDistanceSpecial and bisector. The 'Vector<N-1,T>'
// temporaries in SqrDistanceSpecial become plain number[] of length N-1. The
// dimension aliases are dropped because the runtime-dimension Vector serves
// every N.

import type { DCPQuery } from './DCPQuery.js';
import type { Hyperellipsoid } from './Hyperellipsoid.js';
import { logAssert } from './Logger.js';
import { Vector, add, dot, length, mul, sub } from './Vector.js';

export interface DistPointHyperellipsoidResult {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest
    // hyperellipsoid point.
    closest: [Vector, Vector];
}

// The bisection algorithm to find the unique root of F(t). Only the first
// 'numComponents' entries of e, y and x participate.
//
// Upstream numerical caveat (preserved, not fixed): the root satisfies
// s > -pSqr[numComponents-1] = -1, and the lower bracket is
// smin = z[numComponents-1] - 1. When z[numComponents-1] (the query-point
// coordinate along the smallest-extent axis, divided by that extent) is tiny
// but nonzero, smin rounds to exactly -1 and the root itself is closer to -1
// than the double-precision spacing there. The final division by
// (s + pSqr[i]) then cancels catastrophically: the reported closest point
// leaves the hyperellipsoid and, when s + pSqr[i] underflows to zero, the
// distance is Infinity. A coordinate of exactly zero takes the separate
// SqrDistanceSpecial branch and is handled correctly, so the defect is
// confined to a narrow band around the axis hyperplanes of the smallest
// extent (and hence around the center). Upstream has the same behavior; it
// offers exact arithmetic (T = BSRational) as the remedy. Fixing this in
// double precision would require reformulating the bisection variable, which
// is a change to upstream's algorithm rather than a translation fix.
function bisector(numComponents: number, e: Vector, y: Vector,
    x: Vector): number {
    let sumZSqr = 0;
    const z = new Array<number>(numComponents).fill(0);
    for (let i = 0; i < numComponents; ++i) {
        z[i] = y.values[i] / e.values[i];
        sumZSqr += z[i] * z[i];
    }

    if (sumZSqr === 1) {
        // The point is on the hyperellipsoid.
        for (let i = 0; i < numComponents; ++i) {
            x.values[i] = y.values[i];
        }
        return 0;
    }

    const emin = e.values[numComponents - 1];
    const pSqr = new Vector(y.size);
    const numerator = new Vector(y.size);
    for (let i = 0; i < numComponents; ++i) {
        const p = e.values[i] / emin;
        pSqr.values[i] = p * p;
        numerator.values[i] = pSqr.values[i] * z[i];
    }

    let s = 0;
    let smin = z[numComponents - 1] - 1;
    let smax: number;
    if (sumZSqr < 1) {
        // The point is strictly inside the hyperellipsoid.
        smax = 0;
    }
    else {
        // The point is strictly outside the hyperellipsoid.
        smax = length(numerator, true) - 1;
    }

    const jmax = 2048;
    for (let j = 0; j < jmax; ++j) {
        s = 0.5 * (smin + smax);
        if (s === smin || s === smax) {
            break;
        }

        let g = -1;
        for (let i = 0; i < numComponents; ++i) {
            const ratio = numerator.values[i] / (s + pSqr.values[i]);
            g += ratio * ratio;
        }

        if (g > 0) {
            smin = s;
        }
        else if (g < 0) {
            smax = s;
        }
        else {
            break;
        }
    }

    let sqrDist = 0;
    for (let i = 0; i < numComponents; ++i) {
        x.values[i] = pSqr.values[i] * y.values[i] / (s + pSqr.values[i]);
        const diff = x.values[i] - y.values[i];
        sqrDist += diff * diff;
    }
    return sqrDist;
}

// The hyperellipsoid is sum_{d=0}^{N-1} (x[d]/e[d])^2 = 1 with the e[d]
// positive and nonincreasing: e[d] >= e[d + 1] for all d. The query point is
// (y[0],...,y[N-1]) with y[d] >= 0 for all d. The function returns the
// squared distance from the query point to the hyperellipsoid. It also
// computes the hyperellipsoid point (x[0],...,x[N-1]) that is closest to
// (y[0],...,y[N-1]), where x[d] >= 0 for all d.
function sqrDistanceSpecial(e: Vector, y: Vector, x: Vector): number {
    const n = e.size;
    let sqrDist = 0;

    const ePos = new Vector(n);
    const yPos = new Vector(n);
    const xPos = new Vector(n);
    let numPos = 0;
    for (let i = 0; i < n; ++i) {
        if (y.values[i] > 0) {
            ePos.values[numPos] = e.values[i];
            yPos.values[numPos] = y.values[i];
            ++numPos;
        }
        else {
            x.values[i] = 0;
        }
    }

    if (y.values[n - 1] > 0) {
        sqrDist = bisector(numPos, ePos, yPos, xPos);
    }
    else {
        // y[N-1] = 0
        const numer = new Array<number>(n - 1).fill(0);
        const denom = new Array<number>(n - 1).fill(0);
        const eNm1Sqr = e.values[n - 1] * e.values[n - 1];
        for (let i = 0; i < numPos; ++i) {
            numer[i] = ePos.values[i] * yPos.values[i];
            denom[i] = ePos.values[i] * ePos.values[i] - eNm1Sqr;
        }

        let inSubHyperbox = true;
        for (let i = 0; i < numPos; ++i) {
            if (numer[i] >= denom[i]) {
                inSubHyperbox = false;
                break;
            }
        }

        let inSubHyperellipsoid = false;
        if (inSubHyperbox) {
            // yPos[] is inside the axis-aligned bounding box of the
            // subhyperellipsoid. This intermediate test is designed to guard
            // against the division by zero when ePos[i] == e[N-1] for some i.
            const xde = new Array<number>(n - 1).fill(0);
            let discr = 1;
            for (let i = 0; i < numPos; ++i) {
                xde[i] = numer[i] / denom[i];
                discr -= xde[i] * xde[i];
            }
            if (discr > 0) {
                // yPos[] is inside the subhyperellipsoid. The closest
                // hyperellipsoid point has x[N-1] > 0.
                sqrDist = 0;
                for (let i = 0; i < numPos; ++i) {
                    xPos.values[i] = ePos.values[i] * xde[i];
                    const diff = xPos.values[i] - yPos.values[i];
                    sqrDist += diff * diff;
                }
                x.values[n - 1] = e.values[n - 1] * Math.sqrt(discr);
                sqrDist += x.values[n - 1] * x.values[n - 1];
                inSubHyperellipsoid = true;
            }
        }

        if (!inSubHyperellipsoid) {
            // yPos[] is outside the subhyperellipsoid. The closest
            // hyperellipsoid point has x[N-1] == 0 and is on the
            // domain-boundary hyperellipsoid.
            x.values[n - 1] = 0;
            sqrDist = bisector(numPos, ePos, yPos, xPos);
        }
    }

    // Fill in those x[] values that were not zeroed out initially.
    numPos = 0;
    for (let i = 0; i < n; ++i) {
        if (y.values[i] > 0) {
            x.values[i] = xPos.values[numPos];
            ++numPos;
        }
    }

    return sqrDist;
}

// The hyperellipsoid is sum_{d=0}^{N-1} (x[d]/e[d])^2 = 1 with no constraints
// on the ordering of the e[d]. The query point is (y[0],...,y[N-1]) with no
// constraints on the signs of the components. The function returns the
// squared distance from the query point to the hyperellipsoid. It also
// computes the hyperellipsoid point (x[0],...,x[N-1]) that is closest to
// (y[0],...,y[N-1]).
function sqrDistance(e: Vector, y: Vector, x: Vector): number {
    const n = e.size;

    // Determine negations for y to the first octant.
    const negate = new Array<boolean>(n).fill(false);
    for (let i = 0; i < n; ++i) {
        negate[i] = (y.values[i] < 0);
    }

    // Determine the axis order for decreasing extents. Upstream sorts an
    // array of std::pair<T,int32_t>, which compares the first members and
    // breaks ties with the second members.
    const permute: [number, number][] = [];
    for (let i = 0; i < n; ++i) {
        permute.push([-e.values[i], i]);
    }
    permute.sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[1] - b[1]));

    const invPermute = new Array<number>(n).fill(0);
    for (let i = 0; i < n; ++i) {
        invPermute[permute[i][1]] = i;
    }

    const locE = new Vector(n);
    const locY = new Vector(n);
    for (let i = 0; i < n; ++i) {
        const j = permute[i][1];
        locE.values[i] = e.values[j];
        locY.values[i] = Math.abs(y.values[j]);
    }

    const locX = new Vector(n);
    const result = sqrDistanceSpecial(locE, locY, locX);

    // Restore the axis order and reflections.
    for (let i = 0; i < n; ++i) {
        const j = invPermute[i];
        if (negate[i]) {
            locX.values[j] = -locX.values[j];
        }
        x.values[i] = locX.values[j];
    }

    return result;
}

export class DistPointHyperellipsoid
    implements DCPQuery<Vector, Hyperellipsoid, DistPointHyperellipsoidResult> {
    // The query for any hyperellipsoid.
    compute(point: Vector, hyperellipsoid: Hyperellipsoid):
        DistPointHyperellipsoidResult {
        const n = point.size;
        logAssert(n === hyperellipsoid.center.size,
            'DistPointHyperellipsoid: mismatched dimensions.');

        // Compute the coordinates of Y in the hyperellipsoid coordinate
        // system.
        const diff = sub(point, hyperellipsoid.center);
        const y = new Vector(n);
        for (let i = 0; i < n; ++i) {
            y.values[i] = dot(diff, hyperellipsoid.axis[i]);
        }

        // Compute the closest hyperellipsoid point in the axis-aligned
        // coordinate system.
        const x = new Vector(n);
        const sqrDist = sqrDistance(hyperellipsoid.extent, y, x);

        // Convert back to the original coordinate system.
        let closest1 = hyperellipsoid.center.clone();
        for (let i = 0; i < n; ++i) {
            closest1 = add(closest1, mul(x.values[i], hyperellipsoid.axis[i]));
        }

        return {
            distance: Math.sqrt(sqrDist),
            sqrDistance: sqrDist,
            closest: [point.clone(), closest1]
        };
    }

    // The hyperellipsoid is assumed to be axis-aligned and centered at the
    // origin, so only the extent[] values are used. (The port of the second
    // upstream 'operator()' overload.)
    computeAxisAligned(point: Vector, extent: Vector):
        DistPointHyperellipsoidResult {
        logAssert(point.size === extent.size,
            'DistPointHyperellipsoid: mismatched dimensions.');

        const closest1 = new Vector(point.size);
        const sqrDist = sqrDistance(extent, point, closest1);
        return {
            distance: Math.sqrt(sqrDist),
            sqrDistance: sqrDist,
            closest: [point.clone(), closest1]
        };
    }
}
