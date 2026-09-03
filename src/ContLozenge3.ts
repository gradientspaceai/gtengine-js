// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContLozenge3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for lozenges in 3D.
//
// Port notes: see ContAlignedBox.ts for the Cont* naming precedent
// (getContainerLozenge3, inContainerLozenge3, with the container type as a
// suffix so the exported names are globally unique). Upstream returns a
// vestigial 'true' from GetContainer and fills an output reference; the port
// returns the lozenge. The number-of-points argument is dropped in favor of
// the array length. Upstream copies the fitted OrientedBox3 by value and then
// mutates it; ApprGaussian3::getParameters returns the fitter's own object,
// so the port clones it.

import { ApprGaussian3 } from './ApprGaussian3.js';
import { DistPointRectangle } from './DistPointRectangle.js';
import { logAssert } from './Logger.js';
import { Lozenge3 } from './Lozenge3.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

// Compute the plane of the lozenge rectangle using least-squares fit.
// Parallel planes are chosen close enough together so that all the data
// points lie between them. The radius is half the distance between the two
// planes. The half-cylinder and quarter-cylinder side pieces are chosen
// using a method similar to that used for fitting by capsules.
//
// Upstream reads points[0] unconditionally, so at least one point is
// required.
export function getContainerLozenge3(points: readonly Vector[]): Lozenge3 {
    logAssert(points.length > 0, 'getContainerLozenge3: no points.');
    for (const point of points) {
        logAssert(point.size === 3, 'getContainerLozenge3: points must be 3D.');
    }

    const fitter = new ApprGaussian3();
    fitter.fit(points);
    const box = fitter.getParameters().clone();

    let diff = sub(points[0], box.center);
    let wMin = dot(box.axis[0], diff);
    let wMax = wMin;
    let w: number;
    for (let i = 1; i < points.length; ++i) {
        diff = sub(points[i], box.center);
        w = dot(box.axis[0], diff);
        if (w < wMin) {
            wMin = w;
        }
        else if (w > wMax) {
            wMax = w;
        }
    }

    const radius = 0.5 * (wMax - wMin);
    const rSqr = radius * radius;
    box.center = add(box.center, mul(0.5 * (wMax + wMin), box.axis[0]));

    let aMin = Number.MAX_VALUE;
    let aMax = -aMin;
    let bMin = Number.MAX_VALUE;
    let bMax = -bMin;
    let discr: number, radical: number, u: number, v: number, test: number;
    for (let i = 0; i < points.length; ++i) {
        diff = sub(points[i], box.center);
        u = dot(box.axis[2], diff);
        v = dot(box.axis[1], diff);
        w = dot(box.axis[0], diff);
        discr = rSqr - w * w;
        radical = Math.sqrt(Math.max(discr, 0));

        test = u + radical;
        if (test < aMin) {
            aMin = test;
        }

        test = u - radical;
        if (test > aMax) {
            aMax = test;
        }

        test = v + radical;
        if (test < bMin) {
            bMin = test;
        }

        test = v - radical;
        if (test > bMax) {
            bMax = test;
        }
    }

    // The enclosing region might be a capsule or a sphere.
    if (aMin >= aMax) {
        test = 0.5 * (aMin + aMax);
        aMin = test;
        aMax = test;
    }
    if (bMin >= bMax) {
        test = 0.5 * (bMin + bMax);
        bMin = test;
        bMax = test;
    }

    // Make correction for points inside mitered corner but outside quarter
    // sphere. Upstream selects the extreme values to update through the
    // pointers aExtreme and bExtreme; the port tracks which of the min/max
    // pair is selected with booleans.
    for (let i = 0; i < points.length; ++i) {
        diff = sub(points[i], box.center);
        u = dot(box.axis[2], diff);
        v = dot(box.axis[1], diff);

        let aIsMax: boolean | null = null;
        let bIsMax = false;

        if (u > aMax) {
            if (v > bMax) {
                aIsMax = true;
                bIsMax = true;
            }
            else if (v < bMin) {
                aIsMax = true;
                bIsMax = false;
            }
        }
        else if (u < aMin) {
            if (v > bMax) {
                aIsMax = false;
                bIsMax = true;
            }
            else if (v < bMin) {
                aIsMax = false;
                bIsMax = false;
            }
        }

        if (aIsMax !== null) {
            const aExtreme = aIsMax ? aMax : aMin;
            const bExtreme = bIsMax ? bMax : bMin;
            const deltaU = u - aExtreme;
            const deltaV = v - bExtreme;
            const deltaSumSqr = deltaU * deltaU + deltaV * deltaV;
            w = dot(box.axis[0], diff);
            const wSqr = w * w;
            test = deltaSumSqr + wSqr;
            if (test > rSqr) {
                discr = (rSqr - wSqr) / deltaSumSqr;
                const t = -Math.sqrt(Math.max(discr, 0));
                const newA = u + t * deltaU;
                const newB = v + t * deltaV;
                if (aIsMax) {
                    aMax = newA;
                }
                else {
                    aMin = newA;
                }
                if (bIsMax) {
                    bMax = newB;
                }
                else {
                    bMin = newB;
                }
            }
        }
    }

    const lozenge = new Lozenge3();
    lozenge.radius = radius;
    lozenge.rectangle.axis[0] = box.axis[2].clone();
    lozenge.rectangle.axis[1] = box.axis[1].clone();

    // UPSTREAM BUG (fixed here): in the branches where a direction is
    // non-degenerate, upstream offsets the rectangle center by aMin (and/or
    // bMin) instead of the interval midpoint 0.5*(aMin+aMax) (and/or
    // 0.5*(bMin+bMax)), while still setting the extent to the half-width
    // 0.5*(aMax-aMin). Rectangle<3,Real> is a centered primitive
    // (R(s0,s1) = C + s0*A0 + s1*A1 with |si| <= ei), so the upstream
    // rectangle spans [aMin-halfWidth, aMin+halfWidth] rather than the
    // intended [aMin, aMax] and the returned lozenge fails to contain the
    // input points. (The offsets read like a leftover from Wild Magic's
    // corner-origin rectangle; upstream's own comments in Rectangle.h and
    // Lozenge3.h still call the center an 'origin'.) The port uses the
    // midpoint in every branch, which also makes the degenerate branches
    // agree with the upstream sphere branch.
    if (aMin < aMax) {
        if (bMin < bMax) {
            // Container is a lozenge.
            lozenge.rectangle.extent.values[0] = 0.5 * (aMax - aMin);
            lozenge.rectangle.extent.values[1] = 0.5 * (bMax - bMin);
        }
        else {
            // Container is a capsule.
            lozenge.rectangle.extent.values[0] = 0.5 * (aMax - aMin);
            lozenge.rectangle.extent.values[1] = 0;
        }
    }
    else {
        if (bMin < bMax) {
            // Container is a capsule.
            lozenge.rectangle.extent.values[0] = 0;
            lozenge.rectangle.extent.values[1] = 0.5 * (bMax - bMin);
        }
        else {
            // Container is a sphere.
            lozenge.rectangle.extent.values[0] = 0;
            lozenge.rectangle.extent.values[1] = 0;
        }
    }

    lozenge.rectangle.center = add(box.center,
        add(mul(0.5 * (aMin + aMax), box.axis[2]),
            mul(0.5 * (bMin + bMax), box.axis[1])));

    return lozenge;
}

// Test for containment of a point by a lozenge.
export function inContainerLozenge3(point: Vector, lozenge: Lozenge3): boolean {
    const prQuery = new DistPointRectangle();
    const result = prQuery.compute(point, lozenge.rectangle);
    return result.distance <= lozenge.radius;
}
