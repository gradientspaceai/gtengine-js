// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContCapsule3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for capsules in 3D.
//
// Port notes: see ContAlignedBox.ts for the Cont* naming precedent (each
// GetContainer/InContainer/MergeContainers overload is suffixed with the
// container type so the exported names are globally unique). Upstream has
// three InContainer overloads for Capsule3, distinguished by the type of the
// contained object; the port disambiguates by naming the contained type
// before the container suffix:
//   InContainer(point, capsule)       -> inContainerCapsule3
//   InContainer(sphere, capsule)      -> inContainerSphereCapsule3
//   InContainer(testCapsule, capsule) -> inContainerCapsuleCapsule3
// Upstream returns a vestigial 'true' from GetContainer/MergeContainers and
// fills an output reference; the port returns the capsule. The
// number-of-points argument is dropped in favor of the array length.

import { ApprOrthogonalLine3 } from './ApprOrthogonalLine3';
import { Capsule, type Capsule3 } from './Capsule';
import { DistPointLine } from './DistPointLine';
import { DistPointSegment } from './DistPointSegment';
import { Hypersphere, type Sphere3 } from './Hypersphere';
import { Line } from './Line';
import { logAssert } from './Logger';
import { Segment } from './Segment';
import { Vector, add, dot, mul, normalize, sub } from './Vector';
import { computeOrthogonalComplement3 } from './Vector3';

// Compute the axis of the capsule segment using least-squares fitting. The
// radius is the maximum distance from the points to the axis. Hemispherical
// caps are chosen as close together as possible.
//
// Upstream reads points[0] unconditionally, so at least one point is
// required.
export function getContainerCapsule3(points: readonly Vector[]): Capsule3 {
    logAssert(points.length > 0, 'getContainerCapsule3: no points.');
    for (const point of points) {
        logAssert(point.size === 3, 'getContainerCapsule3: points must be 3D.');
    }

    const fitter = new ApprOrthogonalLine3();
    fitter.fit(points);
    const line = fitter.getParameters();

    const plQuery = new DistPointLine();
    let maxRadiusSqr = 0;
    for (let i = 0; i < points.length; ++i) {
        const result = plQuery.compute(points[i], line);
        if (result.sqrDistance > maxRadiusSqr) {
            maxRadiusSqr = result.sqrDistance;
        }
    }

    const basis: Vector[] = [line.direction.clone(), new Vector(3), new Vector(3)];
    computeOrthogonalComplement3(1, basis);

    let minValue = Number.MAX_VALUE;
    let maxValue = -Number.MAX_VALUE;
    for (let i = 0; i < points.length; ++i) {
        const diff = sub(points[i], line.origin);
        const uDotDiff = dot(diff, basis[1]);
        const vDotDiff = dot(diff, basis[2]);
        const wDotDiff = dot(diff, basis[0]);
        const discr = maxRadiusSqr - (uDotDiff * uDotDiff + vDotDiff * vDotDiff);
        const radical = Math.sqrt(Math.max(discr, 0));

        let test = wDotDiff + radical;
        if (test < minValue) {
            minValue = test;
        }

        test = wDotDiff - radical;
        if (test > maxValue) {
            maxValue = test;
        }
    }

    const center = add(line.origin, mul(0.5 * (minValue + maxValue), line.direction));

    let extent: number;
    if (maxValue > minValue) {
        // Container is a capsule.
        extent = 0.5 * (maxValue - minValue);
    }
    else {
        // Container is a sphere.
        extent = 0;
    }

    const capsule = new Capsule(3);
    capsule.segment = Segment.fromCenteredForm(center, line.direction, extent);
    capsule.radius = Math.sqrt(maxRadiusSqr);
    return capsule;
}

// Test for containment of a point by a capsule.
export function inContainerCapsule3(point: Vector, capsule: Capsule3): boolean {
    const psQuery = new DistPointSegment();
    const result = psQuery.compute(point, capsule.segment);
    return result.distance <= capsule.radius;
}

// Test for containment of a sphere by a capsule.
export function inContainerSphereCapsule3(sphere: Sphere3,
    capsule: Capsule3): boolean {
    const rDiff = capsule.radius - sphere.radius;
    if (rDiff >= 0) {
        const psQuery = new DistPointSegment();
        const result = psQuery.compute(sphere.center, capsule.segment);
        return result.distance <= rDiff;
    }
    return false;
}

// Test for containment of a capsule by a capsule.
export function inContainerCapsuleCapsule3(testCapsule: Capsule3,
    capsule: Capsule3): boolean {
    const spherePosEnd = Hypersphere.fromCenterRadius(
        testCapsule.segment.p[1], testCapsule.radius);
    const sphereNegEnd = Hypersphere.fromCenterRadius(
        testCapsule.segment.p[0], testCapsule.radius);
    return inContainerSphereCapsule3(spherePosEnd, capsule)
        && inContainerSphereCapsule3(sphereNegEnd, capsule);
}

// Compute a capsule that contains the input capsules. The returned capsule
// is not necessarily the one of smallest volume that contains the inputs.
export function mergeContainersCapsule3(capsule0: Capsule3,
    capsule1: Capsule3): Capsule3 {
    if (inContainerCapsuleCapsule3(capsule0, capsule1)) {
        return capsule1.clone();
    }

    if (inContainerCapsuleCapsule3(capsule1, capsule0)) {
        return capsule0.clone();
    }

    const cf0 = capsule0.segment.getCenteredForm();
    const cf1 = capsule1.segment.getCenteredForm();
    const P0 = cf0.center, D0 = cf0.direction;
    const P1 = cf1.center, D1 = cf1.direction;

    // Axis of final capsule.
    const line = new Line(3);

    // Axis center is average of input axis centers.
    line.origin = mul(0.5, add(P0, P1));

    // Axis unit direction is average of input axis unit directions.
    if (dot(D0, D1) >= 0) {
        line.direction = add(D0, D1);
    }
    else {
        line.direction = sub(D0, D1);
    }
    normalize(line.direction);

    // Cylinder with axis 'line' must contain the spheres centered at the
    // endpoints of the input capsules.
    const plQuery = new DistPointLine();
    const posEnd0 = capsule0.segment.p[1];
    let radius = plQuery.compute(posEnd0, line).distance + capsule0.radius;

    const negEnd0 = capsule0.segment.p[0];
    let tmp = plQuery.compute(negEnd0, line).distance + capsule0.radius;
    if (tmp > radius) {
        radius = tmp;
    }

    const posEnd1 = capsule1.segment.p[1];
    tmp = plQuery.compute(posEnd1, line).distance + capsule1.radius;
    if (tmp > radius) {
        radius = tmp;
    }

    const negEnd1 = capsule1.segment.p[0];
    tmp = plQuery.compute(negEnd1, line).distance + capsule1.radius;
    if (tmp > radius) {
        radius = tmp;
    }

    // In the following blocks of code, theoretically k1*k1-k0 >= 0, but
    // numerical rounding errors can make it slightly negative. Guard against
    // this.

    // Process sphere <posEnd0,r0>.
    let rDiff = radius - capsule0.radius;
    let rDiffSqr = rDiff * rDiff;
    let diff = sub(line.origin, posEnd0);
    let k0 = dot(diff, diff) - rDiffSqr;
    let k1 = dot(diff, line.direction);
    let discr = k1 * k1 - k0;
    let root = Math.sqrt(Math.max(discr, 0));
    let tPos = -k1 - root;
    let tNeg = -k1 + root;

    // Process sphere <negEnd0,r0>.
    diff = sub(line.origin, negEnd0);
    k0 = dot(diff, diff) - rDiffSqr;
    k1 = dot(diff, line.direction);
    discr = k1 * k1 - k0;
    root = Math.sqrt(Math.max(discr, 0));
    tmp = -k1 - root;
    if (tmp > tPos) {
        tPos = tmp;
    }
    tmp = -k1 + root;
    if (tmp < tNeg) {
        tNeg = tmp;
    }

    // Process sphere <posEnd1,r1>.
    rDiff = radius - capsule1.radius;
    rDiffSqr = rDiff * rDiff;
    diff = sub(line.origin, posEnd1);
    k0 = dot(diff, diff) - rDiffSqr;
    k1 = dot(diff, line.direction);
    discr = k1 * k1 - k0;
    root = Math.sqrt(Math.max(discr, 0));
    tmp = -k1 - root;
    if (tmp > tPos) {
        tPos = tmp;
    }
    tmp = -k1 + root;
    if (tmp < tNeg) {
        tNeg = tmp;
    }

    // Process sphere <negEnd1,r1>.
    diff = sub(line.origin, negEnd1);
    k0 = dot(diff, diff) - rDiffSqr;
    k1 = dot(diff, line.direction);
    discr = k1 * k1 - k0;
    root = Math.sqrt(Math.max(discr, 0));
    tmp = -k1 - root;
    if (tmp > tPos) {
        tPos = tmp;
    }
    tmp = -k1 + root;
    if (tmp < tNeg) {
        tNeg = tmp;
    }

    const center = add(line.origin, mul(0.5 * (tPos + tNeg), line.direction));

    let extent: number;
    if (tPos > tNeg) {
        // Container is a capsule.
        extent = 0.5 * (tPos - tNeg);
    }
    else {
        // Container is a sphere.
        extent = 0;
    }

    const merge = new Capsule(3);
    merge.segment = Segment.fromCenteredForm(center, line.direction, extent);
    merge.radius = radius;
    return merge;
}
