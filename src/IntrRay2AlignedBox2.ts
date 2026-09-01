// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2AlignedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box to be a solid.
//
// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the four
// edges of the box.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): upstream derives
// the ray queries from the line queries only to reuse the protected DoQuery
// members. In TypeScript the derived query cannot keep the canonical
// test()/find() names while changing the first parameter type, so the port
// reuses the line algorithm through a module-private subclass that exposes
// DoQuery. The upstream Result structs add no members to the line results, so
// the port exports type aliases. The ray-specific DoQuery members remain
// available as the protected doQuery() methods.

import type { AlignedBox } from './AlignedBox';
import type { Ray } from './Ray';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import { Vector, add, mul, sub } from './Vector';
import {
    IntrLine2AlignedBox2TI,
    IntrLine2AlignedBox2FI,
    defaultIntrLine2AlignedBox2TIResult,
    defaultIntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2';
import type {
    IntrLine2AlignedBox2TIResult,
    IntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2';
import { IntrIntervalsFI } from './IntrIntervals';

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrRay2AlignedBox2TI.test. Upstream adds no members to the
// line-box result.
export type IntrRay2AlignedBox2TIResult = IntrLine2AlignedBox2TIResult;

// The result of IntrRay2AlignedBox2FI.find. Upstream adds no members to the
// line-box result.
export type IntrRay2AlignedBox2FIResult = IntrLine2AlignedBox2FIResult;

// Accessors for the protected line-versus-aligned-box DoQuery members.
class LineBoxTIAccess extends IntrLine2AlignedBox2TI {
    run(origin: Vector, direction: Vector, extent: Vector,
        result: IntrLine2AlignedBox2TIResult): void {
        this.doQuery(origin, direction, extent, result);
    }
}

class LineBoxFIAccess extends IntrLine2AlignedBox2FI {
    run(origin: Vector, direction: Vector, extent: Vector,
        result: IntrLine2AlignedBox2FIResult): void {
        this.doQuery(origin, direction, extent, result);
    }
}

// Test-intersection query for a ray and a solid aligned box in 2D.
export class IntrRay2AlignedBox2TI implements
    TIQuery<Ray, AlignedBox, IntrRay2AlignedBox2TIResult> {

    private readonly base = new LineBoxTIAccess();

    test(ray: Ray, box: AlignedBox): IntrRay2AlignedBox2TIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector2::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the ray to the aligned-box coordinate system.
        const rayOrigin = sub(ray.origin, boxCenter);

        const result = defaultIntrLine2AlignedBox2TIResult();
        this.doQuery(rayOrigin, ray.direction, boxExtent, result);
        return result;
    }

    protected doQuery(rayOrigin: Vector, rayDirection: Vector,
        boxExtent: Vector, result: IntrRay2AlignedBox2TIResult): void {
        const o = rayOrigin.values;
        const d = rayDirection.values;
        const e = boxExtent.values;
        for (let i = 0; i < 2; ++i) {
            if (Math.abs(o[i]) > e[i] && o[i] * d[i] >= 0) {
                result.intersect = false;
                return;
            }
        }

        this.base.run(rayOrigin, rayDirection, boxExtent, result);
    }
}

// Find-intersection query for a ray and a solid aligned box in 2D.
export class IntrRay2AlignedBox2FI implements
    FIQuery<Ray, AlignedBox, IntrRay2AlignedBox2FIResult> {

    private readonly base = new LineBoxFIAccess();

    find(ray: Ray, box: AlignedBox): IntrRay2AlignedBox2FIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector2::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the ray to the aligned-box coordinate system.
        const rayOrigin = sub(ray.origin, boxCenter);

        const result = defaultIntrLine2AlignedBox2FIResult();
        this.doQuery(rayOrigin, ray.direction, boxExtent, result);
        for (let i = 0; i < result.numIntersections; ++i) {
            result.point[i] = add(ray.origin,
                mul(result.parameter[i], ray.direction));
        }
        return result;
    }

    protected doQuery(rayOrigin: Vector, rayDirection: Vector,
        boxExtent: Vector, result: IntrRay2AlignedBox2FIResult): void {
        this.base.run(rayOrigin, rayDirection, boxExtent, result);

        if (result.intersect) {
            // The line containing the ray intersects the box; the t-interval
            // is [t0,t1]. The ray intersects the box as long as [t0,t1]
            // overlaps the ray t-interval [0,+infinity).
            const rayInterval: [number, number] = [0, MAX_T];
            const iiQuery = new IntrIntervalsFI();
            const iiResult = iiQuery.find(result.parameter, rayInterval);
            result.intersect = iiResult.intersect;
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
    }
}
