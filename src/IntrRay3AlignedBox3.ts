// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the six
// faces of the box (Liang-Barsky clipping). The queries consider the box to
// be a solid. The algorithms are described in
// https://www.geometrictools.com/Documentation/IntersectionLineBox.pdf
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
    IntrLine3AlignedBox3TI,
    IntrLine3AlignedBox3FI,
    defaultIntrLine3AlignedBox3TIResult,
    defaultIntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3';
import type {
    IntrLine3AlignedBox3TIResult,
    IntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3';
import { IntrIntervalsFI } from './IntrIntervals';

// The result of IntrRay3AlignedBox3TI.test. Upstream adds no members to the
// line-box result.
export type IntrRay3AlignedBox3TIResult = IntrLine3AlignedBox3TIResult;

// The result of IntrRay3AlignedBox3FI.find. Upstream adds no members to the
// line-box result.
export type IntrRay3AlignedBox3FIResult = IntrLine3AlignedBox3FIResult;

// Accessors for the protected line-versus-aligned-box DoQuery members.
class LineBoxTIAccess extends IntrLine3AlignedBox3TI {
    run(origin: Vector, direction: Vector, extent: Vector,
        result: IntrLine3AlignedBox3TIResult): void {
        this.doQuery(origin, direction, extent, result);
    }
}

class LineBoxFIAccess extends IntrLine3AlignedBox3FI {
    run(origin: Vector, direction: Vector, extent: Vector,
        result: IntrLine3AlignedBox3FIResult): void {
        this.doQuery(origin, direction, extent, result);
    }
}

// Test-intersection query for a ray and a solid aligned box in 3D.
export class IntrRay3AlignedBox3TI implements
    TIQuery<Ray, AlignedBox, IntrRay3AlignedBox3TIResult> {

    private readonly base = new LineBoxTIAccess();

    test(ray: Ray, box: AlignedBox): IntrRay3AlignedBox3TIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector3::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the ray to the aligned-box coordinate system.
        const rayOrigin = sub(ray.origin, boxCenter);

        const result = defaultIntrLine3AlignedBox3TIResult();
        this.doQuery(rayOrigin, ray.direction, boxExtent, result);
        return result;
    }

    // The caller must ensure that on entry, 'result' is default constructed
    // as if there is no intersection. If an intersection is found, the
    // 'result' values are modified accordingly.
    protected doQuery(rayOrigin: Vector, rayDirection: Vector,
        boxExtent: Vector, result: IntrRay3AlignedBox3TIResult): void {
        const o = rayOrigin.values;
        const d = rayDirection.values;
        const e = boxExtent.values;
        for (let i = 0; i < 3; ++i) {
            if (Math.abs(o[i]) > e[i] && o[i] * d[i] >= 0) {
                result.intersect = false;
                return;
            }
        }

        this.base.run(rayOrigin, rayDirection, boxExtent, result);
    }
}

// Find-intersection query for a ray and a solid aligned box in 3D.
export class IntrRay3AlignedBox3FI implements
    FIQuery<Ray, AlignedBox, IntrRay3AlignedBox3FIResult> {

    private readonly base = new LineBoxFIAccess();

    find(ray: Ray, box: AlignedBox): IntrRay3AlignedBox3FIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector3::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the ray to the aligned-box coordinate system.
        const rayOrigin = sub(ray.origin, boxCenter);

        const result = defaultIntrLine3AlignedBox3FIResult();
        this.doQuery(rayOrigin, ray.direction, boxExtent, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(ray.origin,
                    mul(result.parameter[i], ray.direction));
            }
        }
        return result;
    }

    // The caller must ensure that on entry, 'result' is default constructed
    // as if there is no intersection. If an intersection is found, the
    // 'result' values are modified accordingly.
    protected doQuery(rayOrigin: Vector, rayDirection: Vector,
        boxExtent: Vector, result: IntrRay3AlignedBox3FIResult): void {
        this.base.run(rayOrigin, rayDirection, boxExtent, result);

        if (result.intersect) {
            // The line containing the ray intersects the box; the t-interval
            // is [t0,t1]. The ray intersects the box as long as [t0,t1]
            // overlaps the ray t-interval [0,+infinity).
            const iiQuery = new IntrIntervalsFI();
            const iiResult = iiQuery.findFiniteSemiInfinite(result.parameter,
                0, true);
            if (iiResult.intersect) {
                result.numIntersections = iiResult.numIntersections;
                result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
            }
            else {
                // The ray does not intersect the box.
                result.intersect = false;
                result.numIntersections = 0;
                result.parameter = [0, 0];
                result.point = [Vector.zero(3), Vector.zero(3)];
            }
        }
    }
}
