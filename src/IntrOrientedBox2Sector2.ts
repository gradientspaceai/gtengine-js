// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox2Sector2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The OrientedBox2 object is considered to be a solid.
//
// The test first checks whether the sector vertex is inside the box. If not,
// the box is tested against the two ray boundaries of the sector using the
// method of separating axes. If the box is not separated, it is clipped
// against the ray boundaries and the clipped polygon is tested against the
// disk of the sector radius centered at the sector vertex.
//
// NOTE (upstream limitation, preserved): the two boundary rays are used as
// halfplanes and the wedge is treated as their intersection. That is valid
// only when the sector half-angle is at most pi/2; for a larger half-angle
// the sector is not convex and is strictly larger than the intersection of
// the two halfplanes, so the query can report "no intersection" for a box
// that does meet the sector. In particular, the Sector2 default half-angle
// of pi describes a full disk but is handled as a halfplane wedge.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only a TIQuery specialization for this pair of primitives, which becomes
// IntrOrientedBox2Sector2TI.

import type { OrientedBox } from './OrientedBox';
import type { Sector2 } from './Sector2';
import type { TIQuery } from './TIQuery';
import { Vector, add, dot, mul, sub, negate } from './Vector';
import { perp } from './Vector2';
import { Halfspace } from './Halfspace';
import { Segment } from './Segment';
import { IntrHalfspace2Polygon2FI } from './IntrHalfspace2Polygon2';
import { DistPointSegment } from './DistPointSegment';

// The result of IntrOrientedBox2Sector2TI queries.
export interface IntrOrientedBox2Sector2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrOrientedBox2Sector2TIResult {
    return { intersect: false };
}

// Clip a convex polygon against a halfspace.
//
// Upstream bug (fixed here): IntrOrientedBox2Sector2.h assigns the clip
// result's polygon unconditionally. IntrHalfspace2Polygon2's FIQuery reports
// intersect = true with an EMPTY polygon when no clipping is necessary
// (the polygon already lies in the closed halfspace), so the unconditional
// assignment throws the polygon away and the sector query then reports "no
// intersection". This is reachable in a whole region of configurations, not
// just at a boundary: the second clip decides whether it must run from the
// projection of the original box, so the polygon produced by the first clip
// can easily lie entirely inside the second halfspace. The port keeps the
// current polygon in that case.
function clip(hpQuery: IntrHalfspace2Polygon2FI, halfspace: Halfspace,
    polygon: Vector[]): Vector[] {
    const hpResult = hpQuery.find(halfspace, polygon);
    if (hpResult.intersect && hpResult.polygon.length === 0) {
        // No clipping was necessary; keep the polygon unchanged.
        return polygon;
    }
    return hpResult.polygon;
}

// Test-intersection query for a solid oriented box and a sector in 2D.
export class IntrOrientedBox2Sector2TI implements
    TIQuery<OrientedBox, Sector2, IntrOrientedBox2Sector2TIResult> {

    test(box: OrientedBox, sector: Sector2): IntrOrientedBox2Sector2TIResult {
        const result = defaultTIResult();

        // Determine whether the vertex is inside the box.
        const CmV = sub(box.center, sector.vertex);
        const P = [dot(box.axis[0], CmV), dot(box.axis[1], CmV)];
        if (Math.abs(P[0]) <= box.extent.values[0] &&
            Math.abs(P[1]) <= box.extent.values[1]) {
            // The vertex is inside the box.
            result.intersect = true;
            return result;
        }

        const d = sector.direction.values;

        // Test whether the box is outside the right ray boundary of the
        // sector.
        const U0 = Vector.fromArray([
            +sector.cosAngle * d[0] + sector.sinAngle * d[1],
            -sector.sinAngle * d[0] + sector.cosAngle * d[1]
        ]);
        const N0 = perp(U0);
        const prjcen0 = dot(N0, CmV);
        const radius0 = box.extent.values[0] * Math.abs(dot(N0, box.axis[0]))
            + box.extent.values[1] * Math.abs(dot(N0, box.axis[1]));
        if (prjcen0 > radius0) {
            result.intersect = false;
            return result;
        }

        // Test whether the box is outside the ray of the left boundary of the
        // sector.
        const U1 = Vector.fromArray([
            +sector.cosAngle * d[0] - sector.sinAngle * d[1],
            +sector.sinAngle * d[0] + sector.cosAngle * d[1]
        ]);
        const N1 = negate(perp(U1));
        const prjcen1 = dot(N1, CmV);
        const radius1 = box.extent.values[0] * Math.abs(dot(N1, box.axis[0]))
            + box.extent.values[1] * Math.abs(dot(N1, box.axis[1]));
        if (prjcen1 > radius1) {
            result.intersect = false;
            return result;
        }

        // Initialize the polygon of intersection to be the box.
        const e0U0 = mul(box.extent.values[0], box.axis[0]);
        const e1U1 = mul(box.extent.values[1], box.axis[1]);
        let polygon: Vector[] = [
            sub(sub(box.center, e0U0), e1U1),
            sub(add(box.center, e0U0), e1U1),
            add(add(box.center, e0U0), e1U1),
            add(sub(box.center, e0U0), e1U1)
        ];

        const hpQuery = new IntrHalfspace2Polygon2FI();

        // Clip the box against the right-ray sector boundary.
        if (prjcen0 >= -radius0) {
            const normal = negate(N0);
            const halfspace = Halfspace.fromNormalConstant(normal,
                dot(normal, sector.vertex));
            polygon = clip(hpQuery, halfspace, polygon);
        }

        // Clip the box against the left-ray sector boundary.
        if (prjcen1 >= -radius1) {
            const normal = negate(N1);
            const halfspace = Halfspace.fromNormalConstant(normal,
                dot(normal, sector.vertex));
            polygon = clip(hpQuery, halfspace, polygon);
        }

        const psQuery = new DistPointSegment();
        const numVertices = polygon.length;
        if (numVertices >= 2) {
            for (let i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++) {
                const segment = Segment.fromEndpoints(polygon[i0], polygon[i1]);
                const psResult = psQuery.compute(sector.vertex, segment);
                if (psResult.distance <= sector.radius) {
                    result.intersect = true;
                    return result;
                }
            }
        }

        result.intersect = false;
        return result;
    }
}
