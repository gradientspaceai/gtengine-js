// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment3CanonicalBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a segment and a solid canonical box in 3D.
//
// The segment is P0 + t * (P1 - P0) for 0 <= t <= 1. The direction D = P1-P0
// is generally not unit length.
//
// The canonical box has center at the origin and is aligned with the
// coordinate axes. The extents are E = (e[0],e[1],e[2]). A box point is
// Y = (y[0],y[1],y[2]) with |y[i]| <= e[i] for all i.
//
// The closest point on the segment is stored in closest[0] with parameter t.
// The closest point on the box is stored in closest[1]. When there are
// infinitely many choices for the pair of closest points, only one of them is
// returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment3<T>, CanonicalBox3<T>>'
// becomes the class DistSegment3CanonicalBox3. As upstream does, the result
// type is the line-box result type, re-exported here as the alias
// DistSegment3CanonicalBox3Result.

import type { CanonicalBox3 } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistLine3CanonicalBox3 } from './DistLine3CanonicalBox3';
import type { DistLine3CanonicalBox3Result } from './DistLine3CanonicalBox3';
import { DistPointCanonicalBox } from './DistPointCanonicalBox';
import { Line } from './Line';
import type { Segment3 } from './Segment';
import { Vector, sub } from './Vector';

// Upstream reuses the line-box result type ('using Result = typename
// LBQuery::Result').
export type DistSegment3CanonicalBox3Result = DistLine3CanonicalBox3Result;

export class DistSegment3CanonicalBox3
    implements DCPQuery<Segment3, CanonicalBox3,
    DistSegment3CanonicalBox3Result> {
    compute(segment: Segment3, box: CanonicalBox3):
        DistSegment3CanonicalBox3Result {
        const result: DistSegment3CanonicalBox3Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            closest: [new Vector(3), new Vector(3)]
        };

        // Compute the distance between the box and the line containing the
        // segment. If the closest line point has a parameter in [0,1], it is
        // on the segment and the line result is the segment result. Otherwise
        // the squared distance is a convex function of the line parameter, so
        // clamping the parameter to the nearest endpoint and re-solving as a
        // point-box query gives the segment result.
        const segDirection = sub(segment.p[1], segment.p[0]);
        const line = Line.fromOriginDirection(segment.p[0], segDirection);
        const lbQuery = new DistLine3CanonicalBox3();
        const lbOutput = lbQuery.compute(line, box);
        if (lbOutput.parameter >= 0) {
            if (lbOutput.parameter <= 1) {
                return lbOutput;
            }

            const pbQuery = new DistPointCanonicalBox();
            const pbOutput = pbQuery.compute(segment.p[1], box);
            result.sqrDistance = pbOutput.sqrDistance;
            result.distance = pbOutput.distance;
            result.parameter = 1;
            result.closest[0] = segment.p[1].clone();
            result.closest[1] = pbOutput.closest[1];
        }
        else {
            const pbQuery = new DistPointCanonicalBox();
            const pbOutput = pbQuery.compute(segment.p[0], box);
            result.sqrDistance = pbOutput.sqrDistance;
            result.distance = pbOutput.distance;
            result.parameter = 0;
            result.closest[0] = segment.p[0].clone();
            result.closest[1] = pbOutput.closest[1];
        }
        return result;
    }
}
