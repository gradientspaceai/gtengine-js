// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrArc2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Find-intersection query for two arcs.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// an FIQuery specialization, so the port has only IntrArc2Arc2FI. The nested
// 'Configuration' enumeration becomes the file-qualified exported enum
// IntrArc2Arc2Configuration. The single-argument upstream Arc2::Contains
// (which assumes the point is on the circle of the arc) is
// Arc2.containsOnCircle in the port, and C++ Vector2 equality becomes
// Vector.equals.

import { Arc2 } from './Arc2.js';
import { Hypersphere } from './Hypersphere.js';
import { Vector } from './Vector.js';
import { IntrCircle2Circle2FI } from './IntrCircle2Circle2.js';
import type { FIQuery } from './FIQuery.js';

// The port of std::numeric_limits<int32_t>::max(), which the circle-circle
// query uses as the 'numIntersections' value meaning "the circles are the
// same".
const INT32_MAX = 2147483647;

// The possible configurations of the intersection set. The comments list the
// valid result members for each configuration.
export enum IntrArc2Arc2Configuration {
    NO_INTERSECTION,
    NONCOCIRCULAR_ONE_POINT,        // point[0]
    NONCOCIRCULAR_TWO_POINTS,       // point[0], point[1]
    COCIRCULAR_ONE_POINT,           // point[0]
    COCIRCULAR_TWO_POINTS,          // point[0], point[1]
    COCIRCULAR_ONE_POINT_ONE_ARC,   // point[0], arc[0]
    COCIRCULAR_ONE_ARC,             // arc[0]
    COCIRCULAR_TWO_ARCS             // arc[0], arc[1]
}

// The result of IntrArc2Arc2FI.find.
export interface IntrArc2Arc2FIResult {
    // True if and only if configuration is not NO_INTERSECTION.
    intersect: boolean;

    // One of the enumerations listed previously.
    configuration: IntrArc2Arc2Configuration;

    point: [Vector, Vector];
    arc: [Arc2, Arc2];
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrArc2Arc2FIResult {
    const zeroArc = (): Arc2 => Arc2.fromCenterRadiusEnds(Vector.zero(2), 0,
        Vector.zero(2), Vector.zero(2));
    return {
        intersect: false,
        configuration: IntrArc2Arc2Configuration.NO_INTERSECTION,
        point: [Vector.zero(2), Vector.zero(2)],
        arc: [zeroArc(), zeroArc()]
    };
}

export class IntrArc2Arc2FI implements
    FIQuery<Arc2, Arc2, IntrArc2Arc2FIResult> {

    find(arc0: Arc2, arc1: Arc2): IntrArc2Arc2FIResult {
        // Assume initially there are no intersections. If we find at least
        // one intersection, we will set result.intersect to true.
        const result = defaultFIResult();

        const circle0 = Hypersphere.fromCenterRadius(arc0.center, arc0.radius);
        const circle1 = Hypersphere.fromCenterRadius(arc1.center, arc1.radius);
        const ccQuery = new IntrCircle2Circle2FI();
        const ccResult = ccQuery.find(circle0, circle1);
        if (!ccResult.intersect) {
            // The arcs do not intersect.
            result.configuration = IntrArc2Arc2Configuration.NO_INTERSECTION;
            return result;
        }

        if (ccResult.numIntersections === INT32_MAX) {
            // The arcs are cocircular. Determine whether they overlap. Let
            // arc0 be <A0,A1> and arc1 be <B0,B1>. The points are ordered
            // counterclockwise around the circle of the arc.
            if (arc1.containsOnCircle(arc0.end[0])) {
                result.intersect = true;
                if (arc1.containsOnCircle(arc0.end[1])) {
                    if (arc0.containsOnCircle(arc1.end[0]) &&
                        arc0.containsOnCircle(arc1.end[1])) {
                        if (arc0.end[0].equals(arc1.end[0]) &&
                            arc0.end[1].equals(arc1.end[1])) {
                            // The arcs are the same.
                            result.configuration =
                                IntrArc2Arc2Configuration.COCIRCULAR_ONE_ARC;
                            result.arc[0] = arc0.clone();
                        } else {
                            // arc0 and arc1 overlap in two disjoint subsets.
                            if (!arc0.end[0].equals(arc1.end[1])) {
                                if (!arc1.end[0].equals(arc0.end[1])) {
                                    // The arcs overlap in two disjoint
                                    // subarcs, each of positive subtended
                                    // angle: <A0,B1>, <A1,B0>
                                    result.configuration =
                                        IntrArc2Arc2Configuration.COCIRCULAR_TWO_ARCS;
                                    result.arc[0] = Arc2.fromCenterRadiusEnds(
                                        arc0.center, arc0.radius, arc0.end[0],
                                        arc1.end[1]);
                                    result.arc[1] = Arc2.fromCenterRadiusEnds(
                                        arc0.center, arc0.radius, arc1.end[0],
                                        arc0.end[1]);
                                } else {
                                    // B0 = A1. The intersection is a point
                                    // {A1} and an arc <A0,B1>.
                                    result.configuration =
                                        IntrArc2Arc2Configuration.COCIRCULAR_ONE_POINT_ONE_ARC;
                                    result.point[0] = arc0.end[1].clone();
                                    result.arc[0] = Arc2.fromCenterRadiusEnds(
                                        arc0.center, arc0.radius, arc0.end[0],
                                        arc1.end[1]);
                                }
                            } else {
                                // A0 = B1.
                                if (!arc1.end[0].equals(arc0.end[1])) {
                                    // The intersection is a point {A0} and an
                                    // arc <A1,B0>.
                                    result.configuration =
                                        IntrArc2Arc2Configuration.COCIRCULAR_ONE_POINT_ONE_ARC;
                                    result.point[0] = arc0.end[0].clone();
                                    result.arc[0] = Arc2.fromCenterRadiusEnds(
                                        arc0.center, arc0.radius, arc1.end[0],
                                        arc0.end[1]);
                                } else {
                                    // The arcs share endpoints, so the union
                                    // is a circle.
                                    result.configuration =
                                        IntrArc2Arc2Configuration.COCIRCULAR_TWO_POINTS;
                                    result.point[0] = arc0.end[0].clone();
                                    result.point[1] = arc0.end[1].clone();
                                }
                            }
                        }
                    } else {
                        // Arc0 inside arc1, <B0,A0,A1,B1>.
                        result.configuration =
                            IntrArc2Arc2Configuration.COCIRCULAR_ONE_ARC;
                        result.arc[0] = arc0.clone();
                    }
                } else {
                    if (!arc0.end[0].equals(arc1.end[1])) {
                        // Arc0 and arc1 overlap, <B0,A0,B1,A1>.
                        result.configuration =
                            IntrArc2Arc2Configuration.COCIRCULAR_ONE_ARC;
                        result.arc[0] = Arc2.fromCenterRadiusEnds(arc0.center,
                            arc0.radius, arc0.end[0], arc1.end[1]);
                    } else {
                        // Arc0 and arc1 share endpoint, <B0,A0,B1,A1> with
                        // A0 = B1.
                        result.configuration =
                            IntrArc2Arc2Configuration.COCIRCULAR_ONE_POINT;
                        result.point[0] = arc0.end[0].clone();
                    }
                }
                return result;
            }

            if (arc1.containsOnCircle(arc0.end[1])) {
                result.intersect = true;
                if (!arc0.end[1].equals(arc1.end[0])) {
                    // Arc0 and arc1 overlap in a single arc, <A0,B0,A1,B1>.
                    result.configuration =
                        IntrArc2Arc2Configuration.COCIRCULAR_ONE_ARC;
                    result.arc[0] = Arc2.fromCenterRadiusEnds(arc0.center,
                        arc0.radius, arc1.end[0], arc0.end[1]);
                } else {
                    // Arc0 and arc1 share endpoint, <A0,B0,A1,B1> with
                    // B0 = A1.
                    result.configuration =
                        IntrArc2Arc2Configuration.COCIRCULAR_ONE_POINT;
                    result.point[0] = arc1.end[0].clone();
                }
                return result;
            }

            if (arc0.containsOnCircle(arc1.end[0])) {
                // Arc1 inside arc0, <A0,B0,B1,A1>.
                result.intersect = true;
                result.configuration =
                    IntrArc2Arc2Configuration.COCIRCULAR_ONE_ARC;
                result.arc[0] = arc1.clone();
            } else {
                // Arcs do not overlap, <A0,A1,B0,B1>.
                result.configuration =
                    IntrArc2Arc2Configuration.NO_INTERSECTION;
            }
            return result;
        }

        // Test whether circle-circle intersection points are on the arcs.
        let numIntersections = 0;
        for (let i = 0; i < ccResult.numIntersections; ++i) {
            if (arc0.containsOnCircle(ccResult.point[i]) &&
                arc1.containsOnCircle(ccResult.point[i])) {
                result.point[numIntersections] = ccResult.point[i];
                ++numIntersections;
                result.intersect = true;
            }
        }

        if (numIntersections === 2) {
            result.configuration =
                IntrArc2Arc2Configuration.NONCOCIRCULAR_TWO_POINTS;
        } else if (numIntersections === 1) {
            result.configuration =
                IntrArc2Arc2Configuration.NONCOCIRCULAR_ONE_POINT;
        } else {
            result.configuration = IntrArc2Arc2Configuration.NO_INTERSECTION;
        }

        return result;
    }
}
