// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrTriangle3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the triangles to be solids.
//
// The test-intersection query (TIQuery) uses the method of separating axes to
// determine whether or not the triangles intersect. See
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// Section 5 describes the finite set of potential separating axes.
//
// The find-intersection query (FIQuery) determines how the two triangles are
// positioned and oriented to each other. The algorithm uses the sign of the
// projections of the vertices of triangle1 onto a normal line that is
// perpendicular to the plane of triangle0. The table of possibilities is
// listed next with n = numNegative, p = numPositive and z = numZero.
//
//   n p z  intersection
//   ------------------------------------
//   0 3 0  none
//   0 2 1  vertex
//   0 1 2  edge
//   0 0 3  coplanar triangles or a triangle is degenerate
//   1 2 0  segment (2 edges clipped)
//   1 1 1  segment (1 edge clipped)
//   1 0 2  edge
//   2 1 0  segment (2 edges clipped)
//   2 0 1  vertex
//   3 0 0  none
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// TIQuery and FIQuery specializations become IntrTriangle3Triangle3TI and
// IntrTriangle3Triangle3FI. Each has two operator() overloads (stationary and
// moving triangles); per PORTING.md only the canonical two-argument query
// keeps 'test'/'find' and the moving-triangle overloads become 'testDynamic'
// and 'findDynamic'. The private static helpers become module-private
// functions. The C++ output reference parameters 'tFirst'/'tLast' become a
// mutable holder object.

import type { FIQuery } from './FIQuery.js';
import { IntrSegment2Triangle2FI } from './IntrSegment2Triangle2.js';
import { IntrTriangle2Triangle2FI } from './IntrTriangle2Triangle2.js';
import { logAssert } from './Logger.js';
import { Segment } from './Segment.js';
import type { Segment3 } from './Segment.js';
import type { TIQuery } from './TIQuery.js';
import { Triangle } from './Triangle.js';
import type { Triangle3 } from './Triangle.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { dotPerp } from './Vector2.js';
import { cross, unitCross } from './Vector3.js';

// The result of IntrTriangle3Triangle3TI.test and .testDynamic.
export interface IntrTriangle3Triangle3TIResult {
    // The contact time is 0 for stationary triangles. It is nonnegative for
    // moving triangles.
    intersect: boolean;
    contactTime: number;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrTriangle3Triangle3TIResult():
    IntrTriangle3Triangle3TIResult {
    return { intersect: false, contactTime: 0 };
}

// The result of IntrTriangle3Triangle3FI.find and .findDynamic.
export interface IntrTriangle3Triangle3FIResult {
    // The contact time is 0 for stationary triangles. It is nonnegative for
    // moving triangles.
    intersect: boolean;
    contactTime: number;

    // The intersection is a point (1 vertex), a segment (2 vertices) or a
    // convex polygon (3 or more vertices).
    intersection: Vector[];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrTriangle3Triangle3FIResult():
    IntrTriangle3Triangle3FIResult {
    return { intersect: false, contactTime: 0, intersection: [] };
}

// The C++ 'T& tFirst, T& tLast' output reference parameters.
interface TimeInterval {
    tFirst: number;
    tLast: number;
}

// Translate the triangles so that inTriangle0.v[0] becomes (0,0,0). The
// return value is the pair of translated triangles.
function translateToOrigin(inTriangle0: Triangle3, inTriangle1: Triangle3):
    { origin: Vector, triangle0: Triangle, triangle1: Triangle } {
    const origin = inTriangle0.v[0];
    const triangle0 = Triangle.fromVertices(
        Vector.zero(3),
        sub(inTriangle0.v[1], origin),
        sub(inTriangle0.v[2], origin));
    const triangle1 = Triangle.fromVertices(
        sub(inTriangle1.v[0], origin),
        sub(inTriangle1.v[1], origin),
        sub(inTriangle1.v[2], origin));
    return { origin, triangle0, triangle1 };
}

// The edge directions <V[1]-V[0], V[2]-V[1], V[0]-V[2]> of a triangle.
function edgesOf(triangle: Triangle3): [Vector, Vector, Vector] {
    return [
        sub(triangle.v[1], triangle.v[0]),
        sub(triangle.v[2], triangle.v[1]),
        sub(triangle.v[0], triangle.v[2])
    ];
}

// The triangle is <V[0],V[1],V[2]>. The line is t*direction, where the origin
// is (0,0,0) and the 'direction' is not zero but not necessarily unit length.
// The projections of the triangle vertices onto the line are
// t[i] = Dot(direction, V[i]). Return the extremes tmin = min(t[0],t[1],t[2])
// and tmax = max(t[0],t[1],t[2]).
function scaleProjectOntoLine(triangle: Triangle3, direction: Vector):
    [number, number] {
    let t = dot(direction, triangle.v[0]);
    const tExtreme: [number, number] = [t, t];
    for (let i = 1; i < 3; ++i) {
        t = dot(direction, triangle.v[i]);
        if (t < tExtreme[0]) {
            tExtreme[0] = t;
        }
        else if (t > tExtreme[1]) {
            tExtreme[1] = t;
        }
    }
    return tExtreme;
}

// The triangles are parallel exactly when Cross(N0,N1) is the zero vector.
//
// Upstream bug (IntrTriangle3Triangle3.h, the moving-triangle overloads of
// both TIQuery::operator() and FIQuery::operator()): the parallel test there
// is 'std::fabs(Dot(N0, N1)) < 1'. N0 and N1 are cross products of triangle
// edges, so they are not unit length and the comparison against 1 is
// meaningless: large non-parallel triangles take the coplanar branch (which
// omits the N1 axis and the edge-edge axes) and small parallel triangles take
// the non-parallel branch. The port fixes this by using the same criterion
// as the upstream stationary query, Dot(Cross(N0,N1),Cross(N0,N1)) > 0.
function trianglesAreParallel(N0: Vector, N1: Vector): boolean {
    const N0xN1 = cross(N0, N1);
    return dot(N0xN1, N0xN1) <= 0;
}

// This is the constant velocity separating axis test.
function testOverlapIntervals(tMax: number, speed: number,
    extreme0: readonly [number, number], extreme1: readonly [number, number],
    interval: TimeInterval): boolean {
    if (extreme1[1] < extreme0[0]) {
        // The interval extreme1 is on the left of the interval extreme0.
        if (speed <= 0) {
            // The interval extreme1 is moving away from the interval
            // extreme0.
            return false;
        }

        // Compute the first time of contact on this axis.
        let t = (extreme0[0] - extreme1[1]) / speed;
        if (t > interval.tFirst) {
            interval.tFirst = t;
        }

        if (interval.tFirst > tMax) {
            // The intersection occurs after the specified maximum time.
            return false;
        }

        // Compute the last time of contact on this axis.
        t = (extreme0[1] - extreme1[0]) / speed;
        if (t < interval.tLast) {
            interval.tLast = t;
        }

        if (interval.tFirst > interval.tLast) {
            // The time interval is invalid, so the objects are not
            // intersecting.
            return false;
        }
    }
    else if (extreme0[1] < extreme1[0]) {
        // The interval extreme1 is on the right of the interval extreme0.
        if (speed >= 0) {
            // The interval extreme1 is moving away from the interval
            // extreme0.
            return false;
        }

        // Compute the first time of contact on this axis.
        let t = (extreme0[1] - extreme1[0]) / speed;
        if (t > interval.tFirst) {
            interval.tFirst = t;
        }

        if (interval.tFirst > tMax) {
            // The intersection occurs after the specified maximum time.
            return false;
        }

        // Compute the last time of contact on this axis.
        t = (extreme0[0] - extreme1[1]) / speed;
        if (t < interval.tLast) {
            interval.tLast = t;
        }

        if (interval.tFirst > interval.tLast) {
            // The time interval is invalid, so the objects are not
            // intersecting.
            return false;
        }
    }
    else {
        // The intervals extreme0 and extreme1 are currently overlapping.
        if (speed > 0) {
            // Compute the last time of contact on this axis.
            const t = (extreme0[1] - extreme1[0]) / speed;
            if (t < interval.tLast) {
                interval.tLast = t;
            }

            if (interval.tFirst > interval.tLast) {
                // The time interval is invalid, so the objects are not
                // intersecting.
                return false;
            }
        }
        else if (speed < 0) {
            // Compute the last time of contact on this axis.
            const t = (extreme0[0] - extreme1[1]) / speed;
            if (t < interval.tLast) {
                interval.tLast = t;
            }

            if (interval.tFirst > interval.tLast) {
                // The time interval is invalid, so the objects are not
                // intersecting.
                return false;
            }
        }
    }
    return true;
}

// A projection wrapper to set up for the separating axis test.
function testOverlap(triangle0: Triangle3, triangle1: Triangle3,
    direction: Vector, tMax: number, velocity: Vector,
    interval: TimeInterval): boolean {
    const extreme0 = scaleProjectOntoLine(triangle0, direction);
    const extreme1 = scaleProjectOntoLine(triangle1, direction);
    const speed = dot(direction, velocity);
    return testOverlapIntervals(tMax, speed, extreme0, extreme1, interval);
}

// Test-intersection query for two triangles in 3D.
export class IntrTriangle3Triangle3TI implements
    TIQuery<Triangle3, Triangle3, IntrTriangle3Triangle3TIResult> {

    // The query is for stationary triangles.
    test(inTriangle0: Triangle3, inTriangle1: Triangle3):
        IntrTriangle3Triangle3TIResult {
        logAssert(inTriangle0.dimension === 3 && inTriangle1.dimension === 3,
            'IntrTriangle3Triangle3TI: mismatched sizes.');
        const result = defaultIntrTriangle3Triangle3TIResult();

        // Translate the triangles so that triangle0.v[0] becomes (0,0,0).
        const { triangle0, triangle1 } =
            translateToOrigin(inTriangle0, inTriangle1);

        // Get edge directions and a normal vector for triangle0.
        const E0 = edgesOf(triangle0);
        const N0 = cross(E0[0], E0[1]);

        // Scale-project triangle1 onto the normal line of triangle0 and test
        // for separation. The translation performed initially ensures that
        // triangle0 projects onto its normal line at t = 0.
        let tExtreme0: [number, number] = [0, 0];
        let tExtreme1 = scaleProjectOntoLine(triangle1, N0);
        if (tExtreme0[1] < tExtreme1[0] || tExtreme1[1] < tExtreme0[0]) {
            return result;
        }

        // Get edge directions and a normal vector for triangle1.
        const E1 = edgesOf(triangle1);
        const N1 = cross(E1[0], E1[1]);

        // Scale-project triangle0 onto the normal line of triangle1 and test
        // for separation.
        const projT0V0 = dot(N1, sub(triangle1.v[0], triangle0.v[0]));
        tExtreme0 = [projT0V0, projT0V0];
        tExtreme1 = scaleProjectOntoLine(triangle0, N1);
        if (tExtreme0[1] < tExtreme1[0] || tExtreme1[1] < tExtreme0[0]) {
            return result;
        }

        // At this time, neither normal line is a separation axis for the
        // triangles. If Cross(N0,N1) != (0,0,0), the planes of the triangles
        // are not parallel and must intersect in a line. If
        // Cross(N0,N1) = (0,0,0), the planes are parallel. In fact they are
        // coplanar; for if they were not coplanar, one of the two previous
        // separating axis tests would have determined this and returned from
        // the function call.

        // The potential separating axes are origin+t*direction, where origin
        // is inTriangle.v[0]. In the translated configuration, the potential
        // separating axes are t*direction.
        let direction: Vector;

        const N0xN1 = cross(N0, N1);
        const sqrLengthN0xN1 = dot(N0xN1, N0xN1);
        if (sqrLengthN0xN1 > 0) {
            // The triangles are not parallel. Test for separation by using
            // directions that are cross products of a pair of triangle edges,
            // one edge from triangle0 and one edge from triangle1.
            for (let i1 = 0; i1 < 3; ++i1) {
                for (let i0 = 0; i0 < 3; ++i0) {
                    direction = cross(E0[i0], E1[i1]);
                    tExtreme0 = scaleProjectOntoLine(triangle0, direction);
                    tExtreme1 = scaleProjectOntoLine(triangle1, direction);
                    if (tExtreme0[1] < tExtreme1[0]
                        || tExtreme1[1] < tExtreme0[0]) {
                        return result;
                    }
                }
            }
        }
        else {
            // The triangles are coplanar. Test for separation by using
            // directions that are cross products of a pair of vectors, one
            // vector a normal of a triangle and the other vector an edge from
            // the other triangle.

            // Directions N0xE0[i0].
            for (let i0 = 0; i0 < 3; ++i0) {
                direction = cross(N0, E0[i0]);
                tExtreme0 = scaleProjectOntoLine(triangle0, direction);
                tExtreme1 = scaleProjectOntoLine(triangle1, direction);
                if (tExtreme0[1] < tExtreme1[0]
                    || tExtreme1[1] < tExtreme0[0]) {
                    return result;
                }
            }

            // Directions N1xE1[i1].
            for (let i1 = 0; i1 < 3; ++i1) {
                direction = cross(N1, E1[i1]);
                tExtreme0 = scaleProjectOntoLine(triangle0, direction);
                tExtreme1 = scaleProjectOntoLine(triangle1, direction);
                if (tExtreme0[1] < tExtreme1[0]
                    || tExtreme1[1] < tExtreme0[0]) {
                    return result;
                }
            }
        }

        result.intersect = true;
        return result;
    }

    // The query is for triangles moving with constant linear velocity during
    // the time interval [0,tMax]. This is the port of the four-argument
    // upstream operator() overload.
    testDynamic(tMax: number, inTriangle0: Triangle3, velocity0: Vector,
        inTriangle1: Triangle3, velocity1: Vector):
        IntrTriangle3Triangle3TIResult {
        logAssert(inTriangle0.dimension === 3 && inTriangle1.dimension === 3,
            'IntrTriangle3Triangle3TI: mismatched sizes.');
        const result = defaultIntrTriangle3Triangle3TIResult();

        // The query determines the interval [tFirst,tLast] over which the
        // triangles are intersecting. Start with time interval [0,+infinity).
        const interval: TimeInterval = { tFirst: 0, tLast: Number.MAX_VALUE };

        // Compute the velocity of inTriangle1 relative to inTriangle0.
        const relVelocity = sub(velocity1, velocity0);

        // Translate the triangles so that triangle0.v[0] becomes (0,0,0).
        const { triangle0, triangle1 } =
            translateToOrigin(inTriangle0, inTriangle1);

        // Get edge directions and a normal vector for triangle0.
        const E0 = edgesOf(triangle0);
        const N0 = cross(E0[0], E0[1]);

        // Test for overlap using the separating axis test in the N0
        // direction.
        if (!testOverlap(triangle0, triangle1, N0, tMax, relVelocity,
            interval)) {
            return result;
        }

        // Get edge directions and a normal vector for triangle1.
        const E1 = edgesOf(triangle1);
        const N1 = cross(E1[0], E1[1]);

        if (!trianglesAreParallel(N0, N1)) {
            // The triangles are not parallel.

            // Test for overlap using the separating axis test in the N1
            // direction.
            if (!testOverlap(triangle0, triangle1, N1, tMax, relVelocity,
                interval)) {
                return result;
            }

            // Test for overlap using the separating axis test in the
            // directions E0[i0]xE1[i1].
            for (let i1 = 0; i1 < 3; ++i1) {
                for (let i0 = 0; i0 < 3; ++i0) {
                    const direction = unitCross(E0[i0], E1[i1]);
                    if (!testOverlap(triangle0, triangle1, direction, tMax,
                        relVelocity, interval)) {
                        return result;
                    }
                }
            }
        }
        else {
            // The triangles are coplanar.

            // Test for overlap using the separating axis test in the
            // directions N0xE0[i0].
            for (let i0 = 0; i0 < 3; ++i0) {
                const direction = unitCross(N0, E0[i0]);
                if (!testOverlap(triangle0, triangle1, direction, tMax,
                    relVelocity, interval)) {
                    return result;
                }
            }

            // Test for overlap using the separating axis test in the
            // directions N1xE1[i1].
            for (let i1 = 0; i1 < 3; ++i1) {
                const direction = unitCross(N1, E1[i1]);
                if (!testOverlap(triangle0, triangle1, direction, tMax,
                    relVelocity, interval)) {
                    return result;
                }
            }
        }

        result.intersect = true;
        result.contactTime = interval.tFirst;
        return result;
    }
}

// Choose the coordinate plane most aligned with the plane normal. The return
// value is the lookup table that maps the two in-plane coordinate indices and
// the discarded coordinate index.
function projectionLookup(normal: Vector): [number, number, number] {
    let maxIndex = 0;
    let cmax = Math.abs(normal.values[0]);
    let cvalue = Math.abs(normal.values[1]);
    if (cvalue > cmax) {
        maxIndex = 1;
        cmax = cvalue;
    }
    cvalue = Math.abs(normal.values[2]);
    if (cvalue > cmax) {
        maxIndex = 2;
    }

    if (maxIndex === 0) {
        // Project onto the yz-plane.
        return [1, 2, 0];
    }
    if (maxIndex === 1) {
        // Project onto the xz-plane.
        return [0, 2, 1];
    }
    // maxIndex = 2: project onto the xy-plane.
    return [0, 1, 2];
}

// Lift a 2D point of the projection plane back to the 3D plane through the
// origin with the specified normal.
function liftToPlane(normal: Vector, lookup: readonly number[],
    x0: number, x1: number): Vector {
    const point = new Vector(3);
    point.values[lookup[0]] = x0;
    point.values[lookup[1]] = x1;
    point.values[lookup[2]] = -(normal.values[lookup[0]] * x0
        + normal.values[lookup[1]] * x1) / normal.values[lookup[2]];
    return point;
}

// Compute the point, segment or polygon of intersection of coplanar
// triangles. The intersection is computed by projecting the triangles onto
// the plane and using a find-intersection query for two triangles in 2D. The
// intersection can be empty.
function getCoplanarIntersection(normal: Vector, triangle0: Triangle3,
    triangle1: Triangle3, result: IntrTriangle3Triangle3FIResult): void {
    // Project the triangles onto the coordinate plane most aligned with the
    // plane normal.
    const lookup = projectionLookup(normal);

    const projTriangle0 = new Triangle(2);
    const projTriangle1 = new Triangle(2);
    for (let i = 0; i < 3; ++i) {
        projTriangle0.v[i] = Vector.fromArray([
            triangle0.v[i].values[lookup[0]],
            triangle0.v[i].values[lookup[1]]
        ]);
        projTriangle1.v[i] = Vector.fromArray([
            triangle1.v[i].values[lookup[0]],
            triangle1.v[i].values[lookup[1]]
        ]);
    }

    // 2D triangle intersection queries require counterclockwise ordering of
    // vertices.
    const maxIndex = lookup[2];
    if (normal.values[maxIndex] < 0) {
        // Triangle0 is clockwise; reorder it.
        const swap = projTriangle0.v[1];
        projTriangle0.v[1] = projTriangle0.v[2];
        projTriangle0.v[2] = swap;
    }

    const edge0 = sub(projTriangle1.v[1], projTriangle1.v[0]);
    const edge1 = sub(projTriangle1.v[2], projTriangle1.v[0]);
    if (dotPerp(edge0, edge1) < 0) {
        // Triangle1 is clockwise; reorder it.
        const swap = projTriangle1.v[1];
        projTriangle1.v[1] = projTriangle1.v[2];
        projTriangle1.v[2] = swap;
    }

    const ttQuery = new IntrTriangle2Triangle2FI();
    const ttResult = ttQuery.find(projTriangle0, projTriangle1);
    const numVertices = ttResult.intersection.length;
    if (numVertices === 0) {
        result.intersect = false;
        result.intersection = [];
        return;
    }

    // Lift the 2D polygon of intersection to the 3D triangle space.
    result.intersect = true;
    result.intersection = [];
    for (let i = 0; i < numVertices; ++i) {
        const src = ttResult.intersection[i];
        result.intersection.push(
            liftToPlane(normal, lookup, src.values[0], src.values[1]));
    }
}

// Compute the point or segment of intersection of the 'triangle' with
// 'normal' vector. The input segment is an edge of the other triangle. The
// intersection can be empty.
function intersectsSegment(normal: Vector, triangle: Triangle3,
    segment: Segment3, result: IntrTriangle3Triangle3FIResult): void {
    // Project the triangle and segment onto the coordinate plane most aligned
    // with the plane normal.
    const lookup = projectionLookup(normal);

    const projTriangle = new Triangle(2);
    for (let i = 0; i < 3; ++i) {
        projTriangle.v[i] = Vector.fromArray([
            triangle.v[i].values[lookup[0]],
            triangle.v[i].values[lookup[1]]
        ]);
    }

    const projSegment = new Segment(2);
    for (let i = 0; i < 2; ++i) {
        projSegment.p[i] = Vector.fromArray([
            segment.p[i].values[lookup[0]],
            segment.p[i].values[lookup[1]]
        ]);
    }

    // Compute the intersection with the coincident edge and the triangle.
    const stQuery = new IntrSegment2Triangle2FI();
    const stResult = stQuery.find(projSegment, projTriangle);
    if (stResult.intersect) {
        result.intersect = true;

        // Lift the 2D intersection points to the 3D triangle space.
        result.intersection = [];
        for (let i = 0; i < stResult.numIntersections; ++i) {
            const src = stResult.point[i];
            result.intersection.push(
                liftToPlane(normal, lookup, src.values[0], src.values[1]));
        }
    }
}

// Determine whether the point is inside or strictly outside the triangle.
function containsPoint(normal: Vector, triangle: Triangle3, point: Vector,
    result: IntrTriangle3Triangle3FIResult): void {
    // Project the triangle and point onto the coordinate plane most aligned
    // with the plane normal.
    const lookup = projectionLookup(normal);

    const projTriangle = new Triangle(2);
    for (let i = 0; i < 3; ++i) {
        projTriangle.v[i] = Vector.fromArray([
            triangle.v[i].values[lookup[0]],
            triangle.v[i].values[lookup[1]]
        ]);
    }

    const projPoint = Vector.fromArray([
        point.values[lookup[0]], point.values[lookup[1]]
    ]);

    // Determine whether the point is inside or strictly outside the triangle.
    // The triangle is counterclockwise ordered when sign is +1 or clockwise
    // ordered when sign is -1.
    const maxIndex = lookup[2];
    const sign = (normal.values[maxIndex] > 0 ? 1 : -1);
    for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
        const diffPV0 = sub(projPoint, projTriangle.v[i0]);
        const diffV1V0 = sub(projTriangle.v[i1], projTriangle.v[i0]);
        if (sign * dotPerp(diffPV0, diffV1V0) > 0) {
            // The point is strictly outside edge <V[i0],V[i1]>.
            result.intersect = false;
            result.intersection = [];
            return;
        }
    }

    // Lift the 2D point of intersection to the 3D triangle space.
    result.intersect = true;
    result.intersection = [
        liftToPlane(normal, lookup, projPoint.values[0], projPoint.values[1])
    ];
}

// Support for the query for moving triangles.
enum ProjectionMap {
    // Initial value for construction of Configuration.
    PM_INVALID,

    // 3 vertices project to the same point (min = max).
    PM_3,

    // 2 vertices project to a point (min) and 1 vertex projects to a point
    // (max).
    PM_21,

    // 1 vertex projects to a point (min) and 2 vertices project to a point
    // (max).
    PM_12,

    // 1 vertex projects to a point (min), 1 vertex projects to a point (max)
    // and 1 vertex projects to a point strictly between the min and max
    // points.
    PM_111
}

// The upstream 'Configuration' struct. The 'map' and 'index' members are
// computed but never read by the upstream query (see the note on
// IntrTriangle3Triangle3FI.findDynamic); they are ported for fidelity.
interface Configuration {
    // This is how the vertices map to the projection interval.
    map: ProjectionMap;

    // The sorted indices of the vertices.
    index: number[];

    // The projection interval [min,max].
    min: number;
    max: number;
}

function defaultConfiguration(): Configuration {
    return {
        map: ProjectionMap.PM_INVALID,
        index: [0, 0, 0, 0, 0, 0, 0, 0],
        min: 0,
        max: 0
    };
}

function copyConfiguration(src: Configuration, trg: Configuration): void {
    trg.map = src.map;
    trg.index = src.index.slice();
    trg.min = src.min;
    trg.max = src.max;
}

enum ContactSide {
    CS_LEFT,
    CS_RIGHT,
    CS_NONE
}

// The upstream 'ContactSide& side' output reference parameter along with the
// contact-time configurations.
interface ContactInfo {
    side: ContactSide;
    tcfg0: Configuration;
    tcfg1: Configuration;
}

// The triangle is <V[0],V[1],V[2]>. The line is t*direction, where the origin
// is (0,0,0) and the 'direction' is not zero but not necessarily unit length.
// The projections of the triangle vertices onto the line are
// t[i] = Dot(direction, V[i]). Compute the configuration of the triangle that
// leads to the extreme interval.
function scaleProjectOntoLineConfig(triangle: Triangle3, direction: Vector,
    cfg: Configuration): void {
    // Find the projections of the triangle vertices onto the potential
    // separating axis.
    const d0 = dot(direction, triangle.v[0]);
    const d1 = dot(direction, triangle.v[1]);
    const d2 = dot(direction, triangle.v[2]);

    // Explicit sort of vertices to construct a Configuration object.
    if (d0 <= d1) {
        if (d1 <= d2) {  // d0 <= d1 <= d2
            if (d0 !== d1) {
                cfg.map = (d1 !== d2 ? ProjectionMap.PM_111
                    : ProjectionMap.PM_12);
            }
            else {  // d0 = d1
                cfg.map = (d1 !== d2 ? ProjectionMap.PM_21
                    : ProjectionMap.PM_3);
            }
            cfg.index[0] = 0;
            cfg.index[1] = 1;
            cfg.index[2] = 2;
            cfg.min = d0;
            cfg.max = d2;
        }
        else if (d0 <= d2) {  // d0 <= d2 < d1
            if (d0 !== d2) {
                cfg.map = ProjectionMap.PM_111;
                cfg.index[0] = 0;
                cfg.index[1] = 2;
                cfg.index[2] = 1;
            }
            else {
                cfg.map = ProjectionMap.PM_21;
                cfg.index[0] = 2;
                cfg.index[1] = 0;
                cfg.index[2] = 1;
            }
            cfg.min = d0;
            cfg.max = d1;
        }
        else {  // d2 < d0 <= d1
            cfg.map = (d0 !== d1 ? ProjectionMap.PM_111 : ProjectionMap.PM_12);
            cfg.index[0] = 2;
            cfg.index[1] = 0;
            cfg.index[2] = 1;
            cfg.min = d2;
            cfg.max = d1;
        }
    }
    else if (d2 <= d1) {  // d2 <= d1 < d0
        if (d2 !== d1) {
            cfg.map = ProjectionMap.PM_111;
            cfg.index[0] = 2;
            cfg.index[1] = 1;
            cfg.index[2] = 0;
        }
        else {
            cfg.map = ProjectionMap.PM_21;
            cfg.index[0] = 1;
            cfg.index[1] = 2;
            cfg.index[2] = 0;
        }
        cfg.min = d2;
        cfg.max = d0;
    }
    else if (d2 <= d0) {  // d1 < d2 <= d0
        cfg.map = (d2 !== d0 ? ProjectionMap.PM_111 : ProjectionMap.PM_12);
        cfg.index[0] = 1;
        cfg.index[1] = 2;
        cfg.index[2] = 0;
        cfg.min = d1;
        cfg.max = d0;
    }
    else {  // d1 < d0 < d2
        cfg.map = ProjectionMap.PM_111;
        cfg.index[0] = 1;
        cfg.index[1] = 0;
        cfg.index[2] = 2;
        cfg.min = d1;
        cfg.max = d2;
    }
}

// This is the constant velocity separating axis test. The cfg0 and cfg1
// inputs are the configurations for the triangles at time 0. The contact.tcfg0
// and contact.tcfg1 are the configurations at contact time.
function findOverlapIntervals(tMax: number, speed: number,
    cfg0: Configuration, cfg1: Configuration, contact: ContactInfo,
    interval: TimeInterval): boolean {
    if (cfg1.max < cfg0.min) {
        // The cfg1 interval is on the left of the cfg0 interval.
        if (speed <= 0) {
            // The cfg1 interval is moving away from the cfg0 interval.
            return false;
        }

        // Compute the first time of contact on this axis.
        let t = (cfg0.min - cfg1.max) / speed;

        if (t > interval.tFirst) {
            // Time t is the new maximum first time of contact.
            interval.tFirst = t;
            contact.side = ContactSide.CS_LEFT;
            copyConfiguration(cfg0, contact.tcfg0);
            copyConfiguration(cfg1, contact.tcfg1);
        }

        if (interval.tFirst > tMax) {
            // The intersection occurs after the specified maximum time.
            return false;
        }

        // Compute the last time of contact on this axis.
        t = (cfg0.max - cfg1.min) / speed;
        if (t < interval.tLast) {
            interval.tLast = t;
        }

        if (interval.tFirst > interval.tLast) {
            // The time interval is invalid, so the objects are not
            // intersecting.
            return false;
        }
    }
    else if (cfg0.max < cfg1.min) {
        // The cfg1 interval is on the right of the cfg0 interval.
        if (speed >= 0) {
            // The cfg1 interval is moving away from the cfg0 interval.
            return false;
        }

        // Compute the first time of contact on this axis.
        let t = (cfg0.max - cfg1.min) / speed;

        if (t > interval.tFirst) {
            // Time t is the new maximum first time of contact.
            interval.tFirst = t;
            contact.side = ContactSide.CS_RIGHT;
            copyConfiguration(cfg0, contact.tcfg0);
            copyConfiguration(cfg1, contact.tcfg1);
        }

        if (interval.tFirst > tMax) {
            // The intersection occurs after the specified maximum time.
            return false;
        }

        // Compute the last time of contact on this axis.
        t = (cfg0.min - cfg1.max) / speed;
        if (t < interval.tLast) {
            interval.tLast = t;
        }

        if (interval.tFirst > interval.tLast) {
            // The time interval is invalid, so the objects are not
            // intersecting.
            return false;
        }
    }
    else {
        // The intervals for cfg0 and cfg1 are currently overlapping.
        if (speed > 0) {
            // Compute the last time of contact on this axis.
            const t = (cfg0.max - cfg1.min) / speed;
            if (t < interval.tLast) {
                interval.tLast = t;
            }

            if (interval.tFirst > interval.tLast) {
                // The time interval is invalid, so the objects are not
                // intersecting.
                return false;
            }
        }
        else if (speed < 0) {
            // Compute the last time of contact on this axis.
            const t = (cfg0.min - cfg1.max) / speed;
            if (t < interval.tLast) {
                interval.tLast = t;
            }

            if (interval.tFirst > interval.tLast) {
                // The time interval is invalid, so the objects are not
                // intersecting.
                return false;
            }
        }
    }
    return true;
}

// A projection wrapper to set up for the separating axis test.
function findOverlap(triangle0: Triangle3, triangle1: Triangle3,
    direction: Vector, tMax: number, velocity: Vector, contact: ContactInfo,
    interval: TimeInterval): boolean {
    const cfg0 = defaultConfiguration();
    const cfg1 = defaultConfiguration();
    scaleProjectOntoLineConfig(triangle0, direction, cfg0);
    scaleProjectOntoLineConfig(triangle1, direction, cfg1);
    const speed = dot(direction, velocity);
    return findOverlapIntervals(tMax, speed, cfg0, cfg1, contact, interval);
}

// Find-intersection query for two triangles in 3D.
export class IntrTriangle3Triangle3FI implements
    FIQuery<Triangle3, Triangle3, IntrTriangle3Triangle3FIResult> {

    // The query is for stationary triangles.
    find(inTriangle0: Triangle3, inTriangle1: Triangle3):
        IntrTriangle3Triangle3FIResult {
        logAssert(inTriangle0.dimension === 3 && inTriangle1.dimension === 3,
            'IntrTriangle3Triangle3FI: mismatched sizes.');
        const result = defaultIntrTriangle3Triangle3FIResult();

        // Translate the triangles so that triangle0.v[0] becomes (0,0,0).
        const { origin, triangle0, triangle1 } =
            translateToOrigin(inTriangle0, inTriangle1);

        // Compute a normal vector for the plane containing triangle0.
        const normal = cross(triangle0.v[1], triangle0.v[2]);

        // Determine where the vertices of triangle1 live relative to the
        // plane of triangle0. The 'distance' values are actually signed and
        // scaled distances, the latter because 'normal' is not necessarily
        // unit length.
        let numPositive = 0;
        let numNegative = 0;
        let numZero = 0;
        const distance: number[] = [0, 0, 0];
        const sign: number[] = [0, 0, 0];
        for (let i = 0; i < 3; ++i) {
            distance[i] = dot(normal, triangle1.v[i]);
            if (distance[i] > 0) {
                sign[i] = 1;
                ++numPositive;
            }
            else if (distance[i] < 0) {
                sign[i] = -1;
                ++numNegative;
            }
            else {
                sign[i] = 0;
                ++numZero;
            }
        }

        if (numZero === 0) {
            if (numPositive > 0 && numNegative > 0) {
                // (n,p,z) is (1,2,0) or (2,1,0).
                const signCompare = (numPositive === 1 ? 1 : -1);
                for (let i0 = 1, i1 = 2, i2 = 0; i2 < 3; i0 = i1, i1 = i2++) {
                    if (sign[i2] === signCompare) {
                        const Vi2 = triangle1.v[i2];
                        const t0 = distance[i2] / (distance[i2] - distance[i0]);
                        const diffVi0Vi2 = sub(triangle1.v[i0], Vi2);
                        const p0 = add(Vi2, mul(t0, diffVi0Vi2));
                        const diffVi1Vi2 = sub(triangle1.v[i1], Vi2);
                        const t1 = distance[i2] / (distance[i2] - distance[i1]);
                        const p1 = add(Vi2, mul(t1, diffVi1Vi2));
                        const segment = Segment.fromEndpoints(p0, p1);
                        intersectsSegment(normal, triangle0, segment, result);
                        break;
                    }
                }
            }
            // else: (n,p,z) is (0,3,0) or (3,0,0) and triangle1 is strictly
            // on one side of the plane of triangle0, so no intersection.
        }
        else if (numZero === 1) {
            if (numPositive === 1) {
                // (n,p,z) is (1,1,1). A single vertex of triangle1 is in the
                // plane of triangle0 and the opposing edge of triangle1
                // intersects the plane transversely.
                for (let i0 = 1, i1 = 2, i2 = 0; i2 < 3; i0 = i1, i1 = i2++) {
                    if (sign[i2] === 0) {
                        const p0 = triangle1.v[i2];
                        const Vi1 = triangle1.v[i1];
                        const t = distance[i1] / (distance[i1] - distance[i0]);
                        const diffVi0Vi1 = sub(triangle1.v[i0], Vi1);
                        const p1 = add(Vi1, mul(t, diffVi0Vi1));
                        const segment = Segment.fromEndpoints(p0, p1);
                        intersectsSegment(normal, triangle0, segment, result);
                        break;
                    }
                }
            }
            else {
                // (n,p,z) is (2,0,1) or (0,2,1). A single vertex of triangle1
                // is in the plane of triangle0.
                for (let i = 0; i < 3; ++i) {
                    if (sign[i] === 0) {
                        containsPoint(normal, triangle0, triangle1.v[i],
                            result);
                        break;
                    }
                }
            }
        }
        else if (numZero === 2) {
            // (n,p,z) is (0,1,2) or (1,0,2). Two vertices are on the plane of
            // triangle0, so the segment connecting the vertices is on the
            // plane.
            for (let i0 = 1, i1 = 2, i2 = 0; i2 < 3; i0 = i1, i1 = i2++) {
                if (sign[i2] !== 0) {
                    const segment = Segment.fromEndpoints(triangle1.v[i0],
                        triangle1.v[i1]);
                    intersectsSegment(normal, triangle0, segment, result);
                    break;
                }
            }
        }
        else {  // numZero == 3
            // (n,p,z) is (0,0,3). Triangle1 is contained in the plane of
            // triangle0.
            getCoplanarIntersection(normal, triangle0, triangle1, result);
        }

        if (result.intersect) {
            // Translate the intersection set back to the original coordinate
            // system.
            for (let i = 0; i < result.intersection.length; ++i) {
                result.intersection[i] = add(result.intersection[i], origin);
            }
        }
        return result;
    }

    // The query is for triangles moving with constant linear velocity during
    // the time interval [0,tMax]. This is the port of the four-argument
    // upstream operator() overload.
    //
    // Note that upstream computes a ContactSide and the contact-time
    // Configuration objects but never reads them; the port keeps them for
    // fidelity but they do not affect the result.
    findDynamic(tMax: number, inTriangle0: Triangle3, velocity0: Vector,
        inTriangle1: Triangle3, velocity1: Vector):
        IntrTriangle3Triangle3FIResult {
        logAssert(inTriangle0.dimension === 3 && inTriangle1.dimension === 3,
            'IntrTriangle3Triangle3FI: mismatched sizes.');
        const result = defaultIntrTriangle3Triangle3FIResult();

        // The query determines the interval [tFirst,tLast] over which the
        // triangles are intersecting. Start with time interval [0,+infinity).
        const interval: TimeInterval = { tFirst: 0, tLast: Number.MAX_VALUE };

        // Compute the velocity of inTriangle1 relative to inTriangle0.
        const relVelocity = sub(velocity1, velocity0);

        // Translate the triangles so that triangle0.v[0] becomes (0,0,0).
        const { triangle0, triangle1 } =
            translateToOrigin(inTriangle0, inTriangle1);

        // Get edge directions and a normal vector for triangle0.
        const E0 = edgesOf(triangle0);
        const N0 = cross(E0[0], E0[1]);

        // Find overlap using the separating axis test in the N0 direction.
        const contact: ContactInfo = {
            side: ContactSide.CS_NONE,
            tcfg0: defaultConfiguration(),
            tcfg1: defaultConfiguration()
        };
        if (!findOverlap(triangle0, triangle1, N0, tMax, relVelocity, contact,
            interval)) {
            return result;
        }

        // Get edge directions and a normal vector for triangle1.
        const E1 = edgesOf(triangle1);
        const N1 = cross(E1[0], E1[1]);

        if (!trianglesAreParallel(N0, N1)) {
            // The triangles are not parallel.

            // Test for overlap using the separating axis test in the N1
            // direction.
            if (!findOverlap(triangle0, triangle1, N1, tMax, relVelocity,
                contact, interval)) {
                return result;
            }

            // Test for overlap using the separating axis test in the
            // directions E0[i0]xE1[i1].
            for (let i1 = 0; i1 < 3; ++i1) {
                for (let i0 = 0; i0 < 3; ++i0) {
                    const direction = unitCross(E0[i0], E1[i1]);
                    if (!findOverlap(triangle0, triangle1, direction, tMax,
                        relVelocity, contact, interval)) {
                        return result;
                    }
                }
            }
        }
        else {
            // The triangles are coplanar.

            // Test for overlap using the separating axis test in the
            // directions N0xE0[i0].
            for (let i0 = 0; i0 < 3; ++i0) {
                const direction = unitCross(N0, E0[i0]);
                if (!findOverlap(triangle0, triangle1, direction, tMax,
                    relVelocity, contact, interval)) {
                    return result;
                }
            }

            // Test for overlap using the separating axis test in the
            // directions N1xE1[i1].
            for (let i1 = 0; i1 < 3; ++i1) {
                const direction = unitCross(N1, E1[i1]);
                if (!findOverlap(triangle0, triangle1, direction, tMax,
                    relVelocity, contact, interval)) {
                    return result;
                }
            }
        }

        // Move the triangles to the first contact before computing the
        // contact set.
        const moved0 = new Triangle(3);
        const moved1 = new Triangle(3);
        for (let i = 0; i < 3; ++i) {
            moved0.v[i] = add(inTriangle0.v[i], mul(interval.tFirst, velocity0));
            moved1.v[i] = add(inTriangle1.v[i], mul(interval.tFirst, velocity1));
        }

        const stationary = this.find(moved0, moved1);
        result.intersect = true;
        result.contactTime = interval.tFirst;
        result.intersection = stationary.intersection;
        return result;
    }
}
