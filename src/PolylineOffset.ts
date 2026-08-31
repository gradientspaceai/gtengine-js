// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PolylineOffset.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The polyline has N vertices. If the polyline is open (with N >= 2), the
// segments are <V[0],V[1]>, <V[1],V[2]>, ..., <V[N-2],V[N-1]>. If the
// polyline is closed (with N >= 3), the segments are those of the open
// polyline and the segment <V[N-1],V[0]>. The geometry of the polyline is not
// taken into account in the algorithm. For example, the algorithm does not
// test whether the segments intersect at interior points. If you want an
// offset for a simple polygon, you must ensure that the incoming points form
// a simple polygon. The offset itself might not be a simple polygon when the
// offset distance is sufficiently large.
//
// The segment <V[i], V[i + 1]> is directed with unit-length direction
//   D = (V[i + 1] - V[i])/|V[i + 1] - V[i]|
// A unit-length normal to the segment is chosen to point to the right of the
// segment,
//   N = Perp(D)
// where Perp(x,y) = (y,-x).
//
// For 3 consecutive vertices <V[i], V[i + 1], V[i + 2]>, it is allowed that
// the directed segments <V[i], V[i + 1]> and <V[i + 1], V[i + 2]> be parallel
// as long as the direction vectors are in the same direction, that is,
// Dot(D[i], D[i + 1]) = 1. It is not allowed that the direction vectors are
// in the opposite direction, Dot(D[i], D[i + 1]) = -1, because then there is
// a singularity in the offset distance at D[i + 1].
//
// To compute the offset polyline in the positive normal direction (offset is
// to the right of segments), set 'offsetRight' to true; the offset polyline
// is the 'rightPolyline' member of the result. To compute the offset polyline
// in the negative normal direction (offset is to the left of segments), set
// 'offsetLeft' to true; the offset polyline is the 'leftPolyline' member. You
// can set both Boolean values to true when you want both polylines.
//
// NOTE: The offset depends on the geometry of the polyline. As the offset
// distance increases, the offset polylines can "fold over". The visualization
// will not look right. This code makes no attempt to determine a maximum
// offset distance for which fold-over occurs once you exceed that maximum.
//
// Port notes:
// - Upstream 'Vector2<T>' is the runtime-sized Vector of size 2 (create with
//   'Vector.fromArray([x, y])').
// - Execute writes the two offset polylines through reference parameters; the
//   port returns them in a PolylineOffsetResult object. A polyline that was
//   not requested is an empty array.
// - Upstream stores the incoming vertices as a const reference (hence the
//   deleted copy/move operations). The port likewise stores the caller's
//   array by reference; the caller must not modify it while the query object
//   is alive, because the directions and normals are precomputed from it in
//   the constructor.

import { logAssert } from './Logger';
import { Vector, add, dot, mul, normalize, sub } from './Vector';
import { perp } from './Vector2';

export interface PolylineOffsetResult {
    // The offset polyline to the right of the segments (positive normal
    // direction), or an empty array when 'offsetRight' was false.
    rightPolyline: Vector[];

    // The offset polyline to the left of the segments (negative normal
    // direction), or an empty array when 'offsetLeft' was false.
    leftPolyline: Vector[];
}

export class PolylineOffset {
    private mVertices: readonly Vector[];
    private mIsOpen: boolean;
    private mDirections: Vector[];
    private mNormals: Vector[];

    constructor(vertices: readonly Vector[], isOpen: boolean) {
        // The validation is performed before any sizing. Upstream sizes
        // mDirections/mNormals to 'vertices.size() - 1' in the constructor
        // initializer list, which underflows for an empty input before
        // LogAssert can fire; see the PR notes.
        const numVertices = vertices.length;
        logAssert(
            numVertices >= (isOpen ? 2 : 3),
            'Invalid number of polyline vertices.');

        this.mVertices = vertices;
        this.mIsOpen = isOpen;

        const numEdges = isOpen ? numVertices - 1 : numVertices;
        this.mDirections = new Array<Vector>(numEdges);
        this.mNormals = new Array<Vector>(numEdges);

        for (let i0 = 0, i1 = 1; i1 < numVertices; i0 = i1++) {
            const V0 = this.mVertices[i0];
            const V1 = this.mVertices[i1];
            const D0 = sub(V1, V0);
            normalize(D0);
            this.mDirections[i0] = D0;
            this.mNormals[i0] = perp(D0);
        }

        if (!this.mIsOpen) {
            const numVerticesM1 = numVertices - 1;
            const V0 = this.mVertices[numVerticesM1];
            const V1 = this.mVertices[0];
            const D0 = sub(V1, V0);
            normalize(D0);
            this.mDirections[numVerticesM1] = D0;
            this.mNormals[numVerticesM1] = perp(D0);
        }
    }

    // Member access for the precomputed per-segment frames. Upstream keeps
    // these private; the port exposes read-only accessors, which is useful
    // for testing and for callers that need the segment normals.
    getDirections(): readonly Vector[] {
        return this.mDirections;
    }

    getNormals(): readonly Vector[] {
        return this.mNormals;
    }

    execute(offsetDistance: number, offsetRight: boolean,
        offsetLeft: boolean): PolylineOffsetResult {
        logAssert(
            offsetDistance > 0,
            'The offset distance must be positive.');

        logAssert(
            offsetRight || offsetLeft,
            'Expecting a directive to compute an offset polyline.');

        const result: PolylineOffsetResult = {
            rightPolyline: [],
            leftPolyline: []
        };

        if (offsetRight) {
            result.rightPolyline = this.computeOffsetPolyline(offsetDistance, 1);
        }

        if (offsetLeft) {
            result.leftPolyline = this.computeOffsetPolyline(offsetDistance, -1);
        }

        return result;
    }

    // Upstream ComputeRightPolyline and ComputeLeftPolyline are identical
    // except that the left version subtracts the offset terms rather than
    // adding them. The port merges them, with 'sign' equal to +1 for the
    // right polyline and -1 for the left polyline.
    private computeOffsetPolyline(distance: number, sign: number): Vector[] {
        const numVertices = this.mVertices.length;
        const numNormals = this.mNormals.length;
        const polyline: Vector[] = [];

        // Process the first endpoint depending on whether the polyline is
        // open or closed.
        const one = 1;
        if (this.mIsOpen) {
            const V0 = this.mVertices[0];
            const N0 = this.mNormals[0];
            polyline.push(add(V0, mul(sign * distance, N0)));
        }
        else {
            const V1 = this.mVertices[0];
            const N0 = this.mNormals[numNormals - 1];
            const N1 = this.mNormals[0];
            polyline.push(add(V1, mul(sign * distance / (one + dot(N0, N1)),
                add(N0, N1))));
        }

        // B = N[i0] + N[i1] is the bisector direction at V[i1] for the two
        // normals of the edges sharing V[i1]. The offset vertex is
        // V[i1] + (d / Dot(N[i0],B)) * B, where d is the segment offset
        // distance. B does not have to be normalized because the offset
        // vertex is independent of the length of B. The offset vertex is
        // therefore
        //   V[i1] + (d / (1 + Dot(N[i0], N[i1]))) * (N[i0] + N[i1])
        for (let i0 = 0, i1 = 1, i2 = 2; i2 < numVertices; i0 = i1, i1 = i2++) {
            const V1 = this.mVertices[i1];
            const N0 = this.mNormals[i0];
            const N1 = this.mNormals[i1];
            polyline.push(add(V1, mul(sign * distance / (one + dot(N0, N1)),
                add(N0, N1))));
        }

        // Process the last endpoint depending on whether the polyline is open
        // or closed.
        if (this.mIsOpen) {
            const V0 = this.mVertices[numVertices - 1];
            const N0 = this.mNormals[numNormals - 1];
            polyline.push(add(V0, mul(sign * distance, N0)));
        }
        else {
            const V1 = this.mVertices[numVertices - 1];
            const N0 = this.mNormals[numVertices - 2];
            const N1 = this.mNormals[numNormals - 1];
            polyline.push(add(V1, mul(sign * distance / (one + dot(N0, N1)),
                add(N0, N1))));
        }

        return polyline;
    }
}
