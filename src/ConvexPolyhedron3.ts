// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ConvexPolyhedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The convex polyhedron represented by this class has triangle faces that are
// counterclockwise ordered when viewed from outside the polyhedron. No attempt
// is made to verify that the polyhedron is convex; the caller is responsible
// for enforcing this.
//
// To support geometric algorithms that are formulated using convex quadratic
// programming (such as computing the distance from a point to a convex
// polyhedron), it is necessary to know the planes of the faces and an
// axis-aligned bounding box. If you want either the faces or the box, pass
// 'true' to the appropriate parameters. When planes are generated, the
// normals are not created to be unit length in order to support queries using
// exact rational arithmetic. If a normal to a face is N = (n0,n1,n2) and V is
// a vertex of the face, the plane is Dot(N,X-V) = 0 and is stored as
// (n0,n1,n2,-Dot(N,V)). The normals are computed to be outer pointing.
//
// Port notes:
// - The upstream constructor moves (does not copy) the input arrays. The port
//   stores the array references, the natural analogue of a move, and the
//   caller should not reuse the arrays afterward. As upstream, construction
//   "succeeds" only when there are at least 4 vertices and at least 12
//   indices; otherwise no move occurs and the member arrays have no elements.
// - The default constructor is 'new ConvexPolyhedron3()'; the four-argument
//   constructor becomes optional parameters on the same constructor.
// - Vector3 becomes a size-3 Vector and Vector4 a size-4 Vector, per the
//   runtime-dimension Vector precedent.

import { AlignedBox } from './AlignedBox.js';
import { Vector, dot, sub, hlift, computeExtremes } from './Vector.js';
import { cross } from './Vector3.js';

export class ConvexPolyhedron3 {
    vertices: Vector[];
    indices: number[];
    planes: Vector[];
    alignedBox: AlignedBox;

    constructor(inVertices?: Vector[], inIndices?: number[],
        wantPlanes: boolean = false, wantAlignedBox: boolean = false) {
        this.vertices = [];
        this.indices = [];
        this.planes = [];
        this.alignedBox = new AlignedBox(3);

        if (inVertices !== undefined && inIndices !== undefined
            && inVertices.length >= 4 && inIndices.length >= 12) {
            this.vertices = inVertices;
            this.indices = inIndices;

            if (wantPlanes) {
                this.generatePlanes();
            }

            if (wantAlignedBox) {
                this.generateAlignedBox();
            }
        }
    }

    // If you modify the vertices or indices and you want the new face planes
    // or aligned box computed, call these functions.
    generatePlanes(): void {
        if (this.vertices.length > 0 && this.indices.length > 0) {
            const numTriangles = Math.floor(this.indices.length / 3);
            this.planes = new Array<Vector>(numTriangles);
            for (let t = 0, i = 0; t < numTriangles; ++t) {
                const V0 = this.vertices[this.indices[i++]];
                const V1 = this.vertices[this.indices[i++]];
                const V2 = this.vertices[this.indices[i++]];
                const E1 = sub(V1, V0);
                const E2 = sub(V2, V0);
                const N = cross(E1, E2);
                this.planes[t] = hlift(N, -dot(N, V0));
            }
        }
    }

    generateAlignedBox(): void {
        if (this.vertices.length > 0 && this.indices.length > 0) {
            const extremes = computeExtremes(this.vertices);
            if (extremes !== null) {
                this.alignedBox.min = extremes.vmin;
                this.alignedBox.max = extremes.vmax;
            }
        }
    }
}
