// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ExtremalQuery3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The extremal queries for convex objects are based on the algorithm
// described in
// https://www.geometrictools.com/Documentation/ExtremalPolytopeQueries.pdf
//
// Port notes:
// - The upstream abstract base class becomes a TypeScript abstract class. The
//   protected constructor is public-by-necessity in TypeScript syntax terms
//   but marked 'protected' so only subclasses can call it, matching upstream.
// - Upstream stores 'Polyhedron3<Real> const&', a reference to a polytope
//   whose lifetime the caller owns. The port stores the object, which is a
//   reference in TypeScript; no copy is made, as upstream.
// - GetExtremeVertices has two output reference parameters; the port returns
//   an object { positiveDirection, negativeDirection } holding the indices of
//   the extreme vertices in the polyhedron vertex array.

import { Polyhedron3 } from './Polyhedron3';
import { Vector, sub } from './Vector';
import { unitCross } from './Vector3';

// The result of ExtremalQuery3.getExtremeVertices: the indices, in the
// polyhedron vertex array, of the vertices extreme in the specified
// direction and in its negation.
export interface ExtremalQuery3Result {
    positiveDirection: number;
    negativeDirection: number;
}

export abstract class ExtremalQuery3 {
    protected mPolytope: Polyhedron3;
    protected mFaceNormals: Vector[];

    // The caller must ensure that the input polyhedron is convex.
    protected constructor(polytope: Polyhedron3) {
        this.mPolytope = polytope;
        this.mFaceNormals = [];

        // Create the face normals.
        const vertexPool = this.mPolytope.getVertices();
        const indices = this.mPolytope.getIndices();
        const numTriangles = Math.floor(indices.length / 3);
        this.mFaceNormals = new Array<Vector>(numTriangles);
        for (let t = 0; t < numTriangles; ++t) {
            const v0 = vertexPool[indices[3 * t + 0]];
            const v1 = vertexPool[indices[3 * t + 1]];
            const v2 = vertexPool[indices[3 * t + 2]];
            const edge1 = sub(v1, v0);
            const edge2 = sub(v2, v0);
            this.mFaceNormals[t] = unitCross(edge1, edge2);
        }
    }

    // Member access.
    getPolytope(): Polyhedron3 {
        return this.mPolytope;
    }

    getFaceNormals(): Vector[] {
        return this.mFaceNormals;
    }

    // Compute the extreme vertices in the specified direction and return the
    // indices of the vertices in the polyhedron vertex array.
    abstract getExtremeVertices(direction: Vector): ExtremalQuery3Result;
}
