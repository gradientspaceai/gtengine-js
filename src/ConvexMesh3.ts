// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ConvexMesh3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream is the class template ConvexMesh3<Real>, a plain
// container of public data members. The port keeps it as such. The member
// type aliases 'Vertex' and 'Triangle' become the exported type aliases
// ConvexMesh3Vertex and ConvexMesh3Triangle, prefixed with the owning class
// name because the library re-exports every symbol from a single flat index.
// A Vertex is the runtime-sized Vector with size 3 (upstream Vector3<Real>).
// Upstream requires Real to be an arbitrary-precision type supporting
// division (contract item 3 below); the port instantiates the class for the
// 'number' Vector as usual, and the requirement is a statement about the
// caller's queries, not about this container, which performs no arithmetic.

import { Vector } from './Vector.js';

// The upstream 'Vertex' alias: a 3D point of the mesh.
export type ConvexMesh3Vertex = Vector;

// The upstream 'Triangle' alias: the three vertex indices of a face.
export type ConvexMesh3Triangle = [number, number, number];

export class ConvexMesh3 {
    // A client of ConvexMesh3 is responsible for populating the vertices
    // and indices so that the resulting mesh represents a convex
    // polyhedron.
    //   1. All elements of 'vertices' must be used by the polyhedron.
    //   2. The triangle faces must have the same chirality when viewed
    //      from outside the polyhedron. They are all counterclockwise
    //      oriented or all clockwise oriented when viewed from outside
    //      the polyhedron.
    //   3. The Real type must be an arbitrary-precision type that
    //      supports division.
    //   4. The polyhedron can be degenerate. All the possibilities are
    //      listed next.
    //        point:
    //          vertices.length === 1, triangles.length === 0
    //
    //        line segment:
    //          vertices.length === 2, triangles.length === 0
    //
    //        convex polygon:
    //          vertices.length >= 3, triangles.length > 0 and the
    //          vertices are coplanar
    //
    //        convex polyhedron:
    //          vertices.length >= 3, triangles.length > 0 and the
    //          vertices are not coplanar

    static readonly CFG_EMPTY = 0x00000000;
    static readonly CFG_POINT = 0x00000001;
    static readonly CFG_SEGMENT = 0x00000002;
    static readonly CFG_POLYGON = 0x00000004;
    static readonly CFG_POLYHEDRON = 0x00000008;

    configuration: number;
    vertices: ConvexMesh3Vertex[];
    triangles: ConvexMesh3Triangle[];

    constructor() {
        this.configuration = ConvexMesh3.CFG_EMPTY;
        this.vertices = [];
        this.triangles = [];
    }
}
