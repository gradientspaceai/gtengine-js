// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PolygonTree.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// These classes are used by the triangulators TriangulateEC (ear clipping)
// and TriangulateCDT (constrained Delaunay triangulation). PolygonTree used
// to be a nested class of those, but it was factored out so that
// applications can use either triangulator without duplicating the trees.
//
// NOTE: The polygon member does not duplicate endpoints. For example, if P[]
// are the point locations and the polygon is a triangle with
// counterclockwise ordering <P[i0],P[i1],P[i2]>, then polygon = {i0,i1,i2}.
// The implication is that there are 3 directed edges: {P[i0],P[i1]},
// {P[i1],P[i2]} and {P[i2],P[i0]}.
//
// Port notes:
//   - The nested class PolygonTreeEx::Node becomes the top-level exported
//     class PolygonTreeExNode, because exported names must be unique across
//     the library.
//   - The three GetContainingTriangle overloads become the distinct methods
//     getContainingTriangle (the tree search), getContainingTriangleInList
//     (a list plus per-triangle node indices) and
//     getContainingTriangleWithChirality (a list with a single chirality,
//     static as upstream's is effectively static).
//   - The sentinel std::numeric_limits<size_t>::max() becomes
//     PolygonTreeEx.INVALID (Number.MAX_SAFE_INTEGER). The two-value returns
//     std::pair<size_t,size_t> become { nIndex, tIndex } objects.
//   - Triangles are 3-element index arrays and points are Vector objects of
//     size 2.

import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

// A tree of nested polygons. The root node corresponds to an outer polygon.
// The children of the root correspond to inner polygons, which are polygons
// strictly contained in the outer polygon. Each inner polygon may itself
// contain an outer polygon which in turn can contain inner polygons, thus
// leading to a hierarchy of polygons. The outer polygons have vertices
// listed in counterclockwise order. The inner polygons have vertices listed
// in clockwise order.
export class PolygonTree {
    polygon: number[];
    child: PolygonTree[];

    constructor() {
        this.polygon = [];
        this.child = [];
    }
}

// A node of PolygonTreeEx.
//
// The chirality (winding ordering of the polygon) is +1 for a
// counterclockwise-ordered polygon or -1 for a clockwise-ordered polygon.
//
// The triangulation is computed by the triangulators and explicitly stored
// per tree node.
export class PolygonTreeExNode {
    polygon: number[];
    chirality: number;
    triangulation: [number, number, number][];
    self: number;
    parent: number;
    minChild: number;
    supChild: number;

    constructor() {
        this.polygon = [];
        this.chirality = 0;
        this.triangulation = [];
        this.self = 0;
        this.parent = 0;
        this.minChild = 0;
        this.supChild = 0;
    }
}

// A tree of nested polygons with extra information about the polygon. The
// point locations are specified separately to the triangulators.
//
// The element nodes[0] is the root of the tree with nodes[0].parent = -1
// (unused). If nodes[0] has C children, then nodes[0].minChild = 1 and
// nodes[0].supChild = 1 + C. Generally, nodes[i] is a node with parent
// nodes[p], where p = nodes[i].parent, and children nodes[c], where
// nodes[i].minChild <= c < nodes[i].supChild. If nodes[i].minChild >=
// nodes[i].supChild, the node has no children.
export class PolygonTreeEx {
    // The port of std::numeric_limits<size_t>::max(), the "no such triangle"
    // sentinel.
    static readonly INVALID = Number.MAX_SAFE_INTEGER;

    // The nodes of the polygon tree, organized based on a breadth-first
    // traversal of the tree.
    nodes: PolygonTreeExNode[];

    // These members support TriangulateCDT. The *NodeIndices members store
    // the indices into 'nodes[]' for the triangles in the *Triangles
    // members. For example, the triangle interiorTriangles[t] comes from
    // nodes[interiorNodeIndices[t]].

    // The triangles in the polygon tree that cover each region bounded by an
    // outer polygon and its contained inner polygons. This set is the
    // equivalent of the output of TriangulateEC that uses ear clipping.
    interiorTriangles: [number, number, number][];
    interiorNodeIndices: number[];

    // The triangles in the polygon tree that cover each region bounded by an
    // inner polygon and its contained outer polygons.
    exteriorTriangles: [number, number, number][];
    exteriorNodeIndices: number[];

    // The triangles inside the polygon tree:
    //   insideTriangles = interiorTriangles + exteriorTriangles
    insideTriangles: [number, number, number][];
    insideNodeIndices: number[];

    // The triangles inside the convex hull of the Delaunay triangles but
    // outside the polygon tree. These triangles are not associated with any
    // 'nodes[]' element.
    outsideTriangles: [number, number, number][];

    // All the triangles:
    //   allTriangles = insideTriangles + outsideTriangles.
    allTriangles: [number, number, number][];

    constructor() {
        this.nodes = [];
        this.interiorTriangles = [];
        this.interiorNodeIndices = [];
        this.exteriorTriangles = [];
        this.exteriorNodeIndices = [];
        this.insideTriangles = [];
        this.insideNodeIndices = [];
        this.outsideTriangles = [];
        this.allTriangles = [];
    }

    // Search the polygon tree for the triangle that contains 'test'. If
    // there is such a triangle, the returned pair (nIndex,tIndex) states
    // that the triangle is nodes[nIndex].triangulation[tIndex]. If there is
    // no such triangle, the returned pair is (INVALID, INVALID). The
    // function is naturally recursive, but simulated recursion is used to
    // avoid a large program stack by instead using the heap.
    getContainingTriangle(test: Vector, points: Vector[]):
        { nIndex: number, tIndex: number } {
        const smax = PolygonTreeEx.INVALID;
        const result = { nIndex: smax, tIndex: smax };

        const stack: number[] = [0];
        while (stack.length > 0 && result.nIndex === smax) {
            const nIndex = stack.pop() as number;
            const node = this.nodes[nIndex];
            for (let c = node.minChild; c < node.supChild; ++c) {
                stack.push(c);
            }

            for (let tIndex = 0; tIndex < node.triangulation.length; ++tIndex) {
                if (PolygonTreeEx.pointInTriangle(test, node.chirality,
                    node.triangulation[tIndex], points)) {
                    result.nIndex = nIndex;
                    result.tIndex = tIndex;
                    break;
                }
            }
        }

        return result;
    }

    // Search the triangles for the triangle that contains 'test'. If there
    // is such a triangle, the returned pair (nIndex,tIndex) states that the
    // triangle is nodes[nIndex].triangulation[tIndex]. If there is no such
    // triangle, the returned pair is (INVALID, INVALID). The function uses a
    // linear search of the input triangles. Typical calls pass
    // (insideTriangles, insideNodeIndices), (interiorTriangles,
    // interiorNodeIndices) or (exteriorTriangles, exteriorNodeIndices).
    getContainingTriangleInList(test: Vector, triangles: [number, number, number][],
        nodeIndices: number[], points: Vector[]): { nIndex: number, tIndex: number } {
        logAssert(triangles.length === nodeIndices.length, 'Invalid argument.');

        const smax = PolygonTreeEx.INVALID;
        const result = { nIndex: smax, tIndex: smax };
        for (let tIndex = 0; tIndex < triangles.length; ++tIndex) {
            const nIndex = nodeIndices[tIndex];
            const node = this.nodes[nIndex];
            if (PolygonTreeEx.pointInTriangle(test, node.chirality,
                triangles[tIndex], points)) {
                result.nIndex = nIndex;
                result.tIndex = tIndex;
                break;
            }
        }
        return result;
    }

    // Search the triangles for the triangle that contains 'test'. If there
    // is such a triangle, the returned t-value is in the range
    // 0 <= t < triangles.length; otherwise, INVALID is returned. The
    // function uses a linear search of the input triangles. No information
    // is available about the 'nodes[]' element corresponding to the
    // containing triangle of the test point.
    static getContainingTriangleWithChirality(test: Vector,
        triangles: [number, number, number][], chirality: number,
        points: Vector[]): number {
        let result = PolygonTreeEx.INVALID;
        for (let tIndex = 0; tIndex < triangles.length; ++tIndex) {
            if (PolygonTreeEx.pointInTriangle(test, chirality, triangles[tIndex], points)) {
                result = tIndex;
                break;
            }
        }
        return result;
    }

    // Determine whether 'test' is inside the triangle whose vertices are
    // points[triangle[0]], points[triangle[1]] and points[triangle[2]]. If
    // the points are counterclockwise ordered, set 'chirality' to +1. If the
    // points are clockwise ordered, set 'chirality' to -1.
    private static pointInTriangle(test: Vector, chirality: number,
        triangle: [number, number, number], points: Vector[]): boolean {
        const sign = chirality;
        for (let i1 = 0, i0 = 2; i1 < 3; i0 = i1++) {
            const p0 = points[triangle[i0]].values;
            const p1 = points[triangle[i1]].values;
            const nx = p1[1] - p0[1];
            const ny = p0[0] - p1[0];
            const dx = test.values[0] - p0[0];
            const dy = test.values[1] - p0[1];
            const sdot = sign * (nx * dx + ny * dy);
            if (sdot > 0) {
                return false;
            }
        }
        return true;
    }
}
