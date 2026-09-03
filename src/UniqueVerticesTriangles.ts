// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) UniqueVerticesTriangles.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// UniqueVerticesTriangles is a helper class that provides support for several
// mesh generation and mesh reduction operations. The vertices have type
// VertexType, for which duplicates must be detectable because duplicate
// vertices are eliminated in the operations.
//
//   1. Generate an indexed triangle representation from an array of triples
//      of VertexType. Each triple represents the vertices of a triangle.
//      Presumably, the triangles share vertices. The output is an array of
//      unique VertexType objects (a vertex pool) and an array of triples of
//      indices into the pool, each triple representing a triangle.
//
//   2. Remove duplicate vertices from a vertex pool used by an indexed
//      triangle representation. A new vertex pool of unique vertices is
//      generated and the indexed triangles are modified to be indices into
//      this vertex pool.
//
//   3. Remove unused vertices from a vertex pool used by an indexed triangle
//      representation. A new vertex pool of unique vertices is generated and
//      the indexed triangles are modified to be indices into the new vertex
//      pool.
//
//   4. Remove duplicate and unused vertices from a vertex pool, a combination
//      of the operations in #2 and #3.
//
// Port notes:
//   - Upstream requires VertexType to have a less-than comparison so that it
//     can be a std::map key. Only the induced equivalence affects the output:
//     the unique vertices are numbered in order of first occurrence and the
//     packing loop scatters them to those indices, so the sorted iteration
//     order of the map is not observable. The port therefore takes a key
//     function that maps a vertex to a primitive key and uses a Map, which
//     reproduces the equivalence without requiring an ordering.
//   - RemoveUnused does depend on ordering: it uses std::set<int32_t> to
//     visit the used vertex indices in increasing order. The port sorts the
//     used indices numerically to reproduce that.
//   - Upstream has pairs of member functions distinguished by whether the
//     indices are a flat array or an array of triples. The port names the
//     latter with the suffix 'Triples'.
//   - The output reference parameters become returned object literals.
//   - The upstream precondition checks are compiled out unless
//     GTL_VALIDATE_UNIQUE_VERTICES_TRIANGLES is defined. The port has the
//     'validate' property, which defaults to false to match upstream.

import { logAssert } from './Logger.js';

// The default vertex key. Arrays (for example [x, y, z] positions) are joined
// with commas; other objects are serialized with JSON.stringify; primitives
// are converted to strings. Pass a custom key function to the constructor
// when this is not the equivalence you need.
export function defaultVertexKey(vertex: unknown): string {
    if (Array.isArray(vertex)) {
        return vertex.join(',');
    }
    if (vertex !== null && typeof vertex === 'object') {
        return JSON.stringify(vertex);
    }
    return String(vertex);
}

export class UniqueVerticesTriangles<VertexType> {
    // Port of the GTL_VALIDATE_UNIQUE_VERTICES_TRIANGLES preprocessor symbol.
    // When true, the preconditions of the member functions are verified.
    validate: boolean = false;

    private keyOf: (vertex: VertexType) => string | number;

    constructor(keyOf: (vertex: VertexType) => string | number = defaultVertexKey) {
        this.keyOf = keyOf;
    }

    // See #1 in the comments at the beginning of this file. The preconditions
    // are
    //   1. inVertices.length is a positive multiple of 3
    // The postconditions are
    //   1. vertices has unique vertices
    //   2. indices.length = inVertices.length
    //   3. 0 <= indices[i] < vertices.length
    generateIndexedTriangles(inVertices: VertexType[]): {
        vertices: VertexType[];
        indices: number[];
    } {
        if (this.validate) {
            logAssert(inVertices.length > 0 && inVertices.length % 3 === 0,
                'Invalid number of vertices.');
        }

        const indices = new Array<number>(inVertices.length).fill(0);
        const vertices = this.removeDuplicates(inVertices, indices);
        return { vertices, indices };
    }

    // See #1 in the comments at the beginning of this file. The preconditions
    // are
    //   1. inVertices.length is a positive multiple of 3
    // The postconditions are
    //   1. vertices has unique vertices
    //   2. triangles.length = inVertices.length / 3
    //   3. 0 <= triangles[i][j] < vertices.length
    generateIndexedTrianglesTriples(inVertices: VertexType[]): {
        vertices: VertexType[];
        triangles: number[][];
    } {
        const result = this.generateIndexedTriangles(inVertices);
        return {
            vertices: result.vertices,
            triangles: UniqueVerticesTriangles.toTriples(result.indices)
        };
    }

    // See #2 in the comments at the beginning of the file. The preconditions
    // are
    //   1. inVertices.length is positive
    //   2. inIndices.length is a positive multiple of 3
    //   3. 0 <= inIndices[i] < inVertices.length
    // The postconditions are
    //   1. vertices has unique vertices
    //   2. indices.length = inIndices.length
    //   3. 0 <= indices[i] < vertices.length
    removeDuplicateVertices(inVertices: VertexType[], inIndices: number[]): {
        vertices: VertexType[];
        indices: number[];
    } {
        if (this.validate) {
            logAssert(inVertices.length > 0, 'Invalid number of vertices.');
            logAssert(inIndices.length > 0 && inIndices.length % 3 === 0,
                'Invalid number of indices.');
            for (const index of inIndices) {
                logAssert(0 <= index && index < inVertices.length, 'Invalid index.');
            }
        }

        const inToOutMapping = new Array<number>(inVertices.length).fill(0);
        const vertices = this.removeDuplicates(inVertices, inToOutMapping);

        const indices = new Array<number>(inIndices.length).fill(0);
        for (let i = 0; i < inIndices.length; ++i) {
            indices[i] = inToOutMapping[inIndices[i]];
        }
        return { vertices, indices };
    }

    // See #2 in the comments at the beginning of the file. The preconditions
    // are
    //   1. inVertices.length is positive
    //   2. inTriangles.length is positive
    //   3. 0 <= inTriangles[i][j] < inVertices.length
    // The postconditions are
    //   1. vertices has unique vertices
    //   2. triangles.length = inTriangles.length
    //   3. 0 <= triangles[i][j] < vertices.length
    removeDuplicateVerticesTriples(inVertices: VertexType[], inTriangles: number[][]): {
        vertices: VertexType[];
        triangles: number[][];
    } {
        if (this.validate) {
            logAssert(inTriangles.length > 0, 'Invalid number of triangles.');
        }

        const result = this.removeDuplicateVertices(inVertices,
            UniqueVerticesTriangles.toFlat(inTriangles));
        return {
            vertices: result.vertices,
            triangles: UniqueVerticesTriangles.toTriples(result.indices)
        };
    }

    // See #3 in the comments at the beginning of the file. The preconditions
    // are
    //   1. inVertices.length is positive
    //   2. inIndices.length is a positive multiple of 3
    //   3. 0 <= inIndices[i] < inVertices.length
    // The postconditions are
    //   1. vertices.length is positive
    //   2. indices.length = inIndices.length
    //   3. 0 <= indices[i] < vertices.length
    //   4. each vertices[j] occurs at least once in indices
    removeUnusedVertices(inVertices: VertexType[], inIndices: number[]): {
        vertices: VertexType[];
        indices: number[];
    } {
        if (this.validate) {
            logAssert(inVertices.length > 0, 'Invalid number of vertices.');
            logAssert(inIndices.length > 0 && inIndices.length % 3 === 0,
                'Invalid number of indices.');
            for (const index of inIndices) {
                logAssert(0 <= index && index < inVertices.length, 'Invalid index.');
            }
        }

        const indices = new Array<number>(inIndices.length).fill(0);
        const vertices = this.removeUnused(inVertices, inIndices, indices);
        return { vertices, indices };
    }

    // See #3 in the comments at the beginning of the file. The preconditions
    // are
    //   1. inVertices.length is positive
    //   2. inTriangles.length is positive
    //   3. 0 <= inTriangles[i][j] < inVertices.length
    // The postconditions are
    //   1. vertices.length is positive
    //   2. triangles.length = inTriangles.length
    //   3. 0 <= triangles[i][j] < vertices.length
    //   4. each vertices[j] occurs at least once in triangles
    removeUnusedVerticesTriples(inVertices: VertexType[], inTriangles: number[][]): {
        vertices: VertexType[];
        triangles: number[][];
    } {
        if (this.validate) {
            logAssert(inTriangles.length > 0, 'Invalid number of triangles.');
        }

        const result = this.removeUnusedVertices(inVertices,
            UniqueVerticesTriangles.toFlat(inTriangles));
        return {
            vertices: result.vertices,
            triangles: UniqueVerticesTriangles.toTriples(result.indices)
        };
    }

    // See #4 and the preconditions for removeDuplicateVertices and for
    // removeUnusedVertices.
    removeDuplicateAndUnusedVertices(inVertices: VertexType[], inIndices: number[]): {
        vertices: VertexType[];
        indices: number[];
    } {
        const temp = this.removeDuplicateVertices(inVertices, inIndices);
        return this.removeUnusedVertices(temp.vertices, temp.indices);
    }

    // See #4 and the preconditions for removeDuplicateVerticesTriples and for
    // removeUnusedVerticesTriples.
    removeDuplicateAndUnusedVerticesTriples(inVertices: VertexType[], inTriangles: number[][]): {
        vertices: VertexType[];
        triangles: number[][];
    } {
        const temp = this.removeDuplicateVerticesTriples(inVertices, inTriangles);
        return this.removeUnusedVerticesTriples(temp.vertices, temp.triangles);
    }

    private removeDuplicates(inVertices: VertexType[], inToOutMapping: number[]): VertexType[] {
        // Construct the unique vertices.
        const numInVertices = inVertices.length;
        let numOutVertices = 0;
        const vmap = new Map<string | number, { vertex: VertexType; index: number }>();
        for (let i = 0; i < numInVertices; ++i) {
            const key = this.keyOf(inVertices[i]);
            const element = vmap.get(key);
            if (element !== undefined) {
                // The vertex is a duplicate of one inserted earlier into the
                // map. Its index i will be modified to that of the
                // first-found vertex.
                inToOutMapping[i] = element.index;
            } else {
                // The vertex occurs for the first time.
                vmap.set(key, { vertex: inVertices[i], index: numOutVertices });
                inToOutMapping[i] = numOutVertices;
                ++numOutVertices;
            }
        }

        // Pack the unique vertices into an array.
        const outVertices = new Array<VertexType>(numOutVertices);
        for (const element of vmap.values()) {
            outVertices[element.index] = element.vertex;
        }
        return outVertices;
    }

    private removeUnused(inVertices: VertexType[], inIndices: number[],
        outIndices: number[]): VertexType[] {
        // The upstream std::set visits the used indices in increasing order.
        const usedIndices = [...new Set<number>(inIndices)].sort((a, b) => a - b);

        // Locate the used vertices and pack them into an array.
        const outVertices = new Array<VertexType>(usedIndices.length);
        let numOutVertices = 0;
        const vmap = new Map<number, number>();
        for (const oldIndex of usedIndices) {
            outVertices[numOutVertices] = inVertices[oldIndex];
            vmap.set(oldIndex, numOutVertices);
            ++numOutVertices;
        }

        // Reassign the old indices to the new indices.
        for (let i = 0; i < inIndices.length; ++i) {
            outIndices[i] = vmap.get(inIndices[i]) as number;
        }
        return outVertices;
    }

    private static toFlat(triangles: number[][]): number[] {
        const indices = new Array<number>(3 * triangles.length);
        for (let t = 0; t < triangles.length; ++t) {
            for (let j = 0; j < 3; ++j) {
                indices[3 * t + j] = triangles[t][j];
            }
        }
        return indices;
    }

    private static toTriples(indices: number[]): number[][] {
        const triangles = new Array<number[]>(indices.length / 3);
        for (let t = 0; t < triangles.length; ++t) {
            triangles[t] = [indices[3 * t], indices[3 * t + 1], indices[3 * t + 2]];
        }
        return triangles;
    }
}
