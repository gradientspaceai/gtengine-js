// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) UniqueVerticesSimplices.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// UniqueVerticesSimplices allows mesh generation and elimination of duplicate
// and/or unused vertices. The mesh can be in any dimension D >= 2. In 2
// dimensions, the mesh is a collection of edges. In 3 dimensions, the mesh is
// a collection of triangles. Generally, the mesh is a collection of
// D-dimensional simplices. The following operations are supported, either for
// a mesh topology consisting of indices or of arrays, each array representing
// a simplex.
//
//   1. Generate an indexed simplex representation from an array of simplices,
//      each simplex represented by D contiguous vertices. Presumably, the
//      simplices share vertices. The output is an array of unique vertices
//      (a vertex pool) and an array of D-element arrays of indices into the
//      pool, each such array representing a simplex.
//
//   2. Remove duplicate vertices from a vertex pool used by an indexed
//      simplex representation. A new vertex pool of unique vertices is
//      generated and the indexed simplices are modified to be indices into
//      this vertex pool.
//
//   3. Remove unused vertices from a vertex pool used by an indexed simplex
//      representation. A new vertex pool of unique vertices is generated and
//      the indexed simplices are modified to be indices into the new vertex
//      pool.
//
//   4. Remove duplicate and unused vertices from a vertex pool, a combination
//      of the operations in #2 and #3.
//
// In the Geometric Tools distribution, the class is used for polygon Boolean
// operations (D = 2) and for compactifying triangle meshes (D = 3).
//
// Port notes: upstream is the class template
// UniqueVerticesSimplices<VertexType, IndexType, Dimension>. The compile-time
// Dimension becomes the runtime constructor argument 'dimension' and the
// IndexType is always 'number' (the checks that IndexType is an integral type
// other than bool have no counterpart). Upstream uses VertexType as the key
// type of a std::map, which requires a less-than comparison; JavaScript Map
// compares keys by identity, so the constructor takes a 'keyOf' function that
// maps a vertex to a primitive key with the property that two vertices are
// duplicates if and only if their keys are equal. The default keyOf handles
// primitives, arrays and objects that expose a numeric 'values' array (the
// port's Vector convention). Only the value equality induced by keyOf is
// observable in the upstream algorithm: the vertex order of the output pool
// is the order in which the vertices are first encountered, not the sorted
// order of the std::map. Upstream distinguishes the flat-index and
// array-of-simplices forms by C++ overloading on the output parameter type;
// the port uses distinct names, with the '...Grouped' suffix for the form
// whose simplices are arrays of 'dimension' indices. The upstream output
// parameters become returned object literals.

import { logAssert, logError } from './Logger.js';

// The type of a primitive key that identifies a vertex value.
export type UVSVertexKey = string | number;

// The default vertex-to-key map. It is total for the vertex types that occur
// in the library: numbers, strings, booleans, bigints, arrays (including
// typed arrays) of those, and objects with a numeric 'values' array such as
// the port's Vector. Anything else requires the caller to supply a keyOf
// function. Note that -0 and 0 produce the same key, as they do for the
// upstream std::map ordering, and that all NaN values produce the same key.
function defaultVertexKey(vertex: unknown): UVSVertexKey {
    const type = typeof vertex;
    if (type === 'number' || type === 'boolean' || type === 'bigint') {
        return type + ':' + String(vertex);
    }
    if (type === 'string') {
        return 'string:' + (vertex as string);
    }
    if (vertex instanceof Array || ArrayBuffer.isView(vertex)) {
        return 'a:' + Array.from(vertex as ArrayLike<unknown>,
            component => defaultVertexKey(component)).join('|');
    }
    if (vertex !== null && type === 'object') {
        const values = (vertex as { values?: unknown }).values;
        if (values instanceof Array || ArrayBuffer.isView(values)) {
            return 'v:' + Array.from(values as ArrayLike<unknown>,
                component => defaultVertexKey(component)).join('|');
        }
    }
    return logError('The vertex type requires a keyOf function.');
}

export class UniqueVerticesSimplices<VertexType> {
    private readonly mDimension: number;
    private readonly mKeyOf: (vertex: VertexType) => UVSVertexKey;

    constructor(dimension: number, keyOf?: (vertex: VertexType) => UVSVertexKey) {
        logAssert(dimension >= 2, 'Invalid dimension.');
        this.mDimension = dimension;
        this.mKeyOf = keyOf ?? (defaultVertexKey as (vertex: VertexType) => UVSVertexKey);
    }

    get dimension(): number {
        return this.mDimension;
    }

    // See #1 in the comments at the beginning of this file. The preconditions
    // are
    //   1. inVertices.length is a positive multiple of dimension
    // The postconditions are
    //   1. vertices has unique vertices
    //   2. indices.length = inVertices.length
    //   3. 0 <= indices[i] < vertices.length
    generateIndexedSimplices(inVertices: readonly VertexType[]):
        { vertices: VertexType[], indices: number[] } {
        logAssert(
            inVertices.length > 0 &&
            inVertices.length % this.mDimension === 0,
            'Invalid number of vertices.');

        const indices: number[] = new Array<number>(inVertices.length);
        const vertices = this.removeDuplicates(inVertices, indices);
        return { vertices, indices };
    }

    // See #1 in the comments at the beginning of this file. The preconditions
    // are
    //   1. inVertices.length is a positive multiple of dimension
    // The postconditions are
    //   1. vertices has unique vertices
    //   2. simplices.length = inVertices.length / dimension
    //   3. 0 <= simplices[s][d] < vertices.length
    generateIndexedSimplicesGrouped(inVertices: readonly VertexType[]):
        { vertices: VertexType[], simplices: number[][] } {
        const { vertices, indices } = this.generateIndexedSimplices(inVertices);
        return { vertices, simplices: this.group(indices) };
    }

    // See #2 in the comments at the beginning of the file. The preconditions
    // are
    //   1. inVertices.length is positive
    //   2. inIndices.length is a positive multiple of dimension
    //   3. 0 <= inIndices[i] < inVertices.length
    // The postconditions are
    //   1. vertices has unique vertices
    //   2. indices.length = inIndices.length
    //   3. 0 <= indices[i] < vertices.length
    removeDuplicateVertices(
        inVertices: readonly VertexType[],
        inIndices: readonly number[]):
        { vertices: VertexType[], indices: number[] } {
        logAssert(
            inVertices.length > 0,
            'Invalid number of vertices.');
        logAssert(
            inIndices.length > 0 &&
            inIndices.length % this.mDimension === 0,
            'Invalid number of indices.');
        const numVertices = inVertices.length;
        for (const index of inIndices) {
            logAssert(
                0 <= index && index < numVertices,
                'Invalid index.');
        }

        const inToOutMapping: number[] = new Array<number>(inVertices.length);
        const vertices = this.removeDuplicates(inVertices, inToOutMapping);

        const indices: number[] = new Array<number>(inIndices.length);
        for (let i = 0; i < inIndices.length; ++i) {
            indices[i] = inToOutMapping[inIndices[i]];
        }
        return { vertices, indices };
    }

    // See #2 in the comments at the beginning of the file. The preconditions
    // are
    //   1. inVertices.length is positive
    //   2. inSimplices.length is positive
    //   3. 0 <= inSimplices[s][d] < inVertices.length
    // The postconditions are
    //   1. vertices has unique vertices
    //   2. simplices.length = inSimplices.length
    //   3. 0 <= simplices[s][d] < vertices.length
    removeDuplicateVerticesGrouped(
        inVertices: readonly VertexType[],
        inSimplices: readonly (readonly number[])[]):
        { vertices: VertexType[], simplices: number[][] } {
        logAssert(
            inVertices.length > 0,
            'Invalid number of vertices.');
        logAssert(
            inSimplices.length > 0,
            'Invalid number of simplices.');
        const { vertices, indices } = this.removeDuplicateVertices(inVertices,
            this.ungroup(inSimplices));
        return { vertices, simplices: this.group(indices) };
    }

    // See #3 in the comments at the beginning of the file. The preconditions
    // are
    //   1. inVertices.length is positive
    //   2. inIndices.length is a positive multiple of dimension
    //   3. 0 <= inIndices[i] < inVertices.length
    // The postconditions are
    //   1. vertices.length is positive
    //   2. indices.length = inIndices.length
    //   3. 0 <= indices[i] < vertices.length
    //   4. each vertices[v] occurs at least once in indices[]
    removeUnusedVertices(
        inVertices: readonly VertexType[],
        inIndices: readonly number[]):
        { vertices: VertexType[], indices: number[] } {
        logAssert(
            inVertices.length > 0,
            'Invalid number of vertices.');
        logAssert(
            inIndices.length > 0 &&
            inIndices.length % this.mDimension === 0,
            'Invalid number of indices.');
        const numVertices = inVertices.length;
        for (const index of inIndices) {
            logAssert(
                0 <= index && index < numVertices,
                'Invalid index.');
        }

        const indices: number[] = new Array<number>(inIndices.length);
        const vertices = this.removeUnused(inVertices, inIndices, indices);
        return { vertices, indices };
    }

    // See #3 in the comments at the beginning of the file. The preconditions
    // are
    //   1. inVertices.length is positive
    //   2. inSimplices.length is positive
    //   3. 0 <= inSimplices[s][d] < inVertices.length
    // The postconditions are
    //   1. vertices.length is positive
    //   2. simplices.length = inSimplices.length
    //   3. 0 <= simplices[s][d] < vertices.length
    //   4. each vertices[v] occurs at least once in simplices[][]
    removeUnusedVerticesGrouped(
        inVertices: readonly VertexType[],
        inSimplices: readonly (readonly number[])[]):
        { vertices: VertexType[], simplices: number[][] } {
        logAssert(
            inVertices.length > 0,
            'Invalid number of vertices.');
        logAssert(
            inSimplices.length > 0,
            'Invalid number of simplices.');
        const { vertices, indices } = this.removeUnusedVertices(inVertices,
            this.ungroup(inSimplices));
        return { vertices, simplices: this.group(indices) };
    }

    // See #4 and the preconditions for removeDuplicateVertices and for
    // removeUnusedVertices.
    removeDuplicateAndUnusedVertices(
        inVertices: readonly VertexType[],
        inIndices: readonly number[]):
        { vertices: VertexType[], indices: number[] } {
        const temp = this.removeDuplicateVertices(inVertices, inIndices);
        return this.removeUnusedVertices(temp.vertices, temp.indices);
    }

    // See #4 and the preconditions for removeDuplicateVerticesGrouped and for
    // removeUnusedVerticesGrouped.
    removeDuplicateAndUnusedVerticesGrouped(
        inVertices: readonly VertexType[],
        inSimplices: readonly (readonly number[])[]):
        { vertices: VertexType[], simplices: number[][] } {
        const temp = this.removeDuplicateVerticesGrouped(inVertices, inSimplices);
        return this.removeUnusedVerticesGrouped(temp.vertices, temp.simplices);
    }

    // The port of the C++ reinterpret_cast between an array of simplices and
    // the flat array of indices that shares its storage.
    private group(indices: readonly number[]): number[][] {
        const simplices: number[][] = [];
        for (let i = 0; i < indices.length; i += this.mDimension) {
            simplices.push(indices.slice(i, i + this.mDimension));
        }
        return simplices;
    }

    private ungroup(simplices: readonly (readonly number[])[]): number[] {
        const indices: number[] = [];
        for (const simplex of simplices) {
            logAssert(simplex.length === this.mDimension, 'Invalid simplex.');
            for (let d = 0; d < this.mDimension; ++d) {
                indices.push(simplex[d]);
            }
        }
        return indices;
    }

    // Store in inToOutMapping[] the index in the returned vertex pool of each
    // input vertex.
    private removeDuplicates(
        inVertices: readonly VertexType[],
        inToOutMapping: number[]): VertexType[] {
        // Construct the unique vertices.
        const numInVertices = inVertices.length;
        let numOutVertices = 0;
        const outVertices: VertexType[] = [];
        const vmap = new Map<UVSVertexKey, number>();
        for (let i = 0; i < numInVertices; ++i) {
            const key = this.mKeyOf(inVertices[i]);
            const value = vmap.get(key);
            if (value !== undefined) {
                // The vertex is a duplicate of one inserted earlier into the
                // map. Its index i will be modified to that of the
                // first-found vertex.
                inToOutMapping[i] = value;
            } else {
                // The vertex occurs for the first time.
                vmap.set(key, numOutVertices);
                inToOutMapping[i] = numOutVertices;
                // Pack the unique vertices into an array. Upstream packs them
                // after the loop by iterating the map, which stores the
                // vertex at the same index computed here.
                outVertices.push(inVertices[i]);
                ++numOutVertices;
            }
        }

        return outVertices;
    }

    // Store in outIndices[] the reassigned indices and return the packed
    // vertex pool.
    private removeUnused(
        inVertices: readonly VertexType[],
        inIndices: readonly number[],
        outIndices: number[]): VertexType[] {
        // Get the unique set of used indices. The upstream container is a
        // std::set<IndexType>, which iterates in increasing order; the port
        // replicates the order by sorting.
        const usedIndices = Array.from(new Set<number>(inIndices));
        usedIndices.sort((i0, i1) => i0 - i1);

        // Locate the used vertices and pack them into an array.
        const outVertices: VertexType[] = new Array<VertexType>(usedIndices.length);
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
}
