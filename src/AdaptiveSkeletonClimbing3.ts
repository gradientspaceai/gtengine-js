// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AdaptiveSkeletonClimbing3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Extract level surfaces using an adaptive approach to reduce the triangle
// count. The implementation is for the algorithm described in the paper
//   Multiresolution Isosurface Extraction with Adaptive Skeleton Climbing
//   Tim Poston, Tien-Tsin Wong and Pheng-Ann Heng
//   Computer Graphics forum, volume 17, issue 3, September 1998
//   pages 137-147
// https://onlinelibrary.wiley.com/doi/abs/10.1111/1467-8659.00261
//
// Port notes:
//  - The upstream template parameters <T, Real> both become number. T must
//    hold integer values (upstream restricts it to int{8,16,32}_t and
//    uint{8,16,32}_t); Real is for extraction to floating-point vertices.
//    The input image is passed as ArrayLike<number>, so number[] and any of
//    the typed arrays (Int32Array, Uint8Array, ...) work; the array is
//    aliased, not copied, exactly as the upstream 'T const*'. This matches
//    the AdaptiveSkeletonClimbing2 port.
//  - The upstream int64_t determinants that disambiguate the
//    four-intersection face configuration are computed with bigint so that
//    products of 32-bit voxel values do not lose precision in IEEE doubles.
//  - Vertex (std::array<Real,3>) becomes the tuple type
//    AdaptiveSkeletonClimbing3Vertex = [number, number, number]. Triangle
//    (TriangleKey<true>) becomes TriangleKey constructed with ordered=true.
//    The unused upstream typedef 'Edge' is not ported.
//  - Out-parameters become returned values: extract() returns
//    { vertices, triangles } and computeNormals() returns the normal array.
//    makeUnique() and orientTriangles() mutate the arrays they are given, as
//    upstream does.
//  - The upstream MakeUnique packs vertices and triangles by their std::map
//    insertion indices, so first-encounter order is preserved; the port
//    replicates this with insertion-ordered Maps keyed by a canonical string.
//  - The private debugging function PrintBoxes (ostream output) is not
//    ported.
//  - The tail of the six Get{X,Y,Z}{Min,Max}EdgesM functions (the
//    subdivision-vertex insertion and the chain of edges) is character for
//    character identical in upstream apart from the coordinate index used in
//    the endpoint-swap test, so it is factored into the module-private
//    helper insertSubdividedFaceEdges. This is a behavior-preserving
//    refactor.

import { logAssert, logError } from './Logger';
import { Array2 } from './Array2';
import { TriangleKey } from './TriangleKey';

// A vertex of the extracted level surface, as (x, y, z) in voxel
// coordinates.
export type AdaptiveSkeletonClimbing3Vertex = [number, number, number];

type Vertex = AdaptiveSkeletonClimbing3Vertex;

// Configuration flags for LinearMergeTree nodes. CFG_MULT is the bitwise OR
// of CFG_INCR and CFG_DECR.
const CFG_NONE = 0;
const CFG_INCR = 1;
const CFG_DECR = 2;
const CFG_MULT = 3;
const CFG_ROOT_MASK = 3;
const CFG_EDGE = 4;
const CFG_ZERO_SUBEDGE = 8;

// Indices into the VETable vertex array. The first twelve are the box edges
// (one interpolated vertex per edge), the last six are the face branch
// points of the plus-sign configuration.
const EI_XMIN_YMIN = 0;
const EI_XMIN_YMAX = 1;
const EI_XMAX_YMIN = 2;
const EI_XMAX_YMAX = 3;
const EI_XMIN_ZMIN = 4;
const EI_XMIN_ZMAX = 5;
const EI_XMAX_ZMIN = 6;
const EI_XMAX_ZMAX = 7;
const EI_YMIN_ZMIN = 8;
const EI_YMIN_ZMAX = 9;
const EI_YMAX_ZMIN = 10;
const EI_YMAX_ZMAX = 11;
const FI_XMIN = 12;
const FI_XMAX = 13;
const FI_YMIN = 14;
const FI_YMAX = 15;
const FI_ZMIN = 16;
const FI_ZMAX = 17;
const I_QUANTITY = 18;

// The bit flags corresponding to the edge indices.
const EB_XMIN_YMIN = 1 << EI_XMIN_YMIN;
const EB_XMIN_YMAX = 1 << EI_XMIN_YMAX;
const EB_XMAX_YMIN = 1 << EI_XMAX_YMIN;
const EB_XMAX_YMAX = 1 << EI_XMAX_YMAX;
const EB_XMIN_ZMIN = 1 << EI_XMIN_ZMIN;
const EB_XMIN_ZMAX = 1 << EI_XMIN_ZMAX;
const EB_XMAX_ZMIN = 1 << EI_XMAX_ZMIN;
const EB_XMAX_ZMAX = 1 << EI_XMAX_ZMAX;
const EB_YMIN_ZMIN = 1 << EI_YMIN_ZMIN;
const EB_YMIN_ZMAX = 1 << EI_YMIN_ZMAX;
const EB_YMAX_ZMIN = 1 << EI_YMAX_ZMIN;
const EB_YMAX_ZMAX = 1 << EI_YMAX_ZMAX;

// The error message that upstream uses when a face configuration is
// impossible for a non-integer level value.
const NON_INTEGER_LEVEL_MESSAGE = 'The level value cannot be an exact integer.';

// A monobox produced by the merging phase. (x0,y0,z0) is the minimum corner,
// (dx,dy,dz) the extents and (LX,LY,LZ) the linear-merge-tree node indices
// of the box edges in the x, y and z directions.
class OctBox {
    x0: number;
    y0: number;
    z0: number;
    x1: number;
    y1: number;
    z1: number;
    dx: number;
    dy: number;
    dz: number;
    LX: number;
    LY: number;
    LZ: number;

    constructor(inX0: number, inY0: number, inZ0: number,
        inDX: number, inDY: number, inDZ: number,
        inLX: number, inLY: number, inLZ: number) {
        this.x0 = inX0;
        this.y0 = inY0;
        this.z0 = inZ0;
        this.x1 = inX0 + inDX;
        this.y1 = inY0 + inDY;
        this.z1 = inZ0 + inDZ;
        this.dx = inDX;
        this.dy = inDY;
        this.dz = inDZ;
        this.LX = inLX;
        this.LY = inLY;
        this.LZ = inLZ;
    }
}

// The mutable per-octant merge state used while combining the eight children
// of an octree node.
class MergeBox {
    xStride: number;
    yStride: number;
    zStride: number;
    valid: boolean;

    constructor(stride: number) {
        this.xStride = stride;
        this.yStride = stride;
        this.zStride = stride;
        this.valid = true;
    }
}

// A complete binary tree over the 2^N intervals of one grid line. Each node
// stores the combined sign-change configuration of its subtree and the index
// of the interval containing the (single) sign change when there is one.
class LinearMergeTree {
    private mTwoPowerN: number;
    private mNodes: number[];
    private mZeroBases: number[];

    constructor(N: number) {
        this.mTwoPowerN = 1 << N;
        this.mNodes = new Array<number>(2 * this.mTwoPowerN - 1).fill(0);
        this.mZeroBases = new Array<number>(2 * this.mTwoPowerN - 1).fill(0);
    }

    isNone(i: number): boolean {
        return (this.mNodes[i] & CFG_ROOT_MASK) === CFG_NONE;
    }

    getRootType(i: number): number {
        return this.mNodes[i] & CFG_ROOT_MASK;
    }

    getZeroBase(i: number): number {
        return this.mZeroBases[i];
    }

    setEdge(i: number): void {
        this.mNodes[i] |= CFG_EDGE;

        // Inform all predecessors whether they have a zero subedge.
        if (this.mZeroBases[i] >= 0) {
            while (i > 0) {
                i = Math.trunc((i - 1) / 2);
                this.mNodes[i] |= CFG_ZERO_SUBEDGE;
            }
        }
    }

    isZeroEdge(i: number): boolean {
        return this.mNodes[i] === (CFG_EDGE | CFG_INCR)
            || this.mNodes[i] === (CFG_EDGE | CFG_DECR);
    }

    hasZeroSubedge(i: number): boolean {
        return (this.mNodes[i] & CFG_ZERO_SUBEDGE) !== 0;
    }

    setLevel(level: number, data: ArrayLike<number>, offset: number, stride: number): void {
        // Assert: The 'level' is not an image value. Because the image values
        // are integers, choose 'level' to be a number that does not represent
        // an integer.

        // Determine the sign changes between pairs of consecutive samples.
        const firstLeaf = this.mTwoPowerN - 1;
        for (let i = 0, leaf = firstLeaf; i < this.mTwoPowerN; ++i, ++leaf) {
            const base = offset + stride * i;
            const value0 = data[base];
            const value1 = data[base + stride];

            if (value0 > level) {
                if (value1 > level) {
                    this.mNodes[leaf] = CFG_NONE;
                    this.mZeroBases[leaf] = -1;
                } else {
                    this.mNodes[leaf] = CFG_DECR;
                    this.mZeroBases[leaf] = i;
                }
            } else { // value0 < level
                if (value1 > level) {
                    this.mNodes[leaf] = CFG_INCR;
                    this.mZeroBases[leaf] = i;
                } else {
                    this.mNodes[leaf] = CFG_NONE;
                    this.mZeroBases[leaf] = -1;
                }
            }
        }

        // Propagate the sign change information up the binary tree.
        for (let i = firstLeaf - 1; i >= 0; --i) {
            const twoIp1 = 2 * i + 1;
            const twoIp2 = twoIp1 + 1;
            const value0 = this.mNodes[twoIp1];
            const value1 = this.mNodes[twoIp2];
            const combine = (value0 | value1);
            this.mNodes[i] = combine;
            if (combine === CFG_INCR) {
                if (value0 === CFG_INCR) {
                    this.mZeroBases[i] = this.mZeroBases[twoIp1];
                } else {
                    this.mZeroBases[i] = this.mZeroBases[twoIp2];
                }
            } else if (combine === CFG_DECR) {
                if (value0 === CFG_DECR) {
                    this.mZeroBases[i] = this.mZeroBases[twoIp1];
                } else {
                    this.mZeroBases[i] = this.mZeroBases[twoIp2];
                }
            } else {
                this.mZeroBases[i] = -1;
            }
        }
    }
}

// One node of the vertex-edge wireframe assembled for a single monobox.
class TVertex {
    position: Vertex;
    adjQuantity: number;
    adjacent: [number, number, number, number];
    valid: boolean;

    constructor(inPosition?: Vertex) {
        if (inPosition === undefined) {
            this.position = [0, 0, 0];
            this.valid = false;
        } else {
            this.position = inPosition;
            this.valid = true;
        }
        this.adjQuantity = 0;
        this.adjacent = [0, 0, 0, 0];
    }
}

// The vertex-edge table for a single monobox. The first I_QUANTITY entries
// are reserved for the box-edge and face-branch-point vertices; additional
// entries are the face subdivision vertices of merged boxes.
class VETable {
    private mVertices: TVertex[];

    constructor() {
        this.mVertices = new Array<TVertex>(I_QUANTITY);
        for (let i = 0; i < I_QUANTITY; ++i) {
            this.mVertices[i] = new TVertex();
        }
    }

    isValidVertex(i: number): boolean {
        return this.mVertices[i].valid;
    }

    getNumVertices(): number {
        return this.mVertices.length;
    }

    getVertex(i: number): Vertex {
        return this.mVertices[i].position;
    }

    // The upstream Insert(i, x, y, z) overload.
    insert(i: number, x: number, y: number, z: number): void {
        const vertex = this.mVertices[i];
        vertex.position = [x, y, z];
        vertex.valid = true;
    }

    // The upstream Insert(position) overload, which appends a new vertex.
    append(position: Vertex): void {
        this.mVertices.push(new TVertex(position));
    }

    insertEdge(v0: number, v1: number): void {
        const vertex0 = this.mVertices[v0];
        const vertex1 = this.mVertices[v1];
        vertex0.adjacent[vertex0.adjQuantity] = v1;
        ++vertex0.adjQuantity;
        vertex1.adjacent[vertex1.adjQuantity] = v0;
        ++vertex1.adjQuantity;
    }

    removeTrianglesEC(positions: Vertex[], triangles: TriangleKey[]): void {
        // Ear-clip the wireframe to get the triangles.
        const triangle: [number, number, number] = [0, 0, 0];
        while (this.removeEC(triangle)) {
            const v0 = positions.length;
            const v1 = v0 + 1;
            const v2 = v1 + 1;
            // Bypassing the TriangleKey ordering, as upstream does; the
            // ordered case would produce the same triple because
            // v0 < v1 < v2.
            const tkey = new TriangleKey(true);
            tkey.V[0] = v0;
            tkey.V[1] = v1;
            tkey.V[2] = v2;
            triangles.push(tkey);
            positions.push(this.mVertices[triangle[0]].position);
            positions.push(this.mVertices[triangle[1]].position);
            positions.push(this.mVertices[triangle[2]].position);
        }
    }

    removeTrianglesSE(positions: Vertex[], triangles: TriangleKey[]): void {
        // Compute centroid of vertices.
        const centroid: Vertex = [0, 0, 0];
        const vmax = this.mVertices.length;
        let quantity = 0;
        for (let i = 0; i < vmax; ++i) {
            const vertex = this.mVertices[i];
            if (vertex.valid) {
                for (let j = 0; j < 3; ++j) {
                    centroid[j] += vertex.position[j];
                }
                ++quantity;
            }
        }
        for (let j = 0; j < 3; ++j) {
            centroid[j] /= quantity;
        }

        const v0 = positions.length;
        positions.push(centroid);

        let i1 = I_QUANTITY;
        let v1 = v0 + 1;
        positions.push(this.mVertices[i1].position);

        let i2 = this.mVertices[i1].adjacent[1];
        let v2 = 0;
        for (let i = 0; i < quantity - 1; ++i) {
            v2 = v1 + 1;
            positions.push(this.mVertices[i2].position);
            // Bypassing the TriangleKey ordering, as upstream does; the
            // ordered case would produce the same triple because
            // v0 < v1 < v2.
            const tkey = new TriangleKey(true);
            tkey.V[0] = v0;
            tkey.V[1] = v1;
            tkey.V[2] = v2;
            triangles.push(tkey);
            if (this.mVertices[i2].adjacent[1] !== i1) {
                i1 = i2;
                i2 = this.mVertices[i2].adjacent[1];
            } else {
                i1 = i2;
                i2 = this.mVertices[i2].adjacent[0];
            }
            v1 = v2;
        }

        v2 = v0 + 1;
        // Unlike the previous constructions, it is not guaranteed that
        // v0 < v1, so the ordering constructor is used.
        triangles.push(new TriangleKey(true, v0, v1, v2));
    }

    private removeVertex(i: number): void {
        const vertex0 = this.mVertices[i];
        const a0 = vertex0.adjacent[0];
        const a1 = vertex0.adjacent[1];
        const adjVertex0 = this.mVertices[a0];
        const adjVertex1 = this.mVertices[a1];

        for (let j = 0; j < adjVertex0.adjQuantity; ++j) {
            if (adjVertex0.adjacent[j] === i) {
                adjVertex0.adjacent[j] = a1;
                break;
            }
        }

        for (let j = 0; j < adjVertex1.adjQuantity; ++j) {
            if (adjVertex1.adjacent[j] === i) {
                adjVertex1.adjacent[j] = a0;
                break;
            }
        }

        vertex0.valid = false;

        if (adjVertex0.adjQuantity === 2) {
            if (adjVertex0.adjacent[0] === adjVertex0.adjacent[1]) {
                adjVertex0.valid = false;
            }
        }

        if (adjVertex1.adjQuantity === 2) {
            if (adjVertex1.adjacent[0] === adjVertex1.adjacent[1]) {
                adjVertex1.valid = false;
            }
        }
    }

    // Ear clipping. The triple is filled with table indices, not with mesh
    // indices, so the upstream TriangleKey ordering is deliberately bypassed.
    private removeEC(triangle: [number, number, number]): boolean {
        const numVertices = this.mVertices.length;
        for (let i = 0; i < numVertices; ++i) {
            const vertex = this.mVertices[i];
            if (vertex.valid && vertex.adjQuantity === 2) {
                triangle[0] = i;
                triangle[1] = vertex.adjacent[0];
                triangle[2] = vertex.adjacent[1];
                this.removeVertex(i);
                return true;
            }
        }

        return false;
    }
}

// The three upstream comparators for the std::set of face subdivision
// vertices. Each is a lexicographic order on a rotation of the coordinates.
// Sort0 orders by (x,y,z), Sort1 by (z,x,y) and Sort2 by (y,z,x).
function sortByCoords(a: Vertex, b: Vertex, k0: number, k1: number, k2: number): number {
    if (a[k0] !== b[k0]) {
        return a[k0] < b[k0] ? -1 : 1;
    }
    if (a[k1] !== b[k1]) {
        return a[k1] < b[k1] ? -1 : 1;
    }
    if (a[k2] !== b[k2]) {
        return a[k2] < b[k2] ? -1 : 1;
    }
    return 0;
}

const sort0 = (a: Vertex, b: Vertex): number => sortByCoords(a, b, 0, 1, 2);
const sort1 = (a: Vertex, b: Vertex): number => sortByCoords(a, b, 2, 0, 1);
const sort2 = (a: Vertex, b: Vertex): number => sortByCoords(a, b, 1, 2, 0);

// The std::set<Vertex, SortK> of the upstream face-subdivision code: unique
// vertices in the comparator's order. Because each comparator is a total
// order on the (x,y,z) triple, set equivalence is exact triple equality.
class VertexSet {
    private mMap: Map<string, Vertex>;
    private mCompare: (a: Vertex, b: Vertex) => number;

    constructor(compare: (a: Vertex, b: Vertex) => number) {
        this.mMap = new Map<string, Vertex>();
        this.mCompare = compare;
    }

    insert(position: Vertex): void {
        const key = position[0] + ',' + position[1] + ',' + position[2];
        if (!this.mMap.has(key)) {
            this.mMap.set(key, position);
        }
    }

    sorted(): Vertex[] {
        const values = Array.from(this.mMap.values());
        values.sort(this.mCompare);
        return values;
    }
}

// The common tail of the six Get{X,Y,Z}{Min,Max}EdgesM functions. The face
// polyline runs from the box-edge vertex 'end0' to 'end1' through the
// subdivision vertices 'vSet' (already sorted); 'k' is the coordinate along
// which the endpoints and the subdivision vertices are compared to decide
// whether the endpoints must be swapped.
function insertSubdividedFaceEdges(table: VETable, end0: number, end1: number,
    vSet: Vertex[], k: number): void {
    // Add subdivision.
    let v0 = table.getNumVertices();
    let v1 = v0;
    if (vSet.length === 0) {
        table.insertEdge(end0, end1);
        return;
    }

    const vk0 = Math.floor(vSet[0][k]);
    const vk1 = Math.floor(vSet[vSet.length - 1][k]);
    const ek0 = Math.floor(table.getVertex(end0)[k]);
    const ek1 = Math.floor(table.getVertex(end1)[k]);
    if (ek1 <= vk0 && vk1 <= ek0) {
        const save = end0;
        end0 = end1;
        end1 = save;
    }

    // Add vertices.
    for (const position of vSet) {
        table.append(position);
    }

    // Add edges.
    table.insertEdge(end0, v1);
    ++v1;
    const imax = vSet.length;
    for (let i = 1; i < imax; ++i, ++v0, ++v1) {
        table.insertEdge(v0, v1);
    }
    table.insertEdge(v0, end1);
}

export class AdaptiveSkeletonClimbing3 {
    // Image data.
    private mTwoPowerN: number;
    private mSize: number;
    private mSizeSqr: number;
    private mInputVoxels: ArrayLike<number>;
    private mLevel: number;

    private mFixBoundary: boolean;

    // Trees for linear merging. mXMerge is indexed by (y,z), mYMerge by
    // (x,z) and mZMerge by (x,y). Upstream accesses them as 'a[i1][i0]';
    // per the Array2 port that is 'a.get(i0, i1)', so the accessors below
    // reverse the index order.
    private mXMerge: Array2<LinearMergeTree>;
    private mYMerge: Array2<LinearMergeTree>;
    private mZMerge: Array2<LinearMergeTree>;

    // Monoboxes.
    private mBoxes: OctBox[];

    // Construction. The input image is assumed to contain
    // (2^N+1)-by-(2^N+1)-by-(2^N+1) integer-valued elements where N >= 1.
    // The organization is lexicographic order for (x,y,z). When
    // 'fixBoundary' is true, image boundary voxels are not allowed to merge
    // with any other voxels. This forces highest level of detail on the
    // boundary. The idea is that an image too large to process by itself can
    // be partitioned into smaller subimages and the adaptive skeleton
    // climbing applied to each subimage. By forcing highest resolution on
    // the boundary, adjacent subimages will not have any cracking problems.
    constructor(N: number, inputVoxels: ArrayLike<number>, fixBoundary: boolean = false) {
        this.mTwoPowerN = 1 << N;
        this.mSize = this.mTwoPowerN + 1;
        this.mSizeSqr = this.mSize * this.mSize;
        this.mInputVoxels = inputVoxels;
        this.mLevel = 0;
        this.mFixBoundary = fixBoundary;
        this.mBoxes = [];

        if (N <= 0 || inputVoxels == null) {
            logError('Invalid input.');
        }

        this.mXMerge = new Array2<LinearMergeTree>(this.mSize, this.mSize);
        this.mYMerge = new Array2<LinearMergeTree>(this.mSize, this.mSize);
        this.mZMerge = new Array2<LinearMergeTree>(this.mSize, this.mSize);
        for (let i = 0; i < this.mSize; ++i) {
            for (let j = 0; j < this.mSize; ++j) {
                this.mXMerge.set(j, i, new LinearMergeTree(N));
                this.mYMerge.set(j, i, new LinearMergeTree(N));
                this.mZMerge.set(j, i, new LinearMergeTree(N));
            }
        }
    }

    // The x-direction merge tree of the grid line through (y, z).
    private xMerge(y: number, z: number): LinearMergeTree {
        return this.mXMerge.get(z, y);
    }

    // The y-direction merge tree of the grid line through (x, z).
    private yMerge(x: number, z: number): LinearMergeTree {
        return this.mYMerge.get(z, x);
    }

    // The z-direction merge tree of the grid line through (x, y).
    private zMerge(x: number, y: number): LinearMergeTree {
        return this.mZMerge.get(y, x);
    }

    // Extract the level surface for the specified 'level', which must not be
    // an integer (the image values are integers, so a non-integer level
    // guarantees the surface misses the sample points). The 'depth' controls
    // the resolution: merging is prevented for the top 'depth' levels of the
    // octree, so depth <= 0 is the fully adaptive (coarsest) extraction and
    // depth = N is the full-resolution extraction. The returned triangles
    // index the returned vertices, which contain duplicates until
    // makeUnique() is called.
    extract(level: number, depth: number):
        { vertices: Vertex[], triangles: TriangleKey[] } {
        const localVertices: Vertex[] = [];
        const localTriangles: TriangleKey[] = [];
        this.mBoxes.length = 0;

        this.mLevel = level;
        this.merge(depth);
        this.tessellate(localVertices, localTriangles);

        return { vertices: localVertices, triangles: localTriangles };
    }

    // Remove duplicate vertices and triangles, remapping the triangle
    // indices. The input arrays are modified in place, as upstream.
    makeUnique(vertices: Vertex[], triangles: TriangleKey[]): void {
        const numVertices = vertices.length;
        const numTriangles = triangles.length;
        if (numVertices === 0 || numTriangles === 0) {
            return;
        }

        // Compute the map of unique vertices and assign to them new and
        // unique indices.
        const vmap = new Map<string, { vertex: Vertex, index: number }>();
        let nextVertex = 0;
        for (let v = 0; v < numVertices; ++v) {
            // Keep only unique vertices.
            const key = vertexKey(vertices[v]);
            if (!vmap.has(key)) {
                vmap.set(key, { vertex: vertices[v], index: nextVertex });
                ++nextVertex;
            }
        }

        // Compute the map of unique triangles and assign to them new and
        // unique indices.
        const tmap = new Map<string, { triangle: TriangleKey, index: number }>();
        let nextTriangle = 0;
        for (let t = 0; t < numTriangles; ++t) {
            const triangle = triangles[t];
            for (let i = 0; i < 3; ++i) {
                const element = vmap.get(vertexKey(vertices[triangle.V[i]]));
                logAssert(element !== undefined, 'Expecting the vertex to be in the vmap.');
                triangle.V[i] = element.index;
            }

            // Keep only unique triangles.
            const key = triangle.V[0] + ',' + triangle.V[1] + ',' + triangle.V[2];
            if (!tmap.has(key)) {
                tmap.set(key, { triangle, index: nextTriangle });
                ++nextTriangle;
            }
        }

        // Pack the vertices into an array.
        vertices.length = vmap.size;
        for (const element of vmap.values()) {
            vertices[element.index] = element.vertex;
        }

        // Pack the triangles into an array.
        triangles.length = tmap.size;
        for (const element of tmap.values()) {
            triangles[element.index] = element.triangle;
        }
    }

    // Reorder the triangles so that their normals are consistent with the
    // image gradient. When 'sameDir' is true the normals point in the
    // direction of increasing image values. The triangles are modified in
    // place, as upstream.
    orientTriangles(vertices: Vertex[], triangles: TriangleKey[], sameDir: boolean): void {
        for (const triangle of triangles) {
            // Get the triangle vertices.
            const v0 = vertices[triangle.V[0]];
            const v1 = vertices[triangle.V[1]];
            const v2 = vertices[triangle.V[2]];

            // Construct the triangle normal based on the current
            // orientation.
            const normal = triangleNormal(v0, v1, v2);

            // Get the image gradient at the vertices.
            const grad0 = this.getGradient(v0);
            const grad1 = this.getGradient(v1);
            const grad2 = this.getGradient(v2);

            // Compute the average gradient.
            const gradAvr: Vertex = [0, 0, 0];
            for (let i = 0; i < 3; ++i) {
                gradAvr[i] = (grad0[i] + grad1[i] + grad2[i]) / 3;
            }

            // Compute the dot product of normal and average gradient.
            const dot = gradAvr[0] * normal[0] + gradAvr[1] * normal[1] + gradAvr[2] * normal[2];

            // Choose triangle orientation based on gradient direction.
            if (sameDir) {
                if (dot < 0) {
                    // Wrong orientation, reorder it.
                    const save = triangle.V[1];
                    triangle.V[1] = triangle.V[2];
                    triangle.V[2] = save;
                }
            } else {
                if (dot > 0) {
                    // Wrong orientation, reorder it.
                    const save = triangle.V[1];
                    triangle.V[1] = triangle.V[2];
                    triangle.V[2] = save;
                }
            }
        }
    }

    // Compute a vertex normal to be the normalized area-weighted sum of the
    // normals of the triangles that share that vertex.
    computeNormals(vertices: Vertex[], triangles: TriangleKey[]): Vertex[] {
        const normals: Vertex[] = new Array<Vertex>(vertices.length);
        for (let i = 0; i < vertices.length; ++i) {
            normals[i] = [0, 0, 0];
        }

        for (const triangle of triangles) {
            // Get the triangle vertices and construct the triangle normal.
            const normal = triangleNormal(vertices[triangle.V[0]],
                vertices[triangle.V[1]], vertices[triangle.V[2]]);

            // Maintain the sum of normals at each vertex.
            for (let i = 0; i < 3; ++i) {
                for (let j = 0; j < 3; ++j) {
                    normals[triangle.V[i]][j] += normal[j];
                }
            }
        }

        // The normal vector storage was used to accumulate the sum of
        // triangle normals. Now these vectors must be rescaled to be unit
        // length.
        for (const normal of normals) {
            const sqrLength = normal[0] * normal[0] + normal[1] * normal[1]
                + normal[2] * normal[2];
            const length = Math.sqrt(sqrLength);
            if (length > 0) {
                for (let i = 0; i < 3; ++i) {
                    normal[i] /= length;
                }
            } else {
                for (let i = 0; i < 3; ++i) {
                    normal[i] = 0;
                }
            }
        }

        return normals;
    }

    // The monoboxes produced by the last extract() call. Exposed in place of
    // the upstream debugging function PrintBoxes.
    getNumBoxes(): number {
        return this.mBoxes.length;
    }

    // ------------------------------------------------------------------
    // Support for merging monoboxes.
    // ------------------------------------------------------------------

    private merge(depth: number): void {
        for (let y = 0; y < this.mSize; ++y) {
            for (let z = 0; z < this.mSize; ++z) {
                const offset = this.mSize * (y + this.mSize * z);
                const stride = 1;
                this.xMerge(y, z).setLevel(this.mLevel, this.mInputVoxels, offset, stride);
            }
        }

        for (let x = 0; x < this.mSize; ++x) {
            for (let z = 0; z < this.mSize; ++z) {
                const offset = x + this.mSizeSqr * z;
                const stride = this.mSize;
                this.yMerge(x, z).setLevel(this.mLevel, this.mInputVoxels, offset, stride);
            }
        }

        for (let x = 0; x < this.mSize; ++x) {
            for (let y = 0; y < this.mSize; ++y) {
                const offset = x + this.mSize * y;
                const stride = this.mSizeSqr;
                this.zMerge(x, y).setLevel(this.mLevel, this.mInputVoxels, offset, stride);
            }
        }

        this.mergeBoxes(0, 0, 0, 0, 0, 0, 0, this.mTwoPowerN, depth);
    }

    // The recursive octree merge. The return value indicates that the
    // subtree merged into a single box, which the caller may merge further.
    private mergeBoxes(v: number, LX: number, LY: number, LZ: number,
        x0: number, y0: number, z0: number, stride: number, depth: number): boolean {
        if (stride > 1) { // internal nodes
            const hStride = Math.trunc(stride / 2);
            const vBase = 8 * v;
            const v000 = vBase + 1;
            const v100 = vBase + 2;
            const v010 = vBase + 3;
            const v110 = vBase + 4;
            const v001 = vBase + 5;
            const v101 = vBase + 6;
            const v011 = vBase + 7;
            const v111 = vBase + 8;
            const LX0 = 2 * LX + 1;
            const LX1 = LX0 + 1;
            const LY0 = 2 * LY + 1;
            const LY1 = LY0 + 1;
            const LZ0 = 2 * LZ + 1;
            const LZ1 = LZ0 + 1;
            const x1 = x0 + hStride;
            const y1 = y0 + hStride;
            const z1 = z0 + hStride;

            const dm1 = depth - 1;
            const m000 = this.mergeBoxes(v000, LX0, LY0, LZ0, x0, y0, z0, hStride, dm1);
            const m100 = this.mergeBoxes(v100, LX1, LY0, LZ0, x1, y0, z0, hStride, dm1);
            const m010 = this.mergeBoxes(v010, LX0, LY1, LZ0, x0, y1, z0, hStride, dm1);
            const m110 = this.mergeBoxes(v110, LX1, LY1, LZ0, x1, y1, z0, hStride, dm1);
            const m001 = this.mergeBoxes(v001, LX0, LY0, LZ1, x0, y0, z1, hStride, dm1);
            const m101 = this.mergeBoxes(v101, LX1, LY0, LZ1, x1, y0, z1, hStride, dm1);
            const m011 = this.mergeBoxes(v011, LX0, LY1, LZ1, x0, y1, z1, hStride, dm1);
            const m111 = this.mergeBoxes(v111, LX1, LY1, LZ1, x1, y1, z1, hStride, dm1);

            const r000 = new MergeBox(hStride);
            const r100 = new MergeBox(hStride);
            const r010 = new MergeBox(hStride);
            const r110 = new MergeBox(hStride);
            const r001 = new MergeBox(hStride);
            const r101 = new MergeBox(hStride);
            const r011 = new MergeBox(hStride);
            const r111 = new MergeBox(hStride);

            if (depth <= 0) {
                if (m000 && m001) {
                    this.doZMerge(r000, r001, x0, y0, LZ);
                }

                if (m100 && m101) {
                    this.doZMerge(r100, r101, x1, y0, LZ);
                }

                if (m010 && m011) {
                    this.doZMerge(r010, r011, x0, y1, LZ);
                }

                if (m110 && m111) {
                    this.doZMerge(r110, r111, x1, y1, LZ);
                }

                if (m000 && m010) {
                    this.doYMerge(r000, r010, x0, LY, z0);
                }

                if (m100 && m110) {
                    this.doYMerge(r100, r110, x1, LY, z0);
                }

                if (m001 && m011) {
                    this.doYMerge(r001, r011, x0, LY, z1);
                }

                if (m101 && m111) {
                    this.doYMerge(r101, r111, x1, LY, z1);
                }

                if (m000 && m100) {
                    this.doXMerge(r000, r100, LX, y0, z0);
                }

                if (m010 && m110) {
                    this.doXMerge(r010, r110, LX, y1, z0);
                }

                if (m001 && m101) {
                    this.doXMerge(r001, r101, LX, y0, z1);
                }

                if (m011 && m111) {
                    this.doXMerge(r011, r111, LX, y1, z1);
                }
            }

            if (depth <= 1) {
                if (r000.valid) {
                    if (r000.xStride === stride) {
                        if (r000.yStride === stride) {
                            if (r000.zStride === stride) {
                                return true;
                            } else {
                                this.addBox(x0, y0, z0, stride, stride, hStride, LX, LY, LZ0);
                            }
                        } else {
                            if (r000.zStride === stride) {
                                this.addBox(x0, y0, z0, stride, hStride, stride, LX, LY0, LZ);
                            } else {
                                this.addBox(x0, y0, z0, stride, hStride, hStride, LX, LY0, LZ0);
                            }
                        }
                    } else {
                        if (r000.yStride === stride) {
                            if (r000.zStride === stride) {
                                this.addBox(x0, y0, z0, hStride, stride, stride, LX0, LY, LZ);
                            } else {
                                this.addBox(x0, y0, z0, hStride, stride, hStride, LX0, LY, LZ0);
                            }
                        } else {
                            if (r000.zStride === stride) {
                                this.addBox(x0, y0, z0, hStride, hStride, stride, LX0, LY0, LZ);
                            } else if (m000) {
                                this.addBox(x0, y0, z0, hStride, hStride, hStride, LX0, LY0, LZ0);
                            }
                        }
                    }
                }

                if (r100.valid) {
                    if (r100.yStride === stride) {
                        if (r100.zStride === stride) {
                            this.addBox(x1, y0, z0, hStride, stride, stride, LX1, LY, LZ);
                        } else {
                            this.addBox(x1, y0, z0, hStride, stride, hStride, LX1, LY, LZ0);
                        }
                    } else {
                        if (r100.zStride === stride) {
                            this.addBox(x1, y0, z0, hStride, hStride, stride, LX1, LY0, LZ);
                        } else if (m100) {
                            this.addBox(x1, y0, z0, hStride, hStride, hStride, LX1, LY0, LZ0);
                        }
                    }
                }

                if (r010.valid) {
                    if (r010.xStride === stride) {
                        if (r010.zStride === stride) {
                            this.addBox(x0, y1, z0, stride, hStride, stride, LX, LY1, LZ);
                        } else {
                            this.addBox(x0, y1, z0, stride, hStride, hStride, LX, LY1, LZ0);
                        }
                    } else {
                        if (r010.zStride === stride) {
                            this.addBox(x0, y1, z0, hStride, hStride, stride, LX0, LY1, LZ);
                        } else if (m010) {
                            this.addBox(x0, y1, z0, hStride, hStride, hStride, LX0, LY1, LZ0);
                        }
                    }
                }

                if (r001.valid) {
                    if (r001.xStride === stride) {
                        if (r001.yStride === stride) {
                            this.addBox(x0, y0, z1, stride, stride, hStride, LX, LY, LZ1);
                        } else {
                            this.addBox(x0, y0, z1, stride, hStride, hStride, LX, LY0, LZ1);
                        }
                    } else {
                        if (r001.yStride === stride) {
                            this.addBox(x0, y0, z1, hStride, stride, hStride, LX0, LY, LZ1);
                        } else if (m001) {
                            this.addBox(x0, y0, z1, hStride, hStride, hStride, LX0, LY0, LZ1);
                        }
                    }
                }

                if (r110.valid) {
                    if (r110.zStride === stride) {
                        this.addBox(x1, y1, z0, hStride, hStride, stride, LX1, LY1, LZ);
                    } else if (m110) {
                        this.addBox(x1, y1, z0, hStride, hStride, hStride, LX1, LY1, LZ0);
                    }
                }

                if (r101.valid) {
                    if (r101.yStride === stride) {
                        this.addBox(x1, y0, z1, hStride, stride, hStride, LX1, LY, LZ1);
                    } else if (m101) {
                        this.addBox(x1, y0, z1, hStride, hStride, hStride, LX1, LY0, LZ1);
                    }
                }

                if (r011.valid) {
                    if (r011.xStride === stride) {
                        this.addBox(x0, y1, z1, stride, hStride, hStride, LX, LY1, LZ1);
                    } else if (m011) {
                        this.addBox(x0, y1, z1, hStride, hStride, hStride, LX0, LY1, LZ1);
                    }
                }

                if (r111.valid && m111) {
                    this.addBox(x1, y1, z1, hStride, hStride, hStride, LX1, LY1, LZ1);
                }
            }
            return false;
        } else { // leaf nodes
            if (this.mFixBoundary) {
                // Do not allow boundary voxels to merge with any other
                // voxels.
                this.addBox(x0, y0, z0, 1, 1, 1, LX, LY, LZ);
                return false;
            }

            // A leaf box is mergeable with neighbors as long as all its
            // faces have 0 or 2 sign changes on the edges. That is, a face
            // may not have sign changes on all four edges. If it does, the
            // resulting box for tessellating is 1x1x1 and is handled
            // separately from boxes of larger dimensions.

            // xmin face
            const z1 = z0 + 1;
            let rt0 = this.yMerge(x0, z0).getRootType(LY);
            let rt1 = this.yMerge(x0, z1).getRootType(LY);
            if ((rt0 | rt1) === CFG_MULT) {
                this.addBox(x0, y0, z0, 1, 1, 1, LX, LY, LZ);
                return false;
            }

            // xmax face
            const x1 = x0 + 1;
            rt0 = this.yMerge(x1, z0).getRootType(LY);
            rt1 = this.yMerge(x1, z1).getRootType(LY);
            if ((rt0 | rt1) === CFG_MULT) {
                this.addBox(x0, y0, z0, 1, 1, 1, LX, LY, LZ);
                return false;
            }

            // ymin face
            rt0 = this.zMerge(x0, y0).getRootType(LZ);
            rt1 = this.zMerge(x1, y0).getRootType(LZ);
            if ((rt0 | rt1) === CFG_MULT) {
                this.addBox(x0, y0, z0, 1, 1, 1, LX, LY, LZ);
                return false;
            }

            // ymax face
            const y1 = y0 + 1;
            rt0 = this.zMerge(x0, y1).getRootType(LZ);
            rt1 = this.zMerge(x1, y1).getRootType(LZ);
            if ((rt0 | rt1) === CFG_MULT) {
                this.addBox(x0, y0, z0, 1, 1, 1, LX, LY, LZ);
                return false;
            }

            // zmin face
            rt0 = this.xMerge(y0, z0).getRootType(LX);
            rt1 = this.xMerge(y1, z0).getRootType(LX);
            if ((rt0 | rt1) === CFG_MULT) {
                this.addBox(x0, y0, z0, 1, 1, 1, LX, LY, LZ);
                return false;
            }

            // zmax face
            rt0 = this.xMerge(y0, z1).getRootType(LX);
            rt1 = this.xMerge(y1, z1).getRootType(LX);
            if ((rt0 | rt1) === CFG_MULT) {
                this.addBox(x0, y0, z0, 1, 1, 1, LX, LY, LZ);
                return false;
            }

            return true;
        }
    }

    private doXMerge(r0: MergeBox, r1: MergeBox, LX: number, y0: number, z0: number): boolean {
        if (!r0.valid || !r1.valid || r0.yStride !== r1.yStride || r0.zStride !== r1.zStride) {
            return false;
        }

        // Boxes are potentially x-mergeable.
        const y1 = y0 + r0.yStride;
        const z1 = z0 + r0.zStride;
        let incr = 0;
        let decr = 0;
        for (let y = y0; y <= y1; ++y) {
            for (let z = z0; z <= z1; ++z) {
                const rootType = this.xMerge(y, z).getRootType(LX);
                if (rootType === CFG_MULT) {
                    return false;
                } else if (rootType === CFG_INCR) {
                    ++incr;
                } else if (rootType === CFG_DECR) {
                    ++decr;
                }
            }
        }

        if (incr !== 0 && decr !== 0) {
            return false;
        }

        // Strongly mono, x-merge the boxes.
        r0.xStride *= 2;
        r1.valid = false;
        return true;
    }

    private doYMerge(r0: MergeBox, r1: MergeBox, x0: number, LY: number, z0: number): boolean {
        if (!r0.valid || !r1.valid || r0.xStride !== r1.xStride || r0.zStride !== r1.zStride) {
            return false;
        }

        // Boxes are potentially y-mergeable.
        const x1 = x0 + r0.xStride;
        const z1 = z0 + r0.zStride;
        let incr = 0;
        let decr = 0;
        for (let x = x0; x <= x1; ++x) {
            for (let z = z0; z <= z1; ++z) {
                const rootType = this.yMerge(x, z).getRootType(LY);
                if (rootType === CFG_MULT) {
                    return false;
                } else if (rootType === CFG_INCR) {
                    ++incr;
                } else if (rootType === CFG_DECR) {
                    ++decr;
                }
            }
        }

        if (incr !== 0 && decr !== 0) {
            return false;
        }

        // Strongly mono, y-merge the boxes.
        r0.yStride *= 2;
        r1.valid = false;
        return true;
    }

    private doZMerge(r0: MergeBox, r1: MergeBox, x0: number, y0: number, LZ: number): boolean {
        if (!r0.valid || !r1.valid || r0.xStride !== r1.xStride || r0.yStride !== r1.yStride) {
            return false;
        }

        // Boxes are potentially z-mergeable.
        const x1 = x0 + r0.xStride;
        const y1 = y0 + r0.yStride;
        let incr = 0;
        let decr = 0;
        for (let x = x0; x <= x1; ++x) {
            for (let y = y0; y <= y1; ++y) {
                const rootType = this.zMerge(x, y).getRootType(LZ);
                if (rootType === CFG_MULT) {
                    return false;
                } else if (rootType === CFG_INCR) {
                    ++incr;
                } else if (rootType === CFG_DECR) {
                    ++decr;
                }
            }
        }

        if (incr !== 0 && decr !== 0) {
            return false;
        }

        // Strongly mono, z-merge the boxes.
        r0.zStride *= 2;
        r1.valid = false;
        return true;
    }

    private addBox(x0: number, y0: number, z0: number, dx: number, dy: number, dz: number,
        LX: number, LY: number, LZ: number): void {
        const box = new OctBox(x0, y0, z0, dx, dy, dz, LX, LY, LZ);
        this.mBoxes.push(box);

        // Mark box edges in linear merge trees. This information will be
        // used later for extraction.
        this.xMerge(box.y0, box.z0).setEdge(box.LX);
        this.xMerge(box.y0, box.z1).setEdge(box.LX);
        this.xMerge(box.y1, box.z0).setEdge(box.LX);
        this.xMerge(box.y1, box.z1).setEdge(box.LX);
        this.yMerge(box.x0, box.z0).setEdge(box.LY);
        this.yMerge(box.x0, box.z1).setEdge(box.LY);
        this.yMerge(box.x1, box.z0).setEdge(box.LY);
        this.yMerge(box.x1, box.z1).setEdge(box.LY);
        this.zMerge(box.x0, box.y0).setEdge(box.LZ);
        this.zMerge(box.x0, box.y1).setEdge(box.LZ);
        this.zMerge(box.x1, box.y0).setEdge(box.LZ);
        this.zMerge(box.x1, box.y1).setEdge(box.LZ);
    }

    // ------------------------------------------------------------------
    // Support for tessellating monoboxes.
    // ------------------------------------------------------------------

    private tessellate(positions: Vertex[], triangles: TriangleKey[]): void {
        for (const box of this.mBoxes) {
            // Get vertices on edges of box.
            const table = new VETable();
            const type = this.getVertices(box, table);
            if (type === 0) {
                continue;
            }

            // Add wireframe edges to table, add face-vertices if necessary.
            if (box.dx > 1 || box.dy > 1 || box.dz > 1) {
                // Box is larger than voxel, each face has at most one edge.
                this.getXMinEdgesM(box, type, table);
                this.getXMaxEdgesM(box, type, table);
                this.getYMinEdgesM(box, type, table);
                this.getYMaxEdgesM(box, type, table);
                this.getZMinEdgesM(box, type, table);
                this.getZMaxEdgesM(box, type, table);

                if (table.getNumVertices() > I_QUANTITY) {
                    table.removeTrianglesSE(positions, triangles);
                } else {
                    table.removeTrianglesEC(positions, triangles);
                }
            } else {
                // The box is a 1x1x1 voxel. Do the full edge analysis but no
                // splitting is required.
                this.getXMinEdgesS(box, type, table);
                this.getXMaxEdgesS(box, type, table);
                this.getYMinEdgesS(box, type, table);
                this.getYMaxEdgesS(box, type, table);
                this.getZMinEdgesS(box, type, table);
                this.getZMaxEdgesS(box, type, table);
                table.removeTrianglesEC(positions, triangles);
            }
        }
    }

    private getXInterp(x: number, y: number, z: number): number {
        let index = x + this.mSize * (y + this.mSize * z);
        const f0 = this.mInputVoxels[index];
        ++index;
        const f1 = this.mInputVoxels[index];
        return x + (this.mLevel - f0) / (f1 - f0);
    }

    private getYInterp(x: number, y: number, z: number): number {
        let index = x + this.mSize * (y + this.mSize * z);
        const f0 = this.mInputVoxels[index];
        index += this.mSize;
        const f1 = this.mInputVoxels[index];
        return y + (this.mLevel - f0) / (f1 - f0);
    }

    private getZInterp(x: number, y: number, z: number): number {
        let index = x + this.mSize * (y + this.mSize * z);
        const f0 = this.mInputVoxels[index];
        index += this.mSizeSqr;
        const f1 = this.mInputVoxels[index];
        return z + (this.mLevel - f0) / (f1 - f0);
    }

    // Insert the interpolated vertex on each of the twelve box edges that
    // the level surface crosses; the return value is the bit flag of the
    // crossed edges.
    private getVertices(box: OctBox, table: VETable): number {
        let root: number;
        let type = 0;

        // xmin-ymin edge
        root = this.zMerge(box.x0, box.y0).getZeroBase(box.LZ);
        if (root !== -1) {
            type |= EB_XMIN_YMIN;
            table.insert(EI_XMIN_YMIN, box.x0, box.y0,
                this.getZInterp(box.x0, box.y0, root));
        }

        // xmin-ymax edge
        root = this.zMerge(box.x0, box.y1).getZeroBase(box.LZ);
        if (root !== -1) {
            type |= EB_XMIN_YMAX;
            table.insert(EI_XMIN_YMAX, box.x0, box.y1,
                this.getZInterp(box.x0, box.y1, root));
        }

        // xmax-ymin edge
        root = this.zMerge(box.x1, box.y0).getZeroBase(box.LZ);
        if (root !== -1) {
            type |= EB_XMAX_YMIN;
            table.insert(EI_XMAX_YMIN, box.x1, box.y0,
                this.getZInterp(box.x1, box.y0, root));
        }

        // xmax-ymax edge
        root = this.zMerge(box.x1, box.y1).getZeroBase(box.LZ);
        if (root !== -1) {
            type |= EB_XMAX_YMAX;
            table.insert(EI_XMAX_YMAX, box.x1, box.y1,
                this.getZInterp(box.x1, box.y1, root));
        }

        // xmin-zmin edge
        root = this.yMerge(box.x0, box.z0).getZeroBase(box.LY);
        if (root !== -1) {
            type |= EB_XMIN_ZMIN;
            table.insert(EI_XMIN_ZMIN, box.x0,
                this.getYInterp(box.x0, root, box.z0), box.z0);
        }

        // xmin-zmax edge
        root = this.yMerge(box.x0, box.z1).getZeroBase(box.LY);
        if (root !== -1) {
            type |= EB_XMIN_ZMAX;
            table.insert(EI_XMIN_ZMAX, box.x0,
                this.getYInterp(box.x0, root, box.z1), box.z1);
        }

        // xmax-zmin edge
        root = this.yMerge(box.x1, box.z0).getZeroBase(box.LY);
        if (root !== -1) {
            type |= EB_XMAX_ZMIN;
            table.insert(EI_XMAX_ZMIN, box.x1,
                this.getYInterp(box.x1, root, box.z0), box.z0);
        }

        // xmax-zmax edge
        root = this.yMerge(box.x1, box.z1).getZeroBase(box.LY);
        if (root !== -1) {
            type |= EB_XMAX_ZMAX;
            table.insert(EI_XMAX_ZMAX, box.x1,
                this.getYInterp(box.x1, root, box.z1), box.z1);
        }

        // ymin-zmin edge
        root = this.xMerge(box.y0, box.z0).getZeroBase(box.LX);
        if (root !== -1) {
            type |= EB_YMIN_ZMIN;
            table.insert(EI_YMIN_ZMIN,
                this.getXInterp(root, box.y0, box.z0), box.y0, box.z0);
        }

        // ymin-zmax edge
        root = this.xMerge(box.y0, box.z1).getZeroBase(box.LX);
        if (root !== -1) {
            type |= EB_YMIN_ZMAX;
            table.insert(EI_YMIN_ZMAX,
                this.getXInterp(root, box.y0, box.z1), box.y0, box.z1);
        }

        // ymax-zmin edge
        root = this.xMerge(box.y1, box.z0).getZeroBase(box.LX);
        if (root !== -1) {
            type |= EB_YMAX_ZMIN;
            table.insert(EI_YMAX_ZMIN,
                this.getXInterp(root, box.y1, box.z0), box.y1, box.z0);
        }

        // ymax-zmax edge
        root = this.xMerge(box.y1, box.z1).getZeroBase(box.LX);
        if (root !== -1) {
            type |= EB_YMAX_ZMAX;
            table.insert(EI_YMAX_ZMAX,
                this.getXInterp(root, box.y1, box.z1), box.y1, box.z1);
        }

        return type;
    }

    // The determinant f00 * f11 - f01 * f10 of the four face corner values,
    // computed with bigint to match the upstream int64_t arithmetic. The
    // sign selects the pairing of the four face intersections.
    private faceDeterminant(i00: number, i10: number, i11: number, i01: number): bigint {
        const f00 = BigInt(this.mInputVoxels[i00]);
        const f10 = BigInt(this.mInputVoxels[i10]);
        const f11 = BigInt(this.mInputVoxels[i11]);
        const f01 = BigInt(this.mInputVoxels[i01]);
        return f00 * f11 - f01 * f10;
    }

    // ------------------------------------------------------------------
    // Edge extraction for single boxes (1x1x1).
    // ------------------------------------------------------------------

    private getXMinEdgesS(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_YMIN) {
            faceType |= 0x01;
        }
        if (type & EB_XMIN_YMAX) {
            faceType |= 0x02;
        }
        if (type & EB_XMIN_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_XMIN_ZMAX) {
            faceType |= 0x08;
        }

        switch (faceType) {
            case 0:
                return;
            case 3:
                table.insertEdge(EI_XMIN_YMIN, EI_XMIN_YMAX);
                break;
            case 5:
                table.insertEdge(EI_XMIN_YMIN, EI_XMIN_ZMIN);
                break;
            case 6:
                table.insertEdge(EI_XMIN_YMAX, EI_XMIN_ZMIN);
                break;
            case 9:
                table.insertEdge(EI_XMIN_YMIN, EI_XMIN_ZMAX);
                break;
            case 10:
                table.insertEdge(EI_XMIN_YMAX, EI_XMIN_ZMAX);
                break;
            case 12:
                table.insertEdge(EI_XMIN_ZMIN, EI_XMIN_ZMAX);
                break;
            case 15: {
                // Four vertices, one per edge, need to disambiguate.
                const i = box.x0 + this.mSize * (box.y0 + this.mSize * box.z0);
                // F(x,y,z), F(x,y+1,z), F(x,y+1,z+1), F(x,y,z+1)
                const det = this.faceDeterminant(i, i + this.mSize,
                    i + this.mSize + this.mSizeSqr, i + this.mSizeSqr);

                if (det > 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_YMIN, EI_XMIN_ZMIN);
                    table.insertEdge(EI_XMIN_YMAX, EI_XMIN_ZMAX);
                } else if (det < 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_YMIN, EI_XMIN_ZMAX);
                    table.insertEdge(EI_XMIN_YMAX, EI_XMIN_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to
                    // tessellation.
                    table.insert(FI_XMIN,
                        table.getVertex(EI_XMIN_ZMIN)[0],
                        table.getVertex(EI_XMIN_ZMIN)[1],
                        table.getVertex(EI_XMIN_YMIN)[2]);

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_YMIN, FI_XMIN);
                    table.insertEdge(EI_XMIN_YMAX, FI_XMIN);
                    table.insertEdge(EI_XMIN_ZMIN, FI_XMIN);
                    table.insertEdge(EI_XMIN_ZMAX, FI_XMIN);
                }
                break;
            }
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }
    }

    private getXMaxEdgesS(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMAX_YMIN) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_YMAX) {
            faceType |= 0x02;
        }
        if (type & EB_XMAX_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_XMAX_ZMAX) {
            faceType |= 0x08;
        }

        switch (faceType) {
            case 0:
                return;
            case 3:
                table.insertEdge(EI_XMAX_YMIN, EI_XMAX_YMAX);
                break;
            case 5:
                table.insertEdge(EI_XMAX_YMIN, EI_XMAX_ZMIN);
                break;
            case 6:
                table.insertEdge(EI_XMAX_YMAX, EI_XMAX_ZMIN);
                break;
            case 9:
                table.insertEdge(EI_XMAX_YMIN, EI_XMAX_ZMAX);
                break;
            case 10:
                table.insertEdge(EI_XMAX_YMAX, EI_XMAX_ZMAX);
                break;
            case 12:
                table.insertEdge(EI_XMAX_ZMIN, EI_XMAX_ZMAX);
                break;
            case 15: {
                // Four vertices, one per edge, need to disambiguate.
                const i = box.x1 + this.mSize * (box.y0 + this.mSize * box.z0);
                // F(x,y,z), F(x,y+1,z), F(x,y+1,z+1), F(x,y,z+1)
                const det = this.faceDeterminant(i, i + this.mSize,
                    i + this.mSize + this.mSizeSqr, i + this.mSizeSqr);

                if (det > 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMAX_YMIN, EI_XMAX_ZMIN);
                    table.insertEdge(EI_XMAX_YMAX, EI_XMAX_ZMAX);
                } else if (det < 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMAX_YMIN, EI_XMAX_ZMAX);
                    table.insertEdge(EI_XMAX_YMAX, EI_XMAX_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to
                    // tessellation.
                    table.insert(FI_XMAX,
                        table.getVertex(EI_XMAX_ZMIN)[0],
                        table.getVertex(EI_XMAX_ZMIN)[1],
                        table.getVertex(EI_XMAX_YMIN)[2]);

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMAX_YMIN, FI_XMAX);
                    table.insertEdge(EI_XMAX_YMAX, FI_XMAX);
                    table.insertEdge(EI_XMAX_ZMIN, FI_XMAX);
                    table.insertEdge(EI_XMAX_ZMAX, FI_XMAX);
                }
                break;
            }
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }
    }

    private getYMinEdgesS(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_YMIN) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_YMIN) {
            faceType |= 0x02;
        }
        if (type & EB_YMIN_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_YMIN_ZMAX) {
            faceType |= 0x08;
        }

        switch (faceType) {
            case 0:
                return;
            case 3:
                table.insertEdge(EI_XMIN_YMIN, EI_XMAX_YMIN);
                break;
            case 5:
                table.insertEdge(EI_XMIN_YMIN, EI_YMIN_ZMIN);
                break;
            case 6:
                table.insertEdge(EI_XMAX_YMIN, EI_YMIN_ZMIN);
                break;
            case 9:
                table.insertEdge(EI_XMIN_YMIN, EI_YMIN_ZMAX);
                break;
            case 10:
                table.insertEdge(EI_XMAX_YMIN, EI_YMIN_ZMAX);
                break;
            case 12:
                table.insertEdge(EI_YMIN_ZMIN, EI_YMIN_ZMAX);
                break;
            case 15: {
                // Four vertices, one per edge, need to disambiguate.
                const i = box.x0 + this.mSize * (box.y0 + this.mSize * box.z0);
                // F(x,y,z), F(x+1,y,z), F(x+1,y,z+1), F(x,y,z+1)
                const det = this.faceDeterminant(i, i + 1,
                    i + 1 + this.mSizeSqr, i + this.mSizeSqr);

                if (det > 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_YMIN, EI_YMIN_ZMIN);
                    table.insertEdge(EI_XMAX_YMIN, EI_YMIN_ZMAX);
                } else if (det < 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_YMIN, EI_YMIN_ZMAX);
                    table.insertEdge(EI_XMAX_YMIN, EI_YMIN_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to
                    // tessellation.
                    table.insert(FI_YMIN,
                        table.getVertex(EI_YMIN_ZMIN)[0],
                        table.getVertex(EI_XMIN_YMIN)[1],
                        table.getVertex(EI_XMIN_YMIN)[2]);

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_YMIN, FI_YMIN);
                    table.insertEdge(EI_XMAX_YMIN, FI_YMIN);
                    table.insertEdge(EI_YMIN_ZMIN, FI_YMIN);
                    table.insertEdge(EI_YMIN_ZMAX, FI_YMIN);
                }
                break;
            }
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }
    }

    private getYMaxEdgesS(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_YMAX) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_YMAX) {
            faceType |= 0x02;
        }
        if (type & EB_YMAX_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_YMAX_ZMAX) {
            faceType |= 0x08;
        }

        switch (faceType) {
            case 0:
                return;
            case 3:
                table.insertEdge(EI_XMIN_YMAX, EI_XMAX_YMAX);
                break;
            case 5:
                table.insertEdge(EI_XMIN_YMAX, EI_YMAX_ZMIN);
                break;
            case 6:
                table.insertEdge(EI_XMAX_YMAX, EI_YMAX_ZMIN);
                break;
            case 9:
                table.insertEdge(EI_XMIN_YMAX, EI_YMAX_ZMAX);
                break;
            case 10:
                table.insertEdge(EI_XMAX_YMAX, EI_YMAX_ZMAX);
                break;
            case 12:
                table.insertEdge(EI_YMAX_ZMIN, EI_YMAX_ZMAX);
                break;
            case 15: {
                // Four vertices, one per edge, need to disambiguate.
                const i = box.x0 + this.mSize * (box.y1 + this.mSize * box.z0);
                // F(x,y,z), F(x+1,y,z), F(x+1,y,z+1), F(x,y,z+1)
                const det = this.faceDeterminant(i, i + 1,
                    i + 1 + this.mSizeSqr, i + this.mSizeSqr);

                if (det > 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_YMAX, EI_YMAX_ZMIN);
                    table.insertEdge(EI_XMAX_YMAX, EI_YMAX_ZMAX);
                } else if (det < 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_YMAX, EI_YMAX_ZMAX);
                    table.insertEdge(EI_XMAX_YMAX, EI_YMAX_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to
                    // tessellation.
                    table.insert(FI_YMAX,
                        table.getVertex(EI_YMAX_ZMIN)[0],
                        table.getVertex(EI_XMIN_YMAX)[1],
                        table.getVertex(EI_XMIN_YMAX)[2]);

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_YMAX, FI_YMAX);
                    table.insertEdge(EI_XMAX_YMAX, FI_YMAX);
                    table.insertEdge(EI_YMAX_ZMIN, FI_YMAX);
                    table.insertEdge(EI_YMAX_ZMAX, FI_YMAX);
                }
                break;
            }
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }
    }

    private getZMinEdgesS(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_ZMIN) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_ZMIN) {
            faceType |= 0x02;
        }
        if (type & EB_YMIN_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_YMAX_ZMIN) {
            faceType |= 0x08;
        }

        switch (faceType) {
            case 0:
                return;
            case 3:
                table.insertEdge(EI_XMIN_ZMIN, EI_XMAX_ZMIN);
                break;
            case 5:
                table.insertEdge(EI_XMIN_ZMIN, EI_YMIN_ZMIN);
                break;
            case 6:
                table.insertEdge(EI_XMAX_ZMIN, EI_YMIN_ZMIN);
                break;
            case 9:
                table.insertEdge(EI_XMIN_ZMIN, EI_YMAX_ZMIN);
                break;
            case 10:
                table.insertEdge(EI_XMAX_ZMIN, EI_YMAX_ZMIN);
                break;
            case 12:
                table.insertEdge(EI_YMIN_ZMIN, EI_YMAX_ZMIN);
                break;
            case 15: {
                // Four vertices, one per edge, need to disambiguate.
                const i = box.x0 + this.mSize * (box.y0 + this.mSize * box.z0);
                // F(x,y,z), F(x+1,y,z), F(x+1,y+1,z), F(x,y+1,z)
                const det = this.faceDeterminant(i, i + 1,
                    i + 1 + this.mSize, i + this.mSize);

                if (det > 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_ZMIN, EI_YMIN_ZMIN);
                    table.insertEdge(EI_XMAX_ZMIN, EI_YMAX_ZMIN);
                } else if (det < 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_ZMIN, EI_YMAX_ZMIN);
                    table.insertEdge(EI_XMAX_ZMIN, EI_YMIN_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to
                    // tessellation.
                    table.insert(FI_ZMIN,
                        table.getVertex(EI_YMIN_ZMIN)[0],
                        table.getVertex(EI_XMIN_ZMIN)[1],
                        table.getVertex(EI_XMIN_ZMIN)[2]);

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_ZMIN, FI_ZMIN);
                    table.insertEdge(EI_XMAX_ZMIN, FI_ZMIN);
                    table.insertEdge(EI_YMIN_ZMIN, FI_ZMIN);
                    table.insertEdge(EI_YMAX_ZMIN, FI_ZMIN);
                }
                break;
            }
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }
    }

    private getZMaxEdgesS(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_ZMAX) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_ZMAX) {
            faceType |= 0x02;
        }
        if (type & EB_YMIN_ZMAX) {
            faceType |= 0x04;
        }
        if (type & EB_YMAX_ZMAX) {
            faceType |= 0x08;
        }

        switch (faceType) {
            case 0:
                return;
            case 3:
                table.insertEdge(EI_XMIN_ZMAX, EI_XMAX_ZMAX);
                break;
            case 5:
                table.insertEdge(EI_XMIN_ZMAX, EI_YMIN_ZMAX);
                break;
            case 6:
                table.insertEdge(EI_XMAX_ZMAX, EI_YMIN_ZMAX);
                break;
            case 9:
                table.insertEdge(EI_XMIN_ZMAX, EI_YMAX_ZMAX);
                break;
            case 10:
                table.insertEdge(EI_XMAX_ZMAX, EI_YMAX_ZMAX);
                break;
            case 12:
                table.insertEdge(EI_YMIN_ZMAX, EI_YMAX_ZMAX);
                break;
            case 15: {
                // Four vertices, one per edge, need to disambiguate.
                const i = box.x0 + this.mSize * (box.y0 + this.mSize * box.z1);
                // F(x,y,z), F(x+1,y,z), F(x+1,y+1,z), F(x,y+1,z)
                const det = this.faceDeterminant(i, i + 1,
                    i + 1 + this.mSize, i + this.mSize);

                if (det > 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_ZMAX, EI_YMIN_ZMAX);
                    table.insertEdge(EI_XMAX_ZMAX, EI_YMAX_ZMAX);
                } else if (det < 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_ZMAX, EI_YMAX_ZMAX);
                    table.insertEdge(EI_XMAX_ZMAX, EI_YMIN_ZMAX);
                } else {
                    // Plus-sign configuration, add branch point to
                    // tessellation.
                    table.insert(FI_ZMAX,
                        table.getVertex(EI_YMIN_ZMAX)[0],
                        table.getVertex(EI_XMIN_ZMAX)[1],
                        table.getVertex(EI_XMIN_ZMAX)[2]);

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_ZMAX, FI_ZMAX);
                    table.insertEdge(EI_XMAX_ZMAX, FI_ZMAX);
                    table.insertEdge(EI_YMIN_ZMAX, FI_ZMAX);
                    table.insertEdge(EI_YMAX_ZMAX, FI_ZMAX);
                }
                break;
            }
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }
    }

    // ------------------------------------------------------------------
    // Edge extraction for merged boxes. Each face has at most one polyline,
    // whose endpoints are box-edge vertices and whose interior vertices come
    // from the grid lines interior to the face.
    // ------------------------------------------------------------------

    private getZMinEdgesM(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_ZMIN) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_ZMIN) {
            faceType |= 0x02;
        }
        if (type & EB_YMIN_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_YMAX_ZMIN) {
            faceType |= 0x08;
        }

        let end0 = 0;
        let end1 = 0;
        switch (faceType) {
            case 0:
                return;
            case 3:
                end0 = EI_XMIN_ZMIN;
                end1 = EI_XMAX_ZMIN;
                break;
            case 5:
                end0 = EI_XMIN_ZMIN;
                end1 = EI_YMIN_ZMIN;
                break;
            case 6:
                end0 = EI_YMIN_ZMIN;
                end1 = EI_XMAX_ZMIN;
                break;
            case 9:
                end0 = EI_XMIN_ZMIN;
                end1 = EI_YMAX_ZMIN;
                break;
            case 10:
                end0 = EI_YMAX_ZMIN;
                end1 = EI_XMAX_ZMIN;
                break;
            case 12:
                end0 = EI_YMIN_ZMIN;
                end1 = EI_YMAX_ZMIN;
                break;
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }

        const vSet = new VertexSet(sort0);

        for (let x = box.x0 + 1; x < box.x1; ++x) {
            const merge = this.yMerge(x, box.z0);
            if (merge.isZeroEdge(box.LY) || merge.hasZeroSubedge(box.LY)) {
                const root = merge.getZeroBase(box.LY);
                vSet.insert([x, this.getYInterp(x, root, box.z0), box.z0]);
            }
        }

        for (let y = box.y0 + 1; y < box.y1; ++y) {
            const merge = this.xMerge(y, box.z0);
            if (merge.isZeroEdge(box.LX) || merge.hasZeroSubedge(box.LX)) {
                const root = merge.getZeroBase(box.LX);
                vSet.insert([this.getXInterp(root, y, box.z0), y, box.z0]);
            }
        }

        insertSubdividedFaceEdges(table, end0, end1, vSet.sorted(), 0);
    }

    private getZMaxEdgesM(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_ZMAX) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_ZMAX) {
            faceType |= 0x02;
        }
        if (type & EB_YMIN_ZMAX) {
            faceType |= 0x04;
        }
        if (type & EB_YMAX_ZMAX) {
            faceType |= 0x08;
        }

        let end0 = 0;
        let end1 = 0;
        switch (faceType) {
            case 0:
                return;
            case 3:
                end0 = EI_XMIN_ZMAX;
                end1 = EI_XMAX_ZMAX;
                break;
            case 5:
                end0 = EI_XMIN_ZMAX;
                end1 = EI_YMIN_ZMAX;
                break;
            case 6:
                end0 = EI_YMIN_ZMAX;
                end1 = EI_XMAX_ZMAX;
                break;
            case 9:
                end0 = EI_XMIN_ZMAX;
                end1 = EI_YMAX_ZMAX;
                break;
            case 10:
                end0 = EI_YMAX_ZMAX;
                end1 = EI_XMAX_ZMAX;
                break;
            case 12:
                end0 = EI_YMIN_ZMAX;
                end1 = EI_YMAX_ZMAX;
                break;
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }

        const vSet = new VertexSet(sort0);

        for (let x = box.x0 + 1; x < box.x1; ++x) {
            const merge = this.yMerge(x, box.z1);
            if (merge.isZeroEdge(box.LY) || merge.hasZeroSubedge(box.LY)) {
                const root = merge.getZeroBase(box.LY);
                vSet.insert([x, this.getYInterp(x, root, box.z1), box.z1]);
            }
        }

        for (let y = box.y0 + 1; y < box.y1; ++y) {
            const merge = this.xMerge(y, box.z1);
            if (merge.isZeroEdge(box.LX) || merge.hasZeroSubedge(box.LX)) {
                const root = merge.getZeroBase(box.LX);
                vSet.insert([this.getXInterp(root, y, box.z1), y, box.z1]);
            }
        }

        insertSubdividedFaceEdges(table, end0, end1, vSet.sorted(), 0);
    }

    private getYMinEdgesM(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_YMIN) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_YMIN) {
            faceType |= 0x02;
        }
        if (type & EB_YMIN_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_YMIN_ZMAX) {
            faceType |= 0x08;
        }

        let end0 = 0;
        let end1 = 0;
        switch (faceType) {
            case 0:
                return;
            case 3:
                end0 = EI_XMIN_YMIN;
                end1 = EI_XMAX_YMIN;
                break;
            case 5:
                end0 = EI_XMIN_YMIN;
                end1 = EI_YMIN_ZMIN;
                break;
            case 6:
                end0 = EI_YMIN_ZMIN;
                end1 = EI_XMAX_YMIN;
                break;
            case 9:
                end0 = EI_XMIN_YMIN;
                end1 = EI_YMIN_ZMAX;
                break;
            case 10:
                end0 = EI_YMIN_ZMAX;
                end1 = EI_XMAX_YMIN;
                break;
            case 12:
                end0 = EI_YMIN_ZMIN;
                end1 = EI_YMIN_ZMAX;
                break;
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }

        const vSet = new VertexSet(sort1);

        for (let x = box.x0 + 1; x < box.x1; ++x) {
            const merge = this.zMerge(x, box.y0);
            if (merge.isZeroEdge(box.LZ) || merge.hasZeroSubedge(box.LZ)) {
                const root = merge.getZeroBase(box.LZ);
                vSet.insert([x, box.y0, this.getZInterp(x, box.y0, root)]);
            }
        }

        for (let z = box.z0 + 1; z < box.z1; ++z) {
            const merge = this.xMerge(box.y0, z);
            if (merge.isZeroEdge(box.LX) || merge.hasZeroSubedge(box.LX)) {
                const root = merge.getZeroBase(box.LX);
                vSet.insert([this.getXInterp(root, box.y0, z), box.y0, z]);
            }
        }

        insertSubdividedFaceEdges(table, end0, end1, vSet.sorted(), 2);
    }

    private getYMaxEdgesM(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_YMAX) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_YMAX) {
            faceType |= 0x02;
        }
        if (type & EB_YMAX_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_YMAX_ZMAX) {
            faceType |= 0x08;
        }

        let end0 = 0;
        let end1 = 0;
        switch (faceType) {
            case 0:
                return;
            case 3:
                end0 = EI_XMIN_YMAX;
                end1 = EI_XMAX_YMAX;
                break;
            case 5:
                end0 = EI_XMIN_YMAX;
                end1 = EI_YMAX_ZMIN;
                break;
            case 6:
                end0 = EI_YMAX_ZMIN;
                end1 = EI_XMAX_YMAX;
                break;
            case 9:
                end0 = EI_XMIN_YMAX;
                end1 = EI_YMAX_ZMAX;
                break;
            case 10:
                end0 = EI_YMAX_ZMAX;
                end1 = EI_XMAX_YMAX;
                break;
            case 12:
                end0 = EI_YMAX_ZMIN;
                end1 = EI_YMAX_ZMAX;
                break;
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }

        const vSet = new VertexSet(sort1);

        for (let x = box.x0 + 1; x < box.x1; ++x) {
            const merge = this.zMerge(x, box.y1);
            if (merge.isZeroEdge(box.LZ) || merge.hasZeroSubedge(box.LZ)) {
                const root = merge.getZeroBase(box.LZ);
                vSet.insert([x, box.y1, this.getZInterp(x, box.y1, root)]);
            }
        }

        for (let z = box.z0 + 1; z < box.z1; ++z) {
            const merge = this.xMerge(box.y1, z);
            if (merge.isZeroEdge(box.LX) || merge.hasZeroSubedge(box.LX)) {
                const root = merge.getZeroBase(box.LX);
                vSet.insert([this.getXInterp(root, box.y1, z), box.y1, z]);
            }
        }

        insertSubdividedFaceEdges(table, end0, end1, vSet.sorted(), 2);
    }

    private getXMinEdgesM(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMIN_YMIN) {
            faceType |= 0x01;
        }
        if (type & EB_XMIN_YMAX) {
            faceType |= 0x02;
        }
        if (type & EB_XMIN_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_XMIN_ZMAX) {
            faceType |= 0x08;
        }

        let end0 = 0;
        let end1 = 0;
        switch (faceType) {
            case 0:
                return;
            case 3:
                end0 = EI_XMIN_YMIN;
                end1 = EI_XMIN_YMAX;
                break;
            case 5:
                end0 = EI_XMIN_YMIN;
                end1 = EI_XMIN_ZMIN;
                break;
            case 6:
                end0 = EI_XMIN_ZMIN;
                end1 = EI_XMIN_YMAX;
                break;
            case 9:
                end0 = EI_XMIN_YMIN;
                end1 = EI_XMIN_ZMAX;
                break;
            case 10:
                end0 = EI_XMIN_ZMAX;
                end1 = EI_XMIN_YMAX;
                break;
            case 12:
                end0 = EI_XMIN_ZMIN;
                end1 = EI_XMIN_ZMAX;
                break;
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }

        const vSet = new VertexSet(sort2);

        for (let z = box.z0 + 1; z < box.z1; ++z) {
            const merge = this.yMerge(box.x0, z);
            if (merge.isZeroEdge(box.LY) || merge.hasZeroSubedge(box.LY)) {
                const root = merge.getZeroBase(box.LY);
                vSet.insert([box.x0, this.getYInterp(box.x0, root, z), z]);
            }
        }

        for (let y = box.y0 + 1; y < box.y1; ++y) {
            const merge = this.zMerge(box.x0, y);
            if (merge.isZeroEdge(box.LZ) || merge.hasZeroSubedge(box.LZ)) {
                const root = merge.getZeroBase(box.LZ);
                vSet.insert([box.x0, y, this.getZInterp(box.x0, y, root)]);
            }
        }

        insertSubdividedFaceEdges(table, end0, end1, vSet.sorted(), 1);
    }

    private getXMaxEdgesM(box: OctBox, type: number, table: VETable): void {
        let faceType = 0;
        if (type & EB_XMAX_YMIN) {
            faceType |= 0x01;
        }
        if (type & EB_XMAX_YMAX) {
            faceType |= 0x02;
        }
        if (type & EB_XMAX_ZMIN) {
            faceType |= 0x04;
        }
        if (type & EB_XMAX_ZMAX) {
            faceType |= 0x08;
        }

        let end0 = 0;
        let end1 = 0;
        switch (faceType) {
            case 0:
                return;
            case 3:
                end0 = EI_XMAX_YMIN;
                end1 = EI_XMAX_YMAX;
                break;
            case 5:
                end0 = EI_XMAX_YMIN;
                end1 = EI_XMAX_ZMIN;
                break;
            case 6:
                end0 = EI_XMAX_ZMIN;
                end1 = EI_XMAX_YMAX;
                break;
            case 9:
                end0 = EI_XMAX_YMIN;
                end1 = EI_XMAX_ZMAX;
                break;
            case 10:
                end0 = EI_XMAX_ZMAX;
                end1 = EI_XMAX_YMAX;
                break;
            case 12:
                end0 = EI_XMAX_ZMIN;
                end1 = EI_XMAX_ZMAX;
                break;
            default:
                logError(NON_INTEGER_LEVEL_MESSAGE);
        }

        const vSet = new VertexSet(sort2);

        for (let z = box.z0 + 1; z < box.z1; ++z) {
            const merge = this.yMerge(box.x1, z);
            if (merge.isZeroEdge(box.LY) || merge.hasZeroSubedge(box.LY)) {
                const root = merge.getZeroBase(box.LY);
                vSet.insert([box.x1, this.getYInterp(box.x1, root, z), z]);
            }
        }

        for (let y = box.y0 + 1; y < box.y1; ++y) {
            const merge = this.zMerge(box.x1, y);
            if (merge.isZeroEdge(box.LZ) || merge.hasZeroSubedge(box.LZ)) {
                const root = merge.getZeroBase(box.LZ);
                vSet.insert([box.x1, y, this.getZInterp(box.x1, y, root)]);
            }
        }

        insertSubdividedFaceEdges(table, end0, end1, vSet.sorted(), 1);
    }

    // ------------------------------------------------------------------
    // Support for normal vector calculations.
    // ------------------------------------------------------------------

    // The gradient of the trilinear interpolant of the image at the given
    // position. Positions outside the voxel grid produce the zero vector.
    private getGradient(position: Vertex): Vertex {
        const vzero: Vertex = [0, 0, 0];
        const x = Math.trunc(position[0]);
        if (x < 0 || x >= this.mTwoPowerN) {
            return vzero;
        }

        const y = Math.trunc(position[1]);
        if (y < 0 || y >= this.mTwoPowerN) {
            return vzero;
        }

        const z = Math.trunc(position[2]);
        if (z < 0 || z >= this.mTwoPowerN) {
            return vzero;
        }

        const i000 = x + this.mSize * (y + this.mSize * z);
        const i100 = i000 + 1;
        const i010 = i000 + this.mSize;
        const i110 = i100 + this.mSize;
        const i001 = i000 + this.mSizeSqr;
        const i101 = i100 + this.mSizeSqr;
        const i011 = i010 + this.mSizeSqr;
        const i111 = i110 + this.mSizeSqr;
        const f000 = this.mInputVoxels[i000];
        const f100 = this.mInputVoxels[i100];
        const f010 = this.mInputVoxels[i010];
        const f110 = this.mInputVoxels[i110];
        const f001 = this.mInputVoxels[i001];
        const f101 = this.mInputVoxels[i101];
        const f011 = this.mInputVoxels[i011];
        const f111 = this.mInputVoxels[i111];

        const fx = position[0] - x;
        const fy = position[1] - y;
        const fz = position[2] - z;
        const oneMinusX = 1 - fx;
        const oneMinusY = 1 - fy;
        const oneMinusZ = 1 - fz;
        const gradient: Vertex = [0, 0, 0];
        let tmp0: number;
        let tmp1: number;

        tmp0 = oneMinusY * (f100 - f000) + fy * (f110 - f010);
        tmp1 = oneMinusY * (f101 - f001) + fy * (f111 - f011);
        gradient[0] = oneMinusZ * tmp0 + fz * tmp1;

        tmp0 = oneMinusX * (f010 - f000) + fx * (f110 - f100);
        tmp1 = oneMinusX * (f011 - f001) + fx * (f111 - f101);
        gradient[1] = oneMinusZ * tmp0 + fz * tmp1;

        tmp0 = oneMinusX * (f001 - f000) + fx * (f101 - f100);
        tmp1 = oneMinusX * (f011 - f010) + fx * (f111 - f110);
        gradient[2] = oneMinusY * tmp0 + fy * tmp1;

        return gradient;
    }
}

// The canonical string key of a vertex, used to replicate the upstream
// std::map<Vertex, int32_t> duplicate detection.
function vertexKey(v: Vertex): string {
    return v[0] + ',' + v[1] + ',' + v[2];
}

// The cross product (v1 - v0) x (v2 - v0).
function triangleNormal(v0: Vertex, v1: Vertex, v2: Vertex): Vertex {
    const edge1: Vertex = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    const edge2: Vertex = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    return [
        edge1[1] * edge2[2] - edge1[2] * edge2[1],
        edge1[2] * edge2[0] - edge1[0] * edge2[2],
        edge1[0] * edge2[1] - edge1[1] * edge2[0]
    ];
}
