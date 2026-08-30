// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AdaptiveSkeletonClimbing2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Extract level curves using an adaptive approach to reduce the segment
// count. The implementation is for the algorithm described in the paper
//   Multiresolution Isosurface Extraction with Adaptive Skeleton Climbing
//   Tim Poston, Tien-Tsin Wong and Pheng-Ann Heng
//   Computer Graphics forum, volume 17, issue 3, September 1998
//   pages 137-147
// https://onlinelibrary.wiley.com/doi/abs/10.1111/1467-8659.00261
//
// Port notes: the upstream template types <T, Real> both become number. The
// image type T must contain integer values (upstream restricts T to the
// 8-, 16- and 32-bit integer types); the Real type is for extraction to
// floating-point vertices. The upstream int64_t determinant used to
// disambiguate the four-intersection configuration is computed with bigint
// so that products of 32-bit pixel values do not lose precision in IEEE
// doubles. Vertex (std::array<Real, 2>) and Edge (std::array<int32_t, 2>)
// become [number, number] tuples. The out-parameters of Extract become the
// return value { vertices, edges }; MakeUnique mutates the arrays it is
// given, as upstream does. The upstream MakeUnique packs vertices and edges
// by their std::map insertion indices, so first-encounter order is
// preserved; the port replicates this with (insertion-ordered) Map keyed by
// a canonical string of the tuple. The private debugging function
// PrintRectangles (ostream output) is not ported.

import { logAssert, logError } from './Logger';

// Configuration flags for LinearMergeTree nodes. CFG_MULT is the bitwise OR
// of CFG_INCR and CFG_DECR.
const CFG_NONE = 0;
const CFG_INCR = 1;
const CFG_DECR = 2;
const CFG_MULT = 3;

// Helper classes for the skeleton climbing (upstream private nested types).
class QuadRectangle {
    xOrigin: number;
    yOrigin: number;
    xStride: number;
    yStride: number;
    valid: boolean;

    constructor() {
        this.xOrigin = 0;
        this.yOrigin = 0;
        this.xStride = 0;
        this.yStride = 0;
        this.valid = false;
    }

    initialize(xOrigin: number, yOrigin: number, xStride: number, yStride: number): void {
        this.xOrigin = xOrigin;
        this.yOrigin = yOrigin;
        this.xStride = xStride;
        this.yStride = yStride;
        this.valid = true;
    }

    copyFrom(other: QuadRectangle): void {
        this.xOrigin = other.xOrigin;
        this.yOrigin = other.yOrigin;
        this.xStride = other.xStride;
        this.yStride = other.yStride;
        this.valid = other.valid;
    }
}

class QuadNode {
    r00: QuadRectangle;
    r10: QuadRectangle;
    r01: QuadRectangle;
    r11: QuadRectangle;

    constructor() {
        // The QuadRectangle members are default-constructed (invalid), as
        // upstream.
        this.r00 = new QuadRectangle();
        this.r10 = new QuadRectangle();
        this.r01 = new QuadRectangle();
        this.r11 = new QuadRectangle();
    }

    initialize(xOrigin: number, yOrigin: number, xNext: number, yNext: number, stride: number): void {
        this.r00.initialize(xOrigin, yOrigin, stride, stride);
        this.r10.initialize(xNext, yOrigin, stride, stride);
        this.r01.initialize(xOrigin, yNext, stride, stride);
        this.r11.initialize(xNext, yNext, stride, stride);
    }

    // Replicates the C++ value copy 'QuadNode node1 = node0'.
    copyFrom(other: QuadNode): void {
        this.r00.copyFrom(other.r00);
        this.r10.copyFrom(other.r10);
        this.r01.copyFrom(other.r01);
        this.r11.copyFrom(other.r11);
    }

    isMono(): boolean {
        return !this.r10.valid && !this.r01.valid && !this.r11.valid;
    }

    getQuantity(): number {
        let quantity = 0;

        if (this.r00.valid) {
            ++quantity;
        }

        if (this.r10.valid) {
            ++quantity;
        }

        if (this.r01.valid) {
            ++quantity;
        }

        if (this.r11.valid) {
            ++quantity;
        }

        return quantity;
    }
}

class LinearMergeTree {
    private mTwoPowerN: number;
    private mNodes: number[];

    constructor(N: number) {
        this.mTwoPowerN = 1 << N;
        this.mNodes = new Array<number>(2 * this.mTwoPowerN - 1).fill(0);
    }

    // Member access.
    getQuantity(): number {
        return 2 * this.mTwoPowerN - 1;
    }

    getNode(i: number): number {
        return this.mNodes[i];
    }

    getEdge(i: number): number {
        // assert: mNodes[i] === CFG_INCR || mNodes[i] === CFG_DECR

        // Traverse binary tree looking for incr or decr leaf node.
        const firstLeaf = this.mTwoPowerN - 1;
        while (i < firstLeaf) {
            i = 2 * i + 1;
            if (this.mNodes[i] === CFG_NONE) {
                ++i;
            }
        }

        return i - firstLeaf;
    }

    setLevel(level: number, data: ArrayLike<number>, offset: number, stride: number): void {
        // Assert: The 'level' is not an image value. Because the image
        // values are integers, choose 'level' to be a number that does not
        // represent an integer.

        // Determine the sign changes between pairs of consecutive samples.
        const firstLeaf = this.mTwoPowerN - 1;
        for (let i = 0, leaf = firstLeaf; i < this.mTwoPowerN; ++i, ++leaf) {
            const base = offset + stride * i;
            const value0 = data[base];
            const value1 = data[base + stride];

            if (value0 > level) {
                if (value1 > level) {
                    this.mNodes[leaf] = CFG_NONE;
                } else {
                    this.mNodes[leaf] = CFG_DECR;
                }
            } else { // value0 < level
                if (value1 > level) {
                    this.mNodes[leaf] = CFG_INCR;
                } else {
                    this.mNodes[leaf] = CFG_NONE;
                }
            }
        }

        // Propagate the sign change information up the binary tree.
        for (let i = firstLeaf - 1; i >= 0; --i) {
            const twoIp1 = 2 * i + 1;
            const child0 = this.mNodes[twoIp1];
            const child1 = this.mNodes[twoIp1 + 1];
            this.mNodes[i] = (child0 | child1);
        }
    }
}

class ASCRectangle {
    xOrigin: number;
    yOrigin: number;
    xStride: number;
    yStride: number;
    yOfXMin: number;
    yOfXMax: number;
    xOfYMin: number;
    xOfYMax: number;

    // A 4-bit flag for how the level set intersects the rectangle
    // boundary.
    //   bit 0 = xmin edge
    //   bit 1 = xmax edge
    //   bit 2 = ymin edge
    //   bit 3 = ymax edge
    // A bit is set if the corresponding edge is intersected by the level
    // set. This information is known from the CFG flags for
    // LinearMergeTree. Intersection occurs whenever the flag is CFG_INCR
    // or CFG_DECR.
    type: number;

    constructor(xOrigin: number, yOrigin: number, xStride: number, yStride: number) {
        this.xOrigin = xOrigin;
        this.yOrigin = yOrigin;
        this.xStride = xStride;
        this.yStride = yStride;
        this.yOfXMin = -1;
        this.yOfXMax = -1;
        this.xOfYMin = -1;
        this.xOfYMax = -1;
        this.type = 0;
    }
}

class AreaMergeTree {
    private mXMerge: LinearMergeTree[];
    private mYMerge: LinearMergeTree[];
    private mNodes: QuadNode[];

    constructor(N: number, xMerge: LinearMergeTree[], yMerge: LinearMergeTree[]) {
        this.mXMerge = xMerge;
        this.mYMerge = yMerge;
        const numNodes = ((1 << (2 * (N + 1))) - 1) / 3;
        this.mNodes = new Array<QuadNode>(numNodes);
        for (let i = 0; i < numNodes; ++i) {
            this.mNodes[i] = new QuadNode();
        }
    }

    constructMono(A: number, LX: number, LY: number, xOrigin: number, yOrigin: number,
        stride: number, depth: number): void {
        if (stride > 1) { // internal nodes
            const hStride = Math.trunc(stride / 2);

            let ABase = 4 * A;
            const A00 = ++ABase;
            const A10 = ++ABase;
            const A01 = ++ABase;
            const A11 = ++ABase;

            let LXBase = 2 * LX;
            const LX0 = ++LXBase;
            const LX1 = ++LXBase;

            let LYBase = 2 * LY;
            const LY0 = ++LYBase;
            const LY1 = ++LYBase;

            const xNext = xOrigin + hStride;
            const yNext = yOrigin + hStride;

            const depthM1 = depth - 1;
            this.constructMono(A00, LX0, LY0, xOrigin, yOrigin, hStride, depthM1);
            this.constructMono(A10, LX1, LY0, xNext, yOrigin, hStride, depthM1);
            this.constructMono(A01, LX0, LY1, xOrigin, yNext, hStride, depthM1);
            this.constructMono(A11, LX1, LY1, xNext, yNext, hStride, depthM1);

            if (depth >= 0) {
                // Merging is prevented above the specified depth in the
                // tree. This allows a single object to produce any
                // resolution isocontour rather than using multiple objects
                // to do so.
                this.mNodes[A].initialize(xOrigin, yOrigin, xNext, yNext, hStride);
                return;
            }

            const mono00 = this.mNodes[A00].isMono();
            const mono10 = this.mNodes[A10].isMono();
            const mono01 = this.mNodes[A01].isMono();
            const mono11 = this.mNodes[A11].isMono();

            const node0 = new QuadNode();
            node0.initialize(xOrigin, yOrigin, xNext, yNext, hStride);
            const node1 = new QuadNode();
            node1.copyFrom(node0);

            // Merge x first, y second.
            if (mono00 && mono10) {
                this.doXMerge(node0.r00, node0.r10, LX, yOrigin);
            }
            if (mono01 && mono11) {
                this.doXMerge(node0.r01, node0.r11, LX, yNext);
            }
            if (mono00 && mono01) {
                this.doYMerge(node0.r00, node0.r01, xOrigin, LY);
            }
            if (mono10 && mono11) {
                this.doYMerge(node0.r10, node0.r11, xNext, LY);
            }

            // Merge y first, x second.
            if (mono00 && mono01) {
                this.doYMerge(node1.r00, node1.r01, xOrigin, LY);
            }
            if (mono10 && mono11) {
                this.doYMerge(node1.r10, node1.r11, xNext, LY);
            }
            if (mono00 && mono10) {
                this.doXMerge(node1.r00, node1.r10, LX, yOrigin);
            }
            if (mono01 && mono11) {
                this.doXMerge(node1.r01, node1.r11, LX, yNext);
            }

            // Choose the merge that produced the smallest number of
            // rectangles.
            if (node0.getQuantity() <= node1.getQuantity()) {
                this.mNodes[A].copyFrom(node0);
            } else {
                this.mNodes[A].copyFrom(node1);
            }
        } else { // leaf nodes
            this.mNodes[A].r00.initialize(xOrigin, yOrigin, 1, 1);
        }
    }

    getRectangles(A: number, LX: number, LY: number, xOrigin: number, yOrigin: number,
        stride: number, rectangles: ASCRectangle[]): void {
        const hStride = Math.trunc(stride / 2);
        let ABase = 4 * A;
        const A00 = ++ABase;
        const A10 = ++ABase;
        const A01 = ++ABase;
        const A11 = ++ABase;
        let LXBase = 2 * LX;
        const LX0 = ++LXBase;
        const LX1 = ++LXBase;
        let LYBase = 2 * LY;
        const LY0 = ++LYBase;
        const LY1 = ++LYBase;
        const xNext = xOrigin + hStride;
        const yNext = yOrigin + hStride;

        const r00 = this.mNodes[A].r00;
        if (r00.valid) {
            if (r00.xStride === stride) {
                if (r00.yStride === stride) {
                    rectangles.push(this.getRectangle(r00, LX, LY));
                } else {
                    rectangles.push(this.getRectangle(r00, LX, LY0));
                }
            } else {
                if (r00.yStride === stride) {
                    rectangles.push(this.getRectangle(r00, LX0, LY));
                } else {
                    this.getRectangles(A00, LX0, LY0, xOrigin, yOrigin, hStride, rectangles);
                }
            }
        }

        const r10 = this.mNodes[A].r10;
        if (r10.valid) {
            if (r10.yStride === stride) {
                rectangles.push(this.getRectangle(r10, LX1, LY));
            } else {
                this.getRectangles(A10, LX1, LY0, xNext, yOrigin, hStride, rectangles);
            }
        }

        const r01 = this.mNodes[A].r01;
        if (r01.valid) {
            if (r01.xStride === stride) {
                rectangles.push(this.getRectangle(r01, LX, LY1));
            } else {
                this.getRectangles(A01, LX0, LY1, xOrigin, yNext, hStride, rectangles);
            }
        }

        const r11 = this.mNodes[A].r11;
        if (r11.valid) {
            this.getRectangles(A11, LX1, LY1, xNext, yNext, hStride, rectangles);
        }
    }

    private doXMerge(r0: QuadRectangle, r1: QuadRectangle, LX: number, yOrigin: number): void {
        if (r0.valid && r1.valid && r0.yStride === r1.yStride) {
            // Rectangles are x-mergeable.
            let incr = 0, decr = 0;
            for (let y = 0; y <= r0.yStride; ++y) {
                switch (this.mXMerge[yOrigin + y].getNode(LX)) {
                    case CFG_MULT:
                        return;
                    case CFG_INCR:
                        ++incr;
                        break;
                    case CFG_DECR:
                        ++decr;
                        break;
                }
            }

            if (incr === 0 || decr === 0) {
                // Strongly mono, x-merge the rectangles.
                r0.xStride *= 2;
                r1.valid = false;
            }
        }
    }

    private doYMerge(r0: QuadRectangle, r1: QuadRectangle, xOrigin: number, LY: number): void {
        if (r0.valid && r1.valid && r0.xStride === r1.xStride) {
            // Rectangles are y-mergeable.
            let incr = 0, decr = 0;
            for (let x = 0; x <= r0.xStride; ++x) {
                switch (this.mYMerge[xOrigin + x].getNode(LY)) {
                    case CFG_MULT:
                        return;
                    case CFG_INCR:
                        ++incr;
                        break;
                    case CFG_DECR:
                        ++decr;
                        break;
                }
            }

            if (incr === 0 || decr === 0) {
                // Strongly mono, y-merge the rectangles.
                r0.yStride *= 2;
                r1.valid = false;
            }
        }
    }

    private getRectangle(qrect: QuadRectangle, LX: number, LY: number): ASCRectangle {
        const rect = new ASCRectangle(qrect.xOrigin, qrect.yOrigin, qrect.xStride, qrect.yStride);

        // xmin edge
        let merge = this.mYMerge[qrect.xOrigin];
        if (merge.getNode(LY) !== CFG_NONE) {
            rect.yOfXMin = merge.getEdge(LY);
            if (rect.yOfXMin !== -1) {
                rect.type |= 0x01;
            }
        }

        // xmax edge
        merge = this.mYMerge[qrect.xOrigin + qrect.xStride];
        if (merge.getNode(LY) !== CFG_NONE) {
            rect.yOfXMax = merge.getEdge(LY);
            if (rect.yOfXMax !== -1) {
                rect.type |= 0x02;
            }
        }

        // ymin edge
        merge = this.mXMerge[qrect.yOrigin];
        if (merge.getNode(LX) !== CFG_NONE) {
            rect.xOfYMin = merge.getEdge(LX);
            if (rect.xOfYMin !== -1) {
                rect.type |= 0x04;
            }
        }

        // ymax edge
        merge = this.mXMerge[qrect.yOrigin + qrect.yStride];
        if (merge.getNode(LX) !== CFG_NONE) {
            rect.xOfYMax = merge.getEdge(LX);
            if (rect.xOfYMax !== -1) {
                rect.type |= 0x08;
            }
        }

        return rect;
    }
}

export class AdaptiveSkeletonClimbing2 {
    // Storage of image data.
    private mTwoPowerN: number;
    private mSize: number;
    private mInputPixels: ArrayLike<number>;

    // Trees for linear merging.
    private mXMerge: LinearMergeTree[];
    private mYMerge: LinearMergeTree[];

    // Tree for area merging.
    private mXYMerge: AreaMergeTree;

    // Construction. The input image is assumed to contain
    // (2^N+1)-by-(2^N+1) elements where N >= 1, with integer-valued
    // pixels. The organization is row-major order for (x,y).
    constructor(N: number, inputPixels: ArrayLike<number>) {
        this.mTwoPowerN = 1 << N;
        this.mSize = this.mTwoPowerN + 1;
        this.mInputPixels = inputPixels;

        if (N <= 0 || inputPixels == null) {
            logError('Invalid input.');
        }

        this.mXMerge = new Array<LinearMergeTree>(this.mSize);
        this.mYMerge = new Array<LinearMergeTree>(this.mSize);
        for (let i = 0; i < this.mSize; ++i) {
            this.mXMerge[i] = new LinearMergeTree(N);
            this.mYMerge[i] = new LinearMergeTree(N);
        }

        this.mXYMerge = new AreaMergeTree(N, this.mXMerge, this.mYMerge);
    }

    // Extract the level curve for the specified 'level', which must not be
    // an integer (image values are integers, so this guarantees the level
    // curve misses the sample points). The 'depth' controls the resolution:
    // merging is prevented above the specified depth in the tree when
    // depth >= 0; pass a negative depth for the fully adaptive (coarsest)
    // extraction. The vertices are (x, y) pairs and the edges are pairs of
    // indices into the vertices array.
    extract(level: number, depth: number):
        { vertices: [number, number][], edges: [number, number][] } {
        const rectangles: ASCRectangle[] = [];
        const localVertices: [number, number][] = [];
        const localEdges: [number, number][] = [];

        this.setLevel(level, depth);
        this.getRectangles(rectangles);
        for (const rectangle of rectangles) {
            if (rectangle.type > 0) {
                this.getComponents(level, rectangle, localVertices, localEdges);
            }
        }

        return { vertices: localVertices, edges: localEdges };
    }

    // Remove duplicate vertices and edges, remapping the edge indices. The
    // input arrays are modified in place, as upstream.
    makeUnique(vertices: [number, number][], edges: [number, number][]): void {
        const numVertices = vertices.length;
        const numEdges = edges.length;
        if (numVertices === 0 || numEdges === 0) {
            return;
        }

        // Compute the map of unique vertices and assign to them new and
        // unique indices.
        const vmap = new Map<string, { vertex: [number, number], index: number }>();
        let nextVertex = 0;
        for (let v = 0; v < numVertices; ++v) {
            // Keep only unique vertices.
            const key = vertices[v][0] + ',' + vertices[v][1];
            if (!vmap.has(key)) {
                vmap.set(key, { vertex: vertices[v], index: nextVertex });
                ++nextVertex;
            }
        }

        // Compute the map of unique edges and assign to them new and
        // unique indices.
        const emap = new Map<string, { edge: [number, number], index: number }>();
        let nextEdge = 0;
        for (let e = 0; e < numEdges; ++e) {
            // Replace old vertex indices by new vertex indices.
            const edge = edges[e];
            for (let i = 0; i < 2; ++i) {
                const vertex = vertices[edge[i]];
                const element = vmap.get(vertex[0] + ',' + vertex[1]);
                logAssert(element !== undefined, 'Expecting the vertex to be in the vmap.');
                edge[i] = element.index;
            }

            // Keep only unique edges.
            const key = edge[0] + ',' + edge[1];
            if (!emap.has(key)) {
                emap.set(key, { edge, index: nextEdge });
                ++nextEdge;
            }
        }

        // Pack the vertices into an array.
        vertices.length = vmap.size;
        for (const element of vmap.values()) {
            vertices[element.index] = element.vertex;
        }

        // Pack the edges into an array.
        edges.length = emap.size;
        for (const element of emap.values()) {
            edges[element.index] = element.edge;
        }
    }

    // Support for extraction of level sets.
    private getInterp(level: number, base: number, index: number, increment: number): number {
        const f0 = this.mInputPixels[index];
        index += increment;
        const f1 = this.mInputPixels[index];
        logAssert((f0 - level) * (f1 - level) < 0, 'Unexpected condition.');
        return base + (level - f0) / (f1 - f0);
    }

    private addEdge(vertices: [number, number][], edges: [number, number][],
        x0: number, y0: number, x1: number, y1: number): void {
        const v0 = vertices.length;
        const v1 = v0 + 1;
        edges.push([v0, v1]);
        vertices.push([x0, y0]);
        vertices.push([x1, y1]);
    }

    private setLevel(level: number, depth: number): void {
        for (let y = 0; y < this.mSize; ++y) {
            const offset = this.mSize * y;
            const stride = 1;
            this.mXMerge[y].setLevel(level, this.mInputPixels, offset, stride);
        }

        for (let x = 0; x < this.mSize; ++x) {
            const offset = x;
            const stride = this.mSize;
            this.mYMerge[x].setLevel(level, this.mInputPixels, offset, stride);
        }

        this.mXYMerge.constructMono(0, 0, 0, 0, 0, this.mTwoPowerN, depth);
    }

    private getRectangles(rectangles: ASCRectangle[]): void {
        this.mXYMerge.getRectangles(0, 0, 0, 0, 0, this.mTwoPowerN, rectangles);
    }

    private getComponents(level: number, rectangle: ASCRectangle,
        vertices: [number, number][], edges: [number, number][]): void {
        let x: number, y: number;
        let x0: number, y0: number, x1: number, y1: number;

        switch (rectangle.type) {
            case 3: // two vertices, on xmin and xmax
                logAssert(rectangle.yOfXMin !== -1, 'Unexpected condition.');
                x = rectangle.xOrigin;
                y = rectangle.yOfXMin;
                x0 = x;
                y0 = this.getInterp(level, y, x + this.mSize * y, this.mSize);

                logAssert(rectangle.yOfXMax !== -1, 'Unexpected condition.');
                x = rectangle.xOrigin + rectangle.xStride;
                y = rectangle.yOfXMax;
                x1 = x;
                y1 = this.getInterp(level, y, x + this.mSize * y, this.mSize);

                this.addEdge(vertices, edges, x0, y0, x1, y1);
                break;
            case 5: // two vertices, on xmin and ymin
                logAssert(rectangle.yOfXMin !== -1, 'Unexpected condition.');
                x = rectangle.xOrigin;
                y = rectangle.yOfXMin;
                x0 = x;
                y0 = this.getInterp(level, y, x + this.mSize * y, this.mSize);

                logAssert(rectangle.xOfYMin !== -1, 'Unexpected condition.');
                x = rectangle.xOfYMin;
                y = rectangle.yOrigin;
                x1 = this.getInterp(level, x, x + this.mSize * y, 1);
                y1 = y;

                this.addEdge(vertices, edges, x0, y0, x1, y1);
                break;
            case 6: // two vertices, on xmax and ymin
                logAssert(rectangle.yOfXMax !== -1, 'Unexpected condition.');
                x = rectangle.xOrigin + rectangle.xStride;
                y = rectangle.yOfXMax;
                x0 = x;
                y0 = this.getInterp(level, y, x + this.mSize * y, this.mSize);

                logAssert(rectangle.xOfYMin !== -1, 'Unexpected condition.');
                x = rectangle.xOfYMin;
                y = rectangle.yOrigin;
                x1 = this.getInterp(level, x, x + this.mSize * y, 1);
                y1 = y;

                this.addEdge(vertices, edges, x0, y0, x1, y1);
                break;
            case 9: // two vertices, on xmin and ymax
                logAssert(rectangle.yOfXMin !== -1, 'Unexpected condition.');
                x = rectangle.xOrigin;
                y = rectangle.yOfXMin;
                x0 = x;
                y0 = this.getInterp(level, y, x + this.mSize * y, this.mSize);

                logAssert(rectangle.xOfYMax !== -1, 'Unexpected condition.');
                x = rectangle.xOfYMax;
                y = rectangle.yOrigin + rectangle.yStride;
                x1 = this.getInterp(level, x, x + this.mSize * y, 1);
                y1 = y;

                this.addEdge(vertices, edges, x0, y0, x1, y1);
                break;
            case 10: // two vertices, on xmax and ymax
                logAssert(rectangle.yOfXMax !== -1, 'Unexpected condition.');
                x = rectangle.xOrigin + rectangle.xStride;
                y = rectangle.yOfXMax;
                x0 = x;
                y0 = this.getInterp(level, y, x + this.mSize * y, this.mSize);

                logAssert(rectangle.xOfYMax !== -1, 'Unexpected condition.');
                x = rectangle.xOfYMax;
                y = rectangle.yOrigin + rectangle.yStride;
                x1 = this.getInterp(level, x, x + this.mSize * y, 1);
                y1 = y;

                this.addEdge(vertices, edges, x0, y0, x1, y1);
                break;
            case 12: // two vertices, on ymin and ymax
                logAssert(rectangle.xOfYMin !== -1, 'Unexpected condition.');
                x = rectangle.xOfYMin;
                y = rectangle.yOrigin;
                x0 = this.getInterp(level, x, x + this.mSize * y, 1);
                y0 = y;

                logAssert(rectangle.xOfYMax !== -1, 'Unexpected condition.');
                x = rectangle.xOfYMax;
                y = rectangle.yOrigin + rectangle.yStride;
                x1 = this.getInterp(level, x, x + this.mSize * y, 1);
                y1 = y;

                this.addEdge(vertices, edges, x0, y0, x1, y1);
                break;
            case 15: { // four vertices, one per edge, need to disambiguate
                logAssert(rectangle.xStride === 1 && rectangle.yStride === 1,
                    'Unexpected condition.');

                logAssert(rectangle.yOfXMin !== -1, 'Unexpected condition.');
                x = rectangle.xOrigin;
                y = rectangle.yOfXMin;
                x0 = x;
                y0 = this.getInterp(level, y, x + this.mSize * y, this.mSize);

                logAssert(rectangle.yOfXMax !== -1, 'Unexpected condition.');
                x = rectangle.xOrigin + rectangle.xStride;
                y = rectangle.yOfXMax;
                x1 = x;
                y1 = this.getInterp(level, y, x + this.mSize * y, this.mSize);

                logAssert(rectangle.xOfYMin !== -1, 'Unexpected condition.');
                x = rectangle.xOfYMin;
                y = rectangle.yOrigin;
                const fx2 = this.getInterp(level, x, x + this.mSize * y, 1);
                const fy2 = y;

                logAssert(rectangle.xOfYMax !== -1, 'Unexpected condition.');
                x = rectangle.xOfYMax;
                y = rectangle.yOrigin + rectangle.yStride;
                const fx3 = this.getInterp(level, x, x + this.mSize * y, 1);
                const fy3 = y;

                // The determinant of integer pixel values is computed with
                // bigint, matching the upstream int64_t arithmetic.
                let index = rectangle.xOrigin + this.mSize * rectangle.yOrigin;
                const i00 = BigInt(this.mInputPixels[index]);
                ++index;
                const i10 = BigInt(this.mInputPixels[index]);
                index += this.mSize;
                const i11 = BigInt(this.mInputPixels[index]);
                --index;
                const i01 = BigInt(this.mInputPixels[index]);

                const det = i00 * i11 - i01 * i10;
                if (det > 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P2> and
                    // <P1,P3>.
                    this.addEdge(vertices, edges, x0, y0, fx2, fy2);
                    this.addEdge(vertices, edges, x1, y1, fx3, fy3);
                } else if (det < 0n) {
                    // Disjoint hyperbolic segments, pair <P0,P3> and
                    // <P1,P2>.
                    this.addEdge(vertices, edges, x0, y0, fx3, fy3);
                    this.addEdge(vertices, edges, x1, y1, fx2, fy2);
                } else {
                    // Plus-sign configuration, add branch point to
                    // tessellation.
                    const fx4 = fx2, fy4 = y0;
                    this.addEdge(vertices, edges, x0, y0, fx4, fy4);
                    this.addEdge(vertices, edges, x1, y1, fx4, fy4);
                    this.addEdge(vertices, edges, fx2, fy2, fx4, fy4);
                    this.addEdge(vertices, edges, fx3, fy3, fx4, fy4);
                }
                break;
            }
            default:
                logError('Unexpected condition.');
        }
    }
}
