// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CurveExtractor.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Abstract base class for extracting level curves from a 2D image. The
// derived classes (CurveExtractorSquares, CurveExtractorTriangles) provide
// the extraction itself; this class provides the rational vertex and edge
// representations, the removal of duplicate vertices and edges, and the
// conversion of rational vertices to floating-point vertices.
//
// Upstream the image type T must be one of the integer types int8_t,
// int16_t, int32_t, uint8_t, uint16_t or uint32_t, and internal integer
// computations are performed using int64_t. The type Real is for extraction
// to floating-point vertices.
//
// Port notes:
//   - The C++ static_assert that constrains T has no TypeScript equivalent;
//     the pixels are numbers that are expected to hold integer values.
//   - The int64_t vertex components become numbers holding integer values.
//     The rational comparisons multiply two components, so a product may
//     exceed 2^53 and lose exactness; compareProducts() below falls back to
//     bigint arithmetic in that case, which reproduces the exact int64_t
//     comparisons of upstream as long as the components themselves are exact
//     integers.
//   - Upstream has two Extract member functions distinguished by their
//     output parameter types. The port names the abstract rational one
//     extract() and the floating-point one extractReal(). Both return the
//     vertices and edges instead of writing to reference parameters, except
//     makeUnique() which modifies the arrays it is given, as upstream does.
//   - The nested structs Vertex and Edge become the top-level classes
//     CurveExtractorVertex and CurveExtractorEdge because exported names must
//     be unique across the library.

import { logAssert } from './Logger.js';

// Compare the products a*b and c*d exactly, returning -1, 0 or +1. The
// double products are used when they are exact; otherwise the comparison is
// redone with bigint. If the true product has magnitude of at least 2^53, so
// does the rounded double product, which makes the test conservative.
function compareProducts(a: number, b: number, c: number, d: number): number {
    const p = a * b;
    const q = c * d;
    if (Math.abs(p) <= Number.MAX_SAFE_INTEGER && Math.abs(q) <= Number.MAX_SAFE_INTEGER) {
        return p < q ? -1 : (p > q ? 1 : 0);
    }
    const bp = BigInt(a) * BigInt(b);
    const bq = BigInt(c) * BigInt(d);
    return bp < bq ? -1 : (bp > bq ? 1 : 0);
}

function gcd(a: number, b: number): number {
    let x = Math.abs(a);
    let y = Math.abs(b);
    while (y !== 0) {
        const t = x % y;
        x = y;
        y = t;
    }
    return x;
}

// A canonical string for the rational pair, used to replicate the equivalence
// relation of the upstream std::map<Vertex, int32_t> keys. Two vertices have
// the same key if and only if neither compares less than the other.
function rationalKey(vertex: CurveExtractorVertex): string {
    const gx = gcd(vertex.xNumer, vertex.xDenom) || 1;
    const gy = gcd(vertex.yNumer, vertex.yDenom) || 1;
    return `${vertex.xNumer / gx}/${vertex.xDenom / gx},${vertex.yNumer / gy}/${vertex.yDenom / gy}`;
}

// The level curves form a graph of vertices and edges. The vertices are
// computed as pairs of nonnegative rational numbers. CurveExtractorVertex
// represents the rational pair (xNumer/xDenom, yNumer/yDenom) as
// (xNumer, xDenom, yNumer, yDenom), where all components are nonnegative.
// The edges connect pairs of vertices, forming a graph that represents the
// level set.
export class CurveExtractorVertex {
    xNumer: number;
    xDenom: number;
    yNumer: number;
    yDenom: number;

    // The vertex generation leads to the numerator and denominator having
    // the same sign. The constructor changes sign to ensure the numerator and
    // denominator are both positive. (The subtractions from zero rather than
    // unary negations keep a zero component as +0.)
    constructor(inXNumer: number = 0, inXDenom: number = 0,
        inYNumer: number = 0, inYDenom: number = 0) {
        if (inXDenom > 0) {
            this.xNumer = inXNumer;
            this.xDenom = inXDenom;
        } else {
            this.xNumer = 0 - inXNumer;
            this.xDenom = 0 - inXDenom;
        }

        if (inYDenom > 0) {
            this.yNumer = inYNumer;
            this.yDenom = inYDenom;
        } else {
            this.yNumer = 0 - inYNumer;
            this.yDenom = 0 - inYDenom;
        }
    }

    // The 4-argument constructor guarantees that xDenom > 0 and yDenom > 0.
    // The following comparisons assume that the denominators are positive.
    equals(other: CurveExtractorVertex): boolean {
        return (
            // xn0 / xd0 == xn1 / xd1
            compareProducts(this.xNumer, other.xDenom, other.xNumer, this.xDenom) === 0
            &&
            // yn0 / yd0 == yn1 / yd1
            compareProducts(this.yNumer, other.yDenom, other.yNumer, this.yDenom) === 0
        );
    }

    lessThan(other: CurveExtractorVertex): boolean {
        const compareX = compareProducts(this.xNumer, other.xDenom, other.xNumer, this.xDenom);
        if (compareX < 0) {
            // xn0/xd0 < xn1/xd1
            return true;
        }
        if (compareX > 0) {
            // xn0/xd0 > xn1/xd1
            return false;
        }

        // yn0/yd0 < yn1/yd1
        return compareProducts(this.yNumer, other.yDenom, other.yNumer, this.yDenom) < 0;
    }
}

// An edge is an unordered pair of vertex indices stored in increasing order.
export class CurveExtractorEdge {
    v: [number, number];

    constructor(v0: number = 0, v1: number = 0) {
        if (v0 < v1) {
            this.v = [v0, v1];
        } else {
            this.v = [v1, v0];
        }
    }

    equals(other: CurveExtractorEdge): boolean {
        return this.v[0] === other.v[0] && this.v[1] === other.v[1];
    }

    lessThan(other: CurveExtractorEdge): boolean {
        for (let i = 0; i < 2; ++i) {
            if (this.v[i] < other.v[i]) {
                return true;
            }
            if (this.v[i] > other.v[i]) {
                return false;
            }
        }
        return false;
    }
}

export abstract class CurveExtractor {
    protected mXBound: number;
    protected mYBound: number;
    protected mInputPixels: ArrayLike<number>;
    protected mPixels: number[];

    // The input is a 2D image with lexicographically ordered pixels (x,y)
    // stored in a linear array. Pixel (x,y) is stored in the array at
    // location index = x + xBound * y. The inputs xBound and yBound must each
    // be 2 or larger so that there is at least one image square to process.
    // The inputPixels must contain at least xBound * yBound elements.
    protected constructor(xBound: number, yBound: number, inputPixels: ArrayLike<number>) {
        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mInputPixels = inputPixels;
        logAssert(xBound > 1 && yBound > 1 && inputPixels !== null &&
            inputPixels.length >= xBound * yBound, 'Invalid input.');
        this.mPixels = new Array<number>(xBound * yBound).fill(0);
    }

    // Extract level curves and return rational vertices.
    abstract extract(level: number): {
        vertices: CurveExtractorVertex[];
        edges: CurveExtractorEdge[];
    };

    // Extract level curves and return floating-point vertices.
    extractReal(level: number, removeDuplicateVertices: boolean): {
        vertices: [number, number][];
        edges: CurveExtractorEdge[];
    } {
        const rational = this.extract(level);
        const edges = rational.edges;
        if (removeDuplicateVertices) {
            this.makeUnique(rational.vertices, edges);
        }
        return { vertices: this.convert(rational.vertices), edges };
    }

    // The extraction has duplicate vertices on edges shared by pixels. This
    // function will eliminate the duplicates. The input arrays are modified
    // in place, as the upstream reference parameters are.
    makeUnique(vertices: CurveExtractorVertex[], edges: CurveExtractorEdge[]): void {
        const numVertices = vertices.length;
        const numEdges = edges.length;
        if (numVertices === 0 || numEdges === 0) {
            return;
        }

        // Compute the map of unique vertices and assign to them new and
        // unique indices.
        const vmap = new Map<string, { vertex: CurveExtractorVertex; index: number }>();
        let nextVertex = 0;
        for (let v = 0; v < numVertices; ++v) {
            // Keep only unique vertices.
            const key = rationalKey(vertices[v]);
            if (!vmap.has(key)) {
                vmap.set(key, { vertex: vertices[v], index: nextVertex });
                ++nextVertex;
            }
        }

        // Compute the map of unique edges and assign to them new and unique
        // indices.
        const emap = new Map<string, { edge: CurveExtractorEdge; index: number }>();
        let nextEdge = 0;
        for (let e = 0; e < numEdges; ++e) {
            // Replace old vertex indices by new vertex indices.
            const edge = edges[e];
            for (let i = 0; i < 2; ++i) {
                const element = vmap.get(rationalKey(vertices[edge.v[i]]));
                logAssert(element !== undefined, 'Expecting the vertex to be in the vmap.');
                edge.v[i] = element.index;
            }

            // Keep only unique edges. NOTE: upstream does not reorder the
            // vertex indices of the edge after the replacement, so an edge
            // whose new indices are decreasing is not recognized as a
            // duplicate of the same undirected edge with increasing indices.
            // The port preserves that behavior.
            const key = `${edge.v[0]},${edge.v[1]}`;
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

    // Convert from CurveExtractorVertex to floating-point pairs.
    convert(input: CurveExtractorVertex[]): [number, number][] {
        const output = new Array<[number, number]>(input.length);
        for (let i = 0; i < input.length; ++i) {
            const rxNumer = input[i].xNumer;
            const rxDenom = input[i].xDenom;
            const ryNumer = input[i].yNumer;
            const ryDenom = input[i].yDenom;
            output[i] = [rxNumer / rxDenom, ryNumer / ryDenom];
        }
        return output;
    }

    protected addVertex(vertices: CurveExtractorVertex[],
        xNumer: number, xDenom: number, yNumer: number, yDenom: number): void {
        vertices.push(new CurveExtractorVertex(xNumer, xDenom, yNumer, yDenom));
    }

    protected addEdge(vertices: CurveExtractorVertex[], edges: CurveExtractorEdge[],
        xNumer0: number, xDenom0: number, yNumer0: number, yDenom0: number,
        xNumer1: number, xDenom1: number, yNumer1: number, yDenom1: number): void {
        const v0 = vertices.length;
        const v1 = v0 + 1;
        edges.push(new CurveExtractorEdge(v0, v1));
        vertices.push(new CurveExtractorVertex(xNumer0, xDenom0, yNumer0, yDenom0));
        vertices.push(new CurveExtractorVertex(xNumer1, xDenom1, yNumer1, yDenom1));
    }
}
