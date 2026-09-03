// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SurfaceExtractor.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Abstract base class for level-surface extraction from a 3D image of
// integer-valued voxels. The image type T must be an integer type in
// upstream (int8_t through uint32_t) with internal computations in int64_t.
//
// Port notes: T and Real are both number. Voxel values are expected to be
// integers; the rational vertex arithmetic (products of numerators and
// denominators) is exact as long as the products stay within 2^53, the
// safe-integer range of IEEE double (upstream uses int64_t, which has the
// analogous constraint at 2^63). The nested structs Vertex and Triangle are
// exported as SurfaceExtractorVertex and SurfaceExtractorTriangle because
// src/index.ts star-exports every file. The C++ comparison operators become
// the methods equals() and lessThan(). The pure virtual
// Extract(level, vertices&, triangles&) is the abstract method
// extractRational(level) returning the arrays, and the non-virtual overload
// Extract(level, removeDuplicateVertices, vertices&, triangles&) is
// extract(level, removeDuplicateVertices). The upstream std::map in
// MakeUnique is used only for duplicate detection (new indices are assigned
// in insertion order, not map order), so the port uses a Map keyed by a
// canonical string; the output is identical.

import { logAssert } from './Logger.js';

// A vertex of the extracted level surface, stored as the rational triple
// (xNumer/xDenom, yNumer/yDenom, zNumer/zDenom) with all components
// nonnegative and denominators positive.
export class SurfaceExtractorVertex {
    xNumer: number;
    xDenom: number;
    yNumer: number;
    yDenom: number;
    zNumer: number;
    zDenom: number;

    constructor(inXNumer: number = 0, inXDenom: number = 1,
        inYNumer: number = 0, inYDenom: number = 1,
        inZNumer: number = 0, inZDenom: number = 1) {
        // The vertex generation leads to the numerator and denominator
        // having the same sign. This constructor changes sign to ensure the
        // numerator and denominator are both positive.
        if (inXDenom > 0) {
            this.xNumer = inXNumer;
            this.xDenom = inXDenom;
        } else {
            this.xNumer = -inXNumer;
            this.xDenom = -inXDenom;
        }

        if (inYDenom > 0) {
            this.yNumer = inYNumer;
            this.yDenom = inYDenom;
        } else {
            this.yNumer = -inYNumer;
            this.yDenom = -inYDenom;
        }

        if (inZDenom > 0) {
            this.zNumer = inZNumer;
            this.zDenom = inZDenom;
        } else {
            this.zNumer = -inZNumer;
            this.zDenom = -inZDenom;
        }
    }

    // The constructor guarantees that xDenom > 0, yDenom > 0 and zDenom > 0.
    // The comparisons assume the denominators are positive.
    equals(other: SurfaceExtractorVertex): boolean {
        // xn0/xd0 == xn1/xd1
        return this.xNumer * other.xDenom === other.xNumer * this.xDenom
            // yn0/yd0 == yn1/yd1
            && this.yNumer * other.yDenom === other.yNumer * this.yDenom
            // zn0/zd0 == zn1/zd1
            && this.zNumer * other.zDenom === other.zNumer * this.zDenom;
    }

    lessThan(other: SurfaceExtractorVertex): boolean {
        const xn0txd1 = this.xNumer * other.xDenom;
        const xn1txd0 = other.xNumer * this.xDenom;
        if (xn0txd1 < xn1txd0) {
            // xn0/xd0 < xn1/xd1
            return true;
        }
        if (xn0txd1 > xn1txd0) {
            // xn0/xd0 > xn1/xd1
            return false;
        }

        const yn0tyd1 = this.yNumer * other.yDenom;
        const yn1tyd0 = other.yNumer * this.yDenom;
        if (yn0tyd1 < yn1tyd0) {
            // yn0/yd0 < yn1/yd1
            return true;
        }
        if (yn0tyd1 > yn1tyd0) {
            // yn0/yd0 > yn1/yd1
            return false;
        }

        const zn0tzd1 = this.zNumer * other.zDenom;
        const zn1tzd0 = other.zNumer * this.zDenom;
        // zn0/zd0 < zn1/zd1
        return zn0tzd1 < zn1tzd0;
    }
}

// A triangle of the extracted level surface, stored as indices into the
// vertex array. The constructor stores the cyclic permutation of (v0,v1,v2)
// with the minimum index first, preserving the winding order.
export class SurfaceExtractorTriangle {
    v: [number, number, number];

    constructor(v0: number = 0, v1: number = 0, v2: number = 0) {
        // After the code is executed, (v[0],v[1],v[2]) is a cyclic
        // permutation of (v0,v1,v2) with v[0] = min{v0,v1,v2}.
        if (v0 < v1) {
            if (v0 < v2) {
                this.v = [v0, v1, v2];
            } else {
                this.v = [v2, v0, v1];
            }
        } else {
            if (v1 < v2) {
                this.v = [v1, v2, v0];
            } else {
                this.v = [v2, v0, v1];
            }
        }
    }

    equals(other: SurfaceExtractorTriangle): boolean {
        return this.v[0] === other.v[0] && this.v[1] === other.v[1]
            && this.v[2] === other.v[2];
    }

    lessThan(other: SurfaceExtractorTriangle): boolean {
        for (let i = 0; i < 3; ++i) {
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

export abstract class SurfaceExtractor {
    protected mXBound: number;
    protected mYBound: number;
    protected mZBound: number;
    protected mXYBound: number;
    protected mInputVoxels: ArrayLike<number>;
    protected mVoxels: number[];

    // The input is a 3D image with lexicographically ordered voxels (x,y,z)
    // stored in a linear array. Voxel (x,y,z) is stored in the array at
    // location index = x + xBound * (y + yBound * z). The inputs xBound,
    // yBound and zBound must each be 2 or larger so that there is at least
    // one image cube to process. The inputVoxels must contain at least
    // xBound * yBound * zBound elements.
    protected constructor(xBound: number, yBound: number, zBound: number,
        inputVoxels: ArrayLike<number>) {
        logAssert(xBound > 1 && yBound > 1 && zBound > 1
            && inputVoxels !== null && inputVoxels !== undefined,
            'Invalid input.');
        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mZBound = zBound;
        this.mXYBound = xBound * yBound;
        this.mInputVoxels = inputVoxels;
        this.mVoxels = new Array<number>(xBound * yBound * zBound).fill(0);
    }

    // Extract level surfaces and return rational vertices. The level
    // surfaces form a graph of vertices, edges and triangles. Derived
    // classes implement this (upstream pure virtual
    // Extract(level, vertices&, triangles&)).
    abstract extractRational(level: number): {
        vertices: SurfaceExtractorVertex[];
        triangles: SurfaceExtractorTriangle[];
    };

    // Extract level surfaces, converting the rational vertices to
    // floating-point triples [x, y, z].
    extract(level: number, removeDuplicateVertices: boolean): {
        vertices: [number, number, number][];
        triangles: SurfaceExtractorTriangle[];
    } {
        const { vertices: rationalVertices, triangles } = this.extractRational(level);
        if (removeDuplicateVertices) {
            this.makeUnique(rationalVertices, triangles);
        }
        return { vertices: this.convert(rationalVertices), triangles };
    }

    // The extraction has duplicate vertices on edges shared by voxels. This
    // function eliminates the duplicates, modifying the arrays in place.
    makeUnique(vertices: SurfaceExtractorVertex[],
        triangles: SurfaceExtractorTriangle[]): void {
        const numVertices = vertices.length;
        const numTriangles = triangles.length;
        if (numVertices === 0 || numTriangles === 0) {
            return;
        }

        // Compute the map of unique vertices and assign to them new and
        // unique indices.
        const vmap = new Map<string, number>();
        const uniqueVertices: SurfaceExtractorVertex[] = [];
        const vertexKeys: string[] = new Array<string>(numVertices);
        let nextVertex = 0;
        for (let v = 0; v < numVertices; ++v) {
            // Keep only unique vertices.
            const key = SurfaceExtractor.vertexKey(vertices[v]);
            vertexKeys[v] = key;
            if (!vmap.has(key)) {
                vmap.set(key, nextVertex);
                uniqueVertices.push(vertices[v]);
                ++nextVertex;
            }
        }

        // Compute the map of unique triangles and assign to them new and
        // unique indices.
        const tmap = new Map<string, number>();
        const uniqueTriangles: SurfaceExtractorTriangle[] = [];
        let nextTriangle = 0;
        for (let t = 0; t < numTriangles; ++t) {
            const triangle = triangles[t];
            for (let i = 0; i < 3; ++i) {
                const newIndex = vmap.get(vertexKeys[triangle.v[i]]);
                logAssert(newIndex !== undefined,
                    'Expecting the vertex to be in the vmap.');
                triangle.v[i] = newIndex;
            }

            // Keep only unique triangles. The Triangle constructor
            // recomputes the canonical cyclic permutation for the remapped
            // indices, matching the upstream comparison semantics.
            const canonical = new SurfaceExtractorTriangle(
                triangle.v[0], triangle.v[1], triangle.v[2]);
            const key = `${canonical.v[0]},${canonical.v[1]},${canonical.v[2]}`;
            if (!tmap.has(key)) {
                tmap.set(key, nextTriangle);
                uniqueTriangles.push(canonical);
                ++nextTriangle;
            }
        }

        // Pack the vertices into the caller's array.
        vertices.length = uniqueVertices.length;
        for (let i = 0; i < uniqueVertices.length; ++i) {
            vertices[i] = uniqueVertices[i];
        }

        // Pack the triangles into the caller's array.
        triangles.length = uniqueTriangles.length;
        for (let i = 0; i < uniqueTriangles.length; ++i) {
            triangles[i] = uniqueTriangles[i];
        }
    }

    // Convert from rational vertices to floating-point triples.
    convert(input: SurfaceExtractorVertex[]): [number, number, number][] {
        const output: [number, number, number][] = new Array(input.length);
        for (let i = 0; i < input.length; ++i) {
            output[i] = [
                input[i].xNumer / input[i].xDenom,
                input[i].yNumer / input[i].yDenom,
                input[i].zNumer / input[i].zDenom
            ];
        }
        return output;
    }

    // The extraction does not use any topological information about the
    // level surfaces. The triangles can be a mixture of clockwise-ordered
    // and counterclockwise-ordered. This function is an attempt to give the
    // triangles a consistent ordering by selecting a normal in approximately
    // the same direction as the average gradient at the vertices (when
    // sameDir is true), or in the opposite direction (when sameDir is
    // false). This might not always produce a consistent order, but is fast.
    orientTriangles(vertices: [number, number, number][],
        triangles: SurfaceExtractorTriangle[], sameDir: boolean): void {
        for (const triangle of triangles) {
            // Get the triangle vertices.
            const v0 = vertices[triangle.v[0]];
            const v1 = vertices[triangle.v[1]];
            const v2 = vertices[triangle.v[2]];

            // Construct the triangle normal based on the current
            // orientation.
            const edge1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            const edge2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            const normal = [
                edge1[1] * edge2[2] - edge1[2] * edge2[1],
                edge1[2] * edge2[0] - edge1[0] * edge2[2],
                edge1[0] * edge2[1] - edge1[1] * edge2[0]
            ];

            // Get the image gradient at the vertices.
            const grad0 = this.getGradient(v0);
            const grad1 = this.getGradient(v1);
            const grad2 = this.getGradient(v2);

            // Compute the average gradient.
            const gradAvr = [
                (grad0[0] + grad1[0] + grad2[0]) / 3,
                (grad0[1] + grad1[1] + grad2[1]) / 3,
                (grad0[2] + grad1[2] + grad2[2]) / 3
            ];

            // Compute the dot product of normal and average gradient.
            const dot = gradAvr[0] * normal[0] + gradAvr[1] * normal[1]
                + gradAvr[2] * normal[2];

            // Choose triangle orientation based on gradient direction.
            if (sameDir) {
                if (dot < 0) {
                    // Wrong orientation, reorder it.
                    const save = triangle.v[1];
                    triangle.v[1] = triangle.v[2];
                    triangle.v[2] = save;
                }
            } else {
                if (dot > 0) {
                    // Wrong orientation, reorder it.
                    const save = triangle.v[1];
                    triangle.v[1] = triangle.v[2];
                    triangle.v[2] = save;
                }
            }
        }
    }

    // Use this function if you want vertex normals for dynamic lighting of
    // the mesh. A vertex normal is the normalized area-weighted sum of the
    // normals to the triangles that share that vertex.
    computeNormals(vertices: [number, number, number][],
        triangles: SurfaceExtractorTriangle[]): [number, number, number][] {
        const normals: [number, number, number][] = new Array(vertices.length);
        for (let i = 0; i < vertices.length; ++i) {
            normals[i] = [0, 0, 0];
        }

        for (const triangle of triangles) {
            // Get the triangle vertices.
            const v0 = vertices[triangle.v[0]];
            const v1 = vertices[triangle.v[1]];
            const v2 = vertices[triangle.v[2]];

            // Construct the triangle normal.
            const edge1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            const edge2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            const normal = [
                edge1[1] * edge2[2] - edge1[2] * edge2[1],
                edge1[2] * edge2[0] - edge1[0] * edge2[2],
                edge1[0] * edge2[1] - edge1[1] * edge2[0]
            ];

            // Maintain the sum of normals at each vertex.
            for (let i = 0; i < 3; ++i) {
                for (let j = 0; j < 3; ++j) {
                    normals[triangle.v[i]][j] += normal[j];
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

    protected abstract getGradient(pos: [number, number, number]): [number, number, number];

    // A canonical string key for a rational vertex: each fraction is reduced
    // to lowest terms (denominators are positive by construction), so two
    // vertices produce the same key exactly when the upstream operator==
    // reports them equal.
    private static vertexKey(vertex: SurfaceExtractorVertex): string {
        const gx = SurfaceExtractor.gcd(vertex.xNumer, vertex.xDenom);
        const gy = SurfaceExtractor.gcd(vertex.yNumer, vertex.yDenom);
        const gz = SurfaceExtractor.gcd(vertex.zNumer, vertex.zDenom);
        return `${vertex.xNumer / gx}/${vertex.xDenom / gx},`
            + `${vertex.yNumer / gy}/${vertex.yDenom / gy},`
            + `${vertex.zNumer / gz}/${vertex.zDenom / gz}`;
    }

    private static gcd(a: number, b: number): number {
        a = Math.abs(a);
        b = Math.abs(b);
        while (b !== 0) {
            const r = a % b;
            a = b;
            b = r;
        }
        return (a !== 0 ? a : 1);
    }
}
