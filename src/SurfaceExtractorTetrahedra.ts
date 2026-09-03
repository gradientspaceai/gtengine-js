// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SurfaceExtractorTetrahedra.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The level set extraction algorithm implemented here is described in the
// document
// https://www.geometrictools.com/Documentation/ExtractLevelSurfaces.pdf
//
// Each image cube is partitioned into 5 tetrahedra, with the partition
// alternating by the parity of (x+y+z) so that adjacent cubes agree on the
// shared faces. The level surface within a tetrahedron is computed by linear
// interpolation, so the extraction is free of the topological ambiguities of
// the cube-based extractors.
//
// Port notes: T and Real are both number; the voxel values are expected to be
// integers, and the internal int64_t computations map onto IEEE doubles (the
// products stay exact while they remain within the safe-integer range). The
// upstream pure virtual Extract(level, vertices&, triangles&) is the base
// class method extractRational(level), which returns the two arrays. The
// nested struct Edge is exported as SurfaceExtractorTetrahedraEdge. The
// upstream std::map<Vertex,int32_t> and std::set<Edge>/std::set<Triangle> are
// replaced by Maps keyed by canonical strings; the vertex indices are
// assigned in insertion order (as upstream), and the edges and triangles are
// sorted at the end by the upstream operator< so the output ordering matches.
// The edge set is computed but, as upstream, not reported by the extraction;
// getEdges() exposes it for inspection.

import {
    SurfaceExtractor,
    SurfaceExtractorTriangle,
    SurfaceExtractorVertex
} from './SurfaceExtractor.js';

// An undirected edge of the extracted level surface, stored as the sorted
// pair of vertex indices.
export class SurfaceExtractorTetrahedraEdge {
    v: [number, number];

    constructor(v0: number = 0, v1: number = 0) {
        if (v0 < v1) {
            this.v = [v0, v1];
        } else {
            this.v = [v1, v0];
        }
    }

    equals(other: SurfaceExtractorTetrahedraEdge): boolean {
        return this.v[0] === other.v[0] && this.v[1] === other.v[1];
    }

    lessThan(other: SurfaceExtractorTetrahedraEdge): boolean {
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

// The greatest common divisor of |a| and |b|, at least 1.
function gcd(a: number, b: number): number {
    a = Math.abs(a);
    b = Math.abs(b);
    while (b !== 0) {
        const r = a % b;
        a = b;
        b = r;
    }
    return (a !== 0 ? a : 1);
}

// A canonical string key for a rational vertex. The SurfaceExtractorVertex
// constructor guarantees positive denominators, so reducing each fraction to
// lowest terms produces the same key exactly when the upstream Vertex
// operator== reports the vertices equal (equivalently, when neither is less
// than the other under the upstream operator<).
function vertexKey(vertex: SurfaceExtractorVertex): string {
    const gx = gcd(vertex.xNumer, vertex.xDenom);
    const gy = gcd(vertex.yNumer, vertex.yDenom);
    const gz = gcd(vertex.zNumer, vertex.zDenom);
    return `${vertex.xNumer / gx}/${vertex.xDenom / gx},`
        + `${vertex.yNumer / gy}/${vertex.yDenom / gy},`
        + `${vertex.zNumer / gz}/${vertex.zDenom / gz}`;
}

export class SurfaceExtractorTetrahedra extends SurfaceExtractor {
    // The unique vertices in insertion order, together with the map from a
    // canonical vertex key to the assigned index.
    protected mVMap: Map<string, number>;
    protected mVList: SurfaceExtractorVertex[];
    protected mESet: Map<string, SurfaceExtractorTetrahedraEdge>;
    protected mTSet: Map<string, SurfaceExtractorTriangle>;
    protected mNextVertex: number;

    // The input is a 3D image with lexicographically ordered voxels (x,y,z)
    // stored in a linear array. Voxel (x,y,z) is stored in the array at
    // location index = x + xBound * (y + yBound * z). The inputs xBound,
    // yBound and zBound must each be 2 or larger so that there is at least
    // one image cube to process. The inputVoxels must contain at least
    // xBound * yBound * zBound elements.
    constructor(xBound: number, yBound: number, zBound: number,
        inputVoxels: ArrayLike<number>) {
        super(xBound, yBound, zBound, inputVoxels);
        this.mVMap = new Map<string, number>();
        this.mVList = [];
        this.mESet = new Map<string, SurfaceExtractorTetrahedraEdge>();
        this.mTSet = new Map<string, SurfaceExtractorTriangle>();
        this.mNextVertex = 0;
    }

    // Extract level surfaces and return rational vertices. Use the
    // base-class extract() if you want real-valued vertices.
    override extractRational(level: number): {
        vertices: SurfaceExtractorVertex[];
        triangles: SurfaceExtractorTriangle[];
    } {
        // Adjust the image so that the level set is F(x,y,z) = 0.
        for (let i = 0; i < this.mVoxels.length; ++i) {
            this.mVoxels[i] = this.mInputVoxels[i] - level;
        }

        this.mVMap.clear();
        this.mVList.length = 0;
        this.mESet.clear();
        this.mTSet.clear();
        this.mNextVertex = 0;

        for (let z = 0, zp = 1; zp < this.mZBound; ++z, ++zp) {
            const zParity = (z & 1);
            for (let y = 0, yp = 1; yp < this.mYBound; ++y, ++yp) {
                const yParity = (y & 1);
                for (let x = 0, xp = 1; xp < this.mXBound; ++x, ++xp) {
                    const xParity = (x & 1);

                    const i000 = x + this.mXBound * (y + this.mYBound * z);
                    const i100 = i000 + 1;
                    const i010 = i000 + this.mXBound;
                    const i110 = i010 + 1;
                    const i001 = i000 + this.mXYBound;
                    const i101 = i001 + 1;
                    const i011 = i001 + this.mXBound;
                    const i111 = i011 + 1;
                    const f000 = this.mVoxels[i000];
                    const f100 = this.mVoxels[i100];
                    const f010 = this.mVoxels[i010];
                    const f110 = this.mVoxels[i110];
                    const f001 = this.mVoxels[i001];
                    const f101 = this.mVoxels[i101];
                    const f011 = this.mVoxels[i011];
                    const f111 = this.mVoxels[i111];

                    if (xParity ^ yParity ^ zParity) {
                        // 1205
                        this.processTetrahedron(
                            xp, y, z, f100,
                            xp, yp, z, f110,
                            x, y, z, f000,
                            xp, y, zp, f101);

                        // 3027
                        this.processTetrahedron(
                            x, yp, z, f010,
                            x, y, z, f000,
                            xp, yp, z, f110,
                            x, yp, zp, f011);

                        // 4750
                        this.processTetrahedron(
                            x, y, zp, f001,
                            x, yp, zp, f011,
                            xp, y, zp, f101,
                            x, y, z, f000);

                        // 6572
                        this.processTetrahedron(
                            xp, yp, zp, f111,
                            xp, y, zp, f101,
                            x, yp, zp, f011,
                            xp, yp, z, f110);

                        // 0752
                        this.processTetrahedron(
                            x, y, z, f000,
                            x, yp, zp, f011,
                            xp, y, zp, f101,
                            xp, yp, z, f110);
                    } else {
                        // 0134
                        this.processTetrahedron(
                            x, y, z, f000,
                            xp, y, z, f100,
                            x, yp, z, f010,
                            x, y, zp, f001);

                        // 2316
                        this.processTetrahedron(
                            xp, yp, z, f110,
                            x, yp, z, f010,
                            xp, y, z, f100,
                            xp, yp, zp, f111);

                        // 5461
                        this.processTetrahedron(
                            xp, y, zp, f101,
                            x, y, zp, f001,
                            xp, yp, zp, f111,
                            xp, y, z, f100);

                        // 7643
                        this.processTetrahedron(
                            x, yp, zp, f011,
                            xp, yp, zp, f111,
                            x, y, zp, f001,
                            x, yp, z, f010);

                        // 6314
                        this.processTetrahedron(
                            xp, yp, zp, f111,
                            x, yp, z, f010,
                            xp, y, z, f100,
                            x, y, zp, f001);
                    }
                }
            }
        }

        // Pack vertices into an array; the index assigned by addVertex is the
        // position in the array.
        const vertices = this.mVList.slice();

        // Pack triangles into an array, ordered by the upstream std::set
        // comparison.
        const triangles = Array.from(this.mTSet.values());
        triangles.sort((t0, t1) => (t0.lessThan(t1) ? -1 : (t1.lessThan(t0) ? 1 : 0)));

        return { vertices, triangles };
    }

    // The edges of the most recent extraction, ordered by the upstream
    // std::set comparison. Upstream computes the edge array in Extract but
    // does not report it to the caller.
    getEdges(): SurfaceExtractorTetrahedraEdge[] {
        const edges = Array.from(this.mESet.values());
        edges.sort((e0, e1) => (e0.lessThan(e1) ? -1 : (e1.lessThan(e0) ? 1 : 0)));
        return edges;
    }

    protected override getGradient(pos: [number, number, number]):
        [number, number, number] {
        const x = Math.trunc(pos[0]);
        if (x < 0 || x + 1 >= this.mXBound) {
            return [0, 0, 0];
        }

        const y = Math.trunc(pos[1]);
        if (y < 0 || y + 1 >= this.mYBound) {
            return [0, 0, 0];
        }

        const z = Math.trunc(pos[2]);
        if (z < 0 || z + 1 >= this.mZBound) {
            return [0, 0, 0];
        }

        // Get image values at corners of voxel.
        const i000 = x + this.mXBound * (y + this.mYBound * z);
        const i100 = i000 + 1;
        const i010 = i000 + this.mXBound;
        const i110 = i010 + 1;
        const i001 = i000 + this.mXYBound;
        const i101 = i001 + 1;
        const i011 = i001 + this.mXBound;
        const i111 = i011 + 1;
        const f000 = this.mVoxels[i000];
        const f100 = this.mVoxels[i100];
        const f010 = this.mVoxels[i010];
        const f110 = this.mVoxels[i110];
        const f001 = this.mVoxels[i001];
        const f101 = this.mVoxels[i101];
        const f011 = this.mVoxels[i011];
        const f111 = this.mVoxels[i111];

        const dx = pos[0] - x;
        const dy = pos[1] - y;
        const dz = pos[2] - z;

        const grad: [number, number, number] = [0, 0, 0];

        if ((x & 1) ^ (y & 1) ^ (z & 1)) {
            // The cube is partitioned into the corner tetrahedra at the
            // corners (1,0,0), (0,1,0), (0,0,1) and (1,1,1) plus the central
            // tetrahedron formed by the remaining corners.
            if (dx - dy - dz >= 0) {
                // 1205
                grad[0] = +f100 - f000;
                grad[1] = -f100 + f110;
                grad[2] = -f100 + f101;
            } else if (dx - dy + dz <= 0) {
                // 3027
                grad[0] = -f010 + f110;
                grad[1] = +f010 - f000;
                grad[2] = -f010 + f011;
            } else if (dx + dy - dz <= 0) {
                // 4750
                grad[0] = -f001 + f101;
                grad[1] = -f001 + f011;
                grad[2] = +f001 - f000;
            } else if (dx + dy + dz >= 2) {
                // Upstream bug: the test is 'dx + dy + dz >= 0', which is
                // always true for a point in the voxel, so the corner
                // tetrahedron at (1,1,1) absorbed the central tetrahedron and
                // the '0752' branch below was unreachable. The plane through
                // (1,1,0), (1,0,1) and (0,1,1) that cuts off the corner
                // (1,1,1) is dx + dy + dz = 2, matching the even-parity
                // branch below where the analogous cutting planes are at 1.
                // 6572
                grad[0] = +f111 - f011;
                grad[1] = +f111 - f101;
                grad[2] = +f111 - f110;
            } else {
                // 0752
                grad[0] = 0.5 * (-f000 - f011 + f101 + f110);
                grad[1] = 0.5 * (-f000 + f011 - f101 + f110);
                grad[2] = 0.5 * (-f000 + f011 + f101 - f110);
            }
        } else {
            // The cube is partitioned into the corner tetrahedra at the
            // corners (0,0,0), (1,1,0), (1,0,1) and (0,1,1) plus the central
            // tetrahedron formed by the remaining corners.
            if (dx + dy + dz <= 1) {
                // 0134
                grad[0] = -f000 + f100;
                grad[1] = -f000 + f010;
                grad[2] = -f000 + f001;
            } else if (dx + dy - dz >= 1) {
                // 2316
                grad[0] = +f110 - f010;
                grad[1] = +f110 - f100;
                grad[2] = -f110 + f111;
            } else if (dx - dy + dz >= 1) {
                // 5461
                grad[0] = +f101 - f001;
                grad[1] = -f101 + f111;
                grad[2] = +f101 - f100;
            } else if (-dx + dy + dz >= 1) {
                // 7643
                grad[0] = -f011 + f111;
                grad[1] = +f011 - f001;
                grad[2] = +f011 - f010;
            } else {
                // 6314
                grad[0] = 0.5 * (f111 - f010 + f100 - f001);
                grad[1] = 0.5 * (f111 + f010 - f100 - f001);
                grad[2] = 0.5 * (f111 - f010 - f100 + f001);
            }
        }

        return grad;
    }

    protected addVertex(v: SurfaceExtractorVertex): number {
        const key = vertexKey(v);
        const existing = this.mVMap.get(key);
        if (existing !== undefined) {
            // Vertex already in map, just return its unique index.
            return existing;
        }
        // Vertex not in map, insert it and assign it a unique index.
        const i = this.mNextVertex++;
        this.mVMap.set(key, i);
        this.mVList.push(v);
        return i;
    }

    protected addEdge(v0: SurfaceExtractorVertex,
        v1: SurfaceExtractorVertex): void {
        const i0 = this.addVertex(v0);
        const i1 = this.addVertex(v1);
        this.insertEdge(i0, i1);
    }

    protected addTriangle(v0: SurfaceExtractorVertex,
        v1: SurfaceExtractorVertex, v2: SurfaceExtractorVertex): void {
        const i0 = this.addVertex(v0);
        const i1 = this.addVertex(v1);
        const i2 = this.addVertex(v2);

        // Nothing to do if triangle already exists.
        const triangle = new SurfaceExtractorTriangle(i0, i1, i2);
        if (this.mTSet.has(SurfaceExtractorTetrahedra.triangleKey(triangle))) {
            return;
        }

        // Prevent double-sided triangles. Swapping the last two vertices of
        // the canonical form yields the canonical form of the reversed
        // triangle, because the minimum index remains first.
        const reversed = new SurfaceExtractorTriangle(
            triangle.v[0], triangle.v[2], triangle.v[1]);
        if (this.mTSet.has(SurfaceExtractorTetrahedra.triangleKey(reversed))) {
            return;
        }

        // At this time the triangle is not double-sided in the mesh.
        this.insertEdge(i0, i1);
        this.insertEdge(i1, i2);
        this.insertEdge(i2, i0);
        this.mTSet.set(SurfaceExtractorTetrahedra.triangleKey(triangle), triangle);
    }

    private insertEdge(i0: number, i1: number): void {
        const edge = new SurfaceExtractorTetrahedraEdge(i0, i1);
        const key = `${edge.v[0]},${edge.v[1]}`;
        if (!this.mESet.has(key)) {
            this.mESet.set(key, edge);
        }
    }

    private static triangleKey(triangle: SurfaceExtractorTriangle): string {
        return `${triangle.v[0]},${triangle.v[1]},${triangle.v[2]}`;
    }

    // Support for extraction with linear interpolation. The tetrahedron has
    // vertices (xi,yi,zi) with function values fi for 0 <= i <= 3. The 81
    // sign patterns of (f0,f1,f2,f3) are handled by reducing to the patterns
    // whose first nonzero value is positive.
    protected processTetrahedron(
        x0: number, y0: number, z0: number, f0: number,
        x1: number, y1: number, z1: number, f1: number,
        x2: number, y2: number, z2: number, f2: number,
        x3: number, y3: number, z3: number, f3: number): void {
        let xn0 = 0, yn0 = 0, zn0 = 0, d0 = 0;
        let xn1 = 0, yn1 = 0, zn1 = 0, d1 = 0;
        let xn2 = 0, yn2 = 0, zn2 = 0, d2 = 0;
        let xn3 = 0, yn3 = 0, zn3 = 0, d3 = 0;

        if (f0 !== 0) {
            // convert to case +***
            if (f0 < 0) {
                f0 = -f0;
                f1 = -f1;
                f2 = -f2;
                f3 = -f3;
            }

            if (f1 > 0) {
                if (f2 > 0) {
                    if (f3 > 0) {
                        // ++++
                        return;
                    } else if (f3 < 0) {
                        // +++-
                        d0 = f0 - f3;
                        xn0 = f0 * x3 - f3 * x0;
                        yn0 = f0 * y3 - f3 * y0;
                        zn0 = f0 * z3 - f3 * z0;
                        d1 = f1 - f3;
                        xn1 = f1 * x3 - f3 * x1;
                        yn1 = f1 * y3 - f3 * y1;
                        zn1 = f1 * z3 - f3 * z1;
                        d2 = f2 - f3;
                        xn2 = f2 * x3 - f3 * x2;
                        yn2 = f2 * y3 - f3 * y2;
                        zn2 = f2 * z3 - f3 * z2;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                    } else {
                        // +++0
                        this.addVertex(
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                } else if (f2 < 0) {
                    d0 = f0 - f2;
                    xn0 = f0 * x2 - f2 * x0;
                    yn0 = f0 * y2 - f2 * y0;
                    zn0 = f0 * z2 - f2 * z0;
                    d1 = f1 - f2;
                    xn1 = f1 * x2 - f2 * x1;
                    yn1 = f1 * y2 - f2 * y1;
                    zn1 = f1 * z2 - f2 * z1;

                    if (f3 > 0) {
                        // ++-+
                        d2 = f3 - f2;
                        xn2 = f3 * x2 - f2 * x3;
                        yn2 = f3 * y2 - f2 * y3;
                        zn2 = f3 * z2 - f2 * z3;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                    } else if (f3 < 0) {
                        // ++--
                        d2 = f0 - f3;
                        xn2 = f0 * x3 - f3 * x0;
                        yn2 = f0 * y3 - f3 * y0;
                        zn2 = f0 * z3 - f3 * z0;
                        d3 = f1 - f3;
                        xn3 = f1 * x3 - f3 * x1;
                        yn3 = f1 * y3 - f3 * y1;
                        zn3 = f1 * z3 - f3 * z1;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn3, d3, yn3, d3, zn3, d3),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                    } else {
                        // ++-0
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                } else {
                    if (f3 > 0) {
                        // ++0+
                        this.addVertex(
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
                    } else if (f3 < 0) {
                        // ++0-
                        d0 = f0 - f3;
                        xn0 = f0 * x3 - f3 * x0;
                        yn0 = f0 * y3 - f3 * y0;
                        zn0 = f0 * z3 - f3 * z0;
                        d1 = f1 - f3;
                        xn1 = f1 * x3 - f3 * x1;
                        yn1 = f1 * y3 - f3 * y1;
                        zn1 = f1 * z3 - f3 * z1;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
                    } else {
                        // ++00
                        this.addEdge(
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1),
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                }
            } else if (f1 < 0) {
                if (f2 > 0) {
                    d0 = f0 - f1;
                    xn0 = f0 * x1 - f1 * x0;
                    yn0 = f0 * y1 - f1 * y0;
                    zn0 = f0 * z1 - f1 * z0;
                    d1 = f2 - f1;
                    xn1 = f2 * x1 - f1 * x2;
                    yn1 = f2 * y1 - f1 * y2;
                    zn1 = f2 * z1 - f1 * z2;

                    if (f3 > 0) {
                        // +-++
                        d2 = f3 - f1;
                        xn2 = f3 * x1 - f1 * x3;
                        yn2 = f3 * y1 - f1 * y3;
                        zn2 = f3 * z1 - f1 * z3;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                    } else if (f3 < 0) {
                        // +-+-
                        d2 = f0 - f3;
                        xn2 = f0 * x3 - f3 * x0;
                        yn2 = f0 * y3 - f3 * y0;
                        zn2 = f0 * z3 - f3 * z0;
                        d3 = f2 - f3;
                        xn3 = f2 * x3 - f3 * x2;
                        yn3 = f2 * y3 - f3 * y2;
                        zn3 = f2 * z3 - f3 * z2;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn3, d3, yn3, d3, zn3, d3),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                    } else {
                        // +-+0
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                } else if (f2 < 0) {
                    d0 = f1 - f0;
                    xn0 = f1 * x0 - f0 * x1;
                    yn0 = f1 * y0 - f0 * y1;
                    zn0 = f1 * z0 - f0 * z1;
                    d1 = f2 - f0;
                    xn1 = f2 * x0 - f0 * x2;
                    yn1 = f2 * y0 - f0 * y2;
                    zn1 = f2 * z0 - f0 * z2;

                    if (f3 > 0) {
                        // +--+
                        d2 = f1 - f3;
                        xn2 = f1 * x3 - f3 * x1;
                        yn2 = f1 * y3 - f3 * y1;
                        zn2 = f1 * z3 - f3 * z1;
                        d3 = f2 - f3;
                        xn3 = f2 * x3 - f3 * x2;
                        yn3 = f2 * y3 - f3 * y2;
                        zn3 = f2 * z3 - f3 * z2;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn3, d3, yn3, d3, zn3, d3),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                    } else if (f3 < 0) {
                        // +---
                        d2 = f3 - f0;
                        xn2 = f3 * x0 - f0 * x3;
                        yn2 = f3 * y0 - f0 * y3;
                        zn2 = f3 * z0 - f0 * z3;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(xn2, d2, yn2, d2, zn2, d2));
                    } else {
                        // +--0
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                } else {
                    d0 = f1 - f0;
                    xn0 = f1 * x0 - f0 * x1;
                    yn0 = f1 * y0 - f0 * y1;
                    zn0 = f1 * z0 - f0 * z1;

                    if (f3 > 0) {
                        // +-0+
                        d1 = f1 - f3;
                        xn1 = f1 * x3 - f3 * x1;
                        yn1 = f1 * y3 - f3 * y1;
                        zn1 = f1 * z3 - f3 * z1;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
                    } else if (f3 < 0) {
                        // +-0-
                        d1 = f3 - f0;
                        xn1 = f3 * x0 - f0 * x3;
                        yn1 = f3 * y0 - f0 * y3;
                        zn1 = f3 * z0 - f0 * z3;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
                    } else {
                        // +-00
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1),
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                }
            } else {
                if (f2 > 0) {
                    if (f3 > 0) {
                        // +0++
                        this.addVertex(
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1));
                    } else if (f3 < 0) {
                        // +0+-
                        d0 = f0 - f3;
                        xn0 = f0 * x3 - f3 * x0;
                        yn0 = f0 * y3 - f3 * y0;
                        zn0 = f0 * z3 - f3 * z0;
                        d1 = f2 - f3;
                        xn1 = f2 * x3 - f3 * x2;
                        yn1 = f2 * y3 - f3 * y2;
                        zn1 = f2 * z3 - f3 * z2;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1));
                    } else {
                        // +0+0
                        this.addEdge(
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                } else if (f2 < 0) {
                    d0 = f2 - f0;
                    xn0 = f2 * x0 - f0 * x2;
                    yn0 = f2 * y0 - f0 * y2;
                    zn0 = f2 * z0 - f0 * z2;

                    if (f3 > 0) {
                        // +0-+
                        d1 = f2 - f3;
                        xn1 = f2 * x3 - f3 * x2;
                        yn1 = f2 * y3 - f3 * y2;
                        zn1 = f2 * z3 - f3 * z2;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1));
                    } else if (f3 < 0) {
                        // +0--
                        d1 = f0 - f3;
                        xn1 = f0 * x3 - f3 * x0;
                        yn1 = f0 * y3 - f3 * y0;
                        zn1 = f0 * z3 - f3 * z0;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1));
                    } else {
                        // +0-0
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                } else {
                    if (f3 > 0) {
                        // +00+
                        this.addEdge(
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
                    } else if (f3 < 0) {
                        // +00-
                        d0 = f0 - f3;
                        xn0 = f0 * x3 - f3 * x0;
                        yn0 = f0 * y3 - f3 * y0;
                        zn0 = f0 * z3 - f3 * z0;
                        this.addTriangle(
                            new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
                    } else {
                        // +000
                        this.addTriangle(
                            new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                            new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1),
                            new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                    }
                }
            }
        } else if (f1 !== 0) {
            // convert to case 0+**
            if (f1 < 0) {
                f1 = -f1;
                f2 = -f2;
                f3 = -f3;
            }

            if (f2 > 0) {
                if (f3 > 0) {
                    // 0+++
                    this.addVertex(
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1));
                } else if (f3 < 0) {
                    // 0++-
                    d0 = f2 - f3;
                    xn0 = f2 * x3 - f3 * x2;
                    yn0 = f2 * y3 - f3 * y2;
                    zn0 = f2 * z3 - f3 * z2;
                    d1 = f1 - f3;
                    xn1 = f1 * x3 - f3 * x1;
                    yn1 = f1 * y3 - f3 * y1;
                    zn1 = f1 * z3 - f3 * z1;
                    this.addTriangle(
                        new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                        new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1));
                } else {
                    // 0++0
                    this.addEdge(
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                        new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                }
            } else if (f2 < 0) {
                d0 = f2 - f1;
                xn0 = f2 * x1 - f1 * x2;
                yn0 = f2 * y1 - f1 * y2;
                zn0 = f2 * z1 - f1 * z2;

                if (f3 > 0) {
                    // 0+-+
                    d1 = f2 - f3;
                    xn1 = f2 * x3 - f3 * x2;
                    yn1 = f2 * y3 - f3 * y2;
                    zn1 = f2 * z3 - f3 * z2;
                    this.addTriangle(
                        new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                        new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1));
                } else if (f3 < 0) {
                    // 0+--
                    d1 = f1 - f3;
                    xn1 = f1 * x3 - f3 * x1;
                    yn1 = f1 * y3 - f3 * y1;
                    zn1 = f1 * z3 - f3 * z1;
                    this.addTriangle(
                        new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                        new SurfaceExtractorVertex(xn1, d1, yn1, d1, zn1, d1),
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1));
                } else {
                    // 0+-0
                    this.addTriangle(
                        new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                        new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                }
            } else {
                if (f3 > 0) {
                    // 0+0+
                    this.addEdge(
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                        new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
                } else if (f3 < 0) {
                    // 0+0-
                    d0 = f1 - f3;
                    xn0 = f1 * x3 - f3 * x1;
                    yn0 = f1 * y3 - f3 * y1;
                    zn0 = f1 * z3 - f3 * z1;
                    this.addTriangle(
                        new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                        new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
                } else {
                    // 0+00
                    this.addTriangle(
                        new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                        new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1),
                        new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
                }
            }
        } else if (f2 !== 0) {
            // convert to case 00+*
            if (f2 < 0) {
                f2 = -f2;
                f3 = -f3;
            }

            if (f3 > 0) {
                // 00++
                this.addEdge(
                    new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                    new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1));
            } else if (f3 < 0) {
                // 00+-
                d0 = f2 - f3;
                xn0 = f2 * x3 - f3 * x2;
                yn0 = f2 * y3 - f3 * y2;
                zn0 = f2 * z3 - f3 * z2;
                this.addTriangle(
                    new SurfaceExtractorVertex(xn0, d0, yn0, d0, zn0, d0),
                    new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                    new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1));
            } else {
                // 00+0
                this.addTriangle(
                    new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                    new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                    new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
            }
        } else if (f3 !== 0) {
            // cases 000+ or 000-
            this.addTriangle(
                new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
        } else {
            // case 0000
            this.addTriangle(
                new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1));
            this.addTriangle(
                new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
            this.addTriangle(
                new SurfaceExtractorVertex(x0, 1, y0, 1, z0, 1),
                new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1),
                new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
            this.addTriangle(
                new SurfaceExtractorVertex(x1, 1, y1, 1, z1, 1),
                new SurfaceExtractorVertex(x2, 1, y2, 1, z2, 1),
                new SurfaceExtractorVertex(x3, 1, y3, 1, z3, 1));
        }
    }
}
