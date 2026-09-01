// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SurfaceExtractorCubes.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The level set extraction algorithm implemented here is described in
// Section 3.2 of the document
// https://www.geometrictools.com/Documentation/LevelSetExtraction.pdf
//
// Port notes: T and Real are both number; the voxel values are expected to be
// integers, and the internal int64_t computations map onto IEEE doubles (the
// products stay exact while they remain within the safe-integer range). The
// upstream pure virtual Extract(level, vertices&, triangles&) is the base
// class method extractRational(level), which returns the two arrays. The
// nested VETable/TVertex classes are module-private. LogError becomes
// logError from src/Logger.ts.

import { logError } from './Logger';
import {
    SurfaceExtractor,
    SurfaceExtractorTriangle,
    SurfaceExtractorVertex
} from './SurfaceExtractor';

// Indices of the 12 voxel edges and the 6 voxel faces in the VETable. An
// edge name lists the two coordinates that are constant on that edge; for
// example, EI_XMIN_YMIN is the edge with x = xmin and y = ymin, along which
// z varies.
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

// The bit masks for the edge indices, used to accumulate the voxel 'type'.
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

// A vertex of the voxel wireframe, with up to 4 adjacent wireframe vertices.
class TVertex {
    pos: SurfaceExtractorVertex;
    numAdjacents: number;
    adj: number[];
    valid: boolean;

    constructor() {
        this.pos = new SurfaceExtractorVertex();
        this.numAdjacents = 0;
        this.adj = [0, 0, 0, 0];
        this.valid = false;
    }
}

// Vertex-edge-triangle table to support mesh topology. There are 12 slots
// for the voxel-edge vertices and 6 slots for the face branch points.
class VETable {
    private mVertex: TVertex[];

    constructor() {
        this.mVertex = new Array<TVertex>(18);
        for (let i = 0; i < 18; ++i) {
            this.mVertex[i] = new TVertex();
        }
    }

    isValidVertex(i: number): boolean {
        return this.mVertex[i].valid;
    }

    getXN(i: number): number {
        return this.mVertex[i].pos.xNumer;
    }

    getXD(i: number): number {
        return this.mVertex[i].pos.xDenom;
    }

    getYN(i: number): number {
        return this.mVertex[i].pos.yNumer;
    }

    getYD(i: number): number {
        return this.mVertex[i].pos.yDenom;
    }

    getZN(i: number): number {
        return this.mVertex[i].pos.zNumer;
    }

    getZD(i: number): number {
        return this.mVertex[i].pos.zDenom;
    }

    // Upstream Insert(i, pos).
    insertVertex(i: number, pos: SurfaceExtractorVertex): void {
        const vertex = this.mVertex[i];
        vertex.pos = pos;
        vertex.valid = true;
    }

    // Upstream Insert(i0, i1).
    insertEdge(i0: number, i1: number): void {
        const vertex0 = this.mVertex[i0];
        const vertex1 = this.mVertex[i1];
        vertex0.adj[vertex0.numAdjacents++] = i1;
        vertex1.adj[vertex1.numAdjacents++] = i0;
    }

    // Ear-clip the wireframe to get the triangles, appending the results to
    // the caller's arrays. Each triangle contributes three new vertices, so
    // the arrays contain duplicates (removed later by makeUnique).
    removeTriangles(vertices: SurfaceExtractorVertex[],
        triangles: SurfaceExtractorTriangle[]): void {
        for (; ;) {
            const triangle = this.remove();
            if (triangle === null) {
                break;
            }
            const v0 = vertices.length;
            const v1 = v0 + 1;
            const v2 = v1 + 1;
            triangles.push(new SurfaceExtractorTriangle(v0, v1, v2));
            vertices.push(this.mVertex[triangle[0]].pos);
            vertices.push(this.mVertex[triangle[1]].pos);
            vertices.push(this.mVertex[triangle[2]].pos);
        }
    }

    private removeVertex(i: number): void {
        const vertex0 = this.mVertex[i];
        // assert: vertex0.numAdjacents == 2
        const a0 = vertex0.adj[0];
        const a1 = vertex0.adj[1];
        const adjVertex0 = this.mVertex[a0];
        const adjVertex1 = this.mVertex[a1];

        for (let j = 0; j < adjVertex0.numAdjacents; ++j) {
            if (adjVertex0.adj[j] === i) {
                adjVertex0.adj[j] = a1;
                break;
            }
        }

        for (let j = 0; j < adjVertex1.numAdjacents; ++j) {
            if (adjVertex1.adj[j] === i) {
                adjVertex1.adj[j] = a0;
                break;
            }
        }

        vertex0.valid = false;

        if (adjVertex0.numAdjacents === 2) {
            if (adjVertex0.adj[0] === adjVertex0.adj[1]) {
                adjVertex0.valid = false;
            }
        }

        if (adjVertex1.numAdjacents === 2) {
            if (adjVertex1.adj[0] === adjVertex1.adj[1]) {
                adjVertex1.valid = false;
            }
        }
    }

    // Upstream Remove(triangle&) returning bool; the port returns the triple
    // of table indices or null when no ear remains.
    private remove(): [number, number, number] | null {
        for (let i = 0; i < 18; ++i) {
            const vertex = this.mVertex[i];
            if (vertex.valid && vertex.numAdjacents === 2) {
                const triangle: [number, number, number] =
                    [i, vertex.adj[0], vertex.adj[1]];
                this.removeVertex(i);
                return triangle;
            }
        }
        return null;
    }
}

export class SurfaceExtractorCubes extends SurfaceExtractor {
    // The input is a 3D image with lexicographically ordered voxels (x,y,z)
    // stored in a linear array. Voxel (x,y,z) is stored in the array at
    // location index = x + xBound * (y + yBound * z). The inputs xBound,
    // yBound and zBound must each be 2 or larger so that there is at least
    // one image cube to process. The inputVoxels must contain at least
    // xBound * yBound * zBound elements.
    constructor(xBound: number, yBound: number, zBound: number,
        inputVoxels: ArrayLike<number>) {
        super(xBound, yBound, zBound, inputVoxels);
    }

    // Extract level surfaces and return rational vertices. Use the
    // base-class extract() if you want real-valued vertices.
    override extractRational(level: number): {
        vertices: SurfaceExtractorVertex[];
        triangles: SurfaceExtractorTriangle[];
    } {
        // Adjust the image so that the level set is F(x,y,z) = 0. The
        // precondition for 'level' is that it is not exactly a voxel value.
        // However, T is an integer type, so we cannot pass in a 'level' that
        // has a fractional value. To circumvent this, the voxel values are
        // doubled so that they are even integers. The level value is doubled
        // and 1 added to obtain an odd integer, guaranteeing 'level' is not a
        // voxel value.
        const levelI64 = 2 * level + 1;
        for (let i = 0; i < this.mVoxels.length; ++i) {
            const inputI64 = 2 * this.mInputVoxels[i];
            this.mVoxels[i] = inputI64 - levelI64;
        }

        const vertices: SurfaceExtractorVertex[] = [];
        const triangles: SurfaceExtractorTriangle[] = [];
        for (let z = 0; z < this.mZBound - 1; ++z) {
            for (let y = 0; y < this.mYBound - 1; ++y) {
                for (let x = 0; x < this.mXBound - 1; ++x) {
                    // Get vertices on edges of box (if any).
                    const table = new VETable();
                    const type = this.getVertices(x, y, z, table);
                    if (type !== 0) {
                        // Get edges on faces of box.
                        this.getXMinEdges(x, y, z, type, table);
                        this.getXMaxEdges(x, y, z, type, table);
                        this.getYMinEdges(x, y, z, type, table);
                        this.getYMaxEdges(x, y, z, type, table);
                        this.getZMinEdges(x, y, z, type, table);
                        this.getZMaxEdges(x, y, z, type, table);

                        // Ear-clip the wireframe mesh.
                        table.removeTriangles(vertices, triangles);
                    }
                }
            }
        }

        return { vertices, triangles };
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
        const oneMX = 1 - dx;
        const oneMY = 1 - dy;
        const oneMZ = 1 - dz;

        const grad: [number, number, number] = [0, 0, 0];

        let tmp0 = oneMY * (f100 - f000) + dy * (f110 - f010);
        let tmp1 = oneMY * (f101 - f001) + dy * (f111 - f011);
        grad[0] = oneMZ * tmp0 + dz * tmp1;

        tmp0 = oneMX * (f010 - f000) + dx * (f110 - f100);
        tmp1 = oneMX * (f011 - f001) + dx * (f111 - f101);
        grad[1] = oneMZ * tmp0 + dz * tmp1;

        tmp0 = oneMX * (f001 - f000) + dx * (f101 - f100);
        tmp1 = oneMX * (f011 - f010) + dx * (f111 - f110);
        grad[2] = oneMY * tmp0 + dy * tmp1;

        return grad;
    }

    // Insert into the table a vertex for each voxel edge whose endpoint
    // values have opposite signs. The return value is the bitmask of the
    // edges that were assigned vertices.
    protected getVertices(x: number, y: number, z: number,
        table: VETable): number {
        let type = 0;

        // Get the image values at the corners of the voxel.
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

        const x0 = x;
        const x1 = x0 + 1;
        const y0 = y;
        const y1 = y0 + 1;
        const z0 = z;
        const z1 = z0 + 1;
        let d = 0;

        // xmin-ymin edge
        if (f000 * f001 < 0) {
            type |= EB_XMIN_YMIN;
            d = f001 - f000;
            table.insertVertex(EI_XMIN_YMIN,
                new SurfaceExtractorVertex(x0, 1, y0, 1, z0 * d - f000, d));
        }

        // xmin-ymax edge
        if (f010 * f011 < 0) {
            type |= EB_XMIN_YMAX;
            d = f011 - f010;
            table.insertVertex(EI_XMIN_YMAX,
                new SurfaceExtractorVertex(x0, 1, y1, 1, z0 * d - f010, d));
        }

        // xmax-ymin edge
        if (f100 * f101 < 0) {
            type |= EB_XMAX_YMIN;
            d = f101 - f100;
            table.insertVertex(EI_XMAX_YMIN,
                new SurfaceExtractorVertex(x1, 1, y0, 1, z0 * d - f100, d));
        }

        // xmax-ymax edge
        if (f110 * f111 < 0) {
            type |= EB_XMAX_YMAX;
            d = f111 - f110;
            table.insertVertex(EI_XMAX_YMAX,
                new SurfaceExtractorVertex(x1, 1, y1, 1, z0 * d - f110, d));
        }

        // xmin-zmin edge
        if (f000 * f010 < 0) {
            type |= EB_XMIN_ZMIN;
            d = f010 - f000;
            table.insertVertex(EI_XMIN_ZMIN,
                new SurfaceExtractorVertex(x0, 1, y0 * d - f000, d, z0, 1));
        }

        // xmin-zmax edge
        if (f001 * f011 < 0) {
            type |= EB_XMIN_ZMAX;
            d = f011 - f001;
            table.insertVertex(EI_XMIN_ZMAX,
                new SurfaceExtractorVertex(x0, 1, y0 * d - f001, d, z1, 1));
        }

        // xmax-zmin edge
        if (f100 * f110 < 0) {
            type |= EB_XMAX_ZMIN;
            d = f110 - f100;
            table.insertVertex(EI_XMAX_ZMIN,
                new SurfaceExtractorVertex(x1, 1, y0 * d - f100, d, z0, 1));
        }

        // xmax-zmax edge
        if (f101 * f111 < 0) {
            type |= EB_XMAX_ZMAX;
            d = f111 - f101;
            table.insertVertex(EI_XMAX_ZMAX,
                new SurfaceExtractorVertex(x1, 1, y0 * d - f101, d, z1, 1));
        }

        // ymin-zmin edge
        if (f000 * f100 < 0) {
            type |= EB_YMIN_ZMIN;
            d = f100 - f000;
            table.insertVertex(EI_YMIN_ZMIN,
                new SurfaceExtractorVertex(x0 * d - f000, d, y0, 1, z0, 1));
        }

        // ymin-zmax edge
        if (f001 * f101 < 0) {
            type |= EB_YMIN_ZMAX;
            d = f101 - f001;
            table.insertVertex(EI_YMIN_ZMAX,
                new SurfaceExtractorVertex(x0 * d - f001, d, y0, 1, z1, 1));
        }

        // ymax-zmin edge
        if (f010 * f110 < 0) {
            type |= EB_YMAX_ZMIN;
            d = f110 - f010;
            table.insertVertex(EI_YMAX_ZMIN,
                new SurfaceExtractorVertex(x0 * d - f010, d, y1, 1, z0, 1));
        }

        // ymax-zmax edge
        if (f011 * f111 < 0) {
            type |= EB_YMAX_ZMAX;
            d = f111 - f011;
            table.insertVertex(EI_YMAX_ZMAX,
                new SurfaceExtractorVertex(x0 * d - f011, d, y1, 1, z1, 1));
        }

        return type;
    }

    protected getXMinEdges(x: number, y: number, z: number, type: number,
        table: VETable): void {
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
                break;
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
                let i = x + this.mXBound * (y + this.mYBound * z);
                // F(x,y,z)
                const f00 = this.mVoxels[i];
                i += this.mXBound;
                // F(x,y+1,z)
                const f10 = this.mVoxels[i];
                i += this.mXYBound;
                // F(x,y+1,z+1)
                const f11 = this.mVoxels[i];
                i -= this.mXBound;
                // F(x,y,z+1)
                const f01 = this.mVoxels[i];
                const det = f00 * f11 - f01 * f10;

                if (det > 0) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_YMIN, EI_XMIN_ZMIN);
                    table.insertEdge(EI_XMIN_YMAX, EI_XMIN_ZMAX);
                } else if (det < 0) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_YMIN, EI_XMIN_ZMAX);
                    table.insertEdge(EI_XMIN_YMAX, EI_XMIN_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to the
                    // tessellation.
                    table.insertVertex(FI_XMIN, new SurfaceExtractorVertex(
                        table.getXN(EI_XMIN_ZMIN), table.getXD(EI_XMIN_ZMIN),
                        table.getYN(EI_XMIN_ZMIN), table.getYD(EI_XMIN_ZMIN),
                        table.getZN(EI_XMIN_YMIN), table.getZD(EI_XMIN_YMIN)));

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_YMIN, FI_XMIN);
                    table.insertEdge(EI_XMIN_YMAX, FI_XMIN);
                    table.insertEdge(EI_XMIN_ZMIN, FI_XMIN);
                    table.insertEdge(EI_XMIN_ZMAX, FI_XMIN);
                }
                break;
            }
            default:
                logError('Unexpected condition.');
        }
    }

    protected getXMaxEdges(x: number, y: number, z: number, type: number,
        table: VETable): void {
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
                break;
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
                let i = (x + 1) + this.mXBound * (y + this.mYBound * z);
                // F(x,y,z)
                const f00 = this.mVoxels[i];
                i += this.mXBound;
                // F(x,y+1,z)
                const f10 = this.mVoxels[i];
                i += this.mXYBound;
                // F(x,y+1,z+1)
                const f11 = this.mVoxels[i];
                i -= this.mXBound;
                // F(x,y,z+1)
                const f01 = this.mVoxels[i];
                const det = f00 * f11 - f01 * f10;

                if (det > 0) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMAX_YMIN, EI_XMAX_ZMIN);
                    table.insertEdge(EI_XMAX_YMAX, EI_XMAX_ZMAX);
                } else if (det < 0) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMAX_YMIN, EI_XMAX_ZMAX);
                    table.insertEdge(EI_XMAX_YMAX, EI_XMAX_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to the
                    // tessellation.
                    table.insertVertex(FI_XMAX, new SurfaceExtractorVertex(
                        table.getXN(EI_XMAX_ZMIN), table.getXD(EI_XMAX_ZMIN),
                        table.getYN(EI_XMAX_ZMIN), table.getYD(EI_XMAX_ZMIN),
                        table.getZN(EI_XMAX_YMIN), table.getZD(EI_XMAX_YMIN)));

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMAX_YMIN, FI_XMAX);
                    table.insertEdge(EI_XMAX_YMAX, FI_XMAX);
                    table.insertEdge(EI_XMAX_ZMIN, FI_XMAX);
                    table.insertEdge(EI_XMAX_ZMAX, FI_XMAX);
                }
                break;
            }
            default:
                logError('Unexpected condition.');
        }
    }

    protected getYMinEdges(x: number, y: number, z: number, type: number,
        table: VETable): void {
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
                break;
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
                let i = x + this.mXBound * (y + this.mYBound * z);
                // F(x,y,z)
                const f00 = this.mVoxels[i];
                ++i;
                // F(x+1,y,z)
                const f10 = this.mVoxels[i];
                i += this.mXYBound;
                // F(x+1,y,z+1)
                const f11 = this.mVoxels[i];
                --i;
                // F(x,y,z+1)
                const f01 = this.mVoxels[i];
                const det = f00 * f11 - f01 * f10;

                if (det > 0) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_YMIN, EI_YMIN_ZMIN);
                    table.insertEdge(EI_XMAX_YMIN, EI_YMIN_ZMAX);
                } else if (det < 0) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_YMIN, EI_YMIN_ZMAX);
                    table.insertEdge(EI_XMAX_YMIN, EI_YMIN_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to the
                    // tessellation.
                    table.insertVertex(FI_YMIN, new SurfaceExtractorVertex(
                        table.getXN(EI_YMIN_ZMIN), table.getXD(EI_YMIN_ZMIN),
                        table.getYN(EI_XMIN_YMIN), table.getYD(EI_XMIN_YMIN),
                        table.getZN(EI_XMIN_YMIN), table.getZD(EI_XMIN_YMIN)));

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_YMIN, FI_YMIN);
                    table.insertEdge(EI_XMAX_YMIN, FI_YMIN);
                    table.insertEdge(EI_YMIN_ZMIN, FI_YMIN);
                    table.insertEdge(EI_YMIN_ZMAX, FI_YMIN);
                }
                break;
            }
            default:
                logError('Unexpected condition.');
        }
    }

    protected getYMaxEdges(x: number, y: number, z: number, type: number,
        table: VETable): void {
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
                break;
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
                let i = x + this.mXBound * ((y + 1) + this.mYBound * z);
                // F(x,y,z)
                const f00 = this.mVoxels[i];
                ++i;
                // F(x+1,y,z)
                const f10 = this.mVoxels[i];
                i += this.mXYBound;
                // F(x+1,y,z+1)
                const f11 = this.mVoxels[i];
                --i;
                // F(x,y,z+1)
                const f01 = this.mVoxels[i];
                const det = f00 * f11 - f01 * f10;

                if (det > 0) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_YMAX, EI_YMAX_ZMIN);
                    table.insertEdge(EI_XMAX_YMAX, EI_YMAX_ZMAX);
                } else if (det < 0) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_YMAX, EI_YMAX_ZMAX);
                    table.insertEdge(EI_XMAX_YMAX, EI_YMAX_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to the
                    // tessellation.
                    table.insertVertex(FI_YMAX, new SurfaceExtractorVertex(
                        table.getXN(EI_YMAX_ZMIN), table.getXD(EI_YMAX_ZMIN),
                        table.getYN(EI_XMIN_YMAX), table.getYD(EI_XMIN_YMAX),
                        table.getZN(EI_XMIN_YMAX), table.getZD(EI_XMIN_YMAX)));

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_YMAX, FI_YMAX);
                    table.insertEdge(EI_XMAX_YMAX, FI_YMAX);
                    table.insertEdge(EI_YMAX_ZMIN, FI_YMAX);
                    table.insertEdge(EI_YMAX_ZMAX, FI_YMAX);
                }
                break;
            }
            default:
                logError('Unexpected condition.');
        }
    }

    protected getZMinEdges(x: number, y: number, z: number, type: number,
        table: VETable): void {
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
                break;
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
                let i = x + this.mXBound * (y + this.mYBound * z);
                // F(x,y,z)
                const f00 = this.mVoxels[i];
                ++i;
                // F(x+1,y,z)
                const f10 = this.mVoxels[i];
                i += this.mXBound;
                // F(x+1,y+1,z)
                const f11 = this.mVoxels[i];
                --i;
                // F(x,y+1,z)
                const f01 = this.mVoxels[i];
                const det = f00 * f11 - f01 * f10;

                if (det > 0) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_ZMIN, EI_YMIN_ZMIN);
                    table.insertEdge(EI_XMAX_ZMIN, EI_YMAX_ZMIN);
                } else if (det < 0) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_ZMIN, EI_YMAX_ZMIN);
                    table.insertEdge(EI_XMAX_ZMIN, EI_YMIN_ZMIN);
                } else {
                    // Plus-sign configuration, add branch point to the
                    // tessellation.
                    table.insertVertex(FI_ZMIN, new SurfaceExtractorVertex(
                        table.getXN(EI_YMIN_ZMIN), table.getXD(EI_YMIN_ZMIN),
                        table.getYN(EI_XMIN_ZMIN), table.getYD(EI_XMIN_ZMIN),
                        table.getZN(EI_XMIN_ZMIN), table.getZD(EI_XMIN_ZMIN)));

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_ZMIN, FI_ZMIN);
                    table.insertEdge(EI_XMAX_ZMIN, FI_ZMIN);
                    table.insertEdge(EI_YMIN_ZMIN, FI_ZMIN);
                    table.insertEdge(EI_YMAX_ZMIN, FI_ZMIN);
                }
                break;
            }
            default:
                logError('Unexpected condition.');
        }
    }

    protected getZMaxEdges(x: number, y: number, z: number, type: number,
        table: VETable): void {
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
                break;
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
                let i = x + this.mXBound * (y + this.mYBound * (z + 1));
                // F(x,y,z)
                const f00 = this.mVoxels[i];
                ++i;
                // F(x+1,y,z)
                const f10 = this.mVoxels[i];
                i += this.mXBound;
                // F(x+1,y+1,z)
                const f11 = this.mVoxels[i];
                --i;
                // F(x,y+1,z)
                const f01 = this.mVoxels[i];
                const det = f00 * f11 - f01 * f10;

                if (det > 0) {
                    // Disjoint hyperbolic segments, pair <P0,P2>, <P1,P3>.
                    table.insertEdge(EI_XMIN_ZMAX, EI_YMIN_ZMAX);
                    table.insertEdge(EI_XMAX_ZMAX, EI_YMAX_ZMAX);
                } else if (det < 0) {
                    // Disjoint hyperbolic segments, pair <P0,P3>, <P1,P2>.
                    table.insertEdge(EI_XMIN_ZMAX, EI_YMAX_ZMAX);
                    table.insertEdge(EI_XMAX_ZMAX, EI_YMIN_ZMAX);
                } else {
                    // Plus-sign configuration, add branch point to the
                    // tessellation.
                    table.insertVertex(FI_ZMAX, new SurfaceExtractorVertex(
                        table.getXN(EI_YMIN_ZMAX), table.getXD(EI_YMIN_ZMAX),
                        table.getYN(EI_XMIN_ZMAX), table.getYD(EI_XMIN_ZMAX),
                        table.getZN(EI_XMIN_ZMAX), table.getZD(EI_XMIN_ZMAX)));

                    // Add edges sharing the branch point.
                    table.insertEdge(EI_XMIN_ZMAX, FI_ZMAX);
                    table.insertEdge(EI_XMAX_ZMAX, FI_ZMAX);
                    table.insertEdge(EI_YMIN_ZMAX, FI_ZMAX);
                    table.insertEdge(EI_YMAX_ZMAX, FI_ZMAX);
                }
                break;
            }
            default:
                logError('Unexpected condition.');
        }
    }
}
