// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TetrahedraRasterizer.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Rasterize a set of tetrahedra into a 3D grid: each grid point is
// classified by the index of a tetrahedron containing it (-1 when no
// tetrahedron contains it).
//
// Port notes: the template parameter T (a floating-point type) is number.
// The constructor takes the vertex and tetrahedron arrays directly (upstream
// passes counts and pointers). The upstream operator() is the method
// rasterize(), which returns the grid as an Int32Array instead of filling a
// caller-passed std::vector<int32_t>&. The multithreaded path
// (numThreads > 0 with std::thread) is omitted because JavaScript is
// single-threaded; rasterize() always uses the single-threaded path, which
// produces identical output. As in upstream, rasterize() clips the stored
// tetrahedra bounding boxes against the region in place, so calling
// rasterize() again with a larger region uses the previously clipped boxes.

import { logError } from './Logger';

export class TetrahedraRasterizer {
    // Constructor arguments.
    private mNumVertices: number;
    private mVertices: ReadonlyArray<readonly number[]>;
    private mNumTetrahedra: number;
    private mTetrahedra: ReadonlyArray<readonly number[]>;

    // Axis-aligned bounding boxes for the tetrahedra.
    private mTetraMin: number[][];
    private mTetraMax: number[][];
    private mValid: boolean[];

    // Vertices and axis-aligned bounding boxes in grid coordinates.
    private mGridVertices: number[][];
    private mGridTetraMin: number[][];
    private mGridTetraMax: number[][];

    // The tetrahedra are stored as indexed primitives where the indices are
    // relative to the vertices[] array. The vertices of tetrahedron t are
    // vertices[tetrahedra[t][i]] for 0 <= i < 4. If v0, v1, v2 and v3 are
    // those vertices, the triangle faces have vertices {v0,v2,v1},
    // {v0,v1,v3}, {v0,v3,v2} and {v1,v2,v3}. The faces are counterclockwise
    // ordered when viewed by an observer outside the tetrahedron. The
    // canonical tetrahedron is {v0,v1,v2,v3} where v0 = (0,0,0),
    // v1 = (1,0,0), v2 = (0,1,0) and v3 = (0,0,1).
    constructor(vertices: ReadonlyArray<readonly number[]>,
        tetrahedra: ReadonlyArray<readonly number[]>) {
        if (vertices.length === 0 || tetrahedra.length === 0) {
            logError('Invalid argument.');
        }

        this.mNumVertices = vertices.length;
        this.mVertices = vertices;
        this.mNumTetrahedra = tetrahedra.length;
        this.mTetrahedra = tetrahedra;

        this.mTetraMin = new Array(this.mNumTetrahedra);
        this.mTetraMax = new Array(this.mNumTetrahedra);
        this.mValid = new Array<boolean>(this.mNumTetrahedra).fill(false);
        this.mGridVertices = new Array(this.mNumVertices);
        for (let v = 0; v < this.mNumVertices; ++v) {
            this.mGridVertices[v] = [0, 0, 0];
        }
        this.mGridTetraMin = new Array(this.mNumTetrahedra);
        this.mGridTetraMax = new Array(this.mNumTetrahedra);
        for (let t = 0; t < this.mNumTetrahedra; ++t) {
            this.mGridTetraMin[t] = [0, 0, 0];
            this.mGridTetraMax[t] = [0, 0, 0];
        }

        this.computeTetrahedraAABBs();
    }

    // Rasterize the tetrahedra into a 3D grid. The input region is a box in
    // the coordinate system of the vertices. The box is associated with a 3D
    // grid of specified bounds. A grid point is a 3-tuple (x,y,z) with
    // integer coordinates satisfying 0 <= x < bound[0], 0 <= y < bound[1]
    // and 0 <= z < bound[2]. The point is classified based on whether or not
    // it is contained by a tetrahedron, and the classification is stored as
    // an integer in grid[i] where the index is
    // i = x + bound[0] * (y + bound[1] * z). If the point is not contained
    // by a tetrahedron, grid[i] is set to -1. If the point is contained by a
    // tetrahedron, grid[i] is set to the tetrahedron index t where
    // 0 <= t < numTetrahedra.
    rasterize(regionMin: readonly number[], regionMax: readonly number[],
        bound: readonly number[]): Int32Array {
        if (bound[0] < 2 || bound[1] < 2 || bound[2] < 2) {
            logError('Invalid argument.');
        }

        // Initialize the grid values to -1. When a grid cell is contained
        // in a tetrahedron, the index of that tetrahedron is stored in
        // grid[i]. All such contained grid[] values are nonnegative.
        const grid = new Int32Array(bound[0] * bound[1] * bound[2]);
        grid.fill(-1);

        // Clip-cull the tetrahedra bounding boxes against the region.
        this.clipCullAABBs(regionMin, regionMax);

        // Transform the vertices and tetrahedra bounding boxes to grid
        // coordinates.
        this.transformToGridCoordinates(regionMin, regionMax, bound);

        for (let t = 0; t < this.mNumTetrahedra; ++t) {
            if (this.mValid[t]) {
                this.rasterizeTetrahedron(t, bound, grid);
            }
        }

        return grid;
    }

    // Compute the axis-aligned bounding boxes of the tetrahedra.
    private computeTetrahedraAABBs(): void {
        for (let t = 0; t < this.mNumTetrahedra; ++t) {
            const tetraMin = [...this.mVertices[this.mTetrahedra[t][0]]];
            const tetraMax = [...tetraMin];
            for (let j = 1; j < 4; ++j) {
                const vertex = this.mVertices[this.mTetrahedra[t][j]];
                for (let i = 0; i < 3; ++i) {
                    if (vertex[i] < tetraMin[i]) {
                        tetraMin[i] = vertex[i];
                    } else if (vertex[i] > tetraMax[i]) {
                        tetraMax[i] = vertex[i];
                    }
                }
            }
            this.mTetraMin[t] = tetraMin;
            this.mTetraMax[t] = tetraMax;
        }
    }

    // Clip-cull the tetrahedra bounding boxes against the region. The
    // mValid[t] is true whenever the box of tetrahedron t intersects the
    // region.
    private clipCullAABBs(regionMin: readonly number[],
        regionMax: readonly number[]): void {
        for (let t = 0; t < this.mNumTetrahedra; ++t) {
            const tetraMin = this.mTetraMin[t];
            const tetraMax = this.mTetraMax[t];
            this.mValid[t] = true;
            for (let i = 0; i < 3; ++i) {
                tetraMin[i] = Math.max(tetraMin[i], regionMin[i]);
                tetraMax[i] = Math.min(tetraMax[i], regionMax[i]);
                if (tetraMin[i] > tetraMax[i]) {
                    this.mValid[t] = false;
                }
            }
        }
    }

    private transformToGridCoordinates(regionMin: readonly number[],
        regionMax: readonly number[], bound: readonly number[]): void {
        const multiplier = [0, 0, 0];
        for (let i = 0; i < 3; ++i) {
            multiplier[i] = (bound[i] - 1) / (regionMax[i] - regionMin[i]);
        }

        for (let v = 0; v < this.mNumVertices; ++v) {
            const vertex = this.mVertices[v];
            const gridVertex = this.mGridVertices[v];
            for (let i = 0; i < 3; ++i) {
                gridVertex[i] = multiplier[i] * (vertex[i] - regionMin[i]);
            }
        }

        for (let t = 0; t < this.mNumTetrahedra; ++t) {
            const tetraMin = this.mTetraMin[t];
            const tetraMax = this.mTetraMax[t];
            const gridTetraMin = this.mGridTetraMin[t];
            const gridTetraMax = this.mGridTetraMax[t];
            for (let i = 0; i < 3; ++i) {
                gridTetraMin[i] = Math.ceil(
                    multiplier[i] * (tetraMin[i] - regionMin[i]));
                gridTetraMax[i] = Math.floor(
                    multiplier[i] * (tetraMax[i] - regionMin[i]));
            }
        }
    }

    private rasterizeTetrahedron(t: number, bound: readonly number[],
        grid: Int32Array): void {
        const imin = this.mGridTetraMin[t];
        const imax = this.mGridTetraMax[t];
        const v0 = this.mGridVertices[this.mTetrahedra[t][0]];
        const v1 = this.mGridVertices[this.mTetrahedra[t][1]];
        const v2 = this.mGridVertices[this.mTetrahedra[t][2]];
        const v3 = this.mGridVertices[this.mTetrahedra[t][3]];
        const gridP = [0, 0, 0];

        for (let i2 = imin[2]; i2 <= imax[2]; ++i2) {
            gridP[2] = i2;
            for (let i1 = imin[1]; i1 <= imax[1]; ++i1) {
                gridP[1] = i1;

                // Find the smallest x-index of a contained grid point on
                // this grid row.
                let i0min: number;
                for (i0min = imin[0]; i0min <= imax[0]; ++i0min) {
                    gridP[0] = i0min;
                    if (TetrahedraRasterizer.pointInTetrahedron(gridP, v0, v1, v2, v3)) {
                        break;
                    }
                }
                if (i0min > imax[0]) {
                    continue;
                }

                // Find the largest x-index of a contained grid point on
                // this grid row. The tetrahedron is convex, so all points
                // between the extremes are contained.
                let i0max: number;
                for (i0max = imax[0]; i0max >= i0min; --i0max) {
                    gridP[0] = i0max;
                    if (TetrahedraRasterizer.pointInTetrahedron(gridP, v0, v1, v2, v3)) {
                        break;
                    }
                }

                const base = bound[0] * (i1 + bound[1] * i2);
                for (let i0 = i0min, j = i0 + base; i0 <= i0max; ++i0, ++j) {
                    grid[j] = t;
                }
            }
        }
    }

    // The point is contained when it is on the negative side of (or exactly
    // on) each of the four face planes, so grid points on the boundary of a
    // tetrahedron are classified as contained.
    private static pointInTetrahedron(P: readonly number[],
        V0: readonly number[], V1: readonly number[],
        V2: readonly number[], V3: readonly number[]): boolean {
        const PmV0 = TetrahedraRasterizer.sub(P, V0);
        const V1mV0 = TetrahedraRasterizer.sub(V1, V0);
        const V2mV0 = TetrahedraRasterizer.sub(V2, V0);
        if (TetrahedraRasterizer.dotCross(PmV0, V2mV0, V1mV0) > 0) {
            return false;
        }

        const V3mV0 = TetrahedraRasterizer.sub(V3, V0);
        if (TetrahedraRasterizer.dotCross(PmV0, V1mV0, V3mV0) > 0) {
            return false;
        }

        if (TetrahedraRasterizer.dotCross(PmV0, V3mV0, V2mV0) > 0) {
            return false;
        }

        const PmV1 = TetrahedraRasterizer.sub(P, V1);
        const V2mV1 = TetrahedraRasterizer.sub(V2, V1);
        const V3mV1 = TetrahedraRasterizer.sub(V3, V1);
        if (TetrahedraRasterizer.dotCross(PmV1, V2mV1, V3mV1) > 0) {
            return false;
        }

        return true;
    }

    private static sub(U: readonly number[], V: readonly number[]): number[] {
        return [U[0] - V[0], U[1] - V[1], U[2] - V[2]];
    }

    private static dot(U: readonly number[], V: readonly number[]): number {
        return U[0] * V[0] + U[1] * V[1] + U[2] * V[2];
    }

    private static cross(U: readonly number[], V: readonly number[]): number[] {
        return [
            U[1] * V[2] - U[2] * V[1],
            U[2] * V[0] - U[0] * V[2],
            U[0] * V[1] - U[1] * V[0]
        ];
    }

    private static dotCross(U: readonly number[], V: readonly number[],
        W: readonly number[]): number {
        return TetrahedraRasterizer.dot(U, TetrahedraRasterizer.cross(V, W));
    }
}
