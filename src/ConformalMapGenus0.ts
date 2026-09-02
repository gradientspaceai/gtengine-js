// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ConformalMapGenus0.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Conformally map a 2-dimensional manifold mesh with the topology of a
// sphere to a sphere. The algorithm is an implementation of the one in the
// paper
//    S. Haker, S. Angenent, A. Tannenbaum, R. Kikinis, G. Sapiro and
//    M. Halle, Conformal surface parameterization for texture mapping,
//    IEEE Transactions on Visualization and Computer Graphics,
//    Volume 6, Number 2, pages 181-189, 2000
// The paper is available at https://ieeexplore.ieee.org/document/856998 but
// is not freely downloadable.
//
// Port notes:
// - 'operator()' becomes 'compute(...)', per PORTING.md. Upstream's
//   (numPositions, positions) and (numTriangles, indices) count-plus-pointer
//   pairs become the arrays alone, whose lengths carry the counts.
// - The sparse matrix A is upstream's std::map keyed by (row,col); the port
//   uses the LinearSystemSparseEntry array of LinearSystem.ts. The entries
//   are appended in the same (row,col) increasing order the std::map would
//   iterate, and LinearSystem.solveSymmetricCGSparse re-sorts them anyway,
//   so the floating-point accumulation order matches upstream.
// - Upstream reads element.first.V[], the vertices of the *unordered edge
//   key*, so v0 < v1. The port's ETManifoldMesh::Edge stores the vertices in
//   the order first encountered, so the port takes min/max explicitly.
// - Upstream dereferences edge->T[1] without checking it, which is a null
//   pointer for a boundary edge. The port asserts that the mesh is closed
//   rather than throwing a TypeError deep inside the loop.
// - The returned Vector2/Vector3 arrays are the port's runtime-dimension
//   Vector of sizes 2 and 3.

import { GTE_C_INV_LN_2, GTE_C_LN_10, GTE_C_PI } from './Constants';
import { ETManifoldMesh } from './ETManifoldMesh';
import { LinearSystem } from './LinearSystem';
import type { LinearSystemSparseEntry } from './LinearSystem';
import { logAssert } from './Logger';
import { Polynomial1 } from './Polynomial1';
import { Vector, dot, length, sub } from './Vector';
import { cross } from './Vector3';

export class ConformalMapGenus0 {
    // Conformal mapping to a plane. The plane's (px,py) points correspond to
    // the mesh's (mx,my,mz) points.
    private mPlaneCoordinates: Vector[];
    private mMinPlaneCoordinate: Vector;
    private mMaxPlaneCoordinate: Vector;

    // Conformal mapping to a sphere. The sphere's (sx,sy,sz) points
    // correspond to the mesh's (mx,my,mz) points.
    private mSphereCoordinates: Vector[];
    private mSphereRadius: number;

    constructor() {
        this.mPlaneCoordinates = [];
        this.mMinPlaneCoordinate = new Vector(2);
        this.mMaxPlaneCoordinate = new Vector(2);
        this.mSphereCoordinates = [];
        this.mSphereRadius = 0;
    }

    // The input mesh should be a closed, manifold surface that has the
    // topology of a sphere (a genus-0 surface). The 'indices' array has
    // 3 * numTriangles elements, three vertex indices per triangle.
    //
    // The returned value is 'true' whenever the conjugate gradient algorithm
    // converged. Even if it did not, the results might still be acceptable.
    compute(positions: readonly Vector[], indices: readonly number[],
        punctureTriangle: number): boolean {
        const numPositions = positions.length;
        const numTriangles = indices.length / 3;
        logAssert(numPositions > 0, 'The mesh must have positions.');
        logAssert(Number.isInteger(numTriangles) && numTriangles > 0,
            'The index array must have 3 indices per triangle.');
        logAssert(Number.isInteger(punctureTriangle) && punctureTriangle >= 0 &&
            punctureTriangle < numTriangles, 'Invalid puncture triangle.');

        let converged = true;
        this.mPlaneCoordinates = new Array<Vector>(numPositions);
        this.mSphereCoordinates = new Array<Vector>(numPositions);
        for (let i = 0; i < numPositions; ++i) {
            this.mPlaneCoordinates[i] = new Vector(2);
            this.mSphereCoordinates[i] = new Vector(3);
        }

        // Construct a triangle-edge representation of the mesh.
        const graph = new ETManifoldMesh();
        for (let t = 0; t < numTriangles; ++t) {
            graph.insert(indices[3 * t], indices[3 * t + 1], indices[3 * t + 2]);
        }
        const emap = graph.getEdges();

        // Construct the nondiagonal entries of the sparse matrix A. The
        // weight of edge <v0,v1> is minus one half of the sum of the
        // cotangents of the two angles opposite the edge.
        const A: LinearSystemSparseEntry[] = [];
        for (const element of emap) {
            const v0 = Math.min(element.V[0], element.V[1]);
            const v1 = Math.max(element.V[0], element.V[1]);

            let value = 0;
            for (let j = 0; j < 2; ++j) {
                const triangle = element.T[j];
                logAssert(triangle !== null,
                    'The mesh must be a closed manifold surface.');
                for (let i = 0; i < 3; ++i) {
                    const v2 = triangle.V[i];
                    if (v2 !== v0 && v2 !== v1) {
                        const E0 = sub(positions[v0], positions[v2]);
                        const E1 = sub(positions[v1], positions[v2]);
                        value += dot(E0, E1) / length(cross(E0, E1));
                    }
                }
            }

            value *= -0.5;
            A.push({ row: v0, col: v1, value });
        }

        // Construct the diagonal entries of the sparse matrix A.
        const tmp = new Array<number>(numPositions).fill(0);
        for (const element of A) {
            tmp[element.row] -= element.value;
            tmp[element.col] -= element.value;
        }
        for (let i = 0; i < numPositions; ++i) {
            A.push({ row: i, col: i, value: tmp[i] });
        }
        logAssert(numPositions + emap.length === A.length, 'Mismatched sizes.');

        // Construct the sparse column vector B.
        const v0 = indices[3 * punctureTriangle];
        const v1 = indices[3 * punctureTriangle + 1];
        const v2 = indices[3 * punctureTriangle + 2];
        const V0 = positions[v0];
        const V1 = positions[v1];
        const V2 = positions[v2];
        const E10 = sub(V1, V0);
        const E20 = sub(V2, V0);
        const E12 = sub(V1, V2);
        const normal = cross(E20, E10);
        const len10 = length(E10);
        const invLen10 = 1 / len10;
        const twoArea = length(normal);
        const invLenNormal = 1 / twoArea;
        const invProd = invLen10 * invLenNormal;
        const re0 = -invLen10;
        const im0 = invProd * dot(E12, E10);
        const re1 = invLen10;
        const im1 = invProd * dot(E20, E10);
        const re2 = 0;
        const im2 = -len10 * invLenNormal;

        // Solve the sparse system for the real parts.
        const maxIterations = 1024;
        const tolerance = 1e-06;
        tmp.fill(0);
        tmp[v0] = re0;
        tmp[v1] = re1;
        tmp[v2] = re2;
        let solution = LinearSystem.solveSymmetricCGSparse(numPositions, A, tmp,
            maxIterations, tolerance);
        if (solution.iterations >= maxIterations) {
            converged = false;
        }
        for (let i = 0; i < numPositions; ++i) {
            this.mPlaneCoordinates[i].values[0] = solution.X[i];
        }

        // Solve the sparse system for the imaginary parts.
        tmp.fill(0);
        tmp[v0] = -im0;
        tmp[v1] = -im1;
        tmp[v2] = -im2;
        solution = LinearSystem.solveSymmetricCGSparse(numPositions, A, tmp,
            maxIterations, tolerance);
        if (solution.iterations >= maxIterations) {
            converged = false;
        }
        for (let i = 0; i < numPositions; ++i) {
            this.mPlaneCoordinates[i].values[1] = solution.X[i];
        }

        // Scale to [-1,1]^2 for numerical conditioning in later steps. The
        // same affine map is applied to both components, so the map is a
        // similarity and the conformality is preserved.
        let fmin = this.mPlaneCoordinates[0].values[0], fmax = fmin;
        for (let i = 0; i < numPositions; ++i) {
            const p = this.mPlaneCoordinates[i].values;
            if (p[0] < fmin) {
                fmin = p[0];
            } else if (p[0] > fmax) {
                fmax = p[0];
            }
            if (p[1] < fmin) {
                fmin = p[1];
            } else if (p[1] > fmax) {
                fmax = p[1];
            }
        }
        const halfRange = 0.5 * (fmax - fmin);
        const invHalfRange = 1 / halfRange;
        for (let i = 0; i < numPositions; ++i) {
            const p = this.mPlaneCoordinates[i].values;
            p[0] = -1 + invHalfRange * (p[0] - fmin);
            p[1] = -1 + invHalfRange * (p[1] - fmin);
        }

        // Map the plane coordinates to the sphere using inverse
        // stereographic projection. The main issue is selecting a
        // translation in (x,y) and a radius of the projection sphere. Both
        // factors strongly influence the final result.

        // Use the average as the south pole. The points tend to be clustered
        // approximately in the middle of the conformally mapped punctured
        // triangle, so the average is a good choice to place the pole.
        let originX = 0, originY = 0;
        for (let i = 0; i < numPositions; ++i) {
            originX += this.mPlaneCoordinates[i].values[0];
            originY += this.mPlaneCoordinates[i].values[1];
        }
        originX /= numPositions;
        originY /= numPositions;
        for (let i = 0; i < numPositions; ++i) {
            this.mPlaneCoordinates[i].values[0] -= originX;
            this.mPlaneCoordinates[i].values[1] -= originY;
        }

        this.mMinPlaneCoordinate = this.mPlaneCoordinates[0].clone();
        this.mMaxPlaneCoordinate = this.mPlaneCoordinates[0].clone();
        for (let i = 1; i < numPositions; ++i) {
            const p = this.mPlaneCoordinates[i].values;
            if (p[0] < this.mMinPlaneCoordinate.values[0]) {
                this.mMinPlaneCoordinate.values[0] = p[0];
            } else if (p[0] > this.mMaxPlaneCoordinate.values[0]) {
                this.mMaxPlaneCoordinate.values[0] = p[0];
            }

            if (p[1] < this.mMinPlaneCoordinate.values[1]) {
                this.mMinPlaneCoordinate.values[1] = p[1];
            } else if (p[1] > this.mMaxPlaneCoordinate.values[1]) {
                this.mMaxPlaneCoordinate.values[1] = p[1];
            }
        }

        // Select the radius of the sphere so that the projected punctured
        // triangle has an area whose fraction of the total spherical area is
        // the same fraction as the area of the punctured triangle to the
        // total area of the original triangle mesh.
        let twoTotalArea = 0;
        for (let t = 0; t < numTriangles; ++t) {
            const P0 = positions[indices[3 * t]];
            const P1 = positions[indices[3 * t + 1]];
            const P2 = positions[indices[3 * t + 2]];
            const E0 = sub(P1, P0);
            const E1 = sub(P2, P0);
            twoTotalArea += length(cross(E0, E1));
        }
        this.computeSphereRadius(v0, v1, v2, twoArea / twoTotalArea);
        const sqrSphereRadius = this.mSphereRadius * this.mSphereRadius;

        // Inverse stereographic projection to obtain sphere coordinates. The
        // sphere is centered at the origin and has radius 1.
        for (let i = 0; i < numPositions; ++i) {
            const p = this.mPlaneCoordinates[i];
            const rSqr = dot(p, p);
            const mult = 1 / (rSqr + sqrSphereRadius);
            const x = 2 * mult * sqrSphereRadius * p.values[0];
            const y = 2 * mult * sqrSphereRadius * p.values[1];
            const z = mult * this.mSphereRadius * (rSqr - sqrSphereRadius);
            const s = this.mSphereCoordinates[i].values;
            s[0] = x / this.mSphereRadius;
            s[1] = y / this.mSphereRadius;
            s[2] = z / this.mSphereRadius;
        }

        return converged;
    }

    // Conformal mapping of the mesh to a plane. The array of coordinates has
    // a one-to-one correspondence with the input vertex array.
    getPlaneCoordinates(): Vector[] {
        return this.mPlaneCoordinates;
    }

    getMinPlaneCoordinate(): Vector {
        return this.mMinPlaneCoordinate;
    }

    getMaxPlaneCoordinate(): Vector {
        return this.mMaxPlaneCoordinate;
    }

    // Conformal mapping of the mesh to a sphere (centered at the origin).
    // The array of coordinates has a one-to-one correspondence with the
    // input vertex array. The coordinates are on the unit sphere; the radius
    // reported by getSphereRadius() is the radius of the sphere used by the
    // inverse stereographic projection, not the radius of the output sphere.
    getSphereCoordinates(): Vector[] {
        return this.mSphereCoordinates;
    }

    getSphereRadius(): number {
        return this.mSphereRadius;
    }

    private computeSphereRadius(v0: number, v1: number, v2: number,
        areaFraction: number): void {
        const V0 = this.mPlaneCoordinates[v0];
        const V1 = this.mPlaneCoordinates[v1];
        const V2 = this.mPlaneCoordinates[v2];

        const r0Sqr = dot(V0, V0);
        const r1Sqr = dot(V1, V1);
        const r2Sqr = dot(V2, V2);
        const diffR10 = r1Sqr - r0Sqr;
        const diffR20 = r2Sqr - r0Sqr;
        const diffX10 = V1.values[0] - V0.values[0];
        const diffY10 = V1.values[1] - V0.values[1];
        const diffX20 = V2.values[0] - V0.values[0];
        const diffY20 = V2.values[1] - V0.values[1];
        const diffRX10 = V1.values[0] * r0Sqr - V0.values[0] * r1Sqr;
        const diffRY10 = V1.values[1] * r0Sqr - V0.values[1] * r1Sqr;
        const diffRX20 = V2.values[0] * r0Sqr - V0.values[0] * r2Sqr;
        const diffRY20 = V2.values[1] * r0Sqr - V0.values[1] * r2Sqr;

        const c0 = diffR20 * diffRY10 - diffR10 * diffRY20;
        const c1 = diffR20 * diffY10 - diffR10 * diffY20;
        const d0 = diffR10 * diffRX20 - diffR20 * diffRX10;
        const d1 = diffR10 * diffX20 - diffR20 * diffX10;
        const e0 = diffRX10 * diffRY20 - diffRX20 * diffRY10;
        const e1 = diffRX10 * diffY20 - diffRX20 * diffY10;
        const e2 = diffX10 * diffY20 - diffX20 * diffY10;

        const poly0 = new Polynomial1(6);
        poly0.set(0, 0);
        poly0.set(1, 0);
        poly0.set(2, e0 * e0);
        poly0.set(3, c0 * c0 + d0 * d0 + 2 * e0 * e1);
        poly0.set(4, 2 * (c0 * c1 + d0 * d1 + e0 * e1) + e1 * e1);
        poly0.set(5, c1 * c1 + d1 * d1 + 2 * e1 * e2);
        poly0.set(6, e2 * e2);

        const qpoly0 = Polynomial1.fromCoefficients([r0Sqr, 1]);
        const qpoly1 = Polynomial1.fromCoefficients([r1Sqr, 1]);
        const qpoly2 = Polynomial1.fromCoefficients([r2Sqr, 1]);

        const tmp = areaFraction * GTE_C_PI;
        const amp = tmp * tmp;

        let poly1 = qpoly0.mul(amp);
        poly1 = poly1.mul(qpoly0);
        poly1 = poly1.mul(qpoly0);
        poly1 = poly1.mul(qpoly0);
        poly1 = poly1.mul(qpoly1);
        poly1 = poly1.mul(qpoly1);
        poly1 = poly1.mul(qpoly2);
        poly1 = poly1.mul(qpoly2);

        const poly2 = poly1.sub(poly0);
        logAssert(poly2.getDegree() <= 8, 'Expecting degree no larger than 8.');

        // Bound a root near zero and apply bisection to find t.
        let tmin = 0, fmin = poly2.evaluate(tmin);
        let tmax = 1;
        const fmaxInitial = poly2.evaluate(tmax);
        logAssert(fmin > 0 && fmaxInitial < 0,
            'Expecting opposite-signed extremes.');

        // Determine the number of iterations to get 'digits' of accuracy.
        const digits = 6;
        const tmp0 = Math.log(tmax - tmin);
        const tmp1 = digits * GTE_C_LN_10;
        const arg = (tmp0 + tmp1) * GTE_C_INV_LN_2;
        const maxIterations = Math.trunc(arg + 0.5);
        let tmid = 0;
        for (let i = 0; i < maxIterations; ++i) {
            tmid = 0.5 * (tmin + tmax);
            const fmid = poly2.evaluate(tmid);
            const product = fmid * fmin;
            if (product < 0) {
                tmax = tmid;
            } else {
                tmin = tmid;
                fmin = fmid;
            }
        }

        this.mSphereRadius = Math.sqrt(tmid);
    }
}
