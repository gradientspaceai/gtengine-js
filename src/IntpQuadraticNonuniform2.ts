// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpQuadraticNonuniform2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Quadratic interpolation of a network of triangles whose vertices are of
// the form (x,y,f(x,y)). This code is an implementation of the algorithm
// found in
//
//   Zoltan J. Cendes and Steven H. Wong,
//   C1 quadratic interpolation over arbitrary point sets,
//   IEEE Computer Graphics & Applications,
//   pp. 8-16, 1987
//
// Port notes:
// * Upstream the mesh is a template parameter constrained only by a
//   duck-typed interface, so the port declares that interface explicitly as
//   IntpQuadraticNonuniform2TriangleMesh (the precedent set by
//   IntpLinearNonuniform2 in B47). Both PlanarMesh and Delaunay2Mesh
//   satisfy it. The C++ per-triangle accessors that return 'bool' and write
//   to a reference parameter,
//     bool GetVertices(int32_t, std::array<Vector2<T>, 3>&)
//     bool GetIndices(int32_t, std::array<int32_t, 3>&)
//     bool GetAdjacencies(int32_t, std::array<int32_t, 3>&)
//     bool GetBarycentrics(int32_t, Vector2<T> const&, std::array<T, 3>&)
//   become getTriangleVertices(t), getTriangleIndices(t),
//   getTriangleAdjacencies(t) and getBarycentrics(t, P), which return the
//   triple or null. The naming matches PlanarMesh and Delaunay2Mesh, whose
//   unindexed GetVertices()/GetIndices()/GetAdjacencies() overloads keep the
//   short names.
// * GetInvalidIndex() is optional on the mesh interface; when it is absent
//   the sentinel is -1 (PlanarMesh reports -1, Delaunay2Mesh reports
//   getInvalidIndex() === -1).
// * The two constructors become the static factories fromSpatialDelta() and
//   fromDerivatives(), following the ambiguous-overload precedent.
// * The C++ operator()(P, F&, FX&, FY&) returning bool becomes evaluate(P),
//   which returns { valid, F, FX, FY }; the values are meaningless when
//   valid is false.

import { AlignedBox } from './AlignedBox';
import { inscribeCircle2 } from './ContScribeCircle2';
import { DistPointAlignedBox } from './DistPointAlignedBox';
import { logAssert } from './Logger';
import { Vector } from './Vector';
import { computeBarycentrics2 } from './Vector2';

// The duck-typed triangle mesh required by the interpolator.
export interface IntpQuadraticNonuniform2TriangleMesh {
    getNumVertices(): number;
    getNumTriangles(): number;

    // The array of mesh vertices, each a 2D Vector.
    getVertices(): readonly Vector[];

    // The flat array of 3 * getNumTriangles() vertex indices.
    getIndices(): readonly number[];

    // The vertices of triangle t, or null when t is invalid.
    getTriangleVertices(t: number): readonly Vector[] | null;

    // The vertex indices of triangle t, or null when t is invalid.
    getTriangleIndices(t: number): readonly number[] | null;

    // The adjacent triangles of triangle t, or null when t is invalid. A
    // negative adjacency indicates a boundary edge.
    getTriangleAdjacencies(t: number): readonly number[] | null;

    // The barycentric coordinates of P with respect to triangle t, or null
    // when the triangle is degenerate.
    getBarycentrics(t: number, p: Vector): readonly number[] | null;

    // The index of the triangle containing P, or the invalid index when P is
    // outside the mesh.
    getContainingTriangle(p: Vector): number;

    // The "no such triangle" sentinel. When absent, -1 is used.
    getInvalidIndex?(): number;
}

export interface IntpQuadraticNonuniform2Result {
    // Valid is true if and only if the input point P is in the convex hull
    // of the input vertices, in which case the interpolation is valid.
    valid: boolean;

    // The interpolated function value and its first-order partial
    // derivatives; meaningful only when valid is true.
    F: number;
    FX: number;
    FY: number;
}

// The per-triangle data computed by the preprocessing pass.
class TriangleData {
    // The center of the inscribed circle of the triangle.
    center: Vector;

    // The cross-edge intersection points, one per edge.
    intersect: [Vector, Vector, Vector];

    // The 19 Bezier control values of the piecewise quadratic.
    coeff: number[];

    constructor() {
        this.center = Vector.zero(2);
        this.intersect = [Vector.zero(2), Vector.zero(2), Vector.zero(2)];
        this.coeff = new Array<number>(19).fill(0);
    }
}

// The port of the private Jet class: a function value and its first-order
// partial derivatives at a sample point.
interface Jet {
    F: number;
    FX: number;
    FY: number;
}

export class IntpQuadraticNonuniform2 {
    private mMesh: IntpQuadraticNonuniform2TriangleMesh;
    private mF: readonly number[];
    private mFX: readonly number[];
    private mFY: readonly number[];
    private mTData: TriangleData[];

    private constructor(mesh: IntpQuadraticNonuniform2TriangleMesh,
        F: readonly number[]) {
        logAssert(F.length >= mesh.getNumVertices(),
            'IntpQuadraticNonuniform2: too few function samples.');
        this.mMesh = mesh;
        this.mF = F;
        this.mFX = [];
        this.mFY = [];
        this.mTData = [];
    }

    // The port of the first constructor. It requires only F and a measure of
    // the rate of change of the function values relative to changes in the
    // spatial variables. The df/dx and df/dy values are estimated at the
    // sample points using mesh normals and spatialDelta.
    static fromSpatialDelta(mesh: IntpQuadraticNonuniform2TriangleMesh,
        F: readonly number[], spatialDelta: number): IntpQuadraticNonuniform2 {
        const intp = new IntpQuadraticNonuniform2(mesh, F);
        intp.estimateDerivatives(spatialDelta);
        intp.processTriangles();
        return intp;
    }

    // The port of the second constructor. It requires you to specify the
    // function values F and the first-order partial derivative values df/dx
    // and df/dy.
    static fromDerivatives(mesh: IntpQuadraticNonuniform2TriangleMesh,
        F: readonly number[], FX: readonly number[], FY: readonly number[]):
        IntpQuadraticNonuniform2 {
        const intp = new IntpQuadraticNonuniform2(mesh, F);
        logAssert(FX.length >= mesh.getNumVertices()
            && FY.length >= mesh.getNumVertices(),
            'IntpQuadraticNonuniform2: too few derivative samples.');
        intp.mFX = FX;
        intp.mFY = FY;
        intp.processTriangles();
        return intp;
    }

    // Quadratic interpolation. The 'valid' field is true if and only if the
    // input point is in the convex hull of the input vertices, in which case
    // the interpolation is valid.
    evaluate(P: Vector): IntpQuadraticNonuniform2Result {
        logAssert(P.size === 2, 'IntpQuadraticNonuniform2: the point must be 2D.');

        const invalidIndex = this.mMesh.getInvalidIndex
            ? this.mMesh.getInvalidIndex() : -1;
        const t = this.mMesh.getContainingTriangle(P);
        if (t === invalidIndex || t < 0 || t >= this.mTData.length) {
            // The point is outside the triangulation.
            return { valid: false, F: 0, FX: 0, FY: 0 };
        }

        // Get the vertices of the triangle.
        const V = this.mMesh.getTriangleVertices(t);
        if (V === null) {
            return { valid: false, F: 0, FX: 0, FY: 0 };
        }

        // Get the additional information for the triangle.
        const tData = this.mTData[t];

        // Determine which of the six subtriangles contains the target point.
        // Theoretically, P must be in one of these subtriangles.
        let sub0 = tData.center;
        let sub1 = Vector.zero(2);
        let sub2 = tData.intersect[2];
        let bary: [number, number, number] = [0, 0, 0];
        let index = 0;

        const barybox = AlignedBox.fromMinMax(Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 1, 1]));
        const pbQuery = new DistPointAlignedBox();
        let minIndex = 0;
        let minDistance = -1;
        let minBary: [number, number, number] = [0, 0, 0];
        let minSub0 = Vector.zero(2);
        let minSub1 = Vector.zero(2);
        let minSub2 = Vector.zero(2);

        for (index = 1; index <= 6; ++index) {
            sub1 = sub2;
            if ((index % 2) !== 0) {
                // index is odd
                sub2 = V[Math.floor(index / 2)];
            } else {
                // index is even
                sub2 = tData.intersect[Math.floor((index - 1) / 2)];
            }

            const local = computeBarycentrics2(P, sub0, sub1, sub2);
            bary = local.bary;
            if (local.valid
                && 0 <= bary[0] && bary[0] <= 1
                && 0 <= bary[1] && bary[1] <= 1
                && 0 <= bary[2] && bary[2] <= 1) {
                // P is in the triangle <sub0,sub1,sub2>.
                break;
            }

            // When computing with floating-point arithmetic, rounding errors
            // can cause us to reach this code when, theoretically, the point
            // is in the subtriangle. Keep track of the (b0,b1,b2) that is
            // closest to the barycentric cube [0,1]^3 and choose the triangle
            // corresponding to it when all 6 tests previously fail.
            const distance =
                pbQuery.compute(Vector.fromArray(bary), barybox).distance;
            if (minIndex === 0 || distance < minDistance) {
                minDistance = distance;
                minIndex = index;
                minBary = [bary[0], bary[1], bary[2]];
                minSub0 = sub0;
                minSub1 = sub1;
                minSub2 = sub2;
            }
        }

        // If the subtriangle was not found, rounding errors caused problems.
        // Choose the barycentric point closest to the box.
        if (index > 6) {
            index = minIndex;
            bary = minBary;
            sub0 = minSub0;
            sub1 = minSub1;
            sub2 = minSub2;
        }

        // Fetch Bezier control points.
        const bez: number[] = [
            tData.coeff[0],
            tData.coeff[12 + index],
            tData.coeff[13 + (index % 6)],
            tData.coeff[index],
            tData.coeff[6 + index],
            tData.coeff[1 + (index % 6)]
        ];

        // Evaluate the Bezier quadratic.
        const F = bary[0] * (bez[0] * bary[0] + bez[1] * bary[1] + bez[2] * bary[2])
            + bary[1] * (bez[1] * bary[0] + bez[3] * bary[1] + bez[4] * bary[2])
            + bary[2] * (bez[2] * bary[0] + bez[4] * bary[1] + bez[5] * bary[2]);

        // Evaluate the barycentric derivatives of F.
        const FU = 2 * (bez[0] * bary[0] + bez[1] * bary[1] + bez[2] * bary[2]);
        const FV = 2 * (bez[1] * bary[0] + bez[3] * bary[1] + bez[4] * bary[2]);
        const FW = 2 * (bez[2] * bary[0] + bez[4] * bary[1] + bez[5] * bary[2]);
        const duw = FU - FW;
        const dvw = FV - FW;

        // Convert back to (x,y) coordinates.
        const m00 = sub0.values[0] - sub2.values[0];
        const m10 = sub0.values[1] - sub2.values[1];
        const m01 = sub1.values[0] - sub2.values[0];
        const m11 = sub1.values[1] - sub2.values[1];
        const inv = 1 / (m00 * m11 - m10 * m01);

        return {
            valid: true,
            F,
            FX: inv * (m11 * duw - m10 * dvw),
            FY: inv * (m00 * dvw - m01 * duw)
        };
    }

    private estimateDerivatives(spatialDelta: number): void {
        const numVertices = this.mMesh.getNumVertices();
        const vertices = this.mMesh.getVertices();
        const numTriangles = this.mMesh.getNumTriangles();
        const indices = this.mMesh.getIndices();

        const FX = new Array<number>(numVertices).fill(0);
        const FY = new Array<number>(numVertices).fill(0);
        const FZ = new Array<number>(numVertices).fill(0);

        // Accumulate normals at spatial locations (averaging process).
        let i = 0;
        for (let t = 0; t < numTriangles; ++t) {
            // Get three vertices of the triangle.
            const v0 = indices[i++];
            const v1 = indices[i++];
            const v2 = indices[i++];

            // Compute the normal vector of the triangle (with positive
            // z-component).
            const p0 = vertices[v0].values;
            const p1 = vertices[v1].values;
            const p2 = vertices[v2].values;
            const dx1 = p1[0] - p0[0];
            const dy1 = p1[1] - p0[1];
            const dz1 = this.mF[v1] - this.mF[v0];
            const dx2 = p2[0] - p0[0];
            const dy2 = p2[1] - p0[1];
            const dz2 = this.mF[v2] - this.mF[v0];
            let nx = dy1 * dz2 - dy2 * dz1;
            let ny = dz1 * dx2 - dz2 * dx1;
            let nz = dx1 * dy2 - dx2 * dy1;
            if (nz < 0) {
                nx = -nx;
                ny = -ny;
                nz = -nz;
            }

            FX[v0] += nx; FY[v0] += ny; FZ[v0] += nz;
            FX[v1] += nx; FY[v1] += ny; FZ[v1] += nz;
            FX[v2] += nx; FY[v2] += ny; FZ[v2] += nz;
        }

        // Scale the normals to form (x,y,-1).
        for (let j = 0; j < numVertices; ++j) {
            if (FZ[j] !== 0) {
                const inv = -spatialDelta / FZ[j];
                FX[j] *= inv;
                FY[j] *= inv;
            } else {
                FX[j] = 0;
                FY[j] = 0;
            }
        }

        this.mFX = FX;
        this.mFY = FY;
    }

    private processTriangles(): void {
        // Add degenerate triangles to boundary triangles so that
        // interpolation at the boundary can be treated in the same way as
        // interpolation in the interior.

        // Compute the centers of the inscribed circles for the triangles.
        const vertices = this.mMesh.getVertices();
        const numTriangles = this.mMesh.getNumTriangles();
        const indices = this.mMesh.getIndices();
        this.mTData = [];
        for (let t = 0; t < numTriangles; ++t) {
            this.mTData.push(new TriangleData());
        }
        let i = 0;
        for (let t = 0; t < numTriangles; ++t) {
            const v0 = indices[i++];
            const v1 = indices[i++];
            const v2 = indices[i++];
            const circle =
                inscribeCircle2(vertices[v0], vertices[v1], vertices[v2]);
            if (circle !== null) {
                this.mTData[t].center = circle.center.clone();
            }
            // Upstream passes a value-initialized Circle2 to Inscribe and
            // ignores the returned 'bool', so a degenerate triangle leaves
            // the center at (0,0). The port keeps that behavior; a Delaunay
            // triangulation has no degenerate triangles.
        }

        // Compute the cross-edge intersections.
        for (let t = 0; t < numTriangles; ++t) {
            this.computeCrossEdgeIntersections(t);
        }

        // Compute the Bezier coefficients.
        for (let t = 0; t < numTriangles; ++t) {
            this.computeCoefficients(t);
        }
    }

    private computeCrossEdgeIntersections(t: number): void {
        // Get the vertices of the triangle.
        const V = this.mMesh.getTriangleVertices(t);
        if (V === null) {
            return;
        }

        // Get the centers of the adjacent triangles.
        const tData = this.mTData[t];
        const adjacencies = this.mMesh.getTriangleAdjacencies(t);
        if (adjacencies === null) {
            return;
        }

        for (let j0 = 2, j1 = 0; j1 < 3; j0 = j1++) {
            const a = adjacencies[j0];
            const V0 = V[j0].values;
            const V1 = V[j1].values;
            if (a >= 0) {
                // Get the center of the adjacent triangle's inscribed circle.
                const U = this.mTData[a].center.values;
                const C = tData.center.values;
                const m00 = V0[1] - V1[1];
                const m01 = V1[0] - V0[0];
                const m10 = C[1] - U[1];
                const m11 = U[0] - C[0];
                const r0 = m00 * V0[0] + m01 * V0[1];
                const r1 = m10 * C[0] + m11 * C[1];
                const invDet = 1 / (m00 * m11 - m01 * m10);
                tData.intersect[j0] = Vector.fromArray([
                    (m11 * r0 - m01 * r1) * invDet,
                    (m00 * r1 - m10 * r0) * invDet
                ]);
            } else {
                // No adjacent triangle, use the center of the edge.
                tData.intersect[j0] = Vector.fromArray([
                    0.5 * (V0[0] + V1[0]),
                    0.5 * (V0[1] + V1[1])
                ]);
            }
        }
    }

    private computeCoefficients(t: number): void {
        // Get the vertices of the triangle.
        const Vt = this.mMesh.getTriangleVertices(t);
        if (Vt === null) {
            return;
        }
        const V = [Vt[0].values, Vt[1].values, Vt[2].values];

        // Get the additional information for the triangle.
        const tData = this.mTData[t];

        // Get the sample data at the main triangle vertices.
        const indices = this.mMesh.getTriangleIndices(t);
        if (indices === null) {
            return;
        }
        const jet: Jet[] = [];
        for (let j = 0; j < 3; ++j) {
            const k = indices[j];
            jet.push({ F: this.mF[k], FX: this.mFX[k], FY: this.mFY[k] });
        }

        // Get the centers of the adjacent triangles.
        const adjacencies = this.mMesh.getTriangleAdjacencies(t);
        if (adjacencies === null) {
            return;
        }
        const U: Vector[] = [Vector.zero(2), Vector.zero(2), Vector.zero(2)];
        for (let j0 = 2, j1 = 0; j1 < 3; j0 = j1++) {
            const a = adjacencies[j0];
            if (a >= 0) {
                // Get the center of the adjacent triangle's inscribed circle.
                U[j0] = this.mTData[a].center;
            } else {
                // No adjacent triangle, use the center of the edge.
                U[j0] = Vector.fromArray([
                    0.5 * (V[j0][0] + V[j1][0]),
                    0.5 * (V[j0][1] + V[j1][1])
                ]);
            }
        }

        // Compute the intermediate terms.
        const cenT = this.mMesh.getBarycentrics(t, tData.center);
        const cen0 = this.mMesh.getBarycentrics(t, U[0]);
        const cen1 = this.mMesh.getBarycentrics(t, U[1]);
        const cen2 = this.mMesh.getBarycentrics(t, U[2]);
        if (cenT === null || cen0 === null || cen1 === null || cen2 === null) {
            return;
        }

        const alpha = (cenT[1] * cen1[0] - cenT[0] * cen1[1]) / (cen1[0] - cenT[0]);
        const beta = (cenT[2] * cen2[1] - cenT[1] * cen2[2]) / (cen2[1] - cenT[1]);
        const gamma = (cenT[0] * cen0[2] - cenT[2] * cen0[0]) / (cen0[2] - cenT[2]);
        const oneMinusAlpha = 1 - alpha;
        const oneMinusBeta = 1 - beta;
        const oneMinusGamma = 1 - gamma;

        const A = new Array<number>(9).fill(0);
        const B = new Array<number>(9).fill(0);

        let tmp = cenT[0] * V[0][0] + cenT[1] * V[1][0] + cenT[2] * V[2][0];
        A[0] = 0.5 * (tmp - V[0][0]);
        A[1] = 0.5 * (tmp - V[1][0]);
        A[2] = 0.5 * (tmp - V[2][0]);
        A[3] = 0.5 * beta * (V[2][0] - V[0][0]);
        A[4] = 0.5 * oneMinusGamma * (V[1][0] - V[0][0]);
        A[5] = 0.5 * gamma * (V[0][0] - V[1][0]);
        A[6] = 0.5 * oneMinusAlpha * (V[2][0] - V[1][0]);
        A[7] = 0.5 * alpha * (V[1][0] - V[2][0]);
        A[8] = 0.5 * oneMinusBeta * (V[0][0] - V[2][0]);

        tmp = cenT[0] * V[0][1] + cenT[1] * V[1][1] + cenT[2] * V[2][1];
        B[0] = 0.5 * (tmp - V[0][1]);
        B[1] = 0.5 * (tmp - V[1][1]);
        B[2] = 0.5 * (tmp - V[2][1]);
        B[3] = 0.5 * beta * (V[2][1] - V[0][1]);
        B[4] = 0.5 * oneMinusGamma * (V[1][1] - V[0][1]);
        B[5] = 0.5 * gamma * (V[0][1] - V[1][1]);
        B[6] = 0.5 * oneMinusAlpha * (V[2][1] - V[1][1]);
        B[7] = 0.5 * alpha * (V[1][1] - V[2][1]);
        B[8] = 0.5 * oneMinusBeta * (V[0][1] - V[2][1]);

        // Compute the Bezier coefficients.
        const coeff = tData.coeff;
        coeff[2] = jet[0].F;
        coeff[4] = jet[1].F;
        coeff[6] = jet[2].F;

        coeff[14] = jet[0].F + A[0] * jet[0].FX + B[0] * jet[0].FY;
        coeff[7] = jet[0].F + A[3] * jet[0].FX + B[3] * jet[0].FY;
        coeff[8] = jet[0].F + A[4] * jet[0].FX + B[4] * jet[0].FY;
        coeff[16] = jet[1].F + A[1] * jet[1].FX + B[1] * jet[1].FY;
        coeff[9] = jet[1].F + A[5] * jet[1].FX + B[5] * jet[1].FY;
        coeff[10] = jet[1].F + A[6] * jet[1].FX + B[6] * jet[1].FY;
        coeff[18] = jet[2].F + A[2] * jet[2].FX + B[2] * jet[2].FY;
        coeff[11] = jet[2].F + A[7] * jet[2].FX + B[7] * jet[2].FY;
        coeff[12] = jet[2].F + A[8] * jet[2].FX + B[8] * jet[2].FY;

        coeff[5] = alpha * coeff[10] + oneMinusAlpha * coeff[11];
        coeff[17] = alpha * coeff[16] + oneMinusAlpha * coeff[18];
        coeff[1] = beta * coeff[12] + oneMinusBeta * coeff[7];
        coeff[13] = beta * coeff[18] + oneMinusBeta * coeff[14];
        coeff[3] = gamma * coeff[8] + oneMinusGamma * coeff[9];
        coeff[15] = gamma * coeff[14] + oneMinusGamma * coeff[16];
        coeff[0] = cenT[0] * coeff[14] + cenT[1] * coeff[16] + cenT[2] * coeff[18];
    }
}
