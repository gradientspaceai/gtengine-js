// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprCylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithm for least-squares fitting of a point set by a cylinder is
// described in
//   https://www.geometrictools.com/Documentation/CylinderFitting.pdf
// This document shows how to compute the cylinder radius r and the cylinder
// axis as a line C + h * W with origin C, unit-length direction W, and any
// real-valued h. The implementation here adds one additional step. It
// projects the point set onto the cylinder axis, computes the bounding
// h-interval [hmin, hmax] for the projections and sets the cylinder center
// to C + ((hmin + hmax) / 2) * W and the cylinder height to hmax - hmin.
//
// Port notes:
// * The four upstream constructor overloads become the named static
//   factories 'fromHemisphereSearch', 'fromEigenIndex' and
//   'fromCylinderAxis' (PORTING.md: ambiguous C++ constructor overloads
//   become static factories). As upstream, 'fromHemisphereSearch' with
//   fitPoints = false selects the fit-to-mesh mode.
// * The two 'operator()' overloads become 'compute' (fit to points, the
//   canonical query) and 'computeMesh' (fit to a triangle mesh).
// * Upstream runs the hemisphere search in std::thread workers when
//   numThreads > 0. The port has no threads; 'numThreads' is preserved
//   because it changes the (theta,phi) sample partition and the order of
//   the min-reduction, and therefore which sample wins a tie. The work is
//   performed sequentially in the same partition order, so the results are
//   deterministic and match a serialized execution of the upstream threads.

import { GTE_C_HALF_PI, GTE_C_TWO_PI } from './Constants';
import { ApprCircle2 } from './ApprCircle2';
import { Cylinder3 } from './Cylinder3';
import { Hypersphere } from './Hypersphere';
import { logAssert } from './Logger';
import {
    Matrix, addMatrix, divMatrix, mulMatrix, negateMatrix, outerProduct,
    subMatrix
} from './Matrix';
import { trace3x3 } from './Matrix3x3';
import { SymmetricEigensolver3x3 } from './SymmetricEigensolver3x3';
import { Vector, add, dot, mul, normalize, sub } from './Vector';
import { dotPerp } from './Vector2';
import { computeOrthogonalComplement3 } from './Vector3';

// The port of the private ApprCylinder3::ConstructorType enum class.
enum ConstructorType {
    FIT_BY_HEMISPHERE_SEARCH,
    FIT_USING_COVARIANCE_EIGENVECTOR,
    FIT_USING_SPECIFIED_AXIS,
    FIT_TO_MESH
}

// The result of the G(W) evaluation: the error, the projected center and
// the squared radius (upstream's output reference parameters).
interface GResult {
    gValue: number;
    PC: Vector;
    rsqr: number;
}

export class ApprCylinder3 {
    private mConstructorType: ConstructorType;

    // Parameters for the hemisphere-search constructor.
    private mNumThreads: number;
    private mNumThetaSamples: number;
    private mNumPhiSamples: number;

    // Parameter for the eigenvector-index constructor.
    private mEigenIndex: number;

    // Parameter for the specified-axis constructor.
    private mCylinderAxis: Vector;

    // A copy of the input points but translated by their average for
    // numerical robustness.
    private mX: Vector[];

    // Preprocessed information that depends only on the sample points. This
    // allows precomputed summations so that G(...) can be evaluated
    // extremely fast.
    private mMu: Vector;      // Vector<6>
    private mF0: Matrix;      // Matrix<3,3>
    private mF1: Matrix;      // Matrix<3,6>
    private mF2: Matrix;      // Matrix<6,6>

    private constructor(constructorType: ConstructorType, numThreads: number,
        numThetaSamples: number, numPhiSamples: number, eigenIndex: number,
        cylinderAxis: Vector) {
        this.mConstructorType = constructorType;
        this.mNumThreads = numThreads;
        this.mNumThetaSamples = numThetaSamples;
        this.mNumPhiSamples = numPhiSamples;
        this.mEigenIndex = eigenIndex;
        this.mCylinderAxis = cylinderAxis;
        this.mX = [];
        this.mMu = new Vector(6);
        this.mF0 = new Matrix(3, 3);
        this.mF1 = new Matrix(3, 6);
        this.mF2 = new Matrix(6, 6);
    }

    // Search the hemisphere for a minimum, choose numThetaSamples and
    // numPhiSamples to be positive (and preferably large). These are used to
    // generate a hemispherical grid of samples to be evaluated to find the
    // cylinder axis-direction W. Upstream runs the search multithreaded when
    // numThreads > 0; see the port notes at the top of this file.
    //
    // Set fitPoints to 'true' to use the algorithm described in the
    // aforementioned PDF file. Set fitPoints to 'false' if you want to fit a
    // cylinder to a triangle mesh; then call 'computeMesh'.
    static fromHemisphereSearch(numThreads: number, numThetaSamples: number,
        numPhiSamples: number, fitPoints: boolean = true): ApprCylinder3 {
        return new ApprCylinder3(
            fitPoints ? ConstructorType.FIT_BY_HEMISPHERE_SEARCH :
                ConstructorType.FIT_TO_MESH,
            numThreads, numThetaSamples, numPhiSamples, 0, new Vector(3));
    }

    // Choose one of the eigenvectors for the covariance matrix as the
    // cylinder axis direction. If eigenIndex is 0, the eigenvector
    // associated with the smallest eigenvalue is chosen. If eigenIndex is 2,
    // the eigenvector associated with the largest eigenvalue is chosen. If
    // eigenIndex is 1, the eigenvector associated with the median eigenvalue
    // is chosen; keep in mind that this could be the minimum or maximum
    // eigenvalue if the eigenspace has dimension 2 or 3.
    static fromEigenIndex(eigenIndex: number): ApprCylinder3 {
        return new ApprCylinder3(
            ConstructorType.FIT_USING_COVARIANCE_EIGENVECTOR,
            0, 0, 0, eigenIndex, new Vector(3));
    }

    // Choose the cylinder axis. If cylinderAxis is not the zero vector, the
    // factory will normalize it.
    static fromCylinderAxis(cylinderAxis: Vector): ApprCylinder3 {
        logAssert(cylinderAxis.size === 3,
            'ApprCylinder3: the cylinder axis must be 3D.');
        const axis = cylinderAxis.clone();
        normalize(axis, true);
        return new ApprCylinder3(ConstructorType.FIT_USING_SPECIFIED_AXIS,
            0, 0, 0, 0, axis);
    }

    // The algorithm must estimate 6 parameters, so the number of points must
    // be at least 6 but preferably larger. The returned value is the
    // root-mean-square of the least-squares error. The 'cylinder' is an
    // output parameter, matching the upstream reference parameter.
    compute(points: readonly Vector[], cylinder: Cylinder3): number {
        logAssert(this.mConstructorType !== ConstructorType.FIT_TO_MESH,
            'Call computeMesh(points, indices, cylinder) for fitting to a mesh.');

        const numPoints = points.length;
        logAssert(numPoints >= 6, 'Fitting requires at least 6 points.');
        for (const point of points) {
            logAssert(point.size === 3, 'ApprCylinder3: points must be 3D.');
        }

        this.mX = [];
        cylinder.axis.origin = new Vector(3);
        cylinder.axis.direction = new Vector(3);
        cylinder.radius = 0;
        cylinder.height = 0;

        const average = this.preprocess(points);

        // Fit the points based on which factory the caller used. The
        // direction is either estimated or selected directly or indirectly
        // by the caller. The center and squared radius are estimated.
        let minPC = new Vector(3);
        let minW = new Vector(3);
        let minRSqr = 0, minError = 0;

        if (this.mConstructorType === ConstructorType.FIT_BY_HEMISPHERE_SEARCH) {
            logAssert(this.mNumThetaSamples > 0 && this.mNumPhiSamples > 0,
                'The number of theta and psi samples must be positive.');

            // Search the hemisphere for the vector that leads to minimum
            // error and use it for the cylinder axis.
            const result = (this.mNumThreads === 0) ?
                this.computeSingleThreaded() : this.computeMultiThreaded();
            minError = result.error;
            minPC = result.PC;
            minW = result.W;
            minRSqr = result.rsqr;
        }
        else if (this.mConstructorType ===
            ConstructorType.FIT_USING_COVARIANCE_EIGENVECTOR) {
            logAssert(this.mEigenIndex < 3 && this.mEigenIndex >= 0,
                'Eigenvector index is out of range.');

            // Use the eigenvector corresponding to mEigenIndex of the
            // eigenvectors of the covariance matrix as the cylinder axis
            // direction. The eigenvectors are sorted from smallest
            // eigenvalue (mEigenIndex = 0) to largest eigenvalue
            // (mEigenIndex = 2).
            const result = this.computeUsingCovariance();
            minError = result.error;
            minPC = result.PC;
            minW = result.W;
            minRSqr = result.rsqr;
        }
        else {  // ConstructorType.FIT_USING_SPECIFIED_AXIS
            logAssert(!this.mCylinderAxis.equals(new Vector(3)),
                'The cylinder axis must be nonzero.');

            const result = this.computeUsingDirection();
            minError = result.error;
            minPC = result.PC;
            minW = result.W;
            minRSqr = result.rsqr;
        }

        // Translate back to the original space by the average of the points.
        cylinder.axis.origin = add(minPC, average);
        cylinder.axis.direction = minW;

        // Compute the cylinder radius.
        cylinder.radius = Math.sqrt(minRSqr);

        // Project the points onto the cylinder axis and choose the cylinder
        // center and cylinder height as described in the comments at the top
        // of this file.
        let hmin = 0, hmax = 0;
        for (let i = 0; i < numPoints; ++i) {
            const h = dot(cylinder.axis.direction,
                sub(points[i], cylinder.axis.origin));
            hmin = Math.min(h, hmin);
            hmax = Math.max(h, hmax);
        }

        const hmid = 0.5 * (hmin + hmax);
        cylinder.axis.origin = add(cylinder.axis.origin,
            mul(hmid, cylinder.axis.direction));
        cylinder.height = hmax - hmin;
        return minError;
    }

    // Use this function for fitting a cylinder to a mesh. For each candidate
    // cylinder axis direction on the hemisphere, project the triangles onto
    // a plane perpendicular to the axis. Compute the sum of 2*area for the
    // projected triangles. The direction that minimizes this measurement is
    // used as the cylinder axis direction. The projected points for this
    // direction are fit with a circle whose center is used to generate the
    // cylinder center and whose radius is used for the cylinder radius. The
    // projections of the points onto the axis are used to determine the
    // minimum and maximum coordinates along the axis, which are then used to
    // compute the cylinder height. The cylinder center is adjusted to be
    // midway between the minimum and maximum coordinates along the axis.
    //
    // The 'indices' array has 3 entries per triangle.
    computeMesh(points: readonly Vector[], indices: readonly number[],
        cylinder: Cylinder3): void {
        logAssert(this.mConstructorType === ConstructorType.FIT_TO_MESH,
            'Call compute(points, cylinder) for fitting to points.');

        const numPoints = points.length;
        const numTriangles = Math.floor(indices.length / 3);
        logAssert(numPoints >= 6 && numTriangles >= 2,
            'Fitting requires at least 6 points and 2 triangles.');
        for (const point of points) {
            logAssert(point.size === 3, 'ApprCylinder3: points must be 3D.');
        }
        logAssert(this.mNumThetaSamples > 0 && this.mNumPhiSamples > 0,
            'The number of theta and psi samples must be positive.');

        this.mX = [];
        cylinder.axis.origin = new Vector(3);
        cylinder.axis.direction = new Vector(3);
        cylinder.radius = 0;
        cylinder.height = 0;

        // Translate the points by translating their average to the origin.
        // This helps with numerical stability.
        let average = new Vector(3);
        for (let i = 0; i < numPoints; ++i) {
            average = add(average, points[i]);
        }
        average = mul(average, 1 / numPoints);
        this.mX = new Array<Vector>(numPoints);
        for (let i = 0; i < numPoints; ++i) {
            this.mX[i] = sub(points[i], average);
        }

        // Search the hemisphere for the vector that leads to minimum
        // measure and use it for the cylinder axis.
        if (this.mNumThreads === 0) {
            this.fitToMeshSingleThreaded(this.mX, numTriangles, indices,
                cylinder);
        }
        else {
            this.fitToMeshMultiThreaded(this.mX, numTriangles, indices,
                cylinder);
        }

        // Translate to the original coordinate system.
        cylinder.axis.origin = add(cylinder.axis.origin, average);
    }

    // Copy the points and translate by the average for numerical
    // robustness, then precompute the summations used by G(...). The
    // average is returned (upstream's output reference parameter).
    private preprocess(points: readonly Vector[]): Vector {
        const numPoints = points.length;
        const rNumPoints = numPoints;
        let average = new Vector(3);
        for (let i = 0; i < numPoints; ++i) {
            average = add(average, points[i]);
        }
        average = mul(average, 1 / rNumPoints);
        this.mX = new Array<Vector>(numPoints);
        for (let i = 0; i < numPoints; ++i) {
            this.mX[i] = sub(points[i], average);
        }

        const products: number[][] = new Array(numPoints);
        this.mMu = new Vector(6);
        for (let i = 0; i < numPoints; ++i) {
            const x = this.mX[i].values;
            const p = [
                x[0] * x[0],
                x[0] * x[1],
                x[0] * x[2],
                x[1] * x[1],
                x[1] * x[2],
                x[2] * x[2]
            ];
            products[i] = p;
            this.mMu.values[0] += p[0];
            this.mMu.values[1] += 2 * p[1];
            this.mMu.values[2] += 2 * p[2];
            this.mMu.values[3] += p[3];
            this.mMu.values[4] += 2 * p[4];
            this.mMu.values[5] += p[5];
        }
        this.mMu = mul(this.mMu, 1 / rNumPoints);

        this.mF0 = new Matrix(3, 3);
        this.mF1 = new Matrix(3, 6);
        this.mF2 = new Matrix(6, 6);
        for (let i = 0; i < numPoints; ++i) {
            const p = products[i];
            const delta = Vector.fromArray([
                p[0] - this.mMu.values[0],
                2 * p[1] - this.mMu.values[1],
                2 * p[2] - this.mMu.values[2],
                p[3] - this.mMu.values[3],
                2 * p[4] - this.mMu.values[4],
                p[5] - this.mMu.values[5]
            ]);
            this.mF0.set(0, 0, this.mF0.get(0, 0) + p[0]);
            this.mF0.set(0, 1, this.mF0.get(0, 1) + p[1]);
            this.mF0.set(0, 2, this.mF0.get(0, 2) + p[2]);
            this.mF0.set(1, 1, this.mF0.get(1, 1) + p[3]);
            this.mF0.set(1, 2, this.mF0.get(1, 2) + p[4]);
            this.mF0.set(2, 2, this.mF0.get(2, 2) + p[5]);
            this.mF1 = addMatrix(this.mF1, outerProduct(this.mX[i], delta));
            this.mF2 = addMatrix(this.mF2, outerProduct(delta, delta));
        }
        this.mF0 = divMatrix(this.mF0, rNumPoints);
        this.mF0.set(1, 0, this.mF0.get(0, 1));
        this.mF0.set(2, 0, this.mF0.get(0, 2));
        this.mF0.set(2, 1, this.mF0.get(1, 2));
        this.mF1 = divMatrix(this.mF1, rNumPoints);
        this.mF2 = divMatrix(this.mF2, rNumPoints);

        return average;
    }

    private computeUsingDirection():
        { error: number, PC: Vector, W: Vector, rsqr: number } {
        const W = this.mCylinderAxis.clone();
        const g = this.G(W);
        return { error: g.gValue, PC: g.PC, W, rsqr: g.rsqr };
    }

    private computeUsingCovariance():
        { error: number, PC: Vector, W: Vector, rsqr: number } {
        let covar = new Matrix(3, 3);  // zero matrix
        for (const X of this.mX) {
            covar = addMatrix(covar, outerProduct(X, X));
        }
        covar = divMatrix(covar, this.mX.length);
        const solver = new SymmetricEigensolver3x3();
        const { evecs } = solver.solve(covar.get(0, 0), covar.get(0, 1),
            covar.get(0, 2), covar.get(1, 1), covar.get(1, 2),
            covar.get(2, 2), true, +1);
        const W = Vector.fromArray(evecs[this.mEigenIndex]);
        const g = this.G(W);
        return { error: g.gValue, PC: g.PC, W, rsqr: g.rsqr };
    }

    private computeSingleThreaded():
        { error: number, PC: Vector, W: Vector, rsqr: number } {
        const iMultiplier = GTE_C_TWO_PI / this.mNumThetaSamples;
        const jMultiplier = GTE_C_HALF_PI / this.mNumPhiSamples;

        // Handle the north pole (0,0,1) separately.
        let minW = Vector.fromArray([0, 0, 1]);
        const northPole = this.G(minW);
        let minError = northPole.gValue;
        let minPC = northPole.PC;
        let minRSqr = northPole.rsqr;

        for (let j = 1; j <= this.mNumPhiSamples; ++j) {
            const phi = jMultiplier * j;  // in [0,pi/2]
            const csphi = Math.cos(phi);
            const snphi = Math.sin(phi);
            for (let i = 0; i < this.mNumThetaSamples; ++i) {
                const theta = iMultiplier * i;  // in [0,2*pi)
                const cstheta = Math.cos(theta);
                const sntheta = Math.sin(theta);
                const W = Vector.fromArray([
                    cstheta * snphi, sntheta * snphi, csphi]);
                const g = this.G(W);
                if (g.gValue < minError) {
                    minError = g.gValue;
                    minRSqr = g.rsqr;
                    minW = W;
                    minPC = g.PC;
                }
            }
        }

        return { error: minError, PC: minPC, W: minW, rsqr: minRSqr };
    }

    // The port of ComputeMultiThreaded. The per-thread partition of the phi
    // samples and the order of the final min-reduction are preserved so that
    // the result matches a serialized execution of the upstream threads.
    private computeMultiThreaded():
        { error: number, PC: Vector, W: Vector, rsqr: number } {
        const iMultiplier = GTE_C_TWO_PI / this.mNumThetaSamples;
        const jMultiplier = GTE_C_HALF_PI / this.mNumPhiSamples;

        // Handle the north pole (0,0,1) separately.
        let minW = Vector.fromArray([0, 0, 1]);
        const northPole = this.G(minW);
        let minError = northPole.gValue;
        let minPC = northPole.PC;
        let minRSqr = northPole.rsqr;

        interface Local {
            error: number;
            rsqr: number;
            W: Vector;
            PC: Vector;
            jmin: number;
            jmax: number;
        }

        const numThreads = this.mNumThreads;
        const local: Local[] = new Array(numThreads);
        const numPhiSamplesPerThread =
            Math.floor(this.mNumPhiSamples / numThreads);
        for (let t = 0; t < numThreads; ++t) {
            local[t] = {
                error: Number.MAX_VALUE,
                rsqr: 0,
                W: new Vector(3),
                PC: new Vector(3),
                jmin: numPhiSamplesPerThread * t,
                jmax: numPhiSamplesPerThread * (t + 1)
            };
        }
        local[numThreads - 1].jmax = this.mNumPhiSamples + 1;

        for (let t = 0; t < numThreads; ++t) {
            for (let j = local[t].jmin; j < local[t].jmax; ++j) {
                // phi in [0,pi/2]
                const phi = jMultiplier * j;
                const csphi = Math.cos(phi);
                const snphi = Math.sin(phi);
                for (let i = 0; i < this.mNumThetaSamples; ++i) {
                    // theta in [0,2*pi)
                    const theta = iMultiplier * i;
                    const cstheta = Math.cos(theta);
                    const sntheta = Math.sin(theta);
                    const W = Vector.fromArray([
                        cstheta * snphi, sntheta * snphi, csphi]);
                    const g = this.G(W);
                    if (g.gValue < local[t].error) {
                        local[t].error = g.gValue;
                        local[t].rsqr = g.rsqr;
                        local[t].W = W;
                        local[t].PC = g.PC;
                    }
                }
            }
        }

        for (let t = 0; t < numThreads; ++t) {
            if (local[t].error < minError) {
                minError = local[t].error;
                minRSqr = local[t].rsqr;
                minW = local[t].W;
                minPC = local[t].PC;
            }
        }

        return { error: minError, PC: minPC, W: minW, rsqr: minRSqr };
    }

    private G(W: Vector): GResult {
        const P = subMatrix(Matrix.identity(3, 3), outerProduct(W, W));
        const w = W.values;
        const S = Matrix.fromArray(3, 3, [
            0, -w[2], w[1],
            w[2], 0, -w[0],
            -w[1], w[0], 0
        ]);

        const A = mulMatrix(mulMatrix(P, this.mF0), P);
        const hatA = negateMatrix(mulMatrix(mulMatrix(S, A), S));
        const hatAA = mulMatrix(hatA, A);
        const trace = trace3x3(hatAA);
        const Q = divMatrix(hatA, trace);
        const pVec = Vector.fromArray([
            P.get(0, 0), P.get(0, 1), P.get(0, 2),
            P.get(1, 1), P.get(1, 2), P.get(2, 2)]);
        const alpha = mulMatrix(this.mF1, pVec);
        const beta = mulMatrix(Q, alpha);
        const term0 = dot(pVec, mulMatrix(this.mF2, pVec));
        const term1 = 4 * dot(alpha, beta);
        const term2 = 4 * dot(beta, mulMatrix(this.mF0, beta));
        const gValue = (term0 - term1 + term2) / this.mX.length;

        const PC = beta;
        const rsqr = dot(pVec, this.mMu) + dot(PC, PC);
        return { gValue, PC, rsqr };
    }

    private fitToMeshSingleThreaded(points: readonly Vector[],
        numTriangles: number, indices: readonly number[],
        cylinder: Cylinder3): void {
        const iMultiplier = GTE_C_TWO_PI / this.mNumThetaSamples;
        const jMultiplier = GTE_C_HALF_PI / this.mNumPhiSamples;

        // Handle the north pole (0,0,1) separately.
        let minDirection = Vector.fromArray([0, 0, 1]);
        let minMeasure = ApprCylinder3.getProjectionMeasure(minDirection,
            points, numTriangles, indices);

        // Process a regular grid of (theta,phi) angles.
        for (let j = 1; j <= this.mNumPhiSamples; ++j) {
            const phi = jMultiplier * j;  // in [0,pi/2]
            const csphi = Math.cos(phi);
            const snphi = Math.sin(phi);
            for (let i = 0; i < this.mNumThetaSamples; ++i) {
                const theta = iMultiplier * i;  // in [0,2*pi)
                const cstheta = Math.cos(theta);
                const sntheta = Math.sin(theta);
                const direction = Vector.fromArray([
                    cstheta * snphi, sntheta * snphi, csphi]);

                const measure = ApprCylinder3.getProjectionMeasure(direction,
                    points, numTriangles, indices);
                if (measure < minMeasure) {
                    minDirection = direction;
                    minMeasure = measure;
                }
            }
        }

        ApprCylinder3.finishCylinder(minDirection, points, cylinder);
    }

    // The port of FitToMeshMultiThreaded; see the note on
    // computeMultiThreaded.
    private fitToMeshMultiThreaded(points: readonly Vector[],
        numTriangles: number, indices: readonly number[],
        cylinder: Cylinder3): void {
        const iMultiplier = GTE_C_TWO_PI / this.mNumThetaSamples;
        const jMultiplier = GTE_C_HALF_PI / this.mNumPhiSamples;

        // Handle the north pole (0,0,1) separately.
        let minDirection = Vector.fromArray([0, 0, 1]);
        let minMeasure = ApprCylinder3.getProjectionMeasure(minDirection,
            points, numTriangles, indices);

        interface Local {
            direction: Vector;
            measure: number;
            jmin: number;
            jmax: number;
        }

        const numThreads = this.mNumThreads;
        const local: Local[] = new Array(numThreads);
        const numPhiSamplesPerThread =
            Math.floor(this.mNumPhiSamples / numThreads);
        for (let t = 0; t < numThreads; ++t) {
            local[t] = {
                direction: new Vector(3),
                measure: Number.MAX_VALUE,
                jmin: numPhiSamplesPerThread * t,
                jmax: numPhiSamplesPerThread * (t + 1)
            };
        }
        local[numThreads - 1].jmax = this.mNumPhiSamples + 1;

        for (let t = 0; t < numThreads; ++t) {
            for (let j = local[t].jmin; j < local[t].jmax; ++j) {
                // phi in [0,pi/2]
                const phi = jMultiplier * j;
                const csphi = Math.cos(phi);
                const snphi = Math.sin(phi);
                for (let i = 0; i < this.mNumThetaSamples; ++i) {
                    // theta in [0,2*pi)
                    const theta = iMultiplier * i;
                    const cstheta = Math.cos(theta);
                    const sntheta = Math.sin(theta);
                    const direction = Vector.fromArray([
                        cstheta * snphi, sntheta * snphi, csphi]);

                    const measure = ApprCylinder3.getProjectionMeasure(
                        direction, points, numTriangles, indices);
                    if (measure < local[t].measure) {
                        local[t].direction = direction;
                        local[t].measure = measure;
                    }
                }
            }
        }

        for (let t = 0; t < numThreads; ++t) {
            if (local[t].measure < minMeasure) {
                minMeasure = local[t].measure;
                minDirection = local[t].direction;
            }
        }

        ApprCylinder3.finishCylinder(minDirection, points, cylinder);
    }

    private static getProjectionMeasure(direction: Vector,
        points: readonly Vector[], numTriangles: number,
        indices: readonly number[]): number {
        // NOTE: 'direction' is cloned because computeOrthogonalComplement3
        // normalizes basis[0] in place; C++ copies the vector into the basis
        // array, TypeScript would alias it.
        const basis: Vector[] = [direction.clone(), new Vector(3),
            new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const numPoints = points.length;
        const projections: Vector[] = new Array(numPoints);
        for (let i = 0; i < numPoints; ++i) {
            projections[i] = Vector.fromArray([
                dot(basis[1], points[i]), dot(basis[2], points[i])]);
        }

        // Add up 2*area of the triangles.
        let measure = 0;
        let k = 0;
        for (let t = 0; t < numTriangles; ++t) {
            const V0 = projections[indices[k++]];
            const V1 = projections[indices[k++]];
            const V2 = projections[indices[k++]];
            const edge10 = sub(V1, V0);
            const edge20 = sub(V2, V0);
            measure += Math.abs(dotPerp(edge10, edge20));
        }
        return measure;
    }

    private static finishCylinder(minDirection: Vector,
        points: readonly Vector[], cylinder: Cylinder3): void {
        const basis: Vector[] = [minDirection.clone(), new Vector(3),
            new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const numPoints = points.length;
        const projections: Vector[] = new Array(numPoints);
        // NOTE (upstream): hmax is initialized to 0 rather than to
        // -max(). Because the points passed to this function have been
        // translated so that their average is the origin, the projections
        // straddle 0 and the initialization is harmless. The port preserves
        // the upstream behavior.
        let hmin = Number.MAX_VALUE;
        let hmax = 0;
        for (let i = 0; i < numPoints; ++i) {
            const h = dot(basis[0], points[i]);
            if (h < hmin) {
                hmin = h;
            }
            if (h > hmax) {
                hmax = h;
            }

            projections[i] = Vector.fromArray([
                dot(basis[1], points[i]), dot(basis[2], points[i])]);
        }
        const fitter = new ApprCircle2();
        const circle = new Hypersphere(2);
        fitter.fitUsingSquaredLengths(projections, circle);

        const minCenter = add(mul(circle.center.values[0], basis[1]),
            mul(circle.center.values[1], basis[2]));
        cylinder.axis.origin = add(minCenter,
            mul(0.5 * (hmax + hmin), minDirection));
        cylinder.axis.direction = minDirection;
        cylinder.radius = circle.radius;
        cylinder.height = hmax - hmin;
    }
}

