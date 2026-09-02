// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RiemannianGeodesic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Computing geodesics on a surface is a differential geometric topic that
// involves Riemannian geometry. The algorithm for constructing geodesics that
// is implemented here uses a multiresolution approach. A description of the
// algorithm is in the document
// https://www.geometrictools.com/Documentation/RiemannianGeodesics.pdf
//
// Port notes:
// - The class is abstract: derived classes must implement computeMetric and
//   computeChristoffel1 (upstream pure virtual functions).
// - Upstream 'GVector<Real>' points become 'GVector'; the metric tensor and
//   the Christoffel symbols are 'GMatrix' objects.
// - 'ComputeGeodesic' writes into 'int32_t& quantity' and a
//   'std::vector<GVector<Real>>&'; the port returns { quantity, path }.
// - 'Subdivide' and 'Refine' take the midpoint by non-const reference and
//   overwrite it; the port returns { mid, changed } with a fresh midpoint
//   vector, leaving the input untouched.
// - 'std::function<void(void)> refineCallback' becomes a '() => void' field.

import { GMatrix } from './GMatrix';
import { GVector } from './GVector';
import { logAssert } from './Logger';
import { addMatrix, inverse, mulMatrix } from './Matrix';
import { Vector, add, dot, mul, sub } from './Vector';

// The result of computeGeodesic.
export interface RiemannianGeodesicResult {
    // The number of polyline vertices, pow(2,subdivisions)+1.
    quantity: number;

    // The polyline approximation to the geodesic curve.
    path: GVector[];
}

// The result of subdivide and refine: the relocated midpoint and whether the
// relocation strictly decreased the sum of the two segment lengths.
export interface RiemannianGeodesicRefineResult {
    mid: GVector;
    changed: boolean;
}

export abstract class RiemannianGeodesic {
    // Tweakable parameters.
    // 1. The integral samples are the number of samples used in the
    //    Trapezoid Rule numerical integrator.
    // 2. The search samples are the number of samples taken along a ray for
    //    the steepest descent algorithm used to refine the vertices of the
    //    polyline approximation to the geodesic curve.
    // 3. The derivative step is the value of h used for centered difference
    //    approximations df/dx = (f(x+h)-f(x-h))/(2*h) in the steepest descent
    //    algorithm.
    // 4. The number of subdivisions indicates how many times the polyline
    //    segments should be subdivided. The number of polyline vertices will
    //    be pow(2,subdivisions)+1.
    // 5. The number of refinements per subdivision. Setting this to a
    //    positive value appears necessary when the geodesic curve has a large
    //    length.
    // 6. The search radius is the distance over which the steepest descent
    //    algorithm searches for a minimum on the line whose direction is the
    //    estimated gradient. The default of 1 means the search interval is
    //    [-L,L], where L is the length of the gradient. If the search radius
    //    is r, then the interval is [-r*L,r*L].
    //
    // The derived quantities mIntegralStep, mSearchStep and mDerivativeFactor
    // are computed by the constructor from the defaults. Upstream does not
    // recompute them when the public parameters are modified, so the port
    // exposes updateDerivedParameters() for callers that change
    // integralSamples, searchSamples or derivativeStep after construction.
    integralSamples: number;  // default = 16
    searchSamples: number;    // default = 32
    derivativeStep: number;   // default = 0.0001
    subdivisions: number;     // default = 7
    refinements: number;      // default = 8
    searchRadius: number;     // default = 1.0

    // A callback that is executed during each call of refine.
    refineCallback: () => void;

    protected mDimension: number;
    protected mMetric: GMatrix;
    protected mMetricInverse: GMatrix;
    protected mChristoffel1: GMatrix[];
    protected mChristoffel2: GMatrix[];
    protected mMetricDerivative: GMatrix[];
    protected mMetricInverseExists: boolean;

    // Progress parameters that are useful to refineCallback.
    protected mSubdivide: number;
    protected mRefine: number;
    protected mCurrentQuantity: number;

    // Derived tweaking parameters.
    protected mIntegralStep: number;      // = 1/(integralSamples-1)
    protected mSearchStep: number;        // = 1/searchSamples
    protected mDerivativeFactor: number;  // = 1/(2*derivativeStep)

    // Construction. The input dimension must be two or larger.
    protected constructor(dimension: number) {
        logAssert(dimension >= 2, 'Dimension must be at least 2.');

        this.integralSamples = 16;
        this.searchSamples = 32;
        this.derivativeStep = 1e-04;
        this.subdivisions = 7;
        this.refinements = 8;
        this.searchRadius = 1;
        this.refineCallback = () => { /* no-op */ };

        this.mDimension = (dimension >= 2 ? dimension : 2);
        this.mMetric = new GMatrix(this.mDimension, this.mDimension);
        this.mMetricInverse = new GMatrix(this.mDimension, this.mDimension);
        this.mChristoffel1 = new Array<GMatrix>(this.mDimension);
        this.mChristoffel2 = new Array<GMatrix>(this.mDimension);
        this.mMetricDerivative = new Array<GMatrix>(this.mDimension);
        this.mMetricInverseExists = true;
        this.mSubdivide = 0;
        this.mRefine = 0;
        this.mCurrentQuantity = 0;
        this.mIntegralStep = 1 / (this.integralSamples - 1);
        this.mSearchStep = 1 / this.searchSamples;
        this.mDerivativeFactor = 0.5 / this.derivativeStep;

        for (let i = 0; i < this.mDimension; ++i) {
            this.mChristoffel1[i] = new GMatrix(this.mDimension, this.mDimension);
            this.mChristoffel2[i] = new GMatrix(this.mDimension, this.mDimension);
            this.mMetricDerivative[i] = new GMatrix(this.mDimension, this.mDimension);
        }
    }

    // Recompute the derived parameters from integralSamples, searchSamples
    // and derivativeStep. Upstream computes them only in the constructor.
    updateDerivedParameters(): void {
        this.mIntegralStep = 1 / (this.integralSamples - 1);
        this.mSearchStep = 1 / this.searchSamples;
        this.mDerivativeFactor = 0.5 / this.derivativeStep;
    }

    // The dimension of the manifold.
    getDimension(): number {
        return this.mDimension;
    }

    // Returns the length of the line segment connecting the points.
    computeSegmentLength(point0: GVector, point1: GVector): number {
        // The Trapezoid Rule is used for integration of the length integral.
        // The computeMetric function internally modifies mMetric, which means
        // the qForm values are actually varying even though diff does not.
        const diff = sub(point1, point0);

        // Evaluate the integrand at point0.
        this.computeMetric(point0);
        let qForm = dot(diff, mulMatrix(this.mMetric, diff));
        logAssert(qForm > 0, 'Unexpected condition.');
        let length = Math.sqrt(qForm);

        // Evaluate the integrand at point1.
        this.computeMetric(point1);
        qForm = dot(diff, mulMatrix(this.mMetric, diff));
        logAssert(qForm > 0, 'Unexpected condition.');
        length += Math.sqrt(qForm);
        length *= 0.5;

        const imax = this.integralSamples - 2;
        for (let i = 1; i <= imax; ++i) {
            // Evaluate the integrand at point0+t*(point1-point0).
            const t = this.mIntegralStep * i;
            const temp = toGVector(add(point0, mul(t, diff)));
            this.computeMetric(temp);
            qForm = dot(diff, mulMatrix(this.mMetric, diff));
            logAssert(qForm > 0, 'Unexpected condition.');
            length += Math.sqrt(qForm);
        }
        length *= this.mIntegralStep;
        return length;
    }

    // Compute the total length of the polyline. The lengths of the segments
    // are computed relative to the metric tensor.
    computeTotalLength(quantity: number, path: readonly GVector[]): number {
        logAssert(quantity >= 2, 'Path must have at least two points.');
        let length = this.computeSegmentLength(path[0], path[1]);
        for (let i = 1, ip1 = 2; ip1 < quantity; ++i, ++ip1) {
            length += this.computeSegmentLength(path[i], path[ip1]);
        }
        return length;
    }

    // Returns a polyline approximation to a geodesic curve connecting the
    // points.
    computeGeodesic(end0: GVector, end1: GVector): RiemannianGeodesicResult {
        logAssert(this.subdivisions < 32, 'Exceeds maximum iterations.');
        const quantity = (1 << this.subdivisions) + 1;
        const path = new Array<GVector>(quantity);
        for (let i = 0; i < quantity; ++i) {
            path[i] = new GVector(this.mDimension);
        }

        this.mCurrentQuantity = 2;
        path[0] = end0.clone();
        path[1] = end1.clone();

        for (this.mSubdivide = 1; this.mSubdivide <= this.subdivisions;
            ++this.mSubdivide) {
            // A subdivision essentially doubles the number of points.
            const newQuantity = 2 * this.mCurrentQuantity - 1;
            logAssert(newQuantity <= quantity, 'Unexpected condition.');

            // Copy the old points so that there are slots for the midpoints
            // during the subdivision, the slots interleaved between the old
            // points.
            for (let i = this.mCurrentQuantity - 1; i > 0; --i) {
                path[2 * i] = path[i].clone();
            }

            // Subdivide the polyline.
            for (let i = 0; i <= this.mCurrentQuantity - 2; ++i) {
                const twoI = 2 * i;
                path[twoI + 1] =
                    this.subdivide(path[twoI], path[twoI + 2]).mid;
            }

            this.mCurrentQuantity = newQuantity;

            // Refine the current polyline vertices.
            for (this.mRefine = 1; this.mRefine <= this.refinements;
                ++this.mRefine) {
                for (let i = 1, im1 = 0, ip1 = 2;
                    i <= this.mCurrentQuantity - 2; ++i, ++im1, ++ip1) {
                    path[i] = this.refine(path[im1], path[i], path[ip1]).mid;
                }
            }
        }

        logAssert(this.mCurrentQuantity === quantity, 'Unexpected condition.');
        this.mSubdivide = 0;
        this.mRefine = 0;
        this.mCurrentQuantity = 0;
        return { quantity, path };
    }

    // Start with the midpoint M of the line segment (E0,E1) and use a
    // steepest descent algorithm to move M so that
    //     Length(E0,M) + Length(M,E1) < Length(E0,E1)
    // This is essentially a relaxation scheme that inserts points into the
    // current polyline approximation to the geodesic curve.
    subdivide(end0: GVector, end1: GVector): RiemannianGeodesicRefineResult {
        const mid = toGVector(mul(0.5, add(end0, end1)));
        const save = this.refineCallback;
        this.refineCallback = () => { /* no-op */ };
        const result = this.refine(end0, mid, end1);
        this.refineCallback = save;
        return result;
    }

    // Apply the steepest descent algorithm to move the midpoint M of the line
    // segment (E0,E1) so that
    //     Length(E0,M) + Length(M,E1) < Length(E0,E1)
    // This is essentially a relaxation scheme that inserts points into the
    // current polyline approximation to the geodesic curve.
    refine(end0: GVector, mid: GVector,
        end1: GVector): RiemannianGeodesicRefineResult {
        // Estimate the gradient vector for the function
        // F(m) = Length(e0,m) + Length(m,e1).
        const temp = mid.clone();
        const gradient = new GVector(this.mDimension);
        for (let i = 0; i < this.mDimension; ++i) {
            temp.values[i] = mid.values[i] + this.derivativeStep;
            let g = this.computeSegmentLength(end0, temp);
            g += this.computeSegmentLength(temp, end1);

            temp.values[i] = mid.values[i] - this.derivativeStep;
            g -= this.computeSegmentLength(end0, temp);
            g -= this.computeSegmentLength(temp, end1);

            temp.values[i] = mid.values[i];
            gradient.values[i] = g * this.mDerivativeFactor;
        }

        // Compute the length sum for the current midpoint.
        let length0 = this.computeSegmentLength(end0, mid);
        let length1 = this.computeSegmentLength(mid, end1);
        const oldLength = length0 + length1;

        const multiplier = this.mSearchStep * this.searchRadius;
        let minLength = oldLength;
        let minPoint = mid.clone();
        for (let i = -this.searchSamples; i <= this.searchSamples; ++i) {
            const tRay = multiplier * i;
            const pRay = toGVector(sub(mid, mul(tRay, gradient)));
            length0 = this.computeSegmentLength(end0, pRay);
            length1 = this.computeSegmentLength(end1, pRay);
            const newLength = length0 + length1;
            if (newLength < minLength) {
                minLength = newLength;
                minPoint = pRay;
            }
        }

        this.refineCallback();
        return { mid: minPoint, changed: minLength < oldLength };
    }

    // Information to be used during the callback.
    getSubdivisionStep(): number {
        return this.mSubdivide;
    }

    getRefinementStep(): number {
        return this.mRefine;
    }

    getCurrentQuantity(): number {
        return this.mCurrentQuantity;
    }

    // Curvature computations to measure how close the approximating polyline
    // is to a geodesic.

    // Returns the total curvature of the line segment connecting the points.
    computeSegmentCurvature(point0: GVector, point1: GVector): number {
        // The Trapezoid Rule is used for integration of the curvature
        // integral. The computeIntegrand function internally modifies
        // mMetric, which means the curvature values are actually varying even
        // though diff does not.
        const diff = sub(point1, point0);

        // Evaluate the integrand at point0.
        let curvature = this.computeIntegrand(point0, diff);

        // Evaluate the integrand at point1.
        curvature += this.computeIntegrand(point1, diff);
        curvature *= 0.5;

        const imax = this.integralSamples - 2;
        for (let i = 1; i <= imax; ++i) {
            // Evaluate the integrand at point0+t*(point1-point0).
            const t = this.mIntegralStep * i;
            const temp = toGVector(add(point0, mul(t, diff)));
            curvature += this.computeIntegrand(temp, diff);
        }
        curvature *= this.mIntegralStep;
        return curvature;
    }

    // Compute the total curvature of the polyline. The curvatures of the
    // segments are computed relative to the metric tensor.
    computeTotalCurvature(quantity: number,
        path: readonly GVector[]): number {
        logAssert(quantity >= 2, 'Path must have at least two points.');
        let curvature = this.computeSegmentCurvature(path[0], path[1]);
        for (let i = 1, ip1 = 2; ip1 < quantity; ++i, ++ip1) {
            curvature += this.computeSegmentCurvature(path[i], path[ip1]);
        }
        return curvature;
    }

    // Support for computeSegmentCurvature.
    protected computeIntegrand(pos: GVector, der: Vector): number {
        this.computeMetric(pos);
        this.computeChristoffel1(pos);
        this.computeMetricInverse();
        this.computeChristoffel2();

        // g_{ij}*der_{i}*der_{j}
        const qForm0 = dot(der, mulMatrix(this.mMetric, der));
        logAssert(qForm0 > 0, 'Unexpected condition.');

        // gamma_{kij}*der_{k}*der_{i}*der_{j}
        let mat = new GMatrix(this.mDimension, this.mDimension);
        for (let k = 0; k < this.mDimension; ++k) {
            mat = toGMatrix(addMatrix(mat,
                mulMatrix(this.mChristoffel1[k], der.values[k])));
        }
        // This product can be negative because mat is not guaranteed to be
        // positive semidefinite. No assertion is added.
        const qForm1 = dot(der, mulMatrix(mat, der));

        const ratio = -qForm1 / qForm0;

        // Compute the acceleration.
        const acc = mul(ratio, der);
        for (let k = 0; k < this.mDimension; ++k) {
            acc.values[k] += dot(der, mulMatrix(this.mChristoffel2[k], der));
        }

        // Compute the curvature.
        return Math.sqrt(dot(acc, mulMatrix(this.mMetric, acc)));
    }

    // Compute the metric tensor for the specified point. Derived classes are
    // responsible for implementing this function.
    protected abstract computeMetric(point: GVector): void;

    // Compute the Christoffel symbols of the first kind for the current
    // point. Derived classes are responsible for implementing this function.
    protected abstract computeChristoffel1(point: GVector): void;

    // Compute the inverse of the current metric tensor. The function returns
    // 'true' iff the inverse exists.
    protected computeMetricInverse(): boolean {
        const result = inverse(this.mMetric);
        this.mMetricInverse = toGMatrix(result.inverse);
        this.mMetricInverseExists = result.invertible;
        return this.mMetricInverseExists;
    }

    // Compute the derivative of the metric tensor for the current state. This
    // is a triply indexed quantity, the values computed using the Christoffel
    // symbols of the first kind.
    protected computeMetricDerivative(): void {
        for (let derivative = 0; derivative < this.mDimension; ++derivative) {
            for (let i0 = 0; i0 < this.mDimension; ++i0) {
                for (let i1 = 0; i1 < this.mDimension; ++i1) {
                    this.mMetricDerivative[derivative].set(i0, i1,
                        this.mChristoffel1[derivative].get(i0, i1) +
                        this.mChristoffel1[derivative].get(i1, i0));
                }
            }
        }
    }

    // Compute the Christoffel symbols of the second kind for the current
    // state. The values depend on the inverse of the metric tensor, so they
    // may be computed only when the inverse exists. The function returns
    // 'true' whenever the inverse metric tensor exists.
    protected computeChristoffel2(): boolean {
        for (let i2 = 0; i2 < this.mDimension; ++i2) {
            for (let i0 = 0; i0 < this.mDimension; ++i0) {
                for (let i1 = 0; i1 < this.mDimension; ++i1) {
                    let value = 0;
                    for (let j = 0; j < this.mDimension; ++j) {
                        value += this.mMetricInverse.get(i2, j) *
                            this.mChristoffel1[j].get(i0, i1);
                    }
                    this.mChristoffel2[i2].set(i0, i1, value);
                }
            }
        }
        return this.mMetricInverseExists;
    }
}

// The Vector free functions return the base Vector type; the geodesic
// algorithm passes the results back into GVector-typed slots, so convert.
function toGVector(v: Vector): GVector {
    return GVector.fromArray(v.values);
}

function toGMatrix(m: { numRows: number, numCols: number, values: number[] } |
    { numRows: number, numCols: number, values: readonly number[] }): GMatrix {
    const result = new GMatrix(m.numRows, m.numCols);
    for (let i = 0; i < result.numElements; ++i) {
        result.values[i] = m.values[i];
    }
    return result;
}
