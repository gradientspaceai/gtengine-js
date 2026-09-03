// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprCone3EllipseAndPoints.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An infinite single-sided cone is fit to a 3D ellipse that is known to be
// the intersection of a plane with the cone. The ellipse itself is not
// enough information, producing the cone vertex and cone direction as a
// function of the cone angle. Additional points on the cone are required
// to determine the cone angle. The algorithm description is
// https://www.geometrictools.com/Documentation/FitConeToEllipseAndPoints.pdf
//
// Port notes:
// * The upstream header contains two classes; both are ported here.
// * The nested 'Control' struct becomes the exported class
//   'ApprCone3EllipseAndPointsControl' (global export uniqueness), whose
//   constructor supplies the upstream default values.
// * The private static 'ComputeCone' is a module-private function.
// * 'ApprCone3ExtractEllipses::Extract' writes the ellipses through an output
//   reference; the port returns them. The upstream 'OBBTree'/'OBBNode' type
//   aliases are dropped in favor of the concrete types.

import { ApprEllipse2 } from './ApprEllipse2.js';
import { ApprGaussian3 } from './ApprGaussian3.js';
import { BVTreeNode } from './BVTree.js';
import { Cone } from './Cone.js';
import type { Cone3 } from './Cone.js';
import { GTE_C_HALF_PI } from './Constants.js';
import { Ellipse3 } from './Ellipse3.js';
import { Hyperellipsoid } from './Hyperellipsoid.js';
import type { Ellipse2 } from './Hyperellipsoid.js';
import { Hyperplane } from './Hyperplane.js';
import type { Plane3 } from './Hyperplane.js';
import { logAssert } from './Logger.js';
import { Minimize1 } from './Minimize1.js';
import type { OrientedBox3 } from './OrientedBox.js';
import { OrientedBoxBV } from './OrientedBoxBV.js';
import { OrientedBoxTreeOfPoints } from './OrientedBoxTreeOfPoints.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

// The default control parameters appear to be reasonable for applications,
// but they are exposed to the caller for tuning.
export class ApprCone3EllipseAndPointsControl {
    // The least-squares error function is updated with the penalty value for
    // a points[i] that is below the plane supporting the cone; that is, when
    // the dot product Dot(coneDirection, points[i] - coneVertex) < 0.
    penalty: number;

    // Parameters for Minimize1.
    maxSubdivisions: number;
    maxBisections: number;
    epsilon: number;
    tolerance: number;

    // Search for the minimum on [0 + padding, pi/2 - padding] to avoid
    // divisions by zero of the least-squares error function at the endpoints
    // of [0,pi/2].
    padding: number;

    constructor() {
        this.penalty = 1;
        this.maxSubdivisions = 8;
        this.maxBisections = 64;
        this.epsilon = 1e-08;
        this.tolerance = 1e-04;
        this.padding = 1e-03;
    }

    validParameters(): boolean {
        return this.penalty > 0
            && this.maxSubdivisions > 0
            && this.maxBisections > 0
            && this.epsilon > 0
            && this.tolerance > 0
            && this.padding > 0;
    }
}

export class ApprCone3EllipseAndPoints {
    // The ellipse must be the intersection of a plane with the cone. In an
    // application, typically the ellipse is estimated from point samples of
    // the intersection which are then fitted with the ellipse.
    static fit(ellipse: Ellipse3, points: readonly Vector[],
        control: ApprCone3EllipseAndPointsControl =
            new ApprCone3EllipseAndPointsControl()): Cone3 {
        logAssert(control.validParameters(), 'Invalid control parameter.');

        // Upstream divides the accumulated error by points.size() and the
        // ellipse extents appear in denominators of ComputeCone, neither of
        // which is guarded upstream. The port asserts, following the
        // established Appr*/Cont* precedent of an explicit empty-input guard.
        logAssert(points.length > 0, 'ApprCone3EllipseAndPoints: no points.');
        logAssert(ellipse.extent.get(0) > 0 && ellipse.extent.get(1) > 0,
            'ApprCone3EllipseAndPoints: ellipse extents must be positive.');

        // The sign pair (sigma0,sigma1) is captured by the error function and
        // is varied between the four minimizer runs, as upstream does with
        // its lambda captures by reference.
        let sigma0 = 0;
        let sigma1 = 0;

        const F = (theta: number): number => {
            const cone = computeCone(theta, sigma0, sigma1, ellipse);

            let error = 0;
            for (const point of points) {
                const diff = sub(point, cone.ray.origin);
                const d = dot(cone.ray.direction, diff);
                if (d >= 0) {
                    const sqrLen = dot(diff, diff);
                    const quad = d * d - cone.cosAngleSqr * sqrLen;
                    error += quad * quad;
                }
                else {
                    error += control.penalty;
                }
            }

            return Math.sqrt(error) / points.length;
        };

        const minimizer = new Minimize1(F, control.maxSubdivisions,
            control.maxBisections, control.epsilon, control.tolerance);
        const t0 = control.padding;
        const t1 = GTE_C_HALF_PI - control.padding;
        let minError = -1;
        let minCone = new Cone(3);

        const signs: ReadonlyArray<readonly [number, number]> =
            [[1, 1], [1, -1], [-1, 1], [-1, -1]];
        for (const [s0, s1] of signs) {
            sigma0 = s0;
            sigma1 = s1;
            const { tMin, fMin } = minimizer.getMinimum(t0, t1);
            if (t0 < tMin && tMin < t1) {
                if (minError === -1 || fMin < minError) {
                    minError = fMin;
                    minCone = computeCone(tMin, sigma0, sigma1, ellipse);
                }
            }
        }

        logAssert(minError !== -1, 'Failed to find fitted cone.');

        return minCone;
    }
}

// For a cone angle theta, the ellipse determines the cone vertex and cone
// axis direction up to the sign choices (sigma0,sigma1) for sin(phi) and
// cos(phi), where phi is the angle between the cone axis and the ellipse
// plane normal. The relationship sin(phi) = e * cos(theta) holds, where e is
// the eccentricity of the ellipse.
function computeCone(theta: number, sigma0: number, sigma1: number,
    ellipse: Ellipse3): Cone3 {
    const C = ellipse.center;
    const N = ellipse.normal;
    const U = ellipse.axis[0];
    const a = ellipse.extent.get(0);
    const b = ellipse.extent.get(1);
    const bDivA = b / a;
    const eSqr = Math.max(0, 1 - bDivA * bDivA);
    const omesqr = 1 - eSqr;
    const e = Math.sqrt(eSqr);

    const snTheta = Math.sin(theta);
    const csTheta = Math.cos(theta);
    const snPhi = sigma0 * e * csTheta;
    const snPhiSqr = snPhi * snPhi;
    const csPhi = sigma1 * Math.sqrt(Math.max(0, 1 - snPhiSqr));
    const h = a * omesqr * csTheta / (snTheta * Math.abs(csPhi));
    const D = add(mul(csPhi, N), mul(snPhi, U));
    const snThetaSqr = snTheta * snTheta;
    const csThetaSqr = csTheta * csTheta;
    const Q = sub(C, mul((h * snPhi * snThetaSqr) / (csThetaSqr - snPhiSqr), U));
    const K = sub(Q, mul(h, D));

    const cone = new Cone(3);
    cone.makeInfiniteCone();
    cone.setAngle(theta);
    cone.ray.origin = K;
    cone.ray.direction = D;
    return cone;
}


// If the points contain only elliptical cross sections of intersection of
// planes with the cone, extract the ellipses so that one of them can be used
// as input to ApprCone3EllipseAndPoints.
export class ApprCone3ExtractEllipses {
    private mBoxExtentEpsilon: number;
    private mCosAngleEpsilon: number;
    private mOBBTree: readonly BVTreeNode<OrientedBoxBV>[];
    private mPlanes: Plane3[];
    private mIndices: number[][];
    private mBoxes: OrientedBox3[];
    private mEllipses: Ellipse3[];

    constructor() {
        this.mBoxExtentEpsilon = 0;
        this.mCosAngleEpsilon = 0;
        this.mOBBTree = [];
        this.mPlanes = [];
        this.mIndices = [];
        this.mBoxes = [];
        this.mEllipses = [];
    }

    // The boxExtentEpsilon determines when a box is deemed "flat." The
    // cosAngleEpsilon is used to decide when two flat boxes are in the same
    // plane. Upstream writes the ellipses through an output reference; the
    // port returns them (they remain available from getEllipses()).
    extract(points: readonly Vector[], boxExtentEpsilon: number,
        cosAngleEpsilon: number): Ellipse3[] {
        this.mBoxExtentEpsilon = Math.max(boxExtentEpsilon, 0);
        this.mCosAngleEpsilon = Math.max(cosAngleEpsilon, 0);
        this.mOBBTree = [];
        this.mPlanes = [];
        this.mIndices = [];
        this.mBoxes = [];
        this.mEllipses = [];

        this.createOBBTree(points);
        this.locatePlanes(0);
        this.associatePointsWithPlanes(points);

        this.mEllipses = new Array<Ellipse3>(this.mIndices.length);
        for (let i = 0; i < this.mIndices.length; ++i) {
            this.mEllipses[i] = computeEllipse(points, this.mIndices[i]);
        }

        return this.mEllipses.slice();
    }

    // Access the ellipses extracted from the input points.
    getEllipses(): readonly Ellipse3[] {
        return this.mEllipses;
    }

    // Accessors for informational purposes.
    getOBBTree(): readonly BVTreeNode<OrientedBoxBV>[] {
        return this.mOBBTree;
    }

    getPlanes(): readonly Plane3[] {
        return this.mPlanes;
    }

    getIndices(): readonly number[][] {
        return this.mIndices;
    }

    getBoxes(): readonly OrientedBox3[] {
        return this.mBoxes;
    }

    private createOBBTree(points: readonly Vector[]): void {
        const creator = new OrientedBoxTreeOfPoints();
        creator.create(points);
        this.mOBBTree = creator.getNodes();
    }

    // A node whose bounding box is flat (some extent is at most
    // mBoxExtentEpsilon) supports a plane of the point set. Descend only into
    // the nodes that are not yet flat.
    private locatePlanes(nodeIndex: number): void {
        const node = this.mOBBTree[nodeIndex];
        if (node.maxIndex >= node.minIndex + 2) {
            const box = node.boundingVolume.box;
            for (let j = 0; j < 3; ++j) {
                if (box.extent.get(j) <= this.mBoxExtentEpsilon) {
                    this.mBoxes.push(box.clone());
                    const plane = Hyperplane.fromNormalOrigin(box.axis[j],
                        box.center);
                    this.processPlane(plane);
                    return;
                }
            }
        }

        if (node.leftChild !== BVTreeNode.invalid) {
            this.locatePlanes(node.leftChild);
        }

        if (node.rightChild !== BVTreeNode.invalid) {
            this.locatePlanes(node.rightChild);
        }
    }

    // Store the plane unless it is effectively the same as one already
    // stored. Planes with nearly opposite normals are compared using the
    // negated constant.
    private processPlane(plane: Plane3): void {
        const epsilon = this.mCosAngleEpsilon;
        const oneMinusEpsilon = 1 - epsilon;

        for (let i = 0; i < this.mPlanes.length; ++i) {
            let cosAngle = dot(plane.normal, this.mPlanes[i].normal);
            let absDiff: number;
            if (cosAngle > 0) {
                absDiff = Math.abs(plane.constant - this.mPlanes[i].constant);
                if (cosAngle >= oneMinusEpsilon && absDiff <= epsilon) {
                    // The planes are effectively the same.
                    return;
                }
            }
            else {
                cosAngle = -cosAngle;
                absDiff = Math.abs(plane.constant + this.mPlanes[i].constant);
                if (cosAngle >= oneMinusEpsilon && absDiff <= epsilon) {
                    // The planes are effectively the same.
                    return;
                }
            }
        }

        this.mPlanes.push(plane);
    }

    // Each point is assigned to the plane nearest to it.
    private associatePointsWithPlanes(points: readonly Vector[]): void {
        // Upstream initializes minJ to size_t(-1) and assigns it only inside
        // the loop over mPlanes. When no flat box was found, mPlanes is empty,
        // the inner loop never executes and 'mIndices[minJ].push_back(i)'
        // indexes out of bounds on the very first point (undefined behavior).
        // The port returns early instead, leaving mIndices empty so that no
        // ellipses are extracted.
        this.mIndices = [];
        if (this.mPlanes.length === 0) {
            return;
        }

        for (let j = 0; j < this.mPlanes.length; ++j) {
            this.mIndices.push([]);
        }

        for (let i = 0; i < points.length; ++i) {
            const point = points[i];
            let minDistance = Number.MAX_VALUE;
            let minJ = 0;
            for (let j = 0; j < this.mPlanes.length; ++j) {
                const plane = this.mPlanes[j];
                const diff = sub(point, plane.origin);
                const distance = Math.abs(dot(plane.normal, diff));
                if (distance < minDistance) {
                    minDistance = distance;
                    minJ = j;
                }
            }

            this.mIndices[minJ].push(i);
        }

        // A node holding exactly three points always has a flat bounding box
        // (three points are coplanar), so LocatePlanes can emit spurious
        // planes that end up with no points of their own. Upstream then fits
        // a Gaussian and an ellipse to an empty index set. The port discards
        // such planes, which keeps mPlanes, mIndices and mEllipses parallel
        // and well defined.
        const keep: number[] = [];
        for (let j = 0; j < this.mIndices.length; ++j) {
            if (this.mIndices[j].length > 0) {
                keep.push(j);
            }
        }

        if (keep.length !== this.mPlanes.length) {
            this.mPlanes = keep.map((j) => this.mPlanes[j]);
            this.mIndices = keep.map((j) => this.mIndices[j]);
        }
    }
}

// Fit the indexed subset of the points with an ellipse in the plane of those
// points.
function computeEllipse(points: readonly Vector[],
    indices: readonly number[]): Ellipse3 {
    // A plane's index set can be empty when every point is closer to some
    // other plane. Upstream then fits a Gaussian and an ellipse to zero
    // points, which produces garbage (and an oriented-box fit of an empty
    // point set inside ApprEllipse2). The port reports the condition,
    // following the established empty-input guard precedent.
    logAssert(indices.length > 0,
        'ApprCone3ExtractEllipses: no points associated with a plane.');

    // Fit the points with a 3D Gaussian distribution. The eigenvalues are
    // computed in nondecreasing order, which means the smallest eigenvalue
    // corresponds to the normal vector gbox.axis[0] of the plane of the
    // points. Use gbox.axis[1] and gbox.axis[2] as the spanners of the plane
    // of the points.
    const gfitter = new ApprGaussian3();
    gfitter.fitIndexed(points, indices);
    const gbox = gfitter.getParameters();

    // Project the points onto the plane as 2-tuples.
    const projections = new Array<Vector>(indices.length);
    for (let i = 0; i < indices.length; ++i) {
        const diff = sub(points[indices[i]], gbox.center);
        projections[i] = Vector.fromArray([
            dot(gbox.axis[1], diff), dot(gbox.axis[2], diff)]);
    }

    // Fit the projected points with a 2D ellipse.
    const efitter = new ApprEllipse2();
    const numIterations = 1024;
    const useEllipseForInitialGuess = false;
    const ellipse2: Ellipse2 = new Hyperellipsoid(2);
    efitter.compute(projections, numIterations, useEllipseForInitialGuess,
        ellipse2);

    // Lift the 2D ellipse to a 3D ellipse.
    const ellipse3 = new Ellipse3();
    ellipse3.center = add(gbox.center,
        add(mul(ellipse2.center.get(0), gbox.axis[1]),
            mul(ellipse2.center.get(1), gbox.axis[2])));
    ellipse3.normal = gbox.axis[0].clone();
    ellipse3.axis[0] = add(mul(ellipse2.axis[0].get(0), gbox.axis[1]),
        mul(ellipse2.axis[0].get(1), gbox.axis[2]));
    ellipse3.axis[1] = add(mul(ellipse2.axis[1].get(0), gbox.axis[1]),
        mul(ellipse2.axis[1].get(1), gbox.axis[2]));
    ellipse3.extent = ellipse2.extent.clone();
    return ellipse3;
}
