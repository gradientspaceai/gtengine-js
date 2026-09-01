// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint3ConvexPolyhedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a convex polyhedron in 3D. The
// algorithm is based on using an LCP solver for the convex quadratic
// programming problem. For details, see
// https://www.geometrictools.com/Documentation/ConvexQuadraticProgramming.pdf
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector3<T>, ConvexPolyhedron3<T>>'
// becomes the class DistPoint3ConvexPolyhedron3 with the result type
// DistPoint3ConvexPolyhedron3Result. The 'std::unique_ptr<LCPSolver<T>>'
// member becomes an 'LCPSolver | null'; the per-query solver of the
// numTriangles == 0 construction is created and dropped in compute(), as
// upstream does. Upstream's 'bool LCPSolver::Solve(q, M, w, z)' returns
// '{ success, result, w, z }' in the port.

import { ConvexPolyhedron3 } from './ConvexPolyhedron3';
import type { DCPQuery } from './DCPQuery';
import { LCPSolver } from './LCPSolver';
import { Vector, dot, hlift, hproject, sub } from './Vector';

export interface DistPoint3ConvexPolyhedron3Result {
    queryIsSuccessful: boolean;

    // These members are valid only when queryIsSuccessful is true; otherwise,
    // they are all set to zero.
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest polyhedron
    // point.
    closest: [Vector, Vector];

    // The number of iterations used by LCPSolver regardless of whether the
    // query is successful.
    numLCPIterations: number;
}

export class DistPoint3ConvexPolyhedron3
    implements DCPQuery<Vector, ConvexPolyhedron3,
    DistPoint3ConvexPolyhedron3Result> {
    private mMaxLCPIterations: number;
    private mLCP: LCPSolver | null;

    // Construction. If you have no knowledge of the number of faces for the
    // convex polyhedra you plan on applying the query to, pass 'numTriangles'
    // of zero. This is a request to compute() to create the LCP solver for
    // each query, and this requires memory allocation and deallocation per
    // query. If you plan on applying the query multiple times to a single
    // polyhedron, even if the vertices of the polyhedron are modified for
    // each query, then pass 'numTriangles' to be the number of triangle faces
    // for that polyhedron. This lets compute() know to create the LCP solver
    // once at construction time, thus avoiding the memory management costs
    // during the query.
    constructor(numTriangles: number = 0) {
        if (numTriangles > 0) {
            const n = numTriangles + 3;
            this.mLCP = new LCPSolver(n);
            this.mMaxLCPIterations = this.mLCP.getMaxIterations();
        }
        else {
            this.mLCP = null;
            this.mMaxLCPIterations = 0;
        }
    }

    // Default maximum iterations is 144 (n = 12, maxIterations = n*n). If the
    // solver fails to converge, try increasing the maximum number of
    // iterations.
    setMaxLCPIterations(maxLCPIterations: number): void {
        this.mMaxLCPIterations = maxLCPIterations;
        if (this.mLCP !== null) {
            this.mLCP.setMaxIterations(this.mMaxLCPIterations);
        }
    }

    compute(point: Vector, polyhedron: ConvexPolyhedron3):
        DistPoint3ConvexPolyhedron3Result {
        const result: DistPoint3ConvexPolyhedron3Result = {
            queryIsSuccessful: false,
            distance: 0,
            sqrDistance: 0,
            closest: [new Vector(3), new Vector(3)],
            numLCPIterations: 0
        };

        const numTriangles = polyhedron.planes.length;
        if (numTriangles === 0) {
            // The polyhedron planes and aligned box need to be created.
            return result;
        }

        const n = numTriangles + 3;

        // Translate the point and convex polyhedron so that the polyhedron is
        // in the first octant. The translation is not explicit; rather, the q
        // and M for the LCP are initialized using the translation
        // information.
        const boxMin = polyhedron.alignedBox.min;
        const hmin = hlift(boxMin, 1);

        const q = new Array<number>(n).fill(0);
        for (let r = 0; r < 3; ++r) {
            q[r] = boxMin.values[r] - point.values[r];
        }
        for (let r = 3, t = 0; r < n; ++r, ++t) {
            q[r] = -dot(polyhedron.planes[t], hmin);
        }

        const M = new Array<number>(n * n).fill(0);
        M[0] = 1;
        M[1] = 0;
        M[2] = 0;
        M[n] = 0;
        M[n + 1] = 1;
        M[n + 2] = 0;
        M[2 * n] = 0;
        M[2 * n + 1] = 0;
        M[2 * n + 2] = 1;
        for (let t = 0, c = 3; t < numTriangles; ++t, ++c) {
            const normal = hproject(polyhedron.planes[t]);
            for (let r = 0; r < 3; ++r) {
                M[c + n * r] = normal.values[r];
                M[r + n * c] = -normal.values[r];
            }
        }
        for (let r = 3; r < n; ++r) {
            for (let c = 3; c < n; ++c) {
                M[c + n * r] = 0;
            }
        }

        const needsLCP = (this.mLCP === null);
        if (needsLCP) {
            this.mLCP = new LCPSolver(n);
            if (this.mMaxLCPIterations > 0) {
                this.mLCP.setMaxIterations(this.mMaxLCPIterations);
            }
        }
        const lcp = this.mLCP as LCPSolver;

        const output = lcp.solve(q, M);
        if (output.success) {
            result.queryIsSuccessful = true;
            result.closest[0] = point.clone();
            for (let i = 0; i < 3; ++i) {
                result.closest[1].values[i] = output.z[i] + boxMin.values[i];
            }

            const diff = sub(result.closest[1], result.closest[0]);
            result.sqrDistance = dot(diff, diff);
            result.distance = Math.sqrt(result.sqrDistance);
        }
        else {
            // If you reach this case, the maximum number of iterations was
            // not specified to be large enough or there is a problem due to
            // floating-point rounding errors.
            result.queryIsSuccessful = false;
        }

        result.numLCPIterations = lcp.getNumIterations();
        if (needsLCP) {
            this.mLCP = null;
        }
        return result;
    }
}
