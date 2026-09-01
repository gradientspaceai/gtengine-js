// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Hyperplane.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The hyperplane is represented as Dot(U, X - P) = 0 where U is a unit-length
// normal vector, P is the hyperplane origin, and X is any point on the
// hyperplane. The user must ensure that the normal vector is unit length. The
// hyperplane constant is c = Dot(U, P) so that Dot(U, X) = c. If P is not
// specified when constructing a hyperplane, it is chosen to be the point on
// the plane closest to the origin, P = c * U.
//
// NOTE: You cannot set 'origin' and 'constant' independently. Use the static
// factories instead.
//
//   // Construct from normal N and constant c.
//   const plane = Hyperplane.fromNormalConstant(N, c);  // origin = c * N
//
//   // Construct from normal N and origin P.
//   const plane = Hyperplane.fromNormalOrigin(N, P);  // constant = Dot(N, P)
//
//   const plane = new Hyperplane(3);  // N = (0,0,1), P = (0,0,0), c = 0
//   plane.normal = ...;
//   plane.constant = 3;
//   // If you consume plane now, the origin and constant are inconsistent.
//   // Instead use
//   plane = Hyperplane.fromNormalConstant(Vector.fromArray([0, 0, 1]), 3);
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Hyperplane(n)' for the default
// constructor, named static factories that copy their Vector arguments,
// comparison methods). The upstream 'Hyperplane(std::array<Vector<N,T>,N>)'
// constructor becomes the static factory fromPoints().

import { logAssert } from './Logger';
import { Matrix } from './Matrix';
import { SingularValueDecomposition } from './SingularValueDecomposition';
import { Vector, dot, mul, sub } from './Vector';
import { computeOrthogonalComplement2 } from './Vector2';
import { unitCross } from './Vector3';

export class Hyperplane {
    // Public member access.
    normal: Vector;
    origin: Vector;
    constant: number;

    // The port of the default constructor, which sets the normal to
    // (0,...,0,1), the origin to (0,...,0) and the constant to zero. The
    // dimension N of the C++ template is a constructor argument here.
    constructor(n: number) {
        logAssert(n >= 2, 'Invalid dimension.');
        this.normal = new Vector(n);
        this.normal.makeUnit(n - 1);
        this.origin = new Vector(n);
        this.constant = 0;
    }

    // The port of 'Hyperplane(inNormal, inConstant)'.
    static fromNormalConstant(inNormal: Vector,
        inConstant: number): Hyperplane {
        const hyperplane = new Hyperplane(inNormal.size);
        hyperplane.normal = inNormal.clone();
        hyperplane.origin = mul(inNormal, inConstant);
        hyperplane.constant = inConstant;
        return hyperplane;
    }

    // The port of 'Hyperplane(inNormal, inOrigin)'.
    static fromNormalOrigin(inNormal: Vector, inOrigin: Vector): Hyperplane {
        logAssert(inNormal.size === inOrigin.size,
            'Hyperplane: mismatched sizes.');
        const hyperplane = new Hyperplane(inNormal.size);
        hyperplane.normal = inNormal.clone();
        hyperplane.origin = inOrigin.clone();
        hyperplane.constant = dot(inNormal, inOrigin);
        return hyperplane;
    }

    // U is a unit-length vector in the orthogonal complement of the set
    // {p[1]-p[0],...,p[n-1]-p[0]} and c = Dot(U,p[0]), where the p[i] are
    // points on the hyperplane. The array must contain exactly N points of
    // dimension N.
    static fromPoints(p: readonly Vector[]): Hyperplane {
        const n = p.length;
        logAssert(n >= 2, 'Invalid dimension.');
        for (let i = 0; i < n; ++i) {
            logAssert(p[i].size === n, 'Hyperplane: mismatched sizes.');
        }

        const hyperplane = new Hyperplane(n);
        if (n === 3) {
            // The port of the Dimension == 3 specialization of
            // ComputeFromPoints.
            const edge0 = sub(p[1], p[0]);
            const edge1 = sub(p[2], p[0]);
            hyperplane.normal = unitCross(edge0, edge1);
        }
        else if (n === 2) {
            // Port fix for an upstream bug. Upstream routes N = 2 through the
            // generic ComputeFromPoints, which constructs a
            // SingularValueDecomposition(2, 1, 32); that constructor has
            // LogAssert(mNumCols >= 2 && ...), so Hyperplane<2,T> cannot be
            // built from points at all upstream. The 1-dimensional orthogonal
            // complement of a single edge in 2D is computed directly with the
            // GTE helper ComputeOrthogonalComplement for Vector2, which sets
            // v[1] = -Perp(v[0]) and orthonormalizes. Degenerate input (the
            // two points coincide) yields the zero normal, matching the
            // UnitCross behavior of the N = 3 path.
            const v: Vector[] = [sub(p[1], p[0]), new Vector(2)];
            computeOrthogonalComplement2(1, v);
            hyperplane.normal = v[1];
        }
        else {
            // The port of the generic ComputeFromPoints (Dimension != 3).
            const edge = Matrix.zero(n, n - 1);
            for (let i0 = 0, i1 = 1; i1 < n; i0 = i1++) {
                edge.setCol(i0, sub(p[i1], p[0]));
            }

            // Compute the 1-dimensional orthogonal complement of the edges of
            // the simplex formed by the points p[]. Port fix for an upstream
            // bug: upstream calls svd.Solve(&edge[0], -1), but the second
            // parameter of Solve is the positive 'multiplier' used to form
            // the convergence threshold (it was a sort type in an older API),
            // and Solve has LogAssert(multiplier > 0). The port uses the
            // documented default multiplier.
            const svd = new SingularValueDecomposition(n, n - 1, 32);
            svd.solve(edge.values);
            hyperplane.normal = Vector.fromArray(svd.getUColumn(n - 1));
        }

        hyperplane.constant = dot(hyperplane.normal, p[0]);
        hyperplane.origin = mul(hyperplane.normal, hyperplane.constant);
        return hyperplane;
    }

    // The dimension N of the hyperplane.
    get dimension(): number {
        return this.normal.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Hyperplane {
        const hyperplane = new Hyperplane(this.dimension);
        hyperplane.normal = this.normal.clone();
        hyperplane.origin = this.origin.clone();
        hyperplane.constant = this.constant;
        return hyperplane;
    }

    // Comparisons to support sorted containers.
    equals(hyperplane: Hyperplane): boolean {
        return this.normal.equals(hyperplane.normal)
            && this.origin.equals(hyperplane.origin)
            && this.constant === hyperplane.constant;
    }

    notEquals(hyperplane: Hyperplane): boolean {
        return !this.equals(hyperplane);
    }

    lessThan(hyperplane: Hyperplane): boolean {
        if (this.normal.lessThan(hyperplane.normal)) {
            return true;
        }

        if (this.normal.greaterThan(hyperplane.normal)) {
            return false;
        }

        if (this.origin.lessThan(hyperplane.origin)) {
            return true;
        }

        if (this.origin.greaterThan(hyperplane.origin)) {
            return false;
        }

        return this.constant < hyperplane.constant;
    }

    lessThanOrEqual(hyperplane: Hyperplane): boolean {
        return !hyperplane.lessThan(this);
    }

    greaterThan(hyperplane: Hyperplane): boolean {
        return hyperplane.lessThan(this);
    }

    greaterThanOrEqual(hyperplane: Hyperplane): boolean {
        return !this.lessThan(hyperplane);
    }
}

// Alias for convenience (the port of the upstream template alias).
export type Plane3 = Hyperplane;
