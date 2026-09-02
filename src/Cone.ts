// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Cone.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An infinite cone is defined by a vertex V, a unit-length direction D and an
// angle A with 0 < A < pi/2. A point X is on the cone when
//   Dot(D, X - V) = |X - V| * cos(A)
// A solid cone includes points on the cone and in the region that contains
// the cone ray V + h * D for h >= 0. It is defined by
//   Dot(D, X - V) >= |X - V| * cos(A)
// The height of any point Y in space relative to the cone is defined by
// h = Dot(D, Y - V), which is the signed length of the projection of X - V
// onto the cone axis. Observe that we have restricted the cone definition to
// an acute angle A, so |X - V| * cos(A) >= 0; therefore, points on or inside
// the cone have nonnegative heights: Dot(D, X - V) >= 0. The infinite solid
// cone is the "positive cone," which means that the non-vertex points inside
// the cone have positive heights. Although rare in computer graphics, one
// might also want to consider the "negative cone," which is defined by
//   -Dot(D, X - V) <= -|X - V| * cos(A)
// The non-vertex points inside this cone have negative heights.
//
// For many of the geometric queries involving cones, we can avoid the square
// root computation implied by |X - V|. The positive cone is defined by
//   Dot(D, X - V)^2 >= |X - V|^2 * cos(A)^2,
// which is a quadratic inequality, but the squaring of the terms leads to an
// inequality that includes points X in the negative cone. When using the
// quadratic inequality for the positive cone, we need to include also the
// constraint Dot(D, X - V) >= 0.
//
// There are four different types of cones. They all involve V, D and A. The
// differences are based on restrictions to the heights of the cone points.
// The height range is defined to be the interval of possible heights, say,
// [hmin,hmax] with 0 <= hmin < hmax <= +infinity.
//     1. infinite cone: hmin = 0, hmax = +infinity
//     2. infinite truncated cone: hmin > 0, hmax = +infinity
//     3. finite cone: hmin >= 0, hmax < +infinity
//     4. frustum of a cone: hmin > 0, hmax < +infinity
// The infinite truncated cone is truncated for h-minimum; the radius of the
// disk at h-minimum is rmin = hmin * tan(A). The finite cone is truncated for
// h-maximum; the radius of the disk at h-maximum is rmax = hmax * tan(A).
// The frustum of a cone is truncated both for h-minimum and h-maximum.
//
// A technical problem when creating a data structure to represent a cone is
// deciding how to represent +infinity in the height range. Upstream cannot
// use a floating-point infinity because the arbitrary precision types
// BSNumber and BSRational have no representation for it, so +infinity is
// represented by setting the maxHeight member to -1. The member functions
// isFinite() and isInfinite() compare maxHeight to -1 and report the correct
// state. This port keeps that convention.
//
// The consequence is that comparisons between heights require extra logic.
// For example, the point-in-cone test using the quadratic inequality is
//   const delta = sub(point, cone.ray.origin);
//   const h = dot(cone.ray.direction, delta);
//   const pointInCone = cone.heightInRange(h) &&
//       h * h >= dot(delta, delta) * cone.cosAngleSqr;
// For a more sophisticated query, such as determining the interval of
// intersection of two height intervals [h0,h1] and [cone.hmin,cone.hmax],
// see IntrIntervals.ts, which supports semi-infinite intervals.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions. The C++ dimension template parameter N becomes a constructor
// argument; createMesh requires N = 3 and asserts it. The CreateMesh output
// reference parameters become a returned object literal.

import { GTE_C_HALF_PI, GTE_C_QUARTER_PI, GTE_C_TWO_PI } from './Constants';
import { logAssert } from './Logger';
import { Matrix, mulMatrix } from './Matrix';
import { Ray } from './Ray';
import { UniqueVerticesSimplices } from './UniqueVerticesSimplices';
import { Vector, add, div, hlift } from './Vector';
import { computeOrthogonalComplement3 } from './Vector3';

export class Cone {
    // The cone axis direction (ray.direction) must be unit length.
    ray: Ray;

    // The angle must be in (0,pi/2). The other members are derived from angle
    // to avoid calling trigonometric functions in geometric queries (for
    // speed). You may set the angle and compute these by calling
    // setAngle(inAngle).
    angle: number;
    cosAngle: number;
    sinAngle: number;
    tanAngle: number;
    cosAngleSqr: number;
    sinAngleSqr: number;
    invSinAngle: number;

    // The heights must satisfy 0 <= minHeight < maxHeight <= +infinity. For
    // an infinite cone, maxHeight is set to -1. For a finite cone, maxHeight
    // is set to a positive number. Be careful not to use maxHeight without
    // understanding this interpretation.
    private mMinHeight: number;
    private mMaxHeight: number;

    // Create an infinite cone with
    //   vertex = (0,...,0)
    //   axis = (0,...,0,1)
    //   angle = pi/4
    //   minimum height = 0
    //   maximum height = +infinity
    // The dimension N of the C++ template is a constructor argument here.
    constructor(n: number) {
        this.ray = new Ray(n);
        this.ray.origin.makeZero();
        this.ray.direction.makeUnit(n - 1);
        this.angle = 0;
        this.cosAngle = 0;
        this.sinAngle = 0;
        this.tanAngle = 0;
        this.cosAngleSqr = 0;
        this.sinAngleSqr = 0;
        this.invSinAngle = 0;
        this.mMinHeight = 0;
        this.mMaxHeight = -1;
        this.setAngle(GTE_C_QUARTER_PI);
        this.makeInfiniteCone();
    }

    // Create an infinite cone with the specified vertex, axis direction,
    // angle and with heights
    //   minimum height = 0
    //   maximum height = +infinity
    static fromRayAngle(inRay: Ray, inAngle: number): Cone {
        const cone = new Cone(inRay.dimension);
        cone.ray = inRay.clone();
        cone.setAngle(inAngle);
        cone.makeInfiniteCone();
        return cone;
    }

    // Create an infinite truncated cone with the specified vertex, axis
    // direction, angle and nonnegative minimum height. The maximum height is
    // +infinity. If you specify a minimum height of 0, you get the equivalent
    // of an infinite cone.
    static fromRayAngleMinHeight(inRay: Ray, inAngle: number,
        inMinHeight: number): Cone {
        const cone = new Cone(inRay.dimension);
        cone.ray = inRay.clone();
        cone.setAngle(inAngle);
        cone.makeInfiniteTruncatedCone(inMinHeight);
        return cone;
    }

    // Create a finite cone or a frustum of a cone with all parameters
    // specified. If you specify a minimum height of 0, you get a finite cone.
    // If you specify a positive minimum height, you get a frustum of a cone.
    static fromRayAngleMinMaxHeight(inRay: Ray, inAngle: number,
        inMinHeight: number, inMaxHeight: number): Cone {
        const cone = new Cone(inRay.dimension);
        cone.ray = inRay.clone();
        cone.setAngle(inAngle);
        cone.makeConeFrustum(inMinHeight, inMaxHeight);
        return cone;
    }

    // The dimension N of the space containing the cone.
    get dimension(): number {
        return this.ray.dimension;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Cone {
        const cone = new Cone(this.dimension);
        cone.ray = this.ray.clone();
        cone.angle = this.angle;
        cone.cosAngle = this.cosAngle;
        cone.sinAngle = this.sinAngle;
        cone.tanAngle = this.tanAngle;
        cone.cosAngleSqr = this.cosAngleSqr;
        cone.sinAngleSqr = this.sinAngleSqr;
        cone.invSinAngle = this.invSinAngle;
        cone.mMinHeight = this.mMinHeight;
        cone.mMaxHeight = this.mMaxHeight;
        return cone;
    }

    // The angle must be in (0,pi/2). The function sets 'angle' and computes
    // 'cosAngle', 'sinAngle', 'tanAngle', 'cosAngleSqr', 'sinAngleSqr' and
    // 'invSinAngle'.
    setAngle(inAngle: number): void {
        logAssert(0 < inAngle && inAngle < GTE_C_HALF_PI, 'Invalid angle.');
        this.angle = inAngle;
        this.cosAngle = Math.cos(this.angle);
        this.sinAngle = Math.sin(this.angle);
        this.tanAngle = Math.tan(this.angle);
        this.cosAngleSqr = this.cosAngle * this.cosAngle;
        this.sinAngleSqr = this.sinAngle * this.sinAngle;
        this.invSinAngle = 1 / this.sinAngle;
    }

    // Set the heights to obtain one of the four types of cones. Be aware that
    // an infinite cone has maxHeight set to -1. Be careful not to use
    // maxHeight without understanding this interpretation.
    makeInfiniteCone(): void {
        this.mMinHeight = 0;
        this.mMaxHeight = -1;
    }

    makeInfiniteTruncatedCone(inMinHeight: number): void {
        logAssert(inMinHeight >= 0, 'Invalid minimum height.');
        this.mMinHeight = inMinHeight;
        this.mMaxHeight = -1;
    }

    makeFiniteCone(inMaxHeight: number): void {
        logAssert(inMaxHeight > 0, 'Invalid maximum height.');
        this.mMinHeight = 0;
        this.mMaxHeight = inMaxHeight;
    }

    makeConeFrustum(inMinHeight: number, inMaxHeight: number): void {
        logAssert(inMinHeight >= 0 && inMaxHeight > inMinHeight,
            'Invalid minimum or maximum height.');
        this.mMinHeight = inMinHeight;
        this.mMaxHeight = inMaxHeight;
    }

    // Get the height extremes. For an infinite cone, maxHeight is -1. For a
    // finite cone, maxHeight is a positive number.
    getMinHeight(): number {
        return this.mMinHeight;
    }

    getMaxHeight(): number {
        return this.mMaxHeight;
    }

    heightInRange(h: number): boolean {
        return this.mMinHeight <= h &&
            (this.mMaxHeight !== -1 ? h <= this.mMaxHeight : true);
    }

    heightLessThanMin(h: number): boolean {
        return h < this.mMinHeight;
    }

    heightGreaterThanMax(h: number): boolean {
        return (this.mMaxHeight !== -1 ? h > this.mMaxHeight : false);
    }

    isFinite(): boolean {
        return this.mMaxHeight !== -1;
    }

    isInfinite(): boolean {
        return this.mMaxHeight === -1;
    }

    // Comparisons to support sorted containers. These are based only on
    // 'ray', 'angle', 'minHeight' and 'maxHeight'.
    equals(cone: Cone): boolean {
        return this.ray.equals(cone.ray)
            && this.angle === cone.angle
            && this.mMinHeight === cone.mMinHeight
            && this.mMaxHeight === cone.mMaxHeight;
    }

    notEquals(cone: Cone): boolean {
        return !this.equals(cone);
    }

    lessThan(cone: Cone): boolean {
        if (this.ray.lessThan(cone.ray)) {
            return true;
        }

        if (this.ray.greaterThan(cone.ray)) {
            return false;
        }

        if (this.angle < cone.angle) {
            return true;
        }

        if (this.angle > cone.angle) {
            return false;
        }

        if (this.mMinHeight < cone.mMinHeight) {
            return true;
        }

        if (this.mMinHeight > cone.mMinHeight) {
            return false;
        }

        return this.mMaxHeight < cone.mMaxHeight;
    }

    lessThanOrEqual(cone: Cone): boolean {
        return !cone.lessThan(this);
    }

    greaterThan(cone: Cone): boolean {
        return cone.lessThan(this);
    }

    greaterThanOrEqual(cone: Cone): boolean {
        return !this.lessThan(cone);
    }

    // Support for visualization. This requires a 3-dimensional cone. The
    // upstream output reference parameters become the returned object.
    createMesh(numMinVertices: number, inscribed: boolean):
        { vertices: Vector[], indices: number[] } {
        logAssert(this.dimension === 3,
            'Meshes can be generated only for 3-dimensional cones.');
        logAssert(this.isFinite(),
            'Meshes can be generated only for finite cones.');
        logAssert(numMinVertices >= 3,
            'At least three vertices are required.');

        const hMin = this.getMinHeight();
        const hMax = this.getMaxHeight();
        const rMin = hMin * this.tanAngle;
        const rMax = hMax * this.tanAngle;
        // TODO (upstream): The next line used to be 0.5*rMax/rMin-1, but if
        // hMin is zero, then rMin is zero and tNumExtra is 'inf'. How was this
        // equation derived? For now, guard against the division by zero using
        // 0.5*(1+rMax)/(1+rMin)-1.
        const tNumExtra = 0.5 * (1 + rMax) / (1 + rMin) - 1;
        let numExtra = 0;
        if (tNumExtra > 0) {
            numExtra = Math.ceil(tNumExtra);
        }
        const numMaxVertices = 2 * numMinVertices * (1 + numExtra);

        let polygonMin: Vector[];
        let polygonMax: Vector[];
        if (inscribed) {
            polygonMin = generateInscribed(numMinVertices, rMin);
            polygonMax = generateInscribed(numMaxVertices, rMax);
        }
        else {
            polygonMin = generateCircumscribed(numMinVertices, rMin);
            polygonMax = generateCircumscribed(numMaxVertices, rMax);
        }

        const { vertices, indices } = createConeFrustumMesh(numMinVertices,
            numMaxVertices, numExtra, hMin, hMax, polygonMin, polygonMax);

        // Transform to the coordinate system of the cone.
        const basis: Vector[] = [this.ray.direction.clone(), new Vector(3),
            new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const rotate = Matrix.zero(3, 3);
        rotate.setCol(0, basis[1]);
        rotate.setCol(1, basis[2]);
        rotate.setCol(2, basis[0]);
        for (let i = 0; i < vertices.length; ++i) {
            vertices[i] = add(mulMatrix(rotate, vertices[i]) as Vector,
                this.ray.origin);
        }

        return { vertices, indices };
    }
}

// The polygon of 'numVertices' vertices inscribed in the circle of the
// specified radius, with the first vertex repeated at the end.
function generateInscribed(numVertices: number, radius: number): Vector[] {
    const theta = GTE_C_TWO_PI / numVertices;
    const polygon: Vector[] = new Array<Vector>(numVertices + 1);
    for (let i = 0; i < numVertices; ++i) {
        const angle = i * theta;
        polygon[i] = Vector.fromArray([radius * Math.cos(angle),
            radius * Math.sin(angle)]);
    }
    polygon[numVertices] = polygon[0].clone();
    return polygon;
}

// The polygon of 'numVertices' vertices circumscribed about the circle of the
// specified radius, with the first vertex repeated at the end.
function generateCircumscribed(numVertices: number, radius: number): Vector[] {
    const theta = GTE_C_TWO_PI / numVertices;
    const inscribed: Vector[] = new Array<Vector>(numVertices + 1);
    for (let i = 0; i < numVertices; ++i) {
        const angle = i * theta;
        inscribed[i] = Vector.fromArray([radius * Math.cos(angle),
            radius * Math.sin(angle)]);
    }
    inscribed[numVertices] = inscribed[0].clone();

    const divisor = 1 + Math.cos(theta);
    const polygon: Vector[] = new Array<Vector>(numVertices + 1);
    for (let i = 0, ip1 = 1; i < numVertices; ++i, ++ip1) {
        polygon[i] = div(add(inscribed[i], inscribed[ip1]), divisor);
    }
    polygon[numVertices] = polygon[0].clone();
    return polygon;
}

// Build the triangle soup of the cone frustum in the cone's local coordinate
// system (the axis is the z-axis and the vertex is the origin), then convert
// it to an indexed mesh with unique vertices.
function createConeFrustumMesh(numMinVertices: number, numMaxVertices: number,
    numExtra: number, hMin: number, hMax: number,
    polygonMin: readonly Vector[], polygonMax: readonly Vector[]):
    { vertices: Vector[], indices: number[] } {
    const vertexPool: Vector[] = [];
    let V0: Vector, V1: Vector, V2: Vector;
    for (let i0 = 0, i1 = 1; i0 < numMinVertices; i0 = i1++) {
        const j0 = 2 * (numExtra + 1) * i0;
        V0 = hlift(polygonMin[i0], hMin);
        for (let k0 = 0, k1 = 1; k0 <= numExtra; k0 = k1++) {
            V1 = hlift(polygonMax[j0 + k1], hMax);
            V2 = hlift(polygonMax[j0 + k0], hMax);
            vertexPool.push(V0.clone(), V1, V2);
        }

        const j1 = 2 * (numExtra + 1) * i1;
        V0 = hlift(polygonMin[i1], hMin);
        for (let k0 = 0, k1 = 1; k0 <= numExtra; k0 = k1++) {
            V1 = hlift(polygonMax[j1 - k0], hMax);
            V2 = hlift(polygonMax[j1 - k1], hMax);
            vertexPool.push(V0.clone(), V1, V2);
        }

        const jmid = j0 + (numExtra + 1);
        V0 = hlift(polygonMax[jmid], hMax);
        V1 = hlift(polygonMin[i0], hMin);
        V2 = hlift(polygonMin[i1], hMin);
        vertexPool.push(V0, V1, V2);
    }

    // The disk at h-minimum.
    V0 = Vector.fromArray([0, 0, hMin]);
    for (let i0 = 0, i1 = 1; i0 < numMinVertices; i0 = i1++) {
        V1 = hlift(polygonMin[i1], hMin);
        V2 = hlift(polygonMin[i0], hMin);
        vertexPool.push(V0.clone(), V1, V2);
    }

    // The disk at h-maximum.
    V0 = Vector.fromArray([0, 0, hMax]);
    for (let i0 = 0, i1 = 1; i0 < numMaxVertices; i0 = i1++) {
        V1 = hlift(polygonMax[i0], hMax);
        V2 = hlift(polygonMax[i1], hMax);
        vertexPool.push(V0.clone(), V1, V2);
    }

    const uvs = new UniqueVerticesSimplices<Vector>(3);
    return uvs.generateIndexedSimplices(vertexPool);
}

// Alias for convenience (the port of the upstream template alias).
export type Cone3 = Cone;
