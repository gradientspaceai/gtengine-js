// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContPointInPolyhedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// This class contains various implementations for point-in-polyhedron
// queries. The planes stored with the faces are used in all cases to reject
// ray-face intersection tests, a quick culling operation.
//
// The algorithm is to cast a ray from the input point P and test for
// intersection against each face of the polyhedron. If the ray only
// intersects faces at interior points (not vertices, not edge points), then
// the point is inside when the number of intersections is odd and the point
// is outside when the number of intersections is even. If the ray intersects
// an edge or a vertex, then the counting must be handled differently. The
// details are tedious. As an alternative, the approach here is to allow you
// to specify 2*N+1 rays, where N >= 0. You should choose these rays randomly.
// Each ray reports "inside" or "outside". Whichever result occurs N+1 or more
// times is the "winner". The number of directions passed to the constructor
// is 2*N+1. If you are feeling lucky, choose a single direction.
//
// Port notes:
//   - The nested enum 'FaceType' becomes the exported, file-qualified enum
//     'PointInPolyhedron3FaceType'; the nested class 'Face' becomes the
//     exported class 'PointInPolyhedron3Face' (a nested class is not
//     expressible in TypeScript).
//   - The (count, pointer) pairs of the constructor become arrays; the
//     'numPoints' parameter has no counterpart because the upstream member
//     'mNumPoints' is never read. The arrays are stored by reference, which
//     has the same aliasing behavior as the upstream 'const*' members, so the
//     caller's data must stay alive and unmodified for the lifetime of the
//     object.
//   - 'Contains' becomes 'contains'.
//   - The upstream 'ContainsFunction' function-pointer type becomes the
//     exported type 'PointInPolyhedron3ContainsFunction'; its 'bool& odd'
//     output parameter becomes a boolean return value (the new value of
//     'odd').

import { PointInPolygon2 } from './ContPointInPolygon2.js';
import { Hyperplane, type Plane3 } from './Hyperplane.js';
import { IntrRay3Plane3FI } from './IntrRay3Plane3.js';
import { IntrRay3Triangle3TI } from './IntrRay3Triangle3.js';
import { logAssert } from './Logger.js';
import { Ray } from './Ray.js';
import { Triangle } from './Triangle.js';
import { Vector, dot, sub } from './Vector.js';
import { computeOrthogonalComplement3 } from './Vector3.js';

export enum PointInPolyhedron3FaceType {
    TRIANGLE = 0,
    CONVEX = 1,
    SIMPLE = 2
}

export class PointInPolyhedron3Face {
    // The members 'indices' and 'plane' are used for triangle faces, for
    // convex polygon faces, and for simple polygon faces. The member
    // 'triangles' is used only for simple faces that are not convex.
    //
    // When you view the face from outside, the vertices are counterclockwise
    // ordered. The indices array stores the indices into the vertex array.
    //
    // NOTE (upstream bug): the upstream member is
    // 'std::array<int32_t, 3> indices', a fixed size of 3, yet the CONVEX and
    // SIMPLE code paths call 'indices.size()' expecting the face's actual
    // vertex count (older GTE releases declared the member as
    // 'std::vector<int32_t>'). As written upstream, a convex or simple face
    // with more than three vertices cannot be represented: ContainsC0 fans
    // only the first triangle and SharedContains projects only the first
    // three vertices, so the point-in-polygon test runs against a truncated
    // polygon and reports wrong containment. The port restores the
    // variable-length array.
    indices: number[];

    // The normal vector is unit length and points to the outside of the
    // polyhedron.
    plane: Plane3;

    // Each simple face may be triangulated. The indices are relative to the
    // vertex array. Each triple of indices represents a triangle in the
    // triangulation.
    triangles: number[];

    constructor() {
        this.indices = [];
        this.plane = new Hyperplane(3);
        this.triangles = [];
    }

    // The port of aggregate initialization of a Face. The arrays and the
    // plane are copied, matching C++ value semantics.
    static fromIndicesPlane(indices: readonly number[], plane: Plane3,
        triangles: readonly number[] = []): PointInPolyhedron3Face {
        const face = new PointInPolyhedron3Face();
        face.indices = indices.slice();
        face.plane = plane.clone();
        face.triangles = triangles.slice();
        return face;
    }
}

// The port of the upstream 'ContainsFunction' typedef. The 'bool& odd' output
// parameter is returned instead.
export type PointInPolyhedron3ContainsFunction = (
    method: number,
    PIP: PointInPolygon2,
    projIntersect: Vector,
    odd: boolean) => boolean;

function containsPointConvex(method: number, PIP: PointInPolygon2,
    projIntersect: Vector, odd: boolean): boolean {
    if (method === 1) {
        if (PIP.containsConvexOrderN(projIntersect)) {
            // The ray intersects the face.
            return !odd;
        }
    } else {
        if (PIP.containsConvexOrderLogN(projIntersect)) {
            // The ray intersects the face.
            return !odd;
        }
    }
    return odd;
}

function containsPointSimple(_method: number, PIP: PointInPolygon2,
    projIntersect: Vector, odd: boolean): boolean {
    if (PIP.contains(projIntersect)) {
        // The ray intersects the face.
        return !odd;
    }
    return odd;
}

// For all types of faces. The ray origin is the test point. The ray direction
// is one of those passed to the constructor. The plane normal is a
// unit-length normal to the face that points outside the polyhedron.
function fastNoIntersect(ray: Ray, plane: Plane3): boolean {
    const planeDistance = dot(plane.normal, ray.origin) - plane.constant;
    const planeAngle = dot(plane.normal, ray.direction);

    if (planeDistance < 0) {
        // The ray origin is on the negative side of the plane.
        if (planeAngle <= 0) {
            // The ray points away from the plane.
            return true;
        }
    }

    if (planeDistance > 0) {
        // The ray origin is on the positive side of the plane.
        if (planeAngle >= 0) {
            // The ray points away from the plane.
            return true;
        }
    }

    return false;
}

export class PointInPolyhedron3 {
    private readonly mType: PointInPolyhedron3FaceType;
    private readonly mPoints: readonly Vector[];
    private readonly mFaces: readonly PointInPolyhedron3Face[];
    private readonly mMethod: number;
    private readonly mDirections: readonly Vector[];

    // The 'contains' query uses the 'method' parameter to determine which
    // point-in-polygon test is used.
    //
    // Triangle faces.
    //   0 : The parameter is unused. Ray-triangle tests are performed without
    //       projection on the planes of the triangles.
    //
    // Convex faces.
    //   0 : Use a triangle fan and perform a ray-triangle intersection query
    //       for each triangle.
    //   1 : Find the point of intersection of ray and plane of polygon. Test
    //       whether that point is inside the convex polygon using an O(N)
    //       test.
    //   2 : Find the point of intersection of ray and plane of polygon. Test
    //       whether that point is inside the convex polygon using an
    //       O(log N) test.
    //
    // Simple faces that are not convex.
    //   0 : Iterate over the triangles of each face and perform a
    //       ray-triangle intersection query for each triangle. This requires
    //       that the Face.triangles array be initialized for each face.
    //   1 : Find the point of intersection of ray and plane of polygon. Test
    //       whether that point is inside the polygon using an O(N) test. The
    //       Face.triangles array is not used for this method, so it does not
    //       have to be initialized for each face.
    constructor(type: PointInPolyhedron3FaceType, points: readonly Vector[],
        faces: readonly PointInPolyhedron3Face[],
        directions: readonly Vector[], method: number) {
        this.mType = type;
        this.mPoints = points;
        this.mFaces = faces;
        this.mMethod = method;
        this.mDirections = directions;
    }

    // The number of rays cast by the query; this is the number of directions
    // passed to the constructor.
    get numRays(): number {
        return this.mDirections.length;
    }

    // This function selects the actual algorithm based on the type of face
    // selected in the constructor. As upstream, an unsupported (type, method)
    // pair reports 'false'.
    contains(p: Vector): boolean {
        if (this.mType === PointInPolyhedron3FaceType.TRIANGLE) {
            return this.containsT0(p);
        }

        if (this.mType === PointInPolyhedron3FaceType.CONVEX) {
            if (this.mMethod === 0) {
                return this.containsC0(p);
            } else {
                // mMethod is 1 or 2
                return this.containsC1C2(p, this.mMethod);
            }
        }

        if (this.mType === PointInPolyhedron3FaceType.SIMPLE) {
            if (this.mMethod === 0) {
                return this.containsS0(p);
            }

            if (this.mMethod === 1) {
                return this.containsS1(p, this.mMethod);
            }
        }

        return false;
    }

    // The majority vote over the rays. The upstream comparison
    // 'insideCount > mNumRays / 2' uses C++ integer division.
    private isInside(insideCount: number): boolean {
        return insideCount > Math.floor(this.mDirections.length / 2);
    }

    // For triangle faces.
    private containsT0(p: Vector): boolean {
        let insideCount = 0;

        const rtQuery = new IntrRay3Triangle3TI();
        const triangle = new Triangle(3);
        const ray = new Ray(3);
        ray.origin = p.clone();

        for (let j = 0; j < this.mDirections.length; ++j) {
            ray.direction = this.mDirections[j].clone();

            // Zero intersections to start with.
            let odd = false;

            for (let i = 0; i < this.mFaces.length; ++i) {
                const face = this.mFaces[i];

                // Attempt to quickly cull the triangle.
                if (fastNoIntersect(ray, face.plane)) {
                    continue;
                }

                // Get the triangle vertices.
                for (let k = 0; k < 3; ++k) {
                    triangle.v[k] = this.mPoints[face.indices[k]];
                }

                // Test for intersection.
                if (rtQuery.test(ray, triangle).intersect) {
                    // The ray intersects the triangle.
                    odd = !odd;
                }
            }

            if (odd) {
                insideCount++;
            }
        }

        return this.isInside(insideCount);
    }

    // For convex faces.
    private containsC0(p: Vector): boolean {
        let insideCount = 0;

        const rtQuery = new IntrRay3Triangle3TI();
        const triangle = new Triangle(3);
        const ray = new Ray(3);
        ray.origin = p.clone();

        for (let j = 0; j < this.mDirections.length; ++j) {
            ray.direction = this.mDirections[j].clone();

            // Zero intersections to start with.
            let odd = false;

            for (let i = 0; i < this.mFaces.length; ++i) {
                const face = this.mFaces[i];

                // Attempt to quickly cull the face.
                if (fastNoIntersect(ray, face.plane)) {
                    continue;
                }

                // Process the triangles in a trifan of the face.
                const numVerticesM1 = face.indices.length - 1;
                triangle.v[0] = this.mPoints[face.indices[0]];
                for (let k = 1; k < numVerticesM1; ++k) {
                    triangle.v[1] = this.mPoints[face.indices[k]];
                    triangle.v[2] = this.mPoints[face.indices[k + 1]];

                    if (rtQuery.test(ray, triangle).intersect) {
                        // The ray intersects the triangle.
                        odd = !odd;
                    }
                }
            }

            if (odd) {
                insideCount++;
            }
        }

        return this.isInside(insideCount);
    }

    private containsC1C2(p: Vector, method: number): boolean {
        return this.sharedContains(p, method, containsPointConvex);
    }

    // For simple faces.
    private containsS0(p: Vector): boolean {
        let insideCount = 0;

        const rtQuery = new IntrRay3Triangle3TI();
        const triangle = new Triangle(3);
        const ray = new Ray(3);
        ray.origin = p.clone();

        for (let j = 0; j < this.mDirections.length; ++j) {
            ray.direction = this.mDirections[j].clone();

            // Zero intersections to start with.
            let odd = false;

            for (let i = 0; i < this.mFaces.length; ++i) {
                const face = this.mFaces[i];

                // Attempt to quickly cull the face.
                if (fastNoIntersect(ray, face.plane)) {
                    continue;
                }

                // The triangulation must exist to use it.
                const numTriangles = Math.floor(face.triangles.length / 3);
                logAssert(numTriangles > 0, 'Triangulation must exist.');

                // Process the triangles in a triangulation of the face.
                let currIndex = 0;
                for (let t = 0; t < numTriangles; ++t) {
                    // Get the triangle vertices.
                    for (let k = 0; k < 3; ++k) {
                        triangle.v[k] = this.mPoints[face.triangles[currIndex++]];
                    }

                    // Test for intersection.
                    if (rtQuery.test(ray, triangle).intersect) {
                        // The ray intersects the triangle.
                        odd = !odd;
                    }
                }
            }

            if (odd) {
                insideCount++;
            }
        }

        return this.isInside(insideCount);
    }

    private containsS1(p: Vector, method: number): boolean {
        return this.sharedContains(p, method, containsPointSimple);
    }

    // Shared code for containsC1C2 and containsS1.
    private sharedContains(p: Vector, method: number,
        containsPoint: PointInPolyhedron3ContainsFunction): boolean {
        let insideCount = 0;

        const rpQuery = new IntrRay3Plane3FI();
        const ray = new Ray(3);
        ray.origin = p.clone();

        for (let j = 0; j < this.mDirections.length; ++j) {
            ray.direction = this.mDirections[j].clone();

            // Zero intersections to start with.
            let odd = false;

            for (let i = 0; i < this.mFaces.length; ++i) {
                const face = this.mFaces[i];

                // Attempt to quickly cull the face.
                if (fastNoIntersect(ray, face.plane)) {
                    continue;
                }

                // Compute the ray-plane intersection.
                const result = rpQuery.find(ray, face.plane);

                // If you trigger this assertion, numerical round-off errors
                // have led to a discrepancy between fastNoIntersect and the
                // find() result.
                logAssert(result.intersect, 'Unexpected condition.');

                // Get a coordinate system for the plane. Use vertex 0 as the
                // origin.
                const V0 = this.mPoints[face.indices[0]];
                const basis: Vector[] = [
                    face.plane.normal.clone(), new Vector(3), new Vector(3)
                ];
                computeOrthogonalComplement3(1, basis);

                // Project the intersection onto the plane.
                let diff = sub(result.point, V0);
                const projIntersect = Vector.fromArray([
                    dot(basis[1], diff), dot(basis[2], diff)
                ]);

                // Project the face vertices onto the plane of the face.
                // Vertex 0 is always the origin.
                const numIndices = face.indices.length;
                const projVertices: Vector[] = new Array<Vector>(numIndices);
                projVertices[0] = Vector.zero(2);
                for (let k = 1; k < numIndices; ++k) {
                    diff = sub(this.mPoints[face.indices[k]], V0);
                    projVertices[k] = Vector.fromArray([
                        dot(basis[1], diff), dot(basis[2], diff)
                    ]);
                }

                // Test whether the intersection point is in the polygon.
                const PIP = new PointInPolygon2(projVertices);
                odd = containsPoint(method, PIP, projIntersect, odd);
            }

            if (odd) {
                insideCount++;
            }
        }

        return this.isInside(insideCount);
    }
}
