import { describe, it, expect } from 'vitest';
import {
    PointInPolyhedron3, PointInPolyhedron3Face, PointInPolyhedron3FaceType
} from '../src/ContPointInPolyhedron3';
import { Hyperplane } from '../src/Hyperplane';
import type { Plane3 } from '../src/Hyperplane';
import { Vector, dot, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

// A set of pseudo-random unit-length directions; the count must be odd so
// that the majority vote is decisive.
function randomDirections(count: number, seed: number): Vector[] {
    const rand = makeRandom(seed);
    const out: Vector[] = [];
    while (out.length < count) {
        const d = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
        const length = Math.hypot(d.get(0), d.get(1), d.get(2));
        if (length < 0.25) {
            continue;
        }
        out.push(Vector.fromArray(d.values.map(c => c / length)));
    }
    return out;
}

// The plane of a face whose vertices are counterclockwise when viewed from
// outside the polyhedron; the resulting normal points outside.
function facePlane(points: readonly Vector[],
    indices: readonly number[]): Plane3 {
    const V0 = points[indices[0]];
    const V1 = points[indices[1]];
    const V2 = points[indices[2]];
    const n = cross(sub(V1, V0), sub(V2, V0));
    const length = Math.hypot(n.get(0), n.get(1), n.get(2));
    const normal = Vector.fromArray(n.values.map(c => c / length));
    return Hyperplane.fromNormalOrigin(normal, V0);
}

// A fan triangulation of a face, from its first vertex.
function fanTriangles(indices: readonly number[]): number[] {
    const out: number[] = [];
    for (let k = 1; k + 1 < indices.length; ++k) {
        out.push(indices[0], indices[k], indices[k + 1]);
    }
    return out;
}

function makeFace(points: readonly Vector[], indices: readonly number[],
    triangles?: readonly number[]): PointInPolyhedron3Face {
    return PointInPolyhedron3Face.fromIndicesPlane(indices,
        facePlane(points, indices),
        triangles === undefined ? fanTriangles(indices) : triangles);
}

// ---------------------------------------------------------------------------
// The unit cube [0,1]^3, with the six quadrilateral faces counterclockwise
// when viewed from outside.
// ---------------------------------------------------------------------------
const cubePoints: Vector[] = [
    v3(0, 0, 0), v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0),
    v3(0, 0, 1), v3(1, 0, 1), v3(1, 1, 1), v3(0, 1, 1)
];
const cubeQuads: number[][] = [
    [0, 3, 2, 1],   // z = 0, normal (0,0,-1)
    [4, 5, 6, 7],   // z = 1, normal (0,0,+1)
    [0, 1, 5, 4],   // y = 0, normal (0,-1,0)
    [1, 2, 6, 5],   // x = 1, normal (+1,0,0)
    [2, 3, 7, 6],   // y = 1, normal (0,+1,0)
    [3, 0, 4, 7]    // x = 0, normal (-1,0,0)
];
const cubeConvexFaces = cubeQuads.map(q => makeFace(cubePoints, q));
const cubeTriangleFaces: PointInPolyhedron3Face[] = [];
for (const q of cubeQuads) {
    cubeTriangleFaces.push(makeFace(cubePoints, [q[0], q[1], q[2]]));
    cubeTriangleFaces.push(makeFace(cubePoints, [q[0], q[2], q[3]]));
}

function cubeContains(p: Vector): boolean {
    for (let k = 0; k < 3; ++k) {
        if (p.get(k) <= 0 || p.get(k) >= 1) {
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// The tetrahedron with vertices (0,0,0), (1,0,0), (0,1,0), (0,0,1).
// ---------------------------------------------------------------------------
const tetraPoints: Vector[] = [
    v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)
];
const tetraIndices: number[][] = [
    [0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]
];
const tetraFaces = tetraIndices.map(f => makeFace(tetraPoints, f));

function tetraContains(p: Vector): boolean {
    const x = p.get(0), y = p.get(1), z = p.get(2);
    return x > 0 && y > 0 && z > 0 && x + y + z < 1;
}

// ---------------------------------------------------------------------------
// A non-convex L-shaped prism: the base polygon
//   (0,0), (2,0), (2,1), (1,1), (1,2), (0,2)
// in the plane z = 0 (counterclockwise in xy), extruded to z = 1. The two
// hexagonal faces are simple but not convex.
// ---------------------------------------------------------------------------
const lBase: [number, number][] = [
    [0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]
];
const lPoints: Vector[] = [
    ...lBase.map(([x, y]) => v3(x, y, 0)),
    ...lBase.map(([x, y]) => v3(x, y, 1))
];
// The bottom face (normal (0,0,-1)) is the base polygon reversed; a fan from
// the base vertex 0 triangulates the L exactly (verified by area: 1 + 1/2 +
// 1/2 + 1 = 3, the area of the L).
const lBottomIndices = [5, 4, 3, 2, 1, 0];
const lBottomTriangles = fanTriangles([0, 1, 2, 3, 4, 5]);
const lTopIndices = [6, 7, 8, 9, 10, 11];
const lFaces: PointInPolyhedron3Face[] = [
    makeFace(lPoints, lBottomIndices, lBottomTriangles),
    makeFace(lPoints, lTopIndices, fanTriangles(lTopIndices))
];
for (let i = 0; i < 6; ++i) {
    const j = (i + 1) % 6;
    lFaces.push(makeFace(lPoints, [i, j, j + 6, i + 6]));
}

function lContains(p: Vector): boolean {
    const x = p.get(0), y = p.get(1), z = p.get(2);
    if (z <= 0 || z >= 1) {
        return false;
    }
    const inArm0 = x > 0 && x < 2 && y > 0 && y < 1;
    const inArm1 = x > 0 && x < 1 && y > 0 && y < 2;
    return inArm0 || inArm1;
}

// ---------------------------------------------------------------------------

const dirs7 = randomDirections(7, 12345);
const dirs1 = [Vector.fromArray([1, 2, 3].map(c => c / Math.sqrt(14)))];

describe('PointInPolyhedron3Face', () => {
    it('default constructs an empty face with the default plane', () => {
        const face = new PointInPolyhedron3Face();
        expect(face.indices).toEqual([]);
        expect(face.triangles).toEqual([]);
        expect(face.plane.normal.values).toEqual([0, 0, 1]);
        expect(face.plane.constant).toBe(0);
    });

    it('copies the arrays and the plane in fromIndicesPlane', () => {
        const indices = [1, 2, 3];
        const triangles = [1, 2, 3];
        const plane = Hyperplane.fromNormalOrigin(v3(0, 0, 1), v3(0, 0, 5));
        const face = PointInPolyhedron3Face.fromIndicesPlane(indices, plane,
            triangles);
        indices[0] = 99;
        triangles[0] = 99;
        plane.constant = -1;
        expect(face.indices).toEqual([1, 2, 3]);
        expect(face.triangles).toEqual([1, 2, 3]);
        expect(face.plane.constant).toBe(5);
    });

    it('holds more than three indices (upstream fixes a fixed-size array)', () => {
        const face = makeFace(cubePoints, cubeQuads[0]);
        expect(face.indices.length).toBe(4);
        expect(face.plane.normal.values).toEqual([0, 0, -1]);
        expect(face.plane.constant).toBeCloseTo(0, 12);
    });
});

describe('PointInPolyhedron3 with triangle faces', () => {
    const cube = new PointInPolyhedron3(PointInPolyhedron3FaceType.TRIANGLE,
        cubePoints, cubeTriangleFaces, dirs7, 0);
    const tetra = new PointInPolyhedron3(PointInPolyhedron3FaceType.TRIANGLE,
        tetraPoints, tetraFaces, dirs7, 0);

    it('classifies known interior and exterior points of a cube', () => {
        expect(cube.contains(v3(0.5, 0.5, 0.5))).toBe(true);
        expect(cube.contains(v3(0.01, 0.99, 0.5))).toBe(true);
        expect(cube.contains(v3(1.5, 0.5, 0.5))).toBe(false);
        expect(cube.contains(v3(-0.5, 0.5, 0.5))).toBe(false);
        expect(cube.contains(v3(0.5, 0.5, 10))).toBe(false);
        expect(cube.contains(v3(5, 5, 5))).toBe(false);
    });

    it('classifies known interior and exterior points of a tetrahedron', () => {
        expect(tetra.contains(v3(0.1, 0.1, 0.1))).toBe(true);
        expect(tetra.contains(v3(0.2, 0.2, 0.5))).toBe(true);
        expect(tetra.contains(v3(0.4, 0.4, 0.4))).toBe(false);
        expect(tetra.contains(v3(-0.1, 0.1, 0.1))).toBe(false);
        expect(tetra.contains(v3(1, 1, 1))).toBe(false);
    });

    it('reports the number of rays it casts', () => {
        expect(cube.numRays).toBe(7);
        expect(new PointInPolyhedron3(PointInPolyhedron3FaceType.TRIANGLE,
            cubePoints, cubeTriangleFaces, dirs1, 0).numRays).toBe(1);
    });

    it('agrees with the analytic test on randomized points', () => {
        const rand = makeRandom(8080);
        for (let trial = 0; trial < 300; ++trial) {
            const p = v3(3 * rand() - 1, 3 * rand() - 1, 3 * rand() - 1);
            expect(cube.contains(p)).toBe(cubeContains(p));
            expect(tetra.contains(p)).toBe(tetraContains(p));
        }
    });

    it('survives grazing rays through vertices and edges by majority vote', () => {
        // From the cube center, (1,1,1)/sqrt(3) exits exactly through the
        // vertex (1,1,1) and (1,1,0)/sqrt(2) exits through an edge; a single
        // such ray can miscount, so the vote is over a set of rays that
        // includes them.
        const s3 = 1 / Math.sqrt(3);
        const s2 = 1 / Math.sqrt(2);
        const grazing = [
            v3(s3, s3, s3), v3(s2, s2, 0), v3(-s2, 0, s2),
            ...randomDirections(4, 999)
        ];
        const grazingCube = new PointInPolyhedron3(
            PointInPolyhedron3FaceType.TRIANGLE, cubePoints, cubeTriangleFaces,
            grazing, 0);
        expect(grazingCube.contains(v3(0.5, 0.5, 0.5))).toBe(true);
        expect(grazingCube.contains(v3(0.25, 0.25, 0.25))).toBe(true);
        expect(grazingCube.contains(v3(2, 2, 2))).toBe(false);
        expect(grazingCube.contains(v3(-1, 0.5, 0.5))).toBe(false);
        expect(grazingCube.contains(v3(0.5, 0.5, 3))).toBe(false);
    });

    it('is unaffected by the ray count when the rays are generic', () => {
        for (const count of [1, 3, 5, 11]) {
            const query = new PointInPolyhedron3(
                PointInPolyhedron3FaceType.TRIANGLE, cubePoints,
                cubeTriangleFaces, randomDirections(count, 31 * count), 0);
            expect(query.contains(v3(0.5, 0.5, 0.5))).toBe(true);
            expect(query.contains(v3(1.5, 0.5, 0.5))).toBe(false);
        }
    });
});

describe('PointInPolyhedron3 with convex faces', () => {
    // Method 0 uses a triangle fan; methods 1 and 2 project onto the plane of
    // the face and run the O(N) and O(log N) convex point-in-polygon tests.
    // All three must agree, and all three need the full four indices of each
    // quadrilateral face.
    const queries = [0, 1, 2].map(method => new PointInPolyhedron3(
        PointInPolyhedron3FaceType.CONVEX, cubePoints, cubeConvexFaces, dirs7,
        method));

    it('classifies known points the same way for methods 0, 1 and 2', () => {
        const inside = [v3(0.5, 0.5, 0.5), v3(0.02, 0.5, 0.97),
            v3(0.9, 0.1, 0.5)];
        const outside = [v3(1.5, 0.5, 0.5), v3(-2, 0.5, 0.5),
            v3(0.5, 0.5, 1.001), v3(3, 3, 3)];
        for (const query of queries) {
            for (const p of inside) {
                expect(query.contains(p)).toBe(true);
            }
            for (const p of outside) {
                expect(query.contains(p)).toBe(false);
            }
        }
    });

    it('agrees with the analytic test on randomized points for all methods', () => {
        const rand = makeRandom(2718);
        for (let trial = 0; trial < 150; ++trial) {
            const p = v3(3 * rand() - 1, 3 * rand() - 1, 3 * rand() - 1);
            const expected = cubeContains(p);
            for (const query of queries) {
                expect(query.contains(p)).toBe(expected);
            }
        }
    });

    it('handles a convex polyhedron with a triangular and a quad face', () => {
        // A square pyramid: a quadrilateral base and four triangular sides.
        const points = [
            v3(-1, -1, 0), v3(1, -1, 0), v3(1, 1, 0), v3(-1, 1, 0),
            v3(0, 0, 2)
        ];
        const faces = [
            makeFace(points, [0, 3, 2, 1]),
            makeFace(points, [0, 1, 4]),
            makeFace(points, [1, 2, 4]),
            makeFace(points, [2, 3, 4]),
            makeFace(points, [3, 0, 4])
        ];
        for (const method of [0, 1, 2]) {
            const query = new PointInPolyhedron3(
                PointInPolyhedron3FaceType.CONVEX, points, faces, dirs7,
                method);
            expect(query.contains(v3(0, 0, 0.5))).toBe(true);
            expect(query.contains(v3(0.4, 0.4, 1))).toBe(true);
            expect(query.contains(v3(0, 0, 2.5))).toBe(false);
            expect(query.contains(v3(0.9, 0.9, 1.5))).toBe(false);
            expect(query.contains(v3(0, 0, -0.5))).toBe(false);
        }
    });
});

describe('PointInPolyhedron3 with simple faces', () => {
    // Method 0 iterates over the triangulation of each face; method 1
    // projects onto the plane of the face and runs the O(N) simple-polygon
    // test. The L-shaped prism has two non-convex hexagonal faces.
    const queries = [0, 1].map(method => new PointInPolyhedron3(
        PointInPolyhedron3FaceType.SIMPLE, lPoints, lFaces, dirs7, method));

    it('classifies known points of the L-shaped prism', () => {
        // In the long arm, in the tall arm and in the shared corner.
        const inside = [v3(1.5, 0.5, 0.5), v3(0.5, 1.5, 0.5),
            v3(0.5, 0.5, 0.5), v3(1.9, 0.1, 0.9)];
        // The notch of the L is outside, as are points beyond every face.
        const outside = [v3(1.5, 1.5, 0.5), v3(1.9, 1.9, 0.5),
            v3(1.5, 0.5, 1.5), v3(-1, 0.5, 0.5), v3(3, 3, 3)];
        for (const query of queries) {
            for (const p of inside) {
                expect(query.contains(p)).toBe(true);
            }
            for (const p of outside) {
                expect(query.contains(p)).toBe(false);
            }
        }
    });

    it('agrees with the analytic test on randomized points for both methods', () => {
        const rand = makeRandom(161803);
        for (let trial = 0; trial < 150; ++trial) {
            const p = v3(4 * rand() - 1, 4 * rand() - 1, 3 * rand() - 1);
            const expected = lContains(p);
            for (const query of queries) {
                expect(query.contains(p)).toBe(expected);
            }
        }
    });

    it('agrees with the triangle-face query on the same triangulated prism', () => {
        // Build the same solid as a pure triangle mesh and compare.
        const triFaces: PointInPolyhedron3Face[] = [];
        for (const face of lFaces) {
            for (let t = 0; t < face.triangles.length; t += 3) {
                triFaces.push(makeFace(lPoints, [face.triangles[t],
                    face.triangles[t + 1], face.triangles[t + 2]]));
            }
        }
        const triQuery = new PointInPolyhedron3(
            PointInPolyhedron3FaceType.TRIANGLE, lPoints, triFaces, dirs7, 0);
        const rand = makeRandom(577);
        for (let trial = 0; trial < 120; ++trial) {
            const p = v3(4 * rand() - 1, 4 * rand() - 1, 3 * rand() - 1);
            expect(triQuery.contains(p)).toBe(lContains(p));
            expect(queries[0].contains(p)).toBe(triQuery.contains(p));
        }
    });

    it('requires a triangulation for method 0', () => {
        const facesNoTriangles = lFaces.map(f =>
            PointInPolyhedron3Face.fromIndicesPlane(f.indices, f.plane, []));
        const query = new PointInPolyhedron3(
            PointInPolyhedron3FaceType.SIMPLE, lPoints, facesNoTriangles,
            dirs7, 0);
        expect(() => query.contains(v3(0.5, 0.5, 0.5)))
            .toThrow('Triangulation must exist.');

        // Method 1 does not use the triangulation.
        const query1 = new PointInPolyhedron3(
            PointInPolyhedron3FaceType.SIMPLE, lPoints, facesNoTriangles,
            dirs7, 1);
        expect(query1.contains(v3(0.5, 0.5, 0.5))).toBe(true);
        expect(query1.contains(v3(1.5, 1.5, 0.5))).toBe(false);
    });

    it('reports false for an unsupported simple-face method, as upstream', () => {
        const query = new PointInPolyhedron3(
            PointInPolyhedron3FaceType.SIMPLE, lPoints, lFaces, dirs7, 2);
        expect(query.contains(v3(0.5, 0.5, 0.5))).toBe(false);
    });
});

describe('PointInPolyhedron3 face-plane culling', () => {
    it('keeps the outward normals consistent with the built faces', () => {
        // Sanity check on the test fixtures themselves: every face normal
        // points away from an interior point of the solid.
        const check = (interior: Vector,
            faces: readonly PointInPolyhedron3Face[]): void => {
            for (const face of faces) {
                const signed = dot(face.plane.normal, interior) -
                    face.plane.constant;
                expect(signed).toBeLessThan(0);
            }
        };
        check(v3(0.5, 0.5, 0.5), cubeConvexFaces);
        check(v3(0.5, 0.5, 0.5), cubeTriangleFaces);
        check(v3(0.2, 0.2, 0.2), tetraFaces);
        check(v3(0.5, 0.5, 0.5), lFaces);
    });

    it('classifies points on and just off a face of the cube', () => {
        const query = new PointInPolyhedron3(
            PointInPolyhedron3FaceType.TRIANGLE, cubePoints, cubeTriangleFaces,
            dirs7, 0);
        // Just inside and just outside the x = 1 face.
        expect(query.contains(v3(1 - 1e-8, 0.5, 0.5))).toBe(true);
        expect(query.contains(v3(1 + 1e-8, 0.5, 0.5))).toBe(false);
        // Just inside and just outside the z = 0 face.
        expect(query.contains(v3(0.5, 0.5, 1e-8))).toBe(true);
        expect(query.contains(v3(0.5, 0.5, -1e-8))).toBe(false);
    });
});
