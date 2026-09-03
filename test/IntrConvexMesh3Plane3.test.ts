import { describe, it, expect } from 'vitest';
import { ConvexMesh3 } from '../src/ConvexMesh3.js';
import type { ConvexMesh3Triangle } from '../src/ConvexMesh3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    IntrConvexMesh3Plane3FI,
    defaultIntrConvexMesh3Plane3FIResult
} from '../src/IntrConvexMesh3Plane3.js';

const Q = IntrConvexMesh3Plane3FI;

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The unit cube [0,1]^3 with outward-facing (counterclockwise from outside)
// triangle faces.
function unitCube(): ConvexMesh3 {
    const mesh = new ConvexMesh3();
    mesh.configuration = ConvexMesh3.CFG_POLYHEDRON;
    mesh.vertices = [
        v3(0, 0, 0), v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0),
        v3(0, 0, 1), v3(1, 0, 1), v3(1, 1, 1), v3(0, 1, 1)
    ];
    mesh.triangles = [
        [0, 3, 2], [0, 2, 1],  // z = 0
        [4, 5, 6], [4, 6, 7],  // z = 1
        [0, 1, 5], [0, 5, 4],  // y = 0
        [2, 3, 7], [2, 7, 6],  // y = 1
        [0, 4, 7], [0, 7, 3],  // x = 0
        [1, 2, 6], [1, 6, 5]   // x = 1
    ];
    return mesh;
}

// The regular octahedron with vertices at +-1 along each axis, faces outward.
function octahedron(): ConvexMesh3 {
    const mesh = new ConvexMesh3();
    mesh.configuration = ConvexMesh3.CFG_POLYHEDRON;
    mesh.vertices = [
        v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
        v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)
    ];
    mesh.triangles = [
        [0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
        [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]
    ];
    return mesh;
}

// A copy of the mesh with its vertex list reindexed by 'perm', where the
// vertex at old index i moves to index perm[i].
function permuteVertices(mesh: ConvexMesh3, perm: readonly number[]):
    ConvexMesh3 {
    const copy = new ConvexMesh3();
    copy.configuration = mesh.configuration;
    copy.vertices = new Array<Vector>(perm.length);
    for (let i = 0; i < perm.length; ++i) {
        copy.vertices[perm[i]] = mesh.vertices[i].clone();
    }
    copy.triangles = mesh.triangles.map(t =>
        [perm[t[0]], perm[t[1]], perm[t[2]]] as ConvexMesh3Triangle);
    return copy;
}

// The volume enclosed by a closed mesh whose triangles face outward.
function volume(mesh: ConvexMesh3): number {
    let sum = 0;
    for (const t of mesh.triangles) {
        const a = mesh.vertices[t[0]];
        const b = mesh.vertices[t[1]];
        const c = mesh.vertices[t[2]];
        sum += dot(a, cross(b, c));
    }
    return sum / 6;
}

// The signed distance of X from the plane.
function signedDistance(plane: Hyperplane, X: Vector): number {
    return dot(plane.normal, X) - plane.constant;
}

// The area of a planar polygon, using the Newell normal.
function polygonArea(polygon: readonly Vector[]): number {
    let n = v3(0, 0, 0);
    for (let i0 = polygon.length - 1, i1 = 0; i1 < polygon.length; i0 = i1++) {
        n = Vector.fromArray([
            n.values[0] + cross(polygon[i0], polygon[i1]).values[0],
            n.values[1] + cross(polygon[i0], polygon[i1]).values[1],
            n.values[2] + cross(polygon[i0], polygon[i1]).values[2]
        ]);
    }
    return 0.5 * Math.sqrt(dot(n, n));
}

// Verify that the polygon is planar with respect to the plane, is traversed
// in one rotational direction and has no repeated vertices.
function expectSimpleConvexPolygon(polygon: readonly Vector[],
    plane: Hyperplane): void {
    expect(polygon.length).toBeGreaterThanOrEqual(3);
    for (const p of polygon) {
        expect(Math.abs(signedDistance(plane, p))).toBeLessThan(1e-12);
    }
    let sign = 0;
    for (let i = 0; i < polygon.length; ++i) {
        const a = polygon[i];
        const b = polygon[(i + 1) % polygon.length];
        const c = polygon[(i + 2) % polygon.length];
        expect(Math.sqrt(dot(sub(b, a), sub(b, a)))).toBeGreaterThan(1e-9);
        const turn = dot(plane.normal, cross(sub(b, a), sub(c, b)));
        if (Math.abs(turn) > 1e-12) {
            const s = turn > 0 ? 1 : -1;
            if (sign === 0) {
                sign = s;
            }
            else {
                expect(s).toBe(sign);
            }
        }
    }
    expect(sign).not.toBe(0);
}

// Every vertex of the mesh is on the given side of the plane (or on it).
function expectOnSide(mesh: ConvexMesh3, plane: Hyperplane,
    side: number): void {
    for (const vertex of mesh.vertices) {
        expect(side * signedDistance(plane, vertex)).toBeGreaterThan(-1e-12);
    }
}

describe('IntrConvexMesh3Plane3', () => {
    it('has the documented configuration and request constants', () => {
        expect(Q.CFG_EMPTY).toBe(0);
        expect(Q.CFG_SPLIT).toBe(48);
        expect(Q.CFG_POS_SIDE_STRICT).toBe(16);
        expect(Q.CFG_POS_SIDE_VERTEX).toBe(17);
        expect(Q.CFG_POS_SIDE_EDGE).toBe(18);
        expect(Q.CFG_POS_SIDE_POLYGON).toBe(20);
        expect(Q.CFG_POS_SIDE_TANGENT).toBe(23);
        expect(Q.CFG_NEG_SIDE_STRICT).toBe(32);
        expect(Q.CFG_NEG_SIDE_VERTEX).toBe(33);
        expect(Q.CFG_NEG_SIDE_EDGE).toBe(34);
        expect(Q.CFG_NEG_SIDE_POLYGON).toBe(36);
        expect(Q.CFG_NEG_SIDE_TANGENT).toBe(39);
        expect(Q.REQ_INTR_BOTH).toBe(3);
        expect(Q.REQ_POLYHEDRON_BOTH).toBe(12);
        expect(Q.REQ_ALL).toBe(15);

        const result = defaultIntrConvexMesh3Plane3FIResult();
        expect(result.configuration).toBe(Q.CFG_EMPTY);
        expect(result.requested).toBe(Q.REQ_CONFIGURATION_ONLY);
        expect(result.intersectionPolygon).toHaveLength(0);
    });

    it('classifies a polyhedron strictly on one side of the plane', () => {
        const query = new Q();
        const cube = unitCube();

        // The cube is strictly below z = 2.
        const above = Hyperplane.fromNormalConstant(v3(0, 0, 1), 2);
        const rAbove = query.find(cube, above, Q.REQ_ALL);
        expect(rAbove.configuration).toBe(Q.CFG_NEG_SIDE_STRICT);
        expect(rAbove.intersectionPolygon).toHaveLength(0);
        expect(rAbove.intersectionMesh.vertices).toHaveLength(0);
        expect(volume(rAbove.negativePolyhedron)).toBeCloseTo(1, 12);
        expect(rAbove.positivePolyhedron.vertices).toHaveLength(0);

        // The cube is strictly above z = -2.
        const below = Hyperplane.fromNormalConstant(v3(0, 0, 1), -2);
        const rBelow = query.find(cube, below, Q.REQ_ALL);
        expect(rBelow.configuration).toBe(Q.CFG_POS_SIDE_STRICT);
        expect(volume(rBelow.positivePolyhedron)).toBeCloseTo(1, 12);
        expect(rBelow.negativePolyhedron.vertices).toHaveLength(0);
    });

    it('reports a single vertex of tangency', () => {
        const query = new Q();
        const cube = unitCube();
        // The plane x+y+z = 0 supports the cube at the origin.
        const plane = Hyperplane.fromNormalConstant(v3(1, 1, 1), 0);
        const result = query.find(cube, plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_POS_SIDE_VERTEX);
        expect(result.intersectionMesh.configuration)
            .toBe(ConvexMesh3.CFG_POINT);
        expect(result.intersectionPolygon).toHaveLength(1);
        expect(result.intersectionPolygon[0].values).toEqual([0, 0, 0]);
        expect(result.intersectionMesh.vertices).toHaveLength(1);
        expect(result.intersectionMesh.triangles).toHaveLength(0);
    });

    it('reports a single edge of tangency', () => {
        const query = new Q();
        const cube = unitCube();
        // The plane x+y = 0 supports the cube along the edge x = y = 0.
        const plane = Hyperplane.fromNormalConstant(v3(1, 1, 0), 0);
        const result = query.find(cube, plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_POS_SIDE_EDGE);
        expect(result.intersectionMesh.configuration)
            .toBe(ConvexMesh3.CFG_SEGMENT);
        expect(result.intersectionPolygon).toHaveLength(2);
        const points = result.intersectionPolygon.map(p => p.values.join(','));
        expect(points.sort()).toEqual(['0,0,0', '0,0,1']);
        expect(result.intersectionMesh.triangles).toHaveLength(0);
    });

    it('reports a coplanar face as an ordered boundary polygon', () => {
        const query = new Q();
        const cube = unitCube();
        // The plane z = 1 contains the top face of the cube.
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 1);
        const result = query.find(cube, plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_NEG_SIDE_POLYGON);
        expect(result.intersectionMesh.configuration)
            .toBe(ConvexMesh3.CFG_POLYGON);
        expect(result.intersectionMesh.vertices).toHaveLength(4);
        expect(result.intersectionMesh.triangles).toHaveLength(2);

        // The polygon must be the boundary square in cyclic order, so each
        // consecutive pair is a unit-length side, not a diagonal. (With the
        // upstream index-map code the order is an arbitrary permutation.)
        const polygon = result.intersectionPolygon;
        expect(polygon).toHaveLength(4);
        expectSimpleConvexPolygon(polygon, plane);
        for (let i = 0; i < 4; ++i) {
            const e = sub(polygon[(i + 1) % 4], polygon[i]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(1, 12);
        }
        expect(polygonArea(polygon)).toBeCloseTo(1, 12);

        // The same face, with the cube vertices relabeled so that the
        // boundary cycle is not monotone in the vertex indices. This is the
        // regression test for the upstream ordering bug: the upstream code
        // reports the square as two sides and two diagonals here.
        const shuffled = permuteVertices(cube, [0, 1, 2, 3, 4, 6, 5, 7]);
        const rShuffled = query.find(shuffled, plane, Q.REQ_ALL);
        expect(rShuffled.intersectionPolygon).toHaveLength(4);
        expectSimpleConvexPolygon(rShuffled.intersectionPolygon, plane);
        for (let i = 0; i < 4; ++i) {
            const e = sub(rShuffled.intersectionPolygon[(i + 1) % 4],
                rShuffled.intersectionPolygon[i]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(1, 12);
        }
        expect(polygonArea(rShuffled.intersectionPolygon)).toBeCloseTo(1, 12);

        // The plane z = 0 contains the bottom face, with the cube above it.
        const bottom = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);
        const rBottom = query.find(cube, bottom, Q.REQ_ALL);
        expect(rBottom.configuration).toBe(Q.CFG_POS_SIDE_POLYGON);
        expect(rBottom.intersectionPolygon).toHaveLength(4);
        expectSimpleConvexPolygon(rBottom.intersectionPolygon, bottom);
        expect(polygonArea(rBottom.intersectionPolygon)).toBeCloseTo(1, 12);
    });

    it('splits the cube by a plane parallel to a face', () => {
        const query = new Q();
        const cube = unitCube();
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0.25);
        const result = query.find(cube, plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_SPLIT);

        expect(volume(result.positivePolyhedron)).toBeCloseTo(0.75, 10);
        expect(volume(result.negativePolyhedron)).toBeCloseTo(0.25, 10);
        expectOnSide(result.positivePolyhedron, plane, +1);
        expectOnSide(result.negativePolyhedron, plane, -1);

        // The polygon of intersection is the unit square at z = 0.25. Each
        // of the four cut side faces of the cube is split into two triangles,
        // so the polygon carries an extra collinear vertex where the plane
        // crosses each face diagonal: eight vertices in all. The
        // triangulation adds the average point as an interior vertex.
        expect(result.intersectionPolygon).toHaveLength(8);
        expectSimpleConvexPolygon(result.intersectionPolygon, plane);
        expect(polygonArea(result.intersectionPolygon)).toBeCloseTo(1, 10);
        expect(result.intersectionMesh.configuration)
            .toBe(ConvexMesh3.CFG_POLYGON);
        expect(result.intersectionMesh.vertices).toHaveLength(9);
        expect(result.intersectionMesh.triangles).toHaveLength(8);
        for (const vertex of result.intersectionMesh.vertices) {
            expect(Math.abs(signedDistance(plane, vertex)))
                .toBeLessThan(1e-12);
        }
    });

    it('splits the cube by a plane that cuts off a corner', () => {
        const query = new Q();
        const cube = unitCube();
        // x + y + z = 0.5 cuts off the corner tetrahedron at the origin.
        const plane = Hyperplane.fromNormalConstant(v3(1, 1, 1), 0.5);
        const result = query.find(cube, plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_SPLIT);

        // The corner piece is a tetrahedron with legs of length 0.5.
        expect(volume(result.negativePolyhedron)).toBeCloseTo(0.5 ** 3 / 6, 10);
        expect(volume(result.positivePolyhedron)).toBeCloseTo(
            1 - 0.5 ** 3 / 6, 10);
        expectOnSide(result.positivePolyhedron, plane, +1);
        expectOnSide(result.negativePolyhedron, plane, -1);

        // The three cut faces each contribute a collinear midpoint, so the
        // triangular cross section is reported with six vertices.
        expect(result.intersectionPolygon).toHaveLength(6);
        expectSimpleConvexPolygon(result.intersectionPolygon, plane);
        // An equilateral triangle with sides of length sqrt(2)/2.
        const side = Math.SQRT2 / 2;
        expect(polygonArea(result.intersectionPolygon)).toBeCloseTo(
            Math.sqrt(3) / 4 * side * side, 10);
    });

    it('honors the requested-information flags', () => {
        const query = new Q();
        const cube = unitCube();
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0.5);

        const only = query.find(cube, plane, Q.REQ_CONFIGURATION_ONLY);
        expect(only.configuration).toBe(Q.CFG_SPLIT);
        expect(only.positivePolyhedron.vertices).toHaveLength(0);
        expect(only.negativePolyhedron.vertices).toHaveLength(0);
        expect(only.intersectionPolygon).toHaveLength(0);
        expect(only.intersectionMesh.vertices).toHaveLength(0);

        const posOnly = query.find(cube, plane, Q.REQ_POLYHEDRON_POS);
        expect(posOnly.positivePolyhedron.vertices.length).toBeGreaterThan(0);
        expect(posOnly.negativePolyhedron.vertices).toHaveLength(0);
        expect(posOnly.intersectionPolygon).toHaveLength(0);

        const polygonOnly = query.find(cube, plane, Q.REQ_INTR_POLYGON);
        expect(polygonOnly.intersectionPolygon).toHaveLength(8);
        expect(polygonOnly.intersectionMesh.vertices).toHaveLength(0);
        expect(polygonOnly.positivePolyhedron.vertices).toHaveLength(0);

        // The default value of 'requested' is REQ_ALL.
        const all = query.find(cube, plane);
        expect(all.requested).toBe(Q.REQ_ALL);
        expect(all.positivePolyhedron.vertices.length).toBeGreaterThan(0);
        expect(all.negativePolyhedron.vertices.length).toBeGreaterThan(0);
    });

    it('does not alias or modify the input mesh', () => {
        const query = new Q();
        const cube = unitCube();
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 2);
        const result = query.find(cube, plane, Q.REQ_ALL);
        expect(result.negativePolyhedron.vertices[0]).not.toBe(
            cube.vertices[0]);
        result.negativePolyhedron.vertices[0].values[0] = 99;
        expect(cube.vertices[0].values[0]).toBe(0);
    });

    it('conserves volume and side classification for random planes', () => {
        // A deterministic linear congruential generator keeps the test
        // reproducible.
        let seed = 424242;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const query = new Q();
        const meshes: [ConvexMesh3, number][] = [
            [unitCube(), 1],
            [octahedron(), 4 / 3]
        ];
        let numSplit = 0;
        for (const [mesh, totalVolume] of meshes) {
            expect(volume(mesh)).toBeCloseTo(totalVolume, 12);
            for (let trial = 0; trial < 120; ++trial) {
                const normal = v3(2 * rand() - 1, 2 * rand() - 1,
                    2 * rand() - 1);
                const lengthN = Math.sqrt(dot(normal, normal));
                if (lengthN < 1e-3) {
                    continue;
                }
                for (let d = 0; d < 3; ++d) {
                    normal.values[d] /= lengthN;
                }
                const constant = 2 * rand() - 1;
                const plane = Hyperplane.fromNormalConstant(normal, constant);
                const result = query.find(mesh, plane, Q.REQ_ALL);

                if (result.configuration !== Q.CFG_SPLIT) {
                    continue;
                }
                ++numSplit;

                const volPos = volume(result.positivePolyhedron);
                const volNeg = volume(result.negativePolyhedron);
                expect(volPos).toBeGreaterThan(0);
                expect(volNeg).toBeGreaterThan(0);
                expect(volPos + volNeg).toBeCloseTo(totalVolume, 8);

                expectOnSide(result.positivePolyhedron, plane, +1);
                expectOnSide(result.negativePolyhedron, plane, -1);

                // The polygon of intersection is planar, simple and convex.
                expectSimpleConvexPolygon(result.intersectionPolygon, plane);

                // The intersection mesh triangulates the same polygon.
                expect(result.intersectionMesh.vertices.length).toBe(
                    result.intersectionPolygon.length + 1);
                expect(result.intersectionMesh.triangles.length).toBe(
                    result.intersectionPolygon.length);
                let meshArea = 0;
                for (const t of result.intersectionMesh.triangles) {
                    const a = result.intersectionMesh.vertices[t[0]];
                    const b = result.intersectionMesh.vertices[t[1]];
                    const c = result.intersectionMesh.vertices[t[2]];
                    const n = cross(sub(b, a), sub(c, a));
                    meshArea += 0.5 * Math.sqrt(dot(n, n));
                }
                expect(meshArea).toBeCloseTo(
                    polygonArea(result.intersectionPolygon), 8);
            }
        }
        expect(numSplit).toBeGreaterThan(100);
    });

    it('detects the impossible all-zero triangle case', () => {
        // A "polyhedron" that is a single triangle in the plane z = 0 with an
        // extra vertex off the plane on each side. The split branch is
        // entered (there are positive and negative vertices), and the
        // coplanar triangle triggers the upstream LogError.
        const mesh = new ConvexMesh3();
        mesh.configuration = ConvexMesh3.CFG_POLYHEDRON;
        mesh.vertices = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1), v3(0, 0, -1)
        ];
        mesh.triangles = [[0, 1, 2] as ConvexMesh3Triangle];
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);
        const query = new Q();
        expect(() => query.find(mesh, plane, Q.REQ_ALL)).toThrow();
    });
});
