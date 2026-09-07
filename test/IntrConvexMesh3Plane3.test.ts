import { describe, it, expect } from 'vitest';
import { ConvexMesh3 } from '../src/ConvexMesh3.js';
import type { ConvexMesh3Triangle } from '../src/ConvexMesh3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, fc, rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';
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

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrConvexMesh3Plane3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrConvexMesh3Plane3 verification', () => {
    // Is the mesh a closed manifold? Every directed edge must appear exactly
    // once, and its reverse must appear exactly once.
    function expectClosedManifold(mesh: ConvexMesh3): void {
        const seen = new Map<string, number>();
        for (const t of mesh.triangles) {
            for (let j0 = 2, j1 = 0; j1 < 3; j0 = j1++) {
                const key = t[j0] + '->' + t[j1];
                seen.set(key, (seen.get(key) ?? 0) + 1);
            }
        }
        for (const [key, count] of seen) {
            expect(count).toBe(1);
            const [a, b] = key.split('->');
            expect(seen.get(b + '->' + a)).toBe(1);
        }
        // Euler characteristic of a sphere: V - E + F = 2.
        const used = new Set<number>();
        for (const t of mesh.triangles) {
            used.add(t[0]);
            used.add(t[1]);
            used.add(t[2]);
        }
        const numE = seen.size / 2;
        expect(used.size - numE + mesh.triangles.length).toBe(2);
    }

    // The convexity/simplicity assertions need a nondegenerate cut: a plane
    // that grazes a cube corner produces a polygon with a vanishing edge.
    function minEdgeLength(polygon: readonly Vector[]): number {
        let m = Number.POSITIVE_INFINITY;
        for (let i = 0; i < polygon.length; ++i) {
            const a = polygon[i];
            const b = polygon[(i + 1) % polygon.length];
            m = Math.min(m, Math.sqrt(dot(sub(b, a), sub(b, a))));
        }
        return m;
    }

    const arbPlaneThroughCube = fc.tuple(
        unitVector(3),
        fc.double({ min: -0.6, max: 0.6, noNaN: true }),
        fc.double({ min: -0.6, max: 0.6, noNaN: true }),
        fc.double({ min: -0.6, max: 0.6, noNaN: true }))
        .map(([n, x, y, z]) => Hyperplane.fromNormalOrigin(n,
            v3(0.5 + x, 0.5 + y, 0.5 + z)));

    it('splits the cube into two closed manifolds that conserve volume', () => {
        const query = new Q();
        check(arbPlaneThroughCube, plane => {
            const cube = unitCube();
            const result = query.find(cube, plane, Q.REQ_ALL);
            if (result.configuration !== Q.CFG_SPLIT) { return; }
            // Upstream requires exact arithmetic (its static_assert); with
            // binary64 a plane that passes within rounding of a mesh vertex
            // splits an edge at a distance of ~1e-16 from it, and the split
            // meshes then carry degenerate triangles. Those cuts are excluded
            // from the manifold assertions and pinned separately.
            if (minEdgeLength(result.intersectionPolygon) < 1e-6) { return; }
            expectClosedManifold(result.positivePolyhedron);
            expectClosedManifold(result.negativePolyhedron);
            expectClose(volume(result.positivePolyhedron)
                + volume(result.negativePolyhedron), 1, 1e-9, 1e-9);
            expect(volume(result.positivePolyhedron)).toBeGreaterThan(-1e-12);
            expect(volume(result.negativePolyhedron)).toBeGreaterThan(-1e-12);
            expectOnSide(result.positivePolyhedron, plane, +1);
            expectOnSide(result.negativePolyhedron, plane, -1);
        }, 100);
    });

    it('the cut polygon is planar, convex and closes the two halves', () => {
        const query = new Q();
        check(arbPlaneThroughCube, plane => {
            const cube = unitCube();
            const result = query.find(cube, plane, Q.REQ_ALL);
            if (result.configuration !== Q.CFG_SPLIT) { return; }
            if (minEdgeLength(result.intersectionPolygon) < 1e-6) { return; }
            expectSimpleConvexPolygon(result.intersectionPolygon, plane);
            // Every polygon vertex is a vertex of both split polyhedra.
            for (const p of result.intersectionPolygon) {
                const inPos = result.positivePolyhedron.vertices.some(
                    v => Math.abs(v.values[0] - p.values[0]) < 1e-12
                        && Math.abs(v.values[1] - p.values[1]) < 1e-12
                        && Math.abs(v.values[2] - p.values[2]) < 1e-12);
                const inNeg = result.negativePolyhedron.vertices.some(
                    v => Math.abs(v.values[0] - p.values[0]) < 1e-12
                        && Math.abs(v.values[1] - p.values[1]) < 1e-12
                        && Math.abs(v.values[2] - p.values[2]) < 1e-12);
                expect(inPos).toBe(true);
                expect(inNeg).toBe(true);
            }
            // The triangulated intersection mesh has the same area as the
            // polygon and lies in the plane.
            let meshArea = 0;
            for (const t of result.intersectionMesh.triangles) {
                const a = result.intersectionMesh.vertices[t[0]];
                const b = result.intersectionMesh.vertices[t[1]];
                const c = result.intersectionMesh.vertices[t[2]];
                meshArea += 0.5 * Math.sqrt(dot(cross(sub(b, a), sub(c, a)),
                    cross(sub(b, a), sub(c, a))));
            }
            expectClose(meshArea, polygonArea(result.intersectionPolygon),
                1e-9, 1e-9);
            for (const v of result.intersectionMesh.vertices) {
                expect(Math.abs(signedDistance(plane, v))).toBeLessThan(1e-11);
            }
        }, 100);
    });

    it('splits the octahedron the same way', () => {
        const query = new Q();
        check(fc.tuple(unitVector(3),
            fc.double({ min: -0.5, max: 0.5, noNaN: true })),
            ([n, c]) => {
                const plane = Hyperplane.fromNormalConstant(n, c);
                const mesh = octahedron();
                const result = query.find(mesh, plane, Q.REQ_ALL);
                if (result.configuration !== Q.CFG_SPLIT) { return; }
                expectClosedManifold(result.positivePolyhedron);
                expectClosedManifold(result.negativePolyhedron);
                expectClose(volume(result.positivePolyhedron)
                    + volume(result.negativePolyhedron), volume(mesh),
                    1e-9, 1e-9);
                if (minEdgeLength(result.intersectionPolygon) < 1e-6) {
                    return;
                }
                expectSimpleConvexPolygon(result.intersectionPolygon, plane);
            }, 100);
    });

    it('the split results do not alias the split-vertex pool or each other', () => {
        // UniqueVerticesSimplices packs references to the vertices it is
        // given, so the three output meshes would otherwise share Vector
        // objects; upstream copies vertices by value.
        const query = new Q();
        const cube = unitCube();
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0.4);
        const result = query.find(cube, plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_SPLIT);

        const posSet = new Set(result.positivePolyhedron.vertices);
        for (const v of result.negativePolyhedron.vertices) {
            expect(posSet.has(v)).toBe(false);
        }
        for (const v of result.intersectionMesh.vertices) {
            expect(posSet.has(v)).toBe(false);
        }
        for (const v of result.intersectionPolygon) {
            expect(posSet.has(v)).toBe(false);
        }
        // Mutating one mesh leaves the others (and the input) untouched.
        const negBefore = result.negativePolyhedron.vertices.map(
            v => v.values.slice());
        const cubeBefore = cube.vertices.map(v => v.values.slice());
        for (const v of result.positivePolyhedron.vertices) {
            v.values[0] = 42;
        }
        expect(result.negativePolyhedron.vertices.map(v => v.values.slice()))
            .toEqual(negBefore);
        expect(cube.vertices.map(v => v.values.slice())).toEqual(cubeBefore);
    });

    it('the coplanar-face polygon does not alias the input mesh', () => {
        const query = new Q();
        const cube = unitCube();
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 1);
        const result = query.find(cube, plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_NEG_SIDE_POLYGON);
        const cubeSet = new Set(cube.vertices);
        for (const v of result.intersectionMesh.vertices) {
            expect(cubeSet.has(v)).toBe(false);
        }
        for (const v of result.intersectionPolygon) {
            expect(cubeSet.has(v)).toBe(false);
        }
        const before = cube.vertices.map(v => v.values.slice());
        for (const v of result.intersectionMesh.vertices) {
            v.values[2] = -7;
        }
        expect(cube.vertices.map(v => v.values.slice())).toEqual(before);
    });

    it('places the triangulation apex at the reciprocal-scaled average', () => {
        // Upstream divides the accumulated polygon-vertex sum with
        // Vector::operator/=, which multiplies by the reciprocal of the
        // vertex count; a componentwise division differs in the last bit for
        // most sums. The cut below is a triangle whose x-coordinates sum to
        // 5, and 5/3 !== 5*(1/3) in binary64.
        expect(5 / 3).not.toBe(5 * (1 / 3));

        const mesh = new ConvexMesh3();
        mesh.configuration = ConvexMesh3.CFG_POLYHEDRON;
        // A tetrahedron whose z = 0 cross-section is the triangle with
        // vertices (1,0,0), (2.5,0,0) and (1.5,1.5,0): the x-coordinates sum
        // to exactly 5.
        mesh.vertices = [
            v3(1, 0, -1), v3(4, 0, -1), v3(2, 3, -1), v3(1, 0, 1)
        ];
        mesh.triangles = [
            [0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]
        ];
        const plane = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);
        const result = new Q().find(mesh, plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_SPLIT);

        const polygon = result.intersectionPolygon;
        expect(polygon.length).toBe(3);
        let sumX = 0;
        let sumY = 0;
        for (const p of polygon) {
            sumX += p.values[0];
            sumY += p.values[1];
        }
        expect(sumX).toBe(5);
        const apexX = sumX * (1 / 3);
        const apexY = sumY * (1 / 3);
        // The apex is the extra vertex of the positive polyhedron that is on
        // the plane but not on the polygon.
        const apex = result.positivePolyhedron.vertices.filter(v =>
            v.values[2] === 0 && !polygon.some(p =>
                p.values[0] === v.values[0] && p.values[1] === v.values[1]));
        expect(apex.length).toBe(1);
        expect(apex[0].values[0]).toBe(apexX);
        expect(apex[0].values[1]).toBe(apexY);
        // The componentwise division would give a different last bit.
        expect(apex[0].values[0]).not.toBe(sumX / 3);
    });

    it('the split is invariant under a relabeling of the input vertices', () => {
        const query = new Q();
        check(fc.tuple(arbPlaneThroughCube,
            fc.constantFrom(
                [0, 1, 2, 3, 4, 5, 6, 7],
                [7, 6, 5, 4, 3, 2, 1, 0],
                [2, 5, 0, 7, 4, 1, 6, 3],
                [3, 0, 6, 1, 7, 2, 5, 4])),
            ([plane, perm]) => {
                const a = query.find(unitCube(), plane, Q.REQ_ALL);
                const b = query.find(permuteVertices(unitCube(), perm), plane,
                    Q.REQ_ALL);
                expect(b.configuration).toBe(a.configuration);
                if (a.configuration !== Q.CFG_SPLIT) { return; }
                expectClose(volume(b.positivePolyhedron),
                    volume(a.positivePolyhedron), 1e-9, 1e-9);
                expectClose(volume(b.negativePolyhedron),
                    volume(a.negativePolyhedron), 1e-9, 1e-9);
                expect(b.intersectionPolygon.length)
                    .toBe(a.intersectionPolygon.length);
                expectClose(polygonArea(b.intersectionPolygon),
                    polygonArea(a.intersectionPolygon), 1e-9, 1e-9);
            }, 60);
    });

    it('is equivariant under a rigid motion of mesh and plane', () => {
        const query = new Q();
        check(fc.tuple(arbPlaneThroughCube, rotationFrame(3),
            wellScaledVector(3, -2, 2)),
            ([plane, R, t]) => {
                const rot = (x: Vector) => Vector.fromArray([
                    dot(R[0], x), dot(R[1], x), dot(R[2], x)]);
                const map = (x: Vector) => Vector.fromArray([
                    dot(R[0], x) + t.values[0], dot(R[1], x) + t.values[1],
                    dot(R[2], x) + t.values[2]]);
                const a = query.find(unitCube(), plane, Q.REQ_ALL);
                const moved = new ConvexMesh3();
                moved.configuration = ConvexMesh3.CFG_POLYHEDRON;
                moved.vertices = unitCube().vertices.map(map);
                moved.triangles = unitCube().triangles;
                const plane2 = Hyperplane.fromNormalOrigin(rot(plane.normal),
                    map(plane.origin));
                const b = query.find(moved, plane2, Q.REQ_ALL);
                if (a.configuration !== Q.CFG_SPLIT
                    || b.configuration !== Q.CFG_SPLIT) {
                    return;
                }
                expectClose(volume(b.positivePolyhedron),
                    volume(a.positivePolyhedron), 1e-9, 1e-9);
                expectClose(polygonArea(b.intersectionPolygon),
                    polygonArea(a.intersectionPolygon), 1e-9, 1e-9);
            }, 60);
    });

    it('reports the tangential configurations with the documented flags', () => {
        const query = new Q();
        const cube = unitCube();
        // A face in the plane.
        const face = query.find(cube,
            Hyperplane.fromNormalConstant(v3(0, 0, 1), 0), Q.REQ_ALL);
        expect(face.configuration).toBe(Q.CFG_POS_SIDE_POLYGON);
        expect(face.intersectionPolygon.length).toBe(4);
        expectSimpleConvexPolygon(face.intersectionPolygon,
            Hyperplane.fromNormalConstant(v3(0, 0, 1), 0));
        // The polygon has unit side lengths (the #301 cycle-order fix).
        for (let i = 0; i < 4; ++i) {
            const a = face.intersectionPolygon[i];
            const b = face.intersectionPolygon[(i + 1) % 4];
            expectClose(Math.sqrt(dot(sub(b, a), sub(b, a))), 1, 1e-12, 1e-12);
        }

        // A single edge in the plane.
        const n = v3(1, 0, 1);
        const edge = query.find(cube,
            Hyperplane.fromNormalConstant(n, 0), Q.REQ_ALL);
        expect(edge.configuration).toBe(Q.CFG_POS_SIDE_EDGE);
        expect(edge.intersectionPolygon.length).toBe(2);

        // A single vertex in the plane.
        const vertexPlane = Hyperplane.fromNormalConstant(v3(1, 1, 1), 0);
        const vertex = query.find(cube, vertexPlane, Q.REQ_ALL);
        expect(vertex.configuration).toBe(Q.CFG_POS_SIDE_VERTEX);
        expect(vertex.intersectionPolygon.length).toBe(1);
        expect(vertex.intersectionPolygon[0].values).toEqual([0, 0, 0]);
        // The reported vertex is a copy.
        expect(vertex.intersectionPolygon[0]).not.toBe(cube.vertices[0]);
    });

    it('pins the floating-point degeneracy of a near-vertex cut', () => {
        // Upstream static_asserts that Real is an arbitrary-precision type,
        // because the classification Dot(N,X) - c > 0 / < 0 / == 0 must be
        // exact. The port instantiates the query for 'number' (see the port
        // notes in the source), so a plane that passes within rounding of a
        // mesh vertex classifies that vertex strictly on one side and splits
        // the incident edges a few 1e-17 away from it. The volumes stay
        // right, but the split meshes carry degenerate triangles (a repeated
        // index) and the polygon repeats a vertex. Pinned so the behaviour is
        // visible rather than surprising.
        const n = v3(-0.7071067811865476, 0, 0.7071067811865475);
        const plane = Hyperplane.fromNormalConstant(n,
            -5.551115123125783e-17);
        // The cube diagonal plane z = x passes through four cube vertices.
        for (const i of [0, 3, 5, 6]) {
            expect(Math.abs(signedDistance(plane, unitCube().vertices[i])))
                .toBeLessThan(1e-16);
        }
        const result = new Q().find(unitCube(), plane, Q.REQ_ALL);
        expect(result.configuration).toBe(Q.CFG_SPLIT);
        // The volumes are still exactly right.
        expectClose(volume(result.positivePolyhedron), 0.5, 1e-12, 1e-12);
        expectClose(volume(result.negativePolyhedron), 0.5, 1e-12, 1e-12);
        // ... but the triangulation is degenerate.
        const degenerate = result.positivePolyhedron.triangles.some(t =>
            t[0] === t[1] || t[1] === t[2] || t[2] === t[0]);
        expect(degenerate).toBe(true);
        expect(minEdgeLength(result.intersectionPolygon)).toBe(0);
    });

    it('honors REQ_CONFIGURATION_ONLY without computing any geometry', () => {
        const query = new Q();
        check(arbPlaneThroughCube, plane => {
            const full = query.find(unitCube(), plane, Q.REQ_ALL);
            const only = query.find(unitCube(), plane,
                Q.REQ_CONFIGURATION_ONLY);
            expect(only.configuration).toBe(full.configuration);
            expect(only.intersectionPolygon.length).toBe(0);
            expect(only.intersectionMesh.vertices.length).toBe(0);
            expect(only.positivePolyhedron.vertices.length).toBe(0);
            expect(only.negativePolyhedron.vertices.length).toBe(0);
        }, 60);
    });
});
