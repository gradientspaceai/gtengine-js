import { describe, it, expect } from 'vitest';
import { ConvexHull3 } from '../src/ConvexHull3';
import { Vector } from '../src/Vector';

const v3 = (x: number, y: number, z: number): Vector =>
    Vector.fromArray([x, y, z]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function sub(a: Vector, b: Vector): [number, number, number] {
    return [a.values[0] - b.values[0], a.values[1] - b.values[1],
        a.values[2] - b.values[2]];
}

function cross(a: readonly number[], b: readonly number[]): [number, number, number] {
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]];
}

function dot(a: readonly number[], b: readonly number[]): number {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Every input point must be on the negative side of (or on) every face plane,
// where the face normal is the outward normal implied by the counterclockwise
// face ordering.
function checkConvexAndContains(points: readonly Vector[],
    hull: readonly number[], tolerance: number): void {
    expect(hull.length % 3).toBe(0);
    expect(hull.length).toBeGreaterThan(0);
    for (let f = 0; f < hull.length; f += 3) {
        const v0 = points[hull[f]];
        const v1 = points[hull[f + 1]];
        const v2 = points[hull[f + 2]];
        const normal = cross(sub(v1, v0), sub(v2, v0));
        const normalLength = Math.sqrt(dot(normal, normal));
        expect(normalLength).toBeGreaterThan(0);
        for (const p of points) {
            const signedDistance = dot(normal, sub(p, v0)) / normalLength;
            expect(signedDistance).toBeLessThanOrEqual(tolerance);
        }
    }
}

// The hull mesh is closed and oriented, so counting the distinct undirected
// edges must give E = 3*T/2 and Euler's formula V - E + T = 2 must hold.
function checkEuler(numVertices: number, hull: readonly number[]): void {
    const numTriangles = hull.length / 3;
    const edges = new Set<string>();
    for (let f = 0; f < hull.length; f += 3) {
        for (let i = 0; i < 3; ++i) {
            const a = hull[f + i];
            const b = hull[f + (i + 1) % 3];
            edges.add(a < b ? `${a},${b}` : `${b},${a}`);
        }
    }
    expect(edges.size).toBe(3 * numTriangles / 2);
    expect(numVertices - edges.size + numTriangles).toBe(2);
}

describe('ConvexHull3', () => {
    it('rejects invalid input', () => {
        const ch = new ConvexHull3();
        expect(() => ch.compute([])).toThrow();
        expect(() => ch.compute([Vector.fromArray([0, 0])])).toThrow();
    });

    it('computes a 0-dimensional hull for coincident points', () => {
        const points = [v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3)];
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(0);
        expect(ch.getHull().length).toBe(1);
        expect(points[ch.getHull()[0]].equals(v3(1, 2, 3))).toBe(true);
        expect(ch.getVertices().length).toBe(1);
    });

    it('computes a 1-dimensional hull for collinear points', () => {
        // Points on the line (t, 2t, 3t), given out of order and with a
        // duplicate.
        const ts = [3, -1, 0, 5, 2, 5];
        const points = ts.map(t => v3(t, 2 * t, 3 * t));
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(1);
        const hull = ch.getHull();
        expect(hull.length).toBe(2);
        const endpointTs = [ts[hull[0]], ts[hull[1]]].sort((a, b) => a - b);
        expect(endpointTs).toEqual([-1, 5]);
        expect(ch.getVertices().length).toBe(2);
    });

    it('computes a 2-dimensional hull for coplanar points', () => {
        // A unit square in the plane z = 4, plus interior and edge-interior
        // points that must not be hull vertices.
        const points = [
            v3(0, 0, 4), v3(1, 0, 4), v3(1, 1, 4), v3(0, 1, 4),
            v3(0.5, 0.5, 4), v3(0.5, 0, 4), v3(0.25, 0.75, 4)
        ];
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(2);
        const hull = ch.getHull();
        expect(hull.length).toBe(4);
        expect([...hull].sort((a, b) => a - b)).toEqual([0, 1, 2, 3]);

        // The polygon must be ordered (consecutive corners of the square).
        for (let i = 0; i < 4; ++i) {
            const a = points[hull[i]];
            const b = points[hull[(i + 1) % 4]];
            const d = sub(b, a);
            expect(Math.abs(dot(d, d) - 1)).toBeLessThan(1e-12);
        }
        expect(ch.getVertices().length).toBe(4);
    });

    it('computes the hull of the corners of a cube', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            points.push(v3(i & 1, (i >> 1) & 1, (i >> 2) & 1));
        }
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        expect(ch.getVertices().length).toBe(8);

        const hull = ch.getHull();
        // Each square face of the cube is split into 2 triangles.
        expect(hull.length).toBe(3 * 12);
        checkConvexAndContains(points, hull, 1e-12);
        checkEuler(8, hull);

        // The hull mesh is a closed, oriented manifold.
        const mesh = ch.getHullMesh();
        expect(mesh.getNumTriangles()).toBe(12);
        expect(mesh.getNumEdges()).toBe(18);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(true);
        expect(mesh.getNumVertices()).toBe(8);
    });

    it('ignores interior points and duplicates', () => {
        // Cube corners, the cube center (interior), a face center (on the
        // boundary but not a vertex) and a repeated corner.
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            points.push(v3(i & 1, (i >> 1) & 1, (i >> 2) & 1));
        }
        points.push(v3(0.5, 0.5, 0.5));
        points.push(v3(0.5, 0.5, 0));
        points.push(v3(1, 1, 1));
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        const vertices = [...ch.getVertices()];
        // All 8 distinct corners are hull vertices.
        for (let i = 0; i < 8; ++i) {
            expect(vertices).toContain(i);
        }
        // The cube center (index 8) is strictly interior, so it is not a hull
        // vertex. The repeated corner (index 10) is dropped by the uniqueness
        // filter in favor of the first occurrence (index 7).
        expect(vertices).not.toContain(8);
        expect(vertices).not.toContain(10);
        // The face center (index 9) lies on the boundary; upstream keeps such
        // a point as a vertex of coplanar hull triangles, which the header
        // documents as expected behavior.
        expect(vertices.length).toBe(9);
        checkConvexAndContains(points, ch.getHull(), 1e-12);
        checkEuler(vertices.length, ch.getHull());
    });

    it('computes the hull of a regular octahedron', () => {
        const points = [
            v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
            v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)
        ];
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        expect(ch.getVertices().length).toBe(6);
        const hull = ch.getHull();
        expect(hull.length).toBe(3 * 8);
        checkConvexAndContains(points, hull, 1e-12);
        checkEuler(6, hull);
    });

    it('computes hulls of random clouds', () => {
        const random = makeRandom(20260902);
        for (let trial = 0; trial < 6; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 30; ++i) {
                points.push(v3(2 * random() - 1, 2 * random() - 1, 2 * random() - 1));
            }
            const ch = new ConvexHull3();
            ch.compute(points);
            expect(ch.getDimension()).toBe(3);

            const hull = ch.getHull();
            checkConvexAndContains(points, hull, 1e-9);
            checkEuler(ch.getVertices().length, hull);

            // The hull vertices are exactly the indices appearing in the
            // triangle list.
            const used = new Set<number>(hull);
            expect([...used].sort((a, b) => a - b))
                .toEqual([...ch.getVertices()].sort((a, b) => a - b));

            const mesh = ch.getHullMesh();
            expect(mesh.isClosed()).toBe(true);
            expect(mesh.isOriented()).toBe(true);
        }
    });

    it('computes hulls of random points on a sphere (all points are vertices)', () => {
        const random = makeRandom(777);
        const points: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            const z = 2 * random() - 1;
            const theta = 2 * Math.PI * random();
            const r = Math.sqrt(1 - z * z);
            points.push(v3(r * Math.cos(theta), r * Math.sin(theta), z));
        }
        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        // Points on a sphere are in convex position, so all are hull vertices.
        expect(ch.getVertices().length).toBe(points.length);
        checkConvexAndContains(points, ch.getHull(), 1e-9);
        checkEuler(points.length, ch.getHull());
    });

    it('is exact on a near-degenerate set (interval fallback to rational)', () => {
        // A grid of points on the plane z = 0 with one point lifted by a
        // quantity so tiny that the interval predicate cannot resolve the
        // sign of the determinant. The exact fallback must classify it as
        // above the plane, making the hull 3-dimensional.
        const points: Vector[] = [];
        for (let x = 0; x <= 2; ++x) {
            for (let y = 0; y <= 2; ++y) {
                points.push(v3(x, y, 0));
            }
        }
        const lift = Math.pow(2, -60);
        points.push(v3(1, 1, lift));

        const ch = new ConvexHull3();
        ch.compute(points);
        expect(ch.getDimension()).toBe(3);
        // The lifted point is a hull vertex even though it is only 2^-60
        // above the plane of the others.
        expect([...ch.getVertices()]).toContain(points.length - 1);
        checkConvexAndContains(points, ch.getHull(), 1e-12);
        checkEuler(ch.getVertices().length, ch.getHull());

        // Without the lift the same set is coplanar and the hull is 2D.
        const flat = points.slice(0, points.length - 1);
        const chFlat = new ConvexHull3();
        chFlat.compute(flat);
        expect(chFlat.getDimension()).toBe(2);
        expect(chFlat.getHull().length).toBe(4);
    });

    it('reuses the functor across data sets', () => {
        const ch = new ConvexHull3();
        ch.compute([v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)]);
        expect(ch.getDimension()).toBe(3);
        expect(ch.getHull().length).toBe(3 * 4);

        ch.compute([v3(5, 5, 5), v3(5, 5, 5)]);
        expect(ch.getDimension()).toBe(0);
        expect(ch.getHull().length).toBe(1);
        expect(ch.getVertices().length).toBe(1);
        expect(ch.getHullMesh().getNumTriangles()).toBe(0);
    });
});
