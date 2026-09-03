import { describe, it, expect } from 'vitest';
import {
    IntrTriangle3Triangle3TI,
    IntrTriangle3Triangle3FI,
    defaultIntrTriangle3Triangle3FIResult,
    defaultIntrTriangle3Triangle3TIResult
} from '../src/IntrTriangle3Triangle3.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function tri(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(v0), Vector.fromArray(v1),
        Vector.fromArray(v2));
}

function translated(T: Triangle, d: Vector): Triangle {
    return Triangle.fromVertices(add(T.v[0], d), add(T.v[1], d),
        add(T.v[2], d));
}

// The signed distance from X to the plane of the triangle.
function planeDistance(T: Triangle, X: Vector): number {
    const n = cross(sub(T.v[1], T.v[0]), sub(T.v[2], T.v[0]));
    return dot(n, sub(X, T.v[0])) / Math.sqrt(dot(n, n));
}

// True when X is in the (closed) triangle, allowing for round-off.
function inTriangle(T: Triangle, X: Vector, tol: number): boolean {
    const n = cross(sub(T.v[1], T.v[0]), sub(T.v[2], T.v[0]));
    const area2 = dot(n, n);
    if (area2 === 0) {
        return false;
    }
    if (Math.abs(planeDistance(T, X)) > tol) {
        return false;
    }
    for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
        const e = sub(T.v[i1], T.v[i0]);
        const d = sub(X, T.v[i0]);
        if (dot(cross(e, d), n) < -tol * area2) {
            return false;
        }
    }
    return true;
}

// Sort the points lexicographically so that a result can be compared to an
// expected point set independently of the ordering the query produces.
function sorted(points: readonly number[][]): number[][] {
    return points.slice().sort(
        (a, b) => (a[0] - b[0]) || (a[1] - b[1]) || (a[2] - b[2]));
}

function toArrays(points: readonly Vector[]): number[][] {
    return points.map((p) => [p.values[0], p.values[1], p.values[2]]);
}

function expectPoints(actual: readonly Vector[], expected: number[][]): void {
    expect(actual.length).toBe(expected.length);
    const a = sorted(toArrays(actual));
    const b = sorted(expected);
    for (let i = 0; i < a.length; ++i) {
        for (let j = 0; j < 3; ++j) {
            expect(a[i][j]).toBeCloseTo(b[i][j], 10);
        }
    }
}

const ti = new IntrTriangle3Triangle3TI();
const fi = new IntrTriangle3Triangle3FI();

// The reference triangle in the plane z = 0.
const base = tri([0, 0, 0], [4, 0, 0], [0, 4, 0]);

describe('IntrTriangle3Triangle3', () => {
    it('default-constructs the results as no intersection', () => {
        expect(defaultIntrTriangle3Triangle3TIResult())
            .toEqual({ intersect: false, contactTime: 0 });
        const r = defaultIntrTriangle3Triangle3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.contactTime).toBe(0);
        expect(r.intersection).toEqual([]);
    });

    it('finds the segment for the (1,2,0) configuration', () => {
        // One vertex of triangle1 is below the plane of triangle0 and two are
        // above, so two edges are clipped.
        const T1 = tri([1, 1, -1], [2, 1, 1], [1, 2, 1]);
        const result = fi.find(base, T1);
        expect(result.intersect).toBe(true);
        expectPoints(result.intersection, [[1.5, 1, 0], [1, 1.5, 0]]);
        expect(ti.test(base, T1).intersect).toBe(true);
    });

    it('finds the segment for the (1,1,1) configuration', () => {
        // One vertex is on the plane of triangle0 and the opposite edge
        // crosses the plane.
        const T1 = tri([1, 1, 0], [1, 2, 1], [2, 1, -1]);
        const result = fi.find(base, T1);
        expect(result.intersect).toBe(true);
        expectPoints(result.intersection, [[1, 1, 0], [1.5, 1.5, 0]]);
        expect(ti.test(base, T1).intersect).toBe(true);
    });

    it('finds the single vertex for the (0,2,1) configuration', () => {
        const T1 = tri([1, 1, 0], [1, 1, 2], [2, 1, 2]);
        const result = fi.find(base, T1);
        expect(result.intersect).toBe(true);
        expectPoints(result.intersection, [[1, 1, 0]]);
        expect(ti.test(base, T1).intersect).toBe(true);
    });

    it('rejects a touching vertex that misses the triangle', () => {
        // The vertex on the plane of triangle0 is outside triangle0.
        const T1 = tri([9, 9, 0], [9, 9, 2], [10, 9, 2]);
        const result = fi.find(base, T1);
        expect(result.intersect).toBe(false);
        expect(result.intersection).toEqual([]);
        expect(ti.test(base, T1).intersect).toBe(false);
    });

    it('finds the coincident edge for the (0,1,2) configuration', () => {
        const T1 = tri([1, 1, 0], [3, 1, 0], [1, 1, 2]);
        const result = fi.find(base, T1);
        expect(result.intersect).toBe(true);
        expectPoints(result.intersection, [[1, 1, 0], [3, 1, 0]]);
        expect(ti.test(base, T1).intersect).toBe(true);
    });

    it('clips a coincident edge that extends beyond the triangle', () => {
        // The edge on the plane runs from inside triangle0 to well outside;
        // only the portion inside triangle0 is reported. Triangle0 is bounded
        // by x + y <= 4, so the clip is at x = 3.
        const T1 = tri([1, 1, 0], [9, 1, 0], [1, 1, 2]);
        const result = fi.find(base, T1);
        expect(result.intersect).toBe(true);
        expectPoints(result.intersection, [[1, 1, 0], [3, 1, 0]]);
    });

    it('reports no intersection when triangle1 is strictly on one side', () => {
        const above = tri([1, 1, 1], [2, 1, 2], [1, 2, 3]);
        const below = tri([1, 1, -1], [2, 1, -2], [1, 2, -3]);
        for (const T1 of [above, below]) {
            expect(fi.find(base, T1).intersect).toBe(false);
            expect(ti.test(base, T1).intersect).toBe(false);
        }
    });

    it('finds the polygon of intersection of coplanar triangles', () => {
        const T1 = tri([1, 1, 0], [5, 1, 0], [1, 5, 0]);
        const result = fi.find(base, T1);
        expect(result.intersect).toBe(true);
        expect(result.intersection.length).toBeGreaterThanOrEqual(3);
        for (const p of result.intersection) {
            expect(p.values[2]).toBeCloseTo(0, 10);
            expect(inTriangle(base, p, 1e-9)).toBe(true);
            expect(inTriangle(T1, p, 1e-9)).toBe(true);
        }
        expect(ti.test(base, T1).intersect).toBe(true);
    });

    it('reports coplanar triangles that do not overlap as disjoint', () => {
        const T1 = tri([9, 9, 0], [12, 9, 0], [9, 12, 0]);
        expect(fi.find(base, T1).intersect).toBe(false);
        expect(ti.test(base, T1).intersect).toBe(false);
    });

    it('reports parallel but non-coplanar triangles as disjoint', () => {
        const T1 = translated(base, vec(0, 0, 1));
        expect(fi.find(base, T1).intersect).toBe(false);
        expect(ti.test(base, T1).intersect).toBe(false);
    });

    it('is invariant under translation of both triangles', () => {
        const T1 = tri([1, 1, -1], [2, 1, 1], [1, 2, 1]);
        const d = vec(-7, 3.5, 11);
        const moved = fi.find(translated(base, d), translated(T1, d));
        const original = fi.find(base, T1);
        expect(moved.intersect).toBe(original.intersect);
        expectPoints(moved.intersection,
            toArrays(original.intersection.map((p) => add(p, d))));
    });

    it('is symmetric in the TI query', () => {
        const others = [
            tri([1, 1, -1], [2, 1, 1], [1, 2, 1]),
            tri([1, 1, 1], [2, 1, 2], [1, 2, 3]),
            tri([1, 1, 0], [5, 1, 0], [1, 5, 0]),
            tri([9, 9, 0], [12, 9, 0], [9, 12, 0]),
            translated(base, vec(0, 0, 1))
        ];
        for (const other of others) {
            expect(ti.test(base, other).intersect)
                .toBe(ti.test(other, base).intersect);
        }
    });

    it('reports contact time zero for already-intersecting moving triangles', () => {
        const T1 = tri([1, 1, -1], [2, 1, 1], [1, 2, 1]);
        const zero = vec(0, 0, 0);
        const result = ti.testDynamic(5, base, zero, T1, zero);
        expect(result.intersect).toBe(true);
        expect(result.contactTime).toBe(0);
    });

    it('computes the contact time of approaching parallel triangles', () => {
        const T0 = tri([-5, -5, 0], [5, -5, 0], [0, 5, 0]);
        const T1 = translated(T0, vec(0, 0, 3));
        const result = ti.testDynamic(10, T0, vec(0, 0, 0), T1, vec(0, 0, -1));
        expect(result.intersect).toBe(true);
        expect(result.contactTime).toBeCloseTo(3, 10);

        // The contact occurs after the maximum time.
        expect(ti.testDynamic(2, T0, vec(0, 0, 0), T1, vec(0, 0, -1)).intersect)
            .toBe(false);
        // Moving apart never contacts.
        expect(ti.testDynamic(100, T0, vec(0, 0, 0), T1, vec(0, 0, 1)).intersect)
            .toBe(false);
    });

    it('computes the contact set of approaching non-parallel triangles', () => {
        // A large vertical triangle sweeping down toward the horizontal one.
        const T0 = tri([-5, -5, 0], [5, -5, 0], [0, 5, 0]);
        const T1 = tri([0, 0, 4], [1, 0, 6], [-1, 0, 6]);
        const result = fi.findDynamic(10, T0, vec(0, 0, 0), T1, vec(0, 0, -1));
        expect(result.intersect).toBe(true);
        expect(result.contactTime).toBeCloseTo(4, 8);
        expect(result.intersection.length).toBeGreaterThan(0);
        for (const p of result.intersection) {
            expect(p.values[2]).toBeCloseTo(0, 8);
        }
    });

    it('reports no dynamic intersection for triangles that never meet', () => {
        const T0 = tri([-1, -1, 0], [1, -1, 0], [0, 1, 0]);
        const T1 = translated(T0, vec(20, 0, 3));
        expect(ti.testDynamic(10, T0, vec(0, 0, 0), T1, vec(0, 0, -1)).intersect)
            .toBe(false);
        expect(fi.findDynamic(10, T0, vec(0, 0, 0), T1, vec(0, 0, -1)).intersect)
            .toBe(false);
    });

    it('agrees between TI, FI and the dynamic queries (randomized)', () => {
        let seed = 60130212;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();
        const zero = vec(0, 0, 0);

        let numHits = 0;
        for (let trial = 0; trial < 1500; ++trial) {
            // Keep the triangles small and close so that intersections are
            // common but not degenerate.
            const T0 = tri(
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)],
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)],
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)]);
            const T1 = tri(
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)],
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)],
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)]);

            const tiResult = ti.test(T0, T1).intersect;
            const fiResult = fi.find(T0, T1);

            // Stationary triangles with zero velocities must reproduce the
            // stationary test.
            const dynamicResult = ti.testDynamic(0, T0, zero, T1, zero);
            expect(dynamicResult.intersect).toBe(tiResult);
            if (dynamicResult.intersect) {
                expect(dynamicResult.contactTime).toBe(0);
            }

            if (fiResult.intersect) {
                ++numHits;
                expect(fiResult.intersection.length).toBeGreaterThan(0);
                // The TI query must agree that they intersect.
                expect(tiResult).toBe(true);
                // Every reported point is on both triangles.
                for (const p of fiResult.intersection) {
                    expect(inTriangle(T0, p, 1e-7)).toBe(true);
                    expect(inTriangle(T1, p, 1e-7)).toBe(true);
                }
                // The centroid of the intersection set is on both triangles
                // (the set is convex).
                let centroid = Vector.zero(3);
                for (const p of fiResult.intersection) {
                    centroid = add(centroid, p);
                }
                centroid = mul(1 / fiResult.intersection.length, centroid);
                expect(inTriangle(T0, centroid, 1e-6)).toBe(true);
                expect(inTriangle(T1, centroid, 1e-6)).toBe(true);
            }
        }
        // Random triangles are never exactly coplanar, so the intersection
        // set here is always a point or a segment; the polygon case is
        // covered by the coplanar tests above.
        expect(numHits).toBeGreaterThan(50);
    });

    it('finds the same contact set as the stationary query at contact time (randomized)', () => {
        let seed = 8081991;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        let numContacts = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const T0 = tri(
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)],
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)],
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)]);
            const T1 = translated(tri(
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)],
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)],
                [rnd(-2, 2), rnd(-2, 2), rnd(-2, 2)]), vec(0, 0, 8));
            const velocity1 = vec(0, 0, -1);
            const result = fi.findDynamic(100, T0, vec(0, 0, 0), T1, velocity1);
            if (!result.intersect) {
                continue;
            }
            ++numContacts;
            expect(result.contactTime).toBeGreaterThanOrEqual(0);
            // Move triangle1 to the contact time and rerun the stationary
            // query; the results must be identical.
            const movedT1 = translated(T1, mul(result.contactTime, velocity1));
            const stationary = fi.find(T0, movedT1);
            expectPoints(result.intersection, toArrays(stationary.intersection));
        }
        expect(numContacts).toBeGreaterThan(50);
    });
});
