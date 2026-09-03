import { describe, it, expect } from 'vitest';
import { Delaunay2, Delaunay2SearchInfo } from '../src/Delaunay2.js';
import { Vector } from '../src/Vector.js';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function cross(a: Vector, b: Vector, c: Vector): number {
    return (b.values[0] - a.values[0]) * (c.values[1] - a.values[1])
        - (b.values[1] - a.values[1]) * (c.values[0] - a.values[0]);
}

// The InCircle determinant. It is positive when d is strictly inside the
// circumcircle of the counterclockwise triangle (a,b,c).
function inCircle(a: Vector, b: Vector, c: Vector, d: Vector): number {
    const ax = a.values[0] - d.values[0], ay = a.values[1] - d.values[1];
    const bx = b.values[0] - d.values[0], by = b.values[1] - d.values[1];
    const cx = c.values[0] - d.values[0], cy = c.values[1] - d.values[1];
    const a2 = ax * ax + ay * ay;
    const b2 = bx * bx + by * by;
    const c2 = cx * cx + cy * cy;
    return a2 * (bx * cy - by * cx)
        - b2 * (ax * cy - ay * cx)
        + c2 * (ax * by - ay * bx);
}

// Structural checks applied to every successful triangulation.
function checkTriangulation(del: Delaunay2, points: Vector[]): void {
    const indices = del.getIndices();
    const adjacencies = del.getAdjacencies();
    const numTriangles = del.getNumTriangles();
    expect(indices.length).toBe(3 * numTriangles);
    expect(adjacencies.length).toBe(3 * numTriangles);
    expect(numTriangles).toBeGreaterThan(0);

    // Every triangle is counterclockwise and nondegenerate.
    for (let t = 0; t < numTriangles; ++t) {
        const [i0, i1, i2] = [indices[3 * t], indices[3 * t + 1], indices[3 * t + 2]];
        expect(cross(points[i0], points[i1], points[i2])).toBeGreaterThan(0);
    }

    // The adjacency relation is symmetric and consistent with shared edges.
    for (let t = 0; t < numTriangles; ++t) {
        for (let j = 0; j < 3; ++j) {
            const adj = adjacencies[3 * t + j];
            if (adj === -1) {
                continue;
            }
            expect(adj).toBeGreaterThanOrEqual(0);
            expect(adj).toBeLessThan(numTriangles);
            const v0 = indices[3 * t + j];
            const v1 = indices[3 * t + ((j + 1) % 3)];
            // The adjacent triangle contains the reversed directed edge.
            let found = false;
            for (let k = 0; k < 3; ++k) {
                if (indices[3 * adj + k] === v1
                    && indices[3 * adj + ((k + 1) % 3)] === v0) {
                    found = true;
                    expect(adjacencies[3 * adj + k]).toBe(t);
                }
            }
            expect(found).toBe(true);
        }
    }

    // Euler-style sanity: each directed edge appears exactly once.
    const directed = new Set<string>();
    for (let t = 0; t < numTriangles; ++t) {
        for (let j = 0; j < 3; ++j) {
            const key = `${indices[3 * t + j]},${indices[3 * t + ((j + 1) % 3)]}`;
            expect(directed.has(key)).toBe(false);
            directed.add(key);
        }
    }
}

// The empty-circumcircle (Delaunay) property, verified brute force with a
// tolerance to allow for cocircular configurations.
function checkEmptyCircumcircle(del: Delaunay2, points: Vector[],
    tolerance = 1e-9): void {
    const indices = del.getIndices();
    const numTriangles = del.getNumTriangles();
    const duplicates = del.getDuplicates();
    for (let t = 0; t < numTriangles; ++t) {
        const a = points[indices[3 * t]];
        const b = points[indices[3 * t + 1]];
        const c = points[indices[3 * t + 2]];
        for (let i = 0; i < points.length; ++i) {
            if (duplicates[i] !== i) {
                continue;
            }
            if (i === indices[3 * t] || i === indices[3 * t + 1]
                || i === indices[3 * t + 2]) {
                continue;
            }
            expect(inCircle(a, b, c, points[i])).toBeLessThanOrEqual(tolerance);
        }
    }
}

// The convex hull of the triangulation, from getHull(), as an ordered
// counterclockwise cycle of vertex indices.
function hullCycle(del: Delaunay2): number[] {
    const hull = del.getHull();
    expect(hull.length % 2).toBe(0);
    const next = new Map<number, number>();
    for (let i = 0; i < hull.length; i += 2) {
        expect(next.has(hull[i])).toBe(false);
        next.set(hull[i], hull[i + 1]);
    }
    const start = hull[0];
    const cycle = [start];
    let current = next.get(start) as number;
    while (current !== start) {
        cycle.push(current);
        current = next.get(current) as number;
        expect(cycle.length).toBeLessThanOrEqual(next.size);
    }
    expect(cycle.length).toBe(next.size);
    return cycle;
}

describe('Delaunay2', () => {
    it('rejects an empty input', () => {
        expect(() => new Delaunay2().compute([])).toThrow();
    });

    it('reports dimension 0 for coincident points', () => {
        const del = new Delaunay2();
        expect(del.compute([v2(2, 3), v2(2, 3), v2(2, 3)])).toBe(false);
        expect(del.getDimension()).toBe(0);
        expect(del.getLine().origin.equals(v2(2, 3))).toBe(true);
        // Port fix relative to upstream: the vertex count is the input count
        // rather than 0.
        expect(del.getNumVertices()).toBe(3);
        expect(del.getNumTriangles()).toBe(0);
        expect(() => del.getHull()).toThrow();
    });

    it('reports dimension 1 for collinear points', () => {
        const points = [v2(0, 0), v2(1, 2), v2(2, 4), v2(3, 6)];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(false);
        expect(del.getDimension()).toBe(1);
        expect(del.getNumVertices()).toBe(4);
        const line = del.getLine();
        expect(line.origin.equals(v2(0, 0))).toBe(true);
        const s = Math.sqrt(5);
        expect(line.direction.values[0]).toBeCloseTo(1 / s, 12);
        expect(line.direction.values[1]).toBeCloseTo(2 / s, 12);
    });

    it('triangulates a single triangle', () => {
        const points = [v2(0, 0), v2(1, 0), v2(0, 1)];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getDimension()).toBe(2);
        expect(del.getNumTriangles()).toBe(1);
        expect(del.getNumUniqueVertices()).toBe(3);
        expect(del.getAdjacencies()).toEqual([-1, -1, -1]);
        checkTriangulation(del, points);
        checkEmptyCircumcircle(del, points);
        expect(hullCycle(del).length).toBe(3);
    });

    it('triangulates a clockwise triangle into a counterclockwise one', () => {
        const points = [v2(0, 0), v2(0, 1), v2(1, 0)];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumTriangles()).toBe(1);
        checkTriangulation(del, points);
    });

    it('triangulates a unit square (two triangles) with a valid diagonal', () => {
        const points = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumTriangles()).toBe(2);
        checkTriangulation(del, points);
        // The four points are cocircular, so both diagonals are Delaunay.
        checkEmptyCircumcircle(del, points, 1e-9);
        expect(hullCycle(del).length).toBe(4);
    });

    it('triangulates a square with its center: four triangles', () => {
        const points = [v2(-1, -1), v2(1, -1), v2(1, 1), v2(-1, 1), v2(0, 0)];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumTriangles()).toBe(4);
        checkTriangulation(del, points);
        checkEmptyCircumcircle(del, points);
        expect(hullCycle(del).length).toBe(4);
    });

    it('triangulates a 5x5 grid: 2*(n-1)^2 triangles and a 16-vertex hull', () => {
        const points: Vector[] = [];
        for (let y = 0; y < 5; ++y) {
            for (let x = 0; x < 5; ++x) {
                points.push(v2(x, y));
            }
        }
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumUniqueVertices()).toBe(25);
        // Euler: for a triangulation of n points with h hull points and no
        // interior collinear degeneracies, numTriangles = 2n - h - 2.
        // n = 25 and the hull is the square boundary with 16 points.
        expect(del.getNumTriangles()).toBe(2 * 25 - 16 - 2);
        checkTriangulation(del, points);
        checkEmptyCircumcircle(del, points, 1e-9);
        expect(hullCycle(del).length).toBe(16);
    });

    it('handles duplicate points and reports them via getDuplicates', () => {
        const points = [
            v2(0, 0), v2(4, 0), v2(0, 4),
            v2(0, 0),           // duplicate of 0
            v2(1, 1),
            v2(4, 0),           // duplicate of 1
            v2(1, 1)            // duplicate of 4
        ];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumVertices()).toBe(7);
        expect(del.getNumUniqueVertices()).toBe(4);
        const dup = del.getDuplicates();
        expect(dup[0]).toBe(0);
        expect(dup[1]).toBe(1);
        expect(dup[2]).toBe(2);
        expect(dup[3]).toBe(0);
        expect(dup[4]).toBe(4);
        expect(dup[5]).toBe(1);
        expect(dup[6]).toBe(4);
        checkTriangulation(del, points);
        checkEmptyCircumcircle(del, points);
        // No duplicate index appears in the triangulation.
        for (const i of del.getIndices()) {
            expect(dup[i]).toBe(i);
        }
    });

    it('triangulates points on a circle plus the center', () => {
        const n = 12;
        const points: Vector[] = [v2(0, 0)];
        for (let i = 0; i < n; ++i) {
            const t = (2 * Math.PI * i) / n;
            points.push(v2(Math.cos(t), Math.sin(t)));
        }
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumTriangles()).toBe(2 * (n + 1) - n - 2);
        checkTriangulation(del, points);
        checkEmptyCircumcircle(del, points, 1e-9);
        expect(hullCycle(del).length).toBe(n);
    });

    it('triangulates a near-degenerate configuration exercising exact arithmetic', () => {
        // A cluster of points that are cocircular or nearly cocircular to
        // within a few ulps, plus points nearly collinear with an existing
        // edge. The interval arithmetic cannot resolve these signs, so the
        // exact BSNumber fallback is exercised.
        const e = Number.EPSILON;
        const points = [
            v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1),        // cocircular square
            v2(0.5, 0.5),                                   // exact center
            v2(0.5 + e, 0.5),
            v2(0.5, 0.5 + 2 * e),
            v2(0.5 - 3 * e, 0.5 - e),
            v2(0.25, 0.25 + e),                             // near diagonal
            v2(0.75, 0.75 - e)
        ];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getNumUniqueVertices()).toBe(10);
        checkTriangulation(del, points);
        // Every input point is a vertex of the triangulation.
        const used = new Set(del.getIndices());
        for (let i = 0; i < points.length; ++i) {
            expect(used.has(i)).toBe(true);
        }
    });

    it('triangulates points with widely varying magnitudes', () => {
        const points = [
            v2(0, 0), v2(1e8, 0), v2(0, 1e8),
            v2(1e-8, 1e-8), v2(1, 1), v2(1e4, 1e4), v2(3e7, 1e7)
        ];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        checkTriangulation(del, points);
        const used = new Set(del.getIndices());
        expect(used.size).toBe(points.length);
    });

    it('satisfies the empty-circumcircle property on random point sets', () => {
        const rnd = makeRandom(1234567);
        for (let trial = 0; trial < 6; ++trial) {
            const n = 20 + Math.floor(rnd() * 40);
            const points: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                points.push(v2(rnd() * 2 - 1, rnd() * 2 - 1));
            }
            const del = new Delaunay2();
            expect(del.compute(points)).toBe(true);
            expect(del.getNumUniqueVertices()).toBe(n);
            checkTriangulation(del, points);
            checkEmptyCircumcircle(del, points, 1e-12);

            // Euler's relation: 2n - h - 2 triangles for points in general
            // position.
            const h = hullCycle(del).length;
            expect(del.getNumTriangles()).toBe(2 * n - h - 2);
        }
    });

    it('satisfies the empty-circumcircle property on random integer point sets', () => {
        // Small integer coordinates produce many cocircular and collinear
        // configurations, which drive the exact-arithmetic path.
        const rnd = makeRandom(2468013);
        for (let trial = 0; trial < 6; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 60; ++i) {
                points.push(v2(Math.floor(rnd() * 9), Math.floor(rnd() * 9)));
            }
            const del = new Delaunay2();
            expect(del.compute(points)).toBe(true);
            checkTriangulation(del, points);
            checkEmptyCircumcircle(del, points, 1e-9);
        }
    });

    it('finds containing triangles and reports -1 outside the hull', () => {
        const points = [v2(-1, -1), v2(1, -1), v2(1, 1), v2(-1, 1), v2(0, 0)];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);

        const info = new Delaunay2SearchInfo();
        // Interior points map to a triangle that contains them.
        const rnd = makeRandom(555);
        for (let k = 0; k < 50; ++k) {
            const p = v2(rnd() * 1.8 - 0.9, rnd() * 1.8 - 0.9);
            const t = del.getContainingTriangle(p, info);
            expect(t).toBeGreaterThanOrEqual(0);
            expect(t).toBeLessThan(del.getNumTriangles());
            const idx = del.getTriangleIndices(t);
            expect(idx).not.toBeNull();
            const [a, b, c] = (idx as [number, number, number])
                .map(i => points[i]) as [Vector, Vector, Vector];
            expect(cross(a, b, p)).toBeGreaterThanOrEqual(-1e-12);
            expect(cross(b, c, p)).toBeGreaterThanOrEqual(-1e-12);
            expect(cross(c, a, p)).toBeGreaterThanOrEqual(-1e-12);
            // The search path is well formed.
            expect(info.numPath).toBeGreaterThan(0);
            expect(info.finalTriangle).toBe(t);
        }

        // Points outside the hull return -1.
        for (const p of [v2(5, 0), v2(0, -5), v2(-3, -3), v2(2, 2)]) {
            expect(del.getContainingTriangle(p, info)).toBe(-1);
        }

        // Reusing finalTriangle as initialTriangle still works.
        info.initialTriangle = info.finalTriangle;
        expect(del.getContainingTriangle(v2(0.1, 0.1), info))
            .toBeGreaterThanOrEqual(0);

        // The search requires dimension 2.
        const collinear = new Delaunay2();
        collinear.compute([v2(0, 0), v2(1, 1), v2(2, 2)]);
        expect(() => collinear.getContainingTriangle(v2(0, 1),
            new Delaunay2SearchInfo())).toThrow();
    });

    it('exposes graph, indices and adjacencies consistently', () => {
        const points = [v2(0, 0), v2(2, 0), v2(2, 2), v2(0, 2), v2(1, 0.5)];
        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);

        const graph = del.getGraph();
        expect(graph.getNumTriangles()).toBe(del.getNumTriangles());

        for (let t = 0; t < del.getNumTriangles(); ++t) {
            const idx = del.getTriangleIndices(t) as [number, number, number];
            const adj = del.getTriangleAdjacencies(t) as [number, number, number];
            expect(idx).toEqual([del.getIndices()[3 * t], del.getIndices()[3 * t + 1],
                del.getIndices()[3 * t + 2]]);
            expect(adj).toEqual([del.getAdjacencies()[3 * t],
                del.getAdjacencies()[3 * t + 1], del.getAdjacencies()[3 * t + 2]]);
            expect(graph.getTriangle(idx[0], idx[1], idx[2])).not.toBeNull();
        }

        expect(del.getTriangleIndices(-1)).toBeNull();
        expect(del.getTriangleIndices(del.getNumTriangles())).toBeNull();
        expect(del.getTriangleAdjacencies(del.getNumTriangles())).toBeNull();
        expect(del.getVertices().length).toBe(5);
    });

    it('reuses the functor for multiple data sets', () => {
        const del = new Delaunay2();
        const square = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
        expect(del.compute(square)).toBe(true);
        expect(del.getNumTriangles()).toBe(2);

        expect(del.compute([v2(0, 0), v2(1, 1), v2(2, 2)])).toBe(false);
        expect(del.getDimension()).toBe(1);
        expect(del.getNumTriangles()).toBe(0);

        expect(del.compute(square)).toBe(true);
        expect(del.getDimension()).toBe(2);
        expect(del.getNumTriangles()).toBe(2);
        checkTriangulation(del, square);
    });

    it('is deterministic across repeated runs of the same input', () => {
        const rnd = makeRandom(31337);
        const points: Vector[] = [];
        for (let i = 0; i < 40; ++i) {
            points.push(v2(rnd() * 10, rnd() * 10));
        }
        const a = new Delaunay2();
        const b = new Delaunay2();
        expect(a.compute(points)).toBe(true);
        expect(b.compute(points)).toBe(true);
        expect(a.getIndices()).toEqual(b.getIndices());
        expect(a.getAdjacencies()).toEqual(b.getAdjacencies());
    });
});
