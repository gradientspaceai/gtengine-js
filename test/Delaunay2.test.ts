import { describe, it, expect } from 'vitest';
import { Delaunay2, Delaunay2SearchInfo } from '../src/Delaunay2.js';
import { Vector } from '../src/Vector.js';
import { SWInterval } from '../src/SWInterval.js';
import { IntrinsicsVector2 } from '../src/Vector2.js';
import { check, fc, latticeVector, wellScaledVector } from './helpers/arbitraries.js';
import { exactDyadic, inCircle2, orient2, twiceSignedAreaExact } from './helpers/exact.js';

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

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks against exact
// bigint arithmetic. Input coordinates are finite doubles, so the coordinates
// of one point set scale by a common power of two into exact integers
// (exactDyadic); the orientation and in-circle determinants are homogeneous,
// so that scale never changes their signs. The reference answers are
// therefore the mathematically exact ones, which is what the interval +
// BSNumber predicates of Delaunay2 are required to reproduce.
// ---------------------------------------------------------------------------

// Exposes the protected sign predicates so they can be compared against the
// exact bigint reference predicates.
class Delaunay2Probe extends Delaunay2 {
    line(pIndex: number, v0Index: number, v1Index: number): number {
        return this.toLine(pIndex, v0Index, v1Index);
    }

    circumcircle(pIndex: number, v0Index: number, v1Index: number,
        v2Index: number): number {
        return this.toCircumcircle(pIndex, v0Index, v1Index, v2Index);
    }
}

function exactOf(points: readonly Vector[]): bigint[][] {
    const flat: number[] = [];
    for (const p of points) {
        flat.push(p.values[0], p.values[1]);
    }
    const s = exactDyadic(flat);
    return points.map((_, i) => [s[2 * i], s[2 * i + 1]]);
}

const keyOfExact = (p: readonly bigint[]): string =>
    String(p[0]) + ',' + String(p[1]);

// Andrew monotone chain on exact integer coordinates, returning indices into
// 'pts'. One index when all points coincide, two (the extremes) when they are
// collinear, otherwise the strictly convex counterclockwise hull.
function exactHullIndices(pts: readonly (readonly bigint[])[]): number[] {
    const idx = pts.map((_, i) => i);
    const cmp = (a: number, b: number): number =>
        pts[a][0] < pts[b][0] ? -1 : pts[a][0] > pts[b][0] ? 1
            : pts[a][1] < pts[b][1] ? -1 : pts[a][1] > pts[b][1] ? 1 : 0;
    idx.sort(cmp);
    const uniq: number[] = [];
    for (const i of idx) {
        if (uniq.length === 0 || cmp(uniq[uniq.length - 1], i) !== 0) {
            uniq.push(i);
        }
    }
    if (uniq.length < 3) {
        return uniq;
    }
    const turn = (o: number, a: number, b: number): number =>
        orient2(pts[o][0], pts[o][1], pts[a][0], pts[a][1], pts[b][0], pts[b][1]);
    const build = (src: readonly number[]): number[] => {
        const out: number[] = [];
        for (const p of src) {
            while (out.length >= 2
                && turn(out[out.length - 2], out[out.length - 1], p) <= 0) {
                out.pop();
            }
            out.push(p);
        }
        out.pop();
        return out;
    };
    return build(uniq).concat(build(uniq.slice().reverse()));
}

// Every observable of a successful triangulation, checked exactly.
function verifyTriangulation(del: Delaunay2, points: readonly Vector[]): void {
    const exact = exactOf(points);
    const indices = del.getIndices();
    const adjacencies = del.getAdjacencies();
    const numTriangles = del.getNumTriangles();
    const duplicates = del.getDuplicates();

    expect(del.getNumVertices()).toBe(points.length);
    expect(indices.length).toBe(3 * numTriangles);
    expect(adjacencies.length).toBe(3 * numTriangles);
    expect(numTriangles).toBeGreaterThan(0);

    // getDuplicates(): the extreme vertices chosen by IntrinsicsVector2 are
    // always the first occurrence of their coordinates (the min/max/distance
    // scans use strict inequalities), so duplicates[i] is the smallest index
    // carrying the same coordinates as vertex i.
    const firstIndex = new Map<string, number>();
    for (let i = 0; i < points.length; ++i) {
        const key = keyOfExact(exact[i]);
        if (!firstIndex.has(key)) {
            firstIndex.set(key, i);
        }
        expect(duplicates[i]).toBe(firstIndex.get(key));
    }
    expect(del.getNumUniqueVertices()).toBe(firstIndex.size);

    // The mesh vertices are exactly the non-duplicate input indices.
    const used = new Set<number>(indices);
    const expectedUsed = new Set<number>();
    for (let i = 0; i < points.length; ++i) {
        if (duplicates[i] === i) {
            expectedUsed.add(i);
        }
    }
    expect(Array.from(used).sort((a, b) => a - b))
        .toEqual(Array.from(expectedUsed).sort((a, b) => a - b));

    // Every triangle is exactly counterclockwise (so nondegenerate).
    for (let t = 0; t < numTriangles; ++t) {
        const a = exact[indices[3 * t]];
        const b = exact[indices[3 * t + 1]];
        const c = exact[indices[3 * t + 2]];
        expect(orient2(a[0], a[1], b[0], b[1], c[0], c[1])).toBe(1);
    }

    // Each directed edge occurs at most once, and adjacency is the symmetric
    // relation induced by shared edges.
    const directed = new Map<string, number>();
    for (let t = 0; t < numTriangles; ++t) {
        for (let j = 0; j < 3; ++j) {
            const key = indices[3 * t + j] + ',' + indices[3 * t + ((j + 1) % 3)];
            expect(directed.has(key)).toBe(false);
            directed.set(key, t);
        }
    }
    for (let t = 0; t < numTriangles; ++t) {
        for (let j = 0; j < 3; ++j) {
            const v0 = indices[3 * t + j];
            const v1 = indices[3 * t + ((j + 1) % 3)];
            const opposite = directed.get(v1 + ',' + v0);
            const adj = adjacencies[3 * t + j];
            expect(adj).toBe(opposite === undefined ? -1 : opposite);
        }
    }

    // The exact empty-circumcircle property: no unique vertex is strictly
    // inside the circumcircle of any triangle.
    for (let t = 0; t < numTriangles; ++t) {
        const i0 = indices[3 * t], i1 = indices[3 * t + 1], i2 = indices[3 * t + 2];
        const a = exact[i0], b = exact[i1], c = exact[i2];
        for (let i = 0; i < points.length; ++i) {
            if (duplicates[i] !== i || i === i0 || i === i1 || i === i2) {
                continue;
            }
            const d = exact[i];
            expect(inCircle2(a[0], a[1], b[0], b[1], c[0], c[1], d[0], d[1]))
                .toBeLessThanOrEqual(0);
        }
    }

    // The triangles tile the convex hull: the exact areas agree. Combined
    // with the counterclockwise orientation and the once-only directed edges,
    // this rules out both gaps and overlaps.
    let twiceArea = 0n;
    for (let t = 0; t < numTriangles; ++t) {
        const a = exact[indices[3 * t]];
        const b = exact[indices[3 * t + 1]];
        const c = exact[indices[3 * t + 2]];
        twiceArea += (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    }
    const hull = exactHullIndices(exact);
    const hullPoly = hull.map(i => exact[i]);
    expect(twiceSignedAreaExact(hullPoly)).toBe(twiceArea);

    // getHull(): the boundary edges form one counterclockwise convex cycle
    // enclosing the same area. It may contain more vertices than the strictly
    // convex hull, because input points lying on a hull edge are vertices of
    // the triangulation.
    const boundary = del.getHull();
    expect(boundary.length % 2).toBe(0);
    const next = new Map<number, number>();
    for (let i = 0; i < boundary.length; i += 2) {
        expect(next.has(boundary[i])).toBe(false);
        next.set(boundary[i], boundary[i + 1]);
    }
    const cycle: number[] = [boundary[0]];
    for (;;) {
        const n = next.get(cycle[cycle.length - 1]);
        expect(n).not.toBe(undefined);
        if (n === boundary[0]) {
            break;
        }
        cycle.push(n as number);
        expect(cycle.length).toBeLessThanOrEqual(next.size);
    }
    expect(cycle.length).toBe(next.size);
    const cyclePoly = cycle.map(i => exact[i]);
    expect(twiceSignedAreaExact(cyclePoly)).toBe(twiceArea);
    for (let i = 0; i < cyclePoly.length; ++i) {
        const a = cyclePoly[i];
        const b = cyclePoly[(i + 1) % cyclePoly.length];
        const c = cyclePoly[(i + 2) % cyclePoly.length];
        // Convex: a left turn, or collinear where a hull edge carries an
        // interior input point.
        expect(orient2(a[0], a[1], b[0], b[1], c[0], c[1]))
            .toBeGreaterThanOrEqual(0);
    }
    // Every strictly convex hull vertex is on the boundary cycle.
    const onCycle = new Set<string>(cycle.map(i => keyOfExact(exact[i])));
    for (const i of hull) {
        expect(onCycle.has(keyOfExact(exact[i]))).toBe(true);
    }
}

// Run Delaunay2 and check the dimension classification and, when the
// dimension is 2, the whole triangulation.
function verifyDelaunay(points: Vector[]): void {
    const exact = exactOf(points);
    const hull = exactHullIndices(exact);
    const expectedDimension = hull.length === 1 ? 0 : hull.length === 2 ? 1 : 2;

    const del = new Delaunay2();
    const is2D = del.compute(points);
    expect(is2D).toBe(expectedDimension === 2);
    expect(del.getDimension()).toBe(expectedDimension);
    expect(del.getNumVertices()).toBe(points.length);
    if (!is2D) {
        expect(del.getNumTriangles()).toBe(0);
        return;
    }
    verifyTriangulation(del, points);
}

describe('Delaunay2 verification', () => {
    it('reproduces the exact Delaunay triangulation of lattice point sets', () => {
        check(fc.array(latticeVector(2, -6, 6), { minLength: 1, maxLength: 14 }),
            (points) => {
                verifyDelaunay(points);
            }, 120);
    }, 30000);

    it('reproduces the exact Delaunay triangulation of clustered lattice points', () => {
        // A coarse lattice makes duplicates, collinear triples and cocircular
        // quadruples common, which is where the exact predicates matter.
        check(fc.array(latticeVector(2, -2, 2), { minLength: 1, maxLength: 14 }),
            (points) => {
                verifyDelaunay(points);
            }, 120);
    }, 30000);

    it('reproduces the exact Delaunay triangulation of well-scaled real points', () => {
        check(fc.array(wellScaledVector(2, -5, 5), { minLength: 1, maxLength: 12 }),
            (points) => {
                verifyDelaunay(points);
            }, 120);
    }, 30000);

    it('classifies exactly collinear input as dimension 0 or 1', () => {
        const arb = fc.tuple(
            latticeVector(2, -10, 10),
            latticeVector(2, -5, 5).filter(d => d.values[0] !== 0 || d.values[1] !== 0),
            fc.array(fc.integer({ min: -20, max: 20 }), { minLength: 1, maxLength: 10 }));
        check(arb, ([base, dir, ts]) => {
            const points = ts.map(t => Vector.fromArray([
                base.values[0] + t * dir.values[0],
                base.values[1] + t * dir.values[1]]));
            const del = new Delaunay2();
            expect(del.compute(points)).toBe(false);
            expect(del.getDimension()).toBe(new Set(ts).size === 1 ? 0 : 1);
            expect(del.getNumTriangles()).toBe(0);
            expect(del.getNumVertices()).toBe(points.length);
        }, 200);
    });

    it('toLine and toCircumcircle agree with the exact bigint predicates', () => {
        // The lattice is coarse so that collinear triples and cocircular
        // quadruples -- the cases where the SWInterval bound straddles zero
        // and the BSNumber fallback decides the sign -- occur constantly.
        const arb = fc.array(latticeVector(2, -3, 3), { minLength: 3, maxLength: 6 });
        check(arb, (points) => {
            const del = new Delaunay2Probe();
            if (!del.compute(points)) {
                return;
            }
            const e = exactOf(points);
            const n = points.length;
            for (let v0 = 0; v0 < n; ++v0) {
                for (let v1 = 0; v1 < n; ++v1) {
                    if (v0 === v1 || keyOfExact(e[v0]) === keyOfExact(e[v1])) {
                        continue;
                    }
                    for (let p = 0; p < n; ++p) {
                        // ToLine returns +1 when P is on the right of the
                        // directed line <V0,V1>, which is the negated 2D
                        // orientation determinant.
                        // '| 0' keeps the negated zero an integer zero.
                        expect(del.line(p, v0, v1)).toBe(
                            -orient2(e[v0][0], e[v0][1], e[v1][0], e[v1][1],
                                e[p][0], e[p][1]) | 0);
                    }
                }
            }
            for (let v0 = 0; v0 < n; ++v0) {
                for (let v1 = 0; v1 < n; ++v1) {
                    for (let v2 = 0; v2 < n; ++v2) {
                        const ccw = orient2(e[v0][0], e[v0][1], e[v1][0], e[v1][1],
                            e[v2][0], e[v2][1]);
                        if (ccw !== 1) {
                            continue;   // documented for CCW triangles only
                        }
                        for (let p = 0; p < n; ++p) {
                            // ToCircumcircle returns +1 outside, -1 inside,
                            // 0 on the circumcircle.
                            expect(del.circumcircle(p, v0, v1, v2)).toBe(
                                -inCircle2(e[v0][0], e[v0][1], e[v1][0], e[v1][1],
                                    e[v2][0], e[v2][1], e[p][0], e[p][1]) | 0);
                        }
                    }
                }
            }
        }, 30);
    }, 30000);

    it('the exact fallback decides signs the interval arithmetic cannot', () => {
        // Four cocircular points: the in-circle determinant is exactly zero,
        // so the SWInterval bound necessarily contains zero and only the
        // BSNumber path can return 0. The same holds for ToLine on three
        // exactly collinear points.
        const points = [v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4), v2(2, 0)];
        const del = new Delaunay2Probe();
        expect(del.compute(points)).toBe(true);

        // Replicate the interval expression tree of ToCircumcircle for
        // (V0,V1,V2) = ((0,0),(4,0),(4,4)) and P = (0,4) to confirm that the
        // interval bound really is indeterminate.
        const iSub = (a: number, b: number): SWInterval => SWInterval.sub(a, b);
        const iAdd = (a: number, b: number): SWInterval => SWInterval.add(a, b);
        const V = [[0, 0], [4, 0], [4, 4]];
        const P = [0, 4];
        const x: SWInterval[] = [], y: SWInterval[] = [], z: SWInterval[] = [];
        for (let i = 0; i < 3; ++i) {
            x.push(iSub(V[i][0], P[0]));
            y.push(iSub(V[i][1], P[1]));
            const s0 = iAdd(V[i][0], P[0]);
            const s1 = iAdd(V[i][1], P[1]);
            z.push(s0.mul(x[i]).add(s1.mul(y[i])));
        }
        const c0 = y[1].mul(z[2]).sub(y[2].mul(z[1]));
        const c1 = y[2].mul(z[0]).sub(y[0].mul(z[2]));
        const c2 = y[0].mul(z[1]).sub(y[1].mul(z[0]));
        const det = x[0].mul(c0).add(x[1].mul(c1)).add(x[2].mul(c2));
        expect(det.get(0)).toBeLessThanOrEqual(0);
        expect(det.get(1)).toBeGreaterThanOrEqual(0);

        // The exact answer: (0,4) is on the circumcircle of the other three.
        expect(del.circumcircle(3, 0, 1, 2)).toBe(0);
        // (2,0) lies exactly on the line through (0,0) and (4,0).
        expect(del.line(4, 0, 1)).toBe(0);

        verifyTriangulation(del, points);
    });

    it('getContainingTriangle agrees with an exact brute-force search', () => {
        const arb = fc.tuple(
            fc.array(latticeVector(2, -4, 4), { minLength: 4, maxLength: 10 }),
            fc.array(fc.tuple(fc.integer({ min: -12, max: 12 }),
                fc.integer({ min: -12, max: 12 })),
            { minLength: 1, maxLength: 6 }));
        check(arb, ([points, rawQueries]) => {
            const del = new Delaunay2();
            if (!del.compute(points)) {
                return;
            }
            // Halved integers so queries land on vertices, on edges and in
            // triangle interiors.
            const queries = rawQueries.map(([a, b]) => v2(a / 2, b / 2));
            const indices = del.getIndices();
            const numTriangles = del.getNumTriangles();
            const all = exactOf(points.concat(queries));
            const e = all.slice(0, points.length);
            const q = all.slice(points.length);

            for (let k = 0; k < queries.length; ++k) {
                const info = new Delaunay2SearchInfo();
                const t = del.getContainingTriangle(queries[k], info);
                const contains = (tri: number): boolean => {
                    for (let j = 0; j < 3; ++j) {
                        const a = e[indices[3 * tri + j]];
                        const b = e[indices[3 * tri + ((j + 1) % 3)]];
                        if (orient2(a[0], a[1], b[0], b[1], q[k][0], q[k][1]) < 0) {
                            return false;
                        }
                    }
                    return true;
                };
                if (t === -1) {
                    // The point is outside every triangle.
                    for (let tri = 0; tri < numTriangles; ++tri) {
                        expect(contains(tri)).toBe(false);
                    }
                } else {
                    expect(t).toBeGreaterThanOrEqual(0);
                    expect(t).toBeLessThan(numTriangles);
                    expect(contains(t)).toBe(true);
                    // The reported path starts at the initial triangle and
                    // ends at the returned one.
                    expect(info.numPath).toBeGreaterThan(0);
                    expect(info.path[0]).toBe(info.initialTriangle);
                    expect(info.path[info.numPath - 1]).toBe(t);
                    expect(info.finalTriangle).toBe(t);
                }
            }
        }, 60);
    }, 30000);

    it('is invariant under permutation of the input (same triangle set)', () => {
        const arb = fc.array(latticeVector(2, -4, 4), { minLength: 3, maxLength: 12 })
            .chain(points => fc.tuple(fc.constant(points),
                fc.shuffledSubarray(points.map((_, i) => i),
                    { minLength: points.length, maxLength: points.length })));
        check(arb, ([points, perm]) => {
            const a = new Delaunay2();
            if (!a.compute(points)) {
                return;
            }
            const permuted = perm.map(i => points[i]);
            const b = new Delaunay2();
            expect(b.compute(permuted)).toBe(true);
            expect(b.getNumUniqueVertices()).toBe(a.getNumUniqueVertices());

            // Compare the triangles as sets of coordinate triples. Delaunay
            // triangulations are unique only when no four points are
            // cocircular, so compare areas always and triangle sets only in
            // the generic case.
            const ea = exactOf(points), eb = exactOf(permuted);
            const asKeys = (del: Delaunay2, e: bigint[][]): string[] => {
                const ind = del.getIndices();
                const out: string[] = [];
                for (let t = 0; t < del.getNumTriangles(); ++t) {
                    const tri = [e[ind[3 * t]], e[ind[3 * t + 1]], e[ind[3 * t + 2]]]
                        .map(keyOfExact).sort();
                    out.push(tri.join('|'));
                }
                return out.sort();
            };
            const ka = asKeys(a, ea), kb = asKeys(b, eb);
            const area = (del: Delaunay2, e: bigint[][]): bigint => {
                const ind = del.getIndices();
                let s = 0n;
                for (let t = 0; t < del.getNumTriangles(); ++t) {
                    const p0 = e[ind[3 * t]], p1 = e[ind[3 * t + 1]],
                        p2 = e[ind[3 * t + 2]];
                    s += (p1[0] - p0[0]) * (p2[1] - p0[1])
                        - (p1[1] - p0[1]) * (p2[0] - p0[0]);
                }
                return s;
            };
            expect(area(b, eb)).toBe(area(a, ea));
            const cocircular = ka.length !== kb.length
                || ka.some((k, i) => k !== kb[i]);
            if (cocircular) {
                // Both must still be valid Delaunay triangulations.
                verifyTriangulation(a, points);
                verifyTriangulation(b, permuted);
            }
        }, 60);
    }, 30000);

    it('reuse does not leak state between data sets', () => {
        const arb = fc.tuple(
            fc.array(latticeVector(2, -4, 4), { minLength: 1, maxLength: 10 }),
            fc.array(latticeVector(2, -4, 4), { minLength: 1, maxLength: 10 }));
        check(arb, ([first, second]) => {
            const shared = new Delaunay2();
            shared.compute(first);
            const sharedResult = shared.compute(second);
            const fresh = new Delaunay2();
            expect(sharedResult).toBe(fresh.compute(second));
            expect(shared.getDimension()).toBe(fresh.getDimension());
            expect(shared.getNumVertices()).toBe(fresh.getNumVertices());
            expect(shared.getNumUniqueVertices()).toBe(fresh.getNumUniqueVertices());
            expect(shared.getNumTriangles()).toBe(fresh.getNumTriangles());
            expect(shared.getIndices()).toEqual(fresh.getIndices());
            expect(shared.getAdjacencies()).toEqual(fresh.getAdjacencies());
            expect(shared.getDuplicates()).toEqual(fresh.getDuplicates());
        }, 60);
    }, 30000);
});

describe('Delaunay2 verification regressions', () => {
    // Upstream bug, fixed in the port. IntrinsicsVector2 measures the
    // perpendicular distances against a normalized direction, so its
    // floating-point roundoff reports dimension 2 for these exactly collinear
    // sets even though Delaunay2<T> passes epsilon = 0. Before the fix,
    // [(-3,-9),(0,0)] produced compute() === true with getIndices() ===
    // [0,1,1] (a degenerate triangle whose second and third vertices are the
    // same), and the longer sets threw 'Unexpected condition.' or 'Unexpected
    // termination of loop while searching for a triangle.' from deep inside
    // the incremental update.
    const collinearCases: number[][][] = [
        [[-3, -9], [0, 0]],
        [[0, 0], [1, 3], [2, 6], [3, 9]],
        [[0, 0], [1, 3], [2, 6], [3, 9], [4, 12]],
        [[0, 0], [3, 1], [6, 2], [9, 3]],
        [[0, 0], [1, 7], [2, 14], [3, 21]]
    ];

    it('the intrinsics computation really does misclassify these sets', () => {
        for (const c of collinearCases) {
            const points = c.map(p => v2(p[0], p[1]));
            expect(new IntrinsicsVector2(points, 0).dimension).toBe(2);
        }
    });

    it('classifies exactly collinear input as dimension 1 despite that', () => {
        for (const c of collinearCases) {
            const points = c.map(p => v2(p[0], p[1]));
            const del = new Delaunay2();
            expect(del.compute(points), JSON.stringify(c)).toBe(false);
            expect(del.getDimension()).toBe(1);
            expect(del.getNumTriangles()).toBe(0);
            expect(del.getIndices().length).toBe(0);
            expect(del.getNumVertices()).toBe(points.length);
            // The reported line contains every input vertex exactly (the
            // points are integer multiples of the direction from the origin).
            const line = del.getLine();
            for (const p of points) {
                const dx = p.values[0] - line.origin.values[0];
                const dy = p.values[1] - line.origin.values[1];
                const perp = dx * line.direction.values[1]
                    - dy * line.direction.values[0];
                expect(Math.abs(perp)).toBeLessThan(1e-12);
            }
        }
    });

    it('repairs a degenerate seed triangle when an off-line vertex exists', () => {
        // (1, 3.0000000000000004) is one ulp off the line through the other
        // three points, which are collinear. The perpendicular distance is
        // below the roundoff of the normalized direction, so the intrinsics
        // computation reports extreme = [0,3,3]: a seed triangle whose second
        // and third vertices are the same. Upstream inserts that degenerate
        // triangle; the port re-selects the third extreme with the exact
        // ToLine predicate and triangulates correctly.
        const points = [v2(0, 0), v2(1, 3.0000000000000004), v2(2, 6), v2(3, 9)];
        const info = new IntrinsicsVector2(points, 0);
        expect(info.dimension).toBe(2);
        expect(info.extreme[1]).toBe(info.extreme[2]);

        const del = new Delaunay2();
        expect(del.compute(points)).toBe(true);
        expect(del.getDimension()).toBe(2);
        expect(del.getNumTriangles()).toBe(2);
        verifyTriangulation(del, points);
    });

    it('toCircumcircle returns integer zero rather than -0', () => {
        // The exact path negates the BSNumber sign; JavaScript's unary minus
        // maps the integer 0 to -0, which is a distinct value under Object.is.
        const points = [v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4)];
        const del = new Delaunay2Probe();
        expect(del.compute(points)).toBe(true);
        const sign = del.circumcircle(3, 0, 1, 2);
        expect(sign).toBe(0);
        expect(Object.is(sign, -0)).toBe(false);
    });
});
