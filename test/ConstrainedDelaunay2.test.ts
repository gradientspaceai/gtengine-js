import { describe, it, expect } from 'vitest';
import { ConstrainedDelaunay2 } from '../src/ConstrainedDelaunay2';
import { Vector } from '../src/Vector';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function signedArea(a: Vector, b: Vector, c: Vector): number {
    return 0.5 * ((b.values[0] - a.values[0]) * (c.values[1] - a.values[1])
        - (b.values[1] - a.values[1]) * (c.values[0] - a.values[0]));
}

// The total signed area of the graph triangles; also asserts that every
// triangle is counterclockwise (positive area), which means the triangulation
// is a valid nonoverlapping tiling of the convex hull.
function totalArea(cdt: ConstrainedDelaunay2): number {
    const vertices = cdt.getVertices();
    let sum = 0;
    for (const tri of cdt.getGraph().getTriangles()) {
        const area = signedArea(vertices[tri.V[0]], vertices[tri.V[1]],
            vertices[tri.V[2]]);
        expect(area).toBeGreaterThan(0);
        sum += area;
    }
    return sum;
}

function triangleKey(v0: number, v1: number, v2: number): string {
    return [v0, v1, v2].sort((a, b) => a - b).join(',');
}

function hasEdge(cdt: ConstrainedDelaunay2, v0: number, v1: number): boolean {
    return cdt.getGraph().getEdge(v0, v1) !== null;
}

// A regular grid of (n x n) points with unit spacing; index (i, j) is
// i + n * j.
function gridPoints(n: number): Vector[] {
    const points: Vector[] = [];
    for (let j = 0; j < n; ++j) {
        for (let i = 0; i < n; ++i) {
            points.push(v2(i, j));
        }
    }
    return points;
}

// The Delaunay empty-circumcircle test for a single triangle against a single
// point, using doubles with a tolerance.
function inCircumcircle(a: Vector, b: Vector, c: Vector, p: Vector): boolean {
    const ax = a.values[0] - p.values[0], ay = a.values[1] - p.values[1];
    const bx = b.values[0] - p.values[0], by = b.values[1] - p.values[1];
    const cx = c.values[0] - p.values[0], cy = c.values[1] - p.values[1];
    const det =
        (ax * ax + ay * ay) * (bx * cy - by * cx)
        - (bx * bx + by * by) * (ax * cy - ay * cx)
        + (cx * cx + cy * cy) * (ax * by - ay * bx);
    // The triangles are counterclockwise, so a positive determinant means the
    // point is strictly inside the circumcircle.
    return det > 1e-9;
}

describe('ConstrainedDelaunay2', () => {
    it('computes the unconstrained Delaunay triangulation', () => {
        const points = gridPoints(4);
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        expect(cdt.getDimension()).toBe(2);
        expect(cdt.getNumVertices()).toBe(16);
        // The convex hull is a 3x3 square.
        expect(totalArea(cdt)).toBeCloseTo(9, 12);
        expect(cdt.getInsertedEdges().length).toBe(0);
    });

    it('inserts an edge that crosses a triangle strip', () => {
        // A 4x4 grid. The segment from (0,0) to (3,1) passes through no other
        // grid vertex, so no partitioning occurs.
        const points = gridPoints(4);
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        const numTrianglesBefore = cdt.getGraph().getNumTriangles();
        const areaBefore = totalArea(cdt);

        const e0 = 0;          // (0, 0)
        const e1 = 3 + 4 * 1;  // (3, 1)
        const partitioned = cdt.insert([e0, e1]);

        expect(partitioned[0]).toBe(e0);
        expect(partitioned[partitioned.length - 1]).toBe(e1);
        expect(partitioned).toEqual([e0, e1]);

        // The constrained edge is now an edge of the triangulation.
        expect(hasEdge(cdt, e0, e1)).toBe(true);

        // The retriangulation preserves the tiling of the convex hull and the
        // triangle count.
        expect(totalArea(cdt)).toBeCloseTo(areaBefore, 12);
        expect(cdt.getGraph().getNumTriangles()).toBe(numTrianglesBefore);

        const inserted = cdt.getInsertedEdges();
        expect(inserted.length).toBe(1);
        expect([inserted[0].V[0], inserted[0].V[1]])
            .toEqual([Math.min(e0, e1), Math.max(e0, e1)]);
    });

    it('partitions an edge that passes through interior vertices', () => {
        // The main diagonal of a 4x4 grid passes through the grid vertices
        // (1,1) and (2,2), so the edge is partitioned into 3 subedges.
        const points = gridPoints(4);
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        const areaBefore = totalArea(cdt);

        const e0 = 0;               // (0, 0)
        const e1 = 3 + 4 * 3;       // (3, 3)
        const m1 = 1 + 4 * 1;       // (1, 1)
        const m2 = 2 + 4 * 2;       // (2, 2)
        const partitioned = cdt.insert([e0, e1]);
        expect(partitioned).toEqual([e0, m1, m2, e1]);

        expect(hasEdge(cdt, e0, m1)).toBe(true);
        expect(hasEdge(cdt, m1, m2)).toBe(true);
        expect(hasEdge(cdt, m2, e1)).toBe(true);

        const inserted = cdt.getInsertedEdges();
        expect(inserted.length).toBe(3);
        expect(inserted.map(k => [k.V[0], k.V[1]])).toEqual([
            [e0, m1], [m1, m2], [m2, e1]
        ]);

        expect(totalArea(cdt)).toBeCloseTo(areaBefore, 12);
    });

    it('returns the input edge when it is already in the triangulation', () => {
        const points = gridPoints(3);
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        const areaBefore = totalArea(cdt);
        const numTrianglesBefore = cdt.getGraph().getNumTriangles();

        // (0,0)-(1,0) is a boundary edge of the grid triangulation.
        const partitioned = cdt.insert([0, 1]);
        expect(partitioned).toEqual([0, 1]);
        expect(cdt.getGraph().getNumTriangles()).toBe(numTrianglesBefore);
        expect(totalArea(cdt)).toBeCloseTo(areaBefore, 12);
        const inserted = cdt.getInsertedEdges();
        expect(inserted.length).toBe(1);
        expect([inserted[0].V[0], inserted[0].V[1]]).toEqual([0, 1]);
    });

    it('forces a non-Delaunay diagonal of a quadrilateral', () => {
        // The Delaunay triangulation of this convex quadrilateral uses the
        // short diagonal <0,2>; the constrained triangulation must use the
        // long diagonal <1,3>.
        const points = [v2(0, 0), v2(1, -3), v2(2, 0), v2(1, 3)];
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        expect(cdt.getGraph().getNumTriangles()).toBe(2);
        expect(hasEdge(cdt, 0, 2)).toBe(true);
        expect(hasEdge(cdt, 1, 3)).toBe(false);

        const partitioned = cdt.insert([1, 3]);
        expect(partitioned).toEqual([1, 3]);
        expect(hasEdge(cdt, 1, 3)).toBe(true);
        expect(hasEdge(cdt, 0, 2)).toBe(false);
        expect(cdt.getGraph().getNumTriangles()).toBe(2);
        expect(totalArea(cdt)).toBeCloseTo(6, 12);
    });

    it('inserts several constrained edges and keeps them all', () => {
        const points = gridPoints(5);
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        const areaBefore = totalArea(cdt);

        // Edges chosen so that no two of them cross in their interiors. The
        // header documents interior crossings as ill-posed input.
        const constraints: [number, number][] = [
            [0, 4 + 5 * 1],          // (0,0) - (4,1)
            [5 * 2, 4 + 5 * 2],      // (0,2) - (4,2), through 3 grid vertices
            [5 * 3, 4 + 5 * 4]       // (0,3) - (4,4)
        ];
        for (const edge of constraints) {
            const partitioned = cdt.insert(edge);
            expect(partitioned[0]).toBe(edge[0]);
            expect(partitioned[partitioned.length - 1]).toBe(edge[1]);
            for (let i = 0; i + 1 < partitioned.length; ++i) {
                expect(hasEdge(cdt, partitioned[i], partitioned[i + 1])).toBe(true);
            }
        }

        // All constrained subedges are still present after later insertions.
        for (const key of cdt.getInsertedEdges()) {
            expect(hasEdge(cdt, key.V[0], key.V[1])).toBe(true);
        }
        // The middle constraint was partitioned at the 3 interior vertices.
        expect(cdt.getInsertedEdges().length).toBe(1 + 4 + 1);
        expect(totalArea(cdt)).toBeCloseTo(areaBefore, 12);

        // The compact arrays can be refreshed from the constrained graph.
        cdt.updateIndicesAdjacencies();
        expect(cdt.getNumTriangles()).toBe(cdt.getGraph().getNumTriangles());
        expect(cdt.getIndices().length).toBe(3 * cdt.getNumTriangles());
    });

    it('keeps the Delaunay property away from the constrained edges', () => {
        // Random points; insert one constrained edge between two of them and
        // verify that every triangle whose edges are all unconstrained is
        // still locally Delaunay.
        const random = makeRandom(424242);
        const points: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            points.push(v2(random() * 4, random() * 4));
        }
        // Two well-separated points to constrain.
        points.push(v2(0.05, 2));
        points.push(v2(3.95, 2));
        const e0 = points.length - 2;
        const e1 = points.length - 1;

        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        const areaBefore = totalArea(cdt);

        // The triangles of the unconstrained Delaunay triangulation.
        const before = new Set<string>();
        for (const tri of cdt.getGraph().getTriangles()) {
            before.add(triangleKey(tri.V[0], tri.V[1], tri.V[2]));
        }

        const partitioned = cdt.insert([e0, e1]);
        expect(partitioned[0]).toBe(e0);
        expect(partitioned[partitioned.length - 1]).toBe(e1);
        for (let i = 0; i + 1 < partitioned.length; ++i) {
            expect(hasEdge(cdt, partitioned[i], partitioned[i + 1])).toBe(true);
        }
        expect(totalArea(cdt)).toBeCloseTo(areaBefore, 10);

        // Local Delaunay test away from the constraint: any pair of adjacent
        // triangles that both survived the retriangulation of the triangle
        // strip is still locally Delaunay. (Inside the strip, upstream
        // retriangulates with a closest-vertex bisection rather than a
        // Delaunay rule, so the strip triangles need not be Delaunay.)
        const vertices = cdt.getVertices();
        let numChecked = 0;
        for (const tri of cdt.getGraph().getTriangles()) {
            const tkey = triangleKey(tri.V[0], tri.V[1], tri.V[2]);
            if (!before.has(tkey)) {
                continue;
            }
            for (let i = 0; i < 3; ++i) {
                const adj = tri.T[i];
                if (adj === null) {
                    continue;
                }
                if (!before.has(triangleKey(adj.V[0], adj.V[1], adj.V[2]))) {
                    continue;
                }
                const opposite = adj.getOppositeVertexOfEdge(
                    tri.V[i], tri.V[(i + 1) % 3]);
                if (!opposite.found) {
                    continue;
                }
                ++numChecked;
                expect(inCircumcircle(vertices[tri.V[0]], vertices[tri.V[1]],
                    vertices[tri.V[2]], vertices[opposite.uOpposite])).toBe(false);
            }
        }
        expect(numChecked).toBeGreaterThan(0);
    });

    it('maps duplicate vertices to their representatives', () => {
        // Index 4 duplicates index 0 and index 5 duplicates index 2.
        const points = [v2(0, 0), v2(2, 0), v2(2, 2), v2(0, 2), v2(0, 0), v2(2, 2)];
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        expect(cdt.getNumVertices()).toBe(6);
        expect(cdt.getNumUniqueVertices()).toBe(4);
        expect(cdt.getDuplicates()[4]).toBe(0);
        expect(cdt.getDuplicates()[5]).toBe(2);

        // Inserting <4,5> is the same as inserting <0,2>.
        const partitioned = cdt.insert([4, 5]);
        expect(partitioned).toEqual([0, 2]);
        expect(hasEdge(cdt, 0, 2)).toBe(true);
        const inserted = cdt.getInsertedEdges();
        expect(inserted.length).toBe(1);
        expect([inserted[0].V[0], inserted[0].V[1]]).toEqual([0, 2]);
    });

    it('rejects invalid edges', () => {
        const points = gridPoints(3);
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(points)).toBe(true);
        expect(() => cdt.insert([2, 2])).toThrow();
        expect(() => cdt.insert([-1, 2])).toThrow();
        expect(() => cdt.insert([0, 9])).toThrow();
        // Duplicate vertices collapse to the same representative.
        const dup = new ConstrainedDelaunay2();
        expect(dup.compute([v2(0, 0), v2(1, 0), v2(0, 1), v2(0, 0)])).toBe(true);
        expect(() => dup.insert([0, 3])).toThrow();
    });

    it('rejects insertion when the intrinsic dimension is less than 2', () => {
        const collinear = new ConstrainedDelaunay2();
        expect(collinear.compute([v2(0, 0), v2(1, 1), v2(2, 2)])).toBe(false);
        expect(collinear.getDimension()).toBe(1);
        expect(() => collinear.insert([0, 2])).toThrow();
    });

    it('clears the inserted edges when recomputing (upstream bug fix)', () => {
        const cdt = new ConstrainedDelaunay2();
        expect(cdt.compute(gridPoints(4))).toBe(true);
        cdt.insert([0, 3 + 4 * 1]);
        expect(cdt.getInsertedEdges().length).toBe(1);

        // Upstream leaves mInsertedEdges populated across data sets; the port
        // clears it.
        expect(cdt.compute(gridPoints(3))).toBe(true);
        expect(cdt.getInsertedEdges().length).toBe(0);
        expect(totalArea(cdt)).toBeCloseTo(4, 12);
    });

    it('inserts random constrained edges into random point sets', () => {
        const random = makeRandom(20260902);
        for (let trial = 0; trial < 4; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 20; ++i) {
                points.push(v2(random() * 10, random() * 10));
            }
            const cdt = new ConstrainedDelaunay2();
            expect(cdt.compute(points)).toBe(true);
            const areaBefore = totalArea(cdt);

            // A single random constrained edge; a second edge could cross it
            // in its interior, which the header documents as ill-posed.
            let e0 = Math.floor(random() * points.length);
            let e1 = Math.floor(random() * points.length);
            if (e0 === e1) {
                e1 = (e1 + 1) % points.length;
            }
            e0 = cdt.getDuplicates()[e0];
            e1 = cdt.getDuplicates()[e1];
            if (e0 === e1) {
                continue;
            }

            const partitioned = cdt.insert([e0, e1]);
            expect(partitioned[0]).toBe(e0);
            expect(partitioned[partitioned.length - 1]).toBe(e1);
            for (let i = 0; i + 1 < partitioned.length; ++i) {
                expect(hasEdge(cdt, partitioned[i], partitioned[i + 1])).toBe(true);
            }
            // The triangulation still tiles the convex hull with
            // counterclockwise triangles (totalArea asserts the latter).
            expect(totalArea(cdt)).toBeCloseTo(areaBefore, 9);
        }
    });
});
