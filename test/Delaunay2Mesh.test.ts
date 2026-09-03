import { describe, it, expect } from 'vitest';
import { Delaunay2 } from '../src/Delaunay2.js';
import { Delaunay2Mesh } from '../src/Delaunay2Mesh.js';
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

function makeTriangulation(points: readonly Vector[]): Delaunay2 {
    const delaunay = new Delaunay2();
    expect(delaunay.compute(points)).toBe(true);
    return delaunay;
}

describe('Delaunay2Mesh', () => {
    it('requires a 2-dimensional triangulation', () => {
        const collinear = new Delaunay2();
        expect(collinear.compute([v2(0, 0), v2(1, 1), v2(2, 2)])).toBe(false);
        expect(() => new Delaunay2Mesh(collinear)).toThrow();
    });

    it('mirrors the accessors of the underlying triangulation', () => {
        const points = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1), v2(0.4, 0.6)];
        const delaunay = makeTriangulation(points);
        const mesh = new Delaunay2Mesh(delaunay);

        expect(mesh.getNumVertices()).toBe(delaunay.getNumVertices());
        expect(mesh.getNumVertices()).toBe(points.length);
        expect(mesh.getNumTriangles()).toBe(delaunay.getNumTriangles());
        expect(mesh.getVertices()).toBe(delaunay.getVertices());
        expect(mesh.getIndices()).toEqual(delaunay.getIndices());
        expect(mesh.getAdjacencies()).toEqual(delaunay.getAdjacencies());
        expect(mesh.getInvalidIndex()).toBe(-1);

        for (let t = 0; t < mesh.getNumTriangles(); ++t) {
            const indices = mesh.getTriangleIndices(t);
            expect(indices).toEqual(delaunay.getTriangleIndices(t));
            expect(indices).not.toBeNull();
            const adjacencies = mesh.getTriangleAdjacencies(t);
            expect(adjacencies).toEqual(delaunay.getTriangleAdjacencies(t));

            const vertices = mesh.getTriangleVertices(t);
            expect(vertices).not.toBeNull();
            for (let i = 0; i < 3; ++i) {
                expect((vertices as Vector[])[i].equals(
                    points[(indices as number[])[i]])).toBe(true);
            }

            // The compact index array agrees with the per-triangle accessor.
            for (let i = 0; i < 3; ++i) {
                expect((indices as number[])[i]).toBe(mesh.getIndices()[3 * t + i]);
                expect((adjacencies as number[])[i])
                    .toBe(mesh.getAdjacencies()[3 * t + i]);
            }
        }
    });

    it('returns null for out-of-range triangle indices', () => {
        const points = [v2(0, 0), v2(1, 0), v2(0, 1)];
        const mesh = new Delaunay2Mesh(makeTriangulation(points));
        expect(mesh.getNumTriangles()).toBe(1);
        for (const t of [-1, 1, 100]) {
            expect(mesh.getTriangleIndices(t)).toBeNull();
            expect(mesh.getTriangleAdjacencies(t)).toBeNull();
            expect(mesh.getTriangleVertices(t)).toBeNull();
            expect(mesh.getBarycentrics(t, v2(0.25, 0.25))).toBeNull();
        }
    });

    it('computes exact barycentric coordinates', () => {
        const points = [v2(0, 0), v2(4, 0), v2(0, 4)];
        const mesh = new Delaunay2Mesh(makeTriangulation(points));
        const indices = mesh.getTriangleIndices(0) as number[];

        // A point whose barycentric coordinates are exact binary fractions.
        const p = v2(1, 2);
        const bary = mesh.getBarycentrics(0, p) as number[];
        expect(bary).not.toBeNull();
        expect(bary[0] + bary[1] + bary[2]).toBe(1);

        // Reconstruct the point from the barycentric coordinates.
        for (let j = 0; j < 2; ++j) {
            let value = 0;
            for (let i = 0; i < 3; ++i) {
                value += bary[i] * points[indices[i]].values[j];
            }
            expect(value).toBeCloseTo(p.values[j], 12);
        }

        // A vertex of the triangle has a unit barycentric coordinate.
        const atVertex = mesh.getBarycentrics(0, points[indices[1]]) as number[];
        expect(atVertex).toEqual([0, 1, 0]);
    });

    it('locates containing triangles and cross-checks barycentrics', () => {
        const random = makeRandom(9091);
        const points: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            points.push(v2(random() * 4 - 2, random() * 4 - 2));
        }
        const delaunay = makeTriangulation(points);
        const mesh = new Delaunay2Mesh(delaunay);

        let numFound = 0;
        for (let trial = 0; trial < 40; ++trial) {
            const p = v2(random() * 4 - 2, random() * 4 - 2);
            const t = mesh.getContainingTriangle(p);
            if (t === mesh.getInvalidIndex()) {
                continue;
            }
            ++numFound;
            expect(t).toBeGreaterThanOrEqual(0);
            expect(t).toBeLessThan(mesh.getNumTriangles());

            const bary = mesh.getBarycentrics(t, p) as number[];
            expect(bary).not.toBeNull();
            // The point is inside the triangle, so the barycentrics are
            // nonnegative and sum to one.
            expect(bary[0] + bary[1] + bary[2]).toBeCloseTo(1, 12);
            for (const b of bary) {
                expect(b).toBeGreaterThanOrEqual(-1e-12);
            }

            // The barycentric combination reproduces the query point.
            const indices = mesh.getTriangleIndices(t) as number[];
            for (let j = 0; j < 2; ++j) {
                let value = 0;
                for (let i = 0; i < 3; ++i) {
                    value += bary[i] * points[indices[i]].values[j];
                }
                expect(value).toBeCloseTo(p.values[j], 10);
            }
        }
        expect(numFound).toBeGreaterThan(0);

        // A point far outside the convex hull is not in any triangle.
        expect(mesh.getContainingTriangle(v2(100, 100)))
            .toBe(mesh.getInvalidIndex());
    });
});
