import { describe, it, expect } from 'vitest';
import { Delaunay3 } from '../src/Delaunay3.js';
import { Delaunay3Mesh } from '../src/Delaunay3Mesh.js';
import { Vector } from '../src/Vector.js';

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

function makeTetrahedralization(points: readonly Vector[]): Delaunay3 {
    const delaunay = new Delaunay3();
    expect(delaunay.compute(points)).toBe(true);
    return delaunay;
}

describe('Delaunay3Mesh', () => {
    it('requires a 3-dimensional tetrahedralization', () => {
        const coplanar = new Delaunay3();
        expect(coplanar.compute(
            [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0)])).toBe(false);
        expect(() => new Delaunay3Mesh(coplanar)).toThrow();
    });

    it('mirrors the accessors of the underlying tetrahedralization', () => {
        const points = [
            v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1),
            v3(1, 1, 1), v3(0.3, 0.3, 0.2)
        ];
        const delaunay = makeTetrahedralization(points);
        const mesh = new Delaunay3Mesh(delaunay);

        expect(mesh.getNumVertices()).toBe(delaunay.getNumVertices());
        expect(mesh.getNumVertices()).toBe(points.length);
        expect(mesh.getNumTetrahedra()).toBe(delaunay.getNumTetrahedra());
        expect(mesh.getNumTetrahedra()).toBeGreaterThan(0);
        expect(mesh.getVertices()).toBe(delaunay.getVertices());
        expect(mesh.getIndices()).toEqual(delaunay.getIndices());
        expect(mesh.getAdjacencies()).toEqual(delaunay.getAdjacencies());
        expect(mesh.getInvalidIndex()).toBe(-1);

        for (let t = 0; t < mesh.getNumTetrahedra(); ++t) {
            const indices = mesh.getTetrahedronIndices(t);
            expect(indices).toEqual(delaunay.getTetrahedronIndices(t));
            expect(indices).not.toBeNull();
            const adjacencies = mesh.getTetrahedronAdjacencies(t);
            expect(adjacencies).toEqual(delaunay.getTetrahedronAdjacencies(t));

            const vertices = mesh.getTetrahedronVertices(t);
            expect(vertices).not.toBeNull();
            for (let i = 0; i < 4; ++i) {
                expect((vertices as Vector[])[i].equals(
                    points[(indices as number[])[i]])).toBe(true);
                expect((indices as number[])[i]).toBe(mesh.getIndices()[4 * t + i]);
                expect((adjacencies as number[])[i])
                    .toBe(mesh.getAdjacencies()[4 * t + i]);
            }
        }
    });

    it('returns null for out-of-range tetrahedron indices', () => {
        const points = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        const mesh = new Delaunay3Mesh(makeTetrahedralization(points));
        expect(mesh.getNumTetrahedra()).toBe(1);
        for (const t of [-1, 1, 100]) {
            expect(mesh.getTetrahedronIndices(t)).toBeNull();
            expect(mesh.getTetrahedronAdjacencies(t)).toBeNull();
            expect(mesh.getTetrahedronVertices(t)).toBeNull();
            expect(mesh.getBarycentrics(t, v3(0.1, 0.1, 0.1))).toBeNull();
        }
    });

    it('computes exact barycentric coordinates', () => {
        const points = [v3(0, 0, 0), v3(4, 0, 0), v3(0, 4, 0), v3(0, 0, 4)];
        const mesh = new Delaunay3Mesh(makeTetrahedralization(points));
        const indices = mesh.getTetrahedronIndices(0) as number[];

        const p = v3(1, 1, 1);
        const bary = mesh.getBarycentrics(0, p) as number[];
        expect(bary).not.toBeNull();
        expect(bary[0] + bary[1] + bary[2] + bary[3]).toBe(1);

        for (let j = 0; j < 3; ++j) {
            let value = 0;
            for (let i = 0; i < 4; ++i) {
                value += bary[i] * points[indices[i]].values[j];
            }
            expect(value).toBeCloseTo(p.values[j], 12);
        }

        const atVertex = mesh.getBarycentrics(0, points[indices[2]]) as number[];
        expect(atVertex).toEqual([0, 0, 1, 0]);
    });

    it('locates containing tetrahedra and cross-checks barycentrics', () => {
        const random = makeRandom(31337);
        const points: Vector[] = [];
        for (let i = 0; i < 20; ++i) {
            points.push(v3(random() * 2 - 1, random() * 2 - 1, random() * 2 - 1));
        }
        const delaunay = makeTetrahedralization(points);
        const mesh = new Delaunay3Mesh(delaunay);

        let numFound = 0;
        for (let trial = 0; trial < 20; ++trial) {
            const p = v3(random() * 2 - 1, random() * 2 - 1, random() * 2 - 1);
            const t = mesh.getContainingTetrahedron(p);
            if (t === mesh.getInvalidIndex()) {
                continue;
            }
            ++numFound;
            expect(t).toBeGreaterThanOrEqual(0);
            expect(t).toBeLessThan(mesh.getNumTetrahedra());

            const bary = mesh.getBarycentrics(t, p) as number[];
            expect(bary).not.toBeNull();
            expect(bary[0] + bary[1] + bary[2] + bary[3]).toBeCloseTo(1, 12);
            for (const b of bary) {
                expect(b).toBeGreaterThanOrEqual(-1e-12);
            }

            const indices = mesh.getTetrahedronIndices(t) as number[];
            for (let j = 0; j < 3; ++j) {
                let value = 0;
                for (let i = 0; i < 4; ++i) {
                    value += bary[i] * points[indices[i]].values[j];
                }
                expect(value).toBeCloseTo(p.values[j], 10);
            }
        }
        expect(numFound).toBeGreaterThan(0);

        expect(mesh.getContainingTetrahedron(v3(50, 50, 50)))
            .toBe(mesh.getInvalidIndex());
    });
});
