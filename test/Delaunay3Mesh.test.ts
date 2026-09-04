import { describe, it, expect } from 'vitest';
import { Delaunay3 } from '../src/Delaunay3.js';
import { Delaunay3Mesh } from '../src/Delaunay3Mesh.js';
import { Vector } from '../src/Vector.js';
import { fc, check, latticeVector, expectClose } from './helpers/arbitraries.js';
import { orient3 } from './helpers/exact.js';

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

// ---------------------------------------------------------------------------
// Verification (V11): property-based cross-checks. All point sets are integer
// lattices, so the containment oracle is exact bigint orientation arithmetic.
// ---------------------------------------------------------------------------

const bigXYZ = (p: Vector): [bigint, bigint, bigint] =>
    [BigInt(p.values[0]), BigInt(p.values[1]), BigInt(p.values[2])];

/**
 * The indices of every tetrahedron of the tetrahedralization that contains P,
 * determined exactly. The test is independent of the vertex ordering
 * convention: P is in the closed tetrahedron when replacing any one vertex by
 * P does not flip the orientation.
 */
function containingTetrahedraExact(mesh: Delaunay3Mesh, p: Vector): number[] {
    const vertices = mesh.getVertices();
    const q = bigXYZ(p);
    const result: number[] = [];
    for (let t = 0; t < mesh.getNumTetrahedra(); ++t) {
        const idx = mesh.getTetrahedronIndices(t) as
            [number, number, number, number];
        const v = idx.map(i => bigXYZ(vertices[i]));
        const det = orient3(v[0], v[1], v[2], v[3]);
        if (det === 0) { continue; }
        let inside = true;
        for (let i = 0; i < 4; ++i) {
            const w = v.slice();
            w[i] = q;
            const s = orient3(w[0], w[1], w[2], w[3]);
            if (s !== 0 && s !== det) {
                inside = false;
                break;
            }
        }
        if (inside) { result.push(t); }
    }
    return result;
}

const latticeCloud3 = (count: number, range: number): fc.Arbitrary<Vector[]> =>
    fc.array(latticeVector(3, -range, range),
        { minLength: count, maxLength: count });

function buildMesh3(points: readonly Vector[]): Delaunay3Mesh | null {
    const delaunay = new Delaunay3();
    if (!delaunay.compute(points) || delaunay.getDimension() !== 3) {
        return null;
    }
    return new Delaunay3Mesh(delaunay);
}

describe('Delaunay3Mesh verification', () => {
    it('locates the containing tetrahedron exactly', () => {
        check(fc.tuple(latticeCloud3(14, 5), latticeVector(3, -2, 2)),
            ([points, query]) => {
                const mesh = buildMesh3(points);
                if (mesh === null) { return true; }
                const expected = containingTetrahedraExact(mesh, query);
                const t = mesh.getContainingTetrahedron(query);
                if (expected.length === 0) {
                    expect(t).toBe(mesh.getInvalidIndex());
                }
                else {
                    // On a shared face, edge or vertex several tetrahedra
                    // contain the point; the search may return any of them.
                    expect(expected).toContain(t);
                }
                return true;
            }, 120);
    }, 30000);

    it('computes barycentrics that reproduce the query point', () => {
        check(fc.tuple(latticeCloud3(14, 5), latticeVector(3, -2, 2)),
            ([points, query]) => {
                const mesh = buildMesh3(points);
                if (mesh === null) { return true; }
                const t = mesh.getContainingTetrahedron(query);
                if (t === mesh.getInvalidIndex()) { return true; }

                const bary = mesh.getBarycentrics(t, query);
                expect(bary).not.toBeNull();
                const b = bary as [number, number, number, number];
                const vertices = mesh.getTetrahedronVertices(t) as
                    [Vector, Vector, Vector, Vector];

                expectClose(b[0] + b[1] + b[2] + b[3], 1, 1e-12, 1e-12);
                for (const bi of b) {
                    expect(bi).toBeGreaterThanOrEqual(-1e-12);
                    expect(bi).toBeLessThanOrEqual(1 + 1e-12);
                }

                // The coordinates are exact rationals rounded once to double,
                // so the reconstruction is accurate to a few ulps of the
                // coordinate range.
                for (let j = 0; j < 3; ++j) {
                    let value = 0;
                    for (let i = 0; i < 4; ++i) {
                        value += b[i] * vertices[i].values[j];
                    }
                    expectClose(value, query.values[j], 1e-11, 1e-11);
                }
                return true;
            }, 120);
    }, 30000);

    it('interpolates a linear function exactly', () => {
        check(fc.tuple(latticeCloud3(14, 5), latticeVector(3, -2, 2),
            latticeVector(4, -4, 4)), ([points, query, coeff]) => {
            const mesh = buildMesh3(points);
            if (mesh === null) { return true; }
            const t = mesh.getContainingTetrahedron(query);
            if (t === mesh.getInvalidIndex()) { return true; }
            const bary = mesh.getBarycentrics(t, query);
            if (bary === null) { return true; }

            const f = (v: Vector): number =>
                coeff.values[0] * v.values[0] + coeff.values[1] * v.values[1]
                + coeff.values[2] * v.values[2] + coeff.values[3];
            const vertices = mesh.getTetrahedronVertices(t) as
                [Vector, Vector, Vector, Vector];
            let interpolated = 0;
            for (let i = 0; i < 4; ++i) {
                interpolated += bary[i] * f(vertices[i]);
            }
            expectClose(interpolated, f(query), 1e-10, 1e-10);
            return true;
        }, 120);
    }, 30000);

    it('keeps the accessors consistent with the tetrahedralization', () => {
        check(latticeCloud3(12, 5), points => {
            const mesh = buildMesh3(points);
            if (mesh === null) { return true; }
            const vertices = mesh.getVertices();
            const indices = mesh.getIndices();
            const adjacencies = mesh.getAdjacencies();
            expect(indices.length).toBe(4 * mesh.getNumTetrahedra());
            expect(adjacencies.length).toBe(indices.length);

            for (let t = 0; t < mesh.getNumTetrahedra(); ++t) {
                const idx = mesh.getTetrahedronIndices(t) as
                    [number, number, number, number];
                const adj = mesh.getTetrahedronAdjacencies(t) as
                    [number, number, number, number];
                const tv = mesh.getTetrahedronVertices(t) as
                    [Vector, Vector, Vector, Vector];
                for (let j = 0; j < 4; ++j) {
                    expect(idx[j]).toBe(indices[4 * t + j]);
                    expect(adj[j]).toBe(adjacencies[4 * t + j]);
                    expect(tv[j].values).toEqual(vertices[idx[j]].values);
                    // The returned vertices are copies, not aliases of the
                    // tetrahedralization storage (upstream returns by value).
                    expect(tv[j]).not.toBe(vertices[idx[j]]);
                    if (adj[j] !== mesh.getInvalidIndex()) {
                        const back = mesh.getTetrahedronAdjacencies(adj[j]) as
                            [number, number, number, number];
                        expect(back).toContain(t);
                    }
                }
            }
            expect(mesh.getTetrahedronIndices(-1)).toBeNull();
            expect(mesh.getTetrahedronIndices(mesh.getNumTetrahedra())).toBeNull();
            expect(mesh.getBarycentrics(mesh.getNumTetrahedra(),
                points[0])).toBeNull();
            return true;
        }, 80);
    }, 30000);
});
