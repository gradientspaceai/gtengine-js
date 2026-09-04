import { describe, it, expect } from 'vitest';
import { ExtremalQuery3 } from '../src/ExtremalQuery3.js';
import type { ExtremalQuery3Result } from '../src/ExtremalQuery3.js';
import { Polyhedron3 } from '../src/Polyhedron3.js';
import { Vector, dot } from '../src/Vector.js';
import { check, fc, rotationFrame, scaled, unitVector } from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A minimal concrete subclass: the extreme vertices are found by an
// exhaustive search over the polyhedron vertices. Later batches port the
// BSP and projection accelerations (ExtremalQuery3BSP, ExtremalQuery3PRJ).
class BruteForceExtremalQuery3 extends ExtremalQuery3 {
    constructor(polytope: Polyhedron3) {
        super(polytope);
    }

    override getExtremeVertices(direction: Vector): ExtremalQuery3Result {
        const vertices = this.mPolytope.getVertices();
        const unique = this.mPolytope.getUniqueIndices();
        let positiveDirection = unique[0];
        let negativeDirection = unique[0];
        let maxDot = dot(direction, vertices[unique[0]]);
        let minDot = maxDot;
        for (const i of unique) {
            const d = dot(direction, vertices[i]);
            if (d > maxDot) {
                maxDot = d;
                positiveDirection = i;
            }
            if (d < minDot) {
                minDot = d;
                negativeDirection = i;
            }
        }
        return { positiveDirection, negativeDirection };
    }
}

// The unit cube [0,1]^3 with counterclockwise-when-viewed-from-outside faces.
const cubeVertices: Vector[] = [
    v3(0, 0, 0), v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0),
    v3(0, 0, 1), v3(1, 0, 1), v3(1, 1, 1), v3(0, 1, 1)
];

const cubeIndices: number[] = [
    0, 3, 2, 0, 2, 1,
    4, 5, 6, 4, 6, 7,
    0, 1, 5, 0, 5, 4,
    3, 7, 6, 3, 6, 2,
    0, 4, 7, 0, 7, 3,
    1, 2, 6, 1, 6, 5
];

const tetraVertices: Vector[] = [
    v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)
];

const tetraIndices: number[] = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];

function makeCube(): Polyhedron3 {
    return new Polyhedron3(cubeVertices, cubeIndices.length, cubeIndices, true);
}

function makeTetra(): Polyhedron3 {
    return new Polyhedron3(tetraVertices, tetraIndices.length, tetraIndices, true);
}

describe('ExtremalQuery3', () => {
    it('stores the polytope by reference', () => {
        const cube = makeCube();
        const query = new BruteForceExtremalQuery3(cube);
        expect(query.getPolytope()).toBe(cube);
    });

    it('computes one unit-length face normal per triangle', () => {
        const query = new BruteForceExtremalQuery3(makeCube());
        const normals = query.getFaceNormals();
        expect(normals.length).toBe(12);
        for (const n of normals) {
            expect(n.size).toBe(3);
            const len = Math.sqrt(dot(n, n));
            expect(len).toBeCloseTo(1, 12);
        }
    });

    it('computes outer-pointing cube face normals', () => {
        const query = new BruteForceExtremalQuery3(makeCube());
        const normals = query.getFaceNormals();
        const expected = [
            [0, 0, -1], [0, 0, -1],
            [0, 0, 1], [0, 0, 1],
            [0, -1, 0], [0, -1, 0],
            [0, 1, 0], [0, 1, 0],
            [-1, 0, 0], [-1, 0, 0],
            [1, 0, 0], [1, 0, 0]
        ];
        for (let t = 0; t < 12; ++t) {
            for (let d = 0; d < 3; ++d) {
                expect(normals[t].values[d]).toBeCloseTo(expected[t][d], 12);
            }
        }
    });

    it('has face normals that point away from the polyhedron centroid', () => {
        const query = new BruteForceExtremalQuery3(makeTetra());
        const normals = query.getFaceNormals();
        expect(normals.length).toBe(4);
        const centroid = v3(0.25, 0.25, 0.25);
        const indices = query.getPolytope().getIndices();
        const vertices = query.getPolytope().getVertices();
        for (let t = 0; t < 4; ++t) {
            const V0 = vertices[indices[3 * t]];
            // Dot(N, centroid - V0) < 0 for an outer-pointing normal.
            const diff = v3(
                centroid.values[0] - V0.values[0],
                centroid.values[1] - V0.values[1],
                centroid.values[2] - V0.values[2]);
            expect(dot(normals[t], diff)).toBeLessThan(0);
        }
    });

    it('normalizes the tetrahedron slanted-face normal', () => {
        const query = new BruteForceExtremalQuery3(makeTetra());
        const n = query.getFaceNormals()[3];
        const s = 1 / Math.sqrt(3);
        expect(n.values[0]).toBeCloseTo(s, 12);
        expect(n.values[1]).toBeCloseTo(s, 12);
        expect(n.values[2]).toBeCloseTo(s, 12);
    });

    it('returns the same face-normal array on repeated access', () => {
        const query = new BruteForceExtremalQuery3(makeCube());
        expect(query.getFaceNormals()).toBe(query.getFaceNormals());
    });

    it('supports the abstract extreme-vertex query through a subclass', () => {
        const query = new BruteForceExtremalQuery3(makeCube());
        const vertices = query.getPolytope().getVertices();

        const along = query.getExtremeVertices(v3(1, 0, 0));
        expect(vertices[along.positiveDirection].values[0]).toBe(1);
        expect(vertices[along.negativeDirection].values[0]).toBe(0);

        const diagonal = query.getExtremeVertices(v3(1, 1, 1));
        expect(vertices[diagonal.positiveDirection].values).toEqual([1, 1, 1]);
        expect(vertices[diagonal.negativeDirection].values).toEqual([0, 0, 0]);
    });

    it('is consistent under negation of the direction', () => {
        const query = new BruteForceExtremalQuery3(makeCube());
        const d = v3(0.3, -0.7, 0.2);
        const r0 = query.getExtremeVertices(d);
        const r1 = query.getExtremeVertices(v3(-0.3, 0.7, -0.2));
        expect(r1.positiveDirection).toBe(r0.negativeDirection);
        expect(r1.negativeDirection).toBe(r0.positiveDirection);
    });

    it('throws when the polytope is invalid (null vertex pool)', () => {
        const invalid = new Polyhedron3(cubeVertices, 6, [0, 1, 2, 0, 2, 3], true);
        expect(invalid.isValid()).toBe(false);
        expect(() => new BruteForceExtremalQuery3(invalid)).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). The base class only computes
// one unit-length outer-pointing normal per triangle, so the properties check
// that construction against an independent normalization of the cross product
// over random affine images of two convex solids.
// ---------------------------------------------------------------------------

const affinePolytope = fc.record({
    which: fc.integer({ min: 0, max: 1 }),
    frame: rotationFrame(3),
    scale: fc.tuple(scaled(0.5, 3, 16), scaled(0.5, 3, 16), scaled(0.5, 3, 16)),
    translate: fc.tuple(scaled(-5, 5, 32), scaled(-5, 5, 32), scaled(-5, 5, 32))
}).map(({ which, frame, scale, translate }) => {
    const base = which === 0
        ? { vertices: cubeVertices, indices: cubeIndices }
        : { vertices: tetraVertices, indices: tetraIndices };
    const map = (p: Vector): Vector => {
        const s = [scale[0] * p.get(0), scale[1] * p.get(1), scale[2] * p.get(2)];
        return v3(
            frame[0].get(0) * s[0] + frame[1].get(0) * s[1] + frame[2].get(0) * s[2] + translate[0],
            frame[0].get(1) * s[0] + frame[1].get(1) * s[1] + frame[2].get(1) * s[2] + translate[1],
            frame[0].get(2) * s[0] + frame[1].get(2) * s[1] + frame[2].get(2) * s[2] + translate[2]);
    };
    const vertices = base.vertices.map(map);
    const indices = base.indices.slice();
    return new Polyhedron3(vertices, indices.length, indices, true);
});

describe('ExtremalQuery3 verification', () => {
    it('face normals are the unit cross products of the face edges', () => {
        check(affinePolytope, (polytope) => {
            const query = new BruteForceExtremalQuery3(polytope);
            const normals = query.getFaceNormals();
            const vertices = polytope.getVertices();
            const indices = polytope.getIndices();
            expect(normals.length).toBe(indices.length / 3);
            for (let t = 0; t < normals.length; ++t) {
                const V0 = vertices[indices[3 * t + 0]];
                const V1 = vertices[indices[3 * t + 1]];
                const V2 = vertices[indices[3 * t + 2]];
                const e1 = [V1.get(0) - V0.get(0), V1.get(1) - V0.get(1), V1.get(2) - V0.get(2)];
                const e2 = [V2.get(0) - V0.get(0), V2.get(1) - V0.get(1), V2.get(2) - V0.get(2)];
                const n = [e1[1] * e2[2] - e1[2] * e2[1],
                    e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0]];
                const len = Math.hypot(n[0], n[1], n[2]);
                for (let i = 0; i < 3; ++i) {
                    // UnitCross normalizes the cross product; the difference
                    // from an independent normalization is a few ulps.
                    expect(Math.abs(normals[t].get(i) - n[i] / len))
                        .toBeLessThanOrEqual(1e-12);
                }
                expect(Math.hypot(normals[t].get(0), normals[t].get(1),
                    normals[t].get(2))).toBeCloseTo(1, 12);
            }
        });
    });

    it('face normals point away from the interior', () => {
        check(affinePolytope, (polytope) => {
            const query = new BruteForceExtremalQuery3(polytope);
            const normals = query.getFaceNormals();
            const vertices = polytope.getVertices();
            const indices = polytope.getIndices();
            const centroid = polytope.computeVertexAverage();
            for (let t = 0; t < normals.length; ++t) {
                const V0 = vertices[indices[3 * t + 0]];
                // Dot(N, centroid - V0) < 0 for an outer-pointing normal.
                const inward = dot(normals[t], v3(centroid.get(0) - V0.get(0),
                    centroid.get(1) - V0.get(1), centroid.get(2) - V0.get(2)));
                expect(inward).toBeLessThan(0);
                // Every vertex is on the inner side of every face plane.
                for (const V of vertices) {
                    const d = dot(normals[t], v3(V.get(0) - V0.get(0),
                        V.get(1) - V0.get(1), V.get(2) - V0.get(2)));
                    expect(d).toBeLessThanOrEqual(1e-12 * Math.max(1, Math.abs(inward)));
                }
            }
        });
    });

    it('holds the polytope by reference and caches the normal array', () => {
        check(affinePolytope, (polytope) => {
            const query = new BruteForceExtremalQuery3(polytope);
            expect(query.getPolytope()).toBe(polytope);
            expect(query.getFaceNormals()).toBe(query.getFaceNormals());
        });
    });

    it('the extreme-vertex query is antisymmetric and supporting', () => {
        check(fc.tuple(affinePolytope, unitVector(3)), ([polytope, direction]) => {
            const query = new BruteForceExtremalQuery3(polytope);
            const result = query.getExtremeVertices(direction);
            const negated = query.getExtremeVertices(
                v3(-direction.get(0), -direction.get(1), -direction.get(2)));
            const vertices = polytope.getVertices();
            // Antisymmetry holds for the projected values (index ties may pick
            // different vertices of the same supporting face).
            expect(dot(direction, vertices[result.positiveDirection]))
                .toBeCloseTo(dot(direction, vertices[negated.negativeDirection]), 12);
            expect(dot(direction, vertices[result.negativeDirection]))
                .toBeCloseTo(dot(direction, vertices[negated.positiveDirection]), 12);
            // No vertex projects beyond the reported extremes.
            const hi = dot(direction, vertices[result.positiveDirection]);
            const lo = dot(direction, vertices[result.negativeDirection]);
            for (const i of polytope.getUniqueIndices()) {
                const d = dot(direction, vertices[i]);
                expect(d).toBeLessThanOrEqual(hi);
                expect(d).toBeGreaterThanOrEqual(lo);
            }
        });
    });
});
