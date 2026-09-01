import { describe, it, expect } from 'vitest';
import { ExtremalQuery3 } from '../src/ExtremalQuery3';
import type { ExtremalQuery3Result } from '../src/ExtremalQuery3';
import { Polyhedron3 } from '../src/Polyhedron3';
import { Vector, dot } from '../src/Vector';

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
