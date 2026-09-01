import { describe, it, expect } from 'vitest';
import { Polyhedron3 } from '../src/Polyhedron3';
import { Vector } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The unit cube [0,1]^3. Vertex i has bit pattern chosen so that the faces
// below are outward-pointing and counterclockwise.
const cubeVertices: Vector[] = [
    v3(0, 0, 0), v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0),
    v3(0, 0, 1), v3(1, 0, 1), v3(1, 1, 1), v3(0, 1, 1)
];

const cubeIndices: number[] = [
    0, 3, 2, 0, 2, 1, // z = 0
    4, 5, 6, 4, 6, 7, // z = 1
    0, 1, 5, 0, 5, 4, // y = 0
    3, 7, 6, 3, 6, 2, // y = 1
    0, 4, 7, 0, 7, 3, // x = 0
    1, 2, 6, 1, 6, 5  // x = 1
];

// The tetrahedron with vertices at the origin and the three unit points.
const tetraVertices: Vector[] = [
    v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)
];

const tetraIndices: number[] = [
    0, 2, 1, // z = 0
    0, 3, 2, // x = 0
    0, 1, 3, // y = 0
    1, 2, 3  // the slanted face
];

describe('Polyhedron3 construction', () => {
    it('accepts a valid cube and copies the indices', () => {
        const indices = cubeIndices.slice();
        const poly = new Polyhedron3(cubeVertices, indices.length, indices,
            true);
        expect(poly.isValid()).toBe(true);
        expect(poly.counterClockwise()).toBe(true);
        expect(poly.getIndices()).toEqual(cubeIndices);
        indices[0] = 99;
        expect(poly.getIndices()[0]).toBe(0);
    });

    it('shares (does not copy) the vertex pool, as upstream does', () => {
        const poly = new Polyhedron3(cubeVertices, cubeIndices.length,
            cubeIndices, true);
        expect(poly.getVertexPool()).toBe(cubeVertices);
        expect(poly.getVertices()).toBe(cubeVertices);
    });

    it('records the unique indices in ascending order', () => {
        // A pool with two unused points; only 0..7 are referenced.
        const pool = cubeVertices.concat([v3(9, 9, 9), v3(8, 8, 8)]);
        const poly = new Polyhedron3(pool, cubeIndices.length, cubeIndices,
            true);
        expect(poly.getUniqueIndices()).toEqual([0, 1, 2, 3, 4, 5, 6, 7]);
    });

    it('rejects invalid inputs and resets to the failed state', () => {
        const check = (poly: Polyhedron3): void => {
            expect(poly.isValid()).toBe(false);
            expect(poly.counterClockwise()).toBe(false);
            expect(poly.getIndices()).toEqual([]);
            expect(poly.getUniqueIndices()).toEqual([]);
            expect(poly.getVertexPool()).toBeNull();
            expect(() => poly.getVertices()).toThrow();
        };

        // Fewer than 12 indices (fewer than four triangles).
        check(new Polyhedron3(tetraVertices, 9, [0, 1, 2, 0, 2, 3, 0, 3, 1],
            true));
        // Not a multiple of 3.
        check(new Polyhedron3(cubeVertices, 13, cubeIndices, true));
        // Null vertex pool.
        check(new Polyhedron3(null, cubeIndices.length, cubeIndices, true));
        // Null indices.
        check(new Polyhedron3(cubeVertices, 12, null, true));
        // Port deviation: an indices array shorter than numIndices.
        check(new Polyhedron3(cubeVertices, 36, [0, 1, 2, 3], true));
    });

    it('an invalid polyhedron reports zero geometric quantities', () => {
        const poly = new Polyhedron3(null, 36, cubeIndices, true);
        expect(poly.computeSurfaceArea()).toBe(0);
        expect(poly.computeVolume()).toBe(0);
        expect(poly.computeVertexAverage().values).toEqual([0, 0, 0]);
    });
});

describe('Polyhedron3 geometric queries: the unit cube', () => {
    const poly = new Polyhedron3(cubeVertices, cubeIndices.length,
        cubeIndices, true);

    it('the vertex average is the cube center', () => {
        expect(poly.computeVertexAverage().values).toEqual([0.5, 0.5, 0.5]);
    });

    it('the surface area is 6', () => {
        expect(poly.computeSurfaceArea()).toBeCloseTo(6, 12);
    });

    it('the volume is 1', () => {
        expect(poly.computeVolume()).toBeCloseTo(1, 12);
    });

    it('the queries ignore unused points in the vertex pool', () => {
        const pool = cubeVertices.concat([v3(100, 100, 100)]);
        const other = new Polyhedron3(pool, cubeIndices.length, cubeIndices,
            true);
        expect(other.computeVertexAverage().values).toEqual([0.5, 0.5, 0.5]);
        expect(other.computeSurfaceArea()).toBeCloseTo(6, 12);
        expect(other.computeVolume()).toBeCloseTo(1, 12);
    });

    it('the volume is independent of the face orientation', () => {
        // Reverse each triangle to flip the orientation; upstream takes the
        // absolute value of the signed volume.
        const reversed: number[] = [];
        for (let i = 0; i < cubeIndices.length; i += 3) {
            reversed.push(cubeIndices[i + 2], cubeIndices[i + 1],
                cubeIndices[i]);
        }
        const other = new Polyhedron3(cubeVertices, reversed.length, reversed,
            false);
        expect(other.counterClockwise()).toBe(false);
        expect(other.computeVolume()).toBeCloseTo(1, 12);
        expect(other.computeSurfaceArea()).toBeCloseTo(6, 12);
    });

    it('the volume is invariant under translation of the whole cube', () => {
        // The divergence-theorem sum uses positions relative to the origin,
        // but a closed surface makes the result translation invariant.
        const shifted = cubeVertices.map(v => v3(v.get(0) + 10,
            v.get(1) - 5, v.get(2) + 3));
        const other = new Polyhedron3(shifted, cubeIndices.length, cubeIndices,
            true);
        expect(other.computeVolume()).toBeCloseTo(1, 10);
        expect(other.computeSurfaceArea()).toBeCloseTo(6, 12);
        expect(other.computeVertexAverage().values)
            .toEqual([10.5, -4.5, 3.5]);
    });

    it('scaling the cube by s scales area by s^2 and volume by s^3', () => {
        const s = 3;
        const scaled = cubeVertices.map(v => v3(s * v.get(0), s * v.get(1),
            s * v.get(2)));
        const other = new Polyhedron3(scaled, cubeIndices.length, cubeIndices,
            true);
        expect(other.computeSurfaceArea()).toBeCloseTo(6 * s * s, 12);
        expect(other.computeVolume()).toBeCloseTo(s * s * s, 12);
    });
});

describe('Polyhedron3 geometric queries: the unit tetrahedron', () => {
    const poly = new Polyhedron3(tetraVertices, tetraIndices.length,
        tetraIndices, true);

    it('is valid with exactly four triangles', () => {
        expect(poly.isValid()).toBe(true);
        expect(poly.getIndices().length).toBe(12);
        expect(poly.getUniqueIndices()).toEqual([0, 1, 2, 3]);
    });

    it('the vertex average is (1/4, 1/4, 1/4)', () => {
        const average = poly.computeVertexAverage();
        expect(average.get(0)).toBeCloseTo(0.25, 12);
        expect(average.get(1)).toBeCloseTo(0.25, 12);
        expect(average.get(2)).toBeCloseTo(0.25, 12);
    });

    it('the surface area is 3/2 + sqrt(3)/2', () => {
        expect(poly.computeSurfaceArea())
            .toBeCloseTo(1.5 + Math.sqrt(3) / 2, 12);
    });

    it('the volume is 1/6', () => {
        expect(poly.computeVolume()).toBeCloseTo(1 / 6, 12);
    });
});
