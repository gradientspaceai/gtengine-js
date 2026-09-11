import { describe, it, expect } from 'vitest';
import { Polyhedron3 } from '../src/Polyhedron3.js';
import { Vector, add, sub, mul, length } from '../src/Vector.js';
import { cross, dotCross } from '../src/Vector3.js';
import { check, expectClose, fc, rotationFrame, vector }
    from './helpers/arbitraries.js';

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

describe('Polyhedron3 verification', () => {
    // A tetrahedron with a comfortably nonzero volume, as a vertex pool for
    // the module-level tetraIndices faces (12 indices, the minimum).
    const tetra = () => fc.tuple(vector(3, -4, 4), vector(3, -4, 4),
        vector(3, -4, 4), vector(3, -4, 4))
        .filter(vs => Math.abs(dotCross(sub(vs[1], vs[0]), sub(vs[2], vs[0]),
            sub(vs[3], vs[0]))) > 1)
        .map(vs => vs.map(v => v.clone()));

    it('the volume of a tetrahedron is |DotCross(e0,e1,e2)|/6', () => {
        check(tetra(), vs => {
            const p = new Polyhedron3(vs, 12, tetraIndices, true);
            expect(p.isValid()).toBe(true);
            const expected = Math.abs(dotCross(sub(vs[1], vs[0]),
                sub(vs[2], vs[0]), sub(vs[3], vs[0]))) / 6;
            expectClose(p.computeVolume(), expected, 1e-9, 1e-9);
        });
    });

    it('the surface area is the sum of the four triangle areas', () => {
        check(tetra(), vs => {
            const p = new Polyhedron3(vs, 12, tetraIndices, true);
            let expected = 0;
            for (let t = 0; t < 4; ++t) {
                const a = vs[tetraIndices[3 * t]];
                const b = vs[tetraIndices[3 * t + 1]];
                const c = vs[tetraIndices[3 * t + 2]];
                expected += 0.5 * length(cross(sub(b, a), sub(c, a)));
            }
            expectClose(p.computeSurfaceArea(), expected, 1e-9, 1e-9);
        });
    });

    it('the volume and area of a closed mesh are rigid-motion invariant',
        () => {
            // ComputeVolume sums DotCross(V0,V1,V2) over the faces, which is
            // the divergence-theorem formula: translation invariant only
            // because the mesh is closed.
            check(fc.tuple(tetra(), rotationFrame(3), vector(3, -20, 20)),
                ([vs, frame, shift]) => {
                    const p = new Polyhedron3(vs, 12, tetraIndices, true);
                    const moved = vs.map(v => add(shift,
                        add(mul(v.get(0), frame[0]),
                            add(mul(v.get(1), frame[1]),
                                mul(v.get(2), frame[2])))));
                    const q = new Polyhedron3(moved, 12, tetraIndices, true);
                    expectClose(q.computeVolume(), p.computeVolume(),
                        1e-9, 1e-9);
                    expectClose(q.computeSurfaceArea(),
                        p.computeSurfaceArea(), 1e-9, 1e-9);
                }, 100);
        });

    it('the vertex average uses the unique indices, not the whole pool', () => {
        check(fc.tuple(tetra(), vector(3, -100, 100)), ([vs, unused]) => {
            const pool = [...vs, unused];
            const p = new Polyhedron3(pool, 12, tetraIndices, true);
            const average = p.computeVertexAverage();
            for (let d = 0; d < 3; ++d) {
                expectClose(average.get(d), 0.25 * (vs[0].get(d)
                    + vs[1].get(d) + vs[2].get(d) + vs[3].get(d)),
                    1e-9, 1e-9);
            }
        });
    });

    it('the unique indices are ascending and duplicate free', () => {
        check(fc.array(fc.integer({ min: 0, max: 5 }),
            { minLength: 12, maxLength: 12 }).map(a => a.map(x => x)),
            indices => {
                const pool = [];
                for (let i = 0; i < 6; ++i) {
                    pool.push(Vector.fromArray([i, 0, 0]));
                }
                const p = new Polyhedron3(pool, 12, indices, true);
                const unique = p.getUniqueIndices();
                for (let i = 1; i < unique.length; ++i) {
                    expect(unique[i]).toBeGreaterThan(unique[i - 1]);
                }
                expect(new Set(indices).size).toBe(unique.length);
            });
    });

    it('the vertex pool is shared and the index array is copied', () => {
        check(tetra(), vs => {
            const indices = [...tetraIndices];
            const p = new Polyhedron3(vs, 12, indices, true);
            indices[0] = 3;
            expect(p.getIndices()[0]).toBe(0);
            // Upstream stores the shared_ptr, so pool edits are visible.
            const before = p.computeVertexAverage().get(0);
            vs[0].set(0, vs[0].get(0) + 4);
            expectClose(p.computeVertexAverage().get(0), before + 1,
                1e-9, 1e-9);
        });
    });

    it('rejects index counts below 12 or not a multiple of 3', () => {
        check(fc.integer({ min: 0, max: 20 }), numIndices => {
            const pool = [];
            for (let i = 0; i < 8; ++i) {
                pool.push(Vector.fromArray([i, 0, 0]));
            }
            const indices = new Array<number>(24).fill(0);
            const p = new Polyhedron3(pool, numIndices, indices, true);
            expect(p.isValid())
                .toBe(numIndices >= 12 && numIndices % 3 === 0);
            if (!p.isValid()) {
                expect(p.counterClockwise()).toBe(false);
                expect(p.getIndices().length).toBe(0);
                expect(p.getUniqueIndices().length).toBe(0);
                expect(p.computeVolume()).toBe(0);
                expect(p.computeSurfaceArea()).toBe(0);
                expect(p.computeVertexAverage().values).toEqual([0, 0, 0]);
            }
        });
    });
});
