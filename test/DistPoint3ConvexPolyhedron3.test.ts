import { describe, expect, it } from 'vitest';
import { ConvexPolyhedron3 } from '../src/ConvexPolyhedron3';
import { DistPoint3ConvexPolyhedron3 }
    from '../src/DistPoint3ConvexPolyhedron3';
import { Vector, length, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// The unit cube [0,1]^3 with faces wound counterclockwise when viewed from
// outside, so the generated plane normals point outward.
function cube(): ConvexPolyhedron3 {
    const vertices = [
        v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(0, 1, 0),
        v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)
    ];
    const indices = [
        0, 3, 2, 0, 2, 1,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5
    ];
    return new ConvexPolyhedron3(vertices, indices, true, true);
}

function tetrahedron(): ConvexPolyhedron3 {
    const vertices = [v(0, 0, 0), v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
    const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
    return new ConvexPolyhedron3(vertices, indices, true, true);
}

// The exact distance from a point to the solid unit cube [0,1]^3.
function cubeDistance(p: Vector): number {
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const d = Math.max(0, -p.values[i], p.values[i] - 1);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

describe('DistPoint3ConvexPolyhedron3', () => {
    it('computes the distance to a cube face', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(2, 0.5, 0.5), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(1, 8);
        expect(result.closest[0].values).toEqual([2, 0.5, 0.5]);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 8);
        expect(result.closest[1].values[1]).toBeCloseTo(0.5, 8);
        expect(result.closest[1].values[2]).toBeCloseTo(0.5, 8);
    });

    it('computes the distance to a cube corner', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(2, 2, 2), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(Math.sqrt(3), 7);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1, 7);
        }
    });

    it('computes the distance to a cube edge', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(2, 2, 0.5), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(Math.SQRT2, 7);
    });

    it('returns zero distance for a point inside', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(0.5, 0.25, 0.75), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(0, 8);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i])
                .toBeCloseTo(result.closest[0].values[i], 8);
        }
    });

    it('returns zero distance for a point on the boundary', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(1, 0.5, 0.5), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(0, 8);
    });

    it('handles a tetrahedron', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const t = tetrahedron();
        // The point (1,1,1) is outside the plane x+y+z=1; the closest point
        // is the face centroid (1/3,1/3,1/3) and the distance is
        // (3-1)/sqrt(3).
        const result = query.compute(v(1, 1, 1), t);
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(2 / Math.sqrt(3), 7);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1 / 3, 6);
        }

        const inside = query.compute(v(0.1, 0.1, 0.1), t);
        expect(inside.distance).toBeCloseTo(0, 8);
    });

    it('reports failure for a polyhedron with no planes', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const empty = new ConvexPolyhedron3();
        const result = query.compute(v(1, 2, 3), empty);
        expect(result.queryIsSuccessful).toBe(false);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expect(result.numLCPIterations).toBe(0);
        expect(result.closest[0].values).toEqual([0, 0, 0]);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('reuses a preallocated LCP solver', () => {
        // The cube has 12 triangle faces.
        const query = new DistPoint3ConvexPolyhedron3(12);
        const c = cube();
        const r0 = query.compute(v(3, 0.5, 0.5), c);
        const r1 = query.compute(v(-2, 0.5, 0.5), c);
        expect(r0.distance).toBeCloseTo(2, 7);
        expect(r1.distance).toBeCloseTo(2, 7);
        expect(r0.numLCPIterations).toBeGreaterThan(0);
    });

    it('honors setMaxLCPIterations', () => {
        const query = new DistPoint3ConvexPolyhedron3(12);
        query.setMaxLCPIterations(500);
        const result = query.compute(v(2, 2, 2), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(Math.sqrt(3), 7);
    });

    it('agrees with the exact cube distance on random points', () => {
        const query = new DistPoint3ConvexPolyhedron3(12);
        const c = cube();
        let seed = 456789123;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const p = v(4 * rand() - 1.5, 4 * rand() - 1.5, 4 * rand() - 1.5);
            const result = query.compute(p, c);
            expect(result.queryIsSuccessful).toBe(true);
            expect(result.distance).toBeCloseTo(cubeDistance(p), 6);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 6);
            // The closest point is in the cube.
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i])
                    .toBeGreaterThanOrEqual(-1e-6);
                expect(result.closest[1].values[i])
                    .toBeLessThanOrEqual(1 + 1e-6);
            }
        }
    });
});
