import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistPlane3AlignedBox3 } from '../src/DistPlane3AlignedBox3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { Vector, dot, normalize, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = v(...normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, v(...origin));
}

// The distance from a plane to a solid box is zero when the box straddles
// the plane; otherwise the closest box point is a vertex, so the distance is
// the smallest |Dot(N,V) - c| over the vertices.
function analyticDistance(p: Hyperplane, b: AlignedBox): number {
    const vertices = b.getVertices();
    let minSigned = Number.MAX_VALUE;
    let maxSigned = -Number.MAX_VALUE;
    for (const vertex of vertices) {
        const s = dot(p.normal, vertex) - p.constant;
        minSigned = Math.min(minSigned, s);
        maxSigned = Math.max(maxSigned, s);
    }
    if (minSigned <= 0 && maxSigned >= 0) {
        return 0;
    }
    return minSigned > 0 ? minSigned : -maxSigned;
}

function expectConsistent(p: Hyperplane, b: AlignedBox,
    result: { distance: number, sqrDistance: number, closest: [Vector, Vector] }):
    void {
    // The reported distance is the length of the segment joining the closest
    // points.
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 10);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 10);

    // closest[0] is on the plane.
    expect(dot(p.normal, result.closest[0])).toBeCloseTo(p.constant, 10);

    // closest[1] is in the box.
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[1].values[i]).toBeGreaterThanOrEqual(
            b.min.values[i] - 1e-10);
        expect(result.closest[1].values[i]).toBeLessThanOrEqual(
            b.max.values[i] + 1e-10);
    }
}

describe('DistPlane3AlignedBox3', () => {
    const query = new DistPlane3AlignedBox3();

    it('computes the distance for an axis-aligned separating plane', () => {
        const b = box([-1, -1, -1], [1, 1, 1]);
        const p = plane([0, 0, 1], [0, 0, 5]);
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.sqrDistance).toBeCloseTo(16, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
        expectConsistent(p, b, result);
    });

    it('reports zero distance when the plane cuts the box', () => {
        const b = box([-1, -2, -3], [4, 5, 6]);
        const p = plane([1, 1, 1], [0, 0, 0]);
        const result = query.compute(p, b);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expectConsistent(p, b, result);
    });

    it('handles a diagonal plane touching a box corner', () => {
        // The plane Dot((1,1,1)/sqrt(3), X) = sqrt(3) passes through the
        // corner (1,1,1) of the unit box.
        const b = box([-1, -1, -1], [1, 1, 1]);
        const p = plane([1, 1, 1], [1, 1, 1]);
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(p, b, result);
    });

    it('is translation covariant', () => {
        const b0 = box([-1, -2, -3], [1, 2, 3]);
        const b1 = box([9, 8, 7], [11, 12, 13]);
        const p0 = plane([1, 2, 3], [0, 0, 8]);
        const p1 = plane([1, 2, 3], [10, 10, 18]);
        const r0 = query.compute(p0, b0);
        const r1 = query.compute(p1, b1);
        expect(r1.distance).toBeCloseTo(r0.distance, 10);
        for (let i = 0; i < 3; ++i) {
            expect(r1.closest[0].values[i] - r0.closest[0].values[i])
                .toBeCloseTo(10, 10);
            expect(r1.closest[1].values[i] - r0.closest[1].values[i])
                .toBeCloseTo(10, 10);
        }
    });

    it('matches the analytic vertex-based distance on random inputs', () => {
        let seed = 987654321;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 400; ++trial) {
            const min = [rand() * 4 - 2, rand() * 4 - 2, rand() * 4 - 2];
            const max = [min[0] + rand() * 3 + 0.1, min[1] + rand() * 3 + 0.1,
                min[2] + rand() * 3 + 0.1];
            const b = box(min, max);
            const p = plane(
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                [rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4]);
            const result = query.compute(p, b);
            expect(result.distance).toBeCloseTo(analyticDistance(p, b), 9);
            expectConsistent(p, b, result);
        }
    });
});
