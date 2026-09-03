import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPoint3Parallelepiped3 } from '../src/DistPoint3Parallelepiped3.js';
import type { DistPoint3Parallelepiped3Result }
    from '../src/DistPoint3Parallelepiped3.js';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox.js';
import { Matrix, mulMatrix } from '../src/Matrix.js';
import { inverse3x3 } from '../src/Matrix3x3.js';
import { Parallelepiped3 } from '../src/Parallelepiped3.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ppd(center: number[], a0: number[], a1: number[], a2: number[]):
    Parallelepiped3 {
    return Parallelepiped3.fromCenterAxis(v(...center),
        [v(...a0), v(...a1), v(...a2)]);
}

// The parallelepiped coordinates of a point: X = C + sum_i s[i] * A[i].
function coordinates(p: Vector, box: Parallelepiped3): Vector {
    const B = new Matrix(3, 3);
    B.setCol(0, box.axis[0]);
    B.setCol(1, box.axis[1]);
    B.setCol(2, box.axis[2]);
    return mulMatrix(inverse3x3(B).inverse, sub(p, box.center));
}

// A grid sampling of the solid parallelepiped gives an upper bound for the
// distance that converges to the true value as n grows.
function sampledDistance(p: Vector, box: Parallelepiped3, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i0 = 0; i0 <= n; ++i0) {
        const s0 = -1 + (2 * i0) / n;
        const t0 = add(box.center, mul(s0, box.axis[0]));
        for (let i1 = 0; i1 <= n; ++i1) {
            const s1 = -1 + (2 * i1) / n;
            const t1 = add(t0, mul(s1, box.axis[1]));
            for (let i2 = 0; i2 <= n; ++i2) {
                const s2 = -1 + (2 * i2) / n;
                const d = sub(p, add(t1, mul(s2, box.axis[2])));
                best = Math.min(best, dot(d, d));
            }
        }
    }
    return Math.sqrt(best);
}

function expectConsistent(p: Vector, box: Parallelepiped3,
    result: DistPoint3Parallelepiped3Result): void {
    expect(result.closest[0].equals(p)).toBe(true);

    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 9);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 9);

    // The closest point lies in the parallelepiped.
    const s = coordinates(result.closest[1], box);
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(s.values[i])).toBeLessThanOrEqual(1 + 1e-9);
    }
}

describe('DistPoint3Parallelepiped3', () => {
    const query = new DistPoint3Parallelepiped3();
    const unit = ppd([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]);

    it('reports zero distance for interior points', () => {
        for (const p of [v(0, 0, 0), v(0.5, -0.5, 0.25), v(1, 1, 1),
            v(-1, 0.3, 0.9)]) {
            const result = query.compute(p, unit);
            expect(result.distance).toBeCloseTo(0, 12);
            expectConsistent(p, unit, result);
        }
    });

    it('computes known exterior distances for the unit cube', () => {
        // Face region.
        let result = query.compute(v(3, 0, 0), unit);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);

        // Edge region.
        result = query.compute(v(3, 4, 0), unit);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9), 12);

        // Vertex region.
        result = query.compute(v(-3, -5, -7), unit);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 16 + 36), 12);
        expect(result.closest[1].values[0]).toBeCloseTo(-1, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(-1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(-1, 12);
    });

    it('agrees with the point-canonical-box query for orthogonal axes', () => {
        const extent = v(2, 3, 0.5);
        const cbox = CanonicalBox.fromExtent(extent);
        const box = ppd([0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 0.5]);
        const cbQuery = new DistPointCanonicalBox();

        let seed = 55555;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 300; ++trial) {
            const p = v(rand() * 12 - 6, rand() * 12 - 6, rand() * 12 - 6);
            const expected = cbQuery.compute(p, cbox);
            const result = query.compute(p, box);
            expect(result.distance).toBeCloseTo(expected.distance, 9);
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i])
                    .toBeCloseTo(expected.closest[1].values[i], 8);
            }
            expectConsistent(p, box, result);
        }
    });

    it('matches a grid sampling for a sheared parallelepiped', () => {
        let seed = 20260901;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const box = ppd([0.5, -0.25, 0.75], [1.5, 0.3, -0.2],
            [0.4, 1.2, 0.5], [-0.3, 0.2, 1.1]);
        for (let trial = 0; trial < 30; ++trial) {
            const p = v(rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4);
            const result = query.compute(p, box);
            const sampled = sampledDistance(p, box, 16);
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-9);
            expect(sampled - result.distance).toBeLessThan(0.25);
            expectConsistent(p, box, result);
        }
    });

    it('is invariant under rigid motion of the point and the solid', () => {
        const c = Math.cos(0.6), s = Math.sin(0.6);
        const rot = (p: Vector) => v(
            c * p.values[0] - s * p.values[2],
            p.values[1],
            s * p.values[0] + c * p.values[2]);
        const shift = v(1, -2, 3);

        const box = ppd([0, 0, 0], [1.2, 0.1, 0], [0.2, 1.1, 0.3],
            [0, 0.1, 0.9]);
        const moved = Parallelepiped3.fromCenterAxis(
            add(rot(box.center), shift),
            [rot(box.axis[0]), rot(box.axis[1]), rot(box.axis[2])]);

        let seed = 3141592;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 100; ++trial) {
            const p = v(rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4);
            const r0 = query.compute(p, box);
            const r1 = query.compute(add(rot(p), shift), moved);
            expect(r1.distance).toBeCloseTo(r0.distance, 9);
        }
    });
});
