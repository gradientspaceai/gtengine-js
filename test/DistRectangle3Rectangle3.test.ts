import { describe, expect, it } from 'vitest';
import { DistRectangle3Rectangle3 } from '../src/DistRectangle3Rectangle3.js';
import type { DistRectangle3Rectangle3Result }
    from '../src/DistRectangle3Rectangle3.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// An orthonormal frame parameterized by two angles.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        v(ca, sa, 0),
        v(-sa * cb, ca * cb, sb),
        v(sa * sb, -ca * sb, cb)
    ];
}

function rect(center: number[], axis: Vector[], e0: number, e1: number):
    Rectangle {
    return Rectangle.fromCenterAxisExtent(v(...center), [axis[0], axis[1]],
        v(e0, e1));
}

// The exact squared distance from a point to a solid rectangle: clamp the
// rectangle coordinates of the point to the extents.
function pointRectangleSqr(p: Vector, r: Rectangle): number {
    const delta = sub(p, r.center);
    let closest = r.center.clone();
    for (let i = 0; i < 2; ++i) {
        const e = r.extent.values[i];
        const s = Math.min(Math.max(dot(r.axis[i], delta), -e), e);
        closest = add(closest, mul(s, r.axis[i]));
    }
    const d = sub(p, closest);
    return dot(d, d);
}

// A grid sampling of rectangle0 combined with the exact point-rectangle
// distance to rectangle1 gives an upper bound whose error is at most the
// grid spacing.
function sampledDistance(r0: Rectangle, r1: Rectangle, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= n; ++i) {
        const s0 = r0.extent.values[0] * (-1 + (2 * i) / n);
        const base = add(r0.center, mul(s0, r0.axis[0]));
        for (let j = 0; j <= n; ++j) {
            const s1 = r0.extent.values[1] * (-1 + (2 * j) / n);
            best = Math.min(best,
                pointRectangleSqr(add(base, mul(s1, r0.axis[1])), r1));
        }
    }
    return Math.sqrt(best);
}

function expectConsistent(r0: Rectangle, r1: Rectangle,
    result: DistRectangle3Rectangle3Result): void {
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 8);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);

    const check = (c: [number, number], r: Rectangle, closest: Vector) => {
        for (let i = 0; i < 2; ++i) {
            expect(Math.abs(c[i]))
                .toBeLessThanOrEqual(r.extent.values[i] + 1e-8);
        }
        const p = add(r.center,
            add(mul(c[0], r.axis[0]), mul(c[1], r.axis[1])));
        for (let i = 0; i < 3; ++i) {
            expect(p.values[i]).toBeCloseTo(closest.values[i], 7);
        }
    };
    check(result.cartesian0, r0, result.closest[0]);
    check(result.cartesian1, r1, result.closest[1]);
}

describe('DistRectangle3Rectangle3', () => {
    const query = new DistRectangle3Rectangle3();

    it('computes the distance between parallel overlapping rectangles', () => {
        const r0 = rect([0, 0, 0], frame(0, 0), 1, 1);
        const r1 = rect([0, 0, 7], frame(0, 0), 1, 1);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(7, 10);
        expect(result.closest[0].values[2]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(7, 10);
        expectConsistent(r0, r1, result);
    });

    it('computes the distance between coplanar separated rectangles', () => {
        const r0 = rect([0, 0, 0], frame(0, 0), 1, 1);
        const r1 = rect([5, 0, 0], frame(0, 0), 1, 1);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.cartesian0[0]).toBeCloseTo(1, 8);
        expect(result.cartesian1[0]).toBeCloseTo(-1, 8);
        expectConsistent(r0, r1, result);
    });

    it('reports zero distance for crossing rectangles', () => {
        const r0 = rect([0, 0, 0], [v(1, 0, 0), v(0, 1, 0)], 2, 2);
        const r1 = rect([0, 0, 0], [v(1, 0, 0), v(0, 0, 1)], 2, 2);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(r0, r1, result);
    });

    it('finds a corner-corner closest pair', () => {
        const r0 = rect([0, 0, 0], frame(0, 0), 1, 1);
        const r1 = rect([5, 5, 0], frame(0, 0), 1, 1);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 9), 8);
        expect(result.cartesian0[0]).toBeCloseTo(1, 8);
        expect(result.cartesian0[1]).toBeCloseTo(1, 8);
        expect(result.cartesian1[0]).toBeCloseTo(-1, 8);
        expect(result.cartesian1[1]).toBeCloseTo(-1, 8);
        expectConsistent(r0, r1, result);
    });

    it('is symmetric in its arguments', () => {
        const r0 = rect([0.3, -1.2, 0.4], frame(0.5, 1.1), 1.5, 0.6);
        const r1 = rect([-3.1, 1.4, 2.9], frame(2.0, 0.3), 0.8, 1.7);
        const r01 = query.compute(r0, r1);
        const r10 = query.compute(r1, r0);
        expect(r10.distance).toBeCloseTo(r01.distance, 9);
        expectConsistent(r0, r1, r01);
        expectConsistent(r1, r0, r10);
    });

    it('matches a dense sampling on random inputs', () => {
        let seed = 135792468;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rv = () => [rand() * 6 - 3, rand() * 6 - 3, rand() * 6 - 3];

        for (let trial = 0; trial < 20; ++trial) {
            const r0 = rect(rv(), frame(rand() * Math.PI, rand() * Math.PI),
                rand() * 2 + 0.1, rand() * 2 + 0.1);
            const r1 = rect(rv(), frame(rand() * Math.PI, rand() * Math.PI),
                rand() * 2 + 0.1, rand() * 2 + 0.1);
            const result = query.compute(r0, r1);
            const sampled = Math.min(sampledDistance(r0, r1, 70),
                sampledDistance(r1, r0, 70));
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-8);
            expect(sampled - result.distance).toBeLessThan(0.2);
            expectConsistent(r0, r1, result);
        }
    });
});
