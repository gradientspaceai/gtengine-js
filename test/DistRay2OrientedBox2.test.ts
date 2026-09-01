import { describe, expect, it } from 'vitest';
import { DistRay2OrientedBox2 } from '../src/DistRay2OrientedBox2';
import { OrientedBox } from '../src/OrientedBox';
import { Ray } from '../src/Ray';
import { Vector, add, dot, length, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(center: number[], axis: number[][],
    extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(v(...center),
        axis.map(a => v(...a)), v(...extent));
}

// The exact squared distance from a point to the solid oriented box, computed
// independently of the library.
function pointBoxSqrDistance(p: Vector, b: OrientedBox): number {
    const delta = sub(p, b.center);
    let sqrDistance = 0;
    for (let i = 0; i < p.size; ++i) {
        const y = dot(delta, b.axis[i]);
        const excess = Math.max(0, Math.abs(y) - b.extent.values[i]);
        sqrDistance += excess * excess;
    }
    return sqrDistance;
}

// Verify that the reported box point is inside the box.
function verifyBoxPoint(b: OrientedBox, q: Vector): void {
    const delta = sub(q, b.center);
    for (let i = 0; i < q.size; ++i) {
        expect(Math.abs(dot(delta, b.axis[i])))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
    }
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

// The distance from a ray point to the solid box is a convex function of the
// ray parameter, so a ternary search over [0,tmax] converges to the true
// minimum independently of the query under test.
function bruteForce(r: Ray, b: OrientedBox, tmax: number): number {
    const f = (t: number): number =>
        pointBoxSqrDistance(add(r.origin, mul(t, r.direction)), b);
    let lo = 0;
    let hi = tmax;
    for (let i = 0; i < 200; ++i) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) < f(m1)) {
            hi = m1;
        }
        else {
            lo = m0;
        }
    }
    return Math.sqrt(f(0.5 * (lo + hi)));
}

// Verify that the reported closest points are consistent with the reported
// distance and lie on their primitives.
function verifyClosest(r: Ray, b: OrientedBox,
    result: { distance: number, parameter: number, closest: [Vector, Vector] }
): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    const onRay = add(r.origin, mul(result.parameter, r.direction));
    for (let i = 0; i < onRay.size; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onRay.values[i], 9);
    }
    verifyBoxPoint(b, result.closest[1]);
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistRay2OrientedBox2', () => {
    const query = new DistRay2OrientedBox2();
    const axisAligned = box([0, 0], [[1, 0], [0, 1]], [1, 1]);
    const c = Math.SQRT1_2;
    const rotated = box([0, 0], [[c, c], [-c, c]], [1, 1]);

    it('matches the aligned result for an axis-aligned oriented box', () => {
        const r = ray([3, 0], [1, 0]);
        const result = query.compute(r, axisAligned);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
    });

    it('handles a box rotated by 45 degrees', () => {
        // The rotated unit box has vertices at distance sqrt(2) from the
        // origin along the coordinate axes; (sqrt(2),0) is a vertex.
        const r = ray([4, 0], [1, 0]);
        const result = query.compute(r, rotated);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 10);
        verifyClosest(r, rotated, result);
    });

    it('reports zero distance when the ray enters the box', () => {
        const r = ray([-5, 0.25], [1, 0]);
        const result = query.compute(r, rotated);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('handles a degenerate box with zero extents', () => {
        const b = box([2, 1], [[1, 0], [0, 1]], [0, 0]);
        const r = ray([0, 1], [-1, 0]);
        const result = query.compute(r, b);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.parameter).toBe(0);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 24681357;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const angle = 2 * Math.PI * rand();
            const ca = Math.cos(angle);
            const sa = Math.sin(angle);
            const b = box([2 * rand() - 1, 2 * rand() - 1],
                [[ca, sa], [-sa, ca]],
                [0.2 + 2 * rand(), 0.2 + 2 * rand()]);
            const r = ray([8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-6) {
                continue;
            }
            const result = query.compute(r, b);
            const brute = bruteForce(r, b, 1e6);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeCloseTo(brute, 6);
            verifyClosest(r, b, result);
        }
    });
});
