import { describe, expect, it } from 'vitest';
import { DistPointHyperellipsoid } from '../src/DistPointHyperellipsoid.js';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ellipsoid(center: number[], axis: number[][],
    extent: number[]): Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(v(...center),
        axis.map(a => v(...a)), v(...extent));
}

function axisAligned(extent: number[]): Hyperellipsoid {
    const n = extent.length;
    const axes: number[][] = [];
    for (let i = 0; i < n; ++i) {
        const a = new Array<number>(n).fill(0);
        a[i] = 1;
        axes.push(a);
    }
    const center = new Array<number>(n).fill(0);
    return ellipsoid(center, axes, extent);
}

// The value of the hyperellipsoid equation at a point in the coordinate
// system of the hyperellipsoid; it is 1 exactly on the surface.
function equationValue(p: Vector, e: Hyperellipsoid): number {
    const delta = sub(p, e.center);
    let sum = 0;
    for (let i = 0; i < e.dimension; ++i) {
        const y = dot(delta, e.axis[i]) / e.extent.values[i];
        sum += y * y;
    }
    return sum;
}

// Minimum distance from a 2D point to an axis-aligned ellipse, computed by a
// dense angular sampling followed by local refinement.
function bruteForce2D(p: Vector, e0: number, e1: number): number {
    const at = (t: number): number =>
        length(sub(p, v(e0 * Math.cos(t), e1 * Math.sin(t))));
    let best = Number.MAX_VALUE;
    let bt = 0;
    const n = 4096;
    for (let i = 0; i < n; ++i) {
        const t = 2 * Math.PI * i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    let h = 2 * Math.PI / n;
    for (let pass = 0; pass < 80; ++pass) {
        for (const s of [1, -1]) {
            const d = at(bt + s * h);
            if (d < best) {
                best = d;
                bt = bt + s * h;
            }
        }
        h *= 0.7;
    }
    return best;
}

describe('DistPointHyperellipsoid', () => {
    const query = new DistPointHyperellipsoid();

    it('reduces to the circle distance when the extents are equal', () => {
        const circle = axisAligned([2, 2]);
        const p = v(5, 0);
        const result = query.compute(p, circle);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
    });

    it('reduces to the sphere distance in 3D', () => {
        const sphere = axisAligned([3, 3, 3]);
        const p = v(0, 0, 10);
        const result = query.compute(p, sphere);
        expect(result.distance).toBeCloseTo(7, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(3, 10);
    });

    it('returns zero distance for a point on the ellipse', () => {
        const e = axisAligned([2, 1]);
        const p = v(2, 0);
        const result = query.compute(p, e);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
    });

    it('returns the minor semi-axis for the ellipse center', () => {
        // The center is equidistant from the two ends of the minor axis; the
        // algorithm reports one of them.
        const e = axisAligned([3, 1]);
        const result = query.compute(v(0, 0), e);
        expect(result.distance).toBeCloseTo(1, 10);
        expect(Math.abs(result.closest[1].values[1])).toBeCloseTo(1, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 10);
    });

    it('returns the smallest semi-axis for the ellipsoid center', () => {
        const e = axisAligned([4, 3, 2]);
        const result = query.compute(v(0, 0, 0), e);
        expect(result.distance).toBeCloseTo(2, 10);
    });

    it('handles a point inside but off center', () => {
        const e = axisAligned([2, 1]);
        const result = query.compute(v(0.5, 0), e);
        // The closest point is on the minor axis at (0,+-1)? No: the closest
        // ellipse point to an interior point near the center is found by the
        // bisection; verify by the equation and by a brute-force minimum.
        expect(equationValue(result.closest[1], e)).toBeCloseTo(1, 9);
        expect(result.distance).toBeCloseTo(bruteForce2D(v(0.5, 0), 2, 1), 8);
    });

    it('agrees with computeAxisAligned for an axis-aligned hyperellipsoid',
        () => {
            const e = axisAligned([3, 1.5, 0.75]);
            const p = v(1, -2, 4);
            const r0 = query.compute(p, e);
            const r1 = query.computeAxisAligned(p, e.extent);
            expect(r1.distance).toBeCloseTo(r0.distance, 12);
            for (let i = 0; i < 3; ++i) {
                expect(r1.closest[1].values[i])
                    .toBeCloseTo(r0.closest[1].values[i], 12);
            }
        });

    it('handles a rotated and translated ellipse', () => {
        const c = Math.SQRT1_2;
        const e = ellipsoid([1, 2], [[c, c], [-c, c]], [3, 1]);
        // A point on the surface: center + 3 * axis[0].
        const onSurface = add(e.center, mul(3, e.axis[0]));
        const r0 = query.compute(onSurface, e);
        expect(r0.distance).toBeCloseTo(0, 10);

        // A point far along axis[1].
        const outside = add(e.center, mul(5, e.axis[1]));
        const r1 = query.compute(outside, e);
        expect(r1.distance).toBeCloseTo(4, 10);
    });

    it('agrees with a dense sampling for random 2D queries', () => {
        let seed = 112233445;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const e0 = 0.3 + 3 * rand();
            const e1 = 0.3 + 3 * rand();
            const e = axisAligned([e0, e1]);
            const p = v(8 * rand() - 4, 8 * rand() - 4);
            const result = query.compute(p, e);
            expect(result.distance).toBeCloseTo(bruteForce2D(p, e0, e1), 6);
            expect(equationValue(result.closest[1], e)).toBeCloseTo(1, 8);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 10);
        }
    });

    it('produces surface points whose normal is parallel to the offset',
        () => {
            let seed = 998877665;
            const rand = (): number => {
                seed = (1103515245 * seed + 12345) % 2147483648;
                return seed / 2147483648;
            };
            for (let trial = 0; trial < 200; ++trial) {
                const e = axisAligned([0.4 + 2 * rand(), 0.4 + 2 * rand(),
                    0.4 + 2 * rand()]);
                const p = v(8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4);
                const result = query.compute(p, e);
                const x = result.closest[1];
                expect(equationValue(x, e)).toBeCloseTo(1, 7);

                // The gradient of the hyperellipsoid equation at x is
                // 2*(x[i]/e[i]^2); it must be parallel to p - x.
                const grad = v(
                    x.values[0] / (e.extent.values[0] * e.extent.values[0]),
                    x.values[1] / (e.extent.values[1] * e.extent.values[1]),
                    x.values[2] / (e.extent.values[2] * e.extent.values[2]));
                const diff = sub(p, x);
                const lg = length(grad);
                const ld = length(diff);
                if (lg > 1e-8 && ld > 1e-6) {
                    const cosine = dot(grad, diff) / (lg * ld);
                    expect(Math.abs(cosine)).toBeCloseTo(1, 5);
                }
                expect(ld).toBeCloseTo(result.distance, 10);
            }
        });

    it('throws for a dimension mismatch', () => {
        const e = axisAligned([1, 2]);
        expect(() => query.compute(v(1, 2, 3), e)).toThrow();
    });
});
