import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import {
    IntrRay2AlignedBox2TI,
    IntrRay2AlignedBox2FI
} from '../src/IntrRay2AlignedBox2';
import {
    IntrRay2OrientedBox2TI,
    IntrRay2OrientedBox2FI,
    defaultIntrRay2OrientedBox2FIResult
} from '../src/IntrRay2OrientedBox2';
import { OrientedBox } from '../src/OrientedBox';
import { Ray } from '../src/Ray';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function rotatedBox(center: Vector, angle: number, extent: Vector):
    OrientedBox {
    const c = Math.cos(angle), s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(center, [vec(c, s), vec(-s, c)],
        extent);
}

// The signed "outside" amount of P relative to the solid box: negative
// strictly inside, zero on the boundary, positive outside.
function boxSignedDepth(box: OrientedBox, P: Vector): number {
    const diff = sub(P, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 2; ++i) {
        const value = Math.abs(dot(diff, box.axis[i])) - box.extent.values[i];
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

describe('IntrRay2OrientedBox2', () => {
    const ti = new IntrRay2OrientedBox2TI();
    const fi = new IntrRay2OrientedBox2FI();
    const unitAxes = [vec(1, 0), vec(0, 1)];

    it('matches the aligned-box query when the axes are standard', () => {
        const center = vec(1, -2);
        const extent = vec(2, 0.5);
        const obox = OrientedBox.fromCenterAxisExtent(center, unitAxes,
            extent);
        const abox = AlignedBox.fromMinMax(sub(center, extent),
            add(center, extent));
        const aTI = new IntrRay2AlignedBox2TI();
        const aFI = new IntrRay2AlignedBox2FI();

        const rays = [
            ray([-10, -2], [1, 0]),
            ray([1, -2], [0, 1]),
            ray([10, -2], [1, 0]),
            ray([-10, 10], [1, 0]),
            ray([-10, -1.5], [1, 0])
        ];
        for (const r of rays) {
            expect(ti.test(r, obox).intersect)
                .toBe(aTI.test(r, abox).intersect);
            const o = fi.find(r, obox);
            const a = aFI.find(r, abox);
            expect(o.intersect).toBe(a.intersect);
            expect(o.numIntersections).toBe(a.numIntersections);
            for (let i = 0; i < o.numIntersections; ++i) {
                expect(o.parameter[i]).toBeCloseTo(a.parameter[i], 12);
            }
        }
    });

    it('finds the entry and exit points of a rotated box', () => {
        // A box rotated 45 degrees about the origin with extents (1,1) is a
        // diamond with vertices at distance sqrt(2) along the axes.
        const box = rotatedBox(vec(0, 0), Math.PI / 4, vec(1, 1));
        const result = fi.find(ray([-5, 0], [1, 0]), box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(5 - Math.SQRT2, 9);
        expect(result.parameter[1]).toBeCloseTo(5 + Math.SQRT2, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-Math.SQRT2, 9);
        expect(result.point[0].values[1]).toBeCloseTo(0, 9);
        expect(result.point[1].values[0]).toBeCloseTo(Math.SQRT2, 9);
    });

    it('clips at the ray origin when the origin is inside', () => {
        const box = rotatedBox(vec(0, 0), Math.PI / 6, vec(2, 1));
        const result = fi.find(ray([0, 0], [1, 0]), box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeGreaterThan(0);
    });

    it('rejects a box behind the ray', () => {
        const box = rotatedBox(vec(0, 0), 0.3, vec(1, 1));
        const r = ray([5, 0], [1, 0]);
        expect(fi.find(r, box).intersect).toBe(false);
        expect(ti.test(r, box).intersect).toBe(false);
    });

    it('has the documented default result', () => {
        const result = defaultIntrRay2OrientedBox2FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('rejects non-2D boxes', () => {
        const box3 = new OrientedBox(3);
        const r = ray([0, 0], [1, 0]);
        expect(() => ti.test(r, box3)).toThrow();
        expect(() => fi.find(r, box3)).toThrow();
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 606060;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const box = rotatedBox(vec(0.5, -0.25), 0.7, vec(1.5, 0.75));
        for (let trial = 0; trial < 400; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, box);
            expect(ti.test(r, box).intersect).toBe(result.intersect);

            let sampledHit = false;
            for (let k = 0; k <= 2000; ++k) {
                const t = k * 0.01;
                if (boxSignedDepth(box,
                    add(r.origin, mul(t, r.direction))) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            for (let i = 0; i < result.numIntersections; ++i) {
                expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-12);
                const P = add(r.origin, mul(result.parameter[i], r.direction));
                expect(sub(P, result.point[i]).values[0]).toBeCloseTo(0, 9);
                expect(sub(P, result.point[i]).values[1]).toBeCloseTo(0, 9);
                expect(boxSignedDepth(box, P)).toBeLessThan(1e-8);
            }
        }
    });
});
