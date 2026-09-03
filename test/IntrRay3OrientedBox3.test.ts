import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrRay3AlignedBox3TI,
    IntrRay3AlignedBox3FI
} from '../src/IntrRay3AlignedBox3.js';
import {
    IntrRay3OrientedBox3TI,
    IntrRay3OrientedBox3FI,
    defaultIntrRay3OrientedBox3FIResult
} from '../src/IntrRay3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function orthonormalFrame(w: Vector): Vector[] {
    const v = [w.clone(), Vector.zero(3), Vector.zero(3)];
    normalize(v[0]);
    computeOrthogonalComplement3(1, v, false);
    return [v[0], v[1], v[2]];
}

function boxSignedDepth(box: OrientedBox, P: Vector): number {
    const diff = sub(P, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 3; ++i) {
        const value = Math.abs(dot(diff, box.axis[i])) - box.extent.values[i];
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

describe('IntrRay3OrientedBox3', () => {
    const ti = new IntrRay3OrientedBox3TI();
    const fi = new IntrRay3OrientedBox3FI();
    const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

    it('matches the aligned-box query when the axes are standard', () => {
        const center = vec(1, -2, 3);
        const extent = vec(2, 0.5, 1);
        const obox = OrientedBox.fromCenterAxisExtent(center, unitAxes,
            extent);
        const abox = AlignedBox.fromMinMax(sub(center, extent),
            add(center, extent));
        const aTI = new IntrRay3AlignedBox3TI();
        const aFI = new IntrRay3AlignedBox3FI();

        const rays = [
            ray([-10, -2, 3], [1, 0, 0]),
            ray([1, -2, 3], [0, 1, 0]),
            ray([10, -2, 3], [1, 0, 0]),
            ray([-10, 10, 3], [1, 0, 0]),
            ray([0, 0, 0], [1, -2, 3])
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

    it('finds entry and exit points of a rotated box', () => {
        // Rotate 45 degrees about z; the box with extents (1,1,1) becomes a
        // diamond cross-section reaching sqrt(2) along x.
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0),
            [vec(c, c, 0), vec(-c, c, 0), vec(0, 0, 1)], vec(1, 1, 1));
        const result = fi.find(ray([-5, 0, 0], [1, 0, 0]), box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(5 - Math.SQRT2, 9);
        expect(result.parameter[1]).toBeCloseTo(5 + Math.SQRT2, 9);
    });

    it('clips at the ray origin when the origin is inside', () => {
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0), unitAxes,
            vec(2, 1, 1));
        const result = fi.find(ray([0, 0, 0], [1, 0, 0]), box);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
    });

    it('rejects a box behind the ray', () => {
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0), unitAxes,
            vec(1, 1, 1));
        const r = ray([5, 0, 0], [1, 0, 0]);
        expect(fi.find(r, box).intersect).toBe(false);
        expect(ti.test(r, box).intersect).toBe(false);
    });

    it('has the documented default result', () => {
        const result = defaultIntrRay3OrientedBox3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('rejects non-3D boxes', () => {
        const box2 = new OrientedBox(2);
        const r = ray([0, 0, 0], [1, 0, 0]);
        expect(() => ti.test(r, box2)).toThrow();
        expect(() => fi.find(r, box2)).toThrow();
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 8080808;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const axes = orthonormalFrame(vec(0.3, 0.8, -0.5));
        const box = OrientedBox.fromCenterAxisExtent(vec(0.5, -0.25, 1), axes,
            vec(1.5, 0.75, 1));

        for (let trial = 0; trial < 300; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, box);
            expect(ti.test(r, box).intersect).toBe(result.intersect);

            let sampledHit = false;
            for (let k = 0; k <= 1500; ++k) {
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

            if (result.intersect) {
                for (let i = 0; i < 2; ++i) {
                    expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-12);
                    const P = add(r.origin,
                        mul(result.parameter[i], r.direction));
                    expect(sub(P, result.point[i]).values[0])
                        .toBeCloseTo(0, 9);
                    expect(boxSignedDepth(box, P)).toBeLessThan(1e-8);
                }
            }
        }
    });
});
