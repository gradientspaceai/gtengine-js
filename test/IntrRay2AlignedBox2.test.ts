import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import { IntrLine2AlignedBox2FI } from '../src/IntrLine2AlignedBox2.js';
import {
    IntrRay2AlignedBox2TI,
    IntrRay2AlignedBox2FI
} from '../src/IntrRay2AlignedBox2.js';
import {
    intrRay2AlignedBox2TIDoQuery,
    intrRay2AlignedBox2FIDoQuery
} from '../src/IntrRay2AlignedBox2.js';
import {
    defaultIntrLine2AlignedBox2TIResult,
    defaultIntrLine2AlignedBox2FIResult
} from '../src/IntrLine2AlignedBox2.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(vec(min), vec(max));
}

function ray(p: number[], d: number[]): Ray {
    const dir = vec(d);
    normalize(dir);
    return Ray.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function insideBox(b: AlignedBox, x: Vector): boolean {
    for (let i = 0; i < 2; ++i) {
        if (x.values[i] < b.min.values[i] - 1e-12 ||
            x.values[i] > b.max.values[i] + 1e-12) {
            return false;
        }
    }
    return true;
}

describe('IntrRay2AlignedBox2', () => {
    const ti = new IntrRay2AlignedBox2TI();
    const fi = new IntrRay2AlignedBox2FI();
    const unit = box([-1, -1], [1, 1]);

    it('finds both crossings for a ray aimed at the box', () => {
        const r = ray([-5, 0], [1, 0]);
        expect(ti.test(r, unit).intersect).toBe(true);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips the near end when the ray origin is inside the box', () => {
        const r = ray([0.5, 0], [1, 0]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
    });

    it('reports no intersection when the ray points away', () => {
        const r = ray([-5, 0], [-1, 0]);
        expect(ti.test(r, unit).intersect).toBe(false);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports a single touching point on the exit edge', () => {
        const r = ray([1, 0], [1, 0]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a ray that grazes a corner', () => {
        // The ray runs along x = 1, touching the right edge of the box.
        const r = ray([1, -5], [0, 1]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
    });

    it('is the line query clipped to t >= 0', () => {
        const rand = makeRandom(20260901);
        const lineFI = new IntrLine2AlignedBox2FI();
        const b = box([-1, 0.5], [2, 3]);
        for (let trial = 0; trial < 400; ++trial) {
            const r = ray([8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1]);
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lineResult = lineFI.find(l, b);
            const rayResult = fi.find(r, b);
            expect(ti.test(r, b).intersect).toBe(rayResult.intersect);

            if (!lineResult.intersect || lineResult.parameter[1] < 0) {
                expect(rayResult.intersect).toBe(false);
            }
            else {
                expect(rayResult.intersect).toBe(true);
                expect(rayResult.parameter[0]).toBeCloseTo(
                    Math.max(lineResult.parameter[0], 0), 12);
                expect(rayResult.parameter[1]).toBeCloseTo(
                    lineResult.parameter[1], 12);
            }
        }
    });

    it('agrees with dense sampling along the ray', () => {
        const rand = makeRandom(31337);
        const b = box([-1, 0.5], [2, 3]);
        for (let trial = 0; trial < 80; ++trial) {
            const r = ray([6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(r, b);

            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = (12 * k) / n;
                const x = add(r.origin, mul(t, r.direction));
                if (insideBox(b, x)) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(2e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(2e-3);
            }
        }
    });
});

describe('intrRay2AlignedBox2 DoQuery helpers', () => {
    // The helpers take the ray in the box-centered coordinate system. The box
    // below is already centered at the origin, so no translation is needed.
    const b = box([-2, -1], [2, 1]);
    const extent = vec([2, 1]);

    it('the TI helper agrees with the class query', () => {
        const cases: Array<[number, number, number, number]> = [
            [-5, 0, 1, 0],
            [5, 0, 1, 0],
            [0, 0, 0, 1],
            [-5, 3, 1, 0],
            [-4, -3, 1, 1]
        ];
        for (const [px, py, dx, dy] of cases) {
            const d = vec([dx, dy]);
            normalize(d);
            const result = defaultIntrLine2AlignedBox2TIResult();
            intrRay2AlignedBox2TIDoQuery(vec([px, py]), d, extent, result);
            const expected = new IntrRay2AlignedBox2TI().test(
                Ray.fromOriginDirection(vec([px, py]), d), b);
            expect(result.intersect).toBe(expected.intersect);
        }
    });

    it('the FI helper clips to the ray t-interval and leaves points alone', () => {
        const result = defaultIntrLine2AlignedBox2FIResult();
        intrRay2AlignedBox2FIDoQuery(vec([0, 0]), vec([1, 0]), extent, result);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The ray starts inside the box, so t0 is clamped to 0.
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('the FI helper rejects a ray pointing away from the box', () => {
        const result = defaultIntrLine2AlignedBox2FIResult();
        intrRay2AlignedBox2FIDoQuery(vec([-5, 0]), vec([-1, 0]), extent,
            result);
        expect(result.intersect).toBe(false);
    });
});
