import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import {
    IntrLine3AlignedBox3FI
} from '../src/IntrLine3AlignedBox3.js';
import { Line } from '../src/Line.js';
import {
    IntrRay3AlignedBox3TI,
    IntrRay3AlignedBox3FI
} from '../src/IntrRay3AlignedBox3.js';

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
    for (let i = 0; i < 3; ++i) {
        if (x.values[i] < b.min.values[i] - 1e-12 ||
            x.values[i] > b.max.values[i] + 1e-12) {
            return false;
        }
    }
    return true;
}

describe('IntrRay3AlignedBox3', () => {
    const ti = new IntrRay3AlignedBox3TI();
    const fi = new IntrRay3AlignedBox3FI();
    const unit = box([-1, -1, -1], [1, 1, 1]);

    it('clips the near end when the ray origin is inside the box', () => {
        const r = ray([0, 0, 0], [0, 0, 1]);
        expect(ti.test(r, unit).intersect).toBe(true);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values[2]).toBeCloseTo(1, 12);
    });

    it('finds both crossings for a ray that starts outside', () => {
        const r = ray([0, 0, -5], [0, 0, 1]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[2]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports no intersection when the ray points away from the box', () => {
        const r = ray([0, 0, -5], [0, 0, -1]);
        expect(ti.test(r, unit).intersect).toBe(false);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('reports the touching point when the ray origin is on the far face', () => {
        // The origin is on the +z face and the ray points outward, so the
        // overlap with [0,+infinity) is the single parameter 0.
        const r = ray([0, 0, 1], [0, 0, 1]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(0, 12);
        expect(ti.test(r, unit).intersect).toBe(true);
    });

    it('handles a ray grazing an edge of the box', () => {
        const r = ray([1, 1, -5], [0, 0, 1]);
        const result = fi.find(r, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
    });

    it('is the line query clipped to t >= 0', () => {
        const rand = makeRandom(4242);
        const lineFI = new IntrLine3AlignedBox3FI();
        const b = box([-1, -2, 0.5], [2, 1, 3]);
        for (let trial = 0; trial < 300; ++trial) {
            const r = ray(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
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
        const rand = makeRandom(99);
        const b = box([-1, -2, 0.5], [2, 1, 3]);
        for (let trial = 0; trial < 60; ++trial) {
            const r = ray(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
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
