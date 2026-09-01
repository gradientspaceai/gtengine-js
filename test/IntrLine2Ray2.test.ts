import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line';
import { Ray } from '../src/Ray';
import { Vector, add, mul, sub, dot, length } from '../src/Vector';
import { IntrLine2Ray2TI, IntrLine2Ray2FI } from '../src/IntrLine2Ray2';

const INT32_MAX = 2147483647;
const MAX_T = Number.MAX_VALUE;

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(px: number, py: number, dx: number, dy: number): Line {
    return Line.fromOriginDirection(v2(px, py), v2(dx, dy));
}

function ray(px: number, py: number, dx: number, dy: number): Ray {
    return Ray.fromOriginDirection(v2(px, py), v2(dx, dy));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrLine2Ray2', () => {
    const ti = new IntrLine2Ray2TI();
    const fi = new IntrLine2Ray2FI();

    it('finds the transverse intersection point and parameters', () => {
        // The x axis and the ray from (2,-1) going up cross at (2,0).
        const result = fi.find(line(0, 0, 1, 0), ray(2, -1, 0, 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.lineParameter[0]).toBeCloseTo(2, 12);
        expect(result.lineParameter[1]).toBe(result.lineParameter[0]);
        expect(result.rayParameter[0]).toBeCloseTo(1, 12);
        expect(result.rayParameter[1]).toBe(result.rayParameter[0]);
        expect(result.point.values[0]).toBeCloseTo(2, 12);
        expect(result.point.values[1]).toBeCloseTo(0, 12);
        expect(ti.test(line(0, 0, 1, 0), ray(2, -1, 0, 1)).numIntersections)
            .toBe(1);
    });

    it('rejects an intersection behind the ray origin', () => {
        // The lines cross at (2,0) but the ray points away from it.
        const l = line(0, 0, 1, 0);
        const r = ray(2, -1, 0, -1);
        expect(fi.find(l, r).intersect).toBe(false);
        expect(fi.find(l, r).numIntersections).toBe(0);
        expect(ti.test(l, r).intersect).toBe(false);
        expect(ti.test(l, r).numIntersections).toBe(0);
    });

    it('reports intersection when the ray origin is on the line', () => {
        const l = line(0, 0, 1, 0);
        const r = ray(3, 0, 0, 1);
        const result = fi.find(l, r);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.rayParameter[0]).toBe(0);
        expect(result.lineParameter[0]).toBeCloseTo(3, 12);
    });

    it('reports no intersection for parallel but distinct components', () => {
        const l = line(0, 0, 1, 0);
        const r = ray(0, 1, 1, 0);
        expect(ti.test(l, r).intersect).toBe(false);
        expect(ti.test(l, r).numIntersections).toBe(0);
        const result = fi.find(l, r);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports the collinear case with the documented parameters', () => {
        const l = line(0, 0, 1, 0);
        const r = ray(3, 0, 1, 0);
        const tiResult = ti.test(l, r);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(INT32_MAX);

        const result = fi.find(l, r);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.lineParameter).toEqual([-MAX_T, MAX_T]);
        expect(result.rayParameter).toEqual([0, MAX_T]);
        expect(result.point.values).toEqual([0, 0]);
    });

    it('reports the collinear case for an opposite ray direction', () => {
        const l = line(0, 0, 1, 0);
        const r = ray(3, 0, -1, 0);
        expect(ti.test(l, r).numIntersections).toBe(INT32_MAX);
        expect(fi.find(l, r).numIntersections).toBe(INT32_MAX);
    });

    it('agrees between TI and FI on random configurations', () => {
        const rnd = makeRandom(8675309);
        let numHit = 0, numMiss = 0;
        for (let k = 0; k < 500; ++k) {
            const l = line(rnd() * 6 - 3, rnd() * 6 - 3,
                rnd() * 2 - 1, rnd() * 2 - 1);
            const r = ray(rnd() * 6 - 3, rnd() * 6 - 3,
                rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(l.direction, l.direction) < 1e-8 ||
                dot(r.direction, r.direction) < 1e-8) {
                continue;
            }
            const tiResult = ti.test(l, r);
            const fiResult = fi.find(l, r);
            expect(tiResult.intersect).toBe(fiResult.intersect);
            expect(tiResult.numIntersections).toBe(fiResult.numIntersections);

            if (fiResult.numIntersections === 1) {
                // The reported point lies on both components.
                const pl = add(l.origin,
                    mul(fiResult.lineParameter[0], l.direction));
                const pr = add(r.origin,
                    mul(fiResult.rayParameter[0], r.direction));
                expect(length(sub(pl, fiResult.point))).toBeLessThan(1e-9);
                expect(length(sub(pr, fiResult.point))).toBeLessThan(1e-9);
                expect(fiResult.rayParameter[0]).toBeGreaterThanOrEqual(0);
                ++numHit;
            } else {
                ++numMiss;
            }
        }
        expect(numHit).toBeGreaterThan(0);
        expect(numMiss).toBeGreaterThan(0);
    });

    it('agrees with a brute-force sampling of the ray', () => {
        const rnd = makeRandom(20250101);
        for (let k = 0; k < 200; ++k) {
            const l = line(rnd() * 4 - 2, rnd() * 4 - 2,
                rnd() * 2 - 1, rnd() * 2 - 1);
            const r = ray(rnd() * 4 - 2, rnd() * 4 - 2,
                rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(l.direction, l.direction) < 1e-8 ||
                dot(r.direction, r.direction) < 1e-8) {
                continue;
            }
            const intersect = ti.test(l, r).intersect;

            // The signed side of the line at ray points t=0 and t=large. A
            // sign change means the ray crosses the line.
            const sideAt = (t: number): number => {
                const p = sub(add(r.origin, mul(t, r.direction)), l.origin);
                return p.values[0] * l.direction.values[1]
                    - p.values[1] * l.direction.values[0];
            };
            const s0 = sideAt(0);
            const s1 = sideAt(1e6);
            if (s0 * s1 < 0) {
                expect(intersect).toBe(true);
            }
        }
    });
});
