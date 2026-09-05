import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Vector, add, mul, sub, length } from '../src/Vector.js';
import { IntrLine2Line2TI, IntrLine2Line2FI } from '../src/IntrLine2Line2.js';
import { check, expectVectorClose, fc } from './helpers/arbitraries.js';

function line(px: number, py: number, dx: number, dy: number): Line {
    return Line.fromOriginDirection(Vector.fromArray([px, py]),
        Vector.fromArray([dx, dy]));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

const INT32_MAX = 2147483647;

describe('IntrLine2Line2', () => {
    const ti = new IntrLine2Line2TI();
    const fi = new IntrLine2Line2FI();

    it('finds the single intersection of two transverse lines', () => {
        // x axis and the line through (2,-1) with direction (0,1).
        const l0 = line(0, 0, 1, 0);
        const l1 = line(2, -1, 0, 1);
        expect(ti.test(l0, l1).numIntersections).toBe(1);

        const result = fi.find(l0, l1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point.values[0]).toBeCloseTo(2, 12);
        expect(result.point.values[1]).toBeCloseTo(0, 12);
        expect(result.line0Parameter).toEqual([2, 2]);
        expect(result.line1Parameter).toEqual([1, 1]);
    });

    it('handles non-unit-length directions', () => {
        // Directions need not be unit length; the parameters scale with them.
        const l0 = line(0, 0, 2, 0);
        const l1 = line(3, -4, 0, 4);
        const result = fi.find(l0, l1);
        expect(result.numIntersections).toBe(1);
        expect(result.line0Parameter[0]).toBeCloseTo(1.5, 12);
        expect(result.line1Parameter[0]).toBeCloseTo(1, 12);
        expect(result.point.values[0]).toBeCloseTo(3, 12);
        expect(result.point.values[1]).toBeCloseTo(0, 12);
    });

    it('reports no intersection for parallel distinct lines', () => {
        const l0 = line(0, 0, 1, 1);
        const l1 = line(0, 1, 2, 2);
        const t = ti.test(l0, l1);
        expect(t.intersect).toBe(false);
        expect(t.numIntersections).toBe(0);

        const f = fi.find(l0, l1);
        expect(f.intersect).toBe(false);
        expect(f.numIntersections).toBe(0);
        expect(f.line0Parameter).toEqual([0, 0]);
        expect(f.point.values).toEqual([0, 0]);
    });

    it('reports infinitely many intersections for identical lines', () => {
        const l0 = line(1, 2, 3, 4);
        // The same line with a different origin and a scaled direction.
        const l1 = line(1 + 6, 2 + 8, -1.5, -2);
        const t = ti.test(l0, l1);
        expect(t.intersect).toBe(true);
        expect(t.numIntersections).toBe(INT32_MAX);

        const f = fi.find(l0, l1);
        expect(f.intersect).toBe(true);
        expect(f.numIntersections).toBe(INT32_MAX);
        expect(f.line0Parameter).toEqual([-Number.MAX_VALUE, Number.MAX_VALUE]);
        expect(f.line1Parameter).toEqual([-Number.MAX_VALUE, Number.MAX_VALUE]);
    });

    it('reports the same line when the two origins coincide', () => {
        const l0 = line(5, 5, 1, 0);
        const l1 = line(5, 5, -3, 0);
        expect(fi.find(l0, l1).numIntersections).toBe(INT32_MAX);
    });

    it('is consistent between TI and FI on random line pairs', () => {
        const rand = makeRandom(1357911);
        let numSingle = 0, numParallel = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const l0 = line(4 * rand() - 2, 4 * rand() - 2,
                2 * rand() - 1, 2 * rand() - 1);
            const l1 = line(4 * rand() - 2, 4 * rand() - 2,
                2 * rand() - 1, 2 * rand() - 1);
            if (length(l0.direction) < 1e-6 || length(l1.direction) < 1e-6) {
                continue;
            }

            const t = ti.test(l0, l1);
            const f = fi.find(l0, l1);
            expect(f.intersect).toBe(t.intersect);
            expect(f.numIntersections).toBe(t.numIntersections);

            if (f.numIntersections === 1) {
                ++numSingle;
                // The point lies on both lines.
                const p0 = add(l0.origin, mul(f.line0Parameter[0], l0.direction));
                const p1 = add(l1.origin, mul(f.line1Parameter[0], l1.direction));
                expect(length(sub(p0, f.point))).toBeLessThan(1e-9);
                expect(length(sub(p1, f.point))).toBeLessThan(1e-9);
            } else {
                ++numParallel;
            }
        }
        expect(numSingle).toBeGreaterThan(100);
        // Random directions are essentially never parallel.
        expect(numParallel).toBe(0);
    });

    it('detects parallel and identical lines built to be so', () => {
        // Integer inputs so the DotPerp comparisons against zero are exact
        // (the queries test for exact equality, as upstream does).
        const rand = makeRandom(24681012);
        let numSame = 0, numDistinct = 0;
        const pick = (n: number) => Math.floor((2 * n + 1) * rand()) - n;
        for (let trial = 0; trial < 300; ++trial) {
            const dx = pick(8), dy = pick(8);
            const px = pick(8), py = pick(8);
            const scale = pick(4);
            if (scale === 0 || (dx === 0 && dy === 0)) {
                continue;
            }
            const l0 = line(px, py, dx, dy);

            // A point further along the same line, with a scaled direction.
            const s = Math.floor(4 * rand());
            const same = line(px + s * dx, py + s * dy, scale * dx, scale * dy);
            expect(fi.find(l0, same).numIntersections).toBe(INT32_MAX);
            expect(ti.test(l0, same).numIntersections).toBe(INT32_MAX);
            ++numSame;

            // A parallel but distinct line: offset along the perpendicular.
            const distinct = line(px + dy, py - dx, scale * dx, scale * dy);
            expect(fi.find(l0, distinct).numIntersections).toBe(0);
            expect(ti.test(l0, distinct).intersect).toBe(false);
            ++numDistinct;
        }
        expect(numSame).toBeGreaterThan(100);
        expect(numDistinct).toBeGreaterThan(100);
    });
});

describe('IntrLine2Line2 verification', () => {
    const ti = new IntrLine2Line2TI();
    const fi = new IntrLine2Line2FI();

    // Integer origins and directions: DotPerp is a difference of products of
    // small integers, hence exact, so the parallel/same/transverse
    // classification below is decided without any round-off.
    const latticeLine = fc.tuple(fc.integer({ min: -5, max: 5 }),
        fc.integer({ min: -5, max: 5 }), fc.integer({ min: -3, max: 3 }),
        fc.integer({ min: -3, max: 3 }))
        .filter(([, , dx, dy]) => dx !== 0 || dy !== 0)
        .map(([px, py, dx, dy]) => line(px, py, dx, dy));

    it('classification matches the exact DotPerp signs on lattice lines', () => {
        check(fc.tuple(latticeLine, latticeLine), ([l0, l1]) => {
            const d0 = l0.direction.values, d1 = l1.direction.values;
            const p0 = l0.origin.values, p1 = l1.origin.values;
            const cross = d0[0] * d1[1] - d0[1] * d1[0];
            const qx = p1[0] - p0[0], qy = p1[1] - p0[1];
            const qPerpD1 = qx * d1[1] - qy * d1[0];
            const t = ti.test(l0, l1);
            const f = fi.find(l0, l1);
            expect(f.intersect).toBe(t.intersect);
            expect(f.numIntersections).toBe(t.numIntersections);
            if (cross !== 0) {
                expect(f.numIntersections).toBe(1);
            } else if (qPerpD1 !== 0) {
                expect(f.intersect).toBe(false);
                expect(f.numIntersections).toBe(0);
            } else {
                expect(f.numIntersections).toBe(INT32_MAX);
                expect(f.line0Parameter).toEqual([-Number.MAX_VALUE,
                    Number.MAX_VALUE]);
                expect(f.line1Parameter).toEqual([-Number.MAX_VALUE,
                    Number.MAX_VALUE]);
            }
        });
    });

    it('the reported point lies on both lines at the reported parameters', () => {
        check(fc.tuple(latticeLine, latticeLine), ([l0, l1]) => {
            const f = fi.find(l0, l1);
            if (f.numIntersections !== 1) {
                return;
            }
            const [s0, s0b] = f.line0Parameter;
            const [s1, s1b] = f.line1Parameter;
            expect(s0b).toBe(s0);
            expect(s1b).toBe(s1);
            expect(Number.isFinite(s0)).toBe(true);
            expect(Number.isFinite(s1)).toBe(true);
            expectVectorClose(f.point, add(l0.origin, mul(s0, l0.direction)),
                0, 0);
            // The two parameters share one exactly computed denominator, so
            // the two reconstructions agree to a few ulps of the coordinates.
            const other = add(l1.origin, mul(s1, l1.direction));
            const scale = 1 + length(l0.origin) + length(l1.origin)
                + Math.abs(s0) + Math.abs(s1);
            expect(length(sub(f.point, other)))
                .toBeLessThanOrEqual(1e-12 * scale);
        });
    });

    it('is symmetric under argument swap', () => {
        check(fc.tuple(latticeLine, latticeLine), ([l0, l1]) => {
            const a = fi.find(l0, l1);
            const b = fi.find(l1, l0);
            expect(b.numIntersections).toBe(a.numIntersections);
            expect(b.intersect).toBe(a.intersect);
            expect(ti.test(l1, l0).numIntersections)
                .toBe(ti.test(l0, l1).numIntersections);
            if (a.numIntersections !== 1) {
                return;
            }
            // The roles of the two parameters swap with the arguments.
            expect(b.line0Parameter[0]).toBeCloseTo(a.line1Parameter[0], 9);
            expect(b.line1Parameter[0]).toBeCloseTo(a.line0Parameter[0], 9);
            const scale = 1 + length(a.point);
            expect(length(sub(a.point, b.point)))
                .toBeLessThanOrEqual(1e-12 * scale);
        });
    });

    it('non-unit-length directions rescale the parameters but not the point', () => {
        check(fc.tuple(latticeLine, latticeLine, fc.integer({ min: 2, max: 5 })),
            ([l0, l1, k]) => {
                const a = fi.find(l0, l1);
                if (a.numIntersections !== 1) {
                    return;
                }
                const scaled = Line.fromOriginDirection(l0.origin,
                    mul(k, l0.direction));
                const b = fi.find(scaled, l1);
                expect(b.numIntersections).toBe(1);
                expect(b.line0Parameter[0] * k).toBeCloseTo(a.line0Parameter[0], 9);
                const scale = 1 + length(a.point);
                expect(length(sub(a.point, b.point)))
                    .toBeLessThanOrEqual(1e-12 * scale);
            });
    });

    it('coincident and parallel-distinct lines are classified exactly', () => {
        const base = line(1, 2, 3, 4);
        // The same line with a shifted origin and a scaled direction.
        const same = line(1 + 6, 2 + 8, -6, -8);
        expect(fi.find(base, same).numIntersections).toBe(INT32_MAX);
        expect(ti.test(base, same).numIntersections).toBe(INT32_MAX);
        const parallel = line(2, 2, 3, 4);
        expect(fi.find(base, parallel).intersect).toBe(false);
        expect(ti.test(base, parallel).intersect).toBe(false);
    });
});
