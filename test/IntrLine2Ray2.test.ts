import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, mul, sub, dot, length } from '../src/Vector.js';
import { IntrLine2Ray2TI, IntrLine2Ray2FI } from '../src/IntrLine2Ray2.js';

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

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, latticeVector, unitVector, wellScaledVector,
    expectVectorClose
} from './helpers/arbitraries.js';
import { IntrLine2Line2FI } from '../src/IntrLine2Line2.js';
import { dotPerp } from '../src/Vector2.js';

const INT32_MAX_V = 2147483647;

// A line and a ray whose directions are transverse with a comfortable margin,
// so the line-line intersection parameter is well conditioned (|t| is bounded
// by roughly |origin difference| / 0.1).
const transverseLineRay = fc.tuple(wellScaledVector(2), unitVector(2),
    wellScaledVector(2), unitVector(2))
    .filter(([, d0, , d1]) => Math.abs(dotPerp(d0, d1)) > 0.1)
    .map(([o0, d0, o1, d1]) => ({
        line: Line.fromOriginDirection(o0, d0),
        ray: Ray.fromOriginDirection(o1, d1)
    }));

describe('IntrLine2Ray2 verification', () => {
    const tiq = new IntrLine2Ray2TI();
    const fiq = new IntrLine2Ray2FI();
    const llq = new IntrLine2Line2FI();

    it('TI and FI agree on intersect and numIntersections', () => {
        check(fc.tuple(wellScaledVector(2), unitVector(2),
            wellScaledVector(2), unitVector(2)), ([o0, d0, o1, d1]) => {
            const l = Line.fromOriginDirection(o0, d0);
            const r = Ray.fromOriginDirection(o1, d1);
            const t = tiq.test(l, r);
            const f = fiq.find(l, r);
            expect(t.intersect).toBe(f.intersect);
            expect(t.numIntersections).toBe(f.numIntersections);
        });
    });

    it('a single-point hit lies on the line and on the ray', () => {
        check(transverseLineRay, ({ line: l, ray: r }) => {
            const f = fiq.find(l, r);
            if (f.numIntersections !== 1) {
                return;
            }
            // Upstream duplicates the single parameter into both entries.
            expect(f.lineParameter[0]).toBe(f.lineParameter[1]);
            expect(f.rayParameter[0]).toBe(f.rayParameter[1]);
            expect(f.rayParameter[0]).toBeGreaterThanOrEqual(0);
            expectVectorClose(f.point,
                add(l.origin, mul(f.lineParameter[0], l.direction)), 0, 0);
            expectVectorClose(f.point,
                add(r.origin, mul(f.rayParameter[0], r.direction)),
                1e-9, 1e-12);
            for (const x of f.point.values) {
                expect(Number.isNaN(x)).toBe(false);
            }
        });
    });

    it('a ray hit is exactly a line-line hit with parameter >= 0', () => {
        check(transverseLineRay, ({ line: l, ray: r }) => {
            const f = fiq.find(l, r);
            const ll = llq.find(l,
                Line.fromOriginDirection(r.origin, r.direction));
            expect(ll.numIntersections).toBe(1);
            const onRay = ll.line1Parameter[0] >= 0;
            expect(f.intersect).toBe(onRay);
            if (onRay) {
                expect(f.lineParameter[0]).toBe(ll.line0Parameter[0]);
                expect(f.rayParameter[0]).toBe(ll.line1Parameter[0]);
            } else {
                // The "no intersection" result keeps the default fields.
                expect(f.numIntersections).toBe(0);
                expect(f.lineParameter).toEqual([0, 0]);
                expect(f.rayParameter).toEqual([0, 0]);
            }
        });
    });

    it('a ray fired at a sampled line point always hits', () => {
        check(fc.tuple(wellScaledVector(2), unitVector(2), unitVector(2),
            fc.double({ min: -5, max: 5, noNaN: true }),
            fc.double({ min: 0.5, max: 5, noNaN: true })),
            ([o, d, u, t, s]) => {
                const l = Line.fromOriginDirection(o, d);
                const target = add(o, mul(t, d));
                // The ray starts away from the target and points at it.
                const r = Ray.fromOriginDirection(sub(target, mul(s, u)), u);
                if (Math.abs(dotPerp(d, u)) <= 0.1) {
                    return;   // near-parallel; the hit is ill conditioned
                }
                const f = fiq.find(l, r);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(1);
                expect(f.rayParameter[0]).toBeGreaterThanOrEqual(0);
                expectVectorClose(f.point, target, 1e-9, 1e-12);
            });
    });

    it('a ray on the line is collinear with the documented parameters', () => {
        // Integer coordinates keep 'o + t*d' exact, so the collinearity test
        // DotPerp(Q, D) == 0 of IntrLine2Line2 is exact as well; with
        // floating-point origins the shifted ray origin is only nearly on the
        // line and upstream (correctly) classifies it as parallel-distinct.
        check(fc.tuple(latticeVector(2), latticeVector(2),
            fc.integer({ min: -5, max: 5 }), fc.boolean())
            .filter(([, d]) => d.values[0] !== 0 || d.values[1] !== 0),
            ([o, d, t, flip]) => {
                const l = Line.fromOriginDirection(o, d);
                const dir = flip ? mul(-1, d) : d;
                const r = Ray.fromOriginDirection(add(o, mul(t, d)), dir);
                const f = fiq.find(l, r);
                const q = tiq.test(l, r);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(INT32_MAX_V);
                expect(q.numIntersections).toBe(INT32_MAX_V);
                expect(f.lineParameter).toEqual([-Number.MAX_VALUE,
                    Number.MAX_VALUE]);
                expect(f.rayParameter).toEqual([0, Number.MAX_VALUE]);
                // point is not filled in the collinear case.
                expect(f.point.values).toEqual([0, 0]);
            });
    });

    it('a ray parallel to but off the line never intersects', () => {
        check(fc.tuple(wellScaledVector(2), unitVector(2),
            fc.double({ min: 0.5, max: 5, noNaN: true }), fc.boolean()),
            ([o, d, offset, flip]) => {
                const l = Line.fromOriginDirection(o, d);
                const n = Vector.fromArray([-d.values[1], d.values[0]]);
                const dir = flip ? mul(-1, d) : d;
                const r = Ray.fromOriginDirection(add(o, mul(offset, n)), dir);
                const f = fiq.find(l, r);
                expect(f.intersect).toBe(false);
                expect(f.numIntersections).toBe(0);
                expect(tiq.test(l, r).intersect).toBe(false);
            });
    });
});
