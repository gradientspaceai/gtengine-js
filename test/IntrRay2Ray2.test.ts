import { describe, it, expect } from 'vitest';
import { Ray } from '../src/Ray.js';
import { Vector, add, mul, normalize, sub, length } from '../src/Vector.js';
import { IntrRay2Ray2TI, IntrRay2Ray2FI } from '../src/IntrRay2Ray2.js';

const INT32_MAX = 2147483647;
const MAX_T = Number.MAX_VALUE;

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
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

describe('IntrRay2Ray2', () => {
    const ti = new IntrRay2Ray2TI();
    const fi = new IntrRay2Ray2FI();

    it('finds the single crossing of two transversal rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([2, -1], [0, 1]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(r0, r1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.ray0Parameter[0]).toBeCloseTo(2, 12);
        expect(result.ray1Parameter[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('rejects a crossing that lies behind one of the rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([-2, -1], [0, 1]);
        expect(ti.test(r0, r1).intersect).toBe(false);
        const result = fi.find(r0, r1);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports no intersection for parallel distinct rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([0, 1], [1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(false);
        expect(tiResult.numIntersections).toBe(0);
        expect(fi.find(r0, r1).intersect).toBe(false);
    });

    it('reports a ray of intersection for collinear same-direction rays', () => {
        // ray1 starts ahead of ray0.
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([3, 0], [1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(INT32_MAX);

        const result = fi.find(r0, r1);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.ray0Parameter).toEqual([3, MAX_T]);
        expect(result.ray1Parameter).toEqual([0, MAX_T]);
        expect(result.point[0].values).toEqual([3, 0]);
    });

    it('handles collinear same-direction rays with ray1 behind ray0', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([-3, 0], [1, 0]);
        const result = fi.find(r0, r1);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.ray0Parameter).toEqual([0, MAX_T]);
        expect(result.ray1Parameter).toEqual([3, MAX_T]);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('reports a segment for overlapping opposite collinear rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([4, 0], [-1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(2);

        const result = fi.find(r0, r1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.ray0Parameter).toEqual([0, 4]);
        expect(result.ray1Parameter).toEqual([0, 4]);
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([4, 0]);
    });

    it('reports no intersection for disjoint opposite collinear rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([-4, 0], [-1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(false);
        expect(tiResult.numIntersections).toBe(0);
        expect(fi.find(r0, r1).intersect).toBe(false);
    });

    it('handles opposite collinear rays that share an origin', () => {
        // Upstream reports numIntersections = 1 from the test query (t == 0)
        // but 2 from the find query, whose t >= 0 branch produces the
        // degenerate segment [0,0]. The port preserves both behaviors.
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([0, 0], [-1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(r0, r1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.ray0Parameter).toEqual([0, 0]);
        expect(result.ray1Parameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('agrees between the test and find queries on random configurations', () => {
        const rand = makeRandom(13579);
        let single = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const r0 = ray([6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1]);
            const r1 = ray([6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1]);
            const tiResult = ti.test(r0, r1);
            const fiResult = fi.find(r0, r1);
            expect(tiResult.intersect).toBe(fiResult.intersect);
            expect(tiResult.numIntersections).toBe(fiResult.numIntersections);

            if (fiResult.numIntersections === 1) {
                ++single;
                expect(fiResult.ray0Parameter[0]).toBeGreaterThanOrEqual(0);
                expect(fiResult.ray1Parameter[0]).toBeGreaterThanOrEqual(0);
                const p0 = add(r0.origin,
                    mul(fiResult.ray0Parameter[0], r0.direction));
                const p1 = add(r1.origin,
                    mul(fiResult.ray1Parameter[0], r1.direction));
                expect(length(sub(p0, p1))).toBeLessThan(1e-9);
                expect(length(sub(p0, fiResult.point[0]))).toBeLessThan(1e-9);
            }
        }
        expect(single).toBeGreaterThan(50);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, latticeVector, unitVector, wellScaledVector, expectVectorClose
} from './helpers/arbitraries.js';
import { IntrLine2Line2FI } from '../src/IntrLine2Line2.js';
import { Line } from '../src/Line.js';
import { dotPerp } from '../src/Vector2.js';

const INT32_MAX_V = 2147483647;

// Two rays whose directions are transverse with a comfortable margin.
const transverseRays = fc.tuple(wellScaledVector(2), unitVector(2),
    wellScaledVector(2), unitVector(2))
    .filter(([, d0, , d1]) => Math.abs(dotPerp(d0, d1)) > 0.1)
    .map(([o0, d0, o1, d1]) => ({
        ray0: Ray.fromOriginDirection(o0, d0),
        ray1: Ray.fromOriginDirection(o1, d1)
    }));

// An exactly representable unit direction, so that collinear configurations
// built from integer offsets stay exactly collinear.
const axisDirection = fc.integer({ min: 0, max: 3 }).map(i =>
    Vector.fromArray([[1, 0], [0, 1], [-1, 0], [0, -1]][i]));

describe('IntrRay2Ray2 verification', () => {
    const tiq = new IntrRay2Ray2TI();
    const fiq = new IntrRay2Ray2FI();
    const llq = new IntrLine2Line2FI();

    it('TI and FI agree for transverse rays', () => {
        check(transverseRays, ({ ray0, ray1 }) => {
            const t = tiq.test(ray0, ray1);
            const f = fiq.find(ray0, ray1);
            expect(t.intersect).toBe(f.intersect);
            expect(t.numIntersections).toBe(f.numIntersections);
        });
    });

    it('a transverse hit lies on both rays at nonnegative parameters', () => {
        check(transverseRays, ({ ray0, ray1 }) => {
            const f = fiq.find(ray0, ray1);
            const ll = llq.find(
                Line.fromOriginDirection(ray0.origin, ray0.direction),
                Line.fromOriginDirection(ray1.origin, ray1.direction));
            expect(ll.numIntersections).toBe(1);
            const onBoth = ll.line0Parameter[0] >= 0
                && ll.line1Parameter[0] >= 0;
            expect(f.intersect).toBe(onBoth);
            if (!onBoth) {
                expect(f.numIntersections).toBe(0);
                return;
            }
            expect(f.numIntersections).toBe(1);
            expect(f.ray0Parameter[0]).toBeGreaterThanOrEqual(0);
            expect(f.ray1Parameter[0]).toBeGreaterThanOrEqual(0);
            expectVectorClose(f.point[0],
                add(ray0.origin, mul(f.ray0Parameter[0], ray0.direction)),
                0, 0);
            expectVectorClose(f.point[0],
                add(ray1.origin, mul(f.ray1Parameter[0], ray1.direction)),
                1e-9, 1e-12);
        });
    });

    it('the query is symmetric under argument swap', () => {
        check(transverseRays, ({ ray0, ray1 }) => {
            const a = fiq.find(ray0, ray1);
            const b = fiq.find(ray1, ray0);
            expect(a.intersect).toBe(b.intersect);
            expect(a.numIntersections).toBe(b.numIntersections);
            if (a.numIntersections === 1) {
                // Swapping negates both DotPerp values in each quotient, so
                // the parameters are bit-identical apart from the sign of a
                // zero; '+ 0' normalizes -0 to +0 for the Object.is compare.
                expect(a.ray0Parameter[0] + 0).toBe(b.ray1Parameter[0] + 0);
                expect(a.ray1Parameter[0] + 0).toBe(b.ray0Parameter[0] + 0);
            }
        });
    });

    it('collinear same-direction rays report a ray of intersections', () => {
        check(fc.tuple(latticeVector(2), axisDirection,
            fc.integer({ min: -5, max: 5 })), ([o, d, k]) => {
            const ray0 = Ray.fromOriginDirection(o, d);
            const ray1 = Ray.fromOriginDirection(add(o, mul(k, d)), d);
            const f = fiq.find(ray0, ray1);
            expect(tiq.test(ray0, ray1).numIntersections).toBe(INT32_MAX_V);
            expect(f.numIntersections).toBe(INT32_MAX_V);
            // t is the parameter of the ray1 origin along ray0.
            const t = k;
            if (t >= 0) {
                expect(f.ray0Parameter).toEqual([t, Number.MAX_VALUE]);
                expect(f.ray1Parameter).toEqual([0, Number.MAX_VALUE]);
                expectVectorClose(f.point[0], ray1.origin, 0, 0);
            } else {
                expect(f.ray0Parameter).toEqual([0, Number.MAX_VALUE]);
                expect(f.ray1Parameter).toEqual([-t, Number.MAX_VALUE]);
                expectVectorClose(f.point[0], ray0.origin, 0, 0);
            }
        });
    });

    it('collinear opposite rays overlap in the segment of the origins', () => {
        check(fc.tuple(latticeVector(2), axisDirection,
            fc.integer({ min: -5, max: 5 })), ([o, d, k]) => {
            const ray0 = Ray.fromOriginDirection(o, d);
            const ray1 = Ray.fromOriginDirection(add(o, mul(k, d)),
                mul(-1, d));
            const f = fiq.find(ray0, ray1);
            const t = tiq.test(ray0, ray1);
            const overlapT = k;   // Dot(D0, O1 - O0)
            if (overlapT > 0) {
                expect(f.numIntersections).toBe(2);
                expect(t.numIntersections).toBe(2);
                expect(f.ray0Parameter).toEqual([0, overlapT]);
                expect(f.ray1Parameter).toEqual([0, overlapT]);
                expectVectorClose(f.point[0], ray0.origin, 0, 0);
                expectVectorClose(f.point[1], ray1.origin, 0, 0);
            } else if (overlapT < 0) {
                expect(f.intersect).toBe(false);
                expect(t.intersect).toBe(false);
            } else {
                // Upstream quirk, preserved: at t == 0 the rays meet in the
                // single common origin, but TI reports numIntersections = 1
                // while FI reports 2 with a degenerate (zero-length) overlap
                // segment. Both agree on the boolean intersect.
                expect(t.intersect).toBe(true);
                expect(f.intersect).toBe(true);
                expect(t.numIntersections).toBe(1);
                expect(f.numIntersections).toBe(2);
                expect(f.ray0Parameter).toEqual([0, 0]);
                expectVectorClose(f.point[0], f.point[1], 0, 0);
            }
        });
    });

    it('the collinear result satisfies the documented reconstruction', () => {
        check(fc.tuple(latticeVector(2), axisDirection,
            fc.integer({ min: 1, max: 5 })), ([o, d, k]) => {
            const ray0 = Ray.fromOriginDirection(o, d);
            const ray1 = Ray.fromOriginDirection(add(o, mul(k, d)),
                mul(-1, d));
            const f = fiq.find(ray0, ray1);
            expect(f.numIntersections).toBe(2);
            expectVectorClose(f.point[1],
                add(ray0.origin, mul(f.ray0Parameter[1], ray0.direction)),
                0, 0);
            expectVectorClose(f.point[0],
                add(ray1.origin, mul(f.ray1Parameter[1], ray1.direction)),
                0, 0);
        });
    });

    it('parallel but distinct rays never intersect', () => {
        check(fc.tuple(wellScaledVector(2), unitVector(2),
            fc.double({ min: 0.5, max: 5, noNaN: true }), fc.boolean()),
            ([o, d, offset, flip]) => {
                const n = Vector.fromArray([-d.values[1], d.values[0]]);
                const ray0 = Ray.fromOriginDirection(o, d);
                const ray1 = Ray.fromOriginDirection(add(o, mul(offset, n)),
                    flip ? mul(-1, d) : d);
                expect(fiq.find(ray0, ray1).intersect).toBe(false);
                expect(tiq.test(ray0, ray1).intersect).toBe(false);
            });
    });
});
