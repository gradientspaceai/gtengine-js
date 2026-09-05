import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Segment } from '../src/Segment.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import {
    IntrHalfspace3Segment3TI,
    IntrHalfspace3Segment3FI
} from '../src/IntrHalfspace3Segment3.js';
import { sub } from '../src/Vector.js';
import { check, fc } from './helpers/arbitraries.js';

function halfspace(nx: number, ny: number, nz: number, c: number): Halfspace {
    const n = Vector.fromArray([nx, ny, nz]);
    normalize(n);
    return Halfspace.fromNormalConstant(n, c);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrHalfspace3Segment3', () => {
    const ti = new IntrHalfspace3Segment3TI();
    const fi = new IntrHalfspace3Segment3FI();
    const zUp = halfspace(0, 0, 1, 0);  // z >= 0

    it('returns the original segment when it is inside the halfspace', () => {
        const s = segment([0, 0, 1], [1, 2, 3]);
        expect(ti.test(zUp, s).intersect).toBe(true);
        const result = fi.find(zUp, s);
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(2);
        expect(result.point[0].values).toEqual([0, 0, 1]);
        expect(result.point[1].values).toEqual([1, 2, 3]);
    });

    it('returns no intersection when the segment is strictly outside', () => {
        const s = segment([0, 0, -1], [1, 2, -3]);
        expect(ti.test(zUp, s).intersect).toBe(false);
        const result = fi.find(zUp, s);
        expect(result.intersect).toBe(false);
        expect(result.numPoints).toBe(0);
    });

    it('returns the crossing point for a straddling segment', () => {
        // (n,p,z) = (1,1,0). Upstream reports the clip point only; see the
        // port notes about the upstream inconsistency with its own table.
        const s = segment([0, 0, -1], [0, 0, 3]);
        expect(ti.test(zUp, s).intersect).toBe(true);
        const result = fi.find(zUp, s);
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(1);
        expect(result.point[0].values[2]).toBeCloseTo(0, 12);

        // The clip point interpolates correctly for an off-axis segment.
        const s2 = segment([0, 0, -1], [4, 8, 1]);
        const r2 = fi.find(zUp, s2);
        expect(r2.point[0].values[0]).toBeCloseTo(2, 12);
        expect(r2.point[0].values[1]).toBeCloseTo(4, 12);
        expect(r2.point[0].values[2]).toBeCloseTo(0, 12);
    });

    it('returns the touching endpoint when one endpoint is on the plane', () => {
        // (n,p,z) = (1,0,1): only the endpoint on the plane is in the closed
        // halfspace.
        const s0 = segment([1, 2, 0], [1, 2, -5]);
        const r0 = fi.find(zUp, s0);
        expect(r0.intersect).toBe(true);
        expect(r0.numPoints).toBe(1);
        expect(r0.point[0].values).toEqual([1, 2, 0]);

        const s1 = segment([1, 2, -5], [1, 2, 0]);
        const r1 = fi.find(zUp, s1);
        expect(r1.numPoints).toBe(1);
        expect(r1.point[0].values).toEqual([1, 2, 0]);
    });

    it('returns the original segment when it lies in the plane', () => {
        // (n,p,z) = (0,0,2).
        const s = segment([0, 0, 0], [1, 1, 0]);
        const result = fi.find(zUp, s);
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(2);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([1, 1, 0]);
    });

    it('handles a degenerate segment (both endpoints equal)', () => {
        expect(fi.find(zUp, segment([1, 1, 1], [1, 1, 1])).numPoints).toBe(2);
        expect(fi.find(zUp, segment([1, 1, -1], [1, 1, -1])).numPoints).toBe(0);
        expect(ti.test(zUp, segment([1, 1, -1], [1, 1, -1])).intersect).toBe(false);
    });

    it('keeps TI and FI consistent and returns points in the halfspace', () => {
        const rand = makeRandom(606060);
        let numHit = 0, numMiss = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const h = halfspace(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1,
                2 * rand() - 1);
            const s = segment(
                [4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2],
                [4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2]);

            const t = ti.test(h, s).intersect;
            const f = fi.find(h, s);
            expect(f.intersect).toBe(t);

            // Independent oracle: the maximum signed distance of the endpoints.
            const d0 = dot(h.normal, s.p[0]) - h.constant;
            const d1 = dot(h.normal, s.p[1]) - h.constant;
            expect(t).toBe(Math.max(d0, d1) >= 0);

            if (f.intersect) {
                ++numHit;
                for (let i = 0; i < f.numPoints; ++i) {
                    expect(dot(h.normal, f.point[i]) - h.constant)
                        .toBeGreaterThan(-1e-12);
                }
            } else {
                ++numMiss;
                expect(f.numPoints).toBe(0);
            }
        }
        expect(numHit).toBeGreaterThan(50);
        expect(numMiss).toBeGreaterThan(50);
    });
});

describe('IntrHalfspace3Segment3 verification', () => {
    const ti = new IntrHalfspace3Segment3TI();
    const fi = new IntrHalfspace3Segment3FI();

    // Integer normals and constants with integer endpoints. The queries only
    // use the signs of dot(N,P) - c and the ratio s0/(s0-s1), both invariant
    // under a positive scaling of (N,c), so an unnormalized integer normal
    // exercises exactly the same code paths while making the sign tests
    // exact. Small integers make the (n,p,z) case with an endpoint exactly on
    // the plane a common draw rather than a measure-zero accident.
    const latticeHs = fc.tuple(fc.integer({ min: -3, max: 3 }),
        fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 }),
        fc.integer({ min: -4, max: 4 }))
        .filter(([a, b, c]) => a !== 0 || b !== 0 || c !== 0)
        .map(([a, b, c, d]) =>
            Halfspace.fromNormalConstant(Vector.fromArray([a, b, c]), d));
    const latticeSeg = fc.tuple(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }))
        .map(([p0, p1]) => segment(p0, p1));

    it('TI and FI agree and match the exact (n,p,z) table', () => {
        check(fc.tuple(latticeHs, latticeSeg), ([h, s]) => {
            const f0 = dot(h.normal, s.p[0]) - h.constant;
            const f1 = dot(h.normal, s.p[1]) - h.constant;
            const numNeg = (f0 < 0 ? 1 : 0) + (f1 < 0 ? 1 : 0);
            const numPos = (f0 > 0 ? 1 : 0) + (f1 > 0 ? 1 : 0);
            const r = fi.find(h, s);
            expect(ti.test(h, s).intersect).toBe(Math.max(f0, f1) >= 0);
            expect(r.intersect).toBe(numNeg < 2);
            if (numNeg === 0) {
                expect(r.numPoints).toBe(2);
            } else if (numNeg === 1) {
                // Upstream reports one point in both sub-cases: the crossing
                // point when numPositive is 1, or the on-plane endpoint.
                expect(r.numPoints).toBe(1);
                expect(numPos + 1).toBeLessThanOrEqual(2);
            } else {
                expect(r.numPoints).toBe(0);
            }
        });
    });

    it('reported points lie on the segment and in the closed halfspace', () => {
        check(fc.tuple(latticeHs, latticeSeg), ([h, s]) => {
            const r = fi.find(h, s);
            const d = sub(s.p[1], s.p[0]);
            const dd = dot(d, d);
            for (let k = 0; k < r.numPoints; ++k) {
                const p = r.point[k];
                expect(Number.isFinite(p.values[0])).toBe(true);
                // The parameter of p along the segment must be in [0,1] and p
                // must reproduce (1-t)P0 + tP1. Everything here is a single
                // ratio of exact integers, so 1e-12 is generous.
                const t = dd > 0 ? dot(sub(p, s.p[0]), d) / dd : 0;
                expect(t).toBeGreaterThanOrEqual(-1e-12);
                expect(t).toBeLessThanOrEqual(1 + 1e-12);
                for (let i = 0; i < 3; ++i) {
                    expect(p.values[i]).toBeCloseTo(
                        s.p[0].values[i] + t * d.values[i], 10);
                }
                expect(dot(h.normal, p) - h.constant)
                    .toBeGreaterThanOrEqual(-1e-12 * (1 + Math.abs(h.constant)));
            }
        });
    });

    it('endpoints are copies, not aliases of the input segment', () => {
        const s = segment([0, 0, 1], [1, 2, 3]);
        const r = fi.find(halfspace(0, 0, 1, 0), s);
        expect(r.numPoints).toBe(2);
        expect(r.point[0]).not.toBe(s.p[0]);
        r.point[0].set(0, 99);
        expect(s.p[0].values[0]).toBe(0);
    });

    it('the clipped-segment case reports one point (upstream behavior)', () => {
        // (n,p,z) = (1,1,0). The intersection is the sub-segment from the
        // positive endpoint to the crossing point, but upstream sets
        // numPoints = 1 and stores only the crossing point. Preserved by the
        // port; see upstream issue #139.
        const r = fi.find(halfspace(0, 0, 1, 0), segment([0, 0, -1], [0, 0, 3]));
        expect(r.intersect).toBe(true);
        expect(r.numPoints).toBe(1);
        expect(r.point[0].values).toEqual([0, 0, 0]);
    });

    it('a degenerate (zero-length) segment behaves like its single point', () => {
        const h = halfspace(0, 0, 1, 0);
        const inside = fi.find(h, segment([1, 1, 2], [1, 1, 2]));
        expect(inside.numPoints).toBe(2);
        const on = fi.find(h, segment([1, 1, 0], [1, 1, 0]));
        expect(on.intersect).toBe(true);
        expect(on.numPoints).toBe(2);
        const outside = fi.find(h, segment([1, 1, -2], [1, 1, -2]));
        expect(outside.intersect).toBe(false);
        expect(outside.numPoints).toBe(0);
    });
});
