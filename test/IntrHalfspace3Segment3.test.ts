import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Segment } from '../src/Segment.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import {
    IntrHalfspace3Segment3TI,
    IntrHalfspace3Segment3FI
} from '../src/IntrHalfspace3Segment3.js';

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
