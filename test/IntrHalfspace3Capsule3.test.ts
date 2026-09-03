import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import { Halfspace } from '../src/Halfspace.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, sub, mul, dot, normalize } from '../src/Vector.js';
import { IntrHalfspace3Capsule3TI } from '../src/IntrHalfspace3Capsule3.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function capsule(p0: Vector, p1: Vector, radius: number): Capsule {
    return Capsule.fromSegmentRadius(Segment.fromEndpoints(p0, p1), radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrHalfspace3Capsule3', () => {
    const ti = new IntrHalfspace3Capsule3TI();

    // The halfspace z >= 0.
    const upper = Halfspace.fromNormalConstant(v3(0, 0, 1), 0);

    it('reports intersection when the capsule is inside the halfspace', () => {
        expect(ti.test(upper, capsule(v3(-1, 0, 3), v3(1, 0, 3), 1))
            .intersect).toBe(true);
    });

    it('reports no intersection when the capsule is strictly below', () => {
        expect(ti.test(upper, capsule(v3(-1, 0, -3), v3(1, 0, -3), 1))
            .intersect).toBe(false);
    });

    it('reports intersection when only the radius reaches the plane', () => {
        // The segment lies at z = -0.5 and the radius is 1, so the capsule
        // reaches z = 0.5 > 0.
        expect(ti.test(upper, capsule(v3(-1, 0, -0.5), v3(1, 0, -0.5), 1))
            .intersect).toBe(true);
    });

    it('reports intersection for exact tangency (closed halfspace)', () => {
        // The capsule touches the plane z = 0 at a single point.
        expect(ti.test(upper, capsule(v3(-1, 0, -1), v3(1, 0, -1), 1))
            .intersect).toBe(true);
    });

    it('uses the plane constant', () => {
        // The halfspace z >= 5.
        const high = Halfspace.fromNormalConstant(v3(0, 0, 1), 5);
        expect(ti.test(high, capsule(v3(0, 0, 3), v3(0, 0, 3.9), 1))
            .intersect).toBe(false);
        expect(ti.test(high, capsule(v3(0, 0, 3), v3(0, 0, 4.5), 1))
            .intersect).toBe(true);
    });

    it('handles a tilted halfspace normal', () => {
        const n = v3(1, 1, 1);
        normalize(n);
        const hs = Halfspace.fromNormalConstant(n, 0);
        // A capsule centered well on the negative side.
        expect(ti.test(hs, capsule(v3(-3, -3, -3), v3(-2, -2, -2), 0.5))
            .intersect).toBe(false);
        expect(ti.test(hs, capsule(v3(-3, -3, -3), v3(1, 1, 1), 0.5))
            .intersect).toBe(true);
    });

    it('agrees with a brute-force maximum over sampled capsule points', () => {
        const rnd = makeRandom(9182);
        for (let k = 0; k < 200; ++k) {
            const n = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(n, n) < 1e-6) {
                continue;
            }
            normalize(n);
            const hs = Halfspace.fromNormalConstant(n, rnd() * 4 - 2);
            const c = capsule(
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rnd() * 1.5);
            const intersect = ti.test(hs, c).intersect;

            // The exact maximum of Dot(N,X)-d over the capsule is
            // max(Dot(N,P0), Dot(N,P1)) - d + radius. Sample points of the
            // capsule and confirm that a positive sample implies intersect.
            let maxSample = -Infinity;
            for (let s = 0; s <= 12; ++s) {
                const base = add(c.segment.p[0],
                    mul(s / 12, sub(c.segment.p[1], c.segment.p[0])));
                // The maximizing surface offset is radius * N.
                const p = add(base, mul(c.radius, n));
                maxSample = Math.max(maxSample, dot(n, p) - hs.constant);
            }
            if (maxSample > 1e-12) {
                expect(intersect).toBe(true);
            }
            if (maxSample < -1e-12) {
                expect(intersect).toBe(false);
            }
        }
    });

    it('handles zero radius and a zero-length segment', () => {
        // A point capsule at the origin touches the closed halfspace z >= 0.
        expect(ti.test(upper, capsule(v3(0, 0, 0), v3(0, 0, 0), 0))
            .intersect).toBe(true);
        expect(ti.test(upper, capsule(v3(0, 0, -1e-9), v3(0, 0, -1e-9), 0))
            .intersect).toBe(false);
    });

    it('throws when the dimensions are not 3', () => {
        const hs2 = Halfspace.fromNormalConstant(Vector.fromArray([0, 1]), 0);
        const c = capsule(v3(0, 0, 0), v3(1, 0, 0), 1);
        expect(() => ti.test(hs2, c)).toThrow();
    });
});
