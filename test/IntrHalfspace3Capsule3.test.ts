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

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, seededRandom
} from './helpers/arbitraries.js';
import { length as vlength } from '../src/Vector.js';

const halfspace3 = fc.tuple(unitVector(3),
    fc.double({ min: -5, max: 5, noNaN: true }))
    .map(([n, c]) => Halfspace.fromNormalConstant(n, c));

const capsule3 = fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
    fc.double({ min: 0.2, max: 4, noNaN: true }),
    fc.double({ min: 0.1, max: 2, noNaN: true }))
    .map(([p0, u, len, r]) => Capsule.fromSegmentRadius(
        Segment.fromEndpoints(p0, add(p0, mul(len, u))), r));

describe('IntrHalfspace3Capsule3 verification', () => {
    const tiq = new IntrHalfspace3Capsule3TI();

    it('equals the support-point test Dot(N, S) >= constant', () => {
        check(fc.tuple(halfspace3, capsule3), ([h, c]) => {
            // The capsule point farthest along the normal is the endpoint
            // with the larger projection, pushed out by the radius.
            const e0 = dot(h.normal, c.segment.p[0]);
            const e1 = dot(h.normal, c.segment.p[1]);
            const support = Math.max(e0, e1) + c.radius;
            expect(tiq.test(h, c).intersect).toBe(support >= h.constant);
        });
    });

    it('agrees with a dense sampling of the capsule', () => {
        const rnd = seededRandom(0x9cd12f70);
        for (let trial = 0; trial < 150; ++trial) {
            const n = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (vlength(n) < 0.3) {
                continue;
            }
            normalize(n);
            const h = Halfspace.fromNormalConstant(n, rnd() * 6 - 3);
            const p0 = Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3,
                rnd() * 6 - 3]);
            const u = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            const c = Capsule.fromSegmentRadius(
                Segment.fromEndpoints(p0, add(p0, mul(0.5 + rnd() * 3, u))),
                0.2 + rnd() * 1.2);
            // Sample the capsule: points on the medial segment pushed out in
            // many directions, which covers the whole solid.
            let maxSigned = -Infinity;
            for (let i = 0; i <= 30; ++i) {
                const p = add(c.segment.p[0],
                    mul(i / 30, sub(c.segment.p[1], c.segment.p[0])));
                for (let k = 0; k < 60; ++k) {
                    const a = (2 * Math.PI * k) / 60;
                    for (const b of [-1, -0.5, 0, 0.5, 1]) {
                        const s = Math.sqrt(Math.max(0, 1 - b * b));
                        const dir = Vector.fromArray([s * Math.cos(a),
                            s * Math.sin(a), b]);
                        const q = add(p, mul(c.radius, dir));
                        maxSigned = Math.max(maxSigned,
                            dot(h.normal, q) - h.constant);
                    }
                }
            }
            const got = tiq.test(h, c).intersect;
            if (maxSigned >= 0) {
                expect(got).toBe(true);
            }
            if (maxSigned < -1e-6) {
                // The sampling is a lower bound on the true maximum, so a
                // negative sample does not prove separation; allow a margin.
                if (maxSigned < -0.2) {
                    expect(got).toBe(false);
                }
            }
        }
    }, 30000);

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(halfspace3, capsule3, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([h, c, R, tr]) => {
            const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
            const xf = (p: Vector): Vector => add(tr, rot(p));
            const n2 = rot(h.normal);
            // Dot(N', X') = Dot(N, X) + Dot(N', T), so the constant shifts.
            const h2 = Halfspace.fromNormalConstant(n2,
                h.constant + dot(n2, tr));
            const c2 = Capsule.fromSegmentRadius(
                Segment.fromEndpoints(xf(c.segment.p[0]), xf(c.segment.p[1])),
                c.radius);
            const e0 = dot(h.normal, c.segment.p[0]);
            const e1 = dot(h.normal, c.segment.p[1]);
            if (Math.abs(Math.max(e0, e1) + c.radius - h.constant) < 1e-9) {
                return;   // exactly tangent; the decision may flip
            }
            expect(tiq.test(h, c).intersect).toBe(tiq.test(h2, c2).intersect);
        });
    });

    it('a zero-radius capsule reduces to the segment test', () => {
        check(fc.tuple(halfspace3, wellScaledVector(3, -4, 4),
            wellScaledVector(3, -4, 4)), ([h, p0, p1]) => {
            const c = Capsule.fromSegmentRadius(
                Segment.fromEndpoints(p0, p1), 0);
            const inside = dot(h.normal, p0) >= h.constant
                || dot(h.normal, p1) >= h.constant;
            expect(tiq.test(h, c).intersect).toBe(inside);
        });
    });

    it('throws when the dimensions are not 3', () => {
        const h2 = Halfspace.fromNormalConstant(Vector.fromArray([1, 0]), 0);
        const c3 = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(Vector.fromArray([0, 0, 0]),
                Vector.fromArray([1, 0, 0])), 1);
        expect(() => tiq.test(h2, c3)).toThrow();
    });
});
