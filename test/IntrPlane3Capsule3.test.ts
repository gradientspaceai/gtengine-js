import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import { DistPointHyperplane } from '../src/DistPointHyperplane.js';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3Capsule3TI,
    defaultIntrPlane3Capsule3TIResult
} from '../src/IntrPlane3Capsule3.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import {
    check, fc, positive, rotationFrame, segment as arbSegment,
    plane as arbPlane, seededRandom, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

function capsule(p0: number[], p1: number[], radius: number): Capsule {
    return Capsule.fromSegmentRadius(
        Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1)),
        radius);
}

// An independent test: sample the capsule segment densely and check whether
// any sample is within 'radius' of the plane.
function bruteForceIntersect(P: Hyperplane, C: Capsule): boolean {
    const query = new DistPointHyperplane();
    const p0 = C.segment.p[0];
    const delta = sub(C.segment.p[1], p0);
    const n = 2000;
    for (let i = 0; i <= n; ++i) {
        const X = add(p0, mul(i / n, delta));
        if (query.compute(X, P).distance <= C.radius) {
            return true;
        }
    }
    return false;
}

const ti = new IntrPlane3Capsule3TI();

describe('IntrPlane3Capsule3', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrPlane3Capsule3TIResult().intersect).toBe(false);
    });

    it('detects endpoints on opposite sides of the plane', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        expect(ti.test(P, capsule([0, 0, -2], [0, 0, 3], 0.25)).intersect)
            .toBe(true);
    });

    it('detects an endpoint exactly on the plane', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        expect(ti.test(P, capsule([0, 0, 0], [0, 0, 3], 0.25)).intersect)
            .toBe(true);
    });

    it('detects a capsule that reaches the plane only with its end sphere', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        // Both endpoints are above the plane; the closer one is at z = 0.5.
        expect(ti.test(P, capsule([0, 0, 0.5], [0, 0, 3], 1)).intersect)
            .toBe(true);
        expect(ti.test(P, capsule([0, 0, 0.5], [0, 0, 3], 0.25)).intersect)
            .toBe(false);
    });

    it('reports the tangent configuration (distance equals radius)', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        expect(ti.test(P, capsule([0, 0, 1], [0, 0, 3], 1)).intersect)
            .toBe(true);
    });

    it('agrees with a dense sampling of the capsule on random inputs', () => {
        let state = 424242;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand(), rand(), rand()]);
            const C = capsule([rand() * 3, rand() * 3, rand() * 3],
                [rand() * 3, rand() * 3, rand() * 3],
                0.2 + Math.abs(rand()));
            const actual = ti.test(P, C).intersect;
            expect(actual).toBe(bruteForceIntersect(P, C));
            if (actual) {
                ++numHits;
            }
        }
        expect(numHits).toBeGreaterThan(20);
        expect(numHits).toBeLessThan(130);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrPlane3Capsule3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrPlane3Capsule3 verification', () => {
    const q = new IntrPlane3Capsule3TI();
    const vp = new DistPointHyperplane();

    const arbCapsule = fc.tuple(arbSegment(3), positive(4))
        .map(([s, r]) => Capsule.fromSegmentRadius(s, r));

    it('agrees with the exact solid-capsule criterion', () => {
        // For a unit-length plane normal the signed distance over the solid
        // capsule sweeps the interval
        //   [min(sd0,sd1) - r, max(sd0,sd1) + r],
        // so the capsule meets the plane exactly when that interval contains
        // zero. This is an independent derivation of upstream's two-branch
        // test.
        check(fc.tuple(arbPlane(3), arbCapsule), ([P, C]) => {
            const sd0 = vp.compute(C.segment.p[0], P).signedDistance;
            const sd1 = vp.compute(C.segment.p[1], P).signedDistance;
            const lo = Math.min(sd0, sd1) - C.radius;
            const hi = Math.max(sd0, sd1) + C.radius;
            expect(q.test(P, C).intersect).toBe(lo <= 0 && hi >= 0);
        });
    });

    it('reports an intersection whenever sampled capsule points straddle the plane', () => {
        // Sample the solid capsule (points within distance r of the segment)
        // and check that a sign change in Dot(N,X) - c forces 'intersect'.
        const rnd = seededRandom(0x5eed01);
        let numStraddle = 0;
        for (let iter = 0; iter < 300; ++iter) {
            const n = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (dot(n, n) < 1e-6) { continue; }
            normalize(n);
            const P = Hyperplane.fromNormalOrigin(n,
                Vector.fromArray([rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2]));
            const p0 = Vector.fromArray([rnd() * 4 - 2, rnd() * 4 - 2,
                rnd() * 4 - 2]);
            const p1 = Vector.fromArray([rnd() * 4 - 2, rnd() * 4 - 2,
                rnd() * 4 - 2]);
            const radius = 0.1 + rnd() * 1.5;
            const C = Capsule.fromSegmentRadius(Segment.fromEndpoints(p0, p1),
                radius);

            let sawPositive = false;
            let sawNegative = false;
            for (let i = 0; i <= 20; ++i) {
                const t = i / 20;
                const base = add(mul(1 - t, p0), mul(t, p1));
                for (let j = 0; j < 12; ++j) {
                    const d = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                        rnd() * 2 - 1]);
                    const len = Math.sqrt(dot(d, d));
                    if (len === 0) { continue; }
                    // Scale to a radius in [0, r] so the sample is inside the
                    // solid capsule.
                    const X = add(base, mul((radius * rnd()) / len, d));
                    const sd = dot(n, X) - P.constant;
                    if (sd > 1e-9) { sawPositive = true; }
                    if (sd < -1e-9) { sawNegative = true; }
                }
            }
            if (sawPositive && sawNegative) {
                ++numStraddle;
                expect(q.test(P, C).intersect).toBe(true);
            }
        }
        expect(numStraddle).toBeGreaterThan(50);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbPlane(3), arbCapsule, rotationFrame(3),
            wellScaledVector(3)),
            ([P, C, R, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(R[0], v), dot(R[1], v), dot(R[2], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const P2 = Hyperplane.fromNormalOrigin(rot(P.normal),
                    map(P.origin));
                const C2 = Capsule.fromSegmentRadius(
                    Segment.fromEndpoints(map(C.segment.p[0]),
                        map(C.segment.p[1])), C.radius);
                const sd0 = vp.compute(C.segment.p[0], P).signedDistance;
                const sd1 = vp.compute(C.segment.p[1], P).signedDistance;
                // Skip configurations within rounding distance of tangency,
                // where the rotated copy can land on the other side of the
                // '<=' comparison.
                const slack = Math.min(
                    Math.abs(Math.abs(sd0) - C.radius),
                    Math.abs(Math.abs(sd1) - C.radius),
                    Math.abs(sd0), Math.abs(sd1));
                if (slack < 1e-8) { return; }
                expect(q.test(P2, C2).intersect)
                    .toBe(q.test(P, C).intersect);
            });
    });

    it('treats a zero-radius capsule as its segment', () => {
        const P = Hyperplane.fromNormalOrigin(Vector.fromArray([0, 0, 1]),
            Vector.zero(3));
        const across = Capsule.fromSegmentRadius(Segment.fromEndpoints(
            Vector.fromArray([0, 0, -1]), Vector.fromArray([0, 0, 1])), 0);
        expect(q.test(P, across).intersect).toBe(true);
        const above = Capsule.fromSegmentRadius(Segment.fromEndpoints(
            Vector.fromArray([0, 0, 1]), Vector.fromArray([0, 0, 2])), 0);
        expect(q.test(P, above).intersect).toBe(false);
        const touching = Capsule.fromSegmentRadius(Segment.fromEndpoints(
            Vector.fromArray([0, 0, 0]), Vector.fromArray([0, 0, 2])), 0);
        expect(q.test(P, touching).intersect).toBe(true);
    });

    it('treats a degenerate (point) capsule as a sphere', () => {
        check(fc.tuple(arbPlane(3), wellScaledVector(3), positive(3)),
            ([P, c, r]) => {
                const C = Capsule.fromSegmentRadius(
                    Segment.fromEndpoints(c, c.clone()), r);
                const d = Math.abs(vp.compute(c, P).signedDistance);
                expect(q.test(P, C).intersect).toBe(d <= r);
            });
    });

    it('detects a capsule that touches the plane only through its end sphere', () => {
        // Both endpoints are strictly above the plane, but the lower endpoint
        // sphere reaches down to it. This is the second branch of the query.
        const P = Hyperplane.fromNormalOrigin(Vector.fromArray([0, 0, 1]),
            Vector.zero(3));
        const C = Capsule.fromSegmentRadius(Segment.fromEndpoints(
            Vector.fromArray([0, 0, 0.5]), Vector.fromArray([0, 0, 3])), 0.75);
        expect(q.test(P, C).intersect).toBe(true);
        const D = Capsule.fromSegmentRadius(Segment.fromEndpoints(
            Vector.fromArray([0, 0, 0.5]), Vector.fromArray([0, 0, 3])), 0.25);
        expect(q.test(P, D).intersect).toBe(false);
    });
});
