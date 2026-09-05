import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, sub, mul, dot, length } from '../src/Vector.js';
import { IntrCapsule3Capsule3TI } from '../src/IntrCapsule3Capsule3.js';

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

// Brute-force containment: is the point inside the solid capsule?
function inCapsule(p: Vector, c: Capsule): boolean {
    const p0 = c.segment.p[0];
    const p1 = c.segment.p[1];
    const dir = sub(p1, p0);
    const dd = dot(dir, dir);
    let t = dd > 0 ? dot(sub(p, p0), dir) / dd : 0;
    t = Math.max(0, Math.min(1, t));
    const closest = add(p0, mul(t, dir));
    return length(sub(p, closest)) <= c.radius;
}

describe('IntrCapsule3Capsule3', () => {
    const ti = new IntrCapsule3Capsule3TI();

    it('reports intersection for overlapping parallel capsules', () => {
        const c0 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 1);
        const c1 = capsule(v3(-1, 1.5, 0), v3(1, 1.5, 0), 1);
        expect(ti.test(c0, c1).intersect).toBe(true);
    });

    it('reports no intersection for separated parallel capsules', () => {
        const c0 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 1);
        const c1 = capsule(v3(-1, 2.5, 0), v3(1, 2.5, 0), 1);
        expect(ti.test(c0, c1).intersect).toBe(false);
    });

    it('reports intersection when the capsules are exactly tangent', () => {
        // The segment distance is exactly 2 = r0 + r1, so the closed solids
        // touch and the query uses <=.
        const c0 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 1);
        const c1 = capsule(v3(-1, 2, 0), v3(1, 2, 0), 1);
        expect(ti.test(c0, c1).intersect).toBe(true);
    });

    it('handles perpendicular crossing capsules', () => {
        const c0 = capsule(v3(-2, 0, 0), v3(2, 0, 0), 0.25);
        const c1 = capsule(v3(0, -2, 0.4), v3(0, 2, 0.4), 0.25);
        // The axes pass within 0.4 of each other, sum of radii is 0.5.
        expect(ti.test(c0, c1).intersect).toBe(true);

        const c2 = capsule(v3(0, -2, 0.6), v3(0, 2, 0.6), 0.25);
        expect(ti.test(c0, c2).intersect).toBe(false);
    });

    it('is symmetric in its arguments', () => {
        const rnd = makeRandom(12345);
        for (let k = 0; k < 200; ++k) {
            const c0 = capsule(
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rnd());
            const c1 = capsule(
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rnd());
            expect(ti.test(c0, c1).intersect).toBe(ti.test(c1, c0).intersect);
        }
    });

    it('agrees with a brute-force sampling of the first capsule', () => {
        // Sampling the medial segment of capsule0 and offsetting by up to
        // radius0 produces points of capsule0. If any lies in capsule1, the
        // query must report an intersection.
        const rnd = makeRandom(777);
        let numFound = 0;
        for (let k = 0; k < 120; ++k) {
            const c0 = capsule(
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                0.3 + rnd() * 0.5);
            const c1 = capsule(
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                0.3 + rnd() * 0.5);
            const intersect = ti.test(c0, c1).intersect;

            let sampleHit = false;
            for (let s = 0; s <= 20 && !sampleHit; ++s) {
                const base = add(c0.segment.p[0],
                    mul(s / 20, sub(c0.segment.p[1], c0.segment.p[0])));
                for (let a = 0; a < 12 && !sampleHit; ++a) {
                    for (let b = 0; b < 6 && !sampleHit; ++b) {
                        const theta = (2 * Math.PI * a) / 12;
                        const phi = (Math.PI * b) / 5;
                        const off = mul(c0.radius, v3(
                            Math.sin(phi) * Math.cos(theta),
                            Math.sin(phi) * Math.sin(theta),
                            Math.cos(phi)));
                        if (inCapsule(add(base, off), c1)) {
                            sampleHit = true;
                        }
                    }
                }
            }

            if (sampleHit) {
                expect(intersect).toBe(true);
                ++numFound;
            }
        }
        expect(numFound).toBeGreaterThan(0);
    });

    it('handles degenerate zero-length segments and zero radii', () => {
        // Two points, treated as spheres of the given radii.
        const p0 = capsule(v3(0, 0, 0), v3(0, 0, 0), 1);
        const p1 = capsule(v3(1.5, 0, 0), v3(1.5, 0, 0), 1);
        expect(ti.test(p0, p1).intersect).toBe(true);

        const p2 = capsule(v3(2.5, 0, 0), v3(2.5, 0, 0), 1);
        expect(ti.test(p0, p2).intersect).toBe(false);

        // Zero radius: the capsules degenerate to their segments.
        const s0 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 0);
        const s1 = capsule(v3(0, -1, 0), v3(0, 1, 0), 0);
        expect(ti.test(s0, s1).intersect).toBe(true);

        const s2 = capsule(v3(0, -1, 1), v3(0, 1, 1), 0);
        expect(ti.test(s0, s2).intersect).toBe(false);
    });

    it('throws when a capsule is not 3-dimensional', () => {
        const c3 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 1);
        const c2 = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(Vector.fromArray([0, 0]),
                Vector.fromArray([1, 0])), 1);
        expect(() => ti.test(c3, c2)).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, seededRandom
} from './helpers/arbitraries.js';
import { DistSegmentSegment } from '../src/DistSegmentSegment.js';

const capsuleArb = fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
    fc.double({ min: 0.2, max: 4, noNaN: true }),
    fc.double({ min: 0.1, max: 2, noNaN: true }))
    .map(([p0, u, len, r]) => Capsule.fromSegmentRadius(
        Segment.fromEndpoints(p0, add(p0, mul(len, u))), r));

// Minimum distance between two segments, by dense sampling. Used only as a
// bracket (an upper bound on the true distance).
function sampledSegmentDistance(a: Segment, b: Segment, n: number): number {
    let best = Infinity;
    for (let i = 0; i <= n; ++i) {
        const p = add(a.p[0], mul(i / n, sub(a.p[1], a.p[0])));
        for (let j = 0; j <= n; ++j) {
            const q = add(b.p[0], mul(j / n, sub(b.p[1], b.p[0])));
            best = Math.min(best, length(sub(p, q)));
        }
    }
    return best;
}

describe('IntrCapsule3Capsule3 verification', () => {
    const tiq = new IntrCapsule3Capsule3TI();

    it('is symmetric under argument swap', () => {
        check(fc.tuple(capsuleArb, capsuleArb), ([c0, c1]) => {
            expect(tiq.test(c0, c1).intersect).toBe(tiq.test(c1, c0).intersect);
        });
    });

    it('agrees with a dense sampling of the two medial segments', () => {
        const rnd = seededRandom(0x71c4e0d3);
        for (let trial = 0; trial < 150; ++trial) {
            const mk = (): Capsule => {
                const p0 = Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3,
                    rnd() * 6 - 3]);
                const u = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                    rnd() * 2 - 1]);
                const p1 = add(p0, mul(0.5 + rnd() * 3, u));
                return Capsule.fromSegmentRadius(
                    Segment.fromEndpoints(p0, p1), 0.2 + rnd() * 1.2);
            };
            const c0 = mk(), c1 = mk();
            const rSum = c0.radius + c1.radius;
            const sampled = sampledSegmentDistance(c0.segment, c1.segment, 40);
            const got = tiq.test(c0, c1).intersect;
            // The sampled distance is an upper bound on the true distance, so
            // 'sampled <= rSum' implies an intersection; the converse needs a
            // margin for the sampling error.
            if (sampled <= rSum) {
                expect(got).toBe(true);
            }
            if (sampled > rSum + 0.5) {
                expect(got).toBe(false);
            }
        }
    }, 30000);

    it('a capsule always intersects itself', () => {
        check(capsuleArb, c => {
            expect(tiq.test(c, c).intersect).toBe(true);
        });
    });

    it('a capsule intersects a transverse capsule crossing its axis', () => {
        // A contained sub-capsule of the same medial segment is *not* used
        // here: that configuration is exactly collinear and hits the upstream
        // DistSegmentSegment parallelism defect pinned below.
        check(fc.tuple(capsuleArb, unitVector(3),
            fc.double({ min: 0.2, max: 0.8, noNaN: true }),
            fc.double({ min: 0.2, max: 2, noNaN: true }),
            fc.double({ min: 0.1, max: 1, noNaN: true })),
            ([c, d, s, len, r]) => {
                const e = sub(c.segment.p[1], c.segment.p[0]);
                const mid = add(c.segment.p[0], mul(s, e));
                if (Math.abs(dot(d, e)) > 0.9 * length(e)) {
                    return;   // nearly collinear; see the pinned defect
                }
                const other = Capsule.fromSegmentRadius(
                    Segment.fromEndpoints(sub(mid, mul(len, d)),
                        add(mid, mul(len, d))), r);
                expect(tiq.test(c, other).intersect).toBe(true);
                expect(tiq.test(other, c).intersect).toBe(true);
            });
    });

    it('inherits the upstream DistSegmentSegment parallelism defect', () => {
        // Upstream issue #418: DCPQuery<Segment,Segment>::operator() decides
        // parallelism with 'det = a*c - b*b > 0', which for mathematically
        // parallel directions can round one ulp above zero. The nonparallel
        // branch then runs on rounding noise. IntrCapsule3Capsule3 calls that
        // (non-robust) query, exactly as upstream does, so a capsule that
        // strictly contains another can be reported as non-intersecting.
        // The behaviour is preserved, not fixed, because the defect belongs
        // to DistSegmentSegment; ComputeRobust gets it right.
        const p0 = Vector.fromArray([0, 0, 0]);
        const p1 = Vector.fromArray([0, -2.8284271276953388,
            -2.828427121774217]);
        const outer = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), 0.1);
        const e = sub(p1, p0);
        const inner = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(add(p0, mul(0.10000000000000005, e)),
                add(p0, mul(0.8749998538754875, e))), 0.05);
        const q = new DistSegmentSegment();
        expect(q.compute(outer.segment, inner.segment).distance)
            .toBeGreaterThan(0.5);
        expect(q.computeRobust(outer.segment, inner.segment).distance)
            .toBe(0);
        // The consequence: a false negative for a contained capsule.
        expect(tiq.test(outer, inner).intersect).toBe(false);
    });

    it('capsules separated by more than the radius sum do not intersect', () => {
        check(fc.tuple(capsuleArb, unitVector(3),
            fc.double({ min: 0.5, max: 5, noNaN: true })),
            ([c, d, extra]) => {
                // Translate a copy far along d: the medial segments are then
                // at least (2*len + rSum + extra) apart in the worst case, so
                // shift by the diameter of the configuration.
                const len = length(sub(c.segment.p[1], c.segment.p[0]));
                const shift = mul(2 * len + 2 * c.radius + extra + 1, d);
                const far = Capsule.fromSegmentRadius(
                    Segment.fromEndpoints(add(c.segment.p[0], shift),
                        add(c.segment.p[1], shift)), c.radius);
                expect(tiq.test(c, far).intersect).toBe(false);
            });
    });

    it('is invariant under a common rigid motion', () => {
        check(fc.tuple(capsuleArb, capsuleArb, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([c0, c1, R, tr]) => {
            const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
            const xf = (c: Capsule): Capsule => Capsule.fromSegmentRadius(
                Segment.fromEndpoints(add(tr, rot(c.segment.p[0])),
                    add(tr, rot(c.segment.p[1]))), c.radius);
            const a = tiq.test(c0, c1).intersect;
            const b = tiq.test(xf(c0), xf(c1)).intersect;
            // Skip configurations right on the touching boundary, where the
            // decision legitimately flips under rounding.
            const rSum = c0.radius + c1.radius;
            const approx = sampledSegmentDistance(c0.segment, c1.segment, 20);
            if (Math.abs(approx - rSum) < 1e-6) {
                return;
            }
            expect(a).toBe(b);
        });
    });

    it('throws for capsules that are not 3-dimensional', () => {
        const c2 = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(Vector.fromArray([0, 0]),
                Vector.fromArray([1, 0])), 1);
        const c3 = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(Vector.fromArray([0, 0, 0]),
                Vector.fromArray([1, 0, 0])), 1);
        expect(() => tiq.test(c2, c3)).toThrow();
    });
});
