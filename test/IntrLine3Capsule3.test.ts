import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, normalize, sub, length } from '../src/Vector.js';
import {
    IntrLine3Capsule3TI,
    IntrLine3Capsule3FI
} from '../src/IntrLine3Capsule3.js';
import {
    intrLine3Capsule3FIDoQuery,
    defaultIntrLine3Capsule3FIResult
} from '../src/IntrLine3Capsule3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function line(p: number[], d: number[]): Line {
    const dir = vec(d);
    normalize(dir);
    return Line.fromOriginDirection(vec(p), dir);
}

function capsule(p0: number[], p1: number[], radius: number): Capsule {
    return Capsule.fromSegmentRadius(
        Segment.fromEndpoints(vec(p0), vec(p1)), radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Distance from a point to a segment, computed independently of the library.
function distanceToSegment(x: Vector, p0: Vector, p1: Vector): number {
    const d = sub(p1, p0);
    const dd = dot(d, d);
    let t = dd > 0 ? dot(sub(x, p0), d) / dd : 0;
    t = Math.min(1, Math.max(0, t));
    return length(sub(x, add(p0, mul(t, d))));
}

describe('IntrLine3Capsule3', () => {
    const ti = new IntrLine3Capsule3TI();
    const fi = new IntrLine3Capsule3FI();

    // Capsule with axis segment from (0,0,-2) to (0,0,2) and radius 1.
    const cap = capsule([0, 0, -2], [0, 0, 2], 1);

    it('finds the chord of a line crossing the capsule wall', () => {
        const l = line([0, 0, 0], [1, 0, 0]);
        expect(ti.test(l, cap).intersect).toBe(true);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a line parallel to the capsule axis', () => {
        const l = line([0.5, 0, -10], [0, 0, 1]);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const zOffset = Math.sqrt(1 - 0.25) + 2;
        expect(result.parameter[0]).toBeCloseTo(10 - zOffset, 12);
        expect(result.parameter[1]).toBeCloseTo(10 + zOffset, 12);
        // Both points are exactly on the capsule boundary.
        for (const p of result.point) {
            expect(distanceToSegment(p, vec([0, 0, -2]), vec([0, 0, 2])))
                .toBeCloseTo(1, 10);
        }
    });

    it('reports no intersection for a parallel line outside the radius', () => {
        const l = line([2, 0, -10], [0, 0, 1]);
        expect(ti.test(l, cap).intersect).toBe(false);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('finds the chord through a hemispherical cap only', () => {
        // The plane z = 2.5 cuts the top cap in a circle of radius
        // sqrt(1 - 0.25).
        const l = line([0, 0, 2.5], [1, 0, 0]);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const half = Math.sqrt(0.75);
        expect(result.parameter[0]).toBeCloseTo(-half, 10);
        expect(result.parameter[1]).toBeCloseTo(half, 10);
    });

    it('reports a single point for a line tangent to the capsule wall', () => {
        const l = line([1, -10, 0], [0, 1, 0]);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(10, 10);
        expect(result.parameter[1]).toBeCloseTo(result.parameter[0], 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 10);
        expect(result.point[0].values[1]).toBeCloseTo(0, 10);
    });

    it('misses a capsule entirely', () => {
        const l = line([5, 5, 0], [1, 0, 0]);
        expect(ti.test(l, cap).intersect).toBe(false);
        expect(fi.find(l, cap).intersect).toBe(false);
    });

    it('degrades on a zero-length capsule segment (upstream limitation)', () => {
        // A zero-length segment has no centered-form direction, so the
        // capsule coordinate frame built by the find query is degenerate. The
        // distance-based test query is still correct, but the find query
        // returns parameters that do not correspond to the sphere. The port
        // preserves upstream behavior; callers must supply a nondegenerate
        // capsule segment.
        const sphere = capsule([1, 2, 3], [1, 2, 3], 2);
        const l = line([1, 2, -10], [0, 0, 1]);
        expect(ti.test(l, sphere).intersect).toBe(true);

        const result = fi.find(l, sphere);
        // The true intersections would be at t = 11 and t = 15.
        expect(result.parameter[0]).not.toBeCloseTo(11, 6);
        expect(result.parameter[1]).not.toBeCloseTo(15, 6);
        for (const p of result.point) {
            expect(length(sub(p, vec([1, 2, 3])))).not.toBeCloseTo(2, 6);
        }
    });

    it('agrees with the distance-based test query and dense sampling', () => {
        const rand = makeRandom(8675309);
        const p0 = vec([-1, 0.5, 0]);
        const p1 = vec([2, -1, 1.5]);
        const cap2 = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), 0.8);
        for (let trial = 0; trial < 120; ++trial) {
            const l = line(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(l, cap2);
            expect(ti.test(l, cap2).intersect).toBe(result.intersect);

            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = -10 + (20 * k) / n;
                const x = add(l.origin, mul(t, l.direction));
                if (distanceToSegment(x, p0, p1) <= 0.8) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(3e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(3e-3);
            }

            if (result.intersect && result.numIntersections === 2) {
                // The endpoints of the chord lie on the capsule boundary.
                for (const p of result.point) {
                    expect(distanceToSegment(p, p0, p1)).toBeCloseTo(0.8, 8);
                }
            }
        }
    });
});

describe('intrLine3Capsule3FIDoQuery', () => {
    const c = capsule([0, 0, -1], [0, 0, 1], 1);

    it('matches the class query but does not compute points', () => {
        const l = line([-5, 0, 0], [1, 0, 0]);
        const result = defaultIntrLine3Capsule3FIResult();
        intrLine3Capsule3FIDoQuery(l.origin, l.direction, c, result);
        const expected = new IntrLine3Capsule3FI().find(l, c);
        expect(result.intersect).toBe(expected.intersect);
        expect(result.numIntersections).toBe(expected.numIntersections);
        expect(result.parameter[0]).toBeCloseTo(expected.parameter[0], 12);
        expect(result.parameter[1]).toBeCloseTo(expected.parameter[1], 12);
        // DoQuery leaves 'point' at its default value.
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([0, 0, 0]);
    });

    it('reports no intersection for a line missing the capsule', () => {
        const result = defaultIntrLine3Capsule3FIResult();
        intrLine3Capsule3FIDoQuery(vec([-5, 3, 0]), vec([1, 0, 0]), c,
            result);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, expectClose,
    expectVectorClose, seededRandom
} from './helpers/arbitraries.js';

// A capsule with a nondegenerate medial segment (the find query builds a
// coordinate frame from the segment direction) and a line with a unit-length
// direction.
const lineCapsule = fc.tuple(wellScaledVector(3, -6, 6), unitVector(3),
    wellScaledVector(3, -4, 4), unitVector(3),
    fc.double({ min: 0.5, max: 5, noNaN: true }),
    fc.double({ min: 0.2, max: 3, noNaN: true }))
    .map(([o, d, p0, u, len, radius]) => ({
        line: Line.fromOriginDirection(o, d),
        capsule: Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, add(p0, mul(len, u))), radius)
    }));

// Distance from a point to the capsule medial segment.
function distToSegment(seg: Segment, p: Vector): number {
    const e = sub(seg.p[1], seg.p[0]);
    const den = dot(e, e);
    let t = den > 0 ? dot(sub(p, seg.p[0]), e) / den : 0;
    t = Math.min(1, Math.max(0, t));
    return length(sub(p, add(seg.p[0], mul(t, e))));
}

describe('IntrLine3Capsule3 verification', () => {
    const tiq = new IntrLine3Capsule3TI();
    const fiq = new IntrLine3Capsule3FI();

    it('TI and FI agree away from tangency', () => {
        check(lineCapsule, ({ line: l, capsule: c }) => {
            const f = fiq.find(l, c);
            const t = tiq.test(l, c);
            if (f.intersect && f.parameter[1] - f.parameter[0] < 1e-6) {
                return;   // grazing; the two formulations may disagree
            }
            if (t.intersect !== f.intersect) {
                // Only accept a disagreement in the near-tangent band.
                const mid = f.intersect
                    ? 0.5 * (f.parameter[0] + f.parameter[1]) : 0;
                const p = add(l.origin, mul(mid, l.direction));
                expect(Math.abs(distToSegment(c.segment, p) - c.radius))
                    .toBeLessThan(1e-6);
                return;
            }
            expect(t.intersect).toBe(f.intersect);
        });
    });

    it('the reported points are on the line and on the capsule', () => {
        check(lineCapsule, ({ line: l, capsule: c }) => {
            const f = fiq.find(l, c);
            if (!f.intersect) {
                return;
            }
            expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
            // Upstream fills both entries whenever intersect is true.
            for (let i = 0; i < 2; ++i) {
                expectVectorClose(f.point[i],
                    add(l.origin, mul(f.parameter[i], l.direction)), 0, 0);
                // A boundary point is at distance exactly radius from the
                // medial segment. The parameters carry a square root, so the
                // residual grows near tangency; 1e-7 covers the sampled band.
                expectClose(distToSegment(c.segment, f.point[i]), c.radius,
                    1e-7, 1e-9);
                expect(Number.isNaN(f.parameter[i])).toBe(false);
            }
        });
    });

    it('a fine sweep of the line agrees with the reported interval', () => {
        const rnd = seededRandom(0x6ad3f21b);
        for (let trial = 0; trial < 150; ++trial) {
            const p0 = Vector.fromArray([rnd() * 4 - 2, rnd() * 4 - 2,
                rnd() * 4 - 2]);
            const u = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (length(u) < 0.3) {
                continue;
            }
            normalize(u);
            const seg = Segment.fromEndpoints(p0,
                add(p0, mul(0.5 + rnd() * 3, u)));
            const c = Capsule.fromSegmentRadius(seg, 0.3 + rnd() * 1.5);
            const d = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (length(d) < 0.3) {
                continue;
            }
            normalize(d);
            const l = Line.fromOriginDirection(
                Vector.fromArray([rnd() * 8 - 4, rnd() * 8 - 4,
                    rnd() * 8 - 4]), d);
            const f = fiq.find(l, c);
            for (let k = 0; k <= 600; ++k) {
                const t = -12 + (24 * k) / 600;
                const p = add(l.origin, mul(t, d));
                if (distToSegment(seg, p) < c.radius - 1e-5) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-6);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-6);
                }
            }
        }
    }, 30000);

    it('a line along the capsule axis spans the whole capsule', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
            fc.double({ min: 0.5, max: 5, noNaN: true }),
            fc.double({ min: 0.2, max: 3, noNaN: true }), fc.boolean()),
            ([p0, u, len, radius, flip]) => {
                const seg = Segment.fromEndpoints(p0, add(p0, mul(len, u)));
                const c = Capsule.fromSegmentRadius(seg, radius);
                const center = mul(0.5, add(seg.p[0], seg.p[1]));
                const dir = flip ? mul(-1, u) : u;
                const l = Line.fromOriginDirection(center, dir);
                const f = fiq.find(l, c);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(2);
                // The hemispherical-cap branch: |t| = extent + radius.
                const half = 0.5 * length(sub(seg.p[1], seg.p[0]));
                expectClose(f.parameter[0], -(half + radius), 1e-9, 1e-9);
                expectClose(f.parameter[1], half + radius, 1e-9, 1e-9);
            });
    });

    it('a line parallel to the axis but outside the wall misses', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
            fc.double({ min: 0.5, max: 5, noNaN: true }),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 1.05, max: 3, noNaN: true })),
            ([p0, R, len, radius, scale]) => {
                const seg = Segment.fromEndpoints(p0, add(p0, mul(len, R[0])));
                const c = Capsule.fromSegmentRadius(seg, radius);
                const off = add(p0, mul(scale * radius, R[1]));
                const l = Line.fromOriginDirection(off, R[0]);
                expect(fiq.find(l, c).intersect).toBe(false);
                expect(tiq.test(l, c).intersect).toBe(false);
            });
    });

    it('a line through a sampled interior point always hits', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
            fc.double({ min: 0.5, max: 5, noNaN: true }),
            fc.double({ min: 0.2, max: 3, noNaN: true }),
            fc.double({ min: 0, max: 1, noNaN: true }), unitVector(3),
            fc.double({ min: 0, max: 0.8, noNaN: true }), unitVector(3)),
            ([p0, u, len, radius, s, off, frac, d]) => {
                const seg = Segment.fromEndpoints(p0, add(p0, mul(len, u)));
                const c = Capsule.fromSegmentRadius(seg, radius);
                const base = add(seg.p[0], mul(s, sub(seg.p[1], seg.p[0])));
                const target = add(base, mul(frac * radius, off));
                const l = Line.fromOriginDirection(target, d);
                const f = fiq.find(l, c);
                expect(tiq.test(l, c).intersect).toBe(true);
                expect(f.intersect).toBe(true);
                expect(f.parameter[0]).toBeLessThanOrEqual(1e-9);
                expect(f.parameter[1]).toBeGreaterThanOrEqual(-1e-9);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineCapsule, rotationFrame(3), wellScaledVector(3,
            -4, 4)), ([{ line: l, capsule: c }, R, tr]) => {
            const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
            const xf = (p: Vector): Vector => add(tr, rot(p));
            const l2 = Line.fromOriginDirection(xf(l.origin), rot(l.direction));
            const c2 = Capsule.fromSegmentRadius(
                Segment.fromEndpoints(xf(c.segment.p[0]), xf(c.segment.p[1])),
                c.radius);
            const f1 = fiq.find(l, c);
            const f2 = fiq.find(l2, c2);
            const graze = (f: typeof f1): boolean => f.intersect
                && f.parameter[1] - f.parameter[0] < 1e-4;
            if (graze(f1) || graze(f2)) {
                return;
            }
            expect(f1.intersect).toBe(f2.intersect);
            if (!f1.intersect) {
                return;
            }
            expectClose(f1.parameter[0], f2.parameter[0], 1e-7, 1e-8);
            expectClose(f1.parameter[1], f2.parameter[1], 1e-7, 1e-8);
        });
    });

    it('the exported DoQuery reproduces the class result', () => {
        check(lineCapsule, ({ line: l, capsule: c }) => {
            const res = defaultIntrLine3Capsule3FIResult();
            intrLine3Capsule3FIDoQuery(l.origin, l.direction, c, res);
            const f = fiq.find(l, c);
            expect(res.intersect).toBe(f.intersect);
            expect(res.numIntersections).toBe(f.numIntersections);
            expect(res.parameter[0]).toBe(f.parameter[0]);
            expect(res.parameter[1]).toBe(f.parameter[1]);
            // The helper leaves point[] untouched; the class fills it.
            expect(res.point[0].values).toEqual([0, 0, 0]);
        });
    });
});
