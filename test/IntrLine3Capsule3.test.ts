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
        // 120 trials x 20001 samples: deterministic but slow enough to
        // exceed the 5 s default when the suite runs in parallel.
    }, 30000);
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
            // A line parallel to the medial segment goes through
            // DistLineSegment, whose parallel case is upstream issue #418:
            // for a collinear line it returns a nonzero distance, so TI
            // reports a miss where the root-based FI correctly reports a hit.
            const w = c.segment.getCenteredForm().direction;
            if (Math.abs(Math.abs(dot(w, l.direction)) - 1) < 1e-9) {
                return;
            }
            const f = fiq.find(l, c);
            const t = tiq.test(l, c);
            if (f.intersect && f.parameter[1] - f.parameter[0] < 1e-6) {
                return;   // grazing; the two formulations may disagree
            }
            if (t.intersect !== f.intersect) {
                // Only accept a disagreement in the near-tangent band, that
                // is, when the minimum distance from the line to the medial
                // segment is at the radius. That minimum is a convex function
                // of the line parameter, so a ternary search finds it; the
                // previous form probed the line origin instead, which says
                // nothing when the find query reports a miss.
                let lo = -1e4, hi = 1e4;
                const dist = (u: number): number => distToSegment(c.segment,
                    add(l.origin, mul(u, l.direction)));
                for (let k = 0; k < 300; ++k) {
                    const a = lo + (hi - lo) / 3, b = hi - (hi - lo) / 3;
                    if (dist(a) < dist(b)) { hi = b; } else { lo = a; }
                }
                expect(Math.abs(dist(0.5 * (lo + hi)) - c.radius))
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
                // The cap-junction plane (s = 0 or s = 1 with d perpendicular
                // to the axis) used to be excluded here: upstream collapsed
                // the interval to a point there (#461 item 3). The port fixes
                // it, so the case is now covered by this property and pinned
                // deterministically below.
                const l = Line.fromOriginDirection(target, d);
                const f = fiq.find(l, c);
                // A line parallel to the medial segment is excluded from the
                // TI comparison only: DistLineSegment's parallel case is
                // upstream issue #418 and reports a nonzero distance for a
                // collinear line.
                if (Math.abs(Math.abs(dot(u, d)) - 1) >= 1e-9) {
                    expect(tiq.test(l, c).intersect).toBe(true);
                }
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

// ---------------------------------------------------------------------------
// Regression tests for the two upstream defects fixed in the port
// (gtengine-js issue #461 items 2 and 3). Every test here fails on the
// unfixed algorithm.
// ---------------------------------------------------------------------------

// An oracle that is independent of the query: it works in world coordinates
// (the query builds its own capsule frame), solves the three boundary
// quadratics -- the infinite cylinder about the medial line and the two end
// spheres -- and keeps every root whose point lies in the solid capsule. The
// solid capsule is the union of the solid finite cylinder and the two solid
// end balls, so its boundary is contained in those three surfaces; the
// capsule is convex, so the line meets it in the interval spanned by the
// accepted roots.
function oracleCapsuleInterval(c: Capsule, l: Line): [number, number] | null {
    const C0 = c.segment.p[0];
    const A = sub(c.segment.p[1], C0);
    const W = mul(1 / length(A), A);
    const r = c.radius;
    const D = l.direction;
    const Q = sub(l.origin, C0);
    const perp = (v: Vector): Vector => sub(v, mul(dot(v, W), W));
    const roots: number[] = [];
    const addRoots = (a2: number, a1: number, a0: number): void => {
        if (a2 === 0) {
            if (a1 !== 0) { roots.push(-0.5 * a0 / a1); }
            return;
        }
        const discr = a1 * a1 - a0 * a2;
        if (discr < 0) { return; }
        const root = Math.sqrt(discr);
        roots.push((-a1 - root) / a2, (-a1 + root) / a2);
    };
    // The infinite cylinder about the medial line.
    const Qp = perp(Q), Dp = perp(D);
    addRoots(dot(Dp, Dp), dot(Qp, Dp), dot(Qp, Qp) - r * r);
    // The two end spheres (|D| = 1, so the leading coefficient is 1).
    for (const C of [C0, c.segment.p[1]]) {
        const E = sub(l.origin, C);
        addRoots(1, dot(E, D), dot(E, E) - r * r);
    }
    let lo = Number.POSITIVE_INFINITY, hi = Number.NEGATIVE_INFINITY;
    for (const t of roots) {
        const p = add(l.origin, mul(t, D));
        if (distToSegment(c.segment, p) <= r * (1 + 1e-9)) {
            if (t < lo) { lo = t; }
            if (t > hi) { hi = t; }
        }
    }
    return lo <= hi ? [lo, hi] : null;
}

describe('IntrLine3Capsule3 upstream-defect regressions', () => {
    const tiq = new IntrLine3Capsule3TI();
    const fiq = new IntrLine3Capsule3FI();

    it('reports a line tangent to a hemispherical cap (#461 item 2)', () => {
        // The line touches the bottom cap at its pole. Upstream accepts the
        // single root in the cap branch but never sets 'intersect' there, so
        // it reported intersect = false with numIntersections = 1 and no
        // points, and the ray and segment clips (guarded by 'intersect')
        // dropped the contact entirely.
        const c = capsule([0, 0, 0], [0, 0, -3], 2);
        const l = line([0, 0, -5], [1, 0, 0]);
        expect(tiq.test(l, c).intersect).toBe(true);
        const f = fiq.find(l, c);
        expect(f.intersect).toBe(true);
        expect(f.numIntersections).toBe(1);
        expect(f.parameter[0] + 0).toBe(0);
        expect(f.parameter[1] + 0).toBe(0);
        // Both points are filled and lie on the cap.
        for (let i = 0; i < 2; ++i) {
            expectVectorClose(f.point[i], vec([0, 0, -5]), 1e-12, 1e-12);
            expectClose(distToSegment(c.segment, f.point[i]), c.radius,
                1e-12, 1e-12);
        }
    });

    it('reports the whole chord in a cap-junction plane (#461 item 3)', () => {
        // The line lies in the plane z = +e that separates the cylinder wall
        // from the top hemisphere, through the medial endpoint. Upstream
        // accepted the same junction point twice -- once by the wall test
        // |z| <= e and once by the cap test z >= e -- and its early return
        // after two accepted roots discarded the true far endpoint, leaving a
        // degenerate interval.
        const s = Segment.fromEndpoints(vec([0, 0, 0]),
            vec([0, -4.999999999999982, -4.999999999999973]));
        const c = capsule([0, 2.9999999999999933, -2.9999999999999933],
            [0, 0, 0], 0.25);
        const cf = s.getCenteredForm();
        const l = Line.fromOriginDirection(cf.center, cf.direction);
        const f = fiq.find(l, c);
        expect(f.intersect).toBe(true);
        expect(f.numIntersections).toBe(2);
        // The chord through the centre of a junction disk has length 2*r.
        expectClose(f.parameter[1] - f.parameter[0], 2 * c.radius, 1e-12,
            1e-12);
        const oracle = oracleCapsuleInterval(c, l);
        expect(oracle).not.toBeNull();
        const o = oracle as [number, number];
        expectClose(f.parameter[0], o[0], 1e-12, 1e-12);
        expectClose(f.parameter[1], o[1], 1e-12, 1e-12);
    });

    it('sweeps lines lying exactly in a cap-junction plane', () => {
        // A capsule of half-length 1 and radius 1/2 in a frame tilted away
        // from the coordinate axes, so that the z-component of the line
        // direction in capsule coordinates is a rounding residual rather than
        // an exact zero: that is what splits the two junction-circle roots
        // between the wall test and the cap test. The junction disk has
        // radius 1/2, so a line in its plane at distance 'off' from the axis
        // cuts a chord of length 2*sqrt(r^2 - off^2).
        const W = vec([1, 2, 3]);
        normalize(W);
        const A = vec([0, -3, 2]);
        normalize(A);
        const B = vec([W.values[1] * A.values[2] - W.values[2] * A.values[1],
            W.values[2] * A.values[0] - W.values[0] * A.values[2],
            W.values[0] * A.values[1] - W.values[1] * A.values[0]]);
        normalize(B);
        const centre = vec([0.25, -0.5, 0.75]);
        const c = Capsule.fromSegmentRadius(Segment.fromEndpoints(
            sub(centre, W), add(centre, W)), 0.5);
        let exercised = 0;
        for (let k = 0; k < 64; ++k) {
            const th = (k / 64) * Math.PI;
            const d = add(mul(Math.cos(th), A), mul(Math.sin(th), B));
            normalize(d);
            const perp = add(mul(-Math.sin(th), A), mul(Math.cos(th), B));
            for (const off of [0, 0.1, 0.25, 0.4, 0.49]) {
                for (const end of [add(centre, W), sub(centre, W)]) {
                    const o = add(end, mul(off, perp));
                    const f = fiq.find(Line.fromOriginDirection(o, d), c);
                    expect(f.intersect).toBe(true);
                    expect(f.numIntersections).toBe(2);
                    expectClose(f.parameter[1] - f.parameter[0],
                        2 * Math.sqrt(0.25 - off * off), 1e-9, 1e-9);
                    ++exercised;
                }
            }
        }
        expect(exercised).toBe(640);
    });

    it('matches the independent world-frame oracle', () => {
        const rnd = seededRandom(0x51c0de3);
        let hits = 0;
        for (let trial = 0; trial < 4000; ++trial) {
            const p0 = Vector.fromArray([rnd() * 4 - 2, rnd() * 4 - 2,
                rnd() * 4 - 2]);
            const u = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (length(u) < 0.3) { continue; }
            normalize(u);
            const seg = Segment.fromEndpoints(p0,
                add(p0, mul(0.5 + rnd() * 3, u)));
            const c = Capsule.fromSegmentRadius(seg, 0.3 + rnd() * 1.5);
            const d = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (length(d) < 0.3) { continue; }
            normalize(d);
            // Aim a third of the lines at the cap-junction circles, where the
            // upstream double count lived: origin on a junction plane and
            // direction perpendicular to the axis.
            let l: Line;
            if (trial % 3 === 0) {
                const axis = sub(seg.p[1], seg.p[0]);
                const w = mul(1 / length(axis), axis);
                const dd = sub(d, mul(dot(d, w), w));
                if (length(dd) < 0.3) { continue; }
                normalize(dd);
                const end = trial % 6 === 0 ? seg.p[0] : seg.p[1];
                l = Line.fromOriginDirection(
                    add(end, mul((rnd() - 0.5) * c.radius, dd)), dd);
            }
            else {
                l = Line.fromOriginDirection(
                    Vector.fromArray([rnd() * 8 - 4, rnd() * 8 - 4,
                        rnd() * 8 - 4]), d);
            }
            const oracle = oracleCapsuleInterval(c, l);
            const f = fiq.find(l, c);
            if (oracle === null) {
                expect(f.intersect).toBe(false);
                continue;
            }
            if (oracle[1] - oracle[0] < 1e-6) {
                continue;   // grazing: entry and exit round together
            }
            ++hits;
            expect(f.intersect).toBe(true);
            expectClose(f.parameter[0], oracle[0], 1e-9, 1e-9);
            expectClose(f.parameter[1], oracle[1], 1e-9, 1e-9);
        }
        expect(hits).toBeGreaterThan(500);
    }, 30000);
});
