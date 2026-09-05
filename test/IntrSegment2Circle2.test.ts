import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    IntrSegment2Circle2TI,
    IntrSegment2Circle2FI,
    defaultIntrSegment2Circle2FIResult,
    intrSegment2Circle2FIDoQuery
} from '../src/IntrSegment2Circle2.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(center), radius);
}

function diskDepth(c: Hypersphere, P: Vector): number {
    const d = sub(P, c.center);
    return dot(d, d) - c.radius * c.radius;
}

describe('IntrSegment2Circle2', () => {
    const ti = new IntrSegment2Circle2TI();
    const fi = new IntrSegment2Circle2FI();

    it('finds both crossings of a long chord through a disk', () => {
        // The centered form has origin (0,0), direction (1,0), extent 5.
        const result = fi.find(segment([-5, 0], [5, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips against the segment endpoints', () => {
        // The segment stops halfway across the disk.
        const result = fi.find(segment([-5, 0], [0, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // Centered form: origin (-2.5,0), direction (1,0), extent 2.5.
        expect(result.parameter[0]).toBeCloseTo(1.5, 12);
        expect(result.parameter[1]).toBeCloseTo(2.5, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0, 12);
    });

    it('reports a segment fully inside the disk as no boundary crossing', () => {
        const s = segment([-0.2, 0], [0.2, 0]);
        const c = circle([0, 0], 1);
        const result = fi.find(s, c);
        // The line-disk t-interval is [-1,1] and the segment interval is
        // [-0.2,0.2]; the overlap is the whole segment.
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-0.2, 12);
        expect(result.parameter[1]).toBeCloseTo(0.2, 12);
        expect(ti.test(s, c).intersect).toBe(true);
    });

    it('rejects a segment that stops short of the disk', () => {
        const s = segment([-5, 0], [-2, 0]);
        const c = circle([0, 0], 1);
        expect(fi.find(s, c).intersect).toBe(false);
        expect(fi.find(s, c).numIntersections).toBe(0);
        expect(ti.test(s, c).intersect).toBe(false);
    });

    it('rejects a segment whose line misses the circle', () => {
        const s = segment([-5, 3], [5, 3]);
        const c = circle([0, 0], 1);
        expect(fi.find(s, c).intersect).toBe(false);
        expect(ti.test(s, c).intersect).toBe(false);
    });

    it('reports a single touching point', () => {
        // The segment endpoint just reaches the circle.
        const result = fi.find(segment([-5, 0], [-1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 9);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const result = defaultIntrSegment2Circle2FIResult();
        intrSegment2Circle2FIDoQuery(vec(0, 0), vec(1, 0), 5,
            circle([0, 0], 1), result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('agrees with brute-force sampling on random segments', () => {
        let seed = 1122334;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = circle([-0.75, 0.5], 1.25);
        for (let trial = 0; trial < 400; ++trial) {
            const s = segment([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 8 - 4, rand() * 8 - 4]);
            const result = fi.find(s, c);
            expect(ti.test(s, c).intersect).toBe(result.intersect);

            let sampledHit = false;
            for (let k = 0; k <= 400; ++k) {
                const t = k / 400;
                const P = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
                if (diskDepth(c, P) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            const { center, direction, extent } = s.getCenteredForm();
            for (let i = 0; i < result.numIntersections; ++i) {
                expect(Math.abs(result.parameter[i]))
                    .toBeLessThanOrEqual(extent + 1e-9);
                const P = add(center, mul(result.parameter[i], direction));
                expect(sub(P, result.point[i]).values[0]).toBeCloseTo(0, 9);
                expect(sub(P, result.point[i]).values[1]).toBeCloseTo(0, 9);
                expect(diskDepth(c, P)).toBeLessThan(1e-8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrSegment2Circle2.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { length } from '../src/Vector.js';
import { IntrLine2Circle2FI } from '../src/IntrLine2Circle2.js';
import { Line } from '../src/Line.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so a direction built from
// (cos a, sin a) never has a subnormal component; a subnormal component makes
// the exact DotPerp tests inside the line-line query underflow.
const angle2 = () => wellScaled(-Math.PI, Math.PI);

const segCircle = fc.tuple(
    wellScaled(-4, 4), wellScaled(-4, 4),
    wellScaled(-4, 4), wellScaled(-4, 4),
    wellScaled(-3, 3), wellScaled(-3, 3),
    fc.double({ min: 0.25, max: 3, noNaN: true, noDefaultInfinity: true })
).filter(([x0, y0, x1, y1]) =>
    (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) > 1e-4)
    .map(([x0, y0, x1, y1, cx, cy, r]) => ({
        segment: Segment.fromEndpoints(vec(x0, y0), vec(x1, y1)),
        circle: Hypersphere.fromCenterRadius(vec(cx, cy), r)
    }));

describe('IntrSegment2Circle2 verification', () => {
    const tiq = new IntrSegment2Circle2TI();
    const fiq = new IntrSegment2Circle2FI();
    const lcq = new IntrLine2Circle2FI();

    it('TI and FI always agree on intersect', () => {
        check(segCircle, ({ segment: s, circle: c }) => {
            expect(tiq.test(s, c).intersect).toBe(fiq.find(s, c).intersect);
        });
    });

    it('the segment hit is the line hit clipped to |t| <= extent', () => {
        check(segCircle, ({ segment: s, circle: c }) => {
            const cf = s.getCenteredForm();
            const line = Line.fromOriginDirection(cf.center, cf.direction);
            const lc = lcq.find(line, c);
            const f = fiq.find(s, c);
            if (!lc.intersect) {
                expect(f.intersect).toBe(false);
                expect(f.numIntersections).toBe(0);
                return;
            }
            const t0 = Math.max(lc.parameter[0], -cf.extent);
            const t1 = Math.min(lc.parameter[1], cf.extent);
            if (t0 > t1) {
                expect(f.intersect).toBe(false);
                return;
            }
            expect(f.intersect).toBe(true);
            expect(f.numIntersections).toBe(t0 < t1 ? 2 : 1);
            expect(f.parameter[0]).toBe(t0);
            if (f.numIntersections === 2) {
                expect(f.parameter[1]).toBe(t1);
            }
        });
    });

    it('reported points lie on the segment and in the closed disk', () => {
        check(segCircle, ({ segment: s, circle: c }) => {
            const cf = s.getCenteredForm();
            const f = fiq.find(s, c);
            for (let i = 0; i < f.numIntersections; ++i) {
                // The reported parameter is in the centered convention.
                expect(Math.abs(f.parameter[i]))
                    .toBeLessThanOrEqual(cf.extent);
                expectVectorClose(f.point[i],
                    add(cf.center, mul(f.parameter[i], cf.direction)), 0, 0);
                expect(diskDepth(c, f.point[i]))
                    .toBeLessThanOrEqual(1e-7 * (1 + c.radius * c.radius));
                // Equivalently, the point is at fraction u in [0,1] of the
                // endpoint form.
                const u = (f.parameter[i] + cf.extent) / (2 * cf.extent);
                expect(u).toBeGreaterThanOrEqual(-1e-12);
                expect(u).toBeLessThanOrEqual(1 + 1e-12);
                expectVectorClose(f.point[i],
                    add(s.p[0], mul(u, sub(s.p[1], s.p[0]))),
                    1e-9, 1e-9);
            }
        });
    });

    it('a segment strictly inside the disk is clipped to its own extent',
        () => {
            check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3),
                fc.double({ min: 1, max: 3, noNaN: true,
                    noDefaultInfinity: true }),
                angle2(), angle2(),
                fc.double({ min: 0, max: 0.4, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: 0, max: 0.4, noNaN: true,
                    noDefaultInfinity: true })),
            ([cx, cy, rad, a0, a1, u0, u1]) => {
                const c = Hypersphere.fromCenterRadius(vec(cx, cy), rad);
                const p0 = add(c.center,
                    vec(rad * u0 * Math.cos(a0), rad * u0 * Math.sin(a0)));
                const p1 = add(c.center,
                    vec(rad * u1 * Math.cos(a1), rad * u1 * Math.sin(a1)));
                if (length(sub(p1, p0)) < 1e-3) { return; }
                const s = Segment.fromEndpoints(p0, p1);
                const cf = s.getCenteredForm();
                const f = fiq.find(s, c);
                // Both endpoints are inside, so the clip is the whole
                // segment. This is the shape of the contained-primitive
                // defect reported in upstream issue #304 for the sphere and
                // ellipsoid TI queries; this query gets it right because it
                // clips an interval rather than testing endpoints.
                expect(tiq.test(s, c).intersect).toBe(true);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(2);
                expectClose(f.parameter[0], -cf.extent, 1e-12, 1e-12);
                expectClose(f.parameter[1], cf.extent, 1e-12, 1e-12);
                expectVectorClose(f.point[0], p0, 1e-9, 1e-9);
                expectVectorClose(f.point[1], p1, 1e-9, 1e-9);
            });
        });

    it('a fine sweep of the segment agrees with the reported interval', () => {
        const rnd = seededRandom(0x2fa41b7);
        for (let trial = 0; trial < 200; ++trial) {
            const c = Hypersphere.fromCenterRadius(
                vec(rnd() * 4 - 2, rnd() * 4 - 2), 0.5 + rnd() * 2);
            const s = Segment.fromEndpoints(
                vec(rnd() * 8 - 4, rnd() * 8 - 4),
                vec(rnd() * 8 - 4, rnd() * 8 - 4));
            if (length(sub(s.p[1], s.p[0])) < 1e-2) { continue; }
            const cf = s.getCenteredForm();
            const f = fiq.find(s, c);
            for (let k = 0; k <= 1000; ++k) {
                const t = -cf.extent + (2 * cf.extent * k) / 1000;
                const p = add(cf.center, mul(t, cf.direction));
                if (diskDepth(c, p) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(
                        f.parameter[f.numIntersections - 1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('the exported DoQuery reproduces the class result', () => {
        check(segCircle, ({ segment: s, circle: c }) => {
            const cf = s.getCenteredForm();
            const result = defaultIntrSegment2Circle2FIResult();
            intrSegment2Circle2FIDoQuery(cf.center, cf.direction, cf.extent, c,
                result);
            const f = fiq.find(s, c);
            expect(result.intersect).toBe(f.intersect);
            expect(result.numIntersections).toBe(f.numIntersections);
            expect(result.parameter).toEqual(f.parameter);
            expect(result.point[0].values).toEqual([0, 0]);
            expect(result.point[1].values).toEqual([0, 0]);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(segCircle, angle2(), wellScaled(-2, 2),
            wellScaled(-2, 2)),
        ([{ segment: s, circle: c }, ang, tx, ty]) => {
            const ca = Math.cos(ang), sa = Math.sin(ang);
            const rot = (v: Vector): Vector => vec(
                ca * v.get(0) - sa * v.get(1),
                sa * v.get(0) + ca * v.get(1));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty));
            const s2 = Segment.fromEndpoints(xf(s.p[0]), xf(s.p[1]));
            const c2 = Hypersphere.fromCenterRadius(xf(c.center), c.radius);
            const f0 = fiq.find(s, c);
            const f1 = fiq.find(s2, c2);
            if (f0.numIntersections !== f1.numIntersections) {
                return;   // grazing configuration
            }
            for (let i = 0; i < f0.numIntersections; ++i) {
                expectVectorClose(f1.point[i], xf(f0.point[i]), 1e-7, 1e-8);
            }
        });
    });

    it('a segment that misses the disk reports the default result', () => {
        const f = fiq.find(segment([-5, 3], [5, 3]), circle([0, 0], 1));
        expect(f.intersect).toBe(false);
        expect(f.numIntersections).toBe(0);
        expect(f.parameter).toEqual([0, 0]);
        expect(f.point[0].values).toEqual([0, 0]);
        expect(f.point[1].values).toEqual([0, 0]);
    });
});
