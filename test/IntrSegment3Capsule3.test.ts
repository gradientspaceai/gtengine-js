import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import {
    IntrSegment3Capsule3TI,
    IntrSegment3Capsule3FI,
    defaultIntrSegment3Capsule3FIResult,
    intrSegment3Capsule3FIDoQuery
} from '../src/IntrSegment3Capsule3.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function capsule(p0: number[], p1: number[], radius: number): Capsule {
    return Capsule.fromSegmentRadius(segment(p0, p1), radius);
}

function capsuleSignedDepth(c: Capsule, P: Vector): number {
    const p0 = c.segment.p[0], p1 = c.segment.p[1];
    const e = sub(p1, p0);
    const d = sub(P, p0);
    const ee = dot(e, e);
    let t = ee > 0 ? dot(d, e) / ee : 0;
    t = Math.max(0, Math.min(1, t));
    const closest = add(p0, mul(t, e));
    const r = sub(P, closest);
    return Math.sqrt(dot(r, r)) - c.radius;
}

describe('IntrSegment3Capsule3', () => {
    const ti = new IntrSegment3Capsule3TI();
    const fi = new IntrSegment3Capsule3FI();

    it('finds both crossings of a capsule along the x-axis', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        // Centered form: origin (0,0,0), direction (1,0,0), extent 10.
        const result = fi.find(segment([-10, 0, 0], [10, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 9);
        expect(result.parameter[1]).toBeCloseTo(2, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-2, 9);
        expect(result.point[1].values[0]).toBeCloseTo(2, 9);
    });

    it('clips against the segment endpoints', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        // Centered form: origin (-5,0,0), direction (1,0,0), extent 5.
        const result = fi.find(segment([-10, 0, 0], [0, 0, 0]), c);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(3, 9);
        expect(result.parameter[1]).toBeCloseTo(5, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-2, 9);
        expect(result.point[1].values[0]).toBeCloseTo(0, 9);
    });

    it('rejects a segment that stops short of the capsule', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const s = segment([-10, 0, 0], [-5, 0, 0]);
        const result = fi.find(s, c);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(ti.test(s, c).intersect).toBe(false);
    });

    it('rejects a segment whose line misses the capsule', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const s = segment([-10, 5, 0], [10, 5, 0]);
        expect(fi.find(s, c).intersect).toBe(false);
        expect(ti.test(s, c).intersect).toBe(false);
    });

    it('reports the test query using the segment-segment distance', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        expect(ti.test(segment([-10, 1, 0], [10, 1, 0]), c).intersect)
            .toBe(true);
        expect(ti.test(segment([-10, 1.0001, 0], [10, 1.0001, 0]), c)
            .intersect).toBe(false);
        // A degenerate segment inside the capsule.
        expect(ti.test(segment([0, 0, 0], [0, 0, 0]), c).intersect).toBe(true);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = defaultIntrSegment3Capsule3FIResult();
        intrSegment3Capsule3FIDoQuery(vec(0, 0, 0), vec(1, 0, 0), 10, c,
            result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 9);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with brute-force sampling on random segments', () => {
        let seed = 31415926;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = capsule([-1, 0.5, 0.25], [1.5, -0.5, 1], 0.9);
        for (let trial = 0; trial < 250; ++trial) {
            const s = segment([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4]);
            const result = fi.find(s, c);

            let sampledHit = false;
            for (let k = 0; k <= 400; ++k) {
                const t = k / 400;
                const P = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
                if (capsuleSignedDepth(c, P) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
                expect(ti.test(s, c).intersect).toBe(true);
            }

            if (result.intersect) {
                const { center, direction, extent } = s.getCenteredForm();
                for (let i = 0; i < 2; ++i) {
                    expect(Math.abs(result.parameter[i]))
                        .toBeLessThanOrEqual(extent + 1e-9);
                    const P = add(center,
                        mul(result.parameter[i], direction));
                    expect(sub(P, result.point[i]).values[0])
                        .toBeCloseTo(0, 9);
                    expect(capsuleSignedDepth(c, P)).toBeLessThan(1e-7);
                }
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrSegment3Capsule3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { IntrLine3Capsule3FI } from '../src/IntrLine3Capsule3.js';
import { Line } from '../src/Line.js';
import { length } from '../src/Vector.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no frame component is a
// subnormal that would underflow when squared.
const angle3 = () => wellScaled(-Math.PI, Math.PI);

const segCapsule = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    fc.double({ min: 0.25, max: 2, noNaN: true, noDefaultInfinity: true })
).filter(([sx0, sy0, sz0, sx1, sy1, sz1, x0, y0, z0, x1, y1, z1]) =>
    (sx1 - sx0) ** 2 + (sy1 - sy0) ** 2 + (sz1 - sz0) ** 2 > 1e-2
    && (x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2 > 1e-2)
    .map(([sx0, sy0, sz0, sx1, sy1, sz1, x0, y0, z0, x1, y1, z1, r]) => ({
        segment: Segment.fromEndpoints(vec(sx0, sy0, sz0), vec(sx1, sy1, sz1)),
        capsule: Capsule.fromSegmentRadius(
            Segment.fromEndpoints(vec(x0, y0, z0), vec(x1, y1, z1)), r)
    }));

// The distance to a convex set is convex, so capsuleSignedDepth is convex
// along the segment; a ternary search finds its minimum reliably. A value
// near zero means the segment grazes the capsule, where the distance-based TI
// and the root-based FI round differently. This guard also keeps the property
// away from the DistSegmentSegment parallelism defect (upstream #418), which
// only bites when the two medial segments are (numerically) parallel and the
// true distance is at the radius.
function minSegCapsuleDepth(c: Capsule, s: Segment): number {
    const e = sub(s.p[1], s.p[0]);
    let lo = 0, hi = 1;
    for (let k = 0; k < 200; ++k) {
        const p = lo + (hi - lo) / 3, q = hi - (hi - lo) / 3;
        const fp = capsuleSignedDepth(c, add(s.p[0], mul(p, e)));
        const fq = capsuleSignedDepth(c, add(s.p[0], mul(q, e)));
        if (fp < fq) { hi = q; } else { lo = p; }
    }
    return capsuleSignedDepth(c, add(s.p[0], mul(0.5 * (lo + hi), e)));
}

// True when the line containing the segment lies in the capsule's
// cap-junction plane z = +-e while being perpendicular to the capsule axis.
// IntrLine3Capsule3 partitions the surface into the wall (|z| <= e) and the
// two hemispheres (z >= e, z <= -e); a line lying exactly in a junction plane
// can have one root accepted by the wall test and the other rejected by both,
// so the query returns a degenerate interval instead of the true chord. See
// the regression test below.
function inJunctionPlane(c: Capsule, s: Segment): boolean {
    const cf = c.segment.getCenteredForm();
    const w = cf.direction;
    const d = s.getCenteredForm().direction;
    const dz = Math.abs(dot(w, d));
    const zc = dot(w, sub(s.getCenteredForm().center, cf.center));
    return dz < 1e-6 && Math.abs(Math.abs(zc) - cf.extent) < 1e-6;
}

describe('IntrSegment3Capsule3 verification', () => {
    const tiq = new IntrSegment3Capsule3TI();
    const fiq = new IntrSegment3Capsule3FI();
    const lcq = new IntrLine3Capsule3FI();

    it('TI and FI agree away from grazing configurations', () => {
        check(segCapsule, ({ segment: s, capsule: c }) => {
            if (Math.abs(minSegCapsuleDepth(c, s)) < 1e-6
                || inJunctionPlane(c, s)) {
                return;
            }
            expect(tiq.test(s, c).intersect).toBe(fiq.find(s, c).intersect);
        });
    });

    it('pins the cap-junction-plane divergence of IntrLine3Capsule3', () => {
        // Upstream robustness defect (preserved): the segment starts exactly
        // at the capsule's medial endpoint and runs perpendicular to the
        // capsule axis, so the line containing it lies in the plane z = +e
        // that separates the cylinder wall from the top hemisphere. The true
        // intersection is a chord of length 2*radius, but IntrLine3Capsule3
        // accepts only one root (the wall test takes |z| <= e, the cap test
        // takes z >= e, and round-off in z puts the exit root on neither
        // side), so the line interval collapses to a point that the segment
        // clip then discards. TI, which is distance based, gets it right.
        const s = segment([0, 0, 0],
            [0, -4.999999999999982, -4.999999999999973]);
        const c = capsule([0, 2.9999999999999933, -2.9999999999999933],
            [0, 0, 0], 0.25);
        // The segment endpoint is 0.25 deep inside the solid capsule.
        expectClose(capsuleSignedDepth(c, s.p[0]), -0.25, 1e-12, 1e-12);
        expect(tiq.test(s, c).intersect).toBe(true);
        const f = fiq.find(s, c);
        expect(f.intersect).toBe(false);

        // The line query is where the collapse happens.
        const cf = s.getCenteredForm();
        const lc = lcq.find(Line.fromOriginDirection(cf.center, cf.direction),
            c);
        expect(lc.intersect).toBe(true);
        expect(lc.parameter[1] - lc.parameter[0]).toBeLessThan(1e-12);
    });

    it('the segment hit is the line hit clipped to |t| <= extent', () => {
        check(segCapsule, ({ segment: s, capsule: c }) => {
            const cf = s.getCenteredForm();
            const line = Line.fromOriginDirection(cf.center, cf.direction);
            const lc = lcq.find(line, c);
            const f = fiq.find(s, c);
            if (!lc.intersect) {
                expect(f.intersect).toBe(false);
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
            // Normalize the -0/+0 tie (toBe uses Object.is).
            expect(f.parameter[0] + 0).toBe(t0 + 0);
        });
    });

    it('reported points lie on the segment and on the capsule surface', () => {
        check(segCapsule, ({ segment: s, capsule: c }) => {
            const cf = s.getCenteredForm();
            const f = fiq.find(s, c);
            if (!f.intersect) {
                // Upstream quirk (preserved): a line tangent to a capsule
                // hemisphere leaves numIntersections = 1 with intersect
                // false. See test/IntrRay3Capsule3.test.ts for the pinned
                // reproduction.
                expect(f.numIntersections === 0 || f.numIntersections === 1)
                    .toBe(true);
                expect(f.point[0].values).toEqual([0, 0, 0]);
                expect(f.point[1].values).toEqual([0, 0, 0]);
                return;
            }
            for (let i = 0; i < 2; ++i) {
                // The parameter is in the centered convention, C + t*D with
                // |t| <= e.
                expect(Math.abs(f.parameter[i]))
                    .toBeLessThanOrEqual(cf.extent + 1e-12);
                expectVectorClose(f.point[i],
                    add(cf.center, mul(f.parameter[i], cf.direction)), 0, 0);
                const depth = capsuleSignedDepth(c, f.point[i]);
                if (Math.abs(Math.abs(f.parameter[i]) - cf.extent) > 1e-12) {
                    // Not a clipped endpoint, so the point is on the surface.
                    expectClose(depth, 0, 1e-7, 1e-8);
                }
                else {
                    expect(depth).toBeLessThanOrEqual(1e-7);
                }
            }
        });
    });

    it('a segment strictly inside the capsule is clipped to its endpoints',
        () => {
            check(fc.tuple(
                wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
                wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
                fc.double({ min: 0.25, max: 2, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: 0, max: 1, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: 0, max: 1, noNaN: true,
                    noDefaultInfinity: true })),
            ([x0, y0, z0, x1, y1, z1, r, u0, u1]) => {
                const a = vec(x0, y0, z0), b = vec(x1, y1, z1);
                if (length(sub(b, a)) < 1e-2) { return; }
                const c = Capsule.fromSegmentRadius(
                    Segment.fromEndpoints(a, b), r);
                // Two points on the medial segment: both strictly inside.
                const p0 = add(a, mul(u0, sub(b, a)));
                const p1 = add(a, mul(u1, sub(b, a)));
                if (length(sub(p1, p0)) < 1e-3) { return; }
                const s = Segment.fromEndpoints(p0, p1);
                const cf = s.getCenteredForm();
                // The contained-segment case that upstream issue #304 found
                // broken for the sphere and ellipsoid TI queries; this query
                // gets it right because its TI is distance based.
                expect(tiq.test(s, c).intersect).toBe(true);
                const f = fiq.find(s, c);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(2);
                expectClose(f.parameter[0], -cf.extent, 1e-12, 1e-12);
                expectClose(f.parameter[1], cf.extent, 1e-12, 1e-12);
                expectVectorClose(f.point[0], p0, 1e-8, 1e-9);
                expectVectorClose(f.point[1], p1, 1e-8, 1e-9);
            });
        });

    it('a fine sweep of the segment agrees with the reported interval', () => {
        const rnd = seededRandom(0x51ce07a);
        for (let trial = 0; trial < 150; ++trial) {
            const c = Capsule.fromSegmentRadius(Segment.fromEndpoints(
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2)),
                0.3 + rnd() * 1.5);
            if (length(sub(c.segment.p[1], c.segment.p[0])) < 1e-2) {
                continue;
            }
            const s = Segment.fromEndpoints(
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5),
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5));
            if (length(sub(s.p[1], s.p[0])) < 1e-2) { continue; }
            const cf = s.getCenteredForm();
            const f = fiq.find(s, c);
            for (let k = 0; k <= 500; ++k) {
                const t = -cf.extent + (2 * cf.extent * k) / 500;
                const p = add(cf.center, mul(t, cf.direction));
                if (capsuleSignedDepth(c, p) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('the exported DoQuery reproduces the class result', () => {
        check(segCapsule, ({ segment: s, capsule: c }) => {
            const cf = s.getCenteredForm();
            const result = defaultIntrSegment3Capsule3FIResult();
            intrSegment3Capsule3FIDoQuery(cf.center, cf.direction, cf.extent,
                c, result);
            const f = fiq.find(s, c);
            expect(result.intersect).toBe(f.intersect);
            expect(result.numIntersections).toBe(f.numIntersections);
            expect(result.parameter).toEqual(f.parameter);
            expect(result.point[0].values).toEqual([0, 0, 0]);
            expect(result.point[1].values).toEqual([0, 0, 0]);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(segCapsule, angle3(), angle3(), angle3(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2)),
        ([{ segment: s, capsule: c }, a1, a2, a3, tx, ty, tz]) => {
            if (Math.abs(minSegCapsuleDepth(c, s)) < 1e-6
                || inJunctionPlane(c, s)) {
                return;   // grazing (tangent) segment, or the junction plane
            }
            const ca = Math.cos(a1), sa = Math.sin(a1);
            const cb = Math.cos(a2), sb = Math.sin(a2);
            const cc = Math.cos(a3), sc = Math.sin(a3);
            const f0v = vec(ca * cb, sa * cb, -sb);
            const f1v = vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc,
                cb * sc);
            const f2v = vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc,
                cb * cc);
            const rot = (v: Vector): Vector => add(mul(v.get(0), f0v),
                add(mul(v.get(1), f1v), mul(v.get(2), f2v)));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty, tz));
            const s2 = Segment.fromEndpoints(xf(s.p[0]), xf(s.p[1]));
            const c2 = Capsule.fromSegmentRadius(Segment.fromEndpoints(
                xf(c.segment.p[0]), xf(c.segment.p[1])), c.radius);
            const g0 = fiq.find(s, c);
            const g1 = fiq.find(s2, c2);
            expect(g1.intersect).toBe(g0.intersect);
            if (!g0.intersect) { return; }
            expect(g1.numIntersections).toBe(g0.numIntersections);
            for (let i = 0; i < 2; ++i) {
                expectVectorClose(g1.point[i], xf(g0.point[i]), 1e-6, 1e-7);
            }
        });
    });

    it('a segment beyond the combined reach misses', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const s = segment([-5, 3, 0], [5, 3, 0]);
        expect(fiq.find(s, c).intersect).toBe(false);
        expect(tiq.test(s, c).intersect).toBe(false);
        const f = fiq.find(s, c);
        expect(f.numIntersections).toBe(0);
        expect(f.parameter).toEqual([0, 0]);
        expect(f.point[0].values).toEqual([0, 0, 0]);
        expect(f.point[1].values).toEqual([0, 0, 0]);
    });
});
