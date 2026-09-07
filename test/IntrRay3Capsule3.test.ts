import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import {
    IntrRay3Capsule3TI,
    IntrRay3Capsule3FI,
    defaultIntrRay3Capsule3FIResult,
    intrRay3Capsule3FIDoQuery
} from '../src/IntrRay3Capsule3.js';
import { Ray } from '../src/Ray.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function capsule(p0: number[], p1: number[], radius: number): Capsule {
    return Capsule.fromSegmentRadius(
        Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1)),
        radius);
}

// Distance from P to the capsule axis segment minus the radius: negative
// strictly inside the solid capsule.
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

describe('IntrRay3Capsule3', () => {
    const ti = new IntrRay3Capsule3TI();
    const fi = new IntrRay3Capsule3FI();

    it('finds both crossings of a capsule along the x-axis', () => {
        // A capsule with axis from (-1,0,0) to (1,0,0) and radius 1 spans
        // x in [-2,2] at y=z=0.
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = fi.find(ray([-10, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(8, 9);
        expect(result.parameter[1]).toBeCloseTo(12, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-2, 9);
        expect(result.point[1].values[0]).toBeCloseTo(2, 9);
    });

    it('crosses the cylindrical wall', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = fi.find(ray([0, -10, 0], [0, 1, 0]), c);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(9, 9);
        expect(result.parameter[1]).toBeCloseTo(11, 9);
    });

    it('clips at the ray origin when the origin is inside', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = fi.find(ray([0, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 9);
    });

    it('rejects a capsule behind the ray', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const r = ray([10, 0, 0], [1, 0, 0]);
        const result = fi.find(r, c);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(ti.test(r, c).intersect).toBe(false);
    });

    it('reports the test query using the ray-segment distance', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        // The ray passes at distance 1 from the axis: touching counts.
        expect(ti.test(ray([-10, 1, 0], [1, 0, 0]), c).intersect).toBe(true);
        expect(ti.test(ray([-10, 1.0001, 0], [1, 0, 0]), c).intersect)
            .toBe(false);
        // A ray pointing away from the capsule.
        expect(ti.test(ray([10, 0, 0], [1, 0, 0]), c).intersect).toBe(false);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = defaultIntrRay3Capsule3FIResult();
        intrRay3Capsule3FIDoQuery(vec(-10, 0, 0), vec(1, 0, 0), c, result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(8, 9);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 5150515;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = capsule([-1, 0.5, 0.25], [1.5, -0.5, 1], 0.9);
        for (let trial = 0; trial < 250; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, c);

            let sampledHit = false;
            for (let k = 0; k <= 1500; ++k) {
                const t = k * 0.01;
                if (capsuleSignedDepth(c,
                    add(r.origin, mul(t, r.direction))) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
                expect(ti.test(r, c).intersect).toBe(true);
            }

            if (result.intersect) {
                for (let i = 0; i < 2; ++i) {
                    expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-9);
                    const P = add(r.origin,
                        mul(result.parameter[i], r.direction));
                    expect(sub(P, result.point[i]).values[0])
                        .toBeCloseTo(0, 9);
                    expect(capsuleSignedDepth(c, P)).toBeLessThan(1e-7);
                }
                expect(result.parameter[0])
                    .toBeLessThanOrEqual(result.parameter[1] + 1e-12);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrRay3Capsule3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { IntrLine3Capsule3FI } from '../src/IntrLine3Capsule3.js';
import { Line } from '../src/Line.js';
import { length } from '../src/Vector.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no direction component is a
// subnormal that would underflow when squared.
const angle3 = () => wellScaled(-Math.PI, Math.PI);

function unitDir(th: number, ph: number): Vector {
    const d = vec(Math.cos(th) * Math.cos(ph), Math.sin(th) * Math.cos(ph),
        Math.sin(ph));
    normalize(d);
    return d;
}

const rayCapsule = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    angle3(), angle3(),
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    fc.double({ min: 0.25, max: 2, noNaN: true, noDefaultInfinity: true })
).filter(([, , , , , x0, y0, z0, x1, y1, z1]) =>
    (x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2 > 1e-2)
    .map(([ox, oy, oz, th, ph, x0, y0, z0, x1, y1, z1, r]) => ({
        ray: Ray.fromOriginDirection(vec(ox, oy, oz), unitDir(th, ph)),
        capsule: Capsule.fromSegmentRadius(
            Segment.fromEndpoints(vec(x0, y0, z0), vec(x1, y1, z1)), r)
    }));

// The distance to a convex set is convex, so capsuleSignedDepth is convex
// along the ray; a ternary search finds its minimum reliably. A value near
// zero means the ray grazes the capsule, where the distance-based TI and the
// root-based FI round differently.
function minRayCapsuleDepth(c: Capsule, r: Ray, tMax = 40): number {
    let lo = 0, hi = tMax;
    for (let k = 0; k < 200; ++k) {
        const p = lo + (hi - lo) / 3, q = hi - (hi - lo) / 3;
        const fp = capsuleSignedDepth(c, add(r.origin, mul(p, r.direction)));
        const fq = capsuleSignedDepth(c, add(r.origin, mul(q, r.direction)));
        if (fp < fq) { hi = q; } else { lo = p; }
    }
    return capsuleSignedDepth(c,
        add(r.origin, mul(0.5 * (lo + hi), r.direction)));
}

// True when the line containing the ray lies in the capsule's cap-junction
// plane z = +-e while being perpendicular to the capsule axis.
// IntrLine3Capsule3 partitions the surface into the wall (|z| <= e) and the
// two hemispheres (z >= e, z <= -e); a line lying exactly in a junction plane
// can have one root accepted by the wall test and the other rejected by both,
// so the query returns a degenerate interval instead of the true chord. See
// the regression test in test/IntrSegment3Capsule3.test.ts.
function inJunctionPlane(c: Capsule, r: Ray): boolean {
    const cf = c.segment.getCenteredForm();
    const w = cf.direction;
    const dz = Math.abs(dot(w, r.direction));
    const zc = dot(w, sub(r.origin, cf.center));
    return dz < 1e-6 && Math.abs(Math.abs(zc) - cf.extent) < 1e-6;
}

describe('IntrRay3Capsule3 verification', () => {
    const tiq = new IntrRay3Capsule3TI();
    const fiq = new IntrRay3Capsule3FI();
    const lcq = new IntrLine3Capsule3FI();

    it('TI and FI agree away from grazing configurations', () => {
        check(rayCapsule, ({ ray: r, capsule: c }) => {
            if (Math.abs(minRayCapsuleDepth(c, r)) < 1e-6
                || inJunctionPlane(c, r)) {
                return;
            }
            expect(tiq.test(r, c).intersect).toBe(fiq.find(r, c).intersect);
        });
    });

    it('the ray hit is the line hit clipped to t >= 0', () => {
        check(rayCapsule, ({ ray: r, capsule: c }) => {
            const line = Line.fromOriginDirection(r.origin, r.direction);
            const lc = lcq.find(line, c);
            const f = fiq.find(r, c);
            if (!lc.intersect) {
                expect(f.intersect).toBe(false);
                return;
            }
            const t0 = Math.max(lc.parameter[0], 0);
            const t1 = lc.parameter[1];
            if (t1 < 0) {
                expect(f.intersect).toBe(false);
                return;
            }
            expect(f.intersect).toBe(true);
            expect(f.numIntersections).toBe(t0 < t1 ? 2 : 1);
            // Normalize the -0/+0 tie (toBe uses Object.is).
            expect(f.parameter[0] + 0).toBe(t0 + 0);
        });
    });

    it('reported points lie on the ray and on the capsule surface', () => {
        check(rayCapsule, ({ ray: r, capsule: c }) => {
            const f = fiq.find(r, c);
            if (!f.intersect) {
                // Upstream quirk (preserved): a line tangent to a capsule
                // hemisphere leaves numIntersections = 1 with intersect
                // false, because only the wall-tangent branch of
                // IntrLine3Capsule3 sets intersect for a single root. See the
                // regression test below.
                expect(f.numIntersections === 0 || f.numIntersections === 1)
                    .toBe(true);
                expect(f.point[0].values).toEqual([0, 0, 0]);
                expect(f.point[1].values).toEqual([0, 0, 0]);
                return;
            }
            // Upstream fills both entries whenever intersect is true.
            for (let i = 0; i < 2; ++i) {
                expect(f.parameter[i]).toBeGreaterThanOrEqual(0);
                expectVectorClose(f.point[i],
                    add(r.origin, mul(f.parameter[i], r.direction)), 0, 0);
                // A clipped endpoint is the ray origin (inside the capsule);
                // an unclipped one is on the surface.
                const depth = capsuleSignedDepth(c, f.point[i]);
                if (f.parameter[i] > 0) {
                    expectClose(depth, 0, 1e-7, 1e-8);
                }
                else {
                    expect(depth).toBeLessThanOrEqual(1e-7);
                }
            }
        });
    });

    it('a fine sweep of the ray agrees with the reported interval', () => {
        const rnd = seededRandom(0x2d81b3c);
        for (let trial = 0; trial < 150; ++trial) {
            const c = Capsule.fromSegmentRadius(Segment.fromEndpoints(
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2)),
                0.3 + rnd() * 1.5);
            if (length(sub(c.segment.p[1], c.segment.p[0])) < 1e-2) {
                continue;
            }
            const d = vec(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d, d) < 1e-4) { continue; }
            normalize(d);
            const r = Ray.fromOriginDirection(
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5), d);
            const f = fiq.find(r, c);
            for (let k = 0; k <= 500; ++k) {
                const t = (18 * k) / 500;
                if (capsuleSignedDepth(c, add(r.origin, mul(t, d))) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('a ray whose origin is inside the capsule starts at parameter 0',
        () => {
            check(fc.tuple(
                wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
                wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
                fc.double({ min: 0.25, max: 2, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: 0, max: 1, noNaN: true,
                    noDefaultInfinity: true }),
                angle3(), angle3()),
            ([x0, y0, z0, x1, y1, z1, r, u, th, ph]) => {
                const p0 = vec(x0, y0, z0), p1 = vec(x1, y1, z1);
                if (length(sub(p1, p0)) < 1e-2) { return; }
                const c = Capsule.fromSegmentRadius(
                    Segment.fromEndpoints(p0, p1), r);
                // A point on the medial segment is strictly inside.
                const o = add(p0, mul(u, sub(p1, p0)));
                const f = fiq.find(
                    Ray.fromOriginDirection(o, unitDir(th, ph)), c);
                expect(f.intersect).toBe(true);
                expect(f.parameter[0]).toBe(0);
                expect(f.numIntersections).toBe(2);
                expectClose(capsuleSignedDepth(c, f.point[1]), 0, 1e-7, 1e-8);
            });
        });

    it('the exported DoQuery reproduces the class result', () => {
        check(rayCapsule, ({ ray: r, capsule: c }) => {
            const result = defaultIntrRay3Capsule3FIResult();
            intrRay3Capsule3FIDoQuery(r.origin, r.direction, c, result);
            const f = fiq.find(r, c);
            expect(result.intersect).toBe(f.intersect);
            expect(result.numIntersections).toBe(f.numIntersections);
            expect(result.parameter).toEqual(f.parameter);
            // DoQuery does not fill the points.
            expect(result.point[0].values).toEqual([0, 0, 0]);
            expect(result.point[1].values).toEqual([0, 0, 0]);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(rayCapsule, angle3(), angle3(), angle3(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2)),
        ([{ ray: r, capsule: c }, a1, a2, a3, tx, ty, tz]) => {
            if (Math.abs(minRayCapsuleDepth(c, r)) < 1e-6
                || inJunctionPlane(c, r)) {
                return;   // grazing (tangent) ray, or the junction plane
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
            const r2 = Ray.fromOriginDirection(xf(r.origin), rot(r.direction));
            const c2 = Capsule.fromSegmentRadius(Segment.fromEndpoints(
                xf(c.segment.p[0]), xf(c.segment.p[1])), c.radius);
            const g0 = fiq.find(r, c);
            const g1 = fiq.find(r2, c2);
            expect(g1.intersect).toBe(g0.intersect);
            if (!g0.intersect) { return; }
            expect(g1.numIntersections).toBe(g0.numIntersections);
            for (let i = 0; i < 2; ++i) {
                expectClose(g1.parameter[i], g0.parameter[i], 1e-6, 1e-7);
                expectVectorClose(g1.point[i], xf(g0.point[i]), 1e-6, 1e-7);
            }
        });
    });

    it('pins the hemisphere-tangency quirk of IntrLine3Capsule3', () => {
        // Upstream defect (preserved): IntrLine3Capsule3 sets intersect =
        // true for a single root only in the wall-tangent branch. When the
        // line is tangent to one of the hemispherical caps, the cap branch
        // records the root and the trailing block duplicates it, but
        // 'intersect' is never set. The ray and segment queries then skip
        // their interval clip entirely (it is guarded by 'intersect') and
        // report intersect = false with numIntersections = 1.
        const c = capsule([0, 0, 0], [0, 0, -3], 2);
        const r = ray([0, 0, -5], [1, 0, 0]);
        // The ray origin is exactly on the cap surface.
        expect(capsuleSignedDepth(c, r.origin)).toBe(0);
        // The distance-based TI reports the contact.
        expect(tiq.test(r, c).intersect).toBe(true);
        const f = fiq.find(r, c);
        expect(f.intersect).toBe(false);
        expect(f.numIntersections).toBe(1);
        // (+0 normalizes the -0/+0 tie that toEqual would reject.)
        expect(f.parameter.map(t => t + 0)).toEqual([0, 0]);
        // No point is filled, because operator() guards on 'intersect'.
        expect(f.point[0].values).toEqual([0, 0, 0]);
        expect(f.point[1].values).toEqual([0, 0, 0]);

        // The line query is where the flag is dropped.
        const lc = lcq.find(
            Line.fromOriginDirection(r.origin, r.direction), c);
        expect(lc.intersect).toBe(false);
        expect(lc.numIntersections).toBe(1);
    });

    it('a ray parallel to the axis and outside the wall misses', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        expect(fiq.find(ray([-5, 2.5, 0], [1, 0, 0]), c).intersect)
            .toBe(false);
        expect(tiq.test(ray([-5, 2.5, 0], [1, 0, 0]), c).intersect)
            .toBe(false);
        // Along the axis the capsule spans [-2, 2].
        const f = fiq.find(ray([-5, 0, 0], [1, 0, 0]), c);
        expect(f.intersect).toBe(true);
        expectClose(f.parameter[0], 3, 1e-12, 1e-12);
        expectClose(f.parameter[1], 7, 1e-12, 1e-12);
    });
});
