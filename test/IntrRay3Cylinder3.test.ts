import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import {
    IntrRay3Cylinder3FI,
    defaultIntrRay3Cylinder3FIResult,
    intrRay3Cylinder3FIDoQuery
} from '../src/IntrRay3Cylinder3.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(Vector.fromArray(origin), d), radius, height);
}

// Negative strictly inside the solid cylinder, zero on the boundary.
function cylinderSignedDepth(c: Cylinder3, P: Vector): number {
    const diff = sub(P, c.axis.origin);
    const z = dot(diff, c.axis.direction);
    const radial = sub(diff, mul(z, c.axis.direction));
    const rDepth = Math.sqrt(dot(radial, radial)) - c.radius;
    const zDepth = Math.abs(z) - 0.5 * c.height;
    return Math.max(rDepth, zDepth);
}

describe('IntrRay3Cylinder3FI', () => {
    const fi = new IntrRay3Cylinder3FI();

    it('finds both crossings of the cylinder wall', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(ray([-10, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(9, 9);
        expect(result.parameter[1]).toBeCloseTo(11, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 9);
        expect(result.point[1].values[0]).toBeCloseTo(1, 9);
    });

    it('finds the crossings of the end disks', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(ray([0, 0, -10], [0, 0, 1]), c);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(8, 9);
        expect(result.parameter[1]).toBeCloseTo(12, 9);
    });

    it('clips at the ray origin when the origin is inside', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(ray([0, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 9);
    });

    it('rejects a cylinder behind the ray', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(ray([10, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('rejects a ray whose line misses the cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        expect(fi.find(ray([-10, 5, 0], [1, 0, 0]), c).intersect).toBe(false);
        // Above the end disk.
        expect(fi.find(ray([-10, 0, 5], [1, 0, 0]), c).intersect).toBe(false);
    });

    it('rejects infinite cylinders', () => {
        const infinite = cylinder([0, 0, 0], [0, 0, 1], 1, -1);
        expect(() => fi.find(ray([-10, 0, 0], [1, 0, 0]), infinite)).toThrow();
    });

    it('exposes the DoQuery helper without computing points', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = defaultIntrRay3Cylinder3FIResult();
        intrRay3Cylinder3FIDoQuery(vec(-10, 0, 0), vec(1, 0, 0), c, result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(9, 9);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 2718281;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = cylinder([0.5, -0.5, 0.25], [0.3, 0.8, -0.5], 1.1, 2.5);
        for (let trial = 0; trial < 250; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, c);

            let sampledHit = false;
            for (let k = 0; k <= 1500; ++k) {
                const t = k * 0.01;
                if (cylinderSignedDepth(c,
                    add(r.origin, mul(t, r.direction))) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            if (result.intersect) {
                for (let i = 0; i < 2; ++i) {
                    expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-9);
                    const P = add(r.origin,
                        mul(result.parameter[i], r.direction));
                    expect(sub(P, result.point[i]).values[0])
                        .toBeCloseTo(0, 9);
                    expect(cylinderSignedDepth(c, P)).toBeLessThan(1e-7);
                }
                expect(result.parameter[0])
                    .toBeLessThanOrEqual(result.parameter[1] + 1e-12);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrRay3Cylinder3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { IntrLine3Cylinder3FI } from '../src/IntrLine3Cylinder3.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no direction component is a
// subnormal that would underflow when squared.
const angle3 = () => wellScaled(-Math.PI, Math.PI);

function unitDir3(th: number, ph: number): Vector {
    const d = vec(Math.cos(th) * Math.cos(ph), Math.sin(th) * Math.cos(ph),
        Math.sin(ph));
    normalize(d);
    return d;
}

const rayCylinder = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    angle3(), angle3(),
    wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2),
    angle3(), angle3(),
    fc.double({ min: 0.25, max: 2, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.5, max: 5, noNaN: true, noDefaultInfinity: true })
).map(([ox, oy, oz, th, ph, cx, cy, cz, ath, aph, r, h]) => ({
    ray: Ray.fromOriginDirection(vec(ox, oy, oz), unitDir3(th, ph)),
    cylinder: Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(vec(cx, cy, cz), unitDir3(ath, aph)), r, h)
}));

// cylinderSignedDepth is a maximum of convex functions, hence convex along
// the ray; a ternary search finds its minimum reliably.
function minRayCylDepth(c: Cylinder3, r: Ray, tMax = 40): number {
    let lo = 0, hi = tMax;
    for (let k = 0; k < 200; ++k) {
        const p = lo + (hi - lo) / 3, q = hi - (hi - lo) / 3;
        const fp = cylinderSignedDepth(c, add(r.origin, mul(p, r.direction)));
        const fq = cylinderSignedDepth(c, add(r.origin, mul(q, r.direction)));
        if (fp < fq) { hi = q; } else { lo = p; }
    }
    return cylinderSignedDepth(c,
        add(r.origin, mul(0.5 * (lo + hi), r.direction)));
}

describe('IntrRay3Cylinder3 verification', () => {
    const fiq = new IntrRay3Cylinder3FI();
    const lcq = new IntrLine3Cylinder3FI();

    it('the ray hit is the line hit clipped to t >= 0', () => {
        check(rayCylinder, ({ ray: r, cylinder: c }) => {
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

    it('reported points lie on the ray and on the cylinder boundary', () => {
        check(rayCylinder, ({ ray: r, cylinder: c }) => {
            const f = fiq.find(r, c);
            if (!f.intersect) {
                expect(f.numIntersections).toBe(0);
                expect(f.point[0].values).toEqual([0, 0, 0]);
                expect(f.point[1].values).toEqual([0, 0, 0]);
                return;
            }
            // Upstream fills both entries whenever intersect is true.
            for (let i = 0; i < 2; ++i) {
                expect(f.parameter[i]).toBeGreaterThanOrEqual(0);
                expect(Number.isFinite(f.parameter[i])).toBe(true);
                expectVectorClose(f.point[i],
                    add(r.origin, mul(f.parameter[i], r.direction)), 0, 0);
                const depth = cylinderSignedDepth(c, f.point[i]);
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
        const rnd = seededRandom(0x77af31b);
        for (let trial = 0; trial < 150; ++trial) {
            const axis = vec(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(axis, axis) < 1e-4) { continue; }
            normalize(axis);
            const c = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(
                    vec(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                    axis),
                0.3 + rnd() * 1.5, 0.5 + rnd() * 4);
            const d = vec(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d, d) < 1e-4) { continue; }
            normalize(d);
            const r = Ray.fromOriginDirection(
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5), d);
            const f = fiq.find(r, c);
            for (let k = 0; k <= 500; ++k) {
                const t = (18 * k) / 500;
                if (cylinderSignedDepth(c, add(r.origin, mul(t, d))) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('a ray whose origin is inside the cylinder starts at parameter 0',
        () => {
            check(fc.tuple(wellScaled(-2, 2), wellScaled(-2, 2),
                wellScaled(-2, 2), angle3(), angle3(),
                fc.double({ min: 0.25, max: 2, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: 0.5, max: 5, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: -0.8, max: 0.8, noNaN: true,
                    noDefaultInfinity: true }),
                angle3(), angle3()),
            ([cx, cy, cz, ath, aph, r, h, u, th, ph]) => {
                const axis = unitDir3(ath, aph);
                const c = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(vec(cx, cy, cz), axis), r, h);
                // A point on the axis, strictly between the end disks.
                const o = add(c.axis.origin, mul(u * 0.5 * h, axis));
                const f = fiq.find(
                    Ray.fromOriginDirection(o, unitDir3(th, ph)), c);
                expect(f.intersect).toBe(true);
                expect(f.parameter[0]).toBe(0);
                expect(f.numIntersections).toBe(2);
                expectClose(cylinderSignedDepth(c, f.point[1]), 0, 1e-7, 1e-8);
            });
        });

    it('the exported DoQuery reproduces the class result', () => {
        check(rayCylinder, ({ ray: r, cylinder: c }) => {
            const result = defaultIntrRay3Cylinder3FIResult();
            intrRay3Cylinder3FIDoQuery(r.origin, r.direction, c, result);
            const f = fiq.find(r, c);
            expect(result.intersect).toBe(f.intersect);
            expect(result.numIntersections).toBe(f.numIntersections);
            expect(result.parameter).toEqual(f.parameter);
            expect(result.point[0].values).toEqual([0, 0, 0]);
            expect(result.point[1].values).toEqual([0, 0, 0]);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(rayCylinder, angle3(), angle3(), angle3(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2)),
        ([{ ray: r, cylinder: c }, a1, a2, a3, tx, ty, tz]) => {
            if (Math.abs(minRayCylDepth(c, r)) < 1e-6) {
                return;   // grazing the wall or a rim
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
            const c2 = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(xf(c.axis.origin),
                    rot(c.axis.direction)), c.radius, c.height);
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

    it('rejects an infinite cylinder', () => {
        // Port deviation, documented in the source: upstream reads
        // cylinder.height directly and would compute nonsense for the
        // height = -1 infinite sentinel.
        const c = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(vec(0, 0, 0), vec(0, 0, 1)), 1, 4);
        c.makeInfiniteCylinder();
        expect(c.isFinite()).toBe(false);
        expect(() => fiq.find(ray([-5, 0, 0], [1, 0, 0]), c)).toThrow(
            'Infinite cylinders are not yet supported.');
    });
});
