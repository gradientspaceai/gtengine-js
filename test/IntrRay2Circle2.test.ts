import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    IntrRay2Circle2TI,
    IntrRay2Circle2FI,
    defaultIntrRay2Circle2FIResult,
    intrRay2Circle2FIDoQuery
} from '../src/IntrRay2Circle2.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(center), radius);
}

// Squared distance from the disk center minus r^2: negative inside.
function diskDepth(c: Hypersphere, P: Vector): number {
    const d = sub(P, c.center);
    return dot(d, d) - c.radius * c.radius;
}

describe('IntrRay2Circle2', () => {
    const ti = new IntrRay2Circle2TI();
    const fi = new IntrRay2Circle2FI();

    it('finds both crossings when the ray starts outside', () => {
        const result = fi.find(ray([-5, 0], [1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips the near crossing when the ray starts inside the disk', () => {
        const result = fi.find(ray([0, 0], [1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
    });

    it('rejects a disk behind the ray origin', () => {
        const r = ray([5, 0], [1, 0]);
        const c = circle([0, 0], 1);
        expect(fi.find(r, c).intersect).toBe(false);
        expect(fi.find(r, c).numIntersections).toBe(0);
        expect(ti.test(r, c).intersect).toBe(false);
    });

    it('reports a single touching point when the ray origin is on the circle', () => {
        const result = fi.find(ray([1, 0], [1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
    });

    it('reports the tangent point', () => {
        const result = fi.find(ray([-5, 1], [1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(5, 9);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('misses a disk that the line also misses', () => {
        const r = ray([-5, 3], [1, 0]);
        const c = circle([0, 0], 1);
        expect(fi.find(r, c).intersect).toBe(false);
        expect(ti.test(r, c).intersect).toBe(false);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const result = defaultIntrRay2Circle2FIResult();
        intrRay2Circle2FIDoQuery(vec(-5, 0), vec(1, 0), circle([0, 0], 1),
            result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 90210;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = circle([0.5, -1], 1.75);
        for (let trial = 0; trial < 400; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, c);
            expect(ti.test(r, c).intersect).toBe(result.intersect);

            // Sample the ray to detect whether it enters the disk.
            let sampledHit = false;
            for (let k = 0; k <= 4000; ++k) {
                const t = k * 0.005;
                if (diskDepth(c, add(r.origin, mul(t, r.direction))) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            for (let i = 0; i < result.numIntersections; ++i) {
                expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-12);
                const P = add(r.origin, mul(result.parameter[i], r.direction));
                expect(sub(P, result.point[i]).values[0]).toBeCloseTo(0, 9);
                expect(sub(P, result.point[i]).values[1]).toBeCloseTo(0, 9);
                expect(diskDepth(c, P)).toBeLessThan(1e-8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrRay2Circle2.h.
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

const rayCircle = fc.tuple(
    wellScaled(-4, 4), wellScaled(-4, 4), angle2(),
    wellScaled(-3, 3), wellScaled(-3, 3),
    fc.double({ min: 0.25, max: 3, noNaN: true, noDefaultInfinity: true })
).map(([ox, oy, a, cx, cy, r]) => ({
    ray: Ray.fromOriginDirection(vec(ox, oy),
        vec(Math.cos(a), Math.sin(a))),
    circle: Hypersphere.fromCenterRadius(vec(cx, cy), r)
}));

// The discriminant of the line-circle quadratic; near zero the line grazes
// the circle and the discrete answers are ambiguous.
function rayDiscriminant(r: Ray, c: Hypersphere): number {
    const diff = sub(r.origin, c.center);
    const a1 = dot(r.direction, diff);
    const a0 = dot(diff, diff) - c.radius * c.radius;
    return a1 * a1 - a0;
}

describe('IntrRay2Circle2 verification', () => {
    const tiq = new IntrRay2Circle2TI();
    const fiq = new IntrRay2Circle2FI();
    const lcq = new IntrLine2Circle2FI();

    it('TI and FI always agree on intersect', () => {
        check(rayCircle, ({ ray: r, circle: c }) => {
            expect(tiq.test(r, c).intersect).toBe(fiq.find(r, c).intersect);
        });
    });

    it('the ray hit is the line hit clipped to t >= 0', () => {
        check(rayCircle, ({ ray: r, circle: c }) => {
            if (Math.abs(rayDiscriminant(r, c)) < 1e-9) {
                return;   // grazing; the clip endpoints are ambiguous
            }
            const line = Line.fromOriginDirection(r.origin, r.direction);
            const lc = lcq.find(line, c);
            const f = fiq.find(r, c);
            if (!lc.intersect) {
                expect(f.intersect).toBe(false);
                expect(f.numIntersections).toBe(0);
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
            expect(f.parameter[0]).toBe(t0);
            if (f.numIntersections === 2) {
                expect(f.parameter[1]).toBe(t1);
            }
        });
    });

    it('reported points lie on the ray and in the closed disk', () => {
        check(rayCircle, ({ ray: r, circle: c }) => {
            const f = fiq.find(r, c);
            for (let i = 0; i < f.numIntersections; ++i) {
                expect(f.parameter[i]).toBeGreaterThanOrEqual(0);
                expectVectorClose(f.point[i],
                    add(r.origin, mul(f.parameter[i], r.direction)), 0, 0);
                // Inside the closed disk, up to the conditioning of the root.
                expect(diskDepth(c, f.point[i]))
                    .toBeLessThanOrEqual(1e-7 * (1 + c.radius * c.radius));
            }
            for (let i = 0; i < f.numIntersections; ++i) {
                expect(Number.isFinite(f.parameter[i])).toBe(true);
            }
            if (f.numIntersections === 2) {
                expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
            }
        });
    });

    it('a fine sweep of the ray agrees with the reported interval', () => {
        const rnd = seededRandom(0x1177aa31);
        for (let trial = 0; trial < 200; ++trial) {
            const c = Hypersphere.fromCenterRadius(
                vec(rnd() * 4 - 2, rnd() * 4 - 2), 0.5 + rnd() * 2);
            const a = rnd() * 2 * Math.PI;
            const d = vec(Math.cos(a), Math.sin(a));
            const r = Ray.fromOriginDirection(
                vec(rnd() * 6 - 3, rnd() * 6 - 3), d);
            const f = fiq.find(r, c);
            for (let k = 0; k <= 1200; ++k) {
                const t = (12 * k) / 1200;
                const p = add(r.origin, mul(t, d));
                // Strictly inside the disk implies the parameter is inside
                // the reported interval.
                if (diskDepth(c, p) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(
                        f.parameter[f.numIntersections - 1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('a ray whose origin is inside the disk starts at parameter 0', () => {
        check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3),
            fc.double({ min: 0.5, max: 3, noNaN: true,
                noDefaultInfinity: true }),
            angle2(),
            fc.double({ min: 0, max: 0.9, noNaN: true,
                noDefaultInfinity: true }),
            angle2()),
        ([cx, cy, rad, ao, u, ad]) => {
            const c = Hypersphere.fromCenterRadius(vec(cx, cy), rad);
            const o = add(c.center,
                vec(rad * u * Math.cos(ao), rad * u * Math.sin(ao)));
            const r = Ray.fromOriginDirection(o,
                vec(Math.cos(ad), Math.sin(ad)));
            const f = fiq.find(r, c);
            expect(f.intersect).toBe(true);
            expect(f.parameter[0]).toBe(0);
            expect(f.numIntersections).toBe(2);
            expect(f.parameter[1]).toBeGreaterThan(0);
            expectClose(length(sub(f.point[1], c.center)), rad, 1e-9, 1e-9);
        });
    });

    it('a ray pointing away from a disk it misses reports nothing', () => {
        check(rayCircle, ({ ray: r, circle: c }) => {
            const back = Ray.fromOriginDirection(r.origin,
                mul(-1, r.direction));
            const f = fiq.find(r, c);
            const b = fiq.find(back, c);
            const line = Line.fromOriginDirection(r.origin, r.direction);
            const lc = lcq.find(line, c);
            // The union of the two opposite rays is the whole line.
            expect(f.intersect || b.intersect).toBe(lc.intersect);
        });
    });

    it('the exported DoQuery reproduces the class result', () => {
        check(rayCircle, ({ ray: r, circle: c }) => {
            const result = defaultIntrRay2Circle2FIResult();
            intrRay2Circle2FIDoQuery(r.origin, r.direction, c, result);
            const f = fiq.find(r, c);
            expect(result.intersect).toBe(f.intersect);
            expect(result.numIntersections).toBe(f.numIntersections);
            expect(result.parameter).toEqual(f.parameter);
            // DoQuery does not fill the points.
            expect(result.point[0].values).toEqual([0, 0]);
            expect(result.point[1].values).toEqual([0, 0]);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(rayCircle, angle2(), wellScaled(-2, 2),
            wellScaled(-2, 2)),
        ([{ ray: r, circle: c }, ang, tx, ty]) => {
            if (Math.abs(rayDiscriminant(r, c)) < 1e-6) {
                return;   // grazing
            }
            const ca = Math.cos(ang), sa = Math.sin(ang);
            const rot = (v: Vector): Vector => vec(
                ca * v.get(0) - sa * v.get(1),
                sa * v.get(0) + ca * v.get(1));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty));
            const r2 = Ray.fromOriginDirection(xf(r.origin), rot(r.direction));
            const c2 = Hypersphere.fromCenterRadius(xf(c.center), c.radius);
            const f0 = fiq.find(r, c);
            const f1 = fiq.find(r2, c2);
            expect(f1.intersect).toBe(f0.intersect);
            expect(f1.numIntersections).toBe(f0.numIntersections);
            for (let i = 0; i < f0.numIntersections; ++i) {
                expectClose(f1.parameter[i], f0.parameter[i], 1e-7, 1e-8);
                expectVectorClose(f1.point[i], xf(f0.point[i]), 1e-7, 1e-8);
            }
        });
    });

    it('a zero-radius circle is met only when the ray passes through it',
        () => {
            const c = circle([1, 0], 0);
            expect(fiq.find(ray([-5, 0], [1, 0]), c).intersect).toBe(true);
            expect(fiq.find(ray([-5, 1], [1, 0]), c).intersect).toBe(false);
            // Pointing away.
            expect(fiq.find(ray([-5, 0], [-1, 0]), c).intersect).toBe(false);
        });
});
