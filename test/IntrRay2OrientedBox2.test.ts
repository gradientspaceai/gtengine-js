import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrRay2AlignedBox2TI,
    IntrRay2AlignedBox2FI
} from '../src/IntrRay2AlignedBox2.js';
import {
    IntrRay2OrientedBox2TI,
    IntrRay2OrientedBox2FI,
    defaultIntrRay2OrientedBox2FIResult
} from '../src/IntrRay2OrientedBox2.js';
import { OrientedBox } from '../src/OrientedBox.js';
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

function rotatedBox(center: Vector, angle: number, extent: Vector):
    OrientedBox {
    const c = Math.cos(angle), s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(center, [vec(c, s), vec(-s, c)],
        extent);
}

// The signed "outside" amount of P relative to the solid box: negative
// strictly inside, zero on the boundary, positive outside.
function boxSignedDepth(box: OrientedBox, P: Vector): number {
    const diff = sub(P, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 2; ++i) {
        const value = Math.abs(dot(diff, box.axis[i])) - box.extent.values[i];
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

describe('IntrRay2OrientedBox2', () => {
    const ti = new IntrRay2OrientedBox2TI();
    const fi = new IntrRay2OrientedBox2FI();
    const unitAxes = [vec(1, 0), vec(0, 1)];

    it('matches the aligned-box query when the axes are standard', () => {
        const center = vec(1, -2);
        const extent = vec(2, 0.5);
        const obox = OrientedBox.fromCenterAxisExtent(center, unitAxes,
            extent);
        const abox = AlignedBox.fromMinMax(sub(center, extent),
            add(center, extent));
        const aTI = new IntrRay2AlignedBox2TI();
        const aFI = new IntrRay2AlignedBox2FI();

        const rays = [
            ray([-10, -2], [1, 0]),
            ray([1, -2], [0, 1]),
            ray([10, -2], [1, 0]),
            ray([-10, 10], [1, 0]),
            ray([-10, -1.5], [1, 0])
        ];
        for (const r of rays) {
            expect(ti.test(r, obox).intersect)
                .toBe(aTI.test(r, abox).intersect);
            const o = fi.find(r, obox);
            const a = aFI.find(r, abox);
            expect(o.intersect).toBe(a.intersect);
            expect(o.numIntersections).toBe(a.numIntersections);
            for (let i = 0; i < o.numIntersections; ++i) {
                expect(o.parameter[i]).toBeCloseTo(a.parameter[i], 12);
            }
        }
    });

    it('finds the entry and exit points of a rotated box', () => {
        // A box rotated 45 degrees about the origin with extents (1,1) is a
        // diamond with vertices at distance sqrt(2) along the axes.
        const box = rotatedBox(vec(0, 0), Math.PI / 4, vec(1, 1));
        const result = fi.find(ray([-5, 0], [1, 0]), box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(5 - Math.SQRT2, 9);
        expect(result.parameter[1]).toBeCloseTo(5 + Math.SQRT2, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-Math.SQRT2, 9);
        expect(result.point[0].values[1]).toBeCloseTo(0, 9);
        expect(result.point[1].values[0]).toBeCloseTo(Math.SQRT2, 9);
    });

    it('clips at the ray origin when the origin is inside', () => {
        const box = rotatedBox(vec(0, 0), Math.PI / 6, vec(2, 1));
        const result = fi.find(ray([0, 0], [1, 0]), box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeGreaterThan(0);
    });

    it('rejects a box behind the ray', () => {
        const box = rotatedBox(vec(0, 0), 0.3, vec(1, 1));
        const r = ray([5, 0], [1, 0]);
        expect(fi.find(r, box).intersect).toBe(false);
        expect(ti.test(r, box).intersect).toBe(false);
    });

    it('has the documented default result', () => {
        const result = defaultIntrRay2OrientedBox2FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('rejects non-2D boxes', () => {
        const box3 = new OrientedBox(3);
        const r = ray([0, 0], [1, 0]);
        expect(() => ti.test(r, box3)).toThrow();
        expect(() => fi.find(r, box3)).toThrow();
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 606060;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const box = rotatedBox(vec(0.5, -0.25), 0.7, vec(1.5, 0.75));
        for (let trial = 0; trial < 400; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, box);
            expect(ti.test(r, box).intersect).toBe(result.intersect);

            let sampledHit = false;
            for (let k = 0; k <= 2000; ++k) {
                const t = k * 0.01;
                if (boxSignedDepth(box,
                    add(r.origin, mul(t, r.direction))) < -1e-6) {
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
                expect(boxSignedDepth(box, P)).toBeLessThan(1e-8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrRay2OrientedBox2.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { IntrLine2OrientedBox2FI } from '../src/IntrLine2OrientedBox2.js';
import { Line } from '../src/Line.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so a direction built from
// (cos a, sin a) never has a subnormal component; a subnormal component makes
// the exact DotPerp tests inside the line-line query underflow.
const angle2 = () => wellScaled(-Math.PI, Math.PI);

const extent2 = () => fc.double(
    { min: 0.25, max: 3, noNaN: true, noDefaultInfinity: true });

const rayBox = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), angle2(),
    wellScaled(-3, 3), wellScaled(-3, 3), angle2(),
    extent2(), extent2()
).map(([ox, oy, ad, cx, cy, ab, e0, e1]) => ({
    ray: Ray.fromOriginDirection(vec(ox, oy), vec(Math.cos(ad), Math.sin(ad))),
    box: rotatedBox(vec(cx, cy), ab, vec(e0, e1))
}));

// The largest signed slab depth: negative strictly inside the box.
function obDepth(box: OrientedBox, P: Vector): number {
    const diff = sub(P, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 2; ++i) {
        const v = Math.abs(dot(diff, box.axis[i])) - box.extent.values[i];
        if (v > worst) { worst = v; }
    }
    return worst;
}

// obDepth is a maximum of convex functions, so it is convex along a ray; a
// ternary search finds its minimum over t in [0, tMax] reliably. A value near
// zero means the ray grazes the box, where the TI slab test and the FI
// interval clip round differently and the discrete answers are ambiguous.
function minRayDepth(box: OrientedBox, r: Ray, tMax = 40): number {
    let lo = 0, hi = tMax;
    for (let k = 0; k < 200; ++k) {
        const a = lo + (hi - lo) / 3, b = hi - (hi - lo) / 3;
        const fa = obDepth(box, add(r.origin, mul(a, r.direction)));
        const fb = obDepth(box, add(r.origin, mul(b, r.direction)));
        if (fa < fb) { hi = b; } else { lo = a; }
    }
    return obDepth(box, add(r.origin, mul(0.5 * (lo + hi), r.direction)));
}

describe('IntrRay2OrientedBox2 verification', () => {
    const tiq = new IntrRay2OrientedBox2TI();
    const fiq = new IntrRay2OrientedBox2FI();
    const lbq = new IntrLine2OrientedBox2FI();

    it('TI and FI always agree on intersect', () => {
        check(rayBox, ({ ray: r, box: b }) => {
            if (Math.abs(minRayDepth(b, r)) < 1e-6) {
                return;   // grazing; the two formulations round differently
            }
            expect(tiq.test(r, b).intersect).toBe(fiq.find(r, b).intersect);
        });
    });

    it('is bit identical to the aligned-box query on the box-frame ray',
        () => {
            check(rayBox, ({ ray: r, box: b }) => {
                const diff = sub(r.origin, b.center);
                const o = vec(dot(diff, b.axis[0]), dot(diff, b.axis[1]));
                const d = vec(dot(r.direction, b.axis[0]),
                    dot(r.direction, b.axis[1]));
                const alignedRay = Ray.fromOriginDirection(o, d);
                const alignedBox = AlignedBox.fromMinMax(
                    mul(-1, b.extent), b.extent);
                const a = new IntrRay2AlignedBox2FI().find(alignedRay,
                    alignedBox);
                const f = fiq.find(r, b);
                expect(f.intersect).toBe(a.intersect);
                expect(f.numIntersections).toBe(a.numIntersections);
                expect(f.parameter).toEqual(a.parameter);
                expect(new IntrRay2OrientedBox2TI().test(r, b).intersect)
                    .toBe(new IntrRay2AlignedBox2TI()
                        .test(alignedRay, alignedBox).intersect);
            });
        });

    it('the ray hit is the line hit clipped to t >= 0', () => {
        check(rayBox, ({ ray: r, box: b }) => {
            const line = Line.fromOriginDirection(r.origin, r.direction);
            const lb = lbq.find(line, b);
            const f = fiq.find(r, b);
            if (!lb.intersect) {
                expect(f.intersect).toBe(false);
                return;
            }
            const t0 = Math.max(lb.parameter[0], 0);
            const t1 = lb.parameter[1];
            if (t1 < 0) {
                expect(f.intersect).toBe(false);
                return;
            }
            expect(f.intersect).toBe(true);
            expect(f.numIntersections).toBe(t0 < t1 ? 2 : 1);
            // Normalize the -0/+0 tie: Math.max(-0, 0) is +0 while the
            // query's clip can produce -0, and toBe uses Object.is.
            expect(f.parameter[0] + 0).toBe(t0 + 0);
        });
    });

    it('reported points lie on the ray and in the closed box', () => {
        check(rayBox, ({ ray: r, box: b }) => {
            const f = fiq.find(r, b);
            for (let i = 0; i < f.numIntersections; ++i) {
                expect(f.parameter[i]).toBeGreaterThanOrEqual(0);
                expectVectorClose(f.point[i],
                    add(r.origin, mul(f.parameter[i], r.direction)), 0, 0);
                expect(obDepth(b, f.point[i])).toBeLessThanOrEqual(1e-9);
            }
        });
    });

    it('a fine sweep of the ray agrees with the reported interval', () => {
        const rnd = seededRandom(0x77c31a5);
        for (let trial = 0; trial < 200; ++trial) {
            const b = rotatedBox(vec(rnd() * 4 - 2, rnd() * 4 - 2),
                rnd() * 2 * Math.PI,
                vec(0.3 + rnd() * 2, 0.3 + rnd() * 2));
            const a = rnd() * 2 * Math.PI;
            const d = vec(Math.cos(a), Math.sin(a));
            const r = Ray.fromOriginDirection(
                vec(rnd() * 8 - 4, rnd() * 8 - 4), d);
            const f = fiq.find(r, b);
            for (let k = 0; k <= 1000; ++k) {
                const t = (14 * k) / 1000;
                const p = add(r.origin, mul(t, d));
                if (obDepth(b, p) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(
                        f.parameter[f.numIntersections - 1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('a ray whose origin is inside the box starts at parameter 0', () => {
        check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3), angle2(),
            extent2(), extent2(),
            fc.double({ min: -0.85, max: 0.85, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -0.85, max: 0.85, noNaN: true,
                noDefaultInfinity: true }),
            angle2()),
        ([cx, cy, ab, e0, e1, u0, u1, ad]) => {
            const b = rotatedBox(vec(cx, cy), ab, vec(e0, e1));
            const o = add(b.center, add(mul(u0 * e0, b.axis[0]),
                mul(u1 * e1, b.axis[1])));
            const r = Ray.fromOriginDirection(o,
                vec(Math.cos(ad), Math.sin(ad)));
            const f = fiq.find(r, b);
            expect(f.intersect).toBe(true);
            expect(f.parameter[0]).toBe(0);
        });
    });

    it('a ray fired at an interior box point hits and its reverse misses',
        () => {
            check(fc.tuple(wellScaled(-2, 2), wellScaled(-2, 2), angle2(),
                extent2(), extent2(),
                fc.double({ min: -0.7, max: 0.7, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: -0.7, max: 0.7, noNaN: true,
                    noDefaultInfinity: true }),
                angle2(),
                fc.double({ min: 6, max: 12, noNaN: true,
                    noDefaultInfinity: true })),
            ([cx, cy, ab, e0, e1, u0, u1, ad, dist]) => {
                const b = rotatedBox(vec(cx, cy), ab, vec(e0, e1));
                const target = add(b.center, add(mul(u0 * e0, b.axis[0]),
                    mul(u1 * e1, b.axis[1])));
                const d = vec(Math.cos(ad), Math.sin(ad));
                const origin = sub(target, mul(dist, d));
                expect(fiq.find(Ray.fromOriginDirection(origin, d), b)
                    .intersect).toBe(true);
                // The origin is at distance >= 6 from a box whose extents are
                // at most 3, so the reversed ray cannot reach the box.
                if (obDepth(b, origin) > 0) {
                    expect(fiq.find(Ray.fromOriginDirection(origin,
                        mul(-1, d)), b).intersect).toBe(false);
                }
            });
        });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(rayBox, angle2(), wellScaled(-2, 2), wellScaled(-2, 2)),
            ([{ ray: r, box: b }, ang, tx, ty]) => {
                const ca = Math.cos(ang), sa = Math.sin(ang);
                const rot = (v: Vector): Vector => vec(
                    ca * v.get(0) - sa * v.get(1),
                    sa * v.get(0) + ca * v.get(1));
                const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty));
                const r2 = Ray.fromOriginDirection(xf(r.origin),
                    rot(r.direction));
                const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                    [rot(b.axis[0]), rot(b.axis[1])], b.extent);
                if (Math.abs(minRayDepth(b, r)) < 1e-6) {
                    return;   // grazing corner or edge
                }
                const f0 = fiq.find(r, b);
                const f1 = fiq.find(r2, b2);
                if (f0.numIntersections !== f1.numIntersections) {
                    return;
                }
                for (let i = 0; i < f0.numIntersections; ++i) {
                    expectClose(f1.parameter[i], f0.parameter[i], 1e-7, 1e-8);
                    expectVectorClose(f1.point[i], xf(f0.point[i]),
                        1e-7, 1e-8);
                }
            });
    });

    it('a zero-extent box is a point the ray must pass through', () => {
        const b = rotatedBox(vec(1, 1), 0.3, vec(0, 0));
        expect(fiq.find(ray([-4, 1], [1, 0]), b).intersect).toBe(true);
        expect(fiq.find(ray([-4, 1.5], [1, 0]), b).intersect).toBe(false);
        const f = fiq.find(ray([-4, 1], [1, 0]), b);
        expect(f.numIntersections).toBe(1);
        expectClose(f.parameter[0], 5, 1e-12, 1e-12);
    });

    it('a miss leaves every field of the result at its default', () => {
        const b = rotatedBox(vec(0, 0), 0.4, vec(1, 1));
        const f = fiq.find(ray([-5, 5], [1, 0]), b);
        expect(f.intersect).toBe(false);
        expect(f.numIntersections).toBe(0);
        expect(f.parameter).toEqual([0, 0]);
        expect(f.point[0].values).toEqual([0, 0]);
        expect(f.point[1].values).toEqual([0, 0]);
        expect(defaultIntrRay2OrientedBox2FIResult().numIntersections).toBe(0);
    });
});
