import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrRay3AlignedBox3TI,
    IntrRay3AlignedBox3FI
} from '../src/IntrRay3AlignedBox3.js';
import {
    IntrRay3OrientedBox3TI,
    IntrRay3OrientedBox3FI,
    defaultIntrRay3OrientedBox3FIResult
} from '../src/IntrRay3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function orthonormalFrame(w: Vector): Vector[] {
    const v = [w.clone(), Vector.zero(3), Vector.zero(3)];
    normalize(v[0]);
    computeOrthogonalComplement3(1, v, false);
    return [v[0], v[1], v[2]];
}

function boxSignedDepth(box: OrientedBox, P: Vector): number {
    const diff = sub(P, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 3; ++i) {
        const value = Math.abs(dot(diff, box.axis[i])) - box.extent.values[i];
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

describe('IntrRay3OrientedBox3', () => {
    const ti = new IntrRay3OrientedBox3TI();
    const fi = new IntrRay3OrientedBox3FI();
    const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

    it('matches the aligned-box query when the axes are standard', () => {
        const center = vec(1, -2, 3);
        const extent = vec(2, 0.5, 1);
        const obox = OrientedBox.fromCenterAxisExtent(center, unitAxes,
            extent);
        const abox = AlignedBox.fromMinMax(sub(center, extent),
            add(center, extent));
        const aTI = new IntrRay3AlignedBox3TI();
        const aFI = new IntrRay3AlignedBox3FI();

        const rays = [
            ray([-10, -2, 3], [1, 0, 0]),
            ray([1, -2, 3], [0, 1, 0]),
            ray([10, -2, 3], [1, 0, 0]),
            ray([-10, 10, 3], [1, 0, 0]),
            ray([0, 0, 0], [1, -2, 3])
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

    it('finds entry and exit points of a rotated box', () => {
        // Rotate 45 degrees about z; the box with extents (1,1,1) becomes a
        // diamond cross-section reaching sqrt(2) along x.
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0),
            [vec(c, c, 0), vec(-c, c, 0), vec(0, 0, 1)], vec(1, 1, 1));
        const result = fi.find(ray([-5, 0, 0], [1, 0, 0]), box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(5 - Math.SQRT2, 9);
        expect(result.parameter[1]).toBeCloseTo(5 + Math.SQRT2, 9);
    });

    it('clips at the ray origin when the origin is inside', () => {
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0), unitAxes,
            vec(2, 1, 1));
        const result = fi.find(ray([0, 0, 0], [1, 0, 0]), box);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
    });

    it('rejects a box behind the ray', () => {
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0), unitAxes,
            vec(1, 1, 1));
        const r = ray([5, 0, 0], [1, 0, 0]);
        expect(fi.find(r, box).intersect).toBe(false);
        expect(ti.test(r, box).intersect).toBe(false);
    });

    it('has the documented default result', () => {
        const result = defaultIntrRay3OrientedBox3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('rejects non-3D boxes', () => {
        const box2 = new OrientedBox(2);
        const r = ray([0, 0, 0], [1, 0, 0]);
        expect(() => ti.test(r, box2)).toThrow();
        expect(() => fi.find(r, box2)).toThrow();
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 8080808;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const axes = orthonormalFrame(vec(0.3, 0.8, -0.5));
        const box = OrientedBox.fromCenterAxisExtent(vec(0.5, -0.25, 1), axes,
            vec(1.5, 0.75, 1));

        for (let trial = 0; trial < 300; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, box);
            expect(ti.test(r, box).intersect).toBe(result.intersect);

            let sampledHit = false;
            for (let k = 0; k <= 1500; ++k) {
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

            if (result.intersect) {
                for (let i = 0; i < 2; ++i) {
                    expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-12);
                    const P = add(r.origin,
                        mul(result.parameter[i], r.direction));
                    expect(sub(P, result.point[i]).values[0])
                        .toBeCloseTo(0, 9);
                    expect(boxSignedDepth(box, P)).toBeLessThan(1e-8);
                }
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrRay3OrientedBox3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { IntrLine3OrientedBox3FI } from '../src/IntrLine3OrientedBox3.js';
import { Line } from '../src/Line.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no frame component is a
// subnormal that would underflow when squared.
const angle3 = () => wellScaled(-Math.PI, Math.PI);
const extent3 = () => fc.double(
    { min: 0.25, max: 3, noNaN: true, noDefaultInfinity: true });

// R = Rz(a)*Ry(b)*Rx(c); the columns are the box axes.
function rotFrame(a: number, b: number, c: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const cc = Math.cos(c), sc = Math.sin(c);
    return [
        vec(ca * cb, sa * cb, -sb),
        vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc, cb * sc),
        vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc, cb * cc)];
}

const rayBox3 = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    angle3(), angle3(),
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    angle3(), angle3(), angle3(),
    extent3(), extent3(), extent3()
).map(([ox, oy, oz, th, ph, cx, cy, cz, a, b, c, e0, e1, e2]) => {
    const d = vec(Math.cos(th) * Math.cos(ph), Math.sin(th) * Math.cos(ph),
        Math.sin(ph));
    normalize(d);
    return {
        ray: Ray.fromOriginDirection(vec(ox, oy, oz), d),
        box: OrientedBox.fromCenterAxisExtent(vec(cx, cy, cz),
            rotFrame(a, b, c), vec(e0, e1, e2))
    };
});

// boxSignedDepth is a maximum of convex functions, hence convex along the
// ray; a ternary search finds its minimum reliably. A value near zero means
// the ray grazes the box, where the TI separating-axis test and the FI
// interval clip round differently and the discrete answers are ambiguous.
function minRayBoxDepth(box: OrientedBox, r: Ray, tMax = 40): number {
    let lo = 0, hi = tMax;
    for (let k = 0; k < 200; ++k) {
        const p = lo + (hi - lo) / 3, q = hi - (hi - lo) / 3;
        const fp = boxSignedDepth(box, add(r.origin, mul(p, r.direction)));
        const fq = boxSignedDepth(box, add(r.origin, mul(q, r.direction)));
        if (fp < fq) { hi = q; } else { lo = p; }
    }
    return boxSignedDepth(box, add(r.origin,
        mul(0.5 * (lo + hi), r.direction)));
}

describe('IntrRay3OrientedBox3 verification', () => {
    const tiq = new IntrRay3OrientedBox3TI();
    const fiq = new IntrRay3OrientedBox3FI();
    const lbq = new IntrLine3OrientedBox3FI();

    it('TI and FI agree away from grazing configurations', () => {
        check(rayBox3, ({ ray: r, box: b }) => {
            if (Math.abs(minRayBoxDepth(b, r)) < 1e-6) {
                return;
            }
            expect(tiq.test(r, b).intersect).toBe(fiq.find(r, b).intersect);
        });
    });

    it('is bit identical to the aligned-box query on the box-frame ray',
        () => {
            check(rayBox3, ({ ray: r, box: b }) => {
                const diff = sub(r.origin, b.center);
                const o = vec(dot(diff, b.axis[0]), dot(diff, b.axis[1]),
                    dot(diff, b.axis[2]));
                const d = vec(dot(r.direction, b.axis[0]),
                    dot(r.direction, b.axis[1]), dot(r.direction, b.axis[2]));
                const alignedRay = Ray.fromOriginDirection(o, d);
                const alignedBox = AlignedBox.fromMinMax(mul(-1, b.extent),
                    b.extent);
                const a = new IntrRay3AlignedBox3FI().find(alignedRay,
                    alignedBox);
                const f = fiq.find(r, b);
                expect(f.intersect).toBe(a.intersect);
                expect(f.numIntersections).toBe(a.numIntersections);
                expect(f.parameter).toEqual(a.parameter);
                expect(tiq.test(r, b).intersect).toBe(
                    new IntrRay3AlignedBox3TI().test(alignedRay, alignedBox)
                        .intersect);
            });
        });

    it('the ray hit is the line hit clipped to t >= 0', () => {
        check(rayBox3, ({ ray: r, box: b }) => {
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
            // Normalize the -0/+0 tie (toBe uses Object.is).
            expect(f.parameter[0] + 0).toBe(t0 + 0);
        });
    });

    it('fills both point entries whenever it intersects', () => {
        check(rayBox3, ({ ray: r, box: b }) => {
            const f = fiq.find(r, b);
            if (!f.intersect) {
                expect(f.point[0].values).toEqual([0, 0, 0]);
                expect(f.point[1].values).toEqual([0, 0, 0]);
                return;
            }
            // Upstream fills point[0] and point[1] whenever intersect is
            // true, duplicating the point when numIntersections is 1.
            for (let i = 0; i < 2; ++i) {
                expect(f.parameter[i]).toBeGreaterThanOrEqual(0);
                expectVectorClose(f.point[i],
                    add(r.origin, mul(f.parameter[i], r.direction)), 0, 0);
                expect(boxSignedDepth(b, f.point[i]))
                    .toBeLessThanOrEqual(1e-9);
            }
            if (f.numIntersections === 1) {
                expect(f.parameter[0]).toBe(f.parameter[1]);
            }
        });
    });

    it('a fine sweep of the ray agrees with the reported interval', () => {
        const rnd = seededRandom(0x31c9f2);
        for (let trial = 0; trial < 150; ++trial) {
            const b = OrientedBox.fromCenterAxisExtent(
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rotFrame(rnd() * 6, rnd() * 6, rnd() * 6),
                vec(0.3 + rnd() * 2, 0.3 + rnd() * 2, 0.3 + rnd() * 2));
            const d = vec(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d, d) < 1e-4) { continue; }
            normalize(d);
            const r = Ray.fromOriginDirection(
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5), d);
            const f = fiq.find(r, b);
            for (let k = 0; k <= 500; ++k) {
                const t = (18 * k) / 500;
                const p = add(r.origin, mul(t, d));
                if (boxSignedDepth(b, p) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('a ray whose origin is inside the box starts at parameter 0', () => {
        check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
            angle3(), angle3(), angle3(),
            extent3(), extent3(), extent3(),
            fc.double({ min: -0.8, max: 0.8, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -0.8, max: 0.8, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -0.8, max: 0.8, noNaN: true,
                noDefaultInfinity: true }),
            angle3(), angle3()),
        ([cx, cy, cz, a, b, c, e0, e1, e2, u0, u1, u2, th, ph]) => {
            const box = OrientedBox.fromCenterAxisExtent(vec(cx, cy, cz),
                rotFrame(a, b, c), vec(e0, e1, e2));
            const o = add(box.center,
                add(mul(u0 * e0, box.axis[0]),
                    add(mul(u1 * e1, box.axis[1]),
                        mul(u2 * e2, box.axis[2]))));
            const d = vec(Math.cos(th) * Math.cos(ph),
                Math.sin(th) * Math.cos(ph), Math.sin(ph));
            normalize(d);
            const f = fiq.find(Ray.fromOriginDirection(o, d), box);
            expect(f.intersect).toBe(true);
            expect(f.parameter[0]).toBe(0);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(rayBox3, angle3(), angle3(), angle3(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2)),
        ([{ ray: r, box: b }, a1, a2, a3, tx, ty, tz]) => {
            if (Math.abs(minRayBoxDepth(b, r)) < 1e-6) {
                return;   // grazing face, edge or corner
            }
            const f = rotFrame(a1, a2, a3);
            const rot = (v: Vector): Vector => add(mul(v.get(0), f[0]),
                add(mul(v.get(1), f[1]), mul(v.get(2), f[2])));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty, tz));
            const r2 = Ray.fromOriginDirection(xf(r.origin), rot(r.direction));
            const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                [rot(b.axis[0]), rot(b.axis[1]), rot(b.axis[2])], b.extent);
            const f0 = fiq.find(r, b);
            const f1 = fiq.find(r2, b2);
            expect(f1.intersect).toBe(f0.intersect);
            if (!f0.intersect) { return; }
            expect(f1.numIntersections).toBe(f0.numIntersections);
            for (let i = 0; i < 2; ++i) {
                expectClose(f1.parameter[i], f0.parameter[i], 1e-7, 1e-8);
                expectVectorClose(f1.point[i], xf(f0.point[i]), 1e-7, 1e-8);
            }
        });
    });

    it('rejects a non-3D box', () => {
        const box2 = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        const r = ray([0, 0, 0], [1, 0, 0]);
        expect(() => fiq.find(r, box2)).toThrow();
        expect(() => tiq.test(r, box2)).toThrow();
    });
});
