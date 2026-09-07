import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import {
    IntrLine3Ellipsoid3TI,
    IntrLine3Ellipsoid3FI,
    defaultIntrLine3Ellipsoid3FIResult,
    intrLine3Ellipsoid3FIDoQuery
} from '../src/IntrLine3Ellipsoid3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function line(origin: number[], direction: number[]): Line {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Line.fromOriginDirection(Vector.fromArray(origin), d);
}

function ellipsoid(center: Vector, axis: Vector[], extent: Vector):
    Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

// (X-C)^T*M*(X-C) - 1; zero on the surface, negative inside.
function quadratic(e: Hyperellipsoid, X: Vector): number {
    const diff = sub(X, e.center);
    let sum = 0;
    for (let d = 0; d < 3; ++d) {
        const t = dot(diff, e.axis[d]) / e.extent.values[d];
        sum += t * t;
    }
    return sum - 1;
}

const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

describe('IntrLine3Ellipsoid3', () => {
    const ti = new IntrLine3Ellipsoid3TI();
    const fi = new IntrLine3Ellipsoid3FI();

    it('finds the two crossings of the unit sphere', () => {
        const sphere = ellipsoid(vec(0, 0, 0), unitAxes, vec(1, 1, 1));
        const result = fi.find(line([0, 0, 0], [1, 0, 0]), sphere);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('finds the two crossings of a stretched ellipsoid', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(3, 2, 1));
        const alongX = fi.find(line([0, 0, 0], [1, 0, 0]), e);
        expect(alongX.parameter[0]).toBeCloseTo(-3, 12);
        expect(alongX.parameter[1]).toBeCloseTo(3, 12);

        const alongZ = fi.find(line([0, 0, 0], [0, 0, 1]), e);
        expect(alongZ.parameter[0]).toBeCloseTo(-1, 12);
        expect(alongZ.parameter[1]).toBeCloseTo(1, 12);
    });

    it('reports a tangent line as a single intersection', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 1, 1));
        const result = fi.find(line([0, 1, 0], [1, 0, 0]), e);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBe(result.parameter[0]);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
        expect(result.point[1].values[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection for a line that misses the ellipsoid', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 1, 1));
        const l = line([0, 2, 0], [1, 0, 0]);
        const result = fi.find(l, e);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(l, e).intersect).toBe(false);
    });

    it('agrees between the test and find queries at tangency', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 1, 1));
        const l = line([0, 1, 0], [1, 0, 0]);
        // The test query uses discr >= 0, so tangency counts.
        expect(ti.test(l, e).intersect).toBe(true);
        expect(fi.find(l, e).intersect).toBe(true);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(1, 1, 1));
        const result = defaultIntrLine3Ellipsoid3FIResult();
        const d = vec(1, 0, 0);
        intrLine3Ellipsoid3FIDoQuery(vec(0, 0, 0), d, e, result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        // The points are left at their default values.
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('rejects non-3D ellipsoids', () => {
        const e2 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(() => ti.test(line([0, 0, 0], [1, 0, 0]), e2)).toThrow();
    });

    it('puts the reported points on the ellipsoid for random lines', () => {
        let seed = 555777;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const w = vec(0.3, -0.7, 0.5);
        normalize(w);
        const basis = [w.clone(), Vector.zero(3), Vector.zero(3)];
        computeOrthogonalComplement3(1, basis, false);
        const e = ellipsoid(vec(0.5, -1, 2), basis, vec(3, 1.5, 0.75));

        let numHits = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const l = line([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(l, e);
            expect(ti.test(l, e).intersect).toBe(result.intersect);
            if (!result.intersect) {
                continue;
            }
            ++numHits;
            for (let i = 0; i < 2; ++i) {
                expect(Math.abs(quadratic(e, result.point[i])))
                    .toBeLessThan(1e-8);
                const onLine = add(l.origin,
                    mul(result.parameter[i], l.direction));
                expect(Math.sqrt(dot(sub(onLine, result.point[i]),
                    sub(onLine, result.point[i])))).toBeLessThan(1e-9);
            }
            if (result.numIntersections === 2) {
                // The midpoint of the chord is strictly inside.
                const mid = add(l.origin, mul(
                    0.5 * (result.parameter[0] + result.parameter[1]),
                    l.direction));
                expect(quadratic(e, mid)).toBeLessThan(0);
                expect(result.parameter[0]).toBeLessThan(result.parameter[1]);
            }
        }
        expect(numHits).toBeGreaterThan(10);
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrLine3Ellipsoid3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { length } from '../src/Vector.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no frame component is a
// subnormal that would underflow when squared.
const angleE = () => wellScaled(-Math.PI, Math.PI);
const extentE = () => fc.double(
    { min: 0.3, max: 3, noNaN: true, noDefaultInfinity: true });

// R = Rz(a)*Ry(b)*Rx(c); the columns are the ellipsoid axes.
function rotFrameE(a: number, b: number, c: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const cc = Math.cos(c), sc = Math.sin(c);
    return [
        vec(ca * cb, sa * cb, -sb),
        vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc, cb * sc),
        vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc, cb * cc)];
}

function unitDirE(th: number, ph: number): Vector {
    const d = vec(Math.cos(th) * Math.cos(ph), Math.sin(th) * Math.cos(ph),
        Math.sin(ph));
    normalize(d);
    return d;
}

const lineEllipsoid = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    angleE(), angleE(),
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    angleE(), angleE(), angleE(),
    extentE(), extentE(), extentE()
).map(([ox, oy, oz, th, ph, cx, cy, cz, a, b, c, e0, e1, e2]) => ({
    line: Line.fromOriginDirection(vec(ox, oy, oz), unitDirE(th, ph)),
    ellipsoid: Hyperellipsoid.fromCenterAxisExtent(vec(cx, cy, cz),
        rotFrameE(a, b, c), vec(e0, e1, e2))
}));

// The discriminant of the quadratic the query solves. Near zero the line
// grazes the ellipsoid and the discrete answers are ambiguous.
function ellipsoidDiscriminant(l: Line, e: Hyperellipsoid): number {
    const diff = sub(l.origin, e.center);
    let a2 = 0, a1 = 0, a0 = -1;
    for (let d = 0; d < 3; ++d) {
        const inv = 1 / (e.extent.values[d] * e.extent.values[d]);
        const dd = dot(l.direction, e.axis[d]);
        const pd = dot(diff, e.axis[d]);
        a2 += dd * dd * inv;
        a1 += dd * pd * inv;
        a0 += pd * pd * inv;
    }
    return a1 * a1 - a0 * a2;
}

describe('IntrLine3Ellipsoid3 verification', () => {
    const tiq = new IntrLine3Ellipsoid3TI();
    const fiq = new IntrLine3Ellipsoid3FI();

    it('TI and FI agree away from tangency', () => {
        check(lineEllipsoid, ({ line: l, ellipsoid: e }) => {
            // TI accepts discr >= 0 and FI accepts discr > 0 or == 0, so the
            // two agree exactly; the tolerance band only guards against the
            // two expressions rounding differently, which they cannot here
            // because both compute a1*a1 - a0*a2 from the same terms.
            expect(tiq.test(l, e).intersect).toBe(fiq.find(l, e).intersect);
        });
    });

    it('reported points lie on the line and on the ellipsoid', () => {
        check(lineEllipsoid, ({ line: l, ellipsoid: e }) => {
            const f = fiq.find(l, e);
            if (!f.intersect) {
                expect(f.numIntersections).toBe(0);
                expect(f.point[0].values).toEqual([0, 0, 0]);
                expect(f.point[1].values).toEqual([0, 0, 0]);
                return;
            }
            // Upstream fills both entries whenever intersect is true.
            for (let i = 0; i < 2; ++i) {
                expectVectorClose(f.point[i],
                    add(l.origin, mul(f.parameter[i], l.direction)), 0, 0);
                // The quadratic residual is conditioned by the square root
                // taken to form the parameter; the ellipsoid extents scale it.
                expectClose(quadratic(e, f.point[i]), 0, 1e-7, 1e-8);
            }
            if (f.numIntersections === 1) {
                expect(f.parameter[0]).toBe(f.parameter[1]);
            }
            else {
                expect(f.numIntersections).toBe(2);
                expect(f.parameter[0]).toBeLessThan(f.parameter[1]);
            }
        });
    });

    it('the roots satisfy the quadratic identities', () => {
        check(lineEllipsoid, ({ line: l, ellipsoid: e }) => {
            const f = fiq.find(l, e);
            if (f.numIntersections !== 2) { return; }
            const diff = sub(l.origin, e.center);
            let a2 = 0, a1 = 0, a0 = -1;
            for (let d = 0; d < 3; ++d) {
                const inv = 1 / (e.extent.values[d] * e.extent.values[d]);
                const dd = dot(l.direction, e.axis[d]);
                const pd = dot(diff, e.axis[d]);
                a2 += dd * dd * inv;
                a1 += dd * pd * inv;
                a0 += pd * pd * inv;
            }
            // t0 + t1 = -2*a1/a2 and t0*t1 = a0/a2.
            expectClose(f.parameter[0] + f.parameter[1], -2 * a1 / a2,
                1e-7, 1e-8);
            expectClose(f.parameter[0] * f.parameter[1], a0 / a2, 1e-7, 1e-8);
        });
    });

    it('a fine sweep of the line agrees with the reported interval', () => {
        const rnd = seededRandom(0x4c17b93);
        for (let trial = 0; trial < 150; ++trial) {
            const e = Hyperellipsoid.fromCenterAxisExtent(
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rotFrameE(rnd() * 6, rnd() * 6, rnd() * 6),
                vec(0.4 + rnd() * 2, 0.4 + rnd() * 2, 0.4 + rnd() * 2));
            const d = vec(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d, d) < 1e-4) { continue; }
            normalize(d);
            const l = Line.fromOriginDirection(
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5), d);
            const f = fiq.find(l, e);
            for (let k = 0; k <= 800; ++k) {
                const t = -20 + (40 * k) / 800;
                if (quadratic(e, add(l.origin, mul(t, d))) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-6);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-6);
                }
            }
        }
    }, 30000);

    it('a line through a sampled surface point always intersects', () => {
        check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
            angleE(), angleE(), angleE(), extentE(), extentE(), extentE(),
            angleE(), angleE(), angleE(), angleE()),
        ([cx, cy, cz, a, b, c, e0, e1, e2, th, ph, dth, dph]) => {
            const e = Hyperellipsoid.fromCenterAxisExtent(vec(cx, cy, cz),
                rotFrameE(a, b, c), vec(e0, e1, e2));
            // A point on the surface in the ellipsoid frame.
            const u = vec(e0 * Math.cos(th) * Math.cos(ph),
                e1 * Math.sin(th) * Math.cos(ph), e2 * Math.sin(ph));
            const p = add(e.center, add(mul(u.get(0), e.axis[0]),
                add(mul(u.get(1), e.axis[1]), mul(u.get(2), e.axis[2]))));
            const l = Line.fromOriginDirection(p, unitDirE(dth, dph));
            expect(tiq.test(l, e).intersect).toBe(true);
            const f = fiq.find(l, e);
            expect(f.intersect).toBe(true);
            let found = false;
            for (let i = 0; i < 2; ++i) {
                if (length(sub(f.point[i], p)) < 1e-6 * (1 + e0 + e1 + e2)) {
                    found = true;
                }
            }
            expect(found).toBe(true);
        });
    });

    it('the exported DoQuery reproduces the class result', () => {
        check(lineEllipsoid, ({ line: l, ellipsoid: e }) => {
            const result = defaultIntrLine3Ellipsoid3FIResult();
            intrLine3Ellipsoid3FIDoQuery(l.origin, l.direction, e, result);
            const f = fiq.find(l, e);
            expect(result.intersect).toBe(f.intersect);
            expect(result.numIntersections).toBe(f.numIntersections);
            expect(result.parameter).toEqual(f.parameter);
            expect(result.point[0].values).toEqual([0, 0, 0]);
            expect(result.point[1].values).toEqual([0, 0, 0]);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(lineEllipsoid, angleE(), angleE(), angleE(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2)),
        ([{ line: l, ellipsoid: e }, a1, a2, a3, tx, ty, tz]) => {
            if (Math.abs(ellipsoidDiscriminant(l, e)) < 1e-8) {
                return;   // tangency
            }
            const fr = rotFrameE(a1, a2, a3);
            const rot = (v: Vector): Vector => add(mul(v.get(0), fr[0]),
                add(mul(v.get(1), fr[1]), mul(v.get(2), fr[2])));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty, tz));
            const l2 = Line.fromOriginDirection(xf(l.origin), rot(l.direction));
            const e2 = Hyperellipsoid.fromCenterAxisExtent(xf(e.center),
                [rot(e.axis[0]), rot(e.axis[1]), rot(e.axis[2])], e.extent);
            const f0 = fiq.find(l, e);
            const f1 = fiq.find(l2, e2);
            expect(f1.intersect).toBe(f0.intersect);
            if (!f0.intersect) { return; }
            expect(f1.numIntersections).toBe(f0.numIntersections);
            for (let i = 0; i < 2; ++i) {
                expectClose(f1.parameter[i], f0.parameter[i], 1e-6, 1e-7);
                expectVectorClose(f1.point[i], xf(f0.point[i]), 1e-6, 1e-7);
            }
        });
    });

    it('reduces to the sphere case for equal extents', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 2, 2));
        const f = fiq.find(line([-10, 0, 0], [1, 0, 0]), e);
        expect(f.intersect).toBe(true);
        expect(f.numIntersections).toBe(2);
        expectClose(f.parameter[0], 8, 1e-12, 1e-12);
        expectClose(f.parameter[1], 12, 1e-12, 1e-12);
    });

    it('rejects a non-3D ellipsoid', () => {
        const e2 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        const l = line([0, 0, 0], [1, 0, 0]);
        expect(() => tiq.test(l, e2)).toThrow();
        expect(() => fiq.find(l, e2)).toThrow();
    });
});
