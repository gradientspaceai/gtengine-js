import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import { IntrHalfspace3Ellipsoid3TI } from '../src/IntrHalfspace3Ellipsoid3.js';
import { Vector, add, dot, mul, normalize } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function halfspace(normal: number[], constant: number): Halfspace {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Halfspace.fromNormalConstant(n, constant);
}

function ellipsoid(center: Vector, axis: Vector[], extent: Vector):
    Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

function axesFromDirection(w: Vector): Vector[] {
    const v = [w.clone(), Vector.zero(3), Vector.zero(3)];
    computeOrthogonalComplement3(1, v, false);
    return [v[0], v[1], v[2]];
}

// The maximum of Dot(N,X) - c over points X on the ellipsoid surface,
// computed by brute-force sampling. The query reports an intersection when
// this maximum is nonnegative.
function sampledMaxSignedDistance(h: Halfspace, e: Hyperellipsoid): number {
    let maxValue = -Number.MAX_VALUE;
    const numTheta = 120, numPhi = 120;
    for (let i = 0; i <= numTheta; ++i) {
        const theta = (2 * Math.PI * i) / numTheta;
        for (let j = 0; j <= numPhi; ++j) {
            const phi = (Math.PI * j) / numPhi;
            const u = [
                Math.sin(phi) * Math.cos(theta),
                Math.sin(phi) * Math.sin(theta),
                Math.cos(phi)
            ];
            let X = e.center.clone();
            for (let d = 0; d < 3; ++d) {
                X = add(X, mul(e.extent.values[d] * u[d], e.axis[d]));
            }
            const value = dot(h.normal, X) - h.constant;
            if (value > maxValue) {
                maxValue = value;
            }
        }
    }
    return maxValue;
}

describe('IntrHalfspace3Ellipsoid3TI', () => {
    const query = new IntrHalfspace3Ellipsoid3TI();
    const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

    it('handles a unit sphere against a moving plane', () => {
        const sphere = ellipsoid(vec(0, 0, 0), unitAxes, vec(1, 1, 1));

        // The projection interval maximum is 1 - c, so the halfspace
        // x >= c intersects the sphere exactly when c <= 1.
        expect(query.test(halfspace([1, 0, 0], 0.5), sphere).intersect)
            .toBe(true);
        expect(query.test(halfspace([1, 0, 0], 1), sphere).intersect)
            .toBe(true);
        expect(query.test(halfspace([1, 0, 0], 1 + 1e-9), sphere).intersect)
            .toBe(false);
        expect(query.test(halfspace([1, 0, 0], -5), sphere).intersect)
            .toBe(true);
    });

    it('uses the ellipsoid extent along the halfspace normal', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 1, 3));

        // Along x the support is 2, along y it is 1, along z it is 3.
        expect(query.test(halfspace([1, 0, 0], 1.999), e).intersect).toBe(true);
        expect(query.test(halfspace([1, 0, 0], 2.001), e).intersect).toBe(false);
        expect(query.test(halfspace([0, 1, 0], 0.999), e).intersect).toBe(true);
        expect(query.test(halfspace([0, 1, 0], 1.001), e).intersect).toBe(false);
        expect(query.test(halfspace([0, 0, 1], 2.999), e).intersect).toBe(true);
        expect(query.test(halfspace([0, 0, 1], 3.001), e).intersect).toBe(false);
    });

    it('accounts for the ellipsoid center', () => {
        const e = ellipsoid(vec(10, 0, 0), unitAxes, vec(1, 1, 1));
        expect(query.test(halfspace([1, 0, 0], 11 - 1e-6), e).intersect)
            .toBe(true);
        expect(query.test(halfspace([1, 0, 0], 11 + 1e-6), e).intersect)
            .toBe(false);
    });

    it('is invariant to the sign of the ellipsoid axes', () => {
        const e0 = ellipsoid(vec(1, 2, 3), unitAxes, vec(2, 1, 3));
        const e1 = ellipsoid(vec(1, 2, 3),
            [vec(-1, 0, 0), vec(0, -1, 0), vec(0, 0, 1)], vec(2, 1, 3));
        const h = halfspace([1, 2, -1], 4);
        expect(query.test(h, e0).intersect).toBe(query.test(h, e1).intersect);
    });

    it('rejects non-3D inputs', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(1, 1, 1));
        const h2 = Halfspace.fromNormalConstant(
            Vector.fromArray([1, 0]), 0);
        expect(() => query.test(h2, e)).toThrow();
    });

    it('agrees with brute-force sampling on random configurations', () => {
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 40; ++trial) {
            const w = vec(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1);
            if (Math.abs(dot(w, w)) < 1e-6) {
                continue;
            }
            normalize(w);
            const axes = axesFromDirection(w);
            const e = ellipsoid(
                vec(rand() * 4 - 2, rand() * 4 - 2, rand() * 4 - 2), axes,
                vec(0.5 + rand() * 2, 0.5 + rand() * 2, 0.5 + rand() * 2));
            const h = halfspace(
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                rand() * 6 - 3);

            const expectedMax = sampledMaxSignedDistance(h, e);
            // Skip near-tangential configurations where the sampled maximum
            // is not accurate enough to decide the sign.
            if (Math.abs(expectedMax) < 1e-2) {
                continue;
            }
            expect(query.test(h, e).intersect).toBe(expectedMax >= 0);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrHalfspace3Ellipsoid3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { sub } from '../src/Vector.js';

function unitAxes3(): Vector[] {
    return [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];
}

// wellScaled snaps |a| < 1e-3 to exactly zero, so no frame component is a
// subnormal that would underflow when squared.
const angleH = () => wellScaled(-Math.PI, Math.PI);
const extentH = () => fc.double(
    { min: 0.3, max: 3, noNaN: true, noDefaultInfinity: true });

// R = Rz(a)*Ry(b)*Rx(c); the columns are the ellipsoid axes.
function rotFrameH(a: number, b: number, c: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const cc = Math.cos(c), sc = Math.sin(c);
    return [
        vec(ca * cb, sa * cb, -sb),
        vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc, cb * sc),
        vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc, cb * cc)];
}

function unitDirH(th: number, ph: number): Vector {
    const d = vec(Math.cos(th) * Math.cos(ph), Math.sin(th) * Math.cos(ph),
        Math.sin(ph));
    normalize(d);
    return d;
}

const halfspaceEllipsoid = fc.tuple(
    angleH(), angleH(), wellScaled(-4, 4),
    wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
    angleH(), angleH(), angleH(),
    extentH(), extentH(), extentH()
).map(([th, ph, k, cx, cy, cz, a, b, c, e0, e1, e2]) => ({
    halfspace: Halfspace.fromNormalConstant(unitDirH(th, ph), k),
    ellipsoid: Hyperellipsoid.fromCenterAxisExtent(vec(cx, cy, cz),
        rotFrameH(a, b, c), vec(e0, e1, e2))
}));

// The exact support value: max over the ellipsoid surface of Dot(N,X) - c is
// Dot(N,C) - c + sqrt(sum_i (e_i * Dot(N,U_i))^2). The query computes the
// same quantity through M^{-1}; this is an independent derivation.
function supportValue(h: Halfspace, e: Hyperellipsoid): number {
    let sum = 0;
    for (let d = 0; d < 3; ++d) {
        const t = e.extent.values[d] * dot(h.normal, e.axis[d]);
        sum += t * t;
    }
    return dot(h.normal, e.center) - h.constant + Math.sqrt(sum);
}

describe('IntrHalfspace3Ellipsoid3 verification', () => {
    const tiq = new IntrHalfspace3Ellipsoid3TI();

    it('agrees with the closed-form support value', () => {
        check(halfspaceEllipsoid, ({ halfspace: h, ellipsoid: e }) => {
            const s = supportValue(h, e);
            if (Math.abs(s) < 1e-9) {
                return;   // exactly on the plane; the two forms may round
            }
            expect(tiq.test(h, e).intersect).toBe(s >= 0);
        });
    });

    it('agrees with a dense sampling of the ellipsoid surface', () => {
        const rnd = seededRandom(0x1a97c3);
        for (let trial = 0; trial < 60; ++trial) {
            const e = Hyperellipsoid.fromCenterAxisExtent(
                vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rotFrameH(rnd() * 6, rnd() * 6, rnd() * 6),
                vec(0.4 + rnd() * 2, 0.4 + rnd() * 2, 0.4 + rnd() * 2));
            const h = Halfspace.fromNormalConstant(
                unitDirH(rnd() * 6, rnd() * 6 - 3), rnd() * 8 - 4);
            const support = supportValue(h, e);
            if (Math.abs(support) < 1e-3) { continue; }

            // The sampled maximum is a lower bound on the true support; when
            // the true support is comfortably positive the sampling finds a
            // point in the halfspace, and when it is comfortably negative no
            // sample can be in it.
            let best = -Number.MAX_VALUE;
            const n = 120;
            for (let i = 0; i <= n; ++i) {
                const th = (2 * Math.PI * i) / n;
                for (let j = 0; j <= n; ++j) {
                    const ph = -Math.PI / 2 + (Math.PI * j) / n;
                    const u = vec(
                        e.extent.values[0] * Math.cos(th) * Math.cos(ph),
                        e.extent.values[1] * Math.sin(th) * Math.cos(ph),
                        e.extent.values[2] * Math.sin(ph));
                    const p = add(e.center, add(mul(u.get(0), e.axis[0]),
                        add(mul(u.get(1), e.axis[1]),
                            mul(u.get(2), e.axis[2]))));
                    const value = dot(h.normal, p) - h.constant;
                    if (value > best) { best = value; }
                }
            }
            const intersect = tiq.test(h, e).intersect;
            if (best >= 0) {
                expect(intersect).toBe(true);
            }
            if (support < 0) {
                expect(intersect).toBe(false);
                expect(best).toBeLessThan(0);
            }
            expectClose(best, support, 2e-3, 2e-3);
        }
    }, 30000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(halfspaceEllipsoid, angleH(), angleH(), angleH(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2)),
        ([{ halfspace: h, ellipsoid: e }, a1, a2, a3, tx, ty, tz]) => {
            if (Math.abs(supportValue(h, e)) < 1e-6) { return; }
            const fr = rotFrameH(a1, a2, a3);
            const rot = (v: Vector): Vector => add(mul(v.get(0), fr[0]),
                add(mul(v.get(1), fr[1]), mul(v.get(2), fr[2])));
            const t = vec(tx, ty, tz);
            const xf = (v: Vector): Vector => add(rot(v), t);
            const n2 = rot(h.normal);
            // Dot(N', X') = Dot(N, X) + Dot(N', T), so the constant shifts.
            const h2 = Halfspace.fromNormalConstant(n2,
                h.constant + dot(n2, t));
            const e2 = Hyperellipsoid.fromCenterAxisExtent(xf(e.center),
                [rot(e.axis[0]), rot(e.axis[1]), rot(e.axis[2])], e.extent);
            expect(tiq.test(h2, e2).intersect).toBe(tiq.test(h, e).intersect);
        });
    });

    it('a degenerate ellipsoid is not supported (zero extents divide by 0)',
        () => {
            // Hyperellipsoid requires positive extents for GetMInverse to be
            // meaningful; a zero extent collapses the ellipsoid to a disk and
            // the support radius along that axis is zero.
            const e = ellipsoid(vec(0, 0, 0), unitAxes3(), vec(2, 3, 1e-8));
            // The plane z = 0.5 misses the (essentially flat) ellipsoid only
            // if its support along z is below 0.5.
            expect(tiq.test(halfspace([0, 0, 1], 0.5), e).intersect)
                .toBe(false);
            expect(tiq.test(halfspace([0, 0, 1], -0.5), e).intersect)
                .toBe(true);
            expect(tiq.test(halfspace([1, 0, 0], 1.5), e).intersect)
                .toBe(true);
            expect(tiq.test(halfspace([1, 0, 0], 2.5), e).intersect)
                .toBe(false);
        });

    it('rejects mismatched dimensions', () => {
        const e2 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(() => tiq.test(halfspace([0, 0, 1], 0), e2)).toThrow();
    });
});
