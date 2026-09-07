import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import {
    IntrEllipsoid3Ellipsoid3TI,
    IntrEllipsoid3Ellipsoid3Classification as C
} from '../src/IntrEllipsoid3Ellipsoid3.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

function ellipsoid(center: number[], axes: Vector[], extent: number[]):
    Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(Vector.fromArray(center), axes,
        Vector.fromArray(extent));
}

function sphere(center: number[], radius: number): Hyperellipsoid {
    return ellipsoid(center, unitAxes, [radius, radius, radius]);
}

function orthonormalFrame(w: Vector): Vector[] {
    const v = [w.clone(), Vector.zero(3), Vector.zero(3)];
    normalize(v[0]);
    computeOrthogonalComplement3(1, v, false);
    return [v[0], v[1], v[2]];
}

// (X-C)^T*M*(X-C) - 1; negative strictly inside.
function quadratic(e: Hyperellipsoid, X: Vector): number {
    const diff = sub(X, e.center);
    let sum = 0;
    for (let d = 0; d < 3; ++d) {
        const t = dot(diff, e.axis[d]) / e.extent.values[d];
        sum += t * t;
    }
    return sum - 1;
}

// Sample the surface of 'e' and return the extreme values of the quadratic of
// 'other' over those samples.
function extremesOnSurface(e: Hyperellipsoid, other: Hyperellipsoid):
    { min: number, max: number } {
    let minValue = Number.MAX_VALUE, maxValue = -Number.MAX_VALUE;
    const numTheta = 90, numPhi = 90;
    for (let i = 0; i < numTheta; ++i) {
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
            const value = quadratic(other, X);
            if (value < minValue) {
                minValue = value;
            }
            if (value > maxValue) {
                maxValue = value;
            }
        }
    }
    return { min: minValue, max: maxValue };
}

describe('IntrEllipsoid3Ellipsoid3TI', () => {
    const query = new IntrEllipsoid3Ellipsoid3TI();

    it('classifies concentric spheres', () => {
        const big = sphere([0, 0, 0], 2);
        const small = sphere([0, 0, 0], 1);
        expect(query.test(big, small)).toEqual({
            intersect: true,
            classification: C.ELLIPSOID0_CONTAINS_ELLIPSOID1
        });
        expect(query.test(small, big)).toEqual({
            intersect: true,
            classification: C.ELLIPSOID1_CONTAINS_ELLIPSOID0
        });
        expect(query.test(big, big)).toEqual({
            intersect: true,
            classification: C.ELLIPSOIDS_INTERSECTING
        });
    });

    it('classifies concentric ellipsoids that cross', () => {
        // A long thin ellipsoid poking out of a sphere.
        const s = sphere([0, 0, 0], 1);
        const needle = ellipsoid([0, 0, 0], unitAxes, [3, 0.2, 0.2]);
        const result = query.test(s, needle);
        expect(result.intersect).toBe(true);
        expect(result.classification).toBe(C.ELLIPSOIDS_INTERSECTING);
    });

    it('classifies separated spheres', () => {
        const a = sphere([0, 0, 0], 1);
        const b = sphere([5, 0, 0], 1);
        const result = query.test(a, b);
        expect(result.intersect).toBe(false);
        expect(result.classification).toBe(C.ELLIPSOIDS_SEPARATED);
    });

    it('classifies overlapping spheres', () => {
        const a = sphere([0, 0, 0], 1);
        const b = sphere([1.5, 0, 0], 1);
        const result = query.test(a, b);
        expect(result.intersect).toBe(true);
        expect(result.classification).toBe(C.ELLIPSOIDS_INTERSECTING);
    });

    it('classifies off-center containment', () => {
        const outer = sphere([0, 0, 0], 5);
        const inner = sphere([1, 0, 0], 1);
        expect(query.test(outer, inner).classification)
            .toBe(C.ELLIPSOID0_CONTAINS_ELLIPSOID1);
        expect(query.test(inner, outer).classification)
            .toBe(C.ELLIPSOID1_CONTAINS_ELLIPSOID0);
    });

    it('classifies containment of a rotated ellipsoid', () => {
        const axes = orthonormalFrame(vec(1, 1, 1));
        const outer = ellipsoid([0, 0, 0], unitAxes, [6, 6, 6]);
        const inner = ellipsoid([0.5, -0.5, 1], axes, [2, 1, 0.5]);
        expect(query.test(outer, inner).classification)
            .toBe(C.ELLIPSOID0_CONTAINS_ELLIPSOID1);
        expect(query.test(inner, outer).classification)
            .toBe(C.ELLIPSOID1_CONTAINS_ELLIPSOID0);
    });

    it('exercises the one-, two- and three-term root solvers', () => {
        // ellipsoid0 is the unit sphere, so the transformed problem has
        // D = diag(1/9,1/4,1) and K equal to the center of ellipsoid1. The
        // number of nonzero K components selects the root solver.
        const e0 = sphere([0, 0, 0], 1);
        const extents: [number, number, number] = [3, 2, 1];

        // One nonzero component: the single-term solver.
        const oneTerm = ellipsoid([4, 0, 0], unitAxes, extents);
        expect(query.test(e0, oneTerm).classification)
            .toBe(C.ELLIPSOIDS_INTERSECTING);
        const oneTermFar = ellipsoid([6, 0, 0], unitAxes, extents);
        expect(query.test(e0, oneTermFar).classification)
            .toBe(C.ELLIPSOIDS_SEPARATED);

        // Two nonzero components: the two-term solver.
        const twoTerm = ellipsoid([2, 1.5, 0], unitAxes, extents);
        const twoResult = query.test(e0, twoTerm);
        const twoExpected = extremesOnSurface(twoTerm, e0);
        expect(twoResult.intersect).toBe(twoExpected.min < 0);

        // Three nonzero components: the three-term solver.
        const threeTerm = ellipsoid([2, 1.5, 0.5], unitAxes, extents);
        const threeResult = query.test(e0, threeTerm);
        const threeExpected = extremesOnSurface(threeTerm, e0);
        expect(threeResult.intersect).toBe(threeExpected.min < 0);
    });

    it('rejects non-3D ellipsoids', () => {
        const e2 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(() => query.test(e2, e2)).toThrow();
    });

    it('agrees with the analytic sphere-sphere classification', () => {
        let seed = 1357911;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 150; ++trial) {
            const r0 = 0.3 + rand() * 2;
            const r1 = 0.3 + rand() * 2;
            const c1 = vec(rand() * 6 - 3, rand() * 6 - 3, rand() * 6 - 3);
            const distance = Math.sqrt(dot(c1, c1));
            // Skip configurations too close to a classification boundary.
            if (Math.abs(distance - (r0 + r1)) < 1e-2 ||
                Math.abs(distance + r1 - r0) < 1e-2 ||
                Math.abs(distance + r0 - r1) < 1e-2 ||
                distance < 1e-6) {
                continue;
            }

            const e0 = sphere([0, 0, 0], r0);
            const e1 = sphere([c1.values[0], c1.values[1], c1.values[2]], r1);
            const result = query.test(e0, e1);

            let expected: C;
            if (distance + r1 <= r0) {
                expected = C.ELLIPSOID0_CONTAINS_ELLIPSOID1;
            } else if (distance + r0 <= r1) {
                expected = C.ELLIPSOID1_CONTAINS_ELLIPSOID0;
            } else if (distance <= r0 + r1) {
                expected = C.ELLIPSOIDS_INTERSECTING;
            } else {
                expected = C.ELLIPSOIDS_SEPARATED;
            }
            expect(result.classification).toBe(expected);
            expect(result.intersect)
                .toBe(expected !== C.ELLIPSOIDS_SEPARATED);
        }
    });

    it('agrees with brute-force sampling for general ellipsoids', () => {
        let seed = 24681012;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        let numTested = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const axes0 = orthonormalFrame(
                vec(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1));
            const axes1 = orthonormalFrame(
                vec(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1));
            const e0 = ellipsoid([0, 0, 0], axes0,
                [0.5 + rand() * 2, 0.5 + rand() * 2, 0.5 + rand() * 2]);
            const e1 = ellipsoid(
                [rand() * 4 - 2, rand() * 4 - 2, rand() * 4 - 2], axes1,
                [0.4 + rand() * 1.5, 0.4 + rand() * 1.5, 0.4 + rand() * 1.5]);

            const on1 = extremesOnSurface(e1, e0);
            const on0 = extremesOnSurface(e0, e1);

            // Skip configurations that are numerically ambiguous.
            const tol = 5e-2;
            if (Math.abs(on1.max) < tol || Math.abs(on1.min) < tol ||
                Math.abs(on0.max) < tol) {
                continue;
            }

            let expected: C;
            if (on1.max < 0) {
                expected = C.ELLIPSOID0_CONTAINS_ELLIPSOID1;
            } else if (on0.max < 0) {
                expected = C.ELLIPSOID1_CONTAINS_ELLIPSOID0;
            } else if (on1.min < 0) {
                expected = C.ELLIPSOIDS_INTERSECTING;
            } else {
                expected = C.ELLIPSOIDS_SEPARATED;
            }

            const result = query.test(e0, e1);
            expect(result.classification).toBe(expected);
            expect(result.intersect)
                .toBe(expected !== C.ELLIPSOIDS_SEPARATED);
            ++numTested;
        }
        expect(numTested).toBeGreaterThan(20);
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrEllipsoid3Ellipsoid3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, seededRandom, wellScaled
} from './helpers/arbitraries.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no frame component is a
// subnormal that would underflow when squared.
const angleQ = () => wellScaled(-Math.PI, Math.PI);
const extentQ = () => fc.double(
    { min: 0.4, max: 2.5, noNaN: true, noDefaultInfinity: true });

// R = Rz(a)*Ry(b)*Rx(c); the columns are the ellipsoid axes.
function rotFrameQ(a: number, b: number, c: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    const cc = Math.cos(c), sc = Math.sin(c);
    return [
        vec(ca * cb, sa * cb, -sb),
        vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc, cb * sc),
        vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc, cb * cc)];
}

const ellipsoidArb = (span: number) => fc.tuple(
    wellScaled(-span, span), wellScaled(-span, span), wellScaled(-span, span),
    angleQ(), angleQ(), angleQ(), extentQ(), extentQ(), extentQ()
).map(([cx, cy, cz, a, b, c, e0, e1, e2]) =>
    Hyperellipsoid.fromCenterAxisExtent(vec(cx, cy, cz),
        rotFrameQ(a, b, c), vec(e0, e1, e2)));

const ellipsoidPair = fc.tuple(ellipsoidArb(2.5), ellipsoidArb(2.5))
    .map(([e0, e1]) => ({ e0, e1 }));

// A dense sampling of the surface of 'e'.
function surfaceSamples(e: Hyperellipsoid, n = 96): Vector[] {
    const out: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        const th = (2 * Math.PI * i) / n;
        for (let j = 0; j <= n; ++j) {
            const ph = (Math.PI * j) / n;
            const u = [Math.sin(ph) * Math.cos(th), Math.sin(ph) * Math.sin(th),
                Math.cos(ph)];
            let x = e.center.clone();
            for (let d = 0; d < 3; ++d) {
                x = add(x, mul(e.extent.values[d] * u[d], e.axis[d]));
            }
            out.push(x);
        }
    }
    return out;
}

// The independent classification: sample the surface of e1 and look at the
// sign of the quadratic of e0 there, then break the "all outside" tie by
// asking whether the center of e0 is inside e1.
//   max < 0            -> e0 contains e1
//   min > 0, C0 inside e1 -> e1 contains e0
//   min > 0, otherwise -> separated
//   min < 0 < max      -> intersecting
// 'margin' is the smallest |value| seen, which says how close the decision is
// to a tie.
function referenceClassify(e0: Hyperellipsoid, e1: Hyperellipsoid):
    { kind: C, margin: number } {
    let lo = Number.MAX_VALUE, hi = -Number.MAX_VALUE;
    for (const x of surfaceSamples(e1)) {
        const q = quadratic(e0, x);
        if (q < lo) { lo = q; }
        if (q > hi) { hi = q; }
    }
    const margin = Math.min(Math.abs(lo), Math.abs(hi));
    if (hi < 0) { return { kind: C.ELLIPSOID0_CONTAINS_ELLIPSOID1, margin }; }
    if (lo > 0) {
        const inside = quadratic(e1, e0.center) < 0;
        return {
            kind: inside ? C.ELLIPSOID1_CONTAINS_ELLIPSOID0
                : C.ELLIPSOIDS_SEPARATED,
            margin: Math.min(margin, Math.abs(quadratic(e1, e0.center)))
        };
    }
    return { kind: C.ELLIPSOIDS_INTERSECTING, margin };
}

describe('IntrEllipsoid3Ellipsoid3 verification', () => {
    const query = new IntrEllipsoid3Ellipsoid3TI();

    // Upstream's GetRoots isolates the roots of f(s) with brackets built from
    // an ad-hoc epsilon = 0.001 (its own comment asks "What role does epsilon
    // play?") and guards them with LogAssert on the sign of F at the bracket
    // endpoints. Those asserts are reachable for perfectly ordinary input --
    // see the pinned regression below -- so the randomized properties skip a
    // case that throws instead of pretending the query answered it.
    function tryTest(a: Hyperellipsoid, b: Hyperellipsoid):
        { intersect: boolean, classification: C } | null {
        try {
            return query.test(a, b);
        }
        catch {
            return null;
        }
    }

    it('pins the reachable GetRoots bracketing assert', () => {
        // Two ordinary ellipsoids: a ball of radius 0.4 at the origin and an
        // ellipsoid of extents (0.49, 0.41, 0.91) whose center is 0.001 away
        // and whose frame is a near-180-degree rotation about x. Upstream
        // (and therefore the port) throws instead of classifying. The failure
        // is erratic in the center offset: 0.002 and 0.0001 classify, 0.001
        // and 1e-6 throw.
        const e0 = ellipsoid([0, 0, 0], unitAxes, [0.4, 0.4, 0.4]);
        const frame = [
            vec(0.9999995000000417, 0, -0.0009999998333333417),
            vec(8.586571855778066e-11, -0.9999999999999963,
                8.586568993587257e-8),
            vec(-0.000999999833333338, -8.586573286873543e-8,
                -0.999999500000038)];
        const e1 = ellipsoid([0.001, 0, 0], frame,
            [0.4904589041303176, 0.41139783664312857, 0.9127878433790204]);
        expect(() => query.test(e0, e1)).toThrow('Unexpected condition.');
    });

    it('matches a brute-force surface sampling classification', () => {
        const rnd = seededRandom(0x2b71fa9);
        let checked = 0;
        for (let trial = 0; trial < 120; ++trial) {
            const mk = (span: number): Hyperellipsoid =>
                Hyperellipsoid.fromCenterAxisExtent(
                    vec(rnd() * 2 * span - span, rnd() * 2 * span - span,
                        rnd() * 2 * span - span),
                    rotFrameQ(rnd() * 6, rnd() * 6, rnd() * 6),
                    vec(0.4 + rnd() * 2, 0.4 + rnd() * 2, 0.4 + rnd() * 2));
            const e0 = mk(2), e1 = mk(2);
            const ref = referenceClassify(e0, e1);
            // The sampled extremes are only as accurate as the 96x97 grid, so
            // a decision that hinges on a value within 1e-2 of zero is a
            // near-tangency the sampling cannot resolve.
            if (ref.margin < 1e-2) { continue; }
            const r = tryTest(e0, e1);
            if (r === null) { continue; }
            expect(r.classification).toBe(ref.kind);
            expect(r.intersect).toBe(ref.kind !== C.ELLIPSOIDS_SEPARATED);
            ++checked;
        }
        expect(checked).toBeGreaterThan(60);
    }, 60000);

    it('never reports INVALID and keeps intersect consistent', () => {
        check(ellipsoidPair, ({ e0, e1 }) => {
            const r = tryTest(e0, e1);
            if (r === null) { return; }
            expect(r.classification).not.toBe(C.INVALID);
            expect(r.intersect)
                .toBe(r.classification !== C.ELLIPSOIDS_SEPARATED);
        }, 60);
    });

    it('classifies a swapped pair with the containment roles exchanged',
        () => {
            const rnd = seededRandom(0x515c0de);
            for (let trial = 0; trial < 120; ++trial) {
                const mk = (): Hyperellipsoid =>
                    Hyperellipsoid.fromCenterAxisExtent(
                        vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                        rotFrameQ(rnd() * 6, rnd() * 6, rnd() * 6),
                        vec(0.4 + rnd() * 2, 0.4 + rnd() * 2, 0.4 + rnd() * 2));
                const e0 = mk(), e1 = mk();
                const ref = referenceClassify(e0, e1);
                if (ref.margin < 1e-2) { continue; }
                const ra = tryTest(e0, e1), rb = tryTest(e1, e0);
                if (ra === null || rb === null) { continue; }
                const a = ra.classification;
                const b = rb.classification;
                const flip = (k: C): C =>
                    k === C.ELLIPSOID0_CONTAINS_ELLIPSOID1
                        ? C.ELLIPSOID1_CONTAINS_ELLIPSOID0
                        : (k === C.ELLIPSOID1_CONTAINS_ELLIPSOID0
                            ? C.ELLIPSOID0_CONTAINS_ELLIPSOID1 : k);
                expect(b).toBe(flip(a));
            }
        }, 60000);

    it('is equivariant under a rigid motion', () => {
        const rnd = seededRandom(0x77e1a3);
        for (let trial = 0; trial < 150; ++trial) {
            const mk = (): Hyperellipsoid =>
                Hyperellipsoid.fromCenterAxisExtent(
                    vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                    rotFrameQ(rnd() * 6, rnd() * 6, rnd() * 6),
                    vec(0.4 + rnd() * 2, 0.4 + rnd() * 2, 0.4 + rnd() * 2));
            const e0 = mk(), e1 = mk();
            const fr = rotFrameQ(rnd() * 6, rnd() * 6, rnd() * 6);
            const t = vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2);
            const rot = (v: Vector): Vector => add(mul(v.get(0), fr[0]),
                add(mul(v.get(1), fr[1]), mul(v.get(2), fr[2])));
            const xf = (e: Hyperellipsoid): Hyperellipsoid =>
                Hyperellipsoid.fromCenterAxisExtent(add(rot(e.center), t),
                    [rot(e.axis[0]), rot(e.axis[1]), rot(e.axis[2])],
                    e.extent);
            const a = tryTest(e0, e1);
            const b = tryTest(xf(e0), xf(e1));
            if (a === null || b === null) { continue; }
            expect(b.intersect).toBe(a.intersect);
            expect(b.classification).toBe(a.classification);
        }
    }, 30000);

    it('is equivariant under a uniform scale about the origin', () => {
        const rnd = seededRandom(0x9c02f1);
        for (let trial = 0; trial < 150; ++trial) {
            const mk = (): Hyperellipsoid =>
                Hyperellipsoid.fromCenterAxisExtent(
                    vec(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                    rotFrameQ(rnd() * 6, rnd() * 6, rnd() * 6),
                    vec(0.4 + rnd() * 2, 0.4 + rnd() * 2, 0.4 + rnd() * 2));
            const e0 = mk(), e1 = mk();
            const s = 0.25 + 3 * rnd();
            const scale = (e: Hyperellipsoid): Hyperellipsoid =>
                Hyperellipsoid.fromCenterAxisExtent(mul(s, e.center),
                    [e.axis[0], e.axis[1], e.axis[2]], mul(s, e.extent));
            const a = tryTest(e0, e1);
            const b = tryTest(scale(e0), scale(e1));
            if (a === null || b === null) { continue; }
            expect(b.intersect).toBe(a.intersect);
            expect(b.classification).toBe(a.classification);
        }
    }, 30000);

    it('never separates ellipsoids that share a center', () => {
        check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
            angleQ(), angleQ(), angleQ(), extentQ(), extentQ(), extentQ(),
            angleQ(), angleQ(), angleQ(), extentQ(), extentQ(), extentQ()),
        ([cx, cy, cz, a0, b0, c0, x0, y0, z0, a1, b1, c1, x1, y1, z1]) => {
            const centre = vec(cx, cy, cz);
            const e0 = Hyperellipsoid.fromCenterAxisExtent(centre,
                rotFrameQ(a0, b0, c0), vec(x0, y0, z0));
            const e1 = Hyperellipsoid.fromCenterAxisExtent(centre.clone(),
                rotFrameQ(a1, b1, c1), vec(x1, y1, z1));
            const r = tryTest(e0, e1);
            if (r === null) { return; }
            expect(r.intersect).toBe(true);
            expect(r.classification).not.toBe(C.ELLIPSOIDS_SEPARATED);
            expect(r.classification).not.toBe(C.INVALID);
        }, 100);
    });

    it('never separates two copies of the same ellipsoid', () => {
        // Coincident centers take the K == 0 branch, where the decision is
        // maxSqrDistance < 1 versus minSqrDistance > 1 with every eigenvalue
        // of M2 equal to 1 in exact arithmetic. Round-off in the eigensolver
        // puts them a few ulps either side of 1, so the exact-coincidence
        // answer is a knife edge between INTERSECTING and a containment; only
        // "not separated" is stable.
        check(ellipsoidArb(2), e => {
            const r = tryTest(e, e.clone());
            if (r === null) { return; }
            expect(r.intersect).toBe(true);
            expect(r.classification).not.toBe(C.ELLIPSOIDS_SEPARATED);
            expect(r.classification).not.toBe(C.INVALID);
        }, 100);
    });

    it('separates ellipsoids moved beyond their combined reach', () => {
        check(fc.tuple(ellipsoidArb(1), ellipsoidArb(1)), ([e0, e1]) => {
            const far = Hyperellipsoid.fromCenterAxisExtent(
                add(e1.center, vec(100, 0, 0)),
                [e1.axis[0], e1.axis[1], e1.axis[2]], e1.extent);
            const r = tryTest(e0, far);
            if (r === null) { return; }
            expect(r.intersect).toBe(false);
            expect(r.classification).toBe(C.ELLIPSOIDS_SEPARATED);
        }, 100);
    });

    it('rejects non-3D ellipsoids', () => {
        const e2 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(() => query.test(e2, e2)).toThrow();
    });
});
