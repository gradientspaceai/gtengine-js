import { describe, it, expect } from 'vitest';
import { Cone } from '../src/Cone.js';
import {
    IntrLine3Cone3FI,
    IntrLine3Cone3FIResultType,
    intrLine3Cone3Convert,
    intrLine3Cone3ConvertPoint
} from '../src/IntrLine3Cone3.js';
import {
    IntrRay3Cone3FI,
    defaultIntrRay3Cone3FIResult
} from '../src/IntrRay3Cone3.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { QFNumber } from '../src/QFNumber.js';
import { cross } from '../src/Vector3.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, rotationFrame, scaled,
    wellScaledVector
} from './helpers/arbitraries.js';

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

// A cone with apex at 'origin', axis 'direction' and the given half-angle.
// A negative maxHeight means the cone is infinite.
function cone(origin: number[], direction: number[], angle: number,
    minHeight: number, maxHeight: number): Cone {
    const C = new Cone(3);
    const d = Vector.fromArray(direction);
    normalize(d);
    C.ray = Ray.fromOriginDirection(Vector.fromArray(origin), d);
    C.setAngle(angle);
    if (maxHeight < 0) {
        if (minHeight > 0) {
            C.makeInfiniteTruncatedCone(minHeight);
        }
        else {
            C.makeInfiniteCone();
        }
    }
    else if (minHeight > 0) {
        C.makeConeFrustum(minHeight, maxHeight);
    }
    else {
        C.makeFiniteCone(maxHeight);
    }
    return C;
}

// True when X is in the solid cone (with a small tolerance).
function inSolidCone(C: Cone, X: Vector, tolerance: number): boolean {
    const diff = sub(X, C.ray.origin);
    const h = dot(C.ray.direction, diff);
    if (h < C.getMinHeight() - tolerance) {
        return false;
    }
    if (C.isFinite() && h > C.getMaxHeight() + tolerance) {
        return false;
    }
    return h >= length(diff) * C.cosAngle - tolerance;
}

const fi = new IntrRay3Cone3FI();
const lineFI = new IntrLine3Cone3FI();

const T = IntrLine3Cone3FIResultType;
const quarterPi = Math.PI / 4;

describe('IntrRay3Cone3', () => {
    it('default-constructs the result as empty', () => {
        const r = defaultIntrRay3Cone3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.type).toBe(T.isEmpty);
    });

    it('finds the known segment for a ray crossing a finite cone', () => {
        // Cone with apex at the origin, axis +z, half-angle pi/4, heights
        // [0,4]. The ray starts at (-5,0,2) heading in +x; the cone cross
        // section at z = 2 is the disk of radius 2, so the segment endpoints
        // are (-2,0,2) and (2,0,2).
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const R = ray([-5, 0, 2], [1, 0, 0]);
        const result = fi.find(R, C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.isSegment);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        expect(P0.values[0]).toBeCloseTo(-2, 10);
        expect(P1.values[0]).toBeCloseTo(2, 10);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(3, 10);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(7, 10);
    });

    it('clips the line result at the ray origin (block 14)', () => {
        // The same configuration but the ray origin is inside the cone, so
        // the line segment [3,7] is clipped to [0,7] in ray parameters.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const R = ray([0, 0, 2], [1, 0, 0]);
        const result = fi.find(R, C);
        expect(result.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(2, 10);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        expect(P0.values[0]).toBeCloseTo(0, 12);
    });

    it('reports empty when the intersection is entirely behind the origin', () => {
        // Block 15: the whole line-cone segment has t < 0.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const R = ray([5, 0, 2], [1, 0, 0]);
        const result = fi.find(R, C);
        expect(result.intersect).toBe(false);
        expect(result.type).toBe(T.isEmpty);
    });

    it('clips to a sliver when the origin is just inside the boundary', () => {
        // The cone cross section at z = 2 is the disk of radius 2, so a ray
        // starting just inside the boundary and heading out exits almost
        // immediately. An exactly-on-the-boundary origin (upstream block 16)
        // is not representable in floating point for this cone.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const inside = fi.find(ray([2 - 1e-6, 0, 2], [1, 0, 0]), C);
        expect(inside.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(inside.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(inside.t[1])).toBeCloseTo(1e-6, 12);

        const outside = fi.find(ray([2 + 1e-6, 0, 2], [1, 0, 0]), C);
        expect(outside.intersect).toBe(false);
        expect(outside.type).toBe(T.isEmpty);
    });

    it('keeps a positive ray for an infinite cone (block 17)', () => {
        // The ray starts on the axis inside an infinite cone, so the whole
        // ray is inside.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);
        const R = ray([0, 0, 1], [0, 0, 1]);
        const result = fi.find(R, C);
        expect(result.type).toBe(T.isRayPositive);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        expect(P0.values[2]).toBeCloseTo(1, 12);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(0, 12);
    });

    it('turns a negative ray into a segment or a point (blocks 18-20)', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);

        // Block 18: the ray travels down the axis from inside the cone, so
        // it exits at the apex.
        const R18 = ray([0, 0, 1], [0, 0, -1]);
        const r18 = fi.find(R18, C);
        expect(r18.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(r18.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(r18.t[1])).toBeCloseTo(1, 10);
        const Q1 = intrLine3Cone3ConvertPoint(r18.P[1]);
        expect(length(Q1)).toBeCloseTo(0, 10);

        // Block 19: below the apex heading further down, nothing is hit.
        const R19 = ray([0, 0, -1], [0, 0, -1]);
        expect(fi.find(R19, C).type).toBe(T.isEmpty);

        // Block 20: starting exactly at the apex heading down touches only
        // the apex.
        const R20 = ray([0, 0, 0], [0, 0, -1]);
        const r20 = fi.find(R20, C);
        expect(r20.type).toBe(T.isPoint);
        const P0 = intrLine3Cone3ConvertPoint(r20.P[0]);
        expect(length(P0)).toBeCloseTo(0, 12);
    });

    it('keeps a ray that starts off-axis inside an infinite cone', () => {
        // The line x = 0.5, y = 0 enters the 45-degree cone at z = 0.5, so a
        // downward ray from (0.5,0,3) is clipped to the segment ending there.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);
        const down = fi.find(ray([0.5, 0, 3], [0, 0, -1]), C);
        expect(down.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(down.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(down.t[1])).toBeCloseTo(2.5, 8);

        // The same line traversed upward from inside is an unbounded ray.
        const up = fi.find(ray([0.5, 0, 3], [0, 0, 1]), C);
        expect(up.type).toBe(T.isRayPositive);
        expect(intrLine3Cone3Convert(up.t[0])).toBeCloseTo(0, 12);
    });

    it('clamps to the cone height range for a cone frustum', () => {
        // A frustum with heights [1,3]. A ray down the axis from above must
        // produce the segment between the caps.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 1, 3);
        const R = ray([0, 0, 10], [0, 0, -1]);
        const result = fi.find(R, C);
        expect(result.type).toBe(T.isSegment);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        const zs = [P0.values[2], P1.values[2]].sort((a, b) => a - b);
        expect(zs[0]).toBeCloseTo(1, 10);
        expect(zs[1]).toBeCloseTo(3, 10);
    });

    it('agrees with the line query on the ray portion (randomized)', () => {
        let seed = 13572468;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (): number => 6 * rand() - 3;

        const cones = [
            cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1),
            cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 3),
            cone([0.5, -0.25, 0], [1, 1, 1], 0.6, 0.5, 2.5),
            cone([0, 0, 0], [0, 1, 0], 0.4, 0.25, -1)
        ];

        let numSegments = 0;
        for (let trial = 0; trial < 2000; ++trial) {
            const C = cones[trial % cones.length];
            const R = ray([rnd(), rnd(), rnd()], [rnd(), rnd(), rnd()]);
            const result = fi.find(R, C);
            const lineResult = lineFI.find(
                Line.fromOriginDirection(R.origin, R.direction), C);

            if (result.intersect) {
                // The ray result is never larger than the line result.
                expect(lineResult.intersect).toBe(true);
            }

            if (result.type === T.isSegment) {
                ++numSegments;
                const t0 = intrLine3Cone3Convert(result.t[0]);
                const t1 = intrLine3Cone3Convert(result.t[1]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                expect(t1).toBeGreaterThanOrEqual(t0 - 1e-12);
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
                expect(inSolidCone(C, P0, 1e-8)).toBe(true);
                expect(inSolidCone(C, P1, 1e-8)).toBe(true);
                // The midpoint is inside the solid cone too (convexity).
                const mid = mul(0.5, add(P0, P1));
                expect(inSolidCone(C, mid, 1e-8)).toBe(true);
                // The endpoints are on the ray.
                for (const [t, P] of [[t0, P0], [t1, P1]] as [number, Vector][]) {
                    const X = add(R.origin, mul(t, R.direction));
                    for (let i = 0; i < 3; ++i) {
                        expect(P.values[i]).toBeCloseTo(X.values[i], 8);
                    }
                }
            }
            else if (result.type === T.isPoint) {
                const t0 = intrLine3Cone3Convert(result.t[0]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                expect(inSolidCone(C, P0, 1e-8)).toBe(true);
            }
            else if (result.type === T.isRayPositive) {
                const t0 = intrLine3Cone3Convert(result.t[0]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                expect(inSolidCone(C, P0, 1e-8)).toBe(true);
                // Walking along the ray direction stays inside.
                const far = add(P0, mul(100, R.direction));
                expect(inSolidCone(C, far, 1e-6)).toBe(true);
            }
            else {
                expect(result.type).toBe(T.isEmpty);
            }
        }
        expect(numSegments).toBeGreaterThan(50);
    });
});

// ---------------------------------------------------------------------------
// Verification (V35): the ray query must be the line query clipped to
// [0,+infinity), and the reported interval must be exactly the set of ray
// parameters whose point is in the solid cone.
// ---------------------------------------------------------------------------

describe('IntrRay3Cone3 verification', () => {
    const fiv = new IntrRay3Cone3FI();
    const lineFIv = new IntrLine3Cone3FI();

    // A cone of one of the four kinds: infinite, infinite truncated, finite
    // and frustum. Coordinates are well scaled so that no dot product
    // underflows.
    // Directions and positions are drawn from a uniform grid rather than
    // fast-check's double(), whose bit-pattern-uniform sampling produces
    // near-zero components (and hence axis-aligned frames and degenerate
    // configurations) far more often than a uniform distribution would.
    const gridPoint = (lo: number, hi: number): fc.Arbitrary<Vector> =>
        fc.tuple(scaled(lo, hi), scaled(lo, hi), scaled(lo, hi))
            .map(a => Vector.fromArray(a));

    const gridDirection = gridPoint(-1, 1)
        .filter(v => length(v) > 0.3)
        .map(v => { const u = v.clone(); normalize(u); return u; });

    const coneArb = fc.tuple(gridPoint(-3, 3), gridDirection,
        fc.double({ min: 0.15, max: 1.35, noNaN: true }),
        fc.integer({ min: 0, max: 3 }),
        fc.double({ min: 0.25, max: 2, noNaN: true }),
        fc.double({ min: 0.25, max: 3, noNaN: true }))
        .map(([v, a, angle, kind, h0, dh]) => {
            const C = new Cone(3);
            C.ray = Ray.fromOriginDirection(v, a);
            C.setAngle(angle);
            if (kind === 0) { C.makeInfiniteCone(); }
            else if (kind === 1) { C.makeInfiniteTruncatedCone(h0); }
            else if (kind === 2) { C.makeFiniteCone(h0 + dh); }
            else { C.makeConeFrustum(h0, h0 + dh); }
            return C;
        });


    // A point of the solid cone, from its (height, radius, angle) parameters.
    function coneSample(C: Cone, s: number, u: number, phi: number): Vector {
        const hmin = C.getMinHeight();
        const hmax = C.isFinite() ? C.getMaxHeight() : hmin + 4;
        const h = hmin + s * (hmax - hmin);
        const r = u * h * C.tanAngle;
        const a = C.ray.direction;
        const k = Math.abs(a.values[0]) < 0.5 ? 0
            : (Math.abs(a.values[1]) < 0.5 ? 1 : 2);
        const e = Vector.unit(3, k);
        const q0 = sub(e, mul(dot(e, a), a));
        normalize(q0);
        const q1 = cross(a, q0);
        return add(add(C.ray.origin, mul(h, a)),
            add(mul(r * Math.cos(phi), q0), mul(r * Math.sin(phi), q1)));
    }

    // Half the rays are aimed at a point of the solid cone so that hits are
    // common; the rest are unconstrained.
    const coneAndRay = coneArb.chain(C => fc.tuple(fc.constant(C),
        gridPoint(-4, 4), gridDirection, fc.boolean(),
        scaled(0, 1), scaled(0, 1), scaled(0, 2 * Math.PI))
        .map(([cone, o, d, aim, s, u, phi]) => {
            if (!aim) {
                return [Ray.fromOriginDirection(o, d), cone] as [Ray, Cone];
            }
            const target = coneSample(cone, s, u, phi);
            const dir = sub(target, o);
            if (length(dir) < 1e-3) {
                return [Ray.fromOriginDirection(o, d), cone] as [Ray, Cone];
            }
            normalize(dir);
            return [Ray.fromOriginDirection(o, dir), cone] as [Ray, Cone];
        }));

    // The line-cone query solves c2*t^2 + 2*c1*t + c0 = 0 and branches on the
    // sign of the discriminant c1^2 - c0*c2 and on an exact test for the line
    // passing through the cone vertex. Both are formed by subtracting nearly
    // equal quantities, so when the discriminant is at the noise level (a line
    // tangent to the cone, or a line through the vertex such as one along the
    // cone axis) the branch taken is decided by rounding and the reported set
    // can be a point, or empty, where the true answer is a ray. Those
    // configurations are excluded from the properties below; see the PR notes
    // for the defect in the shared line-cone query.
    function wellConditioned(origin: Vector, dir: Vector, C: Cone): boolean {
        const PmV = sub(origin, C.ray.origin);
        const UdU = dot(dir, dir);
        const DdU = dot(C.ray.direction, dir);
        const DdPmV = dot(C.ray.direction, PmV);
        const UdPmV = dot(dir, PmV);
        const PmVdPmV = dot(PmV, PmV);
        const k = C.cosAngleSqr;
        const c2 = DdU * DdU - k * UdU;
        const c1 = DdU * DdPmV - k * UdPmV;
        const c0 = DdPmV * DdPmV - k * PmVdPmV;
        const discr = c1 * c1 - c0 * c2;
        const discrScale = Math.max(c1 * c1, Math.abs(c0 * c2));
        if (Math.abs(discr) < 1e-6 * discrScale) {
            return false;
        }
        if (Math.abs(c2) < 1e-6 * Math.max(DdU * DdU, k * UdU)) {
            return false;
        }
        // The height clamping against hmin and hmax is a second family of
        // knife edges; require the roots to be off the caps by a margin.
        const root = [(-c1 - Math.sqrt(Math.max(discr, 0))) / c2,
            (-c1 + Math.sqrt(Math.max(discr, 0))) / c2];
        for (const t of root) {
            const h = t * DdU + DdPmV;
            for (const cap of [0, C.getMinHeight(),
                C.isFinite() ? C.getMaxHeight() : Number.NaN]) {
                if (Number.isFinite(cap) && Math.abs(h - cap) < 1e-9) {
                    return false;
                }
            }
        }
        return true;
    }

    // The reported t-interval as a pair of numbers, with +/-infinity for the
    // unbounded types. The QFNumber sentinel for an infinite endpoint is
    // (+/-1, 0, d), which does not convert to a number, so it is replaced
    // here.
    function intervalOf(result: { type: number, t: [QFNumber, QFNumber] }):
        [number, number] | null {
        switch (result.type) {
            case T.isEmpty:
                return null;
            case T.isPoint: {
                const a = intrLine3Cone3Convert(result.t[0]);
                return [a, a];
            }
            case T.isSegment:
                return [intrLine3Cone3Convert(result.t[0]),
                    intrLine3Cone3Convert(result.t[1])];
            case T.isRayPositive:
                return [intrLine3Cone3Convert(result.t[0]),
                    Number.POSITIVE_INFINITY];
            default:  // isRayNegative
                return [Number.NEGATIVE_INFINITY,
                    intrLine3Cone3Convert(result.t[1])];
        }
    }

    it('is the line result clipped to t >= 0', () => {
        let numTested = 0, numNonEmpty = 0;
        check(coneAndRay, ([R, C]) => {
            if (!wellConditioned(R.origin, R.direction, C)) {
                return;
            }
            ++numTested;
            const lr = lineFIv.find(
                Line.fromOriginDirection(R.origin, R.direction), C);
            const li = intervalOf(lr);
            const rr = fiv.find(R, C);
            const ri = intervalOf(rr);

            if (li === null) {
                expect(ri).toBeNull();
                return;
            }
            // Skip the knife edges where an endpoint of the line interval is
            // numerically indistinguishable from the ray origin: there the
            // exact QFNumber comparison inside the query and the comparison
            // on the converted doubles can disagree.
            if (Math.abs(li[0]) < 1e-9 || Math.abs(li[1]) < 1e-9) {
                return;
            }
            const lo = Math.max(li[0], 0);
            const hi = li[1];
            if (hi < lo) {
                expect(ri).toBeNull();
                return;
            }
            expect(ri).not.toBeNull();
            const got = ri as [number, number];
            if (hi === Number.POSITIVE_INFINITY) {
                expect(rr.type).toBe(T.isRayPositive);
                expectClose(got[0], lo, 1e-9, 1e-9);
                return;
            }
            if (lo === hi) {
                expect(rr.type).toBe(T.isPoint);
            }
            expectClose(got[0], lo, 1e-9, 1e-9);
            expectClose(got[1], hi, 1e-9, 1e-9);
            ++numNonEmpty;
        });
        // The conditioning filter must not reject everything, and both the
        // empty and the nonempty outcome must occur.
        expect(numTested).toBeGreaterThan(100);
        expect(numNonEmpty).toBeGreaterThan(10);
        expect(numTested - numNonEmpty).toBeGreaterThan(10);
    });

    it('reports exactly the ray parameters whose point is in the cone', () => {
        check(coneAndRay, ([R, C]) => {
            if (!wellConditioned(R.origin, R.direction, C)) {
                return;
            }
            const rr = fiv.find(R, C);
            const ri = intervalOf(rr);
            const lo = ri === null ? 0 : ri[0];
            const hi = ri === null ? -1 : Math.min(ri[1], 12);
            const steps = 160;
            for (let i = 0; i <= steps; ++i) {
                const t = (i / steps) * 12;
                const X = add(R.origin, mul(t, R.direction));
                // A point comfortably inside the solid cone must be inside
                // the reported interval.
                if (inSolidCone(C, X, -1e-6)) {
                    expect(ri).not.toBeNull();
                    expect(t).toBeGreaterThan(lo - 1e-5);
                    expect(t).toBeLessThan(hi + 1e-5);
                }
                // A parameter comfortably inside the reported interval must
                // map to a point of the solid cone.
                if (ri !== null && t > lo + 1e-6 && t < hi - 1e-6) {
                    expect(inSolidCone(C, X, 1e-6)).toBe(true);
                }
            }
        }, 60);
    }, 30000);

    it('reports points that lie on the ray and in the cone', () => {
        check(coneAndRay, ([R, C]) => {
            if (!wellConditioned(R.origin, R.direction, C)) {
                return;
            }
            const rr = fiv.find(R, C);
            if (!rr.intersect) {
                return;
            }
            const P0 = intrLine3Cone3ConvertPoint(rr.P[0]);
            expect(inSolidCone(C, P0, 1e-7)).toBe(true);
            const t0 = intrLine3Cone3Convert(rr.t[0]);
            expect(t0).toBeGreaterThan(-1e-12);
            expectVectorClose(P0, add(R.origin, mul(t0, R.direction)),
                1e-8, 1e-8);
            if (rr.type === T.isSegment || rr.type === T.isPoint) {
                const P1 = intrLine3Cone3ConvertPoint(rr.P[1]);
                const t1 = intrLine3Cone3Convert(rr.t[1]);
                expect(inSolidCone(C, P1, 1e-7)).toBe(true);
                expectVectorClose(P1, add(R.origin, mul(t1, R.direction)),
                    1e-8, 1e-8);
                expect(t1).toBeGreaterThan(t0 - 1e-12);
            }
            else if (rr.type === T.isRayPositive) {
                // P[1] carries the ray direction rather than a point.
                const D = intrLine3Cone3ConvertPoint(rr.P[1]);
                expectVectorClose(D, R.direction, 1e-12, 1e-12);
            }
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(coneAndRay, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([[R, C], frame, tr]) => {
                if (!wellConditioned(R.origin, R.direction, C)) {
                    return;
                }
                const xfDir = (v: Vector): Vector => {
                    const w = new Vector(3);
                    for (let i = 0; i < 3; ++i) {
                        w.values[i] = frame[0].values[i] * v.values[0]
                            + frame[1].values[i] * v.values[1]
                            + frame[2].values[i] * v.values[2];
                    }
                    return w;
                };
                const xf = (v: Vector): Vector => add(xfDir(v), tr);
                const D = new Cone(3);
                D.ray = Ray.fromOriginDirection(xf(C.ray.origin),
                    xfDir(C.ray.direction));
                D.setAngle(C.angle);
                if (C.getMinHeight() > 0) {
                    if (C.isFinite()) {
                        D.makeConeFrustum(C.getMinHeight(), C.getMaxHeight());
                    }
                    else {
                        D.makeInfiniteTruncatedCone(C.getMinHeight());
                    }
                }
                else if (C.isFinite()) {
                    D.makeFiniteCone(C.getMaxHeight());
                }
                else {
                    D.makeInfiniteCone();
                }
                const r0 = fiv.find(R, C);
                const r1 = fiv.find(Ray.fromOriginDirection(xf(R.origin),
                    xfDir(R.direction)), D);
                const i0 = intervalOf(r0);
                const i1 = intervalOf(r1);
                if (i0 === null || i1 === null) {
                    // A grazing tangency can be classified either way after
                    // the motion; only assert when both are nonempty or the
                    // interval is not degenerate.
                    if (i0 !== null && i0[1] - i0[0] > 1e-6) {
                        expect(i1).not.toBeNull();
                    }
                    if (i1 !== null && i1[1] - i1[0] > 1e-6) {
                        expect(i0).not.toBeNull();
                    }
                    return;
                }
                expectClose(i0[0], i1[0], 1e-6, 1e-6);
                if (Number.isFinite(i0[1]) && Number.isFinite(i1[1])) {
                    expectClose(i0[1], i1[1], 1e-6, 1e-6);
                }
            });
    });
});
