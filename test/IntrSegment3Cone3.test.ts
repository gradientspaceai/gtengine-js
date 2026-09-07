import { describe, it, expect } from 'vitest';
import { Cone } from '../src/Cone.js';
import {
    IntrLine3Cone3FI,
    IntrLine3Cone3FIResultType,
    intrLine3Cone3Convert,
    intrLine3Cone3ConvertPoint
} from '../src/IntrLine3Cone3.js';
import {
    IntrSegment3Cone3FI,
    defaultIntrSegment3Cone3FIResult
} from '../src/IntrSegment3Cone3.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Segment } from '../src/Segment.js';
import { QFNumber } from '../src/QFNumber.js';
import { cross } from '../src/Vector3.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, scaled,
    wellScaledVector
} from './helpers/arbitraries.js';

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
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

const fi = new IntrSegment3Cone3FI();
const lineFI = new IntrLine3Cone3FI();

const T = IntrLine3Cone3FIResultType;
const quarterPi = Math.PI / 4;

describe('IntrSegment3Cone3', () => {
    it('default-constructs the result as empty', () => {
        const r = defaultIntrSegment3Cone3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.type).toBe(T.isEmpty);
    });

    it('finds the known chord of a finite cone', () => {
        // The cone cross section at z = 2 is the disk of radius 2. The
        // segment from (-5,0,2) to (5,0,2) is parameterized over [0,1], so
        // the chord is t in [0.3,0.7].
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const S = segment([-5, 0, 2], [5, 0, 2]);
        const result = fi.find(S, C);
        expect(result.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(0.3, 10);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(0.7, 10);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        expect(P0.values[0]).toBeCloseTo(-2, 10);
        expect(P1.values[0]).toBeCloseTo(2, 10);
        expect(P0.values[2]).toBeCloseTo(2, 12);
    });

    it('clips both ends to the segment interval (block 24)', () => {
        // The segment is strictly inside the chord, so the result is the
        // whole segment [0,1].
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const S = segment([-1, 0, 2], [1, 0, 2]);
        const result = fi.find(S, C);
        expect(result.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(1, 12);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        expect(P0.values[0]).toBeCloseTo(-1, 12);
        expect(P1.values[0]).toBeCloseTo(1, 12);
    });

    it('reports empty for segments on either side of the chord (block 23)', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        expect(fi.find(segment([-9, 0, 2], [-5, 0, 2]), C).type).toBe(T.isEmpty);
        expect(fi.find(segment([5, 0, 2], [9, 0, 2]), C).type).toBe(T.isEmpty);
    });

    it('reports empty for a segment entirely below the cone', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const result = fi.find(segment([-3, -3, -2], [3, 3, -2]), C);
        expect(result.intersect).toBe(false);
    });

    it('clips an unbounded line result to the segment (blocks 27 and 30)', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);

        // Block 30: the line result is a negative ray ending at the apex. The
        // segment runs downward from (0.5,0,3) past the entry point z = 0.5.
        const down = fi.find(segment([0.5, 0, 3], [0.5, 0, -3]), C);
        expect(down.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(down.t[0])).toBeCloseTo(0, 12);
        // The segment direction has length 6 and the entry is 2.5 below the
        // start, so t = 2.5/6.
        expect(intrLine3Cone3Convert(down.t[1])).toBeCloseTo(2.5 / 6, 8);

        // Block 27: the line result is a positive ray, and the segment runs
        // upward from inside, so the whole segment is the intersection.
        const up = fi.find(segment([0.5, 0, 3], [0.5, 0, 9]), C);
        expect(up.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(up.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(up.t[1])).toBeCloseTo(1, 12);
    });

    it('reports empty when an unbounded line result misses the segment', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);
        // Block 29: the negative-ray line result ends at the apex, well
        // before the segment starts.
        const result = fi.find(segment([0.5, 0, -1], [0.5, 0, -5]), C);
        expect(result.type).toBe(T.isEmpty);
    });

    it('touches the apex in a single point (block 31)', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);
        const result = fi.find(segment([0, 0, 0], [0, 0, -4]), C);
        expect(result.type).toBe(T.isPoint);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        expect(length(P0)).toBeCloseTo(0, 12);
    });

    it('clamps to the cone frustum height range', () => {
        // The frustum has heights [1,3]; a segment down the axis from above
        // must be clipped to the two cap planes.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 1, 3);
        const result = fi.find(segment([0, 0, 10], [0, 0, -10]), C);
        expect(result.type).toBe(T.isSegment);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        const zs = [P0.values[2], P1.values[2]].sort((a, b) => a - b);
        expect(zs[0]).toBeCloseTo(1, 10);
        expect(zs[1]).toBeCloseTo(3, 10);
    });

    it('agrees with the line query and stays in the segment (randomized)', () => {
        let seed = 55501234;
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
            const p0 = [rnd(), rnd(), rnd()];
            const p1 = [rnd(), rnd(), rnd()];
            const S = segment(p0, p1);
            const result = fi.find(S, C);

            // The result is never a ray: the segment is bounded.
            expect(result.type).not.toBe(T.isRayPositive);
            expect(result.type).not.toBe(T.isRayNegative);

            if (result.intersect) {
                // The line query must find at least as much.
                const dir = sub(S.p[1], S.p[0]);
                const lineResult = lineFI.find(
                    Line.fromOriginDirection(S.p[0], dir), C);
                expect(lineResult.intersect).toBe(true);
            }

            if (result.type === T.isSegment) {
                ++numSegments;
                const t0 = intrLine3Cone3Convert(result.t[0]);
                const t1 = intrLine3Cone3Convert(result.t[1]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                expect(t1).toBeLessThanOrEqual(1 + 1e-12);
                expect(t1).toBeGreaterThanOrEqual(t0 - 1e-12);
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
                expect(inSolidCone(C, P0, 1e-8)).toBe(true);
                expect(inSolidCone(C, P1, 1e-8)).toBe(true);
                expect(inSolidCone(C, mul(0.5, add(P0, P1)), 1e-8)).toBe(true);
                // The endpoints are on the segment.
                const dir = sub(S.p[1], S.p[0]);
                for (const [t, P] of [[t0, P0], [t1, P1]] as [number, Vector][]) {
                    const X = add(S.p[0], mul(t, dir));
                    for (let i = 0; i < 3; ++i) {
                        expect(P.values[i]).toBeCloseTo(X.values[i], 8);
                    }
                }
            }
            else if (result.type === T.isPoint) {
                const t0 = intrLine3Cone3Convert(result.t[0]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                expect(t0).toBeLessThanOrEqual(1 + 1e-12);
                expect(inSolidCone(C,
                    intrLine3Cone3ConvertPoint(result.P[0]), 1e-8)).toBe(true);
            }
            else {
                expect(result.type).toBe(T.isEmpty);
            }
        }
        expect(numSegments).toBeGreaterThan(30);
    });
});

// ---------------------------------------------------------------------------
// Verification (V35): the segment query must be the line query clipped to
// [0,1] in the parameterization p[0] + t*(p[1] - p[0]), and the reported
// interval must be exactly the set of parameters whose point is in the solid
// cone.
// ---------------------------------------------------------------------------

describe('IntrSegment3Cone3 verification', () => {
    const fiv = new IntrSegment3Cone3FI();
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

    // Half the segments have one endpoint in the solid cone so that hits are
    // common; the rest are unconstrained.
    const coneAndSegment = coneArb.chain(C => fc.tuple(fc.constant(C),
        gridPoint(-4, 4), gridPoint(-4, 4), fc.boolean(),
        scaled(0, 1), scaled(0, 1), scaled(0, 2 * Math.PI))
        .map(([cone, p0, p1, aim, s, u, phi]) => {
            const q1 = aim ? coneSample(cone, s, u, phi) : p1;
            if (length(sub(q1, p0)) < 0.5) {
                return [Segment.fromEndpoints(p0, add(p0,
                    Vector.fromArray([1, 0, 0]))), cone] as [Segment, Cone];
            }
            return [Segment.fromEndpoints(p0, q1), cone] as [Segment, Cone];
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

    it('is the line result clipped to [0,1]', () => {
        let numTested = 0, numNonEmpty = 0;
        check(coneAndSegment, ([S, C]) => {
            if (!wellConditioned(S.p[0], sub(S.p[1], S.p[0]), C)) {
                return;
            }
            ++numTested;
            // The query parameterizes the segment as p[0] + t*(p[1] - p[0]),
            // so the line query must use the same (non-unit) direction.
            const dir = sub(S.p[1], S.p[0]);
            const lr = lineFIv.find(
                Line.fromOriginDirection(S.p[0], dir), C);
            const li = intervalOf(lr);
            const sr = fiv.find(S, C);
            const si = intervalOf(sr);

            if (li === null) {
                expect(si).toBeNull();
                return;
            }
            // Skip the knife edges where a line-interval endpoint coincides
            // numerically with a segment endpoint.
            for (const v of [li[0], li[1]]) {
                if (Math.abs(v) < 1e-9 || Math.abs(v - 1) < 1e-9) {
                    return;
                }
            }
            const lo = Math.max(li[0], 0);
            const hi = Math.min(li[1], 1);
            if (hi < lo) {
                expect(si).toBeNull();
                return;
            }
            expect(si).not.toBeNull();
            const got = si as [number, number];
            expectClose(got[0], lo, 1e-9, 1e-9);
            expectClose(got[1], hi, 1e-9, 1e-9);
            expect(sr.type).toBe(lo === hi ? T.isPoint : T.isSegment);
            ++numNonEmpty;
        });
        // The conditioning filter must not reject everything, and both the
        // empty and the nonempty outcome must occur.
        expect(numTested).toBeGreaterThan(100);
        expect(numNonEmpty).toBeGreaterThan(10);
        expect(numTested - numNonEmpty).toBeGreaterThan(10);
    });

    it('reports exactly the parameters whose point is in the cone', () => {
        check(coneAndSegment, ([S, C]) => {
            if (!wellConditioned(S.p[0], sub(S.p[1], S.p[0]), C)) {
                return;
            }
            const dir = sub(S.p[1], S.p[0]);
            const sr = fiv.find(S, C);
            const si = intervalOf(sr);
            const lo = si === null ? 0 : si[0];
            const hi = si === null ? -1 : si[1];
            const scale = length(dir);
            const steps = 160;
            for (let i = 0; i <= steps; ++i) {
                const t = i / steps;
                const X = add(S.p[0], mul(t, dir));
                if (inSolidCone(C, X, -1e-6)) {
                    expect(si).not.toBeNull();
                    expect(t).toBeGreaterThan(lo - 1e-5 / scale);
                    expect(t).toBeLessThan(hi + 1e-5 / scale);
                }
                if (si !== null && t > lo + 1e-6 / scale
                    && t < hi - 1e-6 / scale) {
                    expect(inSolidCone(C, X, 1e-6)).toBe(true);
                }
            }
        }, 60);
    }, 30000);

    it('reports points that lie on the segment and in the cone', () => {
        check(coneAndSegment, ([S, C]) => {
            if (!wellConditioned(S.p[0], sub(S.p[1], S.p[0]), C)) {
                return;
            }
            const sr = fiv.find(S, C);
            if (!sr.intersect) {
                return;
            }
            // The segment query never reports an unbounded type.
            expect(sr.type === T.isPoint || sr.type === T.isSegment).toBe(true);
            const dir = sub(S.p[1], S.p[0]);
            for (let k = 0; k < 2; ++k) {
                const t = intrLine3Cone3Convert(sr.t[k]);
                expect(t).toBeGreaterThan(-1e-12);
                expect(t).toBeLessThan(1 + 1e-12);
                const P = intrLine3Cone3ConvertPoint(sr.P[k]);
                expect(inSolidCone(C, P, 1e-7)).toBe(true);
                expectVectorClose(P, add(S.p[0], mul(t, dir)), 1e-8, 1e-8);
            }
            expect(intrLine3Cone3Convert(sr.t[1]))
                .toBeGreaterThan(intrLine3Cone3Convert(sr.t[0]) - 1e-12);
        });
    });

    it('is unchanged when the segment endpoints are swapped', () => {
        check(coneAndSegment, ([S, C]) => {
            if (!wellConditioned(S.p[0], sub(S.p[1], S.p[0]), C)
                || !wellConditioned(S.p[1], sub(S.p[0], S.p[1]), C)) {
                return;
            }
            const a = fiv.find(S, C);
            const b = fiv.find(Segment.fromEndpoints(S.p[1], S.p[0]), C);
            const ia = intervalOf(a);
            const ib = intervalOf(b);
            if (ia === null || ib === null) {
                // Grazing tangency can round either way; assert only when the
                // reported interval has a real length.
                if (ia !== null && ia[1] - ia[0] > 1e-6) {
                    expect(ib).not.toBeNull();
                }
                if (ib !== null && ib[1] - ib[0] > 1e-6) {
                    expect(ia).not.toBeNull();
                }
                return;
            }
            // Reversing the segment maps the parameter t to 1 - t, so the
            // endpoints of the reported interval swap.
            const dir = sub(S.p[1], S.p[0]);
            expectVectorClose(add(S.p[0], mul(ia[0], dir)),
                add(S.p[1], mul(ib[1], sub(S.p[0], S.p[1]))), 1e-7, 1e-7);
            expectVectorClose(add(S.p[0], mul(ia[1], dir)),
                add(S.p[1], mul(ib[0], sub(S.p[0], S.p[1]))), 1e-7, 1e-7);
        });
    });
});
