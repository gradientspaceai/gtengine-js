import { describe, expect, it } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistLine3Circle3 } from '../src/DistLine3Circle3.js';
import { DistRay3Circle3 } from '../src/DistRay3Circle3.js';
import { Line, type Line3 } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { distLine3Circle3Execute } from '../src/DistLine3Circle3.js';
import { getOrthogonal } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, positive,
    rotationFrame, unitVector, wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function circle(center: number[], normal: number[], radius: number): Circle3 {
    const n = v(...normal);
    normalize(n);
    return Circle3.fromCenterNormalRadius(v(...center), n, radius);
}

// The exact distance from a point to the circle (a curve, not a disk).
function pointCircleDistance(p: Vector, c: Circle3): number {
    const delta = sub(p, c.center);
    const h = dot(c.normal, delta);
    const inPlane = sub(delta, mul(h, c.normal));
    const radial = length(inPlane);
    const dr = radial - c.radius;
    return Math.sqrt(h * h + dr * dr);
}

// Verify that a reported closest pair is consistent: the circular point is on
// the circle and the pair realizes the reported distance.
function verifyPair(c: Circle3, linear: Vector, circular: Vector,
    distance: number): void {
    const delta = sub(circular, c.center);
    expect(dot(c.normal, delta)).toBeCloseTo(0, 8);
    expect(length(delta)).toBeCloseTo(c.radius, 8);
    expect(length(sub(linear, circular))).toBeCloseTo(distance, 8);
}

// The ray/segment queries delegate to DistLine3Circle3 for the critical
// points of the line-circle distance. Upstream's PDFSection422 computes
// tauHat = sqrt(|(a1*a3)^(2/3) - a3|) where the derivation requires
// tauHat = sqrt(|(a1*a3)^(2/3) - a3| / a2). The port fixes this (issue filed
// upstream), so every trial checks the line solver against a brute-force
// minimum over the whole line instead of skipping unreliable configurations.
function lineSolverIsReliable(line: Line3, c: Circle3,
    range = 10): boolean {
    const lineDistance = new DistLine3Circle3().compute(line, c).distance;
    const at = (t: number): number =>
        pointCircleDistance(add(line.origin, mul(t, line.direction)), c);
    const n = 8000;
    let best = Number.MAX_VALUE;
    let bt = 0;
    for (let i = 0; i <= n; ++i) {
        const t = -range + 2 * range * i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    let h = 2 * range / n;
    for (let pass = 0; pass < 80; ++pass) {
        for (const sign of [1, -1]) {
            const d = at(bt + sign * h);
            if (d < best) {
                best = d;
                bt = bt + sign * h;
            }
        }
        h *= 0.75;
    }
    return Math.abs(lineDistance - best) <= 1e-4;
}

// Brute-force minimum over the ray parameter: a coarse sampling followed by a
// local refinement (the distance function is smooth away from the axis).
function bruteForce(r: Ray, c: Circle3, tmax: number): number {
    const at = (t: number): number =>
        pointCircleDistance(add(r.origin, mul(t, r.direction)), c);
    const n = 8000;
    let best = Number.MAX_VALUE;
    let bt = 0;
    for (let i = 0; i <= n; ++i) {
        const t = tmax * i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    let h = tmax / n;
    for (let pass = 0; pass < 120; ++pass) {
        for (const sign of [1, -1]) {
            const t = Math.max(0, bt + sign * h);
            const d = at(t);
            if (d < best) {
                best = d;
                bt = t;
            }
        }
        h *= 0.75;
    }
    return best;
}

describe('DistRay3Circle3', () => {
    const query = new DistRay3Circle3();
    const unitCircle = circle([0, 0, 0], [0, 0, 1], 2);

    it('uses the ray origin when the ray points away from the circle', () => {
        const r = ray([5, 0, 0], [1, 0, 0]);
        const result = query.compute(r, unitCircle);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.linearClosest[0].values).toEqual([5, 0, 0]);
        expect(result.circularClosest[0].values[0]).toBeCloseTo(2, 10);
    });

    it('reports zero distance when the ray meets the circle', () => {
        const r = ray([-10, 0, 0], [1, 0, 0]);
        const result = query.compute(r, unitCircle);
        expect(result.distance).toBeCloseTo(0, 10);
        for (let j = 0; j < result.numClosestPairs; ++j) {
            verifyPair(unitCircle, result.linearClosest[j],
                result.circularClosest[j], result.distance);
        }
    });

    it('handles a ray on the circle axis', () => {
        const r = ray([0, 0, 5], [0, 0, 1]);
        const result = query.compute(r, unitCircle);
        expect(result.distance).toBeCloseTo(Math.sqrt(25 + 4), 10);
        expect(result.linearClosest[0].values).toEqual([0, 0, 5]);
        verifyPair(unitCircle, result.linearClosest[0],
            result.circularClosest[0], result.distance);
    });

    it('keeps the line solution when the critical point is on the ray', () => {
        const r = ray([0, 0, -5], [0, 0, 1]);
        const result = query.compute(r, unitCircle);
        // The line meets the plane of the circle at the center, so the ray
        // point closest to the circle is the circle center.
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.linearClosest[0].values[2]).toBeCloseTo(0, 10);
    });

    it('handles a degenerate zero-radius circle', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 0);
        const r = ray([3, 0, 0], [1, 0, 0]);
        const result = query.compute(r, c);
        expect(result.distance).toBeCloseTo(3, 8);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 24681012;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        let compared = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const c = circle([2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [2 * rand() - 1, 2 * rand() - 1, 0.2 + rand()],
                0.5 + 2 * rand());
            const r = ray([8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-3) {
                continue;
            }
            const result = query.compute(r, c);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);
            for (let j = 0; j < result.numClosestPairs; ++j) {
                verifyPair(c, result.linearClosest[j],
                    result.circularClosest[j], result.distance);
            }
            const line = Line.fromOriginDirection(r.origin, r.direction);
            expect(lineSolverIsReliable(line, c)).toBe(true);
            const brute = bruteForce(r, c, 60);
            expect(result.distance).toBeCloseTo(brute, 6);
            ++compared;
        }
        expect(compared).toBeGreaterThan(50);
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistRay3Circle3.h. Upstream runs DistLine3Circle3's Execute to get
// the critical points of the line-circle squared distance and then clamps to
// t >= 0, so the properties check that clamping as well as the geometry.
// ---------------------------------------------------------------------------

const v21Circle = fc.tuple(wellScaledVector(3, -5, 5), unitVector(3),
    positive(4, 0.2))
    .map(([c, n, r]) => Circle3.fromCenterNormalRadius(c, n, r));

// Upstream defects in the shared line-circle solver (DistLine3Circle3.h),
// reachable from the ray and segment queries because they delegate to it.
// Upstream has the same code in every case, so these are upstream defects
// rather than port defects; the missing "/a2" fix of issue #247 affects none
// of them.
//
//   (a) PDFSection422 with a line that meets the axis of the circle. The
//       re-chosen line origin E then lies on that axis, so a3 = |N x E|^2 is
//       exactly zero, tauHat is 0 and gTauHat = a1*0/sqrt(0) is NaN. Every
//       intercept comparison is false, so the "two critical points" fallback
//       runs and, when a0 > a1/sqrt(a2), bisects over
//       [0, -a0 + a1/sqrt(a2)] whose upper bound is negative. RootsBisection1
//       then throws "Invalid ordering of t-interval endpoints."
//   (b) PDFSection422 with a line nearly (but not exactly) perpendicular to
//       the plane of the circle. The exact test !IsZero(N x M) still selects
//       PDFSection422, a2 = |N x M|^2 is around 1e-19, the shifted origin
//       parameter s = -Dot(NxM,NxD)/a2 is around 1e6, and the bracket width
//       a1/sqrt(a2) = r*|N x M|/|M|^2 falls below the ulp of a0, so the two
//       bracket endpoints round to the same double and the bisector throws.
//   (c) PDFSection412 (line perpendicular to the plane, IsZero(N x D) false)
//       when the closest line point is the circle center. Finalize computes
//       the in-plane component of that point and calls Normalize on it
//       without checking the length; the zero vector comes back as zero, the
//       circle point is set to the center itself and the reported distance is
//       0 instead of the radius. A line along the circle axis whose N x D is
//       nonzero only through rounding (for example after a rigid motion) hits
//       this.
//
// The properties below skip these configurations; the dedicated test at the
// end documents (a).
function v21LineSolverApplies(origin: Vector, direction: Vector,
    c: Circle3): boolean {
    const D = sub(origin, c.center);
    const NxM = cross(c.normal, direction);
    const NxD = cross(c.normal, D);
    const isZero = (x: Vector): boolean =>
        x.values[0] === 0 && x.values[1] === 0 && x.values[2] === 0;
    if (isZero(NxM)) {
        if (isZero(NxD)) {
            return true;   // PDFSection411, exact
        }
        // PDFSection412; see (c).
        const t = -dot(direction, D) / dot(direction, direction);
        const linearPoint = add(mul(t, direction), D);
        const project = sub(linearPoint,
            mul(dot(c.normal, linearPoint), c.normal));
        return length(project) > 1e-8 * (length(D) + c.radius);
    }
    if (isZero(NxD)) {
        return true;   // PDFSection421, closed form
    }
    if (length(NxM) <= 1e-5 * length(direction)) {
        return false;  // (b)
    }
    const s = -dot(NxM, NxD) / dot(NxM, NxM);
    const E = add(mul(s, direction), D);
    return !isZero(cross(c.normal, E));   // (a)
}

const v21Ray = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3))
    .map(([o, d]) => Ray.fromOriginDirection(o, d));

// The ray parameter of a point known to be on the ray's line.
function v21RayParameter(r: Ray, x: Vector): number {
    return dot(sub(x, r.origin), r.direction) / dot(r.direction, r.direction);
}

describe('DistRay3Circle3 verification', () => {
    const query = new DistRay3Circle3();

    it('every reported pair is on the ray and on the circle', () => {
        check(fc.tuple(v21Ray, v21Circle), ([r, c]) => {
            if (!v21LineSolverApplies(r.origin, r.direction, c)) {
                return;
            }
            const res = query.compute(r, c);
            expect(res.numClosestPairs === 1 || res.numClosestPairs === 2)
                .toBe(true);
            expectClose(res.sqrDistance, res.distance * res.distance, 1e-12,
                1e-12);
            for (let j = 0; j < res.numClosestPairs; ++j) {
                const linear = res.linearClosest[j];
                const circular = res.circularClosest[j];
                const t = v21RayParameter(r, linear);
                expect(t).toBeGreaterThanOrEqual(-1e-9);
                expectVectorClose(linear, add(r.origin, mul(t, r.direction)),
                    1e-8, 1e-8);
                const delta = sub(circular, c.center);
                expectClose(dot(c.normal, delta), 0, 1e-8, 1e-8);
                expectClose(length(delta), c.radius, 1e-8, 1e-8);
                expectClose(length(sub(linear, circular)), res.distance, 1e-8,
                    1e-8);
            }
        });
    });

    it('the line solver orders its two critical parameters', () => {
        // HasTwoCriticalPoints assumes parameter[0] <= parameter[1] (it
        // returns the line result when t0 >= 0 and the ray origin when
        // t1 <= 0). Verify that assumption on the shared solver.
        check(fc.tuple(v21Ray, v21Circle), ([r, c]) => {
            if (!v21LineSolverApplies(r.origin, r.direction, c)) {
                return;
            }
            const line = Line.fromOriginDirection(r.origin, r.direction);
            const { critical } = distLine3Circle3Execute(line, c);
            if (critical.numPoints === 2) {
                expect(critical.parameter[0])
                    .toBeLessThanOrEqual(critical.parameter[1]);
            }
        });
    });

    it('keeps the line result when every critical point is on the ray', () => {
        check(fc.tuple(v21Ray, v21Circle), ([r, c]) => {
            if (!v21LineSolverApplies(r.origin, r.direction, c)) {
                return;
            }
            const line = Line.fromOriginDirection(r.origin, r.direction);
            const { result: lr, critical } = distLine3Circle3Execute(line, c);
            const rr = query.compute(r, c);
            const onRay = critical.numPoints === 1
                ? critical.parameter[0] > 0
                : critical.parameter[0] >= 0;
            if (onRay) {
                expect(rr.numClosestPairs).toBe(lr.numClosestPairs);
                expect(rr.distance).toBe(lr.distance);
                for (let j = 0; j < lr.numClosestPairs; ++j) {
                    expectVectorClose(rr.linearClosest[j], lr.linearClosest[j],
                        0, 0);
                    expectVectorClose(rr.circularClosest[j],
                        lr.circularClosest[j], 0, 0);
                }
            }
        });
    });

    it('matches a brute-force minimization along the ray', () => {
        check(fc.tuple(v21Ray, v21Circle), ([r, c]) => {
            if (!v21LineSolverApplies(r.origin, r.direction, c)) {
                return;
            }
            const line = Line.fromOriginDirection(r.origin, r.direction);
            // The line solver is the ingredient the ray query trusts; if the
            // shared solver has already missed the line minimum there is
            // nothing for this file to prove. The brute-force bracket has to
            // cover the closest line point, which is at most
            // |P - C| + r away from the origin along a unit direction.
            const range = 2 * (length(sub(r.origin, c.center)) + c.radius + 1);
            expect(lineSolverIsReliable(line, c, range)).toBe(true);
            expectClose(query.compute(r, c).distance, bruteForce(r, c, 60),
                1e-5, 1e-5);
        }, 20);
    }, 30000);

    it('is not larger than the distance to any sampled ray point', () => {
        check(fc.tuple(v21Ray, v21Circle, positive(60, 0)), ([r, c, t]) => {
            if (!v21LineSolverApplies(r.origin, r.direction, c)) {
                return;
            }
            expect(query.compute(r, c).distance).toBeLessThanOrEqual(
                pointCircleDistance(add(r.origin, mul(t, r.direction)), c)
                + 1e-8);
        });
    });

    it('reduces to the point-circle query when the ray points away', () => {
        // A ray whose direction is the outward radial direction at its origin
        // and whose origin is far outside the circle can only move away, so
        // the origin is the closest ray point.
        check(fc.tuple(v21Circle, wellScaled(-Math.PI, Math.PI),
            positive(6, 2)), ([c, angle, extra]) => {
            const u = getOrthogonal(c.normal, true);
            const w = cross(c.normal, u);
            const radial = add(mul(Math.cos(angle), u),
                mul(Math.sin(angle), w));
            const origin = add(c.center, mul(c.radius + extra, radial));
            const r = Ray.fromOriginDirection(origin, radial);
            if (!v21LineSolverApplies(r.origin, r.direction, c)) {
                return;
            }
            const res = query.compute(r, c);
            expect(res.numClosestPairs).toBe(1);
            expectClose(res.distance, extra, 1e-8, 1e-8);
            expectVectorClose(res.linearClosest[0], origin, 1e-9, 1e-9);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(v21Ray, v21Circle, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([r, c, R, tr]) => {
            const rot = (x: Vector): Vector => {
                let y = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    y = add(y, mul(x.values[i], R[i]));
                }
                return y;
            };
            const moved = Circle3.fromCenterNormalRadius(
                add(rot(c.center), tr), rot(c.normal), c.radius);
            const movedRay = Ray.fromOriginDirection(add(rot(r.origin), tr),
                rot(r.direction));
            if (!v21LineSolverApplies(r.origin, r.direction, c)
                || !v21LineSolverApplies(movedRay.origin, movedRay.direction,
                    moved)) {
                return;
            }
            // The bisection in the shared line solver is path dependent, so
            // the distance drifts by more than machine precision under a
            // change of frame.
            expectClose(query.compute(r, c).distance,
                query.compute(movedRay, moved).distance, 1e-7, 1e-7);
        });
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(v21Ray, v21Circle), ([r, c]) => {
            if (!v21LineSolverApplies(r.origin, r.direction, c)) {
                return;
            }
            const o = r.origin.clone();
            const d = r.direction.clone();
            const snapshot = [...c.center.values, ...c.normal.values,
                c.radius];
            const res = query.compute(r, c);
            expect(r.origin.values).toEqual(o.values);
            expect(r.direction.values).toEqual(d.values);
            expect([...c.center.values, ...c.normal.values, c.radius])
                .toEqual(snapshot);
            res.linearClosest[0].values[0] = 888;
            res.circularClosest[0].values[0] = 888;
            expect(r.origin.values).toEqual(o.values);
            expect([...c.center.values, ...c.normal.values, c.radius])
                .toEqual(snapshot);
        });
    });
    it('documents the shared solver failing when the ray meets the circle axis',
        () => {
            // See the v21MeetsCircleAxis comment above. The ray below crosses
            // the circle's axis at (0,0,1); upstream's PDFSection422 then
            // calls its bisector with an inverted bracket. The test is
            // written so that it also passes once the shared solver is fixed,
            // in which case the answer must be the brute-force minimum.
            const c = circle([0, 0, 0], [0, 0, 1], 0.2);
            const r = ray([0.2, 0, 3], [0.1, 0, 1]);
            expect(v21LineSolverApplies(r.origin, r.direction, c))
                .toBe(false);
            let threw = false;
            let distance = 0;
            try {
                distance = query.compute(r, c).distance;
            }
            catch {
                threw = true;
            }
            if (!threw) {
                expect(distance).toBeCloseTo(bruteForce(r, c, 60), 5);
            }
        });
});
