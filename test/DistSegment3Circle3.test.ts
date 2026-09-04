import { describe, expect, it } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistLine3Circle3 } from '../src/DistLine3Circle3.js';
import { DistSegment3Circle3 } from '../src/DistSegment3Circle3.js';
import { Line, type Line3 } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
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

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
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

// Brute-force minimum over the segment parameter: a dense sampling followed
// by a local refinement.
function bruteForce(s: Segment, c: Circle3): number {
    const direction = sub(s.p[1], s.p[0]);
    const at = (t: number): number =>
        pointCircleDistance(add(s.p[0], mul(t, direction)), c);
    const n = 8000;
    let best = Number.MAX_VALUE;
    let bt = 0;
    for (let i = 0; i <= n; ++i) {
        const t = i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    let h = 1 / n;
    for (let pass = 0; pass < 120; ++pass) {
        for (const sign of [1, -1]) {
            const t = Math.min(1, Math.max(0, bt + sign * h));
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

describe('DistSegment3Circle3', () => {
    const query = new DistSegment3Circle3();
    const unitCircle = circle([0, 0, 0], [0, 0, 1], 2);

    it('clamps to the first endpoint', () => {
        const s = segment([5, 0, 0], [10, 0, 0]);
        const result = query.compute(s, unitCircle);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.linearClosest[0].values).toEqual([5, 0, 0]);
    });

    it('clamps to the second endpoint', () => {
        const s = segment([10, 0, 0], [5, 0, 0]);
        const result = query.compute(s, unitCircle);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.linearClosest[0].values).toEqual([5, 0, 0]);
    });

    it('reports zero distance when the segment meets the circle', () => {
        const s = segment([-5, 0, 0], [5, 0, 0]);
        const result = query.compute(s, unitCircle);
        expect(result.distance).toBeCloseTo(0, 10);
        for (let j = 0; j < result.numClosestPairs; ++j) {
            verifyPair(unitCircle, result.linearClosest[j],
                result.circularClosest[j], result.distance);
        }
    });

    it('keeps the interior solution when the critical point is on the segment',
        () => {
            const s = segment([0, 0, -5], [0, 0, 5]);
            const result = query.compute(s, unitCircle);
            expect(result.distance).toBeCloseTo(2, 10);
            expect(result.linearClosest[0].values[2]).toBeCloseTo(0, 10);
        });

    it('handles a nearly degenerate segment', () => {
        // A zero-length segment is not a valid input (the algorithm divides
        // by Dot(M,M) with M = P1 - P0), but a very short segment behaves
        // like a point query.
        const s = segment([5, 0, 0], [5, 1e-9, 0]);
        const result = query.compute(s, unitCircle);
        expect(result.distance).toBeCloseTo(3, 8);
    });

    it('handles a degenerate zero-radius circle', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 0);
        const s = segment([3, 0, 0], [6, 0, 0]);
        const result = query.compute(s, c);
        expect(result.distance).toBeCloseTo(3, 8);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 97531864;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        let compared = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const c = circle([2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [2 * rand() - 1, 2 * rand() - 1, 0.2 + rand()],
                0.5 + 2 * rand());
            const s = segment(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4]);
            const result = query.compute(s, c);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);
            for (let j = 0; j < result.numClosestPairs; ++j) {
                verifyPair(c, result.linearClosest[j],
                    result.circularClosest[j], result.distance);
            }
            const line = Line.fromOriginDirection(s.p[0],
                sub(s.p[1], s.p[0]));
            expect(lineSolverIsReliable(line, c)).toBe(true);
            const brute = bruteForce(s, c);
            expect(result.distance).toBeCloseTo(brute, 6);
            ++compared;
        }
        expect(compared).toBeGreaterThan(50);
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistSegment3Circle3.h. Upstream runs DistLine3Circle3's Execute to
// get the critical points of the line-circle squared distance and then clamps
// to 0 <= t <= 1, so the properties check that clamping as well as the
// geometry.
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
    // (a) and its neighbourhood: a3 = |N x E|^2 that is tiny relative to
    // |E|^2 is the rounding-level version of "the line meets the axis", and
    // tauHat, gTauHat and the resulting bracket lose all their significant
    // digits well before a3 reaches exactly zero.
    return length(cross(c.normal, E)) > 1e-6 * length(E);
}

const v21Segment = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -8, 8))
    .filter(([a, b]) => length(sub(b, a)) > 1e-2)
    .map(([a, b]) => Segment.fromEndpoints(a, b));

// The segment parameter of a point known to be on the segment's line, with
// upstream's convention P0 + t*(P1-P0).
function v21SegParameter(s: Segment, x: Vector): number {
    const dir = sub(s.p[1], s.p[0]);
    return dot(sub(x, s.p[0]), dir) / dot(dir, dir);
}

describe('DistSegment3Circle3 verification', () => {
    const query = new DistSegment3Circle3();

    it('every reported pair is on the segment and on the circle', () => {
        check(fc.tuple(v21Segment, v21Circle), ([s, c]) => {
            if (!v21LineSolverApplies(s.p[0], sub(s.p[1], s.p[0]), c)) {
                return;
            }
            const res = query.compute(s, c);
            expect(res.numClosestPairs === 1 || res.numClosestPairs === 2)
                .toBe(true);
            expectClose(res.sqrDistance, res.distance * res.distance, 1e-12,
                1e-12);
            const dir = sub(s.p[1], s.p[0]);
            for (let j = 0; j < res.numClosestPairs; ++j) {
                const linear = res.linearClosest[j];
                const circular = res.circularClosest[j];
                const t = v21SegParameter(s, linear);
                expect(t).toBeGreaterThanOrEqual(-1e-9);
                expect(t).toBeLessThanOrEqual(1 + 1e-9);
                expectVectorClose(linear, add(s.p[0], mul(t, dir)), 1e-8,
                    1e-8);
                const delta = sub(circular, c.center);
                expectClose(dot(c.normal, delta), 0, 1e-8, 1e-8);
                expectClose(length(delta), c.radius, 1e-8, 1e-8);
                expectClose(length(sub(linear, circular)), res.distance, 1e-8,
                    1e-8);
            }
        });
    });

    it('keeps the line result when every critical point is inside [0,1]',
        () => {
            check(fc.tuple(v21Segment, v21Circle), ([s, c]) => {
                if (!v21LineSolverApplies(s.p[0], sub(s.p[1], s.p[0]), c)) {
                    return;
                }
                const line = Line.fromOriginDirection(s.p[0],
                    sub(s.p[1], s.p[0]));
                const { result: lr, critical } =
                    distLine3Circle3Execute(line, c);
                const inside = critical.numPoints === 1
                    ? critical.parameter[0] > 0 && critical.parameter[0] < 1
                    : critical.parameter[0] >= 0 && critical.parameter[1] <= 1
                        && critical.parameter[0] < 1
                        && critical.parameter[1] > 0;
                const sr = query.compute(s, c);
                if (inside) {
                    expect(sr.numClosestPairs).toBe(lr.numClosestPairs);
                    expect(sr.distance).toBe(lr.distance);
                    for (let j = 0; j < lr.numClosestPairs; ++j) {
                        expectVectorClose(sr.linearClosest[j],
                            lr.linearClosest[j], 0, 0);
                        expectVectorClose(sr.circularClosest[j],
                            lr.circularClosest[j], 0, 0);
                    }
                }
            });
        });

    it('matches a brute-force minimization along the segment', () => {
        check(fc.tuple(v21Segment, v21Circle), ([s, c]) => {
            if (!v21LineSolverApplies(s.p[0], sub(s.p[1], s.p[0]), c)) {
                return;
            }
            // Check the shared line solver on the same line, but with a
            // unit direction and a bracket wide enough to contain the closest
            // line point: the segment direction can be very short, so the
            // line parameter of that point can be in the thousands.
            const unit = sub(s.p[1], s.p[0]);
            normalize(unit);
            const unitLine = Line.fromOriginDirection(s.p[0], unit);
            const range = 2 * (length(sub(s.p[0], c.center)) + c.radius + 1);
            expect(lineSolverIsReliable(unitLine, c, range)).toBe(true);
            expectClose(query.compute(s, c).distance, bruteForce(s, c), 1e-5,
                1e-5);
        }, 20);
    }, 30000);

    it('is not larger than the distance to any sampled segment point', () => {
        check(fc.tuple(v21Segment, v21Circle,
            fc.double({ min: 0, max: 1, noNaN: true })), ([s, c, t]) => {
            if (!v21LineSolverApplies(s.p[0], sub(s.p[1], s.p[0]), c)) {
                return;
            }
            const x = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
            expect(query.compute(s, c).distance)
                .toBeLessThanOrEqual(pointCircleDistance(x, c) + 1e-8);
        });
    });

    it('reduces to the point-circle query for a segment pointing away', () => {
        // Both endpoints on the outward radial ray at the same angle: the
        // distance to the circle grows along the ray, so the nearer endpoint
        // is the closest segment point.
        check(fc.tuple(v21Circle, wellScaled(-Math.PI, Math.PI),
            positive(4, 1), positive(4, 1)), ([c, angle, e0, e1]) => {
            const u = getOrthogonal(c.normal, true);
            const w = cross(c.normal, u);
            const radial = add(mul(Math.cos(angle), u),
                mul(Math.sin(angle), w));
            const near = Math.min(e0, e1);
            const far = Math.max(e0, e1);
            if (far - near < 1e-3) {
                return;
            }
            const s = Segment.fromEndpoints(
                add(c.center, mul(c.radius + far, radial)),
                add(c.center, mul(c.radius + near, radial)));
            if (!v21LineSolverApplies(s.p[0], sub(s.p[1], s.p[0]), c)) {
                return;
            }
            const res = query.compute(s, c);
            expect(res.numClosestPairs).toBe(1);
            expectClose(res.distance, near, 1e-7, 1e-7);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(v21Segment, v21Circle, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([s, c, R, tr]) => {
            const rot = (x: Vector): Vector => {
                let y = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    y = add(y, mul(x.values[i], R[i]));
                }
                return y;
            };
            const moved = Circle3.fromCenterNormalRadius(
                add(rot(c.center), tr), rot(c.normal), c.radius);
            const movedSeg = Segment.fromEndpoints(add(rot(s.p[0]), tr),
                add(rot(s.p[1]), tr));
            if (!v21LineSolverApplies(s.p[0], sub(s.p[1], s.p[0]), c)
                || !v21LineSolverApplies(movedSeg.p[0],
                    sub(movedSeg.p[1], movedSeg.p[0]), moved)) {
                return;
            }
            // The bisection in the shared line solver is path dependent, so
            // the distance drifts by more than machine precision under a
            // change of frame.
            expectClose(query.compute(s, c).distance,
                query.compute(movedSeg, moved).distance, 1e-7, 1e-7);
        });
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(v21Segment, v21Circle), ([s, c]) => {
            if (!v21LineSolverApplies(s.p[0], sub(s.p[1], s.p[0]), c)) {
                return;
            }
            const p0 = s.p[0].clone();
            const p1 = s.p[1].clone();
            const snapshot = [...c.center.values, ...c.normal.values,
                c.radius];
            const res = query.compute(s, c);
            expect(s.p[0].values).toEqual(p0.values);
            expect(s.p[1].values).toEqual(p1.values);
            expect([...c.center.values, ...c.normal.values, c.radius])
                .toEqual(snapshot);
            res.linearClosest[0].values[0] = 888;
            res.circularClosest[0].values[0] = 888;
            expect(s.p[0].values).toEqual(p0.values);
            expect([...c.center.values, ...c.normal.values, c.radius])
                .toEqual(snapshot);
        });
    });
    it('documents the shared solver failing when the segment meets the axis',
        () => {
            // See the v21MeetsCircleAxis comment above. The segment below
            // crosses the circle's axis at (0,0,1); upstream's PDFSection422
            // then calls its bisector with an inverted bracket. The test is
            // written so that it also passes once the shared solver is fixed,
            // in which case the answer must be the brute-force minimum.
            const c = circle([0, 0, 0], [0, 0, 1], 0.2);
            const s = segment([0.2, 0, 3], [-0.2, 0, -1]);
            expect(v21LineSolverApplies(s.p[0], sub(s.p[1], s.p[0]), c))
                .toBe(false);
            let threw = false;
            let distance = 0;
            try {
                distance = query.compute(s, c).distance;
            }
            catch {
                threw = true;
            }
            if (!threw) {
                expect(distance).toBeCloseTo(bruteForce(s, c), 5);
            }
        });
});
