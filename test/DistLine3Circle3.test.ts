import { describe, expect, it } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistLine3Circle3, distLine3Circle3Execute } from '../src/DistLine3Circle3.js';
import { Line } from '../src/Line.js';
import {
    Vector, add, dot, length, mul, normalize, sub
} from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    seededRandom, unitVector, wellScaledVector
} from './helpers/arbitraries.js';
import { cross } from '../src/Vector3.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function circle(center: number[], normal: number[], radius: number): Circle3 {
    const n = v(...normal);
    normalize(n);
    return Circle3.fromCenterNormalRadius(v(...center), n, radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistLine3Circle3', () => {
    const query = new DistLine3Circle3();
    const unit = circle([0, 0, 0], [0, 0, 1], 1);

    it('reports equidistance for the normal line through the center', () => {
        const result = query.compute(line([0, 0, -5], [0, 0, 1]), unit);
        expect(result.equidistant).toBe(true);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(1, 12);
        expect(result.linearClosest[0].values).toEqual([0, 0, 0]);
        const d = sub(result.circularClosest[0], unit.center);
        expect(Math.sqrt(dot(d, d))).toBeCloseTo(1, 10);
    });

    it('handles a line perpendicular to the plane but off the center', () => {
        // The vertical line through (3,0,z); the closest circle point is
        // (1,0,0) and the closest line point is (3,0,0).
        const result = query.compute(line([3, 0, -7], [0, 0, 1]), unit);
        expect(result.equidistant).toBe(false);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.linearClosest[0].values[2]).toBeCloseTo(0, 10);
        expect(result.circularClosest[0].values[0]).toBeCloseTo(1, 10);
    });

    it('reports two closest pairs for a line on the circle axis line',
        () => {
            // A line through the circle center, in the plane of the circle:
            // it meets the circle at two diametrically opposite points.
            const result = query.compute(line([0, 0, 0], [1, 0, 0]), unit);
            expect(result.numClosestPairs).toBe(2);
            expect(result.distance).toBeCloseTo(0, 10);
            expect(result.circularClosest[0].values[0]).toBeCloseTo(-1, 10);
            expect(result.circularClosest[1].values[0]).toBeCloseTo(1, 10);
        });

    it('reports two closest pairs for a symmetric skew line', () => {
        // A line parallel to the x-axis at height z = 2 through the axis of
        // the circle. It is symmetric about the yz-plane, so both (+-1,0,0)
        // are closest.
        const result = query.compute(line([0, 0, 2], [1, 0, 0]), unit);
        expect(result.numClosestPairs).toBe(2);
        expect(result.distance).toBeCloseTo(2, 10);
    });

    it('measures a line in the plane of the circle that misses it', () => {
        const result = query.compute(line([0, 3, 0], [1, 0, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(2, 8);
        expect(result.circularClosest[0].values[1]).toBeCloseTo(1, 8);
    });

    it('reports zero distance for a line through a circle point', () => {
        const result = query.compute(line([1, 0, 0], [0, 0, 1]), unit);
        expect(result.distance).toBeCloseTo(0, 8);
    });

    it('exposes the critical points through the friend-style export', () => {
        const { result, critical } =
            distLine3Circle3Execute(line([3, 0, -7], [0, 0, 1]), unit);
        expect(critical.numPoints).toBe(1);
        expect(result.numClosestPairs).toBe(1);
        expect(critical.distance[0]).toBeCloseTo(result.distance, 12);
        // The critical line point matches the line parameter.
        const ln = line([3, 0, -7], [0, 0, 1]);
        const p = add(ln.origin, mul(critical.parameter[0], ln.direction));
        for (let i = 0; i < 3; ++i) {
            expect(p.values[i]).toBeCloseTo(critical.linearPoint[0].values[i],
                10);
        }
    });

    it('does not require a unit-length line direction', () => {
        const r0 = query.compute(line([3, 0, -7], [0, 0, 1]), unit);
        const r1 = query.compute(line([3, 0, -7], [0, 0, 5]), unit);
        expect(r1.distance).toBeCloseTo(r0.distance, 10);
        for (let i = 0; i < 3; ++i) {
            expect(r1.linearClosest[0].values[i]).toBeCloseTo(
                r0.linearClosest[0].values[i], 8);
        }
    });

    it('agrees with a dense sampling of the circle', () => {
        const rnd = makeRandom(160901);
        const center = v(0.5, -1, 0.25);
        const normal = v(1, 2, -0.5);
        normalize(normal);
        const radius = 1.5;
        const c = Circle3.fromCenterNormalRadius(center, normal, radius);

        const U = v(-normal.values[1], normal.values[0], 0);
        normalize(U);
        const W = cross(normal, U);

        for (let trial = 0; trial < 40; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, c);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);

            for (let j = 0; j < result.numClosestPairs; ++j) {
                // The circle point is on the circle.
                const d = sub(result.circularClosest[j], center);
                expect(Math.sqrt(dot(d, d))).toBeCloseTo(radius, 6);
                expect(Math.abs(dot(d, normal))).toBeLessThan(1e-6);

                // The line point is on the line.
                const w = sub(result.linearClosest[j], ln.origin);
                const f = sub(w, mul(dot(w, ln.direction)
                    / dot(ln.direction, ln.direction), ln.direction));
                expect(Math.sqrt(dot(f, f))).toBeLessThan(1e-6);

                // The pair realizes the reported distance.
                const e = sub(result.linearClosest[j],
                    result.circularClosest[j]);
                expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 5);
            }

            // No sampled circle point is closer to the line.
            const n = 4000;
            const dd = dot(ln.direction, ln.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i < n; ++i) {
                const t = 2 * Math.PI * i / n;
                const q = add(center, add(mul(radius * Math.cos(t), U),
                    mul(radius * Math.sin(t), W)));
                const w = sub(q, ln.origin);
                const s = dot(w, ln.direction) / dd;
                const f = sub(w, mul(s, ln.direction));
                best = Math.min(best, dot(f, f));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-5);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistLine3Circle3.h, including the tauHat fix of
// issue #247 (upstream omits the division by a2 when solving G'(tau) = 1).
// ---------------------------------------------------------------------------

describe('DistLine3Circle3 verification', () => {
    const query = new DistLine3Circle3();

    const circleArb = fc.tuple(wellScaledVector(3, -5, 5), unitVector(3),
        finite(0.25, 4)).map(([c, n, r]) =>
            Circle3.fromCenterNormalRadius(c, n, r));

    const lineArb = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

    // Distance from a point to the line P + t*D (D is not required to be unit
    // length).
    function pointLineDistance(p: Vector, ln: Line): number {
        const diff = sub(p, ln.origin);
        const t = dot(diff, ln.direction) / dot(ln.direction, ln.direction);
        return length(sub(diff, mul(t, ln.direction)));
    }

    // An orthonormal basis {u, w} for the plane of the circle.
    function planeBasis(circle: Circle3): [Vector, Vector] {
        const n = circle.normal;
        const seed = Math.abs(n.values[0]) < 0.9 ? v(1, 0, 0) : v(0, 1, 0);
        const u = sub(seed, mul(dot(seed, n), n));
        normalize(u);
        return [u, cross(n, u)];
    }

    // Independent minimization of the point-line distance over the circle: a
    // dense scan followed by a golden-section refinement of the bracketing
    // interval. The scan makes no assumption about the number of local
    // minima, so it detects a query that converged to the wrong critical
    // point (which is exactly the failure mode of the upstream tauHat).
    function bruteForceDistance(ln: Line, circle: Circle3): number {
        const [u, w] = planeBasis(circle);
        const at = (a: number): Vector => add(circle.center,
            add(mul(circle.radius * Math.cos(a), u),
                mul(circle.radius * Math.sin(a), w)));
        const f = (a: number): number => pointLineDistance(at(a), ln);
        const n = 4096;
        let bestI = 0;
        let best = Number.POSITIVE_INFINITY;
        for (let i = 0; i < n; ++i) {
            const y = f((2 * Math.PI * i) / n);
            if (y < best) { best = y; bestI = i; }
        }
        const h = (2 * Math.PI) / n;
        let lo = ((2 * Math.PI * bestI) / n) - h;
        let hi = lo + 2 * h;
        const phi = (Math.sqrt(5) - 1) / 2;
        for (let i = 0; i < 200; ++i) {
            const m0 = hi - phi * (hi - lo);
            const m1 = lo + phi * (hi - lo);
            if (f(m0) <= f(m1)) { hi = m1; } else { lo = m0; }
        }
        return Math.min(best, f(0.5 * (lo + hi)));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(lineArb, circleArb), ([ln, circle]) => {
            const r = query.compute(ln, circle);
            expect(r.numClosestPairs === 1 || r.numClosestPairs === 2)
                .toBe(true);
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            for (let j = 0; j < r.numClosestPairs; ++j) {
                // The linear closest point is on the line.
                expectClose(pointLineDistance(r.linearClosest[j], ln), 0,
                    1e-7, 1e-7);
                // The circular closest point is on the circle: it satisfies
                // |X-C| = r and lies in the plane of the circle.
                const delta = sub(r.circularClosest[j], circle.center);
                expectClose(length(delta), circle.radius, 1e-7, 1e-7);
                expectClose(dot(circle.normal, delta), 0, 1e-7, 1e-7);
                // The reported distance is the distance of this pair.
                expectClose(length(sub(r.linearClosest[j],
                    r.circularClosest[j])), r.distance, 1e-7, 1e-7);
            }
        }, 150);
    });

    it('matches an independent minimization over the circle', () => {
        check(fc.tuple(lineArb, circleArb), ([ln, circle]) => {
            const r = query.compute(ln, circle);
            const best = bruteForceDistance(ln, circle);
            // The scan-plus-refine reference is accurate well below 1e-9; the
            // tolerance covers the bisection tolerance of the query.
            expectClose(r.distance, best, 1e-7, 1e-7);
        }, 60);
    }, 30000);

    it('matches the brute force for lines that nearly touch the circle',
        () => {
            // Near-tangent configurations put a0 close to +-intercept in
            // PDFSection422, which is where tauHat is used both as a bracket
            // endpoint and as a critical point.
            check(fc.tuple(unitVector(3), finite(0.5, 3), unitVector(3),
                finite(-2, 2), finite(-0.2, 0.2)),
                ([n, radius, dir, along, off]) => {
                    const circle = Circle3.fromCenterNormalRadius(
                        v(0, 0, 0), n, radius);
                    const seed = Math.abs(n.values[0]) < 0.9
                        ? v(1, 0, 0) : v(0, 1, 0);
                    const u = sub(seed, mul(dot(seed, n), n));
                    normalize(u);
                    const origin = add(mul(radius + off, u), mul(along, n));
                    const ln = Line.fromOriginDirection(origin, dir);
                    const r = query.compute(ln, circle);
                    expectClose(r.distance, bruteForceDistance(ln, circle),
                        1e-6, 1e-6);
                }, 60);
        }, 30000);

    it('is equidistant only for the line through the center along N', () => {
        check(circleArb, circle => {
            const ln = Line.fromOriginDirection(circle.center, circle.normal);
            const r = query.compute(ln, circle);
            expect(r.equidistant).toBe(true);
            expect(r.numClosestPairs).toBe(1);
            expectVectorClose(r.linearClosest[0], circle.center, 1e-12, 1e-12);
            expectClose(r.distance, circle.radius, 1e-12, 1e-12);
        });
    });

    it('handles a line perpendicular to the plane off the center', () => {
        check(fc.tuple(circleArb, finite(0.1, 6), finite(-3, 3)),
            ([circle, offset, along]) => {
                const [u] = planeBasis(circle);
                const origin = add(circle.center,
                    add(mul(offset, u), mul(along, circle.normal)));
                const ln = Line.fromOriginDirection(origin, circle.normal);
                const r = query.compute(ln, circle);
                expect(r.equidistant).toBe(false);
                expectClose(r.distance, Math.abs(offset - circle.radius),
                    1e-9, 1e-9);
            });
    });

    it('handles a line through the center that is not perpendicular', () => {
        // PDFSection421: the line origin is on the normal line through the
        // center but the line is not parallel to N.
        check(fc.tuple(circleArb, finite(-3, 3), unitVector(3)),
            ([circle, along, dir]) => {
                const origin = add(circle.center, mul(along, circle.normal));
                const ln = Line.fromOriginDirection(origin, dir);
                const r = query.compute(ln, circle);
                expectClose(r.distance, bruteForceDistance(ln, circle),
                    1e-6, 1e-6);
            }, 60);
    }, 30000);

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineArb, circleArb, rotationFrame(3),
            wellScaledVector(3, -6, 6)), ([ln, circle, frame, shift]) => {
            const rot = (p: Vector): Vector =>
                add(add(mul(p.values[0], frame[0]),
                    mul(p.values[1], frame[1])), mul(p.values[2], frame[2]));
            const movedLine = Line.fromOriginDirection(
                add(shift, rot(ln.origin)), rot(ln.direction));
            const movedCircle = Circle3.fromCenterNormalRadius(
                add(shift, rot(circle.center)), rot(circle.normal),
                circle.radius);
            const r0 = query.compute(ln, circle);
            const r1 = query.compute(movedLine, movedCircle);
            // The bisection is path dependent, so the tolerance covers drift
            // over the iterations.
            expectClose(r0.distance, r1.distance, 1e-7, 1e-7);
        }, 100);
    });

    it('is invariant to the length of the line direction', () => {
        check(fc.tuple(lineArb, circleArb, finite(0.1, 10)),
            ([ln, circle, scale]) => {
                const scaled = Line.fromOriginDirection(ln.origin,
                    mul(scale, ln.direction));
                const r0 = query.compute(ln, circle);
                const r1 = query.compute(scaled, circle);
                expectClose(r0.distance, r1.distance, 1e-7, 1e-7);
            }, 100);
    });

    it('exports the critical points used by ray and segment queries', () => {
        check(fc.tuple(lineArb, circleArb), ([ln, circle]) => {
            const { result, critical } = distLine3Circle3Execute(ln, circle);
            expect(critical.numPoints === 1 || critical.numPoints === 2)
                .toBe(true);
            let best = Number.POSITIVE_INFINITY;
            for (let i = 0; i < critical.numPoints; ++i) {
                // Each critical linear point is the line point at its
                // parameter (in the equidistant case the stored point is the
                // circle center, which is the line point there too).
                expectVectorClose(critical.linearPoint[i],
                    add(ln.origin, mul(critical.parameter[i], ln.direction)),
                    1e-6, 1e-6);
                best = Math.min(best, critical.distance[i]);
            }
            expectClose(result.distance, best, 1e-9, 1e-9);
        }, 100);
    });
});

// ---------------------------------------------------------------------------
// Regression for the a3 = |NxE|^2 = 0 defect described in
// src/DistLine3Circle3.ts: upstream builds a NaN 'intercept' and then calls
// Bisect with tauMin > tauMax, which throws "Invalid ordering of t-interval
// endpoints".
// ---------------------------------------------------------------------------

describe('DistLine3Circle3 lines meeting the circle axis', () => {
    const query = new DistLine3Circle3();
    const unit = Circle3.fromCenterNormalRadius(v(0, 0, 0), v(0, 0, 1), 1);

    function pointLineDistance(p: Vector, ln: Line): number {
        const diff = sub(p, ln.origin);
        const t = dot(diff, ln.direction) / dot(ln.direction, ln.direction);
        return length(sub(diff, mul(t, ln.direction)));
    }

    function bruteForce(ln: Line, circle: Circle3): number {
        // The circle lies in the z = 0 plane in these fixtures.
        const n = 200000;
        let best = Number.POSITIVE_INFINITY;
        for (let i = 0; i < n; ++i) {
            const a = (2 * Math.PI * i) / n;
            const p = add(circle.center,
                v(circle.radius * Math.cos(a), circle.radius * Math.sin(a), 0));
            best = Math.min(best, pointLineDistance(p, ln));
        }
        return best;
    }

    it('does not throw for a line crossing the axis away from the center',
        () => {
            // The line meets the z-axis at (0,0,5); upstream throws here.
            const ln = Line.fromOriginDirection(v(1, 0, 6), v(1, 0, 1));
            const r = query.compute(ln, unit);
            expect(Number.isFinite(r.distance)).toBe(true);
            // sqrt(2)/2 * distance in the x-z plane from (0,0,5) to (1,0,0):
            // compare against a dense sampling of the circle.
            expect(Math.abs(r.distance - bruteForce(ln, unit)))
                .toBeLessThan(1e-6);
        }, 30000);

    it('agrees with brute force for a family of axis-crossing lines', () => {
        // height is where the line meets the z-axis, and the direction is
        // oblique so that a0 != 0 (the sub-branch upstream mishandles).
        for (const height of [-4, -1.5, -0.25, 0.25, 1.5, 4]) {
            for (const dir of [[1, 0, 1], [0, 1, 2], [2, 1, -3], [1, 1, 1]]) {
                const axisPoint = v(0, 0, height);
                const d = v(dir[0], dir[1], dir[2]);
                // Offset the origin along the line so that it is not the
                // axis point itself (otherwise the query takes PDFSection421).
                const origin = add(axisPoint, mul(1.7, d));
                const ln = Line.fromOriginDirection(origin, d);
                const r = query.compute(ln, unit);
                expect(Math.abs(r.distance - bruteForce(ln, unit)))
                    .toBeLessThan(1e-5);
                for (let j = 0; j < r.numClosestPairs; ++j) {
                    expect(Math.abs(length(sub(r.circularClosest[j],
                        unit.center)) - unit.radius)).toBeLessThan(1e-9);
                    expect(Math.abs(length(sub(r.linearClosest[j],
                        r.circularClosest[j])) - r.distance))
                        .toBeLessThan(1e-9);
                }
            }
        }
    }, 30000);

    it('matches PDFSection421 when the line origin is moved onto the axis',
        () => {
            // The same line described by two different origins must give the
            // same distance: one origin on the axis (PDFSection421) and one
            // off it (PDFSection422 with a3 = 0).
            for (const height of [-3, -0.5, 0.5, 3]) {
                for (const dir of [[1, 0, 1], [1, 2, 3], [-2, 1, 1]]) {
                    const d = v(dir[0], dir[1], dir[2]);
                    const onAxis = Line.fromOriginDirection(v(0, 0, height), d);
                    const offAxis = Line.fromOriginDirection(
                        add(v(0, 0, height), mul(2.3, d)), d);
                    const r0 = query.compute(onAxis, unit);
                    const r1 = query.compute(offAxis, unit);
                    expect(Math.abs(r0.distance - r1.distance))
                        .toBeLessThan(1e-9);
                }
            }
        });
});

// ---------------------------------------------------------------------------
// Regression for upstream issue #247: PDFSection422 omits the division by a2
// when solving G'(tau) = 1, so tauHat is off by a factor of sqrt(a2). tauHat
// is a bisection bracket endpoint in the three-critical-point branch
// (|a0| < intercept), and when a2 > 1 the bracket no longer contains the
// intended root. a2 = |NxM|^2 exceeds 1 only for a long direction vector,
// which is exactly how the segment-circle query calls this file (M = P1-P0).
// ---------------------------------------------------------------------------

describe('DistLine3Circle3 tauHat bracket (upstream issue #247)', () => {
    const query = new DistLine3Circle3();

    function pointLineDistance(p: Vector, ln: Line): number {
        const diff = sub(p, ln.origin);
        const t = dot(diff, ln.direction) / dot(ln.direction, ln.direction);
        return length(sub(diff, mul(t, ln.direction)));
    }

    function brute(ln: Line, c: Circle3): number {
        const n = c.normal;
        const seed = Math.abs(n.values[0]) < 0.9 ? v(1, 0, 0) : v(0, 1, 0);
        const u = sub(seed, mul(dot(seed, n), n));
        normalize(u);
        const w = cross(n, u);
        let best = Number.POSITIVE_INFINITY;
        const samples = 200000;
        for (let i = 0; i < samples; ++i) {
            const a = (2 * Math.PI * i) / samples;
            const p = add(c.center, add(mul(c.radius * Math.cos(a), u),
                mul(c.radius * Math.sin(a), w)));
            best = Math.min(best, pointLineDistance(p, ln));
        }
        return best;
    }

    it('computes the near-zero distance upstream reports as 2.92', () => {
        // Found by a randomized search over the three-critical-point branch.
        // The upstream tauHat gives 2.9236 for this configuration; the true
        // distance is about 0.0131.
        const circle = Circle3.fromCenterNormalRadius(v(0, 0, 0),
            v(0.3616963707341183, -0.8460173734208769, 0.3917018499673034),
            3.044704365124926);
        const ln = Line.fromOriginDirection(
            v(3.2986813923344016, 0.9599527320824564, -0.40380136808380485),
            v(-5.10283449324895, 1.0798898539187123, 0.41357166795318245));
        const r = query.compute(ln, circle);
        expect(Math.abs(r.distance - brute(ln, circle))).toBeLessThan(1e-6);
        expect(r.distance).toBeLessThan(0.02);
    }, 30000);

    it('agrees with brute force for long direction vectors', () => {
        // Scaling the direction scales a2 = |NxM|^2 quadratically, so this
        // sweep spends most of its samples where a2 > 1.
        const rnd = seededRandom(0x5eed247);
        for (let iter = 0; iter < 400; ++iter) {
            const nrm = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (length(nrm) < 0.3) { continue; }
            normalize(nrm);
            const circle = Circle3.fromCenterNormalRadius(v(0, 0, 0), nrm,
                0.5 + 2.5 * rnd());
            const d = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (length(d) < 0.3) { continue; }
            normalize(d);
            const ln = Line.fromOriginDirection(
                v((rnd() * 2 - 1) * 4, (rnd() * 2 - 1) * 4,
                    (rnd() * 2 - 1) * 4),
                mul(1 + 8 * rnd(), d));
            const r = query.compute(ln, circle);
            const n = 8192;
            const seed = Math.abs(nrm.values[0]) < 0.9
                ? v(1, 0, 0) : v(0, 1, 0);
            const u = sub(seed, mul(dot(seed, nrm), nrm));
            normalize(u);
            const w = cross(nrm, u);
            let best = Number.POSITIVE_INFINITY;
            for (let i = 0; i < n; ++i) {
                const a = (2 * Math.PI * i) / n;
                const p = add(mul(circle.radius * Math.cos(a), u),
                    mul(circle.radius * Math.sin(a), w));
                best = Math.min(best, pointLineDistance(p, ln));
            }
            expect(r.distance).toBeLessThanOrEqual(best + 1e-9);
            expect(r.distance).toBeGreaterThan(best - 1e-2);
        }
    }, 30000);
});
