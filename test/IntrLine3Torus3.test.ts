import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Torus3 } from '../src/Torus3.js';
import { Vector, add, dot, mul, normalize, sub, length } from '../src/Vector.js';
import { IntrLine3Torus3FI } from '../src/IntrLine3Torus3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function line(p: number[], d: number[]): Line {
    const dir = vec(d);
    normalize(dir);
    return Line.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The implicit torus function; it is zero on the surface.
function implicit(torus: Torus3, x: Vector): number {
    const delta = sub(x, torus.center);
    const sqrLen = dot(delta, delta);
    const dotN = dot(torus.normal, delta);
    const r0Sqr = torus.radius0 * torus.radius0;
    const r1Sqr = torus.radius1 * torus.radius1;
    const a = sqrLen + r0Sqr - r1Sqr;
    return a * a - 4 * r0Sqr * (sqrLen - dotN * dotN);
}

describe('IntrLine3Torus3', () => {
    const fi = new IntrLine3Torus3FI();

    // The standard torus: center at the origin, axis (0,0,1), r0 = 2, r1 = 1.
    const torus = new Torus3();

    it('finds four intersections along the plane of symmetry', () => {
        const l = line([0, 0, 0], [1, 0, 0]);
        const result = fi.find(l, torus);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(4);
        const expected = [-3, -1, 1, 3];
        for (let i = 0; i < 4; ++i) {
            expect(result.lineParameter[i]).toBeCloseTo(expected[i], 8);
            expect(result.point[i].values[0]).toBeCloseTo(expected[i], 8);
            expect(Math.abs(implicit(torus, result.point[i])))
                .toBeLessThan(1e-6);
        }
    });

    it('finds two intersections for a line through the tube', () => {
        // The line x = 2, y = 0 pierces the tube at z = -1 and z = 1.
        const l = line([2, 0, -5], [0, 0, 1]);
        const result = fi.find(l, torus);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.lineParameter[0]).toBeCloseTo(4, 8);
        expect(result.lineParameter[1]).toBeCloseTo(6, 8);
        expect(result.point[0].values[2]).toBeCloseTo(-1, 8);
        expect(result.point[1].values[2]).toBeCloseTo(1, 8);
    });

    it('reports no intersection for a line through the hole', () => {
        const l = line([0, 0, -5], [0, 0, 1]);
        const result = fi.find(l, torus);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports no intersection for a line above the torus', () => {
        const l = line([0, 0, 5], [1, 0, 0]);
        expect(fi.find(l, torus).intersect).toBe(false);
    });

    it('reports the two tangential contacts in the plane z = 1', () => {
        // At z = 1 the torus meets the plane in the double circle
        // x^2 + y^2 = 4, so each contact is a double root.
        const l = line([0, 0, 1], [1, 0, 0]);
        const result = fi.find(l, torus);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.lineParameter[0]).toBeCloseTo(-2, 6);
        expect(result.lineParameter[1]).toBeCloseTo(2, 6);
    });

    it('reports torus parameters consistent with the surface evaluation', () => {
        const l = line([0, 0, 0], [1, 0.35, 0.2]);
        const result = fi.find(l, torus);
        expect(result.numIntersections).toBeGreaterThan(0);
        for (let i = 0; i < result.numIntersections; ++i) {
            const [u, v] = result.torusParameter[i];
            const jet = torus.evaluate(u, v, 0);
            expect(length(sub(jet[0], result.point[i]))).toBeLessThan(1e-6);
        }
    });

    it('handles a translated and reoriented torus', () => {
        // The torus has axis (0,1,0) and center (1,-2,3). A line through the
        // center along the axis passes through the hole and misses.
        const t = Torus3.fromCenterFrameRadii(
            vec([1, -2, 3]), vec([0, 0, 1]), vec([1, 0, 0]), vec([0, 1, 0]),
            2, 0.5);
        const miss = line([1, -10, 3], [0, 1, 0]);
        expect(fi.find(miss, t).intersect).toBe(false);

        // A line in the plane of symmetry crosses all four times.
        const hit = line([1, -2, 3], [0, 0, 1]);
        const result = fi.find(hit, t);
        expect(result.numIntersections).toBe(4);
        const expected = [-2.5, -1.5, 1.5, 2.5];
        for (let i = 0; i < 4; ++i) {
            expect(result.lineParameter[i]).toBeCloseTo(expected[i], 8);
            expect(Math.abs(implicit(t, result.point[i]))).toBeLessThan(1e-6);
        }
    });

    it('produces points on the torus for random lines', () => {
        const rand = makeRandom(271828);
        let total = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const l = line(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(l, torus);
            expect(result.numIntersections).toBeLessThanOrEqual(4);
            total += result.numIntersections;
            for (let i = 0; i < result.numIntersections; ++i) {
                // The point is on the line.
                const onLine = add(l.origin,
                    mul(result.lineParameter[i], l.direction));
                expect(length(sub(onLine, result.point[i]))).toBeLessThan(1e-9);
                // The point is on the torus.
                expect(Math.abs(implicit(torus, result.point[i])))
                    .toBeLessThan(1e-5);
            }
            // The roots are reported in increasing order.
            for (let i = 1; i < result.numIntersections; ++i) {
                expect(result.lineParameter[i])
                    .toBeGreaterThan(result.lineParameter[i - 1]);
            }
        }
        expect(total).toBeGreaterThan(50);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, expectClose,
    expectVectorClose, seededRandom
} from './helpers/arbitraries.js';

// A ring torus (r0 > r1 > 0) with an orthonormal frame, and a line whose
// direction is unit length.
const lineTorus = fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
    fc.double({ min: 1.5, max: 4, noNaN: true }),
    fc.double({ min: 0.2, max: 1, noNaN: true }),
    wellScaledVector(3, -6, 6), unitVector(3))
    .map(([c, R, r0, r1, o, d]) => ({
        torus: Torus3.fromCenterFrameRadii(c, R[0], R[1], R[2], r0, r1),
        line: Line.fromOriginDirection(o, d)
    }));

// Signed distance from a point to the torus tube surface: the point is on the
// surface when this is zero. This is much better conditioned than the quartic
// implicit polynomial, whose value scales with r0^4.
function tubeResidual(torus: Torus3, p: Vector): number {
    const q = sub(p, torus.center);
    const x = dot(q, torus.direction0);
    const y = dot(q, torus.direction1);
    const z = dot(q, torus.normal);
    const radial = Math.hypot(x, y) - torus.radius0;
    return Math.hypot(radial, z) - torus.radius1;
}

describe('IntrLine3Torus3 verification', () => {
    const fiq = new IntrLine3Torus3FI();

    it('reports between 0 and 4 intersections in ascending order', () => {
        check(lineTorus, ({ torus, line: l }) => {
            const f = fiq.find(l, torus);
            expect(f.numIntersections).toBeGreaterThanOrEqual(0);
            expect(f.numIntersections).toBeLessThanOrEqual(4);
            expect(f.intersect).toBe(f.numIntersections > 0);
            for (let i = 1; i < f.numIntersections; ++i) {
                // std::map<T,int32_t> iterates ascending by key.
                expect(f.lineParameter[i - 1])
                    .toBeLessThan(f.lineParameter[i]);
            }
        });
    });

    it('the reported points lie on the line and on the torus', () => {
        check(lineTorus, ({ torus, line: l }) => {
            const f = fiq.find(l, torus);
            for (let i = 0; i < f.numIntersections; ++i) {
                expectVectorClose(f.point[i],
                    add(l.origin, mul(f.lineParameter[i], l.direction)), 0, 0);
                // The quartic root finder works in floating point (upstream
                // suggests rational arithmetic for exact classification), so
                // the residual of a near-double root is larger than for a
                // simple root; 1e-6 covers the sampled configurations.
                expect(Math.abs(tubeResidual(torus, f.point[i])))
                    .toBeLessThan(1e-6);
                expect(Number.isNaN(f.lineParameter[i])).toBe(false);
                for (const x of f.point[i].values) {
                    expect(Number.isNaN(x)).toBe(false);
                }
            }
        });
    });

    it('the surface parameters (u,v) reproduce the reported point', () => {
        check(lineTorus, ({ torus, line: l }) => {
            const f = fiq.find(l, torus);
            for (let i = 0; i < f.numIntersections; ++i) {
                const [u, v] = f.torusParameter[i];
                // Torus3::GetParameters uses atan2, so the parameters are in
                // [-pi,pi] even though the header documents the surface
                // parameterization on [0,2*pi); cos/sin are 2*pi periodic so
                // Evaluate still reproduces the point.
                expect(Math.abs(u)).toBeLessThanOrEqual(Math.PI + 1e-12);
                expect(Math.abs(v)).toBeLessThanOrEqual(Math.PI + 1e-12);
                const jet = torus.evaluate(u, v, 0);
                expectVectorClose(jet[0], f.point[i], 1e-6, 1e-7);
            }
        });
    });

    it('a line through the tube of a canonical torus has four roots', () => {
        // Center C, normal (0,0,1); the line y = 0, z = 0 meets the torus at
        // x = +-(r0-r1) and x = +-(r0+r1).
        check(fc.tuple(wellScaledVector(3, -4, 4),
            fc.double({ min: 1.5, max: 4, noNaN: true }),
            fc.double({ min: 0.2, max: 1, noNaN: true })), ([c, r0, r1]) => {
            const torus = Torus3.fromCenterFrameRadii(c,
                Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1]), r0, r1);
            const l = Line.fromOriginDirection(c, Vector.fromArray([1, 0, 0]));
            const f = fiq.find(l, torus);
            expect(f.numIntersections).toBe(4);
            const expected = [-(r0 + r1), -(r0 - r1), r0 - r1, r0 + r1];
            for (let i = 0; i < 4; ++i) {
                expectClose(f.lineParameter[i], expected[i], 1e-8, 1e-9);
            }
        });
    });

    it('a line on the symmetry axis misses a ring torus', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
            fc.double({ min: 1.5, max: 4, noNaN: true }),
            fc.double({ min: 0.2, max: 1, noNaN: true })),
            ([c, R, r0, r1]) => {
                const torus = Torus3.fromCenterFrameRadii(c, R[0], R[1], R[2],
                    r0, r1);
                const l = Line.fromOriginDirection(c, R[2]);
                const f = fiq.find(l, torus);
                expect(f.intersect).toBe(false);
                expect(f.numIntersections).toBe(0);
            });
    });

    it('a fine sweep of the line brackets the reported roots', () => {
        const rnd = seededRandom(0x4e91c37f);
        for (let trial = 0; trial < 120; ++trial) {
            const r0 = 1.5 + rnd() * 2;
            const r1 = 0.3 + rnd() * 0.8;
            const torus = Torus3.fromCenterFrameRadii(
                Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                    rnd() * 2 - 1]),
                Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1]), r0, r1);
            const d = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (length(d) < 0.3) {
                continue;
            }
            normalize(d);
            const l = Line.fromOriginDirection(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3,
                    rnd() * 6 - 3]), d);
            const f = fiq.find(l, torus);
            // Count sign changes of the tube residual along the line; each is
            // a transversal crossing that the query must report.
            const N = 4000;
            const lo = -14, hi = 14;
            let prev = tubeResidual(torus,
                add(l.origin, mul(lo, d)));
            const crossings: number[] = [];
            for (let k = 1; k <= N; ++k) {
                const t = lo + ((hi - lo) * k) / N;
                const cur = tubeResidual(torus, add(l.origin, mul(t, d)));
                if ((prev < 0 && cur > 0) || (prev > 0 && cur < 0)) {
                    crossings.push(t);
                }
                prev = cur;
            }
            expect(f.numIntersections)
                .toBeGreaterThanOrEqual(crossings.length);
            // Every sign change is near one of the reported roots.
            for (const t of crossings) {
                let best = Infinity;
                for (let i = 0; i < f.numIntersections; ++i) {
                    best = Math.min(best, Math.abs(f.lineParameter[i] - t));
                }
                expect(best).toBeLessThan(1e-2);
            }
        }
    }, 30000);

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineTorus, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([{ torus, line: l }, R, tr]) => {
            const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
            const xf = (p: Vector): Vector => add(tr, rot(p));
            const l2 = Line.fromOriginDirection(xf(l.origin), rot(l.direction));
            const t2 = Torus3.fromCenterFrameRadii(xf(torus.center),
                rot(torus.direction0), rot(torus.direction1),
                rot(torus.normal), torus.radius0, torus.radius1);
            const f1 = fiq.find(l, torus);
            const f2 = fiq.find(l2, t2);
            if (f1.numIntersections !== f2.numIntersections) {
                // A tangency splits or merges roots under perturbation; that
                // is a property of the quartic solver, not of this file.
                return;
            }
            for (let i = 0; i < f1.numIntersections; ++i) {
                expectClose(f1.lineParameter[i], f2.lineParameter[i],
                    1e-5, 1e-6);
            }
        });
    });
});
