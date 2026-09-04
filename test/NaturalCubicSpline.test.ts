import { describe, it, expect } from 'vitest';
import { NaturalCubicSpline } from '../src/NaturalCubicSpline.js';
import { Vector, length, sub } from '../src/Vector.js';

function vec(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// The points of a smooth 3D test curve sampled at the given times.
function samplePoints(times: readonly number[],
    f: (t: number) => number[]): Vector[] {
    return times.map(t => vec(...f(t)));
}

const TIMES = [0, 1, 2.5, 4, 5];

describe('NaturalCubicSpline', () => {
    it('rejects invalid input', () => {
        expect(() => NaturalCubicSpline.createFree(
            [vec(0, 0), vec(1, 1)], [0, 1])).toThrow();
        expect(() => NaturalCubicSpline.createFree(
            [vec(0, 0), vec(1, 1), vec(2, 0)], [0, 1])).toThrow();
        expect(() => NaturalCubicSpline.createFree(
            [vec(0, 0), vec(1, 1, 1), vec(2, 0)], [0, 1, 2])).toThrow();
    });

    it('reports its construction state, domain and segment count', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(t)]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        expect(spline.isConstructed()).toBe(true);
        expect(spline.getDimension()).toBe(2);
        expect(spline.getTMin()).toBe(0);
        expect(spline.getTMax()).toBe(5);
        expect(spline.getNumSegments()).toBe(TIMES.length - 1);
        expect(spline.getPolynomials().length).toBe(TIMES.length - 1);
        for (const poly of spline.getPolynomials()) {
            expect(poly.length).toBe(4);
        }
    });

    it('interpolates the knots for all three boundary conditions', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(t), Math.cos(t)]);
        const d0 = vec(1, 1, 0);
        const d1 = vec(1, Math.cos(5), -Math.sin(5));
        const splines = [
            NaturalCubicSpline.createFree(points, TIMES),
            NaturalCubicSpline.createClosed(points, TIMES),
            NaturalCubicSpline.createClamped(points, TIMES, d0, d1)
        ];
        for (const spline of splines) {
            const jet = spline.createJet();
            for (let i = 0; i < TIMES.length; ++i) {
                spline.evaluate(TIMES[i], 0, jet);
                expect(length(sub(jet[0], points[i]))).toBeLessThan(1e-10);
            }
        }
    });

    it('is C2 continuous at the interior knots', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(2 * t)]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        const jetL = spline.createJet();
        const jetR = spline.createJet();
        const h = 1e-7;
        for (let i = 1; i < TIMES.length - 1; ++i) {
            const t = TIMES[i];
            spline.evaluate(t - h, 2, jetL);
            spline.evaluate(t + h, 2, jetR);
            for (let order = 0; order <= 2; ++order) {
                expect(length(sub(jetL[order], jetR[order]))).toBeLessThan(1e-5);
            }
        }
    });

    it('has zero second derivative at the ends for a free spline', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(2 * t), t * t]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        const jet = spline.createJet();
        spline.evaluate(TIMES[0], 2, jet);
        expect(length(jet[2])).toBeLessThan(1e-10);
        spline.evaluate(TIMES[TIMES.length - 1], 2, jet);
        expect(length(jet[2])).toBeLessThan(1e-10);
    });

    it('matches the requested end derivatives for a clamped spline', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(2 * t)]);
        const d0 = vec(1, 2);
        const d1 = vec(1, 2 * Math.cos(10));
        const spline = NaturalCubicSpline.createClamped(points, TIMES, d0, d1);
        const jet = spline.createJet();
        spline.evaluate(TIMES[0], 1, jet);
        expect(length(sub(jet[1], d0))).toBeLessThan(1e-10);
        spline.evaluate(TIMES[TIMES.length - 1], 1, jet);
        expect(length(sub(jet[1], d1))).toBeLessThan(1e-10);
    });

    it('matches the first and second derivatives at the ends when closed', () => {
        const points = samplePoints(TIMES, t => [Math.cos(t), Math.sin(t)]);
        const spline = NaturalCubicSpline.createClosed(points, TIMES);
        const jet0 = spline.createJet();
        const jet1 = spline.createJet();
        spline.evaluate(TIMES[0], 2, jet0);
        spline.evaluate(TIMES[TIMES.length - 1], 2, jet1);
        expect(length(sub(jet0[1], jet1[1]))).toBeLessThan(1e-10);
        expect(length(sub(jet0[2], jet1[2]))).toBeLessThan(1e-10);
    });

    it('reproduces a cubic polynomial exactly when clamped', () => {
        const P = (t: number) => [
            1 + 2 * t - 0.5 * t * t + 0.25 * t * t * t,
            -3 + t * t * t,
            7 - t
        ];
        const D = (t: number) => [2 - t + 0.75 * t * t, 3 * t * t, -1];
        const points = samplePoints(TIMES, P);
        const spline = NaturalCubicSpline.createClamped(points, TIMES,
            vec(...D(0)), vec(...D(5)));
        const jet = spline.createJet();
        for (const t of [0, 0.3, 1, 1.7, 2.5, 3.3, 4, 4.9, 5]) {
            spline.evaluate(t, 3, jet);
            expect(length(sub(jet[0], vec(...P(t))))).toBeLessThan(1e-10);
            expect(length(sub(jet[1], vec(...D(t))))).toBeLessThan(1e-10);
            // The third derivative of the reproduced cubic is constant.
            expect(jet[3].get(0)).toBeCloseTo(1.5, 8);
            expect(jet[3].get(1)).toBeCloseTo(6, 8);
            expect(jet[3].get(2)).toBeCloseTo(0, 8);
        }
    });

    it('derivatives agree with finite differences of the position', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(2 * t), Math.exp(-t)]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        const jet = spline.createJet();
        const jetP = spline.createJet();
        const jetM = spline.createJet();
        const h = 1e-6;
        for (const t of [0.4, 1.6, 2.9, 3.5, 4.6]) {
            spline.evaluate(t, 2, jet);
            spline.evaluate(t + h, 1, jetP);
            spline.evaluate(t - h, 1, jetM);
            for (let i = 0; i < 3; ++i) {
                const d1 = (jetP[0].get(i) - jetM[0].get(i)) / (2 * h);
                const d2 = (jetP[1].get(i) - jetM[1].get(i)) / (2 * h);
                expect(jet[1].get(i)).toBeCloseTo(d1, 6);
                expect(jet[2].get(i)).toBeCloseTo(d2, 5);
            }
        }
    });

    it('clamps evaluation outside the time interval to the end segments', () => {
        const points = samplePoints(TIMES, t => [t, t * t]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        const jet = spline.createJet();
        spline.evaluate(-10, 0, jet);
        const before = jet[0].clone();
        spline.evaluate(TIMES[0], 0, jet);
        expect(length(sub(before, jet[0]))).toBeLessThan(1e-12);
        spline.evaluate(100, 0, jet);
        const after = jet[0].clone();
        spline.evaluate(TIMES[TIMES.length - 1], 0, jet);
        expect(length(sub(after, jet[0]))).toBeLessThan(1e-12);
    });

    it('interpolates random data with C2 continuity (randomized)', () => {
        let seed = 24680;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const n = 3 + Math.floor(6 * rand());
            const times: number[] = [0];
            for (let i = 1; i < n; ++i) {
                times.push(times[i - 1] + 0.25 + rand());
            }
            const points: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                points.push(vec(4 * rand() - 2, 4 * rand() - 2));
            }
            const spline = NaturalCubicSpline.createFree(points, times);
            const jet = spline.createJet();
            for (let i = 0; i < n; ++i) {
                spline.evaluate(times[i], 0, jet);
                expect(length(sub(jet[0], points[i]))).toBeLessThan(1e-8);
            }
            // Natural boundary conditions.
            spline.evaluate(times[0], 2, jet);
            expect(length(jet[2])).toBeLessThan(1e-8);
            spline.evaluate(times[n - 1], 2, jet);
            expect(length(jet[2])).toBeLessThan(1e-8);
            // C2 across interior knots.
            const h = 1e-7;
            const jetL = spline.createJet();
            const jetR = spline.createJet();
            for (let i = 1; i < n - 1; ++i) {
                spline.evaluate(times[i] - h, 2, jetL);
                spline.evaluate(times[i] + h, 2, jetR);
                expect(length(sub(jetL[2], jetR[2]))).toBeLessThan(1e-4);
            }
        }
    });
});


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). NaturalCubicSpline.h was read
// line by line against src/NaturalCubicSpline.ts. The block row-reduction and
// back-substitution are checked against the conditions that define the spline
// rather than against themselves: on segment i the curve is
// p_i(u) = c0 + u*c1 + u^2*c2 + u^3*c3 with u = (t - times[i])/delta[i], so
// interpolation, C1 and C2 joins and each factory's boundary condition are all
// exact linear identities on the coefficients returned by GetPolynomials.
import {
    check, fc, expectClose, expectVectorClose, wellScaledVector
} from './helpers/arbitraries.js';

interface SplineInput {
    points: Vector[];
    times: number[];
    deltas: number[];
}

const cubicInput = (dim: number, minPoints = 3, maxPoints = 7):
    fc.Arbitrary<SplineInput> =>
    fc.integer({ min: minPoints, max: maxPoints }).chain(n => fc.tuple(
        fc.array(wellScaledVector(dim, -5, 5), { minLength: n, maxLength: n }),
        fc.array(fc.double({ min: 0.3, max: 2, noNaN: true,
            noDefaultInfinity: true }),
        { minLength: n - 1, maxLength: n - 1 })))
        .map(([points, gaps]) => {
            const times = [0];
            for (const g of gaps) { times.push(times[times.length - 1] + g); }
            return { points, times, deltas: gaps };
        });

// The value and the first two derivatives with respect to t of segment i at
// its local parameter u.
function segmentJet(poly: readonly Vector[], delta: number, u: number):
    { p: Vector, d1: Vector, d2: Vector } {
    const n = poly[0].size;
    const p = new Vector(n);
    const d1 = new Vector(n);
    const d2 = new Vector(n);
    for (let k = 0; k < n; ++k) {
        const c0 = poly[0].values[k], c1 = poly[1].values[k];
        const c2 = poly[2].values[k], c3 = poly[3].values[k];
        p.values[k] = c0 + u * (c1 + u * (c2 + u * c3));
        d1.values[k] = (c1 + u * (2 * c2 + u * 3 * c3)) / delta;
        d2.values[k] = (2 * c2 + u * 6 * c3) / (delta * delta);
    }
    return { p, d1, d2 };
}

function expectCubicSplineIdentities(curve: NaturalCubicSpline,
    input: SplineInput): void {
    const polys = curve.getPolynomials();
    const { points, deltas } = input;
    expect(polys.length).toBe(points.length - 1);
    for (let i = 0; i < polys.length; ++i) {
        // Interpolation of the left knot of each segment.
        expectVectorClose(segmentJet(polys[i], deltas[i], 0).p, points[i],
            1e-9, 1e-9);
    }
    // Interpolation of the right end of the last segment.
    const last = polys.length - 1;
    expectVectorClose(segmentJet(polys[last], deltas[last], 1).p,
        points[points.length - 1], 1e-8, 1e-8);

    // C0, C1 and C2 joins at the interior knots.
    for (let i = 0; i + 1 < polys.length; ++i) {
        const left = segmentJet(polys[i], deltas[i], 1);
        const right = segmentJet(polys[i + 1], deltas[i + 1], 0);
        expectVectorClose(left.p, right.p, 1e-8, 1e-8);
        expectVectorClose(left.d1, right.d1, 1e-8, 1e-8);
        expectVectorClose(left.d2, right.d2, 1e-8, 1e-8);
    }
}

describe('NaturalCubicSpline verification', () => {
    it('free splines interpolate, are C2 and have zero end curvature', () => {
        check(cubicInput(3), input => {
            const curve = NaturalCubicSpline.createFree(input.points,
                input.times);
            expectCubicSplineIdentities(curve, input);
            const polys = curve.getPolynomials();
            const last = polys.length - 1;
            const zero = new Vector(3);
            expectVectorClose(segmentJet(polys[0], input.deltas[0], 0).d2,
                zero, 1e-9, 1e-9);
            expectVectorClose(
                segmentJet(polys[last], input.deltas[last], 1).d2, zero,
                1e-8, 1e-8);
        });
    });

    it('clamped splines interpolate, are C2 and match the end tangents', () => {
        check(fc.tuple(cubicInput(2), wellScaledVector(2, -4, 4),
            wellScaledVector(2, -4, 4)), ([input, d0, d1]) => {
            const curve = NaturalCubicSpline.createClamped(input.points,
                input.times, d0, d1);
            expectCubicSplineIdentities(curve, input);
            const polys = curve.getPolynomials();
            const last = polys.length - 1;
            expectVectorClose(segmentJet(polys[0], input.deltas[0], 0).d1, d0,
                1e-9, 1e-9);
            expectVectorClose(
                segmentJet(polys[last], input.deltas[last], 1).d1, d1,
                1e-8, 1e-8);
        });
    });

    it('closed splines interpolate, are C2 and are periodic at the wrap', () => {
        check(cubicInput(2), input => {
            const curve = NaturalCubicSpline.createClosed(input.points,
                input.times);
            expectCubicSplineIdentities(curve, input);
            const polys = curve.getPolynomials();
            const last = polys.length - 1;
            const begin = segmentJet(polys[0], input.deltas[0], 0);
            const end = segmentJet(polys[last], input.deltas[last], 1);
            expectVectorClose(begin.d1, end.d1, 1e-8, 1e-8);
            expectVectorClose(begin.d2, end.d2, 1e-8, 1e-8);
        });
    });

    it('evaluate agrees with the stored polynomial coefficients', () => {
        check(fc.tuple(cubicInput(3), fc.double({ min: 0, max: 1,
            noNaN: true, noDefaultInfinity: true })), ([input, u]) => {
            const curve = NaturalCubicSpline.createFree(input.points,
                input.times);
            const polys = curve.getPolynomials();
            for (let i = 0; i < polys.length; ++i) {
                const t = input.times[i] + u * input.deltas[i];
                const expected = segmentJet(polys[i], input.deltas[i], u);
                const jet = curve.createJet();
                curve.evaluate(t, 2, jet);
                if (u === 1 && i + 1 < polys.length) {
                    // t is exactly the next knot, which selects the next
                    // segment; the values agree by C2 continuity, already
                    // checked above.
                    continue;
                }
                expectVectorClose(jet[0], expected.p, 1e-9, 1e-9);
                expectVectorClose(jet[1], expected.d1, 1e-9, 1e-9);
                expectVectorClose(jet[2], expected.d2, 1e-9, 1e-9);
            }
        });
    });

    it('clamped interpolation reproduces cubic polynomials exactly', () => {
        // A clamped cubic spline is the unique C2 interpolant with the given
        // end tangents, so it must reproduce any cubic sampled at the knots.
        check(fc.tuple(cubicInput(2, 4, 6),
            fc.array(wellScaledVector(2, -2, 2), { minLength: 4, maxLength: 4 }),
            fc.double({ min: 0, max: 1, noNaN: true,
                noDefaultInfinity: true })), ([input, c, s]) => {
            const at = (t: number): Vector => Vector.fromArray([0, 1].map(k =>
                c[0].values[k] + t * (c[1].values[k] + t * (c[2].values[k] +
                    t * c[3].values[k]))));
            const der = (t: number): Vector => Vector.fromArray([0, 1].map(k =>
                c[1].values[k] + t * (2 * c[2].values[k] +
                    3 * t * c[3].values[k])));
            const times = input.times;
            const tmax = times[times.length - 1];
            const curve = NaturalCubicSpline.createClamped(times.map(at),
                times, der(times[0]), der(tmax));
            const t = s * tmax;
            const jet = curve.createJet();
            curve.evaluate(t, 1, jet);
            expectVectorClose(jet[0], at(t), 1e-8, 1e-8);
            expectVectorClose(jet[1], der(t), 1e-8, 1e-8);
        });
    });

    it('is equivariant under translation and scaling of the points', () => {
        check(fc.tuple(cubicInput(2), wellScaledVector(2, -6, 6),
            fc.double({ min: -3, max: 3, noNaN: true,
                noDefaultInfinity: true }).filter(s => Math.abs(s) > 0.1),
            fc.double({ min: 0, max: 1, noNaN: true,
                noDefaultInfinity: true })), ([input, shift, scale, u]) => {
            const base = NaturalCubicSpline.createFree(input.points,
                input.times);
            const moved = NaturalCubicSpline.createFree(
                input.points.map(p => Vector.fromArray([
                    scale * p.values[0] + shift.values[0],
                    scale * p.values[1] + shift.values[1]])), input.times);
            const t = u * input.times[input.times.length - 1];
            const a = base.createJet();
            const b = moved.createJet();
            base.evaluate(t, 1, a);
            moved.evaluate(t, 1, b);
            for (let k = 0; k < 2; ++k) {
                expectClose(b[0].values[k], scale * a[0].values[k] +
                    shift.values[k], 1e-8, 1e-8);
                expectClose(b[1].values[k], scale * a[1].values[k],
                    1e-8, 1e-8);
            }
        });
    });

    it('clamps evaluation outside the time interval', () => {
        check(fc.tuple(cubicInput(2), fc.double({ min: 0.001, max: 50,
            noNaN: true, noDefaultInfinity: true })), ([input, d]) => {
            const curve = NaturalCubicSpline.createFree(input.points,
                input.times);
            const jet = curve.createJet();
            curve.evaluate(input.times[0] - d, 0, jet);
            expectVectorClose(jet[0], input.points[0], 1e-9, 1e-9);
            curve.evaluate(input.times[input.times.length - 1] + d, 0, jet);
            expectVectorClose(jet[0], input.points[input.points.length - 1],
                1e-8, 1e-8);
        });
    });
});
