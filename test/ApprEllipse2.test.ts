import { describe, it, expect } from 'vitest';
import { ApprEllipse2 } from '../src/ApprEllipse2.js';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import { Vector, dot } from '../src/Vector.js';
import {
    check, expectClose, fc, finite, seededRandom
} from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Points that lie exactly on the ellipse with the given center, orientation
// angle and semi-axis lengths (a along the rotated x-axis, b along the
// rotated y-axis).
function ellipsePoints(cx: number, cy: number, a: number, b: number,
    angle: number, count: number, phase: number = 0): Vector[] {
    const ca = Math.cos(angle), sa = Math.sin(angle);
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        const t = phase + 2 * Math.PI * i / count;
        const x = a * Math.cos(t), y = b * Math.sin(t);
        points.push(v2(cx + ca * x - sa * y, cy + sa * x + ca * y));
    }
    return points;
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

// |Dot(u,v)| is 1 when the unit-length u and v are parallel up to sign.
function parallelUpToSign(u: Vector, v: Vector): number {
    return Math.abs(dot(u, v));
}

describe('ApprEllipse2', () => {
    it('recovers an axis-aligned ellipse from exact points', () => {
        const points = ellipsePoints(0, 0, 4, 2, 0, 64);
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(points, 256, false, ellipse);

        expect(error).toBeLessThan(1e-16);
        expect(ellipse.center.values[0]).toBeCloseTo(0, 8);
        expect(ellipse.center.values[1]).toBeCloseTo(0, 8);
        // The eigenvalues of M are sorted increasing, so extent[0] is the
        // major semi-axis length and extent[1] the minor one.
        expect(ellipse.extent.values[0]).toBeCloseTo(4, 8);
        expect(ellipse.extent.values[1]).toBeCloseTo(2, 8);
        expect(parallelUpToSign(ellipse.axis[0], v2(1, 0))).toBeCloseTo(1, 10);
        expect(parallelUpToSign(ellipse.axis[1], v2(0, 1))).toBeCloseTo(1, 10);
    });

    it('recovers a translated and rotated ellipse from exact points', () => {
        const angle = 0.4;
        const points = ellipsePoints(2, -1, 3, 1.5, angle, 64);
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(points, 256, false, ellipse);

        expect(error).toBeLessThan(1e-16);
        expect(ellipse.center.values[0]).toBeCloseTo(2, 8);
        expect(ellipse.center.values[1]).toBeCloseTo(-1, 8);
        expect(ellipse.extent.values[0]).toBeCloseTo(3, 8);
        expect(ellipse.extent.values[1]).toBeCloseTo(1.5, 8);

        const major = v2(Math.cos(angle), Math.sin(angle));
        const minor = v2(-Math.sin(angle), Math.cos(angle));
        expect(parallelUpToSign(ellipse.axis[0], major)).toBeCloseTo(1, 10);
        expect(parallelUpToSign(ellipse.axis[1], minor)).toBeCloseTo(1, 10);

        // The unit-length axes are orthogonal.
        expect(dot(ellipse.axis[0], ellipse.axis[1])).toBeCloseTo(0, 12);
        expect(dot(ellipse.axis[0], ellipse.axis[0])).toBeCloseTo(1, 12);

        // The implicit equation is satisfied by the input points.
        for (const p of points) {
            const d = [p.values[0] - ellipse.center.values[0],
                p.values[1] - ellipse.center.values[1]];
            let q = 0;
            for (let i = 0; i < 2; ++i) {
                const c = d[0] * ellipse.axis[i].values[0] +
                    d[1] * ellipse.axis[i].values[1];
                q += (c * c) / (ellipse.extent.values[i] *
                    ellipse.extent.values[i]);
            }
            expect(q).toBeCloseTo(1, 6);
        }
    });

    it('handles the circle special case', () => {
        const points = ellipsePoints(-2, 5, 3, 3, 0, 48);
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(points, 128, false, ellipse);

        expect(error).toBeLessThan(1e-16);
        expect(ellipse.center.values[0]).toBeCloseTo(-2, 8);
        expect(ellipse.center.values[1]).toBeCloseTo(5, 8);
        expect(ellipse.extent.values[0]).toBeCloseTo(3, 8);
        expect(ellipse.extent.values[1]).toBeCloseTo(3, 8);
        // The axes are still orthonormal even though the axis directions are
        // not unique for a circle.
        expect(dot(ellipse.axis[0], ellipse.axis[1])).toBeCloseTo(0, 12);
    });

    it('continues the minimization when the ellipse is the initial guess', () => {
        const points = ellipsePoints(1, 1, 5, 2, -0.7, 96);
        const fitter = new ApprEllipse2();
        const ellipse = new Hyperellipsoid(2);
        const error0 = fitter.compute(points, 1, false, ellipse);
        const error1 = fitter.compute(points, 3, true, ellipse);
        const error2 = fitter.compute(points, 200, true, ellipse);

        // Restarting from the previous result continues the descent, so the
        // error never increases (up to floating-point round-off once the
        // minimizer has converged).
        expect(error1).toBeLessThanOrEqual(error0 + 1e-20);
        expect(error2).toBeLessThanOrEqual(error1 + 1e-20);
        expect(error2).toBeLessThan(1e-16);
        expect(ellipse.center.values[0]).toBeCloseTo(1, 6);
        expect(ellipse.center.values[1]).toBeCloseTo(1, 6);
        expect(ellipse.extent.values[0]).toBeCloseTo(5, 6);
        expect(ellipse.extent.values[1]).toBeCloseTo(2, 6);
    });

    it('reports the initial-guess error when no iterations are requested', () => {
        const points = ellipsePoints(0, 0, 4, 2, 0, 40);
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(points, 0, false, ellipse);

        // The initial guess is the oriented bounding box of the points, so
        // the error is nonzero but the returned ellipse is well formed.
        expect(error).toBeGreaterThan(0);
        expect(Number.isFinite(error)).toBe(true);
        expect(ellipse.extent.values[0]).toBeGreaterThan(0);
        expect(ellipse.extent.values[1]).toBeGreaterThan(0);
        expect(dot(ellipse.axis[0], ellipse.axis[1])).toBeCloseTo(0, 12);
    });

    it('is robust to small noise', () => {
        const rnd = makeRandom(1234567);
        const exact = ellipsePoints(0.5, -0.25, 3, 2, 0.9, 200);
        const noisy = exact.map(p => v2(
            p.values[0] + 1e-3 * (2 * rnd() - 1),
            p.values[1] + 1e-3 * (2 * rnd() - 1)));
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(noisy, 256, false, ellipse);

        expect(error).toBeLessThan(1e-4);
        expect(ellipse.center.values[0]).toBeCloseTo(0.5, 2);
        expect(ellipse.center.values[1]).toBeCloseTo(-0.25, 2);
        expect(ellipse.extent.values[0]).toBeCloseTo(3, 2);
        expect(ellipse.extent.values[1]).toBeCloseTo(2, 2);
    });

    it('recovers random ellipses (aggregated)', () => {
        const rnd = makeRandom(24680);
        let worstExtentError = 0;
        let worstCenterError = 0;
        for (let trial = 0; trial < 12; ++trial) {
            const cx = 4 * rnd() - 2, cy = 4 * rnd() - 2;
            const a = 1 + 3 * rnd();
            const b = 0.5 + 0.5 * a * rnd();
            const angle = Math.PI * rnd();
            const points = ellipsePoints(cx, cy, a, b, angle, 60, rnd());
            const ellipse = new Hyperellipsoid(2);
            new ApprEllipse2().compute(points, 512, false, ellipse);
            const major = Math.max(a, b), minor = Math.min(a, b);
            worstExtentError = Math.max(worstExtentError,
                Math.abs(ellipse.extent.values[0] - major),
                Math.abs(ellipse.extent.values[1] - minor));
            worstCenterError = Math.max(worstCenterError,
                Math.abs(ellipse.center.values[0] - cx),
                Math.abs(ellipse.center.values[1] - cy));
        }
        expect(worstExtentError).toBeLessThan(1e-5);
        expect(worstCenterError).toBeLessThan(1e-5);
    });

    it('rejects an empty point set', () => {
        const ellipse = new Hyperellipsoid(2);
        expect(() => new ApprEllipse2().compute([], 4, false, ellipse))
            .toThrow(/no points/);
    });

    it('rejects points that are not 2D', () => {
        const ellipse = new Hyperellipsoid(2);
        expect(() => new ApprEllipse2().compute(
            [Vector.fromArray([1, 2, 3])], 4, false, ellipse)).toThrow(/2D/);
    });
});

// ---------------------------------------------------------------------------
// Property-based verification (VERIFYING.md).

describe('ApprEllipse2 verification', () => {
    // The error function F(C,M)/n that the minimizer reports, computed
    // independently of the quartic/quadratic line-minimization polynomials
    // that UpdateCenter and UpdateMatrix use to evaluate it.
    function meanSquaredAlgebraicError(points: readonly Vector[],
        ellipse: Hyperellipsoid): number {
        let error = 0;
        for (const P of points) {
            const d0 = P.values[0] - ellipse.center.values[0];
            const d1 = P.values[1] - ellipse.center.values[1];
            let a = -1;
            for (let k = 0; k < 2; ++k) {
                const u = ellipse.axis[k].values;
                const proj = d0 * u[0] + d1 * u[1];
                a += (proj * proj) /
                    (ellipse.extent.values[k] * ellipse.extent.values[k]);
            }
            error += a * a;
        }
        return error / points.length;
    }

    // Samples of an ellipse, optionally perturbed, as an arbitrary.
    const sampleSet = fc.tuple(finite(-4, 4), finite(-4, 4), finite(1, 4),
        finite(1, 4), finite(0, 3.2), fc.integer({ min: 8, max: 40 }),
        finite(0, 6));

    it('reports the error function value of the ellipse it returns', () => {
        // The value returned by compute() is produced by evaluating the
        // degree-4 (center) or degree-2 (matrix) polynomial that models F on
        // the descent line. That must agree with a direct evaluation of F at
        // the returned (C,M), which is an independent check of every
        // polynomial coefficient.
        const rnd = seededRandom(20260904);
        check(fc.tuple(sampleSet, fc.integer({ min: 0, max: 6 })),
            ([[cx, cy, a, b, angle, count, phase], numIterations]) => {
                const points = ellipsePoints(cx, cy, a, b, angle, count, phase)
                    .map(p => v2(p.values[0] + 0.02 * (2 * rnd() - 1),
                        p.values[1] + 0.02 * (2 * rnd() - 1)));
                const ellipse = new Hyperellipsoid(2);
                const error = new ApprEllipse2().compute(points, numIterations,
                    false, ellipse);
                // Relative tolerance: the polynomial is evaluated by Horner
                // from mean moments while the reference sums the residuals
                // directly, so the two differ by accumulated round-off.
                expectClose(error, meanSquaredAlgebraicError(points, ellipse),
                    1e-14, 1e-7);
            }, 60);
    });

    it('never increases the error as iterations are added', () => {
        // Each of the two half-steps minimizes F along a descent line and
        // falls back to the current F value when no positive root improves
        // it, so the reported error is non-increasing in numIterations.
        const rnd = seededRandom(777);
        check(sampleSet, ([cx, cy, a, b, angle, count, phase]) => {
            const points = ellipsePoints(cx, cy, a, b, angle, count, phase)
                .map(p => v2(p.values[0] + 0.05 * (2 * rnd() - 1),
                    p.values[1] + 0.05 * (2 * rnd() - 1)));
            let previous = Number.MAX_VALUE;
            for (const numIterations of [0, 1, 2, 4, 8]) {
                const ellipse = new Hyperellipsoid(2);
                const error = new ApprEllipse2().compute(points,
                    numIterations, false, ellipse);
                expect(error).toBeLessThanOrEqual(previous * (1 + 1e-9) + 1e-12);
                previous = error;
            }
        }, 40);
    });

    it('keeps an exact ellipse fixed when it is the initial guess', () => {
        check(sampleSet, ([cx, cy, a, b, angle, count, phase]) => {
            const points = ellipsePoints(cx, cy, a, b, angle, count, phase);
            const ca = Math.cos(angle), sa = Math.sin(angle);
            const guess = Hyperellipsoid.fromCenterAxisExtent(v2(cx, cy),
                [v2(ca, sa), v2(-sa, ca)], v2(a, b));
            const error = new ApprEllipse2().compute(points, 4, true, guess);
            // The samples satisfy F = 0 exactly, so the gradients vanish and
            // the minimizer cannot move away from the exact solution.
            expect(error).toBeLessThan(1e-18);
            expectClose(guess.center.values[0], cx, 1e-8, 1e-8);
            expectClose(guess.center.values[1], cy, 1e-8, 1e-8);
            const extents = [guess.extent.values[0], guess.extent.values[1]]
                .sort((u, v) => u - v);
            const expected = [a, b].sort((u, v) => u - v);
            expectClose(extents[0], expected[0], 1e-8, 1e-8);
            expectClose(extents[1], expected[1], 1e-8, 1e-8);
        }, 60);
    });

    it('is equivariant under rigid motions of the samples', () => {
        // Every step of the algorithm - the initial oriented box, the
        // gradients and the eigendecomposition - commutes with a rigid
        // motion, so the fit of the moved samples is the moved fit.
        const rnd = seededRandom(31337);
        check(fc.tuple(sampleSet, finite(0, 3.1), finite(-5, 5), finite(-5, 5)),
            ([[cx, cy, a, b, angle, count, phase], rot, tx, ty]) => {
                const points = ellipsePoints(cx, cy, a, b, angle, count, phase)
                    .map(p => v2(p.values[0] + 0.02 * (2 * rnd() - 1),
                        p.values[1] + 0.02 * (2 * rnd() - 1)));
                const cr = Math.cos(rot), sr = Math.sin(rot);
                const moved = points.map(p => v2(
                    cr * p.values[0] - sr * p.values[1] + tx,
                    sr * p.values[0] + cr * p.values[1] + ty));

                const e0 = new Hyperellipsoid(2);
                const err0 = new ApprEllipse2().compute(points, 4, false, e0);
                const e1 = new Hyperellipsoid(2);
                const err1 = new ApprEllipse2().compute(moved, 4, false, e1);

                // Tolerance: equivariance is exact in real arithmetic but
                // only approximate in floating point. The descent step length
                // is a root of a polynomial built from the sample moments, and
                // the matrix step is halved until Sylvester's criterion holds,
                // so a round-off difference in the rotated moments is
                // amplified over the iterations. 1e-4 relative is far tighter
                // than any structural (mis-ported coefficient) error would be.
                expectClose(err0, err1, 1e-9, 1e-2);
                expectClose(cr * e0.center.values[0] - sr * e0.center.values[1]
                    + tx, e1.center.values[0], 1e-5, 1e-4);
                expectClose(sr * e0.center.values[0] + cr * e0.center.values[1]
                    + ty, e1.center.values[1], 1e-5, 1e-4);
                expectClose(e0.extent.values[0], e1.extent.values[0], 1e-5, 1e-4);
                expectClose(e0.extent.values[1], e1.extent.values[1], 1e-5, 1e-4);
            }, 40);
    });

    it('rejects degenerate inputs', () => {
        const ellipse = new Hyperellipsoid(2);
        expect(() => new ApprEllipse2().compute([], 1, false, ellipse))
            .toThrow(/no points/);
        expect(() => new ApprEllipse2().compute(
            [Vector.fromArray([1, 2, 3])], 1, false, ellipse)).toThrow(/2D/);
    });
});
