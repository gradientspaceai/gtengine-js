import { describe, it, expect } from 'vitest';
import { Minimize1 } from '../src/Minimize1.js';

// A brute-force scan of F on [t0,t1] used as an independent estimate of the
// global minimum.
function scanMinimum(F: (t: number) => number, t0: number, t1: number,
    samples: number): { tMin: number; fMin: number } {
    let tMin = t0;
    let fMin = F(t0);
    for (let i = 1; i <= samples; ++i) {
        const t = t0 + (t1 - t0) * (i / samples);
        const f = F(t);
        if (f < fMin) {
            fMin = f;
            tMin = t;
        }
    }
    return { tMin, fMin };
}

describe('Minimize1', () => {
    it('rejects invalid subdivision and bisection counts', () => {
        expect(() => new Minimize1((t) => t, 0, 8)).toThrow('Invalid argument.');
        expect(() => new Minimize1((t) => t, 8, 0)).toThrow('Invalid argument.');
        expect(() => new Minimize1((t) => t, -1, -1)).toThrow('Invalid argument.');
    });

    it('rejects an initial guess outside the interval', () => {
        const minimizer = new Minimize1((t) => t * t, 8, 64);
        expect(() => minimizer.getMinimum(0, 1, -0.5)).toThrow('Invalid initial t value.');
        expect(() => minimizer.getMinimum(0, 1, 1.5)).toThrow('Invalid initial t value.');
        expect(() => minimizer.getMinimum(0, 1, 0)).not.toThrow();
        expect(() => minimizer.getMinimum(0, 1, 1)).not.toThrow();
    });

    it('clamps epsilon and tolerance to nonnegative values', () => {
        const minimizer = new Minimize1((t) => t * t, 8, 64, -1, -2);
        expect(minimizer.getEpsilon()).toBe(0);
        expect(minimizer.getTolerance()).toBe(0);

        minimizer.setEpsilon(1e-10);
        minimizer.setTolerance(1e-6);
        expect(minimizer.getEpsilon()).toBe(1e-10);
        expect(minimizer.getTolerance()).toBe(1e-6);

        minimizer.setEpsilon(-1e-10);
        minimizer.setTolerance(-1e-6);
        expect(minimizer.getEpsilon()).toBe(0);
        expect(minimizer.getTolerance()).toBe(0);
    });

    it('uses the midpoint of the interval as the default initial guess', () => {
        const samples: number[] = [];
        const minimizer = new Minimize1((t) => {
            samples.push(t);
            return t * t;
        }, 8, 64);
        minimizer.getMinimum(-2, 6);
        // The first three evaluations are at t0, the initial guess and t1.
        expect(samples[0]).toBe(-2);
        expect(samples[1]).toBe(2);
        expect(samples[2]).toBe(6);
    });

    it('finds the vertex of a parabola essentially exactly', () => {
        // F(t) = (t - 1)^2 + 3. The parabolic interpolation of three samples
        // of a quadratic is exact, so the vertex is located immediately.
        const F = (t: number) => (t - 1) * (t - 1) + 3;
        const minimizer = new Minimize1(F, 8, 64, 1e-12, 1e-12);
        const { tMin, fMin } = minimizer.getMinimum(-1, 4, 0.5);
        expect(tMin).toBeCloseTo(1, 10);
        expect(fMin).toBeCloseTo(3, 10);
    });

    it('finds the minimum of a convex quartic', () => {
        // F(t) = (t - 0.3)^4 + 2 has a single (flat) minimum at t = 0.3.
        const F = (t: number) => {
            const d = t - 0.3;
            return d * d * d * d + 2;
        };
        const minimizer = new Minimize1(F, 16, 128, 1e-14, 1e-8);
        const { tMin, fMin } = minimizer.getMinimum(-2, 5);
        expect(tMin).toBeCloseTo(0.3, 4);
        expect(fMin).toBeCloseTo(2, 10);
    });

    it('finds the minimum of a convex function whose minimum is an endpoint', () => {
        // F(t) = t is monotone, so no polyline is V-shaped and the search is
        // a pure subdivision. The smallest sample is the left endpoint.
        const minimizer = new Minimize1((t) => t, 6, 32);
        const { tMin, fMin } = minimizer.getMinimum(0, 1);
        expect(tMin).toBe(0);
        expect(fMin).toBe(0);
    });

    it('handles the symmetric case where the parabola vertex is the midpoint', () => {
        // F(t) = t^4 - 3*t^2 has minima at t = +/-sqrt(1.5) with value
        // -2.25. The initial polyline on [-3,3] is symmetric about t = 0, so
        // the parabola vertex is the middle sample and the neighborhood
        // subdivision branch of the bracketed search is exercised.
        const F = (t: number) => t * t * t * t - 3 * t * t;
        const minimizer = new Minimize1(F, 16, 128, 1e-14, 1e-10);
        const { tMin, fMin } = minimizer.getMinimum(-3, 3);
        expect(Math.abs(tMin)).toBeCloseTo(Math.sqrt(1.5), 6);
        expect(fMin).toBeCloseTo(-2.25, 10);
    });

    it('finds the global minimum of a multimodal quartic', () => {
        // F(t) = t^4 + t^3 - 3*t^2 has a local minimum near t = 0.906 with
        // value about -1.045 and a global minimum near t = -1.656 with value
        // about -5.249.
        const F = (t: number) => t * t * (t * t + t - 3);
        const minimizer = new Minimize1(F, 16, 128, 1e-14, 1e-10);
        const { tMin, fMin } = minimizer.getMinimum(-3, 3);
        const scan = scanMinimum(F, -3, 3, 1000000);
        expect(fMin).toBeLessThanOrEqual(scan.fMin + 1e-9);
        expect(tMin).toBeCloseTo(scan.tMin, 4);
        expect(fMin).toBeCloseTo(F(tMin), 12);
    });

    it('finds the global minimum of a multimodal trigonometric function', () => {
        // F(t) = sin(t) + 0.1*t on [0, 4*pi] has two local minima; the
        // smaller is the one nearer the left end because of the linear term.
        const F = (t: number) => Math.sin(t) + 0.1 * t;
        const minimizer = new Minimize1(F, 16, 128, 1e-14, 1e-10);
        const { tMin, fMin } = minimizer.getMinimum(0, 4 * Math.PI);
        const scan = scanMinimum(F, 0, 4 * Math.PI, 1000000);
        expect(fMin).toBeLessThanOrEqual(scan.fMin + 1e-9);
        expect(tMin).toBeCloseTo(scan.tMin, 4);
        // The derivative cos(t) + 0.1 vanishes at the minimum.
        expect(Math.cos(tMin) + 0.1).toBeCloseTo(0, 6);
    });

    it('honors the bisection budget for a slowly converging search', () => {
        // With a single bisection allowed, the search returns the best of the
        // few samples it has taken rather than a converged minimum.
        const F = (t: number) => (t - 0.3) * (t - 0.3);
        const coarse = new Minimize1(F, 1, 1, 0, 0);
        const fine = new Minimize1(F, 16, 128, 1e-14, 1e-12);
        const coarseResult = coarse.getMinimum(-2, 5);
        const fineResult = fine.getMinimum(-2, 5);
        expect(fineResult.fMin).toBeLessThanOrEqual(coarseResult.fMin);
        expect(fineResult.tMin).toBeCloseTo(0.3, 8);
    });

    it('respects the tolerance-based convergence test', () => {
        // A larger tolerance stops the bisection sooner, so the reported
        // minimum is no better than the one from a tighter tolerance.
        const F = (t: number) => Math.cosh(t - 1.25);
        const loose = new Minimize1(F, 8, 64, 1e-2, 1e-1);
        const tight = new Minimize1(F, 8, 64, 1e-14, 1e-12);
        const looseResult = loose.getMinimum(-4, 4);
        const tightResult = tight.getMinimum(-4, 4);
        expect(tightResult.fMin).toBeLessThanOrEqual(looseResult.fMin + 1e-15);
        expect(tightResult.tMin).toBeCloseTo(1.25, 6);
        expect(tightResult.fMin).toBeCloseTo(1, 12);
    });

    it('locates a minimum at an interior sample of a piecewise function', () => {
        // A V-shaped function with a corner at t = 0.5. The parabolic fit is
        // not exact, so the search relies on the bracketing logic.
        const F = (t: number) => Math.abs(t - 0.5);
        const minimizer = new Minimize1(F, 16, 200, 1e-16, 1e-14);
        const { tMin, fMin } = minimizer.getMinimum(-1, 2);
        expect(tMin).toBeCloseTo(0.5, 6);
        expect(fMin).toBeLessThan(1e-6);
    });
});
