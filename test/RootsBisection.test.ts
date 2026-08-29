import { describe, it, expect } from 'vitest';
import { RootsBisection } from '../src/RootsBisection';

describe('RootsBisection', () => {
    it('finds the root of cos on [0,2] to machine precision', () => {
        const { iterations, root } = RootsBisection.find(Math.cos, 0, 2, 2048);
        expect(iterations).toBeGreaterThanOrEqual(2);
        expect(iterations).toBeLessThanOrEqual(2048);
        expect(Math.abs(root - Math.PI / 2)).toBeLessThanOrEqual(Number.EPSILON * 2);
    });

    it('converges to the floating-point neighbor of 1/3 with large iteration budget', () => {
        const F = (t: number) => t - 1 / 3;
        const { iterations, root } = RootsBisection.find(F, 0, 1, 200);
        // Bisection halts when the midpoint rounds to an interval endpoint
        // or F evaluates to exactly zero.
        expect(iterations).toBeLessThanOrEqual(60);
        expect(Math.abs(root - 1 / 3)).toBeLessThanOrEqual(Number.EPSILON);
    });

    it('returns 1 when the root is at the left bracket endpoint', () => {
        const F = (t: number) => t;
        const { iterations, root } = RootsBisection.find(F, 0, 1, 100);
        expect(iterations).toBe(1);
        expect(root).toBe(0);
    });

    it('returns 1 when the root is at the right bracket endpoint', () => {
        const F = (t: number) => t - 1;
        const { iterations, root } = RootsBisection.find(F, -3, 1, 100);
        expect(iterations).toBe(1);
        expect(root).toBe(1);
    });

    it('returns 0 when there is no sign change', () => {
        const F = (t: number) => t * t + 1;
        const { iterations } = RootsBisection.find(F, -1, 1, 100);
        expect(iterations).toBe(0);
    });

    it('returns 0 for an invalid interval (t0 >= t1)', () => {
        const F = (t: number) => t;
        expect(RootsBisection.find(F, 1, 1, 100).iterations).toBe(0);
        expect(RootsBisection.find(F, 2, 1, 100).iterations).toBe(0);
    });

    it('breaks early when the midpoint evaluates to exactly zero', () => {
        // F(t) = t - 0.25 on [0,1]: midpoints 0.5, then 0.25 where F = 0.
        const F = (t: number) => t - 0.25;
        const { iterations, root } = RootsBisection.find(F, 0, 1, 100);
        expect(root).toBe(0.25);
        expect(iterations).toBe(3);
    });

    it('returns maxIterations+1 when the budget is exhausted', () => {
        const F = (t: number) => t - 1 / 3;
        const { iterations, root } = RootsBisection.find(F, 0, 1, 5);
        expect(iterations).toBe(6);
        // The root estimate is still the best midpoint so far.
        expect(Math.abs(root - 1 / 3)).toBeLessThanOrEqual(1 / 32);
    });

    it('handles a flat region around the root (cubic contact)', () => {
        const F = (t: number) => Math.pow(t - 0.5, 3);
        const { root } = RootsBisection.find(F, 0.1, 1.3, 2048);
        expect(Math.abs(root - 0.5)).toBeLessThanOrEqual(1e-15);
    });

    it('handles a discontinuous sign function', () => {
        const x0 = 0.7234;
        const F = (t: number) => Math.sign(t - x0);
        const { root } = RootsBisection.find(F, 0, 1, 2048);
        expect(Math.abs(root - x0)).toBeLessThanOrEqual(Number.EPSILON);
    });

    describe('overload with known endpoint values', () => {
        it('accepts sign surrogates for infinite endpoint values', () => {
            // F(t) = 1/(t - 0.5) - 1/(t + 0.5) - like function with a root;
            // use F(t) = tan(t) on an interval straddling pi/2 asymmetric?
            // Instead: F(t) = 1/t - 4 on (0, 1]; F(0+) = +infinity. Pass
            // sign surrogates f0 = 1, f1 = F(1) = -3.
            const F = (t: number) => 1 / t - 4;
            const { iterations, root } = RootsBisection.find(F, 0, 1, 1, -3, 2048);
            expect(iterations).toBeGreaterThanOrEqual(2);
            expect(Math.abs(root - 0.25)).toBeLessThanOrEqual(1e-15);
        });

        it('returns 1 immediately when a supplied endpoint value is zero', () => {
            const F = (t: number) => t - 2;
            expect(RootsBisection.find(F, 2, 3, 0, 1, 100)).toEqual({ iterations: 1, root: 2 });
            expect(RootsBisection.find(F, 1, 2, -1, 0, 100)).toEqual({ iterations: 1, root: 2 });
        });

        it('returns 0 when supplied endpoint values have the same sign', () => {
            const F = (t: number) => t * t + 1;
            expect(RootsBisection.find(F, -1, 1, 2, 2, 100).iterations).toBe(0);
        });

        it('returns 0 for an invalid interval', () => {
            const F = (t: number) => t;
            expect(RootsBisection.find(F, 1, 0, -1, 1, 100).iterations).toBe(0);
        });

        it('matches the two-argument overload on a regular bracket', () => {
            const F = (t: number) => Math.exp(t) - 2;
            const a = RootsBisection.find(F, 0, 1, 2048);
            const b = RootsBisection.find(F, 0, 1, F(0), F(1), 2048);
            expect(b.root).toBe(a.root);
            expect(b.iterations).toBe(a.iterations);
            expect(Math.abs(a.root - Math.LN2)).toBeLessThanOrEqual(Number.EPSILON * 2);
        });
    });
});
