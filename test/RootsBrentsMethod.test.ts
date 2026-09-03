import { describe, it, expect } from 'vitest';
import { RootsBrentsMethod } from '../src/RootsBrentsMethod.js';

// Zero tolerances: iterate until F(t) is exactly zero or the bracket
// collapses to consecutive floating-point numbers.
const tight = { maxIter: 200, negF: 0, posF: 0, stepT: 0, convT: 0 } as const;

function findTight(F: (t: number) => number, t0: number, t1: number) {
    return RootsBrentsMethod.find(F, t0, t1, tight.maxIter, tight.negF,
        tight.posF, tight.stepT, tight.convT);
}

describe('RootsBrentsMethod', () => {
    it('finds the root of cos on [0,3] to machine precision', () => {
        const { found, root } = findTight(Math.cos, 0, 3);
        expect(found).toBe(true);
        expect(Math.abs(root - Math.PI / 2)).toBeLessThanOrEqual(Number.EPSILON * 4);
    });

    it('finds the root of a cubic polynomial with tight tolerances', () => {
        // p(t) = (t - 1/3)(t^2 + 1), single real root 1/3.
        const F = (t: number) => (t - 1 / 3) * (t * t + 1);
        const { found, root } = findTight(F, -2, 2);
        expect(found).toBe(true);
        expect(Math.abs(root - 1 / 3)).toBeLessThanOrEqual(Number.EPSILON);
    });

    it('accepts a bracket endpoint that satisfies the function tolerance', () => {
        const F = (t: number) => t * (t - 2);
        // Root exactly at t0 = 0.
        let r = findTight(F, 0, 1);
        expect(r).toEqual({ found: true, root: 0 });
        // Root exactly at t1 = 2.
        r = findTight(F, 1, 2);
        expect(r).toEqual({ found: true, root: 2 });
    });

    it('uses the f-tolerances to accept near-root endpoints', () => {
        const F = (t: number) => t - 1;
        // F(0.9999) = -1e-4 is within [-1e-3, 1e-3].
        const r = RootsBrentsMethod.find(F, 0.9999, 2, 100, -1e-3, 1e-3, 0, 0);
        expect(r.found).toBe(true);
        expect(r.root).toBe(0.9999);
    });

    it('rejects an interval without a sign change', () => {
        const F = (t: number) => t * t + 1;
        expect(findTight(F, -1, 1).found).toBe(false);
    });

    it('rejects invalid parameters', () => {
        const F = (t: number) => t;
        // t1 <= t0.
        expect(RootsBrentsMethod.find(F, 1, 1, 100, 0, 0, 0, 0).found).toBe(false);
        expect(RootsBrentsMethod.find(F, 2, 1, 100, 0, 0, 0, 0).found).toBe(false);
        // maxIterations == 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 0, 0, 0, 0, 0).found).toBe(false);
        // negFTolerance > 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 100, 1e-8, 0, 0, 0).found).toBe(false);
        // posFTolerance < 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 100, 0, -1e-8, 0, 0).found).toBe(false);
        // stepTTolerance < 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 100, 0, 0, -1e-8, 0).found).toBe(false);
        // convTTolerance < 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 100, 0, 0, 0, -1e-8).found).toBe(false);
    });

    it('terminates via the interval-size tolerance', () => {
        const F = (t: number) => t - 1 / 3;
        const r = RootsBrentsMethod.find(F, 0, 1, 200, 0, 0, 0, 1e-6);
        expect(r.found).toBe(true);
        expect(Math.abs(r.root - 1 / 3)).toBeLessThanOrEqual(1e-6);
    });

    it('handles a pathological flat region (ninth-order contact)', () => {
        const F = (t: number) => Math.pow(t - 1, 9);
        const { found, root } = findTight(F, 0, 1.75);
        expect(found).toBe(true);
        expect(Math.abs(root - 1)).toBeLessThanOrEqual(1e-12);
    });

    it('handles a discontinuous sign function', () => {
        const x0 = 0.7234;
        const F = (t: number) => Math.sign(t - x0);
        const { found, root } = findTight(F, 0, 1);
        expect(found).toBe(true);
        expect(Math.abs(root - x0)).toBeLessThanOrEqual(1e-15);
    });

    it('converges much faster than bisection on a smooth function', () => {
        // Count function evaluations; Brent should need far fewer than the
        // ~50+ bisections required for machine-precision convergence.
        let count = 0;
        const F = (t: number) => { ++count; return Math.exp(t) - 2; };
        const { found, root } = findTight(F, 0, 1);
        expect(found).toBe(true);
        expect(Math.abs(root - Math.LN2)).toBeLessThanOrEqual(Number.EPSILON * 2);
        expect(count).toBeLessThan(30);
    });

    it('fails to converge when the iteration budget is too small', () => {
        // One iteration is not enough for exact-zero tolerances on a
        // transcendental function; the algorithm reports failure rather
        // than an insufficiently accurate root.
        const r = RootsBrentsMethod.find(Math.cos, 0, 3, 1, 0, 0, 0, 0);
        expect(r.found).toBe(false);
    });

    it('solves a linear function in a single secant step', () => {
        // The secant iterate is exact for linear F, so even a budget of one
        // iteration finds the root with exact-zero tolerances.
        const F = (t: number) => t - 1 / 3;
        const r = RootsBrentsMethod.find(F, 0, 1, 1, 0, 0, 0, 0);
        expect(r.found).toBe(true);
        expect(r.root).toBe(1 / 3);
    });
});
