import { describe, it, expect } from 'vitest';
import { RootsBisection2 } from '../src/RootsBisection2';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// Enough iterations that each bisector terminates by exhausting the 53-bit
// mantissa rather than by the iteration cap.
const MAX_ITERATIONS = 4096;

// A 2x2 nonlinear system with an exact rational root. The unit circle
// F(x,y) = x^2 + y^2 - 1 meets the line G(x,y) = y - (0.1 x + 0.74) at
// (x,y) = (0.6, 0.8):  0.1 * 0.6 + 0.74 = 0.8 and 0.6^2 + 0.8^2 = 1.
const circle = (x: number, y: number): number => x * x + y * y - 1;
const line = (x: number, y: number): number => y - (0.1 * x + 0.74);

// Rectangle satisfying the bisector's requirements:
//   F(xMin,y) < 0 < F(xMax,y) for all y in [yMin,yMax], because
//   xMin^2 + yMax^2 = 0.8825 < 1 < 1.0525 = xMax^2 + yMin^2.
//   G(x,yMin) < 0 < G(x,yMax) for all x in [xMin,xMax], because the line
//   values 0.78..0.815 lie strictly inside [0.7, 0.85].
const X_MIN = 0.4, X_MAX = 0.75, Y_MIN = 0.7, Y_MAX = 0.85;

describe('RootsBisection2 construction', () => {
    it('rejects nonpositive maximum iteration counts', () => {
        expect(() => new RootsBisection2(0, 8)).toThrow('Invalid maximum iterations.');
        expect(() => new RootsBisection2(8, 0)).toThrow('Invalid maximum iterations.');
        expect(() => new RootsBisection2(-2, -2)).toThrow('Invalid maximum iterations.');
    });

    it('starts with the root bound guarantee intact', () => {
        expect(new RootsBisection2(8, 8).noGuaranteeForRootBound()).toBe(false);
    });
});

describe('RootsBisection2 interval validation', () => {
    it('rejects a misordered x-interval', () => {
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        expect(() => bisector.find(circle, line, X_MAX, X_MIN, Y_MIN, Y_MAX)).toThrow(
            'Invalid ordering of t-interval endpoints.');
    });

    it('rejects a misordered y-interval', () => {
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        expect(() => bisector.find(circle, line, X_MIN, X_MAX, Y_MAX, Y_MIN)).toThrow(
            'Invalid ordering of t-interval endpoints.');
    });
});

describe('RootsBisection2 known root', () => {
    it('finds the circle-line intersection to near machine precision', () => {
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const result = bisector.find(circle, line, X_MIN, X_MAX, Y_MIN, Y_MAX);

        expect(result.iterations).toBeGreaterThan(2);
        expect(result.xRoot).toBeCloseTo(0.6, 12);
        expect(result.yRoot).toBeCloseTo(0.8, 12);
        expect(Math.abs(result.fAtRoot)).toBeLessThan(1e-12);
        expect(Math.abs(result.gAtRoot)).toBeLessThan(1e-12);
        expect(bisector.noGuaranteeForRootBound()).toBe(false);
    });

    it('reports outputs consistent with the returned root', () => {
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const { xRoot, yRoot, fAtRoot, gAtRoot } =
            bisector.find(circle, line, X_MIN, X_MAX, Y_MIN, Y_MAX);
        expect(fAtRoot).toBe(circle(xRoot, yRoot));
        expect(gAtRoot).toBe(line(xRoot, yRoot));
    });

    it('finds the root of a hyperbola-line system', () => {
        // F(x,y) = x + y - 3 and G(x,y) = x*y - 2 meet at (1,2) and (2,1).
        // On [1.5,2.5]x[0.6,1.4] only (2,1) is present, and the branch
        // y = 2/x maps that x-range into [0.8, 4/3], strictly inside the
        // y-band.
        const xMin = 1.5, xMax = 2.5, yMin = 0.6, yMax = 1.4;
        const F = (x: number, y: number): number => x + y - 3;
        const G = (x: number, y: number): number => x * y - 2;

        // F(xMin,y) = y - 1.5 < 0 and F(xMax,y) = y - 0.5 > 0 for every y in
        // the band; G(x,yMin) = 0.6x - 2 < 0 and G(x,yMax) = 1.4x - 2 > 0
        // for every x in the range.
        expect(F(xMin, yMax)).toBeLessThan(0);
        expect(F(xMax, yMin)).toBeGreaterThan(0);
        expect(G(xMax, yMin)).toBeLessThan(0);
        expect(G(xMin, yMax)).toBeGreaterThan(0);

        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const { xRoot, yRoot } = bisector.find(F, G, xMin, xMax, yMin, yMax);

        expect(xRoot).toBeCloseTo(2, 12);
        expect(yRoot).toBeCloseTo(1, 12);
        expect(bisector.noGuaranteeForRootBound()).toBe(false);
    });
});

describe('RootsBisection2 convergence behavior', () => {
    it('improves accuracy as the x-iteration cap grows', () => {
        let previousError = Number.MAX_VALUE;
        for (const cap of [4, 8, 16, 32]) {
            const bisector = new RootsBisection2(cap, MAX_ITERATIONS);
            const { iterations, xRoot } =
                bisector.find(circle, line, X_MIN, X_MAX, Y_MIN, Y_MAX);
            expect(iterations).toBe(cap + 1);
            const error = Math.abs(xRoot - 0.6);
            expect(error).toBeLessThan(previousError);
            // The x-interval has length 0.35 and is halved each iteration
            // after the first, which counts the endpoint evaluations.
            expect(error).toBeLessThanOrEqual(0.35 / Math.pow(2, cap - 1));
            previousError = error;
        }
    });

    it('limits the y-accuracy when the y-iteration cap is small', () => {
        const coarse = new RootsBisection2(MAX_ITERATIONS, 6);
        const fine = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const coarseRoot = coarse.find(circle, line, X_MIN, X_MAX, Y_MIN, Y_MAX);
        const fineRoot = fine.find(circle, line, X_MIN, X_MAX, Y_MIN, Y_MAX);
        expect(Math.abs(coarseRoot.yRoot - 0.8)).toBeGreaterThan(
            Math.abs(fineRoot.yRoot - 0.8));
        expect(Math.abs(coarseRoot.yRoot - 0.8)).toBeLessThan(0.15 / Math.pow(2, 4));
    });
});

describe('RootsBisection2 sign-precondition failures', () => {
    it('reports zero iterations when F has the same sign at both x-ends', () => {
        // Shift the circle so F > 0 on the whole rectangle.
        const F = (x: number, y: number): number => x * x + y * y + 1;
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const result = bisector.find(F, line, X_MIN, X_MAX, Y_MIN, Y_MAX);
        expect(result.iterations).toBe(0);
        expect(result.xRoot).toBe(0);
        expect(result.fAtRoot).toBe(0);
        expect(bisector.noGuaranteeForRootBound()).toBe(true);
    });

    it('reports the failure when G has the same sign at both y-ends', () => {
        // G > 0 throughout, so every nested y-bisection fails.
        const G = (_x: number, y: number): number => y + 1;
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const result = bisector.find(circle, G, X_MIN, X_MAX, Y_MIN, Y_MAX);
        expect(bisector.noGuaranteeForRootBound()).toBe(true);
        // The failed y-bisection reports root 0, so the x-bisection solves
        // F(x,0) = x^2 - 1 = 0 on [0.4,0.75], which also has no sign change.
        expect(result.iterations).toBe(0);
        expect(result.yRoot).toBe(0);
        expect(result.gAtRoot).toBe(0);
    });

    it('still bisects in x when the y-bisection fails but F changes sign', () => {
        // G < 0 throughout, so mYRoot is forced to 0 and the x-bisection
        // solves F(x,0) = x^2 - 0.25 = 0 with root x = 0.5.
        const G = (_x: number, y: number): number => -(y + 1);
        const F = (x: number, y: number): number => x * x + y * y - 0.25;
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const result = bisector.find(F, G, 0.1, 0.9, Y_MIN, Y_MAX);
        // The first midpoint of [0.1,0.9] is exactly 0.5, an exact root, so
        // the bisection stops on its first iteration (numbered 2).
        expect(result.iterations).toBe(2);
        expect(result.xRoot).toBe(0.5);
        expect(result.yRoot).toBe(0);
        expect(bisector.noGuaranteeForRootBound()).toBe(true);
    });
});

describe('RootsBisection2 endpoint and iteration edge cases', () => {
    it('returns 2 iterations with zeroed outputs when xMaxIterations is 1', () => {
        // Upstream leaves the x-outputs unwritten here (see issue #84); the
        // RootsBisection1 port zeroes them, so xRoot and fAtRoot are 0 while
        // the y-outputs come from the last endpoint evaluation, x = xMax.
        const bisector = new RootsBisection2(1, MAX_ITERATIONS);
        const result = bisector.find(circle, line, X_MIN, X_MAX, Y_MIN, Y_MAX);
        expect(result.iterations).toBe(2);
        expect(result.xRoot).toBe(0);
        expect(result.fAtRoot).toBe(0);
        expect(result.yRoot).toBeCloseTo(0.1 * X_MAX + 0.74, 12);
        expect(bisector.noGuaranteeForRootBound()).toBe(false);
    });

    it('returns 2 iterations with a zeroed y-root when yMaxIterations is 1', () => {
        const bisector = new RootsBisection2(MAX_ITERATIONS, 1);
        const result = bisector.find(circle, line, X_MIN, X_MAX, Y_MIN, Y_MAX);
        // Every y-bisection reports 2 iterations with a zeroed root, so the
        // x-bisection solves F(x,0) = x^2 - 1 on [0.4,0.75]: no sign change.
        expect(result.yRoot).toBe(0);
        expect(result.gAtRoot).toBe(0);
        expect(result.iterations).toBe(0);
    });

    it('returns xMin with 1 iteration when F(xMin,*) is exactly zero', () => {
        // Choose the rectangle so the x-bisection's left endpoint is already
        // a root: the line gives y = 0.75 at x = 0.5 and 0.5^2 + 0.75^2 is
        // not 1, so build a tailored pair instead.
        const G = (x: number, y: number): number => y - 0.5 * x;
        const F = (x: number, y: number): number => (x - 1) * (y + 1);
        // On [1,3]x[0.25,1.75]: G(x,0.25) = 0.25 - 0.5x < 0 for x >= 1 and
        // G(x,1.75) = 1.75 - 0.5x > 0 for x <= 3, so each y-bisection works.
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const result = bisector.find(F, G, 1, 3, 0.25, 1.75);
        expect(result.iterations).toBe(1);
        expect(result.xRoot).toBe(1);
        expect(result.fAtRoot).toBe(0);
        // Upstream quirk: the y-outputs come from the last XFunction
        // evaluation, which is at x = xMax, not at the returned xRoot.
        expect(result.yRoot).toBeCloseTo(1.5, 12);
        expect(result.gAtRoot).not.toBe(G(result.xRoot, result.yRoot));
    });

    it('returns xMax with 1 iteration when F(xMax,*) is exactly zero', () => {
        const G = (x: number, y: number): number => y - 0.5 * x;
        const F = (x: number, y: number): number => (x - 3) * (y + 1);
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const result = bisector.find(F, G, 1, 3, 0.25, 1.75);
        expect(result.iterations).toBe(1);
        expect(result.xRoot).toBe(3);
        expect(result.fAtRoot).toBe(0);
        // Here the last evaluation is at x = xMax = xRoot, so the y-outputs
        // are consistent.
        expect(result.yRoot).toBeCloseTo(1.5, 12);
        expect(Math.abs(result.gAtRoot)).toBeLessThan(1e-15);
    });
});

describe('RootsBisection2 root-bound flag semantics', () => {
    it('forgets an early y-failure when a later y-bisection succeeds', () => {
        // Upstream assigns rather than accumulates the y-status, so a
        // y-bisection failure at x = xMin is erased by the success at the
        // subsequent evaluations. The port preserves the quirk.
        let calls = 0;
        const G = (_x: number, y: number): number => {
            // The first two calls (from the y-bisection at x = xMin) report
            // a positive value, so that y-bisection fails.
            ++calls;
            return calls <= 2 ? 1 : y - 0.5;
        };
        const F = (x: number, _y: number): number => x - 1.5;
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const result = bisector.find(F, G, 1, 3, 0.25, 1.75);
        expect(result.xRoot).toBeCloseTo(1.5, 12);
        expect(result.yRoot).toBeCloseTo(0.5, 12);
        expect(bisector.noGuaranteeForRootBound()).toBe(false);
    });

    it('keeps the flag set when the last y-bisection fails', () => {
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);
        const alwaysPositive = (_x: number, y: number): number => y + 1;
        bisector.find(circle, alwaysPositive, X_MIN, X_MAX, Y_MIN, Y_MAX);
        expect(bisector.noGuaranteeForRootBound()).toBe(true);
    });
});

describe('RootsBisection2 randomized cross-check', () => {
    it('agrees with an independent Newton solve on translated circle-line systems', () => {
        const rng = makeRng(20260830);
        const bisector = new RootsBisection2(MAX_ITERATIONS, MAX_ITERATIONS);

        for (let trial = 0; trial < 40; ++trial) {
            // Circle of radius r centered at the origin and the line
            // y = a*x + b, with a small slope so the rectangle preconditions
            // are easy to satisfy.
            const r = 0.8 + 0.4 * rng();
            const a = 0.2 * (2 * rng() - 1);
            // Target root angle in the first quadrant, away from the axes.
            const theta = 0.6 + 0.4 * rng();
            const xStar = r * Math.cos(theta);
            const yStar = r * Math.sin(theta);
            const b = yStar - a * xStar;

            const F = (x: number, y: number): number => x * x + y * y - r * r;
            const G = (x: number, y: number): number => y - (a * x + b);

            // Rectangle bracketing the root; verify the preconditions at the
            // corners, which suffice for these monotone functions.
            const xMin = xStar - 0.1, xMax = xStar + 0.1;
            const yMin = yStar - 0.1, yMax = yStar + 0.1;
            if (F(xMin, yMax) >= 0 || F(xMax, yMin) <= 0) {
                continue;
            }
            if (G(xMin, yMin) >= 0 || G(xMin, yMax) <= 0) {
                continue;
            }
            if (G(xMax, yMin) >= 0 || G(xMax, yMax) <= 0) {
                continue;
            }

            const result = bisector.find(F, G, xMin, xMax, yMin, yMax);
            expect(result.iterations).toBeGreaterThan(1);
            expect(result.xRoot).toBeCloseTo(xStar, 10);
            expect(result.yRoot).toBeCloseTo(yStar, 10);
            expect(Math.abs(result.fAtRoot)).toBeLessThan(1e-9);
            expect(Math.abs(result.gAtRoot)).toBeLessThan(1e-9);
        }
    });
});
