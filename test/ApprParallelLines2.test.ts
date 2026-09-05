import { describe, it, expect } from 'vitest';
import { ApprParallelLines2 } from '../src/ApprParallelLines2.js';
import { Vector, dot } from '../src/Vector.js';
import { check, expectClose, fc, finite, positive, vector } from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Samples that lie exactly on the two lines C -+ r*U + s*V, where
// V = (cos(angle),sin(angle)) and U = (-sin(angle),cos(angle)).
function twoLineSamples(angle: number, cx: number, cy: number, r: number,
    parameters: readonly number[]): Vector[] {
    const V = v2(Math.cos(angle), Math.sin(angle));
    const U = v2(-Math.sin(angle), Math.cos(angle));
    const points: Vector[] = [];
    for (const sign of [-1, 1]) {
        for (const s of parameters) {
            points.push(v2(
                cx + sign * r * U.values[0] + s * V.values[0],
                cy + sign * r * U.values[1] + s * V.values[1]));
        }
    }
    return points;
}

// The quantity the algorithm minimizes: the average of (u^2 - r^2)^2, where
// u is the component of the sample along the normal U of the two lines. It
// is zero exactly when every sample lies on one of the two lines.
function fitError(points: readonly Vector[], center: Vector,
    direction: Vector, radius: number): number {
    const U = v2(-direction.values[1], direction.values[0]);
    let sum = 0;
    for (const p of points) {
        const diff = v2(p.values[0] - center.values[0],
            p.values[1] - center.values[1]);
        const u = dot(diff, U);
        const d = u * u - radius * radius;
        sum += d * d;
    }
    return sum / points.length;
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprParallelLines2', () => {
    it('recovers two axis-aligned parallel lines', () => {
        // The lines y = 3 and y = -3, so V = (1,0), C = (0,0), radius 3.
        const points = twoLineSamples(0, 0, 0, 3, [-4, -2, -1, 0, 1, 2, 4]);
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 128);

        expect(Math.abs(result.direction.values[0])).toBeCloseTo(1, 8);
        expect(result.direction.values[1]).toBeCloseTo(0, 8);
        expect(result.radius).toBeCloseTo(3, 8);
        expect(result.center.values[1]).toBeCloseTo(0, 8);
        expect(fitError(points, result.center, result.direction,
            result.radius)).toBeCloseTo(0, 10);
    });

    it('recovers two rotated parallel lines', () => {
        const angle = 0.7;
        const points = twoLineSamples(angle, 2, -1, 1.5,
            [-3, -2, -1, 0, 1, 2, 3]);
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 128);

        const V = v2(Math.cos(angle), Math.sin(angle));
        expect(Math.abs(dot(result.direction, V))).toBeCloseTo(1, 6);
        expect(result.radius).toBeCloseTo(1.5, 6);
        expect(fitError(points, result.center, result.direction,
            result.radius)).toBeCloseTo(0, 8);
    });

    it('returns a unit-length direction with the center on the normal axis', () => {
        const angle = -0.35;
        const points = twoLineSamples(angle, -4, 5, 0.75, [-2, -1, 0, 1, 2, 5]);
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 128);

        // V is unit length and C has no V-component (only a U-component).
        expect(dot(result.direction, result.direction)).toBeCloseTo(1, 10);
        expect(dot(result.center, result.direction)).toBeCloseTo(0, 8);
        expect(result.radius).toBeGreaterThan(0);
    });

    it('is invariant under translation of the samples along the lines', () => {
        const angle = 0.9;
        const base = twoLineSamples(angle, 0, 0, 2, [-2, -1, 0, 1, 2, 3]);
        const V = v2(Math.cos(angle), Math.sin(angle));
        const shifted = base.map(p => v2(p.values[0] + 7 * V.values[0],
            p.values[1] + 7 * V.values[1]));

        const fitter = new ApprParallelLines2();
        const a = fitter.fit(base, 1024);
        const b = fitter.fit(shifted, 1024);
        expect(b.radius).toBeCloseTo(a.radius, 6);
        expect(Math.abs(dot(a.direction, b.direction))).toBeCloseTo(1, 6);
    });

    it('recovers the separation of two clusters of noisy samples', () => {
        const random = makeRandom(20260830);
        const angle = 0.25;
        const V = v2(Math.cos(angle), Math.sin(angle));
        const U = v2(-Math.sin(angle), Math.cos(angle));
        const points: Vector[] = [];
        for (let i = 0; i < 40; ++i) {
            const sign = (i % 2 === 0) ? 1 : -1;
            const s = 10 * random() - 5;
            const noise = 0.01 * (2 * random() - 1);
            const u = sign * 2 + noise;
            points.push(v2(3 + u * U.values[0] + s * V.values[0],
                1 + u * U.values[1] + s * V.values[1]));
        }
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 2048);
        expect(Math.abs(dot(result.direction, V))).toBeCloseTo(1, 2);
        expect(result.radius).toBeCloseTo(2, 2);
    });

    it('does not do worse than the trivial (0,1) starting estimate', () => {
        // The algorithm seeds the minimization with sigma = 0, gamma = 1 and
        // only replaces it when a candidate has a smaller error, so the
        // reported fit is at least as good as the seed for every input.
        const random = makeRandom(13579);
        for (let trial = 0; trial < 10; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 24; ++i) {
                points.push(v2(10 * random() - 5, 10 * random() - 5));
            }
            const fitter = new ApprParallelLines2();
            const result = fitter.fit(points, 128);

            // The seed fit: direction (1,0) with the center and radius
            // computed from the moments of the mean-centered samples.
            let ax = 0, ay = 0;
            for (const p of points) {
                ax += p.values[0];
                ay += p.values[1];
            }
            ax /= points.length;
            ay /= points.length;
            let Z02 = 0, Z03 = 0;
            for (const p of points) {
                const y = p.values[1] - ay;
                Z02 += y * y;
                Z03 += y * y * y;
            }
            Z02 /= points.length;
            Z03 /= points.length;
            const k = Z03 / (2 * Z02);
            const seedCenter = v2(0, ay + k);
            const seedRadius = Math.sqrt(k * k + Z02);
            const seedError = fitError(points, seedCenter, v2(1, 0), seedRadius);

            const error = fitError(points, result.center, result.direction,
                result.radius);
            expect(error).toBeLessThanOrEqual(seedError + 1e-8);
            expect(dot(result.direction, result.direction)).toBeCloseTo(1, 8);
        }
    });

    it('recovers unevenly sampled lines with nonzero odd moments', () => {
        // Five samples on one line and three on the other, at an angle: the
        // mean-centered data is not centrally symmetric, so the odd moments
        // Z12, Z30, Z21 and Z03 that drive the polynomial F are all nonzero.
        const angle = 0.45;
        const V = v2(Math.cos(angle), Math.sin(angle));
        const U = v2(-Math.sin(angle), Math.cos(angle));
        const points: Vector[] = [];
        const push = (sign: number, s: number) => {
            points.push(v2(
                -2 + sign * 1.25 * U.values[0] + s * V.values[0],
                3 + sign * 1.25 * U.values[1] + s * V.values[1]));
        };
        for (const s of [-2, -1, 0, 1.5, 4]) {
            push(1, s);
        }
        for (const s of [-3, 0.5, 2]) {
            push(-1, s);
        }

        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 2048);
        expect(Math.abs(dot(result.direction, V))).toBeCloseTo(1, 6);
        expect(result.radius).toBeCloseTo(1.25, 6);
        expect(fitError(points, result.center, result.direction,
            result.radius)).toBeCloseTo(0, 6);
    });

    it('handles isotropic samples, for which no direction is preferred', () => {
        // Samples spread uniformly on a circle have Z11 = 0, Z20 = Z02,
        // vanishing odd moments and Z40 = Z04 = 3*Z22, so f0 and f1 are the
        // zero polynomial in exact arithmetic and the error is the same for
        // every direction. Whichever direction is reported, the geometry
        // must still be consistent: r^2 = average(u^2) = 1/2.
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            const angle = 2 * Math.PI * i / 8;
            points.push(v2(Math.cos(angle), Math.sin(angle)));
        }
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 128);
        expect(dot(result.direction, result.direction)).toBeCloseTo(1, 8);
        expect(result.radius).toBeCloseTo(Math.SQRT1_2, 8);
        expect(Number.isFinite(result.center.values[0])).toBe(true);
        expect(Number.isFinite(result.center.values[1])).toBe(true);
    });

    it('handles four points at the corners of a rectangle', () => {
        // The two long sides of the rectangle are the parallel lines.
        const points = [v2(-3, 1), v2(3, 1), v2(-3, -1), v2(3, -1)];
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 128);
        expect(result.radius).toBeCloseTo(1, 8);
        expect(Math.abs(result.direction.values[0])).toBeCloseTo(1, 8);
        expect(fitError(points, result.center, result.direction,
            result.radius)).toBeCloseTo(0, 10);
    });
});

describe('ApprParallelLines2 verification', () => {
    // The error function the algorithm minimizes, evaluated for the line
    // direction V(angle): E = average(((u-k)^2 - r^2)^2) with
    // u = Dot(P - A, U), k = S30/(2*S20) and r^2 = k^2 + S20. Expanding the
    // square gives upstream's S40 - 4*k*S30 + (4*k^2 - S20)*S20.
    const modelError = (points: readonly Vector[], angle: number): number => {
        const n = points.length;
        let ax = 0, ay = 0;
        for (const p of points) { ax += p.get(0); ay += p.get(1); }
        ax /= n; ay /= n;
        const ux = -Math.sin(angle), uy = Math.cos(angle);
        let s20 = 0, s30 = 0, s40 = 0;
        for (const p of points) {
            const u = (p.get(0) - ax) * ux + (p.get(1) - ay) * uy;
            const u2 = u * u;
            s20 += u2; s30 += u2 * u; s40 += u2 * u2;
        }
        s20 /= n; s30 /= n; s40 /= n;
        const k = s30 / (2 * s20);
        return s40 - 4 * k * s30 + (4 * k * k - s20) * s20;
    };

    // The same error for the fitted model, using its own center and radius.
    const resultError = (points: readonly Vector[], center: Vector,
        direction: Vector, radius: number): number => {
        const ux = -direction.get(1), uy = direction.get(0);
        let sum = 0;
        for (const p of points) {
            const u = (p.get(0) - center.get(0)) * ux
                + (p.get(1) - center.get(1)) * uy;
            const d = u * u - radius * radius;
            sum += d * d;
        }
        return sum / points.length;
    };

    // Samples exactly on the two lines C -+ r*U + s*V.
    const twoLines = (angle: number, c: Vector, r: number,
        ss: readonly number[]): Vector[] => {
        const vx = Math.cos(angle), vy = Math.sin(angle);
        const ux = -vy, uy = vx;
        const points: Vector[] = [];
        for (const sign of [-1, 1]) {
            for (const s of ss) {
                points.push(Vector.fromArray([
                    c.get(0) + sign * r * ux + s * vx,
                    c.get(1) + sign * r * uy + s * vy]));
            }
        }
        return points;
    };

    // At least four well-separated positions along the lines. Fewer than
    // three distinct positions would let a second family of parallel lines
    // (perpendicular to V) fit the samples exactly too, so the minimizer
    // would not be unique.
    const positionsArb = fc.array(finite(0.8, 2),
        { minLength: 4, maxLength: 8 })
        .map(gaps => {
            let acc = -3;
            return gaps.map(g => (acc += g));
        });

    // sigma = sin(angle) is produced by a square root, so the fitted
    // direction always has a nonnegative second component; angles are drawn
    // from (0, pi) to make the recovered direction comparable to V(angle).
    const configArb = fc.tuple(
        finite(0.2, Math.PI - 0.2), vector(2, -5, 5), positive(4, 0.5),
        positionsArb);

    it('recovers two parallel lines its samples lie on', () => {
        check(configArb, ([angle, c, r, ss]) => {
            const points = twoLines(angle, c, r, ss);
            const result = new ApprParallelLines2().fit(points, 1024);

            const V = Vector.fromArray([Math.cos(angle), Math.sin(angle)]);
            // RootsPolynomial bisects to a fixed iteration budget, so the
            // recovered root carries a relative error of a few 1e-5.
            expectClose(Math.abs(dot(result.direction, V)), 1, 1e-4, 1e-4);
            expectClose(result.radius, r, 1e-4, 1e-4);

            // The fit is exact, so the residual is zero up to the root
            // accuracy of the bisected degree-16 polynomial (observed up to
            // ~3e-6 relative to the scale).
            const scale = points.reduce((u, p) => u + dot(p, p), 0);
            expect(resultError(points, result.center, result.direction,
                result.radius)).toBeLessThanOrEqual(1e-7 * (1 + scale));
        });
    });

    it('returns a unit direction and a center with no V-component', () => {
        // gamma is computed as -f0(sigma^2)/(sigma*f1(sigma^2)) rather than
        // from gamma^2 + sigma^2 = 1, so V is unit length only to the
        // accuracy of the root of h. That is fine for well-conditioned data
        // such as these two-line samples; near-degenerate samples can make
        // upstream return a badly non-unit direction (see the upstream bug
        // notes for this file).
        check(configArb, ([angle, c, r, ss]) => {
            const result = new ApprParallelLines2()
                .fit(twoLines(angle, c, r, ss), 1024);
            // Observed |V|^2 - 1 up to ~2.2e-6 on well-conditioned samples
            // (the bisection stops at 1024 iterations or the root tolerance),
            // so the bound is 1e-5: still far below the O(1) errors of #380.
            expectClose(dot(result.direction, result.direction), 1,
                1e-5, 1e-5);
            // Center V-component: same root-accuracy bound (observed 1.1e-5).
            expectClose(dot(result.center, result.direction), 0, 1e-4, 1e-4);
            expect(result.radius).toBeGreaterThan(0);
            expect(Number.isFinite(result.radius)).toBe(true);
        });
    });

    it('never returns a worse model than a sampled direction', () => {
        // The candidate roots come from a degree-16 polynomial; the returned
        // (sigma,gamma) must minimize the error at least as well as any other
        // direction does. Upstream's ComputeF sets a30[1] = -3, dropping the
        // factor Z12 that the equivalent evaluation in UpdateParameters uses;
        // that defect (fixed in the port) makes the returned direction miss
        // the minimum, which this property detects.
        check(fc.tuple(configArb, fc.array(finite(-0.4, 0.4),
            { minLength: 18, maxLength: 18 })),
            ([[angle, c, r, ss], noise]) => {
                const points = twoLines(angle, c, r, ss).map((p, i) =>
                    Vector.fromArray([p.get(0) + noise[i % noise.length],
                        p.get(1) + noise[(i + 7) % noise.length]]));

                const result = new ApprParallelLines2().fit(points, 1024);
                const best = resultError(points, result.center,
                    result.direction, result.radius);
                let sampled = Number.MAX_VALUE;
                for (let i = 0; i < 180; ++i) {
                    sampled = Math.min(sampled,
                        modelError(points, (Math.PI * i) / 180));
                }
                // The grid is coarse, so it can only be slightly better than
                // the true minimum; a wrong direction loses by far more.
                expect(best).toBeLessThanOrEqual(
                    sampled * (1 + 1e-3) + 1e-9);
            });
    });

    it('does not modify the input samples', () => {
        check(configArb, ([angle, c, r, ss]) => {
            const points = twoLines(angle, c, r, ss);
            const before = points.map(p => [...p.values]);
            new ApprParallelLines2().fit(points, 1024);
            expect(points.map(p => [...p.values])).toEqual(before);
        });
    });

    it('rejects roots with sigma^2 > 1 for horizontal line pairs', () => {
        // gamma^2 + sigma^2 = 1, so a root of the squared polynomial h with
        // sigma^2 > 1 is spurious. Upstream keeps it, and for samples on two
        // horizontal lines it returns the non-unit direction (0, 1.7320509)
        // with a NaN radius. The port filters those roots out.
        check(fc.tuple(finite(-6, 6), positive(4, 0.5), positionsArb),
            ([offset, r, xs]) => {
                const points: Vector[] = [];
                for (const sign of [-1, 1]) {
                    for (const x of xs) {
                        points.push(Vector.fromArray(
                            [x, offset + sign * r]));
                    }
                }
                const result = new ApprParallelLines2().fit(points, 1024);
                expectClose(dot(result.direction, result.direction), 1,
                    1e-9, 1e-9);
                expect(Number.isFinite(result.radius)).toBe(true);
                expectClose(result.radius, r, 1e-6, 1e-6);
                expectClose(Math.abs(result.direction.get(0)), 1, 1e-6, 1e-6);
            });
    });
});
