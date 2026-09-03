import { describe, it, expect } from 'vitest';
import { ApprGaussian2 } from '../src/ApprGaussian2.js';
import { Vector, dot } from '../src/Vector.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// The covariance matrix (normalized by the number of points, as upstream
// computes it) of the point set, computed independently of the fitter.
function covariance(points: readonly Vector[]): {
    mean: Vector, c00: number, c01: number, c11: number
} {
    const n = points.length;
    const mean = v2(0, 0);
    for (const p of points) {
        mean.values[0] += p.values[0];
        mean.values[1] += p.values[1];
    }
    mean.values[0] /= n;
    mean.values[1] /= n;
    let c00 = 0, c01 = 0, c11 = 0;
    for (const p of points) {
        const d0 = p.values[0] - mean.values[0];
        const d1 = p.values[1] - mean.values[1];
        c00 += d0 * d0;
        c01 += d0 * d1;
        c11 += d1 * d1;
    }
    return { mean: mean, c00: c00 / n, c01: c01 / n, c11: c11 / n };
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprGaussian2', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprGaussian2();
        const box = fitter.getParameters();
        expect(box.center.values).toEqual([0, 0]);
        expect(box.axis[0].values).toEqual([0, 0]);
        expect(box.axis[1].values).toEqual([0, 0]);
        expect(box.extent.values).toEqual([0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
    });

    it('recovers the mean and the axis-aligned covariance', () => {
        // Samples with variance 4 along x and variance 1 along y, centered
        // at (10,-5).
        const points: Vector[] = [];
        for (const dx of [-2, 2]) {
            for (const dy of [-1, 1]) {
                points.push(v2(10 + dx, -5 + dy));
            }
        }
        const fitter = new ApprGaussian2();
        expect(fitter.fit(points)).toBe(true);
        const box = fitter.getParameters();
        expect(box.center.values[0]).toBeCloseTo(10, 12);
        expect(box.center.values[1]).toBeCloseTo(-5, 12);

        // The extents are the eigenvalues in increasing order, so the
        // smaller variance (y) comes first.
        expect(box.extent.values[0]).toBeCloseTo(1, 12);
        expect(box.extent.values[1]).toBeCloseTo(4, 12);
        expect(Math.abs(box.axis[0].values[1])).toBeCloseTo(1, 12);
        expect(Math.abs(box.axis[1].values[0])).toBeCloseTo(1, 12);
    });

    it('recovers the axes of a rotated distribution', () => {
        const angle = 0.4;
        const c = Math.cos(angle), s = Math.sin(angle);
        const points: Vector[] = [];
        for (const u of [-3, -1, 1, 3]) {
            for (const w of [-0.5, 0.5]) {
                points.push(v2(1 + u * c - w * s, 2 + u * s + w * c));
            }
        }
        const fitter = new ApprGaussian2();
        expect(fitter.fit(points)).toBe(true);
        const box = fitter.getParameters();
        expect(box.center.values[0]).toBeCloseTo(1, 12);
        expect(box.center.values[1]).toBeCloseTo(2, 12);

        // Variance along the (c,s) axis is (9+9+1+1)/4 = 5, along the
        // perpendicular axis it is 0.25.
        expect(box.extent.values[0]).toBeCloseTo(0.25, 12);
        expect(box.extent.values[1]).toBeCloseTo(5, 12);

        // axis[1] is parallel to (c,s).
        expect(Math.abs(dot(box.axis[1], v2(c, s)))).toBeCloseTo(1, 12);
        expect(Math.abs(dot(box.axis[0], v2(-s, c)))).toBeCloseTo(1, 12);
    });

    it('satisfies the eigenvalue identity for random point sets', () => {
        const random = makeRandom(2024);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 30; ++i) {
                points.push(v2(10 * random() - 5, 4 * random() - 2));
            }
            const fitter = new ApprGaussian2();
            expect(fitter.fit(points)).toBe(true);
            const box = fitter.getParameters();
            const cov = covariance(points);

            expect(box.center.values[0]).toBeCloseTo(cov.mean.values[0], 10);
            expect(box.center.values[1]).toBeCloseTo(cov.mean.values[1], 10);

            // Increasing order of the eigenvalues.
            expect(box.extent.values[0]).toBeLessThanOrEqual(box.extent.values[1]);

            // Trace and determinant identities.
            expect(box.extent.values[0] + box.extent.values[1])
                .toBeCloseTo(cov.c00 + cov.c11, 10);
            expect(box.extent.values[0] * box.extent.values[1])
                .toBeCloseTo(cov.c00 * cov.c11 - cov.c01 * cov.c01, 10);

            // Orthonormality of the axes.
            for (let i = 0; i < 2; ++i) {
                expect(dot(box.axis[i], box.axis[i])).toBeCloseTo(1, 12);
            }
            expect(dot(box.axis[0], box.axis[1])).toBeCloseTo(0, 12);

            // C * axis[i] = extent[i] * axis[i].
            for (let i = 0; i < 2; ++i) {
                const a0 = box.axis[i].values[0], a1 = box.axis[i].values[1];
                expect(cov.c00 * a0 + cov.c01 * a1)
                    .toBeCloseTo(box.extent.values[i] * a0, 10);
                expect(cov.c01 * a0 + cov.c11 * a1)
                    .toBeCloseTo(box.extent.values[i] * a1, 10);
            }
        }
    });

    it('produces zero extents for coincident points', () => {
        const points = [v2(3, 4), v2(3, 4), v2(3, 4)];
        const fitter = new ApprGaussian2();
        expect(fitter.fit(points)).toBe(true);
        const box = fitter.getParameters();
        expect(box.center.values[0]).toBeCloseTo(3, 12);
        expect(box.center.values[1]).toBeCloseTo(4, 12);
        expect(box.extent.values[0]).toBeCloseTo(0, 12);
        expect(box.extent.values[1]).toBeCloseTo(0, 12);

        // No extent is positive, so the error is zero for every point.
        expect(fitter.error(v2(100, -100))).toBe(0);
    });

    it('has a zero extent for collinear points', () => {
        const points = [v2(0, 0), v2(1, 1), v2(2, 2), v2(-1, -1)];
        const fitter = new ApprGaussian2();
        expect(fitter.fit(points)).toBe(true);
        const box = fitter.getParameters();
        expect(box.extent.values[0]).toBeCloseTo(0, 12);
        expect(box.extent.values[1]).toBeGreaterThan(0);
        expect(Math.abs(dot(box.axis[1], v2(Math.SQRT1_2, Math.SQRT1_2))))
            .toBeCloseTo(1, 12);
    });

    it('computes the Mahalanobis-style error', () => {
        const points: Vector[] = [];
        for (const dx of [-2, 2]) {
            for (const dy of [-1, 1]) {
                points.push(v2(dx, dy));
            }
        }
        const fitter = new ApprGaussian2();
        fitter.fit(points);
        // extent = (1,4) with axis[0] = +-(0,1) and axis[1] = +-(1,0), so
        // the error at (2,1) is (1/1)^2 + (2/4)^2.
        expect(fitter.error(v2(2, 1))).toBeCloseTo(1 + 0.25, 12);
        expect(fitter.error(v2(0, 0))).toBeCloseTo(0, 12);
    });

    it('deep-copies the parameters', () => {
        const source = new ApprGaussian2();
        source.fit([v2(-1, 0), v2(1, 0), v2(0, -2), v2(0, 2)]);
        const target = new ApprGaussian2();
        target.copyParameters(source);
        expect(target.getParameters().extent.values)
            .toEqual(source.getParameters().extent.values);

        // Mutating the source parameters must not affect the copy.
        const before = target.getParameters().center.values[0];
        source.getParameters().center.values[0] = 99;
        expect(target.getParameters().center.values[0]).toBe(before);
    });

    it('fits an indexed subset of the observations', () => {
        const points = [v2(-2, 0), v2(2, 0), v2(1000, 1000)];
        const fitter = new ApprGaussian2();
        expect(fitter.fit(points, [0, 1])).toBe(true);
        expect(fitter.getParameters().center.values[0]).toBeCloseTo(0, 12);
        expect(fitter.getParameters().extent.values[1]).toBeCloseTo(4, 12);
    });
});
