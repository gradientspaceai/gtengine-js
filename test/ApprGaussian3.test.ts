import { describe, it, expect } from 'vitest';
import { ApprGaussian3 } from '../src/ApprGaussian3';
import { Vector, dot } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The covariance matrix (normalized by the number of points, as upstream
// computes it), computed independently of the fitter.
function covariance(points: readonly Vector[]): { mean: Vector, c: number[][] } {
    const n = points.length;
    const mean = v3(0, 0, 0);
    for (const p of points) {
        for (let d = 0; d < 3; ++d) {
            mean.values[d] += p.values[d];
        }
    }
    for (let d = 0; d < 3; ++d) {
        mean.values[d] /= n;
    }
    const c = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const p of points) {
        const d = [
            p.values[0] - mean.values[0],
            p.values[1] - mean.values[1],
            p.values[2] - mean.values[2]
        ];
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                c[i][j] += d[i] * d[j];
            }
        }
    }
    for (let i = 0; i < 3; ++i) {
        for (let j = 0; j < 3; ++j) {
            c[i][j] /= n;
        }
    }
    return { mean: mean, c: c };
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprGaussian3', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprGaussian3();
        const box = fitter.getParameters();
        expect(box.center.values).toEqual([0, 0, 0]);
        for (let i = 0; i < 3; ++i) {
            expect(box.axis[i].values).toEqual([0, 0, 0]);
        }
        expect(box.extent.values).toEqual([0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
    });

    it('recovers the mean and the axis-aligned covariance', () => {
        // Variances 9, 4 and 1 along x, y and z, centered at (1,2,3).
        const points: Vector[] = [];
        for (const dx of [-3, 3]) {
            for (const dy of [-2, 2]) {
                for (const dz of [-1, 1]) {
                    points.push(v3(1 + dx, 2 + dy, 3 + dz));
                }
            }
        }
        const fitter = new ApprGaussian3();
        expect(fitter.fit(points)).toBe(true);
        const box = fitter.getParameters();
        expect(box.center.values[0]).toBeCloseTo(1, 12);
        expect(box.center.values[1]).toBeCloseTo(2, 12);
        expect(box.center.values[2]).toBeCloseTo(3, 12);

        // Increasing order of eigenvalues: z, y, x.
        expect(box.extent.values[0]).toBeCloseTo(1, 10);
        expect(box.extent.values[1]).toBeCloseTo(4, 10);
        expect(box.extent.values[2]).toBeCloseTo(9, 10);
        expect(Math.abs(box.axis[0].values[2])).toBeCloseTo(1, 10);
        expect(Math.abs(box.axis[1].values[1])).toBeCloseTo(1, 10);
        expect(Math.abs(box.axis[2].values[0])).toBeCloseTo(1, 10);
    });

    it('recovers the axes of a rotated distribution', () => {
        // Rotate the frame about the z-axis by 0.6 radians.
        const angle = 0.6;
        const c = Math.cos(angle), s = Math.sin(angle);
        const e0 = v3(c, s, 0), e1 = v3(-s, c, 0), e2 = v3(0, 0, 1);
        const points: Vector[] = [];
        for (const a of [-4, 4]) {
            for (const b of [-2, 2]) {
                for (const d of [-1, 1]) {
                    points.push(v3(
                        a * e0.values[0] + b * e1.values[0] + d * e2.values[0],
                        a * e0.values[1] + b * e1.values[1] + d * e2.values[1],
                        a * e0.values[2] + b * e1.values[2] + d * e2.values[2]));
                }
            }
        }
        const fitter = new ApprGaussian3();
        expect(fitter.fit(points)).toBe(true);
        const box = fitter.getParameters();
        expect(box.extent.values[0]).toBeCloseTo(1, 10);
        expect(box.extent.values[1]).toBeCloseTo(4, 10);
        expect(box.extent.values[2]).toBeCloseTo(16, 10);
        expect(Math.abs(dot(box.axis[0], e2))).toBeCloseTo(1, 10);
        expect(Math.abs(dot(box.axis[1], e1))).toBeCloseTo(1, 10);
        expect(Math.abs(dot(box.axis[2], e0))).toBeCloseTo(1, 10);
    });

    it('satisfies the eigenvalue identities for random point sets', () => {
        const random = makeRandom(31337);
        for (let trial = 0; trial < 15; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 40; ++i) {
                points.push(v3(6 * random() - 3, 4 * random() - 2,
                    2 * random() - 1));
            }
            const fitter = new ApprGaussian3();
            expect(fitter.fit(points)).toBe(true);
            const box = fitter.getParameters();
            const cov = covariance(points);

            for (let d = 0; d < 3; ++d) {
                expect(box.center.values[d]).toBeCloseTo(cov.mean.values[d], 10);
            }

            // Increasing order of the eigenvalues.
            expect(box.extent.values[0]).toBeLessThanOrEqual(box.extent.values[1]);
            expect(box.extent.values[1]).toBeLessThanOrEqual(box.extent.values[2]);

            // The trace identity.
            const trace = cov.c[0][0] + cov.c[1][1] + cov.c[2][2];
            expect(box.extent.values[0] + box.extent.values[1]
                + box.extent.values[2]).toBeCloseTo(trace, 8);

            // Orthonormality of the axes.
            for (let i = 0; i < 3; ++i) {
                expect(dot(box.axis[i], box.axis[i])).toBeCloseTo(1, 10);
                for (let j = i + 1; j < 3; ++j) {
                    expect(dot(box.axis[i], box.axis[j])).toBeCloseTo(0, 10);
                }
            }

            // C * axis[i] = extent[i] * axis[i].
            for (let i = 0; i < 3; ++i) {
                for (let r = 0; r < 3; ++r) {
                    let sum = 0;
                    for (let k = 0; k < 3; ++k) {
                        sum += cov.c[r][k] * box.axis[i].values[k];
                    }
                    expect(sum).toBeCloseTo(
                        box.extent.values[i] * box.axis[i].values[r], 8);
                }
            }
        }
    });

    it('produces zero extents for coincident points', () => {
        const points = [v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3)];
        const fitter = new ApprGaussian3();
        expect(fitter.fit(points)).toBe(true);
        const box = fitter.getParameters();
        expect(box.center.values[0]).toBeCloseTo(1, 12);
        for (let i = 0; i < 3; ++i) {
            expect(box.extent.values[i]).toBeCloseTo(0, 12);
        }
        expect(fitter.error(v3(50, 50, 50))).toBe(0);
    });

    it('has two zero extents for collinear points', () => {
        const points = [v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2), v3(3, 3, 3)];
        const fitter = new ApprGaussian3();
        expect(fitter.fit(points)).toBe(true);
        const box = fitter.getParameters();
        expect(box.extent.values[0]).toBeCloseTo(0, 10);
        expect(box.extent.values[1]).toBeCloseTo(0, 10);
        expect(box.extent.values[2]).toBeGreaterThan(0);
        const unit = 1 / Math.sqrt(3);
        expect(Math.abs(dot(box.axis[2], v3(unit, unit, unit))))
            .toBeCloseTo(1, 10);
    });

    it('computes the Mahalanobis-style error', () => {
        const points: Vector[] = [];
        for (const dx of [-3, 3]) {
            for (const dy of [-2, 2]) {
                for (const dz of [-1, 1]) {
                    points.push(v3(dx, dy, dz));
                }
            }
        }
        const fitter = new ApprGaussian3();
        fitter.fit(points);
        // extent = (1,4,9) along +-z, +-y, +-x, so the error at (9,4,1) is
        // (1/1)^2 + (4/4)^2 + (9/9)^2 = 3.
        expect(fitter.error(v3(9, 4, 1))).toBeCloseTo(3, 8);
        expect(fitter.error(v3(0, 0, 0))).toBeCloseTo(0, 12);
    });

    it('deep-copies the parameters', () => {
        const source = new ApprGaussian3();
        source.fit([v3(-1, 0, 0), v3(1, 0, 0), v3(0, -2, 0), v3(0, 2, 0),
            v3(0, 0, -3), v3(0, 0, 3)]);
        const target = new ApprGaussian3();
        target.copyParameters(source);
        expect(target.getParameters().extent.values)
            .toEqual(source.getParameters().extent.values);
        const before = target.getParameters().axis[0].values[0];
        source.getParameters().axis[0].values[0] = 42;
        expect(target.getParameters().axis[0].values[0]).toBe(before);
    });
});
