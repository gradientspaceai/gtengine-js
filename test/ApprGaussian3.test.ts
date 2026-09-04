import { describe, it, expect } from 'vitest';
import { ApprGaussian3 } from '../src/ApprGaussian3.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { ApprQuery } from '../src/ApprQuery.js';
import { check, expectClose, expectVectorClose, fc, positive, rotationFrame, vector, wellScaledVector } from './helpers/arbitraries.js';

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

describe('ApprGaussian3 verification', () => {
    // The covariance of the points about their mean, computed directly.
    const moments = (points: readonly Vector[]) => {
        const n = points.length;
        const mean = [0, 0, 0];
        for (const p of points) {
            for (let d = 0; d < 3; ++d) { mean[d] += p.get(d); }
        }
        for (let d = 0; d < 3; ++d) { mean[d] /= n; }
        const c = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (const p of points) {
            const df = [p.get(0) - mean[0], p.get(1) - mean[1],
                p.get(2) - mean[2]];
            for (let i = 0; i < 3; ++i) {
                for (let j = 0; j < 3; ++j) { c[i][j] += df[i] * df[j]; }
            }
        }
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) { c[i][j] /= n; }
        }
        return { mean: mean, c: c };
    };

    // wellScaledVector keeps every coordinate either exactly zero or above
    // 1e-3, so the covariance never acquires a subnormal off-diagonal entry.
    // Upstream's GetCosSin computes sqrt(u*u + v*v) without scaling, which
    // loses most of the mantissa when u*u + v*v is subnormal.
    const pointsArb = fc.array(wellScaledVector(3, -10, 10),
        { minLength: 2, maxLength: 12 });

    it('reports the mean, an orthonormal frame and the eigenvalues in '
        + 'increasing order', () => {
            check(pointsArb, points => {
                const fitter = new ApprGaussian3();
                expect(fitter.fitIndexed(points,
                    points.map((_, i) => i))).toBe(true);
                const box = fitter.getParameters();
                const m = moments(points);

                for (let d = 0; d < 3; ++d) {
                    expectClose(box.center.get(d), m.mean[d], 1e-9, 1e-9);
                }

                // The axes are orthonormal.
                for (let i = 0; i < 3; ++i) {
                    for (let j = 0; j < 3; ++j) {
                        expectClose(dot(box.axis[i], box.axis[j]),
                            i === j ? 1 : 0, 1e-10, 1e-10);
                    }
                }

                // The extents are the eigenvalues in increasing order.
                expect(box.extent.get(0)).toBeLessThanOrEqual(box.extent.get(1));
                expect(box.extent.get(1)).toBeLessThanOrEqual(box.extent.get(2));

                // Each (extent, axis) pair is an eigenpair of the covariance.
                for (let i = 0; i < 3; ++i) {
                    const a = box.axis[i];
                    const e = box.extent.get(i);
                    for (let r = 0; r < 3; ++r) {
                        const lhs = m.c[r][0] * a.get(0) + m.c[r][1] * a.get(1)
                            + m.c[r][2] * a.get(2);
                        expectClose(lhs, e * a.get(r), 1e-8, 1e-8);
                    }
                }
            });
        });

    it('recovers the principal frame of an axis-symmetric sample set', () => {
        // The six samples C +/- a*U0, C +/- b*U1, C +/- c*U2 have mean C and
        // covariance diag(a^2/3, b^2/3, c^2/3) in the {U0,U1,U2} basis.
        // wellScaledVector: a center component near 1e-147 puts covariance
        // off-diagonals in the underflow region of SymmetricEigensolver3x3
        // (#379) and the eigenvalues come back wrong; that is upstream's
        // conditioning, not this fitter's.
        check(fc.tuple(wellScaledVector(3, -8, 8), positive(6, 1), positive(6, 1),
            positive(6, 1), rotationFrame(3))
            .filter(([, a, b, c]) => Math.abs(a - b) > 0.25
                && Math.abs(b - c) > 0.25 && Math.abs(a - c) > 0.25),
            ([center, a, b, c, frame]) => {
                const scales = [a, b, c];
                const points: Vector[] = [];
                for (let i = 0; i < 3; ++i) {
                    for (const s of [scales[i], -scales[i]]) {
                        points.push(Vector.fromArray([
                            center.get(0) + s * frame[i].get(0),
                            center.get(1) + s * frame[i].get(1),
                            center.get(2) + s * frame[i].get(2)]));
                    }
                }

                const fitter = new ApprGaussian3();
                expect(fitter.fit(points)).toBe(true);
                const box = fitter.getParameters();
                expectVectorClose(box.center, center, 1e-9, 1e-9);

                const evals = scales.map(s => (s * s) / 3);
                const order = [0, 1, 2].sort((i, j) => evals[i] - evals[j]);
                for (let i = 0; i < 3; ++i) {
                    expectClose(box.extent.get(i), evals[order[i]], 1e-7, 1e-7);
                    expectClose(Math.abs(dot(box.axis[i], frame[order[i]])), 1,
                        1e-6, 1e-6);
                }
            });
    });

    it('the model error is the sum of squared axis ratios and vanishes at '
        + 'the center', () => {
            check(fc.tuple(pointsArb, vector(3, -20, 20)), ([points, p]) => {
                const fitter = new ApprGaussian3();
                expect(fitter.fit(points)).toBe(true);
                const box = fitter.getParameters();
                expect(fitter.error(box.center.clone())).toBe(0);

                const diff = sub(p, box.center);
                let expected = 0;
                for (let i = 0; i < 3; ++i) {
                    if (box.extent.get(i) > 0) {
                        const ratio = dot(diff, box.axis[i]) / box.extent.get(i);
                        expected += ratio * ratio;
                    }
                }
                expect(fitter.error(p)).toBe(expected);
            });
        });

    it('coplanar samples give a zero smallest extent along the plane normal',
        () => {
            // The samples span the plane z = 3 (the triangle (0,0), (a,0),
            // (0,b) guarantees rank 2 in the plane), so the covariance has
            // exactly one zero eigenvalue and its eigenvector is the plane
            // normal.
            check(fc.tuple(positive(8, 1), positive(8, 1),
                fc.array(wellScaledVector(2, -5, 5), { maxLength: 6 })),
                ([a, b, extra]) => {
                    const flat = [Vector.fromArray([0, 0]),
                        Vector.fromArray([a, 0]), Vector.fromArray([0, b]),
                        ...extra];
                    const points = flat.map(
                        q => Vector.fromArray([q.get(0), q.get(1), 3]));
                    const fitter = new ApprGaussian3();
                    expect(fitter.fit(points)).toBe(true);
                    const box = fitter.getParameters();
                    expectClose(box.extent.get(0), 0, 1e-11, 0);
                    expect(box.extent.get(1)).toBeGreaterThan(1e-3);
                    expectClose(Math.abs(box.axis[0].get(2)), 1, 1e-8, 1e-8);
                });
        });

    it('copyParameters deep-copies the oriented box', () => {
        check(pointsArb, points => {
            const source = new ApprGaussian3();
            expect(source.fit(points)).toBe(true);
            const target = new ApprGaussian3();
            target.copyParameters(source);

            const a = source.getParameters();
            const b = target.getParameters();
            expect(b).not.toBe(a);
            expectVectorClose(b.center, a.center, 0, 0);
            expectVectorClose(b.extent, a.extent, 0, 0);

            const centerBefore = [...b.center.values];
            a.center.set(0, 12345);
            a.axis[0].set(0, -999);
            expect([...b.center.values]).toEqual(centerBefore);
            expect(b.axis[0].get(0)).not.toBe(-999);
        });
    });

    it('a failed fit zeroes every parameter', () => {
        ApprQuery.validateIndices = true;
        try {
            const fitter = new ApprGaussian3();
            expect(fitter.fit([Vector.fromArray([1, 2, 3]),
                Vector.fromArray([4, 5, 6])])).toBe(true);
            expect(fitter.fit([Vector.fromArray([1, 2, 3])], [0, 5]))
                .toBe(false);
            const box = fitter.getParameters();
            expect(box.center.values).toEqual([0, 0, 0]);
            expect(box.extent.values).toEqual([0, 0, 0]);
            for (let i = 0; i < 3; ++i) {
                expect(box.axis[i].values).toEqual([0, 0, 0]);
            }
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });
});
