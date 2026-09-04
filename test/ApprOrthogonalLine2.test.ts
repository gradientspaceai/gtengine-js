import { describe, it, expect } from 'vitest';
import { ApprOrthogonalLine2 } from '../src/ApprOrthogonalLine2.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { check, expectClose, fc, rotationFrame, vector, wellScaled, wellScaledVector } from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprOrthogonalLine2', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprOrthogonalLine2();
        const line = fitter.getParameters();
        expect(line.origin.values).toEqual([0, 0]);
        expect(line.direction.values).toEqual([0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
    });

    it('recovers a line from points that lie exactly on it', () => {
        // The line through (1,2) with unit direction (3,4)/5.
        const d = v2(0.6, 0.8);
        const points: Vector[] = [];
        for (const t of [-2, -1, 0, 1, 2]) {
            points.push(v2(1 + t * d.values[0], 2 + t * d.values[1]));
        }
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(true);
        const line = fitter.getParameters();
        expect(line.origin.values[0]).toBeCloseTo(1, 12);
        expect(line.origin.values[1]).toBeCloseTo(2, 12);
        expect(Math.abs(dot(line.direction, d))).toBeCloseTo(1, 12);
        expect(dot(line.direction, line.direction)).toBeCloseTo(1, 12);

        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 14);
        }
    });

    it('recovers a vertical line, which the height fitter cannot fit', () => {
        const points = [v2(4, -3), v2(4, 0), v2(4, 1), v2(4, 6)];
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(true);
        const line = fitter.getParameters();
        expect(line.origin.values[0]).toBeCloseTo(4, 12);
        expect(Math.abs(line.direction.values[1])).toBeCloseTo(1, 12);
        expect(line.direction.values[0]).toBeCloseTo(0, 12);
    });

    it('computes the squared orthogonal distance as the error', () => {
        const points = [v2(0, 0), v2(1, 0), v2(2, 0), v2(3, 0)];
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(true);
        // The fitted line is the x-axis through (1.5,0).
        expect(fitter.error(v2(10, 3))).toBeCloseTo(9, 10);
        expect(fitter.error(v2(-10, 0))).toBeCloseTo(0, 10);
    });

    it('is not unique for a symmetric point set', () => {
        // The covariance matrix is a multiple of the identity, so the
        // maximum eigenvalue has multiplicity 2.
        const points = [v2(1, 0), v2(-1, 0), v2(0, 1), v2(0, -1)];
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(false);
        // The line parameters are still assigned (the fit always succeeds).
        expect(fitter.getParameters().origin.values[0]).toBeCloseTo(0, 12);
        expect(dot(fitter.getParameters().direction,
            fitter.getParameters().direction)).toBeCloseTo(1, 12);
    });

    it('is not unique for coincident points', () => {
        const points = [v2(5, 5), v2(5, 5), v2(5, 5)];
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().origin.values[0]).toBeCloseTo(5, 12);
    });

    it('minimizes the sum of squared orthogonal distances', () => {
        const random = makeRandom(9091);
        for (let trial = 0; trial < 15; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 20; ++i) {
                const t = 6 * random() - 3;
                points.push(v2(t + 0.2 * (2 * random() - 1),
                    2 * t - 1 + 0.2 * (2 * random() - 1)));
            }
            const fitter = new ApprOrthogonalLine2();
            expect(fitter.fit(points)).toBe(true);
            const line = fitter.getParameters();

            let best = 0;
            for (const p of points) {
                best += fitter.error(p);
            }

            // Rotating the fitted direction about the fitted origin must not
            // reduce the sum of squared orthogonal distances.
            for (const delta of [-0.2, -0.05, -0.01, 0.01, 0.05, 0.2]) {
                const angle = Math.atan2(line.direction.values[1],
                    line.direction.values[0]) + delta;
                const u = v2(Math.cos(angle), Math.sin(angle));
                let sum = 0;
                for (const p of points) {
                    const dx = p.values[0] - line.origin.values[0];
                    const dy = p.values[1] - line.origin.values[1];
                    const proj = dx * u.values[0] + dy * u.values[1];
                    sum += dx * dx + dy * dy - proj * proj;
                }
                expect(sum).toBeGreaterThanOrEqual(best - 1e-12);
            }
        }
    });

    it('deep-copies the parameters', () => {
        const source = new ApprOrthogonalLine2();
        source.fit([v2(0, 0), v2(1, 1), v2(2, 2)]);
        const target = new ApprOrthogonalLine2();
        target.copyParameters(source);
        expect(target.getParameters().origin.values[0]).toBeCloseTo(1, 12);
        source.getParameters().origin.values[0] = 55;
        expect(target.getParameters().origin.values[0]).toBeCloseTo(1, 12);
    });
});

describe('ApprOrthogonalLine2 verification', () => {
    const pointsArb = fc.array(wellScaledVector(2, -10, 10),
        { minLength: 2, maxLength: 12 });

    it('reports the mean as the origin and a unit-length direction', () => {
        check(pointsArb, points => {
            const fitter = new ApprOrthogonalLine2();
            fitter.fit(points);
            const line = fitter.getParameters();

            const n = points.length;
            const mean = [0, 0];
            for (const p of points) {
                mean[0] += p.get(0); mean[1] += p.get(1);
            }
            const invN = 1 / n;
            mean[0] *= invN; mean[1] *= invN;
            expectClose(line.origin.get(0), mean[0], 1e-9, 1e-9);
            expectClose(line.origin.get(1), mean[1], 1e-9, 1e-9);
            expectClose(dot(line.direction, line.direction), 1, 1e-12, 1e-12);
        });
    });

    it('the model error is the squared distance to the fitted line', () => {
        check(fc.tuple(pointsArb, vector(2, -20, 20)), ([points, q]) => {
            const fitter = new ApprOrthogonalLine2();
            fitter.fit(points);
            const line = fitter.getParameters();

            // Independent computation: project q onto the line and measure.
            const diff = sub(q, line.origin);
            const t = dot(diff, line.direction);
            const foot = add(line.origin, mul(t, line.direction));
            const d = sub(q, foot);
            const expected = dot(d, d);

            // error() is |diff|^2 - dot^2, which loses digits to
            // cancellation when q is far along the line, so the comparison
            // is relative to |diff|^2.
            const scale = dot(diff, diff);
            expect(Math.abs(fitter.error(q) - expected))
                .toBeLessThanOrEqual(1e-9 + 1e-9 * scale);
        });
    });

    it('recovers a line its samples lie on', () => {
        // wellScaled keeps every coordinate either exactly zero or above
        // 1e-3; a subnormal coordinate would leave the covariance with a
        // subnormal entry, which upstream's unscaled sqrt(u*u + v*v) cannot
        // resolve.
        check(fc.tuple(wellScaledVector(2, -8, 8), rotationFrame(2),
            fc.array(wellScaled(-10, 10), { minLength: 2, maxLength: 10 })
                .filter(ts => Math.max(...ts) - Math.min(...ts) > 0.5)),
            ([origin, frame, ts]) => {
                const dir = frame[0];
                const points = ts.map(t => add(origin, mul(t, dir)));
                const fitter = new ApprOrthogonalLine2();
                expect(fitter.fit(points)).toBe(true);
                const line = fitter.getParameters();

                // The fitted direction is the line direction up to sign.
                expectClose(Math.abs(dot(line.direction, dir)), 1, 1e-8, 1e-8);
                for (const p of points) {
                    expectClose(fitter.error(p), 0, 1e-12, 0);
                }
            });
    });

    it('the total fit error is invariant under a rigid motion of the data',
        () => {
            // sum_i error(P_i) is the least-squares objective, i.e. the
            // smallest eigenvalue of the (unnormalized) covariance. Unlike
            // the fitted direction, which is ill-conditioned when the two
            // eigenvalues are close, that spectral quantity is a rigid
            // invariant.
            check(fc.tuple(pointsArb, rotationFrame(2),
                wellScaledVector(2, -10, 10)), ([points, frame, t]) => {
                    const move = (p: Vector): Vector => add(t,
                        add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])));

                    const f0 = new ApprOrthogonalLine2();
                    const f1 = new ApprOrthogonalLine2();
                    f0.fit(points);
                    f1.fit(points.map(move));

                    let e0 = 0, e1 = 0, scale = 0;
                    for (const p of points) {
                        e0 += f0.error(p);
                        e1 += f1.error(move(p));
                        scale += dot(p, p);
                    }
                    expect(Math.abs(e1 - e0))
                        .toBeLessThanOrEqual(1e-8 + 1e-8 * scale);
                });
        });

    it('reports a non-unique fit for isotropic and coincident data', () => {
        // The four samples (+/-r,0), (0,+/-r) have covariance 2r^2 * I, so
        // the maximum eigenvalue is not simple and upstream returns false.
        // Integer coordinates keep the mean exactly (0,0), which is what
        // makes the two eigenvalues bit-identical.
        check(fc.integer({ min: 1, max: 8 }), r => {
            const isotropic = [Vector.fromArray([r, 0]),
                Vector.fromArray([-r, 0]), Vector.fromArray([0, r]),
                Vector.fromArray([0, -r])];
            const diagonal = [Vector.fromArray([r, r]),
                Vector.fromArray([-r, -r]), Vector.fromArray([r, -r]),
                Vector.fromArray([-r, r])];
            const coincident = [Vector.fromArray([r, r]),
                Vector.fromArray([r, r]), Vector.fromArray([r, r])];
            for (const points of [isotropic, diagonal, coincident]) {
                const fitter = new ApprOrthogonalLine2();
                expect(fitter.fit(points)).toBe(false);
                // Upstream zeroes the parameters only when ValidIndices
                // fails, so the fitted line is still reported.
                expectClose(dot(fitter.getParameters().direction,
                    fitter.getParameters().direction), 1, 1e-12, 1e-12);
            }
        });
    });

    it('leaves the origin at (0,0) for an empty sample set', () => {
        // Regression: upstream computes 'mean /= (Real)numIndices' and GTE's
        // Vector operator/= zeroes the vector when the divisor is zero. A
        // plain division would store (NaN,NaN) here.
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit([])).toBe(false);
        const line = fitter.getParameters();
        expect(line.origin.values).toEqual([0, 0]);
        expect(Number.isNaN(line.origin.get(0))).toBe(false);
    });
});
