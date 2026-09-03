import { describe, it, expect } from 'vitest';
import { ApprOrthogonalLine3 } from '../src/ApprOrthogonalLine3.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { check, expectClose, fc, finite, rotationFrame, vector, wellScaledVector } from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprOrthogonalLine3', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprOrthogonalLine3();
        const line = fitter.getParameters();
        expect(line.origin.values).toEqual([0, 0, 0]);
        expect(line.direction.values).toEqual([0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
    });

    it('recovers a line from points that lie exactly on it', () => {
        const unit = 1 / Math.sqrt(3);
        const d = v3(unit, unit, unit);
        const points: Vector[] = [];
        for (const t of [-3, -1, 0, 1, 3]) {
            points.push(v3(1 + t * d.values[0], -2 + t * d.values[1],
                4 + t * d.values[2]));
        }
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(true);
        const line = fitter.getParameters();
        expect(line.origin.values[0]).toBeCloseTo(1, 12);
        expect(line.origin.values[1]).toBeCloseTo(-2, 12);
        expect(line.origin.values[2]).toBeCloseTo(4, 12);
        expect(Math.abs(dot(line.direction, d))).toBeCloseTo(1, 10);
        expect(dot(line.direction, line.direction)).toBeCloseTo(1, 12);

        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 14);
        }
    });

    it('recovers an axis-aligned line', () => {
        const points = [v3(0, 3, -1), v3(1, 3, -1), v3(5, 3, -1), v3(-4, 3, -1)];
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(true);
        const line = fitter.getParameters();
        expect(Math.abs(line.direction.values[0])).toBeCloseTo(1, 10);
        expect(line.origin.values[1]).toBeCloseTo(3, 12);
        expect(line.origin.values[2]).toBeCloseTo(-1, 12);
    });

    it('computes the squared orthogonal distance as the error', () => {
        const points = [v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0), v3(3, 0, 0)];
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(true);
        // The fitted line is the x-axis through (1.5,0,0).
        expect(fitter.error(v3(0, 3, 4))).toBeCloseTo(25, 8);
        expect(fitter.error(v3(-100, 0, 0))).toBeCloseTo(0, 8);
    });

    it('is not unique for an isotropic point set', () => {
        // The vertices of a regular octahedron have covariance a multiple of
        // the identity, so the maximum eigenvalue has multiplicity 3.
        const points = [v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
            v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)];
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(false);
        expect(dot(fitter.getParameters().direction,
            fitter.getParameters().direction)).toBeCloseTo(1, 10);
    });

    it('is not unique for coincident points', () => {
        const points = [v3(2, 2, 2), v3(2, 2, 2), v3(2, 2, 2)];
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().origin.values[0]).toBeCloseTo(2, 12);
    });

    it('is not unique for a planar isotropic point set', () => {
        // Points on a circle in the z = 0 plane: the two largest eigenvalues
        // are equal.
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            const angle = 2 * Math.PI * i / 8;
            points.push(v3(Math.cos(angle), Math.sin(angle), 0));
        }
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(false);
    });

    it('minimizes the sum of squared orthogonal distances', () => {
        const random = makeRandom(5150);
        for (let trial = 0; trial < 10; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 25; ++i) {
                const t = 8 * random() - 4;
                points.push(v3(
                    t + 0.15 * (2 * random() - 1),
                    2 * t + 1 + 0.15 * (2 * random() - 1),
                    -t + 0.15 * (2 * random() - 1)));
            }
            const fitter = new ApprOrthogonalLine3();
            expect(fitter.fit(points)).toBe(true);
            const line = fitter.getParameters();

            let best = 0;
            for (const p of points) {
                best += fitter.error(p);
            }

            // Perturbing the direction must not reduce the objective.
            for (const eps of [-0.05, -0.01, 0.01, 0.05]) {
                for (let axis = 0; axis < 3; ++axis) {
                    const u = line.direction.clone();
                    u.values[axis] += eps;
                    const len = Math.sqrt(dot(u, u));
                    for (let d = 0; d < 3; ++d) {
                        u.values[d] /= len;
                    }
                    let sum = 0;
                    for (const p of points) {
                        let sqrlen = 0, proj = 0;
                        for (let d = 0; d < 3; ++d) {
                            const diff = p.values[d] - line.origin.values[d];
                            sqrlen += diff * diff;
                            proj += diff * u.values[d];
                        }
                        sum += sqrlen - proj * proj;
                    }
                    expect(sum).toBeGreaterThanOrEqual(best - 1e-10);
                }
            }
        }
    });

    it('deep-copies the parameters', () => {
        const source = new ApprOrthogonalLine3();
        source.fit([v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2)]);
        const target = new ApprOrthogonalLine3();
        target.copyParameters(source);
        expect(target.getParameters().origin.values[0]).toBeCloseTo(1, 12);
        source.getParameters().origin.values[0] = 33;
        expect(target.getParameters().origin.values[0]).toBeCloseTo(1, 12);
    });
});

describe('ApprOrthogonalLine3 verification', () => {
    const pointsArb = fc.array(wellScaledVector(3, -10, 10),
        { minLength: 2, maxLength: 12 });

    it('reports the mean as the origin and a unit-length direction', () => {
        check(pointsArb, points => {
            const fitter = new ApprOrthogonalLine3();
            fitter.fit(points);
            const line = fitter.getParameters();

            const mean = [0, 0, 0];
            for (const p of points) {
                for (let d = 0; d < 3; ++d) { mean[d] += p.get(d); }
            }
            const invN = 1 / points.length;
            for (let d = 0; d < 3; ++d) {
                expectClose(line.origin.get(d), mean[d] * invN, 1e-9, 1e-9);
            }
            expectClose(dot(line.direction, line.direction), 1, 1e-12, 1e-12);
        });
    });

    it('the model error is the squared distance to the fitted line', () => {
        check(fc.tuple(pointsArb, vector(3, -20, 20)), ([points, q]) => {
            const fitter = new ApprOrthogonalLine3();
            fitter.fit(points);
            const line = fitter.getParameters();

            const diff = sub(q, line.origin);
            const t = dot(diff, line.direction);
            const d = sub(q, add(line.origin, mul(t, line.direction)));
            const expected = dot(d, d);

            // error() is |diff|^2 - dot^2, which cancels when q is far along
            // the line, so the comparison is relative to |diff|^2.
            const scale = dot(diff, diff);
            expect(Math.abs(fitter.error(q) - expected))
                .toBeLessThanOrEqual(1e-9 + 1e-9 * scale);
        });
    });

    it('recovers a line its samples lie on', () => {
        check(fc.tuple(vector(3, -8, 8), rotationFrame(3),
            fc.array(finite(-10, 10), { minLength: 2, maxLength: 10 })
                .filter(ts => Math.max(...ts) - Math.min(...ts) > 0.5)),
            ([origin, frame, ts]) => {
                const dir = frame[0];
                const points = ts.map(t => add(origin, mul(t, dir)));
                const fitter = new ApprOrthogonalLine3();
                expect(fitter.fit(points)).toBe(true);
                const line = fitter.getParameters();
                expectClose(Math.abs(dot(line.direction, dir)), 1, 1e-7, 1e-7);
                for (const p of points) {
                    expectClose(fitter.error(p), 0, 1e-12, 0);
                }
            });
    });

    it('the total fit error is invariant under a rigid motion of the data',
        () => {
            // sum_i error(P_i) is the least-squares objective, a spectral
            // quantity of the covariance and therefore a rigid invariant,
            // unlike the fitted direction when the two largest eigenvalues
            // are close.
            check(fc.tuple(pointsArb, rotationFrame(3), vector(3, -10, 10)),
                ([points, frame, t]) => {
                    const move = (p: Vector): Vector => add(t,
                        add(mul(p.get(0), frame[0]),
                            add(mul(p.get(1), frame[1]),
                                mul(p.get(2), frame[2]))));

                    const f0 = new ApprOrthogonalLine3();
                    const f1 = new ApprOrthogonalLine3();
                    f0.fit(points);
                    f1.fit(points.map(move));

                    let e0 = 0, e1 = 0, scale = 0;
                    for (const p of points) {
                        e0 += f0.error(p);
                        e1 += f1.error(move(p));
                        scale += dot(p, p);
                    }
                    expect(Math.abs(e1 - e0))
                        .toBeLessThanOrEqual(1e-7 + 1e-7 * scale);
                });
        });

    it('reports a non-unique fit for isotropic and coincident data', () => {
        // The six samples (+/-r,0,0), (0,+/-r,0), (0,0,+/-r) have covariance
        // (r^2/3) * I, so the maximum eigenvalue is not simple. Integer
        // coordinates keep the mean exactly zero, which is what makes the
        // three eigenvalues bit-identical.
        check(fc.integer({ min: 1, max: 8 }), r => {
            const isotropic: Vector[] = [];
            for (let d = 0; d < 3; ++d) {
                for (const s of [r, -r]) {
                    const p = new Vector(3);
                    p.set(d, s);
                    isotropic.push(p);
                }
            }
            const coincident = [Vector.fromArray([r, r, r]),
                Vector.fromArray([r, r, r])];
            for (const points of [isotropic, coincident]) {
                const fitter = new ApprOrthogonalLine3();
                expect(fitter.fit(points)).toBe(false);
                expectClose(dot(fitter.getParameters().direction,
                    fitter.getParameters().direction), 1, 1e-12, 1e-12);
            }
        });
    });

    it('leaves the origin at the origin for an empty sample set', () => {
        // Upstream multiplies by 1/numIndices here, so an empty set gives
        // (NaN,NaN,NaN); the fit still reports failure. Pinned so the
        // difference from ApprOrthogonalLine2 (which uses operator/=) stays
        // visible.
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit([])).toBe(false);
        expect(fitter.getParameters().origin.values.every(Number.isNaN))
            .toBe(true);
    });
});
