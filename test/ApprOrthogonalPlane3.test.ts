import { describe, it, expect } from 'vitest';
import { ApprOrthogonalPlane3 } from '../src/ApprOrthogonalPlane3.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { check, expectClose, fc, positive, rotationFrame, vector, wellScaledVector } from './helpers/arbitraries.js';

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

describe('ApprOrthogonalPlane3', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprOrthogonalPlane3();
        const p = fitter.getParameters();
        expect(p.origin.values).toEqual([0, 0, 0]);
        expect(p.normal.values).toEqual([0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(3);
    });

    it('recovers a plane from points that lie exactly on it', () => {
        // The plane through (1,1,1) with unit normal (1,1,1)/sqrt(3), using
        // two in-plane basis vectors.
        const unit = 1 / Math.sqrt(3);
        const n = v3(unit, unit, unit);
        const e0 = v3(1, -1, 0), e1 = v3(1, 1, -2);
        const points: Vector[] = [];
        for (const a of [-2, -1, 1, 2]) {
            for (const b of [-1, 0, 1]) {
                points.push(v3(
                    1 + a * e0.values[0] + b * e1.values[0],
                    1 + a * e0.values[1] + b * e1.values[1],
                    1 + a * e0.values[2] + b * e1.values[2]));
            }
        }
        const fitter = new ApprOrthogonalPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(p.origin.values[0]).toBeCloseTo(1, 12);
        expect(p.origin.values[1]).toBeCloseTo(1, 12);
        expect(p.origin.values[2]).toBeCloseTo(1, 12);
        expect(Math.abs(dot(p.normal, n))).toBeCloseTo(1, 10);
        expect(dot(p.normal, p.normal)).toBeCloseTo(1, 12);

        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 12);
        }
    });

    it('recovers a horizontal plane', () => {
        const points = [v3(0, 0, 3), v3(1, 0, 3), v3(0, 1, 3), v3(-2, 5, 3)];
        const fitter = new ApprOrthogonalPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(Math.abs(p.normal.values[2])).toBeCloseTo(1, 12);
        expect(p.origin.values[2]).toBeCloseTo(3, 12);
    });

    it('recovers a vertical plane, which the height fitter cannot fit', () => {
        // The plane x = 2.
        const points = [v3(2, 0, 0), v3(2, 1, 0), v3(2, 0, 1), v3(2, -3, 4)];
        const fitter = new ApprOrthogonalPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(Math.abs(p.normal.values[0])).toBeCloseTo(1, 12);
        expect(p.origin.values[0]).toBeCloseTo(2, 12);
    });

    it('computes the absolute distance to the plane as the error', () => {
        const points = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0)];
        const fitter = new ApprOrthogonalPlane3();
        expect(fitter.fit(points)).toBe(true);
        expect(fitter.error(v3(10, -10, 2.5))).toBeCloseTo(2.5, 10);
        expect(fitter.error(v3(10, -10, -2.5))).toBeCloseTo(2.5, 10);
        expect(fitter.error(v3(3, 4, 0))).toBeCloseTo(0, 10);
    });

    it('returns a normal orthogonal to a collinear point set', () => {
        // Two eigenvalues of the covariance matrix are zero in exact
        // arithmetic, so any unit normal orthogonal to the line minimizes the
        // objective. (Round-off makes the reported uniqueness flag
        // unreliable here, so only the geometry is checked.)
        const points = [v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2), v3(-1, -1, -1)];
        const fitter = new ApprOrthogonalPlane3();
        fitter.fit(points);
        const p = fitter.getParameters();
        expect(dot(p.normal, p.normal)).toBeCloseTo(1, 10);
        const unit = 1 / Math.sqrt(3);
        expect(dot(p.normal, v3(unit, unit, unit))).toBeCloseTo(0, 10);
        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 10);
        }
    });

    it('is not unique for coincident points', () => {
        const points = [v3(7, 8, 9), v3(7, 8, 9), v3(7, 8, 9)];
        const fitter = new ApprOrthogonalPlane3();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().origin.values[0]).toBeCloseTo(7, 12);
    });

    it('minimizes the sum of squared distances to the plane', () => {
        const random = makeRandom(60613);
        for (let trial = 0; trial < 10; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 30; ++i) {
                const x = 4 * random() - 2;
                const y = 4 * random() - 2;
                points.push(v3(x, y, 0.5 * x - 0.25 * y
                    + 0.05 * (2 * random() - 1)));
            }
            const fitter = new ApprOrthogonalPlane3();
            expect(fitter.fit(points)).toBe(true);
            const p = fitter.getParameters();

            let best = 0;
            for (const q of points) {
                const e = fitter.error(q);
                best += e * e;
            }

            // Perturbing the normal must not reduce the objective.
            for (const eps of [-0.05, -0.01, 0.01, 0.05]) {
                for (let axis = 0; axis < 3; ++axis) {
                    const u = p.normal.clone();
                    u.values[axis] += eps;
                    const len = Math.sqrt(dot(u, u));
                    for (let d = 0; d < 3; ++d) {
                        u.values[d] /= len;
                    }
                    let sum = 0;
                    for (const q of points) {
                        let proj = 0;
                        for (let d = 0; d < 3; ++d) {
                            proj += (q.values[d] - p.origin.values[d]) * u.values[d];
                        }
                        sum += proj * proj;
                    }
                    expect(sum).toBeGreaterThanOrEqual(best - 1e-10);
                }
            }
        }
    });

    it('deep-copies the parameters', () => {
        const source = new ApprOrthogonalPlane3();
        source.fit([v3(0, 0, 5), v3(1, 0, 5), v3(0, 1, 5)]);
        const target = new ApprOrthogonalPlane3();
        target.copyParameters(source);
        expect(target.getParameters().origin.values[2]).toBeCloseTo(5, 12);
        source.getParameters().origin.values[2] = 0;
        expect(target.getParameters().origin.values[2]).toBeCloseTo(5, 12);
    });
});

describe('ApprOrthogonalPlane3 verification', () => {
    const pointsArb = fc.array(wellScaledVector(3, -10, 10),
        { minLength: 3, maxLength: 12 });

    it('reports the mean as the origin and a unit-length normal', () => {
        check(pointsArb, points => {
            const fitter = new ApprOrthogonalPlane3();
            fitter.fit(points);
            const p = fitter.getParameters();

            const mean = [0, 0, 0];
            for (const q of points) {
                for (let d = 0; d < 3; ++d) { mean[d] += q.get(d); }
            }
            const invN = 1 / points.length;
            for (let d = 0; d < 3; ++d) {
                expectClose(p.origin.get(d), mean[d] * invN, 1e-9, 1e-9);
            }
            expectClose(dot(p.normal, p.normal), 1, 1e-12, 1e-12);
        });
    });

    it('the model error is the distance to the fitted plane', () => {
        check(fc.tuple(pointsArb, vector(3, -20, 20)), ([points, q]) => {
            const fitter = new ApprOrthogonalPlane3();
            fitter.fit(points);
            const p = fitter.getParameters();

            // Independent computation: drop q onto the plane and measure.
            const diff = sub(q, p.origin);
            const signed = dot(diff, p.normal);
            const foot = sub(q, mul(signed, p.normal));
            const d = sub(q, foot);
            expectClose(fitter.error(q), Math.sqrt(dot(d, d)), 1e-9, 1e-9);
            expect(fitter.error(q)).toBeGreaterThanOrEqual(0);
        });
    });

    it('recovers a plane its samples lie on', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), rotationFrame(3),
            fc.tuple(positive(8, 1), positive(8, 1),
                fc.array(wellScaledVector(2, -5, 5), { maxLength: 6 })),
            ), ([origin, frame, [a, b, extra]]) => {
                const planar = [[0, 0], [a, 0], [0, b],
                    ...extra.map(q => [q.get(0), q.get(1)])];
                const points = planar.map(([u, v]) => add(origin,
                    add(mul(u, frame[0]), mul(v, frame[1]))));

                const fitter = new ApprOrthogonalPlane3();
                expect(fitter.fit(points)).toBe(true);
                const p = fitter.getParameters();

                // The fitted normal is the plane normal up to sign.
                expectClose(Math.abs(dot(p.normal, frame[2])), 1, 1e-7, 1e-7);
                for (const q of points) {
                    expectClose(fitter.error(q), 0, 1e-8, 0);
                }
            });
    });

    it('the total squared fit error is invariant under a rigid motion',
        () => {
            // sum_i error(P_i)^2 is the least-squares objective, a spectral
            // quantity of the covariance and therefore a rigid invariant,
            // unlike the normal itself when the two smallest eigenvalues are
            // close.
            check(fc.tuple(pointsArb, rotationFrame(3),
                wellScaledVector(3, -10, 10)), ([points, frame, t]) => {
                    const move = (p: Vector): Vector => add(t,
                        add(mul(p.get(0), frame[0]),
                            add(mul(p.get(1), frame[1]),
                                mul(p.get(2), frame[2]))));

                    const f0 = new ApprOrthogonalPlane3();
                    const f1 = new ApprOrthogonalPlane3();
                    f0.fit(points);
                    f1.fit(points.map(move));

                    let e0 = 0, e1 = 0, scale = 0;
                    for (const p of points) {
                        e0 += f0.error(p) ** 2;
                        e1 += f1.error(move(p)) ** 2;
                        scale += dot(p, p);
                    }
                    expect(Math.abs(e1 - e0))
                        .toBeLessThanOrEqual(1e-7 + 1e-7 * scale);
                });
        });

    it('reports a non-unique fit for isotropic and coincident data', () => {
        // The six samples (+/-r,0,0), (0,+/-r,0), (0,0,+/-r) have covariance
        // (r^2/3) * I, so the minimum eigenvalue is not simple.
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
                Vector.fromArray([r, r, r]), Vector.fromArray([r, r, r])];
            for (const points of [isotropic, coincident]) {
                const fitter = new ApprOrthogonalPlane3();
                expect(fitter.fit(points)).toBe(false);
                expectClose(dot(fitter.getParameters().normal,
                    fitter.getParameters().normal), 1, 1e-12, 1e-12);
            }
        });
    });

    it('copyParameters deep-copies the parameters', () => {
        check(pointsArb, points => {
            const source = new ApprOrthogonalPlane3();
            source.fit(points);
            const target = new ApprOrthogonalPlane3();
            target.copyParameters(source);

            const s = source.getParameters();
            const t = target.getParameters();
            expect(t.origin).not.toBe(s.origin);
            expect(t.normal).not.toBe(s.normal);
            expect([...t.origin.values]).toEqual([...s.origin.values]);
            expect([...t.normal.values]).toEqual([...s.normal.values]);

            s.origin.set(0, 4321);
            expect(t.origin.get(0)).not.toBe(4321);
        });
    });
});
