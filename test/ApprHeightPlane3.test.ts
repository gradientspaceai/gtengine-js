import { describe, it, expect } from 'vitest';
import { ApprHeightPlane3 } from '../src/ApprHeightPlane3.js';
import { Vector } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, finite, positive, vector, wellScaledVector } from './helpers/arbitraries.js';

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

describe('ApprHeightPlane3', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprHeightPlane3();
        const p = fitter.getParameters();
        expect(p.average.values).toEqual([0, 0, 0]);
        expect(p.coefficients.values).toEqual([0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(3);
    });

    it('recovers a plane from points that lie exactly on it', () => {
        // z = 3x - 2y + 5 sampled on a 4x4 grid.
        const points: Vector[] = [];
        for (let i = 0; i < 4; ++i) {
            for (let j = 0; j < 4; ++j) {
                const x = i - 1.5, y = j - 1.5;
                points.push(v3(x, y, 3 * x - 2 * y + 5));
            }
        }
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(p.average.values[0]).toBeCloseTo(0, 12);
        expect(p.average.values[1]).toBeCloseTo(0, 12);
        expect(p.average.values[2]).toBeCloseTo(5, 12);
        expect(p.coefficients.values[0]).toBeCloseTo(3, 10);
        expect(p.coefficients.values[1]).toBeCloseTo(-2, 10);
        expect(p.coefficients.values[2]).toBe(-1);

        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 18);
        }
    });

    it('recovers a horizontal plane', () => {
        const points = [v3(0, 0, 4), v3(1, 0, 4), v3(0, 1, 4), v3(3, -2, 4)];
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(p.coefficients.values[0]).toBeCloseTo(0, 10);
        expect(p.coefficients.values[1]).toBeCloseTo(0, 10);
        expect(p.average.values[2]).toBeCloseTo(4, 12);
    });

    it('is exact for three noncollinear points', () => {
        const points = [v3(0, 0, 1), v3(1, 0, 4), v3(0, 1, -1)];
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(p.coefficients.values[0]).toBeCloseTo(3, 10);
        expect(p.coefficients.values[1]).toBeCloseTo(-2, 10);
        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 18);
        }
    });

    it('produces least-squares residuals that satisfy the normal equations', () => {
        const random = makeRandom(777);
        for (let trial = 0; trial < 15; ++trial) {
            const points: Vector[] = [];
            const n = 25;
            for (let i = 0; i < n; ++i) {
                const x = 4 * random() - 2;
                const y = 4 * random() - 2;
                points.push(v3(x, y, 1.5 * x + 0.25 * y - 3
                    + 0.5 * (2 * random() - 1)));
            }
            const fitter = new ApprHeightPlane3();
            expect(fitter.fit(points)).toBe(true);
            const p = fitter.getParameters();
            const a = p.coefficients.values[0], b = p.coefficients.values[1];

            // The residual r = a*dx + b*dy - dz must be orthogonal to dx and
            // to dy (the normal equations of the least-squares problem).
            let dotX = 0, dotY = 0, sumR = 0;
            for (const q of points) {
                const dx = q.values[0] - p.average.values[0];
                const dy = q.values[1] - p.average.values[1];
                const dz = q.values[2] - p.average.values[2];
                const r = a * dx + b * dy - dz;
                dotX += r * dx;
                dotY += r * dy;
                sumR += r;
            }
            expect(dotX).toBeCloseTo(0, 8);
            expect(dotY).toBeCloseTo(0, 8);
            expect(sumR).toBeCloseTo(0, 8);
        }
    });

    it('computes the squared vertical error', () => {
        const points = [v3(0, 0, 0), v3(1, 0, 3), v3(0, 1, -2)];
        const fitter = new ApprHeightPlane3();
        fitter.fit(points);
        // The fitted plane is z = 3x - 2y; at (1,1,5) the deviation is 4.
        expect(fitter.error(v3(1, 1, 5))).toBeCloseTo(16, 8);
    });

    it('fails for collinear (x,y) samples', () => {
        const points = [v3(0, 0, 1), v3(1, 1, 2), v3(2, 2, 3), v3(3, 3, 9)];
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().average.values).toEqual([0, 0, 0]);
        expect(fitter.getParameters().coefficients.values).toEqual([0, 0, 0]);
    });

    it('fails for coincident points', () => {
        const points = [v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3)];
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().coefficients.values).toEqual([0, 0, 0]);
    });

    it('deep-copies the parameters', () => {
        const source = new ApprHeightPlane3();
        source.fit([v3(0, 0, 1), v3(1, 0, 4), v3(0, 1, -1)]);
        const target = new ApprHeightPlane3();
        target.copyParameters(source);
        expect(target.getParameters().coefficients.values[0]).toBeCloseTo(3, 10);
        source.getParameters().coefficients.values[0] = 99;
        expect(target.getParameters().coefficients.values[0]).toBeCloseTo(3, 10);
    });
});

describe('ApprHeightPlane3 verification', () => {
    // Domain samples whose (x,y) covariance is nonsingular: the triangle
    // (0,0), (a,0), (0,b) plus optional extra points.
    const domainArb = fc.tuple(positive(8, 1), positive(8, 1),
        fc.array(wellScaledVector(2, -5, 5), { maxLength: 6 }))
        .map(([a, b, extra]) => [[0, 0], [a, 0], [0, b],
            ...extra.map(q => [q.get(0), q.get(1)])]);

    it('recovers the exact plane of height data', () => {
        check(fc.tuple(domainArb, finite(-5, 5), finite(-5, 5), finite(-5, 5)),
            ([domain, a, b, c]) => {
                const points = domain.map(([x, y]) =>
                    Vector.fromArray([x, y, a * x + b * y + c]));
                const fitter = new ApprHeightPlane3();
                expect(fitter.fit(points)).toBe(true);
                const p = fitter.getParameters();

                // The parameters are ((xAvr,yAvr,zAvr),(a,b,-1)).
                expect(p.coefficients.get(2)).toBe(-1);
                expectClose(p.coefficients.get(0), a, 1e-7, 1e-7);
                expectClose(p.coefficients.get(1), b, 1e-7, 1e-7);

                const n = domain.length;
                const xAvr = domain.reduce((u, q) => u + q[0], 0) / n;
                const yAvr = domain.reduce((u, q) => u + q[1], 0) / n;
                expectClose(p.average.get(0), xAvr, 1e-9, 1e-9);
                expectClose(p.average.get(1), yAvr, 1e-9, 1e-9);

                // The error of every sample vanishes.
                for (const q of points) {
                    expectClose(fitter.error(q), 0, 1e-12, 0);
                }
            });
    });

    it('satisfies the two least-squares normal equations', () => {
        check(fc.tuple(domainArb, fc.array(finite(-10, 10),
            { minLength: 10, maxLength: 10 })), ([domain, zs]) => {
                const points = domain.map(([x, y], i) =>
                    Vector.fromArray([x, y, zs[i % zs.length]]));
                const fitter = new ApprHeightPlane3();
                expect(fitter.fit(points)).toBe(true);
                const p = fitter.getParameters();
                const a = p.coefficients.get(0), b = p.coefficients.get(1);

                // d/da and d/db of sum [a*dx + b*dy - dz]^2 vanish.
                let ra = 0, rb = 0, scale = 0;
                for (const q of points) {
                    const dx = q.get(0) - p.average.get(0);
                    const dy = q.get(1) - p.average.get(1);
                    const dz = q.get(2) - p.average.get(2);
                    const r = a * dx + b * dy - dz;
                    ra += r * dx;
                    rb += r * dy;
                    scale += Math.abs(a * dx) + Math.abs(b * dy) + Math.abs(dz);
                }
                expect(Math.abs(ra)).toBeLessThanOrEqual(1e-8 + 1e-8 * scale);
                expect(Math.abs(rb)).toBeLessThanOrEqual(1e-8 + 1e-8 * scale);
            });
    });

    it('is equivariant under a translation of the samples', () => {
        check(fc.tuple(domainArb, finite(-5, 5), finite(-5, 5), finite(-5, 5),
            vector(3, -20, 20)), ([domain, a, b, c, t]) => {
                const points = domain.map(([x, y]) =>
                    Vector.fromArray([x, y, a * x + b * y + c]));
                const shifted = points.map(q => Vector.fromArray([
                    q.get(0) + t.get(0), q.get(1) + t.get(1),
                    q.get(2) + t.get(2)]));

                const f0 = new ApprHeightPlane3();
                const f1 = new ApprHeightPlane3();
                expect(f0.fit(points)).toBe(true);
                expect(f1.fit(shifted)).toBe(true);
                for (let i = 0; i < 2; ++i) {
                    expectClose(f1.getParameters().coefficients.get(i),
                        f0.getParameters().coefficients.get(i), 1e-7, 1e-7);
                }
                for (let i = 0; i < 3; ++i) {
                    expectClose(f1.getParameters().average.get(i),
                        f0.getParameters().average.get(i) + t.get(i),
                        1e-8, 1e-8);
                }
            });
    });

    it('fails and zeroes the parameters for collinear or empty domains',
        () => {
            // The domain samples all lie on the y axis, so the 2x2 domain
            // covariance is singular. Integer coordinates and a power-of-two
            // count keep the averages exact (upstream divides by multiplying
            // with 1/n).
            check(fc.tuple(fc.constantFrom(1, 2, 4, 8),
                fc.array(fc.integer({ min: -8, max: 8 }),
                    { minLength: 8, maxLength: 8 }),
                fc.array(finite(-10, 10), { minLength: 8, maxLength: 8 })),
                ([n, ys, zs]) => {
                    const points = Array.from({ length: n },
                        (_, i) => Vector.fromArray([0, ys[i], zs[i]]));
                    const fitter = new ApprHeightPlane3();
                    expect(fitter.fit(points)).toBe(false);
                    const p = fitter.getParameters();
                    expect(p.average.values).toEqual([0, 0, 0]);
                    expect(p.coefficients.values).toEqual([0, 0, 0]);
                });

            const empty = new ApprHeightPlane3();
            expect(empty.fit([])).toBe(false);
            expect(empty.getParameters().average.values).toEqual([0, 0, 0]);
        });

    it('copyParameters deep-copies the parameters', () => {
        check(fc.tuple(domainArb, finite(-5, 5), finite(-5, 5)),
            ([domain, a, b]) => {
                const source = new ApprHeightPlane3();
                expect(source.fit(domain.map(([x, y]) =>
                    Vector.fromArray([x, y, a * x + b * y])))).toBe(true);
                const target = new ApprHeightPlane3();
                target.copyParameters(source);

                const s = source.getParameters();
                const t = target.getParameters();
                expect(t.average).not.toBe(s.average);
                expect(t.coefficients).not.toBe(s.coefficients);
                expectVectorClose(t.average, s.average, 0, 0);
                expectVectorClose(t.coefficients, s.coefficients, 0, 0);

                s.coefficients.set(0, 4321);
                expect(t.coefficients.get(0)).not.toBe(4321);
            });
    });
});
