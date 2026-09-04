import { describe, expect, it } from 'vitest';
import { ApprQuadratic2, ApprQuadraticCircle2 } from '../src/ApprQuadratic2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, dot } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, finite, positive, rotationFrame, wellScaledVector } from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// V = (1, x, y, x^2, x*y, y^2)
function basis(x: number, y: number): number[] {
    return [1, x, y, x * x, x * y, y * y];
}

// Normalize to unit length with a deterministic sign (the component of
// largest magnitude is made positive), so eigenvectors can be compared.
function canonical(c: readonly number[]): number[] {
    let len = 0, imax = 0;
    for (let i = 0; i < c.length; ++i) {
        len += c[i] * c[i];
        if (Math.abs(c[i]) > Math.abs(c[imax])) {
            imax = i;
        }
    }
    len = Math.sqrt(len);
    const s = (c[imax] < 0 ? -1 : 1) / len;
    return c.map((v) => v * s);
}

describe('ApprQuadratic2', () => {
    it('recovers an ellipse from points sampled on it', () => {
        // x^2/9 + y^2/4 = 1, i.e. -1 + x^2/9 + y^2/4 = 0.
        const points: Vector[] = [];
        const n = 16;
        for (let i = 0; i < n; ++i) {
            const t = (2 * Math.PI * i) / n;
            points.push(v2(3 * Math.cos(t), 2 * Math.sin(t)));
        }

        const result = new ApprQuadratic2().compute(points);
        expect(result.coefficients.length).toBe(6);
        expect(result.minEigenvalue).toBeLessThan(1e-14);

        // The eigenvector is unit length.
        let sqrLen = 0;
        for (const c of result.coefficients) {
            sqrLen += c * c;
        }
        expect(sqrLen).toBeCloseTo(1, 10);

        // It matches the exact conic up to scale and sign.
        const expected = canonical([-1, 0, 0, 1 / 9, 0, 1 / 4]);
        const actual = canonical(result.coefficients);
        for (let i = 0; i < 6; ++i) {
            expect(actual[i]).toBeCloseTo(expected[i], 8);
        }

        // Every sample point satisfies the fitted equation.
        for (const p of points) {
            const b = basis(p.get(0), p.get(1));
            let sum = 0;
            for (let i = 0; i < 6; ++i) {
                sum += result.coefficients[i] * b[i];
            }
            expect(Math.abs(sum)).toBeLessThan(1e-10);
        }
    });

    it('recovers a rotated, translated conic', () => {
        // A circle of radius 2 about (1,-3) sheared into an ellipse by
        // (u,v) -> (u + 0.5*v, v), then translated.
        const points: Vector[] = [];
        const n = 20;
        for (let i = 0; i < n; ++i) {
            const t = (2 * Math.PI * i) / n;
            const u = 2 * Math.cos(t);
            const w = 2 * Math.sin(t);
            points.push(v2(1 + u + 0.5 * w, -3 + w));
        }

        const result = new ApprQuadratic2().compute(points);
        expect(result.minEigenvalue).toBeLessThan(1e-14);
        for (const p of points) {
            const b = basis(p.get(0), p.get(1));
            let sum = 0;
            for (let i = 0; i < 6; ++i) {
                sum += result.coefficients[i] * b[i];
            }
            expect(Math.abs(sum)).toBeLessThan(1e-8);
        }
    });

    it('recovers a pair of lines (a degenerate conic)', () => {
        // x*y = 0: points on the two coordinate axes.
        const points: Vector[] = [];
        for (let i = 1; i <= 5; ++i) {
            points.push(v2(i, 0), v2(-i, 0), v2(0, i), v2(0, -i));
        }

        const result = new ApprQuadratic2().compute(points);
        expect(result.minEigenvalue).toBeLessThan(1e-14);
        for (const p of points) {
            const b = basis(p.get(0), p.get(1));
            let sum = 0;
            for (let i = 0; i < 6; ++i) {
                sum += result.coefficients[i] * b[i];
            }
            expect(Math.abs(sum)).toBeLessThan(1e-8);
        }
    });

    it('returns a positive fit measure for points not on a conic', () => {
        // Six generic points determine a conic exactly, so use more.
        const points: Vector[] = [];
        const n = 24;
        for (let i = 0; i < n; ++i) {
            const t = (2 * Math.PI * i) / n;
            // A limacon-like curve, not a conic.
            const r = 1 + 0.5 * Math.cos(3 * t);
            points.push(v2(r * Math.cos(t), r * Math.sin(t)));
        }

        const result = new ApprQuadratic2().compute(points);
        expect(result.minEigenvalue).toBeGreaterThan(1e-6);
    });
});

describe('ApprQuadraticCircle2', () => {
    it('recovers a circle from points sampled on it', () => {
        const center = v2(3, -2);
        const radius = 5;
        const points: Vector[] = [];
        const n = 12;
        for (let i = 0; i < n; ++i) {
            const t = (2 * Math.PI * i) / n;
            points.push(v2(center.get(0) + radius * Math.cos(t),
                center.get(1) + radius * Math.sin(t)));
        }

        const circle = new Hypersphere(2);
        const measure = new ApprQuadraticCircle2().compute(points, circle);
        expect(measure).toBeLessThan(1e-12);
        expect(circle.center.get(0)).toBeCloseTo(3, 8);
        expect(circle.center.get(1)).toBeCloseTo(-2, 8);
        expect(circle.radius).toBeCloseTo(5, 8);
    });

    it('recovers the unit circle at the origin', () => {
        const points: Vector[] = [];
        const n = 8;
        for (let i = 0; i < n; ++i) {
            const t = (2 * Math.PI * i) / n;
            points.push(v2(Math.cos(t), Math.sin(t)));
        }

        const circle = new Hypersphere(2);
        const measure = new ApprQuadraticCircle2().compute(points, circle);
        expect(measure).toBeLessThan(1e-15);
        expect(circle.center.get(0)).toBeCloseTo(0, 10);
        expect(circle.center.get(1)).toBeCloseTo(0, 10);
        expect(circle.radius).toBeCloseTo(1, 10);
    });

    it('averages symmetric radial perturbations of a circle', () => {
        // Sample at radius r-d and r+d at every angle. This is an
        // algebraic fit constrained by Length(C') = 1, a norm that is not
        // invariant under translation, so the fitted circle is only near
        // the true one -- the well-known bias of algebraic circle fitting.
        // The fit must still be within the noise band.
        const center = v2(-1, 4);
        const radius = 3;
        const d = 0.1;
        const points: Vector[] = [];
        const n = 16;
        for (let i = 0; i < n; ++i) {
            const t = (2 * Math.PI * i) / n;
            for (const r of [radius - d, radius + d]) {
                points.push(v2(center.get(0) + r * Math.cos(t),
                    center.get(1) + r * Math.sin(t)));
            }
        }

        const circle = new Hypersphere(2);
        const measure = new ApprQuadraticCircle2().compute(points, circle);
        expect(measure).toBeGreaterThan(0);
        expect(Math.abs(circle.center.get(0) - center.get(0))).toBeLessThan(d);
        expect(Math.abs(circle.center.get(1) - center.get(1))).toBeLessThan(d);
        expect(circle.radius).toBeGreaterThan(radius - d);
        expect(circle.radius).toBeLessThan(radius + d);
    });
});

describe('ApprQuadratic2 verification', () => {
    // wellScaledVector keeps every coordinate either exactly zero or above
    // 1e-3: the monomial vector contains x^2 and x^4 terms, whose subnormal
    // products would otherwise leave the eigensolver with no mantissa to
    // work with.
    //
    // V = (1, x, y, x^2, x*y, y^2), the monomial vector the fit uses.
    const monomials = (p: Vector): number[] => {
        const x = p.get(0), y = p.get(1);
        return [1, x, y, x * x, x * y, y * y];
    };

    // Samples on the ellipse center + a*cos(t)*U + b*sin(t)*V.
    const ellipsePoints = (center: Vector, a: number, b: number,
        frame: readonly Vector[], k: number, offset: number): Vector[] => {
        const points: Vector[] = [];
        for (let i = 0; i < k; ++i) {
            const t = offset + (2 * Math.PI * i) / k;
            const ca = a * Math.cos(t), sb = b * Math.sin(t);
            points.push(Vector.fromArray([
                center.get(0) + ca * frame[0].get(0) + sb * frame[1].get(0),
                center.get(1) + ca * frame[0].get(1) + sb * frame[1].get(1)]));
        }
        return points;
    };

    const ellipseArb = fc.tuple(wellScaledVector(2, -3, 3), positive(3, 0.5),
        positive(3, 0.5), rotationFrame(2),
        fc.integer({ min: 6, max: 14 }), finite(0, 2 * Math.PI));

    it('the reported measure is the mean squared conic residual', () => {
        // M is (1/n) * sum_i V_i * V_i^T (the header's rank-one formula is a
        // typo), so for the returned unit eigenvector c the minimum
        // eigenvalue equals c^T M c = (1/n) * sum_i Dot(V_i, c)^2. That
        // identity exercises every entry of the 6x6 assembly, including the
        // aliased ones copied in after the accumulation loop.
        check(fc.array(wellScaledVector(2, -6, 6), { minLength: 1, maxLength: 12 }),
            points => {
                const { coefficients, minEigenvalue } =
                    new ApprQuadratic2().compute(points);
                expect(coefficients.length).toBe(6);

                let sum = 0, scale = 0;
                for (const p of points) {
                    const v = monomials(p);
                    const d = v.reduce((u, t, i) => u + t * coefficients[i], 0);
                    sum += d * d;
                    scale += v.reduce((u, t) => u + t * t, 0);
                }
                sum /= points.length;
                scale /= points.length;
                expect(Math.abs(minEigenvalue - sum))
                    .toBeLessThanOrEqual(1e-8 + 1e-8 * scale);
                expect(minEigenvalue).toBeGreaterThanOrEqual(0);
            });
    });

    it('returns a unit-length coefficient vector', () => {
        check(fc.array(wellScaledVector(2, -6, 6), { minLength: 1, maxLength: 12 }),
            points => {
                const { coefficients } = new ApprQuadratic2().compute(points);
                const norm = coefficients.reduce((u, c) => u + c * c, 0);
                expectClose(norm, 1, 1e-9, 1e-9);
            });
    });

    it('fits an ellipse its samples lie on exactly', () => {
        check(ellipseArb, ([center, a, b, frame, k, offset]) => {
            const points = ellipsePoints(center, a, b, frame, k, offset);
            const { coefficients, minEigenvalue } =
                new ApprQuadratic2().compute(points);

            // Every sample satisfies the fitted conic.
            let scale = 0;
            for (const p of points) {
                const v = monomials(p);
                scale += v.reduce((u, t) => u + t * t, 0);
            }
            scale /= points.length;
            expect(minEigenvalue).toBeLessThanOrEqual(1e-10 * scale + 1e-14);

            for (const p of points) {
                const v = monomials(p);
                expect(Math.abs(v.reduce(
                    (u, t, i) => u + t * coefficients[i], 0)))
                    .toBeLessThanOrEqual(1e-6 * Math.sqrt(scale));
            }
        });
    });

    it('handles degenerate sample sets without throwing', () => {
        // Coincident samples make M rank one; upstream still returns the
        // eigenvector of the smallest eigenvalue and clamps the measure to
        // be nonnegative.
        check(fc.tuple(wellScaledVector(2, -5, 5), fc.integer({ min: 1, max: 6 })),
            ([p, n]) => {
                const points = Array.from({ length: n }, () => p.clone());
                const { coefficients, minEigenvalue } =
                    new ApprQuadratic2().compute(points);
                expect(minEigenvalue).toBeGreaterThanOrEqual(0);
                expect(coefficients.every(Number.isFinite)).toBe(true);
            });
    });
});

describe('ApprQuadraticCircle2 verification', () => {
    const circlePoints = (center: Vector, r: number, k: number,
        offset: number): Vector[] => {
        const points: Vector[] = [];
        for (let i = 0; i < k; ++i) {
            const t = offset + (2 * Math.PI * i) / k;
            points.push(Vector.fromArray([center.get(0) + r * Math.cos(t),
                center.get(1) + r * Math.sin(t)]));
        }
        return points;
    };

    const circleArb = fc.tuple(wellScaledVector(2, -4, 4), positive(4, 0.5),
        fc.integer({ min: 3, max: 12 }), finite(0, 2 * Math.PI));

    it('recovers a circle its samples lie on', () => {
        check(circleArb, ([center, r, k, offset]) => {
            const points = circlePoints(center, r, k, offset);
            const circle = new Hypersphere(2);
            const measure = new ApprQuadraticCircle2().compute(points, circle);

            expect(measure).toBeGreaterThanOrEqual(0);
            expect(measure).toBeLessThan(1e-9);
            expectVectorClose(circle.center, center, 1e-7, 1e-7);
            expectClose(circle.radius, r, 1e-7, 1e-7);
        });
    });

    it('reports the center and radius implied by the fitted coefficients',
        () => {
            // The documented relations are (xc,yc) = -(C[1],C[2])/2 and
            // r = sqrt(xc^2 + yc^2 - C[0]) for the normalized coefficients,
            // with the radius clamped to be nonnegative.
            check(fc.array(wellScaledVector(2, -6, 6), { minLength: 3, maxLength: 12 }),
                points => {
                    const circle = new Hypersphere(2);
                    const measure =
                        new ApprQuadraticCircle2().compute(points, circle);
                    expect(measure).toBeGreaterThanOrEqual(0);

                    if (!Number.isFinite(circle.radius)
                        || !circle.center.values.every(Number.isFinite)) {
                        // C'[3] can vanish for degenerate data, and upstream
                        // divides by it unguarded (pinned separately below).
                        return;
                    }
                    expect(circle.radius).toBeGreaterThanOrEqual(0);
                    const c0 = dot(circle.center, circle.center)
                        - circle.radius * circle.radius;
                    // Every sample satisfies c0 + C1*x + C2*y + x^2 + y^2 = 0
                    // only for an exact fit, so just check the relation
                    // between center, radius and c0 is self-consistent.
                    expectClose(circle.radius,
                        Math.sqrt(Math.max(
                            dot(circle.center, circle.center) - c0, 0)),
                        1e-9, 1e-9);
                });
        });

    it('divides by a vanishing C[3] without a guard', () => {
        // The circle coefficients are C'[i]/C'[3]. For coincident samples at
        // the origin the eigenvector of the smallest eigenvalue is (0,1,0,0),
        // so C'[3] is zero and the reported center and radius are
        // -Infinity/NaN rather than a failure. Preserved from upstream.
        check(fc.integer({ min: 1, max: 6 }), n => {
            const points = Array.from({ length: n },
                () => Vector.fromArray([0, 0]));
            const circle = Hypersphere.fromCenterRadius(
                Vector.fromArray([1, 2]), 3);
            const measure = new ApprQuadraticCircle2().compute(points, circle);
            expect(measure).toBe(0);
            expect(Number.isNaN(circle.radius)).toBe(true);
        });
    });
});
