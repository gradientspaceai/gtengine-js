import { describe, expect, it } from 'vitest';
import { ApprQuadratic2, ApprQuadraticCircle2 } from '../src/ApprQuadratic2';
import { Hypersphere } from '../src/Hypersphere';
import { Vector } from '../src/Vector';

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
