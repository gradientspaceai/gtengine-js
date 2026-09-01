import { describe, it, expect } from 'vitest';
import { ImplicitCurve2 } from '../src/ImplicitCurve2';
import { Vector, dot, length } from '../src/Vector';
import { dotPerp } from '../src/Vector2';

// F(x,y) = x^2 + y^2 - r^2, the circle of radius r centered at the origin.
class Circle2Implicit extends ImplicitCurve2 {
    constructor(private r: number) {
        super();
    }

    override f(p: Vector): number {
        return p.get(0) * p.get(0) + p.get(1) * p.get(1) - this.r * this.r;
    }
    override fx(p: Vector): number { return 2 * p.get(0); }
    override fy(p: Vector): number { return 2 * p.get(1); }
    override fxx(): number { return 2; }
    override fxy(): number { return 0; }
    override fyy(): number { return 2; }
}

// F(x,y) = x^2/a^2 + y^2/b^2 - 1, the axis-aligned ellipse.
class Ellipse2Implicit extends ImplicitCurve2 {
    constructor(private a: number, private b: number) {
        super();
    }

    override f(p: Vector): number {
        const x = p.get(0), y = p.get(1);
        return x * x / (this.a * this.a) + y * y / (this.b * this.b) - 1;
    }
    override fx(p: Vector): number { return 2 * p.get(0) / (this.a * this.a); }
    override fy(p: Vector): number { return 2 * p.get(1) / (this.b * this.b); }
    override fxx(): number { return 2 / (this.a * this.a); }
    override fxy(): number { return 0; }
    override fyy(): number { return 2 / (this.b * this.b); }
}

// A line through the origin: F(x,y) = 2x - 3y. The gradient is constant and
// nonzero, so the curvature is zero.
class Line2Implicit extends ImplicitCurve2 {
    constructor() {
        super();
    }

    override f(p: Vector): number { return 2 * p.get(0) - 3 * p.get(1); }
    override fx(): number { return 2; }
    override fy(): number { return -3; }
    override fxx(): number { return 0; }
    override fxy(): number { return 0; }
    override fyy(): number { return 0; }
}

// A degenerate example whose gradient vanishes at the origin:
// F(x,y) = x^2 + y^2.
class Degenerate2Implicit extends ImplicitCurve2 {
    constructor() {
        super();
    }

    override f(p: Vector): number {
        return p.get(0) * p.get(0) + p.get(1) * p.get(1);
    }
    override fx(p: Vector): number { return 2 * p.get(0); }
    override fy(p: Vector): number { return 2 * p.get(1); }
    override fxx(): number { return 2; }
    override fxy(): number { return 0; }
    override fyy(): number { return 2; }
}

function vec2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

describe('ImplicitCurve2', () => {
    it('reports whether a point is on the curve', () => {
        const curve = new Circle2Implicit(2);
        expect(curve.isOnCurve(vec2(2, 0), 1e-12)).toBe(true);
        expect(curve.isOnCurve(vec2(Math.SQRT2, Math.SQRT2), 1e-12)).toBe(true);
        expect(curve.isOnCurve(vec2(2.001, 0), 1e-12)).toBe(false);
        expect(curve.isOnCurve(vec2(2.001, 0), 1e-2)).toBe(true);
    });

    it('matches finite differences of F for the gradient', () => {
        const curve = new Ellipse2Implicit(3, 1.5);
        const h = 1e-5;
        for (const t of [0.3, 1.1, 2.7, 4.4, 5.9]) {
            const p = vec2(3 * Math.cos(t), 1.5 * Math.sin(t));
            const g = curve.getGradient(p);
            const dx = (curve.f(vec2(p.get(0) + h, p.get(1)))
                - curve.f(vec2(p.get(0) - h, p.get(1)))) / (2 * h);
            const dy = (curve.f(vec2(p.get(0), p.get(1) + h))
                - curve.f(vec2(p.get(0), p.get(1) - h))) / (2 * h);
            expect(g.get(0)).toBeCloseTo(dx, 8);
            expect(g.get(1)).toBeCloseTo(dy, 8);
        }
    });

    it('matches finite differences of the gradient for the Hessian', () => {
        const curve = new Ellipse2Implicit(2, 0.75);
        const h = 1e-5;
        const p = vec2(2 * Math.cos(0.8), 0.75 * Math.sin(0.8));
        const hess = curve.getHessian(p);
        const gxp = curve.getGradient(vec2(p.get(0) + h, p.get(1)));
        const gxm = curve.getGradient(vec2(p.get(0) - h, p.get(1)));
        const gyp = curve.getGradient(vec2(p.get(0), p.get(1) + h));
        const gym = curve.getGradient(vec2(p.get(0), p.get(1) - h));
        expect(hess.get(0, 0)).toBeCloseTo((gxp.get(0) - gxm.get(0)) / (2 * h), 8);
        expect(hess.get(0, 1)).toBeCloseTo((gyp.get(0) - gym.get(0)) / (2 * h), 8);
        expect(hess.get(1, 0)).toBeCloseTo((gxp.get(1) - gxm.get(1)) / (2 * h), 8);
        expect(hess.get(1, 1)).toBeCloseTo((gyp.get(1) - gym.get(1)) / (2 * h), 8);
        // The Hessian is symmetric.
        expect(hess.get(0, 1)).toBe(hess.get(1, 0));
    });

    it('computes an orthonormal frame whose normal is the unit gradient', () => {
        const curve = new Circle2Implicit(2.5);
        for (const t of [0, 0.4, 1.9, 3.6, 5.2]) {
            const p = vec2(2.5 * Math.cos(t), 2.5 * Math.sin(t));
            const frame = curve.getFrame(p);
            expect(length(frame.tangent)).toBeCloseTo(1, 12);
            expect(length(frame.normal)).toBeCloseTo(1, 12);
            expect(dot(frame.tangent, frame.normal)).toBeCloseTo(0, 12);
            // The normal is the normalized gradient, which for the circle
            // points radially outward.
            expect(frame.normal.get(0)).toBeCloseTo(Math.cos(t), 12);
            expect(frame.normal.get(1)).toBeCloseTo(Math.sin(t), 12);
            // The tangent is orthogonal to the gradient, so it is tangent to
            // the level set.
            expect(dot(frame.tangent, curve.getGradient(p))).toBeCloseTo(0, 10);
            // The frame is a consistently oriented basis of the plane.
            expect(Math.abs(dotPerp(frame.tangent, frame.normal)))
                .toBeCloseTo(1, 12);
        }
    });

    it('computes curvature -1/r on the circle F = x^2 + y^2 - r^2', () => {
        for (const r of [0.5, 1, 4]) {
            const curve = new Circle2Implicit(r);
            for (const t of [0, 0.9, 2.2, 4.8]) {
                const p = vec2(r * Math.cos(t), r * Math.sin(t));
                const result = curve.getCurvature(p);
                expect(result.valid).toBe(true);
                // The sign follows the orientation of the gradient; the
                // magnitude is 1/r.
                expect(result.curvature).toBeCloseTo(-1 / r, 12);
                expect(Math.abs(result.curvature)).toBeCloseTo(1 / r, 12);
            }
        }
    });

    it('computes the analytic curvature of an ellipse', () => {
        const a = 3, b = 1.25;
        const curve = new Ellipse2Implicit(a, b);
        for (const t of [0, 0.6, 1.5708, 2.9, 4.1, 5.7]) {
            const p = vec2(a * Math.cos(t), b * Math.sin(t));
            const result = curve.getCurvature(p);
            expect(result.valid).toBe(true);
            // kappa = a*b / (a^2 sin^2 t + b^2 cos^2 t)^{3/2}, negated by the
            // gradient orientation of F.
            const denom = Math.pow(a * a * Math.sin(t) * Math.sin(t)
                + b * b * Math.cos(t) * Math.cos(t), 1.5);
            expect(result.curvature).toBeCloseTo(-a * b / denom, 10);
        }
    });

    it('computes zero curvature on a line', () => {
        const curve = new Line2Implicit();
        const result = curve.getCurvature(vec2(3, 2));
        expect(result.valid).toBe(true);
        expect(result.curvature).toBeCloseTo(0, 14);
        const frame = curve.getFrame(vec2(3, 2));
        expect(length(frame.tangent)).toBeCloseTo(1, 12);
        // The line direction is (3,2)/sqrt(13); the tangent is parallel to it.
        expect(Math.abs(dot(frame.tangent,
            vec2(3 / Math.sqrt(13), 2 / Math.sqrt(13))))).toBeCloseTo(1, 12);
    });

    it('reports an invalid curvature where the gradient vanishes', () => {
        const curve = new Degenerate2Implicit();
        const result = curve.getCurvature(vec2(0, 0));
        expect(result.valid).toBe(false);
        expect(result.curvature).toBe(0);
    });

    it('agrees with the parametric curvature formula (randomized)', () => {
        // For F(x,y) = x^2/a^2 + y^2/b^2 - 1, compare the implicit curvature
        // with the parametric formula for X(t) = (a cos t, b sin t):
        // kappa = (x' y'' - y' x'') / (x'^2 + y'^2)^{3/2}.
        let seed = 987654321;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const a = 0.5 + 3 * rand();
            const b = 0.5 + 3 * rand();
            const t = 2 * Math.PI * rand();
            const curve = new Ellipse2Implicit(a, b);
            const p = vec2(a * Math.cos(t), b * Math.sin(t));
            expect(curve.isOnCurve(p, 1e-12)).toBe(true);
            const result = curve.getCurvature(p);
            const xp = -a * Math.sin(t), yp = b * Math.cos(t);
            const xpp = -a * Math.cos(t), ypp = -b * Math.sin(t);
            const kappa = (xp * ypp - yp * xpp)
                / Math.pow(xp * xp + yp * yp, 1.5);
            expect(result.valid).toBe(true);
            // The implicit gradient points outward, the parametric traversal
            // is counterclockwise, so the two differ by a sign.
            expect(result.curvature).toBeCloseTo(-kappa, 9);
        }
    });
});
