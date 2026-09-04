import { describe, it, expect } from 'vitest';
import { FrenetFrame2, FrenetFrame3 } from '../src/FrenetFrame.js';
import { ParametricCurve } from '../src/ParametricCurve.js';
import { Vector, dot, length } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

// A counterclockwise circle of radius r in the plane: X(t) = r*(cos t, sin t).
class Circle2 extends ParametricCurve {
    constructor(private r: number, tmin: number, tmax: number) {
        super(2, tmin, tmax);
        this.mConstructed = true;
    }

    override evaluate(t: number, order: number, jet: Vector[]): void {
        const c = Math.cos(t), s = Math.sin(t);
        const r = this.r;
        jet[0] = Vector.fromArray([r * c, r * s]);
        if (order >= 1) { jet[1] = Vector.fromArray([-r * s, r * c]); }
        if (order >= 2) { jet[2] = Vector.fromArray([-r * c, -r * s]); }
        if (order >= 3) { jet[3] = Vector.fromArray([r * s, -r * c]); }
    }
}

// A straight line in the plane, which has zero curvature.
class Line2 extends ParametricCurve {
    constructor() {
        super(2, 0, 1);
        this.mConstructed = true;
    }

    override evaluate(t: number, order: number, jet: Vector[]): void {
        jet[0] = Vector.fromArray([1 + 2 * t, 3 - t]);
        if (order >= 1) { jet[1] = Vector.fromArray([2, -1]); }
        if (order >= 2) { jet[2] = Vector.fromArray([0, 0]); }
        if (order >= 3) { jet[3] = Vector.fromArray([0, 0]); }
    }
}

// A curve whose derivative vanishes everywhere (a constant point).
class Constant3 extends ParametricCurve {
    constructor() {
        super(3, 0, 1);
        this.mConstructed = true;
    }

    override evaluate(_t: number, _order: number, jet: Vector[]): void {
        for (let i = 0; i < 4; ++i) {
            jet[i] = Vector.fromArray([i === 0 ? 5 : 0, 0, 0]);
        }
    }
}

// The helix X(t) = (a cos t, a sin t, b t), with curvature a/(a^2+b^2) and
// torsion b/(a^2+b^2).
class Helix3 extends ParametricCurve {
    constructor(private a: number, private b: number) {
        super(3, 0, 4 * Math.PI);
        this.mConstructed = true;
    }

    override evaluate(t: number, order: number, jet: Vector[]): void {
        const c = Math.cos(t), s = Math.sin(t);
        const a = this.a, b = this.b;
        jet[0] = Vector.fromArray([a * c, a * s, b * t]);
        if (order >= 1) { jet[1] = Vector.fromArray([-a * s, a * c, b]); }
        if (order >= 2) { jet[2] = Vector.fromArray([-a * c, -a * s, 0]); }
        if (order >= 3) { jet[3] = Vector.fromArray([a * s, -a * c, 0]); }
    }
}

describe('FrenetFrame2', () => {
    it('rejects a curve whose dimension is not 2', () => {
        expect(() => new FrenetFrame2(new Helix3(1, 1)))
            .toThrow(/2-dimensional/);
    });

    it('exposes the curve', () => {
        const curve = new Circle2(2, 0, 2 * Math.PI);
        expect(new FrenetFrame2(curve).getCurve()).toBe(curve);
    });

    it('produces an orthonormal frame on a circle', () => {
        const frame = new FrenetFrame2(new Circle2(3, 0, 2 * Math.PI));
        for (const t of [0, 0.7, 1.9, 3.5, 5.6]) {
            const { position, tangent, normal } = frame.compute(t);
            expect(Math.hypot(tangent.values[0], tangent.values[1]))
                .toBeCloseTo(1, 12);
            expect(Math.hypot(normal.values[0], normal.values[1]))
                .toBeCloseTo(1, 12);
            expect(dot(tangent, normal)).toBeCloseTo(0, 12);
            expect(Math.hypot(position.values[0], position.values[1]))
                .toBeCloseTo(3, 12);
        }
    });

    it('has the outward radial normal on a counterclockwise circle', () => {
        // Perp(u0,u1) = (u1,-u0), so with T = (-sin t, cos t) the normal is
        // (cos t, sin t), the outward radial direction.
        const frame = new FrenetFrame2(new Circle2(2, 0, 2 * Math.PI));
        for (const t of [0.3, 1.2, 4.4]) {
            const { normal } = frame.compute(t);
            expect(normal.values[0]).toBeCloseTo(Math.cos(t), 12);
            expect(normal.values[1]).toBeCloseTo(Math.sin(t), 12);
        }
    });

    it('has curvature 1/r on a circle of radius r', () => {
        for (const r of [0.5, 1, 3, 10]) {
            const frame = new FrenetFrame2(new Circle2(r, 0, 2 * Math.PI));
            for (const t of [0, 1.1, 2.8, 5.0]) {
                expect(frame.getCurvature(t)).toBeCloseTo(1 / r, 11);
            }
        }
    });

    it('has zero curvature on a line', () => {
        const frame = new FrenetFrame2(new Line2());
        expect(frame.getCurvature(0.4)).toBeCloseTo(0, 14);
    });

    it('returns zero curvature when the speed is zero', () => {
        class Point2 extends ParametricCurve {
            constructor() { super(2, 0, 1); this.mConstructed = true; }
            override evaluate(_t: number, _o: number, jet: Vector[]): void {
                for (let i = 0; i < 4; ++i) { jet[i] = new Vector(2); }
            }
        }
        expect(new FrenetFrame2(new Point2()).getCurvature(0.5)).toBe(0);
    });
});

describe('FrenetFrame3', () => {
    it('rejects a curve whose dimension is not 3', () => {
        expect(() => new FrenetFrame3(new Circle2(1, 0, 1)))
            .toThrow(/3-dimensional/);
    });

    it('produces a right-handed orthonormal frame on a helix', () => {
        const frame = new FrenetFrame3(new Helix3(2, 0.75));
        for (const t of [0, 0.9, 2.6, 5.1, 9.4]) {
            const { tangent, normal, binormal } = frame.compute(t);
            for (const v of [tangent, normal, binormal]) {
                expect(Math.sqrt(dot(v, v))).toBeCloseTo(1, 12);
            }
            expect(dot(tangent, normal)).toBeCloseTo(0, 12);
            expect(dot(tangent, binormal)).toBeCloseTo(0, 12);
            expect(dot(normal, binormal)).toBeCloseTo(0, 12);

            // Right-handed: Cross(normal, binormal) = tangent.
            const t2 = cross(normal, binormal);
            for (let k = 0; k < 3; ++k) {
                expect(t2.values[k]).toBeCloseTo(tangent.values[k], 11);
            }
        }
    });

    it('has the inward normal on a helix', () => {
        const frame = new FrenetFrame3(new Helix3(2, 0.75));
        for (const t of [0.4, 2.2, 6.6]) {
            const { normal } = frame.compute(t);
            expect(normal.values[0]).toBeCloseTo(-Math.cos(t), 11);
            expect(normal.values[1]).toBeCloseTo(-Math.sin(t), 11);
            expect(normal.values[2]).toBeCloseTo(0, 11);
        }
    });

    it('has the known curvature and torsion of a helix', () => {
        for (const [a, b] of [[1, 1], [2, 0.5], [3, -1.5]]) {
            const frame = new FrenetFrame3(new Helix3(a, b));
            const denom = a * a + b * b;
            for (const t of [0, 1.3, 4.7, 8.8]) {
                expect(frame.getCurvature(t)).toBeCloseTo(a / denom, 11);
                expect(frame.getTorsion(t)).toBeCloseTo(b / denom, 11);
            }
        }
    });

    it('has zero torsion on a planar circle', () => {
        const frame = new FrenetFrame3(new Helix3(2, 0));
        for (const t of [0.2, 1.8, 3.9]) {
            expect(frame.getCurvature(t)).toBeCloseTo(0.5, 11);
            expect(frame.getTorsion(t)).toBeCloseTo(0, 11);
        }
    });

    it('returns zero curvature and torsion for a degenerate curve', () => {
        const frame = new FrenetFrame3(new Constant3());
        expect(frame.getCurvature(0.5)).toBe(0);
        expect(frame.getTorsion(0.5)).toBe(0);
    });

    it('returns the position from the curve', () => {
        const frame = new FrenetFrame3(new Helix3(2, 1));
        const { position } = frame.compute(1.25);
        expect(position.values[0]).toBeCloseTo(2 * Math.cos(1.25), 12);
        expect(position.values[1]).toBeCloseTo(2 * Math.sin(1.25), 12);
        expect(position.values[2]).toBeCloseTo(1.25, 12);
    });
});


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). FrenetFrame.h was read line by
// line against src/FrenetFrame.ts. Curvature and torsion are geometric
// invariants: they do not depend on the parameterization and they are
// unchanged by a rigid motion (torsion is a signed pseudo-scalar, so it is
// unchanged by a rotation). Those two properties are checked against
// independently reparameterized and rotated copies of random polynomial
// curves, which is a genuine cross-check of the powers and normalizations in
// GetCurvature/GetTorsion that no restatement of the formulas would give.
import {
    check, fc, expectClose, expectVectorClose, wellScaledVector,
    rotationFrame
} from './helpers/arbitraries.js';

// X(t) = C0 + C1 t + C2 t^2 + C3 t^3, with exact derivatives.
class PolyCurve extends ParametricCurve {
    constructor(dimension: number, private c: Vector[]) {
        super(dimension, -1, 1);
        this.mConstructed = true;
    }

    override evaluate(t: number, order: number, jet: Vector[]): void {
        const n = this.mDimension;
        const at = (r: number): Vector => {
            const v = new Vector(n);
            for (let k = 0; k < n; ++k) {
                let sum = 0;
                for (let j = r; j < this.c.length; ++j) {
                    let factor = 1;
                    for (let q = 0; q < r; ++q) { factor *= j - q; }
                    sum += factor * this.c[j].values[k] * Math.pow(t, j - r);
                }
                v.values[k] = sum;
            }
            return v;
        };
        for (let r = 0; r <= Math.min(order, 3); ++r) { jet[r] = at(r); }
    }
}

// Y(u) = X(s(u)) with s(u) = u + beta*u^2, whose derivative 1 + 2*beta*u is
// positive on [-1,1] for |beta| < 1/2. The jet follows from the chain rule.
class ReparamCurve extends ParametricCurve {
    constructor(private base: ParametricCurve, private beta: number) {
        super(base.getDimension(), -0.5, 0.5);
        this.mConstructed = true;
    }

    s(u: number): number {
        return u + this.beta * u * u;
    }

    override evaluate(u: number, order: number, jet: Vector[]): void {
        const inner = this.base.createJet();
        this.base.evaluate(this.s(u), 3, inner);
        const d1 = 1 + 2 * this.beta * u;
        const d2 = 2 * this.beta;
        const n = this.mDimension;
        const make = (f: (k: number) => number): Vector => {
            const v = new Vector(n);
            for (let k = 0; k < n; ++k) { v.values[k] = f(k); }
            return v;
        };
        jet[0] = make(k => inner[0].values[k]);
        if (order >= 1) {
            jet[1] = make(k => inner[1].values[k] * d1);
        }
        if (order >= 2) {
            jet[2] = make(k => inner[2].values[k] * d1 * d1 +
                inner[1].values[k] * d2);
        }
        if (order >= 3) {
            jet[3] = make(k => inner[3].values[k] * d1 * d1 * d1 +
                3 * inner[2].values[k] * d1 * d2);
        }
    }
}

// Coefficients whose leading terms are big enough that the curve is not
// nearly straight: a nearly zero cross product makes curvature and torsion
// ill conditioned and the invariance properties meaningless.
const polyCoefficients = (dim: number): fc.Arbitrary<Vector[]> =>
    fc.array(wellScaledVector(dim, -3, 3), { minLength: 4, maxLength: 4 })
        .filter(c => {
            for (let j = 1; j <= 3; ++j) {
                let sum = 0;
                for (const x of c[j].values) { sum += x * x; }
                if (sum < 0.25) { return false; }
            }
            return true;
        });

const param = fc.double({ min: -0.4, max: 0.4, noNaN: true,
    noDefaultInfinity: true });

// The curvature and torsion formulas divide by |V|^3 and |V x A|^2, so a
// nearly stationary or nearly straight point makes them ill conditioned and
// the invariance comparisons meaningless.
function wellConditioned(curve: ParametricCurve, t: number): boolean {
    const jet = curve.createJet();
    curve.evaluate(t, 2, jet);
    if (length(jet[1]) < 0.3) { return false; }
    if (curve.getDimension() === 3) {
        return length(cross(jet[1], jet[2])) > 0.05;
    }
    return Math.abs(jet[1].values[0] * jet[2].values[1] -
        jet[1].values[1] * jet[2].values[0]) > 0.05;
}

describe('FrenetFrame verification', () => {
    it('FrenetFrame3 curvature and torsion are parameterization invariant', () => {
        check(fc.tuple(polyCoefficients(3), param,
            fc.double({ min: -0.4, max: 0.4, noNaN: true,
                noDefaultInfinity: true })), ([c, u, beta]) => {
            const base = new PolyCurve(3, c);
            const reparam = new ReparamCurve(base, beta);
            const t = reparam.s(u);
            const a = new FrenetFrame3(base);
            const b = new FrenetFrame3(reparam);
            if (!wellConditioned(base, t)) { return; }
            const kappa = a.getCurvature(t);
            // The reparameterization is a smooth monotone change of variable,
            // so the two curvature values are the same number computed along
            // different arithmetic paths.
            expectClose(b.getCurvature(u), kappa, 1e-8, 1e-8);
            expectClose(b.getTorsion(u), a.getTorsion(t), 1e-7, 1e-7);

            // The frame is likewise unchanged (s'(u) > 0 keeps the tangent
            // direction).
            const fa = a.compute(t);
            const fb = b.compute(u);
            expectVectorClose(fb.tangent, fa.tangent, 1e-8, 1e-8);
            expectVectorClose(fb.normal, fa.normal, 1e-7, 1e-7);
            expectVectorClose(fb.binormal, fa.binormal, 1e-7, 1e-7);
        });
    });

    it('FrenetFrame2 curvature is parameterization invariant', () => {
        check(fc.tuple(polyCoefficients(2), param,
            fc.double({ min: -0.4, max: 0.4, noNaN: true,
                noDefaultInfinity: true })), ([c, u, beta]) => {
            const base = new PolyCurve(2, c);
            const reparam = new ReparamCurve(base, beta);
            const t = reparam.s(u);
            const a = new FrenetFrame2(base);
            const b = new FrenetFrame2(reparam);
            if (!wellConditioned(base, t)) { return; }
            expectClose(b.getCurvature(u), a.getCurvature(t), 1e-8, 1e-8);
            const fa = a.compute(t);
            const fb = b.compute(u);
            expectVectorClose(fb.tangent, fa.tangent, 1e-8, 1e-8);
            expectVectorClose(fb.normal, fa.normal, 1e-8, 1e-8);
        });
    });

    it('the 3D frame is right-handed and orthonormal', () => {
        check(fc.tuple(polyCoefficients(3), param), ([c, t]) => {
            const curve = new PolyCurve(3, c);
            if (!wellConditioned(curve, t)) { return; }
            const frame = new FrenetFrame3(curve);
            const { tangent, normal, binormal } = frame.compute(t);
            for (const v of [tangent, normal, binormal]) {
                expectClose(length(v), 1, 1e-9, 1e-9);
            }
            expectClose(dot(tangent, normal), 0, 1e-9, 1e-9);
            expectClose(dot(tangent, binormal), 0, 1e-9, 1e-9);
            expectClose(dot(normal, binormal), 0, 1e-9, 1e-9);
            // Right handed: Cross(normal, binormal) is the tangent.
            expectVectorClose(cross(normal, binormal), tangent, 1e-8, 1e-8);
            // The normal points toward the center of curvature, so it has a
            // positive component along the acceleration.
            const jet = frame.getCurve().createJet();
            frame.getCurve().evaluate(t, 2, jet);
            expect(dot(normal, jet[2])).toBeGreaterThan(-1e-12);
        });
    });

    it('the 2D frame is orthonormal with the clockwise-rotated normal', () => {
        check(fc.tuple(polyCoefficients(2), param), ([c, t]) => {
            const curve = new PolyCurve(2, c);
            if (!wellConditioned(curve, t)) { return; }
            const frame = new FrenetFrame2(curve);
            const { tangent, normal } = frame.compute(t);
            expectClose(length(tangent), 1, 1e-9, 1e-9);
            expectClose(length(normal), 1, 1e-9, 1e-9);
            expectClose(dot(tangent, normal), 0, 1e-12, 1e-12);
            // Perp(x,y) = (y,-x) is a clockwise rotation by pi/2.
            expectClose(normal.values[0], tangent.values[1], 0, 0);
            expectClose(normal.values[1], -tangent.values[0], 0, 0);
        });
    });

    it('curvature and torsion are invariant under rotations', () => {
        check(fc.tuple(polyCoefficients(3), param, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([c, t, R, shift]) => {
            const moved = c.map((v, j) => Vector.fromArray([0, 1, 2].map(k =>
                R[0].values[k] * v.values[0] + R[1].values[k] * v.values[1] +
                R[2].values[k] * v.values[2] + (j === 0 ? shift.values[k] : 0))));
            const a = new FrenetFrame3(new PolyCurve(3, c));
            const b = new FrenetFrame3(new PolyCurve(3, moved));
            if (!wellConditioned(a.getCurve(), t)) { return; }
            const kappa = a.getCurvature(t);
            expectClose(b.getCurvature(t), kappa, 1e-9, 1e-9);
            expectClose(b.getTorsion(t), a.getTorsion(t), 1e-8, 1e-8);
        });
    });

    it('returns zero curvature and torsion at a stationary point', () => {
        // Upstream guards the division with speedSqr > 0 and Dot(cross,cross)
        // > 0 and documents the indeterminate case as zero.
        const curve = new PolyCurve(3, [
            Vector.fromArray([1, 2, 3]), new Vector(3), new Vector(3),
            new Vector(3)]);
        const frame = new FrenetFrame3(curve);
        expect(frame.getCurvature(0.3)).toBe(0);
        expect(frame.getTorsion(0.3)).toBe(0);
        const flat = new FrenetFrame2(new PolyCurve(2, [
            Vector.fromArray([1, 2]), new Vector(2), new Vector(2),
            new Vector(2)]));
        expect(flat.getCurvature(0.3)).toBe(0);
    });
});
