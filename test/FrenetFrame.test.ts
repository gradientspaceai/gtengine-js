import { describe, it, expect } from 'vitest';
import { FrenetFrame2, FrenetFrame3 } from '../src/FrenetFrame';
import { ParametricCurve } from '../src/ParametricCurve';
import { Vector, dot } from '../src/Vector';
import { cross } from '../src/Vector3';

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
