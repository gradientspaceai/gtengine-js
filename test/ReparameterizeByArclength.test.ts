import { describe, it, expect } from 'vitest';
import { ReparameterizeByArclength } from '../src/ReparameterizeByArclength';
import { BezierCurve } from '../src/BezierCurve';
import { ParametricCurve } from '../src/ParametricCurve';
import { Vector } from '../src/Vector';

function vec(...values: number[]): Vector {
    const v = new Vector(values.length);
    for (let i = 0; i < values.length; ++i) {
        v.values[i] = values[i];
    }
    return v;
}

// A circle of radius r traversed on [0, 2*pi]; its speed is the constant r,
// so arclength s corresponds exactly to t = s / r.
class Circle2 extends ParametricCurve {
    constructor(private r: number) {
        super(2, 0, 2 * Math.PI);
        this.mConstructed = true;
    }

    override evaluate(t: number, order: number, jet: Vector[]): void {
        const c = Math.cos(t), s = Math.sin(t), r = this.r;
        jet[0] = vec(r * c, r * s);
        if (order >= 1) { jet[1] = vec(-r * s, r * c); }
        if (order >= 2) { jet[2] = vec(-r * c, -r * s); }
        if (order >= 3) { jet[3] = vec(r * s, -r * c); }
    }
}

describe('ReparameterizeByArclength', () => {
    it('exposes the curve, its domain and its total arclength', () => {
        const curve = new Circle2(3);
        const repar = new ReparameterizeByArclength(curve);
        expect(repar.getCurve()).toBe(curve);
        expect(repar.getTMin()).toBe(0);
        expect(repar.getTMax()).toBeCloseTo(2 * Math.PI, 14);
        expect(repar.getTotalArclength()).toBeCloseTo(6 * Math.PI, 8);
    });

    it('clamps s at or below zero to tMin', () => {
        const repar = new ReparameterizeByArclength(new Circle2(2));
        for (const s of [0, -1, -1e6]) {
            const out = repar.getT(s, true);
            expect(out.t).toBe(0);
            expect(out.f).toBe(0);
            expect(out.numIterations).toBe(0);
        }
    });

    it('clamps s at or above the total arclength to tMax', () => {
        const repar = new ReparameterizeByArclength(new Circle2(2));
        const L = repar.getTotalArclength();
        for (const s of [L, L + 1, 1e6]) {
            const out = repar.getT(s, true);
            expect(out.t).toBeCloseTo(2 * Math.PI, 14);
            expect(out.f).toBe(0);
            expect(out.numIterations).toBe(0);
        }
    });

    it('gives uniform angles for uniform arclength samples on a circle', () => {
        for (const r of [1, 2.5, 7]) {
            const repar = new ReparameterizeByArclength(new Circle2(r));
            const L = repar.getTotalArclength();
            const n = 16;
            for (const useBisection of [true, false]) {
                let previous = -1;
                for (let i = 1; i < n; ++i) {
                    const s = L * i / n;
                    const out = repar.getT(s, useBisection);
                    // Constant speed r means t = s / r.
                    expect(out.t).toBeCloseTo(s / r, 8);
                    expect(out.t).toBeGreaterThan(previous);
                    previous = out.t;
                    expect(Math.abs(out.f)).toBeLessThan(1e-8);
                }
            }
        }
    });

    it('gives equal angular increments between consecutive samples', () => {
        const repar = new ReparameterizeByArclength(new Circle2(4));
        const L = repar.getTotalArclength();
        const n = 12;
        const ts: number[] = [];
        for (let i = 0; i <= n; ++i) {
            ts.push(repar.getT(L * i / n, false).t);
        }
        const expectedStep = 2 * Math.PI / n;
        for (let i = 1; i <= n; ++i) {
            expect(ts[i] - ts[i - 1]).toBeCloseTo(expectedStep, 7);
        }
    });

    it('inverts the arclength of a non-uniform-speed Bezier curve', () => {
        const curve = new BezierCurve(2, 3,
            [vec(0, 0), vec(0.2, 3), vec(2.8, -2), vec(4, 1)]);
        const repar = new ReparameterizeByArclength(curve);
        const L = repar.getTotalArclength();
        expect(L).toBeGreaterThan(0);
        for (const useBisection of [true, false]) {
            for (let i = 1; i < 10; ++i) {
                const s = L * i / 10;
                const out = repar.getT(s, useBisection);
                expect(out.t).toBeGreaterThan(0);
                expect(out.t).toBeLessThan(1);
                // The residual F(t,s) = Arclength(tMin,t) - s is tiny.
                expect(Math.abs(curve.getLength(0, out.t) - s))
                    .toBeLessThan(1e-9);
            }
        }
    });

    it('agrees between bisection and the Newton hybrid', () => {
        const curve = new BezierCurve(3, 4,
            [vec(0, 0, 0), vec(1, 2, -1), vec(3, -1, 2), vec(4, 2, 1),
                vec(6, 0, 3)]);
        const repar = new ReparameterizeByArclength(curve);
        const L = repar.getTotalArclength();
        for (let i = 1; i < 8; ++i) {
            const s = L * i / 8;
            const a = repar.getT(s, true);
            const b = repar.getT(s, false);
            expect(a.t).toBeCloseTo(b.t, 9);
        }
    });

    it('uses far fewer iterations for the Newton hybrid than for bisection',
        () => {
            const curve = new BezierCurve(2, 3,
                [vec(0, 0), vec(1, 3), vec(3, -1), vec(4, 2)]);
            const repar = new ReparameterizeByArclength(curve);
            const s = 0.5 * repar.getTotalArclength();
            const bisect = repar.getT(s, true);
            const newton = repar.getT(s, false);
            expect(bisect.numIterations).toBeGreaterThan(10);
            expect(newton.numIterations).toBeLessThan(bisect.numIterations);
            expect(newton.numIterations).toBeGreaterThan(0);
        });

    it('is monotone in s', () => {
        const curve = new BezierCurve(2, 3,
            [vec(0, 0), vec(1, 3), vec(3, -1), vec(4, 2)]);
        const repar = new ReparameterizeByArclength(curve);
        const L = repar.getTotalArclength();
        let previous = -1;
        for (let i = 0; i <= 20; ++i) {
            const out = repar.getT(L * i / 20, false);
            expect(out.t).toBeGreaterThanOrEqual(previous);
            previous = out.t;
        }
        expect(previous).toBeCloseTo(1, 12);
    });

    it('is consistent with ParametricCurve.getTime', () => {
        const curve = new BezierCurve(2, 3,
            [vec(0, 0), vec(1, 3), vec(3, -1), vec(4, 2)]);
        const repar = new ReparameterizeByArclength(curve);
        const L = repar.getTotalArclength();
        for (let i = 1; i < 6; ++i) {
            const s = L * i / 6;
            expect(repar.getT(s, false).t).toBeCloseTo(curve.getTime(s), 6);
        }
    });
});
