import { describe, it, expect } from 'vitest';
import { NURBSCurve } from '../src/NURBSCurve';
import { BSplineCurve } from '../src/BSplineCurve';
import { BasisFunctionInput, UniqueKnot } from '../src/BasisFunction';
import { Vector } from '../src/Vector';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function vec(...values: number[]): Vector {
    const v = new Vector(values.length);
    for (let i = 0; i < values.length; ++i) {
        v.values[i] = values[i];
    }
    return v;
}

function positionOf(curve: NURBSCurve | BSplineCurve, t: number): Vector {
    const jet = curve.createJet();
    curve.evaluate(t, 0, jet);
    return jet[0];
}

describe('NURBSCurve', () => {
    it('reports the domain of its basis function', () => {
        const input = new BasisFunctionInput(5, 2);
        const curve = new NURBSCurve(2, input);
        expect(curve.isConstructed()).toBe(true);
        expect(curve.getNumControls()).toBe(5);
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBe(1);
        expect(curve.getBasisFunction().getDegree()).toBe(2);
    });

    it('defaults the controls to zero and the weights to zero', () => {
        const curve = new NURBSCurve(3, new BasisFunctionInput(4, 2));
        for (let i = 0; i < 4; ++i) {
            expect(curve.getWeight(i)).toBe(0);
            expect(curve.getControl(i).values).toEqual([0, 0, 0]);
        }
    });

    it('copies controls and weights and supports set/get', () => {
        const c = [vec(0, 0), vec(1, 1), vec(2, 0), vec(3, 3)];
        const w = [1, 2, 3, 4];
        const curve = new NURBSCurve(2, new BasisFunctionInput(4, 2), c, w);
        c[1].values[0] = 100;
        w[1] = 100;
        expect(curve.getControl(1).values[0]).toBe(1);
        expect(curve.getWeight(1)).toBe(2);

        curve.setControl(2, vec(9, 9));
        curve.setWeight(2, 7);
        expect(curve.getControl(2).values).toEqual([9, 9]);
        expect(curve.getWeight(2)).toBe(7);

        // Out-of-range indices are ignored and return element 0.
        curve.setControl(-1, vec(5, 5));
        curve.setWeight(99, 5);
        expect(curve.getControl(-1).values).toEqual([0, 0]);
        expect(curve.getWeight(99)).toBe(1);
    });

    it('reduces to a B-spline curve when all weights are 1', () => {
        const rand = makeRandom(0xabc123);
        for (const degree of [1, 2, 3]) {
            const numControls = degree + 4;
            const controls: Vector[] = [];
            for (let i = 0; i < numControls; ++i) {
                controls.push(vec(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2));
            }
            const weights = new Array<number>(numControls).fill(1);
            const nurbs = new NURBSCurve(3,
                new BasisFunctionInput(numControls, degree), controls, weights);
            const bspline = new BSplineCurve(3,
                new BasisFunctionInput(numControls, degree), controls);
            for (let s = 0; s <= 20; ++s) {
                const t = s / 20;
                const a = positionOf(nurbs, t);
                const b = positionOf(bspline, t);
                for (let k = 0; k < 3; ++k) {
                    expect(a.values[k]).toBeCloseTo(b.values[k], 11);
                }
            }
        }
    });

    it('reduces to a B-spline curve when all weights are equal', () => {
        const controls = [vec(0, 0), vec(1, 2), vec(3, 1), vec(4, -1), vec(6, 2)];
        const weights = new Array<number>(5).fill(3.75);
        const nurbs = new NURBSCurve(2, new BasisFunctionInput(5, 2),
            controls, weights);
        const bspline = new BSplineCurve(2, new BasisFunctionInput(5, 2),
            controls);
        for (let s = 0; s <= 10; ++s) {
            const t = s / 10;
            const a = positionOf(nurbs, t);
            const b = positionOf(bspline, t);
            expect(a.values[0]).toBeCloseTo(b.values[0], 11);
            expect(a.values[1]).toBeCloseTo(b.values[1], 11);
        }
    });

    it('represents a quarter circle exactly (rational quadratic Bezier)', () => {
        // Controls (1,0), (1,1), (0,1) with weights 1, sqrt(2)/2, 1 on the
        // Bezier knot vector {0,0,0,1,1,1}.
        const controls = [vec(1, 0), vec(1, 1), vec(0, 1)];
        const w = Math.SQRT1_2;
        const curve = new NURBSCurve(2, new BasisFunctionInput(3, 2),
            controls, [1, w, 1]);
        for (let s = 0; s <= 20; ++s) {
            const t = s / 20;
            const p = positionOf(curve, t);
            const r = Math.hypot(p.values[0], p.values[1]);
            expect(r).toBeCloseTo(1, 12);
            // The arc stays in the first quadrant and sweeps monotonically.
            expect(p.values[0]).toBeGreaterThanOrEqual(-1e-12);
            expect(p.values[1]).toBeGreaterThanOrEqual(-1e-12);
        }
        const p0 = positionOf(curve, 0);
        const p1 = positionOf(curve, 1);
        expect(p0.values[0]).toBeCloseTo(1, 12);
        expect(p0.values[1]).toBeCloseTo(0, 12);
        expect(p1.values[0]).toBeCloseTo(0, 12);
        expect(p1.values[1]).toBeCloseTo(1, 12);
    });

    it('sweeps the quarter circle with monotone angle', () => {
        const controls = [vec(1, 0), vec(1, 1), vec(0, 1)];
        const w = Math.SQRT1_2;
        const curve = new NURBSCurve(2, new BasisFunctionInput(3, 2),
            controls, [1, w, 1]);
        let previous = -1;
        for (let s = 0; s <= 20; ++s) {
            const p = positionOf(curve, s / 20);
            const angle = Math.atan2(p.values[1], p.values[0]);
            expect(angle).toBeGreaterThan(previous);
            previous = angle;
        }
        expect(previous).toBeCloseTo(Math.PI / 2, 12);
    });

    it('has the arclength of the quarter circle', () => {
        const controls = [vec(1, 0), vec(1, 1), vec(0, 1)];
        const curve = new NURBSCurve(2, new BasisFunctionInput(3, 2),
            controls, [1, Math.SQRT1_2, 1]);
        expect(curve.getTotalLength()).toBeCloseTo(Math.PI / 2, 8);
    });

    it('computes derivatives that match finite differences', () => {
        const controls = [vec(1, 0), vec(1, 1), vec(0, 1)];
        const curve = new NURBSCurve(2, new BasisFunctionInput(3, 2),
            controls, [1, Math.SQRT1_2, 1]);
        const jet = curve.createJet();
        const h = 1e-5;
        for (const t of [0.2, 0.4, 0.6, 0.8]) {
            curve.evaluate(t, 3, jet);
            const p = positionOf(curve, t);
            const pp = positionOf(curve, t + h);
            const pm = positionOf(curve, t - h);
            for (let k = 0; k < 2; ++k) {
                const d1 = (pp.values[k] - pm.values[k]) / (2 * h);
                const d2 = (pp.values[k] - 2 * p.values[k] + pm.values[k])
                    / (h * h);
                expect(jet[1].values[k]).toBeCloseTo(d1, 6);
                expect(jet[2].values[k]).toBeCloseTo(d2, 3);
            }
        }
    });

    it('has the tangent of the circle (perpendicular to the radius)', () => {
        const controls = [vec(1, 0), vec(1, 1), vec(0, 1)];
        const curve = new NURBSCurve(2, new BasisFunctionInput(3, 2),
            controls, [1, Math.SQRT1_2, 1]);
        const jet = curve.createJet();
        for (const t of [0.1, 0.5, 0.9]) {
            curve.evaluate(t, 1, jet);
            const dotRT = jet[0].values[0] * jet[1].values[0]
                + jet[0].values[1] * jet[1].values[1];
            expect(dotRT).toBeCloseTo(0, 10);
        }
    });

    it('supports periodic splines by index wrapping', () => {
        const input = new BasisFunctionInput(4, 2);
        input.periodic = true;
        input.numUniqueKnots = 4 + 2 * 2 - 1;
        input.uniqueKnots = [];
        const n = input.numUniqueKnots;
        for (let i = 0; i < n; ++i) {
            input.uniqueKnots.push(new UniqueKnot(i / (n - 1), 1));
        }
        const controls = [vec(1, 0), vec(0, 1), vec(-1, 0), vec(0, -1)];
        const curve = new NURBSCurve(2, input, controls, [1, 1, 1, 1]);
        expect(curve.isConstructed()).toBe(true);
        // A closed curve: X(tmin) equals X(tmax).
        const a = positionOf(curve, curve.getTMin());
        const b = positionOf(curve, curve.getTMax());
        expect(a.values[0]).toBeCloseTo(b.values[0], 10);
        expect(a.values[1]).toBeCloseTo(b.values[1], 10);
    });

    it('zeroes the jet for order >= SUP_ORDER', () => {
        const controls = [vec(1, 0), vec(1, 1), vec(0, 1)];
        const curve = new NURBSCurve(2, new BasisFunctionInput(3, 2),
            controls, [1, 1, 1]);
        const jet = curve.createJet();
        curve.evaluate(0.5, 4, jet);
        for (let i = 0; i < 4; ++i) {
            expect(jet[i].values).toEqual([0, 0]);
        }
    });
});
