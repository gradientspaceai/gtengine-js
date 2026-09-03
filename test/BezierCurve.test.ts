import { describe, it, expect } from 'vitest';
import { BezierCurve } from '../src/BezierCurve.js';
import { Vector } from '../src/Vector.js';

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

// The de Casteljau evaluation of a Bezier curve, an independent computation.
function deCasteljau(controls: readonly Vector[], t: number): Vector {
    let level = controls.map(c => c.clone());
    while (level.length > 1) {
        const next: Vector[] = [];
        for (let i = 0; i + 1 < level.length; ++i) {
            const v = new Vector(level[i].size);
            for (let k = 0; k < v.size; ++k) {
                v.values[k] = (1 - t) * level[i].values[k]
                    + t * level[i + 1].values[k];
            }
            next.push(v);
        }
        level = next;
    }
    return level[0];
}

function positionOf(curve: BezierCurve, t: number): Vector {
    const jet = curve.createJet();
    curve.evaluate(t, 0, jet);
    return jet[0];
}

describe('BezierCurve', () => {
    it('rejects degree less than 2', () => {
        expect(() => new BezierCurve(2, 1, [vec(0, 0), vec(1, 1)]))
            .toThrow(/Invalid input/);
    });

    it('reports its degree, control count and domain', () => {
        const controls = [vec(0, 0), vec(1, 2), vec(3, 1)];
        const curve = new BezierCurve(2, 2, controls);
        expect(curve.isConstructed()).toBe(true);
        expect(curve.getDegree()).toBe(2);
        expect(curve.getNumControls()).toBe(3);
        expect(curve.getDimension()).toBe(2);
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBe(1);
        expect(curve.getControls().length).toBe(3);
    });

    it('copies the control points (C++ value semantics)', () => {
        const p1 = vec(1, 2);
        const controls = [vec(0, 0), p1, vec(3, 1)];
        const curve = new BezierCurve(2, 2, controls);
        p1.values[0] = 100;
        expect(curve.getControls()[1].values[0]).toBe(1);
    });

    it('interpolates the endpoints of a quadratic', () => {
        const controls = [vec(0, 0), vec(1, 2), vec(3, 1)];
        const curve = new BezierCurve(2, 2, controls);
        const p0 = positionOf(curve, 0);
        const p1 = positionOf(curve, 1);
        expect(p0.values[0]).toBeCloseTo(0, 12);
        expect(p0.values[1]).toBeCloseTo(0, 12);
        expect(p1.values[0]).toBeCloseTo(3, 12);
        expect(p1.values[1]).toBeCloseTo(1, 12);
    });

    it('interpolates the endpoints of a cubic', () => {
        const controls = [vec(-1, 0, 2), vec(1, 2, 0), vec(3, 1, -1),
            vec(4, -2, 5)];
        const curve = new BezierCurve(3, 3, controls);
        const p0 = positionOf(curve, 0);
        const p1 = positionOf(curve, 1);
        for (let k = 0; k < 3; ++k) {
            expect(p0.values[k]).toBeCloseTo(controls[0].values[k], 12);
            expect(p1.values[k]).toBeCloseTo(controls[3].values[k], 12);
        }
    });

    it('matches the Bernstein form for a quadratic', () => {
        const controls = [vec(0, 0), vec(1, 2), vec(3, 1)];
        const curve = new BezierCurve(2, 2, controls);
        for (const t of [0, 0.25, 0.5, 0.75, 1]) {
            const omt = 1 - t;
            const expected = [0, 0];
            for (let k = 0; k < 2; ++k) {
                expected[k] = omt * omt * controls[0].values[k]
                    + 2 * t * omt * controls[1].values[k]
                    + t * t * controls[2].values[k];
            }
            const p = positionOf(curve, t);
            expect(p.values[0]).toBeCloseTo(expected[0], 12);
            expect(p.values[1]).toBeCloseTo(expected[1], 12);
        }
    });

    it('agrees with de Casteljau for degrees 2 through 6', () => {
        const rand = makeRandom(0x51ced1);
        for (let degree = 2; degree <= 6; ++degree) {
            const controls: Vector[] = [];
            for (let i = 0; i <= degree; ++i) {
                controls.push(vec(4 * rand() - 2, 4 * rand() - 2,
                    4 * rand() - 2));
            }
            const curve = new BezierCurve(3, degree, controls);
            for (let s = 0; s <= 10; ++s) {
                const t = s / 10;
                const actual = positionOf(curve, t);
                const expected = deCasteljau(controls, t);
                for (let k = 0; k < 3; ++k) {
                    expect(actual.values[k]).toBeCloseTo(expected.values[k], 11);
                }
            }
        }
    });

    it('computes derivatives that match finite differences', () => {
        const rand = makeRandom(0xbeef01);
        for (let degree = 2; degree <= 5; ++degree) {
            const controls: Vector[] = [];
            for (let i = 0; i <= degree; ++i) {
                controls.push(vec(4 * rand() - 2, 4 * rand() - 2));
            }
            const curve = new BezierCurve(2, degree, controls);
            const jet = curve.createJet();
            const h = 1e-5;
            for (const t of [0.2, 0.4, 0.6, 0.8]) {
                curve.evaluate(t, 3, jet);
                const pPlus = deCasteljau(controls, t + h);
                const pMinus = deCasteljau(controls, t - h);
                const p = deCasteljau(controls, t);
                for (let k = 0; k < 2; ++k) {
                    const d1 = (pPlus.values[k] - pMinus.values[k]) / (2 * h);
                    const d2 = (pPlus.values[k] - 2 * p.values[k]
                        + pMinus.values[k]) / (h * h);
                    expect(jet[1].values[k]).toBeCloseTo(d1, 6);
                    expect(jet[2].values[k]).toBeCloseTo(d2, 3);
                }
            }
        }
    });

    it('has an exact first derivative at the endpoints', () => {
        // X'(0) = d*(C1 - C0) and X'(1) = d*(Cd - C{d-1}).
        const controls = [vec(-1, 0), vec(1, 2), vec(3, 1), vec(4, -2)];
        const curve = new BezierCurve(2, 3, controls);
        const jet = curve.createJet();
        curve.evaluate(0, 1, jet);
        expect(jet[1].values[0]).toBeCloseTo(3 * (1 - (-1)), 12);
        expect(jet[1].values[1]).toBeCloseTo(3 * (2 - 0), 12);
        curve.evaluate(1, 1, jet);
        expect(jet[1].values[0]).toBeCloseTo(3 * (4 - 3), 12);
        expect(jet[1].values[1]).toBeCloseTo(3 * (-2 - 1), 12);
    });

    it('has a constant second derivative for a quadratic', () => {
        const controls = [vec(0, 0), vec(1, 2), vec(3, 1)];
        const curve = new BezierCurve(2, 2, controls);
        const jet = curve.createJet();
        for (const t of [0, 0.3, 0.7, 1]) {
            curve.evaluate(t, 2, jet);
            expect(jet[2].values[0]).toBeCloseTo(2 * (0 - 2 * 1 + 3), 12);
            expect(jet[2].values[1]).toBeCloseTo(2 * (0 - 2 * 2 + 1), 12);
        }
    });

    it('returns a zero third derivative for a quadratic', () => {
        const controls = [vec(0, 0), vec(1, 2), vec(3, 1)];
        const curve = new BezierCurve(2, 2, controls);
        const jet = curve.createJet();
        curve.evaluate(0.5, 3, jet);
        expect(jet[3].values[0]).toBe(0);
        expect(jet[3].values[1]).toBe(0);
    });

    it('has a constant third derivative for a cubic', () => {
        const c = [vec(-1, 0), vec(1, 2), vec(3, 1), vec(4, -2)];
        const curve = new BezierCurve(2, 3, c);
        const jet = curve.createJet();
        for (const t of [0, 0.3, 0.9]) {
            curve.evaluate(t, 3, jet);
            for (let k = 0; k < 2; ++k) {
                const expected = 6 * (-c[0].values[k] + 3 * c[1].values[k]
                    - 3 * c[2].values[k] + c[3].values[k]);
                expect(jet[3].values[k]).toBeCloseTo(expected, 11);
            }
        }
    });

    it('reproduces a straight line when the controls are collinear', () => {
        const controls = [vec(0, 0), vec(1, 1), vec(2, 2), vec(3, 3)];
        const curve = new BezierCurve(2, 3, controls);
        for (let s = 0; s <= 10; ++s) {
            const t = s / 10;
            const p = positionOf(curve, t);
            expect(p.values[0]).toBeCloseTo(3 * t, 11);
            expect(p.values[1]).toBeCloseTo(3 * t, 11);
        }
    });

    it('zeroes the jet for order >= SUP_ORDER', () => {
        const controls = [vec(0, 0), vec(1, 2), vec(3, 1)];
        const curve = new BezierCurve(2, 2, controls);
        const jet = curve.createJet();
        curve.evaluate(0.5, 4, jet);
        for (let i = 0; i < 4; ++i) {
            expect(jet[i].values[0]).toBe(0);
            expect(jet[i].values[1]).toBe(0);
        }
    });

    it('computes the arclength of a degenerate (collinear) cubic', () => {
        const controls = [vec(0, 0), vec(1, 0), vec(2, 0), vec(3, 0)];
        const curve = new BezierCurve(2, 3, controls);
        expect(curve.getTotalLength()).toBeCloseTo(3, 8);
    });
});
