import { describe, it, expect } from 'vitest';
import { NURBSCurve } from '../src/NURBSCurve.js';
import { BSplineCurve } from '../src/BSplineCurve.js';
import { BasisFunctionInput, UniqueKnot } from '../src/BasisFunction.js';
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


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). NURBSCurve.h was read line by
// line against src/NURBSCurve.ts. The file contributes the rational quotient
// and its derivative recursion, so the properties below check the position
// against a directly evaluated rational sum built on a recursive Cox-de Boor
// basis (no shared code) and the derivative recursion against numerical
// differentiation of the analytic lower-order derivative.
import {
    check, fc, expectClose, expectVectorClose, wellScaledVector
} from './helpers/arbitraries.js';

// The textbook Cox-de Boor recursion, written independently of
// BasisFunction.
function nurbsBasis(knots: readonly number[], i: number, d: number,
    t: number): number {
    if (d === 0) {
        return (knots[i] <= t && t < knots[i + 1]) ? 1 : 0;
    }
    const a = knots[i + d] - knots[i];
    const b = knots[i + d + 1] - knots[i + 1];
    let value = 0;
    if (a > 0) {
        value += ((t - knots[i]) / a) * nurbsBasis(knots, i, d - 1, t);
    }
    if (b > 0) {
        value += ((knots[i + d + 1] - t) / b) *
            nurbsBasis(knots, i + 1, d - 1, t);
    }
    return value;
}

// The rational curve evaluated straight from its definition.
function rationalPoint(curve: NURBSCurve, t: number): Vector {
    const basis = curve.getBasisFunction();
    const knots = basis.getKnots();
    const degree = basis.getDegree();
    const numControls = curve.getNumControls();
    const dim = curve.getDimension();
    const x = new Vector(dim);
    let w = 0;
    for (let i = 0; i < basis.getNumControls(); ++i) {
        const j = (i >= numControls ? i - numControls : i);
        const b = nurbsBasis(knots, i, degree, t) * curve.getWeight(j);
        for (let k = 0; k < dim; ++k) {
            x.values[k] += b * curve.getControl(j).values[k];
        }
        w += b;
    }
    for (let k = 0; k < dim; ++k) { x.values[k] /= w; }
    return x;
}

// Open uniform NURBS curves of degree 1..4 with positive weights.
const nurbsData = (dim: number): fc.Arbitrary<{
    curve: NURBSCurve, controls: Vector[], weights: number[] }> =>
    fc.integer({ min: 1, max: 4 }).chain(degree =>
        fc.integer({ min: 1, max: 5 }).chain(extra => {
            const numControls = degree + extra;
            return fc.tuple(
                fc.constant(degree),
                fc.array(wellScaledVector(dim, -6, 6),
                    { minLength: numControls, maxLength: numControls }),
                fc.array(fc.double({ min: 0.2, max: 4, noNaN: true,
                    noDefaultInfinity: true }),
                { minLength: numControls, maxLength: numControls }));
        }))
        .map(([degree, controls, weights]) => ({
            curve: new NURBSCurve(dim,
                new BasisFunctionInput(controls.length, degree),
                controls, weights),
            controls,
            weights
        }));

const inside = fc.double({ min: 0.001, max: 0.999, noNaN: true,
    noDefaultInfinity: true });

describe('NURBSCurve verification', () => {
    it('matches the rational sum of the Cox-de Boor basis', () => {
        check(fc.tuple(nurbsData(3), inside), ([{ curve }, t]) => {
            const jet = curve.createJet();
            curve.evaluate(t, 0, jet);
            expectVectorClose(jet[0], rationalPoint(curve, t), 1e-10, 1e-10);
        });
    });

    it('is unchanged when every weight is scaled by the same factor', () => {
        // The quotient is homogeneous of degree zero in the weights.
        check(fc.tuple(nurbsData(3), inside,
            fc.double({ min: 0.1, max: 8, noNaN: true,
                noDefaultInfinity: true })),
        ([{ curve, controls, weights }, t, s]) => {
            const scaled = new NURBSCurve(3,
                new BasisFunctionInput(controls.length,
                    curve.getBasisFunction().getDegree()),
                controls, weights.map(w => s * w));
            const a = curve.createJet();
            const b = scaled.createJet();
            curve.evaluate(t, 2, a);
            scaled.evaluate(t, 2, b);
            for (let order = 0; order <= 2; ++order) {
                expectVectorClose(b[order], a[order], 1e-8, 1e-8);
            }
        });
    });

    it('reduces to a B-spline curve when the weights are equal', () => {
        check(fc.tuple(nurbsData(3), inside,
            fc.double({ min: 0.5, max: 3, noNaN: true,
                noDefaultInfinity: true })),
        ([{ controls }, t, w]) => {
            const degree = Math.min(3, controls.length - 1);
            const input = new BasisFunctionInput(controls.length, degree);
            const nurbs = new NURBSCurve(3, input, controls,
                new Array<number>(controls.length).fill(w));
            const bspline = new BSplineCurve(3,
                new BasisFunctionInput(controls.length, degree), controls);
            const a = nurbs.createJet();
            const b = bspline.createJet();
            nurbs.evaluate(t, 3, a);
            bspline.evaluate(t, 3, b);
            for (let order = 0; order <= 3; ++order) {
                expectVectorClose(a[order], b[order], 1e-8, 1e-8);
            }
        });
    });

    it('derivatives agree with numerical differentiation', () => {
        // Each order is differentiated from the analytic previous order, so
        // the quotient-rule recursion in Evaluate is checked one step at a
        // time rather than through three chained finite differences.
        check(fc.tuple(nurbsData(2), fc.double({ min: 0.05, max: 0.95,
            noNaN: true, noDefaultInfinity: true })), ([{ curve }, t]) => {
            const h = 2e-3;
            // A degree-d curve is only C^{d-1}, so the derivative of order d
            // jumps at every knot; keep the whole five-point stencil inside
            // one knot span.
            for (const k of curve.getBasisFunction().getKnots()) {
                if (Math.abs(t - k) < 4 * h) { return; }
            }
            const at = (u: number, order: number): Vector => {
                const jet = curve.createJet();
                curve.evaluate(u, 3, jet);
                return jet[order];
            };
            // Five-point centered stencil (error O(step^4)) at two step
            // sizes. A rational curve can have a pole of its analytic
            // continuation just outside the knot span, which makes the high
            // derivatives enormous and the truncation error unpredictable, so
            // the tolerance is calibrated from the observed convergence
            // between the two steps instead of from a fixed constant. A wrong
            // coefficient in the quotient-rule recursion changes the value by
            // an O(1) relative amount, which no calibration can hide.
            const stencil = (step: number, order: number, k: number):
            number => {
                const p1 = at(t + step, order - 1).values[k];
                const m1 = at(t - step, order - 1).values[k];
                const p2 = at(t + 2 * step, order - 1).values[k];
                const m2 = at(t - 2 * step, order - 1).values[k];
                return (8 * (p1 - m1) - (p2 - m2)) / (12 * step);
            };
            for (let order = 1; order <= 3; ++order) {
                const analytic = at(t, order);
                for (let k = 0; k < 2; ++k) {
                    const coarse = stencil(h, order, k);
                    const fine = stencil(0.5 * h, order, k);
                    const drift = Math.abs(coarse - fine);
                    expectClose(fine, analytic.values[k], 1e-8 + drift, 1e-9);
                }
            }
        }, 60);
    });

    it('interpolates the end control points of an open knot vector', () => {
        check(nurbsData(3), ({ curve, controls }) => {
            const jet = curve.createJet();
            curve.evaluate(curve.getTMin(), 0, jet);
            expectVectorClose(jet[0], controls[0], 1e-10, 1e-10);
            curve.evaluate(curve.getTMax(), 0, jet);
            expectVectorClose(jet[0], controls[controls.length - 1],
                1e-10, 1e-10);
        });
    });

    it('stays inside the bounding box of the control points', () => {
        // With positive weights the curve is a convex combination of the
        // control points.
        check(fc.tuple(nurbsData(2), inside), ([{ curve, controls }, t]) => {
            const jet = curve.createJet();
            curve.evaluate(t, 0, jet);
            for (let k = 0; k < 2; ++k) {
                const lo = Math.min(...controls.map(c => c.values[k]));
                const hi = Math.max(...controls.map(c => c.values[k]));
                const span = Math.max(1, hi - lo);
                expect(jet[0].values[k]).toBeGreaterThanOrEqual(lo - 1e-9 * span);
                expect(jet[0].values[k]).toBeLessThanOrEqual(hi + 1e-9 * span);
            }
        });
    });

    it('copies its inputs and clamps the accessor indices', () => {
        check(nurbsData(3), ({ curve, controls, weights }) => {
            const n = controls.length;
            expect(curve.getNumControls()).toBe(n);
            for (let i = 0; i < n; ++i) {
                expect(curve.getControl(i)).not.toBe(controls[i]);
                expectVectorClose(curve.getControl(i), controls[i], 0, 0);
                expect(curve.getWeight(i)).toBe(weights[i]);
            }
            // Out-of-range accesses return element 0 and out-of-range writes
            // are ignored, exactly as upstream.
            expectVectorClose(curve.getControl(-1), curve.getControl(0), 0, 0);
            expectVectorClose(curve.getControl(n), curve.getControl(0), 0, 0);
            expect(curve.getWeight(-1)).toBe(curve.getWeight(0));
            expect(curve.getWeight(n)).toBe(curve.getWeight(0));
            const before = curve.getControl(0).clone();
            curve.setControl(-1, Vector.fromArray([9, 9, 9]));
            curve.setWeight(n, 42);
            expectVectorClose(curve.getControl(0), before, 0, 0);
            expect(curve.getWeight(0)).toBe(weights[0]);
        });
    });
});
