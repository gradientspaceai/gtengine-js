import { describe, it, expect } from 'vitest';
import { BSplineCurve } from '../src/BSplineCurve.js';
import { BasisFunctionInput, UniqueKnot } from '../src/BasisFunction.js';
import { ParametricCurve } from '../src/ParametricCurve.js';
import { Vector } from '../src/Vector.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

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

// The de Boor evaluation of an open uniform B-spline via the Cox-de Boor
// basis functions, computed independently of BasisFunction.
function coxDeBoor(knot: number[], i: number, j: number, t: number): number {
    if (j > 0) {
        let result = 0;
        let denom = knot[i + j] - knot[i];
        if (denom > 0) {
            result += (t - knot[i]) * coxDeBoor(knot, i, j - 1, t) / denom;
        }
        denom = knot[i + j + 1] - knot[i + 1];
        if (denom > 0) {
            result += (knot[i + j + 1] - t) * coxDeBoor(knot, i + 1, j - 1, t) / denom;
        }
        return result;
    }
    return (knot[i] <= t && t < knot[i + 1]) ? 1 : 0;
}

// The open-uniform knot vector on [0,1] used by BasisFunctionInput(n, d).
function openUniformKnots(numControls: number, degree: number): number[] {
    const numKnots = numControls + degree + 1;
    const knot = new Array<number>(numKnots).fill(0);
    let i: number;
    for (i = 0; i <= degree; ++i) {
        knot[i] = 0;
    }
    const factor = 1 / (numControls - degree);
    for (/**/; i < numControls; ++i) {
        knot[i] = (i - degree) * factor;
    }
    for (/**/; i < numKnots; ++i) {
        knot[i] = 1;
    }
    return knot;
}

function referencePosition(controls: Vector[], degree: number, t: number): Vector {
    const knot = openUniformKnots(controls.length, degree);
    const dimension = controls[0].size;
    const result = new Vector(dimension);
    for (let i = 0; i < controls.length; ++i) {
        const n = coxDeBoor(knot, i, degree, t);
        for (let k = 0; k < dimension; ++k) {
            result.values[k] += n * controls[i].values[k];
        }
    }
    return result;
}

describe('BSplineCurve: construction', () => {
    it('takes its domain from the basis function', () => {
        const input = new BasisFunctionInput(5, 2);
        const curve = new BSplineCurve(3, input);
        expect(curve.isConstructed()).toBe(true);
        expect(curve.getDimension()).toBe(3);
        expect(curve.getNumControls()).toBe(5);
        expect(curve.getBasisFunction().getDegree()).toBe(2);
        expect(curve.getTMin()).toBe(curve.getBasisFunction().getMinDomain());
        expect(curve.getTMax()).toBe(curve.getBasisFunction().getMaxDomain());
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBe(1);
    });

    it('zero-fills deferred control points', () => {
        const curve = new BSplineCurve(2, new BasisFunctionInput(4, 1));
        for (let i = 0; i < 4; ++i) {
            expect(curve.getControl(i).values).toEqual([0, 0]);
        }
    });

    it('copies the supplied control points', () => {
        const controls = [vec(0, 0), vec(1, 0), vec(1, 1)];
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 2), controls);
        controls[2].values[1] = -100;
        expect(curve.getControl(2).values).toEqual([1, 1]);
    });

    it('setControl copies and ignores out-of-range indices', () => {
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 1));
        const p = vec(2, 3);
        curve.setControl(1, p);
        p.values[0] = 99;
        expect(curve.getControl(1).values).toEqual([2, 3]);

        // Out-of-range writes are silently ignored, as upstream.
        curve.setControl(-1, vec(7, 7));
        curve.setControl(3, vec(7, 7));
        expect(curve.getControl(0).values).toEqual([0, 0]);
        expect(curve.getControl(2).values).toEqual([0, 0]);
    });

    it('getControl of an out-of-range index returns control 0', () => {
        const controls = [vec(1, 2), vec(3, 4), vec(5, 6)];
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 1), controls);
        expect(curve.getControl(-1).values).toEqual([1, 2]);
        expect(curve.getControl(3).values).toEqual([1, 2]);
    });

    it('getControls aliases the internal storage', () => {
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 1));
        curve.getControls()[0].values[0] = 42;
        expect(curve.getControl(0).values).toEqual([42, 0]);
    });
});

describe('BSplineCurve: evaluation', () => {
    it('interpolates the endpoints of an open curve', () => {
        const controls = [vec(0, 0), vec(1, 2), vec(3, -1), vec(4, 4)];
        const curve = new BSplineCurve(2, new BasisFunctionInput(4, 3), controls);
        expect(curve.getPosition(0).values[0]).toBeCloseTo(0, 12);
        expect(curve.getPosition(0).values[1]).toBeCloseTo(0, 12);
        expect(curve.getPosition(1).values[0]).toBeCloseTo(4, 12);
        expect(curve.getPosition(1).values[1]).toBeCloseTo(4, 12);
    });

    it('degree 1 is the polyline through the control points', () => {
        const controls = [vec(0, 0), vec(2, 1), vec(3, 5)];
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 1), controls);
        // With three control points and degree 1, the knots are
        // {0,0,0.5,1,1}, so the curve is the polyline.
        for (const [t, x, y] of [[0, 0, 0], [0.25, 1, 0.5], [0.5, 2, 1],
            [0.75, 2.5, 3], [1, 3, 5]] as Array<[number, number, number]>) {
            const p = curve.getPosition(t);
            expect(p.values[0]).toBeCloseTo(x, 12);
            expect(p.values[1]).toBeCloseTo(y, 12);
        }
    });

    it('a degree-2 curve with 3 controls is the quadratic Bezier curve', () => {
        const p0 = vec(0, 0), p1 = vec(1, 2), p2 = vec(2, 0);
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 2), [p0, p1, p2]);
        for (let k = 0; k <= 10; ++k) {
            const t = k / 10;
            const s = 1 - t;
            const bx = s * s * 0 + 2 * s * t * 1 + t * t * 2;
            const by = s * s * 0 + 2 * s * t * 2 + t * t * 0;
            const p = curve.getPosition(t);
            expect(p.values[0]).toBeCloseTo(bx, 12);
            expect(p.values[1]).toBeCloseTo(by, 12);

            // The first derivative of the Bezier curve.
            const jet = curve.createJet();
            curve.evaluate(t, 1, jet);
            expect(jet[1].values[0]).toBeCloseTo(2 * s * 1 + 2 * t * 1, 10);
            expect(jet[1].values[1]).toBeCloseTo(2 * s * 2 - 2 * t * 2, 10);
        }
    });

    it('partitions unity: an affine combination of equal controls is fixed', () => {
        // If all control points are identical, the curve is the constant.
        const c = vec(5, -3, 7);
        const curve = new BSplineCurve(3, new BasisFunctionInput(6, 3),
            [c, c.clone(), c.clone(), c.clone(), c.clone(), c.clone()]);
        for (let k = 0; k <= 8; ++k) {
            const p = curve.getPosition(k / 8);
            expect(p.values[0]).toBeCloseTo(5, 10);
            expect(p.values[1]).toBeCloseTo(-3, 10);
            expect(p.values[2]).toBeCloseTo(7, 10);
        }
    });

    it('cross-checks positions against an independent Cox-de Boor evaluation', () => {
        const rand = makeRandom(20260901);
        for (const [numControls, degree] of [[4, 1], [5, 2], [6, 3], [8, 2]]) {
            const controls: Vector[] = [];
            for (let i = 0; i < numControls; ++i) {
                controls.push(vec(rand() * 10 - 5, rand() * 10 - 5, rand() * 10 - 5));
            }
            const curve = new BSplineCurve(3,
                new BasisFunctionInput(numControls, degree), controls);
            for (let k = 0; k < 20; ++k) {
                // Avoid t = 1 where the half-open Cox-de Boor support makes
                // the reference return the zero vector.
                const t = 0.999 * (k / 20);
                const actual = curve.getPosition(t);
                const expected = referencePosition(controls, degree, t);
                for (let d = 0; d < 3; ++d) {
                    expect(actual.values[d]).toBeCloseTo(expected.values[d], 10);
                }
            }
        }
    });

    it('derivatives agree with central differences', () => {
        const rand = makeRandom(4242);
        const controls: Vector[] = [];
        for (let i = 0; i < 7; ++i) {
            controls.push(vec(rand() * 4 - 2, rand() * 4 - 2));
        }
        const curve = new BSplineCurve(2, new BasisFunctionInput(7, 3), controls);
        const h = 1e-5;
        const jet = curve.createJet();
        for (const t of [0.13, 0.27, 0.41, 0.58, 0.72, 0.86]) {
            curve.evaluate(t, 3, jet);
            const pPlus = curve.getPosition(t + h);
            const pMinus = curve.getPosition(t - h);
            const p = curve.getPosition(t);
            for (let d = 0; d < 2; ++d) {
                const d1 = (pPlus.values[d] - pMinus.values[d]) / (2 * h);
                const d2 = (pPlus.values[d] - 2 * p.values[d] + pMinus.values[d]) / (h * h);
                expect(jet[1].values[d]).toBeCloseTo(d1, 6);
                expect(jet[2].values[d]).toBeCloseTo(d2, 3);
            }
        }
    });

    it('the third derivative of a cubic B-spline is piecewise constant', () => {
        const controls = [vec(0, 0), vec(1, 3), vec(2, -1), vec(3, 2), vec(4, 0)];
        const curve = new BSplineCurve(2, new BasisFunctionInput(5, 3), controls);
        const jet = curve.createJet();
        curve.evaluate(0.1, 3, jet);
        const first = [jet[3].values[0], jet[3].values[1]];
        curve.evaluate(0.2, 3, jet);
        expect(jet[3].values[0]).toBeCloseTo(first[0], 10);
        expect(jet[3].values[1]).toBeCloseTo(first[1], 10);
    });

    it('returns a zero jet when the order is at least SUP_ORDER', () => {
        const controls = [vec(1, 1), vec(2, 2), vec(3, 4)];
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 2), controls);
        const jet = curve.createJet();
        for (let i = 0; i < ParametricCurve.SUP_ORDER; ++i) {
            jet[i].values[0] = 123;
        }
        curve.evaluate(0.5, ParametricCurve.SUP_ORDER, jet);
        for (let i = 0; i < ParametricCurve.SUP_ORDER; ++i) {
            expect(jet[i].values).toEqual([0, 0]);
        }
    });

    it('fills only the requested orders', () => {
        const controls = [vec(0, 0), vec(1, 1), vec(2, 0)];
        const curve = new BSplineCurve(2, new BasisFunctionInput(3, 2), controls);
        const jet = curve.createJet();
        const sentinel = jet[1];
        curve.evaluate(0.5, 0, jet);
        // Only jet[0] was replaced; jet[1] is still the zero vector supplied.
        expect(jet[1]).toBe(sentinel);
        expect(jet[1].values).toEqual([0, 0]);
    });
});

describe('BSplineCurve: periodic curves', () => {
    it('wraps the control-point index for a periodic curve', () => {
        const input = new BasisFunctionInput(4, 2);
        input.periodic = true;
        const controls = [vec(1, 0), vec(0, 1), vec(-1, 0), vec(0, -1)];
        const curve = new BSplineCurve(2, input, controls);
        expect(curve.getNumControls()).toBe(4);
        // The basis function stores numControls + degree control indices; the
        // curve wraps them, so evaluation succeeds over the whole domain.
        const tmin = curve.getTMin();
        const tmax = curve.getTMax();
        for (let k = 0; k <= 8; ++k) {
            const t = tmin + (tmax - tmin) * (k / 8);
            const p = curve.getPosition(t);
            expect(Number.isFinite(p.values[0])).toBe(true);
            expect(Number.isFinite(p.values[1])).toBe(true);
            // The curve is inside the convex hull of the control points.
            expect(Math.abs(p.values[0])).toBeLessThanOrEqual(1 + 1e-12);
            expect(Math.abs(p.values[1])).toBeLessThanOrEqual(1 + 1e-12);
        }
    });

    it('a periodic curve of identical controls is constant', () => {
        const input = new BasisFunctionInput(5, 2);
        input.periodic = true;
        const c = vec(2, -4);
        const curve = new BSplineCurve(2, input,
            [c, c.clone(), c.clone(), c.clone(), c.clone()]);
        for (let k = 0; k <= 8; ++k) {
            const t = curve.getTMin()
                + (curve.getTMax() - curve.getTMin()) * (k / 8);
            const p = curve.getPosition(t);
            expect(p.values[0]).toBeCloseTo(2, 10);
            expect(p.values[1]).toBeCloseTo(-4, 10);
        }
    });
});

describe('BSplineCurve: non-uniform knots and inherited quantities', () => {
    it('accepts a non-uniform knot vector', () => {
        const input = new BasisFunctionInput();
        input.numControls = 4;
        input.degree = 2;
        input.uniform = false;
        input.periodic = false;
        input.numUniqueKnots = 3;
        input.uniqueKnots = [new UniqueKnot(), new UniqueKnot(), new UniqueKnot()];
        input.uniqueKnots[0].t = 0;
        input.uniqueKnots[0].multiplicity = 3;
        input.uniqueKnots[1].t = 0.25;
        input.uniqueKnots[1].multiplicity = 1;
        input.uniqueKnots[2].t = 1;
        input.uniqueKnots[2].multiplicity = 3;

        const controls = [vec(0, 0), vec(1, 1), vec(2, -1), vec(3, 0)];
        const curve = new BSplineCurve(2, input, controls);
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBe(1);
        expect(curve.getPosition(0).values[0]).toBeCloseTo(0, 12);
        expect(curve.getPosition(1).values[0]).toBeCloseTo(3, 12);
    });

    it('computes the length of a straight-line curve', () => {
        // Degree 1 with two controls is the segment from p0 to p1.
        const curve = new BSplineCurve(2, new BasisFunctionInput(2, 1),
            [vec(0, 0), vec(3, 4)]);
        expect(curve.getTotalLength()).toBeCloseTo(5, 8);
        expect(curve.getLength(0, 0.5)).toBeCloseTo(2.5, 8);
        const tangent = curve.getTangent(0.5);
        expect(tangent.values[0]).toBeCloseTo(0.6, 10);
        expect(tangent.values[1]).toBeCloseTo(0.8, 10);
        expect(curve.getSpeed(0.5)).toBeCloseTo(5, 10);
    });

    it('getTime inverts getLength on a straight-line curve', () => {
        const curve = new BSplineCurve(2, new BasisFunctionInput(2, 1),
            [vec(0, 0), vec(10, 0)]);
        expect(curve.getTime(2.5)).toBeCloseTo(0.25, 6);
        expect(curve.getTime(0)).toBeCloseTo(0, 8);
        // Lengths beyond the total clamp to the maximum time.
        expect(curve.getTime(100)).toBe(1);
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). Positions and derivatives are
// cross-checked against an independent Cox-de Boor evaluation: the derivative
// of a degree-d B-spline is the degree-(d-1) B-spline on the trimmed knot
// vector whose control points are d*(C[i]-C[i-1])/(t[i+d]-t[i]), which shares
// no code with BasisFunction::Evaluate.
// ---------------------------------------------------------------------------

describe('BSplineCurve verification', () => {
    // The control points and knots of the order-th derivative curve.
    function derivativeData(controls: Vector[], degree: number, order: number):
        { controls: Vector[], knot: number[], degree: number } {
        const dimension = controls[0].size;
        let knot = openUniformKnots(controls.length, degree);
        let ctrl = controls.map(c => c.clone());
        let deg = degree;
        for (let r = 0; r < order; ++r) {
            const next: Vector[] = [];
            for (let i = 1; i < ctrl.length; ++i) {
                const denom = knot[i + deg] - knot[i];
                const q = new Vector(dimension);
                if (denom > 0) {
                    for (let k = 0; k < dimension; ++k) {
                        q.values[k] = deg * (ctrl[i].values[k]
                            - ctrl[i - 1].values[k]) / denom;
                    }
                }
                next.push(q);
            }
            ctrl = next;
            knot = knot.slice(1, knot.length - 1);
            deg = deg - 1;
        }
        return { controls: ctrl, knot, degree: deg };
    }

    function referenceJet(controls: Vector[], degree: number, order: number,
        t: number): Vector {
        const data = derivativeData(controls, degree, order);
        const dimension = controls[0].size;
        const result = new Vector(dimension);
        if (data.degree < 0) { return result; }
        for (let i = 0; i < data.controls.length; ++i) {
            const n = coxDeBoor(data.knot, i, data.degree, t);
            for (let k = 0; k < dimension; ++k) {
                result.values[k] += n * data.controls[i].values[k];
            }
        }
        return result;
    }

    interface Spec {
        numControls: number;
        degree: number;
        dimension: number;
        controls: Vector[];
    }

    const spec: fc.Arbitrary<Spec> =
        fc.tuple(fc.integer({ min: 2, max: 9 }), fc.integer({ min: 1, max: 3 }),
            fc.integer({ min: 1, max: 3 }), fc.integer({ min: 1, max: 1 << 20 }))
            .map(([numControls, degreeRaw, dimension, seed]) => {
                const degree = Math.min(degreeRaw, numControls - 1);
                const rand = makeRandom(seed);
                const controls: Vector[] = [];
                for (let i = 0; i < numControls; ++i) {
                    const values: number[] = [];
                    for (let d = 0; d < dimension; ++d) {
                        values.push(rand() * 8 - 4);
                    }
                    controls.push(vec(...values));
                }
                return { numControls, degree, dimension, controls };
            });

    it('matches an independent Cox-de Boor evaluation of the whole jet', () => {
        check(fc.tuple(spec, fc.integer({ min: 0, max: 999 })),
            ([s, tIndex]) => {
                const curve = new BSplineCurve(s.dimension,
                    new BasisFunctionInput(s.numControls, s.degree), s.controls);
                // Avoid t = 1, where the half-open support of the reference
                // Cox-de Boor recursion returns zero.
                const t = tIndex / 1000;
                // At an interior knot j/(numControls - degree) the degree-th
                // derivative is discontinuous and the reference recursion's
                // half-open convention need not match the curve's; skip.
                const m = s.numControls - s.degree;
                if (Math.abs(t * m - Math.round(t * m)) < 1e-9) {
                    return true;
                }
                const jet = curve.createJet();
                const order = Math.min(3, s.degree);
                curve.evaluate(t, order, jet);
                for (let k = 0; k <= order; ++k) {
                    const expected = referenceJet(s.controls, s.degree, k, t);
                    for (let d = 0; d < s.dimension; ++d) {
                        expectClose(jet[k].values[d], expected.values[d],
                            1e-9, 1e-9);
                    }
                }
            });
    }, 30000);

    it('fills exactly the requested jet slots and no others', () => {
        check(fc.tuple(spec, fc.integer({ min: 0, max: 4 })), ([s, order]) => {
            const curve = new BSplineCurve(s.dimension,
                new BasisFunctionInput(s.numControls, s.degree), s.controls);
            const jet = curve.createJet();
            const sentinels = jet.slice();
            curve.evaluate(0.5, order, jet);
            if (order >= ParametricCurve.SUP_ORDER) {
                // Invalid order: every slot is zeroed in place.
                for (let i = 0; i < ParametricCurve.SUP_ORDER; ++i) {
                    expect(jet[i]).toBe(sentinels[i]);
                    for (let d = 0; d < s.dimension; ++d) {
                        expect(jet[i].values[d]).toBe(0);
                    }
                }
                return;
            }
            for (let i = 0; i < ParametricCurve.SUP_ORDER; ++i) {
                if (i <= order) {
                    // Replaced by a freshly computed vector.
                    expect(jet[i]).not.toBe(sentinels[i]);
                } else {
                    // Untouched, as upstream leaves the higher slots alone.
                    expect(jet[i]).toBe(sentinels[i]);
                }
            }
        });
    });

    it('is a partition of unity and translation equivariant', () => {
        check(fc.tuple(spec, fc.integer({ min: 0, max: 999 }),
            fc.double({ min: -5, max: 5, noNaN: true })),
            ([s, tIndex, shift]) => {
                const t = tIndex / 1000;
                const curve = new BSplineCurve(s.dimension,
                    new BasisFunctionInput(s.numControls, s.degree), s.controls);
                const shifted = new BSplineCurve(s.dimension,
                    new BasisFunctionInput(s.numControls, s.degree),
                    s.controls.map(c => {
                        const q = c.clone();
                        for (let d = 0; d < s.dimension; ++d) {
                            q.values[d] += shift;
                        }
                        return q;
                    }));
                const p = curve.getPosition(t);
                const q = shifted.getPosition(t);
                for (let d = 0; d < s.dimension; ++d) {
                    expectClose(q.values[d], p.values[d] + shift, 1e-9, 1e-9);
                }

                // Identical control points give the constant curve, which is
                // the partition-of-unity statement.
                const constant = new BSplineCurve(s.dimension,
                    new BasisFunctionInput(s.numControls, s.degree),
                    s.controls.map(() => s.controls[0]));
                const c = constant.getPosition(t);
                for (let d = 0; d < s.dimension; ++d) {
                    expectClose(c.values[d], s.controls[0].values[d],
                        1e-9, 1e-9);
                }
            });
    });

    it('copies the control points in and shares them out, as upstream does',
        () => {
            check(fc.tuple(spec, fc.integer({ min: -3, max: 12 })),
                ([s, index]) => {
                    const curve = new BSplineCurve(s.dimension,
                        new BasisFunctionInput(s.numControls, s.degree),
                        s.controls);
                    // The constructor copied the inputs.
                    const saved = s.controls[0].values.slice();
                    s.controls[0].values[0] += 1000;
                    expect(curve.getControl(0).values[0]).toBe(saved[0]);
                    s.controls[0].values[0] = saved[0];

                    // getControl clamps out-of-range indices to control 0.
                    const control = curve.getControl(index);
                    if (index < 0 || index >= s.numControls) {
                        expect(control).toBe(curve.getControl(0));
                    } else {
                        expect(control).toBe(curve.getControls()[index]);
                    }

                    // setControl copies and ignores out-of-range indices.
                    const replacement = new Vector(s.dimension);
                    replacement.values[0] = 42;
                    curve.setControl(index, replacement);
                    if (0 <= index && index < s.numControls) {
                        expect(curve.getControl(index).values[0]).toBe(42);
                        expect(curve.getControl(index)).not.toBe(replacement);
                    } else {
                        expect(curve.getNumControls()).toBe(s.numControls);
                    }
                });
        });

    it('takes its domain from the basis function', () => {
        check(spec, s => {
            const input = new BasisFunctionInput(s.numControls, s.degree);
            const curve = new BSplineCurve(s.dimension, input, s.controls);
            expect(curve.getTMin()).toBe(curve.getBasisFunction().getMinDomain());
            expect(curve.getTMax()).toBe(curve.getBasisFunction().getMaxDomain());
            expect(curve.isConstructed()).toBe(true);
            expect(curve.getNumControls()).toBe(s.numControls);
            // An open curve interpolates its end control points.
            const start = curve.getPosition(curve.getTMin());
            const end = curve.getPosition(curve.getTMax());
            for (let d = 0; d < s.dimension; ++d) {
                expectClose(start.values[d], s.controls[0].values[d], 1e-9, 1e-9);
                expectClose(end.values[d],
                    s.controls[s.numControls - 1].values[d], 1e-9, 1e-9);
            }
        });
    });

    it('zero-fills deferred control points', () => {
        check(spec, s => {
            const curve = new BSplineCurve(s.dimension,
                new BasisFunctionInput(s.numControls, s.degree));
            expect(curve.getNumControls()).toBe(s.numControls);
            for (let i = 0; i < s.numControls; ++i) {
                expect(curve.getControl(i).size).toBe(s.dimension);
                for (let d = 0; d < s.dimension; ++d) {
                    expect(curve.getControl(i).values[d]).toBe(0);
                }
            }
            const p = curve.getPosition(0.5);
            for (let d = 0; d < s.dimension; ++d) {
                expect(p.values[d]).toBe(0);
            }
        });
    });
});
