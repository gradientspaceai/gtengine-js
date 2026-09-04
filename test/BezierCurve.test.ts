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


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). BezierCurve.h was read line by
// line against src/BezierCurve.ts. The Horner-style Compute with its
// combinatorial table is cross-checked against de Casteljau's algorithm for
// the position and against the hodograph (the de Casteljau evaluation of the
// scaled forward differences) for each of the first three derivatives, which
// is what the mControls[order] tables are supposed to encode.
import {
    check, fc, expectClose, expectVectorClose, wellScaledVector
} from './helpers/arbitraries.js';

// The de Casteljau evaluation above (shared with the earlier tests in this
// file) is the independent reference for the position.

// The r-th derivative of a degree-d Bezier curve is the degree-(d-r) Bezier
// curve on the r-th forward differences of the control points, scaled by
// d*(d-1)*...*(d-r+1).
function hodograph(controls: readonly Vector[], order: number): Vector[] {
    let p = controls.map(c => c.clone());
    const degree = controls.length - 1;
    let scale = 1;
    for (let r = 0; r < order; ++r) {
        const q: Vector[] = [];
        for (let i = 0; i + 1 < p.length; ++i) {
            const d = new Vector(p[i].size);
            for (let k = 0; k < d.size; ++k) {
                d.values[k] = p[i + 1].values[k] - p[i].values[k];
            }
            q.push(d);
        }
        p = q;
        scale *= degree - r;
    }
    return p.map(v => {
        const s = new Vector(v.size);
        for (let k = 0; k < s.size; ++k) { s.values[k] = scale * v.values[k]; }
        return s;
    });
}

const bezierData = (dim: number, minDegree = 2, maxDegree = 7):
    fc.Arbitrary<Vector[]> =>
    fc.integer({ min: minDegree, max: maxDegree }).chain(d =>
        fc.array(wellScaledVector(dim, -8, 8),
            { minLength: d + 1, maxLength: d + 1 }));

const inUnit = fc.double({ min: 0, max: 1, noNaN: true,
    noDefaultInfinity: true });

describe('BezierCurve verification', () => {
    it('agrees with de Casteljau for the position', () => {
        check(fc.tuple(bezierData(3), inUnit), ([controls, t]) => {
            const curve = new BezierCurve(3, controls.length - 1, controls);
            const jet = curve.createJet();
            curve.evaluate(t, 0, jet);
            // Both evaluations are numerically stable on [0,1] for control
            // points of this scale; the difference is a few ulps.
            expectVectorClose(jet[0], deCasteljau(controls, t), 1e-11, 1e-11);
        });
    });

    it('agrees with the hodograph for derivatives 1 through 3', () => {
        check(fc.tuple(bezierData(3), inUnit), ([controls, t]) => {
            const degree = controls.length - 1;
            const curve = new BezierCurve(3, degree, controls);
            const jet = curve.createJet();
            curve.evaluate(t, 3, jet);
            for (let order = 1; order <= 3; ++order) {
                if (order > degree) {
                    // Upstream only builds the third-difference table for
                    // degree >= 3 and returns a zero third derivative below
                    // that; a degree-2 curve has no third derivative anyway.
                    expectVectorClose(jet[order], new Vector(3), 1e-12, 1e-12);
                    continue;
                }
                const expected = deCasteljau(hodograph(controls, order), t);
                expectVectorClose(jet[order], expected, 1e-10, 1e-10);
            }
        });
    });

    it('interpolates the endpoints and their tangents exactly', () => {
        check(bezierData(2), controls => {
            const degree = controls.length - 1;
            const curve = new BezierCurve(2, degree, controls);
            const jet = curve.createJet();
            curve.evaluate(0, 1, jet);
            expectVectorClose(jet[0], controls[0], 1e-12, 1e-12);
            for (let k = 0; k < 2; ++k) {
                expectClose(jet[1].values[k], degree *
                    (controls[1].values[k] - controls[0].values[k]),
                1e-10, 1e-10);
            }
            curve.evaluate(1, 1, jet);
            expectVectorClose(jet[0], controls[degree], 1e-12, 1e-12);
            for (let k = 0; k < 2; ++k) {
                expectClose(jet[1].values[k], degree *
                    (controls[degree].values[k] -
                        controls[degree - 1].values[k]), 1e-10, 1e-10);
            }
        });
    });

    it('is equivariant under affine maps of the control points', () => {
        // X(t) = sum B_i(t) C_i with sum B_i(t) = 1, so A*X(t) + b is the
        // curve of the mapped control points.
        check(fc.tuple(bezierData(2), inUnit,
            fc.array(fc.double({ min: -3, max: 3, noNaN: true,
                noDefaultInfinity: true }), { minLength: 6, maxLength: 6 })),
        ([controls, t, m]) => {
            const mapped = controls.map(c => Vector.fromArray([
                m[0] * c.values[0] + m[1] * c.values[1] + m[2],
                m[3] * c.values[0] + m[4] * c.values[1] + m[5]]));
            const degree = controls.length - 1;
            const a = new BezierCurve(2, degree, controls);
            const b = new BezierCurve(2, degree, mapped);
            const ja = a.createJet();
            const jb = b.createJet();
            a.evaluate(t, 1, ja);
            b.evaluate(t, 1, jb);
            expectClose(jb[0].values[0],
                m[0] * ja[0].values[0] + m[1] * ja[0].values[1] + m[2],
                1e-9, 1e-9);
            expectClose(jb[0].values[1],
                m[3] * ja[0].values[0] + m[4] * ja[0].values[1] + m[5],
                1e-9, 1e-9);
            // The derivative transforms by the linear part only.
            expectClose(jb[1].values[0],
                m[0] * ja[1].values[0] + m[1] * ja[1].values[1], 1e-9, 1e-9);
            expectClose(jb[1].values[1],
                m[3] * ja[1].values[0] + m[4] * ja[1].values[1], 1e-9, 1e-9);
        });
    });

    it('stays inside the bounding box of its control polygon', () => {
        // The Bernstein basis is nonnegative and sums to 1 on [0,1], so the
        // curve lies in the convex hull of the control points; the axis
        // aligned box is a cheap consequence that any sign error in the
        // combinatorial table would violate.
        check(fc.tuple(bezierData(2), inUnit), ([controls, t]) => {
            const curve = new BezierCurve(2, controls.length - 1, controls);
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

    it('copies the control points and reports the domain', () => {
        check(bezierData(3), controls => {
            const curve = new BezierCurve(3, controls.length - 1, controls);
            expect(curve.getTMin()).toBe(0);
            expect(curve.getTMax()).toBe(1);
            expect(curve.getNumControls()).toBe(controls.length);
            expect(curve.getDegree()).toBe(controls.length - 1);
            expect(curve.isConstructed()).toBe(true);
            const stored = curve.getControls();
            for (let i = 0; i < controls.length; ++i) {
                expect(stored[i]).not.toBe(controls[i]);
                expectVectorClose(stored[i], controls[i], 0, 0);
            }
        });
    });

    it('zeroes the whole jet when the order is at least SUP_ORDER', () => {
        check(fc.tuple(bezierData(2), inUnit,
            fc.integer({ min: 4, max: 9 })), ([controls, t, order]) => {
            const curve = new BezierCurve(2, controls.length - 1, controls);
            const jet = curve.createJet();
            for (let i = 0; i < jet.length; ++i) { jet[i] = vec(5, 7); }
            curve.evaluate(t, order, jet);
            for (const j of jet) {
                expect(j.values).toEqual([0, 0]);
            }
        });
    });
});
