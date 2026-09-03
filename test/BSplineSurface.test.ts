import { describe, it, expect } from 'vitest';
import { check, fc, finite, expectClose, expectVectorClose } from './helpers/arbitraries.js';
import { BSplineSurface } from '../src/BSplineSurface.js';
import { BSplineCurve } from '../src/BSplineCurve.js';
import { BasisFunctionInput, UniqueKnot } from '../src/BasisFunction.js';
import { ParametricSurface } from '../src/ParametricSurface.js';
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

// A bilinear surface: degree 1 in each direction with two control points per
// direction, so X(u,v) is the bilinear interpolation of the four controls.
function makeBilinear(controls: Vector[]): BSplineSurface {
    const input = [new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1)];
    return new BSplineSurface(controls[0].size, input, controls);
}

describe('BSplineSurface construction', () => {
    it('takes its domain from the basis functions', () => {
        const input = [new BasisFunctionInput(4, 2), new BasisFunctionInput(5, 3)];
        const surface = new BSplineSurface(3, input);
        expect(surface.isConstructed()).toBe(true);
        expect(surface.getUMin()).toBe(0);
        expect(surface.getUMax()).toBe(1);
        expect(surface.getVMin()).toBe(0);
        expect(surface.getVMax()).toBe(1);
        expect(surface.isRectangular()).toBe(true);
        expect(surface.getNumControls(0)).toBe(4);
        expect(surface.getNumControls(1)).toBe(5);
        expect(surface.getBasisFunction(0).getDegree()).toBe(2);
        expect(surface.getBasisFunction(1).getDegree()).toBe(3);
    });

    it('zero-fills deferred control points and copies supplied ones', () => {
        const input = [new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1)];
        const deferred = new BSplineSurface(3, input);
        for (let i1 = 0; i1 < 2; ++i1) {
            for (let i0 = 0; i0 < 2; ++i0) {
                expect(deferred.getControl(i0, i1).values).toEqual([0, 0, 0]);
            }
        }

        const controls = [vec(0, 0, 0), vec(1, 0, 0), vec(0, 1, 0), vec(1, 1, 5)];
        const surface = new BSplineSurface(3, input, controls);
        // The constructor copies, so mutating the input does not change the
        // surface.
        controls[3].values[2] = -100;
        expect(surface.getControl(1, 1).values).toEqual([1, 1, 5]);

        // setControl also copies, and out-of-range indices are ignored.
        const replacement = vec(9, 9, 9);
        surface.setControl(0, 0, replacement);
        replacement.values[0] = 0;
        expect(surface.getControl(0, 0).values).toEqual([9, 9, 9]);
        surface.setControl(5, 0, vec(-1, -1, -1));
        surface.setControl(0, -1, vec(-1, -1, -1));
        expect(surface.getControl(0, 0).values).toEqual([9, 9, 9]);
        // Out-of-range getControl returns the first control point.
        expect(surface.getControl(7, 7).values).toEqual([9, 9, 9]);
    });
});

describe('BSplineSurface evaluation', () => {
    it('interpolates the corner control points of an open spline', () => {
        const input = [new BasisFunctionInput(4, 3), new BasisFunctionInput(5, 2)];
        const rand = makeRandom(11);
        const controls: Vector[] = [];
        for (let i = 0; i < 20; ++i) {
            controls.push(vec(rand(), rand(), rand()));
        }
        const surface = new BSplineSurface(3, input, controls);
        const jet = surface.createJet();

        surface.evaluate(0, 0, 0, jet);
        expect(jet[0].values).toEqual(controls[0].values);
        surface.evaluate(1, 0, 0, jet);
        expect(jet[0].values).toEqual(controls[3].values);
        surface.evaluate(0, 1, 0, jet);
        expect(jet[0].values).toEqual(controls[16].values);
        surface.evaluate(1, 1, 0, jet);
        expect(jet[0].values).toEqual(controls[19].values);
    });

    it('reproduces the control point when all controls are equal (partition of unity)', () => {
        const input = [new BasisFunctionInput(6, 3), new BasisFunctionInput(4, 2)];
        const p = vec(1.5, -2.25, 7);
        const controls: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            controls.push(p.clone());
        }
        const surface = new BSplineSurface(3, input, controls);
        const jet = surface.createJet();
        const rand = makeRandom(77);
        for (let i = 0; i < 20; ++i) {
            surface.evaluate(rand(), rand(), 2, jet);
            for (let k = 0; k < 3; ++k) {
                expect(jet[0].values[k]).toBeCloseTo(p.values[k], 12);
                // The basis derivatives sum to zero.
                expect(jet[1].values[k]).toBeCloseTo(0, 10);
                expect(jet[2].values[k]).toBeCloseTo(0, 10);
                expect(jet[3].values[k]).toBeCloseTo(0, 9);
                expect(jet[4].values[k]).toBeCloseTo(0, 9);
                expect(jet[5].values[k]).toBeCloseTo(0, 9);
            }
        }
    });

    it('matches the analytic bilinear surface and its derivatives', () => {
        // X(u,v) = (1-u)(1-v)P00 + u(1-v)P10 + (1-u)v P01 + uv P11
        const p00 = vec(0, 0, 1), p10 = vec(2, 0, 3), p01 = vec(0, 4, -1), p11 = vec(2, 4, 8);
        const surface = makeBilinear([p00, p10, p01, p11]);
        const jet = surface.createJet();
        const rand = makeRandom(2024);
        for (let i = 0; i < 30; ++i) {
            const u = rand(), v = rand();
            surface.evaluate(u, v, 2, jet);
            for (let k = 0; k < 3; ++k) {
                const a = p00.values[k], b = p10.values[k];
                const c = p01.values[k], d = p11.values[k];
                const position = (1 - u) * (1 - v) * a + u * (1 - v) * b
                    + (1 - u) * v * c + u * v * d;
                const du = (1 - v) * (b - a) + v * (d - c);
                const dv = (1 - u) * (c - a) + u * (d - b);
                const duv = a - b - c + d;
                expect(jet[0].values[k]).toBeCloseTo(position, 12);
                expect(jet[1].values[k]).toBeCloseTo(du, 12);
                expect(jet[2].values[k]).toBeCloseTo(dv, 12);
                // The surface is linear in each variable separately.
                expect(jet[3].values[k]).toBeCloseTo(0, 12);
                expect(jet[4].values[k]).toBeCloseTo(duv, 12);
                expect(jet[5].values[k]).toBeCloseTo(0, 12);
            }
        }
    });

    it('matches finite differences of the position for a cubic-by-quadratic surface', () => {
        const input = [new BasisFunctionInput(6, 3), new BasisFunctionInput(5, 2)];
        const rand = makeRandom(4321);
        const controls: Vector[] = [];
        for (let i = 0; i < 30; ++i) {
            controls.push(vec(rand(), rand(), rand()));
        }
        const surface = new BSplineSurface(3, input, controls);
        const jet = surface.createJet();
        const h = 1e-5;
        for (const [u, v] of [[0.3, 0.4], [0.5, 0.55], [0.72, 0.2]]) {
            surface.evaluate(u, v, 1, jet);
            const du = jet[1].clone(), dv = jet[2].clone();
            const uPlus = surface.getPosition(u + h, v);
            const uMinus = surface.getPosition(u - h, v);
            const vPlus = surface.getPosition(u, v + h);
            const vMinus = surface.getPosition(u, v - h);
            for (let k = 0; k < 3; ++k) {
                expect(du.values[k]).toBeCloseTo(
                    (uPlus.values[k] - uMinus.values[k]) / (2 * h), 6);
                expect(dv.values[k]).toBeCloseTo(
                    (vPlus.values[k] - vMinus.values[k]) / (2 * h), 6);
            }
        }
    });

    it('returns a zero jet when the requested order is too large', () => {
        const surface = makeBilinear([vec(0, 0, 1), vec(1, 0, 2), vec(0, 1, 3), vec(1, 1, 4)]);
        const jet = surface.createJet();
        surface.evaluate(0.5, 0.5, ParametricSurface.SUP_ORDER, jet);
        for (let i = 0; i < ParametricSurface.SUP_ORDER; ++i) {
            expect(jet[i].values).toEqual([0, 0, 0]);
        }
    });

    it('clamps parameters outside the domain of an open spline', () => {
        const surface = makeBilinear([vec(0, 0, 0), vec(1, 0, 0), vec(0, 1, 0), vec(1, 1, 1)]);
        expect(surface.getPosition(-3, -3).values).toEqual([0, 0, 0]);
        expect(surface.getPosition(5, 5).values).toEqual([1, 1, 1]);
    });

    it('wraps the control indices for periodic splines', () => {
        // A periodic spline in u: the control points are not replicated in
        // storage (evaluate wraps the index instead), so the surface must be
        // closed in u.
        const numControls0 = 4, degree0 = 2;
        const input0 = new BasisFunctionInput();
        input0.numControls = numControls0;
        input0.degree = degree0;
        input0.uniform = true;
        input0.periodic = true;
        // The internal knot vector has numControls + 2*degree + 1 knots, all
        // with multiplicity 1 and uniformly spaced.
        const numUnique = numControls0 + 2 * degree0 + 1;
        input0.numUniqueKnots = numUnique;
        input0.uniqueKnots = [];
        for (let i = 0; i < numUnique; ++i) {
            input0.uniqueKnots.push(new UniqueKnot(i / (numUnique - 1), 1));
        }

        const input1 = new BasisFunctionInput(2, 1);
        const rand = makeRandom(5150);
        const controls: Vector[] = [];
        for (let i = 0; i < numControls0 * 2; ++i) {
            controls.push(vec(rand(), rand(), rand()));
        }

        const surface = new BSplineSurface(3, [input0, input1], controls);
        expect(surface.getBasisFunction(0).isPeriodic()).toBe(true);
        expect(surface.getNumControls(0)).toBe(numControls0);

        const umin = surface.getUMin(), umax = surface.getUMax();
        for (const v of [0, 0.3, 1]) {
            const first = surface.getPosition(umin, v);
            const last = surface.getPosition(umax, v);
            for (let k = 0; k < 3; ++k) {
                expect(first.values[k]).toBeCloseTo(last.values[k], 12);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V05). The strongest independent check available for a
// tensor-product surface is to rebuild it as a curve of curves using the
// separately ported BSplineCurve class, which shares no code with
// BSplineSurface::Compute.
// ---------------------------------------------------------------------------
describe('BSplineSurface verification', () => {
    interface Config {
        n0: number; d0: number; n1: number; d1: number; dim: number; controls: Vector[];
    }

    const configArb: fc.Arbitrary<Config> = fc.tuple(
        fc.integer({ min: 2, max: 6 }),   // numControls0
        fc.integer({ min: 1, max: 3 }),   // degree0
        fc.integer({ min: 2, max: 6 }),   // numControls1
        fc.integer({ min: 1, max: 3 }),   // degree1
        fc.integer({ min: 1, max: 3 }),   // dimension
        fc.array(finite(-20, 20), { minLength: 6 * 6 * 3, maxLength: 6 * 6 * 3 })
    ).filter(([n0, d0, n1, d1]) => d0 < n0 && d1 < n1)
        .map(([n0, d0, n1, d1, dim, data]) => {
            const controls: Vector[] = [];
            for (let i = 0; i < n0 * n1; ++i) {
                const c = new Vector(dim);
                for (let k = 0; k < dim; ++k) { c.values[k] = data[i * 3 + k]; }
                controls.push(c);
            }
            return { n0, d0, n1, d1, dim, controls };
        });

    function makeSurface(c: Config): BSplineSurface {
        const input = [new BasisFunctionInput(c.n0, c.d0), new BasisFunctionInput(c.n1, c.d1)];
        return new BSplineSurface(c.dim, input, c.controls);
    }

    /** X(u,v) computed as a B-spline curve whose control points are the
     *  points of the n1 u-direction B-spline curves through the control rows. */
    function nestedCurveEvaluation(c: Config, u: number, v: number): Vector {
        const rowPoints: Vector[] = [];
        for (let i1 = 0; i1 < c.n1; ++i1) {
            const row: Vector[] = [];
            for (let i0 = 0; i0 < c.n0; ++i0) {
                row.push(c.controls[i0 + c.n0 * i1]);
            }
            const curve = new BSplineCurve(c.dim, new BasisFunctionInput(c.n0, c.d0), row);
            rowPoints.push(curve.getPosition(u));
        }
        const outer = new BSplineCurve(c.dim, new BasisFunctionInput(c.n1, c.d1), rowPoints);
        return outer.getPosition(v);
    }

    it('agrees with a nested pair of BSplineCurve evaluations', () => {
        check(fc.tuple(configArb, finite(0, 1), finite(0, 1)), ([c, u, v]) => {
            const surface = makeSurface(c);
            const got = surface.getPosition(u, v);
            const want = nestedCurveEvaluation(c, u, v);
            // Only the summation order differs, so a tight relative tolerance
            // suffices.
            expectVectorClose(got, want, 1e-11, 1e-11);
            return true;
        }, 60);
    });

    it('is invariant under affine maps of the control points', () => {
        check(fc.tuple(configArb, finite(0, 1), finite(0, 1), finite(-3, 3), finite(-5, 5)),
            ([c, u, v, scale, shift]) => {
                const surface = makeSurface(c);
                const moved = c.controls.map(p => {
                    const q = new Vector(c.dim);
                    for (let k = 0; k < c.dim; ++k) { q.values[k] = scale * p.values[k] + shift; }
                    return q;
                });
                const other = makeSurface({ ...c, controls: moved });
                const p = surface.getPosition(u, v);
                const q = other.getPosition(u, v);
                for (let k = 0; k < c.dim; ++k) {
                    expectClose(q.values[k], scale * p.values[k] + shift, 1e-9, 1e-10);
                }
                return true;
            }, 60);
    });

    it('reproduces a repeated control point (partition of unity)', () => {
        check(fc.tuple(configArb, finite(0, 1), finite(0, 1), finite(-5, 5)),
            ([c, u, v, value]) => {
                const controls = c.controls.map(() => {
                    const p = new Vector(c.dim);
                    p.values.fill(value);
                    return p;
                });
                const surface = makeSurface({ ...c, controls });
                const p = surface.getPosition(u, v);
                for (let k = 0; k < c.dim; ++k) {
                    expectClose(p.values[k], value, 1e-12, 1e-12);
                }
                return true;
            }, 60);
    });

    // Degrees of at least 2 make the surface C^1, so a central difference
    // that straddles a knot is still accurate. For degree 1 the first
    // derivative is a step function and finite differences are meaningless
    // near a knot.
    it('first-order jet entries match central differences', () => {
        const smooth = configArb.filter(c => c.d0 >= 2 && c.d1 >= 2);
        check(fc.tuple(smooth, finite(0.2, 0.8), finite(0.2, 0.8)), ([c, u, v]) => {
            const surface = makeSurface(c);
            const jet = surface.createJet();
            surface.evaluate(u, v, 1, jet);
            const h = 1e-5;
            for (let k = 0; k < c.dim; ++k) {
                const du = (surface.getPosition(u + h, v).values[k]
                    - surface.getPosition(u - h, v).values[k]) / (2 * h);
                const dv = (surface.getPosition(u, v + h).values[k]
                    - surface.getPosition(u, v - h).values[k]) / (2 * h);
                // Central differences of a piecewise polynomial: the O(h^2)
                // truncation error plus the O(eps/h) rounding error dominate.
                expectClose(jet[1].values[k], du, 1e-4, 1e-5);
                expectClose(jet[2].values[k], dv, 1e-4, 1e-5);
            }
            return true;
        }, 60);
    });

    it('interpolates the corner control points of an open spline', () => {
        check(configArb, c => {
            const surface = makeSurface(c);
            const bf0 = surface.getBasisFunction(0);
            const bf1 = surface.getBasisFunction(1);
            const corners: Array<[number, number, number]> = [
                [bf0.getMinDomain(), bf1.getMinDomain(), 0],
                [bf0.getMaxDomain(), bf1.getMinDomain(), c.n0 - 1],
                [bf0.getMinDomain(), bf1.getMaxDomain(), c.n0 * (c.n1 - 1)],
                [bf0.getMaxDomain(), bf1.getMaxDomain(), c.n0 * c.n1 - 1]
            ];
            for (const [u, v, index] of corners) {
                expectVectorClose(surface.getPosition(u, v), c.controls[index], 1e-11, 1e-11);
            }
            return true;
        }, 60);
    });

    it('setControl and getControl ignore out-of-range indices', () => {
        const c: Config = {
            n0: 3, d0: 1, n1: 3, d1: 1, dim: 2,
            controls: Array.from({ length: 9 }, (_, i) => {
                const p = new Vector(2);
                p.values[0] = i;
                p.values[1] = -i;
                return p;
            })
        };
        const surface = makeSurface(c);
        const replacement = new Vector(2);
        replacement.values[0] = 42;
        for (const [i0, i1] of [[-1, 0], [0, -1], [3, 0], [0, 3]]) {
            surface.setControl(i0, i1, replacement);
            // Out-of-range reads return the first control point.
            expect(surface.getControl(i0, i1).values[0]).toBe(0);
        }
        // Nothing was written.
        for (let i = 0; i < 9; ++i) {
            expect(surface.getControls()[i].values[0]).toBe(i);
        }
        // setControl copies its argument, as upstream does.
        surface.setControl(1, 1, replacement);
        replacement.values[0] = -99;
        expect(surface.getControl(1, 1).values[0]).toBe(42);
    });

    it('leaves the unused jet entries untouched for low orders', () => {
        const c: Config = {
            n0: 3, d0: 2, n1: 3, d1: 2, dim: 2,
            controls: Array.from({ length: 9 }, (_, i) => {
                const p = new Vector(2);
                p.values[0] = i;
                p.values[1] = i * i;
                return p;
            })
        };
        const surface = makeSurface(c);
        const jet = surface.createJet();
        for (let i = 1; i < ParametricSurface.SUP_ORDER; ++i) {
            jet[i].values[0] = 777;
        }
        surface.evaluate(0.5, 0.5, 0, jet);
        for (let i = 1; i < ParametricSurface.SUP_ORDER; ++i) {
            expect(jet[i].values[0]).toBe(777);
        }
        // An order at or above SUP_ORDER zeroes the whole jet.
        surface.evaluate(0.5, 0.5, ParametricSurface.SUP_ORDER, jet);
        for (let i = 0; i < ParametricSurface.SUP_ORDER; ++i) {
            expect(jet[i].values[0]).toBe(0);
        }
    });
});
