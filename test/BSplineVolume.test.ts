import { describe, it, expect } from 'vitest';
import { check, fc, finite, expectClose, expectVectorClose } from './helpers/arbitraries.js';
import { BSplineVolume } from '../src/BSplineVolume.js';
import { BSplineSurface } from '../src/BSplineSurface.js';
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

// A trilinear volume: degree 1 in each direction with two control points per
// direction, so X(u,v,w) is the trilinear interpolation of the eight
// controls, stored as control[i0 + 2*(i1 + 2*i2)].
function makeTrilinear(controls: Vector[]): BSplineVolume {
    const input = [new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1),
        new BasisFunctionInput(2, 1)];
    return new BSplineVolume(controls[0].size, input, controls);
}

describe('BSplineVolume construction', () => {
    it('reports the basis functions, domains and control counts', () => {
        const input = [new BasisFunctionInput(4, 2), new BasisFunctionInput(5, 3),
            new BasisFunctionInput(3, 1)];
        const volume = new BSplineVolume(3, input);
        expect(volume.isConstructed()).toBe(true);
        expect(volume.getDimension()).toBe(3);
        for (let dim = 0; dim < 3; ++dim) {
            expect(volume.getMinDomain(dim)).toBe(0);
            expect(volume.getMaxDomain(dim)).toBe(1);
        }
        expect(volume.getNumControls(0)).toBe(4);
        expect(volume.getNumControls(1)).toBe(5);
        expect(volume.getNumControls(2)).toBe(3);
        expect(volume.getBasisFunction(2).getDegree()).toBe(1);
    });

    it('zero-fills deferred control points and copies supplied ones', () => {
        const input = [new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1),
            new BasisFunctionInput(2, 1)];
        const deferred = new BSplineVolume(2, input);
        expect(deferred.getControl(1, 1, 1).values).toEqual([0, 0]);

        const controls: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            controls.push(vec(i, -i));
        }
        const volume = new BSplineVolume(2, input, controls);
        controls[7].values[0] = 999;
        expect(volume.getControl(1, 1, 1).values).toEqual([7, -7]);
        expect(volume.getControl(1, 0, 1).values).toEqual([5, -5]);

        const replacement = vec(42, 43);
        volume.setControl(0, 1, 0, replacement);
        replacement.values[0] = 0;
        expect(volume.getControl(0, 1, 0).values).toEqual([42, 43]);
        // Out-of-range setControl is a no-op and out-of-range getControl
        // returns the first control point.
        volume.setControl(2, 0, 0, vec(-1, -1));
        expect(volume.getControl(0, 0, 0).values).toEqual([0, -0]);
        expect(volume.getControl(9, 9, 9).values).toEqual([0, -0]);
    });
});

describe('BSplineVolume evaluation', () => {
    it('interpolates the corner control points of an open spline', () => {
        const controls: Vector[] = [];
        const rand = makeRandom(3);
        for (let i = 0; i < 8; ++i) {
            controls.push(vec(rand(), rand(), rand()));
        }
        const volume = makeTrilinear(controls);
        const jet = volume.createJet();

        const corners: [number, number, number, number][] = [
            [0, 0, 0, 0], [1, 0, 0, 1], [0, 1, 0, 2], [1, 1, 0, 3],
            [0, 0, 1, 4], [1, 0, 1, 5], [0, 1, 1, 6], [1, 1, 1, 7]
        ];
        for (const [u, v, w, index] of corners) {
            volume.evaluate(u, v, w, 0, jet);
            expect(jet[0].values).toEqual(controls[index].values);
        }
    });

    it('matches the analytic trilinear volume and its derivatives', () => {
        const c: Vector[] = [];
        const rand = makeRandom(808);
        for (let i = 0; i < 8; ++i) {
            c.push(vec(rand(), rand(), rand()));
        }
        const volume = makeTrilinear(c);
        const jet = volume.createJet();

        // f(u,v,w) = sum over the 8 corners of the trilinear weights.
        const at = (k: number, i0: number, i1: number, i2: number): number =>
            c[i0 + 2 * (i1 + 2 * i2)].values[k];

        for (let trial = 0; trial < 20; ++trial) {
            const u = rand(), v = rand(), w = rand();
            volume.evaluate(u, v, w, 2, jet);
            for (let k = 0; k < 3; ++k) {
                const g = (a: number, b: number, cc: number): number => at(k, a, b, cc);
                const wu = [1 - u, u], wv = [1 - v, v], ww = [1 - w, w];
                const du = [-1, 1];

                let position = 0, dU = 0, dV = 0, dW = 0;
                let dUV = 0, dUW = 0, dVW = 0;
                for (let i2 = 0; i2 < 2; ++i2) {
                    for (let i1 = 0; i1 < 2; ++i1) {
                        for (let i0 = 0; i0 < 2; ++i0) {
                            const value = g(i0, i1, i2);
                            position += wu[i0] * wv[i1] * ww[i2] * value;
                            dU += du[i0] * wv[i1] * ww[i2] * value;
                            dV += wu[i0] * du[i1] * ww[i2] * value;
                            dW += wu[i0] * wv[i1] * du[i2] * value;
                            dUV += du[i0] * du[i1] * ww[i2] * value;
                            dUW += du[i0] * wv[i1] * du[i2] * value;
                            dVW += wu[i0] * du[i1] * du[i2] * value;
                        }
                    }
                }

                expect(jet[0].values[k]).toBeCloseTo(position, 12);
                expect(jet[1].values[k]).toBeCloseTo(dU, 12);
                expect(jet[2].values[k]).toBeCloseTo(dV, 12);
                expect(jet[3].values[k]).toBeCloseTo(dW, 12);
                // The volume is linear in each variable separately.
                expect(jet[4].values[k]).toBeCloseTo(0, 12);
                expect(jet[5].values[k]).toBeCloseTo(0, 12);
                expect(jet[6].values[k]).toBeCloseTo(0, 12);
                expect(jet[7].values[k]).toBeCloseTo(dUV, 12);
                expect(jet[8].values[k]).toBeCloseTo(dUW, 12);
                expect(jet[9].values[k]).toBeCloseTo(dVW, 12);
            }
        }
    });

    it('reproduces the control point when all controls are equal (partition of unity)', () => {
        const input = [new BasisFunctionInput(5, 3), new BasisFunctionInput(4, 2),
            new BasisFunctionInput(3, 2)];
        const p = vec(-1, 2.5, 0.125);
        const controls: Vector[] = [];
        for (let i = 0; i < 5 * 4 * 3; ++i) {
            controls.push(p.clone());
        }
        const volume = new BSplineVolume(3, input, controls);
        const jet = volume.createJet();
        const rand = makeRandom(19);
        for (let i = 0; i < 20; ++i) {
            volume.evaluate(rand(), rand(), rand(), 1, jet);
            for (let k = 0; k < 3; ++k) {
                expect(jet[0].values[k]).toBeCloseTo(p.values[k], 12);
                expect(jet[1].values[k]).toBeCloseTo(0, 9);
                expect(jet[2].values[k]).toBeCloseTo(0, 9);
                expect(jet[3].values[k]).toBeCloseTo(0, 9);
            }
        }
    });

    it('matches finite differences for a cubic volume', () => {
        const input = [new BasisFunctionInput(6, 3), new BasisFunctionInput(5, 2),
            new BasisFunctionInput(4, 2)];
        const rand = makeRandom(64);
        const controls: Vector[] = [];
        for (let i = 0; i < 6 * 5 * 4; ++i) {
            controls.push(vec(rand(), rand(), rand()));
        }
        const volume = new BSplineVolume(3, input, controls);
        const jet = volume.createJet();
        const h = 1e-5;
        const position = (u: number, v: number, w: number): Vector => {
            const local = volume.createJet();
            volume.evaluate(u, v, w, 0, local);
            return local[0];
        };
        for (const [u, v, w] of [[0.25, 0.5, 0.75], [0.6, 0.2, 0.4]]) {
            volume.evaluate(u, v, w, 1, jet);
            const dU = jet[1].clone(), dV = jet[2].clone(), dW = jet[3].clone();
            const uP = position(u + h, v, w), uM = position(u - h, v, w);
            const vP = position(u, v + h, w), vM = position(u, v - h, w);
            const wP = position(u, v, w + h), wM = position(u, v, w - h);
            for (let k = 0; k < 3; ++k) {
                expect(dU.values[k]).toBeCloseTo((uP.values[k] - uM.values[k]) / (2 * h), 6);
                expect(dV.values[k]).toBeCloseTo((vP.values[k] - vM.values[k]) / (2 * h), 6);
                expect(dW.values[k]).toBeCloseTo((wP.values[k] - wM.values[k]) / (2 * h), 6);
            }
        }
    });

    it('returns a zero jet when the requested order is too large', () => {
        const controls: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            controls.push(vec(i + 1, i + 2, i + 3));
        }
        const volume = makeTrilinear(controls);
        const jet = volume.createJet();
        volume.evaluate(0.5, 0.5, 0.5, BSplineVolume.SUP_ORDER, jet);
        for (let i = 0; i < BSplineVolume.SUP_ORDER; ++i) {
            expect(jet[i].values).toEqual([0, 0, 0]);
        }
    });

    it('wraps the control indices for a periodic direction', () => {
        const numControls = 4, degree = 2;
        const periodic = new BasisFunctionInput();
        periodic.numControls = numControls;
        periodic.degree = degree;
        periodic.uniform = true;
        periodic.periodic = true;
        const numUnique = numControls + 2 * degree + 1;
        periodic.numUniqueKnots = numUnique;
        periodic.uniqueKnots = [];
        for (let i = 0; i < numUnique; ++i) {
            periodic.uniqueKnots.push(new UniqueKnot(i / (numUnique - 1), 1));
        }

        const input = [periodic, new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1)];
        const rand = makeRandom(1212);
        const controls: Vector[] = [];
        for (let i = 0; i < numControls * 2 * 2; ++i) {
            controls.push(vec(rand(), rand(), rand()));
        }
        const volume = new BSplineVolume(3, input, controls);
        expect(volume.getBasisFunction(0).isPeriodic()).toBe(true);

        const jet = volume.createJet();
        const umin = volume.getMinDomain(0), umax = volume.getMaxDomain(0);
        volume.evaluate(umin, 0.5, 0.25, 0, jet);
        const first = jet[0].clone();
        volume.evaluate(umax, 0.5, 0.25, 0, jet);
        for (let k = 0; k < 3; ++k) {
            expect(first.values[k]).toBeCloseTo(jet[0].values[k], 12);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V05). The volume is cross-checked against the separately
// ported BSplineSurface and BSplineCurve classes: X(u,v,w) is the B-spline
// curve in w whose control points are the w-slices evaluated as B-spline
// surfaces in (u,v). Those classes share no code with BSplineVolume::Compute,
// so the check exercises the index arithmetic and the triple summation
// independently.
// ---------------------------------------------------------------------------
describe('BSplineVolume verification', () => {
    interface Config {
        n: [number, number, number];
        d: [number, number, number];
        dim: number;
        controls: Vector[];
    }

    const configArb: fc.Arbitrary<Config> = fc.tuple(
        fc.integer({ min: 2, max: 4 }), fc.integer({ min: 1, max: 3 }),
        fc.integer({ min: 2, max: 4 }), fc.integer({ min: 1, max: 3 }),
        fc.integer({ min: 2, max: 4 }), fc.integer({ min: 1, max: 3 }),
        fc.integer({ min: 1, max: 3 }),
        fc.array(finite(-20, 20), { minLength: 4 * 4 * 4 * 3, maxLength: 4 * 4 * 4 * 3 })
    ).filter(([n0, d0, n1, d1, n2, d2]) => d0 < n0 && d1 < n1 && d2 < n2)
        .map(([n0, d0, n1, d1, n2, d2, dim, data]) => {
            const controls: Vector[] = [];
            for (let i = 0; i < n0 * n1 * n2; ++i) {
                const c = new Vector(dim);
                for (let k = 0; k < dim; ++k) { c.values[k] = data[i * 3 + k]; }
                controls.push(c);
            }
            return { n: [n0, n1, n2] as [number, number, number],
                d: [d0, d1, d2] as [number, number, number], dim, controls };
        });

    const inputsOf = (c: Config): BasisFunctionInput[] => [
        new BasisFunctionInput(c.n[0], c.d[0]),
        new BasisFunctionInput(c.n[1], c.d[1]),
        new BasisFunctionInput(c.n[2], c.d[2])
    ];

    const makeVolume = (c: Config): BSplineVolume =>
        new BSplineVolume(c.dim, inputsOf(c), c.controls);

    function position(volume: BSplineVolume, u: number, v: number, w: number): Vector {
        const jet = volume.createJet();
        volume.evaluate(u, v, w, 0, jet);
        return jet[0];
    }

    /** X(u,v,w) via a stack of BSplineSurface slices and a BSplineCurve in w. */
    function nestedEvaluation(c: Config, u: number, v: number, w: number): Vector {
        const [n0, n1, n2] = c.n;
        const slicePoints: Vector[] = [];
        for (let i2 = 0; i2 < n2; ++i2) {
            const slice: Vector[] = [];
            for (let i1 = 0; i1 < n1; ++i1) {
                for (let i0 = 0; i0 < n0; ++i0) {
                    slice.push(c.controls[i0 + n0 * (i1 + n1 * i2)]);
                }
            }
            const surface = new BSplineSurface(c.dim,
                [new BasisFunctionInput(n0, c.d[0]), new BasisFunctionInput(n1, c.d[1])],
                slice);
            slicePoints.push(surface.getPosition(u, v));
        }
        const curve = new BSplineCurve(c.dim,
            new BasisFunctionInput(n2, c.d[2]), slicePoints);
        return curve.getPosition(w);
    }

    it('agrees with a stack of BSplineSurface slices evaluated as a curve', () => {
        check(fc.tuple(configArb, finite(0, 1), finite(0, 1), finite(0, 1)),
            ([c, u, v, w]) => {
                const volume = makeVolume(c);
                expectVectorClose(position(volume, u, v, w),
                    nestedEvaluation(c, u, v, w), 1e-11, 1e-11);
                return true;
            }, 50);
    });

    it('is invariant under affine maps of the control points', () => {
        check(fc.tuple(configArb, finite(0, 1), finite(0, 1), finite(0, 1),
            finite(-3, 3), finite(-5, 5)), ([c, u, v, w, scale, shift]) => {
                const volume = makeVolume(c);
                const moved = c.controls.map(p => {
                    const q = new Vector(c.dim);
                    for (let k = 0; k < c.dim; ++k) { q.values[k] = scale * p.values[k] + shift; }
                    return q;
                });
                const other = makeVolume({ ...c, controls: moved });
                const p = position(volume, u, v, w);
                const q = position(other, u, v, w);
                for (let k = 0; k < c.dim; ++k) {
                    expectClose(q.values[k], scale * p.values[k] + shift, 1e-9, 1e-10);
                }
                return true;
            }, 50);
    });

    it('interpolates the eight corner control points of an open spline', () => {
        check(configArb, c => {
            const volume = makeVolume(c);
            const [n0, n1, n2] = c.n;
            for (const [a0, a1, a2] of [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0],
                [0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]) {
                const i0 = a0 * (n0 - 1);
                const i1 = a1 * (n1 - 1);
                const i2 = a2 * (n2 - 1);
                const p = position(volume, a0, a1, a2);
                expectVectorClose(p, c.controls[i0 + n0 * (i1 + n1 * i2)], 1e-11, 1e-11);
            }
            return true;
        }, 50);
    });

    it('first-order jet entries match central differences', () => {
        const smooth = configArb.filter(c => c.d[0] >= 2 && c.d[1] >= 2 && c.d[2] >= 2);
        check(fc.tuple(smooth, finite(0.25, 0.75), finite(0.25, 0.75), finite(0.25, 0.75)),
            ([c, u, v, w]) => {
                const volume = makeVolume(c);
                const jet = volume.createJet();
                volume.evaluate(u, v, w, 1, jet);
                const h = 1e-5;
                for (let k = 0; k < c.dim; ++k) {
                    const du = (position(volume, u + h, v, w).values[k]
                        - position(volume, u - h, v, w).values[k]) / (2 * h);
                    const dv = (position(volume, u, v + h, w).values[k]
                        - position(volume, u, v - h, w).values[k]) / (2 * h);
                    const dw = (position(volume, u, v, w + h).values[k]
                        - position(volume, u, v, w - h).values[k]) / (2 * h);
                    // C^1 surfaces: the truncation error is O(h) at a knot
                    // and the rounding error is O(eps/h).
                    expectClose(jet[1].values[k], du, 1e-4, 1e-5);
                    expectClose(jet[2].values[k], dv, 1e-4, 1e-5);
                    expectClose(jet[3].values[k], dw, 1e-4, 1e-5);
                }
                return true;
            }, 40);
    });

    it('setControl and getControl ignore out-of-range indices and copy inputs', () => {
        const c: Config = {
            n: [2, 2, 2], d: [1, 1, 1], dim: 2,
            controls: Array.from({ length: 8 }, (_, i) => {
                const p = new Vector(2);
                p.values[0] = i;
                return p;
            })
        };
        const volume = makeVolume(c);
        const replacement = new Vector(2);
        replacement.values[0] = 42;
        for (const [i0, i1, i2] of [[-1, 0, 0], [0, -1, 0], [0, 0, -1],
            [2, 0, 0], [0, 2, 0], [0, 0, 2]]) {
            volume.setControl(i0, i1, i2, replacement);
            expect(volume.getControl(i0, i1, i2).values[0]).toBe(0);
        }
        for (let i = 0; i < 8; ++i) {
            expect(volume.getControls()[i].values[0]).toBe(i);
        }
        volume.setControl(1, 1, 1, replacement);
        replacement.values[0] = -99;
        expect(volume.getControl(1, 1, 1).values[0]).toBe(42);
    });

    it('leaves the unused jet entries untouched for low orders', () => {
        const c: Config = {
            n: [3, 3, 3], d: [2, 2, 2], dim: 1,
            controls: Array.from({ length: 27 }, (_, i) => {
                const p = new Vector(1);
                p.values[0] = i;
                return p;
            })
        };
        const volume = makeVolume(c);
        const jet = volume.createJet();
        for (let i = 1; i < BSplineVolume.SUP_ORDER; ++i) { jet[i].values[0] = 777; }
        volume.evaluate(0.5, 0.5, 0.5, 0, jet);
        for (let i = 1; i < BSplineVolume.SUP_ORDER; ++i) {
            expect(jet[i].values[0]).toBe(777);
        }
        volume.evaluate(0.5, 0.5, 0.5, 1, jet);
        for (let i = 4; i < BSplineVolume.SUP_ORDER; ++i) {
            expect(jet[i].values[0]).toBe(777);
        }
        volume.evaluate(0.5, 0.5, 0.5, BSplineVolume.SUP_ORDER, jet);
        for (let i = 0; i < BSplineVolume.SUP_ORDER; ++i) {
            expect(jet[i].values[0]).toBe(0);
        }
    });

    // The second-order jet slots are ordered d2X/du2, d2X/dv2, d2X/dw2,
    // d2X/dudv, d2X/dudw, d2X/dvdw -- a different order from BSplineSurface.
    it('places the second-order derivatives in the upstream jet order', () => {
        const c: Config = {
            n: [3, 3, 3], d: [2, 2, 2], dim: 1,
            controls: Array.from({ length: 27 }, (_, i) => {
                const p = new Vector(1);
                // A control lattice whose spline is exactly u^2 + 2*v^2 + 3*w^2
                // is not available in closed form, so use distinct weights per
                // axis and compare with finite differences instead.
                const i0 = i % 3;
                const i1 = Math.floor(i / 3) % 3;
                const i2 = Math.floor(i / 9);
                p.values[0] = i0 + 10 * i1 + 100 * i2 + i0 * i1 * 0.5;
                return p;
            })
        };
        const volume = makeVolume(c);
        const jet = volume.createJet();
        const u = 0.4, v = 0.6, w = 0.55, h = 1e-4;
        volume.evaluate(u, v, w, 2, jet);

        const f = (a: number, b: number, cc: number): number =>
            position(volume, a, b, cc).values[0];
        const duu = (f(u + h, v, w) - 2 * f(u, v, w) + f(u - h, v, w)) / (h * h);
        const dvv = (f(u, v + h, w) - 2 * f(u, v, w) + f(u, v - h, w)) / (h * h);
        const dww = (f(u, v, w + h) - 2 * f(u, v, w) + f(u, v, w - h)) / (h * h);
        const duv = (f(u + h, v + h, w) - f(u + h, v - h, w)
            - f(u - h, v + h, w) + f(u - h, v - h, w)) / (4 * h * h);
        const duw = (f(u + h, v, w + h) - f(u + h, v, w - h)
            - f(u - h, v, w + h) + f(u - h, v, w - h)) / (4 * h * h);
        const dvw = (f(u, v + h, w + h) - f(u, v + h, w - h)
            - f(u, v - h, w + h) + f(u, v - h, w - h)) / (4 * h * h);

        // Second differences with h = 1e-4 carry an O(eps/h^2) rounding error
        // of about 1e-8 relative to the magnitudes here.
        expectClose(jet[4].values[0], duu, 1e-3, 1e-5);
        expectClose(jet[5].values[0], dvv, 1e-3, 1e-5);
        expectClose(jet[6].values[0], dww, 1e-3, 1e-5);
        expectClose(jet[7].values[0], duv, 1e-3, 1e-5);
        expectClose(jet[8].values[0], duw, 1e-3, 1e-5);
        expectClose(jet[9].values[0], dvw, 1e-3, 1e-5);
    });
});
