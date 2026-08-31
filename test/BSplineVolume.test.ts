import { describe, it, expect } from 'vitest';
import { BSplineVolume } from '../src/BSplineVolume';
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
