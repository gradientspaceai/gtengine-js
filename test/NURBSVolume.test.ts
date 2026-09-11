import { describe, expect, it } from 'vitest';
import { NURBSVolume } from '../src/NURBSVolume.js';
import { BSplineVolume } from '../src/BSplineVolume.js';
import { BasisFunctionInput, UniqueKnot } from '../src/BasisFunction.js';
import { Vector, length as vectorLength, sub } from '../src/Vector.js';
import {
    check, fc, finite, invertibleMatrix, wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

const SQRT_HALF = Math.SQRT1_2;

// The standard rational quadratic Bezier representation of a quarter circle of
// radius 1: control points (1,0), (1,1), (0,1) with weights 1, 1/sqrt(2), 1.
const ARC_POINTS: [number, number][] = [[1, 0], [1, 1], [0, 1]];
const ARC_WEIGHTS = [1, SQRT_HALF, 1];

// A quarter of a hollow cylinder: u sweeps the rational quarter circle, v
// interpolates the radius linearly from rMin to rMax and w is the linear
// extrusion along the z-axis to 'height'. The weights depend only on the u
// index, so they factor out of the v and w sums, giving exactly
//   X(u,v,w) = (a(u)*R(v), b(u)*R(v), height*w)
// with a^2+b^2 = 1 and R(v) = rMin + (rMax-rMin)*v.
function makeQuarterTube(rMin: number, rMax: number,
    height: number): NURBSVolume {
    const input = [new BasisFunctionInput(3, 2), new BasisFunctionInput(2, 1),
        new BasisFunctionInput(2, 1)];
    const controls: Vector[] = [];
    const weights: number[] = [];
    for (let i2 = 0; i2 < 2; ++i2) {
        for (let i1 = 0; i1 < 2; ++i1) {
            const radius = (i1 === 0 ? rMin : rMax);
            for (let i0 = 0; i0 < 3; ++i0) {
                const [a, b] = ARC_POINTS[i0];
                controls.push(v3(radius * a, radius * b, i2 * height));
                weights.push(ARC_WEIGHTS[i0]);
            }
        }
    }
    return new NURBSVolume(3, input, controls, weights);
}

// A control net with pseudo-random but deterministic entries.
function makeRandomNet(n0: number, n1: number, n2: number):
    { controls: Vector[], weights: number[] } {
    const controls: Vector[] = [];
    const weights: number[] = [];
    for (let i2 = 0; i2 < n2; ++i2) {
        for (let i1 = 0; i1 < n1; ++i1) {
            for (let i0 = 0; i0 < n0; ++i0) {
                const t = i0 + n0 * (i1 + n1 * i2);
                controls.push(v3(i0 + 0.3 * Math.sin(t),
                    i1 + 0.25 * Math.cos(2 * t),
                    i2 + 0.2 * Math.sin(3 * t + 1)));
                weights.push(1 + 0.5 * (1 + Math.sin(5 * t)));
            }
        }
    }
    return { controls, weights };
}

function positionOf(volume: NURBSVolume, u: number, v: number,
    w: number): Vector {
    const jet = volume.createJet();
    volume.evaluate(u, v, w, 0, jet);
    return jet[0];
}

describe('NURBSVolume construction and member access', () => {
    it('reports the domain, dimension and control counts', () => {
        const volume = makeQuarterTube(1, 2, 3);
        expect(volume.isConstructed()).toBe(true);
        expect(volume.getDimension()).toBe(3);
        for (let dim = 0; dim < 3; ++dim) {
            expect(volume.getMinDomain(dim)).toBe(0);
            expect(volume.getMaxDomain(dim)).toBe(1);
        }
        expect(volume.getNumControls(0)).toBe(3);
        expect(volume.getNumControls(1)).toBe(2);
        expect(volume.getNumControls(2)).toBe(2);
        expect(volume.getBasisFunction(0).getDegree()).toBe(2);
        expect(volume.getBasisFunction(1).getDegree()).toBe(1);
    });

    it('respects a non-unit knot domain', () => {
        const input0 = new BasisFunctionInput(3, 2);
        for (const knot of input0.uniqueKnots) {
            knot.t = -1 + 4 * knot.t;
        }
        const volume = new NURBSVolume(3, [input0,
            new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1)]);
        expect(volume.getMinDomain(0)).toBe(-1);
        expect(volume.getMaxDomain(0)).toBe(3);
        expect(volume.getMinDomain(1)).toBe(0);
        expect(volume.getMaxDomain(1)).toBe(1);
    });

    it('zero-fills controls and weights when they are deferred', () => {
        const volume = new NURBSVolume(3, [new BasisFunctionInput(3, 2),
            new BasisFunctionInput(2, 1), new BasisFunctionInput(4, 2)]);
        expect(volume.getControls().length).toBe(24);
        expect(volume.getWeights().length).toBe(24);
        for (const c of volume.getControls()) {
            expect(c.values).toEqual([0, 0, 0]);
        }
        for (const w of volume.getWeights()) {
            expect(w).toBe(0);
        }
    });

    it('copies the input controls (C++ value semantics)', () => {
        const { controls, weights } = makeRandomNet(3, 2, 2);
        const volume = new NURBSVolume(3, [new BasisFunctionInput(3, 2),
            new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1)],
            controls, weights);
        controls[0].set(0, 99);
        expect(volume.getControl(0, 0, 0).values[0]).not.toBe(99);
    });

    it('stores controls and weights in lexicographical order', () => {
        const n0 = 3, n1 = 2, n2 = 4;
        const { controls, weights } = makeRandomNet(n0, n1, n2);
        const volume = new NURBSVolume(3, [new BasisFunctionInput(n0, 2),
            new BasisFunctionInput(n1, 1), new BasisFunctionInput(n2, 3)],
            controls, weights);
        for (let i2 = 0; i2 < n2; ++i2) {
            for (let i1 = 0; i1 < n1; ++i1) {
                for (let i0 = 0; i0 < n0; ++i0) {
                    const index = i0 + n0 * (i1 + n1 * i2);
                    expect(volume.getControl(i0, i1, i2).values)
                        .toEqual(controls[index].values);
                    expect(volume.getWeight(i0, i1, i2)).toBe(weights[index]);
                }
            }
        }
    });

    it('setControl and setWeight ignore out-of-range indices and getters fall '
        + 'back to element 0', () => {
            const volume = makeQuarterTube(1, 2, 1);
            const before = volume.getControl(0, 0, 0).clone();
            volume.setControl(-1, 0, 0, v3(9, 9, 9));
            volume.setControl(0, 4, 0, v3(9, 9, 9));
            volume.setControl(0, 0, 7, v3(9, 9, 9));
            volume.setWeight(3, 0, 0, 42);
            expect(volume.getControl(0, 0, 0).values).toEqual(before.values);
            expect(volume.getControl(-1, 0, 0).values).toEqual(before.values);
            expect(volume.getWeight(9, 9, 9)).toBe(volume.getWeights()[0]);

            const p = v3(4, 5, 6);
            volume.setControl(1, 1, 1, p);
            p.set(0, 99);
            expect(volume.getControl(1, 1, 1).values).toEqual([4, 5, 6]);
            volume.setWeight(1, 1, 1, 0.25);
            expect(volume.getWeight(1, 1, 1)).toBe(0.25);
        });
});

describe('NURBSVolume exact quadric solids', () => {
    it('reproduces a quarter of a hollow cylinder', () => {
        const rMin = 1, rMax = 2.5, height = 3;
        const volume = makeQuarterTube(rMin, rMax, height);
        for (let iu = 0; iu <= 8; ++iu) {
            for (let iv = 0; iv <= 4; ++iv) {
                for (let iw = 0; iw <= 2; ++iw) {
                    const u = iu / 8, v = iv / 4, w = iw / 2;
                    const x = positionOf(volume, u, v, w);
                    const radius = rMin + (rMax - rMin) * v;
                    expect(Math.hypot(x.values[0], x.values[1]))
                        .toBeCloseTo(radius, 12);
                    expect(x.values[2]).toBeCloseTo(height * w, 12);
                }
            }
        }
        // Exact corners of the solid.
        expect(positionOf(volume, 0, 0, 0).values[0]).toBeCloseTo(rMin, 12);
        expect(positionOf(volume, 1, 1, 1).values[1]).toBeCloseTo(rMax, 12);
        expect(positionOf(volume, 1, 1, 1).values[2]).toBeCloseTo(height, 12);
        const mid = positionOf(volume, 0.5, 0.5, 0.5);
        const r = rMin + 0.5 * (rMax - rMin);
        expect(mid.values[0]).toBeCloseTo(r * SQRT_HALF, 12);
        expect(mid.values[1]).toBeCloseTo(r * SQRT_HALF, 12);
        expect(mid.values[2]).toBeCloseTo(0.5 * height, 12);
    });

    it('has the expected tangent structure on the tube', () => {
        const volume = makeQuarterTube(1, 2, 4);
        const jet = volume.createJet();
        for (const [u, v, w] of [[0.2, 0.3, 0.4], [0.6, 0.8, 0.1]]) {
            volume.evaluate(u, v, w, 1, jet);
            // dX/du is orthogonal to the radial direction (circular sweep).
            const radial = [jet[0].values[0], jet[0].values[1]];
            expect(radial[0] * jet[1].values[0] + radial[1] * jet[1].values[1])
                .toBeCloseTo(0, 10);
            // dX/dv is purely radial: parallel to (x,y,0).
            expect(jet[2].values[0] * radial[1] - jet[2].values[1] * radial[0])
                .toBeCloseTo(0, 10);
            expect(jet[2].values[2]).toBeCloseTo(0, 12);
            // dX/dw is the extrusion direction.
            expect(jet[3].values[0]).toBeCloseTo(0, 12);
            expect(jet[3].values[1]).toBeCloseTo(0, 12);
            expect(jet[3].values[2]).toBeCloseTo(4, 12);
        }
    });
});

describe('NURBSVolume reduction to BSplineVolume', () => {
    it('matches the B-spline volume when all weights are 1', () => {
        const n0 = 4, n1 = 4, n2 = 3;
        const { controls } = makeRandomNet(n0, n1, n2);
        const weights = new Array<number>(controls.length).fill(1);
        const input = () => [new BasisFunctionInput(n0, 3),
            new BasisFunctionInput(n1, 2), new BasisFunctionInput(n2, 2)];
        const nurbs = new NURBSVolume(3, input(), controls, weights);
        const bspline = new BSplineVolume(3, input(), controls);

        const jetN = nurbs.createJet();
        const jetB = bspline.createJet();
        for (let iu = 0; iu <= 4; ++iu) {
            for (let iv = 0; iv <= 4; ++iv) {
                for (let iw = 0; iw <= 4; ++iw) {
                    nurbs.evaluate(iu / 4, iv / 4, iw / 4, 2, jetN);
                    bspline.evaluate(iu / 4, iv / 4, iw / 4, 2, jetB);
                    for (let k = 0; k < NURBSVolume.SUP_ORDER; ++k) {
                        expect(vectorLength(sub(jetN[k], jetB[k])))
                            .toBeLessThan(1e-10);
                    }
                }
            }
        }
    });

    it('matches the B-spline volume when all weights are a common constant',
        () => {
            const { controls } = makeRandomNet(3, 3, 3);
            const weights = new Array<number>(controls.length).fill(4.25);
            const input = () => [new BasisFunctionInput(3, 2),
                new BasisFunctionInput(3, 2), new BasisFunctionInput(3, 2)];
            const nurbs = new NURBSVolume(3, input(), controls, weights);
            const bspline = new BSplineVolume(3, input(), controls);
            const jetN = nurbs.createJet();
            const jetB = bspline.createJet();
            for (const [u, v, w] of [[0.1, 0.4, 0.9], [0.5, 0.5, 0.5],
                [0.83, 0.27, 0.61]]) {
                nurbs.evaluate(u, v, w, 2, jetN);
                bspline.evaluate(u, v, w, 2, jetB);
                for (let k = 0; k < NURBSVolume.SUP_ORDER; ++k) {
                    expect(vectorLength(sub(jetN[k], jetB[k])))
                        .toBeLessThan(1e-10);
                }
            }
        });

    it('wraps the control indices for a periodic direction', () => {
        const numControls = 5, degree = 2;
        const makePeriodic = (): BasisFunctionInput => {
            const input = new BasisFunctionInput();
            input.numControls = numControls;
            input.degree = degree;
            input.uniform = true;
            input.periodic = true;
            const numUnique = numControls + 2 * degree + 1;
            input.numUniqueKnots = numUnique;
            input.uniqueKnots = [];
            for (let i = 0; i < numUnique; ++i) {
                input.uniqueKnots.push(new UniqueKnot(i / (numUnique - 1), 1));
            }
            return input;
        };
        const { controls } = makeRandomNet(numControls, 2, 2);
        const weights = new Array<number>(controls.length).fill(1);
        const input = () => [makePeriodic(), new BasisFunctionInput(2, 1),
            new BasisFunctionInput(2, 1)];
        const nurbs = new NURBSVolume(3, input(), controls, weights);
        const bspline = new BSplineVolume(3, input(), controls);
        expect(nurbs.getBasisFunction(0).isPeriodic()).toBe(true);

        const jetN = nurbs.createJet();
        const jetB = bspline.createJet();
        const umin = nurbs.getMinDomain(0), umax = nurbs.getMaxDomain(0);
        for (let iu = 0; iu <= 9; ++iu) {
            const u = umin + (umax - umin) * (iu / 9);
            nurbs.evaluate(u, 0.4, 0.7, 2, jetN);
            bspline.evaluate(u, 0.4, 0.7, 2, jetB);
            for (let k = 0; k < NURBSVolume.SUP_ORDER; ++k) {
                expect(vectorLength(sub(jetN[k], jetB[k]))).toBeLessThan(1e-10);
            }
        }
    });
});

describe('NURBSVolume derivatives', () => {
    const n0 = 4, n1 = 4, n2 = 4;
    const { controls, weights } = makeRandomNet(n0, n1, n2);
    const volume = new NURBSVolume(3, [new BasisFunctionInput(n0, 3),
        new BasisFunctionInput(n1, 3), new BasisFunctionInput(n2, 3)],
        controls, weights);

    it('matches central differences of the position', () => {
        const h = 1e-5;
        const jet = volume.createJet();
        for (const [u, v, w] of [[0.3, 0.4, 0.5], [0.5, 0.5, 0.5],
            [0.62, 0.17, 0.81]]) {
            volume.evaluate(u, v, w, 1, jet);
            const fd = [
                sub(positionOf(volume, u + h, v, w), positionOf(volume, u - h, v, w)),
                sub(positionOf(volume, u, v + h, w), positionOf(volume, u, v - h, w)),
                sub(positionOf(volume, u, v, w + h), positionOf(volume, u, v, w - h))];
            for (let d = 0; d < 3; ++d) {
                for (let k = 0; k < 3; ++k) {
                    expect(fd[d].values[k] / (2 * h))
                        .toBeCloseTo(jet[1 + d].values[k], 6);
                }
            }
        }
    });

    it('matches central differences of the first derivatives', () => {
        const h = 1e-4;
        // Index of the first-order derivative in a jet: 1 = u, 2 = v, 3 = w.
        const der = (order: number, u: number, v: number, w: number): Vector => {
            const jet = volume.createJet();
            volume.evaluate(u, v, w, 1, jet);
            return jet[order].clone();
        };
        const central = (order: number, axis: number, u: number, v: number,
            w: number): Vector => {
            const plus = [u, v, w], minus = [u, v, w];
            plus[axis] += h;
            minus[axis] -= h;
            const d = sub(der(order, plus[0], plus[1], plus[2]),
                der(order, minus[0], minus[1], minus[2]));
            for (let k = 0; k < 3; ++k) {
                d.values[k] /= 2 * h;
            }
            return d;
        };

        const jet = volume.createJet();
        for (const [u, v, w] of [[0.35, 0.45, 0.55], [0.7, 0.3, 0.6]]) {
            volume.evaluate(u, v, w, 2, jet);
            // d2X/du2, d2X/dv2, d2X/dw2.
            const pure = [[1, 0, 4], [2, 1, 5], [3, 2, 6]];
            for (const [order, axis, slot] of pure) {
                const fd = central(order, axis, u, v, w);
                for (let k = 0; k < 3; ++k) {
                    expect(fd.values[k]).toBeCloseTo(jet[slot].values[k], 5);
                }
            }
            // Mixed partials, taken both ways round to check symmetry:
            // jet[7] = d2X/dudv, jet[8] = d2X/dudw, jet[9] = d2X/dvdw.
            const mixed = [[1, 1, 2, 0, 7], [1, 2, 3, 0, 8], [2, 2, 3, 1, 9]];
            for (const [orderA, axisA, orderB, axisB, slot] of mixed) {
                const fdA = central(orderA, axisA, u, v, w);
                const fdB = central(orderB, axisB, u, v, w);
                for (let k = 0; k < 3; ++k) {
                    expect(fdA.values[k]).toBeCloseTo(jet[slot].values[k], 5);
                    expect(fdB.values[k]).toBeCloseTo(jet[slot].values[k], 5);
                }
            }
        }
    });

    it('leaves higher-order jet slots untouched for lower orders', () => {
        const jet = volume.createJet();
        for (let k = 0; k < NURBSVolume.SUP_ORDER; ++k) {
            jet[k] = v3(-1, -1, -1);
        }
        volume.evaluate(0.5, 0.5, 0.5, 1, jet);
        for (let k = 4; k < NURBSVolume.SUP_ORDER; ++k) {
            expect(jet[k].values).toEqual([-1, -1, -1]);
        }
    });

    it('returns a zero jet when the order is too large', () => {
        const jet = volume.createJet();
        for (let k = 0; k < NURBSVolume.SUP_ORDER; ++k) {
            jet[k] = v3(-1, -1, -1);
        }
        volume.evaluate(0.5, 0.5, 0.5, NURBSVolume.SUP_ORDER, jet);
        for (let k = 0; k < NURBSVolume.SUP_ORDER; ++k) {
            expect(jet[k].values).toEqual([0, 0, 0]);
        }
    });
});

describe('NURBSVolume degenerate and boundary behavior', () => {
    it('interpolates the corner control points of an open uniform solid', () => {
        const n0 = 4, n1 = 3, n2 = 3;
        const { controls, weights } = makeRandomNet(n0, n1, n2);
        const volume = new NURBSVolume(3, [new BasisFunctionInput(n0, 3),
            new BasisFunctionInput(n1, 2), new BasisFunctionInput(n2, 2)],
            controls, weights);
        for (const [u, v, w, i0, i1, i2] of [
            [0, 0, 0, 0, 0, 0], [1, 0, 0, n0 - 1, 0, 0],
            [0, 1, 0, 0, n1 - 1, 0], [0, 0, 1, 0, 0, n2 - 1],
            [1, 1, 1, n0 - 1, n1 - 1, n2 - 1]]) {
            expect(vectorLength(sub(positionOf(volume, u, v, w),
                volume.getControl(i0, i1, i2)))).toBeLessThan(1e-12);
        }
    });

    it('stays inside the bounding box of the control net for positive weights',
        () => {
            const { controls, weights } = makeRandomNet(4, 4, 3);
            const volume = new NURBSVolume(3, [new BasisFunctionInput(4, 3),
                new BasisFunctionInput(4, 2), new BasisFunctionInput(3, 2)],
                controls, weights);
            const lo = [0, 1, 2].map(k => Math.min(...controls.map(c => c.values[k])));
            const hi = [0, 1, 2].map(k => Math.max(...controls.map(c => c.values[k])));
            for (let iu = 0; iu <= 5; ++iu) {
                for (let iv = 0; iv <= 5; ++iv) {
                    for (let iw = 0; iw <= 5; ++iw) {
                        const x = positionOf(volume, iu / 5, iv / 5, iw / 5);
                        for (let k = 0; k < 3; ++k) {
                            expect(x.values[k]).toBeGreaterThanOrEqual(lo[k] - 1e-12);
                            expect(x.values[k]).toBeLessThanOrEqual(hi[k] + 1e-12);
                        }
                    }
                }
            }
        });

    it('clamps parameters outside the domain', () => {
        const volume = makeQuarterTube(1, 2, 1);
        expect(positionOf(volume, -1, -1, -1).values)
            .toEqual(positionOf(volume, 0, 0, 0).values);
        expect(positionOf(volume, 2, 2, 2).values)
            .toEqual(positionOf(volume, 1, 1, 1).values);
    });

    it('works in dimensions other than 3', () => {
        const controls: Vector[] = [];
        const weights: number[] = [];
        for (let i2 = 0; i2 < 2; ++i2) {
            for (let i1 = 0; i1 < 2; ++i1) {
                for (let i0 = 0; i0 < 2; ++i0) {
                    controls.push(Vector.fromArray([i0, i1, i2, i0 + i1 + i2]));
                    weights.push(1);
                }
            }
        }
        const volume = new NURBSVolume(4, [new BasisFunctionInput(2, 1),
            new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1)],
            controls, weights);
        expect(volume.getDimension()).toBe(4);
        // Trilinear interpolation of the unit cube corners.
        const x = positionOf(volume, 0.25, 0.5, 0.75);
        expect(x.values[0]).toBeCloseTo(0.25, 12);
        expect(x.values[1]).toBeCloseTo(0.5, 12);
        expect(x.values[2]).toBeCloseTo(0.75, 12);
        expect(x.values[3]).toBeCloseTo(1.5, 12);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V46): property-based checks against NURBSVolume.h.
// ---------------------------------------------------------------------------

function volumeNetArb(n0: number, n1: number, n2: number) {
    const n = n0 * n1 * n2;
    return fc.record({
        controls: fc.array(wellScaledVector(3, -5, 5),
            { minLength: n, maxLength: n }),
        weights: fc.array(finite(0.25, 4), { minLength: n, maxLength: n })
    });
}

function jetOfVolume(volume: NURBSVolume, u: number, v: number, w: number,
    order: number): Vector[] {
    const jet = volume.createJet();
    volume.evaluate(u, v, w, order, jet);
    return jet;
}


// A B-spline of degree d is only C^(d-1) at an interior knot, so a central
// difference of a derivative straddling a knot has O(h) error instead of
// O(h^2) and the comparison below would be meaningless. An open uniform
// basis with numControls control points of degree 'degree' has its interior
// knots at j/(numControls - degree); keep the sample away from all of them.
function awayFromKnots(numControls: number, degree: number, t: number,
    eps = 0.02): boolean {
    const numSpans = numControls - degree;
    for (let j = 1; j < numSpans; ++j) {
        if (Math.abs(t - j / numSpans) < eps) {
            return false;
        }
    }
    return true;
}

describe('NURBSVolume verification', () => {
    it('reduces to the B-spline volume when all weights are equal', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.integer({ min: 2, max: 4 }), fc.integer({ min: 2, max: 4 }),
            finite(0, 1), finite(0, 1), finite(0, 1), finite(0.25, 4))
            .chain(([n0, n1, n2, u, v, w, h]) => volumeNetArb(n0, n1, n2)
                .map(net => ({ n0, n1, n2, u, v, w, h, net }))),
            ({ n0, n1, n2, u, v, w, h, net }) => {
                const input = [
                    new BasisFunctionInput(n0, Math.min(2, n0 - 1)),
                    new BasisFunctionInput(n1, Math.min(2, n1 - 1)),
                    new BasisFunctionInput(n2, Math.min(2, n2 - 1))];
                const weights = net.weights.map(() => h);
                const nurbs = new NURBSVolume(3, input, net.controls, weights);
                const bspline = new BSplineVolume(3, input, net.controls);
                const jetN = jetOfVolume(nurbs, u, v, w, 2);
                const jetB = bspline.createJet();
                bspline.evaluate(u, v, w, 2, jetB);
                for (let k = 0; k < NURBSVolume.SUP_ORDER; ++k) {
                    for (let i = 0; i < 3; ++i) {
                        const a = jetN[k].values[i], b = jetB[k].values[i];
                        expect(Math.abs(a - b)).toBeLessThanOrEqual(
                            1e-11 * (1 + Math.abs(b)));
                    }
                }
            });
    });

    it('reproduces an affine map of the parameter cube exactly', () => {
        // With degree-1 bases and unit weights, the volume is the trilinear
        // interpolant of its control lattice. When the control points are an
        // affine image of the lattice coordinates i/(n-1), the interpolant is
        // that same affine map of (u,v,w).
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.integer({ min: 2, max: 4 }), fc.integer({ min: 2, max: 4 }),
            finite(0, 1), finite(0, 1), finite(0, 1),
            invertibleMatrix(3, 1e-2), wellScaledVector(3, -3, 3)),
            ([n0, n1, n2, u, v, w, M, t]) => {
                const input = [new BasisFunctionInput(n0, 1),
                    new BasisFunctionInput(n1, 1), new BasisFunctionInput(n2, 1)];
                const apply = (a: number, b: number, c: number): Vector => {
                    const y = new Vector(3);
                    for (let r = 0; r < 3; ++r) {
                        y.values[r] = M.get(r, 0) * a + M.get(r, 1) * b
                            + M.get(r, 2) * c + t.values[r];
                    }
                    return y;
                };
                const controls: Vector[] = [];
                const weights: number[] = [];
                for (let i2 = 0; i2 < n2; ++i2) {
                    for (let i1 = 0; i1 < n1; ++i1) {
                        for (let i0 = 0; i0 < n0; ++i0) {
                            controls.push(apply(i0 / (n0 - 1), i1 / (n1 - 1),
                                i2 / (n2 - 1)));
                            weights.push(1);
                        }
                    }
                }
                const volume = new NURBSVolume(3, input, controls, weights);
                const x = jetOfVolume(volume, u, v, w, 0)[0];
                const expected = apply(u, v, w);
                for (let r = 0; r < 3; ++r) {
                    expect(Math.abs(x.values[r] - expected.values[r]))
                        .toBeLessThanOrEqual(
                            1e-11 * (1 + Math.abs(expected.values[r])));
                }
            });
    });

    it('lies in the bounding box of the control net (positive weights)', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.integer({ min: 2, max: 4 }), fc.integer({ min: 2, max: 4 }),
            finite(0, 1), finite(0, 1), finite(0, 1))
            .chain(([n0, n1, n2, u, v, w]) => volumeNetArb(n0, n1, n2)
                .map(net => ({ n0, n1, n2, u, v, w, net }))),
            ({ n0, n1, n2, u, v, w, net }) => {
                const input = [
                    new BasisFunctionInput(n0, Math.min(2, n0 - 1)),
                    new BasisFunctionInput(n1, Math.min(2, n1 - 1)),
                    new BasisFunctionInput(n2, Math.min(2, n2 - 1))];
                const volume = new NURBSVolume(3, input, net.controls,
                    net.weights);
                const x = jetOfVolume(volume, u, v, w, 0)[0];
                for (let i = 0; i < 3; ++i) {
                    let lo = Number.POSITIVE_INFINITY;
                    let hi = Number.NEGATIVE_INFINITY;
                    for (const c of net.controls) {
                        lo = Math.min(lo, c.values[i]);
                        hi = Math.max(hi, c.values[i]);
                    }
                    const slack = 1e-12 * (1 + hi - lo);
                    expect(x.values[i]).toBeGreaterThanOrEqual(lo - slack);
                    expect(x.values[i]).toBeLessThanOrEqual(hi + slack);
                }
            });
    });

    it('is invariant when every weight is scaled by the same factor', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.integer({ min: 2, max: 4 }), fc.integer({ min: 2, max: 4 }),
            finite(0, 1), finite(0, 1), finite(0, 1), finite(0.25, 4))
            .chain(([n0, n1, n2, u, v, w, s]) => volumeNetArb(n0, n1, n2)
                .map(net => ({ n0, n1, n2, u, v, w, s, net }))),
            ({ n0, n1, n2, u, v, w, s, net }) => {
                const input = [
                    new BasisFunctionInput(n0, Math.min(2, n0 - 1)),
                    new BasisFunctionInput(n1, Math.min(2, n1 - 1)),
                    new BasisFunctionInput(n2, Math.min(2, n2 - 1))];
                const a = new NURBSVolume(3, input, net.controls, net.weights);
                const b = new NURBSVolume(3, input, net.controls,
                    net.weights.map(h => s * h));
                const xa = jetOfVolume(a, u, v, w, 1);
                const xb = jetOfVolume(b, u, v, w, 1);
                for (let k = 0; k < 4; ++k) {
                    for (let i = 0; i < 3; ++i) {
                        expect(Math.abs(xa[k].values[i] - xb[k].values[i]))
                            .toBeLessThanOrEqual(
                                1e-11 * (1 + Math.abs(xa[k].values[i])));
                    }
                }
            });
    });

    it('first-order derivatives match central differences', () => {
        const h = 1e-4;
        check(fc.tuple(fc.integer({ min: 3, max: 4 }),
            fc.integer({ min: 3, max: 4 }), fc.integer({ min: 3, max: 4 }),
            finite(0.2, 0.8), finite(0.2, 0.8), finite(0.2, 0.8))
            .filter(([n0, n1, n2, u, v, w]) => awayFromKnots(n0, 2, u)
                && awayFromKnots(n1, 2, v) && awayFromKnots(n2, 2, w))
            .chain(([n0, n1, n2, u, v, w]) => volumeNetArb(n0, n1, n2)
                .map(net => ({ n0, n1, n2, u, v, w, net }))),
            ({ n0, n1, n2, u, v, w, net }) => {
                const input = [new BasisFunctionInput(n0, 2),
                    new BasisFunctionInput(n1, 2), new BasisFunctionInput(n2, 2)];
                const volume = new NURBSVolume(3, input, net.controls,
                    net.weights);
                const jet = jetOfVolume(volume, u, v, w, 1);
                const args: Array<[number, number, number]> = [
                    [u + h, v, w], [u - h, v, w],
                    [u, v + h, w], [u, v - h, w],
                    [u, v, w + h], [u, v, w - h]];
                const p = args.map(a => jetOfVolume(volume, a[0], a[1], a[2], 0)[0]);
                for (let i = 0; i < 3; ++i) {
                    for (let d = 0; d < 3; ++d) {
                        const fd = (p[2 * d].values[i] - p[2 * d + 1].values[i])
                            / (2 * h);
                        expect(Math.abs(jet[1 + d].values[i] - fd)).toBeLessThan(
                            1e-5 * (1 + Math.abs(fd)));
                    }
                }
            });
    });

    it('second-order derivatives match central differences of the first', () => {
        const h = 1e-4;
        // jet[4..6] are d2/du2, d2/dv2, d2/dw2 and jet[7..9] are d2/dudv,
        // d2/dudw, d2/dvdw. Each is the derivative of jet[1+a] in direction b.
        const pairs: Array<[number, number, number]> = [
            [4, 0, 0], [5, 1, 1], [6, 2, 2], [7, 0, 1], [8, 0, 2], [9, 1, 2]];
        check(fc.tuple(fc.integer({ min: 3, max: 4 }),
            fc.integer({ min: 3, max: 4 }), fc.integer({ min: 3, max: 4 }),
            finite(0.2, 0.8), finite(0.2, 0.8), finite(0.2, 0.8))
            .filter(([n0, n1, n2, u, v, w]) => awayFromKnots(n0, 2, u)
                && awayFromKnots(n1, 2, v) && awayFromKnots(n2, 2, w))
            .chain(([n0, n1, n2, u, v, w]) => volumeNetArb(n0, n1, n2)
                .map(net => ({ n0, n1, n2, u, v, w, net }))),
            ({ n0, n1, n2, u, v, w, net }) => {
                const input = [new BasisFunctionInput(n0, 2),
                    new BasisFunctionInput(n1, 2), new BasisFunctionInput(n2, 2)];
                const volume = new NURBSVolume(3, input, net.controls,
                    net.weights);
                const jet = jetOfVolume(volume, u, v, w, 2);
                const base = [u, v, w];
                const shifted = (d: number, s: number): Vector[] => {
                    const a = base.slice();
                    a[d] += s * h;
                    return jetOfVolume(volume, a[0], a[1], a[2], 1);
                };
                const plus = [shifted(0, 1), shifted(1, 1), shifted(2, 1)];
                const minus = [shifted(0, -1), shifted(1, -1), shifted(2, -1)];
                for (const [slot, a, b] of pairs) {
                    for (let i = 0; i < 3; ++i) {
                        const fd = (plus[b][1 + a].values[i]
                            - minus[b][1 + a].values[i]) / (2 * h);
                        expect(Math.abs(jet[slot].values[i] - fd)).toBeLessThan(
                            1e-4 * (1 + Math.abs(fd)));
                    }
                }
            });
    });

    it('the indexed accessors use i0 + n0*(i1 + n1*i2) and copy the input',
        () => {
            check(fc.tuple(fc.integer({ min: 2, max: 4 }),
                fc.integer({ min: 2, max: 4 }), fc.integer({ min: 2, max: 4 }),
                wellScaledVector(3), finite(0.5, 3)),
                ([n0, n1, n2, p, h]) => {
                    const input = [new BasisFunctionInput(n0, 1),
                        new BasisFunctionInput(n1, 1),
                        new BasisFunctionInput(n2, 1)];
                    const volume = new NURBSVolume(3, input);
                    const i0 = n0 - 1, i1 = n1 - 1, i2 = n2 - 1;
                    const pCopy = p.clone();
                    volume.setControl(i0, i1, i2, p);
                    volume.setWeight(i0, i1, i2, h);
                    p.set(0, p.get(0) + 100);
                    const index = i0 + n0 * (i1 + n1 * i2);
                    expect(volume.getControl(i0, i1, i2).get(0)).toBe(pCopy.get(0));
                    expect(volume.getControls()[index].get(1)).toBe(pCopy.get(1));
                    expect(volume.getWeight(i0, i1, i2)).toBe(h);
                    expect(volume.getWeights()[index]).toBe(h);

                    // Out-of-range writes are ignored and out-of-range reads
                    // return element 0.
                    volume.setControl(n0, 0, 0, Vector.fromArray([9, 9, 9]));
                    volume.setWeight(0, n1, 0, 42);
                    expect(volume.getControl(0, 0, n2).get(0))
                        .toBe(volume.getControls()[0].get(0));
                    expect(volume.getWeight(-1, 0, 0))
                        .toBe(volume.getWeights()[0]);
                    for (const c of volume.getControls()) {
                        expect(c.get(0)).not.toBe(9);
                    }
                });
        });

    it('zeroes the entire jet for an out-of-range order', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 3 }),
            fc.integer({ min: 2, max: 3 }), fc.integer({ min: 2, max: 3 }),
            finite(0, 1), finite(0, 1), finite(0, 1),
            fc.integer({ min: NURBSVolume.SUP_ORDER, max: 16 }))
            .chain(([n0, n1, n2, u, v, w, order]) => volumeNetArb(n0, n1, n2)
                .map(net => ({ n0, n1, n2, u, v, w, order, net }))),
            ({ n0, n1, n2, u, v, w, order, net }) => {
                const input = [new BasisFunctionInput(n0, 1),
                    new BasisFunctionInput(n1, 1), new BasisFunctionInput(n2, 1)];
                const volume = new NURBSVolume(3, input, net.controls,
                    net.weights);
                const jet = volume.createJet();
                for (const e of jet) {
                    e.values[0] = 1;
                    e.values[1] = 2;
                    e.values[2] = 3;
                }
                volume.evaluate(u, v, w, order, jet);
                for (const e of jet) {
                    expect(e.values[0] + 0).toBe(0);
                    expect(e.values[1] + 0).toBe(0);
                    expect(e.values[2] + 0).toBe(0);
                }
            });
    });

    it('reports the basis-function domains and control counts', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 5 }),
            fc.integer({ min: 2, max: 5 }), fc.integer({ min: 2, max: 5 })),
            ([n0, n1, n2]) => {
                const input = [new BasisFunctionInput(n0, Math.min(2, n0 - 1)),
                    new BasisFunctionInput(n1, Math.min(2, n1 - 1)),
                    new BasisFunctionInput(n2, Math.min(2, n2 - 1))];
                const volume = new NURBSVolume(3, input);
                expect(volume.isConstructed()).toBe(true);
                expect(volume.getDimension()).toBe(3);
                const counts = [n0, n1, n2];
                for (let dim = 0; dim < 3; ++dim) {
                    expect(volume.getNumControls(dim)).toBe(counts[dim]);
                    expect(volume.getMinDomain(dim)).toBe(
                        volume.getBasisFunction(dim).getMinDomain());
                    expect(volume.getMaxDomain(dim)).toBe(
                        volume.getBasisFunction(dim).getMaxDomain());
                }
                expect(volume.getControls().length).toBe(n0 * n1 * n2);
                expect(volume.getWeights().length).toBe(n0 * n1 * n2);
            });
    });
});
