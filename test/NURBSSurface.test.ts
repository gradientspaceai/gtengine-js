import { describe, expect, it } from 'vitest';
import { NURBSSurface } from '../src/NURBSSurface.js';
import { BSplineSurface } from '../src/BSplineSurface.js';
import { BasisFunctionInput, UniqueKnot } from '../src/BasisFunction.js';
import { ParametricSurface } from '../src/ParametricSurface.js';
import { Vector, length as vectorLength, sub } from '../src/Vector.js';
import {
    check, fc, finite, invertibleMatrix, wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

const SQRT_HALF = Math.SQRT1_2;

// The standard rational quadratic Bezier representation of a quarter circle of
// radius 1 in a plane: control points (1,0), (1,1), (0,1) with weights
// 1, 1/sqrt(2), 1.
const ARC_POINTS: [number, number][] = [[1, 0], [1, 1], [0, 1]];
const ARC_WEIGHTS = [1, SQRT_HALF, 1];

// A quarter cylinder of radius 1 and height 'height': the u direction is the
// rational quarter circle, the v direction is a linear extrusion.
function makeQuarterCylinder(height: number): NURBSSurface {
    const input = [new BasisFunctionInput(3, 2), new BasisFunctionInput(2, 1)];
    const controls: Vector[] = [];
    const weights: number[] = [];
    for (let i1 = 0; i1 < 2; ++i1) {
        for (let i0 = 0; i0 < 3; ++i0) {
            const [a, b] = ARC_POINTS[i0];
            controls.push(v3(a, b, i1 * height));
            weights.push(ARC_WEIGHTS[i0]);
        }
    }
    return new NURBSSurface(3, input, controls, weights);
}

// An octant of the unit sphere as a surface of revolution: the u direction is
// the quarter-circle profile (r,z) in the xz-plane, the v direction rotates it
// through a quarter turn. The control point is (r_i*a_j, r_i*b_j, z_i) with
// weight w_i*w_j, so x = r(u)*a(v), y = r(u)*b(v), z = z(u) and
// x^2+y^2+z^2 = r^2+z^2 = 1.
function makeSphereOctant(): NURBSSurface {
    const input = [new BasisFunctionInput(3, 2), new BasisFunctionInput(3, 2)];
    const controls: Vector[] = [];
    const weights: number[] = [];
    for (let i1 = 0; i1 < 3; ++i1) {
        const [a, b] = ARC_POINTS[i1];
        for (let i0 = 0; i0 < 3; ++i0) {
            const [r, z] = ARC_POINTS[i0];
            controls.push(v3(r * a, r * b, z));
            weights.push(ARC_WEIGHTS[i0] * ARC_WEIGHTS[i1]);
        }
    }
    return new NURBSSurface(3, input, controls, weights);
}

// A control net with pseudo-random but deterministic entries.
function makeRandomNet(numControls0: number, numControls1: number):
    { controls: Vector[], weights: number[] } {
    const controls: Vector[] = [];
    const weights: number[] = [];
    for (let i1 = 0; i1 < numControls1; ++i1) {
        for (let i0 = 0; i0 < numControls0; ++i0) {
            const t = i0 + numControls0 * i1;
            controls.push(v3(i0 + 0.3 * Math.sin(t), i1 + 0.2 * Math.cos(2 * t),
                0.5 * Math.sin(3 * t + 1)));
            weights.push(1 + 0.5 * (1 + Math.sin(5 * t)));
        }
    }
    return { controls, weights };
}

function positionOf(surface: NURBSSurface, u: number, v: number): Vector {
    const jet = surface.createJet();
    surface.evaluate(u, v, 0, jet);
    return jet[0];
}

describe('NURBSSurface construction and member access', () => {
    it('takes its domain from the basis functions', () => {
        const surface = makeQuarterCylinder(2);
        expect(surface.isConstructed()).toBe(true);
        expect(surface.getDimension()).toBe(3);
        expect(surface.getUMin()).toBe(0);
        expect(surface.getUMax()).toBe(1);
        expect(surface.getVMin()).toBe(0);
        expect(surface.getVMax()).toBe(1);
        expect(surface.isRectangular()).toBe(true);
        expect(surface.getNumControls(0)).toBe(3);
        expect(surface.getNumControls(1)).toBe(2);
        expect(surface.getBasisFunction(0).getDegree()).toBe(2);
        expect(surface.getBasisFunction(1).getDegree()).toBe(1);
    });

    it('respects a non-unit knot domain', () => {
        const input = new BasisFunctionInput(3, 2);
        for (const knot of input.uniqueKnots) {
            knot.t = 2 + 3 * knot.t;
        }
        const surface = new NURBSSurface(3, [input, new BasisFunctionInput(2, 1)]);
        expect(surface.getUMin()).toBe(2);
        expect(surface.getUMax()).toBe(5);
        expect(surface.getVMin()).toBe(0);
        expect(surface.getVMax()).toBe(1);
    });

    it('zero-fills controls and weights when they are deferred', () => {
        const surface = new NURBSSurface(3,
            [new BasisFunctionInput(3, 2), new BasisFunctionInput(4, 2)]);
        expect(surface.getControls().length).toBe(12);
        expect(surface.getWeights().length).toBe(12);
        for (let i1 = 0; i1 < 4; ++i1) {
            for (let i0 = 0; i0 < 3; ++i0) {
                expect(surface.getControl(i0, i1).values).toEqual([0, 0, 0]);
                expect(surface.getWeight(i0, i1)).toBe(0);
            }
        }
    });

    it('copies the input controls (C++ value semantics)', () => {
        const controls = [v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0),
            v3(0, 1, 0), v3(1, 1, 0), v3(2, 1, 0)];
        const surface = new NURBSSurface(3,
            [new BasisFunctionInput(3, 2), new BasisFunctionInput(2, 1)],
            controls, [1, 1, 1, 1, 1, 1]);
        controls[0].set(0, 99);
        expect(surface.getControl(0, 0).values).toEqual([0, 0, 0]);
    });

    it('stores controls and weights in row-major order', () => {
        const { controls, weights } = makeRandomNet(3, 4);
        const surface = new NURBSSurface(3,
            [new BasisFunctionInput(3, 2), new BasisFunctionInput(4, 3)],
            controls, weights);
        for (let i1 = 0; i1 < 4; ++i1) {
            for (let i0 = 0; i0 < 3; ++i0) {
                const index = i0 + 3 * i1;
                expect(surface.getControl(i0, i1).values)
                    .toEqual(controls[index].values);
                expect(surface.getWeight(i0, i1)).toBe(weights[index]);
                expect(surface.getControls()[index].values)
                    .toEqual(controls[index].values);
                expect(surface.getWeights()[index]).toBe(weights[index]);
            }
        }
    });

    it('setControl and setWeight ignore out-of-range indices and getters '
        + 'fall back to element 0', () => {
            const surface = makeQuarterCylinder(1);
            const before = surface.getControl(0, 0).clone();
            surface.setControl(-1, 0, v3(9, 9, 9));
            surface.setControl(0, 5, v3(9, 9, 9));
            surface.setControl(3, 0, v3(9, 9, 9));
            surface.setWeight(-1, 0, 42);
            surface.setWeight(0, 2, 42);
            expect(surface.getControl(0, 0).values).toEqual(before.values);
            expect(surface.getControl(-1, 0).values).toEqual(before.values);
            expect(surface.getControl(0, 7).values).toEqual(before.values);
            expect(surface.getWeight(-1, 0)).toBe(surface.getWeights()[0]);
            expect(surface.getWeight(0, 7)).toBe(surface.getWeights()[0]);
        });

    it('setControl copies its argument and setWeight stores the value', () => {
        const surface = makeQuarterCylinder(1);
        const p = v3(4, 5, 6);
        surface.setControl(1, 1, p);
        p.set(0, 99);
        expect(surface.getControl(1, 1).values).toEqual([4, 5, 6]);
        surface.setWeight(1, 1, 0.25);
        expect(surface.getWeight(1, 1)).toBe(0.25);
    });
});

describe('NURBSSurface exact quadric patches', () => {
    it('reproduces a quarter cylinder of radius 1', () => {
        const height = 2.5;
        const surface = makeQuarterCylinder(height);
        for (let iu = 0; iu <= 10; ++iu) {
            const u = iu / 10;
            for (let iv = 0; iv <= 4; ++iv) {
                const v = iv / 4;
                const x = positionOf(surface, u, v);
                expect(Math.hypot(x.values[0], x.values[1])).toBeCloseTo(1, 12);
                expect(x.values[2]).toBeCloseTo(height * v, 12);
            }
        }
        // The end points and the midpoint of the arc are exact.
        expect(positionOf(surface, 0, 0).values[0]).toBeCloseTo(1, 12);
        expect(positionOf(surface, 0, 0).values[1]).toBeCloseTo(0, 12);
        expect(positionOf(surface, 1, 0).values[0]).toBeCloseTo(0, 12);
        expect(positionOf(surface, 1, 0).values[1]).toBeCloseTo(1, 12);
        const mid = positionOf(surface, 0.5, 0.5);
        expect(mid.values[0]).toBeCloseTo(SQRT_HALF, 12);
        expect(mid.values[1]).toBeCloseTo(SQRT_HALF, 12);
        expect(mid.values[2]).toBeCloseTo(0.5 * height, 12);
    });

    it('the u-tangent of the cylinder is orthogonal to the radius', () => {
        const surface = makeQuarterCylinder(3);
        for (let iu = 0; iu <= 8; ++iu) {
            const u = iu / 8;
            const jet = surface.createJet();
            surface.evaluate(u, 0.5, 1, jet);
            const radial = v3(jet[0].values[0], jet[0].values[1], 0);
            const tangent = jet[1];
            const d = radial.values[0] * tangent.values[0]
                + radial.values[1] * tangent.values[1];
            expect(d).toBeCloseTo(0, 10);
            // The v-tangent is the extrusion direction, of length 'height'.
            expect(jet[2].values[0]).toBeCloseTo(0, 12);
            expect(jet[2].values[1]).toBeCloseTo(0, 12);
            expect(jet[2].values[2]).toBeCloseTo(3, 12);
        }
    });

    it('reproduces an octant of the unit sphere', () => {
        const surface = makeSphereOctant();
        for (let iu = 0; iu <= 8; ++iu) {
            for (let iv = 0; iv <= 8; ++iv) {
                const x = positionOf(surface, iu / 8, iv / 8);
                expect(vectorLength(x)).toBeCloseTo(1, 12);
            }
        }
        // Corners of the patch: (1,0,0), (0,1,0) and the pole (0,0,1).
        expect(positionOf(surface, 0, 0).values[0]).toBeCloseTo(1, 12);
        expect(positionOf(surface, 0, 1).values[1]).toBeCloseTo(1, 12);
        expect(positionOf(surface, 1, 0.37).values[2]).toBeCloseTo(1, 12);
    });

    it('the sphere normal is parallel to the position', () => {
        const surface = makeSphereOctant();
        const jet = surface.createJet();
        for (const [u, v] of [[0.25, 0.25], [0.5, 0.75], [0.8, 0.4]]) {
            surface.evaluate(u, v, 1, jet);
            const x = jet[0];
            // Both tangents are orthogonal to the (unit) position vector.
            for (const t of [jet[1], jet[2]]) {
                const d = x.values[0] * t.values[0] + x.values[1] * t.values[1]
                    + x.values[2] * t.values[2];
                expect(d).toBeCloseTo(0, 10);
            }
        }
    });
});

describe('NURBSSurface reduction to BSplineSurface', () => {
    it('matches the B-spline surface when all weights are 1', () => {
        const numControls0 = 5, numControls1 = 4;
        const { controls } = makeRandomNet(numControls0, numControls1);
        const weights = new Array<number>(controls.length).fill(1);
        const input = () => [new BasisFunctionInput(numControls0, 3),
            new BasisFunctionInput(numControls1, 2)];
        const nurbs = new NURBSSurface(3, input(), controls, weights);
        const bspline = new BSplineSurface(3, input(), controls);

        const jetN = nurbs.createJet();
        const jetB = bspline.createJet();
        for (let iu = 0; iu <= 7; ++iu) {
            for (let iv = 0; iv <= 7; ++iv) {
                const u = iu / 7, v = iv / 7;
                nurbs.evaluate(u, v, 2, jetN);
                bspline.evaluate(u, v, 2, jetB);
                for (let k = 0; k < ParametricSurface.SUP_ORDER; ++k) {
                    expect(vectorLength(sub(jetN[k], jetB[k])))
                        .toBeLessThan(1e-10);
                }
            }
        }
    });

    it('matches the B-spline surface when all weights are a common constant',
        () => {
            // The weights appear in both numerator and denominator, so a
            // common scale factor cancels.
            const { controls } = makeRandomNet(4, 4);
            const weights = new Array<number>(controls.length).fill(7.5);
            const input = () => [new BasisFunctionInput(4, 3),
                new BasisFunctionInput(4, 3)];
            const nurbs = new NURBSSurface(3, input(), controls, weights);
            const bspline = new BSplineSurface(3, input(), controls);
            const jetN = nurbs.createJet();
            const jetB = bspline.createJet();
            for (const [u, v] of [[0.1, 0.9], [0.5, 0.5], [0.73, 0.21]]) {
                nurbs.evaluate(u, v, 2, jetN);
                bspline.evaluate(u, v, 2, jetB);
                for (let k = 0; k < ParametricSurface.SUP_ORDER; ++k) {
                    expect(vectorLength(sub(jetN[k], jetB[k])))
                        .toBeLessThan(1e-10);
                }
            }
        });

    it('matches the periodic B-spline surface when all weights are 1', () => {
        // Exercises the index wrapping in compute(...).
        const makePeriodic = (numControls: number,
            degree: number): BasisFunctionInput => {
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
        const numControls0 = 6, numControls1 = 5;
        const { controls } = makeRandomNet(numControls0, numControls1);
        const weights = new Array<number>(controls.length).fill(1);
        const input = () => [makePeriodic(numControls0, 2),
            new BasisFunctionInput(numControls1, 2)];
        const nurbs = new NURBSSurface(3, input(), controls, weights);
        const bspline = new BSplineSurface(3, input(), controls);
        const jetN = nurbs.createJet();
        const jetB = bspline.createJet();
        expect(nurbs.getBasisFunction(0).isPeriodic()).toBe(true);
        const umin = nurbs.getUMin(), umax = nurbs.getUMax();
        for (let iu = 0; iu <= 9; ++iu) {
            const u = umin + (umax - umin) * (iu / 9);
            for (let iv = 0; iv <= 5; ++iv) {
                nurbs.evaluate(u, iv / 5, 2, jetN);
                bspline.evaluate(u, iv / 5, 2, jetB);
                for (let k = 0; k < ParametricSurface.SUP_ORDER; ++k) {
                    expect(vectorLength(sub(jetN[k], jetB[k])))
                        .toBeLessThan(1e-10);
                }
            }
        }
    });
});

describe('NURBSSurface derivatives', () => {
    it('matches central differences of the position for nonuniform weights',
        () => {
            const numControls0 = 5, numControls1 = 4;
            const { controls, weights } = makeRandomNet(numControls0, numControls1);
            const surface = new NURBSSurface(3,
                [new BasisFunctionInput(numControls0, 3),
                    new BasisFunctionInput(numControls1, 2)], controls, weights);

            const h = 1e-5;
            const jet = surface.createJet();
            // The sample points avoid the interior knots (u = v = 0.5),
            // where a low-order derivative of the basis is discontinuous and
            // a centered difference is only first-order accurate.
            for (const [u, v] of [[0.3, 0.4], [0.55, 0.45], [0.62, 0.17],
                [0.85, 0.9]]) {
                surface.evaluate(u, v, 1, jet);
                const derU = jet[1].clone();
                const derV = jet[2].clone();
                const fdU = sub(positionOf(surface, u + h, v),
                    positionOf(surface, u - h, v));
                const fdV = sub(positionOf(surface, u, v + h),
                    positionOf(surface, u, v - h));
                for (let k = 0; k < 3; ++k) {
                    expect(fdU.values[k] / (2 * h)).toBeCloseTo(derU.values[k], 6);
                    expect(fdV.values[k] / (2 * h)).toBeCloseTo(derV.values[k], 6);
                }
            }
        });

    it('matches central differences of the first derivatives', () => {
        const numControls0 = 5, numControls1 = 5;
        const { controls, weights } = makeRandomNet(numControls0, numControls1);
        const surface = new NURBSSurface(3,
            [new BasisFunctionInput(numControls0, 3),
                new BasisFunctionInput(numControls1, 3)], controls, weights);

        const h = 1e-4;
        const derivative = (order: 1 | 2, u: number, v: number): Vector => {
            const jet = surface.createJet();
            surface.evaluate(u, v, 1, jet);
            return jet[order].clone();
        };

        const jet = surface.createJet();
        // As above, the sample points avoid the interior knot at 0.5.
        for (const [u, v] of [[0.35, 0.45], [0.62, 0.6], [0.72, 0.28]]) {
            surface.evaluate(u, v, 2, jet);
            const fdUU = sub(derivative(1, u + h, v), derivative(1, u - h, v));
            const fdUV = sub(derivative(1, u, v + h), derivative(1, u, v - h));
            const fdVV = sub(derivative(2, u, v + h), derivative(2, u, v - h));
            for (let k = 0; k < 3; ++k) {
                expect(fdUU.values[k] / (2 * h)).toBeCloseTo(jet[3].values[k], 5);
                expect(fdUV.values[k] / (2 * h)).toBeCloseTo(jet[4].values[k], 5);
                expect(fdVV.values[k] / (2 * h)).toBeCloseTo(jet[5].values[k], 5);
            }
            // The mixed partials commute: d2X/dudv computed from dX/dv also
            // matches jet[4].
            const fdVU = sub(derivative(2, u + h, v), derivative(2, u - h, v));
            for (let k = 0; k < 3; ++k) {
                expect(fdVU.values[k] / (2 * h)).toBeCloseTo(jet[4].values[k], 5);
            }
        }
    });

    it('leaves higher-order jet slots untouched for lower orders', () => {
        const surface = makeQuarterCylinder(1);
        const jet = surface.createJet();
        for (let k = 0; k < ParametricSurface.SUP_ORDER; ++k) {
            jet[k] = v3(-1, -1, -1);
        }
        surface.evaluate(0.5, 0.5, 0, jet);
        expect(jet[0].values[0]).toBeCloseTo(SQRT_HALF, 12);
        for (let k = 1; k < ParametricSurface.SUP_ORDER; ++k) {
            expect(jet[k].values).toEqual([-1, -1, -1]);
        }
    });

    it('returns a zero jet when the order is too large', () => {
        const surface = makeQuarterCylinder(1);
        const jet = surface.createJet();
        for (let k = 0; k < ParametricSurface.SUP_ORDER; ++k) {
            jet[k] = v3(-1, -1, -1);
        }
        surface.evaluate(0.5, 0.5, ParametricSurface.SUP_ORDER, jet);
        for (let k = 0; k < ParametricSurface.SUP_ORDER; ++k) {
            expect(jet[k].values).toEqual([0, 0, 0]);
        }
    });
});

describe('NURBSSurface degenerate and boundary behavior', () => {
    it('interpolates the corner control points of an open uniform patch', () => {
        const { controls, weights } = makeRandomNet(4, 3);
        const surface = new NURBSSurface(3,
            [new BasisFunctionInput(4, 3), new BasisFunctionInput(3, 2)],
            controls, weights);
        const corners: [number, number, number, number][] = [
            [0, 0, 0, 0], [1, 0, 3, 0], [0, 1, 0, 2], [1, 1, 3, 2]];
        for (const [u, v, i0, i1] of corners) {
            const x = positionOf(surface, u, v);
            expect(vectorLength(sub(x, surface.getControl(i0, i1))))
                .toBeLessThan(1e-12);
        }
    });

    it('clamps parameters outside the domain', () => {
        const surface = makeQuarterCylinder(2);
        expect(positionOf(surface, -1, -1).values)
            .toEqual(positionOf(surface, 0, 0).values);
        expect(positionOf(surface, 2, 2).values)
            .toEqual(positionOf(surface, 1, 1).values);
    });

    it('supports nonuniform knots', () => {
        const input0 = new BasisFunctionInput();
        input0.numControls = 5;
        input0.degree = 2;
        input0.uniform = false;
        input0.periodic = false;
        input0.numUniqueKnots = 4;
        input0.uniqueKnots = [new UniqueKnot(0, 3), new UniqueKnot(0.25, 1),
            new UniqueKnot(0.8, 1), new UniqueKnot(1, 3)];
        const { controls, weights } = makeRandomNet(5, 3);
        const surface = new NURBSSurface(3,
            [input0, new BasisFunctionInput(3, 2)], controls, weights);
        // The patch still interpolates its corner control points and every
        // point is a convex combination of the control points (positive
        // weights), so it lies in the bounding box of the net.
        expect(vectorLength(sub(positionOf(surface, 0, 0),
            surface.getControl(0, 0)))).toBeLessThan(1e-12);
        expect(vectorLength(sub(positionOf(surface, 1, 1),
            surface.getControl(4, 2)))).toBeLessThan(1e-12);
        for (let iu = 0; iu <= 6; ++iu) {
            for (let iv = 0; iv <= 6; ++iv) {
                const x = positionOf(surface, iu / 6, iv / 6);
                for (let k = 0; k < 3; ++k) {
                    const lo = Math.min(...controls.map(c => c.values[k]));
                    const hi = Math.max(...controls.map(c => c.values[k]));
                    expect(x.values[k]).toBeGreaterThanOrEqual(lo - 1e-12);
                    expect(x.values[k]).toBeLessThanOrEqual(hi + 1e-12);
                }
            }
        }
    });

    it('works in dimensions other than 3', () => {
        const controls = [Vector.fromArray([1, 0]), Vector.fromArray([1, 1]),
            Vector.fromArray([0, 1]), Vector.fromArray([2, 0]),
            Vector.fromArray([2, 2]), Vector.fromArray([0, 2])];
        const surface = new NURBSSurface(2,
            [new BasisFunctionInput(3, 2), new BasisFunctionInput(2, 1)],
            controls, [1, SQRT_HALF, 1, 1, SQRT_HALF, 1]);
        // A quarter annulus: the radius interpolates linearly from 1 to 2.
        for (let iu = 0; iu <= 5; ++iu) {
            for (let iv = 0; iv <= 5; ++iv) {
                const x = positionOf(surface, iu / 5, iv / 5);
                expect(Math.hypot(x.values[0], x.values[1]))
                    .toBeCloseTo(1 + iv / 5, 12);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V46): property-based checks against NURBSSurface.h.
// ---------------------------------------------------------------------------

// A control net of well-scaled points and strictly positive weights.
function netArb(numControls0: number, numControls1: number) {
    const n = numControls0 * numControls1;
    return fc.record({
        controls: fc.array(wellScaledVector(3, -5, 5),
            { minLength: n, maxLength: n }),
        weights: fc.array(finite(0.25, 4), { minLength: n, maxLength: n })
    });
}

// A jet of the requested order, as an array of SUP_ORDER vectors.
function jetOf(surface: NURBSSurface, u: number, v: number,
    order: number): Vector[] {
    const jet = surface.createJet();
    surface.evaluate(u, v, order, jet);
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

describe('NURBSSurface verification', () => {
    it('reduces to the B-spline surface when all weights are equal', () => {
        // The weight cancels from numerator and denominator exactly when it
        // is constant, so the only difference from BSplineSurface is the
        // division by a sum of basis values that is 1 up to rounding.
        check(fc.tuple(fc.integer({ min: 2, max: 5 }),
            fc.integer({ min: 2, max: 5 }), finite(0, 1), finite(0, 1),
            finite(0.25, 4)).chain(([n0, n1, u, v, w]) =>
                netArb(n0, n1).map(net => ({ n0, n1, u, v, w, net }))),
            ({ n0, n1, u, v, w, net }) => {
                const d0 = Math.min(2, n0 - 1), d1 = Math.min(2, n1 - 1);
                const input = [new BasisFunctionInput(n0, d0),
                    new BasisFunctionInput(n1, d1)];
                const weights = net.weights.map(() => w);
                const nurbs = new NURBSSurface(3, input, net.controls, weights);
                const bspline = new BSplineSurface(3, input, net.controls);
                const jetN = jetOf(nurbs, u, v, 2);
                const jetB = bspline.createJet();
                bspline.evaluate(u, v, 2, jetB);
                for (let k = 0; k < 6; ++k) {
                    for (let i = 0; i < 3; ++i) {
                        const a = jetN[k].values[i], b = jetB[k].values[i];
                        expect(Math.abs(a - b)).toBeLessThanOrEqual(
                            1e-11 * (1 + Math.abs(b)));
                    }
                }
            });
    });

    it('lies in the bounding box of the control net (positive weights)', () => {
        // With positive weights the surface point is a convex combination of
        // the control points, so it cannot leave their bounding box.
        check(fc.tuple(fc.integer({ min: 2, max: 5 }),
            fc.integer({ min: 2, max: 5 }), finite(0, 1), finite(0, 1))
            .chain(([n0, n1, u, v]) =>
                netArb(n0, n1).map(net => ({ n0, n1, u, v, net }))),
            ({ n0, n1, u, v, net }) => {
                const input = [new BasisFunctionInput(n0, Math.min(2, n0 - 1)),
                    new BasisFunctionInput(n1, Math.min(2, n1 - 1))];
                const surface = new NURBSSurface(3, input, net.controls,
                    net.weights);
                const x = jetOf(surface, u, v, 0)[0];
                for (let i = 0; i < 3; ++i) {
                    let lo = Number.POSITIVE_INFINITY;
                    let hi = Number.NEGATIVE_INFINITY;
                    for (const c of net.controls) {
                        lo = Math.min(lo, c.values[i]);
                        hi = Math.max(hi, c.values[i]);
                    }
                    // The convex-combination bound holds up to the rounding
                    // of the rational sum.
                    const slack = 1e-12 * (1 + hi - lo);
                    expect(x.values[i]).toBeGreaterThanOrEqual(lo - slack);
                    expect(x.values[i]).toBeLessThanOrEqual(hi + slack);
                }
            });
    });

    it('is invariant when every weight is scaled by the same factor', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 5 }),
            fc.integer({ min: 2, max: 5 }), finite(0, 1), finite(0, 1),
            finite(0.25, 4)).chain(([n0, n1, u, v, s]) =>
                netArb(n0, n1).map(net => ({ n0, n1, u, v, s, net }))),
            ({ n0, n1, u, v, s, net }) => {
                const input = [new BasisFunctionInput(n0, Math.min(2, n0 - 1)),
                    new BasisFunctionInput(n1, Math.min(2, n1 - 1))];
                const a = new NURBSSurface(3, input, net.controls, net.weights);
                const b = new NURBSSurface(3, input, net.controls,
                    net.weights.map(w => s * w));
                const xa = jetOf(a, u, v, 1);
                const xb = jetOf(b, u, v, 1);
                for (let k = 0; k < 3; ++k) {
                    for (let i = 0; i < 3; ++i) {
                        expect(Math.abs(xa[k].values[i] - xb[k].values[i]))
                            .toBeLessThanOrEqual(
                                1e-11 * (1 + Math.abs(xa[k].values[i])));
                    }
                }
            });
    });

    it('is affinely equivariant in the control points', () => {
        // A NURBS surface is a weighted average of its control points, so an
        // affine map of the net maps the surface.
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.integer({ min: 2, max: 4 }), finite(0, 1), finite(0, 1),
            invertibleMatrix(3, 1e-2), wellScaledVector(3, -3, 3))
            .chain(([n0, n1, u, v, M, t]) =>
                netArb(n0, n1).map(net => ({ n0, n1, u, v, M, t, net }))),
            ({ n0, n1, u, v, M, t, net }) => {
                const input = [new BasisFunctionInput(n0, Math.min(2, n0 - 1)),
                    new BasisFunctionInput(n1, Math.min(2, n1 - 1))];
                const a = new NURBSSurface(3, input, net.controls, net.weights);
                const mapped = net.controls.map(c => {
                    const y = new Vector(3);
                    for (let r = 0; r < 3; ++r) {
                        y.values[r] = M.get(r, 0) * c.values[0]
                            + M.get(r, 1) * c.values[1]
                            + M.get(r, 2) * c.values[2] + t.values[r];
                    }
                    return y;
                });
                const b = new NURBSSurface(3, input, mapped, net.weights);
                const xa = jetOf(a, u, v, 0)[0];
                const xb = jetOf(b, u, v, 0)[0];
                for (let r = 0; r < 3; ++r) {
                    const expected = M.get(r, 0) * xa.values[0]
                        + M.get(r, 1) * xa.values[1]
                        + M.get(r, 2) * xa.values[2] + t.values[r];
                    expect(Math.abs(xb.values[r] - expected))
                        .toBeLessThanOrEqual(1e-10 * (1 + Math.abs(expected)));
                }
            });
    });

    it('first-order derivatives match central differences', () => {
        const h = 1e-4;
        check(fc.tuple(fc.integer({ min: 3, max: 5 }),
            fc.integer({ min: 3, max: 5 }), finite(0.2, 0.8), finite(0.2, 0.8))
            .filter(([n0, n1, u, v]) => awayFromKnots(n0, 2, u)
                && awayFromKnots(n1, 2, v))
            .chain(([n0, n1, u, v]) =>
                netArb(n0, n1).map(net => ({ n0, n1, u, v, net }))),
            ({ n0, n1, u, v, net }) => {
                const input = [new BasisFunctionInput(n0, 2),
                    new BasisFunctionInput(n1, 2)];
                const surface = new NURBSSurface(3, input, net.controls,
                    net.weights);
                const jet = jetOf(surface, u, v, 1);
                const pu = jetOf(surface, u + h, v, 0)[0];
                const mu = jetOf(surface, u - h, v, 0)[0];
                const pv = jetOf(surface, u, v + h, 0)[0];
                const mv = jetOf(surface, u, v - h, 0)[0];
                for (let i = 0; i < 3; ++i) {
                    // The central difference has O(h^2) truncation error and
                    // O(eps/h) round-off, so ~1e-7 absolute at h = 1e-4.
                    const du = (pu.values[i] - mu.values[i]) / (2 * h);
                    const dv = (pv.values[i] - mv.values[i]) / (2 * h);
                    expect(Math.abs(jet[1].values[i] - du)).toBeLessThan(
                        1e-5 * (1 + Math.abs(du)));
                    expect(Math.abs(jet[2].values[i] - dv)).toBeLessThan(
                        1e-5 * (1 + Math.abs(dv)));
                }
            });
    });

    it('second-order derivatives match central differences of the first', () => {
        const h = 1e-4;
        check(fc.tuple(fc.integer({ min: 3, max: 5 }),
            fc.integer({ min: 3, max: 5 }), finite(0.2, 0.8), finite(0.2, 0.8))
            .filter(([n0, n1, u, v]) => awayFromKnots(n0, 2, u)
                && awayFromKnots(n1, 2, v))
            .chain(([n0, n1, u, v]) =>
                netArb(n0, n1).map(net => ({ n0, n1, u, v, net }))),
            ({ n0, n1, u, v, net }) => {
                const input = [new BasisFunctionInput(n0, 2),
                    new BasisFunctionInput(n1, 2)];
                const surface = new NURBSSurface(3, input, net.controls,
                    net.weights);
                const jet = jetOf(surface, u, v, 2);
                const pu = jetOf(surface, u + h, v, 1);
                const mu = jetOf(surface, u - h, v, 1);
                const pv = jetOf(surface, u, v + h, 1);
                const mv = jetOf(surface, u, v - h, 1);
                for (let i = 0; i < 3; ++i) {
                    const duu = (pu[1].values[i] - mu[1].values[i]) / (2 * h);
                    const duv = (pv[1].values[i] - mv[1].values[i]) / (2 * h);
                    const dvu = (pu[2].values[i] - mu[2].values[i]) / (2 * h);
                    const dvv = (pv[2].values[i] - mv[2].values[i]) / (2 * h);
                    expect(Math.abs(jet[3].values[i] - duu)).toBeLessThan(
                        1e-4 * (1 + Math.abs(duu)));
                    expect(Math.abs(jet[5].values[i] - dvv)).toBeLessThan(
                        1e-4 * (1 + Math.abs(dvv)));
                    // The mixed partial is symmetric, so both differences
                    // approximate jet[4].
                    expect(Math.abs(jet[4].values[i] - duv)).toBeLessThan(
                        1e-4 * (1 + Math.abs(duv)));
                    expect(Math.abs(jet[4].values[i] - dvu)).toBeLessThan(
                        1e-4 * (1 + Math.abs(dvu)));
                }
            });
    });

    it('setControl copies the input and getControl/getWeight use i0 + n0*i1',
        () => {
            check(fc.tuple(fc.integer({ min: 2, max: 5 }),
                fc.integer({ min: 2, max: 5 }), wellScaledVector(3), finite(0.5, 3)),
                ([n0, n1, p, w]) => {
                    const input = [new BasisFunctionInput(n0, 1),
                        new BasisFunctionInput(n1, 1)];
                    const surface = new NURBSSurface(3, input);
                    const i0 = n0 - 1, i1 = n1 - 1;
                    const pCopy = p.clone();
                    surface.setControl(i0, i1, p);
                    surface.setWeight(i0, i1, w);
                    p.set(0, p.get(0) + 100);
                    expect(surface.getControl(i0, i1).get(0)).toBe(pCopy.get(0));
                    expect(surface.getWeight(i0, i1)).toBe(w);
                    expect(surface.getControls()[i0 + n0 * i1].get(1))
                        .toBe(pCopy.get(1));
                    expect(surface.getWeights()[i0 + n0 * i1]).toBe(w);

                    // Out-of-range reads return element 0 and out-of-range
                    // writes are ignored, as upstream.
                    surface.setControl(-1, 0, Vector.fromArray([9, 9, 9]));
                    surface.setControl(n0, 0, Vector.fromArray([9, 9, 9]));
                    surface.setWeight(0, n1, 42);
                    expect(surface.getControl(-1, 0).get(0))
                        .toBe(surface.getControls()[0].get(0));
                    expect(surface.getWeight(n0, 0))
                        .toBe(surface.getWeights()[0]);
                    for (const c of surface.getControls()) {
                        expect(c.get(0)).not.toBe(9);
                    }
                });
        });

    it('zeroes the entire jet for an out-of-range order', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.integer({ min: 2, max: 4 }), finite(0, 1), finite(0, 1),
            fc.integer({ min: ParametricSurface.SUP_ORDER, max: 12 }))
            .chain(([n0, n1, u, v, order]) =>
                netArb(n0, n1).map(net => ({ n0, n1, u, v, order, net }))),
            ({ n0, n1, u, v, order, net }) => {
                const input = [new BasisFunctionInput(n0, 1),
                    new BasisFunctionInput(n1, 1)];
                const surface = new NURBSSurface(3, input, net.controls,
                    net.weights);
                const jet = surface.createJet();
                // Seed the jet with nonzero values so the zeroing is visible.
                for (const e of jet) {
                    e.values[0] = 1;
                    e.values[1] = 2;
                    e.values[2] = 3;
                }
                surface.evaluate(u, v, order, jet);
                for (const e of jet) {
                    expect(e.values[0] + 0).toBe(0);
                    expect(e.values[1] + 0).toBe(0);
                    expect(e.values[2] + 0).toBe(0);
                }
            });
    });

    it('the domain comes from the basis functions, not the (0,1) defaults',
        () => {
            check(fc.tuple(fc.integer({ min: 2, max: 5 }),
                fc.integer({ min: 2, max: 5 })), ([n0, n1]) => {
                    const input = [new BasisFunctionInput(n0, Math.min(2, n0 - 1)),
                        new BasisFunctionInput(n1, Math.min(2, n1 - 1))];
                    const surface = new NURBSSurface(3, input);
                    expect(surface.getUMin()).toBe(
                        surface.getBasisFunction(0).getMinDomain());
                    expect(surface.getUMax()).toBe(
                        surface.getBasisFunction(0).getMaxDomain());
                    expect(surface.getVMin()).toBe(
                        surface.getBasisFunction(1).getMinDomain());
                    expect(surface.getVMax()).toBe(
                        surface.getBasisFunction(1).getMaxDomain());
                    expect(surface.isRectangular()).toBe(true);
                    expect(surface.isConstructed()).toBe(true);
                    expect(surface.getNumControls(0)).toBe(n0);
                    expect(surface.getNumControls(1)).toBe(n1);
                });
        });
});
