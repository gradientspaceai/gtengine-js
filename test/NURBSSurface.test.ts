import { describe, expect, it } from 'vitest';
import { NURBSSurface } from '../src/NURBSSurface.js';
import { BSplineSurface } from '../src/BSplineSurface.js';
import { BasisFunctionInput, UniqueKnot } from '../src/BasisFunction.js';
import { ParametricSurface } from '../src/ParametricSurface.js';
import { Vector, length as vectorLength, sub } from '../src/Vector.js';

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
