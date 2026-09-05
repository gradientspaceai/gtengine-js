import { describe, it, expect } from 'vitest';
import {
    IntpBSplineUniform, IntpBSplineUniform1, IntpBSplineUniform2,
    IntpBSplineUniform3, IntpBSplineUniformCacheMode, IntpBSplineUniformShared
} from '../src/IntpBSplineUniform.js';
import type { IntpBSplineUniformControls } from '../src/IntpBSplineUniform.js';
import { Vector } from '../src/Vector.js';
import { add as vadd, mul as vmul } from '../src/Vector.js';
import { check, expectClose, fc, finite, seededRandom }
    from './helpers/arbitraries.js';

const NO_CACHING = IntpBSplineUniformCacheMode.NO_CACHING;
const PRE_CACHING = IntpBSplineUniformCacheMode.PRE_CACHING;
const ON_DEMAND_CACHING = IntpBSplineUniformCacheMode.ON_DEMAND_CACHING;
const ALL_MODES = [NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING];

// A controls adapter over a flat array of numbers. The index tuple is
// stored with index 0 varying fastest.
class NumberControls implements IntpBSplineUniformControls<number> {
    constructor(public sizes: number[], public data: number[]) {
    }

    getSize(dimension: number): number {
        return this.sizes[dimension];
    }

    get(indices: readonly number[]): number {
        let index = 0;
        for (let d = this.sizes.length - 1; d >= 0; --d) {
            index = index * this.sizes[d] + indices[d];
        }
        return this.data[index];
    }

    add(c0: number, c1: number): number {
        return c0 + c1;
    }

    mul(c0: number, s: number): number {
        return c0 * s;
    }
}

// A controls adapter whose control points are Vector objects, to exercise
// the generic control-point type.
class VectorControls implements IntpBSplineUniformControls<Vector> {
    constructor(public size: number, public data: Vector[]) {
    }

    getSize(_dimension: number): number {
        return this.size;
    }

    get(indices: readonly number[]): Vector {
        return this.data[indices[0]];
    }

    add(c0: Vector, c1: Vector): Vector {
        return vadd(c0, c1);
    }

    mul(c0: Vector, s: number): Vector {
        return vmul(c0, s);
    }
}

// The uniform B-spline basis functions on [0,1], written directly from the
// standard formulas. Index j is the weight of control point i + j.
function basis(degree: number, u: number): number[] {
    switch (degree) {
        case 1:
            return [1 - u, u];
        case 2:
            return [
                0.5 * (1 - u) * (1 - u),
                0.5 * (1 + 2 * u - 2 * u * u),
                0.5 * u * u
            ];
        case 3:
            return [
                (1 - u) * (1 - u) * (1 - u) / 6,
                (3 * u * u * u - 6 * u * u + 4) / 6,
                (-3 * u * u * u + 3 * u * u + 3 * u + 1) / 6,
                u * u * u / 6
            ];
    }
    throw new Error('unsupported degree');
}

// The independently computed 1-dimensional B-spline value.
function reference1(degree: number, controls: number[], t: number): number {
    const c = controls.length;
    const tmin = -0.5, tmax = c - 0.5;
    const dsdt = (c - degree) / (tmax - tmin);
    let i: number, u: number;
    if (t > tmin) {
        if (t < tmax) {
            const smd = dsdt * (t - tmin);
            i = Math.floor(smd);
            u = smd - i;
        }
        else {
            i = c - 1 - degree;
            u = 1;
        }
    }
    else {
        i = 0;
        u = 0;
    }
    const N = basis(degree, u);
    let sum = 0;
    for (let j = 0; j <= degree; ++j) {
        sum += controls[i + j] * N[j];
    }
    return sum;
}

describe('IntpBSplineUniform static helpers', () => {
    it('computes the blending matrix for degree 0', () => {
        expect(IntpBSplineUniformShared.computeBlendingMatrix(0)).toEqual([1]);
    });

    it('computes the blending matrix for degree 1', () => {
        // Row 0 is Q_{1,1}(s) = 1 - s, row 1 is Q_{1,0}(s) = s.
        const A = IntpBSplineUniformShared.computeBlendingMatrix(1);
        expect(A.length).toBe(4);
        const expected = [1, -1, 0, 1];
        for (let i = 0; i < 4; ++i) {
            expect(A[i]).toBeCloseTo(expected[i], 14);
        }
    });

    it('computes the blending matrix for degree 2', () => {
        // Row 0 is Q_{2,2}(s) = (1 - s)^2/2, row 1 is
        // Q_{2,1}(s) = (1 + 2s - 2s^2)/2, row 2 is Q_{2,0}(s) = s^2/2.
        const A = IntpBSplineUniformShared.computeBlendingMatrix(2);
        expect(A.length).toBe(9);
        const expected = [0.5, -1, 0.5, 0.5, 1, -1, 0, 0, 0.5];
        for (let i = 0; i < 9; ++i) {
            expect(A[i]).toBeCloseTo(expected[i], 14);
        }
    });

    it('computes the blending matrix rows that match the standard basis', () => {
        for (const degree of [1, 2, 3]) {
            const A = IntpBSplineUniformShared.computeBlendingMatrix(degree);
            const degreeP1 = degree + 1;
            for (const u of [0, 0.25, 0.5, 0.75, 1]) {
                const N = basis(degree, u);
                for (let row = 0; row <= degree; ++row) {
                    let value = 0;
                    for (let col = degree; col >= 0; --col) {
                        value = value * u + A[col + degreeP1 * row];
                    }
                    // Row r holds the weight of control point i + r.
                    expect(value).toBeCloseTo(N[row], 13);
                }
            }
        }
    });

    it('computes the derivative coefficients', () => {
        const d2 = IntpBSplineUniformShared.computeDCoefficients(2);
        expect(d2.dCoefficients).toEqual([1, 1, 1, 1, 2, 2]);
        expect(d2.ellMax).toEqual([2, 4, 5]);

        const d0 = IntpBSplineUniformShared.computeDCoefficients(0);
        expect(d0.dCoefficients).toEqual([1]);
        expect(d0.ellMax).toEqual([0]);

        const d3 = IntpBSplineUniformShared.computeDCoefficients(3);
        expect(d3.dCoefficients.length).toBe(10);
        expect(d3.ellMax).toEqual([3, 6, 8, 9]);
    });

    it('computes the powers of ds/dt', () => {
        // numControls = 8, degree = 3, tmin = -0.5, tmax = 7.5, so
        // ds/dt = (8 - 3)/8 = 0.625.
        const p = IntpBSplineUniformShared.computePowers(3, 8, -0.5, 7.5);
        expect(p.length).toBe(4);
        expect(p[0]).toBe(1);
        expect(p[1]).toBeCloseTo(0.625, 15);
        expect(p[2]).toBeCloseTo(0.625 ** 2, 15);
        expect(p[3]).toBeCloseTo(0.625 ** 3, 15);
    });

    it('computes the key for the interval and local parameter', () => {
        const c = 5, degree = 2;
        const tmin = -0.5, tmax = c - 0.5;
        const dsdt = (c - degree) / (tmax - tmin);  // 0.6

        // Below or at tmin.
        expect(IntpBSplineUniformShared.getKey(-10, tmin, tmax, dsdt, c, degree))
            .toEqual({ index: 0, u: 0 });
        expect(IntpBSplineUniformShared.getKey(tmin, tmin, tmax, dsdt, c, degree))
            .toEqual({ index: 0, u: 0 });

        // At or above tmax the domain is extended to the closed support.
        expect(IntpBSplineUniformShared.getKey(tmax, tmin, tmax, dsdt, c, degree))
            .toEqual({ index: c - 1 - degree, u: 1 });
        expect(IntpBSplineUniformShared.getKey(100, tmin, tmax, dsdt, c, degree))
            .toEqual({ index: c - 1 - degree, u: 1 });

        // Interior: smd = 0.6*(t + 0.5).
        const key = IntpBSplineUniformShared.getKey(2, tmin, tmax, dsdt, c, degree);
        expect(key.index).toBe(1);
        expect(key.u).toBeCloseTo(0.5, 14);
    });
});

describe('IntpBSplineUniform1', () => {
    it('throws when the degree is too large for the control count', () => {
        const controls = new NumberControls([3], [1, 2, 3]);
        expect(() => new IntpBSplineUniform1(2, controls, 0, NO_CACHING))
            .toThrow('Incompatible degree and number of controls.');
        expect(() => new IntpBSplineUniform1(3, controls, 0, NO_CACHING))
            .toThrow('Incompatible degree and number of controls.');
        // degree + 1 < numControls is acceptable.
        expect(() => new IntpBSplineUniform1(1, controls, 0, NO_CACHING))
            .not.toThrow();
    });

    it('provides member access', () => {
        const controls = new NumberControls([6], [1, 2, 3, 4, 5, 6]);
        const interp = new IntpBSplineUniform1(3, controls, 0, PRE_CACHING);
        expect(interp.getDegree(0)).toBe(3);
        expect(interp.getNumControls(0)).toBe(6);
        expect(interp.getTMin(0)).toBe(-0.5);
        expect(interp.getTMax(0)).toBe(5.5);
        expect(interp.getCacheMode()).toBe(PRE_CACHING);
    });

    it('matches the standard basis functions for degrees 1, 2 and 3', () => {
        const data = [1, 2, 0, -1, 3, 2.5, -4, 1];
        for (const degree of [1, 2, 3]) {
            for (const mode of ALL_MODES) {
                const controls = new NumberControls([data.length], data);
                const interp = new IntpBSplineUniform1(degree, controls, 0, mode);
                for (let k = 0; k <= 40; ++k) {
                    const t = -0.5 + (k / 40) * data.length;
                    expect(interp.evaluate([0], [t]))
                        .toBeCloseTo(reference1(degree, data, t), 12);
                }
            }
        }
    });

    it('is a partition of unity', () => {
        const n = 7;
        const data = new Array<number>(n).fill(1);
        for (const degree of [1, 2, 3, 4]) {
            for (const mode of ALL_MODES) {
                const controls = new NumberControls([n], data);
                const interp = new IntpBSplineUniform1(degree, controls, 0, mode);
                for (let k = 0; k <= 30; ++k) {
                    const t = -0.5 + (k / 30) * n;
                    expect(interp.evaluate([0], [t])).toBeCloseTo(1, 12);
                }
            }
        }
    });

    it('reproduces linear control data', () => {
        // For controls C[m] = a + b*m, the uniform B-spline of degree d
        // satisfies value(t) = a + b*(s - d + (d-1)/2) where
        // s - d = dsdt*(t - tmin), and the t-derivative is b*dsdt.
        const n = 9;
        const a = -2, b = 0.75;
        const data: number[] = [];
        for (let m = 0; m < n; ++m) {
            data.push(a + b * m);
        }
        for (const degree of [1, 2, 3, 4]) {
            const tmin = -0.5, tmax = n - 0.5;
            const dsdt = (n - degree) / (tmax - tmin);
            for (const mode of ALL_MODES) {
                const controls = new NumberControls([n], data);
                const interp = new IntpBSplineUniform1(degree, controls, 0, mode);
                for (let k = 1; k < 30; ++k) {
                    const t = tmin + (k / 30) * (tmax - tmin);
                    const smd = dsdt * (t - tmin);
                    expect(interp.evaluate([0], [t]))
                        .toBeCloseTo(a + b * (smd + (degree - 1) / 2), 11);
                    expect(interp.evaluate([1], [t])).toBeCloseTo(b * dsdt, 11);
                }
            }
        }
    });

    it('computes derivatives that match finite differences', () => {
        const data = [1, 2, 0, -1, 3, 2.5, -4, 1];
        const controls = new NumberControls([data.length], data);
        const interp = new IntpBSplineUniform1(3, controls, 0, NO_CACHING);
        const h = 1e-6;
        for (const t of [0.3, 1.9, 3.1, 5.5, 6.2]) {
            const d1 = (interp.evaluate([0], [t + h]) - interp.evaluate([0], [t - h]))
                / (2 * h);
            expect(interp.evaluate([1], [t])).toBeCloseTo(d1, 6);
            const d2 = (interp.evaluate([1], [t + h]) - interp.evaluate([1], [t - h]))
                / (2 * h);
            expect(interp.evaluate([2], [t])).toBeCloseTo(d2, 5);
        }
    });

    it('returns the zero control point for out-of-range derivative orders', () => {
        const data = [1, 2, 0, -1, 3];
        for (const mode of ALL_MODES) {
            const controls = new NumberControls([data.length], data);
            const interp = new IntpBSplineUniform1(2, controls, 0, mode);
            expect(interp.evaluate([3], [1])).toBe(0);
            expect(interp.evaluate([-1], [1])).toBe(0);
        }
    });

    it('clamps the parameter to the domain', () => {
        const data = [1, 2, 0, -1, 3];
        const controls = new NumberControls([data.length], data);
        const interp = new IntpBSplineUniform1(2, controls, 0, NO_CACHING);
        expect(interp.evaluate([0], [-100]))
            .toBeCloseTo(interp.evaluate([0], [-0.5]), 14);
        expect(interp.evaluate([0], [100]))
            .toBeCloseTo(interp.evaluate([0], [4.5]), 14);
        // At the ends the degree-2 spline is the average of the two
        // outermost controls.
        expect(interp.evaluate([0], [-100])).toBeCloseTo(0.5 * (data[0] + data[1]), 13);
        expect(interp.evaluate([0], [100])).toBeCloseTo(0.5 * (data[3] + data[4]), 13);
    });

    it('supports a non-numeric control-point type', () => {
        const points = [
            Vector.fromArray([0, 0]),
            Vector.fromArray([1, 2]),
            Vector.fromArray([3, -1]),
            Vector.fromArray([4, 4]),
            Vector.fromArray([6, 0])
        ];
        const controls = new VectorControls(points.length, points);
        const interp = new IntpBSplineUniform1(2, controls,
            Vector.zero(2), NO_CACHING);
        const xData = points.map(p => p.get(0));
        const yData = points.map(p => p.get(1));
        for (let k = 0; k <= 20; ++k) {
            const t = -0.5 + (k / 20) * points.length;
            const P = interp.evaluate([0], [t]);
            expect(P.get(0)).toBeCloseTo(reference1(2, xData, t), 12);
            expect(P.get(1)).toBeCloseTo(reference1(2, yData, t), 12);
        }
    });

    it('agrees with the general-dimension interpolator', () => {
        const data = [1, 2, 0, -1, 3, 2.5, -4, 1];
        for (const degree of [1, 2, 3]) {
            for (const mode of ALL_MODES) {
                const specialized = new IntpBSplineUniform1(degree,
                    new NumberControls([data.length], data), 0, mode);
                const general = new IntpBSplineUniform([degree],
                    new NumberControls([data.length], data), 0, mode);
                for (let k = 0; k <= 25; ++k) {
                    const t = -0.5 + (k / 25) * data.length;
                    for (let order = 0; order <= degree; ++order) {
                        expect(general.evaluate([order], [t]))
                            .toBeCloseTo(specialized.evaluate([order], [t]), 11);
                    }
                }
            }
        }
    });
});

describe('IntpBSplineUniform2', () => {
    it('throws when a degree is too large for the control count', () => {
        const controls = new NumberControls([5, 3], new Array<number>(15).fill(0));
        expect(() => new IntpBSplineUniform2([3, 2], controls, 0, NO_CACHING))
            .toThrow('Incompatible degree and number of controls.');
    });

    it('provides member access', () => {
        const controls = new NumberControls([5, 6], new Array<number>(30).fill(0));
        const interp = new IntpBSplineUniform2([2, 3], controls, 0, NO_CACHING);
        expect(interp.getDegree(0)).toBe(2);
        expect(interp.getDegree(1)).toBe(3);
        expect(interp.getNumControls(0)).toBe(5);
        expect(interp.getNumControls(1)).toBe(6);
        expect(interp.getTMin(1)).toBe(-0.5);
        expect(interp.getTMax(1)).toBe(5.5);
        expect(interp.getCacheMode()).toBe(NO_CACHING);
    });

    it('is the tensor product of two 1-dimensional splines', () => {
        // With C(i,j) = f[i]*g[j] the 2-dimensional interpolator factors
        // into the product of the 1-dimensional interpolators.
        const f = [1, 2, 0, -1, 3, 2.5];
        const g = [-2, 0.5, 1, 4, 0];
        const data: number[] = [];
        for (let j = 0; j < g.length; ++j) {
            for (let i = 0; i < f.length; ++i) {
                data.push(f[i] * g[j]);
            }
        }
        const degrees = [2, 3];
        for (const mode of ALL_MODES) {
            const interp = new IntpBSplineUniform2(degrees,
                new NumberControls([f.length, g.length], data), 0, mode);
            const i0 = new IntpBSplineUniform1(degrees[0],
                new NumberControls([f.length], f), 0, NO_CACHING);
            const i1 = new IntpBSplineUniform1(degrees[1],
                new NumberControls([g.length], g), 0, NO_CACHING);
            for (let a = 0; a <= 8; ++a) {
                for (let b = 0; b <= 8; ++b) {
                    const t0 = -0.5 + (a / 8) * f.length;
                    const t1 = -0.5 + (b / 8) * g.length;
                    for (const order of [[0, 0], [1, 0], [0, 1], [1, 1], [2, 1]]) {
                        const expected = i0.evaluate([order[0]], [t0])
                            * i1.evaluate([order[1]], [t1]);
                        expect(interp.evaluate(order, [t0, t1]))
                            .toBeCloseTo(expected, 10);
                    }
                }
            }
        }
    });

    it('is a partition of unity', () => {
        const data = new Array<number>(5 * 6).fill(1);
        for (const mode of ALL_MODES) {
            const interp = new IntpBSplineUniform2([2, 3],
                new NumberControls([5, 6], data), 0, mode);
            for (let a = 0; a <= 6; ++a) {
                for (let b = 0; b <= 6; ++b) {
                    const t0 = -0.5 + (a / 6) * 5;
                    const t1 = -0.5 + (b / 6) * 6;
                    expect(interp.evaluate([0, 0], [t0, t1])).toBeCloseTo(1, 12);
                }
            }
        }
    });

    it('returns the zero control point for out-of-range derivative orders', () => {
        const data = new Array<number>(5 * 6).fill(1);
        const interp = new IntpBSplineUniform2([2, 3],
            new NumberControls([5, 6], data), 0, NO_CACHING);
        expect(interp.evaluate([3, 0], [1, 1])).toBe(0);
        expect(interp.evaluate([0, 4], [1, 1])).toBe(0);
        expect(interp.evaluate([0, -1], [1, 1])).toBe(0);
    });

    it('agrees with the general-dimension interpolator', () => {
        const sizes = [5, 6];
        const data: number[] = [];
        let seed = 12345;
        const random = () => {
            seed = (seed * 16807) % 2147483647;
            return seed / 2147483647 - 0.5;
        };
        for (let i = 0; i < sizes[0] * sizes[1]; ++i) {
            data.push(random());
        }
        const degrees = [2, 3];
        for (const mode of ALL_MODES) {
            const specialized = new IntpBSplineUniform2(degrees,
                new NumberControls(sizes, data), 0, mode);
            const general = new IntpBSplineUniform(degrees,
                new NumberControls(sizes, data), 0, mode);
            for (let k = 0; k < 40; ++k) {
                const t0 = -0.5 + Math.abs(random() + 0.5) * sizes[0];
                const t1 = -0.5 + Math.abs(random() + 0.5) * sizes[1];
                for (const order of [[0, 0], [1, 0], [0, 2], [2, 3]]) {
                    expect(general.evaluate(order, [t0, t1]))
                        .toBeCloseTo(specialized.evaluate(order, [t0, t1]), 11);
                }
            }
        }
    });
});

describe('IntpBSplineUniform3', () => {
    it('throws when a degree is too large for the control count', () => {
        const controls = new NumberControls([5, 5, 3], new Array<number>(75).fill(0));
        expect(() => new IntpBSplineUniform3([2, 2, 2], controls, 0, NO_CACHING))
            .toThrow('Incompatible degree and number of controls.');
    });

    it('is the tensor product of three 1-dimensional splines', () => {
        const f = [1, 2, 0, -1, 3];
        const g = [-2, 0.5, 1, 4];
        const h = [0.25, -1, 2, 3, 1];
        const sizes = [f.length, g.length, h.length];
        const data: number[] = [];
        for (let k = 0; k < h.length; ++k) {
            for (let j = 0; j < g.length; ++j) {
                for (let i = 0; i < f.length; ++i) {
                    data.push(f[i] * g[j] * h[k]);
                }
            }
        }
        const degrees = [2, 1, 3];
        for (const mode of ALL_MODES) {
            const interp = new IntpBSplineUniform3(degrees,
                new NumberControls(sizes, data), 0, mode);
            const i0 = new IntpBSplineUniform1(degrees[0],
                new NumberControls([f.length], f), 0, NO_CACHING);
            const i1 = new IntpBSplineUniform1(degrees[1],
                new NumberControls([g.length], g), 0, NO_CACHING);
            const i2 = new IntpBSplineUniform1(degrees[2],
                new NumberControls([h.length], h), 0, NO_CACHING);
            for (let a = 0; a <= 4; ++a) {
                for (let b = 0; b <= 4; ++b) {
                    for (let c = 0; c <= 4; ++c) {
                        const t = [-0.5 + (a / 4) * sizes[0],
                            -0.5 + (b / 4) * sizes[1],
                            -0.5 + (c / 4) * sizes[2]];
                        for (const order of [[0, 0, 0], [1, 0, 1], [2, 1, 3]]) {
                            const expected = i0.evaluate([order[0]], [t[0]])
                                * i1.evaluate([order[1]], [t[1]])
                                * i2.evaluate([order[2]], [t[2]]);
                            expect(interp.evaluate(order, t)).toBeCloseTo(expected, 10);
                        }
                    }
                }
            }
        }
    });

    it('agrees with the general-dimension interpolator', () => {
        const sizes = [5, 4, 6];
        const data: number[] = [];
        let seed = 777;
        const random = () => {
            seed = (seed * 16807) % 2147483647;
            return seed / 2147483647 - 0.5;
        };
        for (let i = 0; i < sizes[0] * sizes[1] * sizes[2]; ++i) {
            data.push(random());
        }
        const degrees = [2, 1, 3];
        for (const mode of ALL_MODES) {
            const specialized = new IntpBSplineUniform3(degrees,
                new NumberControls(sizes, data), 0, mode);
            const general = new IntpBSplineUniform(degrees,
                new NumberControls(sizes, data), 0, mode);
            for (let k = 0; k < 25; ++k) {
                const t = [
                    -0.5 + Math.abs(random() + 0.5) * sizes[0],
                    -0.5 + Math.abs(random() + 0.5) * sizes[1],
                    -0.5 + Math.abs(random() + 0.5) * sizes[2]
                ];
                for (const order of [[0, 0, 0], [1, 1, 1], [2, 0, 2]]) {
                    expect(general.evaluate(order, t))
                        .toBeCloseTo(specialized.evaluate(order, t), 11);
                }
            }
        }
    });

    it('returns the zero control point for out-of-range derivative orders', () => {
        const sizes = [5, 4, 6];
        const data = new Array<number>(5 * 4 * 6).fill(1);
        const interp = new IntpBSplineUniform3([2, 1, 3],
            new NumberControls(sizes, data), 0, NO_CACHING);
        expect(interp.evaluate([3, 0, 0], [0, 0, 0])).toBe(0);
        expect(interp.evaluate([0, 2, 0], [0, 0, 0])).toBe(0);
        expect(interp.evaluate([0, 0, -1], [0, 0, 0])).toBe(0);
    });
});

describe('IntpBSplineUniform (general dimension)', () => {
    it('throws when a degree is too large for the control count', () => {
        const controls = new NumberControls([4], [1, 2, 3, 4]);
        expect(() => new IntpBSplineUniform([3], controls, 0, NO_CACHING))
            .toThrow('Incompatible degree and number of controls.');
    });

    it('returns the zero control point when the inputs are too short', () => {
        const controls = new NumberControls([5, 5], new Array<number>(25).fill(1));
        const interp = new IntpBSplineUniform([2, 2], controls, 0, NO_CACHING);
        expect(interp.evaluate([0], [0, 0])).toBe(0);
        expect(interp.evaluate([0, 0], [0])).toBe(0);
        expect(interp.evaluate([0, 0], [0, 0])).toBeCloseTo(1, 12);
    });

    it('provides member access', () => {
        const controls = new NumberControls([5, 6], new Array<number>(30).fill(0));
        const interp = new IntpBSplineUniform([2, 3], controls, 0, ON_DEMAND_CACHING);
        expect(interp.getDegree(1)).toBe(3);
        expect(interp.getNumControls(0)).toBe(5);
        expect(interp.getTMin(0)).toBe(-0.5);
        expect(interp.getTMax(0)).toBe(4.5);
        expect(interp.getCacheMode()).toBe(ON_DEMAND_CACHING);
    });

    it('agrees across the three cache modes in four dimensions', () => {
        const sizes = [4, 4, 4, 4];
        const degrees = [1, 2, 1, 2];
        let seed = 4242;
        const random = () => {
            seed = (seed * 16807) % 2147483647;
            return seed / 2147483647 - 0.5;
        };
        const data: number[] = [];
        for (let i = 0; i < 4 * 4 * 4 * 4; ++i) {
            data.push(random());
        }
        const interps = ALL_MODES.map(mode => new IntpBSplineUniform(degrees,
            new NumberControls(sizes, data), 0, mode));
        for (let k = 0; k < 20; ++k) {
            const t = [
                -0.5 + Math.abs(random() + 0.5) * 4,
                -0.5 + Math.abs(random() + 0.5) * 4,
                -0.5 + Math.abs(random() + 0.5) * 4,
                -0.5 + Math.abs(random() + 0.5) * 4
            ];
            for (const order of [[0, 0, 0, 0], [1, 0, 1, 0], [1, 2, 1, 2]]) {
                const v0 = interps[0].evaluate(order, t);
                expect(interps[1].evaluate(order, t)).toBeCloseTo(v0, 11);
                expect(interps[2].evaluate(order, t)).toBeCloseTo(v0, 11);
            }
        }
    });

    it('is a partition of unity in four dimensions', () => {
        const sizes = [4, 5, 4, 5];
        const degrees = [1, 2, 2, 3];
        const data = new Array<number>(4 * 5 * 4 * 5).fill(1);
        for (const mode of ALL_MODES) {
            const interp = new IntpBSplineUniform(degrees,
                new NumberControls(sizes, data), 0, mode);
            for (let k = 0; k <= 5; ++k) {
                const t = sizes.map(s => -0.5 + (k / 5) * s);
                expect(interp.evaluate([0, 0, 0, 0], t)).toBeCloseTo(1, 11);
            }
        }
    });

    it('supports degree 0, for which the spline is piecewise constant', () => {
        // Degree 0 exercises the special case in the blending-matrix
        // computation; the value is the control point of the containing
        // interval.
        const data = [1, 2, 0, -1, 3];
        const interp = new IntpBSplineUniform([0],
            new NumberControls([data.length], data), 0, NO_CACHING);
        expect(interp.getTMin(0)).toBe(-0.5);
        expect(interp.getTMax(0)).toBe(4.5);
        for (let m = 0; m < data.length; ++m) {
            // ds/dt = 5/5 = 1, so the interval of index m is
            // [-0.5 + m, 0.5 + m).
            expect(interp.evaluate([0], [m])).toBeCloseTo(data[m], 14);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against upstream
// IntpBSplineUniform.h.
// ---------------------------------------------------------------------------

const degreeArb = fc.integer({ min: 1, max: 3 });
const controlArb = fc.integer({ min: -20, max: 20 });

describe('IntpBSplineUniform verification', () => {
    // Row r of the blending matrix holds the coefficients of the local basis
    // polynomial Q_{d,d-r}(s). The B-spline basis is a partition of unity, so
    // the row polynomials sum to one for every s in [0,1], and each of them is
    // nonnegative there.
    it('produces a blending matrix whose rows are a partition of unity', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 5 }), finite(0, 1)),
            ([degree, u]) => {
                const A = IntpBSplineUniformShared.computeBlendingMatrix(degree);
                const degreeP1 = degree + 1;
                expect(A.length).toBe(degreeP1 * degreeP1);
                let sum = 0;
                for (let row = 0; row < degreeP1; ++row) {
                    let value = 0;
                    for (let col = degree; col >= 0; --col) {
                        value = value * u + A[col + degreeP1 * row];
                    }
                    expect(value).toBeGreaterThanOrEqual(-1e-12);
                    sum += value;
                }
                expectClose(sum, 1, 1e-12, 1e-12);
                return true;
            });
    });

    // ComputePowers must fill powerDSDT[i] with (ds/dt)^i. Degree 0 is the
    // case upstream overruns its buffer (issue #135); the port relies on the
    // JavaScript array growing and is exercised by the degree-0 spline below.
    it('produces the powers of ds/dt', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 5 }),
            fc.integer({ min: 6, max: 12 })), ([degree, numControls]) => {
                const tmin = -0.5, tmax = numControls - 0.5;
                const powers = IntpBSplineUniformShared.computePowers(degree,
                    numControls, tmin, tmax);
                const dsdt = (numControls - degree) / (tmax - tmin);
                expect(powers[0]).toBe(1);
                for (let i = 1; i <= degree; ++i) {
                    expectClose(powers[i], Math.pow(dsdt, i), 1e-12, 1e-12);
                }
                return true;
            });
    });

    // GetKey maps t to the interval index and the local parameter. Below tmin
    // it pins (0, 0); at or above tmax it pins the last interval with u = 1,
    // which extends the s-domain [d, c+1) to its support [d, c+1].
    it('produces a key with 0 <= u <= 1 in the correct interval', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 4 }),
            fc.integer({ min: 6, max: 12 }), finite(-20, 20)),
            ([degree, numControls, t]) => {
                const tmin = -0.5, tmax = numControls - 0.5;
                const dsdt = (numControls - degree) / (tmax - tmin);
                const { index, u } = IntpBSplineUniformShared.getKey(t, tmin,
                    tmax, dsdt, numControls, degree);
                expect(index).toBeGreaterThanOrEqual(0);
                expect(index).toBeLessThanOrEqual(numControls - 1 - degree);
                expect(u).toBeGreaterThanOrEqual(0);
                expect(u).toBeLessThanOrEqual(1);
                if (t <= tmin) {
                    expect(index).toBe(0);
                    expect(u).toBe(0);
                } else if (t >= tmax) {
                    expect(index).toBe(numControls - 1 - degree);
                    expect(u).toBe(1);
                } else {
                    // index + u is the continuous s - d coordinate.
                    expectClose(index + u, dsdt * (t - tmin), 1e-9, 1e-12);
                }
                return true;
            });
    });

    // The three cache modes evaluate the same spline. The two caching modes
    // share ComputeTensor and the same accumulation order, so they must agree
    // bit for bit; the non-caching mode sums the same terms in a different
    // order, so it agrees to round-off.
    it('agrees across the three cache modes in one dimension', () => {
        check(fc.tuple(degreeArb,
            fc.array(controlArb, { minLength: 6, maxLength: 10 }),
            fc.integer({ min: 0, max: 3 }), finite(-2, 12)),
            ([degree, data, order, t]) => {
                if (order > degree) {
                    return true;
                }
                const controls = new NumberControls([data.length], data);
                const make = (mode: IntpBSplineUniformCacheMode) =>
                    new IntpBSplineUniform1<number>(degree, controls, 0, mode);
                const none = make(NO_CACHING).evaluate([order], [t]);
                const pre = make(PRE_CACHING).evaluate([order], [t]);
                const demand = make(ON_DEMAND_CACHING).evaluate([order], [t]);
                expect(demand).toBe(pre);
                expectClose(none, pre, 1e-9, 1e-11);
                return true;
            });
    });

    // The hand-optimized one-dimensional specialization and the
    // general-dimension implementation must evaluate the same spline.
    it('the 1D specialization agrees with the general-dimension class', () => {
        check(fc.tuple(degreeArb,
            fc.array(controlArb, { minLength: 6, maxLength: 10 }),
            fc.integer({ min: 0, max: 3 }), finite(-2, 12),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([degree, data, order, t, mode]) => {
                if (order > degree) {
                    return true;
                }
                const controls = new NumberControls([data.length], data);
                const special = new IntpBSplineUniform1<number>(degree, controls,
                    0, mode);
                const general = new IntpBSplineUniform<number>([degree], controls,
                    0, mode);
                expectClose(special.evaluate([order], [t]),
                    general.evaluate([order], [t]), 1e-9, 1e-11);
                return true;
            });
    });

    // The two-dimensional spline is the tensor product of the two
    // one-dimensional splines, so it must agree with the general-dimension
    // class on separable control data as well as on arbitrary data.
    it('the 2D specialization agrees with the general-dimension class', () => {
        check(fc.tuple(degreeArb, degreeArb,
            fc.integer({ min: 5, max: 7 }), fc.integer({ min: 5, max: 7 }),
            fc.integer({ min: 0, max: 0xffff }), finite(-1, 8), finite(-1, 8),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([d0, d1, n0, n1, seed, t0, t1, mode]) => {
                const rand = seededRandom(seed + 1);
                const data: number[] = [];
                for (let i = 0; i < n0 * n1; ++i) {
                    data.push(Math.round(20 * (2 * rand() - 1)));
                }
                const controls = new NumberControls([n0, n1], data);
                const special = new IntpBSplineUniform2<number>([d0, d1],
                    controls, 0, mode);
                const general = new IntpBSplineUniform<number>([d0, d1],
                    controls, 0, mode);
                for (let o0 = 0; o0 <= d0 && o0 <= 1; ++o0) {
                    for (let o1 = 0; o1 <= d1 && o1 <= 1; ++o1) {
                        expectClose(special.evaluate([o0, o1], [t0, t1]),
                            general.evaluate([o0, o1], [t0, t1]), 1e-8, 1e-10);
                    }
                }
                return true;
            }, 60);
    }, 30000);

    it('the 3D specialization agrees with the general-dimension class', () => {
        check(fc.tuple(degreeArb, degreeArb, degreeArb,
            fc.integer({ min: 0, max: 0xffff }),
            finite(-1, 6), finite(-1, 6), finite(-1, 6),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([d0, d1, d2, seed, t0, t1, t2, mode]) => {
                const n = 5;
                const rand = seededRandom(seed + 1);
                const data: number[] = [];
                for (let i = 0; i < n * n * n; ++i) {
                    data.push(Math.round(20 * (2 * rand() - 1)));
                }
                const controls = new NumberControls([n, n, n], data);
                const special = new IntpBSplineUniform3<number>([d0, d1, d2],
                    controls, 0, mode);
                const general = new IntpBSplineUniform<number>([d0, d1, d2],
                    controls, 0, mode);
                expectClose(special.evaluate([0, 0, 0], [t0, t1, t2]),
                    general.evaluate([0, 0, 0], [t0, t1, t2]), 1e-8, 1e-10);
                expectClose(special.evaluate([1, 0, 0], [t0, t1, t2]),
                    general.evaluate([1, 0, 0], [t0, t1, t2]), 1e-7, 1e-9);
                return true;
            }, 40);
    }, 30000);

    // The one-dimensional spline must agree with the standard uniform
    // B-spline basis, written out independently above.
    it('agrees with the standard uniform B-spline basis', () => {
        check(fc.tuple(degreeArb,
            fc.array(controlArb, { minLength: 6, maxLength: 10 }),
            finite(-2, 12),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([degree, data, t, mode]) => {
                const controls = new NumberControls([data.length], data);
                const spline = new IntpBSplineUniform1<number>(degree, controls,
                    0, mode);
                expectClose(spline.evaluate([0], [t]),
                    reference1(degree, data, t), 1e-9, 1e-11);
                return true;
            });
    });

    // The basis is a partition of unity in every dimension, so a constant
    // control field is reproduced everywhere and its derivative vanishes.
    it('reproduces a constant control field in two dimensions', () => {
        check(fc.tuple(degreeArb, degreeArb, fc.integer({ min: -9, max: 9 }),
            finite(-2, 8), finite(-2, 8),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([d0, d1, value, t0, t1, mode]) => {
                const n = 6;
                const controls = new NumberControls([n, n],
                    new Array<number>(n * n).fill(value));
                const spline = new IntpBSplineUniform2<number>([d0, d1],
                    controls, 0, mode);
                expectClose(spline.evaluate([0, 0], [t0, t1]), value,
                    1e-9, 1e-11);
                expectClose(spline.evaluate([1, 0], [t0, t1]), 0, 1e-8, 1e-9);
                return true;
            });
    });

    // The derivative accessor must agree with numerical differentiation. The
    // spline is a polynomial of degree at most three on each interval, so the
    // fourth-order five-point stencil is exact for it up to round-off; the
    // sample point is placed at the middle of the domain and the step is a
    // small fraction of one interval so the stencil stays inside it.
    it('computes derivatives that match a fourth-order stencil', () => {
        check(fc.tuple(degreeArb,
            fc.array(controlArb, { minLength: 8, maxLength: 10 }),
            fc.integer({ min: 1, max: 3 }).map(k => k / 4),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([degree, data, frac, mode]) => {
                const controls = new NumberControls([data.length], data);
                const spline = new IntpBSplineUniform1<number>(degree, controls,
                    0, mode);
                const tmin = spline.getTMin(0), tmax = spline.getTMax(0);
                const numIntervals = data.length - degree;
                const width = (tmax - tmin) / numIntervals;
                // A point strictly inside interval 2 of the s-parameter.
                const t = tmin + width * (2 + frac);
                const h = width / 16;
                const f = (x: number) => spline.evaluate([0], [x]);
                const d1 = (-f(t + 2 * h) + 8 * f(t + h) - 8 * f(t - h)
                    + f(t - 2 * h)) / (12 * h);
                expectClose(spline.evaluate([1], [t]), d1, 1e-6, 1e-7);
                if (degree >= 2) {
                    const d2 = (-f(t + 2 * h) + 16 * f(t + h) - 30 * f(t)
                        + 16 * f(t - h) - f(t - 2 * h)) / (12 * h * h);
                    expectClose(spline.evaluate([2], [t]), d2, 1e-4, 1e-5);
                }
                return true;
            });
    });

    // Degree 0 gives a piecewise-constant spline. Upstream's ComputePowers
    // writes powerDSDT[1] into a buffer sized degree+1 = 1 (issue #135); in
    // the port the array simply grows and degree 0 evaluates correctly.
    it('supports degree 0 in every cache mode', () => {
        check(fc.tuple(fc.array(controlArb, { minLength: 4, maxLength: 8 }),
            finite(-2, 10),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([data, t, mode]) => {
                const controls = new NumberControls([data.length], data);
                const spline = new IntpBSplineUniform1<number>(0, controls, 0,
                    mode);
                const value = spline.evaluate([0], [t]);
                // The value is one of the control points: the degree-0 basis
                // selects the control point of the containing interval.
                expect(data).toContain(value);
                const { index } = IntpBSplineUniformShared.getKey(t, -0.5,
                    data.length - 0.5, 1, data.length, 0);
                expect(value).toBe(data[index]);
                return true;
            });
    });

    // Upstream returns the zero control point when any order is negative or
    // exceeds the degree of its dimension.
    it('returns the zero control point for orders outside [0, degree]', () => {
        check(fc.tuple(degreeArb,
            fc.array(controlArb, { minLength: 6, maxLength: 8 }),
            fc.integer({ min: 4, max: 8 }), finite(-2, 10),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([degree, data, order, t, mode]) => {
                const controls = new NumberControls([data.length], data);
                const spline = new IntpBSplineUniform1<number>(degree, controls,
                    0, mode);
                expect(spline.evaluate([order], [t])).toBe(0);
                const general = new IntpBSplineUniform<number>([degree],
                    controls, 0, mode);
                expect(general.evaluate([order], [t])).toBe(0);
                expect(general.evaluate([-1], [t])).toBe(0);
                return true;
            });
    });

    // The domain of the interpolator is [-1/2, c - 1/2] and the evaluators pin
    // the key outside it, so the spline is constant beyond either end.
    it('is constant outside the parameter domain', () => {
        check(fc.tuple(degreeArb,
            fc.array(controlArb, { minLength: 6, maxLength: 9 }),
            finite(0.01, 50),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([degree, data, delta, mode]) => {
                const controls = new NumberControls([data.length], data);
                const spline = new IntpBSplineUniform1<number>(degree, controls,
                    0, mode);
                const tmin = spline.getTMin(0), tmax = spline.getTMax(0);
                expect(spline.evaluate([0], [tmin - delta]))
                    .toBe(spline.evaluate([0], [tmin]));
                expect(spline.evaluate([0], [tmax + delta]))
                    .toBe(spline.evaluate([0], [tmax]));
                return true;
            });
    });

    // The evaluator is linear in the control points: it is a fixed weighted
    // sum of them for a given (order, t).
    it('is linear in the control points', () => {
        check(fc.tuple(degreeArb,
            fc.array(controlArb, { minLength: 7, maxLength: 9 }),
            fc.array(controlArb, { minLength: 9, maxLength: 9 }),
            fc.integer({ min: -4, max: 4 }), fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: 0, max: 1 }), finite(-1, 9),
            fc.constantFrom(NO_CACHING, PRE_CACHING, ON_DEMAND_CACHING)),
            ([degree, d1, d2raw, a, b, order, t, mode]) => {
                if (order > degree) {
                    return true;
                }
                const d2 = d2raw.slice(0, d1.length);
                const dc = d1.map((v, i) => a * v + b * d2[i]);
                const make = (data: number[]) =>
                    new IntpBSplineUniform1<number>(degree,
                        new NumberControls([data.length], data), 0, mode);
                const v1 = make(d1).evaluate([order], [t]);
                const v2 = make(d2).evaluate([order], [t]);
                const vc = make(dc).evaluate([order], [t]);
                expectClose(a * v1 + b * v2, vc, 1e-8, 1e-10);
                return true;
            });
    });
});
