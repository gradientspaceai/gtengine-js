import { describe, it, expect } from 'vitest';
import {
    IntpBSplineUniform, IntpBSplineUniform1, IntpBSplineUniform2,
    IntpBSplineUniform3, IntpBSplineUniformCacheMode, IntpBSplineUniformShared
} from '../src/IntpBSplineUniform.js';
import type { IntpBSplineUniformControls } from '../src/IntpBSplineUniform.js';
import { Vector } from '../src/Vector.js';
import { add as vadd, mul as vmul } from '../src/Vector.js';

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
