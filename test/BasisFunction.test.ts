import { describe, it, expect } from 'vitest';
import { BasisFunction, BasisFunctionInput, UniqueKnot } from '../src/BasisFunction.js';

// Evaluate the single basis function N_{i,degree} at t (order 'order'),
// honoring local support: indices outside [minIndex, maxIndex] are zero.
function basisValue(bf: BasisFunction, order: number, i: number, t: number): number {
    const { minIndex, maxIndex } = bf.evaluate(t, order);
    return (minIndex <= i && i <= maxIndex) ? bf.getValue(order, i) : 0;
}

// A nonuniform floating (non-open) input: all knots have multiplicity 1.
function makeFloatingInput(numControls: number, degree: number, knots: number[]): BasisFunctionInput {
    const input = new BasisFunctionInput();
    input.numControls = numControls;
    input.degree = degree;
    input.uniform = false;
    input.periodic = false;
    input.numUniqueKnots = knots.length;
    input.uniqueKnots = knots.map((t) => new UniqueKnot(t, 1));
    return input;
}

// A uniform periodic input with integer knots 0, 1, ..., n + 2 * degree.
function makePeriodicInput(numControls: number, degree: number): BasisFunctionInput {
    const input = new BasisFunctionInput();
    input.numControls = numControls;
    input.degree = degree;
    input.uniform = true;
    input.periodic = true;
    input.numUniqueKnots = numControls + 2 * degree + 1;
    input.uniqueKnots = [];
    for (let i = 0; i < input.numUniqueKnots; ++i) {
        input.uniqueKnots.push(new UniqueKnot(i, 1));
    }
    return input;
}

describe('BasisFunctionInput', () => {
    it('default constructor zero-initializes', () => {
        const input = new BasisFunctionInput();
        expect(input.numControls).toBe(0);
        expect(input.degree).toBe(0);
        expect(input.uniform).toBe(false);
        expect(input.periodic).toBe(false);
        expect(input.numUniqueKnots).toBe(0);
        expect(input.uniqueKnots).toEqual([]);
    });

    it('constructs an open uniform configuration on [0,1]', () => {
        const input = new BasisFunctionInput(6, 2);
        expect(input.numControls).toBe(6);
        expect(input.degree).toBe(2);
        expect(input.uniform).toBe(true);
        expect(input.periodic).toBe(false);
        expect(input.numUniqueKnots).toBe(5);
        expect(input.uniqueKnots[0].t).toBe(0);
        expect(input.uniqueKnots[0].multiplicity).toBe(3);
        for (let i = 1; i <= 3; ++i) {
            expect(input.uniqueKnots[i].t).toBeCloseTo(i / 4, 15);
            expect(input.uniqueKnots[i].multiplicity).toBe(1);
        }
        expect(input.uniqueKnots[4].t).toBe(1);
        expect(input.uniqueKnots[4].multiplicity).toBe(3);
    });
});

describe('BasisFunction creation and member access', () => {
    it('builds the expanded knot vector for an open uniform input', () => {
        const bf = new BasisFunction(new BasisFunctionInput(7, 3));
        expect(bf.getNumControls()).toBe(7);
        expect(bf.getDegree()).toBe(3);
        expect(bf.getNumUniqueKnots()).toBe(5);
        expect(bf.getNumKnots()).toBe(11); // n + d + 1
        expect(bf.getMinDomain()).toBe(0);
        expect(bf.getMaxDomain()).toBe(1);
        expect(bf.isOpen()).toBe(true);
        expect(bf.isUniform()).toBe(true);
        expect(bf.isPeriodic()).toBe(false);

        const knots = bf.getKnots();
        const expected = [0, 0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1, 1];
        expect(knots.length).toBe(expected.length);
        for (let i = 0; i < expected.length; ++i) {
            expect(knots[i]).toBeCloseTo(expected[i], 15);
        }

        const unique = bf.getUniqueKnots();
        expect(unique.length).toBe(5);
        expect(unique[0].multiplicity).toBe(4);
        expect(unique[4].multiplicity).toBe(4);
    });

    it('default constructor creates an empty object that can be created later', () => {
        const bf = new BasisFunction();
        expect(bf.getNumControls()).toBe(0);
        bf.create(new BasisFunctionInput(4, 2));
        expect(bf.getNumControls()).toBe(4);
        expect(() => bf.create(new BasisFunctionInput(4, 2))).toThrow('Object already created.');
    });

    it('classifies a floating (non-open) knot vector', () => {
        const input = makeFloatingInput(5, 2, [0, 1, 2.5, 3, 4.5, 5, 6, 7]);
        const bf = new BasisFunction(input);
        expect(bf.isOpen()).toBe(false);
        expect(bf.getMinDomain()).toBe(2.5); // knots[degree]
        expect(bf.getMaxDomain()).toBe(5); // knots[numControls]
    });

    it('validates input', () => {
        expect(() => new BasisFunction(new BasisFunctionInput(1, 1)))
            .toThrow('Invalid number of control points.');
        const badDegree = new BasisFunctionInput(4, 2);
        badDegree.degree = 4;
        expect(() => new BasisFunction(badDegree)).toThrow('Invalid degree.');

        const badUnique = makeFloatingInput(3, 1, [0]);
        badUnique.uniqueKnots[0].multiplicity = 5;
        expect(() => new BasisFunction(badUnique)).toThrow('Invalid number of unique knots.');

        const notIncreasing = makeFloatingInput(5, 2, [0, 1, 1, 2, 3, 4, 5, 6]);
        expect(() => new BasisFunction(notIncreasing))
            .toThrow('Unique knots are not strictly increasing.');

        const badMult0 = makeFloatingInput(5, 2, [0, 1, 2, 3, 4, 5, 6, 7]);
        badMult0.uniqueKnots[0].multiplicity = 4;
        expect(() => new BasisFunction(badMult0)).toThrow('Invalid first multiplicity.');

        const badMult1 = makeFloatingInput(5, 2, [0, 1, 2, 3, 4, 5, 6, 7]);
        badMult1.uniqueKnots[7].multiplicity = 4;
        expect(() => new BasisFunction(badMult1)).toThrow('Invalid last multiplicity.');

        const badMultI = makeFloatingInput(5, 2, [0, 1, 2, 3, 4, 5, 6, 7]);
        badMultI.uniqueKnots[3].multiplicity = 4;
        expect(() => new BasisFunction(badMultI)).toThrow('Invalid interior multiplicity.');
    });
});

describe('BasisFunction hand-computed values', () => {
    it('degree-1 hat functions on knots [0,0,1/2,1,1]', () => {
        const input = new BasisFunctionInput();
        input.numControls = 3;
        input.degree = 1;
        input.uniform = true;
        input.periodic = false;
        input.numUniqueKnots = 3;
        input.uniqueKnots = [new UniqueKnot(0, 2), new UniqueKnot(0.5, 1), new UniqueKnot(1, 2)];
        const bf = new BasisFunction(input);
        expect(bf.isOpen()).toBe(true);

        let r = bf.evaluate(0.3, 0);
        expect(r.minIndex).toBe(0);
        expect(r.maxIndex).toBe(1);
        expect(bf.getValue(0, 0)).toBeCloseTo(0.4, 14); // (1/2 - t) / (1/2)
        expect(bf.getValue(0, 1)).toBeCloseTo(0.6, 14); // t / (1/2)

        r = bf.evaluate(0.75, 0);
        expect(r.minIndex).toBe(1);
        expect(r.maxIndex).toBe(2);
        expect(bf.getValue(0, 1)).toBeCloseTo(0.5, 14); // (1 - t) / (1/2)
        expect(bf.getValue(0, 2)).toBeCloseTo(0.5, 14); // (t - 1/2) / (1/2)
    });

    it('degree-2 open with 3 controls reproduces the Bernstein basis', () => {
        // Knot vector [0,0,0,1,1,1]: N_i = B_i^2 (quadratic Bernstein).
        const bf = new BasisFunction(new BasisFunctionInput(3, 2));
        const t = 0.3;
        const r = bf.evaluate(t, 2);
        expect(r.minIndex).toBe(0);
        expect(r.maxIndex).toBe(2);

        expect(bf.getValue(0, 0)).toBeCloseTo((1 - t) * (1 - t), 14);
        expect(bf.getValue(0, 1)).toBeCloseTo(2 * t * (1 - t), 14);
        expect(bf.getValue(0, 2)).toBeCloseTo(t * t, 14);

        expect(bf.getValue(1, 0)).toBeCloseTo(-2 * (1 - t), 14);
        expect(bf.getValue(1, 1)).toBeCloseTo(2 - 4 * t, 14);
        expect(bf.getValue(1, 2)).toBeCloseTo(2 * t, 14);

        expect(bf.getValue(2, 0)).toBeCloseTo(2, 14);
        expect(bf.getValue(2, 1)).toBeCloseTo(-4, 14);
        expect(bf.getValue(2, 2)).toBeCloseTo(2, 14);
    });

    it('degree-3 open with 4 controls reproduces the cubic Bernstein basis', () => {
        // Knot vector [0,0,0,0,1,1,1,1]: N_i = B_i^3 (cubic Bernstein).
        const bf = new BasisFunction(new BasisFunctionInput(4, 3));
        const t = 0.4;
        const s = 1 - t;
        bf.evaluate(t, 3);

        expect(bf.getValue(0, 0)).toBeCloseTo(s * s * s, 14);
        expect(bf.getValue(0, 1)).toBeCloseTo(3 * t * s * s, 14);
        expect(bf.getValue(0, 2)).toBeCloseTo(3 * t * t * s, 14);
        expect(bf.getValue(0, 3)).toBeCloseTo(t * t * t, 14);

        expect(bf.getValue(1, 0)).toBeCloseTo(-3 * s * s, 14);
        expect(bf.getValue(1, 1)).toBeCloseTo(3 * s * (1 - 3 * t), 14);
        expect(bf.getValue(1, 2)).toBeCloseTo(3 * t * (2 - 3 * t), 14);
        expect(bf.getValue(1, 3)).toBeCloseTo(3 * t * t, 14);

        expect(bf.getValue(2, 0)).toBeCloseTo(6 * s, 13);
        expect(bf.getValue(2, 1)).toBeCloseTo(18 * t - 12, 13);
        expect(bf.getValue(2, 2)).toBeCloseTo(6 - 18 * t, 13);
        expect(bf.getValue(2, 3)).toBeCloseTo(6 * t, 13);

        expect(bf.getValue(3, 0)).toBeCloseTo(-6, 13);
        expect(bf.getValue(3, 1)).toBeCloseTo(18, 13);
        expect(bf.getValue(3, 2)).toBeCloseTo(-18, 13);
        expect(bf.getValue(3, 3)).toBeCloseTo(6, 13);
    });

    it('at the domain endpoints of an open curve, one basis function is 1', () => {
        const bf = new BasisFunction(new BasisFunctionInput(6, 3));
        let r = bf.evaluate(0, 0);
        expect(bf.getValue(0, r.minIndex)).toBeCloseTo(1, 14);
        r = bf.evaluate(1, 0);
        expect(bf.getValue(0, r.maxIndex)).toBeCloseTo(1, 14);
    });
});

describe('BasisFunction partition of unity and nonnegativity', () => {
    it('open uniform, degrees 1 through 4', () => {
        for (let degree = 1; degree <= 4; ++degree) {
            const bf = new BasisFunction(new BasisFunctionInput(degree + 3, degree));
            for (let k = 0; k <= 100; ++k) {
                const t = k / 100;
                const { minIndex, maxIndex } = bf.evaluate(t, 0);
                expect(maxIndex - minIndex).toBe(degree);
                let sum = 0;
                for (let i = minIndex; i <= maxIndex; ++i) {
                    const value = bf.getValue(0, i);
                    expect(value).toBeGreaterThanOrEqual(0);
                    sum += value;
                }
                expect(sum).toBeCloseTo(1, 12);
            }
        }
    });

    it('nonuniform floating knots, degree 2', () => {
        const input = makeFloatingInput(5, 2, [0, 1, 2.5, 3, 4.5, 5, 6, 7]);
        const bf = new BasisFunction(input);
        const tMin = bf.getMinDomain();
        const tMax = bf.getMaxDomain();
        for (let k = 0; k <= 100; ++k) {
            const t = tMin + (k / 100) * (tMax - tMin);
            const { minIndex, maxIndex } = bf.evaluate(t, 0);
            let sum = 0;
            for (let i = minIndex; i <= maxIndex; ++i) {
                const value = bf.getValue(0, i);
                expect(value).toBeGreaterThanOrEqual(0);
                sum += value;
            }
            expect(sum).toBeCloseTo(1, 12);
        }
    });

    it('sum of derivatives of any order is 0 across the basis', () => {
        const input = makeFloatingInput(6, 3, [0, 1, 2, 3.5, 4, 5.5, 6.5, 7, 8, 9]);
        const bf = new BasisFunction(input);
        for (const t of [3.8, 4.6, 5.2, 5.9, 6.3]) {
            const { minIndex, maxIndex } = bf.evaluate(t, 3);
            for (let order = 1; order <= 3; ++order) {
                let sum = 0;
                for (let i = minIndex; i <= maxIndex; ++i) {
                    sum += bf.getValue(order, i);
                }
                expect(Math.abs(sum)).toBeLessThan(1e-10);
            }
        }
    });
});

describe('BasisFunction local support', () => {
    it('basis values vanish outside [minIndex, maxIndex]', () => {
        const bf = new BasisFunction(new BasisFunctionInput(6, 2));
        // Knots [0,0,0,0.25,0.5,0.75,1,1,1]; t = 0.6 lies in [0.5, 0.75).
        const { minIndex, maxIndex } = bf.evaluate(0.6, 0);
        expect(minIndex).toBe(2);
        expect(maxIndex).toBe(4);
        const numColumns = bf.getNumControls() + bf.getDegree();
        for (let i = 0; i < numColumns; ++i) {
            if (i < minIndex || i > maxIndex) {
                expect(bf.getValue(0, i)).toBe(0);
            }
        }
    });
});

describe('BasisFunction derivatives vs finite differences', () => {
    it('orders 1-3 match central differences of the next-lower order', () => {
        const input = makeFloatingInput(6, 3, [0, 1, 2, 3.5, 4, 5.5, 6.5, 7, 8, 9]);
        const bf = new BasisFunction(input);
        const h = 1e-5;
        const n = bf.getNumControls();
        for (const t of [3.8, 4.6, 5.9]) {
            for (let i = 0; i < n + bf.getDegree(); ++i) {
                for (let order = 1; order <= 3; ++order) {
                    const fp = basisValue(bf, order - 1, i, t + h);
                    const fm = basisValue(bf, order - 1, i, t - h);
                    const estimate = (fp - fm) / (2 * h);
                    const exact = basisValue(bf, order, i, t);
                    expect(Math.abs(exact - estimate)).toBeLessThan(1e-4);
                }
            }
        }
    });
});

describe('BasisFunction periodic vs open behavior', () => {
    it('periodic creation appends degree controls and is floating', () => {
        const bf = new BasisFunction(makePeriodicInput(4, 2));
        expect(bf.isPeriodic()).toBe(true);
        expect(bf.isOpen()).toBe(false);
        expect(bf.getNumControls()).toBe(6); // input controls + degree
        expect(bf.getNumKnots()).toBe(9);
        expect(bf.getMinDomain()).toBe(2); // knots[degree]
        expect(bf.getMaxDomain()).toBe(6); // knots[numControls]
    });

    it('periodic evaluation wraps t into the domain', () => {
        const bf = new BasisFunction(makePeriodicInput(4, 2));
        // Domain is [2,6] with length 4; 6.5 wraps to 2.5, -1.5 wraps to 2.5.
        const ref = bf.evaluate(2.5, 1);
        const refValues: number[] = [];
        const refDerivs: number[] = [];
        for (let i = ref.minIndex; i <= ref.maxIndex; ++i) {
            refValues.push(bf.getValue(0, i));
            refDerivs.push(bf.getValue(1, i));
        }

        for (const t of [6.5, -1.5, 10.5]) {
            const r = bf.evaluate(t, 1);
            expect(r.minIndex).toBe(ref.minIndex);
            expect(r.maxIndex).toBe(ref.maxIndex);
            for (let i = r.minIndex; i <= r.maxIndex; ++i) {
                expect(bf.getValue(0, i)).toBeCloseTo(refValues[i - r.minIndex], 12);
                expect(bf.getValue(1, i)).toBeCloseTo(refDerivs[i - r.minIndex], 12);
            }
        }
    });

    it('periodic basis satisfies partition of unity over its domain', () => {
        const bf = new BasisFunction(makePeriodicInput(5, 3));
        const tMin = bf.getMinDomain();
        const tMax = bf.getMaxDomain();
        for (let k = 0; k <= 100; ++k) {
            const t = tMin + (k / 100) * (tMax - tMin);
            const { minIndex, maxIndex } = bf.evaluate(t, 0);
            let sum = 0;
            for (let i = minIndex; i <= maxIndex; ++i) {
                sum += bf.getValue(0, i);
            }
            expect(sum).toBeCloseTo(1, 12);
        }
    });

    it('nonperiodic evaluation clamps t outside the domain', () => {
        const bf = new BasisFunction(new BasisFunctionInput(5, 2));
        const atMin = bf.evaluate(0, 0);
        const v0 = bf.getValue(0, atMin.minIndex);
        const below = bf.evaluate(-0.75, 0);
        expect(below.minIndex).toBe(atMin.minIndex);
        expect(below.maxIndex).toBe(atMin.maxIndex);
        expect(bf.getValue(0, below.minIndex)).toBe(v0);

        const atMax = bf.evaluate(1, 0);
        const v1 = bf.getValue(0, atMax.maxIndex);
        const above = bf.evaluate(2.25, 0);
        expect(above.minIndex).toBe(atMax.minIndex);
        expect(above.maxIndex).toBe(atMax.maxIndex);
        expect(bf.getValue(0, above.maxIndex)).toBe(v1);
    });
});

describe('BasisFunction argument validation', () => {
    it('evaluate rejects order > 3 and getValue rejects bad arguments', () => {
        const bf = new BasisFunction(new BasisFunctionInput(4, 2));
        expect(() => bf.evaluate(0.5, 4)).toThrow('Invalid order.');
        bf.evaluate(0.5, 0);
        expect(() => bf.getValue(4, 0)).toThrow('Invalid order.');
        expect(() => bf.getValue(0, -1)).toThrow('Invalid index.');
        expect(() => bf.getValue(0, bf.getNumControls() + bf.getDegree()))
            .toThrow('Invalid index.');
    });
});
