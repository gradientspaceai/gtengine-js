import { describe, it, expect } from 'vitest';
import { PdeFilter, PdeFilterScaleType } from '../src/PdeFilter.js';
import { check, expectClose, fc, scaled } from './helpers/arbitraries.js';

// Concrete subclass exposing the protected scaling state and recording the
// template-method call order.
class TestFilter extends PdeFilter {
    calls: string[] = [];

    constructor(quantity: number, data: ArrayLike<number>,
        borderValue: number, scaleType: PdeFilterScaleType) {
        super(quantity, data, borderValue, scaleType);
    }

    protected override onPreUpdate(): void {
        this.calls.push('pre');
    }

    protected override onUpdate(): void {
        this.calls.push('update');
    }

    protected override onPostUpdate(): void {
        this.calls.push('post');
    }

    getMin(): number { return this.mMin; }
    getOffset(): number { return this.mOffset; }
    getScale(): number { return this.mScale; }
}

describe('PdeFilter', () => {
    it('provides member access', () => {
        const data = [1, 2, 3, 4];
        const filter = new TestFilter(4, data, 0.5, PdeFilterScaleType.NONE);
        expect(filter.getQuantity()).toBe(4);
        expect(filter.getBorderValue()).toBe(0.5);
        expect(filter.getScaleType()).toBe(PdeFilterScaleType.NONE);
        expect(filter.getTimeStep()).toBe(0);
        filter.setTimeStep(0.125);
        expect(filter.getTimeStep()).toBe(0.125);
    });

    it('update calls onPreUpdate, onUpdate and onPostUpdate in order', () => {
        const filter = new TestFilter(2, [0, 1], 0, PdeFilterScaleType.NONE);
        filter.update();
        expect(filter.calls).toEqual(['pre', 'update', 'post']);
        filter.update();
        expect(filter.calls).toEqual(['pre', 'update', 'post', 'pre', 'update', 'post']);
    });

    it('computes NONE scaling', () => {
        const filter = new TestFilter(3, [2, 6, 4], 0, PdeFilterScaleType.NONE);
        expect(filter.getMin()).toBe(2);
        expect(filter.getOffset()).toBe(0);
        expect(filter.getScale()).toBe(1);
    });

    it('computes UNIT scaling: d\' = (d-min)/(max-min)', () => {
        const filter = new TestFilter(3, [2, 6, 4], 0, PdeFilterScaleType.UNIT);
        expect(filter.getMin()).toBe(2);
        expect(filter.getOffset()).toBe(0);
        expect(filter.getScale()).toBeCloseTo(1 / 4, 15);
    });

    it('computes SYMMETRIC scaling: d\' = -1 + 2(d-min)/(max-min)', () => {
        const filter = new TestFilter(3, [2, 6, 4], 0, PdeFilterScaleType.SYMMETRIC);
        expect(filter.getMin()).toBe(2);
        expect(filter.getOffset()).toBe(-1);
        expect(filter.getScale()).toBeCloseTo(2 / 4, 15);
    });

    it('computes PRESERVE_ZERO scaling with max >= -min', () => {
        const filter = new TestFilter(3, [-2, 6, 4], 0, PdeFilterScaleType.PRESERVE_ZERO);
        expect(filter.getOffset()).toBe(0);
        expect(filter.getScale()).toBeCloseTo(1 / 6, 15);
        // mMin is reset to zero so unscaling preserves zero.
        expect(filter.getMin()).toBe(0);
    });

    it('computes PRESERVE_ZERO scaling with max < -min', () => {
        const filter = new TestFilter(3, [-8, 6, 4], 0, PdeFilterScaleType.PRESERVE_ZERO);
        expect(filter.getOffset()).toBe(0);
        expect(filter.getScale()).toBeCloseTo(-1 / -8, 15);
        expect(filter.getMin()).toBe(0);
    });

    it('uses identity scaling for constant data', () => {
        const filter = new TestFilter(3, [5, 5, 5], 0, PdeFilterScaleType.UNIT);
        expect(filter.getMin()).toBe(5);
        expect(filter.getOffset()).toBe(0);
        expect(filter.getScale()).toBe(1);
    });

    it('treats Number.MAX_VALUE border value as the Neumann marker', () => {
        const filter = new TestFilter(2, [0, 1], Number.MAX_VALUE,
            PdeFilterScaleType.NONE);
        expect(filter.getBorderValue()).toBe(Number.MAX_VALUE);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream PdeFilter.h constructor.
// ---------------------------------------------------------------------------

describe('PdeFilter verification', () => {
    // scaled() draws from a uniform grid, so consecutive distinct samples
    // are never a subnormal apart; finite() would make (max - min) subnormal
    // and 1 / (max - min) infinite, which says nothing about the port.
    const data = fc.array(scaled(-40, 40), { minLength: 1, maxLength: 40 });
    const scaleType = fc.constantFrom(
        PdeFilterScaleType.NONE, PdeFilterScaleType.UNIT,
        PdeFilterScaleType.SYMMETRIC, PdeFilterScaleType.PRESERVE_ZERO);

    // The transform the 2D and 3D filters apply to the image on construction:
    // mOffset + (value - mMin) * mScale (see PdeFilter2.h line 59).
    const transform = (filter: TestFilter, value: number) =>
        filter.getOffset() + (value - filter.getMin()) * filter.getScale();

    it('records the data extremes the way the upstream scan does', () => {
        check(fc.tuple(data, scaleType), ([values, type]) => {
            const filter = new TestFilter(values.length, values, 0, type);
            const minValue = Math.min(...values);
            const maxValue = Math.max(...values);
            if (type === PdeFilterScaleType.PRESERVE_ZERO && minValue !== maxValue) {
                // PRESERVE_ZERO resets mMin to zero after computing mScale.
                expect(filter.getMin()).toBe(0);
            } else {
                expect(filter.getMin()).toBe(minValue);
            }
            expect(filter.getQuantity()).toBe(values.length);
        });
    });

    it('UNIT maps the data range onto [0, 1]', () => {
        check(data.filter(v => Math.min(...v) < Math.max(...v)), values => {
            const filter = new TestFilter(values.length, values, 0,
                PdeFilterScaleType.UNIT);
            const minValue = Math.min(...values);
            const maxValue = Math.max(...values);
            expectClose(transform(filter, minValue), 0);
            expectClose(transform(filter, maxValue), 1, 1e-12, 1e-12);
            for (const value of values) {
                const scaled = transform(filter, value);
                expect(scaled).toBeGreaterThanOrEqual(-1e-12);
                expect(scaled).toBeLessThanOrEqual(1 + 1e-12);
            }
        });
    });

    it('SYMMETRIC maps the data range onto [-1, 1]', () => {
        check(data.filter(v => Math.min(...v) < Math.max(...v)), values => {
            const filter = new TestFilter(values.length, values, 0,
                PdeFilterScaleType.SYMMETRIC);
            const minValue = Math.min(...values);
            const maxValue = Math.max(...values);
            expectClose(transform(filter, minValue), -1, 1e-12, 1e-12);
            expectClose(transform(filter, maxValue), 1, 1e-12, 1e-12);
            // SYMMETRIC is the UNIT map composed with d -> -1 + 2 d.
            const unit = new TestFilter(values.length, values, 0,
                PdeFilterScaleType.UNIT);
            for (const value of values) {
                expectClose(transform(filter, value),
                    -1 + 2 * transform(unit, value), 1e-12, 1e-12);
            }
        });
    });

    it('PRESERVE_ZERO fixes zero and sends the larger magnitude to +-1', () => {
        check(data.filter(v => Math.min(...v) < Math.max(...v)), values => {
            const filter = new TestFilter(values.length, values, 0,
                PdeFilterScaleType.PRESERVE_ZERO);
            const minValue = Math.min(...values);
            const maxValue = Math.max(...values);
            // The name is accurate: the transform of 0 is 0.
            expect(transform(filter, 0)).toBe(0);
            if (maxValue >= -minValue) {
                expect(filter.getScale()).toBe(1 / maxValue);
                expectClose(transform(filter, maxValue), 1, 1e-12, 1e-12);
            } else {
                expect(filter.getScale()).toBe(-1 / minValue);
                expectClose(transform(filter, minValue), -1, 1e-12, 1e-12);
            }
            for (const value of values) {
                expect(Math.abs(transform(filter, value))).toBeLessThanOrEqual(1 + 1e-12);
            }
        });
    });

    it('NONE shifts the data by the data minimum (upstream #60)', () => {
        // Upstream documents ScaleType::NONE as "the data is processed as
        // is", but mMin keeps the data minimum while mOffset is 0 and mScale
        // is 1, so the transform applied by the derived filters is
        // d - min, not d. Preserved faithfully.
        check(data.filter(v => Math.min(...v) < Math.max(...v)), values => {
            const filter = new TestFilter(values.length, values, 0,
                PdeFilterScaleType.NONE);
            expect(filter.getOffset()).toBe(0);
            expect(filter.getScale()).toBe(1);
            const minValue = Math.min(...values);
            expect(filter.getMin()).toBe(minValue);
            for (const value of values) {
                expect(transform(filter, value)).toBe(value - minValue);
            }
        });
    });

    it('constant data uses the identity scale for every scale type', () => {
        check(fc.tuple(scaled(-40, 40), fc.integer({ min: 1, max: 20 }), scaleType),
            ([value, quantity, type]) => {
                const filter = new TestFilter(quantity,
                    new Array<number>(quantity).fill(value), 0, type);
                expect(filter.getOffset()).toBe(0);
                expect(filter.getScale()).toBe(1);
                expect(filter.getMin()).toBe(value);
                // The stored image becomes all zeros, whatever the constant.
                expect(transform(filter, value)).toBe(0);
            });
    });

    it('only the first quantity samples are scanned', () => {
        // Upstream reads data[0 .. mQuantity-1]; a shorter quantity must
        // ignore the tail of the array.
        check(fc.tuple(data, fc.nat({ max: 40 })), ([values, cut]) => {
            const quantity = 1 + (cut % values.length);
            const filter = new TestFilter(quantity, values, 0, PdeFilterScaleType.UNIT);
            const prefix = values.slice(0, quantity);
            expect(filter.getMin()).toBe(Math.min(...prefix));
            const maxValue = Math.max(...prefix);
            if (Math.min(...prefix) < maxValue) {
                expect(filter.getScale()).toBe(1 / (maxValue - Math.min(...prefix)));
            }
        });
    });

    it('the border value is stored verbatim, MAX_VALUE marking Neumann', () => {
        check(fc.tuple(data, scaled(-1e6, 1e6)), ([values, border]) => {
            const dirichlet = new TestFilter(values.length, values, border,
                PdeFilterScaleType.UNIT);
            expect(dirichlet.getBorderValue()).toBe(border);
            expect(dirichlet.getBorderValue() === Number.MAX_VALUE).toBe(false);
            const neumann = new TestFilter(values.length, values, Number.MAX_VALUE,
                PdeFilterScaleType.UNIT);
            expect(neumann.getBorderValue()).toBe(Number.MAX_VALUE);
        });
    });

    it('update runs the template method in order on every call', () => {
        check(fc.tuple(data, fc.integer({ min: 1, max: 6 })), ([values, iterations]) => {
            const filter = new TestFilter(values.length, values, 0,
                PdeFilterScaleType.UNIT);
            for (let i = 0; i < iterations; ++i) {
                filter.update();
            }
            const expected: string[] = [];
            for (let i = 0; i < iterations; ++i) {
                expected.push('pre', 'update', 'post');
            }
            expect(filter.calls).toEqual(expected);
        });
    });
});
