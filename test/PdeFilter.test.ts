import { describe, it, expect } from 'vitest';
import { PdeFilter, PdeFilterScaleType } from '../src/PdeFilter';

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
