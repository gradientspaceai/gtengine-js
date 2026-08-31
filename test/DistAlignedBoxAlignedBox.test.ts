import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { DistAlignedBoxAlignedBox } from '../src/DistAlignedBoxAlignedBox';
import { Vector } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

// The reference distance: for aligned boxes the squared distance is the sum
// of the per-axis interval gaps.
function referenceSqrDistance(b0: AlignedBox, b1: AlignedBox): number {
    let sum = 0;
    for (let i = 0; i < b0.min.size; ++i) {
        const gap = Math.max(0,
            Math.max(b0.min.values[i] - b1.max.values[i],
                b1.min.values[i] - b0.max.values[i]));
        sum += gap * gap;
    }
    return sum;
}

describe('DistAlignedBoxAlignedBox', () => {
    const query = new DistAlignedBoxAlignedBox();

    it('reports zero distance for overlapping boxes', () => {
        const b0 = box([0, 0], [2, 2]);
        const b1 = box([1, 1], [3, 3]);
        const result = query.compute(b0, b1);
        expect(result.distance).toBe(0);
        // The closest sets are the overlap box for both inputs.
        expect(result.closest[0].min.values).toEqual([1, 1]);
        expect(result.closest[0].max.values).toEqual([2, 2]);
        expect(result.closest[1].min.values).toEqual([1, 1]);
        expect(result.closest[1].max.values).toEqual([2, 2]);
    });

    it('reports zero distance for a contained box', () => {
        const b0 = box([0, 0, 0], [10, 10, 10]);
        const b1 = box([2, 3, 4], [5, 6, 7]);
        const result = query.compute(b0, b1);
        expect(result.distance).toBe(0);
        expect(result.closest[0].min.values).toEqual([2, 3, 4]);
        expect(result.closest[0].max.values).toEqual([5, 6, 7]);
    });

    it('measures a one-axis gap', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([4, 0], [5, 1]);
        const result = query.compute(b0, b1);
        expect(result.distance).toBeCloseTo(3, 12);
        // The closest sets degenerate in x and span the overlap in y.
        expect(result.closest[0].min.values).toEqual([1, 0]);
        expect(result.closest[0].max.values).toEqual([1, 1]);
        expect(result.closest[1].min.values).toEqual([4, 0]);
        expect(result.closest[1].max.values).toEqual([4, 1]);
    });

    it('measures a corner-to-corner gap', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([4, 5], [6, 7]);
        const result = query.compute(b0, b1);
        expect(result.sqrDistance).toBeCloseTo(9 + 16, 12);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[0].min.values).toEqual([1, 1]);
        expect(result.closest[0].max.values).toEqual([1, 1]);
        expect(result.closest[1].min.values).toEqual([4, 5]);
        expect(result.closest[1].max.values).toEqual([4, 5]);
    });

    it('reports zero distance for boxes that just touch', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([1, 0], [2, 1]);
        const result = query.compute(b0, b1);
        expect(result.distance).toBe(0);
        expect(result.closest[0].min.values[0]).toBe(1);
        expect(result.closest[0].max.values[0]).toBe(1);
    });

    it('is symmetric', () => {
        const b0 = box([-1, 2, 0], [3, 4, 1]);
        const b1 = box([7, -5, 0.5], [9, -1, 2]);
        const a = query.compute(b0, b1);
        const b = query.compute(b1, b0);
        expect(b.distance).toBeCloseTo(a.distance, 12);
        expect(b.closest[1].min.values).toEqual(a.closest[0].min.values);
        expect(b.closest[0].min.values).toEqual(a.closest[1].min.values);
    });

    it('agrees with the per-axis gap reference and the closest sets realize '
        + 'the distance', () => {
        let seed = 2468;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 10 - 5;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const lo0 = [rand(), rand(), rand()];
            const lo1 = [rand(), rand(), rand()];
            const b0 = box(lo0, lo0.map(x => x + Math.abs(rand()) + 0.1));
            const b1 = box(lo1, lo1.map(x => x + Math.abs(rand()) + 0.1));
            const result = query.compute(b0, b1);
            expect(result.sqrDistance).toBeCloseTo(
                referenceSqrDistance(b0, b1), 9);

            // The centers of the closest sets realize the distance and lie
            // inside their respective input boxes.
            let sum = 0;
            for (let i = 0; i < 3; ++i) {
                const p0 = 0.5 * (result.closest[0].min.values[i]
                    + result.closest[0].max.values[i]);
                const p1 = 0.5 * (result.closest[1].min.values[i]
                    + result.closest[1].max.values[i]);
                expect(p0).toBeGreaterThanOrEqual(b0.min.values[i] - 1e-12);
                expect(p0).toBeLessThanOrEqual(b0.max.values[i] + 1e-12);
                expect(p1).toBeGreaterThanOrEqual(b1.min.values[i] - 1e-12);
                expect(p1).toBeLessThanOrEqual(b1.max.values[i] + 1e-12);
                sum += (p0 - p1) * (p0 - p1);
            }
            expect(sum).toBeCloseTo(result.sqrDistance, 9);
        }
    });
});
