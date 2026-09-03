import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    getContainerAlignedBox,
    inContainerAlignedBox,
    mergeContainersAlignedBox
} from '../src/ContAlignedBox.js';
import { Vector } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

describe('getContainerAlignedBox', () => {
    it('computes the extremes of a 2D point set', () => {
        const points = [v(1, 5), v(-3, 2), v(4, -1), v(0, 0)];
        const result = getContainerAlignedBox(points);
        expect(result).not.toBeNull();
        expect(result!.min.values).toEqual([-3, -1]);
        expect(result!.max.values).toEqual([4, 5]);
    });

    it('computes the extremes of a 3D point set', () => {
        const points = [v(1, 2, 3), v(-1, 7, 0), v(5, -2, 2)];
        const result = getContainerAlignedBox(points)!;
        expect(result.min.values).toEqual([-1, -2, 0]);
        expect(result.max.values).toEqual([5, 7, 3]);
    });

    it('degenerates to a point box for a single point', () => {
        const result = getContainerAlignedBox([v(3, -4)])!;
        expect(result.min.values).toEqual([3, -4]);
        expect(result.max.values).toEqual([3, -4]);
    });

    it('returns null for an empty point set', () => {
        expect(getContainerAlignedBox([])).toBeNull();
    });

    it('contains every input point (randomized)', () => {
        let seed = 12345;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff - 0.5;
        };

        for (let trial = 0; trial < 50; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 20; ++i) {
                points.push(v(10 * rand(), 10 * rand(), 10 * rand()));
            }
            const result = getContainerAlignedBox(points)!;
            for (const point of points) {
                expect(inContainerAlignedBox(point, result)).toBe(true);
            }
            // The box is minimal: every face is touched by some point.
            for (let d = 0; d < 3; ++d) {
                const coords = points.map(p => p.values[d]);
                expect(result.min.values[d]).toBe(Math.min(...coords));
                expect(result.max.values[d]).toBe(Math.max(...coords));
            }
        }
    });
});

describe('inContainerAlignedBox', () => {
    const b = box([-1, -2], [3, 4]);

    it('accepts interior points', () => {
        expect(inContainerAlignedBox(v(0, 0), b)).toBe(true);
        expect(inContainerAlignedBox(v(2.9, 3.9), b)).toBe(true);
    });

    it('accepts boundary points (the boundary is part of the box)', () => {
        expect(inContainerAlignedBox(v(-1, -2), b)).toBe(true);
        expect(inContainerAlignedBox(v(3, 4), b)).toBe(true);
        expect(inContainerAlignedBox(v(-1, 0), b)).toBe(true);
        expect(inContainerAlignedBox(v(0, 4), b)).toBe(true);
    });

    it('rejects points outside along each axis', () => {
        expect(inContainerAlignedBox(v(-1.5, 0), b)).toBe(false);
        expect(inContainerAlignedBox(v(3.5, 0), b)).toBe(false);
        expect(inContainerAlignedBox(v(0, -2.5), b)).toBe(false);
        expect(inContainerAlignedBox(v(0, 4.5), b)).toBe(false);
    });

    it('handles a degenerate (zero-volume) box', () => {
        const flat = box([0, 0], [0, 1]);
        expect(inContainerAlignedBox(v(0, 0.5), flat)).toBe(true);
        expect(inContainerAlignedBox(v(1e-12, 0.5), flat)).toBe(false);
    });

    it('throws on a dimension mismatch', () => {
        expect(() => inContainerAlignedBox(v(0, 0, 0), b))
            .toThrow('inContainerAlignedBox: mismatched dimensions.');
    });
});

describe('mergeContainersAlignedBox', () => {
    it('computes the minimum box containing two disjoint boxes', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([2, -3], [4, 0.5]);
        const merge = mergeContainersAlignedBox(b0, b1);
        expect(merge.min.values).toEqual([0, -3]);
        expect(merge.max.values).toEqual([4, 1]);
    });

    it('is idempotent: merge(b, b) === b', () => {
        const b0 = box([-1, 2, 3], [5, 6, 7]);
        const merge = mergeContainersAlignedBox(b0, b0);
        expect(merge.equals(b0)).toBe(true);
    });

    it('is commutative', () => {
        const b0 = box([-1, 2], [5, 6]);
        const b1 = box([0, -4], [3, 8]);
        expect(mergeContainersAlignedBox(b0, b1)
            .equals(mergeContainersAlignedBox(b1, b0))).toBe(true);
    });

    it('returns the outer box when one box contains the other', () => {
        const outer = box([-10, -10], [10, 10]);
        const inner = box([-1, 2], [3, 4]);
        expect(mergeContainersAlignedBox(outer, inner).equals(outer)).toBe(true);
        expect(mergeContainersAlignedBox(inner, outer).equals(outer)).toBe(true);
    });

    it('is associative', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([-5, 2], [3, 3]);
        const b2 = box([2, -7], [9, 0]);
        const left = mergeContainersAlignedBox(
            mergeContainersAlignedBox(b0, b1), b2);
        const right = mergeContainersAlignedBox(
            b0, mergeContainersAlignedBox(b1, b2));
        expect(left.equals(right)).toBe(true);
    });

    it('contains all vertices of both input boxes', () => {
        const b0 = box([-1, -2, -3], [1, 0, 4]);
        const b1 = box([0, 5, -7], [2, 6, -1]);
        const merge = mergeContainersAlignedBox(b0, b1);
        for (const vertex of [...b0.getVertices(), ...b1.getVertices()]) {
            expect(inContainerAlignedBox(vertex, merge)).toBe(true);
        }
    });

    it('does not alias the inputs', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([2, 2], [3, 3]);
        const merge = mergeContainersAlignedBox(b0, b1);
        merge.min.values[0] = -100;
        expect(b0.min.values[0]).toBe(0);
    });

    it('throws on a dimension mismatch', () => {
        expect(() => mergeContainersAlignedBox(box([0, 0], [1, 1]), new AlignedBox(3)))
            .toThrow('mergeContainersAlignedBox: mismatched dimensions.');
    });
});
