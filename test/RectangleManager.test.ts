import { describe, it, expect } from 'vitest';
import { RectangleManager } from '../src/RectangleManager.js';
import { AlignedBox } from '../src/AlignedBox.js';
import { Vector } from '../src/Vector.js';

function rect(min: [number, number], max: [number, number]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

// A brute-force computation of the overlapping pairs, used as the independent
// cross-check for the sort-and-sweep manager.
function bruteForceOverlap(rectangles: AlignedBox[]): string[] {
    const pairs: string[] = [];
    for (let i = 0; i < rectangles.length; ++i) {
        for (let j = i + 1; j < rectangles.length; ++j) {
            const r0 = rectangles[i];
            const r1 = rectangles[j];
            let overlap = true;
            for (let d = 0; d < 2; ++d) {
                if (r0.max.values[d] < r1.min.values[d]
                    || r0.min.values[d] > r1.max.values[d]) {
                    overlap = false;
                    break;
                }
            }
            if (overlap) {
                pairs.push(`${i},${j}`);
            }
        }
    }
    return pairs;
}

function managerOverlap(manager: RectangleManager): string[] {
    return manager.getOverlap().map(k => `${k.V[0]},${k.V[1]}`);
}

// A small deterministic pseudorandom generator (mulberry32) so the randomized
// cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function randomRectangle(rand: () => number): AlignedBox {
    const min: [number, number] = [0, 0];
    const max: [number, number] = [0, 0];
    for (let d = 0; d < 2; ++d) {
        // Integer-valued coordinates keep all comparisons exact.
        const lo = Math.floor(rand() * 20);
        const len = 1 + Math.floor(rand() * 5);
        min[d] = lo;
        max[d] = lo + len;
    }
    return rect(min, max);
}

describe('RectangleManager', () => {
    it('finds the overlaps of a simple configuration', () => {
        // Rectangles 0 and 1 overlap; rectangle 2 is far away.
        const rectangles = [
            rect([0, 0], [2, 2]),
            rect([1, 1], [3, 3]),
            rect([10, 10], [11, 11])
        ];
        const manager = new RectangleManager(rectangles);
        expect(managerOverlap(manager)).toEqual(['0,1']);
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(rectangles));
    });

    it('treats touching rectangles as overlapping', () => {
        const rectangles = [
            rect([0, 0], [1, 1]),
            rect([1, 0], [2, 1])
        ];
        const manager = new RectangleManager(rectangles);
        expect(managerOverlap(manager)).toEqual(['0,1']);
        expect(bruteForceOverlap(rectangles)).toEqual(['0,1']);
    });

    it('does not report x-only or y-only overlap', () => {
        // The x-intervals overlap but the y-intervals do not.
        const rectangles = [
            rect([0, 0], [4, 1]),
            rect([2, 5], [6, 7])
        ];
        const manager = new RectangleManager(rectangles);
        expect(managerOverlap(manager)).toEqual([]);

        // The y-intervals overlap but the x-intervals do not.
        const rectangles2 = [
            rect([0, 0], [1, 4]),
            rect([5, 2], [7, 6])
        ];
        expect(managerOverlap(new RectangleManager(rectangles2))).toEqual([]);
    });

    it('reports no overlaps for well-separated rectangles', () => {
        const rectangles = [
            rect([0, 0], [1, 1]),
            rect([2, 2], [3, 3]),
            rect([4, 4], [5, 5])
        ];
        const manager = new RectangleManager(rectangles);
        expect(managerOverlap(manager)).toEqual([]);
    });

    it('reports every pair for a fully overlapping configuration', () => {
        const rectangles = [
            rect([0, 0], [10, 10]),
            rect([1, 1], [9, 9]),
            rect([2, 2], [8, 8]),
            rect([3, 3], [7, 7])
        ];
        const manager = new RectangleManager(rectangles);
        expect(managerOverlap(manager))
            .toEqual(['0,1', '0,2', '0,3', '1,2', '1,3', '2,3']);
    });

    it('handles an empty and a single-rectangle collection', () => {
        expect(managerOverlap(new RectangleManager([]))).toEqual([]);
        const one = [rect([0, 0], [1, 1])];
        expect(managerOverlap(new RectangleManager(one))).toEqual([]);
    });

    it('copies rectangles in and out of the manager', () => {
        const rectangles = [rect([0, 0], [1, 1]), rect([5, 5], [6, 6])];
        const manager = new RectangleManager(rectangles);

        const moved = rect([0, 0], [6, 6]);
        manager.setRectangle(0, moved);
        // setRectangle writes through to the caller's array (upstream holds a
        // reference to it) but stores a copy of the input rectangle.
        expect(rectangles[0].min.values).toEqual([0, 0]);
        expect(rectangles[0].max.values).toEqual([6, 6]);
        moved.max.values[0] = 100;
        expect(rectangles[0].max.values[0]).toBe(6);

        const fetched = manager.getRectangle(0);
        fetched.max.values[0] = 55;
        expect(rectangles[0].max.values[0]).toBe(6);
    });

    it('updates the overlap set incrementally after moves', () => {
        const rectangles = [
            rect([0, 0], [2, 2]),
            rect([10, 0], [12, 2]),
            rect([20, 0], [22, 2])
        ];
        const manager = new RectangleManager(rectangles);
        expect(managerOverlap(manager)).toEqual([]);

        // Move rectangle 1 onto rectangle 0: the pair becomes overlapping.
        manager.setRectangle(1, rect([1, 0], [3, 2]));
        manager.update();
        expect(managerOverlap(manager)).toEqual(['0,1']);
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(rectangles));

        // Move it back off: the pair is removed again.
        manager.setRectangle(1, rect([10, 0], [12, 2]));
        manager.update();
        expect(managerOverlap(manager)).toEqual([]);
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(rectangles));

        // Grow rectangle 0 to cover everything.
        manager.setRectangle(0, rect([-1, -1], [30, 30]));
        manager.update();
        expect(managerOverlap(manager)).toEqual(['0,1', '0,2']);
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(rectangles));
    });

    it('separates rectangles that only lose their y-overlap', () => {
        // The x-intervals keep overlapping; only the y-move removes the pair.
        const rectangles = [
            rect([0, 0], [4, 4]),
            rect([2, 2], [6, 6])
        ];
        const manager = new RectangleManager(rectangles);
        expect(managerOverlap(manager)).toEqual(['0,1']);

        manager.setRectangle(1, rect([2, 10], [6, 14]));
        manager.update();
        expect(managerOverlap(manager)).toEqual([]);
    });

    it('matches brute force after initialize for random rectangles', () => {
        const rand = makeRandom(0x2EC7);
        for (let trial = 0; trial < 60; ++trial) {
            const rectangles: AlignedBox[] = [];
            const count = 2 + Math.floor(rand() * 12);
            for (let i = 0; i < count; ++i) {
                rectangles.push(randomRectangle(rand));
            }
            const manager = new RectangleManager(rectangles);
            expect(managerOverlap(manager))
                .toEqual(bruteForceOverlap(rectangles));
        }
    });

    it('matches brute force after many incremental moves', () => {
        const rand = makeRandom(0x5A17);
        for (let trial = 0; trial < 25; ++trial) {
            const rectangles: AlignedBox[] = [];
            const count = 3 + Math.floor(rand() * 8);
            for (let i = 0; i < count; ++i) {
                rectangles.push(randomRectangle(rand));
            }
            const manager = new RectangleManager(rectangles);
            expect(managerOverlap(manager))
                .toEqual(bruteForceOverlap(rectangles));

            for (let move = 0; move < 12; ++move) {
                const i = Math.floor(rand() * count);
                const current = manager.getRectangle(i);
                // Small translations, the intended use of the incremental
                // update (the endpoints stay nearly sorted).
                const dx = Math.floor(rand() * 5) - 2;
                const dy = Math.floor(rand() * 5) - 2;
                const moved = rect(
                    [current.min.values[0] + dx, current.min.values[1] + dy],
                    [current.max.values[0] + dx, current.max.values[1] + dy]);
                manager.setRectangle(i, moved);
                manager.update();
                expect(managerOverlap(manager))
                    .toEqual(bruteForceOverlap(rectangles));
            }
        }
    });

    it('matches brute force after large jumps and re-initialization', () => {
        const rand = makeRandom(0x9C0B);
        const rectangles: AlignedBox[] = [];
        for (let i = 0; i < 10; ++i) {
            rectangles.push(randomRectangle(rand));
        }
        const manager = new RectangleManager(rectangles);
        for (let move = 0; move < 20; ++move) {
            const i = Math.floor(rand() * rectangles.length);
            manager.setRectangle(i, randomRectangle(rand));
            manager.update();
            expect(managerOverlap(manager))
                .toEqual(bruteForceOverlap(rectangles));
        }

        // A full re-initialization reproduces the same overlap set.
        manager.initialize();
        expect(managerOverlap(manager))
            .toEqual(bruteForceOverlap(rectangles));
    });
});
