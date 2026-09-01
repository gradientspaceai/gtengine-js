import { describe, it, expect } from 'vitest';
import { BoxManager } from '../src/BoxManager';
import { AlignedBox } from '../src/AlignedBox';
import { Vector } from '../src/Vector';

function box(min: [number, number, number], max: [number, number, number]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

// A brute-force computation of the overlapping pairs, used as the independent
// cross-check for the sort-and-sweep manager.
function bruteForceOverlap(boxes: AlignedBox[]): string[] {
    const pairs: string[] = [];
    for (let i = 0; i < boxes.length; ++i) {
        for (let j = i + 1; j < boxes.length; ++j) {
            const b0 = boxes[i];
            const b1 = boxes[j];
            let overlap = true;
            for (let d = 0; d < 3; ++d) {
                if (b0.max.values[d] < b1.min.values[d]
                    || b0.min.values[d] > b1.max.values[d]) {
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

function managerOverlap(manager: BoxManager): string[] {
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

function randomBox(rand: () => number): AlignedBox {
    const min: [number, number, number] = [0, 0, 0];
    const max: [number, number, number] = [0, 0, 0];
    for (let d = 0; d < 3; ++d) {
        // Integer-valued coordinates keep all comparisons exact.
        const lo = Math.floor(rand() * 20);
        const len = 1 + Math.floor(rand() * 5);
        min[d] = lo;
        max[d] = lo + len;
    }
    return box(min, max);
}

describe('BoxManager', () => {
    it('finds the overlaps of a simple configuration', () => {
        // Boxes 0 and 1 overlap; box 2 is far away.
        const boxes = [
            box([0, 0, 0], [2, 2, 2]),
            box([1, 1, 1], [3, 3, 3]),
            box([10, 10, 10], [11, 11, 11])
        ];
        const manager = new BoxManager(boxes);
        expect(managerOverlap(manager)).toEqual(['0,1']);
    });

    it('treats touching boxes as overlapping', () => {
        const boxes = [
            box([0, 0, 0], [1, 1, 1]),
            box([1, 0, 0], [2, 1, 1])
        ];
        const manager = new BoxManager(boxes);
        expect(managerOverlap(manager)).toEqual(['0,1']);
        expect(bruteForceOverlap(boxes)).toEqual(['0,1']);
    });

    it('reports no overlaps for well-separated boxes', () => {
        const boxes = [
            box([0, 0, 0], [1, 1, 1]),
            box([2, 2, 2], [3, 3, 3]),
            box([4, 4, 4], [5, 5, 5])
        ];
        const manager = new BoxManager(boxes);
        expect(managerOverlap(manager)).toEqual([]);
    });

    it('separates boxes that overlap in x and y but not z', () => {
        const boxes = [
            box([0, 0, 0], [4, 4, 1]),
            box([1, 1, 5], [3, 3, 6])
        ];
        const manager = new BoxManager(boxes);
        expect(managerOverlap(manager)).toEqual([]);
    });

    it('returns overlap keys sorted lexicographically with V[0] < V[1]', () => {
        const boxes = [
            box([0, 0, 0], [10, 10, 10]),
            box([1, 1, 1], [2, 2, 2]),
            box([3, 3, 3], [4, 4, 4]),
            box([1, 1, 1], [4, 4, 4])
        ];
        const manager = new BoxManager(boxes);
        const keys = manager.getOverlap();
        for (const k of keys) {
            expect(k.V[0]).toBeLessThan(k.V[1]);
        }
        const flat = managerOverlap(manager);
        expect(flat).toEqual([...flat].sort());
        expect(flat).toEqual(bruteForceOverlap(boxes));
    });

    it('matches brute force after initialize on random configurations', () => {
        const rand = makeRandom(12345);
        for (let trial = 0; trial < 20; ++trial) {
            const boxes: AlignedBox[] = [];
            for (let i = 0; i < 12; ++i) {
                boxes.push(randomBox(rand));
            }
            const manager = new BoxManager(boxes);
            expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));
        }
    });

    it('matches brute force after repeated setBox/update moves', () => {
        const rand = makeRandom(98765);
        const boxes: AlignedBox[] = [];
        for (let i = 0; i < 10; ++i) {
            boxes.push(randomBox(rand));
        }
        const manager = new BoxManager(boxes);
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));

        for (let move = 0; move < 40; ++move) {
            const i = Math.floor(rand() * boxes.length);
            // Move the box by a small integer offset in each dimension.
            const b = manager.getBox(i);
            const min: [number, number, number] = [0, 0, 0];
            const max: [number, number, number] = [0, 0, 0];
            for (let d = 0; d < 3; ++d) {
                const delta = Math.floor(rand() * 5) - 2;
                min[d] = b.min.values[d] + delta;
                max[d] = b.max.values[d] + delta;
            }
            manager.setBox(i, box(min, max));
            manager.update();
            expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));
        }
    });

    it('matches brute force when boxes are resized, not only translated', () => {
        const rand = makeRandom(2024);
        const boxes: AlignedBox[] = [];
        for (let i = 0; i < 8; ++i) {
            boxes.push(randomBox(rand));
        }
        const manager = new BoxManager(boxes);
        for (let move = 0; move < 30; ++move) {
            const i = Math.floor(rand() * boxes.length);
            manager.setBox(i, randomBox(rand));
            manager.update();
            expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));
        }
    });

    it('reinitializes after boxes are inserted or removed from the array', () => {
        const rand = makeRandom(555);
        const boxes: AlignedBox[] = [];
        for (let i = 0; i < 6; ++i) {
            boxes.push(randomBox(rand));
        }
        const manager = new BoxManager(boxes);
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));

        // Insert two boxes and reinitialize, as documented by initialize().
        boxes.push(randomBox(rand));
        boxes.push(box([0, 0, 0], [20, 20, 20]));
        manager.initialize();
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));

        // Remove boxes and reinitialize.
        boxes.splice(1, 3);
        manager.initialize();
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));

        // Incremental updates still work after the reinitialization.
        for (let move = 0; move < 10; ++move) {
            const i = Math.floor(rand() * boxes.length);
            manager.setBox(i, randomBox(rand));
            manager.update();
            expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));
        }
    });

    it('writes setBox through to the caller array and copies the box', () => {
        const boxes = [
            box([0, 0, 0], [1, 1, 1]),
            box([5, 5, 5], [6, 6, 6])
        ];
        const manager = new BoxManager(boxes);
        const moved = box([0, 0, 0], [10, 10, 10]);
        manager.setBox(1, moved);
        expect(boxes[1].max.values[0]).toBe(10);

        // The stored box is a copy, so mutating the argument afterward has no
        // effect on the manager.
        moved.max.values[0] = 1000;
        expect(manager.getBox(1).max.values[0]).toBe(10);
        expect(boxes[1].max.values[0]).toBe(10);
    });

    it('returns a copy from getBox', () => {
        const boxes = [box([0, 0, 0], [1, 1, 1]), box([2, 2, 2], [3, 3, 3])];
        const manager = new BoxManager(boxes);
        const b = manager.getBox(0);
        b.min.values[0] = -100;
        expect(boxes[0].min.values[0]).toBe(0);
    });

    it('handles a single box and an empty array', () => {
        const one = [box([0, 0, 0], [1, 1, 1])];
        expect(managerOverlap(new BoxManager(one))).toEqual([]);

        const none: AlignedBox[] = [];
        expect(managerOverlap(new BoxManager(none))).toEqual([]);
    });

    it('handles degenerate boxes with min equal to max', () => {
        const boxes = [
            box([1, 1, 1], [1, 1, 1]),
            box([0, 0, 0], [2, 2, 2]),
            box([5, 5, 5], [5, 5, 5])
        ];
        const manager = new BoxManager(boxes);
        expect(managerOverlap(manager)).toEqual(bruteForceOverlap(boxes));
        expect(managerOverlap(manager)).toEqual(['0,1']);
    });

    it('handles many identical boxes (all pairs overlap)', () => {
        const boxes = [
            box([0, 0, 0], [1, 1, 1]),
            box([0, 0, 0], [1, 1, 1]),
            box([0, 0, 0], [1, 1, 1])
        ];
        const manager = new BoxManager(boxes);
        expect(managerOverlap(manager)).toEqual(['0,1', '0,2', '1,2']);
    });
});
