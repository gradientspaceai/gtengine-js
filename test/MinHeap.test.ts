import { describe, it, expect } from 'vitest';
import { MinHeap, type MinHeapRecord } from '../src/MinHeap.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

describe('MinHeap', () => {
    it('drains inserted elements in ascending value order', () => {
        const heap = new MinHeap<string>(8);
        const entries: Array<[string, number]> = [
            ['e', 5], ['b', 2], ['g', 7], ['a', 1], ['c', 3], ['h', 8], ['f', 6], ['d', 4]
        ];
        for (const [key, value] of entries) {
            expect(heap.insert(key, value)).not.toBeNull();
            expect(heap.isValid()).toBe(true);
        }
        expect(heap.getNumElements()).toBe(8);

        const drained: Array<[string, number]> = [];
        for (let r = heap.remove(); r !== null; r = heap.remove()) {
            drained.push([r.key, r.value]);
            expect(heap.isValid()).toBe(true);
        }
        expect(drained).toEqual([
            ['a', 1], ['b', 2], ['c', 3], ['d', 4], ['e', 5], ['f', 6], ['g', 7], ['h', 8]
        ]);
        expect(heap.getNumElements()).toBe(0);
    });

    it('getMinimum reads the root without removing it', () => {
        const heap = new MinHeap<number>(4);
        expect(heap.getMinimum()).toBeNull();
        heap.insert(10, 3);
        heap.insert(20, 1);
        heap.insert(30, 2);
        expect(heap.getMinimum()).toEqual({ key: 20, value: 1 });
        expect(heap.getNumElements()).toBe(3);
    });

    it('insert returns null when the heap is full', () => {
        const heap = new MinHeap<number>(2);
        expect(heap.insert(1, 1)).not.toBeNull();
        expect(heap.insert(2, 2)).not.toBeNull();
        expect(heap.insert(3, 3)).toBeNull();
        expect(heap.getNumElements()).toBe(2);
    });

    it('remove returns null when the heap is empty', () => {
        const heap = new MinHeap<number>(2);
        expect(heap.remove()).toBeNull();
    });

    it('supports duplicate values', () => {
        const heap = new MinHeap<number>(6);
        for (const v of [2, 1, 2, 1, 2, 1]) {
            heap.insert(v, v);
        }
        const values: number[] = [];
        for (let r = heap.remove(); r !== null; r = heap.remove()) {
            values.push(r.value);
        }
        expect(values).toEqual([1, 1, 1, 2, 2, 2]);
    });

    it('update decreases a value and restores the heap', () => {
        const heap = new MinHeap<string>(4);
        heap.insert('a', 10);
        const rec = heap.insert('b', 20)!;
        heap.insert('c', 30);
        heap.update(rec, 1);
        expect(heap.isValid()).toBe(true);
        expect(heap.getMinimum()).toEqual({ key: 'b', value: 1 });
    });

    it('update increases a value and restores the heap', () => {
        const heap = new MinHeap<string>(4);
        const rec = heap.insert('a', 1)!;
        heap.insert('b', 20);
        heap.insert('c', 30);
        heap.update(rec, 100);
        expect(heap.isValid()).toBe(true);
        expect(heap.getMinimum()).toEqual({ key: 'b', value: 20 });
        const drained: number[] = [];
        for (let r = heap.remove(); r !== null; r = heap.remove()) {
            drained.push(r.value);
        }
        expect(drained).toEqual([20, 30, 100]);
    });

    it('update with an equal value or null record is a no-op', () => {
        const heap = new MinHeap<string>(2);
        const rec = heap.insert('a', 5)!;
        heap.update(rec, 5);
        expect(heap.getMinimum()).toEqual({ key: 'a', value: 5 });
        heap.update(null, 3);
        expect(heap.isValid()).toBe(true);
    });

    it('record identity is stable so keys can be edited in place', () => {
        // The pattern from the upstream polyline-decimation example.
        const heap = new MinHeap<{ label: string }>(3);
        const recA = heap.insert({ label: 'a' }, 3)!;
        heap.insert({ label: 'b' }, 1);
        heap.insert({ label: 'c' }, 2);
        recA.key.label = 'a-edited';
        heap.update(recA, 0);
        expect(heap.getMinimum()!.key.label).toBe('a-edited');
    });

    it('reset clears the heap and restores capacity', () => {
        const heap = new MinHeap<number>(2);
        heap.insert(1, 1);
        heap.insert(2, 2);
        heap.reset(3);
        expect(heap.getNumElements()).toBe(0);
        expect(heap.getMinimum()).toBeNull();
        expect(heap.insert(9, 9)).not.toBeNull();
        expect(heap.insert(8, 8)).not.toBeNull();
        expect(heap.insert(7, 7)).not.toBeNull();
        expect(heap.insert(6, 6)).toBeNull();
        expect(heap.remove()!.value).toBe(7);
    });

    it('accepts a custom comparator (max-heap via reversed order)', () => {
        const heap = new MinHeap<number, number>(4, (a, b) => a > b);
        for (const v of [3, 1, 4, 2]) {
            heap.insert(v, v);
        }
        const drained: number[] = [];
        for (let r = heap.remove(); r !== null; r = heap.remove()) {
            drained.push(r.value);
        }
        expect(drained).toEqual([4, 3, 2, 1]);
    });

    it('maintains the heap invariant under randomized insert/remove/update', () => {
        const rng = makeRng(0x12345678);
        const capacity = 64;
        const heap = new MinHeap<number>(capacity);
        // Model of the live heap contents, keyed by the unique element key.
        const live = new Map<number, { rec: MinHeapRecord<number, number>, value: number }>();
        let nextKey = 0;

        for (let step = 0; step < 3000; ++step) {
            const op = rng();
            if (op < 0.45) {
                const value = Math.floor(rng() * 1000);
                const key = nextKey++;
                const record = heap.insert(key, value);
                if (record !== null) {
                    live.set(key, { rec: record, value });
                } else {
                    expect(live.size).toBe(capacity);
                }
            } else if (op < 0.75) {
                const removed = heap.remove();
                if (live.size === 0) {
                    expect(removed).toBeNull();
                } else {
                    const minValue = Math.min(...[...live.values()].map(e => e.value));
                    expect(removed!.value).toBe(minValue);
                    expect(live.get(removed!.key)!.value).toBe(minValue);
                    live.delete(removed!.key);
                }
            } else if (live.size > 0) {
                const keys = [...live.keys()];
                const key = keys[Math.floor(rng() * keys.length)];
                const entry = live.get(key)!;
                const newValue = Math.floor(rng() * 1000);
                heap.update(entry.rec, newValue);
                entry.value = newValue;
            }

            expect(heap.isValid()).toBe(true);
            expect(heap.getNumElements()).toBe(live.size);
        }

        // Final drain must be sorted and match the model multiset.
        const expected = [...live.values()].map(e => e.value).sort((a, b) => a - b);
        const drained: number[] = [];
        for (let r = heap.remove(); r !== null; r = heap.remove()) {
            drained.push(r.value);
        }
        expect(drained).toEqual(expected);
    });
});
