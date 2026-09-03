import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
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

describe('MinHeap verification', () => {
    type Op =
        | { kind: 'insert'; value: number }
        | { kind: 'remove' }
        | { kind: 'update'; which: number; value: number };

    const opArb: fc.Arbitrary<Op> = fc.oneof(
        fc.record({ kind: fc.constant('insert' as const), value: fc.integer({ min: -50, max: 50 }) }),
        fc.record({ kind: fc.constant('remove' as const) }),
        fc.record({
            kind: fc.constant('update' as const),
            which: fc.integer({ min: 0, max: 1000 }),
            value: fc.integer({ min: -50, max: 50 })
        }));

    /** Run a script of operations against the heap and a sorted-array model. */
    function runScript(capacity: number, ops: Op[]) {
        const heap = new MinHeap<number>(capacity);
        const live: { rec: MinHeapRecord<number, number>; value: number; key: number }[] = [];
        let nextKey = 0;
        for (const op of ops) {
            if (op.kind === 'insert') {
                const key = nextKey++;
                const rec = heap.insert(key, op.value);
                if (live.length === capacity) {
                    expect(rec).toBeNull();
                } else {
                    expect(rec).not.toBeNull();
                    live.push({ rec: rec!, value: op.value, key });
                }
            } else if (op.kind === 'remove') {
                const r = heap.remove();
                if (live.length === 0) {
                    expect(r).toBeNull();
                } else {
                    const minValue = live.reduce((m, e) => Math.min(m, e.value), Infinity);
                    expect(r!.value).toBe(minValue);
                    const i = live.findIndex(e => e.key === r!.key);
                    expect(i).toBeGreaterThanOrEqual(0);
                    expect(live[i]!.value).toBe(minValue);
                    live.splice(i, 1);
                }
            } else if (live.length > 0) {
                const e = live[op.which % live.length]!;
                heap.update(e.rec, op.value);
                e.value = op.value;
            }
            expect(heap.isValid()).toBe(true);
            expect(heap.getNumElements()).toBe(live.length);
            if (live.length > 0) {
                const minValue = live.reduce((m, x) => Math.min(m, x.value), Infinity);
                expect(heap.getMinimum()!.value).toBe(minValue);
            } else {
                expect(heap.getMinimum()).toBeNull();
            }
        }
        return { heap, live };
    }

    it('agrees with a sorted-array model under random op scripts', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 12 }),
            fc.array(opArb, { minLength: 0, maxLength: 60 })),
            ([capacity, ops]) => {
                const { heap, live } = runScript(capacity, ops);
                const expected = live.map(e => e.value).sort((a, b) => a - b);
                const drained: number[] = [];
                for (let r = heap.remove(); r !== null; r = heap.remove()) {
                    drained.push(r.value);
                    expect(heap.isValid()).toBe(true);
                }
                expect(drained).toEqual(expected);
                expect(heap.getNumElements()).toBe(0);
            }, 100);
    });

    it('the record index stays correct: any record can be pulled to the root', () => {
        // update() navigates using record.index, so if the index bookkeeping
        // of insert/remove/update were wrong, lowering a record's value below
        // every other value would fail to make it the minimum.
        check(fc.tuple(fc.integer({ min: 1, max: 10 }),
            fc.array(opArb, { minLength: 1, maxLength: 40 })),
            ([capacity, ops]) => {
                const { heap, live } = runScript(capacity, ops);
                for (const e of live) {
                    const original = e.value;
                    heap.update(e.rec, -1000);
                    expect(heap.isValid()).toBe(true);
                    expect(heap.getMinimum()!.key).toBe(e.key);
                    heap.update(e.rec, original);
                    e.value = original;
                    expect(heap.isValid()).toBe(true);
                }
                // Pushing a record to the top of the range keeps the heap
                // valid too (the sift-down branch of update).
                for (const e of live) {
                    const original = e.value;
                    heap.update(e.rec, 1000);
                    expect(heap.isValid()).toBe(true);
                    heap.update(e.rec, original);
                    expect(heap.isValid()).toBe(true);
                }
                return true;
            }, 60);
    });

    it('drains in ascending order for any multiset of values', () => {
        check(fc.array(finite(-100, 100), { minLength: 0, maxLength: 40 }), values => {
            const heap = new MinHeap<number>(values.length);
            values.forEach((v, i) => { heap.insert(i, v); });
            const drained: number[] = [];
            for (let r = heap.remove(); r !== null; r = heap.remove()) { drained.push(r.value); }
            const expected = [...values].sort((a, b) => a - b);
            return drained.length === expected.length
                && drained.every((v, i) => Object.is(v, expected[i]));
        });
    });

    it('a reversed comparator drains in descending order', () => {
        check(fc.array(finite(-100, 100), { minLength: 1, maxLength: 30 }), values => {
            const heap = new MinHeap<number, number>(values.length, (a, b) => b < a);
            values.forEach((v, i) => { heap.insert(i, v); });
            const drained: number[] = [];
            for (let r = heap.remove(); r !== null; r = heap.remove()) { drained.push(r.value); }
            const expected = [...values].sort((a, b) => b - a);
            return drained.every((v, i) => Object.is(v, expected[i]));
        });
    });

    it('a full heap rejects further inserts without disturbing its contents', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 8 }),
            fc.array(finite(-20, 20), { minLength: 0, maxLength: 20 })),
            ([capacity, values]) => {
                const heap = new MinHeap<number>(capacity);
                const accepted: number[] = [];
                values.forEach((v, i) => {
                    const rec = heap.insert(i, v);
                    if (rec === null) {
                        expect(heap.getNumElements()).toBe(capacity);
                    } else {
                        accepted.push(v);
                    }
                });
                expect(accepted.length).toBe(Math.min(capacity, values.length));
                const drained: number[] = [];
                for (let r = heap.remove(); r !== null; r = heap.remove()) { drained.push(r.value); }
                return drained.length === accepted.length
                    && drained.every((v, i) => Object.is(v, accepted.slice().sort((a, b) => a - b)[i]));
            });
    });

    it('update with the same value leaves the heap ordering unchanged', () => {
        check(fc.array(finite(-50, 50), { minLength: 1, maxLength: 20 }), values => {
            const heap = new MinHeap<number>(values.length);
            const recs = values.map((v, i) => heap.insert(i, v)!);
            const before = recs.map(r => r.index);
            for (const r of recs) { heap.update(r, r.value); }
            return heap.isValid() && recs.every((r, i) => r.index === before[i]);
        });
    });

    it('update(null, value) is a no-op', () => {
        check(fc.tuple(fc.array(finite(-50, 50), { minLength: 0, maxLength: 10 }), finite()),
            ([values, v]) => {
                const heap = new MinHeap<number>(values.length);
                values.forEach((x, i) => { heap.insert(i, x); });
                const n = heap.getNumElements();
                heap.update(null, v);
                return heap.getNumElements() === n && heap.isValid();
            });
    });

    it('reset empties the heap and installs the new capacity', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 10 }), fc.integer({ min: 0, max: 10 })),
            ([capacity0, capacity1]) => {
                const heap = new MinHeap<number>(capacity0);
                for (let i = 0; i < capacity0; ++i) { heap.insert(i, i); }
                heap.reset(capacity1);
                if (heap.getNumElements() !== 0 || heap.getMinimum() !== null) { return false; }
                let accepted = 0;
                for (let i = 0; i < capacity1 + 3; ++i) {
                    if (heap.insert(i, -i) !== null) { ++accepted; }
                }
                return accepted === capacity1 && heap.isValid();
            });
    });

    it('an empty heap answers null to every query', () => {
        check(fc.integer({ min: -5, max: 0 }), capacity => {
            const heap = new MinHeap<number>(capacity);
            return heap.getNumElements() === 0 && heap.getMinimum() === null
                && heap.remove() === null && heap.insert(0, 0) === null && heap.isValid();
        });
    });

    it('removed keys are not corrupted by later inserts reusing the record', () => {
        // Upstream copies the key into the caller's variable; the port
        // returns a reference, but insert rebinds record.key rather than
        // mutating it, so an already-returned key object stays intact.
        check(fc.array(finite(-10, 10), { minLength: 2, maxLength: 8 }), values => {
            const heap = new MinHeap<{ id: number }>(values.length);
            values.forEach((v, i) => { heap.insert({ id: i }, v); });
            const first = heap.remove()!;
            const snapshot = first.key.id;
            heap.insert({ id: 999 }, -1000);
            return first.key.id === snapshot;
        });
    });

    it('the sifted element may sit at the last slot during remove (upstream quirk)', () => {
        // Upstream's Remove loop uses 'child <= last', so the slot holding
        // the element being sifted down can be inspected as a child. The
        // comparison is then against the element itself and the loop breaks,
        // which is why the quirk is harmless. This exercises that path.
        const heap = new MinHeap<number>(4);
        heap.insert(0, 0);
        heap.insert(1, 10);
        heap.insert(2, 20);
        heap.insert(3, 5);
        // Heap array is [0, 5, 20, 10]; removing the root moves 10 to the
        // root with last = 3, and the loop reaches child = 3 == last.
        expect(heap.remove()!.value).toBe(0);
        expect(heap.isValid()).toBe(true);
        expect(heap.getMinimum()!.value).toBe(5);
        const drained: number[] = [];
        for (let r = heap.remove(); r !== null; r = heap.remove()) { drained.push(r.value); }
        expect(drained).toEqual([5, 10, 20]);
    });
});
