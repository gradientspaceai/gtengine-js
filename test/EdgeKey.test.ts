import { describe, it, expect } from 'vitest';
import { EdgeKey } from '../src/EdgeKey.js';
import { FeatureKey } from '../src/FeatureKey.js';
import { check, fc } from './helpers/arbitraries.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

describe('EdgeKey construction', () => {
    it('initializes to invalid indices when the vertices are omitted', () => {
        for (const ordered of [false, true]) {
            const key = new EdgeKey(ordered);
            expect(key.V).toEqual([-1, -1]);
            expect(key.n).toBe(2);
            expect(key.ordered).toBe(ordered);
        }
    });

    it('is a FeatureKey with two vertex indices', () => {
        const key = new EdgeKey(true, 3, 7);
        expect(key).toBeInstanceOf(FeatureKey);
        expect(key.n).toBe(2);
    });
});

describe('EdgeKey ordered', () => {
    it('stores (v0, v1) exactly as given', () => {
        expect(new EdgeKey(true, 3, 7).V).toEqual([3, 7]);
        expect(new EdgeKey(true, 7, 3).V).toEqual([7, 3]);
        expect(new EdgeKey(true, -5, 0).V).toEqual([-5, 0]);
        expect(new EdgeKey(true, 4, 4).V).toEqual([4, 4]);
    });

    it('distinguishes the two orientations of an edge', () => {
        const forward = new EdgeKey(true, 2, 9);
        const reverse = new EdgeKey(true, 9, 2);
        expect(forward.equals(reverse)).toBe(false);
        expect(forward.mapKey()).not.toBe(reverse.mapKey());
        // The reverse of the reverse is the original.
        const twice = new EdgeKey(true, reverse.V[1], reverse.V[0]);
        expect(twice.equals(forward)).toBe(true);
    });
});

describe('EdgeKey unordered', () => {
    it('stores (min, max)', () => {
        expect(new EdgeKey(false, 3, 7).V).toEqual([3, 7]);
        expect(new EdgeKey(false, 7, 3).V).toEqual([3, 7]);
        expect(new EdgeKey(false, 0, -5).V).toEqual([-5, 0]);
        expect(new EdgeKey(false, 4, 4).V).toEqual([4, 4]);
    });

    it('is invariant under swapping the inputs (randomized)', () => {
        const rng = makeRng(12345);
        for (let trial = 0; trial < 200; ++trial) {
            const v0 = Math.floor(rng() * 50);
            const v1 = Math.floor(rng() * 50);
            const key0 = new EdgeKey(false, v0, v1);
            const key1 = new EdgeKey(false, v1, v0);
            expect(key0.V).toEqual(key1.V);
            expect(key0.equals(key1)).toBe(true);
            expect(key0.V[0]).toBe(Math.min(v0, v1));
            expect(key0.V[1]).toBe(Math.max(v0, v1));
            // The stored indices are a permutation of the inputs.
            expect([...key0.V].sort((a, b) => a - b))
                .toEqual([v0, v1].sort((a, b) => a - b));
        }
    });
});

describe('EdgeKey use as a map key', () => {
    it('collapses the two orientations only when unordered', () => {
        const unordered = new Set<string>();
        unordered.add(new EdgeKey(false, 1, 4).mapKey());
        unordered.add(new EdgeKey(false, 4, 1).mapKey());
        expect(unordered.size).toBe(1);

        const ordered = new Set<string>();
        ordered.add(new EdgeKey(true, 1, 4).mapKey());
        ordered.add(new EdgeKey(true, 4, 1).mapKey());
        expect(ordered.size).toBe(2);
    });

    it('hashes equal keys to the same value', () => {
        const key0 = new EdgeKey(false, 8, 2);
        const key1 = new EdgeKey(false, 2, 8);
        expect(FeatureKey.equal(key0, key1)).toBe(true);
        expect(key0.hashValue()).toBe(key1.hashValue());
    });
});

describe('EdgeKey verification', () => {
    const idx = fc.integer({ min: -2, max: 5 });

    it('ordered keys store the inputs verbatim', () => {
        check(fc.tuple(idx, idx), ([v0, v1]) => {
            const key = new EdgeKey(true, v0, v1);
            expect(key.V).toEqual([v0, v1]);
            expect(key.ordered).toBe(true);
            expect(key.n).toBe(2);
        });
    });

    it('unordered keys are sorted and swap-invariant', () => {
        check(fc.tuple(idx, idx), ([v0, v1]) => {
            const a = new EdgeKey(false, v0, v1);
            const b = new EdgeKey(false, v1, v0);
            expect(a.V[0]).toBeLessThanOrEqual(a.V[1]);
            expect(a.V[0]).toBe(Math.min(v0, v1));
            expect(a.V[1]).toBe(Math.max(v0, v1));
            expect(a.equals(b)).toBe(true);
            expect(a.mapKey()).toBe(b.mapKey());
            expect(a.hashValue()).toBe(b.hashValue());
        });
    });

    it('ordered keys distinguish the two orientations', () => {
        check(fc.tuple(idx, idx).filter(([v0, v1]) => v0 !== v1),
            ([v0, v1]) => {
                const a = new EdgeKey(true, v0, v1);
                const b = new EdgeKey(true, v1, v0);
                expect(a.equals(b)).toBe(false);
                // The unordered key of either orientation is the same.
                const u0 = new EdgeKey(false, v0, v1);
                const u1 = new EdgeKey(false, v1, v0);
                expect(u0.equals(u1)).toBe(true);
            });
    });

    it('a Set of unordered edges counts undirected edges', () => {
        check(fc.array(fc.tuple(idx, idx).filter(([a, b]) => a !== b),
            { minLength: 1, maxLength: 12 }), (pairs) => {
            const set = new Set<string>();
            for (const [v0, v1] of pairs) {
                set.add(new EdgeKey(false, v0, v1).mapKey());
            }
            const brute = new Set<string>();
            for (const [v0, v1] of pairs) {
                brute.add(`${Math.min(v0, v1)},${Math.max(v0, v1)}`);
            }
            expect(set.size).toBe(brute.size);
        });
    });
});
