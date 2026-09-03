import { describe, it, expect } from 'vitest';
import { EdgeKey } from '../src/EdgeKey.js';
import { FeatureKey } from '../src/FeatureKey.js';

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
