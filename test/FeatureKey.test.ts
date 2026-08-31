import { describe, it, expect } from 'vitest';
import { FeatureKey } from '../src/FeatureKey';
import { hashCombine } from '../src/HashCombine';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

function makeKey(values: number[], ordered: boolean = false): FeatureKey {
    const key = new FeatureKey(values.length, ordered);
    for (let i = 0; i < values.length; ++i) {
        key.V[i] = values[i];
    }
    return key;
}

describe('FeatureKey construction', () => {
    it('zero-fills V and records n and ordered', () => {
        const key = new FeatureKey(3, true);
        expect(key.V).toEqual([0, 0, 0]);
        expect(key.n).toBe(3);
        expect(key.ordered).toBe(true);

        const unordered = new FeatureKey(2, false);
        expect(unordered.V).toEqual([0, 0]);
        expect(unordered.n).toBe(2);
        expect(unordered.ordered).toBe(false);
    });

    it('rejects a nonpositive number of vertex indices', () => {
        expect(() => new FeatureKey(0, false)).toThrow(
            'The number of vertex indices must be positive.');
        expect(() => new FeatureKey(-1, false)).toThrow();
    });
});

describe('FeatureKey comparisons', () => {
    it('compares equal keys as equal', () => {
        const key0 = makeKey([2, 5, 7]);
        const key1 = makeKey([2, 5, 7]);
        expect(key0.equals(key1)).toBe(true);
        expect(key0.notEqual(key1)).toBe(false);
        expect(key0.lessThan(key1)).toBe(false);
        expect(key0.greaterThan(key1)).toBe(false);
        expect(key0.lessThanOrEqual(key1)).toBe(true);
        expect(key0.greaterThanOrEqual(key1)).toBe(true);
    });

    it('orders keys lexicographically', () => {
        // The first differing component decides the order.
        const key0 = makeKey([1, 9, 9]);
        const key1 = makeKey([2, 0, 0]);
        expect(key0.lessThan(key1)).toBe(true);
        expect(key1.lessThan(key0)).toBe(false);
        expect(key0.greaterThan(key1)).toBe(false);
        expect(key1.greaterThan(key0)).toBe(true);

        // Ties in the leading components are broken by later components.
        const key2 = makeKey([3, 4, 5]);
        const key3 = makeKey([3, 4, 6]);
        expect(key2.lessThan(key3)).toBe(true);
        expect(key3.lessThan(key2)).toBe(false);
        expect(key2.notEqual(key3)).toBe(true);
    });

    it('satisfies the strict-weak-ordering identities on random keys', () => {
        const rng = makeRng(0x5eed1234);
        const keys: FeatureKey[] = [];
        for (let k = 0; k < 60; ++k) {
            keys.push(makeKey([
                Math.floor(3 * rng()),
                Math.floor(3 * rng()),
                Math.floor(3 * rng())
            ]));
        }

        for (const a of keys) {
            // Irreflexivity of '<'.
            expect(a.lessThan(a)).toBe(false);
            for (const b of keys) {
                // Trichotomy.
                const lt = a.lessThan(b);
                const gt = a.greaterThan(b);
                const eq = a.equals(b);
                expect([lt, gt, eq].filter(x => x).length).toBe(1);
                // The derived operators are defined from '<'.
                expect(a.lessThanOrEqual(b)).toBe(lt || eq);
                expect(a.greaterThanOrEqual(b)).toBe(gt || eq);
                expect(a.notEqual(b)).toBe(!eq);
                // Symmetry of the reversed comparison.
                expect(b.lessThan(a)).toBe(gt);
            }
        }
    });
});

describe('FeatureKey hashing and map support', () => {
    it('matches the upstream fold of HashCombine over V', () => {
        const key = makeKey([4, 8, 15]);
        let expected = 0;
        for (const v of key.V) {
            expected = hashCombine(expected, v);
        }
        expect(key.hashValue()).toBe(expected);
        expect(FeatureKey.hashValue(key)).toBe(expected);
    });

    it('gives equal keys equal hashes', () => {
        const rng = makeRng(24680);
        for (let k = 0; k < 50; ++k) {
            const values = [
                Math.floor(100 * rng()),
                Math.floor(100 * rng()),
                Math.floor(100 * rng())
            ];
            const key0 = makeKey(values);
            const key1 = makeKey(values.slice());
            expect(key0.equals(key1)).toBe(true);
            expect(key0.hashValue()).toBe(key1.hashValue());
            expect(FeatureKey.equal(key0, key1)).toBe(true);
        }
    });

    it('produces mapKey values that agree with equals', () => {
        const rng = makeRng(13579);
        const keys: FeatureKey[] = [];
        for (let k = 0; k < 80; ++k) {
            keys.push(makeKey([
                Math.floor(4 * rng()),
                Math.floor(4 * rng()),
                Math.floor(4 * rng())
            ]));
        }
        for (const a of keys) {
            for (const b of keys) {
                expect(a.mapKey() === b.mapKey()).toBe(a.equals(b));
            }
        }

        // Keys of different lengths never collide.
        expect(makeKey([1, 2]).mapKey()).not.toBe(makeKey([1, 2, 0]).mapKey());
        expect(makeKey([1]).mapKey()).not.toBe(makeKey([1, 2]).mapKey());
    });

    it('reports distinct hash values for typical distinct keys', () => {
        // Not a requirement of the algorithm, but a sanity check that the
        // hash is sensitive to the order of the components.
        expect(makeKey([1, 2, 3]).hashValue())
            .not.toBe(makeKey([3, 2, 1]).hashValue());
    });
});

describe('FeatureKey.compare', () => {
    it('sorts keys in the order that std::set<FeatureKey> iterates', () => {
        const keys = [
            makeKey([2, 0, 1]),
            makeKey([0, 3, 4]),
            makeKey([2, 0, 0]),
            makeKey([0, 1, 9]),
            makeKey([1, 1, 1])
        ];
        keys.sort(FeatureKey.compare);
        expect(keys.map(k => k.V)).toEqual([
            [0, 1, 9],
            [0, 3, 4],
            [1, 1, 1],
            [2, 0, 0],
            [2, 0, 1]
        ]);

        // The sort is consistent with lessThan on adjacent elements.
        for (let i = 0; i + 1 < keys.length; ++i) {
            expect(keys[i].lessThan(keys[i + 1])).toBe(true);
        }
    });

    it('returns 0 exactly for equal keys', () => {
        expect(FeatureKey.compare(makeKey([5, 6]), makeKey([5, 6]))).toBe(0);
        expect(FeatureKey.compare(makeKey([5, 6]), makeKey([5, 7]))).toBe(-1);
        expect(FeatureKey.compare(makeKey([5, 7]), makeKey([5, 6]))).toBe(1);
    });
});
