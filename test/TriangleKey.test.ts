import { describe, it, expect } from 'vitest';
import { TriangleKey } from '../src/TriangleKey.js';
import { FeatureKey } from '../src/FeatureKey.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// All permutations of the input positions (0,1,2).
const permutations: number[][] = [
    [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]
];

// Parity of a permutation, computed as the parity of its inversion count.
// An even permutation has parity 1 and corresponds to a permutation matrix
// with determinant 1; an odd permutation has parity -1.
function parity(perm: readonly number[]): number {
    let inversions = 0;
    for (let i = 0; i < perm.length; ++i) {
        for (let j = i + 1; j < perm.length; ++j) {
            if (perm[i] > perm[j]) {
                ++inversions;
            }
        }
    }
    return (inversions % 2 === 0) ? 1 : -1;
}

// The permutation of the positions of 'input' that produces 'output'. The
// values of 'input' must be distinct.
function permutationOf(input: readonly number[],
    output: readonly number[]): number[] {
    return [...output].map(value => input.indexOf(value));
}

describe('TriangleKey construction', () => {
    it('initializes to invalid indices when the vertices are omitted', () => {
        for (const ordered of [false, true]) {
            const key = new TriangleKey(ordered);
            expect(key.V).toEqual([-1, -1, -1]);
            expect(key.n).toBe(3);
            expect(key.ordered).toBe(ordered);
        }
    });

    it('is a FeatureKey with three vertex indices', () => {
        const key = new TriangleKey(true, 5, 1, 9);
        expect(key).toBeInstanceOf(FeatureKey);
        expect(key.n).toBe(3);
    });
});

describe('TriangleKey ordered', () => {
    it('matches the documented canonical forms', () => {
        // (v0, v1, v2) = (4, 7, 9); the cyclic permutations all store
        // (4, 7, 9) and the transpositions all store (4, 9, 7).
        expect(new TriangleKey(true, 4, 7, 9).V).toEqual([4, 7, 9]);
        expect(new TriangleKey(true, 7, 9, 4).V).toEqual([4, 7, 9]);
        expect(new TriangleKey(true, 9, 4, 7).V).toEqual([4, 7, 9]);
        expect(new TriangleKey(true, 4, 9, 7).V).toEqual([4, 9, 7]);
        expect(new TriangleKey(true, 9, 7, 4).V).toEqual([4, 9, 7]);
        expect(new TriangleKey(true, 7, 4, 9).V).toEqual([4, 9, 7]);
    });

    it('stores a cyclic permutation beginning with the minimum', () => {
        const rng = makeRng(9001);
        for (let trial = 0; trial < 200; ++trial) {
            const values = [
                Math.floor(rng() * 1000),
                Math.floor(rng() * 1000),
                Math.floor(rng() * 1000)
            ];
            if (new Set(values).size !== 3) {
                continue;  // require distinct indices for the parity test
            }
            for (const perm of permutations) {
                const input = perm.map(i => values[i]);
                const key = new TriangleKey(true, input[0], input[1], input[2]);

                // V[0] is the minimum of the inputs.
                expect(key.V[0]).toBe(Math.min(...values));

                // V is a permutation of the inputs.
                expect([...key.V].sort((a, b) => a - b))
                    .toEqual([...values].sort((a, b) => a - b));

                // The permutation taking the input to V is even, so the
                // induced permutation matrix has determinant 1.
                expect(parity(permutationOf(input, key.V))).toBe(1);
            }
        }
    });

    it('groups the six input orderings into two orientation classes', () => {
        const values = [11, 3, 25];
        const byParity = new Map<number, Set<string>>();
        for (const perm of permutations) {
            const input = perm.map(i => values[i]);
            const key = new TriangleKey(true, input[0], input[1], input[2]);
            const p = parity(perm);
            if (!byParity.has(p)) {
                byParity.set(p, new Set<string>());
            }
            byParity.get(p)!.add(key.mapKey());
        }
        // Each orientation class collapses to one canonical key, and the two
        // classes are distinct.
        expect(byParity.get(1)!.size).toBe(1);
        expect(byParity.get(-1)!.size).toBe(1);
        const even = [...byParity.get(1)!][0];
        const odd = [...byParity.get(-1)!][0];
        expect(even).not.toBe(odd);
    });

    it('handles repeated indices', () => {
        expect(new TriangleKey(true, 5, 5, 3).V).toEqual([3, 5, 5]);
        expect(new TriangleKey(true, 5, 3, 5).V).toEqual([3, 5, 5]);
        expect(new TriangleKey(true, 3, 5, 5).V).toEqual([3, 5, 5]);
        expect(new TriangleKey(true, 2, 2, 2).V).toEqual([2, 2, 2]);
    });
});

describe('TriangleKey unordered', () => {
    it('sorts the indices increasingly for every input ordering', () => {
        const rng = makeRng(4242);
        for (let trial = 0; trial < 300; ++trial) {
            const values = [
                Math.floor(rng() * 40),
                Math.floor(rng() * 40),
                Math.floor(rng() * 40)
            ];
            const expected = [...values].sort((a, b) => a - b);
            for (const perm of permutations) {
                const input = perm.map(i => values[i]);
                const key = new TriangleKey(false, input[0], input[1], input[2]);
                expect(key.V).toEqual(expected);
            }
        }
    });

    it('collapses all six orderings to one map key', () => {
        const keys = new Set<string>();
        for (const perm of permutations) {
            const input = perm.map(i => [6, 2, 8][i]);
            keys.add(new TriangleKey(false, input[0], input[1], input[2])
                .mapKey());
        }
        expect(keys.size).toBe(1);
        expect([...keys][0]).toBe('2,6,8');
    });

    it('handles repeated and negative indices', () => {
        expect(new TriangleKey(false, 5, 5, 3).V).toEqual([3, 5, 5]);
        expect(new TriangleKey(false, 3, 5, 5).V).toEqual([3, 5, 5]);
        expect(new TriangleKey(false, -1, -4, 2).V).toEqual([-4, -1, 2]);
        expect(new TriangleKey(false, 7, 7, 7).V).toEqual([7, 7, 7]);
    });
});
