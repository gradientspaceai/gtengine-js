import { describe, it, expect } from 'vitest';
import { TetrahedronKey } from '../src/TetrahedronKey.js';
import { FeatureKey } from '../src/FeatureKey.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// All 24 permutations of the input positions (0,1,2,3).
function allPermutations(items: readonly number[]): number[][] {
    if (items.length <= 1) {
        return [[...items]];
    }
    const result: number[][] = [];
    for (let i = 0; i < items.length; ++i) {
        const rest = [...items.slice(0, i), ...items.slice(i + 1)];
        for (const tail of allPermutations(rest)) {
            result.push([items[i], ...tail]);
        }
    }
    return result;
}

const permutations = allPermutations([0, 1, 2, 3]);

// Parity of a permutation, computed as the parity of its inversion count.
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

function permutationOf(input: readonly number[],
    output: readonly number[]): number[] {
    return [...output].map(value => input.indexOf(value));
}

describe('TetrahedronKey construction', () => {
    it('initializes to invalid indices when the vertices are omitted', () => {
        for (const ordered of [false, true]) {
            const key = new TetrahedronKey(ordered);
            expect(key.V).toEqual([-1, -1, -1, -1]);
            expect(key.n).toBe(4);
            expect(key.ordered).toBe(ordered);
        }
    });

    it('is a FeatureKey with four vertex indices', () => {
        const key = new TetrahedronKey(true, 5, 1, 9, 3);
        expect(key).toBeInstanceOf(FeatureKey);
        expect(key.n).toBe(4);
    });
});

describe('TetrahedronKey ordered', () => {
    it('stores one of the twelve documented forms of the input', () => {
        // The header documents that for inputs (v0,v1,v2,v3) the storage
        // (V[0],V[1],V[2],V[3]) is one of these tuples of input slots.
        const documented: string[] = [
            [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2],
            [1, 3, 2, 0], [1, 2, 0, 3], [1, 0, 3, 2],
            [2, 3, 0, 1], [2, 0, 1, 3], [2, 1, 3, 0],
            [3, 1, 0, 2], [3, 0, 2, 1], [3, 2, 1, 0]
        ].map(form => form.join(','));

        const observed = new Set<string>();
        const rng = makeRng(2468);
        for (let trial = 0; trial < 200; ++trial) {
            const values = [0, 0, 0, 0].map(() => Math.floor(rng() * 1000));
            if (new Set(values).size !== 4) {
                continue;
            }
            for (const perm of permutations) {
                const input = perm.map(i => values[i]);
                const key = new TetrahedronKey(true, input[0], input[1],
                    input[2], input[3]);
                const form = permutationOf(input, key.V).join(',');
                expect(documented).toContain(form);
                observed.add(form);
            }
        }
        // Every documented form is realized.
        expect(observed.size).toBe(12);
    });

    it('stores V[0] = min and V[1] = min of the remaining indices', () => {
        const rng = makeRng(777);
        for (let trial = 0; trial < 100; ++trial) {
            const values = [0, 0, 0, 0].map(() => Math.floor(rng() * 1000));
            if (new Set(values).size !== 4) {
                continue;
            }
            for (const perm of permutations) {
                const input = perm.map(i => values[i]);
                const key = new TetrahedronKey(true, input[0], input[1],
                    input[2], input[3]);
                const sorted = [...values].sort((a, b) => a - b);
                expect(key.V[0]).toBe(sorted[0]);
                expect(key.V[1]).toBe(sorted[1]);
                // V is a permutation of the inputs...
                expect([...key.V].sort((a, b) => a - b)).toEqual(sorted);
                // ... and an even one, so the induced permutation matrix has
                // determinant 1.
                expect(parity(permutationOf(input, key.V))).toBe(1);
            }
        }
    });

    it('groups the 24 input orderings into two orientation classes', () => {
        const values = [17, 4, 31, 9];
        const evenKeys = new Set<string>();
        const oddKeys = new Set<string>();
        for (const perm of permutations) {
            const input = perm.map(i => values[i]);
            const key = new TetrahedronKey(true, input[0], input[1], input[2],
                input[3]);
            (parity(perm) === 1 ? evenKeys : oddKeys).add(key.mapKey());
        }
        expect(evenKeys.size).toBe(1);
        expect(oddKeys.size).toBe(1);
        expect([...evenKeys][0]).not.toBe([...oddKeys][0]);
    });

    it('handles repeated indices', () => {
        expect(new TetrahedronKey(true, 5, 5, 5, 5).V).toEqual([5, 5, 5, 5]);
        // The minimum search uses strict comparisons, so the first
        // occurrence of a repeated minimum is the one chosen for V[0].
        expect(new TetrahedronKey(true, 2, 2, 7, 9).V).toEqual([2, 2, 7, 9]);
    });
});

describe('TetrahedronKey unordered', () => {
    it('sorts the indices increasingly for every input ordering', () => {
        const rng = makeRng(31337);
        for (let trial = 0; trial < 100; ++trial) {
            const values = [0, 0, 0, 0].map(() => Math.floor(rng() * 30));
            const expected = [...values].sort((a, b) => a - b);
            for (const perm of permutations) {
                const input = perm.map(i => values[i]);
                const key = new TetrahedronKey(false, input[0], input[1],
                    input[2], input[3]);
                expect(key.V).toEqual(expected);
            }
        }
    });

    it('sorts numerically, not lexicographically as strings', () => {
        // The default JavaScript Array.sort would produce [10,2,3,9].
        expect(new TetrahedronKey(false, 10, 2, 9, 3).V).toEqual([2, 3, 9, 10]);
        expect(new TetrahedronKey(false, -2, 11, -30, 5).V)
            .toEqual([-30, -2, 5, 11]);
    });

    it('collapses all 24 orderings to one map key', () => {
        const keys = new Set<string>();
        for (const perm of permutations) {
            const input = perm.map(i => [6, 2, 8, 1][i]);
            keys.add(new TetrahedronKey(false, input[0], input[1], input[2],
                input[3]).mapKey());
        }
        expect(keys.size).toBe(1);
        expect([...keys][0]).toBe('1,2,6,8');
    });
});

describe('TetrahedronKey.getOppositeFace', () => {
    const oppositeFace = TetrahedronKey.getOppositeFace();

    it('has the upstream table values', () => {
        expect(oppositeFace.map(face => [...face])).toEqual([
            [1, 2, 3],
            [0, 3, 2],
            [0, 1, 3],
            [0, 2, 1]
        ]);
    });

    it('returns the same table on each call', () => {
        expect(TetrahedronKey.getOppositeFace()).toBe(oppositeFace);
    });

    it('lists the three vertices other than the opposite one', () => {
        for (let j = 0; j < 4; ++j) {
            const face = [...oppositeFace[j]];
            expect(face).not.toContain(j);
            expect([...face].sort((a, b) => a - b))
                .toEqual([0, 1, 2, 3].filter(i => i !== j));
        }
    });

    it('orders each face counterclockwise viewed from outside', () => {
        // A positively oriented tetrahedron: the determinant of
        // [p1-p0, p2-p0, p3-p0] is +1.
        const p = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([0, 1, 0]),
            Vector.fromArray([0, 0, 1])
        ];
        for (let j = 0; j < 4; ++j) {
            const [a, b, c] = oppositeFace[j];
            const normal = cross(sub(p[b], p[a]), sub(p[c], p[a]));
            // The normal must point away from the opposite vertex p[j].
            expect(dot(normal, sub(p[a], p[j]))).toBeGreaterThan(0);
        }
    });
});
