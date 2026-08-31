import { describe, it, expect } from 'vitest';
import { UniqueVerticesSimplices } from '../src/UniqueVerticesSimplices';

type P2 = [number, number];
type P3 = [number, number, number];

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// Two triangles sharing the edge from (1,0,0) to (0,1,0), as a soup of six
// vertices with two duplicated pairs.
const triangleSoup: P3[] = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0],
    [1, 0, 0], [1, 1, 0], [0, 1, 0]
];

// A polyline of three segments, as a soup of six endpoints. The shared
// endpoints are duplicated.
const edgeSoup: P2[] = [
    [0, 0], [1, 0],
    [1, 0], [1, 1],
    [1, 1], [0, 0]
];

describe('UniqueVerticesSimplices construction', () => {
    it('requires a dimension of at least 2', () => {
        expect(() => new UniqueVerticesSimplices<P3>(1)).toThrow('Invalid dimension.');
        expect(() => new UniqueVerticesSimplices<P3>(0)).toThrow('Invalid dimension.');
        expect(new UniqueVerticesSimplices<P3>(3).dimension).toBe(3);
        expect(new UniqueVerticesSimplices<P2>(2).dimension).toBe(2);
    });

    it('rejects vertex types that the default key function cannot handle', () => {
        const uvs = new UniqueVerticesSimplices<unknown>(2);
        expect(() => uvs.generateIndexedSimplices([() => 0, () => 1]))
            .toThrow('The vertex type requires a keyOf function.');
    });

    it('uses a caller-supplied key function', () => {
        // Round the coordinates so that nearly equal vertices merge.
        const uvs = new UniqueVerticesSimplices<P2>(2,
            v => v.map(c => Math.round(c)).join(','));
        const { vertices, indices } = uvs.generateIndexedSimplices([
            [0, 0], [1.0001, 0], [0.9998, 0], [2, 2]
        ]);
        expect(vertices.length).toBe(3);
        expect(indices).toEqual([0, 1, 1, 2]);
        // The representative kept is the first occurrence.
        expect(vertices[1]).toEqual([1.0001, 0]);
    });

    it('accepts objects with a numeric values array', () => {
        const uvs = new UniqueVerticesSimplices<{ values: number[] }>(2);
        const { vertices, indices } = uvs.generateIndexedSimplices([
            { values: [1, 2] }, { values: [3, 4] },
            { values: [1, 2] }, { values: [5, 6] }
        ]);
        expect(vertices.length).toBe(3);
        expect(indices).toEqual([0, 1, 0, 2]);
    });
});

describe('UniqueVerticesSimplices.generateIndexedSimplices', () => {
    it('builds an indexed triangle mesh from a soup (D = 3)', () => {
        const uvs = new UniqueVerticesSimplices<P3>(3);
        const { vertices, indices } = uvs.generateIndexedSimplices(triangleSoup);

        // The four distinct positions, in order of first occurrence.
        expect(vertices).toEqual([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]);
        expect(indices).toEqual([0, 1, 2, 1, 3, 2]);
        expect(indices.length).toBe(triangleSoup.length);
        // Postcondition 3: the indices are in range.
        for (const index of indices) {
            expect(index).toBeGreaterThanOrEqual(0);
            expect(index).toBeLessThan(vertices.length);
        }
        // The indexed mesh reconstructs the input soup.
        expect(indices.map(i => vertices[i])).toEqual(triangleSoup);
    });

    it('builds an indexed edge mesh from a soup (D = 2)', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const { vertices, indices } = uvs.generateIndexedSimplices(edgeSoup);
        expect(vertices).toEqual([[0, 0], [1, 0], [1, 1]]);
        expect(indices).toEqual([0, 1, 1, 2, 2, 0]);
    });

    it('groups the indices into simplices', () => {
        const uvs = new UniqueVerticesSimplices<P3>(3);
        const { vertices, simplices } =
            uvs.generateIndexedSimplicesGrouped(triangleSoup);
        expect(vertices).toEqual([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]);
        expect(simplices).toEqual([[0, 1, 2], [1, 3, 2]]);
        expect(simplices.length).toBe(triangleSoup.length / 3);
    });

    it('rejects a vertex count that is not a positive multiple of D', () => {
        const uvs = new UniqueVerticesSimplices<P3>(3);
        expect(() => uvs.generateIndexedSimplices([]))
            .toThrow('Invalid number of vertices.');
        expect(() => uvs.generateIndexedSimplices(triangleSoup.slice(0, 4)))
            .toThrow('Invalid number of vertices.');
    });
});

describe('UniqueVerticesSimplices.removeDuplicateVertices', () => {
    // A pool where vertices 0 and 3 coincide, as do 1 and 4.
    const pool: P2[] = [[0, 0], [1, 0], [1, 1], [0, 0], [1, 0]];
    const indices = [0, 1, 1, 2, 3, 4];

    it('renumbers the indices onto the unique pool', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const result = uvs.removeDuplicateVertices(pool, indices);
        expect(result.vertices).toEqual([[0, 0], [1, 0], [1, 1]]);
        expect(result.indices).toEqual([0, 1, 1, 2, 0, 1]);
        // The referenced positions are unchanged.
        expect(result.indices.map(i => result.vertices[i]))
            .toEqual(indices.map(i => pool[i]));
    });

    it('works on grouped simplices', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const simplices = [[0, 1], [1, 2], [3, 4]];
        const result = uvs.removeDuplicateVerticesGrouped(pool, simplices);
        expect(result.vertices).toEqual([[0, 0], [1, 0], [1, 1]]);
        expect(result.simplices).toEqual([[0, 1], [1, 2], [0, 1]]);
    });

    it('is the identity when the pool is already unique', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const unique: P2[] = [[0, 0], [1, 0], [1, 1]];
        const result = uvs.removeDuplicateVertices(unique, [0, 1, 1, 2]);
        expect(result.vertices).toEqual(unique);
        expect(result.indices).toEqual([0, 1, 1, 2]);
    });

    it('validates its inputs', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        expect(() => uvs.removeDuplicateVertices([], [0, 1]))
            .toThrow('Invalid number of vertices.');
        expect(() => uvs.removeDuplicateVertices(pool, []))
            .toThrow('Invalid number of indices.');
        expect(() => uvs.removeDuplicateVertices(pool, [0, 1, 2]))
            .toThrow('Invalid number of indices.');
        expect(() => uvs.removeDuplicateVertices(pool, [0, 5]))
            .toThrow('Invalid index.');
        expect(() => uvs.removeDuplicateVertices(pool, [0, -1]))
            .toThrow('Invalid index.');
        expect(() => uvs.removeDuplicateVerticesGrouped(pool, []))
            .toThrow('Invalid number of simplices.');
        expect(() => uvs.removeDuplicateVerticesGrouped([], [[0, 1]]))
            .toThrow('Invalid number of vertices.');
    });
});

describe('UniqueVerticesSimplices.removeUnusedVertices', () => {
    // Only vertices 1, 3 and 4 are used.
    const pool: P2[] = [[0, 0], [1, 0], [2, 0], [3, 0], [4, 0]];
    const indices = [3, 1, 1, 4];

    it('packs the used vertices in increasing old-index order', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const result = uvs.removeUnusedVertices(pool, indices);
        // std::set<IndexType> visits 1, 3, 4 in increasing order.
        expect(result.vertices).toEqual([[1, 0], [3, 0], [4, 0]]);
        expect(result.indices).toEqual([1, 0, 0, 2]);
        // The referenced positions are unchanged.
        expect(result.indices.map(i => result.vertices[i]))
            .toEqual(indices.map(i => pool[i]));
        // Postcondition 4: every output vertex is used.
        expect(new Set(result.indices).size).toBe(result.vertices.length);
    });

    it('works on grouped simplices', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const result = uvs.removeUnusedVerticesGrouped(pool, [[3, 1], [1, 4]]);
        expect(result.vertices).toEqual([[1, 0], [3, 0], [4, 0]]);
        expect(result.simplices).toEqual([[1, 0], [0, 2]]);
    });

    it('is the identity when every vertex is used', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const result = uvs.removeUnusedVertices(pool, [0, 1, 2, 3, 4, 0]);
        expect(result.vertices).toEqual(pool);
        expect(result.indices).toEqual([0, 1, 2, 3, 4, 0]);
    });

    it('validates its inputs', () => {
        const uvs = new UniqueVerticesSimplices<P2>(2);
        expect(() => uvs.removeUnusedVertices([], [0, 1]))
            .toThrow('Invalid number of vertices.');
        expect(() => uvs.removeUnusedVertices(pool, []))
            .toThrow('Invalid number of indices.');
        expect(() => uvs.removeUnusedVertices(pool, [0, 1, 2]))
            .toThrow('Invalid number of indices.');
        expect(() => uvs.removeUnusedVertices(pool, [0, 9]))
            .toThrow('Invalid index.');
        expect(() => uvs.removeUnusedVerticesGrouped(pool, []))
            .toThrow('Invalid number of simplices.');
        expect(() => uvs.removeUnusedVerticesGrouped(pool, [[0, 1, 2]]))
            .toThrow('Invalid simplex.');
    });
});

describe('UniqueVerticesSimplices.removeDuplicateAndUnusedVertices', () => {
    it('combines the two reductions', () => {
        // Vertices 0 and 4 coincide; vertices 2 and 5 are unused.
        const pool: P2[] = [[0, 0], [1, 0], [9, 9], [1, 1], [0, 0], [8, 8]];
        const indices = [0, 1, 1, 3, 4, 3];
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const result = uvs.removeDuplicateAndUnusedVertices(pool, indices);

        // The pool has no duplicates and no unused vertices.
        expect(result.vertices).toEqual([[0, 0], [1, 0], [1, 1]]);
        expect(new Set(result.vertices.map(v => v.join(','))).size)
            .toBe(result.vertices.length);
        expect(new Set(result.indices).size).toBe(result.vertices.length);
        // The referenced positions are unchanged.
        expect(result.indices.map(i => result.vertices[i]))
            .toEqual(indices.map(i => pool[i]));
    });

    it('combines the two reductions on grouped simplices', () => {
        const pool: P2[] = [[0, 0], [1, 0], [9, 9], [1, 1], [0, 0], [8, 8]];
        const simplices = [[0, 1], [1, 3], [4, 3]];
        const uvs = new UniqueVerticesSimplices<P2>(2);
        const result = uvs.removeDuplicateAndUnusedVerticesGrouped(pool, simplices);
        expect(result.vertices).toEqual([[0, 0], [1, 0], [1, 1]]);
        expect(result.simplices.map(s => s.map(i => result.vertices[i])))
            .toEqual(simplices.map(s => s.map(i => pool[i])));
    });
});

describe('UniqueVerticesSimplices randomized round trips', () => {
    it('preserves the geometry of a random soup for D = 2, 3 and 4', () => {
        const rng = makeRng(0xc0ffee11);
        for (const dimension of [2, 3, 4]) {
            for (let trial = 0; trial < 12; ++trial) {
                const numSimplices = 1 + Math.floor(8 * rng());
                const soup: number[][] = [];
                for (let k = 0; k < dimension * numSimplices; ++k) {
                    // A small coordinate range guarantees many duplicates.
                    soup.push([Math.floor(3 * rng()), Math.floor(3 * rng())]);
                }

                const uvs = new UniqueVerticesSimplices<number[]>(dimension);
                const { vertices, indices } = uvs.generateIndexedSimplices(soup);

                // The indexed representation reproduces the soup exactly.
                expect(indices.length).toBe(soup.length);
                expect(indices.map(i => vertices[i])).toEqual(soup);
                // The output pool has no duplicates.
                expect(new Set(vertices.map(v => v.join(','))).size)
                    .toBe(vertices.length);
                // Every pool vertex is used.
                expect(new Set(indices).size).toBe(vertices.length);

                // The grouped form is the same data, reshaped.
                const grouped = uvs.generateIndexedSimplicesGrouped(soup);
                expect(grouped.vertices).toEqual(vertices);
                expect(([] as number[]).concat(...grouped.simplices)).toEqual(indices);
                expect(grouped.simplices.length).toBe(numSimplices);

                // Reducing an already reduced mesh changes nothing.
                const reduced =
                    uvs.removeDuplicateAndUnusedVertices(vertices, indices);
                expect(reduced.vertices).toEqual(vertices);
                expect(reduced.indices).toEqual(indices);
            }
        }
    });

    it('agrees with an independent reduction of a random indexed mesh', () => {
        const rng = makeRng(0xfeed5678);
        for (let trial = 0; trial < 20; ++trial) {
            const numVertices = 3 + Math.floor(10 * rng());
            const pool: P2[] = [];
            for (let v = 0; v < numVertices; ++v) {
                pool.push([Math.floor(4 * rng()), Math.floor(4 * rng())]);
            }
            const numIndices = 2 * (1 + Math.floor(6 * rng()));
            const indices: number[] = [];
            for (let i = 0; i < numIndices; ++i) {
                indices.push(Math.floor(numVertices * rng()));
            }

            const uvs = new UniqueVerticesSimplices<P2>(2);
            const result = uvs.removeDuplicateAndUnusedVertices(pool, indices);

            // Independent computation: the used positions in the order in
            // which their canonical representative first occurs in the pool.
            const firstOccurrence = new Map<string, number>();
            for (let v = 0; v < pool.length; ++v) {
                const key = pool[v].join(',');
                if (!firstOccurrence.has(key)) {
                    firstOccurrence.set(key, v);
                }
            }
            const usedKeys = new Set<string>(indices.map(i => pool[i].join(',')));
            const expected = Array.from(usedKeys)
                .sort((a, b) =>
                    (firstOccurrence.get(a) as number) -
                    (firstOccurrence.get(b) as number));

            expect(result.vertices.map(v => v.join(','))).toEqual(expected);
            expect(result.indices.map(i => result.vertices[i].join(',')))
                .toEqual(indices.map(i => pool[i].join(',')));
        }
    });
});
