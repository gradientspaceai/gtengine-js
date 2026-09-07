import { describe, it, expect } from 'vitest';
import { UniqueVerticesSimplices } from '../src/UniqueVerticesSimplices.js';
import { Vector } from '../src/Vector.js';
import { check, fc } from './helpers/arbitraries.js';

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
    it('copies vertices that expose clone() instead of aliasing the input', () => {
        // Upstream's std::vector<VertexType> holds copies; a caller mutating
        // its input must not see the change in the unique-vertex array.
        const a = Vector.fromArray([1, 2, 3]);
        const b = Vector.fromArray([4, 5, 6]);
        const uvs = new UniqueVerticesSimplices<Vector>(3);
        const { vertices } = uvs.generateIndexedSimplices([a, b, a]);
        expect(vertices.length).toBe(2);
        expect(vertices[0]).not.toBe(a);
        a.values[0] = 99;
        expect(vertices[0].values[0]).toBe(1);
    });

});

describe('UniqueVerticesSimplices verification', () => {
    // Small integer vertices so duplicates are frequent.
    const vertexArb = fc.tuple(fc.integer({ min: 0, max: 3 }),
        fc.integer({ min: 0, max: 3 }), fc.integer({ min: 0, max: 3 }))
        .map(([x, y, z]) => [x, y, z] as P3);

    const soupArb = fc.array(vertexArb, { minLength: 3, maxLength: 21 })
        .map((v) => v.slice(0, 3 * Math.floor(v.length / 3)))
        .filter((v) => v.length > 0);

    const key = (v: P3): string => v.join(',');

    const indexedArb = fc.tuple(
        fc.array(vertexArb, { minLength: 1, maxLength: 8 }),
        fc.array(fc.nat(), { minLength: 3, maxLength: 18 }))
        .map(([pool, raw]) => {
            const indices = raw.slice(0, 3 * Math.floor(raw.length / 3))
                .map((i) => i % pool.length);
            return { pool, indices };
        })
        .filter(({ indices }) => indices.length > 0);

    it('generateIndexedSimplices reproduces the input soup exactly', () => {
        check(soupArb, (soup) => {
            const uvs = new UniqueVerticesSimplices<P3>(3);
            const { vertices, indices } = uvs.generateIndexedSimplices(soup);
            expect(indices.length).toBe(soup.length);
            for (let i = 0; i < soup.length; ++i) {
                expect(vertices[indices[i]]).toEqual(soup[i]);
                expect(indices[i]).toBeGreaterThanOrEqual(0);
                expect(indices[i]).toBeLessThan(vertices.length);
            }
            // The pool is unique and in first-encounter order.
            expect(new Set(vertices.map(key)).size).toBe(vertices.length);
            const firstSeen: string[] = [];
            const seen = new Set<string>();
            for (const v of soup) {
                if (!seen.has(key(v))) {
                    seen.add(key(v));
                    firstSeen.push(key(v));
                }
            }
            expect(vertices.map(key)).toEqual(firstSeen);
        });
    });

    it('the grouped form is the flat form cut into simplices', () => {
        check(soupArb, (soup) => {
            const uvs = new UniqueVerticesSimplices<P3>(3);
            const flat = uvs.generateIndexedSimplices(soup);
            const grouped = uvs.generateIndexedSimplicesGrouped(soup);
            expect(grouped.vertices).toEqual(flat.vertices);
            expect(grouped.simplices.length).toBe(soup.length / 3);
            expect(grouped.simplices.flat()).toEqual(flat.indices);
        });
    });

    it('removeDuplicateVertices preserves the geometry of every index', () => {
        check(indexedArb, ({ pool, indices }) => {
            const uvs = new UniqueVerticesSimplices<P3>(3);
            const out = uvs.removeDuplicateVertices(pool, indices);
            expect(out.indices.length).toBe(indices.length);
            expect(new Set(out.vertices.map(key)).size).toBe(out.vertices.length);
            for (let i = 0; i < indices.length; ++i) {
                expect(out.vertices[out.indices[i]]).toEqual(pool[indices[i]]);
            }
            // Running it again is the identity on an already-unique pool.
            const again = uvs.removeDuplicateVertices(out.vertices, out.indices);
            expect(again.vertices).toEqual(out.vertices);
            expect(again.indices).toEqual(out.indices);
        });
    });

    it('removeUnusedVertices keeps used vertices in old-index order', () => {
        check(indexedArb, ({ pool, indices }) => {
            const uvs = new UniqueVerticesSimplices<P3>(3);
            const out = uvs.removeUnusedVertices(pool, indices);
            const used = Array.from(new Set(indices)).sort((a, b) => a - b);
            expect(out.vertices.length).toBe(used.length);
            expect(out.vertices).toEqual(used.map((i) => pool[i]));
            for (let i = 0; i < indices.length; ++i) {
                expect(out.vertices[out.indices[i]]).toEqual(pool[indices[i]]);
            }
            // Postcondition 4: every output vertex is referenced.
            expect(new Set(out.indices).size).toBe(out.vertices.length);
            // Idempotent once nothing is unused.
            const again = uvs.removeUnusedVertices(out.vertices, out.indices);
            expect(again.vertices).toEqual(out.vertices);
            expect(again.indices).toEqual(out.indices);
        });
    });

    it('removeDuplicateAndUnusedVertices equals the composition', () => {
        check(indexedArb, ({ pool, indices }) => {
            const uvs = new UniqueVerticesSimplices<P3>(3);
            const temp = uvs.removeDuplicateVertices(pool, indices);
            const expected = uvs.removeUnusedVertices(temp.vertices, temp.indices);
            const actual = uvs.removeDuplicateAndUnusedVertices(pool, indices);
            expect(actual.vertices).toEqual(expected.vertices);
            expect(actual.indices).toEqual(expected.indices);
            // The result is fully compacted: unique and fully used.
            expect(new Set(actual.vertices.map(key)).size)
                .toBe(actual.vertices.length);
            expect(new Set(actual.indices).size).toBe(actual.vertices.length);
            for (let i = 0; i < indices.length; ++i) {
                expect(actual.vertices[actual.indices[i]])
                    .toEqual(pool[indices[i]]);
            }
        });
    });

    it('the grouped entry points agree with the flat ones', () => {
        check(indexedArb, ({ pool, indices }) => {
            const uvs = new UniqueVerticesSimplices<P3>(3);
            const simplices: number[][] = [];
            for (let i = 0; i < indices.length; i += 3) {
                simplices.push(indices.slice(i, i + 3));
            }

            const dupFlat = uvs.removeDuplicateVertices(pool, indices);
            const dupGrouped = uvs.removeDuplicateVerticesGrouped(pool, simplices);
            expect(dupGrouped.vertices).toEqual(dupFlat.vertices);
            expect(dupGrouped.simplices.flat()).toEqual(dupFlat.indices);

            const unusedFlat = uvs.removeUnusedVertices(pool, indices);
            const unusedGrouped = uvs.removeUnusedVerticesGrouped(pool, simplices);
            expect(unusedGrouped.vertices).toEqual(unusedFlat.vertices);
            expect(unusedGrouped.simplices.flat()).toEqual(unusedFlat.indices);

            const bothFlat = uvs.removeDuplicateAndUnusedVertices(pool, indices);
            const bothGrouped =
                uvs.removeDuplicateAndUnusedVerticesGrouped(pool, simplices);
            expect(bothGrouped.vertices).toEqual(bothFlat.vertices);
            expect(bothGrouped.simplices.flat()).toEqual(bothFlat.indices);
        });
    });

    it('no entry point aliases the caller vertex objects', () => {
        // Upstream stores VertexType by value in std::map and in the output
        // std::vector, so the outputs are copies. This is a regression test
        // for removeUnusedVertices, which packed references before.
        check(indexedArb, ({ pool, indices }) => {
            const uvs = new UniqueVerticesSimplices<P3>(3);
            const results = [
                uvs.removeDuplicateVertices(pool, indices).vertices,
                uvs.removeUnusedVertices(pool, indices).vertices,
                uvs.removeDuplicateAndUnusedVertices(pool, indices).vertices,
                uvs.generateIndexedSimplices(indices.map((i) => pool[i])).vertices
            ];
            for (const vertices of results) {
                for (const v of vertices) {
                    for (const input of pool) {
                        expect(v).not.toBe(input);
                    }
                }
            }
        });
    });

    it('mutating the input after the call leaves the output unchanged', () => {
        const pool: P3[] = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [9, 9, 9]];
        const indices = [0, 1, 2];
        const uvs = new UniqueVerticesSimplices<P3>(3);
        const unused = uvs.removeUnusedVertices(pool, indices);
        const duplicates = uvs.removeDuplicateVertices(pool, indices);
        pool[1][0] = 42;
        expect(unused.vertices[1]).toEqual([1, 0, 0]);
        expect(duplicates.vertices[1]).toEqual([1, 0, 0]);
    });

    it('vertices that expose clone() are copied by every entry point', () => {
        const pool = [Vector.fromArray([1, 2, 3]), Vector.fromArray([4, 5, 6])];
        const indices = [0, 1, 0];
        const uvs = new UniqueVerticesSimplices<Vector>(3);
        const unused = uvs.removeUnusedVertices(pool, indices);
        expect(unused.vertices[0]).not.toBe(pool[0]);
        pool[0].values[0] = 99;
        expect(unused.vertices[0].values[0]).toBe(1);
    });

    it('immutable vertex types pass through unchanged', () => {
        check(fc.array(fc.integer({ min: 0, max: 4 }),
            { minLength: 2, maxLength: 12 })
            .map((v) => v.slice(0, 2 * Math.floor(v.length / 2)))
            .filter((v) => v.length > 0), (soup) => {
            const uvs = new UniqueVerticesSimplices<number>(2);
            const { vertices, indices } = uvs.generateIndexedSimplices(soup);
            expect(vertices).toEqual(Array.from(new Set(soup)));
            for (let i = 0; i < soup.length; ++i) {
                expect(vertices[indices[i]]).toBe(soup[i]);
            }
        });
    });

    it('the preconditions are enforced', () => {
        const uvs = new UniqueVerticesSimplices<P3>(3);
        const pool: P3[] = [[0, 0, 0], [1, 0, 0]];
        expect(() => uvs.generateIndexedSimplices([])).toThrow();
        // Not a multiple of the dimension.
        expect(() => uvs.generateIndexedSimplices([[0, 0, 0], [1, 0, 0]]))
            .toThrow();
        expect(() => uvs.removeDuplicateVertices([], [0, 0, 0])).toThrow();
        expect(() => uvs.removeDuplicateVertices(pool, [])).toThrow();
        expect(() => uvs.removeDuplicateVertices(pool, [0, 1])).toThrow();
        expect(() => uvs.removeDuplicateVertices(pool, [0, 1, 2])).toThrow();
        expect(() => uvs.removeDuplicateVertices(pool, [0, 1, -1])).toThrow();
        expect(() => uvs.removeUnusedVertices(pool, [0, 1, 2])).toThrow();
        expect(() => uvs.removeUnusedVerticesGrouped(pool, [])).toThrow();
        expect(() => uvs.removeDuplicateVerticesGrouped(pool, [[0, 1]]))
            .toThrow();
        expect(() => new UniqueVerticesSimplices<P3>(1)).toThrow();
    });
});
