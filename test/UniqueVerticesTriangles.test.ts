import { describe, it, expect } from 'vitest';
import { check, fc } from './helpers/arbitraries.js';
import { UniqueVerticesTriangles, defaultVertexKey } from '../src/UniqueVerticesTriangles.js';

type P3 = [number, number, number];

// Two triangles sharing the edge from (1,0,0) to (0,1,0), given as six
// vertices with duplicates.
const sharedEdgeSoup: P3[] = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0],
    [1, 0, 0], [1, 1, 0], [0, 1, 0]
];

describe('defaultVertexKey', () => {
    it('joins arrays and distinguishes different tuples', () => {
        expect(defaultVertexKey([1, 2, 3])).toBe('1,2,3');
        expect(defaultVertexKey([1, 2, 3])).toBe(defaultVertexKey([1, 2, 3]));
        expect(defaultVertexKey([1, 2, 3])).not.toBe(defaultVertexKey([1, 3, 2]));
    });

    it('handles primitives and objects', () => {
        expect(defaultVertexKey(7)).toBe('7');
        expect(defaultVertexKey('a')).toBe('a');
        expect(defaultVertexKey({ x: 1 })).toBe('{"x":1}');
    });
});

describe('UniqueVerticesTriangles.generateIndexedTriangles', () => {
    it('builds a vertex pool in first-occurrence order', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const result = uvt.generateIndexedTriangles(sharedEdgeSoup);
        expect(result.vertices).toEqual([
            [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]
        ]);
        expect(result.indices).toEqual([0, 1, 2, 1, 3, 2]);
    });

    it('satisfies the documented postconditions', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const result = uvt.generateIndexedTriangles(sharedEdgeSoup);
        expect(result.indices.length).toBe(sharedEdgeSoup.length);
        for (const index of result.indices) {
            expect(index).toBeGreaterThanOrEqual(0);
            expect(index).toBeLessThan(result.vertices.length);
        }
        // The pool is duplicate free.
        const keys = new Set(result.vertices.map(defaultVertexKey));
        expect(keys.size).toBe(result.vertices.length);
        // Every input vertex is recovered by its index.
        for (let i = 0; i < sharedEdgeSoup.length; ++i) {
            expect(result.vertices[result.indices[i]]).toEqual(sharedEdgeSoup[i]);
        }
    });

    it('keeps every vertex when there are no duplicates', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const soup: P3[] = [[0, 0, 0], [1, 0, 0], [0, 1, 0]];
        const result = uvt.generateIndexedTriangles(soup);
        expect(result.vertices).toEqual(soup);
        expect(result.indices).toEqual([0, 1, 2]);
    });

    it('produces triples with the Triples variant', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const result = uvt.generateIndexedTrianglesTriples(sharedEdgeSoup);
        expect(result.vertices.length).toBe(4);
        expect(result.triangles).toEqual([[0, 1, 2], [1, 3, 2]]);
    });

    it('validates the input count when validation is enabled', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        uvt.validate = true;
        expect(() => uvt.generateIndexedTriangles([])).toThrow('Invalid number of vertices.');
        expect(() => uvt.generateIndexedTriangles([[0, 0, 0], [1, 0, 0]]))
            .toThrow('Invalid number of vertices.');
        // Validation is off by default, matching upstream.
        expect(() => new UniqueVerticesTriangles<P3>().generateIndexedTriangles([]))
            .not.toThrow();
    });

    it('honors a custom key function', () => {
        // Collapse vertices that agree after rounding to one decimal place.
        const uvt = new UniqueVerticesTriangles<P3>(
            v => v.map(c => Math.round(c * 10) / 10).join(','));
        const soup: P3[] = [
            [0, 0, 0], [1, 0, 0], [0, 1, 0],
            [0.001, 0, 0], [1.002, 0, 0], [0, 1.004, 0]
        ];
        const result = uvt.generateIndexedTriangles(soup);
        expect(result.vertices.length).toBe(3);
        expect(result.indices).toEqual([0, 1, 2, 0, 1, 2]);
        // The retained representative is the first occurrence.
        expect(result.vertices[0]).toEqual([0, 0, 0]);
    });
});

describe('UniqueVerticesTriangles.removeDuplicateVertices', () => {
    it('collapses a pool with duplicates and remaps the indices', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [
            [0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]
        ];
        const indices = [0, 1, 2, 3, 4, 2];
        const result = uvt.removeDuplicateVertices(pool, indices);
        expect(result.vertices).toEqual([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]);
        // Index 3 (the duplicate of [1,0,0]) becomes 1.
        expect(result.indices).toEqual([0, 1, 2, 1, 3, 2]);
    });

    it('preserves the geometry referenced by each index', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [
            [2, 0, 0], [0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 0], [1, 1, 1]
        ];
        const indices = [0, 1, 3, 2, 4, 5];
        const result = uvt.removeDuplicateVertices(pool, indices);
        for (let i = 0; i < indices.length; ++i) {
            expect(result.vertices[result.indices[i]]).toEqual(pool[indices[i]]);
        }
    });

    it('works on triples', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0]];
        const result = uvt.removeDuplicateVerticesTriples(pool, [[0, 1, 2], [3, 4, 2]]);
        expect(result.vertices.length).toBe(4);
        expect(result.triangles).toEqual([[0, 1, 2], [1, 3, 2]]);
    });

    it('validates preconditions when enabled', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        uvt.validate = true;
        expect(() => uvt.removeDuplicateVertices([], [0, 0, 0]))
            .toThrow('Invalid number of vertices.');
        expect(() => uvt.removeDuplicateVertices([[0, 0, 0]], [0, 0]))
            .toThrow('Invalid number of indices.');
        expect(() => uvt.removeDuplicateVertices([[0, 0, 0]], [0, 0, 5]))
            .toThrow('Invalid index.');
        expect(() => uvt.removeDuplicateVerticesTriples([[0, 0, 0]], []))
            .toThrow('Invalid number of triangles.');
    });
});

describe('UniqueVerticesTriangles.removeUnusedVertices', () => {
    it('drops the unreferenced vertices and renumbers in increasing order', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [
            [0, 0, 0], [9, 9, 9], [1, 0, 0], [8, 8, 8], [0, 1, 0]
        ];
        // Vertices 1 and 3 are unused.
        const result = uvt.removeUnusedVertices(pool, [0, 2, 4]);
        expect(result.vertices).toEqual([[0, 0, 0], [1, 0, 0], [0, 1, 0]]);
        expect(result.indices).toEqual([0, 1, 2]);
    });

    it('visits the used indices in increasing order, as std::set does', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [
            [0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]
        ];
        // The first reference is to index 4, but the output pool is ordered
        // by the old indices, not by first use.
        const result = uvt.removeUnusedVertices(pool, [4, 2, 0]);
        expect(result.vertices).toEqual([[0, 0, 0], [2, 0, 0], [4, 0, 0]]);
        expect(result.indices).toEqual([2, 1, 0]);
    });

    it('keeps duplicate vertices that are separately referenced', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [[0, 0, 0], [0, 0, 0], [1, 0, 0]];
        const result = uvt.removeUnusedVertices(pool, [0, 1, 2]);
        expect(result.vertices.length).toBe(3);
    });

    it('preserves the geometry referenced by each index', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [
            [0, 0, 0], [9, 9, 9], [1, 0, 0], [8, 8, 8], [0, 1, 0], [7, 7, 7]
        ];
        const indices = [4, 2, 0, 0, 2, 4];
        const result = uvt.removeUnusedVertices(pool, indices);
        for (let i = 0; i < indices.length; ++i) {
            expect(result.vertices[result.indices[i]]).toEqual(pool[indices[i]]);
        }
        // Every output vertex occurs at least once in the output indices.
        for (let j = 0; j < result.vertices.length; ++j) {
            expect(result.indices).toContain(j);
        }
    });

    it('works on triples', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [
            [0, 0, 0], [9, 9, 9], [1, 0, 0], [8, 8, 8], [0, 1, 0]
        ];
        const result = uvt.removeUnusedVerticesTriples(pool, [[0, 2, 4]]);
        expect(result.vertices.length).toBe(3);
        expect(result.triangles).toEqual([[0, 1, 2]]);
    });
});

describe('UniqueVerticesTriangles.removeDuplicateAndUnusedVertices', () => {
    it('applies both operations', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [
            [0, 0, 0],   // used
            [5, 5, 5],   // unused
            [1, 0, 0],   // used
            [0, 0, 0],   // duplicate of index 0, used
            [0, 1, 0],   // used
            [5, 5, 5]    // duplicate of the unused vertex
        ];
        const indices = [0, 2, 4, 3, 2, 4];
        const result = uvt.removeDuplicateAndUnusedVertices(pool, indices);
        expect(result.vertices).toEqual([[0, 0, 0], [1, 0, 0], [0, 1, 0]]);
        expect(result.indices).toEqual([0, 1, 2, 0, 1, 2]);
        for (let i = 0; i < indices.length; ++i) {
            expect(result.vertices[result.indices[i]]).toEqual(pool[indices[i]]);
        }
    });

    it('works on triples', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [
            [0, 0, 0], [5, 5, 5], [1, 0, 0], [0, 0, 0], [0, 1, 0], [5, 5, 5]
        ];
        const result = uvt.removeDuplicateAndUnusedVerticesTriples(pool,
            [[0, 2, 4], [3, 2, 4]]);
        expect(result.vertices).toEqual([[0, 0, 0], [1, 0, 0], [0, 1, 0]]);
        expect(result.triangles).toEqual([[0, 1, 2], [0, 1, 2]]);
    });

    it('is idempotent on an already clean mesh', () => {
        const uvt = new UniqueVerticesTriangles<P3>();
        const pool: P3[] = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]];
        const indices = [0, 1, 2, 1, 3, 2];
        const once = uvt.removeDuplicateAndUnusedVertices(pool, indices);
        expect(once.vertices).toEqual(pool);
        expect(once.indices).toEqual(indices);
        const twice = uvt.removeDuplicateAndUnusedVertices(once.vertices, once.indices);
        expect(twice.vertices).toEqual(once.vertices);
        expect(twice.indices).toEqual(once.indices);
    });
});

describe('UniqueVerticesTriangles randomized round trip', () => {
    it('preserves the triangle geometry through generate and cleanup', () => {
        let state = 24680;
        const next = (m: number) => {
            state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
            return state % m;
        };

        // Build a triangle soup drawn from a small lattice, so duplicates are
        // common.
        const soup: P3[] = [];
        for (let i = 0; i < 300; ++i) {
            soup.push([next(4), next(4), next(4)]);
        }

        const uvt = new UniqueVerticesTriangles<P3>();
        const generated = uvt.generateIndexedTriangles(soup);

        // The pool holds only distinct lattice points.
        expect(generated.vertices.length).toBeLessThanOrEqual(64);
        const keys = new Set(generated.vertices.map(defaultVertexKey));
        expect(keys.size).toBe(generated.vertices.length);

        // The indexed representation reproduces the soup exactly.
        for (let i = 0; i < soup.length; ++i) {
            expect(generated.vertices[generated.indices[i]]).toEqual(soup[i]);
        }

        // The full cleanup does not change an already unique, fully used
        // pool other than possibly renumbering it.
        const cleaned = uvt.removeDuplicateAndUnusedVertices(
            generated.vertices, generated.indices);
        expect(cleaned.vertices.length).toBe(generated.vertices.length);
        for (let i = 0; i < soup.length; ++i) {
            expect(cleaned.vertices[cleaned.indices[i]]).toEqual(soup[i]);
        }
    });
});

describe('UniqueVerticesTriangles verification', () => {
    type V = [number, number, number];

    // Small integer coordinates so duplicates actually occur.
    const vertexArb: fc.Arbitrary<V> = fc.tuple(
        fc.integer({ min: 0, max: 3 }),
        fc.integer({ min: 0, max: 3 }),
        fc.integer({ min: 0, max: 3 })) as fc.Arbitrary<V>;

    /** A soup of triangles: a positive multiple of three vertices. */
    const soupArb = fc.array(vertexArb, { minLength: 3, maxLength: 30 })
        .map(vs => vs.slice(0, 3 * Math.floor(vs.length / 3)))
        .filter(vs => vs.length >= 3);

    const key = (v: V) => defaultVertexKey(v);

    /** An indexed mesh: a vertex pool plus a positive multiple of three indices. */
    const meshArb = fc.array(vertexArb, { minLength: 1, maxLength: 12 })
        .chain(vertices => fc.array(fc.integer({ min: 0, max: vertices.length - 1 }),
            { minLength: 3, maxLength: 30 })
            .map(indices => ({
                vertices,
                indices: indices.slice(0, 3 * Math.floor(indices.length / 3))
            }))
            .filter(m => m.indices.length >= 3));

    it('generateIndexedTriangles meets every documented postcondition', () => {
        check(soupArb, inVertices => {
            const uvt = new UniqueVerticesTriangles<V>();
            const { vertices, indices } = uvt.generateIndexedTriangles(inVertices);
            // 2. indices.length = inVertices.length
            expect(indices.length).toBe(inVertices.length);
            // 3. 0 <= indices[i] < vertices.length
            expect(indices.every(i => 0 <= i && i < vertices.length)).toBe(true);
            // 1. vertices has unique vertices
            expect(new Set(vertices.map(key)).size).toBe(vertices.length);
            // The pool reproduces the input soup exactly.
            expect(indices.map(i => key(vertices[i]!))).toEqual(inVertices.map(key));
            return true;
        });
    });

    it('the unique vertices are numbered in order of first occurrence', () => {
        check(soupArb, inVertices => {
            const uvt = new UniqueVerticesTriangles<V>();
            const { vertices } = uvt.generateIndexedTriangles(inVertices);
            // An independent computation of the same pool.
            const seen = new Set<string | number>();
            const expected: V[] = [];
            for (const v of inVertices) {
                const k = key(v);
                if (!seen.has(k)) { seen.add(k); expected.push(v); }
            }
            return vertices.length === expected.length
                && vertices.every((v, i) => key(v) === key(expected[i]!));
        });
    });

    it('removeDuplicateVertices preserves the geometry each index refers to', () => {
        check(meshArb, ({ vertices: inVertices, indices: inIndices }) => {
            const uvt = new UniqueVerticesTriangles<V>();
            const out = uvt.removeDuplicateVertices(inVertices, inIndices);
            expect(out.indices.length).toBe(inIndices.length);
            expect(new Set(out.vertices.map(key)).size).toBe(out.vertices.length);
            expect(out.indices.every(i => 0 <= i && i < out.vertices.length)).toBe(true);
            for (let i = 0; i < inIndices.length; ++i) {
                expect(key(out.vertices[out.indices[i]!]!))
                    .toBe(key(inVertices[inIndices[i]!]!));
            }
            return true;
        });
    });

    it('removeUnusedVertices keeps exactly the referenced vertices, in old-index order', () => {
        check(meshArb, ({ vertices: inVertices, indices: inIndices }) => {
            const uvt = new UniqueVerticesTriangles<V>();
            const out = uvt.removeUnusedVertices(inVertices, inIndices);
            // Independent computation: the used old indices sorted ascending,
            // matching upstream's std::set iteration order.
            const used = [...new Set(inIndices)].sort((a, b) => a - b);
            expect(out.vertices.length).toBe(used.length);
            expect(out.vertices.map(key)).toEqual(used.map(i => key(inVertices[i]!)));
            expect(out.indices.length).toBe(inIndices.length);
            for (let i = 0; i < inIndices.length; ++i) {
                expect(out.indices[i]).toBe(used.indexOf(inIndices[i]!));
                expect(key(out.vertices[out.indices[i]!]!))
                    .toBe(key(inVertices[inIndices[i]!]!));
            }
            // Postcondition 4: every output vertex is referenced.
            expect(new Set(out.indices).size).toBe(out.vertices.length);
            return true;
        });
    });

    it('removeUnusedVertices keeps duplicates that are separately referenced', () => {
        check(meshArb, ({ vertices: inVertices, indices: inIndices }) => {
            const uvt = new UniqueVerticesTriangles<V>();
            const out = uvt.removeUnusedVertices(inVertices, inIndices);
            // It removes unused slots only; it never merges equal vertices.
            return out.vertices.length === new Set(inIndices).size;
        });
    });

    it('removeDuplicateAndUnusedVertices is the composition of the two steps', () => {
        check(meshArb, ({ vertices: inVertices, indices: inIndices }) => {
            const uvt = new UniqueVerticesTriangles<V>();
            const both = uvt.removeDuplicateAndUnusedVertices(inVertices, inIndices);
            const step1 = uvt.removeDuplicateVertices(inVertices, inIndices);
            const step2 = uvt.removeUnusedVertices(step1.vertices, step1.indices);
            expect(both.indices).toEqual(step2.indices);
            expect(both.vertices.map(key)).toEqual(step2.vertices.map(key));
            // The result is both duplicate-free and fully used, and it still
            // describes the same geometry.
            expect(new Set(both.vertices.map(key)).size).toBe(both.vertices.length);
            expect(new Set(both.indices).size).toBe(both.vertices.length);
            for (let i = 0; i < inIndices.length; ++i) {
                expect(key(both.vertices[both.indices[i]!]!))
                    .toBe(key(inVertices[inIndices[i]!]!));
            }
            return true;
        });
    });

    it('removeDuplicateAndUnusedVertices is idempotent', () => {
        check(meshArb, ({ vertices: inVertices, indices: inIndices }) => {
            const uvt = new UniqueVerticesTriangles<V>();
            const once = uvt.removeDuplicateAndUnusedVertices(inVertices, inIndices);
            const twice = uvt.removeDuplicateAndUnusedVertices(once.vertices, once.indices);
            return JSON.stringify(twice.indices) === JSON.stringify(once.indices)
                && JSON.stringify(twice.vertices) === JSON.stringify(once.vertices);
        });
    });

    it('the Triples variants are exactly the flat variants reshaped', () => {
        check(meshArb, ({ vertices: inVertices, indices: inIndices }) => {
            const uvt = new UniqueVerticesTriangles<V>();
            const triples: number[][] = [];
            for (let t = 0; t < inIndices.length / 3; ++t) {
                triples.push([inIndices[3 * t]!, inIndices[3 * t + 1]!, inIndices[3 * t + 2]!]);
            }
            const reshape = (flat: number[]) => {
                const out: number[][] = [];
                for (let t = 0; t < flat.length / 3; ++t) {
                    out.push([flat[3 * t]!, flat[3 * t + 1]!, flat[3 * t + 2]!]);
                }
                return out;
            };
            const dupFlat = uvt.removeDuplicateVertices(inVertices, inIndices);
            const dupTri = uvt.removeDuplicateVerticesTriples(inVertices, triples);
            expect(dupTri.triangles).toEqual(reshape(dupFlat.indices));
            expect(dupTri.vertices.map(key)).toEqual(dupFlat.vertices.map(key));

            const unusedFlat = uvt.removeUnusedVertices(inVertices, inIndices);
            const unusedTri = uvt.removeUnusedVerticesTriples(inVertices, triples);
            expect(unusedTri.triangles).toEqual(reshape(unusedFlat.indices));
            expect(unusedTri.vertices.map(key)).toEqual(unusedFlat.vertices.map(key));

            const bothFlat = uvt.removeDuplicateAndUnusedVertices(inVertices, inIndices);
            const bothTri = uvt.removeDuplicateAndUnusedVerticesTriples(inVertices, triples);
            expect(bothTri.triangles).toEqual(reshape(bothFlat.indices));
            expect(bothTri.vertices.map(key)).toEqual(bothFlat.vertices.map(key));
            return true;
        });
    });

    it('generateIndexedTrianglesTriples round trips the triangle soup', () => {
        check(soupArb, inVertices => {
            const uvt = new UniqueVerticesTriangles<V>();
            const { vertices, triangles } = uvt.generateIndexedTrianglesTriples(inVertices);
            expect(triangles.length).toBe(inVertices.length / 3);
            for (let t = 0; t < triangles.length; ++t) {
                for (let j = 0; j < 3; ++j) {
                    expect(key(vertices[triangles[t]![j]!]!))
                        .toBe(key(inVertices[3 * t + j]!));
                }
            }
            return true;
        });
    });

    it('validation rejects malformed inputs only when enabled', () => {
        check(fc.array(vertexArb, { minLength: 0, maxLength: 8 })
            .filter(vs => vs.length === 0 || vs.length % 3 !== 0), inVertices => {
                const quiet = new UniqueVerticesTriangles<V>();
                expect(() => quiet.generateIndexedTriangles(inVertices)).not.toThrow();
                const loud = new UniqueVerticesTriangles<V>();
                loud.validate = true;
                expect(() => loud.generateIndexedTriangles(inVertices)).toThrow();
                return true;
            });
    });

    it('a custom key function drives the equivalence used for duplicates', () => {
        // Collapse vertices by their first coordinate only.
        check(soupArb, inVertices => {
            const uvt = new UniqueVerticesTriangles<V>(v => v[0]);
            const { vertices, indices } = uvt.generateIndexedTriangles(inVertices);
            const distinctFirst = new Set(inVertices.map(v => v[0]));
            return vertices.length === distinctFirst.size
                && indices.every((idx, i) => vertices[idx]![0] === inVertices[i]![0]);
        });
    });
});
