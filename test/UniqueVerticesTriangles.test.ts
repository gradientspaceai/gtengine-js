import { describe, it, expect } from 'vitest';
import { UniqueVerticesTriangles, defaultVertexKey } from '../src/UniqueVerticesTriangles';

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
