import { describe, it, expect } from 'vitest';
import { SegmentMesh, SegmentMeshTopology } from '../src/SegmentMesh';
import { Vector } from '../src/Vector';

// A vertex pool V[i] = (i, 2*i) so a vertex is identifiable by its contents.
function pool(count: number): Vector[] {
    const vertices: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        vertices.push(Vector.fromArray([i, 2 * i]));
    }
    return vertices;
}

// Vertex degree induced by the segment index pairs.
function degrees(numVertices: number,
    indices: readonly [number, number][]): number[] {
    const degree = new Array<number>(numVertices).fill(0);
    for (const pair of indices) {
        ++degree[pair[0]];
        ++degree[pair[1]];
    }
    return degree;
}

describe('SegmentMesh construction', () => {
    it('the default constructor has UNKNOWN topology and no data', () => {
        const mesh = new SegmentMesh();
        expect(mesh.getTopology()).toBe(SegmentMeshTopology.UNKNOWN);
        expect(mesh.getVertices()).toEqual([]);
        expect(mesh.getIndices()).toEqual([]);
        expect(mesh.dimension).toBe(0);
    });

    it('exposes the topology enum through the class, as upstream does', () => {
        expect(SegmentMesh.Topology).toBe(SegmentMeshTopology);
        expect(SegmentMesh.Topology.INDEXED)
            .toBe(SegmentMeshTopology.INDEXED);
    });

    it('copies the vertex pool (C++ value semantics)', () => {
        const vertices = pool(4);
        const mesh = SegmentMesh.fromDisjoint(vertices);
        vertices[0].set(0, 99);
        expect(mesh.getVertices()[0].values).toEqual([0, 0]);
        expect(mesh.dimension).toBe(2);
    });
});

describe('SegmentMesh DISJOINT topology', () => {
    it('pairs consecutive vertices: S[i] = {2*i, 2*i+1}', () => {
        const mesh = SegmentMesh.fromDisjoint(pool(6));
        expect(mesh.getTopology()).toBe(SegmentMeshTopology.DISJOINT);
        expect(mesh.getIndices()).toEqual([[0, 1], [2, 3], [4, 5]]);
    });

    it('drops the trailing vertex when the count is odd (L/2 segments)', () => {
        const mesh = SegmentMesh.fromDisjoint(pool(5));
        expect(mesh.getIndices()).toEqual([[0, 1], [2, 3]]);
        expect(mesh.getVertices().length).toBe(5);
    });

    it('gives every vertex of an even pool degree one', () => {
        const mesh = SegmentMesh.fromDisjoint(pool(6));
        expect(degrees(6, mesh.getIndices())).toEqual([1, 1, 1, 1, 1, 1]);
    });

    it('requires at least two vertices', () => {
        expect(() => SegmentMesh.fromDisjoint(pool(1)))
            .toThrow('Invalid number of vertices.');
    });
});

describe('SegmentMesh CONTIGUOUS topology', () => {
    it('open polyline has L-1 segments <V[i],V[i+1]>', () => {
        const mesh = SegmentMesh.fromContiguous(pool(4), true);
        expect(mesh.getTopology()).toBe(SegmentMeshTopology.CONTIGUOUS_OPEN);
        expect(mesh.getIndices()).toEqual([[0, 1], [1, 2], [2, 3]]);
    });

    it('open polyline of two vertices is a single segment', () => {
        const mesh = SegmentMesh.fromContiguous(pool(2), true);
        expect(mesh.getIndices()).toEqual([[0, 1]]);
    });

    it('closed polyline has L segments, the wrap-around stored first', () => {
        const mesh = SegmentMesh.fromContiguous(pool(4), false);
        expect(mesh.getTopology()).toBe(SegmentMeshTopology.CONTIGUOUS_CLOSED);
        // Upstream fills mIndices[i1] = {i0, i1} with i0 the predecessor of
        // i1, so the wrap-around segment <V[L-1],V[0]> is stored first.
        expect(mesh.getIndices()).toEqual([[3, 0], [0, 1], [1, 2], [2, 3]]);
    });

    it('gives every closed-polyline vertex degree two', () => {
        const mesh = SegmentMesh.fromContiguous(pool(5), false);
        expect(degrees(5, mesh.getIndices())).toEqual([2, 2, 2, 2, 2]);
    });

    it('gives open-polyline interior vertices degree two, ends degree one', () => {
        const mesh = SegmentMesh.fromContiguous(pool(5), true);
        expect(degrees(5, mesh.getIndices())).toEqual([1, 2, 2, 2, 1]);
    });

    it('requires at least two vertices', () => {
        expect(() => SegmentMesh.fromContiguous(pool(1), true))
            .toThrow('Invalid number of vertices.');
        expect(() => SegmentMesh.fromContiguous([], false))
            .toThrow('Invalid number of vertices.');
    });
});

describe('SegmentMesh INDEXED topology', () => {
    it('stores a copy of the caller index pairs', () => {
        const indices: [number, number][] = [[0, 2], [2, 1]];
        const mesh = SegmentMesh.fromIndexed(pool(4), indices, true);
        expect(mesh.getTopology()).toBe(SegmentMeshTopology.INDEXED);
        indices[0][0] = 3;
        expect(mesh.getIndices()).toEqual([[0, 2], [2, 1]]);
    });

    it('does not require every vertex to participate', () => {
        const mesh = SegmentMesh.fromIndexed(pool(5), [[0, 1]], true);
        expect(mesh.getVertices().length).toBe(5);
        expect(mesh.getIndices().length).toBe(1);
    });

    it('allows degenerate segments with equal endpoints', () => {
        const mesh = SegmentMesh.fromIndexed(pool(3), [[1, 1]], true);
        expect(mesh.getIndices()).toEqual([[1, 1]]);
    });

    it('validates indices when requested', () => {
        expect(() => SegmentMesh.fromIndexed(pool(3), [[0, 3]], true))
            .toThrow('Invalid index into vertex array.');
        expect(() => SegmentMesh.fromIndexed(pool(3), [[5, 0]], true))
            .toThrow('Invalid index into vertex array.');
    });

    it('skips validation when validateIndices is false', () => {
        const mesh = SegmentMesh.fromIndexed(pool(3), [[0, 7]], false);
        expect(mesh.getIndices()).toEqual([[0, 7]]);
    });

    it('requires at least two vertices and one index pair', () => {
        expect(() => SegmentMesh.fromIndexed(pool(1), [[0, 0]], true))
            .toThrow('Invalid number of vertices or indices.');
        expect(() => SegmentMesh.fromIndexed(pool(3), [], true))
            .toThrow('Invalid number of vertices or indices.');
    });

    it('reproduces the other topologies, as the class documentation says', () => {
        const vertices = pool(6);
        const numVertices = vertices.length;

        const disjointAsIndexed: [number, number][] = [];
        for (let i = 0; i <= (numVertices - 2) / 2; ++i) {
            disjointAsIndexed.push([2 * i, 2 * i + 1]);
        }
        expect(SegmentMesh.fromDisjoint(vertices).getIndices())
            .toEqual(disjointAsIndexed);

        const openAsIndexed: [number, number][] = [];
        for (let i = 0; i <= numVertices - 2; ++i) {
            openAsIndexed.push([i, i + 1]);
        }
        expect(SegmentMesh.fromContiguous(vertices, true).getIndices())
            .toEqual(openAsIndexed);

        // The closed topology has the same segment set as the documented
        // S[i] = {i, (i+1) % L}, but upstream stores it rotated by one.
        const closedAsIndexed: [number, number][] = [];
        for (let i = 0; i <= numVertices - 1; ++i) {
            closedAsIndexed.push([i, (i + 1) % numVertices]);
        }
        const key = (pairs: readonly [number, number][]) =>
            pairs.map(p => p[0] + ',' + p[1]).sort();
        expect(key(SegmentMesh.fromContiguous(vertices, false).getIndices()))
            .toEqual(key(closedAsIndexed));
    });
});

describe('SegmentMesh vertex dimension', () => {
    it('supports 3D vertex pools', () => {
        const vertices = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([1, 1, 1])
        ];
        const mesh = SegmentMesh.fromContiguous(vertices, false);
        expect(mesh.dimension).toBe(3);
        expect(mesh.getIndices()).toEqual([[2, 0], [0, 1], [1, 2]]);
    });

    it('rejects a pool with mixed vertex dimensions', () => {
        const vertices = [
            Vector.fromArray([0, 0]),
            Vector.fromArray([1, 1, 1])
        ];
        expect(() => SegmentMesh.fromContiguous(vertices, true))
            .toThrow('SegmentMesh: mismatched vertex dimensions.');
    });
});
