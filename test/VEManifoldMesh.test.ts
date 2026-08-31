import { describe, it, expect } from 'vitest';
import {
    VEManifoldMesh,
    VEManifoldMeshVertex,
    VEManifoldMeshEdge
} from '../src/VEManifoldMesh';

// Build a closed polyline with the directed edges <0,1>, <1,2>, <2,0>.
function makeTriangleLoop(): VEManifoldMesh {
    const mesh = new VEManifoldMesh();
    expect(mesh.insert(0, 1)).not.toBeNull();
    expect(mesh.insert(1, 2)).not.toBeNull();
    expect(mesh.insert(2, 0)).not.toBeNull();
    return mesh;
}

describe('VEManifoldMesh insertion', () => {
    it('creates the edge and its two vertices', () => {
        const mesh = new VEManifoldMesh();
        const edge = mesh.insert(3, 7);
        expect(edge).not.toBeNull();
        expect(edge!.V).toEqual([3, 7]);
        expect(edge!.E).toEqual([null, null]);

        expect(mesh.getNumVertices()).toBe(2);
        expect(mesh.getNumEdges()).toBe(1);
        expect(mesh.getEdge(3, 7)).toBe(edge);

        // The edge key is the directed pair, so <7,3> is a different edge.
        expect(mesh.getEdge(7, 3)).toBeNull();

        const vertex3 = mesh.getVertex(3);
        const vertex7 = mesh.getVertex(7);
        expect(vertex3).not.toBeNull();
        expect(vertex7).not.toBeNull();
        // A one-edge vertex always has its edge at index zero.
        expect(vertex3!.E).toEqual([edge, null]);
        expect(vertex7!.E).toEqual([edge, null]);
        expect(mesh.isClosed()).toBe(false);
    });

    it('returns null when the directed edge already exists', () => {
        const mesh = new VEManifoldMesh();
        const edge = mesh.insert(0, 1);
        expect(mesh.insert(0, 1)).toBeNull();
        expect(mesh.getNumEdges()).toBe(1);
        expect(mesh.getEdge(0, 1)).toBe(edge);

        // The reversed directed edge is a distinct key and inserts fine, but
        // it makes vertices 0 and 1 shared twice, so the polyline closes.
        expect(mesh.insert(1, 0)).not.toBeNull();
        expect(mesh.getNumEdges()).toBe(2);
        expect(mesh.isClosed()).toBe(true);
    });

    it('links the adjacent edges through the shared vertex', () => {
        const mesh = new VEManifoldMesh();
        const e01 = mesh.insert(0, 1)!;
        const e12 = mesh.insert(1, 2)!;

        // E[i] of an edge points to the edge sharing V[i].
        expect(e01.E).toEqual([null, e12]);
        expect(e12.E).toEqual([e01, null]);

        const vertex1 = mesh.getVertex(1)!;
        expect(vertex1.E).toEqual([e01, e12]);
        expect(mesh.isClosed()).toBe(false);
    });

    it('closes the polyline when every vertex is shared twice', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.getNumVertices()).toBe(3);
        expect(mesh.getNumEdges()).toBe(3);
        expect(mesh.isClosed()).toBe(true);

        const e01 = mesh.getEdge(0, 1)!;
        const e12 = mesh.getEdge(1, 2)!;
        const e20 = mesh.getEdge(2, 0)!;
        expect(e01.E).toEqual([e20, e12]);
        expect(e12.E).toEqual([e01, e20]);
        expect(e20.E).toEqual([e12, e01]);

        for (const vertex of mesh.getVertices()) {
            expect(vertex.E[0]).not.toBeNull();
            expect(vertex.E[1]).not.toBeNull();
        }
    });

    it('throws when the insertion would make the mesh nonmanifold', () => {
        const mesh = new VEManifoldMesh();
        mesh.insert(0, 1);
        mesh.insert(1, 2);
        // Vertex 1 already has two edges.
        expect(() => mesh.insert(1, 3)).toThrow('The mesh must be manifold.');
    });

    it('returns null instead of throwing when the exception is disabled', () => {
        const mesh = new VEManifoldMesh();
        mesh.insert(0, 1);
        mesh.insert(1, 2);
        mesh.throwOnNonmanifoldInsertion(false);
        expect(mesh.insert(1, 3)).toBeNull();

        // Upstream quirk, preserved by the port: the rejected edge was added
        // to the edge map before the manifold test, so it stays there.
        expect(mesh.getEdge(1, 3)).not.toBeNull();
        // Vertex 3 was never reached, so the loop stopped at the failure.
        expect(mesh.getVertex(3)).toBeNull();
        // Vertex 1 is unchanged: it still refers to the first two edges.
        expect(mesh.getVertex(1)!.E).toEqual([
            mesh.getEdge(0, 1), mesh.getEdge(1, 2)
        ]);
    });
});

describe('VEManifoldMesh removal', () => {
    it('returns false for an edge not in the mesh', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.remove(0, 2)).toBe(false);
        expect(mesh.remove(5, 6)).toBe(false);
        expect(mesh.getNumEdges()).toBe(3);
    });

    it('opens a closed polyline and keeps the endpoints', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.remove(1, 2)).toBe(true);
        expect(mesh.getNumEdges()).toBe(2);
        // Vertices 1 and 2 are still used by the remaining edges.
        expect(mesh.getNumVertices()).toBe(3);
        expect(mesh.isClosed()).toBe(false);

        const e01 = mesh.getEdge(0, 1)!;
        const e20 = mesh.getEdge(2, 0)!;
        expect(mesh.getEdge(1, 2)).toBeNull();

        // One-edge vertices always have the edge at index zero.
        expect(mesh.getVertex(1)!.E).toEqual([e01, null]);
        expect(mesh.getVertex(2)!.E).toEqual([e20, null]);
        expect(mesh.getVertex(0)!.E).toEqual([e01, e20]);

        // The adjacency references to the removed edge were cleared.
        expect(e01.E).toEqual([e20, null]);
        expect(e20.E).toEqual([null, e01]);
    });

    it('removes vertices when their last edge is removed', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.remove(0, 1)).toBe(true);
        expect(mesh.remove(1, 2)).toBe(true);
        // Vertex 1 has no remaining edge.
        expect(mesh.getVertex(1)).toBeNull();
        expect(mesh.getNumVertices()).toBe(2);

        expect(mesh.remove(2, 0)).toBe(true);
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        // A mesh with no vertices is vacuously closed.
        expect(mesh.isClosed()).toBe(true);
    });

    it('restores the original state after insert then remove', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.remove(1, 2)).toBe(true);
        const e12 = mesh.insert(1, 2)!;
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.getNumEdges()).toBe(3);
        expect(mesh.getNumVertices()).toBe(3);

        const e01 = mesh.getEdge(0, 1)!;
        const e20 = mesh.getEdge(2, 0)!;
        expect(e12.E).toEqual([e01, e20]);
        expect(e01.E).toEqual([e20, e12]);
        expect(e20.E).toEqual([e12, e01]);
    });
});

describe('VEManifoldMesh iteration order', () => {
    it('returns vertices by increasing index and edges lexicographically', () => {
        const mesh = new VEManifoldMesh();
        // Insert out of order; the accessors must sort as std::map iterates.
        mesh.insert(5, 2);
        mesh.insert(2, 9);
        mesh.insert(9, 5);

        expect(mesh.getVertices().map(v => v.V)).toEqual([2, 5, 9]);
        expect(mesh.getEdges().map(e => e.V)).toEqual([[2, 9], [5, 2], [9, 5]]);
    });
});

describe('VEManifoldMesh creators', () => {
    it('uses the supplied vertex and edge creators', () => {
        class TaggedVertex extends VEManifoldMeshVertex {
            tag: string;
            constructor(v: number) {
                super(v);
                this.tag = 'v' + v;
            }
        }
        class TaggedEdge extends VEManifoldMeshEdge {
            tag: string;
            constructor(v0: number, v1: number) {
                super(v0, v1);
                this.tag = 'e' + v0 + '_' + v1;
            }
        }

        const mesh = new VEManifoldMesh(
            v => new TaggedVertex(v),
            (v0, v1) => new TaggedEdge(v0, v1));
        const edge = mesh.insert(4, 6)!;
        expect(edge).toBeInstanceOf(TaggedEdge);
        expect((edge as TaggedEdge).tag).toBe('e4_6');
        expect(mesh.getVertex(4)).toBeInstanceOf(TaggedVertex);
        expect((mesh.getVertex(6) as TaggedVertex).tag).toBe('v6');
    });
});
