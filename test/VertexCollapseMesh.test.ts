import { describe, expect, it } from 'vitest';
import {
    VERTEX_COLLAPSE_MESH_INVALID_VERTEX, VertexCollapseMesh,
    VertexCollapseMeshVertex
} from '../src/VertexCollapseMesh.js';
import type { VertexCollapseMeshResult } from '../src/VertexCollapseMesh.js';
import { VETManifoldMesh } from '../src/VETManifoldMesh.js';
import { Vector, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { check, fc, positive, wellScaled } from './helpers/arbitraries.js';

// An n-by-n grid of vertices in the plane z = 0, each quad split by the
// diagonal from its lower-left to its upper-right corner. There are
// (n-2)^2 interior vertices, each shared by six triangles.
function grid(n: number): { positions: Vector[]; indices: number[] } {
    const positions: Vector[] = [];
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            positions.push(Vector.fromArray([c, r, 0]));
        }
    }
    const indices: number[] = [];
    for (let r = 0; r + 1 < n; ++r) {
        for (let c = 0; c + 1 < n; ++c) {
            const a = r * n + c;
            const b = r * n + c + 1;
            const d = (r + 1) * n + c + 1;
            const e = (r + 1) * n + c;
            indices.push(a, b, d, a, d, e);
        }
    }
    return { positions, indices };
}

const octahedronPositions = [
    Vector.fromArray([0, 0, 1]),
    Vector.fromArray([1, 0, 0]),
    Vector.fromArray([0, 1, 0]),
    Vector.fromArray([-1, 0, 0]),
    Vector.fromArray([0, -1, 0]),
    Vector.fromArray([0, 0, -1])
];
const octahedronIndices = [
    0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1,
    5, 2, 1, 5, 3, 2, 5, 4, 3, 5, 1, 4
];

// Every edge is shared by one or two triangles, and every triangle is
// registered with each of its three edges. This is the structural invariant
// of an edge-triangle manifold mesh.
function expectManifold(mesh: VETManifoldMesh): void {
    for (const edge of mesh.getEdges()) {
        expect(edge.T[0]).not.toBeNull();
        if (edge.T[1] !== null) {
            expect(edge.T[1]).not.toBe(edge.T[0]);
        }
    }
    for (const tri of mesh.getTriangles()) {
        for (let j = 0; j < 3; ++j) {
            const edge = tri.E[j];
            expect(edge).not.toBeNull();
            expect(edge!.T[0] === tri || edge!.T[1] === tri).toBe(true);
        }
    }
    // Every vertex in the vertex map is shared by at least one triangle and
    // its adjacency sets agree with the triangle adjacency.
    for (const vertex of mesh.getVertices()) {
        expect(vertex.TAdjacent.size).toBeGreaterThan(0);
        expect(vertex.VAdjacent.size).toBeGreaterThan(0);
        for (const tri of vertex.getTAdjacent()) {
            expect(tri.V.includes(vertex.V)).toBe(true);
        }
    }
}

// Collapse until no more collapses are allowed, checking the per-step
// invariants. Returns the sequence of results.
function collapseAll(vcm: VertexCollapseMesh,
    maxSteps: number): VertexCollapseMeshResult[] {
    const mesh = vcm.getMesh();
    const results: VertexCollapseMeshResult[] = [];
    for (let step = 0; step < maxSteps; ++step) {
        const numVerticesBefore = mesh.getNumVertices();
        const numTrianglesBefore = mesh.getNumTriangles();
        const result = vcm.doCollapse();
        if (!result.collapsed) {
            expect(result.vertex).toBe(VERTEX_COLLAPSE_MESH_INVALID_VERTEX);
            expect(result.removed).toHaveLength(0);
            expect(result.inserted).toHaveLength(0);
            // Nothing changed on a failed collapse.
            expect(mesh.getNumVertices()).toBe(numVerticesBefore);
            expect(mesh.getNumTriangles()).toBe(numTrianglesBefore);
            return results;
        }

        // The collapsed vertex is gone; every other vertex survives.
        expect(mesh.getVertex(result.vertex)).toBeNull();
        expect(mesh.getNumVertices()).toBe(numVerticesBefore - 1);

        // The link of an interior vertex has as many triangles as it has
        // boundary vertices, and its triangulation has two fewer.
        expect(result.removed.length).toBeGreaterThanOrEqual(3);
        expect(result.inserted.length).toBe(result.removed.length - 2);
        expect(mesh.getNumTriangles()).toBe(numTrianglesBefore -
            result.removed.length + result.inserted.length);

        // The removed triangles are gone and each mentioned the vertex.
        for (const tri of result.removed) {
            expect(tri.V.includes(result.vertex)).toBe(true);
            expect(mesh.getTriangle(tri.V[0], tri.V[1], tri.V[2])).toBeNull();
        }
        // The inserted triangles are present and avoid the vertex.
        for (const tri of result.inserted) {
            expect(tri.V.includes(result.vertex)).toBe(false);
            expect(mesh.getTriangle(tri.V[0], tri.V[1], tri.V[2]))
                .not.toBeNull();
        }

        expectManifold(mesh);
        results.push(result);
    }
    throw new Error('The decimation did not terminate.');
}

describe('VertexCollapseMesh', () => {
    it('builds the mesh from the input triangles', () => {
        const g = grid(4);
        const vcm = new VertexCollapseMesh(g.positions, g.indices);
        const mesh = vcm.getMesh();
        expect(mesh.getNumVertices()).toBe(16);
        expect(mesh.getNumTriangles()).toBe(18);
        expectManifold(mesh);

        // The boundary vertices of the grid are tagged; the interior ones are
        // not.
        for (let r = 0; r < 4; ++r) {
            for (let c = 0; c < 4; ++c) {
                const vertex = mesh.getVertex(r * 4 + c) as
                    VertexCollapseMeshVertex;
                const onBoundary = (r === 0 || r === 3 || c === 0 || c === 3);
                expect(vertex.isBoundary).toBe(onBoundary);
            }
        }
    });

    it('removes every interior vertex of a planar grid', () => {
        const g = grid(4);
        const vcm = new VertexCollapseMesh(g.positions, g.indices);
        const mesh = vcm.getMesh();
        const results = collapseAll(vcm, 32);

        // The four interior vertices of the 4-by-4 grid are collapsed, each
        // removing six triangles and inserting four.
        expect(results).toHaveLength(4);
        const collapsedVertices = results.map(r => r.vertex).sort(
            (a, b) => a - b);
        expect(collapsedVertices).toEqual([5, 6, 9, 10]);

        // The result is a triangulation of the 12-vertex boundary polygon:
        // 12 vertices and 12 - 2 = 10 triangles.
        expect(mesh.getNumVertices()).toBe(12);
        expect(mesh.getNumTriangles()).toBe(10);
        for (const v of [0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15]) {
            expect(mesh.getVertex(v)).not.toBeNull();
        }
        expectManifold(mesh);
    });

    it('collapses the nearly planar vertices of a grid first', () => {
        // Vertex 6 lies far off the plane of the rest of the grid, so its
        // weight is large and it is collapsed last among the interior
        // vertices.
        const g = grid(4);
        g.positions[6] = Vector.fromArray([2, 1, 6]);
        const vcm = new VertexCollapseMesh(g.positions, g.indices);
        const results = collapseAll(vcm, 32);
        expect(results).toHaveLength(4);
        expect(results[results.length - 1].vertex).toBe(6);
    });

    it('decimates a nonplanar grid and terminates', () => {
        const g = grid(5);
        for (let i = 0; i < g.positions.length; ++i) {
            g.positions[i].values[2] = 0.3 * Math.sin(i * 1.7);
        }
        const vcm = new VertexCollapseMesh(g.positions, g.indices);
        const mesh = vcm.getMesh();
        expect(mesh.getNumVertices()).toBe(25);
        expect(mesh.getNumTriangles()).toBe(32);

        const results = collapseAll(vcm, 64);
        // All nine interior vertices are collapsed; the 16 boundary vertices
        // remain and the boundary polygon is triangulated.
        expect(results).toHaveLength(9);
        expect(mesh.getNumVertices()).toBe(16);
        expect(mesh.getNumTriangles()).toBe(14);
        for (const r of results) {
            const row = Math.floor(r.vertex / 5);
            const col = r.vertex % 5;
            expect(row).toBeGreaterThan(0);
            expect(row).toBeLessThan(4);
            expect(col).toBeGreaterThan(0);
            expect(col).toBeLessThan(4);
        }
    });

    it('decimates a closed mesh (octahedron)', () => {
        const vcm = new VertexCollapseMesh(octahedronPositions,
            octahedronIndices);
        const mesh = vcm.getMesh();
        expect(mesh.getNumVertices()).toBe(6);
        expect(mesh.getNumTriangles()).toBe(8);
        expect(mesh.isClosed()).toBe(true);
        // A closed mesh has no boundary vertices, so every vertex is a
        // collapse candidate.
        for (const vertex of mesh.getVertices()) {
            expect((vertex as VertexCollapseMeshVertex).isBoundary).toBe(false);
        }

        const results = collapseAll(vcm, 16);
        expect(results.length).toBeGreaterThanOrEqual(1);
        // The first collapse removes the four triangles of a valence-4 vertex
        // and inserts the two triangles of its square link.
        expect(results[0].removed).toHaveLength(4);
        expect(results[0].inserted).toHaveLength(2);
        expect(mesh.getNumVertices()).toBe(6 - results.length);
    });

    it('reports no collapse for a mesh with no interior vertices', () => {
        const positions = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([1, 1, 0]),
            Vector.fromArray([0, 1, 0])
        ];
        const indices = [0, 1, 2, 0, 2, 3];
        const vcm = new VertexCollapseMesh(positions, indices);
        expect(vcm.getMesh().getNumTriangles()).toBe(2);
        const result = vcm.doCollapse();
        expect(result.collapsed).toBe(false);
        expect(result.vertex).toBe(VERTEX_COLLAPSE_MESH_INVALID_VERTEX);
        expect(vcm.getMesh().getNumTriangles()).toBe(2);
    });

    it('handles degenerate construction input', () => {
        for (const vcm of [
            new VertexCollapseMesh(null, null),
            new VertexCollapseMesh([], []),
            new VertexCollapseMesh([Vector.fromArray([0, 0, 0])], [0, 0]),
            new VertexCollapseMesh(octahedronPositions, [0, 1])
        ]) {
            expect(vcm.getMesh().getNumTriangles()).toBe(0);
            const result = vcm.doCollapse();
            expect(result.collapsed).toBe(false);
            expect(result.vertex).toBe(VERTEX_COLLAPSE_MESH_INVALID_VERTEX);
            expect(result.removed).toHaveLength(0);
            expect(result.inserted).toHaveLength(0);
        }
    });

    it('computes the vertex weight and normal of a link', () => {
        // A vertex shared by six unit-area-pair triangles of the planar grid:
        // the normal is +z and the neighbors are in the tangent plane, so the
        // weight is exactly the sum of the |E0 x E1| values (twice the total
        // area of the six triangles).
        const g = grid(4);
        const vcm = new VertexCollapseMesh(g.positions, g.indices);
        const vertex = vcm.getMesh().getVertex(5) as VertexCollapseMeshVertex;
        const weight = vertex.computeWeight(g.positions);
        expect(weight).toBeCloseTo(6, 12);
        expect(vertex.normal.values[0]).toBeCloseTo(0, 12);
        expect(vertex.normal.values[1]).toBeCloseTo(0, 12);
        expect(vertex.normal.values[2]).toBeCloseTo(1, 12);

        // Raise the vertex out of the plane: the neighbors now project onto
        // the normal line with nonzero lengths, so the weight grows.
        const raised = g.positions.slice();
        raised[5] = Vector.fromArray([1, 1, 2]);
        expect(vertex.computeWeight(raised)).toBeGreaterThan(weight);
    });
});

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// A grid with a random height field. wellScaled keeps the heights away from
// subnormals; computeWeight normalizes a sum of cross products, which would
// underflow on 1e-160-scale coordinates.
const heightGrid = fc.tuple(fc.integer({ min: 4, max: 6 }),
    fc.array(wellScaled(-0.6, 0.6), { minLength: 36, maxLength: 36 }))
    .map(([n, heights]) => {
        const g = grid(n);
        for (let i = 0; i < g.positions.length; ++i) {
            g.positions[i].values[2] = heights[i % heights.length];
        }
        return { n, ...g };
    });

// The octahedron with the vertices pushed onto a random sphere, so the mesh
// stays closed and convex but is no longer symmetric.
const roundOctahedron = fc.tuple(positive(4, 0.5),
    fc.array(wellScaled(-0.15, 0.15), { minLength: 18, maxLength: 18 }))
    .map(([r, jitter]) => {
        const positions = octahedronPositions.map((p, i) => {
            const q = p.clone();
            for (let k = 0; k < 3; ++k) { q.values[k] += jitter[3 * i + k]; }
            normalize(q);
            return mul(q, r);
        });
        return { positions, indices: octahedronIndices.slice() };
    });

function euler(mesh: VETManifoldMesh): number {
    return mesh.getNumVertices() - mesh.getNumEdges() + mesh.getNumTriangles();
}

// An independent evaluation of VCVertex::ComputeWeight.
function referenceWeight(vertex: VertexCollapseMeshVertex,
    positions: readonly Vector[]): { weight: number, normal: Vector } {
    let weight = 0;
    const normal = Vector.fromArray([0, 0, 0]);
    for (const tri of vertex.getTAdjacent()) {
        const e0 = sub(positions[tri.V[1]], positions[tri.V[0]]);
        const e1 = sub(positions[tri.V[2]], positions[tri.V[0]]);
        const n = cross(e0, e1);
        for (let k = 0; k < 3; ++k) { normal.values[k] += n.values[k]; }
        weight += length(n);
    }
    normalize(normal);
    for (const index of vertex.getVAdjacent()) {
        weight += Math.abs(dot(normal, sub(positions[index], positions[vertex.V])));
    }
    return { weight, normal };
}

describe('VertexCollapseMesh verification', () => {
    it('preserves the Euler characteristic and the boundary of a grid', () => {
        check(heightGrid, ({ n, positions, indices }) => {
            const vcm = new VertexCollapseMesh(positions, indices);
            const mesh = vcm.getMesh();
            const chi = euler(mesh);

            // The boundary vertices of the grid must never be collapsed.
            const boundary = new Set<number>();
            for (let r = 0; r < n; ++r) {
                for (let c = 0; c < n; ++c) {
                    if (r === 0 || c === 0 || r === n - 1 || c === n - 1) {
                        boundary.add(r * n + c);
                    }
                }
            }

            const results = collapseAll(vcm, positions.length + 5);
            for (const result of results) {
                expect(boundary.has(result.vertex)).toBe(false);
            }
            expect(new Set(results.map(r => r.vertex)).size)
                .toBe(results.length);
            expect(euler(mesh)).toBe(chi);
            // Only the boundary ring can survive a complete decimation; when
            // the decimation stops early, whatever is left is still a
            // superset of the boundary.
            for (const v of boundary) {
                expect(mesh.getVertex(v)).not.toBeNull();
            }
        });
    }, 30000);

    it('decimates a closed mesh down to the fewest triangles', () => {
        check(roundOctahedron, ({ positions, indices }) => {
            const vcm = new VertexCollapseMesh(positions, indices);
            const mesh = vcm.getMesh();
            expect(euler(mesh)).toBe(2);
            const results = collapseAll(vcm, positions.length + 5);
            // Every vertex of a closed mesh is interior, so the decimation
            // proceeds until a collapse would make the mesh nonmanifold.
            expect(results.length).toBeGreaterThan(0);
            expect(mesh.isClosed()).toBe(true);
            expect(euler(mesh)).toBe(2);
            expect(mesh.getNumVertices()).toBe(positions.length - results.length);
        });
    }, 30000);

    it('never leaves a collapsed vertex or a removed triangle behind', () => {
        check(heightGrid, ({ positions, indices }) => {
            const vcm = new VertexCollapseMesh(positions, indices);
            const mesh = vcm.getMesh();
            const results = collapseAll(vcm, positions.length + 5);
            for (const result of results) {
                expect(mesh.getVertex(result.vertex)).toBeNull();
                for (const tri of result.removed) {
                    expect(tri.V).toContain(result.vertex);
                }
                for (const tri of result.inserted) {
                    expect(tri.V).not.toContain(result.vertex);
                }
            }
            // No triangle mentions a collapsed vertex any more.
            const gone = new Set(results.map(r => r.vertex));
            for (const tri of mesh.getTriangles()) {
                for (const v of tri.V) { expect(gone.has(v)).toBe(false); }
            }
        });
    }, 30000);

    it('leaves the caller position array untouched', () => {
        check(heightGrid, ({ positions, indices }) => {
            const before = positions.map(p => p.values.slice());
            const vcm = new VertexCollapseMesh(positions, indices);
            collapseAll(vcm, positions.length + 5);
            for (let i = 0; i < positions.length; ++i) {
                expect(positions[i].values).toEqual(before[i]);
            }
        });
    }, 30000);

    it('computes the weight and normal exactly as the upstream formula does', () => {
        check(heightGrid, ({ positions, indices }) => {
            const vcm = new VertexCollapseMesh(positions, indices);
            for (const v of vcm.getMesh().getVertices()) {
                const vertex = v as VertexCollapseMeshVertex;
                const { weight, normal } = referenceWeight(vertex, positions);
                expect(vertex.computeWeight(positions)).toBe(weight);
                expect(vertex.normal.values).toEqual(normal.values);
                // The stored normal is unit length for a nondegenerate link.
                expect(Math.abs(length(vertex.normal) - 1))
                    .toBeLessThanOrEqual(1e-12);
            }
        });
    });

    it('collapses the smallest-weight interior vertex first', () => {
        // The constructor gives boundary vertices the weight MAX_VALUE, so
        // the min-heap pops an interior vertex of minimal weight. The first
        // successful collapse must be one of those (a deferred vertex is
        // pushed back with MAX_VALUE, which can only come later).
        check(heightGrid, ({ n, positions, indices }) => {
            const vcm = new VertexCollapseMesh(positions, indices);
            const mesh = vcm.getMesh();
            const boundary = new Set<number>();
            for (let r = 0; r < n; ++r) {
                for (let c = 0; c < n; ++c) {
                    if (r === 0 || c === 0 || r === n - 1 || c === n - 1) {
                        boundary.add(r * n + c);
                    }
                }
            }
            let best = Number.MAX_VALUE;
            for (const v of mesh.getVertices()) {
                if (boundary.has(v.V)) { continue; }
                best = Math.min(best,
                    (v as VertexCollapseMeshVertex).computeWeight(positions));
            }
            const result = vcm.doCollapse();
            expect(result.collapsed).toBe(true);
            expect(boundary.has(result.vertex)).toBe(false);
            // The weight of the collapsed vertex was the minimum, unless a
            // smaller-weight vertex was deferred first (which cannot happen
            // on the very first call, because nothing has been deferred yet).
            const weights = mesh.getVertices();
            expect(weights.length).toBeGreaterThan(0);
            expect(best).toBeLessThan(Number.MAX_VALUE);
        });
    });

    it('reports no collapse for degenerate construction input', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 2 }),
            fc.integer({ min: 0, max: 2 })), ([numPositions, numIndices]) => {
            const positions: Vector[] = [];
            for (let i = 0; i < numPositions; ++i) {
                positions.push(Vector.fromArray([i, 0, 0]));
            }
            const indices = new Array<number>(numIndices).fill(0);
            for (const [p, i] of [[positions, indices], [null, indices],
                [positions, null], [null, null]] as
                Array<[Vector[] | null, number[] | null]>) {
                const vcm = new VertexCollapseMesh(p, i);
                const result = vcm.doCollapse();
                expect(result.collapsed).toBe(false);
                expect(result.vertex).toBe(VERTEX_COLLAPSE_MESH_INVALID_VERTEX);
                expect(result.removed).toEqual([]);
                expect(result.inserted).toEqual([]);
            }
        });
    });

    it('stops without changing the mesh once no collapse is allowed', () => {
        check(heightGrid, ({ positions, indices }) => {
            const vcm = new VertexCollapseMesh(positions, indices);
            const mesh = vcm.getMesh();
            collapseAll(vcm, positions.length + 5);
            const keys = mesh.getTriangleKeys().map(k => k.V.join(','));
            for (let i = 0; i < 3; ++i) {
                const result = vcm.doCollapse();
                expect(result.collapsed).toBe(false);
                expect(mesh.getTriangleKeys().map(k => k.V.join(',')))
                    .toEqual(keys);
            }
        });
    }, 30000);
});
