import { describe, it, expect } from 'vitest';
import {
    AABBBoundingVolume, AABBBVTreeOfTriangles, aabbBoundingVolumeOps
} from '../src/AABBBVTreeOfTriangles';
import { AlignedBox } from '../src/AlignedBox';
import { BVTree, BVTreeNode } from '../src/BVTree';
import {
    intersectLineTriangle, intersectRayTriangle, intersectSegmentTriangle
} from '../src/BVTreeOfTriangles';
import type { LinearTriangleResult } from '../src/BVTreeOfTriangles';
import { Triangle } from '../src/Triangle';
import { Vector, normalize } from '../src/Vector';

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

type Tri = [number, number, number];

interface Mesh {
    vertices: Vector[];
    triangles: Tri[];
}

// A closed axis-aligned box mesh with each face split into an n-by-n grid of
// quads, each quad split into two triangles. Every edge is shared by exactly
// two triangles and the mesh has 12*n*n triangles.
function makeBoxMesh(n: number, lo: Vector, hi: Vector): Mesh {
    const vertices: Vector[] = [];
    const triangles: Tri[] = [];
    const index = new Map<string, number>();

    const at = (t: number, lo0: number, hi0: number): number =>
        lo0 + (hi0 - lo0) * (t / n);

    const addVertex = (x: number, y: number, z: number): number => {
        const key = x + ',' + y + ',' + z;
        const found = index.get(key);
        if (found !== undefined) {
            return found;
        }
        const i = vertices.length;
        vertices.push(v3(x, y, z));
        index.set(key, i);
        return i;
    };

    for (let axis = 0; axis < 3; ++axis) {
        const a1 = (axis + 1) % 3;
        const a2 = (axis + 2) % 3;
        for (let side = 0; side < 2; ++side) {
            const fixed = side === 0 ? lo.get(axis) : hi.get(axis);
            for (let i = 0; i < n; ++i) {
                for (let j = 0; j < n; ++j) {
                    const corner = (di: number, dj: number): number => {
                        const c = [0, 0, 0];
                        c[axis] = fixed;
                        c[a1] = at(i + di, lo.get(a1), hi.get(a1));
                        c[a2] = at(j + dj, lo.get(a2), hi.get(a2));
                        return addVertex(c[0], c[1], c[2]);
                    };
                    const v00 = corner(0, 0);
                    const v10 = corner(1, 0);
                    const v11 = corner(1, 1);
                    const v01 = corner(0, 1);
                    // Winding is irrelevant here; the queries are two-sided.
                    triangles.push([v00, v10, v11]);
                    triangles.push([v00, v11, v01]);
                }
            }
        }
    }

    return { vertices: vertices, triangles: triangles };
}

// A tetrahedron: the smallest closed mesh.
const tetra: Mesh = {
    vertices: [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)],
    triangles: [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
};

// The component-wise min/max over the vertices of a set of triangles: the
// tight axis-aligned box, computed independently of the port.
function tightBox(mesh: Mesh, triangleIndices: readonly number[]): AlignedBox {
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];
    for (const t of triangleIndices) {
        const tri = mesh.triangles[t];
        for (let j = 0; j < 3; ++j) {
            const p = mesh.vertices[tri[j]];
            for (let k = 0; k < 3; ++k) {
                min[k] = Math.min(min[k], p.get(k));
                max[k] = Math.max(max[k], p.get(k));
            }
        }
    }
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

// A deterministic pseudorandom generator, so failures reproduce.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

const queryFns: ((P: Vector, Q: Vector, t: Triangle) => LinearTriangleResult)[] = [
    intersectLineTriangle, intersectRayTriangle, intersectSegmentTriangle
];

interface Hit {
    triangleIndex: number;
    point: Vector;
    parameter: number;
}

// The brute-force reference for execute(): test every triangle and collect the
// hits with the same std::set-ordered-by-parameter-only container semantics
// that the port reproduces (upstream #167). Insertion order is triangle order.
function bruteForceExecute(mesh: Mesh, queryType: number, P: Vector, Q: Vector): Hit[] {
    const query = queryFns[queryType];
    const out: Hit[] = [];
    for (let t = 0; t < mesh.triangles.length; ++t) {
        const tri = mesh.triangles[t];
        const triangle = Triangle.fromVertices(mesh.vertices[tri[0]],
            mesh.vertices[tri[1]], mesh.vertices[tri[2]]);
        const result = query(P, Q, triangle);
        if (!result.intersect) {
            continue;
        }
        let lo = 0;
        let hi = out.length;
        while (lo < hi) {
            const mid = lo + Math.floor((hi - lo) / 2);
            if (out[mid].parameter < result.parameter) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        if (lo < out.length && !(result.parameter < out[lo].parameter)) {
            // An equivalent element is already a member of the set.
            continue;
        }
        out.splice(lo, 0, {
            triangleIndex: t, point: result.point, parameter: result.parameter
        });
    }
    return out;
}

// The triangles a tree reports as candidates: the union over the reported leaf
// nodes of their partition ranges.
function candidateTriangles(tree: AABBBVTreeOfTriangles,
    nodeIndices: readonly number[]): Set<number> {
    const nodes = tree.getNodes();
    const partition = tree.getPartition();
    const result = new Set<number>();
    for (const nodeIndex of nodeIndices) {
        const node = nodes[nodeIndex];
        for (let i = node.minIndex; i <= node.maxIndex; ++i) {
            result.add(partition[i]);
        }
    }
    return result;
}

// Exposes the protected getLeafIndices, the pattern used by
// test/BVTreeOfTriangles.test.ts.
class TestableTree extends AABBBVTreeOfTriangles {
    leafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }
}

function buildTree(mesh: Mesh, height?: number): TestableTree {
    const tree = new TestableTree();
    tree.createFromTriangles(mesh.vertices, mesh.triangles,
        height === undefined ? BVTree.fullHeight : height);
    return tree;
}

function segmentEnd(P: Vector, D: Vector, scale: number): Vector {
    return Vector.fromArray(P.values.map((c, k) => c + scale * D.get(k)));
}

// ---------------------------------------------------------------------------
// AABBBoundingVolume
// ---------------------------------------------------------------------------

describe('AABBBoundingVolume', () => {
    it('default-constructs the AlignedBox3 default (min = -1, max = +1)', () => {
        const bv = new AABBBoundingVolume();
        expect(bv.box.dimension).toBe(3);
        expect(bv.box.min.values).toEqual([-1, -1, -1]);
        expect(bv.box.max.values).toEqual([1, 1, 1]);
    });

    it('splits along the axis of largest extent, through the box center', () => {
        const bv = new AABBBoundingVolume();

        bv.box = AlignedBox.fromMinMax(v3(-4, 0, 0), v3(6, 1, 2));
        let axis = bv.getSplittingAxis();
        expect(axis.origin.values).toEqual([1, 0.5, 1]);
        expect(axis.direction.values).toEqual([1, 0, 0]);

        bv.box = AlignedBox.fromMinMax(v3(0, -10, 0), v3(1, 10, 2));
        axis = bv.getSplittingAxis();
        expect(axis.origin.values).toEqual([0.5, 0, 1]);
        expect(axis.direction.values).toEqual([0, 1, 0]);

        bv.box = AlignedBox.fromMinMax(v3(0, 0, -3), v3(1, 2, 7));
        axis = bv.getSplittingAxis();
        expect(axis.origin.values).toEqual([0.5, 1, 2]);
        expect(axis.direction.values).toEqual([0, 0, 1]);
    });

    it('breaks extent ties in favor of the smaller axis index', () => {
        // The upstream comparisons are strict, so equal extents keep the
        // earlier axis.
        const bv = new AABBBoundingVolume();
        bv.box = AlignedBox.fromMinMax(v3(0, 0, 0), v3(2, 2, 2));
        expect(bv.getSplittingAxis().direction.values).toEqual([1, 0, 0]);

        bv.box = AlignedBox.fromMinMax(v3(0, 0, 0), v3(1, 2, 2));
        expect(bv.getSplittingAxis().direction.values).toEqual([0, 1, 0]);
    });

    it('handles a degenerate (zero-extent) box', () => {
        const bv = new AABBBoundingVolume();
        bv.box = AlignedBox.fromMinMax(v3(3, 4, 5), v3(3, 4, 5));
        const axis = bv.getSplittingAxis();
        expect(axis.origin.values).toEqual([3, 4, 5]);
        expect(axis.direction.values).toEqual([1, 0, 0]);
    });

    it('tests the linear components against the box', () => {
        const bv = new AABBBoundingVolume();
        bv.box = AlignedBox.fromMinMax(v3(0, 0, 0), v3(1, 1, 1));

        const P = v3(-5, 0.5, 0.5);
        const Q = v3(1, 0, 0);
        expect(AABBBoundingVolume.intersectLine(P, Q, bv)).toBe(true);
        expect(AABBBoundingVolume.intersectRay(P, Q, bv)).toBe(true);
        // The segment endpoints are P and Q, not P and P+Q.
        expect(AABBBoundingVolume.intersectSegment(P, v3(5, 0.5, 0.5), bv)).toBe(true);
        expect(AABBBoundingVolume.intersectSegment(P, v3(-2, 0.5, 0.5), bv)).toBe(false);

        // A ray pointing away from the box misses; the line through the same
        // point still hits.
        expect(AABBBoundingVolume.intersectRay(P, v3(-1, 0, 0), bv)).toBe(false);
        expect(AABBBoundingVolume.intersectLine(P, v3(-1, 0, 0), bv)).toBe(true);

        // A line that misses entirely.
        const M = v3(-5, 3, 0.5);
        expect(AABBBoundingVolume.intersectLine(M, Q, bv)).toBe(false);
        expect(AABBBoundingVolume.intersectRay(M, Q, bv)).toBe(false);
        expect(AABBBoundingVolume.intersectSegment(M, v3(5, 3, 0.5), bv)).toBe(false);
    });

    it('bundles the static operations into aabbBoundingVolumeOps', () => {
        const created = aabbBoundingVolumeOps.create();
        expect(created).toBeInstanceOf(AABBBoundingVolume);
        expect(created.box.max.values).toEqual([1, 1, 1]);

        const bv = new AABBBoundingVolume();
        bv.box = AlignedBox.fromMinMax(v3(0, 0, 0), v3(1, 1, 1));
        const P = v3(-5, 0.5, 0.5);
        const Q = v3(1, 0, 0);
        expect(aabbBoundingVolumeOps.intersectLine(P, Q, bv)).toBe(true);
        expect(aabbBoundingVolumeOps.intersectRay(P, Q, bv)).toBe(true);
        expect(aabbBoundingVolumeOps.intersectSegment(P, v3(5, 0.5, 0.5), bv)).toBe(true);
        expect(aabbBoundingVolumeOps.intersectRay(P, v3(-1, 0, 0), bv)).toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Tree construction
// ---------------------------------------------------------------------------

describe('AABBBVTreeOfTriangles construction', () => {
    it('builds a single-node tree for one triangle with a tight leaf box', () => {
        const tree = buildTree({
            vertices: [v3(0, 0, 0), v3(2, 0, 0), v3(0, 3, 5)],
            triangles: [[0, 1, 2]]
        });
        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        const root = tree.getNodes()[0];
        expect(root.leftChild).toBe(BVTreeNode.invalid);
        expect(root.boundingVolume.box.min.values).toEqual([0, 0, 0]);
        expect(root.boundingVolume.box.max.values).toEqual([2, 3, 5]);
    });

    it('gives the root the tight box of the whole tetrahedron', () => {
        const tree = buildTree(tetra);
        const root = tree.getNodes()[0];
        expect(root.boundingVolume.box.min.values).toEqual([0, 0, 0]);
        expect(root.boundingVolume.box.max.values).toEqual([1, 1, 1]);
    });

    it('gives every visited node the tight box of its triangle range', () => {
        const meshes: Mesh[] = [
            tetra,
            makeBoxMesh(1, v3(-1, -2, -3), v3(4, 1, 2)),
            makeBoxMesh(2, v3(0, 0, 0), v3(1, 1, 1)),
            makeBoxMesh(3, v3(-2, -1, 0), v3(2, 5, 1))
        ];
        for (const mesh of meshes) {
            const tree = buildTree(mesh);
            const partition = tree.getPartition();
            let visited = 0;
            for (const node of tree.getNodes()) {
                if (node.minIndex === BVTreeNode.invalid) {
                    continue;
                }
                ++visited;
                const range: number[] = [];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    range.push(partition[i]);
                }
                const expected = tightBox(mesh, range);
                expect(node.boundingVolume.box.min.values).toEqual(expected.min.values);
                expect(node.boundingVolume.box.max.values).toEqual(expected.max.values);
            }
            expect(visited).toBeGreaterThanOrEqual(mesh.triangles.length);
        }
    });

    it('gives every leaf the tight box of its single triangle', () => {
        const mesh = makeBoxMesh(2, v3(-1, 0, 2), v3(3, 4, 5));
        const tree = buildTree(mesh);
        const partition = tree.getPartition();
        let leafCount = 0;
        for (const node of tree.getNodes()) {
            if (node.minIndex === BVTreeNode.invalid ||
                node.minIndex !== node.maxIndex) {
                continue;
            }
            ++leafCount;
            const expected = tightBox(mesh, [partition[node.minIndex]]);
            expect(node.boundingVolume.box.min.values).toEqual(expected.min.values);
            expect(node.boundingVolume.box.max.values).toEqual(expected.max.values);
        }
        expect(leafCount).toBe(mesh.triangles.length);
    });

    it('nests child boxes inside parent boxes', () => {
        const mesh = makeBoxMesh(3, v3(0, 0, 0), v3(2, 3, 5));
        const tree = buildTree(mesh);
        const nodes = tree.getNodes();
        let interiorCount = 0;
        for (const node of nodes) {
            if (node.leftChild === BVTreeNode.invalid) {
                continue;
            }
            ++interiorCount;
            for (const c of [node.leftChild, node.rightChild]) {
                const child = nodes[c];
                for (let k = 0; k < 3; ++k) {
                    expect(child.boundingVolume.box.min.get(k))
                        .toBeGreaterThanOrEqual(node.boundingVolume.box.min.get(k));
                    expect(child.boundingVolume.box.max.get(k))
                        .toBeLessThanOrEqual(node.boundingVolume.box.max.get(k));
                }
            }
        }
        expect(interiorCount).toBeGreaterThan(0);
    });

    it('respects a user-specified height, leaving multi-triangle leaves tight', () => {
        const mesh = makeBoxMesh(2, v3(0, 0, 0), v3(1, 1, 1));
        const full = buildTree(mesh);
        const shallow = buildTree(mesh, 2);
        expect(shallow.getHeight()).toBe(2);
        expect(full.getHeight()).toBeGreaterThan(2);

        const partition = shallow.getPartition();
        let multiLeaf = 0;
        for (const node of shallow.getNodes()) {
            if (node.minIndex === BVTreeNode.invalid ||
                node.leftChild !== BVTreeNode.invalid) {
                continue;
            }
            if (node.maxIndex > node.minIndex) {
                ++multiLeaf;
            }
            const range: number[] = [];
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                range.push(partition[i]);
            }
            const expected = tightBox(mesh, range);
            expect(node.boundingVolume.box.min.values).toEqual(expected.min.values);
            expect(node.boundingVolume.box.max.values).toEqual(expected.max.values);
        }
        expect(multiLeaf).toBeGreaterThan(0);
    });

    it('rejects degenerate creation input', () => {
        const tree = new AABBBVTreeOfTriangles();
        expect(() => tree.createFromTriangles([v3(0, 0, 0), v3(1, 0, 0)], [[0, 1, 1]]))
            .toThrow();
        expect(() => tree.createFromTriangles(tetra.vertices, [])).toThrow();
    });

    it('copies the input vertices and triangles', () => {
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const triangles: Tri[] = [[0, 1, 2]];
        const tree = buildTree({ vertices: vertices, triangles: triangles });
        vertices[0].set(0, 100);
        triangles[0][0] = 2;
        expect(tree.getVertices()[0].values).toEqual([0, 0, 0]);
        expect(tree.getTriangles()[0]).toEqual([0, 1, 2]);
    });
});

// ---------------------------------------------------------------------------
// Queries
// ---------------------------------------------------------------------------

describe('AABBBVTreeOfTriangles queries', () => {
    it('reports leaves conservatively: every hit triangle is a candidate', () => {
        const mesh = makeBoxMesh(3, v3(-1, -1, -1), v3(1, 1, 1));
        const tree = buildTree(mesh);
        const rand = makeRandom(20250901);
        let checked = 0;
        for (let trial = 0; trial < 250; ++trial) {
            const P = v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2);
            const D = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (normalize(D) === 0) {
                continue;
            }
            for (const queryType of [BVTree.LINE_QUERY, BVTree.RAY_QUERY,
                BVTree.SEGMENT_QUERY]) {
                const Q = queryType === BVTree.SEGMENT_QUERY
                    ? segmentEnd(P, D, 6) : D;
                const candidates = candidateTriangles(tree,
                    tree.leafIndices(queryType, P, Q));
                for (const hit of bruteForceExecute(mesh, queryType, P, Q)) {
                    expect(candidates.has(hit.triangleIndex)).toBe(true);
                    ++checked;
                }
            }
        }
        expect(checked).toBeGreaterThan(200);
    });

    it('prunes: the candidate set is smaller than the whole mesh', () => {
        const mesh = makeBoxMesh(4, v3(-1, -1, -1), v3(1, 1, 1));
        const tree = buildTree(mesh);
        const rand = makeRandom(7);
        let pruned = 0;
        let total = 0;
        for (let trial = 0; trial < 100; ++trial) {
            const P = v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3);
            const D = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (normalize(D) === 0) {
                continue;
            }
            ++total;
            const candidates = candidateTriangles(tree,
                tree.leafIndices(BVTree.LINE_QUERY, P, D));
            if (candidates.size < mesh.triangles.length) {
                ++pruned;
            }
        }
        expect(total).toBeGreaterThan(90);
        expect(pruned).toBe(total);
    });

    it('matches brute force on random lines, rays and segments', () => {
        const mesh = makeBoxMesh(3, v3(-1, -2, -1.5), v3(2, 1, 3));
        const tree = buildTree(mesh);
        const rand = makeRandom(424242);
        let hits = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const P = v3(8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4);
            const D = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (normalize(D) === 0) {
                continue;
            }
            for (const queryType of [BVTree.LINE_QUERY, BVTree.RAY_QUERY,
                BVTree.SEGMENT_QUERY]) {
                const Q = queryType === BVTree.SEGMENT_QUERY
                    ? segmentEnd(P, D, 12) : D;
                const expected = bruteForceExecute(mesh, queryType, P, Q);
                const actual = tree.execute(queryType, P, Q);
                expect(actual.intersections.length).toBe(expected.length);
                for (let i = 0; i < expected.length; ++i) {
                    expect(actual.intersections[i].parameter)
                        .toBeCloseTo(expected[i].parameter, 12);
                    expect(actual.intersections[i].triangleIndex)
                        .toBe(expected[i].triangleIndex);
                    for (let k = 0; k < 3; ++k) {
                        expect(actual.intersections[i].point.get(k))
                            .toBeCloseTo(expected[i].point.get(k), 12);
                    }
                }
                hits += expected.length;
            }
        }
        expect(hits).toBeGreaterThan(100);
    });

    it('matches brute force on a shallow tree with multi-triangle leaves', () => {
        const mesh = makeBoxMesh(2, v3(0, 0, 0), v3(1, 1, 1));
        const tree = buildTree(mesh, 2);
        const rand = makeRandom(999);
        let hits = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const P = v3(3 * rand() - 1, 3 * rand() - 1, 3 * rand() - 1);
            const D = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (normalize(D) === 0) {
                continue;
            }
            const expected = bruteForceExecute(mesh, BVTree.LINE_QUERY, P, D);
            const actual = tree.execute(BVTree.LINE_QUERY, P, D);
            expect(actual.intersections.length).toBe(expected.length);
            for (let i = 0; i < expected.length; ++i) {
                expect(actual.intersections[i].parameter)
                    .toBeCloseTo(expected[i].parameter, 12);
            }
            hits += expected.length;
        }
        expect(hits).toBeGreaterThan(40);
    });

    it('finds the two crossings of a line through a closed mesh', () => {
        const mesh = makeBoxMesh(2, v3(-1, -1, -1), v3(1, 1, 1));
        const tree = buildTree(mesh);
        // (y, z) = (0.25, -0.375) avoids every grid line and every quad
        // diagonal, so the line crosses exactly two faces.
        const P = v3(-10, 0.25, -0.375);
        const Q = v3(1, 0, 0);
        const result = tree.execute(BVTree.LINE_QUERY, P, Q);
        expect(result.intersections.length).toBe(2);
        expect(result.intersections[0].parameter).toBeCloseTo(9, 12);
        expect(result.intersections[1].parameter).toBeCloseTo(11, 12);
        expect(result.intersections[0].point.values[0]).toBeCloseTo(-1, 12);
        expect(result.intersections[1].point.values[0]).toBeCloseTo(1, 12);
    });

    it('reports nothing for a line that misses the mesh', () => {
        const mesh = makeBoxMesh(2, v3(-1, -1, -1), v3(1, 1, 1));
        const tree = buildTree(mesh);
        const result = tree.execute(BVTree.LINE_QUERY, v3(-10, 5, 0), v3(1, 0, 0));
        expect(result.intersections.length).toBe(0);
        // The root box test fails, so the traversal reports no leaves at all.
        expect(result.nodeIndices.length).toBe(0);
    });

    it('distinguishes ray and segment restrictions from the line', () => {
        const mesh = makeBoxMesh(1, v3(-1, -1, -1), v3(1, 1, 1));
        const tree = buildTree(mesh);
        const P = v3(-10, 0.25, -0.375);
        const D = v3(1, 0, 0);

        expect(tree.execute(BVTree.LINE_QUERY, P, D).intersections.length).toBe(2);
        expect(tree.execute(BVTree.RAY_QUERY, P, D).intersections.length).toBe(2);
        expect(tree.execute(BVTree.RAY_QUERY, P, v3(-1, 0, 0))
            .intersections.length).toBe(0);
        // A segment ending inside the mesh crosses only the near face.
        expect(tree.execute(BVTree.SEGMENT_QUERY, P, v3(0, 0.25, -0.375))
            .intersections.length).toBe(1);
        // A segment fully outside crosses nothing.
        expect(tree.execute(BVTree.SEGMENT_QUERY, P, v3(-5, 0.25, -0.375))
            .intersections.length).toBe(0);
    });

    it('reports the root leaf of a height-0 tree for every query (upstream #103)', () => {
        // BVTree.getLeafIndices never tests a leaf's own bounding volume, so
        // the sole node of a single-triangle tree is reported even for a query
        // far away. The exact triangle test then rejects it.
        const tree = buildTree({
            vertices: [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)],
            triangles: [[0, 1, 2]]
        });
        const far = v3(1000, 1000, 1000);
        const result = tree.execute(BVTree.LINE_QUERY, far, v3(1, 0, 0));
        expect(result.nodeIndices).toEqual([0]);
        expect(result.intersections.length).toBe(0);
    });

    it('drops coincident-parameter hits (upstream #167)', () => {
        // The quad diagonals of makeBoxMesh(1) run from the (0,0) corner to
        // the (1,1) corner in face coordinates, so a line along +z at
        // (x, y) = (0.5, 0.5) passes through the shared diagonal of both
        // triangles of the z = 0 face and of the z = 1 face: four triangle
        // hits at two distinct parameters.
        const mesh = makeBoxMesh(1, v3(0, 0, 0), v3(1, 1, 1));
        const tree = buildTree(mesh);
        const P = v3(0.5, 0.5, -3);
        const Q = v3(0, 0, 1);

        // Confirm the setup: four triangles genuinely contain a crossing.
        let rawHits = 0;
        for (const tri of mesh.triangles) {
            const triangle = Triangle.fromVertices(mesh.vertices[tri[0]],
                mesh.vertices[tri[1]], mesh.vertices[tri[2]]);
            if (intersectLineTriangle(P, Q, triangle).intersect) {
                ++rawHits;
            }
        }
        expect(rawHits).toBe(4);

        const result = tree.execute(BVTree.LINE_QUERY, P, Q);
        expect(result.intersections.length).toBe(2);
        const parameters = result.intersections.map(x => x.parameter);
        expect(new Set(parameters).size).toBe(2);
        expect(parameters[0]).toBeCloseTo(3, 12);
        expect(parameters[1]).toBeCloseTo(4, 12);
    });

    it('agrees with brute force on a tetrahedron over many random rays', () => {
        const tree = buildTree(tetra);
        const rand = makeRandom(31337);
        let hits = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const P = v3(3 * rand() - 1, 3 * rand() - 1, 3 * rand() - 1);
            const D = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (normalize(D) === 0) {
                continue;
            }
            const expected = bruteForceExecute(tetra, BVTree.RAY_QUERY, P, D);
            const actual = tree.execute(BVTree.RAY_QUERY, P, D);
            expect(actual.intersections.map(x => x.triangleIndex))
                .toEqual(expected.map(x => x.triangleIndex));
            hits += expected.length;
        }
        expect(hits).toBeGreaterThan(20);
    });

    it('keeps the intersections sorted, on the line and inside their triangle box', () => {
        const mesh = makeBoxMesh(3, v3(-1, -1, -1), v3(1, 1, 1));
        const tree = buildTree(mesh);
        const rand = makeRandom(5150);
        for (let trial = 0; trial < 100; ++trial) {
            const P = v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3);
            const D = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (normalize(D) === 0) {
                continue;
            }
            const result = tree.execute(BVTree.LINE_QUERY, P, D);
            for (let i = 1; i < result.intersections.length; ++i) {
                expect(result.intersections[i - 1].parameter)
                    .toBeLessThan(result.intersections[i].parameter);
            }
            for (const hit of result.intersections) {
                const onLine = segmentEnd(P, D, hit.parameter);
                const box = tightBox(mesh, [hit.triangleIndex]);
                for (let k = 0; k < 3; ++k) {
                    expect(hit.point.get(k)).toBeCloseTo(onLine.get(k), 10);
                    expect(hit.point.get(k))
                        .toBeGreaterThanOrEqual(box.min.get(k) - 1e-9);
                    expect(hit.point.get(k))
                        .toBeLessThanOrEqual(box.max.get(k) + 1e-9);
                }
            }
        }
    });

    it('is unaffected by the order of the input triangles', () => {
        const mesh = makeBoxMesh(2, v3(-1, -1, -1), v3(1, 1, 1));
        const permuted: Mesh = {
            vertices: mesh.vertices,
            triangles: mesh.triangles.slice().reverse()
        };
        const treeA = buildTree(mesh);
        const treeB = buildTree(permuted);
        const rand = makeRandom(64738);
        for (let trial = 0; trial < 60; ++trial) {
            const P = v3(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3);
            const D = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (normalize(D) === 0) {
                continue;
            }
            const a = treeA.execute(BVTree.LINE_QUERY, P, D).intersections
                .map(x => x.parameter);
            const b = treeB.execute(BVTree.LINE_QUERY, P, D).intersections
                .map(x => x.parameter);
            expect(a.length).toBe(b.length);
            for (let i = 0; i < a.length; ++i) {
                expect(a[i]).toBeCloseTo(b[i], 12);
            }
        }
    });
});
