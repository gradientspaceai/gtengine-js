import { describe, it, expect } from 'vitest';
import { AlignedBoxTreeOfTriangles } from '../src/AlignedBoxTreeOfTriangles.js';
import { AlignedBoxBV } from '../src/AlignedBoxBV.js';
import { BVTree, BVTreeNode } from '../src/BVTree.js';
import {
    intersectLineTriangle, intersectRayTriangle, intersectSegmentTriangle
} from '../src/BVTreeOfTriangles.js';
import type { LinearTriangleResult } from '../src/BVTreeOfTriangles.js';
import { Triangle } from '../src/Triangle.js';
import { Vector } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

type Tri = [number, number, number];

interface Mesh {
    vertices: Vector[];
    triangles: Tri[];
}

// A closed axis-aligned box mesh with each face split into an n-by-n grid of
// quads, each quad split into two triangles.
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

const tetra: Mesh = {
    vertices: [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)],
    triangles: [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
};

// The tight axis-aligned bounds over a set of triangles, computed
// independently of the port.
function tightBounds(mesh: Mesh, triangleIndices: readonly number[]):
    { min: number[], max: number[] } {
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];
    for (const t of triangleIndices) {
        for (const j of mesh.triangles[t]) {
            for (let k = 0; k < 3; ++k) {
                min[k] = Math.min(min[k], mesh.vertices[j].get(k));
                max[k] = Math.max(max[k], mesh.vertices[j].get(k));
            }
        }
    }
    return { min: min, max: max };
}

const queryFns: ((P: Vector, Q: Vector, t: Triangle) => LinearTriangleResult)[] = [
    intersectLineTriangle, intersectRayTriangle, intersectSegmentTriangle
];

interface Hit {
    triangleIndex: number;
    parameter: number;
}

// The brute-force reference for execute(): test every triangle and collect
// every hit ordered by (parameter, triangleIndex) -- the port's corrected
// container semantics (upstream #167 is fixed in the port).
function bruteForceExecute(mesh: Mesh, queryType: number, P: Vector,
    Q: Vector): Hit[] {
    const query = queryFns[queryType];
    const out: Hit[] = [];
    for (let t = 0; t < mesh.triangles.length; ++t) {
        const tri = mesh.triangles[t];
        const triangle = Triangle.fromVertices(mesh.vertices[tri[0]],
            mesh.vertices[tri[1]], mesh.vertices[tri[2]]);
        const result = query(P, Q, triangle);
        if (result.intersect) {
            out.push({ triangleIndex: t, parameter: result.parameter });
        }
    }
    out.sort((a, b) => a.parameter - b.parameter ||
        a.triangleIndex - b.triangleIndex);
    return out;
}

class TestableTree extends AlignedBoxTreeOfTriangles {
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

// An independent, recursive reference for the stack-based getLeafIndices.
function referenceLeaves(tree: TestableTree, queryType: number, P: Vector,
    Q: Vector): number[] {
    const nodes = tree.getNodes();
    const hits = (bv: AlignedBoxBV): boolean =>
        queryType === BVTree.LINE_QUERY ? AlignedBoxBV.intersectLine(P, Q, bv)
            : (queryType === BVTree.RAY_QUERY
                ? AlignedBoxBV.intersectRay(P, Q, bv)
                : AlignedBoxBV.intersectSegment(P, Q, bv));

    const out: number[] = [];
    const visit = (nodeIndex: number): void => {
        const node = nodes[nodeIndex];
        if (node.leftChild !== BVTreeNode.invalid &&
            node.rightChild !== BVTreeNode.invalid) {
            if (hits(node.boundingVolume)) {
                visit(node.leftChild);
                visit(node.rightChild);
            }
        } else {
            out.push(nodeIndex);
        }
    };
    visit(0);
    return out;
}

describe('AlignedBoxTreeOfTriangles construction', () => {
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
        const root = buildTree(tetra).getNodes()[0];
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
            let leafCount = 0;
            for (const node of tree.getNodes()) {
                if (node.minIndex === BVTreeNode.invalid) {
                    continue;
                }
                const range: number[] = [];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    range.push(partition[i]);
                }
                const expected = tightBounds(mesh, range);
                expect(node.boundingVolume.box.min.values).toEqual(expected.min);
                expect(node.boundingVolume.box.max.values).toEqual(expected.max);
                if (node.minIndex === node.maxIndex) {
                    ++leafCount;
                }
            }
            expect(leafCount).toBe(mesh.triangles.length);
        }
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
            const expected = tightBounds(mesh, range);
            expect(node.boundingVolume.box.min.values).toEqual(expected.min);
            expect(node.boundingVolume.box.max.values).toEqual(expected.max);
        }
        expect(multiLeaf).toBeGreaterThan(0);
    });

    it('handles a degenerate (zero-area) triangle', () => {
        const mesh: Mesh = {
            vertices: [v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2), v3(0, 1, 0)],
            triangles: [[0, 1, 2], [0, 1, 3]]
        };
        const tree = buildTree(mesh);
        const root = tree.getNodes()[0];
        expect(root.boundingVolume.box.min.values).toEqual([0, 0, 0]);
        expect(root.boundingVolume.box.max.values).toEqual([2, 2, 2]);
    });

    it('does not alias the input vertices in the node boxes', () => {
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const tree = buildTree({ vertices: vertices, triangles: [[0, 1, 2]] });
        tree.getNodes()[0].boundingVolume.box.min.set(0, -100);
        expect(vertices[0].get(0)).toBe(0);
        expect(tree.getVertices()[0].get(0)).toBe(0);
    });

    it('rejects degenerate creation input', () => {
        const tree = new AlignedBoxTreeOfTriangles();
        expect(() => tree.createFromTriangles([v3(0, 0, 0), v3(1, 0, 0)], [[0, 1, 1]]))
            .toThrow();
        expect(() => tree.createFromTriangles(tetra.vertices, [])).toThrow();
    });
});

describe('AlignedBoxTreeOfTriangles queries', () => {
    it('finds the same intersections as brute force (known values)', () => {
        const mesh = makeBoxMesh(1, v3(0, 0, 0), v3(1, 1, 1));
        const tree = buildTree(mesh);

        // A line through the box, entering at x = 0 and exiting at x = 1.
        const P = v3(-2, 0.3, 0.7);
        const Q = v3(1, 0, 0);
        const result = tree.execute(BVTree.LINE_QUERY, P, Q);
        expect(result.intersections.length).toBe(2);
        expect(result.intersections[0].parameter).toBeCloseTo(2, 12);
        expect(result.intersections[1].parameter).toBeCloseTo(3, 12);
        expect(result.intersections[0].point.get(0)).toBeCloseTo(0, 12);
        expect(result.intersections[1].point.get(0)).toBeCloseTo(1, 12);
        expect(result.intersections.map(h => h.triangleIndex))
            .toEqual(bruteForceExecute(mesh, BVTree.LINE_QUERY, P, Q)
                .map(h => h.triangleIndex));

        // A ray starting inside the box hits only the exit face.
        const inside = v3(0.5, 0.3, 0.7);
        const rayResult = tree.execute(BVTree.RAY_QUERY, inside, Q);
        expect(rayResult.intersections.length).toBe(1);
        expect(rayResult.intersections[0].parameter).toBeCloseTo(0.5, 12);

        // A segment that stops short of the box has no intersections.
        expect(tree.execute(BVTree.SEGMENT_QUERY, v3(-2, 0.3, 0.7),
            v3(-1, 0.3, 0.7)).intersections.length).toBe(0);
    });

    it('keeps coincident hits at a shared edge (upstream #167 fixed)', () => {
        // The two triangles of a quad share the diagonal edge from (0,0,0) to
        // (1,1,0). A line through a point of that edge hits both triangles at
        // the same parameter; upstream's std::set ordered by parameter alone
        // would drop one of them.
        const mesh: Mesh = {
            vertices: [v3(0, 0, 0), v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0)],
            triangles: [[0, 1, 2], [0, 2, 3]]
        };
        const tree = buildTree(mesh);
        const P = v3(0.5, 0.5, -1);
        const Q = v3(0, 0, 1);
        const result = tree.execute(BVTree.LINE_QUERY, P, Q);
        expect(result.intersections.length).toBe(2);
        expect(result.intersections[0].parameter)
            .toBe(result.intersections[1].parameter);
        expect(result.intersections.map(h => h.triangleIndex)).toEqual([0, 1]);
    });

    it('agrees with brute force on randomized queries', () => {
        const meshes: Mesh[] = [
            tetra,
            makeBoxMesh(2, v3(-1, -1, -1), v3(1, 1, 1)),
            makeBoxMesh(3, v3(0, 0, 0), v3(2, 1, 3))
        ];
        const rand = makeRandom(2024);
        for (const mesh of meshes) {
            const tree = buildTree(mesh);
            for (const queryType of [BVTree.LINE_QUERY, BVTree.RAY_QUERY,
                BVTree.SEGMENT_QUERY]) {
                for (let trial = 0; trial < 25; ++trial) {
                    const P = v3(6 * (2 * rand() - 1), 6 * (2 * rand() - 1),
                        6 * (2 * rand() - 1));
                    const target = v3(2 * rand() - 1, 2 * rand() - 1,
                        2 * rand() - 1);
                    const d = Vector.fromArray(
                        target.values.map((c, k) => c - P.get(k)));
                    const length = Math.hypot(d.get(0), d.get(1), d.get(2));
                    const Q = queryType === BVTree.SEGMENT_QUERY
                        ? Vector.fromArray(P.values.map((c, k) => c + 3 * d.get(k)))
                        : Vector.fromArray(d.values.map(c => c / length));

                    const result = tree.execute(queryType, P, Q);
                    const expected = bruteForceExecute(mesh, queryType, P, Q);
                    expect(result.intersections.map(h => h.triangleIndex))
                        .toEqual(expected.map(h => h.triangleIndex));
                    for (let i = 0; i < expected.length; ++i) {
                        expect(result.intersections[i].parameter)
                            .toBe(expected[i].parameter);
                    }
                    expect(tree.leafIndices(queryType, P, Q))
                        .toEqual(referenceLeaves(tree, queryType, P, Q));
                }
            }
        }
    });

    it('reports no leaves for a linear component that misses the root box', () => {
        const tree = buildTree(tetra);
        const P = v3(100, 100, 100);
        expect(tree.leafIndices(BVTree.LINE_QUERY, P, v3(1, 0, 0)).length).toBe(0);
        expect(tree.execute(BVTree.LINE_QUERY, P, v3(1, 0, 0)).intersections.length)
            .toBe(0);
    });
});
