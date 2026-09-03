import { describe, it, expect } from 'vitest';
import { OrientedBoxTreeOfTriangles } from '../src/OrientedBoxTreeOfTriangles';
import { OrientedBoxBV } from '../src/OrientedBoxBV';
import { BVTree, BVTreeNode } from '../src/BVTree';
import {
    intersectLineTriangle, intersectRayTriangle, intersectSegmentTriangle
} from '../src/BVTreeOfTriangles';
import type { LinearTriangleResult } from '../src/BVTreeOfTriangles';
import { Triangle } from '../src/Triangle';
import { Vector, dot, sub } from '../src/Vector';

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

function boxContains(bv: OrientedBoxBV, point: Vector, epsilon: number): boolean {
    const box = bv.box;
    const diff = sub(point, box.center);
    for (let k = 0; k < 3; ++k) {
        const y = dot(diff, box.axis[k]);
        if (Math.abs(y) > box.extent.get(k) + epsilon) {
            return false;
        }
    }
    return true;
}

function expectValidBox(bv: OrientedBoxBV): void {
    const box = bv.box;
    for (let i = 0; i < 3; ++i) {
        expect(box.extent.get(i)).toBeGreaterThanOrEqual(0);
        expect(dot(box.axis[i], box.axis[i])).toBeCloseTo(1, 10);
        for (let j = i + 1; j < 3; ++j) {
            expect(dot(box.axis[i], box.axis[j])).toBeCloseTo(0, 10);
        }
    }
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

class TestableTree extends OrientedBoxTreeOfTriangles {
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

function referenceLeaves(tree: TestableTree, queryType: number, P: Vector,
    Q: Vector): number[] {
    const nodes = tree.getNodes();
    const hits = (bv: OrientedBoxBV): boolean =>
        queryType === BVTree.LINE_QUERY ? OrientedBoxBV.intersectLine(P, Q, bv)
            : (queryType === BVTree.RAY_QUERY
                ? OrientedBoxBV.intersectRay(P, Q, bv)
                : OrientedBoxBV.intersectSegment(P, Q, bv));

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

describe('OrientedBoxTreeOfTriangles leaf bounding volumes', () => {
    it('makes a leaf box tight around its triangle with one zero extent', () => {
        // A right triangle in the plane z = 0.
        const mesh: Mesh = {
            vertices: [v3(0, 0, 0), v3(4, 0, 0), v3(0, 2, 0)],
            triangles: [[0, 1, 2]]
        };
        const tree = buildTree(mesh);
        expect(tree.getNodes().length).toBe(1);
        const bv = tree.getNodes()[0].boundingVolume;

        // The zeroed extent is the one perpendicular to the plane of the
        // triangle.
        const zeroIndices = bv.box.extent.values
            .map((e, i) => ({ e: e, i: i })).filter(x => x.e === 0)
            .map(x => x.i);
        expect(zeroIndices.length).toBe(1);
        expect(Math.abs(bv.box.axis[zeroIndices[0]].get(2))).toBeCloseTo(1, 10);

        for (const j of mesh.triangles[0]) {
            expect(boxContains(bv, mesh.vertices[j], 1e-9)).toBe(true);
        }
    });

    it('zeroes the SMALLEST extent, so the leaf box still holds the triangle', () => {
        // Regression for the upstream bug in OrientedBoxTreeOfTriangles.h,
        // where the last comparison of the "find the minimum" scan is
        // 'absExtent > minAbsExtent'. With '>' the largest extent is zeroed
        // whenever extent[2] is not the smallest, collapsing the leaf box
        // along its longest axis so that it no longer contains its triangle.
        const rand = makeRandom(20260902);
        const meshes: Mesh[] = [
            tetra,
            makeBoxMesh(2, v3(-1, -2, -3), v3(4, 1, 2)),
            {
                vertices: Array.from({ length: 30 },
                    () => v3(5 * rand() - 2, 4 * rand(), 3 * rand() - 1)),
                triangles: Array.from({ length: 10 },
                    (_, k): Tri => [3 * k, 3 * k + 1, 3 * k + 2])
            }
        ];
        for (const mesh of meshes) {
            const tree = buildTree(mesh);
            const partition = tree.getPartition();
            let leafCount = 0;
            for (const node of tree.getNodes()) {
                if (node.leftChild !== BVTreeNode.invalid ||
                    node.minIndex === BVTreeNode.invalid) {
                    continue;
                }
                ++leafCount;
                expect(node.minIndex).toBe(node.maxIndex);
                expectValidBox(node.boundingVolume);

                const extents = node.boundingVolume.box.extent.values;
                expect(extents.filter(e => e === 0).length)
                    .toBeGreaterThanOrEqual(1);
                // The zeroed extent was the smallest of the three; every
                // other extent is at least as large as the zeroed one.
                const maxExtent = Math.max(...extents);
                expect(maxExtent).toBeGreaterThan(0);

                const tri = mesh.triangles[partition[node.minIndex]];
                for (const j of tri) {
                    expect(boxContains(node.boundingVolume, mesh.vertices[j],
                        1e-8)).toBe(true);
                }
            }
            expect(leafCount).toBe(mesh.triangles.length);
        }
    });

    it('handles a degenerate (zero-area) triangle', () => {
        const mesh: Mesh = {
            vertices: [v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2), v3(0, 1, 0)],
            triangles: [[0, 1, 2], [0, 1, 3]]
        };
        const tree = buildTree(mesh);
        const partition = tree.getPartition();
        for (const node of tree.getNodes()) {
            if (node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                for (const j of mesh.triangles[partition[i]]) {
                    expect(boxContains(node.boundingVolume, mesh.vertices[j],
                        1e-8)).toBe(true);
                }
            }
        }
    });
});

describe('OrientedBoxTreeOfTriangles interior bounding volumes', () => {
    it('contains every vertex of its range in every node box', () => {
        const meshes: Mesh[] = [
            tetra,
            makeBoxMesh(1, v3(-1, -2, -3), v3(4, 1, 2)),
            makeBoxMesh(2, v3(0, 0, 0), v3(1, 1, 1)),
            makeBoxMesh(3, v3(-2, -1, 0), v3(2, 5, 1))
        ];
        for (const mesh of meshes) {
            const tree = buildTree(mesh);
            const partition = tree.getPartition();
            let interior = 0;
            for (const node of tree.getNodes()) {
                if (node.minIndex === BVTreeNode.invalid) {
                    continue;
                }
                expectValidBox(node.boundingVolume);
                if (node.maxIndex > node.minIndex) {
                    ++interior;
                }
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    for (const j of mesh.triangles[partition[i]]) {
                        expect(boxContains(node.boundingVolume,
                            mesh.vertices[j], 1e-8)).toBe(true);
                    }
                }
            }
            expect(interior).toBeGreaterThan(0);
        }
    });

    it('contains the child ranges in the parent box', () => {
        const mesh = makeBoxMesh(2, v3(0, 0, 0), v3(2, 3, 5));
        const tree = buildTree(mesh);
        const nodes = tree.getNodes();
        const partition = tree.getPartition();
        let interiorCount = 0;
        for (const node of nodes) {
            if (node.leftChild === BVTreeNode.invalid) {
                continue;
            }
            ++interiorCount;
            for (const c of [node.leftChild, node.rightChild]) {
                const child = nodes[c];
                for (let i = child.minIndex; i <= child.maxIndex; ++i) {
                    for (const j of mesh.triangles[partition[i]]) {
                        expect(boxContains(node.boundingVolume,
                            mesh.vertices[j], 1e-8)).toBe(true);
                    }
                }
            }
        }
        expect(interiorCount).toBeGreaterThan(0);
    });

    it('respects a user-specified height, leaving multi-triangle leaves', () => {
        const mesh = makeBoxMesh(2, v3(0, 0, 0), v3(1, 1, 1));
        const shallow = buildTree(mesh, 2);
        expect(shallow.getHeight()).toBe(2);
        expect(buildTree(mesh).getHeight()).toBeGreaterThan(2);

        const partition = shallow.getPartition();
        let multiLeaf = 0;
        for (const node of shallow.getNodes()) {
            if (node.leftChild !== BVTreeNode.invalid ||
                node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            if (node.maxIndex > node.minIndex) {
                ++multiLeaf;
            }
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                for (const j of mesh.triangles[partition[i]]) {
                    expect(boxContains(node.boundingVolume, mesh.vertices[j],
                        1e-8)).toBe(true);
                }
            }
        }
        expect(multiLeaf).toBeGreaterThan(0);
    });

    it('does not alias the input vertices in the node boxes', () => {
        const vertices = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const tree = buildTree({ vertices: vertices, triangles: [[0, 1, 2]] });
        tree.getNodes()[0].boundingVolume.box.center.set(0, -100);
        expect(vertices[0].get(0)).toBe(0);
        expect(tree.getVertices()[0].get(0)).toBe(0);
    });

    it('rejects degenerate creation input', () => {
        const tree = new OrientedBoxTreeOfTriangles();
        expect(() => tree.createFromTriangles([v3(0, 0, 0), v3(1, 0, 0)],
            [[0, 1, 1]])).toThrow();
        expect(() => tree.createFromTriangles(tetra.vertices, [])).toThrow();
    });
});

describe('OrientedBoxTreeOfTriangles queries', () => {
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
        expect(result.intersections.map(h => h.triangleIndex))
            .toEqual(bruteForceExecute(mesh, BVTree.LINE_QUERY, P, Q)
                .map(h => h.triangleIndex));

        // A ray starting inside the box hits only the exit face.
        const rayResult = tree.execute(BVTree.RAY_QUERY, v3(0.5, 0.3, 0.7), Q);
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
        const result = tree.execute(BVTree.LINE_QUERY, v3(0.5, 0.5, -1),
            v3(0, 0, 1));
        expect(result.intersections.length).toBe(2);
        expect(result.intersections[0].parameter)
            .toBe(result.intersections[1].parameter);
        expect(result.intersections.map(h => h.triangleIndex)).toEqual([0, 1]);
    });

    it('agrees with brute force on randomized queries', () => {
        const meshes: Mesh[] = [
            tetra,
            makeBoxMesh(2, v3(-1, -1, -1), v3(1, 1, 1))
        ];
        const rand = makeRandom(2024);
        for (const mesh of meshes) {
            const tree = buildTree(mesh);
            for (const queryType of [BVTree.LINE_QUERY, BVTree.RAY_QUERY,
                BVTree.SEGMENT_QUERY]) {
                for (let trial = 0; trial < 20; ++trial) {
                    const P = v3(6 * (2 * rand() - 1), 6 * (2 * rand() - 1),
                        6 * (2 * rand() - 1));
                    const target = v3(2 * rand() - 1, 2 * rand() - 1,
                        2 * rand() - 1);
                    const d = sub(target, P);
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
