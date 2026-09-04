import { describe, it, expect } from 'vitest';
import { BVTree, BVTreeNode } from '../src/BVTree.js';
import type {
    BVTreeBoundingVolume, BVTreeSplittingAxis, BVTreeVolumeOps
} from '../src/BVTree.js';
import {
    BVTreeOfTriangles, BVTreeOfTrianglesIntersection, intersectLineTriangle,
    intersectRayTriangle, intersectSegmentTriangle
} from '../src/BVTreeOfTriangles.js';
import { Triangle } from '../src/Triangle.js';
import {
    Vector, add, div, dot, length, mul, normalize, sub
} from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, fc, scaled
} from './helpers/arbitraries.js';

// ---------------------------------------------------------------------------
// A concrete BoundingVolume: an axis-aligned bounding box, the same minimal
// harness used by BVTree.test.ts.
// ---------------------------------------------------------------------------

class AABB implements BVTreeBoundingVolume {
    min: number[] = [Infinity, Infinity, Infinity];
    max: number[] = [-Infinity, -Infinity, -Infinity];

    reset(): void {
        this.min = [Infinity, Infinity, Infinity];
        this.max = [-Infinity, -Infinity, -Infinity];
    }

    grow(p: Vector): void {
        for (let j = 0; j < 3; ++j) {
            if (p.values[j] < this.min[j]) {
                this.min[j] = p.values[j];
            }
            if (p.values[j] > this.max[j]) {
                this.max[j] = p.values[j];
            }
        }
    }

    contains(p: Vector, eps: number = 1e-12): boolean {
        for (let j = 0; j < 3; ++j) {
            if (p.values[j] < this.min[j] - eps || p.values[j] > this.max[j] + eps) {
                return false;
            }
        }
        return true;
    }

    getSplittingAxis(): BVTreeSplittingAxis {
        const origin = Vector.fromArray([
            0.5 * (this.min[0] + this.max[0]),
            0.5 * (this.min[1] + this.max[1]),
            0.5 * (this.min[2] + this.max[2])
        ]);
        let jmax = 0;
        let emax = this.max[0] - this.min[0];
        for (let j = 1; j < 3; ++j) {
            const e = this.max[j] - this.min[j];
            if (e > emax) {
                emax = e;
                jmax = j;
            }
        }
        return { origin: origin, direction: Vector.unit(3, jmax) };
    }
}

function intersectSlabs(P: Vector, D: Vector, tmin0: number, tmax0: number,
    box: AABB): boolean {
    let tmin = tmin0;
    let tmax = tmax0;
    for (let j = 0; j < 3; ++j) {
        if (Math.abs(D.values[j]) < 1e-15) {
            if (P.values[j] < box.min[j] || P.values[j] > box.max[j]) {
                return false;
            }
        } else {
            let t0 = (box.min[j] - P.values[j]) / D.values[j];
            let t1 = (box.max[j] - P.values[j]) / D.values[j];
            if (t0 > t1) {
                const t = t0;
                t0 = t1;
                t1 = t;
            }
            tmin = Math.max(tmin, t0);
            tmax = Math.min(tmax, t1);
            if (tmin > tmax) {
                return false;
            }
        }
    }
    return true;
}

const aabbOps: BVTreeVolumeOps<AABB> = {
    create: () => new AABB(),
    intersectLine: (P, Q, bv) => intersectSlabs(P, Q, -Infinity, Infinity, bv),
    intersectRay: (P, Q, bv) => intersectSlabs(P, Q, 0, Infinity, bv),
    intersectSegment: (P, Q, bv) => intersectSlabs(P, sub(Q, P), 0, 1, bv)
};

class AABBTreeOfTriangles extends BVTreeOfTriangles<AABB> {
    constructor() {
        super(aabbOps);
    }

    private growTriangle(bv: AABB, triangleIndex: number): void {
        const tri = this.mTriangles[triangleIndex];
        for (let k = 0; k < 3; ++k) {
            bv.grow(this.mVertices[tri[k]]);
        }
    }

    protected computeInteriorBoundingVolume(i0: number, i1: number, bv: AABB): void {
        bv.reset();
        for (let i = i0; i <= i1; ++i) {
            this.growTriangle(bv, this.mPartition[i]);
        }
    }

    protected computeLeafBoundingVolume(i: number, bv: AABB): void {
        bv.reset();
        this.growTriangle(bv, this.mPartition[i]);
    }
}

function isLeaf(node: BVTreeNode<AABB>): boolean {
    return node.leftChild === BVTreeNode.invalid
        && node.rightChild === BVTreeNode.invalid;
}

// The axis-aligned box [0,2]x[0,3]x[0,5] as 12 triangles.
function makeBox(): {
    vertices: Vector[];
    triangles: [number, number, number][];
} {
    const vertices = [
        Vector.fromArray([0, 0, 0]),
        Vector.fromArray([2, 0, 0]),
        Vector.fromArray([2, 3, 0]),
        Vector.fromArray([0, 3, 0]),
        Vector.fromArray([0, 0, 5]),
        Vector.fromArray([2, 0, 5]),
        Vector.fromArray([2, 3, 5]),
        Vector.fromArray([0, 3, 5])
    ];
    const triangles: [number, number, number][] = [
        [0, 2, 1], [0, 3, 2],   // z = 0
        [4, 5, 6], [4, 6, 7],   // z = 5
        [0, 1, 5], [0, 5, 4],   // y = 0
        [3, 7, 6], [3, 6, 2],   // y = 3
        [0, 4, 7], [0, 7, 3],   // x = 0
        [1, 2, 6], [1, 6, 5]    // x = 2
    ];
    return { vertices: vertices, triangles: triangles };
}

// A brute-force replacement for execute() that visits every triangle and
// keeps every hit, ordered by (parameter, triangleIndex) -- the port's
// corrected container semantics (upstream #167 is fixed in the port).
function bruteForce(queryType: number, P: Vector, Q: Vector,
    vertices: Vector[], triangles: [number, number, number][]):
    BVTreeOfTrianglesIntersection[] {
    const query = [intersectLineTriangle, intersectRayTriangle,
        intersectSegmentTriangle][queryType];
    const intersections: BVTreeOfTrianglesIntersection[] = [];
    for (let t = 0; t < triangles.length; ++t) {
        const tri = triangles[t];
        const triangle = Triangle.fromVertices(vertices[tri[0]], vertices[tri[1]],
            vertices[tri[2]]);
        const output = query(P, Q, triangle);
        if (output.intersect) {
            intersections.push(new BVTreeOfTrianglesIntersection(t, output.point,
                output.parameter));
        }
    }
    intersections.sort((a, b) => a.parameter - b.parameter ||
        a.triangleIndex - b.triangleIndex);
    return intersections;
}

describe('BVTreeOfTriangles', () => {
    it('computes triangle centroids and copies its inputs', () => {
        const tree = new AABBTreeOfTriangles();
        const { vertices, triangles } = makeBox();
        tree.createFromTriangles(vertices, triangles);

        expect(tree.getCentroids().length).toBe(12);
        for (let t = 0; t < triangles.length; ++t) {
            const tri = triangles[t];
            for (let j = 0; j < 3; ++j) {
                const expected = (vertices[tri[0]].values[j] + vertices[tri[1]].values[j]
                    + vertices[tri[2]].values[j]) / 3;
                expect(tree.getCentroids()[t].values[j]).toBeCloseTo(expected, 15);
            }
        }

        // 12 primitives round up to 16 = 2^4 leaves, so the height is 4.
        expect(tree.getHeight()).toBe(4);
        expect(tree.getNodes().length).toBe(31);

        vertices[0].values[0] = 100;
        triangles[0][0] = 5;
        expect(tree.getVertices()[0].values[0]).toBe(0);
        expect(tree.getTriangles()[0][0]).toBe(0);
    });

    it('every node bounding volume contains its triangles', () => {
        const tree = new AABBTreeOfTriangles();
        const { vertices, triangles } = makeBox();
        tree.createFromTriangles(vertices, triangles);

        const partition = tree.getPartition();
        expect([...partition].sort((a, b) => a - b)).toEqual(
            triangles.map((_t, i) => i));

        for (const node of tree.getNodes()) {
            if (node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                const tri = triangles[partition[i]];
                for (let k = 0; k < 3; ++k) {
                    expect(node.boundingVolume.contains(vertices[tri[k]])).toBe(true);
                }
            }
        }
    });

    it('a ray through the box hits the two expected faces', () => {
        const tree = new AABBTreeOfTriangles();
        const { vertices, triangles } = makeBox();
        tree.createFromTriangles(vertices, triangles);

        // A ray from below the box, along +z through (1, 1, *).
        const P = Vector.fromArray([1, 1, -1]);
        const Q = Vector.fromArray([0, 0, 1]);
        const result = tree.execute(BVTree.RAY_QUERY, P, Q);

        expect(result.intersections.length).toBe(2);
        expect(result.intersections[0].parameter).toBeCloseTo(1, 12);
        expect(result.intersections[1].parameter).toBeCloseTo(6, 12);
        expect(result.intersections[0].point.values[0]).toBeCloseTo(1, 12);
        expect(result.intersections[0].point.values[1]).toBeCloseTo(1, 12);
        expect(result.intersections[0].point.values[2]).toBeCloseTo(0, 12);
        expect(result.intersections[1].point.values[2]).toBeCloseTo(5, 12);
        // The entry face is on z = 0 (triangles 0 and 1) and the exit face is
        // on z = 5 (triangles 2 and 3).
        expect([0, 1]).toContain(result.intersections[0].triangleIndex);
        expect([2, 3]).toContain(result.intersections[1].triangleIndex);

        // The reported node indices are distinct leaves.
        expect(new Set(result.nodeIndices).size).toBe(result.nodeIndices.length);
        for (const index of result.nodeIndices) {
            expect(isLeaf(tree.getNodes()[index])).toBe(true);
        }
    });

    it('matches brute force for random rays, lines and segments', () => {
        const tree = new AABBTreeOfTriangles();
        const { vertices, triangles } = makeBox();
        tree.createFromTriangles(vertices, triangles);

        let seed = 4242;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 200; ++trial) {
            // A target strictly inside the box and an origin outside it.
            const target = Vector.fromArray([
                0.1 + 1.8 * rand(),
                0.1 + 2.8 * rand(),
                0.1 + 4.8 * rand()
            ]);
            const origin = Vector.fromArray([
                -3 - 4 * rand(),
                -3 - 4 * rand(),
                -3 - 4 * rand()
            ]);
            const direction = sub(target, origin);

            const rayResult = tree.execute(BVTree.RAY_QUERY, origin, direction);
            const rayExpected = bruteForce(BVTree.RAY_QUERY, origin, direction,
                vertices, triangles);
            expect(rayResult.intersections.map((h) => h.triangleIndex)).toEqual(
                rayExpected.map((h) => h.triangleIndex));
            expect(rayResult.intersections.map((h) => h.parameter)).toEqual(
                rayExpected.map((h) => h.parameter));
            // A ray from outside toward an interior point enters and exits.
            expect(rayResult.intersections.length).toBe(2);

            const lineResult = tree.execute(BVTree.LINE_QUERY, origin, direction);
            const lineExpected = bruteForce(BVTree.LINE_QUERY, origin, direction,
                vertices, triangles);
            expect(lineResult.intersections.map((h) => h.triangleIndex)).toEqual(
                lineExpected.map((h) => h.triangleIndex));

            // A segment from the origin to the interior target: one crossing.
            const segResult = tree.execute(BVTree.SEGMENT_QUERY, origin, target);
            const segExpected = bruteForce(BVTree.SEGMENT_QUERY, origin, target,
                vertices, triangles);
            expect(segResult.intersections.map((h) => h.triangleIndex)).toEqual(
                segExpected.map((h) => h.triangleIndex));
            expect(segResult.intersections.length).toBe(1);

            // The intersection points lie on the linear component and on the
            // reported triangle's plane.
            for (const hit of rayResult.intersections) {
                for (let j = 0; j < 3; ++j) {
                    expect(hit.point.values[j]).toBeCloseTo(
                        origin.values[j] + hit.parameter * direction.values[j], 10);
                }
            }
        }
    });

    it('is sorted by (parameter, triangleIndex) and keeps equal-parameter intersections', () => {
        // Upstream stores the intersections in a std::set ordered by the
        // parameter alone, so two triangles hit at the same parameter keep
        // only a single entry (upstream #167). The port fixes this: ties are
        // broken by triangleIndex and every hit is kept. The ray runs along
        // the shared diagonals of both the z = 0 and z = 5 faces.
        const tree = new AABBTreeOfTriangles();
        const { vertices, triangles } = makeBox();
        tree.createFromTriangles(vertices, triangles);

        const P = Vector.fromArray([1, 1.5, -1]);
        const Q = Vector.fromArray([0, 0, 1]);
        const result = tree.execute(BVTree.RAY_QUERY, P, Q);

        // Triangles 0 and 1 are both hit at parameter 1 and triangles 2 and
        // 3 at parameter 6; all four hits are kept, ties ordered by
        // triangleIndex.
        expect(result.intersections.length).toBe(4);
        expect(result.intersections.map((x) => x.triangleIndex))
            .toEqual([0, 1, 2, 3]);
        expect(result.intersections[0].parameter).toBeCloseTo(1, 12);
        expect(result.intersections[1].parameter).toBeCloseTo(1, 12);
        expect(result.intersections[2].parameter).toBeCloseTo(6, 12);
        expect(result.intersections[3].parameter).toBeCloseTo(6, 12);

        // The parameters are non-decreasing.
        for (let i = 1; i < result.intersections.length; ++i) {
            expect(result.intersections[i - 1].parameter)
                .toBeLessThanOrEqual(result.intersections[i].parameter);
        }
    });

    it('reports nothing for a linear component that misses the mesh', () => {
        const tree = new AABBTreeOfTriangles();
        const { vertices, triangles } = makeBox();
        tree.createFromTriangles(vertices, triangles);

        const result = tree.execute(BVTree.LINE_QUERY,
            Vector.fromArray([100, 100, 100]), Vector.fromArray([0, 0, 1]));
        expect(result.nodeIndices.length).toBe(0);
        expect(result.intersections.length).toBe(0);
    });

    it('handles a single triangle (height 0)', () => {
        const tree = new AABBTreeOfTriangles();
        const vertices = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([0, 1, 0])
        ];
        tree.createFromTriangles(vertices, [[0, 1, 2]]);

        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        expect(isLeaf(tree.getNodes()[0])).toBe(true);

        const result = tree.execute(BVTree.RAY_QUERY,
            Vector.fromArray([0.25, 0.25, -1]), Vector.fromArray([0, 0, 1]));
        expect(result.nodeIndices).toEqual([0]);
        expect(result.intersections.length).toBe(1);
        expect(result.intersections[0].triangleIndex).toBe(0);
        expect(result.intersections[0].parameter).toBeCloseTo(1, 12);
    });

    it('honors a user-specified height', () => {
        const tree = new AABBTreeOfTriangles();
        const { vertices, triangles } = makeBox();
        tree.createFromTriangles(vertices, triangles, 1);

        expect(tree.getHeight()).toBe(1);
        expect(tree.getNodes().length).toBe(3);
        // The children are leaves representing 6 triangles each; execute
        // still finds the correct intersections.
        const P = Vector.fromArray([1, 1, -1]);
        const Q = Vector.fromArray([0, 0, 1]);
        const result = tree.execute(BVTree.RAY_QUERY, P, Q);
        expect(result.intersections.length).toBe(2);
        expect(result.intersections[0].parameter).toBeCloseTo(1, 12);
        expect(result.intersections[1].parameter).toBeCloseTo(6, 12);
    });

    it('the default Intersection is invalid with zero parameter', () => {
        const item = new BVTreeOfTrianglesIntersection();
        expect(item.triangleIndex).toBe(BVTreeOfTrianglesIntersection.invalid);
        expect(item.parameter).toBe(0);
        expect(item.point.values).toEqual([0, 0, 0]);

        const other = new BVTreeOfTrianglesIntersection(3, Vector.zero(3), 1);
        expect(item.lessThan(other)).toBe(true);
        expect(other.lessThan(item)).toBe(false);
    });

    it('throws for invalid inputs', () => {
        const tree = new AABBTreeOfTriangles();
        expect(() => tree.createFromTriangles([], [])).toThrow(
            'Expecting at least 3 vertices and at least 1 triangle.');
        expect(() => tree.createFromTriangles([
            Vector.zero(3), Vector.zero(3), Vector.zero(3)
        ], [])).toThrow('Expecting at least 3 vertices and at least 1 triangle.');
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V07): property-based re-check of BVTreeOfTriangles.h.
//
// The cross-checks below compare the pruned traversal against a brute-force
// sweep of every triangle, and the two must agree exactly. That is only a
// meaningful comparison when the per-triangle query is well conditioned: the
// FIQuery divides the barycentric numerators by Dot(D, N), so a direction
// whose component along the triangle normal is subnormal produces infinite
// parameters and hits that no bounding volume can contain. Coordinates are
// therefore drawn from a uniform grid (scaled(), no dynamic range at all)
// instead of fc.double(), which samples bit patterns and reaches 1e-320
// readily. Exactly-parallel configurations are still generated and are
// harmless: both sides run the same query and both report no intersection.
// ---------------------------------------------------------------------------

const gridCoordinate = scaled(-6, 6, 512);

const gridDirection = fc.tuple(scaled(-1, 1, 64), scaled(-1, 1, 64),
    scaled(-1, 1, 64))
    .map(c => Vector.fromArray([c[0], c[1], c[2]]))
    .filter(v => dot(v, v) > 0.05)
    .map(v => { const u = v.clone(); normalize(u); return u; });

const gridPoint = fc.tuple(gridCoordinate, gridCoordinate, gridCoordinate)
    .map(c => Vector.fromArray([c[0], c[1], c[2]]));

// A triangle soup: independent, non-degenerate triangles with their own
// vertices, so no shared edges muddy the brute-force comparison.
const triangleSoup = fc.array(
    fc.tuple(gridPoint, gridPoint, gridPoint)
        .filter((t) => {
            const n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
            return dot(n, n) > 1;
        }),
    { minLength: 1, maxLength: 9 })
    .map((tris) => {
        const vertices: Vector[] = [];
        const triangles: [number, number, number][] = [];
        for (const t of tris) {
            const b = vertices.length;
            vertices.push(t[0], t[1], t[2]);
            triangles.push([b, b + 1, b + 2]);
        }
        return { vertices: vertices, triangles: triangles };
    });

function bvtReachable(tree: AABBTreeOfTriangles): number[] {
    const nodes = tree.getNodes();
    const list: number[] = [];
    const stack = [0];
    while (stack.length > 0) {
        const index = stack.pop()!;
        list.push(index);
        const node = nodes[index];
        if (!isLeaf(node)) {
            stack.push(node.leftChild, node.rightChild);
        }
    }
    return list;
}

// The point lies in the plane of the triangle and has nonnegative barycentric
// coordinates (it is inside or on the solid triangle).
function expectOnTriangle(point: Vector, tri: [number, number, number],
    vertices: Vector[]): void {
    const v0 = vertices[tri[0]];
    const e1 = sub(vertices[tri[1]], v0);
    const e2 = sub(vertices[tri[2]], v0);
    const n = cross(e1, e2);
    const d = sub(point, v0);
    const nn = dot(n, n);
    // Plane membership, scaled by the triangle size so the tolerance is
    // relative to the data rather than absolute.
    expect(Math.abs(dot(n, d))).toBeLessThanOrEqual(1e-8 * nn);
    const b1 = dot(cross(d, e2), n) / nn;
    const b2 = dot(cross(e1, d), n) / nn;
    expect(b1).toBeGreaterThanOrEqual(-1e-8);
    expect(b2).toBeGreaterThanOrEqual(-1e-8);
    expect(b1 + b2).toBeLessThanOrEqual(1 + 1e-8);
}

describe('BVTreeOfTriangles verification', () => {
    it('computes the exact triangle centroids and copies its inputs', () => {
        check(triangleSoup, (mesh) => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromTriangles(mesh.vertices, mesh.triangles);

            const centroids = tree.getCentroids();
            const third = 1 / 3;
            for (let t = 0; t < mesh.triangles.length; ++t) {
                const tri = mesh.triangles[t];
                for (let j = 0; j < 3; ++j) {
                    // Upstream is (V0 + V1 + V2) / three, and Vector's
                    // operator/= multiplies by the reciprocal rather than
                    // dividing, so the port must too. (v * (1/3) and v / 3
                    // differ by an ulp on many inputs, which is what makes
                    // this an exact rather than an approximate check.)
                    expect(centroids[t].values[j]).toBe(
                        (mesh.vertices[tri[0]].values[j]
                            + mesh.vertices[tri[1]].values[j]
                            + mesh.vertices[tri[2]].values[j]) * third);
                }
            }

            expect(tree.getTriangles().map(t => [t[0], t[1], t[2]]))
                .toEqual(mesh.triangles.map(t => [t[0], t[1], t[2]]));
            for (let i = 0; i < mesh.vertices.length; ++i) {
                expect(tree.getVertices()[i]).not.toBe(mesh.vertices[i]);
            }
        }, 100);
    });

    it('matches brute force over every triangle for all three query types', () => {
        check(fc.tuple(triangleSoup, gridPoint, gridDirection, gridPoint),
            (input) => {
                const mesh = input[0];
                const P = input[1];
                const D = input[2];
                const Q = input[3];
                const tree = new AABBTreeOfTriangles();
                tree.createFromTriangles(mesh.vertices, mesh.triangles);

                const cases: Array<{ queryType: number; A: Vector; B: Vector }> = [
                    { queryType: BVTree.LINE_QUERY, A: P, B: D },
                    { queryType: BVTree.RAY_QUERY, A: P, B: D },
                    { queryType: BVTree.SEGMENT_QUERY, A: P, B: Q }
                ];

                for (const c of cases) {
                    const result = tree.execute(c.queryType, c.A, c.B);
                    const expected = bruteForce(c.queryType, c.A, c.B,
                        mesh.vertices, mesh.triangles);
                    // The tree runs the identical per-triangle query, so the
                    // surviving hits agree bit for bit; only the pruning can
                    // differ, and a correct bounding volume prunes nothing
                    // that the triangle query would have hit.
                    expect(result.intersections.map(x => x.triangleIndex))
                        .toEqual(expected.map(x => x.triangleIndex));
                    for (let k = 0; k < expected.length; ++k) {
                        expect(result.intersections[k].parameter)
                            .toBe(expected[k].parameter);
                        expect([...result.intersections[k].point.values])
                            .toEqual([...expected[k].point.values]);
                    }
                    // The reported leaves are distinct leaves of the tree.
                    const nodes = tree.getNodes();
                    expect(new Set(result.nodeIndices).size)
                        .toBe(result.nodeIndices.length);
                    for (const index of result.nodeIndices) {
                        expect(isLeaf(nodes[index])).toBe(true);
                    }
                    expect(bvtReachable(tree)).toEqual(
                        expect.arrayContaining(result.nodeIndices));
                }
            }, 60);
    });

    it('sorts the hits by (parameter, triangleIndex) and keeps them all', () => {
        check(fc.tuple(triangleSoup, gridPoint, gridDirection), (input) => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromTriangles(input[0].vertices, input[0].triangles);
            const hits = tree.execute(BVTree.LINE_QUERY, input[1],
                input[2]).intersections;
            for (let k = 1; k < hits.length; ++k) {
                expect(hits[k - 1].lessThan(hits[k])).toBe(true);
                expect(hits[k].lessThan(hits[k - 1])).toBe(false);
            }
            // No triangle is reported twice: the leaf ranges are disjoint.
            expect(new Set(hits.map(h => h.triangleIndex)).size)
                .toBe(hits.length);
        }, 100);
    });

    it('places each hit on its triangle and on the linear component', () => {
        check(fc.tuple(triangleSoup, gridPoint, gridDirection, gridPoint),
            (input) => {
                const mesh = input[0];
                const P = input[1];
                const D = input[2];
                const Q = input[3];
                const tree = new AABBTreeOfTriangles();
                tree.createFromTriangles(mesh.vertices, mesh.triangles);

                // Line and ray: the reported parameter is the parameter of
                // P + t * Q directly.
                for (const queryType of [BVTree.LINE_QUERY, BVTree.RAY_QUERY]) {
                    for (const hit of tree.execute(queryType, P, D).intersections) {
                        if (queryType === BVTree.RAY_QUERY) {
                            expect(hit.parameter).toBeGreaterThanOrEqual(0);
                        }
                        for (let j = 0; j < 3; ++j) {
                            expectClose(hit.point.values[j],
                                P.values[j] + hit.parameter * D.values[j],
                                1e-9, 1e-9);
                        }
                        expectOnTriangle(hit.point,
                            mesh.triangles[hit.triangleIndex], mesh.vertices);
                    }
                }

                // Segment: upstream BVTreeOfTriangles.h forwards the FIQuery
                // parameter unchanged, and that is the parameter s of the
                // centered form C + s * D with |s| <= e -- NOT the t in [0,1]
                // of (1-t)*P + t*Q that BVTree.h documents. OBBTreeOfTriangles.h
                // does apply the conversion; see the port notes. The port
                // preserves the upstream value, which this property pins.
                const diff = sub(Q, P);
                const len = length(diff);
                if (len > 1e-3) {
                    const dir = div(diff, len);
                    const mid = mul(0.5, add(P, Q));
                    for (const hit of tree.execute(BVTree.SEGMENT_QUERY, P, Q)
                        .intersections) {
                        expect(Math.abs(hit.parameter))
                            .toBeLessThanOrEqual(0.5 * len + 1e-9);
                        for (let j = 0; j < 3; ++j) {
                            expectClose(hit.point.values[j],
                                mid.values[j] + hit.parameter * dir.values[j],
                                1e-9, 1e-9);
                        }
                        expectOnTriangle(hit.point,
                            mesh.triangles[hit.triangleIndex], mesh.vertices);
                    }
                }
            }, 50);
    });

    it('bounds every node by its own triangles and tiles the partition', () => {
        check(triangleSoup, (mesh) => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromTriangles(mesh.vertices, mesh.triangles);
            const nodes = tree.getNodes();
            const partition = tree.getPartition();

            expect([...partition].sort((a, b) => a - b))
                .toEqual([...Array(mesh.triangles.length).keys()]);

            const covered = new Array<number>(mesh.triangles.length).fill(0);
            for (const index of bvtReachable(tree)) {
                const node = nodes[index];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const tri = mesh.triangles[partition[i]];
                    for (let k = 0; k < 3; ++k) {
                        expect(node.boundingVolume.contains(mesh.vertices[tri[k]]))
                            .toBe(true);
                    }
                }
                if (isLeaf(node)) {
                    for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                        ++covered[i];
                    }
                }
            }
            expect(covered.every(c => c === 1)).toBe(true);
        }, 100);
    });
});
