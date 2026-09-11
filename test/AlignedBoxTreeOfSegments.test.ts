import { describe, it, expect } from 'vitest';
import { AlignedBoxTreeOfSegments } from '../src/AlignedBoxTreeOfSegments.js';
import { AlignedBoxBV } from '../src/AlignedBoxBV.js';
import { BVTree, BVTreeNode } from '../src/BVTree.js';
import { Vector, add, mul } from '../src/Vector.js';
import {
    check, fc, finite, latticeVector, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A deterministic pseudorandom generator, so failures reproduce.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

type Seg = [number, number];

interface Curve {
    vertices: Vector[];
    segments: Seg[];
}

// A closed polyline on a helix; consecutive segments share a vertex.
function makeHelix(n: number, radius: number, pitch: number): Curve {
    const vertices: Vector[] = [];
    const segments: Seg[] = [];
    for (let i = 0; i < n; ++i) {
        const angle = (2 * Math.PI * i) / n;
        vertices.push(v3(radius * Math.cos(angle), radius * Math.sin(angle),
            pitch * i));
    }
    for (let i = 0; i + 1 < n; ++i) {
        segments.push([i, i + 1]);
    }
    segments.push([n - 1, 0]);
    return { vertices: vertices, segments: segments };
}

function makeRandomCurve(seed: number, numVertices: number,
    numSegments: number, scale: number): Curve {
    const rand = makeRandom(seed);
    const vertices: Vector[] = [];
    for (let i = 0; i < numVertices; ++i) {
        vertices.push(v3(scale * (2 * rand() - 1), scale * (2 * rand() - 1),
            scale * (2 * rand() - 1)));
    }
    const segments: Seg[] = [];
    for (let i = 0; i < numSegments; ++i) {
        const i0 = Math.floor(rand() * numVertices);
        let i1 = Math.floor(rand() * numVertices);
        if (i1 === i0) {
            i1 = (i0 + 1) % numVertices;
        }
        segments.push([i0, i1]);
    }
    return { vertices: vertices, segments: segments };
}

// The component-wise min/max over the endpoints of a set of segments,
// computed independently of the port.
function tightBounds(curve: Curve, segmentIndices: readonly number[]):
    { min: number[], max: number[] } {
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];
    for (const s of segmentIndices) {
        for (const j of curve.segments[s]) {
            for (let k = 0; k < 3; ++k) {
                min[k] = Math.min(min[k], curve.vertices[j].get(k));
                max[k] = Math.max(max[k], curve.vertices[j].get(k));
            }
        }
    }
    return { min: min, max: max };
}

// Exposes the protected getLeafIndices.
class TestableTree extends AlignedBoxTreeOfSegments {
    leafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }
}

function buildTree(curve: Curve, height?: number): TestableTree {
    const tree = new TestableTree();
    tree.createFromSegments(curve.vertices, curve.segments,
        height === undefined ? BVTree.fullHeight : height);
    return tree;
}

// The segment indices a tree reports as candidates.
function candidates(tree: TestableTree, nodeIndices: readonly number[]): Set<number> {
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

describe('AlignedBoxTreeOfSegments construction', () => {
    it('builds a single-node tree for one segment with a tight leaf box', () => {
        const tree = buildTree({
            vertices: [v3(2, 5, -1), v3(-3, 1, 4)],
            segments: [[0, 1]]
        });
        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        const root = tree.getNodes()[0];
        expect(root.leftChild).toBe(BVTreeNode.invalid);
        expect(root.boundingVolume.box.min.values).toEqual([-3, 1, -1]);
        expect(root.boundingVolume.box.max.values).toEqual([2, 5, 4]);
    });

    it('gives the root the tight box of all the segment endpoints', () => {
        const curve = makeHelix(8, 2, 0.5);
        const tree = buildTree(curve);
        const all = curve.segments.map((_, i) => i);
        const expected = tightBounds(curve, all);
        const root = tree.getNodes()[0];
        expect(root.boundingVolume.box.min.values).toEqual(expected.min);
        expect(root.boundingVolume.box.max.values).toEqual(expected.max);
    });

    it('gives every visited node the tight box of its segment range', () => {
        const curves: Curve[] = [
            { vertices: [v3(0, 0, 0), v3(1, 0, 0), v3(1, 1, 0)],
                segments: [[0, 1], [1, 2]] },
            makeHelix(9, 3, 0.25),
            makeHelix(16, 1, 2),
            makeRandomCurve(11, 12, 17, 5)
        ];
        for (const curve of curves) {
            const tree = buildTree(curve);
            const partition = tree.getPartition();
            let leafCount = 0;
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
                const expected = tightBounds(curve, range);
                expect(node.boundingVolume.box.min.values).toEqual(expected.min);
                expect(node.boundingVolume.box.max.values).toEqual(expected.max);
                if (node.minIndex === node.maxIndex) {
                    ++leafCount;
                }
            }
            expect(leafCount).toBe(curve.segments.length);
            expect(visited).toBeGreaterThanOrEqual(curve.segments.length);
        }
    });

    it('nests child boxes inside parent boxes', () => {
        const curve = makeRandomCurve(22, 20, 31, 4);
        const tree = buildTree(curve);
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

    it('respects a user-specified height, leaving multi-segment leaves tight', () => {
        const curve = makeHelix(20, 2, 0.3);
        const full = buildTree(curve);
        const shallow = buildTree(curve, 2);
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
            const expected = tightBounds(curve, range);
            expect(node.boundingVolume.box.min.values).toEqual(expected.min);
            expect(node.boundingVolume.box.max.values).toEqual(expected.max);
        }
        expect(multiLeaf).toBeGreaterThan(0);
    });

    it('handles degenerate (zero-length) segments', () => {
        const curve: Curve = {
            vertices: [v3(1, 2, 3), v3(4, 5, 6)],
            segments: [[0, 0], [1, 1], [0, 1]]
        };
        const tree = buildTree(curve);
        const partition = tree.getPartition();
        for (const node of tree.getNodes()) {
            if (node.minIndex === BVTreeNode.invalid ||
                node.minIndex !== node.maxIndex) {
                continue;
            }
            const s = partition[node.minIndex];
            const expected = tightBounds(curve, [s]);
            expect(node.boundingVolume.box.min.values).toEqual(expected.min);
            expect(node.boundingVolume.box.max.values).toEqual(expected.max);
        }
    });

    it('does not alias the input vertices in the node boxes', () => {
        const vertices = [v3(0, 0, 0), v3(1, 1, 1)];
        const tree = buildTree({ vertices: vertices, segments: [[0, 1]] });
        tree.getNodes()[0].boundingVolume.box.min.set(0, -100);
        expect(vertices[0].get(0)).toBe(0);
        expect(tree.getVertices()[0].get(0)).toBe(0);
    });

    it('rejects an empty vertex set', () => {
        const tree = new AlignedBoxTreeOfSegments();
        expect(() => tree.createFromSegments([], [])).toThrow();
    });
});

describe('AlignedBoxTreeOfSegments queries', () => {
    it('matches a recursive reference traversal on random queries', () => {
        const curve = makeRandomCurve(33, 24, 37, 5);
        const tree = buildTree(curve);
        const rand = makeRandom(44);

        let totalReported = 0;
        for (const queryType of [BVTree.LINE_QUERY, BVTree.RAY_QUERY,
            BVTree.SEGMENT_QUERY]) {
            for (let trial = 0; trial < 20; ++trial) {
                const target = curve.vertices[
                    Math.floor(rand() * curve.vertices.length)];
                const P = v3(20 * (2 * rand() - 1), 20 * (2 * rand() - 1),
                    20 * (2 * rand() - 1));
                const dir = Vector.fromArray(
                    target.values.map((c, k) => c - P.get(k)));
                const length = Math.hypot(dir.get(0), dir.get(1), dir.get(2));
                const Q = queryType === BVTree.SEGMENT_QUERY
                    ? Vector.fromArray(P.values.map((c, k) => c + 2 * dir.get(k)))
                    : Vector.fromArray(dir.values.map(c => c / length));

                const leaves = tree.leafIndices(queryType, P, Q);
                expect(leaves).toEqual(referenceLeaves(tree, queryType, P, Q));
                totalReported += candidates(tree, leaves).size;
            }
        }
        expect(totalReported).toBeGreaterThan(0);
    });

    it('reports every segment whose box an axis-parallel query meets', () => {
        // The line x = v.x, y = v.y meets the box of every segment having v as
        // an endpoint, with no rounding in the coordinates. Every segment box
        // that the line meets must be reported, because the boxes of all the
        // ancestors of its leaf contain that box.
        const curve = makeRandomCurve(55, 14, 21, 3);
        const tree = buildTree(curve);
        for (const v of curve.vertices) {
            const P = v3(v.get(0), v.get(1), -100);
            const Q = v3(0, 0, 1);
            const reported = candidates(tree,
                tree.leafIndices(BVTree.LINE_QUERY, P, Q));
            let hitCount = 0;
            for (let s = 0; s < curve.segments.length; ++s) {
                const bounds = tightBounds(curve, [s]);
                const bv = new AlignedBoxBV();
                bv.box.min = Vector.fromArray(bounds.min);
                bv.box.max = Vector.fromArray(bounds.max);
                if (AlignedBoxBV.intersectLine(P, Q, bv)) {
                    ++hitCount;
                    expect(reported.has(s)).toBe(true);
                }
            }
            // Each vertex is an endpoint of at least one segment, whose box
            // therefore contains the line's (x, y).
            expect(hitCount).toBeGreaterThan(0);
        }
    });

    it('reports no leaves for a linear component that misses the root box', () => {
        const tree = buildTree(makeHelix(10, 1, 0.2));
        const P = v3(100, 100, 100);
        expect(tree.leafIndices(BVTree.LINE_QUERY, P, v3(1, 0, 0)).length).toBe(0);
        expect(tree.leafIndices(BVTree.RAY_QUERY, P, v3(1, 0, 0)).length).toBe(0);
        expect(tree.leafIndices(BVTree.SEGMENT_QUERY, P, v3(200, 100, 100)).length)
            .toBe(0);
    });

    it('execute() returns the same leaf indices as getLeafIndices', () => {
        const tree = buildTree(makeHelix(12, 2, 0.5));
        const P = v3(-10, 0, 1);
        const Q = v3(1, 0, 0);
        expect(tree.execute(BVTree.LINE_QUERY, P, Q))
            .toEqual(tree.leafIndices(BVTree.LINE_QUERY, P, Q));
    });
});

// ---------------------------------------------------------------------------
// V43 verification: the tree invariants (tight boxes over segment endpoints,
// nested boxes, partitioned ranges) and the traversal against a recursive one.
// ---------------------------------------------------------------------------

function reachableNodes(tree: TestableTree): number[] {
    const nodes = tree.getNodes();
    const out: number[] = [];
    const visit = (nodeIndex: number): void => {
        out.push(nodeIndex);
        const node = nodes[nodeIndex];
        if (node.leftChild !== BVTreeNode.invalid &&
            node.rightChild !== BVTreeNode.invalid) {
            visit(node.leftChild);
            visit(node.rightChild);
        }
    };
    visit(0);
    return out;
}

describe('AlignedBoxTreeOfSegments verification', () => {
    // Integer endpoints make the tight bounds exactly representable, so the
    // box comparisons need no tolerance.
    const curveArb = fc.array(latticeVector(3, -8, 8),
        { minLength: 2, maxLength: 16 }).chain(vertices =>
            fc.tuple(fc.constant(vertices),
                fc.array(fc.tuple(fc.nat(), fc.nat()),
                    { minLength: 1, maxLength: 16 })))
        .map(([vertices, raw]) => {
            const segments: Seg[] = raw.map(([a, b]) => {
                const i0 = a % vertices.length;
                let i1 = b % vertices.length;
                if (i1 === i0) {
                    i1 = (i0 + 1) % vertices.length;
                }
                return [i0, i1] as Seg;
            });
            return { vertices: vertices, segments: segments } as Curve;
        });

    const queryArb = fc.tuple(fc.integer({ min: 0, max: 2 }),
        wellScaledVector(3, -12, 12), unitVector(3), finite(1, 25));

    it('gives every reachable node the tight box of its range', () => {
        check(fc.tuple(curveArb, fc.option(fc.integer({ min: 0, max: 4 }),
            { nil: undefined })), ([curve, height]) => {
            const tree = buildTree(curve, height);
            const nodes = tree.getNodes();
            const partition = tree.getPartition();
            for (const nodeIndex of reachableNodes(tree)) {
                const node = nodes[nodeIndex];
                const range: number[] = [];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    range.push(partition[i]);
                }
                expect(range.length).toBeGreaterThan(0);
                const bounds = tightBounds(curve, range);
                expect([...node.boundingVolume.box.min.values])
                    .toEqual(bounds.min);
                expect([...node.boundingVolume.box.max.values])
                    .toEqual(bounds.max);
                // The box contains both endpoints of every segment in range.
                for (const s of range) {
                    for (const j of curve.segments[s]) {
                        for (let k = 0; k < 3; ++k) {
                            expect(curve.vertices[j].get(k))
                                .toBeGreaterThanOrEqual(
                                    node.boundingVolume.box.min.get(k));
                            expect(curve.vertices[j].get(k))
                                .toBeLessThanOrEqual(
                                    node.boundingVolume.box.max.get(k));
                        }
                    }
                }
            }
        }, 100);
    });

    it('partitions the segments across the children of every node', () => {
        check(curveArb, curve => {
            const tree = buildTree(curve);
            const nodes = tree.getNodes();
            for (const nodeIndex of reachableNodes(tree)) {
                const node = nodes[nodeIndex];
                if (node.leftChild === BVTreeNode.invalid) {
                    continue;
                }
                const left = nodes[node.leftChild];
                const right = nodes[node.rightChild];
                expect(left.minIndex).toBe(node.minIndex);
                expect(right.maxIndex).toBe(node.maxIndex);
                expect(right.minIndex).toBe(left.maxIndex + 1);
                for (const child of [left, right]) {
                    for (let k = 0; k < 3; ++k) {
                        expect(child.boundingVolume.box.min.get(k))
                            .toBeGreaterThanOrEqual(
                                node.boundingVolume.box.min.get(k));
                        expect(child.boundingVolume.box.max.get(k))
                            .toBeLessThanOrEqual(
                                node.boundingVolume.box.max.get(k));
                    }
                }
            }
            expect([...tree.getPartition()].sort((a, b) => a - b))
                .toEqual(curve.segments.map((_, i) => i));
        }, 100);
    });

    it('matches a recursive traversal and reports only leaves', () => {
        check(fc.tuple(curveArb, queryArb),
            ([curve, [queryType, P, D, len]]) => {
                const tree = buildTree(curve);
                const Q = (queryType === BVTree.SEGMENT_QUERY
                    ? add(P, mul(D, len)) : D);
                const reported = tree.execute(queryType, P, Q);
                expect(reported).toEqual(
                    referenceLeaves(tree, queryType, P, Q));
                expect(new Set(reported).size).toBe(reported.length);
                const nodes = tree.getNodes();
                for (const nodeIndex of reported) {
                    expect(nodes[nodeIndex].leftChild)
                        .toBe(BVTreeNode.invalid);
                }
                // Every leaf whose own box the linear component meets is
                // reported (GetLeafIndices never tests a leaf's own volume,
                // so the reported set is a superset).
                for (const nodeIndex of reachableNodes(tree)) {
                    const node = nodes[nodeIndex];
                    if (node.leftChild !== BVTreeNode.invalid) {
                        continue;
                    }
                    const bv = node.boundingVolume;
                    const hit = queryType === BVTree.LINE_QUERY
                        ? AlignedBoxBV.intersectLine(P, Q, bv)
                        : (queryType === BVTree.RAY_QUERY
                            ? AlignedBoxBV.intersectRay(P, Q, bv)
                            : AlignedBoxBV.intersectSegment(P, Q, bv));
                    if (hit) {
                        expect(reported.includes(nodeIndex)).toBe(true);
                    }
                }
            }, 100);
    });

    it('copies the vertices and segments handed to createFromSegments', () => {
        check(curveArb, curve => {
            const tree = buildTree(curve);
            const before = tree.getVertices()[0].get(0);
            curve.vertices[0].set(0, curve.vertices[0].get(0) + 100);
            expect(tree.getVertices()[0].get(0)).toBe(before);
            curve.vertices[0].set(0, curve.vertices[0].get(0) - 100);
        }, 50);
    });
});
