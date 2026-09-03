import { describe, it, expect } from 'vitest';
import { OrientedBoxBV, orientedBoxBVOps } from '../src/OrientedBoxBV.js';
import { BVTree, BVTreeNode } from '../src/BVTree.js';
import { BVTreeOfPoints } from '../src/BVTreeOfPoints.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, sub } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

const r = Math.SQRT1_2;

// A box rotated by 45 degrees about the z-axis, centered at the origin, with
// extents (2, 1, 1).
function rotatedBox(): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(v3(0, 0, 0),
        [v3(r, r, 0), v3(-r, r, 0), v3(0, 0, 1)], v3(2, 1, 1));
}

function bvOf(box: OrientedBox): OrientedBoxBV {
    const bv = new OrientedBoxBV();
    bv.box = box;
    return bv;
}

describe('OrientedBoxBV', () => {
    it('default-constructs the OrientedBox3 default', () => {
        const bv = new OrientedBoxBV();
        expect(bv.box.dimension).toBe(3);
        expect(bv.box.center.values).toEqual([0, 0, 0]);
        expect(bv.box.extent.values).toEqual([1, 1, 1]);
        expect(bv.box.axis[0].values).toEqual([1, 0, 0]);
        expect(bv.box.axis[1].values).toEqual([0, 1, 0]);
        expect(bv.box.axis[2].values).toEqual([0, 0, 1]);
    });

    it('copies the box in fromBox', () => {
        const box = rotatedBox();
        const bv = OrientedBoxBV.fromBox(box);
        expect(bv.box.center.values).toEqual([0, 0, 0]);
        box.center.set(0, 100);
        expect(bv.box.center.get(0)).toBe(0);
    });

    it('splits along the box axis of largest extent, through the center', () => {
        const box = rotatedBox();
        box.center = v3(3, -4, 5);
        let axis = bvOf(box).getSplittingAxis();
        expect(axis.origin.values).toEqual([3, -4, 5]);
        expect(axis.direction.values).toEqual([r, r, 0]);

        box.extent = v3(1, 7, 2);
        axis = bvOf(box).getSplittingAxis();
        expect(axis.direction.values).toEqual([-r, r, 0]);

        box.extent = v3(1, 2, 9);
        axis = bvOf(box).getSplittingAxis();
        expect(axis.direction.values).toEqual([0, 0, 1]);
    });

    it('breaks extent ties in favor of the smaller axis index', () => {
        // The upstream comparisons are strict, so equal extents keep the
        // earlier axis.
        const box = rotatedBox();
        box.extent = v3(2, 2, 2);
        expect(bvOf(box).getSplittingAxis().direction.values).toEqual([r, r, 0]);
        box.extent = v3(1, 3, 3);
        expect(bvOf(box).getSplittingAxis().direction.values).toEqual([-r, r, 0]);
    });

    it('returns copies from getSplittingAxis, not aliases of the box', () => {
        const box = rotatedBox();
        const bv = bvOf(box);
        const axis = bv.getSplittingAxis();
        axis.origin.set(0, 12);
        axis.direction.set(0, 12);
        expect(bv.box.center.get(0)).toBe(0);
        expect(bv.box.axis[0].get(0)).toBe(r);
    });

    it('tests the linear components against the rotated box', () => {
        const bv = bvOf(rotatedBox());

        // Along the long axis, through the center.
        const P = v3(-10 * r, -10 * r, 0);
        const Q = v3(r, r, 0);
        expect(OrientedBoxBV.intersectLine(P, Q, bv)).toBe(true);
        expect(OrientedBoxBV.intersectRay(P, Q, bv)).toBe(true);
        // A ray pointing away from the box misses; the line does not.
        expect(OrientedBoxBV.intersectRay(P, v3(-r, -r, 0), bv)).toBe(false);
        expect(OrientedBoxBV.intersectLine(P, v3(-r, -r, 0), bv)).toBe(true);

        // The segment endpoints are P and Q, not P and P+Q. This segment
        // stops short of the box.
        expect(OrientedBoxBV.intersectSegment(P, v3(-5 * r, -5 * r, 0), bv))
            .toBe(false);
        expect(OrientedBoxBV.intersectSegment(P, v3(10 * r, 10 * r, 0), bv))
            .toBe(true);

        // Vertical lines at offsets along the second box axis: |y1| <= 1 hits.
        const inside = v3(-0.5 * r, 0.5 * r, 0);
        const outside = v3(-1.5 * r, 1.5 * r, 0);
        expect(OrientedBoxBV.intersectLine(inside, v3(0, 0, 1), bv)).toBe(true);
        expect(OrientedBoxBV.intersectLine(outside, v3(0, 0, 1), bv)).toBe(false);

        // A ray from above the box, aimed down, hits; aimed up, misses.
        const above = v3(0, 0, 10);
        expect(OrientedBoxBV.intersectRay(above, v3(0, 0, -1), bv)).toBe(true);
        expect(OrientedBoxBV.intersectRay(above, v3(0, 0, 1), bv)).toBe(false);
    });

    it('agrees with the box-coordinate containment test on random queries', () => {
        // Randomized cross-check: a segment from an interior point to an
        // exterior point must intersect; a segment between two points on the
        // same side of a slab face must not.
        const box = rotatedBox();
        const bv = bvOf(box);
        let s = 12345 >>> 0;
        const rand = (): number => {
            s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
            return s / 4294967296;
        };

        for (let trial = 0; trial < 200; ++trial) {
            const y = [
                (2 * rand() - 1) * 1.9,
                (2 * rand() - 1) * 0.9,
                (2 * rand() - 1) * 0.9
            ];
            const interior = Vector.fromArray([0, 0, 0]);
            for (let k = 0; k < 3; ++k) {
                for (let j = 0; j < 3; ++j) {
                    interior.set(j, interior.get(j) + y[k] * box.axis[k].get(j));
                }
            }
            // Sanity: the constructed point really is inside.
            for (let k = 0; k < 3; ++k) {
                expect(Math.abs(dot(sub(interior, box.center), box.axis[k])))
                    .toBeLessThanOrEqual(box.extent.get(k) + 1e-12);
            }

            const far = v3(20 + 10 * rand(), 30 * rand(), 40 * rand());
            expect(OrientedBoxBV.intersectSegment(interior, far, bv)).toBe(true);
            expect(OrientedBoxBV.intersectLine(interior, v3(1, 0, 0), bv)).toBe(true);

            // Two points beyond the +z face of the box are separated from it.
            const p0 = v3(10 * rand() - 5, 10 * rand() - 5, 1.5 + rand());
            const p1 = v3(10 * rand() - 5, 10 * rand() - 5, 1.5 + rand());
            expect(OrientedBoxBV.intersectSegment(p0, p1, bv)).toBe(false);
        }
    });

    it('bundles the static operations into orientedBoxBVOps', () => {
        const created = orientedBoxBVOps.create();
        expect(created).toBeInstanceOf(OrientedBoxBV);
        expect(created.box.extent.values).toEqual([1, 1, 1]);

        const bv = bvOf(rotatedBox());
        const P = v3(-10 * r, -10 * r, 0);
        const Q = v3(r, r, 0);
        expect(orientedBoxBVOps.intersectLine(P, Q, bv)).toBe(true);
        expect(orientedBoxBVOps.intersectRay(P, Q, bv)).toBe(true);
        expect(orientedBoxBVOps.intersectRay(P, v3(-r, -r, 0), bv)).toBe(false);
        expect(orientedBoxBVOps.intersectSegment(P, v3(10 * r, 10 * r, 0), bv))
            .toBe(true);
    });
});

// A minimal concrete BVTree that uses OrientedBoxBV as its bounding volume.
// The boxes are axis-aligned oriented boxes (the box axes are the coordinate
// axes), which keeps the expected values independent of any eigensolver.
class AxisAlignedOBBTreeOfPoints extends BVTreeOfPoints<OrientedBoxBV> {
    constructor() {
        super(orientedBoxBVOps);
    }

    leafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }

    private fit(indices: readonly number[], bv: OrientedBoxBV): void {
        const min = [Infinity, Infinity, Infinity];
        const max = [-Infinity, -Infinity, -Infinity];
        for (const i of indices) {
            const p = this.mVertices[i];
            for (let k = 0; k < 3; ++k) {
                min[k] = Math.min(min[k], p.get(k));
                max[k] = Math.max(max[k], p.get(k));
            }
        }
        for (let k = 0; k < 3; ++k) {
            bv.box.center.set(k, 0.5 * (min[k] + max[k]));
            bv.box.extent.set(k, 0.5 * (max[k] - min[k]));
            bv.box.axis[k] = Vector.unit(3, k);
        }
    }

    protected override computeInteriorBoundingVolume(i0: number, i1: number,
        bv: OrientedBoxBV): void {
        const indices: number[] = [];
        for (let i = i0; i <= i1; ++i) {
            indices.push(this.mPartition[i]);
        }
        this.fit(indices, bv);
    }

    protected override computeLeafBoundingVolume(i: number,
        bv: OrientedBoxBV): void {
        this.fit([this.mPartition[i]], bv);
    }
}

describe('OrientedBoxBV as a BVTree bounding volume', () => {
    it('builds a tree whose nodes bound their points and finds hit leaves', () => {
        const points: Vector[] = [];
        let s = 777 >>> 0;
        const rand = (): number => {
            s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
            return s / 4294967296;
        };
        for (let i = 0; i < 32; ++i) {
            points.push(v3(10 * rand() - 5, 10 * rand() - 5, 10 * rand() - 5));
        }

        const tree = new AxisAlignedOBBTreeOfPoints();
        tree.create(points, BVTree.fullHeight);

        // Every node's box contains the points of its range.
        const partition = tree.getPartition();
        for (const node of tree.getNodes()) {
            if (node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            const box = node.boundingVolume.box;
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                const p = points[partition[i]];
                for (let k = 0; k < 3; ++k) {
                    expect(Math.abs(dot(sub(p, box.center), box.axis[k])))
                        .toBeLessThanOrEqual(box.extent.get(k) + 1e-12);
                }
            }
        }

        // A line through a data point reports the leaf holding that point.
        for (let t = 0; t < points.length; t += 5) {
            const nodeIndices = tree.leafIndices(BVTree.LINE_QUERY, points[t],
                v3(0, 0, 1));
            const candidates = new Set<number>();
            for (const nodeIndex of nodeIndices) {
                const node = tree.getNodes()[nodeIndex];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    candidates.add(partition[i]);
                }
            }
            expect(candidates.has(t)).toBe(true);
        }
    });
});
