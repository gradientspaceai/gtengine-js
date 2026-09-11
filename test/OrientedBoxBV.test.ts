import { describe, it, expect } from 'vitest';
import { OrientedBoxBV, orientedBoxBVOps } from '../src/OrientedBoxBV.js';
import { BVTree, BVTreeNode } from '../src/BVTree.js';
import { BVTreeOfPoints } from '../src/BVTreeOfPoints.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import {
    check, expectVectorClose, fc, finite, rotationFrame, unitVector,
    wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// V43 verification: the splitting axis and the three linear-component
// predicates against an independent slab clipper in box coordinates.
// ---------------------------------------------------------------------------

// Clip the parameter interval [t0, t1] of P + t*D against the slabs
// |x[k]| <= extent[k] of the box. Independent of the ported Intr* code.
function obbSlabHit(box: OrientedBox, P: Vector, D: Vector, t0: number,
    t1: number): boolean {
    const delta = sub(P, box.center);
    let tmin = t0;
    let tmax = t1;
    for (let k = 0; k < 3; ++k) {
        const p = dot(delta, box.axis[k]);
        const d = dot(D, box.axis[k]);
        const e = box.extent.get(k);
        if (d !== 0) {
            const a = (-e - p) / d;
            const b = (e - p) / d;
            tmin = Math.max(tmin, Math.min(a, b));
            tmax = Math.min(tmax, Math.max(a, b));
            if (tmin > tmax) {
                return false;
            }
        } else if (p < -e || p > e) {
            return false;
        }
    }
    return true;
}

// The slab answer is trustworthy only away from tangency, so ask it on a
// slightly shrunk and a slightly grown box; when the two agree the
// configuration is robust and the query must return that answer.
function robustObbSlab(box: OrientedBox, P: Vector, D: Vector, t0: number,
    t1: number): boolean | undefined {
    const eps = 1e-9;
    const shrunk = box.clone();
    const grown = box.clone();
    for (let k = 0; k < 3; ++k) {
        shrunk.extent.set(k, box.extent.get(k) - eps);
        grown.extent.set(k, box.extent.get(k) + eps);
    }
    const a = obbSlabHit(shrunk, P, D, t0, t1);
    const b = obbSlabHit(grown, P, D, t0, t1);
    return a === b ? a : undefined;
}

describe('OrientedBoxBV verification', () => {
    const bvArb = fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
        fc.array(finite(0.25, 4), { minLength: 3, maxLength: 3 }))
        .map(([center, axis, e]) => OrientedBoxBV.fromBox(
            OrientedBox.fromCenterAxisExtent(center, axis,
                Vector.fromArray(e))));

    it('splits through the center along the axis of largest extent', () => {
        check(bvArb, bv => {
            const { origin, direction } = bv.getSplittingAxis();
            expectVectorClose(origin, bv.box.center, 0, 0);
            let maxIndex = 0;
            for (let k = 1; k < 3; ++k) {
                if (bv.box.extent.get(k) > bv.box.extent.get(maxIndex)) {
                    maxIndex = k;
                }
            }
            expectVectorClose(direction, bv.box.axis[maxIndex], 0, 0);
            // The axis is returned as a copy of the box axis, not an alias.
            direction.set(0, direction.get(0) + 1);
            expect(bv.box.axis[maxIndex].get(0))
                .not.toBe(direction.get(0));
        });
    });

    it('agrees with the slab clipper for lines, rays and segments', () => {
        check(fc.tuple(bvArb, wellScaledVector(3, -8, 8), unitVector(3),
            finite(0.5, 12)), ([bv, P, D, len]) => {
            const Q = add(P, mul(D, len));
            const line = robustObbSlab(bv.box, P, D, -Infinity, Infinity);
            if (line !== undefined) {
                expect(OrientedBoxBV.intersectLine(P, D, bv)).toBe(line);
            }
            const ray = robustObbSlab(bv.box, P, D, 0, Infinity);
            if (ray !== undefined) {
                expect(OrientedBoxBV.intersectRay(P, D, bv)).toBe(ray);
            }
            const seg = robustObbSlab(bv.box, P, sub(Q, P), 0, 1);
            if (seg !== undefined) {
                expect(OrientedBoxBV.intersectSegment(P, Q, bv)).toBe(seg);
            }
        });
    });

    it('is monotone in segment, ray and line', () => {
        check(fc.tuple(bvArb, wellScaledVector(3, -8, 8), unitVector(3),
            finite(0.5, 12)), ([bv, P, D, len]) => {
            const Q = add(P, mul(D, len));
            if (OrientedBoxBV.intersectSegment(P, Q, bv)) {
                expect(OrientedBoxBV.intersectRay(P, D, bv)).toBe(true);
            }
            if (OrientedBoxBV.intersectRay(P, D, bv)) {
                expect(OrientedBoxBV.intersectLine(P, D, bv)).toBe(true);
            }
        });
    });

    it('hits when an endpoint is inside the box', () => {
        check(fc.tuple(bvArb, wellScaledVector(3, -0.9, 0.9),
            wellScaledVector(3, -8, 8)), ([bv, t, Q]) => {
            let P = bv.box.center.clone();
            for (let k = 0; k < 3; ++k) {
                P = add(P, mul(bv.box.axis[k],
                    t.get(k) * bv.box.extent.get(k)));
            }
            expect(OrientedBoxBV.intersectSegment(P, Q, bv)).toBe(true);
        });
    });

    it('copies the box in fromBox', () => {
        check(bvArb, bv => {
            const copy = OrientedBoxBV.fromBox(bv.box);
            copy.box.center.set(0, copy.box.center.get(0) + 100);
            copy.box.axis[0].set(1, 42);
            expect(bv.box.center.get(0)).not.toBe(copy.box.center.get(0));
            expect(bv.box.axis[0].get(1)).not.toBe(42);
        });
    });
});
