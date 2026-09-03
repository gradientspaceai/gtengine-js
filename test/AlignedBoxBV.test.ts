import { describe, it, expect } from 'vitest';
import { AlignedBoxBV } from '../src/AlignedBoxBV.js';
import { AlignedBox } from '../src/AlignedBox.js';
import { Vector, add, mul, normalize, sub } from '../src/Vector.js';

function V(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function unit(x: number, y: number, z: number): Vector {
    const v = V(x, y, z);
    normalize(v);
    return v;
}

function makeBV(min: [number, number, number],
    max: [number, number, number]): AlignedBoxBV {
    return AlignedBoxBV.fromBox(AlignedBox.fromMinMax(V(...min), V(...max)));
}

describe('AlignedBoxBV', () => {
    it('default-constructs a degenerate box at the origin', () => {
        const bv = new AlignedBoxBV();
        expect(bv.box.min.values).toEqual([0, 0, 0]);
        expect(bv.box.max.values).toEqual([0, 0, 0]);
    });

    it('copies the box in fromBox', () => {
        const box = AlignedBox.fromMinMax(V(-1, -2, -3), V(4, 5, 6));
        const bv = AlignedBoxBV.fromBox(box);
        box.min.values[0] = 100;
        expect(bv.box.min.values).toEqual([-1, -2, -3]);
        expect(bv.box.max.values).toEqual([4, 5, 6]);
    });

    it('splits along the axis of largest extent, breaking ties toward x', () => {
        // Extents (3,1,2): the largest is along x.
        let s = makeBV([-3, -1, -2], [3, 1, 2]).getSplittingAxis();
        expect(s.origin.values).toEqual([0, 0, 0]);
        expect(s.direction.values).toEqual([1, 0, 0]);

        // Extents (1,3,2): the largest is along y. The center is not zero.
        s = makeBV([0, 0, 0], [2, 6, 4]).getSplittingAxis();
        expect(s.origin.values).toEqual([1, 3, 2]);
        expect(s.direction.values).toEqual([0, 1, 0]);

        // Extents (1,2,3): the largest is along z.
        s = makeBV([-1, -2, -3], [1, 2, 3]).getSplittingAxis();
        expect(s.direction.values).toEqual([0, 0, 1]);

        // A cube: the strict comparisons keep the first maximum, x.
        s = makeBV([-1, -1, -1], [1, 1, 1]).getSplittingAxis();
        expect(s.direction.values).toEqual([1, 0, 0]);

        // Ties between y and z: y is chosen because z is not strictly larger.
        s = makeBV([-1, -2, -2], [1, 2, 2]).getSplittingAxis();
        expect(s.direction.values).toEqual([0, 1, 0]);
    });

    it('tests line intersection', () => {
        const bv = makeBV([-1, -1, -1], [1, 1, 1]);

        // A line through the box.
        expect(AlignedBoxBV.intersectLine(V(0, 0, -5), unit(0, 0, 1), bv))
            .toBe(true);
        // The same line's opposite direction: a line is unbounded, so it
        // still intersects.
        expect(AlignedBoxBV.intersectLine(V(0, 0, -5), unit(0, 0, -1), bv))
            .toBe(true);
        // A line missing the box.
        expect(AlignedBoxBV.intersectLine(V(5, 5, 0), unit(0, 0, 1), bv))
            .toBe(false);
        // A line grazing the corner (1,1,1).
        expect(AlignedBoxBV.intersectLine(V(1, 1, 1), unit(1, -1, 0), bv))
            .toBe(true);
    });

    it('tests ray intersection', () => {
        const bv = makeBV([-1, -1, -1], [1, 1, 1]);

        // The ray points at the box.
        expect(AlignedBoxBV.intersectRay(V(0, 0, -5), unit(0, 0, 1), bv))
            .toBe(true);
        // The ray points away from the box.
        expect(AlignedBoxBV.intersectRay(V(0, 0, -5), unit(0, 0, -1), bv))
            .toBe(false);
        // The ray origin is inside the box.
        expect(AlignedBoxBV.intersectRay(V(0, 0, 0), unit(1, 2, 3), bv))
            .toBe(true);
        // The ray misses the box.
        expect(AlignedBoxBV.intersectRay(V(5, 5, -5), unit(0, 0, 1), bv))
            .toBe(false);
    });

    it('tests segment intersection', () => {
        const bv = makeBV([-1, -1, -1], [1, 1, 1]);

        // The segment crosses the box.
        expect(AlignedBoxBV.intersectSegment(V(0, 0, -5), V(0, 0, 5), bv))
            .toBe(true);
        // The segment stops short of the box.
        expect(AlignedBoxBV.intersectSegment(V(0, 0, -5), V(0, 0, -2), bv))
            .toBe(false);
        // The segment is entirely inside the box.
        expect(AlignedBoxBV.intersectSegment(V(-0.5, 0, 0), V(0.5, 0, 0), bv))
            .toBe(true);
        // The segment touches the face x = 1 at one endpoint.
        expect(AlignedBoxBV.intersectSegment(V(1, 0, 0), V(3, 0, 0), bv))
            .toBe(true);
        // A segment that misses.
        expect(AlignedBoxBV.intersectSegment(V(2, 2, 2), V(3, 3, 3), bv))
            .toBe(false);
    });

    it('agrees with brute-force sampling on random segments', () => {
        let seed = 20240613;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const R = (): number => 6 * rand() - 3;

        const bv = makeBV([-1, -0.5, -2], [1.5, 2, 0.5]);
        const box = bv.box;
        const inBox = (p: Vector): boolean => {
            for (let k = 0; k < 3; ++k) {
                if (p.values[k] < box.min.values[k] ||
                    p.values[k] > box.max.values[k]) {
                    return false;
                }
            }
            return true;
        };

        let numHits = 0;
        let numMisses = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const p0 = V(R(), R(), R());
            const p1 = V(R(), R(), R());
            const intersect = AlignedBoxBV.intersectSegment(p0, p1, bv);

            // Sample the segment densely. Any sampled point inside the box
            // proves an intersection; the converse is not guaranteed, so the
            // sampling can only confirm hits.
            let sampledHit = false;
            const delta = sub(p1, p0);
            for (let i = 0; i <= 400; ++i) {
                if (inBox(add(p0, mul(delta, i / 400)))) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(intersect).toBe(true);
                ++numHits;
            }
            else if (!intersect) {
                ++numMisses;
            }

            // The ray from p0 through p1 must intersect whenever the segment
            // does, and the line must intersect whenever the ray does.
            const direction = sub(p1, p0);
            normalize(direction);
            const rayHit = AlignedBoxBV.intersectRay(p0, direction, bv);
            const lineHit = AlignedBoxBV.intersectLine(p0, direction, bv);
            if (intersect) {
                expect(rayHit).toBe(true);
            }
            if (rayHit) {
                expect(lineHit).toBe(true);
            }
        }
        expect(numHits).toBeGreaterThan(5);
        expect(numMisses).toBeGreaterThan(5);
    });
});
