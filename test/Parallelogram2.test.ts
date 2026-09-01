import { describe, it, expect } from 'vitest';
import { Parallelogram2 } from '../src/Parallelogram2';
import { Vector } from '../src/Vector';
import { dotPerp } from '../src/Vector2';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

describe('Parallelogram2 construction', () => {
    it('the default constructor is the [-1,1]^2 square', () => {
        const p = new Parallelogram2();
        expect(p.center.values).toEqual([0, 0]);
        expect(p.axis.length).toBe(2);
        expect(p.axis[0].values).toEqual([1, 0]);
        expect(p.axis[1].values).toEqual([0, 1]);
    });

    it('fromCenterAxis copies the inputs', () => {
        const center = v2(1, 2);
        const axis = [v2(2, 0), v2(1, 3)];
        const p = Parallelogram2.fromCenterAxis(center, axis);
        center.set(0, 99);
        axis[0].set(0, 99);
        expect(p.center.values).toEqual([1, 2]);
        expect(p.axis[0].values).toEqual([2, 0]);
        expect(p.axis[1].values).toEqual([1, 3]);
    });

    it('requires a right-handed basis', () => {
        // DotPerp((1,0),(0,-1)) = -1 < 0, a left-handed basis.
        expect(() => Parallelogram2.fromCenterAxis(v2(0, 0),
            [v2(1, 0), v2(0, -1)])).toThrow(
            'The axes must form a right-handed basis.');
        // Parallel axes are degenerate: DotPerp = 0.
        expect(() => Parallelogram2.fromCenterAxis(v2(0, 0),
            [v2(1, 0), v2(2, 0)])).toThrow();
        // A right-handed non-orthogonal basis is accepted.
        expect(dotPerp(v2(1, 0), v2(1, 1))).toBeGreaterThan(0);
        expect(() => Parallelogram2.fromCenterAxis(v2(0, 0),
            [v2(1, 0), v2(1, 1)])).not.toThrow();
    });

    it('rejects vectors that are not 2D', () => {
        expect(() => Parallelogram2.fromCenterAxis(Vector.fromArray([0, 0, 0]),
            [v2(1, 0), v2(0, 1)])).toThrow();
    });

    it('clone is a deep copy', () => {
        const p = new Parallelogram2();
        const copy = p.clone();
        copy.center.set(0, 5);
        copy.axis[1].set(1, 7);
        expect(p.center.values).toEqual([0, 0]);
        expect(p.axis[1].values).toEqual([0, 1]);
    });
});

describe('Parallelogram2 getVertices', () => {
    it('gives the bit-pattern-ordered corners of the default square', () => {
        const vertices = new Parallelogram2().getVertices();
        expect(vertices.length).toBe(4);
        expect(vertices[0].values).toEqual([-1, -1]);
        expect(vertices[1].values).toEqual([1, -1]);
        expect(vertices[2].values).toEqual([-1, 1]);
        expect(vertices[3].values).toEqual([1, 1]);
    });

    it('handles a sheared, translated parallelogram', () => {
        // center (1,2), axes (2,0) and (1,3).
        const p = Parallelogram2.fromCenterAxis(v2(1, 2),
            [v2(2, 0), v2(1, 3)]);
        const vertices = p.getVertices();
        expect(vertices[0].values).toEqual([-2, -1]);
        expect(vertices[1].values).toEqual([2, -1]);
        expect(vertices[2].values).toEqual([0, 5]);
        expect(vertices[3].values).toEqual([4, 5]);
    });

    it('the area of the parallelogram is 4*|DotPerp(A0,A1)|', () => {
        const p = Parallelogram2.fromCenterAxis(v2(1, 2),
            [v2(2, 0), v2(1, 3)]);
        const vertices = p.getVertices();
        // Triangles (0,1,3) and (0,3,2) tile the parallelogram.
        const area = (v: Vector[], a: number, b: number, c: number): number => {
            const e0 = Vector.fromArray([v[b].get(0) - v[a].get(0),
                v[b].get(1) - v[a].get(1)]);
            const e1 = Vector.fromArray([v[c].get(0) - v[a].get(0),
                v[c].get(1) - v[a].get(1)]);
            return 0.5 * Math.abs(dotPerp(e0, e1));
        };
        const total = area(vertices, 0, 1, 3) + area(vertices, 0, 3, 2);
        expect(total).toBeCloseTo(4 * dotPerp(p.axis[0], p.axis[1]), 12);
        expect(total).toBeCloseTo(24, 12);
    });

    it('the vertices are the center plus the four sign combinations', () => {
        const p = Parallelogram2.fromCenterAxis(v2(-3, 7),
            [v2(1, 1), v2(-1, 2)]);
        const vertices = p.getVertices();
        for (let i = 0; i < 4; ++i) {
            const s0 = 2 * (i & 1) - 1;
            const s1 = 2 * ((i >> 1) & 1) - 1;
            expect(vertices[i].values).toEqual([
                p.center.get(0) + s0 * p.axis[0].get(0)
                    + s1 * p.axis[1].get(0),
                p.center.get(1) + s0 * p.axis[0].get(1)
                    + s1 * p.axis[1].get(1)
            ]);
        }
    });
});

describe('Parallelogram2 comparisons', () => {
    const base = new Parallelogram2();

    it('equals compares the center and both axes', () => {
        expect(base.equals(new Parallelogram2())).toBe(true);
        expect(base.notEquals(new Parallelogram2())).toBe(false);

        const other = base.clone();
        other.axis[1] = v2(0, 2);
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);
    });

    it('lessThan orders by center, then the axis array lexicographically', () => {
        const smallCenter = base.clone();
        smallCenter.center = v2(-1, 0);
        expect(smallCenter.lessThan(base)).toBe(true);

        const smallAxis0 = base.clone();
        smallAxis0.axis[0] = v2(0.5, 0);
        expect(smallAxis0.lessThan(base)).toBe(true);

        // axis[0] equal, axis[1] smaller.
        const smallAxis1 = base.clone();
        smallAxis1.axis[1] = v2(0, 0.5);
        expect(smallAxis1.lessThan(base)).toBe(true);
        expect(base.lessThan(smallAxis1)).toBe(false);
    });

    it('the derived comparisons are consistent', () => {
        const bigger = base.clone();
        bigger.axis[1] = v2(0, 2);
        expect(base.lessThanOrEqual(bigger)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(bigger.greaterThan(base)).toBe(true);
        expect(bigger.greaterThanOrEqual(base)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});
