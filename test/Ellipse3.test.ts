import { describe, it, expect } from 'vitest';
import { Ellipse3 } from '../src/Ellipse3';
import { Vector, dot, sub, add, mul } from '../src/Vector';
import { cross } from '../src/Vector3';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

describe('Ellipse3 construction', () => {
    it('the default constructor is the unit circle in the x-y plane', () => {
        const ellipse = new Ellipse3();
        expect(ellipse.center.values).toEqual([0, 0, 0]);
        expect(ellipse.normal.values).toEqual([0, 0, 1]);
        expect(ellipse.axis.length).toBe(2);
        expect(ellipse.axis[0].values).toEqual([1, 0, 0]);
        expect(ellipse.axis[1].values).toEqual([0, 1, 0]);
        expect(ellipse.extent.values).toEqual([1, 1]);
    });

    it('the default frame is orthonormal and right-handed', () => {
        const e = new Ellipse3();
        expect(dot(e.axis[0], e.axis[1])).toBe(0);
        expect(cross(e.axis[0], e.axis[1]).values).toEqual(e.normal.values);
    });

    it('fromCenterNormalAxisExtent copies the inputs', () => {
        const center = v3(1, 2, 3);
        const normal = v3(0, 1, 0);
        const axis = [v3(0, 0, 1), v3(1, 0, 0)];
        const extent = Vector.fromArray([4, 2]);
        const ellipse = Ellipse3.fromCenterNormalAxisExtent(center, normal,
            axis, extent);
        center.set(0, 99);
        normal.set(1, 99);
        axis[0].set(2, 99);
        extent.set(0, 99);
        expect(ellipse.center.values).toEqual([1, 2, 3]);
        expect(ellipse.normal.values).toEqual([0, 1, 0]);
        expect(ellipse.axis[0].values).toEqual([0, 0, 1]);
        expect(ellipse.axis[1].values).toEqual([1, 0, 0]);
        expect(ellipse.extent.values).toEqual([4, 2]);
    });

    it('rejects mismatched sizes', () => {
        expect(() => Ellipse3.fromCenterNormalAxisExtent(v3(0, 0, 0),
            v3(0, 0, 1), [v3(1, 0, 0)], Vector.fromArray([1, 1]))).toThrow();
        expect(() => Ellipse3.fromCenterNormalAxisExtent(v3(0, 0, 0),
            v3(0, 0, 1), [v3(1, 0, 0), v3(0, 1, 0)],
            Vector.fromArray([1, 1, 1]))).toThrow();
    });

    it('clone is a deep copy', () => {
        const ellipse = new Ellipse3();
        const copy = ellipse.clone();
        copy.axis[0].set(0, 5);
        copy.extent.set(1, 7);
        expect(ellipse.axis[0].values).toEqual([1, 0, 0]);
        expect(ellipse.extent.values).toEqual([1, 1]);
    });
});

describe('Ellipse3 parameterization', () => {
    // X(t) = C + e0*cos(t)*A0 + e1*sin(t)*A1 lies in the plane of the ellipse
    // and satisfies (Dot(A0,X-C)/e0)^2 + (Dot(A1,X-C)/e1)^2 = 1.
    const ellipse = Ellipse3.fromCenterNormalAxisExtent(v3(1, -2, 3),
        v3(0, 1, 0), [v3(0, 0, 1), v3(1, 0, 0)], Vector.fromArray([4, 2]));

    it('parameterized points lie in the plane and on the ellipse', () => {
        const e0 = ellipse.extent.get(0);
        const e1 = ellipse.extent.get(1);
        for (let k = 0; k < 16; ++k) {
            const t = -Math.PI + (k / 16) * 2 * Math.PI;
            const x = add(add(ellipse.center,
                mul(e0 * Math.cos(t), ellipse.axis[0])),
                mul(e1 * Math.sin(t), ellipse.axis[1]));
            const diff = sub(x, ellipse.center);
            expect(dot(ellipse.normal, diff)).toBeCloseTo(0, 12);
            const s0 = dot(ellipse.axis[0], diff) / e0;
            const s1 = dot(ellipse.axis[1], diff) / e1;
            expect(s0 * s0 + s1 * s1).toBeCloseTo(1, 12);
        }
    });

    it('t = 0 gives the major-axis endpoint', () => {
        const x = add(ellipse.center, mul(4, ellipse.axis[0]));
        expect(x.values).toEqual([1, -2, 7]);
    });
});

describe('Ellipse3 comparisons', () => {
    const base = new Ellipse3();

    it('equals compares center, normal, both axes and the extent', () => {
        expect(base.equals(new Ellipse3())).toBe(true);
        expect(base.notEquals(new Ellipse3())).toBe(false);

        const other = base.clone();
        other.axis[1] = v3(0, 2, 0);
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);

        const other2 = base.clone();
        other2.extent = Vector.fromArray([1, 2]);
        expect(base.equals(other2)).toBe(false);
    });

    it('lessThan follows the upstream member order', () => {
        const smallCenter = base.clone();
        smallCenter.center = v3(-1, 0, 0);
        expect(smallCenter.lessThan(base)).toBe(true);

        const smallNormal = base.clone();
        smallNormal.normal = v3(0, 0, -1);
        expect(smallNormal.lessThan(base)).toBe(true);

        const smallAxis0 = base.clone();
        smallAxis0.axis[0] = v3(0.5, 0, 0);
        expect(smallAxis0.lessThan(base)).toBe(true);

        const smallAxis1 = base.clone();
        smallAxis1.axis[1] = v3(0, 0.5, 0);
        expect(smallAxis1.lessThan(base)).toBe(true);

        const smallExtent = base.clone();
        smallExtent.extent = Vector.fromArray([0.5, 1]);
        expect(smallExtent.lessThan(base)).toBe(true);
        expect(base.lessThan(smallExtent)).toBe(false);
    });

    it('the derived comparisons are consistent', () => {
        const bigger = base.clone();
        bigger.extent = Vector.fromArray([2, 1]);
        expect(base.lessThanOrEqual(bigger)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(bigger.greaterThan(base)).toBe(true);
        expect(bigger.greaterThanOrEqual(base)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});
