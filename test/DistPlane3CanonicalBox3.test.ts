import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPlane3CanonicalBox3 } from '../src/DistPlane3CanonicalBox3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { Vector, dot, length, normalize, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(...extent: number[]): CanonicalBox {
    return CanonicalBox.fromExtent(v(...extent));
}

function plane(normal: number[], constant: number): Hyperplane {
    return Hyperplane.fromNormalConstant(v(...normal), constant);
}

function unitPlane(normal: number[], constant: number): Hyperplane {
    const n = v(...normal);
    normalize(n);
    return Hyperplane.fromNormalConstant(n, constant);
}

// The exact plane-box distance for a unit-length normal. The box is centered
// at the origin, so Dot(N,X) ranges over [-R,R] with R = sum_i e[i]*|n[i]|.
function exactDistance(p: Hyperplane, b: CanonicalBox): number {
    let radius = 0;
    for (let i = 0; i < 3; ++i) {
        radius += b.extent.values[i] * Math.abs(p.normal.values[i]);
    }
    return Math.max(0, Math.abs(p.constant) - radius);
}

function verifyClosest(p: Hyperplane, b: CanonicalBox,
    result: { distance: number, closest: [Vector, Vector] }): void {
    // closest[0] is on the plane.
    expect(dot(p.normal, result.closest[0])).toBeCloseTo(p.constant, 9);
    // closest[1] is in the box.
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(result.closest[1].values[i]))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
    }
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistPlane3CanonicalBox3', () => {
    const query = new DistPlane3CanonicalBox3();
    const unitBox = box(1, 1, 1);

    it('computes the distance for an axis-aligned plane', () => {
        const p = plane([0, 0, 1], 2);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(1, 12);
        expect(result.sqrDistance).toBeCloseTo(1, 12);
        verifyClosest(p, unitBox, result);
    });

    it('reports zero distance when the plane cuts the box', () => {
        const p = plane([0, 0, 1], 0.5);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        verifyClosest(p, unitBox, result);
    });

    it('handles a plane touching a box face', () => {
        const p = plane([1, 0, 0], 1);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('computes the distance for an oblique plane', () => {
        // The unit normal (1,1,1)/sqrt(3) has box radius sqrt(3).
        const p = unitPlane([1, 1, 1], 3);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(3 - Math.sqrt(3), 10);
        verifyClosest(p, unitBox, result);
    });

    it('handles a normal with a negative component (reflection path)', () => {
        const p = unitPlane([-1, 1, 0], 4);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 10);
        verifyClosest(p, unitBox, result);
    });

    it('handles the one-nonzero-component normal path', () => {
        const b = box(1, 2, 3);
        const p = plane([1, 0, 0], 5);
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(4, 12);
        // The upstream 1D path reports the plane point (p0,e1,e2) and the box
        // point (clamp(p0),e1,e2).
        expect(result.closest[0].values).toEqual([5, 2, 3]);
        expect(result.closest[1].values).toEqual([1, 2, 3]);
    });

    it('handles a degenerate zero normal', () => {
        // A zero normal is degenerate; upstream falls back to a point-box
        // query using the stored plane origin.
        const p = Hyperplane.fromNormalOrigin(v(0, 0, 0), v(3, 0, 0));
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[0].values).toEqual([3, 0, 0]);
        expect(result.closest[1].values).toEqual([1, 0, 0]);
    });

    it('handles a degenerate box with zero extents', () => {
        const b = box(0, 0, 0);
        const p = plane([0, 0, 1], 2);
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('agrees with the exact formula on random inputs', () => {
        let seed = 1029384756;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 400; ++trial) {
            const n = v(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (dot(n, n) < 1e-6) {
                continue;
            }
            normalize(n);
            const p = Hyperplane.fromNormalConstant(n, 6 * rand() - 3);
            const b = box(0.2 + 2 * rand(), 0.2 + 2 * rand(),
                0.2 + 2 * rand());
            const result = query.compute(p, b);
            expect(result.distance).toBeCloseTo(exactDistance(p, b), 9);
            verifyClosest(p, b, result);
        }
    });

    it('is symmetric under negation of the plane', () => {
        let seed = 5647382910 % 2147483648;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 100; ++trial) {
            const n = v(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (dot(n, n) < 1e-6) {
                continue;
            }
            normalize(n);
            const c = 6 * rand() - 3;
            const b = box(0.5 + rand(), 0.5 + rand(), 0.5 + rand());
            const r0 = query.compute(Hyperplane.fromNormalConstant(n, c), b);
            const negated = v(-n.values[0], -n.values[1], -n.values[2]);
            const r1 = query.compute(
                Hyperplane.fromNormalConstant(negated, -c), b);
            expect(r1.distance).toBeCloseTo(r0.distance, 9);
        }
    });
});
