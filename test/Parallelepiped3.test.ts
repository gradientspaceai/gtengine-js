import { describe, it, expect } from 'vitest';
import { Parallelepiped3 } from '../src/Parallelepiped3.js';
import { Vector, add, mul } from '../src/Vector.js';
import { Polyhedron3 } from '../src/Polyhedron3.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, finite,
    positive, rotationFrame, vector } from './helpers/arbitraries.js';
import { dotCross } from '../src/Vector3.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

describe('Parallelepiped3 construction', () => {
    it('the default constructor is the [-1,1]^3 cube', () => {
        const p = new Parallelepiped3();
        expect(p.center.values).toEqual([0, 0, 0]);
        expect(p.axis.length).toBe(3);
        expect(p.axis[0].values).toEqual([1, 0, 0]);
        expect(p.axis[1].values).toEqual([0, 1, 0]);
        expect(p.axis[2].values).toEqual([0, 0, 1]);
    });

    it('fromCenterAxis copies the inputs', () => {
        const center = v3(1, 2, 3);
        const axis = [v3(2, 0, 0), v3(0, 3, 0), v3(0, 0, 4)];
        const p = Parallelepiped3.fromCenterAxis(center, axis);
        center.set(0, 99);
        axis[2].set(2, 99);
        expect(p.center.values).toEqual([1, 2, 3]);
        expect(p.axis[2].values).toEqual([0, 0, 4]);
    });

    it('requires a right-handed basis', () => {
        // Swapping two axes makes the basis left-handed.
        expect(() => Parallelepiped3.fromCenterAxis(v3(0, 0, 0),
            [v3(0, 1, 0), v3(1, 0, 0), v3(0, 0, 1)])).toThrow(
            'The axes must form a right-handed basis.');
        // Coplanar axes are degenerate: DotCross = 0.
        expect(() => Parallelepiped3.fromCenterAxis(v3(0, 0, 0),
            [v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0)])).toThrow();
        // A right-handed non-orthogonal basis is accepted.
        const sheared = [v3(1, 0, 0), v3(1, 1, 0), v3(1, 1, 1)];
        expect(dotCross(sheared[0], sheared[1], sheared[2])).toBe(1);
        expect(() => Parallelepiped3.fromCenterAxis(v3(0, 0, 0),
            sheared)).not.toThrow();
    });

    it('rejects vectors that are not 3D', () => {
        expect(() => Parallelepiped3.fromCenterAxis(Vector.fromArray([0, 0]),
            [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)])).toThrow();
    });

    it('clone is a deep copy', () => {
        const p = new Parallelepiped3();
        const copy = p.clone();
        copy.center.set(0, 5);
        copy.axis[2].set(2, 7);
        expect(p.center.values).toEqual([0, 0, 0]);
        expect(p.axis[2].values).toEqual([0, 0, 1]);
    });
});

describe('Parallelepiped3 getVertices', () => {
    it('gives the bit-pattern-ordered corners of the default cube', () => {
        const vertices = new Parallelepiped3().getVertices();
        expect(vertices.length).toBe(8);
        expect(vertices[0].values).toEqual([-1, -1, -1]);
        expect(vertices[1].values).toEqual([1, -1, -1]);
        expect(vertices[2].values).toEqual([-1, 1, -1]);
        expect(vertices[3].values).toEqual([1, 1, -1]);
        expect(vertices[4].values).toEqual([-1, -1, 1]);
        expect(vertices[5].values).toEqual([1, -1, 1]);
        expect(vertices[6].values).toEqual([-1, 1, 1]);
        expect(vertices[7].values).toEqual([1, 1, 1]);
    });

    it('handles a scaled, translated box', () => {
        const p = Parallelepiped3.fromCenterAxis(v3(1, 2, 3),
            [v3(2, 0, 0), v3(0, 3, 0), v3(0, 0, 4)]);
        const vertices = p.getVertices();
        expect(vertices[0].values).toEqual([-1, -1, -1]);
        expect(vertices[7].values).toEqual([3, 5, 7]);
        expect(vertices[3].values).toEqual([3, 5, -1]);
        expect(vertices[4].values).toEqual([-1, -1, 7]);
    });

    it('the vertices are the center plus the eight sign combinations', () => {
        const p = Parallelepiped3.fromCenterAxis(v3(-3, 7, 2),
            [v3(1, 0, 0), v3(1, 1, 0), v3(1, 1, 1)]);
        const vertices = p.getVertices();
        for (let i = 0; i < 8; ++i) {
            const s = [2 * (i & 1) - 1, 2 * ((i >> 1) & 1) - 1,
                2 * ((i >> 2) & 1) - 1];
            for (let c = 0; c < 3; ++c) {
                const expected = p.center.get(c) + s[0] * p.axis[0].get(c)
                    + s[1] * p.axis[1].get(c) + s[2] * p.axis[2].get(c);
                expect(vertices[i].get(c)).toBe(expected);
            }
        }
    });

    it('the volume is 8*DotCross(A0,A1,A2)', () => {
        // For the axis-aligned case the box is 4-by-6-by-8 = 192, and
        // 8 * DotCross((2,0,0),(0,3,0),(0,0,4)) = 8 * 24 = 192.
        const axis = [v3(2, 0, 0), v3(0, 3, 0), v3(0, 0, 4)];
        expect(8 * dotCross(axis[0], axis[1], axis[2])).toBe(192);
        const p = Parallelepiped3.fromCenterAxis(v3(1, 2, 3), axis);
        const vertices = p.getVertices();
        // Extents along each axis direction are twice the axis lengths.
        expect(vertices[1].get(0) - vertices[0].get(0)).toBe(4);
        expect(vertices[2].get(1) - vertices[0].get(1)).toBe(6);
        expect(vertices[4].get(2) - vertices[0].get(2)).toBe(8);
    });
});

describe('Parallelepiped3 comparisons', () => {
    const base = new Parallelepiped3();

    it('equals compares the center and all three axes', () => {
        expect(base.equals(new Parallelepiped3())).toBe(true);
        expect(base.notEquals(new Parallelepiped3())).toBe(false);

        const other = base.clone();
        other.axis[2] = v3(0, 0, 2);
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);
    });

    it('lessThan orders by center, then the axis array lexicographically', () => {
        const smallCenter = base.clone();
        smallCenter.center = v3(-1, 0, 0);
        expect(smallCenter.lessThan(base)).toBe(true);

        const smallAxis0 = base.clone();
        smallAxis0.axis[0] = v3(0.5, 0, 0);
        expect(smallAxis0.lessThan(base)).toBe(true);

        const smallAxis2 = base.clone();
        smallAxis2.axis[2] = v3(0, 0, 0.5);
        expect(smallAxis2.lessThan(base)).toBe(true);
        expect(base.lessThan(smallAxis2)).toBe(false);
    });

    it('the derived comparisons are consistent', () => {
        const bigger = base.clone();
        bigger.axis[2] = v3(0, 0, 2);
        expect(base.lessThanOrEqual(bigger)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(bigger.greaterThan(base)).toBe(true);
        expect(bigger.greaterThanOrEqual(base)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});

describe('Parallelepiped3 verification', () => {
    // Right-handed axes: DotCross(A0, A1, A2) > 0 is required by the factory.
    const parallelepiped = () => fc.tuple(vector(3, -5, 5), rotationFrame(3),
        fc.array(positive(4, 0.1), { minLength: 3, maxLength: 3 }),
        fc.array(finite(-0.6, 0.6), { minLength: 3, maxLength: 3 }))
        .map(([c, frame, s, shear]) => {
            const a0 = mul(s[0], frame[0]);
            const a1 = add(mul(s[1], frame[1]), mul(shear[0] * s[1], frame[0]));
            const a2 = add(mul(s[2], frame[2]),
                add(mul(shear[1] * s[2], frame[0]),
                    mul(shear[2] * s[2], frame[1])));
            return Parallelepiped3.fromCenterAxis(c, [a0, a1, a2]);
        });
    const key = (p: Parallelepiped3) => [...p.center.values,
        ...p.axis.flatMap(a => [...a.values])];

    it('getVertices emits bit-pattern order, not the documented CCW order',
        () => {
            // Upstream #155: the comment says "counterclockwise order" but
            // the code emits vertices[i] = C + sum sign[d]*A[d] with
            // i = b[2]b[1]b[0] and sign[d] = 2*b[d] - 1.
            check(parallelepiped(), p => {
                const v = p.getVertices();
                expect(v.length).toBe(8);
                for (let i = 0; i < 8; ++i) {
                    let expected = p.center;
                    for (let d = 0; d < 3; ++d) {
                        const sign = ((i >> d) & 1) === 1 ? 1 : -1;
                        expected = add(expected, mul(sign, p.axis[d]));
                    }
                    for (let d = 0; d < 3; ++d) {
                        expectClose(v[i].get(d), expected.get(d),
                            1e-12, 1e-12);
                    }
                }
            });
        });

    it('the volume of the vertex hull is 8*|DotCross(A0,A1,A2)|', () => {
        // Cross-check against an independent computation: the closed
        // triangulation of the eight corners, measured by Polyhedron3.
        check(parallelepiped(), p => {
            const v = p.getVertices();
            // Faces of the bit-indexed box, each split into two triangles.
            // Winding is irrelevant: ComputeVolume takes the absolute value.
            const faces = [
                [0, 1, 3, 2], [4, 6, 7, 5],   // -A2, +A2
                [0, 4, 5, 1], [2, 3, 7, 6],   // -A1, +A1
                [0, 2, 6, 4], [1, 5, 7, 3]];  // -A0, +A0
            const indices: number[] = [];
            for (const [a, b, c, d] of faces) {
                indices.push(a, b, c, a, c, d);
            }
            const polyhedron = new Polyhedron3(v, indices.length, indices,
                true);
            expectClose(polyhedron.computeVolume(),
                8 * Math.abs(dotCross(p.axis[0], p.axis[1], p.axis[2])),
                1e-9, 1e-9);
        }, 50);
    });

    it('opposite corners average to the center', () => {
        check(parallelepiped(), p => {
            const v = p.getVertices();
            for (let i = 0; i < 8; ++i) {
                for (let d = 0; d < 3; ++d) {
                    expectClose(0.5 * (v[i].get(d) + v[7 - i].get(d)),
                        p.center.get(d), 1e-12, 1e-12);
                }
            }
        });
    });

    it('the factory rejects a left-handed or degenerate basis', () => {
        check(fc.tuple(vector(3, -5, 5), rotationFrame(3)), ([c, frame]) => {
            expect(() => Parallelepiped3.fromCenterAxis(c,
                [frame[1], frame[0], frame[2]]))
                .toThrow('The axes must form a right-handed basis.');
            // A zero axis makes DotCross exactly zero. (A repeated axis
            // would be zero in exact arithmetic but can round to a tiny
            // positive value, which the upstream '> 0' test accepts.)
            expect(() => Parallelepiped3.fromCenterAxis(c,
                [frame[0], frame[1], new Vector(3)]))
                .toThrow('The axes must form a right-handed basis.');
        });
    });

    it('the comparisons follow the (center, axis) member order', () => {
        check(fc.tuple(parallelepiped(), parallelepiped()), ([a, b]) => {
            const cmp = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(cmp < 0);
            expect(a.greaterThan(b)).toBe(cmp > 0);
            expect(a.lessThanOrEqual(b)).toBe(cmp <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(cmp >= 0);
            expect(a.equals(b)).toBe(cmp === 0);
        });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(parallelepiped(), { minLength: 4, maxLength: 5 }),
            ps => {
                expectStrictWeakOrder(ps, (x, y) => x.lessThan(y));
            }, 30);
    });

    it('equals is element equality: a NaN axis breaks self-equality', () => {
        // Regression: upstream operator== is 'center == other.center &&
        // axis == other.axis' with std::array's element-by-element ==.
        const p = new Parallelepiped3();
        p.axis[0] = Vector.fromArray([NaN, 0, 0]);
        expect(p.equals(p)).toBe(false);
        expect(p.notEquals(p)).toBe(true);
        expect(p.lessThan(p)).toBe(false);
        expect(p.lessThanOrEqual(p)).toBe(true);
    });

    it('the factory and clone copy every vector', () => {
        check(parallelepiped(), p => {
            const axis = p.axis.map(a => a.clone());
            const copy = Parallelepiped3.fromCenterAxis(p.center, axis);
            const cloned = copy.clone();
            axis[0].set(0, 999);
            copy.axis[1].set(1, 888);
            expect(copy.axis[0].get(0)).not.toBe(999);
            expect(cloned.axis[1].get(1)).not.toBe(888);
        });
    });
});
