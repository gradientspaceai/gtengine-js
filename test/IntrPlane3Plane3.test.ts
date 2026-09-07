import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3Plane3TI,
    IntrPlane3Plane3FI,
    defaultIntrPlane3Plane3FIResult
} from '../src/IntrPlane3Plane3.js';
import { Vector, add, dot, mul, normalize } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function plane(normal: number[], constant: number): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalConstant(n, constant);
}

describe('IntrPlane3Plane3', () => {
    const ti = new IntrPlane3Plane3TI();
    const fi = new IntrPlane3Plane3FI();

    it('has the documented default result', () => {
        const result = defaultIntrPlane3Plane3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.isLine).toBe(false);
        expect(result.line.origin.values).toEqual([0, 0, 0]);
        expect(result.plane.normal.values).toEqual([0, 0, 0]);
    });

    it('intersects two orthogonal planes in the expected line', () => {
        const p0 = plane([1, 0, 0], 0);
        const p1 = plane([0, 1, 0], 0);
        expect(ti.test(p0, p1).intersect).toBe(true);

        const result = fi.find(p0, p1);
        expect(result.intersect).toBe(true);
        expect(result.isLine).toBe(true);
        expect(result.line.origin.values[0]).toBeCloseTo(0, 12);
        expect(result.line.origin.values[1]).toBeCloseTo(0, 12);
        expect(result.line.origin.values[2]).toBeCloseTo(0, 12);
        // Cross((1,0,0),(0,1,0)) = (0,0,1).
        expect(result.line.direction.values[0]).toBeCloseTo(0, 12);
        expect(result.line.direction.values[1]).toBeCloseTo(0, 12);
        expect(result.line.direction.values[2]).toBeCloseTo(1, 12);
    });

    it('offsets the line origin by the plane constants', () => {
        const p0 = plane([1, 0, 0], 2);
        const p1 = plane([0, 1, 0], 3);
        const result = fi.find(p0, p1);
        expect(result.isLine).toBe(true);
        expect(result.line.origin.values[0]).toBeCloseTo(2, 12);
        expect(result.line.origin.values[1]).toBeCloseTo(3, 12);
        expect(result.line.origin.values[2]).toBeCloseTo(0, 12);
    });

    it('reports coplanar planes with the same normal direction', () => {
        const p0 = plane([0, 0, 1], 4);
        const p1 = plane([0, 0, 1], 4);
        expect(ti.test(p0, p1).intersect).toBe(true);
        const result = fi.find(p0, p1);
        expect(result.intersect).toBe(true);
        expect(result.isLine).toBe(false);
        expect(result.plane.constant).toBe(4);
        expect(result.plane.normal.values).toEqual([0, 0, 1]);
    });

    it('reports coplanar planes with opposite normal directions', () => {
        const p0 = plane([0, 0, 1], 4);
        const p1 = plane([0, 0, -1], -4);
        expect(ti.test(p0, p1).intersect).toBe(true);
        const result = fi.find(p0, p1);
        expect(result.intersect).toBe(true);
        expect(result.isLine).toBe(false);
    });

    it('reports parallel but distinct planes as disjoint', () => {
        const p0 = plane([0, 0, 1], 4);
        const p1 = plane([0, 0, 1], 5);
        expect(ti.test(p0, p1).intersect).toBe(false);
        const result = fi.find(p0, p1);
        expect(result.intersect).toBe(false);
        expect(result.isLine).toBe(false);

        const p2 = plane([0, 0, -1], 5);
        expect(ti.test(p0, p2).intersect).toBe(false);
        expect(fi.find(p0, p2).intersect).toBe(false);
    });

    it('rejects non-3D planes', () => {
        const p2 = Hyperplane.fromNormalConstant(
            Vector.fromArray([1, 0]), 0);
        expect(() => ti.test(p2, p2)).toThrow();
        expect(() => fi.find(p2, p2)).toThrow();
    });

    it('produces a line lying in both planes for random pairs', () => {
        let seed = 777333;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 300; ++trial) {
            const p0 = plane([rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                rand() * 6 - 3);
            const p1 = plane([rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                rand() * 6 - 3);
            const result = fi.find(p0, p1);
            expect(ti.test(p0, p1).intersect).toBe(result.intersect);
            if (!result.isLine) {
                continue;
            }

            // The direction is unit length and orthogonal to both normals.
            expect(dot(result.line.direction, result.line.direction))
                .toBeCloseTo(1, 9);
            expect(dot(result.line.direction, p0.normal)).toBeCloseTo(0, 9);
            expect(dot(result.line.direction, p1.normal)).toBeCloseTo(0, 9);

            // Several points of the line satisfy both plane equations.
            for (const t of [-3, 0, 2.5]) {
                const X = add(result.line.origin, mul(t, result.line.direction));
                expect(dot(p0.normal, X) - p0.constant).toBeCloseTo(0, 8);
                expect(dot(p1.normal, X) - p1.constant).toBeCloseTo(0, 8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrPlane3Plane3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, wellScaled
} from './helpers/arbitraries.js';
import { length, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no normal component is a
// subnormal that would underflow when squared.
const angleP = () => wellScaled(-Math.PI, Math.PI);

function unitDirP(th: number, ph: number): Vector {
    const d = vec(Math.cos(th) * Math.cos(ph), Math.sin(th) * Math.cos(ph),
        Math.sin(ph));
    normalize(d);
    return d;
}

const planePair = fc.tuple(
    angleP(), angleP(), wellScaled(-5, 5),
    angleP(), angleP(), wellScaled(-5, 5)
).map(([t0, p0, c0, t1, p1, c1]) => ({
    plane0: Hyperplane.fromNormalConstant(unitDirP(t0, p0), c0),
    plane1: Hyperplane.fromNormalConstant(unitDirP(t1, p1), c1)
}));

describe('IntrPlane3Plane3 verification', () => {
    const tiq = new IntrPlane3Plane3TI();
    const fiq = new IntrPlane3Plane3FI();

    it('TI and FI always agree on intersect', () => {
        check(planePair, ({ plane0, plane1 }) => {
            expect(tiq.test(plane0, plane1).intersect)
                .toBe(fiq.find(plane0, plane1).intersect);
        });
    });

    it('the reported line lies in both planes', () => {
        check(planePair, ({ plane0, plane1 }) => {
            const f = fiq.find(plane0, plane1);
            if (!f.intersect || !f.isLine) { return; }

            // Two nearly parallel planes make the intersection line
            // ill-conditioned in two independent ways, so the tolerances are
            // scaled by the two condition numbers rather than fixed:
            //  - the direction is Cross(N0,N1) normalized, and |Cross| is the
            //    sine of the angle between the planes, so normalizing
            //    amplifies round-off by 1/sin;
            //  - the origin is (c0 - d*c1)*N0/(1-d^2) + ..., whose magnitude
            //    grows like 1/(1-d^2); the plane-equation residual is then
            //    eps*|origin|.
            const sinTheta = Math.max(
                length(cross(plane0.normal, plane1.normal)), 1e-12);
            const dirTol = 1e-13 / sinTheta;
            // The origin is built as c0'*N0 + c1'*N1 with
            // c0' = (c0 - d*c1)/(1 - d^2); those terms have magnitude
            // ~(|c0|+|c1|)/sin^2 and cancel down to |origin|, so the absolute
            // error of the origin is about eps times that magnitude.
            const cancel = 1 + (Math.abs(plane0.constant)
                + Math.abs(plane1.constant)) / (sinTheta * sinTheta);
            const originTol = 1e-14 * cancel;

            // The direction is unit and orthogonal to both normals.
            expectClose(length(f.line.direction), 1, 1e-12, 1e-12);
            expectClose(dot(f.line.direction, plane0.normal), 0, dirTol, 0);
            expectClose(dot(f.line.direction, plane1.normal), 0, dirTol, 0);

            // Every point of the line satisfies both plane equations.
            for (const t of [-3, -1, 0, 1, 3]) {
                const p = add(f.line.origin, mul(t, f.line.direction));
                expectClose(dot(plane0.normal, p), plane0.constant,
                    originTol, 1e-12);
                expectClose(dot(plane1.normal, p), plane1.constant,
                    originTol, 1e-12);
            }

            // The origin is the point of the line closest to the world
            // origin, because it is the combination c0*N0 + c1*N1 that lies
            // in the plane spanned by the normals.
            expectClose(dot(f.line.origin, f.line.direction), 0,
                originTol, 1e-12);
        });
    });

    it('classifies parallel and coincident planes as upstream documents',
        () => {
            check(fc.tuple(angleP(), angleP(), wellScaled(-5, 5),
                wellScaled(-5, 5), fc.boolean()),
            ([th, ph, c0, c1, flip]) => {
                const n = unitDirP(th, ph);
                const p0 = Hyperplane.fromNormalConstant(n, c0);
                // Same plane, possibly with the normal reversed.
                const same = flip
                    ? Hyperplane.fromNormalConstant(mul(-1, n), -c0)
                    : Hyperplane.fromNormalConstant(n.clone(), c0);
                const f = fiq.find(p0, same);
                expect(f.intersect).toBe(true);
                expect(f.isLine).toBe(false);
                expect(tiq.test(p0, same).intersect).toBe(true);
                // The reported plane is a copy of plane0.
                expectVectorClose(f.plane.normal, p0.normal, 0, 0);
                expect(f.plane.constant).toBe(p0.constant);
                f.plane.normal.set(0, 1234);
                expect(p0.normal.get(0)).not.toBe(1234);

                // A parallel but distinct plane never intersects.
                if (c1 !== c0) {
                    const distinct = flip
                        ? Hyperplane.fromNormalConstant(mul(-1, n), -c1)
                        : Hyperplane.fromNormalConstant(n.clone(), c1);
                    const g = fiq.find(p0, distinct);
                    expect(g.intersect).toBe(false);
                    expect(g.isLine).toBe(false);
                    expect(tiq.test(p0, distinct).intersect).toBe(false);
                }
            });
        });

    it('is symmetric under argument swap', () => {
        check(planePair, ({ plane0, plane1 }) => {
            const f = fiq.find(plane0, plane1);
            const g = fiq.find(plane1, plane0);
            expect(g.intersect).toBe(f.intersect);
            expect(g.isLine).toBe(f.isLine);
            if (!f.intersect || !f.isLine) { return; }
            // Cross(N1,N0) = -Cross(N0,N1), so the directions are opposite.
            expectVectorClose(g.line.direction, mul(-1, f.line.direction),
                1e-12, 1e-12);
            // The origins are the same point of the same line.
            expectVectorClose(g.line.origin, f.line.origin, 1e-9, 1e-10);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(planePair, angleP(), angleP(), angleP(),
            wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3)),
        ([{ plane0, plane1 }, a1, a2, a3, tx, ty, tz]) => {
            const f = fiq.find(plane0, plane1);
            if (!f.intersect || !f.isLine) { return; }
            const ca = Math.cos(a1), sa = Math.sin(a1);
            const cb = Math.cos(a2), sb = Math.sin(a2);
            const cc = Math.cos(a3), sc = Math.sin(a3);
            const b0 = vec(ca * cb, sa * cb, -sb);
            const b1 = vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc,
                cb * sc);
            const b2 = vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc,
                cb * cc);
            const rot = (v: Vector): Vector => add(mul(v.get(0), b0),
                add(mul(v.get(1), b1), mul(v.get(2), b2)));
            const t = vec(tx, ty, tz);
            const xf = (v: Vector): Vector => add(rot(v), t);
            const xfPlane = (p: Hyperplane): Hyperplane => {
                const n = rot(p.normal);
                return Hyperplane.fromNormalConstant(n,
                    p.constant + dot(n, t));
            };
            const g = fiq.find(xfPlane(plane0), xfPlane(plane1));
            expect(g.intersect).toBe(true);
            expect(g.isLine).toBe(true);
            expectVectorClose(g.line.direction, rot(f.line.direction),
                1e-8, 1e-9);
            // The transformed line is the same set of points, but its origin
            // is recomputed as the point closest to the world origin, so only
            // membership can be compared.
            for (const s of [-2, 0, 2]) {
                const p = add(g.line.origin, mul(s, g.line.direction));
                const q = sub(p, xf(f.line.origin));
                expectClose(length(sub(q,
                    mul(dot(q, g.line.direction), g.line.direction))), 0,
                1e-7, 1e-8);
            }
        });
    });

    it('reports no NaN in any field when it intersects', () => {
        check(planePair, ({ plane0, plane1 }) => {
            const f = fiq.find(plane0, plane1);
            if (!f.intersect) { return; }
            const values = f.isLine
                ? f.line.origin.values.concat(f.line.direction.values)
                : f.plane.normal.values.concat([f.plane.constant]);
            for (const v of values) {
                expect(Number.isFinite(v)).toBe(true);
            }
        });
    });

    it('rejects non-3D planes', () => {
        const p2 = Hyperplane.fromNormalConstant(Vector.fromArray([1, 0]), 0);
        const p3 = plane([0, 0, 1], 0);
        expect(() => tiq.test(p2, p3)).toThrow();
        expect(() => fiq.find(p2, p3)).toThrow();
    });
});
