import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrTriangle3AlignedBox3TI,
    IntrTriangle3AlignedBox3FI
} from '../src/IntrTriangle3AlignedBox3.js';
import {
    IntrTriangle3OrientedBox3TI,
    IntrTriangle3OrientedBox3FI,
    defaultIntrTriangle3OrientedBox3FIResult,
    defaultIntrTriangle3OrientedBox3TIResult
} from '../src/IntrTriangle3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function tri(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(v0), Vector.fromArray(v1),
        Vector.fromArray(v2));
}

// A right-handed orthonormal frame from a rotation about the given axis.
function rotationFrame(axis: Vector, angle: number): Vector[] {
    const u = axis.clone();
    normalize(u);
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)].map(
        (e) => add(add(mul(c, e), mul(s, cross(u, e))),
            mul((1 - c) * dot(u, e), u)));
}

// True when the point is in the oriented box (with tolerance).
function inOBB(B: OrientedBox, X: Vector, tol: number): boolean {
    const diff = sub(X, B.center);
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(dot(B.axis[i], diff)) > B.extent.values[i] + tol) {
            return false;
        }
    }
    return true;
}

function planeDistance(T: Triangle, X: Vector): number {
    const n = cross(sub(T.v[1], T.v[0]), sub(T.v[2], T.v[0]));
    return dot(n, sub(X, T.v[0])) / Math.sqrt(dot(n, n));
}

function inTriangle(T: Triangle, X: Vector, tol: number): boolean {
    const n = cross(sub(T.v[1], T.v[0]), sub(T.v[2], T.v[0]));
    const area2 = dot(n, n);
    for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
        const e = sub(T.v[i1], T.v[i0]);
        const d = sub(X, T.v[i0]);
        if (dot(cross(e, d), n) < -tol * area2) {
            return false;
        }
    }
    return true;
}

const ti = new IntrTriangle3OrientedBox3TI();
const fi = new IntrTriangle3OrientedBox3FI();
const abTI = new IntrTriangle3AlignedBox3TI();
const abFI = new IntrTriangle3AlignedBox3FI();

const identityAxes = [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)];
const unitOBB = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0), identityAxes,
    vec(1, 1, 1));

describe('IntrTriangle3OrientedBox3', () => {
    it('default-constructs the results as empty', () => {
        expect(defaultIntrTriangle3OrientedBox3TIResult())
            .toEqual({ intersect: false });
        const r = defaultIntrTriangle3OrientedBox3FIResult();
        expect(r.insidePolygon).toEqual([]);
        expect(r.outsidePolygons).toEqual([]);
    });

    it('keeps a triangle fully inside the box', () => {
        const T = tri([-0.5, -0.5, 0], [0.5, -0.5, 0], [0, 0.5, 0]);
        expect(ti.test(T, unitOBB).intersect).toBe(true);
        const result = fi.find(T, unitOBB);
        expect(result.insidePolygon.length).toBe(3);
        expect(result.outsidePolygons.length).toBe(0);
    });

    it('rejects a triangle entirely outside the box', () => {
        const T = tri([5, 5, 5], [6, 5, 5], [5, 6, 5]);
        expect(ti.test(T, unitOBB).intersect).toBe(false);
        expect(fi.find(T, unitOBB).insidePolygon.length).toBe(0);
    });

    it('clips a large triangle to a rotated box cross section', () => {
        // The box is rotated 45 degrees about z, so the slice at z = 0 is the
        // square of "radius" 1 rotated by 45 degrees, that is, the diamond
        // with vertices at distance sqrt(2) along the axes.
        const box = OrientedBox.fromCenterAxisExtent(vec(0, 0, 0),
            rotationFrame(vec(0, 0, 1), Math.PI / 4), vec(1, 1, 1));
        const T = tri([-10, -10, 0], [10, -10, 0], [0, 10, 0]);
        const result = fi.find(T, box);
        expect(result.insidePolygon.length).toBe(4);
        for (const p of result.insidePolygon) {
            expect(inOBB(box, p, 1e-10)).toBe(true);
            expect(p.values[2]).toBeCloseTo(0, 12);
            expect(Math.abs(p.values[0]) + Math.abs(p.values[1]))
                .toBeCloseTo(Math.SQRT2, 10);
        }
    });

    it('matches the aligned-box query when the axes are the standard basis', () => {
        const B = AlignedBox.fromMinMax(vec(-1, -1, -1), vec(1, 1, 1));
        let seed = 33445566;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        for (let trial = 0; trial < 600; ++trial) {
            const T = tri(
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]);
            expect(ti.test(T, unitOBB).intersect).toBe(abTI.test(T, B).intersect);
            const obbResult = fi.find(T, unitOBB);
            const abResult = abFI.find(T, B);
            expect(obbResult.insidePolygon.length)
                .toBe(abResult.insidePolygon.length);
            for (let i = 0; i < obbResult.insidePolygon.length; ++i) {
                for (let j = 0; j < 3; ++j) {
                    expect(obbResult.insidePolygon[i].values[j])
                        .toBeCloseTo(abResult.insidePolygon[i].values[j], 10);
                }
            }
        }
    });

    it('is invariant under a rigid motion of triangle and box', () => {
        const axes = rotationFrame(vec(2, -1, 0.5), 0.9);
        const translation = vec(3, -2, 4);
        const mapPoint = (p: Vector): Vector => {
            let q = translation.clone();
            for (let i = 0; i < 3; ++i) {
                q = add(q, mul(p.values[i], axes[i]));
            }
            return q;
        };
        const mapDirection = (p: Vector): Vector => {
            let q = Vector.zero(3);
            for (let i = 0; i < 3; ++i) {
                q = add(q, mul(p.values[i], axes[i]));
            }
            return q;
        };

        let seed = 90210777;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        let numInside = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const T = tri(
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]);
            const before = fi.find(T, unitOBB);
            if (before.insidePolygon.length > 0) {
                ++numInside;
            }
            const movedT = Triangle.fromVertices(mapPoint(T.v[0]),
                mapPoint(T.v[1]), mapPoint(T.v[2]));
            const movedBox = OrientedBox.fromCenterAxisExtent(
                mapPoint(unitOBB.center), unitOBB.axis.map(mapDirection),
                unitOBB.extent);
            const after = fi.find(movedT, movedBox);
            expect(after.insidePolygon.length)
                .toBe(before.insidePolygon.length);
            for (let i = 0; i < before.insidePolygon.length; ++i) {
                const expected = mapPoint(before.insidePolygon[i]);
                for (let j = 0; j < 3; ++j) {
                    expect(after.insidePolygon[i].values[j])
                        .toBeCloseTo(expected.values[j], 8);
                }
            }
            expect(ti.test(movedT, movedBox).intersect)
                .toBe(ti.test(T, unitOBB).intersect);
        }
        expect(numInside).toBeGreaterThan(50);
    });

    it('agrees with TI and produces valid geometry (randomized)', () => {
        let seed = 11223344;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        const box = OrientedBox.fromCenterAxisExtent(vec(0.25, -0.5, 0.75),
            rotationFrame(vec(1, 2, -1), 0.8), vec(1, 0.5, 1.5));

        let numInside = 0;
        let numClipped = 0;
        for (let trial = 0; trial < 1200; ++trial) {
            const T = tri(
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]);
            const tiResult = ti.test(T, box).intersect;
            const fiResult = fi.find(T, box);
            if (fiResult.insidePolygon.length > 0) {
                ++numInside;
                expect(tiResult).toBe(true);
                if (fiResult.insidePolygon.length > 3) {
                    ++numClipped;
                }
                for (const p of fiResult.insidePolygon) {
                    expect(inOBB(box, p, 1e-9)).toBe(true);
                    expect(Math.abs(planeDistance(T, p))).toBeLessThan(1e-7);
                    expect(inTriangle(T, p, 1e-7)).toBe(true);
                }
            }
            for (const poly of fiResult.outsidePolygons) {
                for (const p of poly) {
                    expect(Math.abs(planeDistance(T, p))).toBeLessThan(1e-7);
                }
            }
        }
        expect(numInside).toBeGreaterThan(80);
        expect(numClipped).toBeGreaterThan(30);
    });
});
