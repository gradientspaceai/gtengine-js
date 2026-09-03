import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrLine3Plane3TI,
    IntrLine3Plane3FI,
    defaultIntrLine3Plane3FIResult,
    intrLine3Plane3FIDoQuery
} from '../src/IntrLine3Plane3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function line(origin: number[], direction: number[]): Line {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Line.fromOriginDirection(Vector.fromArray(origin), d);
}

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

const ti = new IntrLine3Plane3TI();
const fi = new IntrLine3Plane3FI();

describe('IntrLine3Plane3', () => {
    it('finds a transverse intersection at a known point', () => {
        const L = line([0, 0, 5], [0, 0, -1]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        const result = fi.find(L, P);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.point.values).toEqual([0, 0, 2]);
        expect(ti.test(L, P).intersect).toBe(true);
    });

    it('reports a parallel disjoint line as no intersection', () => {
        const L = line([0, 0, 5], [1, 0, 0]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        const result = fi.find(L, P);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(L, P).intersect).toBe(false);
    });

    it('reports a coincident line with the int32 max sentinel', () => {
        const L = line([1, 2, 2], [1, 1, 0]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        const result = fi.find(L, P);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2147483647);
        expect(result.parameter).toBe(0);
        // The reported point is the line origin.
        expect(result.point.values).toEqual([1, 2, 2]);
        expect(ti.test(L, P).intersect).toBe(true);
    });

    it('exposes the DoQuery helper without computing the point', () => {
        const result = defaultIntrLine3Plane3FIResult();
        intrLine3Plane3FIDoQuery(vec(0, 0, 5), vec(0, 0, -1),
            plane([0, 0, 1], [0, 0, 2]), result);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(3, 12);
        // DoQuery leaves 'point' at its default value.
        expect(result.point.values).toEqual([0, 0, 0]);
    });

    it('agrees with a direct signed-distance computation on random inputs', () => {
        let state = 20250901;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const L = line([rand() * 4, rand() * 4, rand() * 4],
                [rand(), rand(), rand() + 0.001]);
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand() * 3, rand() * 3, rand() * 3]);
            const result = fi.find(L, P);
            expect(ti.test(L, P).intersect).toBe(result.intersect);
            if (result.numIntersections === 1) {
                ++numHits;
                // The reported point lies on the plane and on the line.
                expect(dot(P.normal, result.point) - P.constant)
                    .toBeCloseTo(0, 10);
                const onLine = add(L.origin, mul(result.parameter, L.direction));
                for (let i = 0; i < 3; ++i) {
                    expect(result.point.values[i]).toBeCloseTo(onLine.values[i], 12);
                }
            }
        }
        expect(numHits).toBeGreaterThan(300);
    });
});
