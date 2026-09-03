import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import {
    IntrLine3Ellipsoid3TI,
    IntrLine3Ellipsoid3FI,
    defaultIntrLine3Ellipsoid3FIResult,
    intrLine3Ellipsoid3FIDoQuery
} from '../src/IntrLine3Ellipsoid3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function line(origin: number[], direction: number[]): Line {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Line.fromOriginDirection(Vector.fromArray(origin), d);
}

function ellipsoid(center: Vector, axis: Vector[], extent: Vector):
    Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

// (X-C)^T*M*(X-C) - 1; zero on the surface, negative inside.
function quadratic(e: Hyperellipsoid, X: Vector): number {
    const diff = sub(X, e.center);
    let sum = 0;
    for (let d = 0; d < 3; ++d) {
        const t = dot(diff, e.axis[d]) / e.extent.values[d];
        sum += t * t;
    }
    return sum - 1;
}

const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

describe('IntrLine3Ellipsoid3', () => {
    const ti = new IntrLine3Ellipsoid3TI();
    const fi = new IntrLine3Ellipsoid3FI();

    it('finds the two crossings of the unit sphere', () => {
        const sphere = ellipsoid(vec(0, 0, 0), unitAxes, vec(1, 1, 1));
        const result = fi.find(line([0, 0, 0], [1, 0, 0]), sphere);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('finds the two crossings of a stretched ellipsoid', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(3, 2, 1));
        const alongX = fi.find(line([0, 0, 0], [1, 0, 0]), e);
        expect(alongX.parameter[0]).toBeCloseTo(-3, 12);
        expect(alongX.parameter[1]).toBeCloseTo(3, 12);

        const alongZ = fi.find(line([0, 0, 0], [0, 0, 1]), e);
        expect(alongZ.parameter[0]).toBeCloseTo(-1, 12);
        expect(alongZ.parameter[1]).toBeCloseTo(1, 12);
    });

    it('reports a tangent line as a single intersection', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 1, 1));
        const result = fi.find(line([0, 1, 0], [1, 0, 0]), e);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBe(result.parameter[0]);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
        expect(result.point[1].values[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection for a line that misses the ellipsoid', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 1, 1));
        const l = line([0, 2, 0], [1, 0, 0]);
        const result = fi.find(l, e);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(l, e).intersect).toBe(false);
    });

    it('agrees between the test and find queries at tangency', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 1, 1));
        const l = line([0, 1, 0], [1, 0, 0]);
        // The test query uses discr >= 0, so tangency counts.
        expect(ti.test(l, e).intersect).toBe(true);
        expect(fi.find(l, e).intersect).toBe(true);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(1, 1, 1));
        const result = defaultIntrLine3Ellipsoid3FIResult();
        const d = vec(1, 0, 0);
        intrLine3Ellipsoid3FIDoQuery(vec(0, 0, 0), d, e, result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        // The points are left at their default values.
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('rejects non-3D ellipsoids', () => {
        const e2 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(() => ti.test(line([0, 0, 0], [1, 0, 0]), e2)).toThrow();
    });

    it('puts the reported points on the ellipsoid for random lines', () => {
        let seed = 555777;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const w = vec(0.3, -0.7, 0.5);
        normalize(w);
        const basis = [w.clone(), Vector.zero(3), Vector.zero(3)];
        computeOrthogonalComplement3(1, basis, false);
        const e = ellipsoid(vec(0.5, -1, 2), basis, vec(3, 1.5, 0.75));

        let numHits = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const l = line([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(l, e);
            expect(ti.test(l, e).intersect).toBe(result.intersect);
            if (!result.intersect) {
                continue;
            }
            ++numHits;
            for (let i = 0; i < 2; ++i) {
                expect(Math.abs(quadratic(e, result.point[i])))
                    .toBeLessThan(1e-8);
                const onLine = add(l.origin,
                    mul(result.parameter[i], l.direction));
                expect(Math.sqrt(dot(sub(onLine, result.point[i]),
                    sub(onLine, result.point[i])))).toBeLessThan(1e-9);
            }
            if (result.numIntersections === 2) {
                // The midpoint of the chord is strictly inside.
                const mid = add(l.origin, mul(
                    0.5 * (result.parameter[0] + result.parameter[1]),
                    l.direction));
                expect(quadratic(e, mid)).toBeLessThan(0);
                expect(result.parameter[0]).toBeLessThan(result.parameter[1]);
            }
        }
        expect(numHits).toBeGreaterThan(10);
    });
});
