import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid';
import {
    IntrEllipsoid3Ellipsoid3TI,
    IntrEllipsoid3Ellipsoid3Classification as C
} from '../src/IntrEllipsoid3Ellipsoid3';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';
import { computeOrthogonalComplement3 } from '../src/Vector3';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

function ellipsoid(center: number[], axes: Vector[], extent: number[]):
    Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(Vector.fromArray(center), axes,
        Vector.fromArray(extent));
}

function sphere(center: number[], radius: number): Hyperellipsoid {
    return ellipsoid(center, unitAxes, [radius, radius, radius]);
}

function orthonormalFrame(w: Vector): Vector[] {
    const v = [w.clone(), Vector.zero(3), Vector.zero(3)];
    normalize(v[0]);
    computeOrthogonalComplement3(1, v, false);
    return [v[0], v[1], v[2]];
}

// (X-C)^T*M*(X-C) - 1; negative strictly inside.
function quadratic(e: Hyperellipsoid, X: Vector): number {
    const diff = sub(X, e.center);
    let sum = 0;
    for (let d = 0; d < 3; ++d) {
        const t = dot(diff, e.axis[d]) / e.extent.values[d];
        sum += t * t;
    }
    return sum - 1;
}

// Sample the surface of 'e' and return the extreme values of the quadratic of
// 'other' over those samples.
function extremesOnSurface(e: Hyperellipsoid, other: Hyperellipsoid):
    { min: number, max: number } {
    let minValue = Number.MAX_VALUE, maxValue = -Number.MAX_VALUE;
    const numTheta = 90, numPhi = 90;
    for (let i = 0; i < numTheta; ++i) {
        const theta = (2 * Math.PI * i) / numTheta;
        for (let j = 0; j <= numPhi; ++j) {
            const phi = (Math.PI * j) / numPhi;
            const u = [
                Math.sin(phi) * Math.cos(theta),
                Math.sin(phi) * Math.sin(theta),
                Math.cos(phi)
            ];
            let X = e.center.clone();
            for (let d = 0; d < 3; ++d) {
                X = add(X, mul(e.extent.values[d] * u[d], e.axis[d]));
            }
            const value = quadratic(other, X);
            if (value < minValue) {
                minValue = value;
            }
            if (value > maxValue) {
                maxValue = value;
            }
        }
    }
    return { min: minValue, max: maxValue };
}

describe('IntrEllipsoid3Ellipsoid3TI', () => {
    const query = new IntrEllipsoid3Ellipsoid3TI();

    it('classifies concentric spheres', () => {
        const big = sphere([0, 0, 0], 2);
        const small = sphere([0, 0, 0], 1);
        expect(query.test(big, small)).toEqual({
            intersect: true,
            classification: C.ELLIPSOID0_CONTAINS_ELLIPSOID1
        });
        expect(query.test(small, big)).toEqual({
            intersect: true,
            classification: C.ELLIPSOID1_CONTAINS_ELLIPSOID0
        });
        expect(query.test(big, big)).toEqual({
            intersect: true,
            classification: C.ELLIPSOIDS_INTERSECTING
        });
    });

    it('classifies concentric ellipsoids that cross', () => {
        // A long thin ellipsoid poking out of a sphere.
        const s = sphere([0, 0, 0], 1);
        const needle = ellipsoid([0, 0, 0], unitAxes, [3, 0.2, 0.2]);
        const result = query.test(s, needle);
        expect(result.intersect).toBe(true);
        expect(result.classification).toBe(C.ELLIPSOIDS_INTERSECTING);
    });

    it('classifies separated spheres', () => {
        const a = sphere([0, 0, 0], 1);
        const b = sphere([5, 0, 0], 1);
        const result = query.test(a, b);
        expect(result.intersect).toBe(false);
        expect(result.classification).toBe(C.ELLIPSOIDS_SEPARATED);
    });

    it('classifies overlapping spheres', () => {
        const a = sphere([0, 0, 0], 1);
        const b = sphere([1.5, 0, 0], 1);
        const result = query.test(a, b);
        expect(result.intersect).toBe(true);
        expect(result.classification).toBe(C.ELLIPSOIDS_INTERSECTING);
    });

    it('classifies off-center containment', () => {
        const outer = sphere([0, 0, 0], 5);
        const inner = sphere([1, 0, 0], 1);
        expect(query.test(outer, inner).classification)
            .toBe(C.ELLIPSOID0_CONTAINS_ELLIPSOID1);
        expect(query.test(inner, outer).classification)
            .toBe(C.ELLIPSOID1_CONTAINS_ELLIPSOID0);
    });

    it('classifies containment of a rotated ellipsoid', () => {
        const axes = orthonormalFrame(vec(1, 1, 1));
        const outer = ellipsoid([0, 0, 0], unitAxes, [6, 6, 6]);
        const inner = ellipsoid([0.5, -0.5, 1], axes, [2, 1, 0.5]);
        expect(query.test(outer, inner).classification)
            .toBe(C.ELLIPSOID0_CONTAINS_ELLIPSOID1);
        expect(query.test(inner, outer).classification)
            .toBe(C.ELLIPSOID1_CONTAINS_ELLIPSOID0);
    });

    it('exercises the one-, two- and three-term root solvers', () => {
        // ellipsoid0 is the unit sphere, so the transformed problem has
        // D = diag(1/9,1/4,1) and K equal to the center of ellipsoid1. The
        // number of nonzero K components selects the root solver.
        const e0 = sphere([0, 0, 0], 1);
        const extents: [number, number, number] = [3, 2, 1];

        // One nonzero component: the single-term solver.
        const oneTerm = ellipsoid([4, 0, 0], unitAxes, extents);
        expect(query.test(e0, oneTerm).classification)
            .toBe(C.ELLIPSOIDS_INTERSECTING);
        const oneTermFar = ellipsoid([6, 0, 0], unitAxes, extents);
        expect(query.test(e0, oneTermFar).classification)
            .toBe(C.ELLIPSOIDS_SEPARATED);

        // Two nonzero components: the two-term solver.
        const twoTerm = ellipsoid([2, 1.5, 0], unitAxes, extents);
        const twoResult = query.test(e0, twoTerm);
        const twoExpected = extremesOnSurface(twoTerm, e0);
        expect(twoResult.intersect).toBe(twoExpected.min < 0);

        // Three nonzero components: the three-term solver.
        const threeTerm = ellipsoid([2, 1.5, 0.5], unitAxes, extents);
        const threeResult = query.test(e0, threeTerm);
        const threeExpected = extremesOnSurface(threeTerm, e0);
        expect(threeResult.intersect).toBe(threeExpected.min < 0);
    });

    it('rejects non-3D ellipsoids', () => {
        const e2 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(() => query.test(e2, e2)).toThrow();
    });

    it('agrees with the analytic sphere-sphere classification', () => {
        let seed = 1357911;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 150; ++trial) {
            const r0 = 0.3 + rand() * 2;
            const r1 = 0.3 + rand() * 2;
            const c1 = vec(rand() * 6 - 3, rand() * 6 - 3, rand() * 6 - 3);
            const distance = Math.sqrt(dot(c1, c1));
            // Skip configurations too close to a classification boundary.
            if (Math.abs(distance - (r0 + r1)) < 1e-2 ||
                Math.abs(distance + r1 - r0) < 1e-2 ||
                Math.abs(distance + r0 - r1) < 1e-2 ||
                distance < 1e-6) {
                continue;
            }

            const e0 = sphere([0, 0, 0], r0);
            const e1 = sphere([c1.values[0], c1.values[1], c1.values[2]], r1);
            const result = query.test(e0, e1);

            let expected: C;
            if (distance + r1 <= r0) {
                expected = C.ELLIPSOID0_CONTAINS_ELLIPSOID1;
            } else if (distance + r0 <= r1) {
                expected = C.ELLIPSOID1_CONTAINS_ELLIPSOID0;
            } else if (distance <= r0 + r1) {
                expected = C.ELLIPSOIDS_INTERSECTING;
            } else {
                expected = C.ELLIPSOIDS_SEPARATED;
            }
            expect(result.classification).toBe(expected);
            expect(result.intersect)
                .toBe(expected !== C.ELLIPSOIDS_SEPARATED);
        }
    });

    it('agrees with brute-force sampling for general ellipsoids', () => {
        let seed = 24681012;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        let numTested = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const axes0 = orthonormalFrame(
                vec(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1));
            const axes1 = orthonormalFrame(
                vec(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1));
            const e0 = ellipsoid([0, 0, 0], axes0,
                [0.5 + rand() * 2, 0.5 + rand() * 2, 0.5 + rand() * 2]);
            const e1 = ellipsoid(
                [rand() * 4 - 2, rand() * 4 - 2, rand() * 4 - 2], axes1,
                [0.4 + rand() * 1.5, 0.4 + rand() * 1.5, 0.4 + rand() * 1.5]);

            const on1 = extremesOnSurface(e1, e0);
            const on0 = extremesOnSurface(e0, e1);

            // Skip configurations that are numerically ambiguous.
            const tol = 5e-2;
            if (Math.abs(on1.max) < tol || Math.abs(on1.min) < tol ||
                Math.abs(on0.max) < tol) {
                continue;
            }

            let expected: C;
            if (on1.max < 0) {
                expected = C.ELLIPSOID0_CONTAINS_ELLIPSOID1;
            } else if (on0.max < 0) {
                expected = C.ELLIPSOID1_CONTAINS_ELLIPSOID0;
            } else if (on1.min < 0) {
                expected = C.ELLIPSOIDS_INTERSECTING;
            } else {
                expected = C.ELLIPSOIDS_SEPARATED;
            }

            const result = query.test(e0, e1);
            expect(result.classification).toBe(expected);
            expect(result.intersect)
                .toBe(expected !== C.ELLIPSOIDS_SEPARATED);
            ++numTested;
        }
        expect(numTested).toBeGreaterThan(20);
    });
});
