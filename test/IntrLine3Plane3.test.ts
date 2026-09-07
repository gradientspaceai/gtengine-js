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
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, latticeVector,
    line as arbLine, plane as arbPlane, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrLine3Plane3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrLine3Plane3 verification', () => {
    const tiQ = new IntrLine3Plane3TI();
    const fiQ = new IntrLine3Plane3FI();

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(arbLine(3), arbPlane(3)), ([L, P]) => {
            expect(fiQ.find(L, P).intersect).toBe(tiQ.test(L, P).intersect);
        });
    });

    it('the FI point is on the line and on the plane', () => {
        // The intersection parameter is -signedDistance / Dot(D,N); when the
        // line is nearly parallel to the plane the parameter (and hence the
        // point) is huge and the plane residual grows with it, so the
        // tolerance is relative to |t|.
        check(fc.tuple(arbLine(3), arbPlane(3)), ([L, P]) => {
            const r = fiQ.find(L, P);
            if (r.numIntersections !== 1) { return; }
            const onLine = add(L.origin, mul(r.parameter, L.direction));
            expectVectorClose(r.point, onLine, 1e-12, 1e-12);
            const residual = dot(P.normal, r.point) - P.constant;
            const scale = 1 + Math.abs(r.parameter)
                + Math.abs(dot(P.normal, L.origin));
            expect(Math.abs(residual)).toBeLessThanOrEqual(1e-13 * scale);
        });
    });

    it('reports a line lying in the plane as infinitely many points', () => {
        // Build a direction orthogonal to the plane normal and an origin on
        // the plane, so that Dot(D,N) and the point-plane distance are both
        // exactly zero for lattice inputs.
        check(fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -6, 6),
            latticeVector(3, -6, 6)),
            ([n, u, w]) => {
                const N = cross(u, w);
                if (dot(N, N) === 0) { return; }
                const D = cross(N, u);
                if (dot(D, D) === 0) { return; }
                const P = Hyperplane.fromNormalOrigin(N, n);
                // Move the origin n onto the plane along a plane direction.
                const O = add(n, mul(3, D));
                const L = Line.fromOriginDirection(O, D);
                const r = fiQ.find(L, P);
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(2147483647);
                expect(r.parameter).toBe(0);
                expect(r.point.values).toEqual(O.values);
                expect(tiQ.test(L, P).intersect).toBe(true);
            });
    });

    it('reports a line parallel to and off the plane as empty', () => {
        check(fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -6, 6),
            latticeVector(3, -6, 6)),
            ([n, u, w]) => {
                const N = cross(u, w);
                if (dot(N, N) === 0) { return; }
                const D = cross(N, u);
                if (dot(D, D) === 0) { return; }
                const P = Hyperplane.fromNormalOrigin(N, n);
                // Offset the line origin off the plane along the normal.
                const O = add(n, N);
                const L = Line.fromOriginDirection(O, D);
                const r = fiQ.find(L, P);
                expect(r.intersect).toBe(false);
                expect(r.numIntersections).toBe(0);
                expect(tiQ.test(L, P).intersect).toBe(false);
            });
    });

    it('is equivariant under a rigid motion of line and plane', () => {
        check(fc.tuple(arbLine(3), arbPlane(3), rotationFrame(3),
            wellScaledVector(3)),
            ([L, P, R, t]) => {
                // A line exactly parallel to the plane is a discontinuity of
                // the query (the parallel branch tests an exact zero), and a
                // rotation perturbs Dot(D,N) off zero, so only transverse
                // configurations are compared.
                if (Math.abs(dot(L.direction, P.normal)) < 1e-6) { return; }
                const map = (v: Vector) => add(Vector.fromArray([
                    dot(R[0], v), dot(R[1], v), dot(R[2], v)]), t);
                const rot = (v: Vector) => Vector.fromArray([
                    dot(R[0], v), dot(R[1], v), dot(R[2], v)]);
                const L2 = Line.fromOriginDirection(map(L.origin),
                    rot(L.direction));
                const P2 = Hyperplane.fromNormalOrigin(rot(P.normal),
                    map(P.origin));
                const r0 = fiQ.find(L, P);
                const r1 = fiQ.find(L2, P2);
                expect(r1.intersect).toBe(r0.intersect);
                expect(r1.numIntersections).toBe(r0.numIntersections);
                if (r0.numIntersections === 1) {
                    // The parameter is invariant; the rotated point drifts by
                    // the conditioning of the near-parallel configurations.
                    const scale = 1 + Math.abs(r0.parameter);
                    expectClose(r1.parameter, r0.parameter, 1e-9 * scale,
                        1e-9);
                }
            });
    });

    it('handles a zero-length line direction through the parallel branch', () => {
        const P = Hyperplane.fromNormalOrigin(Vector.fromArray([0, 0, 1]),
            Vector.fromArray([0, 0, 0]));
        const onPlane = Line.fromOriginDirection(Vector.fromArray([1, 2, 0]),
            Vector.zero(3));
        const r0 = fiQ.find(onPlane, P);
        expect(r0.intersect).toBe(true);
        expect(r0.numIntersections).toBe(2147483647);
        expect(r0.point.values).toEqual([1, 2, 0]);

        const offPlane = Line.fromOriginDirection(Vector.fromArray([1, 2, 3]),
            Vector.zero(3));
        const r1 = fiQ.find(offPlane, P);
        expect(r1.intersect).toBe(false);
        expect(r1.numIntersections).toBe(0);
        expect(tiQ.test(offPlane, P).intersect).toBe(false);
    });

    it('the exported DoQuery agrees with find', () => {
        check(fc.tuple(arbLine(3), arbPlane(3)), ([L, P]) => {
            const r = defaultIntrLine3Plane3FIResult();
            intrLine3Plane3FIDoQuery(L.origin, L.direction, P, r);
            const f = fiQ.find(L, P);
            expect(r.intersect).toBe(f.intersect);
            expect(r.numIntersections).toBe(f.numIntersections);
            expect(r.parameter).toBe(f.parameter);
            // DoQuery leaves 'point' at the default value.
            expect(r.point.values).toEqual([0, 0, 0]);
        });
    });
});
