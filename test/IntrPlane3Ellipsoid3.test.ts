import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3Ellipsoid3TI,
    defaultIntrPlane3Ellipsoid3TIResult
} from '../src/IntrPlane3Ellipsoid3.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { mulMatrix } from '../src/Matrix.js';
import {
    check, expectClose, fc, plane as arbPlane, positive, rotationFrame,
    wellScaledVector
} from './helpers/arbitraries.js';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

// An axis-aligned ellipsoid with the given center and semi-axis lengths.
function ellipsoid(center: number[], extent: number[]): Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(Vector.fromArray(center),
        [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
            Vector.fromArray([0, 0, 1])],
        Vector.fromArray(extent));
}

// An independent test: the support of the ellipsoid in the plane-normal
// direction is sqrt(sum_i (e_i * Dot(N, U_i))^2), so the plane and ellipsoid
// intersect when |signed distance of the center| <= that support.
function bruteForceIntersect(P: Hyperplane, E: Hyperellipsoid): boolean {
    let sum = 0;
    for (let i = 0; i < 3; ++i) {
        const term = E.extent.values[i] * dot(P.normal, E.axis[i]);
        sum += term * term;
    }
    const sd = dot(P.normal, E.center) - P.constant;
    return Math.abs(sd) <= Math.sqrt(sum);
}

const ti = new IntrPlane3Ellipsoid3TI();

describe('IntrPlane3Ellipsoid3', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrPlane3Ellipsoid3TIResult().intersect).toBe(false);
    });

    it('uses the semi-axis length along the plane normal', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        // Along z the support is 1.
        expect(ti.test(plane([0, 0, 1], [0, 0, 0.9]), E).intersect).toBe(true);
        expect(ti.test(plane([0, 0, 1], [0, 0, 1.1]), E).intersect).toBe(false);
        // Along x the support is 3.
        expect(ti.test(plane([1, 0, 0], [2.9, 0, 0]), E).intersect).toBe(true);
        expect(ti.test(plane([1, 0, 0], [3.1, 0, 0]), E).intersect).toBe(false);
    });

    it('reports the tangent plane as an intersection', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        expect(ti.test(plane([0, 0, 1], [0, 0, 1]), E).intersect).toBe(true);
    });

    it('reduces to the sphere case for equal extents', () => {
        const E = ellipsoid([1, 2, 3], [2, 2, 2]);
        expect(ti.test(plane([1, 1, 1], [1, 2, 3]), E).intersect).toBe(true);
        // A plane at distance 2.5 from the center misses the radius-2 sphere.
        const n = Vector.fromArray([1, 1, 1]);
        normalize(n);
        const origin = Vector.fromArray([
            1 + 2.5 * n.values[0], 2 + 2.5 * n.values[1], 3 + 2.5 * n.values[2]
        ]);
        expect(ti.test(Hyperplane.fromNormalOrigin(n, origin), E).intersect)
            .toBe(false);
    });

    it('agrees with the support-function test on random inputs', () => {
        let state = 555111;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 400; ++trial) {
            // A rotated ellipsoid: rotate the standard frame about z.
            const angle = rand() * Math.PI;
            const c = Math.cos(angle);
            const s = Math.sin(angle);
            const E = Hyperellipsoid.fromCenterAxisExtent(
                Vector.fromArray([rand() * 3, rand() * 3, rand() * 3]),
                [Vector.fromArray([c, s, 0]), Vector.fromArray([-s, c, 0]),
                    Vector.fromArray([0, 0, 1])],
                Vector.fromArray([0.5 + Math.abs(rand()) * 2,
                    0.5 + Math.abs(rand()) * 2, 0.5 + Math.abs(rand()) * 2]));
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand() * 3, rand() * 3, rand() * 3]);
            const actual = ti.test(P, E).intersect;
            expect(actual).toBe(bruteForceIntersect(P, E));
            if (actual) {
                ++numHits;
            }
        }
        expect(numHits).toBeGreaterThan(50);
        expect(numHits).toBeLessThan(350);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrPlane3Ellipsoid3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrPlane3Ellipsoid3 verification', () => {
    const q = new IntrPlane3Ellipsoid3TI();

    const arbEllipsoid = fc.tuple(wellScaledVector(3), rotationFrame(3),
        fc.tuple(positive(4), positive(4), positive(4)))
        .map(([c, axis, e]) => Hyperellipsoid.fromCenterAxisExtent(c, axis,
            Vector.fromArray([e[0], e[1], e[2]])));

    it('agrees with the exact support point of the ellipsoid', () => {
        // The ellipsoid point maximizing Dot(N, X - C) is
        //   X = C + M^{-1} N / sqrt(N^T M^{-1} N),
        // so the extreme signed distances are Dot(N,C) - c +- sqrt(...). The
        // support point is built here from the port's own getM/getMInverse
        // and checked to be on the ellipsoid, which makes this an independent
        // derivation of the query's comparison.
        check(fc.tuple(arbPlane(3), arbEllipsoid), ([P, E]) => {
            const MInv = E.getMInverse();
            const M = E.getM();
            const MiN = mulMatrix(MInv, P.normal) as Vector;
            const s = dot(P.normal, MiN);
            expect(s).toBeGreaterThan(0);
            const root = Math.sqrt(s);
            const support = add(E.center, mul(1 / root, MiN));
            // The support point is on the ellipsoid: (X-C)^T M (X-C) = 1.
            const d = sub(support, E.center);
            expectClose(dot(d, mulMatrix(M, d) as Vector), 1, 1e-7, 1e-7);
            const center = dot(P.normal, E.center) - P.constant;
            const expected = Math.abs(center) <= root;
            expect(q.test(P, E).intersect).toBe(expected);
        });
    });

    it('reports an intersection whenever sampled surface points straddle the plane', () => {
        check(fc.tuple(arbPlane(3), arbEllipsoid), ([P, E]) => {
            let lo = Number.POSITIVE_INFINITY;
            let hi = Number.NEGATIVE_INFINITY;
            for (let i = 0; i <= 24; ++i) {
                const phi = (Math.PI * i) / 24;
                for (let j = 0; j < 48; ++j) {
                    const theta = (2 * Math.PI * j) / 48;
                    const X = add(E.center, add(
                        mul(E.extent.values[0] * Math.sin(phi)
                            * Math.cos(theta), E.axis[0]),
                        add(mul(E.extent.values[1] * Math.sin(phi)
                            * Math.sin(theta), E.axis[1]),
                            mul(E.extent.values[2] * Math.cos(phi),
                                E.axis[2]))));
                    const sd = dot(P.normal, X) - P.constant;
                    lo = Math.min(lo, sd);
                    hi = Math.max(hi, sd);
                }
            }
            const got = q.test(P, E).intersect;
            // The sampled range is a subset of the true range, so a straddle
            // proves an intersection; the reverse direction is asserted only
            // with a margin that the sampling density can support.
            if (lo < -1e-9 && hi > 1e-9) { expect(got).toBe(true); }
            const maxExtent = Math.max(E.extent.values[0], E.extent.values[1],
                E.extent.values[2]);
            if (lo > 0.05 * maxExtent || hi < -0.05 * maxExtent) {
                expect(got).toBe(false);
            }
        }, 60);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbPlane(3), arbEllipsoid, rotationFrame(3),
            wellScaledVector(3)),
            ([P, E, R, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(R[0], v), dot(R[1], v), dot(R[2], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const P2 = Hyperplane.fromNormalOrigin(rot(P.normal),
                    map(P.origin));
                const E2 = Hyperellipsoid.fromCenterAxisExtent(map(E.center),
                    E.axis.map(rot), E.extent);
                const MiN = mulMatrix(E.getMInverse(), P.normal) as Vector;
                const root = Math.sqrt(dot(P.normal, MiN));
                const center = Math.abs(dot(P.normal, E.center) - P.constant);
                // Skip near-tangency, where rounding of the rotated copy can
                // cross the '<=' comparison.
                if (Math.abs(center - root) < 1e-7 * (1 + root)) { return; }
                expect(q.test(P2, E2).intersect)
                    .toBe(q.test(P, E).intersect);
            });
    });

    it('reduces to the sphere query for equal extents', () => {
        check(fc.tuple(arbPlane(3), wellScaledVector(3), rotationFrame(3),
            positive(4)),
            ([P, c, axis, r]) => {
                const E = Hyperellipsoid.fromCenterAxisExtent(c, axis,
                    Vector.fromArray([r, r, r]));
                const d = Math.abs(dot(P.normal, c) - P.constant);
                // sqrt(N^T M^{-1} N) = r for a sphere and a unit normal.
                const got = q.test(P, E).intersect;
                if (Math.abs(d - r) > 1e-9 * (1 + r)) {
                    expect(got).toBe(d <= r);
                }
            });
    });

    it('detects a plane tangent to an axis-aligned ellipsoid', () => {
        const E = Hyperellipsoid.fromCenterAxisExtent(Vector.zero(3),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([3, 2, 1]));
        const tangent = Hyperplane.fromNormalOrigin(
            Vector.fromArray([1, 0, 0]), Vector.fromArray([3, 0, 0]));
        expect(q.test(tangent, E).intersect).toBe(true);
        const outside = Hyperplane.fromNormalOrigin(
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([3 + 1e-12, 0, 0]));
        expect(q.test(outside, E).intersect).toBe(false);
    });
});
