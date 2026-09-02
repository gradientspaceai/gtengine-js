import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid';
import { Hyperplane } from '../src/Hyperplane';
import {
    IntrPlane3Ellipsoid3TI,
    defaultIntrPlane3Ellipsoid3TIResult
} from '../src/IntrPlane3Ellipsoid3';
import { Vector, dot, normalize } from '../src/Vector';

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
