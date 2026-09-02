import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace';
import { Hyperellipsoid } from '../src/Hyperellipsoid';
import { IntrHalfspace3Ellipsoid3TI } from '../src/IntrHalfspace3Ellipsoid3';
import { Vector, add, dot, mul, normalize } from '../src/Vector';
import { computeOrthogonalComplement3 } from '../src/Vector3';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function halfspace(normal: number[], constant: number): Halfspace {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Halfspace.fromNormalConstant(n, constant);
}

function ellipsoid(center: Vector, axis: Vector[], extent: Vector):
    Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

function axesFromDirection(w: Vector): Vector[] {
    const v = [w.clone(), Vector.zero(3), Vector.zero(3)];
    computeOrthogonalComplement3(1, v, false);
    return [v[0], v[1], v[2]];
}

// The maximum of Dot(N,X) - c over points X on the ellipsoid surface,
// computed by brute-force sampling. The query reports an intersection when
// this maximum is nonnegative.
function sampledMaxSignedDistance(h: Halfspace, e: Hyperellipsoid): number {
    let maxValue = -Number.MAX_VALUE;
    const numTheta = 120, numPhi = 120;
    for (let i = 0; i <= numTheta; ++i) {
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
            const value = dot(h.normal, X) - h.constant;
            if (value > maxValue) {
                maxValue = value;
            }
        }
    }
    return maxValue;
}

describe('IntrHalfspace3Ellipsoid3TI', () => {
    const query = new IntrHalfspace3Ellipsoid3TI();
    const unitAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];

    it('handles a unit sphere against a moving plane', () => {
        const sphere = ellipsoid(vec(0, 0, 0), unitAxes, vec(1, 1, 1));

        // The projection interval maximum is 1 - c, so the halfspace
        // x >= c intersects the sphere exactly when c <= 1.
        expect(query.test(halfspace([1, 0, 0], 0.5), sphere).intersect)
            .toBe(true);
        expect(query.test(halfspace([1, 0, 0], 1), sphere).intersect)
            .toBe(true);
        expect(query.test(halfspace([1, 0, 0], 1 + 1e-9), sphere).intersect)
            .toBe(false);
        expect(query.test(halfspace([1, 0, 0], -5), sphere).intersect)
            .toBe(true);
    });

    it('uses the ellipsoid extent along the halfspace normal', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(2, 1, 3));

        // Along x the support is 2, along y it is 1, along z it is 3.
        expect(query.test(halfspace([1, 0, 0], 1.999), e).intersect).toBe(true);
        expect(query.test(halfspace([1, 0, 0], 2.001), e).intersect).toBe(false);
        expect(query.test(halfspace([0, 1, 0], 0.999), e).intersect).toBe(true);
        expect(query.test(halfspace([0, 1, 0], 1.001), e).intersect).toBe(false);
        expect(query.test(halfspace([0, 0, 1], 2.999), e).intersect).toBe(true);
        expect(query.test(halfspace([0, 0, 1], 3.001), e).intersect).toBe(false);
    });

    it('accounts for the ellipsoid center', () => {
        const e = ellipsoid(vec(10, 0, 0), unitAxes, vec(1, 1, 1));
        expect(query.test(halfspace([1, 0, 0], 11 - 1e-6), e).intersect)
            .toBe(true);
        expect(query.test(halfspace([1, 0, 0], 11 + 1e-6), e).intersect)
            .toBe(false);
    });

    it('is invariant to the sign of the ellipsoid axes', () => {
        const e0 = ellipsoid(vec(1, 2, 3), unitAxes, vec(2, 1, 3));
        const e1 = ellipsoid(vec(1, 2, 3),
            [vec(-1, 0, 0), vec(0, -1, 0), vec(0, 0, 1)], vec(2, 1, 3));
        const h = halfspace([1, 2, -1], 4);
        expect(query.test(h, e0).intersect).toBe(query.test(h, e1).intersect);
    });

    it('rejects non-3D inputs', () => {
        const e = ellipsoid(vec(0, 0, 0), unitAxes, vec(1, 1, 1));
        const h2 = Halfspace.fromNormalConstant(
            Vector.fromArray([1, 0]), 0);
        expect(() => query.test(h2, e)).toThrow();
    });

    it('agrees with brute-force sampling on random configurations', () => {
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 40; ++trial) {
            const w = vec(rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1);
            if (Math.abs(dot(w, w)) < 1e-6) {
                continue;
            }
            normalize(w);
            const axes = axesFromDirection(w);
            const e = ellipsoid(
                vec(rand() * 4 - 2, rand() * 4 - 2, rand() * 4 - 2), axes,
                vec(0.5 + rand() * 2, 0.5 + rand() * 2, 0.5 + rand() * 2));
            const h = halfspace(
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                rand() * 6 - 3);

            const expectedMax = sampledMaxSignedDistance(h, e);
            // Skip near-tangential configurations where the sampled maximum
            // is not accurate enough to decide the sign.
            if (Math.abs(expectedMax) < 1e-2) {
                continue;
            }
            expect(query.test(h, e).intersect).toBe(expectedMax >= 0);
        }
    });
});
