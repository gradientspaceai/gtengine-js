import { describe, it, expect } from 'vitest';
import { ApprSphere3 } from '../src/ApprSphere3.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, positive, rotationFrame, vector } from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Points exactly on the sphere with the given center and radius, spread by a
// deterministic spiral over the sphere.
function spherePoints(cx: number, cy: number, cz: number, r: number,
    count: number): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        const cosPhi = -1 + 2 * (i + 0.5) / count;
        const sinPhi = Math.sqrt(1 - cosPhi * cosPhi);
        const theta = 2.399963229728653 * i;  // the golden angle
        points.push(v3(
            cx + r * sinPhi * Math.cos(theta),
            cy + r * sinPhi * Math.sin(theta),
            cz + r * cosPhi));
    }
    return points;
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprSphere3.fitUsingSquaredLengths', () => {
    it('recovers a sphere from points that lie exactly on it', () => {
        const fitter = new ApprSphere3();
        const sphere = new Hypersphere(3);
        const points = spherePoints(1, -3, 2, 4, 200);
        expect(fitter.fitUsingSquaredLengths(points, sphere)).toBe(true);
        expect(sphere.center.values[0]).toBeCloseTo(1, 8);
        expect(sphere.center.values[1]).toBeCloseTo(-3, 8);
        expect(sphere.center.values[2]).toBeCloseTo(2, 8);
        expect(sphere.radius).toBeCloseTo(4, 8);
    });

    it('fits four noncoplanar points to their circumscribed sphere', () => {
        const fitter = new ApprSphere3();
        const sphere = new Hypersphere(3);
        const points = [v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        expect(fitter.fitUsingSquaredLengths(points, sphere)).toBe(true);
        expect(sphere.center.values[0]).toBeCloseTo(0, 12);
        expect(sphere.center.values[1]).toBeCloseTo(0, 12);
        expect(sphere.center.values[2]).toBeCloseTo(0, 12);
        expect(sphere.radius).toBeCloseTo(1, 12);
    });

    it('fails for coplanar points and zeroes the sphere', () => {
        const fitter = new ApprSphere3();
        const sphere = Hypersphere.fromCenterRadius(v3(9, 9, 9), 9);
        const points = [v3(1, 0, 5), v3(-1, 0, 5), v3(0, 1, 5), v3(0, -1, 5)];
        expect(fitter.fitUsingSquaredLengths(points, sphere)).toBe(false);
        expect(sphere.center.values).toEqual([0, 0, 0]);
        expect(sphere.radius).toBe(0);
    });

    it('fails for collinear points', () => {
        const fitter = new ApprSphere3();
        const sphere = new Hypersphere(3);
        const points = [v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2), v3(-3, -3, -3)];
        expect(fitter.fitUsingSquaredLengths(points, sphere)).toBe(false);
        expect(sphere.center.values).toEqual([0, 0, 0]);
        expect(sphere.radius).toBe(0);
    });

    it('fails for coincident points', () => {
        const fitter = new ApprSphere3();
        const sphere = new Hypersphere(3);
        const points = [v3(2, 3, 4), v3(2, 3, 4), v3(2, 3, 4), v3(2, 3, 4)];
        expect(fitter.fitUsingSquaredLengths(points, sphere)).toBe(false);
        expect(sphere.radius).toBe(0);
    });

    it('is close to the true sphere for noisy samples', () => {
        const random = makeRandom(987);
        const exact = spherePoints(0, 0, 0, 3, 400);
        const points = exact.map(p => {
            const scale = 1 + 0.002 * (2 * random() - 1);
            return v3(5 + scale * p.values[0], -2 + scale * p.values[1],
                7 + scale * p.values[2]);
        });
        const fitter = new ApprSphere3();
        const sphere = new Hypersphere(3);
        expect(fitter.fitUsingSquaredLengths(points, sphere)).toBe(true);
        expect(sphere.center.values[0]).toBeCloseTo(5, 2);
        expect(sphere.center.values[1]).toBeCloseTo(-2, 2);
        expect(sphere.center.values[2]).toBeCloseTo(7, 2);
        expect(sphere.radius).toBeCloseTo(3, 2);
    });
});

describe('ApprSphere3.fitUsingLengths', () => {
    it('converges to the exact sphere when starting from the average', () => {
        const fitter = new ApprSphere3();
        const sphere = new Hypersphere(3);
        const points = spherePoints(1, -3, 2, 4, 200);
        const iterations = fitter.fitUsingLengths(points, 1024, true, sphere, 1e-14);
        expect(iterations).toBeLessThan(1024);
        expect(sphere.center.values[0]).toBeCloseTo(1, 6);
        expect(sphere.center.values[1]).toBeCloseTo(-3, 6);
        expect(sphere.center.values[2]).toBeCloseTo(2, 6);
        expect(sphere.radius).toBeCloseTo(4, 6);
    });

    it('agrees with the squared-length fit for exact data', () => {
        const points = spherePoints(-1, 0.5, 3, 2, 300);
        const fitter = new ApprSphere3();
        const sphereA = new Hypersphere(3);
        fitter.fitUsingSquaredLengths(points, sphereA);
        const sphereB = new Hypersphere(3);
        fitter.fitUsingLengths(points, 4096, true, sphereB, 0);
        for (let d = 0; d < 3; ++d) {
            expect(sphereB.center.values[d]).toBeCloseTo(sphereA.center.values[d], 5);
        }
        expect(sphereB.radius).toBeCloseTo(sphereA.radius, 5);
    });

    it('leaves the incoming sphere untouched for zero iterations', () => {
        const fitter = new ApprSphere3();
        const sphere = Hypersphere.fromCenterRadius(v3(1, 1, 1), 3);
        const points = spherePoints(0, 0, 0, 1, 20);
        expect(fitter.fitUsingLengths(points, 0, false, sphere)).toBe(0);
        expect(sphere.center.values).toEqual([1, 1, 1]);
        expect(sphere.radius).toBe(3);
    });

    it('skips samples that coincide with the current center', () => {
        const fitter = new ApprSphere3();
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 1);
        const points = [
            v3(0, 0, 0),
            v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
            v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)
        ];
        fitter.fitUsingLengths(points, 1, false, sphere);
        // Six unit-length samples and one zero-length sample.
        expect(sphere.radius).toBeCloseTo(6 / 7, 12);
        for (let d = 0; d < 3; ++d) {
            expect(sphere.center.values[d]).toBeCloseTo(0, 12);
        }
    });
});

describe('ApprSphere3 verification', () => {
    // k well-spread samples on the sphere (a Fibonacci spiral). The samples
    // span three dimensions, so the 3x3 system of the algebraic fit is well
    // conditioned.
    const spherePoints = (c: Vector, r: number, k: number): Vector[] => {
        const golden = Math.PI * (3 - Math.sqrt(5));
        const points: Vector[] = [];
        for (let i = 0; i < k; ++i) {
            const z = 1 - (2 * (i + 0.5)) / k;
            const rho = Math.sqrt(Math.max(1 - z * z, 0));
            const phi = golden * i;
            points.push(Vector.fromArray([
                c.get(0) + r * rho * Math.cos(phi),
                c.get(1) + r * rho * Math.sin(phi),
                c.get(2) + r * z]));
        }
        return points;
    };

    // The 6 samples c +/- r*axis[i] of an orthonormal frame. Their unit
    // directions from c sum to zero, which makes c a fixed point of the
    // fitUsingLengths iteration.
    const symmetricPoints = (c: Vector, r: number, axis: Vector[]): Vector[] => {
        const points: Vector[] = [];
        for (const a of axis) {
            for (const s of [r, -r]) {
                points.push(Vector.fromArray([
                    c.get(0) + s * a.get(0),
                    c.get(1) + s * a.get(1),
                    c.get(2) + s * a.get(2)]));
            }
        }
        return points;
    };

    const sphereArb = fc.tuple(vector(3, -5, 5), positive(5, 0.5),
        fc.integer({ min: 8, max: 24 }));

    it('fitUsingSquaredLengths recovers a sphere its samples lie on', () => {
        // sum_i (|X_i-C|^2 - r^2)^2 is zero at the true sphere, the global
        // minimizer that the normal equations solve for.
        check(sphereArb, ([c, r, k]) => {
            const points = spherePoints(c, r, k);
            const sphere = new Hypersphere(3);
            expect(new ApprSphere3().fitUsingSquaredLengths(points, sphere))
                .toBe(true);
            expectVectorClose(sphere.center, c, 1e-8, 1e-8);
            expectClose(sphere.radius, r, 1e-8, 1e-8);
        });
    });

    it('fitUsingSquaredLengths is equivariant under translation', () => {
        check(fc.tuple(sphereArb, vector(3, -20, 20)), ([[c, r, k], t]) => {
            const points = spherePoints(c, r, k);
            const shifted = points.map(p => Vector.fromArray([
                p.get(0) + t.get(0), p.get(1) + t.get(1), p.get(2) + t.get(2)]));

            const a = new Hypersphere(3);
            const b = new Hypersphere(3);
            const fitter = new ApprSphere3();
            expect(fitter.fitUsingSquaredLengths(points, a)).toBe(true);
            expect(fitter.fitUsingSquaredLengths(shifted, b)).toBe(true);
            for (let i = 0; i < 3; ++i) {
                expectClose(b.center.get(i), a.center.get(i) + t.get(i),
                    1e-7, 1e-7);
            }
            expectClose(b.radius, a.radius, 1e-7, 1e-7);
        });
    });

    it('fitUsingSquaredLengths fails and zeroes the sphere for coplanar and '
        + 'coincident samples', () => {
            // Samples in the z = 0 plane give M02 = M12 = M22 = 0 exactly, so
            // the cofactor expansion of the determinant is exactly zero.
            check(fc.tuple(fc.array(vector(2, -10, 10),
                { minLength: 1, maxLength: 8 }), vector(3, -10, 10)),
                ([flat, p]) => {
                    const inPlane = flat.map(
                        q => Vector.fromArray([q.get(0), q.get(1), 0]));
                    const coincident = flat.map(() => p.clone());
                    for (const points of [inPlane, coincident]) {
                        const sphere = Hypersphere.fromCenterRadius(
                            Vector.fromArray([7, 8, 9]), 10);
                        expect(new ApprSphere3()
                            .fitUsingSquaredLengths(points, sphere))
                            .toBe(false);
                        expect(sphere.center.values).toEqual([0, 0, 0]);
                        expect(sphere.radius).toBe(0);
                    }
                });
        });

    it('fitUsingLengths keeps a symmetric exact fit at its fixed point',
        () => {
            check(fc.tuple(vector(3, -5, 5), positive(5, 0.5),
                rotationFrame(3)), ([c, r, axis]) => {
                    const points = symmetricPoints(c, r, axis);
                    const sphere = new Hypersphere(3);
                    const iterations = new ApprSphere3().fitUsingLengths(
                        points, 8, true, sphere);
                    expect(iterations).toBeLessThanOrEqual(8);
                    expectVectorClose(sphere.center, c, 1e-8, 1e-8);
                    expectClose(sphere.radius, r, 1e-8, 1e-8);
                });
        });

    it('fitUsingLengths honors maxIterations and epsilon', () => {
        check(sphereArb, ([c, r, k]) => {
            const points = spherePoints(c, r, k);

            // Zero iterations leaves the incoming sphere untouched.
            const untouched = Hypersphere.fromCenterRadius(
                Vector.fromArray([1, 2, 3]), 4);
            expect(new ApprSphere3().fitUsingLengths(points, 0, false,
                untouched)).toBe(0);
            expect(untouched.center.values).toEqual([1, 2, 3]);
            expect(untouched.radius).toBe(4);

            // A huge epsilon stops after the first update, iteration 0.
            const early = new Hypersphere(3);
            expect(new ApprSphere3().fitUsingLengths(points, 25, true, early,
                1e6)).toBe(0);
        });
    });

    it('neither fit mutates its input samples', () => {
        check(sphereArb, ([c, r, k]) => {
            const points = spherePoints(c, r, k);
            const before = points.map(p => [...p.values]);
            const sphere = new Hypersphere(3);
            const fitter = new ApprSphere3();
            fitter.fitUsingSquaredLengths(points, sphere);
            fitter.fitUsingLengths(points, 4, true, sphere);
            expect(points.map(p => [...p.values])).toEqual(before);
        });
    });
});
