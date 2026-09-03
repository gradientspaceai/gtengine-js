import { describe, it, expect } from 'vitest';
import { ApprSphere3 } from '../src/ApprSphere3.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector } from '../src/Vector.js';

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
