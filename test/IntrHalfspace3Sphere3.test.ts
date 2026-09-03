import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import { IntrHalfspace3Sphere3TI } from '../src/IntrHalfspace3Sphere3.js';

function halfspace(nx: number, ny: number, nz: number, c: number): Halfspace {
    const n = Vector.fromArray([nx, ny, nz]);
    normalize(n);
    return Halfspace.fromNormalConstant(n, c);
}

function sphere(cx: number, cy: number, cz: number, r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray([cx, cy, cz]), r);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrHalfspace3Sphere3', () => {
    const ti = new IntrHalfspace3Sphere3TI();

    it('detects a sphere fully inside the halfspace', () => {
        // Halfspace z >= 0.
        expect(ti.test(halfspace(0, 0, 1, 0), sphere(0, 0, 5, 1)).intersect)
            .toBe(true);
    });

    it('detects a sphere fully outside the halfspace', () => {
        expect(ti.test(halfspace(0, 0, 1, 0), sphere(0, 0, -5, 1)).intersect)
            .toBe(false);
    });

    it('treats tangency as an intersection and reports the exact threshold', () => {
        // Center at z = -1 with radius 1 touches the plane z = 0.
        expect(ti.test(halfspace(0, 0, 1, 0), sphere(0, 0, -1, 1)).intersect)
            .toBe(true);
        expect(ti.test(halfspace(0, 0, 1, 0), sphere(0, 0, -1.0000001, 1)).intersect)
            .toBe(false);
    });

    it('accounts for the plane constant', () => {
        // Halfspace z >= 3.
        const h = halfspace(0, 0, 1, 3);
        expect(ti.test(h, sphere(0, 0, 2, 1)).intersect).toBe(true);
        expect(ti.test(h, sphere(0, 0, 1.5, 1)).intersect).toBe(false);
    });

    it('handles a zero-radius sphere as a point-in-halfspace test', () => {
        const h = halfspace(1, 1, 1, 0);
        expect(ti.test(h, sphere(1, 1, 1, 0)).intersect).toBe(true);
        expect(ti.test(h, sphere(0, 0, 0, 0)).intersect).toBe(true);
        expect(ti.test(h, sphere(-1, -1, -1, 0)).intersect).toBe(false);
    });

    it('agrees with the signed-distance oracle on random configurations', () => {
        const rand = makeRandom(50505);
        let numIn = 0, numOut = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const h = halfspace(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1,
                4 * rand() - 2);
            const s = sphere(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2,
                2 * rand());
            const signedDistance = dot(h.normal, s.center) - h.constant;
            const oracle = signedDistance >= -s.radius;
            const intersect = ti.test(h, s).intersect;
            expect(intersect).toBe(oracle);
            if (intersect) {
                ++numIn;
            } else {
                ++numOut;
            }
        }
        expect(numIn).toBeGreaterThan(50);
        expect(numOut).toBeGreaterThan(50);
    });
});
