import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import { IntrHalfspace3OrientedBox3TI } from '../src/IntrHalfspace3OrientedBox3.js';

function halfspace(nx: number, ny: number, nz: number, c: number): Halfspace {
    const n = Vector.fromArray([nx, ny, nz]);
    normalize(n);
    return Halfspace.fromNormalConstant(n, c);
}

function rotationAxes(ax: number, ay: number, az: number): Vector[] {
    const cx = Math.cos(ax), sx = Math.sin(ax);
    const cy = Math.cos(ay), sy = Math.sin(ay);
    const cz = Math.cos(az), sz = Math.sin(az);
    const m = [
        [cz * cy, cz * sy * sx - sz * cx, cz * sy * cx + sz * sx],
        [sz * cy, sz * sy * sx + cz * cx, sz * sy * cx - cz * sx],
        [-sy, cy * sx, cy * cx]
    ];
    return [
        Vector.fromArray([m[0][0], m[1][0], m[2][0]]),
        Vector.fromArray([m[0][1], m[1][1], m[2][1]]),
        Vector.fromArray([m[0][2], m[1][2], m[2][2]])
    ];
}

function box(center: number[], axes: Vector[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center), axes,
        Vector.fromArray(extent));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrHalfspace3OrientedBox3', () => {
    const ti = new IntrHalfspace3OrientedBox3TI();
    const identity = rotationAxes(0, 0, 0);

    it('detects an axis-aligned box inside, outside and straddling', () => {
        const h = halfspace(0, 0, 1, 0);  // z >= 0
        expect(ti.test(h, box([0, 0, 5], identity, [1, 1, 1])).intersect).toBe(true);
        expect(ti.test(h, box([0, 0, -5], identity, [1, 1, 1])).intersect).toBe(false);
        expect(ti.test(h, box([0, 0, 0], identity, [1, 1, 1])).intersect).toBe(true);
    });

    it('treats face contact as an intersection and finds the exact threshold', () => {
        const h = halfspace(0, 0, 1, 0);
        expect(ti.test(h, box([0, 0, -1], identity, [1, 1, 1])).intersect).toBe(true);
        expect(ti.test(h, box([0, 0, -1.0000001], identity, [1, 1, 1])).intersect)
            .toBe(false);
    });

    it('uses the projected radius of a rotated box', () => {
        // A unit cube rotated 45 degrees about z has half-width sqrt(2) along
        // x, so it reaches the plane x >= 0 when its center is at -sqrt(2).
        const h = halfspace(1, 0, 0, 0);
        const axes = rotationAxes(0, 0, Math.PI / 4);
        expect(ti.test(h, box([-Math.SQRT2 + 1e-12, 0, 0], axes, [1, 1, 1])).intersect)
            .toBe(true);
        expect(ti.test(h, box([-Math.SQRT2 - 1e-6, 0, 0], axes, [1, 1, 1])).intersect)
            .toBe(false);
        // Without the rotation, contact happens at -1 instead.
        expect(ti.test(h, box([-1.2, 0, 0], identity, [1, 1, 1])).intersect)
            .toBe(false);
    });

    it('handles a degenerate box (a point)', () => {
        const h = halfspace(1, 1, 1, 0);
        const axes = rotationAxes(0.2, 0.3, 0.4);
        expect(ti.test(h, box([0, 0, 0], axes, [0, 0, 0])).intersect).toBe(true);
        expect(ti.test(h, box([-1, -1, -1], axes, [0, 0, 0])).intersect).toBe(false);
    });

    it('agrees with a vertex-enumeration oracle on random configurations', () => {
        const rand = makeRandom(112233);
        let numIn = 0, numOut = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const h = halfspace(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1,
                4 * rand() - 2);
            const b = box([4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2],
                rotationAxes(2 * Math.PI * rand(), 2 * Math.PI * rand(),
                    2 * Math.PI * rand()),
                [0.1 + rand(), 0.1 + rand(), 0.1 + rand()]);

            // Oracle: the box intersects the halfspace iff some vertex has a
            // nonnegative signed distance (the box is a convex hull of its
            // vertices).
            let oracle = false;
            for (const v of b.getVertices()) {
                if (dot(h.normal, v) - h.constant >= 0) {
                    oracle = true;
                }
            }

            const intersect = ti.test(h, b).intersect;
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
