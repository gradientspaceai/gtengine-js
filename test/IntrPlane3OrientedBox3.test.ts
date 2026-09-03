import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3OrientedBox3TI,
    defaultIntrPlane3OrientedBox3TIResult
} from '../src/IntrPlane3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, normalize } from '../src/Vector.js';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

// An oriented box rotated by 'angle' about the z-axis.
function box(center: number[], angle: number, extent: number[]): OrientedBox {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center),
        [Vector.fromArray([c, s, 0]), Vector.fromArray([-s, c, 0]),
            Vector.fromArray([0, 0, 1])],
        Vector.fromArray(extent));
}

// An independent test: the box vertices straddle (or touch) the plane.
function bruteForceIntersect(P: Hyperplane, B: OrientedBox): boolean {
    let numNeg = 0;
    let numPos = 0;
    for (const V of B.getVertices()) {
        const sd = dot(P.normal, V) - P.constant;
        if (sd < 0) {
            ++numNeg;
        }
        else if (sd > 0) {
            ++numPos;
        }
        else {
            return true;
        }
    }
    return numNeg > 0 && numPos > 0;
}

const ti = new IntrPlane3OrientedBox3TI();

describe('IntrPlane3OrientedBox3', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrPlane3OrientedBox3TIResult().intersect).toBe(false);
    });

    it('detects a plane cutting through an axis-aligned box', () => {
        const B = box([0, 0, 0], 0, [1, 1, 1]);
        expect(ti.test(plane([0, 0, 1], [0, 0, 0]), B).intersect).toBe(true);
        expect(ti.test(plane([0, 0, 1], [0, 0, 0.5]), B).intersect).toBe(true);
    });

    it('detects the tangent plane at a box face', () => {
        const B = box([0, 0, 0], 0, [1, 1, 1]);
        expect(ti.test(plane([0, 0, 1], [0, 0, 1]), B).intersect).toBe(true);
        expect(ti.test(plane([0, 0, 1], [0, 0, 1.0001]), B).intersect)
            .toBe(false);
    });

    it('accounts for the box orientation', () => {
        // A square rotated 45 degrees about z has half-width sqrt(2) in x.
        const B = box([0, 0, 0], Math.PI / 4, [1, 1, 1]);
        expect(ti.test(plane([1, 0, 0], [1.4, 0, 0]), B).intersect).toBe(true);
        expect(ti.test(plane([1, 0, 0], [1.5, 0, 0]), B).intersect).toBe(false);
        // The unrotated box would not reach x = 1.4.
        const A = box([0, 0, 0], 0, [1, 1, 1]);
        expect(ti.test(plane([1, 0, 0], [1.4, 0, 0]), A).intersect).toBe(false);
    });

    it('agrees with a vertex-sign test on random inputs', () => {
        let state = 13579;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand() * 2, rand() * 2, rand() * 2]);
            const B = box([rand() * 3, rand() * 3, rand() * 3],
                rand() * Math.PI,
                [0.2 + Math.abs(rand()), 0.2 + Math.abs(rand()),
                    0.2 + Math.abs(rand())]);
            const actual = ti.test(P, B).intersect;
            expect(actual).toBe(bruteForceIntersect(P, B));
            if (actual) {
                ++numHits;
            }
        }
        expect(numHits).toBeGreaterThan(50);
        expect(numHits).toBeLessThan(350);
    });
});
