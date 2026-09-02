import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule';
import { DistPointHyperplane } from '../src/DistPointHyperplane';
import { Hyperplane } from '../src/Hyperplane';
import {
    IntrPlane3Capsule3TI,
    defaultIntrPlane3Capsule3TIResult
} from '../src/IntrPlane3Capsule3';
import { Segment } from '../src/Segment';
import { Vector, add, mul, normalize, sub } from '../src/Vector';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

function capsule(p0: number[], p1: number[], radius: number): Capsule {
    return Capsule.fromSegmentRadius(
        Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1)),
        radius);
}

// An independent test: sample the capsule segment densely and check whether
// any sample is within 'radius' of the plane.
function bruteForceIntersect(P: Hyperplane, C: Capsule): boolean {
    const query = new DistPointHyperplane();
    const p0 = C.segment.p[0];
    const delta = sub(C.segment.p[1], p0);
    const n = 2000;
    for (let i = 0; i <= n; ++i) {
        const X = add(p0, mul(i / n, delta));
        if (query.compute(X, P).distance <= C.radius) {
            return true;
        }
    }
    return false;
}

const ti = new IntrPlane3Capsule3TI();

describe('IntrPlane3Capsule3', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrPlane3Capsule3TIResult().intersect).toBe(false);
    });

    it('detects endpoints on opposite sides of the plane', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        expect(ti.test(P, capsule([0, 0, -2], [0, 0, 3], 0.25)).intersect)
            .toBe(true);
    });

    it('detects an endpoint exactly on the plane', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        expect(ti.test(P, capsule([0, 0, 0], [0, 0, 3], 0.25)).intersect)
            .toBe(true);
    });

    it('detects a capsule that reaches the plane only with its end sphere', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        // Both endpoints are above the plane; the closer one is at z = 0.5.
        expect(ti.test(P, capsule([0, 0, 0.5], [0, 0, 3], 1)).intersect)
            .toBe(true);
        expect(ti.test(P, capsule([0, 0, 0.5], [0, 0, 3], 0.25)).intersect)
            .toBe(false);
    });

    it('reports the tangent configuration (distance equals radius)', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        expect(ti.test(P, capsule([0, 0, 1], [0, 0, 3], 1)).intersect)
            .toBe(true);
    });

    it('agrees with a dense sampling of the capsule on random inputs', () => {
        let state = 424242;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand(), rand(), rand()]);
            const C = capsule([rand() * 3, rand() * 3, rand() * 3],
                [rand() * 3, rand() * 3, rand() * 3],
                0.2 + Math.abs(rand()));
            const actual = ti.test(P, C).intersect;
            expect(actual).toBe(bruteForceIntersect(P, C));
            if (actual) {
                ++numHits;
            }
        }
        expect(numHits).toBeGreaterThan(20);
        expect(numHits).toBeLessThan(130);
    });
});
