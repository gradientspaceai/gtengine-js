import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, sub, mul, dot, length } from '../src/Vector.js';
import { IntrCapsule3Capsule3TI } from '../src/IntrCapsule3Capsule3.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function capsule(p0: Vector, p1: Vector, radius: number): Capsule {
    return Capsule.fromSegmentRadius(Segment.fromEndpoints(p0, p1), radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Brute-force containment: is the point inside the solid capsule?
function inCapsule(p: Vector, c: Capsule): boolean {
    const p0 = c.segment.p[0];
    const p1 = c.segment.p[1];
    const dir = sub(p1, p0);
    const dd = dot(dir, dir);
    let t = dd > 0 ? dot(sub(p, p0), dir) / dd : 0;
    t = Math.max(0, Math.min(1, t));
    const closest = add(p0, mul(t, dir));
    return length(sub(p, closest)) <= c.radius;
}

describe('IntrCapsule3Capsule3', () => {
    const ti = new IntrCapsule3Capsule3TI();

    it('reports intersection for overlapping parallel capsules', () => {
        const c0 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 1);
        const c1 = capsule(v3(-1, 1.5, 0), v3(1, 1.5, 0), 1);
        expect(ti.test(c0, c1).intersect).toBe(true);
    });

    it('reports no intersection for separated parallel capsules', () => {
        const c0 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 1);
        const c1 = capsule(v3(-1, 2.5, 0), v3(1, 2.5, 0), 1);
        expect(ti.test(c0, c1).intersect).toBe(false);
    });

    it('reports intersection when the capsules are exactly tangent', () => {
        // The segment distance is exactly 2 = r0 + r1, so the closed solids
        // touch and the query uses <=.
        const c0 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 1);
        const c1 = capsule(v3(-1, 2, 0), v3(1, 2, 0), 1);
        expect(ti.test(c0, c1).intersect).toBe(true);
    });

    it('handles perpendicular crossing capsules', () => {
        const c0 = capsule(v3(-2, 0, 0), v3(2, 0, 0), 0.25);
        const c1 = capsule(v3(0, -2, 0.4), v3(0, 2, 0.4), 0.25);
        // The axes pass within 0.4 of each other, sum of radii is 0.5.
        expect(ti.test(c0, c1).intersect).toBe(true);

        const c2 = capsule(v3(0, -2, 0.6), v3(0, 2, 0.6), 0.25);
        expect(ti.test(c0, c2).intersect).toBe(false);
    });

    it('is symmetric in its arguments', () => {
        const rnd = makeRandom(12345);
        for (let k = 0; k < 200; ++k) {
            const c0 = capsule(
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rnd());
            const c1 = capsule(
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                rnd());
            expect(ti.test(c0, c1).intersect).toBe(ti.test(c1, c0).intersect);
        }
    });

    it('agrees with a brute-force sampling of the first capsule', () => {
        // Sampling the medial segment of capsule0 and offsetting by up to
        // radius0 produces points of capsule0. If any lies in capsule1, the
        // query must report an intersection.
        const rnd = makeRandom(777);
        let numFound = 0;
        for (let k = 0; k < 120; ++k) {
            const c0 = capsule(
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                0.3 + rnd() * 0.5);
            const c1 = capsule(
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                v3(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                0.3 + rnd() * 0.5);
            const intersect = ti.test(c0, c1).intersect;

            let sampleHit = false;
            for (let s = 0; s <= 20 && !sampleHit; ++s) {
                const base = add(c0.segment.p[0],
                    mul(s / 20, sub(c0.segment.p[1], c0.segment.p[0])));
                for (let a = 0; a < 12 && !sampleHit; ++a) {
                    for (let b = 0; b < 6 && !sampleHit; ++b) {
                        const theta = (2 * Math.PI * a) / 12;
                        const phi = (Math.PI * b) / 5;
                        const off = mul(c0.radius, v3(
                            Math.sin(phi) * Math.cos(theta),
                            Math.sin(phi) * Math.sin(theta),
                            Math.cos(phi)));
                        if (inCapsule(add(base, off), c1)) {
                            sampleHit = true;
                        }
                    }
                }
            }

            if (sampleHit) {
                expect(intersect).toBe(true);
                ++numFound;
            }
        }
        expect(numFound).toBeGreaterThan(0);
    });

    it('handles degenerate zero-length segments and zero radii', () => {
        // Two points, treated as spheres of the given radii.
        const p0 = capsule(v3(0, 0, 0), v3(0, 0, 0), 1);
        const p1 = capsule(v3(1.5, 0, 0), v3(1.5, 0, 0), 1);
        expect(ti.test(p0, p1).intersect).toBe(true);

        const p2 = capsule(v3(2.5, 0, 0), v3(2.5, 0, 0), 1);
        expect(ti.test(p0, p2).intersect).toBe(false);

        // Zero radius: the capsules degenerate to their segments.
        const s0 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 0);
        const s1 = capsule(v3(0, -1, 0), v3(0, 1, 0), 0);
        expect(ti.test(s0, s1).intersect).toBe(true);

        const s2 = capsule(v3(0, -1, 1), v3(0, 1, 1), 0);
        expect(ti.test(s0, s2).intersect).toBe(false);
    });

    it('throws when a capsule is not 3-dimensional', () => {
        const c3 = capsule(v3(-1, 0, 0), v3(1, 0, 0), 1);
        const c2 = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(Vector.fromArray([0, 0]),
                Vector.fromArray([1, 0])), 1);
        expect(() => ti.test(c3, c2)).toThrow();
    });
});
