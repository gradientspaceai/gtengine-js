import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import { IntrHalfspace3Sphere3TI } from '../src/IntrHalfspace3Sphere3.js';
import { length } from '../src/Vector.js';
import { check, fc, positive, rotationFrame, unitVector, wellScaled } from './helpers/arbitraries.js';

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

describe('IntrHalfspace3Sphere3 verification', () => {
    const ti = new IntrHalfspace3Sphere3TI();

    const hsArb = fc.tuple(unitVector(3), wellScaled(-6, 6))
        .map(([n, c]) => Halfspace.fromNormalConstant(n, c));
    const sphereArb = fc.tuple(
        fc.array(wellScaled(-6, 6), { minLength: 3, maxLength: 3 }),
        positive(4)).map(([c, r]) => sphere(c[0], c[1], c[2], r));

    it('a sampled ball point inside the halfspace forces intersect = true', () => {
        const rnd = makeRandom(0x7f4a13);
        check(fc.tuple(hsArb, sphereArb), ([h, s]) => {
            const hit = ti.test(h, s).intersect;
            let found = false;
            for (let k = 0; k < 300; ++k) {
                // A uniform-ish direction; only membership matters.
                let dx = 2 * rnd() - 1, dy = 2 * rnd() - 1, dz = 2 * rnd() - 1;
                const len = Math.hypot(dx, dy, dz);
                if (len < 1e-6) {
                    continue;
                }
                const t = s.radius * rnd() / len;
                dx *= t; dy *= t; dz *= t;
                const p = Vector.fromArray([s.center.values[0] + dx,
                    s.center.values[1] + dy, s.center.values[2] + dz]);
                if (dot(h.normal, p) - h.constant >= 0) {
                    found = true;
                    break;
                }
            }
            if (found) {
                expect(hit).toBe(true);
            }
        }, 60);
    }, 30000);

    it('intersect is decided by the extreme point C + r*N', () => {
        check(fc.tuple(hsArb, sphereArb), ([h, s]) => {
            // The support point of the ball in the normal direction is the
            // one that maximizes dot(N, X); the ball meets the halfspace iff
            // that point does. This builds the witness geometrically instead
            // of repeating the query's algebra.
            const p = Vector.fromArray([
                s.center.values[0] + s.radius * h.normal.values[0],
                s.center.values[1] + s.radius * h.normal.values[1],
                s.center.values[2] + s.radius * h.normal.values[2]]);
            const f = dot(h.normal, p) - h.constant;
            if (Math.abs(f) < 1e-12 * (1 + Math.abs(h.constant) + s.radius)) {
                return;    // exactly tangent: sign is not robust
            }
            expect(ti.test(h, s).intersect).toBe(f >= 0);
        });
    });

    it('is equivariant under a rigid motion of the halfspace and sphere', () => {
        check(fc.tuple(hsArb, sphereArb, rotationFrame(3),
            fc.array(wellScaled(-4, 4), { minLength: 3, maxLength: 3 })),
            ([h, s, R, t]) => {
                const apply = (v: Vector): number[] => {
                    const out = [0, 0, 0];
                    for (let d = 0; d < 3; ++d) {
                        out[d] = R[0].values[d] * v.values[0]
                            + R[1].values[d] * v.values[1]
                            + R[2].values[d] * v.values[2];
                    }
                    return out;
                };
                const n2 = apply(h.normal);
                const c2 = apply(s.center);
                // The plane dot(N,X) = c maps to dot(R*N, Y) = c + dot(R*N,t).
                const nv = Vector.fromArray(n2);
                const constant = h.constant + (n2[0] * t[0] + n2[1] * t[1]
                    + n2[2] * t[2]);
                const h2 = Halfspace.fromNormalConstant(nv, constant);
                const s2 = sphere(c2[0] + t[0], c2[1] + t[1], c2[2] + t[2],
                    s.radius);
                const f = dot(h.normal, s.center) - h.constant + s.radius;
                if (Math.abs(f) < 1e-9 * (1 + Math.abs(h.constant)
                    + s.radius + length(s.center))) {
                    return;    // near tangency, the rotated arithmetic differs
                }
                expect(ti.test(h2, s2).intersect).toBe(ti.test(h, s).intersect);
            });
    });

    it('a zero-radius sphere is the point test on its center', () => {
        const h = halfspace(0, 0, 1, 2);   // z >= 2
        expect(ti.test(h, sphere(0, 0, 3, 0)).intersect).toBe(true);
        expect(ti.test(h, sphere(0, 0, 2, 0)).intersect).toBe(true);  // on plane
        expect(ti.test(h, sphere(0, 0, 1.5, 0)).intersect).toBe(false);
    });
});
