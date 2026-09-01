import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3';
import { Halfspace } from '../src/Halfspace';
import { Line } from '../src/Line';
import { Vector, add, mul, dot, normalize } from '../src/Vector';
import { computeOrthogonalComplement3 } from '../src/Vector3';
import { IntrHalfspace3Cylinder3TI } from '../src/IntrHalfspace3Cylinder3';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function cylinder(c: Vector, d: Vector, r: number, h: number): Cylinder3 {
    const dir = d.clone();
    normalize(dir);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(c, dir), r, h);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The exact maximum of Dot(N,X)-d over the solid cylinder.
function exactMax(hs: Halfspace, cyl: Cylinder3): number {
    const center = dot(hs.normal, cyl.axis.origin) - hs.constant;
    const ndw = dot(hs.normal, cyl.axis.direction);
    const absNdW = Math.abs(ndw);
    return center + cyl.radius * Math.sqrt(Math.max(0, 1 - absNdW * absNdW))
        + 0.5 * cyl.height * absNdW;
}

describe('IntrHalfspace3Cylinder3', () => {
    const ti = new IntrHalfspace3Cylinder3TI();
    const upper = Halfspace.fromNormalConstant(v3(0, 0, 1), 0);

    it('reports intersection when the cylinder is inside the halfspace', () => {
        expect(ti.test(upper, cylinder(v3(0, 0, 4), v3(0, 0, 1), 1, 2))
            .intersect).toBe(true);
    });

    it('reports no intersection when the cylinder is strictly below', () => {
        expect(ti.test(upper, cylinder(v3(0, 0, -4), v3(0, 0, 1), 1, 2))
            .intersect).toBe(false);
    });

    it('reports intersection for exact tangency (closed halfspace)', () => {
        // The cap of the cylinder is exactly on the plane z = 0.
        expect(ti.test(upper, cylinder(v3(0, 0, -1), v3(0, 0, 1), 1, 2))
            .intersect).toBe(true);
    });

    it('uses the corrected sqrt(1 - Dot(N,W)^2) clamp', () => {
        // Regression test for the upstream std::max((T)1, ...) typo, which
        // forces root = 1 and inflates the projection interval by the full
        // radius. Here the axis is parallel to the halfspace normal, so the
        // radius contributes nothing: the cylinder spans z in [-2.5,-0.5],
        // entirely below z = 0.
        const cyl = cylinder(v3(0, 0, -1.5), v3(0, 0, 1), 1, 2);
        expect(exactMax(upper, cyl)).toBeCloseTo(-0.5, 12);
        expect(ti.test(upper, cyl).intersect).toBe(false);
    });

    it('accounts for the radius when the axis is parallel to the plane', () => {
        // The axis lies at z = -0.5 and is orthogonal to the normal, so the
        // wall reaches z = 0.5.
        expect(ti.test(upper, cylinder(v3(0, 0, -0.5), v3(1, 0, 0), 1, 4))
            .intersect).toBe(true);
        expect(ti.test(upper, cylinder(v3(0, 0, -1.5), v3(1, 0, 0), 1, 4))
            .intersect).toBe(false);
    });

    it('agrees with a brute-force sampling of the cylinder surface', () => {
        const rnd = makeRandom(4242);
        let numTrue = 0, numFalse = 0;
        for (let k = 0; k < 300; ++k) {
            const n = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(n, n) < 1e-6) {
                continue;
            }
            normalize(n);
            const hs = Halfspace.fromNormalConstant(n, rnd() * 4 - 2);

            const axis = v3(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(axis, axis) < 1e-6) {
                continue;
            }
            const cyl = cylinder(
                v3(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                axis, 0.2 + rnd(), 0.5 + rnd() * 3);
            const intersect = ti.test(hs, cyl).intersect;

            // Sample the boundary of the cylinder and take the maximum of
            // Dot(N,X)-d. It must not exceed the analytic maximum, and a
            // strictly positive sample forces intersect = true.
            const basis = [cyl.axis.direction.clone(), new Vector(3),
                new Vector(3)];
            computeOrthogonalComplement3(1, basis);
            let maxSample = -Infinity;
            for (let a = 0; a < 40; ++a) {
                const theta = (2 * Math.PI * a) / 40;
                const radial = add(mul(Math.cos(theta), basis[1]),
                    mul(Math.sin(theta), basis[2]));
                for (const s of [-0.5, 0, 0.5]) {
                    for (const rr of [0, cyl.radius]) {
                        const p = add(add(cyl.axis.origin,
                            mul(s * cyl.height, cyl.axis.direction)),
                            mul(rr, radial));
                        maxSample = Math.max(maxSample,
                            dot(n, p) - hs.constant);
                    }
                }
            }

            expect(maxSample).toBeLessThanOrEqual(exactMax(hs, cyl) + 1e-9);
            if (maxSample > 1e-9) {
                expect(intersect).toBe(true);
                ++numTrue;
            }
            if (exactMax(hs, cyl) < -1e-9) {
                expect(intersect).toBe(false);
                ++numFalse;
            }
        }
        expect(numTrue).toBeGreaterThan(0);
        expect(numFalse).toBeGreaterThan(0);
    });

    it('handles a zero-radius, zero-height cylinder as a point', () => {
        expect(ti.test(upper, cylinder(v3(0, 0, 0), v3(0, 0, 1), 0, 0))
            .intersect).toBe(true);
        expect(ti.test(upper, cylinder(v3(0, 0, -1), v3(0, 0, 1), 0, 0))
            .intersect).toBe(false);
    });

    it('throws for an infinite cylinder', () => {
        const cyl = cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2);
        cyl.makeInfiniteCylinder();
        expect(() => ti.test(upper, cyl)).toThrow();
    });

    it('throws when the halfspace is not 3-dimensional', () => {
        const hs2 = Halfspace.fromNormalConstant(Vector.fromArray([0, 1]), 0);
        expect(() => ti.test(hs2, cylinder(v3(0, 0, 0), v3(0, 0, 1), 1, 2)))
            .toThrow();
    });
});
