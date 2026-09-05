import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Halfspace } from '../src/Halfspace.js';
import { Line } from '../src/Line.js';
import { Vector, add, mul, dot, normalize } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import { IntrHalfspace3Cylinder3TI } from '../src/IntrHalfspace3Cylinder3.js';

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

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, seededRandom
} from './helpers/arbitraries.js';
import { sub, length as vlength } from '../src/Vector.js';

const halfspaceV31 = fc.tuple(unitVector(3),
    fc.double({ min: -5, max: 5, noNaN: true }))
    .map(([n, c]) => Halfspace.fromNormalConstant(n, c));

const cylinderV31 = fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
    fc.double({ min: 0.1, max: 2, noNaN: true }),
    fc.double({ min: 0.2, max: 5, noNaN: true }))
    .map(([o, w, r, h]) => Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(o, w), r, h));

// The signed distance of the farthest cylinder point along the halfspace
// normal, computed from the geometry rather than from the query's algebra.
function supportSigned(h: Halfspace, c: Cylinder3): number {
    const ndw = dot(h.normal, c.axis.direction);
    const radial = Math.sqrt(Math.max(0, 1 - ndw * ndw));
    return dot(h.normal, c.axis.origin) - h.constant
        + c.radius * radial + 0.5 * c.height * Math.abs(ndw);
}

describe('IntrHalfspace3Cylinder3 verification', () => {
    const tiq = new IntrHalfspace3Cylinder3TI();

    it('equals the geometric support-point test', () => {
        check(fc.tuple(halfspaceV31, cylinderV31), ([h, c]) => {
            expect(tiq.test(h, c).intersect).toBe(supportSigned(h, c) >= 0);
        });
    });

    it('agrees with a dense sampling of the cylinder', () => {
        const rnd = seededRandom(0x30f7b1a2);
        for (let trial = 0; trial < 150; ++trial) {
            const n = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (vlength(n) < 0.3) {
                continue;
            }
            normalize(n);
            const h = Halfspace.fromNormalConstant(n, rnd() * 6 - 3);
            const w = Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1,
                rnd() * 2 - 1]);
            if (vlength(w) < 0.3) {
                continue;
            }
            normalize(w);
            const c = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(Vector.fromArray([rnd() * 6 - 3,
                    rnd() * 6 - 3, rnd() * 6 - 3]), w),
                0.2 + rnd() * 1.2, 0.5 + rnd() * 3);
            const basis = [w.clone(), new Vector(3), new Vector(3)];
            computeOrthogonalComplement3(1, basis);
            const half = 0.5 * c.height;
            let maxSigned = -Infinity;
            for (let i = 0; i <= 20; ++i) {
                const z = -half + (2 * half * i) / 20;
                for (let k = 0; k < 90; ++k) {
                    const a = (2 * Math.PI * k) / 90;
                    const p = add(c.axis.origin, add(mul(z, basis[0]),
                        add(mul(c.radius * Math.cos(a), basis[1]),
                            mul(c.radius * Math.sin(a), basis[2]))));
                    maxSigned = Math.max(maxSigned,
                        dot(h.normal, p) - h.constant);
                }
            }
            const got = tiq.test(h, c).intersect;
            if (maxSigned >= 0) {
                expect(got).toBe(true);
            }
            if (maxSigned < -0.2) {
                expect(got).toBe(false);
            }
        }
    }, 30000);

    it('the fixed clamp matters: an axis-aligned normal uses radial 0', () => {
        // Regression for the upstream 'std::max((T)1, 1 - absNdW^2)' typo.
        // With N == W the true radial term is sqrt(1 - 1) = 0, so the support
        // point is at distance h/2 above the axis origin. Placing the plane
        // between h/2 and h/2 + r separates them; the upstream expression
        // (which always yields root = 1) would report an intersection.
        const w = Vector.fromArray([0, 0, 1]);
        const c = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(Vector.fromArray([0, 0, 0]), w), 3, 2);
        // Support along +z is z = 1; the upstream typo would give 1 + 3 = 4.
        const h = Halfspace.fromNormalConstant(w, 2.5);
        expect(supportSigned(h, c)).toBeLessThan(0);
        expect(tiq.test(h, c).intersect).toBe(false);
        // Just below the true support the query must report an intersection.
        expect(tiq.test(Halfspace.fromNormalConstant(w, 0.9), c).intersect)
            .toBe(true);
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(halfspaceV31, cylinderV31, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([h, c, R, tr]) => {
            const rot = (p: Vector): Vector => add(mul(p.values[0], R[0]),
                add(mul(p.values[1], R[1]), mul(p.values[2], R[2])));
            const n2 = rot(h.normal);
            const h2 = Halfspace.fromNormalConstant(n2,
                h.constant + dot(n2, tr));
            const c2 = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(add(tr, rot(c.axis.origin)),
                    rot(c.axis.direction)), c.radius, c.height);
            if (Math.abs(supportSigned(h, c)) < 1e-9) {
                return;   // exactly tangent
            }
            expect(tiq.test(h, c).intersect).toBe(tiq.test(h2, c2).intersect);
        });
    });

    it('a zero-radius zero-height cylinder is the axis origin', () => {
        check(fc.tuple(halfspaceV31, wellScaledVector(3, -4, 4),
            unitVector(3)), ([h, o, w]) => {
            const c = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(o, w), 0, 0);
            expect(tiq.test(h, c).intersect)
                .toBe(dot(h.normal, o) - h.constant >= 0);
            expect(sub(o, o).values).toEqual([0, 0, 0]);
        });
    });
});
