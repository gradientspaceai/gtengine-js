import { describe, it, expect } from 'vitest';
import { CurvatureFlow3 } from '../src/CurvatureFlow3.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';

const NEUMANN = Number.MAX_VALUE;

function build(bound: number, f: (x: number, y: number, z: number) => number): number[] {
    const data: number[] = [];
    for (let z = 0; z < bound; ++z) {
        for (let y = 0; y < bound; ++y) {
            for (let x = 0; x < bound; ++x) {
                data.push(f(x, y, z));
            }
        }
    }
    return data;
}

describe('CurvatureFlow3', () => {
    it('leaves a linear ramp unchanged (the level surfaces are planes)', () => {
        const bound = 5;
        const filter = new CurvatureFlow3(bound, bound, bound, 1, 1, 1,
            build(bound, (x, y) => x + 2 * y), null, NEUMANN, PdeFilterScaleType.NONE);
        filter.setTimeStep(0.25);
        filter.update();
        for (let z = 1; z <= 3; ++z) {
            for (let y = 1; y <= 3; ++y) {
                for (let x = 1; x <= 3; ++x) {
                    expect(filter.getU(x, y, z)).toBeCloseTo(x + 2 * y, 12);
                }
            }
        }
    });

    it('leaves a constant image unchanged through the zero-gradient branch', () => {
        const filter = new CurvatureFlow3(5, 5, 5, 1, 1, 1,
            new Array<number>(125).fill(-4), null, NEUMANN, PdeFilterScaleType.NONE);
        filter.setTimeStep(0.5);
        filter.update();
        for (let z = 0; z < 5; ++z) {
            for (let y = 0; y < 5; ++y) {
                for (let x = 0; x < 5; ++x) {
                    // Constant data is shifted to zero by the base-class
                    // scaling (upstream issue #60).
                    expect(filter.getU(x, y, z)).toBe(0);
                }
            }
        }
    });

    it('moves a spherical level surface at twice its mean curvature 1/r', () => {
        // For u(x,y,z) = distance to the center, the numerator/denominator is
        // |grad(u)| * div(grad(u)/|grad(u)|) = 2/r on the sphere of radius r.
        const bound = 9, c = 4;
        const radial = (x: number, y: number, z: number) =>
            Math.sqrt((x - c) ** 2 + (y - c) ** 2 + (z - c) ** 2);
        const filter = new CurvatureFlow3(bound, bound, bound, 1, 1, 1,
            build(bound, radial), null, NEUMANN, PdeFilterScaleType.NONE);
        const dt = 0.05;
        filter.setTimeStep(dt);
        filter.update();

        for (const r of [2, 3]) {
            const value = filter.getU(c + r, c, c);
            expect(value).toBeGreaterThan(r);
            // The discrete estimate is within 10% of the exact 2/r.
            const speed = (value - r) / dt;
            expect(Math.abs(speed - 2 / r)).toBeLessThan(0.1 * (2 / r));
        }
        // Isotropy: the six axis directions agree.
        const reference = filter.getU(c + 2, c, c);
        expect(filter.getU(c - 2, c, c)).toBeCloseTo(reference, 12);
        expect(filter.getU(c, c + 2, c)).toBeCloseTo(reference, 12);
        expect(filter.getU(c, c - 2, c)).toBeCloseTo(reference, 12);
        expect(filter.getU(c, c, c + 2)).toBeCloseTo(reference, 12);
        expect(filter.getU(c, c, c - 2)).toBeCloseTo(reference, 12);
    });

    it('agrees with the 2D flow on an image that is constant along z', () => {
        // When u does not depend on z, uz = uxz = uyz = uzz = 0 and the 3D
        // numerator reduces to uy*(uxx*uy - uxy*ux) + ux*(uyy*ux - uxy*uy),
        // which is the 2D numerator with the mixed term weighted by 2 rather
        // than the 0.5 that CurvatureFlow2 uses.
        const bound = 7, c = 3;
        const f = (x: number, y: number) => Math.hypot(x - c, y - c);
        const filter = new CurvatureFlow3(bound, bound, bound, 1, 1, 1,
            build(bound, (x, y) => f(x, y)), null, NEUMANN, PdeFilterScaleType.NONE);
        const dt = 0.05;
        filter.setTimeStep(dt);
        filter.update();

        // The result is independent of z away from the duplicated border.
        for (let z = 1; z <= bound - 2; ++z) {
            expect(filter.getU(c + 2, c, z)).toBeCloseTo(filter.getU(c + 2, c, 1), 12);
        }
        // The cylinder of radius r moves at curvature 1/r.
        expect((filter.getU(c + 2, c, 3) - 2) / dt).toBeCloseTo(0.5, 1);
    });

    it('reproduces the upstream update formula on a hand-computed voxel', () => {
        const bound = 3;
        // A quadratic-ish sample with all 27 neighbors distinct.
        const value = (x: number, y: number, z: number) =>
            1 + x + 2 * y + 3 * z + 0.5 * x * y - 0.25 * y * z + 0.75 * x * z + x * x;
        const data = build(bound, value);
        const filter = new CurvatureFlow3(bound, bound, bound, 1, 1, 1, data, null, 0,
            PdeFilterScaleType.NONE);
        const dt = 0.02;
        filter.setTimeStep(dt);

        // Reproduce the constructor's ScaleType.NONE shift.
        const shift = -Math.min(...data);
        const u = (x: number, y: number, z: number) =>
            data[x + bound * (y + bound * z)] + shift;

        const ux = 0.5 * (u(2, 1, 1) - u(0, 1, 1));
        const uy = 0.5 * (u(1, 2, 1) - u(1, 0, 1));
        const uz = 0.5 * (u(1, 1, 2) - u(1, 1, 0));
        const uxx = u(2, 1, 1) - 2 * u(1, 1, 1) + u(0, 1, 1);
        const uyy = u(1, 2, 1) - 2 * u(1, 1, 1) + u(1, 0, 1);
        const uzz = u(1, 1, 2) - 2 * u(1, 1, 1) + u(1, 1, 0);
        const uxy = 0.25 * (u(0, 0, 1) + u(2, 2, 1) - u(2, 0, 1) - u(0, 2, 1));
        const uxz = 0.25 * (u(0, 1, 0) + u(2, 1, 2) - u(2, 1, 0) - u(0, 1, 2));
        const uyz = 0.25 * (u(1, 0, 0) + u(1, 2, 2) - u(1, 2, 0) - u(1, 0, 2));

        const denom = ux * ux + uy * uy + uz * uz;
        const numer0 = uy * (uxx * uy - uxy * ux) + ux * (uyy * ux - uxy * uy);
        const numer1 = uz * (uxx * uz - uxz * ux) + ux * (uzz * ux - uxz * uz);
        const numer2 = uz * (uyy * uz - uyz * uy) + uy * (uzz * uy - uyz * uz);
        const expected = u(1, 1, 1) + dt * (numer0 + numer1 + numer2) / denom;

        filter.update();
        expect(filter.getU(1, 1, 1)).toBeCloseTo(expected, 12);
    });

    it('produces finite values over many iterations of a blob', () => {
        const bound = 9, c = 4;
        const blob = (x: number, y: number, z: number) =>
            Math.sqrt((x - c) ** 2 + 1.5 * (y - c) ** 2 + 0.5 * (z - c) ** 2) - 3;
        const filter = new CurvatureFlow3(bound, bound, bound, 1, 1, 1,
            build(bound, blob), null, NEUMANN, PdeFilterScaleType.NONE);
        filter.setTimeStep(0.1);
        for (let i = 0; i < 10; ++i) {
            filter.update();
        }
        for (let z = 0; z < bound; ++z) {
            for (let y = 0; y < bound; ++y) {
                for (let x = 0; x < bound; ++x) {
                    expect(Number.isFinite(filter.getU(x, y, z))).toBe(true);
                }
            }
        }
    });
});
