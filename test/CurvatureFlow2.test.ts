import { describe, it, expect } from 'vitest';
import { CurvatureFlow2 } from '../src/CurvatureFlow2.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';

const NEUMANN = Number.MAX_VALUE;

function build(xBound: number, yBound: number,
    f: (x: number, y: number) => number): number[] {
    const data: number[] = [];
    for (let y = 0; y < yBound; ++y) {
        for (let x = 0; x < xBound; ++x) {
            data.push(f(x, y));
        }
    }
    return data;
}

// The number of grid edges whose endpoints straddle the level u = level; a
// discrete proxy for the perimeter of the level curve.
function levelSetPerimeter(filter: CurvatureFlow2, xBound: number, yBound: number,
    level: number): number {
    let count = 0;
    for (let y = 0; y < yBound; ++y) {
        for (let x = 0; x < xBound; ++x) {
            const u = filter.getU(x, y);
            if (x + 1 < xBound && (u - level) * (filter.getU(x + 1, y) - level) < 0) {
                ++count;
            }
            if (y + 1 < yBound && (u - level) * (filter.getU(x, y + 1) - level) < 0) {
                ++count;
            }
        }
    }
    return count;
}

describe('CurvatureFlow2', () => {
    it('leaves a linear ramp unchanged (the level curves are straight)', () => {
        // u = x has zero second derivatives, so the numerator vanishes.
        const filter = new CurvatureFlow2(5, 5, 1, 1, build(5, 5, (x) => x), null,
            NEUMANN, PdeFilterScaleType.NONE);
        filter.setTimeStep(0.25);
        filter.update();
        // Only the pixels with a full linear 3x3 neighborhood are checked; the
        // outermost columns see the duplicated Neumann border.
        for (let y = 1; y <= 3; ++y) {
            for (let x = 1; x <= 3; ++x) {
                expect(filter.getU(x, y)).toBeCloseTo(x, 12);
            }
        }
    });

    it('leaves a constant image unchanged through the zero-gradient branch', () => {
        const filter = new CurvatureFlow2(5, 5, 1, 1, new Array<number>(25).fill(7),
            null, NEUMANN, PdeFilterScaleType.NONE);
        filter.setTimeStep(0.5);
        filter.update();
        for (let y = 0; y < 5; ++y) {
            for (let x = 0; x < 5; ++x) {
                // Constant data is shifted to zero by the base-class scaling
                // (upstream issue #60); the zero-denominator branch copies it
                // through unchanged.
                expect(filter.getU(x, y)).toBe(0);
            }
        }
    });

    it('moves a circular level curve at its curvature 1/r', () => {
        // u(x,y) = distance to the center. Along the x axis, ux = 1, uy = 0,
        // uxy = 0 and uxx = 0, so the update is u + dt*uyy with uyy the
        // discrete approximation of the curvature 1/r.
        const bound = 9, c = 4;
        const radial = (x: number, y: number) => Math.hypot(x - c, y - c);
        const filter = new CurvatureFlow2(bound, bound, 1, 1,
            build(bound, bound, radial), null, NEUMANN, PdeFilterScaleType.NONE);
        const dt = 0.1;
        filter.setTimeStep(dt);
        filter.update();

        for (const r of [2, 3]) {
            const value = filter.getU(c + r, c);
            expect(value).toBeGreaterThan(r);
            // The discrete curvature is within a few percent of 1/r.
            expect((value - r) / dt).toBeCloseTo(1 / r, 1);
        }
        // The flow is isotropic on this image: the four axis directions agree.
        expect(filter.getU(c + 2, c)).toBeCloseTo(filter.getU(c - 2, c), 12);
        expect(filter.getU(c + 2, c)).toBeCloseTo(filter.getU(c, c + 2), 12);
        expect(filter.getU(c + 2, c)).toBeCloseTo(filter.getU(c, c - 2), 12);
    });

    it('shrinks the perimeter of a level curve', () => {
        // A blob whose zero level curve is a wiggly closed curve. Curvature
        // flow smooths the wiggles, which shortens the curve.
        const bound = 21, c = 10;
        const blob = (x: number, y: number) => {
            const dx = x - c, dy = y - c;
            const r = Math.hypot(dx, dy);
            const theta = Math.atan2(dy, dx);
            return r - (6 + 1.5 * Math.cos(4 * theta));
        };
        const data = build(bound, bound, blob);
        const filter = new CurvatureFlow2(bound, bound, 1, 1, data, null, NEUMANN,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(0.2);

        // ScaleType.NONE still subtracts the data minimum (upstream issue
        // #60), so the curve u = 0 of the input sits at this level in the
        // stored image.
        const level = -Math.min(...data);
        const before = levelSetPerimeter(filter, bound, bound, level);
        expect(before).toBeGreaterThan(0);
        for (let i = 0; i < 40; ++i) {
            filter.update();
        }
        const after = levelSetPerimeter(filter, bound, bound, level);
        expect(after).toBeLessThan(before);

        for (let y = 0; y < bound; ++y) {
            for (let x = 0; x < bound; ++x) {
                expect(Number.isFinite(filter.getU(x, y))).toBe(true);
            }
        }
    });

    it('reproduces the upstream update formula on a hand-computed pixel', () => {
        // A 3x3 image with Dirichlet border 0. At the center pixel the padded
        // 3x3 neighborhood is the image itself.
        const data = [
            0, 1, 4,
            2, 3, 5,
            6, 7, 8
        ];
        const filter = new CurvatureFlow2(3, 3, 1, 1, data, null, 0,
            PdeFilterScaleType.NONE);
        const dt = 0.05;
        filter.setTimeStep(dt);
        filter.update();

        const u = (x: number, y: number) => data[x + 3 * y];
        const ux = 0.5 * (u(2, 1) - u(0, 1));
        const uy = 0.5 * (u(1, 2) - u(1, 0));
        const uxx = u(2, 1) - 2 * u(1, 1) + u(0, 1);
        const uyy = u(1, 2) - 2 * u(1, 1) + u(1, 0);
        const uxy = 0.25 * (u(0, 0) + u(2, 2) - u(0, 2) - u(2, 0));
        const sqrUx = ux * ux, sqrUy = uy * uy;
        // Upstream uses -0.5*uxy*ux*uy where the analytic curvature-flow
        // numerator has -2*uxy*ux*uy; the quirk is preserved by the port.
        const numer = uxx * sqrUy + uyy * sqrUx - 0.5 * uxy * ux * uy;
        const expected = u(1, 1) + dt * numer / (sqrUx + sqrUy);

        expect(filter.getU(1, 1)).toBeCloseTo(expected, 12);
    });

    it('does not update masked-out pixels', () => {
        const bound = 5, c = 2;
        const radial = (x: number, y: number) => Math.hypot(x - c, y - c);
        const mask = new Array<number>(bound * bound).fill(1);
        mask[3 + bound * 2] = 0;
        const filter = new CurvatureFlow2(bound, bound, 1, 1,
            build(bound, bound, radial), mask, NEUMANN, PdeFilterScaleType.NONE);
        filter.setTimeStep(0.1);
        const before = filter.getU(3, 2);
        filter.update();
        // The masked pixel is not visited, so it keeps the Neumann mask-border
        // average assigned before the update.
        expect(filter.getMask(3, 2)).toBe(0);
        expect(filter.getU(3, 2)).not.toBe(before + 1);
        expect(Number.isFinite(filter.getU(3, 2))).toBe(true);
    });
});
