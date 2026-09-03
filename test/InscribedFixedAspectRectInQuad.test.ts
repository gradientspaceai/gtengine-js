import { describe, it, expect } from 'vitest';
import { InscribedFixedAspectRectInQuad } from '../src/InscribedFixedAspectRectInQuad.js';
import { Vector } from '../src/Vector.js';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Test whether the point is inside or on the counterclockwise convex quad,
// allowing a small tolerance for floating-point roundoff.
function insideQuad(quad: Vector[], x: number, y: number, epsilon = 1e-9): boolean {
    for (let i = 0; i < 4; ++i) {
        const p0 = quad[i];
        const p1 = quad[(i + 1) % 4];
        const ex = p1.values[0] - p0.values[0];
        const ey = p1.values[1] - p0.values[1];
        const px = x - p0.values[0];
        const py = y - p0.values[1];
        if (ex * py - ey * px < -epsilon) {
            return false;
        }
    }
    return true;
}

function rectCorners(origin: Vector, width: number, height: number):
    Array<[number, number]> {
    const u = origin.values[0];
    const v = origin.values[1];
    return [[u, v], [u + width, v], [u + width, v + height], [u, v + height]];
}

describe('InscribedFixedAspectRectInQuad', () => {
    it('fills an axis-aligned rectangle when the aspect ratios match', () => {
        // The quad is the rectangle [0,4]x[0,2] whose aspect ratio is 2, so
        // the inscribed rectangle is the quad itself.
        const quad = [v2(0, 0), v2(4, 0), v2(4, 2), v2(0, 2)];
        const result = InscribedFixedAspectRectInQuad.execute(quad, 2);
        expect(result.isUnique).toBe(true);
        expect(result.rectWidth).toBeCloseTo(4, 12);
        expect(result.rectHeight).toBeCloseTo(2, 12);
        expect(result.rectOrigin.values[0]).toBeCloseTo(0, 12);
        expect(result.rectOrigin.values[1]).toBeCloseTo(0, 12);
    });

    it('reports infinitely many solutions when the rectangle can slide', () => {
        // The largest square inside [0,4]x[0,2] is 2x2 and can be translated
        // horizontally, so the solution is not unique.
        const quad = [v2(0, 0), v2(4, 0), v2(4, 2), v2(0, 2)];
        const result = InscribedFixedAspectRectInQuad.execute(quad, 1);
        expect(result.isUnique).toBe(false);
        expect(result.rectWidth).toBeCloseTo(2, 12);
        expect(result.rectHeight).toBeCloseTo(2, 12);
        expect(result.rectOrigin.values[1]).toBeCloseTo(0, 12);
    });

    it('inscribes the maximum square in a symmetric trapezoid', () => {
        // The trapezoid has the bottom edge from (0,0) to (4,0) and the top
        // edge from (3,2) to (1,2). The 2x2 square with origin (1,0) has its
        // top corners exactly at the top vertices, so it is the maximum.
        const quad = [v2(0, 0), v2(4, 0), v2(3, 2), v2(1, 2)];
        const result = InscribedFixedAspectRectInQuad.execute(quad, 1);
        expect(result.isUnique).toBe(true);
        expect(result.rectWidth).toBeCloseTo(2, 10);
        expect(result.rectHeight).toBeCloseTo(2, 10);
        expect(result.rectOrigin.values[0]).toBeCloseTo(1, 10);
        expect(result.rectOrigin.values[1]).toBeCloseTo(0, 10);
    });

    it('honors the aspect ratio w/h', () => {
        const quad = [v2(0, 0), v2(6, 0), v2(6, 6), v2(0, 6)];
        for (const aspectRatio of [0.5, 1, 2, 3]) {
            const result = InscribedFixedAspectRectInQuad.execute(quad, aspectRatio);
            expect(result.rectWidth / result.rectHeight).toBeCloseTo(aspectRatio, 10);
            // The limiting dimension is 6.
            expect(Math.max(result.rectWidth, result.rectHeight)).toBeCloseTo(6, 10);
        }
    });

    it('is translation equivariant', () => {
        const quad = [v2(0, 0), v2(5, 1), v2(4, 4), v2(-1, 3)];
        const shift = v2(10, -7);
        const shifted = quad.map(p =>
            v2(p.values[0] + shift.values[0], p.values[1] + shift.values[1]));
        const base = InscribedFixedAspectRectInQuad.execute(quad, 1.5);
        const moved = InscribedFixedAspectRectInQuad.execute(shifted, 1.5);
        expect(moved.rectWidth).toBeCloseTo(base.rectWidth, 9);
        expect(moved.rectHeight).toBeCloseTo(base.rectHeight, 9);
        expect(moved.rectOrigin.values[0])
            .toBeCloseTo(base.rectOrigin.values[0] + shift.values[0], 9);
        expect(moved.rectOrigin.values[1])
            .toBeCloseTo(base.rectOrigin.values[1] + shift.values[1], 9);
    });

    it('scales the rectangle with a uniform scaling of the quad', () => {
        const quad = [v2(0, 0), v2(5, 1), v2(4, 4), v2(-1, 3)];
        const scale = 3;
        const scaled = quad.map(p => v2(scale * p.values[0], scale * p.values[1]));
        const base = InscribedFixedAspectRectInQuad.execute(quad, 2);
        const bigger = InscribedFixedAspectRectInQuad.execute(scaled, 2);
        expect(bigger.rectWidth).toBeCloseTo(scale * base.rectWidth, 9);
        expect(bigger.rectHeight).toBeCloseTo(scale * base.rectHeight, 9);
        expect(bigger.rectOrigin.values[0])
            .toBeCloseTo(scale * base.rectOrigin.values[0], 9);
    });

    it('throws when the quad does not have 4 vertices', () => {
        expect(() => InscribedFixedAspectRectInQuad.execute(
            [v2(0, 0), v2(1, 0), v2(1, 1)], 1)).toThrow();
    });

    it('produces an inscribed rectangle that no larger one can beat', () => {
        const random = makeRandom(160934);
        for (let trial = 0; trial < 100; ++trial) {
            // Build a convex counterclockwise quadrilateral by perturbing the
            // corners of the unit square outward.
            const quad = [
                v2(-1 - random(), -1 - random()),
                v2(1 + random(), -1 - random()),
                v2(1 + random(), 1 + random()),
                v2(-1 - random(), 1 + random())
            ];
            const aspectRatio = 0.4 + 2 * random();
            const result = InscribedFixedAspectRectInQuad.execute(quad, aspectRatio);

            expect(result.rectWidth).toBeGreaterThan(0);
            expect(result.rectHeight).toBeCloseTo(result.rectWidth / aspectRatio, 9);

            // The rectangle is inscribed: all four corners are in the quad.
            for (const [x, y] of rectCorners(result.rectOrigin, result.rectWidth,
                result.rectHeight)) {
                expect(insideQuad(quad, x, y, 1e-8)).toBe(true);
            }

            // Growing the rectangle by 1% about its origin pushes at least one
            // corner outside, so the returned width is maximal.
            const grownWidth = 1.01 * result.rectWidth;
            const grownHeight = grownWidth / aspectRatio;
            let allInside = true;
            for (const [x, y] of rectCorners(result.rectOrigin, grownWidth, grownHeight)) {
                if (!insideQuad(quad, x, y, 1e-8)) {
                    allInside = false;
                    break;
                }
            }
            expect(allInside).toBe(false);
        }
    });

    it('never beats the maximum width found by a brute-force search', () => {
        const random = makeRandom(6022140);
        for (let trial = 0; trial < 20; ++trial) {
            const quad = [
                v2(-1 - random(), -1 - random()),
                v2(1 + random(), -1 - random()),
                v2(1 + random(), 1 + random()),
                v2(-1 - random(), 1 + random())
            ];
            const aspectRatio = 0.5 + random();
            const result = InscribedFixedAspectRectInQuad.execute(quad, aspectRatio);

            // Sample origins on a grid and, for each, bisect on the width to
            // find the largest inscribed rectangle with that origin.
            let best = 0;
            const lo = -2.5;
            const hi = 2.5;
            const steps = 40;
            for (let i = 0; i <= steps; ++i) {
                const u = lo + (hi - lo) * i / steps;
                for (let j = 0; j <= steps; ++j) {
                    const v = lo + (hi - lo) * j / steps;
                    if (!insideQuad(quad, u, v)) {
                        continue;
                    }
                    let low = 0;
                    let high = 6;
                    for (let k = 0; k < 40; ++k) {
                        const mid = 0.5 * (low + high);
                        const corners = rectCorners(v2(u, v), mid, mid / aspectRatio);
                        const fits = corners.every(([x, y]) => insideQuad(quad, x, y));
                        if (fits) {
                            low = mid;
                        }
                        else {
                            high = mid;
                        }
                    }
                    best = Math.max(best, low);
                }
            }
            // The grid search is coarse, so it can only give a lower bound on
            // the optimum. The query result must be at least as large.
            expect(result.rectWidth).toBeGreaterThan(best - 1e-6);
        }
    });
});
