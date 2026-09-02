import { describe, expect, it } from 'vitest';
import { approximateEllipseByArcs } from '../src/ApprEllipseByArcs';
import { Vector, length, sub } from '../src/Vector';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// The implicit ellipse function value; it is 1 on the ellipse.
function implicit(a: number, b: number, p: Vector): number {
    const x = p.values[0] / a;
    const y = p.values[1] / b;
    return x * x + y * y;
}

// The maximum deviation of the arcs from the ellipse, measured by the
// implicit ellipse function at samples along each arc.
function maxArcDeviation(a: number, b: number, numArcs: number,
    samplesPerArc: number): number {
    const result = approximateEllipseByArcs(a, b, numArcs);
    expect(result).not.toBeNull();
    const { points, centers, radii } = result!;

    let maxError = 0;
    for (let i = 0; i < numArcs; ++i) {
        const c = centers[i];
        const d0 = sub(points[i], c);
        const d1 = sub(points[i + 1], c);
        const angle0 = Math.atan2(d0.values[1], d0.values[0]);
        let angle1 = Math.atan2(d1.values[1], d1.values[0]);
        // The arc is traversed counterclockwise from points[i] to
        // points[i+1].
        while (angle1 < angle0) {
            angle1 += 2 * Math.PI;
        }
        for (let j = 0; j <= samplesPerArc; ++j) {
            const t = j / samplesPerArc;
            const angle = angle0 + t * (angle1 - angle0);
            const q = vec(c.values[0] + radii[i] * Math.cos(angle),
                c.values[1] + radii[i] * Math.sin(angle));
            maxError = Math.max(maxError, Math.abs(implicit(a, b, q) - 1));
        }
    }
    return maxError;
}

describe('approximateEllipseByArcs degenerate inputs', () => {
    it('rejects fewer than two arcs and circles', () => {
        expect(approximateEllipseByArcs(2, 1, 1)).toBeNull();
        expect(approximateEllipseByArcs(2, 1, 0)).toBeNull();
        expect(approximateEllipseByArcs(2, 1, -3)).toBeNull();
        expect(approximateEllipseByArcs(1, 1, 4)).toBeNull();
        expect(approximateEllipseByArcs(2.5, 2.5, 8)).toBeNull();
    });
});

describe('approximateEllipseByArcs sample points', () => {
    it('produces the requested container sizes and the axis endpoints', () => {
        const a = 3, b = 1;
        for (const numArcs of [2, 3, 5, 16]) {
            const result = approximateEllipseByArcs(a, b, numArcs);
            expect(result).not.toBeNull();
            const { points, centers, radii } = result!;
            expect(points.length).toBe(numArcs + 1);
            expect(centers.length).toBe(numArcs);
            expect(radii.length).toBe(numArcs);

            expect(points[0].values[0]).toBeCloseTo(a, 14);
            expect(points[0].values[1]).toBeCloseTo(0, 14);
            expect(points[numArcs].values[0]).toBeCloseTo(0, 14);
            expect(points[numArcs].values[1]).toBeCloseTo(b, 14);
        }
    });

    it('places every sample point on the ellipse in the first quadrant', () => {
        for (const [a, b] of [[3, 1], [1, 4], [10, 9.5], [0.5, 2]]) {
            const numArcs = 7;
            const result = approximateEllipseByArcs(a, b, numArcs);
            expect(result).not.toBeNull();
            const points = result!.points;
            for (let i = 0; i <= numArcs; ++i) {
                expect(implicit(a, b, points[i])).toBeCloseTo(1, 12);
                expect(points[i].values[0]).toBeGreaterThanOrEqual(0);
                expect(points[i].values[1]).toBeGreaterThanOrEqual(0);
            }
            // The points are generated counterclockwise from (a,0) to (0,b),
            // so x strictly decreases and y strictly increases.
            for (let i = 1; i <= numArcs; ++i) {
                expect(points[i].values[0])
                    .toBeLessThan(points[i - 1].values[0]);
                expect(points[i].values[1])
                    .toBeGreaterThan(points[i - 1].values[1]);
            }
        }
    });
});

describe('approximateEllipseByArcs arcs', () => {
    it('interpolates the arc endpoints (a C0-continuous chain)', () => {
        const a = 4, b = 1.5, numArcs = 6;
        const result = approximateEllipseByArcs(a, b, numArcs);
        expect(result).not.toBeNull();
        const { points, centers, radii } = result!;
        for (let i = 0; i < numArcs; ++i) {
            expect(radii[i]).toBeGreaterThan(0);
            // Both endpoints of arc i lie on the circle (centers[i],
            // radii[i]), so consecutive arcs share an endpoint exactly.
            expect(length(sub(points[i], centers[i])))
                .toBeCloseTo(radii[i], 10);
            expect(length(sub(points[i + 1], centers[i])))
                .toBeCloseTo(radii[i], 10);
        }
    });

    it('makes the end arcs tangent to the axes by symmetry', () => {
        const a = 2, b = 1, numArcs = 5;
        const result = approximateEllipseByArcs(a, b, numArcs);
        expect(result).not.toBeNull();
        const { centers, radii } = result!;
        // The first arc is fit through the reflection of points[1] about the
        // x-axis, so its center is on the x-axis (tangency at (a,0)).
        expect(centers[0].values[1]).toBeCloseTo(0, 12);
        expect(centers[0].values[0]).toBeCloseTo(a - radii[0], 12);
        // The last arc is fit through the reflection of points[n-1] about the
        // y-axis, so its center is on the y-axis (tangency at (0,b)).
        expect(centers[numArcs - 1].values[0]).toBeCloseTo(0, 12);
        expect(centers[numArcs - 1].values[1])
            .toBeCloseTo(b - radii[numArcs - 1], 12);

        // The osculating circle radii at the axis endpoints are b^2/a and
        // a^2/b. The end arcs are close to (but not exactly) osculating.
        expect(radii[0]).toBeGreaterThan(0);
        expect(radii[0]).toBeLessThan(a);
        expect(radii[numArcs - 1]).toBeGreaterThan(b);
    });

    it('approximates the ellipse with an error that decreases with numArcs', () => {
        const a = 2, b = 1;
        const e2 = maxArcDeviation(a, b, 2, 64);
        const e4 = maxArcDeviation(a, b, 4, 64);
        const e8 = maxArcDeviation(a, b, 8, 64);
        const e16 = maxArcDeviation(a, b, 16, 64);
        expect(e4).toBeLessThan(e2);
        expect(e8).toBeLessThan(e4);
        expect(e16).toBeLessThan(e8);
        // Measured accuracy for a=2, b=1: 8 arcs track the implicit function to
        // better than 1.5%, and 16 arcs to better than 0.5%.
        expect(e8).toBeLessThan(1.5e-2);
        expect(e16).toBeLessThan(5e-3);
    });

    it('handles the tall ellipse (a < b) as well as the wide one', () => {
        // Swapping (a,b) mirrors the configuration about the line y = x.
        const wide = approximateEllipseByArcs(3, 1, 6);
        const tall = approximateEllipseByArcs(1, 3, 6);
        expect(wide).not.toBeNull();
        expect(tall).not.toBeNull();
        const n = 6;
        for (let i = 0; i <= n; ++i) {
            // points are ordered from (a,0) to (0,b) in both cases, so the
            // tall ellipse's point i corresponds to the wide ellipse's point
            // n-i with the coordinates swapped.
            expect(tall!.points[i].values[0])
                .toBeCloseTo(wide!.points[n - i].values[1], 10);
            expect(tall!.points[i].values[1])
                .toBeCloseTo(wide!.points[n - i].values[0], 10);
        }
        expect(maxArcDeviation(1, 3, 8, 64)).toBeLessThan(6e-2);
    });
});
