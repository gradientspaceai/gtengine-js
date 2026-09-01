import { describe, expect, it } from 'vitest';
import { Arc2 } from '../src/Arc2';
import { DistPoint2Arc2 } from '../src/DistPoint2Arc2';
import { Vector, dot, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// An arc of the unit circle centered at C from angle a0 to angle a1
// (counterclockwise).
function arc(center: number[], radius: number, a0: number,
    a1: number): Arc2 {
    const c = v(...center);
    return Arc2.fromCenterRadiusEnds(c, radius,
        v(c.values[0] + radius * Math.cos(a0),
            c.values[1] + radius * Math.sin(a0)),
        v(c.values[0] + radius * Math.cos(a1),
            c.values[1] + radius * Math.sin(a1)));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistPoint2Arc2', () => {
    const query = new DistPoint2Arc2();

    it('uses the circle closest point when it lies on the arc', () => {
        // The upper-right quarter arc of the unit circle.
        const a = arc([0, 0], 1, 0, Math.PI / 2);
        const result = query.compute(v(3, 3), a);
        const s = Math.SQRT1_2;
        expect(result.equidistant).toBe(false);
        expect(result.closest[1].values[0]).toBeCloseTo(s, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(s, 12);
        expect(result.distance).toBeCloseTo(3 * Math.SQRT2 - 1, 12);
    });

    it('falls back to the nearer endpoint when the circle point is off arc',
        () => {
            const a = arc([0, 0], 1, 0, Math.PI / 2);
            // The point is below the arc; the closest circle point is (0,-1),
            // which is not on the arc, so the closest arc point is (1,0).
            const result = query.compute(v(0, -5), a);
            expect(result.equidistant).toBe(false);
            expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
            expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
            expect(result.distance).toBeCloseTo(Math.sqrt(1 + 25), 12);
        });

    it('selects the other endpoint when it is nearer', () => {
        const a = arc([0, 0], 1, 0, Math.PI / 2);
        const result = query.compute(v(-5, 0), a);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(25 + 1), 12);
    });

    it('reports equidistance at the circle center', () => {
        const a = arc([2, -1], 3, 0.3, 2.0);
        const result = query.compute(v(2, -1), a);
        expect(result.equidistant).toBe(true);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(a.end[0].values[0], 12);
        expect(result.closest[1].values[1]).toBeCloseTo(a.end[0].values[1], 12);
    });

    it('reports zero distance for a point on the arc', () => {
        const a = arc([0, 0], 2, 0, Math.PI);
        const result = query.compute(v(0, 2), a);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('ties on equal endpoint distances choose the first endpoint', () => {
        // A symmetric arc, with a point on the axis of symmetry opposite the
        // arc so that both endpoints are equidistant.
        const a = arc([0, 0], 1, Math.PI / 4, 3 * Math.PI / 4);
        const result = query.compute(v(0, -4), a);
        expect(result.closest[1].values[0]).toBeCloseTo(Math.SQRT1_2, 12);
    });

    it('handles a zero-radius arc', () => {
        const c = v(1, 1);
        const a = Arc2.fromCenterRadiusEnds(c, 0, c, c);
        const result = query.compute(v(4, 5), a);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values).toEqual([1, 1]);
    });

    it('agrees with a dense sampling of the arc', () => {
        const rnd = makeRandom(20260901);
        const a0 = -0.4, a1 = 2.2;
        const center = [0.5, -0.75];
        const radius = 1.75;
        const a = arc(center, radius, a0, a1);

        for (let trial = 0; trial < 60; ++trial) {
            const p = v(8 * rnd() - 4, 8 * rnd() - 4);
            const result = query.compute(p, a);

            // The reported closest point is on the circle of the arc.
            const rad = Math.sqrt(dot(sub(result.closest[1], v(...center)),
                sub(result.closest[1], v(...center))));
            expect(rad).toBeCloseTo(radius, 8);

            // The reported closest point realizes the reported distance.
            const d = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(d, d))).toBeCloseTo(result.distance, 10);

            // No sampled arc point is closer.
            const n = 4000;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                const t = a0 + (i / n) * (a1 - a0);
                const q = v(center[0] + radius * Math.cos(t),
                    center[1] + radius * Math.sin(t));
                const e = sub(p, q);
                best = Math.min(best, dot(e, e));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});
