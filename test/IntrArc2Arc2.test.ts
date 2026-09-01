import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2';
import { Vector, sub, length } from '../src/Vector';
import {
    IntrArc2Arc2FI,
    IntrArc2Arc2Configuration as Cfg
} from '../src/IntrArc2Arc2';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// A point on the unit circle at the given angle in degrees. The four axis
// directions are exact so that endpoint equality tests are exact.
function pt(deg: number): Vector {
    switch (((deg % 360) + 360) % 360) {
        case 0: return v2(1, 0);
        case 90: return v2(0, 1);
        case 180: return v2(-1, 0);
        case 270: return v2(0, -1);
        default: {
            const a = (deg * Math.PI) / 180;
            return v2(Math.cos(a), Math.sin(a));
        }
    }
}

// The arc of the unit circle from deg0 counterclockwise to deg1.
function arc(deg0: number, deg1: number): Arc2 {
    return Arc2.fromCenterRadiusEnds(v2(0, 0), 1, pt(deg0), pt(deg1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrArc2Arc2', () => {
    const fi = new IntrArc2Arc2FI();

    it('reports no intersection for arcs on disjoint circles', () => {
        const a0 = arc(0, 90);
        const a1 = Arc2.fromCenterRadiusEnds(v2(10, 0), 1, pt(0), pt(90));
        const result = fi.find(a0, a1);
        expect(result.intersect).toBe(false);
        expect(result.configuration).toBe(Cfg.NO_INTERSECTION);
    });

    it('finds two noncocircular intersection points', () => {
        // Two unit circles centered at (0,0) and (1,1) meet at (1,0) and
        // (0,1). Use arcs that contain both points.
        const a0 = arc(-45, 135);
        const a1 = Arc2.fromCenterRadiusEnds(v2(1, 1), 1,
            v2(1, 1 - 1), v2(1 - 1, 1));  // from (1,0) to (0,1)
        const result = fi.find(a0, a1);
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.NONCOCIRCULAR_TWO_POINTS);
        for (let i = 0; i < 2; ++i) {
            expect(length(result.point[i])).toBeCloseTo(1, 12);
            expect(length(sub(result.point[i], v2(1, 1))))
                .toBeCloseTo(1, 12);
        }
    });

    it('finds one noncocircular intersection point', () => {
        // The arc of the second circle contains only (0,1).
        const a0 = arc(45, 135);
        const a1 = Arc2.fromCenterRadiusEnds(v2(1, 1), 1, v2(0, 1), v2(1, 0));
        const result = fi.find(a0, a1);
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.NONCOCIRCULAR_ONE_POINT);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('reports a single arc when the arcs are identical', () => {
        const a = arc(0, 90);
        const result = fi.find(a, a.clone());
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.COCIRCULAR_ONE_ARC);
        expect(result.arc[0].equals(a)).toBe(true);
    });

    it('reports the containing arc when arc0 is inside arc1', () => {
        const a0 = arc(30, 60);
        const a1 = arc(0, 90);
        const result = fi.find(a0, a1);
        expect(result.configuration).toBe(Cfg.COCIRCULAR_ONE_ARC);
        expect(result.arc[0].equals(a0)).toBe(true);
    });

    it('reports the contained arc when arc1 is inside arc0', () => {
        const a0 = arc(0, 90);
        const a1 = arc(30, 60);
        const result = fi.find(a0, a1);
        expect(result.configuration).toBe(Cfg.COCIRCULAR_ONE_ARC);
        expect(result.arc[0].equals(a1)).toBe(true);
    });

    it('reports the overlap of two partially overlapping cocircular arcs',
        () => {
            // <B0,A0,B1,A1> with arc0 = [0,90] and arc1 = [270,45].
            const a0 = arc(0, 90);
            const a1 = arc(270, 45);
            const result = fi.find(a0, a1);
            expect(result.intersect).toBe(true);
            expect(result.configuration).toBe(Cfg.COCIRCULAR_ONE_ARC);
            expect(result.arc[0].end[0].equals(pt(0))).toBe(true);
            expect(result.arc[0].end[1].equals(pt(45))).toBe(true);

            // <A0,B0,A1,B1> with arc0 = [0,90] and arc1 = [45,135].
            const b = fi.find(arc(0, 90), arc(45, 135));
            expect(b.configuration).toBe(Cfg.COCIRCULAR_ONE_ARC);
            expect(b.arc[0].end[0].equals(pt(45))).toBe(true);
            expect(b.arc[0].end[1].equals(pt(90))).toBe(true);
        });

    it('reports two disjoint overlap arcs', () => {
        // arc0 = [0,270] and arc1 = [180,90] overlap in [0,90] and
        // [180,270].
        const result = fi.find(arc(0, 270), arc(180, 90));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.COCIRCULAR_TWO_ARCS);
        expect(result.arc[0].end[0].equals(pt(0))).toBe(true);
        expect(result.arc[0].end[1].equals(pt(90))).toBe(true);
        expect(result.arc[1].end[0].equals(pt(180))).toBe(true);
        expect(result.arc[1].end[1].equals(pt(270))).toBe(true);
    });

    it('reports two points when the arcs share both endpoints', () => {
        const result = fi.find(arc(0, 180), arc(180, 0));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.COCIRCULAR_TWO_POINTS);
        expect(result.point[0].equals(pt(0))).toBe(true);
        expect(result.point[1].equals(pt(180))).toBe(true);
    });

    it('reports one point when the arcs share a single endpoint', () => {
        // arc0 = [0,90] and arc1 = [90,180] touch at 90 degrees.
        const result = fi.find(arc(0, 90), arc(90, 180));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.COCIRCULAR_ONE_POINT);
        expect(result.point[0].equals(pt(90))).toBe(true);

        // The reversed order touches at 0 degrees.
        const b = fi.find(arc(90, 180), arc(0, 90));
        expect(b.configuration).toBe(Cfg.COCIRCULAR_ONE_POINT);
        expect(b.point[0].equals(pt(90))).toBe(true);
    });

    it('reports one point and one arc', () => {
        // arc0 = [0,270] and arc1 = [270,90] overlap in [0,90] and touch at
        // 270 degrees.
        const result = fi.find(arc(0, 270), arc(270, 90));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.COCIRCULAR_ONE_POINT_ONE_ARC);
        expect(result.point[0].equals(pt(270))).toBe(true);
        expect(result.arc[0].end[0].equals(pt(0))).toBe(true);
        expect(result.arc[0].end[1].equals(pt(90))).toBe(true);
    });

    it('reports no intersection for disjoint cocircular arcs', () => {
        const result = fi.find(arc(0, 45), arc(90, 135));
        expect(result.intersect).toBe(false);
        expect(result.configuration).toBe(Cfg.NO_INTERSECTION);
    });

    it('keeps intersect consistent with the configuration', () => {
        const rnd = makeRandom(1729);
        let numIntersect = 0;
        for (let k = 0; k < 500; ++k) {
            const a0 = arc(rnd() * 360, rnd() * 360);
            const cx = rnd() < 0.4 ? 0 : rnd() * 3 - 1.5;
            const cy = rnd() < 0.4 ? 0 : rnd() * 3 - 1.5;
            const r = (cx === 0 && cy === 0 && rnd() < 0.5) ? 1
                : 0.3 + rnd() * 2;
            const d0 = rnd() * 360, d1 = rnd() * 360;
            const a1 = Arc2.fromCenterRadiusEnds(v2(cx, cy), r,
                v2(cx + r * Math.cos((d0 * Math.PI) / 180),
                    cy + r * Math.sin((d0 * Math.PI) / 180)),
                v2(cx + r * Math.cos((d1 * Math.PI) / 180),
                    cy + r * Math.sin((d1 * Math.PI) / 180)));

            const result = fi.find(a0, a1);
            expect(result.intersect)
                .toBe(result.configuration !== Cfg.NO_INTERSECTION);

            if (result.configuration === Cfg.NONCOCIRCULAR_ONE_POINT ||
                result.configuration === Cfg.NONCOCIRCULAR_TWO_POINTS) {
                const n = result.configuration ===
                    Cfg.NONCOCIRCULAR_TWO_POINTS ? 2 : 1;
                for (let i = 0; i < n; ++i) {
                    const p = result.point[i];
                    expect(length(sub(p, a0.center)))
                        .toBeCloseTo(a0.radius, 8);
                    expect(length(sub(p, a1.center)))
                        .toBeCloseTo(a1.radius, 8);
                    expect(a0.containsOnCircle(p)).toBe(true);
                    expect(a1.containsOnCircle(p)).toBe(true);
                }
                ++numIntersect;
            }
        }
        expect(numIntersect).toBeGreaterThan(0);
    });
});
