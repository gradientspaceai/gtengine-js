import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { Line } from '../src/Line';
import { OrientedBox } from '../src/OrientedBox';
import { Vector, add, mul, sub, dot, normalize } from '../src/Vector';
import {
    IntrLine2AlignedBox2TI,
    IntrLine2AlignedBox2FI
} from '../src/IntrLine2AlignedBox2';
import {
    IntrLine2OrientedBox2TI,
    IntrLine2OrientedBox2FI
} from '../src/IntrLine2OrientedBox2';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(px: number, py: number, dx: number, dy: number): Line {
    return Line.fromOriginDirection(v2(px, py), v2(dx, dy));
}

// An oriented box with the given center, rotation angle and extents.
function obox(cx: number, cy: number, angle: number, e0: number,
    e1: number): OrientedBox {
    const c = Math.cos(angle), s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(v2(cx, cy),
        [v2(c, s), v2(-s, c)], v2(e0, e1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Brute-force containment in the solid oriented box.
function inBox(p: Vector, box: OrientedBox): boolean {
    const diff = sub(p, box.center);
    for (let i = 0; i < 2; ++i) {
        if (Math.abs(dot(diff, box.axis[i])) > box.extent.values[i] + 1e-12) {
            return false;
        }
    }
    return true;
}

describe('IntrLine2OrientedBox2', () => {
    const ti = new IntrLine2OrientedBox2TI();
    const fi = new IntrLine2OrientedBox2FI();
    const abTI = new IntrLine2AlignedBox2TI();
    const abFI = new IntrLine2AlignedBox2FI();

    it('finds the entry and exit parameters for an axis-aligned box', () => {
        // A zero-angle oriented box is the aligned box [-1,1]^2.
        const result = fi.find(line(0, 0, 1, 0), obox(0, 0, 0, 1, 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a rotated box', () => {
        // A unit square rotated 45 degrees has corners at distance sqrt(2)
        // along the coordinate axes.
        const box = obox(0, 0, Math.PI / 4, 1, 1);
        const result = fi.find(line(0, 0, 1, 0), box);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-Math.SQRT2, 12);
        expect(result.parameter[1]).toBeCloseTo(Math.SQRT2, 12);

        // A horizontal line above the top corner misses.
        expect(ti.test(line(0, 1.5, 1, 0), box).intersect).toBe(false);
        expect(ti.test(line(0, 1.4, 1, 0), box).intersect).toBe(true);
    });

    it('reports a corner contact for a tangent line', () => {
        // A rotated box touched by a horizontal line at its topmost corner.
        // The topmost corner of the box with axes (0.6,0.8) and (-0.8,0.6)
        // and unit extents is (-0.2, 1.4).
        const box = OrientedBox.fromCenterAxisExtent(v2(0, 0),
            [v2(0.6, 0.8), v2(-0.8, 0.6)], v2(1, 1));
        const result = fi.find(line(0, 1.4, 1, 0), box);
        expect(result.intersect).toBe(true);
        for (let i = 0; i < result.numIntersections; ++i) {
            expect(result.point[i].values[0]).toBeCloseTo(-0.2, 9);
            expect(result.point[i].values[1]).toBeCloseTo(1.4, 9);
        }

        // Just above the corner there is no intersection.
        expect(ti.test(line(0, 1.4 + 1e-6, 1, 0), box).intersect).toBe(false);
    });

    it('matches the aligned-box query when the axes are the standard basis',
        () => {
            const rnd = makeRandom(5150);
            const ab = AlignedBox.fromMinMax(v2(-1, -2), v2(3, 1));
            const { center, extent } = ab.getCenteredForm();
            const ob = OrientedBox.fromCenterAxisExtent(center,
                [v2(1, 0), v2(0, 1)], extent);
            for (let k = 0; k < 300; ++k) {
                const l = line(rnd() * 8 - 4, rnd() * 8 - 4,
                    rnd() * 2 - 1, rnd() * 2 - 1);
                if (dot(l.direction, l.direction) < 1e-8) {
                    continue;
                }
                expect(ti.test(l, ob).intersect)
                    .toBe(abTI.test(l, ab).intersect);
                const r0 = fi.find(l, ob);
                const r1 = abFI.find(l, ab);
                expect(r0.numIntersections).toBe(r1.numIntersections);
                for (let i = 0; i < r0.numIntersections; ++i) {
                    expect(r0.parameter[i]).toBeCloseTo(r1.parameter[i], 9);
                }
            }
        });

    it('agrees with the TI query and a brute-force sampling', () => {
        const rnd = makeRandom(606);
        let numHit = 0, numMiss = 0;
        for (let k = 0; k < 300; ++k) {
            const box = obox(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * Math.PI,
                0.2 + rnd(), 0.2 + rnd());
            const d = v2(rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(d, d) < 1e-8) {
                continue;
            }
            normalize(d);
            const l = Line.fromOriginDirection(
                v2(rnd() * 6 - 3, rnd() * 6 - 3), d);

            const tiResult = ti.test(l, box).intersect;
            const fiResult = fi.find(l, box);
            expect(tiResult).toBe(fiResult.intersect);

            // Every reported intersection point must be in the box.
            for (let i = 0; i < fiResult.numIntersections; ++i) {
                expect(inBox(fiResult.point[i], box)).toBe(true);
            }

            // Dense sampling of the line: a sample inside the box forces the
            // query to report an intersection.
            let sampleHit = false;
            for (let s = -600; s <= 600 && !sampleHit; ++s) {
                if (inBox(add(l.origin, mul(s * 0.01, l.direction)), box)) {
                    sampleHit = true;
                }
            }
            if (sampleHit) {
                expect(tiResult).toBe(true);
                ++numHit;
            } else if (!tiResult) {
                ++numMiss;
            }
        }
        expect(numHit).toBeGreaterThan(0);
        expect(numMiss).toBeGreaterThan(0);
    });

    it('handles a degenerate box with zero extents', () => {
        // The box degenerates to the point (1,1).
        const box = obox(1, 1, 0.3, 0, 0);
        expect(ti.test(line(0, 0, 1, 1), box).intersect).toBe(true);
        expect(ti.test(line(0, 0.5, 1, 1), box).intersect).toBe(false);
    });
});
