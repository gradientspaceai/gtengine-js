import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrSegment2AlignedBox2TI,
    IntrSegment2AlignedBox2FI
} from '../src/IntrSegment2AlignedBox2.js';
import {
    IntrSegment2OrientedBox2TI,
    IntrSegment2OrientedBox2FI,
    defaultIntrSegment2OrientedBox2FIResult
} from '../src/IntrSegment2OrientedBox2.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function rotatedBox(center: Vector, angle: number, extent: Vector):
    OrientedBox {
    const c = Math.cos(angle), s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(center, [vec(c, s), vec(-s, c)],
        extent);
}

function boxSignedDepth(box: OrientedBox, P: Vector): number {
    const diff = sub(P, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 2; ++i) {
        const value = Math.abs(dot(diff, box.axis[i])) - box.extent.values[i];
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

describe('IntrSegment2OrientedBox2', () => {
    const ti = new IntrSegment2OrientedBox2TI();
    const fi = new IntrSegment2OrientedBox2FI();
    const unitAxes = [vec(1, 0), vec(0, 1)];

    it('matches the aligned-box query when the axes are standard', () => {
        const center = vec(1, -2);
        const extent = vec(2, 0.5);
        const obox = OrientedBox.fromCenterAxisExtent(center, unitAxes,
            extent);
        const abox = AlignedBox.fromMinMax(sub(center, extent),
            add(center, extent));
        const aTI = new IntrSegment2AlignedBox2TI();
        const aFI = new IntrSegment2AlignedBox2FI();

        const segments = [
            segment([-10, -2], [10, -2]),
            segment([1, -10], [1, 10]),
            segment([-10, 10], [10, 10]),
            segment([1, -2], [4, -2]),
            segment([-10, -2], [-5, -2])
        ];
        for (const s of segments) {
            expect(ti.test(s, obox).intersect)
                .toBe(aTI.test(s, abox).intersect);
            const o = fi.find(s, obox);
            const a = aFI.find(s, abox);
            expect(o.intersect).toBe(a.intersect);
            expect(o.numIntersections).toBe(a.numIntersections);
            for (let i = 0; i < o.numIntersections; ++i) {
                expect(o.parameter[i]).toBeCloseTo(a.parameter[i], 12);
                expect(o.cdeParameter[i]).toBeCloseTo(a.cdeParameter[i], 12);
                expect(o.point[i].values[0])
                    .toBeCloseTo(a.point[i].values[0], 12);
                expect(o.point[i].values[1])
                    .toBeCloseTo(a.point[i].values[1], 12);
            }
        }
    });

    it('reports correct world-space points for a rotated box', () => {
        // Regression test for the upstream point-reconstruction bug: with a
        // 45-degree rotation the box is a diamond reaching sqrt(2) along x.
        const box = rotatedBox(vec(0, 0), Math.PI / 4, vec(1, 1));
        const result = fi.find(segment([-5, 0], [5, 0]), box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(-Math.SQRT2, 9);
        expect(result.point[0].values[1]).toBeCloseTo(0, 9);
        expect(result.point[1].values[0]).toBeCloseTo(Math.SQRT2, 9);
        expect(result.point[1].values[1]).toBeCloseTo(0, 9);

        // The endpoint-form parameters recover the same points.
        for (let i = 0; i < 2; ++i) {
            const t = result.parameter[i];
            const Q = add(mul(1 - t, vec(-5, 0)), mul(t, vec(5, 0)));
            expect(Q.values[0]).toBeCloseTo(result.point[i].values[0], 9);
            expect(Q.values[1]).toBeCloseTo(result.point[i].values[1], 9);
        }
    });

    it('reports parameters in [0,1] and the centered form in cdeParameter', () => {
        const box = rotatedBox(vec(1, 1), 0.4, vec(1, 0.5));
        const s = segment([-4, 1], [6, 1]);
        const result = fi.find(s, box);
        expect(result.intersect).toBe(true);
        const { center, direction, extent } = s.getCenteredForm();
        for (let i = 0; i < result.numIntersections; ++i) {
            expect(result.parameter[i]).toBeGreaterThanOrEqual(0);
            expect(result.parameter[i]).toBeLessThanOrEqual(1);
            expect(Math.abs(result.cdeParameter[i]))
                .toBeLessThanOrEqual(extent + 1e-9);
            const P = add(center, mul(result.cdeParameter[i], direction));
            expect(P.values[0]).toBeCloseTo(result.point[i].values[0], 9);
            expect(P.values[1]).toBeCloseTo(result.point[i].values[1], 9);
        }
    });

    it('rejects a segment that stops short of the box', () => {
        const box = rotatedBox(vec(0, 0), 0.3, vec(1, 1));
        const s = segment([-10, 0], [-5, 0]);
        expect(fi.find(s, box).intersect).toBe(false);
        expect(ti.test(s, box).intersect).toBe(false);
    });

    it('handles a degenerate segment (upstream reports 2 intersections)', () => {
        const box = rotatedBox(vec(0, 0), 0.3, vec(1, 1));
        const inside = segment([0.1, -0.2], [0.1, -0.2]);
        const result = fi.find(inside, box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.cdeParameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0.1, -0.2]);
        expect(result.point[1].values).toEqual([0.1, -0.2]);

        const outside = segment([9, 9], [9, 9]);
        expect(fi.find(outside, box).intersect).toBe(false);
    });

    it('has the documented default result', () => {
        const result = defaultIntrSegment2OrientedBox2FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.cdeParameter).toEqual([0, 0]);
    });

    it('rejects non-2D boxes', () => {
        const box3 = new OrientedBox(3);
        const s = segment([0, 0], [1, 1]);
        expect(() => ti.test(s, box3)).toThrow();
        expect(() => fi.find(s, box3)).toThrow();
    });

    it('agrees with brute-force sampling on random segments', () => {
        let seed = 20240901;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const box = rotatedBox(vec(0.5, -0.25), 0.7, vec(1.5, 0.75));
        for (let trial = 0; trial < 400; ++trial) {
            const s = segment([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 8 - 4, rand() * 8 - 4]);
            const result = fi.find(s, box);
            expect(ti.test(s, box).intersect).toBe(result.intersect);

            let sampledHit = false;
            for (let k = 0; k <= 400; ++k) {
                const t = k / 400;
                const P = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
                if (boxSignedDepth(box, P) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            for (let i = 0; i < result.numIntersections; ++i) {
                const t = result.parameter[i];
                expect(t).toBeGreaterThanOrEqual(-1e-9);
                expect(t).toBeLessThanOrEqual(1 + 1e-9);
                const P = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
                expect(P.values[0]).toBeCloseTo(result.point[i].values[0], 8);
                expect(P.values[1]).toBeCloseTo(result.point[i].values[1], 8);
                expect(boxSignedDepth(box, result.point[i]))
                    .toBeLessThan(1e-8);
            }
        }
    });
});
