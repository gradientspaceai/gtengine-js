import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, sub } from '../src/Vector';
import { inContainerAlignedBox } from '../src/ContAlignedBox';
import {
    IntrSegment2AlignedBox2TI,
    IntrSegment2AlignedBox2FI
} from '../src/IntrSegment2AlignedBox2';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrSegment2AlignedBox2', () => {
    const ti = new IntrSegment2AlignedBox2TI();
    const fi = new IntrSegment2AlignedBox2FI();
    const unit = box([-1, -1], [1, 1]);

    it('clips a segment crossing the box', () => {
        const s = segment([-3, 0], [3, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The endpoint-form parameters for the crossing at x = -1 and x = 1.
        expect(result.parameter[0]).toBeCloseTo(1 / 3, 12);
        expect(result.parameter[1]).toBeCloseTo(2 / 3, 12);
        // The centered form of the segment has extent 3 and direction (1,0).
        expect(result.cdeParameter[0]).toBeCloseTo(-1, 12);
        expect(result.cdeParameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports the whole segment when it is inside the box', () => {
        const s = segment([-0.5, -0.25], [0.5, 0.25]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0.5, 12);
    });

    it('clips a segment with one endpoint inside the box', () => {
        const s = segment([0, 0], [4, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports 2 intersections when the segment touches only a corner', () => {
        // The segment touches the box only at the corner (1,1).
        const s = segment([0, 2], [2, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        // Upstream promotes a single-point overlap to 2 intersections so the
        // caller computes both (identical) points.
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('misses when the segment stops short of the box', () => {
        const s = segment([-5, 0], [-2, 0]);
        expect(ti.test(s, unit).intersect).toBe(false);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('misses when the supporting line misses the box', () => {
        const s = segment([-5, 3], [5, 3]);
        expect(ti.test(s, unit).intersect).toBe(false);
        expect(fi.find(s, unit).intersect).toBe(false);
    });

    it('handles a degenerate segment inside the box', () => {
        const s = segment([0.25, -0.5], [0.25, -0.5]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.cdeParameter).toEqual([0, 0]);
        expect(result.point[0].values[0]).toBeCloseTo(0.25, 12);
        expect(result.point[1].values[1]).toBeCloseTo(-0.5, 12);
    });

    it('handles a degenerate segment outside the box', () => {
        const s = segment([5, 5], [5, 5]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('handles a translated box', () => {
        const shifted = box([9, 9], [11, 11]);
        const s = segment([7, 10], [13, 10]);
        expect(ti.test(s, shifted).intersect).toBe(true);
        const result = fi.find(s, shifted);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(9, 12);
        expect(result.point[1].values[0]).toBeCloseTo(11, 12);
        expect(result.point[0].values[1]).toBeCloseTo(10, 12);
    });

    it('agrees with brute-force sampling on random configurations', () => {
        const rnd = makeRandom(31415);
        let tiFiMismatch = 0;
        let sampleMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;
        const samples = 3000;

        for (let trial = 0; trial < 250; ++trial) {
            const cx = 2 * rnd() - 1, cy = 2 * rnd() - 1;
            const ex = 0.2 + rnd(), ey = 0.2 + rnd();
            const b = box([cx - ex, cy - ey], [cx + ex, cy + ey]);
            const p0 = vec(6 * rnd() - 3, 6 * rnd() - 3);
            const p1 = add(p0, add(mul(2, sub(vec(cx, cy), p0)),
                vec(3 * rnd() - 1.5, 3 * rnd() - 1.5)));
            const d = sub(p1, p0);
            if (dot(d, d) < 1e-6) {
                continue;
            }
            const s = Segment.fromEndpoints(p0, p1);

            const tiResult = ti.test(s, b);
            const fiResult = fi.find(s, b);
            if (tiResult.intersect !== fiResult.intersect) {
                ++tiFiMismatch;
            }

            let sampled = false;
            for (let k = 0; k <= samples; ++k) {
                if (inContainerAlignedBox(add(p0, mul(k / samples, d)), b)) {
                    sampled = true;
                    break;
                }
            }
            if (sampled && !fiResult.intersect) {
                ++sampleMismatch;
            }

            if (fiResult.intersect) {
                ++hits;
                for (let i = 0; i < fiResult.numIntersections; ++i) {
                    const t = fiResult.parameter[i];
                    if (t < -1e-9 || t > 1 + 1e-9) {
                        ++pointMismatch;
                    }
                    // The reported point must match the endpoint-form
                    // evaluation of the reported parameter.
                    const expected = add(p0, mul(t, d));
                    const diff = sub(fiResult.point[i], expected);
                    if (Math.sqrt(dot(diff, diff)) > 1e-9) {
                        ++pointMismatch;
                    }
                    // ... and it must be in the box (up to rounding).
                    const grown = AlignedBox.fromMinMax(
                        sub(b.min, vec(1e-9, 1e-9)),
                        add(b.max, vec(1e-9, 1e-9)));
                    if (!inContainerAlignedBox(expected, grown)) {
                        ++pointMismatch;
                    }
                }
            }
        }

        expect(hits).toBeGreaterThan(20);
        expect([tiFiMismatch, sampleMismatch, pointMismatch]).toEqual([0, 0, 0]);
    });
});
