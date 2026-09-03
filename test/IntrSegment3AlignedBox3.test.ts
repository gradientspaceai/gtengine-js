import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { inContainerAlignedBox } from '../src/ContAlignedBox.js';
import {
    IntrSegment3AlignedBox3TI,
    IntrSegment3AlignedBox3FI,
    defaultIntrSegment3AlignedBox3FIResult,
    intrSegment3AlignedBox3FIDoQuery
} from '../src/IntrSegment3AlignedBox3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
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

describe('IntrSegment3AlignedBox3', () => {
    const ti = new IntrSegment3AlignedBox3TI();
    const fi = new IntrSegment3AlignedBox3FI();
    const unit = box([-1, -1, -1], [1, 1, 1]);

    it('clips a segment crossing the box', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The parameters are relative to the centered form, extent 3.
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports the whole segment when it lies inside the box', () => {
        const s = segment([-0.5, -0.25, 0.1], [0.5, 0.25, -0.1]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(-0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0.5, 12);
    });

    it('clips a diagonal segment through the box', () => {
        const s = segment([-3, -3, -3], [3, 3, 3]);
        const result = fi.find(s, unit);
        expect(result.numIntersections).toBe(2);
        for (let j = 0; j < 3; ++j) {
            expect(result.point[0].values[j]).toBeCloseTo(-1, 12);
            expect(result.point[1].values[j]).toBeCloseTo(1, 12);
        }
    });

    it('reports a single point when the segment touches a corner', () => {
        // The segment lies in the plane z = 1 and touches the box only at
        // the corner (1,1,1).
        const s = segment([0, 2, 1], [2, 0, 1]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[2]).toBeCloseTo(1, 12);
    });

    it('misses when the segment stops short of the box', () => {
        const s = segment([-5, 0, 0], [-2, 0, 0]);
        expect(ti.test(s, unit).intersect).toBe(false);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('misses when the supporting line misses the box', () => {
        const s = segment([-5, 3, 0], [5, 3, 0]);
        expect(ti.test(s, unit).intersect).toBe(false);
        expect(fi.find(s, unit).intersect).toBe(false);
    });

    it('handles a degenerate segment', () => {
        // A zero-length segment has zero direction after normalize; the
        // separating-axis test then compares |origin[i]| with boxExtent[i].
        const inside = segment([0.25, -0.5, 0], [0.25, -0.5, 0]);
        expect(ti.test(inside, unit).intersect).toBe(true);
        const outside = segment([5, 5, 5], [5, 5, 5]);
        expect(ti.test(outside, unit).intersect).toBe(false);
        expect(fi.find(outside, unit).intersect).toBe(false);
    });

    it('handles a translated box', () => {
        const shifted = box([9, 9, 9], [11, 11, 11]);
        const s = segment([7, 10, 10], [13, 10, 10]);
        expect(ti.test(s, shifted).intersect).toBe(true);
        const result = fi.find(s, shifted);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(9, 12);
        expect(result.point[1].values[0]).toBeCloseTo(11, 12);
        expect(result.point[0].values[1]).toBeCloseTo(10, 12);
    });

    it('the exported FI DoQuery matches the class query', () => {
        const s = segment([-4, 0.3, -0.2], [4, -0.1, 0.5]);
        const { center: boxCenter, extent: boxExtent } = unit.getCenteredForm();
        const shifted = Segment.fromEndpoints(sub(s.p[0], boxCenter),
            sub(s.p[1], boxCenter));
        const cf = shifted.getCenteredForm();
        const direct = defaultIntrSegment3AlignedBox3FIResult();
        intrSegment3AlignedBox3FIDoQuery(cf.center, cf.direction, cf.extent,
            boxExtent, direct);
        const viaClass = fi.find(s, unit);
        expect(direct.intersect).toBe(viaClass.intersect);
        expect(direct.numIntersections).toBe(viaClass.numIntersections);
        expect(direct.parameter[0]).toBeCloseTo(viaClass.parameter[0], 12);
        expect(direct.parameter[1]).toBeCloseTo(viaClass.parameter[1], 12);
    });

    it('agrees with brute-force sampling on random configurations', () => {
        const rnd = makeRandom(271828);
        let tiFiMismatch = 0;
        let sampleMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;
        const samples = 2000;

        for (let trial = 0; trial < 200; ++trial) {
            const c = vec(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            const e = vec(0.2 + rnd(), 0.2 + rnd(), 0.2 + rnd());
            const b = AlignedBox.fromMinMax(sub(c, e), add(c, e));
            const p0 = vec(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const p1 = add(p0, add(mul(2, sub(c, p0)),
                vec(3 * rnd() - 1.5, 3 * rnd() - 1.5, 3 * rnd() - 1.5)));
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
                const cf = s.getCenteredForm();
                const eps = vec(1e-9, 1e-9, 1e-9);
                const grown = AlignedBox.fromMinMax(sub(b.min, eps),
                    add(b.max, eps));
                for (let i = 0; i < fiResult.numIntersections; ++i) {
                    const t = fiResult.parameter[i];
                    if (Math.abs(t) > cf.extent + 1e-9) {
                        ++pointMismatch;
                    }
                    const expected = add(cf.center, mul(t, cf.direction));
                    const diff = sub(fiResult.point[i], expected);
                    if (Math.sqrt(dot(diff, diff)) > 1e-9) {
                        ++pointMismatch;
                    }
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
