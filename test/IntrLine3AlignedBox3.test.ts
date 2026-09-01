import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { Line } from '../src/Line';
import { Vector, normalize } from '../src/Vector';
import {
    IntrLine3AlignedBox3TI,
    IntrLine3AlignedBox3FI
} from '../src/IntrLine3AlignedBox3';

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

function line(p: number[], d: number[]): Line {
    return Line.fromOriginDirection(Vector.fromArray(p), Vector.fromArray(d));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrLine3AlignedBox3', () => {
    const ti = new IntrLine3AlignedBox3TI();
    const fi = new IntrLine3AlignedBox3FI();
    const unit = box([-1, -1, -1], [1, 1, 1]);

    it('finds the entry and exit of an axis-parallel line', () => {
        const l = line([0, 0, 0], [0, 0, 1]);
        expect(ti.test(l, unit).intersect).toBe(true);
        const result = fi.find(l, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values).toEqual([0, 0, -1]);
        expect(result.point[1].values).toEqual([0, 0, 1]);
    });

    it('finds the crossing of a main diagonal', () => {
        // The diagonal direction (1,1,1) through the origin crosses the unit
        // cube at t = -1/sqrt(3) and +1/sqrt(3) once normalized.
        const d = Vector.fromArray([1, 1, 1]);
        normalize(d);
        const l = Line.fromOriginDirection(Vector.zero(3), d);
        const result = fi.find(l, unit);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[1]).toBeCloseTo(Math.sqrt(3), 12);
        expect(result.parameter[0]).toBeCloseTo(-Math.sqrt(3), 12);
        for (let k = 0; k < 3; ++k) {
            expect(result.point[1].values[k]).toBeCloseTo(1, 12);
            expect(result.point[0].values[k]).toBeCloseTo(-1, 12);
        }
    });

    it('handles a non-centered box and an off-origin line', () => {
        // The line from (0,2,2) along +x through the box [1,3]x[1,3]x[1,3].
        const result = fi.find(line([0, 2, 2], [1, 0, 0]), box([1, 1, 1], [3, 3, 3]));
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(1, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
    });

    it('reports a single point for a line touching an edge of the box', () => {
        // Direction (1,-1,0) through (0,2,0) touches the cube only at the
        // edge point (1,1,z=0).
        const l = line([0, 2, 0], [1, -1, 0]);
        expect(ti.test(l, unit).intersect).toBe(true);
        const result = fi.find(l, unit);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBe(result.parameter[1]);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection for a line that misses the box', () => {
        const l = line([0, 3, 0], [1, 0, 0]);
        expect(ti.test(l, unit).intersect).toBe(false);
        const result = fi.find(l, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('handles degenerate boxes', () => {
        // A flat box (a square in the plane z = 0).
        const flat = box([-1, -1, 0], [1, 1, 0]);
        expect(ti.test(line([0, 0, -5], [0, 0, 1]), flat).intersect).toBe(true);
        const result = fi.find(line([0, 0, -5], [0, 0, 1]), flat);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(5, 12);
        // A point box.
        const point = box([1, 2, 3], [1, 2, 3]);
        const d = Vector.fromArray([1, 2, 3]);
        normalize(d);
        expect(ti.test(Line.fromOriginDirection(Vector.zero(3), d), point).intersect)
            .toBe(true);
        expect(ti.test(line([0, 0, 0], [1, 0, 0]), point).intersect).toBe(false);
    });

    it('agrees with a sampling oracle and keeps TI and FI consistent', () => {
        const rand = makeRandom(20260101);
        let numHit = 0, numMiss = 0;
        for (let trial = 0; trial < 600; ++trial) {
            const lo: number[] = [], hi: number[] = [];
            for (let k = 0; k < 3; ++k) {
                const a = 3 * rand() - 1.5;
                lo.push(a);
                hi.push(a + 0.2 + rand());
            }
            const b = box(lo, hi);
            const d = Vector.fromArray([2 * rand() - 1, 2 * rand() - 1,
                2 * rand() - 1]);
            if (normalize(d) < 1e-6) {
                continue;
            }
            const l = Line.fromOriginDirection(Vector.fromArray(
                [3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5]), d);

            const t = ti.test(l, b).intersect;
            const f = fi.find(l, b);
            expect(f.intersect).toBe(t);

            if (f.intersect) {
                ++numHit;
                expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
                for (let k = 0; k <= 8; ++k) {
                    const s = f.parameter[0] +
                        (k / 8) * (f.parameter[1] - f.parameter[0]);
                    for (let dim = 0; dim < 3; ++dim) {
                        const v = l.origin.values[dim] + s * d.values[dim];
                        expect(v).toBeGreaterThan(b.min.values[dim] - 1e-9);
                        expect(v).toBeLessThan(b.max.values[dim] + 1e-9);
                    }
                }
            } else {
                ++numMiss;
                let anyInside = false;
                for (let k = -400; k <= 400; ++k) {
                    const s = k * 0.02;
                    let inside = true;
                    for (let dim = 0; dim < 3; ++dim) {
                        const v = l.origin.values[dim] + s * d.values[dim];
                        if (v < b.min.values[dim] || v > b.max.values[dim]) {
                            inside = false;
                        }
                    }
                    anyInside = anyInside || inside;
                }
                expect(anyInside).toBe(false);
            }
        }
        expect(numHit).toBeGreaterThan(30);
        expect(numMiss).toBeGreaterThan(50);
    });
});
