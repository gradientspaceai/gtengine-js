import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Line } from '../src/Line.js';
import { Vector, normalize } from '../src/Vector.js';
import {
    IntrLine2AlignedBox2TI,
    IntrLine2AlignedBox2FI
} from '../src/IntrLine2AlignedBox2.js';
import {
    intrLine2AlignedBox2TIDoQuery,
    intrLine2AlignedBox2FIDoQuery,
    defaultIntrLine2AlignedBox2TIResult,
    defaultIntrLine2AlignedBox2FIResult
} from '../src/IntrLine2AlignedBox2.js';
import { add, mul, sub } from '../src/Vector.js';
import { check, expectVectorClose, fc, positive, seededRandom, unitVector, wellScaled } from './helpers/arbitraries.js';

function box(minX: number, minY: number, maxX: number, maxY: number): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray([minX, minY]),
        Vector.fromArray([maxX, maxY]));
}

function line(px: number, py: number, dx: number, dy: number): Line {
    return Line.fromOriginDirection(Vector.fromArray([px, py]),
        Vector.fromArray([dx, dy]));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrLine2AlignedBox2', () => {
    const ti = new IntrLine2AlignedBox2TI();
    const fi = new IntrLine2AlignedBox2FI();
    const unit = box(-1, -1, 1, 1);

    it('finds the entry and exit parameters of an axis-parallel line', () => {
        // The x axis crosses [-1,1]^2 at t = -1 and t = 1.
        const result = fi.find(line(0, 0, 1, 0), unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values).toEqual([-1, 0]);
        expect(result.point[1].values).toEqual([1, 0]);
        expect(ti.test(line(0, 0, 1, 0), unit).intersect).toBe(true);
    });

    it('handles a line whose origin is outside the box', () => {
        // From (-5,0) along +x: enters at t = 4 and exits at t = 6.
        const result = fi.find(line(-5, 0, 1, 0), unit);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
    });

    it('handles a diagonal line through a non-centered box', () => {
        // The line through (0,0) with direction (1,1) crosses [1,3]x[1,3] at
        // t = 1 and t = 3.
        const result = fi.find(line(0, 0, 1, 1), box(1, 1, 3, 3));
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(1, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[1].values).toEqual([3, 3]);
    });

    it('reports a single point for a line touching a corner', () => {
        // Direction (1,-1) through (0,2) touches [-1,1]^2 only at (1,1).
        const l = line(0, 2, 1, -1);
        expect(ti.test(l, unit).intersect).toBe(true);
        const result = fi.find(l, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBe(result.parameter[1]);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('reports a grazing line along an edge as a segment', () => {
        // The line y = 1 lies on the top edge of the box.
        const result = fi.find(line(0, 1, 1, 0), unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection for a line that misses the box', () => {
        const l = line(0, 3, 1, 0);
        expect(ti.test(l, unit).intersect).toBe(false);
        const result = fi.find(l, unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('handles a degenerate (zero-extent) box', () => {
        // A "box" that is the single point (2,3).
        const point = box(2, 3, 2, 3);
        const through = line(0, 0, 2, 3);
        normalize(through.direction);
        expect(ti.test(through, point).intersect).toBe(true);
        expect(ti.test(line(0, 0, 1, 0), point).intersect).toBe(false);
    });

    it('agrees with a sampling oracle and keeps TI and FI consistent', () => {
        const rand = makeRandom(778899);
        let numHit = 0, numMiss = 0;
        for (let trial = 0; trial < 800; ++trial) {
            const x0 = 3 * rand() - 1.5, y0 = 3 * rand() - 1.5;
            const b = box(x0, y0, x0 + 0.2 + rand(), y0 + 0.2 + rand());
            const d = Vector.fromArray([2 * rand() - 1, 2 * rand() - 1]);
            if (normalize(d) < 1e-6) {
                continue;
            }
            const l = Line.fromOriginDirection(
                Vector.fromArray([3 * rand() - 1.5, 3 * rand() - 1.5]), d);

            const t = ti.test(l, b).intersect;
            const f = fi.find(l, b);
            expect(f.intersect).toBe(t);

            if (f.intersect) {
                ++numHit;
                expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
                // Sample points of the reported parameter interval; they must
                // lie in the box (up to round-off).
                for (let k = 0; k <= 8; ++k) {
                    const s = f.parameter[0] +
                        (k / 8) * (f.parameter[1] - f.parameter[0]);
                    for (let dim = 0; dim < 2; ++dim) {
                        const v = l.origin.values[dim] + s * d.values[dim];
                        expect(v).toBeGreaterThan(b.min.values[dim] - 1e-9);
                        expect(v).toBeLessThan(b.max.values[dim] + 1e-9);
                    }
                }
                // The endpoints are on the box boundary.
                for (let i = 0; i < f.numIntersections; ++i) {
                    let onBoundary = false;
                    for (let dim = 0; dim < 2; ++dim) {
                        const v = f.point[i].values[dim];
                        if (Math.abs(v - b.min.values[dim]) < 1e-9 ||
                            Math.abs(v - b.max.values[dim]) < 1e-9) {
                            onBoundary = true;
                        }
                    }
                    expect(onBoundary).toBe(true);
                }
            } else {
                ++numMiss;
                // A dense sampling of the line must find no point in the box.
                let anyInside = false;
                for (let k = -400; k <= 400; ++k) {
                    const s = k * 0.02;
                    let inside = true;
                    for (let dim = 0; dim < 2; ++dim) {
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
        expect(numHit).toBeGreaterThan(50);
        expect(numMiss).toBeGreaterThan(50);
    });
});

describe('intrLine2AlignedBox2 DoQuery helpers', () => {
    // The helpers take the line in the box-centered coordinate system, which
    // is what the classes pass to them. The box below is already centered at
    // the origin, so no translation is needed.
    const b = box(-2, -1, 2, 1);
    const extent = Vector.fromArray([2, 1]);

    it('the TI helper matches the class query', () => {
        const cases: Array<[number, number, number, number]> = [
            [0, 0, 1, 0],
            [0, 5, 1, 0],
            [-4, -3, 1, 1],
            [3, 0, 0, 1],
            [-2, 1, 1, 0]
        ];
        for (const [px, py, dx, dy] of cases) {
            const d = Vector.fromArray([dx, dy]);
            normalize(d);
            const result = defaultIntrLine2AlignedBox2TIResult();
            intrLine2AlignedBox2TIDoQuery(Vector.fromArray([px, py]), d,
                extent, result);
            const expected = new IntrLine2AlignedBox2TI().test(
                line(px, py, d.values[0], d.values[1]), b);
            expect(result.intersect).toBe(expected.intersect);
        }
    });

    it('the FI helper fills parameters but not points', () => {
        const result = defaultIntrLine2AlignedBox2FIResult();
        intrLine2AlignedBox2FIDoQuery(Vector.fromArray([0, 0]),
            Vector.fromArray([1, 0]), extent, result);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
        // DoQuery leaves 'point' at its default value.
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('the FI helper reports no intersection for a missing line', () => {
        const result = defaultIntrLine2AlignedBox2FIResult();
        intrLine2AlignedBox2FIDoQuery(Vector.fromArray([0, 5]),
            Vector.fromArray([1, 0]), extent, result);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });
});

describe('IntrLine2AlignedBox2 verification', () => {
    const ti = new IntrLine2AlignedBox2TI();
    const fi = new IntrLine2AlignedBox2FI();

    const boxArb = fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4),
        positive(4), positive(4))
        .map(([x, y, w, h]) => box(x, y, x + w, y + h));
    const lineArb = fc.tuple(wellScaled(-6, 6), wellScaled(-6, 6),
        unitVector(2))
        .map(([x, y, d]) => line(x, y, d.values[0], d.values[1]));

    function inBox(p: Vector, b: AlignedBox, tol: number): boolean {
        for (let i = 0; i < 2; ++i) {
            if (p.values[i] < b.min.values[i] - tol
                || p.values[i] > b.max.values[i] + tol) {
                return false;
            }
        }
        return true;
    }

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(lineArb, boxArb), ([l, b]) => {
            expect(fi.find(l, b).intersect).toBe(ti.test(l, b).intersect);
        });
    });

    it('FI parameters are ordered and their points lie on the line and in the box', () => {
        check(fc.tuple(lineArb, boxArb), ([l, b]) => {
            const r = fi.find(l, b);
            if (!r.intersect) {
                expect(r.numIntersections).toBe(0);
                return;
            }
            expect(r.parameter[0]).toBeLessThanOrEqual(r.parameter[1]);
            if (r.numIntersections === 1) {
                expect(r.parameter[1]).toBe(r.parameter[0]);
            }
            const scale = 1 + Math.abs(b.min.values[0]) + Math.abs(b.max.values[0])
                + Math.abs(b.min.values[1]) + Math.abs(b.max.values[1]);
            for (let k = 0; k < r.numIntersections; ++k) {
                const p = r.point[k];
                expect(Number.isFinite(p.values[0])).toBe(true);
                expect(Number.isFinite(p.values[1])).toBe(true);
                expectVectorClose(p,
                    add(l.origin, mul(r.parameter[k], l.direction)), 0, 0);
                // Each clip parameter is one division of exactly representable
                // differences, so the point sits on the face to a few ulps.
                expect(inBox(p, b, 1e-12 * scale)).toBe(true);
            }
        });
    });

    it('the reported [t0,t1] is exactly the set of line parameters inside the box', () => {
        const rnd = seededRandom(0x3ac71d);
        check(fc.tuple(lineArb, boxArb), ([l, b]) => {
            const r = fi.find(l, b);
            const scale = 1 + Math.abs(b.min.values[0]) + Math.abs(b.max.values[0])
                + Math.abs(b.min.values[1]) + Math.abs(b.max.values[1]);
            for (let k = 0; k < 400; ++k) {
                const t = 20 * rnd() - 10;
                const p = add(l.origin, mul(t, l.direction));
                if (inBox(p, b, -1e-9 * scale)) {
                    // Strictly inside: the query must have found the hit and
                    // the parameter must lie in the reported interval.
                    expect(r.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(r.parameter[0] - 1e-9 * scale);
                    expect(t).toBeLessThanOrEqual(r.parameter[1] + 1e-9 * scale);
                }
            }
        }, 60);
    }, 30000);

    it('an axis-parallel line hitting a face reports the whole face span', () => {
        const b = box(-1, -1, 1, 1);
        const r = fi.find(line(0, 5, 0, -1), b);
        expect(r.intersect).toBe(true);
        expect(r.numIntersections).toBe(2);
        expect(r.parameter[0]).toBeCloseTo(4, 12);
        expect(r.parameter[1]).toBeCloseTo(6, 12);
    });

    it('a line grazing a corner reports a single point', () => {
        const b = box(0, 0, 2, 2);
        // The line through (2,0) with direction (1,1) touches only the corner.
        const d = Vector.fromArray([1, 1]);
        normalize(d);
        const r = fi.find(Line.fromOriginDirection(Vector.fromArray([2, 0]), d), b);
        expect(r.intersect).toBe(true);
        expect(r.numIntersections).toBe(1);
        expect(r.parameter[0]).toBe(r.parameter[1]);
        expectVectorClose(r.point[0], Vector.fromArray([2, 0]), 1e-12, 1e-12);
    });

    it('a zero-extent (degenerate) box is the point-on-line test', () => {
        const pt = box(1, 1, 1, 1);
        const d = Vector.fromArray([1, 1]);
        normalize(d);
        const through = fi.find(
            Line.fromOriginDirection(Vector.fromArray([0, 0]), d), pt);
        expect(through.intersect).toBe(true);
        expect(through.numIntersections).toBe(1);
        expect(fi.find(line(0, 0, 1, 0), pt).intersect).toBe(false);
    });

    it('the exported DoQuery helpers reproduce the class results', () => {
        check(fc.tuple(lineArb, boxArb), ([l, b]) => {
            const { center, extent } = b.getCenteredForm();
            const origin = sub(l.origin, center);
            const tr = defaultIntrLine2AlignedBox2TIResult();
            intrLine2AlignedBox2TIDoQuery(origin, l.direction, extent, tr);
            expect(tr.intersect).toBe(ti.test(l, b).intersect);

            const fr = defaultIntrLine2AlignedBox2FIResult();
            intrLine2AlignedBox2FIDoQuery(origin, l.direction, extent, fr);
            const expected = fi.find(l, b);
            expect(fr.intersect).toBe(expected.intersect);
            expect(fr.numIntersections).toBe(expected.numIntersections);
            expect(fr.parameter[0]).toBe(expected.parameter[0]);
            expect(fr.parameter[1]).toBe(expected.parameter[1]);
        });
    });
});
