import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Line } from '../src/Line.js';
import { Vector, normalize } from '../src/Vector.js';
import {
    IntrLine3AlignedBox3TI,
    IntrLine3AlignedBox3FI
} from '../src/IntrLine3AlignedBox3.js';
import {
    intrLine3AlignedBox3TIDoQuery,
    intrLine3AlignedBox3FIDoQuery,
    defaultIntrLine3AlignedBox3TIResult,
    defaultIntrLine3AlignedBox3FIResult
} from '../src/IntrLine3AlignedBox3.js';
import { add, mul, sub } from '../src/Vector.js';
import { check, expectVectorClose, fc, positive, seededRandom, unitVector, wellScaled } from './helpers/arbitraries.js';

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

describe('intrLine3AlignedBox3 DoQuery helpers', () => {
    // The helpers take the line in the box-centered coordinate system. The
    // box below is centered at the origin, so no translation is needed.
    const b = box([-2, -1, -3], [2, 1, 3]);
    const extent = Vector.fromArray([2, 1, 3]);

    it('the TI helper matches the class query', () => {
        const cases: Array<[number[], number[]]> = [
            [[0, 0, 0], [1, 0, 0]],
            [[0, 5, 0], [1, 0, 0]],
            [[-4, -3, 0], [1, 1, 0]],
            [[0, 0, 8], [0, 0, -1]],
            [[-2, 1, 3], [1, 0, 0]]
        ];
        for (const [p, d] of cases) {
            const dir = Vector.fromArray(d);
            normalize(dir);
            const result = defaultIntrLine3AlignedBox3TIResult();
            intrLine3AlignedBox3TIDoQuery(Vector.fromArray(p), dir, extent,
                result);
            const expected = new IntrLine3AlignedBox3TI().test(
                line(p, [dir.values[0], dir.values[1], dir.values[2]]), b);
            expect(result.intersect).toBe(expected.intersect);
        }
    });

    it('the FI helper fills parameters but not points', () => {
        const result = defaultIntrLine3AlignedBox3FIResult();
        intrLine3AlignedBox3FIDoQuery(Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]), extent, result);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
        // DoQuery leaves 'point' at its default value.
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([0, 0, 0]);
    });

    it('the FI helper reports no intersection for a missing line', () => {
        const result = defaultIntrLine3AlignedBox3FIResult();
        intrLine3AlignedBox3FIDoQuery(Vector.fromArray([0, 5, 0]),
            Vector.fromArray([1, 0, 0]), extent, result);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });
});

describe('IntrLine3AlignedBox3 verification', () => {
    const ti = new IntrLine3AlignedBox3TI();
    const fi = new IntrLine3AlignedBox3FI();

    const boxArb = fc.tuple(
        fc.array(wellScaled(-4, 4), { minLength: 3, maxLength: 3 }),
        fc.array(positive(4), { minLength: 3, maxLength: 3 }))
        .map(([lo, e]) => box(lo, [lo[0] + e[0], lo[1] + e[1], lo[2] + e[2]]));
    const lineArb = fc.tuple(
        fc.array(wellScaled(-6, 6), { minLength: 3, maxLength: 3 }),
        unitVector(3))
        .map(([p, d]) => line(p, [...d.values]));

    function inBox(p: Vector, b: AlignedBox, tol: number): boolean {
        for (let i = 0; i < 3; ++i) {
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
            let scale = 1;
            for (let i = 0; i < 3; ++i) {
                scale += Math.abs(b.min.values[i]) + Math.abs(b.max.values[i]);
            }
            // Upstream fills both points whenever 'intersect' is true, even
            // for the single-point case (parameter[1] == parameter[0]).
            for (let k = 0; k < 2; ++k) {
                const p = r.point[k];
                expectVectorClose(p,
                    add(l.origin, mul(r.parameter[k], l.direction)), 0, 0);
                expect(inBox(p, b, 1e-12 * scale)).toBe(true);
            }
        });
    });

    it('the reported [t0,t1] contains every line parameter strictly inside the box', () => {
        const rnd = seededRandom(0x77c3ae1);
        check(fc.tuple(lineArb, boxArb), ([l, b]) => {
            const r = fi.find(l, b);
            let scale = 1;
            for (let i = 0; i < 3; ++i) {
                scale += Math.abs(b.min.values[i]) + Math.abs(b.max.values[i]);
            }
            for (let k = 0; k < 400; ++k) {
                const t = 24 * rnd() - 12;
                const p = add(l.origin, mul(t, l.direction));
                if (inBox(p, b, -1e-9 * scale)) {
                    expect(r.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(r.parameter[0] - 1e-9 * scale);
                    expect(t).toBeLessThanOrEqual(r.parameter[1] + 1e-9 * scale);
                }
            }
        }, 60);
    }, 30000);

    it('a robustly interior hit survives a translation of line and box', () => {
        check(fc.tuple(lineArb, boxArb,
            fc.array(wellScaled(-5, 5), { minLength: 3, maxLength: 3 })),
            ([l, b, t]) => {
                const a = fi.find(l, b);
                if (!a.intersect || a.numIntersections !== 2) {
                    return;
                }
                // The clip is done in the box-centered form, so a translation
                // perturbs every difference it forms. Only a hit that passes
                // strictly through the interior is guaranteed to survive; the
                // margin below rejects lines that merely graze a face, which
                // legitimately flip under an arbitrarily small perturbation.
                let scale = 1;
                for (let i = 0; i < 3; ++i) {
                    scale += Math.abs(b.min.values[i]) + Math.abs(b.max.values[i])
                        + Math.abs(t[i]);
                }
                const mid = add(l.origin,
                    mul(0.5 * (a.parameter[0] + a.parameter[1]), l.direction));
                let margin = Infinity;
                for (let i = 0; i < 3; ++i) {
                    margin = Math.min(margin, mid.values[i] - b.min.values[i],
                        b.max.values[i] - mid.values[i]);
                }
                if (margin <= 1e-6 * scale) {
                    return;
                }
                const shift = (v: Vector): Vector => Vector.fromArray(
                    [v.values[0] + t[0], v.values[1] + t[1], v.values[2] + t[2]]);
                const c = fi.find(
                    Line.fromOriginDirection(shift(l.origin), l.direction),
                    AlignedBox.fromMinMax(shift(b.min), shift(b.max)));
                expect(c.intersect).toBe(true);
                expect(Math.abs(a.parameter[0] - c.parameter[0]))
                    .toBeLessThanOrEqual(1e-9 * scale);
                expect(Math.abs(a.parameter[1] - c.parameter[1]))
                    .toBeLessThanOrEqual(1e-9 * scale);
            });
    });

    it('a line grazing a box corner reports a single point', () => {
        const b = box([0, 0, 0], [2, 2, 2]);
        const d = Vector.fromArray([1, 1, 0]);
        normalize(d);
        const r = fi.find(
            Line.fromOriginDirection(Vector.fromArray([2, 0, 1]), d), b);
        expect(r.intersect).toBe(true);
        expect(r.numIntersections).toBe(1);
        expect(r.parameter[0]).toBe(r.parameter[1]);
        expectVectorClose(r.point[0], r.point[1], 0, 0);
    });

    it('the exported DoQuery helpers reproduce the class results', () => {
        check(fc.tuple(lineArb, boxArb), ([l, b]) => {
            const { center, extent } = b.getCenteredForm();
            const origin = sub(l.origin, center);
            const tr = defaultIntrLine3AlignedBox3TIResult();
            intrLine3AlignedBox3TIDoQuery(origin, l.direction, extent, tr);
            expect(tr.intersect).toBe(ti.test(l, b).intersect);

            const fr = defaultIntrLine3AlignedBox3FIResult();
            intrLine3AlignedBox3FIDoQuery(origin, l.direction, extent, fr);
            const expected = fi.find(l, b);
            expect(fr.intersect).toBe(expected.intersect);
            expect(fr.numIntersections).toBe(expected.numIntersections);
            expect(fr.parameter[0]).toBe(expected.parameter[0]);
            expect(fr.parameter[1]).toBe(expected.parameter[1]);
        });
    });
});
